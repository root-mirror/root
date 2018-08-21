// $Id$
// Author: Sergey Linev   20/10/2017

/*************************************************************************
 * Copyright (C) 1995-2013, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "THttpWSHandler.h"

#include "THttpWSEngine.h"
#include "THttpCallArg.h"

#include <thread>

/////////////////////////////////////////////////////////////////////////
///
/// THttpWSHandler
///
/// Class for user-side handling of websocket with THttpServer
/// 1. Create derived from  THttpWSHandler class and implement
///     ProcessWS() method, where all web sockets request handled.
/// 2. Register instance of derived class to running THttpServer
///
///        TUserWSHandler *handler = new TUserWSHandler("name1","title");
///        THttpServer *server = new THttpServer("http:8090");
///        server->Register("/subfolder", handler)
///
/// 3. Now server can accept web socket connection from outside.
///    For instance, from JavaScirpt one can connect to it with code:
///
///        var ws = new WebSocket("ws://hostname:8090/subfolder/name1/root.websocket")
///
/// 4. In the ProcessWS(THttpCallArg *arg) method following code should be implemented:
///
///     if (arg->IsMethod("WS_CONNECT")) {
///         return true;  // to accept incoming request
///      }
///
///      if (arg->IsMethod("WS_READY")) {
///          fWSId = arg->GetWSId(); // fWSId should be member of the user class
///          return true; // connection established
///      }
///
///     if (arg->IsMethod("WS_CLOSE")) {
///         fWSId = 0;
///         return true; // confirm close of socket
///     }
///
///     if (arg->IsMethod("WS_DATA")) {
///         // received data stored as POST data
///         std::string str((const char *)arg->GetPostData(), arg->GetPostDataLength());
///         std::cout << "got string " << str << std::endl;
///         // immediately send data back using websocket id
///         SendCharStarWS(fWSId, "our reply");
///         return true;
///     }
///
///////////////////////////////////////////////////////////////////////////

ClassImp(THttpWSHandler);

////////////////////////////////////////////////////////////////////////////////
/// normal constructor

THttpWSHandler::THttpWSHandler(const char *name, const char *title) : TNamed(name, title)
{
}

////////////////////////////////////////////////////////////////////////////////
/// destructor
/// Delete all websockets handles

THttpWSHandler::~THttpWSHandler()
{
   SetDisabled();
}

/// Returns current number of websocket connections
Int_t THttpWSHandler::GetNumWS()
{
   std::lock_guard<std::mutex> grd(fMutex);
   return fEngines.size();
}

////////////////////////////////////////////////////////////////////////////////
/// Return websocket id with given sequential number
/// Number of websockets returned with GetNumWS() method

UInt_t THttpWSHandler::GetWS(Int_t num)
{
   std::lock_guard<std::mutex> grd(fMutex);
   auto iter = fEngines.begin() + num;
   return (*iter)->GetId();
}

////////////////////////////////////////////////////////////////////////////////
/// Find websocket connection handle with given id
/// If book_send parameter specified, have to book send operation under the mutex

std::shared_ptr<THttpWSEngine> THttpWSHandler::FindEngine(UInt_t wsid, Bool_t book_send)
{
   if (IsDisabled())
      return nullptr;

   std::lock_guard<std::mutex> grd(fMutex);

   for (auto &eng : fEngines)
      if (eng->GetId() == wsid) {

         // not allow to work with disabled engine
         if (eng->fDisabled)
            return nullptr;

         if (book_send) {
            if (eng->fMTSend) {
               Error("FindEngine", "Try to book next send operation before previous completed");
               return nullptr;
            }
            eng->fMTSend = kTRUE;
         }
         return eng;
      }

   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove and destroy WS connection

void THttpWSHandler::RemoveEngine(std::shared_ptr<THttpWSEngine> &engine)
{
   {
      std::lock_guard<std::mutex> grd(fMutex);

      for (auto iter = fEngines.begin(); iter != fEngines.end(); iter++)
         if (*iter == engine) {
            if (engine->fMTSend)
               Error("RemoveEngine", "Trying to remove WS engine during send operation");

            engine->fDisabled = true;
            fEngines.erase(iter);
            break;
         }
   }

   engine->ClearHandle();
}

////////////////////////////////////////////////////////////////////////////////
/// Process request to websocket
/// Different kind of requests coded into THttpCallArg::Method
///  "WS_CONNECT" - connection request
///  "WS_READY" - connection ready
///  "WS_CLOSE" - connection closed
/// All other are normal data, which are delivered to users

Bool_t THttpWSHandler::HandleWS(std::shared_ptr<THttpCallArg> &arg)
{
   if (IsDisabled())
      return kFALSE;

   if (!arg->GetWSId())
      return ProcessWS(arg.get());

   // normally here one accept or reject connection requests
   if (arg->IsMethod("WS_CONNECT"))
      return ProcessWS(arg.get());

   auto engine = FindEngine(arg->GetWSId());

   if (arg->IsMethod("WS_READY")) {

      if (engine) {
         Error("HandleWS", "WS engine with similar id exists %u", arg->GetWSId());
         RemoveEngine(engine);
      }

      engine = arg->TakeWSEngine();
      {
         std::lock_guard<std::mutex> grd(fMutex);
         fEngines.emplace_back(engine);
      }

      if (!ProcessWS(arg.get())) {
         // if connection refused, remove engine again
         RemoveEngine(engine);
         return kFALSE;
      }

      return kTRUE;
   }

   if (arg->IsMethod("WS_CLOSE")) {
      // connection is closed, one can remove handle

      if (engine) {
         engine->ClearHandle();
         RemoveEngine(engine);
      }

      return ProcessWS(arg.get());
   }

   Bool_t check_send  = engine ? engine->PreviewData(arg) : kFALSE;

   Bool_t res = kTRUE;

   if (!check_send) {

      res = ProcessWS(arg.get());

      check_send = engine ? engine->PostProcess(arg) : kFALSE;
   }

   if (check_send)
      PerformSend(engine);


   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Close connection with given websocket id

void THttpWSHandler::CloseWS(UInt_t wsid)
{
   auto engine = FindEngine(wsid);

   if (engine)
      RemoveEngine(engine);
}

////////////////////////////////////////////////////////////////////////////////
/// Send binary data via given websocket id
/// Returns -1 - in case of error
///          0 - when operation was executed immediately
///          1 - when send operation will be performed in different thread

Int_t THttpWSHandler::SendWS(UInt_t wsid, const void *buf, int len)
{
   auto engine = FindEngine(wsid, kTRUE);
   if (!engine) return -1;

   if (!AllowMTSend() && engine->CanSendDirectly()) {
      engine->Send(buf, len);
      engine->fMTSend = false; // probably we do not need to lock mutex to reset flag
      CompleteWSSend(engine->GetId());
      return 0;
   }

   // now we indicate that there is data and any thread can access it
   {
      std::lock_guard<std::mutex> grd(engine->fDataMutex);

      if (engine->fKind != THttpWSEngine::kNone) {
         Error("SendWS", "Data kind is not empty - something screwed up");
         return -1;
      }

      engine->fData.resize(len);
      std::copy((const char *)buf, (const char *)buf + len, engine->fData.begin());

      engine->fDoingSend = false;
      engine->fKind = THttpWSEngine::kData;
   }

   return RunSendingThrd(engine);
}

////////////////////////////////////////////////////////////////////////////////
/// Send data stored in the buffer
/// Returns 0 - when operation was executed immediately
///         1 - when send operation will be performed in different thread

Int_t THttpWSHandler::RunSendingThrd(std::shared_ptr<THttpWSEngine> engine)
{
   // actually lonpoll engine does not require thread to reply data in buffer
   if (!engine->RequireSendThrd()) {

      if (engine->CanSendDirectly())
         return PerformSend(engine);

      // handling will be performed in http request handler
      return 1;
   }

   std::thread thrd([this, engine] {
      PerformSend(engine);
   });

   thrd.detach(); // let continue thread execution without thread handle

   return 1;
}


////////////////////////////////////////////////////////////////////////////////
/// Perform send operation, stored in buffer

Int_t THttpWSHandler::PerformSend(std::shared_ptr<THttpWSEngine> engine)
{
   {
      std::lock_guard<std::mutex> grd(engine->fDataMutex);

      // no need to do somthing - operation was processed already by somebody else
      if (engine->fKind == THttpWSEngine::kNone)
         return 0;

      if (engine->fDoingSend)
         return 1;
      engine->fDoingSend = true;
   }

   if (IsDisabled() || engine->fDisabled)
      return 0;

   switch (engine->fKind) {
   case THttpWSEngine::kData:
      engine->Send(engine->fData.data(), engine->fData.length());
      break;
   case THttpWSEngine::kHeader:
      engine->SendHeader(engine->fHdr.c_str(), engine->fData.data(), engine->fData.length());
      break;
   case THttpWSEngine::kText:
      engine->SendCharStar(engine->fData.c_str());
      break;
   default:
      break;
   }

   engine->fData.clear();
   engine->fHdr.clear();

   {
      std::lock_guard<std::mutex> grd(engine->fDataMutex);
      engine->fDoingSend = false;
      engine->fKind = THttpWSEngine::kNone;
   }

   engine->fMTSend = false; // probably we do not need to lock mutex to reset flag
   CompleteWSSend(engine->GetId());

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Send binary data with text header via given websocket id
/// Returns -1 - in case of error,
///          0 - when operation was executed immediately,
///          1 - when send operation will be performed in different thread,

Int_t THttpWSHandler::SendHeaderWS(UInt_t wsid, const char *hdr, const void *buf, int len)
{
   auto engine = FindEngine(wsid, kTRUE);
   if (!engine) return -1;

   if (!AllowMTSend() && engine->CanSendDirectly()) {
      engine->SendHeader(hdr, buf, len);
      engine->fMTSend = false; // probably we do not need to lock mutex to reset flag
      CompleteWSSend(engine->GetId());
      return 0;
   }


   // now we indicate that there is data and any thread can access it
   {
      std::lock_guard<std::mutex> grd(engine->fDataMutex);

      if (engine->fKind != THttpWSEngine::kNone) {
         Error("SendWS", "Data kind is not empty - something screwed up");
         return -1;
      }

      engine->fHdr = hdr;
      engine->fData.resize(len);
      std::copy((const char *)buf, (const char *)buf + len, engine->fData.begin());

      engine->fDoingSend = false;
      engine->fKind = THttpWSEngine::kHeader;
   }

   return RunSendingThrd(engine);
}

////////////////////////////////////////////////////////////////////////////////
/// Send string via given websocket id
/// Returns -1 - in case of error,
///          0 - when operation was executed immediately,
///          1 - when send operation will be performed in different thread,

Int_t THttpWSHandler::SendCharStarWS(UInt_t wsid, const char *str)
{
   auto engine = FindEngine(wsid, kTRUE);
   if (!engine) return -1;

   if (!AllowMTSend() && engine->CanSendDirectly()) {
      engine->SendCharStar(str);
      engine->fMTSend = false; // probably we do not need to lock mutex to reset flag
      CompleteWSSend(engine->GetId());
      return 0;
   }

   // now we indicate that there is data and any thread can access it
   {
      std::lock_guard<std::mutex> grd(engine->fDataMutex);

      if (engine->fKind != THttpWSEngine::kNone) {
         Error("SendWS", "Data kind is not empty - something screwed up");
         return -1;
      }

      engine->fData = str;

      engine->fDoingSend = false;
      engine->fKind = THttpWSEngine::kText;
   }

   return RunSendingThrd(engine);
}
