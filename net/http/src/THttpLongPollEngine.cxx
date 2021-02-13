// $Id$
// Author: Sergey Linev   8/01/2018

/*************************************************************************
 * Copyright (C) 1995-2013, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "THttpLongPollEngine.h"

#include "TError.h"
#include "THttpCallArg.h"

#include <cstring>
#include <cstdlib>

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// THttpLongPollEngine                                                  //
//                                                                      //
// Emulation of websocket with long poll requests                       //
// Allows to send data from server to client without explicit request   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

const std::string THttpLongPollEngine::gLongPollNope = "<<nope>>";

//////////////////////////////////////////////////////////////////////////
/// constructor

THttpLongPollEngine::THttpLongPollEngine(bool raw) : THttpWSEngine(), fRaw(raw)
{
}

//////////////////////////////////////////////////////////////////////////
/// returns ID of the engine, created from this pointer

UInt_t THttpLongPollEngine::GetId() const
{
   const void *ptr = (const void *)this;
   return TString::Hash((void *)&ptr, sizeof(void *));
}

//////////////////////////////////////////////////////////////////////////
/// clear request, normally called shortly before destructor

void THttpLongPollEngine::ClearHandle(Bool_t)
{
   std::shared_ptr<THttpCallArg> poll;

   {
      std::lock_guard<std::mutex> grd(fMutex);
      poll = std::move(fPoll);
   }

   if (poll) {
      poll->Set404();
      poll->NotifyCondition();
   }
}

//////////////////////////////////////////////////////////////////////////
/// Create raw buffer which should be send as reply
/// For the raw mode all information must be send via binary response

std::string THttpLongPollEngine::MakeBuffer(const void *buf, int len, const char *hdr)
{
   std::string res;

   if (!fRaw) {
      res.resize(len);
      std::copy((const char *)buf, (const char *)buf + len, res.begin());
      return res;
   }

   int hdrlen = hdr ? strlen(hdr) : 0;
   std::string hdrstr = "bin:";
   hdrstr.append(std::to_string(hdrlen));

   while ((hdrstr.length() + 1 + hdrlen) % 8 != 0)
      hdrstr.append(" ");
   hdrstr.append(":");
   if (hdrlen > 0)
      hdrstr.append(hdr);

   res.resize(hdrstr.length() + len);
   std::copy(hdrstr.begin(), hdrstr.begin() + hdrstr.length(), res.begin());
   std::copy((const char *)buf, (const char *)buf + len, res.begin() + hdrstr.length());

   return res;
}

//////////////////////////////////////////////////////////////////////////
/// Send binary data via connection

void THttpLongPollEngine::Send(const void *buf, int len)
{
   std::shared_ptr<THttpCallArg> poll;

   std::string sendbuf = MakeBuffer(buf, len);

   {
      std::lock_guard<std::mutex> grd(fMutex);
      if (fPoll) {
         poll = std::move(fPoll);
      } else if (fBufKind == kNoBuf)  {
         fBufKind = kBinBuf;
         std::swap(fBuf, sendbuf);
         return;
      }
   }

   if(!poll) {
      Error("Send", "Operation invoked before polling request obtained");
      return;
   }

   poll->SetBinaryContent(std::move(sendbuf));
   poll->NotifyCondition();
}

//////////////////////////////////////////////////////////////////////////
/// Send binary data with text header via connection

void THttpLongPollEngine::SendHeader(const char *hdr, const void *buf, int len)
{
   std::shared_ptr<THttpCallArg> poll;

   std::string sendbuf = MakeBuffer(buf, len, hdr);

   {
      std::lock_guard<std::mutex> grd(fMutex);
      if (fPoll) {
         poll = std::move(fPoll);
      } else if (fBufKind == kNoBuf)  {
         fBufKind = kBinBuf;
         if (!fRaw && hdr) fBufHeader = hdr;
         std::swap(fBuf, sendbuf);
         return;
      }
   }

   if(!poll) {
      Error("SendHeader", "Operation invoked before polling request obtained");
      return;
   }

   poll->SetBinaryContent(std::move(sendbuf));
   if (!fRaw)
      poll->SetExtraHeader("LongpollHeader", hdr);
   poll->NotifyCondition();
}

//////////////////////////////////////////////////////////////////////////
/// Send const char data
/// Either do it immediately or keep in internal buffer

void THttpLongPollEngine::SendCharStar(const char *buf)
{
   std::shared_ptr<THttpCallArg> poll;

   std::string sendbuf(fRaw ? "txt:" : "");
   sendbuf.append(buf);

   {
      std::lock_guard<std::mutex> grd(fMutex);
      if (fPoll) {
         poll = std::move(fPoll);
      } else if (fBufKind == kNoBuf) {
         fBufKind = fRaw ? kBinBuf : kTxtBuf;
         std::swap(fBuf, sendbuf);
         return;
      }
   }

   if(!poll) {
      Error("SendCharStart", "Operation invoked before polling request obtained");
      return;
   }

   if (fRaw) poll->SetBinaryContent(std::move(sendbuf));
        else poll->SetTextContent(std::move(sendbuf));
   poll->NotifyCondition();
}

//////////////////////////////////////////////////////////////////////////////
/// Preview data for given socket
/// Method called by WS handler before processing websocket data
/// Returns kTRUE when user should ignore such http request - it is for internal use

Bool_t THttpLongPollEngine::PreProcess(std::shared_ptr<THttpCallArg> &arg)
{
   if (!strstr(arg->GetQuery(), "&dummy"))
      return kFALSE;

   std::shared_ptr<THttpCallArg> poll;

   // if request cannot be postponed just reply it
   if (arg->CanPostpone()) {
      std::lock_guard<std::mutex> grd(fMutex);
      if (fBufKind != kNoBuf) {
         // there is data which can be send directly
         poll = arg;
      } else {
         arg->SetPostponed(); // mark http request as pending, http server should wait for notification
         poll = std::move(fPoll);
         fPoll = arg;         // keep reference on polling request
      }

   } else {
      poll = arg;
   }

   if (poll) {
      // if there buffered data, it will be provided
      PostProcess(poll);
      poll->NotifyCondition();         // inform http server that request is processed
   }

   // if arguments has "&dummy" string, user should not process it
   return kTRUE;
}

//////////////////////////////////////////////////////////////////////////////
/// Normally requests from client does not replied directly for longpoll socket
/// Therefore one can use such request to send data, which was submitted before to the queue

void THttpLongPollEngine::PostProcess(std::shared_ptr<THttpCallArg> &arg)
{
   EBufKind kind = kNoBuf;
   std::string sendbuf, sendhdr;

   {
      std::lock_guard<std::mutex> grd(fMutex);
      if (fBufKind != kNoBuf) {
         kind = fBufKind;
         fBufKind = kNoBuf;
         std::swap(sendbuf, fBuf);
         std::swap(sendhdr, fBufHeader);
      }
   }

   if (kind == kTxtBuf) {
      arg->SetTextContent(std::move(sendbuf));
   } else if (kind == kBinBuf) {
      arg->SetBinaryContent(std::move(sendbuf));
      if (!sendhdr.empty())
         arg->SetExtraHeader("LongpollHeader", sendhdr.c_str());
   } else if (fRaw) {
      arg->SetBinaryContent(std::string("txt:") + gLongPollNope);
   } else {
      arg->SetTextContent(std::string(gLongPollNope));
   }
}

//////////////////////////////////////////////////////////////////////////////
/// Indicate that polling requests is there or buffer empty and can be immediately invoked

Bool_t THttpLongPollEngine::CanSendDirectly()
{
   std::lock_guard<std::mutex> grd(fMutex);
   return fPoll || (fBufKind == kNoBuf) ? kTRUE : kFALSE;
}
