// Author: Sergey Linev <s.linev@gsi.de>
// Date: 2017-10-16
// Warning: This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/RWebWindow.hxx>

#include <ROOT/RWebWindowsManager.hxx>
#include <ROOT/RLogger.hxx>

#include "RWebWindowWSHandler.hxx"
#include "THttpCallArg.h"
#include "TUrl.h"
#include "TROOT.h"

#include <cstring>
#include <cstdlib>
#include <utility>
#include <assert.h>
#include <algorithm>
#include <fstream>

using namespace ROOT::Experimental;
using namespace std::string_literals;

//////////////////////////////////////////////////////////////////////////////////////////
/// Destructor for WebConn
/// Notify special HTTP request which blocks headless browser from exit

RWebWindow::WebConn::~WebConn()
{
   if (fHold) {
      fHold->SetTextContent("console.log('execute holder script');  if (window) setTimeout (window.close, 1000); if (window) window.close();");
      fHold->NotifyCondition();
      fHold.reset();
   }
}


/** \class ROOT::Experimental::RWebWindow
\ingroup webdisplay

Represents web window, which can be shown in web browser or any other supported environment

Window can be configured to run either in the normal or in the batch (headless) mode.
In second case no any graphical elements will be created. For the normal window one can configure geometry
(width and height), which are applied when window shown.

Each window can be shown several times (if allowed) in different places - either as the
CEF (chromium embedded) window or in the standard web browser. When started, window will open and show
HTML page, configured with RWebWindow::SetDefaultPage() method.

Typically (but not necessarily) clients open web socket connection to the window and one can exchange data,
using RWebWindow::Send() method and call-back function assigned via RWebWindow::SetDataCallBack().

*/


//////////////////////////////////////////////////////////////////////////////////////////
/// RWebWindow constructor
/// Should be defined here because of std::unique_ptr<RWebWindowWSHandler>

RWebWindow::RWebWindow() = default;

//////////////////////////////////////////////////////////////////////////////////////////
/// RWebWindow destructor
/// Closes all connections and remove window from manager

RWebWindow::~RWebWindow()
{
   if (fMaster)
      fMaster->RemoveEmbedWindow(fMasterConnId, fMasterChannel);

   if (fWSHandler)
      fWSHandler->SetDisabled();

   if (fMgr) {

      // make copy of all connections
      auto lst = GetConnections();

      {
         // clear connections vector under mutex
         std::lock_guard<std::mutex> grd(fConnMutex);
         fConn.clear();
         fPendingConn.clear();
      }

      for (auto &conn : lst) {
         conn->fActive = false;
         for (auto &elem: conn->fEmbed)
            elem.second->fMaster.reset();
      }

      fMgr->Unregister(*this);
   }
}

//////////////////////////////////////////////////////////////////////////////////////////
/// Configure window to show some of existing JSROOT panels
/// It uses "file:rootui5sys/panel/panel.html" as default HTML page
/// At the moment only FitPanel is existing

void RWebWindow::SetPanelName(const std::string &name)
{
   {
      std::lock_guard<std::mutex> grd(fConnMutex);
      if (!fConn.empty()) {
         R__LOG_ERROR(WebGUILog()) << "Cannot configure panel when connection exists";
         return;
      }
   }

   fPanelName = name;
   SetDefaultPage("file:rootui5sys/panel/panel.html");
}

//////////////////////////////////////////////////////////////////////////////////////////
/// Assigns manager reference, window id and creates websocket handler, used for communication with the clients

std::shared_ptr<RWebWindowWSHandler>
RWebWindow::CreateWSHandler(std::shared_ptr<RWebWindowsManager> mgr, unsigned id, double tmout)
{
   fMgr = mgr;
   fId = id;
   fOperationTmout = tmout;

   fSendMT = fMgr->IsUseSenderThreads();
   fWSHandler = std::make_shared<RWebWindowWSHandler>(*this, Form("win%u", GetId()));

   return fWSHandler;
}

//////////////////////////////////////////////////////////////////////////////////////////
/// Return URL string to access web window
/// If remote flag is specified, real HTTP server will be started automatically

std::string RWebWindow::GetUrl(bool remote)
{
   return fMgr->GetUrl(*this, remote);
}

//////////////////////////////////////////////////////////////////////////////////////////
/// Return THttpServer instance serving requests to the window

THttpServer *RWebWindow::GetServer()
{
   return fMgr->GetServer();
}

//////////////////////////////////////////////////////////////////////////////////////////
/// Show window in specified location
/// See ROOT::Experimental::RWebWindowsManager::Show() docu for more info
/// returns (future) connection id (or 0 when fails)

unsigned RWebWindow::Show(const RWebDisplayArgs &args)
{
   return fMgr->ShowWindow(*this, false, args);
}

//////////////////////////////////////////////////////////////////////////////////////////
/// Create batch job for specified window
/// Normally only single batch job is used, but many can be created
/// See ROOT::Experimental::RWebWindowsManager::Show() docu for more info
/// returns (future) connection id (or 0 when fails)

unsigned RWebWindow::MakeBatch(bool create_new, const RWebDisplayArgs &args)
{
   unsigned connid = 0;
   if (!create_new)
      connid = FindBatch();
   if (!connid)
      connid = fMgr->ShowWindow(*this, true, args);
   return connid;
}

//////////////////////////////////////////////////////////////////////////////////////////
/// Returns connection id of batch job
/// Connection to that job may not be initialized yet
/// If connection does not exists, returns 0

unsigned RWebWindow::FindBatch()
{
   std::lock_guard<std::mutex> grd(fConnMutex);

   for (auto &entry : fPendingConn) {
      if (entry->fBatchMode)
         return entry->fConnId;
   }

   for (auto &conn : fConn) {
      if (conn->fBatchMode)
         return conn->fConnId;
   }

   return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////
/// Returns first connection id where window is displayed
/// It could be that connection(s) not yet fully established - but also not timed out
/// Batch jobs will be ignored here
/// Returns 0 if connection not exists

unsigned RWebWindow::GetDisplayConnection() const
{
   std::lock_guard<std::mutex> grd(fConnMutex);

   for (auto &entry : fPendingConn) {
      if (!entry->fBatchMode)
         return entry->fConnId;
   }

   for (auto &conn : fConn) {
      if (!conn->fBatchMode)
         return conn->fConnId;
   }

   return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////
/// Find connection with given websocket id
/// Connection mutex should be locked before method calling

std::shared_ptr<RWebWindow::WebConn> RWebWindow::FindOrCreateConnection(unsigned wsid, bool make_new, const char *query)
{
   std::lock_guard<std::mutex> grd(fConnMutex);

   for (auto &conn : fConn) {
      if (conn->fWSId == wsid)
         return conn;
   }

   // put code to create new connection here to stay under same locked mutex
   if (make_new) {
      // check if key was registered already

      std::shared_ptr<WebConn> key;
      std::string keyvalue;

      if (query) {
         TUrl url;
         url.SetOptions(query);
         if (url.HasOption("key"))
            keyvalue = url.GetValueFromOptions("key");
      }

      if (!keyvalue.empty())
         for (size_t n = 0; n < fPendingConn.size(); ++n)
            if (fPendingConn[n]->fKey == keyvalue) {
               key = std::move(fPendingConn[n]);
               fPendingConn.erase(fPendingConn.begin() + n);
               break;
            }

      if (key) {
         key->fWSId = wsid;
         key->fActive = true;
         key->ResetStamps(); // TODO: probably, can be moved outside locked area
         fConn.emplace_back(key);
      } else {
         fConn.emplace_back(std::make_shared<WebConn>(++fConnCnt, wsid));
      }
   }

   return nullptr;
}

//////////////////////////////////////////////////////////////////////////////////////////
/// Remove connection with given websocket id

std::shared_ptr<RWebWindow::WebConn> RWebWindow::RemoveConnection(unsigned wsid)
{

   std::shared_ptr<WebConn> res;

   {
      std::lock_guard<std::mutex> grd(fConnMutex);

      for (size_t n = 0; n < fConn.size(); ++n)
         if (fConn[n]->fWSId == wsid) {
            res = std::move(fConn[n]);
            fConn.erase(fConn.begin() + n);
            res->fActive = false;
            break;
         }
   }

   if (res)
      for (auto &elem: res->fEmbed)
         elem.second->fMaster.reset();

   return res;
}

//////////////////////////////////////////////////////////////////////////////////////////
/// Process special http request, used to hold headless browser running
/// Such requests should not be replied for the long time
/// Be aware that function called directly from THttpServer thread, which is not same thread as window

bool RWebWindow::ProcessBatchHolder(std::shared_ptr<THttpCallArg> &arg)
{
   std::string query = arg->GetQuery();

   if (query.compare(0, 4, "key=") != 0)
      return false;

   std::string key = query.substr(4);

   std::shared_ptr<THttpCallArg> prev;

   bool found_key = false;

   // use connection mutex to access hold request
   {
      std::lock_guard<std::mutex> grd(fConnMutex);
      for (auto &entry : fPendingConn) {
         if (entry->fKey == key)  {
            assert(!found_key); // indicate error if many same keys appears
            found_key = true;
            prev = std::move(entry->fHold);
            entry->fHold = arg;
         }
      }


      for (auto &conn : fConn) {
         if (conn->fKey == key) {
            assert(!found_key); // indicate error if many same keys appears
            prev = std::move(conn->fHold);
            conn->fHold = arg;
            found_key = true;
         }
      }
   }

   if (prev) {
      prev->SetTextContent("console.log('execute holder script'); if (window) window.close();");
      prev->NotifyCondition();
   }

   return found_key;
}

//////////////////////////////////////////////////////////////////////////////////////////
/// Provide data to user callback
/// User callback must be executed in the window thread

void RWebWindow::ProvideQueueEntry(unsigned connid, EQueueEntryKind kind, std::string &&arg)
{
   {
      std::lock_guard<std::mutex> grd(fInputQueueMutex);
      fInputQueue.emplace(connid, kind, std::move(arg));
   }

   InvokeCallbacks();
}

//////////////////////////////////////////////////////////////////////////////////////////
/// Invoke callbacks with existing data
/// Must be called from appropriate thread

void RWebWindow::InvokeCallbacks(bool force)
{
   if (fCallbacksThrdIdSet && (fCallbacksThrdId != std::this_thread::get_id()) && !force)
      return;

   while (true) {
      unsigned connid;
      EQueueEntryKind kind;
      std::string arg;

      {
         std::lock_guard<std::mutex> grd(fInputQueueMutex);
         if (fInputQueue.size() == 0)
            return;
         auto &entry = fInputQueue.front();
         connid = entry.fConnId;
         kind = entry.fKind;
         arg = std::move(entry.fData);
         fInputQueue.pop();
      }

      switch (kind) {
      case kind_None: break;
      case kind_Connect:
         if (fConnCallback)
            fConnCallback(connid);
         break;
      case kind_Data:
         if (fDataCallback)
            fDataCallback(connid, arg);
         break;
      case kind_Disconnect:
         if (fDisconnCallback)
            fDisconnCallback(connid);
         break;
      }
   }
}

//////////////////////////////////////////////////////////////////////////////////////////
/// Add display handle and associated key
/// Key is random number generated when starting new window
/// When client is connected, key should be supplied to correctly identify it

unsigned RWebWindow::AddDisplayHandle(bool batch_mode, const std::string &key, std::unique_ptr<RWebDisplayHandle> &handle)
{
   std::lock_guard<std::mutex> grd(fConnMutex);

   ++fConnCnt;

   auto conn = std::make_shared<WebConn>(fConnCnt, batch_mode, key);

   std::swap(conn->fDisplayHandle, handle);

   fPendingConn.emplace_back(conn);

   return fConnCnt;
}

//////////////////////////////////////////////////////////////////////////////////////////
/// Returns true if provided key value already exists (in processes map or in existing connections)

bool RWebWindow::HasKey(const std::string &key) const
{
   std::lock_guard<std::mutex> grd(fConnMutex);

   for (auto &entry : fPendingConn) {
      if (entry->fKey == key)
         return true;
   }

   for (auto &conn : fConn) {
      if (conn->fKey == key)
         return true;
   }

   return false;
}

//////////////////////////////////////////////////////////////////////////////////////////
/// Check if started process(es) establish connection. After timeout such processed will be killed
/// Method invoked from http server thread, therefore appropriate mutex must be used on all relevant data

void RWebWindow::CheckPendingConnections()
{
   if (!fMgr) return;

   timestamp_t stamp = std::chrono::system_clock::now();

   float tmout = fMgr->GetLaunchTmout();

   ConnectionsList_t selected;

   {
      std::lock_guard<std::mutex> grd(fConnMutex);

      auto pred = [&](std::shared_ptr<WebConn> &e) {
         std::chrono::duration<double> diff = stamp - e->fSendStamp;

         if (diff.count() > tmout) {
            R__LOG_DEBUG(0, WebGUILog()) << "Halt process after " << diff.count() << " sec";
            selected.emplace_back(e);
            return true;
         }

         return false;
      };

      fPendingConn.erase(std::remove_if(fPendingConn.begin(), fPendingConn.end(), pred), fPendingConn.end());
   }

}


//////////////////////////////////////////////////////////////////////////////////////////
/// Check if there are connection which are inactive for longer time
/// For instance, batch browser will be stopped if no activity for 30 sec is there

void RWebWindow::CheckInactiveConnections()
{
   timestamp_t stamp = std::chrono::system_clock::now();

   double batch_tmout = 20.;

   std::vector<std::shared_ptr<WebConn>> clr;

   {
      std::lock_guard<std::mutex> grd(fConnMutex);

      auto pred = [&](std::shared_ptr<WebConn> &conn) {
         std::chrono::duration<double> diff = stamp - conn->fSendStamp;
         // introduce large timeout
         if ((diff.count() > batch_tmout) && conn->fBatchMode) {
            conn->fActive = false;
            clr.emplace_back(conn);
            return true;
         }
         return false;
      };

      fConn.erase(std::remove_if(fConn.begin(), fConn.end(), pred), fConn.end());
   }

   for (auto &entry : clr)
      ProvideQueueEntry(entry->fConnId, kind_Disconnect, ""s);

}

/////////////////////////////////////////////////////////////////////////
/// Configure maximal number of allowed connections - 0 is unlimited
/// Will not affect already existing connections
/// Default is 1 - the only client is allowed

void RWebWindow::SetConnLimit(unsigned lmt)
{
   std::lock_guard<std::mutex> grd(fConnMutex);

   fConnLimit = lmt;
}

/////////////////////////////////////////////////////////////////////////
/// returns configured connections limit (0 - default)

unsigned RWebWindow::GetConnLimit() const
{
   std::lock_guard<std::mutex> grd(fConnMutex);

   return fConnLimit;
}

/////////////////////////////////////////////////////////////////////////
/// Configures connection token (default none)
/// When specified, in URL of webpage such token should be provided as &token=value parameter,
/// otherwise web window will refuse connection

void RWebWindow::SetConnToken(const std::string &token)
{
   std::lock_guard<std::mutex> grd(fConnMutex);

   fConnToken = token;
}

/////////////////////////////////////////////////////////////////////////
/// Returns configured connection token

std::string RWebWindow::GetConnToken() const
{
   std::lock_guard<std::mutex> grd(fConnMutex);

   return fConnToken;
}

//////////////////////////////////////////////////////////////////////////////////////////
/// Processing of websockets call-backs, invoked from RWebWindowWSHandler
/// Method invoked from http server thread, therefore appropriate mutex must be used on all relevant data

bool RWebWindow::ProcessWS(THttpCallArg &arg)
{
   if (arg.GetWSId() == 0)
      return true;

   if (arg.IsMethod("WS_CONNECT")) {

      std::lock_guard<std::mutex> grd(fConnMutex);

      if (!fConnToken.empty()) {
         TUrl url;
         url.SetOptions(arg.GetQuery());
         // refuse connection which does not provide proper token
         if (!url.HasOption("token") || (fConnToken != url.GetValueFromOptions("token"))) {
            R__LOG_DEBUG(0, WebGUILog()) << "Refuse connection without proper token";
            return false;
         }
      }

      // refuse connection when number of connections exceed limit
      if (fConnLimit && (fConn.size() >= fConnLimit))
         return false;

      return true;
   }

   if (arg.IsMethod("WS_READY")) {
      auto conn = FindOrCreateConnection(arg.GetWSId(), true, arg.GetQuery());

      if (conn) {
         R__LOG_ERROR(WebGUILog()) << "WSHandle with given websocket id " << arg.GetWSId() << " already exists";
         return false;
      }

      return true;
   }

   if (arg.IsMethod("WS_CLOSE")) {
      // connection is closed, one can remove handle, associated window will be closed

      auto conn = RemoveConnection(arg.GetWSId());

      if (conn)
         ProvideQueueEntry(conn->fConnId, kind_Disconnect, ""s);

      return true;
   }

   if (!arg.IsMethod("WS_DATA")) {
      R__LOG_ERROR(WebGUILog()) << "only WS_DATA request expected!";
      return false;
   }

   auto conn = FindConnection(arg.GetWSId());

   if (!conn) {
      R__LOG_ERROR(WebGUILog()) << "Get websocket data without valid connection - ignore!!!";
      return false;
   }

   if (arg.GetPostDataLength() <= 0)
      return true;

   // here processing of received data should be performed
   // this is task for the implemented windows

   const char *buf = (const char *)arg.GetPostData();
   char *str_end = nullptr;

   unsigned long ackn_oper = std::strtoul(buf, &str_end, 10);
   if (!str_end || *str_end != ':') {
      R__LOG_ERROR(WebGUILog()) << "missing number of acknowledged operations";
      return false;
   }

   unsigned long can_send = std::strtoul(str_end + 1, &str_end, 10);
   if (!str_end || *str_end != ':') {
      R__LOG_ERROR(WebGUILog()) << "missing can_send counter";
      return false;
   }

   unsigned long nchannel = std::strtoul(str_end + 1, &str_end, 10);
   if (!str_end || *str_end != ':') {
      R__LOG_ERROR(WebGUILog()) << "missing channel number";
      return false;
   }

   Long_t processed_len = (str_end + 1 - buf);

   if (processed_len > arg.GetPostDataLength()) {
      R__LOG_ERROR(WebGUILog()) << "corrupted buffer";
      return false;
   }

   std::string cdata(str_end + 1, arg.GetPostDataLength() - processed_len);

   timestamp_t stamp = std::chrono::system_clock::now();

   {
      std::lock_guard<std::mutex> grd(conn->fMutex);

      conn->fSendCredits += ackn_oper;
      conn->fRecvCount++;
      conn->fClientCredits = (int)can_send;
      conn->fRecvStamp = stamp;
   }

   if (fProtocolCnt >= 0)
      if (!fProtocolConnId || (conn->fConnId == fProtocolConnId)) {
         fProtocolConnId = conn->fConnId; // remember connection

         // record send event only for normal channel or very first message via ch0
         if ((nchannel != 0) || (cdata.find("READY=") == 0)) {
            if (fProtocol.length() > 2)
               fProtocol.insert(fProtocol.length() - 1, ",");
            fProtocol.insert(fProtocol.length() - 1, "\"send\"");

            std::ofstream pfs(fProtocolFileName);
            pfs.write(fProtocol.c_str(), fProtocol.length());
            pfs.close();
         }
      }

   if (nchannel == 0) {
      // special system channel
      if ((cdata.find("READY=") == 0) && !conn->fReady) {
         std::string key = cdata.substr(6);

         if (key.empty() && IsNativeOnlyConn()) {
            RemoveConnection(conn->fWSId);
            return false;
         }

         if (!key.empty() && !conn->fKey.empty() && (conn->fKey != key)) {
            R__LOG_ERROR(WebGUILog()) << "Key mismatch after established connection " << key << " != " << conn->fKey;
            RemoveConnection(conn->fWSId);
            return false;
         }

         if (!fPanelName.empty()) {
            // initialization not yet finished, appropriate panel should be started
            Send(conn->fConnId, "SHOWPANEL:"s + fPanelName);
            conn->fReady = 5;
         } else {
            ProvideQueueEntry(conn->fConnId, kind_Connect, ""s);
            conn->fReady = 10;
         }
      } else if (cdata.compare(0,8,"CLOSECH=") == 0) {
         int channel = std::stoi(cdata.substr(8));
         auto iter = conn->fEmbed.find(channel);
         if (iter != conn->fEmbed.end()) {
            iter->second->ProvideQueueEntry(conn->fConnId, kind_Disconnect, ""s);
            conn->fEmbed.erase(iter);
         }
      }
   } else if (fPanelName.length() && (conn->fReady < 10)) {
      if (cdata == "PANEL_READY") {
         R__LOG_DEBUG(0, WebGUILog()) << "Get panel ready " << fPanelName;
         ProvideQueueEntry(conn->fConnId, kind_Connect, ""s);
         conn->fReady = 10;
      } else {
         ProvideQueueEntry(conn->fConnId, kind_Disconnect, ""s);
         RemoveConnection(conn->fWSId);
      }
   } else if (nchannel == 1) {
      ProvideQueueEntry(conn->fConnId, kind_Data, std::move(cdata));
   } else if (nchannel > 1) {
      // process embed window
      auto embed_window = conn->fEmbed[nchannel];
      if (embed_window)
         embed_window->ProvideQueueEntry(conn->fConnId, kind_Data, std::move(cdata));
   }

   CheckDataToSend();

   return true;
}

void RWebWindow::CompleteWSSend(unsigned wsid)
{
   auto conn = FindConnection(wsid);

   if (!conn)
      return;

   {
      std::lock_guard<std::mutex> grd(conn->fMutex);
      conn->fDoingSend = false;
   }

   CheckDataToSend(conn);
}

//////////////////////////////////////////////////////////////////////////////////////////
/// Prepare text part of send data
/// Should be called under locked connection mutex

std::string RWebWindow::_MakeSendHeader(std::shared_ptr<WebConn> &conn, bool txt, const std::string &data, int chid)
{
   std::string buf;

   if (!conn->fWSId || !fWSHandler) {
      R__LOG_ERROR(WebGUILog()) << "try to send text data when connection not established";
      return buf;
   }

   if (conn->fSendCredits <= 0) {
      R__LOG_ERROR(WebGUILog()) << "No credits to send text data via connection";
      return buf;
   }

   if (conn->fDoingSend) {
      R__LOG_ERROR(WebGUILog()) << "Previous send operation not completed yet";
      return buf;
   }

   if (txt)
      buf.reserve(data.length() + 100);

   buf.append(std::to_string(conn->fRecvCount));
   buf.append(":");
   buf.append(std::to_string(conn->fSendCredits));
   buf.append(":");
   conn->fRecvCount = 0; // we confirm how many packages was received
   conn->fSendCredits--;

   buf.append(std::to_string(chid));
   buf.append(":");

   if (txt) {
      buf.append(data);
   } else if (data.length()==0) {
      buf.append("$$nullbinary$$");
   } else {
      buf.append("$$binary$$");
   }

   return buf;
}

//////////////////////////////////////////////////////////////////////////////////////////
/// Checks if one should send data for specified connection
/// Returns true when send operation was performed

bool RWebWindow::CheckDataToSend(std::shared_ptr<WebConn> &conn)
{
   std::string hdr, data;

   {
      std::lock_guard<std::mutex> grd(conn->fMutex);

      if (!conn->fActive || (conn->fSendCredits <= 0) || conn->fDoingSend) return false;

      if (!conn->fQueue.empty()) {
         QueueItem &item = conn->fQueue.front();
         hdr = _MakeSendHeader(conn, item.fText, item.fData, item.fChID);
         if (!hdr.empty() && !item.fText)
            data = std::move(item.fData);
         conn->fQueue.pop();
      } else if ((conn->fClientCredits < 3) && (conn->fRecvCount > 1)) {
         // give more credits to the client
         hdr = _MakeSendHeader(conn, true, "KEEPALIVE", 0);
      }

      if (hdr.empty()) return false;

      conn->fDoingSend = true;
   }

   int res = 0;

   if (data.empty()) {
      res = fWSHandler->SendCharStarWS(conn->fWSId, hdr.c_str());
   } else {
      res = fWSHandler->SendHeaderWS(conn->fWSId, hdr.c_str(), data.data(), data.length());
   }

   // submit operation, will be processed
   if (res >=0) return true;


   // failure, clear sending flag
   std::lock_guard<std::mutex> grd(conn->fMutex);
   conn->fDoingSend = false;
   return false;
}


//////////////////////////////////////////////////////////////////////////////////////////
/// Checks if new data can be send (internal use only)
/// If necessary, provide credits to the client

void RWebWindow::CheckDataToSend(bool only_once)
{
   // make copy of all connections to be independent later, only active connections are checked
   auto arr = GetConnections(0, true);

   do {
      bool isany = false;

      for (auto &conn : arr)
         if (CheckDataToSend(conn))
            isany = true;

      if (!isany) break;

   } while (!only_once);
}

///////////////////////////////////////////////////////////////////////////////////
/// Special method to process all internal activity when window runs in separate thread

void RWebWindow::Sync()
{
   InvokeCallbacks();

   CheckDataToSend();

   CheckPendingConnections();

   CheckInactiveConnections();
}

///////////////////////////////////////////////////////////////////////////////////
/// Returns window address which is used in URL

std::string RWebWindow::GetAddr() const
{
    return fWSHandler->GetName();
}

///////////////////////////////////////////////////////////////////////////////////
/// Returns relative URL address for the specified window
/// Address can be required if one needs to access data from one window into another window
/// Used for instance when inserting panel into canvas

std::string RWebWindow::GetRelativeAddr(const std::shared_ptr<RWebWindow> &win) const
{
   if (fMgr != win->fMgr) {
      R__LOG_ERROR(WebGUILog()) << "Same web window manager should be used";
      return "";
   }

   std::string res("../");
   res.append(win->GetAddr());
   res.append("/");
   return res;
}

/////////////////////////////////////////////////////////////////////////
/// Set client version, used as prefix in scripts URL
/// When changed, web browser will reload all related JS files while full URL will be different
/// Default is empty value - no extra string in URL
/// Version should be string like "1.2" or "ver1.subv2" and not contain any special symbols

void RWebWindow::SetClientVersion(const std::string &vers)
{
   std::lock_guard<std::mutex> grd(fConnMutex);
   fClientVersion = vers;
}

/////////////////////////////////////////////////////////////////////////
/// Returns current client version

std::string RWebWindow::GetClientVersion() const
{
   std::lock_guard<std::mutex> grd(fConnMutex);
   return fClientVersion;
}

/////////////////////////////////////////////////////////////////////////
/// Set arbitrary JSON code, which is accessible via conn.GetUserArgs() method
/// This JSON code injected into main HTML document into JSROOT.connectWebWindow()
/// Must be called before RWebWindow::Show() method is called

void RWebWindow::SetUserArgs(const std::string &args)
{
   std::lock_guard<std::mutex> grd(fConnMutex);
   fUserArgs = args;
}

/////////////////////////////////////////////////////////////////////////
/// Returns configured user arguments for web window
/// See \ref SetUserArgs method for more details

std::string RWebWindow::GetUserArgs() const
{
   std::lock_guard<std::mutex> grd(fConnMutex);
   return fUserArgs;
}

///////////////////////////////////////////////////////////////////////////////////
/// Returns current number of active clients connections

int RWebWindow::NumConnections(bool with_pending) const
{
   std::lock_guard<std::mutex> grd(fConnMutex);
   auto sz = fConn.size();
   if (with_pending)
      sz += fPendingConn.size();
   return sz;
}

///////////////////////////////////////////////////////////////////////////////////
/// Configures recording of communication data in protocol file
/// Provided filename will be used to store JSON array with names of written files - text or binary
/// If data was send from client, "send" entry will be placed. JSON file will look like:
///    ["send","msg0.txt","send","msg1.txt","msg2.txt"]
/// If empty file name is provided, data recording will be disabled
/// Recorded data can be used in JSROOT directly to test client code without running C++ server

void RWebWindow::RecordData(const std::string &fname, const std::string &fprefix)
{
   fProtocolFileName = fname;
   fProtocolCnt = fProtocolFileName.empty() ? -1 : 0;
   fProtocolConnId = fProtocolFileName.empty() ? 0 : GetConnectionId(0);
   fProtocolPrefix = fprefix;
   fProtocol = "[]"; // empty array
}

///////////////////////////////////////////////////////////////////////////////////
/// Returns connection for specified connection number
/// Only active connections are returned - where clients confirms connection
/// Total number of connections can be retrieved with NumConnections() method

unsigned RWebWindow::GetConnectionId(int num) const
{
   std::lock_guard<std::mutex> grd(fConnMutex);
   return ((num >= 0) && (num < (int)fConn.size()) && fConn[num]->fActive) ? fConn[num]->fConnId : 0;
}

///////////////////////////////////////////////////////////////////////////////////
/// returns true if specified connection id exists
/// connid is connection (0 - any)
/// if only_active==false, also inactive connections check or connections which should appear

bool RWebWindow::HasConnection(unsigned connid, bool only_active) const
{
   std::lock_guard<std::mutex> grd(fConnMutex);

   for (auto &conn : fConn) {
      if (connid && (conn->fConnId != connid))
         continue;
      if (conn->fActive || !only_active)
         return true;
   }

   if (!only_active)
      for (auto &conn : fPendingConn) {
         if (!connid || (conn->fConnId == connid))
            return true;
      }

   return false;
}

///////////////////////////////////////////////////////////////////////////////////
/// Closes all connection to clients
/// Normally leads to closing of all correspondent browser windows
/// Some browsers (like firefox) do not allow by default to close window

void RWebWindow::CloseConnections()
{
   SubmitData(0, true, "CLOSE", 0);
}

///////////////////////////////////////////////////////////////////////////////////
/// Close specified connection
/// Connection id usually appears in the correspondent call-backs

void RWebWindow::CloseConnection(unsigned connid)
{
   if (connid)
      SubmitData(connid, true, "CLOSE", 0);
}

///////////////////////////////////////////////////////////////////////////////////
/// returns connection (or all active connections)

RWebWindow::ConnectionsList_t RWebWindow::GetConnections(unsigned connid, bool only_active) const
{
   ConnectionsList_t arr;

   {
      std::lock_guard<std::mutex> grd(fConnMutex);

      for (auto &conn : fConn) {
         if ((conn->fActive || !only_active) && (!connid || (conn->fConnId == connid)))
            arr.push_back(conn);
      }

      if (!only_active)
         for (auto &conn : fPendingConn)
            if (!connid || (conn->fConnId == connid))
               arr.push_back(conn);
   }

   return arr;
}

///////////////////////////////////////////////////////////////////////////////////
/// returns true if sending via specified connection can be performed
/// if direct==true, checks if direct sending (without queuing) is possible
/// if connid==0, all existing connections are checked

bool RWebWindow::CanSend(unsigned connid, bool direct) const
{
   auto arr = GetConnections(connid, direct); // for direct sending connection has to be active

   auto maxqlen = GetMaxQueueLength();

   for (auto &conn : arr) {

      std::lock_guard<std::mutex> grd(conn->fMutex);

      if (direct && (!conn->fQueue.empty() || (conn->fSendCredits == 0) || conn->fDoingSend))
         return false;

      if (conn->fQueue.size() >= maxqlen)
         return false;
   }

   return true;
}

///////////////////////////////////////////////////////////////////////////////////
/// returns send queue length for specified connection
/// if connid==0, maximal value for all connections is returned
/// If wrong connection is specified, -1 is return

int RWebWindow::GetSendQueueLength(unsigned connid) const
{
   int maxq = -1;

   for (auto &conn : GetConnections(connid)) {
      std::lock_guard<std::mutex> grd(conn->fMutex);
      int len = conn->fQueue.size();
      if (len > maxq) maxq = len;
   }

   return maxq;
}


///////////////////////////////////////////////////////////////////////////////////
/// Internal method to send data
/// Allows to specify channel. chid==1 is normal communication, chid==0 for internal with higher priority
/// If connid==0, data will be send to all connections

void RWebWindow::SubmitData(unsigned connid, bool txt, std::string &&data, int chid)
{
   if (fMaster)
      return fMaster->SubmitData(fMasterConnId, txt, std::move(data), fMasterChannel);

   auto arr = GetConnections(connid);
   auto cnt = arr.size();
   auto maxqlen = GetMaxQueueLength();

   timestamp_t stamp = std::chrono::system_clock::now();

   for (auto &conn : arr) {

      if (fProtocolCnt >= 0)
         if (!fProtocolConnId || (conn->fConnId == fProtocolConnId)) {
            fProtocolConnId = conn->fConnId; // remember connection
            std::string fname = fProtocolPrefix;
            fname.append("msg");
            fname.append(std::to_string(fProtocolCnt++));
            fname.append(txt ? ".txt" : ".bin");

            std::ofstream ofs(fname);
            ofs.write(data.c_str(), data.length());
            ofs.close();

            if (fProtocol.length() > 2)
               fProtocol.insert(fProtocol.length() - 1, ",");
            fProtocol.insert(fProtocol.length() - 1, "\""s + fname + "\""s);

            std::ofstream pfs(fProtocolFileName);
            pfs.write(fProtocol.c_str(), fProtocol.length());
            pfs.close();
         }

      conn->fSendStamp = stamp;

      std::lock_guard<std::mutex> grd(conn->fMutex);

      if (conn->fQueue.size() < maxqlen) {
         if (--cnt)
            conn->fQueue.emplace(chid, txt, std::string(data)); // make copy
         else
            conn->fQueue.emplace(chid, txt, std::move(data));  // move content
      } else {
         R__LOG_ERROR(WebGUILog()) << "Maximum queue length achieved";
      }
   }

   CheckDataToSend();
}

///////////////////////////////////////////////////////////////////////////////////
/// Sends data to specified connection
/// If connid==0, data will be send to all connections

void RWebWindow::Send(unsigned connid, const std::string &data)
{
   SubmitData(connid, true, std::string(data), 1);
}

///////////////////////////////////////////////////////////////////////////////////
/// Send binary data to specified connection
/// If connid==0, data will be sent to all connections

void RWebWindow::SendBinary(unsigned connid, std::string &&data)
{
   SubmitData(connid, false, std::move(data), 1);
}

///////////////////////////////////////////////////////////////////////////////////
/// Send binary data to specified connection
/// If connid==0, data will be sent to all connections

void RWebWindow::SendBinary(unsigned connid, const void *data, std::size_t len)
{
   std::string buf;
   buf.resize(len);
   std::copy((const char *)data, (const char *)data + len, buf.begin());
   SubmitData(connid, false, std::move(buf), 1);
}

///////////////////////////////////////////////////////////////////////////////////
/// Assign thread id which has to be used for callbacks

void RWebWindow::AssignCallbackThreadId()
{
   fCallbacksThrdIdSet = true;
   fCallbacksThrdId = std::this_thread::get_id();
   if (!RWebWindowsManager::IsMainThrd()) {
      fProcessMT = true;
   } else if (fMgr->IsUseHttpThread()) {
      // special thread is used by the manager, but main thread used for the canvas - not supported
      R__LOG_ERROR(WebGUILog()) << "create web window from main thread when THttpServer created with special thread - not supported";
   }
}

/////////////////////////////////////////////////////////////////////////////////
/// Set call-back function for data, received from the clients via websocket
///
/// Function should have signature like void func(unsigned connid, const std::string &data)
/// First argument identifies connection (unique for each window), second argument is received data
///
/// At the moment when callback is assigned, RWebWindow working thread is detected.
/// If called not from main application thread, RWebWindow::Run() function must be regularly called from that thread.
///
/// Most simple way to assign call-back - use of c++11 lambdas like:
/// ~~~ {.cpp}
/// auto win = RWebWindow::Create();
/// win->SetDefaultPage("file:./page.htm");
/// win->SetDataCallBack(
///          [](unsigned connid, const std::string &data) {
///                  printf("Conn:%u data:%s\n", connid, data.c_str());
///           }
///       );
/// win->Show();
/// ~~~

void RWebWindow::SetDataCallBack(WebWindowDataCallback_t func)
{
   AssignCallbackThreadId();
   fDataCallback = func;
}

/////////////////////////////////////////////////////////////////////////////////
/// Set call-back function for new connection

void RWebWindow::SetConnectCallBack(WebWindowConnectCallback_t func)
{
   AssignCallbackThreadId();
   fConnCallback = func;
}

/////////////////////////////////////////////////////////////////////////////////
/// Set call-back function for disconnecting

void RWebWindow::SetDisconnectCallBack(WebWindowConnectCallback_t func)
{
   AssignCallbackThreadId();
   fDisconnCallback = func;
}

/////////////////////////////////////////////////////////////////////////////////
/// Set call-backs function for connect, data and disconnect events

void RWebWindow::SetCallBacks(WebWindowConnectCallback_t conn, WebWindowDataCallback_t data, WebWindowConnectCallback_t disconn)
{
   AssignCallbackThreadId();
   fConnCallback = conn;
   fDataCallback = data;
   fDisconnCallback = disconn;
}

/////////////////////////////////////////////////////////////////////////////////
/// Waits until provided check function or lambdas returns non-zero value
/// Check function has following signature: int func(double spent_tm)
/// Waiting will be continued, if function returns zero.
/// Parameter spent_tm is time in seconds, which already spent inside the function
/// First non-zero value breaks loop and result is returned.
/// Runs application mainloop and short sleeps in-between

int RWebWindow::WaitFor(WebWindowWaitFunc_t check)
{
   return fMgr->WaitFor(*this, check);
}

/////////////////////////////////////////////////////////////////////////////////
/// Waits until provided check function or lambdas returns non-zero value
/// Check function has following signature: int func(double spent_tm)
/// Waiting will be continued, if function returns zero.
/// Parameter spent_tm in lambda is time in seconds, which already spent inside the function
/// First non-zero value breaks waiting loop and result is returned (or 0 if time is expired).
/// Runs application mainloop and short sleeps in-between
/// WebGui.OperationTmout rootrc parameter defines waiting time in seconds

int RWebWindow::WaitForTimed(WebWindowWaitFunc_t check)
{
   return fMgr->WaitFor(*this, check, true, GetOperationTmout());
}

/////////////////////////////////////////////////////////////////////////////////
/// Waits until provided check function or lambdas returns non-zero value
/// Check function has following signature: int func(double spent_tm)
/// Waiting will be continued, if function returns zero.
/// Parameter spent_tm in lambda is time in seconds, which already spent inside the function
/// First non-zero value breaks waiting loop and result is returned (or 0 if time is expired).
/// Runs application mainloop and short sleeps in-between
/// duration (in seconds) defines waiting time

int RWebWindow::WaitForTimed(WebWindowWaitFunc_t check, double duration)
{
   return fMgr->WaitFor(*this, check, true, duration);
}


/////////////////////////////////////////////////////////////////////////////////
/// Run window functionality for specified time
/// If no action can be performed - just sleep specified time

void RWebWindow::Run(double tm)
{
   if (!fCallbacksThrdIdSet || (fCallbacksThrdId != std::this_thread::get_id())) {
      R__LOG_WARNING(WebGUILog()) << "Change thread id where RWebWindow is executed";
      fCallbacksThrdIdSet = true;
      fCallbacksThrdId = std::this_thread::get_id();
   }

   if (tm <= 0) {
      Sync();
   } else {
      WaitForTimed([](double) { return 0; }, tm);
   }
}


/////////////////////////////////////////////////////////////////////////////////
/// Add embed window

unsigned RWebWindow::AddEmbedWindow(std::shared_ptr<RWebWindow> window, int channel)
{
   if (channel < 2)
      return 0;

   auto arr = GetConnections(0, true);
   if (arr.size() == 0)
      return 0;

   // check if channel already occupied
   if (arr[0]->fEmbed.find(channel) != arr[0]->fEmbed.end())
      return 0;

   arr[0]->fEmbed[channel] = window;

   return arr[0]->fConnId;
}

/////////////////////////////////////////////////////////////////////////////////
/// Remove RWebWindow associated with the channel

void RWebWindow::RemoveEmbedWindow(unsigned connid, int channel)
{
   auto arr = GetConnections(connid);

   for (auto &conn : arr) {
      auto iter = conn->fEmbed.find(channel);
      if (iter != conn->fEmbed.end())
         conn->fEmbed.erase(iter);
   }
}


/////////////////////////////////////////////////////////////////////////////////
/// Create new RWebWindow
/// Using default RWebWindowsManager

std::shared_ptr<RWebWindow> RWebWindow::Create()
{
   return RWebWindowsManager::Instance()->CreateWindow();
}

/////////////////////////////////////////////////////////////////////////////////
/// Terminate ROOT session
/// Tries to correctly close THttpServer, associated with RWebWindowsManager
/// After that exit from process

void RWebWindow::TerminateROOT()
{

   // workaround to release all connection-specific handles as soon as possible
   // required to work with QWebEngine
   // once problem solved, can be removed here
   ConnectionsList_t arr1, arr2;

   {
      std::lock_guard<std::mutex> grd(fConnMutex);
      std::swap(arr1, fConn);
      std::swap(arr2, fPendingConn);
   }

   fMgr->Terminate();
}

/////////////////////////////////////////////////////////////////////////////////
/// Static method to show web window
/// Has to be used instead of RWebWindow::Show() when window potentially can be embed into other windows
/// Soon RWebWindow::Show() method will be done protected

unsigned RWebWindow::ShowWindow(std::shared_ptr<RWebWindow> window, const RWebDisplayArgs &args)
{
   if (!window)
      return 0;

   if (args.GetBrowserKind() == RWebDisplayArgs::kEmbedded) {
      unsigned connid = args.fMaster ? args.fMaster->AddEmbedWindow(window, args.fMasterChannel) : 0;

      if (connid > 0) {
         window->fMaster = args.fMaster;
         window->fMasterConnId = connid;
         window->fMasterChannel = args.fMasterChannel;

         // inform client that connection is established and window initialized
         args.fMaster->SubmitData(connid, true, "EMBED_DONE"s, args.fMasterChannel);

         // provide call back for window itself that connection is ready
         window->ProvideQueueEntry(connid, kind_Connect, ""s);
      }

      return connid;
   }

   return window->Show(args);
}

