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

#include <ROOT/RWebWindowsManager.hxx>

#include <ROOT/RLogger.hxx>
#include <ROOT/RWebDisplayArgs.hxx>
#include <ROOT/RWebDisplayHandle.hxx>

#include "RWebWindowWSHandler.hxx"

#include "THttpServer.h"

#include "TSystem.h"
#include "TRandom.h"
#include "TString.h"
#include "TApplication.h"
#include "TTimer.h"
#include "TROOT.h"
#include "TEnv.h"

#include <thread>
#include <chrono>

using namespace ROOT::Experimental;

/** \class ROOT::Experimental::RWebWindowsManager
\ingroup webdisplay

Central instance to create and show web-based windows like Canvas or FitPanel.

Manager responsible to creating THttpServer instance, which is used for RWebWindow's
communication with clients.

Method RWebWindowsManager::Show() used to show window in specified location.
*/

//////////////////////////////////////////////////////////////////////////////////////////
/// Returns default window manager
/// Used to display all standard ROOT elements like TCanvas or TFitPanel

std::shared_ptr<RWebWindowsManager> &RWebWindowsManager::Instance()
{
   static std::shared_ptr<RWebWindowsManager> sInstance = std::make_shared<RWebWindowsManager>();
   return sInstance;
}

//////////////////////////////////////////////////////////////////
/// This thread id used to identify main application thread, where ROOT event processing runs
/// To inject code in that thread, one should use TTimer (like THttpServer does)
/// In other threads special run methods have to be invoked like RWebWindow::Run()
///
/// TODO: probably detection of main thread should be delivered by central ROOT instances like gApplication or gROOT
/// Main thread can only make sense if special processing runs there and one can inject own functionality there

static std::thread::id gWebWinMainThrd = std::this_thread::get_id();

//////////////////////////////////////////////////////////////////////////////////////////
/// Returns true when called from main process
/// Main process recognized at the moment when library is loaded

bool RWebWindowsManager::IsMainThrd()
{
   return std::this_thread::get_id() == gWebWinMainThrd;
}

//////////////////////////////////////////////////////////////////////////////////////////
/// window manager constructor
/// Required here for correct usage of unique_ptr<THttpServer>

RWebWindowsManager::RWebWindowsManager() = default;

//////////////////////////////////////////////////////////////////////////////////////////
/// window manager destructor
/// Required here for correct usage of unique_ptr<THttpServer>

RWebWindowsManager::~RWebWindowsManager()
{
   if (gApplication && fServer && !fServer->IsTerminated()) {
      gApplication->Disconnect("Terminate(Int_t)", fServer.get(), "SetTerminate()");
      fServer->SetTerminate();
   }
}

//////////////////////////////////////////////////////////////////////////////////////////
/// Creates http server, if required - with real http engine (civetweb)
/// One could configure concrete HTTP port, which should be used for the server,
/// provide following entry in rootrc file:
///
///      WebGui.HttpPort: 8088
///
/// or specify range of http ports, which can be used:
///
///      WebGui.HttpPortMin: 8800
///      WebGui.HttpPortMax: 9800
///
/// By default range [8800..9800] is used
///
/// One also can bind HTTP server socket to loopback address,
/// In that case only connection from localhost will be available:
///
///      WebGui.HttpLoopback: yes
///
/// Or one could specify hostname which should be used for binding of server socket
///
///      WebGui.HttpBind: hostname | ipaddress
///
/// To use secured protocol, following parameter should be specified
///
///      WebGui.UseHttps: yes
///      WebGui.ServerCert: sertificate_filename.pem
///
/// One also can configure usage of special thread of processing of http server requests
///
///      WebGui.HttpThrd: no
///
/// Extra threads can be used to send data to different clients via websocket (default no)
///
///      WebGui.SenderThrds: no
///
/// If required, one could change websocket timeouts (default is 10000 ms)
///
///      WebGui.HttpWSTmout: 10000
///
/// By default, THttpServer created in restricted mode which only allows websocket handlers
/// and processes only very few other related http requests. For security reasons such mode
/// should be always enabled. Only if it is really necessary to process all other kinds
/// of HTTP requests, one could specify 0 for following parameter (default 1):
///
///      WebGui.WSOnly: 1
///
/// Following parameter controls browser max-age caching parameter for files (default 3600)
///
///      WebGui.HttpMaxAge: 3600
///
/// One also can configure usage of FastCGI server for web windows:
///
///      WebGui.FastCgiPort: 4000
///      WebGui.FastCgiThreads: 10
///
/// To be able start web browser for such windows, one can provide real URL of the
/// web server which will connect with that FastCGI instance:
///
///      WebGui.FastCgiServer: https://your_apache_server.com/root_cgi_path
///

bool RWebWindowsManager::CreateServer(bool with_http)
{
   // explicitly protect server creation
   std::lock_guard<std::recursive_mutex> grd(fMutex);

   if (!fServer) {

      fServer = std::make_unique<THttpServer>("basic_sniffer");

      const char *serv_thrd = gEnv->GetValue("WebGui.HttpThrd", "");
      if (serv_thrd && strstr(serv_thrd, "yes"))
         fUseHttpThrd = true;
      else if (serv_thrd && strstr(serv_thrd, "no"))
         fUseHttpThrd = false;

      const char *send_thrds = gEnv->GetValue("WebGui.SenderThrds", "");
      if (send_thrds && *send_thrds) {
         if (strstr(send_thrds, "yes"))
            fUseSenderThreads = true;
         else if (strstr(send_thrds, "no"))
            fUseSenderThreads = false;
         else
            R__LOG_ERROR(WebGUILog()) << "WebGui.SenderThrds has to be yes or no";
      }

      if (IsUseHttpThread())
         fServer->CreateServerThread();

      if (gApplication)
         gApplication->Connect("Terminate(Int_t)", "THttpServer", fServer.get(), "SetTerminate()");

      Int_t wsonly = gEnv->GetValue("WebGui.WSOnly", 1);
      fServer->SetWSOnly(wsonly != 0);

      // this is location where all ROOT UI5 sources are collected
      // normally it is $ROOTSYS/ui5 or <prefix>/ui5 location
      TString ui5dir = gSystem->Getenv("ROOTUI5SYS");
      if (ui5dir.Length() == 0)
         ui5dir = gEnv->GetValue("WebGui.RootUi5Path","");

      if (ui5dir.Length() == 0)
         ui5dir.Form("%s/ui5", TROOT::GetDataDir().Data());

      if (gSystem->ExpandPathName(ui5dir)) {
         R__LOG_ERROR(WebGUILog()) << "Path to ROOT ui5 sources " << ui5dir << " not found, set ROOTUI5SYS correctly";
         ui5dir = ".";
      }

      fServer->AddLocation("rootui5sys/", ui5dir.Data());
   }

   if (!with_http || fServer->IsAnyEngine())
      return true;

   int http_port = gEnv->GetValue("WebGui.HttpPort", 0);
   int http_min = gEnv->GetValue("WebGui.HttpPortMin", 8800);
   int http_max = gEnv->GetValue("WebGui.HttpPortMax", 9800);
   int http_wstmout = gEnv->GetValue("WebGui.HttpWSTmout", 10000);
   int http_maxage = gEnv->GetValue("WebGui.HttpMaxAge", -1);
   int fcgi_port = gEnv->GetValue("WebGui.FastCgiPort", 0);
   int fcgi_thrds = gEnv->GetValue("WebGui.FastCgiThreads", 10);
   const char *fcgi_serv = gEnv->GetValue("WebGui.FastCgiServer", "");
   fLaunchTmout = gEnv->GetValue("WebGui.LaunchTmout", 30.);
   const char *http_loopback = gEnv->GetValue("WebGui.HttpLoopback", "no");
   const char *http_bind = gEnv->GetValue("WebGui.HttpBind", "");
   const char *http_ssl = gEnv->GetValue("WebGui.UseHttps", "no");
   const char *ssl_cert = gEnv->GetValue("WebGui.ServerCert", "rootserver.pem");

   bool assign_loopback = http_loopback && strstr(http_loopback, "yes");
   bool use_secure = http_ssl && strstr(http_ssl, "yes");
   int ntry = 100;

   if ((http_port < 0) && (fcgi_port <= 0)) {
      R__LOG_ERROR(WebGUILog()) << "Not allowed to create HTTP server, check WebGui.HttpPort variable";
      return false;
   }

   if (http_port < 0) {
      ntry = 0;
   } else {

      if (http_port == 0)
         gRandom->SetSeed(0);

      if (http_max - http_min < ntry)
         ntry = http_max - http_min;
   }

   if (fcgi_port > 0)
      ntry++;

   while (ntry-- >= 0) {
      if ((http_port == 0) && (fcgi_port <= 0)) {
         if ((http_min <= 0) || (http_max <= http_min)) {
            R__LOG_ERROR(WebGUILog()) << "Wrong HTTP range configuration, check WebGui.HttpPortMin/Max variables";
            return false;
         }

         http_port = (int)(http_min + (http_max - http_min) * gRandom->Rndm(1));
      }

      TString engine, url;

      if (fcgi_port > 0) {
         engine.Form("fastcgi:%d?thrds=%d", fcgi_port, fcgi_thrds);
         if (!fServer->CreateEngine(engine)) return false;
         if (fcgi_serv && (strlen(fcgi_serv) > 0)) fAddr = fcgi_serv;
         if (http_port < 0) return true;
         fcgi_port = 0;
      } else if (http_port > 0) {
         url = use_secure ? "https://" : "http://";
         engine.Form("%s:%d?websocket_timeout=%d", (use_secure ? "https" : "http"), http_port, http_wstmout);
         if (assign_loopback) {
            engine.Append("&loopback");
            url.Append("localhost");
         } else if (http_bind && (strlen(http_bind) > 0)) {
            engine.Append("&bind=");
            engine.Append(http_bind);
            url.Append(http_bind);
         } else {
            url.Append("localhost");
         }

         if (http_maxage >= 0)
            engine.Append(TString::Format("&max_age=%d", http_maxage));

         if (use_secure) {
            engine.Append("&ssl_cert=");
            engine.Append(ssl_cert);
         }

         if (fServer->CreateEngine(engine)) {
            fAddr = url.Data();
            fAddr.append(":");
            fAddr.append(std::to_string(http_port));
            return true;
         }
         http_port = 0;
      }
   }

   return false;
}

//////////////////////////////////////////////////////////////////////////////////////////
/// Creates new window
/// To show window, RWebWindow::Show() have to be called

std::shared_ptr<RWebWindow> RWebWindowsManager::CreateWindow()
{

   // we book manager mutex for a longer operation, locked again in server creation
   std::lock_guard<std::recursive_mutex> grd(fMutex);

   if (!CreateServer()) {
      R__LOG_ERROR(WebGUILog()) << "Cannot create server when creating window";
      return nullptr;
   }

   std::shared_ptr<RWebWindow> win = std::make_shared<RWebWindow>();

   if (!win) {
      R__LOG_ERROR(WebGUILog()) << "Fail to create RWebWindow instance";
      return nullptr;
   }

   double dflt_tmout = gEnv->GetValue("WebGui.OperationTmout", 50.);

   auto wshandler = win->CreateWSHandler(Instance(), ++fIdCnt, dflt_tmout);

   if (gEnv->GetValue("WebGui.RecordData", 0) > 0) {
      std::string fname, prefix;
      if (fIdCnt > 1) {
         prefix = std::string("f") + std::to_string(fIdCnt) + "_";
         fname = std::string("protcol") + std::to_string(fIdCnt) + ".json";
      } else {
         fname = "protocol.json";
      }
      win->RecordData(fname, prefix);
   }

   const char *token = gEnv->GetValue("WebGui.ConnToken", "");
   if (token && *token)
      win->SetConnToken(token);

   fServer->RegisterWS(wshandler);

   return win;
}

//////////////////////////////////////////////////////////////////////////////////////////
/// Release all references to specified window
/// Called from RWebWindow destructor

void RWebWindowsManager::Unregister(RWebWindow &win)
{
   if (win.fWSHandler)
      fServer->UnregisterWS(win.fWSHandler);
}

//////////////////////////////////////////////////////////////////////////
/// Provide URL address to access specified window from inside or from remote

std::string RWebWindowsManager::GetUrl(const RWebWindow &win, bool remote)
{
   if (!fServer) {
      R__LOG_ERROR(WebGUILog()) << "Server instance not exists when requesting window URL";
      return "";
   }

   std::string addr = "/";

   addr.append(win.fWSHandler->GetName());

   addr.append("/");

   if (remote) {
      if (!CreateServer(true)) {
         R__LOG_ERROR(WebGUILog()) << "Fail to start real HTTP server when requesting URL";
         return "";
      }

      addr = fAddr + addr;
   }

   return addr;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// Show web window in specified location.
///
/// \param batch_mode indicates that browser will run in headless mode
/// \param user_args specifies where and how display web window
///
/// As display args one can use string like "firefox" or "chrome" - these are two main supported web browsers.
/// See RWebDisplayArgs::SetBrowserKind() for all available options. Default value for the browser can be configured
/// when starting root with --web argument like: "root --web=chrome"
///
///  If allowed, same window can be displayed several times (like for TCanvas)
///
///  Following parameters can be configured in rootrc file:
///
///   WebGui.Chrome:  full path to Google Chrome executable
///   WebGui.ChromeBatch: command to start chrome in batch
///   WebGui.ChromeInteractive: command to start chrome in interactive mode
///   WebGui.Firefox: full path to Mozialla Firefox executable
///   WebGui.FirefoxBatch: command to start Firefox in batch mode
///   WebGui.FirefoxInteractive: command to start Firefox in interactive mode
///   WebGui.FirefoxProfile: name of Firefox profile to use
///   WebGui.FirefoxProfilePath: file path to Firefox profile
///   WebGui.FirefoxRandomProfile: usage of random Firefox profile -1 never, 0 - only for batch mode (dflt), 1 - always
///   WebGui.LaunchTmout: time required to start process in seconds (default 30 s)
///   WebGui.OperationTmout: time required to perform WebWindow operation like execute command or update drawings
///   WebGui.RecordData: if specified enables data recording for each web window 0 - off, 1 - on
///   WebGui.JsonComp: compression factor for JSON conversion, if not specified - each widget uses own default values
///   WebGui.ForceHttp: 0 - off (default), 1 - always create real http server to run web window
///   WebGui.Console: -1 - output only console.error(), 0 - add console.warn(), 1  - add console.log() output
///   WebGui.ConnCredits: 10 - number of packets which can be send by server or client without acknowledge from receiving side
///   WebGui.openui5src:   alternative location for openui5 like https://openui5.hana.ondemand.com/
///   WebGui.openui5libs:  list of pre-loaded ui5 libs like sap.m, sap.ui.layout, sap.ui.unified
///   WebGui.openui5theme:  openui5 theme like sap_belize (default) or sap_fiori_3
///
///   HTTP-server related parameters documented in RWebWindowsManager::CreateServer() method

unsigned RWebWindowsManager::ShowWindow(RWebWindow &win, bool batch_mode, const RWebDisplayArgs &user_args)
{
   // silently ignore regular Show() calls in batch mode
   if (!batch_mode && gROOT->IsWebDisplayBatch())
      return 0;

   // for embedded window no any browser need to be started
   if (user_args.GetBrowserKind() == RWebDisplayArgs::kEmbedded)
      return 0;

   // place here while involves conn mutex
   auto token = win.GetConnToken();

   // we book manager mutex for a longer operation,
   std::lock_guard<std::recursive_mutex> grd(fMutex);

   if (!fServer) {
      R__LOG_ERROR(WebGUILog()) << "Server instance not exists to show window";
      return 0;
   }

   std::string key;
   int ntry = 100000;

   do {
      key = std::to_string(gRandom->Integer(0x100000));
   } while ((--ntry > 0) && win.HasKey(key));
   if (ntry == 0) {
      R__LOG_ERROR(WebGUILog()) << "Fail to create unique key for the window";
      return 0;
   }

   RWebDisplayArgs args(user_args);

   if (batch_mode && !args.IsSupportHeadless()) {
      R__LOG_ERROR(WebGUILog()) << "Cannot use batch mode with " << args.GetBrowserName();
      return 0;
   }

   args.SetHeadless(batch_mode);
   if (args.GetWidth() <= 0) args.SetWidth(win.GetWidth());
   if (args.GetHeight() <= 0) args.SetHeight(win.GetHeight());

   bool normal_http = !args.IsLocalDisplay();
   if (!normal_http && (gEnv->GetValue("WebGui.ForceHttp",0) == 1))
      normal_http = true;

   std::string url = GetUrl(win, normal_http);
   if (url.empty()) {
      R__LOG_ERROR(WebGUILog()) << "Cannot create URL for the window";
      return 0;
   }
   if (normal_http && fAddr.empty()) {
      R__LOG_WARNING(WebGUILog()) << "Full URL cannot be produced for window " << url << " to start web browser";
      return 0;
   }

   args.SetUrl(url);

   args.AppendUrlOpt(std::string("key=") + key);
   if (batch_mode) args.AppendUrlOpt("batch_mode");
   if (!token.empty())
      args.AppendUrlOpt(std::string("token=") + token);

   if (!normal_http)
      args.SetHttpServer(GetServer());

   auto handle = RWebDisplayHandle::Display(args);

   if (!handle) {
      R__LOG_ERROR(WebGUILog()) << "Cannot display window in " << args.GetBrowserName();
      return 0;
   }

   return win.AddDisplayHandle(batch_mode, key, handle);
}

//////////////////////////////////////////////////////////////////////////
/// Waits until provided check function or lambdas returns non-zero value
/// Regularly calls WebWindow::Sync() method to let run event loop
/// If call from the main thread, runs system events processing
/// Check function has following signature: int func(double spent_tm)
/// Parameter spent_tm is time in seconds, which already spent inside function
/// Waiting will be continued, if function returns zero.
/// First non-zero value breaks waiting loop and result is returned (or 0 if time is expired).
/// If parameter timed is true, timelimit (in seconds) defines how long to wait

int RWebWindowsManager::WaitFor(RWebWindow &win, WebWindowWaitFunc_t check, bool timed, double timelimit)
{
   int res = 0;
   int cnt = 0;
   double spent = 0;

   auto start = std::chrono::high_resolution_clock::now();

   win.Sync(); // in any case call sync once to ensure

   while ((res = check(spent)) == 0) {

      if (IsMainThrd())
         gSystem->ProcessEvents();

      win.Sync();

      std::this_thread::sleep_for(std::chrono::milliseconds(1));

      std::chrono::duration<double, std::milli> elapsed = std::chrono::high_resolution_clock::now() - start;

      spent = elapsed.count() * 1e-3; // use ms precision

      if (timed && (spent > timelimit))
         return -3;

      cnt++;
   }

   return res;
}

//////////////////////////////////////////////////////////////////////////
/// Terminate http server and ROOT application

void RWebWindowsManager::Terminate()
{
   if (fServer)
      fServer->SetTerminate();

   if (gApplication)
      TTimer::SingleShot(100, "TApplication",  gApplication, "Terminate()");
}
