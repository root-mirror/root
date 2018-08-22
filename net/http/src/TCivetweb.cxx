// $Id$
// Author: Sergey Linev   21/12/2013

/*************************************************************************
 * Copyright (C) 1995-2013, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TCivetweb.h"

#include "../civetweb/civetweb.h"

#include <stdlib.h>
#include <string.h>

#include "THttpServer.h"
#include "THttpWSEngine.h"
#include "TUrl.h"

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TCivetwebWSEngine                                                    //
//                                                                      //
// Implementation of THttpWSEngine for Civetweb                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

class TCivetwebWSEngine : public THttpWSEngine {
protected:
   struct mg_connection *fWSconn;

   /// True websocket requires extra thread to parallelize sending
   virtual Bool_t SupportSendThrd() const { return kTRUE; }

public:
   TCivetwebWSEngine(struct mg_connection *conn) : THttpWSEngine(), fWSconn(conn) {}

   virtual ~TCivetwebWSEngine()
   {
      TCivetwebWSEngine::ClearHandle(kTRUE);
   }

   virtual UInt_t GetId() const { return TString::Hash((void *)&fWSconn, sizeof(void *)); }

   virtual void ClearHandle(Bool_t terminate) override
   {
      if (fWSconn && terminate)
         mg_websocket_write(fWSconn, MG_WEBSOCKET_OPCODE_CONNECTION_CLOSE, nullptr, 0);
      fWSconn = nullptr;
   }

   virtual void Send(const void *buf, int len)
   {
      if (fWSconn)
         mg_websocket_write(fWSconn, MG_WEBSOCKET_OPCODE_BINARY, (const char *)buf, len);
   }

   /////////////////////////////////////////////////////////
   /// Special method to send binary data with text header
   /// For normal websocket it is two separated operation, for other engines could be combined together,
   /// but emulates as two messages on client side
   virtual void SendHeader(const char *hdr, const void *buf, int len)
   {
      if (fWSconn) {
         mg_websocket_write(fWSconn, MG_WEBSOCKET_OPCODE_TEXT, hdr, strlen(hdr));
         mg_websocket_write(fWSconn, MG_WEBSOCKET_OPCODE_BINARY, (const char *)buf, len);
      }
   }

   virtual void SendCharStar(const char *str)
   {
      if (fWSconn)
         mg_websocket_write(fWSconn, MG_WEBSOCKET_OPCODE_TEXT, str, strlen(str));
   }
};

//////////////////////////////////////////////////////////////////////////

int websocket_connect_handler(const struct mg_connection *conn, void *)
{
   const struct mg_request_info *request_info = mg_get_request_info(conn);
   if (!request_info)
      return 1;

   TCivetweb *engine = (TCivetweb *)request_info->user_data;
   if (!engine || engine->IsTerminating())
      return 1;
   THttpServer *serv = engine->GetServer();
   if (!serv)
      return 1;

   auto arg = std::make_shared<THttpCallArg>();
   arg->SetPathAndFileName(request_info->local_uri); // path and file name
   arg->SetQuery(request_info->query_string);        // query arguments
   arg->SetWSId(TString::Hash((void *)&conn, sizeof(void *)));
   arg->SetMethod("WS_CONNECT");

   Bool_t execres = serv->ExecuteWS(arg, kTRUE, kTRUE);

   return execres && !arg->Is404() ? 0 : 1;
}

//////////////////////////////////////////////////////////////////////////

void websocket_ready_handler(struct mg_connection *conn, void *)
{
   const struct mg_request_info *request_info = mg_get_request_info(conn);

   TCivetweb *engine = (TCivetweb *)request_info->user_data;
   if (!engine || engine->IsTerminating())
      return;
   THttpServer *serv = engine->GetServer();
   if (!serv)
      return;

   auto arg = std::make_shared<THttpCallArg>();
   arg->SetPathAndFileName(request_info->local_uri); // path and file name
   arg->SetQuery(request_info->query_string);        // query arguments
   arg->SetMethod("WS_READY");

   // delegate ownership to the arg, id will be automatically set
   arg->CreateWSEngine<TCivetwebWSEngine>(conn);

   serv->ExecuteWS(arg, kTRUE, kTRUE);
}

//////////////////////////////////////////////////////////////////////////

int websocket_data_handler(struct mg_connection *conn, int, char *data, size_t len, void *)
{
   const struct mg_request_info *request_info = mg_get_request_info(conn);

   // do not handle empty data
   if (len == 0)
      return 1;

   TCivetweb *engine = (TCivetweb *)request_info->user_data;
   if (!engine || engine->IsTerminating())
      return 1;
   THttpServer *serv = engine->GetServer();
   if (!serv)
      return 1;

   // seems to be, appears when connection is broken
   if ((len == 2) && ((int)data[0] == 3) && ((int)data[1] == -23))
      return 0;

   auto arg = std::make_shared<THttpCallArg>();
   arg->SetPathAndFileName(request_info->local_uri); // path and file name
   arg->SetQuery(request_info->query_string);        // query arguments
   arg->SetWSId(TString::Hash((void *)&conn, sizeof(void *)));
   arg->SetMethod("WS_DATA");

   arg->SetPostData(std::string(data,len));

   serv->ExecuteWS(arg, kTRUE, kTRUE);

   return 1;
}

//////////////////////////////////////////////////////////////////////////

void websocket_close_handler(const struct mg_connection *conn, void *)
{
   const struct mg_request_info *request_info = mg_get_request_info(conn);

   TCivetweb *engine = (TCivetweb *)request_info->user_data;
   if (!engine || engine->IsTerminating())
      return;
   THttpServer *serv = engine->GetServer();
   if (!serv)
      return;

   auto arg = std::make_shared<THttpCallArg>();
   arg->SetPathAndFileName(request_info->local_uri); // path and file name
   arg->SetQuery(request_info->query_string);        // query arguments
   arg->SetWSId(TString::Hash((void *)&conn, sizeof(void *)));
   arg->SetMethod("WS_CLOSE");

   serv->ExecuteWS(arg, kTRUE, kFALSE); // do not wait for result of execution
}

//////////////////////////////////////////////////////////////////////////

static int log_message_handler(const struct mg_connection *conn, const char *message)
{
   const struct mg_context *ctx = mg_get_context(conn);

   TCivetweb *engine = (TCivetweb *)mg_get_user_data(ctx);

   if (engine)
      return engine->ProcessLog(message);

   // provide debug output
   if ((gDebug > 0) || (strstr(message, "cannot bind to") != 0))
      fprintf(stderr, "Error in <TCivetweb::Log> %s\n", message);

   return 0;
}

//////////////////////////////////////////////////////////////////////////

static int begin_request_handler(struct mg_connection *conn, void *)
{
   const struct mg_request_info *request_info = mg_get_request_info(conn);

   TCivetweb *engine = (TCivetweb *)request_info->user_data;
   if (!engine || engine->IsTerminating())
      return 0;
   THttpServer *serv = engine->GetServer();
   if (!serv)
      return 0;

   auto arg = std::make_shared<THttpCallArg>();

   TString filename;

   Bool_t execres = kTRUE, debug = engine->IsDebugMode();

   if (!debug && serv->IsFileRequested(request_info->local_uri, filename)) {
      if ((filename.Index(".js") != kNPOS) || (filename.Index(".css") != kNPOS)) {
         std::string buf = THttpServer::ReadFileContent(filename.Data());
         if (buf.empty()) {
            arg->Set404();
         } else {
            arg->SetContentType(THttpServer::GetMimeType(filename.Data()));
            arg->SetContent(std::move(buf));
            arg->AddHeader("Cache-Control", "max-age=3600");
            arg->SetZipping();
         }
      } else {
         arg->SetFile(filename.Data());
      }
   } else {
      arg->SetPathAndFileName(request_info->local_uri); // path and file name
      arg->SetQuery(request_info->query_string);        // query arguments
      arg->SetTopName(engine->GetTopName());
      arg->SetMethod(request_info->request_method); // method like GET or POST
      if (request_info->remote_user)
         arg->SetUserName(request_info->remote_user);

      TString header;
      for (int n = 0; n < request_info->num_headers; n++)
         header.Append(
            TString::Format("%s: %s\r\n", request_info->http_headers[n].name, request_info->http_headers[n].value));
      arg->SetRequestHeader(header);

      const char *len = mg_get_header(conn, "Content-Length");
      Int_t ilen = len ? TString(len).Atoi() : 0;

      if (ilen > 0) {
         std::string buf;
         buf.resize(ilen);
         Int_t iread = mg_read(conn, (void *) buf.data(), ilen);
         if (iread == ilen)
            arg->SetPostData(std::move(buf));
      }

      if (debug) {
         TString cont;
         cont.Append("<title>Civetweb echo</title>");
         cont.Append("<h1>Civetweb echo</h1>\n");

         static int count = 0;

         cont.Append(TString::Format("Request %d:<br/>\n<pre>\n", ++count));
         cont.Append(TString::Format("  Method   : %s\n", arg->GetMethod()));
         cont.Append(TString::Format("  PathName : %s\n", arg->GetPathName()));
         cont.Append(TString::Format("  FileName : %s\n", arg->GetFileName()));
         cont.Append(TString::Format("  Query    : %s\n", arg->GetQuery()));
         cont.Append(TString::Format("  PostData : %ld\n", arg->GetPostDataLength()));
         if (arg->GetUserName())
            cont.Append(TString::Format("  User     : %s\n", arg->GetUserName()));

         cont.Append("</pre><p>\n");

         cont.Append("Environment:<br/>\n<pre>\n");
         for (int n = 0; n < request_info->num_headers; n++)
            cont.Append(
               TString::Format("  %s = %s\n", request_info->http_headers[n].name, request_info->http_headers[n].value));
         cont.Append("</pre><p>\n");

         arg->SetContentType("text/html");

         arg->SetContent(cont);

      } else {
         execres = serv->ExecuteHttp(arg);
      }
   }

   if (!execres || arg->Is404()) {
      std::string hdr = arg->FillHttpHeader("HTTP/1.1");
      mg_printf(conn, "%s", hdr.c_str());
   } else if (arg->IsFile()) {
      mg_send_file(conn, (const char *)arg->GetContent());
   } else {

      Bool_t dozip = kFALSE;
      switch (arg->GetZipping()) {
      case THttpCallArg::kNoZip: dozip = kFALSE; break;
      case THttpCallArg::kZipLarge:
         if (arg->GetContentLength() < 10000) break;
      case THttpCallArg::kZip:
         // check if request header has Accept-Encoding
         for (int n = 0; n < request_info->num_headers; n++) {
            TString name = request_info->http_headers[n].name;
            if (name.Index("Accept-Encoding", 0, TString::kIgnoreCase) != 0)
               continue;
            TString value = request_info->http_headers[n].value;
            dozip = (value.Index("gzip", 0, TString::kIgnoreCase) != kNPOS);
            break;
         }

         break;
      case THttpCallArg::kZipAlways: dozip = kTRUE; break;
      }

      if (dozip)
         arg->CompressWithGzip();

      std::string hdr = arg->FillHttpHeader("HTTP/1.1");
      mg_printf(conn, "%s", hdr.c_str());

      if (arg->GetContentLength() > 0)
         mg_write(conn, arg->GetContent(), (size_t)arg->GetContentLength());
   }

   // Returning non-zero tells civetweb that our function has replied to
   // the client, and civetweb should not send client any more data.
   return 1;
}

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TCivetweb                                                            //
//                                                                      //
// http server implementation, based on civetweb embedded server        //
// It is default kind of engine, created for THttpServer                //
// Currently v1.8 from https://github.com/civetweb/civetweb is used     //
//                                                                      //
// Following additional options can be specified:                       //
//    top=foldername - name of top folder, seen in the browser          //
//    thrds=N - use N threads to run civetweb server (default 5)        //
//    auth_file - global authentication file                            //
//    auth_domain - domain name, used for authentication                //
//                                                                      //
// Example:                                                             //
//    new THttpServer("http:8080?top=MyApp&thrds=3");                   //
//                                                                      //
// Authentication:                                                      //
//    When auth_file and auth_domain parameters are specified, access   //
//    to running http server will be possible only after user           //
//    authentication, using so-call digest method. To generate          //
//    authentication file, htdigest routine should be used:             //
//                                                                      //
//        [shell] htdigest -c .htdigest domain_name user                //
//                                                                      //
//    When creating server, parameters should be:                       //
//                                                                      //
//       new THttpServer("http:8080?auth_file=.htdigets&auth_domain=domain_name");  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// constructor

TCivetweb::TCivetweb(Bool_t only_secured)
   : THttpEngine("civetweb", "compact embedded http server"), fCtx(nullptr), fCallbacks(nullptr), fTopName(),
     fDebug(kFALSE), fTerminating(kFALSE), fOnlySecured(only_secured)
{
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TCivetweb::~TCivetweb()
{
   if (fCtx && !fTerminating)
      mg_stop((struct mg_context *)fCtx);
   if (fCallbacks)
      free(fCallbacks);
}

////////////////////////////////////////////////////////////////////////////////
/// process civetweb log message, can be used to detect critical errors

Int_t TCivetweb::ProcessLog(const char *message)
{
   if ((gDebug > 0) || (strstr(message, "cannot bind to") != 0))
      Error("Log", "%s", message);

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates embedded civetweb server
/// As main argument, http port should be specified like "8090".
/// Or one can provide combination of ipaddress and portnumber like 127.0.0.1:8090
/// Extra parameters like in URL string could be specified after '?' mark:
///    thrds=N   - there N is number of threads used by the civetweb (default is 10)
///    top=name  - configure top name, visible in the web browser
///    ssl_certificate=filename - SSL certificate, see docs/OpenSSL.md from civetweb
///    auth_file=filename  - authentication file name, created with htdigets utility
///    auth_domain=domain   - authentication domain
///    websocket_timeout=tm  - set web sockets timeout in seconds (default 300)
///    websocket_disable - disable web sockets handling (default enabled)
///    bind - ip address to bind server socket
///    loopback  - bind specified port to loopback 127.0.0.1 address
///    debug   - enable debug mode, server always returns html page with request info
///    log=filename  - configure civetweb log file
///  Examples:
///     http:8080?websocket_disable
///     http:7546?thrds=30&websocket_timeout=20

Bool_t TCivetweb::Create(const char *args)
{
   fCallbacks = malloc(sizeof(struct mg_callbacks));
   memset(fCallbacks, 0, sizeof(struct mg_callbacks));
   //((struct mg_callbacks *) fCallbacks)->begin_request = begin_request_handler;
   ((struct mg_callbacks *)fCallbacks)->log_message = log_message_handler;
   TString sport = IsSecured() ? "8480s" : "8080", num_threads = "10", websocket_timeout = "300000";
   TString auth_file, auth_domain, log_file, ssl_cert;
   Bool_t use_ws = kTRUE;

   // extract arguments
   if (args && (strlen(args) > 0)) {

      // first extract port number
      sport = "";
      while ((*args != 0) && (*args != '?') && (*args != '/'))
         sport.Append(*args++);
      if (IsSecured() && (sport.Index("s")==kNPOS)) sport.Append("s");

      // than search for extra parameters
      while ((*args != 0) && (*args != '?'))
         args++;

      if (*args == '?') {
         TUrl url(TString::Format("http://localhost/folder%s", args));

         if (url.IsValid()) {
            url.ParseOptions();

            const char *top = url.GetValueFromOptions("top");
            if (top)
               fTopName = top;

            const char *log = url.GetValueFromOptions("log");
            if (log)
               log_file = log;

            Int_t thrds = url.GetIntValueFromOptions("thrds");
            if (thrds > 0)
               num_threads.Form("%d", thrds);

            const char *afile = url.GetValueFromOptions("auth_file");
            if (afile)
               auth_file = afile;

            const char *adomain = url.GetValueFromOptions("auth_domain");
            if (adomain)
               auth_domain = adomain;

            const char *sslc = url.GetValueFromOptions("ssl_cert");
            if (sslc)
               ssl_cert = sslc;

            Int_t wtmout = url.GetIntValueFromOptions("websocket_timeout");
            if (wtmout > 0) {
               websocket_timeout.Format("%d", wtmout * 1000);
               use_ws = kTRUE;
            }

            if (url.HasOption("websocket_disable"))
               use_ws = kFALSE;

            if (url.HasOption("debug"))
               fDebug = kTRUE;

            if (url.HasOption("loopback") && (sport.Index(":") == kNPOS))
               sport = TString("127.0.0.1:") + sport;

            if (url.HasOption("bind") && (sport.Index(":") == kNPOS)) {
               const char *addr = url.GetValueFromOptions("bind");
               if (addr && strlen(addr))
                  sport = TString(addr) + ":" + sport;
            }

            if (GetServer() && url.HasOption("cors")) {
               const char *cors = url.GetValueFromOptions("cors");
               GetServer()->SetCors(cors && *cors ? cors : "*");
            }
         }
      }
   }

   const char *options[20];
   int op(0);

   Info("Create", "Starting HTTP server on port %s", sport.Data());

   options[op++] = "listening_ports";
   options[op++] = sport.Data();
   options[op++] = "num_threads";
   options[op++] = num_threads.Data();

   if (use_ws) {
      options[op++] = "websocket_timeout_ms";
      options[op++] = websocket_timeout.Data();
   }

   if ((auth_file.Length() > 0) && (auth_domain.Length() > 0)) {
      options[op++] = "global_auth_file";
      options[op++] = auth_file.Data();
      options[op++] = "authentication_domain";
      options[op++] = auth_domain.Data();
   }

   if (log_file.Length() > 0) {
      options[op++] = "error_log_file";
      options[op++] = log_file.Data();
   }

   if (ssl_cert.Length() > 0) {
      options[op++] = "ssl_certificate";
      options[op++] = ssl_cert.Data();
   } else if (IsSecured())
      Error("Create", "No SSL certificate file configured");

   options[op++] = nullptr;

   // Start the web server.
   fCtx = mg_start((struct mg_callbacks *)fCallbacks, this, options);

   if (!fCtx)
      return kFALSE;

   mg_set_request_handler((struct mg_context *)fCtx, "/", begin_request_handler, nullptr);

   if (use_ws)
      mg_set_websocket_handler((struct mg_context *)fCtx, "**root.websocket$", websocket_connect_handler,
                               websocket_ready_handler, websocket_data_handler, websocket_close_handler, nullptr);

   return kTRUE;
}
