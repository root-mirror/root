// $Id$
// Author: Sergey Linev   28/12/2013

/*************************************************************************
 * Copyright (C) 1995-2013, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TFastCgi.h"

#include "TThread.h"
#include "TUrl.h"
#include "THttpServer.h"

#include "ROOT/RMakeUnique.hxx"

#include <cstring>

#ifdef WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

////////////////////////////////////////////////////////////////////////////////

class TFastCgiCallArg : public THttpCallArg {

   bool fCanPostpone{false};

protected:

   void CheckWSPageContent(THttpWSHandler *) override
   {
      std::string search = "JSROOT.connectWebWindow({";
      std::string replace = search + "socket_kind:\"longpoll\",";
      ReplaceAllinContent(search, replace, true);
   }

public:
   TFastCgiCallArg(bool can_postpone) : THttpCallArg(), fCanPostpone(can_postpone) {};

   /** All FastCGI requests should be immediately replied to get slot for next */
   Bool_t CanPostpone() const override { return fCanPostpone; }
};

////////////////////////////////////////////////////////////////////////////////


#ifndef HTTP_WITHOUT_FASTCGI

#include "fcgiapp.h"

#include <stdlib.h>

void FCGX_ROOT_send_file(FCGX_Request *request, const char *fname)
{
   std::string buf = THttpServer::ReadFileContent(fname);

   if (buf.empty()) {
      FCGX_FPrintF(request->out,
                   "Status: 404 Not Found\r\n"
                   "Content-Length: 0\r\n" // Always set Content-Length
                   "Connection: close\r\n\r\n");
   } else {

      FCGX_FPrintF(request->out,
                   "Status: 200 OK\r\n"
                   "Content-Type: %s\r\n"
                   "Content-Length: %d\r\n" // Always set Content-Length
                   "\r\n",
                   THttpServer::GetMimeType(fname), (int) buf.length());

      FCGX_PutStr(buf.c_str(), buf.length(), request->out);
   }
}


void process_request(TFastCgi *engine, FCGX_Request *request, bool can_postpone)
{
   int count = 0;
   count++;  // simple static request counter

   const char *inp_path = FCGX_GetParam("PATH_INFO", request->envp);
   if (!inp_path) inp_path = FCGX_GetParam("SCRIPT_FILENAME", request->envp);
   const char *inp_query = FCGX_GetParam("QUERY_STRING", request->envp);
   const char *inp_method = FCGX_GetParam("REQUEST_METHOD", request->envp);
   const char *inp_length = FCGX_GetParam("CONTENT_LENGTH", request->envp);

   auto arg = std::make_shared<TFastCgiCallArg>(can_postpone);
   if (inp_path)
      arg->SetPathAndFileName(inp_path);
   if (inp_query)
      arg->SetQuery(inp_query);
   if (inp_method)
      arg->SetMethod(inp_method);
   if (engine->GetTopName())
      arg->SetTopName(engine->GetTopName());
   int len = 0;
   if (inp_length)
      len = strtol(inp_length, nullptr, 10);
   if (len > 0) {
      std::string buf;
      buf.resize(len);
      int nread = FCGX_GetStr((char *)buf.data(), len, request->in);
      if (nread == len)
         arg->SetPostData(std::move(buf));
   }

   TString header;
   for (char **envp = request->envp; *envp != nullptr; envp++) {
      TString entry = *envp;
      for (Int_t n = 0; n < entry.Length(); n++)
         if (entry[n] == '=') {
            entry[n] = ':';
            break;
         }
      header.Append(entry);
      header.Append("\r\n");
   }
   arg->SetRequestHeader(header);

   TString username = arg->GetRequestHeader("REMOTE_USER");
   if ((username.Length() > 0) && (arg->GetRequestHeader("AUTH_TYPE").Length() > 0))
      arg->SetUserName(username);

   if (engine->IsDebugMode()) {
      FCGX_FPrintF(request->out, "Status: 200 OK\r\n"
                                "Content-type: text/html\r\n"
                                "\r\n"
                                "<title>FastCGI echo</title>"
                                "<h1>FastCGI echo</h1>\n");

      FCGX_FPrintF(request->out, "Request %d:<br/>\n<pre>\n", count);
      FCGX_FPrintF(request->out, "  Method   : %s\n", arg->GetMethod());
      FCGX_FPrintF(request->out, "  PathName : %s\n", arg->GetPathName());
      FCGX_FPrintF(request->out, "  FileName : %s\n", arg->GetFileName());
      FCGX_FPrintF(request->out, "  Query    : %s\n", arg->GetQuery());
      FCGX_FPrintF(request->out, "  PostData : %ld\n", arg->GetPostDataLength());
      FCGX_FPrintF(request->out, "</pre><p>\n");

      FCGX_FPrintF(request->out, "Environment:<br/>\n<pre>\n");
      for (char **envp = request->envp; *envp != nullptr; envp++)
         FCGX_FPrintF(request->out, "  %s\n", *envp);
      FCGX_FPrintF(request->out, "</pre><p>\n");

      return;
   }

   TString fname;

   if (engine->GetServer()->IsFileRequested(inp_path, fname)) {
      FCGX_ROOT_send_file(request, fname.Data());
      return;
   }

   if (!engine->GetServer()->ExecuteHttp(arg) || arg->Is404()) {
      std::string hdr = arg->FillHttpHeader("Status:");
      FCGX_FPrintF(request->out, hdr.c_str());
   } else if (arg->IsFile()) {
      FCGX_ROOT_send_file(request, (const char *)arg->GetContent());
   } else {

      // TODO: check in request header that gzip encoding is supported
      if (arg->GetZipping() != THttpCallArg::kNoZip)
         arg->CompressWithGzip();

      std::string hdr = arg->FillHttpHeader("Status:");
      FCGX_FPrintF(request->out, hdr.c_str());

      FCGX_PutStr((const char *)arg->GetContent(), (int)arg->GetContentLength(), request->out);
   }
}

void run_multi_threads(TFastCgi *engine, Int_t nthrds)
{
   std::condition_variable cond; ///<! condition used to wait for processing
   std::mutex m;
   std::unique_ptr<FCGX_Request> arg;
   int nwaiting = 0;

   auto worker_func = [engine, &cond, &m, &arg, &nwaiting]() {

      while (!engine->IsTerminating()) {

         std::unique_ptr<FCGX_Request> request;

         bool can_postpone = false;

         {
            std::unique_lock<std::mutex> lk(m);
            nwaiting++;
            cond.wait(lk);
            nwaiting--;
            can_postpone = (nwaiting > 5);
            std::swap(arg, request);
         }

         if (request) {
            process_request(engine, request.get(), can_postpone);

            FCGX_Finish_r(request.get());
         }
      }

   };

   // start N workers
   std::vector<std::thread> workers;
   for (int n=0; n< nthrds; ++n)
      workers.emplace_back(worker_func);

   while (!engine->IsTerminating()) {
      auto request = std::make_unique<FCGX_Request>();

      FCGX_InitRequest(request.get(), engine->GetSocket(), 0);

      int rc = FCGX_Accept_r(request.get());

      if (rc != 0)
         continue;

      {
         std::lock_guard<std::mutex> lk(m);
         if (nwaiting > 0)
            std::swap(request, arg);
      }

      if (!request) {
         // notify thread to process request
         cond.notify_one();
      } else {
         // process request ourselfs
         process_request(engine, request.get(), false);
         FCGX_Finish_r(request.get());
      }
   }

   // ensure that all threads are waked up
   cond.notify_all();

   // join all workers
   for (auto & thrd : workers)
      thrd.join();

}

// simple run function to process all requests in same thread

void run_single_thread(TFastCgi *engine)
{

   FCGX_Request request;

   FCGX_InitRequest(&request, engine->GetSocket(), 0);

   while (!engine->IsTerminating()) {

      int rc = FCGX_Accept_r(&request);

      if (rc != 0)
         continue;

      process_request(engine, &request, false);

      FCGX_Finish_r(&request);
   }

}

#endif


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TFastCgi                                                             //
//                                                                      //
// http engine implementation, based on fastcgi package                 //
// Allows to redirect http requests from normal web server like         //
// Apache or lighttpd                                                   //
//                                                                      //
// Configuration example for lighttpd                                   //
//                                                                      //
// server.modules += ( "mod_fastcgi" )                                  //
// fastcgi.server = (                                                   //
//   "/remote_scripts/" =>                                              //
//     (( "host" => "192.168.1.11",                                     //
//        "port" => 9000,                                               //
//        "check-local" => "disable",                                   //
//        "docroot" => "/"                                              //
//     ))                                                               //
// )                                                                    //
//                                                                      //
// When creating THttpServer, one should specify:                       //
//                                                                      //
//  THttpServer* serv = new THttpServer("fastcgi:9000");                //
//                                                                      //
// In this case, requests to lighttpd server will be                    //
// redirected to ROOT session. Like:                                    //
//    http://lighttpdhost/remote_scripts/root.cgi/                      //
//                                                                      //
// Following additional options can be specified                        //
//    top=foldername - name of top folder, seen in the browser          //
//    thrds=N - run N worker threads to process requests, default 10    //
//    debug=1 - run fastcgi server in debug mode                        //
// Example:                                                             //
//    serv->CreateEngine("fastcgi:9000?top=fastcgiserver");             //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// normal constructor

TFastCgi::TFastCgi()
   : THttpEngine("fastcgi", "fastcgi interface to webserver")
{
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TFastCgi::~TFastCgi()
{
   fTerminating = kTRUE;

   // running thread will stopped
   if (fThrd)
      fThrd->join();

   if (fSocket > 0) {
      // close opened socket
      close(fSocket);
      fSocket = 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// initializes fastcgi variables and start thread,
/// which will process incoming http requests

Bool_t TFastCgi::Create(const char *args)
{
#ifndef HTTP_WITHOUT_FASTCGI
   FCGX_Init();

   TString sport = ":9000";
   Int_t nthrds = 10;

   if ((args != 0) && (strlen(args) > 0)) {

      // first extract port number
      sport = ":";
      while ((*args != 0) && (*args >= '0') && (*args <= '9'))
         sport.Append(*args++);

      // than search for extra parameters
      while ((*args != 0) && (*args != '?'))
         args++;

      if (*args == '?') {
         TUrl url(TString::Format("http://localhost/folder%s", args));

         if (url.IsValid()) {

            url.ParseOptions();

            if (url.GetValueFromOptions("debug") != 0)
               fDebugMode = kTRUE;

            if (url.HasOption("thrds"))
               nthrds = url.GetIntValueFromOptions("thrds");

            const char *top = url.GetValueFromOptions("top");
            if (top != 0)
               fTopName = top;
         }
      }
   }

   Info("Create", "Starting FastCGI server on port %s", sport.Data() + 1);

   fSocket = FCGX_OpenSocket(sport.Data(), 10);
   if (!fSocket) return kFALSE;

   if (nthrds > 0)
      fThrd = std::make_unique<std::thread>(run_multi_threads, this, nthrds);
   else
      fThrd = std::make_unique<std::thread>(run_single_thread, this);

   return kTRUE;
#else
   (void) args;
   Error("Create", "ROOT compiled without fastcgi support");
   return kFALSE;
#endif
}
