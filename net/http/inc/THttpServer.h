// $Id$
// Author: Sergey Linev   21/12/2013

/*************************************************************************
 * Copyright (C) 1995-2013, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_THttpServer
#define ROOT_THttpServer

#include "TNamed.h"
#include "TList.h"
#include "THttpCallArg.h"

#include <mutex>
#include <map>
#include <string>
#include <memory>
#include <queue>
#include <thread>
#include <vector>

class THttpEngine;
class THttpTimer;
class TRootSniffer;

class THttpServer : public TNamed {

protected:
   TList fEngines;                      ///<! engines which runs http server
   std::unique_ptr<THttpTimer> fTimer;   ///<! timer used to access main thread
   std::unique_ptr<TRootSniffer> fSniffer; ///<! sniffer provides access to ROOT objects hierarchy
   Bool_t fTerminated{kFALSE};          ///<! termination flag, disables all requests processing
   Long_t fMainThrdId{0};               ///<! id of the thread for processing requests
   Bool_t fOwnThread{kFALSE};           ///<! true when specialized thread allocated for processing requests
   std::thread fThrd;                   ///<! own thread
   Bool_t fWSOnly{kFALSE};              ///<! when true, handle only websockets / longpoll engine

   TString fJSROOTSYS;       ///<! location of local JSROOT files
   TString fTopName{"ROOT"}; ///<! name of top folder, default - "ROOT"
   TString fJSROOT;          ///<! location of external JSROOT files

   std::map<std::string, std::string> fLocations; ///<! list of local directories, which could be accessed via server

   std::string fDefaultPage;     ///<! file name for default page name
   std::string fDefaultPageCont; ///<! content of default html page
   std::string fDrawPage;        ///<! file name for drawing of single element
   std::string fDrawPageCont;    ///<! content of draw html page
   std::string fCors;            ///<! CORS: sets Access-Control-Allow-Origin header for ProcessRequest responses

   std::mutex fMutex;                                        ///<! mutex to protect list with arguments
   std::queue<std::shared_ptr<THttpCallArg>> fArgs;          ///<! submitted arguments

   std::mutex fWSMutex;                                      ///<! mutex to protect WS handler lists
   std::vector<std::shared_ptr<THttpWSHandler>> fWSHandlers; ///<! list of WS handlers

   virtual void MissedRequest(THttpCallArg *arg);

   virtual void ProcessRequest(std::shared_ptr<THttpCallArg> arg);

   virtual void ProcessBatchHolder(std::shared_ptr<THttpCallArg> &arg);

   void StopServerThread();

   std::string BuildWSEntryPage();

   void ReplaceJSROOTLinks(std::shared_ptr<THttpCallArg> &arg);

   static Bool_t VerifyFilePath(const char *fname);

   THttpServer(const THttpServer &) = delete;
   THttpServer &operator=(const THttpServer &) = delete;

public:
   THttpServer(const char *engine = "http:8080");
   virtual ~THttpServer();

   Bool_t CreateEngine(const char *engine);

   Bool_t IsAnyEngine() const { return fEngines.GetSize() > 0; }

   /** returns pointer on objects sniffer */
   TRootSniffer *GetSniffer() const { return fSniffer.get(); }

   void SetSniffer(TRootSniffer *sniff);

   Bool_t IsReadOnly() const;

   void SetReadOnly(Bool_t readonly = kTRUE);

   Bool_t IsWSOnly() const;

   void SetWSOnly(Bool_t on = kTRUE);

   /** set termination flag, no any further requests will be processed */
   void SetTerminate();

   /** returns kTRUE, if server was terminated */
   Bool_t IsTerminated() const { return fTerminated; }

   /** Enable CORS header to ProcessRequests() responses
    * Specified location (typically "*") add as "Access-Control-Allow-Origin" header */
   void SetCors(const std::string &domain = "*") { fCors = domain; }

   /** Returns kTRUE if CORS was configured */
   Bool_t IsCors() const { return !fCors.empty(); }

   /** Returns specified CORS domain */
   const char *GetCors() const { return fCors.c_str(); }

   /** set name of top item in objects hierarchy */
   void SetTopName(const char *top) { fTopName = top; }

   /** returns name of top item in objects hierarchy */
   const char *GetTopName() const { return fTopName.Data(); }

   void SetJSROOT(const char *location);

   void AddLocation(const char *prefix, const char *path);

   void SetDefaultPage(const std::string &filename = "");

   void SetDrawPage(const std::string &filename = "");

   void SetTimer(Long_t milliSec = 100, Bool_t mode = kTRUE);

   void CreateServerThread();

   /** Check if file is requested, thread safe */
   Bool_t IsFileRequested(const char *uri, TString &res) const;

   /** Execute HTTP request */
   Bool_t ExecuteHttp(std::shared_ptr<THttpCallArg> arg);

   /** Submit HTTP request */
   Bool_t SubmitHttp(std::shared_ptr<THttpCallArg> arg, Bool_t can_run_immediately = kFALSE);

   /** Process submitted requests, must be called from appropriate thread */
   Int_t ProcessRequests();

   /** Register object in subfolder */
   Bool_t Register(const char *subfolder, TObject *obj);

   /** Unregister object */
   Bool_t Unregister(TObject *obj);

   /** Register WS handler*/
   void RegisterWS(std::shared_ptr<THttpWSHandler> ws);

   /** Unregister WS handler*/
   void UnregisterWS(std::shared_ptr<THttpWSHandler> ws);

   /** Find web-socket handler with given name */
   std::shared_ptr<THttpWSHandler> FindWS(const char *name);

   /** Execute WS request */
   Bool_t ExecuteWS(std::shared_ptr<THttpCallArg> &arg, Bool_t external_thrd = kFALSE, Bool_t wait_process = kFALSE);

   /** Restrict access to specified object */
   void Restrict(const char *path, const char *options);

   Bool_t RegisterCommand(const char *cmdname, const char *method, const char *icon = nullptr);

   Bool_t Hide(const char *fullname, Bool_t hide = kTRUE);

   Bool_t SetIcon(const char *fullname, const char *iconname);

   Bool_t CreateItem(const char *fullname, const char *title);

   Bool_t SetItemField(const char *fullname, const char *name, const char *value);

   const char *GetItemField(const char *fullname, const char *name);

   /** Guess mime type base on file extension */
   static const char *GetMimeType(const char *path);

   /** Reads content of file from the disk */
   static char *ReadFileContent(const char *filename, Int_t &len);

   /** Reads content of file from the disk, use std::string in return value */
   static std::string ReadFileContent(const std::string &filename);

   ClassDefOverride(THttpServer, 0) // HTTP server for ROOT analysis
};

#endif
