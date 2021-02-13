// $Id$
// Author: Sergey Linev   21/05/2015

/*************************************************************************
 * Copyright (C) 1995-2013, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_THttpCallArg
#define ROOT_THttpCallArg

#include "TObject.h"

#include "TString.h"

#include <condition_variable>
#include <string>
#include <memory>

class THttpServer;
class THttpWSEngine;
class THttpWSHandler;

class THttpCallArg : public TObject {

   friend class THttpServer;
   friend class THttpWSEngine;
   friend class THttpWSHandler;

public:
   enum {
      kNoZip     = 0,             // no zipping
      kZip       = 1,             // zip content if "Accept-Encoding" header contains "gzip"
      kZipLarge  = 2,             // zip if content larger than 10K and "Accept-Encoding" contains "gzip"
      kZipAlways = 3              // zip always
   };

protected:
   TString fTopName;              ///<! top item name
   TString fMethod;               ///<! request method like GET or POST
   TString fPathName;             ///<! item path
   TString fFileName;             ///<! file name
   TString fUserName;             ///<! authenticated user name (if any)
   TString fQuery;                ///<! additional arguments

   UInt_t fWSId{0};               ///<! websocket identifier, used in web-socket related operations

   std::condition_variable fCond; ///<! condition used to wait for processing

   TString fContentType;          ///<! type of content
   TString fRequestHeader;        ///<! complete header, provided with request
   TString fHeader;               ///<! response header like ContentEncoding, Cache-Control and so on
   Int_t fZipping{kNoZip};        ///<! indicate if and when content should be compressed

   Bool_t fNotifyFlag{kFALSE};    ///<!  indicate that notification called

   TString AccessHeader(TString &buf, const char *name, const char *value = nullptr, Bool_t doing_set = kFALSE);

   TString CountHeader(const TString &buf, Int_t number = -1111) const;

   /** Method used to modify content of web page used by web socket handler */
   virtual void CheckWSPageContent(THttpWSHandler *) {}

private:
   std::shared_ptr<THttpWSEngine> fWSEngine; ///<!  web-socket engine, which supplied to run created web socket

   std::string  fContent;  ///<! content - text or binary
   std::string  fPostData; ///<! data received with post request - text - or binary

   void AssignWSId();
   std::shared_ptr<THttpWSEngine> TakeWSEngine();

public:
   explicit THttpCallArg() {} // NOLINT: not allowed to use = default because of TObject::kIsOnHeap detection, see ROOT-10300
   virtual ~THttpCallArg();

   // these methods used to set http request arguments

   /** set request method kind like GET or POST */
   void SetMethod(const char *method) { fMethod = method; }

   /** set engine-specific top-name */
   void SetTopName(const char *topname) { fTopName = topname; }

   void SetPathAndFileName(const char *fullpath);

   /** set request path name */
   void SetPathName(const char *p) { fPathName = p; }

   /** set request file name */
   void SetFileName(const char *f) { fFileName = f; }

   /** set name of authenticated user */
   void SetUserName(const char *n) { fUserName = n; }

   /** set request query */
   void SetQuery(const char *q) { fQuery = q; }

   void SetPostData(void *data, Long_t length, Bool_t make_copy = kFALSE);

   void SetPostData(std::string &&data);

   /** set web-socket id */
   void SetWSId(UInt_t id) { fWSId = id; }

   /** get web-socket id */
   UInt_t GetWSId() const { return fWSId; }

   /** set full set of request header */
   void SetRequestHeader(const char *h) { fRequestHeader = (h ? h : ""); }

   /** returns number of fields in request header */
   Int_t NumRequestHeader() const { return CountHeader(fRequestHeader).Atoi(); }

   /** returns field name in request header */
   TString GetRequestHeaderName(Int_t number) const { return CountHeader(fRequestHeader, number); }

   /** get named field from request header */
   TString GetRequestHeader(const char *name) { return AccessHeader(fRequestHeader, name); }

   /** returns engine-specific top-name */
   const char *GetTopName() const { return fTopName.Data(); }

   /** returns request method like GET or POST */
   const char *GetMethod() const { return fMethod.Data(); }

   /** returns kTRUE if post method is used */
   Bool_t IsMethod(const char *name) const { return fMethod.CompareTo(name) == 0; }

   /** returns kTRUE if post method is used */
   Bool_t IsPostMethod() const { return IsMethod("POST"); }

   /** return pointer on posted with request data */
   const void *GetPostData() const { return fPostData.data(); }

   /** return length of posted with request data */
   Long_t GetPostDataLength() const { return (Long_t) fPostData.length(); }

   /** returns path name from request URL */
   const char *GetPathName() const { return fPathName.Data(); }

   /** returns file name from request URL */
   const char *GetFileName() const { return fFileName.Data(); }

   /** return authenticated user name (0 - when no authentication) */
   const char *GetUserName() const { return fUserName.Length() > 0 ? fUserName.Data() : nullptr; }

   /** returns request query (string after ? in request URL) */
   const char *GetQuery() const { return fQuery.Data(); }

   // these methods used in THttpServer to set results of request processing

   /** set content type like "text/xml" or "application/json" */
   void SetContentType(const char *typ) { fContentType = typ; }

   /** mark reply as 404 error - page/request not exists or refused */
   void Set404() { SetContentType("_404_"); }

   /** Return true if reply can be postponed by server  */
   virtual Bool_t CanPostpone() const { return kTRUE; }

   /** mark as postponed - reply will not be send to client immediately */
   void SetPostponed()
   {
      if (CanPostpone())
         SetContentType("_postponed_");
      else
         Set404();
   }

   /** indicate that http request should response with file content */
   void SetFile(const char *filename = nullptr)
   {
      SetContentType("_file_");
      if (filename)
         fContent = filename;
   }

   void SetText();
   void SetTextContent(std::string &&txt);

   void SetXml();
   void SetXmlContent(std::string &&xml);

   void SetJson();
   void SetJsonContent(std::string &&json);

   void SetBinary();
   void SetBinaryContent(std::string &&bin);

   void AddHeader(const char *name, const char *value);

   void AddNoCacheHeader();

   /** returns number of fields in header */
   Int_t NumHeader() const { return CountHeader(fHeader).Atoi(); }

   /** returns field name in header */
   TString GetHeaderName(Int_t number) const { return CountHeader(fHeader, number); }

   TString GetHeader(const char *name);

   /** Set Content-Encoding header like gzip */
   void SetEncoding(const char *typ) { AccessHeader(fHeader, "Content-Encoding", typ, kTRUE); }

   void SetContent(const char *cont);
   void SetContent(std::string &&cont);
   void ReplaceAllinContent(const std::string &from, const std::string &to, bool once = false);

   Bool_t CompressWithGzip();

   void SetZipping(Int_t mode = kZipLarge) { fZipping = mode; }
   Int_t GetZipping() const { return fZipping; }

   /** add extra http header value to the reply */
   void SetExtraHeader(const char *name, const char *value) { AddHeader(name, value); }

   std::string FillHttpHeader(const char *header = nullptr);

   // these methods used to return results of http request processing

   Bool_t IsContentType(const char *typ) const { return fContentType == typ; }
   const char *GetContentType() const { return fContentType.Data(); }

   Bool_t Is404() const { return IsContentType("_404_"); }
   Bool_t IsFile() const { return IsContentType("_file_"); }
   Bool_t IsPostponed() const { return IsContentType("_postponed_"); }
   Bool_t IsText() const { return IsContentType("text/plain"); }
   Bool_t IsXml() const { return IsContentType("text/xml"); }
   Bool_t IsJson() const { return IsContentType("application/json"); }
   Bool_t IsBinary() const { return IsContentType("application/x-binary"); }

   Long_t GetContentLength() const { return (Long_t) fContent.length(); }
   const void *GetContent() const { return fContent.data(); }

   void NotifyCondition();

   virtual void HttpReplied();

   template <class T, typename... Args>
   void CreateWSEngine(Args... args)
   {
      fWSEngine = std::make_shared<T>(args...);
      AssignWSId();
   }

   ClassDef(THttpCallArg, 0) // Arguments for single HTTP call
};

#endif
