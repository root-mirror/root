// Author: Sergey Linev <s.linev@gsi.de>
// Date: 2018-08-20
// Warning: This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RWebWindowWSHandler
#define ROOT7_RWebWindowWSHandler

#include "THttpWSHandler.h"
#include "TEnv.h"
#include "TUrl.h"

#include <ROOT/RWebWindow.hxx>

#include <string>

using namespace std::string_literals;

namespace ROOT {
namespace Experimental {

/// just wrapper to deliver websockets call-backs to the RWebWindow class

class RWebWindowWSHandler : public THttpWSHandler {

protected:
   Bool_t ProcessBatchHolder(std::shared_ptr<THttpCallArg> &arg) override
   {
      return IsDisabled() ? kFALSE : fWindow.ProcessBatchHolder(arg);
   }

   void VerifyDefaultPageContent(std::shared_ptr<THttpCallArg> &arg) override
   {
      auto token = fWindow.GetConnToken();
      if (!token.empty()) {
         TUrl url;
         url.SetOptions(arg->GetQuery());
         // refuse connection which does not provide proper token
         if (!url.HasOption("token") || (token != url.GetValueFromOptions("token"))) {
            // refuce loading of default web page without token
            arg->SetContent("refused");
            arg->Set404();
            return;
         }
      }

      auto version = fWindow.GetClientVersion();
      if (!version.empty()) {
         std::string search = "jsrootsys/scripts/JSRoot.core."s;
         std::string replace = version + "/jsrootsys/scripts/JSRoot.core."s;
         // replace link to JSROOT main script to emulate new version
         arg->ReplaceAllinContent(search, replace, true);
         arg->AddNoCacheHeader();
      }

      std::string more_args;
      const char *ui5source = gEnv->GetValue("WebGui.openui5src","");
      if (ui5source && *ui5source)
         more_args.append("openui5src: \""s + ui5source + "\","s);
      const char *ui5libs = gEnv->GetValue("WebGui.openui5libs","");
      if (ui5libs && *ui5libs)
         more_args.append("openui5libs: \""s + ui5libs + "\","s);
      const char *ui5theme = gEnv->GetValue("WebGui.openui5theme","");
      if (ui5theme && *ui5theme)
         more_args.append("openui5theme: \""s + ui5theme + "\","s);
      int credits = gEnv->GetValue("WebGui.ConnCredits", 10);
      if ((credits > 0) && (credits != 10))
         more_args.append("credits: "s + std::to_string(credits) + ","s);
      auto user_args = fWindow.GetUserArgs();
      if (!user_args.empty())
         more_args.append("user_args: "s + user_args + ","s);
      if (!more_args.empty()) {
         std::string search = "JSROOT.connectWebWindow({"s;
         std::string replace = search + more_args;
         arg->ReplaceAllinContent(search, replace, true);
         arg->AddNoCacheHeader();
      }
   }

public:
   RWebWindow &fWindow; ///<! window reference

   /// constructor
   RWebWindowWSHandler(RWebWindow &wind, const char *name)
      : THttpWSHandler(name, "RWebWindow websockets handler", kFALSE), fWindow(wind)
   {
   }

   virtual ~RWebWindowWSHandler() = default;

   /// returns content of default web-page
   /// THttpWSHandler interface
   TString GetDefaultPageContent() override { return IsDisabled() ? "" : fWindow.fDefaultPage.c_str(); }

   /// returns true when window allowed to serve files relative to default page
   Bool_t CanServeFiles() const override { return !IsDisabled(); }

   /// Process websocket request - called from THttpServer thread
   /// THttpWSHandler interface
   Bool_t ProcessWS(THttpCallArg *arg) override { return arg && !IsDisabled() ? fWindow.ProcessWS(*arg) : kFALSE; }

   /// Allow processing of WS actions in arbitrary thread
   Bool_t AllowMTProcess() const override { return fWindow.fProcessMT; }

   /// Allows usage of special threads for send operations
   Bool_t AllowMTSend() const override { return fWindow.fSendMT; }

   /// React on completion of multi-threaded send operation
   void CompleteWSSend(UInt_t wsid) override { if (!IsDisabled()) fWindow.CompleteWSSend(wsid); }
};

} // namespace Experimental
} // namespace ROOT

#endif
