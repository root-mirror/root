// Author: Sergey Linev <S.Linev@gsi.de>
// Date: 2017-06-29
// Warning: This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!

// Copyright (c) 2013 The Chromium Embedded Framework Authors. All rights
// reserved. Use of this source code is governed by a BSD-style license that
// can be found in the LICENSE file.

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_cef_simple_app
#define ROOT_cef_simple_app

#include "include/cef_app.h"

#include <string>

#include "gui_handler.h"

class THttpServer;

class RCefWebDisplayHandle;

// Implement application-level callbacks for the browser process.
class SimpleApp : public CefApp,
#if defined(OS_LINUX)
                  public CefPrintHandler,
#endif
                  /*, public CefRenderProcessHandler */
                  public CefBrowserProcessHandler {
protected:
   bool fUseViewes{false};  ///<! is views framework used
   THttpServer *fFirstServer; ///<! first server
   std::string fFirstUrl;   ///<! first URL to open
   std::string fFirstContent; ///<! first page content open
   CefRect fFirstRect;      ///<! original width
   bool fFirstHeadless{false}; ///<! is first window is headless
   RCefWebDisplayHandle *fNextHandle{nullptr}; ///< next handle where browser will be created

   CefRefPtr<GuiHandler> fGuiHandler; ///<! normal handler

public:
   SimpleApp(bool use_viewes,
             THttpServer *serv = nullptr, const std::string &url = "", const std::string &cont = "",
             int width = 0, int height = 0, bool headless = false);

   void SetNextHandle(RCefWebDisplayHandle *handle);

   // CefApp methods:
   virtual CefRefPtr<CefBrowserProcessHandler> GetBrowserProcessHandler() OVERRIDE { return this; }

#if defined(OS_LINUX)
   // only on Linux special print handler is required to return PDF size
   virtual CefRefPtr<CefPrintHandler> GetPrintHandler() OVERRIDE { return this; }
#endif
   // virtual CefRefPtr<CefRenderProcessHandler> GetRenderProcessHandler() OVERRIDE { return this; }

   virtual void OnRegisterCustomSchemes(CefRawPtr<CefSchemeRegistrar> registrar) OVERRIDE;

   // CefBrowserProcessHandler methods:
   virtual void OnContextInitialized() OVERRIDE;

   virtual void
   OnBeforeCommandLineProcessing(const CefString &process_type, CefRefPtr<CefCommandLine> command_line) OVERRIDE;

   virtual void OnBeforeChildProcessLaunch(CefRefPtr<CefCommandLine> command_line) OVERRIDE;

#if defined(OS_LINUX)
   // CefPrintHandler methods
   virtual CefSize GetPdfPaperSize(int device_units_per_inch) OVERRIDE { return CefSize(device_units_per_inch*8.25, device_units_per_inch*11.75); }
   virtual bool OnPrintDialog( CefRefPtr< CefBrowser > browser, bool has_selection, CefRefPtr< CefPrintDialogCallback > callback ) OVERRIDE { return false; }
   virtual bool OnPrintJob( CefRefPtr< CefBrowser > browser, const CefString& document_name, const CefString& pdf_file_path, CefRefPtr< CefPrintJobCallback > callback ) OVERRIDE { return false; }
   virtual void OnPrintReset( CefRefPtr< CefBrowser > browser ) OVERRIDE {}
   virtual void OnPrintSettings( CefRefPtr< CefBrowser > browser, CefRefPtr< CefPrintSettings > settings, bool get_defaults ) OVERRIDE {}
   virtual void OnPrintStart( CefRefPtr< CefBrowser > browser ) OVERRIDE {}
#endif


   void StartWindow(THttpServer *serv, const std::string &url, const std::string &cont, CefRect &rect);

   // CefRenderProcessHandler methods
   // virtual void OnContextCreated(CefRefPtr<CefBrowser> browser, CefRefPtr<CefFrame> frame,
   //                              CefRefPtr<CefV8Context> context) OVERRIDE;

private:
   // Include the default reference counting implementation.
   IMPLEMENT_REFCOUNTING(SimpleApp);
   DISALLOW_COPY_AND_ASSIGN(SimpleApp);
};

#endif // ROOT_cef_simple_app
