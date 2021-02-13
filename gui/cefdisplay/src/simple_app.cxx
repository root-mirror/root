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


#if !defined(_MSC_VER)
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wshadow"
#endif

#include "simple_app.h"

#include <string>
#include <cstdio>
#include <memory>

#include "include/cef_browser.h"
#include "include/cef_version.h"
#include "include/views/cef_browser_view.h"
#include "include/views/cef_window.h"
#include "include/wrapper/cef_helpers.h"

#include "THttpServer.h"
#include "TTimer.h"
#include <ROOT/RLogger.hxx>

#include "RCefWebDisplayHandle.hxx"

namespace {

// When using the Views framework this object provides the delegate
// implementation for the CefWindow that hosts the Views-based browser.
class SimpleWindowDelegate : public CefWindowDelegate {
   CefRefPtr<CefBrowserView> fBrowserView;
   int fWidth{800};  ///< preferred window width
   int fHeight{600}; ///< preferred window height
public:
   explicit SimpleWindowDelegate(CefRefPtr<CefBrowserView> browser_view, int width = 800, int height = 600)
      : fBrowserView(browser_view), fWidth(width), fHeight(height)
   {
   }

   void OnWindowCreated(CefRefPtr<CefWindow> window) OVERRIDE
   {
      // Add the browser view and show the window.
      window->AddChildView(fBrowserView);
      window->Show();

      // Give keyboard focus to the browser view.
      fBrowserView->RequestFocus();
   }

   void OnWindowDestroyed(CefRefPtr<CefWindow> window) OVERRIDE { fBrowserView = nullptr; }

   bool CanClose(CefRefPtr<CefWindow> window) OVERRIDE
   {
      // Allow the window to close if the browser says it's OK.
      CefRefPtr<CefBrowser> browser = fBrowserView->GetBrowser();
      if (browser)
         return browser->GetHost()->TryCloseBrowser();
      return true;
   }

   CefSize GetPreferredSize(CefRefPtr<CefView> view) OVERRIDE
   {
     return CefSize(fWidth, fHeight);
   }


private:

   IMPLEMENT_REFCOUNTING(SimpleWindowDelegate);
   DISALLOW_COPY_AND_ASSIGN(SimpleWindowDelegate);
};


class SimpleBrowserViewDelegate : public CefBrowserViewDelegate {
 public:
  SimpleBrowserViewDelegate() {}

  bool OnPopupBrowserViewCreated(CefRefPtr<CefBrowserView> browser_view,
                                 CefRefPtr<CefBrowserView> popup_browser_view,
                                 bool is_devtools) OVERRIDE {
    // Create a new top-level Window for the popup. It will show itself after
    // creation.
    CefWindow::CreateTopLevelWindow(new SimpleWindowDelegate(popup_browser_view));

    // We created the Window.
    return true;
  }

 private:
  IMPLEMENT_REFCOUNTING(SimpleBrowserViewDelegate);
  DISALLOW_COPY_AND_ASSIGN(SimpleBrowserViewDelegate);
};


} // namespace

SimpleApp::SimpleApp(bool use_viewes,
                     THttpServer *serv, const std::string &url, const std::string &cont,
                     int width, int height, bool headless)
   : CefApp(), CefBrowserProcessHandler(), fUseViewes(use_viewes), fFirstServer(serv), fFirstUrl(url), fFirstContent(cont), fFirstHeadless(headless)
{
   fFirstRect.Set(0, 0, width, height);

#if defined(OS_WIN) || defined(OS_LINUX)
   // Create the browser using the Views framework if "--use-views" is specified
   // via the command-line. Otherwise, create the browser using the native
   // platform framework. The Views framework is currently only supported on
   // Windows and Linux.
#else
   if (fUseViewes) {
      R__LOG_ERROR(CefWebDisplayLog()) << "view framework does not supported by CEF on the platform, switching off";
      fUseViewes = false;
   }
#endif

}


void SimpleApp::SetNextHandle(RCefWebDisplayHandle *handle)
{
   fNextHandle = handle;
}


void SimpleApp::OnRegisterCustomSchemes(CefRawPtr<CefSchemeRegistrar> registrar)
{
   // registrar->AddCustomScheme("rootscheme", true, true, true, true, true, true);
  // registrar->AddCustomScheme("rootscheme", true, false, false, true, false, false);
}

void SimpleApp::OnBeforeCommandLineProcessing(const CefString &process_type, CefRefPtr<CefCommandLine> command_line)
{
//   command_line->AppendSwitch("allow-file-access-from-files");
//   command_line->AppendSwitch("disable-web-security");
//   if (fBatch) {
//      command_line->AppendSwitch("disable-gpu");
//      command_line->AppendSwitch("disable-gpu-compositing");
//      command_line->AppendSwitch("disable-gpu-sandbox");
//   }
}

void SimpleApp::OnBeforeChildProcessLaunch(CefRefPtr<CefCommandLine> command_line)
{
//   command_line->AppendSwitch("allow-file-access-from-files");
//   command_line->AppendSwitch("disable-web-security");
//   if (fLastBatch) {
//      command_line->AppendSwitch("disable-webgl");
//      command_line->AppendSwitch("disable-gpu");
//      command_line->AppendSwitch("disable-gpu-compositing");
//   }
}

void SimpleApp::OnContextInitialized()
{
   CEF_REQUIRE_UI_THREAD();

   if (!fFirstUrl.empty() || !fFirstContent.empty()) {
      StartWindow(fFirstServer, fFirstUrl, fFirstContent, fFirstRect);
      fFirstUrl.clear();
      fFirstContent.clear();
   }
}


void SimpleApp::StartWindow(THttpServer *serv, const std::string &addr, const std::string &cont, CefRect &rect)
{
   CEF_REQUIRE_UI_THREAD();

   if (!fGuiHandler)
      fGuiHandler = new GuiHandler(fUseViewes);

   std::string url;

   //bool is_batch = false;

   if(addr.empty() && !cont.empty()) {
      url = fGuiHandler->AddBatchPage(cont);
      // is_batch = true;
   } else if (serv) {
      url = fGuiHandler->MakePageUrl(serv, addr);
   } else {
      url = addr;
   }

   // Specify CEF browser settings here.
   CefBrowserSettings browser_settings;
   // browser_settings.plugins = STATE_DISABLED;
   // browser_settings.file_access_from_file_urls = STATE_ENABLED;
   // browser_settings.universal_access_from_file_urls = STATE_ENABLED;
   // browser_settings.web_security = STATE_DISABLED;

   if (fUseViewes) {
      // Create the BrowserView.
      CefRefPtr<CefBrowserView> browser_view =
         CefBrowserView::CreateBrowserView(fGuiHandler, url, browser_settings, nullptr, nullptr, new SimpleBrowserViewDelegate());

      // Create the Window. It will show itself after creation.
      CefWindow::CreateTopLevelWindow(new SimpleWindowDelegate(browser_view, rect.width, rect.height));

      if (fNextHandle) {
         fNextHandle->SetBrowser(browser_view->GetBrowser());
         fNextHandle = nullptr; // used only once
      }

   } else {

      CefWindowInfo window_info;

      // TODO: Seems to be, to configure window_info.SetAsWindowless,
      // one should implement CefRenderHandler

      #if defined(OS_WIN)
         RECT wnd_rect = {rect.x, rect.y, rect.x + rect.width, rect.y + rect.height};
         if (!rect.IsEmpty()) window_info.SetAsChild(0, wnd_rect);
         // On Windows we need to specify certain flags that will be passed to
         // CreateWindowEx().
         window_info.SetAsPopup(0, "cefsimple");
         //if (is_batch)
         //   window_info.SetAsWindowless(GetDesktopWindow());
      #elif defined(OS_LINUX)
         if (!rect.IsEmpty()) window_info.SetAsChild(0, rect);
         //if (is_batch)
         //   window_info.SetAsWindowless(kNullWindowHandle);
      #else
         if (!rect.IsEmpty())
            window_info.SetAsChild(0, rect.x, rect.y, rect.width, rect.height );
         //if (is_batch)
         //   window_info.SetAsWindowless(kNullWindowHandle);
      #endif

      // Create the first browser window.
      auto browser = CefBrowserHost::CreateBrowserSync(window_info, fGuiHandler, url, browser_settings, nullptr, nullptr);

      if (fNextHandle) {
         fNextHandle->SetBrowser(browser);
         fNextHandle = nullptr; // used only once
      }

   }

}

