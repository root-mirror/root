// Author: Sergey Linev <S.Linev@gsi.de>
// Date: 2020-08-21
// Warning: This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


#if !defined(_MSC_VER)
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wshadow"
#endif


#include "RCefWebDisplayHandle.hxx"

#include "TTimer.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TEnv.h"
#include "TApplication.h"

#include <ROOT/RMakeUnique.hxx>
#include <ROOT/RLogger.hxx>


class TCefTimer : public TTimer {
public:
   TCefTimer(Long_t milliSec, Bool_t mode) : TTimer(milliSec, mode) {}
   void Timeout() override
   {
      // just let run loop
      CefDoMessageLoopWork();

   }
};


class FrameSourceVisitor : public CefStringVisitor {

   RCefWebDisplayHandle *fHandle{nullptr};

public:

   FrameSourceVisitor(RCefWebDisplayHandle *handle) : CefStringVisitor(), fHandle(handle) {}

   virtual ~FrameSourceVisitor() = default;

   void Visit(const CefString &str) override
   {
      if (fHandle && fHandle->IsValid())
         fHandle->SetContent(str.ToString());
   }

private:
   // Include the default reference counting implementation.
   IMPLEMENT_REFCOUNTING(FrameSourceVisitor);
   DISALLOW_COPY_AND_ASSIGN(FrameSourceVisitor);
};


class HeadlessPrintCallback : public CefPdfPrintCallback {
   bool *fFlag{nullptr};
public:
   HeadlessPrintCallback(bool *flag) : CefPdfPrintCallback(), fFlag(flag) {}
   virtual ~HeadlessPrintCallback() = default;

   void OnPdfPrintFinished(const CefString&, bool ok ) override
   {
      if (fFlag) *fFlag = true;
   }
private:
   // Include the default reference counting implementation.
   IMPLEMENT_REFCOUNTING(HeadlessPrintCallback);
   DISALLOW_COPY_AND_ASSIGN(HeadlessPrintCallback);
};

std::unique_ptr<ROOT::Experimental::RWebDisplayHandle> RCefWebDisplayHandle::CefCreator::Display(const ROOT::Experimental::RWebDisplayArgs &args)
{

   auto handle = std::make_unique<RCefWebDisplayHandle>(args.GetFullUrl());

   if (fCefApp) {
      fCefApp->SetNextHandle(handle.get());

      CefRect rect((args.GetX() > 0) ? args.GetX() : 0, (args.GetY() > 0) ? args.GetY() : 0,
            (args.GetWidth() > 0) ? args.GetWidth() : 800, (args.GetHeight() > 0) ? args.GetHeight() : 600);

      fCefApp->StartWindow(args.GetHttpServer(), args.GetFullUrl(), args.GetPageContent(), rect);

      if (args.IsHeadless())
         handle->WaitForContent(30, args.GetExtraArgs()); // 30 seconds

      return handle;
   }

   bool use_views = GuiHandler::PlatformInit();

#ifdef OS_WIN
   CefMainArgs main_args(GetModuleHandle(nullptr));
#else
   TApplication *root_app = gROOT->GetApplication();

   int cef_argc = 1;
   const char *arg2 = nullptr, *arg3 = nullptr;
   if (args.IsHeadless()) {
      // arg2 = "--allow-file-access-from-files";
      arg2 = "--disable-web-security";
      cef_argc++;
      if (use_views) {
         arg3 = "--ozone-platform=headless";
         cef_argc++;
      }
   }
   char *cef_argv[] = {root_app->Argv(0), (char *) arg2, (char *) arg3, nullptr};

   CefMainArgs main_args(cef_argc, cef_argv);
#endif

   // CEF applications have multiple sub-processes (render, plugin, GPU, etc)
   // that share the same executable. This function checks the command-line and,
   // if this is a sub-process, executes the appropriate logic.

   /*         int exit_code = CefExecuteProcess(main_args, nullptr, nullptr);
        if (exit_code >= 0) {
          // The sub-process has completed so return here.
          return exit_code;
        }
    */

   // Install xlib error handlers so that the application won't be terminated
   // on non-fatal errors.
   //         XSetErrorHandler(XErrorHandlerImpl);
   //         XSetIOErrorHandler(XIOErrorHandlerImpl);

   // Specify CEF global settings here.
   CefSettings settings;

   TString cef_main = TROOT::GetBinDir() + "/cef_main";
   cef_string_ascii_to_utf16(cef_main.Data(), cef_main.Length(), &settings.browser_subprocess_path);

#ifdef OS_LINUX
   // on linux resource directory copied to lib/
   TString path2 = TROOT::GetLibDir() + "/locales";
   cef_string_ascii_to_utf16(path2.Data(), path2.Length(), &settings.locales_dir_path);
   TString path3 = TROOT::GetLibDir();
   cef_string_ascii_to_utf16(path3.Data(), path3.Length(), &settings.resources_dir_path);
#endif

#ifdef OS_WIN
   // on windows resource directory copied to bin/
   TString path2 = TROOT::GetBinDir() + "/locales";
   cef_string_ascii_to_utf16(path2.Data(), path2.Length(), &settings.locales_dir_path);
   TString path3 = TROOT::GetBinDir();
   cef_string_ascii_to_utf16(path3.Data(), path3.Length(), &settings.resources_dir_path);
#endif

#ifdef OS_MACOSX
   // on mac there is framework directory, where resources and libs are combined together
   TString path = TROOT::GetDataDir() + "/Frameworks/Chromium Embedded Framework.framework";
   cef_string_ascii_to_utf16(path.Data(), path.Length(), &settings.framework_dir_path);
#endif

   settings.no_sandbox = true;
   // if (gROOT->IsWebDisplayBatch()) settings.single_process = true;

   // if (batch_mode)
   // settings.windowless_rendering_enabled = true;

   // settings.external_message_pump = true;
   // settings.multi_threaded_message_loop = false;

   std::string plog = "cef.log";
   cef_string_ascii_to_utf16(plog.c_str(), plog.length(), &settings.log_file);

   settings.log_severity = LOGSEVERITY_ERROR; // LOGSEVERITY_VERBOSE, LOGSEVERITY_INFO, LOGSEVERITY_WARNING,
   // LOGSEVERITY_ERROR, LOGSEVERITY_DISABLE
   // settings.uncaught_exception_stack_size = 100;
   // settings.ignore_certificate_errors = true;

   // settings.remote_debugging_port = 7890;

   // SimpleApp implements application-level callbacks for the browser process.
   // It will create the first browser instance in OnContextInitialized() after
   // CEF has initialized.
   fCefApp = new SimpleApp(use_views, args.GetHttpServer(), args.GetFullUrl(), args.GetPageContent(), args.GetWidth(), args.GetHeight(), args.IsHeadless());

   fCefApp->SetNextHandle(handle.get());

   // Initialize CEF for the browser process.
   CefInitialize(main_args, settings, fCefApp.get(), nullptr);

   if (args.IsHeadless()) {
      handle->WaitForContent(30, args.GetExtraArgs()); // 30 seconds
   } else {
      // Create timer to let run CEF message loop together with ROOT event loop
      Int_t interval = gEnv->GetValue("WebGui.CefTimer", 10);
      TCefTimer *timer = new TCefTimer((interval > 0) ? interval : 10, kTRUE);
      timer->TurnOn();
   }

   // window not yet exists here
   return handle;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// Destructor
/// Closes browser window if any

RCefWebDisplayHandle::~RCefWebDisplayHandle()
{
   fValid = kInvalid;

   CloseBrowser();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// Closes associated browser window


void RCefWebDisplayHandle::CloseBrowser()
{
   if (fBrowser) {
      auto host = fBrowser->GetHost();
      if (host) host->CloseBrowser(true);
      fBrowser = nullptr;
   }
}



///////////////////////////////////////////////////////////////////////////////////////////////////
/// Process system events until browser content is available
/// Used in headless mode for batch production like chrome --dump-dom is doing

bool RCefWebDisplayHandle::WaitForContent(int tmout_sec, const std::string &extra_args)
{
   int expired = tmout_sec * 100;
   bool did_try = false, print_finished = false;
   std::string pdffile;
   if (!extra_args.empty() && (extra_args.find("--print-to-pdf=")==0))
      pdffile = extra_args.substr(15);

   while ((--expired > 0) && GetContent().empty() && !print_finished) {

      if (gSystem->ProcessEvents()) break; // interrupted, has to return

      CefDoMessageLoopWork();

      if (fBrowser && !did_try && fBrowser->HasDocument() && !fBrowser->IsLoading() && fBrowser->GetMainFrame()) {
         did_try = true;
         if (pdffile.empty()) {
            fBrowser->GetMainFrame()->GetSource(new FrameSourceVisitor(this));
         } else {
            CefPdfPrintSettings settings;
            fBrowser->GetHost()->PrintToPDF(pdffile, settings, new HeadlessPrintCallback(&print_finished));
         }
      }

      gSystem->Sleep(10); // only 10 ms sleep
   }

   CloseBrowser();

   // call it once here to complete browser window closing, timer is not installed in batch mode
   CefDoMessageLoopWork();

   return !GetContent().empty();
}


void RCefWebDisplayHandle::AddCreator()
{
   auto &entry = FindCreator("cef");
   if (!entry)
      GetMap().emplace("cef", std::make_unique<CefCreator>());
}


struct RCefCreatorReg {
   RCefCreatorReg() { RCefWebDisplayHandle::AddCreator(); }
} newRCefCreatorReg;

