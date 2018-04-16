// Copyright (c) 2013 The Chromium Embedded Framework Authors. All rights
// reserved. Use of this source code is governed by a BSD-style license that
// can be found in the LICENSE file.

#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wshadow"

#include "include/base/cef_build.h"
#include "include/cef_app.h"

#if !defined(OS_WIN)
#include <unistd.h>
#endif

#include "include/cef_render_process_handler.h"
#include "include/base/cef_logging.h"

// Implement application-level callbacks for the browser process.
class MyRendererProcessApp : public CefApp, public CefRenderProcessHandler {

public:
   MyRendererProcessApp() : CefApp(), CefRenderProcessHandler() {}
   virtual ~MyRendererProcessApp() {}

   virtual CefRefPtr< CefRenderProcessHandler > GetRenderProcessHandler() { return this; }

   void OnRegisterCustomSchemes(CefRawPtr<CefSchemeRegistrar> registrar) OVERRIDE
   {
      // registrar->AddCustomScheme("rootscheme", true, true, true, true, true, true);
   }

   virtual void OnWebKitInitialized()
   {
      // Register the extension.
      // CefRegisterExtension("v8/jsroot", "var jsroot_cef3_identifier = true;", NULL);
   }


private:
   // Include the default reference counting implementation.
   IMPLEMENT_REFCOUNTING(MyRendererProcessApp);
   DISALLOW_COPY_AND_ASSIGN(MyRendererProcessApp);
};

// Entry point function for all processes.
int main(int argc, char *argv[])
{
   // printf("Starting CEF_MAIN ARGC %d\n", argc);
   // for (int n = 1; n < argc; n++) printf("ARGV[%d] = %s\n", n, argv[n]);

#if defined(OS_WIN)
   CefMainArgs main_args(::GetModuleHandle(NULL));
#else
   // Provide CEF with command-line arguments.
   CefMainArgs main_args(argc, argv);
#endif

   CefRefPtr<CefApp> app = new MyRendererProcessApp();

   // CEF applications have multiple sub-processes (render, plugin, GPU, etc)
   // that share the same executable. This function checks the command-line and,
   // if this is a sub-process, executes the appropriate logic.
   int exit_code = CefExecuteProcess(main_args, app, NULL);
   if (exit_code >= 0) {
      // The sub-process has completed so return here.
      return exit_code;
   }

   return 0;
}
