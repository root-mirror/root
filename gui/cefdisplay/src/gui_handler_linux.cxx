// Author: Sergey Linev <S.Linev@gsi.de>
// Date: 2017-06-29
// Warning: This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!

/*************************************************************************
 * Copyright (C) 1995-2017, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#if !defined(_MSC_VER)
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wshadow"
#endif

#include "gui_handler.h"

#include "include/cef_config.h"

#ifdef CEF_X11

#include <X11/Xatom.h>
#include <X11/Xlib.h>
#include <string>

#include "include/base/cef_logging.h"
#include "include/cef_browser.h"

int x11_errhandler( Display *dpy, XErrorEvent *err )
{
  // special for modality usage: XGetWindowProperty + XQueryTree()
  if (err->error_code == BadWindow) {
     // if ( err->request_code == 25 && qt_xdnd_handle_badwindow() )
     return 0;
  } else if (err->error_code == BadMatch && err->request_code == 42) {
     //  special case for  X_SetInputFocus
     return 0;
  } else if (err->error_code == BadDrawable && err->request_code == 14) {
     return 0;
  }

  // here XError are forwarded
  char errstr[512];
  XGetErrorText( dpy, err->error_code, errstr, sizeof(errstr) );
  printf( "X11 Error: %d opcode: %d info: %s\n", err->error_code, err->request_code, errstr );
  return 0;
}

bool GuiHandler::PlatformInit()
{
   // install custom X11 error handler to avoid application exit in case of X11 failure
   XSetErrorHandler( x11_errhandler );

   return false; // do not use view framework
}



void GuiHandler::PlatformTitleChange(CefRefPtr<CefBrowser> browser, const CefString &title)
{

   std::string titleStr(title);

   // Retrieve the X11 display shared with Chromium.
   ::Display *display = cef_get_xdisplay();
   DCHECK(display);

   // Retrieve the X11 window handle for the browser.
   ::Window window = browser->GetHost()->GetWindowHandle();
   DCHECK(window != kNullWindowHandle);

   // Retrieve the atoms required by the below XChangeProperty call.
   const char *kAtoms[] = {"_NET_WM_NAME", "UTF8_STRING"};
   Atom atoms[2];
   int result = XInternAtoms(display, const_cast<char **>(kAtoms), 2, false, atoms);
   if (!result) NOTREACHED();

   // Set the window title.
   XChangeProperty(display, window, atoms[0], atoms[1], 8, PropModeReplace,
                   reinterpret_cast<const unsigned char *>(titleStr.c_str()), titleStr.size());

   // TODO(erg): This is technically wrong. So XStoreName and friends expect
   // this in Host Portable Character Encoding instead of UTF-8, which I believe
   // is Compound Text. This shouldn't matter 90% of the time since this is the
   // fallback to the UTF8 property above.
   XStoreName(display, browser->GetHost()->GetWindowHandle(), titleStr.c_str());
}

#else

bool GuiHandler::PlatformInit()
{
   return true; // use view framework
}

void GuiHandler::PlatformTitleChange(CefRefPtr<CefBrowser>, const CefString &)
{
   // do nothing
}


#endif

