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

#include "gui_handler.h"

#include <windows.h>
#include <string>

#include "include/cef_browser.h"
#include "include/cef_config.h"

bool GuiHandler::PlatformInit()
{
#ifdef CEF_X11
   return false; // compiled without ozone support
#else
   return true; // compiled with ozone support
#endif
}

void GuiHandler::PlatformTitleChange(CefRefPtr<CefBrowser> browser, const CefString &title)
{
   CefWindowHandle hwnd = browser->GetHost()->GetWindowHandle();
   SetWindowText(hwnd, std::string(title).c_str());
}
