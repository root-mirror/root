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

#ifndef ROOT7_RCefWebDisplayHandle
#define ROOT7_RCefWebDisplayHandle

#include <ROOT/RWebDisplayHandle.hxx>

#include "simple_app.h"

/** \class RCefWebDisplayHandle
\ingroup cefwebdisplay
*/

class RCefWebDisplayHandle : public ROOT::Experimental::RWebDisplayHandle {
protected:
   class CefCreator : public Creator {

      CefRefPtr<SimpleApp> fCefApp;
   public:
      CefCreator() = default;
      virtual ~CefCreator() = default;

      std::unique_ptr<ROOT::Experimental::RWebDisplayHandle> Display(const ROOT::Experimental::RWebDisplayArgs &args) override;
   };

   enum EValidValues { kValid = 0x3C3C3C3C, kInvalid = 0x92929292 };

   unsigned fValid{kValid};  ///< used to verify if instance valid or not

   CefRefPtr<CefBrowser> fBrowser; ///< associated browser

public:
   RCefWebDisplayHandle(const std::string &url) : ROOT::Experimental::RWebDisplayHandle(url) {}

   virtual ~RCefWebDisplayHandle();

   bool IsValid() const { return fValid == kValid; }

   void SetBrowser(CefRefPtr<CefBrowser> br) { if (IsValid()) fBrowser = br; }

   void CloseBrowser();

   bool WaitForContent(int tmout_sec, const std::string &extra_args);

   static void AddCreator();

};


#endif
