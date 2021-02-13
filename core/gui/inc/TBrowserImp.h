// @(#)root/base:$Id$
// Author: Fons Rademakers   15/11/95

/*************************************************************************
 * Copyright (C) 1995-2021, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


#ifndef ROOT_TBrowserImp
#define ROOT_TBrowserImp

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TBrowserImp                                                          //
//                                                                      //
// ABC describing GUI independent browser implementation protocol.      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TObject.h"

class TBrowser;
class TGMainFrame;

class TBrowserImp {

protected:
   TBrowser  *fBrowser{nullptr};     ///< TBrowser associated with this implementation
   Bool_t     fShowCycles{kFALSE};   ///< Show object cycle numbers in browser

   TBrowserImp(const TBrowserImp&) = delete;
   TBrowserImp &operator=(const TBrowserImp& br) = delete;

public:
   TBrowserImp(TBrowser *b = nullptr);
   TBrowserImp(TBrowser *b, const char *title, UInt_t width, UInt_t height, Option_t *opt = "");
   TBrowserImp(TBrowser *b, const char *title, Int_t x, Int_t y, UInt_t width, UInt_t height, Option_t *opt = "");
   virtual ~TBrowserImp() = default;

   virtual void      Add(TObject *, const char *, Int_t) { }
   virtual void      AddCheckBox(TObject *, Bool_t = kFALSE) { }
   virtual void      CheckObjectItem(TObject *, Bool_t = kFALSE) { }
   virtual void      RemoveCheckBox(TObject *) { }
   virtual void      BrowseObj(TObject *) { }
   TBrowser         *Browser() const { return fBrowser; }
   virtual void      CloseTabs() { }
   virtual void      ExecuteDefaultAction(TObject *) { }
   virtual void      Iconify() { }
   virtual void      RecursiveRemove(TObject *) { }
   virtual void      Refresh(Bool_t = kFALSE) { }
   virtual void      Show() { }
   virtual void      SetDrawOption(Option_t * = "") { }
   virtual Option_t *GetDrawOption() const { return nullptr; }

   virtual Long_t    ExecPlugin(const char *, const char *, const char *, Int_t, Int_t) { return 0; }
   virtual void      SetStatusText(const char *, Int_t) { }
   virtual void      StartEmbedding(Int_t, Int_t) { }
   virtual void      StopEmbedding(const char *) { }

   virtual TGMainFrame *GetMainFrame() const { return nullptr; }

   virtual TBrowser *GetBrowser() const      { return fBrowser; }
   virtual void      SetBrowser(TBrowser *b) { fBrowser = b; }

   ClassDef(TBrowserImp,0)  //ABC describing browser implementation protocol
};

#endif
