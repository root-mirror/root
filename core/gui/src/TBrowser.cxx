// @(#)root/base:$Id$
// Author: Fons Rademakers   25/10/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TBrowser
\ingroup Base

Using a TBrowser one can browse all ROOT objects. It shows in a list
on the left side of the window all browsable ROOT classes. Selecting
one of the classes displays, in the icon-box on the right side, all
objects in the class. Selecting one of the objects in the icon-box,
will place all browsable objects in a new list and draws the
contents of the selected class in the icon-box. And so on....

\image html base_browser.png
*/

#include "TBrowser.h"
#include "TGuiFactory.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TTimer.h"
#include "TContextMenu.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TClass.h"
#include "TApplication.h"

/** \class TBrowserTimer
*/

class TBrowserTimer : public TTimer {

protected:
   TBrowser *fBrowser;
   Bool_t    fActivate;

public:
   TBrowserTimer(TBrowser *b, Long_t ms = 1000)
      : TTimer(ms, kTRUE), fBrowser(b), fActivate(kFALSE) { }
   Bool_t Notify();
};

/** \class TBrowserObject
This class is designed to wrap a Foreign object in order to
inject it into the Browse sub-system.
*/

class TBrowserObject : public TNamed
{

public:

   TBrowserObject(void *obj, TClass *cl, const char *brname);
   ~TBrowserObject(){;}

   void    Browse(TBrowser* b);
   Bool_t  IsFolder() const;
   TClass *IsA() const { return fClass; }

private:
   void     *fObj;   //! pointer to the foreign object
   TClass   *fClass; //! pointer to class of the foreign object

};


ClassImp(TBrowser);

////////////////////////////////////////////////////////////////////////////////
// Make sure the application environment exists and the GUI libs are loaded

Bool_t TBrowser::InitGraphics()
{
   // Make sure the application environment exists. It is need for graphics
   // (colors are initialized in the TApplication ctor).
   if (!gApplication)
      TApplication::CreateApplication();
   // make sure that the Gpad and GUI libs are loaded
   TApplication::NeedGraphicsLibs();
   if (gApplication)
      gApplication->InitializeGraphics();
   if (gROOT->IsBatch()) {
      Warning("TBrowser", "The ROOT browser cannot run in batch mode");
      return kFALSE;
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a new browser with a name, title. Width and height are by
/// default set to 640x400 and (optionally) adjusted by the screen factor
/// (depending on Rint.Canvas.UseScreenFactor to be true or false, default
/// is true).

TBrowser::TBrowser(const char *name, const char *title, TBrowserImp *extimp,
                   Option_t *opt)
   : TNamed(name, title), fLastSelectedObject(0), fImp(extimp), fTimer(0),
     fContextMenu(0), fNeedRefresh(kFALSE)
{
   if (!InitGraphics())
      return;
   if (TClass::IsCallingNew() != TClass::kRealNew) {
      fImp = 0;
   } else {
      Float_t cx = gStyle->GetScreenFactor();
      UInt_t w = UInt_t(cx*800);
      UInt_t h = UInt_t(cx*500);
      if (!fImp) fImp = gGuiFactory->CreateBrowserImp(this, title, w, h, opt);
      Create();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Create a new browser with a name, title, width and height.

TBrowser::TBrowser(const char *name, const char *title, UInt_t width,
                   UInt_t height, TBrowserImp *extimp, Option_t *opt)
   : TNamed(name, title), fLastSelectedObject(0), fImp(extimp), fTimer(0), fContextMenu(0),
     fNeedRefresh(kFALSE)
{
   if (!InitGraphics())
      return;
   if (!fImp) fImp = gGuiFactory->CreateBrowserImp(this, title, width, height, opt);
   Create();
}

////////////////////////////////////////////////////////////////////////////////
/// Create a new browser with a name, title, position, width and height.

TBrowser::TBrowser(const char *name, const char *title, Int_t x, Int_t y,
                   UInt_t width, UInt_t height, TBrowserImp *extimp, Option_t *opt)
   : TNamed(name, title), fLastSelectedObject(0), fImp(extimp), fTimer(0), fContextMenu(0),
     fNeedRefresh(kFALSE)
{
   if (!InitGraphics())
      return;
   fImp = gGuiFactory->CreateBrowserImp(this, title, x, y, width, height, opt);
   Create();
}

////////////////////////////////////////////////////////////////////////////////
/// Create a new browser with a name, title, width and height for TObject *obj.

TBrowser::TBrowser(const char *name, TObject *obj, const char *title, Option_t *opt)
   : TNamed(name, title), fLastSelectedObject(0), fImp(0), fTimer(0), fContextMenu(0),
     fNeedRefresh(kFALSE)
{
   if (!InitGraphics())
      return;
   Float_t cx = gStyle->GetScreenFactor();
   UInt_t w = UInt_t(cx*800);
   UInt_t h = UInt_t(cx*500);

   if (!fImp) fImp = gGuiFactory->CreateBrowserImp(this, title, w, h, opt);
   Create(obj);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a new browser with a name, title, width and height for TObject *obj.

TBrowser::TBrowser(const char *name, TObject *obj, const char *title,
                   UInt_t width, UInt_t height, Option_t *opt)
   : TNamed(name, title), fLastSelectedObject(0), fTimer(0), fContextMenu(0),
     fNeedRefresh(kFALSE)
{
   if (!InitGraphics())
      return;
   fImp = gGuiFactory->CreateBrowserImp(this, title, width, height, opt);
   Create(obj);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a new browser with a name, title, width and height for TObject *obj.

TBrowser::TBrowser(const char *name, TObject *obj, const char *title,
                   Int_t x, Int_t y,
                   UInt_t width, UInt_t height, Option_t *opt)
   : TNamed(name, title), fLastSelectedObject(0), fTimer(0), fContextMenu(0),
     fNeedRefresh(kFALSE)
{
   if (!InitGraphics())
      return;
   fImp = gGuiFactory->CreateBrowserImp(this, title, x, y, width, height, opt);
   Create(obj);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a new browser with a name, title, width and height for TObject *obj.

TBrowser::TBrowser(const char *name, void *obj, TClass *cl,
                   const char *objname, const char *title, Option_t *opt)
   : TNamed(name, title), fLastSelectedObject(0), fTimer(0), fContextMenu(0),
     fNeedRefresh(kFALSE)
{
   if (!InitGraphics())
      return;
   Float_t cx = gStyle->GetScreenFactor();
   UInt_t w = UInt_t(cx*800);
   UInt_t h = UInt_t(cx*500);

   fImp = gGuiFactory->CreateBrowserImp(this, title, w, h, opt);

   Create(new TBrowserObject(obj,cl,objname));
}

////////////////////////////////////////////////////////////////////////////////
/// Create a new browser with a name, title, width and height for TObject *obj.

TBrowser::TBrowser(const char *name, void *obj, TClass *cl,
                   const char *objname, const char *title,
                   UInt_t width, UInt_t height, Option_t *opt)
   : TNamed(name, title), fLastSelectedObject(0), fTimer(0), fContextMenu(0),
     fNeedRefresh(kFALSE)
{
   if (!InitGraphics())
      return;
   fImp = gGuiFactory->CreateBrowserImp(this, title, width, height, opt);
   Create(new TBrowserObject(obj,cl,objname));
}

////////////////////////////////////////////////////////////////////////////////
/// Create a new browser with a name, title, width and height for TObject *obj.

TBrowser::TBrowser(const char *name,void *obj,  TClass *cl,
                   const char *objname, const char *title,
                   Int_t x, Int_t y,
                   UInt_t width, UInt_t height, Option_t *opt)
   : TNamed(name, title), fLastSelectedObject(0), fTimer(0), fContextMenu(0),
     fNeedRefresh(kFALSE)
{
   if (!InitGraphics())
      return;
   fImp = gGuiFactory->CreateBrowserImp(this, title, x, y, width, height, opt);
   Create(new TBrowserObject(obj,cl,objname));
}

////////////////////////////////////////////////////////////////////////////////
/// Delete the browser.

TBrowser::~TBrowser()
{
   Destructor();
}

////////////////////////////////////////////////////////////////////////////////
/// Actual browser destructor.

void TBrowser::Destructor()
{
   if (fImp) fImp->CloseTabs();
   R__LOCKGUARD(gROOTMutex);
   gROOT->GetListOfBrowsers()->Remove(this);
   delete fContextMenu;
   delete fTimer;
   if (fImp) delete fImp;
}

////////////////////////////////////////////////////////////////////////////////
/// Add object with name to browser. If name not set the objects GetName()
/// is used. If check < 0 (default) no check box is drawn, if 0 then
/// unchecked checkbox is added, if 1 checked checkbox is added.

void TBrowser::Add(TObject *obj, const char *name, Int_t check)
{
   if (obj && fImp) {
      fImp->Add(obj, name, check);
      obj->SetBit(kMustCleanup);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Add foreign object with name to browser.
/// 'cl' is the type use to store the value of obj.
/// So literally the following pseudo code should be correct:
///~~~ {.cpp}
///    `cl->GetName()` * ptr = (`cl->GetName()`*) obj;
///~~~
/// and the value of obj is not necessarily the start of the object.
/// If check < 0 (default) no check box is drawn, if 0 then
/// unchecked checkbox is added, if 1 checked checkbox is added.

void TBrowser::Add(void *obj, TClass *cl, const char *name, Int_t check)
{
   if (!obj || !cl) return;
   TObject *to;
   if (cl->IsTObject()) to = (TObject*)cl->DynamicCast(TObject::Class(),obj,kTRUE);
   else                 to = new TBrowserObject(obj,cl,name);

   if (!to) return;
   Add(to,name,check);
}

////////////////////////////////////////////////////////////////////////////////
/// Add checkbox for this item.

void TBrowser::AddCheckBox(TObject *obj, Bool_t check)
{
   if (obj && fImp) {
      fImp->AddCheckBox(obj, check);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Change status of checkbox for this item.

void TBrowser::CheckObjectItem(TObject *obj, Bool_t check)
{
   if (obj && fImp) {
      fImp->CheckObjectItem(obj, check);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Remove checkbox for this item.

void TBrowser::RemoveCheckBox(TObject *obj)
{
   if (obj && fImp) {
      fImp->RemoveCheckBox(obj);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Create the browser, called by the ctors.

void TBrowser::Create(TObject *obj)
{
   fNeedRefresh = kFALSE;

   fTimer = new TBrowserTimer(this);
   gSystem->AddTimer(fTimer);

   R__LOCKGUARD(gROOTMutex);
   gROOT->GetListOfBrowsers()->Add(this);

   // Get the list of globals
   gROOT->GetListOfGlobals(kTRUE);
   gROOT->GetListOfGlobalFunctions(kTRUE);

   fContextMenu = new TContextMenu("BrowserContextMenu") ;

   // Fill the first list from the present TObject obj
   if (obj) {
      Add(obj);
      if (fImp) fImp->BrowseObj(obj);
   } else if (fImp) {
      // Fill the first list with all browsable classes from TROOT
      fImp->BrowseObj(gROOT);
   }

   // The first list will be filled by TWin32BrowserImp ctor
   // with all browsable classes from TROOT
}

////////////////////////////////////////////////////////////////////////////////
/// Execute default action for selected object (action is specified
/// in the `$HOME/.root.mimes` or `$ROOTSYS/etc/root.mimes file`).

void TBrowser::ExecuteDefaultAction(TObject *obj)
{
   if (obj && fImp)
      fImp->ExecuteDefaultAction(obj);
}

////////////////////////////////////////////////////////////////////////////////
/// Recursively remove obj from browser.

void TBrowser::RecursiveRemove(TObject *obj)
{
   if (fImp && obj) {
      fImp->RecursiveRemove(obj);
      fNeedRefresh = kTRUE;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Refresh browser contents.

void TBrowser::Refresh()
{
   fNeedRefresh = kTRUE;
   if (fImp) fImp->Refresh();
   fNeedRefresh = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Assign the last selected object.

void TBrowser::SetSelected(TObject *clickedObject)
{
   fLastSelectedObject = clickedObject;
}

/** \class  TBrowserTimer
Called whenever timer times out.
*/

Bool_t TBrowserTimer::Notify()
{
   if (fBrowser) {
      if (fBrowser->GetRefreshFlag()) {
         fBrowser->SetRefreshFlag(kFALSE);
         fActivate = kTRUE;
      } else if (fActivate) {
         fActivate = kFALSE;
         fBrowser->Refresh();
      }
   }
   Reset();

   return kFALSE;
}

/** \class TBrowserObject
This is a wrapper class to emulate the TObject interface
around an object of a non-TObject class
*/

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TBrowserObject::TBrowserObject(void *obj, TClass *cl, const char *brname)
   : TNamed(brname, cl ? cl->GetName() : ""), fObj(obj), fClass(cl)
{
   if (cl==0) Fatal("Constructor","Class parameter should not be null");
   SetBit(kCanDelete);
}

////////////////////////////////////////////////////////////////////////////////
/// Return kTRUE if the object is a folder (contains browsable objects).

Bool_t TBrowserObject::IsFolder() const
{
   return fClass->IsFolder(fObj);
}

////////////////////////////////////////////////////////////////////////////////
/// Browse the content of the underlying object.

void TBrowserObject::Browse(TBrowser* b)
{
   fClass->Browse(fObj, b);
}
