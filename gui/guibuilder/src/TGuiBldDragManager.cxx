// @(#)root/guibuilder:$Id: 99541fed972db7e4279fadcd4c0a0f0ae28a0a4d $
// Author: Valeriy Onuchin   12/09/04

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


#include "TGuiBldDragManager.h"
#include "TGuiBldEditor.h"
#include "TRootGuiBuilder.h"

#include "TTimer.h"
#include "TList.h"
#include "TClass.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TColor.h"
#include "TImage.h"

#include "TError.h"
#include "TClassMenuItem.h"
#include "TMethod.h"
#include "TMethodCall.h"
#include "TBaseClass.h"
#include "TDataMember.h"
#include "TMethodArg.h"
#include "TToggle.h"
#include "TDataType.h"
#include "TObjString.h"
#include "TInterpreter.h"

#include "KeySymbols.h"
#include "TGResourcePool.h"
#include "TGMenu.h"
#include "TGFileDialog.h"
#include "TGMsgBox.h"
#include "TRandom.h"
#include "TGButton.h"
#include "TGMdi.h"
#include "TGTextEntry.h"
#include "TGDockableFrame.h"
#include "TGColorDialog.h"
#include "TGFontDialog.h"
#include "TGComboBox.h"
#include "TGCanvas.h"
#include "TGLabel.h"
#include "TGProgressBar.h"
#include "TGScrollBar.h"
#include "TVirtualX.h"
#include "strlcpy.h"
#include "snprintf.h"

#undef DEBUG_LOCAL

//_____________________________________________________________________________
//
// TGuiBldDragManager
//
// Drag and drop manager used by the ROOT GUI Builder.
//_____________________________________________________________________________

ClassImp(TGuiBldDragManager);

static UInt_t gGridStep = 8;
static TGuiBldDragManager *gGuiBldDragManager = 0;
TGColorDialog *TGuiBldDragManager::fgGlobalColorDialog = 0;
TGFontDialog *TGuiBldDragManager::fgGlobalFontDialog = 0;

static const char *gSaveMacroTypes[] = {
   "Macro files", "*.C",
   "All files",   "*",
   0,             0
};

static const char *gImageTypes[] = {
   "All files", "*",
   "XPM",     "*.xpm",
   "GIF",     "*.gif",
   "PNG",     "*.png",
   "JPEG",    "*.jpg",
   "TARGA",   "*.tga",
   "BMP",     "*.bmp",
   "ICO",     "*.ico",
   "XCF",     "*.xcf",
   "CURSORS", "*.cur",
   "PPM",     "*.ppm",
   "PNM",     "*.pnm",
   "XBM",     "*.xbm",
   "TIFF",    "*.tiff",
   "Enacapsulated PostScript", "*.eps",
   "PostScript", "*.ps",
   "PDF",        "*.pdf",
   "ASImage XML","*.xml",
    0,             0
};


////////////////////////////////////////////////////////////////////////////////


class TGuiBldMenuDialog : public TGTransientFrame {

friend class TGuiBldDragManager;

public:
   TGButton       *fOK;       // OK button
   //TGButton     *fApply;    // apply button
   TGButton       *fCancel;   // cancel button
   TObject        *fObject;   // selected object/frame
   TMethod        *fMethod;   // method to be applied
   TGLayoutHints  *fL1;       // internally used layout hints
   TGLayoutHints  *fL2;       // internally used layout hints
   TList          *fWidgets;  // list of widgets

public:
   virtual ~TGuiBldMenuDialog();
   TGuiBldMenuDialog(const TGWindow *main, TObject *obj, TMethod *method);

   const char *GetParameters();
   void CloseWindow();
   void ConnectButtonSignals();
   void Build();
   void Popup();
   void ApplyMethod();
   void Add(const char *argname, const char *value, const char *type);

};

static TGuiBldMenuDialog *gMenuDialog = 0;


////////////////////////////////////////////////////////////////////////////////
/// ctor.

TGuiBldMenuDialog::TGuiBldMenuDialog(const TGWindow *main, TObject *obj, TMethod *method) :
    TGTransientFrame(gClient->GetDefaultRoot(), main, 200, 100)
{
   fObject = obj;
   fMethod = method;
   if (!obj) return; // zombie

   fWidgets = new TList();

   fL1 = new TGLayoutHints(kLHintsTop | kLHintsCenterX, 0, 0, 5, 0);
   fL2 = new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5);

   TString title = obj->ClassName();
   title += "::";
   title += method->GetName();

   Build();
   ConnectButtonSignals();

   SetWindowName(title);
   SetIconName(title);
   SetEditDisabled(kEditDisable);

   //TRootGuiBuilder::PropagateBgndColor(this, TRootGuiBuilder::GetBgnd());
}

////////////////////////////////////////////////////////////////////////////////
/// dtor.

TGuiBldMenuDialog::~TGuiBldMenuDialog()
{
   fWidgets->Delete();
   delete fWidgets;
   delete fL1;
   delete fL2;
}

////////////////////////////////////////////////////////////////////////////////
/// Connect buttons signals

void TGuiBldMenuDialog::ConnectButtonSignals()
{
   fOK->Connect("Pressed()", "TGuiBldDragManager", gGuiBldDragManager, "DoDialogOK()");
   //fApply->Connect("Pressed()", "TGuiBldDragManager", gGuiBldDragManager, "DoDialogApply()");
   fCancel->Connect("Pressed()", "TGuiBldDragManager", gGuiBldDragManager, "DoDialogCancel()");
}

////////////////////////////////////////////////////////////////////////////////
/// execute method for obejct with input args

void TGuiBldMenuDialog::ApplyMethod()
{
   const char *params = GetParameters();
   fObject->Execute(fMethod->GetName(), params);
}

////////////////////////////////////////////////////////////////////////////////
/// Return input parameters as single string.

const char *TGuiBldMenuDialog::GetParameters()
{
   static char params[1024];
   char param[256];

   TObjString   *str;
   TObject      *obj;

   Int_t selfobjpos;
//   if (fMenu->GetContextMenu()->GetSelectedMenuItem())
//      selfobjpos = fMenu->GetContextMenu()->GetSelectedMenuItem()->GetSelfObjectPos();
//   else
   selfobjpos = -1;

   params[0] = 0;
   TIter next(fWidgets);
   Int_t nparam = 0;

   while ((obj = next())) {        // first element is label, skip...
      if (obj->IsA() != TGLabel::Class()) break;
      obj = next();                // get either TGTextEntry or TGComboBox
      str = (TObjString *) next(); // get type string

      nparam++;

      const char *type = str->GetString().Data();
      const char *data = 0;

      if (obj->IsA() == TGTextEntry::Class())
         data = ((TGTextEntry *) obj)->GetBuffer()->GetString();

      // TODO: Combobox...

      // if necessary, replace the selected object by it's address
      if (selfobjpos == nparam-1) {
         if (params[0]) strlcat(params, ",", 1024-strlen(params));
         snprintf(param, 255, "(TObject*)0x%lx", (Long_t)fObject);
         strlcat(params, param, 1024-strlen(params));
      }

      if (params[0]) strlcat(params, ",", 1024-strlen(params));
      if (data) {
         if (!strncmp(type, "char*", 5))
            snprintf(param, 255, "\"%s\"", data);
         else
            strlcpy(param, data, sizeof(param));
      } else
         strlcpy(param, "0", sizeof(param));

      strlcat(params, param, 1024-strlen(params));
   }

   // if selected object is the last argument, have to insert it here
   if (selfobjpos == nparam) {
      if (params[0]) strlcat(params, ",", 1024-strlen(params));
      snprintf(param, 255, "(TObject*)0x%lx", (Long_t)fObject);
      strlcat(params, param, 1024-strlen(params));
   }

   return params;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a string describing method argument.

static TString CreateArgumentTitle(TMethodArg *argument)
{
   static TString ret;

   if (argument) {
      ret.Form("(%s)  %s", argument->GetTitle(), argument->GetName());
      if (argument->GetDefault() && *(argument->GetDefault())) {
         ret += "  [default: ";
         ret += argument->GetDefault();
         ret += "]";
      }
   }

   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// Add a label and text input field.

void TGuiBldMenuDialog::Add(const char *argname, const char *value, const char *type)
{
   TGLabel      *l = new TGLabel(this, argname);
   TGTextBuffer *b = new TGTextBuffer(20);
   b->AddText(0, value);
   TGTextEntry  *t = new TGTextEntry(this, b);

   t->Connect("ReturnPressed()", "TGuiBldDragManager", gGuiBldDragManager, "DoDialogOK()");
   t->Resize(260, t->GetDefaultHeight());
   AddFrame(l, fL1);
   AddFrame(t, fL2);

   fWidgets->Add(l);
   fWidgets->Add(t);
   fWidgets->Add(new TObjString(type));
}

////////////////////////////////////////////////////////////////////////////////
/// Close window

void TGuiBldMenuDialog::CloseWindow()
{
   gGuiBldDragManager->DoDialogCancel();
}

////////////////////////////////////////////////////////////////////////////////
/// Build dialog

void TGuiBldMenuDialog::Build()
{
   TMethodArg *argument = 0;
   Int_t selfobjpos = -1;

   TIter next(fMethod->GetListOfMethodArgs());
   Int_t argpos = 0;

   while ((argument = (TMethodArg *) next())) {
      // Do not input argument for self object
      if (selfobjpos != argpos) {
         TString arg_name = CreateArgumentTitle(argument);
         const char *argname = arg_name.Data();
         const char *type = argument->GetTypeName();
         TDataType  *datatype = gROOT->GetType(type);
         const char *charstar = "char*";
         char basictype[32];

         if (datatype) {
            strlcpy(basictype, datatype->GetTypeName(), sizeof(basictype));
         } else {
            TClass *cl = TClass::GetClass(type);
            if (strncmp(type, "enum", 4) && (cl && !(cl->Property() & kIsEnum)))
               Warning("Dialog", "data type is not basic type, assuming (int)");
            strlcpy(basictype, "int", sizeof(basictype));
         }

         if (strchr(argname, '*')) {
            strlcat(basictype, "*", 32-strlen(basictype));
            type = charstar;
         }

         TDataMember *m = argument->GetDataMember();
         if (m && m->GetterMethod(fObject->IsA())) {

            // Get the current value and form it as a text:
            char val[256];

            if (!strncmp(basictype, "char*", 5)) {
               char *tdefval = 0;
               m->GetterMethod()->Execute(fObject, "", &tdefval);
               if (tdefval && strlen(tdefval))
                  strlcpy(val, tdefval, sizeof(val));
            } else if (!strncmp(basictype, "float", 5) ||
                       !strncmp(basictype, "double", 6)) {
               Double_t ddefval = 0.0;
               m->GetterMethod()->Execute(fObject, "", ddefval);
               snprintf(val, 255, "%g", ddefval);
            } else if (!strncmp(basictype, "char", 4) ||
                       !strncmp(basictype, "bool", 4) ||
                       !strncmp(basictype, "int", 3)  ||
                       !strncmp(basictype, "long", 4) ||
                       !strncmp(basictype, "short", 5)) {
               Long_t ldefval = 0L;
               m->GetterMethod()->Execute(fObject, "", ldefval);
               snprintf(val, 255, "%li", ldefval);
            }

            // Find out whether we have options ...

            TList *opt;
            // coverity[returned_pointer]: keep for later use
            if ((opt = m->GetOptions())) {
               Warning("Dialog", "option menu not yet implemented");

            } else {
               // we haven't got options - textfield ...
               Add(argname, val, type);
            }
         } else {    // if m not found ...

            char val[256] = "";
            const char *tval = argument->GetDefault();
            if (tval) strlcpy(val, tval, sizeof(val));
            Add(argname, val, type);
         }
      }
      argpos++;
   }

   // add OK, Apply, Cancel buttons
   TGHorizontalFrame *hf = new TGHorizontalFrame(this, 60, 20, kFixedWidth);
   TGLayoutHints     *l1 = new TGLayoutHints(kLHintsCenterY | kLHintsExpandX, 5, 5, 0, 0);
   UInt_t  width = 0, height = 0;

   fWidgets->Add(l1);

   fOK = new TGTextButton(hf, "&OK", 1);
   hf->AddFrame(fOK, l1);
   fWidgets->Add(fOK);
   height = fOK->GetDefaultHeight();
   width  = TMath::Max(width, fOK->GetDefaultWidth());

/*
   fApply = new TGTextButton(hf, "&Apply", 2);
   hf->AddFrame(fApply, l1);
   fWidgets->Add(fApply);
   height = fApply->GetDefaultHeight();
   width  = TMath::Max(width, fApply->GetDefaultWidth());
*/

   fCancel = new TGTextButton(hf, "&Cancel", 3);
   hf->AddFrame(fCancel, l1);
   fWidgets->Add(fCancel);
   height = fCancel->GetDefaultHeight();
   width  = TMath::Max(width, fCancel->GetDefaultWidth());

   // place buttons at the bottom
   l1 = new TGLayoutHints(kLHintsBottom | kLHintsCenterX, 0, 0, 5, 5);
   AddFrame(hf, l1);
   fWidgets->Add(l1);
   fWidgets->Add(hf);

   hf->Resize((width + 20) * 3, height);

   // map all widgets and calculate size of dialog
   MapSubwindows();
}

////////////////////////////////////////////////////////////////////////////////
/// Popup dialog.

void TGuiBldMenuDialog::Popup()
{
   UInt_t width  = GetDefaultWidth();
   UInt_t height = GetDefaultHeight();

   Resize(width, height);

   Window_t wdummy;
   Int_t x = (Int_t)((TGFrame*)fMain)->GetWidth();
   Int_t y = (Int_t)((TGFrame*)fMain)->GetHeight();
   gVirtualX->TranslateCoordinates(fMain->GetId(), fClient->GetDefaultRoot()->GetId(),
                                   x, y, x, y, wdummy);

   x += 10;
   y += 10;

   // make the message box non-resizable
   SetWMSize(width, height);
   SetWMSizeHints(width, height, width, height, 0, 0);

   SetMWMHints(kMWMDecorAll | kMWMDecorResizeH  | kMWMDecorMaximize |
                              kMWMDecorMinimize | kMWMDecorMenu,
               kMWMFuncAll  | kMWMFuncResize    | kMWMFuncMaximize |
                              kMWMFuncMinimize,
               kMWMInputModeless);

   Move(x, y);
   SetWMPosition(x, y);
   MapRaised();
   fClient->WaitFor(this);
}


///////////////////////// auxilary static functions ///////////////////////////
////////////////////////////////////////////////////////////////////////////////
/// Helper. Return a window located at point x,y (in screen coordinates)

static Window_t GetWindowFromPoint(Int_t x, Int_t y)
{
   Window_t src, dst, child;
   Window_t ret = 0;
   Int_t xx = x;
   Int_t yy = y;

   if (!gGuiBldDragManager || gGuiBldDragManager->IsStopped() ||
       !gClient->IsEditable()) return 0;

   dst = src = child = gVirtualX->GetDefaultRootWindow();

   while (child && dst) {
      src = dst;
      dst = child;
      gVirtualX->TranslateCoordinates(src, dst, xx, yy, xx, yy, child);
      ret = dst;
   }
   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// Helper to layout

static void layoutFrame(TGFrame *frame)
{
   if (!frame || !frame->InheritsFrom(TGCompositeFrame::Class())) {
      return;
   }

   TGCompositeFrame *comp = (TGCompositeFrame*)frame;

   if (comp->GetLayoutManager()) {
      comp->GetLayoutManager()->Layout();
   } else {
      comp->Layout();
   }
   gClient->NeedRedraw(comp);

   TIter next(comp->GetList());
   TGFrameElement *fe;

   while ((fe = (TGFrameElement*)next())) {
      layoutFrame(fe->fFrame);
      gClient->NeedRedraw(fe->fFrame);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Our own error handler (not used yet)

static void GuiBldErrorHandler(Int_t /*level*/, Bool_t /*abort*/,
                                 const char * /*location*/, const char * /*msg*/)
{
}

////////////////////////////////////////////////////////////////////////////////
class TGuiBldDragManagerGrid {

public:
   static UInt_t   fgStep;
   static ULong_t  fgPixel;
   static TGGC    *fgBgnd;

   Pixmap_t    fPixmap;
   TGWindow   *fWindow;
   Int_t       fWinId;

   TGuiBldDragManagerGrid();
   ~TGuiBldDragManagerGrid();
   void  Draw();
   void  SetStep(UInt_t step);
   void  InitPixmap();
   void  InitBgnd();
};

UInt_t   TGuiBldDragManagerGrid::fgStep = gGridStep;
ULong_t  TGuiBldDragManagerGrid::fgPixel = 0;
TGGC    *TGuiBldDragManagerGrid::fgBgnd = 0;

////////////////////////////////////////////////////////////////////////////////
/// Create a grid background for the selected window

TGuiBldDragManagerGrid::TGuiBldDragManagerGrid()
{
   fPixmap = 0;
   fWindow = 0;
   fWinId = 0;

   if (!fgBgnd) {
      InitBgnd();
   }
   SetStep(fgStep);
}

////////////////////////////////////////////////////////////////////////////////
/// ctor.

TGuiBldDragManagerGrid::~TGuiBldDragManagerGrid()
{
   fWindow = gClient->GetWindowById(fWinId);

   if (fWindow) {
      fWindow->SetBackgroundPixmap(0);
      fWindow->SetBackgroundColor(((TGFrame*)fWindow)->GetBackground());
      gClient->NeedRedraw(fWindow, kTRUE);
   }
   if (fPixmap) {
      gVirtualX->DeletePixmap(fPixmap);
   }

   fPixmap = 0;
   fWindow = 0;
   fWinId = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the grid step

void TGuiBldDragManagerGrid::SetStep(UInt_t step)
{
   if (!gClient || !gClient->IsEditable()) {
      return;
   }

   fWindow = (TGWindow*)gClient->GetRoot();
   fWinId = fWindow->GetId();
   fgStep = step;
   InitPixmap();

}

////////////////////////////////////////////////////////////////////////////////
/// Create grid background.

void TGuiBldDragManagerGrid::InitBgnd()
{
   if (fgBgnd) {
      return;
   }

   fgBgnd = new TGGC(TGFrame::GetBckgndGC());

   Float_t r, g, b;

   r = 232./255;
   g = 232./255;
   b = 226./255;

   fgPixel = TColor::RGB2Pixel(r, g, b);
   fgBgnd->SetForeground(fgPixel);
}

////////////////////////////////////////////////////////////////////////////////
/// Create grid background pixmap

void TGuiBldDragManagerGrid::InitPixmap()
{
   if (fPixmap) {
      gVirtualX->DeletePixmap(fPixmap);
   }

   fPixmap = gVirtualX->CreatePixmap(gClient->GetDefaultRoot()->GetId(), fgStep, fgStep);
   gVirtualX->FillRectangle(fPixmap, fgBgnd->GetGC(), 0, 0, fgStep, fgStep);

   if(fgStep > 2) {
      gVirtualX->FillRectangle(fPixmap, TGFrame::GetShadowGC()(),
                               fgStep - 1, fgStep - 1, 1, 1);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw grid over editted frame.

void TGuiBldDragManagerGrid::Draw()
{
   if (!gClient || !gClient->IsEditable()) {
      return;
   }

   fWindow = gClient->GetWindowById(fWinId);

   if (fWindow && (fWindow != gClient->GetRoot())) {
      fWindow->SetBackgroundPixmap(0);
      fWindow->SetBackgroundColor(((TGFrame*)fWindow)->GetBackground());
      gClient->NeedRedraw(fWindow);
   }

   if (!fPixmap) {
      InitPixmap();
   }

   fWindow = (TGWindow*)gClient->GetRoot();
   fWinId = fWindow->GetId();
   fWindow->SetBackgroundPixmap(fPixmap);

   gClient->NeedRedraw(fWindow);
}


////////////////////////////////////////////////////////////////////////////////
class TGuiBldDragManagerRepeatTimer : public TTimer {

private:
   TGuiBldDragManager *fManager;   // back pointer

public:
   TGuiBldDragManagerRepeatTimer(TGuiBldDragManager *m, Long_t ms) :
                                 TTimer(ms, kTRUE) { fManager = m; }
   Bool_t Notify() { fManager->HandleTimer(this); Reset(); return kFALSE; }
};


////////////////////////////////////////////////////////////////////////////////
class TGGrabRect : public TGFrame {

private:
   Pixmap_t    fPixmap;
   ECursor     fType;

public:
   TGGrabRect(Int_t type);
   ~TGGrabRect() {}

   Bool_t HandleButton(Event_t *ev);
   ECursor GetType() const { return fType; }
};

////////////////////////////////////////////////////////////////////////////////
/// ctor.

TGGrabRect::TGGrabRect(Int_t type) :
      TGFrame(gClient->GetDefaultRoot(), 8, 8, kTempFrame)
{
   fType = kTopLeft;

   switch (type) {
      case 0:
         fType = kTopLeft;
         break;
      case 1:
         fType = kTopSide;
         break;
      case 2:
         fType = kTopRight;
         break;
      case 3:
         fType = kBottomLeft;
         break;
      case 4:
         fType = kLeftSide;
         break;
      case 5:
         fType = kRightSide;
         break;
      case 6:
         fType = kBottomSide;
         break;
      case 7:
         fType = kBottomRight;
         break;
   }

   SetWindowAttributes_t attr;
   attr.fMask             = kWAOverrideRedirect | kWASaveUnder;
   attr.fOverrideRedirect = kTRUE;
   attr.fSaveUnder        = kTRUE;

   gVirtualX->ChangeWindowAttributes(fId, &attr);

   fPixmap = gVirtualX->CreatePixmap(gVirtualX->GetDefaultRootWindow(), 8, 8);
   const TGGC *bgc = TRootGuiBuilder::GetPopupHlghtGC();
//fClient->GetResourcePool()->GetSelectedBckgndGC();
   TGGC *gc = new TGGC(TGFrame::GetBckgndGC());

   Pixel_t back;
   fClient->GetColorByName("black", back);
   gc->SetBackground(back);
   gc->SetForeground(back);

   gVirtualX->FillRectangle(fPixmap, bgc->GetGC(), 0, 0, 7, 7);
   gVirtualX->DrawRectangle(fPixmap, gc->GetGC(), 0, 0, 7, 7);
//   gVirtualX->DrawRectangle(fPixmap, fClient->GetResourcePool()->GetSelectedBckgndGC()->GetGC(), 1, 1, 5, 5);

   AddInput(kButtonPressMask);
   SetBackgroundPixmap(fPixmap);

   gVirtualX->SetCursor(fId, gVirtualX->CreateCursor(fType));
}

////////////////////////////////////////////////////////////////////////////////
/// Handle button press event.

Bool_t TGGrabRect::HandleButton(Event_t *ev)
{
   gGuiBldDragManager->CheckDragResize(ev);
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
class TGAroundFrame : public TGFrame {

public:
   TGAroundFrame();
   ~TGAroundFrame() {}
};

////////////////////////////////////////////////////////////////////////////////
/// ctor.

TGAroundFrame::TGAroundFrame() : TGFrame(gClient->GetDefaultRoot(), 1, 1,
                                         kTempFrame | kOwnBackground)
{
   SetWindowAttributes_t attr;
   attr.fMask             = kWAOverrideRedirect | kWASaveUnder;
   attr.fOverrideRedirect = kTRUE;
   attr.fSaveUnder        = kTRUE;

   gVirtualX->ChangeWindowAttributes(fId, &attr);
   ULong_t blue;
   fClient->GetColorByName("blue", blue);
   SetBackgroundColor(blue);
}


////////////////////////////////////////////////////////////////////////////////
class TGuiBldDragManagerPimpl {

friend class TGuiBldDragManager;

private:
   TGuiBldDragManager     *fManager;   // drag and drop manager
   TTimer            *fRepeatTimer;    // repeat rate timer (when mouse stays pressed)
   TGFrame           *fGrab;           // grabbed/selected frame
   TGLayoutHints     *fGrabLayout;     // layout of grabbed frame
   TGFrame           *fSaveGrab;       // used during context menu handling
   TGFrame           *fClickFrame;     // last clicked frame
   TGuiBldDragManagerGrid *fGrid;      //
   ECursor            fResizeType;     // defines resize type
   Int_t              fX0, fY0;        // initial drag position in pixels
   Int_t              fX, fY;          // current drag position in pixels
   Int_t              fXf, fYf;        // offset of inititial position inside frame
   Int_t              fGrabX, fGrabY;  //
   const TGWindow    *fGrabParent;     // parent of the grabbed/seleceted frame
   Int_t              fLastPopupAction;//
   Bool_t             fReplaceOn;
   TGGrabRect        *fGrabRect[8];    // small rectangles drawn over grabbed/selected frame
   TGFrame           *fAroundFrame[4]; // red lines drawn over layouted frame
   Bool_t             fGrabRectHidden;
   TGFrameElement    *fGrabListPosition;
   Bool_t             fButtonPressed;  //
   Bool_t             fCompacted;   //
   TGFrame           *fPlane; // highlighted plain composite frame when mose is moving
   TGFrame           *fSpacePressedFrame; // frame which was grabbed via spacebar pressed
   Bool_t             fPlacePopup; // kTRUE is menu fo frame was placed
   TList             *fFrameMenuTrash; // trash list
   TGFrame           *fMenuObject;     // object/frame for which context menu is created

public:
   TGuiBldDragManagerPimpl(TGuiBldDragManager *m) {
      fManager = m;
      fRepeatTimer = new TGuiBldDragManagerRepeatTimer(m, 100);

      int i = 0;
      for (i = 0; i <8; i++) {
         fGrabRect[i] = new TGGrabRect(i);
      }
      for (i = 0; i <4; i++) {
         fAroundFrame[i] = new TGAroundFrame();
      }

      fFrameMenuTrash = new TList();

      ResetParams();
   }
   void ResetParams() {
      fGrab = 0;
      fSaveGrab = 0;
      fClickFrame = 0;
      fGrid = 0;
      fX0 = fY0 = fX = fY = fXf = fYf = fGrabX = fGrabY = 0;
      fGrabParent = 0;
      fResizeType = kPointer;
      fLastPopupAction = kNoneAct;
      fReplaceOn = kFALSE;
      fGrabLayout = 0;
      fGrabRectHidden = kFALSE;
      fGrabListPosition = 0;
      fButtonPressed = kFALSE;
      fCompacted = kFALSE;
      fPlane = 0;
      fSpacePressedFrame = 0;
      fPlacePopup = kFALSE;
      fFrameMenuTrash->Delete();
      fMenuObject = 0;
   }

   ~TGuiBldDragManagerPimpl() {
      int i;
      for (i = 0; i <8; i++) {
         delete fGrabRect[i];
      }
      for (i = 0; i <4; i++) {
         delete fAroundFrame[i];
      }

      delete fRepeatTimer;
      delete fGrab;
      fFrameMenuTrash->Delete();
      delete fFrameMenuTrash;

      if (fPlane) {
         fPlane->ChangeOptions(fPlane->GetOptions() & ~kRaisedFrame);
         gClient->NeedRedraw(fPlane, kTRUE);
         fPlane = 0;
      }
   }
};


////////////////////////////////////////////////////////////////////////////////
TGuiBldDragManager::TGuiBldDragManager() : TVirtualDragManager() ,
                    TGFrame(gClient->GetDefaultRoot(), 1, 1)
{
   // Constructor. Create "fantom window".

   SetWindowAttributes_t attr;
   attr.fMask             = kWAOverrideRedirect | kWASaveUnder;
   attr.fOverrideRedirect = kTRUE;
   attr.fSaveUnder        = kTRUE;

   gVirtualX->ChangeWindowAttributes(fId, &attr);

   gGuiBldDragManager = this;
   fPimpl = new TGuiBldDragManagerPimpl(this);

   fSelectionIsOn = kFALSE;
   fFrameMenu     = 0;
   fLassoMenu     = 0;
   fEditor        = 0;
   fBuilder       = 0;
   fLassoDrawn    = kFALSE;
   fDropStatus    = kFALSE;
   fStop          = kTRUE;
   fSelected      = 0;
   fListOfDialogs = 0;

   Reset1();
   CreateListOfDialogs();

   TString tmpfile = gSystem->TempDirectory();
   char *s = gSystem->ConcatFileName(tmpfile.Data(),
               TString::Format("RootGuiBldClipboard%d.C", gSystem->GetPid()));
   fPasteFileName = s;
   delete [] s;

   s = gSystem->ConcatFileName(tmpfile.Data(),
               TString::Format("RootGuiBldTmpFile%d.C", gSystem->GetPid()));
   fTmpBuildFile = s;
   delete [] s;

   fName = "Gui Builder Drag Manager";
   SetWindowName(fName.Data());

   // let's try to solve the problems by myself
   SetErrorHandler(GuiBldErrorHandler);

   fClient->UnregisterWindow(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TGuiBldDragManager::~TGuiBldDragManager()
{
   SetEditable(kFALSE);

   delete fPimpl;

   delete fBuilder;
   fBuilder = 0;

//   delete fEditor;
//   fEditor = 0;

   delete fFrameMenu;
   fFrameMenu =0;

   delete fLassoMenu;
   fLassoMenu = 0;

   if (!gSystem->AccessPathName(fPasteFileName.Data())) {
      gSystem->Unlink(fPasteFileName.Data());
   }

   delete fListOfDialogs;

   gGuiBldDragManager = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Reset some parameters

void TGuiBldDragManager::Reset1()
{
   TVirtualDragManager::Init();
   fTargetId = 0;
   fPimpl->fPlacePopup = kFALSE;
   SetCursorType(kPointer);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a list of dialog methods

void TGuiBldDragManager::CreateListOfDialogs()
{
   fListOfDialogs = new TList();

   TList *methodList = IsA()->GetListOfMethods();
   TIter next(methodList);
   TString str;
   TMethod *method;

   while ((method = (TMethod*) next())) {
      str = method->GetCommentString();
      if (str.Contains("*DIALOG")) {
         fListOfDialogs->Add(method);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw grid on editable frame and restore background on previuosly editted one

void TGuiBldDragManager::Snap2Grid()
{
   if (fStop) {
      return;
   }

   delete fPimpl->fGrid;

   fPimpl->fGrid = new TGuiBldDragManagerGrid();
   fPimpl->fGrid->Draw();
}

////////////////////////////////////////////////////////////////////////////////
/// Return the grid step

UInt_t TGuiBldDragManager::GetGridStep()
{
   return fPimpl->fGrid ? fPimpl->fGrid->fgStep : 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the grid step

void TGuiBldDragManager::SetGridStep(UInt_t step)
{
   fPimpl->fGrid->SetStep(step);
}

////////////////////////////////////////////////////////////////////////////////
/// Return kTRUE if event is rejected for processing by drag manager

Bool_t TGuiBldDragManager::IgnoreEvent(Event_t *event)
{
   if (fStop || !fClient || !fClient->IsEditable()) return kTRUE;
   if (event->fType == kClientMessage) return kFALSE;
   if (event->fType == kDestroyNotify) return kFALSE;

   TGWindow *w = fClient->GetWindowById(event->fWindow);

   if (w) {
      if (IsEditDisabled(w)) {
         w = GetEditableParent((TGFrame*)w);
         return !w;
      }
   } else {
      return kTRUE;
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return a pointer to the parent window (which is being editted)

TGFrame* TGuiBldDragManager::InEditable(Window_t id)
{
   if (fStop || !id) {
      return 0;
   }

   Window_t preparent = id;
   Window_t parent = (Window_t)gVirtualX->GetParent(id);

   while (!parent || (parent != fClient->GetDefaultRoot()->GetId())) {
      if (parent == fClient->GetRoot()->GetId()) {
         TGWindow *w = fClient->GetWindowById(preparent);
         return (w ? (TGFrame*)w : 0);
      }
      preparent = parent;
      parent = gVirtualX->GetParent(parent);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Find the first composite parent of window

TGCompositeFrame *TGuiBldDragManager::FindCompositeFrame(Window_t id)
{
   if (fStop || !id) {
      return 0;
   }

   Window_t parent = id;

   while (!parent || (parent != fClient->GetDefaultRoot()->GetId())) {
      TGWindow *w = fClient->GetWindowById(parent);
      if (w) {
         if (w->InheritsFrom(TGCompositeFrame::Class())) {
            return (TGCompositeFrame*)w;
         }
      }
      parent = gVirtualX->GetParent(parent);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Set cursor for selcted/grabbed frame.

void TGuiBldDragManager::SetCursorType(Int_t cur)
{
   if (fStop) {
      return;
   }

   static UInt_t gid = 0;
   static UInt_t rid = 0;

   if (fPimpl->fGrab && (gid != fPimpl->fGrab->GetId())) {
      gVirtualX->SetCursor(fPimpl->fGrab->GetId(),
                           gVirtualX->CreateCursor((ECursor)cur));
      gid = fPimpl->fGrab->GetId();
   }
   if (fClient->IsEditable() && (rid != fClient->GetRoot()->GetId())) {
      gVirtualX->SetCursor(fClient->GetRoot()->GetId(),
                           gVirtualX->CreateCursor((ECursor)cur));
      rid = fClient->GetRoot()->GetId();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Check resize type event.

Bool_t TGuiBldDragManager::CheckDragResize(Event_t *event)
{
   if (fStop) {
      return kFALSE;
   }

   Bool_t ret = kFALSE;
   fPimpl->fResizeType = kPointer;

   for (int i = 0; i < 8; i++) {
      if (fPimpl->fGrabRect[i]->GetId() == event->fWindow) {
         fPimpl->fResizeType = fPimpl->fGrabRect[i]->GetType();
         ret = kTRUE;
      }
   }

   if ((event->fType == kButtonPress) && (fPimpl->fResizeType != kPointer)) {
      fDragType = kDragResize;
      ret = kTRUE;
   }

   SetCursorType(ret ? fPimpl->fResizeType : kPointer);
   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// Redraw the editted window

void TGuiBldDragManager::DoRedraw()
{
   if (fStop || !fClient || !fClient->IsEditable()) {
      return;
   }

   TGWindow *root = (TGWindow*)fClient->GetRoot();

   fClient->NeedRedraw(root, kTRUE);

   if (fBuilder) {
      fClient->NeedRedraw(fBuilder, kTRUE);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Switch editable

void TGuiBldDragManager::SwitchEditable(TGFrame *frame)
{
   if (fStop || !frame) {
      return;
   }

   TGCompositeFrame *comp = 0;

   if (frame->InheritsFrom(TGCompositeFrame::Class()) && CanChangeLayout(frame)) {
      comp = (TGCompositeFrame *)frame;
   } else if (frame->GetParent()->InheritsFrom(TGCompositeFrame::Class())) {
      comp = (TGCompositeFrame *)frame->GetParent();
   }

   if (!comp) {
      return;
   }

   TString str = comp->ClassName();
   str += "::";
   str += comp->GetName();

   if (IsEditDisabled(comp)) {
      if (fBuilder) {
         str += " cannot be editted.";
         fBuilder->UpdateStatusBar(str.Data());
      }
      return;
   }

   if (frame != comp) {
      SelectFrame(frame);
   }

   if (comp->IsEditable()) {
      return;
   }

   RaiseMdiFrame(comp);
   comp->SetEditable(kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Grab/Select frame

void TGuiBldDragManager::SelectFrame(TGFrame *frame, Bool_t add)
{
   if (fStop || !frame || (frame->GetParent() == fClient->GetDefaultRoot()) ||
       !fClient->IsEditable()) {
      return;
   }

   TString str = frame->ClassName();
   str += "::";
   str += frame->GetName();

   if (IsGrabDisabled(frame)) {
      if (fBuilder) {
         str += "can not be selected";
         fBuilder->UpdateStatusBar(str.Data());
      }
      return;
   }

   // do not grab mdi frames (quick hack)
   if (fBuilder && frame->InheritsFrom(TGMdiFrame::Class())) {
      return;
   }


   static Int_t x, x0, y, y0, xx, yy;
   Window_t c;

   RaiseMdiFrame(FindMdiFrame(frame));
   frame->MapRaised();

   if (!add) {

      fDragType = (fDragType != kDragCopy ) ?  kDragMove : fDragType;

      gVirtualX->TranslateCoordinates(frame->GetId(),
                                      fClient->GetDefaultRoot()->GetId(),
                                      0, 0, x0, y0, c);

      x = x0 + frame->GetWidth();
      y = y0 + frame->GetHeight();

      if (fBuilder) {
         str += " selected";
         str += (IsEditDisabled(frame) || IsFixedLayout(frame)  ? ". This frame cannot be editted." :
                 " ");
         str += " Press SpaceBar to unselect the frame.";
         if (IsFixedSize(frame)) str += " This frame cannot be resized.";

         fBuilder->UpdateStatusBar(str.Data());
      }

   } else { //shift mask is on

      gVirtualX->TranslateCoordinates(frame->GetId(),
                                      fClient->GetDefaultRoot()->GetId(),
                                      0, 0, xx, yy, c);

      fDragType = kDragLasso;
      fPimpl->fX0 = x0 = TMath::Min(x0, xx);
      fPimpl->fX = x = TMath::Max(x, xx + (Int_t)frame->GetWidth());
      fPimpl->fY0 = y0 = TMath::Min(y0, yy);
      fPimpl->fY = y = TMath::Max(y, yy + (Int_t)frame->GetHeight());

      DrawLasso();
   }

   fFrameUnder = fPimpl->fGrab = frame;
   fPimpl->fGrab->RequestFocus();

   // quick hack. the special case for TGCanvases
   if (frame->InheritsFrom(TGCanvas::Class())) {
      fSelected = ((TGCanvas*)frame)->GetContainer();

      if (!IsEditDisabled(fSelected)) {
         fSelected->SetEditable(kTRUE);
         if (fBuilder && fBuilder->GetAction()) {
            PlaceFrame((TGFrame*)fBuilder->ExecuteAction(), 0);
         }
      }
   } else {
      fSelected = fPimpl->fGrab;
   }
   ChangeSelected(fPimpl->fGrab);

   SetCursorType(kMove);

   SetLassoDrawn(kFALSE);
   DrawGrabRectangles(fPimpl->fGrab);
}

////////////////////////////////////////////////////////////////////////////////
/// Inform outside wold that selected frame was changed

void TGuiBldDragManager::ChangeSelected(TGFrame *fr)
{
   if (fStop) {
      return;
   }

   TGFrame *sel = fr;

   if (fBuilder && (sel == fBuilder->GetMdiMain()->GetCurrent())) {
      sel = 0;
   }

   if (!fr) {
      UngrabFrame();
   }

   if (fEditor) {
      fEditor->ChangeSelected(sel);
   }

   if (fBuilder) {
      fBuilder->ChangeSelected(sel);
      //fBuilder->Update();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// grab frame (see SelectFrame)

void TGuiBldDragManager::GrabFrame(TGFrame *frame)
{
   if (fStop || !frame || !fClient->IsEditable()) {
      return;
   }

   fPimpl->fGrabParent = frame->GetParent();
   fPimpl->fGrabX = frame->GetX();
   fPimpl->fGrabY = frame->GetY();

   Window_t c;

   gVirtualX->TranslateCoordinates(frame->GetId(),
                                   fClient->GetDefaultRoot()->GetId(),
                                   0, 0, fPimpl->fX0, fPimpl->fY0, c);

   fPimpl->fX = fPimpl->fX0;
   fPimpl->fY = fPimpl->fY0;

   if (frame->GetFrameElement() && frame->GetFrameElement()->fLayout) {
      fPimpl->fGrabLayout = frame->GetFrameElement()->fLayout;
   }

   if (fPimpl->fGrabParent && frame->GetFrameElement() &&
       fPimpl->fGrabParent->InheritsFrom(TGCompositeFrame::Class())) {
      TList *li = ((TGCompositeFrame*)fPimpl->fGrabParent)->GetList();
      fPimpl->fGrabListPosition = (TGFrameElement*)li->Before(frame->GetFrameElement());
      ((TGCompositeFrame*)fPimpl->fGrabParent)->RemoveFrame(frame);
   }

   SetWindowAttributes_t attr;
   attr.fMask             = kWAOverrideRedirect | kWASaveUnder;
   attr.fOverrideRedirect = kTRUE;
   attr.fSaveUnder        = kTRUE;

   gVirtualX->ChangeWindowAttributes(frame->GetId(), &attr);

   frame->UnmapWindow();
   frame->ReparentWindow(fClient->GetDefaultRoot(), fPimpl->fX0, fPimpl->fY0);
   gVirtualX->Update(1);
   frame->Move(fPimpl->fX0, fPimpl->fY0);
   frame->MapRaised();

   if (fBuilder) {
      //fBuilder->Update();
      TString str = frame->ClassName();
      str += "::";
      str += frame->GetName();
      str += " is grabbed";

      fBuilder->UpdateStatusBar(str.Data());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Ungrab/Unselect selected/grabbed frame.

void TGuiBldDragManager::UngrabFrame()
{
   if (fStop || !fPimpl->fGrab) {
      return;
   }

   SetCursorType(kPointer);
   HideGrabRectangles();

   DoRedraw();

   if (fBuilder) {
      //fBuilder->Update();
      TString str = fPimpl->fGrab->ClassName();
      str += "::";
      str += fPimpl->fGrab->GetName();
      str += " ungrabbed";
      fBuilder->UpdateStatusBar(str.Data());
   }
   fSelected = fPimpl->fGrab = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Helper for IsPointVisible

static Bool_t IsParentOfGrab(Window_t id, const TGWindow *grab)
{
   const TGWindow *parent = grab;

   while (parent && (parent != gClient->GetDefaultRoot())) {
      if (parent->GetId() == id) {
         return kTRUE;
      }
      parent = parent->GetParent();
   }

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Helper function for IsSelectedWindow method

Bool_t TGuiBldDragManager::IsPointVisible(Int_t xi, Int_t yi)
{
   Window_t w = gVirtualX->GetDefaultRootWindow();
   Window_t src, dst, child;
   Int_t x = xi;
   Int_t y = yi;
   Bool_t ret = kFALSE;

   gVirtualX->TranslateCoordinates(fPimpl->fGrab->GetId(), w, x, y, x, y, child);

   dst = src = child = w;

   while (child) {
      src = dst;
      dst = child;
      gVirtualX->TranslateCoordinates(src, dst, x, y, x, y, child);

      if (IsParentOfGrab(child, fPimpl->fGrab)) {
         return kTRUE;
      }
   }

   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// Return kTRUE if grabbed/selected frame is not overlapped by other windows.

Bool_t TGuiBldDragManager::IsSelectedVisible()
{
   if (fStop || !fPimpl->fGrab || !fClient->IsEditable()) {
      return kFALSE;
   }

   if (fBuilder) {
      TGMdiFrame *mdi = fBuilder->FindEditableMdiFrame(fPimpl->fGrab);
      if (mdi && (mdi != fBuilder->GetMdiMain()->GetCurrent())) {
         return kFALSE;
      }
   }

   // popup menu was placed
   if (fPimpl->fPlacePopup) {
      return kTRUE;
   }

   static Long64_t was = gSystem->Now();
   static Bool_t visible = kFALSE;

   Long64_t now = gSystem->Now();

   if (now-was < 100) {
      return visible;
   }
   was = now;

   visible = kFALSE;

   if (!IsPointVisible(2, 2)) {
      return visible;
   }

   if (!IsPointVisible(2, fPimpl->fGrab->GetHeight()-2)) {
      return visible;
   }

   if (!IsPointVisible(fPimpl->fGrab->GetWidth()-2, 2)) {
      return visible;
   }

   if (!IsPointVisible(fPimpl->fGrab->GetWidth()-2,
                       fPimpl->fGrab->GetHeight()-2)) {
      return visible;
   }

   visible = kTRUE;

   return visible;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw small grab rectangles around grabbed/selected/frame

void TGuiBldDragManager::DrawGrabRectangles(TGWindow *win)
{
   if (fStop) {
      return;
   }

   TGFrame *frame = win ? (TGFrame *)win : fPimpl->fGrab;

   if (!frame || !fClient->IsEditable() || fPimpl->fPlacePopup) {
      return;
   }

   Window_t w = gVirtualX->GetDefaultRootWindow();
   Window_t c; Int_t x, y;

   gVirtualX->TranslateCoordinates(frame->GetId(), w,  0, 0, x, y, c);

   if (frame->InheritsFrom(TGCompositeFrame::Class()) &&
       CanChangeLayout(frame) && !frame->IsLayoutBroken()) {
      fPimpl->fAroundFrame[0]->MoveResize(x-3, y-3, frame->GetWidth()+6, 2);
      fPimpl->fAroundFrame[0]->MapRaised();
      fPimpl->fAroundFrame[1]->MoveResize(x+frame->GetWidth()+3, y-3, 2, frame->GetHeight()+6);
      fPimpl->fAroundFrame[1]->MapRaised();
      fPimpl->fAroundFrame[2]->MoveResize(x-3, y+frame->GetHeight()+2, frame->GetWidth()+6, 2);
      fPimpl->fAroundFrame[2]->MapRaised();
      fPimpl->fAroundFrame[3]->MoveResize(x-3, y-3, 2, frame->GetHeight()+6);
      fPimpl->fAroundFrame[3]->MapRaised();
   } else {
      for (int i = 0; i < 4; i++) fPimpl->fAroundFrame[i]->UnmapWindow();
   }

   // draw rectangles
   DrawGrabRect(0, x - 6, y - 6);
   DrawGrabRect(1, x + frame->GetWidth()/2 - 3, y - 6);
   DrawGrabRect(2, x + frame->GetWidth(), y - 6);
   DrawGrabRect(3, x - 6, y + frame->GetHeight());
   DrawGrabRect(4, x - 6, y + frame->GetHeight()/2 - 3);
   DrawGrabRect(5, x + frame->GetWidth(), y + frame->GetHeight()/2 - 3);
   DrawGrabRect(6, x + frame->GetWidth()/2 - 3, y + frame->GetHeight());
   DrawGrabRect(7, x + frame->GetWidth(), y + frame->GetHeight());

   fPimpl->fGrabRectHidden = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Helper method to draw grab rectangle at position x,y

void TGuiBldDragManager::DrawGrabRect(Int_t i, Int_t x, Int_t y)
{
   if (fStop) {
      return;
   }

   fPimpl->fGrabRect[i]->Move(x, y);
   fPimpl->fGrabRect[i]->MapRaised();
}

////////////////////////////////////////////////////////////////////////////////
/// Raise composite frame when mouse is moving over it.
/// That allows to highlight position of "plain" composite frames.

void TGuiBldDragManager::HighlightCompositeFrame(Window_t win)
{
   static Window_t gw = 0;

   if (fStop || !win || (win == gw)) {
      return;
   }

   TGWindow *w = fClient->GetWindowById(win);

   if (!w || (w == fPimpl->fPlane) || w->GetEditDisabled() || w->IsEditable() ||
       !w->InheritsFrom(TGCompositeFrame::Class())) {
      return;
   }

   TGFrame *frame = (TGFrame*)w;
   UInt_t opt = frame->GetOptions();

   if ((opt & kRaisedFrame) || (opt & kSunkenFrame)) {
      return;
   }

   gw = win;
   if (fPimpl->fPlane) {
      fPimpl->fPlane->ChangeOptions(fPimpl->fPlane->GetOptions() & ~kRaisedFrame);
      fClient->NeedRedraw(fPimpl->fPlane, kTRUE);
   }
   fPimpl->fPlane = frame;
   fPimpl->fPlane->ChangeOptions(opt | kRaisedFrame);
   fClient->NeedRedraw(fPimpl->fPlane, kTRUE);

   if (fBuilder) {
      TString str = frame->ClassName();
      str += "::";
      str += frame->GetName();
      fBuilder->UpdateStatusBar(str.Data());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// The main event loop is  originated here
/// It repeadeatly queries pointer state and position on the screen.
/// From this info an Event_t structure is built.

Bool_t TGuiBldDragManager::HandleTimer(TTimer *t)
{
   return HandleTimerEvent(0, t);
}

////////////////////////////////////////////////////////////////////////////////
/// Handle timer events or events coming from the recorder.

Bool_t TGuiBldDragManager::HandleTimerEvent(Event_t *e, TTimer *t)
{
   static Int_t gy = 0;
   static Int_t gx = 0;
   static UInt_t gstate = 0;
   static Window_t gw = 0;

   Bool_t ret = kTRUE;

   // if nothing is editted stop timer and reset everything
   if (!fClient || !fClient->IsEditable()) {
      SetEditable(kFALSE);
      return kFALSE;
   }
   if (!IsSelectedVisible()) {
      HideGrabRectangles();
   }
   if (e) {
      if (fPimpl->fRepeatTimer) {
         // we are replaying events from the recorder...
         fPimpl->fRepeatTimer->Reset();
         fPimpl->fRepeatTimer->Remove();
      }
      if (e->fType == kButtonPress)
         return HandleButtonPress(e);
      else if (e->fType == kButtonRelease)
         return HandleButtonRelease(e);
      else if (e->fState & kButton1Mask)
         return HandleMotion(e);
      return kTRUE;
   }
   Window_t dum;
   Event_t ev;
   ev.fCode = kButton1;
   ev.fType = kMotionNotify;
   ev.fState = 0;

   gVirtualX->QueryPointer(gVirtualX->GetDefaultRootWindow(), dum, dum,
                           ev.fXRoot, ev.fYRoot, ev.fX, ev.fY, ev.fState);

   ev.fWindow = GetWindowFromPoint(ev.fXRoot, ev.fYRoot);

   if (ev.fWindow && (gw == ev.fWindow) && (gstate == ev.fState) &&
       (ev.fYRoot == gy) && (ev.fXRoot == gx)) {
      return kFALSE;
   }

   gw = ev.fWindow;
   gstate = ev.fState;
   ev.fState &= ~16;    // ignore "num lock" pressed
   ev.fState &= ~2;     // ignore "caps lock" pressed

   if (!fDragging && !fMoveWaiting && !fPimpl->fButtonPressed &&
       ((ev.fState == kButton1Mask) || (ev.fState == kButton3Mask) ||
        (ev.fState == (kButton1Mask | kKeyShiftMask)) ||
        (ev.fState == (kButton1Mask | kKeyControlMask)))) {

      if (ev.fState & kButton1Mask) ev.fCode = kButton1;
      if (ev.fState & kButton3Mask) ev.fCode = kButton3;

      ev.fType = kButtonPress;
      t->SetTime(40);

      if (fPimpl->fPlane && fClient->GetWindowById(fPimpl->fPlane->GetId())) {
         fPimpl->fPlane->ChangeOptions(fPimpl->fPlane->GetOptions() & ~kRaisedFrame);
         fClient->NeedRedraw(fPimpl->fPlane, kTRUE);
      } else {
         fPimpl->fPlane = 0;
      }

      ret = HandleButtonPress(&ev);
      TimerEvent(&ev);
      return ret;
   }

   if ((fDragging || fMoveWaiting) && (!ev.fState || (ev.fState == kKeyShiftMask)) &&
       fPimpl->fButtonPressed) {

      ev.fType = kButtonRelease;
      t->SetTime(100);

      ret = HandleButtonRelease(&ev);
      TimerEvent(&ev);
      return ret;
   }

   fPimpl->fButtonPressed = (ev.fState & kButton1Mask) ||
                            (ev.fState & kButton2Mask) ||
                            (ev.fState & kButton3Mask);

   if ((ev.fYRoot == gy) && (ev.fXRoot == gx)) return kFALSE;

   gy = ev.fYRoot;
   gx = ev.fXRoot;

   if (!fMoveWaiting && !fDragging && !ev.fState) {
      if (!CheckDragResize(&ev) && fClient->GetWindowById(ev.fWindow)) {
         HighlightCompositeFrame(ev.fWindow);
      }
   } else if (ev.fState & kButton1Mask) {
      HandleMotion(&ev);
      TimerEvent(&ev);
   }
   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// Recognize what was done when mouse button pressed

Bool_t TGuiBldDragManager::RecognizeGesture(Event_t *event, TGFrame *frame)
{
   if (fStop) {
      return kFALSE;
   }

   if (((event->fCode != kButton1) && (event->fCode != kButton3)) ||
       !frame || !fClient->IsEditable()) {
      return kFALSE;
   }

   TGFrame *context_fr = 0;
   Bool_t mdi = kFALSE;

   // hack for editable mdi frames
   if (frame->IsEditable() && frame->InheritsFrom(TGMdiFrame::Class())) {
      context_fr = frame;
      mdi = kTRUE;
   }

   // handle context menu
   if (event->fCode == kButton3) {
      if (!fPimpl->fSpacePressedFrame) {
         if (!mdi) {
            SelectFrame(frame);
            context_fr = fSelected;
         }
      } else {
         context_fr = fPimpl->fSpacePressedFrame;
      }

      HandleButon3Pressed(event, context_fr);
      return kTRUE;
   }

   fDragType = kDragNone;

   if (!fSelectionIsOn) {
      fPimpl->fX0 = event->fXRoot;
      fPimpl->fY0 = event->fYRoot;
   }

   //HideGrabRectangles();
   fPimpl->fClickFrame = frame;

   if (fBuilder && fBuilder->IsExecutable() &&
       frame->InheritsFrom(TGCompositeFrame::Class())) {
      UngrabFrame();
      frame->SetEditable(kTRUE);
      fSource = 0;
      fDragType = kDragLasso;
      goto out;
   }

   if (event->fState & kKeyShiftMask) {
      // drag grabbed frame with shift key pressed should create a copy of grabbed frame
      // move a copy of editable and selected frame
      if (frame == fPimpl->fGrab) {
         fSource = frame;
         fDragType = kDragCopy;
         gVirtualX->SetCursor(frame->GetId(), gVirtualX->CreateCursor(kMove));
         goto out;
      }

      // otherwise do lasso selection
      if (!fSelectionIsOn) {
         fSelectionIsOn = kTRUE;
      } else {
         fPimpl->fX = event->fXRoot;
         fPimpl->fY = event->fYRoot;
         fDragType = kDragLasso;
         DrawLasso();
         return kTRUE;
      }
   }

   CheckDragResize(event);

   if (frame->IsEditable()) {
      fSource = 0;

      if (fDragType != kDragResize) {

         // move editable and selected frame
         if (frame == fPimpl->fGrab) {
            fSource = frame;
            fDragType = kDragMove;
            gVirtualX->SetCursor(frame->GetId(), gVirtualX->CreateCursor(kMove));
            goto out;
         }

         fDragType = kDragLasso;
      }
   } else if ((fDragType != kDragResize) && !fPimpl->fSpacePressedFrame) {

      // special case of TGCanvas
      if (!fPimpl->fGrab && frame->InheritsFrom(TGCanvas::Class()))  {
         TGFrame *cont = ((TGCanvas*)frame)->GetContainer();

         if (!cont->IsEditable()) {
            cont->SetEditable(kTRUE);
            fDragType = kDragLasso;
            goto out;
         }
      }

      fSource = frame;
      SelectFrame(frame, event->fState & kKeyShiftMask);
   }

   if ((fDragType == kDragNone) && !fPimpl->fSpacePressedFrame) {
      SwitchEditable(frame);
      fSource = 0;

      // try again
      CheckDragResize(event);

      if (fDragType == kDragNone) {
         return kFALSE;
      }
   }

out:
   Window_t c;

   gVirtualX->TranslateCoordinates(fClient->GetDefaultRoot()->GetId(),
                                   frame->GetId(),
                                   event->fXRoot, event->fYRoot,
                                   fPimpl->fXf, fPimpl->fYf, c);
   fPimpl->fX = event->fXRoot;
   fPimpl->fY = event->fYRoot;

   fMoveWaiting = kTRUE;
   DoRedraw();

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle 3d mouse pressed (popup context menu)

void TGuiBldDragManager::HandleButon3Pressed(Event_t *event, TGFrame *frame)
{
   if (fStop || !frame) {
      return;
   }

   if (fClient->GetWaitForEvent() == kUnmapNotify) {
      return;
   }

   if (frame == fSelected) {
      Menu4Frame(frame, event->fXRoot, event->fYRoot);
   } else if (frame->IsEditable())  {
      if (fLassoDrawn) {
         Menu4Lasso(event->fXRoot, event->fYRoot);
      } else {
         Menu4Frame(frame, event->fXRoot, event->fYRoot);
      }
   } else {
      TGFrame *base = InEditable(frame->GetId());
      if (base) {
         //SelectFrame(base);
         Menu4Frame(base, event->fXRoot, event->fYRoot);
      } else {
         Menu4Frame(frame, event->fXRoot, event->fYRoot);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Handle button event occured in some ROOT frame

Bool_t TGuiBldDragManager::HandleButton(Event_t *event)
{
   if (fStop) {
      return kFALSE;
   }

   if (event->fCode != kButton3) {
      CloseMenus();
   }

   if (event->fType == kButtonPress) {
      return HandleButtonPress(event);
   } else {
      return HandleButtonRelease(event);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Resize events

Bool_t TGuiBldDragManager::HandleConfigureNotify(Event_t *event)
{
   if (fStop) {
      return kFALSE;
   }

   TGWindow *w = fClient->GetWindowById(event->fWindow);

   if (!w) {
      return kFALSE;
   }

   fPimpl->fCompacted = kFALSE;
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle repaint event

Bool_t TGuiBldDragManager::HandleExpose(Event_t *event)
{
   if (fStop) {
      return kFALSE;
   }

   static Long64_t was = gSystem->Now();
   static Window_t win = 0;
   Long64_t now = gSystem->Now();

   if (event->fCount || (win == event->fWindow) || (now-was < 50) || fDragging) {
      if (fDragging) {
         HideGrabRectangles();
      }
      return kFALSE;
   }

   if (gMenuDialog) {
      HideGrabRectangles();
      gMenuDialog->RaiseWindow();
      return kFALSE;
   }

   if (fLassoDrawn) {
      DrawLasso();
   } else {
      if (IsSelectedVisible()) {
         DrawGrabRectangles();
      }
   }

   win = event->fWindow;
   was = now;

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle all events.

Bool_t TGuiBldDragManager::HandleEvent(Event_t *event)
{
   if (fStop) {
      return kFALSE;
   }

   if (IgnoreEvent(event)) {
      return kFALSE;
   }

   switch (event->fType) {

      case kExpose:
         return HandleExpose(event);

      case kConfigureNotify:
         while (gVirtualX->CheckEvent(fId, kConfigureNotify, *event))
            ;
         return HandleConfigureNotify(event);

      case kGKeyPress:
      case kKeyRelease:
         return HandleKey(event);

      case kFocusIn:
      case kFocusOut:
         //HandleFocusChange(event);
         break;

      case kButtonPress:
         {
            Int_t dbl_clk = kFALSE;

            static Window_t gDbw = 0;
            static Long_t gLastClick = 0;
            static UInt_t gLastButton = 0;
            static Int_t gDbx = 0;
            static Int_t gDby = 0;

            if ((event->fTime - gLastClick < 350) &&
                (event->fCode == gLastButton) &&
                (TMath::Abs(event->fXRoot - gDbx) < 6) &&
                (TMath::Abs(event->fYRoot - gDby) < 6) &&
                (event->fWindow == gDbw)) {
               dbl_clk = kTRUE;
            }

            if (dbl_clk) {
               if (event->fState & kKeyControlMask) {
                  HandleAction(kEndEditAct);
                  return kTRUE;
               } else if (!(event->fState & 0xFF)) {
                  TGFrame *w = (TGFrame*)fClient->GetWindowById(event->fWindow);

                  if (w && (w->GetEditDisabled() & kEditDisableBtnEnable)) {
                     return w->HandleDoubleClick(event);
                  }
                  if (SaveFrame(fTmpBuildFile.Data())) {
                     gROOT->Macro(fTmpBuildFile.Data());
                  }
                  // an easy way to start editting
                  if (fBuilder) fBuilder->HandleMenu(kGUIBLD_FILE_START);
                  return kTRUE;
               }
            } else {
               gDbw = event->fWindow;
               gLastClick = event->fTime;
               gLastButton = event->fCode;
               gDbx = event->fXRoot;
               gDby = event->fYRoot;

               Bool_t ret = HandleButtonPress(event);
               return ret;
            }

            return kFALSE;
         }

      case kButtonRelease:
         return HandleButtonRelease(event);

      case kEnterNotify:
      case kLeaveNotify:
         //HandleCrossing(event);
         break;

      case kMotionNotify:
         while (gVirtualX->CheckEvent(fId, kMotionNotify, *event))
            ;
         return HandleMotion(event);

      case kClientMessage:
         return HandleClientMessage(event);

      case kDestroyNotify:
         return HandleDestroyNotify(event);

      case kSelectionNotify:
         //HandleSelection(event);
         break;

      case kSelectionRequest:
         //HandleSelectionRequest(event);
         break;

      case kSelectionClear:
         //HandleSelectionClear(event);
         break;

      case kColormapNotify:
         //HandleColormapChange(event);
         break;

      default:
         //Warning("HandleEvent", "unknown event (%#x) for (%#x)", event->fType, fId);
         break;
   }

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Mouse double click handler (never should happen)

Bool_t TGuiBldDragManager::HandleDoubleClick(Event_t *)
{
   if (fStop) {
      return kFALSE;
   }

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return a parent which can handle button evevents.

TGFrame *TGuiBldDragManager::GetBtnEnableParent(TGFrame *fr)
{
   TGWindow *parent = fr;

   while (parent && (parent != fClient->GetDefaultRoot())) {
      if (parent->GetEditDisabled() & kEditDisableBtnEnable) {
         return (TGFrame*)parent;
      }
      parent = (TGWindow*)parent->GetParent();
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Unmap all popups

void TGuiBldDragManager::UnmapAllPopups()
{
   TList *li = fClient->GetListOfPopups();
   if (!li->GetEntries()) {
      return;
   }

   TGPopupMenu *pup;
   TIter next(li);

   while ((pup = (TGPopupMenu*)next())) {
      pup->UnmapWindow();
      fClient->ResetWaitFor(pup);
   }
   gVirtualX->GrabPointer(0, 0, 0, 0, kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Handle button press event

Bool_t TGuiBldDragManager::HandleButtonPress(Event_t *event)
{
   if (fStop) {
      return kFALSE;
   }

   fPimpl->fButtonPressed = kTRUE;
   fPimpl->fPlacePopup = kFALSE;

   if (fPimpl->fPlane) {
      fPimpl->fPlane->ChangeOptions(fPimpl->fPlane->GetOptions() & ~kRaisedFrame);
      fClient->NeedRedraw(fPimpl->fPlane, kTRUE);
   }

   if (gMenuDialog) { // keep editor on the top
      gMenuDialog->RaiseWindow();
   }

   // keep undocked toolbar on the top
   //(but under win32 key handling will be broken : todo)
   if (gVirtualX->InheritsFrom("TGX11") && fBuilder &&
       fBuilder->GetToolDock()->IsUndocked()) {
      fBuilder->GetToolDock()->GetUndocked()->RaiseWindow();
   }

   // keep color dialog on the top
   if (fgGlobalColorDialog && fgGlobalColorDialog->IsMapped()) {
      fgGlobalColorDialog->RaiseWindow();
      return kFALSE;
   }

   if ( ((event->fCode != kButton1) && (event->fCode != kButton3)) ||
        (event->fType != kButtonPress) || IgnoreEvent(event)) {
      return kFALSE;
   }

   Reset1();
   //HideGrabRectangles();

   Window_t w = GetWindowFromPoint(event->fXRoot, event->fYRoot);
   TGFrame *fr = 0;

   if (w) {
      fr = (TGFrame*)fClient->GetWindowById(w);
      if (!fr) {
         return kFALSE;
      }

      //fr->HandleButton(event);
      if (!IsEventsDisabled(fr)) {
         TGFrame *btnframe = GetBtnEnableParent(fr);
         if (btnframe) {
            event->fUser[0] = fr->GetId();
            btnframe->HandleButton(event);
         }
      }

      if (IsGrabDisabled(fr)) {
         fr = GetEditableParent(fr);
      }

      if (!fr) {
         return kFALSE;
      }
   } else {
      return kFALSE;
   }

   return RecognizeGesture(event, fr);
}

////////////////////////////////////////////////////////////////////////////////
/// Handle button release event

Bool_t TGuiBldDragManager::HandleButtonRelease(Event_t *event)
{
   if (fStop) {
      return kFALSE;
   }

   // unmap all waiting popups
   if (fClient->GetWaitForEvent() == kUnmapNotify) {
      UnmapAllPopups();
   }

   TGWindow *w = fClient->GetWindowById(event->fWindow);

   if (w && !IsEventsDisabled(w)) {
      TGFrame *btnframe = GetBtnEnableParent((TGFrame*)w);
      if (btnframe) {
         event->fUser[0] = w->GetId();
         btnframe->HandleButton(event);
      }
   }

   fPimpl->fButtonPressed = kFALSE;
   gVirtualX->SetCursor(fClient->GetRoot()->GetId(), gVirtualX->CreateCursor(kPointer));
   EndDrag();
   fSelectionIsOn &= (event->fState & kKeyShiftMask);

   if (fLassoDrawn) {
      DrawLasso();
      return kTRUE;
   }

   if (fPimpl->fClickFrame && !fSelectionIsOn) {

      // make editable the clicked frame if no lasso was drawn
      if ((fPimpl->fClickFrame == fPimpl->fGrab) && (fSelected == fPimpl->fGrab) &&
           !fPimpl->fGrab->IsEditable()) {
         SwitchEditable(fPimpl->fClickFrame);
         return kTRUE;

      // select/grab clicked frame if there was no grab frame
      } else if (!fPimpl->fGrab || ((fPimpl->fClickFrame != fPimpl->fGrab) &&
                                    (fPimpl->fClickFrame != fSelected))) {
         SelectFrame(fPimpl->fClickFrame);
         return kTRUE;
      }

   }

   SelectFrame(fPimpl->fGrab);

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle key event

Bool_t TGuiBldDragManager::HandleKey(Event_t *event)
{
   if (fStop) {
      return kFALSE;
   }

   char   tmp[10];
   UInt_t keysym;
   Bool_t ret = kFALSE;
   TGFileInfo fi;
   static TString dir(".");
   static Bool_t overwr = kFALSE;
   TString fname;

   TGWindow *w = fClient->GetWindowById(GetWindowFromPoint(event->fXRoot, event->fYRoot));

   if (!w || !fPimpl) {
      return kFALSE;
   }

   if (w->GetEditDisabled() & kEditDisableKeyEnable) {
      return ((TGFrame*)w)->HandleKey(event);
   }

   if (event->fType != kGKeyPress) {
      return kFALSE;
   }

   if (IsEditDisabled(w)) {
      TGFrame *parent = GetEditableParent((TGFrame*)w);
      if (parent) {
         event->fWindow = parent->GetId();
         parent->HandleKey(event);
      } else {
         return ((TGFrame*)w)->HandleKey(event);
      }
   }

   fPimpl->fSpacePressedFrame = 0;

   if (fPimpl->fPlane) {
      fPimpl->fPlane->ChangeOptions(fPimpl->fPlane->GetOptions() & ~kRaisedFrame);
      fClient->NeedRedraw(fPimpl->fPlane, kTRUE);
   }

   CloseMenus();

   fi.fFileTypes = gSaveMacroTypes;
   fi.SetIniDir(dir);
   fi.fOverwrite = overwr;

   gVirtualX->LookupString(event, tmp, sizeof(tmp), keysym);

   if (event->fState & kKeyControlMask) {

      switch ((EKeySym)keysym & ~0x20) {
         case kKey_Return:
         case kKey_Enter:
            HandleReturn(kTRUE);
            ret = kTRUE;
            break;
         case kKey_X:
            HandleCut();
            ret = kTRUE;
            break;
         case kKey_C:
            HandleCopy();
            ret = kTRUE;
            break;
         case kKey_V:
            if (fPimpl->fClickFrame && !fPimpl->fClickFrame->IsEditable()) {
               fPimpl->fClickFrame->SetEditable(kTRUE);
            }
            HandlePaste();
            ret = kTRUE;
            break;
         case kKey_B:
         {
            if (fPimpl->fGrab ) {
               BreakLayout();
            }
            ret = kTRUE;
            break;
         }
         case kKey_L:
         {
            if (fPimpl->fGrab && (fPimpl->fClickFrame != fClient->GetRoot())) {
               Compact(kFALSE);
            } else {
               Compact(kTRUE);
            }
            ret = kTRUE;
            break;
         }
         case kKey_R:
            HandleReplace();
            ret = kTRUE;
            break;
         case kKey_S:
            Save();
            ret = kTRUE;
            break;
         case kKey_G:
            HandleGrid();
            ret = kTRUE;
            break;
         case kKey_H:
            SwitchLayout();
            ret = kTRUE;
            break;
         case kKey_N:
            if (fBuilder) {
               fBuilder->NewProject();
            } else {
               TGMainFrame *main = new TGMainFrame(fClient->GetDefaultRoot(), 300, 300);
               main->MapRaised();
               main->SetEditable(kTRUE);
            }
            ret = kTRUE;
            break;
         case kKey_O:
            if (fBuilder) {
               fBuilder->NewProject();
            } else {
               TGMainFrame *main = new TGMainFrame(fClient->GetDefaultRoot(), 300, 300);
               main->MapRaised();
               main->SetEditable(kTRUE);
            }
            new TGFileDialog(fClient->GetDefaultRoot(), this, kFDSave, &fi);

            if (!fi.fFilename) return kTRUE;
            dir = fi.fIniDir;
            overwr = fi.fOverwrite;
            fname = gSystem->BaseName(gSystem->UnixPathName(fi.fFilename));

            if (fname.EndsWith(".C")) {
               gROOT->Macro(fname.Data());
            } else {
               Int_t retval;
               new TGMsgBox(fClient->GetDefaultRoot(), this, "Error...",
                            TString::Format("file (%s) must have extension .C", fname.Data()),
                            kMBIconExclamation, kMBRetry | kMBCancel, &retval);
               if (retval == kMBRetry) {
                  HandleKey(event);
               }
            }
            ret = kTRUE;
            break;
         default:
            break;
      }
   } else {
      switch ((EKeySym)keysym) {
         case kKey_Delete:
         case kKey_Backspace:
            HandleDelete(event->fState & kKeyShiftMask);
            ret = kTRUE;
            break;
         case kKey_Return:
         case kKey_Enter:
            //UnmapAllPopups();
            HandleReturn(kFALSE);
            ret = kTRUE;
            break;
         case kKey_Left:
         case kKey_Right:
         case kKey_Up:
         case kKey_Down:
            if (fLassoDrawn) {
               HandleAlignment(keysym, event->fState & kKeyShiftMask);
            } else if (fPimpl->fGrab) {
               HandleLayoutOrder((keysym == kKey_Right) || (keysym == kKey_Down));
            }
            ret = kTRUE;
            break;
         case kKey_Space:
            //UnmapAllPopups();
            if (fPimpl->fGrab) {
               SwitchEditable(fPimpl->fGrab);

               TGFrame *p = (TGFrame*)GetEditableParent(fPimpl->fGrab);

               if (p) {
                  if (p == fBuilder->GetMdiMain()->GetCurrent()) {
                     UngrabFrame();
                  } else {
                     SelectFrame(p);
                     fSource = fPimpl->fSpacePressedFrame = p;
                  }
               }
            }
            ret = kTRUE;
            break;
         default:
            break;
      }
   }
   if (fBuilder) {
      fBuilder->SetAction(0);
      //fBuilder->Update();
   }

   if (fLassoDrawn) {
      DrawLasso();
   }

   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// Reparent frames

void TGuiBldDragManager::ReparentFrames(TGFrame *newfr, TGCompositeFrame *oldfr)
{
   if (fStop || !fClient->IsEditable() || (newfr == fClient->GetDefaultRoot())) {
      return;
   }

   Int_t x0, y0, xx, yy;
   Window_t c;
   static TGLayoutHints *hints = new TGLayoutHints(kLHintsNormal, 2, 2, 2, 2);

   if (!newfr || !newfr->GetId() || !oldfr || !oldfr->GetId()) return;

   gVirtualX->TranslateCoordinates(newfr->GetId(), oldfr->GetId(),
                                   0, 0, x0, y0, c);

   x0 = x0 < 0 ? 0 : x0;
   y0 = y0 < 0 ? 0 : y0;
   Int_t x = x0 + newfr->GetWidth();
   Int_t y = y0 + newfr->GetHeight();

   TGCompositeFrame *comp = 0;

   if (newfr->InheritsFrom(TGCompositeFrame::Class())) {
      comp = (TGCompositeFrame*)newfr;
      comp->SetLayoutBroken();
   }

   TIter next(oldfr->GetList());
   TGFrameElement *el;

   while ((el = (TGFrameElement*)next())) {
      TGFrame *frame = el->fFrame;

      if ((frame->GetX() >= x0) && (frame->GetY() >= y0) &&
          (frame->GetX() + (Int_t)frame->GetWidth() <= x) &&
          (frame->GetY() + (Int_t)frame->GetHeight() <= y)) {

         if (frame == fPimpl->fGrab) {
            UngrabFrame();
         }

         oldfr->RemoveFrame(frame);

         gVirtualX->TranslateCoordinates(oldfr->GetId(), newfr->GetId(),
                                        frame->GetX(), frame->GetY(), xx, yy, c);

         frame->ReparentWindow(newfr, xx, yy);

         if (comp) {
            comp->AddFrame(frame, hints); // el->fLayout);
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return the list of frames inside of some area

TList *TGuiBldDragManager::GetFramesInside(Int_t x0, Int_t y0, Int_t x, Int_t y)
{
   if (fStop) {
      return 0;
   }

   Int_t xx, yy;

   if (!fClient->GetRoot()->InheritsFrom(TGCompositeFrame::Class())) {
      return 0;
   }

   TList *list = new TList();

   xx = x0; yy = y0;
   x0 = TMath::Min(xx, x); x = TMath::Max(xx, x);
   y0 = TMath::Min(yy, y); y = TMath::Max(yy, y);

   TIter next(((TGCompositeFrame*)fClient->GetRoot())->GetList());
   TGFrameElement *el;

   while ((el = (TGFrameElement*)next())) {
      if ((el->fFrame->GetX() >= x0) && (el->fFrame->GetY() >= y0) &&
          (el->fFrame->GetX() + (Int_t)el->fFrame->GetWidth() <= x) &&
          (el->fFrame->GetY() + (Int_t)el->fFrame->GetHeight() <= y)) {
         list->Add(el->fFrame);
      }
   }
   if (list->IsEmpty()) {
      delete list;
      return 0;
   }
   return list;
}

////////////////////////////////////////////////////////////////////////////////
/// Drop canvas container

void TGuiBldDragManager::DropCanvas(TGCanvas *canvas)
{
   if (fStop) {
      return;
   }

   TGCompositeFrame *comp = (TGCompositeFrame*)canvas->GetParent();
   comp->SetEditable(kTRUE);

   TGCompositeFrame *cont = (TGCompositeFrame*)canvas->GetContainer();
   Int_t x = canvas->GetX();
   Int_t y = canvas->GetY();

   cont->SetEditDisabled(cont->GetEditDisabled() & ~kEditDisableGrab);
   cont->ReparentWindow(comp, x, y);
   canvas->SetContainer(0);
   comp->AddFrame(cont);
   DeleteFrame(canvas);

   if (fBuilder) {
      TString str = cont->ClassName();
      str += "::";
      str += cont->GetName();
      str += " dropped.";
      fBuilder->UpdateStatusBar(str.Data());
   }
   SelectFrame(cont);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a new TGCanvas and place container into it

void TGuiBldDragManager::PutToCanvas(TGCompositeFrame *cont)
{
   if (fStop || !cont) {
      return;
   }

   TGCompositeFrame *comp = (TGCompositeFrame*)cont->GetParent();
   comp->SetEditable(kTRUE);

   UInt_t w = cont->GetWidth()/2;
   UInt_t h = cont->GetHeight()/2;

   w = w < 100 ? 100 : w;
   h = h < 100 ? 100 : h;

   TGCanvas *canvas = new TGCanvas(comp, w, h);
   canvas->Move(cont->GetX(), cont->GetY());
   comp->RemoveFrame(cont);
   comp->AddFrame(canvas);
   cont->ReparentWindow(canvas->GetViewPort());
   canvas->SetContainer(cont);
   cont->SetCleanup(kDeepCleanup);
   canvas->MapSubwindows();
   canvas->MapWindow();
   SelectFrame(canvas);

   if (fBuilder) {
      fBuilder->UpdateStatusBar("Grab action performed. Presss Cntrl-Return to Drop grabbed frame.");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Handling of  return/enter key pressing
///
/// If on is kFALSE:
///    If Return or Enter key was pressed - Grab Act
///    If lasso is  drawn - new composite frame is created and
///                         all frames inside lasso adopted as childrens.
///    If lasso is not drawn and selected frame is composite one,
///                       - new TGCanvas widget is created and selcted frmae became
///                         container for this canvas.
///
/// If on is kTRUE:
///    If Return or Enter key was pressed with Control Key - Drop Act,
///    The opposite action to the Grab Act.
///    If selected/grabbed frame is not a TGCanvas widget -
///          all frames inside the grabbed/selected frame are "dropped" into
///          the underlying frame and the grabbed frame is deleted.
///
///    If selected/grabbed frame is a TGCanvas widget -
///          container frame "dropped" to editable frame

void TGuiBldDragManager::HandleReturn(Bool_t on)
{
   if (fStop) {
      return;
   }

   Int_t x0, y0, x, y, xx, yy;
   Window_t c;
   TGCompositeFrame *parent = 0;
   TList *li = 0;

   if (!fClient->GetRoot()->InheritsFrom(TGCompositeFrame::Class()) ||
       !fClient->IsEditable()) {
      return;
   }

   // if grabbed frame is editable - we need to switch edit to parent
   if (fPimpl->fGrab && fPimpl->fGrab->IsEditable()) {
      ((TGFrame*)fPimpl->fGrab->GetParent())->SetEditable(kTRUE);
   }

   if (fPimpl->fGrab && !fLassoDrawn) {
      if (!on) {
         if (fPimpl->fGrab->InheritsFrom(TGCompositeFrame::Class()) &&
            !fPimpl->fGrab->InheritsFrom(TGCanvas::Class()) &&
            !fPimpl->fGrab->InheritsFrom(TGContainer::Class()) &&
            CanChangeLayout(fPimpl->fGrab) &&
            CanChangeLayout((TGWindow*)fPimpl->fGrab->GetParent())) {
            PutToCanvas((TGCompositeFrame*)fPimpl->fGrab);
            return;
         }
      } else {

         if ((fPimpl->fGrab->IsA() == TGCanvas::Class()) &&
             !((TGCanvas*)fPimpl->fGrab)->GetContainer()->InheritsFrom(TGContainer::Class()) &&
             CanChangeLayout((TGWindow*)fPimpl->fGrab->GetParent())) {
            DropCanvas((TGCanvas*)fPimpl->fGrab);
            return;
         }
      }
   }

   TGCompositeFrame *comp = (TGCompositeFrame*)fClient->GetRoot();

   if (fLassoDrawn) {

      gVirtualX->TranslateCoordinates(fClient->GetDefaultRoot()->GetId(),
                                      fClient->GetRoot()->GetId(),
                                      fPimpl->fX, fPimpl->fY, x, y, c);
      gVirtualX->TranslateCoordinates(fClient->GetDefaultRoot()->GetId(),
                                      fClient->GetRoot()->GetId(),
                                      fPimpl->fX0, fPimpl->fY0, x0, y0, c);

      xx = x0; yy = y0;
      x0 = TMath::Min(xx, x); x = TMath::Max(xx, x);
      y0 = TMath::Min(yy, y); y = TMath::Max(yy, y);

      li = GetFramesInside(x0, y0, x, y);

      if (!on && li) {
         parent = new TGCompositeFrame(comp, x - x0, y - y0);
         parent->MoveResize(x0, y0, x - x0, y - y0);
         ReparentFrames(parent, comp);

         comp->AddFrame(parent);
         parent->MapWindow();
         SetLassoDrawn(kFALSE);
         SelectFrame(parent);

         if (fBuilder) {
            TString str = "Grab action performed.";
            str += " Press Cntrl-Return to Drop grabbed frames.";
            str += " Presss Return for TCanvas Grab";
            fBuilder->UpdateStatusBar(str.Data());
         }
      }
   } else if (on && fPimpl->fGrab) {

      // check if it is forbidden
      if (!CanChangeLayout(fPimpl->fGrab) ||
          !CanChangeLayout((TGWindow*)fPimpl->fGrab->GetParent())) {
         if (fBuilder) {
            fBuilder->UpdateStatusBar("Drop action disabled");
         }
         return;
      }

      if (fPimpl->fGrab->InheritsFrom(TGCompositeFrame::Class())) {
         parent = (TGCompositeFrame*)fPimpl->fGrab;
      } else {
         //parent = (TGCompositeFrame*)fPimpl->fGrab->GetParent();
      }
      if (parent) {
         ReparentFrames(comp, parent);
         DeleteFrame(fPimpl->fGrab);
         UngrabFrame();
         ChangeSelected(0);   //update editors

         if (fBuilder) {
            fBuilder->UpdateStatusBar("Drop action performed");
         }
      }
   }
   delete li;
}

////////////////////////////////////////////////////////////////////////////////
/// Align frames located inside lasso area.

void TGuiBldDragManager::HandleAlignment(Int_t to, Bool_t lineup)
{
   if (fStop) {
      return;
   }

   Int_t x0, y0, x, y, xx, yy;
   Window_t c;
   TGCompositeFrame *comp = 0;

   if (!fClient->GetRoot()->InheritsFrom(TGCompositeFrame::Class()) ||
       !fClient->IsEditable()) {
      return;
   }

   if (fLassoDrawn) {
      gVirtualX->TranslateCoordinates(fClient->GetDefaultRoot()->GetId(),
                                      fClient->GetRoot()->GetId(),
                                      fPimpl->fX, fPimpl->fY, x, y, c);
      gVirtualX->TranslateCoordinates(fClient->GetDefaultRoot()->GetId(),
                                      fClient->GetRoot()->GetId(),
                                      fPimpl->fX0, fPimpl->fY0, x0, y0, c);

      xx = x0; yy = y0;
      x0 = TMath::Min(xx, x); x = TMath::Max(xx, x);
      y0 = TMath::Min(yy, y); y = TMath::Max(yy, y);

      comp = (TGCompositeFrame*)fClient->GetRoot();

      ToGrid(x, y);
      ToGrid(x0, y0);

      TIter next(comp->GetList());
      TGFrameElement *el;
      TGFrame *prev = 0;

      while ((el = (TGFrameElement*)next())) {
         TGFrame *fr = el->fFrame;

         if ((fr->GetX() >= x0) && (fr->GetY() >= y0) &&
             (fr->GetX() + (Int_t)fr->GetWidth() <= x) &&
             (fr->GetY() + (Int_t)fr->GetHeight() <= y)) {

            switch ((EKeySym)to) {
               case kKey_Left:
                  fr->Move(x0, fr->GetY());
                  if (lineup) {
                     // coverity[dead_error_line]
                     if (prev) fr->Move(fr->GetX(), prev->GetY() + prev->GetHeight());
                     else fr->Move(x0, y0);
                  }
                  break;
               case kKey_Right:
                  fr->Move(x - fr->GetWidth(), fr->GetY());
                  if (lineup) {
                     if (prev) fr->Move(fr->GetX(), prev->GetY() + prev->GetHeight());
                     else fr->Move(x - fr->GetWidth(), y0);
                  }
                  break;
               case kKey_Up:
                  fr->Move(fr->GetX(), y0);
                  if (lineup) {
                     if (prev) fr->Move(prev->GetX() + prev->GetWidth(), fr->GetY());
                     else fr->Move(x0, y0);
                  }
                  break;
               case kKey_Down:
                  fr->Move(fr->GetX(), y - fr->GetHeight());
                  if (lineup) {
                     if (prev) fr->Move(prev->GetX() + prev->GetWidth(), fr->GetY());
                     else fr->Move(x0, y - fr->GetHeight());
                  }
                  break;
               default:
                  break;
            }
            prev = fr;
         }
      }
   }
   if (fLassoDrawn) {
      DrawLasso();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Handle delete or crop action
///
/// crop is kFALSE - delete action
///   - if lasso is drawn -> all frames inside lasso area are deleted
///   - if frame is grabbed/selected -> the frame is deleted
/// crop is kTRUE - crop action
///   - if lasso is drawn -> all frames outside of lasso area are deleted
///   - if frame is grabbed/selected -> all frames except the grabbed frame are deleted
///     In both cases the main frame is shrinked to the size of crop area.

void TGuiBldDragManager::HandleDelete(Bool_t crop)
{
   if (fStop) {
      return;
   }

   Int_t x0, y0, x, y, xx, yy, w, h;
   Window_t c;

   if (!fClient->GetRoot()->InheritsFrom(TGCompositeFrame::Class()) ||
       !fClient->IsEditable()) {
      return;
   }

   TGCompositeFrame *comp = 0;
   Bool_t fromGrab = kFALSE;
   TGFrame *frame = fPimpl->fGrab;

   if (fBuilder && crop) {
      comp = fBuilder->FindEditableMdiFrame(fClient->GetRoot());
   } else {
      comp = (TGCompositeFrame*)fClient->GetRoot();
   }

   if (frame && !CanChangeLayout((TGWindow*)frame->GetParent())) {
      frame = GetMovableParent(frame);

      if (!frame) {
         TString str = fPimpl->fGrab->ClassName();
         str += "::";
         str += fPimpl->fGrab->GetName();
         str += " cannot be deleted";

         if (fBuilder) {
            fBuilder->UpdateStatusBar(str.Data());
         }
         return;
      }
   }

   // prepare to crop grabbed frame
   if (frame && !fLassoDrawn && crop) {
      gVirtualX->TranslateCoordinates(frame->GetId(),
                                      fClient->GetDefaultRoot()->GetId(),
                                      -2, -2,
                                      fPimpl->fX0, fPimpl->fY0, c);

      fPimpl->fX = fPimpl->fX0 + frame->GetWidth()+4;
      fPimpl->fY = fPimpl->fY0 + frame->GetHeight()+4;
      fromGrab = kTRUE;
   }

   x0 = fPimpl->fX0; y0 = fPimpl->fY0;
   x  = fPimpl->fX;  y  = fPimpl->fY;
   if (comp) {
      gVirtualX->TranslateCoordinates(fClient->GetDefaultRoot()->GetId(),
                                      comp->GetId(),
                                      fPimpl->fX, fPimpl->fY, x, y, c);
      gVirtualX->TranslateCoordinates(fClient->GetDefaultRoot()->GetId(),
                                      comp->GetId(),
                                      fPimpl->fX0, fPimpl->fY0, x0, y0, c);
   }

   xx = x0; yy = y0;
   x0 = TMath::Min(xx, x); x = TMath::Max(xx, x);
   y0 = TMath::Min(yy, y); y = TMath::Max(yy, y);
   w = x - x0;
   h = y - y0;

   if (fLassoDrawn || fromGrab) {
      if (comp) {
         TIter next(comp->GetList());
         TGFrameElement *el;

         while ((el = (TGFrameElement*)next())) {
            TGFrame *fr = el->fFrame;

            if ((fr->GetX() >= x0) && (fr->GetY() >= y0) &&
                (fr->GetX() + (Int_t)fr->GetWidth() <= x) &&
                (fr->GetY() + (Int_t)fr->GetHeight() <= y)) {
               if (!crop) {
                  DeleteFrame(fr);
               } else {
                  fr->Move(fr->GetX() - x0, fr->GetY() - y0);
               }
            } else {
               if (crop) {
                  DeleteFrame(fr);
               }
            }
         }
         if (crop && comp) {
            gVirtualX->TranslateCoordinates(comp->GetId(), comp->GetParent()->GetId(),
                                             x0, y0, xx, yy, c);

            comp->MoveResize(xx, yy, w, h);

            if (comp->GetParent()->InheritsFrom(TGMdiDecorFrame::Class())) {
               TGMdiDecorFrame *decor = (TGMdiDecorFrame *)comp->GetParent();

               gVirtualX->TranslateCoordinates(decor->GetId(), decor->GetParent()->GetId(),
                                               xx, yy, xx, yy, c);

               Int_t b = 2 * decor->GetBorderWidth();
               decor->MoveResize(xx, yy, comp->GetWidth() + b,
                                 comp->GetHeight() + b + decor->GetTitleBar()->GetDefaultHeight());
            }
         }
      }
   } else { //  no lasso drawn -> delete selected frame
      if (frame)
         DeleteFrame(frame);
      UngrabFrame();
      ChangeSelected(0);   //update editors
   }
   SetLassoDrawn(kFALSE);

   if (fBuilder) {
      //fBuilder->Update();
      fBuilder->UpdateStatusBar(crop ? "Crop action performed" : "Delete action performed");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Delete frame

void TGuiBldDragManager::DeleteFrame(TGFrame *frame)
{
   if (fStop || !frame) {
      return;
   }

   // remove the frame from the list tree and reset the editor...
   fEditor->RemoveFrame(frame);

   frame->UnmapWindow();

   TGCompositeFrame *comp = 0;

   if (frame->GetParent()->InheritsFrom(TGCompositeFrame::Class())) {
      comp = (TGCompositeFrame*)frame->GetParent();
   }

   if (comp) {
      comp->RemoveFrame(frame);
   }

   if (frame == fPimpl->fGrab) {
      UngrabFrame();
   }

   fClient->UnregisterWindow(frame);

   // mem.leak paid for robustness (with possibility "undelete")
   frame->ReparentWindow(fClient->GetDefaultRoot());
}

////////////////////////////////////////////////////////////////////////////////
/// Handle cut action

void TGuiBldDragManager::HandleCut()
{
   if (fStop || !fPimpl->fGrab) {
      return;
   }

   //
   fPimpl->fGrab = GetMovableParent(fPimpl->fGrab);
   HandleCopy();
   DeleteFrame(fPimpl->fGrab);
   ChangeSelected(0);   //update editors
}

////////////////////////////////////////////////////////////////////////////////
/// Handle copy. This method is also used by SaveFrame method.
/// In later case  brk_layout == kFALSE

void TGuiBldDragManager::HandleCopy(Bool_t brk_layout)
{
   if (fStop || !fPimpl->fGrab) {
      return;
   }

   TGMainFrame *tmp = new TGMainFrame(fClient->GetDefaultRoot(),
                                      fPimpl->fGrab->GetWidth(),
                                      fPimpl->fGrab->GetHeight());

   // save coordinates
   Int_t x0 = fPimpl->fGrab->GetX();
   Int_t y0 = fPimpl->fGrab->GetY();

   // save parent name
   TString name = fPimpl->fGrab->GetParent()->GetName();

   ((TGWindow*)fPimpl->fGrab->GetParent())->SetName(tmp->GetName());

   fPimpl->fGrab->SetX(0);
   fPimpl->fGrab->SetY(0);

   TGFrameElement *fe = fPimpl->fGrab->GetFrameElement();

   if (fe) {
      tmp->GetList()->Add(fe);
   }

   tmp->SetLayoutBroken(brk_layout);

   if (!brk_layout) { //save frame
      tmp->SetMWMHints(kMWMDecorAll, kMWMFuncAll, kMWMInputFullApplicationModal);
      tmp->SetWMSize(tmp->GetWidth(), tmp->GetHeight());
      tmp->SetWMSizeHints(tmp->GetDefaultWidth(), tmp->GetDefaultHeight(), 10000, 10000, 0, 0);
      const char *short_name = gSystem->BaseName(fPasteFileName.Data());
      tmp->SetWindowName(short_name);
      tmp->SetIconName(short_name);
      tmp->SetClassHints(short_name, short_name);
      // some problems here under win32
      if (gVirtualX->InheritsFrom("TGX11")) tmp->SetIconPixmap("bld_rgb.xpm");
   }
   Bool_t quite =  brk_layout || (fPasteFileName == fTmpBuildFile);
   tmp->SaveSource(fPasteFileName.Data(), quite ? "keep_names quiet" : "keep_names");
   tmp->GetList()->Remove(fe);

   fPimpl->fGrab->SetX(x0);
   fPimpl->fGrab->SetY(y0);

   ((TGWindow*)fPimpl->fGrab->GetParent())->SetName(name.Data());

   if (fBuilder) {
      TString str = fPimpl->fGrab->ClassName();
      str += "::";
      str += fPimpl->fGrab->GetName();
      str += " copied to clipboard";
      fBuilder->UpdateStatusBar(str.Data());
   }

   delete tmp;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle paste action.

void TGuiBldDragManager::HandlePaste()
{
   if (fStop) {
      return;
   }

   Int_t xp = 0;
   Int_t yp = 0;

   if (gSystem->AccessPathName(fPasteFileName.Data())) {
      return;
   }

   fPasting = kTRUE;
   gROOT->Macro(fPasteFileName.Data());

   Window_t c;
   TGFrame *root = (TGFrame*)fClient->GetRoot();

   if (!fPimpl->fReplaceOn) {
      gVirtualX->TranslateCoordinates(fClient->GetDefaultRoot()->GetId(),
                                      root->GetId(),
                                      fPimpl->fX0, fPimpl->fY0, xp, yp, c);
      ToGrid(xp, yp);

      // fPasteFrame is defined in TVirtualDragManager.h
      // fPasteFrame is a TGMainFrame consisting "the frame to paste"
      // into the editable frame (aka fClient->GetRoot())

      if (fPasteFrame) {
         TGMainFrame *main = (TGMainFrame*)fPasteFrame;
         TGFrame *paste = ((TGFrameElement*)main->GetList()->First())->fFrame;

         UInt_t w = paste->GetWidth();
         UInt_t h = paste->GetHeight();

         if (xp + w > root->GetWidth()) {
            w = root->GetWidth() - xp -1;
         }
         if (yp + h > root->GetHeight()) {
            h = root->GetHeight() - yp -1;
         }

         paste->Resize(w, h);
         fPasteFrame->Move(xp, yp);
         fPimpl->fGrab = fPasteFrame;
         HandleReturn(1);  // drop
      }
   }

   fPasting = kFALSE;

   if (fBuilder) {
      fBuilder->UpdateStatusBar("Paste action performed");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Replace frame (doesn't work yet properly)

void TGuiBldDragManager::DoReplace(TGFrame *frame)
{
   if (fStop || !frame || !fPimpl->fGrab || !fPimpl->fReplaceOn) {
      return;
   }

   Int_t w = fPimpl->fGrab->GetWidth();
   Int_t h = fPimpl->fGrab->GetHeight();
   Int_t x = fPimpl->fGrab->GetX();
   Int_t y = fPimpl->fGrab->GetY();

   if (fBuilder) {
      TString str = fPimpl->fGrab->ClassName();
      str += "::";
      str += fPimpl->fGrab->GetName();
      str += " replaced by ";
      str += frame->ClassName();
      str += "::";
      str += frame->GetName();
      fBuilder->UpdateStatusBar(str.Data());
   }

   TGFrameElement *fe = fPimpl->fGrab->GetFrameElement();

   if (fe) {
      fe->fFrame = 0;
      fPimpl->fGrab->DestroyWindow();
      delete fPimpl->fGrab;
      fPimpl->fGrab = 0;

      fe->fFrame = frame;
      frame->MoveResize(x, y, w, h);
      frame->MapRaised();
      frame->SetFrameElement(fe);
   }

   SelectFrame(frame);
   fPimpl->fReplaceOn = kFALSE;

   TGWindow *root = (TGWindow *)fClient->GetRoot();
   root->SetEditable(kFALSE);
   DoRedraw();
   root->SetEditable(kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Handle replace

void TGuiBldDragManager::HandleReplace()
{
   if (fStop || !fPimpl->fGrab) {
      return;
   }

   fPimpl->fReplaceOn = kTRUE;
   TGFrame *frame = 0;

   if (fBuilder && fBuilder->IsExecutable())  {
      frame = (TGFrame *)fBuilder->ExecuteAction();
   } else {
      HandlePaste();
      frame = fPasteFrame;
   }
   DoReplace(frame);
   fPimpl->fReplaceOn = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a frame which is the same as currently editted frame

void TGuiBldDragManager::CloneEditable()
{
   if (fStop) {
      return;
   }

   TString tmpfile = gSystem->TempDirectory();
   char *s = gSystem->ConcatFileName(tmpfile.Data(), TString::Format("tmp%d.C",
                                     gRandom->Integer(100)));
   Save(s);
   gROOT->Macro(s);
   gSystem->Unlink(s);
   delete [] s;

   if (fClient->GetRoot()->InheritsFrom(TGFrame::Class())) {
      TGFrame *f = (TGFrame *)fClient->GetRoot();
      f->Resize(f->GetWidth() + 10, f->GetHeight() + 10);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Save an editted frame to the file

Bool_t TGuiBldDragManager::Save(const char *file)
{
   if (fStop || !fClient->GetRoot() || !fClient->IsEditable()) {
      return kFALSE;
   }

   TGMainFrame *main = (TGMainFrame*)fClient->GetRoot()->GetMainFrame();
   TGWindow *root = (TGWindow*)fClient->GetRoot();
   TString fname = file;

   root->SetEditable(kFALSE);

   static TImage *img = 0;

   if (!img) {
      img = TImage::Create();
   }
   img->FromWindow(main->GetId());

   if (!file || !file[0]) {
      static TString dir(".");
      static Bool_t overwr = kFALSE;
      TGFileInfo fi;

      fi.fFileTypes = gSaveMacroTypes;
      fi.SetIniDir(dir);
      fi.fOverwrite = overwr;
      new TGFileDialog(fClient->GetDefaultRoot(), this, kFDSave, &fi);

      if (!fi.fFilename) goto out;
      dir = fi.fIniDir;
      overwr = fi.fOverwrite;
      fname = gSystem->BaseName(gSystem->UnixPathName(fi.fFilename));
   }

   if (fname.EndsWith(".C")) {
      main->SetMWMHints(kMWMDecorAll, kMWMFuncAll, kMWMInputFullApplicationModal);
      main->SetWMSize(main->GetWidth(), main->GetHeight());
      main->SetWMSizeHints(main->GetDefaultWidth(), main->GetDefaultHeight(), 10000, 10000, 0, 0);
      main->SetWindowName(fname.Data());
      main->SetIconName(fname.Data());
      main->SetClassHints(fname.Data(), fname.Data());
      // some problems here under win32
      if (gVirtualX->InheritsFrom("TGX11")) main->SetIconPixmap("bld_rgb.xpm");
      main->SaveSource(fname.Data(), file ? "keep_names quiet" : "keep_names");

      fBuilder->AddMacro(fname.Data(), img);

   } else {
      Int_t retval;
      TString msg = TString::Format("file (%s) must have extension .C", fname.Data());

      new TGMsgBox(fClient->GetDefaultRoot(), main, "Error...", msg.Data(),
                   kMBIconExclamation, kMBRetry | kMBCancel, &retval);

      if (retval == kMBRetry) {
         return Save();
      }
   }

out:
   main->RaiseWindow();
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Save composite frame as macro

Bool_t TGuiBldDragManager::SaveFrame(const char *file)
{
   if (fStop || !fClient->GetRoot() || !fClient->IsEditable() ||
       !fPimpl->fGrab || !fPimpl->fGrab->InheritsFrom(TGCompositeFrame::Class())) {
      return kFALSE;
   }

   TString fname = file;

   TGFrame *frame = fPimpl->fGrab;
   SetEditable(kFALSE);

   static TImage *img = 0;

   if (!img) {
      img = TImage::Create();
   }
   img->FromWindow(frame->GetId());

   static TString dir(".");
   static Bool_t overwr = kFALSE;

   TString sav = fPasteFileName;

   if (!file) {
      TGFileInfo fi;

      fi.fFileTypes = gSaveMacroTypes;
      fi.SetIniDir(dir);
      fi.fOverwrite = overwr;
      new TGFileDialog(fClient->GetDefaultRoot(), frame, kFDSave, &fi);

      if (!fi.fFilename) {
         goto out;
      }

      dir = fi.fIniDir;
      overwr = fi.fOverwrite;
      fname = gSystem->BaseName(gSystem->UnixPathName(fi.fFilename));
   }

   if (fname.EndsWith(".C")) {
      fPasteFileName = fname;
      fPimpl->fGrab = frame;
      fStop = kFALSE;
      TGFrameElement *fe = frame->GetFrameElement();

      if (!fe) { // should never happen
         fe = new TGFrameElement();
         fe->fFrame = frame;
         fe->fState = kIsMapped;
         frame->SetFrameElement(fe);
         TGCompositeFrame *comp = (TGCompositeFrame*)frame->GetParent();
         comp->GetList()->Add(fe);
      }
      delete fe->fLayout;
      fe->fLayout = new TGLayoutHints(kLHintsExpandX | kLHintsExpandY);

      HandleCopy(kFALSE);
      fStop = kTRUE;

      fBuilder->AddMacro(fname.Data(), img);
   } else {
      Int_t retval;
      TString msg = TString::Format("file (%s) must have extension .C", fname.Data());

      new TGMsgBox(fClient->GetDefaultRoot(), frame, "Error...", msg.Data(),
                   kMBIconExclamation, kMBRetry | kMBCancel, &retval);

      if (retval == kMBRetry) {
         return SaveFrame();
      }
   }

out:
   fPasteFileName = sav;
   return kTRUE;
}
/*
////////////////////////////////////////////////////////////////////////////////
/// Not used yet. Return 0 if all child frames are inside area x,y, w,h

static Int_t canResize(TGFrame *frame, Int_t x, Int_t y, UInt_t &w, UInt_t &h)
{
   if (frame->InheritsFrom(TGCompositeFrame::Class())) return 0;

   TGCompositeFrame *comp = (TGCompositeFrame*)frame;

   TIter next(comp->GetList());
   TGFrameElement *fe;
   Int_t d = gGuiBldDragManager->GetGridStep();
   Int_t ret = 0;

   while ((fe = (TGFrameElement*)next())) {
      if (x + fe->fFrame->GetX() + fe->fFrame->GetWidth() > w) {
         w = fe->fFrame->GetX() + x  + fe->fFrame->GetWidth();
         ret |= 4;
      }
      if (y + fe->fFrame->GetY() + fe->fFrame->GetHeight() > h) {
         h = fe->fFrame->GetY() + y + fe->fFrame->GetHeight();
         ret |= 8;
      }
   }
   return ret;
}
*/

////////////////////////////////////////////////////////////////////////////////
/// handle resize

void TGuiBldDragManager::DoResize()
{
   if (fStop || !fClient->IsEditable()) {
      return;
   }

   TGFrame *fr = fPimpl->fGrab;

   if (!fr || IsFixedSize(fr) ||
        IsFixedLayout((TGWindow*)fr->GetParent())) {

      fr = (TGFrame*)GetResizableParent(fr);

      if (!fr ) {
         return;
      }
   }

   TGCompositeFrame *comp = 0;

   if (fr->InheritsFrom(TGCompositeFrame::Class())) {
      comp = (TGCompositeFrame*)fr;
   }

   Window_t c;
   Int_t x = fPimpl->fX;
   Int_t y = fPimpl->fY;
   UInt_t w = 0;
   UInt_t h = 0;
   UInt_t wp = ((TGFrame*)fr->GetParent())->GetWidth() - 2;
   UInt_t hp = ((TGFrame*)fr->GetParent())->GetHeight() - 2;

   gVirtualX->TranslateCoordinates(fClient->GetDefaultRoot()->GetId(),
                                   fr->GetId(), x, y, x, y, c);

   ToGrid(x, y);
   HighlightCompositeFrame(fr->GetParent()->GetId());

   switch (fPimpl->fResizeType) {
      case kTopLeft:
         if ((((int)fr->GetWidth() > x) || (x < 0)) &&
             (((int)fr->GetHeight() > y) || (y < 0))) {

            if (fr->GetY() + y < 2) {
               y = 2 - fr->GetY();
            }
            if (fr->GetX() + x < 2) {
               x = 2 - fr->GetX();
            }
            h = fr->GetHeight() - y;
            w = fr->GetWidth() - x;
            x = fr->GetX() + x;
            y = fr->GetY() + y;

            if (!IsFixedH(fr) && !IsFixedW(fr)) {
               fr->MoveResize(x, y, w, h);
               break;
            }
            if (IsFixedH(fr)) {
               fr->MoveResize(x, fr->GetY(), w, fr->GetDefaultHeight());
               break;
            }
            if (IsFixedW(fr)) {
               fr->MoveResize(fr->GetX(), y, fr->GetDefaultWidth(), h);
               break;
            }
         }
         break;
      case kTopRight:
         if ((x > 0) && (((int)fr->GetHeight() > y) || (y < 0))) {

            if (fr->GetY() + y < 2) {
               y = 2 - fr->GetY();
            }
            h = fr->GetHeight() - y;

            if (IsFixedW(fr)) {
               w =  fr->GetDefaultWidth();
            } else {
               w = fr->GetX() + x > Int_t(wp) ? wp - fr->GetX() : UInt_t(x);
            }
            x = fr->GetX();
            y = fr->GetY() + y;

            if (!IsFixedH(fr)) {
               fr->MoveResize(x, y, w, h);
            } else {
               fr->Resize(x, fr->GetDefaultHeight());
            }
         }
         break;
      case kTopSide:
         if (((int)fr->GetHeight() > y) || (y < 0)) {
            if (IsFixedH(fr)) {
               break;
            }

            if (fr->GetY() + y < 2) {
               y = 2 - fr->GetY();
            }
            h = fr->GetHeight() - y;
            w = fr->GetWidth();
            x = fr->GetX();
            y = fr->GetY() + y;

            fr->MoveResize(x, y, w, h);
         }
         break;
      case kBottomLeft:
         if ((((int)fr->GetWidth() > x) || (x < 0)) && (y > 0)) {

            if (fr->GetX() + x < 2) {
               x = 2 - fr->GetX();
            }
            h = fr->GetY() + y > Int_t(hp) ? hp - fr->GetY() : UInt_t(y);
            w = fr->GetWidth() - x;
            x = fr->GetX() + x;

            if (!IsFixedH(fr) && !IsFixedW(fr)) {
               fr->MoveResize(x, fr->GetY(), w, h);
               break;
            }
            if (IsFixedH(fr)) {
               fr->MoveResize(x, fr->GetY(), w, fr->GetDefaultHeight());
               break;
            }
            if (IsFixedW(fr)) {
               fr->MoveResize(fr->GetX(), fr->GetY(),
                              fr->GetDefaultWidth(), h);
               break;
            }
         }
         break;
      case kBottomRight:
         if ((x > 0) && (y > 0)) {
            w = !IsFixedW(fr) ? UInt_t(x) : fr->GetDefaultWidth();
            h = !IsFixedH(fr) ? UInt_t(y) : fr->GetDefaultHeight();

            h = fr->GetY() + h > hp ? hp - fr->GetY() : h;
            w = fr->GetX() + w > wp ? wp - fr->GetX() : w;

            //canResize(comp, 0, 0, w, h);
            fr->Resize(w, h);
         }
         break;
      case kBottomSide:
         if (y > 0) {
            if (IsFixedH(fr)) {
               break;
            }

            w = fr->GetWidth();
            h = fr->GetY() + y > (Int_t)hp ? hp - fr->GetY() : UInt_t(y);

            //canResize(comp, 0, 0, w, h);
            fr->Resize(w, h);
         }
         break;
      case kLeftSide:
         if ((int)fr->GetWidth() > x ) {
            if (IsFixedW(fr)) {
               break;
            }

            if (fr->GetX() + x < 2) {
               x = 2 - fr->GetX();
            }
            w = fr->GetWidth() - x;
            h = fr->GetHeight();
            y = fr->GetY();
            x = fr->GetX() + x;

            //canResize(comp, x, y, w, h);
            fr->MoveResize(x, y, w, h);
         }
         break;
      case kRightSide:
         if (x > 0) {
            if (IsFixedW(fr)) {
               break;
            }

            h = fr->GetHeight();
            w = fr->GetX() + x > (Int_t)wp ? wp - fr->GetX() : UInt_t(x);
            //canResize(comp, 0, 0, w, h);
            fr->Resize(w, h);
         }
         break;
      default:
         break;
   }
   if (comp && (!comp->IsLayoutBroken() || IsFixedLayout(comp))) {
      layoutFrame(comp);
   }

   gVirtualX->SetCursor(fClient->GetRoot()->GetId(),
                        gVirtualX->CreateCursor(fPimpl->fResizeType));
   w = fr->GetWidth();
   h = fr->GetHeight();

   if (fBuilder) {
      TString str = fr->ClassName();
      str += "::";
      str += fr->GetName();
      str += " resized   ";
      str += TString::Format("(%d x %d)", w, h);
      fBuilder->UpdateStatusBar(str.Data());
   }

   fClient->NeedRedraw(fr, kTRUE);
   DoRedraw();
   fEditor->ChangeSelected(fr); //to update the geometry frame after drag resize
}

////////////////////////////////////////////////////////////////////////////////
/// Handle move

void TGuiBldDragManager::DoMove()
{
   if (fStop || !fPimpl->fGrab || !fClient->IsEditable()) {
      return;
   }

   TGWindow *parent = (TGWindow*)fPimpl->fGrab->GetParent();

   // do not remove frame from fixed layout or non-editable parent
   if (IsFixedLayout(parent) || IsEditDisabled(parent)) {
      return;
   }

   Int_t x = fPimpl->fX - fPimpl->fXf;
   Int_t y = fPimpl->fY - fPimpl->fYf;

   static Int_t qq;
   static UInt_t w = 0;
   static UInt_t h = 0;

   if (w == 0) {
      gVirtualX->GetWindowSize(gVirtualX->GetDefaultRootWindow(), qq, qq, w, h);
   }

   //
   Bool_t move = (x > 0) && (y > 0) && ((x + fPimpl->fGrab->GetWidth()) < (w - 0)) &&
                 ((y + fPimpl->fGrab->GetHeight()) < (h - 30));


   // we are out of "win32 world"
   if (!move && !gVirtualX->InheritsFrom("TGX11")) {
      EndDrag();
      return;
   }

   fPimpl->fGrab->Move(x, y);

   if (fBuilder) {
      //fBuilder->Update();
      TString str = fPimpl->fGrab->ClassName();
      str += "::";
      str += fPimpl->fGrab->GetName();
      str += " is moved to absolute position   ";
      str += TString::Format("(%d , %d)", x, y);
      fBuilder->UpdateStatusBar(str.Data());
   }

   CheckTargetUnderGrab();
}

////////////////////////////////////////////////////////////////////////////////
/// Return a pointer to the parent mdi frame

TGFrame *TGuiBldDragManager::FindMdiFrame(TGFrame *in)
{
   if (fStop || !in) {
      return 0;
   }

   TGFrame *p = in;

   while (p && (p != fClient->GetDefaultRoot()) &&
         !p->InheritsFrom(TGMainFrame::Class())) {
      if (p->InheritsFrom(TGMdiFrame::Class())) {
         return p;
      }
      p = (TGFrame*)p->GetParent();
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Raise guibuilder's mdi frame.

void TGuiBldDragManager::RaiseMdiFrame(TGFrame *comp)
{
   if (fStop || !comp) {
      return;
   }

   if (comp && comp->InheritsFrom(TGMdiFrame::Class()) && fBuilder) {
      TGFrame *mdi = fBuilder->FindEditableMdiFrame(comp);
      if (mdi) {
         // dragged frame is taken from some main frame
         //if (fPimpl->fGrab && fClient->GetRoot()->InheritsFrom(TGMainFrame::Class())) {
         //   fBuilder->MapRaised();
         //}
      }
      if (fBuilder->GetMdiMain()->GetCurrent() != comp) {
         fBuilder->GetMdiMain()->SetCurrent((TGMdiFrame*)comp);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Look for the drop target under grabbed/selected frame while moving

void TGuiBldDragManager::CheckTargetUnderGrab()
{
   if (fStop || !fPimpl->fGrab ) {
      return;
   }

   Int_t x = fPimpl->fGrab->GetX();
   Int_t y = fPimpl->fGrab->GetY();
   UInt_t w = fPimpl->fGrab->GetWidth();
   UInt_t h = fPimpl->fGrab->GetHeight();

   Bool_t ok = CheckTargetAtPoint(x - 1, y - 1);

   if (!ok) {
      ok = CheckTargetAtPoint(x + w + 1, y + h + 1);
   }

   if (!ok) {
      ok = CheckTargetAtPoint(x + w + 1, y - 1);
   }

   if (!ok) {
      ok = CheckTargetAtPoint(x - 1, y + h + 1);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Helper. Look for the drop target under grabbed/selected frame while moving.

Bool_t TGuiBldDragManager::CheckTargetAtPoint(Int_t x, Int_t y)
{
   if (fStop || !fPimpl->fGrab) {
      return kFALSE;
   }

   UInt_t ww = fPimpl->fGrab->GetWidth();
   UInt_t hh = fPimpl->fGrab->GetHeight();
   Bool_t ret = kFALSE;
   Window_t c;
   TGWindow *win = 0;

   Window_t w = GetWindowFromPoint(x, y);

   if (w && (w != gVirtualX->GetDefaultRootWindow())) {
      win = fClient->GetWindowById(w);
      TGCompositeFrame *comp = 0;

      if (!win) {
         goto out;
      }

      if (win->InheritsFrom(TGCompositeFrame::Class())) {
         comp = (TGCompositeFrame *)win;
      } else if (win->GetParent() != fClient->GetDefaultRoot()) {
         comp = (TGCompositeFrame *)win->GetParent();
      }

      if (comp) {
         gVirtualX->TranslateCoordinates(fClient->GetDefaultRoot()->GetId(),
                                         comp->GetId(), x, y, x, y, c);

         RaiseMdiFrame(comp);

         if ((comp != fPimpl->fGrab) && (x >= 0) && (y >= 0) &&
             (x + ww <= comp->GetWidth()) &&
             (y + hh <= comp->GetHeight())) {

            if (comp != fTarget) {
               comp->HandleDragEnter(fPimpl->fGrab);

               if (fTarget) fTarget->HandleDragLeave(fPimpl->fGrab);

               else Snap2Grid();
            } else {
               if (fTarget) {
                  fTarget->HandleDragMotion(fPimpl->fGrab);
               }
            }

            fTarget = comp;
            fTargetId = comp->GetId();
            ret = kTRUE;
            return ret;

         } else {
            if (fTarget) {
               fTarget->HandleDragLeave(fPimpl->fGrab);
            }
            fTarget = 0;
            fTargetId = 0;
         }
      }
   }

out:
   if (fTarget) {
      fTarget->HandleDragLeave(fPimpl->fGrab);
   }

   if (!w || !win) {
      fTarget = 0;
      fTargetId = 0;
   }
   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle motion event

Bool_t TGuiBldDragManager::HandleMotion(Event_t *event)
{
   if (fStop) {
      return kFALSE;
   }

   static Long64_t was = gSystem->Now();
   static Int_t gy = event->fYRoot;
   static Int_t gx = event->fXRoot;

   Long64_t now = gSystem->Now();

   if ((now-was < 100) || !(event->fState & kButton1Mask) ||
       ((event->fYRoot == gy) && (event->fXRoot == gx))) {
      return kFALSE;
   }

   was = now;
   gy = event->fYRoot;
   gx = event->fXRoot;

   if (!fDragging) {
      if (fMoveWaiting && ((TMath::Abs(fPimpl->fX - event->fXRoot) > 10) ||
          (TMath::Abs(fPimpl->fY - event->fYRoot) > 10))) {

         return StartDrag(fSource, event->fXRoot, event->fYRoot);
      }
   } else {
      fPimpl->fX = event->fXRoot;
      fPimpl->fY = event->fYRoot;

      switch (fDragType) {
         case kDragLasso:
            DrawLasso();
            fSelectionIsOn = event->fState & kKeyShiftMask;
            break;
         case kDragMove:
         case kDragCopy:
         case kDragLink:
            DoMove();
            break;
         case kDragResize:
            DoResize();
            break;
         default:
            break;
      }
   }
   DrawGrabRectangles();
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Put created frame at position of the last mouse click

void TGuiBldDragManager::PlaceFrame(TGFrame *frame, TGLayoutHints *hints)
{
   Int_t x0, y0, x, y;
   Window_t c;

   if (fStop || !frame || !fClient->IsEditable()) {
      return;
   }

   frame->MapSubwindows();
   TGFrame *root = (TGFrame*)fClient->GetRoot();

   gVirtualX->TranslateCoordinates(fClient->GetDefaultRoot()->GetId(),
                                   root->GetId(),
                                   fPimpl->fX0 , fPimpl->fY0, x0, y0, c);
   gVirtualX->TranslateCoordinates(fClient->GetDefaultRoot()->GetId(),
                                   root->GetId(),
                                   fPimpl->fX , fPimpl->fY, x, y, c);

   ToGrid(x, y);
   ToGrid(x0, y0);

   UInt_t w = TMath::Abs(x - x0);
   UInt_t h = TMath::Abs(y - y0);
   x = x > x0 ? x0 : x;
   y = y > y0 ? y0 : y;

   // do not create frame with size smaller when default size
   w = w < frame->GetDefaultWidth() + 2 ? frame->GetDefaultWidth() + 2 : w;
   h = h < frame->GetDefaultHeight() + 2 ? frame->GetDefaultHeight() + 2 : h;

   // do not create frame out of editable space
   x = x + w > root->GetWidth() ? Int_t(root->GetWidth() - w) : x;
   y = y + h > root->GetHeight() ? Int_t(root->GetHeight() - h) : y;

   frame->Move(x, y);

   UInt_t grid = GetGridStep();

   if (IsFixedW(frame) || IsFixedH(frame) || IsFixedSize(frame)) {
      w = IsFixedW(frame) ? frame->GetDefaultWidth() : w;
      h = IsFixedH(frame) ? frame->GetDefaultHeight() : h;
      frame->Resize(w < grid ? grid : w, h < grid ? grid : h);
   } else {
      if (frame->InheritsFrom(TGVerticalFrame::Class())) {
         frame->Resize(w < grid ? 15*grid : w, h < grid ?  30*grid : h);
      } else if (frame->InheritsFrom(TGHorizontalFrame::Class())) {
         frame->Resize(w < grid ? 30*grid : w, h < grid ?  15*grid : h);
      }
      else frame->Resize(w < 2*grid ? 2*grid : w, h < 2*grid ?  2*grid : h);
   }

   frame->MapRaised();
   frame->SetCleanup(kDeepCleanup);
   frame->AddInput(kButtonPressMask);

   if (fClient->GetRoot()->InheritsFrom(TGCompositeFrame::Class())) {
      TGCompositeFrame *edit = (TGCompositeFrame*)fClient->GetRoot();
      edit->SetCleanup(kDeepCleanup);
      ReparentFrames(frame, edit);
      frame->MapRaised();
      //edit->SetLayoutBroken();
      UInt_t g = 2;
      // temporary hack for status bar
      if (frame->InheritsFrom("TGStatusBar")) {
         edit->AddFrame(frame, new TGLayoutHints(kLHintsBottom | kLHintsExpandX));
      }
      else {
         edit->AddFrame(frame, hints ? hints : new TGLayoutHints(kLHintsNormal, g, g, g, g));
      }

      if (hints && !edit->IsLayoutBroken()) {
         edit->GetLayoutManager()->Layout();
      } else {
         edit->Layout();
      }
   }
   if (fBuilder) {
      TString str = frame->ClassName();
      str += "::";
      str += frame->GetName();
      str += " created";
      fBuilder->UpdateStatusBar(str.Data());
   }

   if (frame->InheritsFrom(TGCanvas::Class())) {
      frame = ((TGCanvas*)frame)->GetContainer();
   }

   SelectFrame(frame);

}
////////////////////////////////////////////////////////////////////////////////
/// Draw lasso for allocation new object

void TGuiBldDragManager::DrawLasso()
{
   if (fStop || !fClient->IsEditable()) {
      return;
   }

   UngrabFrame();

   Int_t x0, y0, x, y;
   Window_t c;
   TGFrame *root = (TGFrame*)fClient->GetRoot();

   gVirtualX->TranslateCoordinates(fClient->GetDefaultRoot()->GetId(), root->GetId(),
                                   fPimpl->fX0 , fPimpl->fY0, x0, y0, c);
   gVirtualX->TranslateCoordinates(fClient->GetDefaultRoot()->GetId(), root->GetId(),
                                   fPimpl->fX , fPimpl->fY, x, y, c);

   UInt_t w, h;
   Bool_t xswap = kFALSE;
   Bool_t yswap = kFALSE;

   // check limits

   if ((x == x0) || ( y==y0 )) return; //lasso is not rectangle -> do not draw it

   if (x > x0) {
      x0 = x0 < 0 ? 0 : x0;
      w = x - x0;
   } else {
      x = x < 0 ? 0 : x;
      w = x0 - x;
      x0 = x;
      xswap = kTRUE;
   }

   if (y > y0) {
      y0 = y0 < 0 ? 0 : y0;
      h = y - y0;
   } else {
      y = y < 0 ? 0 : y;
      h = y0 - y;
      y0 = y;
      yswap = kTRUE;
   }

   w = x0 + w > root->GetWidth() ? root->GetWidth() - x0 : w;
   h = y0 + h > root->GetHeight() ? root->GetHeight() - y0 : h;
   x = x0 + w;
   y = y0 + h;

   ToGrid(x, y);
   ToGrid(x0, y0);

   // correct fPimpl->fX0 , fPimpl->fY0 , fPimpl->fX , fPimpl->fY
   gVirtualX->TranslateCoordinates(root->GetId(), fClient->GetDefaultRoot()->GetId(),
                                   xswap ? x : x0, yswap ? y : y0,
                                   fPimpl->fX0 , fPimpl->fY0,  c);
   gVirtualX->TranslateCoordinates(root->GetId(), fClient->GetDefaultRoot()->GetId(),
                                   xswap ? x0 : x, yswap ? y0 : y,
                                   fPimpl->fX , fPimpl->fY,  c);
   DoRedraw();

   gVirtualX->DrawRectangle(fClient->GetRoot()->GetId(),
                            GetBlackGC()(), x0, y0, w, h);
   gVirtualX->DrawRectangle(fClient->GetRoot()->GetId(),
                            GetBlackGC()(), x0+1, y0+1, w-2, h-2);

   gVirtualX->SetCursor(fId, gVirtualX->CreateCursor(kCross));
   gVirtualX->SetCursor(fClient->GetRoot()->GetId(), gVirtualX->CreateCursor(kCross));

   SetLassoDrawn(kTRUE);
   root->RequestFocus();

   if (fBuilder) {
      TString str = "Lasso drawn. Align frames inside or presss Return key to grab frames.";
      fBuilder->UpdateStatusBar(str.Data());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Handle client message

Bool_t TGuiBldDragManager::HandleClientMessage(Event_t *event)
{
   if (fStop) {
      return kFALSE;
   }

   if ((event->fFormat == 32) && ((Atom_t)event->fUser[0] == gWM_DELETE_WINDOW) &&
       (event->fHandle != gROOT_MESSAGE)) {

      if (fPimpl->fPlane && (fPimpl->fPlane->GetId() == event->fWindow)) {
         fPimpl->fPlane = 0;
      }

      TGWindow *root = (TGWindow*)fClient->GetRoot();
      if (!root || (root == fClient->GetDefaultRoot())) {
         SetEditable(kFALSE);
         return kTRUE;
      }
      TGMainFrame *main = (TGMainFrame*)root->GetMainFrame();

      if (event->fWindow == main->GetId()) {
         if (main != fBuilder) {
            if (fEditor && !fEditor->IsEmbedded()) {
               delete fEditor;
               fEditor = 0;
            }

            SetEditable(kFALSE);
            return kTRUE;
         }

         delete fFrameMenu;
         fFrameMenu =0;

         delete fLassoMenu;
         fLassoMenu = 0;

         delete fPimpl->fGrid;
         fPimpl->fGrid = 0;
         Reset1();

      } else if (fBuilder && (event->fWindow == fBuilder->GetId())) {
         fBuilder->CloseWindow();

      } else if (fEditor && (event->fWindow == fEditor->GetMainFrame()->GetId())) {
         TQObject::Disconnect(fEditor);
         fEditor = 0;
      }

      // to avoid segv. stop editting
      SetEditable(kFALSE);
   }

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle destroy notify

Bool_t TGuiBldDragManager::HandleDestroyNotify(Event_t *event)
{
   if (fPimpl->fPlane && (fPimpl->fPlane->GetId() == event->fWindow)) {
      fPimpl->fPlane = 0;
   }

   return kFALSE;
}


////////////////////////////////////////////////////////////////////////////////
/// not used yet.

Bool_t TGuiBldDragManager::HandleSelection(Event_t *)
{
   if (fStop) {
      return kFALSE;
   }

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// not used yet.

Bool_t TGuiBldDragManager::HandleSelectionRequest(Event_t *)
{
   if (fStop) {
      return kFALSE;
   }

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Find parent frame which can be dragged

TGFrame *TGuiBldDragManager::GetMovableParent(TGWindow *p)
{
   if (fStop) {
      return 0;
   }

   TGFrame *ret = (TGFrame*)p;
   TGWindow *parent = (TGWindow*)ret->GetParent();

   while (parent && (parent != fClient->GetDefaultRoot())) {
      if (!IsFixedLayout(parent) && !IsEditDisabled(parent)) {
         return ret;
      }
      ret = (TGFrame*)parent;
      parent = (TGWindow*)ret->GetParent();
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Find parent frame which can be resized

TGWindow *TGuiBldDragManager::GetResizableParent(TGWindow *p)
{
   if (fStop) {
      return 0;
   }

   TGWindow *parent = p;

   while (parent && (parent != fClient->GetDefaultRoot())) {
      if (!IsFixedSize(parent) &&
          !IsFixedLayout((TGWindow*)parent->GetParent()) &&
          !IsEditDisabled((TGWindow*)parent->GetParent())) {
         return parent;
      }
      parent = (TGWindow*)parent->GetParent();
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Start dragging.

Bool_t TGuiBldDragManager::StartDrag(TGFrame *src, Int_t x, Int_t y)
{
   if (fStop || fDragging) {
      return kFALSE;
   }

   TGFrame *mov = src;

   // special case when frame was grabbed via spacebar pressing
   if (fPimpl->fSpacePressedFrame) {
      if (fDragType == kDragNone) {
         fDragType = kDragMove;
         mov = fPimpl->fSpacePressedFrame;
      } else {
         fPimpl->fSpacePressedFrame = 0;
      }
   }

   TGWindow *parent = (TGWindow*)(mov ? mov->GetParent() : 0);

   // do not remove frame from fixed layout or non-editable parent
   // try to drag "draggable parent"
   if (parent && (IsFixedLayout(parent) || IsEditDisabled(parent))) {
      mov = GetMovableParent(parent);
      if (!mov) {
         return kFALSE;
      }
   }

   SetEditable(kTRUE);  // grab server

   fPimpl->fX = x;
   fPimpl->fY = y;
   fSelectionIsOn = kFALSE;

   fPimpl->fRepeatTimer->Reset();
   gSystem->AddTimer(fPimpl->fRepeatTimer);

   fMoveWaiting = kFALSE;
   fDragging = kTRUE;
   if (src) gVirtualX->SetCursor(src->GetId(), gVirtualX->CreateCursor(kMove));

   switch (fDragType) {
      case kDragCopy:
         HandleCopy();
         HandlePaste();
         GrabFrame(fPimpl->fGrab);
         break;
      case kDragMove:
         fPimpl->fGrab = mov;
         GrabFrame(fPimpl->fGrab);
         break;
      default:
         //fPimpl->fGrab = 0;
         break;
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// End dragging.

Bool_t TGuiBldDragManager::EndDrag()
{
   TGFrame *frame = 0;
   Bool_t ret = kFALSE;

   if (fStop) {
      return kFALSE;
   }

   fMoveWaiting = kFALSE;  // for sanity check

   if (fPimpl->fGrab && (fDragType >= kDragMove) && (fDragType <= kDragLink)) {

      ret = Drop();

   } else if (fBuilder && fBuilder->IsExecutable() &&
              (fDragType == kDragLasso) && !fSelectionIsOn) {

      frame = (TGFrame*)fBuilder->ExecuteAction();
      PlaceFrame(frame, fBuilder->GetAction()->fHints);
      SetLassoDrawn(kFALSE);
      ret = kTRUE;
      //return ret;
   } else if ((fDragType == kDragLasso) && fSelectionIsOn) {

      HandleReturn(kFALSE);
      ret = kTRUE;
   }

   if (!fLassoDrawn) {
      DoRedraw();
   }

   Reset1();
   fPimpl->fSpacePressedFrame = 0;

   if (fBuilder) {
      fBuilder->SetAction(0);
   }

   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// Do cancel action.

Bool_t TGuiBldDragManager::Cancel(Bool_t /*delSrc*/)
{
   if (fStop) {
      return kFALSE;
   }

   fTarget = 0;
   EndDrag();
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Drop grabbed frame

Bool_t TGuiBldDragManager::Drop()
{
   if (fStop || !fDragging || !fPimpl->fGrab ||
       !((fDragType >= kDragMove) && (fDragType <= kDragLink))) {
      return kFALSE;
   }

   fDropStatus = kFALSE;
   TGFrame *frame = 0;
   TGFrame *parent = 0;
   Int_t x, y;
   Window_t c;

   switch (fDragType) {
      case kDragCopy:
      case kDragMove:
         frame = (TGFrame*)fPimpl->fGrab;
         break;
      default:
         break;
   }

   TGWindow *w = fClient->GetWindowById(fTargetId);

   if (fTarget && fPimpl->fGrab && (w == fTarget) &&  w &&
       (w != fClient->GetDefaultRoot())) {
      parent = fTarget;

      gVirtualX->TranslateCoordinates(fClient->GetDefaultRoot()->GetId(),
                                      fTarget->GetId(),
                                      fPimpl->fGrab->GetX(),
                                      fPimpl->fGrab->GetY(), x, y, c);
      fTarget->HandleDragLeave(fPimpl->fGrab);
   } else {
      parent = (TGFrame*)fPimpl->fGrabParent;
      x = fPimpl->fGrabX;
      y = fPimpl->fGrabY;
   }

   //reject move if layout is on
   if (parent && !parent->IsLayoutBroken() && (parent == fPimpl->fGrabParent) ) {
      fDropStatus = 0;
   } else if (parent && frame && (parent != fClient->GetDefaultRoot()) ) {
      ToGrid(x, y);
      fDropStatus = parent->HandleDragDrop(frame, x, y, fPimpl->fGrabLayout);

      // drop was rejected
      if (!fDropStatus) {
         if (fDragType == kDragMove) {  // return dragged frame to initial position
            parent = (TGFrame*)fPimpl->fGrabParent;
            x = fPimpl->fGrabX;
            y = fPimpl->fGrabY;
            frame = fPimpl->fGrab;

            if (parent && frame && (parent != fClient->GetDefaultRoot())) {
               fDropStatus = parent->HandleDragDrop(frame, x, y, fPimpl->fGrabLayout);
            }
         } else { // (fDragType == kDragCopy) - delete it
            DeleteFrame(frame);
         }
      }
   }

   if (fDropStatus) {
      //do not break layout of the new parent if layout there is enabled
      if (parent && !parent->IsLayoutBroken()) {
         parent->Layout();
      }

      if (fBuilder) {
         TString str = frame->ClassName();
         str += "::";
         str += frame->GetName();
         str += " dropped into ";
         str += parent->ClassName();
         str += "::";
         str += parent->GetName();
         str += " at position  ";
         str += TString::Format("(%d , %d)", x, y);
         fBuilder->UpdateStatusBar(str.Data());
      }
      fTarget = 0;
      fTargetId = 0;

      if (parent && (parent == fPimpl->fGrabParent) && fPimpl->fGrabListPosition &&
          frame && parent->InheritsFrom(TGCompositeFrame::Class())) {

         TList *li = ((TGCompositeFrame*)parent)->GetList();
         li->Remove(frame->GetFrameElement());
         li->AddAfter(fPimpl->fGrabListPosition, frame->GetFrameElement());
      }
   } else { // grab frame cannot be dropped
//      if (fDragType == kDragCopy) { // dosn't work (point is not reached ???)
//         HandleDelete(kFALSE);
//      }

      if (fPimpl->fGrab && fPimpl->fGrabParent) {
         fPimpl->fGrab->ReparentWindow(fPimpl->fGrabParent, fPimpl->fGrabX, fPimpl->fGrabY);
         ((TGCompositeFrame*)fPimpl->fGrabParent)->AddFrame(fPimpl->fGrab);
      }
   }

   fPimpl->fGrabParent = 0;
   fPimpl->fGrabX = 0;
   fPimpl->fGrabY = 0;
   fPimpl->fGrabListPosition = 0;

   return fDropStatus;
}

////////////////////////////////////////////////////////////////////////////////
/// Waits for either the mouse move from the given initial ButtonPress location
/// or for the mouse button to be released. If mouse moves away from the initial
/// ButtonPress location before the mouse button is released "IsMoveWaiting"
/// returns kTRUE. If the mouse button released before the mose moved from the
/// initial ButtonPress location, "IsMoveWaiting" returns kFALSE.

Bool_t TGuiBldDragManager::IsMoveWaiting() const
{
   return fMoveWaiting;
}

////////////////////////////////////////////////////////////////////////////////
/// Layout and Resize frame.
/// If global is kFALSE - compact selected frame
/// If global is kFALSE - compact main frame of selected frame

void TGuiBldDragManager::Compact(Bool_t global)
{
   TGCompositeFrame *comp = 0;
   TGFrameElement *fe;

   if (fStop || !fClient || !fClient->IsEditable() || !fPimpl->fGrab) {
      return;
   }

   TGWindow *parent = (TGWindow*)fPimpl->fGrab->GetParent();

   if (global) {
      if (!fBuilder) {
         comp = (TGCompositeFrame*)fClient->GetRoot()->GetMainFrame();
      } else {
         comp = fBuilder->FindEditableMdiFrame(fClient->GetRoot());
         if (!comp) {
            comp = (TGCompositeFrame*)fClient->GetRoot()->GetMainFrame();
         }
      }
   } else {
      if (fPimpl->fGrab &&
          fPimpl->fGrab->InheritsFrom(TGCompositeFrame::Class())) {
         comp = (TGCompositeFrame*)fPimpl->fGrab;
      } else {
         comp = (TGCompositeFrame*)parent;
      }
   }

   if (!comp || IsFixedLayout(comp)  || IsFixedLayout(parent) ||
       IsFixedSize(comp) || IsFixedH(comp) || IsFixedW(comp)) return;

   comp->SetLayoutBroken(kFALSE);

   TIter next(comp->GetList());

   TGFrame *root = (TGFrame *)fClient->GetRoot();
   root->SetEditable(kFALSE);

   TGDimension d;

   if (global) {
      while ((fe = (TGFrameElement*)next())) {
         if (IsFixedLayout(fe->fFrame) || IsFixedSize(fe->fFrame) ||
             IsFixedH(fe->fFrame) || IsFixedW(fe->fFrame)) continue;

         fe->fFrame->SetLayoutBroken(kFALSE);
         d = fe->fFrame->GetDefaultSize();

         // avoid "to point" resizing
         if ((d.fWidth > 10) && (d.fHeight > 10)) {
            fe->fFrame->Resize();
         } else if (d.fWidth > 10) {
            fe->fFrame->Resize(d.fWidth, 10);
         } else if (d.fHeight > 10) {
            fe->fFrame->Resize(10, d.fHeight);
         } else {
            fe->fFrame->Resize(10, 10);
         }
         fClient->NeedRedraw(fe->fFrame);
      }
      if (!IsFixedLayout(root)) {
         root->SetLayoutBroken(kFALSE);
      }
      fPimpl->fCompacted = kTRUE;
   }

   if (!IsFixedLayout(comp)) {
      comp->SetLayoutBroken(kFALSE);
      d = comp->GetDefaultSize();

      // avoid "to point" resizing
      if ((d.fWidth > 10) && (d.fHeight > 10)) {
         comp->Resize();
      } else if (d.fWidth > 10) {
         comp->Resize(d.fWidth, 10);
      } else if (d.fHeight > 10) {
         comp->Resize(10, d.fHeight);
      } else {
         comp->Resize(10, 10);
      }
      layoutFrame(comp);
   }

   if (comp->GetParent()->InheritsFrom(TGMdiDecorFrame::Class())) {
      TGMdiDecorFrame *decor = (TGMdiDecorFrame *)comp->GetParent();
      Int_t b = 2 * decor->GetBorderWidth();
      decor->MoveResize(decor->GetX(), decor->GetY(), comp->GetDefaultWidth() + b,
                        comp->GetDefaultHeight() + b + decor->GetTitleBar()->GetDefaultHeight());
   }

   root->SetEditable(kTRUE);

   fClient->NeedRedraw(comp);
   SelectFrame(comp);
   DoRedraw();
}

////////////////////////////////////////////////////////////////////////////////
/// Grab server.

void TGuiBldDragManager::SetEditable(Bool_t on)
{
   static Bool_t gon = kFALSE;
   static const TGWindow *gw = 0;

   if ((gon == on) && (fClient->GetRoot() == gw)) {
      return;
   }

   gon = on;  gw = fClient->GetRoot();

   if (on) {
      fStop = kFALSE;

      if (fPimpl->fRepeatTimer) {
         fPimpl->fRepeatTimer->Reset();
      } else {
         fPimpl->fRepeatTimer = new TGuiBldDragManagerRepeatTimer(this, 100);
      }
      gSystem->AddTimer(fPimpl->fRepeatTimer);
      ((TGFrame*)fClient->GetRoot())->AddInput(kKeyPressMask | kButtonPressMask);

      Snap2Grid();
   } else {
      HideGrabRectangles();

      if (fPimpl->fRepeatTimer) {
         fPimpl->fRepeatTimer->Remove();
      }

      fSelected = fPimpl->fGrab = 0;

      delete fPimpl->fGrid;
      fPimpl->fGrid = 0;

      fPimpl->ResetParams();

      TGWindow *root = (TGWindow*)fClient->GetRoot();
      if (root) {
         fClient->SetRoot(0);
      }

      if (!gSystem->AccessPathName(fPasteFileName.Data())) {
         gSystem->Unlink(fPasteFileName.Data());
      }

      if (!gSystem->AccessPathName(fTmpBuildFile.Data())) {
         gSystem->Unlink(fTmpBuildFile.Data());
      }

      if (fBuilder) {
         fBuilder->Update();
      }
      //CloseMenus();

      fStop = kTRUE;
   }

   if (on && fClient->IsEditable()) {
      gVirtualX->SetCursor(fClient->GetRoot()->GetId(),
                           gVirtualX->CreateCursor(kPointer));
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return grid coordinates which are close to given

void TGuiBldDragManager::ToGrid(Int_t &x, Int_t &y)
{
   UInt_t step = GetGridStep();
   x = x - x%step;
   y = y - y%step;
}

////////////////////////////////////////////////////////////////////////////////
/// Main handler of actions

void TGuiBldDragManager::HandleAction(Int_t act)
{
   fPimpl->fLastPopupAction = act;

   switch ((EActionType)act) {
      case kPropertyAct:
         CreatePropertyEditor();
         break;
      case kEditableAct:
         if (fPimpl->fSaveGrab) fPimpl->fSaveGrab->SetEditable(kTRUE);
         if (fBuilder) {
            fBuilder->HandleMenu(kGUIBLD_FILE_START);
         }
         break;
      case kCutAct:
         HandleCut();
         break;
      case kCopyAct:
         HandleCopy();
         break;
      case kPasteAct:
         HandlePaste();
         break;
      case kCropAct:
         HandleDelete(kTRUE);
         break;
      case kCompactAct:
         Compact(kFALSE);
         break;
      case kCompactGlobalAct:
         Compact(kTRUE);
         break;
      case kDropAct:
         HandleReturn(kTRUE);
         break;
      case kLayUpAct:
         HandleLayoutOrder(kFALSE);
         break;
      case kLayDownAct:
         HandleLayoutOrder(kTRUE);
         break;
      case kCloneAct:
         CloneEditable();
         break;
      case kGrabAct:
         HandleReturn(kFALSE);
         break;
      case kDeleteAct:
         HandleDelete(kFALSE);
         break;
      case kLeftAct:
         HandleAlignment(kKey_Left);
         break;
      case kRightAct:
         HandleAlignment(kKey_Right);
         break;
      case kUpAct:
         HandleAlignment(kKey_Up);
         break;
      case kDownAct:
         HandleAlignment(kKey_Down);
         break;
      case kEndEditAct:
         if (fBuilder) {
            fBuilder->HandleMenu(kGUIBLD_FILE_STOP);
         }
         SetEditable(kFALSE);
         break;
      case kReplaceAct:
         HandleReplace();
         break;
      case kGridAct:
         HandleGrid();
         break;
      case kBreakLayoutAct:
         BreakLayout();
         break;
      case kSwitchLayoutAct:
      case kLayoutVAct:
      case kLayoutHAct:
         SwitchLayout();
         break;
      case kNewAct:
         if (fBuilder) {
            fBuilder->NewProject();
         } else {
            TGMainFrame *main = new TGMainFrame(fClient->GetDefaultRoot(), 300, 300);
            main->MapRaised();
            main->SetEditable(kTRUE);
         }
         break;
      case kOpenAct:
         if (fBuilder) {
            fBuilder->OpenProject();
         } else {
            TGMainFrame *main = new TGMainFrame(fClient->GetDefaultRoot(), 300, 300);
            main->MapRaised();
            main->SetEditable(kTRUE);
         }
         break;
      case kSaveAct:
         if (fBuilder) {
            if (fBuilder->FindEditableMdiFrame(fClient->GetRoot()) ||
                (!fClient->IsEditable() && fBuilder->GetMdiMain()->GetCurrent())) {
               fBuilder->SaveProject();
            } else {
               Save();
            }
         } else {
            Save();
         }
         break;
      case kSaveFrameAct:
         SaveFrame();
         break;
      default:
         break;
   }

   fPimpl->fPlacePopup = kFALSE;

   if (fBuilder) {
      fBuilder->SetAction(0);
      //fBuilder->Update();
   }

   if (fPimpl->fSaveGrab) {
      fClient->NeedRedraw(fPimpl->fSaveGrab, kTRUE);
   }

   DoRedraw();
}

////////////////////////////////////////////////////////////////////////////////
///  kTRUE - if it's possible to switch disable/enable layout

Bool_t TGuiBldDragManager::CanChangeLayout(TGWindow *w) const
{
   return (!(w->GetEditDisabled() & kEditDisable) &&
           !IsFixedLayout(w) && w->InheritsFrom(TGCompositeFrame::Class()));
}

////////////////////////////////////////////////////////////////////////////////
/// kTRUE - if it's possible to change layout order in the parent's layout of window w

Bool_t TGuiBldDragManager::CanChangeLayoutOrder(TGWindow *w) const
{
   return (w->GetParent()->InheritsFrom(TGCompositeFrame::Class()) &&
           !((TGCompositeFrame*)w->GetParent())->IsLayoutBroken() &&
           !IsFixedLayout((TGWindow*)w->GetParent()));
}

////////////////////////////////////////////////////////////////////////////////
/// kTRUE is frame could be compacted/"layouted"

Bool_t TGuiBldDragManager::CanCompact(TGWindow *w) const
{
   return CanChangeLayout(w);
/*
   return (!IsFixedLayout(w) &&
           w->InheritsFrom(TGCompositeFrame::Class()) &&
           ((TGCompositeFrame*)w)->IsLayoutBroken() &&
           !IsEditDisabled((TGWindow*)w->GetParent()) &&
           !IsFixedLayout((TGWindow*)w->GetParent()));
*/
}

////////////////////////////////////////////////////////////////////////////////
/// Create widget property editor (it could be located outside of guibuilder)

void TGuiBldDragManager::CreatePropertyEditor()
{
//   if (!fPimpl->fClickFrame) return;

   TGWindow *root = (TGWindow*)fClient->GetRoot();
   root->SetEditable(kFALSE);

   fBuilder = (TRootGuiBuilder*)TRootGuiBuilder::Instance();

   fBuilder->Move(fPimpl->fX0, fPimpl->fY0);
   fBuilder->SetWMPosition(fPimpl->fX0, fPimpl->fY0);
   SetPropertyEditor(fBuilder->GetEditor());

   root->SetEditable(kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Helper method

void TGuiBldDragManager::SetPropertyEditor(TGuiBldEditor *e)
{
   fEditor = e;

   if (!fEditor) {
      return;
   }

   ChangeSelected(fPimpl->fClickFrame);
   fEditor->Connect("UpdateSelected(TGFrame*)", "TGuiBldDragManager", this,
                    "HandleUpdateSelected(TGFrame*)");
}

////////////////////////////////////////////////////////////////////////////////
/// Change layout order

void TGuiBldDragManager::HandleLayoutOrder(Bool_t forward)
{
   if (fStop || !fPimpl->fGrab || !fPimpl->fGrab->GetFrameElement() ||
       !CanChangeLayoutOrder(fPimpl->fGrab)) {
      return;
   }

   TGCompositeFrame *comp = (TGCompositeFrame*)fPimpl->fGrab->GetParent();
   TList *li = comp->GetList();
   TGFrameElement *fe = fPimpl->fGrab->GetFrameElement();

   if (!fe) { // sanity check
      return;
   }

   TGFrame *frame;
   TGFrameElement *el;

   if (forward) {
      el = (TGFrameElement *)li->After(fe);
      if (!el) return;
      frame = el->fFrame;

      el->fFrame = fPimpl->fGrab;
      fPimpl->fGrab->SetFrameElement(el);
      fe->fFrame = frame;
      frame->SetFrameElement(fe);
   } else {
      el = (TGFrameElement *)li->Before(fe);

      if (!el) {
         return;
      }
      frame = el->fFrame;

      el->fFrame = fPimpl->fGrab;
      fPimpl->fGrab->SetFrameElement(el);
      fe->fFrame = frame;
      frame->SetFrameElement(fe);
   }

   Bool_t sav = comp->IsLayoutBroken();
   comp->SetLayoutBroken(kFALSE);
   TGWindow *root = (TGWindow *)fClient->GetRoot();
   root->SetEditable(kFALSE);
   comp->Layout();
   DoRedraw();
   root->SetEditable(kTRUE);

   if (sav) {
      comp->SetLayoutBroken(kTRUE);
   }
   SelectFrame(el->fFrame);
}

////////////////////////////////////////////////////////////////////////////////
/// Switch on/of grid drawn.

void TGuiBldDragManager::HandleGrid()
{
   if (fStop) {
      return;
   }

   TGWindow *root = (TGWindow*)fClient->GetRoot();

   if (!root || (root == fClient->GetDefaultRoot())) {
      return;
   }

   if (fPimpl->fGrid->fgStep > 1) {
      fPimpl->fGrid->SetStep(1);
      if (fBuilder) {
         fBuilder->UpdateStatusBar("Grid switched OFF");
      }
   } else {
      fPimpl->fGrid->SetStep(gGridStep);

      if (fBuilder) {
         fBuilder->UpdateStatusBar("Grid switched ON");
      }

      if (root->InheritsFrom(TGCompositeFrame::Class())) {
         TGCompositeFrame *comp = (TGCompositeFrame*)root;
         TIter next(comp->GetList());
         TGFrameElement *fe;
         Int_t x, y, w, h;

         while ((fe = (TGFrameElement*)next())) {
            x = fe->fFrame->GetX();
            y = fe->fFrame->GetY();
            w = fe->fFrame->GetWidth();
            h = fe->fFrame->GetHeight();
            ToGrid(x, y);
            ToGrid(w, h);
            fe->fFrame->MoveResize(x, y, w, h);
         }
      }
   }

   Snap2Grid();
   DrawGrabRectangles();
}

////////////////////////////////////////////////////////////////////////////////
/// Helper to find a frame which can be layouted

TGCompositeFrame *TGuiBldDragManager::FindLayoutFrame(TGFrame *f)
{
   if (fStop || !f) {
      return 0;
   }

   const TGWindow *parent = f->GetParent();
   TGCompositeFrame *ret = 0;

   while (parent && (parent != fClient->GetDefaultRoot())) {
      ret = (TGCompositeFrame*)parent;
      if (parent->InheritsFrom(TGMdiFrame::Class())) return ret;
      parent = parent->GetParent();
   }
   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// When selected frame was changed by guibuilder editor -> update its appearence

void TGuiBldDragManager::HandleUpdateSelected(TGFrame *f)
{
   if (fStop || !f) {
      return;
   }

   TGCompositeFrame *parent = 0;
   if (f->GetParent() &&
       f->GetParent()->InheritsFrom(TGCompositeFrame::Class())) {
      parent = (TGCompositeFrame*)f->GetParent();
   }

   if (!parent || !CanChangeLayout(parent)) {
      return;
   }

   Bool_t sav = parent->IsLayoutBroken();
   parent->SetLayoutBroken(kFALSE);

   if ((parent->GetWidth() < parent->GetDefaultWidth()) ||
        (parent->GetHeight() < parent->GetDefaultHeight())) {
      parent->Resize(parent->GetDefaultSize());
   } else {
      parent->Layout();
      if (f->InheritsFrom(TGCompositeFrame::Class())) {
         layoutFrame(f);
      }
   }
   fClient->NeedRedraw(parent, kTRUE);
   fClient->NeedRedraw(f);

   if (sav) parent->SetLayoutBroken(kTRUE);

   SelectFrame(f);
}

////////////////////////////////////////////////////////////////////////////////
/// Hide/Unmap grab rectangles.

void TGuiBldDragManager::HideGrabRectangles()
{
   static Bool_t first = kFALSE;

   if (fPimpl->fGrabRectHidden) {
      return;
   }
   // skip very first event
   if (!first) {
      first = kTRUE;
      return;
   }
   int i = 0;
   for (i = 0; i < 8; i++) fPimpl->fGrabRect[i]->UnmapWindow();
   for (i = 0; i < 4; i++) fPimpl->fAroundFrame[i]->UnmapWindow();
   fPimpl->fGrabRectHidden = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Delete widget property editor.

void TGuiBldDragManager::DeletePropertyEditor()
{
   if (fStop || !fEditor) {
      return;
   }

   TQObject::Disconnect(fEditor);

   delete fEditor;
   fEditor = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the X coordinate where drag started

Int_t TGuiBldDragManager::GetStrartDragX() const
{
   return fPimpl->fX0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the Y coordinate where drag started

Int_t TGuiBldDragManager::GetStrartDragY() const
{
   return fPimpl->fY0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the current X coordinate of the dragged frame

Int_t TGuiBldDragManager::GetEndDragX() const
{
   return fPimpl->fY;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the current Y coordinate of the dragged frame

Int_t TGuiBldDragManager::GetEndDragY() const
{
   return fPimpl->fY;
}

////////////////////////////////////////////////////////////////////////////////
/// Disable/Enable layout for selected/grabbed composite frame.

void TGuiBldDragManager::BreakLayout()
{
   if (fStop) {
      return;
   }

   TGFrame *frame = fSelected;

   if (!frame) {
      return;
   }

   TString str = frame->ClassName();
   str += "::";
   str += frame->GetName();

   if (IsFixedLayout(frame)) {
      if (fBuilder) {
         str += " layout cannot be broken";
         fBuilder->UpdateStatusBar(str.Data());
      }
      return;
   }

   frame->SetLayoutBroken(!frame->IsLayoutBroken());
   DrawGrabRectangles();

   if (fBuilder) {
      str += (frame->IsLayoutBroken() ? " Disable Layout" : " Enable Layout");
      fBuilder->UpdateStatusBar(str.Data());
   }
   if (fPimpl->fGrab && (fPimpl->fGrab->IsA() == TGCanvas::Class())) {
      fPimpl->fGrab->Layout();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Switch Horizontal/Vertical layout of selected/grabbed composite frame

void TGuiBldDragManager::SwitchLayout()
{
   if (fStop || !fPimpl->fGrab) {
      return;
   }

   TGCompositeFrame *comp = (TGCompositeFrame*)fSelected;

   comp->SetLayoutBroken(kFALSE);

   UInt_t opt = comp->GetOptions();
   TGLayoutManager *m = comp->GetLayoutManager();

   if (!m) {
      return;
   }

   if (m->InheritsFrom(TGHorizontalLayout::Class())) {
      opt &= ~kHorizontalFrame;
      opt |= kVerticalFrame;

      if (fBuilder) {
         TString str = comp->ClassName();
         str += "::";
         str += comp->GetName();
         str += " Vertical Layout ON";
         fBuilder->UpdateStatusBar(str.Data());
      }
   } else if (m->InheritsFrom(TGVerticalLayout::Class())) {
      opt &= ~kVerticalFrame;
      opt |= kHorizontalFrame;

      if (fBuilder) {
         TString str = comp->ClassName();
         str += "::";
         str += comp->GetName();
         str += " Horizontal Layout ON";
         fBuilder->UpdateStatusBar(str.Data());
      }
   }

   comp->ChangeOptions(opt);
   if (!IsFixedSize(comp)) {
      comp->Resize();
   }

   if (fPimpl->fGrab && (fPimpl->fGrab->IsA() == TGCanvas::Class())) {
      fPimpl->fGrab->Layout();
   }

   fClient->NeedRedraw(comp);
   SelectFrame(comp);
}

////////////////////////////////////////////////////////////////////////////////
/// Return the current grabbed/selected frame.

TGFrame *TGuiBldDragManager::GetSelected() const
{
   return fSelected;
}

////////////////////////////////////////////////////////////////////////////////
/// Helper to close all menus

void TGuiBldDragManager::CloseMenus()
{
   void *ud;

   if (fFrameMenu) {
      fFrameMenu->EndMenu(ud);
   }
   if (fLassoMenu) {
      fLassoMenu->EndMenu(ud);
   }
   //UnmapAllPopups();
}

////////////////////////////////////////////////////////////////////////////////
/// Return the parent frame which can be editted.

TGFrame *TGuiBldDragManager::GetEditableParent(TGFrame *fr)
{
   if (!fr || (fr == fClient->GetDefaultRoot())) {
      return 0;
   }

   TGWindow *parent = (TGWindow*)fr->GetParent();

   while (parent && (parent != fClient->GetDefaultRoot())) {
      if (!IsEditDisabled(parent) && !IsGrabDisabled(parent)) {
         return (TGFrame*)parent;
      }
      parent = (TGWindow*)parent->GetParent();
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return a name of icon

static TString FindMenuIconName(TString &in)
{
   Int_t p1 = in.Index("*icon=", 1);
   if (p1 == kNPOS) return "";
   p1 += 6;
   Int_t p2 = in.Index("*", p1);

   if (p2 == kNPOS) return "";

   return in(p1, p2-p1);
}

////////////////////////////////////////////////////////////////////////////////
/// Helper

static Bool_t containBaseClass(const char *somestring, TClass *cl)
{
   TString str = somestring;

   if (str.Contains(cl->GetName())) {
      return kTRUE;
   }

   TIter nextBaseClass(cl->GetListOfBases());
   TBaseClass *bc;

   while ((bc = (TBaseClass*)nextBaseClass())) {
      if (!bc->GetClassPointer()) {
         continue;
      }
      if (containBaseClass(somestring, bc->GetClassPointer())) {
         return kTRUE;
      }
   }

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Add DIALOG entries to the selected frame popup menu

void TGuiBldDragManager::AddDialogMethods(TGPopupMenu *menu, TObject *object)
{
   if (!menu || !object) {
      return;
   }

   TMethod *method;
   TIter next(fListOfDialogs);
   TString str;
   TString pname;
   const TGPicture *pic;
   TClass *cl = object->IsA();
   TString ename;

   while ((method = (TMethod*) next())) {
      ename = method->GetName();
      ename += "...";
      if (menu->GetEntry(ename.Data())) {
         continue;
      }
      if (!containBaseClass(method->GetSignature(), cl)) {
         continue;
      }

      str = method->GetCommentString();
      pname = FindMenuIconName(str);
      pic = fClient->GetPicture(pname.Data());
      menu->AddEntry(ename.Data(), kMethodMenuAct, method, pic);
   }
   menu->AddSeparator();
}

////////////////////////////////////////////////////////////////////////////////
/// Add entries with class //*MENU* methods

void TGuiBldDragManager::AddClassMenuMethods(TGPopupMenu *menu, TObject *object)
{
   if (!menu || !object) {
      return;
   }

   TList *menuItemList;
   TClassMenuItem *menuItem;
   TString str;
   TString pname;
   const TGPicture *pic;
   TMethod  *method;
   TClass   *classPtr = 0;
   TList    *methodList;
   EMenuItemKind menuKind;
   TDataMember *m;

   AddDialogMethods(menu, object);

   menuItemList = object->IsA()->GetMenuList();
   TIter nextItem(menuItemList);

   fPimpl->fMenuObject = (TGFrame*)object;
   nextItem.Reset();

   while ((menuItem = (TClassMenuItem*) nextItem())) {
      switch (menuItem->GetType()) {
         case TClassMenuItem::kPopupStandardList:
            {
               // Standard list of class methods. Rebuild from scratch.
               // Get linked list of objects menu items (i.e. member functions
               // with the token *MENU in their comment fields.
               methodList = new TList;
               object->IsA()->GetMenuItems(methodList);

               TIter next(methodList);

               while ((method = (TMethod*) next())) {
                  if (classPtr != method->GetClass()) {
//                     menu->AddSeparator();
                     classPtr = method->GetClass();
                  }

                  menuKind = method->IsMenuItem();

                  switch (menuKind) {
                     case kMenuDialog:
                     {
                        str = method->GetCommentString();
                        pname = FindMenuIconName(str);
                        pic = fClient->GetPicture(pname.Data());
                        menu->AddEntry(method->GetName(), kMethodMenuAct, method, pic);
                        break;
                     }

                     case kMenuSubMenu:
                        if ((m = method->FindDataMember())) {
                           if (m->GetterMethod()) {
                              TGPopupMenu *r = TRootGuiBuilder::CreatePopup();
                              menu->AddPopup(method->GetName(), r);
                              fPimpl->fFrameMenuTrash->Add(r);
                              TIter nxt(m->GetOptions());
                              TOptionListItem *it;

                              while ((it = (TOptionListItem*) nxt())) {
                                 const char  *name  = it->fOptName;
                                 Long_t val   = it->fValue;

                                 TToggle *t = new TToggle;
                                 t->SetToggledObject(object, method);
                                 t->SetOnValue(val);
                                 fPimpl->fFrameMenuTrash->Add(t);

                                 //r->AddSeparator();
                                 r->AddEntry(name, kToggleMenuAct, t);
                                 if (t->GetState()) r->CheckEntryByData(t);
                              }
                           } else {
                              menu->AddEntry(method->GetName(), kMethodMenuAct, method);
                           }
                        }
                        break;

                     case kMenuToggle:
                        {
                           TToggle *t = new TToggle;
                           t->SetToggledObject(object, method);
                           t->SetOnValue(1);
                           fPimpl->fFrameMenuTrash->Add(t);
                           menu->AddEntry(method->GetName(), kToggleMenuAct, t);
                           if (t->GetState()) menu->CheckEntryByData(t);
                        }
                        break;

                     default:
                        break;
                  }
               }
               delete methodList;
            }
            break;
         case TClassMenuItem::kPopupUserFunction:
            {
               if (menuItem->IsToggle()) {
                  TMethod* method2 =
                        object->IsA()->GetMethodWithPrototype(menuItem->GetFunctionName(),
                                                              menuItem->GetArgs());
                  if (method2) {
                     TToggle *t = new TToggle;
                     t->SetToggledObject(object, method2);
                     t->SetOnValue(1);
                     fPimpl->fFrameMenuTrash->Add(t);

                     menu->AddEntry(method2->GetName(), kToggleMenuAct, t);
                     if (t->GetState()) menu->CheckEntryByData(t);
                  }
               } else {
                  const char* menuItemTitle = menuItem->GetTitle();
                  if (strlen(menuItemTitle)==0) menuItemTitle = menuItem->GetFunctionName();
                  menu->AddEntry(menuItemTitle, kMethodMenuAct, menuItem);
               }
            }

            break;
         default:
            break;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Process a method chosen via frame context menu

void TGuiBldDragManager::DoClassMenu(Int_t id)
{
   if (!fFrameMenu || ((id != kMethodMenuAct) && (id != kToggleMenuAct))) {
      return;
   }

   TGMenuEntry *me = 0;

   if (id == kMethodMenuAct) {
      delete gMenuDialog;
      me = fFrameMenu->GetCurrent();

      if (!me || !fPimpl->fMenuObject) {
         return;
      }
      TMethod *method = (TMethod*)me->GetUserData();
      TString str = method->GetCommentString();

      if (str.Contains("*DIALOG")) {
         TString str2;
         str2.Form("((TGuiBldDragManager*)0x%lx)->%s((%s*)0x%lx)", (ULong_t)this, method->GetName(),
                  fPimpl->fMenuObject->ClassName(), (ULong_t)fPimpl->fMenuObject);
         gCling->Calc((char *)str2.Data());
         //delete fFrameMenu;  // suicide (BB)?
         //fFrameMenu = 0;
         return;
      }
      gMenuDialog = new TGuiBldMenuDialog(fPimpl->fMenuObject, fPimpl->fMenuObject, method);
      gMenuDialog->Popup();

   } else if (id == kToggleMenuAct) {
      me = fFrameMenu->GetCurrent();
      if (!me) {
         return;
      }
      TGPopupMenu *menu = me->GetPopup();
      TToggle *toggle = 0;

      if (menu) {    //process submenu
         toggle = (TToggle*)menu->GetCurrent()->GetUserData();
      } else {    //process check entry
         toggle = (TToggle*)fFrameMenu->GetCurrent()->GetUserData();
      }
      if (toggle) {
         toggle->Toggle();
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Delete dialog and trash

void TGuiBldDragManager::DeleteMenuDialog()
{
   fPimpl->fFrameMenuTrash->Delete();
   gMenuDialog->DeleteWindow();
   gMenuDialog = 0;
   fPimpl->fMenuObject = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Process dialog OK button pressed

void TGuiBldDragManager::DoDialogOK()
{
   gMenuDialog->ApplyMethod();
   DoRedraw();
   DeleteMenuDialog();
   gMenuDialog = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Process dialog Apply button pressed

void TGuiBldDragManager::DoDialogApply()
{
   gMenuDialog->ApplyMethod();
}

////////////////////////////////////////////////////////////////////////////////
/// Process dialog Cancel button pressed

void TGuiBldDragManager::DoDialogCancel()
{
   DeleteMenuDialog();
   gMenuDialog = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Create and  place context menu for selected frame

void TGuiBldDragManager::Menu4Frame(TGFrame *frame, Int_t x, Int_t y)
{
   if (fStop) {
      return;
   }

   fPimpl->fSaveGrab = fPimpl->fGrab;
   fPimpl->fX0 = x;
   fPimpl->fY0 = y;
   fPimpl->fClickFrame = frame;

   Bool_t composite = frame->InheritsFrom(TGCompositeFrame::Class());
   Bool_t compar = frame->GetParent()->InheritsFrom(TGCompositeFrame::Class());

   TGCompositeFrame *cfr = 0;
   TGCompositeFrame *cfrp = 0;
   TGLayoutManager *lm = 0;

   if (composite)  {
      cfr = (TGCompositeFrame *)frame;
      lm = cfr->GetLayoutManager();
   }
   if (compar)  {
      cfrp = (TGCompositeFrame *)frame->GetParent();
   }

   delete fFrameMenu;

   fFrameMenu = TRootGuiBuilder::CreatePopup();
   fFrameMenu->Connect("Activated(Int_t)", "TGuiBldDragManager", this, "DoClassMenu(Int_t)");

   TString title = frame->ClassName();
   title += "::";
   title += frame->GetName();
   fFrameMenu->AddLabel(title.Data());
   fFrameMenu->AddSeparator();

   // special case - menu for editable Mdi frame
   if (fBuilder && (frame == fBuilder->GetMdiMain()->GetCurrent())) {
      if (!gSystem->AccessPathName(fPasteFileName.Data())) {
         fFrameMenu->AddEntry("Paste\tCtrl+V", kPasteAct,
                               0, fClient->GetPicture("bld_paste.png"));
      }
      fFrameMenu->AddEntry("Compact\tCtrl+L", kCompactAct,
                               0, fClient->GetPicture("bld_compact.png"));
      fFrameMenu->AddEntry("Grid On/Off\tCtrl+G", kGridAct,
                              0, fClient->GetPicture("bld_grid.png"));
      fFrameMenu->AddEntry("Save As ...\tCtrl+S", kSaveAct,
                              0, fClient->GetPicture("bld_save.png"));
      fFrameMenu->AddEntry("End Edit\tCtrl+DblClick", kEndEditAct,
                              0, fClient->GetPicture("bld_stop.png"));
      goto out;
   }

   AddClassMenuMethods(fFrameMenu, frame);

   if (!fBuilder) {
      fFrameMenu->AddEntry("Gui Builder", kPropertyAct);
      fFrameMenu->AddSeparator();
   }
/*
   if (!frame->IsEditable() && !InEditable(frame->GetId())) {
      fPimpl->fSaveGrab = frame;
      goto out;
   }
*/
   if (!IsEditDisabled(cfrp)) {
      fFrameMenu->AddSeparator();

      if (composite && !IsFixedLayout(frame) && cfr->GetList()->GetEntries()) {
         fFrameMenu->AddEntry("Drop\tCtrl+Return", kDropAct);
      }

      if (!IsFixedLayout(cfrp)) {
         fFrameMenu->AddEntry("Cut\tCtrl+X", kCutAct,
                               0, fClient->GetPicture("bld_cut.png"));
      }
      //
      fFrameMenu->AddEntry("Copy\tCtrl+C", kCopyAct,
                            0, fClient->GetPicture("bld_copy.png"));

      if (frame->IsEditable() && !IsFixedLayout(frame) &&
          !gSystem->AccessPathName(fPasteFileName.Data())) {
         fFrameMenu->AddEntry("Paste\tCtrl+V", kPasteAct,
                               0, fClient->GetPicture("bld_paste.png"));
      }

      if (!IsFixedLayout(cfrp)) {
         fFrameMenu->AddEntry("Delete\tDel", kDeleteAct,
                              0, fClient->GetPicture("bld_delete.png"));
      }

      if (!IsFixedLayout(cfrp)) {
         fFrameMenu->AddEntry("Crop\tShift+Del", kCropAct,
                               0, fClient->GetPicture("bld_crop.png"));
      }

//      if (!IsFixedLayout(cfrp) && !gSystem->AccessPathName(fPasteFileName.Data())) {
//         fFrameMenu->AddEntry("Replace\tCtrl+R", kReplaceAct,
//                               0, fClient->GetPicture("bld_paste_into.png"));
//      }

      fFrameMenu->AddSeparator();
   } else {
      if (!gSystem->AccessPathName(fPasteFileName.Data()) && !IsFixedLayout(frame)) {
         fFrameMenu->AddEntry("Paste\tCtrl+V", kPasteAct,
                               0, fClient->GetPicture("bld_paste.png"));
      }
      if (frame->GetMainFrame() == frame) {
         fFrameMenu->AddEntry("Clone\tCtrl+A", kCloneAct);
      }
      fFrameMenu->AddSeparator();
   }

   if (CanChangeLayout(frame)) {
      const char *label = (frame->IsLayoutBroken() ? "Allow Layout\tCtrl+B" :
                                                     "Break Layout\tCtrl+B");
      fFrameMenu->AddEntry(label, kBreakLayoutAct,
                            0, fClient->GetPicture("bld_break.png"));
   }

   if (composite && !cfr->GetList()->IsEmpty()) {
      if (CanCompact(frame)) {
         if (!frame->IsEditable()) {
            fFrameMenu->AddEntry("Compact\tCtrl+L", kCompactAct,
                                  0, fClient->GetPicture("bld_compact.png"));
         } else {
            fFrameMenu->AddEntry("Compact\tCtrl+L", kCompactGlobalAct,
                                  0, fClient->GetPicture("bld_compact.png"));
         }
      }

      if (lm && ((lm->IsA() == TGVerticalLayout::Class()) ||
           (lm->IsA() == TGHorizontalLayout::Class())) && !IsFixedLayout(frame)) {

         if (lm->IsA() == TGVerticalLayout::Class()) {
            fFrameMenu->AddEntry("Horizontal\tCtrl+H", kSwitchLayoutAct,
                                 0, fClient->GetPicture("bld_hbox.png"));
         } else if (lm->IsA() == TGHorizontalLayout::Class()) {
            fFrameMenu->AddEntry("Vertical\tCtrl+H", kSwitchLayoutAct,
                                 0, fClient->GetPicture("bld_vbox.png"));
         }
      }
   }

   if (compar && (cfrp->GetList()->GetSize() > 1) && CanChangeLayoutOrder(frame)) {
      if (cfrp->GetList()->First() != frame->GetFrameElement()) {
         fFrameMenu->AddEntry("Lay Up\tUp/Left", kLayUpAct);
      }
      if (cfrp->GetList()->Last() != frame->GetFrameElement()) {
         fFrameMenu->AddEntry("Lay Down\tDown/Right", kLayDownAct);
      }
      fFrameMenu->AddSeparator();
   }

   if (frame->IsEditable()) {
      fFrameMenu->AddEntry("Grid On/Off\tCtrl+G", kGridAct,
                            0, fClient->GetPicture("bld_grid.png"));
   }
   if (composite && !cfr->GetList()->IsEmpty()) {
      fPimpl->fSaveGrab = frame;
      fFrameMenu->AddEntry("Save As ...       ", kSaveFrameAct,
                            0, fClient->GetPicture("bld_save.png"));
   }

out:
   fFrameMenu->Connect("Activated(Int_t)", "TGuiBldDragManager", this, "HandleAction(Int_t)");

   fPimpl->fLastPopupAction = kNoneAct;
   fPimpl->fPlacePopup = kTRUE;

   fFrameMenu->PlaceMenu(x, y, kTRUE, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Create context menu for lasso actions.

void TGuiBldDragManager::Menu4Lasso(Int_t x, Int_t y)
{
   if (fStop || !fLassoDrawn) {
      return;
   }

   DrawLasso();

   delete fLassoMenu;

   fLassoMenu = TRootGuiBuilder::CreatePopup();
   fLassoMenu->AddLabel("Edit actions");
   fLassoMenu->AddSeparator();
   fLassoMenu->AddEntry("Grab\tReturn", kGrabAct);
   fLassoMenu->AddSeparator();
   fLassoMenu->AddEntry("Delete\tDelete", kDeleteAct,
                        0, fClient->GetPicture("bld_delete.png"));
   fLassoMenu->AddEntry("Crop\tShift+Delete", kCropAct,
                        0, fClient->GetPicture("bld_crop.png"));
   fLassoMenu->AddSeparator();
   fLassoMenu->AddEntry("Align Left\tLeft Key", kLeftAct,
                        0, fClient->GetPicture("bld_AlignLeft.png"));
   fLassoMenu->AddEntry("Align Right\tRight Key", kRightAct,
                        0, fClient->GetPicture("bld_AlignRight.png"));
   fLassoMenu->AddEntry("Align Up\tUp Key", kUpAct,
                        0, fClient->GetPicture("bld_AlignTop.png"));
   fLassoMenu->AddEntry("Align Down\tDown Key", kDownAct,
                        0, fClient->GetPicture("bld_AlignBtm.png"));

   fLassoMenu->Connect("Activated(Int_t)", "TGuiBldDragManager", this, "HandleAction(Int_t)");

   fPimpl->fLastPopupAction = kNoneAct;
   fPimpl->fPlacePopup = kTRUE;
   fLassoMenu->PlaceMenu(x, y, kTRUE, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Return kTRUE if paste frame exist.

Bool_t TGuiBldDragManager::IsPasteFrameExist()
{
   return !gSystem->AccessPathName(fPasteFileName.Data());
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to global color dialog. If dialog is not yet created
/// and input parameter is kTRUE - the dialog will be created.

TGColorDialog *TGuiBldDragManager::GetGlobalColorDialog(Bool_t create)
{
   static Int_t retc;
   static Pixel_t color;

   if (!fgGlobalColorDialog && create) {
      fgGlobalColorDialog = new TGColorDialog(gClient->GetDefaultRoot(), 0,
                                              &retc, &color, kFALSE);
      int i = 0;
      for (i = 0; i < 10; i++) {
         fgGlobalColorDialog->GetCustomPalette()->SetColor(i, TColor::Number2Pixel(i));
      }
      for (i = 0; i < 10; i++) {
         fgGlobalColorDialog->GetCustomPalette()->SetColor(10+i, TColor::Number2Pixel(180+i));
      }
   }
   return fgGlobalColorDialog;
}

////////////////////////////////////////////////////////////////////////////////
/// Create global font dialog.

TGFontDialog *TGuiBldDragManager::GetGlobalFontDialog()
{
   static TGFontDialog::FontProp_t prop;

   if (!fgGlobalFontDialog) {
     fgGlobalFontDialog = new TGFontDialog(gClient->GetDefaultRoot(), 0, &prop, "", 0, kFALSE);
   }
   return fgGlobalFontDialog;
}

////////////////////////////////////////////////////////////////////////////////
/// Map dialog and place it relative to selected frame.

void TGuiBldDragManager::MapGlobalDialog(TGMainFrame *dialog, TGFrame *fr)
{
   Int_t x = 0, y = 0;
   Window_t wdummy;
   UInt_t dw = gClient->GetDisplayWidth() - 20;
   UInt_t dh = gClient->GetDisplayHeight() - 50;

   TGFrame *parent = (TGFrame*)fr->GetParent();
   gVirtualX->TranslateCoordinates(parent->GetId(), gClient->GetDefaultRoot()->GetId(),
                                   fr->GetX() + fr->GetWidth(),
                                   fr->GetY() + fr->GetHeight(), x, y, wdummy);

   if (x + dialog->GetWidth() > dw) {
      x = dw - dialog->GetWidth();
   }

   if (y + dialog->GetHeight() > dh) {
      y = dh - dialog->GetHeight();
   }

   dialog->Move(x, y);
   dialog->SetWMPosition(x, y);
   dialog->MapRaised();
}

////////////////////////////////////////////////////////////////////////////////
/// Change background color via context menu.

void TGuiBldDragManager::ChangeBackgroundColor(TGFrame *fr)
{
   TGColorDialog *cd = GetGlobalColorDialog();
   cd->SetCurrentColor(fr->GetBackground());
   cd->Connect("ColorSelected(Pixel_t)", "TGFrame", fr, "ChangeBackground(Pixel_t)");
   MapGlobalDialog(cd, fr);
   fClient->WaitForUnmap(cd);
   TQObject::Disconnect(cd);
}

////////////////////////////////////////////////////////////////////////////////
/// Change background color via context menu for this frame and all subframes.
/// This method is activated via context menu during guibuilding.

void TGuiBldDragManager::ChangeBackgroundColor(TGCompositeFrame *fr)
{
   TGColorDialog *cd = GetGlobalColorDialog();
   cd->SetCurrentColor(fr->GetBackground());
   cd->Connect("ColorSelected(Pixel_t)", "TGCompositeFrame", fr,
               "ChangeSubframesBackground(Pixel_t)");
   MapGlobalDialog(cd, fr);
   fClient->WaitForUnmap(cd);
   TQObject::Disconnect(cd);
}

////////////////////////////////////////////////////////////////////////////////
/// Change text color via color selection dialog. This method is activated
/// via context menu during guibuilding.

void TGuiBldDragManager::ChangeTextColor(TGGroupFrame *fr)
{
   TGGC *gc = fClient->GetResourcePool()->GetGCPool()->FindGC(fr->GetNormGC());

   if (!gc) {
      return;
   }
   ULong_t color = gc->GetForeground();

   TGColorDialog *cd = GetGlobalColorDialog();
   cd->SetCurrentColor(color);
   cd->Connect("ColorSelected(Pixel_t)", "TGGroupFrame", fr, "SetTextColor(Pixel_t)");
   MapGlobalDialog(cd, fr);
   fClient->WaitForUnmap(cd);
   TQObject::Disconnect(cd);
}

////////////////////////////////////////////////////////////////////////////////
/// Change text font via font selection dialog. This method is activated
/// via context menu during guibuilding.

void TGuiBldDragManager::ChangeTextFont(TGGroupFrame *fr)
{
   TGFontDialog *fd = GetGlobalFontDialog();

   TGGC *gc = fClient->GetResourcePool()->GetGCPool()->FindGC(fr->GetNormGC());

   if (!gc) {
      return;
   }

   TGFont *font = fClient->GetResourcePool()->GetFontPool()->FindFont(fr->GetFontStruct());

   if (!font) {
      return;
   }
   fd->SetColor(gc->GetForeground());
   fd->SetFont(font);
   fd->EnableAlign(kFALSE);
   fd->Connect("FontSelected(char*)", "TGGroupFrame", fr, "SetTextFont(char*)");
   fd->Connect("ColorSelected(Pixel_t)", "TGGroupFrame", fr, "SetTextColor(Pixel_t)");

   MapGlobalDialog(fd, fr);
   fClient->WaitForUnmap(fd);
   TQObject::Disconnect(fd);
}

////////////////////////////////////////////////////////////////////////////////
/// Edit properties via font selection dialog. This method is activated
/// via context menu during guibuilding.

void TGuiBldDragManager::ChangeProperties(TGTextButton *fr)
{
   TGFontDialog *fd = GetGlobalFontDialog();

   TGGC *gc = fClient->GetResourcePool()->GetGCPool()->FindGC(fr->GetNormGC());
   if (!gc) {
      return;
   }

   TGFont *font = fClient->GetResourcePool()->GetFontPool()->FindFont(fr->GetFontStruct());

   if (!font) {
      return;
   }
   fd->SetColor(gc->GetForeground());
   fd->SetFont(font);
   fd->SetAlign(fr->GetTextJustify());

   fd->Connect("FontSelected(char*)", "TGTextButton", fr, "SetFont(char*)");
   fd->Connect("ColorSelected(Pixel_t)", "TGTextButton", fr, "SetTextColor(Pixel_t)");
   fd->Connect("AlignSelected(Int_t)", "TGTextButton", fr, "SetTextJustify(Int_t)");

   MapGlobalDialog(fd, fr);
   fClient->WaitForUnmap(fd);
   TQObject::Disconnect(fd);
}

////////////////////////////////////////////////////////////////////////////////
/// Change text color via color selection dialog. This method is activated
/// via context menu during guibuilding.

void TGuiBldDragManager::ChangeTextColor(TGTextButton *fr)
{
   TGGC *gc = gClient->GetResourcePool()->GetGCPool()->FindGC(fr->GetNormGC());

   if (!gc) {
      return;
   }
   ULong_t color = gc->GetForeground();

   TGColorDialog *cd = GetGlobalColorDialog();
   cd->SetCurrentColor(color);
   cd->Connect("ColorSelected(Pixel_t)", "TGTextButton", fr, "SetTextColor(Pixel_t)");

   MapGlobalDialog(cd, fr);
   fClient->WaitForUnmap(cd);
   TQObject::Disconnect(cd);
}

////////////////////////////////////////////////////////////////////////////////
/// Invoke file dialog to assign a new picture.
/// This method is activated via context menu during guibuilding.

void TGuiBldDragManager::ChangePicture(TGPictureButton *fr)
{
   static TGFileInfo fi;
   static TString dir(".");
   static Bool_t overwr = kFALSE;
   TString fname;

   fi.fFileTypes = gImageTypes;
   fi.SetIniDir(dir);
   fi.fOverwrite = overwr;

   TGWindow *root = (TGWindow*)fClient->GetRoot();
   SetEditable(kFALSE);

   new TGFileDialog(fClient->GetDefaultRoot(), fr, kFDOpen, &fi);

   if (!fi.fFilename) {
      root->SetEditable(kTRUE);
      SetEditable(kTRUE);
      return;
   }

   dir    = fi.fIniDir;
   overwr = fi.fOverwrite;
   fname  = fi.fFilename;

   const TGPicture *pic = fClient->GetPicture(fname.Data());

   if (!pic) {
      Int_t retval;
      new TGMsgBox(fClient->GetDefaultRoot(), fr, "Error...",
                   TString::Format("Cannot read image file (%s)", fname.Data()),
                   kMBIconExclamation, kMBRetry | kMBCancel, &retval);

      if (retval == kMBRetry) {
         ChangePicture(fr);
      }
   } else {
      const TGPicture *tmp = fr->GetPicture();
      if (tmp) fClient->FreePicture(tmp);

      fr->SetPicture(pic);

      // not clear how to do at this point
      tmp = fr->GetDisabledPicture();
      if (tmp) fClient->FreePicture(tmp);
   }
   root->SetEditable(kTRUE);
   SetEditable(kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Change background color via context menu

void TGuiBldDragManager::ChangeBackgroundColor(TGCanvas *fr)
{
   TGColorDialog *cd = GetGlobalColorDialog();
   cd->SetCurrentColor(fr->GetBackground());
   cd->Connect("ColorSelected(Pixel_t)", "TGFrame", fr, "ChangeBackground(Pixel_t)");
   cd->Connect("ColorSelected(Pixel_t)", "TGScrollBar", fr->GetHScrollbar(), "ChangeBackground(Pixel_t)");
   cd->Connect("ColorSelected(Pixel_t)", "TGScrollBar", fr->GetVScrollbar(), "ChangeBackground(Pixel_t)");

   MapGlobalDialog(cd, fr);
   fClient->WaitForUnmap(cd);
   TQObject::Disconnect(cd);
}

////////////////////////////////////////////////////////////////////////////////
/// Change background color for list box entries. This method is invoked
/// via context menu during guibuilding.

void TGuiBldDragManager::ChangeBackgroundColor(TGComboBox *fr)
{
   Pixel_t color = TGFrame::GetWhitePixel();

   TGColorDialog *cd = GetGlobalColorDialog();
   cd->SetCurrentColor(color);

   cd->Connect("ColorSelected(Pixel_t)", "TGListBox", fr->GetListBox(),
               "ChangeBackground(Pixel_t)");

   TGLBEntry *se = fr->GetSelectedEntry();

   if (se) {
      cd->Connect("ColorSelected(Pixel_t)", "TGLBEntry", se,
                  "SetBackgroundColor(Pixel_t)");
   }

   TGTextEntry *te = fr->GetTextEntry();

   if (te) {
      cd->Connect("ColorSelected(Pixel_t)", "TGTextEntry", te,
                  "SetBackgroundColor(Pixel_t)");
   }

   MapGlobalDialog(cd, fr);
   fClient->WaitForUnmap(cd);
   TQObject::Disconnect(cd);

   if (se) {
      fClient->NeedRedraw(se, kTRUE); // force redraw
   }

   if (te) {
      fClient->NeedRedraw(te, kTRUE);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Edit properties via font selection dialog. This method is activated
/// via context menu during guibuilding.

void TGuiBldDragManager::ChangeProperties(TGLabel *fr)
{
   TGFontDialog *fd = GetGlobalFontDialog();

   TGGC *gc = fClient->GetResourcePool()->GetGCPool()->FindGC(fr->GetNormGC());

   if (!gc) {
      return;
   }

   TGFont *font = fClient->GetResourcePool()->GetFontPool()->FindFont(fr->GetFontStruct());

   if (!font) {
      return;
   }

   fd->SetColor(gc->GetForeground());
   fd->SetFont(font);
   fd->SetAlign(fr->GetTextJustify());

   fd->Connect("FontSelected(char*)", "TGLabel", fr, "SetTextFont(char*)");
   fd->Connect("ColorSelected(Pixel_t)", "TGLabel", fr, "SetTextColor(Pixel_t)");
   fd->Connect("AlignSelected(Int_t)", "TGLabel", fr, "SetTextJustify(Int_t)");

   MapGlobalDialog(fd, fr);
   fClient->WaitForUnmap(fd);
   TQObject::Disconnect(fd);
}

////////////////////////////////////////////////////////////////////////////////
/// Change text color via color selection dialog. This method is activated
/// via context menu during guibuilding.

void TGuiBldDragManager::ChangeTextColor(TGLabel *fr)
{
   TGGC *gc = gClient->GetResourcePool()->GetGCPool()->FindGC(fr->GetNormGC());

   if (!gc) {
      return;
   }

   ULong_t color = gc->GetForeground();

   TGColorDialog *cd = GetGlobalColorDialog();
   cd->SetCurrentColor(color);
   cd->Connect("ColorSelected(Pixel_t)", "TGLabel", fr, "SetTextColor(Pixel_t)");

   MapGlobalDialog(cd, fr);
   fClient->WaitForUnmap(cd);
   TQObject::Disconnect(cd);
}

////////////////////////////////////////////////////////////////////////////////
/// Set background color for list box entries. This method is invoked
/// via context menu during  guibuilding.

void TGuiBldDragManager::ChangeBackgroundColor(TGListBox *fr)
{
   Pixel_t color = TGFrame::GetWhitePixel();

   TGColorDialog *cd = GetGlobalColorDialog();
   cd->SetCurrentColor(color);
   cd->Connect("ColorSelected(Pixel_t)", "TGListBox", fr, "ChangeBackground(Pixel_t)");

   MapGlobalDialog(cd, fr);
   fClient->WaitForUnmap(cd);
   TQObject::Disconnect(cd);
}

////////////////////////////////////////////////////////////////////////////////
/// Set progress bar color via TGColorDialog.
/// This method is activated via context menu during guibuilding.

void TGuiBldDragManager::ChangeBarColor(TGProgressBar *fr)
{
   ULong_t color = fr->GetBarColor();

   TGColorDialog *cd = GetGlobalColorDialog();

   cd->SetCurrentColor(color);
   cd->Connect("ColorSelected(Pixel_t)", "TGProgressBar", fr,  "SetBarColor(Pixel_t)");

   MapGlobalDialog(cd, fr);
   fClient->WaitForUnmap(cd);
   TQObject::Disconnect(cd);
}

////////////////////////////////////////////////////////////////////////////////
/// Change text color which displays position.

void TGuiBldDragManager::ChangeTextColor(TGProgressBar *fr)
{
   TGGC *gc = gClient->GetResourcePool()->GetGCPool()->FindGC(fr->GetNormGC());

   if (!gc) {
      return;
   }

   Pixel_t pixel = gc->GetForeground();
   TGColorDialog *cd = GetGlobalColorDialog();

   cd->SetCurrentColor(pixel);
   cd->Connect("ColorSelected(Pixel_t)", "TGProgressBar", fr,
               "SetForegroundColor(Pixel_t)");

   MapGlobalDialog(cd, fr);
   fClient->WaitForUnmap(cd);
   TQObject::Disconnect(cd);
}

////////////////////////////////////////////////////////////////////////////////
/// Set text color. This method is invoked
/// via context menu during  guibuilding.

void TGuiBldDragManager::ChangeTextColor(TGTextEntry *fr)
{
   Pixel_t color = fr->GetTextColor();

   TGColorDialog *cd = GetGlobalColorDialog();
   cd->SetCurrentColor(color);
   cd->Connect("ColorSelected(Pixel_t)", "TGTextEntry", fr, "SetTextColor(Pixel_t)");

   MapGlobalDialog(cd, fr);
   fClient->WaitForUnmap(cd);
   TQObject::Disconnect(cd);
}

////////////////////////////////////////////////////////////////////////////////
/// Change text font via font selection dialog. This method is activated
/// via context menu during guibuilding.

void TGuiBldDragManager::ChangeTextFont(TGTextEntry *fr)
{
   TGFontDialog *fd = GetGlobalFontDialog();

   fd->SetColor(fr->GetTextColor());
   FontStruct_t fs = fr->GetFontStruct();
   TGFont *font = fClient->GetResourcePool()->GetFontPool()->FindFont(fs);

   if (font) {
      fd->SetFont(font);
   }

   fd->EnableAlign(kFALSE);
   fd->Connect("FontSelected(char*)", "TGTextEntry", fr, "SetFont(char*)");
   fd->Connect("ColorSelected(Pixel_t)", "TGTextEntry", fr, "SetTextColor(Pixel_t)");

   MapGlobalDialog(fd, fr);
   fClient->WaitForUnmap(fd);
   TQObject::Disconnect(fd);

   int tw, max_ascent, max_descent;
   tw = gVirtualX->TextWidth(fs, fr->GetText(), fr->GetBuffer()->GetTextLength());

   if (tw < 1) {
      TString dummy('w', fr->GetBuffer()->GetBufferLength());
      tw = gVirtualX->TextWidth(fs, dummy.Data(), dummy.Length());
   }

   gVirtualX->GetFontProperties(fs, max_ascent, max_descent);
   fr->Resize(tw + 8, max_ascent + max_descent + 7);
}

////////////////////////////////////////////////////////////////////////////////
/// Invoke file dialog to assign a new image.
/// This method is activated via context menu during guibuilding.

void TGuiBldDragManager::ChangeImage(TGIcon *fr)
{
   static TGFileInfo fi;
   static TString dir(".");
   static Bool_t overwr = kFALSE;
   TString fname;

   fi.fFileTypes = gImageTypes;
   fi.SetIniDir(dir);
   fi.fOverwrite = overwr;

   TGWindow *root = (TGWindow*)fClient->GetRoot();
   SetEditable(kFALSE);

   new TGFileDialog(fClient->GetDefaultRoot(), fr, kFDOpen, &fi);

   if (!fi.fFilename) {
      root->SetEditable(kTRUE);
      gDragManager->SetEditable(kTRUE);
      return;
   }

   dir    = fi.fIniDir;
   overwr = fi.fOverwrite;
   fname  = fi.fFilename;

   TImage *img = TImage::Open(fname.Data());

   if (!img) {
      Int_t retval;
      new TGMsgBox(fClient->GetDefaultRoot(), fr, "Error...",
                   TString::Format("Cannot read image file (%s)", fname.Data()),
                   kMBIconExclamation, kMBRetry | kMBCancel, &retval);

      if (retval == kMBRetry) {
         ChangeImage(fr);
      }
   } else {
      fr->SetImage(img);
      fr->SetImagePath(gSystem->GetDirName(fname.Data()));
   }

   root->SetEditable(kTRUE);
   SetEditable(kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Set lasso drawn flag

void TGuiBldDragManager::SetLassoDrawn(Bool_t on)
{
   if (fLassoDrawn == on) {
      return;
   }

   fLassoDrawn = on;

   if (fBuilder) {
      if (on) {
         fBuilder->EnableEditButtons(kFALSE);
      }

      fBuilder->EnableLassoButtons(on);
   }
}
