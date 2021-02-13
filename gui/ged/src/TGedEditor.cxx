// @(#)root/ged:$Id$
// Author: Marek Biskup, Ilka Antcheva 02/08/2003

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGedEditor                                                           //
//                                                                      //
// The main class of ROOT graphics editor. It manages the appearance    //
// of objects editors according to the selected object in the canvas    //
// (an object became selected after the user click on it using the      //
// left-mouse button).                                                  //
//                                                                      //
// Every object editor provides an object specific GUI and follows a    //
// simple naming convention: it has as a name the object class name     //
// concatinated with 'Editor' (e.g. for TGraph objects the object       //
// editor is TGraphEditor).                                             //
//                                                                      //
// The ROOT graphics editor can be activated by selecting 'Editor'      //
// from the View canvas menu, or SetLine/Fill/Text/MarkerAttributes     //
// from the context menu. The algorithm in use is simple: according to  //
// the selected object <obj> in the canvas it looks for a class name    //
// <obj>Editor. If a class with this name exists, the editor verifies   //
// that this class derives from the base editor class TGedFrame.        //
// It makes an instance of the object editor, scans all object base     //
// classes searching the corresponding object editors and makes an      //
// instance of the base class editor too. Once the object editor is in  //
// place, it sets the user interface elements according to the object   //
// state and is ready for interactions. When a new object of a          //
// different class is selected, a new object editor is loaded in the    //
// editor frame. The old one is cached in memory for potential reuse.   //
//                                                                      //
// Any created canvas will be shown with the editor if you have a       //
// .rootrc file in your working directory containing the the line:      //
// Canvas.ShowEditor:      true                                         //
//                                                                      //
// An created object can be set as selected in a macro by:              //
// canvas->Selected(parent_pad_of_object, object, 1);                   //
// The first parameter can be the canvas itself or the pad containing   //
// 'object'.                                                            //
//                                                                      //
// Begin_Html                                                           //
/*
<img src="gif/TGedEditor.gif">
*/
//End_Html
//////////////////////////////////////////////////////////////////////////

#include "TGedEditor.h"
#include "TCanvas.h"
#include "TGCanvas.h"
#include "TGTab.h"
#include "TGedFrame.h"
#include "TROOT.h"
#include "TClass.h"
#include "TBaseClass.h"
#include "TVirtualX.h"


class TGedTabInfo : public TObject {
   // Helper class for managing visibility and order of created tabs.
public:
   TGTabElement      *fElement;
   TGCompositeFrame  *fContainer;

   TGedTabInfo(TGTabElement* el, TGCompositeFrame* f) :
      fElement(el), fContainer(f) {}
};


ClassImp(TGedEditor);

TGedEditor* TGedEditor::fgFrameCreator = 0;

////////////////////////////////////////////////////////////////////////////////
/// Returns TGedEditor that currently creates TGedFrames.

TGedEditor* TGedEditor::GetFrameCreator()
{
   return fgFrameCreator;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the TGedEditor that currently creates TGedFrames.

void TGedEditor::SetFrameCreator(TGedEditor* e)
{
   fgFrameCreator = e;
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor of graphics editor.

TGedEditor::TGedEditor(TCanvas* canvas, UInt_t width, UInt_t height) :
   TGMainFrame(gClient->GetRoot(), width, height),
   fCan          (0),
   fTab          (0),
   fTabContainer (0),
   fModel        (0),
   fPad          (0),
   fCanvas       (0),
   fClass        (0),
   fGlobal       (kTRUE)
{
   fCan = new TGCanvas(this, 170, 10, kFixedWidth);
   AddFrame(fCan, new TGLayoutHints(kLHintsExpandY | kLHintsExpandX));

   fTab = new TGTab(fCan->GetViewPort(), 10, 10);
   fTab->Associate(fCan);
   fTab->SetCleanup(kDeepCleanup);
   fCan->SetContainer(fTab);

   fTabContainer = GetEditorTab("Style");

   gROOT->GetListOfCleanups()->Add(this);

   SetCanvas(canvas);
   if (fCanvas) {
      UInt_t ch = fCanvas->GetWindowHeight();
      if (ch)
         Resize(GetWidth(), ch > 700 ? 700 : ch);
      else
         Resize(GetWidth(), fCanvas->GetWh()<450 ? 450 : fCanvas->GetWh() + 4);
   } else {
      Resize(width, height);
   }

   MapSubwindows();
   MapWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Editor destructor.

TGedEditor::~TGedEditor()
{
   Hide();

   if(fGlobal){
      TQObject::Disconnect("TCanvas", "Selected(TVirtualPad *, TObject *, Int_t)");
      TQObject::Disconnect("TCanvas", "Closed()");
   }

   // delete class editors
   TIter next(fFrameMap.GetTable());
   TPair* pair;
   while ((pair = (TPair*) next())) {
      if (pair->Value() != 0) {
         TGedFrame* frame  = (TGedFrame*) pair->Value();
         delete frame;
      }
   }

   TGedTabInfo* ti;
   TIter it1(&fCreatedTabs);
   while ((ti = (TGedTabInfo*) it1())) {
      fTab->AddFrame(ti->fElement,0);
      fTab->AddFrame(ti->fContainer,0);
   }

   delete fTab;
   delete ((TGFrameElement*)fList->First())->fLayout;
   delete fCan;
}

////////////////////////////////////////////////////////////////////////////////
/// Virtual method that is called on any change in the dependent frames.
/// This implementation simply calls fPad Modified()/Update().

void TGedEditor::Update(TGedFrame* /*frame*/)
{
   if (fPad) {
      fPad->Modified();
      fPad->Update();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Find or create tab with name.

TGCompositeFrame* TGedEditor::GetEditorTab(const char* name)
{
   return GetEditorTabInfo(name)->fContainer;
}

////////////////////////////////////////////////////////////////////////////////
/// Find or create tab with name.

TGedTabInfo* TGedEditor::GetEditorTabInfo(const char* name)
{
   // look in list of created tabs
   if ( ! fCreatedTabs.IsEmpty()) {
      TIter next(&fCreatedTabs);
      TGedTabInfo* ti;
      while ((ti = (TGedTabInfo *) next())) {
         if (*ti->fElement->GetText() == name)
            return ti;
      }
   }

   // create tab
   TGCompositeFrame* tc = fTab->AddTab(new TGString(name));

   // remove created frame end tab element from the fTab frame
   TGTabElement *te = fTab->GetTabTab(fTab->GetNumberOfTabs() - 1);
   fTab->RemoveFrame(tc);
   fTab->RemoveFrame(te);

   // create a title frame for each tab
   TGedFrame* nf = CreateNameFrame(tc, name);
   if (nf) {
      nf->SetGedEditor(this);
      nf->SetModelClass(0);
      tc->AddFrame(nf, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 2, 2, 2, 2));
   }
   // add to list of created tabs
   TGedTabInfo* ti = new TGedTabInfo(te, tc);
   fCreatedTabs.Add(ti);

   return ti;
}

////////////////////////////////////////////////////////////////////////////////
/// Called when closed via WM close button. Calls Hide().

void TGedEditor::CloseWindow()
{
   Hide();
}

////////////////////////////////////////////////////////////////////////////////
/// Clears windows in editor tab.
/// Unmap and withdraw currently shown frames and thus prepare for
/// construction of a new class layout or destruction.

void TGedEditor::ReinitWorkspace()
{
   TIter next(&fVisibleTabs);
   TGedTabInfo* ti;
   while ((ti = (TGedTabInfo*)next())) {
      TGTabElement     *te = ti->fElement;
      TGCompositeFrame *tc = ti->fContainer;

      fTab->RemoveFrame(te);
      fTab->RemoveFrame(tc);

      TIter frames(tc->GetList());
      frames(); // skip name-frame
      TGFrameElement* fr;
      while ((fr = (TGFrameElement *) frames()) != 0) {
         TGFrame *f = fr->fFrame;
         tc->RemoveFrame(f);
         f->UnmapWindow();
         te->UnmapWindow();
         tc->UnmapWindow();
      }
      fVisibleTabs.Remove(ti);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set editor global.

void TGedEditor::SetGlobal(Bool_t global)
{
   fGlobal = global;
   if (fGlobal) {
      TQObject::Connect("TCanvas", "Selected(TVirtualPad *, TObject *, Int_t)",
                        "TGedEditor", this, "GlobalSetModel(TVirtualPad *, TObject *, Int_t)");

      TQObject::Connect("TCanvas", "Closed()",
                        "TGedEditor", this, "GlobalClosed()");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Delete global editor if no canvas exists.

void TGedEditor::GlobalClosed()
{
   if (gROOT->GetListOfCanvases()->IsEmpty())
      TVirtualPadEditor::Terminate();
}

////////////////////////////////////////////////////////////////////////////////
/// Set canvas to global editor.

void TGedEditor::GlobalSetModel(TVirtualPad *pad, TObject * obj, Int_t ev)
{
   if ((ev != kButton1Down) || !IsMapped() ||
       (obj && obj->InheritsFrom("TColorWheel")))
      return;

   TCanvas* can = pad->GetCanvas();
   // Do nothing if canvas is the same as before or
   // local editor of the canvas is active.
   if (!can || (can == fCanvas || can->GetShowEditor()))
      return;

   Show();
}

////////////////////////////////////////////////////////////////////////////////
/// Connect this editor to the Selected signal of canvas 'c'.

void TGedEditor::ConnectToCanvas(TCanvas *c)
{
   c->Connect("Selected(TVirtualPad*,TObject*,Int_t)", "TGedEditor",
              this, "SetModel(TVirtualPad*,TObject*,Int_t)");
}

////////////////////////////////////////////////////////////////////////////////
/// Disconnect this editor from the Selected signal of fCanvas.

void TGedEditor::DisconnectFromCanvas()
{
   if (fCanvas)
      Disconnect(fCanvas, "Selected(TVirtualPad*,TObject*,Int_t)", this, "SetModel(TVirtualPad*,TObject*,Int_t)");
}

////////////////////////////////////////////////////////////////////////////////
/// Change connection to another canvas.

void TGedEditor::SetCanvas(TCanvas *newcan)
{
   if (fCanvas == newcan) return;

   DisconnectFromCanvas();
   fCanvas = newcan;

   if (!newcan) return;

   SetWindowName(Form("%s_Editor", fCanvas->GetName()));
   fPad = fCanvas->GetSelectedPad();
   if (fPad == 0) fPad = fCanvas;
   ConnectToCanvas(fCanvas);
}

////////////////////////////////////////////////////////////////////////////////
/// Activate object editors according to the selected object.

void TGedEditor::SetModel(TVirtualPad* pad, TObject* obj, Int_t event, Bool_t force)
{
   if ((event != kButton1Down) || (obj && obj->InheritsFrom("TColorWheel")))
      return;

   if (gPad && gPad->GetVirtCanvas()) gPad->GetVirtCanvas()->SetCursor(kWatch);
   gVirtualX->SetCursor(GetId(), gVirtualX->CreateCursor(kWatch));

   fPad = pad;
   if (obj == 0) obj = fPad;

   // keep selected by name
   TGTabElement* seltab = fTab->GetCurrentTab();

   Bool_t mapTabs = kFALSE;
   if (fModel != obj || force) {
      fModel = obj;
      if (fModel == 0 || fModel->IsA() != fClass) {
         ReinitWorkspace();
         mapTabs = kTRUE;
         // add Sytle tab to list of visible tabs
         fVisibleTabs.Add(fCreatedTabs.First());
         if (fModel) {
            fClass = fModel->IsA();
            // build a list of editors
            ActivateEditor(fClass, kTRUE);
         } else {
            fClass = 0;
         }

         // add class editors to fTabContainer
         TGedFrame* gfr;
         TIter ngf(&fGedFrames);
         while ((gfr = (TGedFrame*) ngf()))
            fTabContainer->AddFrame(gfr, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 2, 2, 2, 2));

         fExclMap.Clear();
         fGedFrames.Clear();

         // add visible tabs in fTab
         TIter next(&fVisibleTabs);
         TGedTabInfo* ti;
         while ((ti = (TGedTabInfo *) next())) {
            fTab->AddFrame(ti->fElement,0);
            fTab->AddFrame(ti->fContainer,0);
         }
      }
      ConfigureGedFrames(kTRUE);
   } else {
      ConfigureGedFrames(kFALSE);
   } // end fModel != obj

   if (mapTabs) { // selected object is different class
      TGedTabInfo* ti;
      TIter next(&fVisibleTabs);
      while ((ti = (TGedTabInfo *) next())) {
         ti->fElement->MapWindow();
         ti->fContainer->MapWindow();
      }
      if (seltab == 0 || fTab->SetTab(seltab->GetString(), kFALSE) == kFALSE)
         fTab->SetTab(0, kFALSE);
   }

   if (fGlobal)
      Layout();
   else
      ((TGMainFrame*)GetMainFrame())->Layout();

   if (gPad && gPad->GetVirtCanvas()) gPad->GetVirtCanvas()->SetCursor(kPointer);
   gVirtualX->SetCursor(GetId(), gVirtualX->CreateCursor(kPointer));
}

////////////////////////////////////////////////////////////////////////////////
/// Show editor.

void TGedEditor::Show()
{
   // gPad is setup properly in calling code for global and canvas editor.
   if (gPad) SetCanvas(gPad->GetCanvas());

   if (fCanvas && fGlobal) {
      SetModel(fCanvas->GetClickSelectedPad(), fCanvas->GetClickSelected(), kButton1Down);

      if (fCanvas->GetShowEditor())
         fCanvas->ToggleEditor();

      UInt_t dw = fClient->GetDisplayWidth();
      UInt_t cw = fCanvas->GetWindowWidth();
      UInt_t ch = fCanvas->GetWindowHeight();
      UInt_t cx = (UInt_t)fCanvas->GetWindowTopX();
      UInt_t cy = (UInt_t)fCanvas->GetWindowTopY();
      if (!ch)
         cy = cy + 20;      // embeded canvas protection

      Int_t gedx = 0, gedy = 0;

      if (cw + GetWidth() > dw) {
         gedx = cx + cw - GetWidth();
         gedy = ch - GetHeight();
      } else {
         if (cx > GetWidth())
            gedx = cx - GetWidth() - 20;
         else
            gedx = cx + cw + 10;
         gedy = cy - 20;
      }
      MoveResize(gedx, gedy, GetWidth(), ch > 700 ? 700 : ch);
      SetWMPosition(gedx, gedy);
   } else if (fCanvas) {
      SetModel(fCanvas, fCanvas, kButton1Down);
   }
   MapWindow();
   gVirtualX->RaiseWindow(GetId());

   if (!gROOT->GetListOfCleanups()->FindObject(this))
      gROOT->GetListOfCleanups()->Add(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Hide editor. The editor is put into non-active state.

void TGedEditor::Hide()
{
   UnmapWindow();
   ReinitWorkspace();
   fModel = 0; fClass = 0;
   DisconnectFromCanvas();
   fCanvas = 0; fPad = 0;
   gROOT->GetListOfCleanups()->Remove(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Remove references to fModel in case the fModel is being deleted.
/// Deactivate attribute frames if they point to obj.

void TGedEditor::RecursiveRemove(TObject* obj)
{
   if (obj == fPad) {
      // printf("TGedEditor::RecursiveRemove: %s - pad deleted.\n", locglob);
      SetModel(fCanvas, fCanvas, kButton1Down);
      return;
   }

   if (obj == fModel) {
      // printf("TGedEditor::RecursiveRemove: %s - model deleted.\n", locglob);
      SetModel(fPad, fPad, kButton1Down);
      return;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Searches for GedFrames for given class. In recursive mode look for class
/// editor in its list of bases.

void TGedEditor::ActivateEditor(TClass* cl, Bool_t recurse)
{
   TPair     *pair = (TPair*) fFrameMap.FindObject(cl);
   TClass    *edClass = 0;
   TGedFrame *frame = 0;

   if (pair == 0) {
      edClass = TClass::GetClass(Form("%sEditor", cl->GetName()));

      if (edClass && edClass->InheritsFrom(TGedFrame::Class())) {
         TGWindow *exroot = (TGWindow*) fClient->GetRoot();
         fClient->SetRoot(fTabContainer);
         fgFrameCreator = this;
         frame = reinterpret_cast<TGedFrame*>(edClass->New());
         frame->SetModelClass(cl);
         fgFrameCreator = 0;
         fClient->SetRoot(exroot);
      }
      fFrameMap.Add(cl, frame);
   } else {
      frame =  (TGedFrame*)pair->Value();
   }

   Bool_t exclfr    = kFALSE;
   Bool_t exclbases = kFALSE;

   if (frame) {
      TPair* exclpair = (TPair*) fExclMap.FindObject(cl);
      if (exclpair) {
         exclfr = kTRUE;
         exclbases = (exclpair->Value() != 0);
      }

      if (!exclfr && frame->AcceptModel(fModel)){
         // handle extra tabs in the gedframe
         if (frame->GetExtraTabs()) {
            TIter next(frame->GetExtraTabs());
            TGedFrame::TGedSubFrame* subf;
            while ((subf = (TGedFrame::TGedSubFrame*)next())) {
               // locate the composite frame on created tabs
               TGedTabInfo* ti = GetEditorTabInfo(subf->fName);
               ti->fContainer->AddFrame(subf->fFrame, new TGLayoutHints(kLHintsNormal | kLHintsExpandX));
               if (fVisibleTabs.FindObject(ti) == 0)
                  fVisibleTabs.Add(ti);
            }
         }
         InsertGedFrame(frame);
      }
   }

   if (recurse && !exclbases) {
      if (frame)
         frame->ActivateBaseClassEditors(cl);
      else
         ActivateEditors(cl->GetListOfBases(), recurse);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Searches GedFrames for classes in the given list.

void TGedEditor::ActivateEditors(TList* bcl, Bool_t recurse)
{
   TBaseClass *base;
   TIter next(bcl);

   while ((base = (TBaseClass*) next())) {
      ActivateEditor(base->GetClassPointer(), recurse);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Exclude editor for class cl from current construction.
/// If recurse is true the base-class editors of cl are also excluded.

void  TGedEditor::ExcludeClassEditor(TClass* cl, Bool_t recurse)
{
   TPair* pair = (TPair*) fExclMap.FindObject(cl);
   if (pair) {
      if (recurse && pair->Value() == 0)
         pair->SetValue((TObject*)(Long_t)1); // hack, reuse TObject as Bool_t
   } else {
      fExclMap.Add(cl, (TObject*)(Long_t)(recurse ? 1 : 0));
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Insert GedFrame in fGedFrames list according to priorities.

void TGedEditor::InsertGedFrame(TGedFrame* f)
{
   // printf("%s %s  insert gedframe %s \n", fModel->GetName(), fModel->IsA()->GetName(),f->GetModelClass()->GetName());
   TObjLink* lnk = fGedFrames.FirstLink();
   if (lnk == 0) {
      fGedFrames.Add(f);
      return;
   }
   TGedFrame* cf;
   while (lnk) {
      cf = (TGedFrame*) lnk->GetObject();
      if (f->GetPriority() < cf->GetPriority()) {
         fGedFrames.AddBefore(lnk, f);
         return;
      }
      lnk = lnk->Next();
   }
   fGedFrames.Add(f);
}

////////////////////////////////////////////////////////////////////////////////
/// Call SetModel in class editors.

void TGedEditor::ConfigureGedFrames(Bool_t objChanged)
{
   TGFrameElement *el;

   // Call SetModel on TGedNameFrames (first in the container list)
   // and map extra-tabs.
   TIter vistabs(&fVisibleTabs);
   vistabs(); // skip Style tab
   TGedTabInfo* ti;
   while ((ti = (TGedTabInfo *) vistabs())) {
      TIter fr(ti->fContainer->GetList());
      el = (TGFrameElement*) fr();
      if (el) {
         ((TGedFrame*) el->fFrame)->SetModel(fModel);
         if (objChanged) {
            do {
               el->fFrame->MapSubwindows();
               el->fFrame->Layout();
               el->fFrame->MapWindow();
            } while((el = (TGFrameElement *) fr()));
         }
      }
      ti->fContainer->Layout();
   }

   TIter next(fTabContainer->GetList());
   while ((el = (TGFrameElement *) next())) {
      if ((el->fFrame)->InheritsFrom(TGedFrame::Class())) {
         if (objChanged) {
            el->fFrame->MapSubwindows();
            ((TGedFrame *)(el->fFrame))->SetModel(fModel);
            el->fFrame->Layout();
            el->fFrame->MapWindow();
         } else {
            ((TGedFrame *)(el->fFrame))->SetModel(fModel);
         }
      }
   }
   fTabContainer->Layout();
}

////////////////////////////////////////////////////////////////////////////////
/// Virtual function for creation of top name-frame in each tab.

TGedFrame* TGedEditor::CreateNameFrame(const TGWindow* parent, const char* /*tab_name*/)
{
   return new TGedNameFrame(parent);
}

////////////////////////////////////////////////////////////////////////////////
/// Print contents of fFrameMap.

void TGedEditor::PrintFrameStat()
{
   printf("TGedEditor::PrintFrameStat()\n");
   Int_t sum = 0;
   TIter next(fFrameMap.GetTable());
   TPair* pair;
   while ((pair = (TPair*) next())) {
      if (pair->Value() != 0) {
         TClass* cl  = (TClass*) pair->Key();
         printf("TGedFrame created for %s \n", cl->GetName());
         sum ++;
      }
   }
   printf("SUMMARY: %d editors stored in the local map.\n", sum);
}
