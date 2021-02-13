// @(#)root/eve:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TEveGedEditor.h"
#include "TEveElement.h"
#include "TEveManager.h"

#include "TGToolTip.h"
#include "TGDNDManager.h"
#include "TGMsgBox.h"

#include "TClass.h"
#include "TContextMenu.h"
#include "TVirtualX.h"

/** \class TEveGedEditor
\ingroup TEve
Specialization of TGedEditor for proper update propagation to TEveManager.
*/

ClassImp(TEveGedEditor);

Int_t   TEveGedEditor::fgMaxExtraEditors = 10;
TList  *TEveGedEditor::fgExtraEditors    = new TList;

TContextMenu *TEveGedEditor::fgContextMenu = nullptr;

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TEveGedEditor::TEveGedEditor(TCanvas* canvas, UInt_t width, UInt_t height) :
   TGedEditor(canvas, width, height),
   fElement  (0),
   fObject   (0)
{
   // Remove old name-frame -- it is created in TGedEditor constructor
   // so virtuals are not active yet.
   fTabContainer->RemoveAll();
   TGedFrame* nf = CreateNameFrame(fTabContainer, "Style");
   nf->SetGedEditor(this);
   nf->SetModelClass(0);
   fTabContainer->AddFrame(nf, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 2, 2, 2, 2));

   // Fix priority for TAttMarkerEditor.
   TClass* amClass = TClass::GetClass("TAttMarker");
   TClass* edClass = TClass::GetClass("TAttMarkerEditor");
   TGWindow *exroot = (TGWindow*) fClient->GetRoot();
   fClient->SetRoot(fTabContainer);
   SetFrameCreator(this);
   TGedFrame *frame = reinterpret_cast<TGedFrame*>(edClass->New());
   frame->SetModelClass(amClass);
   {
      Int_t off = edClass->GetDataMemberOffset("fPriority");
      if (off == 0)
         Warning("TEveGedEditor::TEveGedEditor", "Can't fix priority for TAttMarkerEditor.\n");
      else
         * (Int_t*) (((char*)frame) + off) = 1;
   }
   SetFrameCreator(0);
   fClient->SetRoot(exroot);
   fFrameMap.Add(amClass, frame);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TEveGedEditor::~TEveGedEditor()
{
   if (gDebug > 0)
      Info("TEveGedEditor::~TEveGedEditor", "%p going down.", this);
}

////////////////////////////////////////////////////////////////////////////////
/// Called from window-manger close button.
/// Unregister from global list and delete the window.

void TEveGedEditor::CloseWindow()
{
   if (gDebug > 0)
      Info("TEveGedEditor::CloseWindow", "%p closing.", this);

   fgExtraEditors->Remove(this);

   DisplayElement(0);

   if (gDNDManager) {
      if (gDNDManager->GetMainFrame() == this)
         gDNDManager->SetMainFrame(0);
   }
   DeleteWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Create name-frame for a tab.

TGedFrame* TEveGedEditor::CreateNameFrame(const TGWindow* parent, const char* /*tab_name*/)
{
   return new TEveGedNameFrame(parent);
}

////////////////////////////////////////////////////////////////////////////////
/// Return eve-element if it is the model object.

TEveElement* TEveGedEditor::GetEveElement() const
{
   return (fModel == fObject) ? fElement : 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Show a TEveElement in editor.

void TEveGedEditor::DisplayElement(TEveElement* re)
{
   static const TEveException eh("TEveGedEditor::DisplayElement ");

   fElement = re;
   fObject  = fElement ? fElement->GetEditorObject(eh) : 0;
   TGedEditor::SetModel(fPad, fObject, kButton1Down, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Show a TObject in editor.

void TEveGedEditor::DisplayObject(TObject* obj)
{
   fElement = dynamic_cast<TEveElement*>(obj);
   fObject  = obj;
   TGedEditor::SetModel(fPad, obj, kButton1Down, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Set model object.

void TEveGedEditor::SetModel(TVirtualPad* pad, TObject* obj, Int_t event, Bool_t force)
{
   fElement = dynamic_cast<TEveElement*>(obj);
   fObject  = obj;
   TGedEditor::SetModel(pad, obj, event, force);
}

////////////////////////////////////////////////////////////////////////////////
/// Virtual method from TGedEditor, called on every change.
/// Propagates changes to TEveElement and TEveManager.

void TEveGedEditor::Update(TGedFrame* /*gframe*/)
{
   if (fElement)
   {
      fElement->ElementChanged();
      fElement->PropagateVizParamsToProjecteds();
   }

   gEve->Redraw3D();
}

////////////////////////////////////////////////////////////////////////////////
/// Static function to create a new extra editor.

void TEveGedEditor::SpawnNewEditor(TObject* obj)
{
   if (fgExtraEditors->GetSize() >= fgMaxExtraEditors)
   {
      new TGMsgBox(gClient->GetDefaultRoot(), gEve->GetMainWindow(),
                   "Clutter warning",
                   "Maximum number of extra editors reached.",
                   kMBIconStop, kMBOk);
   }

   if (obj)
   {
      TEveGedEditor *ed = new TEveGedEditor();
      ed->DisplayObject(obj);
      ed->SetWindowName(Form("GED %s", obj->GetName()));

      fgExtraEditors->Add(ed);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Element was changed. Update editors showing it.

void TEveGedEditor::ElementChanged(TEveElement* el)
{
   TObject *eobj = el->GetEditorObject("TEveGedEditor::ElementChanged ");
   TObjLink *lnk = fgExtraEditors->FirstLink();
   while (lnk)
   {
      TEveGedEditor *ed = (TEveGedEditor*) lnk->GetObject();
      if (ed->GetModel() == eobj)
         ed->DisplayElement(el);
      lnk = lnk->Next();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Element is being deleted. Close editors showing it.

void TEveGedEditor::ElementDeleted(TEveElement* el)
{
   TObject *eobj = el->GetEditorObject("TEveGedEditor::ElementChanged ");
   TObjLink *lnk = fgExtraEditors->FirstLink();
   while (lnk)
   {
      TEveGedEditor *ed = (TEveGedEditor*) lnk->GetObject();
      if (ed->GetModel() == eobj)
      {
         TObjLink *next = lnk->Next();
         ed->DeleteWindow();
         fgExtraEditors->Remove(lnk);
         lnk = next;
      }
      else
      {
         lnk = lnk->Next();
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Destroys all editors. Called from EVE termination.

void TEveGedEditor::DestroyEditors()
{
   while ( ! fgExtraEditors->IsEmpty())
   {
      TEveGedEditor *ed = (TEveGedEditor*) fgExtraEditors->First();
      ed->DeleteWindow();
      fgExtraEditors->RemoveFirst();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return context menu object shared among eve-ged-editors.

TContextMenu* TEveGedEditor::GetContextMenu()
{
   if (fgContextMenu == 0)
      fgContextMenu = new TContextMenu("", "");
   return fgContextMenu;
}

/** \class TEveGedNameFrame
\ingroup TEve
Specialization of TGedNameFrame used in EVE.
It provides the ability to undock given editor for easier use.
Support for that is also provided from the TEveManager.
*/

ClassImp(TEveGedNameFrame);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TEveGedNameFrame::TEveGedNameFrame(const TGWindow *p, Int_t width, Int_t height,
                                   UInt_t options) :
   TGedFrame(p, width, height, options),
   fNCButton(0)
{
   fNCButton = new TEveGedNameTextButton(this);
   fNCButton->SetTextColor(0x0020a0);
   AddFrame(fNCButton, new TGLayoutHints(kLHintsNormal | kLHintsExpandX));
   fNCButton->Connect("Clicked()", "TEveGedNameFrame", this, "SpawnEditorClone()");
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TEveGedNameFrame::~TEveGedNameFrame()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Set model object.

void TEveGedNameFrame::SetModel(TObject* obj)
{
   if (obj)
   {
      fNCButton->SetText(Form("%s [%s]", obj->GetName(), obj->ClassName()));
      fNCButton->SetToolTipText(obj->GetTitle());
      fNCButton->SetEnabled(kTRUE);
   }
   else
   {
      fNCButton->SetText("No object selected.");
      fNCButton->SetToolTipText(0);
      fNCButton->SetEnabled(kFALSE);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Create a new floating editor with current object.

void TEveGedNameFrame::SpawnEditorClone()
{
   TEveGedEditor::SpawnNewEditor(fGedEditor->GetModel());
}

/** \class TEveGedNameTextButton
\ingroup TEve
Specialization of TGTextButton for EVE name frame.
It opens a context-menu on right-click.
*/

ClassImp(TEveGedNameTextButton);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TEveGedNameTextButton::TEveGedNameTextButton(TEveGedNameFrame* p) :
   TGTextButton(p, ""),
   fFrame(p)
{
   gVirtualX->GrabButton(fId, kAnyButton, kAnyModifier,
                         kButtonPressMask | kButtonReleaseMask,
                         kNone, kNone);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TEveGedNameTextButton::~TEveGedNameTextButton()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Handle button.

Bool_t TEveGedNameTextButton::HandleButton(Event_t* event)
{
   static const TEveException eh("TEveGedNameTextButton::HandleButton ");

   if (fTip) fTip->Hide();
   if (fState == kButtonDisabled) return kTRUE;

   if (event->fCode == kButton3 && event->fType == kButtonPress)
   {
      TEveGedEditor *eged = (TEveGedEditor*) fFrame->GetGedEditor();
      TEveElement   *el   = eged->GetEveElement();
      if (el)
         TEveGedEditor::GetContextMenu()->Popup(event->fXRoot, event->fYRoot,
                                                el->GetObject(eh));
      return 1;
   }
   else if (event->fCode == kButton1)
   {
      return TGTextButton::HandleButton(event);
   }
   else
   {
      return 0;
   }
}
