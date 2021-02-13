// @(#)root/gui:$Id$
// Author: Fons Rademakers   20/02/98

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TRootDialog                                                          //
//                                                                      //
// A TRootDialog is used to prompt for the arguments of an object's     //
// member function. A TRootDialog is created via the context menu's     //
// when selecting a member function taking arguments.                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TRootDialog.h"
#include "TRootContextMenu.h"
#include "TContextMenu.h"
#include "TClassMenuItem.h"
#include "TList.h"
#include "TGLabel.h"
#include "TGTextEntry.h"
#include "TGButton.h"
#include "TObjString.h"
#include "KeySymbols.h"
#include "TVirtualX.h"

extern TGTextEntry *gBlinkingEntry;

ClassImp(TRootDialog);

////////////////////////////////////////////////////////////////////////////////
/// Create a method argument prompt dialog.

TRootDialog::TRootDialog(TRootContextMenu *cmenu, const TGWindow *main,
    const char *title, Bool_t okB, Bool_t cancelB, Bool_t applyB,
    Bool_t helpB) : TGTransientFrame(gClient->GetRoot(), main, 200, 100)
{
   fMenu   = cmenu;

   fOk     = okB;
   fCancel = cancelB;
   fApply  = applyB;
   fHelp   = helpB;

   fWidgets = new TList;

   fL1 = new TGLayoutHints(kLHintsTop | kLHintsCenterX, 0, 0, 5, 0);
   fL2 = new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5);

   SetWindowName(title);
   SetIconName(title);
   SetEditDisabled(kEditDisable);

   AddInput(kKeyPressMask | kEnterWindowMask | kLeaveWindowMask);
}

////////////////////////////////////////////////////////////////////////////////
/// Delete the dialog.

TRootDialog::~TRootDialog()
{
   fWidgets->Delete();
   delete fWidgets;
   delete fL1;
   delete fL2;
}

////////////////////////////////////////////////////////////////////////////////
/// Add a label and text input field.

void TRootDialog::Add(const char *argname, const char *value, const char *type)
{
   TGLabel      *l = new TGLabel(this, argname);
   TString svalue(value);
   // keep double backslashes (e.g. in case of LateX formatting, like \\gamma)
   svalue.ReplaceAll("\\", "\\\\");
   TGTextBuffer *b = new TGTextBuffer(20); b->AddText(0, svalue.Data());
   TGTextEntry  *t = new TGTextEntry(this, b);

   t->Connect("TabPressed()", "TRootDialog", this, "TabPressed()");

   t->Associate(fMenu);
   t->Resize(260, t->GetDefaultHeight());
   AddFrame(l, fL1);
   AddFrame(t, fL2);

   fWidgets->Add(l);
   fWidgets->Add(t);   // TGTextBuffer will be deleted by TGTextEntry
   fWidgets->Add(new TObjString(type));
}

////////////////////////////////////////////////////////////////////////////////
/// Get parameter string (called by contextmenu after OK or Apply has
/// been selected).

const char *TRootDialog::GetParameters()
{
   static TString params;
   TString param;

   TObjString   *str;
   TObject      *obj;

   Int_t selfobjpos;
   if (fMenu->GetContextMenu()->GetSelectedMenuItem())
      selfobjpos = fMenu->GetContextMenu()->GetSelectedMenuItem()->GetSelfObjectPos();
   else
      selfobjpos = -1;

   params.Clear();
   TIter next(fWidgets);
   Int_t nparam = 0;

   while ((obj = next())) {        // first element is label, skip...
      if (obj->IsA() != TGLabel::Class()) break;
      obj = next();                // get either TGTextEntry or TGComboBox
      str = (TObjString *) next(); // get type string

      nparam++;

      const char *type = str ? str->GetString().Data() : 0;
      const char *data = 0;

      if (obj && obj->IsA() == TGTextEntry::Class())
         data = ((TGTextEntry *) obj)->GetBuffer()->GetString();

      // TODO: Combobox...

      // if necessary, replace the selected object by it's address
      if (selfobjpos == nparam-1) {
         if (params.Length()) params += ",";
         param = TString::Format("(TObject*)0x%lx",
               (Long_t)fMenu->GetContextMenu()->GetSelectedObject());
         params += param;
      }

      if (params.Length()) params += ",";
      if (type && data) {
         if (!strncmp(type, "char*", 5))
            param = TString::Format("\"%s\"", data);
         else
            param = TString::Format("(%s)%s", type, data);
      } else
         param = "0";

      params += param;
   }

   // if selected object is the last argument, have to insert it here
   if (selfobjpos == nparam) {
      if (params.Length()) params += ",";
      param = TString::Format("(TObject*)0x%lx",
            (Long_t)fMenu->GetContextMenu()->GetSelectedObject());
      params += param;
   }

   return params.Data();
}

////////////////////////////////////////////////////////////////////////////////
/// Popup dialog.

void TRootDialog::Popup()
{
   //--- create the OK, Apply and Cancel buttons

   UInt_t  nb = 0, width = 0, height = 0;

   TGHorizontalFrame *hf = new TGHorizontalFrame(this, 60, 20, kFixedWidth);
   TGLayoutHints     *l1 = new TGLayoutHints(kLHintsCenterY | kLHintsExpandX, 5, 5, 0, 0);

   // put hf as last in the list to be deleted
   fWidgets->Add(l1);

   TGTextButton *b;
   if (fOk) {
      b = new TGTextButton(hf, "&OK", 1);
      fWidgets->Add(b);
      b->Associate(fMenu);
      hf->AddFrame(b, l1);
      height = b->GetDefaultHeight();
      width  = TMath::Max(width, b->GetDefaultWidth()); ++nb;
   }
   if (fApply) {
      b = new TGTextButton(hf, "&Apply", 2);
      fWidgets->Add(b);
      b->Associate(fMenu);
      hf->AddFrame(b, l1);
      height = b->GetDefaultHeight();
      width  = TMath::Max(width, b->GetDefaultWidth()); ++nb;
   }
   if (fCancel) {
      b = new TGTextButton(hf, "&Cancel", 3);
      fWidgets->Add(b);
      b->Associate(fMenu);
      hf->AddFrame(b, l1);
      height = b->GetDefaultHeight();
      width  = TMath::Max(width, b->GetDefaultWidth()); ++nb;
   }
   if (fHelp) {
      b = new TGTextButton(hf, "Online &Help", 4);
      fWidgets->Add(b);
      b->Associate(fMenu);
      hf->AddFrame(b, l1);
      height = b->GetDefaultHeight();
      width  = TMath::Max(width, b->GetDefaultWidth()); ++nb;
   }

   // place buttons at the bottom
   l1 = new TGLayoutHints(kLHintsBottom | kLHintsCenterX, 0, 0, 5, 5);
   fWidgets->Add(l1);
   fWidgets->Add(hf);

   AddFrame(hf, l1);

   // keep the buttons centered and with the same width
   hf->Resize((width + 20) * nb, height);

   // map all widgets and calculate size of dialog
   MapSubwindows();

   width  = GetDefaultWidth();
   height = GetDefaultHeight();

   Resize(width, height);

   // position relative to the parent's window
   CenterOnParent();

   // make the message box non-resizable
   SetWMSize(width, height);
   SetWMSizeHints(width, height, width, height, 0, 0);

   SetMWMHints(kMWMDecorAll | kMWMDecorResizeH  | kMWMDecorMaximize |
                              kMWMDecorMinimize | kMWMDecorMenu,
               kMWMFuncAll  | kMWMFuncResize    | kMWMFuncMaximize |
                              kMWMFuncMinimize,
               kMWMInputModeless);

   MapWindow();
   fClient->WaitFor(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Called when closed via window manager action.

void TRootDialog::CloseWindow()
{
   // Send Cancel button message to context menu eventhandler
   SendMessage(fMenu, MK_MSG(kC_COMMAND, kCM_BUTTON), 3, 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Handle Tab keyboard navigation in this dialog.

void TRootDialog::TabPressed()
{
   Bool_t setNext = kFALSE;
   TGTextEntry *entry;
   TIter next(fWidgets);

   while ( TObject* obj = next() ) {
      if ( obj->IsA() == TGTextEntry::Class() ) {
         entry = (TGTextEntry*) obj;
         if ( entry == gBlinkingEntry ) {
            setNext = kTRUE;
         } else if ( setNext ) {
            entry->SetFocus();
            entry->End();
            return;
         }
      }
   }

   next.Reset();
   while ( TObject* obj = next() ) {
      if ( obj->IsA() == TGTextEntry::Class() ) {
         entry = (TGTextEntry*) obj;
         entry->SetFocus();
         entry->End();
         return;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// The key press event handler in this dialog.

Bool_t TRootDialog::HandleKey(Event_t* event)
{
   char   tmp[10];
   UInt_t keysym;
   gVirtualX->LookupString(event, tmp, sizeof(tmp), keysym);
   if ((EKeySym)keysym  == kKey_Tab) {

      TGTextEntry *entry;
      TIter next(fWidgets);

      while ( TObject* obj = next() ) {
         if ( obj->IsA() == TGTextEntry::Class() ) {
            entry = (TGTextEntry*) obj;
            entry->TabPressed();
            return kTRUE;
         }
      }
   }

   return TGMainFrame::HandleKey(event);
}
