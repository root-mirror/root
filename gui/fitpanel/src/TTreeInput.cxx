// @(#)root/gui:$Id$
// Author: David Gonzalez Maline  21/10/2008

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Tree Input Widget                                                    //
//                                                                      //
// An dialog box that asks the user for the variables and cuts          //
// of the selected tree in the fitpanel.                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TTreeInput.h"
#include "TGButton.h"
#include "TGLabel.h"
#include "TGTextEntry.h"
#include "strlcpy.h"

enum ETreeInput {
   kTI_TEVARS, kTI_TECUTS
};

ClassImp(TTreeInput);

////////////////////////////////////////////////////////////////////////////////
/// Create simple input dialog.

   TTreeInput::TTreeInput(const TGWindow *p, const TGWindow *main,
                          char *strvars, char *strcuts):
      TGTransientFrame(p, main, 10, 10, kVerticalFrame),
      fStrvars(strvars),
      fStrcuts(strcuts)
{
   if (!p && !main) {
      MakeZombie();
      return;
   }
   SetCleanup(kDeepCleanup);

   TGLabel *label = new TGLabel(this, "Selected Variables: ");
   AddFrame(label, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 0));

   TGTextBuffer *tbuf = new TGTextBuffer(256);  //will be deleted by TGtextEntry
   fTEVars = new TGTextEntry(this, tbuf, kTI_TEVARS);
   fTEVars->Resize(260, fTEVars->GetDefaultHeight());
   AddFrame(fTEVars, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5));

   label = new TGLabel(this, "Selected Cuts: ");
   AddFrame(label, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 0));

   tbuf = new TGTextBuffer(256);  //will be deleted by TGtextEntry
   fTECuts = new TGTextEntry(this, tbuf, kTI_TECUTS);
   fTECuts->Resize(260, fTECuts->GetDefaultHeight());
   AddFrame(fTECuts, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5));

   // create frame and layout hints for Ok and Cancel buttons
   TGHorizontalFrame *hf = new TGHorizontalFrame(this, 60, 20, kFixedWidth);
   hf->SetCleanup(kDeepCleanup);

   // create OK and Cancel buttons in their own frame (hf)
   UInt_t  width = 0, height = 0;

   fOk = new TGTextButton(hf, "&Ok", 1);
   fOk->Associate(this);
   hf->AddFrame(fOk, new TGLayoutHints(kLHintsCenterY | kLHintsExpandX, 5, 5, 0, 0));
   height = fOk->GetDefaultHeight();
   width  = TMath::Max(width, fOk->GetDefaultWidth());

   fCancel = new TGTextButton(hf, "&Cancel", 2);
   fCancel->Associate(this);
   hf->AddFrame(fCancel, new TGLayoutHints(kLHintsCenterY | kLHintsExpandX, 5, 5, 0, 0));
   height = fCancel->GetDefaultHeight();
   width  = TMath::Max(width, fCancel->GetDefaultWidth());

   // place button frame (hf) at the bottom
   AddFrame(hf, new TGLayoutHints(kLHintsBottom | kLHintsCenterX, 0, 0, 5, 5));

   // keep buttons centered and with the same width
   hf->Resize((width + 20) * 2, height);

   // set dialog title
   SetWindowName("Get Input");

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

   // popup dialog and wait till user replies
   MapWindow();
   fTEVars->SetFocus();

   gClient->WaitFor(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Cleanup dialog.

TTreeInput::~TTreeInput()
{
   Cleanup();
}

////////////////////////////////////////////////////////////////////////////////
/// Handle button and text enter events

Bool_t TTreeInput::ProcessMessage(Long_t msg, Long_t parm1, Long_t)
{
   switch (GET_MSG(msg)) {
      case kC_COMMAND:
         switch (GET_SUBMSG(msg)) {
            case kCM_BUTTON:
               switch (parm1) {
                  case 1:
                     // here copy the string from text buffer to return variable
                     // see TFitEditor.cxx for the maximum length:
                     // char variables[256] = {0}; char cuts[256] = {0};
                     strlcpy(fStrvars, fTEVars->GetBuffer()->GetString(), 256);
                     strlcpy(fStrcuts, fTECuts->GetBuffer()->GetString(), 256);
                     delete this;
                     break;
                  case 2:
                     fStrvars[0] = 0;
                     fStrcuts[0] = 0;
                     delete this;
                     break;
               }
               default:
                  break;
         }
         break;

      case kC_TEXTENTRY:
         switch (GET_SUBMSG(msg)) {
            case kTE_ENTER:
               // here copy the string from text buffer to return variable
               // see TFitEditor.cxx for the maximum length:
               // char variables[256] = {0}; char cuts[256] = {0};
               strlcpy(fStrvars, fTEVars->GetBuffer()->GetString(), 256);
               strlcpy(fStrcuts, fTECuts->GetBuffer()->GetString(), 256);
               delete this;
               break;
            case kTE_TAB:
               if ( parm1 == kTI_TEVARS )
                  fTECuts->SetFocus();
               else if ( parm1 == kTI_TECUTS )
                  fTEVars->SetFocus();
               break;
            default:
               break;
         }
         break;

      default:
         break;
   }
   return kTRUE;
}
