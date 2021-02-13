// @(#)root/gui:$Id$
// Author: Fons Rademakers   10/7/2000

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/
/**************************************************************************

    This source is based on Xclass95, a Win95-looking GUI toolkit.
    Copyright (C) 1996, 1997 David Barth, Ricky Ralston, Hector Peraza.

    Xclass95 is free software; you can redistribute it and/or
    modify it under the terms of the GNU Library General Public
    License as published by the Free Software Foundation; either
    version 2 of the License, or (at your option) any later version.

**************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGTextEditDialogs                                                    //
//                                                                      //
// This file defines several dialogs that are used by the TGTextEdit    //
// widget via its associated context popup menu.                        //
// The following dialogs are available: TGSearchDialog, TGGotoDialog    //
// and TGPrintDialog.                                                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TGTextEditDialogs.h"
#include "TGButton.h"
#include "TGLabel.h"
#include "TGTextEntry.h"
#include "TGIcon.h"
#include "TGMsgBox.h"
#include "TGComboBox.h"
#include "TSystem.h"
#include "TObjArray.h"
#include "TVirtualX.h"
#include "strlcpy.h"
#include "snprintf.h"

#include <cstdlib>


ClassImp(TGSearchDialog);
ClassImp(TGPrintDialog);
ClassImp(TGGotoDialog);

static TString gLastSearchString;
TGSearchDialog *TGSearchDialog::fgSearchDialog = nullptr;

////////////////////////////////////////////////////////////////////////////////
/// Create a search dialog box. Used to get from the user the required
/// search instructions. Ret_code is kTRUE when sstruct has been set,
/// kFALSE otherwise (like when dialog was canceled).

TGSearchDialog::TGSearchDialog(const TGWindow *p, const TGWindow *main,
                               UInt_t w, UInt_t h, TGSearchType *sstruct,
                               Int_t *ret_code, UInt_t options) :
     TGTransientFrame(p, main, w, h, options)
{
   if (!p && !main) {
      MakeZombie();
      // coverity [uninit_ctor]
      return;
   }
   fRetCode = ret_code;
   fType = sstruct;

   ChangeOptions((GetOptions() & ~kVerticalFrame) | kHorizontalFrame);

   fF1 = new TGCompositeFrame(this, 60, 20, kVerticalFrame | kFixedWidth);
   fF2 = new TGCompositeFrame(this, 60, 20, kVerticalFrame);
   fF3 = new TGCompositeFrame(fF2, 60, 20, kHorizontalFrame);
   fF4 = new TGCompositeFrame(fF2, 60, 20, kHorizontalFrame);

   fSearchButton = new TGTextButton(fF1, new TGHotString("&Search"), 1);
   fCancelButton = new TGTextButton(fF1, new TGHotString("&Cancel"), 2);
   fF1->Resize(fSearchButton->GetDefaultWidth()+40, GetDefaultHeight());

   fSearchButton->Associate(this);
   fCancelButton->Associate(this);

   fL1 = new TGLayoutHints(kLHintsTop | kLHintsExpandX, 2, 2, 3, 0);
   fL2 = new TGLayoutHints(kLHintsTop | kLHintsRight | kLHintsExpandX,
                           2, 5, 0, 2);
   fL21 = new TGLayoutHints(kLHintsTop | kLHintsRight, 2, 5, 10, 0);

   fF1->AddFrame(fSearchButton, fL1);
   fF1->AddFrame(fCancelButton, fL1);

   AddFrame(fF1, fL21);

   fLSearch = new TGLabel(fF3, new TGHotString("Search &for:"));

   fCombo = new TGComboBox(fF3, "");
   fSearch = fCombo->GetTextEntry();
   fBSearch = fSearch->GetBuffer();
   if (sstruct && sstruct->fBuffer)
      fBSearch->AddText(0, sstruct->fBuffer);
   else if (!gLastSearchString.IsNull())
      fBSearch->AddText(0, gLastSearchString.Data());
   else
      fSearchButton->SetState(kButtonDisabled);
   fSearch->Associate(this);
   fCombo->Resize(220, fSearch->GetDefaultHeight());
   fSearch->SelectAll();

   fL5 = new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 3, 5, 0, 0);
   fL6 = new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 0, 2, 0, 0);

   fF3->AddFrame(fLSearch, fL5);
   fF3->AddFrame(fCombo, fL6);

   fG2 = new TGGroupFrame(fF4, new TGString("Direction"), kHorizontalFrame);

   fL3 = new TGLayoutHints(kLHintsTop | kLHintsRight, 2, 2, 2, 2);
   fL9 = new TGLayoutHints(kLHintsBottom | kLHintsLeft, 0, 0, 0, 0);
   fL4 = new TGLayoutHints(kLHintsBottom | kLHintsLeft, 0, 0, 5, 0);
   fL10 = new TGLayoutHints(kLHintsBottom | kLHintsRight, 0, 0, 5, 0);

   fCaseCheck = new TGCheckButton(fF4, new TGHotString("&Case sensitive"), 1);
   fCaseCheck->Associate(this);
   fF4->AddFrame(fCaseCheck, fL9);

   fDirectionRadio[0] = new TGRadioButton(fG2, new TGHotString("Forward"),  1);
   fDirectionRadio[1] = new TGRadioButton(fG2, new TGHotString("Backward"), 2);

   fG2->AddFrame(fDirectionRadio[0], fL4);
   fG2->AddFrame(fDirectionRadio[1], fL10);
   fDirectionRadio[0]->Associate(this);
   fDirectionRadio[1]->Associate(this);

   if (fType->fCaseSensitive == kFALSE)
      fCaseCheck->SetState(kButtonUp);
   else
      fCaseCheck->SetState(kButtonDown);

   if (fType->fDirection)
      fDirectionRadio[0]->SetState(kButtonDown);
   else
      fDirectionRadio[1]->SetState(kButtonDown);

   fF4->AddFrame(fG2, fL3);

   fF2->AddFrame(fF3, fL1);
   fF2->AddFrame(fF4, fL1);

   AddFrame(fF2, fL2);

   MapSubwindows();
   Resize(GetDefaultSize());
   SetEditDisabled(kEditDisable);

   CenterOnParent();

   SetWindowName("Search");
   SetIconName("Search");

   SetMWMHints(kMWMDecorAll | kMWMDecorMaximize | kMWMDecorMenu,
               kMWMFuncAll | kMWMFuncMaximize | kMWMFuncResize,
               kMWMInputModeless);

   if (fType->fClose) {
      MapWindow();
      fSearch->RequestFocus();
      fClient->WaitFor(this);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Clean up search dialog.

TGSearchDialog::~TGSearchDialog()
{
   if (IsZombie()) return;
   delete fSearchButton;
   delete fCancelButton;
   delete fDirectionRadio[0]; delete fDirectionRadio[1];
   delete fCaseCheck;
   delete fCombo;
   delete fLSearch;
   delete fG2;
   delete fF1; delete fF2; delete fF3; delete fF4;
   delete fL1; delete fL2; delete fL3; delete fL4; delete fL5; delete fL6;
   delete fL21;delete fL9; delete fL10;
}

////////////////////////////////////////////////////////////////////////////////
/// Close the dialog. On close the dialog will be deleted and cannot be
/// re-used.

void TGSearchDialog::CloseWindow()
{
   if (fType->fClose) {
      DeleteWindow();
   } else {
      UnmapWindow();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// emit signal when search text entered

void TGSearchDialog::TextEntered(const char *text)
{
   Emit("TextEntered(const char *)", text);
}

////////////////////////////////////////////////////////////////////////////////
/// Process search dialog widget messages.

Bool_t TGSearchDialog::ProcessMessage(Long_t msg, Long_t parm1, Long_t)
{
   const char *string;

   switch (GET_MSG(msg)) {
      case kC_COMMAND:
         switch (GET_SUBMSG(msg)) {
            case kCM_BUTTON:
               switch (parm1) {
                  case 1:
                     string = fBSearch->GetString();
                     if (fType->fBuffer)
                        delete [] fType->fBuffer;
                     fType->fBuffer = StrDup(string);
                     gLastSearchString = string;
                     *fRetCode = kTRUE;
                     TextEntered(string);
                     fCombo->ReturnPressed();
                     if (fType->fClose) CloseWindow();
                     break;
                  case 2:
                     *fRetCode = kFALSE;
                     CloseWindow();
                     break;
               }
               break;

            case kCM_CHECKBUTTON:
               fType->fCaseSensitive = !fType->fCaseSensitive;
               break;

            case kCM_RADIOBUTTON:
               switch (parm1) {
                  case 1:
                     fType->fDirection = kTRUE;
                     fDirectionRadio[1]->SetState(kButtonUp);
                     break;
                  case 2:
                     fType->fDirection = kFALSE;
                     fDirectionRadio[0]->SetState(kButtonUp);
                     break;
               }
               break;

            default:
               break;
         }
         break;

      case kC_TEXTENTRY:
         switch (GET_SUBMSG(msg)) {
            case kTE_TEXTCHANGED:
               string = fBSearch->GetString();
               if (!string[0]) {
                  fSearchButton->SetState(kButtonDisabled);
               } else {
                  fSearchButton->SetState(kButtonUp);
               }
               break;
            case kTE_ENTER:
               string = fBSearch->GetString();
               if (fType->fBuffer)
                  delete [] fType->fBuffer;
               fType->fBuffer = StrDup(string);
               gLastSearchString = string;
               *fRetCode = kTRUE;
               TextEntered(string);
               if (fType->fClose) CloseWindow();
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

////////////////////////////////////////////////////////////////////////////////
/// Return global search dialog.

TGSearchDialog *&TGSearchDialog::SearchDialog()
{
   return fgSearchDialog;
}


////////////////////////////////////////////////////////////////////////////////
/// Create the printer dialog box. Returns kTRUE in ret_code when
/// printerName and printProg have been set and cancel was not pressed,
/// kFALSE otherwise.

TGPrintDialog::TGPrintDialog(const TGWindow *p, const TGWindow *main,
                             UInt_t w, UInt_t h, char **printerName,
                             char **printProg, Int_t *ret_code,
                             UInt_t options) :
   TGTransientFrame(p, main, w, h, options)
{
   if (!p && !main) {
      MakeZombie();
      // coverity [uninit_ctor]
      return;
   }
   fPrinter      = printerName;
   fPrintCommand = printProg;
   fRetCode      = ret_code;
   fEditDisabled = kEditDisable;

   ChangeOptions((GetOptions() & ~kVerticalFrame) | kHorizontalFrame);

   fF1 = new TGCompositeFrame(this, 60, 20, kVerticalFrame | kFixedWidth);
   fF5 = new TGCompositeFrame(this, 60, 20, kHorizontalFrame);
   fF2 = new TGCompositeFrame(fF5,  60, 20, kVerticalFrame);
   fF3 = new TGCompositeFrame(fF2,  60, 20, kHorizontalFrame);
   fF4 = new TGCompositeFrame(fF2,  60, 20, kHorizontalFrame);

   fPrintButton  = new TGTextButton(fF1, new TGHotString("&Print"), 1);
   fCancelButton = new TGTextButton(fF1, new TGHotString("&Cancel"), 2);
   fF1->Resize(fPrintButton->GetDefaultWidth()+40, GetDefaultHeight());

   fPrintButton->Associate(this);
   fCancelButton->Associate(this);

   fL1 = new TGLayoutHints(kLHintsTop | kLHintsExpandX, 2, 2, 2, 2);
   fL2 = new TGLayoutHints(kLHintsTop | kLHintsRight | kLHintsExpandX,
                           2, 5, 0, 2);
   fL3 = new TGLayoutHints(kLHintsTop | kLHintsRight, 2, 2, 4, 4);
   fL5 = new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 3, 5, 0, 0);
   fL6 = new TGLayoutHints(kLHintsExpandX | kLHintsCenterY, 0, 2, 0, 0);
   fL7 = new TGLayoutHints(kLHintsLeft, 10, 10, 10, 10);

   fL21 = new TGLayoutHints(kLHintsTop | kLHintsRight, 2, 5, 10, 0);

   fF1->AddFrame(fPrintButton, fL1);
   fF1->AddFrame(fCancelButton, fL1);
   AddFrame(fF1, fL3);


   fLPrintCommand = new TGLabel(fF3, new TGHotString("Print command:"));
   fBPrintCommand = new TGTextBuffer(50);
   if ((printProg) && (*printProg))
      fBPrintCommand->AddText(0, *printProg);
   fPrintCommandEntry = new TGTextEntry(fF3, fBPrintCommand);
   fPrintCommandEntry->Associate(this);
   fPrintCommandEntry->Resize(150, fPrintCommandEntry->GetDefaultHeight());

   fF3->AddFrame(fLPrintCommand, fL5);
   fF3->AddFrame(fPrintCommandEntry, fL6);

   fLPrinter = new TGLabel(fF4, new TGHotString("Printer:"));
   if ((printerName) && (*printerName))
      fPrinterEntry = new TGComboBox(fF4, *printerName);
   fBPrinter = fPrinterEntry->GetTextEntry()->GetBuffer();
   fPrinterEntry->Resize(150, fPrinterEntry->GetTextEntry()->GetDefaultHeight());
   fF4->AddFrame(fLPrinter, fL5);
   fF4->AddFrame(fPrinterEntry, fL6);

   fF2->AddFrame(fF3, fL1);
   fF2->AddFrame(fF4, fL1);

   const TGPicture *printerPicture = fClient->GetPicture("printer_s.xpm");
   if (!printerPicture) {
      Error("TGPrintDialog", "printer_s.xpm not found");
      fPrinterIcon = 0;
   } else {
      fPrinterIcon = new TGIcon(fF5, printerPicture, 32, 32);
      fF5->AddFrame(fPrinterIcon, fL7);
   }
   fF5->AddFrame(fF2, fL1);
   AddFrame(fF5, fL1);

   MapSubwindows();
   Resize(GetDefaultSize());

   GetPrinters();
   CenterOnParent();

   SetWindowName("Print");
   SetIconName("Print");

   SetMWMHints(kMWMDecorAll | kMWMDecorMaximize | kMWMDecorMenu,
               kMWMFuncAll | kMWMFuncMaximize | kMWMFuncResize,
               kMWMInputModeless);

   fPrinterEntry->RequestFocus();
   MapWindow();
   fClient->WaitFor(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Clean up print dialog.

TGPrintDialog::~TGPrintDialog()
{
   if (IsZombie()) return;
   delete fPrinterIcon;
   delete fPrintButton;
   delete fCancelButton;
   delete fPrinterEntry;       // deletes also fBPrinter
   delete fPrintCommandEntry;  // deletes also fBPrintCommand
   delete fLPrinter; delete fLPrintCommand;
   delete fF1; delete fF2; delete fF3; delete fF4; delete fF5;
   delete fL1; delete fL2; delete fL3; delete fL5; delete fL6; delete fL7;
   delete fL21;
}

////////////////////////////////////////////////////////////////////////////////
/// Close the dialog. On close the dialog will be deleted and cannot be
/// re-used.

void TGPrintDialog::CloseWindow()
{
   DeleteWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Ask the system fo the list of available printers and populate the combo
/// box. If there is a default printer, select it in the list.

void TGPrintDialog::GetPrinters()
{
   TObject *obj;
   Int_t idx = 1, dflt =1;

   if (gVirtualX->InheritsFrom("TGX11") || gVirtualX->InheritsFrom("TGCocoa")) {
      char *lpstat = gSystem->Which(gSystem->Getenv("PATH"), "lpstat",
                                    kExecutePermission);
      if (lpstat == 0) return;
      TString defaultprinter = gSystem->GetFromPipe("lpstat -d");
      TString printerlist = gSystem->GetFromPipe("lpstat -v");
      TObjArray *tokens = printerlist.Tokenize("\n");
      TIter iter(tokens);
      while((obj = iter())) {
         TString line = obj->GetName();
         TObjArray *tk = line.Tokenize(" ");
         TString pname = ((TObject*)tk->At(2))->GetName();
         if (pname.EndsWith(":")) pname.Remove(pname.Last(':'));
         //if (pname.Contains(":")) pname.Remove(pname.Last(':'));
         if (defaultprinter.Contains(pname)) {
            dflt = idx;
            fPrinterEntry->GetTextEntry()->SetText(pname.Data(), kFALSE);
         }
         fPrinterEntry->AddEntry(pname.Data(), idx++);
      }
      delete [] lpstat;
   }
   else {
      TString defaultprinter = gSystem->GetFromPipe("WMIC Path Win32_Printer where Default=TRUE Get DeviceID");
      TString printerlist = gSystem->GetFromPipe("WMIC Path Win32_Printer Get DeviceID");
      defaultprinter.Remove(0, defaultprinter.First('\n')); // remove "Default"
      printerlist.Remove(0, printerlist.First('\n')); // remove "DeviceID"
      printerlist.ReplaceAll("\r", "");
      TObjArray *tokens = printerlist.Tokenize("\n");
      TIter iter(tokens);
      while((obj = iter())) {
         TString pname = obj->GetName();
         pname.Remove(TString::kTrailing, ' ');
         if (defaultprinter.Contains(pname)) {
            dflt = idx;
            fPrinterEntry->GetTextEntry()->SetText(pname.Data(), kFALSE);
         }
         fPrinterEntry->AddEntry(pname.Data(), idx++);
      }
   }
   fPrinterEntry->Select(dflt, kFALSE);
   fPrinterEntry->Layout();
}

////////////////////////////////////////////////////////////////////////////////
/// Process print dialog widget messages.

Bool_t TGPrintDialog::ProcessMessage(Long_t msg, Long_t parm1, Long_t)
{
   const char *string, *txt;

   switch (GET_MSG(msg)) {
      case kC_COMMAND:
         switch (GET_SUBMSG(msg)) {
            case kCM_BUTTON:
               switch (parm1) {
                  case 1:
                     *fRetCode = kTRUE;
                     {
                        string = fBPrinter->GetString();
                        delete [] *fPrinter;
                        const size_t prSize = strlen(string) + 1;
                        *fPrinter = new char[prSize];
                        strlcpy(*fPrinter, string, prSize);
                     }
                     {
                        string = fBPrintCommand->GetString();
                        delete [] *fPrintCommand;
                        const size_t cmdSize = strlen(string) + 1;
                        *fPrintCommand = new char[cmdSize];
                        strlcpy(*fPrintCommand, string, cmdSize);
                     }

                     if (fBPrintCommand->GetTextLength() == 0) {
                        txt = "Please provide print command or use \"Cancel\"";
                        new TGMsgBox(fClient->GetRoot(), GetMainFrame(),
                                     "Missing Print Parameters", txt, kMBIconExclamation,
                                      kMBOk);
                        return kTRUE;
                     }
                     CloseWindow();
                     break;
                  case 2:
                     *fRetCode = kFALSE;
                     CloseWindow();
                     break;
               }
               break;
         }
         break;

      default:
         break;
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a dialog to GoTo a specific line number. Returns -1 in
/// ret_code in case no valid line number was given or in case
/// cancel was pressed. If on input *ret_code is > 0 then this value
/// will be used as default value.

TGGotoDialog::TGGotoDialog(const TGWindow *p, const TGWindow *main,
                           UInt_t w, UInt_t h, Long_t *ret_code,
                           UInt_t options) :
   TGTransientFrame(p, main, w, h, options)
{
   if (!p && !main) {
      MakeZombie();
      // coverity [uninit_ctor]
      return;
   }
   fRetCode = ret_code;
   fEditDisabled = kEditDisable;

   ChangeOptions((GetOptions() & ~kVerticalFrame) | kHorizontalFrame);

   fF1 = new TGCompositeFrame(this, 60, 20, kVerticalFrame | kFixedWidth);
   fF2 = new TGCompositeFrame(this, 60, 20, kHorizontalFrame);

   fGotoButton = new TGTextButton(fF1, new TGHotString("&Goto"), 1);
   fCancelButton = new TGTextButton(fF1, new TGHotString("&Cancel"), 2);
   fF1->Resize(fGotoButton->GetDefaultWidth()+40, GetDefaultHeight());

   fGotoButton->Associate(this);
   fCancelButton->Associate(this);

   fL1 = new TGLayoutHints(kLHintsCenterY | kLHintsExpandX, 2, 2, 3, 0);
   fL21 = new TGLayoutHints(kLHintsCenterY | kLHintsRight, 2, 5, 10, 0);

   fF1->AddFrame(fGotoButton, fL1);
   fF1->AddFrame(fCancelButton, fL1);
   AddFrame(fF1, fL21);

   fLGoTo = new TGLabel(fF2, new TGHotString("&Goto Line:"));

   fBGoTo = new TGTextBuffer(50);
   if (*fRetCode > 0) {
      char curline[32];
      snprintf(curline, 32, "%ld", *fRetCode);
      fBGoTo->AddText(0, curline);
   } else
      fGotoButton->SetState(kButtonDisabled);
   fGoTo = new TGTextEntry(fF2, fBGoTo);
   fGoTo->Associate(this);
   fGoTo->Resize(220, fGoTo->GetDefaultHeight());
   fGoTo->SelectAll();

   fL5 = new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 3, 5, 0, 0);
   fL6 = new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 0, 2, 0, 0);

   fF2->AddFrame(fLGoTo, fL5);
   fF2->AddFrame(fGoTo, fL5);
   AddFrame(fF2, fL1);

   MapSubwindows();
   Resize(GetDefaultSize());

   CenterOnParent();

   SetWindowName("Goto");
   SetIconName("Print");

   SetMWMHints(kMWMDecorAll | kMWMDecorMaximize | kMWMDecorMenu,
               kMWMFuncAll | kMWMFuncMaximize | kMWMFuncResize,
               kMWMInputModeless);

   MapWindow();
   fGoTo->RequestFocus();
   fClient->WaitFor(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Clean up goto dialog

TGGotoDialog::~TGGotoDialog()
{
   if (IsZombie()) return;
   delete fGotoButton;
   delete fCancelButton;
   delete fGoTo;
   delete fLGoTo;
   delete fF1; delete fF2;
   delete fL1; delete fL5; delete fL6; delete fL21;
}

////////////////////////////////////////////////////////////////////////////////
/// Close the dialog. On close the dialog will be deleted and cannot be
/// re-used.

void TGGotoDialog::CloseWindow()
{
   DeleteWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Process goto dialog widget messages.

Bool_t TGGotoDialog::ProcessMessage(Long_t msg, Long_t parm1, Long_t)
{
   const char *string;

   switch (GET_MSG(msg)) {
      case kC_COMMAND:
         switch (GET_SUBMSG(msg)) {
            case kCM_BUTTON:
               switch (parm1) {
                  case 1:
                     string = fBGoTo->GetString();
                     *fRetCode = (Long_t) atof(string);
                     CloseWindow();
                     break;
                  case 2:
                     *fRetCode = -1;
                     CloseWindow();
                     break;
               }
               break;

            default:
               break;
         }
         break;

      case kC_TEXTENTRY:
         switch (GET_SUBMSG(msg)) {
            case kTE_TEXTCHANGED:
               string = fBGoTo->GetString();
               if (!string[0])
                  fGotoButton->SetState(kButtonDisabled);
               else
                  fGotoButton->SetState(kButtonUp);
               break;
            case kTE_ENTER:
               string = fBGoTo->GetString();
               *fRetCode = (Long_t) atof(string);
               CloseWindow();
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

