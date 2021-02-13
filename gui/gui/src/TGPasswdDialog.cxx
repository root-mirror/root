// @(#)root/gui:$Id$
// Author: G. Ganis  10/10/2005

/*************************************************************************
 * Copyright (C) 1995-2005, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGPasswdDialog                                                       //
//                                                                      //
// Graphic dialog to enter passwords                                    //
//                                                                      //
// Usage:                                                               //
//                                                                      //
// {                                                                    //
//   // Buffer for the passwd                                           //
//   char pwdbuf[128]                                                   //
//                                                                      //
//   Open the dialog box                                                //
//   TGPasswdDialog dialog("My prompt", pwdbuf, 128);                   //
//                                                                      //
//   // Wait until the user is done                                     //
//   while (gROOT->IsInterrupted())                                     //
//      gSystem->DispatchOneEvent(kFALSE);                              //
//                                                                      //
//   // Password is now in pwdbuf                                       //
//   ...                                                                //
//                                                                      //
// }                                                                    //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TGPasswdDialog.h"

#include "TError.h"
#include "TGFrame.h"
#include "TGButton.h"
#include "TGLabel.h"
#include "TGTextEntry.h"
#include "TGTextBuffer.h"
#include "TGString.h"
#include "TROOT.h"
#include "TVirtualX.h"


ClassImp(TGPasswdDialog);

////////////////////////////////////////////////////////////////////////////////
/// Create an editor in a dialog.

TGPasswdDialog::TGPasswdDialog(const char *prompt, char *pwdbuf, Int_t pwdlenmax,
                               UInt_t w, UInt_t h)
{
   fPwdBuf = pwdbuf;
   fPwdLenMax = pwdlenmax;

   const TGWindow *mainw = gClient->GetRoot();
   fDialog = new TGTransientFrame(mainw, mainw, w, h);
   fDialog->Connect("CloseWindow()", "TGPasswdDialog", this, "CloseWindow()");

   // Prompt
   fDialog->AddFrame(new TGLabel(fDialog, prompt),
                     new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 5, 5, 10, 5));

   // Passwd
   fPasswdText = new TGTextBuffer(40);
   fPasswd = new TGTextEntry(fDialog, fPasswdText);
   fPasswd->SetCursorPosition(0);
   fPasswd->Resize(300, fPasswd->GetDefaultHeight());
   fPasswd->SetEchoMode(TGTextEntry::kPassword);
   fPasswd->Connect("ReturnPressed()", "TGPasswdDialog", this, "ReturnPressed()");

   fDialog->AddFrame(fPasswd, new TGLayoutHints(kLHintsCenterY |
                                                kLHintsLeft | kLHintsExpandX,
                                                5, 5, 5, 5));
   // Ok button
   fOk = new TGTextButton(fDialog, "     &Ok     ");
   fOk->Connect("Clicked()", "TGPasswdDialog", this, "ReturnPressed()");
   fDialog->AddFrame(fOk, new TGLayoutHints(kLHintsBottom | kLHintsCenterX, 0, 0, 5, 5));
   // set window title and icon name
   fDialog->SetWindowName("Password dialog");
   fDialog->SetIconName("Password dialog");

   fDialog->MapSubwindows();

   Int_t width  = fDialog->GetDefaultWidth();
   Int_t height = fDialog->GetDefaultHeight();

   fDialog->Resize(width, height);

   fPasswd->SetFocus();
   // position relative to the parent window (which is the root window)
   Window_t wdum;
   int      ax, ay;
   Int_t    mw = ((TGFrame *) mainw)->GetWidth();
   Int_t    mh = ((TGFrame *) mainw)->GetHeight();

   gVirtualX->TranslateCoordinates(mainw->GetId(), mainw->GetId(),
                          (mw - width) >> 1, (mh - height) >> 1, ax, ay, wdum);
   fDialog->Move(ax, ay);
   fDialog->SetWMPosition(ax, ay);

   // make the message box non-resizable
   fDialog->SetWMSize(width, height);
   fDialog->SetWMSizeHints(width, height, width, height, 0, 0);

   // Now we wait for the user
   gROOT->SetInterrupt(kTRUE);

   fDialog->MapWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Delete log window.

TGPasswdDialog::~TGPasswdDialog()
{
   DoClose();
   delete fDialog;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle close button.

void TGPasswdDialog::DoClose()
{
   fDialog->SendCloseMessage();
}

////////////////////////////////////////////////////////////////////////////////
/// Called when closed via window manager action.

void TGPasswdDialog::CloseWindow()
{
   delete this;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle return

void TGPasswdDialog::ReturnPressed()
{
   if (fPwdBuf) {
      Int_t len = strlen(fPasswdText->GetString());
      len = (len < (fPwdLenMax - 1)) ? len : fPwdLenMax - 1;
      memcpy(fPwdBuf, fPasswdText->GetString(), len);
      fPwdBuf[len] = 0;
      fPasswdText->Clear();
   } else
      Error("ReturnPressed", "passwd buffer undefined");

   // We are done
   gROOT->SetInterrupt(kFALSE);

   // Close window
   fDialog->UnmapWindow();
}
