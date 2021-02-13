// @(#)root/ged:$Id$
// Author: Denis Favre-Miville   08/09/05

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  TStyleDialog                                                        //
//                                                                      //
//  This small class is useful to ask the user for a name and a title,  //
//       in order to rename a style, create a new style or import a     //
//       style from a canvas.                                           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TStyleDialog.h"
#include "TStyleManager.h"

#include <TCanvas.h>
#include <TGButton.h>
#include <TGLabel.h>
#include <TGLayout.h>
#include <TGTextEntry.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TVirtualMutex.h>

#include <cstdlib>

ClassImp(TStyleDialog);

enum EStyleDialogWid {
   kName,
   kTitle,
   kButOK,
   kButCancel
};

////////////////////////////////////////////////////////////////////////////////
///  Constructor. Create the dialog window and draw it centered over the
/// main window 'mf'. A pointer to the style to copy or rename is given
/// by 'cur' and the parameter 'mode' specify the mode:
///    1 = copy  |  2 = rename  |  3 = import from canvas.

TStyleDialog::TStyleDialog(TStyleManager *sm, TStyle *cur, Int_t mode,
                              TVirtualPad *currentPad)
                     : TGTransientFrame(0, sm)
{
   fStyleManager = sm;

   // Create the main frame.
   SetCleanup(kNoCleanup);
   fCurStyle = cur;
   fMode = mode;
   fCurPad = currentPad;

   switch (fMode) {
      case 1:  SetWindowName("Create a New Style");            break;
      case 2:  SetWindowName("Rename the Selected Style");     break;
      case 3:  SetWindowName("Import a New Style from Canvas");
   }

   // Create the trash lists to have an effective deletion of every object.
   fTrashListLayout = new TList();
   fTrashListFrame = new TList();

   // Create the layouts and add them to the layout trash list.
   TGLayoutHints *layoutNameLabel = new TGLayoutHints(kLHintsNormal, 0, 70, 3);
   fTrashListLayout->Add(layoutNameLabel);
   TGLayoutHints *layoutTitleLabel = new TGLayoutHints(kLHintsNormal, 0, 39, 3);
   fTrashListLayout->Add(layoutTitleLabel);
   TGLayoutHints *layoutWarningLabel = new TGLayoutHints(kLHintsExpandX);
   fTrashListLayout->Add(layoutWarningLabel);
   TGLayoutHints *layoutOKButton = new TGLayoutHints(kLHintsExpandX, 0, 5);
   fTrashListLayout->Add(layoutOKButton);
   TGLayoutHints *layoutCancelButton = new TGLayoutHints(kLHintsExpandX, 5);
   fTrashListLayout->Add(layoutCancelButton);
   TGLayoutHints *layoutH1 = new TGLayoutHints(kLHintsExpandX, 10, 10, 10, 5);
   fTrashListLayout->Add(layoutH1);
   TGLayoutHints *layoutH2 = new TGLayoutHints(kLHintsExpandX, 10, 10, 5, 5);
   fTrashListLayout->Add(layoutH2);
   TGLayoutHints *layoutH4 = new TGLayoutHints(kLHintsExpandX, 10, 10, 5, 10);
   fTrashListLayout->Add(layoutH4);

   // Create and place the widgets in the main window.
   // Every frame created here must be added to the frame trash list.
   TGHorizontalFrame *h1 = new TGHorizontalFrame(this);
   fTrashListFrame->Add(h1);
   fNameLabel = new TGLabel(h1, "Name:");
   h1->AddFrame(fNameLabel, layoutNameLabel);

   if (fMode == 1) {
      TString newName;
      newName.Form("%s_1", fCurStyle->GetName());
      fName = new TGTextEntry(h1, newName.Data(), kName);
   } else if (fMode == 2) {
      // The names of the 5 basics styles can not be modified.
      fName = new TGTextEntry(h1, fCurStyle->GetName(), kName);
      if ((!strcmp(fName->GetText(), "Default"))
       || (!strcmp(fName->GetText(), "Plain"  ))
       || (!strcmp(fName->GetText(), "Bold"   ))
       || (!strcmp(fName->GetText(), "Video"  ))
       || (!strcmp(fName->GetText(), "Pub"    ))) fName->SetEnabled(kFALSE);
   } else
      fName = new TGTextEntry(h1, "Imported_Style", kName);
   fName->Associate(this);
   fName->Resize(200, 22);
   h1->AddFrame(fName);
   AddFrame(h1, layoutH1);

   TGHorizontalFrame *h2 = new TGHorizontalFrame(this);
   fTrashListFrame->Add(h2);
   fTitleLabel = new TGLabel(h2, "Description:");
   h2->AddFrame(fTitleLabel, layoutTitleLabel);
   switch (fMode) {
      case 1:
      case 2:
         fTitle = new TGTextEntry(h2, fCurStyle->GetTitle(), kTitle);
         break;
      case 3: {
         TString newTitle("Imported from canvas ");
         if (fCurPad->GetCanvas())
            newTitle += fCurPad->GetCanvas()->GetName();
         fTitle = new TGTextEntry(h2, newTitle.Data(), kTitle);
      }
   }
   fTitle->Associate(this);
   fTitle->Resize(200, 22);
   h2->AddFrame(fTitle);
   fTitle->Associate(h2);
   AddFrame(h2, layoutH2);

   TGHorizontalFrame *h3 = new TGHorizontalFrame(this);
   fTrashListFrame->Add(h3);
   fWarnLabel = new TGLabel(h3);
   Pixel_t red;
   gClient->GetColorByName("#FF0000", red);
   fWarnLabel->SetTextColor(red, kFALSE);
   fWarnLabel->Resize(200, 22);
   h3->AddFrame(fWarnLabel, layoutWarningLabel);
   AddFrame(h3, layoutH2);

   TGHorizontalFrame *h4 = new TGHorizontalFrame(this);
   fTrashListFrame->Add(h4);
   fOK = new TGTextButton(h4, "&OK", kButOK);
   fOK->Associate(this);
   h4->AddFrame(fOK, layoutOKButton);
   fOK->Associate(h4);
   fCancel = new TGTextButton(h4, "&Cancel", kButCancel);
   fCancel->Associate(this);
   h4->AddFrame(fCancel, layoutCancelButton);
   fCancel->Associate(h4);
   AddFrame(h4, layoutH4);

   // Refresh the warning message.
   DoUpdate();

   Resize();
   CenterOnParent();
   MapSubwindows();
   Int_t w = GetDefaultWidth();
   Int_t h = GetDefaultHeight();
   SetWMSizeHints(w, h, w, h, 0, 0);
   MapWindow();

   switch (fMode) {
      case 1:
         fOK->SetToolTipText("Create this new style");
         fCancel->SetToolTipText("Cancel the creation ");
         break;
      case 2:
         fOK->SetToolTipText("Rename the selected style");
         fCancel->SetToolTipText("Cancel the rename ");
         break;
      case 3:
         fOK->SetToolTipText("Import this new style from the canvas");
         fCancel->SetToolTipText("Cancel the import");
         break;
   }

   Connect("CloseWindow()", "TStyleDialog", this, "DoCloseWindow()");
   fName->Connect("TextChanged(const char *)", "TStyleDialog", this, "DoUpdate()");
   fOK->Connect("Clicked()", "TStyleDialog", this, "DoOK()");
   fCancel->Connect("Clicked()", "TStyleDialog", this, "DoCancel()");

   gClient->WaitFor(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TStyleDialog::~TStyleDialog()
{
   Disconnect("DoCloseWindow()");
   fName->Disconnect("TextChanged(const char *)");
   fOK->Disconnect("Clicked()");
   fCancel->Disconnect("Clicked()");

   delete fName;
   delete fNameLabel;
   delete fTitle;
   delete fTitleLabel;
   delete fWarnLabel;
   delete fOK;
   delete fCancel;

   TObject *obj1;
   TObject *obj2;

   obj1 = fTrashListFrame->First();
   while (obj1) {
      obj2 = fTrashListFrame->After(obj1);
      fTrashListFrame->Remove(obj1);
      delete obj1;
      obj1 = obj2;
   }
   delete fTrashListFrame;

   obj1 = fTrashListLayout->First();
   while (obj1) {
      obj2 = fTrashListLayout->After(obj1);
      fTrashListLayout->Remove(obj1);
      delete obj1;
      obj1 = obj2;
   }
   delete fTrashListLayout;
}

////////////////////////////////////////////////////////////////////////////////
///  Slot called when the Cancel button is clicked. Close the window
/// without saving submitted changes.

void TStyleDialog::DoCancel()
{
   fStyleManager->SetLastChoice(kFALSE);

   SendCloseMessage();
}

////////////////////////////////////////////////////////////////////////////////
///  Slot called when the window is closed via the window manager.
/// Close the window without saving submitted changes.

void TStyleDialog::DoCloseWindow()
{
   delete this;
}

////////////////////////////////////////////////////////////////////////////////
///  Slot called when the OK button is clicked. Rename or create the style
/// before closing the window.

void TStyleDialog::DoOK()
{
   if (fMode == 2) {
      // Update the name and the title of the style.
      fCurStyle->SetName(fName->GetText());
      fCurStyle->SetTitle(fTitle->GetText());
   } else {
      // Create a new style (copy of fCurStyle), with the given name and title.
      TStyle *tmpStyle = new TStyle(*fCurStyle);
      tmpStyle->SetName(fName->GetText());
      tmpStyle->SetTitle(fTitle->GetText());
      {
         R__LOCKGUARD(gROOTMutex);
         gROOT->GetListOfStyles()->Add(tmpStyle);
      }
      if (fMode == 3) {
         // Import the properties of the canvas.
         TStyle *tmp = gStyle;
         gStyle = tmpStyle;
         gStyle->SetIsReading(kFALSE);
         if (fCurPad->GetCanvas())
            fCurPad->GetCanvas()->UseCurrentStyle();
         gStyle->SetIsReading(kTRUE);
         gStyle = tmp;
      }
   }

   fStyleManager->SetLastChoice(kTRUE);

   SendCloseMessage();
}

////////////////////////////////////////////////////////////////////////////////
///  Slot called every time the name is changed. Provide some protection
/// to avoid letting the user use an empty name or an already used one.
///  A warning message can be shown and the OK button disabled.

void TStyleDialog::DoUpdate()
{
   if (!strlen(fName->GetText())) {
      fWarnLabel->SetText("That name is empty");
      fOK->SetEnabled(kFALSE);
      return;
   }

   if (strstr(fName->GetText(), " ") != 0) {
      fWarnLabel->SetText("That name contains some spaces");
      fOK->SetEnabled(kFALSE);
      return;
   }

   switch (fMode) {
      case 1:
      case 3:
         if (gROOT->GetStyle(fName->GetText())) {
            fWarnLabel->SetText("That name is already used by another style.");
            fOK->SetEnabled(kFALSE);
            return;
         }
         break;
      case 2:
         TStyle *tmp = gROOT->GetStyle(fName->GetText());
         if (tmp && (tmp != fCurStyle)) {
            fWarnLabel->SetText("That name is already used by another style.");
            fOK->SetEnabled(kFALSE);
            return;
         }
   }

   fWarnLabel->SetText("");
   fOK->SetEnabled(kTRUE);
}
