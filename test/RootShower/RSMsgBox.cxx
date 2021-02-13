// Author: Bertrand Bellenot   22/08/02

/*************************************************************************
 * Copyright (C) 1995-2002, Bertrand Bellenot.                           *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see the LICENSE file.                         *
 *************************************************************************/

#include <TSystem.h>
#include <TVirtualX.h>
#include <TRootHelpDialog.h>

#include "RSMsgBox.h"
#include "RootShower.h"
#include "RSHelpText.h"

////////////////////////////////////////////////////////////////////////////////

RootShowerMsgBox::RootShowerMsgBox(const TGWindow *p, const TGWindow *main,
                       UInt_t w, UInt_t h, UInt_t options) :
     TGTransientFrame(p, main, w, h, options)
{
   TGPicture *fIconPicture;
   TGIcon *fIcon;
   UInt_t wh1 = (UInt_t)(0.6 * h);
   UInt_t wh2 = h - wh1;

   fVFrame  = new TGVerticalFrame(this, w, wh1, 0);

   TString theLogoFilename = gProgPath;
   theLogoFilename.Append("/icons/mclogo01.xpm");
   fIconPicture = (TGPicture *)gClient->GetPicture(theLogoFilename);
   fIcon = new TGIcon(this, fIconPicture,
                      fIconPicture->GetWidth(),
                      fIconPicture->GetHeight());
   fLogoLayout = new TGLayoutHints(kLHintsCenterX, 0, 0, 0, 0);
   AddFrame(fIcon, fLogoLayout);

   fLabel1 = new TGLabel(fVFrame, "Appears that You have changed some parameters");
   fLabel2 = new TGLabel(fVFrame, "Would You wish to SAVE changes ?");

   fBly  = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX, 0, 0, 5, 5);
   fBfly = new TGLayoutHints(kLHintsTop | kLHintsRight| kLHintsExpandX, 0, 0, 5, 5);

   fVFrame->AddFrame(fLabel1,fBly);
   fVFrame->AddFrame(fLabel2,fBly);

//------------------------------------------------------------------------------
// OK Cancel Buttons in Horizontal frame
//------------------------------------------------------------------------------

   fHFrame  = new TGHorizontalFrame(this, w, wh2, 0);

   fOkButton = new TGTextButton(fHFrame, "&Yes", 1);
   fOkButton->Associate(this);

   fHelpButton = new TGTextButton(fHFrame, "&Help", 2);
   fHelpButton->Associate(this);

   fCancelButton = new TGTextButton(fHFrame, "&No", 3);
   fCancelButton->Associate(this);

   fL1 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,
                           2, 2, 2, 2);
   fL2 = new TGLayoutHints(kLHintsBottom | kLHintsLeft | kLHintsExpandX,
                           10, 10, 5, 10);

   fHFrame->AddFrame(fOkButton,     fL1);
   fHFrame->AddFrame(fHelpButton,   fL1);
   fHFrame->AddFrame(fCancelButton, fL1);

   fHFrame->Resize(w, fOkButton->GetDefaultHeight());

   AddFrame(fVFrame, fBfly);
   AddFrame(fHFrame, fL2);

   SetWindowName("Root's Monte Carlo Message Box");
   TGDimension size = GetDefaultSize();
   Resize(size);

   SetWMSize(size.fWidth, size.fHeight);
   SetWMSizeHints(size.fWidth, size.fHeight, size.fWidth, size.fHeight, 0, 0);
   SetMWMHints(kMWMDecorAll | kMWMDecorResizeH  | kMWMDecorMaximize |
               kMWMDecorMinimize | kMWMDecorMenu, kMWMFuncAll |
               kMWMFuncResize    | kMWMFuncMaximize | kMWMFuncMinimize,
               kMWMInputModeless);

   // position relative to the parent's window
   Window_t wdummy;
   Int_t ax, ay;
   gVirtualX->TranslateCoordinates(main->GetId(), GetParent()->GetId(),
                 (Int_t)(((TGFrame *) main)->GetWidth() - fWidth) >> 1,
                 (Int_t)(((TGFrame *) main)->GetHeight() - fHeight) >> 1,
                 ax, ay, wdummy);
   Move(ax, ay);

   MapSubwindows();
   MapWindow();

   fClient->WaitFor(this);
}

////////////////////////////////////////////////////////////////////////////////

RootShowerMsgBox::~RootShowerMsgBox()
{
   delete fLogoLayout;
   delete fCancelButton;
   delete fHelpButton;
   delete fOkButton;
   delete fVFrame;
   delete fHFrame;
   delete fLabel1;
   delete fLabel2;
   delete fBly;
}

////////////////////////////////////////////////////////////////////////////////
/// Close dialog in response to window manager close.

void RootShowerMsgBox::CloseWindow()
{
   DeleteWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Process messages sent to this dialog.

Bool_t RootShowerMsgBox::ProcessMessage(Long_t msg, Long_t parm1, Long_t /*parm2*/)
{
   TRootHelpDialog* hd;

   switch (GET_MSG(msg)) {

      case kC_COMMAND:
         switch (GET_SUBMSG(msg)) {
            case kCM_BUTTON:
               switch(parm1) {
                  case 1:
                     gRootShower->SetOk();
                     CloseWindow();
                     break;
                  case 2:
                     hd = new TRootHelpDialog(this, "Help on Message Box", 590, 200);
                     hd->SetText(gRSHelpMsgBox);
                     hd->Popup();
                     fClient->WaitFor(hd);
                     break;
                  case 3:                  // Cancel button
                     gRootShower->SetOk(false);
                     CloseWindow();
                     break;
                  default:
                     break;
               }
               break;
            default:
               break;
         }
         break;
   }
   return kTRUE;
}

