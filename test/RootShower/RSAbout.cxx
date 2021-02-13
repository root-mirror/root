// Author: Bertrand Bellenot   22/08/02

/*************************************************************************
 * Copyright (C) 1995-2002, Bertrand Bellenot.                           *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see the LICENSE file.                         *
 *************************************************************************/

#include <TSystem.h>
#include <TVirtualX.h>
#include <TROOT.h>
#include <TRootHelpDialog.h>

#include "RSAbout.h"
#include "RSVersion.h"


////////////////////////////////////////////////////////////////////////////////

RootShowerAbout::RootShowerAbout(const TGWindow *p, const TGWindow *main,
                       UInt_t w, UInt_t h, UInt_t options) :
     TGTransientFrame(p, main, w, h, options)
{
   Int_t iday,imonth,iyear;
   Char_t message1[80];
   static const Char_t *months[] = {"January","February", "March","April","May","June","July",
                                    "August", "September","October","November","December"};
   TGPicture *fIconPicture;
   TGIcon *fIcon;
   UInt_t wh1 = (UInt_t)(0.6 * h);
   UInt_t wh2 = h - wh1;

   const Char_t *root_version = gROOT->GetVersion();
   Int_t idatqq = gROOT->GetVersionDate();
   iday   = idatqq%100;
   imonth = (idatqq/100)%100;
   iyear  = (idatqq/10000);
   Char_t *root_date = Form("%s %d %4d",months[imonth-1],iday,iyear);
   fVFrame  = new TGVerticalFrame(this, w, wh1, 0);

   TString theLogoFilename = gProgPath;
   theLogoFilename.Append("/icons/mclogo01.xpm");

   fIconPicture = (TGPicture *)gClient->GetPicture(theLogoFilename);
   fIcon = new TGIcon(this, fIconPicture,
                      fIconPicture->GetWidth(),
                      fIconPicture->GetHeight());
   fLogoLayout = new TGLayoutHints(kLHintsCenterX, 0, 0, 0, 0);
   AddFrame(fIcon, fLogoLayout);

   sprintf(message1,"   ROOT Shower Monte Carlo / Event Display v %s   ",
           ROOTSHOWER_RELEASE);
   fLabel1 = new TGLabel(fVFrame, message1);
   sprintf(message1,"   Compiled with Root version %s, release date : %s   ",
           root_version, root_date);
   fLabel2 = new TGLabel(fVFrame, message1);
   fLabel4 = new TGLabel(fVFrame, "   (c) Bertrand Bellenot 1994 - 2008   ");

   fBly  = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX, 0, 0, 5, 5);
   fBfly = new TGLayoutHints(kLHintsTop | kLHintsRight| kLHintsExpandX, 0, 0, 5, 5);

   fVFrame->AddFrame(fLabel1,fBly);
   fVFrame->AddFrame(fLabel2,fBly);
   fVFrame->AddFrame(fLabel4,fBly);

   //------------------------------------------------------------------------------
   // OK Cancel Buttons in Horizontal frame
   //------------------------------------------------------------------------------

   fHFrame  = new TGHorizontalFrame(this, w, wh2, 0);

   fOkButton = new TGTextButton(fHFrame, "&Ok", 1);
   fOkButton->Resize(100, fOkButton->GetDefaultHeight());
   fOkButton->Associate(this);

   fL1 = new TGLayoutHints(kLHintsCenterX | kLHintsExpandX, 2, 2, 2, 2);
   fL2 = new TGLayoutHints(kLHintsBottom | kLHintsCenterX | kLHintsExpandX, 150, 150, 2, 10);

   fHFrame->AddFrame(fOkButton,     fL1);

   fHFrame->Resize(100, fOkButton->GetDefaultHeight());

   AddFrame(fVFrame, fBfly);
   AddFrame(fHFrame, fL2);

   SetWindowName("About RootShower...");
   TGDimension size = GetDefaultSize();
   Resize(size);

   SetWMSize(size.fWidth, size.fHeight);
   SetWMSizeHints(size.fWidth, size.fHeight, size.fWidth, size.fHeight, 0, 0);
   SetMWMHints(kMWMDecorAll | kMWMDecorResizeH  | kMWMDecorMaximize |
               kMWMDecorMinimize | kMWMDecorMenu, kMWMFuncAll |
               kMWMFuncResize    | kMWMFuncMaximize | kMWMFuncMinimize,
               kMWMInputModeless);

   // position relative to the parent's window
   Int_t      ax, ay;
   Window_t wdummy;
   gVirtualX->TranslateCoordinates(main->GetId(), GetParent()->GetId(),
                                 (Int_t)(((TGFrame *) main)->GetWidth() - fWidth) >> 1,
                                 (Int_t)(((TGFrame *) main)->GetHeight() - fHeight) >> 1,
                                 ax, ay, wdummy);
   if (ax < 0) ax = 10;
   if (ay < 0) ay = 10;
   Move(ax, ay);
   SetWMPosition(ax, ay);

   MapSubwindows();
   MapWindow();

   fClient->WaitFor(this);
}

////////////////////////////////////////////////////////////////////////////////

RootShowerAbout::~RootShowerAbout()
{
   delete fOkButton;
   delete fVFrame;
   delete fHFrame;
   delete fLabel1;
   delete fLabel2;
   delete fLabel4;
   delete fBly;
}

////////////////////////////////////////////////////////////////////////////////
/// Close dialog in response to window manager close.

void RootShowerAbout::CloseWindow()
{
   DeleteWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Process messages sent to this dialog.

Bool_t RootShowerAbout::ProcessMessage(Long_t msg, Long_t /*parm1*/,
                                       Long_t /*parm2*/)
{
   switch (GET_MSG(msg)) {
      case kC_COMMAND:
         switch (GET_SUBMSG(msg)) {
            case kCM_BUTTON:
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

