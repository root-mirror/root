// Author: Bertrand Bellenot   22/08/02

/*************************************************************************
 * Copyright (C) 1995-2002, Bertrand Bellenot.                           *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see the LICENSE file.                         *
 *************************************************************************/

#include <cstdlib>
#include <TGLabel.h>
#include <TVirtualX.h>
#include <TRootHelpDialog.h>
#include <TGMsgBox.h>

#include "SettingsDlg.h"
#include "RootShower.h"
#include "constants.h"
#include "RSHelpText.h"

// definition of structure used to retrieve
// PDG number and particle name into the list
// of available primaries

typedef struct {
    Int_t       pdg_code;
    const char *pdg_name;
} str_choice_def;

str_choice_def choice_def[] = {
      {    22,  "gamma"       },
      {    11,  "e-"          },
      {   -11,  "e+"          },
      {    13,  "mu-"         },
      {   -13,  "mu+"         },
      {    15,  "tau-"        },
      {   -15,  "tau+"        },
      {   111,  "pi0"         },
      {   211,  "pi+"         },
      {  -211,  "pi-"         },
      {   221,  "Eta"         },
      {   321,  "K+"          },
      {  -321,  "K-"          },
      {   130,  "K(L)0"       },
      {   310,  "K(S)0"       },
      {   113,  "rho(770)0"   },
      {   213,  "rho(770)+"   },
      {  -213,  "rho(770)-"   },
      {   223,  "omega(782)0" },
      {   333,  "phi(1020)0"  },
      {   443,  "J/psi(1S)0"  },
      {   511,  "B0"          },
      {   513,  "B*0"         },
      {   521,  "B+"          },
      {  -521,  "B-"          },
      {   523,  "B*+"         },
      {  -523,  "B*-"         },
      {   531,  "B(s)0"       },
      {   411,  "D+"          },
      {  -411,  "D-"          },
      {   421,  "D0"          },
      {   431,  "D(s)+"       },
      {  -431,  "D(s)-"       },
      {   433,  "D(s)*+"      },
      {  -433,  "D(s)*-"      },
      {    24,  "W+"          },
      {   -24,  "W-"          },
      {    23,  "Z0"          },
      {     0,    0           }
};

//_________________________________________________
// SettingsDialog
//
// SettingsDialog is a dialog accessing the RootShowerhower parameters.

enum RootShowerSettingsTypes {
    Id1,
    Id2,
    Id3,
    Id4,
    Id5
};


////////////////////////////////////////////////////////////////////////////////
/// Create a dialog window. A dialog window pops up with respect to its
/// "main" window.

SettingsDialog::SettingsDialog(const TGWindow *p, const TGWindow *main, UInt_t w,
                       UInt_t h, UInt_t options)
    : TGTransientFrame(p, main, w, h, options)
{
   Int_t i;
   Char_t tmp[20];
   UInt_t wh1 = (UInt_t)(0.6 * h);
   UInt_t wh2 = h - wh1;

   fFrame1  = new TGHorizontalFrame(this, w, wh2, 0);

   fOkButton = new TGTextButton(fFrame1, "     &Ok     ", 1);
   fOkButton->Associate(this);
   fCancelButton = new TGTextButton(fFrame1, "   &Cancel   ", 2);
   fCancelButton->Associate(this);
   fHelpButton = new TGTextButton(fFrame1, "    &Help    ", 3);
   fHelpButton->Associate(this);

   fL1 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,
                           2, 2, 2, 2);
   fL2 = new TGLayoutHints(kLHintsBottom | kLHintsRight | kLHintsExpandX, 2, 2, 5, 1);

   fFrame1->AddFrame(fOkButton, fL1);
   fFrame1->AddFrame(fHelpButton, fL1);
   fFrame1->AddFrame(fCancelButton, fL1);

   fFrame1->Resize(150, fOkButton->GetDefaultHeight());
   AddFrame(fFrame1, fL2);

   //--------- create Tab widget and some composite frames for Tab testing
   fTab = new TGTab(this, 300, 300);
   fL3 = new TGLayoutHints(kLHintsTop | kLHintsLeft, 10, 10, 10, 10);

   TGCompositeFrame *tf = fTab->AddTab("Physics settings");
   tf->SetLayoutManager(new TGHorizontalLayout(tf));

   fF3 = new TGGroupFrame(tf, "Particle", kVerticalFrame);
   tf->AddFrame(fF3, fL3);
   fF3->SetLayoutManager(new TGMatrixLayout(fF3, 0, 1, 10));

   fF3->AddFrame(fListBox = new TGListBox(fF3, 89));

   for (i = 0; choice_def[i].pdg_name; i++) {
      fListBox->AddEntry(choice_def[i].pdg_name, i);
   }
   fFirstEntry = 0;
   fLastEntry  = 30;
   fListBox->Resize(120, 110);

   fF3->Resize(fF3->GetDefaultSize());

   fF4 = new TGGroupFrame(tf, "E0 / B", kVerticalFrame);
   fF4->SetTitlePos(TGGroupFrame::kRight); // right aligned
   tf->AddFrame(fF4, fL3);
   fF4->SetLayoutManager(new TGMatrixLayout(fF4, 0, 2, 10));

   fF4->AddFrame(new TGLabel(fF4, "E0 [GeV]"));
   fF4->AddFrame(fTxt4 = new TGTextEntry(fF4, new TGTextBuffer(100), Id4));
   fF4->AddFrame(new TGLabel(fF4, "B [kGauss]"));
   fF4->AddFrame(fTxt5 = new TGTextEntry(fF4, new TGTextBuffer(100), Id5));
   fTxt4->Associate(this);
   fTxt5->Associate(this);
   fTxt4->Resize(65, fTxt4->GetDefaultHeight());
   fTxt5->Resize(65, fTxt5->GetDefaultHeight());
   sprintf(tmp,"%1.4f",gRootShower->fE0);
   fTxt4->SetText(tmp);
   sprintf(tmp,"%1.4f",gRootShower->fB);
   fTxt5->SetText(tmp);

   TGLayoutHints *fL5 = new TGLayoutHints(kLHintsBottom | kLHintsExpandX |
                                          kLHintsExpandY, 2, 2, 5, 1);
   AddFrame(fTab, fL5);

   MapSubwindows();
   Resize(GetDefaultSize());

   fF4->Resize(fF4->GetDefaultWidth(),fF3->GetDefaultHeight());

   for (i = 0; choice_def[i].pdg_name; i++) {
      if (gRootShower->fFirstParticle == choice_def[i].pdg_code) {
         fListBox->Select(i);
         fListBox->GetScrollBar()->SetPosition(i-2);
         break;
      }
   }
   fListBox->MapSubwindows();
   fListBox->Layout();
   // position relative to the parent's window
   Window_t wdum;
   Int_t ax, ay;
   gVirtualX->TranslateCoordinates(main->GetId(), GetParent()->GetId(),
                          (Int_t)(((TGFrame *) main)->GetWidth() - fWidth) >> 1,
                          (Int_t)(((TGFrame *) main)->GetHeight() - fHeight) >> 1,
                          ax, ay, wdum);
   Move(ax, ay);

   SetWindowName("Shower Settings");

   // make the message box non-resizable
   UInt_t width  = GetDefaultWidth();
   UInt_t height = GetDefaultHeight();
   SetWMSize(width, height);
   SetWMSizeHints(width, height, width, height, 0, 0);

   SetMWMHints(kMWMDecorAll | kMWMDecorResizeH  | kMWMDecorMaximize |
               kMWMDecorMinimize | kMWMDecorMenu,
               kMWMFuncAll  | kMWMFuncResize    | kMWMFuncMaximize |
               kMWMFuncMinimize, kMWMInputModeless);
   MapWindow();
   fClient->WaitFor(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Delete test dialog widgets.

SettingsDialog::~SettingsDialog()
{
   delete fOkButton;
   delete fCancelButton;
   delete fHelpButton;
   delete fFrame1;
   delete fListBox;
   delete fF3;
   delete fF4;
   delete fTxt4;
   delete fTxt5;
   delete fTab;
   delete fL3;
   delete fL2;
   delete fL1;
}

////////////////////////////////////////////////////////////////////////////////
/// Called when window is closed via the window manager.

void SettingsDialog::CloseWindow()
{
   DeleteWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Process messages coming from widgets associated with the dialog.

Bool_t SettingsDialog::ProcessMessage(Long_t msg, Long_t parm1, Long_t)
{
   Int_t Selection;
   Int_t retval;
   TRootHelpDialog* hd;

   const Char_t *buf_tmp;

   switch (GET_MSG(msg)) {
      case kC_COMMAND:

         switch (GET_SUBMSG(msg)) {
            case kCM_BUTTON:
               switch (parm1) {
                  case 1:
                     Selection = fListBox->GetSelected();
                     if (Selection > 37) {
                        new TGMsgBox(fClient->GetRoot(), this, "Particle selection",
                                     "This particle is not implemented yet !",
                                     kMBIconExclamation, kMBOk, &retval);
                        fListBox->Select(2);
                        break;
                     }
                     gRootShower->fFirstParticle = choice_def[Selection].pdg_code;
                     buf_tmp = fTxt4->GetBuffer()->GetString();
                     gRootShower->fE0 = atof(buf_tmp);
                     buf_tmp = fTxt5->GetBuffer()->GetString();
                     gRootShower->fB = atof(buf_tmp);
                     gRootShower->SettingsModified();
                  case 2:
                     CloseWindow();
                     break;
                  case 3:
                     hd = new TRootHelpDialog(this, "Help on Settings Dialog", 560, 400);
                     hd->SetText(gSettingsHelp);
                     hd->Popup();
                     fClient->WaitFor(hd);
                     break;
                  default:
                     break;
               }
               break;
            case kCM_TAB:
               break;
            default:
               break;
         }
         break;
      case kC_TEXTENTRY:
         switch (GET_SUBMSG(msg)) {
            case kTE_ENTER:
               switch (parm1) {
                  case Id4:
                     fTxt5->SetFocus();
                     break;
                  case Id5:
                     fTxt4->SetFocus();
                     break;
               }
               break;
            case kTE_TAB:
               switch (parm1) {
                  case Id4:
                     fTxt5->SetFocus();
                     break;
                  case Id5:
                     fTxt4->SetFocus();
                     break;
               }
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
