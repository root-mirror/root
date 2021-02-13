//----------------------------------------------------
//
//    A small panel destined to provide graphical interface
//    for "Canvas Save" action
//
//    Author: Dmitry Vasiliev (LNS, Catania)
//
//----------------------------------------------------
//
#include "canvsave.h"
#include <TCollection.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <snprintf.h>

ClassImp(CanvSave);


CanvSave::CanvSave(const TGWindow *p, const TGWindow *main, UInt_t w,
                   UInt_t h, UInt_t options) :
      TGTransientFrame(p, main, w, h, options)
{
   /*
   //--------------------------------------------------------
   //
   //     Constructor for the class CanvSave
   //
   //-------------------------------------------------------
   //
   //
   //                   -------------------------
   //                   |                       |
   //                   |-----------------------|
   //                   |         fFrame        |
   //                   |-----------------------|
   //                   |                       |
   //                   -------------------------
   //                   /         |              \
   //                  /          |               \
   //          ------------    ----------       ----------------
   //          |          |    |        |       |      |       |
   //          |   fF1    |    |  fF2   |       |     fF3      |
   //          |          |    |        |       |      |       |
   //          ------------    ----------       ----------------
   //               |              |             /            \
   //          ------------    ----------    ----------    -----------
   //          |          |    |        |    |        |    |         |
   //          |  fCombo  |    |  fText |    | fSave  |    | fCancel |
   //          |          |    |        |    |        |    |         |
   //          ------------    ----------    ----------    -----------
   //
   //
   /////////////////////////////////////////////////////////////////////////////

   */

   fFrame = new TGCompositeFrame(this, 400, 200, kVerticalFrame);
   AddFrame(fFrame, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));

   fF1 = new TGCompositeFrame(fFrame, 400, 50, kHorizontalFrame);
   fF1->SetLayoutManager(new TGMatrixLayout(fF1, 0, 2, 10));
   fFrame->AddFrame(fF1, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));

   fF2 = new TGCompositeFrame(fFrame, 400, 50, kHorizontalFrame);
   fF2->SetLayoutManager(new TGMatrixLayout(fF2, 0, 2, 10));
   fFrame->AddFrame(fF2, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));

   fF3 = new TGCompositeFrame(fFrame, 400, 100, kHorizontalFrame);
   fFrame->AddFrame(fF3, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));

   fF1->AddFrame(new TGLabel(fF1, new TGHotString("Save as")));
   fCombo = new TGComboBox(fF1, 10);
   fCombo->Associate(this);
   char tmp[20];
   snprintf(tmp,20, "%s", ".ps");
   fCombo->AddEntry(tmp, 1);
   snprintf(tmp,20, "%s", ".gif");
   fCombo->AddEntry(tmp, 2);
   fCombo->Resize(100, 20);
   fCombo->Select(1);
   fF1->AddFrame(fCombo);

   fF2->AddFrame(new TGLabel(fF2, new TGHotString("Enter file name")));
   fBuf = new TGTextBuffer(100);
   fText = new TGTextEntry(fF2, fBuf, -1);
   fText->Resize(100, 20);
   fBuf->AddText(0, "histo.ps");
   fF2->AddFrame(fText);

   fSave = new TGTextButton(fF3, " Save  ", 100);
   fSave->Associate(this);
   fSave->Resize(60, 20);
   fF3->AddFrame(fSave, new TGLayoutHints(kLHintsLeft | kLHintsTop, 60, 10, 2, 2));

   fCancel = new TGTextButton(fF3, "Cancel", 101);
   fCancel->Associate(this);
   fCancel->Resize(60, 20);
   fF3->AddFrame(fCancel, new TGLayoutHints(kLHintsLeft | kLHintsTop, 10, 5, 2, 2));

   MapSubwindows();

   Window_t wdum;
   Int_t ax, ay;

   //printf("GetWidth = %ld\tfWidth = %ld\n",((TGFrame*)main)->GetWidth(),fWidth);

   gVirtualX->TranslateCoordinates(main->GetId(), GetParent()->GetId(),
                                   (((TGFrame *) main)->GetWidth() - fWidth) >> 1,
                                   (((TGFrame *) main)->GetHeight() - fHeight) >> 1,
                                   ax, ay, wdum);
   Move(ax, ay);

   SetWindowName("SAVE PANEL");

   Resize(GetDefaultSize());
   MapWindow();
   fClient->WaitFor(this); //otherwise the current directory in the file view
   //window is not always updated when a picture is saved
}

CanvSave::~CanvSave()
{
   //-----------------------------------------------------
   //
   //      Destructor for the class CanvSave
   //
   //----------------------------------------------------

   delete fCancel;
   delete fSave;
   delete fText;
   delete fCombo;
   delete fF3;
   delete fF2;
   delete fF1;
   delete fFrame;
}

void CanvSave::CloseWindow()
{
   //---------------------------------------------------
   //
   //                CloseWindow()
   //
   //    Closes the window "SAVE PANEL" (apparently)
   //
   //---------------------------------------------------

   delete this;
}

Bool_t CanvSave::ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2)
{
   //------------------------------------------------------------
   //
   //     ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2)
   //
   //   Processes information from GUI items of the panel
   //
   //------------------------------------------------------------

   TIter it(gROOT->GetListOfCanvases());
   TCanvas *c;
   while ((c = (TCanvas*) it()))  {
      if (!strcmp("canvasA", c->GetName())) break;
   }

   switch (GET_MSG(msg)) {
      case kC_COMMAND:

         switch (GET_SUBMSG(msg)) {

            case kCM_BUTTON:

               switch (parm1) {
                  case 100:

                     //'Save' button is clicked

                     if (!strcmp("canvasA", c->GetName())) {
                        c->cd();
                        c->Print(fBuf->GetString());
                     }
                     CloseWindow();
                     break;

                  case 101:
                     CloseWindow();
                     break;

                  default:
                     break;
               }

            case kCM_COMBOBOX:

               // Process Combo box

               switch (parm1) {
                  case 10:
                     if (parm2 == 1) {
                        fBuf->Clear();
                        fBuf->AddText(0, "hist.ps");
                        fClient->NeedRedraw(fText);
                     }
                     if (parm2 == 2) {
                        fBuf->Clear();
                        fBuf->AddText(0, "histo.gif");
                        fClient->NeedRedraw(fText);
                     }
                     break;

                  default:
                     break;
               }

            default:
               break;
         }

      default:
         break;
   }

   return kTRUE;
}
