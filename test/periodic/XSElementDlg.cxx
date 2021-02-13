/*
 * $Header$
 * $Log$
 */

#include <ctype.h>
#include <stdlib.h>

#include <TGTab.h>
#include <TGFrame.h>
#include <TGButton.h>
#include <TGWindow.h>
#include <TGButton.h>
#include <TGMsgBox.h>
#include <TVirtualX.h>

#include "XSElementList.h"
#include "XSPeriodicTable.h"
#include "XSElementDlg.h"

//ClassImp(XSElementDlg);

enum    XSElementDlgEnum {
      ELEM_OK=900,
      ELEM_CANCEL
   };

static  Int_t   activeTab=0;            // Remember last active Tab
static   Int_t   LastWinX = -1;
static   Int_t   LastWinY = -1;


/* ----- XSElementDlg ----- */
XSElementDlg::XSElementDlg( const TGWindow *p, const TGWindow *main,
      UInt_t *retZ, UInt_t w, UInt_t h)
   : TGTransientFrame(p,main,w,h)
{
   selectedZ = retZ;   // Save Address

   lHints = new TGLayoutHints( kLHintsTop | kLHintsLeft |
         kLHintsExpandX | kLHintsExpandY,
         2, 2, 2, 2);

   buttonFrame = new TGHorizontalFrame(this, 120, 20, kFixedWidth);
   okButton = new TGTextButton(buttonFrame, "&Ok", ELEM_OK);
   okButton->Associate(this);
   closeButton = new TGTextButton(buttonFrame, "&Close", ELEM_CANCEL);
   closeButton->Associate(this);

   buttonLayoutHints = new TGLayoutHints(
            kLHintsTop | kLHintsLeft | kLHintsExpandX,
            2, 2, 5, 1);
   buttonFrame->AddFrame(okButton, buttonLayoutHints);
   buttonFrame->AddFrame(closeButton, buttonLayoutHints);

   frameLayoutHints = new TGLayoutHints(kLHintsBottom | kLHintsRight,
         2, 2, 2, 2);
   AddFrame(buttonFrame, frameLayoutHints);

   // ---------- Add Periodic table -------------
   TGCompositeFrame *tf;
   tabMenu = new TGTab(this, w, h);

   tf = tabMenu->AddTab("Periodic Table");

   pTable = new XSPeriodicTable(this,tf,w,h);
   pTable->SelectZ(*selectedZ);
   tf->AddFrame(pTable, lHints);

   // ------------- Table Sorted by Name ------------
   tf = tabMenu->AddTab("Name");

   nameFrame = new TGCompositeFrame(tf,w,h, kVerticalFrame);
   nameListBox = new XSElementList(nameFrame,XSEL_SORTBY_NAME);
   nameListBox->SelectZ(*selectedZ);
   nameFrame->AddFrame(nameListBox,lHints);

   tf->AddFrame(nameFrame,lHints);

   // ------------- Table Sorted by Mnemonic -------------
   tf = tabMenu->AddTab("Mnemonic");

   mnemonicFrame = new TGCompositeFrame(tf,w,h,kVerticalFrame);
   mnemonicListBox = new XSElementList(mnemonicFrame,XSEL_SORTBY_MNEMONIC);
   mnemonicListBox->SelectZ(*selectedZ);
   mnemonicFrame->AddFrame(mnemonicListBox, lHints);

   tf->AddFrame(mnemonicFrame, lHints);

   // ------------- Table Sorted by Charge -------------
   tf = tabMenu->AddTab("Z (Charge)");

   zFrame = new TGCompositeFrame(tf,w,h,kVerticalFrame);
   zListBox = new XSElementList(zFrame,XSEL_SORTBY_Z);
   zListBox->SelectZ(*selectedZ);
   zFrame->AddFrame(zListBox, lHints);
   tf->AddFrame(zFrame, lHints);

   // ------------- Finalise settings -----------------
   lHints2 = new TGLayoutHints(
                   kLHintsBottom | kLHintsExpandX | kLHintsExpandY,
         2, 2, 5, 1);
   AddFrame(tabMenu, lHints2);

   tabMenu->SetTab(activeTab);

   /* --------- Set Windows Position --------- */
   int     ax, ay;
   if (LastWinX == -1) {   // Go to the middle of parent window
      Window_t wdum;
      gVirtualX->TranslateCoordinates(main->GetId(), GetParent()->GetId(),
         (((TGFrame *) main)->GetWidth() - fWidth) >> 1,
         (((TGFrame *) main)->GetHeight() - fHeight) >> 1,
         ax, ay, wdum);
   } else {
      ax = LastWinX;
      ay = LastWinY;
   }
   Move(ax,ay);
   SetWMPosition(ax,ay);

   MapSubwindows();
   Resize(GetDefaultSize());
   SetWindowName("Select Element");
   MapWindow();

   fClient->WaitFor(this);
} // XSElementDlg

/* ----- ~XSElementDlg ----- */
XSElementDlg::~XSElementDlg()
{
   activeTab = tabMenu->GetCurrent();

   delete   closeButton;
   delete   buttonFrame;
   delete   buttonLayoutHints;
   delete   frameLayoutHints;
   delete   lHints;
   delete   tabMenu;
   delete   pTable;
   delete   nameListBox;
   delete   nameFrame;
   delete   mnemonicListBox;
   delete   mnemonicFrame;
   delete   zListBox;
   delete   zFrame;
   delete   lHints2;
} // ~XSElementDlg

/* ----- CloseWindow ----- */
void
XSElementDlg::CloseWindow()
{
   // --- Remember old position ---
   Window_t wdum;
   gVirtualX->TranslateCoordinates(GetId(), GetParent()->GetId(),
      0, 0, LastWinX, LastWinY, wdum);

   delete this;
} // CloseWindow

/* ----- ProcessButton ----- */
Bool_t
XSElementDlg::ProcessButton(Long_t param)
{
   Int_t   activeTab;

   switch (param) {
      case ELEM_OK:
         // Find active Tab
         activeTab = tabMenu->GetCurrent();
         if (activeTab==0)
            goto LERROR;
         else {
            switch (activeTab) {
               case 1:
                  *selectedZ =
                     nameListBox->CurrentZ();
                  break;
               case 2:
                  *selectedZ =
                     mnemonicListBox->CurrentZ();
                  break;
               case 3:
                  *selectedZ =
                     zListBox->CurrentZ();
                  break;
            }
            if (*selectedZ==0)
               goto LERROR;
         }
         break;

      case ELEM_CANCEL:
         *selectedZ = 0;
         break;
      default:
         *selectedZ = param;
         break;
   }
   CloseWindow();
   return kTRUE;

LERROR:
   // Issue a message to select an Element
   int     retval;
   new TGMsgBox(fClient->GetRoot(), this, "Info",
      "Please Select one Element!",
      kMBIconAsterisk, kMBOk, &retval );
   return kFALSE;
} // ProcessButton

/* ----- ProcessMessage ----- */
Bool_t
XSElementDlg::ProcessMessage(Long_t msg, Long_t param1, Long_t /*param2*/)
{
   switch (GET_MSG(msg)) {
      case kC_COMMAND:
         switch (GET_SUBMSG(msg)) {
            case kCM_BUTTON:
               ProcessButton(param1);
               break;

            default:
               break;
         }
      default:
         break;
   }
   return kTRUE;
} // ProcessMessage
