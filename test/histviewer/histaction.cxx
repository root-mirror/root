//---------------------------------------------------------------
//
//      Control Panel for histogram operations
//      Demo version
//
//      Author: Dmitry Vasiliev (LNS, Catania)
//
//---------------------------------------------------------------

#include "histaction.h"
#include <TCollection.h>
#include <TFile.h>
#include <TKey.h>
#include <TROOT.h>
#include <TGMsgBox.h>
#include <TSystem.h>
#include <TH1.h>
#include "canvsave.h"
#include <TCanvas.h>
#include <snprintf.h>

#if !defined S_ISDIR
#define S_ISDIR(m) (((m)&(0170000)) == (0040000))
#endif

ClassImp(HistAction);


Bool_t kSTATUS; //to suppress call of drawHist() when no histo is highlighted
                //in ListBoxA
TCanvas *cA;    //canvas in fCanvasA
TCanvas *cB;    //canvas in fCanvasB
TH1F *hist;     //clone histo in Scan panel
TH1F *temph;    //bin highlight in canvas cB
Long_t cursA, cursB; //current position in ListBoxA and ListBoxB
Int_t totalA, totalB;  //total number of entries in ListBoxA and ListBoxB
Int_t *array;     //array of non empty bins of a histo under scan
Int_t count;      //number of user highlighted buttons in Display Layout
Int_t indi[2];    //survice arrays for algorithm which calculates
Int_t indj[2];    //display layout


HistAction::HistAction(const TGWindow *p, UInt_t w, UInt_t h) :
      TGMainFrame(p, w, h)
{
   //--------------------------------------------------------------------
   // Constructor for the class HistAction. Draws a control panel, divides it
   // into sub areas and maps all subwindows
   //--------------------------------------------------------------------


   Int_t i;

   for (i = 0; i < 16; i++) {
      pads[i] = NULL;
      histInd[i] = -1;
   }
   fHisto = new TObjArray(kMaxHist);
   position = 0;
   totalA = 0;
   resetIter();
   resetFlags();


   /*

                                --------------------------
                                |                        |
                                |                        |
                                -------------fF0----------
                                |                        |
                                --------------------------
                               /                         \
                       ---------------               -------------------
                       |      |      |               |Import \ Scan    |
                       |     fFA     |               |        ---------|
                       |      |      |               |                 |
                       |      |      |               |      fTab       |
                       ---------------               -------------------

   */


   fF0 = new TGCompositeFrame(this, 200, 20, kVerticalFrame);
   AddFrame(fF0, new TGLayoutHints(kLHintsLeft | kLHintsTop));


   fFA = new TGCompositeFrame(fF0, 200, 20, kHorizontalFrame);
   fF0->AddFrame(fFA, new TGLayoutHints(kLHintsTop | kLHintsLeft));


   fTab = new TGTab(fF0, 200, 20);
   fF0->AddFrame(fTab, new TGLayoutHints(kLHintsLeft | kLHintsTop));

   /*
   //-------------------------------------------------------------------
   //
   //  Panel A (main presentation canvas, list box with the names of
   //           histogrammes, control buttons, display layout, close button)
   //
   //-------------------------------------------------------------------
   //
   //
   //               --------------------
   //               |      |           |
   //               |     fFA          |
   //               |      |           |
   //               --------------------
   //              /                    \
   //   ---------------            ----------------
   //   |             |            |              |
   //   |             |            |----fA1-------|
   //   |  fCanvasA   |            |              |
   //   |             |            |              |
   //   ---------------            ----------------
   //                              /               \
   //                             /                 \
   //                      ---------------       ----------------
   //                      |             |       |      |       |
   //              Control |     fA2     |       |      |       |
   //              Buttons |             |       |     fA3      |
   //                      |Matrix Layout|       |      |       |
   //                      |             |       |      |       |
   //                      ---------------       ----------------
   //                                           /               \
   //                                          /                 \  fA4
   //                                 --------------         --------------
   //                                 |            |         |fMultiButton|
   //                                 |            |         |------------|
   //                                 | fListBoxA  |         |fPrevButtonA|
   //                                 |            |         |------------|
   //                                 |            |         |fNextButtonA|
   //                                 --------------         |------------|
   //                                                        |fA5(Display |
   //                                                        |    Layout) |
   //                                                        |------------|
   //                                                        |fCloseButton|
   //                                                        --------------
   //
   /////////////////////////////////////////////////////////////////////////////

   */

   fCanvasA = new TRootEmbeddedCanvas("canvasA", fFA, 400, 400);
   fFA->AddFrame(fCanvasA, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5));
   cA = fCanvasA->GetCanvas();
   //cA->SetCanvasSize(396,396);
   cA->SetFillColor(10);


   fA1 = new TGCompositeFrame(fFA, 100, 20, kVerticalFrame);
   fFA->AddFrame(fA1, new TGLayoutHints(kLHintsTop | kLHintsLeft));


   fA2 = new TGGroupFrame(fA1, "Control Buttons", kHorizontalFrame);
   fA2->SetLayoutManager(new TGMatrixLayout(fA2, 0, 3, 8));
   fA1->AddFrame(fA2, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5));


   fA3 = new TGCompositeFrame(fA1, 200, 50, kHorizontalFrame);
   fA1->AddFrame(fA3, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5));


   fListBoxA = new TGListBox(fA3, M_LIST_A);
   fListBoxA->Associate(this);
   fA3->AddFrame(fListBoxA, new TGLayoutHints(kLHintsTop | kLHintsLeft));
   fListBoxA->Resize(140, 305);
   fListBoxA->Connect("DoubleClicked(char*)", "HistAction", this, "doubleclickedBoxA(char*)");


   fA4 = new TGCompositeFrame(fA3, 100, 20, kVerticalFrame);
   fA3->AddFrame(fA4, new TGLayoutHints(kLHintsTop | kLHintsLeft));


   fMultiButton = new TGCheckButton(fA4, "Multiple selection", M_MULTI);
   fMultiButton->Associate(this);
   fMultiButton->SetState(kButtonUp);
   fA4->AddFrame(fMultiButton, new TGLayoutHints(kLHintsTop | kLHintsLeft,
                 5, 5, 5, 5));


   fPrevButtonA = new TGTextButton(fA4, "Previous", M_PREV_A);
   fPrevButtonA->Associate(this);
   fA4->AddFrame(fPrevButtonA, new TGLayoutHints(kLHintsLeft | kLHintsTop,
                 10, 5, 5, 5));


   fNextButtonA = new TGTextButton(fA4, "   Next   ", M_NEXT_A);
   fNextButtonA->Associate(this);
   fA4->AddFrame(fNextButtonA, new TGLayoutHints(kLHintsLeft | kLHintsTop,
                 10, 5, 5, 5));


   fA5 = new TGGroupFrame(fA4, "Display Layout", kVerticalFrame);
   fA4->AddFrame(fA5, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5));
   fA5->SetLayoutManager(new TGMatrixLayout(fA5, 0, 4, kVerticalFrame));


   for (i = 0;i < 16;i++) {
      fLayoutButton[i] = new TGTextButton(fA5, "     ", 500 + i);
      fLayoutButton[i]->Associate(this);
      fA5->AddFrame(fLayoutButton[i], new TGLayoutHints(kLHintsTop | kLHintsLeft,
                    5, 5, 5, 5));
   }


   fDrawButton = new TGTextButton(fA2, "Draw", M_DRAW);
   fDrawButton->Associate(this);
   fDrawButton->Resize(65, fDrawButton->GetDefaultHeight());
   fA2->AddFrame(fDrawButton);


   fSelectButton = new TGTextButton(fA2, "Select All", M_SELECT);
   fSelectButton->Associate(this);
   fSelectButton->Resize(65, fSelectButton->GetDefaultHeight());
   fA2->AddFrame(fSelectButton);


   fClearButtonA = new TGTextButton(fA2, "Clear", M_CLEAR_A);
   fClearButtonA->Associate(this);
   fClearButtonA->Resize(65, fClearButtonA->GetDefaultHeight());
   fA2->AddFrame(fClearButtonA);


   fSaveButton = new TGTextButton(fA2, "Save Pic", M_SAVE);
   fSaveButton->Associate(this);
   fSaveButton->Resize(65, fSaveButton->GetDefaultHeight());
   fA2->AddFrame(fSaveButton);


   fEditButton = new TGTextButton(fA2, "Edit Pic", M_EDIT);
   fEditButton->Associate(this);
   fEditButton->Resize(65, fEditButton->GetDefaultHeight());
   fA2->AddFrame(fEditButton);


   fCloseButton = new TGTextButton(fA4, "Close Window", M_CLOSE);
   fCloseButton->Associate(this);
   SetWindowAttributes_t wattr;
   wattr.fMask = kWABackPixel;
   gClient->GetColorByName("red", wattr.fBackgroundPixel);
   gVirtualX->ChangeWindowAttributes(fCloseButton->GetId(), &wattr);
   fA4->AddFrame(fCloseButton, new TGLayoutHints(kLHintsTop | kLHintsLeft |
                 kLHintsCenterX, 0, 0, 40, 5));


   /*
   //---------------------------------------------------------------------
   //
   // Panel C (File Manager)
   //
   //---------------------------------------------------------------------
   //
   //
   //                   ----------------------------------
   //                   |                                |
   //                   |  tf = fTab->AddTab("Import")   |
   //                   |                                |
   //                   ----------------------------------
   //                    /                             \
   //                   /                               \
   //          ----------------------              ----------------------
   //          |  |  |    |    |    |              |                    |
   //          |  |  | fC1|    |    |              |        fC2         |
   //          |  |  |    |    |    |              |                    |
   //          ----------------------              ----------------------
   //         /    |    |     |      \                         |
   //   ------  -----  -----  -----  ---------     ----------------------
   //   |    |  |   |  |   |  |   |  |       |     |                    |
   //   |fDir|  |   |  |   |  |   |  |fFilter|     |       fFileView    |
   //   |    |  |   |  |   |  |   |  |  Box  |     |                    |
   //   ------  -----  -----  -----  ---------     ----------------------
   //           fCdup  fList  fDetail
   //                   Mode   Mode
   //
   //
   /////////////////////////////////////////////////////////////////////////////

   */

   TGCompositeFrame *tf = fTab->AddTab("Import");
   fC1 = new TGCompositeFrame(tf, 500, 50, kHorizontalFrame);
   tf->AddFrame(fC1, new TGLayoutHints(kLHintsLeft | kLHintsTop, 5, 5, 5, 5));


   fC2 = new TGCompositeFrame(tf, 500, 300, kVerticalFrame);
   tf->AddFrame(fC2, new TGLayoutHints(kLHintsLeft | kLHintsTop, 5, 5, 5, 5));


   fDir = new TGListBox(fC1, -1);
   //fDir->Associate(this);
   fDir->Resize(400, 20);
   TString tmp;
   tmp.Form("%s", gSystem->WorkingDirectory());
   fDir->AddEntry(tmp, 1);
   fC1->AddFrame(fDir, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2, 2, 2, 2));


   fPcdup = fClient->GetPicture("tb_uplevel.xpm");
   fCdup = new TGPictureButton(fC1, fPcdup, M_CDUP);
   fCdup->SetToolTipText("Up One Level");
   fCdup->Associate(this);
   fC1->AddFrame(fCdup, new TGLayoutHints(kLHintsLeft | kLHintsTop, 5, 5, 5, 5));


   fPlist = fClient->GetPicture("tb_list.xpm");
   fListMode = new TGPictureButton(fC1, fPlist, M_LIST_MODE);
   fListMode->SetToolTipText("List Mode");
   fListMode->Associate(this);
   fListMode->SetState(kButtonUp);
   fListMode->AllowStayDown(kTRUE);
   fC1->AddFrame(fListMode, new TGLayoutHints(kLHintsLeft | kLHintsTop, 5, 5, 5, 5));


   fPdetail = fClient->GetPicture("tb_details.xpm");
   fDetailMode = new TGPictureButton(fC1, fPdetail, M_DETAIL_MODE);
   fDetailMode->SetToolTipText("Details Mode");
   fDetailMode->Associate(this);
   fDetailMode->SetState(kButtonEngaged);
   fDetailMode->AllowStayDown(kTRUE);
   fC1->AddFrame(fDetailMode, new TGLayoutHints(kLHintsLeft | kLHintsTop,
                 5, 5, 5, 5));


   fC1->AddFrame(new TGLabel(fC1, new TGHotString("Filter:")),
                 new TGLayoutHints(kLHintsLeft | kLHintsTop, 10, 5, 5, 5));


   fFilterBox = new TGComboBox(fC1, M_FILTER);
   fFilterBox->Associate(this);
   char filt[50];
   snprintf(filt,50, "%s", "All files  (*)");
   fFilterBox->AddEntry(filt, 1);
   snprintf(filt,50, "%s", "ROOT files  (*.root)");
   fFilterBox->AddEntry(filt, 2);
   fFilterBox->Resize(150, 20);
   fFilterBox->Select(2);
   fC1->AddFrame(fFilterBox, new TGLayoutHints(kLHintsLeft | kLHintsTop, 5, 5, 5, 5));


   fFileView = new TGListView(fC2, 540, 180);
   fFileCont = new TGFileContainer(fFileView->GetViewPort(), 100, 100,
                                   kVerticalFrame, GetWhitePixel());
   fFileCont->Associate(this);
   fFileView->GetViewPort()->SetBackgroundColor(GetWhitePixel());
   fFileView->SetContainer(fFileCont);
   fFileCont->SetFilter("*.root");
   fFileCont->ChangeDirectory(".");
   fFileView->SetViewMode(kLVDetails);
   fFileCont->Sort(kSortByName);
   fC2->AddFrame(fFileView, new TGLayoutHints(kLHintsLeft | kLHintsTop,
                 2, 2, 2, 2));


   /*
   //---------------------------------------------------------------------
   //
   //  Panel B (histo presentation canvas, list box with the indices of non
   //           empty bins, control buttons, fields for displaying digital
   //           information about a histogram)
   //
   //---------------------------------------------------------------------
   //
   //
   //                        ------------------------------
   //                        |                            |
   //                        |  tf = fTab->AddTab("Scan") |
   //                        |                            |
   //                        ------------------------------
   //                                      |
   //                        ------------------------------
   //                        |         |         |        |
   //                        |         |   fFB   |        |
   //                        |         |         |        |
   //                        ------------------------------
   //                       /              |               \
   //                 fB1  /               |                \
   //              ------------      -------------     -----------------
   //              | fName    |      |           |     |               |
   //              |----------|      |  fCanvasB |     |-----fB2-------|
   //              | fTitle   |      |           |     |               |
   //      Matrix  |----------|      -------------     -----------------
   //      Layout  | fChan    |                       /                \
   //              |----------|                      /                  \
   //              | fRange   |             ---------------         -------------
   //              |----------|             |     |       |         |           |
   //              | fEntries |             |    fB4      |         |----fB3----|
   //              ------------             |     |       |         |           |
   //                                       ---------------         -------------
   //                                      /           |           /          |
   //                                     /            |          /           |
   //                               -----------   ----------   --------   -------
   //                               |         |   |        |   |      |   |     |
   //                               | fList   |   |   fB5  |   | fBin |   |fBin |
   //                               |  BoxB   |   |        |   | Cont |   |Range|
   //                               |         |   | Matrix |   |      |   |     |
   //                               -----------   | Layout |   --------   -------
   //                                             |        |
   //                                             ----------
   //                                               Control
   //                                               Buttons
   //
   //
   /////////////////////////////////////////////////////////////////////////////

   */

   tf = fTab->AddTab("Scan");
   fFB = new TGGroupFrame(tf, "Scan Panel", kHorizontalFrame);
   tf->AddFrame(fFB, new TGLayoutHints(kLHintsTop | kLHintsLeft));


   fB1 = new TGGroupFrame(fFB, "General Histo Info", kVerticalFrame);
   fB1->SetLayoutManager(new TGMatrixLayout(fB1, 0, 2, 10));
   fFB->AddFrame(fB1, new TGLayoutHints(kLHintsTop | kLHintsLeft |
                                        kLHintsExpandY));


   fCanvasB = new TRootEmbeddedCanvas("canvasB", fFB, 200, 200);
   fFB->AddFrame(fCanvasB, new TGLayoutHints(kLHintsTop | kLHintsLeft |
                 kLHintsExpandX | kLHintsExpandY,
                 5, 5, 5, 5));
   cB = fCanvasB->GetCanvas();
   cB->SetFillColor(10);


   fB2 = new TGCompositeFrame(fFB, 100, 100, kVerticalFrame);
   fFB->AddFrame(fB2, new TGLayoutHints(kLHintsTop | kLHintsLeft |
                                        kLHintsExpandY));


   fB2->AddFrame(new TGLabel(fB2, new TGHotString("Non empty bins")),
                 new TGLayoutHints(kLHintsLeft | kLHintsTop));


   fB4 = new TGCompositeFrame(fB2, 220, 100, kHorizontalFrame);
   fB2->AddFrame(fB4, new TGLayoutHints(kLHintsTop | kLHintsLeft));


   fListBoxB = new TGListBox(fB4, M_LIST_B);
   fListBoxB->Associate(this);
   fB4->AddFrame(fListBoxB, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 0, 5, 0));
   fListBoxB->Resize(70, 100);


   fB5 = new TGGroupFrame(fB4, "Control Buttons", kVerticalFrame);
   fB5->SetLayoutManager(new TGMatrixLayout(fB5, 0, 2, 10));
   fB4->AddFrame(fB5, new TGLayoutHints(kLHintsTop | kLHintsLeft, 20, 0, 0, 0));


   fPrevButtonB = new TGTextButton(fB5, "Previous", M_PREV_B);
   fPrevButtonB->Associate(this);
   fPrevButtonB->Resize(55, fPrevButtonB->GetDefaultHeight());
   fB5->AddFrame(fPrevButtonB, new TGLayoutHints(kLHintsLeft | kLHintsTop,
                 10, 1, 5, 5));


   fGetButton = new TGTextButton(fB5, "Import", M_IMPORT);
   fGetButton->Associate(this);
   fGetButton->Resize(55, fGetButton->GetDefaultHeight());
   fB5->AddFrame(fGetButton, new TGLayoutHints(kLHintsTop | kLHintsLeft,
                 1, 2, 5, 5));


   fNextButtonB = new TGTextButton(fB5, "Next", M_NEXT_B);
   fNextButtonB->Associate(this);
   fNextButtonB->Resize(55, fNextButtonB->GetDefaultHeight());
   fB5->AddFrame(fNextButtonB, new TGLayoutHints(kLHintsLeft | kLHintsTop,
                 10, 1, 5, 5));


   fClearButtonB = new TGTextButton(fB5, "Clear", M_CLEAR_B);
   fClearButtonB->Associate(this);
   fClearButtonB->Resize(55, fClearButtonB->GetDefaultHeight());
   fB5->AddFrame(fClearButtonB, new TGLayoutHints(kLHintsTop | kLHintsLeft,
                 1, 2, 5, 5));


   fB3 = new TGCompositeFrame(fB2, 100, 20, kVerticalFrame);
   fB3->SetLayoutManager(new TGMatrixLayout(fB3, 0, 2, 10));
   fB2->AddFrame(fB3, new TGLayoutHints(kLHintsTop | kLHintsLeft));


   fB1->AddFrame(new TGLabel(fB1, new TGHotString("Name")),
                 new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX));


   fName = new TGTextEntry(fB1, fNameBuf = new TGTextBuffer(100), -1);
   fB1->AddFrame(fName, new TGLayoutHints(kLHintsTop | kLHintsLeft, 10, 2, 2, 2));
   fName->Resize(100, fName->GetDefaultHeight());


   fB1->AddFrame(new TGLabel(fB1, new TGHotString("Title")),
                 new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX));


   fTitle = new TGTextEntry(fB1, fTitleBuf = new TGTextBuffer(100), -1);
   fB1->AddFrame(fTitle, new TGLayoutHints(kLHintsTop | kLHintsLeft, 10, 2, 2, 2));
   fTitle->Resize(100, fTitle->GetDefaultHeight());


   fB1->AddFrame(new TGLabel(fB1, new TGHotString("Channels")),
                 new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX));


   fChan = new TGTextEntry(fB1, fChanBuf = new TGTextBuffer(100), -1);
   fB1->AddFrame(fChan, new TGLayoutHints(kLHintsTop | kLHintsLeft, 10, 2, 2, 2));
   fChan->Resize(100, fChan->GetDefaultHeight());


   fB1->AddFrame(new TGLabel(fB1, new TGHotString("Range")),
                 new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX));


   fRange = new TGTextEntry(fB1, fRangeBuf = new TGTextBuffer(100), -1);
   fB1->AddFrame(fRange, new TGLayoutHints(kLHintsTop | kLHintsLeft, 10, 2, 2, 2));
   fRange->Resize(100, fRange->GetDefaultHeight());


   fB1->AddFrame(new TGLabel(fB1, new TGHotString("Entries")),
                 new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX));


   fEntries = new TGTextEntry(fB1, fEntriesBuf = new TGTextBuffer(100), -1);
   fB1->AddFrame(fEntries, new TGLayoutHints(kLHintsTop | kLHintsLeft, 10, 2, 2, 2));
   fEntries->Resize(100, fEntries->GetDefaultHeight());


   fB3->AddFrame(new TGLabel(fB3, new TGHotString("Bin content")),
                 new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX));


   fBinCont = new TGTextEntry(fB3, fBinContBuf = new TGTextBuffer(100), -1);
   fB3->AddFrame(fBinCont, new TGLayoutHints(kLHintsTop | kLHintsLeft,
                 0, 2, 2, 2));
   fBinCont->Resize(110, fBinCont->GetDefaultHeight());


   fB3->AddFrame(new TGLabel(fB3, new TGHotString("Bin range")),
                 new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX));


   fBinRange = new TGTextEntry(fB3, fBinRangeBuf = new TGTextBuffer(100), -1);
   fB3->AddFrame(fBinRange, new TGLayoutHints(kLHintsTop | kLHintsLeft,
                 0, 2, 2, 2));
   fBinRange->Resize(110, fBinRange->GetDefaultHeight());


   //----------------------------------------------------------------

   MapSubwindows();

   SetWindowName("Histogram Viewer");

   Resize(GetDefaultSize());
   MapWindow();
}

HistAction::~HistAction()
{
   //----------------------------------------------------------
   //
   //    Destructor for the class HistAction
   //
   //----------------------------------------------------------

   delete fCloseButton;
   delete fDrawButton;
   delete fGetButton;
   delete fSelectButton;
   delete fSaveButton;
   delete fEditButton;
   delete fClearButtonB;
   delete fClearButtonA;
   delete fPrevButtonA;
   delete fNextButtonA;
   delete fPrevButtonB;
   delete fNextButtonB;
   delete fMultiButton;
   delete fCdup;
   delete fListMode;
   delete fDetailMode;
   delete fLayoutButton[0];
   delete fLayoutButton[1];
   delete fLayoutButton[2];
   delete fLayoutButton[3];
   delete fLayoutButton[4];
   delete fLayoutButton[5];
   delete fLayoutButton[6];
   delete fLayoutButton[7];
   delete fLayoutButton[8];
   delete fLayoutButton[9];
   delete fLayoutButton[10];
   delete fLayoutButton[11];
   delete fLayoutButton[12];
   delete fLayoutButton[13];
   delete fLayoutButton[14];
   delete fLayoutButton[15];
   delete fListBoxA;
   delete fListBoxB;
   delete fFilterBox;
   delete fDir;
   delete fName;
   delete fTitle;
   delete fChan;
   delete fRange;
   delete fEntries;
   delete fBinCont;
   delete fBinRange;
   delete fC1;
   delete fC2;
   delete fTab;
   delete fA5;
   delete fA4;
   delete fA3;
   delete fA2;
   delete fCanvasA;
   delete fA1;
   delete fB5;
   delete fB4;
   delete fB3;
   delete fB2;
   delete fCanvasB;
   delete fB1;
   delete fFA;
   delete fFB;
   delete fF0;
}

void HistAction::CloseWindow()
{
   //-----------------------------------------------------
   //
   //           CloseWindow()
   //
   //     Closes the panel "Histogram Viewer"
   //
   //-----------------------------------------------------

   delete this;

}

Bool_t HistAction::importHist(const char *name)
{
   //------------------------------------------------------------------
   //
   //               importHist(const char *name)
   //
   // Allows to import an existing histogram from the memory
   // (needs histo name as an input parameter).
   // This function is not used in GUI at the moment, one can call it from
   // the interpreter.
   //
   // Example:
   //           gui = new HistAction(gClient->GetRoot(),1,1);
   //           TH1F *hist = new TH1F("myhisto","bla-bla-bla",100,0,100);
   //           {
   //           ...
   //           hist->Fill(x);
   //           ...
   //           }
   //           gui->importHist(hist->GetName());
   //
   //-------------------------------------------------------------------

   TH1F *h;
   h = (TH1F*) gROOT->FindObject(name);
   if (!h) return kFALSE;
   if (position == kMaxHist) return kFALSE;
   fHisto->AddAt(h, position++);
   fListBoxA->AddEntry(h->GetName(), ++totalA);
   fListBoxA->MapSubwindows();
   fListBoxA->Layout();
   return kTRUE;
}

Bool_t HistAction::importFromFile(const char *filename)
{
   //-------------------------------------------------------------------
   //
   //                 importFromFile(const char *filename)
   //
   // Imports histogrammes from a file with the name "filename".
   // Opens the file, scans it, if finds an object of the class TH1F or TH2F,
   // imports it.
   // All the other objects are ignored.
   // In case of not a ROOT file returns an error message and takes no further
   // action.
   // This function is called when a user doubly clicks on the file icon
   // in the file list view.
   //
   //--------------------------------------------------------------------

   Int_t retval;
   Int_t l;
   TFile *f;
   TH1F *fH;
   f = new TFile(filename);
   if (!f) return kFALSE;
   if (f->IsZombie()) {
      new TGMsgBox(fClient->GetRoot(), this, "Error Message",
                   "You have chosen not a ROOT file. Please, be attentive.",
                   kMBIconExclamation, kMBOk, &retval);
      delete f;
      f = NULL;
      return kFALSE;
   }
   TKey *key;
   TIter it(f->GetListOfKeys());
   while ((key = (TKey*) it())) {
      if (!strcmp(key->GetClassName(), "TH1F") ||
            !strcmp(key->GetClassName(), "TH2F")) {
         fH = (TH1F*) key->ReadObj();
         if (fH && position < kMaxHist) fHisto->AddAt(fH, position++);
      }
   }
   for (l = totalA; l < position; l++) {
      fListBoxA->AddEntry(fHisto->At(l)->GetName(), l + 1);
   }
   fListBoxA->MapSubwindows();
   fListBoxA->Layout();
   totalA = position;
   return kTRUE;
}

int HistAction::getNextTrueIndex()
{
   //---------------------------------------------------------------
   //
   //                 getNextTrueIndex()
   //
   // Iterates over array "flags", returns the next "TRUE index".
   // In case of no "TRUE index" found returns -1.
   // "TRUE index" means the index of a histogram currently highlighted in
   // the large list box (ListBoxA).
   //
   //---------------------------------------------------------------

   while (cursorIter < kMaxHist) {
      cursorIter++;
      if (flags[cursorIter]) return cursorIter;
   }
   return -1;
}

void HistAction::setCanvasDivision(Int_t number)
{
   //---------------------------------------------------------------
   //
   //              setCanvasDivision(Int_t number)
   //
   // Calculates the xDiv and yDiv parameters which are used to divide
   // the main canvas (CanvasA) into subpads.
   // This function is used in case of the automatic display layout
   // (checkbutton "Multiple selection" is engaged).
   // The function takes as an input parameter the total number of histogrammes
   // to be displayed.
   // Called from drawHist()
   //
   //---------------------------------------------------------------

   Int_t i, j, k;
   for (i = 1; i < 50; i++) {
      k = 0;
      for (j = i - 1; j <= i; j++) {
         if (number <= i*j) {
            k = j;
            break;
         }
      }
      if (number <= i*k) break;
   }
   xDiv = i;
   yDiv = k;
}

void HistAction::drawHist()
{
   //-----------------------------------------------------------------
   //
   //                drawHist()
   //
   // Draws a set of histogrammes in the canvas cA in case of the automatic
   // display layout (checkbutton "Multiple selection" is engaged).
   // Called when the button "Draw" is clicked.
   //
   //----------------------------------------------------------------

   Int_t number;  //number of highlighted histos in ListBoxA
   Int_t i;
   number = -1;
   resetIter();
   while (getNextTrueIndex() != -1) number++;
   setCanvasDivision(number);
   cA->Clear();
   cA->Divide(xDiv, yDiv);
   resetIter();
   for (i = 0; i < number; i++) {
      cA->cd(i + 1);
      ((TH1F*) fHisto->At(getNextTrueIndex()))->Draw();
   }
   cA->cd();
   cA->Modified();
   cA->Update();
}

void HistAction::toScan()
{
   //---------------------------------------------------------------
   //
   //                     toScan()
   //
   // Makes a clone of a selected histo.
   // Puts the clone object to the scan panel (panel B).
   // Called when the button "Import" is clicked.
   //
   //---------------------------------------------------------------

   Int_t retval;
   Int_t chan;
   Float_t from, to;
   Int_t entries;
   TObject *p;
   TString tmp;
   Int_t i;
   Float_t bcont;

   clearScan();
   TIter it(gPad->GetListOfPrimitives());
   while ((p = it())) {
      if (!strcmp(p->ClassName(), "TH1F")) hist = (TH1F*) p->Clone();
   }
   if (!hist) {
      new TGMsgBox(fClient->GetRoot(), this, "Help message",
                   "Choose a histogram by clicking on its area the middle button of the mouse.",
                   kMBIconExclamation, kMBOk, &retval);
      return;
   }
   cB->cd();
   if (hist) {

      hist->Draw();
      chan = hist->GetXaxis()->GetNbins();
      from = hist->GetXaxis()->GetXmin();
      to = hist->GetXaxis()->GetXmax();
      entries = (Int_t) hist->GetEntries();

      tmp.Form("%s", hist->GetName());
      fNameBuf->Clear();
      fNameBuf->AddText(0, tmp);
      fClient->NeedRedraw(fName);

      tmp.Form("%s", hist->GetTitle());
      fTitleBuf->Clear();
      fTitleBuf->AddText(0, tmp);
      fClient->NeedRedraw(fTitle);

      tmp.Form("%i", chan);
      fChanBuf->Clear();
      fChanBuf->AddText(0, tmp);
      fClient->NeedRedraw(fChan);

      tmp.Form("%.1f%s%.1f", from, "...", to);
      fRangeBuf->Clear();
      fRangeBuf->AddText(0, tmp);
      fClient->NeedRedraw(fRange);

      tmp.Form("%i", entries);
      fEntriesBuf->Clear();
      fEntriesBuf->AddText(0, tmp);
      fClient->NeedRedraw(fEntries);

      array = new Int_t[chan];
      for (i = 1; i <= chan; i++) {
         bcont = (float) hist->GetBinContent(i);
         if (bcont > 0.0) {
            tmp.Form("%i", i);
            array[totalB] = i;
            fListBoxB->AddEntry(tmp, ++totalB);
         }
      }
      fListBoxB->MapSubwindows();
      fListBoxB->Layout();
   }
   cB->Update();
}

void HistAction::processBoxB(Int_t par)
{
   //-----------------------------------------------------------
   //
   //        processBoxB(Int_t par)
   //
   // Processes information from the listbox "Non empty bins" (ListBoxB).
   // Updates the fields "Bin content" and "Bin range".
   // Higlights the selected bin on the canvas cB.
   // Called when an item in the listbox "Non empty bins" is highlighted.
   // Input parameter - bin index.
   //
   //------------------------------------------------------------

   Int_t chan;
   Float_t from, to;
   Float_t bfrom, bto;
   TString tmp;
   Float_t bcont;

   if (temph) {
      delete temph;
      temph = NULL;
   }
   if (hist) {
      chan = hist->GetXaxis()->GetNbins();
      from = hist->GetXaxis()->GetXmin();
      to = hist->GetXaxis()->GetXmax();
      temph = new TH1F("h", "h", chan, from, to);
      temph->SetLineColor(kRed);
      temph->SetFillColor(kRed);
      temph->SetBinContent(array[par-1], hist->GetBinContent(array[par-1]));
      cB->cd();
      temph->Draw("SAME");

      bcont = (float) hist->GetBinContent(array[par-1]);
      tmp.Form("%.2f", bcont);
      fBinContBuf->Clear();
      fBinContBuf->AddText(0, tmp);
      fClient->NeedRedraw(fBinCont);

      bfrom = hist->GetBinLowEdge(array[par-1]);
      bto = bfrom + hist->GetBinWidth(array[par-1]);
      tmp.Form("%.2f%s%.2f", bfrom, "...", bto);
      fBinRangeBuf->Clear();
      fBinRangeBuf->AddText(0, tmp);
      fClient->NeedRedraw(fBinRange);
   }
   cB->Update();
}

void HistAction::clearScan()
{
   //-----------------------------------------------------------
   //
   //              clearScan()
   //
   // Clears all the fields in the scan panel (panel B).
   //
   //-----------------------------------------------------------

   TVirtualPad *padsav = gPad;

   fNameBuf->Clear();
   fClient->NeedRedraw(fName);

   fTitleBuf->Clear();
   fClient->NeedRedraw(fTitle);

   fChanBuf->Clear();
   fClient->NeedRedraw(fChan);

   fRangeBuf->Clear();
   fClient->NeedRedraw(fRange);

   fEntriesBuf->Clear();
   fClient->NeedRedraw(fEntries);

   fBinContBuf->Clear();
   fClient->NeedRedraw(fBinCont);

   fBinRangeBuf->Clear();
   fClient->NeedRedraw(fBinRange);

   if (hist) {
      delete hist;
      hist = NULL;
   }

   if (temph) {
      delete temph;
      temph = NULL;
   }

   if (totalB) {
      fListBoxB->RemoveEntries(1, totalB);
      totalB = 0;
      fListBoxB->MapSubwindows();
      fListBoxB->Layout();
   }

   cursB = 0;

   if (array) {
      delete [] array;
      array = NULL;
   }

   cB->Clear();
   cB->Update();

   gPad = padsav;
}

Bool_t HistAction::toDefault(Window_t id)
{
   //------------------------------------------------------------------
   //
   //      toDefault(Window_t id)
   //
   // Changes the color attributes of a window to default values (gray).
   // Used to change the colors of the buttons in the panel "Display Layout".
   // Does not redraw the button.
   // So to visually change the color needs to be followed by
   // the function NeedRedraw(TGWindow *).
   // Input parameter - ID of the button to which the action must be applied.
   //
   //------------------------------------------------------------------

   SetWindowAttributes_t wattr;
   wattr.fMask = kWABackPixel;
   wattr.fBackgroundPixel = GetDefaultFrameBackground();
   gVirtualX->ChangeWindowAttributes(id, &wattr);
   return kTRUE;

}

Bool_t HistAction::toGreen(Window_t id)
{
   //-------------------------------------------------------------
   //
   //       toGreen(Window_t id)
   //
   // The same as above except changing the color to green.
   //
   //-------------------------------------------------------------

   SetWindowAttributes_t wattr;
   wattr.fMask = kWABackPixel;
   gClient->GetColorByName("green", wattr.fBackgroundPixel);
   gVirtualX->ChangeWindowAttributes(id, &wattr);
   return kTRUE;
}

Bool_t HistAction::isOverlap()
{
   //-------------------------------------------------------------
   //
   //         isOverlap()
   //
   // Checks if a selected display layout overlaps with already existing
   // pads in the canvas cA.
   //
   //-------------------------------------------------------------

   Int_t i, j;
   Int_t tmpIndex;

   for (i = 0; i < 4; i++) {
      for (j = 0; j < 4; j++) {
         if (verLay[i] && horLay[j]) {
            tmpIndex = 4 * i + j;
            if (histInd[tmpIndex] != -1) return kTRUE;
         }
      }
   }
   return kFALSE;

}

Bool_t HistAction::isLayout()
{
   //-------------------------------------------------------
   //
   //        isLayout()
   //
   // Checks if display layout is set.
   //
   //--------------------------------------------------------

   Int_t i;
   for (i = 0;i < 4;i++) {
      if (horLay[i] != 0) return kTRUE;
      if (verLay[i] != 0) return kTRUE;
   }
   return kFALSE;

}

void HistAction::paintHist()
{
   //--------------------------------------------------------------------
   //
   //              paintHist()
   //
   // Draws a histo in the canvas cA in case of the user defined display layout.
   // The latest display layout has the highest priority. If an overlap
   // with existing pads is detected, they are deleted from cA.
   // Algorithm virtually divides cA into subpads with the matrix layout (4x4).
   // A real pad in which histo will be drawn is constructed from virtual subpads.
   // The number of virtual subpads for the real pad can change in the range 1-16.
   // Arrays histInd[16] and pads[16] keep the "id" of the histo and the
   // address of the real pad
   //
   //            -----------------
   //            |   |   |   |   |
   //            | 1 | 2 | 3 | 4 |
   //            |---|---|---|---|
   //            |   |   |   |   |
   //            | 5 | 6 | 7 | 8 |
   //            |---|---|---|---|
   //            |   |   |   |   |
   //            | 9 | 10| 11| 12|
   //            |---|---|---|---|
   //            |   |   |   |   |
   //            | 13| 14| 15| 16|
   //            -----------------
   //
   //
   // If a histo with id=20 must be drawn in a pad which embraces virtual subpads
   // 1,2,5,6 then
   //              histInd[0] = 20        pads[0] = address of the real pad
   //              histInd[1] = 20        pads[1] = NULL
   //              histInd[4] = 20        pads[4] = NULL
   //              histInd[5] = 20        pads[5] = NULL
   //
   // To search for the pads to be deleted the algorithm uses only array
   // histInd[].
   // Only one of the virtual subpads of the real pad keeps the address
   // to avoid double deleting of the same object.
   // If there is an overlap between the pads which contain the histo with
   // the same "id", then only the latest version is drawn.
   // All the other pads with this histo (even non overlapping with the current
   // one) will be deleted from the canvas.
   // To have several versions of the same histo drawn in the canvas one has
   // to avoid pads overlapping when setting display layout.
   //--------------------------------------------------------------------

   Int_t retval;
   Float_t xmin = 0.0F;
   Float_t xmax = 0.0F;
   Float_t ymin = 0.0F;
   Float_t ymax = 0.0F;
   Int_t i, j, countLocal;
   Int_t ind;
   Int_t tempind;
   TPad *pad;
   const Float_t ratio = 0.25;

   if (!isLayout()) {
      new TGMsgBox(fClient->GetRoot(), this, "Message",
                   "Set Display Layout.",
                   kMBIconExclamation, kMBOk, &retval);
      return;
   }
   resetIter();
   ind = getNextTrueIndex();
   for (i = 0; i < 4; i++) {
      if (horLay[i] && (xmin == 0.0)) xmin = i * ratio + 0.01;
      if (horLay[i] && (xmin != 0.0)) xmax = (i + 1) * ratio - 0.01;
   }
   for (i = 3; i > -1; i--) {
      if (verLay[i] && (ymin == 0.0)) ymin = (3 - i) * ratio + 0.01;
      if (verLay[i] && (ymin != 0.0)) ymax = (4 - i) * ratio - 0.01;
   }
   if (isOverlap()) {
      for (i = 0; i < 16; i++) {
         if (verLay[i/4] && horLay[i%4]) {
            tempind = histInd[i];
            for (j = 0; j < 16; j++) {
               if (histInd[j] == tempind) {
                  histInd[j] = -1;
                  if (pads[j]) {
                     delete pads[j];
                     pads[j] = NULL;
                  }
               }
            }
         }
      }
   }
   pad = new TPad("pad", "pad", xmin, ymin, xmax, ymax);
   pad->SetFillColor(10);
   cA->cd();
   pad->Draw();
   pad->cd();
   if (fHisto->At(ind))((TH1F*) fHisto->At(ind))->Draw();
   cA->cd();
   cA->Modified();
   cA->Update();

   countLocal = 0;
   for (i = 0; i < 4; i++) {
      for (j = 0; j < 4; j++) {
         if (verLay[i] && horLay[j]) {
            countLocal++;
            histInd[4*i+j] = ind;
            if (countLocal == 1) pads[4*i+j] = pad;
            else pads[4*i+j] = NULL;
         }
      }
   }
   return;
}

Bool_t HistAction::ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2)
{
   //------------------------------------------------------------------
   //
   //    ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2)
   //
   // Processes information from all GUI items.
   // Selecting an item usually generates an event with 4 parameters.
   // The first two are packed into msg (first and second bytes).
   // The other two are parm1 and parm2.
   //
   //------------------------------------------------------------------

   Int_t retval; //for class TGMsgBox
   Int_t buttons;//used to construct message panel when Close button is clicked
   Int_t numb;//to update layout of list boxes
   Int_t i, j;
   Int_t imin, imax;//to calculate display layout
   Int_t jmin, jmax;//to calculate display layout
   TString tmp, filename, ext;//to process information in file container
   TGFileItem *f;//the same as above
   void *p = 0;//the same as above
   TString command;//the same as above

   switch (GET_MSG(msg)) {
      case kC_COMMAND:

         switch (GET_SUBMSG(msg)) {

            case kCM_BUTTON:

               switch (parm1) {
                  case M_DRAW:
                     if (!totalA) {
                        new TGMsgBox(fClient->GetRoot(), this, "Error message",
                                     "Histo container is empty.",
                                     kMBIconExclamation, kMBOk, &retval);
                        break;
                     }

                     if (!kSTATUS) {
                        new TGMsgBox(fClient->GetRoot(), this, "Help message",
                                     "Highlight the name of the histogram to be displayed.",
                                     kMBIconExclamation, kMBOk, &retval);
                        break;
                     }

                     if (fMultiButton->GetState()) drawHist();//Automatic display layout
                     else {
                        paintHist(); //User defined display layout

                        // Total number of buttons which can be set in "Display Layout" panel by a user
                        count = 2;
                     }
                     break;

                  case M_CLEAR_A:
                     for (int k = 0; k < 16; k++) {
                        histInd[k] = -1;
                        if (pads[k]) delete pads[k];
                        pads[k] = NULL;
                     }
                     cA->cd();
                     cA->Clear();
                     cA->Update();
                     break;

                  case M_SAVE:

                     // Ask confirmation to close the window

                     new CanvSave(fClient->GetRoot(), this, 400, 200);
                     fFileCont->DisplayDirectory();
                     break;

                  case M_EDIT:
                     cA->cd();
                     cA->EditorBar();
                     break;

                  case M_IMPORT:

                     // Put a histo to the Scan panel

                     toScan();
                     break;

                  case M_CLEAR_B:
                     clearScan();
                     break;

                  case M_CLOSE:
                     retval = 0;
                     buttons = 0;
                     buttons |= kMBYes;
                     buttons |= kMBNo;
                     new TGMsgBox(fClient->GetRoot(), this, "Confirm action",
                                  "Close Panel 'Histogram Viewer' ?",
                                  kMBIconQuestion, buttons, &retval);
                     if (retval == 1) CloseWindow();
                     break;

                  case M_PREV_B:

                     // One bin up in ListBoxB (non empty bins)

                     if (totalB > 0) {
                        if ((cursB > 0) && (cursB <= totalB)) cursB--;
                        if (cursB < 1) cursB = 1;
                        if (cursB > totalB) cursB = totalB;
                        fListBoxB->Select(cursB);
                        numb = cursB;
                        while (((--numb) % 5) != 0) { }
                        fListBoxB->SetTopEntry(++numb);
                        SendMessage(this, MK_MSG(kC_COMMAND, kCM_LISTBOX), M_LIST_B, cursB);
                     }
                     break;

                  case M_NEXT_B:

                     // One bin down in ListBoxB (non empty bins)

                     if (totalB > 0) {
                        if ((cursB > 0) && (cursB <= totalB)) cursB++;
                        if (cursB < 1) cursB = 1;
                        if (cursB > totalB) cursB = totalB;
                        fListBoxB->Select(cursB);
                        numb = cursB;
                        while (((--numb) % 5) != 0) { }
                        fListBoxB->SetTopEntry(++numb);
                        SendMessage(this, MK_MSG(kC_COMMAND, kCM_LISTBOX), M_LIST_B, cursB);
                     }
                     break;

                  case M_PREV_A:

                     // One histo up in ListBoxA

                     if (!fMultiButton->GetState() && (totalA > 0)) {
                        if ((cursA > 0) && (cursA <= totalA)) cursA--;
                        if (cursA < 1) cursA = 1;
                        if (cursA > totalA) cursA = totalA;
                        fListBoxA->Select(cursA);
                        numb = cursA;
                        while (((--numb) % 14) != 0) { }
                        fListBoxA->SetTopEntry(++numb);
                        SendMessage(this, MK_MSG(kC_COMMAND, kCM_LISTBOX), M_LIST_A, cursA);
                        SendMessage(this, MK_MSG(kC_COMMAND, kCM_BUTTON), M_DRAW, 0);
                     }
                     break;

                  case M_NEXT_A:

                     // One histo down in ListBoxA

                     if (!fMultiButton->GetState() && (totalA > 0)) {
                        if ((cursA > 0) && (cursA <= totalA)) cursA++;
                        if (cursA < 1) cursA = 1;
                        if (cursA > totalA) cursA = totalA;
                        fListBoxA->Select(cursA);
                        numb = cursA;
                        while (((--numb) % 14) != 0) { }
                        fListBoxA->SetTopEntry(++numb);
                        SendMessage(this, MK_MSG(kC_COMMAND, kCM_LISTBOX), M_LIST_A, cursA);
                        SendMessage(this, MK_MSG(kC_COMMAND, kCM_BUTTON), M_DRAW, 0);
                     }
                     break;

                  case M_SELECT:

                     // "Select All" button is clicked

                     fMultiButton->SetState(kButtonDown);
                     SendMessage(this, MK_MSG(kC_COMMAND, kCM_CHECKBUTTON), M_MULTI, 0);

                     // Call twice SetMultipleSelections(kTRUE), otherwise items in the list box
                     // are not highlighted (though proper functionality remains)

                     fListBoxA->SetMultipleSelections(kTRUE);
                     for (i = 1; i <= totalA; i++) {
                        fListBoxA->Select(i);
                        SendMessage(this, MK_MSG(kC_COMMAND, kCM_LISTBOX), M_LIST_A, i);
                     }
                     break;

                  case M_CDUP:

                     // To the parent directory

                     gSystem->ChangeDirectory("..");
                     tmp.Form("%s", gSystem->WorkingDirectory());
                     fDir->RemoveEntry(1);
                     fDir->AddEntry(tmp, 1);
                     fDir->MapSubwindows();
                     fDir->Layout();
                     fFileCont->DisplayDirectory();//otherwise the current directory is
                     //not always updated
                     break;

                  case M_LIST_MODE:

                     // "List" mode of file view

                     fFileView->SetViewMode(kLVList);
                     fDetailMode->SetState(kButtonUp);
                     break;

                  case M_DETAIL_MODE:

                     // "Details" mode of file view

                     fFileView->SetViewMode(kLVDetails);
                     fListMode->SetState(kButtonUp);
                     break;

                  default:

                     // Process the panel "Display Layout"

                     if (parm1 >= 500 && parm1 <= 515 && !fMultiButton->GetState()) {
                        if (count == 2) count = 0;
                        if (count < 2) {
                           toGreen(fLayoutButton[parm1-500]->GetId());
                           verLay[(parm1-500)/4] = 1;
                           horLay[(parm1-500)%4] = 1;
                           fClient->NeedRedraw(fLayoutButton[parm1-500]);
                           indi[count] = (parm1 - 500) / 4;
                           indj[count] = (parm1 - 500) % 4;
                           count++;
                           if (count == 2) {
                              imin = (indi[0] < indi[1]) ? indi[0] : indi[1];
                              imax = (indi[0] > indi[1]) ? indi[0] : indi[1];
                              jmin = (indj[0] < indj[1]) ? indj[0] : indj[1];
                              jmax = (indj[0] > indj[1]) ? indj[0] : indj[1];
                              for (i = 0;i < 4;i++) {
                                 for (j = 0;j < 4;j++) {
                                    if (i >= imin && i <= imax && j >= jmin && j <= jmax) {
                                       toGreen(fLayoutButton[4*i+j]->GetId());
                                       verLay[i] = 1;
                                       horLay[j] = 1;
                                    } else {
                                       toDefault(fLayoutButton[4*i+j]->GetId());
                                       if (i < imin || i > imax) verLay[i] = 0;
                                       if (j < jmin || j > jmax) horLay[j] = 0;
                                    }
                                    fClient->NeedRedraw(fLayoutButton[4*i+j]);
                                 }
                              }
                           }
                           if (count == 1) {
                              for (i = 0;i < 16;i++) {
                                 if (i != (parm1 - 500)) {
                                    toDefault(fLayoutButton[i]->GetId());
                                    if (i / 4 != (parm1 - 500) / 4) verLay[i/4] = 0;
                                    if (i % 4 != (parm1 - 500) % 4) horLay[i%4] = 0;
                                    fClient->NeedRedraw(fLayoutButton[i]);
                                 }
                              }
                           }
                        }
                     }
                     break;
               }

            case kCM_CHECKBUTTON:

               // Multiple selection

               switch (parm1) {
                  case M_MULTI:
                     if (!fListBoxA->GetMultipleSelections()) {
                        if (fListBoxA->GetSelectedEntry())
                           fListBoxA->GetSelectedEntry()->Activate(kFALSE);
                     }
                     fListBoxA->SetMultipleSelections(fMultiButton->GetState());

                     cursA = 0;
                     cA->Clear();
                     cA->Update();
                     for (i = 0; i < 16; i++) {
                        toDefault(fLayoutButton[i]->GetId());
                        fClient->NeedRedraw(fLayoutButton[i]);
                        verLay[i/4] = 0;
                        horLay[i%4] = 0;
                     }
                     count = 0;
                     for (j = 0; j < 16; j++) {
                        pads[j] = NULL;
                        histInd[j] = -1;
                     }
                     resetFlags();
                     kSTATUS = kFALSE;
                     break;

                  default:
                     break;
               }

            case kCM_LISTBOX:

               switch (parm1) {
                  case M_LIST_A:

                     // ListBoxA

                     cursA = parm2; //necessary for "Previous", "Next" buttons in case of
                     //random jumps in list box window

                     if (!fListBoxA->GetMultipleSelections()) {
                        resetFlags();
                        flags[parm2-1] = kTRUE;
                     } else
                        flags[parm2-1] = !flags[parm2-1];
                     kSTATUS = kTRUE;
                     break;

                  case M_LIST_B:

                     // ListBoxB (non empty bins)

                     cursB = parm2; //for "Previous", "Next" buttons
                     processBoxB(parm2);
                     break;

                  default:
                     break;
               }

            case kCM_COMBOBOX:
               switch (parm1) {

                  case M_FILTER:

                     // Set filter on the file type

                     if (parm2 == 1) fFileCont->SetFilter("*");
                     if (parm2 == 2) fFileCont->SetFilter("*.root");
                     fFileCont->DisplayDirectory();
                     break;

                  default:
                     break;
               }

            default:
               break;
         }

      case kC_CONTAINER:
         switch (GET_SUBMSG(msg)) {

            case kCT_ITEMDBLCLICK:

               // Process mouse double clicking in file view container

               if (parm1 == kButton1) {
                  if (fFileCont->NumSelected() == 1) {
                     f = (TGFileItem *) fFileCont->GetNextSelected(&p);
                     if (S_ISDIR(f->GetType())) {
                        fFileCont->ChangeDirectory(f->GetItemName()->GetString());
                        tmp.Form("%s", gSystem->WorkingDirectory());
                        fDir->RemoveEntry(1);
                        fDir->AddEntry(tmp, 1);
                        fDir->MapSubwindows();
                        fDir->Layout();
                     } else {
                        filename.Form("%s", f->GetItemName()->GetString());
                        ext = filename(filename.Last('.')+1, filename.Length());

                        // Call gv for postscript and pdf files

                        if (ext == "ps" || ext == "PDF" || ext == "pdf") {
                           command.Form("%s%s%s%s%s", "gv ",
                                        gSystem->WorkingDirectory(),
                                        "/",
                                        filename.Data(),
                                        "&");
                           gSystem->Exec(command);
                           break;
                        }

                        // Call xv for jpg, gif and bmp files

                        if (ext == "gif" || ext == "jpg" || ext == "bmp") {
                           command.Form("%s%s%s%s%s", "xv ",
                                        gSystem->WorkingDirectory(),
                                        "/",
                                        filename.Data(),
                                        "&");
                           gSystem->Exec(command);
                           break;
                        }

                        // Import root file

                        filename.Form("%s%s%s",
                                gSystem->WorkingDirectory(),
                                "/",
                                f->GetItemName()->GetString());
                        importFromFile(filename);
                     }
                  }
               }
               break;

            default:
               break;
         }

      default:
         break;
   }

   return kTRUE;

}

void HistAction::doubleclickedBoxA(const char * /*text*/)
{
   //------------------------------------------------------------------
   //
   //    doubleclickBoxA(const char *)
   //
   // Handle double click events in fListBoxA.
   // Double clicking in the list of histograms will draw the selected
   // histogram.
   //
   //------------------------------------------------------------------

   paintHist(); //User defined display layout
   // Total number of buttons which can be set in "Display Layout" panel by a user
   count = 2;
}
