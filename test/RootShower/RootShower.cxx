// Author: Bertrand Bellenot   22/08/02

/*************************************************************************
 * Copyright (C) 1995-2002, Bertrand Bellenot.                           *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see the LICENSE file.                         *
 *************************************************************************/

#include <ctime>
#include <iostream>
#include <string>

#include <TROOT.h>
#include <TStyle.h>
#include <TRint.h>
#include <TVirtualX.h>
#include <TEnv.h>
#include <KeySymbols.h>

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TFrame.h>
#include <TH1.h>
#include <TF1.h>

#include <TGMenu.h>
#include <TGFileDialog.h>
#include <TGTextEdit.h>
#include <TGToolTip.h>
#include <TG3DLine.h>
#include <TRootEmbeddedCanvas.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <TSystem.h>
#include <TRootHelpDialog.h>
#include <TGStatusBar.h>
#include <TBrowser.h>
#include <TParticle.h>
#include <TContextMenu.h>
#include "RootShower.h"
#include "MyParticle.h"
#include "GTitleFrame.h"
#include "GButtonFrame.h"
#include "RSMsgBox.h"
#include "RSAbout.h"
#include "SettingsDlg.h"
#include "RSHelpText.h"
#include "MyEvent.h"

#include <TGeoManager.h>
#include <TView.h>
#include <TGToolBar.h>
#include <TGSplitter.h>
#include <TColor.h>
#include <THtml.h>

#include <TPluginManager.h>
#include <TVirtualGL.h>

#include "constants.h"

enum RootShowerMessageTypes {
   M_FILE_OPEN,
   M_FILE_SAVEAS,
   M_FILE_HTML,
   M_FILE_EXIT,
   M_EVENT_NEXT,
   M_EVENT_SELECT,
   M_INTERRUPT_SIMUL,
   M_ZOOM_PLUS,
   M_ZOOM_MOINS,
   M_ZOOM_PLUS2,
   M_ZOOM_MOINS2,

   M_SHOW_PROCESS,
   M_ANIMATE_GIF,
   M_SETTINGS_DLG,
   M_SETTINGS_SAVE,
   M_SHOW_INFOS,
   M_SHOW_3D,
   M_SHOW_TRACK,

   M_VIEW_TOOLBAR,
   M_INSPECT_BROWSER,

   M_HELP_PHYSICS,
   M_HELP_SIMULATION,
   M_HELP_LICENSE,
   M_HELP_ABOUT
};

const char *xpm_names[] = {
   "open.xpm",
   "save.xpm",
   "",
   "settings.xpm",
   "",
   "infos.xpm",
   "view3d.xpm",
   "",
   "browser.xpm",
   "",
   "manual.xpm",
   "help.xpm",
   "license.xpm",
   "about.xpm",
   "",
   "quit.xpm",
   0
};

ToolBarData_t tb_data[] = {
   { "", "Open Root event file",     kFALSE, M_FILE_OPEN,        NULL },
   { "", "Save event in Root file",  kFALSE, M_FILE_SAVEAS,      NULL },
   { "",              0,             0,      -1,                 NULL },
   { "", "Event settings",           kFALSE, M_SETTINGS_DLG,     NULL },
   { "",              0,             0,      -1,                 NULL },
   { "", "Infos on current event",   kFALSE, M_SHOW_INFOS,       NULL },
   { "", "Open 3D viewer",           kFALSE, M_SHOW_3D,          NULL },
   { "",              0,             0,      -1,                 NULL },
   { "", "Start Root browser",       kFALSE, M_INSPECT_BROWSER,  NULL },
   { "",              0,             0,      -1,                 NULL },
   { "", "Physics recalls",          kFALSE, M_HELP_PHYSICS,     NULL },
   { "", "RootShower help",          kFALSE, M_HELP_SIMULATION,  NULL },
   { "", "Display license",          kFALSE, M_HELP_LICENSE,     NULL },
   { "", "About RootShower",         kFALSE, M_HELP_ABOUT,       NULL },
   { "",              0,             0,      -1,                 NULL },
   { "", "Exit Application",         kFALSE, M_FILE_EXIT,        NULL },
   { NULL,            NULL,          0,      0,                  NULL }
};

RootShower      *gRootShower;
Int_t            gColIndex;
TGListTree      *gEventListTree; // event selection TGListTree
TGListTreeItem  *gBaseLTI;
TGListTreeItem  *gTmpLTI;
TGListTreeItem  *gLTI[MAX_PARTICLE];

const TGPicture *bpic, *bspic;
const TGPicture *lpic, *lspic;

const Char_t *filetypes[] = {
   "ROOT files",    "*.root",
   "ROOT macros",   "*.C",
   "GIF  files",    "*.gif",
   "PS   files",    "*.ps",
   "EPS  files",    "*.eps",
   "All files",     "*",
   0,               0
};

enum EGeometrySettingsDialogMessageTypes {
   kM_BUTTON_OK,
   kM_BUTTON_CANCEL,
   kM_COMBOBOX_CHANNELID,
   kM_COMBOBOX_TDC
};


//_________________________________________________
// RootShower
//

Int_t RootShower::fgDefaultXPosition = 20;
Int_t RootShower::fgDefaultYPosition = 20;


////////////////////////////////////////////////////////////////////////////////
/// Create (the) Event Display.
///
/// p = pointer to GMainFrame (not owner)
/// w = width of RootShower frame
/// h = width of RootShower frame

RootShower::RootShower(const TGWindow *p, UInt_t w, UInt_t h):
  TGMainFrame(p, w, h)
{
   fOk                 = kFALSE;
   fModified           = kFALSE;
   fSettingsModified   = kFALSE;
   fIsRunning          = kFALSE;
   fShowProcess        = kFALSE;
   fCreateGIFs         = kFALSE;
   fTimer              = 0;
   fPicIndex           = 1;

   fRootShowerEnv = new TEnv(".rootshowerrc");

   fFirstParticle = fRootShowerEnv->GetValue("RootShower.fFirstParticle", PHOTON);
   fE0            = fRootShowerEnv->GetValue("RootShower.fE0", 10.0);
   fB             = fRootShowerEnv->GetValue("RootShower.fB", 20.000);
   fPicNumber     = fRootShowerEnv->GetValue("RootShower.fPicNumber", 24);
   fPicDelay      = fRootShowerEnv->GetValue("RootShower.fPicDelay", 100);
   fPicReset      = fRootShowerEnv->GetValue("RootShower.fPicReset", 1);

   fEventNr = 0;
   fNRun    = 0;

   bpic = gClient->GetPicture("branch_t.xpm");
   bspic = gClient->GetPicture("branch_t.xpm");

   lpic = gClient->GetPicture("leaf_t.xpm");
   lspic = gClient->GetPicture("leaf_t.xpm");

   // Create menubar and popup menus.
   MakeMenuBarFrame();

   //---- toolbar

   int spacing = 8;
   fToolBar = new TGToolBar(this, 60, 20, kHorizontalFrame | kRaisedFrame);
   for (int i = 0; xpm_names[i]; i++) {
      TString iconname(gProgPath);
#ifdef R__WIN32
      iconname += "\\icons\\";
#else
      iconname += "/icons/";
#endif
      iconname += xpm_names[i];
      tb_data[i].fPixmap = iconname.Data();
      if (strlen(xpm_names[i]) == 0) {
         fToolBar->AddFrame(new TGVertical3DLine(fToolBar), new TGLayoutHints(kLHintsExpandY, 4, 4));
         continue;
      }
      const TGPicture *pic = fClient->GetPicture(tb_data[i].fPixmap);
      TGPictureButton *pb = new TGPictureButton(fToolBar, pic, tb_data[i].fId);
      pb->SetToolTipText(tb_data[i].fTipText);
      tb_data[i].fButton = pb;

      fToolBar->AddButton(this, pb, spacing);
      spacing = 0;
   }
   AddFrame(fToolBar, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 0, 0, 0, 0));
   fToolBar->GetButton(M_SHOW_3D)->SetState(kButtonDisabled);
   fToolBar->GetButton(M_FILE_SAVEAS)->SetState(kButtonDisabled);

   // Layout hints
   fL1 = new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 0, 0, 0);
   fL2 = new TGLayoutHints(kLHintsCenterX | kLHintsExpandX, 0, 0, 0, 0);
   fL3 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY,
                           0, 0, 0, 0);
   fL4 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandY, 5, 5, 2, 2);
   fL5 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY,
                           2, 2, 2, 2);
   fL6 = new TGLayoutHints(kLHintsBottom| kLHintsExpandX, 0, 0, 0, 0);
   fL7 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY,
                           5, 5, 2, 2);
   fL8 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX, 0, 0, 0, 0);

   // CREATE TITLE FRAME
   fTitleFrame = new GTitleFrame(this, "ROOT Shower Monte Carlo", "Event Display", 100, 100);
   AddFrame(fTitleFrame, fL2);

   // CREATE MAIN FRAME
   fMainFrame = new TGCompositeFrame(this, 100, 100, kHorizontalFrame | kRaisedFrame);

   TGVerticalFrame *fV1 = new TGVerticalFrame(fMainFrame, 150, 10, kSunkenFrame | kFixedWidth);
   TGVerticalFrame *fV2 = new TGVerticalFrame(fMainFrame, 10, 10, kSunkenFrame);

   TGLayoutHints *lo;

   lo = new TGLayoutHints(kLHintsLeft | kLHintsExpandY,2,0,2,2);
   fMainFrame->AddFrame(fV1, lo);

   TGVSplitter *splitter = new TGVSplitter(fMainFrame, 5);
   splitter->SetFrame(fV1, kTRUE);
   lo = new TGLayoutHints(kLHintsLeft | kLHintsExpandY, 0, 0 ,0, 0);
   fMainFrame->AddFrame(splitter, lo);

   lo = new TGLayoutHints(kLHintsRight | kLHintsExpandX | kLHintsExpandY,0,2,2,2);
   fMainFrame->AddFrame(fV2, lo);


   // Create Selection frame (i.e. with buttons and geometry selection widgets)
   fSelectionFrame = new TGCompositeFrame(fV1, 100, 100, kVerticalFrame);
   // create button frame
   fButtonFrame = new GButtonFrame (fSelectionFrame, this, M_EVENT_NEXT,
                                    M_EVENT_SELECT, M_INTERRUPT_SIMUL);
   lo = new TGLayoutHints(kLHintsTop | kLHintsCenterX | kLHintsExpandX, 2, 5, 1, 2);
   fSelectionFrame->AddFrame(fButtonFrame, lo);

   fTreeView = new TGCanvas(fSelectionFrame, 150, 10, kSunkenFrame | kDoubleBorder);
   fEventListTree = new TGListTree(fTreeView->GetViewPort(), 10, 10, kHorizontalFrame);
   gEventListTree = fEventListTree;
   fEventListTree->SetCanvas(fTreeView);
   fEventListTree->Associate(this);
   BuildEventTree();
   fTreeView->SetContainer(fEventListTree);
   fSelectionFrame->AddFrame(fTreeView, fL5);

   lo = new TGLayoutHints(kLHintsExpandX | kLHintsExpandY);
   fV1->AddFrame(fSelectionFrame, lo);

   fContextMenu = new TContextMenu("RSContextMenu");

   //__________________________________________________________________________________

   // Create Display frame
   fDisplayFrame = new TGTab(fV2, 580, 360);

   // Create Display Canvas Tab (where the actual main event is displayed)
   TGCompositeFrame *tFrame = fDisplayFrame->AddTab("Main Event (Shower)");

   // Create Layout hints
   fZoomButtonsLayout = new TGLayoutHints(kLHintsBottom | kLHintsLeft | kLHintsExpandX, 5, 2, 2, 2);

   fHFrame = new TGHorizontalFrame(tFrame,0,0,0);
   tFrame->AddFrame(fHFrame, new TGLayoutHints(kLHintsBottom | kLHintsLeft | kLHintsExpandX, 5, 5, 5, 5));
   // Create Zoom Buttons
   fZoomPlusButton = new TGTextButton(fHFrame, "&Zoom Forward", M_ZOOM_PLUS);
   fZoomPlusButton->Associate(this);
   fZoomPlusButton->SetToolTipText("Zoom forward event view");
   fHFrame->AddFrame(fZoomPlusButton, fZoomButtonsLayout);
   fZoomMoinsButton = new TGTextButton(fHFrame, "Zoom &Backward", M_ZOOM_MOINS);
   fZoomMoinsButton->Associate(this);
   fZoomMoinsButton->SetToolTipText("Zoom backward event view");
   fHFrame->AddFrame(fZoomMoinsButton, fZoomButtonsLayout);

   fEmbeddedCanvas = new TRootEmbeddedCanvas("fEmbeddedCanvas", tFrame, 580, 360);
   tFrame->AddFrame(fEmbeddedCanvas, fL5);
   fEmbeddedCanvas->GetCanvas()->SetBorderMode(0);
   fCA = fEmbeddedCanvas->GetCanvas();
   fCA->SetFillColor(1);

   // Create Display Canvas Tab (where the selected event is displayed)
   TGCompositeFrame *tFrame2 = fDisplayFrame->AddTab("Selected Track");

   fHFrame2 = new TGHorizontalFrame(tFrame2,0,0,0);
   tFrame2->AddFrame(fHFrame2, new TGLayoutHints(kLHintsBottom | kLHintsLeft | kLHintsExpandX, 5, 5, 5, 5));
   // Create Zoom Buttons
   fZoomPlusButton2 = new TGTextButton(fHFrame2, "&Zoom Forward", M_ZOOM_PLUS2);
   fZoomPlusButton2->Associate(this);
   fZoomPlusButton2->SetToolTipText("Zoom forward event view");
   fHFrame2->AddFrame(fZoomPlusButton2, fZoomButtonsLayout);
   fZoomMoinsButton2 = new TGTextButton(fHFrame2, "Zoom &Backward", M_ZOOM_MOINS2);
   fZoomMoinsButton2->Associate(this);
   fZoomMoinsButton2->SetToolTipText("Zoom backward event view");
   fHFrame2->AddFrame(fZoomMoinsButton2, fZoomButtonsLayout);

   fEmbeddedCanvas2 = new TRootEmbeddedCanvas("fEmbeddedCanvas2", tFrame2, 580, 360);
   tFrame2->AddFrame(fEmbeddedCanvas2, fL5);
   fEmbeddedCanvas2->GetCanvas()->SetBorderMode(0);
   fCB = fEmbeddedCanvas2->GetCanvas();
   fCB->SetFillColor(1);

   // Create Display Canvas Tab (where the histogram is displayed)
   TGCompositeFrame *tFrame3 = fDisplayFrame->AddTab("Statistics");

   fEmbeddedCanvas3 = new TRootEmbeddedCanvas("fEmbeddedCanvas3", tFrame3, 580, 360);
   tFrame3->AddFrame(fEmbeddedCanvas3, fL5);
   fEmbeddedCanvas3->GetCanvas()->SetBorderMode(0);
   fCC = fEmbeddedCanvas3->GetCanvas();
   fCC->SetFillColor(10);
   fCC->cd();
   fPadC = new TPad("fPadC","Histogram",0.0,0.0,1.0,1.0,10,3,1);
   fPadC->SetFillColor(10);
   fPadC->SetBorderMode(0);
   fPadC->SetBorderSize(0);
   fPadC->Draw();
   // Creation of histogram for particle's energy loss
   fHisto_dEdX = new TH1F("Statistics","Energy loss for each particle",100,0,0.025); // Max at 25 MeV
   fHisto_dEdX->SetFillColor(38);
   fHisto_dEdX->SetStats(kTRUE);
   fHisto_dEdX->SetXTitle("Energy Loss [GeV]");
   fHisto_dEdX->SetLabelFont(42,"X");
   fHisto_dEdX->SetLabelSize(0.03f, "X");
   fHisto_dEdX->GetXaxis()->SetTitleFont(42);
   fHisto_dEdX->SetYTitle("Number");
   fHisto_dEdX->SetLabelFont(42,"Y");
   fHisto_dEdX->SetLabelSize(0.03f, "Y");
   fHisto_dEdX->GetYaxis()->SetTitleFont(42);

   fCC->Update();

   // Create text display Tab
   tFrame = fDisplayFrame->AddTab("PDG Table");
   fTextView = new TGTextEdit(tFrame, 300, 100, kSunkenFrame | kDoubleBorder);
   tFrame->AddFrame(fTextView, fL5);
   TString pdgFilename = gSystem->Getenv("ROOTSYS");
   pdgFilename.Append("/etc/pdg_table.txt");

   fTextView->LoadFile(pdgFilename);

   lo = new TGLayoutHints(kLHintsExpandX | kLHintsExpandY);
   fV2->AddFrame(fDisplayFrame, lo);

   AddFrame(fMainFrame, lo);

   // Create status bar
   Int_t parts[] = {45, 45, 10};
   fStatusBar = new TGStatusBar(this, 50, 10, kHorizontalFrame);
   fStatusBar->SetParts(parts, 3);
   AddFrame(fStatusBar, fL6);
   fStatusBar->SetText("Waiting to start simulation...",0);

   // Finish RootShower for display...
   SetWindowName("Root Shower Event Display");
   SetIconName("Root Shower Event Display");
   MapSubwindows();
   Resize(GetDefaultSize()); // this is used here to init layout algoritme
   MapWindow();
   fEvent = new MyEvent();
   fEvent->GetDetector()->Init();
   fEvent->Init(0, fFirstParticle, fE0, fB);
   Initialize(1);
   gROOT->GetListOfBrowsables()->Add(fEvent,"RootShower Event");
   gSystem->Load("libTreeViewer");
   AddInput(kKeyPressMask | kKeyReleaseMask);
   gVirtualX->SetInputFocus(GetId());
   gRootShower = this;
}


////////////////////////////////////////////////////////////////////////////////
/// Create menubar and popup menus.

void RootShower::MakeMenuBarFrame()
{
   // layout hint items
   fMenuBarLayout = new TGLayoutHints(kLHintsTop| kLHintsLeft | kLHintsExpandX,
                                      0, 0, 0, 0);
   fMenuBarItemLayout = new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0);
   fMenuBarHelpLayout = new TGLayoutHints(kLHintsTop | kLHintsRight);

   fMenuBar = new TGMenuBar(this, 1, 1, kHorizontalFrame | kRaisedFrame);

   // file popup menu
   fMenuFile = new TGPopupMenu(gClient->GetRoot());
   fMenuFile->AddEntry("&Open...\tCtrl+O", M_FILE_OPEN);
   fMenuFile->AddEntry("S&ave as...\tCtrl+A", M_FILE_SAVEAS);
   fMenuFile->AddEntry("&Close", -1);
   fMenuFile->AddSeparator();
   fMenuFile->AddEntry("&Print", -1);
   fMenuFile->AddEntry("P&rint setup...", -1);
   fMenuFile->AddSeparator();
   fMenuFile->AddEntry("E&xit\tCtrl+Q", M_FILE_EXIT);
   fMenuFile->DisableEntry(M_FILE_SAVEAS);
   fMenuFile->Associate(this);

   fMenuEvent = new TGPopupMenu(gClient->GetRoot());
   fMenuEvent->AddLabel("Simulation Settings...");
   fMenuEvent->AddSeparator();
   fMenuEvent->AddEntry("&Settings...", M_SETTINGS_DLG);
   fMenuEvent->AddEntry("Save &Parameters", M_SETTINGS_SAVE);
   fMenuEvent->AddEntry("Show &Process", M_SHOW_PROCESS);
   fMenuEvent->AddEntry("Animated &GIF", M_ANIMATE_GIF);
   fMenuEvent->AddEntry("&Infos...\tCtrl+I", M_SHOW_INFOS);
   fMenuEvent->AddSeparator();
   fMenuEvent->AddEntry("&3D View", M_SHOW_3D);
   fMenuEvent->AddEntry("&Show Selected  Track", M_SHOW_TRACK);
   fMenuEvent->DisableEntry(M_SHOW_INFOS);
   fMenuEvent->DisableEntry(M_SHOW_3D);
   fMenuEvent->DisableEntry(M_SHOW_TRACK);

   fMenuEvent->DisableEntry(M_SHOW_PROCESS);
   fMenuEvent->DisableEntry(M_ANIMATE_GIF);

   fMenuEvent->Associate(this);

   fMenuTools = new TGPopupMenu(gClient->GetRoot());
   fMenuTools->AddLabel("Simulation Tools...");
   fMenuTools->AddSeparator();
   fMenuTools->AddEntry("Start &Browser\tCtrl+B", M_INSPECT_BROWSER);
   fMenuTools->AddEntry("&Create Html Doc", M_FILE_HTML);
   fMenuTools->Associate(this);

   fMenuView = new TGPopupMenu(gClient->GetRoot());
   fMenuView->AddEntry("&Toolbar", M_VIEW_TOOLBAR);
   fMenuView->Associate(this);
   fMenuView->CheckEntry(M_VIEW_TOOLBAR);

   fMenuHelp = new TGPopupMenu(gClient->GetRoot());
   fMenuHelp->AddEntry("&Physics", M_HELP_PHYSICS);
   fMenuHelp->AddEntry("&Simulation", M_HELP_SIMULATION);
   fMenuHelp->AddSeparator();
   fMenuHelp->AddEntry("&License...", M_HELP_LICENSE);
   fMenuHelp->AddEntry("&About...", M_HELP_ABOUT);
   fMenuHelp->Associate(this);

   fMenuBar->AddPopup("&File", fMenuFile, fMenuBarItemLayout);
   fMenuBar->AddPopup("&Event", fMenuEvent, fMenuBarItemLayout);
   fMenuBar->AddPopup("&Tools", fMenuTools, fMenuBarItemLayout);
   fMenuBar->AddPopup("&View", fMenuView, fMenuBarItemLayout);
   fMenuBar->AddPopup("&Help", fMenuHelp, fMenuBarHelpLayout);

   AddFrame(fMenuBar, fMenuBarLayout);
}


////////////////////////////////////////////////////////////////////////////////
/// Destroy menubar and popup menus.

void RootShower::CloseMenuBarFrame()
{
   delete fMenuHelp;
   delete fMenuEvent;
   delete fMenuTools;
   delete fMenuFile;

   delete fMenuBarItemLayout;
   delete fMenuBarHelpLayout;
   delete fMenuBar;
   delete fMenuBarLayout;
}

////////////////////////////////////////////////////////////////////////////////
/// Show or hide toolbar.

void RootShower::ShowToolBar(Bool_t show)
{
   if (show) {
      ShowFrame(fToolBar);
      fMenuView->CheckEntry(M_VIEW_TOOLBAR);
   } else {
      HideFrame(fToolBar);
      fMenuView->UnCheckEntry(M_VIEW_TOOLBAR);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Destroy RootShower object. Delete all created widgets
/// GUI MEMBERS

RootShower::~RootShower()
{
   CloseMenuBarFrame();

   delete fContextMenu;
   delete fZoomPlusButton2;
   delete fZoomMoinsButton2;
   delete fZoomPlusButton;
   delete fZoomMoinsButton;
   delete fHFrame;
   delete fHFrame2;
   delete fZoomButtonsLayout;

   delete fEmbeddedCanvas;
   delete fTextView;
   delete fDisplayFrame;
   delete fEventListTree;
   delete fTreeView;
   delete fButtonFrame;
   delete fSelectionFrame;
   delete fMainFrame;
   delete fTitleFrame;

   delete fL8;
   delete fL7;
   delete fL6;
   delete fL5;
   delete fL4;
   delete fL3;
   delete fL2;
   delete fL1;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the default position on the screen of new RootShower instances.

void RootShower::setDefaultPosition(Int_t x, Int_t y)
{
   fgDefaultXPosition = x;
   fgDefaultYPosition = y;
}

////////////////////////////////////////////////////////////////////////////////
/// Apply layout on the main frame.

void RootShower::Layout()
{
   TGMainFrame::Layout();
}


////////////////////////////////////////////////////////////////////////////////
/// Got close message for this RootShower. The EventDislay and the
/// application will be terminated.

void RootShower::CloseWindow()
{
   if (fModified) {
      new RootShowerMsgBox(gClient->GetRoot(),this, 400, 200);
      if ( fOk ) {
         fRootShowerEnv->SetValue("RootShower.fFirstParticle",fFirstParticle);
         fRootShowerEnv->SetValue("RootShower.fE0",fE0);
         fRootShowerEnv->SetValue("RootShower.fB",fB);
         fRootShowerEnv->SaveLevel(kEnvLocal);
         std::cout << " Saving stuff .... " << std::endl;
#ifdef R__WIN32
         gSystem->Exec("del .rootshowerrc");
         gSystem->Rename(".rootshowerrc.new",".rootshowerrc");
#endif
      }
   }
   std::cout << "Terminating RootShower" << std::endl;
   DeleteWindow();
   gApplication->Terminate(0);
}

////////////////////////////////////////////////////////////////////////////////
/// This event is generated when the frame is resized.

Bool_t RootShower::HandleConfigureNotify(Event_t *event)
{
   TGFrame* f = (TGFrame*)this;
   if ((event->fWidth != f->GetWidth()) || (event->fHeight != f->GetHeight())) {
      UInt_t w = event->fWidth;
      UInt_t h = event->fHeight;
      f->Resize(w,h);
      f->Layout();
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle messages send to the RootShower object.

Bool_t RootShower::ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2)
{
   Window_t wdummy;
   int ax, ay;
   TRootHelpDialog *hd;
   TGListTreeItem *item;
   TGFileInfo fi;
   Char_t  strtmp[250];

   switch (GET_MSG(msg)) {

      case kC_COMMAND:

         switch (GET_SUBMSG(msg)) {

            case kCM_BUTTON:
            case kCM_MENU:
               switch (parm1) {

                  case M_EVENT_NEXT:
                     if (fDisplayFrame->GetCurrent() != 0)
                        fDisplayFrame->SetTab(0);
                     Initialize(0);
                     fStatusBar->SetText("Simulation running, please wait...",0);
                     fButtonFrame->SetState(GButtonFrame::kNoneActive);
                     fMenuEvent->DisableEntry(M_SETTINGS_DLG);
                     OnShowerProduce();
                     fEventListTree->ClearViewPort();
                     fClient->NeedRedraw(fEventListTree);
                     fButtonFrame->SetState(GButtonFrame::kAllActive);
                     fMenuEvent->EnableEntry(M_SETTINGS_DLG);
                     sprintf(strtmp,"Done - Total particles : %d - Waiting for next simulation",
                             fEvent->GetTotal());
                     fStatusBar->SetText(strtmp,0);
                     break;
                  case M_EVENT_SELECT:
                     if (fDisplayFrame->GetCurrent() != 1)
                        fDisplayFrame->SetTab(1);
                     if ((item = fEventListTree->GetSelected()) != 0)
                        OnShowSelected(item);
                     break;
                  case M_INTERRUPT_SIMUL:
                     Interrupt();
                     break;

                  case M_ZOOM_PLUS:
                     fCA->cd();
                     fCA->GetView()->ZoomView(0, 1.25);
                     fCA->Modified();
                     fCA->Update();
                     break;

                  case M_ZOOM_MOINS:
                     fCA->cd();
                     fCA->GetView()->UnzoomView(0, 1.25);
                     fCA->Modified();
                     fCA->Update();
                     break;

                  case M_ZOOM_PLUS2:
                     fCB->cd();
                     fCB->GetView()->ZoomView(0, 1.25);
                     fCB->Modified();
                     fCB->Update();
                     break;

                  case M_ZOOM_MOINS2:
                     fCB->cd();
                     fCB->GetView()->UnzoomView(0, 1.25);
                     fCB->Modified();
                     fCB->Update();
                     break;

                  case M_FILE_OPEN:
                     if (fIsRunning) break;
                     fi.fFileTypes = filetypes;
                     new TGFileDialog(fClient->GetRoot(), this, kFDOpen,&fi);
                     if (!fi.fFilename) return kTRUE;
                     OnOpenFile(fi.fFilename);
                     break;

                  case M_FILE_HTML:
                     {
                        THtml html;
                        html.SetInputDir(gProgPath);
                        html.MakeClass("MyParticle");
                        html.MakeClass("MyDetector");
                        html.MakeClass("EventHeader");
                        html.MakeClass("MyEvent");
                        html.MakeIndex();
                     }
                     break;

                  case M_FILE_SAVEAS:
                     if (fIsRunning) break;
                     fi.fFileTypes = filetypes;
                     new TGFileDialog(fClient->GetRoot(), this, kFDSave,&fi);
                     if (!fi.fFilename) return kTRUE;
                     OnSaveFile(fi.fFilename);
                     break;

                  case M_FILE_EXIT:
                     CloseWindow();   // this also terminates theApp
                     break;

                  case M_SHOW_PROCESS:
                     if (fShowProcess) {
                        fMenuEvent->UnCheckEntry(M_SHOW_PROCESS);
                        fShowProcess = kFALSE;
                     }
                     else {
                        fMenuEvent->CheckEntry(M_SHOW_PROCESS);
                        fShowProcess = kTRUE;
                     }
                     break;

                  case M_ANIMATE_GIF:
                     if (fCreateGIFs) {
                        fMenuEvent->UnCheckEntry(M_ANIMATE_GIF);
                        fCreateGIFs = kFALSE;
                     }
                     else {
                        fMenuEvent->CheckEntry(M_ANIMATE_GIF);
                        fCreateGIFs = kTRUE;
                     }
                     break;

                  case M_SETTINGS_DLG:
                     if (fIsRunning) break;
                     new SettingsDialog(fClient->GetRoot(), this, 400, 200);
                     if (fSettingsModified) {
                        fEvent->Init(0, fFirstParticle, fE0, fB);
                        Initialize(0);
                        gRootShower->Modified();
                        gRootShower->SettingsModified(kFALSE);
                     }
                     break;

                  case M_SETTINGS_SAVE:
                     fRootShowerEnv->SetValue("RootShower.fFirstParticle",fFirstParticle);
                     fRootShowerEnv->SetValue("RootShower.fE0",fE0);
                     fRootShowerEnv->SetValue("RootShower.fB",fB);
                     fRootShowerEnv->SaveLevel(kEnvLocal);
#ifdef R__WIN32
                     gSystem->Exec("del .rootshowerrc");
                     gSystem->Rename(".rootshowerrc.new",".rootshowerrc");
#endif
                     gRootShower->Modified(kFALSE);
                     break;

                  case M_SHOW_INFOS:
                     if (fIsRunning) break;
                     ShowInfos();
                     break;

                  case M_INSPECT_BROWSER:
                     new TBrowser;
                     break;

                  case M_VIEW_TOOLBAR:
                     if (fMenuView->IsEntryChecked(M_VIEW_TOOLBAR))
                        ShowToolBar(kFALSE);
                     else
                        ShowToolBar();
                     break;

                  case M_HELP_PHYSICS:
#ifdef R__WIN32
                     sprintf(strtmp, "start winhlp32 %s\\Physics.hlp",gProgPath);
                     gSystem->Exec(strtmp);
#else
                     sprintf(strtmp, "Help on Physics");
                     hd = new TRootHelpDialog(this, strtmp, 620, 350);
                     hd->SetText(gPhysicsHelpText);
                     gVirtualX->TranslateCoordinates(GetId(), GetParent()->GetId(),
                                                     (Int_t)(GetWidth() - 620) >> 1,
                                                     (Int_t)(GetHeight() - 350) >> 1,
                                                     ax, ay, wdummy);
                     hd->Move(ax, ay);
                     hd->Popup();
                     fClient->WaitFor(hd);
#endif
                     break;

                  case M_HELP_SIMULATION:
                     sprintf(strtmp, "Help on Simulation");
                     hd = new TRootHelpDialog(this, strtmp, 620, 350);
                     hd->SetText(gSimulationHelpText);
                     gVirtualX->TranslateCoordinates(GetId(), GetParent()->GetId(),
                                                     (Int_t)(GetWidth() - 620) >> 1,
                                                     (Int_t)(GetHeight() - 350) >> 1,
                                                     ax, ay, wdummy);
                     hd->Move(ax, ay);
                     hd->Popup();
                     fClient->WaitFor(hd);
                     break;

                  case M_HELP_LICENSE:
                     sprintf(strtmp, "RootShower License");
                     hd = new TRootHelpDialog(this, strtmp, 640, 380);
                     hd->SetText(gHelpLicense);
                     gVirtualX->TranslateCoordinates(GetId(), GetParent()->GetId(),
                                                    (Int_t)(GetWidth() - 640) >> 1,
                                                    (Int_t)(GetHeight() - 380) >> 1,
                                                    ax, ay, wdummy);
                     hd->Move(ax, ay);
                     hd->Popup();
                     fClient->WaitFor(hd);
                     break;

                  case M_HELP_ABOUT:
                     new RootShowerAbout(gClient->GetRoot(),this, 400, 200);
                     break;

                  case M_SHOW_3D:
                     if (gGLManager) {
                        if (fIsRunning) break;
                        TString cmd;
                        fCA->cd();
                        TVirtualViewer3D *viewer3D = fCA->GetViewer3D("ogl");
                        cmd.Form("((TGLViewer *)0x%lx)->SetCurrentCamera(TGLViewer::kCameraPerspXOY)",
                                 (ULong_t)viewer3D);
                        gROOT->ProcessLine(cmd);
                        cmd.Form("((TGLViewer *)0x%lx)->CurrentCamera().RotateRad(0.0, TMath::Pi())",
                                 (ULong_t)viewer3D);
                        gROOT->ProcessLine(cmd);
                        cmd.Form("((TGLViewer *)0x%lx)->CurrentCamera().Dolly(-100, 0, 0)",
                                 (ULong_t)viewer3D);
                        gROOT->ProcessLine(cmd);
                     }
                     break;

                  case M_SHOW_TRACK:
                     if (fIsRunning) break;
                     if (fDisplayFrame->GetCurrent() != 1)
                        fDisplayFrame->SetTab(1);
                     if ((item = fEventListTree->GetSelected()) != 0)
                        OnShowSelected(item);
                     break;


               } // switch parm1
               break; // M_MENU

            } // switch submsg
            break; // case kC_COMMAND

         case kC_LISTTREE:

            switch (GET_SUBMSG(msg)) {

               case kCT_ITEMDBLCLICK:
                  if (parm1 == kButton1) {
                     if (fEventListTree->GetSelected()) {
                        fEventListTree->ClearViewPort();
                        fClient->NeedRedraw(fEventListTree);
                     }
                  }
                  break;

               case kCT_ITEMCLICK:
                  if (parm1 == kButton3) {
                     if (fEventListTree->GetSelected()) {
                        Int_t x = (Int_t)(parm2 & 0xffff);
                        Int_t y = (Int_t)((parm2 >> 16) & 0xffff);
                        Clicked(fEventListTree->GetSelected(), x, y);
                     }
                  }
                  break;

            } // switch submsg
            break; // case kC_LISTTREE
   } // switch msg

   return kTRUE;
}


////////////////////////////////////////////////////////////////////////////////
/// Add item to the TGListTree of the event display. It will be connected
/// to the current TGListTreeItem (i.e. fCurEventListItem)

TGListTreeItem* RootShower::AddToTree(const char *name)
{
   TGListTreeItem *e = 0;
   e = fEventListTree->AddItem(fCurListItem, name);
   return e;
}

////////////////////////////////////////////////////////////////////////////////
/// Add recursively stations and layers (and cells) in TGListTree.

void RootShower::BuildEventTree()
{
   fCurListItem = 0;
   TGListTreeItem *eventLTItem = AddToTree("Event");
   fCurListItem = eventLTItem;
   gBaseLTI = eventLTItem;
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize RootShower display.

void RootShower::Initialize(Int_t set_angles)
{
   Interrupt(kFALSE);
   fEventListTree->DeleteChildren(fCurListItem);
   fEventListTree->ClearViewPort();
   fClient->NeedRedraw(fEventListTree);

   fCB->cd();
   fCB->SetFillColor(1);
   fCB->Clear();
   gGeoManager->GetTopVolume()->Draw();
   fCB->GetView()->SetPerspective();
   if (set_angles)
      fCB->GetView()->SideView();
   gGeoManager->GetTopVolume()->Draw();
   fCB->cd();
   fCB->Update();

   fCA->cd();
   fCA->SetFillColor(1);
   fCA->Clear();
   gGeoManager->GetTopVolume()->Draw();
   fCA->GetView()->SetPerspective();
   if (set_angles)
      fCA->GetView()->SideView();
   gGeoManager->GetTopVolume()->Draw();
   fCA->cd();
   fCA->Update();
   fStatusBar->SetText("",1);
}

////////////////////////////////////////////////////////////////////////////////
/// Produce (generate) one event.

void RootShower::Produce()
{
   Int_t     local_num,local_last,local_end;
   Int_t     old_num;
   Bool_t    first_pass;
   Char_t    strtmp[80];

   // Check if some Event parameters have changed
   if ((fEvent->GetHeader()->GetDate() != fEventTime) ||
       (fEvent->GetHeader()->GetPrimary() != fFirstParticle) ||
       (fEvent->GetHeader()->GetEnergy() != fE0) ||
       (fEvent->GetB() != fB)) {
      fEventNr++;
      fNRun = 0;
   }
   fEvent->SetHeader(fEventNr, fNRun++, fEventTime, fFirstParticle, fE0);
   fEvent->Init(0, fFirstParticle, fE0, fB);

   fMenuFile->DisableEntry(M_FILE_SAVEAS);
   fMenuEvent->DisableEntry(M_SHOW_3D);
   fToolBar->GetButton(M_SHOW_3D)->SetState(kButtonDisabled);
   fToolBar->GetButton(M_FILE_SAVEAS)->SetState(kButtonDisabled);
   Interrupt(kFALSE);
   first_pass = kTRUE;
   old_num = -1;
   // loop events until user interrupt or until all particles are dead
   while ((!IsInterrupted()) && (fEvent->GetNAlives() > 0)) {
      if (first_pass && fEvent->GetTotal() > 1) {
         fEventListTree->OpenItem(gBaseLTI);
         fEventListTree->OpenItem(gLTI[0]);
         fEventListTree->ClearViewPort();
         fClient->NeedRedraw(fEventListTree);
         first_pass = kFALSE;
      }
      if (fEvent->GetTotal() > old_num) {
         sprintf(strtmp,"Simulation running, particles : %4d, please wait...",fEvent->GetTotal());
         old_num = fEvent->GetTotal();
         fStatusBar->SetText(strtmp,0);
         // Update display here to not slow down too much...
         gSystem->ProcessEvents();
      }
      local_last = fEvent->GetLast();
      local_num = 0;
      local_end = kFALSE;
      while ((!IsInterrupted()) && (local_end == kFALSE) && (local_num < (local_last + 1))) {
         // Update display here if fast machine...
         if (fEvent->GetParticle(local_num)->GetStatus() != DEAD) {
            gSystem->ProcessEvents();
            if (fEvent->Action(local_num) == DEAD)
               local_end = kTRUE;
            if (fEvent->GetParticle(local_num)->GetStatus() == CREATED)
               fEvent->GetParticle(local_num)->SetStatus(ALIVE);
         }
         local_num ++;
      }
   }
   fMenuEvent->EnableEntry(M_SHOW_INFOS);
   if (!IsInterrupted()) {
      if (gGLManager) {
         fMenuEvent->EnableEntry(M_SHOW_3D);
         fToolBar->GetButton(M_SHOW_3D)->SetState(kButtonUp);
      }
      fToolBar->GetButton(M_FILE_SAVEAS)->SetState(kButtonUp);
      fMenuFile->EnableEntry(M_FILE_SAVEAS);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize and generate one event.

void RootShower::OnShowerProduce()
{
   Int_t i, j;
   fStatusBar->SetText("",1);

   SetWindowName("Root Shower Event Display");

   // animation logo handling
   if (fPicReset > 0) fPicIndex = 1;
   // animation timer
   if (!fTimer) fTimer = new TTimer(this, fPicDelay);
   fTimer->Reset();
   fTimer->TurnOn();
   fEventTime.Set();

   fIsRunning = kTRUE;
   fHisto_dEdX->Reset();
   Produce();
   Interrupt(kFALSE);
   for (i=0;i<=fEvent->GetTotal();i++) {
      gSystem->ProcessEvents();  // handle GUI events
      if (IsInterrupted()) break;
      // if particle has no child, represent it by a leaf,
      // otherwise by a branch
      if (fEvent->GetParticle(i)->GetChildId(0) == 0) {
         lpic = gClient->GetPicture("leaf_t.xpm");
         lspic = gClient->GetPicture("leaf_t.xpm");
         gLTI[i]->SetPictures(lpic, lspic);
      }
      else {
         bpic = gClient->GetPicture("branch_t.xpm");
         bspic = gClient->GetPicture("branch_t.xpm");
         gLTI[i]->SetPictures(bpic, bspic);
      }
      // Show only charged and massive particles...
      if ((fEvent->GetParticle(i)->GetPdgCode() != PHOTON) &&
          (fEvent->GetParticle(i)->GetPdgCode() != NEUTRINO_E) &&
          (fEvent->GetParticle(i)->GetPdgCode() != NEUTRINO_MUON) &&
          (fEvent->GetParticle(i)->GetPdgCode() != NEUTRINO_TAU) &&
          (fEvent->GetParticle(i)->GetPdgCode() != ANTINEUTRINO_E) &&
          (fEvent->GetParticle(i)->GetPdgCode() != ANTINEUTRINO_MUON) &&
          (fEvent->GetParticle(i)->GetPdgCode() != ANTINEUTRINO_TAU) ) {
         // Fill histogram for particle's energy loss
         fHisto_dEdX->Fill(fEvent->GetParticle(i)->GetELoss());
         for (j=0;j<fEvent->GetParticle(i)->GetNTracks();j++)
            fEvent->GetParticle(i)->GetTrack(j)->Draw();
         // show track by track if "show process" has been chosen
         // into the menu
         if (fShowProcess) {
            fCA->Modified();
            fCA->Update();
            // create one gif image by step if "Animated GIF"
            // has been chosen into the menu
            if (fCreateGIFs) {
               fCA->SaveAs("RSEvent.gif+");
            }
         }
      }
   }
   AppendPad();
   fCA->GetView()->SetPerspective();
   fCA->cd();
   fCA->Modified();
   fCA->Update();
   fPadC->cd();
   // do not fit if not enough particles
   if (fHisto_dEdX->GetEntries() > 10) {
      fHisto_dEdX->Fit("landau","L");
      TF1 *f1 = fHisto_dEdX->GetFunction("landau");
      //delete fit function is fit is a non sense
      if (f1 && f1->GetNDF() > 0) {
         f1->SetLineColor(kRed);
         f1->SetLineWidth(1);
      } else {
         delete f1;
      }
   }
   fHisto_dEdX->Draw();
   fPadC->Modified();
   fPadC->Update();
   fCC->Update();
   fPadC->cd();
   fPadC->SetFillColor(16);
   fPadC->GetFrame()->SetFillColor(10);
   fPadC->Draw();
   fPadC->Update();

   // Open first list tree items
   fEventListTree->OpenItem(gBaseLTI);
   fEventListTree->OpenItem(gLTI[0]);
   fTimer->TurnOff();
   fIsRunning = kFALSE;
   if (fPicReset > 0)
      fTitleFrame->ChangeRightLogo(1);
}

////////////////////////////////////////////////////////////////////////////////
/// No comment...

void RootShower::HighLight(TGListTreeItem * /*item*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Shows track which has been selected into the list tree

void RootShower::OnShowSelected(TGListTreeItem *item)
{
   Int_t i, j, retval;

   fCB->cd();
   fCB->SetFillColor(1);
   fCB->SetBorderMode(0);
   fCB->SetBorderSize(0);
   fCB->Clear();
   fCB->cd();
   // draw geometry
   gGeoManager->GetTopVolume()->Draw();
   fCB->GetView()->SetPerspective();
   fCB->cd();
   fCB->Update();
   retval = -1;
   for (i=0;i<=fEvent->GetTotal();i++) {
      if (gLTI[i] == item) {
         retval = i;
         break;
      }
   }
   if ((retval > -1) &&
       (fEvent->GetParticle(i)->GetPdgCode() != PHOTON) &&
       (fEvent->GetParticle(i)->GetPdgCode() != NEUTRINO_E) &&
       (fEvent->GetParticle(i)->GetPdgCode() != NEUTRINO_MUON) &&
       (fEvent->GetParticle(i)->GetPdgCode() != NEUTRINO_TAU) &&
       (fEvent->GetParticle(i)->GetPdgCode() != ANTINEUTRINO_E) &&
       (fEvent->GetParticle(i)->GetPdgCode() != ANTINEUTRINO_MUON) &&
       (fEvent->GetParticle(i)->GetPdgCode() != ANTINEUTRINO_TAU) ) {
      for (j=0;j<fEvent->GetParticle(retval)->GetNTracks();j++)
         fEvent->GetParticle(retval)->GetTrack(j)->Draw();
   }
   fCB->GetView()->SetPerspective();
   fCB->cd();
   fCB->Modified();
   fCB->Update();
}

////////////////////////////////////////////////////////////////////////////////
/// Opens a root file into which a previous event has been saved.

void RootShower::OnOpenFile(const Char_t *filename)
{
   char   strtmp[256];
   Int_t  i,j;
   TFile *f = new TFile(filename);
   TTree *tree;
   TBranch *branch;
   fStatusBar->SetText("",1);

   fEvent->Init(0, fFirstParticle, fE0, fB);
   fHisto_dEdX->Reset();
   tree = (TTree *)f->Get("RootShower");
   if (tree == NULL) return;
   branch = tree->GetBranch("Event");
   branch->SetAddress(&fEvent);
   tree->GetEntry(0, 1);
   f->Close();

   // take back detector dimensions for selection geometry
   gGeoManager->Import(filename, "detector");
   Initialize(1);

   for (i=0;i<=fEvent->GetTotal();i++) {
      gTmpLTI = fEventListTree->AddItem(gBaseLTI, fEvent->GetParticle(i)->GetName());
      gTmpLTI->SetUserData(fEvent->GetParticle(i));
      sprintf(strtmp,"%1.2f GeV",fEvent->GetParticle(i)->Energy());
      fEventListTree->SetToolTipItem(gTmpLTI, strtmp);
      gLTI[i] = gTmpLTI;

      if (fEvent->GetParticle(i)->GetChildId(0) == 0) {
         lpic = gClient->GetPicture("leaf_t.xpm");
         lspic = gClient->GetPicture("leaf_t.xpm");
         gLTI[i]->SetPictures(lpic, lspic);
      }
      else {
         bpic = gClient->GetPicture("branch_t.xpm");
         bspic = gClient->GetPicture("branch_t.xpm");
         gLTI[i]->SetPictures(bpic, bspic);
      }

      if ((fEvent->GetParticle(i)->GetPdgCode() != PHOTON) &&
          (fEvent->GetParticle(i)->GetPdgCode() != NEUTRINO_E) &&
          (fEvent->GetParticle(i)->GetPdgCode() != NEUTRINO_MUON) &&
          (fEvent->GetParticle(i)->GetPdgCode() != NEUTRINO_TAU) &&
          (fEvent->GetParticle(i)->GetPdgCode() != ANTINEUTRINO_E) &&
          (fEvent->GetParticle(i)->GetPdgCode() != ANTINEUTRINO_MUON) &&
          (fEvent->GetParticle(i)->GetPdgCode() != ANTINEUTRINO_TAU) ) {
         // Fill histogram for particle's energy loss
         fHisto_dEdX->Fill(fEvent->GetParticle(i)->GetELoss());
         for (j=0;j<fEvent->GetParticle(i)->GetNTracks();j++)
           fEvent->GetParticle(i)->GetTrack(j)->Draw();
      }
   }
   // Reparent each list tree item regarding the
   // corresponding particle relations
   for (i=1;i<=fEvent->GetTotal();i++) {
      fEventListTree->Reparent(gLTI[i],
            gLTI[fEvent->GetParticle(i)->GetFirstMother()]);
   }
   fEventListTree->OpenItem(gBaseLTI);
   fEventListTree->OpenItem(gLTI[0]);
   fEventListTree->ClearViewPort();
   fClient->NeedRedraw(fEventListTree);
   AppendPad();

   sprintf(strtmp,"Done - Total particles : %d - Waiting for next simulation",
                   fEvent->GetTotal());
   fStatusBar->SetText(strtmp,0);
   fPadC->cd();
   // do not fit if not enough particles
   if (fHisto_dEdX->GetEntries() > 10) {
      fHisto_dEdX->Fit("landau","L");
      TF1 *f1 = fHisto_dEdX->GetFunction("landau");
      //delete fit function is fit is a non sense
      if (f1 && f1->GetNDF() > 0) {
         f1->SetLineColor(kRed);
         f1->SetLineWidth(1);
      } else {
         delete f1;
      }
   }
   fHisto_dEdX->Draw();
   fPadC->Modified();
   fPadC->Update();
   fCC->Update();
   fPadC->cd();
   fPadC->SetFillColor(16);
   fPadC->GetFrame()->SetFillColor(10);
   fPadC->Draw();
   fPadC->Update();

   fCA->cd();
   fCA->Modified();
   fCA->Update();
   fMenuEvent->EnableEntry(M_SHOW_INFOS);
   if (gGLManager) {
      fMenuEvent->EnableEntry(M_SHOW_3D);
      fToolBar->GetButton(M_SHOW_3D)->SetState(kButtonUp);
   }
   fToolBar->GetButton(M_FILE_SAVEAS)->SetState(kButtonUp);
   fMenuFile->EnableEntry(M_FILE_SAVEAS);
   fButtonFrame->SetState(GButtonFrame::kAllActive);
   sprintf(strtmp,"Root Shower Event Display - %s",filename);
   SetWindowName(strtmp);
}

////////////////////////////////////////////////////////////////////////////////
/// Saves current event into a Root file

void RootShower::OnSaveFile(const Char_t *filename)
{
   TFile *hfile;
   char  strtmp[256];
   gGeoManager->Export(filename, "detector");
   hfile = new TFile(filename,"UPDATE","Root Shower file");
   hfile->SetCompressionLevel(9);
   TTree *hTree = new TTree("RootShower","Root Shower tree");
   hTree->Branch("Event", "MyEvent", &fEvent, 8000, 2);
   hTree->Fill();  //fill the tree
   hTree->Write();
   hTree->Print();
   hfile->Close();
   sprintf(strtmp,"Root Shower Event Display - %s",filename);
   SetWindowName(strtmp);
}

////////////////////////////////////////////////////////////////////////////////
/// Gives infos on current event

void RootShower::ShowInfos()
{
   Window_t wdummy;
   int ax, ay;
   TRootHelpDialog *hd;
   Char_t str[32];
   Char_t Msg[500];
   Double_t dimx,dimy,dimz;

   fEvent->GetDetector()->GetDimensions(&dimx, &dimy, &dimz);

   sprintf(Msg, "  Some information about the current shower\n");
   sprintf(Msg, "%s  Dimensions of the target\n", Msg);
   sprintf(Msg, "%s  X .................... : %1.2e [cm]    \n", Msg, dimx);
   sprintf(Msg, "%s  Y .................... : %1.2e [cm]    \n", Msg, dimy);
   sprintf(Msg, "%s  Z .................... : %1.2e [cm]    \n", Msg, dimz);
   sprintf(Msg, "%s  Magnetic field ....... : %1.2e [kGauss]\n", Msg,
           fEvent->GetB());
   sprintf(Msg, "%s  Initial particle ..... : %s \n", Msg,
           fEvent->GetParticle(0)->GetName());
   sprintf(Msg, "%s  Initial energy ....... : %1.2e [GeV] \n", Msg,
           fEvent->GetHeader()->GetEnergy());
   sprintf(Msg, "%s  Total Energy loss .... : %1.2e [GeV]", Msg,
           fEvent->GetDetector()->GetTotalELoss());

   sprintf(str, "Infos on current shower");
   hd = new TRootHelpDialog(this, str, 420, 155);
   hd->SetText(Msg);
   gVirtualX->TranslateCoordinates( GetId(), GetParent()->GetId(),
              (Int_t)(GetWidth() - 420) >> 1,(Int_t)(GetHeight() - 155) >> 1,
              ax, ay, wdummy);
   hd->Move(ax, ay);
   hd->Popup();
   fClient->WaitFor(hd);
}

////////////////////////////////////////////////////////////////////////////////
/// Handle keyboard events.

Bool_t RootShower::HandleKey(Event_t *event)
{
   char   input[10];
   UInt_t keysym;

   if (event->fType == kGKeyPress) {
      gVirtualX->LookupString(event, input, sizeof(input), keysym);

      switch ((EKeySym)keysym) {   // ignore these keys
         case kKey_Shift:
         case kKey_Control:
         case kKey_Meta:
         case kKey_Alt:
         case kKey_CapsLock:
         case kKey_NumLock:
         case kKey_ScrollLock:
            return kTRUE;
         case kKey_F1:
            SendMessage(this, MK_MSG(kC_COMMAND, kCM_MENU),
                        M_HELP_SIMULATION, 0);
            return kTRUE;
         default:
            break;
      }
      if (event->fState & kKeyControlMask) {   // Cntrl key modifier pressed
         switch ((EKeySym)keysym & ~0x20) {   // treat upper and lower the same
            case kKey_A:
               SendMessage(this, MK_MSG(kC_COMMAND, kCM_MENU),
                           M_FILE_SAVEAS, 0);
               return kTRUE;
            case kKey_B:
               SendMessage(this, MK_MSG(kC_COMMAND, kCM_MENU),
                           M_INSPECT_BROWSER, 0);
               return kTRUE;
            case kKey_I:
               SendMessage(this, MK_MSG(kC_COMMAND, kCM_MENU),
                           M_SHOW_INFOS, 0);
               return kTRUE;
            case kKey_O:
               SendMessage(this, MK_MSG(kC_COMMAND, kCM_MENU),
                           M_FILE_OPEN, 0);
               return kTRUE;
            case kKey_Q:
               SendMessage(this, MK_MSG(kC_COMMAND, kCM_MENU),
                           M_FILE_EXIT, 0);
               return kTRUE;
            default:
               break;
         }
      }
   }
   return TGMainFrame::HandleKey(event);
}

////////////////////////////////////////////////////////////////////////////////
/// Logo animation timer handling.

Bool_t RootShower::HandleTimer(TTimer *)
{
   if (fPicIndex > fPicNumber) fPicIndex = 1;
   fTitleFrame->ChangeRightLogo(fPicIndex);
   fPicIndex++;
   fTimer->Reset();
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute distance from point px,py to objects in event

Int_t RootShower::DistancetoPrimitive(Int_t px, Int_t py)
{
   Int_t i, j;
   Int_t dist = 9999;

   if (fEvent->GetTotal() <= 0) return 0;
   // Browse every track and get related particle infos.
   for (i=0;i<fEvent->GetTotal();i++) {
      for (j=0;j<fEvent->GetParticle(i)->GetNTracks();j++) {
         dist = fEvent->GetParticle(i)->GetTrack(j)->DistancetoPrimitive(px, py);
         if (dist < 2) {
            gPad->SetSelected((TObject*)fEvent->GetParticle(i));
            fStatusBar->SetText(fEvent->GetParticle(i)->GetObjectInfo(px, py),1);
            gPad->SetCursor(kPointer);
            return 0;
         }
      }
   }
   gPad->SetSelected((TObject*)gPad->GetView());
   return gPad->GetView()->DistancetoPrimitive(px,py);
}

////////////////////////////////////////////////////////////////////////////////
/// Process mouse clicks in TGListTree.

void RootShower::Clicked(TGListTreeItem *item, Int_t x, Int_t y)
{
   MyParticle *part = (MyParticle *) item->GetUserData();
   if (part) {
      fContextMenu->Popup(x, y, part);
   }
   fEventListTree->ClearViewPort();
}

////////////////////////////////////////////////////////////////////////////////
/// Main (entry point).

int main(int argc, char **argv)
{
   Bool_t rint = kFALSE;
   for (int i = 0; i < argc; i++) {
      if (!strcmp(argv[i], "-d")) rint = kTRUE;
      if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "-?")) {
         printf("Usage: %s [-d] [-h | -?]\n", argv[0]);
         printf("  -d:     debug and inspect mode via ROOT prompt\n");
         printf("  -h, -?: this message\n");
         return 0;
      }
   }

   TApplication *theApp;
   if (rint)
      theApp = new TRint("App", &argc, argv);
   else
      theApp = new TApplication("App", &argc, argv);

   gStyle->SetOptStat(1111);
   gStyle->SetOptFit(1111);
   gStyle->SetStatFont(42);

   if (!gGLManager) {
      TString x = "win32";
      if (gVirtualX->InheritsFrom("TGX11"))
         x = "x11";
      TPluginHandler *ph = gROOT->GetPluginManager()->FindHandler("TGLManager", x);
      if (ph && ph->LoadPlugin() != -1)
         ph->ExecPlugin(0);
   }

   gRandom->SetSeed( (UInt_t)time( NULL ) );
   const Int_t NRGBs = 5;
   Double_t Stops[NRGBs] = { 0.00, 0.50, 0.75, 0.875, 1.00 };
   Double_t Red[NRGBs] = { 1.00, 1.00, 1.00, 1.00, 1.00 };
   Double_t Green[NRGBs] = { 1.00, 0.75, 0.50, 0.25, 0.00 };
   Double_t Blue[NRGBs] = { 0.00, 0.00, 0.00, 0.00, 0.00 };
   gColIndex = TColor::CreateGradientColorTable(NRGBs, Stops, Red, Green, Blue, 17);

   // Create RootShower
   RootShower theShower(gClient->GetRoot(), 400, 200);

   // run ROOT application
   theApp->Run();

   // pro forma, never reached
   delete theApp;

   return 0;
}
