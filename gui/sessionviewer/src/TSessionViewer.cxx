// @(#)root/sessionviewer:$Id$
// Author: Marek Biskup, Jakub Madejczyk, Bertrand Bellenot 10/08/2005

/*************************************************************************
 * Copyright (C) 1995-2005, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TSessionViewer                                                       //
//                                                                      //
// Widget used to manage PROOF or local sessions, PROOF connections,    //
// queries construction and results handling.                           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TApplication.h"
#include "TROOT.h"
#include "THashList.h"
#include "TClass.h"
#include "TSystem.h"
#include "TGFileDialog.h"
#include "TBrowser.h"
#include "TGButton.h"
#include "TGLayout.h"
#include "TGListTree.h"
#include "TGCanvas.h"
#include "TGLabel.h"
#include "TGTextEntry.h"
#include "TGNumberEntry.h"
#include "TGTableLayout.h"
#include "TGListBox.h"
#include "TGSplitter.h"
#include "TGProgressBar.h"
#include "TGListView.h"
#include "TGMsgBox.h"
#include "TGMenu.h"
#include "TGStatusBar.h"
#include "TGIcon.h"
#include "TChain.h"
#include "TDSet.h"
#include "TFileInfo.h"
#include "TObjString.h"
#include "TProof.h"
#include "TRandom.h"
#include "TSessionViewer.h"
#include "TSessionLogView.h"
#include "TQueryResult.h"
#include "TGTextView.h"
#include "TGToolBar.h"
#include "TGTab.h"
#include "TRootEmbeddedCanvas.h"
#include "TCanvas.h"
#include "TGMimeTypes.h"
#include "TInterpreter.h"
#include "TContextMenu.h"
#include "TG3DLine.h"
#include "TSessionDialogs.h"
#include "TEnv.h"
#include "TH2.h"
#include "TVirtualTreePlayer.h"
#include "TSelector.h"
#include "TFileCollection.h"
#include "TVirtualX.h"
#include "snprintf.h"
#ifdef WIN32
#include "TWin32SplashThread.h"
#endif

#include <cstdlib>

TSessionViewer *gSessionViewer = nullptr;

const char *kConfigFile = ".proofgui.conf";

ClassImp(TQueryDescription);
ClassImp(TSessionDescription);
ClassImp(TSessionServerFrame);
ClassImp(TSessionFrame);
ClassImp(TSessionQueryFrame);
ClassImp(TSessionOutputFrame);
ClassImp(TSessionInputFrame);
ClassImp(TSessionViewer);

const char *xpm_names[] = {
    "monitor01.xpm",
    "monitor02.xpm",
    "monitor03.xpm",
    "monitor04.xpm",
    0
};

const char *conftypes[] = {
   "Config files",  "*.conf",
   "All files",     "*.*",
    0,               0
};

const char *pkgtypes[] = {
   "Package files", "*.par",
   "All files",     "*.*",
    0,               0
};

const char *macrotypes[] = {
   "C files",       "*.[C|c]*",
   "All files",     "*",
   0,               0
};


const char *kFeedbackHistos[] = {
   "PROOF_PacketsHist",
   "PROOF_EventsHist",
   "PROOF_NodeHist",
   "PROOF_LatencyHist",
   "PROOF_ProcTimeHist",
   "PROOF_CpuTimeHist",
   0
};

const char* const kSession_RedirectFile = ".templog";
const char* const kSession_RedirectCmd = ".tempcmd";

// Menu command id's
enum ESessionViewerCommands {
   kFileLoadConfig,
   kFileSaveConfig,
   kFileCloseViewer,
   kFileQuit,

   kSessionNew,
   kSessionAdd,
   kSessionDelete,
   kSessionGetQueries,

   kSessionConnect,
   kSessionDisconnect,
   kSessionShutdown,
   kSessionCleanup,
   kSessionBrowse,
   kSessionShowStatus,
   kSessionReset,

   kQueryNew,
   kQueryEdit,
   kQueryDelete,
   kQuerySubmit,
   kQueryStartViewer,

   kOptionsAutoSave,
   kOptionsStatsHist,
   kOptionsStatsTrace,
   kOptionsSlaveStatsTrace,
   kOptionsFeedback,

   kHelpAbout
};

const char *xpm_toolbar[] = {
    "fileopen.xpm",
    "filesaveas.xpm",
    "",
    "connect.xpm",
    "disconnect.xpm",
    "",
    "query_new.xpm",
    "query_submit.xpm",
    "",
    "about.xpm",
    "",
    "quit.xpm",
    0
};

ToolBarData_t tb_data[] = {
  { "", "Open Config File",     kFALSE, kFileLoadConfig,    0 },
  { "", "Save Config File",     kFALSE, kFileSaveConfig,    0 },
  { "", 0,                      0,      -1,                 0 },
  { "", "Connect",              kFALSE, kSessionConnect,    0 },
  { "", "Disconnect",           kFALSE, kSessionDisconnect, 0 },
  { "", 0,                      0,      -1,                 0 },
  { "", "New Query",            kFALSE, kQueryNew,          0 },
  { "", "Submit Query",         kFALSE, kQuerySubmit,       0 },
  { "", 0,                      0,      -1,                 0 },
  { "", "About Root",           kFALSE, kHelpAbout,         0 },
  { "", 0,                      0,      -1,                 0 },
  { "", "Exit Root",            kFALSE, kFileQuit,          0 },
  { 0,  0,                      0,      0,                  0 }
};


////////////////////////////////////////////////////////////////////////////////
// Server Frame

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TSessionServerFrame::TSessionServerFrame(TGWindow* p, Int_t w, Int_t h) :
   TGCompositeFrame(p, w, h), fFrmNewServer(0), fTxtName(0), fTxtAddress(0),
   fNumPort(0), fLogLevel(0), fTxtConfig(0), fTxtUsrName(0), fSync(0),
   fViewer(0), fBtnAdd(0), fBtnConnect(0)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TSessionServerFrame::~TSessionServerFrame()
{
   Cleanup();
}

////////////////////////////////////////////////////////////////////////////////
/// Build server configuration frame.

void TSessionServerFrame::Build(TSessionViewer *gui)
{
   SetLayoutManager(new TGVerticalLayout(this));

   SetCleanup(kDeepCleanup);

   fViewer = gui;
   fFrmNewServer = new TGGroupFrame(this, "New Session");
   fFrmNewServer->SetCleanup(kDeepCleanup);

   AddFrame(fFrmNewServer, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

   fFrmNewServer->SetLayoutManager(new TGMatrixLayout(fFrmNewServer, 0, 2, 8));

   fFrmNewServer->AddFrame(new TGLabel(fFrmNewServer, "Session Name:"),
                           new TGLayoutHints(kLHintsLeft, 3, 3, 3, 3));
   fFrmNewServer->AddFrame(fTxtName = new TGTextEntry(fFrmNewServer,
                           (const char *)0, 1), new TGLayoutHints());
   fTxtName->Resize(156, fTxtName->GetDefaultHeight());
   fTxtName->Associate(this);
   fFrmNewServer->AddFrame(new TGLabel(fFrmNewServer, "Server name:"),
                           new TGLayoutHints(kLHintsLeft, 3, 3, 3, 3));
   fFrmNewServer->AddFrame(fTxtAddress = new TGTextEntry(fFrmNewServer,
                           (const char *)0, 2), new TGLayoutHints());
   fTxtAddress->Resize(156, fTxtAddress->GetDefaultHeight());
   fTxtAddress->Associate(this);
   fFrmNewServer->AddFrame(new TGLabel(fFrmNewServer, "Port (default: 1093):"),
                           new TGLayoutHints(kLHintsLeft, 3, 3, 3, 3));
   fFrmNewServer->AddFrame(fNumPort = new TGNumberEntry(fFrmNewServer, 1093, 5,
            3, TGNumberFormat::kNESInteger,TGNumberFormat::kNEANonNegative,
            TGNumberFormat::kNELLimitMinMax, 0, 65535),new TGLayoutHints());
   fNumPort->Associate(this);
   fFrmNewServer->AddFrame(new TGLabel(fFrmNewServer, "Configuration File:"),
                           new TGLayoutHints(kLHintsLeft, 3, 3, 3, 3));
   fFrmNewServer->AddFrame(fTxtConfig = new TGTextEntry(fFrmNewServer,
                           (const char *)0, 4), new TGLayoutHints());
   fTxtConfig->Resize(156, fTxtConfig->GetDefaultHeight());
   fTxtConfig->Associate(this);
   fFrmNewServer->AddFrame(new TGLabel(fFrmNewServer, "Log Level:"),
                           new TGLayoutHints(kLHintsLeft, 3, 3, 3, 3));

   fFrmNewServer->AddFrame(fLogLevel = new TGNumberEntry(fFrmNewServer, 0, 5, 5,
                           TGNumberFormat::kNESInteger,
                           TGNumberFormat::kNEANonNegative,
                           TGNumberFormat::kNELLimitMinMax, 0, 5),
                           new TGLayoutHints(kLHintsLeft, 3, 3, 3, 3));
   fLogLevel->Associate(this);

   fFrmNewServer->AddFrame(new TGLabel(fFrmNewServer, "User Name:"),
            new TGLayoutHints(kLHintsLeft, 3, 3, 3, 3));
   fFrmNewServer->AddFrame(fTxtUsrName = new TGTextEntry(fFrmNewServer,
                           (const char *)0, 6), new TGLayoutHints());
   fTxtUsrName->Resize(156, fTxtUsrName->GetDefaultHeight());
   fTxtUsrName->Associate(this);

   fFrmNewServer->AddFrame(new TGLabel(fFrmNewServer, "Process mode :"),
            new TGLayoutHints(kLHintsLeft | kLHintsBottom | kLHintsExpandX,
            3, 3, 3, 3));
   fFrmNewServer->AddFrame(fSync = new TGCheckButton(fFrmNewServer,
      "&Synchronous"), new TGLayoutHints(kLHintsLeft | kLHintsBottom |
      kLHintsExpandX, 3, 3, 3, 3));
   fSync->SetToolTipText("Default Process Mode");
   fSync->SetState(kButtonDown);

   AddFrame(fBtnAdd = new TGTextButton(this, "             Save             "),
            new TGLayoutHints(kLHintsTop | kLHintsCenterX, 5, 5, 15, 5));
   fBtnAdd->SetToolTipText("Add server to the list");
   fBtnAdd->Connect("Clicked()", "TSessionServerFrame", this,
                    "OnBtnAddClicked()");
   AddFrame(fBtnConnect = new TGTextButton(this, "          Connect          "),
                 new TGLayoutHints(kLHintsTop | kLHintsCenterX, 5, 5, 15, 5));
   fBtnConnect->Connect("Clicked()", "TSessionServerFrame", this,
                        "OnBtnConnectClicked()");
   fBtnConnect->SetToolTipText("Connect to the selected server");

   fTxtConfig->Connect("DoubleClicked()", "TSessionServerFrame", this,
                       "OnConfigFileClicked()");

   fTxtName->Connect("TextChanged(char*)", "TSessionServerFrame", this,
                     "SettingsChanged()");
   fTxtAddress->Connect("TextChanged(char*)", "TSessionServerFrame", this,
                       "SettingsChanged()");
   fTxtConfig->Connect("TextChanged(char*)", "TSessionServerFrame", this,
                       "SettingsChanged()");
   fTxtUsrName->Connect("TextChanged(char*)", "TSessionServerFrame", this,
                        "SettingsChanged()");
   fSync->Connect("Clicked()", "TSessionServerFrame", this,
                  "SettingsChanged()");
   fLogLevel->Connect("ValueChanged(Long_t)", "TSessionServerFrame", this,
                      "SettingsChanged()");
   fLogLevel->Connect("ValueSet(Long_t)", "TSessionServerFrame", this,
                      "SettingsChanged()");
   fNumPort->Connect("ValueChanged(Long_t)", "TSessionServerFrame", this,
                     "SettingsChanged()");
   fNumPort->Connect("ValueSet(Long_t)", "TSessionServerFrame", this,
                     "SettingsChanged()");

}

////////////////////////////////////////////////////////////////////////////////
/// Settings have changed, update GUI accordingly.

void TSessionServerFrame::SettingsChanged()
{
   TGTextEntry *sender = dynamic_cast<TGTextEntry*>((TQObject*)gTQSender);
   Bool_t issync = (fSync->GetState() == kButtonDown);
   if ((fViewer->GetActDesc()->fLocal) ||
       (strcmp(fViewer->GetActDesc()->GetName(), fTxtName->GetText())) ||
       (strcmp(fViewer->GetActDesc()->fAddress.Data(), fTxtAddress->GetText())) ||
       (strcmp(fViewer->GetActDesc()->fConfigFile.Data(), fTxtConfig->GetText())) ||
       (strcmp(fViewer->GetActDesc()->fUserName.Data(), fTxtUsrName->GetText())) ||
       (fViewer->GetActDesc()->fLogLevel != fLogLevel->GetIntNumber()) ||
       (fViewer->GetActDesc()->fPort != fNumPort->GetIntNumber()) ||
       (fViewer->GetActDesc()->fSync != issync)) {
      ShowFrame(fBtnAdd);
      HideFrame(fBtnConnect);
   }
   else {
      HideFrame(fBtnAdd);
      ShowFrame(fBtnConnect);
   }
   if (sender) {
      sender->SetFocus();
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Handle expose event in server frame.

Bool_t TSessionServerFrame::HandleExpose(Event_t * /*event*/)
{
   //fTxtName->SelectAll();
   //fTxtName->SetFocus();
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Browse configuration files.

void TSessionServerFrame::OnConfigFileClicked()
{
   // do nothing if connection in progress
   if (fViewer->IsBusy())
      return;
   TGFileInfo fi;
   fi.fFileTypes = conftypes;
   new TGFileDialog(fClient->GetRoot(), fViewer, kFDOpen, &fi);
   if (!fi.fFilename) return;
   fTxtConfig->SetText(gSystem->BaseName(fi.fFilename));
}

////////////////////////////////////////////////////////////////////////////////
/// Delete selected session configuration (remove it from the list).

void TSessionServerFrame::OnBtnDeleteClicked()
{
   // do nothing if connection in progress
   if (fViewer->IsBusy())
      return;
   TString name(fTxtName->GetText());
   TIter next(fViewer->GetSessions());
   TSessionDescription *desc = fViewer->GetActDesc();

   if (desc->fLocal) {
      Int_t retval;
      new TGMsgBox(fClient->GetRoot(), this, "Error Deleting Session",
                   "Deleting Local Sessions is not allowed !",
                    kMBIconExclamation,kMBOk,&retval);
      return;
   }
   // ask for confirmation
   TString m;
   m.Form("Are you sure to delete the server \"%s\"",
          desc->fName.Data());
   Int_t result;
   new TGMsgBox(fClient->GetRoot(), this, "", m.Data(), 0,
                kMBOk | kMBCancel, &result);
   // if confirmed, delete it
   if (result == kMBOk) {
      // remove the Proof session from gROOT list of Proofs
      if (desc->fConnected && desc->fAttached && desc->fProof) {
         desc->fProof->Detach("S");
      }
      // remove it from our sessions list
      fViewer->GetSessions()->Remove((TObject *)desc);
      // update configuration file
      TGListTreeItem *item = fViewer->GetSessionHierarchy()->GetSelected();
      fViewer->GetSessionHierarchy()->DeleteItem(item);

      TObject *obj = fViewer->GetSessions()->Last();
      item = fViewer->GetSessionHierarchy()->FindChildByData(
               fViewer->GetSessionItem(), (void *)obj);
      if (item) {
         fViewer->GetSessionHierarchy()->ClearHighlighted();
         fViewer->GetSessionHierarchy()->OpenItem(item);
         fViewer->GetSessionHierarchy()->HighlightItem(item);
         fViewer->GetSessionHierarchy()->SetSelected(item);
         fViewer->GetSessionHierarchy()->ClearViewPort();
         fClient->NeedRedraw(fViewer->GetSessionHierarchy());
         fViewer->OnListTreeClicked(item, 1, 0, 0);
      }
   }
   if (fViewer->IsAutoSave())
      fViewer->WriteConfiguration();
}

////////////////////////////////////////////////////////////////////////////////
/// Connect to selected server.

void TSessionServerFrame::OnBtnConnectClicked()
{
   // do nothing if connection in progress
   if (fViewer->IsBusy())
      return;

   if (!fViewer->GetSessions()->FindObject(fTxtName->GetText())) {
      OnBtnAddClicked();
   }
   else {
      fViewer->GetActDesc()->fAddress = fTxtAddress->GetText();
      fViewer->GetActDesc()->fPort = fNumPort->GetIntNumber();
      if (strlen(fTxtConfig->GetText()) > 1)
         fViewer->GetActDesc()->fConfigFile = TString(fTxtConfig->GetText());
      else
         fViewer->GetActDesc()->fConfigFile = "";
      fViewer->GetActDesc()->fLogLevel = fLogLevel->GetIntNumber();
      fViewer->GetActDesc()->fUserName = fTxtUsrName->GetText();
      fViewer->GetActDesc()->fSync = (fSync->GetState() == kButtonDown);
      if (fViewer->IsAutoSave())
         fViewer->WriteConfiguration();
   }
   // set flag busy
   fViewer->SetBusy();
   // avoid input events in list tree while connecting
   fViewer->GetSessionHierarchy()->RemoveInput(kPointerMotionMask |
         kEnterWindowMask | kLeaveWindowMask | kKeyPressMask);
   gVirtualX->GrabButton(fViewer->GetSessionHierarchy()->GetId(), kAnyButton,
         kAnyModifier, kButtonPressMask | kButtonReleaseMask, kNone, kNone, kFALSE);
   // set watch cursor to indicate connection in progress
   gVirtualX->SetCursor(fViewer->GetSessionHierarchy()->GetId(),
         gVirtualX->CreateCursor(kWatch));
   gVirtualX->SetCursor(GetId(),gVirtualX->CreateCursor(kWatch));
   // display connection progress bar in first part of status bar
   fViewer->GetStatusBar()->GetBarPart(0)->ShowFrame(fViewer->GetConnectProg());
   // connect to proof startup message (to update progress bar)
   TQObject::Connect("TProof", "StartupMessage(char *,Bool_t,Int_t,Int_t)",
         "TSessionViewer", fViewer, "StartupMessage(char *,Bool_t,Int_t,Int_t)");
   // collect and set-up configuration
   TString address = fTxtAddress->GetText();
   TString url;
   if (address == "lite://") {
      url = address;
   }
   else {
      // collect and set-up configuration
      url = fTxtUsrName->GetText();
      url += "@"; url += address.Data();
      if (fNumPort->GetIntNumber() > 0) {
         url += ":";
         url += fNumPort->GetIntNumber();
      }
   }

   TProofDesc *desc;
   fViewer->GetActDesc()->fProofMgr = TProofMgr::Create(url);
   if (!fViewer->GetActDesc()->fProofMgr ||
       !fViewer->GetActDesc()->fProofMgr->IsValid()) {
      // hide connection progress bar from status bar
      fViewer->GetStatusBar()->GetBarPart(0)->HideFrame(fViewer->GetConnectProg());
      // release busy flag
      fViewer->SetBusy(kFALSE);
      // restore cursors and input
      gVirtualX->SetCursor(GetId(), 0);
      gVirtualX->GrabButton(fViewer->GetSessionHierarchy()->GetId(), kAnyButton,
            kAnyModifier, kButtonPressMask | kButtonReleaseMask, kNone, kNone);
      fViewer->GetSessionHierarchy()->AddInput(kPointerMotionMask |
            kEnterWindowMask | kLeaveWindowMask | kKeyPressMask);
      gVirtualX->SetCursor(fViewer->GetSessionHierarchy()->GetId(), 0);
      return;
   }
   fViewer->UpdateListOfSessions();
   // check if the session already exist before to recreate it
   TList *sessions = fViewer->GetActDesc()->fProofMgr->QuerySessions("");
   if (sessions) {
      TIter nextp(sessions);
      // loop over existing Proof sessions
      while ((desc = (TProofDesc *)nextp())) {
         if ((desc->GetName() == fViewer->GetActDesc()->fTag) ||
             (desc->GetTitle() == fViewer->GetActDesc()->fName)) {
            fViewer->GetActDesc()->fProof =
               fViewer->GetActDesc()->fProofMgr->AttachSession(desc->GetLocalId(), kTRUE);
            fViewer->GetActDesc()->fTag = desc->GetName();
            fViewer->GetActDesc()->fProof->SetAlias(fViewer->GetActDesc()->fName);
            fViewer->GetActDesc()->fConnected = kTRUE;
            fViewer->GetActDesc()->fAttached = kTRUE;

            if (fViewer->GetOptionsMenu()->IsEntryChecked(kOptionsFeedback)) {
               Int_t i = 0;
               // browse list of feedback histos and check user's selected ones
               while (kFeedbackHistos[i]) {
                  if (fViewer->GetCascadeMenu()->IsEntryChecked(41+i)) {
                     fViewer->GetActDesc()->fProof->AddFeedback(kFeedbackHistos[i]);
                     fViewer->GetActDesc()->fNbHistos++;
                  }
                  i++;
               }
               // connect feedback signal
               fViewer->GetActDesc()->fProof->Connect("Feedback(TList *objs)",
                           "TSessionQueryFrame", fViewer->GetQueryFrame(),
                           "Feedback(TList *objs)");
               gROOT->Time();
            }
            else {
               // if feedback option not selected, clear Proof's feedback option
               fViewer->GetActDesc()->fProof->ClearFeedback();
            }

            break;
         }
      }
   }
   if (fViewer->GetActDesc()->fProof == 0) {
      if (fViewer->GetActDesc()->fProofMgr->IsValid()) {
         fViewer->GetActDesc()->fProof = fViewer->GetActDesc()->fProofMgr->CreateSession(
         fViewer->GetActDesc()->fConfigFile);
         desc = 0;
         sessions = fViewer->GetActDesc()->fProofMgr->QuerySessions("");
         if (sessions)
            desc = (TProofDesc *)sessions->Last();
         if (desc) {
            fViewer->GetActDesc()->fProof->SetAlias(fViewer->GetActDesc()->fName);
            fViewer->GetActDesc()->fTag = desc->GetName();
            fViewer->GetActDesc()->fConnected = kTRUE;
            fViewer->GetActDesc()->fAttached = kTRUE;
         }
      }
   }
   if (fViewer->GetActDesc()->fProof) {
      fViewer->GetActDesc()->fConfigFile = fViewer->GetActDesc()->fProof->GetConfFile();
      fViewer->GetActDesc()->fUserName   = fViewer->GetActDesc()->fProof->GetUser();
      fViewer->GetActDesc()->fPort       = fViewer->GetActDesc()->fProof->GetPort();
      fViewer->GetActDesc()->fLogLevel   = fViewer->GetActDesc()->fProof->GetLogLevel();
      if (fViewer->GetActDesc()->fLogLevel < 0)
         fViewer->GetActDesc()->fLogLevel = 0;
      if (fViewer->GetActDesc()->fProof->IsLite())
         fViewer->GetActDesc()->fAddress = "lite://";
      else
         fViewer->GetActDesc()->fAddress = fViewer->GetActDesc()->fProof->GetMaster();
      fViewer->GetActDesc()->fConnected = kTRUE;
      fViewer->GetActDesc()->fProof->SetBit(TProof::kUsingSessionGui);
   }
   fViewer->UpdateListOfSessions();

   // check if connected and valid
   if (fViewer->GetActDesc()->fProof &&
       fViewer->GetActDesc()->fProof->IsValid()) {
      // set log level
      fViewer->GetActDesc()->fProof->SetLogLevel(fViewer->GetActDesc()->fLogLevel);
      // set query type (synch / asynch)
      fViewer->GetActDesc()->fProof->SetQueryMode(fViewer->GetActDesc()->fSync ?
                             TProof::kSync : TProof::kAsync);
      // set connected flag
      fViewer->GetActDesc()->fConnected = kTRUE;
      // change list tree item picture to connected pixmap
      TGListTreeItem *item = fViewer->GetSessionHierarchy()->FindChildByData(
                              fViewer->GetSessionItem(),fViewer->GetActDesc());
      if (item) {
         item->SetPictures(fViewer->GetProofConPict(), fViewer->GetProofConPict());
         // update viewer
         fViewer->OnListTreeClicked(item, 1, 0, 0);
         fViewer->GetSessionHierarchy()->ClearViewPort();
         fClient->NeedRedraw(fViewer->GetSessionHierarchy());
      }
      // connect to progress related signals
      fViewer->GetActDesc()->fProof->Connect("Progress(Long64_t,Long64_t)",
                 "TSessionQueryFrame", fViewer->GetQueryFrame(),
                 "Progress(Long64_t,Long64_t)");
      fViewer->GetActDesc()->fProof->Connect("Progress(Long64_t,Long64_t,Long64_t,Float_t,Float_t,Float_t,Float_t)",
                 "TSessionQueryFrame", fViewer->GetQueryFrame(),
                 "Progress(Long64_t,Long64_t,Long64_t,Float_t,Float_t,Float_t,Float_t)");
      fViewer->GetActDesc()->fProof->Connect("Progress(Long64_t,Long64_t,Long64_t,Float_t,Float_t,Float_t,Float_t,Int_t,Int_t,Float_t)",
                 "TSessionQueryFrame", fViewer->GetQueryFrame(),
                 "Progress(Long64_t,Long64_t,Long64_t,Float_t,Float_t,Float_t,Float_t,Int_t,Int_t,Float_t)");
      fViewer->GetActDesc()->fProof->Connect("StopProcess(Bool_t)",
                 "TSessionQueryFrame", fViewer->GetQueryFrame(),
                 "IndicateStop(Bool_t)");
      fViewer->GetActDesc()->fProof->Connect(
                 "ResetProgressDialog(const char*,Int_t,Long64_t,Long64_t)",
                 "TSessionQueryFrame", fViewer->GetQueryFrame(),
                 "ResetProgressDialog(const char*,Int_t,Long64_t,Long64_t)");
      // enable timer used for status bar icon's animation
      fViewer->EnableTimer();
      // change status bar right icon to connected pixmap
      fViewer->ChangeRightLogo("monitor01.xpm");
      // do not animate yet
      fViewer->SetChangePic(kFALSE);
      // connect to signal "query result ready"
      fViewer->GetActDesc()->fProof->Connect("QueryResultReady(char *)",
                       "TSessionViewer", fViewer, "QueryResultReady(char *)");
      // display connection information on status bar
      TString msg;
      msg.Form("PROOF Cluster %s ready", fViewer->GetActDesc()->fName.Data());
      fViewer->GetStatusBar()->SetText(msg.Data(), 1);
      fViewer->GetSessionFrame()->ProofInfos();
      fViewer->UpdateListOfPackages();
      fViewer->GetSessionFrame()->UpdateListOfDataSets();
      // Enable previously uploaded packages if in auto-enable mode
      if (fViewer->GetActDesc()->fAutoEnable) {
         TPackageDescription *package;
         TIter next(fViewer->GetActDesc()->fPackages);
         while ((package = (TPackageDescription *)next())) {
            if (!package->fEnabled) {
               if (fViewer->GetActDesc()->fProof->EnablePackage(package->fName) != 0)
                  Error("Submit", "Enable package failed");
               else {
                  package->fEnabled = kTRUE;
                  fViewer->GetSessionFrame()->UpdatePackages();
               }
            }
         }
      }
   }
   // hide connection progress bar from status bar
   fViewer->GetStatusBar()->GetBarPart(0)->HideFrame(fViewer->GetConnectProg());
   // release busy flag
   fViewer->SetBusy(kFALSE);
   // restore cursors and input
   gVirtualX->SetCursor(GetId(), 0);
   gVirtualX->GrabButton(fViewer->GetSessionHierarchy()->GetId(), kAnyButton,
         kAnyModifier, kButtonPressMask | kButtonReleaseMask, kNone, kNone);
   fViewer->GetSessionHierarchy()->AddInput(kPointerMotionMask |
         kEnterWindowMask | kLeaveWindowMask | kKeyPressMask);
   gVirtualX->SetCursor(fViewer->GetSessionHierarchy()->GetId(), 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Reset server configuration fields.

void TSessionServerFrame::OnBtnNewServerClicked()
{
   // do nothing if connection in progress
   if (fViewer->IsBusy())
      return;
   fViewer->GetSessionHierarchy()->ClearHighlighted();
   fViewer->GetSessionHierarchy()->OpenItem(fViewer->GetSessionItem());
   fViewer->GetSessionHierarchy()->HighlightItem(fViewer->GetSessionItem());
   fViewer->GetSessionHierarchy()->SetSelected(fViewer->GetSessionItem());
   fViewer->OnListTreeClicked(fViewer->GetSessionItem(), 1, 0, 0);
   fViewer->GetSessionHierarchy()->ClearViewPort();
   fClient->NeedRedraw(fViewer->GetSessionHierarchy());
   fTxtName->SetText("");
   fTxtAddress->SetText("");
   fTxtConfig->SetText("");
   fNumPort->SetIntNumber(1093);
   fLogLevel->SetIntNumber(0);
   fTxtUsrName->SetText("");
}

////////////////////////////////////////////////////////////////////////////////
/// Add newly created session configuration in the list of sessions.

void TSessionServerFrame::OnBtnAddClicked()
{
   Int_t retval;
   Bool_t newSession = kTRUE;
   TSessionDescription* desc = 0;
   // do nothing if connection in progress
   if (fViewer->IsBusy())
      return;

   if ((!fTxtName->GetBuffer()->GetTextLength()) ||
       (!fTxtAddress->GetBuffer()->GetTextLength()) ||
       (!fTxtUsrName->GetBuffer()->GetTextLength())) {
      new TGMsgBox(fClient->GetRoot(), fViewer, "Error Adding Session",
                   "At least one required field is empty !",
                    kMBIconExclamation, kMBOk, &retval);
      return;
   }
   TObject *obj = fViewer->GetSessions()->FindObject(fTxtName->GetText());
   if (obj)
      desc = dynamic_cast<TSessionDescription*>(obj);
   if (desc) {
      new TGMsgBox(fClient->GetRoot(), fViewer, "Adding Session",
                   TString::Format("The session \"%s\" already exists ! Overwrite ?",
                   fTxtName->GetText()), kMBIconQuestion, kMBYes | kMBNo |
                   kMBCancel, &retval);
      if (retval != kMBYes)
         return;
      newSession = kFALSE;
   }
   if (newSession) {
      desc = new TSessionDescription();
      desc->fName = fTxtName->GetText();
      desc->fTag = "";
      desc->fQueries = new TList();
      desc->fPackages = new TList();
      desc->fActQuery = 0;
      desc->fProof = 0;
      desc->fProofMgr = 0;
      desc->fAutoEnable = kFALSE;
      desc->fAddress = fTxtAddress->GetText();
      desc->fPort = fNumPort->GetIntNumber();
      desc->fConnected = kFALSE;
      desc->fAttached = kFALSE;
      desc->fLocal = kFALSE;
      if (strlen(fTxtConfig->GetText()) > 1)
         desc->fConfigFile = TString(fTxtConfig->GetText());
      else
         desc->fConfigFile = "";
      desc->fLogLevel = fLogLevel->GetIntNumber();
      desc->fUserName = fTxtUsrName->GetText();
      desc->fSync = (fSync->GetState() == kButtonDown);
      // add newly created session config to our session list
      fViewer->GetSessions()->Add((TObject *)desc);
      // save into configuration file
      TGListTreeItem *item = fViewer->GetSessionHierarchy()->AddItem(
            fViewer->GetSessionItem(), desc->fName.Data(),
            fViewer->GetProofDisconPict(), fViewer->GetProofDisconPict());
      fViewer->GetSessionHierarchy()->SetToolTipItem(item, "Proof Session");
      item->SetUserData(desc);
      fViewer->GetSessionHierarchy()->ClearHighlighted();
      fViewer->GetSessionHierarchy()->OpenItem(fViewer->GetSessionItem());
      fViewer->GetSessionHierarchy()->OpenItem(item);
      fViewer->GetSessionHierarchy()->HighlightItem(item);
      fViewer->GetSessionHierarchy()->SetSelected(item);
      fViewer->GetSessionHierarchy()->ClearViewPort();
      fClient->NeedRedraw(fViewer->GetSessionHierarchy());
      fViewer->OnListTreeClicked(item, 1, 0, 0);
   }
   else {
      fViewer->GetActDesc()->fName = fTxtName->GetText();
      fViewer->GetActDesc()->fAddress = fTxtAddress->GetText();
      fViewer->GetActDesc()->fPort = fNumPort->GetIntNumber();
      if (strlen(fTxtConfig->GetText()) > 1)
         fViewer->GetActDesc()->fConfigFile = TString(fTxtConfig->GetText());
      fViewer->GetActDesc()->fLogLevel = fLogLevel->GetIntNumber();
      fViewer->GetActDesc()->fUserName = fTxtUsrName->GetText();
      fViewer->GetActDesc()->fSync = (fSync->GetState() == kButtonDown);
      TGListTreeItem *item2 = fViewer->GetSessionHierarchy()->GetSelected();
      item2->SetUserData(fViewer->GetActDesc());
      fViewer->OnListTreeClicked(fViewer->GetSessionHierarchy()->GetSelected(),
                                 1, 0, 0);
   }
   HideFrame(fBtnAdd);
   ShowFrame(fBtnConnect);
   if (fViewer->IsAutoSave())
      fViewer->WriteConfiguration();
}

////////////////////////////////////////////////////////////////////////////////
/// Update fields with values from session description desc.

void TSessionServerFrame::Update(TSessionDescription* desc)
{
   if (desc->fLocal) {
      fTxtName->SetText("");
      fTxtAddress->SetText("");
      fNumPort->SetIntNumber(1093);
      fTxtConfig->SetText("");
      fTxtUsrName->SetText("");
      fLogLevel->SetIntNumber(0);
      return;
   }

   fTxtName->SetText(desc->fName);
   fTxtAddress->SetText(desc->fAddress);
   fNumPort->SetIntNumber(desc->fPort);
   fLogLevel->SetIntNumber(desc->fLogLevel);

   if (desc->fConfigFile.Length() > 1) {
      fTxtConfig->SetText(desc->fConfigFile);
   }
   else {
      fTxtConfig->SetText("");
   }
   fTxtUsrName->SetText(desc->fUserName);
}

////////////////////////////////////////////////////////////////////////////////
/// Process messages for session server frame.
/// Used to navigate between text entry fields.

Bool_t TSessionServerFrame::ProcessMessage(Long_t msg, Long_t parm1, Long_t)
{
   switch (GET_MSG(msg)) {
      case kC_TEXTENTRY:
         switch (GET_SUBMSG(msg)) {
            case kTE_ENTER:
            case kTE_TAB:
               switch (parm1) {
                  case 1: // session name
                     fTxtAddress->SelectAll();
                     fTxtAddress->SetFocus();
                     break;
                  case 2: // server address
                     fNumPort->GetNumberEntry()->SelectAll();
                     fNumPort->GetNumberEntry()->SetFocus();
                     break;
                  case 3: // port number
                     fTxtConfig->SelectAll();
                     fTxtConfig->SetFocus();
                     break;
                  case 4: // configuration file
                     fLogLevel->GetNumberEntry()->SelectAll();
                     fLogLevel->GetNumberEntry()->SetFocus();
                     break;
                  case 5: // log level
                     fTxtUsrName->SelectAll();
                     fTxtUsrName->SetFocus();
                     break;
                  case 6: // user name
                     fTxtName->SelectAll();
                     fTxtName->SetFocus();
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

////////////////////////////////////////////////////////////////////////////////
// Session Frame

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TSessionFrame::TSessionFrame(TGWindow* p, Int_t w, Int_t h) :
   TGCompositeFrame(p, w, h), fTab(0), fFA(0), fFB(0), fFC(0), fFD(0), fFE(0),
   fCommandTxt(0), fCommandBuf(0), fInfoTextView(0), fClearCheck(0),
   fBtnShowLog(0), fBtnNewQuery(0), fBtnGetQueries(0), fLBPackages(0),
   fBtnAdd(0), fBtnRemove(0), fBtnUp(0), fBtnDown(0), fBtnShow(0),
   fBtnShowEnabled(0), fChkMulti(0), fChkEnable(0), fBtnUpload(0),
   fBtnEnable(0), fBtnClear(0), fBtnDisable(0), fDSetView(0), fDataSetTree(0),
   fBtnUploadDSet(0), fBtnRemoveDSet(0), fBtnVerifyDSet(0), fBtnRefresh(0),
   fTxtParallel(0), fLogLevel(0), fApplyLogLevel(0), fApplyParallel(0),
   fViewer(0)
{
   for (int i=0;i<19;++i) fInfoLine[i] = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TSessionFrame::~TSessionFrame()
{
   Cleanup();
}

////////////////////////////////////////////////////////////////////////////////
/// Build session frame.

void TSessionFrame::Build(TSessionViewer *gui)
{
   SetLayoutManager(new TGVerticalLayout(this));
   SetCleanup(kDeepCleanup);
   fViewer  = gui;
   Int_t i,j;

   // main session tab
   fTab = new TGTab(this, 200, 200);
   AddFrame(fTab, new TGLayoutHints(kLHintsTop | kLHintsExpandX |
         kLHintsExpandY, 2, 2, 2, 2));

   // add "Status" tab element
   TGCompositeFrame *tf = fTab->AddTab("Status");
   fFA = new TGCompositeFrame(tf, 100, 100, kVerticalFrame);
   tf->AddFrame(fFA, new TGLayoutHints(kLHintsTop | kLHintsLeft |
         kLHintsExpandX | kLHintsExpandY));

   // add first session information line
   fInfoLine[0] = new TGLabel(fFA, " ");
   fFA->AddFrame(fInfoLine[0], new TGLayoutHints(kLHintsCenterX |
         kLHintsExpandX, 5, 5, 15, 5));

   TGCompositeFrame* frmInfos = new TGHorizontalFrame(fFA, 350, 100);
   frmInfos->SetLayoutManager(new TGTableLayout(frmInfos, 9, 2));

   // add session information lines
   j = 0;
   for (i=0;i<17;i+=2) {
      fInfoLine[i+1] = new TGLabel(frmInfos, " ");
      frmInfos->AddFrame(fInfoLine[i+1], new TGTableLayoutHints(0, 1, j, j+1,
            kLHintsLeft | kLHintsCenterY, 5, 5, 2, 2));
      fInfoLine[i+2] = new TGLabel(frmInfos, " ");
      frmInfos->AddFrame(fInfoLine[i+2], new TGTableLayoutHints(1, 2, j, j+1,
            kLHintsLeft | kLHintsCenterY, 5, 5, 2, 2));
      j++;
   }
   fFA->AddFrame(frmInfos, new TGLayoutHints(kLHintsLeft | kLHintsTop |
         kLHintsExpandX  | kLHintsExpandY, 5, 5, 5, 5));

   // add "new query" and "get queries" buttons
   TGCompositeFrame* frmBut1 = new TGHorizontalFrame(fFA, 350, 100);
   frmBut1->SetCleanup(kDeepCleanup);
   frmBut1->AddFrame(fBtnNewQuery = new TGTextButton(frmBut1, "New Query..."),
         new TGLayoutHints(kLHintsLeft | kLHintsExpandX, 5, 5, 5, 5));
   fBtnNewQuery->SetToolTipText("Open New Query Dialog");
   frmBut1->AddFrame(fBtnGetQueries = new TGTextButton(frmBut1, " Get Queries "),
         new TGLayoutHints(kLHintsLeft | kLHintsExpandX, 5, 5, 5, 5));
   fBtnGetQueries->SetToolTipText("Get List of Queries from the server");
   fBtnShowLog = new TGTextButton(frmBut1, "Show log...");
   fBtnShowLog->SetToolTipText("Show Session log (opens log window)");
   frmBut1->AddFrame(fBtnShowLog, new TGLayoutHints(kLHintsLeft |
         kLHintsExpandX, 5, 5, 5, 5));
   fFA->AddFrame(frmBut1, new TGLayoutHints(kLHintsLeft | kLHintsBottom |
         kLHintsExpandX));

   // add "Commands" tab element
   tf = fTab->AddTab("Commands");
   fFC = new TGCompositeFrame(tf, 100, 100, kVerticalFrame);
   tf->AddFrame(fFC, new TGLayoutHints(kLHintsTop | kLHintsLeft |
         kLHintsExpandX | kLHintsExpandY));

   // add command line label and text entry
   TGCompositeFrame* frmCmd = new TGHorizontalFrame(fFC, 350, 100);
   frmCmd->SetCleanup(kDeepCleanup);
   frmCmd->AddFrame(new TGLabel(frmCmd, "Command Line :"),
         new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 5, 5, 15, 5));
   fCommandBuf = new TGTextBuffer(120);
   frmCmd->AddFrame(fCommandTxt = new TGTextEntry(frmCmd,
         fCommandBuf ),new TGLayoutHints(kLHintsLeft | kLHintsCenterY |
         kLHintsExpandX, 5, 5, 15, 5));
   fFC->AddFrame(frmCmd, new TGLayoutHints(kLHintsExpandX, 5, 5, 10, 5));
   // connect command line text entry to "return pressed" signal
   fCommandTxt->Connect("ReturnPressed()", "TSessionFrame", this,
         "OnCommandLine()");

   // check box for option "clear view"
   fClearCheck = new TGCheckButton(fFC, "Clear view after each command");
   fFC->AddFrame(fClearCheck,new TGLayoutHints(kLHintsLeft | kLHintsTop,
         10, 5, 5, 5));
   fClearCheck->SetState(kButtonUp);
   // add text view for redirected output
   fFC->AddFrame(new TGLabel(fFC, "Output :"),
      new TGLayoutHints(kLHintsLeft | kLHintsTop, 10, 5, 5, 5));
   fInfoTextView = new TGTextView(fFC, 330, 150, "", kSunkenFrame |
         kDoubleBorder);
   fFC->AddFrame(fInfoTextView, new TGLayoutHints(kLHintsLeft |
         kLHintsTop | kLHintsExpandX | kLHintsExpandY, 10, 10, 5, 5));

   // add "Packages" tab element
   tf = fTab->AddTab("Packages");
   fFB = new TGCompositeFrame(tf, 100, 100, kVerticalFrame);
   tf->AddFrame(fFB, new TGLayoutHints(kLHintsTop | kLHintsLeft |
         kLHintsExpandX | kLHintsExpandY));

   // new frame containing packages listbox and control buttons
   TGCompositeFrame* frmcanvas = new TGHorizontalFrame(fFB, 350, 100);

   // packages listbox
   fLBPackages = new TGListBox(frmcanvas);
   fLBPackages->Resize(80,150);
   fLBPackages->SetMultipleSelections(kFALSE);
   frmcanvas->AddFrame(fLBPackages, new TGLayoutHints(kLHintsExpandX |
         kLHintsExpandY, 5, 5, 5, 5));
   // control buttons frame
   TGCompositeFrame* frmBut2 = new TGVerticalFrame(frmcanvas, 150, 100);

   fChkMulti = new TGCheckButton(frmBut2, "Multiple Selection");
   fChkMulti->SetToolTipText("Enable multiple selection in the package list");
   frmBut2->AddFrame(fChkMulti, new TGLayoutHints(kLHintsLeft, 5, 5, 5, 5));

   fBtnAdd = new TGTextButton(frmBut2, "     Add...     ");
   fBtnAdd->SetToolTipText("Add a package to the list");
   frmBut2->AddFrame(fBtnAdd,new TGLayoutHints(kLHintsTop | kLHintsLeft |
         kLHintsExpandX, 5, 5, 5, 5));
   fBtnRemove = new TGTextButton(frmBut2, "Remove");
   fBtnRemove->SetToolTipText("Remove package from the list");
   frmBut2->AddFrame(fBtnRemove,new TGLayoutHints(kLHintsTop | kLHintsLeft |
         kLHintsExpandX, 5, 5, 5, 5));
   fBtnUp = new TGTextButton(frmBut2, "Move Up");
   fBtnUp->SetToolTipText("Move package one step upward in the list");
   frmBut2->AddFrame(fBtnUp,new TGLayoutHints(kLHintsTop | kLHintsLeft |
         kLHintsExpandX, 5, 5, 5, 5));
   fBtnDown = new TGTextButton(frmBut2, "Move Down");
   fBtnDown->SetToolTipText("Move package one step downward in the list");
   frmBut2->AddFrame(fBtnDown,new TGLayoutHints(kLHintsTop | kLHintsLeft |
         kLHintsExpandX, 5, 5, 5, 5));
   frmcanvas->AddFrame(frmBut2, new TGLayoutHints(kLHintsLeft | kLHintsCenterY |
         kLHintsExpandY));
   fFB->AddFrame(frmcanvas, new TGLayoutHints(kLHintsLeft | kLHintsTop |
         kLHintsExpandX | kLHintsExpandY));

   TGCompositeFrame* frmLeg = new TGHorizontalFrame(fFB, 300, 100);
   frmLeg->SetCleanup(kDeepCleanup);
   TGPicture *pic1 = (TGPicture *)fClient->GetPicture("package.xpm");
   TGIcon *icn1 = new TGIcon(frmLeg, pic1, pic1->GetWidth(), pic1->GetHeight());
   frmLeg->AddFrame(icn1, new TGLayoutHints(kLHintsLeft | kLHintsTop,
         5, 5, 0, 5));
   frmLeg->AddFrame(new TGLabel(frmLeg, ": Local"),
         new TGLayoutHints(kLHintsLeft | kLHintsTop, 0, 10, 0, 5));

   TGPicture *pic2 = (TGPicture *)fClient->GetPicture("package_delete.xpm");
   TGIcon *icn2 = new TGIcon(frmLeg, pic2, pic2->GetWidth(), pic2->GetHeight());
   frmLeg->AddFrame(icn2, new TGLayoutHints(kLHintsLeft | kLHintsTop,
         5, 5, 0, 5));
   frmLeg->AddFrame(new TGLabel(frmLeg, ": Uploaded"),
         new TGLayoutHints(kLHintsLeft | kLHintsTop, 0, 10, 0, 5));

   TGPicture *pic3 = (TGPicture *)fClient->GetPicture("package_add.xpm");
   TGIcon *icn3 = new TGIcon(frmLeg, pic3, pic3->GetWidth(), pic3->GetHeight());
   frmLeg->AddFrame(icn3, new TGLayoutHints(kLHintsLeft | kLHintsTop,
         5, 5, 0, 5));
   frmLeg->AddFrame(new TGLabel(frmLeg, ": Enabled"),
         new TGLayoutHints(kLHintsLeft | kLHintsTop, 0, 10, 0, 5));
   fFB->AddFrame(frmLeg, new TGLayoutHints(kLHintsLeft | kLHintsTop |
         kLHintsExpandX, 0, 0, 0, 0));

   TGCompositeFrame* frmBtn = new TGHorizontalFrame(fFB, 300, 100);
   frmBtn->SetCleanup(kDeepCleanup);
   frmBtn->AddFrame(fBtnUpload = new TGTextButton(frmBtn,
         " Upload "), new TGLayoutHints(kLHintsLeft | kLHintsExpandX |
         kLHintsCenterY, 5, 5, 5, 5));
   fBtnUpload->SetToolTipText("Upload selected package(s) to the server");
   frmBtn->AddFrame(fBtnEnable = new TGTextButton(frmBtn,
         " Enable "), new TGLayoutHints(kLHintsLeft | kLHintsExpandX |
         kLHintsCenterY, 5, 5, 5, 5));
   fBtnEnable->SetToolTipText("Enable selected package(s) on the server");
   frmBtn->AddFrame(fBtnDisable = new TGTextButton(frmBtn,
         " Disable "), new TGLayoutHints(kLHintsLeft | kLHintsExpandX |
         kLHintsCenterY, 5, 5, 5, 5));
   fBtnDisable->SetToolTipText("Disable selected package(s) on the server");
   frmBtn->AddFrame(fBtnClear = new TGTextButton(frmBtn,
         " Clear "), new TGLayoutHints(kLHintsLeft | kLHintsExpandX |
         kLHintsCenterY, 5, 5, 5, 5));
   fBtnClear->SetToolTipText("Clear all packages on the server");
   fFB->AddFrame(frmBtn, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));

   fBtnClear->SetEnabled(kFALSE);

   TGCompositeFrame* frmBtn3 = new TGHorizontalFrame(fFB, 300, 100);
   frmBtn3->SetCleanup(kDeepCleanup);
   fBtnShow = new TGTextButton(frmBtn3, "Show packages");
   fBtnShow->SetToolTipText("Show (list) available packages on the server");
   frmBtn3->AddFrame(fBtnShow,new TGLayoutHints(kLHintsCenterY | kLHintsLeft |
         kLHintsExpandX, 5, 5, 5, 5));
   fBtnShowEnabled = new TGTextButton(frmBtn3, "Show Enabled");
   fBtnShowEnabled->SetToolTipText("Show (list) enabled packages on the server");
   frmBtn3->AddFrame(fBtnShowEnabled,new TGLayoutHints(kLHintsCenterY | kLHintsLeft |
         kLHintsExpandX, 5, 5, 5, 5));
   fFB->AddFrame(frmBtn3, new TGLayoutHints(kLHintsExpandX, 0, 0, 0, 0));

   fChkEnable = new TGCheckButton(fFB, "Enable at session startup");
   fChkEnable->SetToolTipText("Enable packages on the server at startup time");
   fFB->AddFrame(fChkEnable, new TGLayoutHints(kLHintsLeft, 5, 5, 5, 5));

   // add "DataSets" tab element
   tf = fTab->AddTab("DataSets");
   fFE = new TGCompositeFrame(tf, 100, 100, kVerticalFrame);
   tf->AddFrame(fFE, new TGLayoutHints(kLHintsTop | kLHintsLeft |
         kLHintsExpandX | kLHintsExpandY));

   // new frame containing datasets treeview and control buttons
   TGCompositeFrame* frmdataset = new TGHorizontalFrame(fFE, 350, 100);

   // datasets list tree
   fDSetView = new TGCanvas(frmdataset, 200, 200, kSunkenFrame | kDoubleBorder);
   frmdataset->AddFrame(fDSetView, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY,
         5, 5, 5, 5));
   fDataSetTree = new TGListTree(fDSetView, kHorizontalFrame);
   fDataSetTree->AddItem(0, "DataSets");

   // control buttons frame
   TGCompositeFrame* frmBut3 = new TGVerticalFrame(frmdataset, 150, 100);

   fBtnUploadDSet = new TGTextButton(frmBut3, "     Upload...     ");
   fBtnUploadDSet->SetToolTipText("Upload a dataset to the cluster");
   frmBut3->AddFrame(fBtnUploadDSet, new TGLayoutHints(kLHintsTop | kLHintsLeft |
         kLHintsExpandX, 5, 5, 5, 5));
   fBtnRemoveDSet = new TGTextButton(frmBut3, "Remove");
   fBtnRemoveDSet->SetToolTipText("Remove dataset from the cluster");
   frmBut3->AddFrame(fBtnRemoveDSet,new TGLayoutHints(kLHintsTop | kLHintsLeft |
         kLHintsExpandX, 5, 5, 5, 5));
   fBtnVerifyDSet = new TGTextButton(frmBut3, "Verify");
   fBtnVerifyDSet->SetToolTipText("Verify dataset on the cluster");
   frmBut3->AddFrame(fBtnVerifyDSet,new TGLayoutHints(kLHintsTop | kLHintsLeft |
         kLHintsExpandX, 5, 5, 5, 5));
   fBtnRefresh = new TGTextButton(frmBut3, "Refresh List");
   fBtnRefresh->SetToolTipText("Refresh List of DataSet/Files present on the cluster");
   frmBut3->AddFrame(fBtnRefresh,new TGLayoutHints(kLHintsTop | kLHintsLeft |
         kLHintsExpandX, 5, 5, 15, 5));

   frmdataset->AddFrame(frmBut3, new TGLayoutHints(kLHintsLeft | kLHintsCenterY |
         kLHintsExpandY, 5, 5, 5, 0));

   fFE->AddFrame(frmdataset, new TGLayoutHints(kLHintsLeft | kLHintsTop |
         kLHintsExpandX | kLHintsExpandY));

   // add "Options" tab element
   tf = fTab->AddTab("Options");
   fFD = new TGCompositeFrame(tf, 100, 100, kVerticalFrame);
   tf->AddFrame(fFD, new TGLayoutHints(kLHintsTop | kLHintsLeft |
         kLHintsExpandX | kLHintsExpandY));

   // add Log Level label and text entry
   TGCompositeFrame* frmLog = new TGHorizontalFrame(fFD, 310, 100, kFixedWidth);
   frmLog->SetCleanup(kDeepCleanup);
   frmLog->AddFrame(fApplyLogLevel = new TGTextButton(frmLog,
         "        Apply        "), new TGLayoutHints(kLHintsRight |
         kLHintsCenterY, 10, 5, 5, 5));
   fApplyLogLevel->SetToolTipText("Apply currently selected log level");
   fLogLevel = new TGNumberEntry(frmLog, 0, 5, 5, TGNumberFormat::kNESInteger,
         TGNumberFormat::kNEANonNegative, TGNumberFormat::kNELLimitMinMax, 0, 5);
   frmLog->AddFrame(fLogLevel, new TGLayoutHints(kLHintsRight |
         kLHintsCenterY, 5, 5, 5, 5));
   frmLog->AddFrame(new TGLabel(frmLog, "Log Level :"),
         new TGLayoutHints(kLHintsRight | kLHintsCenterY, 5, 5, 5, 5));
   fFD->AddFrame(frmLog, new TGLayoutHints(kLHintsLeft, 5, 5, 15, 5));

   // add Parallel Nodes label and text entry
   TGCompositeFrame* frmPar = new TGHorizontalFrame(fFD, 310, 100, kFixedWidth);
   frmPar->SetCleanup(kDeepCleanup);
   frmPar->AddFrame(fApplyParallel = new TGTextButton(frmPar,
         "        Apply        "), new TGLayoutHints(kLHintsRight |
         kLHintsCenterY, 10, 5, 5, 5));
   fApplyParallel->SetToolTipText("Apply currently selected parallel nodes");
   fTxtParallel = new TGTextEntry(frmPar);
   fTxtParallel->SetAlignment(kTextRight);
   fTxtParallel->SetText("99999");
   fTxtParallel->Resize(fLogLevel->GetDefaultWidth(), fTxtParallel->GetDefaultHeight());
   frmPar->AddFrame(fTxtParallel, new TGLayoutHints(kLHintsRight |
         kLHintsCenterY, 5, 5, 5, 5));
   frmPar->AddFrame(new TGLabel(frmPar, "Set Parallel Nodes :"),
         new TGLayoutHints(kLHintsRight | kLHintsCenterY, 5, 5, 5, 5));
   fFD->AddFrame(frmPar, new TGLayoutHints(kLHintsLeft, 5, 5, 5, 5));

   // connect button actions to functions
   fBtnShowLog->Connect("Clicked()", "TSessionFrame", this,
         "OnBtnShowLogClicked()");
   fBtnNewQuery->Connect("Clicked()", "TSessionFrame", this,
         "OnBtnNewQueryClicked()");
   fBtnGetQueries->Connect("Clicked()", "TSessionFrame", this,
         "OnBtnGetQueriesClicked()");

   fChkEnable->Connect("Toggled(Bool_t)", "TSessionFrame", this,
         "OnStartupEnable(Bool_t)");
   fChkMulti->Connect("Toggled(Bool_t)", "TSessionFrame", this,
         "OnMultipleSelection(Bool_t)");
   fBtnAdd->Connect("Clicked()", "TSessionFrame", this,
         "OnBtnAddClicked()");
   fBtnRemove->Connect("Clicked()", "TSessionFrame", this,
         "OnBtnRemoveClicked()");
   fBtnUp->Connect("Clicked()", "TSessionFrame", this,
         "OnBtnUpClicked()");
   fBtnDown->Connect("Clicked()", "TSessionFrame", this,
         "OnBtnDownClicked()");
   fApplyLogLevel->Connect("Clicked()", "TSessionFrame", this,
         "OnApplyLogLevel()");
   fApplyParallel->Connect("Clicked()", "TSessionFrame", this,
         "OnApplyParallel()");
   fBtnUpload->Connect("Clicked()", "TSessionFrame", this,
         "OnUploadPackages()");
   fBtnEnable->Connect("Clicked()", "TSessionFrame", this,
         "OnEnablePackages()");
   fBtnDisable->Connect("Clicked()", "TSessionFrame", this,
         "OnDisablePackages()");
   fBtnClear->Connect("Clicked()", "TSessionFrame", this,
         "OnClearPackages()");
   fBtnShowEnabled->Connect("Clicked()", "TSessionViewer", fViewer,
         "ShowEnabledPackages()");
   fBtnShow->Connect("Clicked()", "TSessionViewer", fViewer,
         "ShowPackages()");

   fBtnUploadDSet->Connect("Clicked()", "TSessionFrame", this,
         "OnBtnUploadDSet()");
   fBtnRemoveDSet->Connect("Clicked()", "TSessionFrame", this,
         "OnBtnRemoveDSet()");
   fBtnVerifyDSet->Connect("Clicked()", "TSessionFrame", this,
         "OnBtnVerifyDSet()");
   fBtnRefresh->Connect("Clicked()", "TSessionFrame", this,
         "UpdateListOfDataSets()");
}

////////////////////////////////////////////////////////////////////////////////
/// Display information on current session.

void TSessionFrame::ProofInfos()
{
   TString buf;

   // if local session
   if (fViewer->GetActDesc()->fLocal) {
      buf.Form("*** Local Session on %s ***", gSystem->HostName());
      fInfoLine[0]->SetText(buf.Data());
      UserGroup_t *userGroup = gSystem->GetUserInfo();
      fInfoLine[1]->SetText("User :");
      if (userGroup) buf.Form("%s", userGroup->fRealName.Data());
      fInfoLine[2]->SetText(buf.Data());
      fInfoLine[3]->SetText("Working directory :");
      buf.Form("%s", gSystem->WorkingDirectory());
      fInfoLine[4]->SetText(buf.Data());
      fInfoLine[5]->SetText(" ");
      fInfoLine[6]->SetText(" ");
      fInfoLine[7]->SetText(" ");
      fInfoLine[8]->SetText(" ");
      fInfoLine[9]->SetText(" ");
      fInfoLine[10]->SetText(" ");
      fInfoLine[11]->SetText(" ");
      fInfoLine[12]->SetText(" ");
      fInfoLine[13]->SetText(" ");
      fInfoLine[14]->SetText(" ");
      fInfoLine[15]->SetText(" ");
      fInfoLine[16]->SetText(" ");
      fInfoLine[17]->SetText(" ");
      fInfoLine[18]->SetText(" ");
      if (userGroup) delete userGroup;
      Layout();
      Resize(GetDefaultSize());
      return;
   }
   // return if not a valid Proof session
   if (!fViewer->GetActDesc()->fConnected ||
       !fViewer->GetActDesc()->fAttached ||
       !fViewer->GetActDesc()->fProof ||
       !fViewer->GetActDesc()->fProof->IsValid())
      return;

   if (!fViewer->GetActDesc()->fProof->IsMaster()) {
      if (fViewer->GetActDesc()->fProof->IsParallel())
         buf.Form("*** Connected to %s (parallel mode, %d workers) ***",
               fViewer->GetActDesc()->fProof->GetMaster(),
               fViewer->GetActDesc()->fProof->GetParallel());
      else
         buf.Form("*** Connected to %s (sequential mode) ***",
               fViewer->GetActDesc()->fProof->GetMaster());
      fInfoLine[0]->SetText(buf.Data());
      fInfoLine[1]->SetText("Port number : ");
      buf.Form("%d", fViewer->GetActDesc()->fProof->GetPort());
      fInfoLine[2]->SetText(buf.Data());
      fInfoLine[3]->SetText("User : ");
      buf.Form("%s", fViewer->GetActDesc()->fProof->GetUser());
      fInfoLine[4]->SetText(buf.Data());
      fInfoLine[5]->SetText("Client protocol version : ");
      buf.Form("%d", fViewer->GetActDesc()->fProof->GetClientProtocol());
      fInfoLine[6]->SetText(buf.Data());
      fInfoLine[7]->SetText("Remote protocol version : ");
      buf.Form("%d", fViewer->GetActDesc()->fProof->GetRemoteProtocol());
      fInfoLine[8]->SetText(buf.Data());
      fInfoLine[9]->SetText("Log level : ");
      buf.Form("%d", fViewer->GetActDesc()->fProof->GetLogLevel());
      fInfoLine[10]->SetText(buf.Data());
      fInfoLine[11]->SetText("Session unique tag : ");
      buf.Form("%s", fViewer->GetActDesc()->fProof->IsValid() ?
            fViewer->GetActDesc()->fProof->GetSessionTag() : " ");
      fInfoLine[12]->SetText(buf.Data());
      fInfoLine[13]->SetText("Total MB's processed :");
      buf.Form("%.2f", float(fViewer->GetActDesc()->fProof->GetBytesRead())/(1024*1024));
      fInfoLine[14]->SetText(buf.Data());
      fInfoLine[15]->SetText("Total real time used (s) :");
      buf.Form("%.3f", fViewer->GetActDesc()->fProof->GetRealTime());
      fInfoLine[16]->SetText(buf.Data());
      fInfoLine[17]->SetText("Total CPU time used (s) :");
      buf.Form("%.3f", fViewer->GetActDesc()->fProof->GetCpuTime());
      fInfoLine[18]->SetText(buf.Data());
   }
   else {
      if (fViewer->GetActDesc()->fProof->IsParallel())
         buf.Form("*** Master server %s (parallel mode, %d workers) ***",
               fViewer->GetActDesc()->fProof->GetMaster(),
               fViewer->GetActDesc()->fProof->GetParallel());
      else
         buf.Form("*** Master server %s (sequential mode) ***",
               fViewer->GetActDesc()->fProof->GetMaster());
      fInfoLine[0]->SetText(buf.Data());
      fInfoLine[1]->SetText("Port number : ");
      buf.Form("%d", fViewer->GetActDesc()->fProof->GetPort());
      fInfoLine[2]->SetText(buf.Data());
      fInfoLine[3]->SetText("User : ");
      buf.Form("%s", fViewer->GetActDesc()->fProof->GetUser());
      fInfoLine[4]->SetText(buf.Data());
      fInfoLine[5]->SetText("Protocol version : ");
      buf.Form("%d", fViewer->GetActDesc()->fProof->GetClientProtocol());
      fInfoLine[6]->SetText(buf.Data());
      fInfoLine[7]->SetText("Image name : ");
      buf.Form("%s",fViewer->GetActDesc()->fProof->GetImage());
      fInfoLine[8]->SetText(buf.Data());
      fInfoLine[9]->SetText("Config directory : ");
      buf.Form("%s", fViewer->GetActDesc()->fProof->GetConfDir());
      fInfoLine[10]->SetText(buf.Data());
      fInfoLine[11]->SetText("Config file : ");
      buf.Form("%s", fViewer->GetActDesc()->fProof->GetConfFile());
      fInfoLine[12]->SetText(buf.Data());
      fInfoLine[13]->SetText("Total MB's processed :");
      buf.Form("%.2f", float(fViewer->GetActDesc()->fProof->GetBytesRead())/(1024*1024));
      fInfoLine[14]->SetText(buf.Data());
      fInfoLine[15]->SetText("Total real time used (s) :");
      buf.Form("%.3f", fViewer->GetActDesc()->fProof->GetRealTime());
      fInfoLine[16]->SetText(buf.Data());
      fInfoLine[17]->SetText("Total CPU time used (s) :");
      buf.Form("%.3f", fViewer->GetActDesc()->fProof->GetCpuTime());
      fInfoLine[18]->SetText(buf.Data());
   }
   Layout();
   Resize(GetDefaultSize());
}

////////////////////////////////////////////////////////////////////////////////
/// Open Upload Dataset dialog.

void TSessionFrame::OnBtnUploadDSet()
{
   if (fViewer->IsBusy())
      return;
   if (fViewer->GetActDesc()->fLocal) return;
   new TUploadDataSetDlg(fViewer, 450, 360);
}

////////////////////////////////////////////////////////////////////////////////
/// Update list of dataset present on the cluster.

void TSessionFrame::UpdateListOfDataSets()
{
   TObjString *dsetname;
   TFileInfo  *dsetfilename;
   // cleanup the list
   fDataSetTree->DeleteChildren(fDataSetTree->GetFirstItem());
   if (fViewer->GetActDesc()->fConnected && fViewer->GetActDesc()->fAttached &&
       fViewer->GetActDesc()->fProof && fViewer->GetActDesc()->fProof->IsValid() &&
       fViewer->GetActDesc()->fProof->IsParallel()) {

      const TGPicture *dseticon = fClient->GetPicture("rootdb_t.xpm");
      // ask for the list of datasets

      // TODO: is now returning a TMap; viewer has to be adapted
      TList *dsetlist = 0; //fViewer->GetActDesc()->fProof->GetDataSets();
      // coverity[dead_error_condition]: to be changed for TMap usage
      if (dsetlist) {
         TGListTreeItem *dsetitem;
         fDataSetTree->OpenItem(fDataSetTree->GetFirstItem());
         TIter nextdset(dsetlist);
         while ((dsetname = (TObjString *)nextdset())) {
            if (!fDataSetTree->FindItemByObj(fDataSetTree->GetFirstItem(), dsetname)) {
               // add the dataset in the tree
               dsetitem = fDataSetTree->AddItem(fDataSetTree->GetFirstItem(),
                                        dsetname->GetName(), dsetname);
               // ask for the list of files in the dataset
               TList *dsetfilelist = fViewer->GetActDesc()->fProof->GetDataSet(
                                                            dsetname->GetName())->GetList();
               if(dsetfilelist) {
                  TIter nextdsetfile(dsetfilelist);
                  while ((dsetfilename = (TFileInfo *)nextdsetfile())) {
                     if (! fDataSetTree->FindItemByObj(dsetitem, dsetfilename)) {
                        // if not already in, add the file name in the tree
                        fDataSetTree->AddItem(dsetitem,
                           dsetfilename->GetFirstUrl()->GetUrl(),
                           dsetfilename, dseticon, dseticon);
                     }
                  }
                  // open the dataset item in order to show the files
                  fDataSetTree->OpenItem(dsetitem);
               }
            }
         }
      }
   }
   // refresh list tree
   fClient->NeedRedraw(fDataSetTree);
}

////////////////////////////////////////////////////////////////////////////////
/// Remove dataset from the list and from the cluster.

void TSessionFrame::OnBtnRemoveDSet()
{
   TGListTreeItem *item;
   TObjString *obj = 0;
   if (fViewer->GetActDesc()->fLocal) return;

   item = fDataSetTree->GetSelected();
   if (!item) return;
   if (item->GetParent() == 0) return;
   if (item->GetParent() == fDataSetTree->GetFirstItem()) {
      // Dataset itself
      obj = (TObjString *)item->GetUserData();
   }
   else if (item->GetParent()->GetParent() == fDataSetTree->GetFirstItem()) {
      // One file of the dataset
      obj = (TObjString *)item->GetParent()->GetUserData();
   }

   // if valid Proof session, set parallel slaves
   if (obj && fViewer->GetActDesc()->fProof &&
      fViewer->GetActDesc()->fProof->IsValid()) {
      fViewer->GetActDesc()->fProof->RemoveDataSet(obj->GetName());
      UpdateListOfDataSets();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Verify that the files in the selected dataset are present on the cluster.

void TSessionFrame::OnBtnVerifyDSet()
{
   TGListTreeItem *item;
   TObjString *obj = 0;
   if (fViewer->GetActDesc()->fLocal) return;

   item = fDataSetTree->GetSelected();
   if (!item) return;
   if (item->GetParent() == 0) return;
   if (item->GetParent() == fDataSetTree->GetFirstItem()) {
      // Dataset itself
      obj = (TObjString *)item->GetUserData();
   }
   else if (item->GetParent()->GetParent() == fDataSetTree->GetFirstItem()) {
      // One file of the dataset
      obj = (TObjString *)item->GetParent()->GetUserData();
   }

   // if valid Proof session, set parallel slaves
   if (obj && fViewer->GetActDesc()->fProof &&
      fViewer->GetActDesc()->fProof->IsValid()) {
      fViewer->GetActDesc()->fProof->VerifyDataSet(obj->GetName());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Apply selected log level on current session.

void TSessionFrame::OnApplyLogLevel()
{
   // if local session, do nothing
   if (fViewer->GetActDesc()->fLocal) return;
   // if valid Proof session, set log level
   if (fViewer->GetActDesc()->fProof &&
      fViewer->GetActDesc()->fProof->IsValid()) {
      fViewer->GetActDesc()->fLogLevel = fLogLevel->GetIntNumber();
      fViewer->GetActDesc()->fProof->SetLogLevel(fViewer->GetActDesc()->fLogLevel);
   }
   fViewer->GetSessionFrame()->ProofInfos();
}

////////////////////////////////////////////////////////////////////////////////
/// Apply selected number of workers on current Proof session.

void TSessionFrame::OnApplyParallel()
{
   // if local session, do nothing
   if (fViewer->GetActDesc()->fLocal) return;
   // if valid Proof session, set parallel slaves
   if (fViewer->GetActDesc()->fProof &&
      fViewer->GetActDesc()->fProof->IsValid()) {
      Int_t nodes = atoi(fTxtParallel->GetText());
      fViewer->GetActDesc()->fProof->SetParallel(nodes);
   }
   fViewer->GetSessionFrame()->ProofInfos();
}

////////////////////////////////////////////////////////////////////////////////
/// Handle multiple selection check button.

void TSessionFrame::OnMultipleSelection(Bool_t on)
{
   fLBPackages->SetMultipleSelections(on);
}

////////////////////////////////////////////////////////////////////////////////
/// Handle multiple selection check button.

void TSessionFrame::OnStartupEnable(Bool_t on)
{
   if (fViewer->GetActDesc())
      fViewer->GetActDesc()->fAutoEnable = on;
}

////////////////////////////////////////////////////////////////////////////////
/// Update list of packages.

void TSessionFrame::UpdatePackages()
{
   TPackageDescription *package;
   const TGPicture *pict;
   fLBPackages->RemoveEntries(0, fLBPackages->GetNumberOfEntries());
   TIter next(fViewer->GetActDesc()->fPackages);
   while ((package = (TPackageDescription *)next())) {
      if (package->fEnabled)
         pict = fClient->GetPicture("package_add.xpm");
      else if (package->fUploaded)
         pict = fClient->GetPicture("package_delete.xpm");
      else
         pict = fClient->GetPicture("package.xpm");
      TGIconLBEntry *entry = new TGIconLBEntry(fLBPackages->GetContainer(),
                                    package->fId, package->fPathName, pict);
      fLBPackages->AddEntry(entry, new TGLayoutHints(kLHintsExpandX | kLHintsTop));
   }
   fLBPackages->Layout();
   fClient->NeedRedraw(fLBPackages->GetContainer());
}

////////////////////////////////////////////////////////////////////////////////
/// Upload selected package(s) to the current session.

void TSessionFrame::OnUploadPackages()
{
   // if local session, do nothing
   if (fViewer->GetActDesc()->fLocal) return;
   // if valid Proof session, upload packages
   if (fViewer->GetActDesc()->fProof &&
      fViewer->GetActDesc()->fProof->IsValid()) {
      TObject *obj;
      TList selected;
      fLBPackages->GetSelectedEntries(&selected);
      TIter next(&selected);
      while ((obj = next())) {
         TString name = obj->GetTitle();
         if (fViewer->GetActDesc()->fProof->UploadPackage(name.Data()) != 0)
            Error("Submit", "Upload package failed");
         else {
            TObject *o = fViewer->GetActDesc()->fPackages->FindObject(gSystem->BaseName(name));
            if (!o) continue;
            TPackageDescription *package =
               dynamic_cast<TPackageDescription *>(o);
            if (package) {
               package->fUploaded = kTRUE;
               ((TGIconLBEntry *)obj)->SetPicture(
                     fClient->GetPicture("package_delete.xpm"));
            }
         }
      }
      UpdatePackages();
   }
   fLBPackages->Layout();
   fClient->NeedRedraw(fLBPackages->GetContainer());
}

////////////////////////////////////////////////////////////////////////////////
/// Enable selected package(s) in the current session.

void TSessionFrame::OnEnablePackages()
{
   // if local session, do nothing
   if (fViewer->GetActDesc()->fLocal) return;
   // if valid Proof session, enable packages
   if (fViewer->GetActDesc()->fProof &&
      fViewer->GetActDesc()->fProof->IsValid()) {
      TObject *obj;
      TList selected;
      fBtnEnable->SetState(kButtonDisabled);
      fLBPackages->GetSelectedEntries(&selected);
      TIter next(&selected);
      while ((obj = next())) {
         TString name = obj->GetTitle();
         TObject *o = fViewer->GetActDesc()->fPackages->FindObject(gSystem->BaseName(name));
         if (!o) continue;
         TPackageDescription *package =
            dynamic_cast<TPackageDescription *>(o);
         if (package) {
            if (!package->fUploaded) {
               if (fViewer->GetActDesc()->fProof->UploadPackage(name.Data()) != 0)
                  Error("Submit", "Upload package failed");
               else {
                  package->fUploaded = kTRUE;
                  ((TGIconLBEntry *)obj)->SetPicture(
                        fClient->GetPicture("package_delete.xpm"));
               }
            }
         }
         if (fViewer->GetActDesc()->fProof->EnablePackage(name) != 0)
            Error("Submit", "Enable package failed");
         else {
            package->fEnabled = kTRUE;
            ((TGIconLBEntry *)obj)->SetPicture(fClient->GetPicture("package_add.xpm"));
         }
      }
      UpdatePackages();
      fBtnEnable->SetState(kButtonUp);
   }
   fLBPackages->Layout();
   fClient->NeedRedraw(fLBPackages->GetContainer());
}

////////////////////////////////////////////////////////////////////////////////
/// Disable selected package(s) in the current session.

void TSessionFrame::OnDisablePackages()
{
   // if local session, do nothing
   if (fViewer->GetActDesc()->fLocal) return;
   // if valid Proof session, disable (clear) packages
   if (fViewer->GetActDesc()->fProof &&
      fViewer->GetActDesc()->fProof->IsValid()) {
      TObject *obj;
      TList selected;
      fLBPackages->GetSelectedEntries(&selected);
      TIter next(&selected);
      while ((obj = next())) {
         TString name = obj->GetTitle();
         if (fViewer->GetActDesc()->fProof->ClearPackage(name) != 0)
            Error("Submit", "Clear package failed");
         else {
            TObject *o = fViewer->GetActDesc()->fPackages->FindObject(gSystem->BaseName(name));
            if (!o) continue;
            TPackageDescription *package =
               dynamic_cast<TPackageDescription *>(o);
            if (package) {
               package->fEnabled = kFALSE;
               package->fUploaded = kFALSE;
               ((TGIconLBEntry *)obj)->SetPicture(fClient->GetPicture("package.xpm"));
            }
         }
      }
      UpdatePackages();
   }
   fLBPackages->Layout();
   fClient->NeedRedraw(fLBPackages->GetContainer());
}

////////////////////////////////////////////////////////////////////////////////
/// Clear (disable) all packages in the current session.

void TSessionFrame::OnClearPackages()
{
   TPackageDescription *package;
   // if local session, do nothing
   if (fViewer->GetActDesc()->fLocal) return;
   // if valid Proof session, clear packages
   if (fViewer->GetActDesc()->fProof &&
      fViewer->GetActDesc()->fProof->IsValid()) {
      if (fViewer->GetActDesc()->fProof->ClearPackages() != 0)
         Error("Submit", "Clear packages failed");
      else {
         TIter next(fViewer->GetActDesc()->fPackages);
         while ((package = (TPackageDescription *)next())) {
            package->fEnabled = kFALSE;
         }
      }
   }
   fLBPackages->Layout();
   fClient->NeedRedraw(fLBPackages->GetContainer());
}

////////////////////////////////////////////////////////////////////////////////
/// Open file dialog and add selected package file to the list.

void TSessionFrame::OnBtnAddClicked()
{
   if (fViewer->IsBusy())
      return;
   TGFileInfo fi;
   TPackageDescription *package;
   TGIconLBEntry *entry;
   fi.fFileTypes = pkgtypes;
   new TGFileDialog(fClient->GetRoot(), fViewer, kFDOpen, &fi);
   if (fi.fMultipleSelection && fi.fFileNamesList) {
      TObjString *el;
      TIter next(fi.fFileNamesList);
      while ((el = (TObjString *) next())) {
         package = new TPackageDescription;
         package->fName = gSystem->BaseName(gSystem->UnixPathName(el->GetString()));
         package->fPathName = gSystem->UnixPathName(el->GetString());
         package->fId   = fViewer->GetActDesc()->fPackages->GetEntries();
         package->fUploaded = kFALSE;
         package->fEnabled = kFALSE;
         fViewer->GetActDesc()->fPackages->Add((TObject *)package);
         entry = new TGIconLBEntry(fLBPackages->GetContainer(),
                                   package->fId, package->fPathName,
                                   fClient->GetPicture("package.xpm"));
         fLBPackages->AddEntry(entry, new TGLayoutHints(kLHintsExpandX | kLHintsTop));
      }
   }
   else if (fi.fFilename) {
      package = new TPackageDescription;
      package->fName = gSystem->BaseName(gSystem->UnixPathName(fi.fFilename));
      package->fPathName = gSystem->UnixPathName(fi.fFilename);
      package->fId   = fViewer->GetActDesc()->fPackages->GetEntries();
      package->fUploaded = kFALSE;
      package->fEnabled = kFALSE;
      fViewer->GetActDesc()->fPackages->Add((TObject *)package);
      entry = new TGIconLBEntry(fLBPackages->GetContainer(),
                                package->fId, package->fPathName,
                                fClient->GetPicture("package.xpm"));
      fLBPackages->AddEntry(entry, new TGLayoutHints(kLHintsExpandX | kLHintsTop));
   }
   fLBPackages->Layout();
   fClient->NeedRedraw(fLBPackages->GetContainer());
}

////////////////////////////////////////////////////////////////////////////////
/// Remove selected package from the list.

void TSessionFrame::OnBtnRemoveClicked()
{
   TPackageDescription *package;
   const TGPicture *pict;
   Int_t pos = fLBPackages->GetSelected();
   fLBPackages->RemoveEntries(0, fLBPackages->GetNumberOfEntries());
   fViewer->GetActDesc()->fPackages->Remove(
         fViewer->GetActDesc()->fPackages->At(pos));
   Int_t id = 0;
   TIter next(fViewer->GetActDesc()->fPackages);
   while ((package = (TPackageDescription *)next())) {
      package->fId = id;
      id++;
      if (package->fEnabled)
         pict = fClient->GetPicture("package_add.xpm");
      else if (package->fUploaded)
         pict = fClient->GetPicture("package_delete.xpm");
      else
         pict = fClient->GetPicture("package.xpm");
      TGIconLBEntry *entry = new TGIconLBEntry(fLBPackages->GetContainer(),
                                    package->fId, package->fPathName, pict);
      fLBPackages->AddEntry(entry, new TGLayoutHints(kLHintsExpandX | kLHintsTop));
   }
   fLBPackages->Layout();
   fClient->NeedRedraw(fLBPackages->GetContainer());
}

////////////////////////////////////////////////////////////////////////////////
/// Move selected package entry one position up in the list.

void TSessionFrame::OnBtnUpClicked()
{
   TPackageDescription *package;
   const TGPicture *pict;
   Int_t pos = fLBPackages->GetSelected();
   if (pos <= 0) return;
   fLBPackages->RemoveEntries(0, fLBPackages->GetNumberOfEntries());
   package = (TPackageDescription *)fViewer->GetActDesc()->fPackages->At(pos);
   fViewer->GetActDesc()->fPackages->Remove(
         fViewer->GetActDesc()->fPackages->At(pos));
   package->fId -= 1;
   fViewer->GetActDesc()->fPackages->AddAt(package, package->fId);
   Int_t id = 0;
   TIter next(fViewer->GetActDesc()->fPackages);
   while ((package = (TPackageDescription *)next())) {
      package->fId = id;
      id++;
      if (package->fEnabled)
         pict = fClient->GetPicture("package_add.xpm");
      else if (package->fUploaded)
         pict = fClient->GetPicture("package_delete.xpm");
      else
         pict = fClient->GetPicture("package.xpm");
      TGIconLBEntry *entry = new TGIconLBEntry(fLBPackages->GetContainer(),
                                    package->fId, package->fPathName, pict);
      fLBPackages->AddEntry(entry, new TGLayoutHints(kLHintsExpandX | kLHintsTop));
   }
   fLBPackages->Select(pos-1);
   fLBPackages->Layout();
   fClient->NeedRedraw(fLBPackages->GetContainer());
}

////////////////////////////////////////////////////////////////////////////////
/// Move selected package entry one position down in the list.

void TSessionFrame::OnBtnDownClicked()
{
   TPackageDescription *package;
   const TGPicture *pict;
   Int_t pos = fLBPackages->GetSelected();
   if (pos == -1 || pos == fViewer->GetActDesc()->fPackages->GetEntries()-1)
      return;
   fLBPackages->RemoveEntries(0, fLBPackages->GetNumberOfEntries());
   package = (TPackageDescription *)fViewer->GetActDesc()->fPackages->At(pos);
   fViewer->GetActDesc()->fPackages->Remove(
         fViewer->GetActDesc()->fPackages->At(pos));
   package->fId += 1;
   fViewer->GetActDesc()->fPackages->AddAt(package, package->fId);
   Int_t id = 0;
   TIter next(fViewer->GetActDesc()->fPackages);
   while ((package = (TPackageDescription *)next())) {
      package->fId = id;
      id++;
      if (package->fEnabled)
         pict = fClient->GetPicture("package_add.xpm");
      else if (package->fUploaded)
         pict = fClient->GetPicture("package_delete.xpm");
      else
         pict = fClient->GetPicture("package.xpm");
      TGIconLBEntry *entry = new TGIconLBEntry(fLBPackages->GetContainer(),
                                    package->fId, package->fPathName, pict);
      fLBPackages->AddEntry(entry, new TGLayoutHints(kLHintsExpandX | kLHintsTop));
   }
   fLBPackages->Select(pos+1);
   fLBPackages->Layout();
   fClient->NeedRedraw(fLBPackages->GetContainer());
}

////////////////////////////////////////////////////////////////////////////////
/// Disconnect from current Proof session.

void TSessionFrame::OnBtnDisconnectClicked()
{
   // if local session, do nothing
   if (fViewer->GetActDesc()->fLocal) return;
   // if valid Proof session, disconnect (close)
   if (fViewer->GetActDesc()->fAttached &&
       fViewer->GetActDesc()->fProof &&
       fViewer->GetActDesc()->fProof->IsValid()) {
      fViewer->GetActDesc()->fProof->Detach();
   }
   // reset connected flag
   fViewer->GetActDesc()->fAttached = kFALSE;
   fViewer->GetActDesc()->fProof = 0;
   // disable animation timer
   fViewer->DisableTimer();
   // change list tree item picture to disconnected pixmap
   TGListTreeItem *item = fViewer->GetSessionHierarchy()->FindChildByData(
                           fViewer->GetSessionItem(), fViewer->GetActDesc());
   if (item) {
      item->SetPictures(fViewer->GetProofDisconPict(),
                        fViewer->GetProofDisconPict());
   }
   // update viewer
   fViewer->OnListTreeClicked(fViewer->GetSessionHierarchy()->GetSelected(),
                              1, 0, 0);
   fViewer->GetSessionHierarchy()->ClearViewPort();
   fClient->NeedRedraw(fViewer->GetSessionHierarchy());
   fViewer->GetStatusBar()->SetText("", 1);
}

////////////////////////////////////////////////////////////////////////////////
/// Show session log.

void TSessionFrame::OnBtnShowLogClicked()
{
   fViewer->ShowLog(0);
}

////////////////////////////////////////////////////////////////////////////////
/// Call "New Query" Dialog.

void TSessionFrame::OnBtnNewQueryClicked()
{
   TNewQueryDlg *dlg = new TNewQueryDlg(fViewer, 350, 310);
   dlg->Popup();
}

////////////////////////////////////////////////////////////////////////////////
/// Get list of queries from current Proof server and populate the list tree.

void TSessionFrame::OnBtnGetQueriesClicked()
{
   TList *lqueries = 0;
   TQueryResult *query = 0;
   TQueryDescription *newquery = 0, *lquery = 0;
   if (fViewer->GetActDesc()->fProof &&
       fViewer->GetActDesc()->fProof->IsValid()) {
      lqueries = fViewer->GetActDesc()->fProof->GetListOfQueries();
   }
   if (lqueries) {
      TIter nextp(lqueries);
      // loop over list of queries received from Proof server
      while ((query = (TQueryResult *)nextp())) {
         // create new query description
         newquery = new TQueryDescription();
         newquery->fReference = TString::Format("%s:%s", query->GetTitle(),
                                                query->GetName());
         // check in our tree if it is already there
         TGListTreeItem *item =
            fViewer->GetSessionHierarchy()->FindChildByData(
                  fViewer->GetSessionItem(), fViewer->GetActDesc());
         // if already there, skip
         if (fViewer->GetSessionHierarchy()->FindChildByName(item,
             newquery->fReference.Data()))
            continue;
         // check also in our query description list
         Bool_t found = kFALSE;
         TIter nextp2(fViewer->GetActDesc()->fQueries);
         while ((lquery = (TQueryDescription *)nextp2())) {
            if (lquery->fReference.CompareTo(newquery->fReference) == 0) {
               found = kTRUE;
               break;
            }
         }
         if (found) continue;
         // build new query description with infos from Proof
         newquery->fStatus = query->IsFinalized() ?
               TQueryDescription::kSessionQueryFinalized :
               (TQueryDescription::ESessionQueryStatus)query->GetStatus();
         newquery->fSelectorString  = query->GetSelecImp()->GetName();
         newquery->fQueryName       = TString::Format("%s:%s", query->GetTitle(),
                                                      query->GetName());
         newquery->fOptions         = query->GetOptions();
         newquery->fEventList       = "";
         newquery->fNbFiles         = 0;
         newquery->fNoEntries       = query->GetEntries();
         newquery->fFirstEntry      = query->GetFirst();
         newquery->fResult          = query;
         newquery->fChain           = 0;
         fViewer->GetActDesc()->fQueries->Add((TObject *)newquery);
         TGListTreeItem *item2 = fViewer->GetSessionHierarchy()->AddItem(item,
                  newquery->fQueryName, fViewer->GetQueryConPict(),
                  fViewer->GetQueryConPict());
         item2->SetUserData(newquery);
         if (query->GetInputList())
            fViewer->GetSessionHierarchy()->AddItem(item2, "InputList");
         if (query->GetOutputList())
            fViewer->GetSessionHierarchy()->AddItem(item2, "OutputList");
      }
   }
   // at the end, update list tree
   fViewer->GetSessionHierarchy()->ClearViewPort();
   fClient->NeedRedraw(fViewer->GetSessionHierarchy());
}

////////////////////////////////////////////////////////////////////////////////
/// Command line handling.

void TSessionFrame::OnCommandLine()
{
   // get command string
   const char *cmd = fCommandTxt->GetText();
   char opt[2];
   // form temporary file path
   TString pathtmp = TString::Format("%s/%s", gSystem->TempDirectory(),
                                     kSession_RedirectCmd);
   // if check box "clear view" is checked, open temp file in write mode
   // (overwrite), in append mode otherwise.
   if (fClearCheck->IsOn())
      snprintf(opt, 2, "w");
   else
      snprintf(opt, 2, "a");

   // if valid Proof session, pass the command to Proof
   if (fViewer->GetActDesc()->fProof &&
       fViewer->GetActDesc()->fProof->IsValid()) {
      // redirect stdout/stderr to temp file
      if (gSystem->RedirectOutput(pathtmp.Data(), opt) != 0) {
         Error("ShowStatus", "stdout/stderr redirection failed; skipping");
         return;
      }
      // execute command line
      fViewer->GetActDesc()->fProof->Exec(cmd);
      // restore back stdout/stderr
      if (gSystem->RedirectOutput(0) != 0) {
         Error("ShowStatus", "stdout/stderr retore failed; skipping");
         return;
      }
      // if check box "clear view" is checked, clear text view
      if (fClearCheck->IsOn())
         fInfoTextView->Clear();
      // load (display) temp file in text view
      fInfoTextView->LoadFile(pathtmp.Data());
      // set focus to "command line" text entry
      fCommandTxt->SetFocus();
   }
   else {
      // if no Proof session, or Proof session not valid,
      // lets execute command line by TApplication

      // redirect stdout/stderr to temp file
      if (gSystem->RedirectOutput(pathtmp.Data(), opt) != 0) {
         Error("ShowStatus", "stdout/stderr redirection failed; skipping");
      }
      // execute command line
      gApplication->ProcessLine(cmd);
      // restore back stdout/stderr
      if (gSystem->RedirectOutput(0) != 0) {
         Error("ShowStatus", "stdout/stderr retore failed; skipping");
      }
      // if check box "clear view" is checked, clear text view
      if (fClearCheck->IsOn())
         fInfoTextView->Clear();
      // load (display) temp file in text view
      fInfoTextView->LoadFile(pathtmp.Data());
      // set focus to "command line" text entry
      fCommandTxt->SetFocus();
   }
   // display bottom of text view
   fInfoTextView->ShowBottom();
}

////////////////////////////////////////////////////////////////////////////////
/// Switch widgets status/visibility for local/remote sessions.

void TSessionFrame::SetLocal(Bool_t local)
{
   if (local) {
      fBtnGetQueries->SetState(kButtonDisabled);
      fBtnShowLog->SetState(kButtonDisabled);
      fTab->HideFrame(fTab->GetTabTab("Options"));
      fTab->HideFrame(fTab->GetTabTab("Packages"));
      fTab->HideFrame(fTab->GetTabTab("DataSets"));
   }
   else {
      fBtnGetQueries->SetState(kButtonUp);
      fBtnShowLog->SetState(kButtonUp);
      fTab->ShowFrame(fTab->GetTabTab("Options"));
      fTab->ShowFrame(fTab->GetTabTab("Packages"));
      fTab->ShowFrame(fTab->GetTabTab("DataSets"));
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Shutdown current session.

void TSessionFrame::ShutdownSession()
{
   // do nothing if connection in progress
   if (fViewer->IsBusy())
      return;

   if (fViewer->GetActDesc()->fLocal) {
      Int_t retval;
      new TGMsgBox(fClient->GetRoot(), this, "Error Shutting down Session",
                   "Shutting down Local Sessions is not allowed !",
                    kMBIconExclamation,kMBOk,&retval);
      return;
   }
   if (!fViewer->GetActDesc()->fAttached ||
       !fViewer->GetActDesc()->fProof ||
       !fViewer->GetActDesc()->fProof->IsValid())
      return;
   // ask for confirmation
   TString m;
   m.Form("Are you sure to shutdown the session \"%s\"",
          fViewer->GetActDesc()->fName.Data());
   Int_t result;
   new TGMsgBox(fClient->GetRoot(), this, "", m.Data(), 0,
                kMBOk | kMBCancel, &result);
   // if confirmed, delete it
   if (result != kMBOk)
      return;
   // remove the Proof session from gROOT list of Proofs
   fViewer->GetActDesc()->fProof->Detach("S");
   // reset connected flag
   fViewer->GetActDesc()->fAttached = kFALSE;
   fViewer->GetActDesc()->fProof = 0;
   // disable animation timer
   fViewer->DisableTimer();
   // change list tree item picture to disconnected pixmap
   TGListTreeItem *item = fViewer->GetSessionHierarchy()->FindChildByData(
                          fViewer->GetSessionItem(), fViewer->GetActDesc());
   if (item) {
      item->SetPictures(fViewer->GetProofDisconPict(),
                        fViewer->GetProofDisconPict());
   }
   // update viewer
   fViewer->OnListTreeClicked(fViewer->GetSessionHierarchy()->GetSelected(),
                              1, 0, 0);
   fViewer->GetSessionHierarchy()->ClearViewPort();
   fClient->NeedRedraw(fViewer->GetSessionHierarchy());
   fViewer->GetStatusBar()->SetText("", 1);
}

//////////////////////////////////////////////////////////////////////////
// Edit Query Frame
//////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// Create a new Query dialog, used by the Session Viewer, to Edit a Query if
/// the editmode flag is set, or to create a new one if not set.

TEditQueryFrame::TEditQueryFrame(TGWindow* p, Int_t w, Int_t h) :
   TGCompositeFrame(p, w, h, kVerticalFrame), fFrmMore(0), fBtnMore(0),
   fTxtQueryName(0), fTxtChain(0), fTxtSelector(0), fTxtOptions(0),
   fNumEntries(0), fNumFirstEntry(0), fTxtParFile(0), fTxtEventList(0),
   fViewer(0), fQuery(0), fChain(0)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Delete query dialog.

TEditQueryFrame::~TEditQueryFrame()
{
   Cleanup();
}

////////////////////////////////////////////////////////////////////////////////
/// Build the "new query" dialog.

void TEditQueryFrame::Build(TSessionViewer *gui)
{
   TGButton *btnTmp;
   fViewer = gui;
   SetCleanup(kDeepCleanup);
   SetLayoutManager(new TGTableLayout(this, 6, 5));

   // add "Query Name" label and text entry
   AddFrame(new TGLabel(this, "Query Name :"),
         new TGTableLayoutHints(0, 1, 0, 1, kLHintsCenterY, 5, 5, 4, 0));
   AddFrame(fTxtQueryName = new TGTextEntry(this,
         (const char *)0, 1), new TGTableLayoutHints(1, 2, 0, 1,
         kLHintsCenterY, 5, 5, 4, 0));

   // add "TChain" label and text entry
   AddFrame(new TGLabel(this, "TChain :"),
         new TGTableLayoutHints(0, 1, 1, 2, kLHintsCenterY, 5, 5, 4, 0));
   AddFrame(fTxtChain = new TGTextEntry(this,
         (const char *)0, 2), new TGTableLayoutHints(1, 2, 1, 2,
         kLHintsCenterY, 5, 5, 4, 0));
   fTxtChain->SetToolTipText("Specify TChain or TDSet from memory or file");
   fTxtChain->SetEnabled(kFALSE);
   // add "Browse" button
   AddFrame(btnTmp = new TGTextButton(this, "Browse..."),
         new TGTableLayoutHints(2, 3, 1, 2, kLHintsCenterY, 5, 0, 4, 8));
   btnTmp->Connect("Clicked()", "TEditQueryFrame", this, "OnBrowseChain()");

   // add "Selector" label and text entry
   AddFrame(new TGLabel(this, "Selector :"),
         new TGTableLayoutHints(0, 1, 2, 3, kLHintsCenterY, 5, 5, 0, 0));
   AddFrame(fTxtSelector = new TGTextEntry(this,
         (const char *)0, 3), new TGTableLayoutHints(1, 2, 2, 3,
         kLHintsCenterY, 5, 5, 0, 0));
   // add "Browse" button
   AddFrame(btnTmp = new TGTextButton(this, "Browse..."),
         new TGTableLayoutHints(2, 3, 2, 3, kLHintsCenterY, 5, 0, 0, 8));
   btnTmp->Connect("Clicked()", "TEditQueryFrame", this, "OnBrowseSelector()");

   // add "Less <<" ("More >>") button
   AddFrame(fBtnMore = new TGTextButton(this, " Less << "),
         new TGTableLayoutHints(2, 3, 4, 5, kLHintsCenterY, 5, 5, 4, 0));
   fBtnMore->Connect("Clicked()", "TEditQueryFrame", this, "OnNewQueryMore()");

   // add (initially hidden) options frame
   fFrmMore = new TGCompositeFrame(this, 200, 200);
   fFrmMore->SetCleanup(kDeepCleanup);

   AddFrame(fFrmMore, new TGTableLayoutHints(0, 3, 5, 6,
         kLHintsExpandX | kLHintsExpandY));
   fFrmMore->SetLayoutManager(new TGTableLayout(fFrmMore, 4, 3));

   // add "Options" label and text entry
   fFrmMore->AddFrame(new TGLabel(fFrmMore, "Options :"),
         new TGTableLayoutHints(0, 1, 0, 1, kLHintsCenterY, 5, 5, 0, 0));
   fFrmMore->AddFrame(fTxtOptions = new TGTextEntry(fFrmMore,
         (const char *)0, 4), new TGTableLayoutHints(1, 2, 0, 1, 0, 17,
         0, 0, 8));
   //fTxtOptions->SetText("ASYN");
   fTxtOptions->SetText("");

   // add "Nb Entries" label and number entry
   fFrmMore->AddFrame(new TGLabel(fFrmMore, "Nb Entries :"),
         new TGTableLayoutHints(0, 1, 1, 2, kLHintsCenterY, 5, 5, 0, 0));
   fFrmMore->AddFrame(fNumEntries = new TGNumberEntry(fFrmMore, 0, 5, -1,
         TGNumberFormat::kNESInteger, TGNumberFormat::kNEAAnyNumber,
         TGNumberFormat::kNELNoLimits), new TGTableLayoutHints(1, 2, 1, 2,
         0, 17, 0, 0, 8));
   // coverity[negative_returns]: no problem with -1, the format is kNESInteger
   fNumEntries->SetIntNumber(-1);
   // add "First Entry" label and number entry
   fFrmMore->AddFrame(new TGLabel(fFrmMore, "First entry :"),
         new TGTableLayoutHints(0, 1, 2, 3, kLHintsCenterY, 5, 5, 0, 0));
   fFrmMore->AddFrame(fNumFirstEntry = new TGNumberEntry(fFrmMore, 0, 5, -1,
         TGNumberFormat::kNESInteger, TGNumberFormat::kNEANonNegative,
         TGNumberFormat::kNELNoLimits), new TGTableLayoutHints(1, 2, 2, 3, 0,
         17, 0, 0, 8));

   // add "Event list" label and text entry
   fFrmMore->AddFrame(new TGLabel(fFrmMore, "Event list :"),
         new TGTableLayoutHints(0, 1, 3, 4, kLHintsCenterY, 5, 5, 0, 0));
   fFrmMore->AddFrame(fTxtEventList = new TGTextEntry(fFrmMore,
         (const char *)0, 6), new TGTableLayoutHints(1, 2, 3, 4, 0, 17,
         5, 0, 0));
   // add "Browse" button
   fFrmMore->AddFrame(btnTmp = new TGTextButton(fFrmMore, "Browse..."),
         new TGTableLayoutHints(2, 3, 3, 4, 0, 6, 0, 0, 8));
   btnTmp->Connect("Clicked()", "TEditQueryFrame", this, "OnBrowseEventList()");

   fTxtQueryName->Associate(this);
   fTxtChain->Associate(this);
   fTxtSelector->Associate(this);
   fTxtOptions->Associate(this);
   fNumEntries->Associate(this);
   fNumFirstEntry->Associate(this);
   fTxtEventList->Associate(this);

   fTxtQueryName->Connect("TextChanged(char*)", "TEditQueryFrame", this,
                        "SettingsChanged()");
   fTxtChain->Connect("TextChanged(char*)", "TEditQueryFrame", this,
                        "SettingsChanged()");
   fTxtSelector->Connect("TextChanged(char*)", "TEditQueryFrame", this,
                        "SettingsChanged()");
   fTxtOptions->Connect("TextChanged(char*)", "TEditQueryFrame", this,
                        "SettingsChanged()");
   fNumEntries->Connect("ValueChanged(Long_t)", "TEditQueryFrame", this,
                        "SettingsChanged()");
   fNumFirstEntry->Connect("ValueChanged(Long_t)", "TEditQueryFrame", this,
                        "SettingsChanged()");
   fTxtEventList->Connect("TextChanged(char*)", "TEditQueryFrame", this,
                        "SettingsChanged()");
}

////////////////////////////////////////////////////////////////////////////////
/// Show/hide options frame and update button text accordingly.

void TEditQueryFrame::OnNewQueryMore()
{
   if (IsVisible(fFrmMore)) {
      HideFrame(fFrmMore);
      fBtnMore->SetText(" More >> ");
   }
   else {
      ShowFrame(fFrmMore);
      fBtnMore->SetText(" Less << ");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Call new chain dialog.

void TEditQueryFrame::OnBrowseChain()
{
   TNewChainDlg *dlg = new TNewChainDlg(fClient->GetRoot(), this);
   dlg->Connect("OnElementSelected(TObject *)", "TEditQueryFrame",
         this, "OnElementSelected(TObject *)");
}

////////////////////////////////////////////////////////////////////////////////
/// Handle OnElementSelected signal coming from new chain dialog.

void TEditQueryFrame::OnElementSelected(TObject *obj)
{
   if (obj) {
      fChain = obj;
      if (obj->IsA() == TChain::Class())
         fTxtChain->SetText(((TChain *)fChain)->GetName());
      else if (obj->IsA() == TDSet::Class())
         fTxtChain->SetText(((TDSet *)fChain)->GetObjName());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Open file browser to choose selector macro.

void TEditQueryFrame::OnBrowseSelector()
{
   TGFileInfo fi;
   fi.fFileTypes = macrotypes;
   new TGFileDialog(fClient->GetRoot(), this, kFDOpen, &fi);
   if (!fi.fFilename) return;
   fTxtSelector->SetText(gSystem->UnixPathName(fi.fFilename));
}

////////////////////////////////////////////////////////////////////////////////
///Browse event list

void TEditQueryFrame::OnBrowseEventList()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Save current settings in main session viewer.

void TEditQueryFrame::OnBtnSave()
{
   // if we are in edition mode and query description is valid,
   // use it, otherwise create a new one
   TQueryDescription *newquery;
   if (fQuery)
      newquery = fQuery;
   else
      newquery = new TQueryDescription();

   // update query description fields
   newquery->fSelectorString  = fTxtSelector->GetText();
   if (fChain) {
      newquery->fTDSetString  = fChain->GetName();
      newquery->fChain        = fChain;
   }
   else {
      newquery->fTDSetString = "";
      newquery->fChain       = 0;
   }
   newquery->fQueryName      = fTxtQueryName->GetText();
   newquery->fOptions        = fTxtOptions->GetText();
   newquery->fNoEntries      = fNumEntries->GetIntNumber();
   newquery->fFirstEntry     = fNumFirstEntry->GetIntNumber();
   newquery->fNbFiles        = 0;
   newquery->fResult         = 0;

   if (newquery->fChain) {
      if (newquery->fChain->IsA() == TChain::Class())
         newquery->fNbFiles = ((TChain *)newquery->fChain)->GetListOfFiles()->GetEntriesFast();
      else if (newquery->fChain->IsA() == TDSet::Class())
         newquery->fNbFiles = ((TDSet *)newquery->fChain)->GetListOfElements()->GetSize();
   }
   // update user data with modified query description
   TGListTreeItem *item = fViewer->GetSessionHierarchy()->GetSelected();
   fViewer->GetSessionHierarchy()->RenameItem(item, newquery->fQueryName);
   item->SetUserData(newquery);
   // update list tree
   fViewer->GetSessionHierarchy()->ClearViewPort();
   fClient->NeedRedraw(fViewer->GetSessionHierarchy());
   fTxtQueryName->SelectAll();
   fTxtQueryName->SetFocus();
   fViewer->WriteConfiguration();
   fViewer->GetQueryFrame()->Modified(kFALSE);
   if (fViewer->GetActDesc()->fLocal ||
      (fViewer->GetActDesc()->fConnected &&
       fViewer->GetActDesc()->fAttached &&
       fViewer->GetActDesc()->fProof &&
       fViewer->GetActDesc()->fProof->IsValid())) {
      fViewer->GetQueryFrame()->GetTab()->SetTab("Status");
      fViewer->GetQueryFrame()->OnBtnSubmit();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Settings have changed, update GUI accordingly.

void TEditQueryFrame::SettingsChanged()
{
   if (fQuery) {
      if ((strcmp(fQuery->fSelectorString.Data(), fTxtSelector->GetText())) ||
          (strcmp(fQuery->fQueryName.Data(), fTxtQueryName->GetText())) ||
          (strcmp(fQuery->fOptions.Data(), fTxtOptions->GetText())) ||
          (fQuery->fNoEntries  != fNumEntries->GetIntNumber()) ||
          (fQuery->fFirstEntry != fNumFirstEntry->GetIntNumber()) ||
          (fQuery->fChain != fChain)) {
         fViewer->GetQueryFrame()->Modified(kTRUE);
      }
      else {
         fViewer->GetQueryFrame()->Modified(kFALSE);
      }
   }
   else {
      if ((fTxtQueryName->GetText()) &&
         ((fTxtQueryName->GetText()) ||
          (fTxtChain->GetText())))
         fViewer->GetQueryFrame()->Modified(kTRUE);
      else
         fViewer->GetQueryFrame()->Modified(kFALSE);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Update entry fields with query description values.

void TEditQueryFrame::UpdateFields(TQueryDescription *desc)
{
   fChain = 0;
   fQuery = desc;
   fTxtChain->SetText("");
   if (desc->fChain) {
      fChain = desc->fChain;
      fTxtChain->SetText(desc->fTDSetString);
   }
   fTxtQueryName->SetText(desc->fQueryName);
   fTxtSelector->SetText(desc->fSelectorString);
   fTxtOptions->SetText(desc->fOptions);
   fNumEntries->SetIntNumber(desc->fNoEntries);
   fNumFirstEntry->SetIntNumber(desc->fFirstEntry);
   fTxtEventList->SetText(desc->fEventList);
}

////////////////////////////////////////////////////////////////////////////////
// Query Frame

////////////////////////////////////////////////////////////////////////////////
/// Constructor

TSessionQueryFrame::TSessionQueryFrame(TGWindow* p, Int_t w, Int_t h) :
   TGCompositeFrame(p, w, h), fBtnSubmit(0), fBtnFinalize(0), fBtnStop(0),
   fBtnAbort(0), fBtnShowLog(0), fBtnRetrieve(0), fBtnSave(0), fInfoTextView(0),
   fModified(0), fFiles(0), fFirst(0), fEntries(0), fPrevTotal(0),
   fPrevProcessed(0), fLabInfos(0), fLabStatus(0), fTotal(0), fRate(0),
   fStatus(kStopped), fTab(0), fFA(0), fFB(0), fFC(0), fFD(0), frmProg(0),
   fECanvas(0), fStatsCanvas(0), fViewer(0), fDesc(0)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TSessionQueryFrame::~TSessionQueryFrame()
{
   Cleanup();
}

////////////////////////////////////////////////////////////////////////////////
/// Build query information frame.

void TSessionQueryFrame::Build(TSessionViewer *gui)
{
   SetLayoutManager(new TGVerticalLayout(this));
   SetCleanup(kDeepCleanup);
   fFirst = fEntries = fPrevTotal = 0;
   fPrevProcessed = 0;
   fViewer = gui;
   fModified = kFALSE;

   // main query tab
   fTab = new TGTab(this, 200, 200);
   AddFrame(fTab, new TGLayoutHints(kLHintsTop | kLHintsExpandX |
         kLHintsExpandY, 2, 2, 2, 2));

   // add "Status" tab element
   TGCompositeFrame *tf = fTab->AddTab("Status");
   fFB = new TGCompositeFrame(tf, 100, 100, kVerticalFrame);
   tf->AddFrame(fFB, new TGLayoutHints(kLHintsTop | kLHintsLeft |
         kLHintsExpandX | kLHintsExpandY));

   // new frame containing control buttons and feedback histos canvas
   TGCompositeFrame* frmcanvas = new TGHorizontalFrame(fFB, 350, 100);
   // control buttons frame
   TGCompositeFrame* frmBut2 = new TGVerticalFrame(frmcanvas, 150, 100);
   fBtnSubmit = new TGTextButton(frmBut2, "        Submit        ");
   fBtnSubmit->SetToolTipText("Submit (process) selected query");
   frmBut2->AddFrame(fBtnSubmit,new TGLayoutHints(kLHintsCenterY | kLHintsLeft |
         kLHintsExpandX, 5, 5, 5, 5));
   fBtnStop = new TGTextButton(frmBut2, "Stop");
   fBtnStop->SetToolTipText("Stop processing query");
   frmBut2->AddFrame(fBtnStop,new TGLayoutHints(kLHintsCenterY | kLHintsLeft |
         kLHintsExpandX, 5, 5, 5, 5));
   fBtnAbort = new TGTextButton(frmBut2, "Abort");
   fBtnAbort->SetToolTipText("Abort processing query");
   frmBut2->AddFrame(fBtnAbort,new TGLayoutHints(kLHintsCenterY | kLHintsLeft |
         kLHintsExpandX, 5, 5, 5, 5));
   frmcanvas->AddFrame(frmBut2, new TGLayoutHints(kLHintsLeft | kLHintsCenterY |
         kLHintsExpandY));
   // feedback histos embedded canvas
   fECanvas = new TRootEmbeddedCanvas("fECanvas", frmcanvas, 400, 150);
   fStatsCanvas = fECanvas->GetCanvas();
   fStatsCanvas->SetFillColor(10);
   fStatsCanvas->SetBorderMode(0);
   frmcanvas->AddFrame(fECanvas, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY,
            4, 4, 4, 4));
   fFB->AddFrame(frmcanvas, new TGLayoutHints(kLHintsLeft | kLHintsTop |
         kLHintsExpandX | kLHintsExpandY));

   // progress infos label
   fLabInfos = new TGLabel(fFB, "                                  ");
   fFB->AddFrame(fLabInfos, new TGLayoutHints(kLHintsLeft, 5, 5, 5, 5));
   // progress status label
   fLabStatus = new TGLabel(fFB, "                                  ");
   fFB->AddFrame(fLabStatus, new TGLayoutHints(kLHintsLeft, 5, 5, 5, 5));

   //progress bar
   frmProg = new TGHProgressBar(fFB, TGProgressBar::kFancy, 350 - 20);
   frmProg->ShowPosition();
   frmProg->SetBarColor("green");
   fFB->AddFrame(frmProg, new TGLayoutHints(kLHintsExpandX, 5, 5, 5, 5));
   // total progress infos
   fFB->AddFrame(fTotal = new TGLabel(fFB,
         " Estimated time left : 0 sec (--- events of --- processed) "),
         new TGLayoutHints(kLHintsLeft, 5, 5, 5, 5));
   // progress rate infos
   fFB->AddFrame(fRate = new TGLabel(fFB,
         " Processing Rate : -- events/sec    "),
         new TGLayoutHints(kLHintsLeft, 5, 5, 5, 5));

   // add "Results" tab element
   tf = fTab->AddTab("Results");
   fFC = new TGCompositeFrame(tf, 100, 100, kVerticalFrame);
   tf->AddFrame(fFC, new TGLayoutHints(kLHintsTop | kLHintsLeft |
         kLHintsExpandX | kLHintsExpandY));
   // query result (header) information text view
   fInfoTextView = new TGTextView(fFC, 330, 185, "", kSunkenFrame |
         kDoubleBorder);
   fFC->AddFrame(fInfoTextView, new TGLayoutHints(kLHintsTop | kLHintsLeft |
         kLHintsExpandY | kLHintsExpandX, 5, 5, 10, 10));

   // add "Retrieve", "Finalize" and "Show Log" buttons
   TGCompositeFrame* frmBut3 = new TGHorizontalFrame(fFC, 350, 100);
   fBtnRetrieve = new TGTextButton(frmBut3, "Retrieve");
   fBtnRetrieve->SetToolTipText("Retrieve query results");
   frmBut3->AddFrame(fBtnRetrieve,new TGLayoutHints(kLHintsTop | kLHintsLeft |
         kLHintsExpandX, 5, 5, 10, 10));
   fBtnFinalize = new TGTextButton(frmBut3, "Finalize");
   fBtnFinalize->SetToolTipText("Finalize query");
   frmBut3->AddFrame(fBtnFinalize,new TGLayoutHints(kLHintsTop | kLHintsLeft |
         kLHintsExpandX, 5, 5, 10, 10));
   fBtnShowLog = new TGTextButton(frmBut3, "Show Log");
   fBtnShowLog->SetToolTipText("Show query log (open log window)");
   frmBut3->AddFrame(fBtnShowLog,new TGLayoutHints(kLHintsTop | kLHintsLeft |
         kLHintsExpandX, 5, 5, 10, 10));
   fFC->AddFrame(frmBut3, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX));

   // add "Results" tab element
   tf = fTab->AddTab("Edit Query");
   fFD = new TEditQueryFrame(tf, 100, 100);
   fFD->Build(fViewer);
   tf->AddFrame(fFD, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 10, 0));
   TString btntxt;
   if (fViewer->GetActDesc()->fLocal ||
      (fViewer->GetActDesc()->fProof &&
       fViewer->GetActDesc()->fProof->IsValid())) {
      btntxt = "         Submit         ";
   }
   else {
      btntxt = "     Apply changes      ";
   }
   tf->AddFrame(fBtnSave = new TGTextButton(tf, btntxt),
                new TGLayoutHints(kLHintsTop | kLHintsLeft, 10, 5, 25, 5));

   // connect button actions to functions
   fBtnSave->Connect("Clicked()", "TEditQueryFrame", fFD,
         "OnBtnSave()");
   fBtnSubmit->Connect("Clicked()", "TSessionQueryFrame", this,
         "OnBtnSubmit()");
   fBtnFinalize->Connect("Clicked()", "TSessionQueryFrame", this,
         "OnBtnFinalize()");
   fBtnStop->Connect("Clicked()", "TSessionQueryFrame", this,
         "OnBtnStop()");
   fBtnAbort->Connect("Clicked()", "TSessionQueryFrame", this,
         "OnBtnAbort()");
   fBtnShowLog->Connect("Clicked()", "TSessionQueryFrame", this,
         "OnBtnShowLog()");
   fBtnRetrieve->Connect("Clicked()", "TSessionQueryFrame", this,
         "OnBtnRetrieve()");
//   fBtnSave->SetState(kButtonDisabled);
   Resize(350, 310);
}

////////////////////////////////////////////////////////////////////////////////
/// Notify changes in query editor settings.

void TSessionQueryFrame::Modified(Bool_t mod)
{
   fModified = mod;
   if (fModified) {
      fBtnSave->SetState(kButtonUp);
   }
   else {
      fBtnSave->SetState(kButtonDisabled);
   }
   if (fViewer->GetActDesc()->fLocal ||
      (fViewer->GetActDesc()->fProof &&
       fViewer->GetActDesc()->fProof->IsValid()))
      fBtnSave->SetState(kButtonUp);
}

////////////////////////////////////////////////////////////////////////////////
/// Feedback function connected to Feedback signal.
/// Used to update feedback histograms.

void TSessionQueryFrame::Feedback(TList *objs)
{
   // if no actual session, just return
   if (!fViewer->GetActDesc()->fAttached)
      return;
   if (!fViewer->GetActDesc()->fProof)
      return;
   if ((fViewer->GetActDesc()->fActQuery) &&
       (fViewer->GetActDesc()->fActQuery->fStatus !=
        TQueryDescription::kSessionQuerySubmitted) &&
       (fViewer->GetActDesc()->fActQuery->fStatus !=
        TQueryDescription::kSessionQueryRunning) )
      return;
   TProof *sender = dynamic_cast<TProof*>((TQObject*)gTQSender);
   // if Proof sender match actual session one, update feedback histos
   if (sender && (sender == fViewer->GetActDesc()->fProof))
      UpdateHistos(objs);
}

////////////////////////////////////////////////////////////////////////////////
/// Update feedback histograms.

void TSessionQueryFrame::UpdateHistos(TList *objs)
{
   TVirtualPad *save = gPad;
   TObject *o;
   Int_t pos = 1;
   Int_t i = 0;
   while (kFeedbackHistos[i]) {
      // check if user has selected this histogram in the option menu
      if (fViewer->GetCascadeMenu()->IsEntryChecked(41+i)) {
         if ( (o = objs->FindObject(kFeedbackHistos[i]))) {
            fStatsCanvas->cd(pos);
            gPad->SetEditable(kTRUE);
            if (TH1 *h = dynamic_cast<TH1*>(o)) {
               h->SetStats(0);
               h->SetBarWidth(0.75);
               h->SetBarOffset(0.125);
               h->SetFillColor(9);
               h->Draw("bar");
               pos++;
            }
            else if (TH2 *h2 = dynamic_cast<TH2*>(o)) {
               h2->Draw();
               pos++;
            }
            gPad->Modified();
         }
      }
      i++;
   }
   // update canvas
   fStatsCanvas->cd();
   fStatsCanvas->Modified();
   fStatsCanvas->Update();
   if (save != 0) {
      save->cd();
   } else {
      gPad = 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Update progress bar and status labels.

void TSessionQueryFrame::Progress(Long64_t total, Long64_t processed)
{
   Long_t tt;
   UInt_t hh=0, mm=0, ss=0;
   TString stm;
   // if no actual session, just return
   if (!fViewer->GetActDesc()->fProof)
      return;
   // if Proof sender does't match actual session one, return
   TProof *sender = dynamic_cast<TProof*>((TQObject*)gTQSender);
   if (!sender || (sender != fViewer->GetActDesc()->fProof))
      return;

   if ((fViewer->GetActDesc()->fActQuery) &&
       (fViewer->GetActDesc()->fActQuery->fStatus !=
        TQueryDescription::kSessionQuerySubmitted) &&
       (fViewer->GetActDesc()->fActQuery->fStatus !=
        TQueryDescription::kSessionQueryRunning) ) {
      fTotal->SetText(" Estimated time left : 0 sec (0 events of 0 processed)        ");
      fRate->SetText(" Processing Rate : 0.0f events/sec   ");
      frmProg->Reset();
      fFB->Layout();
      return;
   }

   if (total < 0)
      total = fPrevTotal;
   else
      fPrevTotal = total;

   // if no change since last call, just return
   if (fPrevProcessed == processed)
      return;
   TString buf;

   // Update information at first call
   if (fEntries != total) {
      buf.Form("PROOF cluster : \"%s\" - %d worker nodes",
               fViewer->GetActDesc()->fProof->GetMaster(),
               fViewer->GetActDesc()->fProof->GetParallel());
      fLabInfos->SetText(buf);

      fEntries = total;
      buf.Form(" %d files, %lld events, starting event %lld",
               fFiles, fEntries, fFirst);
      fLabStatus->SetText(buf);
   }

   // compute progress bar position and update
   Float_t pos = (Float_t)((Double_t)(processed * 100)/(Double_t)total);
   frmProg->SetPosition(pos);
   // if 100%, stop animation and set icon to "connected"
   if (pos >= 100.0) {
      fViewer->SetChangePic(kFALSE);
      fViewer->ChangeRightLogo("monitor01.xpm");
   }

   // get current time
   if ((fViewer->GetActDesc()->fActQuery->fStatus ==
        TQueryDescription::kSessionQueryRunning) ||
       (fViewer->GetActDesc()->fActQuery->fStatus ==
        TQueryDescription::kSessionQuerySubmitted))
      fViewer->GetActDesc()->fActQuery->fEndTime = gSystem->Now();
   TTime tdiff = fViewer->GetActDesc()->fActQuery->fEndTime -
                 fViewer->GetActDesc()->fActQuery->fStartTime;
   Float_t eta = 0;
   if (processed)
      eta = ((Float_t)((Long64_t)tdiff)*total/Float_t(processed) -
            Long64_t(tdiff))/1000.;

   tt = (Long_t)eta;
   if (tt > 0) {
      hh = (UInt_t)(tt / 3600);
      mm = (UInt_t)((tt % 3600) / 60);
      ss = (UInt_t)((tt % 3600) % 60);
   }
   if (hh)
      stm.Form("%d h %d min %d sec", hh, mm, ss);
   else if (mm)
      stm.Form("%d min %d sec", mm, ss);
   else
      stm.Form("%d sec", ss);
   if (processed == total) {
      // finished
      tt = (Long_t) Long64_t(tdiff)/1000;
      if (tt > 0) {
         hh = (UInt_t)(tt / 3600);
         mm = (UInt_t)((tt % 3600) / 60);
         ss = (UInt_t)((tt % 3600) % 60);
      }
      if (hh)
         stm.Form("%d h %d min %d sec", hh, mm, ss);
      else if (mm)
         stm.Form("%d min %d sec", mm, ss);
      else
         stm.Form("%d sec", ss);
      buf.Form(" Processed : %lld events in %s", total, stm.Data());
      fTotal->SetText(buf.Data());
   } else {
      // update status infos
      buf.Form(" Estimated time left : %s (%lld events of %lld processed)        ",
               stm.Data(), processed, total);
      fTotal->SetText(buf.Data());
   }
   if (processed > 0 && (Long64_t)tdiff > 0) {
      buf.Form(" Processing Rate : %.1f events/sec   ",
               (Float_t)processed/(Long64_t)tdiff*1000.);
      fRate->SetText(buf);
   }
   fPrevProcessed = processed;

   fFB->Layout();
}

////////////////////////////////////////////////////////////////////////////////
/// New version of Progress (just forward to the old version
/// for the time being).

void TSessionQueryFrame::Progress(Long64_t total, Long64_t processed,
                                  Long64_t /*bytesread*/, Float_t /*initTime*/,
                                  Float_t /*procTime*/, Float_t /*evtrti*/,
                                  Float_t /*mbrti*/, Int_t /*actw*/,
                                  Int_t /*tses*/, Float_t /*eses*/)
{
   Progress(total, processed);
}

////////////////////////////////////////////////////////////////////////////////
/// Update progress bar and status labels.

void TSessionQueryFrame::ProgressLocal(Long64_t total, Long64_t processed)
{
   Long_t tt;
   UInt_t hh=0, mm=0, ss=0;
   TString stm;
   TString cproc;
   Int_t status;

   switch (fViewer->GetActDesc()->fActQuery->fStatus) {

      case TQueryDescription::kSessionQueryAborted:
         cproc = " - ABORTED";
         status = kAborted;
         break;
      case TQueryDescription::kSessionQueryStopped:
         cproc = " - STOPPED";
         status = kStopped;
         break;
      case TQueryDescription::kSessionQueryRunning:
         cproc = " ";
         status = kRunning;
         break;
      case TQueryDescription::kSessionQueryCompleted:
      case TQueryDescription::kSessionQueryFinalized:
         cproc = " ";
         status = kDone;
         break;
      default:
         status = -1;
         break;
   }
   if (processed < 0) processed = 0;

   frmProg->SetBarColor("green");
   if (status == kAborted)
      frmProg->SetBarColor("red");
   else if (status == kStopped)
      frmProg->SetBarColor("yellow");
   else if (status == -1 ) {
      fTotal->SetText(" Estimated time left : 0 sec (0 events of 0 processed)        ");
      fRate->SetText(" Processing Rate : 0.0f events/sec   ");
      frmProg->Reset();
      fFB->Layout();
      return;
   }

   if (total < 0)
      total = fPrevTotal;
   else
      fPrevTotal = total;

   // if no change since last call, just return
   TString buf;

   // Update information at first call
   if (fEntries != total) {
      fLabInfos->SetText("Local Session");

      fEntries = total;
      buf.Form(" %d files, %lld events, starting event %lld",
               fFiles, fEntries, fFirst);
      fLabStatus->SetText(buf.Data());
   }

   // compute progress bar position and update
   Float_t pos = 0.0;
   if (processed > 0 && total > 0)
      pos = (Float_t)((Double_t)(processed * 100)/(Double_t)total);
   frmProg->SetPosition(pos);
   // if 100%, stop animation and set icon to "connected"
   if (pos >= 100.0) {
      fViewer->SetChangePic(kFALSE);
      fViewer->ChangeRightLogo("monitor01.xpm");
   }

   // get current time
   if (status == kRunning)
      fViewer->GetActDesc()->fActQuery->fEndTime = gSystem->Now();
   TTime tdiff = fViewer->GetActDesc()->fActQuery->fEndTime -
                 fViewer->GetActDesc()->fActQuery->fStartTime;
   Float_t eta = 0;
   if (processed)
      eta = ((Float_t)((Long64_t)tdiff)*total/(Float_t)(processed) -
            (Long64_t)(tdiff))/1000.;

   tt = (Long_t)eta;
   if (tt > 0) {
      hh = (UInt_t)(tt / 3600);
      mm = (UInt_t)((tt % 3600) / 60);
      ss = (UInt_t)((tt % 3600) % 60);
   }
   if (hh)
      stm = TString::Format("%d h %d min %d sec", hh, mm, ss);
   else if (mm)
      stm = TString::Format("%d min %d sec", mm, ss);
   else
      stm = TString::Format("%d sec", ss);
   if ((processed != total) && (status == kRunning)) {
      // update status infos
      buf.Form(" Estimated time left : %s (%lld events of %lld processed)        ",
               stm.Data(), processed, total);
      fTotal->SetText(buf);
   } else {
      tt = (Long_t) Long64_t(tdiff)/1000;
      if (tt > 0) {
         hh = (UInt_t)(tt / 3600);
         mm = (UInt_t)((tt % 3600) / 60);
         ss = (UInt_t)((tt % 3600) % 60);
      }
      if (hh)
         stm = TString::Format("%d h %d min %d sec", hh, mm, ss);
      else if (mm)
         stm = TString::Format("%d min %d sec", mm, ss);
      else
         stm = TString::Format("%d sec", ss);
      buf.Form(" Processed : %ld events in %s",
               (Long_t)processed, stm.Data());
      buf += cproc;
      fTotal->SetText(buf.Data());
   }
   if (processed > 0 && (Long64_t)tdiff > 0) {
      buf.Form(" Processing Rate : %.1f events/sec   ",
               (Float_t)processed/(Long64_t)tdiff*1000.);
      fRate->SetText(buf.Data());
   }
   fPrevProcessed = processed;

   fFB->Layout();
}

////////////////////////////////////////////////////////////////////////////////
/// Indicate that Cancel or Stop was clicked.

void TSessionQueryFrame::IndicateStop(Bool_t aborted)
{
   if (aborted == kTRUE) {
      // Aborted
      frmProg->SetBarColor("red");
   }
   else {
      // Stopped
      frmProg->SetBarColor("yellow");
   }
   // disconnect progress related signals
   if (fViewer->GetActDesc()->fProof &&
       fViewer->GetActDesc()->fProof->IsValid()) {
      fViewer->GetActDesc()->fProof->Disconnect("Progress(Long64_t,Long64_t)",
               this, "Progress(Long64_t,Long64_t)");
      fViewer->GetActDesc()->fProof->Disconnect("Progress(Long64_t,Long64_t,Long64_t,Float_t,Float_t,Float_t,Float_t)",
               this, "Progress(Long64_t,Long64_t,Long64_t,Float_t,Float_t,Float_t,Float_t)");
      fViewer->GetActDesc()->fProof->Disconnect("Progress(Long64_t,Long64_t,Long64_t,Float_t,Float_t,Float_t,Float_t,Int_t,Int_t,Float_t)",
               this, "Progress(Long64_t,Long64_t,Long64_t,Float_t,Float_t,Float_t,Float_t,Int_t,Int_t,Float_t)");
      fViewer->GetActDesc()->fProof->Disconnect("StopProcess(Bool_t)", this,
               "IndicateStop(Bool_t)");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Reset progress frame information fields.

void TSessionQueryFrame::ResetProgressDialog(const char * /*selector*/, Int_t files,
                                        Long64_t first, Long64_t entries)
{
   TString buf;
   fFiles         = files > 0 ? files : 0;
   fFirst         = first;
   fEntries       = entries;
   fPrevProcessed = 0;
   fPrevTotal     = 0;

   if (!fViewer->GetActDesc()->fLocal) {
      frmProg->SetBarColor("green");
      frmProg->Reset();
   }

   buf.Form("%0d files, %0lld events, starting event %0lld",
            fFiles > 0 ? fFiles : 0, fEntries > 0 ? fEntries : 0,
            fFirst >= 0 ? fFirst : 0);
   fLabStatus->SetText(buf.Data());
   // Reconnect the slots
   if (fViewer->GetActDesc()->fProof &&
       fViewer->GetActDesc()->fProof->IsValid()) {
      fViewer->GetActDesc()->fProof->Connect("Progress(Long64_t,Long64_t)",
               "TSessionQueryFrame", this, "Progress(Long64_t,Long64_t)");
      fViewer->GetActDesc()->fProof->Connect("Progress(Long64_t,Long64_t,Long64_t,Float_t,Float_t,Float_t,Float_t)",
               "TSessionQueryFrame", this,
               "Progress(Long64_t,Long64_t,Long64_t,Float_t,Float_t,Float_t,Float_t)");
      fViewer->GetActDesc()->fProof->Connect("StopProcess(Bool_t)",
               "TSessionQueryFrame", this, "IndicateStop(Bool_t)");
      buf.Form("PROOF cluster : \"%s\" - %d worker nodes",
                fViewer->GetActDesc()->fProof->GetMaster(),
                fViewer->GetActDesc()->fProof->GetParallel());
      fLabInfos->SetText(buf.Data());
   }
   else if (fViewer->GetActDesc()->fLocal) {
      fStatsCanvas->Clear();
      fLabInfos->SetText("Local Session");
      fLabStatus->SetText(" ");
   }
   else {
      fLabInfos->SetText(" ");
      fLabStatus->SetText(" ");
   }
   fFB->Layout();
}

////////////////////////////////////////////////////////////////////////////////
/// Finalize query.

void TSessionQueryFrame::OnBtnFinalize()
{
   // check if Proof is valid
   if (fViewer->GetActDesc()->fProof &&
       fViewer->GetActDesc()->fProof->IsValid()) {
      gPad->SetEditable(kFALSE);
      TGListTreeItem *item = fViewer->GetSessionHierarchy()->GetSelected();
      if (!item) return;
      TObject *obj = (TObject *)item->GetUserData();
      if ((obj) && (obj->IsA() == TQueryDescription::Class())) {
         // as it can take time, set watch cursor
         gVirtualX->SetCursor(GetId(),gVirtualX->CreateCursor(kWatch));
         TQueryDescription *query = (TQueryDescription *)obj;
         fViewer->GetActDesc()->fProof->Finalize(query->fReference);
         UpdateButtons(query);
         // restore cursor
         gVirtualX->SetCursor(GetId(), 0);
      }
   }
   if (fViewer->GetActDesc()->fLocal) {
      gPad->SetEditable(kFALSE);
      TChain *chain = (TChain *)fViewer->GetActDesc()->fActQuery->fChain;
      if (chain && chain->GetPlayer())
         chain->GetPlayer()->GetSelectorFromFile()->Terminate();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Stop processing query.

void TSessionQueryFrame::OnBtnStop()
{
   // check for proof validity
   if (fViewer->GetActDesc()->fProof &&
       fViewer->GetActDesc()->fProof->IsValid()) {
      fViewer->GetActDesc()->fProof->StopProcess(kFALSE);
   }
   if (fViewer->GetActDesc()->fLocal) {
      gROOT->SetInterrupt();
      fViewer->GetActDesc()->fActQuery->fStatus =
         TQueryDescription::kSessionQueryStopped;
   }
   // stop icon animation and set connected icon
   fViewer->ChangeRightLogo("monitor01.xpm");
   fViewer->SetChangePic(kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Show query log.

void TSessionQueryFrame::OnBtnShowLog()
{
   TGListTreeItem *item = fViewer->GetSessionHierarchy()->GetSelected();
   if (!item) return;
   TObject *obj = (TObject *)item->GetUserData();
   if ((!obj) || (obj->IsA() != TQueryDescription::Class()))
      return;
   TQueryDescription *query = (TQueryDescription *)obj;
   fViewer->ShowLog(query->fReference.Data());
}

////////////////////////////////////////////////////////////////////////////////
/// Retrieve query.

void TSessionQueryFrame::OnBtnRetrieve()
{
   // check for proof validity
   if (fViewer->GetActDesc()->fAttached &&
       fViewer->GetActDesc()->fProof &&
       fViewer->GetActDesc()->fProof->IsValid()) {
      TGListTreeItem *item = fViewer->GetSessionHierarchy()->GetSelected();
      if (!item) return;
      TObject *obj = (TObject *)item->GetUserData();
      if (obj && obj->IsA() == TQueryDescription::Class()) {
         // as it can take time, set watch cursor
         gVirtualX->SetCursor(GetId(), gVirtualX->CreateCursor(kWatch));
         TQueryDescription *query = (TQueryDescription *)obj;
         Int_t rc = fViewer->GetActDesc()->fProof->Retrieve(query->fReference);
         if (rc == 0)
            fViewer->OnCascadeMenu();
         // restore cursor
         gVirtualX->SetCursor(GetId(), 0);
      }
   }
   if (fViewer->GetActDesc()->fLocal) {
      TGListTreeItem *item=0, *item2=0;
      item = fViewer->GetSessionHierarchy()->FindItemByObj(fViewer->GetSessionItem(),
                                                           fViewer->GetActDesc());
      if (item) {
         item2 = fViewer->GetSessionHierarchy()->FindItemByObj(item,
                                    fViewer->GetActDesc()->fActQuery);
      }
      if (item2) {
         // add input and output list entries
         TChain *chain = (TChain *)fViewer->GetActDesc()->fActQuery->fChain;
         if (chain && chain->GetPlayer()) {
            TSelector *selector = chain->GetPlayer()->GetSelectorFromFile();
            if (selector) {
               TList *objlist = selector->GetOutputList();
               if (objlist)
                  if (!fViewer->GetSessionHierarchy()->FindChildByName(item2, "OutputList"))
                     fViewer->GetSessionHierarchy()->AddItem(item2, "OutputList");
            }
         }
      }
      // update list tree, query frame information, and buttons state
      fViewer->GetSessionHierarchy()->ClearViewPort();
      fClient->NeedRedraw(fViewer->GetSessionHierarchy());
      UpdateInfos();
      UpdateButtons(fViewer->GetActDesc()->fActQuery);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Abort processing query.

void TSessionQueryFrame::OnBtnAbort()
{
   // check for proof validity
   if (fViewer->GetActDesc()->fProof &&
       fViewer->GetActDesc()->fProof->IsValid()) {
      fViewer->GetActDesc()->fProof->StopProcess(kTRUE);
   }
   if (fViewer->GetActDesc()->fLocal) {
      gROOT->SetInterrupt();
      fViewer->GetActDesc()->fActQuery->fStatus =
         TQueryDescription::kSessionQueryAborted;
   }
   // stop icon animation and set connected icon
   fViewer->ChangeRightLogo("monitor01.xpm");
   fViewer->SetChangePic(kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Submit query.

void TSessionQueryFrame::OnBtnSubmit()
{
   Int_t retval;
   Long64_t id = 0;
   TGListTreeItem *item = fViewer->GetSessionHierarchy()->GetSelected();
   if (!item) return;
   // retrieve query description attached to list tree item
   TObject *obj = (TObject *)item->GetUserData();
   if (!obj || obj->IsA() != TQueryDescription::Class())
      return;
   TQueryDescription *newquery = (TQueryDescription *)obj;
   // reset progress information
   ResetProgressDialog(newquery->fSelectorString,
         newquery->fNbFiles, newquery->fFirstEntry, newquery->fNoEntries);
   // set query start time
   newquery->fStartTime = gSystem->Now();
   fViewer->GetActDesc()->fNbHistos = 0;
   // check for proof validity
   if (fViewer->GetActDesc()->fProof &&
       fViewer->GetActDesc()->fProof->IsValid()) {
      fViewer->GetActDesc()->fProof->SetBit(TProof::kUsingSessionGui);
      // set query description status to submitted
      newquery->fStatus = TQueryDescription::kSessionQuerySubmitted;
      // if feedback option selected
      if (fViewer->GetOptionsMenu()->IsEntryChecked(kOptionsFeedback)) {
         Int_t i = 0;
         // browse list of feedback histos and check user's selected ones
         while (kFeedbackHistos[i]) {
            if (fViewer->GetCascadeMenu()->IsEntryChecked(41+i)) {
               fViewer->GetActDesc()->fProof->AddFeedback(kFeedbackHistos[i]);
               fViewer->GetActDesc()->fNbHistos++;
            }
            i++;
         }
         // connect feedback signal
         fViewer->GetActDesc()->fProof->Connect("Feedback(TList *objs)",
                           "TSessionQueryFrame", fViewer->GetQueryFrame(),
                           "Feedback(TList *objs)");
         gROOT->Time();
      }
      else {
         // if feedback option not selected, clear Proof's feedback option
         fViewer->GetActDesc()->fProof->ClearFeedback();
      }
      // set current proof session
      fViewer->GetActDesc()->fProof->cd();
      // check if parameter file has been specified
      if (newquery->fChain) {
         if (fViewer->GetActDesc()->fProof->IsLite()) {
            newquery->fOptions = "";
         }
         // set query reference id
         newquery->fReference= TString::Format("session-%s:q%d",
                            fViewer->GetActDesc()->fProof->GetSessionTag(),
                            fViewer->GetActDesc()->fProof->GetSeqNum()+1);
         if (newquery->fChain->IsA() == TChain::Class()) {
            // TChain case
            newquery->fStatus = TQueryDescription::kSessionQuerySubmitted;
            ((TChain *)newquery->fChain)->SetProof(fViewer->GetActDesc()->fProof);
            id = ((TChain *)newquery->fChain)->Process(newquery->fSelectorString,
                    newquery->fOptions,
                    newquery->fNoEntries > 0 ? newquery->fNoEntries : 1234567890,
                    newquery->fFirstEntry);
         }
         else if (newquery->fChain->IsA() == TDSet::Class()) {
            // TDSet case
            newquery->fStatus = TQueryDescription::kSessionQuerySubmitted;
            id = ((TDSet *)newquery->fChain)->Process(newquery->fSelectorString,
                  newquery->fOptions,
                  newquery->fNoEntries,
                  newquery->fFirstEntry);
         }
      }
      else {
         Error("Submit", "No TChain defined; skipping");
         newquery->fStatus = TQueryDescription::kSessionQueryCreated;
         return;
      }
      // set query reference id to unique identifier
      newquery->fReference= TString::Format("session-%s:q%lld",
                            fViewer->GetActDesc()->fProof->GetSessionTag(), id);
      // start icon animation
      fViewer->SetChangePic(kTRUE);
   }
   else if (fViewer->GetActDesc()->fLocal) { // local session case
      // if feedback option selected
      if (fViewer->GetOptionsMenu()->IsEntryChecked(kOptionsFeedback)) {
         Int_t i = 0;
         // browse list of feedback histos and check user's selected ones
         while (kFeedbackHistos[i]) {
            if (fViewer->GetCascadeMenu()->IsEntryChecked(41+i)) {
               fViewer->GetActDesc()->fNbHistos++;
            }
            i++;
         }
      }
      if (newquery->fChain) {
         if (newquery->fChain->IsA() == TChain::Class()) {
            // TChain case
            newquery->fStatus = TQueryDescription::kSessionQueryRunning;
            fViewer->EnableTimer();
            UpdateButtons(newquery);
            gPad->SetEditable(kFALSE);
            ((TChain *)newquery->fChain)->SetTimerInterval(100);
            id = ((TChain *)newquery->fChain)->Process(newquery->fSelectorString,
                  newquery->fOptions,
                  newquery->fNoEntries > 0 ? newquery->fNoEntries : 1234567890,
                  newquery->fFirstEntry);
            ((TChain *)newquery->fChain)->SetTimerInterval(0);
            OnBtnRetrieve();
            TChain *chain = (TChain *)newquery->fChain;
            ProgressLocal(chain->GetEntries(),
                          chain->GetReadEntry()+1);
            if ((newquery->fStatus != TQueryDescription::kSessionQueryAborted) &&
                (newquery->fStatus != TQueryDescription::kSessionQueryStopped))
               newquery->fStatus = TQueryDescription::kSessionQueryCompleted;
            UpdateButtons(newquery);
         }
         else {
            new TGMsgBox(fClient->GetRoot(), this, "Error Submitting Query",
                         "Only TChains are allowed in Local Session (no TDSet) !",
                          kMBIconExclamation,kMBOk,&retval);
         }
      }
      else {
         Error("Submit", "No TChain defined; skipping");
         newquery->fStatus = TQueryDescription::kSessionQueryCreated;
         return;
      }
      // set query reference id to unique identifier
      newquery->fReference = TString::Format("local-session-%s:q%lld", newquery->fQueryName.Data(), id);
   }
   // update buttons state
   UpdateButtons(newquery);
}

////////////////////////////////////////////////////////////////////////////////
/// Update buttons state for the current query status.

void TSessionQueryFrame::UpdateButtons(TQueryDescription *desc)
{
   TGListTreeItem *item = fViewer->GetSessionHierarchy()->GetSelected();
   if (!item) return;
   // retrieve query description attached to list tree item
   TObject *obj = (TObject *)item->GetUserData();
   if (!obj || obj->IsA() != TQueryDescription::Class())
      return;
   TQueryDescription *query = (TQueryDescription *)obj;
   if (desc != query) return;

   Bool_t submit_en = kFALSE;
   if ((fViewer->GetActDesc()->fProof &&
       fViewer->GetActDesc()->fProof->IsValid()) ||
       fViewer->GetActDesc()->fLocal)
      submit_en = kTRUE;

   switch (desc->fStatus) {
      case TQueryDescription::kSessionQueryFromProof:
         fBtnSubmit->SetEnabled(submit_en);
         fBtnFinalize->SetEnabled(kTRUE);
         fBtnStop->SetEnabled(kFALSE);
         fBtnAbort->SetEnabled(kFALSE);
         fBtnShowLog->SetEnabled(kTRUE);
         fBtnRetrieve->SetEnabled(kTRUE);
         break;

      case TQueryDescription::kSessionQueryCompleted:
         fBtnSubmit->SetEnabled(submit_en);
         fBtnFinalize->SetEnabled(kTRUE);
         if (((desc->fResult == 0) || (desc->fResult &&
              (desc->fResult->IsFinalized() ||
              (desc->fResult->GetInputObject("TDSet") == 0)))) &&
              !(fViewer->GetActDesc()->fLocal))
            fBtnFinalize->SetEnabled(kFALSE);
         fBtnStop->SetEnabled(kFALSE);
         fBtnAbort->SetEnabled(kFALSE);
         fBtnShowLog->SetEnabled(kTRUE);
         fBtnRetrieve->SetEnabled(kTRUE);
         break;

      case TQueryDescription::kSessionQueryCreated:
         fBtnSubmit->SetEnabled(submit_en);
         fBtnFinalize->SetEnabled(kFALSE);
         fBtnStop->SetEnabled(kFALSE);
         fBtnAbort->SetEnabled(kFALSE);
         fBtnShowLog->SetEnabled(kTRUE);
         fBtnRetrieve->SetEnabled(kFALSE);
         break;

      case TQueryDescription::kSessionQuerySubmitted:
         fBtnSubmit->SetEnabled(kFALSE);
         fBtnFinalize->SetEnabled(kFALSE);
         fBtnStop->SetEnabled(kTRUE);
         fBtnAbort->SetEnabled(kTRUE);
         fBtnShowLog->SetEnabled(kTRUE);
         fBtnRetrieve->SetEnabled(kFALSE);
         break;

      case TQueryDescription::kSessionQueryRunning:
         fBtnSubmit->SetEnabled(kFALSE);
         fBtnFinalize->SetEnabled(kFALSE);
         fBtnStop->SetEnabled(kTRUE);
         fBtnAbort->SetEnabled(kTRUE);
         fBtnShowLog->SetEnabled(kTRUE);
         fBtnRetrieve->SetEnabled(kFALSE);
         break;

      case TQueryDescription::kSessionQueryStopped:
         fBtnSubmit->SetEnabled(submit_en);
         fBtnFinalize->SetEnabled(kTRUE);
         fBtnStop->SetEnabled(kFALSE);
         fBtnAbort->SetEnabled(kFALSE);
         fBtnShowLog->SetEnabled(kTRUE);
         fBtnRetrieve->SetEnabled(kTRUE);
         break;

      case TQueryDescription::kSessionQueryAborted:
         fBtnSubmit->SetEnabled(submit_en);
         fBtnFinalize->SetEnabled(kFALSE);
         fBtnStop->SetEnabled(kFALSE);
         fBtnAbort->SetEnabled(kFALSE);
         fBtnShowLog->SetEnabled(kTRUE);
         fBtnRetrieve->SetEnabled(kFALSE);
         break;

      case TQueryDescription::kSessionQueryFinalized:
         fBtnSubmit->SetEnabled(submit_en);
         fBtnFinalize->SetEnabled(kFALSE);
         fBtnStop->SetEnabled(kFALSE);
         fBtnAbort->SetEnabled(kFALSE);
         fBtnShowLog->SetEnabled(kTRUE);
         fBtnRetrieve->SetEnabled(kFALSE);
         break;

      default:
         break;
   }
   if (fViewer->GetActDesc()->fLocal &&
       !(fViewer->GetActDesc()->fActQuery->fChain)) {
      fBtnFinalize->SetEnabled(kFALSE);
      fBtnRetrieve->SetEnabled(kFALSE);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Update query information (header) text view.

void TSessionQueryFrame::UpdateInfos()
{
   TString buffer;
   const char *qst[] = {"aborted  ", "submitted", "running  ",
                        "stopped  ", "completed"};

   if (fViewer->GetActDesc()->fActQuery)
      fFD->UpdateFields(fViewer->GetActDesc()->fActQuery);

   if (fViewer->GetActDesc()->fLocal ||
      (fViewer->GetActDesc()->fConnected &&
       fViewer->GetActDesc()->fAttached &&
       fViewer->GetActDesc()->fProof &&
       fViewer->GetActDesc()->fProof->IsValid())) {
      fBtnSave->SetText("         Submit         ");
   }
   else {
      fBtnSave->SetText("     Apply changes      ");
   }
   fClient->NeedRedraw(fBtnSave);
   fInfoTextView->Clear();
   if (!fViewer->GetActDesc()->fActQuery ||
       !fViewer->GetActDesc()->fActQuery->fResult) {
      ResetProgressDialog("", 0, 0, 0);
      if (fViewer->GetActDesc()->fLocal) {
         if (fViewer->GetActDesc()->fActQuery) {
            TChain *chain = (TChain *)fViewer->GetActDesc()->fActQuery->fChain;
            if (chain) {
               ProgressLocal(chain->GetEntries(),
                             chain->GetReadEntry()+1);
            }
            else {
               ProgressLocal(0, 0);
            }
            UpdateButtons(fViewer->GetActDesc()->fActQuery);
         }
      }
      else {
         fTotal->SetText(" Estimated time left : 0 sec (0 events of 0 processed)        ");
         fRate->SetText(" Processing Rate : 0.0f events/sec   ");
         frmProg->Reset();
         fFB->Layout();
      }
      return;
   }
   TQueryResult *result = fViewer->GetActDesc()->fActQuery->fResult;

   // Status label
   Int_t st = (result->GetStatus() > 0 && result->GetStatus() <=
               TQueryResult::kCompleted) ? result->GetStatus() : 0;

   Int_t qry = result->GetSeqNum();

   buffer = TString::Format("------------------------------------------------------\n");
   // Print header
   if (!result->IsDraw()) {
      const char *fin = result->IsFinalized() ? "finalized" : qst[st];
      const char *arc = result->IsArchived() ? "(A)" : "";
      buffer = TString::Format(" Query No  : %d\n", qry);
      buffer += TString::Format(" Ref       : \"%s:%s\"\n", result->GetTitle(),
                                result->GetName());
      buffer += TString::Format(" Selector  : %s\n",
                                result->GetSelecImp()->GetTitle());
      buffer += TString::Format("Status    : %9s%s\n", fin, arc);
      buffer += TString::Format("------------------------------------------------------\n");
   } else {
      buffer += TString::Format(" Query No  : %d\n", qry);
      buffer += TString::Format(" Ref       : \"%s:%s\"\n", result->GetTitle(),
                                result->GetName());
      buffer += TString::Format(" Selector  : %s\n",
                                result->GetSelecImp()->GetTitle());
      buffer += TString::Format("------------------------------------------------------\n");
   }

   // Time information
   Int_t elapsed = (Int_t)(result->GetEndTime().Convert() -
                           result->GetStartTime().Convert());
   buffer += TString::Format(" Started   : %s\n",
                             result->GetStartTime().AsString());
   buffer += TString::Format(" Real time : %d sec (CPU time: %.1f sec)\n",
                             elapsed, result->GetUsedCPU());

   // Number of events processed, rate, size
   Double_t rate = 0.0;
   if (result->GetEntries() > -1 && elapsed > 0)
      rate = result->GetEntries() / (Double_t)elapsed ;
   Float_t size = ((Float_t)result->GetBytes())/(1024*1024);
   buffer += TString::Format(" Processed : %lld events (size: %.3f MBs)\n",
                             result->GetEntries(), size);
   buffer += TString::Format(" Rate      : %.1f evts/sec\n", rate);

   // Package information
   if (strlen(result->GetParList()) > 1) {
      buffer += TString::Format(" Packages  :  %s\n", result->GetParList());
   }

   // Result information
   TString res = result->GetResultFile();
   if (!result->IsArchived()) {
      Int_t dq = res.Index("queries");
      if (dq > -1) {
         res.Remove(0,res.Index("queries"));
         res.Insert(0,"<PROOF_SandBox>/");
      }
      if (res.BeginsWith("-")) {
         res = (result->GetStatus() == TQueryResult::kAborted) ?
               "not available" : "sent to client";
      }
   }
   if (res.Length() > 1) {
      buffer += TString::Format("------------------------------------------------------\n");
      buffer += TString::Format(" Results   : %s\n", res.Data());
   }

   if (result->GetOutputList() && result->GetOutputList()->GetSize() > 0) {
      buffer += TString::Format(" Outlist   : %d objects\n",
                                result->GetOutputList()->GetSize());
      buffer += TString::Format("------------------------------------------------------\n");
   }
   fInfoTextView->LoadBuffer(buffer.Data());

   //Float_t pos = Float_t((Double_t)(result->GetEntries() * 100)/(Double_t)total);
   if (result->GetStatus() == TQueryResult::kAborted)
      frmProg->SetBarColor("red");
   else if (result->GetStatus() == TQueryResult::kStopped)
      frmProg->SetBarColor("yellow");
   else
      frmProg->SetBarColor("green");

   frmProg->SetPosition(100.0);

   buffer = TString::Format(" Processed : %lld events in %.1f sec", result->GetEntries(),
                           (Float_t)elapsed);
   fTotal->SetText(buffer.Data());
   buffer = TString::Format(" Processing Rate : %.1f events/sec   ", rate);
   fRate->SetText(buffer.Data());
   fFB->Layout();
}

//////////////////////////////////////////////////////////////////////////////////////////
// Output frame

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TSessionOutputFrame::TSessionOutputFrame(TGWindow* p, Int_t w, Int_t h) :
   TGCompositeFrame(p, w, h), fEntryTmp(0), fLVContainer(0), fViewer(0)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TSessionOutputFrame::~TSessionOutputFrame()
{
   delete fLVContainer; // this container is inside the TGListView and is not
                        // deleted automatically
   Cleanup();
}

////////////////////////////////////////////////////////////////////////////////
/// Build query output information frame.

void TSessionOutputFrame::Build(TSessionViewer *gui)
{
   fViewer = gui;
   SetLayoutManager(new TGVerticalLayout(this));
   SetCleanup(kDeepCleanup);

   // Container of object TGListView
   TGListView *frmListView = new TGListView(this, 340, 190);
   fLVContainer = new TGLVContainer(frmListView, kSunkenFrame, GetWhitePixel());
   fLVContainer->Associate(frmListView);
   fLVContainer->SetCleanup(kDeepCleanup);
   AddFrame(frmListView, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY,
         4, 4, 4, 4));

   frmListView->Connect("Clicked(TGLVEntry*, Int_t, Int_t, Int_t)",
         "TSessionOutputFrame", this,
         "OnElementClicked(TGLVEntry* ,Int_t, Int_t, Int_t)");
   frmListView->Connect("DoubleClicked(TGLVEntry*, Int_t, Int_t, Int_t)",
         "TSessionOutputFrame", this,
         "OnElementDblClicked(TGLVEntry* ,Int_t, Int_t, Int_t)");
}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse clicks on list view items.

void TSessionOutputFrame::OnElementClicked(TGLVEntry* entry, Int_t btn, Int_t x,
                                           Int_t y)
{
   TObject *obj = (TObject *)entry->GetUserData();
   if ((obj) && (btn ==3)) {
      // if right button, popup context menu
      fViewer->GetContextMenu()->Popup(x, y, obj, (TBrowser *)0);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Handle double-clicks on list view items.

void TSessionOutputFrame::OnElementDblClicked(TGLVEntry* entry, Int_t , Int_t, Int_t)
{
   char action[512];
   TString act;
   TObject *obj = (TObject *)entry->GetUserData();
   if (!obj) return;
   TString ext = obj->GetName();
   gPad->SetEditable(kFALSE);
   // check default action from root.mimes
   if (fClient->GetMimeTypeList()->GetAction(obj->IsA()->GetName(), action)) {
      act = TString::Format("((%s*)0x%lx)%s", obj->IsA()->GetName(), (Long_t)obj, action);
      if (act[0] == '!') {
         act.Remove(0, 1);
         gSystem->Exec(act.Data());
      } else {
         // do not allow browse
         if (!act.Contains("Browse"))
            gROOT->ProcessLine(act.Data());
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Add object to output list view.

void TSessionOutputFrame::AddObject(TObject *obj)
{
   TGLVEntry *item;
   if (obj) {
      item = new TGLVEntry(fLVContainer, obj->GetName(), obj->IsA()->GetName());
      item->SetUserData(obj);
      fLVContainer->AddItem(item);
   }
}

//////////////////////////////////////////////////////////////////////////////////////////
// Input Frame

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TSessionInputFrame::TSessionInputFrame(TGWindow* p, Int_t w, Int_t h) :
   TGCompositeFrame(p, w, h), fViewer(0), fLVContainer(0)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TSessionInputFrame::~TSessionInputFrame()
{
   delete fLVContainer; // this container is inside the TGListView and is not
                        // deleted automatically
   Cleanup();
}

////////////////////////////////////////////////////////////////////////////////
/// Build query input information frame.

void TSessionInputFrame::Build(TSessionViewer *gui)
{
   fViewer = gui;
   SetLayoutManager(new TGVerticalLayout(this));
   SetCleanup(kDeepCleanup);

   // Container of object TGListView
   TGListView *frmListView = new TGListView(this, 340, 190);
   fLVContainer = new TGLVContainer(frmListView, kSunkenFrame, GetWhitePixel());
   fLVContainer->Associate(frmListView);
   fLVContainer->SetCleanup(kDeepCleanup);
   AddFrame(frmListView, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY,
         4, 4, 4, 4));
}

////////////////////////////////////////////////////////////////////////////////
/// Add object to input list view.

void TSessionInputFrame::AddObject(TObject *obj)
{
   TGLVEntry *item;
   if (obj) {
      item = new TGLVEntry(fLVContainer, obj->GetName(), obj->IsA()->GetName());
      item->SetUserData(obj);
      fLVContainer->AddItem(item);
   }
}

//////////////////////////////////////////////////////////////////////////////////////////
// Session Viewer Main Frame

////////////////////////////////////////////////////////////////////////////////
/// Main Session viewer constructor.

TSessionViewer::TSessionViewer(const char *name, UInt_t w, UInt_t h) :
   TGMainFrame(gClient->GetRoot(), w, h), fSessionHierarchy(0), fSessionItem(0)
{
   // only one session viewer allowed
   if (gSessionViewer)
      // coverity[uninit_member]: already done
      return;
   Build();
   SetWindowName(name);
   Resize(w, h);
   gSessionViewer = this;
}

////////////////////////////////////////////////////////////////////////////////
/// Main Session viewer constructor.

TSessionViewer::TSessionViewer(const char *name, Int_t x, Int_t y, UInt_t w,
                              UInt_t h) : TGMainFrame(gClient->GetRoot(), w, h),
                              fSessionHierarchy(0), fSessionItem(0)
{
   // only one session viewer allowed
   if (gSessionViewer)
      // coverity[uninit_member]: already done
      return;
   Build();
   SetWindowName(name);
   Move(x, y);
   Resize(w, h);
   gSessionViewer = this;
}

////////////////////////////////////////////////////////////////////////////////
/// Read configuration file and populate list of sessions
/// list of queries and list of packages.
/// Read and set also global options as feedback histos.

void TSessionViewer::ReadConfiguration(const char *filename)
{
   if (fViewerEnv)
      delete fViewerEnv;
   fViewerEnv = new TEnv();

   const char *fn = (filename && strlen(filename)) ? filename : fConfigFile.Data();

   fViewerEnv->ReadFile(fn, kEnvUser);

   Bool_t bval = (Bool_t)fViewerEnv->GetValue("Option.Feedback", 1);
   if (bval)
      fOptionsMenu->CheckEntry(kOptionsFeedback);
   else
      fOptionsMenu->UnCheckEntry(kOptionsFeedback);

   bval = (Bool_t)fViewerEnv->GetValue("Option.MasterHistos", 1);
   if (bval) {
      fOptionsMenu->CheckEntry(kOptionsStatsHist);
      gEnv->SetValue("Proof.StatsHist", 1);
   }
   else {
      fOptionsMenu->UnCheckEntry(kOptionsStatsHist);
      gEnv->SetValue("Proof.StatsHist", 0);
   }
   bval = (Bool_t)fViewerEnv->GetValue("Option.MasterEvents", 0);
   if (bval)
      fOptionsMenu->CheckEntry(kOptionsStatsTrace);
   else
      fOptionsMenu->UnCheckEntry(kOptionsStatsTrace);
   bval = (Bool_t)fViewerEnv->GetValue("Option.WorkerEvents", 0);
   if (bval)
      fOptionsMenu->CheckEntry(kOptionsSlaveStatsTrace);
   else
      fOptionsMenu->UnCheckEntry(kOptionsSlaveStatsTrace);

   Int_t i = 0;
   while (kFeedbackHistos[i]) {
      bval = (Bool_t)fViewerEnv->GetValue(Form("Option.%s",kFeedbackHistos[i]),
                                          i == 1 ? 1 : 0);
      if (bval)
         fCascadeMenu->CheckEntry(41+i);
      else
         fCascadeMenu->UnCheckEntry(41+i);
      i++;
   }
   TSessionDescription *proofDesc;
   fSessions->Delete();
   if (fSessionItem)
      fSessionHierarchy->DeleteChildren(fSessionItem);
   else
      fSessionItem = fSessionHierarchy->AddItem(0, "Sessions", fBaseIcon,
            fBaseIcon);
   // add local session description
   TGListTreeItem *item = fSessionHierarchy->AddItem(fSessionItem, "Local",
                                                     fLocal, fLocal);
   fSessionHierarchy->SetToolTipItem(item, "Local Session");
   TSessionDescription *localdesc = new TSessionDescription();
   localdesc->fTag = "";
   localdesc->fName = "Local";
   localdesc->fAddress = "Local";
   localdesc->fPort = 0;
   localdesc->fConfigFile = "";
   localdesc->fLogLevel = 0;
   localdesc->fUserName = "";
   localdesc->fQueries = new TList();
   localdesc->fPackages = new TList();
   localdesc->fActQuery = 0;
   localdesc->fProof = 0;
   localdesc->fProofMgr = 0;
   localdesc->fAttached = kFALSE;
   localdesc->fConnected = kFALSE;
   localdesc->fLocal = kTRUE;
   localdesc->fSync = kTRUE;
   localdesc->fAutoEnable = kFALSE;
   localdesc->fNbHistos = 0;
   item->SetUserData(localdesc);
   fSessions->Add((TObject *)localdesc);
   fActDesc = localdesc;

   SysInfo_t info;
   gSystem->GetSysInfo(&info);
   // if the machine has more than one CPU, add one PROOF lite session
   // (not supported on Windows yet)
   if (!info.fOS.Contains("Microsoft") && info.fCpus > 1) {
      // add proof lite session description
      item = fSessionHierarchy->AddItem(fSessionItem, "Lite",
                                        fProofDiscon, fProofDiscon);
      fSessionHierarchy->SetToolTipItem(item, "PROOF Lite");
      TSessionDescription *litedesc = new TSessionDescription();
      litedesc->fTag = "";
      litedesc->fName = "PROOF Lite";
      litedesc->fAddress = "lite://";
      litedesc->fPort = 0;
      litedesc->fConfigFile = "";
      litedesc->fLogLevel = 0;
      litedesc->fUserName = "";
      litedesc->fQueries = new TList();
      litedesc->fPackages = new TList();
      litedesc->fActQuery = 0;
      litedesc->fProof = 0;
      litedesc->fProofMgr = 0;
      litedesc->fAttached = kFALSE;
      litedesc->fConnected = kFALSE;
      litedesc->fLocal = kFALSE;
      litedesc->fSync = kTRUE;
      litedesc->fAutoEnable = kFALSE;
      litedesc->fNbHistos = 0;
      item->SetUserData(litedesc);
      fSessions->Add((TObject *)litedesc);
      fActDesc = litedesc;
   }
   TIter next(fViewerEnv->GetTable());
   TEnvRec *er;
   while ((er = (TEnvRec*) next())) {
      const char *s;
      if ((s = strstr(er->GetName(), "SessionDescription."))) {
         const char *val = fViewerEnv->GetValue(s, (const char*)0);
         if (val) {
            Int_t cnt = 0;
            char *v = StrDup(val);
            s += 7;
            while (1) {
               TString name = strtok(!cnt ? v : 0, ";");
               if (name.IsNull()) break;
               TString sessiontag = strtok(0, ";");
               TString address = strtok(0, ";");
               if (address.IsNull()) break;
               TString port = strtok(0, ";");
               if (port.IsNull()) break;
               TString loglevel = strtok(0, ";");
               if (loglevel.IsNull()) break;
               TString configfile = strtok(0, ";");
               TString user = strtok(0, ";");
               if (user.IsNull()) break;
               TString sync = strtok(0, ";");
               TString autoen = strtok(0, ";");

               // build session description
               proofDesc = new TSessionDescription();
               proofDesc->fTag = sessiontag.Length() > 2 ? sessiontag.Data() : "";
               proofDesc->fName = name;
               proofDesc->fAddress = address;
               proofDesc->fPort = atoi(port);
               proofDesc->fConfigFile = configfile.Length() > 2 ? configfile.Data() : "";
               proofDesc->fLogLevel = atoi(loglevel);
               proofDesc->fConnected = kFALSE;
               proofDesc->fAttached = kFALSE;
               proofDesc->fLocal = kFALSE;
               proofDesc->fQueries = new TList();
               proofDesc->fPackages = new TList();
               proofDesc->fActQuery = 0;
               proofDesc->fProof = 0;
               proofDesc->fProofMgr = 0;
               proofDesc->fSync = (Bool_t)(atoi(sync));
               proofDesc->fAutoEnable = (Bool_t)(atoi(autoen));
               proofDesc->fUserName = user;
               fSessions->Add((TObject *)proofDesc);
               item = fSessionHierarchy->AddItem(
                     fSessionItem, proofDesc->fName.Data(),
                     fProofDiscon, fProofDiscon);
               fSessionHierarchy->SetToolTipItem(item, "Proof Session");
               item->SetUserData(proofDesc);
               fActDesc = proofDesc;
               cnt++;
            }
            delete [] v;
         }
      }
      if ((s = strstr(er->GetName(), "QueryDescription."))) {
         const char *val = fViewerEnv->GetValue(s, (const char*)0);
         if (val) {
            Int_t cnt = 0;
            char *v = StrDup(val);
            s += 7;
            while (1) {

               TString status = strtok(!cnt ? v : 0, ";");
               if (status.IsNull()) break;

               TString reference = strtok(0, ";");
               if (reference.IsNull()) break;
               TString queryname = strtok(0, ";");
               if (queryname.IsNull()) break;
               TString selector = strtok(0, ";");
               if (selector.IsNull()) break;
               TString dset = strtok(0, ";");
               TString options = strtok(0, ";");
               TString eventlist = strtok(0, ";");
               TString nbfiles = strtok(0, ";");
               TString nbentries = strtok(0, ";");
               TString firstentry = strtok(0, ";");

               TQueryDescription *newquery = new TQueryDescription();
               newquery->fStatus =
                  (TQueryDescription::ESessionQueryStatus)(atoi(status));
               newquery->fSelectorString  = selector.Length() > 2 ? selector.Data() : "";
               newquery->fReference       = reference.Length() > 2 ? reference.Data() : "";
               newquery->fTDSetString     = dset.Length() > 2 ? dset.Data() : "";
               newquery->fQueryName       = queryname.Length() > 2 ? queryname.Data() : "";
               newquery->fOptions         = options.Length() > 2 ? options.Data() : "";
               newquery->fEventList       = eventlist.Length() > 2 ? eventlist.Data() : "";
               newquery->fNbFiles         = atoi(nbfiles);
               newquery->fNoEntries       = atoi(nbentries);
               newquery->fFirstEntry      = atoi(firstentry);
               newquery->fResult          = 0;
               newquery->fChain           = 0;
               fActDesc->fQueries->Add((TObject *)newquery);
               cnt++;
               TGListTreeItem *item1 = fSessionHierarchy->FindChildByData(
                     fSessionItem, fActDesc);
               TGListTreeItem *item2 = fSessionHierarchy->AddItem(
                     item1, newquery->fQueryName, fQueryCon, fQueryCon);
               item2->SetUserData(newquery);
            }
            delete [] v;
         }
      }
   }
   fSessionHierarchy->ClearHighlighted();
   fSessionHierarchy->OpenItem(fSessionItem);
   if (fActDesc == localdesc) {
      fSessionHierarchy->HighlightItem(fSessionItem);
      fSessionHierarchy->SetSelected(fSessionItem);
   }
   else {
      fSessionHierarchy->OpenItem(item);
      fSessionHierarchy->HighlightItem(item);
      fSessionHierarchy->SetSelected(item);
   }
   fSessionHierarchy->ClearViewPort();
   fClient->NeedRedraw(fSessionHierarchy);
}

////////////////////////////////////////////////////////////////////////////////
/// Update list of existing Proof sessions.

void TSessionViewer::UpdateListOfProofs()
{
   // get list of proof sessions
   Bool_t found  = kFALSE;
   Bool_t exists = kFALSE;
   TGListTreeItem *item = 0;
   TSeqCollection *proofs = gROOT->GetListOfProofs();
   TSessionDescription *desc = 0;
   TSessionDescription *newdesc;
   if (proofs) {
      TObject *o = proofs->First();
      if (o && dynamic_cast<TProofMgr *>(o)) {
         TProofMgr *mgr = dynamic_cast<TProofMgr *>(o);
         if (mgr && mgr->QuerySessions("L")) {
            TIter nxd(mgr->QuerySessions("L"));
            TProofDesc *d = 0;
            TProof *p = 0;
            while ((d = (TProofDesc *)nxd())) {
               TIter nextfs(fSessions);
               // check if session exists in the list
               exists = kFALSE;
               while ((desc = (TSessionDescription *)nextfs())) {
                  if ((desc->fTag == d->GetName()) ||
                      (desc->fName == d->GetTitle())) {
                     exists = kTRUE;
                     break;
                  }
               }
               TIter nexts(fSessions);
               found = kFALSE;
               p = d->GetProof();
               while ((desc = (TSessionDescription *)nexts())) {
                  if (desc->fConnected && desc->fAttached)
                     continue;
                  if (p && ((exists && ((desc->fTag == d->GetName()) ||
                      (desc->fName == d->GetTitle()))) ||
                      (!exists && (desc->fAddress == p->GetMaster())))) {
                     desc->fConnected  = kTRUE;
                     desc->fAttached   = kTRUE;
                     desc->fProof      = p;
                     desc->fProofMgr   = mgr;
                     desc->fTag        = d->GetName();
                     item = fSessionHierarchy->FindChildByData(fSessionItem,
                                                               desc);
                     if (item) {
                        item->SetPictures(fProofCon, fProofCon);
                        if (item == fSessionHierarchy->GetSelected()) {
                           fActDesc->fProof->Connect("Progress(Long64_t,Long64_t)",
                                    "TSessionQueryFrame", fQueryFrame,
                                    "Progress(Long64_t,Long64_t)");
                           fActDesc->fProof->Connect("Progress(Long64_t,Long64_t,Long64_t,Float_t,Float_t,Float_t,Float_t)",
                                    "TSessionQueryFrame", fQueryFrame,
                                    "Progress(Long64_t,Long64_t,Long64_t,Float_t,Float_t,Float_t,Float_t)");
                           fActDesc->fProof->Connect("StopProcess(Bool_t)",
                                    "TSessionQueryFrame", fQueryFrame,
                                    "IndicateStop(Bool_t)");
                           fActDesc->fProof->Connect(
                              "ResetProgressDialog(const char*, Int_t,Long64_t,Long64_t)",
                              "TSessionQueryFrame", fQueryFrame,
                              "ResetProgressDialog(const char*,Int_t,Long64_t,Long64_t)");
                           // enable timer used for status bar icon's animation
                           EnableTimer();
                           // change status bar right icon to connected pixmap
                           ChangeRightLogo("monitor01.xpm");
                           // do not animate yet
                           SetChangePic(kFALSE);
                           // connect to signal "query result ready"
                           fActDesc->fProof->Connect("QueryResultReady(char *)",
                                    "TSessionViewer", this, "QueryResultReady(char *)");
                           // display connection information on status bar
                           TString msg;
                           msg.Form("PROOF Cluster %s ready", fActDesc->fName.Data());
                           fStatusBar->SetText(msg.Data(), 1);
                           UpdateListOfPackages();
                           fSessionFrame->UpdatePackages();
                           fSessionFrame->UpdateListOfDataSets();
                           fPopupSrv->DisableEntry(kSessionConnect);
                           fSessionMenu->DisableEntry(kSessionConnect);
                           fPopupSrv->EnableEntry(kSessionDisconnect);
                           fSessionMenu->EnableEntry(kSessionDisconnect);
                           if (fToolBar->GetButton(kSessionDisconnect))
                              fToolBar->GetButton(kSessionDisconnect)->SetState(kButtonUp);
                           if (fToolBar->GetButton(kSessionConnect))
                              fToolBar->GetButton(kSessionConnect)->SetState(kButtonDisabled);
                           fSessionFrame->SetLogLevel(fActDesc->fLogLevel);
                           // update session information frame
                           fSessionFrame->ProofInfos();
                           fSessionFrame->SetLocal(kFALSE);
                           if (fActFrame != fSessionFrame) {
                              fV2->HideFrame(fActFrame);
                              fV2->ShowFrame(fSessionFrame);
                              fActFrame = fSessionFrame;
                           }
                        }
                     }
                     if (desc->fLogLevel < 0)
                        desc->fLogLevel = 0;
                     found = kTRUE;
                     break;
                  }
               }
               if (found) continue;
               newdesc = new TSessionDescription();
               // and fill information from Proof session
               newdesc->fTag       = d->GetName();
               newdesc->fName      = d->GetTitle();
               newdesc->fAddress   = d->GetTitle();
               newdesc->fConnected = kFALSE;
               newdesc->fAttached  = kFALSE;
               newdesc->fProofMgr  = mgr;
               p = d->GetProof();
               if (p) {
                  newdesc->fConnected  = kTRUE;
                  newdesc->fAttached   = kTRUE;
                  if (p->IsLite())
                     newdesc->fAddress = "lite://";
                  else
                     newdesc->fAddress = p->GetMaster();
                  newdesc->fConfigFile = p->GetConfFile();
                  newdesc->fUserName   = p->GetUser();
                  newdesc->fPort       = p->GetPort();
                  newdesc->fLogLevel   = p->GetLogLevel();
                  newdesc->fProof      = p;
                  newdesc->fProof->Connect("Progress(Long64_t,Long64_t)",
                           "TSessionQueryFrame", fQueryFrame,
                           "Progress(Long64_t,Long64_t)");
                  newdesc->fProof->Connect("Progress(Long64_t,Long64_t,Long64_t,Float_t,Float_t,Float_t,Float_t)",
                           "TSessionQueryFrame", fQueryFrame,
                           "Progress(Long64_t,Long64_t,Long64_t,Float_t,Float_t,Float_t,Float_t)");
                  newdesc->fProof->Connect("StopProcess(Bool_t)",
                           "TSessionQueryFrame", fQueryFrame,
                           "IndicateStop(Bool_t)");
                  newdesc->fProof->Connect(
                           "ResetProgressDialog(const char*, Int_t,Long64_t,Long64_t)",
                           "TSessionQueryFrame", fQueryFrame,
                           "ResetProgressDialog(const char*,Int_t,Long64_t,Long64_t)");
                  // enable timer used for status bar icon's animation
                  EnableTimer();
                  // change status bar right icon to connected pixmap
                  ChangeRightLogo("monitor01.xpm");
                  // do not animate yet
                  SetChangePic(kFALSE);
                  // connect to signal "query result ready"
                  newdesc->fProof->Connect("QueryResultReady(char *)",
                           "TSessionViewer", this, "QueryResultReady(char *)");
               }
               newdesc->fQueries    = new TList();
               newdesc->fPackages   = new TList();
               if (newdesc->fLogLevel < 0)
                  newdesc->fLogLevel = 0;
               newdesc->fActQuery   = 0;
               newdesc->fLocal = kFALSE;
               newdesc->fSync = kFALSE;
               newdesc->fAutoEnable = kFALSE;
               newdesc->fNbHistos = 0;
               // add new session description in list tree
               if (p)
                  item = fSessionHierarchy->AddItem(fSessionItem, newdesc->fName.Data(),
                           fProofCon, fProofCon);
               else
                  item = fSessionHierarchy->AddItem(fSessionItem, newdesc->fName.Data(),
                           fProofDiscon, fProofDiscon);
               fSessionHierarchy->SetToolTipItem(item, "Proof Session");
               item ->SetUserData(newdesc);
               // and in our session description list
               fSessions->Add(newdesc);
            }
         }
         return;
      }
      TIter nextp(proofs);
      TProof *proof;
      // loop over existing Proof sessions
      while ((proof = (TProof *)nextp())) {
         TIter nexts(fSessions);
         found = kFALSE;
         // check if session is already in the list
         while ((desc = (TSessionDescription *)nexts())) {
            if (desc->fProof == proof) {
               desc->fConnected = kTRUE;
               desc->fAttached = kTRUE;
               found = kTRUE;
               break;
            }
         }
         if (found) continue;
         // create new session description
         newdesc = new TSessionDescription();
         // and fill information from Proof session
         newdesc->fName       = proof->GetMaster();
         newdesc->fConfigFile = proof->GetConfFile();
         newdesc->fUserName   = proof->GetUser();
         newdesc->fPort       = proof->GetPort();
         newdesc->fLogLevel   = proof->GetLogLevel();
         if (newdesc->fLogLevel < 0)
            newdesc->fLogLevel = 0;
         if (proof->IsLite())
            newdesc->fAddress = "lite://";
         else
            newdesc->fAddress = proof->GetMaster();
         newdesc->fQueries    = new TList();
         newdesc->fPackages   = new TList();
         newdesc->fProof      = proof;
         newdesc->fActQuery   = 0;
         newdesc->fConnected = kTRUE;
         newdesc->fAttached = kTRUE;
         newdesc->fLocal = kFALSE;
         newdesc->fSync = kFALSE;
         newdesc->fAutoEnable = kFALSE;
         newdesc->fNbHistos = 0;
         // add new session description in list tree
         item = fSessionHierarchy->AddItem(fSessionItem, newdesc->fName.Data(),
                  fProofCon, fProofCon);
         fSessionHierarchy->SetToolTipItem(item, "Proof Session");
         item ->SetUserData(newdesc);
         // and in our session description list
         fSessions->Add(newdesc);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Update list of existing Proof sessions.

void TSessionViewer::UpdateListOfSessions()
{
   // get list of proof sessions
   TGListTreeItem *item;
   TList *sessions = fActDesc->fProofMgr->QuerySessions("");
   if (sessions) {
      TIter nextp(sessions);
      TProofDesc *pdesc;
      TProof *proof;
      TSessionDescription *newdesc;
      // loop over existing Proof sessions
      while ((pdesc = (TProofDesc *)nextp())) {
         TIter nexts(fSessions);
         TSessionDescription *desc = 0;
         Bool_t found = kFALSE;
         // check if session is already in the list
         while ((desc = (TSessionDescription *)nexts())) {
            if ((desc->fTag == pdesc->GetName()) ||
                (desc->fName == pdesc->GetTitle())) {
               desc->fConnected = kTRUE;
               found = kTRUE;
               break;
            }
         }
         if (found) continue;
         // create new session description
         newdesc = new TSessionDescription();
         // and fill information from Proof session
         newdesc->fTag        = pdesc->GetName();
         newdesc->fName       = pdesc->GetTitle();
         proof = pdesc->GetProof();
         if (proof) {
            newdesc->fConfigFile = proof->GetConfFile();
            newdesc->fUserName   = proof->GetUser();
            newdesc->fPort       = proof->GetPort();
            newdesc->fLogLevel   = proof->GetLogLevel();
            if (newdesc->fLogLevel < 0)
               newdesc->fLogLevel = 0;
            if (proof->IsLite())
               newdesc->fAddress = "lite://";
            else
               newdesc->fAddress = proof->GetMaster();
            newdesc->fProof      = proof;
         }
         else {
            newdesc->fProof      = 0;
            newdesc->fConfigFile = "";
            newdesc->fUserName   = fActDesc->fUserName;
            newdesc->fPort       = fActDesc->fPort;
            newdesc->fLogLevel   = 0;
            newdesc->fAddress    = fActDesc->fAddress;
         }
         newdesc->fQueries    = new TList();
         newdesc->fPackages   = new TList();
         newdesc->fProofMgr   = fActDesc->fProofMgr;
         newdesc->fActQuery   = 0;
         newdesc->fConnected  = kTRUE;
         newdesc->fAttached   = kFALSE;
         newdesc->fLocal      = kFALSE;
         newdesc->fSync       = kFALSE;
         newdesc->fAutoEnable = kFALSE;
         newdesc->fNbHistos   = 0;
         // add new session description in list tree
         item = fSessionHierarchy->AddItem(fSessionItem, newdesc->fName.Data(),
                  fProofDiscon, fProofDiscon);
         fSessionHierarchy->SetToolTipItem(item, "Proof Session");
         item ->SetUserData(newdesc);
         // and in our session description list
         fSessions->Add(newdesc);
         // set actual description to the last one
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Save actual configuration in config file "filename".

void TSessionViewer::WriteConfiguration(const char *filename)
{
   TSessionDescription *session;
   TQueryDescription *query;
   Int_t scnt = 0, qcnt = 1;
   const char *fname = filename ? filename : fConfigFile.Data();

   delete fViewerEnv;
   gSystem->Unlink(fname);
   fViewerEnv = new TEnv();

   fViewerEnv->SetValue("Option.Feedback",
         (Int_t)fOptionsMenu->IsEntryChecked(kOptionsFeedback));
   fViewerEnv->SetValue("Option.MasterHistos",
         (Int_t)fOptionsMenu->IsEntryChecked(kOptionsStatsHist));
   fViewerEnv->SetValue("Option.MasterEvents",
         (Int_t)fOptionsMenu->IsEntryChecked(kOptionsStatsTrace));
   fViewerEnv->SetValue("Option.WorkerEvents",
         (Int_t)fOptionsMenu->IsEntryChecked(kOptionsSlaveStatsTrace));

   Int_t i = 0;
   // browse list of feedback histos and check user's selected ones
   while (kFeedbackHistos[i]) {
      fViewerEnv->SetValue(Form("Option.%s",kFeedbackHistos[i]),
         (Int_t)fCascadeMenu->IsEntryChecked(41+i));
      i++;
   }

   TIter snext(fSessions);
   while ((session = (TSessionDescription *) snext())) {
      if ((scnt > 0) && ((session->fAddress.Length() < 3) ||
           session->fUserName.Length() < 2)) {
         // skip gROOT's list of sessions
         continue;
      }
      if ((scnt > 0) && (session->fName == session->fAddress)) {
         // skip gROOT's list of proofs
         continue;
      }
      TString sessionstring;
      sessionstring += session->fName;
      sessionstring += ";";
      sessionstring += session->fTag.Length() > 1 ? session->fTag.Data() : " ";
      sessionstring += ";";
      sessionstring += session->fAddress;
      sessionstring += ";";
      sessionstring += TString::Format("%d", session->fPort);
      sessionstring += ";";
      sessionstring += TString::Format("%d", session->fLogLevel);
      sessionstring += ";";
      sessionstring += session->fConfigFile.Length() > 1 ? session->fConfigFile.Data() : " ";
      sessionstring += ";";
      sessionstring += session->fUserName;
      sessionstring += ";";
      sessionstring += TString::Format("%d", session->fSync);
      sessionstring += ";";
      sessionstring += TString::Format("%d", session->fAutoEnable);
      if (scnt > 0) // skip local session
         fViewerEnv->SetValue(Form("SessionDescription.%d",scnt), sessionstring);
      scnt++;

      TIter qnext(session->fQueries);
      while ((query = (TQueryDescription *) qnext())) {
         TString querystring;
         querystring += TString::Format("%d", query->fStatus);
         querystring += ";";
         querystring += query->fReference.Length() > 1 ? query->fReference.Data() : " ";
         querystring += ";";
         querystring += query->fQueryName;
         querystring += ";";
         querystring += query->fSelectorString.Length() > 1 ? query->fSelectorString.Data() : " ";
         querystring += ";";
         querystring += query->fTDSetString.Length() > 1 ? query->fTDSetString.Data() : " ";
         querystring += ";";
         querystring += query->fOptions.Length() > 1 ? query->fOptions.Data() : " ";
         querystring += ";";
         querystring += query->fEventList.Length() > 1 ? query->fEventList.Data() : " ";
         querystring += ";";
         querystring += TString::Format("%d",query->fNbFiles);
         querystring += ";";
         querystring += TString::Format("%lld",query->fNoEntries);
         querystring += ";";
         querystring += TString::Format("%lld",query->fFirstEntry);
         fViewerEnv->SetValue(Form("QueryDescription.%d",qcnt), querystring);
         qcnt++;
      }
   }

   fViewerEnv->WriteFile(fname);
}

////////////////////////////////////////////////////////////////////////////////
/// Build main session viewer frame and subframes.

void TSessionViewer::Build()
{
   TString buf;
   fActDesc = 0;
   fActFrame = 0;
   fLogWindow = 0;
   fBusy = kFALSE;
   fAutoSave = kTRUE;
   fChangePic = kFALSE;
   fStart = fElapsed = 0;

   SetCleanup(kDeepCleanup);
   // set minimun size
   SetWMSizeHints(400 + 200, 370+50, 2000, 1000, 1, 1);

   // collect icons
   fLocal = fClient->GetPicture("local_session.xpm");
   fProofCon = fClient->GetPicture("proof_connected.xpm");
   fProofDiscon = fClient->GetPicture("proof_disconnected.xpm");
   fQueryCon = fClient->GetPicture("query_connected.xpm");
   fQueryDiscon = fClient->GetPicture("query_disconnected.xpm");
   fBaseIcon = fClient->GetPicture("proof_base.xpm");

   //--- File menu
   fFileMenu = new TGPopupMenu(fClient->GetDefaultRoot());
   fFileMenu->AddEntry("&Load Config...", kFileLoadConfig);
   fFileMenu->AddEntry("&Save Config...", kFileSaveConfig);
   fFileMenu->AddSeparator();
   fFileMenu->AddEntry("&Close Viewer",    kFileCloseViewer);
   fFileMenu->AddSeparator();
   fFileMenu->AddEntry("&Quit ROOT",       kFileQuit);

   //--- Session menu
   fSessionMenu = new TGPopupMenu(gClient->GetDefaultRoot());
   fSessionMenu->AddLabel("Session Management");
   fSessionMenu->AddSeparator();
   fSessionMenu->AddEntry("&New Session", kSessionNew);
   fSessionMenu->AddEntry("&Add to the list", kSessionAdd);
   fSessionMenu->AddEntry("De&lete", kSessionDelete);
   fSessionMenu->AddSeparator();
   fSessionMenu->AddEntry("&Connect...", kSessionConnect);
   fSessionMenu->AddEntry("&Disconnect", kSessionDisconnect);
   fSessionMenu->AddEntry("Shutdo&wn",  kSessionShutdown);
   fSessionMenu->AddEntry("&Show status",kSessionShowStatus);
   fSessionMenu->AddEntry("&Get Queries",kSessionGetQueries);
   fSessionMenu->AddSeparator();
   fSessionMenu->AddEntry("&Cleanup", kSessionCleanup);
   fSessionMenu->AddEntry("&Reset",kSessionReset);
   fSessionMenu->DisableEntry(kSessionAdd);

   //--- Query menu
   fQueryMenu = new TGPopupMenu(gClient->GetDefaultRoot());
   fQueryMenu->AddLabel("Query Management");
   fQueryMenu->AddSeparator();
   fQueryMenu->AddEntry("&New...", kQueryNew);
   fQueryMenu->AddEntry("&Edit", kQueryEdit);
   fQueryMenu->AddEntry("&Submit", kQuerySubmit);
   fQueryMenu->AddSeparator();
   fQueryMenu->AddEntry("Start &Viewer", kQueryStartViewer);
   fQueryMenu->AddSeparator();
   fQueryMenu->AddEntry("&Delete", kQueryDelete);

   fViewerEnv = 0;
#ifdef WIN32
   fConfigFile = TString::Format("%s\\%s", gSystem->HomeDirectory(), kConfigFile);
#else
   fConfigFile = TString::Format("%s/%s", gSystem->HomeDirectory(), kConfigFile);
#endif

   fCascadeMenu = new TGPopupMenu(fClient->GetDefaultRoot());
   Int_t i = 0;
   while (kFeedbackHistos[i]) {
      fCascadeMenu->AddEntry(kFeedbackHistos[i], 41+i);
      i++;
   }
   fCascadeMenu->AddEntry("User defined...", 50);
   // disable it for now (until implemented)
   fCascadeMenu->DisableEntry(50);

   //--- Options menu
   fOptionsMenu = new TGPopupMenu(fClient->GetDefaultRoot());
   fOptionsMenu->AddLabel("Global Options");
   fOptionsMenu->AddSeparator();
   fOptionsMenu->AddEntry("&Autosave Config", kOptionsAutoSave);
   fOptionsMenu->AddSeparator();
   fOptionsMenu->AddEntry("Master &Histos", kOptionsStatsHist);
   fOptionsMenu->AddEntry("&Master Events", kOptionsStatsTrace);
   fOptionsMenu->AddEntry("&Worker Events", kOptionsSlaveStatsTrace);
   fOptionsMenu->AddSeparator();
   fOptionsMenu->AddEntry("Feedback &Active", kOptionsFeedback);
   fOptionsMenu->AddSeparator();
   fOptionsMenu->AddPopup("&Feedback Histos", fCascadeMenu);
   fOptionsMenu->CheckEntry(kOptionsAutoSave);

   //--- Help menu
   fHelpMenu = new TGPopupMenu(gClient->GetDefaultRoot());
   fHelpMenu->AddEntry("&About ROOT...",  kHelpAbout);

   fFileMenu->Associate(this);
   fSessionMenu->Associate(this);
   fQueryMenu->Associate(this);
   fOptionsMenu->Associate(this);
   fCascadeMenu->Associate(this);
   fHelpMenu->Associate(this);

   //--- create menubar and add popup menus
   fMenuBar = new TGMenuBar(this, 1, 1, kHorizontalFrame);

   fMenuBar->AddPopup("&File", fFileMenu, new TGLayoutHints(kLHintsTop |
         kLHintsLeft, 0, 4, 0, 0));
   fMenuBar->AddPopup("&Session", fSessionMenu, new TGLayoutHints(kLHintsTop |
         kLHintsLeft, 0, 4, 0, 0));
   fMenuBar->AddPopup("&Query",  fQueryMenu, new TGLayoutHints(kLHintsTop |
         kLHintsLeft, 0, 4, 0, 0));
   fMenuBar->AddPopup("&Options",  fOptionsMenu, new TGLayoutHints(kLHintsTop |
         kLHintsLeft, 0, 4, 0, 0));
   fMenuBar->AddPopup("&Help", fHelpMenu, new TGLayoutHints(kLHintsTop |
         kLHintsRight));

   TGHorizontal3DLine *toolBarSep = new TGHorizontal3DLine(this);
   AddFrame(toolBarSep, new TGLayoutHints(kLHintsTop | kLHintsExpandX));

   AddFrame(fMenuBar, new TGLayoutHints(kLHintsTop | kLHintsLeft |
         kLHintsExpandX, 0, 0, 1, 1));

   toolBarSep = new TGHorizontal3DLine(this);
   AddFrame(toolBarSep, new TGLayoutHints(kLHintsTop | kLHintsExpandX));

   //---- toolbar

   int spacing = 8;
   fToolBar = new TGToolBar(this, 60, 20, kHorizontalFrame);
   for (int ii = 0; xpm_toolbar[ii]; ii++) {
      tb_data[ii].fPixmap = xpm_toolbar[ii];
      if (strlen(xpm_toolbar[ii]) == 0) {
         spacing = 8;
         continue;
      }
      fToolBar->AddButton(this, &tb_data[ii], spacing);
      spacing = 0;
   }
   AddFrame(fToolBar, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 0, 0, 0, 0));
   toolBarSep = new TGHorizontal3DLine(this);
   AddFrame(toolBarSep, new TGLayoutHints(kLHintsTop | kLHintsExpandX));
   if (fToolBar->GetButton(kQuerySubmit))
      fToolBar->GetButton(kQuerySubmit)->SetState(kButtonDisabled);

   fPopupSrv = new TGPopupMenu(fClient->GetDefaultRoot());
   fPopupSrv->AddEntry("Connect",kSessionConnect);
   fPopupSrv->AddEntry("Disconnect",kSessionDisconnect);
   fPopupSrv->AddEntry("Shutdown",kSessionShutdown);
   fPopupSrv->AddEntry("Browse",kSessionBrowse);
   fPopupSrv->AddEntry("Show status",kSessionShowStatus);
   fPopupSrv->AddEntry("Delete", kSessionDelete);
   fPopupSrv->AddEntry("Get Queries",kSessionGetQueries);
   fPopupSrv->AddSeparator();
   fPopupSrv->AddEntry("Cleanup", kSessionCleanup);
   fPopupSrv->AddEntry("Reset",kSessionReset);
   fPopupSrv->Connect("Activated(Int_t)","TSessionViewer", this,
         "MyHandleMenu(Int_t)");

   fPopupQry = new TGPopupMenu(fClient->GetDefaultRoot());
   fPopupQry->AddEntry("Edit",kQueryEdit);
   fPopupQry->AddEntry("Submit",kQuerySubmit);
   fPopupQry->AddSeparator();
   fPopupQry->AddEntry("Start &Viewer", kQueryStartViewer);
   fPopupQry->AddSeparator();
   fPopupQry->AddEntry("Delete",kQueryDelete);
   fPopupQry->Connect("Activated(Int_t)","TSessionViewer", this,
         "MyHandleMenu(Int_t)");


   fSessionMenu->DisableEntry(kSessionGetQueries);
   fSessionMenu->DisableEntry(kSessionShowStatus);
   fPopupSrv->DisableEntry(kSessionGetQueries);
   fPopupSrv->DisableEntry(kSessionShowStatus);
   fPopupSrv->DisableEntry(kSessionDisconnect);
   fPopupSrv->DisableEntry(kSessionShutdown);
   fPopupSrv->DisableEntry(kSessionCleanup);
   fPopupSrv->DisableEntry(kSessionReset);
   fSessionMenu->DisableEntry(kSessionDisconnect);
   fSessionMenu->DisableEntry(kSessionShutdown);
   fSessionMenu->DisableEntry(kSessionCleanup);
   fSessionMenu->DisableEntry(kSessionReset);
   if (fToolBar->GetButton(kSessionDisconnect))
      fToolBar->GetButton(kSessionDisconnect)->SetState(kButtonDisabled);

   //--- Horizontal mother frame -----------------------------------------------
   fHf = new TGHorizontalFrame(this, 10, 10);
   fHf->SetCleanup(kDeepCleanup);

   //--- fV1 -------------------------------------------------------------------
   fV1 = new TGVerticalFrame(fHf, 100, 100, kFixedWidth);
   fV1->SetCleanup(kDeepCleanup);

   fTreeView = new TGCanvas(fV1, 100, 200, kSunkenFrame | kDoubleBorder);
   fV1->AddFrame(fTreeView, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY,
         2, 0, 0, 0));
   fSessionHierarchy = new TGListTree(fTreeView, kHorizontalFrame);
   fSessionHierarchy->DisableOpen();
   fSessionHierarchy->Connect("Clicked(TGListTreeItem*,Int_t,Int_t,Int_t)",
         "TSessionViewer", this,
         "OnListTreeClicked(TGListTreeItem*, Int_t, Int_t, Int_t)");
   fSessionHierarchy->Connect("DoubleClicked(TGListTreeItem*,Int_t)",
         "TSessionViewer", this,
         "OnListTreeDoubleClicked(TGListTreeItem*, Int_t)");
   fV1->Resize(fTreeView->GetDefaultWidth()+100, fV1->GetDefaultHeight());

   //--- fV2 -------------------------------------------------------------------
   fV2 = new TGVerticalFrame(fHf, 350, 310);
   fV2->SetCleanup(kDeepCleanup);

   //--- Server Frame ----------------------------------------------------------
   fServerFrame = new TSessionServerFrame(fV2, 350, 310);
   fSessions = new TList;
   ReadConfiguration();
   fServerFrame->Build(this);
   fV2->AddFrame(fServerFrame, new TGLayoutHints(kLHintsTop | kLHintsExpandX |
         kLHintsExpandY, 2, 0, 1, 2));

   //--- Session Frame ---------------------------------------------------------
   fSessionFrame = new TSessionFrame(fV2, 350, 310);
   fSessionFrame->Build(this);
   fV2->AddFrame(fSessionFrame, new TGLayoutHints(kLHintsTop | kLHintsExpandX |
         kLHintsExpandY, 2, 0, 1, 2));

   //--- Query Frame -----------------------------------------------------------
   fQueryFrame = new TSessionQueryFrame(fV2, 350, 310);
   fQueryFrame->Build(this);
   fV2->AddFrame(fQueryFrame, new TGLayoutHints(kLHintsTop | kLHintsExpandX |
         kLHintsExpandY, 2, 0, 1, 2));

   //--- Output Frame ----------------------------------------------------------
   fOutputFrame = new TSessionOutputFrame(fV2, 350, 310);
   fOutputFrame->Build(this);
   fV2->AddFrame(fOutputFrame, new TGLayoutHints(kLHintsTop | kLHintsExpandX |
         kLHintsExpandY, 2, 0, 1, 2));

   //--- Input Frame -----------------------------------------------------------
   fInputFrame = new TSessionInputFrame(fV2, 350, 310);
   fInputFrame->Build(this);
   fV2->AddFrame(fInputFrame, new TGLayoutHints(kLHintsTop | kLHintsExpandX |
         kLHintsExpandY, 2, 0, 1, 2));

   fHf->AddFrame(fV1, new TGLayoutHints(kLHintsLeft | kLHintsExpandY));

   // add vertical splitter between list tree and frames
   TGVSplitter *splitter = new TGVSplitter(fHf, 4);
   splitter->SetFrame(fV1, kTRUE);
   fHf->AddFrame(splitter,new TGLayoutHints(kLHintsLeft | kLHintsExpandY));
   fHf->AddFrame(new TGVertical3DLine(fHf), new TGLayoutHints(kLHintsLeft |
         kLHintsExpandY));

   fHf->AddFrame(fV2, new TGLayoutHints(kLHintsRight | kLHintsExpandX |
         kLHintsExpandY));

   AddFrame(fHf, new TGLayoutHints(kLHintsRight | kLHintsExpandX |
         kLHintsExpandY));

   // if description available, update server infos frame
   if (fActDesc) {
      if (!fActDesc->fLocal) {
         fServerFrame->Update(fActDesc);
      }
      else {
         fServerFrame->SetAddEnabled();
         fServerFrame->SetConnectEnabled(kFALSE);
      }
   }

   //--- Status Bar ------------------------------------------------------------
   int parts[] = { 36, 49, 15 };
   fStatusBar = new TGStatusBar(this, 10, 10);
   fStatusBar->SetCleanup(kDeepCleanup);
   fStatusBar->SetParts(parts, 3);
   for (int p = 0; p < 3; ++p)
      fStatusBar->GetBarPart(p)->SetCleanup(kDeepCleanup);
   AddFrame(fStatusBar, new TGLayoutHints(kLHintsTop | kLHintsLeft |
         kLHintsExpandX, 0, 0, 1, 1));

   // connection icon (animation) and time info
   fStatusBar->SetText("      00:00:00", 2);
   TGCompositeFrame *leftpart = fStatusBar->GetBarPart(2);
   fRightIconPicture = (TGPicture *)fClient->GetPicture("proof_disconnected.xpm");
   fRightIcon = new TGIcon(leftpart, fRightIconPicture,
         fRightIconPicture->GetWidth(),fRightIconPicture->GetHeight());
   leftpart->AddFrame(fRightIcon, new TGLayoutHints(kLHintsLeft, 2, 0, 0, 0));

   // connection progress bar
   TGCompositeFrame *rightpart = fStatusBar->GetBarPart(0);
   fConnectProg = new TGHProgressBar(rightpart, TGProgressBar::kStandard, 100);
   fConnectProg->ShowPosition();
   fConnectProg->SetBarColor("green");
   rightpart->AddFrame(fConnectProg, new TGLayoutHints(kLHintsExpandX, 1, 1, 1, 1));

   // add user info
   fUserGroup = gSystem->GetUserInfo();
   buf.Form("User : %s - %s", fUserGroup->fRealName.Data(),
             fUserGroup->fGroup.Data());
   fStatusBar->SetText(buf.Data(), 1);

   fTimer = 0;

   // create context menu
   fContextMenu = new TContextMenu("SessionViewerContextMenu") ;

   SetWindowName("ROOT Session Viewer");
   MapSubwindows();
   MapWindow();

   // hide frames
   fServerFrame->SetAddEnabled(kFALSE);
   fStatusBar->GetBarPart(0)->HideFrame(fConnectProg);
   fV2->HideFrame(fSessionFrame);
   fV2->HideFrame(fQueryFrame);
   fV2->HideFrame(fOutputFrame);
   fV2->HideFrame(fInputFrame);
   fQueryFrame->GetQueryEditFrame()->OnNewQueryMore();
   fActFrame = fServerFrame;
   UpdateListOfProofs();
   Resize(610, 420);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TSessionViewer::~TSessionViewer()
{
   delete fUserGroup;
   if (gSessionViewer == this)
      gSessionViewer = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse clicks in list tree.

void TSessionViewer::OnListTreeClicked(TGListTreeItem *entry, Int_t btn,
                                       Int_t x, Int_t y)
{
   TList *objlist;
   TObject *obj;
   TString msg;

   fSessionMenu->DisableEntry(kSessionAdd);
   if (fToolBar->GetButton(kQuerySubmit))
      fToolBar->GetButton(kQuerySubmit)->SetState(kButtonDisabled);
   if (entry->GetParent() == 0) {  // PROOF
      // switch frames only if actual one doesn't match
      if (fActFrame != fServerFrame) {
         fV2->HideFrame(fActFrame);
         fV2->ShowFrame(fServerFrame);
         fActFrame = fServerFrame;
      }
      fSessionMenu->DisableEntry(kSessionDelete);
      fSessionMenu->EnableEntry(kSessionAdd);
      fServerFrame->SetAddEnabled();
      fServerFrame->SetConnectEnabled(kFALSE);
      fPopupSrv->DisableEntry(kSessionConnect);
      fSessionMenu->DisableEntry(kSessionConnect);
      if (fToolBar->GetButton(kSessionConnect))
         fToolBar->GetButton(kSessionConnect)->SetState(kButtonDisabled);
   }
   else if (entry->GetParent()->GetParent() == 0) { // Server
      if (entry->GetUserData()) {
         obj = (TObject *)entry->GetUserData();
         if (!obj || obj->IsA() != TSessionDescription::Class())
            return;
         // update server frame information
         fServerFrame->Update((TSessionDescription *)obj);
         fActDesc = (TSessionDescription*)obj;
         // if Proof valid, update connection infos
         if (fActDesc->fConnected && fActDesc->fAttached &&
             fActDesc->fProof && fActDesc->fProof->IsValid()) {
            fActDesc->fProof->cd();
            msg.Form("PROOF Cluster %s ready", fActDesc->fName.Data());
         }
         else {
            msg.Form("PROOF Cluster %s not connected", fActDesc->fName.Data());
         }
         fStatusBar->SetText(msg.Data(), 1);
      }
      if ((fActDesc->fConnected) && (fActDesc->fAttached)) {
         fPopupSrv->DisableEntry(kSessionConnect);
         fSessionMenu->DisableEntry(kSessionConnect);
         if (fToolBar->GetButton(kSessionConnect))
            fToolBar->GetButton(kSessionConnect)->SetState(kButtonDisabled);
         UpdateListOfPackages();
         fSessionFrame->UpdateListOfDataSets();
      }
      else {
         fPopupSrv->EnableEntry(kSessionConnect);
         fSessionMenu->EnableEntry(kSessionConnect);
         if (fToolBar->GetButton(kSessionConnect))
            fToolBar->GetButton(kSessionConnect)->SetState(kButtonUp);
      }
      // local session
      if (fActDesc->fLocal) {
         if (fActFrame != fSessionFrame) {
            fV2->HideFrame(fActFrame);
            fV2->ShowFrame(fSessionFrame);
            fActFrame = fSessionFrame;
            UpdateListOfPackages();
            fSessionFrame->UpdateListOfDataSets();
         }
         fSessionFrame->SetLocal();
         fServerFrame->SetAddEnabled();
         fServerFrame->SetConnectEnabled(kFALSE);
      }
      // proof session not connected
      if ((!fActDesc->fLocal) && (!fActDesc->fAttached) &&
          (fActFrame != fServerFrame)) {
         fV2->HideFrame(fActFrame);
         fV2->ShowFrame(fServerFrame);
         fActFrame = fServerFrame;
      }
      // proof session connected
      if ((!fActDesc->fLocal) && (fActDesc->fConnected) &&
          (fActDesc->fAttached)) {
         if (fActFrame != fSessionFrame) {
            fV2->HideFrame(fActFrame);
            fV2->ShowFrame(fSessionFrame);
            fActFrame = fSessionFrame;
         }
         fSessionFrame->SetLocal(kFALSE);
      }
      fSessionFrame->SetLogLevel(fActDesc->fLogLevel);
      fServerFrame->SetLogLevel(fActDesc->fLogLevel);
      if (fActDesc->fAutoEnable)
         fSessionFrame->CheckAutoEnPack(kTRUE);
      else
         fSessionFrame->CheckAutoEnPack(kFALSE);
      // update session information frame
      fSessionFrame->ProofInfos();
      fSessionFrame->UpdatePackages();
      fServerFrame->SetAddEnabled(kFALSE);
      fServerFrame->SetConnectEnabled();
   }
   else if (entry->GetParent()->GetParent()->GetParent() == 0) { // query
      obj = (TObject *)entry->GetParent()->GetUserData();
      if (obj && obj->IsA() == TSessionDescription::Class()) {
         fActDesc = (TSessionDescription *)obj;
      }
      obj = (TObject *)entry->GetUserData();
      if (obj && obj->IsA() == TQueryDescription::Class()) {
         fActDesc->fActQuery = (TQueryDescription *)obj;
      }
      // update query information and buttons state
      fQueryFrame->UpdateInfos();
      fQueryFrame->UpdateButtons(fActDesc->fActQuery);
      if (fActFrame != fQueryFrame) {
         fV2->HideFrame(fActFrame);
         fV2->ShowFrame(fQueryFrame);
         fActFrame = fQueryFrame;
      }
      if ((fActDesc->fConnected) && (fActDesc->fAttached) &&
          (fActDesc->fActQuery->fStatus != TQueryDescription::kSessionQueryRunning) &&
          (fActDesc->fActQuery->fStatus != TQueryDescription::kSessionQuerySubmitted) &&
          (fToolBar->GetButton(kQuerySubmit)) )
         fToolBar->GetButton(kQuerySubmit)->SetState(kButtonUp);
      // trick to update feedback histos
      OnCascadeMenu();
   }
   else {   // a list (input, output)
      obj = (TObject *)entry->GetParent()->GetParent()->GetUserData();
      if (obj && obj->IsA() == TSessionDescription::Class()) {
         fActDesc = (TSessionDescription *)obj;
      }
      obj = (TObject *)entry->GetParent()->GetUserData();
      if (obj && obj->IsA() == TQueryDescription::Class()) {
         fActDesc->fActQuery = (TQueryDescription *)obj;
      }
      if (fActDesc->fActQuery) {
         // update input/output list views
         fInputFrame->RemoveAll();
         fOutputFrame->RemoveAll();
         if (fActDesc->fActQuery->fResult) {
            objlist = fActDesc->fActQuery->fResult->GetOutputList();
            if (objlist) {
               TIter nexto(objlist);
               while ((obj = (TObject *) nexto())) {
                  fOutputFrame->AddObject(obj);
               }
            }
            objlist = fActDesc->fActQuery->fResult->GetInputList();
            if (objlist) {
               TIter nexti(objlist);
               while ((obj = (TObject *) nexti())) {
                  fInputFrame->AddObject(obj);
               }
            }
         }
         else {
            TChain *chain = (TChain *)fActDesc->fActQuery->fChain;
            if (chain && chain->GetPlayer()) {
               objlist = chain->GetPlayer()->GetSelectorFromFile()->GetOutputList();
               if (objlist) {
                  TIter nexto(objlist);
                  while ((obj = (TObject *) nexto())) {
                     fOutputFrame->AddObject(obj);
                  }
               }
            }
         }
         fInputFrame->Resize();
         fOutputFrame->Resize();
         fClient->NeedRedraw(fOutputFrame->GetLVContainer());
         fClient->NeedRedraw(fInputFrame->GetLVContainer());
      }
      // switch frames
      if (strstr(entry->GetText(),"Output")) {
         if (fActFrame != fOutputFrame) {
            fV2->HideFrame(fActFrame);
            fV2->ShowFrame(fOutputFrame);
            fActFrame = fOutputFrame;
         }
      }
      else if (strstr(entry->GetText(),"Input")) {
         if (fActFrame != fInputFrame) {
            fV2->HideFrame(fActFrame);
            fV2->ShowFrame(fInputFrame);
            fActFrame = fInputFrame;
         }
      }
   }
   if (btn == 3) { // right button
      // place popup menus
      TGListTreeItem *item = fSessionHierarchy->GetSelected();
      if (!item) return;
      obj = (TObject *)item->GetUserData();
      if (obj && obj->IsA() == TQueryDescription::Class()) {
         fPopupQry->PlaceMenu(x, y, 1, 1);
      }
      else if (obj && obj->IsA() == TSessionDescription::Class()) {
         if (!fActDesc->fLocal)
            fPopupSrv->PlaceMenu(x, y, 1, 1);
      }
   }
   // enable / disable menu entries
   if (fActDesc->fConnected && fActDesc->fAttached) {
      fSessionMenu->EnableEntry(kSessionGetQueries);
      fSessionMenu->EnableEntry(kSessionShowStatus);
      fPopupSrv->EnableEntry(kSessionGetQueries);
      fPopupSrv->EnableEntry(kSessionShowStatus);
      fPopupSrv->EnableEntry(kSessionDisconnect);
      fPopupSrv->EnableEntry(kSessionShutdown);
      fPopupSrv->EnableEntry(kSessionCleanup);
      fPopupSrv->EnableEntry(kSessionReset);
      fSessionMenu->EnableEntry(kSessionDisconnect);
      fSessionMenu->EnableEntry(kSessionShutdown);
      fSessionMenu->EnableEntry(kSessionCleanup);
      fSessionMenu->EnableEntry(kSessionReset);
      if (fToolBar->GetButton(kSessionDisconnect))
         fToolBar->GetButton(kSessionDisconnect)->SetState(kButtonUp);
      fQueryMenu->EnableEntry(kQuerySubmit);
      fPopupQry->EnableEntry(kQuerySubmit);
   }
   else {
      fSessionMenu->DisableEntry(kSessionGetQueries);
      fSessionMenu->DisableEntry(kSessionShowStatus);
      fPopupSrv->DisableEntry(kSessionGetQueries);
      fPopupSrv->DisableEntry(kSessionShowStatus);
      if (entry->GetParent() != 0)
         fSessionMenu->EnableEntry(kSessionDelete);
      fPopupSrv->EnableEntry(kSessionDelete);
      fPopupSrv->DisableEntry(kSessionDisconnect);
      fPopupSrv->DisableEntry(kSessionShutdown);
      fPopupSrv->DisableEntry(kSessionCleanup);
      fPopupSrv->DisableEntry(kSessionReset);
      fSessionMenu->DisableEntry(kSessionDisconnect);
      fSessionMenu->DisableEntry(kSessionShutdown);
      fSessionMenu->DisableEntry(kSessionCleanup);
      fSessionMenu->DisableEntry(kSessionReset);
      if (fToolBar->GetButton(kSessionDisconnect))
         fToolBar->GetButton(kSessionDisconnect)->SetState(kButtonDisabled);
      fQueryMenu->DisableEntry(kQuerySubmit);
      fPopupQry->DisableEntry(kQuerySubmit);
   }
   if (fActDesc->fLocal) {
      fSessionMenu->DisableEntry(kSessionDelete);
      fSessionMenu->DisableEntry(kSessionConnect);
      fSessionMenu->DisableEntry(kSessionDisconnect);
      fSessionMenu->DisableEntry(kSessionShutdown);
      fSessionMenu->DisableEntry(kSessionCleanup);
      fSessionMenu->DisableEntry(kSessionReset);
      if (fToolBar->GetButton(kSessionDisconnect))
         fToolBar->GetButton(kSessionDisconnect)->SetState(kButtonDisabled);
      if (fToolBar->GetButton(kSessionConnect))
         fToolBar->GetButton(kSessionConnect)->SetState(kButtonDisabled);
      fQueryMenu->EnableEntry(kQuerySubmit);
      fPopupQry->EnableEntry(kQuerySubmit);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse double clicks in list tree (connect to server).

void TSessionViewer::OnListTreeDoubleClicked(TGListTreeItem *entry, Int_t /*btn*/)
{
   if (entry == fSessionItem)
      return;
   if (entry->GetParent()->GetParent() == 0) { // Server
      if (entry->GetUserData()) {
         TObject *obj = (TObject *)entry->GetUserData();
         if ((!obj) || (obj->IsA() != TSessionDescription::Class()))
            return;
         fActDesc = (TSessionDescription*)obj;
         // if Proof valid, update connection infos
      }
      if ((!fActDesc->fLocal) && ((!fActDesc->fConnected) ||
          (!fActDesc->fAttached))) {
         fServerFrame->OnBtnConnectClicked();
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Terminate Session : save configuration, clean temporary files and close
/// Proof connections.

void TSessionViewer::Terminate()
{
   // clean-up temporary files
   TString pathtmp;
   pathtmp = TString::Format("%s/%s", gSystem->TempDirectory(), kSession_RedirectFile);
   if (!gSystem->AccessPathName(pathtmp)) {
      gSystem->Unlink(pathtmp);
   }
   pathtmp = TString::Format("%s/%s", gSystem->TempDirectory(), kSession_RedirectCmd);
   if (!gSystem->AccessPathName(pathtmp)) {
      gSystem->Unlink(pathtmp);
   }
   // close opened Proof sessions (if any)
   TIter next(fSessions);
   TSessionDescription *desc = 0;
   while ((desc = (TSessionDescription *)next())) {
      if (desc->fAttached && desc->fProof &&
          desc->fProof->IsValid())
         desc->fProof->Detach();
   }
   // Save configuration
   if (fAutoSave)
      WriteConfiguration();
}

////////////////////////////////////////////////////////////////////////////////
/// Close main Session Viewer window.

void TSessionViewer::CloseWindow()
{
   // clean-up temporary files
   TString pathtmp;
   pathtmp = TString::Format("%s/%s", gSystem->TempDirectory(), kSession_RedirectFile);
   if (!gSystem->AccessPathName(pathtmp)) {
      gSystem->Unlink(pathtmp);
   }
   pathtmp = TString::Format("%s/%s", gSystem->TempDirectory(), kSession_RedirectCmd);
   if (!gSystem->AccessPathName(pathtmp)) {
      gSystem->Unlink(pathtmp);
   }
   // Save configuration
   if (fAutoSave)
      WriteConfiguration();
   Cleanup();
   fSessions->Delete();
   if (fSessionItem)
      fSessionHierarchy->DeleteChildren(fSessionItem);
   delete fSessionHierarchy; // this has been put int TGCanvas which isn't a
                             // TGComposite frame and doesn't do cleanups.
   fClient->FreePicture(fLocal);
   fClient->FreePicture(fProofCon);
   fClient->FreePicture(fProofDiscon);
   fClient->FreePicture(fQueryCon);
   fClient->FreePicture(fQueryDiscon);
   fClient->FreePicture(fBaseIcon);
   delete fTimer;
   DeleteWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Change the right logo (used for animation).

void TSessionViewer::ChangeRightLogo(const char *name)
{
   fClient->FreePicture(fRightIconPicture);
   fRightIconPicture = (TGPicture *)fClient->GetPicture(name);
   fRightIcon->SetPicture(fRightIconPicture);
}

////////////////////////////////////////////////////////////////////////////////
/// Enable animation timer.

void TSessionViewer::EnableTimer()
{
   if (!fTimer) fTimer = new TTimer(this, 500);
   fTimer->Reset();
   fTimer->TurnOn();
   time( &fStart );
}

////////////////////////////////////////////////////////////////////////////////
/// Disable animation timer.

void TSessionViewer::DisableTimer()
{
   if (fTimer)
      fTimer->TurnOff();
   ChangeRightLogo("proof_disconnected.xpm");
}

////////////////////////////////////////////////////////////////////////////////
/// Handle animation timer.

Bool_t TSessionViewer::HandleTimer(TTimer *)
{
   TString buf;
   struct tm *connected;
   Int_t count = gRandom->Integer(4);
   if (count > 3) {
      count = 0;
   }
   if (fChangePic)
      ChangeRightLogo(xpm_names[count]);
   time( &fElapsed );
   time_t elapsed_time = (time_t)difftime( fElapsed, fStart );
   connected = gmtime( &elapsed_time );
   if (connected) {
      buf.Form("      %02d:%02d:%02d", connected->tm_hour,
               connected->tm_min, connected->tm_sec);
      fStatusBar->SetText(buf.Data(), 2);
   }
   else {
      fStatusBar->SetText("      00:00:00", 2);
   }

   if (fActDesc->fLocal) {
      if ((fActDesc->fActQuery) &&
         (fActDesc->fActQuery->fStatus ==
         TQueryDescription::kSessionQueryRunning)) {
         TChain *chain = (TChain *)fActDesc->fActQuery->fChain;
         if (chain)
            fQueryFrame->ProgressLocal(chain->GetEntries(),
                                       chain->GetReadEntry()+1);
      }
   }

   fTimer->Reset();
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Load/append a log msg in the log frame.

void TSessionViewer::LogMessage(const char *msg, Bool_t all)
{
   if (fLogWindow) {
      if (all) {
         // load buffer
         fLogWindow->LoadBuffer(msg);
      } else {
         // append
         fLogWindow->AddBuffer(msg);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Handle signal "query result ready" coming from Proof session.

void TSessionViewer::QueryResultReady(char *query)
{
   TString strtmp;
   strtmp.Form("Query Result Ready for %s", query);
   // show information on status bar
   ShowInfo(strtmp.Data());
   TGListTreeItem *item=0, *item2=0;
   TQueryDescription *lquery = 0;
   // loop over actual queries to find which one is ready

   TIter nexts(fSessions);
   TSessionDescription *desc = 0;
   // check if session is already in the list
   while ((desc = (TSessionDescription *)nexts())) {
      if (desc && !desc->fAttached)
         continue;
      TIter nextp(desc->fQueries);
      while ((lquery = (TQueryDescription *)nextp())) {
         if (lquery->fReference.Contains(query)) {
            // results are ready for this query
            lquery->fResult = desc->fProof->GetQueryResult(query);
            lquery->fStatus = TQueryDescription::kSessionQueryFromProof;
            if (!lquery->fResult)
               break;
            // get query status
            lquery->fStatus = lquery->fResult->IsFinalized() ?
               TQueryDescription::kSessionQueryFinalized :
               (TQueryDescription::ESessionQueryStatus)lquery->fResult->GetStatus();
            // get data set
            TObject *o = lquery->fResult->GetInputObject("TDSet");
            if (o)
               lquery->fChain = (TDSet *) o;
            item = fSessionHierarchy->FindItemByObj(fSessionItem, desc);
            if (item) {
               item2 = fSessionHierarchy->FindItemByObj(item, lquery);
            }
            if (item2) {
               // add input and output list entries
               if (lquery->fResult->GetInputList())
                  if (!fSessionHierarchy->FindChildByName(item2, "InputList"))
                     fSessionHierarchy->AddItem(item2, "InputList");
               if (lquery->fResult->GetOutputList())
                  if (!fSessionHierarchy->FindChildByName(item2, "OutputList"))
                     fSessionHierarchy->AddItem(item2, "OutputList");
            }
            // update list tree, query frame information, and buttons state
            fSessionHierarchy->ClearViewPort();
            fClient->NeedRedraw(fSessionHierarchy);
            fQueryFrame->UpdateInfos();
            fQueryFrame->UpdateButtons(lquery);
            break;
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Clean-up Proof session.

void TSessionViewer::CleanupSession()
{
   TGListTreeItem *item = fSessionHierarchy->GetSelected();
   if (!item) return;
   TObject *obj = (TObject *)item->GetUserData();
   if (!obj || obj->IsA() != TSessionDescription::Class()) return;
   if (!fActDesc->fProof || !fActDesc->fProof->IsValid()) return;
   TString m;
   m.Form("Are you sure to cleanup the session \"%s::%s\"",
         fActDesc->fName.Data(), fActDesc->fTag.Data());
   Int_t result;
   new TGMsgBox(fClient->GetRoot(), this, "", m.Data(), 0,
         kMBYes | kMBNo | kMBCancel, &result);
   if (result == kMBYes) {
      // send cleanup request for the session specified by the tag reference
      TString sessiontag;
      sessiontag.Form("session-%s",fActDesc->fTag.Data());
      fActDesc->fProof->CleanupSession(sessiontag.Data());
      // clear the list of queries
      fActDesc->fQueries->Clear();
      fSessionHierarchy->DeleteChildren(item);
      fSessionFrame->OnBtnGetQueriesClicked();
      if (fAutoSave)
         WriteConfiguration();
   }
   // update list tree
   fSessionHierarchy->ClearViewPort();
   fClient->NeedRedraw(fSessionHierarchy);
}

////////////////////////////////////////////////////////////////////////////////
/// Reset Proof session.

void TSessionViewer::ResetSession()
{
   TGListTreeItem *item = fSessionHierarchy->GetSelected();
   if (!item) return;
   TObject *obj = (TObject *)item->GetUserData();
   if (!obj || obj->IsA() != TSessionDescription::Class()) return;
   if (!fActDesc->fProof || !fActDesc->fProof->IsValid()) return;
   TString m;
   m.Form("Do you really want to reset the session \"%s::%s\"",
         fActDesc->fName.Data(), fActDesc->fAddress.Data());
   Int_t result;
   new TGMsgBox(fClient->GetRoot(), this, "", m.Data(), 0,
         kMBYes | kMBNo | kMBCancel, &result);
   if (result == kMBYes) {
      // reset the session
      TProofMgr *mgr = TProof::Mgr(fActDesc->fAddress);
      if (mgr && mgr->IsValid()) {
         mgr->Reset(fActDesc->fUserName);
      }
      // reset connected flag
      fActDesc->fAttached = kFALSE;
      fActDesc->fProof = 0;
      // disable animation timer
      DisableTimer();
      // change list tree item picture to disconnected pixmap
      TGListTreeItem *item2 = fSessionHierarchy->FindChildByData(
                              fSessionItem, fActDesc);
      if (item2) item2->SetPictures(fProofDiscon, fProofDiscon);

      OnListTreeClicked(fSessionHierarchy->GetSelected(), 1, 0, 0);
      fSessionHierarchy->ClearViewPort();
      fClient->NeedRedraw(fSessionHierarchy);
      fStatusBar->SetText("", 1);
   }
   // update list tree
   fSessionHierarchy->ClearViewPort();
   fClient->NeedRedraw(fSessionHierarchy);
}

////////////////////////////////////////////////////////////////////////////////
/// Delete query from list tree and ask user if they want to delete it also
/// from server.

void TSessionViewer::DeleteQuery()
{
   TGListTreeItem *item = fSessionHierarchy->GetSelected();
   if (!item) return;
   TObject *obj = (TObject *)item->GetUserData();
   if (!obj || obj->IsA() != TQueryDescription::Class()) return;
   TQueryDescription *query = (TQueryDescription *)obj;
   TString m;
   Int_t result = 0;

   if (fActDesc->fAttached && fActDesc->fProof && fActDesc->fProof->IsValid()) {
      if ((fActDesc->fActQuery->fStatus == TQueryDescription::kSessionQuerySubmitted) ||
          (fActDesc->fActQuery->fStatus == TQueryDescription::kSessionQueryRunning) ) {
         new TGMsgBox(fClient->GetRoot(), this, "Delete Query",
                      "Deleting running queries is not allowed", kMBIconExclamation,
                      kMBOk, &result);
         return;
      }
      m.Form("Do you want to delete query \"%s\" from server too ?",
            query->fQueryName.Data());
      new TGMsgBox(fClient->GetRoot(), this, "", m.Data(), kMBIconQuestion,
            kMBYes | kMBNo | kMBCancel, &result);
   }
   else {
      m.Form("Dou you really want to delete query \"%s\" ?",
            query->fQueryName.Data());
      new TGMsgBox(fClient->GetRoot(), this, "", m.Data(), kMBIconQuestion,
            kMBOk | kMBCancel, &result);
   }
   if (result == kMBYes) {
      fActDesc->fProof->Remove(query->fReference.Data());
      fActDesc->fQueries->Remove((TObject *)query);
      fSessionHierarchy->DeleteItem(item);
      delete query;
   }
   else if (result == kMBNo || result == kMBOk) {
      fActDesc->fQueries->Remove((TObject *)query);
      fSessionHierarchy->DeleteItem(item);
      delete query;
   }
   fSessionHierarchy->ClearViewPort();
   fClient->NeedRedraw(fSessionHierarchy);
   if (fAutoSave)
      WriteConfiguration();
}

////////////////////////////////////////////////////////////////////////////////
/// Edit currently selected query.

void TSessionViewer::EditQuery()
{
   TGListTreeItem *item = fSessionHierarchy->GetSelected();
   if (!item) return;
   TObject *obj = (TObject *)item->GetUserData();
   if (!obj || obj->IsA() != TQueryDescription::Class()) return;
   TQueryDescription *query = (TQueryDescription *)obj;
   TNewQueryDlg *dlg = new TNewQueryDlg(this, 350, 310, query, kTRUE);
   dlg->Popup();
}

////////////////////////////////////////////////////////////////////////////////
/// Start TreeViewer from selected TChain.

void TSessionViewer::StartViewer()
{
   TGListTreeItem *item = fSessionHierarchy->GetSelected();
   if (!item) return;
   TObject *obj = (TObject *)item->GetUserData();
   if (!obj || obj->IsA() != TQueryDescription::Class()) return;
   TQueryDescription *query = (TQueryDescription *)obj;
   if (!query->fChain && query->fResult &&
      (obj = query->fResult->GetInputObject("TDSet"))) {
      query->fChain = (TDSet *) obj;
   }
   if (!query->fChain) return;
   if (query->fChain->IsA() == TChain::Class())
      ((TChain *)query->fChain)->StartViewer();
   else if (query->fChain->IsA() == TDSet::Class())
      ((TDSet *)query->fChain)->StartViewer();
}

////////////////////////////////////////////////////////////////////////////////
/// Query the list of uploaded packages from proof and display it
/// into a new text window.

void TSessionViewer::ShowPackages()
{
   Window_t wdummy;
   Int_t  ax, ay;

   if (fActDesc->fLocal) return;
   if (!fActDesc->fProof || !fActDesc->fProof->IsValid())
      return;
   TString pathtmp = TString::Format("%s/%s", gSystem->TempDirectory(),
            kSession_RedirectFile);
   // redirect stdout/stderr to temp file
   if (gSystem->RedirectOutput(pathtmp.Data(), "w") != 0) {
      Error("ShowStatus", "stdout/stderr redirection failed; skipping");
      return;
   }
   fActDesc->fProof->ShowPackages(kTRUE);
   // restore stdout/stderr
   if (gSystem->RedirectOutput(0) != 0) {
      Error("ShowStatus", "stdout/stderr retore failed; skipping");
      return;
   }
   if (!fLogWindow) {
      fLogWindow = new TSessionLogView(this, 700, 100);
   } else {
      // Clear window
      fLogWindow->Clear();
   }
   fLogWindow->LoadFile(pathtmp.Data());
   gVirtualX->TranslateCoordinates(GetId(), fClient->GetDefaultRoot()->GetId(),
                                    0, 0, ax, ay, wdummy);
   fLogWindow->Move(ax, ay + GetHeight() + 35);
   fLogWindow->Popup();
}

////////////////////////////////////////////////////////////////////////////////
/// Update the list of packages.

void TSessionViewer::UpdateListOfPackages()
{
   TObjString *packname;
   TPackageDescription *package;
   if (fActDesc->fConnected && fActDesc->fAttached &&
       fActDesc->fProof && fActDesc->fProof->IsValid() &&
       fActDesc->fProof->IsParallel()) {
      //fActDesc->fPackages->Clear();
      TList *packlist = fActDesc->fProof->GetListOfEnabledPackages();
      if(packlist) {
         TIter nextenabled(packlist);
         while ((packname = (TObjString *)nextenabled())) {
            package = new TPackageDescription;
            package->fName = packname->GetName();
            package->fName += ".par";
            package->fPathName = package->fName;
            package->fId   = fActDesc->fPackages->GetEntries();
            package->fUploaded = kTRUE;
            package->fEnabled = kTRUE;
            if (!fActDesc->fPackages->FindObject(package->fName)) {
               fActDesc->fPackages->Add((TObject *)package);
            }
         }
      }
      packlist = fActDesc->fProof->GetListOfPackages();
      if(packlist) {
         TIter nextpack(packlist);
         while ((packname = (TObjString *)nextpack())) {
            package = new TPackageDescription;
            package->fName = packname->GetName();
            package->fName += ".par";
            package->fPathName = package->fName;
            package->fId   = fActDesc->fPackages->GetEntries();
            package->fUploaded = kTRUE;
            package->fEnabled = kFALSE;
            if (!fActDesc->fPackages->FindObject(package->fName)) {
               fActDesc->fPackages->Add((TObject *)package);
            }
         }
      }
   }
//   fSessionFrame->UpdatePackages();
}

////////////////////////////////////////////////////////////////////////////////
/// Query list of enabled packages from proof and display it
/// into a new text window.

void TSessionViewer::ShowEnabledPackages()
{
   Window_t wdummy;
   Int_t  ax, ay;

   if (fActDesc->fLocal) return;
   if (!fActDesc->fProof || !fActDesc->fProof->IsValid())
      return;
   TString pathtmp = TString::Format("%s/%s", gSystem->TempDirectory(),
         kSession_RedirectFile);
   // redirect stdout/stderr to temp file
   if (gSystem->RedirectOutput(pathtmp.Data(), "w") != 0) {
      Error("ShowStatus", "stdout/stderr redirection failed; skipping");
      return;
   }
   fActDesc->fProof->ShowEnabledPackages(kTRUE);
   // restore stdout/stderr
   if (gSystem->RedirectOutput(0) != 0) {
      Error("ShowStatus", "stdout/stderr retore failed; skipping");
      return;
   }
   if (!fLogWindow) {
      fLogWindow = new TSessionLogView(this, 700, 100);
   } else {
      // Clear window
      fLogWindow->Clear();
   }
   fLogWindow->LoadFile(pathtmp.Data());
   gVirtualX->TranslateCoordinates(GetId(), fClient->GetDefaultRoot()->GetId(),
                                    0, 0, ax, ay, wdummy);
   fLogWindow->Move(ax, ay + GetHeight() + 35);
   fLogWindow->Popup();
}

////////////////////////////////////////////////////////////////////////////////
/// Display the content of the temporary log file for queryref
/// into a new text window.

void TSessionViewer::ShowLog(const char *queryref)
{
   Window_t wdummy;
   Int_t  ax, ay;

   if (fActDesc->fProof) {
      gVirtualX->SetCursor(GetId(),gVirtualX->CreateCursor(kWatch));
      if (!fLogWindow) {
         fLogWindow = new TSessionLogView(this, 700, 100);
      } else {
         // Clear window
         fLogWindow->Clear();
      }
      fActDesc->fProof->Connect("LogMessage(const char*,Bool_t)",
            "TSessionViewer", this, "LogMessage(const char*,Bool_t)");
      Bool_t logonly = fActDesc->fProof->SendingLogToWindow();
      fActDesc->fProof->SendLogToWindow(kTRUE);
      if (queryref)
         fActDesc->fProof->ShowLog(queryref);
      else
         fActDesc->fProof->ShowLog(0);
      fActDesc->fProof->SendLogToWindow(logonly);
      // set log window position at the bottom of Session Viewer
      gVirtualX->TranslateCoordinates(GetId(),
            fClient->GetDefaultRoot()->GetId(), 0, 0, ax, ay, wdummy);
      fLogWindow->Move(ax, ay + GetHeight() + 35);
      fLogWindow->Popup();
      gVirtualX->SetCursor(GetId(), 0);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Display text in status bar.

void TSessionViewer::ShowInfo(const char *txt)
{
   fStatusBar->SetText(txt,0);
   fClient->NeedRedraw(fStatusBar);
   gSystem->ProcessEvents();
}

////////////////////////////////////////////////////////////////////////////////
/// Retrieve and display Proof status.

void TSessionViewer::ShowStatus()
{
   Window_t wdummy;
   Int_t  ax, ay;

   if (!fActDesc->fProof || !fActDesc->fProof->IsValid())
      return;
   TString pathtmp = TString::Format("%s/%s", gSystem->TempDirectory(),
            kSession_RedirectFile);
   // redirect stdout/stderr to temp file
   if (gSystem->RedirectOutput(pathtmp.Data(), "w") != 0) {
      Error("ShowStatus", "stdout/stderr redirection failed; skipping");
      return;
   }
   fActDesc->fProof->GetStatus();
   // restore stdout/stderr
   if (gSystem->RedirectOutput(0) != 0) {
      Error("ShowStatus", "stdout/stderr retore failed; skipping");
      return;
   }
   if (!fLogWindow) {
      fLogWindow = new TSessionLogView(this, 700, 100);
   } else {
      // Clear window
      fLogWindow->Clear();
   }
   fLogWindow->LoadFile(pathtmp.Data());
   gVirtualX->TranslateCoordinates(GetId(), fClient->GetDefaultRoot()->GetId(),
                                    0, 0, ax, ay, wdummy);
   fLogWindow->Move(ax, ay + GetHeight() + 35);
   fLogWindow->Popup();
}

////////////////////////////////////////////////////////////////////////////////
/// Handle startup message (connection progress) coming from Proof session.

void TSessionViewer::StartupMessage(char *msg, Bool_t, Int_t done, Int_t total)
{
   Float_t pos = Float_t(Double_t(done * 100)/Double_t(total));
   fConnectProg->SetPosition(pos);
   fStatusBar->SetText(msg, 1);
}

////////////////////////////////////////////////////////////////////////////////
/// Handle session viewer custom popup menus.

void TSessionViewer::MyHandleMenu(Int_t id)
{
   switch (id) {

      case kSessionDelete:
         fServerFrame->OnBtnDeleteClicked();
         break;

      case kSessionConnect:
         fServerFrame->OnBtnConnectClicked();
         break;

      case kSessionDisconnect:
         fSessionFrame->OnBtnDisconnectClicked();
         break;

      case kSessionShutdown:
         fSessionFrame->ShutdownSession();
         break;

      case kSessionCleanup:
         CleanupSession();
         break;

      case kSessionReset:
         ResetSession();
         break;

      case kSessionBrowse:
         if (fActDesc->fProof && fActDesc->fProof->IsValid()) {
            TBrowser *b = new TBrowser();
            fActDesc->fProof->Browse(b);
         }
         break;

      case kSessionShowStatus:
         ShowStatus();
         break;

      case kSessionGetQueries:
         fSessionFrame->OnBtnGetQueriesClicked();
         break;

      case kQueryEdit:
         EditQuery();
         break;

      case kQueryDelete:
         DeleteQuery();
         break;

      case kQueryStartViewer:
         StartViewer();
         break;

      case kQuerySubmit:
         fQueryFrame->OnBtnSubmit();
         break;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Handle feedback histograms configuration menu.

void TSessionViewer::OnCascadeMenu()
{
   // divide stats canvas by number of selected feedback histos
   fQueryFrame->GetStatsCanvas()->cd();
   fQueryFrame->GetStatsCanvas()->Clear();
   fQueryFrame->GetStatsCanvas()->Modified();
   fQueryFrame->GetStatsCanvas()->Update();
   if (!fActDesc || !fActDesc->fActQuery) return;
   fActDesc->fNbHistos = 0;
   Int_t i = 0;

   if (fActDesc->fAttached && fActDesc->fProof &&
       fActDesc->fProof->IsValid()) {
      if (fOptionsMenu->IsEntryChecked(kOptionsFeedback)) {
         // browse list of feedback histos and check user's selected ones
         while (kFeedbackHistos[i]) {
            if (fCascadeMenu->IsEntryChecked(41+i)) {
               fActDesc->fProof->AddFeedback(kFeedbackHistos[i]);
            }
            i++;
         }
      }
      else {
         // if feedback option not selected, clear Proof's feedback option
         fActDesc->fProof->ClearFeedback();
      }
   }

   i = 0;
   // loop over feedback histo list
   while (kFeedbackHistos[i]) {
      // check if user has selected this histogram in the option menu
      if (fCascadeMenu->IsEntryChecked(41+i))
         fActDesc->fNbHistos++;
      i++;
   }
   fQueryFrame->GetStatsCanvas()->SetEditable(kTRUE);
   fQueryFrame->GetStatsCanvas()->Clear();
   if (fActDesc->fNbHistos == 4)
      fQueryFrame->GetStatsCanvas()->Divide(2, 2);
   else if (fActDesc->fNbHistos > 4)
      fQueryFrame->GetStatsCanvas()->Divide(3, 2);
   else
      fQueryFrame->GetStatsCanvas()->Divide(fActDesc->fNbHistos, 1);

   // if actual query has results, update feedback histos
   if (fActDesc->fActQuery && fActDesc->fActQuery->fResult &&
       fActDesc->fActQuery->fResult->GetOutputList()) {
      fQueryFrame->UpdateHistos(fActDesc->fActQuery->fResult->GetOutputList());
      fQueryFrame->ResetProgressDialog("", 0, 0, 0);
   }
   else if (fActDesc->fActQuery) {
      fQueryFrame->ResetProgressDialog(fActDesc->fActQuery->fSelectorString,
                                       fActDesc->fActQuery->fNbFiles,
                                       fActDesc->fActQuery->fFirstEntry,
                                       fActDesc->fActQuery->fNoEntries);
   }
   fQueryFrame->UpdateInfos();
}
////////////////////////////////////////////////////////////////////////////////
/// Handle messages send to the TSessionViewer object. E.g. all menu entries
/// messages.

Bool_t TSessionViewer::ProcessMessage(Long_t msg, Long_t parm1, Long_t)
{
   TNewQueryDlg *dlg;

   switch (GET_MSG(msg)) {
      case kC_COMMAND:
         switch (GET_SUBMSG(msg)) {
            case kCM_BUTTON:
            case kCM_MENU:
               switch (parm1) {

                  case kFileCloseViewer:
                     CloseWindow();
                     break;

                  case kFileLoadConfig:
                     {
                        TGFileInfo fi;
                        fi.SetFilename(gSystem->BaseName(fConfigFile));
                        fi.SetIniDir(gSystem->HomeDirectory());
                        fi.fFileTypes = conftypes;
                        new TGFileDialog(fClient->GetRoot(), this, kFDOpen, &fi);
                        if (fi.fFilename) {
                           fConfigFile = fi.fFilename;
                           ReadConfiguration(fConfigFile);
                           OnListTreeClicked(fSessionHierarchy->GetSelected(), 1, 0, 0);
                        }
                     }
                     break;

                  case kFileSaveConfig:
                     {
                        TGFileInfo fi;
                        fi.SetFilename(gSystem->BaseName(fConfigFile));
                        fi.SetIniDir(gSystem->HomeDirectory());
                        fi.fFileTypes = conftypes;
                        new TGFileDialog(fClient->GetRoot(), this, kFDSave, &fi);
                        if (fi.fFilename) {
                           fConfigFile = fi.fFilename;
                           WriteConfiguration(fConfigFile);
                        }
                     }
                     break;

                  case kFileQuit:
                     Terminate();
                     if (!gApplication->ReturnFromRun())
                        delete this;
                     gApplication->Terminate(0);
                     break;

                  case kSessionNew:
                     fServerFrame->OnBtnNewServerClicked();
                     break;

                  case kSessionAdd:
                     fServerFrame->OnBtnAddClicked();
                     break;

                  case kSessionDelete:
                     fServerFrame->OnBtnDeleteClicked();
                     break;

                  case kSessionCleanup:
                     CleanupSession();
                     break;

                  case kSessionReset:
                     ResetSession();
                     break;

                  case kSessionConnect:
                     fServerFrame->OnBtnConnectClicked();
                     break;

                  case kSessionDisconnect:
                     fSessionFrame->OnBtnDisconnectClicked();
                     break;

                  case kSessionShutdown:
                     fSessionFrame->ShutdownSession();
                     break;

                  case kSessionShowStatus:
                     ShowStatus();
                     break;

                  case kSessionGetQueries:
                     fSessionFrame->OnBtnGetQueriesClicked();
                     break;

                  case kQueryNew:
                     dlg = new TNewQueryDlg(this, 350, 310);
                     dlg->Popup();
                     break;

                  case kQueryEdit:
                     EditQuery();
                     break;

                  case kQueryDelete:
                     DeleteQuery();
                     break;

                  case kQueryStartViewer:
                     StartViewer();
                     break;

                  case kQuerySubmit:
                     fQueryFrame->OnBtnSubmit();
                     break;

                  case kOptionsAutoSave:
                     if(fOptionsMenu->IsEntryChecked(kOptionsAutoSave)) {
                        fOptionsMenu->UnCheckEntry(kOptionsAutoSave);
                        fAutoSave = kFALSE;
                     }
                     else {
                        fOptionsMenu->CheckEntry(kOptionsAutoSave);
                        fAutoSave = kTRUE;
                     }
                     break;

                  case kOptionsStatsHist:
                     if(fOptionsMenu->IsEntryChecked(kOptionsStatsHist)) {
                        fOptionsMenu->UnCheckEntry(kOptionsStatsHist);
                        gEnv->SetValue("Proof.StatsHist", 0);
                     }
                     else {
                        fOptionsMenu->CheckEntry(kOptionsStatsHist);
                        gEnv->SetValue("Proof.StatsHist", 1);
                     }
                     break;

                  case kOptionsStatsTrace:
                     if(fOptionsMenu->IsEntryChecked(kOptionsStatsTrace)) {
                        fOptionsMenu->UnCheckEntry(kOptionsStatsTrace);
                        gEnv->SetValue("Proof.StatsTrace", 0);
                     }
                     else {
                        fOptionsMenu->CheckEntry(kOptionsStatsTrace);
                        gEnv->SetValue("Proof.StatsTrace", 1);
                     }
                     break;

                  case kOptionsSlaveStatsTrace:
                     if(fOptionsMenu->IsEntryChecked(kOptionsSlaveStatsTrace)) {
                        fOptionsMenu->UnCheckEntry(kOptionsSlaveStatsTrace);
                        gEnv->SetValue("Proof.SlaveStatsTrace", 0);
                     }
                     else {
                        fOptionsMenu->CheckEntry(kOptionsSlaveStatsTrace);
                        gEnv->SetValue("Proof.SlaveStatsTrace", 1);
                     }
                     break;

                  case kOptionsFeedback:
                     if(fOptionsMenu->IsEntryChecked(kOptionsFeedback)) {
                        fOptionsMenu->UnCheckEntry(kOptionsFeedback);
                     }
                     else {
                        fOptionsMenu->CheckEntry(kOptionsFeedback);
                     }
                     break;

                  case 41:
                  case 42:
                  case 43:
                  case 44:
                  case 45:
                  case 46:
                     if (fCascadeMenu->IsEntryChecked(parm1)) {
                        fCascadeMenu->UnCheckEntry(parm1);
                     }
                     else {
                        fCascadeMenu->CheckEntry(parm1);
                     }
                     OnCascadeMenu();
                     break;

                  case 50:
                     if (fCascadeMenu->IsEntryChecked(parm1)) {
                        fCascadeMenu->UnCheckEntry(parm1);
                     }
                     else {
                        fCascadeMenu->CheckEntry(parm1);
                     }
                     OnCascadeMenu();
                     break;

                  case kHelpAbout:
                     {
#ifdef R__UNIX
		         TString rootx = TROOT::GetBinDir() + "/root -a &";
                        gSystem->Exec(rootx);
#else
#ifdef WIN32
                        new TWin32SplashThread(kTRUE);
#else
                        char str[32];
                        sprintf(str, "About ROOT %s...", gROOT->GetVersion());
                        TRootHelpDialog *hd = new TRootHelpDialog(this, str, 600, 400);
                        hd->SetText(gHelpAbout);
                        hd->Popup();
#endif
#endif
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
