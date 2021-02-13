// @(#)root/gui:$Id$
// Author: Katerina Opocenska   11/09/2008

/*************************************************************************
* Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
* All rights reserved.                                                  *
*                                                                       *
* For the licensing terms see $ROOTSYS/LICENSE.                         *
* For the list of contributors see $ROOTSYS/README/CREDITS.             *
*************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  ROOT EVENT RECORDING SYSTEM                                         //
// ==================================================================   //
//                                                                      //
//  TRecorder class provides interface for recording and replaying      //
//  events in ROOT.                                                     //
//  Recorded events are:                                                //
//  - Commands typed by user in commandline ('new TCanvas')             //
//  - GUI events (mouse movement, button clicks, ...)                   //
//                                                                      //
//  All the recorded events from one session are stored in one TFile    //
//  and can be replayed again anytime.                                  //
//                                                                      //
//  Recording                                                           //
//  ==================================================================  //
//                                                                      //
//  1] To start recording                                               //
//                                                                      //
//    TRecorder r(const char *filename, "NEW")                          //
//    TRecorder r(const char *filename, "RECREATE")                     //
//                                                                      //
//    or:                                                               //
//                                                                      //
//    TRecorder *recorder = new TRecorder;                              //
//    recorder->Start(const char *filename, ...)                        //
//                                                                      //
//    -filename      Name of ROOT file in which to save                 //
//                   recorded events.                                   //
//                                                                      //
//  2] To stop recording                                                //
//                                                                      //
//    recorder->Stop()                                                  //
//                                                                      //
//                                                                      //
//  IMPORTANT:                                                          //
//  State capturing is part of recording. It means that if you want to  //
//  record events for some object (window), creation of this object     //
//  must be also recorded.                                              //
//                                                                      //
//    Example:                                                          //
//    --------                                                          //
//    t = new TRecorder();          // Create a new recorder            //
//    t->Start("logfile.root");     // ! Start recording first          //
//                                                                      //
//    c = new TCanvas();            // ! Then, create an object         //
//    c->Dump();                    // Work with that object            //
//                                                                      //
//    t->Stop();                    // Stop recording                   //
//                                                                      //
//  It is strongly recommended to start recording with empty ROOT       //
//  environment, at least with no previously created ROOT GUI.          //
//  This ensures that only events for well known windows are stored.    //
//  Events for windows, which were not created during recording,        //
//  cannot be replayed.                                                 //
//                                                                      //
//  Replaying                                                           //
//  =================================================================== //
//                                                                      //
//  1] To start replaying                                               //
//                                                                      //
//    TRecorder r(const char *filename)                                 //
//    TRecorder r(const char *filename, "READ")                         //
//                                                                      //
//    or:                                                               //
//                                                                      //
//    TRecorder *recorder = new TRecorder;                              //
//    recorder->Replay(const char *filename,                            //
//                      Bool_t showMouseCursor = kTRUE);                //
//                                                                      //
//    -filename         A name of file with recorded events             //
//                      previously created with TRecorder::Start        //
//                                                                      //
//    -showMouseCursor  If kTRUE, mouse cursor is replayed as well.     //
//                      In that case it is not recommended to use mouse //
//                      during replaying.                               //
//                                                                      //
//  In general, it is not recommended to use mouse to change positions  //
//  and states of ROOT windows during replaying.                        //
//                                                                      //
//  IMPORTANT:                                                          //
//  The state of ROOT environment before replaying of some events       //
//  must be exactly the same as before recording them.                  //
//  Therefore it is strongly recommended to start both recording        //
//  and replaying with empty ROOT environment.                          //
//                                                                      //
//  2] To pause replaying                                               //
//                                                                      //
//    recorder->Pause()                                                 //
//                                                                      //
//    Replaying is stopped until recorder->Resume() is called.          //
//                                                                      //
//                                                                      //
//  3] To resume paused replaying                                       //
//                                                                      //
//    recorder->Resume()                                                //
//                                                                      //
//    Resumes previously stopped replaying.                             //
//                                                                      //
//                                                                      //
//  4] To stop replaying before its end                                 //
//                                                                      //
//    recorder->Stop()                                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TRecorder.h"

#include "TROOT.h"
#include "TFile.h"
#include "TTimer.h"
#include "TTree.h"
#include "TSystem.h"
#include "TMutex.h"
#include "TGButton.h"
#include "TGFileDialog.h"
#include "TGLabel.h"
#include "TGWindow.h"
#include "Buttons.h"
#include "TKey.h"
#include "TPaveLabel.h"
#include "TLatex.h"
#include "TVirtualDragManager.h"
#include "TGPicture.h"
#include "KeySymbols.h"
#include "TVirtualX.h"

#include <iomanip>

// Names of ROOT GUI events. Used for listing event logs.
const char *kRecEventNames[] = {
   "KeyPress",
   "KeyRelease",
   "ButtonPress",
   "ButtonRelease",
   "MotionNotify",
   "EnterNotify",
   "LeaveNotify",
   "FocusIn",
   "FocusOut",
   "Expose",
   "ConfigureNotify",
   "MapNotify",
   "UnmapNotify",
   "DestroyNotify",
   "ClientMessage",
   "SelectionClear",
   "SelectionRequest",
   "SelectionNotify",
   "ColormapNotify",
   "ButtonDoubleClick",
   "OtherEvent"
};

// Names of TTrees in the TFile with recorded events
const char *kCmdEventTree   = "CmdEvents";   // Name of TTree with commandline events
const char *kGuiEventTree   = "GuiEvents";   // Name of TTree with GUI events
const char *kWindowsTree    = "WindowsTree"; // Name of TTree with window IDs
const char *kExtraEventTree = "ExtraEvents"; // Name of TTree with extra events (PaveLabels and Texts)
const char *kBranchName     = "MainBranch";  // Name of the main branch in all TTress

ClassImp(TRecorder);


//_____________________________________________________________________________
//
// TGCursorWindow
//
// Window used as fake mouse cursor wile replaying events.
//_____________________________________________________________________________

class TGCursorWindow : public TGFrame {

protected:
   Pixmap_t fPic, fMask;            // Pixmaps used as Window shape

public:
   TGCursorWindow();
   virtual ~TGCursorWindow();
};

static TGCursorWindow *gCursorWin = 0;
static Int_t gDecorWidth  = 0;
static Int_t gDecorHeight = 0;

////////////////////////////////////////////////////////////////////////////////
/// TGCursorWindow constructor.

TGCursorWindow::TGCursorWindow() :
      TGFrame(gClient->GetDefaultRoot(), 32, 32, kTempFrame)
{
   SetWindowAttributes_t wattr;
   const TGPicture *pbg = fClient->GetPicture("recursor.png");
   fPic  = pbg->GetPicture();
   fMask = pbg->GetMask();

   gVirtualX->ShapeCombineMask(fId, 0, 0, fMask);
   SetBackgroundPixmap(fPic);

   wattr.fMask = kWAOverrideRedirect | kWASaveUnder;
   wattr.fSaveUnder = kTRUE;
   wattr.fOverrideRedirect = kTRUE;

   gVirtualX->ChangeWindowAttributes(fId, &wattr);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TGCursorWindow::~TGCursorWindow()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Creates initial INACTIVE state for the recorder

TRecorder::TRecorder()
{
   fFilename = "";
   fRecorderState = new TRecorderInactive();
}

////////////////////////////////////////////////////////////////////////////////
/// Creates a recorder with filename to replay or to record,
/// depending on option (NEW or RECREATE will start recording,
/// READ will start replaying)

TRecorder::TRecorder(const char *filename, Option_t *option)
{
   TString opt(option);
   fFilename = "";
   fRecorderState = new TRecorderInactive();
   if ((opt == "NEW") || (opt == "RECREATE"))
      Start(filename, option);
   else
      Replay(filename);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TRecorder::~TRecorder()
{
   delete fRecorderState;
}

////////////////////////////////////////////////////////////////////////////////
/// Browse the recorder from a ROOT file. This allows to replay a
/// session from the browser.

void TRecorder::Browse(TBrowser *)
{
   Replay(fFilename);
}

////////////////////////////////////////////////////////////////////////////////
/// Starts recording events

void TRecorder::Start(const char *filename, Option_t *option, Window_t *w,
                      Int_t winCount)
{
   fRecorderState->Start(this, filename, option, w, winCount);
}

////////////////////////////////////////////////////////////////////////////////
/// Stopps recording events

void TRecorder::Stop(Bool_t guiCommand)
{
   fRecorderState->Stop(this, guiCommand);
}

////////////////////////////////////////////////////////////////////////////////
/// Replays events from 'filename'

Bool_t TRecorder::Replay(const char *filename, Bool_t showMouseCursor,
                         TRecorder::EReplayModes mode)
{
   return fRecorderState->Replay(this, filename, showMouseCursor, mode);
}

////////////////////////////////////////////////////////////////////////////////
/// Pauses replaying

void TRecorder::Pause()
{
   fRecorderState->Pause(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Resumes replaying

void TRecorder::Resume()
{
   fRecorderState->Resume(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Cancells replaying

void TRecorder::ReplayStop()
{
   fRecorderState->ReplayStop(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Prints out recorded commandline events

void TRecorder::ListCmd(const char *filename)
{
   fRecorderState->ListCmd(filename);
}

////////////////////////////////////////////////////////////////////////////////
/// Prints out recorded GUI events

void TRecorder::ListGui(const char *filename)
{
   fRecorderState->ListGui(filename);
}

////////////////////////////////////////////////////////////////////////////////
/// Changes state from the current to the passed one (newstate)
/// Deletes the old state if delPreviousState = KTRUE

void TRecorder::ChangeState(TRecorderState *newstate, Bool_t delPreviousState)
{
   if (delPreviousState)
      delete fRecorderState;

   fRecorderState = newstate;
}

////////////////////////////////////////////////////////////////////////////////
/// Get current state of recorder.

TRecorder::ERecorderState TRecorder::GetState() const
{
   return fRecorderState->GetState();
}


////////////////////////////////////////////////////////////////////////////////
/// Save previous canvases in a .root file

void TRecorder::PrevCanvases(const char *filename, Option_t *option)
{
   fRecorderState->PrevCanvases(filename,option);
}

//______________________________________________________________________________
// Represents state of TRecorder when replaying

ClassImp(TRecorderReplaying);

////////////////////////////////////////////////////////////////////////////////
/// Allocates all necessary data structures used for replaying
/// What is allocated here is deleted in destructor

TRecorderReplaying::TRecorderReplaying(const char *filename)
{
   fCanv = 0;
   fCmdTree = 0;
   fCmdTreeCounter = 0;
   fEventReplayed = kTRUE;
   fExtraTree = 0;
   fExtraTreeCounter = 0;
   fFilterStatusBar = kFALSE;
   fGuiTree = 0;
   fGuiTreeCounter = 0;
   fNextEvent = 0;
   fRecorder = 0;
   fRegWinCounter = 0;
   fShowMouseCursor = kTRUE;
   fWaitingForWindow = kFALSE;
   fWin = 0;
   fWinTree = 0;
   fWinTreeEntries = 0;
   fFile       = TFile::Open(filename);
   fCmdEvent   = new TRecCmdEvent();
   fGuiEvent   = new TRecGuiEvent();
   fExtraEvent = new TRecExtraEvent();
   fWindowList = new TList();
   fTimer      = new TTimer();
   fMutex      = new TMutex(kFALSE);
   if (!gCursorWin)
      gCursorWin  = new TGCursorWindow();
}

////////////////////////////////////////////////////////////////////////////////
/// Closes all signal-slot connections
/// Frees all memory allocated in contructor.

TRecorderReplaying::~TRecorderReplaying()
{
   fTimer->Disconnect(fTimer, "Timeout()", this, "ReplayRealtime()");
   fTimer->TurnOff();
   // delete fTimer;

   gClient->Disconnect(gClient, "RegisteredWindow(Window_t)", this,
                       "RegisterWindow(Window_t)");

   if (fFile) {
      fFile->Close();
      delete fFile;
   }

   delete fWindowList;
   delete fCmdEvent;
   delete fGuiEvent;
   delete fExtraEvent;
   delete fMutex;
   if (gCursorWin)
      gCursorWin->DeleteWindow();
   gCursorWin = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Initialization of data structures for replaying.
/// Start of replaying.
///
/// Return value:
///  - kTRUE  = everything is OK and replaying has begun
///  - kFALSE = non existing or invalid log file, replaying has not started

Bool_t TRecorderReplaying::Initialize(TRecorder *r, Bool_t showMouseCursor,
                                      TRecorder::EReplayModes)
{
   fWin              = 0;
   fGuiTreeCounter   = 0;
   fCmdTreeCounter   = 0;
   fExtraTreeCounter = 0;
   fRegWinCounter    = 0;
   fRecorder         = 0;

   fFilterStatusBar  = kFALSE;

   fWaitingForWindow = kFALSE;

   fEventReplayed    = 1;

   fRecorder = r;
   fShowMouseCursor = showMouseCursor;

   if (!fFile || fFile->IsZombie() || !fFile->IsOpen())
      return kFALSE;

   fCmdTree   = (TTree*) fFile->Get(kCmdEventTree);
   fWinTree   = (TTree*) fFile->Get(kWindowsTree);
   fGuiTree   = (TTree*) fFile->Get(kGuiEventTree);
   fExtraTree = (TTree*) fFile->Get(kExtraEventTree);

   if (!fCmdTree || !fWinTree || ! fGuiTree || ! fExtraTree) {
      Error("TRecorderReplaying::Initialize",
            "The ROOT file is not valid event logfile.");
      return kFALSE;
   }

   try {
      fCmdTree->SetBranchAddress(kBranchName, &fCmdEvent);
      fWinTree->SetBranchAddress(kBranchName, &fWin);
      fGuiTree->SetBranchAddress(kBranchName, &fGuiEvent);
      fExtraTree->SetBranchAddress(kBranchName, &fExtraEvent);
   }
   catch(...) {
      Error("TRecorderReplaying::Initialize",
            "The ROOT file is not valid event logfile");
      return kFALSE;
   }

   // No event to replay in given ROOT file
   if (!PrepareNextEvent()) {
      Info("TRecorderReplaying::Initialize",
           "Log file empty. No event to replay.");
      return kFALSE;
   }

   // Number of registered windows during recording
   fWinTreeEntries = fWinTree->GetEntries();

   // When a window is registered during replaying,
   // TRecorderReplaying::RegisterWindow(Window_t) is called
   gClient->Connect("RegisteredWindow(Window_t)", "TRecorderReplaying",
                    this, "RegisterWindow(Window_t)");

   Info("TRecorderReplaying::Initialize", "Replaying of file %s started",
        fFile->GetName());

   TFile *f = TFile::Open(fFile->GetName());
   if (f && !f->IsZombie()) {
      TIter nextkey(f->GetListOfKeys());
      TKey *key;
      TObject *obj;
      while ((key = (TKey*)nextkey())) {
         fFilterStatusBar = kTRUE;
         obj = key->ReadObj();
         if (!obj->InheritsFrom("TCanvas"))
            continue;
         fCanv = (TCanvas*) obj;
         fCanv->Draw();
      }
      TCanvas *canvas;
      TIter nextc(gROOT->GetListOfCanvases());
      while ((canvas = (TCanvas*)nextc())) {
         canvas->SetWindowSize(canvas->GetWindowWidth(),
                               canvas->GetWindowHeight());
      }
      fFilterStatusBar = kFALSE;
      f->Close();
   }

   gPad = 0;
   // Starts replaying
   fTimer->Connect("Timeout()","TRecorderReplaying",this,"ReplayRealtime()");
   fTimer->Start(0);

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Creates mapping for the newly registered window w and adds this
/// mapping to fWindowList
///
/// Called by signal whenever a new window is registered during replaying.
///
/// The new window ID is mapped to the old one with the same number in the
/// list of registered windows.
/// It means that 1st new window is mapped to the 1st original,
/// 2nd to the 2nd, Nth new to the Nth original.

void TRecorderReplaying::RegisterWindow(Window_t w)
{
   if (fFilterStatusBar) {
      TGWindow *win = gClient->GetWindowById(w);
      if (win) {
         if (win->GetParent()->InheritsFrom("TGStatusBar")) {
            fFilterStatusBar = kFALSE;
            return;
         }
      }
   }

   // Get original window ID that was registered as 'fRegWinCounter'th
   if (fWinTreeEntries > fRegWinCounter) {
      fWinTree->GetEntry(fRegWinCounter);
   }
   else {
      // More windows registered when replaying then when recording.
      // Cannot continue
      Error("TRecorderReplaying::RegisterWindow",
            "More windows registered than expected");
      //ReplayStop(fRecorder);
      return;
   }

   if ((gDebug > 0) && (fWaitingForWindow)) {
      std::ios::fmtflags f = std::cout.flags(); // store flags
      std::cout << " Window registered: new ID: " << std::hex << w <<
              "  previous ID: " << fWin << std::dec << std::endl;
      std::cout.flags( f ); // restore flags (reset std::hex)
   }

   // Lock mutex for guarding access to fWindowList
   fMutex->Lock();

   // Increases counter of registered windows
   fRegWinCounter++;

   // Creates new mapping of original window (fWin) and a new one (w)
   TRecWinPair *ids = new TRecWinPair(fWin, w);
   // Saves the newly created mapping
   fWindowList->Add(ids);

   // If we are waiting for this window to be registered
   // (Replaying was stopped because of that)
   if (fWaitingForWindow && fGuiEvent->fWindow == fWin) {

      if (gDebug > 0) {
         std::ios::fmtflags f = std::cout.flags(); // store flags
         std::cout << " Window " << std::hex << fGuiEvent->fWindow <<
                 " registered." << std::dec << std::endl;
         std::cout.flags( f ); // restore flags (reset std::hex)
      }

      fNextEvent = fGuiEvent;
      // Sets that we do not wait for this window anymore
      fWaitingForWindow = kFALSE;

      // Start replaying of events again
      fTimer->Start(25);
   }
   fMutex->UnLock();
}

////////////////////////////////////////////////////////////////////////////////
/// All references to the old windows (IDs) in fNextEvent are replaced by
/// new ones according to the mappings in fWindowList

Bool_t TRecorderReplaying::RemapWindowReferences()
{
   // Lock mutex for guarding access to fWindowList
   fMutex->Lock();

   TRecWinPair *ids;
   TListIter it(fWindowList);

   Bool_t found = kFALSE;

   // Iterates through the whole list of mappings
   while ((ids = (TRecWinPair*)it.Next())) {
      // Window that the event belongs to
      if (!found && fGuiEvent->fWindow == 0) {
         fGuiEvent->fWindow = gVirtualX->GetDefaultRootWindow();
         found = kTRUE;
      }
      else if (!found && ids->fKey == fGuiEvent->fWindow) {
         fGuiEvent->fWindow = ids->fValue;
         found = kTRUE;
      }
      for (Int_t i = 0; i < 5; ++i) {
         if ((Long_t) ids->fKey == fGuiEvent->fUser[i])
            fGuiEvent->fUser[i] = ids->fValue;
      }
      if (fGuiEvent->fMasked && ids->fKey == fGuiEvent->fMasked) {
         fGuiEvent->fMasked = ids->fValue;
      }
   }

   if (!found && fGuiEvent->fWindow == 0) {
      fGuiEvent->fWindow = gVirtualX->GetDefaultRootWindow();
      found = kTRUE;
   }
   // Mapping for the event found
   if (found) {
      fMutex->UnLock();
      return kTRUE;
   }

   if (gDebug > 0) {
      // save actual formatting flags
      std::ios_base::fmtflags org_flags = std::cout.flags();
      std::cout << "fGuiTreeCounter = " << std::dec << fGuiTreeCounter <<
              " No mapping found for ID " << std::hex << fGuiEvent->fWindow << std::endl;
      TRecorderInactive::DumpRootEvent(fGuiEvent,0);
      // restore original formatting flags
      std::cout.flags(org_flags);
   }

   // Stopps timer and waits for the appropriate window to be registered
   fTimer->Stop();
   fWaitingForWindow = kTRUE;

   fMutex->UnLock();
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////

Bool_t TRecorderReplaying::FilterEvent(TRecGuiEvent *e)
{
   // Not all the recorded events are replayed.
   // Some of them are generated automatically in ROOT
   // as a consequence of other events.
   //
   // RETURN VALUE:
   //    -  kTRUE  = passed TRecGuiEvent *e should be filtered
   //                (should not be replayed)
   //    -  kFALSE = passed TRecGuiEvent *e should not be filtered
   //                (should be replayed)

   // We do not replay any client messages except closing of windows
   if (e->fType == kClientMessage) {
      if ((e->fFormat == 32) && (e->fHandle != TRecGuiEvent::kROOT_MESSAGE)
          && ((Atom_t)e->fUser[0] == TRecGuiEvent::kWM_DELETE_WINDOW))
         return kFALSE;
      else
         return kTRUE;
   }

   // See TRecorderRecording::SetTypeOfConfigureNotify to get know
   // which kConfigureNotify events are filtered
   if (e->fType == kConfigureNotify && e->fUser[4] == TRecGuiEvent::kCNFilter) {
      return kTRUE;
   }

   if (e->fType == kOtherEvent) {
      if (e->fFormat >= kGKeyPress && e->fFormat < kOtherEvent)
         return kFALSE;
      return kTRUE;
   }

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Finds the next event in log file to replay and sets it to fNextEvent
///
/// Reads both from CmdTree and GuiTree and chooses that event that becomes
/// earlier
/// - fCmdTreeCounter determines actual position in fCmdTree
/// - fGuiTreeCounter determines actual position in fCmdTree
///
/// If GUI event should be replayed, we must first make sure that there is
/// appropriate mapping for this event
///
///  RETURN VALUE:
///  kFALSE = there is no event to be replayed
///  kTRUE  = there is still at least one event to be replayed. Cases:
///             - fNextEvent  = 0 => We are waiting for the appropriate
///                                  window to be registered
///             - fNextEvent != 0 => fNextEvent can be replayed (windows are
///                                  ready)

Bool_t TRecorderReplaying::PrepareNextEvent()
{
   fCmdEvent   =  0;
   fGuiEvent   =  0;
   fExtraEvent =  0;
   fNextEvent  =  0;

   // Reads the next unreplayed commandline event to fCmdEvent
   if (fCmdTree->GetEntries() > fCmdTreeCounter)
      fCmdTree->GetEntry(fCmdTreeCounter);

   // Reads the next unreplayed extra event to fExtraEvent
   if (fExtraTree->GetEntries() > fExtraTreeCounter)
      fExtraTree->GetEntry(fExtraTreeCounter);

   // Reads the next unreplayed GUI event to fGuiEvent
   // Skips GUI events that should not be replayed (FilterEvent call)
   while (fGuiTree->GetEntries() > fGuiTreeCounter) {
      fGuiTree->GetEntry(fGuiTreeCounter);
      if (!fGuiEvent || !FilterEvent(fGuiEvent))
         break;
      fGuiTreeCounter++;
   }

   // Chooses which one will be fNextEvent (the next event to be replayed)
   if (fCmdEvent && fGuiEvent && fExtraEvent) {
      // If there are all uf them, compares their times and chooses the
      // earlier one
      if ((fCmdEvent->GetTime() <= fGuiEvent->GetTime()) &&
          (fCmdEvent->GetTime() <= fExtraEvent->GetTime()))
         fNextEvent = fCmdEvent;
      else {
         if (fGuiEvent->GetTime() <= fExtraEvent->GetTime())
            fNextEvent = fGuiEvent;
         else
            fNextEvent = fExtraEvent;
      }
   }
   else if (fCmdEvent && fGuiEvent) {
      // If there are both of them, compares their times and chooses the
      // earlier one
      if (fCmdEvent->GetTime() <= fGuiEvent->GetTime())
         fNextEvent = fCmdEvent;
      else
         fNextEvent = fGuiEvent;
   }
   else if (fCmdEvent && fExtraEvent ) {
      // If there are both of them, compares their times and chooses the
      // earlier one
      if (fCmdEvent->GetTime() <= fExtraEvent->GetTime())
         fNextEvent = fCmdEvent;
      else
         fNextEvent = fExtraEvent;
   }
   else if (fGuiEvent && fExtraEvent) {
      // If there are both of them, compares their times and chooses the
      // earlier one
      if (fExtraEvent->GetTime() <= fGuiEvent->GetTime())
         fNextEvent = fExtraEvent;
      else
         fNextEvent = fGuiEvent;
   }

   // Nor commandline neither event to replay
   else if (!fCmdEvent && !fGuiEvent && !fExtraEvent)
      fNextEvent = 0;
   // Only GUI event to replay
   else if (fGuiEvent)
      fNextEvent = fGuiEvent;
   // Only commandline event to replay
   else if (fCmdEvent)
      fNextEvent = fCmdEvent;
   else
      fNextEvent = fExtraEvent;

   // Nothing to replay
   if (fNextEvent == 0)
      return kFALSE;

   // Commandline event to replay
   if (fNextEvent == fCmdEvent)
      fCmdTreeCounter++;

   // Extra event to replay
   if (fNextEvent == fExtraEvent)
      fExtraTreeCounter++;

   // GUI event to replay
   if (fNextEvent == fGuiEvent) {
      // We have the new window to send this event to
      if (RemapWindowReferences())
         fGuiTreeCounter++;
      // We do not have it yet (waiting for registraion)
      else
         fNextEvent = 0;
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// ButtonPress and ButtonRelease must be sometimes replayed more times
/// Example: pressing of a button opens small window and user chooses
/// something from that window (color)
/// Window must be opened while user is choosing

Bool_t TRecorderReplaying::CanOverlap()
{
   if (!fGuiEvent) {
      Error("TRecorderReplaying::CanOverlap()", "fGuiEvent = 0");
      return kFALSE;
   }

   // only GUI events overlapping is allowed
   if (fNextEvent->GetType() != TRecEvent::kGuiEvent)
      return kFALSE;


   if (gDebug > 0) {
      std::cout << "Event overlapping " <<
              kRecEventNames[((TRecGuiEvent*)fNextEvent)->fType] << std::endl;
      TRecorderInactive::DumpRootEvent(((TRecGuiEvent*)fNextEvent), 0);
   }

   // GUI event
   TRecGuiEvent *e  = (TRecGuiEvent*) fNextEvent;

   // Overlapping allowed for ButtonPress, ButtonRelease and MotionNotify
   if (e->fType == kButtonPress || e->fType == kButtonRelease ||
       e->fType == kMotionNotify)
      return kTRUE;

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Replays the next event.
///
/// It is called when fTimer times out.
/// Every time fTimer is set again to time equal to time difference between
/// current two events being replayed.
///
/// It can happen that execution of an event lasts different time during the
/// recording and during the replaying.
/// If fTimer times out too early and the previous event has not been yet
/// replayed, it is usually postponed in order
/// to keep events execution in the right order.
/// The excpetions are determined by TRecorderReplaying::CanOverlap()
///

void TRecorderReplaying::ReplayRealtime()
{
   UInt_t keysym;
   char str[2];

   if ((gROOT->GetEditorMode() == kText) ||
       (gROOT->GetEditorMode() == kPaveLabel)){
      gROOT->SetEditorMode();
   }

   // If there are automatically generated ROOT events in the queue, they
   // are let to be handled first
   if (gVirtualX->EventsPending()) {
      gSystem->ProcessEvents();
      return;
   }

   // Previous event has not been replayed yet and it is not allowed for
   // this event to be replayed more times
   if (!fEventReplayed && !CanOverlap())
      return;

   // Event to replay prepared
   if (fNextEvent) {
      // Sets that fNextEvent has not been replayed yet
      fEventReplayed = 0;

      // Remembers its execution time to compute time difference with
      // the next event
      fPreviousEventTime = fNextEvent->GetTime();

      // Special execution of events causing potential deadlocks
      if (fNextEvent->GetType() == TRecEvent::kGuiEvent) {
         TRecGuiEvent *ev = (TRecGuiEvent *)fNextEvent;
         if (ev->fType == kGKeyPress && ev->fState & kKeyControlMask) {
            Event_t *e = ev->CreateEvent(ev);
            gVirtualX->LookupString(e, str, sizeof(str), keysym);
            // catch the ctrl-s event
            if ((keysym & ~0x20) == kKey_S) {
               fEventReplayed = 1;
               PrepareNextEvent();
               ev->ReplayEvent(fShowMouseCursor);
               return;
            }
         }
      }

      // REPLAYS CURRENT EVENT
      fNextEvent->ReplayEvent(fShowMouseCursor);

      // Sets that fNextEvent has been replayed
      fEventReplayed = 1;
   }

   // Prepares new event for replaying
   if (!PrepareNextEvent()) {
      // No more events to be replayed (replaying has finished).

      // Switches recorder back to INACTIVE state
      Info("TRecorderReplaying::ReplayRealtime", "Replaying finished");
      fRecorder->ChangeState(new TRecorderInactive());
      return;
   }
   else {
      // We have event to replay here.

      // It will be replayed with the same time difference to the previous
      // one as when recording.
      // After given time, timer will call this method again
      if (fNextEvent)
         fTimer->Start(Long_t(fNextEvent->GetTime() - fPreviousEventTime));
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Pauses replaying

void TRecorderReplaying::Pause(TRecorder *r)
{
   fTimer->Stop();
   r->ChangeState(new TRecorderPaused(this), kFALSE);
   Info("TRecorderReplaying::Pause", "Replaying paused.");
}

////////////////////////////////////////////////////////////////////////////////
/// Cancels replaying

void TRecorderReplaying::ReplayStop(TRecorder *r)
{
   Info("TRecorderReplaying::ReplayStop", "Replaying cancelled");
   r->ChangeState(new TRecorderInactive());
}

////////////////////////////////////////////////////////////////////////////////
/// Continues previously paused replaying

void TRecorderReplaying::Continue()
{
   if (fNextEvent)
      fTimer->Start(Long_t(fNextEvent->GetTime() - fPreviousEventTime));
}

//______________________________________________________________________________
// Represents state of TRecorder after its creation

ClassImp(TRecorderInactive);

////////////////////////////////////////////////////////////////////////////////
/// Switches from INACTIVE state to RECORDING and starts recording

void TRecorderInactive::Start(TRecorder *r, const char *filename,
                              Option_t *option, Window_t *w, Int_t winCount)
{
   // const char *filename = name of ROOT file where to store recorded events
   // Option_t *option     = option for creation of ROOT file
   // Window_t *w          = list of IDs of recorder windows (if GUI for
   //                        recorder is used) [0 by default]
   // Int_t winCount       = number of IDs it this list [0 by default]

   TRecorderRecording *rec = new TRecorderRecording(r, filename, option, w, winCount);
   if (rec->StartRecording()) {
      r->ChangeState(rec);
      r->fFilename = gSystem->BaseName(filename);
   }
   else
      delete rec;
}

////////////////////////////////////////////////////////////////////////////////
/// Switches from INACTIVE state of recorder to REPLAYING
/// Return kTRUE if replaying has started or kFALSE if it is not possible
/// (bad file etc.)

Bool_t TRecorderInactive::Replay(TRecorder *r, const char *filename,
                                 Bool_t showMouseCursor,
                                 TRecorder::EReplayModes mode)
{
   // const char *filename = name of ROOT file from where to replay recorded
   // events
   // TRecorder::EReplayModes mode     = mode of replaying

   TRecorderReplaying *replay = new TRecorderReplaying(filename);

   if (replay->Initialize(r, showMouseCursor, mode)) {
      r->ChangeState(replay);
      r->fFilename = gSystem->BaseName(filename);
      return kTRUE;
   }
   else {
      delete replay;
      return kFALSE;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Prints out commandline events recorded in given file

void TRecorderInactive::ListCmd(const char *filename)
{
   /*
   if (!TClassTable::GetDict(" TRecCmdEvent")) {
      Error("TRecorderInactive::List", " TRecCmdEvent not in dictionary.");
      return;
   }*/

   TFile *file = TFile::Open(filename);
   if (!file) return;
   if (file->IsZombie() || !file->IsOpen()) {
      delete file;
      return;
   }
   TTree *t1 = (TTree*)file->Get(kCmdEventTree);

   if (!t1) {
      Error("TRecorderInactive::List",
            "The ROOT file is not valid event logfile.");
      delete file;
      return;
   }

   TRecCmdEvent *fCmdEvent = new  TRecCmdEvent();
   t1->SetBranchAddress(kBranchName, &fCmdEvent);

   Int_t entries = t1->GetEntries();
   for (Int_t i = 0; i < entries; ++i) {
      t1->GetEntry(i);
      std::cout << "[" << i << "] " << "fTime=" <<
             (ULong64_t) fCmdEvent->GetTime() << " fText=" <<
             fCmdEvent->GetText() << std::endl;
   }
   std::cout << std::endl;

   delete fCmdEvent;
   delete file;
}

////////////////////////////////////////////////////////////////////////////////
/// Prints out GUI events recorded in given file

void TRecorderInactive::ListGui(const char *filename)
{
   /*
   if (!TClassTable::GetDict("TRecGuiEvent")) {
      Error("TRecorderInactive::ListGui",
            "TRecGuiEvent not in the dictionary.");
      return;
   }*/

   TFile *file = TFile::Open(filename);
   if (!file) return;
   if (file->IsZombie() || !file->IsOpen()) {
      delete file;
      return;
   }
   TTree *t1 = (TTree*)file->Get(kGuiEventTree);

   if (!t1) {
      Error("TRecorderInactive::ListGui",
            "The ROOT file is not valid event logfile.");
      delete file;
      return;
   }

   TRecGuiEvent *guiEvent = new TRecGuiEvent();
   t1->SetBranchAddress(kBranchName, &guiEvent);

   Int_t entries = t1->GetEntries();

   for (Int_t i = 0; i < entries ; ++i) {
      t1->GetEntry(i);
      DumpRootEvent(guiEvent, i);
   }

   delete file;
   delete guiEvent;
}

////////////////////////////////////////////////////////////////////////////////
/// Prints out attributes of one GUI event TRecGuiEvent *e
/// Int_n n is number of event if called in cycle

void TRecorderInactive::DumpRootEvent(TRecGuiEvent *e, Int_t n)
{
   std::ios::fmtflags f = std::cout.flags(); // store flags
   std::cout << "[" << n << "] " << std::dec <<  std::setw(10)
      << e->GetTime().AsString() << std::setw(15) << kRecEventNames[e->fType]
      << " fW:"   << std::hex << e->fWindow
      << " t:"    << std::dec << e->fTime
      << " x:"    << DisplayValid(e->fX)
      << " y:"    << DisplayValid(e->fY)
      << " fXR:"  << DisplayValid(e->fXRoot)
      << " fYR:"  << DisplayValid(e->fYRoot)
      << " c:"    << DisplayValid(e->fCode)
      << " s:"    << DisplayValid(e->fState)
      << " w:"    << DisplayValid(e->fWidth)
      << " h:"    << DisplayValid(e->fHeight)
      << " cnt:"  << DisplayValid(e->fCount)
      << " se:"   << e->fSendEvent
      << " h:"    << e->fHandle
      << " fF:"   << DisplayValid(e->fFormat)
      << " | ";

   for (Int_t i=0; i<5; ++i)
      if (DisplayValid(e->fUser[i]) != -1)
         std::cout << "[" << i << "]=" << DisplayValid(e->fUser[i]);

   if (e->fMasked)
      std::cout << " | fM:" << std::hex << e->fMasked;

   std::cout << std::endl;
   std::cout.flags( f ); // restore flags (reset std::hex)
}

////////////////////////////////////////////////////////////////////////////////
/// Save previous canvases in a .root file

void TRecorderInactive::PrevCanvases(const char *filename, Option_t *option)
{
   fCollect = gROOT->GetListOfCanvases();
   TFile *f = TFile::Open(filename, option);
   if (f && !f->IsZombie()) {
      fCollect->Write();
      f->Close();
      delete f;
   }
}

//______________________________________________________________________________
// Represents state of TRecorder when paused

   ClassImp(TRecorderPaused);

////////////////////////////////////////////////////////////////////////////////
/// Rememeber the recorder state that is paused

TRecorderPaused::TRecorderPaused(TRecorderReplaying *state)
{
   fReplayingState = state;
}

////////////////////////////////////////////////////////////////////////////////
/// Continues replaying

void TRecorderPaused::Resume(TRecorder *r)
{
   fReplayingState->Continue();
   Info("TRecorderPaused::Resume", "Replaying resumed");

   // Switches back to the previous replaying state
   r->ChangeState(fReplayingState);
}

////////////////////////////////////////////////////////////////////////////////
/// Replaying is cancelled

void TRecorderPaused::ReplayStop(TRecorder *r)
{
   delete fReplayingState;

   Info("TRecorderReplaying::ReplayStop", "Reaplying cancelled");
   r->ChangeState(new TRecorderInactive());
}


//______________________________________________________________________________
// Represents state of TRecorder when recording events

ClassImp(TRecorderRecording);

////////////////////////////////////////////////////////////////////////////////
/// Initializes TRecorderRecording for recording
/// What is allocated here is deleted in destructor

TRecorderRecording::TRecorderRecording(TRecorder *r, const char *filename,
                                       Option_t *option, Window_t *w,
                                       Int_t winCount)
{
   fRecorder = r;
   fBeginPave = 0;

   // Remember window IDs of GUI recorder (appropriate events are
   // filtered = not recorded)
   fFilteredIdsCount = winCount;
   fFilteredIds = new Window_t[fFilteredIdsCount];
   for(Int_t i=0; i < fFilteredIdsCount; ++i)
      fFilteredIds[i] = w[i];

   // No unhandled commandline event in the beginning
   fCmdEventPending = kFALSE;

   // Filer pave events (mouse button move)
   fFilterEventPave = kFALSE;

   // No registered windows in the beginning
   fRegWinCounter = 0;

   // New timer for recording
   fTimer      = new TTimer(25, kTRUE);

   fMouseTimer = new TTimer(50, kTRUE);
   fMouseTimer->Connect("Timeout()", "TRecorderRecording", this,
                        "RecordMousePosition()");

   // File where store recorded events
   fFile       = TFile::Open(filename, option);

   // TTrees with windows, commandline events and GUi events
   fWinTree   = new TTree(kWindowsTree,    "Windows");
   fCmdTree   = new TTree(kCmdEventTree,   "Commandline events");
   fGuiTree   = new TTree(kGuiEventTree,   "GUI events");
   fExtraTree = new TTree(kExtraEventTree, "Extra events");

   fWin        = 0;
   fCmdEvent   = new TRecCmdEvent();
   fGuiEvent   = new TRecGuiEvent();
   fExtraEvent = new TRecExtraEvent();
}

////////////////////////////////////////////////////////////////////////////////
/// Freeing of allocated memory

TRecorderRecording::~TRecorderRecording()
{
   delete[] fFilteredIds;

   if (fFile)
      delete fFile;
   delete fMouseTimer;
   delete fTimer;
   delete fCmdEvent;
   delete fGuiEvent;
   delete fExtraEvent;
}

////////////////////////////////////////////////////////////////////////////////
/// Connects appropriate signals and slots in order to gain all registered
/// windows and processed events in ROOT and then starts recording

Bool_t TRecorderRecording::StartRecording()
{
   if (!fFile || fFile->IsZombie() || !fFile->IsOpen())
      return kFALSE;

   // When user types something in the commandline,
   // TRecorderRecording::RecordCmdEvent(const char *line) is called
   gApplication->Connect("LineProcessed(const char*)", "TRecorderRecording",
                         this, "RecordCmdEvent(const char*)");

   // When a new window in ROOT is registered,
   // TRecorderRecording::RegisterWindow(Window_t) is called
   gClient->Connect("RegisteredWindow(Window_t)", "TRecorderRecording", this,
                    "RegisterWindow(Window_t)");

   // When a GUI event (different from kConfigureNotify) is processed in
   // TGClient::HandleEvent or in TGClient::HandleMaskEvent,
   // TRecorderRecording::RecordGuiEvent(Event_t*, Window_t) is called
   gClient->Connect("ProcessedEvent(Event_t*, Window_t)", "TRecorderRecording",
                    this, "RecordGuiEvent(Event_t*, Window_t)");

   // When a kConfigureNotify event is processed in TGFrame::HandleEvent,
   // TRecorderRecording::RecordGuiCNEvent(Event_t*) is called
   TQObject::Connect("TGFrame", "ProcessedConfigure(Event_t*)",
                     "TRecorderRecording", this, "RecordGuiCNEvent(Event_t*)");

   // When a PaveLabel is created, TRecorderRecording::RecordPave(TObject*)
   // is called
   TQObject::Connect("TPad", "RecordPave(const TObject*)", "TRecorderRecording",
                     this, "RecordPave(const TObject*)");

   // When a Text is created, TRecorderRecording::RecordText() is called
   TQObject::Connect("TPad", "RecordLatex(const TObject*)",
                     "TRecorderRecording", this, "RecordText(const TObject*)");

   // When a PaveLabel is created, TRecorderRecording::FilterEventPave()
   // is called to filter mouse clicks events.
   TQObject::Connect("TPad", "EventPave()", "TRecorderRecording", this,
                     "FilterEventPave()");

   // When starting editing a TLatex or a TPaveLabel, StartEditing()
   // is called to memorize edition starting time.
   TQObject::Connect("TPad", "StartEditing()", "TRecorderRecording", this,
                     "StartEditing()");

   // Gui Builder specific events.
   TQObject::Connect("TGuiBldDragManager", "TimerEvent(Event_t*)",
                     "TRecorderRecording", this, "RecordGuiBldEvent(Event_t*)");

   // Creates in TTrees appropriate branches to store registered windows,
   // commandline events and GUI events
   fWinTree->Branch(kBranchName, &fWin, "fWin/l");
   fCmdTree->Branch(kBranchName, " TRecCmdEvent", &fCmdEvent);
   fGuiTree->Branch(kBranchName, "TRecGuiEvent", &fGuiEvent);
   fExtraTree->Branch(kBranchName, "TRecExtraEvent", &fExtraEvent);

   Int_t numCanvases = gROOT->GetListOfCanvases()->LastIndex();

   if (numCanvases >= 0){

      TIter nextwindow (gClient->GetListOfWindows());
      TGWindow *twin;
      Window_t  twin2;
      Int_t cnt = 0;
      while ((twin = (TGWindow*) nextwindow())) {
         twin2 = (Window_t) twin->GetId();
         if (IsFiltered(twin2)) {
            if (gDebug > 0) {
               std::cout << "WindowID "<< twin2 << " filtered" << std::endl;
            }
         }
         else if (twin != gClient->GetRoot()) {
            RegisterWindow(twin2);
         }
         cnt++;
      }
      //Info("TRecorderRecording::StartRecording", "Previous Canvases");
   }

   // Starts the timer for recording
   fTimer->TurnOn();

   // start mouse events recording timer
   fMouseTimer->Start(50);

   Info("TRecorderRecording::StartRecording", "Recording started. Log file: %s",
        fFile->GetName());

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Disconnects all slots and stopps recording.

void TRecorderRecording::Stop(TRecorder *, Bool_t guiCommand)
{
   TQObject::Disconnect("TGuiBldDragManager", "TimerEvent(Event_t*)", this,
                        "RecordGuiBldEvent(Event_t*)");
   TQObject::Disconnect("TGFrame", "ProcessedConfigure(Event_t*)", this,
                        "RecordGuiCNEvent(Event_t*)");
   TQObject::Disconnect("TPad", "RecordPave(const TObject*)", this,
                        "RecordPave(const TObject*)");
   TQObject::Disconnect("TPad", "RecordLatex(const TObject*)", this,
                        "RecordText(const TObject*)");
   TQObject::Disconnect("TPad", "EventPave()", this, "FilterEventPave()");
   TQObject::Disconnect("TPad", "StartEditing()", this, "StartEditing()");
   gClient->Disconnect(gClient, "ProcessedEvent(Event_t*, Window_t)", this,
                       "RecordGuiEvent(Event_t*, Window_t)");
   gClient->Disconnect(gClient, "RegisteredWindow(Window_t)", this,
                       "RegisterWindow(Window_t)");
   gApplication->Disconnect(gApplication, "LineProcessed(const char*)", this,
                            "RecordCmdEvent(const char*)");

   // Decides if to store the last event. It is stored if GUI recorder is used,
   // otherwise it is 'TEventRecorded::Stop' and should not be stored
   if (fCmdEventPending && guiCommand)
      fCmdTree->Fill();

   fRecorder->Write("recorder");
   fFile->Write();
   fFile->Close();
   fTimer->TurnOff();

   fMouseTimer->TurnOff();

   Info("TRecorderRecording::Stop", "Recording finished.");

   fRecorder->ChangeState(new TRecorderInactive());
}

////////////////////////////////////////////////////////////////////////////////
/// This method is called when RegisteredWindow(Window_t) is emitted from
/// TGClient.

void TRecorderRecording::RegisterWindow(Window_t w)
{
   // Stores ID of the registered window in appropriate TTree
   fWin = (ULong64_t) w;
   fWinTree->Fill();
}

////////////////////////////////////////////////////////////////////////////////
/// Records commandline event (text and time) ans saves the previous
/// commandline event
/// This 1 event delay in saving ensures that the last commandline events
/// 'TRecorder::Stop' will be not stored

void TRecorderRecording::RecordCmdEvent(const char *line)
{
   // If there is some previously recorded event, saves it in TTree now
   if (fCmdEventPending)
      fCmdTree->Fill();

   // Fill information about this new commandline event: command text and
   // time of event execution
   fCmdEvent->SetTime(fTimer->GetAbsTime());
   fCmdEvent->SetText((char*)line);

   // This event will be stored next time (if it is not the last one
   // 'TRecorder::Stop')
   fCmdEventPending = kTRUE;
   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Records GUI Event_t *e different from kConfigureNotify (they are
/// recorded in TRecorderRecording::RecordGuiCNEvent)
///
/// It is called via signal-slot when an event is processed in
/// TGClient::HandleEvent(Event_t *event)
/// or in TGClient::HandleMaskEvent(Event_t *event, Window_t wid)
///
/// If signal is emitted from TGClient::HandleEvent(Event_t *event),
/// then wid = 0

void TRecorderRecording::RecordGuiEvent(Event_t *e, Window_t wid)
{
   // If this event is caused by a recorder itself (GUI recorder),
   // it is not recorded
   if (fFilteredIdsCount && IsFiltered(e->fWindow))
      return;

   // Doesn't record the mouse clicks when a pavelabel is recorded
   if  (fFilterEventPave && (e->fCode == 1)) {
      fFilterEventPave = kFALSE;
      return;
   }
   fFilterEventPave = kFALSE;

   // don't record any copy/paste event, as event->fUser[x] parameters
   // will be invalid when replaying on a different OS
   if (e->fType == kSelectionClear || e->fType == kSelectionRequest ||
       e->fType == kSelectionNotify)
      return;

   // Copies all items of e to fGuiEvent
   CopyEvent(e, wid);

   // Saves time of recording
   fGuiEvent->SetTime(fTimer->GetAbsTime());

   // Saves recorded event itself in TTree
   fGuiTree->Fill();
}

////////////////////////////////////////////////////////////////////////////////
/// Special case for the gui builder, having a timer handling some of the
/// events.

void TRecorderRecording::RecordGuiBldEvent(Event_t *e)
{
   e->fFormat = e->fType;
   e->fType = kOtherEvent;

   // Copies all items of e to fGuiEvent
   CopyEvent(e, 0);

   // Saves time of recording
   fGuiEvent->SetTime(fTimer->GetAbsTime());

   // Saves recorded event itself in TTree
   fGuiTree->Fill();
}

////////////////////////////////////////////////////////////////////////////////
/// Try to record all mouse moves...

void TRecorderRecording::RecordMousePosition()
{
   Window_t dum;
   Event_t ev;
   ev.fCode = 0;
   ev.fType = kMotionNotify;
   ev.fState = 0;
   ev.fWindow = 0;
   ev.fUser[0] = ev.fUser[1] = ev.fUser[2] = ev.fUser[3] = ev.fUser[4] = 0;
   ev.fCount = 0;
   ev.fFormat = 0;
   ev.fHandle = 0;
   ev.fHeight = 0;
   ev.fSendEvent = 0;
   ev.fTime = 0;
   ev.fWidth = 0;

   gVirtualX->QueryPointer(gVirtualX->GetDefaultRootWindow(), dum, dum,
                           ev.fXRoot, ev.fYRoot, ev.fX, ev.fY, ev.fState);
   ev.fXRoot -= gDecorWidth;
   ev.fYRoot -= gDecorHeight;

   RecordGuiEvent(&ev, 0);
   fMouseTimer->Reset();
}

////////////////////////////////////////////////////////////////////////////////
/// Records GUI Event_t *e of type kConfigureNotify.
/// It is called via signal-slot when an kConfigureNotify event is processed
/// in TGFrame::HandleEvent

void TRecorderRecording::RecordGuiCNEvent(Event_t *e)
{
   // If this event is caused by a recorder itself, it is not recorded
   if (fFilteredIdsCount && IsFiltered(e->fWindow))
      return;

   // Sets fUser[4] value to one of EConfigureNotifyType
   // According to this value, event is or is not replayed in the future
   SetTypeOfConfigureNotify(e);

   // Copies all items of e to fGuiEvent
   CopyEvent(e, 0);

   // Saves time of recording
   fGuiEvent->SetTime(fTimer->GetAbsTime());

   // Saves recorded event itself in TTree
   fGuiTree->Fill();
}

////////////////////////////////////////////////////////////////////////////////
/// Records TPaveLabel object created in TCreatePrimitives::Pave()

void TRecorderRecording::RecordPave(const TObject *obj)
{
   Long64_t extratime = fBeginPave;
   Long64_t interval = (Long64_t)fTimer->GetAbsTime() - fBeginPave;
   TPaveLabel *pavel = (TPaveLabel *) obj;
   const char *label;
   label = pavel->GetLabel();
   TString aux = "";
   TString cad = "";
   cad = "TPaveLabel *p = new TPaveLabel(";
   cad += pavel->GetX1();
   cad += ",";
   cad += pavel->GetY1();
   cad += ",";
   cad += pavel->GetX2();
   cad += ",";
   cad += pavel->GetY2();
   cad += ",\"\"); p->Draw(); gPad->Modified(); gPad->Update();";
   Int_t i, len = (Int_t)strlen(label);
   interval /= (len + 2);
   RecordExtraEvent(cad, extratime);
   for (i=0; i < len; ++i) {
      cad = "p->SetLabel(\"";
      cad += (aux += label[i]);
      cad += "\"); ";
#ifndef R__WIN32
      cad += " p->SetTextFont(83); p->SetTextSizePixels(14); ";
#endif
      cad += " gPad->Modified(); gPad->Update();";
      extratime += interval;
      RecordExtraEvent(cad, extratime);
   }
   cad  = "p->SetTextFont(";
   cad += pavel->GetTextFont();
   cad += "); p->SetTextSize(";
   cad += pavel->GetTextSize();
   cad += "); gPad->Modified(); gPad->Update();";
   extratime += interval;
   RecordExtraEvent(cad, extratime);
}

////////////////////////////////////////////////////////////////////////////////
/// Records TLatex object created in TCreatePrimitives::Text()

void TRecorderRecording::RecordText(const TObject *obj)
{
   Long64_t extratime = fBeginPave;
   Long64_t interval = (Long64_t)fTimer->GetAbsTime() - fBeginPave;
   TLatex *texto = (TLatex *) obj;
   const char *label;
   label = texto->GetTitle();
   TString aux = "";
   TString cad = "";
   cad = "TLatex *l = new TLatex(";
   cad += texto->GetX();
   cad += ",";
   cad += texto->GetY();
   cad += ",\"\"); l->Draw(); gPad->Modified(); gPad->Update();";
   Int_t i, len = (Int_t)strlen(label);
   interval /= (len + 2);
   RecordExtraEvent(cad, extratime);
   for (i=0; i < len; ++i) {
      cad = "l->SetTitle(\"";
      cad += (aux += label[i]);
      cad += "\"); ";
#ifndef R__WIN32
      cad += " l->SetTextFont(83); l->SetTextSizePixels(14); ";
#endif
      cad += " gPad->Modified(); gPad->Update();";
      extratime += interval;
      RecordExtraEvent(cad, extratime);
   }
   cad  = "l->SetTextFont(";
   cad += texto->GetTextFont();
   cad += "); l->SetTextSize(";
   cad += texto->GetTextSize();
   cad += "); gPad->Modified(); gPad->Update();";
   cad += " TVirtualPad *spad = gPad->GetCanvas()->GetSelectedPad();";
   cad += " gPad->GetCanvas()->Selected(spad, l, kButton1Down);";
   extratime += interval;
   RecordExtraEvent(cad, extratime);
}

////////////////////////////////////////////////////////////////////////////////
/// Change the state of the flag to kTRUE when you are recording a pavelabel.

void TRecorderRecording::FilterEventPave()
{
   fFilterEventPave = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Memorize the starting time of editinga TLatex or a TPaveLabel

void TRecorderRecording::StartEditing()
{
   fBeginPave = fTimer->GetAbsTime();
}

////////////////////////////////////////////////////////////////////////////////
/// Records TLatex or TPaveLabel object created in TCreatePrimitives,
/// ExtTime is needed for the correct replay of these events.

void TRecorderRecording::RecordExtraEvent(TString line, TTime extTime)
{
   fExtraEvent->SetTime(extTime);
   fExtraEvent->SetText(line);
   fExtraTree->Fill();
}

////////////////////////////////////////////////////////////////////////////////
/// Copies all items of given event to fGuiEvent

void TRecorderRecording::CopyEvent(Event_t *e, Window_t wid)
{
   fGuiEvent->fType     = e->fType;
   fGuiEvent->fWindow   = e->fWindow;
   fGuiEvent->fTime     = e->fTime;

   fGuiEvent->fX        = e->fX;
   fGuiEvent->fY        = e->fY;
   fGuiEvent->fXRoot    = e->fXRoot;
   fGuiEvent->fYRoot    = e->fYRoot;

   fGuiEvent->fCode     = e->fCode;
   fGuiEvent->fState    = e->fState;

   fGuiEvent->fWidth    = e->fWidth;
   fGuiEvent->fHeight   = e->fHeight;

   fGuiEvent->fCount       = e->fCount;
   fGuiEvent->fSendEvent   = e->fSendEvent;
   fGuiEvent->fHandle      = e->fHandle;
   fGuiEvent->fFormat      = e->fFormat;

   if (fGuiEvent->fHandle == gROOT_MESSAGE)
      fGuiEvent->fHandle = TRecGuiEvent::kROOT_MESSAGE;

   for(Int_t i=0; i<5; ++i)
      fGuiEvent->fUser[i] = e->fUser[i];

   if (fGuiEvent->fUser[0] == (Int_t)gWM_DELETE_WINDOW)
      fGuiEvent->fUser[0] = TRecGuiEvent::kWM_DELETE_WINDOW;

   if (e->fType == kGKeyPress || e->fType == kKeyRelease) {
      char tmp[10] = {0};
      UInt_t keysym = 0;
      gVirtualX->LookupString(e, tmp, sizeof(tmp), keysym);
      fGuiEvent->fCode = keysym;
   }

   fGuiEvent->fMasked  = wid;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns kTRUE if passed id belongs to window IDs of recorder GUI itself

Bool_t TRecorderRecording::IsFiltered(Window_t id)
{
   for(Int_t i=0; i < fFilteredIdsCount; ++i)
      if (id == fFilteredIds[i])
         return kTRUE;

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Sets type of kConfigureNotify event to one of EConfigureNotify
///
/// On Linux paremeters of GUI event kConfigureNotify are different
/// than parameters of the same event executed on Windows.
/// Therefore we need to distinguish [on Linux], if the event is movement
/// or resize event.
/// On Windows, we do not need to distinguish them.

void TRecorderRecording::SetTypeOfConfigureNotify(Event_t *e)
{
   // On both platforms, we mark the events matching the criteria
   // (automatically generated in ROOT) as events that should be filtered
   // when replaying (TRecGuiEvent::kCNFilter)
   if ((e->fX == 0 && e->fY == 0)) { // || e->fFormat == 32 ) {
      e->fUser[4] = TRecGuiEvent::kCNFilter;
      return;
   }

#ifdef WIN32

   // No need to distinguish between move and resize on Windows
   e->fUser[4] = TRecGuiEvent::kCNMoveResize;

#else

   TGWindow *w = gClient->GetWindowById(e->fWindow);
   if (w) {
      TGFrame *t = (TGFrame *)w;

      // If this event does not cause any change in position or size ->
      // automatically generated event
      if (t->GetWidth() == e->fWidth && t->GetHeight() == e->fHeight &&
          e->fX == t->GetX() && e->fY == t->GetY()) {
         e->fUser[4] = TRecGuiEvent::kCNFilter;
      }
      else {
         // Size of the window did not change -> move
         if (t->GetWidth() == e->fWidth && t->GetHeight() == e->fHeight) {
            e->fUser[4] = TRecGuiEvent::kCNMove;
         }
         // Size of the window changed -> resize
         else {
            e->fUser[4] = TRecGuiEvent::kCNResize;
         }
      }
   }

#endif
}



//______________________________________________________________________________
// The GUI for the recorder

ClassImp(TGRecorder);

////////////////////////////////////////////////////////////////////////////////
/// The GUI for the recorder

TGRecorder::TGRecorder(const TGWindow *p, UInt_t w, UInt_t h) :
   TGMainFrame(p ? p : gClient->GetRoot(), w, h)
{
   TGHorizontalFrame *hframe;
   TGVerticalFrame *vframe;
   SetCleanup(kDeepCleanup);
   fRecorder = new TRecorder();
   fFilteredIds[0] = GetId();

   // Create a horizontal frame widget with buttons
   hframe = new TGHorizontalFrame(this, 200, 75, kChildFrame | kFixedHeight,
                                  (Pixel_t)0x000000);
   fFilteredIds[1] = hframe->GetId();

   // LABEL WITH TIME

   vframe = new TGVerticalFrame(hframe, 200, 75, kChildFrame | kFixedHeight,
                                (Pixel_t)0x000000);
   fFilteredIds[2] = vframe->GetId();

   TGLabel *fStatusLabel = new TGLabel(vframe, "Status:");
   fStatusLabel->SetTextColor(0x7cffff);
   fStatusLabel->SetBackgroundColor((Pixel_t)0x000000);
   vframe->AddFrame(fStatusLabel, new TGLayoutHints(kLHintsLeft | kLHintsTop,
                                                    2, 2, 2, 2));
   fFilteredIds[3] = fStatusLabel->GetId();

   TGLabel *fTimeLabel = new TGLabel(vframe, "Time: ");
   fTimeLabel->SetTextColor(0x7cffff);
   fTimeLabel->SetBackgroundColor((Pixel_t)0x000000);
   vframe->AddFrame(fTimeLabel, new TGLayoutHints(kLHintsLeft | kLHintsTop,
                                                  2, 2, 13, 2));
   fFilteredIds[4] = fTimeLabel->GetId();

   hframe->AddFrame(vframe, new TGLayoutHints(kLHintsLeft | kLHintsExpandY));

   vframe = new TGVerticalFrame(hframe, 200, 75, kChildFrame | kFixedHeight,
                                (Pixel_t)0x000000);
   fFilteredIds[5] = vframe->GetId();

   fStatus = new TGLabel(vframe, "Inactive");
   fStatus->SetTextColor(0x7cffff);
   fStatus->SetBackgroundColor((Pixel_t)0x000000);
   vframe->AddFrame(fStatus, new TGLayoutHints(kLHintsLeft | kLHintsTop,
                                               2, 2, 2, 2));
   fFilteredIds[6] = fStatus->GetId();

   fTimeDisplay = new TGLabel(vframe, "00:00:00");
   fTimeDisplay->SetTextColor(0x7cffff);
   fTimeDisplay->SetTextFont("Helvetica -34", kFALSE);
   fTimeDisplay->SetBackgroundColor((Pixel_t)0x000000);
   vframe->AddFrame(fTimeDisplay, new TGLayoutHints(kLHintsLeft | kLHintsTop,
                                                    2, 2, 2, 2));
   fFilteredIds[7] = fTimeDisplay->GetId();

   hframe->AddFrame(vframe, new TGLayoutHints(kLHintsLeft | kLHintsExpandY,
                                              10, 0, 0, 0));
   AddFrame(hframe, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

   // Create a horizontal frame widget with buttons
   hframe = new TGHorizontalFrame(this, 200, 200);
   fFilteredIds[8] = hframe->GetId();

   // START-STOP button
   fStartStop = new TGPictureButton(hframe,gClient->GetPicture("record.png"));
   fStartStop->SetStyle(gClient->GetStyle());
   fStartStop->Connect("Clicked()","TGRecorder",this,"StartStop()");
   hframe->AddFrame(fStartStop, new TGLayoutHints(kLHintsLeft | kLHintsTop,
                                                  2, 2, 2, 2));
   fStartStop->Resize(40,40);
   fFilteredIds[9] = fStartStop->GetId();

   // REPLAY button
   fReplay = new TGPictureButton(hframe,gClient->GetPicture("replay.png"));
   fReplay->SetStyle(gClient->GetStyle());
   fReplay->Connect("Clicked()","TGRecorder",this,"Replay()");
   hframe->AddFrame(fReplay, new TGLayoutHints(kLHintsLeft | kLHintsTop,
                                               2, 2, 2, 2));
   fReplay->Resize(40,40);
   fFilteredIds[10] = fReplay->GetId();

   // MOUSE CURSOR CHECKBOX
   fCursorCheckBox = new TGCheckButton(this,"Show mouse cursor");
   AddFrame(fCursorCheckBox, new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2));
   fFilteredIds[11] = fCursorCheckBox->GetId();

   // Timer
   fTimer = new TTimer(25);
   fTimer->Connect("Timeout()", "TGRecorder", this, "Update()");

   AddFrame(hframe, new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2));

   SetEditDisabled(kEditDisable | kEditDisableGrab);
   SetWindowName("ROOT Event Recorder");
   MapSubwindows();
   Layout();
   MapWindow();

   SetDefault();
}

////////////////////////////////////////////////////////////////////////////////
/// Sets GUI to the default inactive state

void TGRecorder::SetDefault()
{
   fTimeDisplay->SetText("00:00:00");

   fReplay->SetPicture(gClient->GetPicture("replay.png"));
   fReplay->SetEnabled(kTRUE);

   fCursorCheckBox->SetEnabled(kTRUE);
   fCursorCheckBox->SetOn(kTRUE);

   fStartStop->SetPicture(gClient->GetPicture("record.png"));
   fStartStop->SetEnabled(kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Called when fTimer timeouts (every 0.025 second)
/// Updates GUI of recorder

void TGRecorder::Update()
{
   struct tm *running;
   static int cnt = 0;
   TString stime;
   time( &fElapsed );
   time_t elapsed_time = (time_t)difftime( fElapsed, fStart );
   running = gmtime( &elapsed_time );

   switch(fRecorder->GetState()) {

      // When recording or replaying, updates timer
      // and displays new value of seconds counter
      case TRecorder::kRecording:
      case TRecorder::kReplaying:

         // Every whole second: updates timer and displays new value
         if (cnt >= 10) {
            if (fRecorder->GetState() == TRecorder::kReplaying)
               fStatus->SetText("Replaying");
            else
               fStatus->SetText("Recording");
            stime.Form("%02d:%02d:%02d", running->tm_hour,
                        running->tm_min, running->tm_sec);
            fTimeDisplay->SetText(stime.Data());

            cnt = 0;
            if (gVirtualX->EventsPending()) {
               fStatus->SetText("Waiting...");
               fStatus->SetTextColor((Pixel_t)0xff0000);
            }
            else {
               fStatus->SetTextColor((Pixel_t)0x7cffff);
            }
            fStatus->Resize();
            fTimeDisplay->Resize();
         }
         else
            ++cnt;

         // Changes background color according to the queue of pending events
         fTimer->Reset();
         break;

      // End of replaying or recording. Sets recorder GUI to default state
      case TRecorder::kInactive:
         fStatus->SetText("Inactive");
         fStatus->SetTextColor((Pixel_t)0x7cffff);
         fStatus->Resize();
         fTimer->TurnOff();
         SetDefault();
         break;

      default:
         break;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Handles push of the fStartStop button
/// according to the current recorder state

void TGRecorder::StartStop()
{
   static const char *gFiletypes[] = {
      "All files", "*", "Text files", "*.txt", "ROOT files", "*.root", 0, 0
   };
   TGFileInfo fi;

   switch(fRecorder->GetState()) {

      // Starts recording
      case TRecorder::kInactive:

         fi.fFileTypes = gFiletypes;
         fi.fOverwrite = kFALSE;

         new TGFileDialog(gClient->GetDefaultRoot(),
                          gClient->GetDefaultRoot(),
                          kFDSave,&fi);

         if (fi.fFilename && strlen(fi.fFilename)) {

            if (!gROOT->GetListOfCanvases()->IsEmpty()) {
               fRecorder->PrevCanvases(fi.fFilename, "RECREATE");
               fRecorder->Start(fi.fFilename, "UPDATE", fFilteredIds,
                                fgWidgetsCount);
            }
            else {
               fRecorder->Start(fi.fFilename, "RECREATE", fFilteredIds,
                                fgWidgetsCount);
            }
            fCursorCheckBox->SetDisabledAndSelected(kTRUE);
            fStartStop->SetPicture(gClient->GetPicture("stop.png"));
            fReplay->SetEnabled(kFALSE);
            fTimer->TurnOn();
            time( &fStart );
         }
         break;

      // Stops recording
      case TRecorder::kRecording:
         fRecorder->Stop(kTRUE);
         break;

      // Pauses replaying
      case TRecorder::kReplaying:
         fRecorder->Pause();
         fStartStop->SetPicture(gClient->GetPicture("replay.png"));
         break;

      // Resumes replaying
      case TRecorder::kPaused:
         fRecorder->Resume();
         fStartStop->SetPicture(gClient->GetPicture("pause.png"));
         break;

      default:
         break;
   } // switch
}

////////////////////////////////////////////////////////////////////////////////
/// Handles push of fReplay button
/// according to the current recorder state

void TGRecorder::Replay()
{
   TGFileInfo fi;

   switch(fRecorder->GetState()) {

      // Starts replaying
      case TRecorder::kInactive:

         new TGFileDialog(gClient->GetDefaultRoot(),
                          gClient->GetDefaultRoot(),
                          kFDOpen, &fi);

         if (fi.fFilename && strlen(fi.fFilename)) {
            if (fRecorder->Replay(fi.fFilename, fCursorCheckBox->IsOn())) {

               fTimer->TurnOn();
               time( &fStart );

               fReplay->SetPicture(gClient->GetPicture("stop.png"));
               fStartStop->SetPicture(gClient->GetPicture("pause.png"));

               if (fCursorCheckBox->IsOn())
                  fStartStop->SetEnabled(kFALSE);

               fCursorCheckBox->SetEnabled(kFALSE);
            }
         }
         break;

      // Stops replaying
      case TRecorder::kReplaying:
      case TRecorder::kPaused:
         fRecorder->ReplayStop();
         break;

      default:
         break;

   } // switch
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor. Cleanup the GUI.

TGRecorder::~TGRecorder()
{
   fTimer->TurnOff();
   delete fTimer;
   Cleanup();
}

//______________________________________________________________________________
// Helper class

ClassImp(TRecCmdEvent);
ClassImp(TRecGuiEvent);

////////////////////////////////////////////////////////////////////////////////
/// Replays stored GUI event

void TRecGuiEvent::ReplayEvent(Bool_t showMouseCursor)
{
   Int_t    px, py, dx, dy;
   Window_t wtarget;
   Event_t *e = CreateEvent(this);

   // don't try to replay any copy/paste event, as event->fUser[x]
   // parameters are invalid on different OSes
   if (e->fType == kSelectionClear || e->fType == kSelectionRequest ||
       e->fType == kSelectionNotify) {
      delete e;
      return;
   }

   // Replays movement/resize event
   if (e->fType == kConfigureNotify) {
      TGWindow *w = gClient->GetWindowById(e->fWindow);

      // Theoretically, w should always exist (we found the right mapping,
      // otherwise we would not get here).
      // Anyway, it can happen that it was destroyed by some earlier ROOT event
      // We give higher priority to automatically generated
      // ROOT events in TRecorderReplaying::ReplayRealtime.

      if (w) {
         WindowAttributes_t attr;
         if (e->fUser[4] == TRecGuiEvent::kCNMove) {
            // Linux: movement of the window
            // first get window attribute to compensate the border size
            gVirtualX->GetWindowAttributes(e->fWindow, attr);
            if ((e->fX - attr.fX > 0) && (e->fY - attr.fY > 0))
               w->Move(e->fX - attr.fX, e->fY - attr.fY);
         }
         else {
            if (e->fUser[4] == TRecGuiEvent::kCNResize) {
               // Linux: resize of the window
               w->Resize(e->fWidth, e->fHeight);
            }
            else {
               if (e->fUser[4] == TRecGuiEvent::kCNMoveResize) {
                  // Windows: movement or resize of the window
                  w->MoveResize(e->fX, e->fY, e->fWidth, e->fHeight);
               }
               else {
                  if (gDebug > 0)
                     Error("TRecGuiEvent::ReplayEvent",
                           "kConfigureNotify: Unknown value: fUser[4] = %ld ",
                           e->fUser[4]);
               }
            }
         }
      }
      else {
         // w = 0
         if (gDebug > 0)
            Error("TRecGuiEvent::ReplayEvent",
                  "kConfigureNotify: Window does not exist anymore ");
      }
      delete e;
      return;

   } // kConfigureNotify

   if (showMouseCursor && e->fType == kButtonPress) {
      gVirtualX->TranslateCoordinates(e->fWindow, gVirtualX->GetDefaultRootWindow(),
                                      e->fX, e->fY, px, py, wtarget);
      dx = px - gCursorWin->GetX();
      dy = py - gCursorWin->GetY();
      if (TMath::Abs(dx) > 5) gDecorWidth += dx;
      if (TMath::Abs(dy) > 5) gDecorHeight += dy;
   }
   // Displays fake mouse cursor for MotionNotify event
   if (showMouseCursor && e->fType == kMotionNotify) {
      if (gCursorWin && e->fWindow == gVirtualX->GetDefaultRootWindow()) {
         if (!gCursorWin->IsMapped()) {
            gCursorWin->MapRaised();
         }
         if (gVirtualX->GetDrawMode() == TVirtualX::kCopy) {
//#ifdef R__MACOSX
            // this may have side effects (e.g. stealing focus)
            gCursorWin->RaiseWindow();
//#endif
            gCursorWin->Move(e->fXRoot + gDecorWidth, e->fYRoot + gDecorHeight);
         }
      }
   }

   // Lets all the other events to be handled the same way as when recording
   // first, special case for the gui builder, having a timer handling
   // some of the events
   if (e->fType == kOtherEvent && e->fFormat >= kGKeyPress &&
       e->fFormat < kOtherEvent) {
      e->fType = (EGEventType)e->fFormat;
      if (gDragManager)
         gDragManager->HandleTimerEvent(e, 0);
      delete e;
      return;
   }
   else { // then the normal cases
      if (!fMasked)
         gClient->HandleEvent(e);
      else
         gClient->HandleMaskEvent(e, fMasked);
   }
   delete e;
}

////////////////////////////////////////////////////////////////////////////////
/// Converts TRecGuiEvent type to Event_t type

Event_t *TRecGuiEvent::CreateEvent(TRecGuiEvent *ge)
{
   Event_t *e = new Event_t();

   // Copies all data items

   e->fType   = ge->fType;
   e->fWindow = ge->fWindow;
   e->fTime   = ge->fTime;

   e->fX = ge->fX;
   e->fY = ge->fY;
   e->fXRoot = ge->fXRoot;
   e->fYRoot = ge->fYRoot;

   e->fCode   = ge->fCode;
   e->fState  = ge->fState;

   e->fWidth  = ge->fWidth;
   e->fHeight = ge->fHeight;

   e->fCount  = ge->fCount;
   e->fSendEvent = ge->fSendEvent;

   e->fHandle = ge->fHandle;
   e->fFormat = ge->fFormat;

   if (e->fHandle == TRecGuiEvent::kROOT_MESSAGE)
      e->fHandle = gROOT_MESSAGE;

   for(Int_t i=0; i<5; ++i)
      e->fUser[i] = ge->fUser[i];

   if (e->fUser[0] == TRecGuiEvent::kWM_DELETE_WINDOW)
      e->fUser[0] = gWM_DELETE_WINDOW;

   if (ge->fType == kGKeyPress || ge->fType == kKeyRelease) {
      e->fCode    = gVirtualX->KeysymToKeycode(ge->fCode);
#ifdef R__WIN32
      e->fUser[1] = 1;
      e->fUser[2] = e->fCode;
#endif
   }

   return e;
}

ClassImp(TRecWinPair);
