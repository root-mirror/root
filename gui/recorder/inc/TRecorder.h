// @(#)root/gui:$Id$
// Author: Katerina Opocenska   11/09/2008

/*************************************************************************
* Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
* All rights reserved.                                                  *
*                                                                       *
* For the licensing terms see $ROOTSYS/LICENSE.                         *
* For the list of contributors see $ROOTSYS/README/CREDITS.             *
*************************************************************************/

#ifndef ROOT_TRecorder
#define ROOT_TRecorder

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

#include "TApplication.h"
#include "TError.h"
#include "TTimer.h"
#include "TGClient.h"
#include "TGFrame.h"
#include "TCanvas.h"
#include "THashList.h"

#include <ctime>
#include <iostream>

class TMutex;
class TTree;
class TFile;
class TGPictureButton;
class TGCheckButton;
class TGLabel;
class TRecorderState;

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  TRecEvent                                                           //
//                                                                      //
//  Abstract class that defines interface for a class storing           //
//  information about 1 ROOT event.                                     //
//  Time of event is stored and this event can be replayed.             //
//  Classes TRecCmdEvent and TRecGuiEvent implements this interface     //
//  for command line and GUI events respectively.                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

class TRecEvent : public TObject
{
private:
   TTime   fEventTime;          // Time of original event execution

public:
   //---- Types of events recorded in ROOT.
   enum ERecEventType {
      kCmdEvent,     // Commandline event
      kGuiEvent,    // GUI event
      kExtraEvent
   };

   // Replays (executes) the stored event again
   virtual void ReplayEvent(Bool_t showMouseCursor = kTRUE) = 0;

   // Returns what kind of event it stores
   virtual ERecEventType GetType() const = 0;

   virtual TTime GetTime() const {
      // Returns time of original execution of stored event
      return fEventTime;
   }

   virtual void SetTime(TTime t) {
      // Sets time of event execution
      fEventTime = t;
   }

   ClassDef(TRecEvent,1) // Abstract class. Defines basic interface for storing information about ROOT events
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  TRecCmdEvent                                                        //
//                                                                      //
//  Class used for storing information about 1 commandline event.       //
//  It means 1 command typed in by user in the commandline,             //
//  e.g 'new TCanvas'.                                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

class TRecCmdEvent : public TRecEvent
{
private:
   TString fText;             // Text of stored command

public:
   TRecCmdEvent() {
      // Creates new empty  TRecCmdEvent
   }

   void SetText(const char *text) {
      // Saves text of a command
      fText = text;
   }

   const char *GetText() const {
      // Returns stored text of the command
      return fText.Data();
   }

   virtual ERecEventType GetType() const {
      // Returns what kind of event it stores (commandline event)
      return TRecEvent::kCmdEvent;
   }

   virtual void ReplayEvent(Bool_t) {
      // Stored command is executed again
      std::cout << GetText() << std::endl;
      gApplication->ProcessLine(GetText());
   }

   ClassDef(TRecCmdEvent,1) // Class stores information about 1 commandline event (= 1 command typed by user in commandline)
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  TRecExtraEvent                                                      //
//                                                                      //
//  Class used for storing information about 1 extra event.             //
//  It means 1 TPaveLabel or 1 TLatex event produced in the Canvas      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class TRecExtraEvent : public TRecEvent
{
private:
   TString fText;             // Text of stored command

public:
   TRecExtraEvent() {
      // Creates new empty  TRecExtraEvent
   }

   void SetText(TString text) {
      // Saves text of a command (PaveLabel or Text)
      fText = text;
   }

   TString GetText() const {
      // Returns stored text of the command
      return fText;
   }

   virtual ERecEventType GetType() const {
      // Returns what kind of event it stores (Especial event)
      return TRecEvent::kExtraEvent;
   }

   virtual void ReplayEvent(Bool_t) {
      // Stored event is executed again

      gApplication->ProcessLine(GetText());
   }

   ClassDef(TRecExtraEvent,1) // Class stores information about extra events
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  TRecGuiEvent                                                        //
//                                                                      //
//  Class used for storing information about 1 GUI event in ROOT.       //
//  For list of possible GUI events see EGEventType.                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

class TRecGuiEvent : public TRecEvent
{
protected:
   friend class TRecorderInactive;
   friend class TRecorderPaused;
   friend class TRecorderRecording;
   friend class TRecorderReplaying;

   EGEventType    fType;            // Type of event (see EGEventType)
   Window_t       fWindow;          // Window ID which reported event is relative to
   Time_t         fTime;            // Time event occured in ms
   Int_t          fX;               // Pointer x coordinate in event window
   Int_t          fY;               // Pointer y coordinate in event window
   Int_t          fXRoot;           // x coordinate relative to root
   Int_t          fYRoot;           // y coordinate relative to root
   UInt_t         fCode;            // Key or button code
   UInt_t         fState;           // Key or button mask
   UInt_t         fWidth;           // Width of exposed area
   UInt_t         fHeight;          // Height of exposed area
   Int_t          fCount;           // If non-zero, at least this many more exposes
   Bool_t         fSendEvent;       // True if event came from SendEvent
   Handle_t       fHandle;          // General resource handle (used for atoms or windows)
   Int_t          fFormat;          // Next fields only used by kClientMessageEvent
   Long_t         fUser[5];         // 5 longs can be used by client message events
                                    // NOTE: only [0], [1] and [2] may be used.
                                    // [1] and [2] may contain > 32 bit quantities
                                    // (i.e. pointers on 64 bit machines)
   Window_t       fMasked;          // If non-zero, event recorded in HandleMaskEvent()

public:
   //---- Types of kConfigureNotify GUI event
   enum EConfigureNotifyType {
      kCNMove       = 0,      // Movement of a window (Linux)
      kCNResize     = 1,      // Resize of a window (Linux)
      kCNMoveResize = 2,      // Movement, resize or both (Windows)
      kCNFilter     = 3       // Not replaybale (filtered event).
   };
   //---- Aliases for non cross-platform atoms.
   enum ERootAtoms {
      kWM_DELETE_WINDOW = 10001,
      kROOT_MESSAGE     = 10002
   };

   virtual ERecEventType GetType() const {
      // Returns what kind of event it stores (GUI event)
      return TRecEvent::kGuiEvent;
   }

   virtual void    ReplayEvent(Bool_t showMouseCursor = kTRUE);
   static Event_t *CreateEvent(TRecGuiEvent *ge);

   ClassDef(TRecGuiEvent,1) // Class stores information about 1 GUI event in ROOT
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  TRecWinPair                                                         //
//                                                                      //
//  Class used for storing of window IDs mapping.                       //
//  Remapping of window IDs is needed for replaying events.             //
//  - ID of original window is stored in fKey.                          //
//  - ID of a new window is stored in fValue.                           //
//                                                                      //
//  Whenever an event is replayed, its referenced window ID is changed  //
//  from original to a new one according to the appropriate mapping.    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

class TRecWinPair : public TObject
{
protected:
   friend class TRecorderReplaying;

   Window_t    fKey;    // ID of original window (for which an event was originally recorded)
   Window_t    fValue;  // ID of a new window (for which an event is being replayed)

public:
   // Creates a new key-value mapping of window IDs
   TRecWinPair(Window_t key, Window_t value): fKey(key), fValue(value) {}

   ClassDef(TRecWinPair,1) // Class used for storing of window IDs mapping. Needed for replaying events.
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  TRecorder                                                           //
//                                                                      //
//  Class provides direct recorder/replayer interface for a user.       //
//  See 'ROOT EVENT RECORDING SYSTEM' for more information about usage. //
//                                                                      //
//  Implementation uses C++ design pattern State. Functionality of      //
//  recorder is divided into 4 classes according to the current         //
//  state of recorder.                                                  //
//                                                                      //
//  Internally, there is a pointer to TRecorderState object.            //
//  This object changes whenever state of recorder is changed.          //
//  States of recorder are the following:                               //
//                                                                      //
//  - INACTIVE  Implemented in TRecorderInactive class.                 //
//              Default state after TRecorder object is created.        //
//                                                                      //
//  - RECORDING Implemented in TRecorderRecording class.                //
//                                                                      //
//  - REPLAYING Implemented in TRecorderReplaying class.                //
//                                                                      //
//  - PAUSED    Implemented in TRecorderPause class.                    //
//              Pause of replaying.                                     //
//                                                                      //
//  Every command for TRecorder is just passed                          //
//  to TRecordeState object.                                            //
//  Depending on the current state of recorder, this command is passed  //
//  to some of the above mentioned classes and if valid, handled there. //
//                                                                      //
//  [TRecorder.JPG]                                                     //
//                                                                      //
//  Switching between states is not possible from outside. States are   //
//  switched directly by state objects via:                             //
//                                                                      //
//  ChangeState(TRecorderState* newstate, Bool_t deletePreviousState);  //
//                                                                      //
//  When recorder is switched to a new state, the old state object is   //
//  typically deleted. The only exception is switching from REPLAYING   //
//  state to PAUSED state. The previous state (REPLAYING) is not        //
//  deleted in order to be used again after TRecorder::Resume call.     //
//                                                                      //
//  STATE TRANSITIONS:                                                  //
//  ------------------                                                  //
//                                                                      //
//  INACTIVE  -> RECORDING via TRecorder::Start (Starts recording)      //
//  RECORDING -> INACTIVE  via TRecorder::Stop  (Stops recording)       //
//                                                                      //
//  INACTIVE  -> REPLAYING via TRecorder::Replay     (Starts replaying) //
//  REPLAYING -> INACTIVE  via TRecorder::ReplayStop (Stops replaying)  //
//                                                                      //
//  REPLAYING -> PAUSED    via TRecorder::Pause  (Pause replaying)      //
//  PAUSED    -> REPLAYING via TRecorder::Resume (Resumes replaying)    //
//                                                                      //
//  PAUSED    -> INACTIVE  via TRecorder::ReplayStop (Stops paused      //
//                                                    replaying)        //
//                                                                      //
// [TRecorderStates.JPG]                                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class TRecorder : public TObject
{
private:
   TRecorderState *fRecorderState;   //! Current state of recorder

   TRecorder(const TRecorder&);            // Not implemented.
   TRecorder &operator=(const TRecorder&); // Not implemented.

protected:
   friend class TRecorderState;
   friend class TRecorderInactive;
   friend class TRecorderPaused;
   friend class TRecorderRecording;
   friend class TRecorderReplaying;

   TString      fFilename;           // Events file name
   // Changes state to the new one.
   // See class documentation for information about state changing.
   void  ChangeState(TRecorderState *newstate, Bool_t deletePreviousState = kTRUE);

public:
   //---- Modes of replaying. Only kRealtime implemented so far
   enum EReplayModes {
      kRealtime
   };
   //---- States of recorder. In every moment, recorder is in right
   // one of these states.
   enum ERecorderState {
      kInactive,
      kRecording,
      kPaused,
      kReplaying
   };

   // Creates recorder and sets its state as INACTIVE
   TRecorder();
   TRecorder(const char *filename, Option_t *option = "READ");

   // Deletes recorder together with its current state
   virtual ~TRecorder();

   void Browse(TBrowser *);

   // Starts recording of events to the given file
   void Start(const char *filename, Option_t *option = "RECREATE", Window_t *w = 0, Int_t winCount = 0);

   // Stops recording of events
   void Stop(Bool_t guiCommand = kFALSE);

   // Replays recorded events from given file
   Bool_t Replay(const char *filename, Bool_t showMouseCursor = kTRUE, TRecorder::EReplayModes mode = kRealtime);

   // Replays recorded events from current file
   void Replay() { Replay(fFilename); }   // *MENU*

   // Pauses replaying
   void Pause();

   // Resumes paused replaying
   void Resume();

   // Stops (cancels) replaying
   void ReplayStop();

   // Prints out the list of recorded commandline events
   void ListCmd(const char *filename);

   // Prints out the list of recorded GUI events
   void ListGui(const char *filename);

   // Gets current state of recorder
   virtual TRecorder::ERecorderState GetState() const;

   // Saves all the canvases previous to the TRecorder
   void PrevCanvases(const char *filename, Option_t *option);

   ClassDef(TRecorder,2) // Class provides direct recorder/replayer interface for a user.
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  TRecorderState                                                      //
//                                                                      //
//  Abstract class that defines interface for a state of recorder.      //
//  Inherited classes are:                                              //
//  - TRecorderInactive                                                 //
//  - TRecorderRecording                                                //
//  - TRecorderReplaying                                                //
//  - TRecorderPaused                                                   //
//                                                                      //
//  See TRecorder for more information about creating, using,           //
//  changing and deleting states.                                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class TRecorderState
{
protected:
   friend class TRecorder;
   void ChangeState(TRecorder *r, TRecorderState *s, Bool_t deletePreviousState) { r->ChangeState(s, deletePreviousState); }

public:
   virtual ~TRecorderState() {}
   virtual void   Start(TRecorder *, const char *, Option_t *, Window_t *, Int_t) {}
   virtual void   Stop(TRecorder *, Bool_t ) {}
   virtual Bool_t Replay(TRecorder *, const char *, Bool_t, TRecorder::EReplayModes) { return false; }
   virtual void   Pause(TRecorder *) {}
   virtual void   Resume(TRecorder *) {}
   virtual void   ReplayStop(TRecorder *) {}

   virtual void   ListCmd(const char *) {}
   virtual void   ListGui(const char *) {}

   virtual void   PrevCanvases(const char *, Option_t *) {}

   virtual TRecorder::ERecorderState GetState() const = 0;

   ClassDef(TRecorderState, 0) // Abstract class that defines interface for a state of recorder
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  TRecorderReplaying                                                  //
//                                                                      //
//  Represents state of TRecorder when replaying previously recorded    //
//  events.                                                             //
//                                                                      //
//  Not intended to be used by a user directly.                         //
//  [Replaying.JPG]                                                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class TRecorderReplaying : public TRecorderState
{
private:
   virtual  ~TRecorderReplaying();
   Bool_t   PrepareNextEvent();
   Bool_t   RemapWindowReferences();
   Bool_t   CanOverlap();

   Bool_t   FilterEvent(TRecGuiEvent *e);

   TRecorder  *fRecorder;  // Reference to recorder (owner of this state) is kept in order to switch
                           // recorder to INACTIVE state after replaying is finished

   TFile      *fFile;      // ROOT file which the recorded events are being read from


   TCanvas    *fCanv;      // Used to record the previous canvases


   TTimer     *fTimer;     // Timer used for replaying

   TTree      *fWinTree;   // TTree with recorded windows (=registered during recording)
   TTree      *fGuiTree;   // TTree with recorded GUI events
   TTree      *fCmdTree;   // TTree with recorded commandline events
   TTree      *fExtraTree; // TTree with recorded extra events (PaveLabels and Texts)

   ULong64_t       fWin;            // Window ID being currenty mapped
   TRecGuiEvent   *fGuiEvent;       // GUI event being currently replayed
   TRecCmdEvent   *fCmdEvent;       // Commandline event being currently replayed
   TRecExtraEvent *fExtraEvent;     // Extra event being currently replayed

   Int_t       fRegWinCounter;      // Counter of registered windows when replaying
   Int_t       fGuiTreeCounter;     // Counter of GUI events that have been replayed
   Int_t       fCmdTreeCounter;     // Counter of commandline events that have been replayed
   Int_t       fExtraTreeCounter;   // Counter of extra events that have been replayed

   Int_t       fWinTreeEntries;     // Number of registered windows during _recording_

   TMutex      *fMutex;

   TList      *fWindowList;         // List of TRecWinPair objects. Mapping of window IDs is stored here.

   TRecEvent  *fNextEvent;          // The next event that is going to be replayed (GUI event or commandline)

   TTime       fPreviousEventTime;  // Execution time of the previously replayed event.
                                    // It is used for computing time difference between two events.

   Bool_t      fWaitingForWindow;   // Signalizes that we wait for a window to be registered in order
                                    // to replay the next event fNextEvent.
                                    // Registraion of windows can last different time when recording and replaying.
                                    // If there is an event ready to be replayed but the corresponding windows has not been yet
                                    // registered, we wait (postopone fNextEvent) until it is registered.

   Bool_t      fEventReplayed;      // Signalizes that the last event sent to the replaying has been already replayed.
                                    // Sometimes an execution of an event can take more time than during recording.
                                    // This ensures that the next event is sent to replaying AFTER
                                    // the replaying of the previous one finishes and not earlier.
                                    // Exceptions: ButtonPress and ButtonRelease events (See TRecorderReplaying::CanBeOverlapped)

   Bool_t      fShowMouseCursor;    // Specifies if mouse cursor should be also replayed

   Bool_t      fFilterStatusBar;    // Special flag to filter status bar element

protected:
   friend class TRecorderInactive;
   friend class TRecorderPaused;

   TRecorderReplaying(const char *filename);
   Bool_t     Initialize(TRecorder *r, Bool_t showMouseCursor, TRecorder::EReplayModes mode);

public:
   virtual TRecorder::ERecorderState GetState() const { return TRecorder::kReplaying; }

   virtual void   Pause(TRecorder *r);
   virtual void   Continue();
   virtual void   ReplayStop(TRecorder *r);

   void           RegisterWindow(Window_t w);   //SLOT
   void           ReplayRealtime();             //SLOT

   ClassDef(TRecorderReplaying, 0) // Represents state of TRecorder when replaying
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  TRecorderRecording                                                  //
//                                                                      //
//  Represents state of TRecorder when recording events.                //
//                                                                      //
//  Not intended to be used by a user directly.                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class TRecorderRecording: public TRecorderState
{
private:
   virtual ~TRecorderRecording();
   Bool_t  IsFiltered(Window_t id);
   void    SetTypeOfConfigureNotify(Event_t *e);
   void    CopyEvent(Event_t *e, Window_t wid);

   TRecorder          *fRecorder;         // Reference to recorder (owner of this state) is kept in order to switch
                                          // recorder back to INACTIVE state after recording is finished

   TFile              *fFile;             // ROOT file to store recorded events in
   TTimer             *fTimer;            // Timer used for recording
   TTimer             *fMouseTimer;       // Timer used for recording mouse position
   Long64_t            fBeginPave;        // TLatex/TPaveLabel edition starting time

   TTree              *fWinTree;          // TTree with registered windows
   TTree              *fGuiTree;          // TTree with recorded GUI events
   TTree              *fCmdTree;          // TTree with recorded commandline events
   TTree              *fExtraTree;        // TTree with recorded extra events (PaveLabels and Texts)

   ULong64_t           fWin;              // The newest registered window to be stored in TTree
   TRecGuiEvent       *fGuiEvent;         // The newest GUI event to be stored in TTree
   TRecCmdEvent       *fCmdEvent;         // The newest commandline event to be stored in TTree
   TRecExtraEvent     *fExtraEvent;       // The newest extra event to be stored in TTree

   Bool_t              fCmdEventPending;  // Indication if there is a still pending commandline event that should be stored.
                                          // Commandline events are stored with 1 event delay to ensure skipping
                                          // the last event 'TRecorder::Stop' that is not supposed to be recorded

   Int_t               fRegWinCounter;    // Counter of registered ROOT windows.
                                          // It is increased every time when a new window is registered

   Int_t               fFilteredIdsCount; // Only when GUI for recorder is used: Count of windows in GUI recorder
   Window_t           *fFilteredIds;      // Only when GUI for recorer is used: IDs of windows that creates that GUI.
                                          // Events for GUI recorder are not recorded.
   Bool_t              fFilterEventPave;  // Special flag to filter events during the pave recording

protected:
   friend class TRecorderInactive;
   TRecorderRecording(TRecorder *r, const char *filename, Option_t *option, Window_t *w, Int_t winCount);

   Bool_t StartRecording();

public:
   virtual TRecorder::ERecorderState GetState() const { return TRecorder::kRecording; }

   virtual void Stop(TRecorder *r, Bool_t guiCommand);

   void  RegisterWindow(Window_t w);               //SLOT
   void  RecordCmdEvent(const char *line);         //SLOT
   void  RecordGuiEvent(Event_t *e, Window_t wid); //SLOT
   void  RecordGuiBldEvent(Event_t *e);            //SLOT
   void  RecordGuiCNEvent(Event_t *e);             //SLOT
   void  RecordMousePosition();
   void  RecordPave(const TObject *obj);           //SLOT
   void  RecordText(const TObject *obj);           //SLOT
   void  FilterEventPave();                        //SLOT
   void  StartEditing();                           //SLOT

   void  RecordExtraEvent(TString line, TTime extTime);

   ClassDef(TRecorderRecording, 0) // Represents state of TRecorder when recording events
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  TRecorderInactive                                                   //
//                                                                      //
//  Represents state of TRecorder just after its creation.              //
//  Nor recording neither replaying is being executed in this state.    //
//                                                                      //
//  Not intended to be used by a user directly.                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class TRecorderInactive : public TRecorderState
{

private:
   TSeqCollection *fCollect;

public:
   virtual        ~TRecorderInactive() {}
   TRecorderInactive() : fCollect(0) {}

   virtual void   ListCmd(const char *filename);
   virtual void   ListGui(const char *filename);

   virtual void   Start(TRecorder *r, const char *filename, Option_t *option, Window_t *w = 0, Int_t winCount = 0);
   virtual Bool_t Replay(TRecorder *r, const char *filename, Bool_t showMouseCursor, TRecorder::EReplayModes mode);

   virtual TRecorder::ERecorderState GetState() const { return TRecorder::kInactive; }

   static void    DumpRootEvent(TRecGuiEvent *e, Int_t n);
   static long    DisplayValid(Long_t n) { return ( n < 0 ? -1 : n); }

   void PrevCanvases(const char *filename, Option_t *option);

   ClassDef(TRecorderInactive, 0) // Represents state of TRecorder after its creation
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  TRecorderPaused                                                     //
//                                                                      //
//  Represents state of TRecorder when replaying was paused             //
//  by a user.                                                          //
//  The paused replaying is remembered and after Resume call can        //
//  be continued again.                                                 //
//                                                                      //
//  Not intended to be used by a user directly.                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class TRecorderPaused: public TRecorderState
{
private:
   virtual ~TRecorderPaused() {}

   TRecorderReplaying       *fReplayingState;      // Replaying that is paused

protected:
   friend class TRecorderReplaying;
   TRecorderPaused(TRecorderReplaying *state);

public:
   virtual TRecorder::ERecorderState GetState() const { return TRecorder::kPaused; }

   virtual void Resume(TRecorder *r);
   virtual void ReplayStop(TRecorder *r);

   ClassDef(TRecorderPaused, 0) // Represents state of TRecorder when paused
};


//////////////////////////////////////////////////////////////////////////
//                                                                      //
//  TGRecorder                                                          //
//                                                                      //
//  Provides GUI for TRecorder class.                                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class TGRecorder : public TGMainFrame
{
private:
   TRecorder          *fRecorder;          // Recorder

   TGPictureButton    *fStartStop;         // Button for start and stop of recording
   TGPictureButton    *fReplay;            // Button for start of replaying

   TGLabel            *fStatus;            // Label with actual status
   TGLabel            *fTimeDisplay;       // Label with time counter
   TGCheckButton      *fCursorCheckBox;    // Check box "Show mouse cursor" for replaying

   TTimer             *fTimer;             // Timer for handling GUI of recorder
   time_t              fStart, fElapsed;   // playing/recording time

   static const Int_t  fgWidgetsCount = 12;            // Number of windows in GUI recorder
   Window_t            fFilteredIds[fgWidgetsCount];   // IDs of these windows in GUI recorder

   void                SetDefault();

public:
   TGRecorder(const TGWindow *p = 0, UInt_t w = 230, UInt_t h = 150);
   virtual ~TGRecorder();

   void StartStop();
   void Update();
   void Replay();

   ClassDef(TGRecorder,0) // GUI class of the event recorder.
};

#endif // ROOT_TRecorder
