// @(#)root/proofx:$Id$
// Author: Gerardo Ganis  12/12/2005

/*************************************************************************
 * Copyright (C) 1995-2005, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TProofServLite
\ingroup proofkernel

Version of the PROOF worker server for local running. The client starts
directly the desired number of these workers; the master and daemons are
eliminated, optimizing the number of messages exchanged and created / destroyed.

*/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TProofServLite                                                       //
//                                                                      //
// TProofServLite is the version of the PROOF worker server for local   //
// running. The client starts directly the desired number of these      //
// workers; the master and daemons are eliminated, optimizing the number//
// of messages exchanged and created / destroyed.                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "RConfigure.h"
#include <ROOT/RConfig.hxx>
#include "Riostream.h"

#ifdef WIN32
   #include <io.h>
   typedef long off_t;
   #include <snprintf.h>
#else
#include <netinet/in.h>
#endif
#include <sys/types.h>
#include <cstdlib>

#include "TProofServLite.h"
#include "TEnv.h"
#include "TError.h"
#include "TException.h"
#include "THashList.h"
#include "TInterpreter.h"
#include "TMessage.h"
#include "TProofDebug.h"
#include "TProof.h"
#include "TVirtualProofPlayer.h"
#include "TProofQueryResult.h"
#include "TRegexp.h"
#include "TClass.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TPluginManager.h"
#include "TSocket.h"
#include "TTimeStamp.h"
#include "compiledata.h"

using namespace std;

// debug hook
static volatile Int_t gProofServDebug = 1;

//----- Interrupt signal handler -----------------------------------------------
////////////////////////////////////////////////////////////////////////////////

class TProofServLiteInterruptHandler : public TSignalHandler {
   TProofServLite  *fServ;
public:
   TProofServLiteInterruptHandler(TProofServLite *s)
      : TSignalHandler(kSigUrgent, kFALSE) { fServ = s; }
   Bool_t  Notify();
};

////////////////////////////////////////////////////////////////////////////////
/// Handle urgent data

Bool_t TProofServLiteInterruptHandler::Notify()
{
   fServ->HandleUrgentData();
   if (TROOT::Initialized()) {
      Throw(GetSignal());
   }
   return kTRUE;
}

//----- SigPipe signal handler -------------------------------------------------
////////////////////////////////////////////////////////////////////////////////

class TProofServLiteSigPipeHandler : public TSignalHandler {
   TProofServLite  *fServ;
public:
   TProofServLiteSigPipeHandler(TProofServLite *s) : TSignalHandler(kSigPipe, kFALSE)
      { fServ = s; }
   Bool_t  Notify();
};

////////////////////////////////////////////////////////////////////////////////
/// Handle sig pipe

Bool_t TProofServLiteSigPipeHandler::Notify()
{
   fServ->HandleSigPipe();
   return kTRUE;
}

//----- Termination signal handler ---------------------------------------------
////////////////////////////////////////////////////////////////////////////////

class TProofServLiteTerminationHandler : public TSignalHandler {
   TProofServLite  *fServ;
public:
   TProofServLiteTerminationHandler(TProofServLite *s)
      : TSignalHandler(kSigTermination, kFALSE) { fServ = s; }
   Bool_t  Notify();
};

////////////////////////////////////////////////////////////////////////////////
/// Handle termination

Bool_t TProofServLiteTerminationHandler::Notify()
{
   Printf("TProofServLiteTerminationHandler::Notify: wake up!");

   fServ->HandleTermination();
   return kTRUE;
}

//----- Seg violation signal handler ---------------------------------------------
////////////////////////////////////////////////////////////////////////////////

class TProofServLiteSegViolationHandler : public TSignalHandler {
   TProofServLite  *fServ;
public:
   TProofServLiteSegViolationHandler(TProofServLite *s)
      : TSignalHandler(kSigSegmentationViolation, kFALSE) { fServ = s; }
   Bool_t  Notify();
};

////////////////////////////////////////////////////////////////////////////////
/// Handle seg violation

Bool_t TProofServLiteSegViolationHandler::Notify()
{
   Printf("**** ");
   Printf("**** Segmentation violation: terminating ****");
   Printf("**** ");
   fServ->HandleTermination();
   return kTRUE;
}

//----- Input handler for messages from parent or master -----------------------
////////////////////////////////////////////////////////////////////////////////

class TProofServLiteInputHandler : public TFileHandler {
   TProofServLite  *fServ;
public:
   TProofServLiteInputHandler(TProofServLite *s, Int_t fd) : TFileHandler(fd, 1)
      { fServ = s; }
   Bool_t Notify();
   Bool_t ReadNotify() { return Notify(); }
};

////////////////////////////////////////////////////////////////////////////////
/// Handle input on the socket

Bool_t TProofServLiteInputHandler::Notify()
{
   fServ->HandleSocketInput();

   return kTRUE;
}

ClassImp(TProofServLite);

// Hook to the constructor. This is needed to avoid using the plugin manager
// which may create problems in multi-threaded environments.
extern "C" {
   TApplication *GetTProofServLite(Int_t *argc, char **argv, FILE *flog)
   { return new TProofServLite(argc, argv, flog); }
}

////////////////////////////////////////////////////////////////////////////////
/// Main constructor

TProofServLite::TProofServLite(Int_t *argc, char **argv, FILE *flog)
            : TProofServ(argc, argv, flog)
{
   fInterruptHandler = 0;
   fTerminated = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Finalize the server setup. If master, create the TProof instance to talk
/// the worker or submaster nodes.
/// Return 0 on success, -1 on error

Int_t TProofServLite::CreateServer()
{
   if (gProofDebugLevel > 0)
      Info("CreateServer", "starting server creation");

   // Get file descriptor for log file
   if (fLogFile) {
      // Use the file already open by pmain
      if ((fLogFileDes = fileno(fLogFile)) < 0) {
         Error("CreateServer", "resolving the log file description number");
         return -1;
      }
   }

   // Get socket to be used to call back our xpd
   fSockPath = gEnv->GetValue("ProofServ.OpenSock", "");
   if (fSockPath.Length() <= 0) {
      Error("CreateServer", "Socket setup by xpd undefined");
      return -1;
   }
   TString entity = gEnv->GetValue("ProofServ.Entity", "");
   if (entity.Length() > 0)
      fSockPath.Insert(0,TString::Format("%s/", entity.Data()));

   // Call back the client
   fSocket = new TSocket(fSockPath);
   if (!fSocket || !(fSocket->IsValid())) {
      Error("CreateServer", "Failed to open connection to the client");
      return -1;
   }

   // Send our ordinal, to allow the client to identify us
   TMessage msg;
   msg << fOrdinal;
   fSocket->Send(msg);

   // Get socket descriptor
   Int_t sock = fSocket->GetDescriptor();

   // Install interrupt and message input handlers
   fInterruptHandler = new TProofServLiteInterruptHandler(this);
   gSystem->AddSignalHandler(fInterruptHandler);
   gSystem->AddFileHandler(new TProofServLiteInputHandler(this, sock));

   // Wait (loop) in worker node to allow debugger to connect
   if (gEnv->GetValue("Proof.GdbHook",0) == 2) {
      while (gProofServDebug)
         ;
   }

   if (gProofDebugLevel > 0)
      Info("CreateServer", "Service: %s, ConfDir: %s, IsMaster: %d",
           fService.Data(), fConfDir.Data(), (Int_t)fMasterServ);

   if (Setup() == -1) {
      // Setup failure
      Terminate(0);
      SendLogFile();
      return -1;
   }

   if (!fLogFile) {
      RedirectOutput();
      // If for some reason we failed setting a redirection file for the logs
      // we cannot continue
      if (!fLogFile || (fLogFileDes = fileno(fLogFile)) < 0) {
         Terminate(0);
         SendLogFile(-98);
         return -1;
      }
   }

   // Everybody expects std::iostream to be available, so load it...
   ProcessLine("#include <iostream>", kTRUE);
   ProcessLine("#include <string>",kTRUE); // for std::string std::iostream.

   // Load user functions
   const char *logon;
   logon = gEnv->GetValue("Proof.Load", (char *)0);
   if (logon) {
      char *mac = gSystem->Which(TROOT::GetMacroPath(), logon, kReadPermission);
      if (mac)
         ProcessLine(TString::Format(".L %s", logon), kTRUE);
      delete [] mac;
   }

   // Execute logon macro
   logon = gEnv->GetValue("Proof.Logon", (char *)0);
   if (logon && !NoLogOpt()) {
      char *mac = gSystem->Which(TROOT::GetMacroPath(), logon, kReadPermission);
      if (mac)
         ProcessFile(logon);
      delete [] mac;
   }

   // Save current interpreter context
   gInterpreter->SaveContext();
   gInterpreter->SaveGlobalsContext();

   // Avoid spurious messages at first action
   FlushLogFile();

   // Done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Cleanup. Not really necessary since after this dtor there is no
/// live anyway.

TProofServLite::~TProofServLite()
{
   delete fSocket;
}

////////////////////////////////////////////////////////////////////////////////
/// Called when the client is not alive anymore; terminate the session.

void TProofServLite::HandleSigPipe()
{
   Terminate(0);  // will not return from here....
}

////////////////////////////////////////////////////////////////////////////////
/// Called when the client is not alive anymore; terminate the session.

void TProofServLite::HandleTermination()
{
   Terminate(0);  // will not return from here....
}

////////////////////////////////////////////////////////////////////////////////
/// Print the ProofServ logo on standard output.
/// Return 0 on success, -1 on error

Int_t TProofServLite::Setup()
{
   char str[512];

   if (IsMaster()) {
      snprintf(str, 512, "**** Welcome to the PROOF server @ %s ****", gSystem->HostName());
   } else {
      snprintf(str, 512, "**** PROOF worker server @ %s started ****", gSystem->HostName());
   }

   if (fSocket->Send(str) != 1+static_cast<Int_t>(strlen(str))) {
      Error("Setup", "failed to send proof server startup message");
      return -1;
   }

   // Get client protocol
   if ((fProtocol = gEnv->GetValue("ProofServ.ClientVersion", -1)) < 0) {
      Error("Setup", "remote proof protocol missing");
      return -1;
   }

   // The local user
   UserGroup_t *pw = gSystem->GetUserInfo();
   if (pw) {
      fUser = pw->fUser;
      delete pw;
   }

   // Work dir and ...
   fWorkDir = gEnv->GetValue("ProofServ.Sandbox", TString::Format("~/%s", kPROOF_WorkDir));
   Info("Setup", "fWorkDir: %s", fWorkDir.Data());

   // Get Session tags
   fTopSessionTag = gEnv->GetValue("ProofServ.SessionTag", "-1");
   fSessionTag.Form("%s-%s-%ld-%d", fOrdinal.Data(), gSystem->HostName(),
                                    (Long_t)TTimeStamp().GetSec(), gSystem->GetPid());
   if (gProofDebugLevel > 0)
      Info("Setup", "session tag is %s", fSessionTag.Data());
   if (fTopSessionTag.IsNull()) fTopSessionTag = fSessionTag;

   // Send session tag to client
   TMessage m(kPROOF_SESSIONTAG);
   m << fSessionTag;
   fSocket->Send(m);

   // Get Session dir (sandbox)
   if ((fSessionDir = gEnv->GetValue("ProofServ.SessionDir", "-1")) == "-1") {
      Error("Setup", "Session dir missing");
      return -1;
   }

   // Link the session tag to the log file
   if (gSystem->Getenv("ROOTPROOFLOGFILE")) {
      TString logfile = gSystem->Getenv("ROOTPROOFLOGFILE");
      Int_t iord = logfile.Index(TString::Format("-%s", fOrdinal.Data()));
      if (iord != kNPOS) logfile.Remove(iord);
      logfile += TString::Format("-%s.log", fSessionTag.Data());
      gSystem->Symlink(gSystem->Getenv("ROOTPROOFLOGFILE"), logfile);
   }

   // Goto to the main PROOF working directory
   gSystem->ExpandPathName(fWorkDir);
   if (gProofDebugLevel > 0)
      Info("Setup", "working directory set to %s", fWorkDir.Data());

   // Common setup
   if (SetupCommon() != 0) {
      Error("Setup", "common setup failed");
      return -1;
   }

   // Check every two hours if client is still alive
   fSocket->SetOption(kKeepAlive, 1);

   // Install SigPipe handler to handle kKeepAlive failure
   gSystem->AddSignalHandler(new TProofServLiteSigPipeHandler(this));

   // Install Termination handler
   gSystem->AddSignalHandler(new TProofServLiteTerminationHandler(this));

   // Install seg violation handler
   gSystem->AddSignalHandler(new TProofServLiteSegViolationHandler(this));

   // Done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Terminate the proof server.

void TProofServLite::Terminate(Int_t status)
{
   if (fTerminated)
      // Avoid doubling the exit operations
      exit(1);
   fTerminated = kTRUE;

   // Notify
   Info("Terminate", "starting session termination operations ...");

   // Cleanup session directory
   if (status == 0) {
      // make sure we remain in a "connected" directory
      gSystem->ChangeDirectory("/");
      // needed in case fSessionDir is on NFS ?!
      gSystem->MakeDirectory(fSessionDir+"/.delete");
      gSystem->Exec(TString::Format("%s %s", kRM, fSessionDir.Data()));
   }

   // Cleanup data directory if empty
   if (!fDataDir.IsNull() && !gSystem->AccessPathName(fDataDir, kWritePermission)) {
     if (UnlinkDataDir(fDataDir))
        Info("Terminate", "data directory '%s' has been removed", fDataDir.Data());
   }

   // Remove input and signal handlers to avoid spurious "signals"
   // for closing activities executed upon exit()
   gSystem->RemoveSignalHandler(fInterruptHandler);

   // Stop processing events (set a flag to exit the event loop)
   gSystem->ExitLoop();

   // Notify
   Printf("Terminate: termination operations ended: quitting!");
}

////////////////////////////////////////////////////////////////////////////////
/// Cloning itself via fork.

void TProofServLite::HandleFork(TMessage *mess)
{
   if (!mess) {
      Error("HandleFork", "empty message!");
      return;
   }

   // Extract the ordinals of the clones
   TString clones;
   (*mess) >> clones;
   PDB(kGlobal, 1)
      Info("HandleFork", "cloning to %s", clones.Data());

   TString clone;
   Int_t from = 0;
   while (clones.Tokenize(clone, from, " ")) {

      Int_t rc = 0;
      // Fork
      if ((rc = Fork()) < 0) {
         Error("HandleFork", "failed to fork %s", clone.Data());
         return;
      }

      // If the child, finalize the setup and return
      if (rc == 0) {
         SetupOnFork(clone.Data());
         return;
      }
   }

   // Done
   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Finalize the server setup afetr forking.
/// Return 0 on success, -1 on error

Int_t TProofServLite::SetupOnFork(const char *ord)
{
   if (gProofDebugLevel > 0)
      Info("SetupOnFork", "finalizing setup of %s", ord);

   // Set the ordinal
   fOrdinal = ord;
   TString sord;
   sord.Form("-%s", fOrdinal.Data());

   // Close the current log file
   if (fLogFile) {
      fclose(fLogFile);
      fLogFileDes = -1;
   }

   TString sdir = gSystem->GetDirName(fSessionDir.Data());
   RedirectOutput(sdir.Data(), "a");
   // If for some reason we failed setting a redirection file for the logs
   // we cannot continue
   if (!fLogFile || (fLogFileDes = fileno(fLogFile)) < 0) {
      Terminate(0);
      return -1;
   }
   FlushLogFile();

   // Eliminate existing symlink
   void *dirp = gSystem->OpenDirectory(sdir);
   if (dirp) {
      TString ent;
      const char *e = 0;
      while ((e = gSystem->GetDirEntry(dirp))) {
         ent.Form("%s/%s", sdir.Data(), e);
         FileStat_t st;
         if (gSystem->GetPathInfo(ent.Data(), st) == 0) {
            if (st.fIsLink && ent.Contains(sord)) {
               PDB(kGlobal, 1)
                  Info("SetupOnFork","unlinking: %s", ent.Data());
               gSystem->Unlink(ent);
            }
         }
      }
      gSystem->FreeDirectory(dirp);
   }

   // The session tag
   fSessionTag.Form("%s-%d-%d", gSystem->HostName(), (int)time(0), gSystem->GetPid());

   // Create new symlink
   TString logfile = gSystem->Getenv("ROOTPROOFLOGFILE");
   logfile.ReplaceAll("-0.0", sord.Data());
   gSystem->Setenv("ROOTPROOFLOGFILE", logfile);
   Int_t iord = logfile.Index(sord.Data());
   if (iord != kNPOS) logfile.Remove(iord + sord.Length());
   logfile += TString::Format("-%s.log", fSessionTag.Data());
   gSystem->Symlink(gSystem->Getenv("ROOTPROOFLOGFILE"), logfile);

   // Get socket to be used to call back our xpd
   fSockPath = gEnv->GetValue("ProofServ.OpenSock", "");
   if (fSockPath.Length() <= 0) {
      Error("CreateServer", "Socket setup by xpd undefined");
      return -1;
   }
   TString entity = gEnv->GetValue("ProofServ.Entity", "");
   if (entity.Length() > 0)
      fSockPath.Insert(0, TString::Format("%s/", entity.Data()));

   // Call back the client
   fSocket = new TSocket(fSockPath);
   if (!fSocket || !(fSocket->IsValid())) {
      Error("CreateServer", "Failed to open connection to the client");
      return -1;
   }

   // Send our ordinal, to allow the client to identify us
   TMessage msg;
   msg << fOrdinal;
   fSocket->Send(msg);

   // Get socket descriptor
   Int_t sock = fSocket->GetDescriptor();

   // Install interrupt and message input handlers
   fInterruptHandler = new TProofServLiteInterruptHandler(this);
   gSystem->AddSignalHandler(fInterruptHandler);
   gSystem->AddFileHandler(new TProofServLiteInputHandler(this, sock));

   // Wait (loop) in worker node to allow debugger to connect
   if (gEnv->GetValue("Proof.GdbHook",0) == 2) {
      while (gProofServDebug)
         ;
   }

   if (gProofDebugLevel > 0)
      Info("SetupOnFork", "Service: %s, ConfDir: %s, IsMaster: %d",
           fService.Data(), fConfDir.Data(), (Int_t)fMasterServ);

   if (Setup() == -1) {
      // Setup failure
      Terminate(0);
      SendLogFile();
      return -1;
   }

   // Disallow the interpretation of Rtypes.h, TError.h and TGenericClassInfo.h
   ProcessLine("#define ROOT_Rtypes 0", kTRUE);
   ProcessLine("#define ROOT_TError 0", kTRUE);
   ProcessLine("#define ROOT_TGenericClassInfo 0", kTRUE);

   // Save current interpreter context
   gInterpreter->SaveContext();
   gInterpreter->SaveGlobalsContext();

   // Done
   return 0;
}
