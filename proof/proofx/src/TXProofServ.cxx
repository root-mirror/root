// @(#)root/proofx:$Id$
// Author: Gerardo Ganis  12/12/2005

/*************************************************************************
 * Copyright (C) 1995-2005, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TXProofServ
\ingroup proofx

This class implements the XProofD version of TProofServ, with respect to which it differs only
for the underlying connection technology.

*/

#include "RConfigure.h"
#include <ROOT/RConfig.hxx>
#include "Riostream.h"

#ifdef WIN32
   #include <io.h>
   typedef long off_t;
#endif
#include <sys/types.h>
#include <netinet/in.h>
#include <utime.h>

#include "TXProofServ.h"
#include "TObjString.h"
#include "TEnv.h"
#include "TError.h"
#include "TException.h"
#include "THashList.h"
#include "TInterpreter.h"
#include "TParameter.h"
#include "TProofDebug.h"
#include "TProof.h"
#include "TVirtualProofPlayer.h"
#include "TQueryResultManager.h"
#include "TRegexp.h"
#include "TClass.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TPluginManager.h"
#include "TXSocketHandler.h"
#include "TXUnixSocket.h"
#include "compiledata.h"
#include "TProofNodeInfo.h"
#include "XProofProtocol.h"
#include "snprintf.h"

#include <XrdClient/XrdClientConst.hh>
#include <XrdClient/XrdClientEnv.hh>


// debug hook
static volatile Int_t gProofServDebug = 1;

//----- SigPipe signal handler -------------------------------------------------
////////////////////////////////////////////////////////////////////////////////

class TXProofServSigPipeHandler : public TSignalHandler {
   TXProofServ  *fServ;
public:
   TXProofServSigPipeHandler(TXProofServ *s) : TSignalHandler(kSigInterrupt, kFALSE)
      { fServ = s; }
   Bool_t  Notify();
};

////////////////////////////////////////////////////////////////////////////////

Bool_t TXProofServSigPipeHandler::Notify()
{
   fServ->HandleSigPipe();
   return kTRUE;
}

//----- Termination signal handler ---------------------------------------------
////////////////////////////////////////////////////////////////////////////////

class TXProofServTerminationHandler : public TSignalHandler {
   TXProofServ  *fServ;
public:
   TXProofServTerminationHandler(TXProofServ *s)
      : TSignalHandler(kSigTermination, kFALSE) { fServ = s; }
   Bool_t  Notify();
};

////////////////////////////////////////////////////////////////////////////////

Bool_t TXProofServTerminationHandler::Notify()
{
   Printf("Received SIGTERM: terminating");

   fServ->HandleTermination();
   return kTRUE;
}

//----- Seg violation signal handler ---------------------------------------------
////////////////////////////////////////////////////////////////////////////////

class TXProofServSegViolationHandler : public TSignalHandler {
   TXProofServ  *fServ;
public:
   TXProofServSegViolationHandler(TXProofServ *s)
      : TSignalHandler(kSigSegmentationViolation, kFALSE) { fServ = s; }
   Bool_t  Notify();
};

////////////////////////////////////////////////////////////////////////////////

Bool_t TXProofServSegViolationHandler::Notify()
{
   Printf("**** ");
   Printf("**** Segmentation violation: terminating ****");
   Printf("**** ");
   fServ->HandleTermination();
   return kTRUE;
}

//----- Input handler for messages from parent or master -----------------------
////////////////////////////////////////////////////////////////////////////////

class TXProofServInputHandler : public TFileHandler {
   TXProofServ  *fServ;
public:
   TXProofServInputHandler(TXProofServ *s, Int_t fd) : TFileHandler(fd, 1)
      { fServ = s; }
   Bool_t Notify();
   Bool_t ReadNotify() { return Notify(); }
};

////////////////////////////////////////////////////////////////////////////////

Bool_t TXProofServInputHandler::Notify()
{
   fServ->HandleSocketInput();
   // This request has been completed: remove the client ID from the pipe
   ((TXUnixSocket *) fServ->GetSocket())->RemoveClientID();
   return kTRUE;
}

ClassImp(TXProofServ);

// Hook to the constructor. This is needed to avoid using the plugin manager
// which may create problems in multi-threaded environments.
extern "C" {
   TApplication *GetTXProofServ(Int_t *argc, char **argv, FILE *flog)
   { return new TXProofServ(argc, argv, flog); }
}

////////////////////////////////////////////////////////////////////////////////
/// Main constructor

TXProofServ::TXProofServ(Int_t *argc, char **argv, FILE *flog)
            : TProofServ(argc, argv, flog)
{
   fInterruptHandler = 0;
   fInputHandler = 0;
   fTerminated = kFALSE;

   // TODO:
   //    Int_t useFIFO = 0;
/*   if (GetParameter(fProof->GetInputList(), "PROOF_UseFIFO", useFIFO) != 0) {
      if (useFIFO == 1)
         Info("", "enablig use of FIFO (if allowed by the server)");
      else
         Warning("", "unsupported strategy index (%d): ignore", strategy);
   }
*/
}

////////////////////////////////////////////////////////////////////////////////
/// Finalize the server setup. If master, create the TProof instance to talk
/// the worker or submaster nodes.
/// Return 0 on success, -1 on error

Int_t TXProofServ::CreateServer()
{
   Bool_t xtest = (Argc() > 3 && !strcmp(Argv(3), "test")) ? kTRUE : kFALSE;

   if (gProofDebugLevel > 0)
      Info("CreateServer", "starting%s server creation", (xtest ? " test" : ""));

   // Get file descriptor for log file
   if (fLogFile) {
      // Use the file already open by pmain
      if ((fLogFileDes = fileno(fLogFile)) < 0) {
         Error("CreateServer", "resolving the log file description number");
         return -1;
      }
      // Hide the session start-up logs unless we are in verbose mode
      if (gProofDebugLevel <= 0)
         lseek(fLogFileDes, (off_t) 0, SEEK_END);
   }

   // Global location string in TXSocket
   TXSocket::SetLocation((IsMaster()) ? "master" : "slave");

   // Set debug level in XrdClient
   EnvPutInt(NAME_DEBUG, gEnv->GetValue("XNet.Debug", 0));

   // Get socket to be used to call back our xpd
   if (xtest) {
      // test session, just send the protocol version on the open pipe
      // and exit
      if (!(fSockPath = gSystem->Getenv("ROOTOPENSOCK"))) {
         Error("CreateServer", "test: socket setup by xpd undefined");
         return -1;
      }
      Int_t fpw = (Int_t) strtol(fSockPath.Data(), 0, 10);
      int proto = htonl(kPROOF_Protocol);
      fSockPath = "";
      if (write(fpw, &proto, sizeof(proto)) != sizeof(proto)) {
         Error("CreateServer", "test: sending protocol number");
         return -1;
      }
      exit(0);
   } else {
      fSockPath = gEnv->GetValue("ProofServ.OpenSock", "");
      if (fSockPath.Length() <= 0) {
         Error("CreateServer", "socket setup by xpd undefined");
         return -1;
      }
      TString entity = gEnv->GetValue("ProofServ.Entity", "");
      if (entity.Length() > 0)
         fSockPath.Insert(0,Form("%s/", entity.Data()));
   }

   // Get open socket descriptor, if any
   Int_t sockfd = -1;
   const char *opensock = gSystem->Getenv("ROOTOPENSOCK");
   if (opensock && strlen(opensock) > 0) {
      TSystem::ResetErrno();
      sockfd = (Int_t) strtol(opensock, 0, 10);
      if (TSystem::GetErrno() == ERANGE) {
         sockfd = -1;
         Warning("CreateServer", "socket descriptor: wrong conversion from '%s'", opensock);
      }
      if (sockfd > 0 && gProofDebugLevel > 0)
         Info("CreateServer", "using open connection (descriptor %d)", sockfd);
   }

   // Get the sessions ID
   Int_t psid = gEnv->GetValue("ProofServ.SessionID", -1);
   if (psid < 0) {
     Error("CreateServer", "Session ID undefined");
     return -1;
   }

   // Call back the server
   fSocket = new TXUnixSocket(fSockPath, psid, -1, this, sockfd);
   if (!fSocket || !(fSocket->IsValid())) {
      Error("CreateServer", "Failed to open connection to XrdProofd coordinator");
      return -1;
   }
   // Set compression level, if any
   fSocket->SetCompressionSettings(fCompressMsg);

   // Set the title for debugging
   TString tgt("client");
   if (fOrdinal != "0") {
      tgt = fOrdinal;
      if (tgt.Last('.') != kNPOS) tgt.Remove(tgt.Last('.'));
   }
   fSocket->SetTitle(tgt);

   // Set the this as reference of this socket
   ((TXSocket *)fSocket)->fReference = this;

   // Get socket descriptor
   Int_t sock = fSocket->GetDescriptor();

   // Install message input handlers
   fInputHandler =
      TXSocketHandler::GetSocketHandler(new TXProofServInputHandler(this, sock), fSocket);
   gSystem->AddFileHandler(fInputHandler);

   // Get the client ID
   Int_t cid = gEnv->GetValue("ProofServ.ClientID", -1);
   if (cid < 0) {
     Error("CreateServer", "Client ID undefined");
     SendLogFile();
     return -1;
   }
   ((TXSocket *)fSocket)->SetClientID(cid);

   // debug hooks
   if (IsMaster()) {
      // wait (loop) in master to allow debugger to connect
      if (gEnv->GetValue("Proof.GdbHook",0) == 1) {
         while (gProofServDebug)
            ;
      }
   } else {
      // wait (loop) in slave to allow debugger to connect
      if (gEnv->GetValue("Proof.GdbHook",0) == 2) {
         while (gProofServDebug)
            ;
      }
   }

   if (gProofDebugLevel > 0)
      Info("CreateServer", "Service: %s, ConfDir: %s, IsMaster: %d",
           fService.Data(), fConfDir.Data(), (Int_t)fMasterServ);

   if (Setup() == -1) {
      // Setup failure
      LogToMaster();
      SendLogFile();
      Terminate(0);
      return -1;
   }

   if (!fLogFile) {
      RedirectOutput();
      // If for some reason we failed setting a redirection file for the logs
      // we cannot continue
      if (!fLogFile || (fLogFileDes = fileno(fLogFile)) < 0) {
         LogToMaster();
         SendLogFile(-98);
         Terminate(0);
         return -1;
      }
   }

   // Send message of the day to the client
   if (IsMaster()) {
      if (CatMotd() == -1) {
         LogToMaster();
         SendLogFile(-99);
         Terminate(0);
         return -1;
      }
   }

   // Everybody expects iostream to be available, so load it...
   ProcessLine("#include <iostream>", kTRUE);
   ProcessLine("#include <string>",kTRUE); // for std::string iostream.

   // Load user functions
   const char *logon;
   logon = gEnv->GetValue("Proof.Load", (char *)0);
   if (logon) {
      char *mac = gSystem->Which(TROOT::GetMacroPath(), logon, kReadPermission);
      if (mac)
         ProcessLine(Form(".L %s", logon), kTRUE);
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

   // if master, start slave servers
   if (IsMaster()) {
      TString master;

      if (fConfFile.BeginsWith("lite:")) {
         master = "lite://";
      } else {
         master.Form("proof://%s@__master__", fUser.Data());

         // Add port, if defined
         Int_t port = gEnv->GetValue("ProofServ.XpdPort", -1);
         if (port > -1) {
            master += ":";
            master += port;
         }
      }

      // Make sure that parallel startup via threads is not active
      // (it is broken for xpd because of the locks on gInterpreterMutex)
      gEnv->SetValue("Proof.ParallelStartup", 0);

      // Get plugin manager to load appropriate TProof from
      TPluginManager *pm = gROOT->GetPluginManager();
      if (!pm) {
         Error("CreateServer", "no plugin manager found");
         SendLogFile(-99);
         Terminate(0);
         return -1;
      }

      // Find the appropriate handler
      TPluginHandler *h = pm->FindHandler("TProof", fConfFile);
      if (!h) {
         Error("CreateServer", "no plugin found for TProof with a"
                             " config file of '%s'", fConfFile.Data());
         SendLogFile(-99);
         Terminate(0);
         return -1;
      }

      // load the plugin
      if (h->LoadPlugin() == -1) {
         Error("CreateServer", "plugin for TProof could not be loaded");
         SendLogFile(-99);
         Terminate(0);
         return -1;
      }

      // Make instance of TProof
      if (fConfFile.BeginsWith("lite:")) {
         // Remove input and signal handlers to avoid spurious "signals"
         // during startup
         gSystem->RemoveFileHandler(fInputHandler);
         fProof = reinterpret_cast<TProof*>(h->ExecPlugin(6, master.Data(),
                                                            0, 0,
                                                            fLogLevel,
                                                            fSessionDir.Data(), 0));
         // Re-enable input and signal handlers
         gSystem->AddFileHandler(fInputHandler);
      } else {
         fProof = reinterpret_cast<TProof*>(h->ExecPlugin(5, master.Data(),
                                                            fConfFile.Data(),
                                                            fConfDir.Data(),
                                                            fLogLevel,
                                                            fTopSessionTag.Data()));
      }

      // Save worker info
      if (fProof) fProof->SaveWorkerInfo();

      if (!fProof || (fProof && !fProof->IsValid())) {
         Error("CreateServer", "plugin for TProof could not be executed");
         FlushLogFile();
         delete fProof;
         fProof = 0;
         SendLogFile(-99);
         Terminate(0);
         return -1;
      }
      // Find out if we are a master in direct contact only with workers
      fEndMaster = fProof->IsEndMaster();

      SendLogFile();
   }

   // Setup the shutdown timer
   if (!fShutdownTimer) {
      // Check activity on socket every 5 mins
      fShutdownTimer = new TShutdownTimer(this, 300000);
      fShutdownTimer->Start(-1, kFALSE);
   }

   // Check if schema evolution is effective: clients running versions <=17 do not
   // support that: send a warning message
   if (fProtocol <= 17) {
      TString msg;
      msg.Form("Warning: client version is too old: automatic schema evolution is ineffective.\n"
               "         This may generate compatibility problems between streamed objects.\n"
               "         The advise is to move to ROOT >= 5.21/02 .");
      SendAsynMessage(msg.Data());
   }

   // Setup the idle timer
   if (IsMaster() && !fIdleTOTimer) {
      // Check activity on socket every 5 mins
      Int_t idle_to = gEnv->GetValue("ProofServ.IdleTimeout", -1);
      if (idle_to > 0) {
         fIdleTOTimer = new TIdleTOTimer(this, idle_to * 1000);
         fIdleTOTimer->Start(-1, kTRUE);
         if (gProofDebugLevel > 0)
            Info("CreateServer", " idle timer started (%d secs)", idle_to);
      } else if (gProofDebugLevel > 0) {
         Info("CreateServer", " idle timer not started (no idle timeout requested)");
      }
   }

   // Done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Cleanup. Not really necessary since after this dtor there is no
/// live anyway.

TXProofServ::~TXProofServ()
{
   delete fSocket;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle high priority data sent by the master or client.

void TXProofServ::HandleUrgentData()
{
   // Real-time notification of messages
   TProofServLogHandlerGuard hg(fLogFile, fSocket, "", fRealTimeLog);

   // Get interrupt
   Bool_t fw = kFALSE;
   Int_t iLev = ((TXSocket *)fSocket)->GetInterrupt(fw);
   if (iLev < 0) {
      Error("HandleUrgentData", "error receiving interrupt");
      return;
   }

   PDB(kGlobal, 2)
      Info("HandleUrgentData", "got interrupt: %d\n", iLev);

   if (fProof)
      fProof->SetActive();

   switch (iLev) {

      case TProof::kPing:
         PDB(kGlobal, 2)
            Info("HandleUrgentData", "*** Ping");

         // If master server, propagate interrupt to slaves
         if (fw && IsMaster()) {
            Int_t nbad = fProof->fActiveSlaves->GetSize() - fProof->Ping();
            if (nbad > 0) {
               Info("HandleUrgentData","%d slaves did not reply to ping",nbad);
            }
         }

         // Touch the admin path to show we are alive
         if (fAdminPath.IsNull()) {
            fAdminPath = gEnv->GetValue("ProofServ.AdminPath", "");
         }

         if (!fAdminPath.IsNull()) {
            if (!fAdminPath.EndsWith(".status")) {
               // Update file time stamps
               if (utime(fAdminPath.Data(), 0) != 0)
                  Info("HandleUrgentData", "problems touching path: %s", fAdminPath.Data());
               else
                  PDB(kGlobal, 2)
                     Info("HandleUrgentData", "touching path: %s", fAdminPath.Data());
            } else {
               // Update the status in the file
               //     0     idle
               //     1     running
               //     2     being terminated  (currently unused)
               //     3     queued
               //     4     idle timed-out
               Int_t uss_rc = UpdateSessionStatus(-1);
               if (uss_rc != 0)
                  Error("HandleUrgentData", "problems updating status path: %s (errno: %d)", fAdminPath.Data(), -uss_rc);
            }
         } else {
            Info("HandleUrgentData", "admin path undefined");
         }

         break;

      case TProof::kHardInterrupt:
         Info("HandleUrgentData", "*** Hard Interrupt");

         // If master server, propagate interrupt to slaves
         if (fw && IsMaster())
            fProof->Interrupt(TProof::kHardInterrupt);

         // Flush input socket
         ((TXSocket *)fSocket)->Flush();

         if (IsMaster())
            SendLogFile();

         break;

      case TProof::kSoftInterrupt:
         Info("HandleUrgentData", "Soft Interrupt");

         // If master server, propagate interrupt to slaves
         if (fw && IsMaster())
            fProof->Interrupt(TProof::kSoftInterrupt);

         Interrupt();

         if (IsMaster())
            SendLogFile();

         break;


      case TProof::kShutdownInterrupt:
         Info("HandleUrgentData", "Shutdown Interrupt");

         // When returning for here connection are closed
         HandleTermination();

         break;

      default:
         Error("HandleUrgentData", "unexpected type: %d", iLev);
         break;
   }


   if (fProof) fProof->SetActive(kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Called when the client is not alive anymore; terminate the session.

void TXProofServ::HandleSigPipe()
{
   // Real-time notification of messages

   Info("HandleSigPipe","got sigpipe ... do nothing");
}

////////////////////////////////////////////////////////////////////////////////
/// Called when the client is not alive anymore; terminate the session.

void TXProofServ::HandleTermination()
{
   // If master server, propagate interrupt to slaves
   // (shutdown interrupt send internally).
   if (IsMaster()) {

      // If not idle, try first to stop processing
      if (!fIdle) {
         // Remove pending requests
         fWaitingQueries->Delete();
         // Interrupt the current monitor
         fProof->InterruptCurrentMonitor();
         // Do not wait for ever, but al least 20 seconds
         Long_t timeout = gEnv->GetValue("Proof.ShutdownTimeout", 60);
         timeout = (timeout > 20) ? timeout : 20;
         // Processing will be aborted
         fProof->StopProcess(kTRUE, (Long_t) (timeout / 2));
         // Receive end-of-processing messages, but do not wait for ever
         fProof->Collect(TProof::kActive, timeout);
         // Still not idle
         if (!fIdle)
            Warning("HandleTermination","processing could not be stopped");
      }
      // Close the session
      if (fProof)
         fProof->Close("S");
   }

   Terminate(0);  // will not return from here....
}

////////////////////////////////////////////////////////////////////////////////
/// Print the ProofServ logo on standard output.
/// Return 0 on success, -1 on error

Int_t TXProofServ::Setup()
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
   fUser = gEnv->GetValue("ProofServ.Entity", "");
   if (fUser.Length() >= 0) {
      if (fUser.Contains(":"))
         fUser.Remove(fUser.Index(":"));
      if (fUser.Contains("@"))
         fUser.Remove(fUser.Index("@"));
   } else {
      UserGroup_t *pw = gSystem->GetUserInfo();
      if (pw) {
         fUser = pw->fUser;
         delete pw;
      }
   }

   // Work dir and ...
   if (IsMaster()) {
      TString cf = gEnv->GetValue("ProofServ.ProofConfFile", "");
      if (cf.Length() > 0)
         fConfFile = cf;
   }
   fWorkDir = gEnv->GetValue("ProofServ.Sandbox", Form("~/%s", kPROOF_WorkDir));

   // Get Session tag
   if ((fSessionTag = gEnv->GetValue("ProofServ.SessionTag", "-1")) == "-1") {
      Error("Setup", "Session tag missing");
      return -1;
   }
   // Get top session tag, i.e. the tag of the PROOF session
   if ((fTopSessionTag = gEnv->GetValue("ProofServ.TopSessionTag", "-1")) == "-1") {
      fTopSessionTag = "";
      // Try to extract it from log file path (for backward compatibility)
      if (gSystem->Getenv("ROOTPROOFLOGFILE")) {
         fTopSessionTag = gSystem->GetDirName(gSystem->Getenv("ROOTPROOFLOGFILE"));
         Ssiz_t lstl;
         if ((lstl = fTopSessionTag.Last('/')) != kNPOS) fTopSessionTag.Remove(0, lstl + 1);
         if (fTopSessionTag.BeginsWith("session-")) {
            fTopSessionTag.Remove(0, strlen("session-"));
         } else {
            fTopSessionTag = "";
         }
      }
      if (fTopSessionTag.IsNull()) {
         Error("Setup", "top session tag missing");
         return -1;
      }
   }

   // Make sure the process ID is in the tag
   TString spid = Form("-%d", gSystem->GetPid());
   if (!fSessionTag.EndsWith(spid)) {
      Int_t nd = 0;
      if ((nd = fSessionTag.CountChar('-')) >= 2) {
         Int_t id = fSessionTag.Index("-", fSessionTag.Index("-") + 1);
         if (id != kNPOS) fSessionTag.Remove(id);
      } else if (nd != 1) {
         Warning("Setup", "Wrong number of '-' in session tag: protocol error? %s", fSessionTag.Data());
      }
      // Add this process ID
      fSessionTag += spid;
   }
   if (gProofDebugLevel > 0)
      Info("Setup", "session tags: %s, %s", fTopSessionTag.Data(), fSessionTag.Data());

   // Get Session dir (sandbox)
   if ((fSessionDir = gEnv->GetValue("ProofServ.SessionDir", "-1")) == "-1") {
      Error("Setup", "Session dir missing");
      return -1;
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

   // Send packages off immediately to reduce latency
   fSocket->SetOption(kNoDelay, 1);

   // Check every two hours if client is still alive
   fSocket->SetOption(kKeepAlive, 1);

   // Install SigPipe handler to handle kKeepAlive failure
   gSystem->AddSignalHandler(new TXProofServSigPipeHandler(this));

   // Install Termination handler
   gSystem->AddSignalHandler(new TXProofServTerminationHandler(this));

   // Install seg violation handler
   gSystem->AddSignalHandler(new TXProofServSegViolationHandler(this));

   if (gProofDebugLevel > 0)
      Info("Setup", "successfully completed");

   // Done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Get list of workers to be used from now on.
/// The list must be provided by the caller.

TProofServ::EQueryAction TXProofServ::GetWorkers(TList *workers,
                                                 Int_t & /* prioritychange */,
                                                 Bool_t resume)
{
   TProofServ::EQueryAction rc = kQueryStop;

   // User config files, when enabled, override cluster-wide configuration
   if (gEnv->GetValue("ProofServ.UseUserCfg", 0) != 0) {
      Int_t pc = 1;
      return TProofServ::GetWorkers(workers, pc);
   }

   // seqnum of the query for which we call getworkers
   Bool_t dynamicStartup = gEnv->GetValue("Proof.DynamicStartup", kFALSE);
   TString seqnum = (dynamicStartup) ? "" : XPD_GW_Static;
   if (!fWaitingQueries->IsEmpty()) {
      if (resume) {
         seqnum += ((TProofQueryResult *)(fWaitingQueries->First()))->GetSeqNum();
      } else {
         seqnum += ((TProofQueryResult *)(fWaitingQueries->Last()))->GetSeqNum();
      }
   }
   // Send request to the coordinator
   TObjString *os = 0;
   if (dynamicStartup) {
      // We wait dynto seconds for the first worker to come; -1 means forever
      Int_t dynto = gEnv->GetValue("Proof.DynamicStartupTimeout", -1);
      Bool_t doto = (dynto > 0) ? kTRUE : kFALSE;
      while (!(os = ((TXSocket *)fSocket)->SendCoordinator(kGetWorkers, seqnum.Data()))) {
         if (doto > 0 && --dynto < 0) break;
         // Another second
         gSystem->Sleep(1000);
      }
   } else {
      os = ((TXSocket *)fSocket)->SendCoordinator(kGetWorkers, seqnum.Data());
   }

   // The reply contains some information about the master (image, workdir)
   // followed by the information about the workers; the tokens for each node
   // are separated by '&'
   if (os) {
      TString fl(os->GetName());
      if (fl.BeginsWith(XPD_GW_QueryEnqueued)) {
         SendAsynMessage("+++ Query cannot be processed now: enqueued");
         return kQueryEnqueued;
      }

      // Honour a max number of workers request (typically when running in valgrind)
      Int_t nwrks = -1;
      Bool_t pernode = kFALSE;
      if (gSystem->Getenv("PROOF_NWORKERS")) {
         TString s(gSystem->Getenv("PROOF_NWORKERS"));
         if (s.EndsWith("x")) {
            pernode = kTRUE;
            s.ReplaceAll("x", "");
         }
         if (s.IsDigit()) {
            nwrks = s.Atoi();
            if (!dynamicStartup && (nwrks > 0)) {
               // Notify, except in dynamic workers mode to avoid flooding
               TString msg;
               if (pernode) {
                  msg.Form("+++ Starting max %d workers per node following the setting of PROOF_NWORKERS", nwrks);
               } else {
                  msg.Form("+++ Starting max %d workers following the setting of PROOF_NWORKERS", nwrks);
               }
               SendAsynMessage(msg);
            } else {
               nwrks = -1;
            }
         } else {
            pernode = kFALSE;
         }
      }

      TString tok;
      Ssiz_t from = 0;
      TList *nodecnt = (pernode) ? new TList : 0 ;
      if (fl.Tokenize(tok, from, "&")) {
         if (!tok.IsNull()) {
            TProofNodeInfo *master = new TProofNodeInfo(tok);
            if (!master) {
               Error("GetWorkers", "no appropriate master line got from coordinator");
               return kQueryStop;
            } else {
               // Set image if not yet done and available
               if (fImage.IsNull() && strlen(master->GetImage()) > 0)
                  fImage = master->GetImage();
               SafeDelete(master);
            }
            // Now the workers
            while (fl.Tokenize(tok, from, "&")) {
               if (!tok.IsNull()) {
                  if (nwrks == -1 || nwrks > 0) {
                     // We have the minimal set of information to start
                     rc = kQueryOK;
                     if (pernode && nodecnt) {
                        TProofNodeInfo *ni = new TProofNodeInfo(tok);
                        TParameter<Int_t> *p = 0;
                        Int_t nw = 0;
                        if (!(p = (TParameter<Int_t> *) nodecnt->FindObject(ni->GetNodeName().Data()))) {
                           p = new TParameter<Int_t>(ni->GetNodeName().Data(), nw);
                           nodecnt->Add(p);
                        }
                        nw = p->GetVal();
                        if (gDebug > 0)
                           Info("GetWorkers","%p: name: %s (%s) val: %d (nwrks: %d)",
                                             p, p->GetName(), ni->GetNodeName().Data(),  nw, nwrks);
                        if (nw < nwrks) {
                           if (workers) workers->Add(ni);
                           nw++;
                           p->SetVal(nw);
                        } else {
                           // Two many workers on this machine already
                           SafeDelete(ni);
                        }
                     } else {
                        if (workers)
                           workers->Add(new TProofNodeInfo(tok));
                        // Count down
                        if (nwrks != -1) nwrks--;
                     }
                  } else {
                     // Release this worker (to cleanup the session list in the coordinator and get a fresh
                     // and correct list next call)
                     TProofNodeInfo *ni = new TProofNodeInfo(tok);
                     ReleaseWorker(ni->GetOrdinal().Data());
                  }
               }
            }
         }
      }
      // Cleanup
      if (nodecnt) {
         nodecnt->SetOwner(kTRUE);
         SafeDelete(nodecnt);
      }
   }

   // We are done
   return rc;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle error on the input socket

Bool_t TXProofServ::HandleError(const void *)
{
   // Try reconnection
   if (fSocket && !fSocket->IsValid()) {

      fSocket->Reconnect();
      if (fSocket && fSocket->IsValid()) {
         if (gDebug > 0)
            Info("HandleError",
               "%p: connection to local coordinator re-established", this);
         FlushLogFile();
         return kFALSE;
      }
   }
   Printf("TXProofServ::HandleError: %p: got called ...", this);

   // If master server, propagate interrupt to slaves
   // (shutdown interrupt send internally).
   if (IsMaster())
      fProof->Close("S");

   // Avoid communicating back anything to the coordinator (it is gone)
   if (fSocket) ((TXSocket *)fSocket)->SetSessionID(-1);

   Terminate(0);

   Printf("TXProofServ::HandleError: %p: DONE ... ", this);

   // We are done
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle asynchronous input on the input socket

Bool_t TXProofServ::HandleInput(const void *in)
{
   if (gDebug > 2)
      Printf("TXProofServ::HandleInput %p, in: %p", this, in);

   XHandleIn_t *hin = (XHandleIn_t *) in;
   Int_t acod = (hin) ? hin->fInt1 : kXPD_msg;

   // Act accordingly
   if (acod == kXPD_ping || acod == kXPD_interrupt) {
      // Interrupt or Ping
      HandleUrgentData();

   } else if (acod == kXPD_flush) {
      // Flush stdout, so that we can access the full log file
      Info("HandleInput","kXPD_flush: flushing log file (stdout)");
      fflush(stdout);

   } else if (acod == kXPD_urgent) {
      // Get type
      Int_t type = hin->fInt2;
      switch (type) {
      case TXSocket::kStopProcess:
         {
            // Abort or Stop ?
            Bool_t abort = (hin->fInt3 != 0) ? kTRUE : kFALSE;
            // Timeout
            Int_t timeout = hin->fInt4;
            // Act now
            if (fProof)
               fProof->StopProcess(abort, timeout);
            else
               if (fPlayer)
                  fPlayer->StopProcess(abort, timeout);
         }
         break;
      default:
         Info("HandleInput","kXPD_urgent: unknown type: %d", type);
      }

   } else if (acod == kXPD_inflate) {

      // Obsolete type
      Warning("HandleInput", "kXPD_inflate: obsolete message type");

   } else if (acod == kXPD_priority) {

      // The factor is the priority to be propagated
      fGroupPriority = hin->fInt2;
      if (fProof)
         fProof->BroadcastGroupPriority(fGroup, fGroupPriority);
      // Notify
      Info("HandleInput", "kXPD_priority: group %s priority set to %f",
           fGroup.Data(), (Float_t) fGroupPriority / 100.);

   } else if (acod == kXPD_clusterinfo) {

      // Information about the cluster status
      fTotSessions     = hin->fInt2;
      fActSessions     = hin->fInt3;
      fEffSessions     = (hin->fInt4)/1000.;
      // Notify
      Info("HandleInput", "kXPD_clusterinfo: tot: %d, act: %d, eff: %f",
           fTotSessions, fActSessions, fEffSessions);

   } else {
      // Standard socket input
      HandleSocketInput();
      // This request has been completed: remove the client ID from the pipe
      ((TXSocket *)fSocket)->RemoveClientID();
   }

   // We are done
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Disable read timeout on the underlying socket

void TXProofServ::DisableTimeout()
{
   if (fSocket)
     ((TXSocket *)fSocket)->DisableTimeout();
}

////////////////////////////////////////////////////////////////////////////////
/// Enable read timeout on the underlying socket

void TXProofServ::EnableTimeout()
{
   if (fSocket)
     ((TXSocket *)fSocket)->EnableTimeout();
}

////////////////////////////////////////////////////////////////////////////////
/// Terminate the proof server.

void TXProofServ::Terminate(Int_t status)
{
   if (fTerminated)
      // Avoid doubling the exit operations
      exit(1);
   fTerminated = kTRUE;

   // Notify
   Info("Terminate", "starting session termination operations ...");
   if (fgLogToSysLog > 0) {
      TString s;
      s.Form("%s -1 %.3f %.3f", fgSysLogEntity.Data(), fRealTime, fCpuTime);
      gSystem->Syslog(kLogNotice, s.Data());
   }

   // Notify the memory footprint
   ProcInfo_t pi;
   if (!gSystem->GetProcInfo(&pi)){
      Info("Terminate", "process memory footprint: %ld/%ld kB virtual, %ld/%ld kB resident ",
                        pi.fMemVirtual, fgVirtMemMax, pi.fMemResident, fgResMemMax);
   }

   // Deactivate current monitor, if any
   if (fProof)
      fProof->SetMonitor(0, kFALSE);

   // Cleanup session directory
   if (status == 0) {
      // make sure we remain in a "connected" directory
      gSystem->ChangeDirectory("/");
      // needed in case fSessionDir is on NFS ?!
      gSystem->MakeDirectory(fSessionDir+"/.delete");
      gSystem->Exec(Form("%s %s", kRM, fSessionDir.Data()));
   }

   // Cleanup queries directory if empty
   if (IsMaster()) {
      if (!(fQMgr && fQMgr->Queries() && fQMgr->Queries()->GetSize())) {
         // make sure we remain in a "connected" directory
         gSystem->ChangeDirectory("/");
         // needed in case fQueryDir is on NFS ?!
         gSystem->MakeDirectory(fQueryDir+"/.delete");
         gSystem->Exec(Form("%s %s", kRM, fQueryDir.Data()));
         // Remove lock file
         if (fQueryLock)
            gSystem->Unlink(fQueryLock->GetName());
       }

      // Unlock the query dir owned by this session
      if (fQueryLock)
         fQueryLock->Unlock();
   } else {
      // Try to stop processing if any
      Bool_t abort = (status == 0) ? kFALSE : kTRUE;
      if (!fIdle && fPlayer)
         fPlayer->StopProcess(abort,1);
      gSystem->Sleep(2000);
   }

   // Cleanup data directory if empty
   if (!fDataDir.IsNull() && !gSystem->AccessPathName(fDataDir, kWritePermission)) {
     if (UnlinkDataDir(fDataDir))
        Info("Terminate", "data directory '%s' has been removed", fDataDir.Data());
   }

   // Remove input and signal handlers to avoid spurious "signals"
   // for closing activities executed upon exit()
   gSystem->RemoveFileHandler(fInputHandler);

   // Stop processing events (set a flag to exit the event loop)
   gSystem->ExitLoop();

   // We post the pipe once to wake up the main thread which is waiting for
   // activity on this socket; this fake activity will make it return and
   // eventually exit the loop.
   TXSocket::fgPipe.Post((TXSocket *)fSocket);

   // Notify
   Printf("Terminate: termination operations ended: quitting!");
}

////////////////////////////////////////////////////////////////////////////////
/// Try locking query area of session tagged sessiontag.
/// The id of the locking file is returned in fid and must be
/// unlocked via UnlockQueryFile(fid).

Int_t TXProofServ::LockSession(const char *sessiontag, TProofLockPath **lck)
{
   // We do not need to lock our own session
   if (strstr(sessiontag, fTopSessionTag))
      return 0;

   if (!lck) {
      Info("LockSession","locker space undefined");
      return -1;
   }
   *lck = 0;

   // Check the format
   TString stag = sessiontag;
   TRegexp re("session-.*-.*-.*");
   Int_t i1 = stag.Index(re);
   if (i1 == kNPOS) {
      Info("LockSession","bad format: %s", sessiontag);
      return -1;
   }
   stag.ReplaceAll("session-","");

   // Drop query number, if any
   Int_t i2 = stag.Index(":q");
   if (i2 != kNPOS)
      stag.Remove(i2);

   // Make sure that parent process does not exist anylonger
   TString parlog = fSessionDir;
   parlog = parlog.Remove(parlog.Index("master-")+strlen("master-"));
   parlog += stag;
   if (!gSystem->AccessPathName(parlog)) {
      Info("LockSession","parent still running: do nothing");
      return -1;
   }

   // Lock the query lock file
   TString qlock = fQueryLock->GetName();
   qlock.ReplaceAll(fTopSessionTag, stag);

   if (!gSystem->AccessPathName(qlock)) {
      *lck = new TProofLockPath(qlock);
      if (((*lck)->Lock()) < 0) {
         Info("LockSession","problems locking query lock file");
         SafeDelete(*lck);
         return -1;
      }
   }

   // We are done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Send message to intermediate coordinator to release worker of last ordinal
/// ord.

void TXProofServ::ReleaseWorker(const char *ord)
{
   if (gDebug > 2) Info("ReleaseWorker","releasing: %s", ord);

   ((TXSocket *)fSocket)->SendCoordinator(kReleaseWorker, ord);
}
