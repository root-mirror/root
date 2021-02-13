// @(#)root/proofx:$Id$
// Author: Gerardo Ganis  12/12/2005

/*************************************************************************
 * Copyright (C) 1995-2005, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TXSlave
\ingroup proofx

This is the version of TSlave for workers servers based on XProofD.
See TSlave and TXSocket for details.

*/

#include "TXSlave.h"
#include "TProof.h"
#include "TProofServ.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TUrl.h"
#include "TMessage.h"
#include "TMonitor.h"
#include "TError.h"
#include "TSysEvtHandler.h"
#include "TVirtualMutex.h"
#include "TXSocket.h"
#include "TXSocketHandler.h"
#include "Varargs.h"
#include "XProofProtocol.h"

#include <mutex>

ClassImp(TXSlave);

//______________________________________________________________________________

//---- Hook to the constructor -------------------------------------------------
//---- This is needed to avoid using the plugin manager which may create -------
//---- problems in multi-threaded environments. --------------------------------
TSlave *GetTXSlave(const char *url, const char *ord, Int_t perf,
                     const char *image, TProof *proof, Int_t stype,
                     const char *workdir, const char *msd, Int_t nwk)
{
   return ((TSlave *)(new TXSlave(url, ord, perf, image,
                                    proof, stype, workdir, msd, nwk)));
}

class XSlaveInit {
 public:
   XSlaveInit() {
      TSlave::SetTXSlaveHook(&GetTXSlave);
}};
static XSlaveInit xslave_init;

//______________________________________________________________________________

//---- error handling ----------------------------------------------------------
//---- Needed to avoid blocking on the CINT mutex in printouts -----------------

////////////////////////////////////////////////////////////////////////////////
/// Interface to ErrorHandler (protected).

void TXSlave::DoError(int level, const char *location, const char *fmt, va_list va) const
{
   ::ErrorHandler(level, Form("TXSlave::%s", location), fmt, va);
}

//
// Specific Interrupt signal handler
//
class TXSlaveInterruptHandler : public TSignalHandler {
private:
   TXSocket *fSocket;
public:
   TXSlaveInterruptHandler(TXSocket *s = 0)
      : TSignalHandler(kSigInterrupt, kFALSE), fSocket(s) { }
   Bool_t Notify();
};

////////////////////////////////////////////////////////////////////////////////
/// TXSlave interrupt handler.

Bool_t TXSlaveInterruptHandler::Notify()
{
   Info("Notify","Processing interrupt signal ...");

   // Handle also interrupt condition on socket(s)
   if (fSocket)
      fSocket->SetInterrupt();

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a PROOF slave object. Called via the TProof ctor.

TXSlave::TXSlave(const char *url, const char *ord, Int_t perf,
               const char *image, TProof *proof, Int_t stype,
               const char *workdir, const char *msd, Int_t nwk) : TSlave()
{
   fImage = image;
   fProofWorkDir = workdir;
   fWorkDir = workdir;
   fOrdinal = ord;
   fPerfIdx = perf;
   fProof = proof;
   fSlaveType = (ESlaveType)stype;
   fMsd = msd;
   fNWrks = nwk;
   fIntHandler = 0;
   fValid = kFALSE;

   // Instance of the socket input handler to monitor all the XPD sockets
   TXSocketHandler *sh = TXSocketHandler::GetSocketHandler();
   gSystem->AddFileHandler(sh);

   TXSocket::SetLocation((fProof->IsMaster()) ? "master" : "client");

   Init(url, stype);
}

////////////////////////////////////////////////////////////////////////////////
/// Init a PROOF slave object. Called via the TXSlave ctor.
/// The Init method is technology specific and is overwritten by derived
/// classes.

void TXSlave::Init(const char *host, Int_t stype)
{
   // Url string with host, port information; 'host' may contain 'user' information
   // in the usual form 'user@host'

   // Auxilliary url
   TUrl url(host);
   url.SetProtocol(fProof->fUrl.GetProtocol());
   // Check port
   if (url.GetPort() == TUrl("a").GetPort()) {
      // We use 'rootd' service as default.
      Int_t port = gSystem->GetServiceByName("proofd");
      if (port < 0) {
         if (gDebug > 0)
            Info("Init","service 'proofd' not found by GetServiceByName"
                        ": using default IANA assigned tcp port 1093");
         port = 1093;
      } else {
         if (gDebug > 1)
            Info("Init","port from GetServiceByName: %d", port);
      }
      url.SetPort(port);
   }

   // Fill members
   fName = url.GetHostFQDN();
   fPort = url.GetPort(); // We get the right default if the port is not specified
   // Group specification , if any, uses the password field, i.e. user[:group]
   fGroup = url.GetPasswd();

   // The field 'psid' is interpreted as session ID when we are attaching
   // to an existing session (ID passed in the options field of the url) or
   // to our PROOF protocl version when we are creating a new session
   TString opts(url.GetOptions());
   Bool_t attach = (opts.Length() > 0 && opts.IsDigit()) ? kTRUE : kFALSE;
   Int_t psid = (attach) ? opts.Atoi() : kPROOF_Protocol;

   // Add information about our status (Client or Master)
   TString iam;
   Char_t mode = 's';
   TString alias = fProof->GetTitle();
   if (fProof->IsMaster() && stype == kSlave) {
      iam = "Master";
      mode = 's';
      // Send session tag of the closest master to the slaves
      alias.Form("session-%s|ord:%s", fProof->GetName(), fOrdinal.Data());
   } else if (fProof->IsMaster() && stype == kMaster) {
      iam = "Master";
      mode = 'm';
      // Send session tag of the closest master to the slaves
      if (fNWrks > 1) {
         alias.Form("session-%s|ord:%s|plite:%d", fProof->GetName(), fOrdinal.Data(), fNWrks);
         mode = 'L';
      } else {
         alias.Form("session-%s|ord:%s", fProof->GetName(), fOrdinal.Data());
      }
   } else if (!fProof->IsMaster() && stype == kMaster) {
      iam = "Local Client";
      mode = (attach) ? 'A' : 'M';
   } else {
      Error("Init","Impossible PROOF <-> SlaveType Configuration Requested");
      R__ASSERT(0);
   }

   // Add conf file, if required
   if (fProof->fConfFile.Length() > 0 && fNWrks <= 1)
      alias += Form("|cf:%s",fProof->fConfFile.Data());

   // Send over env variables (may not be supported remotely)
   TString envlist;
   if (!fProof->GetManager() ||
        fProof->GetManager()->GetRemoteProtocol() > 1001) {
         // Check if the user forced locally a given authentication protocol:
         // we need to do the same remotely to get the right credentials
         if (gSystem->Getenv("XrdSecPROTOCOL")) {
            TProof::DelEnvVar("XrdSecPROTOCOL");
            TProof::AddEnvVar("XrdSecPROTOCOL", gSystem->Getenv("XrdSecPROTOCOL"));
         }
         const TList *envs = TProof::GetEnvVars();
         if (envs != 0 ) {
            TIter next(envs);
            for (TObject *o = next(); o != 0; o = next()) {
               TNamed *env = dynamic_cast<TNamed*>(o);
               if (env != 0) {
                  if (!envlist.IsNull())
                     envlist += ",";
                  envlist += Form("%s=%s", env->GetName(), env->GetTitle());
               }
            }
         }
   } else {
      if (fProof->GetManager() && TProof::GetEnvVars())
         Info("Init", "** NOT ** sending user envs - RemoteProtocol : %d",
                      fProof->GetManager()->GetRemoteProtocol());
   }

   // Add to the buffer
   if (!envlist.IsNull())
      alias += Form("|envs:%s", envlist.Data());

   // Open connection to a remote XrdPROOF slave server.
   // Login and authentication are dealt with at this level, if required.
   if (!(fSocket = new TXSocket(url.GetUrl(kTRUE), mode, psid,
                                -1, alias, fProof->GetLogLevel(), this))) {
      ParseBuffer(); // For the log path
      Error("Init", "while opening the connection to %s - exit", url.GetUrl(kTRUE));
      return;
   }

   // The socket may not be valid
   if (!(fSocket->IsValid())) {
      // Notify only if verbosity is on: most likely the failure has already been notified
      PDB(kGlobal,1)
         Error("Init", "some severe error occurred while opening "
                       "the connection at %s - exit", url.GetUrl(kTRUE));
      ParseBuffer(); // For the log path
      // Fill some useful info
      fUser = ((TXSocket *)fSocket)->fUser;
      PDB(kGlobal,3) Info("Init","%s: fUser is .... %s", iam.Data(), fUser.Data());
      SafeDelete(fSocket);
      return;
   }

   // Set the ordinal in the title for debugging
   fSocket->SetTitle(fOrdinal);

   // Check if the remote server supports user envs setting
   if (!fProof->GetManager() && !envlist.IsNull() &&
      ((TXSocket *)fSocket)->GetXrdProofdVersion() <= 1001) {
      Info("Init","user envs setting sent but unsupported remotely - RemoteProtocol : %d",
                     ((TXSocket *)fSocket)->GetXrdProofdVersion());
   }

   // Set the reference to TProof
   ((TXSocket *)fSocket)->fReference = fProof;

   // Protocol run by remote PROOF server
   fProtocol = fSocket->GetRemoteProtocol();

   // Set server type
   fProof->fServType = TProofMgr::kXProofd;

   // Set remote session ID
   fProof->fSessionID = ((TXSocket *)fSocket)->GetSessionID();

   // Extract the log file path and, if any, set URL entry point for the default data pool
   ParseBuffer();

   // Remove socket from global TROOT socket list. Only the TProof object,
   // representing all slave sockets, will be added to this list. This will
   // ensure the correct termination of all proof servers in case the
   // root session terminates.
   {
      R__LOCKGUARD(gROOTMutex);
      gROOT->GetListOfSockets()->Remove(fSocket);
   }

   // Fill some useful info
   fUser = ((TXSocket *)fSocket)->fUser;
   PDB(kGlobal,3) {
      Info("Init","%s: fUser is .... %s", iam.Data(), fUser.Data());
   }

   // Set valid
   fValid = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Parse fBuffer after a connection attempt

void TXSlave::ParseBuffer()
{
   // Set URL entry point for the default data pool
   TString buffer(((TXSocket *)fSocket)->fBuffer);
   if (!buffer.IsNull()) {
      Ssiz_t ilog = buffer.Index("|log:");
      if (ilog != 0) {
         // Extract the pool URL (on master)
         TString dpu = (ilog != kNPOS) ? buffer(0, ilog) : buffer;
         if (dpu.Length() > 0) fProof->SetDataPoolUrl(dpu);
      }
      if (ilog != kNPOS) {
         // The rest, if any, if the log file path from which we extract the working dir
         buffer.Remove(0, ilog + sizeof("|log:") - 1);
         fWorkDir = buffer;
         if ((ilog = fWorkDir.Last('.')) !=  kNPOS) fWorkDir.Remove(ilog);
         if (gDebug > 2)
            Info("ParseBuffer", "workdir is: %s", fWorkDir.Data());
      } else if (fProtocol > 31) {
         Warning("ParseBuffer", "expected log path not found in received startup buffer!");
      }
   }
   // Done
   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Init a PROOF slave object. Called via the TXSlave ctor.
/// The Init method is technology specific and is overwritten by derived
/// classes.

Int_t TXSlave::SetupServ(Int_t, const char *)
{
   // Get back startup message of proofserv (we are now talking with
   // the real proofserver and not anymore with the proofd front-end)
   Int_t what;
   char buf[512];
   if (fSocket->Recv(buf, sizeof(buf), what) <= 0) {
      Error("SetupServ", "failed to receive slave startup message");
      Close("S");
      SafeDelete(fSocket);
      fValid = kFALSE;
      return -1;
   }

   if (what == kMESS_NOTOK) {
      SafeDelete(fSocket);
      fValid = kFALSE;
      return -1;
   }

   // protocols less than 4 are incompatible
   if (fProtocol < 4) {
      Error("SetupServ", "incompatible PROOF versions (remote version "
                         "must be >= 4, is %d)", fProtocol);
      SafeDelete(fSocket);
      fValid = kFALSE;
      return -1;
   }

   fProof->fProtocol   = fProtocol;   // protocol of last slave on master

   // set some socket options
   fSocket->SetOption(kNoDelay, 1);

   // We are done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Destroy slave.

TXSlave::~TXSlave()
{
   Close();
}

////////////////////////////////////////////////////////////////////////////////
/// Close slave socket.

void TXSlave::Close(Option_t *opt)
{
   if (fSocket)
      // Closing socket ...
      fSocket->Close(opt);

   SafeDelete(fInput);
   SafeDelete(fSocket);
}

////////////////////////////////////////////////////////////////////////////////
/// Ping the remote master or slave servers.
/// Returns 0 if ok, -1 if it did not ping or in case of error

Int_t TXSlave::Ping()
{
   if (!IsValid()) return -1;

   return (((TXSocket *)fSocket)->Ping(GetOrdinal()) ? 0 : -1);
}

////////////////////////////////////////////////////////////////////////////////
/// Touch the client admin file to proof we are alive.

void TXSlave::Touch()
{
   if (!IsValid()) return;

   ((TXSocket *)fSocket)->RemoteTouch();
   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Send interrupt to master or slave servers.
/// Returns 0 if ok, -1 in case of error

void TXSlave::Interrupt(Int_t type)
{
   if (!IsValid()) return;

   if (type == TProof::kLocalInterrupt) {

      // Deactivate and flush the local socket (we are not - yet - closing
      // the session, so we do less things that in case of an error ...)
      if (fProof) {

         // Attach to the monitor instance, if any
         TMonitor *mon = fProof->fCurrentMonitor;
         TList *al = mon ? mon->GetListOfActives() : nullptr;
         if (mon && fSocket && al && al->FindObject(fSocket)) {
            // Synchronous collection in TProof
            if (gDebug > 2)
               Info("Interrupt", "%p: deactivating from monitor %p", this, mon);
            mon->DeActivate(fSocket);
         }
         delete al;
      } else {
         Warning("Interrupt", "%p: reference to PROOF missing", this);
      }

      // Post semaphore to wake up anybody waiting
      if (fSocket) ((TXSocket *)fSocket)->PostSemAll();

      return;
   }

   if (fSocket) ((TXSocket *)fSocket)->SendInterrupt(type);
   Info("Interrupt","Interrupt of type %d sent", type);
}

////////////////////////////////////////////////////////////////////////////////
/// Sent stop/abort request to PROOF server. It will be
/// processed asynchronously by a separate thread.

void TXSlave::StopProcess(Bool_t abort, Int_t timeout)
{
   if (!IsValid()) return;

   ((TXSocket *)fSocket)->SendUrgent(TXSocket::kStopProcess, (Int_t)abort, timeout);
   if (gDebug > 0)
      Info("StopProcess", "Request of type %d sent over", abort);
}

////////////////////////////////////////////////////////////////////////////////
/// Find out the remote proofd protocol version.
/// Returns -1 in case of error.

Int_t TXSlave::GetProofdProtocol(TSocket *s)
{
   Int_t rproto = -1;

   UInt_t cproto = 0;
   Int_t len = sizeof(cproto);
   memcpy((char *)&cproto,
      Form(" %d", TSocket::GetClientProtocol()),len);
   Int_t ns = s->SendRaw(&cproto, len);
   if (ns != len) {
      ::Error("TXSlave::GetProofdProtocol",
              "sending %d bytes to proofd server [%s:%d]",
              len, (s->GetInetAddress()).GetHostName(), s->GetPort());
      return -1;
   }

   // Get the remote protocol
   Int_t ibuf[2] = {0};
   len = sizeof(ibuf);
   Int_t nr = s->RecvRaw(ibuf, len);
   if (nr != len) {
      ::Error("TXSlave::GetProofdProtocol",
              "reading %d bytes from proofd server [%s:%d]",
              len, (s->GetInetAddress()).GetHostName(), s->GetPort());
      return -1;
   }
   Int_t kind = net2host(ibuf[0]);
   if (kind == kROOTD_PROTOCOL) {
      rproto = net2host(ibuf[1]);
   } else {
      kind = net2host(ibuf[1]);
      if (kind == kROOTD_PROTOCOL) {
         len = sizeof(rproto);
         nr = s->RecvRaw(&rproto, len);
         if (nr != len) {
            ::Error("TXSlave::GetProofdProtocol",
                    "reading %d bytes from proofd server [%s:%d]",
                    len, (s->GetInetAddress()).GetHostName(), s->GetPort());
            return -1;
         }
         rproto = net2host(rproto);
      }
   }
   if (gDebug > 2)
      ::Info("TXSlave::GetProofdProtocol",
             "remote proofd: buf1: %d, buf2: %d rproto: %d",
             net2host(ibuf[0]),net2host(ibuf[1]),rproto);

   // We are done
   return rproto;
}

////////////////////////////////////////////////////////////////////////////////
/// Send message to intermediate coordinator.
/// If any output is due, this is returned as a generic message

TObjString *TXSlave::SendCoordinator(Int_t kind, const char *msg, Int_t int2)
{
   return ((TXSocket *)fSocket)->SendCoordinator(kind, msg, int2);
}

////////////////////////////////////////////////////////////////////////////////
/// Set an alias for this session. If reconnection is supported, the alias
/// will be communicated to the remote coordinator so that it can be recovered
/// when reconnecting

void TXSlave::SetAlias(const char *alias)
{
   // Nothing to do if not in contact with coordinator
   if (!IsValid()) return;

   ((TXSocket *)fSocket)->SendCoordinator(kSessionAlias, alias);

   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Communicate to the coordinator the priprity of the group to which the
/// user belongs
/// Return 0 on success

Int_t TXSlave::SendGroupPriority(const char *grp, Int_t priority)
{
   // Nothing to do if not in contact with coordinator
   if (!IsValid()) return -1;

   ((TXSocket *)fSocket)->SendCoordinator(kGroupProperties, grp, priority);

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle error on the input socket

Bool_t TXSlave::HandleError(const void *in)
{
   XHandleErr_t *herr = in ? (XHandleErr_t *)in : 0;

   // Try reconnection
   if (fSocket && herr && (herr->fOpt == 1)) {

      ((TXSocket *)fSocket)->Reconnect();
      if (fSocket && fSocket->IsValid()) {
         if (gDebug > 0) {
            if (!strcmp(GetOrdinal(), "0")) {
               Printf("Proof: connection to master at %s:%d re-established",
                     GetName(), GetPort());
            } else {
               Printf("Proof: connection to node '%s' at %s:%d re-established",
                     GetOrdinal(), GetName(), GetPort());
            }
         }
         return kFALSE;
      }
   }

   // This seems a real error: notify the interested parties
   Info("HandleError", "%p:%s:%s got called ... fProof: %p, fSocket: %p (valid: %d)",
                       this, fName.Data(), fOrdinal.Data(), fProof, fSocket,
                       (fSocket ? (Int_t)fSocket->IsValid() : -1));

   // Remove interrupt handler (avoid affecting other clients of the underlying physical
   // connection)
   SetInterruptHandler(kFALSE);

   if (fProof) {

      // Remove PROOF signal handler
      if (fProof->fIntHandler)
         fProof->fIntHandler->Remove();

      Info("HandleError", "%p: proof: %p", this, fProof);

      if (fSocket) {
         // This is need to skip contacting the remote server upon close
         ((TXSocket *)fSocket)->SetSessionID(-1);
         // This is need to interrupt possible pickup waiting status
         ((TXSocket *)fSocket)->SetInterrupt();
         // Synchronous collection in TProof: post fatal message; this will
         // mark the worker as bad and update the internal lists accordingly
         ((TXSocket *)fSocket)->PostMsg(kPROOF_FATAL);
      }

      // On masters we notify clients of the problem occured
      if (fProof->IsMaster()) {
         TString msg(Form("Worker '%s-%s' has been removed from the active list",
                          fName.Data(), fOrdinal.Data()));
         TMessage m(kPROOF_MESSAGE);
         m << msg;
         if (gProofServ)
            gProofServ->GetSocket()->Send(m);
         else
            Warning("HandleError", "%p: global reference to TProofServ missing", this);
      }
   } else {
      Warning("HandleError", "%p: reference to PROOF missing", this);
   }

   Printf("TXSlave::HandleError: %p: DONE ... ", this);

   // We are done
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle asynchronous input on the socket

Bool_t TXSlave::HandleInput(const void *)
{
   if (fProof) {

      // Attach to the monitor instance, if any
      TMonitor *mon = fProof->fCurrentMonitor;

      if (gDebug > 2)
         Info("HandleInput", "%p: %s: proof: %p, mon: %p",
                             this, GetOrdinal(), fProof, mon);

      if (mon && mon->IsActive(fSocket)) {
         // Synchronous collection in TProof
         if (gDebug > 2)
            Info("HandleInput","%p: %s: posting monitor %p", this, GetOrdinal(), mon);
         mon->SetReady(fSocket);
      } else {
         // Asynchronous collection in TProof
         if (gDebug > 2) {
            if (mon) {
               Info("HandleInput", "%p: %s: not active in current monitor"
                                   " - calling TProof::CollectInputFrom",
                                   this, GetOrdinal());
            } else {
               Info("HandleInput", "%p: %s: calling TProof::CollectInputFrom",
                                   this, GetOrdinal());
            }
         }
         if (fProof->CollectInputFrom(fSocket) < 0)
            // Something wrong on the line: flush it
            FlushSocket();
      }
   } else {
      Warning("HandleInput", "%p: %s: reference to PROOF missing", this, GetOrdinal());
      return kFALSE;
   }

   // We are done
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set/Unset the interrupt handler

void TXSlave::SetInterruptHandler(Bool_t on)
{
   if (gDebug > 1)
      Info("SetInterruptHandler", "enter: %d", on);

   if (on) {
      if (!fIntHandler)
         fIntHandler = new TXSlaveInterruptHandler((TXSocket *)fSocket);
      fIntHandler->Add();
   } else {
      if (fIntHandler)
         fIntHandler->Remove();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Clean any input on the socket

void TXSlave::FlushSocket()
{
   if (gDebug > 1)
      Info("FlushSocket", "enter: %p", fSocket);

   if (fSocket)
      TXSocket::fgPipe.Flush(fSocket);
}
