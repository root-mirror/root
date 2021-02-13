// @(#)root/proof:$Id$
// Author: Fons Rademakers   14/02/97

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TSlave
\ingroup proofkernel

Class describing a PROOF worker server. It contains information like the
workers host name, ordinal number, performance index, socket, etc.
Objects of this class can only be created via TProof member functions.

*/

#include <stdlib.h>

#include "RConfigure.h"
#include "TApplication.h"
#include "TSlave.h"
#include "TSlaveLite.h"
#include "TProof.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TUrl.h"
#include "TMessage.h"
#include "TError.h"
#include "TVirtualMutex.h"
#include "TSocket.h"
#include "TObjString.h"

ClassImp(TSlave);

// Hook for the TXSlave constructor
TSlave_t TSlave::fgTXSlaveHook = 0;

////////////////////////////////////////////////////////////////////////////////
/// Create a PROOF slave object. Called via the TProof ctor.

TSlave::TSlave(const char *url, const char *ord, Int_t perf,
               const char *image, TProof *proof, Int_t stype,
               const char *workdir, const char *msd, Int_t)
  : fImage(image), fProofWorkDir(workdir),
    fWorkDir(workdir), fPort(-1),
    fOrdinal(ord), fPerfIdx(perf),
    fProtocol(0), fSocket(0), fProof(proof),
    fInput(0), fBytesRead(0), fRealTime(0),
    fCpuTime(0), fSlaveType((ESlaveType)stype), fStatus(TSlave::kInvalid),
    fParallel(0), fMsd(msd)
{
   fName = TUrl(url).GetHostFQDN();
   fPort = TUrl(url).GetPort();

   Init(url, -1, stype);
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor used by derived classes

TSlave::TSlave()
{
   fPort      = -1;
   fOrdinal   = "-1";
   fPerfIdx   = -1;
   fProof     = 0;
   fSlaveType = kMaster;
   fProtocol  = 0;
   fSocket    = 0;
   fInput     = 0;
   fBytesRead = 0;
   fRealTime  = 0;
   fCpuTime   = 0;
   fStatus    = kInvalid;
   fParallel  = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Init a PROOF slave object. Called via the TSlave ctor.
/// The Init method is technology specific and is overwritten by derived
/// classes.

void TSlave::Init(const char *host, Int_t port, Int_t stype)
{
   // The url contains information about the server type: make sure
   // it is 'proofd' or alike
   TString proto = fProof->fUrl.GetProtocol();
   proto.Insert(5, 'd');

   TUrl hurl(host);
   hurl.SetProtocol(proto);
   if (port > 0)
      hurl.SetPort(port);

   // Add information about our status (Client or Master)
   TString iam;
   if (fProof->IsMaster() && stype == kSlave) {
      iam = "Master";
      hurl.SetOptions("SM");
   } else if (fProof->IsMaster() && stype == kMaster) {
      iam = "Master";
      hurl.SetOptions("MM");
   } else if (!fProof->IsMaster() && stype == kMaster) {
      iam = "Local Client";
      hurl.SetOptions("MC");
   } else {
      Error("Init","Impossible PROOF <-> SlaveType Configuration Requested");
      R__ASSERT(0);
   }

   // Open authenticated connection to remote PROOF slave server.
   // If a connection was already open (fSocket != 0), re-use it
   // to perform authentication (optimization needed to avoid a double
   // opening in case this is called by TXSlave).
   Int_t wsize = 65536;
   fSocket = TSocket::CreateAuthSocket(hurl.GetUrl(), 0, wsize, fSocket);

   if (!fSocket || !fSocket->IsAuthenticated()) {
      SafeDelete(fSocket);
      return;
   }

   // Remove socket from global TROOT socket list. Only the TProof object,
   // representing all slave sockets, will be added to this list. This will
   // ensure the correct termination of all proof servers in case the
   // root session terminates.
   {
      R__LOCKGUARD(gROOTMutex);
      gROOT->GetListOfSockets()->Remove(fSocket);
   }

   // Fill some useful info
   fUser              = fSocket->GetSecContext()->GetUser();
   PDB(kGlobal,3) {
      Info("Init","%s: fUser is .... %s", iam.Data(), fUser.Data());
   }

   if (fSocket->GetRemoteProtocol() >= 14 ) {
      TMessage m(kPROOF_SETENV);

      const TList *envs = TProof::GetEnvVars();
      if (envs != 0 ) {
         TIter next(envs);
         for (TObject *o = next(); o != 0; o = next()) {
            TNamed *env = dynamic_cast<TNamed*>(o);
            if (env != 0) {
               TString def = Form("%s=%s", env->GetName(), env->GetTitle());
               const char *p = def.Data();
               m << p;
            }
         }
      }
      fSocket->Send(m);
   } else {
      Info("Init","** NOT ** Sending kPROOF_SETENV RemoteProtocol : %d",
         fSocket->GetRemoteProtocol());
   }

   char buf[512];
   fSocket->Recv(buf, sizeof(buf));
   if (strcmp(buf, "Okay")) {
      Printf("%s", buf);
      SafeDelete(fSocket);
      return;
   }

}

////////////////////////////////////////////////////////////////////////////////
/// Init a PROOF slave object. Called via the TSlave ctor.
/// The Init method is technology specific and is overwritten by derived
/// classes.

Int_t TSlave::SetupServ(Int_t stype, const char *conffile)
{
   // get back startup message of proofserv (we are now talking with
   // the real proofserver and not anymore with the proofd front-end)
   Int_t what;
   char buf[512];
   if (fSocket->Recv(buf, sizeof(buf), what) <= 0) {
      Error("SetupServ", "failed to receive slave startup message");
      SafeDelete(fSocket);
      return -1;
   }

   if (what == kMESS_NOTOK) {
      SafeDelete(fSocket);
      return -1;
   }

   // exchange protocol level between client and master and between
   // master and slave
   if (fSocket->Send(kPROOF_Protocol, kROOTD_PROTOCOL) != 2*sizeof(Int_t)) {
      Error("SetupServ", "failed to send local PROOF protocol");
      SafeDelete(fSocket);
      return -1;
   }

   if (fSocket->Recv(fProtocol, what) != 2*sizeof(Int_t)) {
      Error("SetupServ", "failed to receive remote PROOF protocol");
      SafeDelete(fSocket);
      return -1;
   }

   // protocols less than 4 are incompatible
   if (fProtocol < 4) {
      Error("SetupServ", "incompatible PROOF versions (remote version"
                      " must be >= 4, is %d)", fProtocol);
      SafeDelete(fSocket);
      return -1;
   }

   fProof->fProtocol   = fProtocol;   // protocol of last slave on master

   if (fProtocol < 5) {
      //
      // Setup authentication related stuff for ald versions
      Bool_t isMaster = (stype == kMaster);
      TString wconf = isMaster ? TString(conffile) : fProofWorkDir;
      if (OldAuthSetup(isMaster, wconf) != 0) {
         Error("SetupServ", "OldAuthSetup: failed to setup authentication");
         SafeDelete(fSocket);
         return -1;
      }
   } else {
      //
      // Send ordinal (and config) info to slave (or master)
      TMessage mess;
      if (stype == kMaster)
         mess << fUser << fOrdinal << TString(conffile);
      else
         mess << fUser << fOrdinal << fProofWorkDir;

      if (fSocket->Send(mess) < 0) {
         Error("SetupServ", "failed to send ordinal and config info");
         SafeDelete(fSocket);
         return -1;
      }
   }

   // set some socket options
   fSocket->SetOption(kNoDelay, 1);

   // Set active state
   fStatus = kActive;

   // We are done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Init a PROOF slave object using the connection opened via s. Used to
/// avoid double opening when an attempt via TXSlave found a remote proofd.

void TSlave::Init(TSocket *s, Int_t stype)
{
   fSocket = s;
   TSlave::Init(s->GetInetAddress().GetHostName(), s->GetPort(), stype);
}

////////////////////////////////////////////////////////////////////////////////
/// Destroy slave.

TSlave::~TSlave()
{
   Close();
}

////////////////////////////////////////////////////////////////////////////////
/// Close slave socket.

void TSlave::Close(Option_t *opt)
{
   if (fSocket) {

      // If local client ...
      if (!(fProof->IsMaster()) && !strncasecmp(opt,"S",1)) {
         // ... tell master and slaves to stop
         Interrupt(TProof::kShutdownInterrupt);
      }

      // deactivate used sec context if talking to proofd daemon running
      // an old protocol (sec context disactivated remotely)
      TSecContext *sc = fSocket->GetSecContext();
      if (sc && sc->IsActive()) {
         TIter last(sc->GetSecContextCleanup(), kIterBackward);
         TSecContextCleanup *nscc = 0;
         while ((nscc = (TSecContextCleanup *)last())) {
            if (nscc->GetType() == TSocket::kPROOFD &&
                nscc->GetProtocol() < 9) {
               sc->DeActivate("");
               break;
            }
         }
      }
   }

   SafeDelete(fInput);
   SafeDelete(fSocket);
}

////////////////////////////////////////////////////////////////////////////////
/// Used to sort slaves by performance index.

Int_t TSlave::Compare(const TObject *obj) const
{
   const TSlave *sl = dynamic_cast<const TSlave*>(obj);

   if (!sl) {
      Error("Compare", "input is not a TSlave object");
      return 0;
   }

   if (fPerfIdx > sl->GetPerfIdx()) return 1;
   if (fPerfIdx < sl->GetPerfIdx()) return -1;
   const char *myord = GetOrdinal();
   const char *otherord = sl->GetOrdinal();
   while (myord && otherord) {
      Int_t myval = atoi(myord);
      Int_t otherval = atoi(otherord);
      if (myval < otherval) return 1;
      if (myval > otherval) return -1;
      myord = strchr(myord, '.');
      if (myord) myord++;
      otherord = strchr(otherord, '.');
      if (otherord) otherord++;
   }
   if (myord) return -1;
   if (otherord) return 1;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Printf info about slave.

void TSlave::Print(Option_t *) const
{
   TString sc;

   const char *sst[] = { "invalid" , "valid", "inactive" };
   Int_t st = fSocket ? ((fStatus == kInactive) ? 2 : 1) : 0;

   Printf("*** Worker %s  (%s)", fOrdinal.Data(), sst[st]);
   Printf("    Host name:               %s", GetName());
   Printf("    Port number:             %d", GetPort());
   Printf("    Worker session tag:      %s", GetSessionTag());
   Printf("    ROOT version|rev|tag:    %s", GetROOTVersion());
   Printf("    Architecture-Compiler:   %s", GetArchCompiler());
   if (fSocket) {
      if (strlen(GetGroup()) > 0) {
         Printf("    User/Group:              %s/%s", GetUser(), GetGroup());
      } else {
         Printf("    User:                    %s", GetUser());
      }
      if (fSocket->GetSecContext())
         Printf("    Security context:        %s", fSocket->GetSecContext()->AsString(sc));
      Printf("    Proofd protocol version: %d", fSocket->GetRemoteProtocol());
      Printf("    Image name:              %s", GetImage());
      Printf("    Working directory:       %s", GetWorkDir());
      Printf("    Performance index:       %d", GetPerfIdx());
      Printf("    MB's processed:          %.2f", float(GetBytesRead())/(1024*1024));
      Printf("    MB's sent:               %.2f", float(fSocket->GetBytesRecv())/(1024*1024));
      Printf("    MB's received:           %.2f", float(fSocket->GetBytesSent())/(1024*1024));
      Printf("    Real time used (s):      %.3f", GetRealTime());
      Printf("    CPU time used (s):       %.3f", GetCpuTime());
   } else {
      if (strlen(GetGroup()) > 0) {
         Printf("    User/Group:              %s/%s", GetUser(), GetGroup());
      } else {
         Printf("    User:                    %s", GetUser());
      }
      Printf("    Security context:");
      Printf("    Proofd protocol version:");
      Printf("    Image name:              %s", GetImage());
      Printf("    Working directory:       %s", GetWorkDir());
      Printf("    Performance index:       %d", GetPerfIdx());
      Printf("    MB's processed:          %.2f", float(GetBytesRead())/(1024*1024));
      Printf("    MB's sent:");
      Printf("    MB's received:");
      Printf("    Real time used (s):      %.3f", GetRealTime());
      Printf("    CPU time used (s):       %.3f", GetCpuTime());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Adopt and register input handler for this slave. Handler will be deleted
/// by the slave.

void TSlave::SetInputHandler(TFileHandler *ih)
{
   fInput = ih;
   fInput->Add();
}

////////////////////////////////////////////////////////////////////////////////
/// Setup authentication related stuff for old versions.
/// Provided for backward compatibility.

Int_t TSlave::OldAuthSetup(Bool_t master, TString wconf)
{
   static OldSlaveAuthSetup_t oldAuthSetupHook = 0;

   if (!oldAuthSetupHook) {
      // Load libraries needed for (server) authentication ...
      TString authlib = "libRootAuth";
      char *p = 0;
      // The generic one
      if ((p = gSystem->DynamicPathName(authlib, kTRUE))) {
         delete[] p;
         if (gSystem->Load(authlib) == -1) {
            Error("OldAuthSetup", "can't load %s",authlib.Data());
            return kFALSE;
         }
      } else {
         Error("OldAuthSetup", "can't locate %s",authlib.Data());
         return -1;
      }
      //
      // Locate OldSlaveAuthSetup
      Func_t f = gSystem->DynFindSymbol(authlib,"OldSlaveAuthSetup");
      if (f)
         oldAuthSetupHook = (OldSlaveAuthSetup_t)(f);
      else {
         Error("OldAuthSetup", "can't find OldSlaveAuthSetup");
         return -1;
      }
   }
   //
   // Setup
   return (*oldAuthSetupHook)(fSocket, master, fOrdinal, wconf);
}

////////////////////////////////////////////////////////////////////////////////
/// Static method returning the appropriate TSlave object for the remote
/// server.

TSlave *TSlave::Create(const char *url, const char *ord, Int_t perf,
                       const char *image, TProof *proof, Int_t stype,
                       const char *workdir, const char *msd, Int_t nwk)
{
   TSlave *s = 0;

   // Check if we are setting up a lite version
   if (!strcmp(url, "lite")) {
      return new TSlaveLite(ord, perf, image, proof, stype, workdir, msd);
   }

   // No need to try a XPD connection in some well defined cases
   Bool_t tryxpd = kTRUE;
   if (!(proof->IsMaster())) {
      if (proof->IsProofd())
         tryxpd = kFALSE;
   } else {
      if (gApplication && (gApplication->Argc() < 3 ||
                          (gApplication->Argc() > 2 && gApplication->Argv(2) &&
                           strncmp(gApplication->Argv(2),"xpd",3))))
         tryxpd = kFALSE;
   }

   // We do this without the plugin manager because it blocks the CINT mutex
   // breaking the parallel startup
   if (!fgTXSlaveHook) {

      // Load the library containing TXSlave ...
      TString proofxlib = "libProofx";
      char *p = 0;
      if ((p = gSystem->DynamicPathName(proofxlib, kTRUE))) {
         delete[] p;
         if (gSystem->Load(proofxlib) == -1)
            ::Error("TSlave::Create", "can't load %s", proofxlib.Data());
      } else
         ::Error("TSlave::Create", "can't locate %s", proofxlib.Data());
   }

   // Load the right class
   if (fgTXSlaveHook && tryxpd) {
      s = (*fgTXSlaveHook)(url, ord, perf, image, proof, stype, workdir, msd, nwk);
   } else {
      s = new TSlave(url, ord, perf, image, proof, stype, workdir, msd);
   }

   return s;
}

////////////////////////////////////////////////////////////////////////////////
/// Ping the remote master or slave servers.
/// Returns 0 if ok, -1 in case of error

Int_t TSlave::Ping()
{
   if (!IsValid()) return -1;

   TMessage mess(kPROOF_PING | kMESS_ACK);
   fSocket->Send(mess);
   if (fSocket->Send(mess) == -1) {
      Warning("Ping","%s: acknowledgement not received", GetOrdinal());
      return -1;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Send interrupt OOB byte to master or slave servers.
/// Returns 0 if ok, -1 in case of error

void TSlave::Interrupt(Int_t type)
{
   if (!IsValid()) return;

   char oobc = (char) type;
   const int kBufSize = 1024;
   char waste[kBufSize];

   // Send one byte out-of-band message to server
   if (fSocket->SendRaw(&oobc, 1, kOob) <= 0) {
      Error("Interrupt", "error sending oobc to slave %s", GetOrdinal());
      return;
   }

   if (type == TProof::kHardInterrupt) {
      char  oob_byte;
      int   n, nch, nbytes = 0, nloop = 0;

      // Receive the OOB byte
      while ((n = fSocket->RecvRaw(&oob_byte, 1, kOob)) < 0) {
         if (n == -2) {   // EWOULDBLOCK
            //
            // The OOB data has not yet arrived: flush the input stream
            //
            // In some systems (Solaris) regular recv() does not return upon
            // receipt of the oob byte, which makes the below call to recv()
            // block indefinitely if there are no other data in the queue.
            // FIONREAD ioctl can be used to check if there are actually any
            // data to be flushed.  If not, wait for a while for the oob byte
            // to arrive and try to read it again.
            //
            fSocket->GetOption(kBytesToRead, nch);
            if (nch == 0) {
               gSystem->Sleep(1000);
               continue;
            }

            if (nch > kBufSize) nch = kBufSize;
            n = fSocket->RecvRaw(waste, nch);
            if (n <= 0) {
               Error("Interrupt", "error receiving waste from slave %s",
                     GetOrdinal());
               break;
            }
            nbytes += n;
         } else if (n == -3) {   // EINVAL
            //
            // The OOB data has not arrived yet
            //
            gSystem->Sleep(100);
            if (++nloop > 100) {  // 10 seconds time-out
               Error("Interrupt", "server %s does not respond", GetOrdinal());
               break;
            }
         } else {
            Error("Interrupt", "error receiving OOB from server %s",
                  GetOrdinal());
            break;
         }
      }

      //
      // Continue flushing the input socket stream until the OOB
      // mark is reached
      //
      while (1) {
         int atmark;

         fSocket->GetOption(kAtMark, atmark);

         if (atmark)
            break;

         // find out number of bytes to read before atmark
         fSocket->GetOption(kBytesToRead, nch);
         if (nch == 0) {
            gSystem->Sleep(1000);
            continue;
         }

         if (nch > kBufSize) nch = kBufSize;
         n = fSocket->RecvRaw(waste, nch);
         if (n <= 0) {
            Error("Interrupt", "error receiving waste (2) from slave %s",
                  GetOrdinal());
            break;
         }
         nbytes += n;
      }
      if (nbytes > 0) {
         if (fProof->IsMaster())
            Info("Interrupt", "slave %s:%s synchronized: %d bytes discarded",
                 GetName(), GetOrdinal(), nbytes);
         else
            Info("Interrupt", "PROOF synchronized: %d bytes discarded", nbytes);
      }

      // Get log file from master or slave after a hard interrupt
      fProof->Collect(this);

   } else if (type == TProof::kSoftInterrupt) {

      // Get log file from master or slave after a soft interrupt
      fProof->Collect(this);

   } else if (type == TProof::kShutdownInterrupt) {

      ; // nothing expected to be returned

   } else {

      // Unexpected message, just receive log file
      fProof->Collect(this);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Sent stop/abort request to PROOF server.

void TSlave::StopProcess(Bool_t abort, Int_t timeout)
{
   // Notify the remote counterpart
   TMessage msg(kPROOF_STOPPROCESS);
   msg << abort;
   if (fProof->fProtocol > 9)
      msg << timeout;
   fSocket->Send(msg);
}

////////////////////////////////////////////////////////////////////////////////
/// Send message to intermediate coordinator. Only meaningful when there is one,
/// i.e. in XPD framework

TObjString *TSlave::SendCoordinator(Int_t, const char *, Int_t)
{
   if (gDebug > 0)
      Info("SendCoordinator","method not implemented for this communication layer");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Set an alias for this session. If reconnection is supported, the alias
/// will be communicated to the remote coordinator so that it can be recovered
/// when reconnecting

void TSlave::SetAlias(const char *)
{
   if (gDebug > 0)
      Info("SetAlias","method not implemented for this communication layer");
   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Set hook to TXSlave ctor

void TSlave::SetTXSlaveHook(TSlave_t xslavehook)
{
   fgTXSlaveHook = xslavehook;
}
