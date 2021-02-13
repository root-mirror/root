// @(#)root/proof:$Id$
// Author: Gerardo Ganis  March 2008

/*************************************************************************
 * Copyright (C) 1995-2005, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TSlaveLite
\ingroup proofkernel

Version of TSlave for local worker servers.
See TSlave for details.

*/



//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TSlaveLite                                                           //
//                                                                      //
// This is the version of TSlave for local worker servers.              //
// See TSlave for details.                                              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "RConfigure.h"
#include "TSlaveLite.h"
#include "TProof.h"
#include "TProofServ.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TUrl.h"
#include "TMessage.h"
#include "TMonitor.h"
#include "TError.h"
#include "TSocket.h"
#include "TSysEvtHandler.h"
#include "TVirtualMutex.h"

ClassImp(TSlaveLite);

//______________________________________________________________________________
//---- error handling ----------------------------------------------------------
//---- Needed to avoid blocking on the CINT mutex in printouts -----------------

////////////////////////////////////////////////////////////////////////////////
/// Interface to ErrorHandler (protected).

void TSlaveLite::DoError(int level, const char *location,
                                    const char *fmt, va_list va) const
{
   ::ErrorHandler(level, Form("TSlaveLite::%s", location), fmt, va);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a PROOF slave object. Called via the TProof ctor.

TSlaveLite::TSlaveLite(const char *ord, Int_t perf,
               const char *image, TProof *proof, Int_t stype,
               const char *workdir, const char *msd, Int_t) : TSlave()
{
   fName = ord;  // Need this during the setup phase; see end of SetupServ
   fImage = image;
   fProofWorkDir = workdir;
   fWorkDir = workdir;
   fOrdinal = ord;
   fPerfIdx = perf;
   fProof = proof;
   fSlaveType = (ESlaveType)stype;
   fMsd = msd;
   fIntHandler = 0;
   fValid = kFALSE;
   fProtocol = kPROOF_Protocol;

   if (fPerfIdx > 0) Init();
}

////////////////////////////////////////////////////////////////////////////////
/// Init a PROOF worker object. Called via the TSlaveLite ctor.

void TSlaveLite::Init()
{
   // Command to be executed
   TString cmd;
   cmd.Form(". %s/worker-%s.env; export ROOTBINDIR=\"%s\"; %s/proofserv proofslave lite %d %d 0&",
            fWorkDir.Data(), fOrdinal.Data(), TROOT::GetBinDir().Data(), TROOT::GetBinDir().Data(),
            gSystem->GetPid(), gDebug);
   // Execute
   if (gSystem->Exec(cmd) != 0) {
      Error("Init", "an error occured while executing 'proofserv'");
      SetBit(kInvalidObject);
      return;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Init a PROOF slave object. Called via the TSlaveLite ctor.
/// The Init method is technology specific and is overwritten by derived
/// classes.

Int_t TSlaveLite::SetupServ(Int_t, const char *)
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

   // Receive the unique tag and save it as name of this object
   TMessage *msg = 0;
   if (fSocket->Recv(msg) <= 0 || !msg || msg->What() != kPROOF_SESSIONTAG) {
      Error("SetupServ", "failed to receive unique session tag");
      Close("S");
      SafeDelete(fSocket);
      fValid = kFALSE;
      return -1;
   }
   // Extract the unique tag
   (*msg) >> fSessionTag;

   // Set the real name (temporarly set to ordinal for the setup)
   fName = gSystem->HostName();

   // We are done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Destroy slave.

TSlaveLite::~TSlaveLite()
{
   Close();
}

////////////////////////////////////////////////////////////////////////////////
/// Close slave socket.

void TSlaveLite::Close(Option_t *opt)
{
   if (fSocket)
      // Closing socket ...
      fSocket->Close(opt);

   SafeDelete(fInput);
   SafeDelete(fSocket);
}

////////////////////////////////////////////////////////////////////////////////
/// Printf info about slave.

void TSlaveLite::Print(Option_t *) const
{
   const char *sst[] = { "invalid" , "valid", "inactive" };
   Int_t st = fSocket ? ((fStatus == kInactive) ? 2 : 1) : 0;

   Printf("*** Worker %s  (%s)", fOrdinal.Data(), sst[st]);
   Printf("    Worker session tag:      %s", GetSessionTag());
   Printf("    ROOT version|rev|tag:    %s", GetROOTVersion());
   Printf("    Architecture-Compiler:   %s", GetArchCompiler());
   if (fSocket) {
      Printf("    Working directory:       %s", GetWorkDir());
      Printf("    MB's processed:          %.2f", float(GetBytesRead())/(1024*1024));
      Printf("    MB's sent:               %.2f", float(fSocket->GetBytesRecv())/(1024*1024));
      Printf("    MB's received:           %.2f", float(fSocket->GetBytesSent())/(1024*1024));
      Printf("    Real time used (s):      %.3f", GetRealTime());
      Printf("    CPU time used (s):       %.3f", GetCpuTime());
   }
}
