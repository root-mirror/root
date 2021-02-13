// @(#)root/proof:$Id$
// Author: Maarten Ballintijn   06/12/03

/*************************************************************************
 * Copyright (C) 1995-2001, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TCondor                                                              //
//                                                                      //
// Interface to the Condor system. TCondor provides a (partial) API for //
// querying and controlling the Condor system, including experimental   //
// extensions like COD (computing on demand)                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <stdlib.h>

#include "TCondor.h"
#include "TList.h"
#include "TSystem.h"
#include "TObjString.h"
#include "TRegexp.h"
#include "TProofDebug.h"
#include "Riostream.h"
#include "TEnv.h"
#include "TClass.h"

ClassImp(TCondorSlave);
ClassImp(TCondor);


////////////////////////////////////////////////////////////////////////////////
/// Create Condor interface object. Uses Condor apps since there is no
/// API yet.

TCondor::TCondor(const char *pool) : fPool(pool), fState(kFree)
{
   fClaims = new TList;

   // Setup Condor

   TString condorHome = gEnv->GetValue("Proof.CondorHome", (char*)0);
   if (condorHome != "") {
      TString path = gSystem->Getenv("PATH");
      path = condorHome + "/bin:" + path;
      gSystem->Setenv("PATH",path);
   }

   TString condorConf = gEnv->GetValue("Proof.CondorConfig", (char*)0);
   if (condorConf != "") {
      gSystem->Setenv("CONDOR_CONFIG",condorConf);
   }

   char *loc = gSystem->Which(gSystem->Getenv("PATH"), "condor_cod",
                                                kExecutePermission);

   if (loc) {
      fValid = kTRUE;
      delete [] loc;
   } else {
      fValid = kFALSE;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Cleanup Condor interface.

TCondor::~TCondor()
{
   PDB(kCondor,1) Info("~TCondor","fState %d", fState );

   if (fState != kFree) {
      Release();
   }
   delete fClaims;
}


////////////////////////////////////////////////////////////////////////////////
/// Print master status

void TCondor::Print(Option_t * opt) const
{
   std::cout << "OBJ: " << IsA()->GetName()
      << "\tPool: \"" << fPool << "\""
      << "\tState: " << fState << std::endl;
   fClaims->Print(opt);
}


////////////////////////////////////////////////////////////////////////////////
/// Claim a VirtualMachine for PROOF usage.

TCondorSlave *TCondor::ClaimVM(const char *vm, const char *cmd)
{
   Int_t port = 0;

   TString claimCmd = Form("condor_cod request -name %s -timeout 10 2>>%s/condor.proof.%d",
                           vm, gSystem->TempDirectory(), gSystem->GetUid() );

   PDB(kCondor,2) Info("ClaimVM","command: %s", claimCmd.Data());
   FILE  *pipe = gSystem->OpenPipe(claimCmd, "r");

   if (!pipe) {
      SysError("ClaimVM","cannot run command: %s", claimCmd.Data());
      return 0;
   }

   TString claimId;
   TString line;
   while (line.Gets(pipe)) {
      PDB(kCondor,3) Info("ClaimVM","line = %s", line.Data());

      if (line.BeginsWith("ClaimId = \"")) {
         line.Remove(0, line.Index("\"")+1);
         line.Chop(); // remove trailing "
         claimId = line;
         PDB(kCondor,1) Info("ClaimVM","claim = '%s'", claimId.Data());
         TRegexp r("[0-9]*$");
         TString num = line(r);
         port = 37000 + atoi(num.Data());
         PDB(kCondor,1) Info("ClaimVM","port = %d", port);
      }
   }

   Int_t r = gSystem->ClosePipe(pipe);
   if (r) {
      Error("ClaimVM","command: %s returned %d", claimCmd.Data(), r);
      return 0;
   } else {
      PDB(kCondor,1) Info("ClaimVM","command: %s returned %d", claimCmd.Data(), r);
   }

   TString jobad("jobad");
   FILE *jf = gSystem->TempFileName(jobad);

   if (jf == 0) return 0;

   TString str(cmd);
   str.ReplaceAll("$(Port)", Form("%d", port));
   fputs(str, jf);

   fclose(jf);

   TString activateCmd = Form("condor_cod activate -id '%s' -jobad %s",
                              claimId.Data(), jobad.Data() );

   PDB(kCondor,2) Info("ClaimVM","command: %s", activateCmd.Data());
   pipe = gSystem->OpenPipe(activateCmd, "r");

   if (!pipe) {
      SysError("ClaimVM","cannot run command: %s", activateCmd.Data());
      return 0;
   }

   while (line.Gets(pipe)) {
      PDB(kCondor,3) Info("ClaimVM","Activate: line = %s", line.Data());
   }

   r = gSystem->ClosePipe(pipe);
   if (r) {
      Error("ClaimVM","command: %s returned %d", activateCmd.Data(), r);
   } else {
      PDB(kCondor,1) Info("ClaimVM","command: %s returned %d", activateCmd.Data(), r);
   }

   gSystem->Unlink(jobad);

   // TODO: get info at the start for all nodes ...
   TCondorSlave *claim = new TCondorSlave;
   claim->fClaimID = claimId;
   TString node(vm);
   node = node.Remove(0, node.Index("@")+1);
   claim->fHostname = node;
   claim->fPort = port;
   claim->fPerfIdx = 100; //set performance index to 100 by default
   claim->fImage = node; //set image to hostname by default

   return claim;
}


////////////////////////////////////////////////////////////////////////////////
/// Get the names of the virtual machines in the pool.
/// Return a TList of TObjString or 0 in case of failure

TList *TCondor::GetVirtualMachines() const
{
   TString poolopt = fPool ? Form("-pool %s", fPool.Data()) : "";
   TString cmd = Form("condor_status %s -format \"%%s\\n\" Name", poolopt.Data());

   PDB(kCondor,2) Info("GetVirtualMachines","command: %s", cmd.Data());

   FILE  *pipe = gSystem->OpenPipe(cmd, "r");

   if (!pipe) {
      SysError("GetVirtualMachines","cannot run command: %s", cmd.Data());
      return 0;
   }

   TString line;
   TList *l = new TList;
   while (line.Gets(pipe)) {
      PDB(kCondor,3) Info("GetVirtualMachines","line = %s", line.Data());
      if (line != "") l->Add(new TObjString(line));
   }

   Int_t r = gSystem->ClosePipe(pipe);
   if (r) {
      delete l;
      Error("GetVirtualMachines","command: %s returned %d", cmd.Data(), r);
      return 0;
   } else {
      PDB(kCondor,1) Info("GetVirtualMachines","command: %s returned %d", cmd.Data(), r);
   }

   return l;
}


////////////////////////////////////////////////////////////////////////////////
/// Claim n virtual machines
/// This function figures out the image and performance index before returning
/// the list of condor slaves

TList *TCondor::Claim(Int_t n, const char *cmd)
{
   if (fState != kFree) {
      Error("Claim","not in state Free");
      return 0;
   }

   TList *vms = GetVirtualMachines();
   TIter next(vms);
   TObjString *vm;
   for(Int_t i=0; i < n && (vm = (TObjString*) next()) != 0; i++ ) {
      TCondorSlave *claim = ClaimVM(vm->GetName(), cmd);
      if (claim != 0) {
         if ( !GetVmInfo(vm->GetName(), claim->fImage, claim->fPerfIdx) ) {
            // assume vm is gone
            delete claim;
         } else {
            fClaims->Add(claim);
            fState = kActive;
         }
      }
   }

   vms->Delete();
   delete vms;

   return fClaims;
}


////////////////////////////////////////////////////////////////////////////////
/// Claim virtual machine with name vmname
/// This function does not figure out the image and performance index before
/// returning the condor slave

TCondorSlave *TCondor::Claim(const char *vmname, const char *cmd)
{
   if (fState != kFree && fState != kActive) {
      Error("Claim","not in state Free or Active");
      return 0;
   }

   TCondorSlave *claim = ClaimVM(vmname, cmd);
   if (claim != 0) {
      fClaims->Add(claim);
      fState = kActive;
   }

   return claim;
}


////////////////////////////////////////////////////////////////////////////////
/// Set the state of workers

Bool_t TCondor::SetState(EState state)
{
   PDB(kCondor,1) Info("SetState","state: %s (%lld)",
                       state == kSuspended ? "kSuspended" : "kActive", Long64_t(gSystem->Now()));
   TIter next(fClaims);
   TCondorSlave *claim;
   while((claim = (TCondorSlave*) next()) != 0) {
      TString cmd = Form("condor_cod %s -id '%s'",
                         state == kSuspended ? "suspend" : "resume",
                         claim->fClaimID.Data());

      PDB(kCondor,2) Info("SetState","command: %s", cmd.Data());
      FILE  *pipe = gSystem->OpenPipe(cmd, "r");

      if (!pipe) {
         SysError("SetState","cannot run command: %s", cmd.Data());
         return kFALSE;
      }

      TString line;
      while (line.Gets(pipe)) {
         PDB(kCondor,3) Info("SetState","line = %s", line.Data());
      }

      Int_t r = gSystem->ClosePipe(pipe);
      if (r) {
         Error("SetState","command: %s returned %d", cmd.Data(), r);
         return kFALSE;
      } else {
         PDB(kCondor,1) Info("SetState","command: %s returned %d", cmd.Data(), r);
      }
   }

   fState = state;
   return kTRUE;
}


////////////////////////////////////////////////////////////////////////////////
/// Suspend worker

Bool_t TCondor::Suspend()
{
   if (fState != kActive) {
      Error("Suspend","not in state Active");
      return kFALSE;
   }

   return SetState(kSuspended);
}


////////////////////////////////////////////////////////////////////////////////
/// Resume worker

Bool_t TCondor::Resume()
{
   if (fState != kSuspended) {
      Error("Suspend","not in state Suspended");
      return kFALSE;
   }

   return SetState(kActive);
}


////////////////////////////////////////////////////////////////////////////////
/// Release worker

Bool_t TCondor::Release()
{
   if (fState == kFree) {
      Error("Suspend","not in state Active or Suspended");
      return kFALSE;
   }

   TCondorSlave *claim;
   while((claim = (TCondorSlave*) fClaims->First()) != 0) {
      TString cmd = Form("condor_cod release -id '%s'", claim->fClaimID.Data());

      PDB(kCondor,2) Info("SetState","command: %s", cmd.Data());
      FILE  *pipe = gSystem->OpenPipe(cmd, "r");

      if (!pipe) {
         SysError("Release","cannot run command: %s", cmd.Data());
         return kFALSE;
      }

      TString line;
      while (line.Gets(pipe)) {
         PDB(kCondor,3) Info("Release","line = %s", line.Data());
      }

      Int_t r = gSystem->ClosePipe(pipe);
      if (r) {
         Error("Release","command: %s returned %d", cmd.Data(), r);
         return kFALSE;
      } else {
         PDB(kCondor,1) Info("Release","command: %s returned %d", cmd.Data(), r);
      }

      fClaims->Remove(claim);
      delete claim;
   }

   fState = kFree;
   return kTRUE;
}


////////////////////////////////////////////////////////////////////////////////
/// Get info about worker status

Bool_t TCondor::GetVmInfo(const char *vm, TString &image, Int_t &perfidx) const
{
   TString cmd = Form("condor_status -format \"%%d:\" Mips -format \"%%s\\n\" FileSystemDomain "
                      "-const 'Name==\"%s\"'", vm);

   PDB(kCondor,2) Info("GetVmInfo","command: %s", cmd.Data());
   FILE  *pipe = gSystem->OpenPipe(cmd, "r");

   if (!pipe) {
      SysError("GetVmInfo","cannot run command: %s", cmd.Data());
      return kFALSE;
   }

   TString line;
   while (line.Gets(pipe)) {
      PDB(kCondor,3) Info("GetVmInfo","line = %s", line.Data());
      if (line != "") {
         TString amips = line(TRegexp("^[0-9]*"));
         perfidx = atoi(amips);
         image = line(TRegexp("[^:]+$"));
         break;
      }
   }

   Int_t r = gSystem->ClosePipe(pipe);
   if (r) {
      Error("GetVmInfo","command: %s returned %d", cmd.Data(), r);
      return kFALSE;
   } else {
      PDB(kCondor,1) Info("GetVmInfo","command: %s returned %d", cmd.Data(), r);
   }

   return kTRUE;
}


////////////////////////////////////////////////////////////////////////////////
/// Get image of the worker

TString TCondor::GetImage(const char *host) const
{
   TString cmd = Form("condor_status -direct %s -format \"Image:%%s\\n\" "
                      "FileSystemDomain", host);

   PDB(kCondor,2) Info("GetImage","command: %s", cmd.Data());

   FILE  *pipe = gSystem->OpenPipe(cmd, "r");

   if (!pipe) {
      SysError("GetImage","cannot run command: %s", cmd.Data());
      return "";
   }

   TString image;
   TString line;
   while (line.Gets(pipe)) {
      PDB(kCondor,3) Info("GetImage","line = %s", line.Data());
      if (line != "") {
         image = line(TRegexp("[^:]+$"));
         break;
      }
   }

   Int_t r = gSystem->ClosePipe(pipe);
   if (r) {
      Error("GetImage","command: %s returned %d", cmd.Data(), r);
      return "";
   } else {
      PDB(kCondor,1) Info("GetImage","command: %s returned %d", cmd.Data(), r);
   }

   return image;
}


////////////////////////////////////////////////////////////////////////////////
/// Print worker status

void TCondorSlave::Print(Option_t * /*opt*/ ) const
{
   std::cout << "OBJ: " << IsA()->GetName()
      << " " << fHostname << ":" << fPort
      << "  Perf: " << fPerfIdx
      << "  Image: " << fImage << std::endl;
}
