// @(#)root/proof:$Id$
// Author: Sangsu Ryu 28/06/2011

/*************************************************************************
 * Copyright (C) 1995-2005, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TSelVerifyDataSet
\ingroup proofkernel

Selector to verify dataset in parallel on workers

*/

#define TSelVerifyDataSet_cxx

#include "TSelVerifyDataSet.h"
#include "TDataSetManager.h"
#include "TDSet.h"
#include "TParameter.h"
#include "TTree.h"
#include "TFile.h"
#include "TNamed.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TFileStager.h"
#include "TProofDebug.h"
#include "TProofServ.h"
#include "TFileCollection.h"
#include "TFileInfo.h"

ClassImp(TSelVerifyDataSet);

////////////////////////////////////////////////////////////////////////////////
/// Constructor

TSelVerifyDataSet::TSelVerifyDataSet(TTree *)
{
   InitMembers();
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor

TSelVerifyDataSet::TSelVerifyDataSet()
{
   InitMembers();
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize members

void TSelVerifyDataSet::InitMembers()
{
   fFopt = -1;
   fSopt = 0;
   fRopt = 0;

   fAllf = 0;
   fCheckstg = 0;
   fNonStgf = 0;
   fReopen = 0;
   fTouch = 0;
   fStgf = 0;
   fNoaction = 0;
   fFullproc = 0;
   fLocateonly = 0;
   fStageonly = 0;
   fDoall       = 0;
   fGetlistonly = 0;
   fScanlist    = 0;
   fDbg = 0;

   fChangedDs = kFALSE;
   fTouched = 0;
   fOpened = 0;
   fDisappeared = 0;
   fSubDataSet = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Worker Begin

void TSelVerifyDataSet::SlaveBegin(TTree *)
{
   TString dsname, opts;

   TNamed* par = dynamic_cast<TNamed*>(fInput->FindObject("PROOF_VerifyDataSet"));
   if (par) {
      dsname = par->GetTitle();
   } else {
      Abort("cannot find dataset name: cannot continue", kAbortProcess);
      return;
   }

   par = dynamic_cast<TNamed*>(fInput->FindObject("PROOF_VerifyDataSetOption"));
   if (par) {
      opts = par->GetTitle();
   } else {
      Abort("cannot find verify options: cannot continue", kAbortProcess);
      return;
   }

   par = dynamic_cast<TNamed*>(fInput->FindObject("PROOF_MSS"));
   if (par) {
      fMss = par->GetTitle();
      PDB(kSelector, 2) Info("SlaveBegin", "dataset MSS: '%s'", fMss.Data());
   }

   par = dynamic_cast<TNamed*>(fInput->FindObject("PROOF_StageOption"));
   if (par) {
      fStageopts = par->GetTitle();
      PDB(kSelector, 2) Info("SlaveBegin", "dataset stage options: '%s'", fStageopts.Data());
   }

   // Extract the directives
   UInt_t o = 0;
   if (!opts.IsNull()) {
      // Selection options
      if (strstr(opts, "allfiles:") || strchr(opts, 'A'))
         o |= TDataSetManager::kAllFiles;
      else if (strstr(opts, "staged:") || strchr(opts, 'D'))
         o |= TDataSetManager::kStagedFiles;
      // Pre-action options
      if (strstr(opts, "open:") || strchr(opts, 'O'))
         o |= TDataSetManager::kReopen;
      if (strstr(opts, "touch:") || strchr(opts, 'T'))
         o |= TDataSetManager::kTouch;
      if (strstr(opts, "nostagedcheck:") || strchr(opts, 'I'))
         o |= TDataSetManager::kNoStagedCheck;
      // Process options
      if (strstr(opts, "noaction:") || strchr(opts, 'N'))
         o |= TDataSetManager::kNoAction;
      if (strstr(opts, "locateonly:") || strchr(opts, 'L'))
         o |= TDataSetManager::kLocateOnly;
      if (strstr(opts, "stageonly:") || strchr(opts, 'S'))
         o |= TDataSetManager::kStageOnly;
      // Auxilliary options
      if (strstr(opts, "verbose:") || strchr(opts, 'V'))
         o |= TDataSetManager::kDebug;
   } else {
      // Default
      o = TDataSetManager::kReopen | TDataSetManager::kDebug;
   }

   PDB(kSelector, 1) Info("SlaveBegin", "o=%d", o);
   // File selection
   fFopt = ((o & TDataSetManager::kAllFiles)) ? -1 : 0;
   if (fFopt >= 0) {
      if ((o & TDataSetManager::kStagedFiles)) {
         fFopt = 10;
      } else {
         if ((o & TDataSetManager::kReopen)) fFopt++;
         if ((o & TDataSetManager::kTouch)) fFopt++;
      }
      if ((o & TDataSetManager::kNoStagedCheck)) fFopt += 100;
   } else {
      if ((o & TDataSetManager::kStagedFiles) || (o & TDataSetManager::kReopen) || (o & TDataSetManager::kTouch)) {
         Warning("SlaveBegin", "kAllFiles mode: ignoring kStagedFiles or kReopen"
                                " or kTouch requests");
      }
      if ((o & TDataSetManager::kNoStagedCheck)) fFopt -= 100;
   }
   PDB(kSelector, 1) Info("SlaveBegin", "fFopt=%d", fFopt);

   // Type of action
   fSopt = ((o & TDataSetManager::kNoAction)) ? -1 : 0;
   if (fSopt >= 0) {
      if ((o & TDataSetManager::kLocateOnly) && (o & TDataSetManager::kStageOnly)) {
         Error("SlaveBegin", "kLocateOnly and kStageOnly cannot be processed concurrently");
         return;
      }
      if ((o & TDataSetManager::kLocateOnly)) fSopt = 1;
      if ((o & TDataSetManager::kStageOnly)) fSopt = 2;
   } else if ((o & TDataSetManager::kLocateOnly) || (o & TDataSetManager::kStageOnly)) {
      Warning("SlaveBegin", "kNoAction mode: ignoring kLocateOnly or kStageOnly requests");
   }
   PDB(kSelector, 1) Info("SlaveBegin", "fSopt=%d", fSopt);

   fDbg = ((o & TDataSetManager::kDebug)) ? kTRUE : kFALSE;

   // File selection, Reopen and Touch options
   fAllf     = (fFopt == -1)               ? kTRUE : kFALSE;
   fCheckstg = (fFopt >= 100 || fFopt < -1) ? kFALSE : kTRUE;
   if (fFopt >= 0) fFopt %= 100;
   fNonStgf  = (fFopt >= 0 && fFopt < 10)   ? kTRUE : kFALSE;
   fReopen   = (fFopt >= 1 && fFopt < 10)   ? kTRUE : kFALSE;
   fTouch    = (fFopt >= 2 && fFopt < 10)   ? kTRUE : kFALSE;
   fStgf     = (fFopt == 10)               ? kTRUE : kFALSE;

   PDB(kSelector, 1) Info("SlaveBegin",
                          "fAllf=%d fCheckstg=%d fNonStgf=%d fReopen=%d fTouch=%d fStgf=%d",
                          fAllf, fCheckstg, fNonStgf, fReopen, fTouch, fStgf);

   // File processing options
   fNoaction   = (fSopt == -1) ? kTRUE : kFALSE;
   fFullproc   = (fSopt == 0)  ? kTRUE : kFALSE;
   fLocateonly = (fSopt == 1)  ? kTRUE : kFALSE;
   fStageonly  = (fSopt == 2)  ? kTRUE : kFALSE;

   PDB(kSelector, 1) Info("SlaveBegin",
                          "fNoaction=%d fFullproc=%d fLocateonly=%d fStageonly=%d",
                          fNoaction, fFullproc, fLocateonly, fStageonly);

   // Run options
   fDoall       = (fRopt == 0) ? kTRUE : kFALSE;
   fGetlistonly = (fRopt == 1) ? kTRUE : kFALSE;
   fScanlist    = (fRopt == 2) ? kTRUE : kFALSE;

   PDB(kSelector, 1) Info("SlaveBegin",
                          "fDoall=%d fGetlistonly=%d fScanlist=%d",
                          fDoall, fGetlistonly, fScanlist);

   TString hostname(TUrl(gSystem->HostName()).GetHostFQDN());
   TString thisordinal = gProofServ ? gProofServ->GetOrdinal() : "n.d";
   TString title =
      TString::Format("TSelVerifyDataSet_%s_%s", hostname.Data(), thisordinal.Data());
   fSubDataSet= new TFileCollection(dsname, title);
}

////////////////////////////////////////////////////////////////////////////////
/// Process a single entry

Bool_t TSelVerifyDataSet::Process(Long64_t entry)
{
   TDSetElement *fCurrent = 0;
   TPair *elemPair = 0;
   if (fInput && (elemPair = dynamic_cast<TPair *>
                      (fInput->FindObject("PROOF_CurrentElement")))) {
      if ((fCurrent = dynamic_cast<TDSetElement *>(elemPair->Value())))
         Info("Process", "entry %lld: file: '%s'", entry, fCurrent->GetName());
   }
   if (!fCurrent) {
      Error("Process", "entry %lld: current element not found!", entry);
      return kFALSE;
   }

   TFileInfo *fileInfo = dynamic_cast<TFileInfo*>(fCurrent->GetAssocObj(0));
   if (!fileInfo) {
      Error("Process", "can not get TFileInfo; returning");
      return kFALSE;
   }

   PDB(kSelector, 1) {
      Info("Process", "input fileinfo: ");
      fileInfo->Print("L");
   }

   TFileStager *stager = 0;
   Bool_t createStager = kFALSE;

   TFileInfo* newfileinfo = new TFileInfo(*fileInfo);
   newfileinfo->SetIndex(fileInfo->GetIndex());

   if (fDoall || fGetlistonly) {

      stager = (fMss && strlen(fMss) > 0) ? TFileStager::Open(fMss) : 0;
      createStager = (stager) ? kFALSE : kTRUE;

      // Check which files have been staged, this can be replaced by a bulk command,
      // once it exists in the xrdclient

      // For real time monitoring
      gSystem->DispatchOneEvent(kTRUE);

      Bool_t changed = kFALSE;
      Bool_t touched = kFALSE;
      Bool_t disappeared = kFALSE;

      TDataSetManager::CheckStagedStatus(newfileinfo, fFopt, -1, 0, stager, createStager,
                                         fDbg, changed, touched, disappeared);

      if (changed) fChangedDs = kTRUE;
      if (touched) fTouched++;
      if (disappeared) fDisappeared++;

      SafeDelete(stager);

      PDB(kSelector, 1) Info("Process",
                             "fChangedDs = %d, fTouched = %d disappeared = %d",
                              fChangedDs, fTouched, fDisappeared);

      // If required to only get the list we are done
      if (fGetlistonly) {
         Info("Process", "updated fileinfo: ");
         newfileinfo->Print("F");
         fSubDataSet->Add(newfileinfo);
         return kTRUE;
      }
   }

   if (!fNoaction && (fDoall || fScanlist)) {

      // Point to the fileinfo
      //newStagedFiles = (!fDoall && fScanlist && flist) ? flist : newStagedFiles;
      if (!fDoall && fScanlist) {
          SafeDelete(newfileinfo);
          newfileinfo = new TFileInfo(*fileInfo);
          newfileinfo->SetIndex(fileInfo->GetIndex());
      }

      // Loop over now staged files
      PDB(kSelector, 1) Info("Process",
                             "file appear to be newly staged; %s",
                             newfileinfo->GetFirstUrl()->GetUrl());

      // If staging files, prepare the stager
      if (fLocateonly || fStageonly) {
         stager = (fMss && strlen(fMss) > 0) ? TFileStager::Open(fMss) : 0;
         createStager = (stager) ? kFALSE : kTRUE;
      }

      // Process the file
      Bool_t changed = kFALSE;
      Bool_t opened = kFALSE;
      TDataSetManager::ProcessFile(newfileinfo, fSopt, fCheckstg, fDoall, stager, createStager, fStageopts,
                                   fDbg, changed, opened);

      if (changed) fChangedDs = kTRUE;
      if (opened) fOpened++;
   }

   PDB(kSelector, 1) {
      Info("Process", "updated fileinfo: ");
      newfileinfo->Print("L");
   }
   fSubDataSet->Add(newfileinfo);

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Worker Terminate

void TSelVerifyDataSet::SlaveTerminate()
{
   if (fSubDataSet) {
      fSubDataSet->Update();
      if (fSubDataSet->GetNFiles() > 0) {
         fOutput->Add(fSubDataSet);
         Info("SlaveTerminate",
              "sub-dataset '%s' added to the output list (%lld files)",
              fSubDataSet->GetTitle(), fSubDataSet->GetNFiles());
      }
      // Add information for registration
      fOutput->Add(new TNamed(TString::Format("DATASET_%s", fSubDataSet->GetName()).Data(),"OT:sortidx:"));
      fOutput->Add(new TNamed("PROOFSERV_RegisterDataSet", ""));
   }

   // Send the number of files disppeared, opened and mark 'changed'if any fileinfo in the dataset has changed
   TString hostname(TUrl(gSystem->HostName()).GetHostFQDN());
   TString thisordinal = gProofServ ? gProofServ->GetOrdinal() : "n.d";
   TString sfdisppeared= TString::Format("PROOF_NoFilesDisppeared_%s_%s", hostname.Data(), thisordinal.Data());
   fOutput->Add(new TParameter<Int_t>(sfdisppeared.Data(), fDisappeared));
   TString sfOpened= TString::Format("PROOF_NoFilesOpened_%s_%s", hostname.Data(), thisordinal.Data());
   fOutput->Add(new TParameter<Int_t>(sfOpened.Data(), fOpened));
   TString sfTouched = TString::Format("PROOF_NoFilesTouched_%s_%s", hostname.Data(), thisordinal.Data());
   fOutput->Add(new TParameter<Int_t>(sfTouched.Data(), fTouched));
   TString schanged= TString::Format("PROOF_DataSetChanged_%s_%s", hostname.Data(), thisordinal.Data());
   fOutput->Add(new TParameter<Bool_t>(schanged.Data(), fChangedDs));
}
