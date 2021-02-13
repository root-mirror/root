// @(#)root/proofplayer:$Id$
// Author: Maarten Ballintijn   07/01/02

/*************************************************************************
 * Copyright (C) 1995-2001, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TProofPlayer
\ingroup proofkernel

Internal class steering processing in PROOF.
Instances of the TProofPlayer class are created on the worker nodes
per session and do the processing.
Instances of its subclass - TProofPlayerRemote are created per each
query on the master(s) and on the client. On the master(s),
TProofPlayerRemote coordinate processing, check the dataset, create
the packetizer and take care of merging the results of the workers.
The instance on the client collects information on the input
(dataset and selector), it invokes the Begin() method and finalizes
the query by calling Terminate().

*/

#include "TProofDraw.h"
#include "TProofPlayer.h"
#include "THashList.h"
#include "TEnv.h"
#include "TEventIter.h"
#include "TVirtualPacketizer.h"
#include "TSelector.h"
#include "TSocket.h"
#include "TProofServ.h"
#include "TProof.h"
#include "TProofOutputFile.h"
#include "TProofSuperMaster.h"
#include "TSlave.h"
#include "TClass.h"
#include "TROOT.h"
#include "TError.h"
#include "TException.h"
#include "MessageTypes.h"
#include "TMessage.h"
#include "TDSetProxy.h"
#include "TString.h"
#include "TSystem.h"
#include "TFile.h"
#include "TFileCollection.h"
#include "TFileInfo.h"
#include "TFileMerger.h"
#include "TProofDebug.h"
#include "TTimer.h"
#include "TMap.h"
#include "TPerfStats.h"
#include "TStatus.h"
#include "TEventList.h"
#include "TProofLimitsFinder.h"
#include "TSortedList.h"
#include "TTree.h"
#include "TEntryList.h"
#include "TDSet.h"
#include "TDrawFeedback.h"
#include "TNamed.h"
#include "TObjString.h"
#include "TQueryResult.h"
#include "TMD5.h"
#include "TMethodCall.h"
#include "TObjArray.h"
#include "TH1.h"
#include "TVirtualMonitoring.h"
#include "TParameter.h"
#include "TOutputListSelectorDataMap.h"
#include "TStopwatch.h"

// Timeout exception
#define kPEX_STOPPED  1001
#define kPEX_ABORTED  1002

// To flag an abort condition: use a local static variable to avoid
// warnings about problems with longjumps
static Bool_t gAbort = kFALSE;

class TAutoBinVal : public TNamed {
private:
   Double_t fXmin, fXmax, fYmin, fYmax, fZmin, fZmax;

public:
   TAutoBinVal(const char *name, Double_t xmin, Double_t xmax, Double_t ymin,
               Double_t ymax, Double_t zmin, Double_t zmax) : TNamed(name,"")
   {
      fXmin = xmin; fXmax = xmax;
      fYmin = ymin; fYmax = ymax;
      fZmin = zmin; fZmax = zmax;
   }
   void GetAll(Double_t& xmin, Double_t& xmax, Double_t& ymin,
               Double_t& ymax, Double_t& zmin, Double_t& zmax)
   {
      xmin = fXmin; xmax = fXmax;
      ymin = fYmin; ymax = fYmax;
      zmin = fZmin; zmax = fZmax;
   }

};

//
// Special timer to dispatch pending events while processing
////////////////////////////////////////////////////////////////////////////////

class TDispatchTimer : public TTimer {
private:
   TProofPlayer    *fPlayer;

public:
   TDispatchTimer(TProofPlayer *p) : TTimer(1000, kFALSE), fPlayer(p) { }

   Bool_t Notify();
};
////////////////////////////////////////////////////////////////////////////////
/// Handle expiration of the timer associated with dispatching pending
/// events while processing. We must act as fast as possible here, so
/// we just set a flag submitting a request for dispatching pending events

Bool_t TDispatchTimer::Notify()
{
   if (gDebug > 0) printf("TDispatchTimer::Notify: called!\n");

   fPlayer->SetBit(TProofPlayer::kDispatchOneEvent);

   // Needed for the next shot
   Reset();
   return kTRUE;
}

//
// Special timer to notify reach of max packet proc time
////////////////////////////////////////////////////////////////////////////////

class TProctimeTimer : public TTimer {
private:
   TProofPlayer    *fPlayer;

public:
   TProctimeTimer(TProofPlayer *p, Long_t to) : TTimer(to, kFALSE), fPlayer(p) { }

   Bool_t Notify();
};
////////////////////////////////////////////////////////////////////////////////
/// Handle expiration of the timer associated with dispatching pending
/// events while processing. We must act as fast as possible here, so
/// we just set a flag submitting a request for dispatching pending events

Bool_t TProctimeTimer::Notify()
{
   if (gDebug > 0) printf("TProctimeTimer::Notify: called!\n");

   fPlayer->SetBit(TProofPlayer::kMaxProcTimeReached);

   // One shot only
   return kTRUE;
}

//
// Special timer to handle stop/abort request via exception raising
////////////////////////////////////////////////////////////////////////////////

class TStopTimer : public TTimer {
private:
   Bool_t           fAbort;
   TProofPlayer    *fPlayer;

public:
   TStopTimer(TProofPlayer *p, Bool_t abort, Int_t to);

   Bool_t Notify();
};

////////////////////////////////////////////////////////////////////////////////
/// Constructor for the timer to stop/abort processing.
/// The 'timeout' is in seconds.
/// Make sure that 'to' make sense, i.e. not larger than 10 days;
/// the minimum value is 10 ms (0 does not seem to start the timer ...).

TStopTimer::TStopTimer(TProofPlayer *p, Bool_t abort, Int_t to)
           : TTimer(((to <= 0 || to > 864000) ? 10 : to * 1000), kFALSE)
{
   if (gDebug > 0)
      Info ("TStopTimer","enter: %d, timeout: %d", abort, to);

   fPlayer = p;
   fAbort = abort;

   if (gDebug > 1)
      Info ("TStopTimer","timeout set to %s ms", fTime.AsString());
}

////////////////////////////////////////////////////////////////////////////////
/// Handle the signal coming from the expiration of the timer
/// associated with an abort or stop request.
/// We raise an exception which will be processed in the
/// event loop.

Bool_t TStopTimer::Notify()
{
   if (gDebug > 0) printf("TStopTimer::Notify: called!\n");

   if (fAbort)
      Throw(kPEX_ABORTED);
   else
      Throw(kPEX_STOPPED);

   return kTRUE;
}

//------------------------------------------------------------------------------

ClassImp(TProofPlayer);

THashList *TProofPlayer::fgDrawInputPars = 0;

////////////////////////////////////////////////////////////////////////////////
/// Default ctor.

TProofPlayer::TProofPlayer(TProof *)
   : fAutoBins(0), fOutput(0), fSelector(0), fCreateSelObj(kTRUE), fSelectorClass(0),
     fFeedbackTimer(0), fFeedbackPeriod(2000),
     fEvIter(0), fSelStatus(0),
     fTotalEvents(0), fReadBytesRun(0), fReadCallsRun(0), fProcessedRun(0),
     fQueryResults(0), fQuery(0), fPreviousQuery(0), fDrawQueries(0),
     fMaxDrawQueries(1), fStopTimer(0), fDispatchTimer(0),
     fProcTimeTimer(0), fProcTime(0),
     fOutputFile(0),
     fSaveMemThreshold(-1), fSavePartialResults(kFALSE), fSaveResultsPerPacket(kFALSE)
{
   fInput         = new TList;
   fExitStatus    = kFinished;
   fProgressStatus = new TProofProgressStatus();
   ResetBit(TProofPlayer::kDispatchOneEvent);
   ResetBit(TProofPlayer::kIsProcessing);
   ResetBit(TProofPlayer::kMaxProcTimeReached);
   ResetBit(TProofPlayer::kMaxProcTimeExtended);

   static Bool_t initLimitsFinder = kFALSE;
   if (!initLimitsFinder && gProofServ && !gProofServ->IsMaster()) {
      THLimitsFinder::SetLimitsFinder(new TProofLimitsFinder);
      initLimitsFinder = kTRUE;
   }

}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TProofPlayer::~TProofPlayer()
{
   fInput->Clear("nodelete");
   SafeDelete(fInput);
   // The output list is owned by fSelector and destroyed in there
   SafeDelete(fSelector);
   SafeDelete(fFeedbackTimer);
   SafeDelete(fEvIter);
   SafeDelete(fQueryResults);
   SafeDelete(fDispatchTimer);
   SafeDelete(fProcTimeTimer);
   SafeDelete(fProcTime);
   SafeDelete(fStopTimer);
}

////////////////////////////////////////////////////////////////////////////////
/// Set processing bit according to 'on'

void TProofPlayer::SetProcessing(Bool_t on)
{
   if (on)
      SetBit(TProofPlayer::kIsProcessing);
   else
      ResetBit(TProofPlayer::kIsProcessing);
}

////////////////////////////////////////////////////////////////////////////////
/// Stop the process after this event. If timeout is positive, start
/// a timer firing after timeout seconds to hard-stop time-expensive
/// events.

void TProofPlayer::StopProcess(Bool_t abort, Int_t timeout)
{
   if (gDebug > 0)
      Info ("StopProcess","abort: %d, timeout: %d", abort, timeout);

   if (fEvIter != 0)
      fEvIter->StopProcess(abort);
   Long_t to = 1;
   if (abort == kTRUE) {
      fExitStatus = kAborted;
   } else {
      fExitStatus = kStopped;
      to = timeout;
   }
   // Start countdown, if needed
   if (to > 0)
      SetStopTimer(kTRUE, abort, to);
}

////////////////////////////////////////////////////////////////////////////////
/// Enable/disable the timer to dispatch pening events while processing.

void TProofPlayer::SetDispatchTimer(Bool_t on)
{
   SafeDelete(fDispatchTimer);
   ResetBit(TProofPlayer::kDispatchOneEvent);
   if (on) {
      fDispatchTimer = new TDispatchTimer(this);
      fDispatchTimer->Start();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Enable/disable the timer to stop/abort processing.
/// The 'timeout' is in seconds.

void TProofPlayer::SetStopTimer(Bool_t on, Bool_t abort, Int_t timeout)
{
   std::lock_guard<std::mutex> lock(fStopTimerMtx);

   // Clean-up the timer
   SafeDelete(fStopTimer);
   if (on) {
      // create timer
      fStopTimer = new TStopTimer(this, abort, timeout);
      // Start the countdown
      fStopTimer->Start();
      if (gDebug > 0)
         Info ("SetStopTimer", "%s timer STARTED (timeout: %d)",
                               (abort ? "ABORT" : "STOP"), timeout);
   } else {
      if (gDebug > 0)
         Info ("SetStopTimer", "timer STOPPED");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Add query result to the list, making sure that there are no
/// duplicates.

void TProofPlayer::AddQueryResult(TQueryResult *q)
{
   if (!q) {
      Warning("AddQueryResult","query undefined - do nothing");
      return;
   }

   // Treat differently normal and draw queries
   if (!(q->IsDraw())) {
      if (!fQueryResults) {
         fQueryResults = new TList;
         fQueryResults->Add(q);
      } else {
         TIter nxr(fQueryResults);
         TQueryResult *qr = 0;
         TQueryResult *qp = 0;
         while ((qr = (TQueryResult *) nxr())) {
            // If same query, remove old version and break
            if (*qr == *q) {
               fQueryResults->Remove(qr);
               delete qr;
               break;
            }
            // Record position according to start time
            if (qr->GetStartTime().Convert() <= q->GetStartTime().Convert())
               qp = qr;
         }

         if (!qp) {
            fQueryResults->AddFirst(q);
         } else {
            fQueryResults->AddAfter(qp, q);
         }
      }
   } else if (IsClient()) {
      // If max reached, eliminate first the oldest one
      if (fDrawQueries == fMaxDrawQueries && fMaxDrawQueries > 0) {
         TIter nxr(fQueryResults);
         TQueryResult *qr = 0;
         while ((qr = (TQueryResult *) nxr())) {
            // If same query, remove old version and break
            if (qr->IsDraw()) {
               fDrawQueries--;
               fQueryResults->Remove(qr);
               delete qr;
               break;
            }
         }
      }
      // Add new draw query
      if (fDrawQueries >= 0 && fDrawQueries < fMaxDrawQueries) {
         fDrawQueries++;
         if (!fQueryResults)
            fQueryResults = new TList;
         fQueryResults->Add(q);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Remove all query result instances referenced 'ref' from
/// the list of results.

void TProofPlayer::RemoveQueryResult(const char *ref)
{
   if (fQueryResults) {
      TIter nxq(fQueryResults);
      TQueryResult *qr = 0;
      while ((qr = (TQueryResult *) nxq())) {
         if (qr->Matches(ref)) {
            fQueryResults->Remove(qr);
            delete qr;
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get query result instances referenced 'ref' from
/// the list of results.

TQueryResult *TProofPlayer::GetQueryResult(const char *ref)
{
   if (fQueryResults) {
      if (ref && strlen(ref) > 0) {
         TIter nxq(fQueryResults);
         TQueryResult *qr = 0;
         while ((qr = (TQueryResult *) nxq())) {
            if (qr->Matches(ref))
               return qr;
         }
      } else {
         // Get last
         return (TQueryResult *) fQueryResults->Last();
      }
   }

   // Nothing found
   return (TQueryResult *)0;
}

////////////////////////////////////////////////////////////////////////////////
/// Set current query and save previous value.

void TProofPlayer::SetCurrentQuery(TQueryResult *q)
{
   fPreviousQuery = fQuery;
   fQuery = q;
}

////////////////////////////////////////////////////////////////////////////////
/// Add object to input list.

void TProofPlayer::AddInput(TObject *inp)
{
   fInput->Add(inp);
}

////////////////////////////////////////////////////////////////////////////////
/// Clear input list.

void TProofPlayer::ClearInput()
{
   fInput->Clear();
}

////////////////////////////////////////////////////////////////////////////////
/// Get output object by name.

TObject *TProofPlayer::GetOutput(const char *name) const
{
   if (fOutput)
      return fOutput->FindObject(name);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Get output list.

TList *TProofPlayer::GetOutputList() const
{
   TList *ol = fOutput;
   if (!ol && fQuery)
      ol = fQuery->GetOutputList();
   return ol;
}

////////////////////////////////////////////////////////////////////////////////
/// Reinitialize fSelector using the selector files in the query result.
/// Needed when Finalize is called after a Process execution for the same
/// selector name.

Int_t TProofPlayer::ReinitSelector(TQueryResult *qr)
{
   Int_t rc = 0;

   // Make sure we have a query
   if (!qr) {
      Info("ReinitSelector", "query undefined - do nothing");
      return -1;
   }

   // Selector name
   TString selec = qr->GetSelecImp()->GetName();
   if (selec.Length() <= 0) {
      Info("ReinitSelector", "selector name undefined - do nothing");
      return -1;
   }

   // Find out if this is a standard selection used for Draw actions
   Bool_t stdselec = TSelector::IsStandardDraw(selec);

   // Find out if this is a precompiled selector: in such a case we do not
   // have the code in TMacros, so we must rely on local libraries
   Bool_t compselec = (selec.Contains(".") || stdselec) ? kFALSE : kTRUE;

   // If not, find out if it needs to be expanded
   TString ipathold;
   if (!stdselec && !compselec) {
      // Check checksums for the versions of the selector files
      Bool_t expandselec = kTRUE;
      TString dir, ipath;
      char *selc = gSystem->Which(TROOT::GetMacroPath(), selec, kReadPermission);
      if (selc) {
         // Check checksums
         TMD5 *md5icur = 0, *md5iold = 0, *md5hcur = 0, *md5hold = 0;
         // Implementation files
         md5icur = TMD5::FileChecksum(selc);
         md5iold = qr->GetSelecImp()->Checksum();
         // Header files
         TString selh(selc);
         Int_t dot = selh.Last('.');
         if (dot != kNPOS) selh.Remove(dot);
         selh += ".h";
         if (!gSystem->AccessPathName(selh, kReadPermission))
            md5hcur = TMD5::FileChecksum(selh);
         md5hold = qr->GetSelecHdr()->Checksum();

         // If nothing has changed nothing to do
         if (md5hcur && md5hold && md5icur && md5iold)
            if (*md5hcur == *md5hold && *md5icur == *md5iold)
               expandselec = kFALSE;

         SafeDelete(md5icur);
         SafeDelete(md5hcur);
         SafeDelete(md5iold);
         SafeDelete(md5hold);
         if (selc) delete [] selc;
      }

      Bool_t ok = kTRUE;
      // Expand selector files, if needed
      if (expandselec) {

         ok = kFALSE;
         // Expand files in a temporary directory
         TUUID u;
         dir = Form("%s/%s",gSystem->TempDirectory(),u.AsString());
         if (!(gSystem->MakeDirectory(dir))) {

            // Export implementation file
            selec = Form("%s/%s",dir.Data(),selec.Data());
            qr->GetSelecImp()->SaveSource(selec);

            // Export header file
            TString seleh = Form("%s/%s",dir.Data(),qr->GetSelecHdr()->GetName());
            qr->GetSelecHdr()->SaveSource(seleh);

            // Adjust include path
            ipathold = gSystem->GetIncludePath();
            ipath = Form("-I%s %s", dir.Data(), gSystem->GetIncludePath());
            gSystem->SetIncludePath(ipath.Data());

            ok = kTRUE;
         }
      }
      TString opt(qr->GetOptions());
      Ssiz_t id = opt.Last('#');
      if (id != kNPOS && id < opt.Length() - 1)
         selec += opt(id + 1, opt.Length());

      if (!ok) {
         Info("ReinitSelector", "problems locating or exporting selector files");
         return -1;
      }
   }

   // Cleanup previous stuff
   SafeDelete(fSelector);
   fSelectorClass = 0;

   // Init the selector now
   Int_t iglevelsave = gErrorIgnoreLevel;
   if (compselec)
      // Silent error printout on first attempt
      gErrorIgnoreLevel = kBreak;

   if ((fSelector = TSelector::GetSelector(selec))) {
      if (compselec)
         gErrorIgnoreLevel = iglevelsave; // restore ignore level
      fSelectorClass = fSelector->IsA();
      fSelector->SetOption(qr->GetOptions());

   } else {
      if (compselec) {
         gErrorIgnoreLevel = iglevelsave; // restore ignore level
         // Retry by loading first the libraries listed in TQueryResult, if any
         if (strlen(qr->GetLibList()) > 0) {
            TString sl(qr->GetLibList());
            TObjArray *oa = sl.Tokenize(" ");
            if (oa) {
               Bool_t retry = kFALSE;
               TIter nxl(oa);
               TObjString *os = 0;
               while ((os = (TObjString *) nxl())) {
                  TString lib = gSystem->BaseName(os->GetName());
                  if (lib != "lib") {
                     lib.ReplaceAll("-l", "lib");
                     if (gSystem->Load(lib) == 0)
                        retry = kTRUE;
                  }
               }
               // Retry now, if the case
               if (retry)
                  fSelector = TSelector::GetSelector(selec);
            }
         }
      }
      if (!fSelector) {
         if (compselec)
            Info("ReinitSelector", "compiled selector re-init failed:"
                                   " automatic reload unsuccessful:"
                                   " please load manually the correct library");
         rc = -1;
      }
   }
   if (fSelector) {
      // Draw needs to reinit temp histos
      fSelector->SetInputList(qr->GetInputList());
      if (stdselec) {
         ((TProofDraw *)fSelector)->DefVar();
      } else {
         // variables may have been initialized in Begin()
         fSelector->Begin(0);
      }
   }

   // Restore original include path, if needed
   if (ipathold.Length() > 0)
      gSystem->SetIncludePath(ipathold.Data());

   return rc;
}

////////////////////////////////////////////////////////////////////////////////
/// Incorporate output object (may not be used in this class).

Int_t TProofPlayer::AddOutputObject(TObject *)
{
   MayNotUse("AddOutputObject");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Incorporate output list (may not be used in this class).

void TProofPlayer::AddOutput(TList *)
{
   MayNotUse("AddOutput");
}

////////////////////////////////////////////////////////////////////////////////
/// Store output list (may not be used in this class).

void TProofPlayer::StoreOutput(TList *)
{
   MayNotUse("StoreOutput");
}

////////////////////////////////////////////////////////////////////////////////
/// Store feedback list (may not be used in this class).

void TProofPlayer::StoreFeedback(TObject *, TList *)
{
   MayNotUse("StoreFeedback");
}

////////////////////////////////////////////////////////////////////////////////
/// Report progress (may not be used in this class).

void TProofPlayer::Progress(Long64_t /*total*/, Long64_t /*processed*/)
{
   MayNotUse("Progress");
}

////////////////////////////////////////////////////////////////////////////////
/// Report progress (may not be used in this class).

void TProofPlayer::Progress(Long64_t /*total*/, Long64_t /*processed*/,
                            Long64_t /*bytesread*/,
                            Float_t /*evtRate*/, Float_t /*mbRate*/,
                            Float_t /*evtrti*/, Float_t /*mbrti*/)
{
   MayNotUse("Progress");
}

////////////////////////////////////////////////////////////////////////////////
/// Report progress (may not be used in this class).

void TProofPlayer::Progress(TProofProgressInfo * /*pi*/)
{
   MayNotUse("Progress");
}

////////////////////////////////////////////////////////////////////////////////
/// Set feedback list (may not be used in this class).

void TProofPlayer::Feedback(TList *)
{
   MayNotUse("Feedback");
}

////////////////////////////////////////////////////////////////////////////////
/// Draw feedback creation proxy. When accessed via TProof avoids
/// link dependency on libProofPlayer.

TDrawFeedback *TProofPlayer::CreateDrawFeedback(TProof *p)
{
   return new TDrawFeedback(p);
}

////////////////////////////////////////////////////////////////////////////////
/// Set draw feedback option.

void TProofPlayer::SetDrawFeedbackOption(TDrawFeedback *f, Option_t *opt)
{
   if (f)
      f->SetOption(opt);
}

////////////////////////////////////////////////////////////////////////////////
/// Delete draw feedback object.

void TProofPlayer::DeleteDrawFeedback(TDrawFeedback *f)
{
   delete f;
}

////////////////////////////////////////////////////////////////////////////////
/// Save the partial results of this query to a dedicated file under the user
/// data directory. The file name has the form
///         <session_tag>.q<query_seq_num>.root
/// The file pat and the file are created if not existing already.
/// Only objects in the outputlist not being TProofOutputFile are saved.
/// The packets list 'packets' is saved if given.
/// Trees not attached to any file are attached to the open file.
/// If 'queryend' is kTRUE evrything is written out (TTrees included).
/// The actual saving action is controlled by 'force' and by fSavePartialResults /
/// fSaveResultsPerPacket:
///
///    fSavePartialResults = kFALSE/kTRUE  no-saving/saving
///    fSaveResultsPerPacket = kFALSE/kTRUE  save-per-query/save-per-packet
///
/// The function CheckMemUsage sets fSavePartialResults = 1 if fSaveMemThreshold > 0 and
/// ProcInfo_t::fMemResident >= fSaveMemThreshold: from that point on partial results
/// are always saved and expensive calls to TSystem::GetProcInfo saved.
/// The switch fSaveResultsPerPacket is instead controlled by the user or admin
/// who can also force saving in all cases; parameter PROOF_SavePartialResults or
/// RC env ProofPlayer.SavePartialResults .
/// However, if 'force' is kTRUE, fSavePartialResults and fSaveResultsPerPacket
/// are ignored.
/// Return -1 in case of problems, 0 otherwise.

Int_t TProofPlayer::SavePartialResults(Bool_t queryend, Bool_t force)
{
   Bool_t save = (force || (fSavePartialResults &&
                 (queryend || fSaveResultsPerPacket))) ? kTRUE : kFALSE;
   if (!save) {
      PDB(kOutput, 2)
         Info("SavePartialResults", "partial result saving disabled");
      return 0;
   }

   // Sanity check
   if (!gProofServ) {
      Error("SavePartialResults", "gProofServ undefined: something really wrong going on!!!");
      return -1;
   }
   if (!fOutput) {
      Error("SavePartialResults", "fOutput undefined: something really wrong going on!!!");
      return -1;
   }

   PDB(kOutput, 1)
      Info("SavePartialResults", "start saving partial results {%d,%d,%d,%d}",
                                 queryend, force, fSavePartialResults, fSaveResultsPerPacket);

   // Get list of processed packets from the iterator
   PDB(kOutput, 2) Info("SavePartialResults", "fEvIter: %p", fEvIter);

   TList *packets = (fEvIter) ? fEvIter->GetPackets() : 0;
   PDB(kOutput, 2) Info("SavePartialResults", "list of packets: %p, sz: %d",
                                              packets, (packets ? packets->GetSize(): -1));

   // Open the file
   const char *oopt = "UPDATE";
   // Check if the file has already been defined
   TString baseName(fOutputFilePath);
   if (fOutputFilePath.IsNull()) {
      baseName.Form("output-%s.q%d.root", gProofServ->GetTopSessionTag(), gProofServ->GetQuerySeqNum());
      if (gProofServ->GetDataDirOpts() && strlen(gProofServ->GetDataDirOpts()) > 0) {
         fOutputFilePath.Form("%s/%s?%s", gProofServ->GetDataDir(), baseName.Data(),
                                          gProofServ->GetDataDirOpts());
      } else {
         fOutputFilePath.Form("%s/%s", gProofServ->GetDataDir(), baseName.Data());
      }
      Info("SavePartialResults", "file with (partial) output: '%s'", fOutputFilePath.Data());
      oopt = "RECREATE";
   }
   // Open the file in write mode
   if (!(fOutputFile = TFile::Open(fOutputFilePath, oopt)) ||
         (fOutputFile && fOutputFile->IsZombie())) {
      Error("SavePartialResults", "cannot open '%s' for writing", fOutputFilePath.Data());
      SafeDelete(fOutputFile);
      return -1;
   }

   // Save current directory
   TDirectory *curdir = gDirectory;
   fOutputFile->cd();

   // Write first the packets list, if required
   if (packets) {
      TDirectory *packetsDir = fOutputFile->mkdir("packets");
      if (packetsDir) packetsDir->cd();
      packets->Write(0, TObject::kSingleKey | TObject::kOverwrite);
      fOutputFile->cd();
   }

   Bool_t notempty = kFALSE;
   // Write out the output list
   TList torm;
   TIter nxo(fOutput);
   TObject *o = 0;
   while ((o = nxo())) {
      // Skip output file drivers
      if (o->InheritsFrom(TProofOutputFile::Class())) continue;
      // Skip control objets
      if (!strncmp(o->GetName(), "PROOF_", 6)) continue;
      // Skip data members mapping
      if (o->InheritsFrom(TOutputListSelectorDataMap::Class())) continue;
      // Skip missing file info
      if (!strcmp(o->GetName(), "MissingFiles")) continue;
      // Trees need a special treatment
      if (o->InheritsFrom("TTree")) {
         TTree *t = (TTree *) o;
         TDirectory *d = t->GetDirectory();
         // If the tree is not attached to any file ...
         if (!d || (d  && !d->InheritsFrom("TFile"))) {
            // ... we attach it
            t->SetDirectory(fOutputFile);
         }
         if (t->GetDirectory() == fOutputFile) {
            if (queryend) {
               // ... we write it out
               o->Write(0, TObject::kOverwrite);
               // At least something in the file
               notempty = kTRUE;
               // Flag for removal from the outputlist
               torm.Add(o);
               // Prevent double-deletion attempts
               t->SetDirectory(0);
            } else {
               // ... or we set in automatic flush mode
               t->SetAutoFlush();
            }
         }
      } else if (queryend || fSaveResultsPerPacket) {
         // Save overwriting what's already there
         o->Write(0, TObject::kOverwrite);
         // At least something in the file
         notempty = kTRUE;
         // Flag for removal from the outputlist
         if (queryend) torm.Add(o);
      }
   }

   // Restore previous directory
   gDirectory = curdir;

   // Close the file if required
   if (notempty) {
      if (!fOutput->FindObject(baseName)) {
         TProofOutputFile *po = 0;
         // Get directions
         TNamed *nm = (TNamed *) fInput->FindObject("PROOF_DefaultOutputOption");
         TString oname = (nm) ? nm->GetTitle() : fOutputFilePath.Data();
         if (nm && oname.BeginsWith("ds:")) {
            oname.Replace(0, 3, "");
            TString qtag =
               TString::Format("%s_q%d", gProofServ->GetTopSessionTag(), gProofServ->GetQuerySeqNum());
            oname.ReplaceAll("<qtag>", qtag);
            // Create the TProofOutputFile for dataset creation
            po = new TProofOutputFile(baseName, "DRO", oname.Data());
         } else {
            Bool_t hasddir = kFALSE;
            // Create the TProofOutputFile for automatic merging
            po = new TProofOutputFile(baseName, "M");
            if (oname.BeginsWith("of:")) oname.Replace(0, 3, "");
            if (gProofServ->IsTopMaster()) {
               if (!strcmp(TUrl(oname, kTRUE).GetProtocol(), "file")) {
                  TString dsrv;
                  TProofServ::GetLocalServer(dsrv);
                  TProofServ::FilterLocalroot(oname, dsrv);
                  oname.Insert(0, dsrv);
               }
            } else {
               if (nm) {
                  // The name has been sent by the client: resolve local place holders
                  oname.ReplaceAll("<file>", baseName);
               } else {
                  // We did not get any indication; the final file will be in the datadir on
                  // the top master and it will be resolved there
                  oname.Form("<datadir>/%s", baseName.Data());
                  hasddir = kTRUE;
               }
            }
            po->SetOutputFileName(oname.Data());
            if (hasddir)
               // Reset the bit, so that <datadir> has a chance to be resolved in AddOutputObject
               po->ResetBit(TProofOutputFile::kOutputFileNameSet);
            po->SetName(gSystem->BaseName(oname.Data()));
         }
         po->AdoptFile(fOutputFile);
         fOutput->Add(po);
         // Flag the nature of this file
         po->SetBit(TProofOutputFile::kSwapFile);
      }
   }
   fOutputFile->Close();
   SafeDelete(fOutputFile);

   // If last call, cleanup the output list from objects saved to file
   if (queryend && torm.GetSize() > 0) {
      TIter nxrm(&torm);
      while ((o = nxrm())) { fOutput->Remove(o); }
   }
   torm.SetOwner(kFALSE);

   PDB(kOutput, 1)
      Info("SavePartialResults", "partial results saved to file");
   // We are done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Make sure that a valid selector object
/// Return -1 in case of problems, 0 otherwise

Int_t TProofPlayer::AssertSelector(const char *selector_file)
{
   if (selector_file && strlen(selector_file)) {
      SafeDelete(fSelector);

      // Get selector files from cache
      TString ocwd = gSystem->WorkingDirectory();
      if (gProofServ) {
         gProofServ->GetCacheLock()->Lock();
         gSystem->ChangeDirectory(gProofServ->GetCacheDir());
      }

      fSelector = TSelector::GetSelector(selector_file);

      if (gProofServ) {
         gSystem->ChangeDirectory(ocwd);
         gProofServ->GetCacheLock()->Unlock();
      }

      if (!fSelector) {
         Error("AssertSelector", "cannot load: %s", selector_file );
        return -1;
      }

      fCreateSelObj = kTRUE;
      Info("AssertSelector", "Processing via filename (%s)", selector_file);
   } else if (!fSelector) {
      Error("AssertSelector", "no TSelector object define : cannot continue!");
      return -1;
   } else {
      Info("AssertSelector", "Processing via TSelector object");
   }
   // Done
   return 0;
}
////////////////////////////////////////////////////////////////////////////////
/// Update fProgressStatus

void TProofPlayer::UpdateProgressInfo()
{
   if (fProgressStatus) {
      fProgressStatus->IncEntries(fProcessedRun);
      fProgressStatus->SetBytesRead(TFile::GetFileBytesRead()-fReadBytesRun);
      fProgressStatus->SetReadCalls(TFile::GetFileReadCalls()-fReadCallsRun);
      fProgressStatus->SetLastUpdate();
      if (gMonitoringWriter)
         gMonitoringWriter->SendProcessingProgress(fProgressStatus->GetEntries(),
                                                   fReadBytesRun, kFALSE);
      fProcessedRun = 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Process specified TDSet on PROOF worker.
/// The return value is -1 in case of error and TSelector::GetStatus()
/// in case of success.

Long64_t TProofPlayer::Process(TDSet *dset, const char *selector_file,
                               Option_t *option, Long64_t nentries,
                               Long64_t first)
{
   PDB(kGlobal,1) Info("Process","Enter");

   fExitStatus = kFinished;
   fOutput = 0;

   TCleanup clean(this);

   fSelectorClass = 0;
   TString wmsg;
   TRY {
      if (AssertSelector(selector_file) != 0 || !fSelector) {
         Error("Process", "cannot assert the selector object");
         return -1;
      }

      fSelectorClass = fSelector->IsA();
      Int_t version = fSelector->Version();
      if (version == 0 && IsClient()) fSelector->GetOutputList()->Clear();

      fOutput = (THashList *) fSelector->GetOutputList();

      if (gProofServ)
         TPerfStats::Start(fInput, fOutput);

      fSelStatus = new TStatus;
      fOutput->Add(fSelStatus);

      fSelector->SetOption(option);
      fSelector->SetInputList(fInput);

      // If in sequential (0-PROOF) mode validate the data set to get
      // the number of entries
      fTotalEvents = nentries;
      if (fTotalEvents < 0 && gProofServ &&
         gProofServ->IsMaster() && !gProofServ->IsParallel()) {
         dset->Validate();
         dset->Reset();
         TDSetElement *e = 0;
         while ((e = dset->Next())) {
            fTotalEvents += e->GetNum();
         }
      }

      dset->Reset();

      // Set parameters controlling the iterator behaviour
      Int_t useTreeCache = 1;
      if (TProof::GetParameter(fInput, "PROOF_UseTreeCache", useTreeCache) == 0) {
         if (useTreeCache > -1 && useTreeCache < 2)
            gEnv->SetValue("ProofPlayer.UseTreeCache", useTreeCache);
      }
      Long64_t cacheSize = -1;
      if (TProof::GetParameter(fInput, "PROOF_CacheSize", cacheSize) == 0) {
         TString sz = TString::Format("%lld", cacheSize);
         gEnv->SetValue("ProofPlayer.CacheSize", sz.Data());
      }
      // Parallel unzipping
      Int_t useParallelUnzip = 0;
      if (TProof::GetParameter(fInput, "PROOF_UseParallelUnzip", useParallelUnzip) == 0) {
         if (useParallelUnzip > -1 && useParallelUnzip < 2)
            gEnv->SetValue("ProofPlayer.UseParallelUnzip", useParallelUnzip);
      }
      // OS file caching (Mac Os X only)
      Int_t dontCacheFiles = 0;
      if (TProof::GetParameter(fInput, "PROOF_DontCacheFiles", dontCacheFiles) == 0) {
         if (dontCacheFiles == 1)
            gEnv->SetValue("ProofPlayer.DontCacheFiles", 1);
      }
      fEvIter = TEventIter::Create(dset, fSelector, first, nentries);

      // Control file object swap
      //     <how>*10 + <force>
      //     <how> =  0       end of run
      //              1       after each packet
      //     <force> = 0      no, swap only if memory threshold is reached
      //               1      swap in all cases, accordingly to <how>
      Int_t opt = 0;
      if (TProof::GetParameter(fInput, "PROOF_SavePartialResults", opt) != 0) {
         opt = gEnv->GetValue("ProofPlayer.SavePartialResults", 0);
      }
      fSaveResultsPerPacket = (opt >= 10) ? kTRUE : kFALSE;
      fSavePartialResults = (opt%10 > 0) ? kTRUE : kFALSE;
      Info("Process", "save partial results? %d  per-packet? %d", fSavePartialResults, fSaveResultsPerPacket);

      // Memory threshold for file object swap
      Float_t memfrac = gEnv->GetValue("ProofPlayer.SaveMemThreshold", -1.);
      if (memfrac > 0.) {
         // The threshold is per core
         SysInfo_t si;
         if (gSystem->GetSysInfo(&si) == 0) {
            fSaveMemThreshold = (Long_t) ((memfrac * si.fPhysRam * 1024.) / si.fCpus);
            Info("Process", "memory threshold for saving objects to file set to %ld kB",
                                 fSaveMemThreshold);
         } else {
            Error("Process", "cannot get SysInfo_t (!)");
         }
      }

      if (version == 0) {
         PDB(kLoop,1) Info("Process","Call Begin(0)");
         fSelector->Begin(0);
      } else {
         if (IsClient()) {
            // on client (for local run)
            PDB(kLoop,1) Info("Process","Call Begin(0)");
            fSelector->Begin(0);
         }
         if (!fSelStatus->TestBit(TStatus::kNotOk)) {
            PDB(kLoop,1) Info("Process","Call SlaveBegin(0)");
            fSelector->SlaveBegin(0);  // Init is called explicitly
                                       // from GetNextEvent()
         }
      }

   } CATCH(excode) {
      ResetBit(TProofPlayer::kIsProcessing);
      Error("Process","exception %d caught", excode);
      gProofServ->GetCacheLock()->Unlock();
      return -1;
   } ENDTRY;

   // Save the results, if needed, closing the file
   if (SavePartialResults(kFALSE) < 0)
      Warning("Process", "problems seetting up file-object swapping");

   // Create feedback lists, if required
   SetupFeedback();

   if (gMonitoringWriter)
      gMonitoringWriter->SendProcessingStatus("STARTED",kTRUE);

   PDB(kLoop,1)
      Info("Process","Looping over Process()");

   // get the byte read counter at the beginning of processing
   fReadBytesRun = TFile::GetFileBytesRead();
   fReadCallsRun = TFile::GetFileReadCalls();
   fProcessedRun = 0;
   // force the first monitoring info
   if (gMonitoringWriter)
      gMonitoringWriter->SendProcessingProgress(0,0,kTRUE);

   // Start asynchronous timer to dispatch pending events
   SetDispatchTimer(kTRUE);

   // Loop over range
   gAbort = kFALSE;
   Long64_t entry;
   fProgressStatus->Reset();
   if (gProofServ) gProofServ->ResetBit(TProofServ::kHighMemory);

   TRY {

      Int_t mrc = -1;
      // Get the frequency for checking memory consumption and logging information
      Long64_t memlogfreq = -1;
      if (((mrc = TProof::GetParameter(fInput, "PROOF_MemLogFreq", memlogfreq))) != 0) memlogfreq = -1;
      Long64_t singleshot = 1;
      Bool_t warnHWMres = kTRUE, warnHWMvir = kTRUE;
      TString lastMsg("(unfortunately no detailed info is available about current packet)");

      // Initial memory footprint
      if (!CheckMemUsage(singleshot, warnHWMres, warnHWMvir, wmsg)) {
         Error("Process", "%s", wmsg.Data());
         wmsg.Insert(0, TString::Format("ERROR:%s, after SlaveBegin(), ", gProofServ ? gProofServ->GetOrdinal() : "gProofServ is nullptr"));
         fSelStatus->Add(wmsg.Data());
         if (gProofServ) {
            gProofServ->SendAsynMessage(wmsg.Data());
            gProofServ->SetBit(TProofServ::kHighMemory);
         }
         fExitStatus = kStopped;
         ResetBit(TProofPlayer::kIsProcessing);
      } else if (!wmsg.IsNull()) {
         Warning("Process", "%s", wmsg.Data());
      }

      TPair *currentElem = 0;
      // The event loop on the worker
      Long64_t fst = -1, num;
      Long_t maxproctime = -1;
      Bool_t newrun = kFALSE;
      while ((fEvIter->GetNextPacket(fst, num) != -1) &&
              !fSelStatus->TestBit(TStatus::kNotOk) &&
              fSelector->GetAbort() == TSelector::kContinue) {
         // This is needed by the inflate infrastructure to calculate
         // sleeping times
         SetBit(TProofPlayer::kIsProcessing);

         // Give the possibility to the selector to access additional info in the
         // incoming packet
         if (dset->Current()) {
            if (!currentElem) {
               currentElem = new TPair(new TObjString("PROOF_CurrentElement"), dset->Current());
               fInput->Add(currentElem);
            } else {
               if (currentElem->Value() != dset->Current()) {
                  currentElem->SetValue(dset->Current());
               } else if (dset->Current()->TestBit(TDSetElement::kNewRun)) {
                  dset->Current()->ResetBit(TDSetElement::kNewRun);
               }
            }
            if (dset->Current()->TestBit(TDSetElement::kNewPacket)) {
               if (dset->TestBit(TDSet::kEmpty)) {
                  lastMsg = "check logs for possible stacktrace - last cycle:";
               } else {
                  TDSetElement *elem = dynamic_cast<TDSetElement *>(currentElem->Value());
                  TString fn = (elem) ? elem->GetFileName() : "<undef>";
                  lastMsg.Form("while processing dset:'%s', file:'%s'"
                              " - check logs for possible stacktrace - last event:", dset->GetName(), fn.Data());
               }
               TProofServ::SetLastMsg(lastMsg);
            }
            // Set the max proc time, if any
            if (dset->Current()->GetMaxProcTime() >= 0.)
               maxproctime = (Long_t) (1000 * dset->Current()->GetMaxProcTime());
            newrun = (dset->Current()->TestBit(TDSetElement::kNewPacket)) ? kTRUE : kFALSE;
         }

         ResetBit(TProofPlayer::kMaxProcTimeReached);
         ResetBit(TProofPlayer::kMaxProcTimeExtended);
         // Setup packet proc time measurement
         if (maxproctime > 0) {
            if (!fProcTimeTimer) fProcTimeTimer = new TProctimeTimer(this, maxproctime);
            fProcTimeTimer->Start(maxproctime, kTRUE); // One shot
            if (!fProcTime) fProcTime = new TStopwatch();
            fProcTime->Reset();                        // Reset counters
         }
         Long64_t refnum = num;
         if (refnum < 0 && maxproctime <= 0) {
            wmsg.Form("neither entries nor max proc time specified:"
                      " risk of infinite loop: processing aborted");
            Error("Process", "%s", wmsg.Data());
            if (gProofServ) {
               wmsg.Insert(0, TString::Format("ERROR:%s, entry:%lld, ",
                                             gProofServ->GetOrdinal(), fProcessedRun));
               gProofServ->SendAsynMessage(wmsg.Data());
            }
            fExitStatus = kAborted;
            ResetBit(TProofPlayer::kIsProcessing);
            break;
         }
         while (refnum < 0 || num--) {

            // Did we use all our time?
            if (TestBit(TProofPlayer::kMaxProcTimeReached)) {
               fProcTime->Stop();
               if (!newrun && !TestBit(TProofPlayer::kMaxProcTimeExtended) && refnum > 0) {
                  // How much are we left with?
                  Float_t xleft = (refnum > num) ? (Float_t) num / (Float_t) (refnum) : 1.;
                  if (xleft < 0.2) {
                     // Give another try, 1.5 times the remaining measured expected time
                     Long_t mpt = (Long_t) (1500 * num / ((Double_t)(refnum - num) / fProcTime->RealTime()));
                     SetBit(TProofPlayer::kMaxProcTimeExtended);
                     fProcTimeTimer->Start(mpt, kTRUE); // One shot
                     ResetBit(TProofPlayer::kMaxProcTimeReached);
                  }
               }
               if (TestBit(TProofPlayer::kMaxProcTimeReached)) {
                  Info("Process", "max proc time reached (%ld msecs): packet processing stopped:\n%s",
                                  maxproctime, lastMsg.Data());

                  break;
               }
            }

            if (!(!fSelStatus->TestBit(TStatus::kNotOk) &&
                   fSelector->GetAbort() == TSelector::kContinue)) break;

            // Get the netry number, taking into account entry or event lists
            entry = fEvIter->GetEntryNumber(fst);
            fst++;

            // Set the last entry
            TProofServ::SetLastEntry(entry);

            if (fSelector->Version() == 0) {
               PDB(kLoop,3)
                  Info("Process","Call ProcessCut(%lld)", entry);
               if (fSelector->ProcessCut(entry)) {
                  PDB(kLoop,3)
                     Info("Process","Call ProcessFill(%lld)", entry);
                  fSelector->ProcessFill(entry);
               }
            } else {
               PDB(kLoop,3)
                  Info("Process","Call Process(%lld)", entry);
               fSelector->Process(entry);
               if (fSelector->GetAbort() == TSelector::kAbortProcess) {
                  ResetBit(TProofPlayer::kIsProcessing);
                  break;
               } else if (fSelector->GetAbort() == TSelector::kAbortFile) {
                  Info("Process", "packet processing aborted following the"
                                  " selector settings:\n%s", lastMsg.Data());
                  fEvIter->InvalidatePacket();
                  fProgressStatus->SetBit(TProofProgressStatus::kFileCorrupted);
               }
            }
            if (!fSelStatus->TestBit(TStatus::kNotOk)) fProcessedRun++;

            // Check the memory footprint, if required
            if (memlogfreq > 0 && (GetEventsProcessed() + fProcessedRun)%memlogfreq == 0) {
               if (!CheckMemUsage(memlogfreq, warnHWMres, warnHWMvir, wmsg)) {
                  Error("Process", "%s", wmsg.Data());
                  if (gProofServ) {
                     wmsg.Insert(0, TString::Format("ERROR:%s, entry:%lld, ",
                                                   gProofServ->GetOrdinal(), entry));
                     gProofServ->SendAsynMessage(wmsg.Data());
                  }
                  fExitStatus = kStopped;
                  ResetBit(TProofPlayer::kIsProcessing);
                  if (gProofServ) gProofServ->SetBit(TProofServ::kHighMemory);
                  break;
               } else {
                  if (!wmsg.IsNull()) {
                     Warning("Process", "%s", wmsg.Data());
                     if (gProofServ) {
                        wmsg.Insert(0, TString::Format("WARNING:%s, entry:%lld, ",
                                                      gProofServ->GetOrdinal(), entry));
                        gProofServ->SendAsynMessage(wmsg.Data());
                     }
                  }
               }
            }
            if (TestBit(TProofPlayer::kDispatchOneEvent)) {
               gSystem->DispatchOneEvent(kTRUE);
               ResetBit(TProofPlayer::kDispatchOneEvent);
            }
            ResetBit(TProofPlayer::kIsProcessing);
            if (fSelStatus->TestBit(TStatus::kNotOk) || gROOT->IsInterrupted()) break;

            // Make sure that the selector abort status is reset
            if (fSelector->GetAbort() == TSelector::kAbortFile)
               fSelector->Abort("status reset", TSelector::kContinue);
         }
      }

   } CATCH(excode) {
      if (excode == kPEX_STOPPED) {
         Info("Process","received stop-process signal");
         fExitStatus = kStopped;
      } else if (excode == kPEX_ABORTED) {
         gAbort = kTRUE;
         Info("Process","received abort-process signal");
         fExitStatus = kAborted;
      } else {
         Error("Process","exception %d caught", excode);
         // Perhaps we need a dedicated status code here ...
         gAbort = kTRUE;
         fExitStatus = kAborted;
      }
      ResetBit(TProofPlayer::kIsProcessing);
   } ENDTRY;

   // Clean-up the envelop for the current element
   TPair *currentElem = 0;
   if ((currentElem = (TPair *) fInput->FindObject("PROOF_CurrentElement"))) {
      if ((currentElem = (TPair *) fInput->Remove(currentElem))) {
         delete currentElem->Key();
         delete currentElem;
      }
   }

   // Final memory footprint
   Long64_t singleshot = 1;
   Bool_t warnHWMres = kTRUE, warnHWMvir = kTRUE;
   Bool_t shrc = CheckMemUsage(singleshot, warnHWMres, warnHWMvir, wmsg);
   if (!wmsg.IsNull()) Warning("Process", "%s (%s)", wmsg.Data(), shrc ? "warn" : "hwm");

   PDB(kGlobal,2)
      Info("Process","%lld events processed", fProgressStatus->GetEntries());

   if (gMonitoringWriter) {
      gMonitoringWriter->SendProcessingProgress(fProgressStatus->GetEntries(),
                                                TFile::GetFileBytesRead()-fReadBytesRun, kFALSE);
      gMonitoringWriter->SendProcessingStatus("DONE");
   }

   // Stop active timers
   SetDispatchTimer(kFALSE);
   if (fStopTimer != 0)
      SetStopTimer(kFALSE, gAbort);
   if (fFeedbackTimer != 0)
      HandleTimer(0);

   StopFeedback();

   // Save the results, if needed, closing the file
   if (SavePartialResults(kTRUE) < 0)
      Warning("Process", "problems saving the results to file");

   SafeDelete(fEvIter);

   // Finalize

   if (fExitStatus != kAborted) {

      TIter nxo(GetOutputList());
      TObject *o = 0;
      while ((o = nxo())) {
         // Special treatment for files
         if (o->IsA() == TProofOutputFile::Class()) {
            TProofOutputFile *of = (TProofOutputFile *)o;
            of->Print();
            of->SetWorkerOrdinal(gProofServ->GetOrdinal());
            const char *dir = of->GetDir();
            if (!dir || (dir && strlen(dir) <= 0)) {
               of->SetDir(gProofServ->GetSessionDir());
            } else if (dir && strlen(dir) > 0) {
               TUrl u(dir);
               if (!strcmp(u.GetHost(), "localhost") || !strcmp(u.GetHost(), "127.0.0.1") ||
                   !strcmp(u.GetHost(), "localhost.localdomain")) {
                  u.SetHost(TUrl(gSystem->HostName()).GetHostFQDN());
                  of->SetDir(u.GetUrl(kTRUE));
               }
               of->Print();
            }
         }
      }

      MapOutputListToDataMembers();

      if (!fSelStatus->TestBit(TStatus::kNotOk)) {
         if (fSelector->Version() == 0) {
            PDB(kLoop,1) Info("Process","Call Terminate()");
            fSelector->Terminate();
         } else {
            PDB(kLoop,1) Info("Process","Call SlaveTerminate()");
            fSelector->SlaveTerminate();
            if (IsClient() && !fSelStatus->TestBit(TStatus::kNotOk)) {
               PDB(kLoop,1) Info("Process","Call Terminate()");
               fSelector->Terminate();
            }
         }
      }

      // Add Selector status in the output list so it can be returned to the client as done
      // by Tree::Process (see ROOT-748). The status from the various workers will be added.
      fOutput->Add(new TParameter<Long64_t>("PROOF_SelectorStatus", (Long64_t) fSelector->GetStatus()));

      if (gProofServ && !gProofServ->IsParallel()) {  // put all the canvases onto the output list
         TIter nxc(gROOT->GetListOfCanvases());
         while (TObject *c = nxc())
            fOutput->Add(c);
      }
   }

   if (gProofServ)
      TPerfStats::Stop();

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Process specified TDSet on PROOF worker with TSelector object
/// The return value is -1 in case of error and TSelector::GetStatus()
/// in case of success.

Long64_t TProofPlayer::Process(TDSet *dset, TSelector *selector,
                               Option_t *option, Long64_t nentries,
                               Long64_t first)
{
   if (!selector) {
      Error("Process", "selector object undefiend!");
      return -1;
   }

   SafeDelete(fSelector);
   fSelector = selector;
   fCreateSelObj = kFALSE;
   return Process(dset, (const char *)0, option, nentries, first);
}

////////////////////////////////////////////////////////////////////////////////
/// Not implemented: meaningful only in the remote player. Returns kFALSE.

Bool_t TProofPlayer::JoinProcess(TList *)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Check the memory usage, if requested.
/// Return kTRUE if OK, kFALSE if above 95% of at least one between virtual or
/// resident limits are depassed.

Bool_t TProofPlayer::CheckMemUsage(Long64_t &mfreq, Bool_t &w80r,
                                   Bool_t &w80v, TString &wmsg)
{
   Long64_t processed = GetEventsProcessed() + fProcessedRun;
   if (mfreq > 0 && processed%mfreq == 0) {
      // Record the memory information
      ProcInfo_t pi;
      if (!gSystem->GetProcInfo(&pi)){
         wmsg = "";
         if (gProofServ)
            Info("CheckMemUsage|Svc", "Memory %ld virtual %ld resident event %lld",
                                      pi.fMemVirtual, pi.fMemResident, processed);
         // Save info in TStatus
         fSelStatus->SetMemValues(pi.fMemVirtual, pi.fMemResident);
         // Apply limit on virtual memory, if any: warn if above 80%, stop if above 95% of max
         if (TProofServ::GetVirtMemMax() > 0) {
            if (pi.fMemVirtual > TProofServ::GetMemStop() * TProofServ::GetVirtMemMax()) {
               wmsg.Form("using more than %d%% of allowed virtual memory (%ld kB)"
                         " - STOP processing", (Int_t) (TProofServ::GetMemStop() * 100), pi.fMemVirtual);
               return kFALSE;
            } else if (pi.fMemVirtual > TProofServ::GetMemHWM() * TProofServ::GetVirtMemMax() && w80v) {
               // Refine monitoring
               mfreq = 1;
               wmsg.Form("using more than %d%% of allowed virtual memory (%ld kB)",
                         (Int_t) (TProofServ::GetMemHWM() * 100), pi.fMemVirtual);
               w80v = kFALSE;
            }
         }
         // Apply limit on resident memory, if any: warn if above 80%, stop if above 95% of max
         if (TProofServ::GetResMemMax() > 0) {
            if (pi.fMemResident > TProofServ::GetMemStop() * TProofServ::GetResMemMax()) {
               wmsg.Form("using more than %d%% of allowed resident memory (%ld kB)"
                         " - STOP processing", (Int_t) (TProofServ::GetMemStop() * 100), pi.fMemResident);
               return kFALSE;
            } else if (pi.fMemResident > TProofServ::GetMemHWM() * TProofServ::GetResMemMax() && w80r) {
               // Refine monitoring
               mfreq = 1;
               if (wmsg.Length() > 0) {
                  wmsg.Form("using more than %d%% of allowed both virtual and resident memory ({%ld,%ld} kB)",
                            (Int_t) (TProofServ::GetMemHWM() * 100), pi.fMemVirtual, pi.fMemResident);
               } else {
                  wmsg.Form("using more than %d%% of allowed resident memory (%ld kB)",
                            (Int_t) (TProofServ::GetMemHWM() * 100), pi.fMemResident);
               }
               w80r = kFALSE;
            }
         }
         // In saving-partial-results mode flag the saving regime when reached to save expensive calls
         // to TSystem::GetProcInfo in SavePartialResults
         if (fSaveMemThreshold > 0 && pi.fMemResident >= fSaveMemThreshold) fSavePartialResults = kTRUE;
      }
   }
   // Done
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Finalize query (may not be used in this class).

Long64_t TProofPlayer::Finalize(Bool_t, Bool_t)
{
   MayNotUse("Finalize");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Finalize query (may not be used in this class).

Long64_t TProofPlayer::Finalize(TQueryResult *)
{
   MayNotUse("Finalize");
   return -1;
}
////////////////////////////////////////////////////////////////////////////////
/// Merge output (may not be used in this class).

void TProofPlayer::MergeOutput(Bool_t)
{
   MayNotUse("MergeOutput");
   return;
}

////////////////////////////////////////////////////////////////////////////////

void TProofPlayer::MapOutputListToDataMembers() const
{
   TOutputListSelectorDataMap* olsdm = new TOutputListSelectorDataMap(fSelector);
   fOutput->Add(olsdm);
}

////////////////////////////////////////////////////////////////////////////////
/// Update automatic binning parameters for given object "name".

void TProofPlayer::UpdateAutoBin(const char *name,
                                 Double_t& xmin, Double_t& xmax,
                                 Double_t& ymin, Double_t& ymax,
                                 Double_t& zmin, Double_t& zmax)
{
   if ( fAutoBins == 0 ) {
      fAutoBins = new THashList;
   }

   TAutoBinVal *val = (TAutoBinVal*) fAutoBins->FindObject(name);

   if ( val == 0 ) {
      //look for info in higher master
      if (gProofServ && !gProofServ->IsTopMaster()) {
         TString key = name;
         TProofLimitsFinder::AutoBinFunc(key,xmin,xmax,ymin,ymax,zmin,zmax);
      }

      val = new TAutoBinVal(name,xmin,xmax,ymin,ymax,zmin,zmax);
      fAutoBins->Add(val);
   } else {
      val->GetAll(xmin,xmax,ymin,ymax,zmin,zmax);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get next packet (may not be used in this class).

TDSetElement *TProofPlayer::GetNextPacket(TSlave *, TMessage *)
{
   MayNotUse("GetNextPacket");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Set up feedback (may not be used in this class).

void TProofPlayer::SetupFeedback()
{
   MayNotUse("SetupFeedback");
}

////////////////////////////////////////////////////////////////////////////////
/// Stop feedback (may not be used in this class).

void TProofPlayer::StopFeedback()
{
   MayNotUse("StopFeedback");
}

////////////////////////////////////////////////////////////////////////////////
/// Draw (may not be used in this class).

Long64_t TProofPlayer::DrawSelect(TDSet * /*set*/, const char * /*varexp*/,
                                  const char * /*selection*/, Option_t * /*option*/,
                                  Long64_t /*nentries*/, Long64_t /*firstentry*/)
{
   MayNotUse("DrawSelect");
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle tree header request.

void TProofPlayer::HandleGetTreeHeader(TMessage *)
{
   MayNotUse("HandleGetTreeHeader|");
}

////////////////////////////////////////////////////////////////////////////////
/// Receive histo from slave.

void TProofPlayer::HandleRecvHisto(TMessage *mess)
{
   TObject *obj = mess->ReadObject(mess->GetClass());
   if (obj->InheritsFrom(TH1::Class())) {
      TH1 *h = (TH1*)obj;
      h->SetDirectory(0);
      TH1 *horg = (TH1*)gDirectory->GetList()->FindObject(h->GetName());
      if (horg)
         horg->Add(h);
      else
         h->SetDirectory(gDirectory);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw the object if it is a canvas.
/// Return 0 in case of success, 1 if it is not a canvas or libProofDraw
/// is not available.

Int_t TProofPlayer::DrawCanvas(TObject *obj)
{
   static Int_t (*gDrawCanvasHook)(TObject *) = 0;

   // Load the library the first time
   if (!gDrawCanvasHook) {
      // Load library needed for graphics ...
      TString drawlib = "libProofDraw";
      char *p = 0;
      if ((p = gSystem->DynamicPathName(drawlib, kTRUE))) {
         delete[] p;
         if (gSystem->Load(drawlib) != -1) {
            // Locate DrawCanvas
            Func_t f = 0;
            if ((f = gSystem->DynFindSymbol(drawlib,"DrawCanvas")))
               gDrawCanvasHook = (Int_t (*)(TObject *))(f);
            else
               Warning("DrawCanvas", "can't find DrawCanvas");
         } else
            Warning("DrawCanvas", "can't load %s", drawlib.Data());
      } else
         Warning("DrawCanvas", "can't locate %s", drawlib.Data());
   }
   if (gDrawCanvasHook && obj)
      return (*gDrawCanvasHook)(obj);
   // No drawing hook or object undefined
   return 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Parse the arguments from var, sel and opt and fill the selector and
/// object name accordingly.
/// Return 0 in case of success, 1 if libProofDraw is not available.

Int_t TProofPlayer::GetDrawArgs(const char *var, const char *sel, Option_t *opt,
                                TString &selector, TString &objname)
{
   static Int_t (*gGetDrawArgsHook)(const char *, const char *, Option_t *,
                                    TString &, TString &) = 0;

   // Load the library the first time
   if (!gGetDrawArgsHook) {
      // Load library needed for graphics ...
      TString drawlib = "libProofDraw";
      char *p = 0;
      if ((p = gSystem->DynamicPathName(drawlib, kTRUE))) {
         delete[] p;
         if (gSystem->Load(drawlib) != -1) {
            // Locate GetDrawArgs
            Func_t f = 0;
            if ((f = gSystem->DynFindSymbol(drawlib,"GetDrawArgs")))
               gGetDrawArgsHook = (Int_t (*)(const char *, const char *, Option_t *,
                                             TString &, TString &))(f);
            else
               Warning("GetDrawArgs", "can't find GetDrawArgs");
         } else
            Warning("GetDrawArgs", "can't load %s", drawlib.Data());
      } else
         Warning("GetDrawArgs", "can't locate %s", drawlib.Data());
   }
   if (gGetDrawArgsHook)
      return (*gGetDrawArgsHook)(var, sel, opt, selector, objname);
   // No parser hook or object undefined
   return 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Create/destroy a named canvas for feedback

void TProofPlayer::FeedBackCanvas(const char *name, Bool_t create)
{
   static void (*gFeedBackCanvasHook)(const char *, Bool_t) = 0;

   // Load the library the first time
   if (!gFeedBackCanvasHook) {
      // Load library needed for graphics ...
      TString drawlib = "libProofDraw";
      char *p = 0;
      if ((p = gSystem->DynamicPathName(drawlib, kTRUE))) {
         delete[] p;
         if (gSystem->Load(drawlib) != -1) {
            // Locate FeedBackCanvas
            Func_t f = 0;
            if ((f = gSystem->DynFindSymbol(drawlib,"FeedBackCanvas")))
               gFeedBackCanvasHook = (void (*)(const char *, Bool_t))(f);
            else
               Warning("FeedBackCanvas", "can't find FeedBackCanvas");
         } else
            Warning("FeedBackCanvas", "can't load %s", drawlib.Data());
      } else
         Warning("FeedBackCanvas", "can't locate %s", drawlib.Data());
   }
   if (gFeedBackCanvasHook) (*gFeedBackCanvasHook)(name, create);
   // No parser hook or object undefined
   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the size in bytes of the cache

Long64_t TProofPlayer::GetCacheSize()
{
   if (fEvIter) return fEvIter->GetCacheSize();
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the number of entries in the learning phase

Int_t TProofPlayer::GetLearnEntries()
{
   if (fEvIter) return fEvIter->GetLearnEntries();
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Switch on/off merge timer

void TProofPlayerRemote::SetMerging(Bool_t on)
{
   if (on) {
      if (!fMergeSTW) fMergeSTW = new TStopwatch();
      PDB(kGlobal,1)
         Info("SetMerging", "ON: mergers: %d", fProof->fMergersCount);
      if (fNumMergers <= 0 && fProof->fMergersCount > 0)
         fNumMergers = fProof->fMergersCount;
   } else if (fMergeSTW) {
      fMergeSTW->Stop();
      Float_t rt = fMergeSTW->RealTime();
      PDB(kGlobal,1)
         Info("SetMerging", "OFF: rt: %f, mergers: %d", rt, fNumMergers);
      if (fQuery) {
         if (!fProof->TestBit(TProof::kIsClient) || fProof->IsLite()) {
            // On the master (or in Lite()) we set the merging time and the numebr of mergers
            fQuery->SetMergeTime(rt);
            fQuery->SetNumMergers(fNumMergers);
         } else {
            // In a standard client we save the transfer-to-client time
            fQuery->SetRecvTime(rt);
         }
         PDB(kGlobal,2) fQuery->Print("F");
      }
   }
}

//------------------------------------------------------------------------------

ClassImp(TProofPlayerLocal);

////////////////////////////////////////////////////////////////////////////////
/// Process the specified TSelector object 'nentries' times.
/// Used to test the PROOF interator mechanism for cycle-driven selectors in a
/// local session.
/// The return value is -1 in case of error and TSelector::GetStatus()
/// in case of success.

Long64_t TProofPlayerLocal::Process(TSelector *selector,
                                    Long64_t nentries, Option_t *option)
{
   if (!selector) {
      Error("Process", "selector object undefiend!");
      return -1;
   }

   TDSetProxy *set = new TDSetProxy("", "", "");
   set->SetBit(TDSet::kEmpty);
   set->SetBit(TDSet::kIsLocal);
   Long64_t rc = Process(set, selector, option, nentries);
   SafeDelete(set);

   // Done
   return rc;
}

////////////////////////////////////////////////////////////////////////////////
/// Process the specified TSelector file 'nentries' times.
/// Used to test the PROOF interator mechanism for cycle-driven selectors in a
/// local session.
/// Process specified TDSet on PROOF worker with TSelector object
/// The return value is -1 in case of error and TSelector::GetStatus()
/// in case of success.

Long64_t TProofPlayerLocal::Process(const char *selector,
                                    Long64_t nentries, Option_t *option)
{
   TDSetProxy *set = new TDSetProxy("", "", "");
   set->SetBit(TDSet::kEmpty);
   set->SetBit(TDSet::kIsLocal);
   Long64_t rc = Process(set, selector, option, nentries);
   SafeDelete(set);

   // Done
   return rc;
}


//------------------------------------------------------------------------------

ClassImp(TProofPlayerRemote);

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TProofPlayerRemote::~TProofPlayerRemote()
{
   SafeDelete(fOutput);      // owns the output list
   SafeDelete(fOutputLists);

   // Objects stored in maps are already deleted when merging the feedback
   SafeDelete(fFeedbackLists);
   SafeDelete(fPacketizer);

   SafeDelete(fProcessMessage);
}

////////////////////////////////////////////////////////////////////////////////
/// Init the packetizer
/// Return 0 on success (fPacketizer is correctly initialized), -1 on failure.

Int_t TProofPlayerRemote::InitPacketizer(TDSet *dset, Long64_t nentries,
                                         Long64_t first, const char *defpackunit,
                                         const char *defpackdata)
{
   SafeDelete(fPacketizer);
   PDB(kGlobal,1) Info("Process","Enter");
   fDSet = dset;
   fExitStatus = kFinished;

   // This is done here to pickup on the fly changes
   Int_t honebyone = 1;
   if (TProof::GetParameter(fInput, "PROOF_MergeTH1OneByOne", honebyone) != 0)
      honebyone = gEnv->GetValue("ProofPlayer.MergeTH1OneByOne", 1);
   fMergeTH1OneByOne = (honebyone == 1) ? kTRUE : kFALSE;

   Bool_t noData = dset->TestBit(TDSet::kEmpty) ? kTRUE : kFALSE;

   TString packetizer;
   TList *listOfMissingFiles = 0;

   TMethodCall callEnv;
   TClass *cl;
   noData = dset->TestBit(TDSet::kEmpty) ? kTRUE : kFALSE;

   if (noData) {

      if (TProof::GetParameter(fInput, "PROOF_Packetizer", packetizer) != 0)
         packetizer = defpackunit;
      else
         Info("InitPacketizer", "using alternate packetizer: %s", packetizer.Data());

      // Get linked to the related class
      cl = TClass::GetClass(packetizer);
      if (cl == 0) {
         Error("InitPacketizer", "class '%s' not found", packetizer.Data());
         fExitStatus = kAborted;
         return -1;
      }

      // Init the constructor
      callEnv.InitWithPrototype(cl, cl->GetName(),"TList*,Long64_t,TList*,TProofProgressStatus*");
      if (!callEnv.IsValid()) {
         Error("InitPacketizer",
               "cannot find correct constructor for '%s'", cl->GetName());
         fExitStatus = kAborted;
         return -1;
      }
      callEnv.ResetParam();
      callEnv.SetParam((Long_t) fProof->GetListOfActiveSlaves());
      callEnv.SetParam((Long64_t) nentries);
      callEnv.SetParam((Long_t) fInput);
      callEnv.SetParam((Long_t) fProgressStatus);

   } else if (dset->TestBit(TDSet::kMultiDSet)) {

      // We have to process many datasets in one go, keeping them separate
      if (fProof->GetRunStatus() != TProof::kRunning) {
         // We have been asked to stop
         Error("InitPacketizer", "received stop/abort request");
         fExitStatus = kAborted;
         return -1;
      }

      // The multi packetizer
      packetizer = "TPacketizerMulti";

      // Get linked to the related class
      cl = TClass::GetClass(packetizer);
      if (cl == 0) {
         Error("InitPacketizer", "class '%s' not found", packetizer.Data());
         fExitStatus = kAborted;
         return -1;
      }

      // Init the constructor
      callEnv.InitWithPrototype(cl, cl->GetName(),"TDSet*,TList*,Long64_t,Long64_t,TList*,TProofProgressStatus*");
      if (!callEnv.IsValid()) {
         Error("InitPacketizer", "cannot find correct constructor for '%s'", cl->GetName());
         fExitStatus = kAborted;
         return -1;
      }
      callEnv.ResetParam();
      callEnv.SetParam((Long_t) dset);
      callEnv.SetParam((Long_t) fProof->GetListOfActiveSlaves());
      callEnv.SetParam((Long64_t) first);
      callEnv.SetParam((Long64_t) nentries);
      callEnv.SetParam((Long_t) fInput);
      callEnv.SetParam((Long_t) fProgressStatus);

      // We are going to test validity during the packetizer initialization
      dset->SetBit(TDSet::kValidityChecked);
      dset->ResetBit(TDSet::kSomeInvalid);

   } else {

      // Lookup - resolve the end-point urls to optmize the distribution.
      // The lookup was previously called in the packetizer's constructor.
      // A list for the missing files may already have been added to the
      // output list; otherwise, if needed it will be created inside
      if ((listOfMissingFiles = (TList *)fInput->FindObject("MissingFiles"))) {
         // Move it to the output list
         fInput->Remove(listOfMissingFiles);
      } else {
         listOfMissingFiles = new TList;
      }
      // Do the lookup; we only skip it if explicitly requested so.
      TString lkopt;
      if (TProof::GetParameter(fInput, "PROOF_LookupOpt", lkopt) != 0 || lkopt != "none")
         dset->Lookup(kTRUE, &listOfMissingFiles);

      if (fProof->GetRunStatus() != TProof::kRunning) {
         // We have been asked to stop
         Error("InitPacketizer", "received stop/abort request");
         fExitStatus = kAborted;
         return -1;
      }

      if (!(dset->GetListOfElements()) ||
          !(dset->GetListOfElements()->GetSize())) {
         if (gProofServ)
            gProofServ->SendAsynMessage("InitPacketizer: No files from the data set were found - Aborting");
         Error("InitPacketizer", "No files from the data set were found - Aborting");
         fExitStatus = kAborted;
         if (listOfMissingFiles) {
            listOfMissingFiles->SetOwner();
            fOutput->Remove(listOfMissingFiles);
            SafeDelete(listOfMissingFiles);
         }
         return -1;
      }

      if (TProof::GetParameter(fInput, "PROOF_Packetizer", packetizer) != 0)
         // Using standard packetizer TAdaptivePacketizer
         packetizer = defpackdata;
      else
         Info("InitPacketizer", "using alternate packetizer: %s", packetizer.Data());

      // Get linked to the related class
      cl = TClass::GetClass(packetizer);
      if (cl == 0) {
         Error("InitPacketizer", "class '%s' not found", packetizer.Data());
         fExitStatus = kAborted;
         return -1;
      }

      // Init the constructor
      callEnv.InitWithPrototype(cl, cl->GetName(),"TDSet*,TList*,Long64_t,Long64_t,TList*,TProofProgressStatus*");
      if (!callEnv.IsValid()) {
         Error("InitPacketizer", "cannot find correct constructor for '%s'", cl->GetName());
         fExitStatus = kAborted;
         return -1;
      }
      callEnv.ResetParam();
      callEnv.SetParam((Long_t) dset);
      callEnv.SetParam((Long_t) fProof->GetListOfActiveSlaves());
      callEnv.SetParam((Long64_t) first);
      callEnv.SetParam((Long64_t) nentries);
      callEnv.SetParam((Long_t) fInput);
      callEnv.SetParam((Long_t) fProgressStatus);

      // We are going to test validity during the packetizer initialization
      dset->SetBit(TDSet::kValidityChecked);
      dset->ResetBit(TDSet::kSomeInvalid);
   }

   // Get an instance of the packetizer
   Long_t ret = 0;
   callEnv.Execute(ret);
   if ((fPacketizer = (TVirtualPacketizer *)ret) == 0) {
      Error("InitPacketizer", "cannot construct '%s'", cl->GetName());
      fExitStatus = kAborted;
      return -1;
   }

   if (!fPacketizer->IsValid()) {
      Error("InitPacketizer",
            "instantiated packetizer object '%s' is invalid", cl->GetName());
      fExitStatus = kAborted;
      SafeDelete(fPacketizer);
      return -1;
   }

   // In multi mode retrieve the list of missing files
   if (!noData && dset->TestBit(TDSet::kMultiDSet)) {
      if ((listOfMissingFiles = (TList *) fInput->FindObject("MissingFiles"))) {
         // Remove it; it will be added to the output list
         fInput->Remove(listOfMissingFiles);
      }
   }

   if (!noData) {
      // Add invalid elements to the list of missing elements
      TDSetElement *elem = 0;
      if (dset->TestBit(TDSet::kSomeInvalid)) {
         TIter nxe(dset->GetListOfElements());
         while ((elem = (TDSetElement *)nxe())) {
            if (!elem->GetValid()) {
               if (!listOfMissingFiles)
                  listOfMissingFiles = new TList;
               listOfMissingFiles->Add(elem->GetFileInfo(dset->GetType()));
               dset->Remove(elem, kFALSE);
            }
         }
         // The invalid elements have been removed
         dset->ResetBit(TDSet::kSomeInvalid);
      }

      // Record the list of missing or invalid elements in the output list
      if (listOfMissingFiles && listOfMissingFiles->GetSize() > 0) {
         TIter missingFiles(listOfMissingFiles);
         TString msg;
         if (gDebug > 0) {
            TFileInfo *fi = 0;
            while ((fi = (TFileInfo *) missingFiles.Next())) {
               if (fi->GetCurrentUrl()) {
                  msg = Form("File not found: %s - skipping!",
                                                fi->GetCurrentUrl()->GetUrl());
               } else {
                  msg = Form("File not found: %s - skipping!", fi->GetName());
               }
               if (gProofServ) gProofServ->SendAsynMessage(msg.Data());
            }
         }
         // Make sure it will be sent back
         if (!GetOutput("MissingFiles")) {
            listOfMissingFiles->SetName("MissingFiles");
            AddOutputObject(listOfMissingFiles);
         }
         TStatus *tmpStatus = (TStatus *)GetOutput("PROOF_Status");
         if (!tmpStatus) AddOutputObject((tmpStatus = new TStatus()));

         // Estimate how much data are missing
         Int_t ngood = dset->GetListOfElements()->GetSize();
         Int_t nbad = listOfMissingFiles->GetSize();
         Double_t xb = Double_t(nbad) / Double_t(ngood + nbad);
         msg = Form(" About %.2f %c of the requested files (%d out of %d) were missing or unusable; details in"
                    " the 'missingFiles' list", xb * 100., '%', nbad, nbad + ngood);
         tmpStatus->Add(msg.Data());
         msg = Form(" +++\n"
                    " +++ About %.2f %c of the requested files (%d out of %d) are missing or unusable; details in"
                    " the 'MissingFiles' list\n"
                    " +++", xb * 100., '%', nbad, nbad + ngood);
         if (gProofServ) gProofServ->SendAsynMessage(msg.Data());
      } else {
         // Cleanup
         SafeDelete(listOfMissingFiles);
      }
   }

   // Done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Process specified TDSet on PROOF.
/// This method is called on client and on the PROOF master.
/// The return value is -1 in case of an error and TSelector::GetStatus() in
/// in case of success.

Long64_t TProofPlayerRemote::Process(TDSet *dset, const char *selector_file,
                                     Option_t *option, Long64_t nentries,
                                     Long64_t first)
{
   PDB(kGlobal,1) Info("Process", "Enter");

   fDSet = dset;
   fExitStatus = kFinished;

   if (!fProgressStatus) {
      Error("Process", "No progress status");
      return -1;
   }
   fProgressStatus->Reset();

   //   delete fOutput;
   if (!fOutput)
      fOutput = new THashList;
   else
      fOutput->Clear();

   SafeDelete(fFeedbackLists);

   if (fProof->IsMaster()){
      TPerfStats::Start(fInput, fOutput);
   } else {
      TPerfStats::Setup(fInput);
   }

   TStopwatch elapsed;

   // Define filename
   TString fn;
   fSelectorFileName = selector_file;

   if (fCreateSelObj) {
      if(!SendSelector(selector_file)) return -1;
      fn = gSystem->BaseName(selector_file);
   } else {
      fn = selector_file;
   }

   TMessage mesg(kPROOF_PROCESS);

   // Parse option
   Bool_t sync = (fProof->GetQueryMode(option) == TProof::kSync);

   TList *inputtmp = 0;  // List of temporary input objects
   TDSet *set = dset;
   if (fProof->IsMaster()) {

      PDB(kPacketizer,1) Info("Process","Create Proxy TDSet");
      set = new TDSetProxy( dset->GetType(), dset->GetObjName(),
                            dset->GetDirectory() );
      if (dset->TestBit(TDSet::kEmpty))
         set->SetBit(TDSet::kEmpty);

      if (InitPacketizer(dset, nentries, first, "TPacketizerUnit", "TPacketizer") != 0) {
         Error("Process", "cannot init the packetizer");
         fExitStatus = kAborted;
         return -1;
      }

      // Reset start, this is now managed by the packetizer
      first = 0;

      // Negative memlogfreq disable checks.
      // If 0 is passed we try to have 100 messages about memory
      // Otherwise we use the frequency passed.
      Int_t mrc = -1;
      Long64_t memlogfreq = -1, mlf;
      if (gSystem->Getenv("PROOF_MEMLOGFREQ")) {
         TString clf(gSystem->Getenv("PROOF_MEMLOGFREQ"));
         if (clf.IsDigit()) { memlogfreq = clf.Atoi(); mrc = 0; }
      }
      if ((mrc = TProof::GetParameter(fProof->GetInputList(), "PROOF_MemLogFreq", mlf)) == 0) memlogfreq = mlf;
      if (memlogfreq == 0) {
         memlogfreq = fPacketizer->GetTotalEntries()/(fProof->GetParallel()*100);
         if (memlogfreq <= 0) memlogfreq = 1;
      }
      if (mrc == 0) fProof->SetParameter("PROOF_MemLogFreq", memlogfreq);


      // Send input data, if any
      TString emsg;
      if (TProof::SendInputData(fQuery, fProof, emsg) != 0)
         Warning("Process", "could not forward input data: %s", emsg.Data());

      // Attach to the transient histogram with the assigned packets, if required
      if (fInput->FindObject("PROOF_StatsHist") != 0) {
         if (!(fProcPackets = (TH1I *) fOutput->FindObject("PROOF_ProcPcktHist"))) {
            Warning("Process", "could not attach to histogram 'PROOF_ProcPcktHist'");
         } else {
            PDB(kLoop,1)
               Info("Process", "attached to histogram 'PROOF_ProcPcktHist' to record"
                               " packets being processed");
         }
      }

   } else {

      // Check whether we have to enforce the use of submergers
      if (gEnv->Lookup("Proof.UseMergers") && !fInput->FindObject("PROOF_UseMergers")) {
         Int_t smg = gEnv->GetValue("Proof.UseMergers",-1);
         if (smg >= 0) {
            fInput->Add(new TParameter<Int_t>("PROOF_UseMergers", smg));
            if (gEnv->Lookup("Proof.MergersByHost")) {
               Int_t mbh = gEnv->GetValue("Proof.MergersByHost",0);
               if (mbh != 0) {
                  // Administrator settings have the priority
                  TObject *o = 0;
                  if ((o = fInput->FindObject("PROOF_MergersByHost"))) { fInput->Remove(o); delete o; }
                  fInput->Add(new TParameter<Int_t>("PROOF_MergersByHost", mbh));
               }
            }
         }
      }

      // For a new query clients should make sure that the temporary
      // output list is empty
      if (fOutputLists) {
         fOutputLists->Delete();
         delete fOutputLists;
         fOutputLists = 0;
      }

      if (!sync) {
         gSystem->RedirectOutput(fProof->fLogFileName);
         Printf(" ");
         Info("Process","starting new query");
      }

      // Define fSelector in Client if processing with filename
      if (fCreateSelObj) {
         SafeDelete(fSelector);
         if (!(fSelector = TSelector::GetSelector(selector_file))) {
            if (!sync)
               gSystem->RedirectOutput(0);
            return -1;
         }
      }

      fSelectorClass = 0;
      fSelectorClass = fSelector->IsA();

      // Add fSelector to inputlist if processing with object
      if (!fCreateSelObj) {
         // In any input list was set into the selector move it to the PROOF
         // input list, because we do not want to stream the selector one
         if (fSelector->GetInputList() && fSelector->GetInputList()->GetSize() > 0) {
            TIter nxi(fSelector->GetInputList());
            TObject *o = 0;
            while ((o = nxi())) {
               if (!fInput->FindObject(o)) {
                  fInput->Add(o);
                  if (!inputtmp) {
                     inputtmp = new TList;
                     inputtmp->SetOwner(kFALSE);
                  }
                  inputtmp->Add(o);
               }
            }
         }
         fInput->Add(fSelector);
      }
      // Set the input list for initialization
      fSelector->SetInputList(fInput);
      fSelector->SetOption(option);
      if (fSelector->GetOutputList()) fSelector->GetOutputList()->Clear();

      PDB(kLoop,1) Info("Process","Call Begin(0)");
      fSelector->Begin(0);

      // Reset the input list to avoid double streaming and related problems (saving
      // the TQueryResult)
      if (!fCreateSelObj) fSelector->SetInputList(0);

      // Send large input data objects, if any
      fProof->SendInputDataFile();

      if (!sync)
         gSystem->RedirectOutput(0);
   }

   TCleanup clean(this);
   SetupFeedback();

   TString opt = option;

   // Old servers need a dedicated streamer
   if (fProof->fProtocol < 13)
      dset->SetWriteV3(kTRUE);

   // Workers will get the entry ranges from the packetizer
   Long64_t num = (gProofServ && gProofServ->IsMaster() && gProofServ->IsParallel()) ? -1 : nentries;
   Long64_t fst = (gProofServ && gProofServ->IsMaster() && gProofServ->IsParallel()) ? -1 : first;

   // Entry- or Event- list ?
   TEntryList *enl = (!fProof->IsMaster()) ? dynamic_cast<TEntryList *>(set->GetEntryList())
                                           : (TEntryList *)0;
   TEventList *evl = (!fProof->IsMaster() && !enl) ? dynamic_cast<TEventList *>(set->GetEntryList())
                                           : (TEventList *)0;
   if (fProof->fProtocol > 14) {
      if (fProcessMessage) delete fProcessMessage;
      fProcessMessage = new TMessage(kPROOF_PROCESS);
      mesg << set << fn << fInput << opt << num << fst << evl << sync << enl;
      (*fProcessMessage) << set << fn << fInput << opt << num << fst << evl << sync << enl;
   } else {
      mesg << set << fn << fInput << opt << num << fst << evl << sync;
      if (enl)
         // Not supported remotely
         Warning("Process","entry lists not supported by the server");
   }

   // Reset the merging progress information
   fProof->ResetMergePrg();

   Int_t nb = fProof->Broadcast(mesg);
   PDB(kGlobal,1) Info("Process", "Broadcast called: %d workers notified", nb);
   if (fProof->IsLite()) fProof->fNotIdle += nb;

   // Reset streamer choice
   if (fProof->fProtocol < 13)
      dset->SetWriteV3(kFALSE);

   // Redirect logs from master to special log frame
   if (IsClient())
      fProof->fRedirLog = kTRUE;

   if (!IsClient()){
      // Signal the start of finalize for the memory log grepping
      Info("Process|Svc", "Start merging Memory information");
   }

   if (!sync) {
      if (IsClient()) {
         // Asynchronous query: just make sure that asynchronous input
         // is enabled and return the prompt
         PDB(kGlobal,1) Info("Process","Asynchronous processing:"
                                       " activating CollectInputFrom");
         fProof->Activate();

         // Receive the acknowledgement and query sequential number
         fProof->Collect();

         return fProof->fSeqNum;

      } else {
         PDB(kGlobal,1) Info("Process","Calling Collect");
         fProof->Collect();

         HandleTimer(0); // force an update of final result
         // This forces a last call to TPacketizer::HandleTimer via the second argument
         // (the first is ignored). This is needed when some events were skipped so that
         // the total number of entries is not the one requested. The packetizer has no
         // way in such a case to understand that processing is finished: it must be told.
         if (fPacketizer) {
            fPacketizer->StopProcess(kFALSE, kTRUE);
            // The progress timer will now stop itself at the next call
            fPacketizer->SetBit(TVirtualPacketizer::kIsDone);
            // Store process info
            elapsed.Stop();
            if (fQuery)
               fQuery->SetProcessInfo(0, 0., fPacketizer->GetBytesRead(),
                                             fPacketizer->GetInitTime(),
                                             elapsed.RealTime());
         }
         StopFeedback();

         return Finalize(kFALSE,sync);
      }
   } else {

      PDB(kGlobal,1) Info("Process","Synchronous processing: calling Collect");
      fProof->Collect();
      if (!(fProof->IsSync())) {
         // The server required to switch to asynchronous mode
         Info("Process", "switching to the asynchronous mode ...");
         return fProof->fSeqNum;
      }

      // Restore prompt logging, for clients (Collect leaves things as they were
      // at the time it was called)
      if (IsClient())
         fProof->fRedirLog = kFALSE;

      if (!IsClient()) {
         // Force an update of final result
         HandleTimer(0);
         // This forces a last call to TPacketizer::HandleTimer via the second argument
         // (the first is ignored). This is needed when some events were skipped so that
         // the total number of entries is not the one requested. The packetizer has no
         // way in such a case to understand that processing is finished: it must be told.
         if (fPacketizer) {
            fPacketizer->StopProcess(kFALSE, kTRUE);
            // The progress timer will now stop itself at the next call
            fPacketizer->SetBit(TVirtualPacketizer::kIsDone);
            // Store process info
            if (fQuery)
               fQuery->SetProcessInfo(0, 0., fPacketizer->GetBytesRead(),
                                             fPacketizer->GetInitTime(),
                                             fPacketizer->GetProcTime());
         }
      } else {
         // Set the input list: maybe required at termination
         if (!fCreateSelObj) fSelector->SetInputList(fInput);
      }
      StopFeedback();

      Long64_t rc = -1;
      if (!IsClient() || GetExitStatus() != TProofPlayer::kAborted)
         rc = Finalize(kFALSE,sync);

      // Remove temporary input objects, if any
      if (inputtmp) {
         TIter nxi(inputtmp);
         TObject *o = 0;
         while ((o = nxi())) fInput->Remove(o);
         SafeDelete(inputtmp);
      }

      // Done
      return rc;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Process specified TDSet on PROOF.
/// This method is called on client and on the PROOF master.
/// The return value is -1 in case of an error and TSelector::GetStatus() in
/// in case of success.

Long64_t TProofPlayerRemote::Process(TDSet *dset, TSelector *selector,
                                     Option_t *option, Long64_t nentries,
                                     Long64_t first)
{
   if (!selector) {
      Error("Process", "selector object undefined");
      return -1;
   }

   // Define fSelector in Client
   if (IsClient() && (selector != fSelector)) {
      SafeDelete(fSelector);
      fSelector = selector;
   }

   fCreateSelObj = kFALSE;
   Long64_t rc = Process(dset, selector->ClassName(), option, nentries, first);
   fCreateSelObj = kTRUE;

   // Done
   return rc;
}

////////////////////////////////////////////////////////////////////////////////
/// Prepares the given list of new workers to join a progressing process.
/// Returns kTRUE on success, kFALSE otherwise.

Bool_t TProofPlayerRemote::JoinProcess(TList *workers)
{
   if (!fProcessMessage || !fProof || !fPacketizer) {
      Error("Process", "Should not happen: fProcessMessage=%p fProof=%p fPacketizer=%p",
         fProcessMessage, fProof, fPacketizer);
      return kFALSE;
   }

   if (!workers || !fProof->IsMaster()) {
      Error("Process", "Invalid call");
      return kFALSE;
   }

   PDB(kGlobal, 1)
      Info("Process", "Preparing %d new worker(s) to process", workers->GetEntries());

   // Sends the file associated to the TSelector, if necessary
   if (fCreateSelObj) {
      PDB(kGlobal, 2)
         Info("Process", "Sending selector file %s", fSelectorFileName.Data());
      if(!SendSelector(fSelectorFileName.Data())) {
         Error("Process", "Problems in sending selector file %s", fSelectorFileName.Data());
         return kFALSE;
      }
   }

   if (fProof->IsLite()) fProof->fNotIdle += workers->GetSize();

   PDB(kGlobal, 2)
      Info("Process", "Adding new workers to the packetizer");
   if (fPacketizer->AddWorkers(workers) == -1) {
      Error("Process", "Cannot add new workers to the packetizer!");
      return kFALSE;  // TODO: make new wrks inactive
   }

   PDB(kGlobal, 2)
      Info("Process", "Broadcasting process message to new workers");
   fProof->Broadcast(*fProcessMessage, workers);

   // Don't call Collect(): we came here from a global Collect() already which
   // will take care of new workers as well

   return kTRUE;

}

////////////////////////////////////////////////////////////////////////////////
/// Merge output in files

Bool_t TProofPlayerRemote::MergeOutputFiles()
{
   PDB(kOutput,1) Info("MergeOutputFiles", "enter: fOutput size: %d", fOutput->GetSize());
   PDB(kOutput,2) fOutput->ls();

   TList *rmList = 0;
   if (fMergeFiles) {
      TIter nxo(fOutput);
      TObject *o = 0;
      TProofOutputFile *pf = 0;
      while ((o = nxo())) {
         if ((pf = dynamic_cast<TProofOutputFile*>(o))) {

            PDB(kOutput,2) pf->Print();

            if (pf->IsMerge()) {

               // Point to the merger
               Bool_t localMerge = (pf->GetTypeOpt() == TProofOutputFile::kLocal) ? kTRUE : kFALSE;
               TFileMerger *filemerger = pf->GetFileMerger(localMerge);
               if (!filemerger) {
                  Error("MergeOutputFiles", "file merger is null in TProofOutputFile! Protocol error?");
                  pf->Print();
                  continue;
               }
               // If only one instance the list in the merger is not yet created: do it now
               if (!pf->IsMerged()) {
                  PDB(kOutput,2) pf->Print();
                  TString fileLoc = TString::Format("%s/%s", pf->GetDir(), pf->GetFileName());
                  filemerger->AddFile(fileLoc);
               }
               // Datadir
               TString ddir, ddopts;
               if (gProofServ) {
                  ddir.Form("%s/", gProofServ->GetDataDir());
                  if (gProofServ->GetDataDirOpts()) ddopts= gProofServ->GetDataDirOpts();
               }
               // Set the output file
               TString outfile(pf->GetOutputFileName());
               if (outfile.Contains("<datadir>/")) {
                  outfile.ReplaceAll("<datadir>/", ddir.Data());
                  if (!ddopts.IsNull())
                     outfile += TString::Format("?%s", ddopts.Data());
                  pf->SetOutputFileName(outfile);
               }
               if ((gProofServ && gProofServ->IsTopMaster()) || (fProof && fProof->IsLite())) {
                  TFile::EFileType ftyp = TFile::kLocal;
                  TString srv;
                  TProofServ::GetLocalServer(srv);
                  TUrl usrv(srv);
                  Bool_t localFile = kFALSE;
                  if (pf->IsRetrieve()) {
                     // This file will be retrieved by the client: we created it in the data dir
                     // and save the file URL on the client in the title
                     if (outfile.BeginsWith("client:")) outfile.Replace(0, 7, "");
                     TString bn = gSystem->BaseName(TUrl(outfile.Data(), kTRUE).GetFile());
                     // The output file path on the master
                     outfile.Form("%s%s", ddir.Data(), bn.Data());
                     // Save the client path in the title if not defined yet
                     if (strlen(pf->GetTitle()) <= 0) pf->SetTitle(bn);
                     // The file is local
                     localFile = kTRUE;
                  } else {
                     // Check if the file is on the master or elsewhere
                     if (outfile.BeginsWith("master:")) outfile.Replace(0, 7, "");
                     // Check locality
                     TUrl uof(outfile.Data(), kTRUE);
                     TString lfn;
                     ftyp = TFile::GetType(uof.GetUrl(), "RECREATE", &lfn);
                     if (ftyp == TFile::kLocal && !srv.IsNull()) {
                        // Check if is a different server
                        if (uof.GetPort() > 0 && usrv.GetPort() > 0 &&
                            usrv.GetPort() != uof.GetPort()) ftyp = TFile::kNet;
                     }
                     // If it is really local set the file name
                     if (ftyp == TFile::kLocal) outfile = lfn;
                     // The file maybe local
                     if (ftyp == TFile::kLocal || ftyp == TFile::kFile) localFile = kTRUE;
                  }
                  // The remote output file name (the one to be used by the client)
                  TString outfilerem(outfile);
                  // For local files we add the local server
                  if (localFile) {
                     // Remove prefix, if any, if included and if Xrootd
                     TProofServ::FilterLocalroot(outfilerem, srv);
                     outfilerem.Insert(0, srv);
                  }
                  // Save the new remote output filename
                  pf->SetOutputFileName(outfilerem);
                  // Align the filename
                  pf->SetFileName(gSystem->BaseName(outfilerem));
               }
               if (!filemerger->OutputFile(outfile)) {
                  Error("MergeOutputFiles", "cannot open the output file");
                  continue;
               }
               // Merge
               PDB(kSubmerger,2) filemerger->PrintFiles("");
               if (!filemerger->Merge()) {
                  Error("MergeOutputFiles", "cannot merge the output files");
                  continue;
               }
               // Remove the files
               TList *fileList = filemerger->GetMergeList();
               if (fileList) {
                  TIter next(fileList);
                  TObjString *url = 0;
                  while((url = (TObjString*)next())) {
                     TUrl u(url->GetName());
                     if (!strcmp(u.GetProtocol(), "file")) {
                        gSystem->Unlink(u.GetFile());
                     } else {
                        gSystem->Unlink(url->GetName());
                     }
                  }
               }
               // Reset the merger
               filemerger->Reset();

            } else {

               // If not yet merged (for example when having only 1 active worker,
               // we need to create the dataset by calling Merge on an effectively empty list
               if (!pf->IsMerged()) {
                  TList dumlist;
                  dumlist.Add(new TNamed("dum", "dum"));
                  dumlist.SetOwner(kTRUE);
                  pf->Merge(&dumlist);
               }
               // Point to the dataset
               TFileCollection *fc = pf->GetFileCollection();
               if (!fc) {
                  Error("MergeOutputFiles", "file collection is null in TProofOutputFile! Protocol error?");
                  pf->Print();
                  continue;
               }
               // Add the collection to the output list for registration and/or to be returned
               // to the client
               fOutput->Add(fc);
               // Do not cleanup at destruction
               pf->ResetFileCollection();
               // Tell the main thread to register this dataset, if needed
               if (pf->IsRegister()) {
                  TString opt;
                  if ((pf->GetTypeOpt() & TProofOutputFile::kOverwrite)) opt += "O";
                  if ((pf->GetTypeOpt() & TProofOutputFile::kVerify)) opt += "V";
                  if (!fOutput->FindObject("PROOFSERV_RegisterDataSet"))
                     fOutput->Add(new TNamed("PROOFSERV_RegisterDataSet", ""));
                  TString tag = TString::Format("DATASET_%s", pf->GetTitle());
                  fOutput->Add(new TNamed(tag, opt));
               }
               // Remove this object from the output list and schedule it for distruction
               fOutput->Remove(pf);
               if (!rmList) rmList = new TList;
               rmList->Add(pf);
               PDB(kOutput,2) fOutput->Print();
            }
         }
      }
   }

   // Remove objects scheduled for removal
   if (rmList && rmList->GetSize() > 0) {
      TIter nxo(rmList);
      TObject *o = 0;
      while((o = nxo())) {
         fOutput->Remove(o);
      }
      rmList->SetOwner(kTRUE);
      delete rmList;
   }

   PDB(kOutput,1) Info("MergeOutputFiles", "done!");

   // Done
   return kTRUE;
}


////////////////////////////////////////////////////////////////////////////////
/// Set the selector's data members:
/// find the mapping of data members to otuput list entries in the output list
/// and apply it.

void TProofPlayerRemote::SetSelectorDataMembersFromOutputList()
{
   TOutputListSelectorDataMap* olsdm
      = TOutputListSelectorDataMap::FindInList(fOutput);
   if (!olsdm) {
      PDB(kOutput,1) Warning("SetSelectorDataMembersFromOutputList",
                             "failed to find map object in output list!");
      return;
   }

   olsdm->SetDataMembers(fSelector);
}

////////////////////////////////////////////////////////////////////////////////

Long64_t TProofPlayerRemote::Finalize(Bool_t force, Bool_t sync)
{
   // Finalize a query.
   // Returns -1 in case of an error, 0 otherwise.

   if (IsClient()) {
      if (fOutputLists == 0) {
         if (force)
            if (fQuery)
               return fProof->Finalize(Form("%s:%s", fQuery->GetTitle(),
                                                     fQuery->GetName()), force);
      } else {
         // Make sure the all objects are in the output list
         PDB(kGlobal,1) Info("Finalize","Calling Merge Output to finalize the output list");
         MergeOutput();
      }
   }

   Long64_t rv = 0;
   if (fProof->IsMaster()) {

      // Fill information for monitoring and stop it
      TStatus *status = (TStatus *) fOutput->FindObject("PROOF_Status");
      if (!status) {
         // The query was aborted: let's add some info in the output list
         status = new TStatus();
         fOutput->Add(status);
         TString emsg = TString::Format("Query aborted after %lld entries", GetEventsProcessed());
         status->Add(emsg);
      }
      status->SetExitStatus((Int_t) GetExitStatus());

      PDB(kOutput,1) Info("Finalize","Calling Merge Output");
      // Some objects (e.g. histos in autobin) may not have been merged yet
      // do it now
      MergeOutput();

      fOutput->SetOwner();

      // Add the active-wrks-vs-proctime info from the packetizer
      if (fPacketizer) {
         TObject *pperf = (TObject *) fPacketizer->GetProgressPerf(kTRUE);
         if (pperf) fOutput->Add(pperf);
         TList *parms = fPacketizer->GetConfigParams(kTRUE);
         if (parms) {
            TIter nxo(parms);
            TObject *o = 0;
            while ((o = nxo())) fOutput->Add(o);
         }

         // If other invalid elements were found during processing, add them to the
         // list of missing elements
         TDSetElement *elem = 0;
         if (fPacketizer->GetFailedPackets()) {
            TString type = (fPacketizer->TestBit(TVirtualPacketizer::kIsTree)) ? "TTree" : "";
            TList *listOfMissingFiles = (TList *) fOutput->FindObject("MissingFiles");
            if (!listOfMissingFiles) {
               listOfMissingFiles = new TList;
               listOfMissingFiles->SetName("MissingFiles");
            }
            TIter nxe(fPacketizer->GetFailedPackets());
            while ((elem = (TDSetElement *)nxe()))
               listOfMissingFiles->Add(elem->GetFileInfo(type));
            if (!fOutput->FindObject(listOfMissingFiles)) fOutput->Add(listOfMissingFiles);
         }
      }

      TPerfStats::Stop();
      // Save memory usage on master
      Long_t vmaxmst, rmaxmst;
      TPerfStats::GetMemValues(vmaxmst, rmaxmst);
      status->SetMemValues(vmaxmst, rmaxmst, kTRUE);

      SafeDelete(fSelector);

   } else {
      if (fExitStatus != kAborted) {

         if (!sync) {
            // Reinit selector (with multi-sessioning we must do this until
            // TSelector::GetSelector() is optimized to i) avoid reloading of an
            // unchanged selector and ii) invalidate existing instances of
            // reloaded selector)
            if (ReinitSelector(fQuery) == -1) {
               Info("Finalize", "problems reinitializing selector \"%s\"",
                    fQuery->GetSelecImp()->GetName());
               return -1;
            }
         }

         if (fPacketizer)
            if (TList *failedPackets = fPacketizer->GetFailedPackets()) {
               fPacketizer->SetFailedPackets(0);
               failedPackets->SetName("FailedPackets");
               AddOutputObject(failedPackets);

               TStatus *status = (TStatus *)GetOutput("PROOF_Status");
               if (!status) AddOutputObject((status = new TStatus()));
               status->Add("Some packets were not processed! Check the the"
                           " 'FailedPackets' list in the output list");
            }

         // Some input parameters may be needed in Terminate
         fSelector->SetInputList(fInput);

         TList *output = fSelector->GetOutputList();
         if (output) {
            TIter next(fOutput);
            while(TObject* obj = next()) {
               if (fProof->IsParallel() || DrawCanvas(obj) == 1)
                  // Either parallel or not a canvas or not able to display it:
                  // just add to the list
                  output->Add(obj);
            }
         } else {
            Warning("Finalize", "undefined output list in the selector! Protocol error?");
         }

         // We need to do this because the output list can be modified in TSelector::Terminate
         // in a way to invalidate existing objects; so we clean the links when still valid and
         // we re-copy back later
         fOutput->SetOwner(kFALSE);
         fOutput->Clear("nodelete");

         // Map output objects to selector members
         SetSelectorDataMembersFromOutputList();

         PDB(kLoop,1) Info("Finalize","Call Terminate()");
         // This is the end of merging
         SetMerging(kFALSE);
         // We measure the merge time
         fProof->fQuerySTW.Reset();
         // Call Terminate now
         fSelector->Terminate();

         rv = fSelector->GetStatus();

         // Copy the output list back and clean the selector's list
         TIter it(output);
         while(TObject* o = it()) {
            fOutput->Add(o);
         }

         // Save the output list in the current query, if any
         if (fQuery) {
            fQuery->SetOutputList(fOutput);
            // Set in finalized state (cannot be done twice)
            fQuery->SetFinalized();
         } else {
            Warning("Finalize","current TQueryResult object is undefined!");
         }

         if (!fCreateSelObj) {
            fInput->Remove(fSelector);
            fOutput->Remove(fSelector);
            if (output) output->Remove(fSelector);
            fSelector = 0;
         }

         // We have transferred copy of the output objects in TQueryResult,
         // so now we can cleanup the selector, making sure that we do not
         // touch the output objects
         if (output) { output->SetOwner(kFALSE); output->Clear("nodelete"); }
         SafeDelete(fSelector);

         // Delete fOutput (not needed anymore, cannot be finalized twice),
         // making sure that the objects saved in TQueryResult are not deleted
         fOutput->SetOwner(kFALSE);
         fOutput->Clear("nodelete");
         SafeDelete(fOutput);

      } else {

         // Cleanup
         fOutput->SetOwner();
         SafeDelete(fSelector);
         if (!fCreateSelObj) fSelector = 0;
      }
   }
   PDB(kGlobal,1) Info("Process","exit");

   if (!IsClient()) {
      Info("Finalize", "finalization on %s finished", gProofServ->GetPrefix());
   }
   fProof->FinalizationDone();

   return rv;
}

////////////////////////////////////////////////////////////////////////////////
/// Finalize the results of a query already processed.

Long64_t TProofPlayerRemote::Finalize(TQueryResult *qr)
{
   PDB(kGlobal,1) Info("Finalize(TQueryResult *)","Enter");

   if (!IsClient()) {
      Info("Finalize(TQueryResult *)",
           "method to be executed only on the clients");
      return -1;
   }

   if (!qr) {
      Info("Finalize(TQueryResult *)", "query undefined");
      return -1;
   }

   if (qr->IsFinalized()) {
      Info("Finalize(TQueryResult *)", "query already finalized");
      return -1;
   }

   // Reset the list
   if (!fOutput)
      fOutput = new THashList;
   else
      fOutput->Clear();

   // Make sure that the temporary output list is empty
   if (fOutputLists) {
      fOutputLists->Delete();
      delete fOutputLists;
      fOutputLists = 0;
   }

   // Re-init the selector
   gSystem->RedirectOutput(fProof->fLogFileName);

   // Import the output list
   TList *tmp = (TList *) qr->GetOutputList();
   if (!tmp) {
      gSystem->RedirectOutput(0);
      Info("Finalize(TQueryResult *)", "outputlist is empty");
      return -1;
   }
   TList *out = fOutput;
   if (fProof->fProtocol < 11)
      out = new TList;
   TIter nxo(tmp);
   TObject *o = 0;
   while ((o = nxo()))
      out->Add(o->Clone());

   // Adopts the list
   if (fProof->fProtocol < 11) {
      out->SetOwner();
      StoreOutput(out);
   }
   gSystem->RedirectOutput(0);

   SetSelectorDataMembersFromOutputList();

   // Finalize it
   SetCurrentQuery(qr);
   Long64_t rc = Finalize();
   RestorePreviousQuery();

   return rc;
}

////////////////////////////////////////////////////////////////////////////////
/// Send the selector file(s) to master or worker nodes.

Bool_t TProofPlayerRemote::SendSelector(const char* selector_file)
{
   // Check input
   if (!selector_file) {
      Info("SendSelector", "Invalid input: selector (file) name undefined");
      return kFALSE;
   }

   if (!strchr(gSystem->BaseName(selector_file), '.')) {
      if (gDebug > 1)
         Info("SendSelector", "selector name '%s' does not contain a '.':"
              " nothing to send, it will be loaded from a library", selector_file);
      return kTRUE;
   }

   // Extract the fine name first
   TString selec = selector_file;
   TString aclicMode;
   TString arguments;
   TString io;
   selec = gSystem->SplitAclicMode(selec, aclicMode, arguments, io);

   // Expand possible envs or '~'
   gSystem->ExpandPathName(selec);

   // Update the macro path
   TString mp(TROOT::GetMacroPath());
   TString np = gSystem->GetDirName(selec);
   if (!np.IsNull()) {
      np += ":";
      if (!mp.BeginsWith(np) && !mp.Contains(":"+np)) {
         Int_t ip = (mp.BeginsWith(".:")) ? 2 : 0;
         mp.Insert(ip, np);
         TROOT::SetMacroPath(mp);
         if (gDebug > 0)
            Info("SendSelector", "macro path set to '%s'", TROOT::GetMacroPath());
      }
   }

   // Header file
   TString header = selec;
   header.Remove(header.Last('.'));
   header += ".h";
   if (gSystem->AccessPathName(header, kReadPermission)) {
      TString h = header;
      header.Remove(header.Last('.'));
      header += ".hh";
      if (gSystem->AccessPathName(header, kReadPermission)) {
         Info("SendSelector",
              "header file not found: tried: %s %s", h.Data(), header.Data());
         return kFALSE;
      }
   }

   // Send files now
   if (fProof->SendFile(selec, (TProof::kBinary | TProof::kForward | TProof::kCp | TProof::kCpBin)) == -1) {
      Info("SendSelector", "problems sending implementation file %s", selec.Data());
      return kFALSE;
   }
   if (fProof->SendFile(header, (TProof::kBinary | TProof::kForward | TProof::kCp)) == -1) {
      Info("SendSelector", "problems sending header file %s", header.Data());
      return kFALSE;
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Merge objects in output the lists.

void TProofPlayerRemote::MergeOutput(Bool_t saveMemValues)
{
   PDB(kOutput,1) Info("MergeOutput","Enter");

   TObject *obj = 0;
   if (fOutputLists) {

      TIter next(fOutputLists);

      TList *list;
      while ( (list = (TList *) next()) ) {

         if (!(obj = fOutput->FindObject(list->GetName()))) {
            obj = list->First();
            list->Remove(obj);
            fOutput->Add(obj);
         }

         if ( list->IsEmpty() ) continue;

         TMethodCall callEnv;
         if (obj->IsA())
            callEnv.InitWithPrototype(obj->IsA(), "Merge", "TCollection*");
         if (callEnv.IsValid()) {
            callEnv.SetParam((Long_t) list);
            callEnv.Execute(obj);
         } else {
            // No Merge interface, return individual objects
            while ( (obj = list->First()) ) {
               fOutput->Add(obj);
               list->Remove(obj);
            }
         }
      }
      SafeDelete(fOutputLists);

   } else {

      PDB(kOutput,1) Info("MergeOutput","fOutputLists empty");
   }

   if (!IsClient() || fProof->IsLite()) {
      // Merge the output files created on workers, if any
      MergeOutputFiles();
   }

   // If there are TProofOutputFile objects we have to make sure that the internal
   // information is consistent for the cases where this object is going to be merged
   // again (e.g. when using submergers or in a multi-master setup). This may not be
   // the case because the first coming in is taken as reference and it has the
   // internal dir and raw dir of the originating worker.
   TString key;
   TNamed *nm = 0;
   TList rmlist;
   TIter nxo(fOutput);
   while ((obj = nxo())) {
      TProofOutputFile *pf = dynamic_cast<TProofOutputFile *>(obj);
      if (pf) {
         if (gProofServ) {
            PDB(kOutput,2) Info("MergeOutput","found TProofOutputFile '%s'", obj->GetName());
            TString dir(pf->GetOutputFileName());
            PDB(kOutput,2) Info("MergeOutput","outputfilename: '%s'", dir.Data());
            // The dir
            if (dir.Last('/') != kNPOS) dir.Remove(dir.Last('/')+1);
            PDB(kOutput,2) Info("MergeOutput","dir: '%s'", dir.Data());
            pf->SetDir(dir);
            // The raw dir; for xrootd based system we include the 'localroot', if any
            TUrl u(dir);
            dir = u.GetFile();
            TString pfx  = gEnv->GetValue("Path.Localroot","");
            if (!pfx.IsNull() &&
               (!strcmp(u.GetProtocol(), "root") || !strcmp(u.GetProtocol(), "xrd")))
               dir.Insert(0, pfx);
            PDB(kOutput,2) Info("MergeOutput","rawdir: '%s'", dir.Data());
            pf->SetDir(dir, kTRUE);
            // The worker ordinal
            pf->SetWorkerOrdinal(gProofServ ? gProofServ->GetOrdinal() : "0");
            // The saved output file name, if any
            key.Form("PROOF_OutputFileName_%s", pf->GetFileName());
            if ((nm = (TNamed *) fOutput->FindObject(key.Data()))) {
               pf->SetOutputFileName(nm->GetTitle());
               rmlist.Add(nm);
            } else if (TestBit(TVirtualProofPlayer::kIsSubmerger)) {
               pf->SetOutputFileName(0);
               pf->ResetBit(TProofOutputFile::kOutputFileNameSet);
            }
            // The filename (order is important to exclude '.merger' from the key)
            dir = pf->GetFileName();
            if (TestBit(TVirtualProofPlayer::kIsSubmerger)) {
               dir += ".merger";
               pf->SetMerged(kFALSE);
            } else {
               if (dir.EndsWith(".merger")) dir.Remove(dir.Last('.'));
            }
            pf->SetFileName(dir);
         } else if (fProof->IsLite()) {
            // The ordinal
            pf->SetWorkerOrdinal("0");
            // The dir
            pf->SetDir(gSystem->GetDirName(pf->GetOutputFileName()));
            // The filename and raw dir
            TUrl u(pf->GetOutputFileName(), kTRUE);
            pf->SetFileName(gSystem->BaseName(u.GetFile()));
            pf->SetDir(gSystem->GetDirName(u.GetFile()), kTRUE);
            // Notify the output path
            Printf("\nOutput file: %s", pf->GetOutputFileName());
         }
      } else {
         PDB(kOutput,2) Info("MergeOutput","output object '%s' is not a TProofOutputFile", obj->GetName());
      }
   }

   // Remove temporary objects from fOutput
   if (rmlist.GetSize() > 0) {
      TIter nxrm(&rmlist);
      while ((obj = nxrm()))
         fOutput->Remove(obj);
      rmlist.SetOwner(kTRUE);
   }

   // If requested (typically in case of submerger to count possible side-effects in that process)
   // save the measured memory usage
   if (saveMemValues) {
      TPerfStats::Stop();
      // Save memory usage on master
      Long_t vmaxmst, rmaxmst;
      TPerfStats::GetMemValues(vmaxmst, rmaxmst);
      TStatus *status = (TStatus *) fOutput->FindObject("PROOF_Status");
      if (status) status->SetMemValues(vmaxmst, rmaxmst, kFALSE);
   }

   PDB(kOutput,1) fOutput->Print();
   PDB(kOutput,1) Info("MergeOutput","leave (%d object(s))", fOutput->GetSize());
}

////////////////////////////////////////////////////////////////////////////////
/// Progress signal.

void TProofPlayerRemote::Progress(Long64_t total, Long64_t processed)
{
   if (IsClient()) {
      fProof->Progress(total, processed);
   } else {
      // Send to the previous tier
      TMessage m(kPROOF_PROGRESS);
      m << total << processed;
      gProofServ->GetSocket()->Send(m);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Progress signal.

void TProofPlayerRemote::Progress(Long64_t total, Long64_t processed,
                                  Long64_t bytesread,
                                  Float_t initTime, Float_t procTime,
                                  Float_t evtrti, Float_t mbrti)
{
   PDB(kGlobal,1)
      Info("Progress","%lld %lld %lld %f %f %f %f", total, processed, bytesread,
                                             initTime, procTime, evtrti, mbrti);

   if (IsClient()) {
      fProof->Progress(total, processed, bytesread, initTime, procTime, evtrti, mbrti);
   } else {
      // Send to the previous tier
      TMessage m(kPROOF_PROGRESS);
      m << total << processed << bytesread << initTime << procTime << evtrti << mbrti;
      gProofServ->GetSocket()->Send(m);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Progress signal.

void TProofPlayerRemote::Progress(TProofProgressInfo *pi)
{
   if (pi) {
      PDB(kGlobal,1)
         Info("Progress","%lld %lld %lld %f %f %f %f %d %f", pi->fTotal, pi->fProcessed, pi->fBytesRead,
                           pi->fInitTime, pi->fProcTime, pi->fEvtRateI, pi->fMBRateI,
                           pi->fActWorkers, pi->fEffSessions);

      if (IsClient()) {
         fProof->Progress(pi->fTotal, pi->fProcessed, pi->fBytesRead,
                           pi->fInitTime, pi->fProcTime,
                           pi->fEvtRateI, pi->fMBRateI,
                           pi->fActWorkers, pi->fTotSessions, pi->fEffSessions);
      } else {
         // Send to the previous tier
         TMessage m(kPROOF_PROGRESS);
         m << pi;
         gProofServ->GetSocket()->Send(m);
      }
   } else {
      Warning("Progress","TProofProgressInfo object undefined!");
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Feedback signal.

void TProofPlayerRemote::Feedback(TList *objs)
{
   fProof->Feedback(objs);
}

////////////////////////////////////////////////////////////////////////////////
/// Stop process after this event.

void TProofPlayerRemote::StopProcess(Bool_t abort, Int_t)
{
   if (fPacketizer != 0)
      fPacketizer->StopProcess(abort, kFALSE);
   if (abort == kTRUE)
      fExitStatus = kAborted;
   else
      fExitStatus = kStopped;
}

////////////////////////////////////////////////////////////////////////////////
/// Incorporate the received object 'obj' into the output list fOutput.
/// The latter is created if not existing.
/// This method short cuts 'StoreOutput + MergeOutput' optimizing the memory
/// consumption.
/// Returns -1 in case of error, 1 if the object has been merged into another
/// one (so that its ownership has not been taken and can be deleted), and 0
/// otherwise.

Int_t TProofPlayerRemote::AddOutputObject(TObject *obj)
{
   PDB(kOutput,1)
      Info("AddOutputObject","Enter: %p (%s)", obj, obj ? obj->ClassName() : "undef");

   // We must something to process
   if (!obj) {
      PDB(kOutput,1) Info("AddOutputObject","Invalid input (obj == 0x0)");
      return -1;
   }

   // Create the output list, if not yet done
   if (!fOutput)
      fOutput = new THashList;

   // Flag about merging
   Bool_t merged = kTRUE;

   // Process event lists first
   TList *elists = dynamic_cast<TList *> (obj);
   if (elists && !strcmp(elists->GetName(), "PROOF_EventListsList")) {

      // Create a global event list, result of merging the event lists
      // coresponding to the various data set elements
      TEventList *evlist = new TEventList("PROOF_EventList");

      // Iterate the list of event list segments
      TIter nxevl(elists);
      TEventList *evl = 0;
      while ((evl = dynamic_cast<TEventList *> (nxevl()))) {

         // Find the file offset (fDSet is the current TDSet instance)
         // locating the element by name
         TIter nxelem(fDSet->GetListOfElements());
         TDSetElement *elem = 0;
         while ((elem = dynamic_cast<TDSetElement *> (nxelem()))) {
            if (!strcmp(elem->GetFileName(), evl->GetName()))
               break;
         }
         if (!elem) {
            Error("AddOutputObject", "Found an event list for %s, but no object with"
                                     " the same name in the TDSet", evl->GetName());
            continue;
         }
         Long64_t offset = elem->GetTDSetOffset();

         // Shift the list by the number of first event in that file
         Long64_t *arr = evl->GetList();
         Int_t num = evl->GetN();
         if (arr && offset > 0)
            for (Int_t i = 0; i < num; i++)
               arr[i] += offset;

         // Add to the global event list
         evlist->Add(evl);
      }

      // Incorporate the resulting global list in fOutput
      SetLastMergingMsg(evlist);
      Incorporate(evlist, fOutput, merged);
      NotifyMemory(evlist);

      // Delete the global list if merged
      if (merged)
         SafeDelete(evlist);

      // The original object has been transformed in something else; we do
      // not have ownership on it
      return 1;
   }

   // Check if we need to merge files
   TProofOutputFile *pf = dynamic_cast<TProofOutputFile*>(obj);
   if (pf) {
      fMergeFiles = kTRUE;
      if (!IsClient() || fProof->IsLite()) {
         if (pf->IsMerge()) {
            Bool_t hasfout = (pf->GetOutputFileName() &&
                              strlen(pf->GetOutputFileName()) > 0 &&
                              pf->TestBit(TProofOutputFile::kOutputFileNameSet)) ? kTRUE : kFALSE;
            Bool_t setfout = (!hasfout || TestBit(TVirtualProofPlayer::kIsSubmerger)) ? kTRUE : kFALSE;
            if (setfout) {

               TString ddir, ddopts;
               if (gProofServ) {
                  ddir.Form("%s/", gProofServ->GetDataDir());
                  if (gProofServ->GetDataDirOpts()) ddopts = gProofServ->GetDataDirOpts();
               }
               // Set the output file
               TString outfile(pf->GetOutputFileName());
               outfile.ReplaceAll("<datadir>/", ddir.Data());
               if (!ddopts.IsNull()) outfile += TString::Format("?%s", ddopts.Data());
               pf->SetOutputFileName(outfile);

               if (gProofServ) {
                  // If submerger, save first the existing filename, if any
                  if (TestBit(TVirtualProofPlayer::kIsSubmerger) && hasfout) {
                     TString key = TString::Format("PROOF_OutputFileName_%s", pf->GetFileName());
                     if (!fOutput->FindObject(key.Data()))
                        fOutput->Add(new TNamed(key.Data(), pf->GetOutputFileName()));
                  }
                  TString of;
                  TProofServ::GetLocalServer(of);
                  if (of.IsNull()) {
                     // Assume an xroot server running on the machine
                     of.Form("root://%s/", gSystem->HostName());
                     if (gSystem->Getenv("XRDPORT")) {
                        TString sp(gSystem->Getenv("XRDPORT"));
                        if (sp.IsDigit())
                           of.Form("root://%s:%s/", gSystem->HostName(), sp.Data());
                     }
                  }
                  TString sessionPath(gProofServ->GetSessionDir());
                  TProofServ::FilterLocalroot(sessionPath, of);
                  of += TString::Format("%s/%s", sessionPath.Data(), pf->GetFileName());
                  if (TestBit(TVirtualProofPlayer::kIsSubmerger)) {
                     if (!of.EndsWith(".merger")) of += ".merger";
                  } else {
                     if (of.EndsWith(".merger")) of.Remove(of.Last('.'));
                  }
                  pf->SetOutputFileName(of);
               }
            }
            // Notify
            PDB(kOutput, 1) pf->Print();
         }
      } else {
         // On clients notify the output path
         Printf("Output file: %s", pf->GetOutputFileName());
      }
   }

   // For other objects we just run the incorporation procedure
   SetLastMergingMsg(obj);
   Incorporate(obj, fOutput, merged);
   NotifyMemory(obj);

   // We are done
   return (merged ? 1 : 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Control output redirection to TProof::fLogFileW

void TProofPlayerRemote::RedirectOutput(Bool_t on)
{
   if (on && fProof && fProof->fLogFileW) {
      TProofServ::SetErrorHandlerFile(fProof->fLogFileW);
      fErrorHandler = SetErrorHandler(TProofServ::ErrorHandler);
   } else if (!on) {
      if (fErrorHandler) {
         TProofServ::SetErrorHandlerFile(0);
         SetErrorHandler(fErrorHandler);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Incorporate the content of the received output list 'out' into the final
/// output list fOutput. The latter is created if not existing.
/// This method short cuts 'StoreOutput + MergeOutput' limiting the memory
/// consumption.

void TProofPlayerRemote::AddOutput(TList *out)
{
   PDB(kOutput,1) Info("AddOutput","Enter");

   // We must something to process
   if (!out) {
      PDB(kOutput,1) Info("AddOutput","Invalid input (out == 0x0)");
      return;
   }

   // Create the output list, if not yet done
   if (!fOutput)
      fOutput = new THashList;

   // Process event lists first
   Bool_t merged = kTRUE;
   TList *elists = dynamic_cast<TList *> (out->FindObject("PROOF_EventListsList"));
   if (elists) {

      // Create a global event list, result of merging the event lists
      // corresponding to the various data set elements
      TEventList *evlist = new TEventList("PROOF_EventList");

      // Iterate the list of event list segments
      TIter nxevl(elists);
      TEventList *evl = 0;
      while ((evl = dynamic_cast<TEventList *> (nxevl()))) {

         // Find the file offset (fDSet is the current TDSet instance)
         // locating the element by name
         TIter nxelem(fDSet->GetListOfElements());
         TDSetElement *elem = 0;
         while ((elem = dynamic_cast<TDSetElement *> (nxelem()))) {
            if (!strcmp(elem->GetFileName(), evl->GetName()))
               break;
         }
         if (!elem) {
            Error("AddOutput", "Found an event list for %s, but no object with"
                               " the same name in the TDSet", evl->GetName());
            continue;
         }
         Long64_t offset = elem->GetTDSetOffset();

         // Shift the list by the number of first event in that file
         Long64_t *arr = evl->GetList();
         Int_t num = evl->GetN();
         if (arr && offset > 0)
            for (Int_t i = 0; i < num; i++)
               arr[i] += offset;

         // Add to the global event list
         evlist->Add(evl);
      }

      // Remove and delete the events lists object to avoid spoiling iteration
      // during next steps
      out->Remove(elists);
      delete elists;

      // Incorporate the resulting global list in fOutput
      SetLastMergingMsg(evlist);
      Incorporate(evlist, fOutput, merged);
      NotifyMemory(evlist);
   }

   // Iterate on the remaining objects in the received list
   TIter nxo(out);
   TObject *obj = 0;
   while ((obj = nxo())) {
      SetLastMergingMsg(obj);
      Incorporate(obj, fOutput, merged);
      // If not merged, drop from the temporary list, as the ownership
      // passes to fOutput
      if (!merged)
         out->Remove(obj);
      NotifyMemory(obj);
   }

   // Done
   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Printout the memory record after merging object 'obj'
/// This record is used by the memory monitor

void TProofPlayerRemote::NotifyMemory(TObject *obj)
{
   if (fProof && (!IsClient() || fProof->IsLite())){
      ProcInfo_t pi;
      if (!gSystem->GetProcInfo(&pi)){
         // For PROOF-Lite we redirect this output to a the open log file so that the
         // memory monitor can pick these messages up
         RedirectOutput(fProof->IsLite());
         Info("NotifyMemory|Svc", "Memory %ld virtual %ld resident after merging object %s",
                                  pi.fMemVirtual, pi.fMemResident, obj->GetName());
         RedirectOutput(0);
      }
      // Record also values for monitoring
      TPerfStats::SetMemValues();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the message to be notified in case of exception

void TProofPlayerRemote::SetLastMergingMsg(TObject *obj)
{
   TString lastMsg = TString::Format("while merging object '%s'", obj->GetName());
   TProofServ::SetLastMsg(lastMsg);
}

////////////////////////////////////////////////////////////////////////////////
/// Incorporate object 'newobj' in the list 'outlist'.
/// The object is merged with an object of the same name already existing in
/// the list, or just added.
/// The boolean merged is set to kFALSE when the object is just added to 'outlist';
/// this happens if the Merge() method does not exist or if a object named as 'obj'
/// is not already in the list. If the obj is not 'merged' than it should not be
/// deleted, unless outlist is not owner of its objects.
/// Return 0 on success, -1 on error.

Int_t TProofPlayerRemote::Incorporate(TObject *newobj, TList *outlist, Bool_t &merged)
{
   merged = kTRUE;

   PDB(kOutput,1)
      Info("Incorporate", "enter: obj: %p (%s), list: %p",
                          newobj, newobj ? newobj->ClassName() : "undef", outlist);

   // The object and list must exist
   if (!newobj || !outlist) {
      Error("Incorporate","Invalid inputs: obj: %p, list: %p", newobj, outlist);
      return -1;
   }

   // Special treatment for histograms in autobin mode
   Bool_t specialH =
      (!fProof || !fProof->TestBit(TProof::kIsClient) || fProof->IsLite()) ? kTRUE : kFALSE;
   if (specialH && newobj->InheritsFrom(TH1::Class())) {
      if (!HandleHistogram(newobj, merged)) {
         if (merged) {
            PDB(kOutput,1) Info("Incorporate", "histogram object '%s' merged", newobj->GetName());
         } else {
            PDB(kOutput,1) Info("Incorporate", "histogram object '%s' added to the"
                                " appropriate list for delayed merging", newobj->GetName());
         }
         return 0;
      }
   }

   // Check if an object with the same name exists already
   TObject *obj = outlist->FindObject(newobj->GetName());

   // If no, add the new object and return
   if (!obj) {
      outlist->Add(newobj);
      merged = kFALSE;
      // Done
      return 0;
   }

   // Locate the Merge(TCollection *) method
   TMethodCall callEnv;
   if (obj->IsA())
      callEnv.InitWithPrototype(obj->IsA(), "Merge", "TCollection*");
   if (callEnv.IsValid()) {
      // Found: put the object in a one-element list
      static TList *xlist = new TList;
      xlist->Add(newobj);
      // Call the method
      callEnv.SetParam((Long_t) xlist);
      callEnv.Execute(obj);
      // Ready for next call
      xlist->Clear();
   } else {
      // Not found: return individual objects
      outlist->Add(newobj);
      merged = kFALSE;
   }

   // Done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Low statistic histograms need a special treatment when using autobin

TObject *TProofPlayerRemote::HandleHistogram(TObject *obj, Bool_t &merged)
{
   TH1 *h = dynamic_cast<TH1 *>(obj);
   if (!h) {
      // Not an histo
      return obj;
   }

   // This is only used if we return (TObject *)0 and there is only one case
   // when we set this to kTRUE
   merged = kFALSE;

   // Does is still needs binning ?
   Bool_t tobebinned = (h->GetBuffer()) ? kTRUE : kFALSE;

   // Number of entries
   Int_t nent = h->GetBufferLength();
   PDB(kOutput,2) Info("HandleHistogram", "h:%s ent:%d, buffer size: %d",
                       h->GetName(), nent, h->GetBufferSize());

   // Attach to the list in the outputlists, if any
   TList *list = 0;
   if (!fOutputLists) {
      PDB(kOutput,2) Info("HandleHistogram", "create fOutputLists");
      fOutputLists = new TList;
      fOutputLists->SetOwner();
   }
   list = (TList *) fOutputLists->FindObject(h->GetName());

   TH1 *href = 0;
   if (tobebinned) {

      // The histogram needs to be projected in a reasonable range: we
      // do this at the end with all the histos, so we need to create
      // a list here
      if (!list) {
         // Create the list
         list = new TList;
         list->SetName(h->GetName());
         list->SetOwner();
         fOutputLists->Add(list);
         // Move in it any previously merged object from the output list
         if (fOutput && (href = (TH1 *) fOutput->FindObject(h->GetName()))) {
            fOutput->Remove(href);
            list->Add(href);
         }
      }
      TIter nxh(list);
      while ((href = (TH1 *) nxh())) {
         if (href->GetBuffer() && href->GetBufferLength() < nent) break;
      }
      if (href) {
         list->AddBefore(href, h);
      } else {
         list->Add(h);
      }
      // Done
      return (TObject *)0;

   } else {

      if (list) {
         TIter nxh(list);
         while ((href = (TH1 *) nxh())) {
            if (href->GetBuffer() || href->GetEntries() < nent) break;
         }
         if (href) {
            list->AddBefore(href, h);
         } else {
            list->Add(h);
         }
         // Done
         return (TObject *)0;

      } else {
         // Check if we can 'Add' the histogram to an existing one; this is more efficient
         // then using Merge
         TH1 *hout = (TH1*) fOutput->FindObject(h->GetName());
         if (hout) {
            // Remove the existing histo from the output list ...
            fOutput->Remove(hout);
            // ... and create either the list to merge in one-go at the end
            // (more efficient than merging one by one) or, if too big, merge
            // these two and start the 'one-by-one' technology
            Int_t hsz = h->GetNbinsX() * h->GetNbinsY() * h->GetNbinsZ();
            if (fMergeTH1OneByOne || (gProofServ && hsz > gProofServ->GetMsgSizeHWM())) {
               list = new TList;
               list->Add(hout);
               h->Merge(list);
               list->SetOwner();
               delete list;
               return h;
            } else {
               list = new TList;
               list->SetName(h->GetName());
               list->SetOwner();
               fOutputLists->Add(list);
               // Add the existing and the incoming histos
               list->Add(hout);
               list->Add(h);
               // Done
               return (TObject *)0;
            }
         } else {
            // This is the first one; add it to the output list
            fOutput->Add(h);
            return (TObject *)0;
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return kTRUE is the histograms 'h0' and 'h1' have the same binning and ranges
/// on the axis (i.e. if they can be just Add-ed for merging).

Bool_t TProofPlayerRemote::HistoSameAxis(TH1 *h0, TH1 *h1)
{
   Bool_t rc = kFALSE;
   if (!h0 || !h1) return rc;

   TAxis *a0 = 0, *a1 = 0;

   // Check X
   a0 = h0->GetXaxis();
   a1 = h1->GetXaxis();
   if (a0->GetNbins() == a1->GetNbins())
      if (TMath::Abs(a0->GetXmax() - a1->GetXmax()) < 1.e-9)
         if (TMath::Abs(a0->GetXmin() - a1->GetXmin()) < 1.e-9) rc = kTRUE;

   // Check Y, if needed
   if (h0->GetDimension() > 1) {
      rc = kFALSE;
      a0 = h0->GetYaxis();
      a1 = h1->GetYaxis();
      if (a0->GetNbins() == a1->GetNbins())
         if (TMath::Abs(a0->GetXmax() - a1->GetXmax()) < 1.e-9)
            if (TMath::Abs(a0->GetXmin() - a1->GetXmin()) < 1.e-9) rc = kTRUE;
   }

   // Check Z, if needed
   if (h0->GetDimension() > 2) {
      rc = kFALSE;
      a0 = h0->GetZaxis();
      a1 = h1->GetZaxis();
      if (a0->GetNbins() == a1->GetNbins())
         if (TMath::Abs(a0->GetXmax() - a1->GetXmax()) < 1.e-9)
            if (TMath::Abs(a0->GetXmin() - a1->GetXmin()) < 1.e-9) rc = kTRUE;
   }

   // Done
   return rc;
}

////////////////////////////////////////////////////////////////////////////////
/// Store received output list.

void TProofPlayerRemote::StoreOutput(TList *out)
{
   PDB(kOutput,1) Info("StoreOutput","Enter");

   if ( out == 0 ) {
      PDB(kOutput,1) Info("StoreOutput","Leave (empty)");
      return;
   }

   TIter next(out);
   out->SetOwner(kFALSE);  // take ownership of the contents

   if (fOutputLists == 0) {
      PDB(kOutput,2) Info("StoreOutput","Create fOutputLists");
      fOutputLists = new TList;
      fOutputLists->SetOwner();
   }
   // process eventlists first
   TList* lists = dynamic_cast<TList*> (out->FindObject("PROOF_EventListsList"));
   if (lists) {
      out->Remove(lists);
      TEventList *mainList = new TEventList("PROOF_EventList");
      out->Add(mainList);
      TIter it(lists);
      TEventList *aList;
      while ( (aList = dynamic_cast<TEventList*> (it())) ) {
         // find file offset
         TIter nxe(fDSet->GetListOfElements());
         TDSetElement *elem;
         while ( (elem = dynamic_cast<TDSetElement*> (nxe())) ) {
            if (strcmp(elem->GetFileName(), aList->GetName()) == 0)
               break;
         }
         if (!elem) {
            Error("StoreOutput", "found the EventList for %s, but no object with that name "
                                 "in the TDSet", aList->GetName());
            continue;
         }
         Long64_t offset = elem->GetTDSetOffset();

         // shift the list by the number of first event in that file
         Long64_t *arr = aList->GetList();
         Int_t num = aList->GetN();
         if (arr && offset)
            for (int i = 0; i < num; i++)
               arr[i] += offset;

         mainList->Add(aList);           // add to the main list
      }
      delete lists;
   }

   TObject *obj;
   while( (obj = next()) ) {
      PDB(kOutput,2) Info("StoreOutput","find list for '%s'", obj->GetName() );

      TList *list = (TList *) fOutputLists->FindObject( obj->GetName() );
      if ( list == 0 ) {
         PDB(kOutput,2) Info("StoreOutput", "list for '%s' not found (creating)", obj->GetName());
         list = new TList;
         list->SetName( obj->GetName() );
         list->SetOwner();
         fOutputLists->Add( list );
      }
      list->Add( obj );
   }

   delete out;
   PDB(kOutput,1) Info("StoreOutput", "leave");
}

////////////////////////////////////////////////////////////////////////////////
/// Merge feedback lists.

TList *TProofPlayerRemote::MergeFeedback()
{
   PDB(kFeedback,1)
      Info("MergeFeedback","Enter");

   if ( fFeedbackLists == 0 ) {
      PDB(kFeedback,1)
         Info("MergeFeedback","Leave (no output)");
      return 0;
   }

   TList *fb = new TList;   // collection of feedback objects
   fb->SetOwner();

   TIter next(fFeedbackLists);

   TMap *map;
   while ( (map = (TMap*) next()) ) {

      PDB(kFeedback,2)
         Info("MergeFeedback", "map %s size: %d", map->GetName(), map->GetSize());

      // turn map into list ...

      TList *list = new TList;
      TIter keys(map);

#ifndef R__TH1MERGEFIXED
      Int_t nbmx = -1;
      TObject *oref = 0;
#endif
      while ( TObject *key = keys() ) {
         TObject *o = map->GetValue(key);
         TH1 *h = dynamic_cast<TH1 *>(o);
#ifndef R__TH1MERGEFIXED
         // Temporary fix for to cope with the problem in TH1::Merge.
         // We need to use a reference histo the one with the largest number
         // of bins so that the histos from all submasters can be correctly
         // fit in
         if (h && !strncmp(o->GetName(),"PROOF_",6)) {
            if (h->GetNbinsX() > nbmx) {
               nbmx=  h->GetNbinsX();
               oref = o;
            }
         }
#endif
         if (h) {
            TIter nxh(list);
            TH1 *href= 0;
            while ((href = (TH1 *)nxh())) {
               if (h->GetBuffer()) {
                  if (href->GetBuffer() && href->GetBufferLength() < h->GetBufferLength()) break;
               } else {
                  if (href->GetBuffer() || href->GetEntries() < h->GetEntries()) break;
               }
            }
            if (href) {
               list->AddBefore(href, h);
            } else {
               list->Add(h);
            }
         } else {
            list->Add(o);
         }
      }

      // clone first object, remove from list
#ifdef R__TH1MERGEFIXED
      TObject *obj = list->First();
#else
      TObject *obj = (oref) ? oref : list->First();
#endif
      list->Remove(obj);
      obj = obj->Clone();
      fb->Add(obj);

      if ( list->IsEmpty() ) {
         delete list;
         continue;
      }

      // merge list with clone
      TMethodCall callEnv;
      if (obj->IsA())
         callEnv.InitWithPrototype(obj->IsA(), "Merge", "TCollection*");
      if (callEnv.IsValid()) {
         callEnv.SetParam((Long_t) list);
         callEnv.Execute(obj);
      } else {
         // No Merge interface, return copy of individual objects
         while ( (obj = list->First()) ) {
            fb->Add(obj->Clone());
            list->Remove(obj);
         }
      }

      delete list;
   }

   PDB(kFeedback,1)
      Info("MergeFeedback","Leave (%d object(s))", fb->GetSize());

   return fb;
}

////////////////////////////////////////////////////////////////////////////////
/// Store feedback results from the specified slave.

void TProofPlayerRemote::StoreFeedback(TObject *slave, TList *out)
{
   PDB(kFeedback,1)
      Info("StoreFeedback","Enter");

   if ( out == 0 ) {
      PDB(kFeedback,1)
         Info("StoreFeedback","Leave (empty)");
      return;
   }

   if ( IsClient() ) {
      // in client
      Feedback(out);
      delete out;
      return;
   }

   if (fFeedbackLists == 0) {
      PDB(kFeedback,2) Info("StoreFeedback","Create fFeedbackLists");
      fFeedbackLists = new TList;
      fFeedbackLists->SetOwner();
   }

   TIter next(out);
   out->SetOwner(kFALSE);  // take ownership of the contents

   const char *ord = ((TSlave*) slave)->GetOrdinal();

   TObject *obj;
   while( (obj = next()) ) {
      PDB(kFeedback,2)
         Info("StoreFeedback","%s: Find '%s'", ord, obj->GetName() );
      TMap *map = (TMap*) fFeedbackLists->FindObject(obj->GetName());
      if ( map == 0 ) {
         PDB(kFeedback,2)
            Info("StoreFeedback", "%s: map for '%s' not found (creating)", ord, obj->GetName());
         // Map must not be owner (ownership is with regards to the keys (only))
         map = new TMap;
         map->SetName(obj->GetName());
         fFeedbackLists->Add(map);
      } else {
         PDB(kFeedback,2)
            Info("StoreFeedback","%s: removing previous value", ord);
         if (map->GetValue(slave))
            delete map->GetValue(slave);
         map->Remove(slave);
      }
      map->Add(slave, obj);
      PDB(kFeedback,2)
         Info("StoreFeedback","%s: %s, size: %d", ord, obj->GetName(), map->GetSize());
   }

   delete out;
   PDB(kFeedback,1)
      Info("StoreFeedback","Leave");
}

////////////////////////////////////////////////////////////////////////////////
/// Setup reporting of feedback objects.

void TProofPlayerRemote::SetupFeedback()
{
   if (IsClient()) return; // Client does not need timer

   fFeedback = (TList*) fInput->FindObject("FeedbackList");

   PDB(kFeedback,1) Info("SetupFeedback","\"FeedbackList\" %sfound",
      fFeedback == 0 ? "NOT ":"");

   if (fFeedback == 0 || fFeedback->GetSize() == 0) return;

   // OK, feedback was requested, setup the timer
   SafeDelete(fFeedbackTimer);
   fFeedbackPeriod = 2000;
   TProof::GetParameter(fInput, "PROOF_FeedbackPeriod", fFeedbackPeriod);
   fFeedbackTimer = new TTimer;
   fFeedbackTimer->SetObject(this);
   fFeedbackTimer->Start(fFeedbackPeriod, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Stop reporting of feedback objects.

void TProofPlayerRemote::StopFeedback()
{
   if (fFeedbackTimer == 0) return;

   PDB(kFeedback,1) Info("StopFeedback","Stop Timer");

   SafeDelete(fFeedbackTimer);
}

////////////////////////////////////////////////////////////////////////////////
/// Send feedback objects to client.

Bool_t TProofPlayerRemote::HandleTimer(TTimer *)
{
   PDB(kFeedback,2) Info("HandleTimer","Entry");

   if (fFeedbackTimer == 0) return kFALSE; // timer already switched off

   // process local feedback objects

   TList *fb = new TList;
   fb->SetOwner();

   TIter next(fFeedback);
   while( TObjString *name = (TObjString*) next() ) {
      TObject *o = fOutput->FindObject(name->GetName());
      if (o != 0) {
         fb->Add(o->Clone());
         // remove the corresponding entry from the feedback list
         TMap *m = 0;
         if (fFeedbackLists &&
            (m = (TMap *) fFeedbackLists->FindObject(name->GetName()))) {
            fFeedbackLists->Remove(m);
            m->DeleteValues();
            delete m;
         }
      }
   }

   if (fb->GetSize() > 0) {
      StoreFeedback(this, fb); // adopts fb
   } else {
      delete fb;
   }

   if (fFeedbackLists == 0) {
      fFeedbackTimer->Start(fFeedbackPeriod, kTRUE);   // maybe next time
      return kFALSE;
   }

   fb = MergeFeedback();

   PDB(kFeedback,2) Info("HandleTimer","Sending %d objects", fb->GetSize());

   TMessage m(kPROOF_FEEDBACK);
   m << fb;

   // send message to client;
   gProofServ->GetSocket()->Send(m);

   delete fb;

   fFeedbackTimer->Start(fFeedbackPeriod, kTRUE);

   return kFALSE; // ignored?
}

////////////////////////////////////////////////////////////////////////////////
/// Get next packet for specified slave.

TDSetElement *TProofPlayerRemote::GetNextPacket(TSlave *slave, TMessage *r)
{
   // The first call to this determines the end of initialization
   SetInitTime();

   if (fProcPackets) {
      Int_t bin = fProcPackets->GetXaxis()->FindBin(slave->GetOrdinal());
      if (bin >= 0) {
         if (fProcPackets->GetBinContent(bin) > 0)
            fProcPackets->Fill(slave->GetOrdinal(), -1);
      }
   }

   TDSetElement *e = fPacketizer->GetNextPacket( slave, r );

   if (e == 0) {
      PDB(kPacketizer,2)
         Info("GetNextPacket","%s: done!", slave->GetOrdinal());
   } else if (e == (TDSetElement*) -1) {
      PDB(kPacketizer,2)
         Info("GetNextPacket","%s: waiting ...", slave->GetOrdinal());
   } else {
      PDB(kPacketizer,2)
         Info("GetNextPacket","%s (%s): '%s' '%s' '%s' %lld %lld",
              slave->GetOrdinal(), slave->GetName(), e->GetFileName(),
              e->GetDirectory(), e->GetObjName(), e->GetFirst(), e->GetNum());
      if (fProcPackets) fProcPackets->Fill(slave->GetOrdinal(), 1);
   }

   return e;
}

////////////////////////////////////////////////////////////////////////////////
/// Is the player running on the client?

Bool_t TProofPlayerRemote::IsClient() const
{
   return fProof ? fProof->TestBit(TProof::kIsClient) : kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw (support for TChain::Draw()).
/// Returns -1 in case of error or number of selected events in case of success.

Long64_t TProofPlayerRemote::DrawSelect(TDSet *set, const char *varexp,
                                        const char *selection, Option_t *option,
                                        Long64_t nentries, Long64_t firstentry)
{
   if (!fgDrawInputPars) {
      fgDrawInputPars = new THashList;
      fgDrawInputPars->Add(new TObjString("FeedbackList"));
      fgDrawInputPars->Add(new TObjString("PROOF_ChainWeight"));
      fgDrawInputPars->Add(new TObjString("PROOF_LineColor"));
      fgDrawInputPars->Add(new TObjString("PROOF_LineStyle"));
      fgDrawInputPars->Add(new TObjString("PROOF_LineWidth"));
      fgDrawInputPars->Add(new TObjString("PROOF_MarkerColor"));
      fgDrawInputPars->Add(new TObjString("PROOF_MarkerStyle"));
      fgDrawInputPars->Add(new TObjString("PROOF_MarkerSize"));
      fgDrawInputPars->Add(new TObjString("PROOF_FillColor"));
      fgDrawInputPars->Add(new TObjString("PROOF_FillStyle"));
      fgDrawInputPars->Add(new TObjString("PROOF_ListOfAliases"));
   }

   TString selector, objname;
   if (GetDrawArgs(varexp, selection, option, selector, objname) != 0) {
      Error("DrawSelect", "parsing arguments");
      return -1;
   }

   TNamed *varexpobj = new TNamed("varexp", varexp);
   TNamed *selectionobj = new TNamed("selection", selection);

   // Save the current input list
   TObject *o = 0;
   TList *savedInput = new TList;
   TIter nxi(fInput);
   while ((o = nxi())) {
      savedInput->Add(o);
      TString n(o->GetName());
      if (fgDrawInputPars &&
          !fgDrawInputPars->FindObject(o->GetName()) &&
          !n.BeginsWith("alias:")) fInput->Remove(o);
   }

   fInput->Add(varexpobj);
   fInput->Add(selectionobj);

   // Make sure we have an object name
   if (objname == "") objname = "htemp";

   fProof->AddFeedback(objname);
   Long64_t r = Process(set, selector, option, nentries, firstentry);
   fProof->RemoveFeedback(objname);

   fInput->Remove(varexpobj);
   fInput->Remove(selectionobj);
   if (TNamed *opt = dynamic_cast<TNamed*> (fInput->FindObject("PROOF_OPTIONS"))) {
      fInput->Remove(opt);
      delete opt;
   }

   delete varexpobj;
   delete selectionobj;

   // Restore the input list
   fInput->Clear();
   TIter nxsi(savedInput);
   while ((o = nxsi()))
      fInput->Add(o);
   savedInput->SetOwner(kFALSE);
   delete savedInput;

   return r;
}

////////////////////////////////////////////////////////////////////////////////
/// Set init time

void TProofPlayerRemote::SetInitTime()
{
   if (fPacketizer)
      fPacketizer->SetInitTime();
}

//------------------------------------------------------------------------------


ClassImp(TProofPlayerSlave);

////////////////////////////////////////////////////////////////////////////////
/// Setup feedback.

void TProofPlayerSlave::SetupFeedback()
{
   TList *fb = (TList*) fInput->FindObject("FeedbackList");
   if (fb) {
      PDB(kFeedback,1)
         Info("SetupFeedback","\"FeedbackList\" found: %d objects", fb->GetSize());
   } else {
      PDB(kFeedback,1)
         Info("SetupFeedback","\"FeedbackList\" NOT found");
   }

   if (fb == 0 || fb->GetSize() == 0) return;

   // OK, feedback was requested, setup the timer

   SafeDelete(fFeedbackTimer);
   fFeedbackPeriod = 2000;
   TProof::GetParameter(fInput, "PROOF_FeedbackPeriod", fFeedbackPeriod);
   fFeedbackTimer = new TTimer;
   fFeedbackTimer->SetObject(this);
   fFeedbackTimer->Start(fFeedbackPeriod, kTRUE);

   fFeedback = fb;
}

////////////////////////////////////////////////////////////////////////////////
/// Stop feedback.

void TProofPlayerSlave::StopFeedback()
{
   if (fFeedbackTimer == 0) return;

   PDB(kFeedback,1) Info("StopFeedback","Stop Timer");

   SafeDelete(fFeedbackTimer);
}

////////////////////////////////////////////////////////////////////////////////
/// Handle timer event.

Bool_t TProofPlayerSlave::HandleTimer(TTimer *)
{
   PDB(kFeedback,2) Info("HandleTimer","Entry");

   // If in sequential (0-slave-PROOF) mode we do not have a packetizer
   // so we also send the info to update the progress bar.
   if (gProofServ) {
      Bool_t sendm = kFALSE;
      TMessage m(kPROOF_PROGRESS);
      if (gProofServ->IsMaster() && !gProofServ->IsParallel()) {
         sendm = kTRUE;
         if (gProofServ->GetProtocol() > 25) {
            m << GetProgressStatus();
         } else if (gProofServ->GetProtocol() > 11) {
            TProofProgressStatus *ps = GetProgressStatus();
            m << fTotalEvents << ps->GetEntries() << ps->GetBytesRead()
              << (Float_t) -1. << (Float_t) ps->GetProcTime()
              << (Float_t) ps->GetRate() << (Float_t) -1.;
         } else {
            m << fTotalEvents << GetEventsProcessed();
         }
      }
      if (sendm) gProofServ->GetSocket()->Send(m);
   }

   if (fFeedback == 0) return kFALSE;

   TList *fb = new TList;
   fb->SetOwner(kFALSE);

   if (fOutput == 0) {
      fOutput = (THashList *) fSelector->GetOutputList();
   }

   if (fOutput) {
      TIter next(fFeedback);
      while( TObjString *name = (TObjString*) next() ) {
         // TODO: find object in memory ... maybe allow only in fOutput ?
         TObject *o = fOutput->FindObject(name->GetName());
         if (o != 0) fb->Add(o);
      }
   }

   PDB(kFeedback,2) Info("HandleTimer","Sending %d objects", fb->GetSize());

   TMessage m(kPROOF_FEEDBACK);
   m << fb;

   // send message to client;
   if (gProofServ) gProofServ->GetSocket()->Send(m);

   delete fb;

   fFeedbackTimer->Start(fFeedbackPeriod, kTRUE);

   return kFALSE; // ignored?
}

////////////////////////////////////////////////////////////////////////////////
/// Handle tree header request.

void TProofPlayerSlave::HandleGetTreeHeader(TMessage *mess)
{
   TMessage answ(kPROOF_GETTREEHEADER);

   TDSet *dset;
   (*mess) >> dset;
   dset->Reset();
   TDSetElement *e = dset->Next();
   Long64_t entries = 0;
   TFile *f = 0;
   TTree *t = 0;
   if (!e) {
      PDB(kGlobal, 1) Info("HandleGetTreeHeader", "empty TDSet");
   } else {
      f = TFile::Open(e->GetFileName());
      t = 0;
      if (f) {
         t = (TTree*) f->Get(e->GetObjName());
         if (t) {
            t->SetMaxVirtualSize(0);
            t->DropBaskets();
            entries = t->GetEntries();

            // compute #entries in all the files
            while ((e = dset->Next()) != 0) {
               TFile *f1 = TFile::Open(e->GetFileName());
               if (f1) {
                  TTree *t1 = (TTree*) f1->Get(e->GetObjName());
                  if (t1) {
                     entries += t1->GetEntries();
                     delete t1;
                  }
                  delete f1;
               }
            }
            t->SetMaxEntryLoop(entries);   // this field will hold the total number of entries ;)
         }
      }
   }
   if (t)
      answ << TString("Success") << t;
   else
      answ << TString("Failed") << t;

   fSocket->Send(answ);

   SafeDelete(t);
   SafeDelete(f);
}


//------------------------------------------------------------------------------

ClassImp(TProofPlayerSuperMaster);

////////////////////////////////////////////////////////////////////////////////
/// Process specified TDSet on PROOF. Runs on super master.
/// The return value is -1 in case of error and TSelector::GetStatus() in
/// in case of success.

Long64_t TProofPlayerSuperMaster::Process(TDSet *dset, const char *selector_file,
                                          Option_t *option, Long64_t nentries,
                                          Long64_t first)
{
   fProgressStatus->Reset();
   PDB(kGlobal,1) Info("Process","Enter");

   TProofSuperMaster *proof = dynamic_cast<TProofSuperMaster*>(GetProof());
   if (!proof) return -1;

   delete fOutput;
   fOutput = new THashList;

   TPerfStats::Start(fInput, fOutput);

   if (!SendSelector(selector_file)) {
      Error("Process", "sending selector %s", selector_file);
      return -1;
   }

   TCleanup clean(this);
   SetupFeedback();

   if (proof->IsMaster()) {

      // make sure the DSet is valid
      if (!dset->ElementsValid()) {
         proof->ValidateDSet(dset);
         if (!dset->ElementsValid()) {
            Error("Process", "could not validate TDSet");
            return -1;
         }
      }

      TList msds;
      msds.SetOwner(); // This will delete TPairs

      TList keyholder; // List to clean up key part of the pairs
      keyholder.SetOwner();
      TList valueholder; // List to clean up value part of the pairs
      valueholder.SetOwner();

      // Construct msd list using the slaves
      TIter nextslave(proof->GetListOfActiveSlaves());
      while (TSlave *sl = dynamic_cast<TSlave*>(nextslave())) {
         TList *submasters = 0;
         TPair *msd = dynamic_cast<TPair*>(msds.FindObject(sl->GetMsd()));
         if (!msd) {
            submasters = new TList;
            submasters->SetName(sl->GetMsd());
            keyholder.Add(submasters);
            TList *setelements = new TSortedList(kSortDescending);
            setelements->SetName(TString(sl->GetMsd())+"_Elements");
            valueholder.Add(setelements);
            msds.Add(new TPair(submasters, setelements));
         } else {
            submasters = dynamic_cast<TList*>(msd->Key());
         }
         if (submasters) submasters->Add(sl);
      }

      // Add TDSetElements to msd list
      Long64_t cur = 0; //start of next element
      TIter nextelement(dset->GetListOfElements());
      while (TDSetElement *elem = dynamic_cast<TDSetElement*>(nextelement())) {

         if (elem->GetNum()<1) continue; // get rid of empty elements

         if (nentries !=-1 && cur>=first+nentries) {
            // we are done
            break;
         }

         if (cur+elem->GetNum()-1<first) {
            //element is before first requested entry
            cur+=elem->GetNum();
            continue;
         }

         if (cur<first) {
            //modify element to get proper start
            elem->SetNum(elem->GetNum()-(first-cur));
            elem->SetFirst(elem->GetFirst()+first-cur);
            cur=first;
         }

         if (nentries==-1 || cur+elem->GetNum()<=first+nentries) {
            cur+=elem->GetNum();
         } else {
            //modify element to get proper end
            elem->SetNum(first+nentries-cur);
            cur=first+nentries;
         }

         TPair *msd = dynamic_cast<TPair*>(msds.FindObject(elem->GetMsd()));
         if (!msd) {
            Error("Process", "data requires mass storage domain '%s'"
                  " which is not accessible in this proof session",
                  elem->GetMsd());
            return -1;
         } else {
            TList *elements = dynamic_cast<TList*>(msd->Value());
            if (elements) elements->Add(elem);
         }
      }

      TList usedmasters;
      TIter nextmsd(msds.MakeIterator());
      while (TPair *msd = dynamic_cast<TPair*>(nextmsd())) {
         TList *submasters = dynamic_cast<TList*>(msd->Key());
         TList *setelements = dynamic_cast<TList*>(msd->Value());

         // distribute elements over the masters
         Int_t nmasters = submasters ? submasters->GetSize() : -1;
         Int_t nelements = setelements ? setelements->GetSize() : -1;
         for (Int_t i=0; i<nmasters; i++) {

            Long64_t nent = 0;
            TDSet set(dset->GetType(), dset->GetObjName(),
                      dset->GetDirectory());
            for (Int_t j = (i*nelements)/nmasters;
                       j < ((i+1)*nelements)/nmasters;
                       j++) {
               TDSetElement *elem = setelements ?
                  dynamic_cast<TDSetElement*>(setelements->At(j)) : (TDSetElement *)0;
               if (elem) {
                  set.Add(elem->GetFileName(), elem->GetObjName(),
                        elem->GetDirectory(), elem->GetFirst(),
                        elem->GetNum(), elem->GetMsd());
                  nent += elem->GetNum();
               } else {
                  Warning("Process", "not a TDSetElement object");
               }
            }

            if (set.GetListOfElements()->GetSize()>0) {
               TMessage mesg(kPROOF_PROCESS);
               TString fn(gSystem->BaseName(selector_file));
               TString opt = option;
               mesg << &set << fn << fInput << opt << Long64_t(-1) << Long64_t(0);

               TSlave *sl = dynamic_cast<TSlave*>(submasters->At(i));
               if (sl) {
                  PDB(kGlobal,1) Info("Process",
                                    "Sending TDSet with %d elements to submaster %s",
                                    set.GetListOfElements()->GetSize(),
                                    sl->GetOrdinal());
                  sl->GetSocket()->Send(mesg);
                  usedmasters.Add(sl);

                  // setup progress info
                  fSlaves.AddLast(sl);
                  fSlaveProgress.Set(fSlaveProgress.GetSize()+1);
                  fSlaveProgress[fSlaveProgress.GetSize()-1] = 0;
                  fSlaveTotals.Set(fSlaveTotals.GetSize()+1);
                  fSlaveTotals[fSlaveTotals.GetSize()-1] = nent;
                  fSlaveBytesRead.Set(fSlaveBytesRead.GetSize()+1);
                  fSlaveBytesRead[fSlaveBytesRead.GetSize()-1] = 0;
                  fSlaveInitTime.Set(fSlaveInitTime.GetSize()+1);
                  fSlaveInitTime[fSlaveInitTime.GetSize()-1] = -1.;
                  fSlaveProcTime.Set(fSlaveProcTime.GetSize()+1);
                  fSlaveProcTime[fSlaveProcTime.GetSize()-1] = -1.;
                  fSlaveEvtRti.Set(fSlaveEvtRti.GetSize()+1);
                  fSlaveEvtRti[fSlaveEvtRti.GetSize()-1] = -1.;
                  fSlaveMBRti.Set(fSlaveMBRti.GetSize()+1);
                  fSlaveMBRti[fSlaveMBRti.GetSize()-1] = -1.;
                  fSlaveActW.Set(fSlaveActW.GetSize()+1);
                  fSlaveActW[fSlaveActW.GetSize()-1] = 0;
                  fSlaveTotS.Set(fSlaveTotS.GetSize()+1);
                  fSlaveTotS[fSlaveTotS.GetSize()-1] = 0;
                  fSlaveEffS.Set(fSlaveEffS.GetSize()+1);
                  fSlaveEffS[fSlaveEffS.GetSize()-1] = 0.;
               } else {
                  Warning("Process", "not a TSlave object");
               }
            }
         }
      }

      if ( !IsClient() ) HandleTimer(0);
      PDB(kGlobal,1) Info("Process","Calling Collect");
      proof->Collect(&usedmasters);
      HandleTimer(0);

   }

   StopFeedback();

   PDB(kGlobal,1) Info("Process","Calling Merge Output");
   MergeOutput();

   TPerfStats::Stop();

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Report progress.

void TProofPlayerSuperMaster::Progress(TSlave *sl, Long64_t total, Long64_t processed)
{
   Int_t idx = fSlaves.IndexOf(sl);
   fSlaveProgress[idx] = processed;
   if (fSlaveTotals[idx] != total)
      Warning("Progress", "total events has changed for slave %s", sl->GetName());
   fSlaveTotals[idx] = total;

   Long64_t tot = 0;
   Int_t i;
   for (i = 0; i < fSlaveTotals.GetSize(); i++) tot += fSlaveTotals[i];
   Long64_t proc = 0;
   for (i = 0; i < fSlaveProgress.GetSize(); i++) proc += fSlaveProgress[i];

   Progress(tot, proc);
}

////////////////////////////////////////////////////////////////////////////////
/// Report progress.

void TProofPlayerSuperMaster::Progress(TSlave *sl, Long64_t total,
                                       Long64_t processed, Long64_t bytesread,
                                       Float_t initTime, Float_t procTime,
                                       Float_t evtrti, Float_t mbrti)
{
   PDB(kGlobal,2)
      Info("Progress","%s: %lld %lld %f %f %f %f", sl->GetName(),
                      processed, bytesread, initTime, procTime, evtrti, mbrti);

   Int_t idx = fSlaves.IndexOf(sl);
   if (fSlaveTotals[idx] != total)
      Warning("Progress", "total events has changed for slave %s", sl->GetName());
   fSlaveTotals[idx] = total;
   fSlaveProgress[idx] = processed;
   fSlaveBytesRead[idx] = bytesread;
   fSlaveInitTime[idx] = (initTime > -1.) ? initTime : fSlaveInitTime[idx];
   fSlaveProcTime[idx] = (procTime > -1.) ? procTime : fSlaveProcTime[idx];
   fSlaveEvtRti[idx] = (evtrti > -1.) ? evtrti : fSlaveEvtRti[idx];
   fSlaveMBRti[idx] = (mbrti > -1.) ? mbrti : fSlaveMBRti[idx];

   Int_t i;
   Long64_t tot = 0;
   Long64_t proc = 0;
   Long64_t bytes = 0;
   Float_t init = -1.;
   Float_t ptime = -1.;
   Float_t erti = 0.;
   Float_t srti = 0.;
   Int_t nerti = 0;
   Int_t nsrti = 0;
   for (i = 0; i < fSlaveTotals.GetSize(); i++) {
      tot += fSlaveTotals[i];
      if (i < fSlaveProgress.GetSize())
         proc += fSlaveProgress[i];
      if (i < fSlaveBytesRead.GetSize())
         bytes += fSlaveBytesRead[i];
      if (i < fSlaveInitTime.GetSize())
         if (fSlaveInitTime[i] > -1. && (init < 0. || fSlaveInitTime[i] < init))
            init = fSlaveInitTime[i];
      if (i < fSlaveProcTime.GetSize())
         if (fSlaveProcTime[i] > -1. && (ptime < 0. || fSlaveProcTime[i] > ptime))
            ptime = fSlaveProcTime[i];
      if (i < fSlaveEvtRti.GetSize())
         if (fSlaveEvtRti[i] > -1.) {
            erti += fSlaveEvtRti[i];
            nerti++;
         }
      if (i < fSlaveMBRti.GetSize())
         if (fSlaveMBRti[i] > -1.) {
            srti += fSlaveMBRti[i];
            nsrti++;
         }
   }
   srti = (nsrti > 0) ? srti / nerti : 0.;

   Progress(tot, proc, bytes, init, ptime, erti, srti);
}

////////////////////////////////////////////////////////////////////////////////
/// Progress signal.

void TProofPlayerSuperMaster::Progress(TSlave *wrk, TProofProgressInfo *pi)
{
   if (pi) {
      PDB(kGlobal,2)
         Info("Progress","%s: %lld %lld %lld %f %f %f %f %d %f", wrk->GetOrdinal(),
                         pi->fTotal, pi->fProcessed, pi->fBytesRead,
                         pi->fInitTime, pi->fProcTime, pi->fEvtRateI, pi->fMBRateI,
                         pi->fActWorkers, pi->fEffSessions);

      Int_t idx = fSlaves.IndexOf(wrk);
      if (fSlaveTotals[idx] != pi->fTotal)
         Warning("Progress", "total events has changed for worker %s", wrk->GetName());
      fSlaveTotals[idx] = pi->fTotal;
      fSlaveProgress[idx] = pi->fProcessed;
      fSlaveBytesRead[idx] = pi->fBytesRead;
      fSlaveInitTime[idx] = (pi->fInitTime > -1.) ? pi->fInitTime : fSlaveInitTime[idx];
      fSlaveProcTime[idx] = (pi->fProcTime > -1.) ? pi->fProcTime : fSlaveProcTime[idx];
      fSlaveEvtRti[idx] = (pi->fEvtRateI > -1.) ? pi->fEvtRateI : fSlaveEvtRti[idx];
      fSlaveMBRti[idx] = (pi->fMBRateI > -1.) ? pi->fMBRateI : fSlaveMBRti[idx];
      fSlaveActW[idx] = (pi->fActWorkers > -1) ? pi->fActWorkers : fSlaveActW[idx];
      fSlaveTotS[idx] = (pi->fTotSessions > -1) ? pi->fTotSessions : fSlaveTotS[idx];
      fSlaveEffS[idx] = (pi->fEffSessions > -1.) ? pi->fEffSessions : fSlaveEffS[idx];

      Int_t i;
      Int_t nerti = 0;
      Int_t nsrti = 0;
      TProofProgressInfo pisum(0, 0, 0, -1., -1., 0., 0., 0, 0, 0.);
      for (i = 0; i < fSlaveTotals.GetSize(); i++) {
         pisum.fTotal += fSlaveTotals[i];
         if (i < fSlaveProgress.GetSize())
            pisum.fProcessed += fSlaveProgress[i];
         if (i < fSlaveBytesRead.GetSize())
            pisum.fBytesRead += fSlaveBytesRead[i];
         if (i < fSlaveInitTime.GetSize())
            if (fSlaveInitTime[i] > -1. && (pisum.fInitTime < 0. || fSlaveInitTime[i] < pisum.fInitTime))
               pisum.fInitTime = fSlaveInitTime[i];
         if (i < fSlaveProcTime.GetSize())
            if (fSlaveProcTime[i] > -1. && (pisum.fProcTime < 0. || fSlaveProcTime[i] > pisum.fProcTime))
               pisum.fProcTime = fSlaveProcTime[i];
         if (i < fSlaveEvtRti.GetSize())
            if (fSlaveEvtRti[i] > -1.) {
               pisum.fEvtRateI += fSlaveEvtRti[i];
               nerti++;
            }
         if (i < fSlaveMBRti.GetSize())
            if (fSlaveMBRti[i] > -1.) {
               pisum.fMBRateI += fSlaveMBRti[i];
               nsrti++;
            }
         if (i < fSlaveActW.GetSize())
            pisum.fActWorkers += fSlaveActW[i];
         if (i < fSlaveTotS.GetSize())
            if (fSlaveTotS[i] > -1 && (pisum.fTotSessions < 0. || fSlaveTotS[i] > pisum.fTotSessions))
               pisum.fTotSessions = fSlaveTotS[i];
         if (i < fSlaveEffS.GetSize())
            if (fSlaveEffS[i] > -1. && (pisum.fEffSessions < 0. || fSlaveEffS[i] > pisum.fEffSessions))
               pisum.fEffSessions = fSlaveEffS[i];
      }
      pisum.fMBRateI = (nsrti > 0) ? pisum.fMBRateI / nerti : 0.;

      Progress(&pisum);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Send progress and feedback to client.

Bool_t TProofPlayerSuperMaster::HandleTimer(TTimer *)
{
   if (fFeedbackTimer == 0) return kFALSE; // timer stopped already

   Int_t i;
   Long64_t tot = 0;
   Long64_t proc = 0;
   Long64_t bytes = 0;
   Float_t init = -1.;
   Float_t ptime = -1.;
   Float_t erti = 0.;
   Float_t srti = 0.;
   Int_t nerti = 0;
   Int_t nsrti = 0;
   for (i = 0; i < fSlaveTotals.GetSize(); i++) {
      tot += fSlaveTotals[i];
      if (i < fSlaveProgress.GetSize())
         proc += fSlaveProgress[i];
      if (i < fSlaveBytesRead.GetSize())
         bytes += fSlaveBytesRead[i];
      if (i < fSlaveInitTime.GetSize())
         if (fSlaveInitTime[i] > -1. && (init < 0. || fSlaveInitTime[i] < init))
            init = fSlaveInitTime[i];
      if (i < fSlaveProcTime.GetSize())
         if (fSlaveProcTime[i] > -1. && (ptime < 0. || fSlaveProcTime[i] > ptime))
            ptime = fSlaveProcTime[i];
      if (i < fSlaveEvtRti.GetSize())
         if (fSlaveEvtRti[i] > -1.) {
            erti += fSlaveEvtRti[i];
            nerti++;
         }
      if (i < fSlaveMBRti.GetSize())
         if (fSlaveMBRti[i] > -1.) {
            srti += fSlaveMBRti[i];
            nsrti++;
         }
   }
   erti = (nerti > 0) ? erti / nerti : 0.;
   srti = (nsrti > 0) ? srti / nerti : 0.;

   TMessage m(kPROOF_PROGRESS);
   if (gProofServ->GetProtocol() > 25) {
      // Fill the message now
      TProofProgressInfo pi(tot, proc, bytes, init, ptime,
                            erti, srti, -1,
                            gProofServ->GetTotSessions(), gProofServ->GetEffSessions());
      m << &pi;
   } else {

      m << tot << proc << bytes << init << ptime << erti << srti;
   }

   // send message to client;
   gProofServ->GetSocket()->Send(m);

   if (fReturnFeedback)
      return TProofPlayerRemote::HandleTimer(0);
   else
      return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Setup reporting of feedback objects and progress messages.

void TProofPlayerSuperMaster::SetupFeedback()
{
   if (IsClient()) return; // Client does not need timer

   TProofPlayerRemote::SetupFeedback();

   if (fFeedbackTimer) {
      fReturnFeedback = kTRUE;
      return;
   } else {
      fReturnFeedback = kFALSE;
   }

   // setup the timer for progress message
   SafeDelete(fFeedbackTimer);
   fFeedbackPeriod = 2000;
   TProof::GetParameter(fInput, "PROOF_FeedbackPeriod", fFeedbackPeriod);
   fFeedbackTimer = new TTimer;
   fFeedbackTimer->SetObject(this);
   fFeedbackTimer->Start(fFeedbackPeriod, kTRUE);
}
