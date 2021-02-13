// @(#)root/proofplayer:$Id$
// Author: Long Tran-Thanh    22/07/07
// Revised: G. Ganis, May 2011

/*************************************************************************
 * Copyright (C) 1995-2002, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TPacketizerUnit
\ingroup proofkernel

This packetizer generates packets of generic units, representing the
number of times an operation cycle has to be repeated by the worker
node, e.g. the number of Monte carlo events to be generated.
Packets sizes are generated taking into account the performance of
worker nodes, based on the time needed to process previous packets,
with the goal of having all workers ending at the same time.

*/


#include "TPacketizerUnit.h"

#include "Riostream.h"
#include "TDSet.h"
#include "TError.h"
#include "TEventList.h"
#include "TMap.h"
#include "TMessage.h"
#include "TMonitor.h"
#include "TNtupleD.h"
#include "TObject.h"
#include "TParameter.h"
#include "TPerfStats.h"
#include "TProofDebug.h"
#include "TProof.h"
#include "TProofPlayer.h"
#include "TProofServ.h"
#include "TSlave.h"
#include "TSocket.h"
#include "TStopwatch.h"
#include "TTimer.h"
#include "TUrl.h"
#include "TClass.h"
#include "TMath.h"


using namespace TMath;
//
// The following utility class manage the state of the
// work to be performed and the slaves involved in the process.
//
// The list of TSlaveStat(s) keep track of the work (being) done
// by each slave
//

//------------------------------------------------------------------------------

class TPacketizerUnit::TSlaveStat : public TVirtualPacketizer::TVirtualSlaveStat {

friend class TPacketizerUnit;

private:
   Long64_t  fLastProcessed; // Number of processed entries of the last packet
   Double_t  fRate;         // Estimated processing rate averaged over circularity
   Double_t  fTimeInstant;   // Starting time of the current packet
   TNtupleD *fCircNtp;       // Keeps circular info for speed calculations
   Long_t    fCircLvl;       // Circularity level

public:
   TSlaveStat(TSlave *sl, TList *input);
   ~TSlaveStat();

//   void        GetCurrentTime();

   void        UpdatePerformance(Double_t time);
   TProofProgressStatus *AddProcessed(TProofProgressStatus *st);

//   ClassDef(TPacketizerUnit::TSlaveStat, 0);
};

////////////////////////////////////////////////////////////////////////////////
/// Main constructor

TPacketizerUnit::TSlaveStat::TSlaveStat(TSlave *slave, TList *input)
                            : fLastProcessed(0),
                              fRate(0), fTimeInstant(0), fCircLvl(5)
{
   // Initialize the circularity ntple for speed calculations
   fCircNtp = new TNtupleD("Speed Circ Ntp", "Circular process info","tm:ev");
   fCircNtp->SetDirectory(0);
   TProof::GetParameter(input, "PROOF_TPacketizerUnitCircularity", fCircLvl);
   fCircLvl = (fCircLvl > 0) ? fCircLvl : 5;
   fCircNtp->SetCircular(fCircLvl);
   fSlave = slave;
   fStatus = new TProofProgressStatus();
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TPacketizerUnit::TSlaveStat::~TSlaveStat()
{
   SafeDelete(fCircNtp);
}

////////////////////////////////////////////////////////////////////////////////
/// Update the circular ntple

void TPacketizerUnit::TSlaveStat::UpdatePerformance(Double_t time)
{
   Double_t ttot = time;
   Double_t *ar = fCircNtp->GetArgs();
   Int_t ne = fCircNtp->GetEntries();
   if (ne <= 0) {
      // First call: just fill one ref entry and return
      fCircNtp->Fill(0., 0);
      fRate = 0.;
      return;
   }
   // Fill the entry
   fCircNtp->GetEntry(ne-1);
   ttot = ar[0] + time;
   fCircNtp->Fill(ttot, GetEntriesProcessed());

   // Calculate the speed
   fCircNtp->GetEntry(0);
   Double_t dtime = (ttot > ar[0]) ? ttot - ar[0] : ne+1 ;
   Long64_t nevts = GetEntriesProcessed() - (Long64_t)ar[1];
   fRate = nevts / dtime;
   PDB(kPacketizer,2)
      Info("UpdatePerformance", "time:%f, dtime:%f, nevts:%lld, speed: %f",
                                time, dtime, nevts, fRate);

}

////////////////////////////////////////////////////////////////////////////////
/// Update the status info to the 'st'.
/// return the difference (*st - *fStatus)

TProofProgressStatus *TPacketizerUnit::TSlaveStat::AddProcessed(TProofProgressStatus *st)
{
   if (st) {
      // The entriesis not correct in 'st'
      Long64_t lastEntries = st->GetEntries() - fStatus->GetEntries();
      // The last proc time should not be added
      fStatus->SetLastProcTime(0.);
      // Get the diff
      TProofProgressStatus *diff = new TProofProgressStatus(*st - *fStatus);
      *fStatus += *diff;
      // Set the correct value
      fStatus->SetLastEntries(lastEntries);
      return diff;
   } else {
      Error("AddProcessed", "status arg undefined");
      return 0;
   }
}

//------------------------------------------------------------------------------

ClassImp(TPacketizerUnit);

////////////////////////////////////////////////////////////////////////////////
/// Constructor

TPacketizerUnit::TPacketizerUnit(TList *slaves, Long64_t num, TList *input,
                                 TProofProgressStatus *st)
                : TVirtualPacketizer(input, st)
{
   PDB(kPacketizer,1) Info("TPacketizerUnit", "enter (num %lld)", num);

   // Init pointer members
   fWrkStats = 0;
   fPackets = 0;
   fInput = input;

   fFixedNum = kFALSE;
   Int_t fixednum = -1;
   if (TProof::GetParameter(input, "PROOF_PacketizerFixedNum", fixednum) != 0 || fixednum <= 0) {
      fFixedNum = kFALSE;
   }
   else {
      Info("TPacketizerUnit", "forcing the same cycles on each worker");
      fFixedNum = kTRUE;
   }

   fCalibFrac = 0.01;
   if (TProof::GetParameter(input, "PROOF_PacketizerCalibFrac", fCalibFrac) != 0 || fCalibFrac <= 0)
      fCalibFrac = 0.01;
   PDB(kPacketizer,1)
      Info("TPacketizerUnit", "size of the calibration packets: %.2f %% of average number per worker", fCalibFrac);

   fMaxPacketTime = 3.;
   Double_t timeLimit = -1;
   if (TProof::GetParameter(input, "PROOF_PacketizerTimeLimit", timeLimit) == 0) {
      fMaxPacketTime = timeLimit;
      Warning("TPacketizerUnit", "PROOF_PacketizerTimeLimit is deprecated: use PROOF_MaxPacketTime instead");
   }
   PDB(kPacketizer,1)
      Info("TPacketizerUnit", "time limit is %lf", fMaxPacketTime);

   // Different default for min packet time
   fMinPacketTime = 1;
   Double_t minPacketTime = 0;
   if (TProof::GetParameter(input, "PROOF_MinPacketTime", minPacketTime) == 0) fMinPacketTime = minPacketTime;
   TParameter<Double_t> *mpt = (TParameter<Double_t> *) fConfigParams->FindObject("PROOF_MinPacketTime");
   if (mpt) {
      mpt->SetVal(fMinPacketTime);
   } else {
      fConfigParams->Add(new TParameter<Double_t>("PROOF_MinPacketTime", fMinPacketTime));
   }

   fProcessing = 0;
   fAssigned = 0;
   fPacketSeq = 0;

   fStopwatch = new TStopwatch();

   fPackets = new TList;
   fPackets->SetOwner();

   fWrkStats = new TMap;
   fWrkStats->SetOwner(kFALSE);
   fWrkExcluded = 0;

   TSlave *slave;
   TIter si(slaves);
   while ((slave = (TSlave*) si.Next())) {
      if (slave->GetParallel() > 0) {
         fWrkStats->Add(slave, new TSlaveStat(slave, input));
      } else {
         if (!fWrkExcluded) {
            fWrkExcluded = new TList;
            fWrkExcluded->SetOwner(kFALSE);
         }
         PDB(kPacketizer,2)
            Info("TPacketizerUnit", "node '%s' has NO active worker: excluded from work distribution", slave->GetOrdinal());
         fWrkExcluded->Add(slave);
      }
   }

   fTotalEntries = 0;
   fNumPerWorker = -1;
   if (num > 0 && AssignWork(0,0,num) != 0)
      Warning("TPacketizerUnit", "some problems assigning work");

   // Save the config parameters in the dedicated list so that they will be saved
   // in the outputlist and therefore in the relevant TQueryResult
   fConfigParams->Add(new TParameter<Float_t>("PROOF_PacketizerCalibFrac", fCalibFrac));

   fStopwatch->Start();
   PDB(kPacketizer,1) Info("TPacketizerUnit", "return");
}

////////////////////////////////////////////////////////////////////////////////
/// Assign work to be done to this packetizer

Int_t TPacketizerUnit::AssignWork(TDSet *, Long64_t, Long64_t num)
{
   if (num < 0) {
      Error("AssignWork", "assigned a negative number (%lld) of cycles - protocol error?", num);
      return -1;
   }

   fTotalEntries += num;
   PDB(kPacketizer,1)
      Info("AssignWork", "assigned %lld additional cycles (new total: %lld)", num, fTotalEntries);

   // Update fixed number counter
   if (fFixedNum && fWrkStats->GetSize() > 0) {
      // Approximate number: the exact number is determined in GetNextPacket
      fNumPerWorker = fTotalEntries / fWrkStats->GetSize();
      if (fNumPerWorker == 0) fNumPerWorker = 1;
   }

   // Update/Save the config parameters in the dedicated list so that they will be saved
   // in the outputlist and therefore in the relevant TQueryResult
   TParameter<Long64_t> *fn =
      (TParameter<Long64_t> *) fConfigParams->FindObject("PROOF_PacketizerFixedNum");
   if (fn) {
      fn->SetVal(fNumPerWorker);
   } else {
      fConfigParams->Add(new TParameter<Long64_t>("PROOF_PacketizerFixedNum", fNumPerWorker));
   }

   // Done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TPacketizerUnit::~TPacketizerUnit()
{
   if (fWrkStats)
      fWrkStats->DeleteValues();
   SafeDelete(fWrkStats);
   SafeDelete(fWrkExcluded);
   SafeDelete(fPackets);
   SafeDelete(fStopwatch);
}

////////////////////////////////////////////////////////////////////////////////
/// Get current time

Double_t TPacketizerUnit::GetCurrentTime()
{
   Double_t retValue = fStopwatch->RealTime();
   fStopwatch->Continue();
   return retValue;
}

////////////////////////////////////////////////////////////////////////////////
/// Get Estimation of the current rate; just summing the current rates of
/// the active workers

Float_t TPacketizerUnit::GetCurrentRate(Bool_t &all)
{
   all = kTRUE;
   // Loop over the workers
   Float_t currate = 0.;
   if (fWrkStats && fWrkStats->GetSize() > 0) {
      TIter nxw(fWrkStats);
      TObject *key;
      while ((key = nxw()) != 0) {
         TSlaveStat *slstat = (TSlaveStat *) fWrkStats->GetValue(key);
         if (slstat && slstat->GetProgressStatus() && slstat->GetEntriesProcessed() > 0) {
            // Sum-up the current rates
            currate += slstat->GetProgressStatus()->GetCurrentRate();
         } else {
            all = kFALSE;
         }
      }
   }
   // Done
   return currate;
}

////////////////////////////////////////////////////////////////////////////////
/// Get next packet

TDSetElement *TPacketizerUnit::GetNextPacket(TSlave *sl, TMessage *r)
{
   if (!fValid)
      return 0;

   // Find slave
   TSlaveStat *slstat = (TSlaveStat*) fWrkStats->GetValue(sl);
   if (!slstat) {
      Warning("GetNextPacket", "Received a packet request from an unknown slave: %s:%s",
         sl->GetName(), sl->GetOrdinal());
      return 0;
   }

   PDB(kPacketizer,2)
      Info("GetNextPacket","worker-%s: fAssigned %lld\t", sl->GetOrdinal(), fAssigned);

   // Update stats & free old element
   Double_t latency = 0., proctime = 0., proccpu = 0.;
   Long64_t bytesRead = -1;
   Long64_t totalEntries = -1; // used only to read an old message type
   Long64_t totev = 0;
   Long64_t numev = -1;

   TProofProgressStatus *status = 0;
   if (sl->GetProtocol() > 18) {
      (*r) >> latency;
      (*r) >> status;

      // Calculate the progress made in the last packet
      TProofProgressStatus *progress = 0;
      if (status) {
         // update the worker status
         numev = status->GetEntries() - slstat->GetEntriesProcessed();
         progress = slstat->AddProcessed(status);
         if (progress) {
            // (*fProgressStatus) += *progress;
            proctime = progress->GetProcTime();
            proccpu  = progress->GetCPUTime();
            totev  = status->GetEntries(); // for backward compatibility
            bytesRead  = progress->GetBytesRead();
            delete progress;
         }
         delete status;
      } else
          Error("GetNextPacket", "no status came in the kPROOF_GETPACKET message");
   } else {

      (*r) >> latency >> proctime >> proccpu;

      // only read new info if available
      if (r->BufferSize() > r->Length()) (*r) >> bytesRead;
      if (r->BufferSize() > r->Length()) (*r) >> totalEntries;
      if (r->BufferSize() > r->Length()) (*r) >> totev;

      numev = totev - slstat->GetEntriesProcessed();
      slstat->GetProgressStatus()->IncEntries(numev);
      slstat->GetProgressStatus()->SetLastUpdate();
   }

   fProgressStatus->IncEntries(numev);
   fProgressStatus->SetLastUpdate();

   fProcessing = 0;

   PDB(kPacketizer,2)
      Info("GetNextPacket","worker-%s (%s): %lld %7.3lf %7.3lf %7.3lf %lld",
                           sl->GetOrdinal(), sl->GetName(),
                           numev, latency, proctime, proccpu, bytesRead);

   if (gPerfStats != 0) {
      gPerfStats->PacketEvent(sl->GetOrdinal(), sl->GetName(), "", numev,
                              latency, proctime, proccpu, bytesRead);
   }

   if (fNumPerWorker > 0 && slstat->GetEntriesProcessed() >= fNumPerWorker) {
      PDB(kPacketizer,2)
         Info("GetNextPacket","worker-%s (%s) is done (%lld cycles)",
                           sl->GetOrdinal(), sl->GetName(), slstat->GetEntriesProcessed());
      return 0;
   }

   if (fAssigned == fTotalEntries) {
      Bool_t done = kTRUE;
      // If we are on a submaster, check if there is something else to do
      if (gProofServ && gProofServ->IsMaster() && !gProofServ->IsTopMaster()) {
         TDSetElement *nxe = gProofServ->GetNextPacket();
         if (nxe) {
            if (AssignWork(0,0,nxe->GetNum()) == 0) {
               if (fAssigned < fTotalEntries) done = kFALSE;
            } else {
               Error("GetNextPacket", "problems assigning additional work: stop");
            }
            SafeDelete(nxe);
         }
      }
      if (done) {
         // Send last timer message
         HandleTimer(0);
         return 0;
      }
   }

   if (fStop) {
      // Send last timer message
      HandleTimer(0);
      return 0;
   }


   Long64_t num;

   // Get the current time
   Double_t cTime = GetCurrentTime();

   if (slstat->fCircNtp->GetEntries() <= 0) {
      // The calibration phase
      Long64_t avg = fTotalEntries / fWrkStats->GetSize();
      num = (Long64_t) (fCalibFrac * avg);
      if (num < 1) num = (avg >= 1) ? avg : 1;
      PDB(kPacketizer,2)
         Info("GetNextPacket", "calibration: total entries %lld, workers %d, frac: %.1f %%, raw num: %lld",
                               fTotalEntries, fWrkStats->GetSize(), fCalibFrac * 100., num);

      // Create a reference entry
      slstat->UpdatePerformance(0.);

   } else {

      if (fNumPerWorker < 0) {

         // Schedule tasks for workers based on the currently estimated processing speeds

         // Update performances
         // slstat->fStatus was updated before;
         slstat->UpdatePerformance(proctime);

         // We need to estimate the total instantaneous rate: for the workers not having yet
         // one we assume the average of those having a measurement
         // The optimal number for worker j is
         //
         //                      n_j =  r_j / Sum r_i * N_left
         //

         Int_t nrm = 0;
         Double_t sumRate = 0.;
         TIter nxwrk(fWrkStats);
         TSlaveStat *wrkStat = 0;
         TSlave *tmpWrk = 0;
         while ((tmpWrk = (TSlave *)nxwrk())) {
            if ((wrkStat = dynamic_cast<TSlaveStat *>(fWrkStats->GetValue(tmpWrk)))) {
               if (wrkStat->fRate > 0) {
                  nrm++;
                  sumRate += wrkStat->fRate;
               }
               PDB(kPacketizer,3)
                  Info("GetNextPacket", "%d: worker-%s: rate %lf /s (sum: %lf /s)",
                                          nrm, tmpWrk->GetOrdinal(), wrkStat->fRate, sumRate);
            } else {
               Warning("GetNextPacket", "dynamic_cast<TSlaveStat *> failing on value for '%s (%s)'! Skipping",
                                        tmpWrk->GetName(), tmpWrk->GetOrdinal());
            }
         }

         // Check consistency
         if (nrm <= 0) {
            Error("GetNextPacket", "no worker has consistent information: stop processing!");
            return (TDSetElement *)0;
         }

         Double_t avgRate = sumRate / nrm;
         // Check if all workers had meaningful rate information
         if (nrm < fWrkStats->GetSize()) {
            // For some workers the measurement is missing: use the average
            sumRate += (fWrkStats->GetSize() - nrm) * avgRate;
         }
         PDB(kPacketizer,2)
            Info("GetNextPacket", "rate: avg: %lf /s/wrk - sum: %lf /s (measurements %d out of %d)",
                                   avgRate, sumRate, nrm, fWrkStats->GetSize());

         // Packet size for this worker
         Double_t wrkRate = (slstat->fRate > 0.) ? slstat->fRate : avgRate ;
         num = (Long64_t) ((fTotalEntries - fAssigned) * wrkRate / sumRate);
         PDB(kPacketizer,2)
            Info("GetNextPacket", "worker-%s (%s): raw packet size: %lld", sl->GetOrdinal(), sl->GetName(), num);

         // Apply time-per-packet limits
         Double_t packTime = num / wrkRate;
         if (fMaxPacketTime > 0. && packTime > fMaxPacketTime) {
            num = (Long64_t) (fMaxPacketTime * wrkRate) ;
            packTime = fMaxPacketTime;
            PDB(kPacketizer,2)
               Info("GetNextPacket", "worker-%s (%s): time-limited packet size: %lld (upper limit: %.2f secs)",
                                     sl->GetOrdinal(), sl->GetName(), num, fMaxPacketTime);
         }
         if (fMinPacketTime > 0. && packTime < fMinPacketTime) {
            num = (Long64_t) (fMinPacketTime * wrkRate);
            PDB(kPacketizer,2)
               Info("GetNextPacket", "worker-%s (%s): time-limited packet size: %lld (lower limit: %.2f secs)",
                                     sl->GetOrdinal(), sl->GetName(), num, fMinPacketTime);
         }

      } else {
         // Fixed number of cycles per worker
         num = fNumPerWorker - slstat->fLastProcessed;
         if (num > 1 && slstat->fRate > 0 && num / slstat->fRate > fMaxPacketTime) {
            num = (Long64_t) (slstat->fRate * fMaxPacketTime);
         }
      }
   }
   // Minimum packet size
   num = (num > 1) ? num : 1;
   fProcessing = (num < (fTotalEntries - fAssigned)) ? num
                                                     : (fTotalEntries - fAssigned);

   // Set the information of the current slave
   slstat->fLastProcessed = fProcessing;
   // Set the start time of the current packet
   slstat->fTimeInstant = cTime;

   // Update the sequential number
   fPacketSeq++;
   TString sseq = TString::Format("p%lld", fPacketSeq);

   PDB(kPacketizer,2)
      Info("GetNextPacket", "worker-%s: num %lld, processing %lld, remaining %lld",sl->GetOrdinal(),
                            num, fProcessing, (fTotalEntries - fAssigned - fProcessing));
   TDSetElement *elem = new TDSetElement(sseq, sseq, "", fAssigned, fProcessing);
   elem->SetBit(TDSetElement::kEmpty);

   // Update the total counter
   fAssigned += slstat->fLastProcessed;

   return elem;
}

////////////////////////////////////////////////////////////////////////////////
/// Adds new workers. Returns the number of workers added, or -1 on failure.

Int_t TPacketizerUnit::AddWorkers(TList *workers)
{
   if (!workers) {
      Error("AddWorkers", "Null list of new workers!");
      return -1;
   }

   Int_t curNumOfWrks = fWrkStats->GetEntries();

   TSlave *sl;
   TIter next(workers);
   while (( sl = dynamic_cast<TSlave*>(next()) ))
      fWrkStats->Add(sl, new TSlaveStat(sl, fInput));

   fNumPerWorker = -1;
   if (fFixedNum && fWrkStats->GetSize() > 0) {
      // Approximate number: the exact number is determined in GetNextPacket
      fNumPerWorker = (fNumPerWorker * curNumOfWrks) / fWrkStats->GetSize();
      if (fNumPerWorker == 0) fNumPerWorker = 1;
   }

   fConfigParams->Add(new TParameter<Long64_t>("PROOF_PacketizerFixedNum", fNumPerWorker));

   return fWrkStats->GetEntries();
}
