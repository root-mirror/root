// @(#)root/proofplayer:$Id$
// Author: Maarten Ballintijn    9/7/2002

/*************************************************************************
 * Copyright (C) 1995-2002, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TVirtualPacketizer
#define ROOT_TVirtualPacketizer

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TVirtualPacketizer                                                   //
//                                                                      //
// Packetizer is a load balancing object created for each query.        //
// It generates packets to be processed on PROOF worker servers.        //
// A packet is an event range (begin entry and number of entries) or    //
// object range (first object and number of objects) in a TTree         //
// (entries) or a directory (objects) in a file.                        //
// Packets are generated taking into account the performance of the     //
// remote machine, the time it took to process a previous packet on     //
// the remote machine, the locality of the database files, etc.         //
//                                                                      //
// TVirtualPacketizer includes common parts of PROOF packetizers.       //
// Look in subclasses for details.                                      //
// The default packetizer is TPacketizerAdaptive.                       //
// To use an alternative one, for instance - the TPacketizer, call:     //
// proof->SetParameter("PROOF_Packetizer", "TPacketizer");              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "TSlave.h"
#include "TProofProgressStatus.h"
#include "TTime.h"


class TDSet;
class TDSetElement;
class TList;
class TMap;
class TMessage;
class TNtuple;
class TNtupleD;
class TProofProgressInfo;
class TSlave;


class TVirtualPacketizer : public TObject {

public:              // public because of Sun CC bug
   class TVirtualSlaveStat;

protected:
   enum EUseEstOpt {        // Option for usage of estimated values
      kEstOff     = 0,
      kEstCurrent = 1,
      kEstAverage = 2
   };

   // General configuration parameters
   Double_t  fMinPacketTime; // minimum packet time
   Double_t  fMaxPacketTime; // maximum packet time
   TList    *fConfigParams;  // List of configuration parameters

   TMap     *fSlaveStats;   // slave status, keyed by correspondig TSlave

   TProofProgressStatus *fProgressStatus; // pointer to status in the player.
   TTimer   *fProgress;     // progress updates timer

   Long64_t  fTotalEntries; // total number of entries to be distributed;
                            // not used in the progressive packetizer
   TList    *fFailedPackets;// a list of packets that failed while processing

   // Members for progress info
   TTime     fStartTime;    // time offset
   Float_t   fInitTime;     // time before processing
   Float_t   fProcTime;     // time since start of processing
   Float_t   fTimeUpdt;     // time between updates
   TNtupleD *fCircProg;     // Keeps circular info for "instantenous"
                            // rate calculations
   Long_t    fCircN;        // Circularity

   TNtuple  *fProgressPerf; // {Active workers, evt rate, MBs read} as a function of processing time
   Float_t   fProcTimeLast; // Time of the last measurement
   Int_t     fActWrksLast;  // Active workers at fProcTimeLast
   Float_t   fEvtRateLast;  // Evt rate at fProcTimeLast
   Float_t   fMBsReadLast;  // MBs read at fProcTimeLast
   Float_t   fEffSessLast;  // Number of effective sessions at fProcTimeLast
   Bool_t    fAWLastFill;   // Whether to fill the last measurement
   Float_t   fReportPeriod; // Time between reports if nothing changes (estimated proc time / 100)

   EUseEstOpt fUseEstOpt;   // Control usage of estimated values for the progress info

   Bool_t   fValid;           // Constructed properly?
   Bool_t   fStop;            // Termination of Process() requested?

   TString  fDataSet;         // Name of the dataset being processed (for dataset-driven runs)

   TList      *fInput;        // Input list

   TVirtualPacketizer(TList *input, TProofProgressStatus *st = 0);
   TVirtualPacketizer(const TVirtualPacketizer &);  // no implementation, will generate
   void operator=(const TVirtualPacketizer &);      // error on accidental usage

   TDSetElement  *CreateNewPacket(TDSetElement* base, Long64_t first, Long64_t num);
   Long64_t       GetEntries(Bool_t tree, TDSetElement *e); // Num of entries or objects
   virtual Bool_t HandleTimer(TTimer *timer);

public:
   enum EStatusBits { kIsInitializing = BIT(16), kIsDone = BIT(17), kIsTree = BIT(18) };
   virtual ~TVirtualPacketizer();

   virtual Int_t           AssignWork(TDSet* /*dset*/, Long64_t /*first*/, Long64_t /*num*/) { return -1; }
   Bool_t                  IsValid() const { return fValid; }
   Long64_t                GetEntriesProcessed() const { return (fProgressStatus? fProgressStatus->GetEntries() : 0); }
   virtual Int_t           GetEstEntriesProcessed(Float_t, Long64_t &ent, Long64_t &bytes, Long64_t &calls)
                           { ent = GetEntriesProcessed(); bytes = GetBytesRead(); calls = GetReadCalls(); return 0; }
   virtual Float_t         GetCurrentRate(Bool_t &all) { all = kTRUE; return (fProgressStatus? fProgressStatus->GetCurrentRate() : 0.); }
   Long64_t                GetTotalEntries() const { return fTotalEntries; }
   virtual TDSetElement   *GetNextPacket(TSlave *sl, TMessage *r);
   virtual void            SetInitTime();
   virtual void            StopProcess(Bool_t abort, Bool_t stoptimer = kFALSE);
   TList                  *GetFailedPackets() { return fFailedPackets; }
   void                    SetFailedPackets(TList *list) { fFailedPackets = list; }
   virtual Int_t           AddWorkers(TList *workers);

   Long64_t      GetBytesRead() const { return (fProgressStatus? fProgressStatus->GetBytesRead() : 0); }
   Long64_t      GetReadCalls() const { return (fProgressStatus? fProgressStatus->GetReadCalls() : 0); }
   Double_t      GetCumProcTime() const { return fProgressStatus->GetProcTime(); }
   Float_t       GetInitTime() const { return fInitTime; }
   Float_t       GetProcTime() const { return fProcTime; }
   TNtuple      *GetProgressPerf(Bool_t steal = kFALSE) { if (steal) { TNtuple *n = fProgressPerf; fProgressPerf = 0; return n;
                                                                } else { return fProgressPerf;} }
   TList        *GetConfigParams(Bool_t steal = kFALSE) { if (steal) { TList *l = fConfigParams; fConfigParams = 0; return l;
                                                                } else { return fConfigParams;} }
   virtual void  MarkBad(TSlave * /*s*/, TProofProgressStatus * /*status*/, TList ** /*missingFiles*/) { return; }
   virtual Int_t AddProcessed(TSlave * /*sl*/, TProofProgressStatus * /*st*/,
                    Double_t /*lat*/, TList ** /*missingFiles*/) { return 0; }
   TProofProgressStatus *GetStatus() { return fProgressStatus; }
   void          SetProgressStatus(TProofProgressStatus *st) { fProgressStatus = st; }
   void          SetTotalEntries(Long64_t ent) { fTotalEntries = ent; }

   TMap         *GetSlaveStats() const { return fSlaveStats; }

   virtual Int_t GetActiveWorkers() { return -1; }

   ClassDef(TVirtualPacketizer,0)  //Generate work packets for parallel processing
};

//------------------------------------------------------------------------------

class TVirtualPacketizer::TVirtualSlaveStat : public TObject {

friend class TPacketizerAdaptive;
friend class TPacketizer;

protected:
   TString        fWrkFQDN;       // Worker FQDN
   TSlave        *fSlave;         // corresponding TSlave record
   TProofProgressStatus *fStatus; // status as of the last finished packet

public:
   const char *GetName() const { return fWrkFQDN.Data(); }
   const char *GetOrdinal() const { return fSlave->GetOrdinal(); }
   Long64_t    GetEntriesProcessed() const { return fStatus?fStatus->GetEntries():-1; }
   Double_t    GetProcTime() const { return fStatus?fStatus->GetProcTime():-1; }
   Float_t     GetAvgRate() { return fStatus->GetRate(); }
   TProofProgressStatus *GetProgressStatus() { return fStatus; }
   virtual TProofProgressStatus *AddProcessed(TProofProgressStatus *st) = 0;
};

#endif
