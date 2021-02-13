// @(#)root/proofplayer:$Id$
// Author: Kristjan Gulbrandsen   11/05/04

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TPerfStats
#define ROOT_TPerfStats

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TPerfStats                                                           //
//                                                                      //
// Provides the interface for the PROOF internal performance measurment //
// and event tracing.                                                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include "TObject.h"
#include "TObjArray.h"
#include "TTimeStamp.h"
#include "TString.h"
#include "TVirtualPerfStats.h"


class TDSet;
class TH1D;
class TH1I;
class TH2D;
class TList;
class TTree;

class TPerfEvent : public TObject {

public:
   TString                       fEvtNode;   // node on which the event was generated
   TTimeStamp                    fTimeStamp; // time offset from start of run
   TVirtualPerfStats::EEventType fType;
   TString                       fSlaveName;
   TString                       fNodeName;
   TString                       fFileName;
   TString                       fFileClass;
   TString                       fSlave;
   Long64_t                      fEventsProcessed;
   Long64_t                      fBytesRead;
   Long64_t                      fLen;
   Double_t                      fLatency;
   Double_t                      fProcTime;
   Double_t                      fCpuTime;
   Bool_t                        fIsStart;
   Bool_t                        fIsOk;

   TPerfEvent(TTimeStamp *offset = 0);
   virtual ~TPerfEvent() {}

   Bool_t   IsSortable() const { return kTRUE; }
   Int_t    Compare(const TObject *obj) const;
   void     Print(Option_t *option="") const;

   ClassDef(TPerfEvent,3) // Class holding TProof Event Info
};


class TPerfStats : public TVirtualPerfStats {

friend class TProofMonSender;

private:
   TTree         *fTrace;        //!TTree with trace events
   TTimeStamp     fTzero;        //!start time of this run
   TPerfEvent    *fPerfEvent;    //!TPerfEvent used to fill tree
   TH1D          *fPacketsHist;  //!histogram of packets processed per slave
   TH1I          *fProcPcktHist; //!histogram of packets being processed per slave
   TH1D          *fEventsHist;   //!histogram of events processed per slave
   TH1D          *fNodeHist;     //!histogram of slaves per file serving node
   TH2D          *fLatencyHist;  //!histogram of latency due to packet requests
   TH2D          *fProcTimeHist; //!histogram of real time spent processing packets
   TH2D          *fCpuTimeHist;  //!histogram of cpu time spent processing packets
   Long64_t       fBytesRead;    //!track bytes read of main file
   Double_t       fTotCpuTime;   //!total cpu time of all slaves
   Long64_t       fTotBytesRead; //!total bytes read on all slaves
   Long64_t       fTotEvents;    //!total number of events processed
   Long64_t       fNumEvents;    //!total number of events to be processed
   Int_t          fSlaves;       //!number of active slaves

   Bool_t         fDoHist;       //!Fill histos
   Bool_t         fDoTrace;      //!Trace details in master
   Bool_t         fDoTraceRate;  //!Trace processing rate in master
   Bool_t         fDoSlaveTrace; //!Full tracing in workers
   Bool_t         fDoQuota;      //!Save stats on SQL server for quota management

   Bool_t         fMonitorPerPacket; //!Whether to send the full entry per each packet

   TObjArray      fMonSenders;   //!Monitoring engines

   TString        fDataSet;      //!Dataset string
   Int_t          fDataSetLen;   //!Maximum size of the dataset string fDataSet
   Int_t          fDataSetSize;  //!# of files in the dataset
   TDSet         *fDSet;         //!Saved pointer to the TDSet object
   TList         *fOutput;       //!Saved pointer to the output list

   static Long_t  fgVirtMemMax;   //! Max virtual memory used by this process
   static Long_t  fgResMemMax;    //! Max resident memory used by this process

   TPerfStats(TList *input, TList *output);
   void WriteQueryLog();

   void SetFile(TFile *) {}

public:
   virtual ~TPerfStats();

   void SimpleEvent(EEventType type);
   void PacketEvent(const char *slave, const char *slavename, const char *filename,
                    Long64_t eventsprocessed, Double_t latency,
                    Double_t proctime, Double_t cputime, Long64_t bytesRead);
   void FileEvent(const char *slave, const char *slavename, const char *nodename, const char *filename,
                  Bool_t isStart);

   void FileOpenEvent(TFile *file, const char *filename, Double_t start);
   void FileReadEvent(TFile *file, Int_t len, Double_t start);
   void UnzipEvent(TObject *tree, Long64_t pos, Double_t start, Int_t complen, Int_t objlen);
   void RateEvent(Double_t proctime, Double_t deltatime,
                  Long64_t eventsprocessed, Long64_t bytesRead);
   void SetBytesRead(Long64_t num);
   Long64_t GetBytesRead() const;
   void SetNumEvents(Long64_t num) { fNumEvents = num; }
   Long64_t GetNumEvents() const { return fNumEvents; }

   void        PrintBasketInfo(Option_t * = "") const {}
   void        SetLoaded(TBranch *, size_t) {}
   void        SetLoaded(size_t, size_t) {}
   void        SetLoadedMiss(TBranch *, size_t) {}
   void        SetLoadedMiss(size_t, size_t) {}
   void        SetMissed(TBranch *, size_t) {}
   void        SetMissed(size_t, size_t) {}
   void        SetUsed(TBranch *, size_t) {}
   void        SetUsed(size_t, size_t) {}
   void        UpdateBranchIndices(TObjArray *) {}

   static void Start(TList *input, TList *output);
   static void Stop();
   static void Setup(TList *input);
   static void SetMemValues();
   static void GetMemValues(Long_t &vmax, Long_t &rmax);

   ClassDef(TPerfStats,0)  // Class for collecting PROOF statistics
};


#endif
