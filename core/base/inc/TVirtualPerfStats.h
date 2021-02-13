// @(#)root/base:$Id$
// Author: Kristjan Gulbrandsen   11/05/04

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TVirtualPerfStats
#define ROOT_TVirtualPerfStats

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TVirtualPerfStats                                                    //
//                                                                      //
// Provides the interface for the PROOF internal performance measurment //
// and event tracing.                                                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include "TObject.h"


class TFile;
class TBranch;
class TChain;

class TVirtualPerfStats : public TObject {

private:
   friend class TChain;

   virtual void SetFile(TFile *) = 0;

public:
   virtual ~TVirtualPerfStats() {}

   enum EEventType {
      kUnDefined,
      kPacket,       //info of single packet processing
      kStart,        //begin of run
      kStop,         //end of run
      kFile,         //file started/finished in packetizer
      kFileOpen,     //opening data file statistics
      kFileRead,     //reading data file event
      kRate,         //processing {evt, MB} rates
      kNumEventType  //number of entries, must be last
   };

   static TVirtualPerfStats *&CurrentPerfStats();  // Return the current perfStats for this thread.

   virtual void SimpleEvent(EEventType type) = 0;

   virtual void PacketEvent(const char *slave, const char *slavename, const char *filename,
                            Long64_t eventsprocessed, Double_t latency,
                            Double_t proctime, Double_t cputime,
                            Long64_t bytesRead) = 0;

   virtual void FileEvent(const char *slave, const char *slavename, const char *nodename,
                          const char *filename, Bool_t isStart) = 0;

   virtual void FileOpenEvent(TFile *file, const char *filename, Double_t start) = 0;

   virtual void FileReadEvent(TFile *file, Int_t len, Double_t start) = 0;

   virtual void UnzipEvent(TObject *tree, Long64_t pos, Double_t start, Int_t complen, Int_t objlen) = 0;

   virtual void RateEvent(Double_t proctime, Double_t deltatime,
                          Long64_t eventsprocessed, Long64_t bytesRead) = 0;

   virtual void SetBytesRead(Long64_t num) = 0;
   virtual Long64_t GetBytesRead() const = 0;
   virtual void SetNumEvents(Long64_t num) = 0;
   virtual Long64_t GetNumEvents() const = 0;

   virtual void PrintBasketInfo(Option_t *option = "") const = 0;
   virtual void SetLoaded(TBranch *b, size_t basketNumber) = 0;
   virtual void SetLoaded(size_t bi, size_t basketNumber) = 0;
   virtual void SetLoadedMiss(TBranch *b, size_t basketNumber) = 0;
   virtual void SetLoadedMiss(size_t bi, size_t basketNumber) = 0;
   virtual void SetMissed(TBranch *b, size_t basketNumber) = 0;
   virtual void SetMissed(size_t bi, size_t basketNumber) = 0;
   virtual void SetUsed(TBranch *b, size_t basketNumber) = 0;
   virtual void SetUsed(size_t bi, size_t basketNumber) = 0;
   virtual void UpdateBranchIndices(TObjArray *branches) = 0;

   static const char *EventType(EEventType type);

   ClassDef(TVirtualPerfStats,0)  // ABC for collecting PROOF statistics
};


#ifndef __CINT__
#define gPerfStats (TVirtualPerfStats::CurrentPerfStats())

#elif defined(__MAKECINT__)
// To properly handle the use of gPerfStats in header files (in static declarations)
R__EXTERN TVirtualPerfStats *gPerfStats;
#endif


#endif
