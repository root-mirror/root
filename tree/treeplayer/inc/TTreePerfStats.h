// @(#)root/treeplayer:$Id$
// Author: Rene Brun 29/10/09

/*************************************************************************
 * Copyright (C) 1995-2009, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TTreePerfStats
#define ROOT_TTreePerfStats

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TTreePerfStats                                                       //
//                                                                      //
// TTree I/O performance measurement                                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include "TVirtualPerfStats.h"
#include "TString.h"
#include <vector>
#include <unordered_map>

class TBrowser;
class TFile;
class TTree;
class TStopwatch;
class TPaveText;
class TGraphErrors;
class TGaxis;
class TText;

class TTreePerfStats : public TVirtualPerfStats {

public:
   struct BasketInfo {
      UInt_t fUsed = {0};       // Number of times the basket was requested from the disk.
      UInt_t fLoaded = {0};     // Number of times the basket was put in the primary TTreeCache
      UInt_t fLoadedMiss = {0}; // Number of times the basket was put in the secondary cache
      UInt_t fMissed = {0};     // Number of times the basket was read directly from the file.
   };

   using BasketList_t = std::vector<std::pair<TBranch*, std::vector<size_t>>>;

protected:
   Int_t         fTreeCacheSize; //TTreeCache buffer size
   Int_t         fNleaves;       //Number of leaves in the tree
   Int_t         fReadCalls;     //Number of read calls
   Int_t         fReadaheadSize; //Readahead cache size
   Long64_t      fBytesRead;     //Number of bytes read
   Long64_t      fBytesReadExtra;//Number of bytes (overhead) of the readahead cache
   Double_t      fRealNorm;      //Real time scale factor for fGraphTime
   Double_t      fRealTime;      //Real time
   Double_t      fCpuTime;       //Cpu time
   Double_t      fDiskTime;      //Time spent in pure raw disk IO
   Double_t      fUnzipTime;     //Time spent uncompressing the data.
   Long64_t      fUnzipInputSize;//Compressed bytes seen by the decompressor.
   Long64_t      fUnzipObjSize;  //Uncompressed bytes produced by the decompressor.
   Double_t      fCompress;      //Tree compression factor
   TString       fName;          //name of this TTreePerfStats
   TString       fHostInfo;      //name of the host system, ROOT version and date
   TFile        *fFile;          //!pointer to the file containing the Tree
   TTree        *fTree;          //!pointer to the Tree being monitored
   TGraphErrors *fGraphIO ;      //pointer to the graph with IO data
   TGraphErrors *fGraphTime ;    //pointer to the graph with timestamp info
   TPaveText    *fPave;          //pointer to annotation pavetext
   TStopwatch   *fWatch;         //TStopwatch pointer
   TGaxis       *fRealTimeAxis;  //pointer to TGaxis object showing real-time
   TText        *fHostInfoText;  //Graphics Text object with the fHostInfo data

   std::unordered_map<TBranch*, size_t>  fBranchIndexCache; // Cache the index of the branch in the cache's array.
   std::vector<std::vector<BasketInfo> > fBasketsInfo;      // Details on which baskets was used, cached, 'miss-cached' or read uncached.Browse

   BasketInfo &GetBasketInfo(TBranch *b, size_t basketNumber);
   BasketInfo &GetBasketInfo(size_t bi, size_t basketNumber);

   virtual void SetFile(TFile *newfile) {
      fFile = newfile;
   }

public:
   TTreePerfStats();
   TTreePerfStats(const char *name, TTree *T);
   virtual ~TTreePerfStats();
   virtual void     Browse(TBrowser *b);
   virtual Int_t    DistancetoPrimitive(Int_t px, Int_t py);
   virtual void     Draw(Option_t *option="");
   virtual void     ExecuteEvent(Int_t event, Int_t px, Int_t py);
   virtual void     Finish();
   virtual Long64_t GetBytesRead() const {return fBytesRead;}
   virtual Long64_t GetBytesReadExtra() const {return fBytesReadExtra;}
   virtual Double_t GetCpuTime()   const {return fCpuTime;}
   virtual Double_t GetDiskTime()  const {return fDiskTime;}
   TGraphErrors    *GetGraphIO()     {return fGraphIO;}
   TGraphErrors    *GetGraphTime()   {return fGraphTime;}
   const char      *GetHostInfo() const{return fHostInfo.Data();}
   const char      *GetName()    const{return fName.Data();}
   virtual Int_t    GetNleaves() const {return fNleaves;}
   virtual Long64_t GetNumEvents() const {return 0;}
   TPaveText       *GetPave()      {return fPave;}
   virtual Int_t    GetReadaheadSize() const {return fReadaheadSize;}
   virtual Int_t    GetReadCalls() const {return fReadCalls;}
   virtual Double_t GetRealTime()  const {return fRealTime;}
   TStopwatch      *GetStopwatch() const {return fWatch;}
   virtual Int_t    GetTreeCacheSize() const {return fTreeCacheSize;}
   virtual Double_t GetUnzipTime() const {return fUnzipTime; }
   virtual void     Paint(Option_t *chopt="");
   virtual void     Print(Option_t *option="") const;

   virtual void     SimpleEvent(EEventType) {}
   virtual void     PacketEvent(const char *, const char *, const char *,
                            Long64_t , Double_t ,Double_t , Double_t ,Long64_t ) {}
   virtual void     FileEvent(const char *, const char *, const char *, const char *, Bool_t) {}
   virtual void     FileOpenEvent(TFile *, const char *, Double_t) {}
   virtual void     FileReadEvent(TFile *file, Int_t len, Double_t start);
   virtual void     UnzipEvent(TObject *tree, Long64_t pos, Double_t start, Int_t complen, Int_t objlen);
   virtual void     RateEvent(Double_t , Double_t , Long64_t , Long64_t) {}

   virtual void     SaveAs(const char *filename="",Option_t *option="") const;
   virtual void     SavePrimitive(std::ostream &out, Option_t *option = "");
   virtual void     SetBytesRead(Long64_t nbytes) {fBytesRead = nbytes;}
   virtual void     SetBytesReadExtra(Long64_t nbytes) {fBytesReadExtra = nbytes;}
   virtual void     SetCompress(Double_t cx) {fCompress = cx;}
   virtual void     SetDiskTime(Double_t t) {fDiskTime = t;}
   virtual void     SetNumEvents(Long64_t) {}
   virtual void     SetCpuTime(Double_t cptime) {fCpuTime = cptime;}
   virtual void     SetGraphIO(TGraphErrors *gr) {fGraphIO = gr;}
   virtual void     SetGraphTime(TGraphErrors *gr) {fGraphTime = gr;}
   virtual void     SetHostInfo(const char *info) {fHostInfo = info;}
   virtual void     SetName(const char *name) {fName = name;}
   virtual void     SetNleaves(Int_t nleaves) {fNleaves = nleaves;}
   virtual void     SetReadaheadSize(Int_t nbytes) {fReadaheadSize = nbytes;}
   virtual void     SetReadCalls(Int_t ncalls) {fReadCalls = ncalls;}
   virtual void     SetRealNorm(Double_t rnorm) {fRealNorm = rnorm;}
   virtual void     SetRealTime(Double_t rtime) {fRealTime = rtime;}
   virtual void     SetTreeCacheSize(Int_t nbytes) {fTreeCacheSize = nbytes;}
   virtual void     SetUnzipTime(Double_t uztime) {fUnzipTime = uztime;}

   virtual void     PrintBasketInfo(Option_t *option = "") const;
   virtual void     SetLoaded(TBranch *b, size_t basketNumber) { ++GetBasketInfo(b, basketNumber).fLoaded; }
   virtual void     SetLoaded(size_t bi, size_t basketNumber) { ++GetBasketInfo(bi, basketNumber).fLoaded; }
   virtual void     SetLoadedMiss(TBranch *b, size_t basketNumber) { ++GetBasketInfo(b, basketNumber).fLoadedMiss; }
   virtual void     SetLoadedMiss(size_t bi, size_t basketNumber) { ++GetBasketInfo(bi, basketNumber).fLoadedMiss; }
   virtual void     SetMissed(TBranch *b, size_t basketNumber) { ++GetBasketInfo(b, basketNumber).fMissed; }
   virtual void     SetMissed(size_t bi, size_t basketNumber) { ++GetBasketInfo(bi, basketNumber).fMissed; }
   virtual void     SetUsed(TBranch *b, size_t basketNumber) { ++GetBasketInfo(b, basketNumber).fUsed; }
   virtual void     SetUsed(size_t bi, size_t basketNumber) { ++GetBasketInfo(bi, basketNumber).fUsed; }
   virtual void     UpdateBranchIndices(TObjArray *branchNames);

   BasketList_t     GetDuplicateBasketCache() const;

   ClassDef(TTreePerfStats, 8) // TTree I/O performance measurement
};

#endif
