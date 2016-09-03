// @(#)root/tree:$Id$
// Author: Rene Brun   04/06/2006

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TTreeCache
\ingroup tree

A specialized TFileCacheRead object for a TTree.

This class acts as a file cache, registering automatically the baskets from
the branches being processed (TTree::Draw or TTree::Process and TSelectors)
when in the learning phase. The learning phase is by default 100 entries.
It can be changed via TTreeCache::SetLearnEntries.

This cache speeds-up considerably the performance, in particular
when the Tree is accessed remotely via a high latency network.

The default cache size (10 Mbytes) may be changed via the function
TTree::SetCacheSize

Only the baskets for the requested entry range are put in the cache

For each Tree being processed a TTreeCache object is created.
This object is automatically deleted when the Tree is deleted or
when the file is deleted.

- Special case of a TChain
  Once the training is done on the first Tree, the list of branches
  in the cache is kept for the following files.

- Special case of a TEventlist
  if the Tree or TChain has a TEventlist, only the buffers
  referenced by the list are put in the cache.

The learning period is started or restarted when:
   - TTree automatically creates a cache. This feature can be
     controlled with an env. variable or the TTreeCache.Size option.
   - TTree::SetCacheSize is called with a non-zero size and a cache
     did not previously exist
   - TTreeCache::StartLearningPhase is called.
   - TTreeCache::SetEntryRange is called
        * and the learning is not yet finished
        * and has not been set to manual
        * and the new minimun entry is different.

The learning period is stopped (and prefetching is actually started) when:
   - TTreeCache::StopLearningPhase is called.
   - An entry outside the 'learning' range is requested
     The 'learning range is from fEntryMin (default to 0) to
     fEntryMin + fgLearnEntries (default to 100).
   - A 'cached' TChain switches over to a new file.

Further, the TreeCache can optimize its behavior on a cache miss.  When
miss optimization is enabled, it will track all branches utilized after
the learning phase (those that cause a cache miss).  When one cache miss
occurs, then all the utilized branches will be prefetched for that event.
This optimization utilizes the observation that infrequently accessed
branches are often accessed together.  For example, this will greatly speed
up an analysis where the results of a trigger are read out for every branch,
but the majority of event collections are read only when the trigger results
pass a set of filters.  NOTE - when this mode is enabled, the memory dedicated
to the cache will up to double in the case of cache miss.

## WHY DO WE NEED the TreeCache when doing data analysis?

When writing a TTree, the branch buffers are kept in memory.
A typical branch buffersize (before compression) is typically 32 KBytes.
After compression, the zipped buffer may be just a few Kbytes.
The branch buffers cannot be much larger in case of Trees with several
hundred or thousand branches.

When writing, this does not generate a performance problem because branch
buffers are always written sequentially and the OS is in general clever enough
to flush the data to the output file when a few MBytes of data have to be written.
When reading at the contrary, one may hit a performance problem when reading
across a network (LAN or WAN) and the network latency is high.
For example in a WAN with 10ms latency, reading 1000 buffers of 10 KBytes each
with no cache will imply 10s penalty where a local read of the 10 MBytes would
take about 1 second.

The TreeCache will try to prefetch all the buffers for the selected branches
such that instead of transferring 1000 buffers of 10 Kbytes, it will be able
to transfer one single large buffer of 10 Mbytes in one single transaction.
Not only the TreeCache minimizes the number of transfers, but in addition
it can sort the blocks to be read in increasing order such that the file
is read sequentially.

Systems like xrootd, dCache or httpd take advantage of the TreeCache in
reading ahead as much data as they can and return to the application
the maximum data specified in the cache and have the next chunk of data ready
when the next request comes.

## HOW TO USE the TreeCache

A few use cases are discussed below. A cache may be created with automatic sizing
when a TTree is used:

Caches are created and automatically sized for TTrees when TTreeCache.Size or
the environment variable ROOT_TTREECACHE_SIZE is set to a sizing factor.

But there are many possible configurations where manual control may be wanted.
In some applications you know a priori the list of branches to read. In other
applications the analysis loop calls several layers of user functions where it
is impossible to predict a priori which branches will be used. This
is probably the most frequent case. In this case ROOT I/O will flag used
branches automatically when a branch buffer is read during the learning phase.
The TreeCache interface provides functions to instruct the cache about the used
branches if they are known a priori. In the examples below, portions of analysis
code are shown. The few statements involving the TreeCache are marked with `//<<<`

### 1. with TTree::Draw

the TreeCache is automatically used by TTree::Draw. The function knows
which branches are used in the query and it puts automatically these branches
in the cache. The entry range is also known automatically.

### 2. with TTree::Process and TSelectors

You must enable the cache and tell the system which branches to cache
and also specify the entry range. It is important to specify the entry range
in case you process only a subset of the events, otherwise you run the risk
to store in the cache entries that you do not need.

#### example 2a
~~~ {.cpp}
    TTree *T = (TTree*)f->Get("mytree");
    Long64_t nentries = T->GetEntries();
    Int_t cachesize = 10000000; //10 MBytes
    T->SetCacheSize(cachesize); //<<<
    T->AddBranchToCache("*",kTRUE);    //<<< add all branches to the cache
    T->Process('myselector.C+");
    //in the TSelector::Process function we read all branches
    T->GetEntry(i);
    ... here you process your entry
~~~
#### example 2b

in the Process function we read a subset of the branches.
Only the branches used in the first entry will be put in the cache
~~~ {.cpp}
    TTree *T = (TTree*)f->Get("mytree");
    //we want to process only the 200 first entries
    Long64_t nentries=200;
    int efirst= 0;
    int elast = efirst+nentries;
    Int_t cachesize = 10000000; //10 MBytes
    TTreeCache::SetLearnEntries(1);  //<<< we can take the decision after 1 entry
    T->SetCacheSize(cachesize);      //<<<
    T->SetCacheEntryRange(efirst,elast); //<<<
    T->Process('myselector.C+","",nentries,efirst);
    // in the TSelector::Process we read only 2 branches
    TBranch *b1 = T->GetBranch("branch1");
    b1->GetEntry(i);
    if (somecondition) return;
    TBranch *b2 = T->GetBranch("branch2");
    b2->GetEntry(i);
    ... here you process your entry
~~~
### 3. with your own event loop

#### example 3a

in your analysis loop, you always use 2 branches. You want to prefetch
the branch buffers for these 2 branches only.
~~~ {.cpp}
    TTree *T = (TTree*)f->Get("mytree");
    TBranch *b1 = T->GetBranch("branch1");
    TBranch *b2 = T->GetBranch("branch2");
    Long64_t nentries = T->GetEntries();
    Int_t cachesize = 10000000; //10 MBytes
    T->SetCacheSize(cachesize);     //<<<
    T->AddBranchToCache(b1,kTRUE);  //<<<add branch1 and branch2 to the cache
    T->AddBranchToCache(b2,kTRUE);  //<<<
    T->StopCacheLearningPhase();    //<<<
    for (Long64_t i=0;i<nentries;i++) {
       T->LoadTree(i); //<<< important call when calling TBranch::GetEntry after
       b1->GetEntry(i);
       if (some condition not met) continue;
       b2->GetEntry(i);
       if (some condition not met) continue;
       //here we read the full event only in some rare cases.
       //there is no point in caching the other branches as it might be
       //more economical to read only the branch buffers really used.
       T->GetEntry(i);
       .. process the rare but interesting cases.
       ... here you process your entry
    }
~~~
#### example 3b

in your analysis loop, you always use 2 branches in the main loop.
you also call some analysis functions where a few more branches will be read.
but you do not know a priori which ones. There is no point in prefetching
branches that will be used very rarely.
~~~ {.cpp}
    TTree *T = (TTree*)f->Get("mytree");
    Long64_t nentries = T->GetEntries();
    Int_t cachesize = 10000000;   //10 MBytes
    T->SetCacheSize(cachesize);   //<<<
    T->SetCacheLearnEntries(5);   //<<< we can take the decision after 5 entries
    TBranch *b1 = T->GetBranch("branch1");
    TBranch *b2 = T->GetBranch("branch2");
    for (Long64_t i=0;i<nentries;i++) {
       T->LoadTree(i);
       b1->GetEntry(i);
       if (some condition not met) continue;
       b2->GetEntry(i);
       //at this point we may call a user function where a few more branches
       //will be read conditionally. These branches will be put in the cache
       //if they have been used in the first 10 entries
       if (some condition not met) continue;
       //here we read the full event only in some rare cases.
       //there is no point in caching the other branches as it might be
       //more economical to read only the branch buffers really used.
       T->GetEntry(i);
       .. process the rare but interesting cases.
       ... here you process your entry
    }
~~~
## SPECIAL CASES WHERE TreeCache should not be activated

When reading only a small fraction of all entries such that not all branch
buffers are read, it might be faster to run without a cache.

## HOW TO VERIFY That the TreeCache has been used and check its performance

Once your analysis loop has terminated, you can access/print the number
of effective system reads for a given file with a code like
(where TFile* f is a pointer to your file)
~~~ {.cpp}
    printf("Reading %lld bytes in %d transactions\n",f->GetBytesRead(),  f->GetReadCalls());
~~~
*/

#include "TSystem.h"
#include "TEnv.h"
#include "TTreeCache.h"
#include "TChain.h"
#include "TList.h"
#include "TBranch.h"
#include "TEventList.h"
#include "TObjString.h"
#include "TRegexp.h"
#include "TLeaf.h"
#include "TFriendElement.h"
#include "TFile.h"
#include "TMath.h"
#include <limits.h>

Int_t TTreeCache::fgLearnEntries = 100;

ClassImp(TTreeCache);

////////////////////////////////////////////////////////////////////////////////
/// Default Constructor.

TTreeCache::TTreeCache() : TFileCacheRead(),
   fEntryMin(0),
   fEntryMax(1),
   fEntryCurrent(-1),
   fEntryNext(-1),
   fNbranches(0),
   fNReadOk(0),
   fNReadMiss(0),
   fNReadPref(0),
   fBranches(0),
   fBrNames(0),
   fTree(0),
   fIsLearning(kTRUE),
   fIsManual(kFALSE),
   fFirstBuffer(kTRUE),
   fOneTime(kFALSE),
   fReverseRead(0),
   fFillTimes(0),
   fFirstTime(kTRUE),
   fFirstEntry(-1),
   fReadDirectionSet(kFALSE),
   fEnabled(kTRUE),
   fPrefillType(GetConfiguredPrefillType()),
   fAutoCreated(kFALSE)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TTreeCache::TTreeCache(TTree *tree, Int_t buffersize) : TFileCacheRead(tree->GetCurrentFile(),buffersize,tree),
   fEntryMin(0),
   fEntryMax(tree->GetEntriesFast()),
   fEntryCurrent(-1),
   fEntryNext(0),
   fNbranches(0),
   fNReadOk(0),
   fNReadMiss(0),
   fNReadPref(0),
   fBranches(0),
   fBrNames(new TList),
   fTree(tree),
   fIsLearning(kTRUE),
   fIsManual(kFALSE),
   fFirstBuffer(kTRUE),
   fOneTime(kFALSE),
   fReverseRead(0),
   fFillTimes(0),
   fFirstTime(kTRUE),
   fFirstEntry(-1),
   fReadDirectionSet(kFALSE),
   fEnabled(kTRUE),
   fPrefillType(GetConfiguredPrefillType()),
   fAutoCreated(kFALSE)
{
   fEntryNext = fEntryMin + fgLearnEntries;
   Int_t nleaves = tree->GetListOfLeaves()->GetEntries();
   fBranches = new TObjArray(nleaves);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor. (in general called by the TFile destructor)

TTreeCache::~TTreeCache()
{
   // Informe the TFile that we have been deleted (in case
   // we are deleted explicitly by legacy user code).
   if (fFile) fFile->SetCacheRead(0, fTree);

   delete fBranches;
   if (fBrNames) {fBrNames->Delete(); delete fBrNames; fBrNames=0;}
}

////////////////////////////////////////////////////////////////////////////////
/// Add a branch to the list of branches to be stored in the cache
/// this function is called by TBranch::GetBasket
/// Returns:
///  - 0 branch added or already included
///  - -1 on error

Int_t TTreeCache::AddBranch(TBranch *b, Bool_t subbranches /*= kFALSE*/)
{
   if (!fIsLearning) {
      return -1;
   }

   // Reject branch that are not from the cached tree.
   if (!b || fTree->GetTree() != b->GetTree()) return -1;

   // Is this the first addition of a branch (and we are learning and we are in
   // the expected TTree), then prefill the cache.  (We expect that in future
   // release the Prefill-ing will be the default so we test for that inside the
   // LearnPrefill call).
   if (fNbranches == 0 && fEntryMin >= 0 && b->GetReadEntry() == fEntryMin) LearnPrefill();

   //Is branch already in the cache?
   Bool_t isNew = kTRUE;
   for (int i=0;i<fNbranches;i++) {
      if (fBranches->UncheckedAt(i) == b) {isNew = kFALSE; break;}
   }
   if (isNew) {
      fTree = b->GetTree();
      fBranches->AddAtAndExpand(b, fNbranches);
      fBrNames->Add(new TObjString(b->GetName()));
      fNbranches++;
      if (gDebug > 0) printf("Entry: %lld, registering branch: %s\n",b->GetTree()->GetReadEntry(),b->GetName());
   }

   // process subbranches
   Int_t res = 0;
   if (subbranches) {
      TObjArray *lb = b->GetListOfBranches();
      Int_t nb = lb->GetEntriesFast();
      for (Int_t j = 0; j < nb; j++) {
         TBranch* branch = (TBranch*) lb->UncheckedAt(j);
         if (!branch) continue;
         if (AddBranch(branch, subbranches)<0) {
            res = -1;
         }
      }
   }
   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Add a branch to the list of branches to be stored in the cache
/// this is to be used by user (thats why we pass the name of the branch).
/// It works in exactly the same way as TTree::SetBranchStatus so you
/// probably want to look over there for details about the use of bname
/// with regular expressions.
/// The branches are taken with respect to the Owner of this TTreeCache
/// (i.e. the original Tree)
/// NB: if bname="*" all branches are put in the cache and the learning phase stopped
/// Returns:
///  - 0 branch added or already included
///  - -1 on error

Int_t TTreeCache::AddBranch(const char *bname, Bool_t subbranches /*= kFALSE*/)
{
   TBranch *branch, *bcount;
   TLeaf *leaf, *leafcount;

   Int_t i;
   Int_t nleaves = (fTree->GetListOfLeaves())->GetEntriesFast();
   TRegexp re(bname,kTRUE);
   Int_t nb = 0;
   Int_t res = 0;

   // first pass, loop on all branches
   // for leafcount branches activate/deactivate in function of status
   Bool_t all = kFALSE;
   if (!strcmp(bname,"*")) all = kTRUE;
   for (i=0;i<nleaves;i++)  {
      leaf = (TLeaf*)(fTree->GetListOfLeaves())->UncheckedAt(i);
      branch = (TBranch*)leaf->GetBranch();
      TString s = branch->GetName();
      if (!all) { //Regexp gives wrong result for [] in name
         TString longname;
         longname.Form("%s.%s",fTree->GetName(),branch->GetName());
         if (strcmp(bname,branch->GetName())
             && longname != bname
             && s.Index(re) == kNPOS) continue;
      }
      nb++;
      if (AddBranch(branch, subbranches)<0) {
         res = -1;
      }
      leafcount = leaf->GetLeafCount();
      if (leafcount && !all) {
         bcount = leafcount->GetBranch();
         if (AddBranch(bcount, subbranches)<0) {
            res = -1;
         }
      }
   }
   if (nb==0 && strchr(bname,'*')==0) {
      branch = fTree->GetBranch(bname);
      if (branch) {
         if (AddBranch(branch, subbranches)<0) {
            res = -1;
         }
         ++nb;
      }
   }

   //search in list of friends
   UInt_t foundInFriend = 0;
   if (fTree->GetListOfFriends()) {
      TIter nextf(fTree->GetListOfFriends());
      TFriendElement *fe;
      TString name;
      while ((fe = (TFriendElement*)nextf())) {
         TTree *t = fe->GetTree();
         if (t==0) continue;

         // If the alias is present replace it with the real name.
         char *subbranch = (char*)strstr(bname,fe->GetName());
         if (subbranch!=bname) subbranch = 0;
         if (subbranch) {
            subbranch += strlen(fe->GetName());
            if ( *subbranch != '.' ) subbranch = 0;
            else subbranch ++;
         }
         if (subbranch) {
            name.Form("%s.%s",t->GetName(),subbranch);
            if (AddBranch(name, subbranches)<0) {
               res = -1;
            }
            ++foundInFriend;
         }
      }
   }
   if (!nb && !foundInFriend) {
      if (gDebug > 0) printf("AddBranch: unknown branch -> %s \n", bname);
      Error("AddBranch", "unknown branch -> %s", bname);
      return -1;
   }
   //if all branches are selected stop the learning phase
   if (*bname == '*') {
      fEntryNext = -1; // We are likely to have change the set of branches, so for the [re-]reading of the cluster.
      StopLearningPhase();
   }
   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove a branch to the list of branches to be stored in the cache
/// this function is called by TBranch::GetBasket.
/// Returns:
///  - 0 branch dropped or not in cache
///  - -1 on error

Int_t TTreeCache::DropBranch(TBranch *b, Bool_t subbranches /*= kFALSE*/)
{
   if (!fIsLearning) {
      return -1;
   }

   // Reject branch that are not from the cached tree.
   if (!b || fTree->GetTree() != b->GetTree()) return -1;

   //Is branch already in the cache?
   if (fBranches->Remove(b)) {
      --fNbranches;
      if (gDebug > 0) printf("Entry: %lld, un-registering branch: %s\n",b->GetTree()->GetReadEntry(),b->GetName());
   }
   delete fBrNames->Remove(fBrNames->FindObject(b->GetName()));

   // process subbranches
   Int_t res = 0;
   if (subbranches) {
      TObjArray *lb = b->GetListOfBranches();
      Int_t nb = lb->GetEntriesFast();
      for (Int_t j = 0; j < nb; j++) {
         TBranch* branch = (TBranch*) lb->UncheckedAt(j);
         if (!branch) continue;
         if (DropBranch(branch, subbranches)<0) {
            res = -1;
         }
      }
   }
   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove a branch to the list of branches to be stored in the cache
/// this is to be used by user (thats why we pass the name of the branch).
/// It works in exactly the same way as TTree::SetBranchStatus so you
/// probably want to look over there for details about the use of bname
/// with regular expressions.
/// The branches are taken with respect to the Owner of this TTreeCache
/// (i.e. the original Tree)
/// NB: if bname="*" all branches are put in the cache and the learning phase stopped
/// Returns:
///  - 0 branch dropped or not in cache
///  - -1 on error

Int_t TTreeCache::DropBranch(const char *bname, Bool_t subbranches /*= kFALSE*/)
{
   TBranch *branch, *bcount;
   TLeaf *leaf, *leafcount;

   Int_t i;
   Int_t nleaves = (fTree->GetListOfLeaves())->GetEntriesFast();
   TRegexp re(bname,kTRUE);
   Int_t nb = 0;
   Int_t res = 0;

   // first pass, loop on all branches
   // for leafcount branches activate/deactivate in function of status
   Bool_t all = kFALSE;
   if (!strcmp(bname,"*")) all = kTRUE;
   for (i=0;i<nleaves;i++)  {
      leaf = (TLeaf*)(fTree->GetListOfLeaves())->UncheckedAt(i);
      branch = (TBranch*)leaf->GetBranch();
      TString s = branch->GetName();
      if (!all) { //Regexp gives wrong result for [] in name
         TString longname;
         longname.Form("%s.%s",fTree->GetName(),branch->GetName());
         if (strcmp(bname,branch->GetName())
             && longname != bname
             && s.Index(re) == kNPOS) continue;
      }
      nb++;
      if (DropBranch(branch, subbranches)<0) {
         res = -1;
      }
      leafcount = leaf->GetLeafCount();
      if (leafcount && !all) {
         bcount = leafcount->GetBranch();
         if (DropBranch(bcount, subbranches)<0) {
            res = -1;
         }
      }
   }
   if (nb==0 && strchr(bname,'*')==0) {
      branch = fTree->GetBranch(bname);
      if (branch) {
         if (DropBranch(branch, subbranches)<0) {
            res = -1;
         }
         ++nb;
      }
   }

   //search in list of friends
   UInt_t foundInFriend = 0;
   if (fTree->GetListOfFriends()) {
      TIter nextf(fTree->GetListOfFriends());
      TFriendElement *fe;
      TString name;
      while ((fe = (TFriendElement*)nextf())) {
         TTree *t = fe->GetTree();
         if (t==0) continue;

         // If the alias is present replace it with the real name.
         char *subbranch = (char*)strstr(bname,fe->GetName());
         if (subbranch!=bname) subbranch = 0;
         if (subbranch) {
            subbranch += strlen(fe->GetName());
            if ( *subbranch != '.' ) subbranch = 0;
            else subbranch ++;
         }
         if (subbranch) {
            name.Form("%s.%s",t->GetName(),subbranch);
            if (DropBranch(name, subbranches)<0) {
               res = -1;
            }
            ++foundInFriend;
         }
      }
   }
   if (!nb && !foundInFriend) {
      if (gDebug > 0) printf("DropBranch: unknown branch -> %s \n", bname);
      Error("DropBranch", "unknown branch -> %s", bname);
      return -1;
   }
   //if all branches are selected stop the learning phase
   if (*bname == '*') {
      fEntryNext = -1; // We are likely to have change the set of branches, so for the [re-]reading of the cluster.
   }
   return res;
}

/**
 * Enable / disable the miss cache.
 */
void TTreeCache::SetOptimizeMisses(bool opt) {

   if (opt && !fMissCache) {
      ResetMissCache();
   }
   fOptimizeMisses = opt;
}

/**
 * Reset all the training
 */
void TTreeCache::ResetMissCache() {

   fLastMiss = -1;
   fFirstMiss = -1;
   // For the most part, we keep the allocated memory around to prevent memory churn.
   if (!fMissBranches) {
      fMissBranches.reset(new std::vector<TBranch*>());
   } else {
      fMissBranches->clear();
   }
   if (!fMissCache) {
      fMissCache.reset(new std::vector<char>());
   } else {
      fMissCache->clear();
   }
   if (!fEntries) {
      fEntries.reset(new std::vector<std::pair<ULong64_t, UInt_t>>());
   } else {
      fEntries->clear();
   }
   if (!fEntryOffsets) {
      fEntryOffsets.reset(new std::vector<size_t>());
   } else {
      fEntryOffsets->clear();
   }
}

std::pair<ULong64_t, UInt_t> TTreeCache::FindBranchBasket(TBranch &b) {
   if (b.GetDirectory() == 0) {
      //printf("Branch at %p has no valid directory.\n", &b);
      return std::make_pair(0, 0);
   }
   if (b.GetDirectory()->GetFile() != fFile) {
      //printf("Branch at %p is in wrong file (branch file %p, my file %p).\n", &b, b.GetDirectory()->GetFile(), fFile);
      return std::make_pair(0, 0);
   }

   //printf("Trying to find a basket for branch %p\n", &b);
   // Pull in metadata about branch; make sure it is valid
   Int_t *lbaskets   = b.GetBasketBytes();
   Long64_t *entries = b.GetBasketEntry();
   if (!lbaskets || !entries) {
      //printf("No baskets or entries.\n");
      return std::make_pair(0, 0);
   }
   //Int_t blistsize = b.GetListOfBaskets()->GetSize();
   Int_t blistsize = b.GetWriteBasket();
   if (blistsize <= 0) {
      //printf("Basket list is size 0.\n");
      return std::make_pair(0, 0);
   }

   // Search for the basket that contains the event of interest.  Unlike the primary cache, we
   // are only interested in a single basket per branch - we don't try to fill the cache.
   Long64_t basketOffset = TMath::BinarySearch(blistsize, entries, fTree->GetReadEntry());
   if (basketOffset < 0) { // No entry found.
      //printf("No entry offset found for entry %ld\n", fTree->GetReadEntry());
      return std::make_pair(0, 0);
   }

   // Check to see if there's already a copy of this basket in memory.  If so, don't fetch it
   if ((basketOffset < blistsize) &&
           b.GetListOfBaskets()->UncheckedAt(basketOffset)) {

       //printf("Basket is already in memory.\n");
       return std::make_pair(0, 0);
   }

   Long64_t pos = b.GetBasketSeek(basketOffset);
   Int_t len = lbaskets[basketOffset];
   if (pos <= 0 || len <= 0) {
      /*printf("Basket returned was invalid (basketOffset=%ld, pos=%ld, len=%d).\n", basketOffset, pos, len);
      for (int idx=0; idx<blistsize; idx++) {
         printf("Basket entry %d, first event %d, pos %ld\n", idx, entries[idx], b.GetBasketSeek(idx));
      }*/
      return std::make_pair(0, 0);
    } // Sanity check
   // Do not cache a basket if it is bigger than the cache size!
   if (len > fBufferSizeMin) {
      //printf("Basket size is greater than the cache size.\n");
      return std::make_pair(0, 0);
   }

   return std::make_pair(pos, len);
}


TBranch *TTreeCache::CalculateMissEntries(Long64_t pos, int len, bool all) {
   if ((pos < 0) || (len < 0)) {  // TODO: unlikely
      return nullptr;
   }
   int count = all ? (fTree->GetListOfLeaves())->GetEntriesFast() : fMissBranches->size();
   fEntries->reserve(count);
   fEntries->clear();
   bool found_request = false;
   TBranch *resultBranch = nullptr;
   //printf("Will search %d branches for basket at %ld.\n", count, pos);
   for (int i=0; i<count; i++) {
      TBranch *b = all ? static_cast<TBranch*>(static_cast<TLeaf*>((fTree->GetListOfLeaves())->UncheckedAt(i))->GetBranch()) : (*fMissBranches)[i];
      std::pair<ULong64_t, UInt_t> iopos = FindBranchBasket(*b);
      if (iopos.second == 0) {  // Error indicator
         continue;
      }
      if (iopos.first == static_cast<ULong64_t>(pos) && iopos.second == static_cast<UInt_t>(len)) {
         found_request = true;
         resultBranch = b;
         // Note that we continue to iterate; fills up the rest of the entries in the cache.
      }
      // At this point, we are ready to push back a new offset
      fEntries->emplace_back(std::move(iopos));
   }
   if (!found_request) {
      // We have gone through all the branches in this file and the requested basket
      // doesn't appear to be in any of them.  Likely a logic error / bug.
      fEntries->clear();
   }
   return resultBranch;
}


/**
 * Process a cache miss; (pos, len) isn't in the buffer.
 *
 * The first time we have a miss, we buffer as many baskets we can (up to the maximum
 * size of the TTreeCache) in memory from all branches that are not in the prefetch list.
 * 
 * Subsequent times, we fetch all the buffers corresponding to branches that had previously
 * seen misses.  If it turns out the (pos, len) isn't in the list of branches, we treat this
 * as if it was the first miss.
 */
bool TTreeCache::ProcessMiss(Long64_t pos, int len) {

   bool firstMiss = false;
   if (fFirstMiss == -1) {
      fFirstMiss = fEntryCurrent;
      firstMiss = true;
   }
   fLastMiss = fEntryCurrent;
   // The first time this is executed, we try to pull in as much data as we can.
   TBranch *b = CalculateMissEntries(pos, len, firstMiss);
   if (!b) {
      if (!firstMiss) {
         // TODO: this recalculates for *all* branches, throwing away the above work.
         b = CalculateMissEntries(pos, len, true);
      }
      if (!b) {
         //printf("ProcessMiss: pos %ld does not appear to correspond to a buffer in this file.\n", pos);
         // We have gone through all the branches in this file and the requested basket
         // doesn't appear to be in any of them.  Likely a logic error / bug.
         fEntries->clear();
         return false;
      }
   }
   // TODO: this should be a set.
   fMissBranches->push_back(b);

   // OK, sort the entries
   std::sort(fEntries->begin(), fEntries->end(), [](const std::pair<ULong64_t, UInt_t> &_a, const std::pair<ULong64_t, UInt_t> &_b) {
      return _a.first < _b.first;
   });
   // Calculate the buffer offsets.  The data corresponding to read at fEntries[i]
   // is at offset fEntryOffsets[i] in the fMissCache.
   fEntryOffsets->reserve(fEntries->size());

   // Now, fetch the buffer.
   std::vector<Long64_t> positions; positions.reserve(fEntries->size());
   std::vector<Int_t> lengths; lengths.reserve(fEntries->size());
   int idx = 0;
   ULong64_t cumulative = 0;
   for (const auto &iopos : *fEntries) {
         positions.push_back(iopos.first);
         lengths.push_back(iopos.second);
         (*fEntryOffsets)[idx++] = cumulative;
         cumulative += iopos.second;
   }
   fMissCache->reserve(cumulative);
   //printf("Reading %lu bytes into miss cache for %lu entries.\n", cumulative, fEntries->size());
   fNMissReadPref += fEntries->size();
   fFile->ReadBuffers(&((*fMissCache)[0]), &positions[0], &lengths[0], fEntries->size());
   fFirstMiss = fLastMiss = fEntryCurrent;

   return true;
}

bool TTreeCache::CheckMissCache(char *buf, Long64_t pos, int len) {

   if (!fOptimizeMisses) {
      return false;
   }
   if ((pos < 0) || (len < 0)) {
      return false;
   }

   //printf("Checking the miss cache for offset=%ld, length=%d\n", pos, len);

   // First, binary search to see if the desired basket is already cached.
   auto iter = std::lower_bound(fEntries->begin(), fEntries->end(), std::make_pair(static_cast<ULong64_t>(pos), static_cast<UInt_t>(len)), [](const std::pair<ULong64_t, UInt_t> &_a, const std::pair<ULong64_t, UInt_t> &_b) {
      return _a.first < _b.first;
   });
   if (iter != fEntries->end()) {
      // TODO: What if there is a length mismatch?
      auto buffer_entry_in_cache = iter - fEntries->begin();
      auto offset = (*fEntryOffsets)[buffer_entry_in_cache];
      memcpy(buf, &((*fMissCache)[offset]), len);
      //printf("Returning data from pos=%ld in miss cache.\n", offset);
      fNMissReadOk++;
      return true;
   }

   //printf("Data not in miss cache.\n");

   // Update the cache, looking for this (pos, len).
   if (!ProcessMiss(pos, len)) {
      //printf("Unable to pull data into miss cache.\n");
      fNMissReadMiss++;
      return false;
   }

   // OK, we updated the cache with as much information as possible.  Seach again for
   // the entry we want.  TODO: ProcessMiss could return an offset.
   iter = std::lower_bound(fEntries->begin(), fEntries->end(), std::make_pair(pos, len), [](const std::pair<ULong64_t, UInt_t> &_a, const std::pair<ULong64_t, UInt_t> &_b) {
      return _a.first < _b.first;
   });
   if (iter != fEntries->end()) {
      auto buffer_entry_in_cache = iter - fEntries->begin();
      auto offset = (*fEntryOffsets)[buffer_entry_in_cache];
      //printf("Expecting data at offset %ld in miss cache.\n", offset);
      memcpy(buf, &((*fMissCache)[offset]), len);
      fNMissReadOk++;
      return true;
   }

   // This must be a logic bug.  ProcessMiss should return false if (pos, len) 
   // wasn't put into fEntries.
   fNMissReadMiss++;
   return false;
}


////////////////////////////////////////////////////////////////////////////////
/// Fill the cache buffer with the branches in the cache.

Bool_t TTreeCache::FillBuffer() {

   if (fNbranches <= 0) return kFALSE;
   TTree *tree = ((TBranch*)fBranches->UncheckedAt(0))->GetTree();
   Long64_t entry = tree->GetReadEntry();
   Long64_t fEntryCurrentMax = 0;

   if (fEnablePrefetching) { // Prefetching mode
      if (fIsLearning) { // Learning mode
         if (fEntryNext >= 0 && entry >= fEntryNext) {
            // entry is outside the learn range, need to stop the learning
            // phase. Doing so may trigger a recursive call to FillBuffer in
            // the process of filling both prefetching buffers
            StopLearningPhase();
            fIsManual = kFALSE;
         }
      }
      if (fIsLearning) { //  Learning mode
         entry = 0;
      }
      if (fFirstTime) {
         //try to detect if it is normal or reverse read
         fFirstEntry = entry;
      }
      else {
         if (fFirstEntry == entry) return kFALSE;
         // Set the read direction
         if (!fReadDirectionSet) {
            if (entry < fFirstEntry) {
               fReverseRead = kTRUE;
               fReadDirectionSet = kTRUE;
            }
            else if (entry > fFirstEntry) {
               fReverseRead =kFALSE;
               fReadDirectionSet = kTRUE;
            }
         }

         if (fReverseRead) {
            // Reverse reading with prefetching
            if (fEntryCurrent >0 && entry < fEntryNext) {
               // We can prefetch the next buffer
               if (entry >= fEntryCurrent) {
                  entry = fEntryCurrent - tree->GetAutoFlush() * fFillTimes;
               }
               if (entry < 0) entry = 0;
            }
            else if (fEntryCurrent >= 0) {
               // We are still reading from the oldest buffer, no need to prefetch a new one
               return kFALSE;
            }
            if (entry < 0) return kFALSE;
            fFirstBuffer = !fFirstBuffer;
         }
         else {
            // Normal reading with prefetching
            if (fEnablePrefetching) {
               if (entry < 0 && fEntryNext > 0) {
                  entry = fEntryCurrent;
               } else if (entry >= fEntryCurrent) {
                  if (entry < fEntryNext) {
                     entry = fEntryNext;
                  }
               }
               else {
                  // We are still reading from the oldest buffer,
                  // no need to prefetch a new one
                  return kFALSE;
               }
               fFirstBuffer = !fFirstBuffer;
            }
         }
      }
   }

   // If the entry is in the range we previously prefetched, there is
   // no point in retrying.   Note that this will also return false
   // during the training phase (fEntryNext is then set intentional to
   // the end of the training phase).
   if (fEntryCurrent <= entry && entry < fEntryNext) return kFALSE;

   // Triggered by the user, not the learning phase
   if (entry == -1)  entry = 0;

   fEntryCurrentMax = fEntryCurrent;
   TTree::TClusterIterator clusterIter = tree->GetClusterIterator(entry);
   fEntryCurrent = clusterIter();
   fEntryNext = clusterIter.GetNextEntry();

   if (fEntryCurrent < fEntryMin) fEntryCurrent = fEntryMin;
   if (fEntryMax <= 0) fEntryMax = tree->GetEntries();
   if (fEntryNext > fEntryMax) fEntryNext = fEntryMax;

   if ( fEnablePrefetching ) {
      if ( entry == fEntryMax ) {
         // We are at the end, no need to do anything else
         return kFALSE;
      }
   }

   // Check if owner has a TEventList set. If yes we optimize for this
   // Special case reading only the baskets containing entries in the
   // list.
   TEventList *elist = fTree->GetEventList();
   Long64_t chainOffset = 0;
   if (elist) {
      if (fTree->IsA() ==TChain::Class()) {
         TChain *chain = (TChain*)fTree;
         Int_t t = chain->GetTreeNumber();
         chainOffset = chain->GetTreeOffset()[t];
      }
   }

   //clear cache buffer
   Int_t fNtotCurrentBuf = 0;
   if (fEnablePrefetching){ //prefetching mode
      if (fFirstBuffer) {
         TFileCacheRead::Prefetch(0,0);
         fNtotCurrentBuf = fNtot;
      }
      else {
         TFileCacheRead::SecondPrefetch(0,0);
         fNtotCurrentBuf = fBNtot;
      }
   }
   else {
      TFileCacheRead::Prefetch(0,0);
      fNtotCurrentBuf = fNtot;
   }

   //store baskets
   Int_t clusterIterations = 0;
   Long64_t minEntry = fEntryCurrent;
   Int_t prevNtot;
   Int_t minBasket = 0;  // We will use this to avoid re-checking the first baskets in the 2nd (or more) run in the while loop.
   Long64_t maxReadEntry = minEntry; // If we are stopped before the end of the 2nd pass, this marker will where we need to start next time.
   do {
      prevNtot = fNtotCurrentBuf;
      Int_t nextMinBasket = INT_MAX;
      UInt_t pass = 0;
      while (pass < 2) {
         // The first pass we add one basket per branches.
         // then in the second pass we add the other baskets of the cluster.
         // This is to support the case where the cache is too small to hold a full cluster.
         ++pass;
         for (Int_t i=0;i<fNbranches;i++) {
            TBranch *b = (TBranch*)fBranches->UncheckedAt(i);
            if (b->GetDirectory()==0) continue;
            if (b->GetDirectory()->GetFile() != fFile) continue;
            Int_t nb = b->GetMaxBaskets();
            Int_t *lbaskets   = b->GetBasketBytes();
            Long64_t *entries = b->GetBasketEntry();
            if (!lbaskets || !entries) continue;
            //we have found the branch. We now register all its baskets
            //from the requested offset to the basket below fEntrymax
            Int_t blistsize = b->GetListOfBaskets()->GetSize();
            Int_t j=minBasket;  // We need this out of the loop so we can find out how far we went.
            Bool_t firstBasketSeen = kFALSE;
            for (;j<nb;j++) {
               // This basket has already been read, skip it
               if (j<blistsize && b->GetListOfBaskets()->UncheckedAt(j)) continue;

               Long64_t pos = b->GetBasketSeek(j);
               Int_t len = lbaskets[j];
               if (pos <= 0 || len <= 0) continue;
               if (len > fBufferSizeMin) {
                  // Do not cache a basket if it is bigger than the cache size!
                  continue;
               }
               //important: do not try to read fEntryNext, otherwise you jump to the next autoflush
               if (entries[j] >= fEntryNext) break; // break out of the for each branch loop.
               if (entries[j] < minEntry && (j<nb-1 && entries[j+1] <= minEntry)) continue;
               if (elist) {
                  Long64_t emax = fEntryMax;
                  if (j<nb-1) emax = entries[j+1]-1;
                  if (!elist->ContainsRange(entries[j]+chainOffset,emax+chainOffset)) continue;
               }
               if (pass==2 && !firstBasketSeen) {
                  // Okay, this has already been requested in the first pass.
                  firstBasketSeen = kTRUE;
                  continue;
               }
               fNReadPref++;

               if ( (fNtotCurrentBuf+len) > fBufferSizeMin ) {
                  // Humm ... we are going to go over the requested size.
                  if (clusterIterations > 0) {
                     // We already have a full cluster and now we would go over the requested
                     // size, let's stop caching (and make sure we start next time from the
                     // end of the previous cluster).
                     if (gDebug > 5) {
                        Info("FillBuffer","Breaking early because %d is greater than %d at cluster iteration %d will restart at %lld",(fNtotCurrentBuf+len), fBufferSizeMin, clusterIterations,minEntry);
                     }
                     fEntryNext = minEntry;
                     break;
                  } else {
                     if (pass == 1) {
                        if ( (fNtotCurrentBuf+len) > 4*fBufferSizeMin ) {
                           // Okay, so we have not even made one pass and we already have
                           // accumulated request for more than twice the memory size ...
                           // So stop for now, and will restart at the same point, hoping
                           // that the basket will still be in memory and not asked again ..
                           fEntryNext = maxReadEntry;
                           if (gDebug > 5) {
                              Info("FillBuffer","Breaking early because %d is greater than 2*%d at cluster iteration %d pass %d will restart at %lld",(fNtotCurrentBuf+len), fBufferSizeMin, clusterIterations,pass,fEntryNext);
                           }
                           break;
                        }
                     } else {
                        // We have made one pass through the branches and thus already
                        // requested one basket per branch, let's stop prefetching
                        // now.
                        if ( (fNtotCurrentBuf+len) > 2*fBufferSizeMin ) {
                           fEntryNext = maxReadEntry;
                           if (gDebug > 5) {
                              Info("FillBuffer","Breaking early because %d is greater than 2*%d at cluster iteration %d pass %d will restart at %lld",(fNtotCurrentBuf+len), fBufferSizeMin, clusterIterations,pass,fEntryNext);
                           }
                           break;
                        }
                     }
                  }
               }
               if (fEnablePrefetching){
                  if (fFirstBuffer) {
                     TFileCacheRead::Prefetch(pos,len);
                     fNtotCurrentBuf = fNtot;
                  }
                  else {
                     TFileCacheRead::SecondPrefetch(pos,len);
                     fNtotCurrentBuf = fBNtot;
                  }
               }
               else {
                  TFileCacheRead::Prefetch(pos,len);
                  fNtotCurrentBuf = fNtot;
               }
               if ( ( j < (nb-1) ) && entries[j+1] > maxReadEntry ) {
                  maxReadEntry = entries[j+1];
               }
               if (fNtotCurrentBuf > 4*fBufferSizeMin) {
                  // Humm something wrong happened.
                  Warning("FillBuffer","There is more data in this cluster (starting at entry %lld to %lld, current=%lld) than usual ... with %d %.3f%% of the branches we already have %d bytes (instead of %d)",
                          fEntryCurrent,fEntryNext, entries[j], i, (100.0*i) / ((float)fNbranches), fNtotCurrentBuf,fBufferSizeMin);

               }
               if (pass==1) {
                  // In the first pass, we record one basket per branch and move on to the next branch.
                  break;
               }
            }

            if (j < nextMinBasket) nextMinBasket = j;
            if (gDebug > 0) printf("Entry: %lld, registering baskets branch %s, fEntryNext=%lld, fNseek=%d, fNtotCurrentBuf=%d\n",minEntry,((TBranch*)fBranches->UncheckedAt(i))->GetName(),fEntryNext,fNseek,fNtotCurrentBuf);
         }
      } // loop for the 2 passes.
      clusterIterations++;

      minEntry = clusterIter.Next();
      if (fIsLearning) {
         fFillTimes++;
      }

      // Continue as long as we still make progress (prevNtot < fNtotCurrentBuf), that the next entry range to be looked at,
      // which start at 'minEntry', is not past the end of the requested range (minEntry < fEntryMax)
      // and we guess that we not going to go over the requested amount of memory by asking for another set
      // of entries (fBufferSizeMin > ((Long64_t)fNtotCurrentBuf*(clusterIterations+1))/clusterIterations).
      // fNtotCurrentBuf / clusterIterations is the average size we are accumulated so far at each loop.
      // and thus (fNtotCurrentBuf / clusterIterations) * (clusterIterations+1) is a good guess at what the next total size
      // would be if we run the loop one more time.   fNtotCurrentBuf and clusterIterations are Int_t but can sometimes
      // be 'large' (i.e. 30Mb * 300 intervals) and can overflow the numerical limit of Int_t (i.e. become
      // artificially negative).   To avoid this issue we promote fNtotCurrentBuf to a long long (64 bits rather than 32 bits)
      if (!((fBufferSizeMin > ((Long64_t)fNtotCurrentBuf*(clusterIterations+1))/clusterIterations) && (prevNtot < fNtotCurrentBuf) && (minEntry < fEntryMax)))
         break;

      //for the reverse reading case
      if (!fIsLearning && fReverseRead){
         if (clusterIterations >= fFillTimes)
            break;
         if (minEntry >= fEntryCurrentMax && fEntryCurrentMax >0)
            break;
      }
      minBasket = nextMinBasket;
      fEntryNext = clusterIter.GetNextEntry();
      if (fEntryNext > fEntryMax) fEntryNext = fEntryMax;
   } while (kTRUE);

   if (fEnablePrefetching) {
      if (fIsLearning) {
         fFirstBuffer = !fFirstBuffer;
      }
      if (!fIsLearning && fFirstTime){
         // First time we add autoFlush entries , after fFillTimes * autoFlush
         // only in reverse prefetching mode
         fFirstTime = kFALSE;
      }
   }
   fIsLearning = kFALSE;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the desired prefill type from the environment or resource variable
/// - 0 - No prefill
/// - 1 - All branches

TTreeCache::EPrefillType TTreeCache::GetConfiguredPrefillType() const
{
   const char *stcp;
   Int_t s = 0;

   if (!(stcp = gSystem->Getenv("ROOT_TTREECACHE_PREFILL")) || !*stcp) {
      s = gEnv->GetValue("TTreeCache.Prefill", 1);
   } else {
      s = TString(stcp).Atoi();
   }

   return static_cast<TTreeCache::EPrefillType>(s);
}

////////////////////////////////////////////////////////////////////////////////
/// Give the total efficiency of the primary cache... defined as the ratio
/// of blocks found in the cache vs. the number of blocks prefetched
/// ( it could be more than 1 if we read the same block from the cache more
///   than once )
///
/// Note: This should eb used at the end of the processing or we will
/// get incomplete stats

Double_t TTreeCache::GetEfficiency() const
{
   if ( !fNReadPref )
      return 0;

   return ((Double_t)fNReadOk / (Double_t)fNReadPref);
}

/**
 * The total efficiency of the 'miss cache' - defined as the ratio
 * of blocks found in the cache versus the number of blocks prefetched
 */
double TTreeCache::GetMissEfficiency() const
{
   if ( !fNMissReadPref ) {return 0;}
   return static_cast<double>(fNMissReadOk) / static_cast<double>(fNMissReadPref);
}

////////////////////////////////////////////////////////////////////////////////
/// This will indicate a sort of relative efficiency... a ratio of the
/// reads found in the cache to the number of reads so far

Double_t TTreeCache::GetEfficiencyRel() const
{
   if ( !fNReadOk && !fNReadMiss )
      return 0;

   return ((Double_t)fNReadOk / (Double_t)(fNReadOk + fNReadMiss));
}

/**
 * Relative efficiency of the 'miss cache' - ratio of the reads found in cache
 * to the number of reads so far.
 */
double TTreeCache::GetMissEfficiencyRel() const
{
   if (!fNMissReadOk && !fNMissReadMiss) {return 0;}

   return static_cast<double>(fNMissReadOk) / static_cast<double>(fNMissReadOk + fNMissReadMiss);
}

////////////////////////////////////////////////////////////////////////////////
/// Static function returning the number of entries used to train the cache
/// see SetLearnEntries

Int_t TTreeCache::GetLearnEntries()
{
   return fgLearnEntries;
}

////////////////////////////////////////////////////////////////////////////////
/// Print cache statistics. Like:
///
/// ~~~ {.cpp}
///    ******TreeCache statistics for file: cms2.root ******
///    Number of branches in the cache ...: 1093
///    Cache Efficiency ..................: 0.997372
///    Cache Efficiency Rel...............: 1.000000
///    Learn entries......................: 100
///    Reading............................: 72761843 bytes in 7 transactions
///    Readahead..........................: 256000 bytes with overhead = 0 bytes
///    Average transaction................: 10394.549000 Kbytes
///    Number of blocks in current cache..: 210, total size: 6280352
/// ~~~
///
/// - if option = "a" the list of blocks in the cache is printed
///   see also class TTreePerfStats.
/// - if option contains 'cachedbranches', the list of branches being
///   cached is printed.

void TTreeCache::Print(Option_t *option) const
{
   TString opt = option;
   opt.ToLower();
   printf("******TreeCache statistics for tree: %s in file: %s ******\n",fTree ? fTree->GetName() : "no tree set",fFile ? fFile->GetName() : "no file set");
   if (fNbranches <= 0) return;
   printf("Number of branches in the cache ...: %d\n",fNbranches);
   printf("Cache Efficiency ..................: %f\n",GetEfficiency());
   printf("Cache Efficiency Rel...............: %f\n",GetEfficiencyRel());
   printf("Secondary Efficiency ..............: %f\n", GetMissEfficiency());
   printf("Secondary Efficiency Rel ..........: %f\n", GetMissEfficiencyRel());
   printf("Learn entries......................: %d\n",TTreeCache::GetLearnEntries());
   if ( opt.Contains("cachedbranches") ) {
      opt.ReplaceAll("cachedbranches","");
      printf("Cached branches....................:\n");
      const TObjArray *cachedBranches = this->GetCachedBranches();
      Int_t nbranches = cachedBranches->GetEntriesFast();
      for (Int_t i = 0; i < nbranches; ++i) {
         TBranch* branch = (TBranch*) cachedBranches->UncheckedAt(i);
         printf("Branch name........................: %s\n",branch->GetName());
      }
   }
   TFileCacheRead::Print(opt);
}

////////////////////////////////////////////////////////////////////////////////
/// Old method ReadBuffer before the addition of the prefetch mechanism.

Int_t TTreeCache::ReadBufferNormal(char *buf, Long64_t pos, Int_t len){
   //Is request already in the cache?
   if (TFileCacheRead::ReadBuffer(buf,pos,len) == 1){
      fNReadOk++;
      return 1;
   }

   //not found in cache. Do we need to fill the cache?
   Bool_t bufferFilled = FillBuffer();
   if (bufferFilled) {
      Int_t res = TFileCacheRead::ReadBuffer(buf,pos,len);

      if (res == 1)
         fNReadOk++;
      else if (res == 0)
         fNReadMiss++;

      return res;
   }
   if (CheckMissCache(buf, pos, len)) {return 1;}

   fNReadMiss++;

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Used to read a chunk from a block previously fetched. It will call FillBuffer
/// even if the cache lookup succeeds, because it will try to prefetch the next block
/// as soon as we start reading from the current block.

Int_t TTreeCache::ReadBufferPrefetch(char *buf, Long64_t pos, Int_t len){
   if (TFileCacheRead::ReadBuffer(buf, pos, len) == 1){
      //call FillBuffer to prefetch next block if necessary
      //(if we are currently reading from the last block available)
      FillBuffer();
      fNReadOk++;
      return 1;
   }

   //keep on prefetching until request is satisfied
   // try to prefetch a couple of times and if request is still not satisfied then
   // fall back to normal reading without prefetching for the current request
   Int_t counter = 0;
   while (1) {
      if(TFileCacheRead::ReadBuffer(buf, pos, len)) {
         break;
      }
      FillBuffer();
      fNReadMiss++;
      counter++;
      if (counter>1) {
        return 0;
      }
   }

   fNReadOk++;
   return 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Read buffer at position pos if the request is in the list of
/// prefetched blocks read from fBuffer.
/// Otherwise try to fill the cache from the list of selected branches,
/// and recheck if pos is now in the list.
/// Returns:
///  - -1 in case of read failure,
///  - 0 in case not in cache,
///  - 1 in case read from cache.
/// This function overloads TFileCacheRead::ReadBuffer.

Int_t TTreeCache::ReadBuffer(char *buf, Long64_t pos, Int_t len)
{
   if (!fEnabled) return 0;

   if (fEnablePrefetching)
      return TTreeCache::ReadBufferPrefetch(buf, pos, len);
   else
      return TTreeCache::ReadBufferNormal(buf, pos, len);
}

////////////////////////////////////////////////////////////////////////////////
/// This will simply clear the cache

void TTreeCache::ResetCache()
{
   TFileCacheRead::Prefetch(0,0);

   if (fEnablePrefetching) {
      fFirstTime = kTRUE;
      TFileCacheRead::SecondPrefetch(0, 0);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Change the underlying buffer size of the cache.
/// If the change of size means some cache content is lost, or if the buffer
/// is now larger, setup for a cache refill the next time there is a read
/// Returns:
///  - 0 if the buffer content is still available
///  - 1 if some or all of the buffer content has been made unavailable
///  - -1 on error

Int_t TTreeCache::SetBufferSize(Int_t buffersize)
{
   Int_t prevsize = GetBufferSize();
   Int_t res = TFileCacheRead::SetBufferSize(buffersize);
   if (res < 0) {
      return res;
   }

   if (res == 0 && buffersize <= prevsize) {
      return res;
   }

   // if content was removed from the buffer, or the buffer was enlarged then
   // empty the prefetch lists and prime to fill the cache again

   TFileCacheRead::Prefetch(0,0);
   if (fEnablePrefetching) {
      TFileCacheRead::SecondPrefetch(0, 0);
   }

   fEntryCurrent = -1;
   if (!fIsLearning) {
      fEntryNext = -1;
   }

   return 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the minimum and maximum entry number to be processed
/// this information helps to optimize the number of baskets to read
/// when prefetching the branch buffers.

void TTreeCache::SetEntryRange(Long64_t emin, Long64_t emax)
{
   // This is called by TTreePlayer::Process in an automatic way...
   // don't restart it if the user has specified the branches.
   Bool_t needLearningStart = (fEntryMin != emin) && fIsLearning && !fIsManual;

   fEntryMin  = emin;
   fEntryMax  = emax;
   fEntryNext  = fEntryMin + fgLearnEntries * (fIsLearning && !fIsManual);
   if (gDebug > 0)
      Info("SetEntryRange", "fEntryMin=%lld, fEntryMax=%lld, fEntryNext=%lld",
                             fEntryMin, fEntryMax, fEntryNext);

   if (needLearningStart) {
      // Restart learning
      StartLearningPhase();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Overload to make sure that the object specific

void TTreeCache::SetFile(TFile *file, TFile::ECacheAction action)
{
   // The infinite recursion is 'broken' by the fact that
   // TFile::SetCacheRead remove the entry from fCacheReadMap _before_
   // calling SetFile (and also by setting fFile to zero before the calling).
   if (fFile) {
      TFile *prevFile = fFile;
      fFile = 0;
      prevFile->SetCacheRead(0, fTree, action);
   }
   TFileCacheRead::SetFile(file, action);
}

////////////////////////////////////////////////////////////////////////////////
/// Static function to set the number of entries to be used in learning mode
/// The default value for n is 10. n must be >= 1

void TTreeCache::SetLearnEntries(Int_t n)
{
   if (n < 1) n = 1;
   fgLearnEntries = n;
}

////////////////////////////////////////////////////////////////////////////////
/// Set whether the learning period is started with a prefilling of the
/// cache and which type of prefilling is used.
/// The two value currently supported are:
///  - TTreeCache::kNoPrefill    disable the prefilling
///  - TTreeCache::kAllBranches  fill the cache with baskets from all branches.
/// The default prefilling behavior can be controlled by setting
/// TTreeCache.Prefill or the environment variable ROOT_TTREECACHE_PREFILL.

void TTreeCache::SetLearnPrefill(TTreeCache::EPrefillType type /* = kNoPrefill */)
{
   fPrefillType = type;
}

////////////////////////////////////////////////////////////////////////////////
/// The name should be enough to explain the method.
/// The only additional comments is that the cache is cleaned before
/// the new learning phase.

void TTreeCache::StartLearningPhase()
{
   fIsLearning = kTRUE;
   fIsManual = kFALSE;
   fNbranches  = 0;
   if (fBrNames) fBrNames->Delete();
   fIsTransferred = kFALSE;
   fEntryCurrent = -1;
}

////////////////////////////////////////////////////////////////////////////////
/// This is the counterpart of StartLearningPhase() and can be used to stop
/// the learning phase. It's useful when the user knows exactly what branches
/// they are going to use.
/// For the moment it's just a call to FillBuffer() since that method
/// will create the buffer lists from the specified branches.

void TTreeCache::StopLearningPhase()
{
   if (fIsLearning) {
      // This will force FillBuffer to read the buffers.
      fEntryNext = -1;
      fIsLearning = kFALSE;
   }
   fIsManual = kTRUE;

   //fill the buffers only once during learning
   if (fEnablePrefetching && !fOneTime) {
      fIsLearning = kTRUE;
      FillBuffer();
      fOneTime = kTRUE;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Update pointer to current Tree and recompute pointers to the branches in the cache.

void TTreeCache::UpdateBranches(TTree *tree)
{

   fTree = tree;

   fEntryMin  = 0;
   fEntryMax  = fTree->GetEntries();

   fEntryCurrent = -1;

   if (fBrNames->GetEntries() == 0 && fIsLearning) {
      // We still need to learn.
      fEntryNext = fEntryMin + fgLearnEntries;
   } else {
      // We learnt from a previous file.
      fIsLearning = kFALSE;
      fEntryNext = -1;
   }
   fNbranches = 0;

   TIter next(fBrNames);
   TObjString *os;
   while ((os = (TObjString*)next())) {
      TBranch *b = fTree->GetBranch(os->GetName());
      if (!b) {
         continue;
      }
      fBranches->AddAt(b, fNbranches);
      fNbranches++;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Perform an initial prefetch, attempting to read as much of the learning
/// phase baskets for all branches at once

void TTreeCache::LearnPrefill()
{
   // This is meant for the learning phase
   if (!fIsLearning) return;

   // This should be called before reading entries, otherwise we'll
   // always exit here, since TBranch adds itself before reading
   if (fNbranches > 0) return;

   // Is the LearnPrefill enabled (using an Int_t here to allow for future
   // extension to alternative Prefilling).
   if (fPrefillType == kNoPrefill) return;

   // Force only the learn entries to be cached by temporarily setting min/max
   // to the learning phase entry range
   // But save all the old values, so we can restore everything to how it was
   Long64_t eminOld = fEntryMin;
   Long64_t emaxOld = fEntryMax;
   Long64_t ecurrentOld = fEntryCurrent;
   Long64_t enextOld = fEntryNext;

   fEntryMin = fEntryCurrent;
   fEntryMax = fEntryNext;

   // Add all branches to be cached. This also sets fIsManual, stops learning,
   // and makes fEntryNext = -1 (which forces a cache fill, which is good)
   AddBranch("*");
   fIsManual = kFALSE; // AddBranch sets fIsManual, so we reset it

   // Now, fill the buffer with the learning phase entry range
   FillBuffer();

   // Leave everything the way we found it
   fIsLearning = kTRUE;
   DropBranch("*"); // This doesn't work unless we're already learning

   // Restore entry values
   fEntryMin = eminOld;
   fEntryMax = emaxOld;
   fEntryCurrent = ecurrentOld;
   fEntryNext = enextOld;
}
