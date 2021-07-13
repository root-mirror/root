// @(#)root/tree:$Id$
// Author: Rene Brun   04/06/2006

/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TTreeCache
\ingroup tree
\brief A cache to speed-up the reading of ROOT datasets

# A cache to speed-up the reading of ROOT datasets

## Table of Contents
- [Motivation](\ref motivation)
- [General Description](\ref description)
- [Changes in behaviour](\ref changesbehaviour)
- [Self-optimization](\ref cachemisses)
- [Examples of usage](\ref examples)
- [Check performance and stats](\ref checkPerf)

\anchor motivation
## Motivation: why having a cache is needed?

When writing a TTree, the branch buffers are kept in memory.
A typical branch buffersize (before compression) is typically 32 KBytes.
After compression, the zipped buffer may be just a few Kbytes.
The branch buffers cannot be much larger in case of TTrees with several
hundred or thousand branches.

When writing, this does not generate a performance problem because branch
buffers are always written sequentially and, thanks to OS optimisations,
content is flushed to the output file when a few MBytes of data are available.
On the other hand, when reading, one may hit performance problems because of
latencies e.g imposed by network.
For example in a WAN with 10ms latency, reading 1000 buffers of 10 KBytes each
with no cache will imply 10s penalty where a local read of the 10 MBytes would
take about 1 second.

The TreeCache tries to prefetch all the buffers for the selected branches
in order to transfer a few multi-Megabytes large buffers instead of many
multi-kilobytes small buffers. In addition, TTreeCache can sort the blocks to
be read in increasing order such that the file is read sequentially.

Systems like xrootd, dCache or httpd take advantage of the TTreeCache in
reading ahead as much data as they can and return to the application
the maximum data specified in the cache and have the next chunk of data ready
when the next request comes.

### Are there cases for which the usage of TTreeCache is detrimental for performance?
Yes, some corner cases. For example, when reading only a small fraction of all
entries such that not all branch buffers are read.

\anchor description
## General Description
This class acts as a file cache, registering automatically the baskets from
the branches being processed via direct manipulation of TTrees or with tools
such as TTree::Draw, TTree::Process, TSelector, TTreeReader and RDataFrame
when in the learning phase. The learning phase is by default 100 entries.
It can be changed via TTreeCache::SetLearnEntries.

The usage of a TTreeCache can considerably improve the runtime performance at
the price of a modest investment in memory, in particular when the TTree is
accessed remotely, e.g. via a high latency network.

For each TTree being processed a TTreeCache object is created.
This object is automatically deleted when the Tree is deleted or
when the file is deleted.
The user can change the size of the cache with the TTree::SetCacheSize method
(by default the size is 30 Megabytes). This feature can be controlled with the
environment variable `ROOT_TTREECACHE_SIZE` or the TTreeCache.Size option.
The entry range for which the cache is active can also be set with the
SetEntryRange method.

\anchor changesbehaviour
## Changes of behavior when using TChain and TEventList

The usage of TChain or TEventList have influence on the behaviour of the cache:

- Special case of a TChain
  Once the training is done on the first Tree, the list of branches
  in the cache is kept for the following files.

- Special case of a TEventlist
  if the Tree or TChain has a TEventlist, only the buffers
  referenced by the list are put in the cache.

The learning phase is started or restarted when:
   - TTree automatically creates a cache.
   - TTree::SetCacheSize is called with a non-zero size and a cache
     did not previously exist
   - TTreeCache::StartLearningPhase is called.
   - TTreeCache::SetEntryRange is called
        * and the learning is not yet finished
        * and has not been set to manual
        * and the new minimun entry is different.

The learning period is stopped (and prefetching is started) when:
   - TTreeCache::StopLearningPhase is called.
   - An entry outside the 'learning' range is requested
     The 'learning range is from fEntryMin (default to 0) to
     fEntryMin + fgLearnEntries.
   - A 'cached' TChain switches over to a new file.


\anchor cachemisses
## Self-optimization in presence of cache misses

The TTreeCache can optimize its behavior on a cache miss. When
miss optimization is enabled (see the SetOptimizeMisses method),
it tracks all branches utilized after the learning phase which caused a cache
miss.
When one cache miss occurs, all the utilized branches are be prefetched
for that event. This optimization utilizes the observation that infrequently
accessed branches are often accessed together.
An example scenario where such behavior is desirable, is an analysis where
a set of collections are read only for a few events in which a certain
condition is respected, e.g. a trigger fired.

### Additional memory and CPU usage when optimizing for cache misses
When this mode is enabled, the memory dedicated to the cache can increase
by at most a factor two in the case of cache miss.
Additionally, on the first miss of an event, we must iterate through all the
"active branches" for the miss cache and find the correct basket.
This can be potentially a CPU-expensive operation compared to, e.g., the
latency of a SSD.  This is why the miss cache is currently disabled by default.

\anchor examples
## Example usages of TTreeCache

A few use cases are discussed below. A cache may be created with automatic
sizing when a TTree is used:

In some applications, e.g. central processing workflows of experiments, the list
of branches to read is known a priori. For these cases, the TTreeCache can be
instructed about the branches which will be read via explicit calls to the TTree
or TTreeCache interfaces.
In less streamlined applications such as analysis, predicting the branches which
will be read can be difficult. In such cases, ROOT I/O flags used branches
automatically when a branch buffer is read during the learning phase.

In the examples below, portions of analysis code are shown.
The few statements involving the TreeCache are marked with `//<<<`

### ROOT::RDataFrame and TTreeReader Examples

If you use RDataFrame or TTreeReader, the system will automatically cache the
best set of branches: no action is required by the user.

### TTree::Draw Example

The TreeCache is automatically used by TTree::Draw. The method knows
which branches are used in the query and it puts automatically these branches
in the cache. The entry range is also inferred automatically.

### TTree::Process and TSelectors Examples

The user must enable the cache and tell the system which branches to cache
and also specify the entry range. It is important to specify the entry range
in case only a subset of the events is processed to avoid wasteful caching.

#### Reading all branches

~~~ {.cpp}
    TTree *T;
    f->GetObject(T, "mytree");
    auto nentries = T->GetEntries();
    auto cachesize = 10000000U; // 10 MBytes
    T->SetCacheSize(cachesize); //<<<
    T->AddBranchToCache("*", true);    //<<< add all branches to the cache
    T->Process("myselector.C+");
    // In the TSelector::Process function we read all branches
    T->GetEntry(i);
    // ... Here the entry is processed
~~~

#### Reading a subset of all branches

In the Process function we read a subset of the branches.
Only the branches used in the first entry will be put in the cache
~~~ {.cpp}
    TTree *T;
    f->GetObject(T, "mytree");
    // We want to process only the 200 first entries
    auto nentries=200UL;
    auto efirst = 0;
    auto elast = efirst+nentries;
    auto cachesize = 10000000U; // 10 MBytes
    TTreeCache::SetLearnEntries(1);  //<<< we can take the decision after 1 entry
    T->SetCacheSize(cachesize);      //<<<
    T->SetCacheEntryRange(efirst,elast); //<<<
    T->Process("myselector.C+","",nentries,efirst);
    // In the TSelector::Process we read only 2 branches
    auto b1 = T->GetBranch("branch1");
    b1->GetEntry(i);
    if (somecondition) return;
    auto b2 = T->GetBranch("branch2");
    b2->GetEntry(i);
    ... Here the entry is processed
~~~
### Custom event loop

#### Always using the same two branches

In this example, exactly two branches are always used: those need to be
prefetched.
~~~ {.cpp}
    TTree *T;
    f->GetObject(T, "mytree");
    auto b1 = T->GetBranch("branch1");
    auto b2 = T->GetBranch("branch2");
    auto nentries = T->GetEntries();
    auto cachesize = 10000000U; //10 MBytes
    T->SetCacheSize(cachesize);     //<<<
    T->AddBranchToCache(b1, true);  //<<< add branch1 and branch2 to the cache
    T->AddBranchToCache(b2, true);  //<<<
    T->StopCacheLearningPhase();    //<<< we do not need the system to guess anything
    for (auto i : TSeqL(nentries)) {
       T->LoadTree(i); //<<< important call when calling TBranch::GetEntry after
       b1->GetEntry(i);
       if (some condition not met) continue;
       b2->GetEntry(i);
       if (some condition not met) continue;
       // Here we read the full event only in some rare cases.
       // There is no point in caching the other branches as it might be
       // more economical to read only the branch buffers really used.
       T->GetEntry(i);
       ... Here the entry is processed
    }
~~~
#### Always using at least the same two branches

In this example, two branches are always used: in addition, some analysis
functions are invoked and those may trigger the reading of other branches which
are a priori not known.
There is no point in prefetching branches that will be used very rarely: we can
rely on the system to cache the right branches.
~~~ {.cpp}
    TTree *T;
    f->GetObject(T, "mytree");
    auto nentries = T->GetEntries();
    auto cachesize = 10000000;   //10 MBytes
    T->SetCacheSize(cachesize);   //<<<
    T->SetCacheLearnEntries(5);   //<<< we can take the decision after 5 entries
    auto b1 = T->GetBranch("branch1");
    auto b2 = T->GetBranch("branch2");
    for (auto i : TSeqL(nentries)) {
       T->LoadTree(i);
       b1->GetEntry(i);
       if (some condition not met) continue;
       b2->GetEntry(i);
       // At this point we may call a user function where a few more branches
       // will be read conditionally. These branches will be put in the cache
       // if they have been used in the first 10 entries
       if (some condition not met) continue;
       // Here we read the full event only in some rare cases.
       // There is no point in caching the other branches as it might be
       // more economical to read only the branch buffers really used.
       T->GetEntry(i);
       .. process the rare but interesting cases.
       ... Here the entry is processed
    }
~~~

\anchor checkPerf
## How can the usage and performance of TTreeCache be verified?

Once the event loop terminated, the number of effective system reads for a
given file can be checked with a code like the following:
~~~ {.cpp}
    printf("Reading %lld bytes in %d transactions\n",myTFilePtr->GetBytesRead(),  f->GetReadCalls());
~~~

Another handy command is:
~~~ {.cpp}
myTreeOrChain.GetTree()->PrintCacheStats();
~~~

*/

#include "TSystem.h"
#include "TEnv.h"
#include "TTreeCache.h"
#include "TChain.h"
#include "TList.h"
#include "TBranch.h"
#include "TBranchElement.h"
#include "TEventList.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TRegexp.h"
#include "TLeaf.h"
#include "TFriendElement.h"
#include "TFile.h"
#include "TMath.h"
#include "TBranchCacheInfo.h"
#include "TVirtualPerfStats.h"
#include <limits.h>

Int_t TTreeCache::fgLearnEntries = 100;

ClassImp(TTreeCache);

////////////////////////////////////////////////////////////////////////////////
/// Default Constructor.

TTreeCache::TTreeCache() : TFileCacheRead(), fPrefillType(GetConfiguredPrefillType())
{
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TTreeCache::TTreeCache(TTree *tree, Int_t buffersize)
   : TFileCacheRead(tree->GetCurrentFile(), buffersize, tree), fEntryMax(tree->GetEntriesFast()), fEntryNext(0),
     fBrNames(new TList), fTree(tree), fPrefillType(GetConfiguredPrefillType())
{
   fEntryNext = fEntryMin + fgLearnEntries;
   Int_t nleaves = tree->GetListOfLeaves()->GetEntriesFast();
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
/// Add a branch discovered by actual usage to the list of branches to be stored
/// in the cache this function is called by TBranch::GetBasket
/// If we are not longer in the training phase this is an error.
/// Returns:
///  - 0 branch added or already included
///  - -1 on error

Int_t TTreeCache::LearnBranch(TBranch *b, Bool_t subbranches /*= kFALSE*/)
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
   if (!fLearnPrefilling && fNbranches == 0) LearnPrefill();

   return AddBranch(b, subbranches);
}

////////////////////////////////////////////////////////////////////////////////
/// Add a branch to the list of branches to be stored in the cache
/// this function is called by the user via TTree::AddBranchToCache.
/// The branch is added even if we are outside of the training phase.
/// Returns:
///  - 0 branch added or already included
///  - -1 on error

Int_t TTreeCache::AddBranch(TBranch *b, Bool_t subbranches /*= kFALSE*/)
{
   // Reject branch that are not from the cached tree.
   if (!b || fTree->GetTree() != b->GetTree()) return -1;

   //Is branch already in the cache?
   Bool_t isNew = kTRUE;
   for (int i=0;i<fNbranches;i++) {
      if (fBranches->UncheckedAt(i) == b) {isNew = kFALSE; break;}
   }
   if (isNew) {
      fTree = b->GetTree();
      fBranches->AddAtAndExpand(b, fNbranches);
      const char *bname = b->GetName();
      if (fTree->IsA() == TChain::Class()) {
         // If we have a TChain, we will need to use the branch name
         // and we better disambiguate them (see atlasFlushed.root for example)
         // in order to cache all the requested branches.
         // We do not do this all the time as GetMother is slow (it contains
         // a linear search from list of top level branch).
         TString build;
         const char *mothername = b->GetMother()->GetName();
         if (b != b->GetMother() && mothername[strlen(mothername)-1] != '.') {
            // Maybe we ought to prefix the name to avoid ambiguity.
            auto bem = dynamic_cast<TBranchElement*>(b->GetMother());
            if (bem->GetType() < 3) {
               // Not a collection.
               build = mothername;
               build.Append(".");
               if (strncmp(bname,build.Data(),build.Length()) != 0) {
                  build.Append(bname);
                  bname = build.Data();
               }
            }
         }
      }
      fBrNames->Add(new TObjString(bname));
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
            if (name != bname && AddBranch(name, subbranches)<0) {
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

////////////////////////////////////////////////////////////////////////////////
/// Start of methods for the miss cache.
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// Enable / disable the miss cache.
///
/// The first time this is called on a TTreeCache object, the corresponding
/// data structures will be allocated.  Subsequent enable / disables will
/// simply turn the functionality on/off.
void TTreeCache::SetOptimizeMisses(Bool_t opt)
{

   if (opt && !fMissCache) {
      ResetMissCache();
   }
   fOptimizeMisses = opt;
}

////////////////////////////////////////////////////////////////////////////////
/// Reset all the miss cache training.
///
/// The contents of the miss cache will be emptied as well as the list of
/// branches used.
void TTreeCache::ResetMissCache()
{

   fLastMiss = -1;
   fFirstMiss = -1;

   if (!fMissCache) {
      fMissCache.reset(new MissCache());
   }
   fMissCache->clear();
}

////////////////////////////////////////////////////////////////////////////////
/// For the event currently being fetched into the miss cache, find the IO
/// (offset / length tuple) to pull in the current basket for a given branch.
///
/// Returns:
/// - IOPos describing the IO operation necessary for the basket on this branch
/// - On failure, IOPos.length will be set to 0.
TTreeCache::IOPos TTreeCache::FindBranchBasketPos(TBranch &b, Long64_t entry)
{
   if (R__unlikely(b.GetDirectory() == 0)) {
      // printf("Branch at %p has no valid directory.\n", &b);
      return IOPos{0, 0};
   }
   if (R__unlikely(b.GetDirectory()->GetFile() != fFile)) {
      // printf("Branch at %p is in wrong file (branch file %p, my file %p).\n", &b, b.GetDirectory()->GetFile(),
      // fFile);
      return IOPos{0, 0};
   }

   // printf("Trying to find a basket for branch %p\n", &b);
   // Pull in metadata about branch; make sure it is valid
   Int_t *lbaskets = b.GetBasketBytes();
   Long64_t *entries = b.GetBasketEntry();
   if (R__unlikely(!lbaskets || !entries)) {
      // printf("No baskets or entries.\n");
      return IOPos{0, 0};
   }
   // Int_t blistsize = b.GetListOfBaskets()->GetSize();
   Int_t blistsize = b.GetWriteBasket();
   if (R__unlikely(blistsize <= 0)) {
      // printf("Basket list is size 0.\n");
      return IOPos{0, 0};
   }

   // Search for the basket that contains the event of interest.  Unlike the primary cache, we
   // are only interested in a single basket per branch - we don't try to fill the cache.
   Long64_t basketOffset = TMath::BinarySearch(blistsize, entries, entry);
   if (basketOffset < 0) { // No entry found.
      // printf("No entry offset found for entry %ld\n", fTree->GetReadEntry());
      return IOPos{0, 0};
   }

   // Check to see if there's already a copy of this basket in memory.  If so, don't fetch it
   if ((basketOffset < blistsize) && b.GetListOfBaskets()->UncheckedAt(basketOffset)) {

      // printf("Basket is already in memory.\n");
      return IOPos{0, 0};
   }

   Long64_t pos = b.GetBasketSeek(basketOffset);
   Int_t len = lbaskets[basketOffset];
   if (R__unlikely(pos <= 0 || len <= 0)) {
      /*printf("Basket returned was invalid (basketOffset=%ld, pos=%ld, len=%d).\n", basketOffset, pos, len);
      for (int idx=0; idx<blistsize; idx++) {
         printf("Basket entry %d, first event %d, pos %ld\n", idx, entries[idx], b.GetBasketSeek(idx));
      }*/
      return IOPos{0, 0};
   } // Sanity check
   // Do not cache a basket if it is bigger than the cache size!
   if (R__unlikely(len > fBufferSizeMin)) {
      // printf("Basket size is greater than the cache size.\n");
      return IOPos{0, 0};
   }

   return {pos, len};
}

////////////////////////////////////////////////////////////////////////////////
/// Given a particular IO description (offset / length) representing a 'miss' of
/// the TTreeCache's primary cache, calculate all the corresponding IO that
/// should be performed.
///
/// `all` indicates that this function should search the set of _all_ branches
/// in this TTree.  When set to false, we only search through branches that
/// have previously incurred a miss.
///
/// Returns:
/// - TBranch pointer corresponding to the basket that will be retrieved by
///   this IO operation.
/// - If no corresponding branch could be found (or an error occurs), this
///   returns nullptr.
TBranch *TTreeCache::CalculateMissEntries(Long64_t pos, Int_t len, Bool_t all)
{
   if (R__unlikely((pos < 0) || (len < 0))) {
      return nullptr;
   }

   int count = all ? (fTree->GetListOfLeaves())->GetEntriesFast() : fMissCache->fBranches.size();
   fMissCache->fEntries.reserve(count);
   fMissCache->fEntries.clear();
   Bool_t found_request = kFALSE;
   TBranch *resultBranch = nullptr;
   Long64_t entry = fTree->GetReadEntry();

   std::vector<std::pair<size_t, Int_t>> basketsInfo;
   auto perfStats = GetTree()->GetPerfStats();

   // printf("Will search %d branches for basket at %ld.\n", count, pos);
   for (int i = 0; i < count; i++) {
      TBranch *b =
         all ? static_cast<TBranch *>(static_cast<TLeaf *>((fTree->GetListOfLeaves())->UncheckedAt(i))->GetBranch())
             : fMissCache->fBranches[i];
      IOPos iopos = FindBranchBasketPos(*b, entry);
      if (iopos.fLen == 0) { // Error indicator
         continue;
      }
      if (iopos.fPos == pos && iopos.fLen == len) {
         found_request = kTRUE;
         resultBranch = b;
         // Note that we continue to iterate; fills up the rest of the entries in the cache.
      }
      // At this point, we are ready to push back a new offset
      fMissCache->fEntries.emplace_back(std::move(iopos));

      if (R__unlikely(perfStats)) {
         Int_t blistsize = b->GetWriteBasket();
         Int_t basketNumber = -1;
         for (Int_t bn = 0; bn < blistsize; ++bn) {
            if (iopos.fPos == b->GetBasketSeek(bn)) {
               basketNumber = bn;
               break;
            }
         }
         if (basketNumber >= 0)
            basketsInfo.emplace_back((size_t)i, basketNumber);
      }
   }
   if (R__unlikely(!found_request)) {
      // We have gone through all the branches in this file and the requested basket
      // doesn't appear to be in any of them.  Likely a logic error / bug.
      fMissCache->fEntries.clear();
   }
   if (R__unlikely(perfStats)) {
      for (auto &info : basketsInfo) {
         perfStats->SetLoadedMiss(info.first, info.second);
      }
   }
   return resultBranch;
}

////////////////////////////////////////////////////////////////////////////////
///
/// Process a cache miss; (pos, len) isn't in the buffer.
///
/// The first time we have a miss, we buffer as many baskets we can (up to the
/// maximum size of the TTreeCache) in memory from all branches that are not in
/// the prefetch list.
///
/// Subsequent times, we fetch all the buffers corresponding to branches that
/// had previously seen misses.  If it turns out the (pos, len) isn't in the
/// list of branches, we treat this as if it was the first miss.
///
/// Returns true if we were able to pull the data into the miss cache.
///
Bool_t TTreeCache::ProcessMiss(Long64_t pos, int len)
{

   Bool_t firstMiss = kFALSE;
   if (fFirstMiss == -1) {
      fFirstMiss = fEntryCurrent;
      firstMiss = kTRUE;
   }
   fLastMiss = fEntryCurrent;
   // The first time this is executed, we try to pull in as much data as we can.
   TBranch *b = CalculateMissEntries(pos, len, firstMiss);
   if (!b) {
      if (!firstMiss) {
         // TODO: this recalculates for *all* branches, throwing away the above work.
         b = CalculateMissEntries(pos, len, kTRUE);
      }
      if (!b) {
         // printf("ProcessMiss: pos %ld does not appear to correspond to a buffer in this file.\n", pos);
         // We have gone through all the branches in this file and the requested basket
         // doesn't appear to be in any of them.  Likely a logic error / bug.
         fMissCache->fEntries.clear();
         return kFALSE;
      }
   }
   // TODO: this should be a set.
   fMissCache->fBranches.push_back(b);

   // OK, sort the entries
   std::sort(fMissCache->fEntries.begin(), fMissCache->fEntries.end());

   // Now, fetch the buffer.
   std::vector<Long64_t> positions;
   positions.reserve(fMissCache->fEntries.size());
   std::vector<Int_t> lengths;
   lengths.reserve(fMissCache->fEntries.size());
   ULong64_t cumulative = 0;
   for (auto &mcentry : fMissCache->fEntries) {
      positions.push_back(mcentry.fIO.fPos);
      lengths.push_back(mcentry.fIO.fLen);
      mcentry.fIndex = cumulative;
      cumulative += mcentry.fIO.fLen;
   }
   fMissCache->fData.reserve(cumulative);
   // printf("Reading %lu bytes into miss cache for %lu entries.\n", cumulative, fEntries->size());
   fNMissReadPref += fMissCache->fEntries.size();
   fFile->ReadBuffers(&(fMissCache->fData[0]), &(positions[0]), &(lengths[0]), fMissCache->fEntries.size());
   fFirstMiss = fLastMiss = fEntryCurrent;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Given an IO operation (pos, len) that was a cache miss in the primary TTC,
/// try the operation again with the miss cache.
///
/// Returns true if the IO operation was successful and the contents of buf
/// were populated with the requested data.
///
Bool_t TTreeCache::CheckMissCache(char *buf, Long64_t pos, int len)
{

   if (!fOptimizeMisses) {
      return kFALSE;
   }
   if (R__unlikely((pos < 0) || (len < 0))) {
      return kFALSE;
   }

   // printf("Checking the miss cache for offset=%ld, length=%d\n", pos, len);

   // First, binary search to see if the desired basket is already cached.
   MissCache::Entry mcentry{IOPos{pos, len}};
   auto iter = std::lower_bound(fMissCache->fEntries.begin(), fMissCache->fEntries.end(), mcentry);

   if (iter != fMissCache->fEntries.end()) {
      if (len > iter->fIO.fLen) {
         ++fNMissReadMiss;
         return kFALSE;
      }
      auto offset = iter->fIndex;
      memcpy(buf, &(fMissCache->fData[offset]), len);
      // printf("Returning data from pos=%ld in miss cache.\n", offset);
      ++fNMissReadOk;
      return kTRUE;
   }

   // printf("Data not in miss cache.\n");

   // Update the cache, looking for this (pos, len).
   if (!ProcessMiss(pos, len)) {
      // printf("Unable to pull data into miss cache.\n");
      ++fNMissReadMiss;
      return kFALSE;
   }

   // OK, we updated the cache with as much information as possible.  Search again for
   // the entry we want.
   iter = std::lower_bound(fMissCache->fEntries.begin(), fMissCache->fEntries.end(), mcentry);

   if (iter != fMissCache->fEntries.end()) {
      auto offset = iter->fIndex;
      // printf("Expecting data at offset %ld in miss cache.\n", offset);
      memcpy(buf, &(fMissCache->fData[offset]), len);
      ++fNMissReadOk;
      return kTRUE;
   }

   // This must be a logic bug.  ProcessMiss should return false if (pos, len)
   // wasn't put into fEntries.
   ++fNMissReadMiss;
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// End of methods for miss cache.
////////////////////////////////////////////////////////////////////////////////

namespace {
struct BasketRanges {
   struct Range {
      Long64_t fMin; ///< Inclusive minimum
      Long64_t fMax; ///< Inclusive maximum

      Range() : fMin(-1), fMax(-1) {}

      void UpdateMin(Long64_t min)
      {
         if (fMin == -1 || min < fMin)
            fMin = min;
      }

      void UpdateMax(Long64_t max)
      {
         if (fMax == -1 || fMax < max)
            fMax = max;
      }

      Bool_t Contains(Long64_t entry) { return (fMin <= entry && entry <= fMax); }
   };

   std::vector<Range> fRanges;
   std::map<Long64_t,size_t> fMinimums;
   std::map<Long64_t,size_t> fMaximums;

   BasketRanges(size_t nBranches) { fRanges.resize(nBranches); }

   void Update(size_t branchNumber, Long64_t min, Long64_t max)
   {
      Range &range = fRanges.at(branchNumber);
      auto old(range);

      range.UpdateMin(min);
      range.UpdateMax(max);

      if (old.fMax != range.fMax) {
         if (old.fMax != -1) {
            auto maxIter = fMaximums.find(old.fMax);
            if (maxIter != fMaximums.end()) {
               if (maxIter->second == 1) {
                  fMaximums.erase(maxIter);
               } else {
                  --(maxIter->second);
               }
            }
         }
         ++(fMaximums[max]);
      }
   }

   void Update(size_t branchNumber, size_t basketNumber, Long64_t *entries, size_t nb, size_t max)
   {
      Update(branchNumber, entries[basketNumber],
             (basketNumber < (nb - 1)) ? (entries[basketNumber + 1] - 1) : max - 1);
   }

   // Check that fMaximums and fMinimums are properly set
   bool CheckAllIncludeRange()
   {
      Range result;
      for (const auto &r : fRanges) {
         if (result.fMin == -1 || result.fMin < r.fMin) {
            if (r.fMin != -1)
               result.fMin = r.fMin;
         }
         if (result.fMax == -1 || r.fMax < result.fMax) {
            if (r.fMax != -1)
               result.fMax = r.fMax;
         }
      }
      // if (result.fMax < result.fMin) {
      //    // No overlapping range.
      // }

      Range allIncludedRange(AllIncludedRange());

      return (result.fMin == allIncludedRange.fMin && result.fMax == allIncludedRange.fMax);
   }

   // This returns a Range object where fMin is the maximum of all the minimun entry
   // number loaded for each branch and fMax is the minimum of all the maximum entry
   // number loaded for each branch.
   // As such it is valid to have fMin > fMax, this is the case where there
   // are no overlap between the branch's range.  For example for 2 branches
   // where we have for one the entry [50,99] and for the other [0,49] then
   // we will have fMin = max(50,0) = 50 and fMax = min(99,49) = 49
   Range AllIncludedRange()
   {
      Range result;
      if (!fMinimums.empty())
         result.fMin = fMinimums.rbegin()->first;
      if (!fMaximums.empty())
         result.fMax = fMaximums.begin()->first;
      return result;
   }

   // Returns the number of branches with at least one baskets registered.
   UInt_t BranchesRegistered()
   {
      UInt_t result = 0;
      for (const auto &r : fRanges) {
         if (r.fMin != -1 && r.fMax != -1)
            ++result;
      }
      return result;
   }

   // Returns true if at least one of the branch's range contains
   // the entry.
   Bool_t Contains(Long64_t entry)
   {
      for (const auto &r : fRanges) {
         if (r.fMin != -1 && r.fMax != -1)
            if (r.fMin <= entry && entry <= r.fMax)
               return kTRUE;
      }
      return kFALSE;
   }

   void Print()
   {
      for (size_t i = 0; i < fRanges.size(); ++i) {
         if (fRanges[i].fMin != -1 || fRanges[i].fMax != -1)
            Printf("Range #%zu : %lld to %lld", i, fRanges[i].fMin, fRanges[i].fMax);
      }
   }
};
} // Anonymous namespace.

////////////////////////////////////////////////////////////////////////////////
/// Fill the cache buffer with the branches in the cache.

Bool_t TTreeCache::FillBuffer()
{

   if (fNbranches <= 0) return kFALSE;
   TTree *tree = ((TBranch*)fBranches->UncheckedAt(0))->GetTree();
   Long64_t entry = tree->GetReadEntry();
   Long64_t fEntryCurrentMax = 0;

   if (entry != -1 && (entry < fEntryMin || fEntryMax < entry))
      return kFALSE;

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
         // The learning phase should start from the minimum entry in the cache
         entry = fEntryMin;
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

   // Set to true to enable all debug output without having to set gDebug
   // Replace this once we have a per module and/or per class debugging level/setting.
   static constexpr bool showMore = kFALSE;

   static const auto PrintAllCacheInfo = [](TObjArray *branches) {
      for (Int_t i = 0; i < branches->GetEntries(); i++) {
         TBranch *b = (TBranch *)branches->UncheckedAt(i);
         b->PrintCacheInfo();
      }
   };

   if (showMore || gDebug > 6)
      Info("FillBuffer", "***** Called for entry %lld", entry);

   if (!fIsLearning && fEntryCurrent <= entry && entry < fEntryNext) {
      // Check if all the basket in the cache have already be used and
      // thus we can reuse the cache.
      Bool_t allUsed = kTRUE;
      for (Int_t i = 0; i < fNbranches; ++i) {
         TBranch *b = (TBranch *)fBranches->UncheckedAt(i);
         if (!b->fCacheInfo.AllUsed()) {
            allUsed = kFALSE;
            break;
         }
      }
      if (allUsed) {
         fEntryNext = entry;
         if (showMore || gDebug > 5)
            Info("FillBuffer", "All baskets used already, so refresh the cache early at entry %lld", entry);
      }
      if (gDebug > 8)
         PrintAllCacheInfo(fBranches);
   }

   // If the entry is in the range we previously prefetched, there is
   // no point in retrying.   Note that this will also return false
   // during the training phase (fEntryNext is then set intentional to
   // the end of the training phase).
   if (fEntryCurrent <= entry && entry < fEntryNext) return kFALSE;

   // Triggered by the user, not the learning phase
   if (entry == -1)
      entry = 0;

   Bool_t resetBranchInfo = kFALSE;
   if (entry < fCurrentClusterStart || fNextClusterStart <= entry) {
      // We are moving on to another set of clusters.
      resetBranchInfo = kTRUE;
      if (showMore || gDebug > 6)
         Info("FillBuffer", "*** Will reset the branch information about baskets");
   } else if (showMore || gDebug > 6) {
      Info("FillBuffer", "*** Info we have on the set of baskets");
      PrintAllCacheInfo(fBranches);
   }

   fEntryCurrentMax = fEntryCurrent;
   TTree::TClusterIterator clusterIter = tree->GetClusterIterator(entry);

   auto entryCurrent = clusterIter();
   auto entryNext    = clusterIter.GetNextEntry();

   if (entryNext < fEntryMin || fEntryMax < entryCurrent) {
      // There is no overlap between the cluster we found [entryCurrent, entryNext[
      // and the authorized range [fEntryMin, fEntryMax]
      // so we have nothing to do
      return kFALSE;
   }

   // If there is overlap between the found cluster and the authorized range
   // update the cache data members with the information about the current cluster.
   fEntryCurrent = entryCurrent;
   fEntryNext = entryNext;

   auto firstClusterEnd = fEntryNext;
   if (showMore || gDebug > 6)
      Info("FillBuffer", "Looking at cluster spanning from %lld to %lld", fEntryCurrent, fEntryNext);

   if (0 <= fNextClusterStart && fEntryCurrent < fNextClusterStart && !(fEntryCurrent < fNextClusterStart))
      fEntryCurrent = fNextClusterStart;
   if (fEntryCurrent < fEntryMin) fEntryCurrent = fEntryMin;
   if (fEntryMax <= 0) fEntryMax = tree->GetEntries();
   if (fEntryNext > fEntryMax) fEntryNext = fEntryMax;

   if ( fEnablePrefetching ) {
      if ( entry == fEntryMax ) {
         // We are at the end, no need to do anything else
         return kFALSE;
      }
   }

   if (resetBranchInfo) {
      // We earlier thought we were onto the next set of clusters.
      if (fCurrentClusterStart != -1 || fNextClusterStart != -1) {
         if (!(fEntryCurrent < fCurrentClusterStart || fEntryCurrent >= fNextClusterStart)) {
            Error("FillBuffer", "Inconsistency: fCurrentClusterStart=%lld fEntryCurrent=%lld fNextClusterStart=%lld "
                                "but fEntryCurrent should not be in between the two",
                  fCurrentClusterStart, fEntryCurrent, fNextClusterStart);
         }
      }

      // Start the next cluster set.
      fCurrentClusterStart = fEntryCurrent;
      fNextClusterStart = firstClusterEnd;
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
   Int_t ntotCurrentBuf = 0;
   if (fEnablePrefetching){ //prefetching mode
      if (fFirstBuffer) {
         TFileCacheRead::Prefetch(0,0);
         ntotCurrentBuf = fNtot;
      }
      else {
         TFileCacheRead::SecondPrefetch(0,0);
         ntotCurrentBuf = fBNtot;
      }
   }
   else {
      TFileCacheRead::Prefetch(0,0);
      ntotCurrentBuf = fNtot;
   }

   //store baskets
   BasketRanges ranges((showMore || gDebug > 6) ? fNbranches : 0);
   BasketRanges reqRanges(fNbranches);
   BasketRanges memRanges((showMore || gDebug > 6) ? fNbranches : 0);
   Int_t clusterIterations = 0;
   Long64_t minEntry = fEntryCurrent;
   Int_t prevNtot;
   Long64_t maxReadEntry = minEntry; // If we are stopped before the end of the 2nd pass, this marker will where we need to start next time.
   Int_t nReadPrefRequest = 0;
   auto perfStats = GetTree()->GetPerfStats();
   do {
      prevNtot = ntotCurrentBuf;
      Long64_t lowestMaxEntry = fEntryMax; // The lowest maximum entry in the TTreeCache for each branch for each pass.

      struct collectionInfo {
         Int_t fClusterStart{-1}; // First basket belonging to the current cluster
         Int_t fCurrent{0};       // Currently visited basket
         Bool_t fLoadedOnce{kFALSE};

         void Rewind() { fCurrent = (fClusterStart >= 0) ? fClusterStart : 0; }
      };
      std::vector<collectionInfo> cursor(fNbranches);
      Bool_t reachedEnd = kFALSE;
      Bool_t skippedFirst = kFALSE;
      Bool_t oncePerBranch = kFALSE;
      Int_t nDistinctLoad = 0;
      Bool_t progress = kTRUE;
      enum ENarrow {
         kFull = 0,
         kNarrow = 1
      };
      enum EPass {
         kStart = 1,
         kRegular = 2,
         kRewind = 3
      };

      auto CollectBaskets = [this, elist, chainOffset, entry, clusterIterations, resetBranchInfo, perfStats,
       &cursor, &lowestMaxEntry, &maxReadEntry, &minEntry,
       &reachedEnd, &skippedFirst, &oncePerBranch, &nDistinctLoad, &progress,
       &ranges, &memRanges, &reqRanges,
       &ntotCurrentBuf, &nReadPrefRequest](EPass pass, ENarrow narrow, Long64_t maxCollectEntry) {
         // The first pass we add one basket per branches around the requested entry
         // then in the second pass we add the other baskets of the cluster.
         // This is to support the case where the cache is too small to hold a full cluster.
         Int_t nReachedEnd = 0;
         Int_t nSkipped = 0;
         auto oldnReadPrefRequest = nReadPrefRequest;
         std::vector<Int_t> potentialVetoes;

         if (showMore || gDebug > 7)
            Info("CollectBaskets", "Called with pass=%d narrow=%d maxCollectEntry=%lld", pass, narrow, maxCollectEntry);

         Bool_t filled = kFALSE;
         for (Int_t i = 0; i < fNbranches; ++i) {
            TBranch *b = (TBranch*)fBranches->UncheckedAt(i);
            if (b->GetDirectory()==0 || b->TestBit(TBranch::kDoNotProcess))
               continue;
            if (b->GetDirectory()->GetFile() != fFile)
               continue;
            potentialVetoes.clear();
            if (pass == kStart && !cursor[i].fLoadedOnce && resetBranchInfo) {
               // First check if we have any cluster that is currently in the
               // cache but was not used and would be reloaded in the next
               // cluster.
               b->fCacheInfo.GetUnused(potentialVetoes);
               if (showMore || gDebug > 7) {
                  TString vetolist;
                  for(auto v : potentialVetoes) {
                     vetolist += v;
                     vetolist.Append(' ');
                  }
                  if (!potentialVetoes.empty())
                     Info("FillBuffer", "*** Potential Vetos for branch #%d: %s", i, vetolist.Data());
               }
               b->fCacheInfo.Reset();
            }
            Int_t nb = b->GetMaxBaskets();
            Int_t *lbaskets   = b->GetBasketBytes();
            Long64_t *entries = b->GetBasketEntry();
            if (!lbaskets || !entries)
               continue;
            //we have found the branch. We now register all its baskets
            // from the requested offset to the basket below fEntryMax
            Int_t blistsize = b->GetListOfBaskets()->GetSize();

            auto maxOfBasket = [this, nb, entries](int j) {
               return ((j < (nb - 1)) ? (entries[j + 1] - 1) : fEntryMax - 1);
            };

            if (pass == kRewind)
               cursor[i].Rewind();
            for (auto &j = cursor[i].fCurrent; j < nb; j++) {
               // This basket has already been read, skip it

               if (j < blistsize && b->GetListOfBaskets()->UncheckedAt(j)) {

                  if (showMore || gDebug > 6) {
                     ranges.Update(i, entries[j], maxOfBasket(j));
                     memRanges.Update(i, entries[j], maxOfBasket(j));
                  }
                  if (entries[j] <= entry && entry <= maxOfBasket(j)) {
                     b->fCacheInfo.SetIsInCache(j);
                     b->fCacheInfo.SetUsed(j);
                     if (narrow) {
                        // In narrow mode, we would select 'only' this basket,
                        // so we are done for this round, let's 'consume' this
                        // basket and go.
                        ++nReachedEnd;
                        ++j;
                        break;
                     }
                  }
                  continue;
               }

               // Important: do not try to read maxCollectEntry, otherwise we might jump to the next autoflush
               if (entries[j] >= maxCollectEntry) {
                  ++nReachedEnd;
                  break; // break out of the for each branch loop.
               }

               Long64_t pos = b->GetBasketSeek(j);
               Int_t len = lbaskets[j];
               if (pos <= 0 || len <= 0)
                  continue;
               if (len > fBufferSizeMin) {
                  // Do not cache a basket if it is bigger than the cache size!
                  if ((showMore || gDebug > 7) &&
                      (!(entries[j] < minEntry && (j < nb - 1 && entries[j + 1] <= minEntry))))
                     Info("FillBuffer", "Skipping branch %s basket %d is too large for the cache: %d > %d",
                          b->GetName(), j, len, fBufferSizeMin);
                  continue;
               }

               if (nReadPrefRequest && entries[j] > (reqRanges.AllIncludedRange().fMax + 1)) {
                  // There is a gap between this basket and the max of the 'lowest' already loaded basket
                  // If we are tight in memory, reading this basket may prevent reading the basket (for the other branches)
                  // that covers this gap, forcing those baskets to be read uncached (because the cache wont be reloaded
                  // until we use this basket).
                  // eg. We could end up with the cache contain
                  //   b1: [428, 514[ // 'this' basket and we can assume [321 to 428[ is already in memory
                  //   b2: [400, 424[
                  // and when reading entry 425 we will read b2's basket uncached.

                  if (showMore || gDebug > 8)
                     Info("FillBuffer", "Skipping for now due to gap %d/%d with %lld > %lld", i, j, entries[j],
                          (reqRanges.AllIncludedRange().fMax + 1));
                  break; // Without consuming the basket.
               }

               if (entries[j] < minEntry && (j<nb-1 && entries[j+1] <= minEntry))
                  continue;

               // We are within the range
               if (cursor[i].fClusterStart == -1)
                  cursor[i].fClusterStart = j;

               if (elist) {
                  Long64_t emax = fEntryMax;
                  if (j<nb-1)
                     emax = entries[j + 1] - 1;
                  if (!elist->ContainsRange(entries[j]+chainOffset,emax+chainOffset))
                     continue;
               }

               if (b->fCacheInfo.HasBeenUsed(j) || b->fCacheInfo.IsInCache(j) || b->fCacheInfo.IsVetoed(j)) {
                  // We already cached and used this basket during this cluster range,
                  // let's not redo it
                  if (showMore || gDebug > 7)
                     Info("FillBuffer", "Skipping basket to avoid redo: %d/%d veto: %d", i, j, b->fCacheInfo.IsVetoed(j));
                  continue;
               }

               if (std::find(std::begin(potentialVetoes), std::end(potentialVetoes), j) != std::end(potentialVetoes)) {
                  // This basket was in the previous cache/cluster and was not used,
                  // let's not read it again. I.e. we bet that it will continue to not
                  // be used.  At worst it will be used and thus read by itself.
                  // Usually in this situation the basket is large so the penalty for
                  // (re)reading it uselessly is high and the penalty to read it by
                  // itself is 'small' (i.e. size bigger than latency).
                  b->fCacheInfo.Veto(j);
                  if (showMore || gDebug > 7)
                     Info("FillBuffer", "Veto-ing cluster %d [%lld,%lld[ in branch %s #%d", j, entries[j],
                          maxOfBasket(j) + 1, b->GetName(), i);
                  continue;
               }

               if (narrow) {
                  if ((((entries[j] > entry)) || (j < nb - 1 && entries[j + 1] <= entry))) {
                     // Keep only the basket that contains the entry
                     if (j == cursor[i].fClusterStart && entry > entries[j])
                        ++nSkipped;
                     if (entries[j] > entry)
                        break;
                     else
                        continue;
                  }
               }

               if ((ntotCurrentBuf + len) > fBufferSizeMin) {
                  // Humm ... we are going to go over the requested size.
                  if (clusterIterations > 0 && cursor[i].fLoadedOnce) {
                     // We already have a full cluster and now we would go over the requested
                     // size, let's stop caching (and make sure we start next time from the
                     // end of the previous cluster).
                     if (showMore || gDebug > 5) {
                        Info(
                           "FillBuffer",
                           "Breaking early because %d is greater than %d at cluster iteration %d will restart at %lld",
                           (ntotCurrentBuf + len), fBufferSizeMin, clusterIterations, minEntry);
                     }
                     fEntryNext = minEntry;
                     filled = kTRUE;
                     break;
                  } else {
                     if (pass == kStart || !cursor[i].fLoadedOnce) {
                        if ((ntotCurrentBuf + len) > 4 * fBufferSizeMin) {
                           // Okay, so we have not even made one pass and we already have
                           // accumulated request for more than twice the memory size ...
                           // So stop for now, and will restart at the same point, hoping
                           // that the basket will still be in memory and not asked again ..
                           fEntryNext = maxReadEntry;

                           if (showMore || gDebug > 5) {
                              Info("FillBuffer", "Breaking early because %d is greater than 4*%d at cluster iteration "
                                                 "%d pass %d will restart at %lld",
                                   (ntotCurrentBuf + len), fBufferSizeMin, clusterIterations, pass, fEntryNext);
                           }
                           filled = kTRUE;
                           break;
                        }
                     } else {
                        // We have made one pass through the branches and thus already
                        // requested one basket per branch, let's stop prefetching
                        // now.
                        if ((ntotCurrentBuf + len) > 2 * fBufferSizeMin) {
                           fEntryNext = maxReadEntry;
                           if (showMore || gDebug > 5) {
                              Info("FillBuffer", "Breaking early because %d is greater than 2*%d at cluster iteration "
                                                 "%d pass %d will restart at %lld",
                                   (ntotCurrentBuf + len), fBufferSizeMin, clusterIterations, pass, fEntryNext);
                           }
                           filled = kTRUE;
                           break;
                        }
                     }
                  }
               }

               ++nReadPrefRequest;

               reqRanges.Update(i, j, entries, nb, fEntryMax);
               if (showMore || gDebug > 6)
                  ranges.Update(i, j, entries, nb, fEntryMax);

               b->fCacheInfo.SetIsInCache(j);

               if (showMore || gDebug > 6)
                  Info("FillBuffer", "*** Registering branch %d basket %d %s", i, j, b->GetName());

               if (!cursor[i].fLoadedOnce) {
                  cursor[i].fLoadedOnce = kTRUE;
                  ++nDistinctLoad;
               }
               if (R__unlikely(perfStats)) {
                  perfStats->SetLoaded(i, j);
               }

               // Actual registering the basket for loading from the file.
               if (fEnablePrefetching){
                  if (fFirstBuffer) {
                     TFileCacheRead::Prefetch(pos,len);
                     ntotCurrentBuf = fNtot;
                  }
                  else {
                     TFileCacheRead::SecondPrefetch(pos,len);
                     ntotCurrentBuf = fBNtot;
                  }
               }
               else {
                  TFileCacheRead::Prefetch(pos,len);
                  ntotCurrentBuf = fNtot;
               }

               if ( ( j < (nb-1) ) && entries[j+1] > maxReadEntry ) {
                  // Info("FillBuffer","maxCollectEntry incremented from %lld to %lld", maxReadEntry, entries[j+1]);
                  maxReadEntry = entries[j+1];
               }
               if (ntotCurrentBuf > 4 * fBufferSizeMin) {
                  // Humm something wrong happened.
                  Warning("FillBuffer", "There is more data in this cluster (starting at entry %lld to %lld, "
                                        "current=%lld) than usual ... with %d %.3f%% of the branches we already have "
                                        "%d bytes (instead of %d)",
                          fEntryCurrent, fEntryNext, entries[j], i, (100.0 * i) / ((float)fNbranches), ntotCurrentBuf,
                          fBufferSizeMin);
               }
               if (pass == kStart) {
                  // In the first pass, we record one basket per branch and move on to the next branch.
                  auto high = maxOfBasket(j);
                  if (high < lowestMaxEntry)
                     lowestMaxEntry = high;
                  // 'Consume' the baskets (i.e. avoid looking at it during a subsequent pass)
                  ++j;
                  break;
               } else if ((j + 1) == nb || entries[j + 1] >= maxReadEntry || entries[j + 1] >= lowestMaxEntry) {
                  // In the other pass, load the baskets until we get to the maximum loaded so far.
                  auto high = maxOfBasket(j);
                  if (high < lowestMaxEntry)
                     lowestMaxEntry = high;
                  // 'Consume' the baskets (i.e. avoid looking at it during a subsequent pass)
                  ++j;
                  break;
               }
            }

            if (cursor[i].fCurrent == nb) {
               ++nReachedEnd;
            }

            if (gDebug > 0)
               Info("CollectBaskets",
                    "Entry: %lld, registering baskets branch %s, fEntryNext=%lld, fNseek=%d, ntotCurrentBuf=%d",
                    minEntry, ((TBranch *)fBranches->UncheckedAt(i))->GetName(), fEntryNext, fNseek, ntotCurrentBuf);
         }
         reachedEnd = (nReachedEnd == fNbranches);
         skippedFirst = (nSkipped > 0);
         oncePerBranch = (nDistinctLoad == fNbranches);
         progress = nReadPrefRequest - oldnReadPrefRequest;
         return filled;
      };

      // First collect all the basket containing the request entry.
      bool full = kFALSE;

      full = CollectBaskets(kStart, kNarrow, fEntryNext);

      // Then fill out from all but the 'largest' branch to even out
      // the range across branches;
      while (!full && !reachedEnd && progress) { // used to be restricted to !oncePerBranch
         full = CollectBaskets(kStart, kFull, std::min(maxReadEntry, fEntryNext));
      }

      resetBranchInfo = kFALSE; // Make sure the 2nd cluster iteration does not erase the info.

      // Then fill out to the end of the cluster.
      if (!full && !fReverseRead) {
         do {
            full = CollectBaskets(kRegular, kFull, fEntryNext);
         } while (!full && !reachedEnd && progress);
      }

      // The restart from the start of the cluster.
      if (!full && skippedFirst) {
         full = CollectBaskets(kRewind, kFull, fEntryNext);
         while (!full && !reachedEnd && progress) {
            full = CollectBaskets(kRegular, kFull, fEntryNext);
         }
      }

      clusterIterations++;

      minEntry = clusterIter.Next();
      if (fIsLearning) {
         fFillTimes++;
      }

      // Continue as long as we still make progress (prevNtot < ntotCurrentBuf), that the next entry range to be looked
      // at,
      // which start at 'minEntry', is not past the end of the requested range (minEntry < fEntryMax)
      // and we guess that we not going to go over the requested amount of memory by asking for another set
      // of entries (fBufferSizeMin > ((Long64_t)ntotCurrentBuf*(clusterIterations+1))/clusterIterations).
      // ntotCurrentBuf / clusterIterations is the average size we are accumulated so far at each loop.
      // and thus (ntotCurrentBuf / clusterIterations) * (clusterIterations+1) is a good guess at what the next total
      // size
      // would be if we run the loop one more time.   ntotCurrentBuf and clusterIterations are Int_t but can sometimes
      // be 'large' (i.e. 30Mb * 300 intervals) and can overflow the numerical limit of Int_t (i.e. become
      // artificially negative).   To avoid this issue we promote ntotCurrentBuf to a long long (64 bits rather than 32
      // bits)
      if (!((fBufferSizeMin > ((Long64_t)ntotCurrentBuf * (clusterIterations + 1)) / clusterIterations) &&
            (prevNtot < ntotCurrentBuf) && (minEntry < fEntryMax))) {
         if (showMore || gDebug > 6)
            Info("FillBuffer", "Breaking because %d <= %lld || (%d >= %d) || %lld >= %lld", fBufferSizeMin,
                 ((Long64_t)ntotCurrentBuf * (clusterIterations + 1)) / clusterIterations, prevNtot, ntotCurrentBuf,
                 minEntry, fEntryMax);
         break;
      }

      //for the reverse reading case
      if (!fIsLearning && fReverseRead) {
         if (clusterIterations >= fFillTimes)
            break;
         if (minEntry >= fEntryCurrentMax && fEntryCurrentMax > 0)
            break;
      }
      fEntryNext = clusterIter.GetNextEntry();
      if (fEntryNext > fEntryMax) fEntryNext = fEntryMax;
      fNextClusterStart = fEntryNext;
   } while (kTRUE);

   if (showMore || gDebug > 6) {
      Info("FillBuffer", "Mem ranges");
      memRanges.Print();
      Info("FillBuffer", "Combined ranges");
      ranges.Print();
      Info("FillBuffer", "Requested ranges");
      reqRanges.Print();
      PrintAllCacheInfo(fBranches);
   }

   if (nReadPrefRequest == 0) {
      // Nothing was added in the cache.  This usually indicates that the baskets
      // contains the requested entry are either already in memory or are too large
      // on their own to fit in the cache.
      if (showMore || gDebug > 5) {
         Info("FillBuffer", "For entry %lld, nothing was added to the cache.", entry);
      }
   } else if (fEntryNext < firstClusterEnd && !reqRanges.Contains(entry)) {
      // Something went very wrong and even-though we searched for the baskets
      // holding 'entry' we somehow ended up with a range of entries that does
      // validate.  So we must have been unable to find or fit the needed basket.
      // And thus even-though, we know the corresponding baskets wont be in the cache,
      // Let's make it official that 'entry' is within the range of this TTreeCache ('s search.)

      // Without this, the next read will be flagged as 'out-of-range' and then we start at
      // the exact same point as this FillBuffer execution resulting in both the requested
      // entry still not being part of the cache **and** the beginning of the cluster being
      // read **again**.

      if (showMore || gDebug > 5) {
         Error("FillBuffer", "Reset the next entry because the currently loaded range does not contains the request "
                             "entry: %lld.  fEntryNext updated from %lld to %lld. %d",
               entry, fEntryNext, firstClusterEnd, nReadPrefRequest);
         reqRanges.Print();
      }

      fEntryNext = firstClusterEnd;
   } else {
      if (showMore || gDebug > 5) {
         Info("FillBuffer", "Complete adding %d baskets from %d branches taking in memory %d out of %d",
              nReadPrefRequest, reqRanges.BranchesRegistered(), ntotCurrentBuf, fBufferSizeMin);
      }
   }

   fNReadPref += nReadPrefRequest;
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

////////////////////////////////////////////////////////////////////////////////
/// The total efficiency of the 'miss cache' - defined as the ratio
/// of blocks found in the cache versus the number of blocks prefetched

Double_t TTreeCache::GetMissEfficiency() const
{
   if (!fNMissReadPref) {
      return 0;
   }
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

////////////////////////////////////////////////////////////////////////////////
/// Relative efficiency of the 'miss cache' - ratio of the reads found in cache
/// to the number of reads so far.

Double_t TTreeCache::GetMissEfficiencyRel() const
{
   if (!fNMissReadOk && !fNMissReadMiss) {
      return 0;
   }

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

   static const auto recordMiss = [](TVirtualPerfStats *perfStats, TObjArray *branches, Bool_t bufferFilled,
                                     Long64_t basketpos) {
      if (gDebug > 6)
         ::Info("TTreeCache::ReadBufferNormal", "Cache miss after an %s FillBuffer: pos=%lld",
                bufferFilled ? "active" : "inactive", basketpos);
      for (Int_t i = 0; i < branches->GetEntries(); ++i) {
         TBranch *b = (TBranch *)branches->UncheckedAt(i);
         Int_t blistsize = b->GetListOfBaskets()->GetSize();
         for (Int_t j = 0; j < blistsize; ++j) {
            if (basketpos == b->GetBasketSeek(j)) {
               if (gDebug > 6)
                  ::Info("TTreeCache::ReadBufferNormal", "   Missing basket: %d for %s", j, b->GetName());
               perfStats->SetMissed(i, j);
            }
         }
      }
   };

   //not found in cache. Do we need to fill the cache?
   Bool_t bufferFilled = FillBuffer();
   if (bufferFilled) {
      Int_t res = TFileCacheRead::ReadBuffer(buf,pos,len);

      if (res == 1)
         fNReadOk++;
      else if (res == 0) {
         fNReadMiss++;
         auto perfStats = GetTree()->GetPerfStats();
         if (perfStats)
            recordMiss(perfStats, fBranches, bufferFilled, pos);
      }

      return res;
   }

   if (CheckMissCache(buf, pos, len)) {
      return 1;
   }

   fNReadMiss++;
   auto perfStats = GetTree()->GetPerfStats();
   if (perfStats)
      recordMiss(perfStats, fBranches, bufferFilled, pos);

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Used to read a chunk from a block previously fetched. It will call FillBuffer
/// even if the cache lookup succeeds, because it will try to prefetch the next block
/// as soon as we start reading from the current block.

Int_t TTreeCache::ReadBufferPrefetch(char *buf, Long64_t pos, Int_t len)
{
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
   for (Int_t i = 0; i < fNbranches; ++i) {
      TBranch *b = (TBranch*)fBranches->UncheckedAt(i);
      if (b->GetDirectory()==0 || b->TestBit(TBranch::kDoNotProcess))
         continue;
      if (b->GetDirectory()->GetFile() != fFile)
         continue;
      b->fCacheInfo.Reset();
   }
   fEntryCurrent = -1;
   fEntryNext = -1;
   fCurrentClusterStart = -1;
   fNextClusterStart = -1;

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

   auto perfStats = GetTree()->GetPerfStats();
   if (perfStats)
      perfStats->UpdateBranchIndices(fBranches);

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

   auto perfStats = GetTree()->GetPerfStats();
   if (perfStats)
      perfStats->UpdateBranchIndices(fBranches);
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

   Long64_t entry = fTree ? fTree->GetReadEntry() : 0;

   // Return early if we are out of the requested range.
   if (entry < fEntryMin || entry > fEntryMax) return;

   fLearnPrefilling = kTRUE;


   // Force only the learn entries to be cached by temporarily setting min/max
   // to the learning phase entry range
   // But save all the old values, so we can restore everything to how it was
   Long64_t eminOld = fEntryMin;
   Long64_t emaxOld = fEntryMax;
   Long64_t ecurrentOld = fEntryCurrent;
   Long64_t enextOld = fEntryNext;
   auto currentClusterStartOld = fCurrentClusterStart;
   auto nextClusterStartOld = fNextClusterStart;

   fEntryMin = std::max(fEntryMin, fEntryCurrent);
   fEntryMax = std::min(fEntryMax, fEntryNext);

   // We check earlier that we are within the authorized range, but
   // we might still be out of the (default) learning range and since
   // this is called before any branch is added to the cache, this means
   // that the user's first GetEntry is this one which is outside of the
   // learning range ... so let's do something sensible-ish.
   // Note: we probably should also fix the learning range but we may
   // or may not have enough information to know if we can move it
   // (for example fEntryMin (eminOld right now) might be the default or user provided)
   if (entry < fEntryMin) fEntryMin = entry;
   if (entry > fEntryMax) fEntryMax = entry;

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
   fCurrentClusterStart = currentClusterStartOld;
   fNextClusterStart = nextClusterStartOld;

   fLearnPrefilling = kFALSE;
}
