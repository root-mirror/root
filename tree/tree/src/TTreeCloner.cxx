// @(#)root/tree:$Id$
// Author: Philippe Canal 07/11/2005

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TTreeCloner
\ingroup tree

Class implementing or helping  the various TTree cloning method
*/

#include "TBasket.h"
#include "TBranch.h"
#include "TBranchClones.h"
#include "TBranchElement.h"
#include "TStreamerInfo.h"
#include "TBranchRef.h"
#include "TError.h"
#include "TProcessID.h"
#include "TTree.h"
#include "TTreeCloner.h"
#include "TFile.h"
#include "TLeafB.h"
#include "TLeafI.h"
#include "TLeafL.h"
#include "TLeafS.h"
#include "TLeafO.h"
#include "TLeafC.h"
#include "TFileCacheRead.h"
#include "TTreeCache.h"
#include "snprintf.h"

#include <algorithm>

////////////////////////////////////////////////////////////////////////////////

Bool_t TTreeCloner::CompareSeek::operator()(UInt_t i1, UInt_t i2)
{
   if (fObject->fBasketSeek[i1] ==  fObject->fBasketSeek[i2]) {
      if (fObject->fBasketEntry[i1] ==  fObject->fBasketEntry[i2]) {
         return i1 < i2;
      }
      return  fObject->fBasketEntry[i1] <  fObject->fBasketEntry[i2];
   }
   return fObject->fBasketSeek[i1] <  fObject->fBasketSeek[i2];
}

////////////////////////////////////////////////////////////////////////////////

Bool_t TTreeCloner::CompareEntry::operator()(UInt_t i1, UInt_t i2)
{
   if (fObject->fBasketEntry[i1] ==  fObject->fBasketEntry[i2]) {
      return i1 < i2;
   }
   return  fObject->fBasketEntry[i1] <  fObject->fBasketEntry[i2];
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.  This object would transfer the data from
/// 'from' to 'to' using the method indicated in method.
///
/// The value of the parameter 'method' determines in which
/// order the branches' baskets are written to the output file.
///
/// When a TTree is filled the data is stored in the individual
/// branches' basket.  Each basket is written individually to
/// the disk as soon as it is full.  In consequence the baskets
/// of branches that contain 'large' data chunk are written to
/// the disk more often.
///
/// There is currently 3 supported sorting order:
///
///     SortBasketsByOffset (the default)
///     SortBasketsByBranch
///     SortBasketsByEntry
///
/// When using SortBasketsByOffset the baskets are written in
/// the output file in the same order as in the original file
/// (i.e. the basket are sorted on their offset in the original
/// file; Usually this also means that the baskets are sorted
/// on the index/number of the _last_ entry they contain)
///
/// When using SortBasketsByBranch all the baskets of each
/// individual branches are stored contiguously.  This tends to
/// optimize reading speed when reading a small number (1->5) of
/// branches, since all their baskets will be clustered together
/// instead of being spread across the file.  However it might
/// decrease the performance when reading more branches (or the full
/// entry).
///
/// When using SortBasketsByEntry the baskets with the lowest
/// starting entry are written first.  (i.e. the baskets are
/// sorted on the index/number of the first entry they contain).
/// This means that on the file the baskets will be in the order
/// in which they will be needed when reading the whole tree
/// sequentially.

TTreeCloner::TTreeCloner(TTree *from, TTree *to, Option_t *method, UInt_t options) :
   fWarningMsg(),
   fIsValid(kTRUE),
   fNeedConversion(kFALSE),
   fOptions(options),
   fFromTree(from),
   fToTree(to),
   fMethod(method),
   fFromBranches( from ? from->GetListOfLeaves()->GetEntries()+1 : 0),
   fToBranches( to ? to->GetListOfLeaves()->GetEntries()+1 : 0),
   fMaxBaskets(CollectBranches()),
   fBasketBranchNum(new UInt_t[fMaxBaskets]),
   fBasketNum(new UInt_t[fMaxBaskets]),
   fBasketSeek(new Long64_t[fMaxBaskets]),
   fBasketEntry(new Long64_t[fMaxBaskets]),
   fBasketIndex(new UInt_t[fMaxBaskets]),
   fPidOffset(0),
   fCloneMethod(TTreeCloner::kDefault),
   fToStartEntries(0),
   fCacheSize(0LL),
   fFileCache(nullptr),
   fPrevCache(nullptr)
{
   TString opt(method);
   opt.ToLower();
   if (opt.Contains("sortbasketsbybranch")) {
      //::Info("TTreeCloner::TTreeCloner","use: kSortBasketsByBranch");
      fCloneMethod = TTreeCloner::kSortBasketsByBranch;
   } else if (opt.Contains("sortbasketsbyentry")) {
      //::Info("TTreeCloner::TTreeCloner","use: kSortBasketsByEntry");
      fCloneMethod = TTreeCloner::kSortBasketsByEntry;
   } else {
      //::Info("TTreeCloner::TTreeCloner","use: kSortBasketsByOffset");
      fCloneMethod = TTreeCloner::kSortBasketsByOffset;
   }
   if (fToTree) fToStartEntries = fToTree->GetEntries();

   if (fFromTree == nullptr) {
      if (to)
         fWarningMsg.Form("An input TTree is required (cloning to %s).",
                          to->GetName());
      else
         fWarningMsg.Form("An input and output TTree are required.");
      if (!(fOptions & kNoWarnings)) {
         Warning("TTreeCloner::TTreeCloner", "%s", fWarningMsg.Data());
      }
      fIsValid = kFALSE;
   }
   if (fToTree == nullptr) {
      fWarningMsg.Form("An output TTree is required (cloning %s).",
                       from ? from->GetName() : "no tree");
      if (!(fOptions & kNoWarnings)) {
         Warning("TTreeCloner::TTreeCloner", "%s", fWarningMsg.Data());
      }
      fIsValid = kFALSE;
   } else if (fToTree->GetDirectory() == nullptr) {
      fWarningMsg.Form("The output TTree (%s) must be associated with a directory.",
                       fToTree->GetName());
      if (!(fOptions & kNoWarnings)) {
         Warning("TTreeCloner::TTreeCloner", "%s", fWarningMsg.Data());
      }
      fIsValid = kFALSE;
   } else if (fToTree->GetCurrentFile() == nullptr) {
      fWarningMsg.Form("The output TTree (%s) must be associated with a directory (%s) that is in a file.",
                       fToTree->GetName(),fToTree->GetDirectory()->GetName());
      if (!(fOptions & kNoWarnings)) {
         Warning("TTreeCloner::TTreeCloner", "%s", fWarningMsg.Data());
      }
      fIsValid = kFALSE;
   } else if (! fToTree->GetDirectory()->IsWritable()) {
      if (fToTree->GetDirectory()==fToTree->GetCurrentFile()) {
         fWarningMsg.Form("The output TTree (%s) must be associated with a writable file (%s).",
                          fToTree->GetName(),fToTree->GetCurrentFile()->GetName());
      } else {
         fWarningMsg.Form("The output TTree (%s) must be associated with a writable directory (%s in %s).",
                          fToTree->GetName(),fToTree->GetDirectory()->GetName(),fToTree->GetCurrentFile()->GetName());
      }
      if (!(fOptions & kNoWarnings)) {
         Warning("TTreeCloner::TTreeCloner", "%s", fWarningMsg.Data());
      }
      fIsValid = kFALSE;
   }

   if (fIsValid && (!(fOptions & kNoFileCache))) {
      fCacheSize = fFromTree->GetCacheAutoSize();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Execute the cloning.

Bool_t TTreeCloner::Exec()
{
   if (!IsValid()) {
      return kFALSE;
   }
   CreateCache();
   ImportClusterRanges();
   CopyStreamerInfos();
   CopyProcessIds();
   CloseOutWriteBaskets();
   CollectBaskets();
   SortBaskets();
   WriteBaskets();
   CopyMemoryBaskets();
   RestoreCache();

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// TTreeCloner destructor

TTreeCloner::~TTreeCloner()
{
   // The file cache was restored to its previous value at the end of Exec,
   // we can safely delete our cache.
   delete fFileCache;
   delete [] fBasketBranchNum;
   delete [] fBasketNum;
   delete [] fBasketSeek;
   delete [] fBasketEntry;
   delete [] fBasketIndex;
}

////////////////////////////////////////////////////////////////////////////////
/// Before we can start adding new basket, we need to flush to
/// disk the partially filled baskets (the WriteBasket)

void TTreeCloner::CloseOutWriteBaskets()
{
   for(Int_t i=0; i<fToBranches.GetEntries(); ++i) {
      TBranch *to = (TBranch*)fToBranches.UncheckedAt(i);
      to->FlushOneBasket(to->GetWriteBasket());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Fill the array of branches, adding the branch 'from' and 'to',
/// and matching the sub-branches of the 'from' and 'to' branches.
/// Returns the total number of baskets in all the from branch and
/// it sub-branches.

UInt_t TTreeCloner::CollectBranches(TBranch *from, TBranch *to) {
   // Since this is called from the constructor, this can not be a virtual function

   UInt_t numBaskets = 0;
   if (from->InheritsFrom(TBranchClones::Class())) {
      TBranchClones *fromclones = (TBranchClones*) from;
      TBranchClones *toclones = (TBranchClones*) to;
      numBaskets += CollectBranches(fromclones->fBranchCount, toclones->fBranchCount);

   } else if (from->InheritsFrom(TBranchElement::Class())) {
      Int_t nb = from->GetListOfLeaves()->GetEntries();
      Int_t fnb = to->GetListOfLeaves()->GetEntries();
      if (nb != fnb && (nb == 0 || fnb == 0)) {
         // We might be in the case where one branch is split
         // while the other is not split.  We must reject this match.
         fWarningMsg.Form("The export branch and the import branch do not have the same split level. (The branch name is %s.)",
                          from->GetName());
         if (!(fOptions & kNoWarnings)) {
            Warning("TTreeCloner::CollectBranches", "%s", fWarningMsg.Data());
         }
         fNeedConversion = kTRUE;
         fIsValid = kFALSE;
         return 0;
      }
      if (((TBranchElement*) from)->GetStreamerType() != ((TBranchElement*) to)->GetStreamerType()) {
         fWarningMsg.Form("The export branch and the import branch do not have the same streamer type. (The branch name is %s.)",
                          from->GetName());
         if (!(fOptions & kNoWarnings)) {
            Warning("TTreeCloner::CollectBranches", "%s", fWarningMsg.Data());
         }
         fIsValid = kFALSE;
         return 0;
      }
      TBranchElement *fromelem = (TBranchElement*) from;
      TBranchElement *toelem = (TBranchElement*) to;
      if (fromelem->fMaximum > toelem->fMaximum) toelem->fMaximum = fromelem->fMaximum;
   } else {

      Int_t nb = from->GetListOfLeaves()->GetEntries();
      Int_t fnb = to->GetListOfLeaves()->GetEntries();
      if (nb != fnb) {
         fWarningMsg.Form("The export branch and the import branch (%s) do not have the same number of leaves (%d vs %d)",
                          from->GetName(), fnb, nb);
         if (!(fOptions & kNoWarnings)) {
            Error("TTreeCloner::CollectBranches", "%s", fWarningMsg.Data());
         }
         fIsValid = kFALSE;
         return 0;
      }
      for (Int_t i=0;i<nb;i++)  {

         TLeaf *fromleaf = (TLeaf*)from->GetListOfLeaves()->At(i);
         TLeaf *toleaf = (TLeaf*)to->GetListOfLeaves()->At(i);
         if (toleaf->IsA() != fromleaf->IsA() ) {
            // The data type do not match, we can not do a fast merge.
            fWarningMsg.Form("The export leaf and the import leaf (%s.%s) do not have the same data type (%s vs %s)",
                              from->GetName(),fromleaf->GetName(),fromleaf->GetTypeName(),toleaf->GetTypeName());
            if (! (fOptions & kNoWarnings) ) {
               Warning("TTreeCloner::CollectBranches", "%s", fWarningMsg.Data());
            }
            fIsValid = kFALSE;
            fNeedConversion = kTRUE;
            return 0;
         }
         toleaf->IncludeRange( fromleaf );
      }

   }

   fFromBranches.AddLast(from);
   if (!from->TestBit(TBranch::kDoNotUseBufferMap)) {
      // Make sure that we reset the Buffer's map if needed.
      to->ResetBit(TBranch::kDoNotUseBufferMap);
   }
   fToBranches.AddLast(to);

   numBaskets += from->GetWriteBasket();
   numBaskets += CollectBranches(from->GetListOfBranches(),to->GetListOfBranches());

   return numBaskets;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill the array of branches, matching the branches of the 'from' and 'to' arrays.
/// Returns the total number of baskets in all the branches.

UInt_t TTreeCloner::CollectBranches(TObjArray *from, TObjArray *to)
{
   // Since this is called from the constructor, this can not be a virtual function

   Int_t fnb = from->GetEntries();
   Int_t tnb = to->GetEntries();
   if (!fnb || !tnb) {
      return 0;
   }

   UInt_t numBasket = 0;
   Int_t fi = 0;
   Int_t ti = 0;
   while (ti < tnb) {
      TBranch* fb = (TBranch*) from->UncheckedAt(fi);
      TBranch* tb = (TBranch*) to->UncheckedAt(ti);
      Int_t firstfi = fi;
      while (strcmp(fb->GetName(), tb->GetName())) {
         ++fi;
         if (fi >= fnb) {
            // continue at the beginning
            fi = 0;
         }
         if (fi==firstfi) {
            // We tried all the branches and there is not match.
            fb = 0;
            break;
         }
         fb = (TBranch*) from->UncheckedAt(fi);
      }
      if (fb) {
         numBasket += CollectBranches(fb, tb);
         ++fi;
         if (fi >= fnb) {
            fi = 0;
         }
      } else {
         if (tb->GetMother()==tb) {
            // Top level branch.
            if (!(fOptions & kIgnoreMissingTopLevel)) {
               fWarningMsg.Form("One of the export top level branches (%s) is not present in the import TTree.",
                                tb->GetName());
               if (!(fOptions & kNoWarnings)) {
                  Error("TTreeCloner::CollectBranches", "%s", fWarningMsg.Data());
               }
               fIsValid = kFALSE;
            }
         } else {
            fWarningMsg.Form("One of the export sub-branches (%s) is not present in the import TTree.",
                             tb->GetName());
            if (!(fOptions & kNoWarnings)) {
               Error("TTreeCloner::CollectBranches", "%s", fWarningMsg.Data());
            }
            fIsValid = kFALSE;
         }
      }
      ++ti;
   }
   return numBasket;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill the array of branches, matching the branches of the 'from' and 'to' TTrees
/// Returns the total number of baskets in all the branches.

UInt_t TTreeCloner::CollectBranches()
{
   // Since this is called from the constructor, this can not be a virtual function

   if (!fFromTree || !fToTree) {
      return 0;
   }
   UInt_t numBasket = CollectBranches(fFromTree->GetListOfBranches(),
                                      fToTree->GetListOfBranches());

   if (fFromTree->GetBranchRef()) {
      fToTree->BranchRef();
      numBasket += CollectBranches(fFromTree->GetBranchRef(),fToTree->GetBranchRef());
   }
   return numBasket;
}

////////////////////////////////////////////////////////////////////////////////
/// Collect the information about the on-file basket that need
/// to be copied.

void TTreeCloner::CollectBaskets()
{
   UInt_t len = fFromBranches.GetEntries();

   for(UInt_t i=0,bi=0; i<len; ++i) {
      TBranch *from = (TBranch*)fFromBranches.UncheckedAt(i);
      for(Int_t b=0; b<from->GetWriteBasket(); ++b,++bi) {
         fBasketBranchNum[bi] = i;
         fBasketNum[bi] = b;
         fBasketSeek[bi] = from->GetBasketSeek(b);
         //fprintf(stderr,"For %s %d %lld\n",from->GetName(),bi,fBasketSeek[bi]);
         fBasketEntry[bi] = from->GetBasketEntry()[b];
         fBasketIndex[bi] = bi;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Make sure that all the needed TStreamerInfo are
/// present in the output file

void TTreeCloner::CopyStreamerInfos()
{
   TFile *fromFile = fFromTree->GetDirectory()->GetFile();
   TFile *toFile = fToTree->GetDirectory()->GetFile();
   TList *l = fromFile->GetStreamerInfoList();
   TIter next(l);
   TStreamerInfo *oldInfo;
   while ( (oldInfo = (TStreamerInfo*)next()) ) {
      if (oldInfo->IsA() != TStreamerInfo::Class()) {
         continue;
      }
      TStreamerInfo *curInfo = 0;
      TClass *cl = TClass::GetClass(oldInfo->GetName());

      if (!cl->IsLoaded() || cl->GetNew()) {
         // Insure that the TStreamerInfo is loaded
         curInfo = (TStreamerInfo*)cl->GetStreamerInfo(oldInfo->GetClassVersion());
         if (oldInfo->GetClassVersion()==1) {
            // We may have a Foreign class let's look using the
            // checksum:
            TStreamerInfo *matchInfo = (TStreamerInfo*)cl->FindStreamerInfo(oldInfo->GetCheckSum());
            if (matchInfo) {
               curInfo = matchInfo;
            }
         }
         curInfo->ForceWriteInfo(toFile);
      } else {
         // If there is no default constructor the GetStreamerInfo
         // will not work. It also means (hopefully) that an
         // inheriting class has a streamerInfo in the list (which
         // will induces the setting of this streamerInfo)

         oldInfo->ForceWriteInfo(toFile);
      }
   }
   delete l;
}

////////////////////////////////////////////////////////////////////////////////
/// Transfer the basket from the input file to the output file

void TTreeCloner::CopyMemoryBaskets()
{
   TBasket *basket = 0;
   for(Int_t i=0; i<fToBranches.GetEntries(); ++i) {
      TBranch *from = (TBranch*)fFromBranches.UncheckedAt( i );
      TBranch *to   = (TBranch*)fToBranches.UncheckedAt( i );

      basket = from->GetListOfBaskets()->GetEntries() ? from->GetBasket(from->GetWriteBasket()) : 0;
      if (basket) {
         basket = (TBasket*)basket->Clone();
         basket->SetBranch(to);
         to->AddBasket(*basket, kFALSE, fToStartEntries+from->GetBasketEntry()[from->GetWriteBasket()]);
      } else {
         to->AddLastBasket(  fToStartEntries+from->GetBasketEntry()[from->GetWriteBasket()] );
      }
      // In older files, if the branch is a TBranchElement non-terminal 'object' branch, it's basket will contain 0
      // events, in newer file in the same case, the write basket will be missing.
      if (from->GetEntries()!=0 && from->GetWriteBasket()==0 && (basket==0 || basket->GetNevBuf()==0)) {
         to->SetEntries(to->GetEntries()+from->GetEntries());
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Make sure that all the needed TStreamerInfo are
/// present in the output file

void TTreeCloner::CopyProcessIds()
{
   // NOTE: We actually need to merge the ProcessId somehow :(

   TFile *fromfile = fFromTree->GetDirectory()->GetFile();
   TFile *tofile = fToTree->GetDirectory()->GetFile();

   fPidOffset = tofile->GetNProcessIDs();

   TIter next(fromfile->GetListOfKeys());
   TKey *key;
   TDirectory::TContext cur(fromfile);
   while ((key = (TKey*)next())) {
      if (!strcmp(key->GetClassName(),"TProcessID")) {
         TProcessID *pid = (TProcessID*)key->ReadObjectAny(0);
         if (!pid) continue;

         //UShort_t out = TProcessID::WriteProcessID(id,tofile);
         UShort_t out = 0;
         TObjArray *pids = tofile->GetListOfProcessIDs();
         Int_t npids = tofile->GetNProcessIDs();
         Bool_t wasIn = kFALSE;
         for (Int_t i=0;i<npids;i++) {
            if (pids->At(i) == pid) {out = (UShort_t)i; wasIn = kTRUE; break;}
         }

         if (!wasIn) {
            TDirectory *dirsav = gDirectory;
            tofile->cd();
            tofile->SetBit(TFile::kHasReferences);
            pids->AddAtAndExpand(pid,npids);
            pid->IncrementCount();
            char name[32];
            snprintf(name,32,"ProcessID%d",npids);
            pid->Write(name);
            tofile->IncrementProcessIDs();
            if (gDebug > 0) {
               Info("WriteProcessID", "name=%s, file=%s", name, tofile->GetName());
            }
            if (dirsav) dirsav->cd();
            out = (UShort_t)npids;
         }
         if (out<fPidOffset) {
            Error("CopyProcessIDs","Copied %s from %s might already exist!\n",
                  pid->GetName(),fromfile->GetName());
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Create a TFileCacheRead if it was requested.

void TTreeCloner::CreateCache()
{
   if (fCacheSize && fFromTree->GetCurrentFile()) {
      TFile *f = fFromTree->GetCurrentFile();
      auto prev = fFromTree->GetReadCache(f);
      if (fFileCache && prev == fFileCache) {
         return;
      }
      fPrevCache = prev;
      // Remove the previous cache if any.
      if (prev) f->SetCacheRead(nullptr, fFromTree);
      // The constructor attach the new cache.
      fFileCache = new TFileCacheRead(f, fCacheSize, fFromTree);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Restore the TFileCacheRead to its previous value.

void TTreeCloner::RestoreCache() {
   if (IsValid() && fFileCache && fFromTree->GetCurrentFile()) {
      TFile *f = fFromTree->GetCurrentFile();
      f->SetCacheRead(nullptr,fFromTree); // Remove our file cache.
      f->SetCacheRead(fPrevCache, fFromTree);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the entries and import the cluster range of the

void TTreeCloner::ImportClusterRanges()
{
   // First undo, the external call to SetEntries
   // We could improve the interface to optional tell the TTreeCloner that the
   // SetEntries was not done.
   fToTree->SetEntries(fToTree->GetEntries() - fFromTree->GetTree()->GetEntries());

   fToTree->ImportClusterRanges( fFromTree->GetTree() );

   // This is only updated by TTree::Fill upon seeing a Flush event in TTree::Fill
   // So we need to propagate (this has also the advantage of turning on the
   // history recording feature of SetAutoFlush for the next iteration)
   fToTree->fFlushedBytes += fFromTree->fFlushedBytes;

   fToTree->SetEntries(fToTree->GetEntries() + fFromTree->GetTree()->GetEntries());
}

////////////////////////////////////////////////////////////////////////////////
/// Set the TFile cache size to be used.
/// Note that the default is to use the same size as the default TTreeCache for
/// the input tree.
/// \param size Size of the cache. Zero disable the use of the cache.

void TTreeCloner::SetCacheSize(Int_t size)
{
   fCacheSize = size;
   if (IsValid() && fFileCache) {
      if (fCacheSize == 0 || fCacheSize != fFileCache->GetBufferSize()) {
         TFile *f = fFromTree->GetCurrentFile();
         f->SetCacheRead(nullptr,fFromTree);
         delete fFileCache;
         fFileCache = nullptr;
      }
   }
   // Note if the TFile cache is needed, it will be created at the
   // beginning of Exec.
}

////////////////////////////////////////////////////////////////////////////////
/// Sort the basket according to the user request.

void TTreeCloner::SortBaskets()
{
   // Currently this sort __has to__ preserve the order
   // of basket for each individual branch.

   switch (fCloneMethod) {
      case kSortBasketsByBranch:
         // nothing to do, it is already sorted.
         break;
      case kSortBasketsByEntry: {
         for(UInt_t i = 0; i < fMaxBaskets; ++i) { fBasketIndex[i] = i; }
         std::sort(fBasketIndex, fBasketIndex+fMaxBaskets, CompareEntry( this) );
         break;
      }
      case kSortBasketsByOffset:
      default: {
         for(UInt_t i = 0; i < fMaxBaskets; ++i) { fBasketIndex[i] = i; }
         std::sort(fBasketIndex, fBasketIndex+fMaxBaskets, CompareSeek( this) );
         break;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Fill the file cache with the next set of basket.
///
/// \param from index of the first lement of fFromBranches to start caching
/// \return The index of first element of fFromBranches that is not in the cache
UInt_t TTreeCloner::FillCache(UInt_t from)
{
   if (!fFileCache) return 0;
   // Reset the cache
   fFileCache->Prefetch(0, 0);
   Long64_t size = 0;
   for (UInt_t j = from; j < fMaxBaskets; ++j) {
      TBranch *frombr = (TBranch *) fFromBranches.UncheckedAt(fBasketBranchNum[fBasketIndex[j]]);


      Int_t index = fBasketNum[ fBasketIndex[j] ];
      Long64_t pos = frombr->GetBasketSeek(index);
      Int_t len = frombr->GetBasketBytes()[index];
      if (pos && len) {
         size += len;
         if (size > fFileCache->GetBufferSize()) {
            return j;
         }
         fFileCache->Prefetch(pos,len);
      }
   }
   return fMaxBaskets;
}

////////////////////////////////////////////////////////////////////////////////
/// Transfer the basket from the input file to the output file

void TTreeCloner::WriteBaskets()
{
   TBasket *basket = new TBasket();
   for(UInt_t j = 0, notCached = 0; j<fMaxBaskets; ++j) {
      TBranch *from = (TBranch*)fFromBranches.UncheckedAt( fBasketBranchNum[ fBasketIndex[j] ] );
      TBranch *to   = (TBranch*)fToBranches.UncheckedAt( fBasketBranchNum[ fBasketIndex[j] ] );

      TFile *tofile = to->GetFile(0);
      TFile *fromfile = from->GetFile(0);

      Int_t index = fBasketNum[ fBasketIndex[j] ];

      Long64_t pos = from->GetBasketSeek(index);
      if (pos!=0) {
         if (fFileCache && j >= notCached) {
            notCached = FillCache(notCached);
         }
         if (from->GetBasketBytes()[index] == 0) {
            from->GetBasketBytes()[index] = basket->ReadBasketBytes(pos, fromfile);
         }
         Int_t len = from->GetBasketBytes()[index];

         basket->LoadBasketBuffers(pos,len,fromfile,fFromTree);
         basket->IncrementPidOffset(fPidOffset);
         basket->CopyTo(tofile);
         to->AddBasket(*basket,kTRUE,fToStartEntries + from->GetBasketEntry()[index]);
      } else {
         TBasket *frombasket = from->GetBasket( index );
         if (frombasket && frombasket->GetNevBuf()>0) {
            TBasket *tobasket = (TBasket*)frombasket->Clone();
            tobasket->SetBranch(to);
            to->AddBasket(*tobasket, kFALSE, fToStartEntries+from->GetBasketEntry()[index]);
            to->FlushOneBasket(to->GetWriteBasket());
         }
      }
   }
   delete basket;
}
