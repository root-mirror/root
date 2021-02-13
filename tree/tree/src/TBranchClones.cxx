// @(#)root/tree:$Id$
// Author: Rene Brun   11/02/96

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TBranchClones.h"

#include "TBasket.h"
#include "TBuffer.h"
#include "TClass.h"
#include "TClonesArray.h"
#include "TDataMember.h"
#include "TDataType.h"
#include "TFile.h"
#include "TLeafI.h"
#include "TRealData.h"
#include "TTree.h"

#include <cstring>

ClassImp(TBranchClones);

/** \class TBranchClones
\ingroup tree

A Branch for the case of an array of clone objects.

See TTree.
*/

////////////////////////////////////////////////////////////////////////////////
/// Default and i/o constructor.

TBranchClones::TBranchClones()
: TBranch()
, fList(0)
, fRead(0)
, fN(0)
, fNdataMax(0)
, fBranchCount(0)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TBranchClones::TBranchClones(TTree *tree, const char* name, void* pointer, Int_t basketsize, Int_t compress, Int_t splitlevel)
: TBranch()
, fList(0)
, fRead(0)
, fN(0)
, fNdataMax(0)
, fBranchCount(0)
{
   Init(tree,0,name,pointer,basketsize,compress,splitlevel);
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TBranchClones::TBranchClones(TBranch *parent, const char* name, void* pointer, Int_t basketsize, Int_t compress, Int_t splitlevel)
: TBranch()
, fList(0)
, fRead(0)
, fN(0)
, fNdataMax(0)
, fBranchCount(0)
{
   Init(0,parent,name,pointer,basketsize,compress,splitlevel);
}

////////////////////////////////////////////////////////////////////////////////
/// Initialization (non-virtual, to be called from constructor).

void TBranchClones::Init(TTree *tree, TBranch *parent, const char* name, void* pointer, Int_t basketsize, Int_t compress, Int_t splitlevel)
{
   if (tree==0 && parent!=0) tree = parent->GetTree();
   fTree   = tree;
   fMother = parent ? parent->GetMother() : this;
   fParent = parent;

   TString leaflist;
   TString branchname;
   TString branchcount;
   SetName(name);
   if ((compress == -1) && tree->GetDirectory()) {
      TFile* bfile = 0;
      if (tree->GetDirectory()) {
         bfile = tree->GetDirectory()->GetFile();
      }
      if (bfile) {
         compress = bfile->GetCompressionSettings();
      }
   }
   char* cpointer = (char*) pointer;
   char** ppointer = (char**) pointer;
   fList = (TClonesArray*) *ppointer;
   fAddress = cpointer;
   TClass* cl = fList->GetClass();
   if (!cl) {
      return;
   }
   tree->BuildStreamerInfo(cl);
   fClassName = cl->GetName();
   fSplitLevel = splitlevel;

   // Create a branch to store the array count.
   if (basketsize < 100) {
      basketsize = 100;
   }
   leaflist.Form("%s_/I", name);
   branchcount.Form("%s_", name);
   fBranchCount = new TBranch(this, branchcount, &fN, leaflist, basketsize);
   fBranchCount->SetBit(kIsClone);
   TLeaf* leafcount = (TLeaf*) fBranchCount->GetListOfLeaves()->UncheckedAt(0);
   fDirectory = fTree->GetDirectory();
   fFileName = "";

   // Loop on all public data members of the class and its base classes.
   const char* itype = 0;
   TRealData* rd = 0;
   TIter next(cl->GetListOfRealData());
   while ((rd = (TRealData *) next())) {
      if (rd->TestBit(TRealData::kTransient)) continue;

      if (rd->IsObject()) {
         continue;
      }
      TDataMember* member = rd->GetDataMember();
      if (!member->IsPersistent()) {
         // -- Skip non-persistent members.
         continue;
      }
      if (!member->IsBasic() || member->IsaPointer()) {
         Warning("BranchClones", "Cannot process: %s::%s", cl->GetName(), member->GetName());
         continue;
      }
      // Forget TObject part if splitlevel = 2.
      if ((splitlevel > 1) || fList->TestBit(TClonesArray::kForgetBits) || cl->CanIgnoreTObjectStreamer()) {
         if (!std::strcmp(member->GetName(), "fBits")) {
            continue;
         }
         if (!std::strcmp(member->GetName(), "fUniqueID")) {
            continue;
         }
      }
      tree->BuildStreamerInfo(TClass::GetClass(member->GetFullTypeName()));
      TDataType* membertype = member->GetDataType();
      Int_t type = membertype->GetType();
      if (!type) {
         Warning("BranchClones", "Cannot process: %s::%s of type zero!", cl->GetName(), member->GetName());
         continue;
      }

      if (type ==  1) {
         itype = "B";
      } else if (type == 2) {
         itype = "S";
      } else if (type == 3) {
         itype = "I";
      } else if (type == 5) {
         itype = "F";
      } else if (type == 8) {
         itype = "D";
      } else if (type == 9) {
         itype = "D";
      } else if (type == 11) {
         itype = "b";
      } else if (type == 12) {
         itype = "s";
      } else if (type == 13) {
         itype = "i";
      }

      leaflist.Form("%s[%s]/%s", member->GetName(), branchcount.Data(), itype);
      Int_t comp = compress;
      branchname.Form("%s.%s", name, rd->GetName());
      TBranch* branch  = new TBranch(this, branchname, this, leaflist, basketsize, comp);
      branch->SetBit(kIsClone);
      TObjArray* leaves = branch->GetListOfLeaves();
      TLeaf* leaf = (TLeaf*) leaves->UncheckedAt(0);
      leaf->SetOffset(rd->GetThisOffset());
      leaf->SetLeafCount(leafcount);
      Int_t arraydim = member->GetArrayDim();
      if (arraydim) {
         Int_t maxindex = 1;
         while (arraydim) {
            maxindex *= member->GetMaxIndex(--arraydim);
         }
         leaf->SetLen(maxindex);
      }
      fBranches.Add(branch);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TBranchClones::~TBranchClones()
{
   delete fBranchCount;
   fBranchCount = 0;
   fBranches.Delete();
   // FIXME: We might own this, possible memory leak.
   fList = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Browse this branch.

void TBranchClones::Browse(TBrowser* b)
{
   fBranches.Browse(b);
}

////////////////////////////////////////////////////////////////////////////////
/// Loop on all branches and fill Basket buffer.

Int_t TBranchClones::FillImpl(ROOT::Internal::TBranchIMTHelper *imtHelper)
{
   Int_t i = 0;
   Int_t nbytes = 0;
   Int_t nbranches = fBranches.GetEntriesFast();
   char** ppointer = (char**) fAddress;
   if (!ppointer) {
      return 0;
   }
   fList = (TClonesArray*) *ppointer;
   fN = fList->GetEntriesFast();
   fEntries++;
   if (fN > fNdataMax) {
      fNdataMax = fList->GetSize();
      TString branchcount;
      branchcount.Form("%s_", GetName());
      TLeafI* leafi = (TLeafI*) fBranchCount->GetLeaf(branchcount);
      leafi->SetMaximum(fNdataMax);
      for (i = 0; i < nbranches; i++) {
         TBranch* branch = (TBranch*) fBranches.UncheckedAt(i);
         TObjArray* leaves = branch->GetListOfLeaves();
         TLeaf* leaf = (TLeaf*) leaves->UncheckedAt(0);
         leaf->SetAddress();
      }
   }
   nbytes += fBranchCount->FillImpl(imtHelper);
   for (i = 0; i < nbranches; i++)  {
      TBranch* branch = (TBranch*) fBranches.UncheckedAt(i);
      TObjArray* leaves = branch->GetListOfLeaves();
      TLeaf* leaf = (TLeaf*) leaves->UncheckedAt(0);
      leaf->Import(fList, fN);
      nbytes += branch->FillImpl(imtHelper);
   }
   return nbytes;
}

////////////////////////////////////////////////////////////////////////////////
/// Read all branches and return total number of bytes read.

Int_t TBranchClones::GetEntry(Long64_t entry, Int_t getall)
{
   if (TestBit(kDoNotProcess) && !getall) {
      return 0;
   }
   Int_t nbytes = fBranchCount->GetEntry(entry, getall);
   TLeaf* leafcount = (TLeaf*) fBranchCount->GetListOfLeaves()->UncheckedAt(0);
   fN = Int_t(leafcount->GetValue());
   if (fN <= 0) {
      if (fList) {
         fList->Clear();
      }
      return 0;
   }
   TBranch* branch = 0;
   Int_t nbranches = fBranches.GetEntriesFast();
   // If fList exists, create clones array objects.
   if (fList) {
      fList->Clear();
      fList->ExpandCreateFast(fN);
      for (Int_t i = 0; i < nbranches; i++)  {
         branch = (TBranch*) fBranches.UncheckedAt(i);
         if (((TLeaf*) branch->GetListOfLeaves()->UncheckedAt(0))->GetOffset() < 0) {
            continue;
         }
         nbytes += branch->GetEntryExport(entry, getall, fList, fN);
      }
   } else {
      for (Int_t i = 0; i < nbranches; i++)  {
         branch = (TBranch*) fBranches.UncheckedAt(i);
         nbytes += branch->GetEntry(entry, getall);
      }
   }
   return nbytes;
}

////////////////////////////////////////////////////////////////////////////////
/// Print branch parameters.

void TBranchClones::Print(Option_t *option) const
{
   fBranchCount->Print(option);
   Int_t nbranches = fBranches.GetEntriesFast();
   for (Int_t i = 0; i < nbranches; i++) {
      TBranch* branch = (TBranch*) fBranches.At(i);
      branch->Print(option);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Reset branch.
///
/// - Existing buffers are deleted
/// - Entries, max and min are reset

void TBranchClones::Reset(Option_t* option)
{
   fEntries = 0;
   fTotBytes = 0;
   fZipBytes = 0;
   Int_t nbranches = fBranches.GetEntriesFast();
   for (Int_t i = 0; i < nbranches; i++) {
      TBranch* branch = (TBranch*) fBranches.At(i);
      branch->Reset(option);
   }
   fBranchCount->Reset();
}

////////////////////////////////////////////////////////////////////////////////
/// Reset branch after a merge.
///
/// - Existing buffers are deleted
/// - Entries, max and min are reset

void TBranchClones::ResetAfterMerge(TFileMergeInfo *info)
{
   fEntries  = 0;
   fTotBytes = 0;
   fZipBytes = 0;
   Int_t nbranches = fBranches.GetEntriesFast();
   for (Int_t i = 0; i < nbranches; i++) {
      TBranch* branch = (TBranch*) fBranches.At(i);
      branch->ResetAfterMerge(info);
   }
   fBranchCount->ResetAfterMerge(info);
}

////////////////////////////////////////////////////////////////////////////////
/// Set address of this branch.

void TBranchClones::SetAddress(void* addr)
{
   fReadEntry = -1;
   fAddress = (char*) addr;
   char** pp= (char**) fAddress;
   if (pp && (*pp == 0)) {
      // We've been asked to allocate an object for the user.
      *pp= (char*) new TClonesArray(fClassName);
   }
   fList = 0;
   if (pp) {
      fList = (TClonesArray*) *pp;
   }
   fBranchCount->SetAddress(&fN);
}

////////////////////////////////////////////////////////////////////////////////
/// Reset basket size for all sub-branches.

void TBranchClones::SetBasketSize(Int_t buffsize)
{
   TBranch::SetBasketSize(buffsize);

   Int_t nbranches = fBranches.GetEntriesFast();
   for (Int_t i = 0; i < nbranches; i++)  {
      TBranch* branch = (TBranch*) fBranches[i];
      branch->SetBasketSize(fBasketSize);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Serialize/Deserialize from a buffer.

void TBranchClones::Streamer(TBuffer& b)
{
   UInt_t R__s, R__c;
   if (b.IsReading()) {
      b.ReadVersion(&R__s, &R__c);
      TNamed::Streamer(b);
      b >> fCompress;
      b >> fBasketSize;
      b >> fEntryOffsetLen;
      b >> fMaxBaskets;
      b >> fWriteBasket;
      b >> fEntryNumber;
      b >> fEntries;
      b >> fTotBytes;
      b >> fZipBytes;
      b >> fOffset;
      b >> fBranchCount;
      fClassName.Streamer(b);
      fBranches.Streamer(b);
      fTree = 0;
      TBranch* branch = 0;
      TLeaf* leaf = 0;
      Int_t nbranches = fBranches.GetEntriesFast();
      for (Int_t i = 0; i < nbranches; i++) {
         branch = (TBranch*) fBranches[i];
         branch->SetBit(kIsClone);
         leaf = (TLeaf*) branch->GetListOfLeaves()->UncheckedAt(0);
         leaf->SetOffset(-1);
      }
      fRead = 1;
      TClass* cl = TClass::GetClass((const char*) fClassName);
      if (!cl) {
         Warning("Streamer", "Unknown class: %s. Cannot read BranchClones: %s", fClassName.Data(), GetName());
         SetBit(kDoNotProcess);
         return;
      }
      if (!cl->GetListOfRealData()) {
         cl->BuildRealData();
      }
      TString branchname;
      TRealData* rd = 0;
      TIter next(cl->GetListOfRealData());
      while ((rd = (TRealData*) next())) {
         if (rd->TestBit(TRealData::kTransient)) continue;

         TDataMember* member = rd->GetDataMember();
         if (!member || !member->IsBasic() || !member->IsPersistent()) {
            continue;
         }
         TDataType* membertype = member->GetDataType();
         if (!membertype->GetType()) {
            continue;
         }
         branchname.Form("%s.%s", GetName(), rd->GetName());
         branch = (TBranch*) fBranches.FindObject(branchname);
         if (!branch) {
            continue;
         }
         TObjArray* leaves = branch->GetListOfLeaves();
         leaf = (TLeaf*) leaves->UncheckedAt(0);
         leaf->SetOffset(rd->GetThisOffset());
      }
      b.CheckByteCount(R__s, R__c, TBranchClones::IsA());
   } else {
      R__c = b.WriteVersion(TBranchClones::IsA(), kTRUE);
      TNamed::Streamer(b);
      b << fCompress;
      b << fBasketSize;
      b << fEntryOffsetLen;
      b << fMaxBaskets;
      b << fWriteBasket;
      b << fEntryNumber;
      b << fEntries;
      b << fTotBytes;
      b << fZipBytes;
      b << fOffset;
      b << fBranchCount;
      fClassName.Streamer(b);
      fBranches.Streamer(b);
      b.SetByteCount(R__c, kTRUE);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Refresh the value of fDirectory (i.e. where this branch writes/reads its buffers)
/// with the current value of fTree->GetCurrentFile unless this branch has been
/// redirected to a different file.  Also update the sub-branches.

void TBranchClones::UpdateFile()
{
   fBranchCount->UpdateFile();
   TBranch::UpdateFile();
}
