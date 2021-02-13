// @(#)root/tree:$Id$
// Author: Rene Brun   19/08/2004

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TBranchRef
\ingroup tree

A branch containing and managing a TRefTable for TRef autoloading.
It loads the TBranch containing a referenced object when requested
by TRef::GetObject(), so the reference can be resolved. The
information which branch to load is stored by TRefTable. Once a
TBranch has read the TBranchRef's current entry it will not be told
to re-read, in case the use has changed objects read from the
branch.

### LIMITATION :
Note that this does NOT allow for autoloading of references spanning
different entries. The TBranchRef's current entry has to correspond
to the entry of the TBranch containing the referenced object.

The TRef cannot be stored in a top-level branch which is a
TBranchObject for the auto-loading to work. E.g. you cannot store
the TRefs in TObjArray, and create a top-level branch storing this
TObjArray.
*/

#include "TBranchRef.h"
#include "TTree.h"
#include "TBasket.h"
#include "TRefTable.h"
#include "TFile.h"
#include "TFriendElement.h"

ClassImp(TBranchRef);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor.

TBranchRef::TBranchRef(): TBranch(), fRequestedEntry(-1), fRefTable(0)
{
   fReadLeaves = (ReadLeaves_t)&TBranchRef::ReadLeavesImpl;
   fFillLeaves = (FillLeaves_t)&TBranchRef::FillLeavesImpl;
}

////////////////////////////////////////////////////////////////////////////////
/// Main constructor called by TTree::BranchRef.

TBranchRef::TBranchRef(TTree *tree)
    : TBranch(), fRequestedEntry(-1), fRefTable(0)
{
   if (!tree) return;
   SetName("TRefTable");
   SetTitle("List of branch numbers with referenced objects");
   fRefTable = new TRefTable(this,100);

   fCompress       = 1;
   fBasketSize     = 32000;
   fAddress        = 0;
   fBasketBytes    = new Int_t[fMaxBaskets];
   fBasketEntry    = new Long64_t[fMaxBaskets];
   fBasketSeek     = new Long64_t[fMaxBaskets];

   for (Int_t i=0;i<fMaxBaskets;i++) {
      fBasketBytes[i] = 0;
      fBasketEntry[i] = 0;
      fBasketSeek[i]  = 0;
   }

   fTree       = tree;
   fMother     = this;
   fDirectory  = fTree->GetDirectory();
   fFileName   = "";
   fReadLeaves = (ReadLeaves_t)&TBranchRef::ReadLeavesImpl;
   fFillLeaves = (FillLeaves_t)&TBranchRef::FillLeavesImpl;
}

////////////////////////////////////////////////////////////////////////////////
/// Typical destructor.

TBranchRef::~TBranchRef()
{
   delete fRefTable;
}

////////////////////////////////////////////////////////////////////////////////
/// Clear entries in the TRefTable.

void TBranchRef::Clear(Option_t *option)
{
   if (fRefTable) fRefTable->Clear(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Fill the branch basket with the referenced objects parent numbers.

Int_t TBranchRef::FillImpl(ROOT::Internal::TBranchIMTHelper *imtHelper)
{
   Int_t nbytes = TBranch::FillImpl(imtHelper);
   return nbytes;
}

////////////////////////////////////////////////////////////////////////////////
/// This function is called by TRefTable::Notify, itself called by
/// TRef::GetObject.
/// The function reads the branch containing the object referenced
/// by the TRef.

Bool_t TBranchRef::Notify()
{
   if (!fRefTable) fRefTable = new TRefTable(this,100);
   UInt_t uid = fRefTable->GetUID();
   TProcessID* context = fRefTable->GetUIDContext();
   if (fReadEntry != fRequestedEntry) {
      // Load the RefTable if we need to.
      GetEntry(fRequestedEntry);
   }
   TBranch *branch = (TBranch*)fRefTable->GetParent(uid, context);
   if (branch) {
      // don't re-read, the user might have changed some object
      if (branch->GetReadEntry() != fRequestedEntry)
         branch->GetEntry(fRequestedEntry);
   } else {
      //scan the TRefTable of possible friend Trees
      TList *friends = fTree->GetListOfFriends();
      if (!friends) return kTRUE;
      TObjLink *lnk = friends->FirstLink();
      while (lnk) {
         TFriendElement* elem = (TFriendElement*)lnk->GetObject();
         TTree *tree = elem->GetTree();
         TBranchRef *bref = tree->GetBranchRef();
         if (bref) {
            if (bref->GetReadEntry() != fRequestedEntry) {
               bref->GetEntry(fRequestedEntry);
            }
            branch = (TBranch*)bref->GetRefTable()->GetParent(uid, context);
            if (branch) {
               // don't re-read, the user might have changed some object
               if (branch->GetReadEntry() != fRequestedEntry)
                  branch->GetEntry(fRequestedEntry);
               return kTRUE;
            }
         }
         lnk = lnk->Next();
      }
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Print the TRefTable branch.

void TBranchRef::Print(Option_t *option) const
{
   TBranch::Print(option);
}

////////////////////////////////////////////////////////////////////////////////
/// This function called by TBranch::GetEntry overloads TBranch::ReadLeaves.

void TBranchRef::ReadLeavesImpl(TBuffer &b)
{
   if (!fRefTable) fRefTable = new TRefTable(this,100);
   fRefTable->ReadBuffer(b);
}

////////////////////////////////////////////////////////////////////////////////
/// This function called by TBranch::Fill overloads TBranch::FillLeaves.

void TBranchRef::FillLeavesImpl(TBuffer &b)
{
   if (!fRefTable) fRefTable = new TRefTable(this,100);
   fRefTable->FillBuffer(b);
}

////////////////////////////////////////////////////////////////////////////////
/// - Existing buffers are deleted
/// - Entries, max and min are reset
/// - TRefTable is cleared.

void TBranchRef::Reset(Option_t *option)
{
   TBranch::Reset(option);
   if (!fRefTable) fRefTable = new TRefTable(this,100);
   fRefTable->Reset();
}

////////////////////////////////////////////////////////////////////////////////
/// Reset a Branch after a Merge operation (drop data but keep customizations)
/// TRefTable is cleared.

void TBranchRef::ResetAfterMerge(TFileMergeInfo *info)
{
   TBranch::ResetAfterMerge(info);
   if (!fRefTable) fRefTable = new TRefTable(this,100);
   fRefTable->Reset();
}

////////////////////////////////////////////////////////////////////////////////
/// Set the current parent branch.
///
/// This function is called by TBranchElement::GetEntry()
/// and TBranchElement::Fill() when reading or writing
/// branches that may contain referenced objects.

Int_t TBranchRef::SetParent(const TObject* object, Int_t branchID)
{
   if (!fRefTable) {
      fRefTable = new TRefTable(this, 100);
   }
   TRefTable::SetRefTable(fRefTable);
   return fRefTable->SetParent(object, branchID);
}

