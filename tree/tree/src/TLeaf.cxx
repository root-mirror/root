// @(#)root/tree:$Id$
// Author: Rene Brun   12/01/96

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TLeaf
\ingroup tree

A TLeaf describes individual elements of a TBranch
See TBranch structure in TTree.
*/

#include "TLeaf.h"
#include "TBranch.h"
#include "TBuffer.h"
#include "TTree.h"
#include "TVirtualPad.h"
#include "TBrowser.h"

#include <cctype>

ClassImp(TLeaf);

////////////////////////////////////////////////////////////////////////////////

TLeaf::TLeaf()
   : TNamed()
   , fNdata(0)
   , fLen(0)
   , fLenType(0)
   , fOffset(0)
   , fIsRange(kFALSE)
   , fIsUnsigned(kFALSE)
   , fLeafCount(0)
   , fBranch(0)
   , fLeafCountValues(0)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Create a Leaf.
///
/// See the TTree and TBranch constructors for explanation of parameters.

TLeaf::TLeaf(TBranch *parent, const char* name, const char *)
   : TNamed(name, name)
   , fNdata(0)
   , fLen(0)
   , fLenType(4)
   , fOffset(0)
   , fIsRange(kFALSE)
   , fIsUnsigned(kFALSE)
   , fLeafCount(0)
   , fBranch(parent)
   , fLeafCountValues(0)
{
   fLeafCount = GetLeafCounter(fLen);

   if (fLen == -1) {
      MakeZombie();
      return;
   }

   const char *bracket = strchr(name, '[');
   if (bracket) fName.ReplaceAll(bracket,"");
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.

TLeaf::TLeaf(const TLeaf& lf) :
  TNamed(lf),
  fNdata(lf.fNdata),
  fLen(lf.fLen),
  fLenType(lf.fLenType),
  fOffset(lf.fOffset),
  fIsRange(lf.fIsRange),
  fIsUnsigned(lf.fIsUnsigned),
  fLeafCount(lf.fLeafCount),
  fBranch(lf.fBranch),
  fLeafCountValues(nullptr)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Assignment operator.

TLeaf& TLeaf::operator=(const TLeaf& lf)
{
   if(this!=&lf) {
      TNamed::operator=(lf);
      fNdata=lf.fNdata;
      fLen=lf.fLen;
      fLenType=lf.fLenType;
      fOffset=lf.fOffset;
      fIsRange=lf.fIsRange;
      fIsUnsigned=lf.fIsUnsigned;
      fLeafCount=lf.fLeafCount;
      fBranch=lf.fBranch;
      if (fLeafCountValues) {
         fLeafCountValues->fStartEntry = -1;
         fLeafCountValues->fValues.resize(0);
      }
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TLeaf::~TLeaf()
{
   if (fBranch) {
      TTree* tree = fBranch->GetTree();
      fBranch = 0;
      if (tree) {
         TObjArray *lst = tree->GetListOfLeaves();
         if (lst->GetLast()!=-1) lst->Remove(this);
      }
   }
   fLeafCount = 0;
   delete fLeafCountValues;
}

////////////////////////////////////////////////////////////////////////////////
/// Browse the content of this leaf.

void TLeaf::Browse(TBrowser* b)
{
   if (strchr(GetName(), '.')) {
      fBranch->GetTree()->Draw(GetName(), "", b ? b->GetDrawOption() : "");
   } else {
      if ((fBranch->GetListOfLeaves()->GetEntries() > 1) ||
          (strcmp(fBranch->GetName(), GetName()) != 0)) {
         TString name(fBranch->GetName());
         if (!name.EndsWith(".")) name += ".";
         name += GetName();
         fBranch->GetTree()->Draw(name, "", b ? b->GetDrawOption() : "");
      } else {
         fBranch->GetTree()->Draw(GetName(), "", b ? b->GetDrawOption() : "");
      }
   }
   if (gPad) {
      gPad->Update();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Pack leaf elements in Basket output buffer.

void TLeaf::FillBasket(TBuffer &)
{
}

////////////////////////////////////////////////////////////////////////////////
/// If the class supports it, generate an offset array base.
///
/// This class only returns `nullptr` on error.
Int_t *TLeaf::GenerateOffsetArrayBase(Int_t base, Int_t events) const
{
   // In order to avoid a virtual call, we assume ROOT developers will override
   // the default GenerateOffsetArray for cases where this function does not apply.

   Int_t *retval = new Int_t[events];
   if (R__unlikely(!retval || !fLeafCount)) {
      delete [] retval;
      return nullptr;
   }

   Long64_t orig_entry = std::max(fBranch->GetReadEntry(), 0LL); // -1 indicates to start at the beginning
   const std::vector<Int_t> *countValues = fLeafCount->GetLeafCountValues(orig_entry, events);

   if (!countValues || ((Int_t)countValues->size()) < events) {
      Error("GenerateOffsetArrayBase", "The leaf %s could not retrieve enough entries from its branch count (%s), ask for %d and got %ld",
            GetName(), fLeafCount->GetName(), events, (long)(countValues ? countValues->size() : -1));
      delete [] retval;
      return nullptr;
   }

   Int_t header = GetOffsetHeaderSize();
   Int_t len = 0;
   for (Int_t idx = 0, offset = base; idx < events; idx++) {
      retval[idx] = offset;
      len = (*countValues)[idx];
      offset += fLenType * len + header;
   }

   return retval;
}


////////////////////////////////////////////////////////////////////////////////
/// Return the full name (including the parent's branch names) of the leaf.

TString TLeaf::GetFullName() const
{
   TString branchname = GetBranch()->GetFullName();
   if (branchname.Length() && (branchname[branchname.Length()-1] == '.'))
      return branchname + GetName();
   else
      return branchname + "." + GetName();
}

////////////////////////////////////////////////////////////////////////////////
/// Return a pointer to the counter of this leaf (if any) or store the number of elements that the leaf contains in
/// countval.
///
/// - If leaf name has the form var[nelem], where nelem is alphanumeric, then
///     if nelem is a leaf name, return countval = 1 and the pointer to
///          the leaf named nelem, otherwise return 0.
/// - If leaf name has the form var[nelem], where nelem is a non-negative integer, then
///     return countval = nelem and a null pointer.
/// - If leaf name has the form of a multi-dimensional array (e.g. var[nelem][nelem2]
///     where nelem and nelem2 are non-negative integers) then
///     return countval = product of all dimension sizes and a null pointer.
/// - If leaf name has the form var[... (and does not match the previous 2
///     cases) return countval = -1 and null pointer;
/// - Otherwise return countval = 1 and a null pointer.

TLeaf* TLeaf::GetLeafCounter(Int_t& countval) const
{
   countval = 1;
   const char* name = GetTitle();
   char* bleft = (char*) strchr(name, '[');
   if (!bleft) {
      return 0;
   }
   bleft++;
   Int_t nch = strlen(bleft);
   char* countname = new char[nch+1];
   strcpy(countname, bleft);
   char* bright = (char*) strchr(countname, ']');
   if (!bright) {
      delete[] countname;
      countname = 0;
      countval = -1;
      return 0;
   }
   char *bleft2 = (char*) strchr(countname, '[');
   *bright = 0;
   nch = strlen(countname);

   // Now search a branch name with a leaf name = countname
   if (fBranch == 0) {
      Error("GetLeafCounter","TLeaf %s is not setup properly, fBranch is null.",GetName());
      delete[] countname;
      return 0;
   }
   if (fBranch->GetTree() == 0) {
      Error("GetLeafCounter","For Leaf %s, the TBranch %s is not setup properly, fTree is null.",GetName(),fBranch->GetName());
      delete[] countname;
      return 0;
   }
   TTree* pTree = fBranch->GetTree();

   TLeaf* leaf = (TLeaf*) GetBranch()->GetListOfLeaves()->FindObject(countname);
   if (leaf == 0) {
      // Try outside the branch:
      leaf = (TLeaf*) pTree->GetListOfLeaves()->FindObject(countname);
   }
   //if not found, make one more trial in case the leaf name has a "."
   if (!leaf && strchr(GetName(), '.')) {
      char* withdot = new char[strlen(GetName())+strlen(countname)+1];
      strcpy(withdot, GetName());
      char* lastdot = strrchr(withdot, '.');
      strcpy(lastdot, countname);
      leaf = (TLeaf*) pTree->GetListOfLeaves()->FindObject(countname);
      delete[] withdot;
      withdot = 0;
   }
   if (!leaf && strchr(countname,'.')) {
      // Not yet found and the countname has a dot in it, let's try
      // to find the leaf using its full name
      leaf = pTree->FindLeaf(countname);
   }
   Int_t i = 0;
   if (leaf) {
      countval = 1;
      leaf->SetRange();
      if (bleft2) {
         sscanf(bleft2, "[%d]", &i);
         countval *= i;
      }
      bleft = bleft2;
      while (bleft) {
         bleft2++;
         bleft = (char*) strchr(bleft2, '[');
         if (!bleft) {
            break;
         }
         sscanf(bleft, "[%d]", &i);
         countval *= i;
         bleft2 = bleft;
      }
      delete[] countname;
      countname = 0;
      return leaf;
   }
   // not found in a branch/leaf. Is it a numerical value?
   for (i = 0; i < nch; i++) {
      if (!isdigit(countname[i])) {
         delete[] countname;
         countname = 0;
         countval = -1;
         return 0;
      }
   }
   sscanf(countname, "%d", &countval);
   if (bleft2) {
      sscanf(bleft2, "[%d]", &i);
      countval *= i;
   }
   bleft = bleft2;
   while (bleft) {
      bleft2++;
      bleft = (char*) strchr(bleft2, '[');
      if (!bleft) {
         break;
      }
      sscanf(bleft, "[%d]", &i);
      countval *= i;
      bleft2 = bleft;
   }

   delete[] countname;
   countname = 0;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// If this branch is a branch count, return the set of collection size for
/// the entry range requested
/// start: first entry to read and return information about
/// len: number of entries to read.
const TLeaf::Counts_t *TLeaf::GetLeafCountValues(Long64_t start, Long64_t len)
{
   if (len <= 0 || !IsRange())
     return nullptr;

   if (fLeafCountValues) {
      if (fLeafCountValues->fStartEntry == start && len < (Long64_t)fLeafCountValues->fValues.size())
      {
         return &fLeafCountValues->fValues;
      }
      if (start >= fLeafCountValues->fStartEntry &&
          (start+len) <= (Long64_t)(fLeafCountValues->fStartEntry + fLeafCountValues->fValues.size()))
      {
         auto &values(fLeafCountValues->fValues);
         values.erase(values.begin(), values.begin() + start-fLeafCountValues->fStartEntry);
         return &values;
      }
   } else {
      fLeafCountValues = new LeafCountValues();
   }


   fLeafCountValues->fValues.clear();
   fLeafCountValues->fValues.reserve(len);
   fLeafCountValues->fStartEntry = start;

   auto branch = GetBranch();
   Long64_t orig_leaf_entry = branch->GetReadEntry();
   for (Long64_t idx = 0; idx < len; ++idx) {
       branch->GetEntry(start + idx);
       auto size = static_cast<Int_t>(GetValue());
       fLeafCountValues->fValues.push_back( size );
   }
   branch->GetEntry(orig_leaf_entry);
   return &(fLeafCountValues->fValues);
}

////////////////////////////////////////////////////////////////////////////////
/// Return the number of effective elements of this leaf, for the current entry.

Int_t TLeaf::GetLen() const
{
   if (fLeafCount) {
      // -- We are a varying length array.
      Int_t len = Int_t(fLeafCount->GetValue());
      if (len > fLeafCount->GetMaximum()) {
         Error("GetLen", "Leaf counter is greater than maximum!  leaf: '%s' len: %d max: %d", GetName(), len, fLeafCount->GetMaximum());
         len = fLeafCount->GetMaximum();
      }
      return len * fLen;
   } else {
      // -- We are a fixed size thing.
      return fLen;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Helper routine for TLeafX::SetAddress.
///
/// The return value is non-zero if we owned the old
/// value buffer and must delete it now.  The size
/// of the value buffer is recalculated and stored,
/// and a decision is made whether or not we own the
/// new value buffer.

Int_t TLeaf::ResetAddress(void* addr, Bool_t calledFromDestructor)
{
   // The kNewValue bit records whether or not we own
   // the current value buffer or not.  If we own it,
   // then we are responsible for deleting it.
   Bool_t deleteValue = kFALSE;
   if (TestBit(kNewValue)) {
      deleteValue = kTRUE;
   }
   // If we are not being called from a destructor,
   // recalculate the value buffer size and decide
   // whether or not we own the new value buffer.
   if (!calledFromDestructor) {
      // -- Recalculate value buffer size and decide ownership of value.
      if (fLeafCount) {
         // -- Varying length array data member.
         fNdata = (fLeafCount->GetMaximum() + 1) * fLen;
      } else {
         // -- Fixed size data member.
         fNdata = fLen;
      }
      // If we were provided an address, then we do not own
      // the value, otherwise we do and must delete it later,
      // keep track of this with bit kNewValue.
      if (addr) {
         ResetBit(kNewValue);
      } else {
         SetBit(kNewValue);
      }
   }
   return deleteValue;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the leaf count of this leaf.

void TLeaf::SetLeafCount(TLeaf *leaf)
{
   if (IsZombie() && (fLen == -1) && leaf) {
      // The constructor noted that it could not find the
      // leafcount.  Now that we did find it, let's remove
      // the side-effects.
      ResetBit(kZombie);
      fLen = 1;
   }
   fLeafCount = leaf;
}

////////////////////////////////////////////////////////////////////////////////
/// Stream a class object.

void TLeaf::Streamer(TBuffer &b)
{
   if (b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = b.ReadVersion(&R__s, &R__c);
      if (R__v > 1) {
         b.ReadClassBuffer(TLeaf::Class(), this, R__v, R__s, R__c);
      } else {
         // -- Process old versions before automatic schema evolution.
         TNamed::Streamer(b);
         b >> fLen;
         b >> fLenType;
         b >> fOffset;
         b >> fIsRange;
         b >> fIsUnsigned;
         b >> fLeafCount;
         b.CheckByteCount(R__s, R__c, TLeaf::IsA());
      }
      if (!fLen) {
         fLen = 1;
      }
      // We do not own the value buffer right now.
      ResetBit(kNewValue);
      SetAddress();
   } else {
      b.WriteClassBuffer(TLeaf::Class(), this);
   }
}

