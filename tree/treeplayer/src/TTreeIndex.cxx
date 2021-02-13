// @(#)root/tree:$Id$
// Author: Rene Brun   05/07/2004

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TTreeIndex
A Tree Index with majorname and minorname.
*/

#include "TTreeIndex.h"

#include "TTreeFormula.h"
#include "TTree.h"
#include "TBuffer.h"
#include "TMath.h"

ClassImp(TTreeIndex);


struct IndexSortComparator {

  IndexSortComparator(Long64_t *major, Long64_t *minor)
        : fValMajor(major), fValMinor(minor)
  {}

   template<typename Index>
   bool operator()(Index i1, Index i2) {
      if( *(fValMajor + i1) == *(fValMajor + i2) )
         return *(fValMinor + i1) < *(fValMinor + i2);
      else
         return *(fValMajor + i1) < *(fValMajor + i2);
   }

  // pointers to the start of index values tables keeping uppder 64bit and lower 64bit
  // of combined indexed 128bit value
  Long64_t *fValMajor, *fValMinor;
};


////////////////////////////////////////////////////////////////////////////////
/// Default constructor for TTreeIndex

TTreeIndex::TTreeIndex(): TVirtualIndex()
{
   fTree               = 0;
   fN                  = 0;
   fIndexValues        = 0;
   fIndexValuesMinor   = 0;
   fIndex              = 0;
   fMajorFormula       = 0;
   fMinorFormula       = 0;
   fMajorFormulaParent = 0;
   fMinorFormulaParent = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Normal constructor for TTreeIndex
///
/// Build an index table using the leaves of Tree T with  major & minor names
/// The index is built with the expressions given in "majorname" and "minorname".
///
/// a Long64_t array fIndexValues is built with:
///
/// -  major = the value of majorname converted to an integer
/// -  minor = the value of minorname converted to an integer
/// -  fIndexValues[i] = major<<31 + minor
///
/// This array is sorted. The sorted fIndex[i] contains the serial number
/// in the Tree corresponding to the pair "major,minor" in fIndexvalues[i].
///
///  Once the index is computed, one can retrieve one entry via
/// ~~~{.cpp}
///     T->GetEntryWithIndex(majornumber, minornumber)
/// ~~~
/// Example:
/// ~~~{.cpp}
///  tree.BuildIndex("Run","Event"); //creates an index using leaves Run and Event
///  tree.GetEntryWithIndex(1234,56789); // reads entry corresponding to
///                                      // Run=1234 and Event=56789
/// ~~~
/// Note that majorname and minorname may be expressions using original
/// Tree variables eg: "run-90000", "event +3*xx". However the result
/// must be integer.
///
/// In case an expression is specified, the equivalent expression must be computed
/// when calling GetEntryWithIndex.
///
/// To build an index with only majorname, specify minorname="0" (default)
///
/// ## TreeIndex and Friend Trees
///
/// Assuming a parent Tree T and a friend Tree TF, the following cases are supported:
/// -  CASE 1: T->GetEntry(entry) is called
///            In this case, the serial number entry is used to retrieve
///            the data in both Trees.
/// -  CASE 2: T->GetEntry(entry) is called, TF has a TreeIndex
///            the expressions given in major/minorname of TF are used
///            to compute the value pair major,minor with the data in T.
///         TF->GetEntryWithIndex(major,minor) is then called (tricky case!)
/// -  CASE 3: T->GetEntryWithIndex(major,minor) is called.
///            It is assumed that both T and TF have a TreeIndex built using
///            the same major and minor name.
///
/// ## Saving the TreeIndex
///
/// Once the index is built, it can be saved with the TTree object
/// with tree.Write(); (if the file has been open in "update" mode).
///
/// The most convenient place to create the index is at the end of
/// the filling process just before saving the Tree header.
/// If a previous index was computed, it is redefined by this new call.
///
/// Note that this function can also be applied to a TChain.
///
/// The return value is the number of entries in the Index (< 0 indicates failure)
///
/// It is possible to play with different TreeIndex in the same Tree.
/// see comments in TTree::SetTreeIndex.

TTreeIndex::TTreeIndex(const TTree *T, const char *majorname, const char *minorname)
           : TVirtualIndex()
{
   fTree               = (TTree*)T;
   fN                  = 0;
   fIndexValues        = 0;
   fIndexValuesMinor   = 0;
   fIndex              = 0;
   fMajorFormula       = 0;
   fMinorFormula       = 0;
   fMajorFormulaParent = 0;
   fMinorFormulaParent = 0;
   fMajorName          = majorname;
   fMinorName          = minorname;
   if (!T) return;
   fN = T->GetEntries();
   if (fN <= 0) {
      MakeZombie();
      Error("TreeIndex","Cannot build a TreeIndex with a Tree having no entries");
      return;
   }

   GetMajorFormula();
   GetMinorFormula();
   if (!fMajorFormula || !fMinorFormula) {
      MakeZombie();
      Error("TreeIndex","Cannot build the index with major=%s, minor=%s",fMajorName.Data(), fMinorName.Data());
      return;
   }
   if ((fMajorFormula->GetNdim() != 1) || (fMinorFormula->GetNdim() != 1)) {
      MakeZombie();
      Error("TreeIndex","Cannot build the index with major=%s, minor=%s",fMajorName.Data(), fMinorName.Data());
      return;
   }
   // accessing array elements should be OK
   //if ((fMajorFormula->GetMultiplicity() != 0) || (fMinorFormula->GetMultiplicity() != 0)) {
   //   MakeZombie();
   //   Error("TreeIndex","Cannot build the index with major=%s, minor=%s that cannot be arrays",fMajorName.Data(), fMinorName.Data());
   //   return;
   //}

   Long64_t *tmp_major = new Long64_t[fN];
   Long64_t *tmp_minor = new Long64_t[fN];
   Long64_t i;
   Long64_t oldEntry = fTree->GetReadEntry();
   Int_t current = -1;
   for (i=0;i<fN;i++) {
      Long64_t centry = fTree->LoadTree(i);
      if (centry < 0) break;
      if (fTree->GetTreeNumber() != current) {
         current = fTree->GetTreeNumber();
         fMajorFormula->UpdateFormulaLeaves();
         fMinorFormula->UpdateFormulaLeaves();
      }
      tmp_major[i] = (Long64_t) fMajorFormula->EvalInstance<LongDouble_t>();
      tmp_minor[i] = (Long64_t) fMinorFormula->EvalInstance<LongDouble_t>();
   }
   fIndex = new Long64_t[fN];
   for(i = 0; i < fN; i++) { fIndex[i] = i; }
   std::sort(fIndex, fIndex + fN, IndexSortComparator(tmp_major, tmp_minor) );
   //TMath::Sort(fN,w,fIndex,0);
   fIndexValues = new Long64_t[fN];
   fIndexValuesMinor = new Long64_t[fN];
   for (i=0;i<fN;i++) {
      fIndexValues[i] = tmp_major[fIndex[i]];
      fIndexValuesMinor[i] = tmp_minor[fIndex[i]];
   }

   delete [] tmp_major;
   delete [] tmp_minor;
   fTree->LoadTree(oldEntry);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TTreeIndex::~TTreeIndex()
{
   if (fTree && fTree->GetTreeIndex() == this) fTree->SetTreeIndex(0);
   delete [] fIndexValues;      fIndexValues = 0;
   delete [] fIndexValuesMinor;      fIndexValuesMinor = 0;
   delete [] fIndex;            fIndex = 0;
   delete fMajorFormula;        fMajorFormula  = 0;
   delete fMinorFormula;        fMinorFormula  = 0;
   delete fMajorFormulaParent;  fMajorFormulaParent = 0;
   delete fMinorFormulaParent;  fMinorFormulaParent = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Append 'add' to this index.  Entry 0 in add will become entry n+1 in this.
/// If delaySort is true, do not sort the value, then you must call
/// Append(0,kFALSE);

void TTreeIndex::Append(const TVirtualIndex *add, Bool_t delaySort )
{

   if (add && add->GetN()) {
      // Create new buffer (if needed)

      const TTreeIndex *ti_add = dynamic_cast<const TTreeIndex*>(add);
      if (ti_add == 0) {
         Error("Append","Can only Append a TTreeIndex to a TTreeIndex but got a %s",
               add->IsA()->GetName());
      }

      Long64_t oldn = fN;
      fN += add->GetN();

      Long64_t *oldIndex = fIndex;
      Long64_t *oldValues = GetIndexValues();
      Long64_t *oldValues2 = GetIndexValuesMinor();

      fIndex = new Long64_t[fN];
      fIndexValues = new Long64_t[fN];
      fIndexValuesMinor = new Long64_t[fN];

      // Copy data
      Long_t size = sizeof(Long64_t) * oldn;
      Long_t add_size = sizeof(Long64_t) * add->GetN();

      memcpy(fIndex,oldIndex, size);
      memcpy(fIndexValues,oldValues, size);
      memcpy(fIndexValuesMinor,oldValues2, size);

      Long64_t *addIndex = ti_add->GetIndex();
      Long64_t *addValues = ti_add->GetIndexValues();
      Long64_t *addValues2 = ti_add->GetIndexValuesMinor();

      memcpy(fIndex + oldn, addIndex, add_size);
      memcpy(fIndexValues + oldn, addValues, add_size);
      memcpy(fIndexValuesMinor + oldn, addValues2, add_size);
      for(Int_t i = 0; i < add->GetN(); i++) {
         fIndex[oldn + i] += oldn;
      }

      delete [] oldIndex;
      delete [] oldValues;
      delete [] oldValues2;
   }

   // Sort.
   if (!delaySort) {
      Long64_t *addValues = GetIndexValues();
      Long64_t *addValues2 = GetIndexValuesMinor();
      Long64_t *ind = fIndex;
      Long64_t *conv = new Long64_t[fN];

      for(Long64_t i = 0; i < fN; i++) { conv[i] = i; }
      std::sort(conv, conv+fN, IndexSortComparator(addValues, addValues2) );
      //Long64_t *w = fIndexValues;
      //TMath::Sort(fN,w,conv,0);

      fIndex = new Long64_t[fN];
      fIndexValues = new Long64_t[fN];
      fIndexValuesMinor = new Long64_t[fN];

      for (Int_t i=0;i<fN;i++) {
         fIndex[i] = ind[conv[i]];
         fIndexValues[i] = addValues[conv[i]];
         fIndexValuesMinor[i] = addValues2[conv[i]];
      }
      delete [] addValues;
      delete [] addValues2;
      delete [] ind;
      delete [] conv;
   }
}



////////////////////////////////////////////////////////////////////////////////
/// conversion from old 64bit indexes
/// return true if index was converted

bool TTreeIndex::ConvertOldToNew()
{
   if( !fIndexValuesMinor && fN ) {
      fIndexValuesMinor = new Long64_t[fN];
      for(int i=0; i<fN; i++) {
         fIndexValuesMinor[i] = (fIndexValues[i] & 0x7fffffff);
         fIndexValues[i] >>= 31;
      }
      return true;
   }
   return false;
}



////////////////////////////////////////////////////////////////////////////////
/// Returns the entry number in this (friend) Tree corresponding to entry in
/// the master Tree 'parent'.
/// In case this (friend) Tree and 'master' do not share an index with the same
/// major and minor name, the entry serial number in the (friend) tree
/// and in the master Tree are assumed to be the same

Long64_t TTreeIndex::GetEntryNumberFriend(const TTree *parent)
{
   if (!parent) return -3;
   // We reached the end of the parent tree
   Long64_t pentry = parent->GetReadEntry();
   if (pentry >= parent->GetEntries())
      return -2;
   GetMajorFormulaParent(parent);
   GetMinorFormulaParent(parent);
   if (!fMajorFormulaParent || !fMinorFormulaParent) return -1;
   if (!fMajorFormulaParent->GetNdim() || !fMinorFormulaParent->GetNdim()) {
      // The Tree Index in the friend has a pair majorname,minorname
      // not available in the parent Tree T.
      // if the friend Tree has less entries than the parent, this is an error
      if (pentry >= fTree->GetEntries()) return -2;
      // otherwise we ignore the Tree Index and return the entry number
      // in the parent Tree.
      return pentry;
   }

   // majorname, minorname exist in the parent Tree
   // we find the current values pair majorv,minorv in the parent Tree
   Double_t majord = fMajorFormulaParent->EvalInstance();
   Double_t minord = fMinorFormulaParent->EvalInstance();
   Long64_t majorv = (Long64_t)majord;
   Long64_t minorv = (Long64_t)minord;
   // we check if this pair exist in the index.
   // if yes, we return the corresponding entry number
   // if not the function returns -1
   return fTree->GetEntryNumberWithIndex(majorv,minorv);
}


////////////////////////////////////////////////////////////////////////////////
/// find position where major|minor values are in the IndexValues tables
/// this is the index in IndexValues table, not entry# !
/// use lower_bound STD algorithm.

Long64_t TTreeIndex::FindValues(Long64_t major, Long64_t minor) const
{
   Long64_t mid, step, pos = 0, count = fN;
   // find lower bound using bisection
   while( count > 0 ) {
      step = count / 2;
      mid = pos + step;
      // check if *mid < major|minor
      if( fIndexValues[mid] < major
          || ( fIndexValues[mid] == major &&  fIndexValuesMinor[mid] < minor ) ) {
         pos = mid+1;
         count -= step + 1;
      } else
         count = step;
   }
   return pos;
}


////////////////////////////////////////////////////////////////////////////////
/// Return entry number corresponding to major and minor number.
/// Note that this function returns only the entry number, not the data
/// To read the data corresponding to an entry number, use TTree::GetEntryWithIndex
/// the BuildIndex function has created a table of Double_t* of sorted values
/// corresponding to val = major<<31 + minor;
/// The function performs binary search in this sorted table.
/// If it finds a pair that maches val, it returns directly the
/// index in the table.
/// If an entry corresponding to major and minor is not found, the function
/// returns the index of the major,minor pair immediatly lower than the
/// requested value, ie it will return -1 if the pair is lower than
/// the first entry in the index.
///
/// See also GetEntryNumberWithIndex

Long64_t TTreeIndex::GetEntryNumberWithBestIndex(Long64_t major, Long64_t minor) const
{
   if (fN == 0) return -1;

   Long64_t pos = FindValues(major, minor);
   if( pos < fN && fIndexValues[pos] == major && fIndexValuesMinor[pos] == minor )
      return fIndex[pos];
   if( --pos < 0 )
      return -1;
   return fIndex[pos];
}


////////////////////////////////////////////////////////////////////////////////
/// Return entry number corresponding to major and minor number.
/// Note that this function returns only the entry number, not the data
/// To read the data corresponding to an entry number, use TTree::GetEntryWithIndex
/// the BuildIndex function has created a table of Double_t* of sorted values
/// corresponding to val = major<<31 + minor;
/// The function performs binary search in this sorted table.
/// If it finds a pair that maches val, it returns directly the
/// index in the table, otherwise it returns -1.
///
/// See also GetEntryNumberWithBestIndex

Long64_t TTreeIndex::GetEntryNumberWithIndex(Long64_t major, Long64_t minor) const
{
   if (fN == 0) return -1;

   Long64_t pos = FindValues(major, minor);
   if( pos < fN && fIndexValues[pos] == major && fIndexValuesMinor[pos] == minor )
      return fIndex[pos];
   return -1;
}


////////////////////////////////////////////////////////////////////////////////

Long64_t* TTreeIndex::GetIndexValuesMinor()  const
{
   return fIndexValuesMinor;
}



////////////////////////////////////////////////////////////////////////////////
/// Return a pointer to the TreeFormula corresponding to the majorname.

TTreeFormula *TTreeIndex::GetMajorFormula()
{
   if (!fMajorFormula) {
      fMajorFormula = new TTreeFormula("Major",fMajorName.Data(),fTree);
      fMajorFormula->SetQuickLoad(kTRUE);
   }
   return fMajorFormula;
}

////////////////////////////////////////////////////////////////////////////////
/// Return a pointer to the TreeFormula corresponding to the minorname.

TTreeFormula *TTreeIndex::GetMinorFormula()
{
   if (!fMinorFormula) {
      fMinorFormula = new TTreeFormula("Minor",fMinorName.Data(),fTree);
      fMinorFormula->SetQuickLoad(kTRUE);
   }
   return fMinorFormula;
}

////////////////////////////////////////////////////////////////////////////////
/// Return a pointer to the TreeFormula corresponding to the majorname in parent tree.

TTreeFormula *TTreeIndex::GetMajorFormulaParent(const TTree *parent)
{
   if (!fMajorFormulaParent) {
      // Prevent TTreeFormula from finding any of the branches in our TTree even if it
      // is a friend of the parent TTree.
      TTree::TFriendLock friendlock(fTree, TTree::kFindLeaf | TTree::kFindBranch | TTree::kGetBranch | TTree::kGetLeaf);
      fMajorFormulaParent = new TTreeFormula("MajorP",fMajorName.Data(),const_cast<TTree*>(parent));
      fMajorFormulaParent->SetQuickLoad(kTRUE);
   }
   if (fMajorFormulaParent->GetTree() != parent) {
      fMajorFormulaParent->SetTree(const_cast<TTree*>(parent));
      fMajorFormulaParent->UpdateFormulaLeaves();
   }
   return fMajorFormulaParent;
}

////////////////////////////////////////////////////////////////////////////////
/// Return a pointer to the TreeFormula corresponding to the minorname in parent tree.

TTreeFormula *TTreeIndex::GetMinorFormulaParent(const TTree *parent)
{
   if (!fMinorFormulaParent) {
      // Prevent TTreeFormula from finding any of the branches in our TTree even if it
      // is a friend of the parent TTree.
      TTree::TFriendLock friendlock(fTree, TTree::kFindLeaf | TTree::kFindBranch | TTree::kGetBranch | TTree::kGetLeaf);
      fMinorFormulaParent = new TTreeFormula("MinorP",fMinorName.Data(),const_cast<TTree*>(parent));
      fMinorFormulaParent->SetQuickLoad(kTRUE);
   }
   if (fMinorFormulaParent->GetTree() != parent) {
      fMinorFormulaParent->SetTree(const_cast<TTree*>(parent));
      fMinorFormulaParent->UpdateFormulaLeaves();
   }
   return fMinorFormulaParent;
}

////////////////////////////////////////////////////////////////////////////////
/// Return kTRUE if index can be applied to the TTree

Bool_t TTreeIndex::IsValidFor(const TTree *parent)
{
   auto *majorFormula = GetMajorFormulaParent(parent);
   auto *minorFormula = GetMinorFormulaParent(parent);
   if ((majorFormula == nullptr || majorFormula->GetNdim() == 0) ||
       (minorFormula == nullptr || minorFormula->GetNdim() == 0))
         return kFALSE;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Print the table with : serial number, majorname, minorname.
/// -  if option = "10" print only the first 10 entries
/// -  if option = "100" print only the first 100 entries
/// -  if option = "1000" print only the first 1000 entries

void TTreeIndex::Print(Option_t * option) const
{
   TString opt = option;
   Bool_t printEntry = kFALSE;
   Long64_t n = fN;
   if (opt.Contains("10"))   n = 10;
   if (opt.Contains("100"))  n = 100;
   if (opt.Contains("1000")) n = 1000;
   if (opt.Contains("all")) {
      printEntry = kTRUE;
   }

   if (printEntry) {
      Printf("\n*****************************************************************");
      Printf("*    Index of Tree: %s/%s",fTree->GetName(),fTree->GetTitle());
      Printf("*****************************************************************");
      Printf("%8s : %16s : %16s : %16s","serial",fMajorName.Data(),fMinorName.Data(),"entry number");
      Printf("*****************************************************************");
      for (Long64_t i=0;i<n;i++) {
         Printf("%8lld :         %8lld :         %8lld :         %8lld",
                i, fIndexValues[i], GetIndexValuesMinor()[i], fIndex[i]);
      }

   } else {
      Printf("\n**********************************************");
      Printf("*    Index of Tree: %s/%s",fTree->GetName(),fTree->GetTitle());
      Printf("**********************************************");
      Printf("%8s : %16s : %16s","serial",fMajorName.Data(),fMinorName.Data());
      Printf("**********************************************");
      for (Long64_t i=0;i<n;i++) {
         Printf("%8lld :         %8lld :         %8lld",
                i, fIndexValues[i],GetIndexValuesMinor()[i]);
     }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TTreeIndex.
/// Note that this Streamer should be changed to an automatic Streamer
/// once TStreamerInfo supports an index of type Long64_t

void TTreeIndex::Streamer(TBuffer &R__b)
{
   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TVirtualIndex::Streamer(R__b);
      fMajorName.Streamer(R__b);
      fMinorName.Streamer(R__b);
      R__b >> fN;
      fIndexValues = new Long64_t[fN];
      R__b.ReadFastArray(fIndexValues,fN);
      if( R__v > 1 ) {
         fIndexValuesMinor = new Long64_t[fN];
         R__b.ReadFastArray(fIndexValuesMinor,fN);
      } else {
         ConvertOldToNew();
      }
      fIndex      = new Long64_t[fN];
      R__b.ReadFastArray(fIndex,fN);
      R__b.CheckByteCount(R__s, R__c, TTreeIndex::IsA());
   } else {
      R__c = R__b.WriteVersion(TTreeIndex::IsA(), kTRUE);
      TVirtualIndex::Streamer(R__b);
      fMajorName.Streamer(R__b);
      fMinorName.Streamer(R__b);
      R__b << fN;
      R__b.WriteFastArray(fIndexValues, fN);
      R__b.WriteFastArray(fIndexValuesMinor, fN);
      R__b.WriteFastArray(fIndex, fN);
      R__b.SetByteCount(R__c, kTRUE);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Called by TChain::LoadTree when the parent chain changes it's tree.

void TTreeIndex::UpdateFormulaLeaves(const TTree *parent)
{
   if (fMajorFormula)       { fMajorFormula->UpdateFormulaLeaves();}
   if (fMinorFormula)       { fMinorFormula->UpdateFormulaLeaves();}
   if (fMajorFormulaParent) {
      if (parent) fMajorFormulaParent->SetTree(const_cast<TTree*>(parent));
      fMajorFormulaParent->UpdateFormulaLeaves();
   }
   if (fMinorFormulaParent) {
      if (parent) fMinorFormulaParent->SetTree(const_cast<TTree*>(parent));
      fMinorFormulaParent->UpdateFormulaLeaves();
   }
}
////////////////////////////////////////////////////////////////////////////////
/// this function is called by TChain::LoadTree and TTreePlayer::UpdateFormulaLeaves
/// when a new Tree is loaded.
/// Because Trees in a TChain may have a different list of leaves, one
/// must update the leaves numbers in the TTreeFormula used by the TreeIndex.

void TTreeIndex::SetTree(const TTree *T)
{
   fTree = (TTree*)T;
}

