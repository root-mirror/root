// @(#)root/tree:$Id$
// Author: Anna Kreshuk 27/10/2006

/*************************************************************************
 * Copyright (C) 1995-2006, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TEntryList
\ingroup tree

A List of entry numbers in a TTree or TChain.

There are two types of entry lists:

#### 1.
    for a TTree (fBlocks data member is non-zero)
    Entry numbers are stored in TEntryListBlocks, which, in their turn, are stored
    in the TObjArray fBlocks. The range of the entry numbers is cut into intervals
    of kBlockSize entries (currently 64000), so that the first block contains
    information which entries out of the first 64000 pass the selection, the second
    block - which entries out of the 64000-127999 interval pass the selection, etc.
    Some blocks, obviously, might be empty. The internal representation of entry
    numbers in the blocks is described in the TEntryListBlock class description, and
    this representation might be changed by calling OptimizeStorage() function
    (when the list is filled via the Enter() function, this is done automatically,
    except for the last block).
    Individual entry lists can be merged (functions Merge() and Add())
    to make an entry list for a TChain of corresponding TTrees.
Begin_Macro
entrylist_figure1.C
End_Macro

#### 2.
    for a TChain (fLists data member is non-zero)
    It contains a TList of sub-lists (TEntryList objects, corresponding to each TTree)
    Trees and lists are matched by the TTree name and its file name (full path).
    All sub-lists are returned by the GetLists() function and individual lists are
    returned by GetEntryList() function. Such lists are no different from the lists for
    TTrees, described above.
Begin_Macro
entrylist_figure2.C
End_Macro

## Operations on entry lists

- __Add__() - if the lists are for the same tree, adds all the entries of the second list
          to the first list. If the lists are for different trees, creates a TEntryList
          with 2 sublists for each TTree. If the lists are for TChains, merges the ones
          for the same trees and adds new sublists for the TTrees that were not included
          in the first TEntryList
- __Subtract__() - if the lists are for the same TTree, removes the entries of the second
               list from the first list. If the lists are for TChains, loops over all
               sub-lists
- __GetEntry(n)__ - returns the n-th entry number
- __Next__()      - returns next entry number. Note, that this function is
                much faster than GetEntry, and it's called when GetEntry() is called
                for 2 or more indices in a row.

## TTree::Draw() and TChain::Draw()

Use option __entrylist__ to write the results of TTree::Draw and TChain::Draw into
an entry list. Example:
~~~ {.cpp}
     tree->Draw(">>elist", "x<0 && y>0", "entrylist");
     TEntryList *elist = (TEntryList*)gDirectory->Get("elist");
~~~
## Example of Loop on TEntryList with a TChain
~~~ {.cpp}
     void loopChain() {
        TFile *fe = TFile::Open("myelist.root");
        TEntryList *myelist = (TEntryList*)fe->Get("myelist");
        TChain *chain = new TChain("ntuple");
        chain->Add("hsimple.root");
        chain->Add("hsimple2.root");
        Long64_t listEntries = myelist->GetN();
        Long64_t chainEntries = chain->GetEntries();
        Int_t treenum = 0;
        chain->SetEntryList(myelist);

        for (entry=start;entry < end;entry++) {
           entryNumber = chain->GetEntryNumber(entry);
           if (entryNumber < 0) break;
           localEntry = chain->LoadTree(entryNumber);
           if (localEntry < 0) break;
           ....
           then either call
              branch->GetEntry(localEntry);
           or
              chain->GetEntry(entryNumber);
           In the later case the LoadTree is then somewhat redudant.
           ...
        }
     }
~~~
When using the TEntryList interface directly, you can get the 'tree number' and entry in
the current tree (i.e. value similar to the return value of LoadTree) from calling
TEntryList::GetEntryAndTree:
~~~ {.cpp}
     Long64_t treeEntry = myelist->GetEntryAndTree(el,treenum);
~~~
to obtain the entry number within the chain you need to add to it the value of
`treeEntry+ch->GetTreeOffset()[treenum]`
such that the loop in the previous example can also be written as:
~~~ {.cpp}
     for (Long64_t el = 0; el < listEntries; el++) {
        Long64_t treeEntry = myelist->GetEntryAndTree(el,treenum);
        Long64_t chainEntry = treeEntry+ch->GetTreeOffset()[treenum];
        printf("el=%lld, treeEntry=%lld, chainEntry=%lld, treenum=%d\n", el, treeEntry, chainEntry, treenum);

        ch->LoadTree(chainEntry); // this also returns treeEntry
        needed_branch->GetEntry(treeEntry);
     }
~~~
## TSelectors

To fill an TEntryList from a TSelector correctly, one must add the TEntryList object
to the output list of the selector (TSelector::fOutput). This is the only way to
make the sub-lists of the TEntryList switch when the current tree of the TChain is
changed.

## Using a TEntryList as input (TTree::SetEntryList() and TChain::SetEntryList())

while the TTree::SetEntryList() function is only setting the TTree::fEntryList
data member, the same function in TChain also finds correspondence between
the TTrees of this TChain and the sub-lists of this TEntryList.

## TEntryList and the current directory

TEntryList objects are automatically added to the current directory (like TTrees).
However, in case of a TEntryList for a chain, only the top-level entry list is added,
not the sub-lists for specific trees. Placing entry lists in the current directory
allows calling them as a part of a TTreeFormula expression, so if the user wants
to extract a sublist from a TChain entry list via the GetEntryList() or some other
function, they have to add it to the current directory to be able to use it in
TTreeFormula expressions.

## TEntryList and TEventList

TTree::SetEventList() and TChain::SetEventList() transform a TEventList into a TEntryList
See comments to those functions for more details
*/

#include "TEntryList.h"
#include "TEntryListBlock.h"
#include "TError.h"
#include "TKey.h"
#include "TBuffer.h"
#include "TTree.h"
#include "TFile.h"
#include "TRegexp.h"
#include "TSystem.h"
#include "TObjString.h"

ClassImp(TEntryList);

////////////////////////////////////////////////////////////////////////////////
/// default c-tor

TEntryList::TEntryList() : fEntriesToProcess(0)
{
   fLists = 0;
   fCurrent = 0;
   fBlocks = 0;
   fN = 0;
   fNBlocks = 0;
   fTreeName = "";
   fFileName = "";
   fStringHash = 0;
   fTreeNumber = -1;
   fDirectory = 0;
   fReapply = kFALSE;
   fLastIndexQueried = -1;
   fLastIndexReturned = 0;
   fShift = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// c-tor with name and title

TEntryList::TEntryList(const char *name, const char *title) :
   TNamed(name, title),
   fEntriesToProcess(0)
{
   fLists = 0;
   fCurrent = 0;
   fBlocks = 0;
   fN = 0;
   fNBlocks = 0;
   fTreeName = "";
   fFileName = "";
   fStringHash = 0;
   fTreeNumber = -1;
   fReapply = kFALSE;

   fDirectory  = gDirectory;
   if (fDirectory) fDirectory->Append(this);

   fLastIndexQueried = -1;
   fLastIndexReturned = 0;
   fShift = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// constructor with name and title, which also sets the tree

TEntryList::TEntryList(const char *name, const char *title, const TTree *tree):TNamed(name, title)
{
   fLists = 0;
   fCurrent = 0;
   fBlocks = 0;
   fN = 0;
   fNBlocks = 0;
   fTreeNumber = -1;
   TEntryList::SetTree(tree);
   fReapply = kFALSE;

   fDirectory  = gDirectory;
   if (fDirectory) fDirectory->Append(this);

   fLastIndexQueried = -1;
   fLastIndexReturned = 0;
   fShift = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// c-tor with name and title, which also sets the treename and the filename

TEntryList::TEntryList(const char *name, const char *title, const char *treename, const char *filename) : TNamed(name, title),fEntriesToProcess(0)
{
   fLists = 0;
   fCurrent = 0;
   fBlocks = 0;
   fNBlocks = 0;
   fN = 0;
   SetTree(treename, filename);
   fTreeNumber = -1;
   fReapply = kFALSE;

   fDirectory  = gDirectory;
   if (fDirectory) fDirectory->Append(this);

   fLastIndexQueried = -1;
   fLastIndexReturned = 0;
   fShift = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// c-tor, which sets the tree

TEntryList::TEntryList(const TTree *tree) : fEntriesToProcess(0)
{
   fLists = 0;
   fCurrent = 0;
   fBlocks = 0;
   fNBlocks = 0;
   fN = 0;

   SetTree(tree);
   fTreeNumber = -1;

   fReapply = kFALSE;
   fDirectory  = gDirectory;
   if (fDirectory) fDirectory->Append(this);

   fLastIndexQueried = -1;
   fLastIndexReturned = 0;
   fShift = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// copy c-tor

TEntryList::TEntryList(const TEntryList &elist) : TNamed(elist)
{
   fNBlocks = elist.fNBlocks;
   fTreeName = elist.fTreeName;
   fFileName = elist.fFileName;
   fStringHash = elist.fStringHash;
   fTreeNumber = elist.fTreeNumber;
   fLastIndexQueried = -1;
   fLastIndexReturned = 0;
   fN = elist.fN;
   fShift = elist.fShift;
   fLists = 0;
   fBlocks = 0;
   fReapply = elist.fReapply;
   fCurrent = 0;
   fEntriesToProcess = elist.fEntriesToProcess;
   if (elist.fLists){
      fLists = new TList();
      TEntryList *el1 = 0;
      TEntryList *el2 = 0;
      TIter next(elist.fLists);
      while((el1 = (TEntryList*)next())){
         el2 = new TEntryList(*el1);
         if (el1==elist.fCurrent)
            fCurrent = el2;
         fLists->Add(el2);
      }
   } else {
      if (elist.fBlocks){
         TEntryListBlock *block1 = 0;
         TEntryListBlock *block2 = 0;
         //or just copy it as a TObjArray??
         fBlocks = new TObjArray();
         for (Int_t i=0; i<fNBlocks; i++){
            block1 = (TEntryListBlock*)elist.fBlocks->UncheckedAt(i);
            block2 = new TEntryListBlock(*block1);
            fBlocks->Add(block2);
         }
      }
      fCurrent = this;
   }
   fDirectory  = 0;

}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TEntryList::~TEntryList()
{
   if (fBlocks){
      fBlocks->Delete();
      delete fBlocks;
   }
   fBlocks = 0;
   if (fLists){
      fLists->Delete();
      delete fLists;
   }

   fLists = 0;

   if (fDirectory) fDirectory->Remove(this);
   fDirectory  = 0;

}

////////////////////////////////////////////////////////////////////////////////
/// Add 2 entry lists

void TEntryList::Add(const TEntryList *elist)
{
   if (fN==0){
      if (!fLists && fTreeName=="" && fFileName==""){
         //this list is empty. copy the other list completely
         fNBlocks = elist->fNBlocks;
         fTreeName = elist->fTreeName;
         fFileName = elist->fFileName;
         fStringHash = elist->fStringHash;
         fTreeNumber = elist->fTreeNumber;
         fLastIndexQueried = -1;
         fLastIndexReturned = 0;
         fN = elist->fN;
         if (elist->fLists){
            fLists = new TList();
            TEntryList *el1 = 0;
            TEntryList *el2 = 0;
            TIter next(elist->fLists);
            while((el1 = (TEntryList*)next())){
               el2 = new TEntryList(*el1);
               if (el1==elist->fCurrent)
                  fCurrent = el2;
               fLists->Add(el2);
            }
         } else {
            if (elist->fBlocks){
               TEntryListBlock *block1 = 0;
               TEntryListBlock *block2 = 0;
               fBlocks = new TObjArray();
               for (Int_t i=0; i<fNBlocks; i++){
                  block1 = (TEntryListBlock*)elist->fBlocks->UncheckedAt(i);
                  block2 = new TEntryListBlock(*block1);
                  fBlocks->Add(block2);
               }
            }
            fCurrent = 0;
         }
         return;
      }
   }

   if (!fLists){
      if (!elist->fLists){
         if (!strcmp(elist->fTreeName.Data(),fTreeName.Data()) && !strcmp(elist->fFileName.Data(),fFileName.Data())){
            //entry lists are for the same tree
            if (!elist->fBlocks)
               //the other list is empty list
               return;
            if (!fBlocks){
               //this entry list is empty
               TEntryListBlock *block1 = 0;
               TEntryListBlock *block2 = 0;
               fNBlocks = elist->fNBlocks;
               fN = elist->fN;
               fBlocks = new TObjArray();
               for (Int_t i=0; i<fNBlocks; i++){
                  block1 = (TEntryListBlock*)elist->fBlocks->UncheckedAt(i);
                  block2 = new TEntryListBlock(*block1);
               fBlocks->Add(block2);
               }
               return;
            }
            //both not empty, merge block by block
            TEntryListBlock *block1=0;
            TEntryListBlock *block2=0;
            Int_t i;
            Int_t nmin = TMath::Min(fNBlocks, elist->fNBlocks);
            Long64_t nnew, nold;
            for (i=0; i<nmin; i++){
               block1 = (TEntryListBlock*)fBlocks->UncheckedAt(i);
               block2 = (TEntryListBlock*)elist->fBlocks->UncheckedAt(i);
               nold = block1->GetNPassed();
               nnew = block1->Merge(block2);
               fN = fN - nold + nnew;
            }
            if (fNBlocks<elist->fNBlocks){
               Int_t nmax = elist->fNBlocks;
               for (i=nmin; i<nmax; i++){
                  block2 = (TEntryListBlock*)elist->fBlocks->UncheckedAt(i);
                  block1 = new TEntryListBlock(*block2);
                  fBlocks->Add(block1);
                  fN+=block1->GetNPassed();
                  fNBlocks++;
               }
            }
            fLastIndexQueried = -1;
            fLastIndexReturned = 0;
         } else {
            //entry lists are for different trees. create a chain entry list with
            //2 sub lists for the first and second entry lists
            fLastIndexQueried = -1;
            fLastIndexReturned = 0;
            fLists = new TList();
            TEntryList *el = new TEntryList();
            el->fTreeName = fTreeName;
            el->fFileName = fFileName;
            el->fBlocks = fBlocks;
            fBlocks = 0;
            el->fNBlocks = fNBlocks;
            el->fN = fN;
            el->fLastIndexQueried = -1;
            el->fLastIndexReturned = 0;
            fLists->Add(el);
            el = new TEntryList(*elist);
            el->fLastIndexQueried = -1;
            el->fLastIndexReturned = 0;
            fLists->Add(el);
            fN+=el->GetN();
            fCurrent = 0;
         }
      } else {
         //second list already has sublists. add one by one
         TEntryList *el = 0;
         TIter next(elist->fLists);
         while ((el = (TEntryList*)next())){
            Add(el);
         }
         fCurrent = 0;
      }
   } else {
      //there are already some sublists in this list, just add another one
      if (!elist->fLists){
         //the other list doesn't have sublists
         TIter next(fLists);
         TEntryList *el = 0;
         Bool_t found = kFALSE;
         while ((el = (TEntryList*)next())){
            if (!strcmp(el->fTreeName.Data(), elist->fTreeName.Data()) &&
                !strcmp(el->fFileName.Data(), elist->fFileName.Data())){
            // if (el->fStringHash == elist->fStringHash){
               //found a list for the same tree
               Long64_t oldn = el->GetN();
               el->Add(elist);
               found = kTRUE;
               fN = fN - oldn + el->GetN();
               break;
            }
         }
         if (!found){
            el = new TEntryList(*elist);
            el->fLastIndexQueried = -1;
            el->fLastIndexReturned = 0;
            fLists->Add(el);
            fN+=el->GetN();
         }
      } else {
         //add all sublists from the other list
         TEntryList *el = 0;
         TIter next(elist->fLists);
         while ((el = (TEntryList*)next())){
            Add(el);
         }
         fCurrent = 0;
      }
      if (fCurrent){
         if (fCurrent->fBlocks){
            Int_t currentblock = (fCurrent->fLastIndexReturned)/kBlockSize;
            TEntryListBlock *block = (TEntryListBlock*)fCurrent->fBlocks->UncheckedAt(currentblock);
            block->ResetIndices();
            fCurrent->fLastIndexReturned = 0;
            fCurrent->fLastIndexQueried = -1;
         }
      }
      fCurrent = 0;
   }

}

////////////////////////////////////////////////////////////////////////////////
/// - When tree = 0, returns from the current list
/// - When tree != 0, finds the list, corresponding to this tree
/// - When tree is a chain, the entry is assumed to be global index and the local
/// entry is recomputed from the treeoffset information of the chain

Int_t TEntryList::Contains(Long64_t entry, TTree *tree)
{
   if (!tree){
      if (fBlocks) {
         //this entry list doesn't contain any sub-lists
         TEntryListBlock *block = 0;
         Int_t nblock = entry/kBlockSize;
         if (nblock >= fNBlocks) return 0;
         block = (TEntryListBlock*)fBlocks->UncheckedAt(nblock);
         return block->Contains(entry-nblock*kBlockSize);
      }
      if (fLists) {
         if (!fCurrent) fCurrent = (TEntryList*)fLists->First();
         return fCurrent->Contains(entry);
      }
      return 0;
   } else {
      Long64_t localEntry = tree->LoadTree(entry);
      SetTree(tree->GetTree());
      if (fCurrent)
         return fCurrent->Contains(localEntry);
   }
   return 0;

}

////////////////////////////////////////////////////////////////////////////////
/// Called by TKey and others to automatically add us to a directory when we are read from a file.

void TEntryList::DirectoryAutoAdd(TDirectory* dir)
{
   SetDirectory(dir);
}

////////////////////////////////////////////////////////////////////////////////
/// Add entry \#entry to the list
/// - When tree = 0, adds to the current list
/// - When tree != 0, finds the list, corresponding to this tree
/// - When tree is a chain, the entry is assumed to be global index and the local
/// entry is recomputed from the treeoffset information of the chain

Bool_t TEntryList::Enter(Long64_t entry, TTree *tree)
{
   if (!tree){
      if (!fLists) {
         if (!fBlocks) fBlocks = new TObjArray();
         TEntryListBlock *block = 0;
         Long64_t nblock = entry/kBlockSize;
         if (nblock >= fNBlocks) {
            if (fNBlocks>0){
               block = (TEntryListBlock*)fBlocks->UncheckedAt(fNBlocks-1);
               if (!block) return 0;
               block->OptimizeStorage();
            }
            for (Int_t i=fNBlocks; i<=nblock; i++){
               block = new TEntryListBlock();
               fBlocks->Add(block);
            }
         fNBlocks = nblock+1;
         }
         block = (TEntryListBlock*)fBlocks->UncheckedAt(nblock);
         if (block->Enter(entry-nblock*kBlockSize)) {
            fN++;
            return 1;
         }
      } else {
         //the entry in the current entry list
         if (!fCurrent) fCurrent = (TEntryList*)fLists->First();
         if (fCurrent->Enter(entry)) {
            if (fLists)
               fN++;
            return 1;
         }
      }
   } else {
      Long64_t localentry = tree->LoadTree(entry);
      SetTree(tree->GetTree());
      if (fCurrent){
         if (fCurrent->Enter(localentry)) {
            if (fLists)
               fN++;
            return 1;
         }
      }
   }
   return 0;

}

////////////////////////////////////////////////////////////////////////////////
/// Remove entry \#entry from the list
/// - When tree = 0, removes from the current list
/// - When tree != 0, finds the list, corresponding to this tree
/// - When tree is a chain, the entry is assumed to be global index and the local
/// entry is recomputed from the treeoffset information of the chain

Bool_t TEntryList::Remove(Long64_t entry, TTree *tree)
{
   if (entry < 0)
     return kFALSE;
   if (!tree) {
      if (!fLists) {
         if (!fBlocks) return 0;
         TEntryListBlock *block = 0;
         Long64_t nblock = entry/kBlockSize;
         block = (TEntryListBlock*)fBlocks->UncheckedAt(nblock);
         if (!block) return 0;
         Long64_t blockindex = entry - nblock*kBlockSize;
         if (block->Remove(blockindex)){
            fN--;
            return 1;
         }
      } else {
         if (!fCurrent) fCurrent = (TEntryList*)fLists->First();
         if (fCurrent->Remove(entry)){
            if (fLists)
               fN--;
            return 1;
         }
      }
   } else {
      Int_t localentry = tree->LoadTree(entry);
      SetTree(tree->GetTree());
      if (fCurrent){
         if (fCurrent->Remove(localentry)) {
            if (fLists)
               fN--;
            return 1;
         }
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the number of the entry \#index of this TEntryList in the TTree or TChain
/// See also Next().

Long64_t TEntryList::GetEntry(Int_t index)
{

   if (index>=fN){
      return -1;
   }
   if (index==fLastIndexQueried+1){
      //in a loop
      return Next();
   } else {
      if (fBlocks) {
         TEntryListBlock *block = 0;
         Long64_t total_passed = 0;
         Int_t i=0;
         while (total_passed<=index && i<fNBlocks){
            block=(TEntryListBlock*)fBlocks->UncheckedAt(i);
            total_passed+=block->GetNPassed();
            i++;
         }
         i--;
         total_passed-=block->GetNPassed();
         if (i!=fLastIndexReturned/kBlockSize){
            block = (TEntryListBlock*)fBlocks->UncheckedAt(fLastIndexReturned/kBlockSize);
            block->ResetIndices();
            block = (TEntryListBlock*)fBlocks->UncheckedAt(i);
         }

         Long64_t localindex = index - total_passed;
         Long64_t blockindex = block->GetEntry(localindex);
         if (blockindex < 0) return -1;
         Long64_t res = i*kBlockSize + blockindex;
         fLastIndexQueried = index;
         fLastIndexReturned = res;
         return res;
      } else {
         //find the corresponding list
         if (!fCurrent) fCurrent = (TEntryList*)fLists->First();
         TIter next(fLists);
         TEntryList *templist;
         Long64_t ntotal = 0;
         if (fCurrent){
            //reset all indices of the current list
            if (fCurrent->fBlocks){
               Int_t currentblock = (fCurrent->fLastIndexReturned)/kBlockSize;
               TEntryListBlock *block = (TEntryListBlock*)fCurrent->fBlocks->UncheckedAt(currentblock);
               block->ResetIndices();
               fCurrent->fLastIndexReturned = 0;
               fCurrent->fLastIndexQueried = -1;
            }
         }
         while ((templist = (TEntryList*)next())){
            if (!fShift){
               ntotal += templist->GetN();
            } else {
               if (templist->GetTreeNumber() >= 0)
                  ntotal += templist->GetN();
            }
            if (ntotal > index)
               break;
         }
         fCurrent = templist;
         if (!fCurrent) return -1;
         Long64_t localentry = index - (ntotal - fCurrent->GetN());
         fLastIndexQueried = index;
         fLastIndexReturned = fCurrent->GetEntry(localentry);
         return fLastIndexReturned;
      }

   }
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the index of "index"-th non-zero entry in the TTree or TChain
/// and the # of the corresponding tree in the chain

Long64_t TEntryList::GetEntryAndTree(Int_t index, Int_t &treenum)
{
//If shift is true, then when the requested entry is found in an entry list,
//for which there is no corresponding tree in the chain, this list is not
//taken into account, and entry from the next list with a tree is returned.
//Example:
//First sublist - 20 entries, second sublist - 5 entries, third sublist - 10 entries
//Second sublist doesn't correspond to any trees of the chain
//Then, when GetEntryAndTree(21, treenum, kTRUE) is called, first entry of the
//third sublist will be returned

   Long64_t result = GetEntry(index);
   if (result < 0) {
      treenum = -1;
      return result;
   }
   R__ASSERT(fLists == nullptr || (fLists != nullptr && fCurrent != nullptr));
   if (fCurrent)
      treenum = fCurrent->fTreeNumber;
   else
      treenum = fTreeNumber;
   if (treenum < 0)
      return -1;

   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// To be able to re-localize the entry-list we identify the file by just the
/// name and the anchor, i.e. we drop protocol, host, options, ...
/// The result in the form 'file#anchor' (or 'file', if no anchor is present)
/// is saved in 'fn'.
/// The function optionally (is 'local' is defined) checks file locality (i.e.
/// protocol 'file://') returning the result in '*local' .

void TEntryList::GetFileName(const char *filename, TString &fn, Bool_t *local)
{
   TUrl u(filename, kTRUE);
   if (local) *local = (!strcmp(u.GetProtocol(), "file")) ? kTRUE : kFALSE;
   if (strlen(u.GetAnchor()) > 0) {
      fn.Form("%s#%s", u.GetFile(), u.GetAnchor());
   } else {
      fn = u.GetFile();
   }
   // Done
   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the entry list, corresponding to treename and filename
/// By default, the filename is first tried as is, and then, if the corresponding list
/// is not found, the filename is expanded to the absolute path, and compared again.
/// To avoid it, use option "ne"

TEntryList *TEntryList::GetEntryList(const char *treename, const char *filename, Option_t *opt)
{
   if (gDebug > 1)
      Info("GetEntryList","tree: %s, file: %s",
                          (treename ? treename : "-"), (filename ? filename : "-"));

   if (!treename || !filename) return 0;
   TString option = opt;
   option.ToUpper();
   Bool_t nexp = option.Contains("NE");

   TString fn;
   Bool_t local;
   GetFileName(filename, fn, &local);
   if (nexp) local = kFALSE;

   if (gDebug > 1)
      Info("GetEntryList", "file: %s, local? %d", filename, local);

   if (!fLists){
      //there are no sublists
      if (!strcmp(treename, fTreeName.Data()) && !(strcmp(fn.Data(), fFileName.Data()))){
         return this;
      } else {
         //if the file is local, try the full name, unless "ne" option was specified
         if (!nexp && local){
            gSystem->ExpandPathName(fn);
            if (!gSystem->IsAbsoluteFileName(fn))
               gSystem->PrependPathName(gSystem->pwd(), fn);
            fn = gSystem->UnixPathName(fn);
            if (!strcmp(treename, fTreeName.Data()) && !(strcmp(fn.Data(), fFileName.Data())))
               return this;
         }
         return 0;
      }
   }

   TString stotal = treename;
   stotal.Append(fn);
   ULong_t newhash = stotal.Hash();

   TIter next(fLists);
   TEntryList *templist;
   while ((templist = (TEntryList*)next())){
      if (templist->fStringHash==0){
         stotal = templist->fTreeName + templist->fFileName;
         templist->fStringHash = stotal.Hash();
      }
      if (gDebug > 1)
         Info("GetEntryList", "file: %s (fn: %s), hash: %lu, element hash: %lu",
                              filename, fn.Data(), newhash, templist->fStringHash);
      if (newhash == templist->fStringHash){
         if (!strcmp(templist->GetTreeName(), treename) && !strcmp(templist->GetFileName(), fn.Data())){
            return templist;
         }
      }
   }

   //didn't find anything for this filename, try the full name too
   if (!nexp && local){
      TString longname = fn;
      gSystem->ExpandPathName(longname);
      if (!gSystem->IsAbsoluteFileName(longname))
         gSystem->PrependPathName(gSystem->pwd(), longname);
      longname = gSystem->UnixPathName(longname);
      stotal = treename;
      stotal.Append(longname);
      newhash = stotal.Hash();
      next.Reset();
      while ((templist = (TEntryList*)next())){
         if (templist->fStringHash==0){
            stotal = templist->fTreeName + templist->fFileName;
            templist->fStringHash = stotal.Hash();
         }
         if (gDebug > 1)
            Info("GetEntryList", "file: %s (longname: %s), hash: %lu, element hash: %lu",
                                 filename, longname.Data(), newhash, templist->fStringHash);
         if (newhash == templist->fStringHash){
            if (templist->fTreeName == treename && templist->fFileName == longname){
               return templist;
            }
         }
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Merge this list with the lists from the collection

Int_t TEntryList::Merge(TCollection *list)
{
   if (!list) return -1;
   TIter next(list);
   TEntryList *elist = 0;
   while ((elist = (TEntryList*)next())) {
      if (!elist->InheritsFrom(TEntryList::Class())) {
         Error("Add","Attempt to add object of class: %s to a %s",elist->ClassName(),this->ClassName());
         return -1;
      }
      Add(elist);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the next non-zero entry index (next after fLastIndexQueried)
/// this function is faster than GetEntry()

Long64_t TEntryList::Next()
{
   Long64_t result;
   if (fN == fLastIndexQueried+1 || fN==0){
      return -1;
   }
   if (fBlocks){
      Int_t iblock = fLastIndexReturned/kBlockSize;
      TEntryListBlock *current_block = (TEntryListBlock*)fBlocks->UncheckedAt(iblock);
      result = current_block->Next();
      if (result>=0) {
         fLastIndexQueried++;
         fLastIndexReturned = result+kBlockSize*iblock;
         return fLastIndexReturned;
      }
      else {
         while (result<0 && iblock<fNBlocks-1) {
            current_block->ResetIndices();
            iblock++;
            current_block = (TEntryListBlock*)fBlocks->UncheckedAt(iblock);
            current_block->ResetIndices();
            result = current_block->Next();
         }
         if (result<0) {
            fLastIndexQueried = -1;
            fLastIndexReturned = 0;
            return -1;
         }
         fLastIndexQueried++;
         fLastIndexReturned = result+kBlockSize*iblock;

         return fLastIndexReturned;
      }
   } else {
      if (!fCurrent) {
         fCurrent = (TEntryList*)fLists->First();
         if (!fCurrent) return 0;
         if (fShift) {
            while (fCurrent->GetTreeNumber()<0) {
               fCurrent = (TEntryList*)fLists->After(fCurrent);
               if (!fCurrent) return 0;
            }
         }
      }
      result = fCurrent->Next();
      if (result>=0) {
         fLastIndexQueried++;
         fLastIndexReturned = result;
         return result;
      } else {
         if (fCurrent){
            //reset all indices of the current list
            if (fCurrent->fBlocks){
               Int_t currentblock = (fCurrent->fLastIndexReturned)/kBlockSize;
               TEntryListBlock *block = (TEntryListBlock*)fCurrent->fBlocks->UncheckedAt(currentblock);
               block->ResetIndices();
               fCurrent->fLastIndexReturned = 0;
               fCurrent->fLastIndexQueried = -1;
            }
         }

         //find the list with the next non-zero entry
         while (result<0 && fCurrent!=((TEntryList*)fLists->Last())){
            if (!fCurrent) return 0;
            fCurrent->fLastIndexQueried = -1;
            fCurrent->fLastIndexReturned = 0;
            fCurrent = (TEntryList*)fLists->After(fCurrent);
            // fCurrent is guarantee to be non-zero because it is not the 'last'
            // element of the list.
            if (!fCurrent) return 0;
            if (!fShift)
               result = fCurrent->Next();
            else {
               if (fCurrent->GetTreeNumber() >= 0)
                  result = fCurrent->Next();
            }
         }
         fLastIndexQueried++;
         fLastIndexReturned = result;
         return result;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Checks if the array representation is more economical and if so, switches to it

void TEntryList::OptimizeStorage()
{
   if (fBlocks){
      TEntryListBlock *block = 0;
      for (Int_t i=0; i<fNBlocks; i++){
         block = (TEntryListBlock*)fBlocks->UncheckedAt(i);
         block->OptimizeStorage();
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Print this list
/// - option = "" - default - print the name of the tree and file
/// - option = "all" - print all the entry numbers

void TEntryList::Print(const Option_t* option) const
{
   TString opt = option;
   opt.ToUpper();
   if (fBlocks) {
      Printf("%s %s %lld", fTreeName.Data(), fFileName.Data(), fN);
      if (opt.Contains("A")){
         TEntryListBlock* block = 0;
         for (Int_t i=0; i<fNBlocks; i++){
            block = (TEntryListBlock*)fBlocks->UncheckedAt(i);
            Int_t shift = i*kBlockSize;
            block->PrintWithShift(shift);
         }
      }
   }
   else {
      TEntryList *elist = 0;
      if (fN>0){
         TIter next(fLists);
         while((elist = (TEntryList*)next())){
            elist->Print(option);
         }
      } else {
         if (!fLists) Printf("%s %s %lld", fTreeName.Data(), fFileName.Data(), fN);
         else {
            TIter next(fLists);
            while ((elist = (TEntryList*)next())){
               Printf("%s %s %lld", elist->GetTreeName(), elist->GetFileName(), elist->GetN());
            }
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Reset this list

void TEntryList::Reset()
{
   //Maybe not delete, but just reset the number of blocks to 0????

   if (fBlocks){
      fBlocks->Delete();
      delete fBlocks;
      fBlocks = 0;
   }
   if (fLists){
      if (!((TEntryList*)fLists->First())->GetDirectory()){
         fLists->Delete();
      }
      delete fLists;
      fLists = 0;
   }
   fCurrent = 0;
   fBlocks = 0;
   fNBlocks = 0;
   fN = 0;
   fTreeName = "";
   fFileName = "";
   fStringHash = 0;
   fTreeNumber = -1;
   fLastIndexQueried = -1;
   fLastIndexReturned = 0;
   fReapply = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Add reference to directory dir. dir can be 0.

void TEntryList::SetDirectory(TDirectory *dir)
{
   if (fDirectory == dir) return;
   if (fDirectory) fDirectory->Remove(this);
   fDirectory = dir;
   if (fDirectory) fDirectory->Append(this);
}

////////////////////////////////////////////////////////////////////////////////
/// If a list for a tree with such name and filename exists, sets it as the current sublist
/// If not, creates this list and sets it as the current sublist
///
/// ! the filename is taken as provided, no extensions to full path or url !

void TEntryList::SetTree(const char *treename, const char *filename)
{
   TEntryList *elist = 0;

   TString fn;
   GetFileName(filename, fn);

   TString stotal = treename;
   stotal.Append(fn.Data());
   //printf("setting tree %s\n", stotal.Data());
   ULong_t newhash = stotal.Hash();
   if (fLists) {
      //find the corresponding entry list and make it current
      if (!fCurrent) fCurrent = (TEntryList*)fLists->First();
      if (fCurrent->fStringHash == 0){
         stotal = fCurrent->fTreeName + fCurrent->fFileName;
         fCurrent->fStringHash = stotal.Hash();
      }
      if (newhash == fCurrent->fStringHash){
         //this list is current
         if (!strcmp(fCurrent->fTreeName, treename) && !strcmp(fCurrent->fFileName, fn.Data())){
            return;
         }
      }
      TIter next(fLists);
      while ((elist = (TEntryList*)next())){
         if (newhash == elist->fStringHash){
            if (elist->fTreeName == treename && elist->fFileName == fn.Data()) {
               //the current entry list was changed. reset the fLastIndexQueried,
               //so that Next() doesn't start with the wrong current list
               //Also, reset those indices in the previously current list
               if (fCurrent->fBlocks){
                  Int_t currentblock = (fCurrent->fLastIndexReturned)/kBlockSize;
                  TEntryListBlock *block = (TEntryListBlock*)fCurrent->fBlocks->UncheckedAt(currentblock);
                  block->ResetIndices();
                  fCurrent->fLastIndexReturned = 0;
                  fCurrent->fLastIndexQueried = -1;
               }
               fCurrent = elist;
               fLastIndexQueried = -3;
               return;
            }
         }
      }
      //didn't find an entry list for this tree, create a new one
      elist = new TEntryList("", "", treename, fn.Data());
      if (elist->GetDirectory()) {
         //sub lists are not added to the current directory
         elist->GetDirectory()->Remove(elist);
         elist->SetDirectory(0);
      }
      fLists->Add(elist);
      fCurrent = elist;
      return;
   } else {
      if (fN==0 && fTreeName=="" && fFileName==""){
         //this is the first tree set to this list
         fTreeName = treename;
         fFileName = fn;
         stotal = fTreeName + fFileName;
         //fStringHash = stotal.Hash();
         fStringHash = newhash;
         fCurrent = this;
      } else {
         if (fStringHash == 0){
            stotal = fTreeName + fFileName;
            fStringHash = stotal.Hash();
         }
         if (newhash != fStringHash){
            //we have a chain and already have an entry list for the first tree
            //move the first entry list to the fLists
            fLists = new TList();
            elist = new TEntryList();
            elist->fTreeName = fTreeName;
            elist->fFileName = fFileName;
            elist->fStringHash = fStringHash;
            elist->fN = fN;
            elist->fTreeNumber = fTreeNumber;
            elist->fBlocks = fBlocks;
            fBlocks = 0;
            elist->fNBlocks = fNBlocks;
            fLists->Add(elist);
            elist = new TEntryList("", "", treename, fn.Data());
            if (elist->GetDirectory()) {
               //sub lists are not added to the current directory
               elist->GetDirectory()->Remove(elist);
               elist->SetDirectory(0);
            }
            fLists->Add(elist);
            fCurrent = elist;
            //the current entry list was changed. reset the fLastIndexQueried,
            //so that Next() doesn't start with the wrong current list
            fLastIndexQueried = -3;

         }
         else {
            //same tree as in the current entry list, don't do anything
            return;
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// If a list for a tree with such name and filename exists, sets it as the current sublist
/// If not, creates this list and sets it as the current sublist
/// The name of the file, where the tree is, is taken as
/// `tree->GetTree()->GetCurrentFile()->GetName()`, and then expanded either to the absolute path,
/// or to full url. If, for some reason, you want to provide
/// the filename in a different format, use SetTree(const char *treename, const char *filename),
/// where the filename is taken "as is".

void TEntryList::SetTree(const TTree *tree)
{
   if (!tree) return;
   auto thisTree = tree->GetTree();
   if (!thisTree) return;

   TString treename;
   if (tree->GetDirectory()->InheritsFrom("TFile")) {
      treename = thisTree->GetName();
   } else {
      treename = TString::Format("%s/%s",tree->GetDirectory()->GetName(),thisTree->GetName());
   }

   TString filename;
   if (tree->GetTree()->GetCurrentFile()){
      filename = tree->GetTree()->GetCurrentFile()->GetName();
      TUrl url(filename.Data(), kTRUE);
      if (!strcmp(url.GetProtocol(), "file")){
         gSystem->ExpandPathName(filename);
         if (!gSystem->IsAbsoluteFileName(filename))
            gSystem->PrependPathName(gSystem->pwd(), filename);
         filename = gSystem->UnixPathName(filename);
         url.SetFile(filename);
      }
      filename = url.GetUrl();
   } else {
      //memory-resident
      filename = "";
   }
   SetTree(treename, filename);

}

////////////////////////////////////////////////////////////////////////////////
/// Remove all the entries of this entry list, that are contained in elist

void TEntryList::Subtract(const TEntryList *elist)
{
   TEntryList *templist = 0;
   if (!fLists){
      if (!fBlocks) return;
      //check if lists are for the same tree
      if (!elist->fLists){
         //second list is also only for 1 tree
         if (!strcmp(elist->fTreeName.Data(),fTreeName.Data()) &&
             !strcmp(elist->fFileName.Data(),fFileName.Data())){
            //same tree
            Long64_t n2 = elist->GetN();
            Long64_t entry;
            for (Int_t i=0; i<n2; i++){
               entry = (const_cast<TEntryList*>(elist))->GetEntry(i);
               Remove(entry);
            }
         } else {
            //different trees
            return;
         }
      } else {
         //second list has sublists, try to find one for the same tree as this list
         TIter next1(elist->GetLists());
         templist = 0;
         Bool_t found = kFALSE;
         while ((templist = (TEntryList*)next1())){
            if (!strcmp(templist->fTreeName.Data(),fTreeName.Data()) &&
                !strcmp(templist->fFileName.Data(),fFileName.Data())){
               found = kTRUE;
               break;
            }
         }
         if (found) {
            Subtract(templist);
         }
      }
   } else {
      //this list has sublists
      TIter next2(fLists);
      templist = 0;
      Long64_t oldn=0;
      while ((templist = (TEntryList*)next2())){
         oldn = templist->GetN();
         templist->Subtract(elist);
         fN = fN - oldn + templist->GetN();
      }
   }
   return;
}

////////////////////////////////////////////////////////////////////////////////

TEntryList operator||(TEntryList &elist1, TEntryList &elist2)
{
   TEntryList eresult = elist1;
   //eresult = elist1;
   // printf("internal in operator1\n");
   eresult.Print("all");
   eresult.Add(&elist2);
   // printf("internal in operator2\n");
   eresult.Print("all");

   return eresult;
}

////////////////////////////////////////////////////////////////////////////////
/// Relocate the file paths.
/// If `oldroot` is defined, replace `oldroot` with `newroot` in all file names,
/// i.e. `oldroot/re/st/of/the/path` will become `newroot`/re/st/of/the/path`.
/// If `oldroot` is null, the new path will be just `newroot/path`.
/// Relocation is mandatory to use the entry-list with the same dataset at a different
/// location (i.e. on a different cluster, machine or disks).

Int_t TEntryList::RelocatePaths(const char *newroot, const char *oldroot)
{
   // At least newroot must be given
   if (!newroot || (newroot && strlen(newroot) <= 0)) {
      Warning("RelocatePaths", "the new location must be given!");
      return -1;
   }

   if (strlen(GetName()) > 0)
      Info("RelocatePaths", "'%s': relocating paths '%s' to '%s'",
                            GetName(), oldroot ? oldroot : "*", newroot);

   Int_t nrl = 0, xnrl = 0;
   // Apply to all underlying lists, if any
   if (fLists) {
      TIter nxl(fLists);
      TEntryList *enl = 0;
      while ((enl = (TEntryList *) nxl())) {
         if ((xnrl = enl->RelocatePaths(newroot, oldroot)) < 0) {
            Warning("RelocatePaths", "problems relocating '%s'", enl->GetName());
         } else {
            nrl += xnrl;
         }
      }
   }
   // Apply to ourselves
   TString temp;
   Ssiz_t lo = 0;
   if (oldroot && (lo = strlen(oldroot)) > 0) {
      if (fFileName.BeginsWith(oldroot)) {
         fFileName.Replace(0, lo, newroot);
         nrl++;
      }
   } else {
      Ssiz_t ilst = fFileName.Last('/');
      if (ilst != kNPOS) {
         fFileName.Replace(0, ilst, newroot);
      } else {
         fFileName.Insert(0, TString::Format("%s/", newroot));
      }
      nrl++;
   }
   if (fStringHash != 0) {
      temp.Form("%s%s", fTreeName.Data(), fFileName.Data());
      fStringHash = temp.Hash();
   }

   // Done
   return nrl;
}

////////////////////////////////////////////////////////////////////////////////
/// Relocate entry list 'enlnm' in file 'fn' replacing 'oldroot' with 'newroot' in
/// filenames. If 'enlnm' is null or '*' all entry lists in the file are relocated.
/// Relocation is mandatory to use the entry-list with the same dataset at a different
/// location (i.e. on a different cluster, machine or disks).
/// This function can be called as many times as need to reach the desired result.
/// The existing 'locations' can be checked qith TEntryList::Scan .

Int_t TEntryList::Relocate(const char *fn,
                          const char *newroot, const char *oldroot, const char *enlnm)
{
   // Open the file for updating
   TFile *fl = TFile::Open(fn, "UPDATE");
   if (!fl || (fl&& fl->IsZombie())) {
      ::Error("TEntryList::Relocate", "file '%s' cannot be open for updating", fn);
      return -1;
   }

   Int_t nrl = 0;
   // Read the lists
   TString nm(enlnm);
   if (nm.IsNull()) nm = "*";
   TRegexp nmrg(nm, kTRUE);
   TIter nxk(fl->GetListOfKeys());
   TKey *key = 0;
   while ((key = (TKey *) nxk())) {
      if (!strcmp(key->GetClassName(), "TEntryList")) {
         TString knm(key->GetName());
         if (knm.Index(nmrg) != kNPOS) {
            TEntryList *enl = dynamic_cast<TEntryList *>(fl->Get(knm));
            if (enl) {
               Int_t xnrl = enl->RelocatePaths(newroot, oldroot);
               if (xnrl >= 0) {
                  enl->Write(knm, TObject::kOverwrite);
                  nrl += xnrl;
               } else {
                  ::Error("TEntryList::Relocate", "problems relocating '%s' ...", enl->GetName());
               }
            }
         }
      }
   }
   // Close the file
   fl->Close();
   delete fl;
   // Done
   return nrl;
}

////////////////////////////////////////////////////////////////////////////////
/// Get in 'c' the string in common at the beginning of 'a' and 'b'
///
/// Return:
/// - 0         a and b are not contained in each other, i.e. c != a && c != b
/// - 1         a is contained in b, i.e. c == a (includes a == empty)
/// - 2         b is contained in a, i.e. c == b (includes b == empty)
/// - 3         b is a, i.e. c == b == a (includes a == b == empty)
/// Auxiliary function for path scans.

static Int_t GetCommonString(TString a, TString b, TString &c)
{
   if (a == b) {
      c = a;
      return 3;
   }
   if (a.IsNull()) {
      c = "";
      return 1;
   }
   if (b.IsNull()) {
      c = "";
      return 2;
   }
   Bool_t ashort = (a.Length() > b.Length()) ? kFALSE : kTRUE;
   Ssiz_t len = (ashort) ? a.Length() : b.Length();
   Int_t lcom = 0;
   for (Int_t i = 0; i < len; i++) {
      if (a[i] != b[i]) break;
      lcom++;
   }
   if (lcom == len) {
      c = ashort ? a : b;
      return ashort ? 1 : 2;
   }
   c = a(0,lcom);
   // Done
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Scan the paths to find the common roots. If 'roots' is defined, add
/// the found roots to the list as TObjStrings.
/// Return the number of roots found.

Int_t TEntryList::ScanPaths(TList *roots, Bool_t notify)
{
   TList *xrl = roots ? roots : new TList;

   Int_t nrl = 0;
   // Apply to all underlying lists, if any
   if (fLists) {
      TIter nxl(fLists);
      TEntryList *enl = 0;
      while ((enl = (TEntryList *) nxl()))
         nrl += enl->ScanPaths(xrl, kFALSE);
   }
   // Apply to ourselves
   Bool_t newobjs = kTRUE;
   TString path = gSystem->GetDirName(fFileName), com;
   TObjString *objs = nullptr;
   TIter nxr(xrl);
   while ((objs = (TObjString *) nxr())) {
      Int_t rc = 0;
      if ((rc = GetCommonString(path, objs->GetString(), com)) != 2) {
         TUrl ucom(com);
         if (strlen(ucom.GetFile()) > 0 && strcmp(ucom.GetFile(), "/")) {
            objs->SetString(com.Data());
            newobjs = kFALSE;
            break;
         }
      }
   }
   if (newobjs) xrl->Add(new TObjString(path));

   // Done
   nrl = xrl->GetSize();
   if (notify) {
      Printf(" * ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ *");
      Printf(" * Entry-list: %s", GetName());
      Printf(" *   %d common root paths found", nrl);
      nxr.Reset();
      while ((objs = (TObjString *) nxr())) {
         Printf(" *     %s", objs->GetName());
      }
      Printf(" * ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ *");
   }

   if (xrl != roots) {
      xrl->SetOwner(kTRUE);
      SafeDelete(xrl);
   }

   // Done
   return nrl;
}

////////////////////////////////////////////////////////////////////////////////
/// Scan TEntryList in 'fn' to find the common parts of paths.
/// If 'roots' is defined, add the found roots to the list as TObjStrings.
/// Return the number of common root paths found.

Int_t TEntryList::Scan(const char *fn, TList *roots)
{
   // Open the file for updating
   TFile *fl = TFile::Open(fn);
   if (!fl || (fl&& fl->IsZombie())) {
      ::Error("TEntryList::Relocate", "file '%s' cannot be open for reading", fn);
      return -1;
   }

   Int_t nrs = 0;
   // Read the lists
   TIter nxk(fl->GetListOfKeys());
   TKey *key = 0;
   while ((key = (TKey *) nxk())) {
      if (!strcmp(key->GetClassName(), "TEntryList")) {
         TEntryList *enl = dynamic_cast<TEntryList *>(fl->Get(key->GetName()));
         if (enl) {
            nrs += enl->ScanPaths(roots);
         } else {
            ::Error("TEntryList::Scan", "object entry-list '%s' not found or not loadable!", key->GetName());
         }
      }
   }
   // Close the file
   fl->Close();
   delete fl;

   // Done
   return nrs;
}

////////////////////////////////////////////////////////////////////////////////
/// Custom streamer for class TEntryList to handle the different interpretation
/// of fFileName between version 1 and >1 .

void TEntryList::Streamer(TBuffer &b)
{
   if (b.IsReading()) {
      UInt_t R__s, R__c;
      Version_t R__v = b.ReadVersion(&R__s, &R__c);
      b.ReadClassBuffer(TEntryList::Class(), this, R__v, R__s, R__c);
      if (R__v <= 1) {
         // The filename contained also the protocol and host: this was dropped
         // in version > 1 to allow re-localization
         GetFileName(fFileName.Data(), fFileName);
      }
   } else {
      b.WriteClassBuffer(TEntryList::Class(), this);
   }
}
