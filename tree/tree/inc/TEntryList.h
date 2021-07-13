// @(#)root/tree:$Id$
// Author: Anna Kreshuk 27/10/2006

/*************************************************************************
 * Copyright (C) 1995-2006, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TEntryList
#define ROOT_TEntryList

#include "TNamed.h"

class TTree;
class TDirectory;
class TObjArray;
class TString;

class TList;
class TCollection;

class TEntryList: public TNamed
{
 private:
   TEntryList& operator=(const TEntryList&); // Not implemented

 protected:
   TList      *fLists;                  ///<  a list of underlying entry lists for each tree of a chain
   TEntryList *fCurrent;                ///<! currently filled entry list

   Int_t            fNBlocks;           ///<  number of TEntryListBlocks
   TObjArray       *fBlocks;            ///<  blocks with indices of passing events (TEntryListBlocks)
   Long64_t         fN;                 ///<  number of entries in the list
   Long64_t         fEntriesToProcess;  ///<  used on proof to set the number of entries to process in a packet
   TString          fTreeName;          ///<  name of the tree
   TString          fFileName;          ///<  name of the file, where the tree is
   ULong_t          fStringHash;        ///<! Hash value of a string of treename and filename
   Int_t            fTreeNumber;        ///<! the index of the tree in the chain (used when the entry
                                        ///<  list is used as input (TTree::SetEntryList())

   Long64_t         fLastIndexQueried;  ///<! used to optimize GetEntry() function from a loop
   Long64_t         fLastIndexReturned; ///<! used to optimize GetEntry() function from a loop
   Bool_t           fShift;             ///<! true when some sub-lists don't correspond to trees
                                        ///<  (when the entry list is used as input in TChain)
   TDirectory      *fDirectory;         ///<! Pointer to directory holding this tree
   Bool_t           fReapply;           ///<  If true, TTree::Draw will 'reapply' the original cut

   void             GetFileName(const char *filename, TString &fn, Bool_t * = 0);

 public:
   enum {kBlockSize = 64000}; //number of entries in each block (not the physical size).

   TEntryList();
   TEntryList(const char *name, const char *title);
   TEntryList(const char *name, const char *title, const TTree *tree);
   TEntryList(const char *name, const char *title, const char *treename, const char *filename);
   TEntryList(const TTree *tree);
   TEntryList(const TEntryList& elist);
   virtual ~TEntryList();

   virtual void        Add(const TEntryList *elist);
   void                AddSubList(TEntryList *elist);
   virtual Int_t       Contains(Long64_t entry, TTree *tree = 0);
   virtual void        DirectoryAutoAdd(TDirectory *);
   virtual Bool_t      Enter(Long64_t entry, TTree *tree = 0);
   virtual TEntryList *GetCurrentList() const { return fCurrent; };
   virtual TEntryList *GetEntryList(const char *treename, const char *filename, Option_t *opt="");
   virtual Long64_t    GetEntry(Int_t index);
   virtual Long64_t    GetEntryAndTree(Int_t index, Int_t &treenum);
   virtual Long64_t    GetEntriesToProcess() const {return fEntriesToProcess;}
   virtual TList      *GetLists() const { return fLists; }
   virtual TDirectory *GetDirectory() const { return fDirectory; }
   virtual Long64_t    GetN() const { return fN; }
   virtual const char *GetTreeName() const { return fTreeName.Data(); }
   virtual const char *GetFileName() const { return fFileName.Data(); }
   virtual Int_t       GetTreeNumber() const { return fTreeNumber; }
   virtual Bool_t      GetReapplyCut() const { return fReapply; };

   Bool_t IsValid() const
   {
      if ((fLists || fBlocks)) return kTRUE;
      return kFALSE;
   }

   virtual Int_t       Merge(TCollection *list);

   virtual Long64_t    Next();
   virtual void        OptimizeStorage();
   virtual Int_t       RelocatePaths(const char *newloc, const char *oldloc = 0);
   virtual Bool_t      Remove(Long64_t entry, TTree *tree = 0);
   virtual void        Reset();
   virtual Int_t       ScanPaths(TList *roots, Bool_t notify = kTRUE);

   virtual void        Print(const Option_t* option = "") const;
   virtual void        SetDirectory(TDirectory *dir);
   virtual void        SetEntriesToProcess(Long64_t nen) { fEntriesToProcess = nen; }
   virtual void        SetShift(Bool_t shift) { fShift = shift; };
   virtual void        SetTree(const TTree *tree);
   virtual void        SetTree(const char *treename, const char *filename);
   virtual void        SetTreeName(const char *treename){ fTreeName = treename; };
   virtual void        SetFileName(const char *filename){ fFileName = filename; };
   virtual void        SetTreeNumber(Int_t index) { fTreeNumber=index;  }
   virtual void        SetReapplyCut(Bool_t apply = kFALSE) {fReapply = apply;}; // *TOGGLE* *GETTER=GetReapplyCut
   virtual void        Subtract(const TEntryList *elist);

   static  Int_t       Relocate(const char *fn,
                                const char *newroot, const char *oldroot = 0, const char *enlnm = 0);
   static  Int_t       Scan(const char *fn, TList *roots);

// Preventing warnings with -Weffc++ in GCC since the overloading of the || operator was a design choice.
#if (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__) >= 40600
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#endif
   friend TEntryList operator||(TEntryList& elist1, TEntryList& elist2);
#if (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__) >= 40600
#pragma GCC diagnostic pop
#endif

   ClassDef(TEntryList, 2);  //A list of entries in a TTree
};
#endif
