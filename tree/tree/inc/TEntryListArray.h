// @(#)root/tree:$Id$
// Author: Bruno Lenzi 12/07/2011

/*************************************************************************
 * Copyright (C) 1995-2006, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TEntryListArray
#define ROOT_TEntryListArray

#include "TEntryList.h"

class TTree;
class TDirectory;
class TObjArray;
class TString;

class TList;
class TCollection;
class TIter;

class TEntryListArray : public TEntryList {

private:
   TEntryListArray& operator=(const TEntryListArray&); // Not implemented

protected:
   TList *fSubLists;                     ///<  a list of underlying entry lists for each event of a TEntryList
   Long64_t         fEntry;              ///<  the entry number, when the list is used for subentries
   TEntryListArray *fLastSubListQueried; ///<! last sublist checked by GetSubListForEntry
   TIter *fSubListIter;                  ///<! to iterate over fSubLists and keep last one checked

   void Init();
   virtual void AddEntriesAndSubLists(const TEntryList *elist);
   virtual void ConvertToTEntryListArray(TEntryList *e);
//    virtual TList* GetSubLists() const {
//       return fSubLists;
//    };
   virtual Bool_t      RemoveSubList(TEntryListArray *e, TTree *tree = 0);
   virtual Bool_t      RemoveSubListForEntry(Long64_t entry, TTree *tree = 0);
   virtual TEntryListArray* SetEntry(Long64_t entry, TTree *tree = 0);


public:
   TEntryListArray();
   TEntryListArray(const char *name, const char *title);
   TEntryListArray(const char *name, const char *title, const TTree *tree);
   TEntryListArray(const char *name, const char *title, const char *treename, const char *filename);
   TEntryListArray(const TTree *tree);
   TEntryListArray(const TEntryListArray& elist);
   TEntryListArray(const TEntryList& elist); // to convert TEL to TELA
   virtual ~TEntryListArray();

   virtual void        Add(const TEntryList *elist);
   virtual Int_t       Contains(Long64_t entry, TTree *tree, Long64_t subentry);
   virtual Int_t       Contains(Long64_t entry, TTree *tree = 0) {
      return TEntryList::Contains(entry, tree);
   };
   virtual Bool_t      Enter(Long64_t entry, TTree *tree, Long64_t subentry);
   virtual Bool_t      Enter(Long64_t entry, TTree *tree = 0) {
      return Enter(entry, tree, -1);
   };
//    virtual Bool_t      Enter(Long64_t entry, TTree *tree, const TEntryList *e);
   virtual TEntryListArray* GetSubListForEntry(Long64_t entry, TTree *tree = 0);
   virtual void        Print(const Option_t* option = "") const;
   virtual Bool_t      Remove(Long64_t entry, TTree *tree, Long64_t subentry);
   virtual Bool_t      Remove(Long64_t entry, TTree *tree = 0) {
      return Remove(entry, tree, -1);
   };
   virtual void        Reset();

   virtual void        SetTree(const char *treename, const char *filename);
   virtual void        SetTree(const TTree *tree) {
      TEntryList::SetTree(tree);   // will take treename and filename from the tree and call the method above
   }
   virtual void        Subtract(const TEntryList *elist);
   virtual TList* GetSubLists() const {
      return fSubLists;
   };

   ClassDef(TEntryListArray, 1);  //A list of entries and subentries in a TTree
};
#endif
