// @(#)root/cont:$Id$
// Author: Fons Rademakers   27/09/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_THashTable
#define ROOT_THashTable


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// THashTable                                                           //
//                                                                      //
// THashTable implements a hash table to store TObject's. The hash      //
// value is calculated using the value returned by the TObject's        //
// Hash() function. Each class inheriting from TObject can override     //
// Hash() as it sees fit.                                               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TCollection.h"
#include "TString.h"

class TList;
class TListIter;
class THashTableIter;


class THashTable : public TCollection {

friend class  THashTableIter;

private:
   TList     **fCont;          //Hash table (table of lists)
   Int_t       fEntries;       //Number of objects in table
   Int_t       fUsedSlots;     //Number of used slots
   Int_t       fRehashLevel;   //Average collision rate which triggers rehash

   Int_t       GetCheckedHashValue(TObject *obj) const;
   Int_t       GetHashValue(const TObject *obj) const;
   Int_t       GetHashValue(TString &s) const { return s.Hash() % fSize; }
   Int_t       GetHashValue(const char *str) const { return ::Hash(str) % fSize; }

   void        AddImpl(Int_t slot, TObject *object);

   THashTable(const THashTable&) = delete;
   THashTable& operator=(const THashTable&) = delete;

public:
   THashTable(Int_t capacity = TCollection::kInitHashTableCapacity, Int_t rehash = 0);
   virtual       ~THashTable();
   void          Add(TObject *obj);
   void          AddBefore(const TObject *before, TObject *obj);
   virtual void  AddAll(const TCollection *col);
   Float_t       AverageCollisions() const;
   void          Clear(Option_t *option="");
   Int_t         Collisions(const char *name) const;
   Int_t         Collisions(TObject *obj) const;
   void          Delete(Option_t *option="");
   Bool_t        Empty() const { return fEntries == 0; }
   TObject      *FindObject(const char *name) const;
   TObject      *FindObject(const TObject *obj) const;
   const TList  *GetListForObject(const char *name) const;
   const TList  *GetListForObject(const TObject *obj) const;
   TObject     **GetObjectRef(const TObject *obj) const;
   Int_t         GetRehashLevel() const { return fRehashLevel; }
   Int_t         GetSize() const { return fEntries; }
   TIterator    *MakeIterator(Bool_t dir = kIterForward) const;
   void          Print(Option_t *option, Int_t recurse) const;
   using TCollection::Print;
   void          Rehash(Int_t newCapacity, Bool_t checkObjValidity = kTRUE);
   TObject      *Remove(TObject *obj);
   TObject      *RemoveSlow(TObject *obj);
   void          SetRehashLevel(Int_t rehash) { fRehashLevel = rehash; }

   ClassDef(THashTable,0)  //A hash table
};

inline Float_t THashTable::AverageCollisions() const
{
   if (fUsedSlots)
      return ((Float_t)fEntries)/((Float_t)fUsedSlots);
   else
      return 0.0;
}

inline Int_t THashTable::GetCheckedHashValue(TObject *obj) const
{
   Int_t i = Int_t(obj->CheckedHash() % fSize); // need intermediary i for Linux g++
   return i;
}

inline Int_t THashTable::GetHashValue(const TObject *obj) const
{
   Int_t i = Int_t(obj->Hash() % fSize);  // need intermediary i for Linux g++
   return i;
}


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// THashTableIter                                                       //
//                                                                      //
// Iterator of hash table.                                              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

class THashTableIter : public TIterator {

private:
   const THashTable *fTable;       //hash table being iterated
   Int_t             fCursor;      //current position in table
   TListIter        *fListCursor;  //current position in collision list
   Bool_t            fDirection;   //iteration direction

   THashTableIter() : fTable(nullptr), fCursor(0), fListCursor(nullptr), fDirection(kIterForward) { }
   Int_t             NextSlot();

public:
   THashTableIter(const THashTable *ht, Bool_t dir = kIterForward);
   THashTableIter(const THashTableIter &iter);
   ~THashTableIter();
   TIterator      &operator=(const TIterator &rhs);
   THashTableIter &operator=(const THashTableIter &rhs);

   const TCollection *GetCollection() const { return fTable; }
   TObject           *Next();
   void               Reset();
   Bool_t             operator!=(const TIterator &aIter) const;
   Bool_t             operator!=(const THashTableIter &aIter) const;
   TObject           *operator*() const;

   ClassDef(THashTableIter,0)  //Hash table iterator
};

#endif
