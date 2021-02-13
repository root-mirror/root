// @(#)root/cont:$Id$
// Author: Fons Rademakers   13/09/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TOrdCollection
#define ROOT_TOrdCollection


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TOrdCollection                                                       //
//                                                                      //
// Ordered collection.                                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TSeqCollection.h"

#include <iterator>


class TOrdCollectionIter;


class TOrdCollection : public TSeqCollection {

friend class  TOrdCollectionIter;

private:
   TObject  **fCont;
   Int_t      fCapacity;
   Int_t      fGapStart;
   Int_t      fGapSize;

   Int_t      PhysIndex(Int_t idx) const;
   Int_t      LogIndex(Int_t idx) const;
   void       MoveGapTo(Int_t newGapStart);
   Bool_t     IllegalIndex(const char *method, Int_t idx) const;
   void       Init(Int_t capacity);
   Bool_t     LowWaterMark() const;
   void       SetCapacity(Int_t newCapacity);

   TOrdCollection(const TOrdCollection&) = delete;
   TOrdCollection& operator=(const TOrdCollection&) = delete;

public:
   enum { kDefaultCapacity = 1, kMinExpand = 8, kShrinkFactor = 2 };

   typedef TOrdCollectionIter Iterator_t;

   TOrdCollection(Int_t capacity = kDefaultCapacity);
   ~TOrdCollection();
   void          Clear(Option_t *option="");
   void          Delete(Option_t *option="");
   TObject     **GetObjectRef(const TObject *obj) const;
   Int_t         IndexOf(const TObject *obj) const;
   TIterator    *MakeIterator(Bool_t dir = kIterForward) const;

   void          AddFirst(TObject *obj);
   void          AddLast(TObject *obj);
   void          AddAt(TObject *obj, Int_t idx);
   void          AddAfter(const TObject *after, TObject *obj);
   void          AddBefore(const TObject *before, TObject *obj);
   void          PutAt(TObject *obj, Int_t idx);
   TObject      *RemoveAt(Int_t idx);
   TObject      *Remove(TObject *obj);

   TObject      *At(Int_t idx) const;
   TObject      *Before(const TObject *obj) const;
   TObject      *After(const TObject *obj) const;
   TObject      *First() const;
   TObject      *Last() const;

   void          Sort();
   Int_t         BinarySearch(TObject *obj);

   ClassDef(TOrdCollection,0)  //An ordered collection
};


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TOrdCollectionIter                                                   //
//                                                                      //
// Iterator of ordered collection.                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

class TOrdCollectionIter : public TIterator,
                           public std::iterator<std::bidirectional_iterator_tag,
                                                TObject*, std::ptrdiff_t,
                                                const TObject**, const TObject*&> {

private:
   const TOrdCollection  *fCol;       //collection being iterated
   Int_t                  fCurCursor; //current position in collection
   Int_t                  fCursor;    //next position in collection
   Bool_t                 fDirection; //iteration direction

   TOrdCollectionIter() : fCol(nullptr), fCurCursor(0), fCursor(0), fDirection(kIterForward) { }

public:
   TOrdCollectionIter(const TOrdCollection *col, Bool_t dir = kIterForward);
   TOrdCollectionIter(const TOrdCollectionIter &iter);
   ~TOrdCollectionIter() { }
   TIterator          &operator=(const TIterator &rhs);
   TOrdCollectionIter &operator=(const TOrdCollectionIter &rhs);

   const TCollection *GetCollection() const { return fCol; }
   TObject           *Next();
   void               Reset();
   Bool_t             operator!=(const TIterator &aIter) const;
   Bool_t             operator!=(const TOrdCollectionIter &aIter) const;
   TObject           *operator*() const;

   ClassDef(TOrdCollectionIter,0)  //Ordered collection iterator
};


//---- inlines -----------------------------------------------------------------

inline Bool_t TOrdCollection::LowWaterMark() const
{
   return (fSize < (fCapacity / 4) && fSize > TCollection::kInitCapacity);
}

inline Int_t TOrdCollection::PhysIndex(Int_t idx) const
   { return (idx < fGapStart) ? idx : idx + fGapSize; }

inline Int_t TOrdCollection::LogIndex(Int_t idx) const
   { return (idx < fGapStart) ? idx : idx - fGapSize; }

#endif
