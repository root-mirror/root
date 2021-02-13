// @(#)root/eve:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TEveChunkManager
#define ROOT_TEveChunkManager

#include "TEveUtil.h"

#include "TArrayC.h"

#include <vector>
#include <set>

/******************************************************************************/
// TEveChunkManager
/******************************************************************************/

class TEveChunkManager
{
private:
   TEveChunkManager(const TEveChunkManager&) = delete;
   TEveChunkManager& operator=(const TEveChunkManager&) = delete;

protected:
   Int_t fS;        // Size of atom
   Int_t fN;        // Number of atoms in a chunk

   Int_t fSize;     // Size of container, number of atoms
   Int_t fVecSize;  // Number of allocated chunks
   Int_t fCapacity; // Available capacity within the chunks

   std::vector<TArrayC*> fChunks; // Memory blocks

   void ReleaseChunks();

public:
   TEveChunkManager();
   TEveChunkManager(Int_t atom_size, Int_t chunk_size);
   virtual ~TEveChunkManager();

   void Reset(Int_t atom_size, Int_t chunk_size);
   void Refit();

   Int_t    S() const { return fS; }
   Int_t    N() const { return fN; }

   Int_t    Size()     const { return fSize; }
   Int_t    VecSize()  const { return fVecSize; }
   Int_t    Capacity() const { return fCapacity; }

   Char_t* Atom(Int_t idx)   const { return fChunks[idx/fN]->fArray + idx%fN*fS; }
   Char_t* Chunk(Int_t chk)  const { return fChunks[chk]->fArray; }
   Int_t   NAtoms(Int_t chk) const { return (chk < fVecSize-1) ? fN : (fSize-1)%fN + 1; }

   Char_t* NewAtom();
   Char_t* NewChunk();


   // Iterator

   struct iterator
   {
      TEveChunkManager *fPlex;
      Char_t           *fCurrent;
      Int_t             fAtomIndex;
      Int_t             fNextChunk;
      Int_t             fAtomsToGo;

      const std::set<Int_t>           *fSelection;
      std::set<Int_t>::const_iterator  fSelectionIterator;

      iterator(TEveChunkManager* p) :
         fPlex(p), fCurrent(0), fAtomIndex(-1),
         fNextChunk(0), fAtomsToGo(0), fSelection(0), fSelectionIterator() {}
      iterator(TEveChunkManager& p) :
         fPlex(&p), fCurrent(0), fAtomIndex(-1),
         fNextChunk(0), fAtomsToGo(0), fSelection(0), fSelectionIterator() {}
      iterator(const iterator& i) :
         fPlex(i.fPlex), fCurrent(i.fCurrent), fAtomIndex(i.fAtomIndex),
         fNextChunk(i.fNextChunk), fAtomsToGo(i.fAtomsToGo),
         fSelection(i.fSelection), fSelectionIterator(i.fSelectionIterator) {}

      iterator& operator=(const iterator& i) {
         fPlex = i.fPlex; fCurrent = i.fCurrent; fAtomIndex = i.fAtomIndex;
         fNextChunk = i.fNextChunk; fAtomsToGo = i.fAtomsToGo;
         fSelection = i.fSelection; fSelectionIterator = i.fSelectionIterator;
         return *this;
      }

      Bool_t  next();
      void    reset() { fCurrent = 0; fAtomIndex = -1; fNextChunk = fAtomsToGo = 0; }

      Char_t* operator()() { return fCurrent; }
      Char_t* operator*()  { return fCurrent; }
      Int_t   index()      { return fAtomIndex; }
   };

   ClassDef(TEveChunkManager, 1); // Vector-like container with chunked memory allocation.
};


//______________________________________________________________________________
inline Char_t* TEveChunkManager::NewAtom()
{
   Char_t *a = (fSize >= fCapacity) ? NewChunk() : Atom(fSize);
   ++fSize;
   return a;
}


/******************************************************************************/
// Templated some-class TEveChunkVector
/******************************************************************************/

template<class T>
class TEveChunkVector : public TEveChunkManager
{
private:
   TEveChunkVector(const TEveChunkVector&);            // Not implemented
   TEveChunkVector& operator=(const TEveChunkVector&); // Not implemented

public:
   TEveChunkVector()                 : TEveChunkManager() {}
   TEveChunkVector(Int_t chunk_size) : TEveChunkManager(sizeof(T), chunk_size) {}
   virtual ~TEveChunkVector() {}

   void Reset(Int_t chunk_size) { Reset(sizeof(T), chunk_size); }

   T* At(Int_t idx)  { return reinterpret_cast<T*>(Atom(idx)); }
   T& Ref(Int_t idx) { return *At(idx); }

   ClassDef(TEveChunkVector, 1); // Templated class for specific atom classes (given as template argument).
};

#endif
