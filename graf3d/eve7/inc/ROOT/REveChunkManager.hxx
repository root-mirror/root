// @(#)root/eve7:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2018

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_REveChunkManager
#define ROOT_REveChunkManager

#include <ROOT/REveUtil.hxx>

#include "TArrayC.h"

#include <vector>

namespace ROOT {
namespace Experimental {

/******************************************************************************/
// REveChunkManager
/******************************************************************************/

class REveChunkManager
{
private:
   REveChunkManager(const REveChunkManager&) = delete;
   REveChunkManager& operator=(const REveChunkManager&) = delete;

protected:
   Int_t fS;        // Size of atom
   Int_t fN;        // Number of atoms in a chunk

   Int_t fSize;     // Size of container, number of atoms
   Int_t fVecSize;  // Number of allocated chunks
   Int_t fCapacity; // Available capacity within the chunks

   std::vector<TArrayC*> fChunks; // Memory blocks

   void ReleaseChunks();

public:
   REveChunkManager();
   REveChunkManager(Int_t atom_size, Int_t chunk_size);
   virtual ~REveChunkManager();

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
      REveChunkManager &fPlex;
      Char_t           *fCurrent{nullptr};
      Int_t             fAtomIndex{-1};
      Int_t             fNextChunk{0};
      Int_t             fAtomsToGo{0};

      const std::set<Int_t>           *fSelection{nullptr};
      std::set<Int_t>::const_iterator  fSelectionIterator;

      iterator(REveChunkManager &p) : fPlex(p) {}

/*      iterator(const iterator &i) :
         fPlex(i.fPlex), fCurrent(i.fCurrent), fAtomIndex(i.fAtomIndex),
         fNextChunk(i.fNextChunk), fAtomsToGo(i.fAtomsToGo),
         fSelection(i.fSelection), fSelectionIterator(i.fSelectionIterator) {}

      iterator& operator=(const iterator& i) {
         fPlex = i.fPlex; fCurrent = i.fCurrent; fAtomIndex = i.fAtomIndex;
         fNextChunk = i.fNextChunk; fAtomsToGo = i.fAtomsToGo;
         fSelection = i.fSelection; fSelectionIterator = i.fSelectionIterator;
         return *this;
      }
*/

      Bool_t  next();
      void    reset() { fCurrent = nullptr; fAtomIndex = -1; fNextChunk = fAtomsToGo = 0; }

      Char_t* operator()() { return fCurrent; }
      Char_t* operator*()  { return fCurrent; }
      Int_t   index()      { return fAtomIndex; }
   };

};


//______________________________________________________________________________
inline Char_t* REveChunkManager::NewAtom()
{
   Char_t *a = (fSize >= fCapacity) ? NewChunk() : Atom(fSize);
   ++fSize;
   return a;
}


/******************************************************************************/
// Templated some-class REveChunkVector
/******************************************************************************/

template<class T>
class REveChunkVector : public REveChunkManager
{
private:
   REveChunkVector(const REveChunkVector&) = delete;
   REveChunkVector& operator=(const REveChunkVector&) = delete;

public:
   REveChunkVector() = default;
   REveChunkVector(Int_t chunk_size) : REveChunkManager(sizeof(T), chunk_size) {}
   virtual ~REveChunkVector() {}

   void Reset(Int_t chunk_size) { Reset(sizeof(T), chunk_size); }

   T *At(Int_t idx)  { return reinterpret_cast<T*>(Atom(idx)); }
   T &Ref(Int_t idx) { return *At(idx); }
};

} // namespace Experimental
} // namespace ROOT

#endif
