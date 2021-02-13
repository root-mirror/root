// @(#)root/roofit:$Id$
// Author: Stephan Hageboeck, CERN, 10/2018
/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** Memory pool for RooArgSet and RooDataSet.
 * \class MemPoolForRooSets
 * \ingroup roofitcore
 * RooArgSet and RooDataSet were using a mempool that guarantees that allocating,
 * de-allocating and re-allocating a set does not yield the same pointer. Since
 * both were using the same logic, the functionality has been put in this class.
 * This class solves RooFit's static destruction order problems by intentionally leaking
 * arenas of the mempool that still contain live objects at the end of the program.
 *
 * When the set types are compared based on a unique ID instead of their pointer,
 * one can go back to normal memory management, and this class becomes obsolete.
 */

#ifndef ROOFIT_ROOFITCORE_SRC_MEMPOOLFORROOSETS_H_
#define ROOFIT_ROOFITCORE_SRC_MEMPOOLFORROOSETS_H_

#include "TStorage.h"

#include <vector>
#include <algorithm>
#include <set>

template <class RooSet_t, std::size_t POOLSIZE>
class MemPoolForRooSets {

  struct Arena {
    Arena()
      : ownedMemory{static_cast<RooSet_t *>(TStorage::ObjectAlloc(POOLSIZE * sizeof(RooSet_t)))},
        memBegin{ownedMemory}, nextItem{ownedMemory},
        memEnd{memBegin + POOLSIZE}, refCount{0}
    {
    }



    Arena(const Arena &) = delete;
    Arena(Arena && other)
      : ownedMemory{other.ownedMemory},
        memBegin{other.memBegin}, nextItem{other.nextItem}, memEnd{other.memEnd},
        refCount{other.refCount}
#ifndef NDEBUG
      , deletedElements { std::move(other.deletedElements) }
#endif
    {
      // Needed for unique ownership
      other.ownedMemory = nullptr;
      other.refCount = 0;
    }



    Arena & operator=(const Arena &) = delete;
    Arena & operator=(Arena && other)
    {
      ownedMemory = other.ownedMemory;
      memBegin = other.memBegin;
      nextItem = other.nextItem;
      memEnd   = other.memEnd;
#ifndef NDEBUG
      deletedElements = std::move(other.deletedElements);
#endif
      refCount     = other.refCount;

      other.ownedMemory = nullptr;
      other.refCount = 0;

      return *this;
    }



    // If there is any user left, the arena shouldn't be deleted.
    // If this happens, nevertheless, one has an order of destruction problem.
    ~Arena()
    {
      if (!ownedMemory) return;

      if (refCount != 0) {
        std::cerr << __FILE__ << ":" << __LINE__ << "Deleting arena " << ownedMemory << " with use count " << refCount
                  << std::endl;
        assert(false);
      }

      ::operator delete(ownedMemory);
    }


    bool inPool(const RooSet_t * const ptr) const {
      return memBegin <= ptr && ptr < memEnd;
    }

    bool inPool(const void * const ptr) const
    {
      return inPool(static_cast<const RooSet_t * const>(ptr));
    }

    bool hasSpace() const { return ownedMemory && nextItem < memEnd; }
    bool empty() const { return refCount == 0; }

    void tryFree(bool freeNonFull) {
      if (ownedMemory && empty() && (!hasSpace() || freeNonFull) ) {
        ::operator delete(ownedMemory);
        ownedMemory = nullptr;
      }
    }

    void * tryAllocate()
    {
      if (!hasSpace()) return nullptr;

      ++refCount;
      return nextItem++;
    }

    bool tryDeallocate(void * ptr)
    {
      if (inPool(ptr)) {
        --refCount;
#ifndef NDEBUG
        const std::size_t index = static_cast<RooSet_t *>(ptr) - memBegin;
        if (deletedElements.count(index) != 0) {
          std::cerr << "Double delete of " << ptr << " at index " << index << " in Arena with refCount " << refCount
              << ".\n\tArena: |" << memBegin << "\t" << ptr << "\t" << memEnd << "|" << std::endl;
          throw;
        }
        deletedElements.insert(index);
#endif
        return true;
      } else
        return false;
    }

    bool memoryOverlaps(const Arena& other) const {
      //Need the reinterpret_cast to correctly check for non-overlap on the last byte of the last element
      return inPool(other.memBegin) || inPool(reinterpret_cast<const char*>(other.memEnd)-1);
    }

    RooSet_t * ownedMemory;
    const RooSet_t * memBegin;
    RooSet_t * nextItem;
    const RooSet_t * memEnd;
    std::size_t refCount;
#ifndef NDEBUG
    std::set<std::size_t> deletedElements;
#endif
  };


  public:
  /// Create empty mem pool.
  MemPoolForRooSets() : fArenas{} {}

  MemPoolForRooSets(const MemPoolForRooSets &) = delete;
  MemPoolForRooSets(MemPoolForRooSets &&)      = delete;
  MemPoolForRooSets & operator=(const MemPoolForRooSets &) = delete;
  MemPoolForRooSets & operator=(MemPoolForRooSets &&) = delete;

  /// Destructor. Should not be called when RooArgSets or RooDataSets are still alive.
  ~MemPoolForRooSets()
  {
    if (!empty()) {
#ifndef _MSC_VER
      std::cerr << __PRETTY_FUNCTION__;
#endif
      std::cerr << " The mem pool being deleted is not empty. This will lead to crashes."
                << std::endl;
      assert(false);
    }
  }



  /// Allocate memory for the templated set type. Fails if bytes != sizeof(RooSet_t).
  void * allocate(std::size_t bytes)
  {
    if (bytes != sizeof(RooSet_t))
      throw std::bad_alloc();

    if (fArenas.empty() || !fArenas.back().hasSpace()) {
      newArena();
      prune();
    }

    void * ptr = fArenas.back().tryAllocate();
    assert(ptr != nullptr);

    return ptr;
  }



  /// Deallocate memory for the templated set type if in pool.
  /// \return True if element was in pool.
  bool deallocate(void * ptr)
  {
    bool deallocSuccess = false;

    if (std::any_of(fArenas.begin(), fArenas.end(),
      [ptr](Arena& arena){return arena.tryDeallocate(ptr);})) {
        deallocSuccess = true;
    }

    if (fTeardownMode) {
      // Try pruning after each dealloc because we are tearing down
      prune();
    }

    return deallocSuccess;
  }



  ////////////////////////////////////////////////////////////////////////////////
  /// Free memory in arenas that don't have space and no users.
  /// In fTeardownMode, it will also delete the arena that still has space.
  ///
  void prune()
  {
    for (auto & arena : fArenas) {
      arena.tryFree(fTeardownMode);
    }

    if (fTeardownMode) {
      fArenas.erase(
          std::remove_if(fArenas.begin(), fArenas.end(), [](Arena& ar){return ar.ownedMemory == nullptr;}),
          fArenas.end());
    }
  }



  /// Test if pool is empty.
  bool empty() const
  {
    return std::all_of(fArenas.begin(), fArenas.end(), [](const Arena & ar) { return ar.empty(); });
  }



  /// Set pool to teardown mode (at program end).
  /// Will prune all empty arenas. Non-empty arenas will survive until all contained elements
  /// are deleted. They may therefore leak if not all elements are destructed.
  void teardown()
  {
    fTeardownMode = true;

    prune();
  }


  private:

  ////////////////////////////////////////////////////////////////////////////////////
  /// RooFit relies on unique pointers for RooArgSets. Here, memory
  /// has to be allocated until a completely new chunk of memory is encountered.
  /// As soon as RooXXXSets can be identified with a unique ID, this becomes obsolete.
  void newArena() {
    std::vector<Arena> failedAllocs;
    while (true) {
      Arena ar;
      if (std::none_of(fArenas.begin(), fArenas.end(),
          [&ar](Arena& other){return ar.memoryOverlaps(other);})) {
        fArenas.push_back(std::move(ar));
        break;
      }
      else {
        failedAllocs.push_back(std::move(ar));
      }
    }
  }



  std::vector<Arena> fArenas;
  bool               fTeardownMode{false};
};

#endif /* ROOFIT_ROOFITCORE_SRC_MEMPOOLFORROOSETS_H_ */
