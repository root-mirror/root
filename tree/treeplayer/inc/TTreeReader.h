// @(#)root/tree:$Id$
// Author: Axel Naumann, 2010-08-02

/*************************************************************************
 * Copyright (C) 1995-2013, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TTreeReader
#define ROOT_TTreeReader


////////////////////////////////////////////////////////////////////////////
//                                                                        //
// TTreeReader                                                            //
//                                                                        //
// A simple interface for reading trees or chains.                        //
//                                                                        //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef ROOT_THashTable
#include "THashTable.h"
#endif
#ifndef ROOT_TTree
#include "TTree.h"
#endif
#ifndef ROOT_TTreeReaderUtils
#include "TTreeReaderUtils.h"
#endif

#include <deque>
#include <iterator>

class TDictionary;
class TDirectory;
class TFileCollection;

namespace ROOT {
namespace Internal {
   class TBranchProxyDirector;
}
}

class TTreeReader: public TObject {
public:

   ///\class TTreeReader::Iterator_t
   /// Iterate through the entries of a TTree.
   ///
   /// This iterator drives the associated TTreeReader; its
   /// dereferencing (and actually even the iteration) will
   /// set the entry number represented by this iterator.
   /// It does not really represent a data element; it simply
   /// returns the entry number (or -1 once the end of the tree
   /// is reached).
   class Iterator_t:
      public std::iterator<std::input_iterator_tag, const Long64_t, Long64_t> {
   private:
      Long64_t fEntry; ///< Entry number of the tree referenced by this iterator; -1 is invalid.
      TTreeReader* fReader; ///< The reader we select the entries on.

      /// Whether the iterator points to a valid entry.
      bool IsValid() const { return fEntry >= 0; }

   public:
      /// Default-initialize the iterator as "past the end".
      Iterator_t(): fEntry(-1), fReader() {}

      /// Initialize the iterator with the reader it steers and a
      /// tree entry number; -1 is invalid.
      Iterator_t(TTreeReader& reader, Long64_t entry):
         fEntry(entry), fReader(&reader) {}

      /// Compare two iterators for equality.
      bool operator==(const Iterator_t& lhs) const {
         // From C++14: value initialized (past-end) it compare equal.
         if (!IsValid() && !lhs.IsValid()) return true;
         return fEntry == lhs.fEntry && fReader == lhs.fReader;
      }

      /// Compare two iterators for inequality.
      bool operator!=(const Iterator_t& lhs) const {
         return !(*this == lhs);
      }

      /// Increment the iterator (postfix i++).
      Iterator_t operator++(int) {
         Iterator_t ret = *this;
         this->operator++();
         return ret;
      }

      /// Increment the iterator (prefix ++i).
      Iterator_t& operator++() {
         if (IsValid()) {
            ++fEntry;
            // Force validity check of new fEntry.
            this->operator*();
            // Don't set the old entry: op* will if needed, and
            // in most cases it just adds a lot of spinning back
            // and forth: in most cases teh sequence is ++i; *i.
         }
         return *this;
      }

      /// Set the entry number in the reader and return it.
      const Long64_t& operator*() {
         if (IsValid()) {
            // If we cannot access that entry, mark the iterator invalid.
            if (fReader->SetEntry(fEntry) != kEntryValid) {
               fEntry = -1;
            }
         }
         // There really is no data in this iterator; return the number.
         return fEntry;
      }

      const Long64_t& operator*() const {
         return **const_cast<Iterator_t*>(this);
      }
   };

   typedef Iterator_t iterator;

   enum EEntryStatus {
      kEntryValid = 0, ///< data read okay
      kEntryNotLoaded, ///< no entry has been loaded yet
      kEntryNoTree, ///< the tree does not exist
      kEntryNotFound, ///< the tree entry number does not exist
      kEntryChainSetupError, ///< problem in accessing a chain element, e.g. file without the tree
      kEntryChainFileError, ///< problem in opening a chain's file
      kEntryDictionaryError, ///< problem reading dictionary info from tree
      kEntryLast, ///< last entry was reached
   };

   TTreeReader():
      fTree(0),
      fDirectory(0),
      fEntryStatus(kEntryNoTree),
      fMostRecentTreeNumber(-1),
      fDirector(0),
      fLastEntry(-1),
      fProxiesSet(false)
   {}

   TTreeReader(TTree* tree);
   TTreeReader(const char* keyname, TDirectory* dir = NULL );

   ~TTreeReader();

   void SetTree(TTree* tree);
   void SetTree(const char* /*keyname*/, TDirectory* /*dir = NULL*/ ) { Error("SetTree()", "Not Implemented!");};
   void SetChain(const char* /*keyname*/, TFileCollection* /*files*/ ) { Error("SetChain()", "Not Implemented!");};

   Bool_t IsChain() const { return TestBit(kBitIsChain); }

   Bool_t Next() { return SetEntry(GetCurrentEntry() + 1) == kEntryValid; }
   EEntryStatus SetEntry(Long64_t entry) { return SetEntryBase(entry, kFALSE); }
   EEntryStatus SetLocalEntry(Long64_t entry) { return SetEntryBase(entry, kTRUE); }
   void SetLastEntry(Long64_t entry) { fLastEntry = entry; }
   EEntryStatus SetEntriesRange(Long64_t first, Long64_t last);
   void Restart();

   EEntryStatus GetEntryStatus() const { return fEntryStatus; }

   TTree* GetTree() const { return fTree; }
   Long64_t GetEntries(Bool_t force) const { return fTree ? (force ? fTree->GetEntries() : fTree->GetEntriesFast() ) : -1; }
   Long64_t GetCurrentEntry() const;

   /// Return an iterator to the 0th TTree entry.
   Iterator_t begin() {
      return Iterator_t(*this, 0);
   }
   Iterator_t end() const { return Iterator_t(); }

protected:
   void Initialize();
   ROOT::Internal::TNamedBranchProxy* FindProxy(const char* branchname) const {
      return (ROOT::Internal::TNamedBranchProxy*) fProxies.FindObject(branchname); }
   TCollection* GetProxies() { return &fProxies; }

   Bool_t RegisterValueReader(ROOT::Internal::TTreeReaderValueBase* reader);
   void DeregisterValueReader(ROOT::Internal::TTreeReaderValueBase* reader);

   EEntryStatus SetEntryBase(Long64_t entry, Bool_t local);

private:

   enum EPropertyBits {
      kBitIsChain = BIT(14) ///< our tree is a chain
   };

   TTree* fTree; ///< tree that's read
   TDirectory* fDirectory; ///< directory (or current file for chains)
   EEntryStatus fEntryStatus; ///< status of most recent read request
   Int_t fMostRecentTreeNumber; ///< TTree::GetTreeNumber() of the most recent tree
   ROOT::Internal::TBranchProxyDirector* fDirector; ///< proxying director, owned
   std::deque<ROOT::Internal::TTreeReaderValueBase*> fValues; ///< readers that use our director
   THashTable   fProxies; ///< attached ROOT::TNamedBranchProxies; owned

   /// The last entry to be processed. When set (i.e. >= 0), it provides a way
   /// to stop looping over the TTree when we reach a certain entry: Next()
   /// returns kEntryLast when GetCurrentEntry() reaches fLastEntry
   Long64_t fLastEntry;
   Bool_t fProxiesSet; ///< True if the proxies have been set, false otherwise

   friend class ROOT::Internal::TTreeReaderValueBase;
   friend class ROOT::Internal::TTreeReaderArrayBase;

   ClassDef(TTreeReader, 0); // A simple interface to read trees
};

#endif // defined TTreeReader
