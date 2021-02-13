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

#include "TTree.h"
#include "TTreeReaderUtils.h"
#include "TNotifyLink.h"

#include <deque>
#include <iterator>
#include <unordered_map>
#include <string>

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
      kEntryBeyondEnd, ///< last entry loop has reached its end
      kEntryBadReader, ///< One of the readers was not successfully initialized.
      kEntryUnknownError ///< LoadTree return less than -4, likely a 'newer' error code.
   };

   enum ELoadTreeStatus {
      kNoTree = 0,       ///< default state, no TTree is connected (formerly 'Zombie' state)
      kLoadTreeNone,     ///< Notify has not been called yet.
      kInternalLoadTree, ///< Notify/LoadTree was last called from SetEntryBase
      kExternalLoadTree  ///< User code called LoadTree directly.
   };

   static constexpr const char * const fgEntryStatusText[kEntryBeyondEnd + 1] = {
      "valid entry",
      "the tree does not exist",
      "the tree entry number does not exist",
      "cannot access chain element",
      "problem in opening a chain's file",
      "problem reading dictionary info from tree",
      "last entry loop has reached its end"
   };

   TTreeReader();

   TTreeReader(TTree* tree, TEntryList* entryList = nullptr);
   TTreeReader(const char* keyname, TDirectory* dir, TEntryList* entryList = nullptr);
   TTreeReader(const char *keyname, TEntryList *entryList = nullptr) : TTreeReader(keyname, nullptr, entryList) {}

   ~TTreeReader();

   void SetTree(TTree* tree, TEntryList* entryList = nullptr);
   void SetTree(const char* keyname, TEntryList* entryList = nullptr) {
      SetTree(keyname, nullptr, entryList);
   }
   void SetTree(const char* keyname, TDirectory* dir, TEntryList* entryList = nullptr);

   Bool_t IsChain() const { return TestBit(kBitIsChain); }

   Bool_t IsInvalid() const { return fLoadTreeStatus == kNoTree; }

   TTree* GetTree() const { return fTree; }
   TEntryList* GetEntryList() const { return fEntryList; }

   ///\{ \name Entry setters

   /// Move to the next entry (or index of the TEntryList if that is set).
   ///
   /// \return false if the previous entry was already the last entry. This allows
   ///   the function to be used in `while (reader.Next()) { ... }`
   Bool_t Next() {
      return SetEntry(GetCurrentEntry() + 1) == kEntryValid;
   }

   /// Set the next entry (or index of the TEntryList if that is set).
   ///
   /// \param entry If not TEntryList is set, the entry is a global entry (i.e.
   /// not the entry number local to the chain's current tree).
   /// \returns the `entry`'s read status, i.e. whether the entry is available.
   EEntryStatus SetEntry(Long64_t entry) { return SetEntryBase(entry, kFALSE); }

   /// Set the next local tree entry. If a TEntryList is set, this function is
   /// equivalent to `SetEntry()`.
   ///
   /// \param entry Entry number of the TChain's current TTree. This is the
   /// entry number passed for instance by `TSelector::Process(entry)`, i.e.
   /// within `TSelector::Process()` always use `SetLocalEntry()` and not
   /// `SetEntry()`!
   /// \return the `entry`'s read status, i.e. whether the entry is available.
   EEntryStatus SetLocalEntry(Long64_t entry) { return SetEntryBase(entry, kTRUE); }

   ///  Set the begin and end entry numbers
   ///
   /// \param beginEntry The first entry that `Next()` will load.
   /// \param endEntry The entry that `Next()` will return `kFALSE` on (i.e. not
   ///   load anymore).
   EEntryStatus SetEntriesRange(Long64_t beginEntry, Long64_t endEntry);

   ///  Get the begin and end entry numbers
   ///
   /// \return a pair contained the begin and end entry numbers.
   std::pair<Long64_t, Long64_t> GetEntriesRange() const { return std::make_pair(fBeginEntry, fEndEntry); }

   /// Restart a Next() loop from entry 0 (of TEntryList index 0 of fEntryList is set).
   void Restart();

   ///\}

   EEntryStatus GetEntryStatus() const { return fEntryStatus; }

   Long64_t GetEntries() const;
   Long64_t GetEntries(Bool_t force);

   /// Returns the index of the current entry being read.
   ///
   /// If `IsChain()`, the returned index corresponds to the global entry number
   /// (i.e. not the entry number local to the chain's current tree).
   /// If `fEntryList`, the returned index corresponds to an index in the
   /// TEntryList; to translate to the TChain's / TTree's entry number pass it
   /// through `reader.GetEntryList()->GetEntry(reader.GetCurrentEntry())`.
   Long64_t GetCurrentEntry() const { return fEntry; }

   Bool_t Notify();

   /// Return an iterator to the 0th TTree entry.
   Iterator_t begin() {
      return Iterator_t(*this, 0);
   }
   /// Return an iterator beyond the last TTree entry.
   Iterator_t end() const { return Iterator_t(); }

protected:
   using NamedProxies_t = std::unordered_map<std::string, std::unique_ptr<ROOT::Internal::TNamedBranchProxy>>;
   void Initialize();
   ROOT::Internal::TNamedBranchProxy* FindProxy(const char* branchname) const
   {
      const auto proxyIt = fProxies.find(branchname);
      return fProxies.end() != proxyIt ? proxyIt->second.get() : nullptr;
   }

   void AddProxy(ROOT::Internal::TNamedBranchProxy *p)
   {
      auto bpName = p->GetName();
#ifndef NDEBUG
      if (fProxies.end() != fProxies.find(bpName)) {
         std::string err = "A proxy with key " + std::string(bpName) + " was already stored!";
         throw std::runtime_error(err);
      }
#endif

      fProxies[bpName].reset(p);
   }

   Bool_t RegisterValueReader(ROOT::Internal::TTreeReaderValueBase* reader);
   void DeregisterValueReader(ROOT::Internal::TTreeReaderValueBase* reader);

   EEntryStatus SetEntryBase(Long64_t entry, Bool_t local);

   Bool_t SetProxies();

private:

   std::string GetProxyKey(const char *branchname)
   {
      std::string key(branchname);
      //key += reinterpret_cast<std::uintptr_t>(fTree);
      return key;
   }

   enum EStatusBits {
      kBitIsChain = BIT(14), ///< our tree is a chain
      kBitHaveWarnedAboutEntryListAttachedToTTree = BIT(15), ///< the tree had a TEntryList and we have warned about that
      kBitSetEntryBaseCallingLoadTree = BIT(16) ///< SetEntryBase is in the process of calling TChain/TTree::LoadTree.
   };

   TTree* fTree = nullptr; ///< tree that's read
   TEntryList* fEntryList = nullptr; ///< entry list to be used
   EEntryStatus fEntryStatus = kEntryNotLoaded; ///< status of most recent read request
   ELoadTreeStatus fLoadTreeStatus = kNoTree;   ///< Indicator on how LoadTree was called 'last' time.
   TNotifyLink<TTreeReader> fNotify; // Callback object used by the TChain to update this proxy
   ROOT::Internal::TBranchProxyDirector* fDirector = nullptr; ///< proxying director, owned
   std::deque<ROOT::Internal::TFriendProxy*> fFriendProxies; ///< proxying for friend TTrees, owned
   std::deque<ROOT::Internal::TTreeReaderValueBase*> fValues; ///< readers that use our director
   NamedProxies_t fProxies; ///< attached ROOT::TNamedBranchProxies; owned

   Long64_t fEntry = -1; ///< Current (non-local) entry of fTree or of fEntryList if set.

   /// The end of the entry loop. When set (i.e. >= 0), it provides a way
   /// to stop looping over the TTree when we reach a certain entry: Next()
   /// returns kFALSE when GetCurrentEntry() reaches fEndEntry.
   Long64_t fEndEntry = -1LL;
   Long64_t fBeginEntry = 0LL; ///< This allows us to propagate the range to the TTreeCache
   Bool_t fProxiesSet = kFALSE; ///< True if the proxies have been set, false otherwise
   Bool_t fSetEntryBaseCallingLoadTree = kFALSE; ///< True if during the LoadTree execution triggered by SetEntryBase.

   friend class ROOT::Internal::TTreeReaderValueBase;
   friend class ROOT::Internal::TTreeReaderArrayBase;

   ClassDef(TTreeReader, 0); // A simple interface to read trees
};

#endif // defined TTreeReader
