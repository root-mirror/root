// @(#)root/treeplayer:$Id$
// Author: Axel Naumann, 2011-09-21

/*************************************************************************
 * Copyright (C) 1995-2013, Rene Brun and Fons Rademakers and al.        *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TTreeReader.h"

#include "TChain.h"
#include "TDirectory.h"
#include "TEntryList.h"
#include "TTreeCache.h"
#include "TTreeReaderValue.h"
#include "TFriendProxy.h"


// clang-format off
/**
 \class TTreeReader
 \ingroup treeplayer
 \brief A simple, robust and fast interface to read values from ROOT columnar datasets such as TTree, TChain or TNtuple

 TTreeReader is associated to TTreeReaderValue and TTreeReaderArray which are handles to concretely
 access the information in the dataset.

 Example code can be found in
  - tutorials/tree/hsimpleReader.C
  - tutorials/tree/h1analysisTreeReader.C
  - <a href="http://root.cern.ch/gitweb?p=roottest.git;a=tree;f=root/tree/reader;hb=HEAD">This example</a>

 You can generate a skeleton of `TTreeReaderValue<T>` and `TTreeReaderArray<T>` declarations
 for all of a tree's branches using `TTree::MakeSelector()`.

 Roottest contains an
 <a href="http://root.cern.ch/gitweb?p=roottest.git;a=tree;f=root/tree/reader;hb=HEAD">example</a>
 showing the full power.

A simpler analysis example can be found below: it histograms a function of the px and py branches.

~~~{.cpp}
// A simple TTreeReader use: read data from hsimple.root (written by hsimple.C)

#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

void hsimpleReader() {
   // Create a histogram for the values we read.
   TH1F("h1", "ntuple", 100, -4, 4);

   // Open the file containing the tree.
   TFile *myFile = TFile::Open("$ROOTSYS/tutorials/hsimple.root");

   // Create a TTreeReader for the tree, for instance by passing the
   // TTree's name and the TDirectory / TFile it is in.
   TTreeReader myReader("ntuple", myFile);

   // The branch "px" contains floats; access them as myPx.
   TTreeReaderValue<Float_t> myPx(myReader, "px");
   // The branch "py" contains floats, too; access those as myPy.
   TTreeReaderValue<Float_t> myPy(myReader, "py");

   // Loop over all entries of the TTree or TChain.
   while (myReader.Next()) {
      // Just access the data as if myPx and myPy were iterators (note the '*'
      // in front of them):
      myHist->Fill(*myPx + *myPy);
   }

   myHist->Draw();
}
~~~

A more complete example including error handling and a few combinations of
TTreeReaderValue and TTreeReaderArray would look like this:

~~~{.cpp}
#include <TFile.h>
#include <TH1.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include "TriggerInfo.h"
#include "Muon.h"
#include "Tau.h"

#include <vector>
#include <iostream>

bool CheckValue(ROOT::Internal::TTreeReaderValueBase& value) {
   if (value.GetSetupStatus() < 0) {
      std::cerr << "Error " << value.GetSetupStatus()
                << "setting up reader for " << value.GetBranchName() << '\n';
      return false;
   }
   return true;
}


// Analyze the tree "MyTree" in the file passed into the function.
// Returns false in case of errors.
bool analyze(TFile* file) {
   // Create a TTreeReader named "MyTree" from the given TDirectory.
   // The TTreeReader gives access to the TTree to the TTreeReaderValue and
   // TTreeReaderArray objects. It knows the current entry number and knows
   // how to iterate through the TTree.
   TTreeReader reader("MyTree", file);

   // Read a single float value in each tree entries:
   TTreeReaderValue<float> weight(reader, "event.weight");

   // Read a TriggerInfo object from the tree entries:
   TTreeReaderValue<TriggerInfo> triggerInfo(reader, "triggerInfo");

   //Read a vector of Muon objects from the tree entries:
   TTreeReaderValue<std::vector<Muon>> muons(reader, "muons");

   //Read the pT for all jets in the tree entry:
   TTreeReaderArray<double> jetPt(reader, "jets.pT");

   // Read the taus in the tree entry:
   TTreeReaderArray<Tau> taus(reader, "taus");


   // Now iterate through the TTree entries and fill a histogram.

   TH1F("hist", "TTreeReader example histogram", 10, 0., 100.);

   bool firstEntry = true;
   while (reader.Next()) {
      if (firstEntry) {
         // Check that branches exist and their types match our expectation.
         if (!CheckValue(weight)) return false;
         if (!CheckValue(triggerInfo)) return false;
         if (!CheckValue(muons)) return false;
         if (!CheckValue(jetPt)) return false;
         if (!CheckValue(taus)) return false;
         firstentry = false;
      }

      // Access the TriggerInfo object as if it's a pointer.
      if (!triggerInfo->hasMuonL1())
         continue;

      // Ditto for the vector<Muon>.
      if (!muons->size())
         continue;

      // Access the jetPt as an array, whether the TTree stores this as
      // a std::vector, std::list, TClonesArray or Jet* C-style array, with
      // fixed or variable array size.
      if (jetPt.GetSize() < 2 || jetPt[0] < 100)
         continue;

      // Access the array of taus.
      if (!taus.IsEmpty()) {
         // Access a float value - need to dereference as TTreeReaderValue
         // behaves like an iterator
         float currentWeight = *weight;
         for (const Tau& tau: taus) {
            hist->Fill(tau.eta(), currentWeight);
         }
      }
   } // TTree entry / event loop

   // Return true if we have iterated through all entries.
   return reader.GetEntryStatus() == TTreeReader::kEntryBeyondEnd;
}
~~~
*/
// clang-format on

ClassImp(TTreeReader);

using namespace ROOT::Internal;

// Provide some storage for the poor little symbol.
constexpr const char * const TTreeReader::fgEntryStatusText[TTreeReader::kEntryBeyondEnd + 1];

////////////////////////////////////////////////////////////////////////////////
/// Default constructor.  Call SetTree to connect to a TTree.

TTreeReader::TTreeReader() : fNotify(this) {}

////////////////////////////////////////////////////////////////////////////////
/// Access data from tree.
///
/// \param tree The TTree or TChain to read from
/// \param entryList It can be a single TEntryList with global entry numbers (supported, as
///                  an extension, also in the case of a TChain) or, if the first parameter
///                  is a TChain, a TEntryList with sub-TEntryLists with local entry numbers.
///                  In the latter case, the TEntryList must be associated to the TChain, as
///                  per chain.SetEntryList(&entryList).

TTreeReader::TTreeReader(TTree* tree, TEntryList* entryList /*= nullptr*/):
   fTree(tree),
   fEntryList(entryList),
   fNotify(this)
{
   if (!fTree) {
      ::Error("TTreeReader::TTreeReader", "TTree is NULL!");
   } else {
      Initialize();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Access data from the tree called keyname in the directory (e.g. TFile)
/// dir, or the current directory if dir is NULL. If keyname cannot be
/// found, or if it is not a TTree, IsInvalid() will return true.
///
/// \param keyname The name of the TTree to read from file
/// \param dir The TDirectory to read keyname from
/// \param entryList It can be a single TEntryList with global entry numbers (supported, as
///                  an extension, also in the case of a TChain) or, if the first parameter
///                  is a TChain, a TEntryList with sub-TEntryLists with local entry numbers.
///                  In the latter case, the TEntryList must be associated to the TChain, as
///                  per chain.SetEntryList(&entryList).

TTreeReader::TTreeReader(const char* keyname, TDirectory* dir, TEntryList* entryList /*= nullptr*/):
   fEntryList(entryList),
   fNotify(this)
{
   if (!dir) dir = gDirectory;
   dir->GetObject(keyname, fTree);
   if (!fTree) {
      std::string msg = "No TTree called ";
      msg += keyname;
      msg += " was found in the selected TDirectory.";
      Error("TTreeReader", "%s", msg.c_str());
   }
   Initialize();
}

////////////////////////////////////////////////////////////////////////////////
/// Tell all value readers that the tree reader does not exist anymore.

TTreeReader::~TTreeReader()
{
   for (std::deque<ROOT::Internal::TTreeReaderValueBase*>::const_iterator
           i = fValues.begin(), e = fValues.end(); i != e; ++i) {
      (*i)->MarkTreeReaderUnavailable();
   }
   if (fTree && fNotify.IsLinked())
      fNotify.RemoveLink(*fTree);

   // Need to clear the map of proxies before deleting the director otherwise
   // they will have a dangling pointer.
   fProxies.clear();

   for (auto feproxy: fFriendProxies) {
      delete feproxy;
   }
   fFriendProxies.clear();

   delete fDirector;
}

////////////////////////////////////////////////////////////////////////////////
/// Initialization of the director.

void TTreeReader::Initialize()
{
   fEntry = -1;
   if (!fTree) {
      fEntryStatus = kEntryNoTree;
      fLoadTreeStatus = kNoTree;
      return;
   }

   fLoadTreeStatus = kLoadTreeNone;
   if (fTree->InheritsFrom(TChain::Class())) {
      SetBit(kBitIsChain);
   } else if (fEntryList && fEntryList->GetLists()) {
      Error("Initialize", "We are not processing a TChain but the TEntryList contains sublists. Please "
                          "provide a simple TEntryList with no sublists instead.");
      fEntryStatus = kEntryNoTree;
      fLoadTreeStatus = kNoTree;
      return;
   }

   fDirector = new ROOT::Internal::TBranchProxyDirector(fTree, -1);

   if (!fNotify.IsLinked()) {
      fNotify.PrependLink(*fTree);

      if (fTree->GetTree()) {
         // The current TTree is already available.
         fSetEntryBaseCallingLoadTree = kTRUE;
         Notify();
         fSetEntryBaseCallingLoadTree = kFALSE;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Callback from TChain and TTree's LoadTree.

Bool_t TTreeReader::Notify()
{

   if (fSetEntryBaseCallingLoadTree) {
      if (fLoadTreeStatus == kExternalLoadTree) {
         // This can happen if someone switched trees behind us.
         // Likely cause: a TChain::LoadTree() e.g. from TTree::Process().
         // This means that "local" should be set!
         // There are two entities switching trees which is bad.
         Warning("SetEntryBase()",
                  "The current tree in the TChain %s has changed (e.g. by TTree::Process) "
                  "even though TTreeReader::SetEntry() was called, which switched the tree "
                  "again. Did you mean to call TTreeReader::SetLocalEntry()?",
                  fTree->GetName());
      }
      fLoadTreeStatus = kInternalLoadTree;
   } else {
      fLoadTreeStatus = kExternalLoadTree;
   }

   if (!fEntryList && fTree->GetEntryList() && !TestBit(kBitHaveWarnedAboutEntryListAttachedToTTree)) {
      Warning("SetEntryBase()",
              "The TTree / TChain has an associated TEntryList. "
              "TTreeReader ignores TEntryLists unless you construct the TTreeReader passing a TEntryList.");
      SetBit(kBitHaveWarnedAboutEntryListAttachedToTTree);
   }

   if (!fDirector->Notify()) {
      Error("SetEntryBase()", "There was an error while notifying the proxies.");
      return false;
   }

   if (fProxiesSet) {
      for (auto value: fValues) {
         value->NotifyNewTree(fTree->GetTree());
      }
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Tell readers we now have a tree.
/// fValues gets insertions during this loop (when parameterized arrays are read),
/// invalidating iterators. Use old-school counting instead.

Bool_t TTreeReader::SetProxies() {

   for (size_t i = 0; i < fValues.size(); ++i) {
      ROOT::Internal::TTreeReaderValueBase* reader = fValues[i];
      reader->CreateProxy();
      if (!reader->GetProxy()){
         return kFALSE;
      }

   }
   // If at least one proxy was there and no error occurred, we assume the proxies to be set.
   fProxiesSet = !fValues.empty();

   // Now we need to properly set the TTreeCache. We do this in steps:
   // 1. We set the entry range according to the entry range of the TTreeReader
   // 2. We add to the cache the branches identifying them by the name the user provided
   //    upon creation of the TTreeReader{Value, Array}s
   // 3. We stop the learning phase.
   // Operations 1, 2 and 3 need to happen in this order. See: https://sft.its.cern.ch/jira/browse/ROOT-9773?focusedCommentId=87837
   if (fProxiesSet) {
      const auto curFile = fTree->GetCurrentFile();
      if (curFile && fTree->GetTree()->GetReadCache(curFile, true)) {
         if (!(-1LL == fEndEntry && 0ULL == fBeginEntry)) {
            // We need to avoid to pass -1 as end entry to the SetCacheEntryRange method
            const auto lastEntry = (-1LL == fEndEntry) ? fTree->GetEntriesFast() : fEndEntry;
            fTree->SetCacheEntryRange(fBeginEntry, lastEntry);
         }
         for (auto value: fValues) {
            fTree->AddBranchToCache(value->GetProxy()->GetBranchName(), true);
         }
         fTree->StopCacheLearningPhase();
      }
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the range of entries to be loaded by `Next()`; end will not be loaded.
///
/// If end <= begin, `end` is ignored (set to `-1`, i.e. will run on all entries from `begin` onwards).
///
/// Example:
///
/// ~~~ {.cpp}
/// reader.SetEntriesRange(3, 5);
/// while (reader.Next()) {
///   // Will load entries 3 and 4.
/// }
/// ~~~
///
/// Note that if a TEntryList is present, beginEntry and endEntry refer to the beginEntry-th/endEntry-th entries of the
/// TEntryList (or the main TEntryList in case it has sub-entrylists). In other words, SetEntriesRange can
/// be used to only loop over part of the TEntryList, but not to further restrict the actual TTree/TChain entry numbers
/// considered.
///
/// \param beginEntry The first entry to be loaded by `Next()`.
/// \param endEntry   The entry where `Next()` will return kFALSE, not loading it.

TTreeReader::EEntryStatus TTreeReader::SetEntriesRange(Long64_t beginEntry, Long64_t endEntry)
{
   if (beginEntry < 0)
      return kEntryNotFound;
   // Complain if the entries number is larger than the tree's / chain's / entry
   // list's number of entries, unless it's a TChain and "max entries" is
   // uninitialized (i.e. TTree::kMaxEntries).
   if (beginEntry >= GetEntries(false) && !(IsChain() && GetEntries(false) == TTree::kMaxEntries)) {
      Error("SetEntriesRange()", "first entry out of range 0..%lld", GetEntries(false));
      return kEntryNotFound;
   }

   // Update data members to correctly reflect the defined range
   if (endEntry > beginEntry)
      fEndEntry = endEntry;
   else
      fEndEntry = -1;

   fBeginEntry = beginEntry;

   if (beginEntry - 1 < 0) 
      // Reset the cache if reading from the first entry of the tree
      Restart();
   else {
      // Load the first entry in the range. SetEntry() will also call SetProxies(),
      // thus adding all the branches to the cache and triggering the learning phase.
      EEntryStatus es = SetEntry(beginEntry - 1);
      if (es != kEntryValid) {
         Error("SetEntriesRange()", "Error setting first entry %lld: %s",
               beginEntry, fgEntryStatusText[(int)es]);
         return es;
      }
   }

   return kEntryValid;
}

void TTreeReader::Restart() {
   fDirector->SetReadEntry(-1);
   fProxiesSet = false; // we might get more value readers, meaning new proxies.
   fEntry = -1;
   if (const auto curFile = fTree->GetCurrentFile()) {
      if (auto tc = fTree->GetTree()->GetReadCache(curFile, true)) {
         tc->DropBranch("*", true);
         tc->ResetCache();
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the number of entries of the TEntryList if one is provided, else
/// of the TTree / TChain, independent of a range set by SetEntriesRange()
/// by calling TTree/TChain::GetEntriesFast.


Long64_t TTreeReader::GetEntries() const {
   if (fEntryList)
      return fEntryList->GetN();
   if (!fTree)
      return -1;
   return fTree->GetEntriesFast();
}


////////////////////////////////////////////////////////////////////////////////
/// Returns the number of entries of the TEntryList if one is provided, else
/// of the TTree / TChain, independent of a range set by SetEntriesRange().
///
/// \param force If `IsChain()` and `force`, determines whether all TFiles of
///   this TChain should be opened to determine the exact number of entries
/// of the TChain. If `!IsChain()`, `force` is ignored.

Long64_t TTreeReader::GetEntries(Bool_t force)  {
   if (fEntryList)
      return fEntryList->GetN();
   if (!fTree)
      return -1;
   if (force) {
      fSetEntryBaseCallingLoadTree = kTRUE;
      auto res = fTree->GetEntries();
      // Go back to where we were:
      fTree->LoadTree(GetCurrentEntry());
      fSetEntryBaseCallingLoadTree = kFALSE;
      return res;
   }
   return fTree->GetEntriesFast();
}



////////////////////////////////////////////////////////////////////////////////
/// Load an entry into the tree, return the status of the read.
/// For chains, entry is the global (i.e. not tree-local) entry number, unless
/// `local` is `true`, in which case `entry` specifies the entry number within
/// the current tree. This is needed for instance for TSelector::Process().

TTreeReader::EEntryStatus TTreeReader::SetEntryBase(Long64_t entry, Bool_t local)
{
   if (IsInvalid()) {
      fEntryStatus = kEntryNoTree;
      fEntry = -1;
      return fEntryStatus;
   }

   fEntry = entry;

   Long64_t entryAfterList = entry;
   if (fEntryList) {
      if (entry >= fEntryList->GetN()) {
         // Passed the end of the chain, Restart() was not called:
         // don't try to load entries anymore. Can happen in these cases:
         // while (tr.Next()) {something()};
         // while (tr.Next()) {somethingelse()}; // should not be calling somethingelse().
         fEntryStatus = kEntryNotFound;
         return fEntryStatus;
      }
      if (entry >= 0) {
         if (fEntryList->GetLists()) {
            R__ASSERT(IsChain());
            int treenum = -1;
            entryAfterList = fEntryList->GetEntryAndTree(entry, treenum);
            entryAfterList += static_cast<TChain *>(fTree)->GetTreeOffset()[treenum];
            // We always translate local entry numbers to global entry numbers for TChain+TEntryList with sublists
            local = false;
         } else {
            // Could be a TTree or a TChain (TTreeReader also supports single TEntryLists for TChains).
            // In both cases, we are using the global entry numbers coming from the single TEntryList.
            entryAfterList = fEntryList->GetEntry(entry);
         }
      }
   }

   TTree* treeToCallLoadOn = local ? fTree->GetTree() : fTree;

   fSetEntryBaseCallingLoadTree = kTRUE;
   const Long64_t loadResult = treeToCallLoadOn->LoadTree(entryAfterList);
   fSetEntryBaseCallingLoadTree = kFALSE;

   if (loadResult < 0) {
      // ROOT-9628 We cover here the case when:
      // - We deal with a TChain
      // - The last file is opened
      // - The TTree is not correctly loaded
      // The system is robust against issues with TTrees associated to the chain
      // when they are not at the end of it.
      if (loadResult == -3 && TestBit(kBitIsChain) && !fTree->GetTree()) {
         fDirector->Notify();
         if (fProxiesSet) {
            for (auto value: fValues) {
               value->NotifyNewTree(fTree->GetTree());
            }
         }
         Warning("SetEntryBase()",
               "There was an issue opening the last file associated to the TChain "
               "being processed.");
         fEntryStatus = kEntryChainFileError;
         return fEntryStatus;
      }

      if (loadResult == -2) {
         fDirector->Notify();
         if (fProxiesSet) {
            for (auto value: fValues) {
               value->NotifyNewTree(fTree->GetTree());
            }
         }
         fEntryStatus = kEntryNotFound;
         return fEntryStatus;
      }

      if (loadResult == -1) {
         // The chain is empty
         fEntryStatus = kEntryNotFound;
         return fEntryStatus;
      }

      if (loadResult == -4) {
         // The TChainElement corresponding to the entry is missing or
         // the TTree is missing from the file.
         fDirector->Notify();
         if (fProxiesSet) {
            for (auto value: fValues) {
               value->NotifyNewTree(fTree->GetTree());
            }
         }
         fEntryStatus = kEntryNotFound;
         return fEntryStatus;
      }

      Warning("SetEntryBase()",
              "Unexpected error '%lld' in %s::LoadTree", loadResult,
              treeToCallLoadOn->IsA()->GetName());

      fEntryStatus = kEntryUnknownError;
      return fEntryStatus;
   }

   if (!fProxiesSet) {
      if (!SetProxies()) {
         fEntryStatus = kEntryDictionaryError;
         return fEntryStatus;
      }
   }

   if ((fEndEntry >= 0 && entry >= fEndEntry) || (fEntry >= fTree->GetEntriesFast())) {
      fEntryStatus = kEntryBeyondEnd;
      return fEntryStatus;
   }
   fDirector->SetReadEntry(loadResult);
   fEntryStatus = kEntryValid;
   return fEntryStatus;
}

////////////////////////////////////////////////////////////////////////////////
/// Set (or update) the which tree to read from. `tree` can be
/// a TTree or a TChain.

void TTreeReader::SetTree(TTree* tree, TEntryList* entryList /*= nullptr*/)
{
   fTree = tree;
   fEntryList = entryList;
   fEntry = -1;

   if (fTree) {
      fLoadTreeStatus = kLoadTreeNone;
      SetBit(kBitIsChain, fTree->InheritsFrom(TChain::Class()));
   } else {
      fLoadTreeStatus = kNoTree;
   }

   if (!fDirector) {
      Initialize();
   }
   else {
      fDirector->SetTree(fTree);
      fDirector->SetReadEntry(-1);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set (or update) the which tree to read from, passing the name of a tree in a
/// directory.
///
/// \param keyname - name of the tree in `dir`
/// \param dir - the `TDirectory` to load `keyname` from (or gDirectory if `nullptr`)
/// \param entryList - the `TEntryList` to attach to the `TTreeReader`.

void TTreeReader::SetTree(const char* keyname, TDirectory* dir, TEntryList* entryList /*= nullptr*/)
{
   TTree* tree = nullptr;
   if (!dir)
      dir = gDirectory;
   dir->GetObject(keyname, tree);
   SetTree(tree, entryList);
}

////////////////////////////////////////////////////////////////////////////////
/// Add a value reader for this tree.

Bool_t TTreeReader::RegisterValueReader(ROOT::Internal::TTreeReaderValueBase* reader)
{
   if (fProxiesSet) {
      Error("RegisterValueReader",
            "Error registering reader for %s: TTreeReaderValue/Array objects must be created before the call to Next() / SetEntry() / SetLocalEntry(), or after TTreeReader::Restart()!",
            reader->GetBranchName());
      return false;
   }
   fValues.push_back(reader);
   return true;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove a value reader for this tree.

void TTreeReader::DeregisterValueReader(ROOT::Internal::TTreeReaderValueBase* reader)
{
   std::deque<ROOT::Internal::TTreeReaderValueBase*>::iterator iReader
      = std::find(fValues.begin(), fValues.end(), reader);
   if (iReader == fValues.end()) {
      Error("DeregisterValueReader", "Cannot find reader of type %s for branch %s", reader->GetDerivedTypeName(), reader->fBranchName.Data());
      return;
   }
   fValues.erase(iReader);
}
