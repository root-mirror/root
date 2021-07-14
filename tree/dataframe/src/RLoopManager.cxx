/*************************************************************************
 * Copyright (C) 1995-2021, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "RConfigure.h" // R__USE_IMT
#include "ROOT/RDataSource.hxx"
#include "ROOT/RDF/GraphNode.hxx"
#include "ROOT/RDF/RActionBase.hxx"
#include "ROOT/RDF/RFilterBase.hxx"
#include "ROOT/RDF/RLoopManager.hxx"
#include "ROOT/RDF/RRangeBase.hxx"
#include "ROOT/RDF/RSlotStack.hxx"
#include "ROOT/RLogger.hxx"
#include "RtypesCore.h" // Long64_t
#include "TStopwatch.h"
#include "TBranchElement.h"
#include "TBranchObject.h"
#include "TChain.h"
#include "TEntryList.h"
#include "TFile.h"
#include "TFriendElement.h"
#include "TInterpreter.h"
#include "TROOT.h" // IsImplicitMTEnabled
#include "TTreeReader.h"
#include "TTree.h" // For MaxTreeSizeRAII. Revert when #6640 will be solved.

#ifdef R__USE_IMT
#include "ROOT/TThreadExecutor.hxx"
#include "ROOT/TTreeProcessorMT.hxx"
#endif

#include <algorithm>
#include <atomic>
#include <exception>
#include <functional>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <sstream>
#include <thread>
#include <unordered_map>
#include <vector>
#include <set>
#include <limits> // For MaxTreeSizeRAII. Revert when #6640 will be solved.

using namespace ROOT::Detail::RDF;
using namespace ROOT::Internal::RDF;

namespace {
/// A helper function that returns all RDF code that is currently scheduled for just-in-time compilation.
/// This allows different RLoopManager instances to share these data.
/// We want RLoopManagers to be able to add their code to a global "code to execute via cling",
/// so that, lazily, we can jit everything that's needed by all RDFs in one go, which is potentially
/// much faster than jitting each RLoopManager's code separately.
static std::string &GetCodeToJit()
{
   static std::string code;
   return code;
}

static bool ContainsLeaf(const std::set<TLeaf *> &leaves, TLeaf *leaf)
{
   return (leaves.find(leaf) != leaves.end());
}

///////////////////////////////////////////////////////////////////////////////
/// This overload does not check whether the leaf/branch is already in bNamesReg. In case this is a friend leaf/branch,
/// `allowDuplicates` controls whether we add both `friendname.bname` and `bname` or just the shorter version.
static void InsertBranchName(std::set<std::string> &bNamesReg, ColumnNames_t &bNames, const std::string &branchName,
                             const std::string &friendName, bool allowDuplicates)
{
   if (!friendName.empty()) {
      // In case of a friend tree, users might prepend its name/alias to the branch names
      const auto friendBName = friendName + "." + branchName;
      if (bNamesReg.insert(friendBName).second)
         bNames.push_back(friendBName);
   }

   if (allowDuplicates || friendName.empty()) {
      if (bNamesReg.insert(branchName).second)
         bNames.push_back(branchName);
   }
}

///////////////////////////////////////////////////////////////////////////////
/// This overload makes sure that the TLeaf has not been already inserted.
static void InsertBranchName(std::set<std::string> &bNamesReg, ColumnNames_t &bNames, const std::string &branchName,
                             const std::string &friendName, std::set<TLeaf *> &foundLeaves, TLeaf *leaf,
                             bool allowDuplicates)
{
   const bool canAdd = allowDuplicates ? true : !ContainsLeaf(foundLeaves, leaf);
   if (!canAdd) {
      return;
   }

   InsertBranchName(bNamesReg, bNames, branchName, friendName, allowDuplicates);

   foundLeaves.insert(leaf);
}

static void ExploreBranch(TTree &t, std::set<std::string> &bNamesReg, ColumnNames_t &bNames, TBranch *b,
                          std::string prefix, std::string &friendName, bool allowDuplicates)
{
   for (auto sb : *b->GetListOfBranches()) {
      TBranch *subBranch = static_cast<TBranch *>(sb);
      auto subBranchName = std::string(subBranch->GetName());
      auto fullName = prefix + subBranchName;

      std::string newPrefix;
      if (!prefix.empty())
         newPrefix = fullName + ".";

      ExploreBranch(t, bNamesReg, bNames, subBranch, newPrefix, friendName, allowDuplicates);

      auto branchDirectlyFromTree = t.GetBranch(fullName.c_str());
      if (!branchDirectlyFromTree)
         branchDirectlyFromTree = t.FindBranch(fullName.c_str()); // try harder
      if (branchDirectlyFromTree)
         InsertBranchName(bNamesReg, bNames, std::string(branchDirectlyFromTree->GetFullName()), friendName,
                          allowDuplicates);

      if (t.GetBranch(subBranchName.c_str()))
         InsertBranchName(bNamesReg, bNames, subBranchName, friendName, allowDuplicates);
   }
}

static void GetBranchNamesImpl(TTree &t, std::set<std::string> &bNamesReg, ColumnNames_t &bNames,
                               std::set<TTree *> &analysedTrees, std::string &friendName, bool allowDuplicates)
{
   std::set<TLeaf *> foundLeaves;
   if (!analysedTrees.insert(&t).second) {
      return;
   }

   const auto branches = t.GetListOfBranches();
   // Getting the branches here triggered the read of the first file of the chain if t is a chain.
   // We check if a tree has been successfully read, otherwise we throw (see ROOT-9984) to avoid further
   // operations
   if (!t.GetTree()) {
      std::string err("GetBranchNames: error in opening the tree ");
      err += t.GetName();
      throw std::runtime_error(err);
   }
   if (branches) {
      for (auto b : *branches) {
         TBranch *branch = static_cast<TBranch *>(b);
         const auto branchName = std::string(branch->GetName());
         if (branch->IsA() == TBranch::Class()) {
            // Leaf list
            auto listOfLeaves = branch->GetListOfLeaves();
            if (listOfLeaves->GetEntriesUnsafe() == 1) {
               auto leaf = static_cast<TLeaf *>(listOfLeaves->UncheckedAt(0));
               InsertBranchName(bNamesReg, bNames, branchName, friendName, foundLeaves, leaf, allowDuplicates);
            }

            for (auto leaf : *listOfLeaves) {
               auto castLeaf = static_cast<TLeaf *>(leaf);
               const auto leafName = std::string(leaf->GetName());
               const auto fullName = branchName + "." + leafName;
               InsertBranchName(bNamesReg, bNames, fullName, friendName, foundLeaves, castLeaf, allowDuplicates);
            }
         } else if (branch->IsA() == TBranchObject::Class()) {
            // TBranchObject
            ExploreBranch(t, bNamesReg, bNames, branch, branchName + ".", friendName, allowDuplicates);
            InsertBranchName(bNamesReg, bNames, branchName, friendName, allowDuplicates);
         } else {
            // TBranchElement
            // Check if there is explicit or implicit dot in the name

            bool dotIsImplied = false;
            auto be = dynamic_cast<TBranchElement *>(b);
            if (!be)
               throw std::runtime_error("GetBranchNames: unsupported branch type");
            // TClonesArray (3) and STL collection (4)
            if (be->GetType() == 3 || be->GetType() == 4)
               dotIsImplied = true;

            if (dotIsImplied || branchName.back() == '.')
               ExploreBranch(t, bNamesReg, bNames, branch, "", friendName, allowDuplicates);
            else
               ExploreBranch(t, bNamesReg, bNames, branch, branchName + ".", friendName, allowDuplicates);

            InsertBranchName(bNamesReg, bNames, branchName, friendName, allowDuplicates);
         }
      }
   }

   // The list of friends needs to be accessed via GetTree()->GetListOfFriends()
   // (and not via GetListOfFriends() directly), otherwise when `t` is a TChain we
   // might not recover the list correctly (https://github.com/root-project/root/issues/6741).
   auto friendTrees = t.GetTree()->GetListOfFriends();

   if (!friendTrees)
      return;

   for (auto friendTreeObj : *friendTrees) {
      auto friendTree = ((TFriendElement *)friendTreeObj)->GetTree();

      std::string frName;
      auto alias = t.GetFriendAlias(friendTree);
      if (alias != nullptr)
         frName = std::string(alias);
      else
         frName = std::string(friendTree->GetName());

      GetBranchNamesImpl(*friendTree, bNamesReg, bNames, analysedTrees, frName, allowDuplicates);
   }
}

static void ThrowIfNSlotsChanged(unsigned int nSlots)
{
   const auto currentSlots = RDFInternal::GetNSlots();
   if (currentSlots != nSlots) {
      std::string msg = "RLoopManager::Run: when the RDataFrame was constructed the number of slots required was " +
                        std::to_string(nSlots) + ", but when starting the event loop it was " +
                        std::to_string(currentSlots) + ".";
      if (currentSlots > nSlots)
         msg += " Maybe EnableImplicitMT() was called after the RDataFrame was constructed?";
      else
         msg += " Maybe DisableImplicitMT() was called after the RDataFrame was constructed?";
      throw std::runtime_error(msg);
   }
}

/**
\struct MaxTreeSizeRAII
\brief Scope-bound change of `TTree::fgMaxTreeSize`.

This RAII object stores the current value result of `TTree::GetMaxTreeSize`,
changes it to maximum at construction time and restores it back at destruction
time. Needed for issue #6523 and should be reverted when #6640 will be solved.
*/
struct MaxTreeSizeRAII {
   Long64_t fOldMaxTreeSize;

   MaxTreeSizeRAII() : fOldMaxTreeSize(TTree::GetMaxTreeSize())
   {
      TTree::SetMaxTreeSize(std::numeric_limits<Long64_t>::max());
   }

   ~MaxTreeSizeRAII() { TTree::SetMaxTreeSize(fOldMaxTreeSize); }
};

struct DatasetLogInfo {
   std::string fDataSet;
   ULong64_t fRangeStart;
   ULong64_t fRangeEnd;
   unsigned int fSlot;
};

std::string LogRangeProcessing(const DatasetLogInfo &info)
{
   std::stringstream msg;
   msg << "Processing " << info.fDataSet << ": entry range [" << info.fRangeStart << "," << info.fRangeEnd - 1
       << "], using slot " << info.fSlot << " in thread " << std::this_thread::get_id() << '.';
   return msg.str();
}

DatasetLogInfo TreeDatasetLogInfo(const TTreeReader &r, unsigned int slot)
{
   const auto tree = r.GetTree();
   const auto chain = dynamic_cast<TChain *>(tree);
   std::string what;
   if (chain) {
      auto files = chain->GetListOfFiles();
      std::vector<std::string> treeNames;
      std::vector<std::string> fileNames;
      for (TObject *f : *files) {
         treeNames.emplace_back(f->GetName());
         fileNames.emplace_back(f->GetTitle());
      }
      what = "trees {";
      for (const auto &t : treeNames) {
         what += t + ",";
      }
      what.back() = '}';
      what += " in files {";
      for (const auto &f : fileNames) {
         what += f + ",";
      }
      what.back() = '}';
   } else {
      const auto treeName = tree->GetName();
      what = std::string("tree \"") + treeName + "\"";
      const auto file = tree->GetCurrentFile();
      if (file)
         what += std::string(" in file \"") + file->GetName() + "\"";
   }
   const auto entryRange = r.GetEntriesRange();
   const ULong64_t end = entryRange.second == -1ll ? tree->GetEntries() : entryRange.second;
   return {std::move(what), static_cast<ULong64_t>(entryRange.first), end, slot};
}

} // anonymous namespace

namespace ROOT {
namespace Detail {
namespace RDF {

/// A RAII object that calls RLoopManager::CleanUpTask at destruction
struct RCallCleanUpTask {
   RLoopManager &fLoopManager;
   unsigned int fArg;
   TTreeReader *fReader;

   RCallCleanUpTask(RLoopManager &lm, unsigned int arg = 0u, TTreeReader *reader = nullptr)
      : fLoopManager(lm), fArg(arg), fReader(reader)
   {
   }
   ~RCallCleanUpTask() { fLoopManager.CleanUpTask(fReader, fArg); }
};

} // namespace RDF
} // namespace Detail
} // namespace ROOT

///////////////////////////////////////////////////////////////////////////////
/// Get all the branches names, including the ones of the friend trees
ColumnNames_t ROOT::Internal::RDF::GetBranchNames(TTree &t, bool allowDuplicates)
{
   std::set<std::string> bNamesSet;
   ColumnNames_t bNames;
   std::set<TTree *> analysedTrees;
   std::string emptyFrName = "";
   GetBranchNamesImpl(t, bNamesSet, bNames, analysedTrees, emptyFrName, allowDuplicates);
   return bNames;
}

RLoopManager::RLoopManager(TTree *tree, const ColumnNames_t &defaultBranches)
   : fTree(std::shared_ptr<TTree>(tree, [](TTree *) {})), fDefaultColumns(defaultBranches),
     fNSlots(RDFInternal::GetNSlots()),
     fLoopType(ROOT::IsImplicitMTEnabled() ? ELoopType::kROOTFilesMT : ELoopType::kROOTFiles),
     fDataBlockNotifier(fNSlots)
{
}

RLoopManager::RLoopManager(ULong64_t nEmptyEntries)
   : fNEmptyEntries(nEmptyEntries), fNSlots(RDFInternal::GetNSlots()),
     fLoopType(ROOT::IsImplicitMTEnabled() ? ELoopType::kNoFilesMT : ELoopType::kNoFiles), fDataBlockNotifier(fNSlots)
{
}

RLoopManager::RLoopManager(std::unique_ptr<RDataSource> ds, const ColumnNames_t &defaultBranches)
   : fDefaultColumns(defaultBranches), fNSlots(RDFInternal::GetNSlots()),
     fLoopType(ROOT::IsImplicitMTEnabled() ? ELoopType::kDataSourceMT : ELoopType::kDataSource),
     fDataSource(std::move(ds)), fDataBlockNotifier(fNSlots)
{
   fDataSource->SetNSlots(fNSlots);
}

struct RSlotRAII {
   RSlotStack &fSlotStack;
   unsigned int fSlot;
   RSlotRAII(RSlotStack &slotStack) : fSlotStack(slotStack), fSlot(slotStack.GetSlot()) {}
   ~RSlotRAII() { fSlotStack.ReturnSlot(fSlot); }
};

/// Run event loop with no source files, in parallel.
void RLoopManager::RunEmptySourceMT()
{
#ifdef R__USE_IMT
   RSlotStack slotStack(fNSlots);
   // Working with an empty tree.
   // Evenly partition the entries according to fNSlots. Produce around 2 tasks per slot.
   const auto nEntriesPerSlot = fNEmptyEntries / (fNSlots * 2);
   auto remainder = fNEmptyEntries % (fNSlots * 2);
   std::vector<std::pair<ULong64_t, ULong64_t>> entryRanges;
   ULong64_t start = 0;
   while (start < fNEmptyEntries) {
      ULong64_t end = start + nEntriesPerSlot;
      if (remainder > 0) {
         ++end;
         --remainder;
      }
      entryRanges.emplace_back(start, end);
      start = end;
   }

   // Each task will generate a subrange of entries
   auto genFunction = [this, &slotStack](const std::pair<ULong64_t, ULong64_t> &range) {
      RSlotRAII slotRAII(slotStack);
      auto slot = slotRAII.fSlot;
      RCallCleanUpTask cleanup(*this, slot);
      InitNodeSlots(nullptr, slot);
      R__LOG_INFO(RDFLogChannel()) << LogRangeProcessing({"an empty source", range.first, range.second, slot});
      try {
         for (auto currEntry = range.first; currEntry < range.second; ++currEntry) {
            RunAndCheckFilters(slot, currEntry);
         }
      } catch (...) {
         // Error might throw in experiment frameworks like CMSSW
         std::cerr << "RDataFrame::Run: event loop was interrupted\n";
         throw;
      }
   };

   ROOT::TThreadExecutor pool;
   pool.Foreach(genFunction, entryRanges);

#endif // not implemented otherwise
}

/// Run event loop with no source files, in sequence.
void RLoopManager::RunEmptySource()
{
   InitNodeSlots(nullptr, 0);
   R__LOG_INFO(RDFLogChannel()) << LogRangeProcessing({"an empty source", 0, fNEmptyEntries, 0u});
   RCallCleanUpTask cleanup(*this);
   try {
      for (ULong64_t currEntry = 0; currEntry < fNEmptyEntries && fNStopsReceived < fNChildren; ++currEntry) {
         RunAndCheckFilters(0, currEntry);
      }
   } catch (...) {
      std::cerr << "RDataFrame::Run: event loop was interrupted\n";
      throw;
   }
}

/// Run event loop over one or multiple ROOT files, in parallel.
void RLoopManager::RunTreeProcessorMT()
{
#ifdef R__USE_IMT
   RSlotStack slotStack(fNSlots);
   const auto &entryList = fTree->GetEntryList() ? *fTree->GetEntryList() : TEntryList();
   auto tp = std::make_unique<ROOT::TTreeProcessorMT>(*fTree, entryList, fNSlots);

   std::atomic<ULong64_t> entryCount(0ull);

   tp->Process([this, &slotStack, &entryCount](TTreeReader &r) -> void {
      RSlotRAII slotRAII(slotStack);
      auto slot = slotRAII.fSlot;
      RCallCleanUpTask cleanup(*this, slot, &r);
      InitNodeSlots(&r, slot);
      R__LOG_INFO(RDFLogChannel()) << LogRangeProcessing(TreeDatasetLogInfo(r, slot));
      const auto entryRange = r.GetEntriesRange(); // we trust TTreeProcessorMT to call SetEntriesRange
      const auto nEntries = entryRange.second - entryRange.first;
      auto count = entryCount.fetch_add(nEntries);
      try {
         // recursive call to check filters and conditionally execute actions
         while (r.Next()) {
            RunAndCheckFilters(slot, count++);
         }
      } catch (...) {
         std::cerr << "RDataFrame::Run: event loop was interrupted\n";
         throw;
      }
   });
#endif // no-op otherwise (will not be called)
}

/// Run event loop over one or multiple ROOT files, in sequence.
void RLoopManager::RunTreeReader()
{
   TTreeReader r(fTree.get(), fTree->GetEntryList());
   if (0 == fTree->GetEntriesFast())
      return;
   RCallCleanUpTask cleanup(*this, 0u, &r);
   InitNodeSlots(&r, 0);
   R__LOG_INFO(RDFLogChannel()) << LogRangeProcessing(TreeDatasetLogInfo(r, 0u));

   // recursive call to check filters and conditionally execute actions
   // in the non-MT case processing can be stopped early by ranges, hence the check on fNStopsReceived
   try {
      while (r.Next() && fNStopsReceived < fNChildren) {
         RunAndCheckFilters(0, r.GetCurrentEntry());
      }
   } catch (...) {
      std::cerr << "RDataFrame::Run: event loop was interrupted\n";
      throw;
   }
   if (r.GetEntryStatus() != TTreeReader::kEntryNotFound && fNStopsReceived < fNChildren) {
      // something went wrong in the TTreeReader event loop
      throw std::runtime_error("An error was encountered while processing the data. TTreeReader status code is: " +
                               std::to_string(r.GetEntryStatus()));
   }
}

/// Run event loop over data accessed through a DataSource, in sequence.
void RLoopManager::RunDataSource()
{
   R__DEBUG_ASSERT(fDataSource != nullptr);
   fDataSource->Initialise();
   auto ranges = fDataSource->GetEntryRanges();
   while (!ranges.empty() && fNStopsReceived < fNChildren) {
      InitNodeSlots(nullptr, 0u);
      fDataSource->InitSlot(0u, 0ull);
      RCallCleanUpTask cleanup(*this);
      try {
         for (const auto &range : ranges) {
            const auto start = range.first;
            const auto end = range.second;
            R__LOG_INFO(RDFLogChannel()) << LogRangeProcessing({fDataSource->GetLabel(), start, end, 0u});
            for (auto entry = start; entry < end && fNStopsReceived < fNChildren; ++entry) {
               if (fDataSource->SetEntry(0u, entry)) {
                  RunAndCheckFilters(0u, entry);
               }
            }
         }
      } catch (...) {
         std::cerr << "RDataFrame::Run: event loop was interrupted\n";
         throw;
      }
      fDataSource->FinaliseSlot(0u);
      ranges = fDataSource->GetEntryRanges();
   }
   fDataSource->Finalise();
}

/// Run event loop over data accessed through a DataSource, in parallel.
void RLoopManager::RunDataSourceMT()
{
#ifdef R__USE_IMT
   R__DEBUG_ASSERT(fDataSource != nullptr);
   RSlotStack slotStack(fNSlots);
   ROOT::TThreadExecutor pool;

   // Each task works on a subrange of entries
   auto runOnRange = [this, &slotStack](const std::pair<ULong64_t, ULong64_t> &range) {
      RSlotRAII slotRAII(slotStack);
      const auto slot = slotRAII.fSlot;
      InitNodeSlots(nullptr, slot);
      RCallCleanUpTask cleanup(*this, slot);
      fDataSource->InitSlot(slot, range.first);
      const auto start = range.first;
      const auto end = range.second;
      R__LOG_INFO(RDFLogChannel()) << LogRangeProcessing({fDataSource->GetLabel(), start, end, slot});
      try {
         for (auto entry = start; entry < end; ++entry) {
            if (fDataSource->SetEntry(slot, entry)) {
               RunAndCheckFilters(slot, entry);
            }
         }
      } catch (...) {
         std::cerr << "RDataFrame::Run: event loop was interrupted\n";
         throw;
      }
      fDataSource->FinaliseSlot(slot);
   };

   fDataSource->Initialise();
   auto ranges = fDataSource->GetEntryRanges();
   while (!ranges.empty()) {
      pool.Foreach(runOnRange, ranges);
      ranges = fDataSource->GetEntryRanges();
   }
   fDataSource->Finalise();
#endif // not implemented otherwise (never called)
}

/// Execute actions and make sure named filters are called for each event.
/// Named filters must be called even if the analysis logic would not require it, lest they report confusing results.
void RLoopManager::RunAndCheckFilters(unsigned int slot, Long64_t entry)
{
   // data-block callbacks run before the rest of the graph
   if (fDataBlockNotifier.CheckFlag(slot)) {
      for (auto &callback : fDataBlockCallbacks) {
         callback(slot);
      }
      fDataBlockNotifier.UnsetFlag(slot);
   }

   for (auto &actionPtr : fBookedActions)
      actionPtr->Run(slot, entry);
   for (auto &namedFilterPtr : fBookedNamedFilters)
      namedFilterPtr->CheckFilters(slot, entry);
   for (auto &callback : fCallbacks)
      callback(slot);
}

/// Build TTreeReaderValues for all nodes
/// This method loops over all filters, actions and other booked objects and
/// calls their `InitSlot` method, to get them ready for running a task.
void RLoopManager::InitNodeSlots(TTreeReader *r, unsigned int slot)
{
   SetupDataBlockCallbacks(r, slot);
   for (auto &ptr : fBookedActions)
      ptr->InitSlot(r, slot);
   for (auto &ptr : fBookedFilters)
      ptr->InitSlot(r, slot);
   for (auto &callback : fCallbacksOnce)
      callback(slot);
}

void RLoopManager::SetupDataBlockCallbacks(TTreeReader *r, unsigned int slot) {
   if (r != nullptr) {
      // we need to set a notifier so that we run the callbacks every time we switch to a new TTree
      // `PrependLink` inserts this notifier into the TTree/TChain's linked list of notifiers
      fDataBlockNotifier.GetChainNotifyLink(slot).PrependLink(*r->GetTree());
   }
   // Whatever the data source, initially set the "new data block" flag:
   // - for TChains, this ensures that we don't skip the first data block because
   //   the correct tree is already loaded
   // - for RDataSources and empty sources, which currently don't have data blocks, this
   //   ensures that we run once per task
   fDataBlockNotifier.SetFlag(slot);
}

/// Initialize all nodes of the functional graph before running the event loop.
/// This method is called once per event-loop and performs generic initialization
/// operations that do not depend on the specific processing slot (i.e. operations
/// that are common for all threads).
void RLoopManager::InitNodes()
{
   EvalChildrenCounts();
   for (auto &filter : fBookedFilters)
      filter->InitNode();
   for (auto &range : fBookedRanges)
      range->InitNode();
   for (auto &ptr : fBookedActions)
      ptr->Initialize();
}

/// Perform clean-up operations. To be called at the end of each event loop.
void RLoopManager::CleanUpNodes()
{
   fMustRunNamedFilters = false;

   // forget RActions and detach TResultProxies
   for (auto &ptr : fBookedActions)
      ptr->Finalize();

   fRunActions.insert(fRunActions.begin(), fBookedActions.begin(), fBookedActions.end());
   fBookedActions.clear();

   // reset children counts
   fNChildren = 0;
   fNStopsReceived = 0;
   for (auto &ptr : fBookedFilters)
      ptr->ResetChildrenCount();
   for (auto &ptr : fBookedRanges)
      ptr->ResetChildrenCount();

   fCallbacks.clear();
   fCallbacksOnce.clear();
   fDataBlockCallbacks.clear();
}

/// Perform clean-up operations. To be called at the end of each task execution.
void RLoopManager::CleanUpTask(TTreeReader *r, unsigned int slot)
{
   if (r != nullptr)
      fDataBlockNotifier.GetChainNotifyLink(slot).RemoveLink(*r->GetTree());
   for (auto &ptr : fBookedActions)
      ptr->FinalizeSlot(slot);
   for (auto &ptr : fBookedFilters)
      ptr->FinaliseSlot(slot);
}

/// Add RDF nodes that require just-in-time compilation to the computation graph.
/// This method also clears the contents of GetCodeToJit().
void RLoopManager::Jit()
{
   // TODO this should be a read lock unless we find GetCodeToJit non-empty
   R__LOCKGUARD(gROOTMutex);

   const std::string code = std::move(GetCodeToJit());
   if (code.empty()) {
      R__LOG_INFO(RDFLogChannel()) << "Nothing to jit and execute.";
      return;
   }

   TStopwatch s;
   s.Start();
   RDFInternal::InterpreterCalc(code, "RLoopManager::Run");
   s.Stop();
   R__LOG_INFO(RDFLogChannel()) << "Just-in-time compilation phase completed"
                                << (s.RealTime() > 1e-3 ? " in " + std::to_string(s.RealTime()) + " seconds." : ".");
}

/// Trigger counting of number of children nodes for each node of the functional graph.
/// This is done once before starting the event loop. Each action sends an `increase children count` signal
/// upstream, which is propagated until RLoopManager. Each time a node receives the signal, in increments its
/// children counter. Each node only propagates the signal once, even if it receives it multiple times.
/// Named filters also send an `increase children count` signal, just like actions, as they always execute during
/// the event loop so the graph branch they belong to must count as active even if it does not end in an action.
void RLoopManager::EvalChildrenCounts()
{
   for (auto &actionPtr : fBookedActions)
      actionPtr->TriggerChildrenCount();
   for (auto &namedFilterPtr : fBookedNamedFilters)
      namedFilterPtr->TriggerChildrenCount();
}

/// Start the event loop with a different mechanism depending on IMT/no IMT, data source/no data source.
/// Also perform a few setup and clean-up operations (jit actions if necessary, clear booked actions after the loop...).
void RLoopManager::Run()
{
   // Change value of TTree::GetMaxTreeSize only for this scope. Revert when #6640 will be solved.
   MaxTreeSizeRAII ctxtmts;

   R__LOG_INFO(RDFLogChannel()) << "Starting event loop number " << fNRuns << '.';

   ThrowIfNSlotsChanged(GetNSlots());

   Jit();

   InitNodes();

   TStopwatch s;
   s.Start();
   switch (fLoopType) {
   case ELoopType::kNoFilesMT: RunEmptySourceMT(); break;
   case ELoopType::kROOTFilesMT: RunTreeProcessorMT(); break;
   case ELoopType::kDataSourceMT: RunDataSourceMT(); break;
   case ELoopType::kNoFiles: RunEmptySource(); break;
   case ELoopType::kROOTFiles: RunTreeReader(); break;
   case ELoopType::kDataSource: RunDataSource(); break;
   }
   s.Stop();

   CleanUpNodes();

   fNRuns++;

   R__LOG_INFO(RDFLogChannel()) << "Finished event loop number " << fNRuns - 1 << " (" << s.CpuTime() << "s CPU, "
                                << s.RealTime() << "s elapsed).";
}

/// Return the list of default columns -- empty if none was provided when constructing the RDataFrame
const ColumnNames_t &RLoopManager::GetDefaultColumnNames() const
{
   return fDefaultColumns;
}

TTree *RLoopManager::GetTree() const
{
   return fTree.get();
}

void RLoopManager::Book(RDFInternal::RActionBase *actionPtr)
{
   fBookedActions.emplace_back(actionPtr);
}

void RLoopManager::Deregister(RDFInternal::RActionBase *actionPtr)
{
   RDFInternal::Erase(actionPtr, fRunActions);
   RDFInternal::Erase(actionPtr, fBookedActions);
}

void RLoopManager::Book(RFilterBase *filterPtr)
{
   fBookedFilters.emplace_back(filterPtr);
   if (filterPtr->HasName()) {
      fBookedNamedFilters.emplace_back(filterPtr);
      fMustRunNamedFilters = true;
   }
}

void RLoopManager::Deregister(RFilterBase *filterPtr)
{
   RDFInternal::Erase(filterPtr, fBookedFilters);
   RDFInternal::Erase(filterPtr, fBookedNamedFilters);
}

void RLoopManager::Book(RRangeBase *rangePtr)
{
   fBookedRanges.emplace_back(rangePtr);
}

void RLoopManager::Deregister(RRangeBase *rangePtr)
{
   RDFInternal::Erase(rangePtr, fBookedRanges);
}

// dummy call, end of recursive chain of calls
bool RLoopManager::CheckFilters(unsigned int, Long64_t)
{
   return true;
}

/// Call `FillReport` on all booked filters
void RLoopManager::Report(ROOT::RDF::RCutFlowReport &rep) const
{
   for (const auto &fPtr : fBookedNamedFilters)
      fPtr->FillReport(rep);
}

void RLoopManager::ToJitExec(const std::string &code) const
{
   R__LOCKGUARD(gROOTMutex);
   GetCodeToJit().append(code);
}

void RLoopManager::RegisterCallback(ULong64_t everyNEvents, std::function<void(unsigned int)> &&f)
{
   if (everyNEvents == 0ull)
      fCallbacksOnce.emplace_back(std::move(f), fNSlots);
   else
      fCallbacks.emplace_back(everyNEvents, std::move(f), fNSlots);
}

std::vector<std::string> RLoopManager::GetFiltersNames()
{
   std::vector<std::string> filters;
   for (auto &filter : fBookedFilters) {
      auto name = (filter->HasName() ? filter->GetName() : "Unnamed Filter");
      filters.push_back(name);
   }
   return filters;
}

std::vector<RNodeBase *> RLoopManager::GetGraphEdges() const
{
   std::vector<RNodeBase *> nodes(fBookedFilters.size() + fBookedRanges.size());
   auto it = std::copy(fBookedFilters.begin(), fBookedFilters.end(), nodes.begin());
   std::copy(fBookedRanges.begin(), fBookedRanges.end(), it);
   return nodes;
}

std::vector<RDFInternal::RActionBase *> RLoopManager::GetAllActions() const
{
   std::vector<RDFInternal::RActionBase *> actions(fBookedActions.size() + fRunActions.size());
   auto it = std::copy(fBookedActions.begin(), fBookedActions.end(), actions.begin());
   std::copy(fRunActions.begin(), fRunActions.end(), it);
   return actions;
}

std::shared_ptr<ROOT::Internal::RDF::GraphDrawing::GraphNode> RLoopManager::GetGraph()
{
   std::string name;
   if (fDataSource) {
      name = fDataSource->GetLabel();
   } else if (fTree) {
      name = fTree->GetName();
   } else {
      name = std::to_string(fNEmptyEntries);
   }

   auto thisNode = std::make_shared<ROOT::Internal::RDF::GraphDrawing::GraphNode>(name);
   thisNode->SetRoot();
   thisNode->SetCounter(0);
   return thisNode;
}

////////////////////////////////////////////////////////////////////////////
/// Return all valid TTree::Branch names (caching results for subsequent calls).
/// Never use fBranchNames directy, always request it through this method.
const ColumnNames_t &RLoopManager::GetBranchNames()
{
   if (fValidBranchNames.empty() && fTree) {
      fValidBranchNames = RDFInternal::GetBranchNames(*fTree, /*allowRepetitions=*/true);
   }
   return fValidBranchNames;
}

bool RLoopManager::HasDSValuePtrs(const std::string &col) const
{
   return fDSValuePtrMap.find(col) != fDSValuePtrMap.end();
}

void RLoopManager::AddDSValuePtrs(const std::string &col, const std::vector<void *> ptrs)
{
   fDSValuePtrMap[col] = ptrs;
}

void RLoopManager::AddDataBlockCallback(std::function<void(unsigned int)> &&callback)
{
   if (callback)
      fDataBlockCallbacks.emplace_back(std::move(callback));
}
