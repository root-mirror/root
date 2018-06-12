// Author: Enrico Guiraud, Danilo Piparo CERN  03/2017

/*************************************************************************
 * Copyright (C) 1995-2016, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_RDFNODES
#define ROOT_RDFNODES

#include "ROOT/RCutFlowReport.hxx"
#include "ROOT/RDataSource.hxx"
#include "ROOT/RDFNodesUtils.hxx"
#include "ROOT/RDFBookedCustomColumns.hxx"
#include "ROOT/RDFUtils.hxx"
#include "ROOT/GraphNode.hxx"
#include "ROOT/RIntegerSequence.hxx"
#include "ROOT/RMakeUnique.hxx"
#include "ROOT/RVec.hxx"
#include "ROOT/TRWSpinLock.hxx"
#include "ROOT/TypeTraits.hxx"
#include "TError.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"

#include <deque> // std::vector substitute in case of vector<bool>
#include <functional>
#include <limits>
#include <map>
#include <numeric> // std::accumulate (FillReport), std::iota (TSlotStack)
#include <stack>
#include <string>
#include <thread>
#include <tuple>
#include <type_traits>
#include <vector>

namespace ROOT {
namespace RDF{
template <typename Proxied, typename DataSource>
class RInterface;
template <typename T>
class RResultPtr;

}// Namespace RDF

namespace Internal {
namespace RDF {
class RActionBase;
class GraphCreatorHelper;

// This is an helper class to allow to pick a slot resorting to a map
// indexed by thread ids.
// WARNING: this class does not work as a regular stack. The size is
// fixed at construction time and no blocking is foreseen.
class TSlotStack {
private:
   unsigned int &GetCount();
   unsigned int &GetIndex();
   unsigned int fCursor;
   std::vector<unsigned int> fBuf;
   ROOT::TRWSpinLock fRWLock;

public:
   TSlotStack() = delete;
   TSlotStack(unsigned int size) : fCursor(size), fBuf(size) { std::iota(fBuf.begin(), fBuf.end(), 0U); }
   void ReturnSlot(unsigned int slotNumber);
   unsigned int GetSlot();
   std::map<std::thread::id, unsigned int> fCountMap;
   std::map<std::thread::id, unsigned int> fIndexMap;
};
} // namespace RDF
} // namespace Internal



namespace Detail {
namespace RDF {
class RCustomColumnBase;

using namespace ROOT::TypeTraits;
namespace RDFInternal = ROOT::Internal::RDF;

// forward declarations for RLoopManager
class RCustomColumnBase;
class RFilterBase;
class RRangeBase;

class RLoopManager;

/// Base class for non-leaf nodes of the computational graph.
/// It only exposes the bare minimum interface required to work as a generic part of the computation graph.
/// RDataFrames and results of transformations can be cast to this type via ROOT::RDF::ToCommonNodeType.
class RNode {
protected:
   RLoopManager *fLoopManager;
   unsigned int fNChildren{0};      ///< Number of nodes of the functional graph hanging from this object
   unsigned int fNStopsReceived{0}; ///< Number of times that a children node signaled to stop processing entries.

public:
   RNode(RLoopManager *lm = nullptr) : fLoopManager(lm) {}
   virtual ~RNode() {}
   virtual bool CheckFilters(unsigned int, Long64_t) = 0;
   virtual void Report(ROOT::RDF::RCutFlowReport &) const = 0;
   virtual void PartialReport(ROOT::RDF::RCutFlowReport &) const = 0;
   virtual void IncrChildrenCount() = 0;
   virtual void StopProcessing() = 0;
   virtual void AddFilterName(std::vector<std::string> &filters) = 0;

   virtual void ResetChildrenCount()
   {
      fNChildren = 0;
      fNStopsReceived = 0;
   }

   virtual RLoopManager *GetLoopManagerUnchecked() { return fLoopManager; }
};

class RLoopManager : public RNode {
   friend class ROOT::Internal::RDF::GraphDrawing::GraphCreatorHelper;
   using RDataSource = ROOT::RDF::RDataSource;
   enum class ELoopType { kROOTFiles, kROOTFilesMT, kNoFiles, kNoFilesMT, kDataSource, kDataSourceMT };
   using Callback_t = std::function<void(unsigned int)>;
   class TCallback {
      const Callback_t fFun;
      const ULong64_t fEveryN;
      std::vector<ULong64_t> fCounters;

   public:
      TCallback(ULong64_t everyN, Callback_t &&f, unsigned int nSlots)
         : fFun(std::move(f)), fEveryN(everyN), fCounters(nSlots, 0ull)
      {
      }

      void operator()(unsigned int slot)
      {
         auto &c = fCounters[slot];
         ++c;
         if (c == fEveryN) {
            c = 0ull;
            fFun(slot);
         }
      }
   };

   class TOneTimeCallback {
      const Callback_t fFun;
      std::vector<int> fHasBeenCalled; // std::vector<bool> is thread-unsafe for our purposes (and generally evil)

   public:
      TOneTimeCallback(Callback_t &&f, unsigned int nSlots) : fFun(std::move(f)), fHasBeenCalled(nSlots, 0) {}

      void operator()(unsigned int slot)
      {
         if (fHasBeenCalled[slot] == 1)
            return;
         fFun(slot);
         fHasBeenCalled[slot] = 1;
      }
   };

   std::vector<RDFInternal::RActionBase *> fBookedActions; ///< Non-owning pointers to actions to be run
   std::vector<RDFInternal::RActionBase *> fRunActions; ///< Non-owning pointers to actions already run
   std::vector<RFilterBase *> fBookedFilters;
   std::vector<RFilterBase *> fBookedNamedFilters; ///< Contains a subset of fBookedFilters, i.e. only the named filters
   std::vector<RRangeBase *> fBookedRanges;

   /// Shared pointer to the input TTree. It does not delete the pointee if the TTree/TChain was passed directly as an
   /// argument to RDataFrame's ctor (in which case we let users retain ownership).
   std::shared_ptr<TTree> fTree{nullptr};
   const ColumnNames_t fDefaultColumns;
   const ULong64_t fNEmptyEntries{0};
   const unsigned int fNSlots{1};
   bool fMustRunNamedFilters{true};
   const ELoopType fLoopType; ///< The kind of event loop that is going to be run (e.g. on ROOT files, on no files)
   std::string fToJit;        ///< code that should be jitted and executed right before the event loop
   const std::unique_ptr<RDataSource> fDataSource; ///< Owning pointer to a data-source object. Null if no data-source
   std::map<std::string, std::string> fAliasColumnNameMap; ///< ColumnNameAlias-columnName pairs
   std::vector<TCallback> fCallbacks;                      ///< Registered callbacks
   std::vector<TOneTimeCallback> fCallbacksOnce; ///< Registered callbacks to invoke just once before running the loop
   /// A unique ID that identifies the computation graph that starts with this RLoopManager.
   /// Used, for example, to jit objects in a namespace reserved for this computation graph
   const unsigned int fID = GetNextID();

   std::vector<RCustomColumnBase* > fCustomColumns; ///< The loopmanager tracks all columns created, without owning them.

   void RunEmptySourceMT();
   void RunEmptySource();
   void RunTreeProcessorMT();
   void RunTreeReader();
   void RunDataSourceMT();
   void RunDataSource();
   void RunAndCheckFilters(unsigned int slot, Long64_t entry);
   void InitNodeSlots(TTreeReader *r, unsigned int slot);
   void InitNodes();
   void CleanUpNodes();
   void CleanUpTask(unsigned int slot);
   void EvalChildrenCounts();
   unsigned int GetNextID() const;

public:
   RLoopManager(TTree *tree, const ColumnNames_t &defaultBranches);
   RLoopManager(ULong64_t nEmptyEntries);
   RLoopManager(std::unique_ptr<RDataSource> ds, const ColumnNames_t &defaultBranches);
   RLoopManager(const RLoopManager &) = delete;
   RLoopManager &operator=(const RLoopManager &) = delete;

   void BuildJittedNodes();
   RLoopManager *GetLoopManagerUnchecked() final { return this; }
   void Run();
   const ColumnNames_t &GetDefaultColumnNames() const;
   TTree *GetTree() const;
   ::TDirectory *GetDirectory() const;
   ULong64_t GetNEmptyEntries() const { return fNEmptyEntries; }
   RDataSource *GetDataSource() const { return fDataSource.get(); }
   void Book(RDFInternal::RActionBase *actionPtr);
   void Deregister(RDFInternal::RActionBase *actionPtr);
   void Book(RFilterBase *filterPtr);
   void Deregister(RFilterBase *filterPtr);
   void Book(RRangeBase *rangePtr);
   void Deregister(RRangeBase *rangePtr);
   bool CheckFilters(unsigned int, Long64_t) final;
   unsigned int GetNSlots() const { return fNSlots; }
   bool MustRunNamedFilters() const { return fMustRunNamedFilters; }
   void Report(ROOT::RDF::RCutFlowReport &rep) const final;
   /// End of recursive chain of calls, does nothing
   void PartialReport(ROOT::RDF::RCutFlowReport &) const final {}
   void SetTree(const std::shared_ptr<TTree> &tree) { fTree = tree; }
   void IncrChildrenCount() final { ++fNChildren; }
   void StopProcessing() final { ++fNStopsReceived; }
   void ToJit(const std::string &s) { fToJit.append(s); }
   void AddColumnAlias(const std::string &alias, const std::string &colName) { fAliasColumnNameMap[alias] = colName; }
   const std::map<std::string, std::string> &GetAliasMap() const { return fAliasColumnNameMap; }
   void RegisterCallback(ULong64_t everyNEvents, std::function<void(unsigned int)> &&f);
   unsigned int GetID() const { return fID; }

   /// End of recursive chain of calls, does nothing
   void AddFilterName(std::vector<std::string> &) {}
   /// For each booked filter, returns either the name or "Unnamed Filter"
   std::vector<std::string> GetFiltersNames();

   /// For all the actions, either booked or run
   std::vector<RDFInternal::RActionBase *> GetAllActions();

   void RegisterCustomColumn(RCustomColumnBase *column){
      fCustomColumns.push_back(column);
   }

   void DeRegisterCustomColumn(RCustomColumnBase *column){
      fCustomColumns.erase(std::remove(fCustomColumns.begin(), fCustomColumns.end(), column), fCustomColumns.end());
   }

   std::vector<RDFInternal::RActionBase *> GetBookedActions(){
      return fBookedActions;
   }
   std::shared_ptr<ROOT::Internal::RDF::GraphDrawing::GraphNode> GetGraph();

};

} //namespace RDF
} //namespace Detail

namespace Internal {
namespace RDF {
namespace GraphDrawing {
// Forward declarations for all nodes. Putting them here because RFilter, RRange, and RCustomColumn have been already
// declared.
std::shared_ptr<GraphNode>
CreateDefineNode(const std::string &columnName, const ROOT::Detail::RDF::RCustomColumnBase *columnPtr);

std::shared_ptr<GraphNode> CreateFilterNode(const ROOT::Detail::RDF::RFilterBase *filterPtr);

std::shared_ptr<GraphNode> CreateRangeNode(const ROOT::Detail::RDF::RRangeBase *rangePtr);

bool CheckIfDefaultOrDSColumn(const std::string &name,
                              const std::shared_ptr<ROOT::Detail::RDF::RCustomColumnBase> &column);
} // namespace GraphDrawing
} // namespace RDF
} // namespace Internal

namespace Detail {
namespace RDF {
class RCustomColumnBase {
protected:
   RLoopManager *fLoopManager; ///< A raw pointer to the RLoopManager at the root of this functional graph. It is only
/// guaranteed to contain a valid address during an event loop.
   const std::string fName;
   unsigned int fNChildren{0};      ///< number of nodes of the functional graph hanging from this object
   unsigned int fNStopsReceived{0}; ///< number of times that a children node signaled to stop processing entries.
   const unsigned int fNSlots;      ///< number of thread slots used by this node, inherited from parent node.
   const bool fIsDataSourceColumn; ///< does the custom column refer to a data-source column? (or a user-define column?)
   std::vector<Long64_t> fLastCheckedEntry;

   RDFInternal::RBookedCustomColumns fCustomColumns;

public:
   RCustomColumnBase(RLoopManager *lm, std::string_view name, const unsigned int nSlots, const bool isDSColumn,
                     const RDFInternal::RBookedCustomColumns &customColumns);

   RCustomColumnBase &operator=(const RCustomColumnBase &) = delete;
   virtual ~RCustomColumnBase(); // outlined defaulted.
   virtual void InitSlot(TTreeReader *r, unsigned int slot) = 0;
   virtual void *GetValuePtr(unsigned int slot) = 0;
   virtual const std::type_info &GetTypeId() const = 0;
   RLoopManager *GetLoopManagerUnchecked() const;
   std::string GetName() const;
   virtual void Update(unsigned int slot, Long64_t entry) = 0;
   virtual void ClearValueReaders(unsigned int slot) = 0;
   bool IsDataSourceColumn() const { return fIsDataSourceColumn; }
   virtual void InitNode();
};

/// A wrapper around a concrete RCustomColumn, which forwards all calls to it
/// RJittedCustomColumn is a placeholder that is put in the collection of custom columns in place of a RCustomColumn
/// that will be just-in-time compiled. Jitted code will assign the concrete RCustomColumn to this RJittedCustomColumn
/// before the event-loop starts.
class RJittedCustomColumn : public RCustomColumnBase {
   std::unique_ptr<RCustomColumnBase> fConcreteCustomColumn = nullptr;

public:
   RJittedCustomColumn(RLoopManager *lm, std::string_view name, unsigned int nSlots)
      : RCustomColumnBase(lm, name, nSlots, /*isDSColumn=*/false, RDFInternal::RBookedCustomColumns())
   {
   }

   void SetCustomColumn(std::unique_ptr<RCustomColumnBase> c) { fConcreteCustomColumn = std::move(c); }

   void InitSlot(TTreeReader *r, unsigned int slot) final;
   void *GetValuePtr(unsigned int slot) final;
   const std::type_info &GetTypeId() const final;
   void Update(unsigned int slot, Long64_t entry) final;
   void ClearValueReaders(unsigned int slot) final;
   void InitNode() final;
};

// clang-format off
namespace CustomColExtraArgs {
struct None{};
struct Slot{};
struct SlotAndEntry{};
}
// clang-format on

} // namespace RDF
} // namespace Detail

namespace Internal {
namespace RDF {
using namespace ROOT::Detail::RDF;
namespace RDFGraphDrawing = ROOT::Internal::RDF::GraphDrawing;

/**
\class ROOT::Internal::RDF::TColumnValue
\ingroup dataframe
\brief Helper class that updates and returns TTree branches as well as RDataFrame temporary columns
\tparam T The type of the column

RDataFrame nodes must access two different types of values during the event loop:
values of real branches, for which TTreeReader{Values,Arrays} act as proxies, or
temporary columns whose values are generated on the fly. While the type of the
value is known at compile time (or just-in-time), it is only at runtime that nodes
can check whether a certain value is generated on the fly or not.

TColumnValue abstracts this difference by providing the same interface for
both cases and handling the reading or generation of new values transparently.
Only one of the two data members fReaderProxy or fValuePtr will be non-null
for a given TColumnValue, depending on whether the value comes from a real
TTree branch or from a temporary column respectively.

RDataFrame nodes can store tuples of TColumnValues and retrieve an updated
value for the column via the `Get` method.
**/
template <typename T>
class TColumnValue {

   using MustUseRVec_t = IsRVec_t<T>;

   // ColumnValue_t is the type of the column or the type of the elements of an array column
   using ColumnValue_t = typename std::conditional<MustUseRVec_t::value, TakeFirstParameter_t<T>, T>::type;
   using TreeReader_t = typename std::conditional<MustUseRVec_t::value, TTreeReaderArray<ColumnValue_t>,
                                                  TTreeReaderValue<ColumnValue_t>>::type;

   /// TColumnValue has a slightly different behaviour whether the column comes from a TTreeReader, a RDataFrame Define
   /// or a RDataSource. It stores which it is as an enum.
   enum class EColumnKind { kTree, kCustomColumn, kDataSource, kInvalid };
   // Set to the correct value by MakeProxy or SetTmpColumn
   EColumnKind fColumnKind = EColumnKind::kInvalid;
   /// The slot this value belongs to. Only needed when querying custom column values, it is set in `SetTmpColumn`.
   unsigned int fSlot = std::numeric_limits<unsigned int>::max();

   // Each element of the following stacks will be in use by a _single task_.
   // Each task will push one element when it starts and pop it when it ends.
   // Stacks will typically be very small (1-2 elements typically) and will only grow over size 1 in case of interleaved
   // task execution i.e. when more than one task needs readers in this worker thread.

   /// Owning ptrs to a TTreeReaderValue or TTreeReaderArray. Only used for Tree columns.
   std::stack<std::unique_ptr<TreeReader_t>> fTreeReaders;
   /// Non-owning ptrs to the value of a custom column.
   std::stack<T *> fCustomValuePtrs;
   /// Non-owning ptrs to the value of a data-source column.
   std::stack<T **> fDSValuePtrs;
   /// Non-owning ptrs to the node responsible for the custom column. Needed when querying custom values.
   std::stack<RCustomColumnBase *> fCustomColumns;
   /// Enumerator for the different properties of the branch storage in memory
   enum class EStorageType : char { kContiguous, kUnknown, kSparse };
   /// Signal whether we ever checked that the branch we are reading with a TTreeReaderArray stores array elements
   /// in contiguous memory. Only used when T == RVec<U>.
   EStorageType fStorageType = EStorageType::kUnknown;
   /// If MustUseRVec, i.e. we are reading an array, we return a reference to this RVec to clients
   RVec<ColumnValue_t> fRVec;
   bool fCopyWarningPrinted = false;

public:
   TColumnValue(){};

   void SetTmpColumn(unsigned int slot, RCustomColumnBase *tmpColumn);

   void MakeProxy(TTreeReader *r, const std::string &bn)
   {
      fColumnKind = EColumnKind::kTree;
      fTreeReaders.emplace(std::make_unique<TreeReader_t>(*r, bn.c_str()));
   }

   /// This overload is used to return scalar quantities (i.e. types that are not read into a RVec)
   template <typename U = T, typename std::enable_if<!TColumnValue<U>::MustUseRVec_t::value, int>::type = 0>
   T &Get(Long64_t entry);

   /// This overload is used to return arrays (i.e. types that are read into a RVec).
   /// In this case the returned T is always a RVec<ColumnValue_t>.
   template <typename U = T, typename std::enable_if<TColumnValue<U>::MustUseRVec_t::value, int>::type = 0>
   T &Get(Long64_t entry);

   void Reset()
   {
      switch (fColumnKind) {
      case EColumnKind::kTree: fTreeReaders.pop(); break;
      case EColumnKind::kCustomColumn:
         fCustomColumns.pop();
         fCustomValuePtrs.pop();
         break;
      case EColumnKind::kDataSource:
         fCustomColumns.pop();
         fDSValuePtrs.pop();
         break;
      case EColumnKind::kInvalid: throw std::runtime_error("ColumnKind not set for this TColumnValue");
      }
   }
};

// Some extern instaniations to speed-up compilation/interpretation time
// These are not active if c++17 is enabled because of a bug in our clang
// See ROOT-9499.
#if __cplusplus < 201703L
extern template class TColumnValue<int>;
extern template class TColumnValue<unsigned int>;
extern template class TColumnValue<char>;
extern template class TColumnValue<unsigned char>;
extern template class TColumnValue<float>;
extern template class TColumnValue<double>;
extern template class TColumnValue<Long64_t>;
extern template class TColumnValue<ULong64_t>;
extern template class TColumnValue<std::vector<int>>;
extern template class TColumnValue<std::vector<unsigned int>>;
extern template class TColumnValue<std::vector<char>>;
extern template class TColumnValue<std::vector<unsigned char>>;
extern template class TColumnValue<std::vector<float>>;
extern template class TColumnValue<std::vector<double>>;
extern template class TColumnValue<std::vector<Long64_t>>;
extern template class TColumnValue<std::vector<ULong64_t>>;
#endif

template <typename T>
struct TRDFValueTuple {
};

template <typename... BranchTypes>
struct TRDFValueTuple<TypeList<BranchTypes...>> {
   using type = std::tuple<TColumnValue<BranchTypes>...>;
};

template <typename BranchType>
using RDFValueTuple_t = typename TRDFValueTuple<BranchType>::type;

/// Clear the proxies of a tuple of TColumnValues
template <typename ValueTuple, std::size_t... S>
void ResetRDFValueTuple(ValueTuple &values, std::index_sequence<S...>)
{
   // hack to expand a parameter pack without c++17 fold expressions.
   std::initializer_list<int> expander{(std::get<S>(values).Reset(), 0)...};
   (void)expander; // avoid "unused variable" warnings
}

class RActionBase {
protected:
   RLoopManager *fLoopManager; ///< A raw pointer to the RLoopManager at the root of this functional
                               /// graph. It is only guaranteed to contain a valid address during an
                               /// event loop.
   const unsigned int fNSlots; ///< Number of thread slots used by this node.

   RBookedCustomColumns fCustomColumns;

public:
   RActionBase(RLoopManager *implPtr, const unsigned int nSlots, const RBookedCustomColumns &customColumns);

   RActionBase(const RActionBase &) = delete;
   RActionBase &operator=(const RActionBase &) = delete;
   virtual ~RActionBase() { fLoopManager->Deregister(this); }

   virtual void Run(unsigned int slot, Long64_t entry) = 0;
   virtual void Initialize() = 0;
   virtual void InitSlot(TTreeReader *r, unsigned int slot) = 0;
   virtual void TriggerChildrenCount() = 0;
   virtual void ClearValueReaders(unsigned int slot) = 0;
   virtual void FinalizeSlot(unsigned int) = 0;
   virtual void Finalize() = 0;
   /// This method is invoked to update a partial result during the event loop, right before passing the result to a
   /// user-defined callback registered via RResultPtr::RegisterCallback
   virtual void *PartialUpdate(unsigned int slot) = 0;
   virtual bool HasRun() const = 0;

   virtual std::shared_ptr< ROOT::Internal::RDF::GraphDrawing::GraphNode> GetGraph() = 0;
};

class RJittedAction : public RActionBase {
private:
   std::unique_ptr<RActionBase> fConcreteAction;

public:
   RJittedAction(RLoopManager &lm) : RActionBase(&lm, lm.GetNSlots(), RDFInternal::RBookedCustomColumns()) {}

   void SetAction(std::unique_ptr<RActionBase> a) { fConcreteAction = std::move(a); }

   void Run(unsigned int slot, Long64_t entry) final;
   void Initialize() final;
   void InitSlot(TTreeReader *r, unsigned int slot) final;
   void TriggerChildrenCount() final;
   void FinalizeSlot(unsigned int) final;
   void Finalize() final;
   void *PartialUpdate(unsigned int slot) final;
   bool HasRun() const final;
   void ClearValueReaders(unsigned int slot) final;

   std::shared_ptr< ROOT::Internal::RDF::GraphDrawing::GraphNode> GetGraph();
};

template <typename Helper, typename PrevDataFrame, typename ColumnTypes_t = typename Helper::ColumnTypes_t>
class RAction final : public RActionBase {
   using TypeInd_t = std::make_index_sequence<ColumnTypes_t::list_size>;

   Helper fHelper;
   const ColumnNames_t fBranches;
   const std::shared_ptr<PrevDataFrame> fPrevDataPtr;
   PrevDataFrame &fPrevData;
   std::vector<RDFValueTuple_t<ColumnTypes_t>> fValues;
   bool fHasRun = false;

public:
   RAction(Helper &&h, const ColumnNames_t &bl, std::shared_ptr<PrevDataFrame> pd, const RBookedCustomColumns &customColumns)
      : RActionBase(pd->GetLoopManagerUnchecked(), pd->GetLoopManagerUnchecked()->GetNSlots(), customColumns), fHelper(std::move(h)),
        fBranches(bl), fPrevDataPtr(std::move(pd)), fPrevData(*fPrevDataPtr), fValues(fNSlots)
   {
   }

   RAction(const RAction &) = delete;
   RAction &operator=(const RAction &) = delete;

   void Initialize() final {
      fHelper.Initialize();
   }

   void InitSlot(TTreeReader *r, unsigned int slot) final
   {
      for (auto &bookedBranch : fCustomColumns.GetColumns())
         bookedBranch.second->InitSlot(r, slot);

      InitRDFValues(slot, fValues[slot], r, fBranches, fCustomColumns, TypeInd_t());
      fHelper.InitTask(r, slot);
   }

   void Run(unsigned int slot, Long64_t entry) final
   {
      // check if entry passes all filters
      if (fPrevData.CheckFilters(slot, entry))
         Exec(slot, entry, TypeInd_t());
   }

   template <std::size_t... S>
   void Exec(unsigned int slot, Long64_t entry, std::index_sequence<S...>)
   {
      (void)entry; // avoid bogus 'unused parameter' warning in gcc4.9
      fHelper.Exec(slot, std::get<S>(fValues[slot]).Get(entry)...);
   }

   void TriggerChildrenCount() final { fPrevData.IncrChildrenCount(); }

   void FinalizeSlot(unsigned int slot) final
   {
      ClearValueReaders(slot);
      for (auto &column : fCustomColumns.GetColumns()) {
         column.second->ClearValueReaders(slot);
      }
      fHelper.CallFinalizeTask(slot);
   }

   void ClearValueReaders(unsigned int slot) { ResetRDFValueTuple(fValues[slot], TypeInd_t()); }

   void Finalize() final
   {
      fHelper.Finalize();
      fHasRun = true;
   }

   std::shared_ptr<RDFGraphDrawing::GraphNode> GetGraph()
   {
      auto prevNode = fPrevData.GetGraph();
      auto prevColumns = prevNode->GetDefinedColumns();

      // Action nodes do not need to ask an helper to create the graph nodes. They are never common nodes between
      // multiple branches
      auto thisNode = std::make_shared< RDFGraphDrawing::GraphNode>(fHelper.GetActionName());
      auto evaluatedNode = thisNode;
      for (auto &column : fCustomColumns.GetColumns()) {
         /* Each column that this node has but the previous hadn't has been defined in between,
          * so it has to be built and appended. */
         if (RDFGraphDrawing::CheckIfDefaultOrDSColumn(column.first, column.second))
            continue;
         if (std::find(prevColumns.begin(), prevColumns.end(), column.first) == prevColumns.end()) {
            auto defineNode = RDFGraphDrawing::CreateDefineNode(column.first, column.second.get());
            evaluatedNode->SetPrevNode(defineNode);
            evaluatedNode = defineNode;
         }
      }

      thisNode->AddDefinedColumns(fCustomColumns.GetNames());
      thisNode->SetAction(HasRun());
      evaluatedNode->SetPrevNode(prevNode);
      return thisNode;
   }

   /// This method is invoked to update a partial result during the event loop, right before passing the result to a
   /// user-defined callback registered via RResultPtr::RegisterCallback
   void *PartialUpdate(unsigned int slot) final { return PartialUpdateImpl(slot); }

   bool HasRun() const { return fHasRun; }

private:
   // this overload is SFINAE'd out if Helper does not implement `PartialUpdate`
   // the template parameter is required to defer instantiation of the method to SFINAE time
   template <typename H = Helper>
   auto PartialUpdateImpl(unsigned int slot) -> decltype(std::declval<H>().PartialUpdate(slot), (void *)(nullptr))
   {
      return &fHelper.PartialUpdate(slot);
   }
   // this one is always available but has lower precedence thanks to `...`
   void *PartialUpdateImpl(...) { throw std::runtime_error("This action does not support callbacks yet!"); }
};

} // namespace RDF
} // namespace Internal

namespace Detail {
namespace RDF {
namespace RDFGraphDrawing = ROOT::Internal::RDF::GraphDrawing;

template <typename F, typename ExtraArgsTag = CustomColExtraArgs::None>
class RCustomColumn final : public RCustomColumnBase {
   // shortcuts
   using NoneTag = CustomColExtraArgs::None;
   using SlotTag = CustomColExtraArgs::Slot;
   using SlotAndEntryTag = CustomColExtraArgs::SlotAndEntry;
   // other types
   using FunParamTypes_t = typename CallableTraits<F>::arg_types;
   using ColumnTypesTmp_t =
      RDFInternal::RemoveFirstParameterIf_t<std::is_same<ExtraArgsTag, SlotTag>::value, FunParamTypes_t>;
   using ColumnTypes_t =
      RDFInternal::RemoveFirstTwoParametersIf_t<std::is_same<ExtraArgsTag, SlotAndEntryTag>::value, ColumnTypesTmp_t>;
   using TypeInd_t = std::make_index_sequence<ColumnTypes_t::list_size>;
   using ret_type = typename CallableTraits<F>::ret_type;
   // Avoid instantiating vector<bool> as `operator[]` returns temporaries in that case. Use std::deque instead.
   using ValuesPerSlot_t =
      typename std::conditional<std::is_same<ret_type, bool>::value, std::deque<ret_type>, std::vector<ret_type>>::type;

   F fExpression;
   const ColumnNames_t fBranches;
   ValuesPerSlot_t fLastResults;

   std::vector<RDFInternal::RDFValueTuple_t<ColumnTypes_t>> fValues;

public:
   RCustomColumn(RLoopManager *lm, std::string_view name, F &&expression, const ColumnNames_t &bl, unsigned int nSlots,
                 const RDFInternal::RBookedCustomColumns &customColumns, bool isDSColumn = false)
      : RCustomColumnBase(lm, name, nSlots, isDSColumn, customColumns), fExpression(std::move(expression)), fBranches(bl),
        fLastResults(fNSlots), fValues(fNSlots)
   {
   }

   RCustomColumn(const RCustomColumn &) = delete;
   RCustomColumn &operator=(const RCustomColumn &) = delete;


   void InitSlot(TTreeReader *r, unsigned int slot) final
   {
       //TODO: Each node calls this method for each column it uses. Multiple nodes may share the same columns, and this would lead to this method being called multiple times.
      RDFInternal::InitRDFValues(slot, fValues[slot], r, fBranches, fCustomColumns, TypeInd_t());
   }

   void *GetValuePtr(unsigned int slot) final { return static_cast<void *>(&fLastResults[slot]); }

   void Update(unsigned int slot, Long64_t entry) final
   {
      if (entry != fLastCheckedEntry[slot]) {
         // evaluate this filter, cache the result
         UpdateHelper(slot, entry, TypeInd_t(), ColumnTypes_t(), ExtraArgsTag{});
         fLastCheckedEntry[slot] = entry;
      }
   }

   const std::type_info &GetTypeId() const
   {
      return fIsDataSourceColumn ? typeid(typename std::remove_pointer<ret_type>::type) : typeid(ret_type);
   }

   template <std::size_t... S, typename... BranchTypes>
   void UpdateHelper(unsigned int slot, Long64_t entry, std::index_sequence<S...>, TypeList<BranchTypes...>, NoneTag)
   {
      fLastResults[slot] = fExpression(std::get<S>(fValues[slot]).Get(entry)...);
      // silence "unused parameter" warnings in gcc
      (void)slot;
      (void)entry;
   }

   template <std::size_t... S, typename... BranchTypes>
   void UpdateHelper(unsigned int slot, Long64_t entry, std::index_sequence<S...>, TypeList<BranchTypes...>, SlotTag)
   {
      fLastResults[slot] = fExpression(slot, std::get<S>(fValues[slot]).Get(entry)...);
      // silence "unused parameter" warnings in gcc
      (void)slot;
      (void)entry;
   }

   template <std::size_t... S, typename... BranchTypes>
   void
   UpdateHelper(unsigned int slot, Long64_t entry, std::index_sequence<S...>, TypeList<BranchTypes...>, SlotAndEntryTag)
   {
      fLastResults[slot] = fExpression(slot, entry, std::get<S>(fValues[slot]).Get(entry)...);
      // silence "unused parameter" warnings in gcc
      (void)slot;
      (void)entry;
   }

   void ClearValueReaders(unsigned int slot) final
   {
      //TODO: Each node calls this method for each column it uses. Multiple nodes may share the same columns, and this would lead to this method being called multiple times.
      RDFInternal::ResetRDFValueTuple(fValues[slot], TypeInd_t());
   }
};

class RFilterBase : public RNode {
protected:
   std::vector<Long64_t> fLastCheckedEntry;
   std::vector<int> fLastResult = {true}; // std::vector<bool> cannot be used in a MT context safely
   std::vector<ULong64_t> fAccepted = {0};
   std::vector<ULong64_t> fRejected = {0};
   const std::string fName;
   const unsigned int fNSlots;      ///< Number of thread slots used by this node, inherited from parent node.

   RDFInternal::RBookedCustomColumns fCustomColumns;

public:
   RFilterBase(RLoopManager *df, std::string_view name, const unsigned int nSlots,
               const RDFInternal::RBookedCustomColumns &customColumns);
   RFilterBase &operator=(const RFilterBase &) = delete;
   virtual ~RFilterBase() { fLoopManager->Deregister(this); }

   virtual void InitSlot(TTreeReader *r, unsigned int slot) = 0;
   bool HasName() const;
   std::string GetName() const;
   virtual void FillReport(ROOT::RDF::RCutFlowReport &) const;
   virtual void TriggerChildrenCount() = 0;
   virtual void ResetReportCount()
   {
      R__ASSERT(!fName.empty()); // this method is to only be called on named filters
      // fAccepted and fRejected could be different than 0 if this is not the first event-loop run using this filter
      std::fill(fAccepted.begin(), fAccepted.end(), 0);
      std::fill(fRejected.begin(), fRejected.end(), 0);
   }
   virtual void ClearValueReaders(unsigned int slot) = 0;
   virtual void ClearTask(unsigned int slot) = 0;
   virtual void InitNode();
   virtual void AddFilterName(std::vector<std::string> &filters) = 0;
   virtual std::shared_ptr<RDFGraphDrawing::GraphNode> GetGraph() = 0;
};

/// A wrapper around a concrete RFilter, which forwards all calls to it
/// RJittedFilter is the type of the node returned by jitted Filter calls: the concrete filter can be created and set
/// at a later time, from jitted code.
class RJittedFilter final : public RFilterBase {
   std::unique_ptr<RFilterBase> fConcreteFilter = nullptr;

public:
   RJittedFilter(RLoopManager *lm, std::string_view name)
      : RFilterBase(lm, name, lm->GetNSlots(), RDFInternal::RBookedCustomColumns())
   {
   }

   void SetFilter(std::unique_ptr<RFilterBase> f);

   void InitSlot(TTreeReader *r, unsigned int slot) final;
   bool CheckFilters(unsigned int slot, Long64_t entry) final;
   void Report(ROOT::RDF::RCutFlowReport &) const final;
   void PartialReport(ROOT::RDF::RCutFlowReport &) const final;
   void FillReport(ROOT::RDF::RCutFlowReport &) const final;
   void IncrChildrenCount() final;
   void StopProcessing() final;
   void ResetChildrenCount() final;
   void TriggerChildrenCount() final;
   void ResetReportCount() final;
   void ClearValueReaders(unsigned int slot) final;
   void InitNode() final;
   void AddFilterName(std::vector<std::string> &filters) final;
   void ClearTask(unsigned int slot) final;

   std::shared_ptr<RDFGraphDrawing::GraphNode> GetGraph(){
      if(fConcreteFilter != nullptr ){
         //Here the filter exists, so it can be served
         return fConcreteFilter->GetGraph();
      }
      throw std::runtime_error("The Jitting should have been invoked before this method.");
   }

};

template <typename FilterF, typename PrevDataFrame>
class RFilter final : public RFilterBase {
   using ColumnTypes_t = typename CallableTraits<FilterF>::arg_types;
   using TypeInd_t = std::make_index_sequence<ColumnTypes_t::list_size>;

   FilterF fFilter;
   const ColumnNames_t fBranches;
   const std::shared_ptr<PrevDataFrame> fPrevDataPtr;
   PrevDataFrame &fPrevData;
   std::vector<RDFInternal::RDFValueTuple_t<ColumnTypes_t>> fValues;

public:
   RFilter(FilterF &&f, const ColumnNames_t &bl, std::shared_ptr<PrevDataFrame> pd,
           const RDFInternal::RBookedCustomColumns &customColumns, std::string_view name = "")
      : RFilterBase(pd->GetLoopManagerUnchecked(), name, pd->GetLoopManagerUnchecked()->GetNSlots(), customColumns),
        fFilter(std::move(f)), fBranches(bl), fPrevDataPtr(std::move(pd)), fPrevData(*fPrevDataPtr), fValues(fNSlots)
   {
   }

   RFilter(const RFilter &) = delete;
   RFilter &operator=(const RFilter &) = delete;

   bool CheckFilters(unsigned int slot, Long64_t entry) final
   {
      if (entry != fLastCheckedEntry[slot]) {
         if (!fPrevData.CheckFilters(slot, entry)) {
            // a filter upstream returned false, cache the result
            fLastResult[slot] = false;
         } else {
            // evaluate this filter, cache the result
            auto passed = CheckFilterHelper(slot, entry, TypeInd_t());
            passed ? ++fAccepted[slot] : ++fRejected[slot];
            fLastResult[slot] = passed;
         }
         fLastCheckedEntry[slot] = entry;
      }
      return fLastResult[slot];
   }

   template <std::size_t... S>
   bool CheckFilterHelper(unsigned int slot, Long64_t entry, std::index_sequence<S...>)
   {
      return fFilter(std::get<S>(fValues[slot]).Get(entry)...);
      // silence "unused parameter" warnings in gcc
      (void)slot;
      (void)entry;
   }

   void InitSlot(TTreeReader *r, unsigned int slot) final
   {
      for (auto &bookedBranch : fCustomColumns.GetColumns())
         bookedBranch.second->InitSlot(r, slot);
      RDFInternal::InitRDFValues(slot, fValues[slot], r, fBranches, fCustomColumns, TypeInd_t());
   }

   // recursive chain of `Report`s
   void Report(ROOT::RDF::RCutFlowReport &rep) const final { PartialReport(rep); }

   void PartialReport(ROOT::RDF::RCutFlowReport &rep) const final
   {
      fPrevData.PartialReport(rep);
      FillReport(rep);
   }

   void StopProcessing() final
   {
      ++fNStopsReceived;
      if (fNStopsReceived == fNChildren)
         fPrevData.StopProcessing();
   }

   void IncrChildrenCount() final
   {
      ++fNChildren;
      // propagate "children activation" upstream. named filters do the propagation via `TriggerChildrenCount`.
      if (fNChildren == 1 && fName.empty())
         fPrevData.IncrChildrenCount();
   }

   void TriggerChildrenCount() final
   {
      R__ASSERT(!fName.empty()); // this method is to only be called on named filters
      fPrevData.IncrChildrenCount();
   }

   virtual void ClearValueReaders(unsigned int slot) final
   {
      RDFInternal::ResetRDFValueTuple(fValues[slot], TypeInd_t());
   }

   void AddFilterName(std::vector<std::string> &filters)
   {
      fPrevData.AddFilterName(filters);
      auto name = (HasName() ? fName : "Unnamed Filter");
      filters.push_back(name);
   }

   virtual void ClearTask(unsigned int slot) final
   {
      for (auto &column : fCustomColumns.GetColumns()) {
         column.second->ClearValueReaders(slot);
      }

      ClearValueReaders(slot);
   }

   std::shared_ptr<RDFGraphDrawing::GraphNode> GetGraph(){
      // Recursively call for the previous node.
      auto prevNode = fPrevData.GetGraph();
      auto prevColumns = prevNode->GetDefinedColumns();

      auto thisNode = RDFGraphDrawing::CreateFilterNode(this);

      /* If the returned node is not new, there is no need to perform any other operation.
       * This is a likely scenario when building the entire graph in which branches share
       * some nodes. */
      if(!thisNode->GetIsNew()){
         return thisNode;
      }

      auto evaluatedNode = thisNode;
      /* Each column that this node has but the previous hadn't has been defined in between,
       * so it has to be built and appended. */

      for (auto &column: fCustomColumns.GetColumns()){
         // Even if treated as custom columns by the Dataframe, datasource columns must not be in the graph.
         if(RDFGraphDrawing::CheckIfDefaultOrDSColumn(column.first, column.second))
            continue;
         if(std::find(prevColumns.begin(), prevColumns.end(), column.first) == prevColumns.end()){
            auto defineNode = RDFGraphDrawing::CreateDefineNode(column.first, column.second.get());
            evaluatedNode->SetPrevNode(defineNode);
            evaluatedNode = defineNode;
         }
      }

      // Keep track of the columns defined up to this point.
      thisNode->AddDefinedColumns(fCustomColumns.GetNames());

      evaluatedNode->SetPrevNode(prevNode);
      return thisNode;
   }
};

class RRangeBase : public RNode {
protected:
   unsigned int fStart;
   unsigned int fStop;
   unsigned int fStride;
   Long64_t fLastCheckedEntry{-1};
   bool fLastResult{true};
   ULong64_t fNProcessedEntries{0};
   bool fHasStopped{false};         ///< True if the end of the range has been reached
   const unsigned int fNSlots;      ///< Number of thread slots used by this node, inherited from parent node.

   void ResetCounters();

public:
   RRangeBase(RLoopManager *implPtr, unsigned int start, unsigned int stop, unsigned int stride,
              const unsigned int nSlots);

   RRangeBase &operator=(const RRangeBase &) = delete;
   virtual ~RRangeBase() { fLoopManager->Deregister(this); }

   void InitNode() { ResetCounters(); }
   virtual std::shared_ptr<RDFGraphDrawing::GraphNode> GetGraph() = 0;
};

template <typename PrevData>
class RRange final : public RRangeBase {
   const std::shared_ptr<PrevData> fPrevDataPtr;
   PrevData &fPrevData;

public:
   RRange(unsigned int start, unsigned int stop, unsigned int stride, std::shared_ptr<PrevData> pd)
      : RRangeBase(pd->GetLoopManagerUnchecked(), start, stop, stride, pd->GetLoopManagerUnchecked()->GetNSlots()),
        fPrevDataPtr(std::move(pd)), fPrevData(*fPrevDataPtr)
   {
   }

   RRange(const RRange &) = delete;
   RRange &operator=(const RRange &) = delete;

   /// Ranges act as filters when it comes to selecting entries that downstream nodes should process
   bool CheckFilters(unsigned int slot, Long64_t entry) final
   {
      if (entry != fLastCheckedEntry) {
         if (fHasStopped)
            return false;
         if (!fPrevData.CheckFilters(slot, entry)) {
            // a filter upstream returned false, cache the result
            fLastResult = false;
         } else {
            // apply range filter logic, cache the result
            ++fNProcessedEntries;
            if (fNProcessedEntries <= fStart || (fStop > 0 && fNProcessedEntries > fStop) ||
                (fStride != 1 && fNProcessedEntries % fStride != 0))
               fLastResult = false;
            else
               fLastResult = true;
            if (fNProcessedEntries == fStop) {
               fHasStopped = true;
               fPrevData.StopProcessing();
            }
         }
         fLastCheckedEntry = entry;
      }
      return fLastResult;
   }

   // recursive chain of `Report`s
   // RRange simply forwards these calls to the previous node
   void Report(ROOT::RDF::RCutFlowReport &rep) const final { fPrevData.PartialReport(rep); }

   void PartialReport(ROOT::RDF::RCutFlowReport &rep) const final { fPrevData.PartialReport(rep); }

   void StopProcessing() final
   {
      ++fNStopsReceived;
      if (fNStopsReceived == fNChildren && !fHasStopped)
         fPrevData.StopProcessing();
   }

   void IncrChildrenCount() final
   {
      ++fNChildren;
      // propagate "children activation" upstream
      if (fNChildren == 1)
         fPrevData.IncrChildrenCount();
   }

   /// This function must be defined by all nodes, but only the filters will add their name
   void AddFilterName(std::vector<std::string> &filters) { fPrevData.AddFilterName(filters); }
   std::shared_ptr<RDFGraphDrawing::GraphNode> GetGraph()
   {
      // TODO: Ranges node have no information about custom columns, hence it is not possible now
      // if defines have been used before.
      auto prevNode = fPrevData.GetGraph();
      auto prevColumns = prevNode->GetDefinedColumns();

      auto thisNode = RDFGraphDrawing::CreateRangeNode(this);

      /* If the returned node is not new, there is no need to perform any other operation.
       * This is a likely scenario when building the entire graph in which branches share
       * some nodes. */
      if (!thisNode->GetIsNew()) {
         return thisNode;
      }
      thisNode->SetPrevNode(prevNode);

      // If there have been some defines before it, this node won't detect them.
      thisNode->AddDefinedColumns(prevColumns);

      return thisNode;
   }
};

} // namespace RDF
} // namespace Detail

// method implementations
namespace Internal {
namespace RDF {

template <typename T>
void TColumnValue<T>::SetTmpColumn(unsigned int slot, ROOT::Detail::RDF::RCustomColumnBase *customColumn)
{
   fCustomColumns.emplace(customColumn);
   // Here we compare names and not typeinfos since they may come from two different contexts: a compiled
   // and a jitted one.
   if (0 != strcmp(customColumn->GetTypeId().name(), typeid(T).name()))
      throw std::runtime_error(
         std::string("TColumnValue: type specified for column \"" + customColumn->GetName() + "\" is ") +
         TypeID2TypeName(typeid(T)) + " but temporary column has type " + TypeID2TypeName(customColumn->GetTypeId()));

   if (customColumn->IsDataSourceColumn()) {
      fColumnKind = EColumnKind::kDataSource;
      fDSValuePtrs.emplace(static_cast<T **>(customColumn->GetValuePtr(slot)));
   } else {
      fColumnKind = EColumnKind::kCustomColumn;
      fCustomValuePtrs.emplace(static_cast<T *>(customColumn->GetValuePtr(slot)));
   }
   fSlot = slot;
}

// This method is executed inside the event-loop, many times per entry
// If need be, the if statement can be avoided using thunks
// (have both branches inside functions and have a pointer to the branch to be executed)
template <typename T>
template <typename U, typename std::enable_if<!TColumnValue<U>::MustUseRVec_t::value, int>::type>
T &TColumnValue<T>::Get(Long64_t entry)
{
   if (fColumnKind == EColumnKind::kTree) {
      return *(fTreeReaders.top()->Get());
   } else {
      fCustomColumns.top()->Update(fSlot, entry);
      return fColumnKind == EColumnKind::kCustomColumn ? *fCustomValuePtrs.top() : **fDSValuePtrs.top();
   }
}

/// This overload is used to return arrays (i.e. types that are read into a RVec)
template <typename T>
template <typename U, typename std::enable_if<TColumnValue<U>::MustUseRVec_t::value, int>::type>
T &TColumnValue<T>::Get(Long64_t entry)
{
   if (fColumnKind == EColumnKind::kTree) {
      auto &readerArray = *fTreeReaders.top();
      // We only use TTreeReaderArrays to read columns that users flagged as type `RVec`, so we need to check
      // that the branch stores the array as contiguous memory that we can actually wrap in an `RVec`.
      // Currently we need the first entry to have been loaded to perform the check
      // TODO Move check to `MakeProxy` once Axel implements this kind of check in TTreeReaderArray using
      // TBranchProxy

      if (EStorageType::kUnknown == fStorageType && readerArray.GetSize() > 1) {
         // We can decide since the array is long enough
         fStorageType = (1 == (&readerArray[1] - &readerArray[0])) ? EStorageType::kContiguous : EStorageType::kSparse;
      }

      const auto readerArraySize = readerArray.GetSize();
      if (EStorageType::kContiguous == fStorageType ||
          (EStorageType::kUnknown == fStorageType && readerArray.GetSize() < 2)) {
         if (readerArraySize > 0) {
            // trigger loading of the contens of the TTreeReaderArray
            // the address of the first element in the reader array is not necessarily equal to
            // the address returned by the GetAddress method
            auto readerArrayAddr = &readerArray.At(0);
            T tvec(readerArrayAddr, readerArraySize);
            swap(fRVec, tvec);
         } else {
            T emptyVec{};
            swap(fRVec, emptyVec);
         }
      } else {
// The storage is not contiguous or we don't know yet: we cannot but copy into the tvec
#ifndef NDEBUG
         if (!fCopyWarningPrinted) {
            Warning("TColumnValue::Get", "Branch %s hangs from a non-split branch. A copy is being performed in order "
                                         "to properly read the content.",
                    readerArray.GetBranchName());
            fCopyWarningPrinted = true;
         }
#else
         (void)fCopyWarningPrinted;
#endif
         if (readerArraySize > 0) {
            (void)readerArray.At(0); // trigger deserialisation
            T tvec(readerArray.begin(), readerArray.end());
            swap(fRVec, tvec);
         } else {
            T emptyVec{};
            swap(fRVec, emptyVec);
         }
      }
      return fRVec;

   } else {
      fCustomColumns.top()->Update(fSlot, entry);
      return fColumnKind == EColumnKind::kCustomColumn ? *fCustomValuePtrs.top() : **fDSValuePtrs.top();
   }
}

} // namespace RDF
} // namespace Internal
} // namespace ROOT
#endif // ROOT_RDFNODES
