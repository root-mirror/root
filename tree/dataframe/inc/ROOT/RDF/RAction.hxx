// Author: Enrico Guiraud, Danilo Piparo CERN  09/2018

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_RACTION
#define ROOT_RACTION

#include "ROOT/RDF/ColumnReaderUtils.hxx"
#include "ROOT/RDF/GraphNode.hxx"
#include "ROOT/RDF/RActionBase.hxx"
#include "ROOT/RDF/RColumnReaderBase.hxx"
#include "ROOT/RDF/Utils.hxx" // ColumnNames_t, IsInternalColumn
#include "ROOT/RDF/RLoopManager.hxx"

#include <array>
#include <cstddef> // std::size_t
#include <memory>
#include <string>
#include <vector>

namespace ROOT {
namespace Internal {
namespace RDF {

namespace RDFDetail = ROOT::Detail::RDF;
namespace RDFGraphDrawing = ROOT::Internal::RDF::GraphDrawing;

namespace GraphDrawing {
std::shared_ptr<GraphNode> AddDefinesToGraph(std::shared_ptr<GraphNode> node,
                                             const RDFInternal::RBookedDefines &defines,
                                             const std::vector<std::string> &prevNodeDefines);
} // namespace GraphDrawing

// clang-format off
/**
 * \class ROOT::Internal::RDF::RAction
 * \ingroup dataframe
 * \brief A RDataFrame node that produces a result
 * \tparam Helper The action helper type, which implements the concrete action logic (e.g. FillHelper, SnapshotHelper)
 * \tparam PrevDataFrame The type of the parent node in the computation graph
 * \tparam ColumnTypes_t A TypeList with the types of the input columns
 *
 */
// clang-format on
template <typename Helper, typename PrevDataFrame, typename ColumnTypes_t = typename Helper::ColumnTypes_t>
class R__CLING_PTRCHECK(off) RAction : public RActionBase {
   using TypeInd_t = std::make_index_sequence<ColumnTypes_t::list_size>;

   Helper fHelper;
   const std::shared_ptr<PrevDataFrame> fPrevDataPtr;
   PrevDataFrame &fPrevData;
   /// Column readers per slot and per input column
   std::vector<std::array<std::unique_ptr<RColumnReaderBase>, ColumnTypes_t::list_size>> fValues;

   /// The nth flag signals whether the nth input column is a custom column or not.
   std::array<bool, ColumnTypes_t::list_size> fIsDefine;

public:
   RAction(Helper &&h, const ColumnNames_t &columns, std::shared_ptr<PrevDataFrame> pd, const RBookedDefines &defines)
      : RActionBase(pd->GetLoopManagerUnchecked(), columns, defines), fHelper(std::forward<Helper>(h)),
        fPrevDataPtr(std::move(pd)), fPrevData(*fPrevDataPtr), fValues(GetNSlots()), fIsDefine()
   {
      const auto nColumns = columns.size();
      const auto &customCols = GetDefines();
      for (auto i = 0u; i < nColumns; ++i)
         fIsDefine[i] = customCols.HasName(columns[i]);
   }

   RAction(const RAction &) = delete;
   RAction &operator=(const RAction &) = delete;
   // must call Deregister here, before fPrevDataFrame is destroyed,
   // otherwise if fPrevDataFrame is fLoopManager we get a use after delete
   ~RAction() { fLoopManager->Deregister(this); }

   /**
      Retrieve a wrapper to the result of the action that knows how to merge
      with others of the same type.
   */
   std::unique_ptr<RDFDetail::RMergeableValueBase> GetMergeableValue() const final
   {
      return fHelper.GetMergeableValue();
   }

   void Initialize() final { fHelper.Initialize(); }

   void InitSlot(TTreeReader *r, unsigned int slot) final
   {
      for (auto &bookedBranch : GetDefines().GetColumns())
         bookedBranch.second->InitSlot(r, slot);
      RDFInternal::RColumnReadersInfo info{RActionBase::GetColumnNames(), RActionBase::GetDefines(), fIsDefine.data(),
                                           fLoopManager->GetDSValuePtrs(), fLoopManager->GetDataSource()};
      fValues[slot] = RDFInternal::MakeColumnReaders(slot, r, ColumnTypes_t{}, info);
      fHelper.InitTask(r, slot);
   }

   template <typename... ColTypes, std::size_t... S>
   void CallExec(unsigned int slot, Long64_t entry, TypeList<ColTypes...>, std::index_sequence<S...>)
   {
      fHelper.Exec(slot, fValues[slot][S]->template Get<ColTypes>(entry)...);
      (void)entry; // avoid "unused parameter" warnings
   }

   void Run(unsigned int slot, Long64_t entry) final
   {
      // check if entry passes all filters
      if (fPrevData.CheckFilters(slot, entry))
         CallExec(slot, entry, ColumnTypes_t{}, TypeInd_t{});
   }

   void TriggerChildrenCount() final { fPrevData.IncrChildrenCount(); }

   /// Clean-up operations to be performed at the end of a task.
   void FinalizeSlot(unsigned int slot) final
   {
      for (auto &column : GetDefines().GetColumns())
         column.second->FinaliseSlot(slot);
      for (auto &v : fValues[slot])
         v.reset();
      fHelper.CallFinalizeTask(slot);
   }

   /// Clean-up and finalize the action result (e.g. merging slot-local results).
   /// It invokes the helper's Finalize method.
   void Finalize() final
   {
      fHelper.Finalize();
      SetHasRun();
   }

   std::shared_ptr<RDFGraphDrawing::GraphNode> GetGraph()
   {
      auto prevNode = fPrevData.GetGraph();
      auto prevColumns = prevNode->GetDefinedColumns();

      // Action nodes do not need to go through CreateFilterNode: they are never common nodes between multiple branches
      auto thisNode = std::make_shared<RDFGraphDrawing::GraphNode>(fHelper.GetActionName());

      auto upmostNode = AddDefinesToGraph(thisNode, GetDefines(), prevColumns);

      thisNode->AddDefinedColumns(GetDefines().GetNames());
      thisNode->SetAction(HasRun());
      upmostNode->SetPrevNode(prevNode);
      return thisNode;
   }

   /// This method is invoked to update a partial result during the event loop, right before passing the result to a
   /// user-defined callback registered via RResultPtr::RegisterCallback
   void *PartialUpdate(unsigned int slot) final { return PartialUpdateImpl(slot); }

private:
   // this overload is SFINAE'd out if Helper does not implement `PartialUpdate`
   // the template parameter is required to defer instantiation of the method to SFINAE time
   template <typename H = Helper>
   auto PartialUpdateImpl(unsigned int slot) -> decltype(std::declval<H>().PartialUpdate(slot), (void *)(nullptr))
   {
      return &fHelper.PartialUpdate(slot);
   }

   // this one is always available but has lower precedence thanks to `...`
   void *PartialUpdateImpl(...) { throw std::runtime_error("This action does not support callbacks!"); }
};

} // namespace RDF
} // namespace Internal
} // namespace ROOT

#endif // ROOT_RACTION
