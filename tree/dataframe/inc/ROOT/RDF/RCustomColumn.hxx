// Author: Enrico Guiraud, Danilo Piparo CERN  09/2018

/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_RCUSTOMCOLUMN
#define ROOT_RCUSTOMCOLUMN

#include "ROOT/RDF/NodesUtils.hxx"
#include "ROOT/RDF/RColumnValue.hxx"
#include "ROOT/RDF/RCustomColumnBase.hxx"
#include "ROOT/RDF/Utils.hxx"
#include "ROOT/RIntegerSequence.hxx"
#include "ROOT/RStringView.hxx"
#include "ROOT/TypeTraits.hxx"
#include "RtypesCore.h"

#include <deque>
#include <type_traits>
#include <vector>

class TTreeReader;

namespace ROOT {
namespace Detail {
namespace RDF {

using namespace ROOT::TypeTraits;

// clang-format off
namespace CustomColExtraArgs {
struct None{};
struct Slot{};
struct SlotAndEntry{};
}
// clang-format on

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

public:
   RCustomColumn(RLoopManager *lm, std::string_view name, F &&expression, const ColumnNames_t &bl, unsigned int nSlots,
                 const RDFInternal::RBookedCustomColumns &customColumns, bool isDSColumn = false)
      : RCustomColumnBase(lm, name, nSlots, isDSColumn, customColumns), fExpression(std::forward<F>(expression)),
        fBranches(bl), fLastResults(fNSlots), fValues(fNSlots)
   {
   }

   RCustomColumn(const RCustomColumn &) = delete;
   RCustomColumn &operator=(const RCustomColumn &) = delete;

   void InitSlot(TTreeReader *r, unsigned int slot) final
   {
      // TODO: Each node calls this method for each column it uses. Multiple nodes may share the same columns, and this
      // would lead to this method being called multiple times.
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
};

} // ns RDF
} // ns Detail
} // ns ROOT

#endif // ROOT_RCUSTOMCOLUMN
