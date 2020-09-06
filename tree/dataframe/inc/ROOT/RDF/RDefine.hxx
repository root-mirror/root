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

#include "ROOT/RDF/ColumnReaders.hxx"
#include "ROOT/RDF/RDefineBase.hxx"
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
class RDefine final : public RDefineBase {
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
   const ColumnNames_t fColumnNames;
   ValuesPerSlot_t fLastResults;

   std::vector<RDFInternal::RDFValueTuple_t> fValues;

   /// The nth flag signals whether the nth input column is a custom column or not.
   std::array<bool, ColumnTypes_t::list_size> fIsDefine;

   template <typename... ColTypes, std::size_t... S>
   void UpdateHelper(unsigned int slot, Long64_t entry, TypeList<ColTypes...>, std::index_sequence<S...>, NoneTag)
   {
      fLastResults[slot] = fExpression(fValues[slot][S]->template Get<ColTypes>(entry)...);
      // silence "unused parameter" warnings in gcc
      (void)slot;
      (void)entry;
   }

   template <typename... ColTypes, std::size_t... S>
   void UpdateHelper(unsigned int slot, Long64_t entry, TypeList<ColTypes...>, std::index_sequence<S...>, SlotTag)
   {
      fLastResults[slot] = fExpression(slot, fValues[slot][S]->template Get<ColTypes>(entry)...);
      // silence "unused parameter" warnings in gcc
      (void)slot;
      (void)entry;
   }

   template <typename... ColTypes, std::size_t... S>
   void
   UpdateHelper(unsigned int slot, Long64_t entry, TypeList<ColTypes...>, std::index_sequence<S...>, SlotAndEntryTag)
   {
      fLastResults[slot] = fExpression(slot, entry, fValues[slot][S]->template Get<ColTypes>(entry)...);
      // silence "unused parameter" warnings in gcc
      (void)slot;
      (void)entry;
   }

public:
   RDefine(std::string_view name, std::string_view type, F expression, const ColumnNames_t &columns,
                 unsigned int nSlots, const RDFInternal::RBookedDefines &defines,
                 const std::map<std::string, std::vector<void *>> &DSValuePtrs)
      : RDefineBase(name, type, nSlots, defines, DSValuePtrs), fExpression(std::move(expression)),
        fColumnNames(columns), fLastResults(fNSlots), fValues(fNSlots), fIsDefine()
   {
      const auto nColumns = fColumnNames.size();
      for (auto i = 0u; i < nColumns; ++i)
         fIsDefine[i] = fDefines.HasName(fColumnNames[i]);
   }

   RDefine(const RDefine &) = delete;
   RDefine &operator=(const RDefine &) = delete;

   void InitSlot(TTreeReader *r, unsigned int slot) final
   {
      if (!fIsInitialized[slot]) {
         fIsInitialized[slot] = true;
         RDFInternal::RColumnReadersInfo info{fColumnNames, fDefines, fIsDefine.data(), fDSValuePtrs};
         RDFInternal::InitColumnReaders(slot, fValues[slot], r, ColumnTypes_t{}, info);
         fLastCheckedEntry[slot] = -1;
      }
   }

   void *GetValuePtr(unsigned int slot) final { return static_cast<void *>(&fLastResults[slot]); }

   void Update(unsigned int slot, Long64_t entry) final
   {
      if (entry != fLastCheckedEntry[slot]) {
         // evaluate this filter, cache the result
         UpdateHelper(slot, entry, ColumnTypes_t{}, TypeInd_t{}, ExtraArgsTag{});
         fLastCheckedEntry[slot] = entry;
      }
   }

   const std::type_info &GetTypeId() const { return typeid(ret_type); }

   void ClearValueReaders(unsigned int slot) final
   {
      if (fIsInitialized[slot]) {
         for (auto &v : fValues[slot])
            v->Reset();
         fIsInitialized[slot] = false;
         fValues[slot].clear();
      }
   }
};

} // ns RDF
} // ns Detail
} // ns ROOT

#endif // ROOT_RCUSTOMCOLUMN
