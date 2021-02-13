// Author: Stefan Wunsch CERN  04/2019

/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

// Include Python.h first before any standard header
#include "Python.h"

#include "ROOT/RIntegerSequence.hxx"
#include "ROOT/RMakeUnique.hxx"
#include "ROOT/RDataSource.hxx"
#include "ROOT/TSeq.hxx"
#include "ROOT/RVec.hxx"

#include <algorithm>
#include <map>
#include <tuple>
#include <string>
#include <typeinfo>
#include <vector>

#ifndef ROOT_RNUMPYDS
#define ROOT_RNUMPYDS

namespace ROOT {

namespace Internal {

namespace RDF {

////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief A RDataSource implementation which takes a collection of RVecs, which
/// are able to adopt data from Numpy arrays
///
/// This component allows to create a data source on a set of columns with data
/// coming from RVecs. The adoption of externally provided data, e.g., via Numpy
/// arrays, with RVecs allows to read arbitrary data from memory.
/// In addition, the data source has to keep a reference on the Python owned data
/// so that the lifetime of the data is tied to the datasource.
template <typename... ColumnTypes>
class RNumpyDS final : public ROOT::RDF::RDataSource {
   using PointerHolderPtrs_t = std::vector<ROOT::Internal::TDS::TPointerHolder *>;

   std::tuple<ROOT::RVec<ColumnTypes>*...> fColumns;
   const std::vector<std::string> fColNames;
   const std::map<std::string, std::string> fColTypesMap;
   // The role of the fPointerHoldersModels is to be initialised with the pack
   // of arguments in the constrcutor signature at construction time
   // Once the number of slots is known, the fPointerHolders are initialised
   // according to the models.
   const PointerHolderPtrs_t fPointerHoldersModels;
   std::vector<PointerHolderPtrs_t> fPointerHolders;
   std::vector<std::pair<ULong64_t, ULong64_t>> fEntryRanges{};
   unsigned int fNSlots{0};
   // Pointer to PyObject holding RVecs
   // The RVecs itself hold a reference to the associated Numpy arrays so that
   // the data cannot go out of scope as long as the datasource survives.
   PyObject* fPyRVecs;

   Record_t GetColumnReadersImpl(std::string_view colName, const std::type_info &id)
   {
      auto colNameStr = std::string(colName);
      // This could be optimised and done statically
      const auto idName = ROOT::Internal::RDF::TypeID2TypeName(id);
      auto it = fColTypesMap.find(colNameStr);
      if (fColTypesMap.end() == it) {
         std::string err = "The specified column name, \"" + colNameStr + "\" is not known to the data source.";
         throw std::runtime_error(err);
      }

      const auto colIdName = it->second;
      if (colIdName != idName) {
         std::string err = "Column " + colNameStr + " has type " + colIdName +
                           " while the id specified is associated to type " + idName;
         throw std::runtime_error(err);
      }

      const auto colBegin = fColNames.begin();
      const auto colEnd = fColNames.end();
      const auto namesIt = std::find(colBegin, colEnd, colName);
      const auto index = std::distance(colBegin, namesIt);

      Record_t ret(fNSlots);
      for (auto slot : ROOT::TSeqU(fNSlots)) {
         ret[slot] = fPointerHolders[index][slot]->GetPointerAddr();
      }
      return ret;
   }

   size_t GetEntriesNumber() { return std::get<0>(fColumns)->size(); }
   template <std::size_t... S>
   void SetEntryHelper(unsigned int slot, ULong64_t entry, std::index_sequence<S...>)
   {
      std::initializer_list<int> expander{
         (*static_cast<ColumnTypes *>(fPointerHolders[S][slot]->GetPointer()) = (*std::get<S>(fColumns))[entry], 0)...};
      (void)expander; // avoid unused variable warnings
   }

   template <std::size_t... S>
   void ColLenghtChecker(std::index_sequence<S...>)
   {
      if (sizeof...(S) < 2)
         return;

      const std::vector<size_t> colLengths{std::get<S>(fColumns)->size()...};
      const auto expectedLen = colLengths[0];
      std::string err;
      for (auto i : TSeqI(1, colLengths.size())) {
         if (expectedLen != colLengths[i]) {
            err += "Column \"" + fColNames[i] + "\" and column \"" + fColNames[0] +
                   "\" have different lengths: " + std::to_string(expectedLen) + " and " +
                   std::to_string(colLengths[i]);
         }
      }
      if (!err.empty()) {
         throw std::runtime_error(err);
      }
   }

protected:
   std::string AsString() { return "Numpy data source"; };

public:
   RNumpyDS(PyObject* pyRVecs,
                  std::pair<std::string, ROOT::RVec<ColumnTypes>*>... colsNameVals)
      : fColumns(std::tuple<ROOT::RVec<ColumnTypes>*...>(colsNameVals.second...)),
        fColNames({colsNameVals.first...}),
        fColTypesMap({{colsNameVals.first, ROOT::Internal::RDF::TypeID2TypeName(typeid(ColumnTypes))}...}),
        fPointerHoldersModels({new ROOT::Internal::TDS::TTypedPointerHolder<ColumnTypes>(new ColumnTypes())...}),
        fPyRVecs(pyRVecs)
   {
      // Take a reference to the data associated with this data source
      Py_INCREF(fPyRVecs);
   }

   ~RNumpyDS()
   {
      for (auto &&ptrHolderv : fPointerHolders) {
         for (auto &&ptrHolder : ptrHolderv) {
            delete ptrHolder;
         }
      }
      // Release the data associated to this data source
      Py_DECREF(fPyRVecs);
   }

   const std::vector<std::string> &GetColumnNames() const { return fColNames; }

   std::vector<std::pair<ULong64_t, ULong64_t>> GetEntryRanges()
   {
      auto entryRanges(std::move(fEntryRanges)); // empty fEntryRanges
      return entryRanges;
   }

   std::string GetTypeName(std::string_view colName) const
   {
      const auto key = std::string(colName);
      return fColTypesMap.at(key);
   }

   bool HasColumn(std::string_view colName) const
   {
      const auto key = std::string(colName);
      const auto endIt = fColTypesMap.end();
      return endIt != fColTypesMap.find(key);
   }

   bool SetEntry(unsigned int slot, ULong64_t entry)
   {
      SetEntryHelper(slot, entry, std::index_sequence_for<ColumnTypes...>());
      return true;
   }

   void SetNSlots(unsigned int nSlots)
   {
      fNSlots = nSlots;
      const auto nCols = fColNames.size();
      fPointerHolders.resize(nCols); // now we need to fill it with the slots, all of the same type
      auto colIndex = 0U;
      for (auto &&ptrHolderv : fPointerHolders) {
         for (auto slot : ROOT::TSeqI(fNSlots)) {
            auto ptrHolder = fPointerHoldersModels[colIndex]->GetDeepCopy();
            ptrHolderv.emplace_back(ptrHolder);
            (void)slot;
         }
         colIndex++;
      }
      for (auto &&ptrHolder : fPointerHoldersModels)
         delete ptrHolder;
   }

   void Initialise()
   {
      ColLenghtChecker(std::index_sequence_for<ColumnTypes...>());
      const auto nEntries = GetEntriesNumber();
      const auto nEntriesInRange = nEntries / fNSlots; // between integers. Should make smaller?
      auto reminder = 1U == fNSlots ? 0 : nEntries % fNSlots;
      fEntryRanges.resize(fNSlots);
      auto init = 0ULL;
      auto end = 0ULL;
      for (auto &&range : fEntryRanges) {
         end = init + nEntriesInRange;
         if (0 != reminder) { // Distribute the reminder among the first chunks
            reminder--;
            end += 1;
         }
         range.first = init;
         range.second = end;
         init = end;
      }
   }

   std::string GetLabel() { return "RNumpyDS"; }
};

// Factory to create datasource able to read Numpy arrays through RVecs
// Note that we have to return the object on the heap so that the interpreter
// does not clean it up during shutdown and causes a double delete.
template <typename... ColumnTypes>
RDataFrame* MakeNumpyDataFrame(PyObject* pyRVecs,
                              std::pair<std::string, ROOT::RVec<ColumnTypes>*> &&... colNameProxyPairs)
{
   return new RDataFrame(std::make_unique<RNumpyDS<ColumnTypes...>>(
      std::forward<PyObject*>(pyRVecs),
      std::forward<std::pair<std::string, ROOT::RVec<ColumnTypes>*>>(colNameProxyPairs)...));
}

} // namespace RDF
} // namespace Internal
} // namespace ROOT

#endif // ROOT_RNUMPYDS
