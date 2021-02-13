// Author: Enrico Guiraud CERN 09/2020

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_RDF_COLUMNREADERUTILS
#define ROOT_RDF_COLUMNREADERUTILS

#include "RColumnReaderBase.hxx"
#include "RDefineBase.hxx"
#include "RDefineReader.hxx"
#include "RDSColumnReader.hxx"
#include "RTreeColumnReader.hxx"

#include <ROOT/RDataSource.hxx>
#include <ROOT/TypeTraits.hxx>
#include <TError.h> // R__ASSERT
#include <TTreeReader.h>

#include <array>
#include <map>
#include <memory>
#include <string>
#include <typeinfo> // for typeid
#include <vector>

namespace ROOT {
namespace Internal {
namespace RDF {

using namespace ROOT::TypeTraits;
namespace RDFDetail = ROOT::Detail::RDF;

template <typename T>
std::unique_ptr<RDFDetail::RColumnReaderBase>
MakeColumnReader(unsigned int slot, RDFDetail::RDefineBase *define, TTreeReader *r, ROOT::RDF::RDataSource *ds,
                 const std::vector<void *> *DSValuePtrsPtr, const std::string &colName)
{
   using Ret_t = std::unique_ptr<RDFDetail::RColumnReaderBase>;

   if (define != nullptr) // it's a Define'd column
      return Ret_t(new RDefineReader(slot, *define, typeid(T)));

   if (DSValuePtrsPtr != nullptr) {
      // reading from a RDataSource with the old column reader interface
      auto &DSValuePtrs = *DSValuePtrsPtr;
      return Ret_t(new RDSColumnReader<T>(DSValuePtrs[slot]));
   }

   if (ds != nullptr) {
      // reading from a RDataSource with the new column reader interface
      return ds->GetColumnReaders(slot, colName, typeid(T));
   }

   // reading from a TTree
   return Ret_t(new RTreeColumnReader<T>(*r, colName));
}

template <typename T>
std::unique_ptr<RDFDetail::RColumnReaderBase>
MakeColumnReadersHelper(unsigned int slot, RDFDetail::RDefineBase *define,
                        const std::map<std::string, std::vector<void *>> &DSValuePtrsMap, TTreeReader *r,
                        ROOT::RDF::RDataSource *ds, const std::string &colName)
{
   const auto DSValuePtrsIt = DSValuePtrsMap.find(colName);
   const std::vector<void *> *DSValuePtrsPtr = DSValuePtrsIt != DSValuePtrsMap.end() ? &DSValuePtrsIt->second : nullptr;
   R__ASSERT(define != nullptr || r != nullptr || DSValuePtrsPtr != nullptr || ds != nullptr);
   return MakeColumnReader<T>(slot, define, r, ds, DSValuePtrsPtr, colName);
}

/// This type aggregates some of the arguments passed to InitColumnReaders.
/// We need to pass a single RColumnReadersInfo object rather than each argument separately because with too many
/// arguments passed, gcc 7.5.0 and cling disagree on the ABI, which leads to the last function argument being read
/// incorrectly from a compiled InitColumnReaders symbols when invoked from a jitted symbol.
struct RColumnReadersInfo {
   const std::vector<std::string> &fColNames;
   const RBookedDefines &fCustomCols;
   const bool *fIsDefine;
   const std::map<std::string, std::vector<void *>> &fDSValuePtrsMap;
   ROOT::RDF::RDataSource *fDataSource;
};

/// Create a group of column readers, one per type in the parameter pack.
/// colInfo.fColNames and colInfo.fIsDefine are expected to have size equal to the parameter pack, and elements ordered
/// accordingly, i.e. fIsDefine[0] refers to fColNames[0] which is of type "ColTypes[0]".
template <typename... ColTypes>
std::array<std::unique_ptr<RDFDetail::RColumnReaderBase>, sizeof...(ColTypes)>
MakeColumnReaders(unsigned int slot, TTreeReader *r, TypeList<ColTypes...>, const RColumnReadersInfo &colInfo)
{
   // see RColumnReadersInfo for why we pass these arguments like this rather than directly as function arguments
   const auto &colNames = colInfo.fColNames;
   const auto &customCols = colInfo.fCustomCols;
   const bool *isDefine = colInfo.fIsDefine;
   const auto &DSValuePtrsMap = colInfo.fDSValuePtrsMap;
   auto *ds = colInfo.fDataSource;

   const auto &customColMap = customCols.GetColumns();

   int i = -1;
   std::array<std::unique_ptr<RDFDetail::RColumnReaderBase>, sizeof...(ColTypes)> ret{
      {{(++i, MakeColumnReadersHelper<ColTypes>(slot, isDefine[i] ? customColMap.at(colNames[i]).get() : nullptr,
                                                DSValuePtrsMap, r, ds, colNames[i]))}...}};
   return ret;

   // avoid bogus "unused variable" warnings
   (void)ds;
   (void)slot;
   (void)r;
}

} // namespace RDF
} // namespace Internal
} // namespace ROOT

#endif // ROOT_RDF_COLUMNREADERS
