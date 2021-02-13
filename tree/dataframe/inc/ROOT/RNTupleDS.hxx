/// \file RNTupleDS.hxx
/// \ingroup NTuple ROOT7
/// \author Jakob Blomer <jblomer@cern.ch>
/// \author Enrico Guiraud <enrico.guiraud@cern.ch>
/// \date 2018-10-04
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_RNTupleDS
#define ROOT_RNTupleDS

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDataSource.hxx>
#include <ROOT/RNTupleUtil.hxx>
#include <ROOT/RStringView.hxx>

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace ROOT {
namespace Experimental {

class RNTupleDescriptor;

namespace Detail {
class RFieldBase;
class RFieldValue;
class RPageSource;

} // namespace Detail

class RNTupleDS final : public ROOT::RDF::RDataSource {
   /// Clones of the first source, one for each slot
   std::vector<std::unique_ptr<ROOT::Experimental::Detail::RPageSource>> fSources;

   std::vector<std::string> fColumnNames;
   std::vector<std::string> fColumnTypes;
   std::vector<size_t> fActiveColumns;

   unsigned fNSlots = 0;
   bool fHasSeenAllRanges = false;

   void AddFields(const RNTupleDescriptor &desc, DescriptorId_t parentId);

public:
   explicit RNTupleDS(std::unique_ptr<ROOT::Experimental::Detail::RPageSource> pageSource);
   ~RNTupleDS() = default;
   void SetNSlots(unsigned int nSlots) final;
   const std::vector<std::string> &GetColumnNames() const final { return fColumnNames; }
   bool HasColumn(std::string_view colName) const final;
   std::string GetTypeName(std::string_view colName) const final;
   std::vector<std::pair<ULong64_t, ULong64_t>> GetEntryRanges() final;

   bool SetEntry(unsigned int slot, ULong64_t entry) final;

   void Initialise() final;
   void Finalise() final;

   std::unique_ptr<ROOT::Detail::RDF::RColumnReaderBase>
   GetColumnReaders(unsigned int /*slot*/, std::string_view /*name*/, const std::type_info &) final;

protected:
   Record_t GetColumnReadersImpl(std::string_view name, const std::type_info &) final;
};

RDataFrame MakeNTupleDataFrame(std::string_view ntupleName, std::string_view fileName);

} // ns Experimental
} // ns ROOT

#endif
