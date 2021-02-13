/// \file RNTupleDS.cxx
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

#include <ROOT/RDF/RColumnReaderBase.hxx>
#include <ROOT/RField.hxx>
#include <ROOT/RFieldValue.hxx>
#include <ROOT/RNTupleDescriptor.hxx>
#include <ROOT/RNTupleDS.hxx>
#include <ROOT/RNTupleUtil.hxx>
#include <ROOT/RPageStorage.hxx>
#include <ROOT/RStringView.hxx>

#include <TError.h>

#include <string>
#include <vector>
#include <typeinfo>
#include <utility>

namespace ROOT {
namespace Experimental {
namespace Detail {
class RNTupleColumnReader final : public ROOT::Detail::RDF::RColumnReaderBase {
   using RFieldBase = ROOT::Experimental::Detail::RFieldBase;
   using RFieldValue = ROOT::Experimental::Detail::RFieldValue;
   using RPageSource = ROOT::Experimental::Detail::RPageSource;

   std::unique_ptr<RFieldBase> fField;
   RFieldValue fValue;
   Long64_t fLastEntry; ///< Last entry number that was read

   std::unique_ptr<RFieldBase> MakeField(const std::string &colName, RPageSource &source)
   {
      const auto &descriptor = source.GetDescriptor();
      const auto fieldId = descriptor.FindFieldId(colName);
      const auto &fieldDescriptor = descriptor.GetFieldDescriptor(fieldId);
      const auto typeName = fieldDescriptor.GetTypeName();
      auto fieldBasePtr = Detail::RFieldBase::Create(fieldDescriptor.GetFieldName(), typeName).Unwrap();
      Detail::RFieldFuse::ConnectRecursively(fieldId, source, *fieldBasePtr);
      return fieldBasePtr;
   }

public:
   RNTupleColumnReader(const std::string &colName, RPageSource &source)
      : fField(MakeField(colName, source)), fValue(fField->GenerateValue()), fLastEntry(-1)
   {
   }

   void *GetImpl(Long64_t entry) final
   {
      if (entry != fLastEntry) {
         fField->Read(entry, &fValue);
         fLastEntry = entry;
      }
      return fValue.GetRawPtr();
   }
};
} // namespace Detail

void RNTupleDS::AddFields(const RNTupleDescriptor &desc, DescriptorId_t parentId)
{
   for (const auto& f : desc.GetFieldRange(parentId)) {
      fColumnNames.emplace_back(desc.GetQualifiedFieldName(f.GetId()));
      fColumnTypes.emplace_back(f.GetTypeName());
      if (f.GetStructure() == ENTupleStructure::kRecord)
         AddFields(desc, f.GetId());
   }
}


RNTupleDS::RNTupleDS(std::unique_ptr<Detail::RPageSource> pageSource)
{
   pageSource->Attach();
   const auto &descriptor = pageSource->GetDescriptor();

   AddFields(descriptor, descriptor.GetFieldZeroId());

   fSources.emplace_back(std::move(pageSource));
}

RDF::RDataSource::Record_t RNTupleDS::GetColumnReadersImpl(std::string_view /* name */, const std::type_info & /* ti */)
{
   // This datasource uses the GetColumnReaders2 API instead (better name in the works)
   return {};
}

std::unique_ptr<ROOT::Detail::RDF::RColumnReaderBase>
RNTupleDS::GetColumnReaders(unsigned int slot, std::string_view name, const std::type_info & /*tid*/)
{
   return std::make_unique<ROOT::Experimental::Detail::RNTupleColumnReader>(std::string(name), *fSources[slot]);
}

bool RNTupleDS::SetEntry(unsigned int, ULong64_t)
{
   return true;
}

std::vector<std::pair<ULong64_t, ULong64_t>> RNTupleDS::GetEntryRanges()
{
   // TODO(jblomer): use cluster boundaries for the entry ranges
   std::vector<std::pair<ULong64_t, ULong64_t>> ranges;
   if (fHasSeenAllRanges) return ranges;

   auto nEntries = fSources[0]->GetNEntries();
   const auto chunkSize = nEntries / fNSlots;
   const auto reminder = 1U == fNSlots ? 0 : nEntries % fNSlots;
   auto start = 0UL;
   auto end = 0UL;
   for (auto i : ROOT::TSeqU(fNSlots)) {
      start = end;
      end += chunkSize;
      ranges.emplace_back(start, end);
      (void)i;
   }
   ranges.back().second += reminder;
   fHasSeenAllRanges = true;
   return ranges;
}


std::string RNTupleDS::GetTypeName(std::string_view colName) const
{
   const auto index = std::distance(
      fColumnNames.begin(), std::find(fColumnNames.begin(), fColumnNames.end(), colName));
   return fColumnTypes[index];
}


bool RNTupleDS::HasColumn(std::string_view colName) const
{
   return std::find(fColumnNames.begin(), fColumnNames.end(), colName) !=
          fColumnNames.end();
}


void RNTupleDS::Initialise()
{
   fHasSeenAllRanges = false;
}


void RNTupleDS::Finalise()
{
}


void RNTupleDS::SetNSlots(unsigned int nSlots)
{
   R__ASSERT(fNSlots == 0);
   R__ASSERT(nSlots > 0);
   fNSlots = nSlots;

   for (unsigned int i = 1; i < fNSlots; ++i) {
      fSources.emplace_back(fSources[0]->Clone());
      R__ASSERT(i == (fSources.size() - 1));
      fSources[i]->Attach();
   }
}
} // ns Experimental
} // ns ROOT


ROOT::RDataFrame ROOT::Experimental::MakeNTupleDataFrame(std::string_view ntupleName, std::string_view fileName)
{
   auto pageSource = ROOT::Experimental::Detail::RPageSource::Create(ntupleName, fileName);
   ROOT::RDataFrame rdf(std::make_unique<RNTupleDS>(std::move(pageSource)));
   return rdf;
}
