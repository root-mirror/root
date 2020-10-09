/// \file RNTupleDS.cxx
/// \ingroup NTuple ROOT7
/// \author Jakob Blomer <jblomer@cern.ch>
/// \date 2018-10-04
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

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

void RNTupleDS::AddFields(const RNTupleDescriptor &desc, DescriptorId_t parentId)
{
   for (const auto& f : desc.GetFieldRange(parentId)) {
      fColumnNames.emplace_back(desc.GetQualifiedFieldName(f.GetId()));
      fColumnTypes.emplace_back(f.GetTypeName());
      std::string normalized;
      TClassEdit::GetNormalizedName(normalized, f.GetTypeName());
      fNormalizedColumnTypes.emplace_back(normalized);
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


RDF::RDataSource::Record_t RNTupleDS::GetColumnReadersImpl(std::string_view name, const std::type_info& ti)
{
   const auto colIdx = std::distance(
      fColumnNames.begin(), std::find(fColumnNames.begin(), fColumnNames.end(), name));

   std::string demangled = ROOT::Internal::RDF::DemangleTypeIdName(ti);
   std::string normalized;
   TClassEdit::GetNormalizedName(normalized, demangled.c_str());
   if (normalized != fNormalizedColumnTypes[colIdx]) {
      std::string err = "The type of column \"";
      err += name;
      err += "\" is ";
      err += fColumnTypes[colIdx];
      if (fColumnTypes[colIdx] != fNormalizedColumnTypes[colIdx])
         err += " (= " + fNormalizedColumnTypes[colIdx] + ")";
      err += " but ";
      err += demangled;
      if (demangled != normalized)
         err += " (= " + normalized + ")";
      err += " has been selected";
      throw std::runtime_error(err);
   }

   std::vector<void*> ptrs;
   for (unsigned int slot = 0; slot < fNSlots; ++slot) {
      if (!fValuePtrs[slot][colIdx]) {
         const auto &descriptor = fSources[slot]->GetDescriptor();
         auto colName = fColumnNames[colIdx];
         auto typeName = fColumnTypes[colIdx];
         auto fieldId = descriptor.FindFieldId(colName);
         fFields[slot][colIdx] = std::unique_ptr<Detail::RFieldBase>(
            Detail::RFieldBase::Create(descriptor.GetFieldDescriptor(fieldId).GetFieldName(), typeName));
         Detail::RFieldFuse::ConnectRecursively(fieldId, *fSources[slot], *fFields[slot][colIdx]);
         fValues[slot][colIdx] = fFields[slot][colIdx]->GenerateValue();
         fValuePtrs[slot][colIdx] = fValues[slot][colIdx].GetRawPtr();
         if (slot == 0)
            fActiveColumns.emplace_back(colIdx);
      }
      ptrs.push_back(&fValuePtrs[slot][colIdx]);
   }

   return ptrs;
}

bool RNTupleDS::SetEntry(unsigned int slot, ULong64_t entryIndex)
{
   for (auto colIdx : fActiveColumns) {
      fFields[slot][colIdx]->Read(entryIndex, &fValues[slot][colIdx]);
   }
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

   fFields.resize(fNSlots);
   fValues.resize(fNSlots);
   fValuePtrs.resize(fNSlots);
   for (unsigned int i = 1; i < fNSlots; ++i) {
      fSources.emplace_back(fSources[0]->Clone());
      R__ASSERT(i == (fSources.size() - 1));
      fSources[i]->Attach();
   }

   auto nColumns = fColumnNames.size();
   for (unsigned int i = 0; i < fNSlots; ++i) {
      fFields[i].resize(nColumns);
      fValues[i].resize(nColumns);
      fValuePtrs[i].resize(nColumns, nullptr);
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
