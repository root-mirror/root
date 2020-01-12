/// \file ROOT/RPageStorageChain.cxx
/// \ingroup NTuple ROOT7
/// \author Simon Leisibach <simon.leisibach@gmail.com>
/// \date 2019-09-09
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/RNTupleOptions.hxx>
#include <ROOT/RLogger.hxx>
#include <ROOT/RPage.hxx>
#include <ROOT/RPageStorage.hxx>
#include <ROOT/RPageStorageChain.hxx>
#include <ROOT/RPageStorageRaw.hxx>
#include <ROOT/RPageStorageRoot.hxx>

#include <TError.h>

#include <iostream>
#include <memory>
#include <unordered_map>

ROOT::Experimental::Detail::RPageSourceChain::RPageSourceChain(std::string_view ntupleName,
                                                               std::vector<std::string> locationVec,
                                                               const RNTupleReadOptions &options)
   : RPageSource{ntupleName, options}, fMetrics{"RPageSourceChain"}
{
   // No need to check if locationVec is empty. It's already tested in RNTupleReader::Open

   for (const auto &location : locationVec) {
      fSources.emplace_back(Create(ntupleName, location, options));
      fSources.back()->Attach();
   }

   fUnsafe = CompareFileMetaData();
   InitializeMemberVariables();
}

ROOT::Experimental::Detail::RPageSourceChain::RPageSourceChain(std::string_view ntupleName,
                                                               std::vector<RPageSource *> sources,
                                                               const RNTupleReadOptions &options)
   : RPageSource{ntupleName, options}, fMetrics{"RPageSourceChain"}
{
   for (const auto &s : sources) {
      fSources.emplace_back(s->Clone());
      fSources.back()->Attach();
   }

   fUnsafe = CompareFileMetaData();
   InitializeMemberVariables();
}

ROOT::Experimental::Detail::RPageSourceChain::RPageSourceChain(std::string_view ntupleName,
                                                               std::vector<std::unique_ptr<RPageSource>> &&sources,
                                                               const RNTupleReadOptions &options)
   : RPageSource{ntupleName, options}, fSources{std::move(sources)}, fMetrics{"RPageSourceChain"}
{
   fUnsafe = CompareFileMetaData();
   InitializeMemberVariables();
}

bool ROOT::Experimental::Detail::RPageSourceChain::CompareFileMetaData()
{
   for (std::size_t i = 1; i < fSources.size(); ++i) {
      // checks only number of fields and columns
      if ((fSources.at(0)->GetDescriptor().GetNFields() != fSources.at(i)->GetDescriptor().GetNFields()) ||
          (fSources.at(0)->GetDescriptor().GetNColumns() != fSources.at(i)->GetDescriptor().GetNColumns())) {
         R__WARNING_HERE("NTuple")
            << "The meta-data (number of fields and columns) of the files do not match. A nullptr was returned. Please "
               "make sure the number of fields and columns in all files are identical!";
         return true;
      }
      // compares all fieldDescriptors
      for (std::size_t j = 0; j < fSources.at(0)->GetDescriptor().GetNFields(); ++j) {
         if (!(fSources.at(0)->GetDescriptor().GetFieldDescriptor(j) ==
               fSources.at(i)->GetDescriptor().GetFieldDescriptor(j))) {
            R__WARNING_HERE("NTuple")
               << "The meta-data of the fields of the files do not match. A nullptr was returned. Please make sure the "
                  "metadata of the fields (fieldName, field order, etc.) is the same across all files!";
            return true;
         }
      }
      // compares all columnDescriptors
      for (std::size_t j = 0; j < fSources.at(0)->GetDescriptor().GetNColumns(); ++j) {
         if (!(fSources.at(0)->GetDescriptor().GetColumnDescriptor(j) ==
               fSources.at(i)->GetDescriptor().GetColumnDescriptor(j))) {
            R__WARNING_HERE("NTuple")
               << "The meta-data of the columns of the files do not match. A nullptr was returned. Please make sure "
                  "the metadata of columns are the same across all files!";
            return true;
         }
      }
   }
   return false;
}

void ROOT::Experimental::Detail::RPageSourceChain::InitializeMemberVariables()
{
   // initialize fNEntryPerSource and fNClusterPerSource
   fNEntryPerSource.resize(fSources.size() + 1);
   fNClusterPerSource.resize(fSources.size() + 1);
   fNEntryPerSource.at(0) = 0;
   fNClusterPerSource.at(0) = 0;
   for (std::size_t i = 0; i < fSources.size(); ++i) {
      fNEntryPerSource.at(i + 1) = fNEntryPerSource.at(i) + fSources.at(i)->GetNEntries();
      fNClusterPerSource.at(i + 1) = fNClusterPerSource.at(i) + fSources.at(i)->GetDescriptor().GetNClusters();
   }

   // initialize fNElementsPerColumnPerSource
   auto nColumns{fSources.at(0)->GetDescriptor().GetNColumns()};
   fNElementsPerColumnPerSource.resize(fSources.size() + 1);
   for (auto &fN : fNElementsPerColumnPerSource) {
      fN.resize(nColumns);
   }
   for (std::size_t i = 0; i < nColumns; ++i) {
      fNElementsPerColumnPerSource.at(0).at(i) = 0;
   }
   for (std::size_t i = 0; i < fSources.size(); ++i) {
      for (std::size_t j = 0; j < nColumns; ++j) {
         auto lastClusterId = fSources.at(i)->GetDescriptor().GetNClusters() - 1;
         auto columnDescript = fSources.at(i)->GetDescriptor().GetClusterDescriptor(lastClusterId).GetColumnRange(j);
         fNElementsPerColumnPerSource.at(i + 1).at(j) =
            columnDescript.fFirstElementIndex + columnDescript.fNElements + fNElementsPerColumnPerSource.at(i).at(j);
      }
   }
}

ROOT::Experimental::RNTupleDescriptor ROOT::Experimental::Detail::RPageSourceChain::DoAttach()
{
   RNTupleDescriptorBuilder descBuilder;
   descBuilder.SetNTuple(fSources.at(0)->GetDescriptor());
   descBuilder.AddFieldsAndColumnsFromDescriptor(fSources.at(0)->GetDescriptor());

   for (std::size_t i = 0; i < fSources.size(); ++i) {
      descBuilder.AddClustersFromDescriptor(fSources.at(i)->GetDescriptor());
   }
   return descBuilder.MoveDescriptor();
}

std::unique_ptr<ROOT::Experimental::Detail::RPageSource> ROOT::Experimental::Detail::RPageSourceChain::Clone() const
{
   std::vector<RPageSource *> sourceVec;
   for (const auto &s : fSources) {
      sourceVec.emplace_back(s.get());
   }
   return std::make_unique<RPageSourceChain>(fNTupleName, sourceVec, fOptions);
}

ROOT::Experimental::Detail::RPage
ROOT::Experimental::Detail::RPageSourceChain::PopulatePage(ColumnHandle_t columnHandle, NTupleSize_t globalIndex)
{
   std::size_t sourceIndex{0};
   for (; sourceIndex < fSources.size(); ++sourceIndex) {
      if (globalIndex < fNEntryPerSource.at(sourceIndex + 1))
         break;
   }
   R__ASSERT(sourceIndex != fSources.size() && "globalIndex is bigger than total number of entries.");

   RPage page{fSources.at(sourceIndex)->PopulatePage(columnHandle, globalIndex - fNEntryPerSource.at(sourceIndex))};

   if (fPageMapper.find(page.GetBuffer()) == fPageMapper.end()) {
      fPageMapper.emplace(page.GetBuffer(), PageInfoChain{sourceIndex, 1});
   } else {
      fPageMapper.at(page.GetBuffer()).fNSamePagePopulated += 1;
   }

   auto clusterId = fDescriptor.FindClusterId(columnHandle.fId, globalIndex);
   auto selfOffset = fDescriptor.GetClusterDescriptor(clusterId).GetColumnRange(columnHandle.fId).fFirstElementIndex;
   auto clusterInfo = RPage::RClusterInfo(clusterId, selfOffset);

   page.SetWindow(page.GetGlobalRangeFirst() + fNElementsPerColumnPerSource.at(sourceIndex).at(columnHandle.fId),
                  clusterInfo);
   return page;
}

ROOT::Experimental::Detail::RPage
ROOT::Experimental::Detail::RPageSourceChain::PopulatePage(ColumnHandle_t columnHandle,
                                                           const RClusterIndex &clusterIndex)
{
   auto clusterId = clusterIndex.GetClusterId();
   std::size_t sourceIndex{0};
   for (; sourceIndex < fSources.size(); ++sourceIndex) {
      if (clusterId < fNClusterPerSource.at(sourceIndex + 1))
         break;
   }
   R__ASSERT(sourceIndex != fSources.size() && "clusterIndex is bigger than total number of entries.");

   RClusterIndex newClusterIndex(clusterIndex.GetClusterId() - fNClusterPerSource.at(sourceIndex),
                                 clusterIndex.GetIndex());
   RPage page{fSources.at(sourceIndex)->PopulatePage(columnHandle, newClusterIndex)};

   if (fPageMapper.find(page.GetBuffer()) == fPageMapper.end()) {
      fPageMapper.emplace(page.GetBuffer(), PageInfoChain{sourceIndex, 1});
   } else {
      fPageMapper.at(page.GetBuffer()).fNSamePagePopulated += 1;
   }

   auto selfOffset = fDescriptor.GetClusterDescriptor(clusterId).GetColumnRange(columnHandle.fId).fFirstElementIndex;
   auto clusterInfo = RPage::RClusterInfo(clusterId, selfOffset);

   page.SetWindow(page.GetGlobalRangeFirst() + fNElementsPerColumnPerSource.at(sourceIndex).at(columnHandle.fId),
                  clusterInfo);
   return page;
}

void ROOT::Experimental::Detail::RPageSourceChain::ReleasePage(RPage &page)
{
   if (page.IsNull())
      return;
   auto mapIterator = fPageMapper.find(page.GetBuffer());
   if (mapIterator == fPageMapper.end())
      R__ASSERT(false && "Page could not be assigned to source and released.");

   fSources.at(mapIterator->second.fSourceId)->ReleasePage(page);
   // Sometimes malloc allocates the same memory location as a previous malloc.
   // To avoid a collision of keys for such cases, the entry is deleted from std::unordered_map when released.
   if (mapIterator->second.fNSamePagePopulated-- == 1)
      fPageMapper.erase(mapIterator);
}
