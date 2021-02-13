/// \file RPageStorage.cxx
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

#include <ROOT/RPageStorage.hxx>
#include <ROOT/RPageStorageFile.hxx>
#include <ROOT/RColumn.hxx>
#include <ROOT/RField.hxx>
#include <ROOT/RNTupleDescriptor.hxx>
#include <ROOT/RNTupleMetrics.hxx>
#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RPagePool.hxx>
#include <ROOT/RPageStorageFile.hxx>
#include <ROOT/RStringView.hxx>

#include <Compression.h>
#include <TError.h>

#include <unordered_map>
#include <utility>


ROOT::Experimental::Detail::RPageStorage::RPageStorage(std::string_view name) : fNTupleName(name)
{
}

ROOT::Experimental::Detail::RPageStorage::~RPageStorage()
{
}

ROOT::Experimental::Detail::RNTupleMetrics &ROOT::Experimental::Detail::RPageStorage::GetMetrics()
{
   static RNTupleMetrics metrics("");
   return metrics;
}


//------------------------------------------------------------------------------


ROOT::Experimental::Detail::RPageSource::RPageSource(std::string_view name, const RNTupleReadOptions &options)
   : RPageStorage(name), fOptions(options)
{
}

ROOT::Experimental::Detail::RPageSource::~RPageSource()
{
}

std::unique_ptr<ROOT::Experimental::Detail::RPageSource> ROOT::Experimental::Detail::RPageSource::Create(
   std::string_view ntupleName, std::string_view location, const RNTupleReadOptions &options)
{
   return std::make_unique<RPageSourceFile>(ntupleName, location, options);
}

ROOT::Experimental::Detail::RPageStorage::ColumnHandle_t
ROOT::Experimental::Detail::RPageSource::AddColumn(DescriptorId_t fieldId, const RColumn &column)
{
   R__ASSERT(fieldId != kInvalidDescriptorId);
   auto columnId = fDescriptor.FindColumnId(fieldId, column.GetIndex());
   R__ASSERT(columnId != kInvalidDescriptorId);
   fActiveColumns.emplace(columnId);
   return ColumnHandle_t{columnId, &column};
}

void ROOT::Experimental::Detail::RPageSource::DropColumn(ColumnHandle_t columnHandle)
{
   fActiveColumns.erase(columnHandle.fId);
}

ROOT::Experimental::NTupleSize_t ROOT::Experimental::Detail::RPageSource::GetNEntries()
{
   return fDescriptor.GetNEntries();
}

ROOT::Experimental::NTupleSize_t ROOT::Experimental::Detail::RPageSource::GetNElements(ColumnHandle_t columnHandle)
{
   return fDescriptor.GetNElements(columnHandle.fId);
}

ROOT::Experimental::ColumnId_t ROOT::Experimental::Detail::RPageSource::GetColumnId(ColumnHandle_t columnHandle)
{
   // TODO(jblomer) distinguish trees
   return columnHandle.fId;
}

void ROOT::Experimental::Detail::RPageSource::UnzipCluster(RCluster *cluster)
{
   if (fTaskScheduler)
      UnzipClusterImpl(cluster);
}


//------------------------------------------------------------------------------


ROOT::Experimental::Detail::RPageSink::RPageSink(std::string_view name, const RNTupleWriteOptions &options)
   : RPageStorage(name), fOptions(options)
{
}

ROOT::Experimental::Detail::RPageSink::~RPageSink()
{
}

std::unique_ptr<ROOT::Experimental::Detail::RPageSink> ROOT::Experimental::Detail::RPageSink::Create(
   std::string_view ntupleName, std::string_view location, const RNTupleWriteOptions &options)
{
   return std::make_unique<RPageSinkFile>(ntupleName, location, options);
}

ROOT::Experimental::Detail::RPageStorage::ColumnHandle_t
ROOT::Experimental::Detail::RPageSink::AddColumn(DescriptorId_t fieldId, const RColumn &column)
{
   auto columnId = fLastColumnId++;
   fDescriptorBuilder.AddColumn(columnId, fieldId, column.GetVersion(), column.GetModel(), column.GetIndex());
   return ColumnHandle_t{columnId, &column};
}


void ROOT::Experimental::Detail::RPageSink::Create(RNTupleModel &model)
{
   fDescriptorBuilder.SetNTuple(fNTupleName, model.GetDescription(), "undefined author",
                                model.GetVersion(), model.GetUuid());

   std::unordered_map<const RFieldBase *, DescriptorId_t> fieldPtr2Id; // necessary to find parent field ids
   const auto &fieldZero = *model.GetFieldZero();
   fDescriptorBuilder.AddField(
      RDanglingFieldDescriptor::FromField(fieldZero)
         .FieldId(fLastFieldId)
         .MakeDescriptor()
         .Unwrap()
   );
   fieldPtr2Id[&fieldZero] = fLastFieldId++;
   for (auto& f : *model.GetFieldZero()) {
      fDescriptorBuilder.AddField(
         RDanglingFieldDescriptor::FromField(f)
            .FieldId(fLastFieldId)
            .MakeDescriptor()
            .Unwrap()
      );
      fDescriptorBuilder.AddFieldLink(fieldPtr2Id[f.GetParent()], fLastFieldId);
      Detail::RFieldFuse::Connect(fLastFieldId, *this, f); // issues in turn one or several calls to AddColumn()
      fieldPtr2Id[&f] = fLastFieldId++;
   }

   auto nColumns = fLastColumnId;
   for (DescriptorId_t i = 0; i < nColumns; ++i) {
      RClusterDescriptor::RColumnRange columnRange;
      columnRange.fColumnId = i;
      columnRange.fFirstElementIndex = 0;
      columnRange.fNElements = 0;
      columnRange.fCompressionSettings = fOptions.GetCompression();
      fOpenColumnRanges.emplace_back(columnRange);
      RClusterDescriptor::RPageRange pageRange;
      pageRange.fColumnId = i;
      fOpenPageRanges.emplace_back(std::move(pageRange));
   }

   CreateImpl(model);
}


void ROOT::Experimental::Detail::RPageSink::CommitPage(ColumnHandle_t columnHandle, const RPage &page)
{
   auto locator = CommitPageImpl(columnHandle, page);

   auto columnId = columnHandle.fId;
   fOpenColumnRanges[columnId].fNElements += page.GetNElements();
   RClusterDescriptor::RPageRange::RPageInfo pageInfo;
   pageInfo.fNElements = page.GetNElements();
   pageInfo.fLocator = locator;
   fOpenPageRanges[columnId].fPageInfos.emplace_back(pageInfo);
}


void ROOT::Experimental::Detail::RPageSink::CommitCluster(ROOT::Experimental::NTupleSize_t nEntries)
{
   auto locator = CommitClusterImpl(nEntries);

   R__ASSERT((nEntries - fPrevClusterNEntries) < ClusterSize_t(-1));
   fDescriptorBuilder.AddCluster(fLastClusterId, RNTupleVersion(), fPrevClusterNEntries,
                                 ClusterSize_t(nEntries - fPrevClusterNEntries));
   fDescriptorBuilder.SetClusterLocator(fLastClusterId, locator);
   for (auto &range : fOpenColumnRanges) {
      fDescriptorBuilder.AddClusterColumnRange(fLastClusterId, range);
      range.fFirstElementIndex += range.fNElements;
      range.fNElements = 0;
   }
   for (auto &range : fOpenPageRanges) {
      RClusterDescriptor::RPageRange fullRange;
      std::swap(fullRange, range);
      range.fColumnId = fullRange.fColumnId;
      fDescriptorBuilder.AddClusterPageRange(fLastClusterId, std::move(fullRange));
   }
   ++fLastClusterId;
   fPrevClusterNEntries = nEntries;
}
