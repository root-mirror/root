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


std::unique_ptr<unsigned char []> ROOT::Experimental::Detail::RPageSource::UnsealPage(
   const RSealedPage &sealedPage, const RColumnElementBase &element)
{
   const auto bytesPacked = element.GetPackedSize(sealedPage.fNElements);
   const auto pageSize = element.GetSize() * sealedPage.fNElements;

   // TODO(jblomer): We might be able to do better memory handling for unsealing pages than a new malloc for every
   // new page.
   auto pageBuffer = std::make_unique<unsigned char[]>(bytesPacked);
   if (sealedPage.fSize != bytesPacked) {
      fDecompressor->Unzip(sealedPage.fBuffer, sealedPage.fSize, bytesPacked, pageBuffer.get());
   } else {
      // We cannot simply map the sealed page as we don't know its life time. Specialized page sources
      // may decide to implement to not use UnsealPage but to custom mapping / decompression code.
      // Note that usually pages are compressed.
      memcpy(pageBuffer.get(), sealedPage.fBuffer, bytesPacked);
   }

   if (!element.IsMappable()) {
      auto unpackedBuffer = new unsigned char[pageSize];
      element.Unpack(unpackedBuffer, pageBuffer.get(), sealedPage.fNElements);
      pageBuffer = std::unique_ptr<unsigned char []>(unpackedBuffer);
   }

   return pageBuffer;
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

   auto &fieldZero = *model.GetFieldZero();
   fDescriptorBuilder.AddField(
      RDanglingFieldDescriptor::FromField(fieldZero)
         .FieldId(fLastFieldId)
         .MakeDescriptor()
         .Unwrap()
   );
   fieldZero.SetOnDiskId(fLastFieldId);
   for (auto& f : *model.GetFieldZero()) {
      fLastFieldId++;
      fDescriptorBuilder.AddField(
         RDanglingFieldDescriptor::FromField(f)
            .FieldId(fLastFieldId)
            .MakeDescriptor()
            .Unwrap()
      );
      fDescriptorBuilder.AddFieldLink(f.GetParent()->GetOnDiskId(), fLastFieldId);
      f.SetOnDiskId(fLastFieldId);
      f.ConnectPageStorage(*this); // issues in turn one or several calls to AddColumn()
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
   fOpenColumnRanges.at(columnHandle.fId).fNElements += page.GetNElements();

   RClusterDescriptor::RPageRange::RPageInfo pageInfo;
   pageInfo.fNElements = page.GetNElements();
   pageInfo.fLocator = CommitPageImpl(columnHandle, page);
   fOpenPageRanges.at(columnHandle.fId).fPageInfos.emplace_back(pageInfo);
}


void ROOT::Experimental::Detail::RPageSink::CommitSealedPage(
   ROOT::Experimental::DescriptorId_t columnId,
   const ROOT::Experimental::Detail::RPageStorage::RSealedPage &sealedPage)
{
   fOpenColumnRanges.at(columnId).fNElements += sealedPage.fNElements;

   RClusterDescriptor::RPageRange::RPageInfo pageInfo;
   pageInfo.fNElements = sealedPage.fNElements;
   pageInfo.fLocator = CommitSealedPageImpl(columnId, sealedPage);
   fOpenPageRanges.at(columnId).fPageInfos.emplace_back(pageInfo);
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


ROOT::Experimental::Detail::RPageStorage::RSealedPage
ROOT::Experimental::Detail::RPageSink::SealPage(
   const RPage &page, const RColumnElementBase &element, int compressionSetting)
{
   unsigned char *buffer = reinterpret_cast<unsigned char *>(page.GetBuffer());
   bool isAdoptedBuffer = true;
   auto packedBytes = page.GetSize();

   if (!element.IsMappable()) {
      packedBytes = element.GetPackedSize(page.GetNElements());
      buffer = new unsigned char[packedBytes];
      isAdoptedBuffer = false;
      element.Pack(buffer, page.GetBuffer(), page.GetNElements());
   }
   auto zippedBytes = packedBytes;

   if ((compressionSetting != 0) || !element.IsMappable()) {
      zippedBytes = fCompressor->Zip(buffer, packedBytes, fOptions.GetCompression());
      if (!isAdoptedBuffer)
         delete[] buffer;
      buffer = const_cast<unsigned char *>(reinterpret_cast<const unsigned char *>(fCompressor->GetZipBuffer()));
      isAdoptedBuffer = true;
   }

   R__ASSERT(isAdoptedBuffer);

   return RSealedPage{buffer, zippedBytes, page.GetNElements()};
}
