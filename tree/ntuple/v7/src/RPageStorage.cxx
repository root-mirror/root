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
#include <ROOT/RPageSinkBuf.hxx>
#include <ROOT/RPageStorageFile.hxx>
#include <ROOT/RStringView.hxx>
#ifdef R__ENABLE_DAOS
# include <ROOT/RPageStorageDaos.hxx>
#endif
#ifdef R__ENABLE_DAVIX
#include <ROOT/RPageStorageS3.hxx>
#endif

#include <Compression.h>
#include <TError.h>

#include <utility>


ROOT::Experimental::Detail::RPageStorage::RPageStorage(std::string_view name) : fNTupleName(name)
{
}

ROOT::Experimental::Detail::RPageStorage::~RPageStorage()
{
}


//------------------------------------------------------------------------------


ROOT::Experimental::Detail::RPageSource::RPageSource(std::string_view name, const RNTupleReadOptions &options)
   : RPageStorage(name), fMetrics(""), fOptions(options)
{
}

ROOT::Experimental::Detail::RPageSource::~RPageSource()
{
}

std::unique_ptr<ROOT::Experimental::Detail::RPageSource> ROOT::Experimental::Detail::RPageSource::Create(
   std::string_view ntupleName, std::string_view location, const RNTupleReadOptions &options)
{
   if (ntupleName.empty()) {
      throw RException(R__FAIL("empty RNTuple name"));
   }
   if (location.empty()) {
      throw RException(R__FAIL("empty storage location"));
   }
   if (location.find("daos://") == 0)
#ifdef R__ENABLE_DAOS
      return std::make_unique<RPageSourceDaos>(ntupleName, location, options);
#else
      throw RException(R__FAIL("This RNTuple build does not support DAOS."));
#endif
   if (location.find("s3://") == 0) {
#ifdef R__ENABLE_DAVIX
      return std::make_unique<RPageSourceS3>(ntupleName, location, options);
#else
      throw RException(R__FAIL("This RNTuple build does not support S3."));
#endif
   }
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

void ROOT::Experimental::Detail::RPageSource::EnableDefaultMetrics(const std::string &prefix)
{
   fMetrics = RNTupleMetrics(prefix);
   fCounters = std::unique_ptr<RCounters>(new RCounters{
      *fMetrics.MakeCounter<RNTupleAtomicCounter*>("nReadV", "", "number of vector read requests"),
      *fMetrics.MakeCounter<RNTupleAtomicCounter*>("nRead", "", "number of byte ranges read"),
      *fMetrics.MakeCounter<RNTupleAtomicCounter*>("szReadPayload", "B", "volume read from storage (required)"),
      *fMetrics.MakeCounter<RNTupleAtomicCounter*>("szReadOverhead", "B", "volume read from storage (overhead)"),
      *fMetrics.MakeCounter<RNTupleAtomicCounter*>("szUnzip", "B", "volume after unzipping"),
      *fMetrics.MakeCounter<RNTupleAtomicCounter*>("nClusterLoaded", "",
                                                   "number of partial clusters preloaded from storage"),
      *fMetrics.MakeCounter<RNTupleAtomicCounter*>("nPageLoaded", "", "number of pages loaded from storage"),
      *fMetrics.MakeCounter<RNTupleAtomicCounter*>("nPagePopulated", "", "number of populated pages"),
      *fMetrics.MakeCounter<RNTupleAtomicCounter*>("timeWallRead", "ns", "wall clock time spent reading"),
      *fMetrics.MakeCounter<RNTupleAtomicCounter*>("timeWallUnzip", "ns", "wall clock time spent decompressing"),
      *fMetrics.MakeCounter<RNTupleTickCounter<RNTupleAtomicCounter>*>("timeCpuRead", "ns", "CPU time spent reading"),
      *fMetrics.MakeCounter<RNTupleTickCounter<RNTupleAtomicCounter>*> ("timeCpuUnzip", "ns",
                                                                        "CPU time spent decompressing"),
      *fMetrics.MakeCounter<RNTupleCalcPerf*> ("bwRead", "MB/s", "bandwidth compressed bytes read per second",
         fMetrics, [](const RNTupleMetrics &metrics) -> std::pair<bool, double> {
            if (const auto szReadPayload = metrics.GetLocalCounter("szReadPayload")) {
               if (const auto szReadOverhead = metrics.GetLocalCounter("szReadOverhead")) {
                  if (const auto timeWallRead = metrics.GetLocalCounter("timeWallRead")) {
                     if (auto walltime = timeWallRead->GetValueAsInt()) {
                        double payload = szReadPayload->GetValueAsInt();
                        double overhead = szReadOverhead->GetValueAsInt();
                        // unit: bytes / nanosecond = GB/s
                        return {true, (1000. * (payload + overhead) / walltime)};
                     }
                  }
               }
            }
            return {false, -1.};
         }
      ),
      *fMetrics.MakeCounter<RNTupleCalcPerf*> ("bwReadUnzip", "MB/s", "bandwidth uncompressed bytes read per second",
         fMetrics, [](const RNTupleMetrics &metrics) -> std::pair<bool, double> {
            if (const auto szUnzip = metrics.GetLocalCounter("szUnzip")) {
               if (const auto timeWallRead = metrics.GetLocalCounter("timeWallRead")) {
                  if (auto walltime = timeWallRead->GetValueAsInt()) {
                     double unzip = szUnzip->GetValueAsInt();
                     // unit: bytes / nanosecond = GB/s
                     return {true, 1000. * unzip / walltime};
                  }
               }
            }
            return {false, -1.};
         }
      ),
      *fMetrics.MakeCounter<RNTupleCalcPerf*> ("bwUnzip", "MB/s", "decompression bandwidth of uncompressed bytes per second",
         fMetrics, [](const RNTupleMetrics &metrics) -> std::pair<bool, double> {
            if (const auto szUnzip = metrics.GetLocalCounter("szUnzip")) {
               if (const auto timeWallUnzip = metrics.GetLocalCounter("timeWallUnzip")) {
                  if (auto walltime = timeWallUnzip->GetValueAsInt()) {
                     double unzip = szUnzip->GetValueAsInt();
                     // unit: bytes / nanosecond = GB/s
                     return {true, 1000. * unzip / walltime};
                  }
               }
            }
            return {false, -1.};
         }
      ),
      *fMetrics.MakeCounter<RNTupleCalcPerf*> ("rtReadEfficiency", "", "ratio of payload over all bytes read",
         fMetrics, [](const RNTupleMetrics &metrics) -> std::pair<bool, double> {
            if (const auto szReadPayload = metrics.GetLocalCounter("szReadPayload")) {
               if (const auto szReadOverhead = metrics.GetLocalCounter("szReadOverhead")) {
                  if (auto payload = szReadPayload->GetValueAsInt()) {
                     // r/(r+o) = 1/((r+o)/r) = 1/(1 + o/r)
                     return {true, 1./(1. + (1. * szReadOverhead->GetValueAsInt()) / payload)};
                  }
               }
            }
            return {false, -1.};
         }
      ),
      *fMetrics.MakeCounter<RNTupleCalcPerf*> ("rtCompression", "", "ratio of compressed bytes / uncompressed bytes",
         fMetrics, [](const RNTupleMetrics &metrics) -> std::pair<bool, double> {
            if (const auto szReadPayload = metrics.GetLocalCounter("szReadPayload")) {
               if (const auto szUnzip = metrics.GetLocalCounter("szUnzip")) {
                  if (auto unzip = szUnzip->GetValueAsInt()) {
                     return {true, (1. * szReadPayload->GetValueAsInt()) / unzip};
                  }
               }
            }
            return {false, -1.};
         }
      )
   });
}


//------------------------------------------------------------------------------


ROOT::Experimental::Detail::RPageSink::RPageSink(std::string_view name, const RNTupleWriteOptions &options)
   : RPageStorage(name), fMetrics(""), fOptions(options.Clone())
{
}

ROOT::Experimental::Detail::RPageSink::~RPageSink()
{
}

std::unique_ptr<ROOT::Experimental::Detail::RPageSink> ROOT::Experimental::Detail::RPageSink::Create(
   std::string_view ntupleName, std::string_view location, const RNTupleWriteOptions &options)
{
   if (ntupleName.empty()) {
      throw RException(R__FAIL("empty RNTuple name"));
   }
   if (location.empty()) {
      throw RException(R__FAIL("empty storage location"));
   }
   std::unique_ptr<ROOT::Experimental::Detail::RPageSink> realSink;
   if (location.find("daos://") == 0) {
#ifdef R__ENABLE_DAOS
      realSink = std::make_unique<RPageSinkDaos>(ntupleName, location, options);
#else
      throw RException(R__FAIL("This RNTuple build does not support DAOS."));
#endif
   } else if (location.find("s3://") == 0) {
#ifdef R__ENABLE_DAVIX
      return std::make_unique<RPageSinkS3>(ntupleName, location, options);
#else
      throw RException(R__FAIL("This RNTuple build does not support S3."));
#endif
   } else {
      realSink = std::make_unique<RPageSinkFile>(ntupleName, location, options);
   }

   if (options.GetUseBufferedWrite())
      return std::make_unique<RPageSinkBuf>(std::move(realSink));
   return realSink;
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
      f.ConnectPageSink(*this); // issues in turn one or several calls to AddColumn()
   }

   auto nColumns = fLastColumnId;
   for (DescriptorId_t i = 0; i < nColumns; ++i) {
      RClusterDescriptor::RColumnRange columnRange;
      columnRange.fColumnId = i;
      columnRange.fFirstElementIndex = 0;
      columnRange.fNElements = 0;
      columnRange.fCompressionSettings = GetWriteOptions().GetCompression();
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
ROOT::Experimental::Detail::RPageSink::SealPage(const RPage &page,
   const RColumnElementBase &element, int compressionSetting, void *buf)
{
   unsigned char *pageBuf = reinterpret_cast<unsigned char *>(page.GetBuffer());
   bool isAdoptedBuffer = true;
   auto packedBytes = page.GetSize();

   if (!element.IsMappable()) {
      packedBytes = element.GetPackedSize(page.GetNElements());
      pageBuf = new unsigned char[packedBytes];
      isAdoptedBuffer = false;
      element.Pack(pageBuf, page.GetBuffer(), page.GetNElements());
   }
   auto zippedBytes = packedBytes;

   if ((compressionSetting != 0) || !element.IsMappable()) {
      zippedBytes = RNTupleCompressor::Zip(pageBuf, packedBytes, compressionSetting, buf);
      if (!isAdoptedBuffer)
         delete[] pageBuf;
      pageBuf = reinterpret_cast<unsigned char *>(buf);
      isAdoptedBuffer = true;
   }

   R__ASSERT(isAdoptedBuffer);

   return RSealedPage{pageBuf, zippedBytes, page.GetNElements()};
}

ROOT::Experimental::Detail::RPageStorage::RSealedPage
ROOT::Experimental::Detail::RPageSink::SealPage(
   const RPage &page, const RColumnElementBase &element, int compressionSetting)
{
   R__ASSERT(fCompressor);
   return SealPage(page, element, compressionSetting, fCompressor->GetZipBuffer());
}

void ROOT::Experimental::Detail::RPageSink::EnableDefaultMetrics(const std::string &prefix)
{
   fMetrics = RNTupleMetrics(prefix);
   fCounters = std::unique_ptr<RCounters>(new RCounters{
      *fMetrics.MakeCounter<RNTupleAtomicCounter*>("nPageCommitted", "", "number of pages committed to storage"),
      *fMetrics.MakeCounter<RNTupleAtomicCounter*>("szWritePayload", "B", "volume written for committed pages"),
      *fMetrics.MakeCounter<RNTupleAtomicCounter*>("szZip", "B", "volume before zipping"),
      *fMetrics.MakeCounter<RNTupleAtomicCounter*>("timeWallWrite", "ns", "wall clock time spent writing"),
      *fMetrics.MakeCounter<RNTupleAtomicCounter*>("timeWallZip", "ns", "wall clock time spent compressing"),
      *fMetrics.MakeCounter<RNTupleTickCounter<RNTupleAtomicCounter>*>("timeCpuWrite", "ns", "CPU time spent writing"),
      *fMetrics.MakeCounter<RNTupleTickCounter<RNTupleAtomicCounter>*> ("timeCpuZip", "ns",
                                                                        "CPU time spent compressing")
   });
}
