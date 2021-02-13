/// \file RPageStorageFile.cxx
/// \ingroup NTuple ROOT7
/// \author Jakob Blomer <jblomer@cern.ch>
/// \date 2019-11-25
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/RCluster.hxx>
#include <ROOT/RClusterPool.hxx>
#include <ROOT/RField.hxx>
#include <ROOT/RLogger.hxx>
#include <ROOT/RNTupleDescriptor.hxx>
#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RNTupleZip.hxx>
#include <ROOT/RPage.hxx>
#include <ROOT/RPageAllocator.hxx>
#include <ROOT/RPagePool.hxx>
#include <ROOT/RPageStorageFile.hxx>
#include <ROOT/RRawFile.hxx>

#include <RVersion.h>
#include <TError.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <utility>

#include <atomic>
#include <condition_variable>
#include <functional>
#include <mutex>
#include <thread>
#include <queue>

ROOT::Experimental::Detail::RPageSinkFile::RPageSinkFile(std::string_view ntupleName, std::string_view path,
   const RNTupleWriteOptions &options)
   : RPageSink(ntupleName, options)
   , fMetrics("RPageSinkRoot")
   , fPageAllocator(std::make_unique<RPageAllocatorHeap>())
{
   R__LOG_WARNING(NTupleLog()) << "The RNTuple file format will change. " <<
      "Do not store real data with this version of RNTuple!";

   fWriter = std::unique_ptr<Internal::RNTupleFileWriter>(Internal::RNTupleFileWriter::Recreate(
      ntupleName, path, options.GetCompression(), options.GetContainerFormat()));
}


ROOT::Experimental::Detail::RPageSinkFile::RPageSinkFile(std::string_view ntupleName, TFile &file,
   const RNTupleWriteOptions &options)
   : RPageSink(ntupleName, options)
   , fMetrics("RPageSinkRoot")
   , fPageAllocator(std::make_unique<RPageAllocatorHeap>())
{
   R__LOG_WARNING(NTupleLog()) << "The RNTuple file format will change. " <<
      "Do not store real data with this version of RNTuple!";

   fWriter = std::unique_ptr<Internal::RNTupleFileWriter>(Internal::RNTupleFileWriter::Append(ntupleName, file));
}


ROOT::Experimental::Detail::RPageSinkFile::RPageSinkFile(std::string_view ntupleName, std::string_view path,
   const RNTupleWriteOptions &options, std::unique_ptr<TFile> &file)
   : RPageSink(ntupleName, options)
   , fMetrics("RPageSinkRoot")
   , fPageAllocator(std::make_unique<RPageAllocatorHeap>())
{
   R__LOG_WARNING(NTupleLog()) << "The RNTuple file format will change. " <<
      "Do not store real data with this version of RNTuple!";
   fWriter = std::unique_ptr<Internal::RNTupleFileWriter>(
      Internal::RNTupleFileWriter::Recreate(ntupleName, path, file));
}


ROOT::Experimental::Detail::RPageSinkFile::~RPageSinkFile()
{
}


void ROOT::Experimental::Detail::RPageSinkFile::CreateImpl(const RNTupleModel & /* model */)
{
   const auto &descriptor = fDescriptorBuilder.GetDescriptor();
   auto szHeader = descriptor.GetHeaderSize();
   auto buffer = std::make_unique<unsigned char[]>(szHeader);
   descriptor.SerializeHeader(buffer.get());

   auto zipBuffer = std::make_unique<unsigned char[]>(szHeader);
   auto szZipHeader = fCompressor(buffer.get(), szHeader, fOptions.GetCompression(),
      [&zipBuffer](const void *b, size_t n, size_t o){ memcpy(zipBuffer.get() + o, b, n); } );
   fWriter->WriteNTupleHeader(zipBuffer.get(), szZipHeader, szHeader);
}


ROOT::Experimental::RClusterDescriptor::RLocator
ROOT::Experimental::Detail::RPageSinkFile::CommitPageImpl(ColumnHandle_t columnHandle, const RPage &page)
{
   unsigned char *buffer = reinterpret_cast<unsigned char *>(page.GetBuffer());
   bool isAdoptedBuffer = true;
   auto packedBytes = page.GetSize();
   auto element = columnHandle.fColumn->GetElement();
   const auto isMappable = element->IsMappable();

   if (!isMappable) {
      packedBytes = (page.GetNElements() * element->GetBitsOnStorage() + 7) / 8;
      buffer = new unsigned char[packedBytes];
      isAdoptedBuffer = false;
      element->Pack(buffer, page.GetBuffer(), page.GetNElements());
   }
   auto zippedBytes = packedBytes;

   if (fOptions.GetCompression() != 0) {
      zippedBytes = fCompressor(buffer, packedBytes, fOptions.GetCompression());
      if (!isAdoptedBuffer)
         delete[] buffer;
      buffer = const_cast<unsigned char *>(reinterpret_cast<const unsigned char *>(fCompressor.GetZipBuffer()));
      isAdoptedBuffer = true;
   }

   auto offsetData = fWriter->WriteBlob(buffer, zippedBytes, packedBytes);
   fClusterMinOffset = std::min(offsetData, fClusterMinOffset);
   fClusterMaxOffset = std::max(offsetData + zippedBytes, fClusterMaxOffset);

   if (!isAdoptedBuffer)
      delete[] buffer;

   RClusterDescriptor::RLocator result;
   result.fPosition = offsetData;
   result.fBytesOnStorage = zippedBytes;
   return result;
}


ROOT::Experimental::RClusterDescriptor::RLocator
ROOT::Experimental::Detail::RPageSinkFile::CommitClusterImpl(ROOT::Experimental::NTupleSize_t /* nEntries */)
{
   RClusterDescriptor::RLocator result;
   result.fPosition = fClusterMinOffset;
   result.fBytesOnStorage = fClusterMaxOffset - fClusterMinOffset;
   fClusterMinOffset = std::uint64_t(-1);
   fClusterMaxOffset = 0;
   return result;
}


void ROOT::Experimental::Detail::RPageSinkFile::CommitDatasetImpl()
{
   const auto &descriptor = fDescriptorBuilder.GetDescriptor();
   auto szFooter = descriptor.GetFooterSize();
   auto buffer = std::make_unique<unsigned char []>(szFooter);
   descriptor.SerializeFooter(buffer.get());

   auto zipBuffer = std::make_unique<unsigned char []>(szFooter);
   auto szZipFooter = fCompressor(buffer.get(), szFooter, fOptions.GetCompression(),
      [&zipBuffer](const void *b, size_t n, size_t o){ memcpy(zipBuffer.get() + o, b, n); } );
   fWriter->WriteNTupleFooter(zipBuffer.get(), szZipFooter, szFooter);
   fWriter->Commit();
}


ROOT::Experimental::Detail::RPage
ROOT::Experimental::Detail::RPageSinkFile::ReservePage(ColumnHandle_t columnHandle, std::size_t nElements)
{
   if (nElements == 0)
      nElements = kDefaultElementsPerPage;
   auto elementSize = columnHandle.fColumn->GetElement()->GetSize();
   return fPageAllocator->NewPage(columnHandle.fId, elementSize, nElements);
}

void ROOT::Experimental::Detail::RPageSinkFile::ReleasePage(RPage &page)
{
   fPageAllocator->DeletePage(page);
}


////////////////////////////////////////////////////////////////////////////////


ROOT::Experimental::Detail::RPage ROOT::Experimental::Detail::RPageAllocatorFile::NewPage(
   ColumnId_t columnId, void *mem, std::size_t elementSize, std::size_t nElements)
{
   RPage newPage(columnId, mem, elementSize * nElements, elementSize);
   newPage.TryGrow(nElements);
   return newPage;
}

void ROOT::Experimental::Detail::RPageAllocatorFile::DeletePage(const RPage& page)
{
   if (page.IsNull())
      return;
   delete[] reinterpret_cast<unsigned char *>(page.GetBuffer());
}


////////////////////////////////////////////////////////////////////////////////


ROOT::Experimental::Detail::RPageSourceFile::RPageSourceFile(std::string_view ntupleName,
   const RNTupleReadOptions &options)
   : RPageSource(ntupleName, options)
   , fMetrics("RPageSourceFile")
   , fPageAllocator(std::make_unique<RPageAllocatorFile>())
   , fPagePool(std::make_shared<RPagePool>())
   , fClusterPool(std::make_unique<RClusterPool>(*this))
{
   fCounters = std::unique_ptr<RCounters>(new RCounters{
      *fMetrics.MakeCounter<RNTupleAtomicCounter*>("nReadV", "", "number of vector read requests"),
      *fMetrics.MakeCounter<RNTupleAtomicCounter*>("nRead", "", "number of byte ranges read"),
      *fMetrics.MakeCounter<RNTupleAtomicCounter*>("szReadPayload", "B", "volume read from file (required)"),
      *fMetrics.MakeCounter<RNTupleAtomicCounter*>("szReadOverhead", "B", "volume read from file (overhead)"),
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
            if (const auto szReadPayload = metrics.GetCounter("RPageSourceFile.szReadPayload")) {
               if (const auto szReadOverhead = metrics.GetCounter("RPageSourceFile.szReadOverhead")) {
                  if (const auto timeWallRead = metrics.GetCounter("RPageSourceFile.timeWallRead")) {
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
            if (const auto szUnzip = metrics.GetCounter("RPageSourceFile.szUnzip")) {
               if (const auto timeWallRead = metrics.GetCounter("RPageSourceFile.timeWallRead")) {
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
            if (const auto szUnzip = metrics.GetCounter("RPageSourceFile.szUnzip")) {
               if (const auto timeWallUnzip = metrics.GetCounter("RPageSourceFile.timeWallUnzip")) {
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
            if (const auto szReadPayload = metrics.GetCounter("RPageSourceFile.szReadPayload")) {
               if (const auto szReadOverhead = metrics.GetCounter("RPageSourceFile.szReadOverhead")) {
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
            if (const auto szReadPayload = metrics.GetCounter("RPageSourceFile.szReadPayload")) {
               if (const auto szUnzip = metrics.GetCounter("RPageSourceFile.szUnzip")) {
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


ROOT::Experimental::Detail::RPageSourceFile::RPageSourceFile(std::string_view ntupleName, std::string_view path,
   const RNTupleReadOptions &options)
   : RPageSourceFile(ntupleName, options)
{
   fFile = ROOT::Internal::RRawFile::Create(path);
   R__ASSERT(fFile);
   fReader = Internal::RMiniFileReader(fFile.get());
}


ROOT::Experimental::Detail::RPageSourceFile::~RPageSourceFile()
{
}


ROOT::Experimental::RNTupleDescriptor ROOT::Experimental::Detail::RPageSourceFile::AttachImpl()
{
   RNTupleDescriptorBuilder descBuilder;
   auto ntpl = fReader.GetNTuple(fNTupleName).Unwrap();

   auto buffer = std::make_unique<unsigned char[]>(ntpl.fLenHeader);
   auto zipBuffer = std::make_unique<unsigned char[]>(ntpl.fNBytesHeader);
   fReader.ReadBuffer(zipBuffer.get(), ntpl.fNBytesHeader, ntpl.fSeekHeader);
   fDecompressor(zipBuffer.get(), ntpl.fNBytesHeader, ntpl.fLenHeader, buffer.get());
   descBuilder.SetFromHeader(buffer.get());

   buffer = std::make_unique<unsigned char[]>(ntpl.fLenFooter);
   zipBuffer = std::make_unique<unsigned char[]>(ntpl.fNBytesFooter);
   fReader.ReadBuffer(zipBuffer.get(), ntpl.fNBytesFooter, ntpl.fSeekFooter);
   fDecompressor(zipBuffer.get(), ntpl.fNBytesFooter, ntpl.fLenFooter, buffer.get());
   descBuilder.AddClustersFromFooter(buffer.get());

   return descBuilder.MoveDescriptor();
}


ROOT::Experimental::Detail::RPage ROOT::Experimental::Detail::RPageSourceFile::PopulatePageFromCluster(
   ColumnHandle_t columnHandle, const RClusterDescriptor &clusterDescriptor, ClusterSize_t::ValueType idxInCluster)
{
   const auto columnId = columnHandle.fId;
   const auto clusterId = clusterDescriptor.GetId();
   const auto &pageRange = clusterDescriptor.GetPageRange(columnId);

   // TODO(jblomer): binary search
   RClusterDescriptor::RPageRange::RPageInfo pageInfo;
   decltype(idxInCluster) firstInPage = 0;
   NTupleSize_t pageNo = 0;
   for (const auto &pi : pageRange.fPageInfos) {
      if (firstInPage + pi.fNElements > idxInCluster) {
         pageInfo = pi;
         break;
      }
      firstInPage += pi.fNElements;
      ++pageNo;
   }
   R__ASSERT(firstInPage <= idxInCluster);
   R__ASSERT((firstInPage + pageInfo.fNElements) > idxInCluster);

   const auto element = columnHandle.fColumn->GetElement();
   const auto elementSize = element->GetSize();

   const auto bytesOnStorage = pageInfo.fLocator.fBytesOnStorage;
   const auto bytesPacked = (element->GetBitsOnStorage() * pageInfo.fNElements + 7) / 8;
   const auto pageSize = elementSize * pageInfo.fNElements;

   auto pageBuffer = new unsigned char[bytesPacked];
   if (fOptions.GetClusterCache() == RNTupleReadOptions::EClusterCache::kOff) {
      fReader.ReadBuffer(pageBuffer, bytesOnStorage, pageInfo.fLocator.fPosition);
      fCounters->fNPageLoaded.Inc();
   } else {
      //printf("GETTING CLUSTER %ld\n", clusterId);
      if (!fCurrentCluster || (fCurrentCluster->GetId() != clusterId) || !fCurrentCluster->ContainsColumn(columnId))
         fCurrentCluster = fClusterPool->GetCluster(clusterId, fActiveColumns);
      R__ASSERT(fCurrentCluster->ContainsColumn(columnId));

      auto cachedPage = fPagePool->GetPage(columnId, RClusterIndex(clusterId, idxInCluster));
      if (!cachedPage.IsNull())
         return cachedPage;

      ROnDiskPage::Key key(columnId, pageNo);
      //printf("POPULATE cluster %ld column %ld page %ld\n", clusterId, columnId, pageNo);
      auto onDiskPage = fCurrentCluster->GetOnDiskPage(key);
      R__ASSERT(onDiskPage);
      R__ASSERT(bytesOnStorage == onDiskPage->GetSize());
      memcpy(pageBuffer, onDiskPage->GetAddress(), onDiskPage->GetSize());
   }

   if (bytesOnStorage != bytesPacked) {
      RNTupleAtomicTimer timer(fCounters->fTimeWallUnzip, fCounters->fTimeCpuUnzip);
      fDecompressor(pageBuffer, bytesOnStorage, bytesPacked);
      fCounters->fSzUnzip.Add(bytesPacked);
   }

   if (!element->IsMappable()) {
      auto unpackedBuffer = new unsigned char[pageSize];
      element->Unpack(unpackedBuffer, pageBuffer, pageInfo.fNElements);
      delete[] pageBuffer;
      pageBuffer = unpackedBuffer;
   }

   const auto indexOffset = clusterDescriptor.GetColumnRange(columnId).fFirstElementIndex;
   auto newPage = fPageAllocator->NewPage(columnId, pageBuffer, elementSize, pageInfo.fNElements);
   newPage.SetWindow(indexOffset + firstInPage, RPage::RClusterInfo(clusterId, indexOffset));
   fPagePool->RegisterPage(newPage,
      RPageDeleter([](const RPage &page, void * /*userData*/)
      {
         RPageAllocatorFile::DeletePage(page);
      }, nullptr));
   fCounters->fNPagePopulated.Inc();
   return newPage;
}


ROOT::Experimental::Detail::RPage ROOT::Experimental::Detail::RPageSourceFile::PopulatePage(
   ColumnHandle_t columnHandle, NTupleSize_t globalIndex)
{
   const auto columnId = columnHandle.fId;
   auto cachedPage = fPagePool->GetPage(columnId, globalIndex);
   if (!cachedPage.IsNull())
      return cachedPage;

   const auto clusterId = fDescriptor.FindClusterId(columnId, globalIndex);
   R__ASSERT(clusterId != kInvalidDescriptorId);
   const auto &clusterDescriptor = fDescriptor.GetClusterDescriptor(clusterId);
   const auto selfOffset = clusterDescriptor.GetColumnRange(columnId).fFirstElementIndex;
   R__ASSERT(selfOffset <= globalIndex);
   return PopulatePageFromCluster(columnHandle, clusterDescriptor, globalIndex - selfOffset);
}


ROOT::Experimental::Detail::RPage ROOT::Experimental::Detail::RPageSourceFile::PopulatePage(
   ColumnHandle_t columnHandle, const RClusterIndex &clusterIndex)
{
   const auto clusterId = clusterIndex.GetClusterId();
   const auto idxInCluster = clusterIndex.GetIndex();
   const auto columnId = columnHandle.fId;
   auto cachedPage = fPagePool->GetPage(columnId, clusterIndex);
   if (!cachedPage.IsNull())
      return cachedPage;

   R__ASSERT(clusterId != kInvalidDescriptorId);
   const auto &clusterDescriptor = fDescriptor.GetClusterDescriptor(clusterId);
   return PopulatePageFromCluster(columnHandle, clusterDescriptor, idxInCluster);
}

void ROOT::Experimental::Detail::RPageSourceFile::ReleasePage(RPage &page)
{
   fPagePool->ReturnPage(page);
}

std::unique_ptr<ROOT::Experimental::Detail::RPageSource> ROOT::Experimental::Detail::RPageSourceFile::Clone() const
{
   auto clone = new RPageSourceFile(fNTupleName, fOptions);
   clone->fFile = fFile->Clone();
   clone->fReader = Internal::RMiniFileReader(clone->fFile.get());
   return std::unique_ptr<RPageSourceFile>(clone);
}

std::unique_ptr<ROOT::Experimental::Detail::RCluster>
ROOT::Experimental::Detail::RPageSourceFile::LoadCluster(DescriptorId_t clusterId, const ColumnSet_t &columns)
{
   fCounters->fNClusterLoaded.Inc();

   const auto &clusterDesc = GetDescriptor().GetClusterDescriptor(clusterId);
   auto clusterLocator = clusterDesc.GetLocator();
   auto clusterSize = clusterLocator.fBytesOnStorage;
   R__ASSERT(clusterSize > 0);

   struct ROnDiskPageLocator {
      ROnDiskPageLocator() = default;
      ROnDiskPageLocator(DescriptorId_t c, NTupleSize_t p, std::uint64_t o, std::uint64_t s)
         : fColumnId(c), fPageNo(p), fOffset(o), fSize(s) {}
      DescriptorId_t fColumnId = 0;
      NTupleSize_t fPageNo = 0;
      std::uint64_t fOffset = 0;
      std::uint64_t fSize = 0;
      std::size_t fBufPos = 0;
   };

   // Collect the page necessary page meta-data and sum up the total size of the compressed and packed pages
   std::vector<ROnDiskPageLocator> onDiskPages;
   auto activeSize = 0;
   for (auto columnId : columns) {
      const auto &pageRange = clusterDesc.GetPageRange(columnId);
      NTupleSize_t pageNo = 0;
      for (const auto &pageInfo : pageRange.fPageInfos) {
         const auto &pageLocator = pageInfo.fLocator;
         activeSize += pageLocator.fBytesOnStorage;
         onDiskPages.emplace_back(ROnDiskPageLocator(
            columnId, pageNo, pageLocator.fPosition, pageLocator.fBytesOnStorage));
         ++pageNo;
      }
   }

   // Linearize the page requests by file offset
   std::sort(onDiskPages.begin(), onDiskPages.end(),
      [](const ROnDiskPageLocator &a, const ROnDiskPageLocator &b) {return a.fOffset < b.fOffset;});

   // In order to coalesce close-by pages, we collect the sizes of the gaps between pages on disk.  We then order
   // the gaps by size, sum them up and find a cutoff for the largest gap that we tolerate when coalescing pages.
   // The size of the cutoff is given by the fraction of extra bytes we are willing to read in order to reduce
   // the number of read requests.  We thus schedule the lowest number of requests given a tolerable fraction
   // of extra bytes.
   // TODO(jblomer): Eventually we may want to select the parameter at runtime according to link latency and speed,
   // memory consumption, device block size.
   float maxOverhead = 0.25 * float(activeSize);
   std::vector<std::size_t> gaps;
   for (unsigned i = 1; i < onDiskPages.size(); ++i) {
      gaps.emplace_back(onDiskPages[i].fOffset - (onDiskPages[i-1].fSize + onDiskPages[i-1].fOffset));
   }
   std::sort(gaps.begin(), gaps.end());
   std::size_t gapCut = 0;
   float szExtra = 0.0;
   for (auto g : gaps) {
      szExtra += g;
      if (szExtra  > maxOverhead)
         break;
      gapCut = g;
   }

   // Prepare the input vector for the RRawFile::ReadV() call
   struct RReadRequest {
      RReadRequest() = default;
      RReadRequest(std::size_t b, std::uint64_t o, std::uint64_t s) : fBufPos(b), fOffset(o), fSize(s) {}
      std::size_t fBufPos = 0;
      std::uint64_t fOffset = 0;
      std::uint64_t fSize = 0;
   };
   std::vector<ROOT::Internal::RRawFile::RIOVec> readRequests;

   ROOT::Internal::RRawFile::RIOVec req;
   std::size_t szPayload = 0;
   std::size_t szOverhead = 0;
   for (auto &s : onDiskPages) {
      R__ASSERT(s.fSize > 0);
      auto readUpTo = req.fOffset + req.fSize;
      R__ASSERT(s.fOffset >= readUpTo);
      auto overhead = s.fOffset - readUpTo;
      szPayload += s.fSize;
      if (overhead <= gapCut) {
         szOverhead += overhead;
         s.fBufPos = reinterpret_cast<intptr_t>(req.fBuffer) + req.fSize + overhead;
         req.fSize += overhead + s.fSize;
         continue;
      }

      // close the current request and open new one
      if (req.fSize > 0)
         readRequests.emplace_back(req);

      req.fBuffer = reinterpret_cast<unsigned char *>(req.fBuffer) + req.fSize;
      s.fBufPos = reinterpret_cast<intptr_t>(req.fBuffer);

      req.fOffset = s.fOffset;
      req.fSize = s.fSize;
   }
   readRequests.emplace_back(req);
   fCounters->fSzReadPayload.Add(szPayload);
   fCounters->fSzReadOverhead.Add(szOverhead);

   // Register the on disk pages in a page map
   auto buffer = new unsigned char[reinterpret_cast<intptr_t>(req.fBuffer) + req.fSize];
   auto pageMap = std::make_unique<ROnDiskPageMapHeap>(std::unique_ptr<unsigned char []>(buffer));
   for (const auto &s : onDiskPages) {
      ROnDiskPage::Key key(s.fColumnId, s.fPageNo);
      pageMap->Register(key, ROnDiskPage(buffer + s.fBufPos, s.fSize));
   }
   fCounters->fNPageLoaded.Add(onDiskPages.size());
   for (auto &r : readRequests) {
      r.fBuffer = buffer + reinterpret_cast<intptr_t>(r.fBuffer);
   }

   auto nReqs = readRequests.size();
   {
      RNTupleAtomicTimer timer(fCounters->fTimeWallRead, fCounters->fTimeCpuRead);
      fFile->ReadV(&readRequests[0], nReqs);
   }
   fCounters->fNReadV.Inc();
   fCounters->fNRead.Add(nReqs);

   auto cluster = std::make_unique<RCluster>(clusterId);
   cluster->Adopt(std::move(pageMap));
   for (auto colId : columns)
      cluster->SetColumnAvailable(colId);
   return cluster;
}


void ROOT::Experimental::Detail::RPageSourceFile::UnzipClusterImpl(RCluster *cluster)
{
   RNTupleAtomicTimer timer(fCounters->fTimeWallUnzip, fCounters->fTimeCpuUnzip);
   fTaskScheduler->Reset();

   const auto clusterId = cluster->GetId();
   const auto &clusterDescriptor = fDescriptor.GetClusterDescriptor(clusterId);

   std::vector<std::unique_ptr<RColumnElementBase>> allElements;

   const auto &columnsInCluster = cluster->GetAvailColumns();
   for (const auto columnId : columnsInCluster) {
      const auto &columnDesc = fDescriptor.GetColumnDescriptor(columnId);

      allElements.emplace_back(RColumnElementBase::Generate(columnDesc.GetModel().GetType()));

      const auto &pageRange = clusterDescriptor.GetPageRange(columnId);
      std::uint64_t pageNo = 0;
      std::uint64_t firstInPage = 0;
      for (const auto &pi : pageRange.fPageInfos) {
         ROnDiskPage::Key key(columnId, pageNo);
         auto onDiskPage = cluster->GetOnDiskPage(key);
         R__ASSERT(onDiskPage);
         R__ASSERT(onDiskPage->GetSize() == pi.fLocator.fBytesOnStorage);

         auto taskFunc =
            [this, columnId, clusterId, firstInPage, onDiskPage,
             element = allElements.back().get(),
             nElements = pi.fNElements,
             indexOffset = clusterDescriptor.GetColumnRange(columnId).fFirstElementIndex
            ] () {
               const auto bytesPacked = (element->GetBitsOnStorage() * nElements + 7) / 8;
               const auto pageSize = element->GetSize() * nElements;

               auto pageBufferPacked = new unsigned char[bytesPacked];
               if (onDiskPage->GetSize() != bytesPacked) {
                  fDecompressor(onDiskPage->GetAddress(), onDiskPage->GetSize(), bytesPacked, pageBufferPacked);
                  fCounters->fSzUnzip.Add(bytesPacked);
               } else {
                  // We cannot simply map the onDiskPage because the cluster pool and the page pool have different
                  // life times
                  memcpy(pageBufferPacked, onDiskPage->GetAddress(), bytesPacked);
               }

               auto pageBuffer = pageBufferPacked;
               if (!element->IsMappable()) {
                  pageBuffer = new unsigned char[pageSize];
                  element->Unpack(pageBuffer, pageBufferPacked, nElements);
                  delete[] pageBufferPacked;
               }

               auto newPage = fPageAllocator->NewPage(columnId, pageBuffer, element->GetSize(), nElements);
               newPage.SetWindow(indexOffset + firstInPage, RPage::RClusterInfo(clusterId, indexOffset));
               fPagePool->PreloadPage(newPage,
                  RPageDeleter([](const RPage &page, void * /*userData*/)
                  {
                     RPageAllocatorFile::DeletePage(page);
                  }, nullptr));
            };

         fTaskScheduler->AddTask(taskFunc);

         firstInPage += pi.fNElements;
         pageNo++;
      } // for all pages in column
   } // for all columns in cluster

   fCounters->fNPagePopulated.Add(cluster->GetNOnDiskPages());

   fTaskScheduler->Wait();
}
