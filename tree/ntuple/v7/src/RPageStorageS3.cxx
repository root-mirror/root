/// \file RPageStorageS3.cxx
/// \ingroup NTuple ROOT7
/// \author Max Orok <maxwellorok@gmail.com>
/// \date 2021-06-03
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2021, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/RCluster.hxx>
#include <ROOT/RClusterPool.hxx>
#include <ROOT/RColumn.hxx>
#include <ROOT/RPageAllocator.hxx>
#include <ROOT/RPagePool.hxx>
#include <ROOT/RPageStorageS3.hxx>

#include <iostream>
#include <string>

#include "davix.hpp"

namespace {
   static const std::string kOidAnchor = std::to_string(std::uint64_t(-1));
   static const std::string kOidHeader = std::to_string(std::uint64_t(-2));
   static const std::string kOidFooter = std::to_string(std::uint64_t(-3));
} // anonymous namespace

namespace ROOT {
namespace Experimental {
namespace Detail {

class RS3Handle {
private:
   std::string fUri;
   Davix::Context fCtx;
   Davix::RequestParams fReqParams;
public:
   RS3Handle() = default;
   explicit RS3Handle(std::string_view ntupleName, std::string_view uri) {
      R__ASSERT(!ntupleName.empty() && !uri.empty());
      fUri += uri;
      if (uri.back() != '/') {
         fUri += '/';
      }
      fUri += ntupleName;
      if (ntupleName.back() != '/') {
         fUri += '/';
      }
      const char *s3_sec = getenv("S3_SECRET_KEY");
      const char *s3_acc = getenv("S3_ACCESS_KEY");
      if (s3_sec && s3_acc) {
         fReqParams.setAwsAuthorizationKeys(s3_sec, s3_acc);
      }
      const char *s3_reg = getenv("S3_REGION");
      if (s3_reg) {
         fReqParams.setAwsRegion(s3_reg);
      }
   }
   RS3Handle(const RS3Handle&) = delete;
   RS3Handle& operator=(const RS3Handle&) = delete;
   RS3Handle(RS3Handle&&) = default;
   RS3Handle& operator=(RS3Handle&&) = default;
   ~RS3Handle() = default;

   void ReadObject(std::string_view key, std::vector<char> &buf) {
      Davix::Uri uri(fUri + std::string(key));
      Davix::DavFile obj(fCtx, fReqParams, uri);
      // Danger: there are no size limits on the amount of data read into the buffer
      try {
         obj.get(nullptr, buf);
      } catch (const Davix::DavixException &err) {
         throw ROOT::Experimental::RException(R__FAIL(std::string("S3 read error: ") + err.what()));
      }
   }
   void WriteObject(std::string_view key, const void *buffer, std::size_t nbytes) {
      Davix::Uri uri(fUri + std::string(key));
      Davix::DavFile target(fCtx, fReqParams, uri);
      try {
         target.put(nullptr, static_cast<const char*>(buffer), nbytes);
      } catch (const Davix::DavixException &err) {
         throw ROOT::Experimental::RException(R__FAIL(std::string("S3 write error: ") + err.what()));
      }
   }
};

} // namespace Detail
} // namespace Experimental
} // namespace ROOT

std::uint32_t
ROOT::Experimental::Detail::RS3NTupleAnchor::Serialize(void *buffer) const
{
   using namespace ROOT::Experimental::Internal::RNTupleSerialization;
   if (buffer != nullptr) {
      auto bytes = reinterpret_cast<unsigned char *>(buffer);
      bytes += SerializeUInt32(fVersion, bytes);
      bytes += SerializeUInt32(fNBytesHeader, bytes);
      bytes += SerializeUInt32(fLenHeader, bytes);
      bytes += SerializeUInt32(fNBytesFooter, bytes);
      bytes += SerializeUInt32(fLenFooter, bytes);
   }
   return 20;
}

std::uint32_t
ROOT::Experimental::Detail::RS3NTupleAnchor::Deserialize(const void *buffer)
{
   using namespace ROOT::Experimental::Internal::RNTupleSerialization;
   auto bytes = reinterpret_cast<const unsigned char *>(buffer);
   bytes += DeserializeUInt32(bytes, &fVersion);
   bytes += DeserializeUInt32(bytes, &fNBytesHeader);
   bytes += DeserializeUInt32(bytes, &fLenHeader);
   bytes += DeserializeUInt32(bytes, &fNBytesFooter);
   bytes += DeserializeUInt32(bytes, &fLenFooter);
   return 20;
}


////////////////////////////////////////////////////////////////////////////////


ROOT::Experimental::Detail::RPageSinkS3::RPageSinkS3(std::string_view ntupleName, std::string_view uri,
   const RNTupleWriteOptions &options)
   : RPageSink(ntupleName, options)
   , fPageAllocator(std::make_unique<RPageAllocatorHeap>())
   , fUri(uri)
{
   R__LOG_WARNING(NTupleLog()) << "The S3 backend is experimental and still under development. " <<
      "Do not store real data with this version of RNTuple!";
   fCompressor = std::make_unique<RNTupleCompressor>();
   EnableDefaultMetrics("RPageSinkS3");
   fS3Handle = std::make_unique<RS3Handle>(ntupleName, fUri);
}

ROOT::Experimental::Detail::RPageSinkS3::~RPageSinkS3() = default;

void ROOT::Experimental::Detail::RPageSinkS3::CreateImpl(const RNTupleModel & /* model */)
{
   const auto &descriptor = fDescriptorBuilder.GetDescriptor();
   auto szHeader = descriptor.GetHeaderSize();
   auto buffer = std::make_unique<unsigned char[]>(szHeader);
   descriptor.SerializeHeader(buffer.get());

   auto zipBuffer = std::make_unique<unsigned char[]>(szHeader);
   auto szZipHeader = fCompressor->Zip(buffer.get(), szHeader, fOptions->GetCompression(),
      [&zipBuffer](const void *b, size_t n, size_t o){ memcpy(zipBuffer.get() + o, b, n); } );
   WriteNTupleHeader(zipBuffer.get(), szZipHeader, szHeader);
}


ROOT::Experimental::RClusterDescriptor::RLocator
ROOT::Experimental::Detail::RPageSinkS3::CommitPageImpl(ColumnHandle_t columnHandle, const RPage &page)
{
   auto element = columnHandle.fColumn->GetElement();
   RPageStorage::RSealedPage sealedPage;
   {
      RNTupleAtomicTimer timer(fCounters->fTimeWallZip, fCounters->fTimeCpuZip);
      sealedPage = SealPage(page, *element, fOptions->GetCompression());
   }

   fCounters->fSzZip.Add(page.GetSize());
   return CommitSealedPageImpl(columnHandle.fId, sealedPage);
}

ROOT::Experimental::RClusterDescriptor::RLocator
ROOT::Experimental::Detail::RPageSinkS3::CommitSealedPageImpl(
   DescriptorId_t /*columnId*/, const RPageStorage::RSealedPage &sealedPage)
{
   auto offsetData = fOid.fetch_add(1);
   {
      RNTupleAtomicTimer timer(fCounters->fTimeWallWrite, fCounters->fTimeCpuWrite);
      fS3Handle->WriteObject(std::to_string(offsetData), sealedPage.fBuffer, sealedPage.fSize);
   }

   RClusterDescriptor::RLocator result;
   result.fPosition = offsetData;
   result.fBytesOnStorage = sealedPage.fSize;
   fCounters->fNPageCommitted.Inc();
   fCounters->fSzWritePayload.Add(sealedPage.fSize);
   return result;
}

// todo(max) figure out we should implement CommitClusterImpl
ROOT::Experimental::RClusterDescriptor::RLocator
ROOT::Experimental::Detail::RPageSinkS3::CommitClusterImpl(ROOT::Experimental::NTupleSize_t /* nEntries */)
{
   return {};
}

void ROOT::Experimental::Detail::RPageSinkS3::CommitDatasetImpl()
{
   const auto &descriptor = fDescriptorBuilder.GetDescriptor();
   auto szFooter = descriptor.GetFooterSize();
   auto buffer = std::make_unique<unsigned char []>(szFooter);
   descriptor.SerializeFooter(buffer.get());

   auto zipBuffer = std::make_unique<unsigned char []>(szFooter);
   auto szZipFooter = fCompressor->Zip(buffer.get(), szFooter, fOptions->GetCompression(),
      [&zipBuffer](const void *b, size_t n, size_t o){ memcpy(zipBuffer.get() + o, b, n); } );
   WriteNTupleFooter(zipBuffer.get(), szZipFooter, szFooter);
   WriteNTupleAnchor();
}

void ROOT::Experimental::Detail::RPageSinkS3::WriteNTupleHeader(
		const void *data, size_t nbytes, size_t lenHeader)
{
   fS3Handle->WriteObject(kOidHeader, data, nbytes);
   fNTupleAnchor.fLenHeader = lenHeader;
   fNTupleAnchor.fNBytesHeader = nbytes;
}

void ROOT::Experimental::Detail::RPageSinkS3::WriteNTupleFooter(
		const void *data, size_t nbytes, size_t lenFooter)
{
   fS3Handle->WriteObject(kOidFooter, data, nbytes);
   fNTupleAnchor.fLenFooter = lenFooter;
   fNTupleAnchor.fNBytesFooter = nbytes;
}

void ROOT::Experimental::Detail::RPageSinkS3::WriteNTupleAnchor() {
   const auto ntplSize = RS3NTupleAnchor::GetSize();
   auto buffer = std::make_unique<unsigned char[]>(ntplSize);
   fNTupleAnchor.Serialize(buffer.get());
   fS3Handle->WriteObject(kOidAnchor, buffer.get(), ntplSize);
}

ROOT::Experimental::Detail::RPage
ROOT::Experimental::Detail::RPageSinkS3::ReservePage(ColumnHandle_t columnHandle, std::size_t nElements)
{
   if (nElements == 0)
      nElements = fOptions->GetNElementsPerPage();
   auto elementSize = columnHandle.fColumn->GetElement()->GetSize();
   return fPageAllocator->NewPage(columnHandle.fId, elementSize, nElements);
}

void ROOT::Experimental::Detail::RPageSinkS3::ReleasePage(RPage &page)
{
   fPageAllocator->DeletePage(page);
}


////////////////////////////////////////////////////////////////////////////////


ROOT::Experimental::Detail::RPage ROOT::Experimental::Detail::RPageAllocatorS3::NewPage(
   ColumnId_t columnId, void *mem, std::size_t elementSize, std::size_t nElements)
{
   RPage newPage(columnId, mem, elementSize * nElements, elementSize);
   newPage.TryGrow(nElements);
   return newPage;
}

void ROOT::Experimental::Detail::RPageAllocatorS3::DeletePage(const RPage& page)
{
   if (page.IsNull())
      return;
   delete[] reinterpret_cast<unsigned char *>(page.GetBuffer());
}


////////////////////////////////////////////////////////////////////////////////


ROOT::Experimental::Detail::RPageSourceS3::RPageSourceS3(
   std::string_view ntupleName, std::string_view uri, const RNTupleReadOptions &options)
   : RPageSource(ntupleName, options)
   , fPageAllocator(std::make_unique<RPageAllocatorS3>())
   , fPagePool(std::make_shared<RPagePool>())
   , fUri(uri)
{
   RPageSource::fDecompressor = std::make_unique<RNTupleDecompressor>();
   EnableDefaultMetrics("RPageSourceS3");
   fS3Handle = std::make_unique<RS3Handle>(ntupleName, fUri);
}

ROOT::Experimental::Detail::RPageSourceS3::~RPageSourceS3() = default;

ROOT::Experimental::RNTupleDescriptor ROOT::Experimental::Detail::RPageSourceS3::AttachImpl()
{
   R__ASSERT(fS3Handle);
   RNTupleDescriptorBuilder descBuilder;

   RS3NTupleAnchor ntpl;
   std::vector<char> buf;
   buf.reserve(RS3NTupleAnchor::GetSize());
   fS3Handle->ReadObject(kOidAnchor, buf);
   if (buf.size() != RS3NTupleAnchor::GetSize()) {
      throw ROOT::Experimental::RException(R__FAIL("error reading RNTuple anchor"));
   }
   ntpl.Deserialize(static_cast<void*>(buf.data()));

   buf.clear();
   buf.reserve(ntpl.fNBytesHeader);
   fS3Handle->ReadObject(kOidHeader, buf);
   if (buf.size() != ntpl.fNBytesHeader) {
      throw ROOT::Experimental::RException(R__FAIL("error reading RNTuple header"));
   }
   auto unzipBuf = std::make_unique<unsigned char[]>(ntpl.fLenHeader);
   fDecompressor->Unzip(buf.data(), ntpl.fNBytesHeader, ntpl.fLenHeader, unzipBuf.get());
   descBuilder.SetFromHeader(unzipBuf.get());

   buf.clear();
   buf.reserve(ntpl.fNBytesFooter);
   fS3Handle->ReadObject(kOidFooter, buf);
   if (buf.size() != ntpl.fNBytesFooter) {
      throw ROOT::Experimental::RException(R__FAIL("error reading RNTuple footer"));
   }
   unzipBuf = std::make_unique<unsigned char[]>(ntpl.fLenFooter);
   fDecompressor->Unzip(buf.data(), ntpl.fNBytesFooter, ntpl.fLenFooter, unzipBuf.get());
   descBuilder.AddClustersFromFooter(unzipBuf.get());

   return descBuilder.MoveDescriptor();
}

// todo(max) implement LoadSealedPage
void ROOT::Experimental::Detail::RPageSourceS3::LoadSealedPage(
   DescriptorId_t columnId, const RClusterIndex &clusterIndex, RSealedPage &sealedPage)
{
   (void)columnId;
   (void)clusterIndex;
   (void)sealedPage;
}

ROOT::Experimental::Detail::RPage ROOT::Experimental::Detail::RPageSourceS3::PopulatePageFromCluster(
   ColumnHandle_t columnHandle, const RClusterDescriptor &clusterDescriptor, ClusterSize_t::ValueType idxInCluster)
{
   const auto columnId = columnHandle.fId;
   const auto clusterId = clusterDescriptor.GetId();

   auto pageInfo = clusterDescriptor.GetPageRange(columnId).Find(idxInCluster);

   const auto element = columnHandle.fColumn->GetElement();
   const auto elementSize = element->GetSize();
   const auto bytesOnStorage = pageInfo.fLocator.fBytesOnStorage;

   const void *sealedPageBuffer = nullptr; // points either to directReadBuffer or to a read-only page in the cluster
   std::unique_ptr<unsigned char []> directReadBuffer; // only used if cluster pool is turned off

   // todo(max) enable cluster cache
   std::vector<char> buf;
   buf.reserve(bytesOnStorage);
   fS3Handle->ReadObject(std::to_string(pageInfo.fLocator.fPosition), buf);
   R__ASSERT(buf.size() == bytesOnStorage);
   directReadBuffer = std::make_unique<unsigned char[]>(bytesOnStorage);
   memcpy(directReadBuffer.get(), buf.data(), bytesOnStorage);
   fCounters->fNPageLoaded.Inc();
   fCounters->fNRead.Inc();
   fCounters->fSzReadPayload.Add(bytesOnStorage);
   sealedPageBuffer = directReadBuffer.get();

   std::unique_ptr<unsigned char []> pageBuffer;
   {
      RNTupleAtomicTimer timer(fCounters->fTimeWallUnzip, fCounters->fTimeCpuUnzip);
      pageBuffer = UnsealPage({sealedPageBuffer, bytesOnStorage, pageInfo.fNElements}, *element);
      fCounters->fSzUnzip.Add(elementSize * pageInfo.fNElements);
   }

   const auto indexOffset = clusterDescriptor.GetColumnRange(columnId).fFirstElementIndex;
   auto newPage = fPageAllocator->NewPage(columnId, pageBuffer.release(), elementSize, pageInfo.fNElements);
   newPage.SetWindow(indexOffset + pageInfo.fFirstInPage, RPage::RClusterInfo(clusterId, indexOffset));
   fPagePool->RegisterPage(newPage,
      RPageDeleter([](const RPage &page, void * /*userData*/)
      {
         RPageAllocatorS3::DeletePage(page);
      }, nullptr));
   fCounters->fNPagePopulated.Inc();
   return newPage;
}

ROOT::Experimental::Detail::RPage ROOT::Experimental::Detail::RPageSourceS3::PopulatePage(
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

ROOT::Experimental::Detail::RPage ROOT::Experimental::Detail::RPageSourceS3::PopulatePage(
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

void ROOT::Experimental::Detail::RPageSourceS3::ReleasePage(RPage &page)
{
   fPagePool->ReturnPage(page);
}

std::unique_ptr<ROOT::Experimental::Detail::RPageSource>
ROOT::Experimental::Detail::RPageSourceS3::Clone() const {
   return std::make_unique<RPageSourceS3>(fNTupleName, fUri, fOptions);
}

std::unique_ptr<ROOT::Experimental::Detail::RCluster>
ROOT::Experimental::Detail::RPageSourceS3::LoadCluster(DescriptorId_t clusterId, const ColumnSet_t &columns)
{
   (void)columns;
   return std::make_unique<RCluster>(clusterId);
}
