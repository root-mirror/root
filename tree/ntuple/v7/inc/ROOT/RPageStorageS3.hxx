/// \file ROOT/RPageStorageS3.hxx
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

#ifndef ROOT7_RPageStorageS3
#define ROOT7_RPageStorageS3

#include <ROOT/RPageStorage.hxx>
#include <ROOT/RNTupleMetrics.hxx>
#include <ROOT/RNTupleZip.hxx>
#include <ROOT/RStringView.hxx>

#include <array>
#include <cstdio>
#include <memory>
#include <string>

namespace ROOT {
namespace Experimental {
namespace Detail {

class RCluster;
class RPageAllocatorHeap;
class RPagePool;
class RS3Handle;

// clang-format off
/**
\class ROOT::Experimental::Detail::RS3NTupleAnchor
\ingroup NTuple
\brief Entry point for an RNTuple in an S3 bucket. It encodes essential
information to read the ntuple; currently, it contains (un)compressed size of
the header/footer blobs.
*/
// clang-format on
struct RS3NTupleAnchor {
   /// Allows for evolving the struct in future versions
   std::uint32_t fVersion = 0;
   /// The size of the compressed ntuple header
   std::uint32_t fNBytesHeader = 0;
   /// The size of the uncompressed ntuple header
   std::uint32_t fLenHeader = 0;
   /// The size of the compressed ntuple footer
   std::uint32_t fNBytesFooter = 0;
   /// The size of the uncompressed ntuple footer
   std::uint32_t fLenFooter = 0;

   bool operator ==(const RS3NTupleAnchor &other) const {
      return fVersion == other.fVersion &&
         fNBytesHeader == other.fNBytesHeader &&
         fLenHeader == other.fLenHeader &&
         fNBytesFooter == other.fNBytesFooter &&
         fLenFooter == other.fLenFooter;
   }

   std::uint32_t Serialize(void *buffer) const;
   std::uint32_t Deserialize(const void *buffer);

   static std::uint32_t GetSize()
   { return RS3NTupleAnchor().Serialize(nullptr); }
};

// clang-format off
/**
\class ROOT::Experimental::Detail::RPageSinkS3
\ingroup NTuple
\brief Storage provider that writes %RNTuple pages to an S3 bucket at the specified path

Currently, an object is allocated for each page + 3 additional objects (anchor/header/footer).
*/
// clang-format on
class RPageSinkS3 : public RPageSink {
private:
   std::unique_ptr<RPageAllocatorHeap> fPageAllocator;
   /// Handle to the S3 location.
   std::unique_ptr<RS3Handle> fS3Handle;
   /// A URI to an S3 location
   std::string fUri;
   /// Object ID for the next committed page; it is automatically incremented in CommitSealedPageImpl()`
   std::atomic<std::uint64_t> fOid{0};

   RS3NTupleAnchor fNTupleAnchor;

protected:
   void CreateImpl(const RNTupleModel &model) final;
   RClusterDescriptor::RLocator CommitPageImpl(ColumnHandle_t columnHandle, const RPage &page) final;
   RClusterDescriptor::RLocator CommitSealedPageImpl(DescriptorId_t columnId,
                                                     const RPageStorage::RSealedPage &sealedPage) final;
   RClusterDescriptor::RLocator CommitClusterImpl(NTupleSize_t nEntries) final;
   void CommitDatasetImpl() final;
   void WriteNTupleHeader(const void *data, size_t nbytes, size_t lenHeader);
   void WriteNTupleFooter(const void *data, size_t nbytes, size_t lenFooter);
   void WriteNTupleAnchor();

public:
   RPageSinkS3(std::string_view ntupleName, std::string_view uri, const RNTupleWriteOptions &options);
   virtual ~RPageSinkS3();

   RPage ReservePage(ColumnHandle_t columnHandle, std::size_t nElements = 0) final;
   void ReleasePage(RPage &page) final;
};

// clang-format off
/**
\class ROOT::Experimental::Detail::RPageAllocatorS3
\ingroup NTuple
\brief Manages pages read from S3
*/
// clang-format on
class RPageAllocatorS3 {
public:
   static RPage NewPage(ColumnId_t columnId, void *mem, std::size_t elementSize, std::size_t nElements);
   static void DeletePage(const RPage& page);
};

// clang-format off
/**
\class ROOT::Experimental::Detail::RPageSourceS3
\ingroup NTuple
\brief Storage provider that reads %RNTuple pages from an S3 bucket

The S3 configuration is specified using the following environment variables:
- `S3_ACCESS_KEY`
- `S3_SECRET_KEY`
- `S3_REGION`

*/
// clang-format on
class RPageSourceS3 : public RPageSource {
private:
   std::unique_ptr<RPageAllocatorS3> fPageAllocator;
   std::shared_ptr<RPagePool> fPagePool;
   std::unique_ptr<RS3Handle> fS3Handle;
   std::string fUri;

   RPage PopulatePageFromCluster(ColumnHandle_t columnHandle,
      const RClusterDescriptor &clusterDescriptor, ClusterSize_t::ValueType idxInCluster);

protected:
   RNTupleDescriptor AttachImpl() final;

public:
   RPageSourceS3(std::string_view ntupleName, std::string_view uri, const RNTupleReadOptions &options);
   std::unique_ptr<RPageSource> Clone() const final;
   virtual ~RPageSourceS3();

   RPage PopulatePage(ColumnHandle_t columnHandle, NTupleSize_t globalIndex) final;
   RPage PopulatePage(ColumnHandle_t columnHandle, const RClusterIndex &clusterIndex) final;
   void ReleasePage(RPage &page) final;

   void LoadSealedPage(DescriptorId_t columnId, const RClusterIndex &clusterIndex,
                       RSealedPage &sealedPage) final;

   std::unique_ptr<RCluster> LoadCluster(DescriptorId_t clusterId, const ColumnSet_t &columns) final;
};

} // namespace Detail

} // namespace Experimental
} // namespace ROOT

#endif
