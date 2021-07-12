/// \file ROOT/RPageStorageFile.hxx
/// \ingroup NTuple ROOT7
/// \author Jakob Blomer <jblomer@cern.ch>
/// \date 2019-11-21
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RPageStorageFile
#define ROOT7_RPageStorageFile

#include <ROOT/RMiniFile.hxx>
#include <ROOT/RNTupleZip.hxx>
#include <ROOT/RPageStorage.hxx>
#include <ROOT/RStringView.hxx>

#include <array>
#include <cstdio>
#include <memory>
#include <string>
#include <utility>

class TFile;

namespace ROOT {

namespace Internal {
class RRawFile;
}

namespace Experimental {
namespace Detail {

class RClusterPool;
class RPageAllocatorHeap;
class RPagePool;


// clang-format off
/**
\class ROOT::Experimental::Detail::RPageSinkFile
\ingroup NTuple
\brief Storage provider that write ntuple pages into a file

The written file can be either in ROOT format or in RNTuple bare format.
*/
// clang-format on
class RPageSinkFile : public RPageSink {
private:
   std::unique_ptr<RPageAllocatorHeap> fPageAllocator;

   std::unique_ptr<Internal::RNTupleFileWriter> fWriter;
   /// Byte offset of the first page of the current cluster
   std::uint64_t fClusterMinOffset = std::uint64_t(-1);
   /// Byte offset of the end of the last page of the current cluster
   std::uint64_t fClusterMaxOffset = 0;
   RPageSinkFile(std::string_view ntupleName, const RNTupleWriteOptions &options);

   RClusterDescriptor::RLocator WriteSealedPage(const RPageStorage::RSealedPage &sealedPage,
                                                std::size_t bytesPacked);

protected:
   void CreateImpl(const RNTupleModel &model) final;
   RClusterDescriptor::RLocator CommitPageImpl(ColumnHandle_t columnHandle, const RPage &page) final;
   RClusterDescriptor::RLocator CommitSealedPageImpl(DescriptorId_t columnId,
                                                     const RPageStorage::RSealedPage &sealedPage) final;
   RClusterDescriptor::RLocator CommitClusterImpl(NTupleSize_t nEntries) final;
   void CommitDatasetImpl() final;

public:
   RPageSinkFile(std::string_view ntupleName, std::string_view path, const RNTupleWriteOptions &options);
   RPageSinkFile(std::string_view ntupleName, std::string_view path, const RNTupleWriteOptions &options,
                 std::unique_ptr<TFile> &file);
   RPageSinkFile(std::string_view ntupleName, TFile &file, const RNTupleWriteOptions &options);
   RPageSinkFile(const RPageSinkFile&) = delete;
   RPageSinkFile& operator=(const RPageSinkFile&) = delete;
   RPageSinkFile(RPageSinkFile&&) = default;
   RPageSinkFile& operator=(RPageSinkFile&&) = default;
   virtual ~RPageSinkFile();

   RPage ReservePage(ColumnHandle_t columnHandle, std::size_t nElements = 0) final;
   void ReleasePage(RPage &page) final;
};

class RPageCache;

// clang-format off
/**
\class ROOT::Experimental::Detail::RPageAllocatorFile
\ingroup NTuple
\brief Manages pages read from the file
*/
// clang-format on
class RPageAllocatorFile {
private:
   std::unique_ptr<RPageCache> fPageCache;
public:
   RPageAllocatorFile();
   static RPage NewPage(ColumnId_t columnId, void *mem, std::size_t elementSize, std::size_t nElements);
   RPage NewPage(ColumnId_t columnId, std::size_t elementSize, std::size_t nElements);
   void DeletePage(const RPage& page);
};


// clang-format off
/**
\class ROOT::Experimental::Detail::RPageSourceFile
\ingroup NTuple
\brief Storage provider that reads ntuple pages from a file
*/
// clang-format on
class RPageSourceFile : public RPageSource {
public:
   /// Cannot process pages larger than 1MB
   static constexpr std::size_t kMaxPageSize = 1024 * 1024;

private:
   /// Populated pages might be shared; there memory buffer is managed by the RPageAllocatorFile
   std::unique_ptr<RPageAllocatorFile> fPageAllocator;
   /// The page pool might, at some point, be used by multiple page sources
   std::shared_ptr<RPagePool> fPagePool;
   /// The last cluster from which a page got populated.  Points into fClusterPool->fPool
   RCluster *fCurrentCluster = nullptr;
   /// An RRawFile is used to request the necessary byte ranges from a local or a remote file
   std::unique_ptr<ROOT::Internal::RRawFile> fFile;
   /// Takes the fFile to read ntuple blobs from it
   Internal::RMiniFileReader fReader;
   /// The cluster pool asynchronously preloads the next few clusters
   std::unique_ptr<RClusterPool> fClusterPool;

   RPageSourceFile(std::string_view ntupleName, const RNTupleReadOptions &options);
   RPage PopulatePageFromCluster(ColumnHandle_t columnHandle, const RClusterDescriptor &clusterDescriptor,
                                 ClusterSize_t::ValueType idxInCluster);

protected:
   RNTupleDescriptor AttachImpl() final;
   void UnzipClusterImpl(RCluster *cluster) final;

public:
   RPageSourceFile(std::string_view ntupleName, std::string_view path, const RNTupleReadOptions &options);
   /// The cloned page source creates a new raw file and reader and opens its own file descriptor to the data.
   /// The meta-data (header and footer) is reread and parsed by the clone.
   std::unique_ptr<RPageSource> Clone() const final;

   RPageSourceFile(const RPageSourceFile&) = delete;
   RPageSourceFile& operator=(const RPageSourceFile&) = delete;
   RPageSourceFile(RPageSourceFile&&) = default;
   RPageSourceFile& operator=(RPageSourceFile&&) = default;
   virtual ~RPageSourceFile();

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
