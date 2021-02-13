/// \file ROOT/RPageStorage.hxx
/// \ingroup NTuple ROOT7
/// \author Jakob Blomer <jblomer@cern.ch>
/// \date 2018-07-19
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RPageStorage
#define ROOT7_RPageStorage

#include <ROOT/RNTupleDescriptor.hxx>
#include <ROOT/RNTupleOptions.hxx>
#include <ROOT/RNTupleUtil.hxx>
#include <ROOT/RPage.hxx>
#include <ROOT/RPageAllocator.hxx>
#include <ROOT/RStringView.hxx>

#include <atomic>
#include <cstddef>
#include <functional>
#include <memory>
#include <unordered_set>

namespace ROOT {
namespace Experimental {

class RNTupleModel;
// TODO(jblomer): factory methods to create tree sinks and sources outside Detail namespace

namespace Detail {

class RCluster;
class RColumn;
class RPagePool;
class RFieldBase;
class RNTupleMetrics;

enum class EPageStorageType {
   kSink,
   kSource,
};

// clang-format off
/**
\class ROOT::Experimental::Detail::RPageStorage
\ingroup NTuple
\brief Common functionality of an ntuple storage for both reading and writing

The RPageStore provides access to a storage container that keeps the bits of pages and clusters comprising
an ntuple.  Concrete implementations can use a TFile, a raw file, an object store, and so on.
*/
// clang-format on
class RPageStorage {
public:
   /// The interface of a task scheduler to schedule page (de)compression tasks
   class RTaskScheduler {
   public:
      /// Start a new set of tasks
      virtual void Reset() = 0;
      /// Take a callable that represents a task
      virtual void AddTask(const std::function<void(void)> &taskFunc) = 0;
      /// Blocks until all scheduled tasks finished
      virtual void Wait() = 0;
   };

protected:
   std::string fNTupleName;
   RTaskScheduler *fTaskScheduler = nullptr;

public:
   explicit RPageStorage(std::string_view name);
   RPageStorage(const RPageStorage &other) = delete;
   RPageStorage& operator =(const RPageStorage &other) = delete;
   RPageStorage(RPageStorage &&other) = default;
   RPageStorage& operator =(RPageStorage &&other) = default;
   virtual ~RPageStorage();

   /// Whether the concrete implementation is a sink or a source
   virtual EPageStorageType GetType() = 0;

   struct RColumnHandle {
      DescriptorId_t fId = kInvalidDescriptorId;
      const RColumn *fColumn = nullptr;

      /// Returns true for a valid column handle; fColumn and fId should always either both
      /// be valid or both be invalid.
      operator bool() const { return fId != kInvalidDescriptorId && fColumn; }
   };
   /// The column handle identifies a column with the current open page storage
   using ColumnHandle_t = RColumnHandle;

   /// Register a new column.  When reading, the column must exist in the ntuple on disk corresponding to the meta-data.
   /// When writing, every column can only be attached once.
   virtual ColumnHandle_t AddColumn(DescriptorId_t fieldId, const RColumn &column) = 0;
   /// Unregisters a column.  A page source decreases the reference counter for the corresponding active column.
   /// For a page sink, dropping columns is currently a no-op.
   virtual void DropColumn(ColumnHandle_t columnHandle) = 0;

   /// Every page store needs to be able to free pages it handed out.  But Sinks and sources have different means
   /// of allocating pages.
   virtual void ReleasePage(RPage &page) = 0;

   /// Returns an empty metrics.  Page storage implementations usually have their own metrics.
   virtual RNTupleMetrics &GetMetrics();

   void SetTaskScheduler(RTaskScheduler *taskScheduler) { fTaskScheduler = taskScheduler; }
};

// clang-format off
/**
\class ROOT::Experimental::Detail::RPageSink
\ingroup NTuple
\brief Abstract interface to write data into an ntuple

The page sink takes the list of columns and afterwards a series of page commits and cluster commits.
The user is responsible to commit clusters at a consistent point, i.e. when all pages corresponding to data
up to the given entry number are committed.
*/
// clang-format on
class RPageSink : public RPageStorage {
protected:
   RNTupleWriteOptions fOptions;

   /// Building the ntuple descriptor while writing is done in the same way for all the storage sink implementations.
   /// Field, column, cluster ids and page indexes per cluster are issued sequentially starting with 0
   DescriptorId_t fLastFieldId = 0;
   DescriptorId_t fLastColumnId = 0;
   DescriptorId_t fLastClusterId = 0;
   NTupleSize_t fPrevClusterNEntries = 0;
   /// Keeps track of the number of elements in the currently open cluster. Indexed by column id.
   std::vector<RClusterDescriptor::RColumnRange> fOpenColumnRanges;
   /// Keeps track of the written pages in the currently open cluster. Indexed by column id.
   std::vector<RClusterDescriptor::RPageRange> fOpenPageRanges;
   RNTupleDescriptorBuilder fDescriptorBuilder;

   virtual void CreateImpl(const RNTupleModel &model) = 0;
   virtual RClusterDescriptor::RLocator CommitPageImpl(ColumnHandle_t columnHandle, const RPage &page) = 0;
   virtual RClusterDescriptor::RLocator CommitClusterImpl(NTupleSize_t nEntries) = 0;
   virtual void CommitDatasetImpl() = 0;

public:
   RPageSink(std::string_view ntupleName, const RNTupleWriteOptions &options);

   RPageSink(const RPageSink&) = delete;
   RPageSink& operator=(const RPageSink&) = delete;
   RPageSink(RPageSink&&) = default;
   RPageSink& operator=(RPageSink&&) = default;
   virtual ~RPageSink();

   /// Guess the concrete derived page source from the file name (location)
   static std::unique_ptr<RPageSink> Create(std::string_view ntupleName, std::string_view location,
                                            const RNTupleWriteOptions &options = RNTupleWriteOptions());
   EPageStorageType GetType() final { return EPageStorageType::kSink; }

   ColumnHandle_t AddColumn(DescriptorId_t fieldId, const RColumn &column) final;
   void DropColumn(ColumnHandle_t /*columnHandle*/) final {}

   /// Physically creates the storage container to hold the ntuple (e.g., a keys a TFile or an S3 bucket)
   /// To do so, Create() calls CreateImpl() after updating the descriptor.
   /// Create() associates column handles to the columns referenced by the model
   void Create(RNTupleModel &model);
   /// Write a page to the storage. The column must have been added before.
   void CommitPage(ColumnHandle_t columnHandle, const RPage &page);
   /// Finalize the current cluster and create a new one for the following data.
   void CommitCluster(NTupleSize_t nEntries);
   /// Finalize the current cluster and the entrire data set.
   void CommitDataset() { CommitDatasetImpl(); }

   /// Get a new, empty page for the given column that can be filled with up to nElements.  If nElements is zero,
   /// the page sink picks an appropriate size.
   virtual RPage ReservePage(ColumnHandle_t columnHandle, std::size_t nElements = 0) = 0;
};

// clang-format off
/**
\class ROOT::Experimental::Detail::RPageSource
\ingroup NTuple
\brief Abstract interface to read data from an ntuple

The page source is initialized with the columns of interest. Pages from those columns can then be
mapped into memory. The page source also gives access to the ntuple's meta-data.
*/
// clang-format on
class RPageSource : public RPageStorage {
public:
   /// Derived from the model (fields) that are actually being requested at a given point in time
   using ColumnSet_t = std::unordered_set<DescriptorId_t>;

protected:
   RNTupleReadOptions fOptions;
   RNTupleDescriptor fDescriptor;
   /// The active columns are implicitly defined by the model fields or views
   ColumnSet_t fActiveColumns;

   virtual RNTupleDescriptor AttachImpl() = 0;
   // Only called if a task scheduler is set. No-op be default.
   virtual void UnzipClusterImpl(RCluster * /* cluster */)
      { }

public:
   RPageSource(std::string_view ntupleName, const RNTupleReadOptions &fOptions);
   RPageSource(const RPageSource&) = delete;
   RPageSource& operator=(const RPageSource&) = delete;
   RPageSource(RPageSource&&) = default;
   RPageSource& operator=(RPageSource&&) = default;
   virtual ~RPageSource();
   /// Guess the concrete derived page source from the file name (location)
   static std::unique_ptr<RPageSource> Create(std::string_view ntupleName, std::string_view location,
                                              const RNTupleReadOptions &options = RNTupleReadOptions());
   /// Open the same storage multiple time, e.g. for reading in multiple threads
   virtual std::unique_ptr<RPageSource> Clone() const = 0;

   EPageStorageType GetType() final { return EPageStorageType::kSource; }
   const RNTupleDescriptor &GetDescriptor() const { return fDescriptor; }
   ColumnHandle_t AddColumn(DescriptorId_t fieldId, const RColumn &column) final;
   void DropColumn(ColumnHandle_t columnHandle) final;

   /// Open the physical storage container for the tree
   void Attach() { fDescriptor = AttachImpl(); }
   NTupleSize_t GetNEntries();
   NTupleSize_t GetNElements(ColumnHandle_t columnHandle);
   ColumnId_t GetColumnId(ColumnHandle_t columnHandle);

   /// Allocates and fills a page that contains the index-th element
   virtual RPage PopulatePage(ColumnHandle_t columnHandle, NTupleSize_t globalIndex) = 0;
   /// Another version of PopulatePage that allows to specify cluster-relative indexes
   virtual RPage PopulatePage(ColumnHandle_t columnHandle, const RClusterIndex &clusterIndex) = 0;

   /// Populates all the pages of the given cluster id and columns; it is possible that some columns do not
   /// contain any pages.  The pages source may load more columns than the minimal necessary set from `columns`.
   /// To indicate which columns have been loaded, LoadCluster() must mark them with SetColumnAvailable().
   /// That includes the ones from the `columns` that don't have pages; otherwise subsequent requests
   /// for the cluster would assume an incomplete cluster and trigger loading again.
   /// LoadCluster() is typically called from the I/O thread of a cluster pool, i.e. the method runs
   /// concurrently to other methods of the page source.
   virtual std::unique_ptr<RCluster> LoadCluster(DescriptorId_t clusterId, const ColumnSet_t &columns) = 0;

   /// Parallel decompression and unpacking of the pages in the given cluster. The unzipped pages are supposed
   /// to be preloaded in a page pool attached to the source. The method is triggered by the cluster pool's
   /// unzip thread. It is an optional optimization, the method can safely do nothing. In particular, the
   /// actual implementation will only run if a task scheduler is set. In practice, a task scheduler is set
   /// if implicit multi-threading is turned on.
   void UnzipCluster(RCluster *cluster);
};

} // namespace Detail

} // namespace Experimental
} // namespace ROOT

#endif
