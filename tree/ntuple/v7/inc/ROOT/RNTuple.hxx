/// \file ROOT/RNTuple.hxx
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

#ifndef ROOT7_RNTuple
#define ROOT7_RNTuple

#include <ROOT/RConfig.hxx> // for R__unlikely
#include <ROOT/RNTupleMetrics.hxx>
#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RNTupleOptions.hxx>
#include <ROOT/RNTupleUtil.hxx>
#include <ROOT/RNTupleView.hxx>
#include <ROOT/RPageStorage.hxx>
#include <ROOT/RStringView.hxx>

#include <iterator>
#include <memory>
#include <sstream>
#include <utility>

class TFile;

namespace ROOT {
namespace Experimental {

class REntry;
class RNTupleModel;
class TTaskGroup;

namespace Detail {
class RPageSink;
class RPageSource;
}


/**
 * Listing of the different options that can be printed by RNTupleReader::GetInfo()
 */
enum class ENTupleInfo {
   kSummary,  // The ntuple name, description, number of entries
   kStorageDetails, // size on storage, page sizes, compression factor, etc.
   kMetrics, // internals performance counters, requires that EnableMetrics() was called
};

/**
 * Listing of the different entry output formats of RNTupleReader::Show()
 */
enum class ENTupleShowFormat {
   kCurrentModelJSON, // prints a single entry/row with the current active model in JSON format.
   kCompleteJSON,  // prints a single entry/row with all the fields in JSON format.
};


class RNTupleImtTaskScheduler : public Detail::RPageStorage::RTaskScheduler {
private:
   std::unique_ptr<TTaskGroup> fTaskGroup;
public:
   void Reset() final;
   void AddTask(const std::function<void(void)> &taskFunc) final;
   void Wait() final;
};

// clang-format off
/**
\class ROOT::Experimental::RNTupleReader
\ingroup NTuple
\brief An RNTuple that is used to read data from storage

An input ntuple provides data from storage as C++ objects. The ntuple model can be created from the data on storage
or it can be imposed by the user. The latter case allows users to read into a specialized ntuple model that covers
only a subset of the fields in the ntuple. The ntuple model is used when reading complete entries.
Individual fields can be read as well by instantiating a tree view.
*/
// clang-format on
class RNTupleReader {
private:
   std::unique_ptr<Detail::RPageSource> fSource;
   /// Needs to be destructed before fSource
   std::unique_ptr<RNTupleModel> fModel;
   /// We use a dedicated on-demand reader for Show() and Scan(). Printing data uses all the fields
   /// from the full model even if the analysis code uses only a subset of fields. The display reader
   /// is a clone of the original reader.
   std::unique_ptr<RNTupleReader> fDisplayReader;
   Detail::RNTupleMetrics fMetrics;
   /// Set as the page source's scheduler for parallel page decompression if IMT is on
   RNTupleImtTaskScheduler fUnzipTasks;

   void ConnectModel(const RNTupleModel &model);
   RNTupleReader *GetDisplayReader();
   void InitPageSource();

public:
   // Browse through the entries
   class RIterator {
   private:
      NTupleSize_t fIndex = kInvalidNTupleIndex;
   public:
      using iterator = RIterator;
      using iterator_category = std::forward_iterator_tag;
      using value_type = NTupleSize_t;
      using difference_type = NTupleSize_t;
      using pointer = NTupleSize_t*;
      using reference = NTupleSize_t&;

      RIterator() = default;
      explicit RIterator(NTupleSize_t index) : fIndex(index) {}
      ~RIterator() = default;

      iterator  operator++(int) /* postfix */        { auto r = *this; fIndex++; return r; }
      iterator& operator++()    /* prefix */         { ++fIndex; return *this; }
      reference operator* ()                         { return fIndex; }
      pointer   operator->()                         { return &fIndex; }
      bool      operator==(const iterator& rh) const { return fIndex == rh.fIndex; }
      bool      operator!=(const iterator& rh) const { return fIndex != rh.fIndex; }
   };


   static std::unique_ptr<RNTupleReader> Open(std::unique_ptr<RNTupleModel> model,
                                              std::string_view ntupleName,
                                              std::string_view storage,
                                              const RNTupleReadOptions &options = RNTupleReadOptions());
   static std::unique_ptr<RNTupleReader> Open(std::string_view ntupleName,
                                              std::string_view storage,
                                              const RNTupleReadOptions &options = RNTupleReadOptions());

   /// The user imposes an ntuple model, which must be compatible with the model found in the data on storage
   RNTupleReader(std::unique_ptr<RNTupleModel> model, std::unique_ptr<Detail::RPageSource> source);
   /// The model is generated from the ntuple metadata on storage
   explicit RNTupleReader(std::unique_ptr<Detail::RPageSource> source);
   std::unique_ptr<RNTupleReader> Clone() { return std::make_unique<RNTupleReader>(fSource->Clone()); }
   ~RNTupleReader();

   RNTupleModel *GetModel();
   NTupleSize_t GetNEntries() const { return fSource->GetNEntries(); }
   const RNTupleDescriptor &GetDescriptor() const { return fSource->GetDescriptor(); }

   /// Prints a detailed summary of the ntuple, including a list of fields.
   void PrintInfo(const ENTupleInfo what = ENTupleInfo::kSummary, std::ostream &output = std::cout);

   /// Shows the values of the i-th entry/row, starting with 0 for the first entry. By default,
   /// prints the output in JSON format.
   /// Uses the visitor pattern to traverse through each field of the given entry.
   void Show(NTupleSize_t index, const ENTupleShowFormat format = ENTupleShowFormat::kCurrentModelJSON,
             std::ostream &output = std::cout);

   /// Analogous to Fill(), fills the default entry of the model. Returns false at the end of the ntuple.
   /// On I/O errors, raises an exception.
   void LoadEntry(NTupleSize_t index) { LoadEntry(index, *fModel->GetDefaultEntry()); }
   /// Fills a user provided entry after checking that the entry has been instantiated from the ntuple model
   void LoadEntry(NTupleSize_t index, REntry &entry) {
      // TODO(jblomer): can be templated depending on the factory method / constructor
      if (R__unlikely(!fModel)) {
         fModel = fSource->GetDescriptor().GenerateModel();
         ConnectModel(*fModel);
      }

      for (auto& value : entry) {
         value.GetField()->Read(index, &value);
      }
   }

   RNTupleGlobalRange GetEntryRange() { return RNTupleGlobalRange(0, GetNEntries()); }

   /// Provides access to an individual field that can contain either a scalar value or a collection, e.g.
   /// GetView<double>("particles.pt") or GetView<std::vector<double>>("particle").  It can as well be the index
   /// field of a collection itself, like GetView<NTupleSize_t>("particle")
   template <typename T>
   RNTupleView<T> GetView(std::string_view fieldName) {
      auto fieldId = fSource->GetDescriptor().FindFieldId(fieldName);
      return RNTupleView<T>(fieldId, fSource.get());
   }
   RNTupleViewCollection GetViewCollection(std::string_view fieldName) {
      auto fieldId = fSource->GetDescriptor().FindFieldId(fieldName);
      return RNTupleViewCollection(fieldId, fSource.get());
   }

   RIterator begin() { return RIterator(0); }
   RIterator end() { return RIterator(GetNEntries()); }

   void EnableMetrics() { fMetrics.Enable(); }
   const Detail::RNTupleMetrics &GetMetrics() const { return fMetrics; }
};

// clang-format off
/**
\class ROOT::Experimental::RNTupleWriter
\ingroup NTuple
\brief An RNTuple that gets filled with entries (data) and writes them to storage

An output ntuple can be filled with entries. The caller has to make sure that the data that gets filled into an ntuple
is not modified for the time of the Fill() call. The fill call serializes the C++ object into the column format and
writes data into the corresponding column page buffers.  Writing of the buffers to storage is deferred and can be
triggered by Flush() or by destructing the ntuple.  On I/O errors, an exception is thrown.
*/
// clang-format on
class RNTupleWriter {
private:
   static constexpr NTupleSize_t kDefaultClusterSizeEntries = 64000;
   std::unique_ptr<Detail::RPageSink> fSink;
   /// Needs to be destructed before fSink
   std::unique_ptr<RNTupleModel> fModel;
   NTupleSize_t fClusterSizeEntries;
   NTupleSize_t fLastCommitted;
   NTupleSize_t fNEntries;

public:
   static std::unique_ptr<RNTupleWriter> Recreate(std::unique_ptr<RNTupleModel> model,
                                                  std::string_view ntupleName,
                                                  std::string_view storage,
                                                  const RNTupleWriteOptions &options = RNTupleWriteOptions());
   static std::unique_ptr<RNTupleWriter> Append(std::unique_ptr<RNTupleModel> model,
                                                std::string_view ntupleName,
                                                TFile &file,
                                                const RNTupleWriteOptions &options = RNTupleWriteOptions());
   RNTupleWriter(std::unique_ptr<RNTupleModel> model, std::unique_ptr<Detail::RPageSink> sink);
   RNTupleWriter(const RNTupleWriter&) = delete;
   RNTupleWriter& operator=(const RNTupleWriter&) = delete;
   ~RNTupleWriter();

   /// The simplest user interface if the default entry that comes with the ntuple model is used
   void Fill() { Fill(*fModel->GetDefaultEntry()); }
   /// Multiple entries can have been instantiated from the tnuple model.  This method will perform
   /// a light check whether the entry comes from the ntuple's own model
   void Fill(REntry &entry) {
      for (auto& value : entry) {
         value.GetField()->Append(value);
      }
      fNEntries++;
      if ((fNEntries % fClusterSizeEntries) == 0)
         CommitCluster();
   }
   /// Ensure that the data from the so far seen Fill calls has been written to storage
   void CommitCluster();
};

// clang-format off
/**
\class ROOT::Experimental::RCollectionNTuple
\ingroup NTuple
\brief A virtual ntuple for collections that can be used to some extent like a real ntuple
*
* This class is between a field and a ntuple.  It carries the offset column for the collection and the default entry
* taken from the collection model.  It does not, however, have a tree model because the collection model has been merged
* into the larger ntuple model.
*/
// clang-format on
class RCollectionNTuple {
private:
   ClusterSize_t fOffset;
   std::unique_ptr<REntry> fDefaultEntry;
public:
   explicit RCollectionNTuple(std::unique_ptr<REntry> defaultEntry);
   RCollectionNTuple(const RCollectionNTuple&) = delete;
   RCollectionNTuple& operator=(const RCollectionNTuple&) = delete;
   ~RCollectionNTuple() = default;

   void Fill() { Fill(fDefaultEntry.get()); }
   void Fill(REntry *entry) {
      for (auto& treeValue : *entry) {
         treeValue.GetField()->Append(treeValue);
      }
      fOffset++;
   }

   ClusterSize_t* GetOffsetPtr() { return &fOffset; }
};

} // namespace Experimental
} // namespace ROOT

#endif
