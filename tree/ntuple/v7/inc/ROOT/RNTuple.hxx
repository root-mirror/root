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
#include <ROOT/RSpan.hxx>
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


#ifdef R__USE_IMT
class TTaskGroup;
class RNTupleImtTaskScheduler : public Detail::RPageStorage::RTaskScheduler {
private:
   std::unique_ptr<TTaskGroup> fTaskGroup;
public:
   RNTupleImtTaskScheduler();
   virtual ~RNTupleImtTaskScheduler() = default;
   void Reset() final;
   void AddTask(const std::function<void(void)> &taskFunc) final;
   void Wait() final;
};
#endif

// clang-format off
/**
\class ROOT::Experimental::RNTupleReader
\ingroup NTuple
\brief An RNTuple that is used to read data from storage

An input ntuple provides data from storage as C++ objects. The ntuple model can be created from the data on storage
or it can be imposed by the user. The latter case allows users to read into a specialized ntuple model that covers
only a subset of the fields in the ntuple. The ntuple model is used when reading complete entries.
Individual fields can be read as well by instantiating a tree view.

~~~ {.cpp}
#include <ROOT/RNTuple.hxx>
using ROOT::Experimental::RNTupleReader;

#include <iostream>

auto ntuple = RNTupleReader::Open("myNTuple", "some/file.root");
std::cout << "myNTuple has " << ntuple->GetNEntries() << " entries\n";
~~~
*/
// clang-format on
class RNTupleReader {
private:
   /// Set as the page source's scheduler for parallel page decompression if IMT is on
   /// Needs to be destructed after the pages source is destructed (an thus be declared before)
   std::unique_ptr<Detail::RPageStorage::RTaskScheduler> fUnzipTasks;

   std::unique_ptr<Detail::RPageSource> fSource;
   /// Needs to be destructed before fSource
   std::unique_ptr<RNTupleModel> fModel;
   /// We use a dedicated on-demand reader for Show() and Scan(). Printing data uses all the fields
   /// from the full model even if the analysis code uses only a subset of fields. The display reader
   /// is a clone of the original reader.
   std::unique_ptr<RNTupleReader> fDisplayReader;
   Detail::RNTupleMetrics fMetrics;

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

   /// Used to specify the underlying RNTuples in OpenFriends()
   struct ROpenSpec {
      std::string fNTupleName;
      std::string fStorage;
      RNTupleReadOptions fOptions;

      ROpenSpec() = default;
      ROpenSpec(std::string_view n, std::string_view s) : fNTupleName(n), fStorage(s) {}
   };

   /// Throws an exception if the model is null.
   static std::unique_ptr<RNTupleReader> Open(std::unique_ptr<RNTupleModel> model,
                                              std::string_view ntupleName,
                                              std::string_view storage,
                                              const RNTupleReadOptions &options = RNTupleReadOptions());
   /// Open an RNTuple for reading.
   ///
   /// Throws an RException if there is no RNTuple with the given name.
   ///
   /// **Example: open an RNTuple and print the number of entries**
   /// ~~~ {.cpp}
   /// #include <ROOT/RNTuple.hxx>
   /// using ROOT::Experimental::RNTupleReader;
   ///
   /// #include <iostream>
   ///
   /// auto ntuple = RNTupleReader::Open("myNTuple", "some/file.root");
   /// std::cout << "myNTuple has " << ntuple->GetNEntries() << " entries\n";
   /// ~~~
   static std::unique_ptr<RNTupleReader> Open(std::string_view ntupleName,
                                              std::string_view storage,
                                              const RNTupleReadOptions &options = RNTupleReadOptions());
   /// Open RNTuples as one virtual, horizontally combined ntuple.  The underlying RNTuples must
   /// have an identical number of entries.  Fields in the combined RNTuple are named with the ntuple name
   /// as a prefix, e.g. myNTuple1.px and myNTuple2.pt (see tutorial ntpl006_friends)
   static std::unique_ptr<RNTupleReader> OpenFriends(std::span<ROpenSpec> ntuples);

   /// The user imposes an ntuple model, which must be compatible with the model found in the data on
   /// storage.
   ///
   /// Throws an exception if the model or the source is null.
   RNTupleReader(std::unique_ptr<RNTupleModel> model, std::unique_ptr<Detail::RPageSource> source);
   /// The model is generated from the ntuple metadata on storage
   ///
   /// Throws an exception if the source is null.
   explicit RNTupleReader(std::unique_ptr<Detail::RPageSource> source);
   std::unique_ptr<RNTupleReader> Clone() { return std::make_unique<RNTupleReader>(fSource->Clone()); }
   ~RNTupleReader();

   RNTupleModel *GetModel();
   NTupleSize_t GetNEntries() const { return fSource->GetNEntries(); }
   const RNTupleDescriptor &GetDescriptor() const { return fSource->GetDescriptor(); }

   /// Prints a detailed summary of the ntuple, including a list of fields.
   ///
   /// **Example: print summary information to stdout**
   /// ~~~ {.cpp}
   /// #include <ROOT/RNTuple.hxx>
   /// using ROOT::Experimental::ENTupleInfo;
   /// using ROOT::Experimental::RNTupleReader;
   ///
   /// #include <iostream>
   ///
   /// auto ntuple = RNTupleReader::Open("myNTuple", "some/file.root");
   /// ntuple->PrintInfo();
   /// // or, equivalently:
   /// ntuple->PrintInfo(ENTupleInfo::kSummary, std::cout);
   /// ~~~
   /// **Example: print detailed column storage data to stderr**
   /// ~~~ {.cpp}
   /// #include <ROOT/RNTuple.hxx>
   /// using ROOT::Experimental::ENTupleInfo;
   /// using ROOT::Experimental::RNTupleReader;
   ///
   /// #include <iostream>
   ///
   /// auto ntuple = RNTupleReader::Open("myNTuple", "some/file.root");
   /// ntuple->PrintInfo(ENTupleInfo::kStorageDetails, std::cerr);
   /// ~~~
   ///
   /// For use of ENTupleInfo::kMetrics, see #EnableMetrics.
   void PrintInfo(const ENTupleInfo what = ENTupleInfo::kSummary, std::ostream &output = std::cout);

   /// Shows the values of the i-th entry/row, starting with 0 for the first entry. By default,
   /// prints the output in JSON format.
   /// Uses the visitor pattern to traverse through each field of the given entry.
   void Show(NTupleSize_t index, const ENTupleShowFormat format = ENTupleShowFormat::kCurrentModelJSON,
             std::ostream &output = std::cout);

   /// Analogous to Fill(), fills the default entry of the model. Returns false at the end of the ntuple.
   /// On I/O errors, raises an exception.
   void LoadEntry(NTupleSize_t index) {
      // TODO(jblomer): can be templated depending on the factory method / constructor
      if (R__unlikely(!fModel)) {
         fModel = fSource->GetDescriptor().GenerateModel();
         ConnectModel(*fModel);
      }
      LoadEntry(index, *fModel->GetDefaultEntry());
   }
   /// Fills a user provided entry after checking that the entry has been instantiated from the ntuple model
   void LoadEntry(NTupleSize_t index, REntry &entry) {
      for (auto& value : entry) {
         value.GetField()->Read(index, &value);
      }
   }

   /// Returns an iterator over the entry indices of the RNTuple.
   ///
   /// **Example: iterate over all entries and print each entry in JSON format**
   /// ~~~ {.cpp}
   /// #include <ROOT/RNTuple.hxx>
   /// using ROOT::Experimental::ENTupleShowFormat;
   /// using ROOT::Experimental::RNTupleReader;
   ///
   /// #include <iostream>
   ///
   /// auto ntuple = RNTupleReader::Open("myNTuple", "some/file.root");
   /// for (auto i : ntuple->GetEntryRange()) {
   ///    ntuple->Show(i, ENTupleShowFormat::kCompleteJSON);
   /// }
   /// ~~~
   RNTupleGlobalRange GetEntryRange() { return RNTupleGlobalRange(0, GetNEntries()); }

   /// Provides access to an individual field that can contain either a scalar value or a collection, e.g.
   /// GetView<double>("particles.pt") or GetView<std::vector<double>>("particle").  It can as well be the index
   /// field of a collection itself, like GetView<NTupleSize_t>("particle").
   ///
   /// Raises an exception if there is no field with the given name.
   ///
   /// **Example: iterate over a field named "pt" of type `float`**
   /// ~~~ {.cpp}
   /// #include <ROOT/RNTuple.hxx>
   /// using ROOT::Experimental::RNTupleReader;
   ///
   /// #include <iostream>
   ///
   /// auto ntuple = RNTupleReader::Open("myNTuple", "some/file.root");
   /// auto pt = ntuple->GetView<float>("pt");
   ///
   /// for (auto i : ntuple->GetEntryRange()) {
   ///    std::cout << i << ": " << pt(i) << "\n";
   /// }
   /// ~~~
   template <typename T>
   RNTupleView<T> GetView(std::string_view fieldName) {
      auto fieldId = fSource->GetDescriptor().FindFieldId(fieldName);
      if (fieldId == kInvalidDescriptorId) {
         throw RException(R__FAIL("no field named '" + std::string(fieldName) + "' in RNTuple '"
            + fSource->GetDescriptor().GetName() + "'"
         ));
      }
      return RNTupleView<T>(fieldId, fSource.get());
   }

   /// Raises an exception if:
   /// * there is no field with the given name or,
   /// * the field is not a collection
   RNTupleViewCollection GetViewCollection(std::string_view fieldName) {
      auto fieldId = fSource->GetDescriptor().FindFieldId(fieldName);
      if (fieldId == kInvalidDescriptorId) {
         throw RException(R__FAIL("no field named '" + std::string(fieldName) + "' in RNTuple '"
            + fSource->GetDescriptor().GetName() + "'"
         ));
      }
      return RNTupleViewCollection(fieldId, fSource.get());
   }

   RIterator begin() { return RIterator(0); }
   RIterator end() { return RIterator(GetNEntries()); }

   /// Enable performance measurements (decompression time, bytes read from storage, etc.)
   ///
   /// **Example: inspect the reader metrics after loading every entry**
   /// ~~~ {.cpp}
   /// #include <ROOT/RNTuple.hxx>
   /// using ROOT::Experimental::ENTupleInfo;
   /// using ROOT::Experimental::RNTupleReader;
   ///
   /// #include <iostream>
   ///
   /// auto ntuple = RNTupleReader::Open("myNTuple", "some/file.root");
   /// // metrics must be turned on beforehand
   /// ntuple->EnableMetrics();
   ///
   /// for (auto i : ntuple->GetEntryRange()) {
   ///    ntuple->LoadEntry(i);
   /// }
   /// ntuple->PrintInfo(ENTupleInfo::kMetrics);
   /// ~~~
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
   /// The page sink's parallel page compression scheduler if IMT is on.
   /// Needs to be destructed after the page sink is destructed and so declared before.
   std::unique_ptr<Detail::RPageStorage::RTaskScheduler> fZipTasks;
   std::unique_ptr<Detail::RPageSink> fSink;
   /// Needs to be destructed before fSink
   std::unique_ptr<RNTupleModel> fModel;
   Detail::RNTupleMetrics fMetrics;
   NTupleSize_t fLastCommitted;
   NTupleSize_t fNEntries;

public:
   /// Throws an exception if the model is null.
   static std::unique_ptr<RNTupleWriter> Recreate(std::unique_ptr<RNTupleModel> model,
                                                  std::string_view ntupleName,
                                                  std::string_view storage,
                                                  const RNTupleWriteOptions &options = RNTupleWriteOptions());
   /// Throws an exception if the model is null.
   static std::unique_ptr<RNTupleWriter> Append(std::unique_ptr<RNTupleModel> model,
                                                std::string_view ntupleName,
                                                TFile &file,
                                                const RNTupleWriteOptions &options = RNTupleWriteOptions());
   /// Throws an exception if the model or the sink is null.
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
      if ((fNEntries % fSink->GetWriteOptions().GetNEntriesPerCluster()) == 0)
         CommitCluster();
   }
   /// Ensure that the data from the so far seen Fill calls has been written to storage
   void CommitCluster();

   std::unique_ptr<REntry> CreateEntry() { return fModel->CreateEntry(); }

   void EnableMetrics() { fMetrics.Enable(); }
   const Detail::RNTupleMetrics &GetMetrics() const { return fMetrics; }
};

// clang-format off
/**
\class ROOT::Experimental::RCollectionNTuple
\ingroup NTuple
\brief A virtual ntuple used for writing untyped collections that can be used to some extent like an RNTupleWriter
*
* This class is between a field and a ntuple.  It carries the offset column for the collection and the default entry
* taken from the collection model.  It does not, however, own an ntuple model because the collection model has been
* merged into the larger ntuple model.
*/
// clang-format on
class RCollectionNTupleWriter {
private:
   ClusterSize_t fOffset;
   std::unique_ptr<REntry> fDefaultEntry;
public:
   explicit RCollectionNTupleWriter(std::unique_ptr<REntry> defaultEntry);
   RCollectionNTupleWriter(const RCollectionNTupleWriter&) = delete;
   RCollectionNTupleWriter& operator=(const RCollectionNTupleWriter&) = delete;
   ~RCollectionNTupleWriter() = default;

   void Fill() { Fill(fDefaultEntry.get()); }
   void Fill(REntry *entry) {
      for (auto &value : *entry) {
         value.GetField()->Append(value);
      }
      fOffset++;
   }

   ClusterSize_t *GetOffsetPtr() { return &fOffset; }
};

} // namespace Experimental
} // namespace ROOT

#endif
