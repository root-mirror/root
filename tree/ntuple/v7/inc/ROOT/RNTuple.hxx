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

#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RNTupleUtil.hxx>
#include <ROOT/RNTupleView.hxx>
#include <ROOT/RStringView.hxx>

#include <iterator>
#include <memory>
#include <utility>

namespace ROOT {
namespace Experimental {

class REntry;
class RNTupleModel;

namespace Detail {
class RPageSink;
class RPageSource;
}

namespace Detail {

// clang-format off
/**
\class ROOT::Experimental::RNTuple
\ingroup NTuple
\brief The RNTuple represents a live dataset, whose structure is defined by an RNTupleModel

RNTuple connects the static information of the RNTupleModel to a source or sink on physical storage.
Reading and writing requires use of the corresponding derived class RInputForest or ROutputForest.
RNTuple writes only complete entries (rows of the data set).  The entry itself is not kept within the
RNTuple, which allows for multiple concurrent entries for the same RNTuple.  Besides reading an entire entry,
the RNTuple can expose views that read only specific fields.
*/
// clang-format on
class RNTuple {
protected:
   std::unique_ptr<RNTupleModel> fModel;
   /// The number of entries is constant for reading and reflects the sum of Fill() operations when writing
   ForestSize_t fNEntries;

   /// Only the derived RInputForest and ROutputForest can be instantiated
   explicit RNTuple(std::unique_ptr<RNTupleModel> model);

public:
   RNTuple(const RNTuple&) = delete;
   RNTuple& operator =(const RNTuple&) = delete;
   ~RNTuple();

   RNTupleModel* GetModel() { return fModel.get(); }
}; // RNTuple

} // namespace Detail


/**
 * Listing of the different options that can be returned by RInputForest::GetInfo()
 */
enum class ENTupleInfo {
   kSummary,  // The ntuple name, description, number of entries
};


// clang-format off
/**
\class ROOT::Experimental::RInputForest
\ingroup NTuple
\brief An RNTuple that is used to read data from storage

An input ntuple provides data from storage as C++ objects. The ntuple model can be created from the data on storage
or it can be imposed by the user. The latter case allows users to read into a specialized ntuple model that covers
only a subset of the fields in the ntuple. The ntuple model is used when reading complete entries.
Individual fields can be read as well by instantiating a tree view.
*/
// clang-format on
class RInputForest : public Detail::RNTuple {
private:
   std::unique_ptr<Detail::RPageSource> fSource;

public:
   // Browse through the entries
   class RIterator : public std::iterator<std::forward_iterator_tag, ForestSize_t> {
   private:
      using iterator = RIterator;
      ForestSize_t fIndex = kInvalidForestIndex;
   public:
      RIterator() = default;
      explicit RIterator(ForestSize_t index) : fIndex(index) {}
      ~RIterator() = default;

      iterator  operator++(int) /* postfix */        { auto r = *this; ++fIndex; return r; }
      iterator& operator++()    /* prefix */         { ++fIndex; return *this; }
      reference operator* ()                         { return fIndex; }
      pointer   operator->()                         { return &fIndex; }
      bool      operator==(const iterator& rh) const { return fIndex == rh.fIndex; }
      bool      operator!=(const iterator& rh) const { return fIndex != rh.fIndex; }
   };


   static std::unique_ptr<RInputForest> Open(std::unique_ptr<RNTupleModel> model,
                                             std::string_view forestName,
                                             std::string_view storage);
   static std::unique_ptr<RInputForest> Open(std::string_view forestName, std::string_view storage);

   /// The user imposes an ntuple model, which must be compatible with the model found in the data on storage
   RInputForest(std::unique_ptr<RNTupleModel> model, std::unique_ptr<Detail::RPageSource> source);
   /// The model is generated from the ntuple metadata on storage
   RInputForest(std::unique_ptr<Detail::RPageSource> source);
   ~RInputForest();

   ForestSize_t GetNEntries() { return fNEntries; }

   std::string GetInfo(const ENTupleInfo what = ENTupleInfo::kSummary);

   /// Analogous to Fill(), fills the default entry of the model. Returns false at the end of the ntuple.
   /// On I/O errors, raises an expection.
   void LoadEntry(ForestSize_t index) { LoadEntry(index, fModel->GetDefaultEntry()); }
   /// Fills a user provided entry after checking that the entry has been instantiated from the ntuple model
   void LoadEntry(ForestSize_t index, REntry* entry) {
      for (auto& value : *entry) {
         value.GetField()->Read(index, &value);
      }
   }

   RNTupleViewRange GetViewRange() { return RNTupleViewRange(0, fNEntries); }

   /// Provides access to an individual field that can contain either a skalar value or a collection, e.g.
   /// GetView<double>("particles.pt") or GetView<std::vector<double>>("particle").  It can as well be the index
   /// field of a collection itself, like GetView<ForestSize_t>("particle")
   template <typename T>
   RNTupleView<T> GetView(std::string_view fieldName) { return RNTupleView<T>(fieldName, fSource.get()); }
   RNTupleViewCollection GetViewCollection(std::string_view fieldName) {
      return RNTupleViewCollection(fieldName, fSource.get());
   }

   RIterator begin() { return RIterator(0); }
   RIterator end() { return RIterator(fNEntries); }
};

// clang-format off
/**
\class ROOT::Experimental::ROutputForest
\ingroup NTuple
\brief An RNTuple that gets filled with entries (data) and writes them to storage

An output ntuple can be filled with entries. The caller has to make sure that the data that gets filled into an ntuple
is not modified for the time of the Fill() call. The fill call serializes the C++ object into the column format and
writes data into the corresponding column page buffers.  Writing of the buffers to storage is deferred and can be
triggered by Flush() or by destructing the ntuple.  On I/O errors, an exception is thrown.
*/
// clang-format on
class ROutputForest : public Detail::RNTuple {
private:
   static constexpr ForestSize_t kDefaultClusterSizeEntries = 8192;
   std::unique_ptr<Detail::RPageSink> fSink;
   ForestSize_t fClusterSizeEntries;
   ForestSize_t fLastCommitted;

public:
   static std::unique_ptr<ROutputForest> Recreate(std::unique_ptr<RNTupleModel> model,
                                                  std::string_view forestName,
                                                  std::string_view storage);
   ROutputForest(std::unique_ptr<RNTupleModel> model, std::unique_ptr<Detail::RPageSink> sink);
   ROutputForest(const ROutputForest&) = delete;
   ROutputForest& operator=(const ROutputForest&) = delete;
   ~ROutputForest();

   /// The simplest user interface if the default entry that comes with the ntuple model is used
   void Fill() { Fill(fModel->GetDefaultEntry()); }
   /// Multiple entries can have been instantiated from the tnuple model.  This method will perform
   /// a light check whether the entry comes from the ntuple's own model
   void Fill(REntry *entry) {
      for (auto& value : *entry) {
         value.GetField()->Append(value);
      }
      fNEntries++;
      if ((fNEntries % fClusterSizeEntries) == 0) CommitCluster();
   }
   /// Ensure that the data from the so far seen Fill calls has been written to storage
   void CommitCluster();
};

// clang-format off
/**
\class ROOT::Experimental::RCollectionForest
\ingroup NTuple
\brief A virtual ntuple for collections that can be used to some extent like a real ntuple
*
* This class is between a field and a ntuple.  It carries the offset column for the collection and the default entry
* taken from the collection model.  It does not, however, have a tree model because the collection model has been merged
* into the larger ntuple model.
*/
// clang-format on
class RCollectionForest {
private:
   ClusterSize_t fOffset;
   std::unique_ptr<REntry> fDefaultEntry;
public:
   explicit RCollectionForest(std::unique_ptr<REntry> defaultEntry);
   RCollectionForest(const RCollectionForest&) = delete;
   RCollectionForest& operator=(const RCollectionForest&) = delete;
   ~RCollectionForest() = default;

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
