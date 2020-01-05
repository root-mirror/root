/// \file ROOT/RNTupleView.hxx
/// \ingroup NTuple ROOT7
/// \author Jakob Blomer <jblomer@cern.ch>
/// \date 2018-10-05
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RNTupleView
#define ROOT7_RNTupleView

#include <ROOT/RField.hxx>
#include <ROOT/RNTupleUtil.hxx>
#include <ROOT/RStringView.hxx>

#include <iterator>
#include <memory>
#include <utility>
#include <unordered_map>

namespace ROOT {
namespace Experimental {


// clang-format off
/**
\class ROOT::Experimental::RNTupleGlobalRange
\ingroup NTuple
\brief Used to loop over indexes (entries or collections) between start and end
*/
// clang-format on
class RNTupleGlobalRange {
private:
   const NTupleSize_t fStart;
   const NTupleSize_t fEnd;
public:
   class RIterator : public std::iterator<std::forward_iterator_tag, NTupleSize_t> {
   private:
      using iterator = RIterator;
      NTupleSize_t fIndex = kInvalidNTupleIndex;
   public:
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

   RNTupleGlobalRange(NTupleSize_t start, NTupleSize_t end) : fStart(start), fEnd(end) {}
   RIterator begin() { return RIterator(fStart); }
   RIterator end() { return RIterator(fEnd); }
};


// clang-format off
/**
\class ROOT::Experimental::RNTupleClusterRange
\ingroup NTuple
\brief Used to loop over entries of collections in a single cluster
*/
// clang-format on
class RNTupleClusterRange {
private:
   const DescriptorId_t fClusterId;
   const ClusterSize_t::ValueType fStart;
   const ClusterSize_t::ValueType fEnd;
public:
   class RIterator : public std::iterator<std::forward_iterator_tag, RClusterIndex> {
   private:
      using iterator = RIterator;
      RClusterIndex fIndex;
   public:
      RIterator() = default;
      explicit RIterator(const RClusterIndex &index) : fIndex(index) {}
      ~RIterator() = default;

      iterator  operator++(int) /* postfix */        { auto r = *this; fIndex++; return r; }
      iterator& operator++()    /* prefix */         { fIndex++; return *this; }
      reference operator* ()                         { return fIndex; }
      pointer   operator->()                         { return &fIndex; }
      bool      operator==(const iterator& rh) const { return fIndex == rh.fIndex; }
      bool      operator!=(const iterator& rh) const { return fIndex != rh.fIndex; }
   };

   RNTupleClusterRange(DescriptorId_t clusterId, ClusterSize_t::ValueType start, ClusterSize_t::ValueType end)
      : fClusterId(clusterId), fStart(start), fEnd(end) {}
   RIterator begin() { return RIterator(RClusterIndex(fClusterId, fStart)); }
   RIterator end() { return RIterator(RClusterIndex(fClusterId, fEnd)); }
};


// clang-format off
/**
\class ROOT::Experimental::RNTupleView
\ingroup NTuple
\brief An RNTupleView provides read-only access to a single field of the ntuple

(NB(jblomer): The ntuple view is very close to TTreeReader. Do we simply want to teach TTreeReader to deal with
RNTuple?)

The view owns a field and its underlying columns in order to fill an ntuple value object with data. Data can be
accessed by index. For top level fields, the index refers to the entry number. Fields that are part of
nested collections have global index numbers that are derived from their parent indexes.

The RNTupleView object is an iterable. That means, all field values in the tree can be sequentially read from begin()
to end().

For simple types, template specializations let the reading become a pure mapping into a page buffer.
*/
// clang-format on
template <typename T, typename T2=void>
class RNTupleView {
   friend class RNTupleReader;
   friend class RNTupleViewCollection;

protected:
   /**
    * fFieldId has fParent always set to null; views access nested fields without looking at the parent
    */
   RField<T, T2> fField;
   Detail::RFieldValue fValue;

   RNTupleView(DescriptorId_t fieldId, Detail::RPageSource* pageSource)
     : fField(pageSource->GetDescriptor().GetFieldDescriptor(fieldId).GetFieldName()), fValue(fField.GenerateValue())
   {
      Detail::RFieldFuse::Connect(fieldId, *pageSource, fField);
      std::unordered_map<const Detail::RFieldBase *, DescriptorId_t> field2Id;
      field2Id[&fField] = fieldId;
      for (auto &f : fField) {
         auto subFieldId = pageSource->GetDescriptor().FindFieldId(f.GetName(), field2Id[f.GetParent()]);
         Detail::RFieldFuse::Connect(subFieldId, *pageSource, f);
         field2Id[&f] = subFieldId;
      }
   }

public:
   RNTupleView(const RNTupleView& other) = delete;
   RNTupleView(RNTupleView&& other) = default;
   RNTupleView& operator=(const RNTupleView& other) = delete;
   RNTupleView& operator=(RNTupleView&& other) = default;
   ~RNTupleView() { fField.DestroyValue(fValue); }

   const T& operator()(NTupleSize_t globalIndex) {
      fField.Read(globalIndex, &fValue);
      return *fValue.Get<T>();
   }

   const T& operator()(const RClusterIndex &clusterIndex) {
      fField.Read(clusterIndex, &fValue);
      return *fValue.Get<T>();
   }
};

// Template specializations in order to directly map simple types into the page pool

template <>
class RNTupleView<float> {
   friend class RNTupleReader;
   friend class RNTupleViewCollection;

protected:
   RField<float> fField;
   RNTupleView(DescriptorId_t fieldId, Detail::RPageSource* pageSource)
      : fField(pageSource->GetDescriptor().GetFieldDescriptor(fieldId).GetFieldName())
   {
      Detail::RFieldFuse::Connect(fieldId, *pageSource, fField);
   }

public:
   RNTupleView(const RNTupleView& other) = delete;
   RNTupleView(RNTupleView&& other) = default;
   RNTupleView& operator=(const RNTupleView& other) = delete;
   RNTupleView& operator=(RNTupleView&& other) = default;
   ~RNTupleView() = default;

   float operator()(NTupleSize_t globalIndex) { return *fField.Map(globalIndex); }
   float operator()(const RClusterIndex &clusterIndex) { return *fField.Map(clusterIndex); }
};

/// Used when reading from RField<float, RCustomSizedFloat>
template<>
class RNTupleView<float, RCustomSizedFloat> {
   friend class RNTupleReader;
   friend class RNTupleViewCollection;

protected:
   RField<float, RCustomSizedFloat> fField;
   RNTupleView(DescriptorId_t fieldId, Detail::RPageSource* pageSource, std::size_t nBits, std::int64_t min,
               std::int64_t max)
      : fField(pageSource->GetDescriptor().GetFieldDescriptor(fieldId).GetFieldName(), nBits, min, max)
   {
      Detail::RFieldFuse::Connect(fieldId, *pageSource, fField);
   }

public:
   RNTupleView(const RNTupleView& other) = delete;
   RNTupleView(RNTupleView&& other) = default;
   RNTupleView& operator=(const RNTupleView& other) = delete;
   RNTupleView& operator=(RNTupleView&& other) = default;
   ~RNTupleView() = default;

   double operator()(NTupleSize_t globalIndex) { return *fField.Map(globalIndex); }
   double operator()(const RClusterIndex &clusterIndex) { return *fField.Map(clusterIndex); }
};

template <>
class RNTupleView<double> {
   friend class RNTupleReader;
   friend class RNTupleViewCollection;

protected:
   RField<double> fField;
   RNTupleView(DescriptorId_t fieldId, Detail::RPageSource* pageSource)
      : fField(pageSource->GetDescriptor().GetFieldDescriptor(fieldId).GetFieldName())
   {
      Detail::RFieldFuse::Connect(fieldId, *pageSource, fField);
   }

public:
   RNTupleView(const RNTupleView& other) = delete;
   RNTupleView(RNTupleView&& other) = default;
   RNTupleView& operator=(const RNTupleView& other) = delete;
   RNTupleView& operator=(RNTupleView&& other) = default;
   ~RNTupleView() = default;

   double operator()(NTupleSize_t globalIndex) { return *fField.Map(globalIndex); }
   double operator()(const RClusterIndex &clusterIndex) { return *fField.Map(clusterIndex); }
};

/// Used when reading from RField<double, RCustomSizedFloat>
template<>
class RNTupleView<double, RCustomSizedFloat> {
   friend class RNTupleReader;
   friend class RNTupleViewCollection;

protected:
   RField<double, RCustomSizedFloat> fField;
   RNTupleView(DescriptorId_t fieldId, Detail::RPageSource* pageSource, std::size_t nBits, std::int64_t min,
               std::int64_t max)
      : fField(pageSource->GetDescriptor().GetFieldDescriptor(fieldId).GetFieldName(), nBits, min, max)
   {
      Detail::RFieldFuse::Connect(fieldId, *pageSource, fField);
   }

public:
   RNTupleView(const RNTupleView& other) = delete;
   RNTupleView(RNTupleView&& other) = default;
   RNTupleView& operator=(const RNTupleView& other) = delete;
   RNTupleView& operator=(RNTupleView&& other) = default;
   ~RNTupleView() = default;

   double operator()(NTupleSize_t globalIndex) { return *fField.Map(globalIndex); }
   double operator()(const RClusterIndex &clusterIndex) { return *fField.Map(clusterIndex); }
};

template <>
class RNTupleView<std::int32_t> {
   friend class RNTupleReader;
   friend class RNTupleViewCollection;

protected:
   RField<std::int32_t> fField;
   RNTupleView(DescriptorId_t fieldId, Detail::RPageSource* pageSource)
      : fField(pageSource->GetDescriptor().GetFieldDescriptor(fieldId).GetFieldName())
   {
      Detail::RFieldFuse::Connect(fieldId, *pageSource, fField);
   }

public:
   RNTupleView(const RNTupleView& other) = delete;
   RNTupleView(RNTupleView&& other) = default;
   RNTupleView& operator=(const RNTupleView& other) = delete;
   RNTupleView& operator=(RNTupleView&& other) = default;
   ~RNTupleView() = default;

   std::int32_t operator()(NTupleSize_t globalIndex) { return *fField.Map(globalIndex); }
   std::int32_t operator()(const RClusterIndex &clusterIndex) { return *fField.Map(clusterIndex); }
};

template <>
class RNTupleView<ClusterSize_t> {
   friend class RNTupleReader;
   friend class RNTupleViewCollection;

protected:
   RField<ClusterSize_t> fField;
   RNTupleView(DescriptorId_t fieldId, Detail::RPageSource* pageSource)
      : fField(pageSource->GetDescriptor().GetFieldDescriptor(fieldId).GetFieldName())
   {
      Detail::RFieldFuse::Connect(fieldId, *pageSource, fField);
   }

public:
   RNTupleView(const RNTupleView& other) = delete;
   RNTupleView(RNTupleView&& other) = default;
   RNTupleView& operator=(const RNTupleView& other) = delete;
   RNTupleView& operator=(RNTupleView&& other) = default;
   ~RNTupleView() = default;

   ClusterSize_t operator()(NTupleSize_t globalIndex) { return *fField.Map(globalIndex); }
   ClusterSize_t operator()(const RClusterIndex &clusterIndex) { return *fField.Map(clusterIndex); }
};

/// Used when reading from RField<std::vector<T>, RCustomSizedFloat> or RField<std::array<T, N>, RCustomSizedFloat>
template <typename T>
class RNTupleView<T, RCustomSizedFloat> {
   friend class RNTupleReader;
   friend class RNTupleViewCollection;

protected:
   RField<T, RCustomSizedFloat> fField;
   Detail::RFieldValue fValue;

   RNTupleView(DescriptorId_t fieldId, Detail::RPageSource* pageSource, std::size_t nBits, std::int64_t min, std::int64_t max)
      : fField(pageSource->GetDescriptor().GetFieldDescriptor(fieldId).GetFieldName(), nBits, min, max),
      fValue(fField.GenerateValue())
   {
      Detail::RFieldFuse::Connect(fieldId, *pageSource, fField);
      std::unordered_map<const Detail::RFieldBase *, DescriptorId_t> field2Id;
      field2Id[&fField] = fieldId;
      for (auto &f : fField) {
         auto subFieldId = pageSource->GetDescriptor().FindFieldId(f.GetName(), field2Id[f.GetParent()]);
         Detail::RFieldFuse::Connect(subFieldId, *pageSource, f);
         field2Id[&f] = subFieldId;
      }
   }

public:
   RNTupleView(const RNTupleView& other) = delete;
   RNTupleView(RNTupleView&& other) = default;
   RNTupleView& operator=(const RNTupleView& other) = delete;
   RNTupleView& operator=(RNTupleView&& other) = default;
   ~RNTupleView() = default;

   const T& operator()(NTupleSize_t globalIndex) {
      fField.Read(globalIndex, &fValue);
      return *fValue.Get<T>();
   }

   const T& operator()(const RClusterIndex &clusterIndex) {
      fField.Read(clusterIndex, &fValue);
      return *fValue.Get<T>();
   }
};


// clang-format off
/**
\class ROOT::Experimental::RNTupleViewCollection
\ingroup NTuple
\brief A view for a collection, that can itself generate new ntuple views for its nested fields.
*/
// clang-format on
class RNTupleViewCollection : public RNTupleView<ClusterSize_t> {
    friend class RNTupleReader;

private:
   Detail::RPageSource* fSource;
   DescriptorId_t fCollectionFieldId;

   RNTupleViewCollection(DescriptorId_t fieldId, Detail::RPageSource* source)
      : RNTupleView<ClusterSize_t>(fieldId, source)
      , fSource(source)
      , fCollectionFieldId(fieldId)
   {}

public:
   RNTupleViewCollection(const RNTupleViewCollection& other) = delete;
   RNTupleViewCollection(RNTupleViewCollection&& other) = default;
   RNTupleViewCollection& operator=(const RNTupleViewCollection& other) = delete;
   RNTupleViewCollection& operator=(RNTupleViewCollection&& other) = default;
   ~RNTupleViewCollection() = default;

   RNTupleClusterRange GetViewRange(NTupleSize_t globalIndex) {
      ClusterSize_t size;
      RClusterIndex collectionStart;
      fField.GetCollectionInfo(globalIndex, &collectionStart, &size);
      return RNTupleClusterRange(collectionStart.GetClusterId(), collectionStart.GetIndex(),
                                 collectionStart.GetIndex() + size);
   }
   RNTupleClusterRange GetViewRange(const RClusterIndex &clusterIndex) {
      ClusterSize_t size;
      RClusterIndex collectionStart;
      fField.GetCollectionInfo(clusterIndex, &collectionStart, &size);
      return RNTupleClusterRange(collectionStart.GetClusterId(), collectionStart.GetIndex(),
                                 collectionStart.GetIndex() + size);
   }

   template <typename T>
   RNTupleView<T> GetView(std::string_view fieldName) {
      auto fieldId = fSource->GetDescriptor().FindFieldId(fieldName, fCollectionFieldId);
      return RNTupleView<T>(fieldId, fSource);
   }
   /// ViewCreator for floating point fields where the number of bits used to store data in storage is user-defined.
   template <typename T, std::size_t nBits, std::int64_t min, std::int64_t max>
   RNTupleView<T, RCustomSizedFloat> GetView(std::string_view fieldName) {
      auto fieldId = fSource->GetDescriptor().FindFieldId(fieldName, fCollectionFieldId);
      return RNTupleView<T, RCustomSizedFloat>(fieldId, fSource, nBits, min, max);
   }
   RNTupleViewCollection GetViewCollection(std::string_view fieldName) {
      auto fieldId = fSource->GetDescriptor().FindFieldId(fieldName, fCollectionFieldId);
      return RNTupleViewCollection(fieldId, fSource);
   }

   ClusterSize_t operator()(NTupleSize_t globalIndex) {
      ClusterSize_t size;
      RClusterIndex collectionStart;
      fField.GetCollectionInfo(globalIndex, &collectionStart, &size);
      return size;
   }
   ClusterSize_t operator()(const RClusterIndex &clusterIndex) {
      ClusterSize_t size;
      RClusterIndex collectionStart;
      fField.GetCollectionInfo(clusterIndex, &collectionStart, &size);
      return size;
   }
};

} // namespace Experimental
} // namespace ROOT

#endif
