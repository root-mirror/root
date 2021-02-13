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
#include <type_traits>
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
   class RIterator {
   private:
      RClusterIndex fIndex;
   public:
      using iterator = RIterator;
      using iterator_category = std::forward_iterator_tag;
      using value_type = RClusterIndex;
      using difference_type = RClusterIndex;
      using pointer = RClusterIndex*;
      using reference = RClusterIndex&;

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


namespace Internal {

template <class FieldT>
class IsMappable {
public:
   using RSuccess = char;
   struct RFailure { char x[2]; };

   template<class C, typename ... ArgsT>
   using MapOverloadT = decltype(std::declval<C>().Map(std::declval<ArgsT>() ...)) (C::*)(ArgsT ...);

   template <class C> static RSuccess Test(MapOverloadT<C, NTupleSize_t>);
   template <class C> static RFailure Test(...);

public:
   static constexpr bool value = sizeof(Test<FieldT>(0)) == sizeof(RSuccess);
};

} // namespace Internal


// clang-format off
/**
\class ROOT::Experimental::RNTupleView
\ingroup NTuple
\brief An RNTupleView provides read-only access to a single field of the ntuple

The view owns a field and its underlying columns in order to fill an ntuple value object with data. Data can be
accessed by index. For top-level fields, the index refers to the entry number. Fields that are part of
nested collections have global index numbers that are derived from their parent indexes.

Fields of simple types with a Map() method will use that and thus expose zero-copy access.
*/
// clang-format on
template <typename T>
class RNTupleView {
   friend class RNTupleReader;
   friend class RNTupleViewCollection;

   using FieldT = RField<T>;

private:
   /// fFieldId has fParent always set to null; views access nested fields without looking at the parent
   FieldT fField;
   /// Used as a Read() destination for fields that are not mappable
   Detail::RFieldValue fValue;

   RNTupleView(DescriptorId_t fieldId, Detail::RPageSource* pageSource)
     : fField(pageSource->GetDescriptor().GetFieldDescriptor(fieldId).GetFieldName()), fValue(fField.GenerateValue())
   {
      Detail::RFieldFuse::ConnectRecursively(fieldId, *pageSource, fField);
   }

public:
   RNTupleView(const RNTupleView& other) = delete;
   RNTupleView(RNTupleView&& other) = default;
   RNTupleView& operator=(const RNTupleView& other) = delete;
   RNTupleView& operator=(RNTupleView&& other) = default;
   ~RNTupleView() { fField.DestroyValue(fValue); }

   RNTupleGlobalRange GetFieldRange() const { return RNTupleGlobalRange(0, fField.GetNElements()); }

   template <typename C = T>
   typename std::enable_if_t<Internal::IsMappable<FieldT>::value, const C&>
   operator()(NTupleSize_t globalIndex) { return *fField.Map(globalIndex); }

   template <typename C = T>
   typename std::enable_if_t<!Internal::IsMappable<FieldT>::value, const C&>
   operator()(NTupleSize_t globalIndex) {
      fField.Read(globalIndex, &fValue);
      return *fValue.Get<T>();
   }

   template <typename C = T>
   typename std::enable_if_t<Internal::IsMappable<FieldT>::value, const C&>
   operator()(const RClusterIndex &clusterIndex) { return *fField.Map(clusterIndex); }

   template <typename C = T>
   typename std::enable_if_t<!Internal::IsMappable<FieldT>::value, const C&>
   operator()(const RClusterIndex &clusterIndex) {
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

   RNTupleClusterRange GetCollectionRange(NTupleSize_t globalIndex) {
      ClusterSize_t size;
      RClusterIndex collectionStart;
      fField.GetCollectionInfo(globalIndex, &collectionStart, &size);
      return RNTupleClusterRange(collectionStart.GetClusterId(), collectionStart.GetIndex(),
                                 collectionStart.GetIndex() + size);
   }
   RNTupleClusterRange GetCollectionRange(const RClusterIndex &clusterIndex) {
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
