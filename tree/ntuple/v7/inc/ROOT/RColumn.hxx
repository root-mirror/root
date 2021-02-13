/// \file ROOT/RColumn.hxx
/// \ingroup NTuple ROOT7
/// \author Jakob Blomer <jblomer@cern.ch>
/// \date 2018-10-09
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RColumn
#define ROOT7_RColumn

#include <ROOT/RColumnElement.hxx>
#include <ROOT/RColumnModel.hxx>
#include <ROOT/RNTupleUtil.hxx>
#include <ROOT/RPage.hxx>
#include <ROOT/RPageStorage.hxx>

#include <TError.h>

#include <memory>

namespace ROOT {
namespace Experimental {
namespace Detail {

// clang-format off
/**
\class ROOT::Experimental::RColumn
\ingroup NTuple
\brief A column is a storage-backed array of a simple, fixed-size type, from which pages can be mapped into memory.

On the primitives data layer, the RColumn and RColumnElement are the equivalents to RField and RTreeValue on the
logical data layer.
*/
// clang-format on
class RColumn {
private:
   RColumnModel fModel;
   /**
    * Columns belonging to the same field are distinguished by their order.  E.g. for an std::string field, there is
    * the offset column with index 0 and the character value column with index 1.
    */
   std::uint32_t fIndex;
   RPageSink *fPageSink;
   RPageSource *fPageSource;
   RPageStorage::ColumnHandle_t fHandleSink;
   RPageStorage::ColumnHandle_t fHandleSource;
   /// Open page into which new elements are being written
   RPage fHeadPage;
   /// The number of elements written resp. available in the column
   NTupleSize_t fNElements;
   /// The currently mapped page for reading
   RPage fCurrentPage;
   /// The column id is used to find matching pages with content when reading
   ColumnId_t fColumnIdSource;
   /// Used to pack and unpack pages on writing/reading
   std::unique_ptr<RColumnElementBase> fElement;

   RColumn(const RColumnModel &model, std::uint32_t index);

public:
   template <typename CppT, EColumnType ColumnT>
   static RColumn *Create(const RColumnModel &model, std::uint32_t index) {
      R__ASSERT(model.GetType() == ColumnT);
      auto column = new RColumn(model, index);
      column->fElement = std::unique_ptr<RColumnElementBase>(new RColumnElement<CppT, ColumnT>(nullptr));
      return column;
   }

   RColumn(const RColumn&) = delete;
   RColumn &operator =(const RColumn&) = delete;
   ~RColumn();

   void Connect(DescriptorId_t fieldId, RPageStorage *pageStorage);

   void Append(const RColumnElementBase &element) {
      void *dst = fHeadPage.TryGrow(1);
      if (dst == nullptr) {
         Flush();
         dst = fHeadPage.TryGrow(1);
         R__ASSERT(dst != nullptr);
      }
      element.WriteTo(dst, 1);
      fNElements++;
   }

   void AppendV(const RColumnElementBase &elemArray, std::size_t count) {
      void *dst = fHeadPage.TryGrow(count);
      if (dst == nullptr) {
         for (unsigned i = 0; i < count; ++i) {
            Append(RColumnElementBase(elemArray, i));
         }
         return;
      }
      elemArray.WriteTo(dst, count);
      fNElements += count;
   }

   void Read(const NTupleSize_t globalIndex, RColumnElementBase *element) {
      if (!fCurrentPage.Contains(globalIndex)) {
         MapPage(globalIndex);
         R__ASSERT(fCurrentPage.Contains(globalIndex));
      }
      void *src = static_cast<unsigned char *>(fCurrentPage.GetBuffer()) +
                  (globalIndex - fCurrentPage.GetGlobalRangeFirst()) * element->GetSize();
      element->ReadFrom(src, 1);
   }

   void Read(const RClusterIndex &clusterIndex, RColumnElementBase *element) {
      if (!fCurrentPage.Contains(clusterIndex)) {
         MapPage(clusterIndex);
      }
      void *src = static_cast<unsigned char *>(fCurrentPage.GetBuffer()) +
                  (clusterIndex.GetIndex() - fCurrentPage.GetClusterRangeFirst()) * element->GetSize();
      element->ReadFrom(src, 1);
   }

   void ReadV(const NTupleSize_t globalIndex, const ClusterSize_t::ValueType count, RColumnElementBase *elemArray) {
      R__ASSERT(count > 0);
      if (!fCurrentPage.Contains(globalIndex)) {
         MapPage(globalIndex);
      }
      NTupleSize_t idxInPage = globalIndex - fCurrentPage.GetGlobalRangeFirst();

      void *src = static_cast<unsigned char *>(fCurrentPage.GetBuffer()) + idxInPage * elemArray->GetSize();
      if (globalIndex + count <= fCurrentPage.GetGlobalRangeLast() + 1) {
         elemArray->ReadFrom(src, count);
      } else {
         ClusterSize_t::ValueType nBatch = fCurrentPage.GetNElements() - idxInPage;
         elemArray->ReadFrom(src, nBatch);
         RColumnElementBase elemTail(*elemArray, nBatch);
         ReadV(globalIndex + nBatch, count - nBatch, &elemTail);
      }
   }

   void ReadV(const RClusterIndex &clusterIndex, const ClusterSize_t::ValueType count, RColumnElementBase *elemArray)
   {
      if (!fCurrentPage.Contains(clusterIndex)) {
         MapPage(clusterIndex);
      }
      NTupleSize_t idxInPage = clusterIndex.GetIndex() - fCurrentPage.GetClusterRangeFirst();

      void* src = static_cast<unsigned char *>(fCurrentPage.GetBuffer()) + idxInPage * elemArray->GetSize();
      if (clusterIndex.GetIndex() + count <= fCurrentPage.GetClusterRangeLast() + 1) {
         elemArray->ReadFrom(src, count);
      } else {
         ClusterSize_t::ValueType nBatch = fCurrentPage.GetNElements() - idxInPage;
         elemArray->ReadFrom(src, nBatch);
         RColumnElementBase elemTail(*elemArray, nBatch);
         ReadV(RClusterIndex(clusterIndex.GetClusterId(), clusterIndex.GetIndex() + nBatch), count - nBatch, &elemTail);
      }
   }

   template <typename CppT, EColumnType ColumnT>
   CppT *Map(const NTupleSize_t globalIndex) {
      if (!fCurrentPage.Contains(globalIndex)) {
         MapPage(globalIndex);
      }
      return reinterpret_cast<CppT*>(
         static_cast<unsigned char *>(fCurrentPage.GetBuffer()) +
         (globalIndex - fCurrentPage.GetGlobalRangeFirst()) * RColumnElement<CppT, ColumnT>::kSize);
   }

   template <typename CppT, EColumnType ColumnT>
   CppT *Map(const RClusterIndex &clusterIndex) {
      if (!fCurrentPage.Contains(clusterIndex)) {
         MapPage(clusterIndex);
      }
      return reinterpret_cast<CppT*>(
         static_cast<unsigned char *>(fCurrentPage.GetBuffer()) +
         (clusterIndex.GetIndex() - fCurrentPage.GetClusterRangeFirst()) * RColumnElement<CppT, ColumnT>::kSize);
   }

   NTupleSize_t GetGlobalIndex(const RClusterIndex &clusterIndex) {
      if (!fCurrentPage.Contains(clusterIndex)) {
         MapPage(clusterIndex);
      }
      return fCurrentPage.GetClusterInfo().GetIndexOffset() + clusterIndex.GetIndex();
   }

   RClusterIndex GetClusterIndex(NTupleSize_t globalIndex) {
      if (!fCurrentPage.Contains(globalIndex)) {
         MapPage(globalIndex);
      }
      return RClusterIndex(fCurrentPage.GetClusterInfo().GetId(),
                           globalIndex - fCurrentPage.GetClusterInfo().GetIndexOffset());
   }

   /// For offset columns only, look at the two adjacent values that define a collection's coordinates
   void GetCollectionInfo(const NTupleSize_t globalIndex, RClusterIndex *collectionStart, ClusterSize_t *collectionSize)
   {
      auto idxStart = (globalIndex == 0) ? 0 : *Map<ClusterSize_t, EColumnType::kIndex>(globalIndex - 1);
      auto idxEnd = *Map<ClusterSize_t, EColumnType::kIndex>(globalIndex);
      auto selfOffset = fCurrentPage.GetClusterInfo().GetIndexOffset();
      if (globalIndex == selfOffset) {
         // Passed cluster boundary
         idxStart = 0;
      }
      *collectionSize = idxEnd - idxStart;
      *collectionStart = RClusterIndex(fCurrentPage.GetClusterInfo().GetId(), idxStart);
   }

   void GetCollectionInfo(const RClusterIndex &clusterIndex,
                          RClusterIndex *collectionStart, ClusterSize_t *collectionSize)
   {
      auto index = clusterIndex.GetIndex();
      auto idxStart = (index == 0) ? 0 : *Map<ClusterSize_t, EColumnType::kIndex>(clusterIndex - 1);
      auto idxEnd = *Map<ClusterSize_t, EColumnType::kIndex>(clusterIndex);
      *collectionSize = idxEnd - idxStart;
      *collectionStart = RClusterIndex(clusterIndex.GetClusterId(), idxStart);
   }

   /// Get the currently active cluster id
   void GetSwitchInfo(NTupleSize_t globalIndex, RClusterIndex *varIndex, std::uint32_t *tag) {
      auto varSwitch = Map<RColumnSwitch, EColumnType::kSwitch>(globalIndex);
      *varIndex = RClusterIndex(fCurrentPage.GetClusterInfo().GetId(), varSwitch->GetIndex());
      *tag = varSwitch->GetTag();
   }

   void Flush();
   void MapPage(const NTupleSize_t index);
   void MapPage(const RClusterIndex &clusterIndex);
   NTupleSize_t GetNElements() const { return fNElements; }
   RColumnElementBase *GetElement() const { return fElement.get(); }
   const RColumnModel &GetModel() const { return fModel; }
   std::uint32_t GetIndex() const { return fIndex; }
   ColumnId_t GetColumnIdSource() const { return fColumnIdSource; }
   RPageSource *GetPageSource() const { return fPageSource; }
   RPageStorage::ColumnHandle_t GetHandleSource() const { return fHandleSource; }
   RPageStorage::ColumnHandle_t GetHandleSink() const { return fHandleSink; }
   RNTupleVersion GetVersion() const { return RNTupleVersion(); }
};

} // namespace Detail

} // namespace Experimental
} // namespace ROOT

#endif
