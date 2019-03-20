/// \file ROOT/RColumn.hxx
/// \ingroup Forest ROOT7
/// \author Jakob Blomer <jblomer@cern.ch>
/// \date 2018-10-09
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2015, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RColumn
#define ROOT7_RColumn

#include <ROOT/RColumnElement.hxx>
#include <ROOT/RColumnModel.hxx>
#include <ROOT/RForestUtil.hxx>
#include <ROOT/RPage.hxx>
#include <ROOT/RPageStorage.hxx>

#include <TError.h>

#include <memory>
#include <vector>

namespace ROOT {
namespace Experimental {

class RColumnModel;

namespace Detail {

// clang-format off
/**
\class ROOT::Experimental::RColumn
\ingroup Forest
\brief A column is a storage-backed array of a simple, fixed-size type, from which pages can be mapped into memory.

On the primitives data layer, the RColumn and RColumnElement are the equivalents to RField and RTreeValue on the
logical data layer.
*/
// clang-format on
class RColumn {
private:
   RColumnModel fModel;
   RPageSink* fPageSink;
   RPageSource* fPageSource;
   RPageStorage::ColumnHandle_t fHandleSink;
   RPageStorage::ColumnHandle_t fHandleSource;
   /// Open page into which new elements are being written
   RPage fHeadPage;
   /// The number of elements written resp. available in the column
   ForestIndex_t fNElements;
   /// The currently mapped page for reading
   RPage fCurrentPage;
   /// The column id is used to find matching pages with content when reading
   ColumnId_t fColumnIdSource;

public:
   RColumn(const RColumnModel& model);
   RColumn(const RColumn&) = delete;
   RColumn& operator =(const RColumn&) = delete;
   ~RColumn() = default;

   void Connect(RPageStorage* pageStorage);

   void Append(const RColumnElementBase& element) {
      void* dst = fHeadPage.Reserve(1);
      if (dst == nullptr) {
         Flush();
         dst = fHeadPage.Reserve(1);
         R__ASSERT(dst != nullptr);
      }
      element.Serialize(dst, 1);
      fNElements++;
   }

   void AppendV(const RColumnElementBase& elemArray, std::size_t count) {
      void* dst = fHeadPage.Reserve(count);
      if (dst == nullptr) {
         for (unsigned i = 0; i < count; ++i) {
            Append(RColumnElementBase(elemArray, i));
         }
         return;
      }
      elemArray.Serialize(dst, count);
      fNElements += count;
   }

   void Read(const ForestIndex_t index, RColumnElementBase* element) {
      if (!fCurrentPage.Contains(index)) {
         MapPage(index);
      }
      void* src = static_cast<unsigned char *>(fCurrentPage.GetBuffer()) +
                  (index - fCurrentPage.GetRangeFirst()) * element->GetSize();
      element->Deserialize(src, 1);
   }

   void ReadV(const ForestIndex_t index, const ForestIndex_t count, RColumnElementBase* elemArray) {
      if (!fCurrentPage.Contains(index)) {
         MapPage(index);
      }
      ForestIndex_t idxInPage = index - fCurrentPage.GetRangeFirst();

      void* src = static_cast<unsigned char *>(fCurrentPage.GetBuffer()) + idxInPage * elemArray->GetSize();
      if (index + count <= fCurrentPage.GetRangeLast() + 1) {
         elemArray->Deserialize(src, count);
      } else {
         ForestIndex_t nBatch = fCurrentPage.GetRangeLast() - idxInPage;
         elemArray->Deserialize(src, nBatch);
         RColumnElementBase elemTail(*elemArray, nBatch);
         ReadV(index + nBatch, count - nBatch, &elemTail);
      }
   }

   /// Map may fall back to Read() and therefore requires a valid element
   template <typename CppT, EColumnType ColumnT>
   CppT* Map(const ForestIndex_t index, RColumnElementBase* element) {
      if (!RColumnElement<CppT, ColumnT>::kIsMappable) {
         Read(index, element);
         return static_cast<CppT*>(element->GetRawContent());
      }

      if (!fCurrentPage.Contains(index)) {
         MapPage(index);
      }
      return reinterpret_cast<CppT*>(
         static_cast<unsigned char *>(fCurrentPage.GetBuffer()) +
         (index - fCurrentPage.GetRangeFirst()) * RColumnElement<CppT, ColumnT>::kSize);
   }

   /// MapV may fail if there are less than count consecutive elements or if the type pair is not mappable
   template <typename CppT, EColumnType ColumnT>
   void* MapV(const ForestIndex_t index, const ForestIndex_t count) {
      if (!RColumnElement<CppT, ColumnT>::kIsMappable) return nullptr;
      if (!fCurrentPage.Contains(index)) {
         MapPage(index);
      }
      if (index + count > fCurrentPage.GetRangeLast() + 1) return nullptr;
      return static_cast<unsigned char *>(fCurrentPage.GetBuffer()) +
             (index - fCurrentPage.GetRangeFirst()) * kColumnElementSizes[static_cast<int>(ColumnT)];
   }

   void Flush();
   void MapPage(const ForestIndex_t index);
   ForestIndex_t GetNElements() { return fNElements; }
   const RColumnModel& GetModel() const { return fModel; }
   ColumnId_t GetColumnIdSource() const { return fColumnIdSource; }
   RPageSource* GetPageSource() const { return fPageSource; }
   RPageStorage::ColumnHandle_t GetHandleSource() const { return fHandleSource; }
};

} // namespace Detail

} // namespace Experimental
} // namespace ROOT

#endif
