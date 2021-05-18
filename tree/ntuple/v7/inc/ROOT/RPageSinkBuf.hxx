/// \file ROOT/RPageSinkBuf.hxx
/// \ingroup NTuple ROOT7
/// \author Jakob Blomer <jblomer@cern.ch>
/// \author Max Orok <maxwellorok@gmail.com>
/// \date 2021-03-17
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2021, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RPageSinkBuf
#define ROOT7_RPageSinkBuf

#include <ROOT/RPageStorage.hxx>

#include <iterator>
#include <list>
#include <memory>

namespace ROOT {
namespace Experimental {
namespace Detail {

// clang-format off
/**
\class ROOT::Experimental::Detail::RPageSinkBuf
\ingroup NTuple
\brief Wrapper sink that coalesces cluster column page writes
*/
// clang-format on
class RPageSinkBuf : public RPageSink {
private:
   /// A buffered column. The column is not responsible for RPage memory management (i.e.
   /// ReservePage/ReleasePage), which is handled by the enclosing RPageSinkBuf.
   class RColumnBuf {
   public:
      struct RPageZipItem {
         RPage fPage;
         // Compression scratch buffer for fSealedPage.
         std::unique_ptr<unsigned char[]> fBuf;
         RPageStorage::RSealedPage fSealedPage;
         explicit RPageZipItem(RPage page)
            : fPage(page), fBuf(nullptr) {}
         bool IsSealed() {
            return fSealedPage.fBuffer != nullptr;
         }
         void AllocateSealedPageBuf() {
            fBuf = std::make_unique<unsigned char[]>(fPage.GetSize());
         }
      };
   public:
      RColumnBuf() = default;
      RColumnBuf(const RColumnBuf&) = delete;
      RColumnBuf& operator=(const RColumnBuf&) = delete;
      RColumnBuf(RColumnBuf&&) = default;
      RColumnBuf& operator=(RColumnBuf&&) = default;
      ~RColumnBuf() = default;
      /// Returns an iterator to the newly buffered page. The iterator remains
      /// valid until the return value of DrainBufferedPages() is destroyed.
      std::list<RPageZipItem>::iterator BufferPage(
         RPageStorage::ColumnHandle_t columnHandle, const RPage &page)
      {
         if (!fCol) {
            fCol = columnHandle;
         }
         fBufferedPages.push_back(RPageZipItem(page));
         return std::prev(fBufferedPages.end());
      }
      const RPageStorage::ColumnHandle_t &GetHandle() const { return fCol; }
      // When the return value of DrainBufferedPages() is destroyed, all iterators
      // returned by GetBuffer are invalidated.
      std::list<RPageZipItem> DrainBufferedPages() {
         std::list<RPageZipItem> drained;
         std::swap(fBufferedPages, drained);
         return drained;
      }
   private:
      RPageStorage::ColumnHandle_t fCol;
      // Using a linked list guarantees that references to list elements are
      // never invalidated by appends in BufferPage.
      std::list<RPageZipItem> fBufferedPages;
   };

private:
   /// The inner sink, responsible for actually performing I/O.
   std::unique_ptr<RPageSink> fInnerSink;
   /// The buffered page sink maintains a copy of the RNTupleModel for the inner sink.
   /// For the unbuffered case, the RNTupleModel is instead managed by a RNTupleWriter.
   std::unique_ptr<RNTupleModel> fInnerModel;
   /// Vector of buffered column pages. Indexed by column id.
   std::vector<RColumnBuf> fBufferedColumns;

protected:
   void CreateImpl(const RNTupleModel &model) final;
   RClusterDescriptor::RLocator CommitPageImpl(ColumnHandle_t columnHandle, const RPage &page) final;
   RClusterDescriptor::RLocator CommitSealedPageImpl(DescriptorId_t columnId, const RSealedPage &sealedPage) final;
   RClusterDescriptor::RLocator CommitClusterImpl(NTupleSize_t nEntries) final;
   void CommitDatasetImpl() final;

public:
   explicit RPageSinkBuf(std::unique_ptr<RPageSink> inner);
   RPageSinkBuf(const RPageSink&) = delete;
   RPageSinkBuf& operator=(const RPageSinkBuf&) = delete;
   RPageSinkBuf(RPageSinkBuf&&) = default;
   RPageSinkBuf& operator=(RPageSinkBuf&&) = default;
   virtual ~RPageSinkBuf() = default;

   RPage ReservePage(ColumnHandle_t columnHandle, std::size_t nElements = 0) final;
   void ReleasePage(RPage &page) final;

   RNTupleMetrics &GetMetrics() final { return fInnerSink->GetMetrics(); }
};

} // namespace Detail
} // namespace Experimental
} // namespace ROOT

#endif
