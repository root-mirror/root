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
   private:
      std::pair<RPageStorage::ColumnHandle_t, std::vector<RPage>> fBuf;
   public:
      void BufferPage(RPageStorage::ColumnHandle_t columnHandle, const RPage &page) {
         if (!fBuf.first) {
            fBuf.first = columnHandle;
         }
         fBuf.second.push_back(page);
      }
      const RPageStorage::ColumnHandle_t &GetHandle() const { return fBuf.first; }
      std::vector<RPage> DrainBufferedPages() {
         std::vector<RPage> drained;
         std::swap(fBuf.second, drained);
         return drained;
      }
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
