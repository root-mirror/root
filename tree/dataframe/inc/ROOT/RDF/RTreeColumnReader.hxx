// Author: Enrico Guiraud CERN 09/2020

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_RDF_RTREECOLUMNREADER
#define ROOT_RDF_RTREECOLUMNREADER

#include "RColumnReaderBase.hxx"
#include <ROOT/RMakeUnique.hxx>
#include <ROOT/RVec.hxx>
#include <Rtypes.h>  // Long64_t, R__CLING_PTRCHECK
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include <memory>
#include <string>

namespace ROOT {
namespace Internal {
namespace RDF {

/// RTreeColumnReader specialization for TTree values read via TTreeReaderValues
template <typename T>
class R__CLING_PTRCHECK(off) RTreeColumnReader final : public ROOT::Detail::RDF::RColumnReaderBase {
   std::unique_ptr<TTreeReaderValue<T>> fTreeValue;

   void *GetImpl(Long64_t) final { return fTreeValue->Get(); }
public:
   /// Construct the RTreeColumnReader. Actual initialization is performed lazily by the Init method.
   RTreeColumnReader(TTreeReader &r, const std::string &colName)
      : fTreeValue(std::make_unique<TTreeReaderValue<T>>(r, colName.c_str()))
   {
   }

   /// The dtor resets the TTreeReaderValue object.
   //
   // Otherwise a race condition is present in which a TTreeReader
   // and its TTreeReader{Value,Array}s can be deleted concurrently:
   // - Thread #1) a task ends and pushes back processing slot
   // - Thread #2) a task starts and overwrites thread-local TTreeReaderValues
   // - Thread #1) first task deletes TTreeReader
   // See https://github.com/root-project/root/commit/26e8ace6e47de6794ac9ec770c3bbff9b7f2e945
   ~RTreeColumnReader() { fTreeValue.reset(); }
};

/// RTreeColumnReader specialization for TTree values read via TTreeReaderArrays.
///
/// TTreeReaderArrays are used whenever the RDF column type is RVec<T>.
template <typename T>
class R__CLING_PTRCHECK(off) RTreeColumnReader<RVec<T>> final : public ROOT::Detail::RDF::RColumnReaderBase {
   std::unique_ptr<TTreeReaderArray<T>> fTreeArray;

   /// Enumerator for the memory layout of the branch
   enum class EStorageType : char { kContiguous, kUnknown, kSparse };

   /// We return a reference to this RVec to clients, to guarantee a stable address and contiguous memory layout.
   RVec<T> fRVec;

   /// Signal whether we ever checked that the branch we are reading with a TTreeReaderArray stores array elements
   /// in contiguous memory.
   EStorageType fStorageType = EStorageType::kUnknown;

   /// Whether we already printed a warning about performing a copy of the TTreeReaderArray contents
   bool fCopyWarningPrinted = false;

   void *GetImpl(Long64_t) final
   {
      auto &readerArray = *fTreeArray;
      // We only use TTreeReaderArrays to read columns that users flagged as type `RVec`, so we need to check
      // that the branch stores the array as contiguous memory that we can actually wrap in an `RVec`.
      // Currently we need the first entry to have been loaded to perform the check
      // TODO Move check to constructor once ROOT-10823 is fixed and TTreeReaderArray itself exposes this information
      const auto arrSize = readerArray.GetSize();
      if (EStorageType::kUnknown == fStorageType && arrSize > 1) {
         // We can decide since the array is long enough
         fStorageType = EStorageType::kContiguous;
         for (auto i = 0u; i < arrSize - 1; ++i) {
            if ((char *)&readerArray[i + 1] - (char *)&readerArray[i] != sizeof(T)) {
               fStorageType = EStorageType::kSparse;
               break;
            }
         }
      }

      const auto readerArraySize = readerArray.GetSize();
      if (EStorageType::kContiguous == fStorageType ||
          (EStorageType::kUnknown == fStorageType && readerArray.GetSize() < 2)) {
         if (readerArraySize > 0) {
            // trigger loading of the contents of the TTreeReaderArray
            // the address of the first element in the reader array is not necessarily equal to
            // the address returned by the GetAddress method
            auto readerArrayAddr = &readerArray.At(0);
            RVec<T> rvec(readerArrayAddr, readerArraySize);
            std::swap(fRVec, rvec);
         } else {
            RVec<T> emptyVec{};
            std::swap(fRVec, emptyVec);
         }
      } else {
         // The storage is not contiguous or we don't know yet: we cannot but copy into the rvec
#ifndef NDEBUG
         if (!fCopyWarningPrinted) {
            Warning("RTreeColumnReader::Get",
                    "Branch %s hangs from a non-split branch. A copy is being performed in order "
                    "to properly read the content.",
                    readerArray.GetBranchName());
            fCopyWarningPrinted = true;
         }
#else
         (void)fCopyWarningPrinted;
#endif
         if (readerArraySize > 0) {
            RVec<T> rvec(readerArray.begin(), readerArray.end());
            std::swap(fRVec, rvec);
         } else {
            RVec<T> emptyVec{};
            std::swap(fRVec, emptyVec);
         }
      }
      return &fRVec;
   }

public:
   RTreeColumnReader(TTreeReader &r, const std::string &colName)
      : fTreeArray(std::make_unique<TTreeReaderArray<T>>(r, colName.c_str()))
   {
   }

   /// See the other class template specializations for an explanation.
   ~RTreeColumnReader() { fTreeArray.reset(); }
};

/// RTreeColumnReader specialization for arrays of boolean values read via TTreeReaderArrays.
///
/// TTreeReaderArray<bool> is used whenever the RDF column type is RVec<bool>.
template <>
class R__CLING_PTRCHECK(off) RTreeColumnReader<RVec<bool>> final : public ROOT::Detail::RDF::RColumnReaderBase {

   std::unique_ptr<TTreeReaderArray<bool>> fTreeArray;

   /// We return a reference to this RVec to clients, to guarantee a stable address and contiguous memory layout
   RVec<bool> fRVec;

   // We always copy the contents of TTreeReaderArray<bool> into an RVec<bool> (never take a view into the memory
   // buffer) because the underlying memory buffer might be the one of a std::vector<bool>, which is not a contiguous
   // slab of bool values.
   // Note that this also penalizes the case in which the column type is actually bool[], but the possible performance
   // gains in this edge case is probably not worth the extra complication required to differentiate the two cases.
   void *GetImpl(Long64_t) final
   {
      auto &readerArray = *fTreeArray;
      const auto readerArraySize = readerArray.GetSize();
      if (readerArraySize > 0) {
         // always perform a copy
         RVec<bool> rvec(readerArray.begin(), readerArray.end());
         std::swap(fRVec, rvec);
      } else {
         RVec<bool> emptyVec{};
         std::swap(fRVec, emptyVec);
      }
      return &fRVec;
   }

public:
   RTreeColumnReader(TTreeReader &r, const std::string &colName)
      : fTreeArray(std::make_unique<TTreeReaderArray<bool>>(r, colName.c_str()))
   {
   }

   /// See the other class template specializations for an explanation.
   ~RTreeColumnReader() { fTreeArray.reset(); }
};

} // namespace RDF
} // namespace Internal
} // namespace ROOT

#endif
