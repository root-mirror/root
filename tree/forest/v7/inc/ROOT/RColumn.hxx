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
#include <ROOT/RForestUtil.hxx>
#include <ROOT/RPage.hxx>

#include <TError.h>

#include <memory>
#include <vector>

namespace ROOT {
namespace Experimental {

class RColumnModel;

namespace Detail {

class RPageStorage;
class RPageSink;
class RPageSource;

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
   friend class RPageSink; // to let it set fHeadPage

private:
   RColumnModel fModel;
   RPageSink* fPageSink;
   RPageSource* fPageSource;
   /// Open page into which new elements are being written
   RPage fHeadPage;
   TreeIndex_t fRangeEnd;

public:
   RColumn(const RColumnModel &model, RPageStorage *pageStorage);
   RColumn(const RColumn&) = delete;
   RColumn& operator =(const RColumn&) = delete;
   ~RColumn() = default;

   void Append(const RColumnElementBase& element) {
      void* dst = fHeadPage.Reserve(1);
      if (dst == nullptr) {
         Flush();
         dst = fHeadPage.Reserve(1);
         R__ASSERT(dst != nullptr);
      }
      element.Serialize(dst);
      fRangeEnd++;
   }
   void Flush();

   void Read(const TreeIndex_t /*index*/, RColumnElementBase* /*element*/) {/*...*/}
   void Map(const TreeIndex_t /*index*/, void ** /*dst*/) {/*...*/}

   // Returns the number of mapped values
   TreeIndex_t MapV(const TreeIndex_t /*index*/, const TreeIndex_t /*count*/, void ** /*dst*/) {return 0;/*...*/}
   void ReadV(const TreeIndex_t /*index*/, const TreeIndex_t /*count*/, void * /*dst*/) {/*...*/}

   TreeIndex_t GetNElements();
   const RColumnModel& GetModel() const { return fModel; }
};

} // namespace Detail

} // namespace Experimental
} // namespace ROOT

#endif
