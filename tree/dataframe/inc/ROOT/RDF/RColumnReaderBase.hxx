// Author: Enrico Guiraud CERN 09/2020

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_INTERNAL_RDF_RCOLUMNREADERBASE
#define ROOT_INTERNAL_RDF_RCOLUMNREADERBASE

#include <Rtypes.h>

namespace ROOT {
namespace Detail {
namespace RDF {

/**
\class ROOT::Internal::RDF::RColumnReaderBase
\ingroup dataframe
\brief Pure virtual base class for all column reader types

This pure virtual class provides a common base class for the different column reader types, e.g. RTreeColumnReader and
RDSColumnReader.
**/
class R__CLING_PTRCHECK(off) RColumnReaderBase {
public:
   virtual ~RColumnReaderBase() = default;

   /// Return the column value for the given entry. Called at most once per entry.
   /// \tparam T The column type
   /// \param entry The entry number
   template <typename T>
   T &Get(Long64_t entry)
   {
      return *static_cast<T *>(GetImpl(entry));
   }

private:
   virtual void *GetImpl(Long64_t entry) = 0;
};

} // namespace RDF
} // namespace Detail
} // namespace ROOT

#endif
