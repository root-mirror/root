// Author: Enrico Guiraud CERN 09/2020

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_RDF_RDSCOLUMNREADER
#define ROOT_RDF_RDSCOLUMNREADER

#include "RColumnReaderBase.hxx"
#include <Rtypes.h>  // Long64_t, R__CLING_PTRCHECK

namespace ROOT {
namespace Internal {
namespace RDF {

/// Column reader type that deals with values read from RDataSources.
template <typename T>
class R__CLING_PTRCHECK(off) RDSColumnReader final : public ROOT::Detail::RDF::RColumnReaderBase {
   T **fDSValuePtr = nullptr;

   void *GetImpl(Long64_t) final { return *fDSValuePtr; }

public:
   RDSColumnReader(void *DSValuePtr) : fDSValuePtr(static_cast<T **>(DSValuePtr)) {}
};

} // namespace RDF
} // namespace Internal
} // namespace ROOT

#endif
