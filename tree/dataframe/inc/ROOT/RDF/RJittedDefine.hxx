// Author: Enrico Guiraud, Danilo Piparo CERN  09/2018

/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_RJITTEDCUSTOMCOLUMN
#define ROOT_RJITTEDCUSTOMCOLUMN

#include "ROOT/RDF/RDefineBase.hxx"
#include "ROOT/RStringView.hxx"
#include "RtypesCore.h"

#include <memory>
#include <type_traits>

class TTreeReader;

namespace ROOT {
namespace Detail {
namespace RDF {

/// A wrapper around a concrete RDefine, which forwards all calls to it
/// RJittedDefine is a placeholder that is put in the collection of custom columns in place of a RDefine
/// that will be just-in-time compiled. Jitted code will assign the concrete RDefine to this RJittedDefine
/// before the event-loop starts.
class RJittedDefine : public RDefineBase {
   std::unique_ptr<RDefineBase> fConcreteDefine = nullptr;

public:
   RJittedDefine(std::string_view name, std::string_view type, unsigned int nSlots,
                       const std::map<std::string, std::vector<void *>> &DSValuePtrs)
      : RDefineBase(name, type, nSlots, RDFInternal::RBookedDefines(), DSValuePtrs, nullptr)
   {
   }

   void SetDefine(std::unique_ptr<RDefineBase> c) { fConcreteDefine = std::move(c); }

   void InitSlot(TTreeReader *r, unsigned int slot) final;
   void *GetValuePtr(unsigned int slot) final;
   const std::type_info &GetTypeId() const final;
   void Update(unsigned int slot, Long64_t entry) final;
   void FinaliseSlot(unsigned int slot) final;
};

} // ns RDF
} // ns Detail
} // ns ROOT

#endif // ROOT_RJITTEDCUSTOMCOLUMN
