// Author: Enrico Guiraud, Danilo Piparo CERN  09/2018

/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/RDF/RJittedDefine.hxx>
#include <TError.h> // R__ASSERT

using namespace ROOT::Detail::RDF;

void RJittedDefine::InitSlot(TTreeReader *r, unsigned int slot)
{
   R__ASSERT(fConcreteDefine != nullptr);
   fConcreteDefine->InitSlot(r, slot);
}

void *RJittedDefine::GetValuePtr(unsigned int slot)
{
   R__ASSERT(fConcreteDefine != nullptr);
   return fConcreteDefine->GetValuePtr(slot);
}

const std::type_info &RJittedDefine::GetTypeId() const
{
   R__ASSERT(fConcreteDefine != nullptr);
   return fConcreteDefine->GetTypeId();
}

void RJittedDefine::Update(unsigned int slot, Long64_t entry)
{
   R__ASSERT(fConcreteDefine != nullptr);
   fConcreteDefine->Update(slot, entry);
}

void RJittedDefine::FinaliseSlot(unsigned int slot)
{
   R__ASSERT(fConcreteDefine != nullptr);
   fConcreteDefine->FinaliseSlot(slot);
}
