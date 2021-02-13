// Author: Enrico Guiraud, Danilo Piparo CERN  09/2018

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "ROOT/RDF/RBookedDefines.hxx"
#include "ROOT/RDF/RLoopManager.hxx"
#include "ROOT/RDF/RJittedAction.hxx"
// Avoid error: invalid application of ‘sizeof’ to incomplete type in RJittedAction::GetMergeableValue
#include "ROOT/RDF/RMergeableValue.hxx"
#include "TError.h"

using ROOT::Internal::RDF::RJittedAction;
using ROOT::Detail::RDF::RLoopManager;

RJittedAction::RJittedAction(RLoopManager &lm) : RActionBase(&lm, {}, ROOT::Internal::RDF::RBookedDefines{}) {}

void RJittedAction::Run(unsigned int slot, Long64_t entry)
{
   R__ASSERT(fConcreteAction != nullptr);
   fConcreteAction->Run(slot, entry);
}

void RJittedAction::Initialize()
{
   R__ASSERT(fConcreteAction != nullptr);
   fConcreteAction->Initialize();
}

void RJittedAction::InitSlot(TTreeReader *r, unsigned int slot)
{
   R__ASSERT(fConcreteAction != nullptr);
   fConcreteAction->InitSlot(r, slot);
}

void RJittedAction::TriggerChildrenCount()
{
   R__ASSERT(fConcreteAction != nullptr);
   fConcreteAction->TriggerChildrenCount();
}

void RJittedAction::FinalizeSlot(unsigned int slot)
{
   R__ASSERT(fConcreteAction != nullptr);
   fConcreteAction->FinalizeSlot(slot);
}

void RJittedAction::Finalize()
{
   R__ASSERT(fConcreteAction != nullptr);
   fConcreteAction->Finalize();
}

void *RJittedAction::PartialUpdate(unsigned int slot)
{
   R__ASSERT(fConcreteAction != nullptr);
   return fConcreteAction->PartialUpdate(slot);
}

bool RJittedAction::HasRun() const
{
   if (fConcreteAction != nullptr) {
      return fConcreteAction->HasRun();
   } else {
      // The action has not been JITted. This means that it has not run.
      return false;
   }
}

void RJittedAction::SetHasRun()
{
   R__ASSERT(fConcreteAction != nullptr);
   return fConcreteAction->SetHasRun();
}

std::shared_ptr<ROOT::Internal::RDF::GraphDrawing::GraphNode> RJittedAction::GetGraph()
{
   R__ASSERT(fConcreteAction != nullptr);
   return fConcreteAction->GetGraph();
}

/**
   Retrieve a wrapper to the result of the action that knows how to merge
   with others of the same type.
*/
std::unique_ptr<ROOT::Detail::RDF::RMergeableValueBase> RJittedAction::GetMergeableValue() const
{
   R__ASSERT(fConcreteAction != nullptr);
   return fConcreteAction->GetMergeableValue();
}
