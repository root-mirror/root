// Author: Enrico Guiraud, Danilo Piparo CERN  09/2018

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_RJITTEDACTION
#define ROOT_RJITTEDACTION

#include "ROOT/RDF/RActionBase.hxx"
#include "ROOT/RDF/RLoopManager.hxx"
#include "RtypesCore.h"

#include <memory>

class TTreeReader;

namespace ROOT {
namespace Detail {
namespace RDF {
class RMergeableValueBase;
} // namespace RDF
} // namespace Detail
} // namespace ROOT

namespace ROOT {
namespace Internal {
namespace RDF {

// fwd decl
namespace GraphDrawing {
class GraphNode;
} // ns GraphDrawing

class RJittedAction : public RActionBase {
private:
   std::unique_ptr<RActionBase> fConcreteAction;

public:
   RJittedAction(RLoopManager &lm);
   ~RJittedAction() { fLoopManager->Deregister(this); }

   void SetAction(std::unique_ptr<RActionBase> a) { fConcreteAction = std::move(a); }

   void Run(unsigned int slot, Long64_t entry) final;
   void Initialize() final;
   void InitSlot(TTreeReader *r, unsigned int slot) final;
   void TriggerChildrenCount() final;
   void FinalizeSlot(unsigned int) final;
   void Finalize() final;
   void *PartialUpdate(unsigned int slot) final;
   bool HasRun() const final;
   void SetHasRun() final;

   std::shared_ptr<GraphDrawing::GraphNode> GetGraph();

   // Helper for RMergeableValue
   std::unique_ptr<ROOT::Detail::RDF::RMergeableValueBase> GetMergeableValue() const final;
};

} // ns RDF
} // ns Internal
} // ns ROOT

#endif // ROOT_RJITTEDACTION
