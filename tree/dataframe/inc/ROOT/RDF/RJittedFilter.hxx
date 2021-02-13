// Author: Enrico Guiraud, Danilo Piparo CERN  09/2018

/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_RJITTEDFILTER
#define ROOT_RJITTEDFILTER

#include "ROOT/RDF/GraphNode.hxx"
#include "ROOT/RDF/RFilterBase.hxx"
#include "ROOT/RDF/RLoopManager.hxx"
#include "ROOT/RStringView.hxx"
#include "RtypesCore.h"

#include <memory>
#include <string>
#include <vector>

class TTreeReader;

namespace ROOT {
namespace RDF {
class RCutFlowReport;
}

namespace Detail {
namespace RDF {

namespace RDFGraphDrawing = ROOT::Internal::RDF::GraphDrawing;

/// A wrapper around a concrete RFilter, which forwards all calls to it
/// RJittedFilter is the type of the node returned by jitted Filter calls: the concrete filter can be created and set
/// at a later time, from jitted code.
class RJittedFilter final : public RFilterBase {
   std::unique_ptr<RFilterBase> fConcreteFilter = nullptr;

public:
   RJittedFilter(RLoopManager *lm, std::string_view name);
   ~RJittedFilter() { fLoopManager->Deregister(this); }

   void SetFilter(std::unique_ptr<RFilterBase> f);

   void InitSlot(TTreeReader *r, unsigned int slot) final;
   bool CheckFilters(unsigned int slot, Long64_t entry) final;
   void Report(ROOT::RDF::RCutFlowReport &) const final;
   void PartialReport(ROOT::RDF::RCutFlowReport &) const final;
   void FillReport(ROOT::RDF::RCutFlowReport &) const final;
   void IncrChildrenCount() final;
   void StopProcessing() final;
   void ResetChildrenCount() final;
   void TriggerChildrenCount() final;
   void ResetReportCount() final;
   void InitNode() final;
   void AddFilterName(std::vector<std::string> &filters) final;
   void FinaliseSlot(unsigned int slot) final;
   std::shared_ptr<RDFGraphDrawing::GraphNode> GetGraph();
};

} // ns RDF
} // ns Detail
} // ns ROOT

#endif // ROOT_RJITTEDFILTER
