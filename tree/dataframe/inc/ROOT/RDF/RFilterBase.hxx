// Author: Enrico Guiraud, Danilo Piparo CERN  09/2018

/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_RFILTERBASE
#define ROOT_RFILTERBASE

#include "ROOT/RDF/RBookedDefines.hxx"
#include "ROOT/RDF/RNodeBase.hxx"
#include "RtypesCore.h"
#include "TError.h" // R_ASSERT

#include <string>
#include <vector>

class TTreeReader;

namespace ROOT {

namespace RDF {
class RCutFlowReport;
} // ns RDF

namespace Detail {
namespace RDF {
namespace RDFInternal = ROOT::Internal::RDF;

class RLoopManager;

class RFilterBase : public RNodeBase {
protected:
   std::vector<Long64_t> fLastCheckedEntry;
   std::vector<int> fLastResult = {true}; // std::vector<bool> cannot be used in a MT context safely
   std::vector<ULong64_t> fAccepted = {0};
   std::vector<ULong64_t> fRejected = {0};
   const std::string fName;
   const unsigned int fNSlots; ///< Number of thread slots used by this node, inherited from parent node.

   RDFInternal::RBookedDefines fDefines;

public:
   RFilterBase(RLoopManager *df, std::string_view name, const unsigned int nSlots,
               const RDFInternal::RBookedDefines &defines);
   RFilterBase &operator=(const RFilterBase &) = delete;

   virtual ~RFilterBase();

   virtual void InitSlot(TTreeReader *r, unsigned int slot) = 0;
   bool HasName() const;
   std::string GetName() const;
   virtual void FillReport(ROOT::RDF::RCutFlowReport &) const;
   virtual void TriggerChildrenCount() = 0;
   virtual void ResetReportCount()
   {
      R__ASSERT(!fName.empty()); // this method is to only be called on named filters
      // fAccepted and fRejected could be different than 0 if this is not the first event-loop run using this filter
      std::fill(fAccepted.begin(), fAccepted.end(), 0);
      std::fill(fRejected.begin(), fRejected.end(), 0);
   }
   /// Clean-up operations to be performed at the end of a task.
   virtual void FinaliseSlot(unsigned int slot) = 0;
   virtual void InitNode();
   virtual void AddFilterName(std::vector<std::string> &filters) = 0;
};

} // ns RDF
} // ns Detail
} // ns ROOT

#endif // ROOT_RFILTERBASE
