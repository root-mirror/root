// @(#)root/tree:$Id$
// Author: Brian Bockelman, 2017

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TBranchIMTHelper
#define ROOT_TBranchIMTHelper

#include "RtypesCore.h"

#ifdef R__USE_IMT
#include "ROOT/TTaskGroup.hxx"
#endif

/// A helper class for managing IMT work during TTree:Fill operations.
///
namespace ROOT {
namespace Internal {

class TBranchIMTHelper {

#ifdef R__USE_IMT
using TaskGroup_t = ROOT::Experimental::TTaskGroup;
#endif

public:
   template<typename FN> void Run(const FN &lambda) {
#ifdef R__USE_IMT
      if (!fGroup) { fGroup.reset(new TaskGroup_t()); }
      fGroup->Run( [=]() {
         auto nbytes = lambda();
         if (nbytes >= 0) {
            fBytes += nbytes;
         } else {
            ++fNerrors;
         }
      });
#else
      (void)lambda;
#endif
   }

   void Wait() {
#ifdef R__USE_IMT
      if (fGroup) fGroup->Wait();
#endif
   }

   Long64_t GetNbytes() { return fBytes; }
   Long64_t GetNerrors() {  return fNerrors; }

private:
   std::atomic<Long64_t> fBytes{0};   // Total number of bytes written by this helper.
   std::atomic<Int_t>    fNerrors{0}; // Total error count of all tasks done by this helper.
#ifdef R__USE_IMT
   std::unique_ptr<TaskGroup_t> fGroup;
#endif
};

} // Internal
} // ROOT

#endif
