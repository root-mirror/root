// Author: Stefan Wunsch, Enrico Guiraud CERN  09/2020

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "ROOT/RDFHelpers.hxx"
#include "TROOT.h"      // IsImplicitMTEnabled
#include "TError.h"     // Warning
#include "RConfigure.h" // R__USE_IMT
#ifdef R__USE_IMT
#include "ROOT/TThreadExecutor.hxx"
#endif // R__USE_IMT

#include <set>

using ROOT::RDF::RResultHandle;

void ROOT::RDF::RunGraphs(std::vector<RResultHandle> handles)
{
   if (handles.empty()) {
      Warning("RunGraphs", "Got an empty list of handles");
      return;
   }

   // Check that there are results which have not yet been run
   unsigned int nNotRun = 0;
   for (const auto &h : handles) {
      if (!h.IsReady())
         nNotRun++;
   }
   if (nNotRun < handles.size()) {
      Warning("RunGraphs", "Got %lu handles from which %lu link to results which are already ready.", handles.size(),
              handles.size() - nNotRun);
      return;
   }
   if (nNotRun == 0)
      return;

   // Find the unique event loops
   auto sameGraph = [](const RResultHandle &a, const RResultHandle &b) { return a.fLoopManager < b.fLoopManager; };
   std::set<RResultHandle, decltype(sameGraph)> s(handles.begin(), handles.end(), sameGraph);
   std::vector<RResultHandle> uniqueLoops(s.begin(), s.end());

   // Trigger the unique event loops
   auto run = [](RResultHandle &h) { h.fLoopManager->Run(); };
#ifdef R__USE_IMT
   if (ROOT::IsImplicitMTEnabled()) {
      ROOT::TThreadExecutor{}.Foreach(run, uniqueLoops);
      return;
   }
#endif // R__USE_IMT
   for (auto &h : uniqueLoops)
      run(h);
}
