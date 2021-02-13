// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_MinimumBuilder
#define ROOT_Minuit2_MinimumBuilder

#include "Minuit2/MnTraceObject.h"
#include "Minuit2/MnPrint.h"

namespace ROOT {

namespace Minuit2 {

class FunctionMinimum;
class MnFcn;
class GradientCalculator;
class MinimumSeed;
class MinimumState;
class MnStrategy;

class MinimumBuilder {

public:
   MinimumBuilder();

   virtual ~MinimumBuilder() {}

   virtual FunctionMinimum Minimum(const MnFcn &, const GradientCalculator &, const MinimumSeed &, const MnStrategy &,
                                   unsigned int, double) const = 0;

   int StorageLevel() const { return fStorageLevel; }
   int PrintLevel() const { return fPrintLevel; }

   bool TraceIter() const { return (fTracer); }
   MnTraceObject *TraceObject() const { return (fTracer); }

   virtual void SetPrintLevel(int level) { fPrintLevel = level; }
   virtual void SetStorageLevel(int level) { fStorageLevel = level; }

   // set trace object (user manages it)
   virtual void SetTraceObject(MnTraceObject &obj) { fTracer = &obj; }

   void TraceIteration(int iter, const MinimumState &state) const
   {
      if (fTracer)
         (*fTracer)(iter, state);
   }

private:
   int fPrintLevel;
   int fStorageLevel;

   MnTraceObject *fTracer; //! tracer object (it is managed by user)
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_MinimumBuilder
