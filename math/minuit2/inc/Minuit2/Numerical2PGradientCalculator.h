// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_Numerical2PGradientCalculator
#define ROOT_Minuit2_Numerical2PGradientCalculator

#include "Minuit2/MnConfig.h"

#include "Minuit2/GradientCalculator.h"

#include <vector>

namespace ROOT {

namespace Minuit2 {

class MnFcn;
class MnUserTransformation;
class MnMachinePrecision;
class MnStrategy;

/**
   class performing the numerical gradient calculation
 */

class Numerical2PGradientCalculator : public GradientCalculator {

public:
   Numerical2PGradientCalculator(const MnFcn &fcn, const MnUserTransformation &par, const MnStrategy &stra)
      : fFcn(fcn), fTransformation(par), fStrategy(stra)
   {
   }

   virtual ~Numerical2PGradientCalculator() {}

   virtual FunctionGradient operator()(const MinimumParameters &) const;

   virtual FunctionGradient operator()(const std::vector<double> &params) const;

   virtual FunctionGradient operator()(const MinimumParameters &, const FunctionGradient &) const;

   const MnFcn &Fcn() const { return fFcn; }
   const MnUserTransformation &Trafo() const { return fTransformation; }
   const MnMachinePrecision &Precision() const;
   const MnStrategy &Strategy() const { return fStrategy; }

   unsigned int Ncycle() const;
   double StepTolerance() const;
   double GradTolerance() const;

private:
   const MnFcn &fFcn;
   const MnUserTransformation &fTransformation;
   const MnStrategy &fStrategy;
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_Numerical2PGradientCalculator
