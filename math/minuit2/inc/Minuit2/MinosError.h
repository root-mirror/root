// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_MinosError
#define ROOT_Minuit2_MinosError

#include "Minuit2/MnCross.h"
#include <utility>

namespace ROOT {

namespace Minuit2 {

//____________________________________________________________________________________
/**
   Class holding the result of Minos (lower and upper values) for a specific parameter
 */

class MinosError {

public:
   MinosError() : fParameter(0), fMinParValue(0.), fUpper(MnCross()), fLower(MnCross()) {}

   MinosError(unsigned int par, double value, const MnCross &low, const MnCross &up)
      : fParameter(par), fMinParValue(value), fUpper(up), fLower(low)
   {
   }

   ~MinosError() {}

   MinosError(const MinosError &err)
      : fParameter(err.fParameter), fMinParValue(err.fMinParValue), fUpper(err.fUpper), fLower(err.fLower)
   {
   }

   MinosError &operator=(const MinosError &) = default;

   MinosError &operator()(const MinosError &err)
   {
      fParameter = err.fParameter;
      fMinParValue = err.fMinParValue;
      fUpper = err.fUpper;
      fLower = err.fLower;
      return *this;
   }

   std::pair<double, double> operator()() const { return std::pair<double, double>(Lower(), Upper()); }
   double Lower() const
   {
      if (AtLowerLimit())
         return LowerState().Parameter(Parameter()).LowerLimit() - fMinParValue;
      if (LowerValid())
         return -1. * LowerState().Error(Parameter()) * (1. + fLower.Value());
      // return Hessian Error in case is invalid
      return -LowerState().Error(Parameter());
   }
   double Upper() const
   {
      if (AtUpperLimit())
         return UpperState().Parameter(Parameter()).UpperLimit() - fMinParValue;
      if (UpperValid())
         return UpperState().Error(Parameter()) * (1. + fUpper.Value());
      // return Hessian Error in case is invalid
      return UpperState().Error(Parameter());
   }
   unsigned int Parameter() const { return fParameter; }
   const MnUserParameterState &LowerState() const { return fLower.State(); }
   const MnUserParameterState &UpperState() const { return fUpper.State(); }
   bool IsValid() const { return fLower.IsValid() && fUpper.IsValid(); }
   bool LowerValid() const { return fLower.IsValid(); }
   bool UpperValid() const { return fUpper.IsValid(); }
   bool AtLowerLimit() const { return fLower.AtLimit(); }
   bool AtUpperLimit() const { return fUpper.AtLimit(); }
   bool AtLowerMaxFcn() const { return fLower.AtMaxFcn(); }
   bool AtUpperMaxFcn() const { return fUpper.AtMaxFcn(); }
   bool LowerNewMin() const { return fLower.NewMinimum(); }
   bool UpperNewMin() const { return fUpper.NewMinimum(); }
   unsigned int NFcn() const { return fUpper.NFcn() + fLower.NFcn(); }
   // return parameter value at the minimum
   double Min() const { return fMinParValue; }

private:
   unsigned int fParameter;
   double fMinParValue;
   MnCross fUpper;
   MnCross fLower;
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_MinosError
