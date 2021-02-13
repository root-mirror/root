// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "Minuit2/SimplexParameters.h"

namespace ROOT {

namespace Minuit2 {

void SimplexParameters::Update(double y, const MnAlgebraicVector &p)
{
   // update the SimplexParameter object with a new value y = FCN(p)
   fSimplexParameters[Jh()] = std::pair<double, MnAlgebraicVector>(y, p);
   if (y < fSimplexParameters[Jl()].first)
      fJLow = Jh();

   unsigned int jh = 0;
   for (unsigned int i = 1; i < fSimplexParameters.size(); i++) {
      if (fSimplexParameters[i].first > fSimplexParameters[jh].first)
         jh = i;
   }
   fJHigh = jh;

   return;
}

MnAlgebraicVector SimplexParameters::Dirin() const
{
   // find simplex direction (vector from big to smaller parameter points)
   MnAlgebraicVector dirin(fSimplexParameters.size() - 1);
   for (unsigned int i = 0; i < fSimplexParameters.size() - 1; i++) {
      double pbig = fSimplexParameters[0].second(i), plit = pbig;
      for (unsigned int j = 0; j < fSimplexParameters.size(); j++) {
         if (fSimplexParameters[j].second(i) < plit)
            plit = fSimplexParameters[j].second(i);
         if (fSimplexParameters[j].second(i) > pbig)
            pbig = fSimplexParameters[j].second(i);
      }
      dirin(i) = pbig - plit;
   }

   return dirin;
}

} // namespace Minuit2

} // namespace ROOT
