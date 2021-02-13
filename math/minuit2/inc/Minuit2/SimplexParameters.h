// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_SimplexParameters
#define ROOT_Minuit2_SimplexParameters

#include <cassert>

#include "Minuit2/MnMatrix.h"

#include <vector>
#include <utility>

namespace ROOT {

namespace Minuit2 {

/**
   class describing the simplex set of points (f(x), x )  which evolve during the minimization
   iteration process.
 */

class SimplexParameters {

public:
   SimplexParameters(const std::vector<std::pair<double, MnAlgebraicVector>> &simpl, unsigned int jh, unsigned int jl)
      : fSimplexParameters(simpl), fJHigh(jh), fJLow(jl)
   {
   }

   ~SimplexParameters() {}

   void Update(double, const MnAlgebraicVector &);

   const std::vector<std::pair<double, MnAlgebraicVector>> &Simplex() const { return fSimplexParameters; }

   const std::pair<double, MnAlgebraicVector> &operator()(unsigned int i) const
   {
      assert(i < fSimplexParameters.size());
      return fSimplexParameters[i];
   }

   unsigned int Jh() const { return fJHigh; }
   unsigned int Jl() const { return fJLow; }
   double Edm() const { return fSimplexParameters[Jh()].first - fSimplexParameters[Jl()].first; }
   MnAlgebraicVector Dirin() const;

private:
   std::vector<std::pair<double, MnAlgebraicVector>> fSimplexParameters;
   unsigned int fJHigh;
   unsigned int fJLow;
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_SimplexParameters
