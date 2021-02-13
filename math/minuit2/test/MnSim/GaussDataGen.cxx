// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "GaussDataGen.h"
#include "GaussFunction.h"
#include "GaussRandomGen.h"
#include "FlatRandomGen.h"

namespace ROOT {

namespace Minuit2 {

GaussDataGen::GaussDataGen(unsigned int n)
{

   // create random generator for mean values of Gaussian [-50, 50)
   FlatRandomGen rand_mean(0., 50.);

   // create random generator for sigma values of Gaussian [1., 11.)
   FlatRandomGen rand_var(6., 5.);

   // errors of measurements (Gaussian, mean=0., sig = 0.01)
   double mvariance = 0.01 * 0.01;
   GaussRandomGen rand_mvar(0., 0.01);

   // simulate data
   fSimMean = rand_mean();
   fSimVar = rand_var();
   double sim_sig = std::sqrt(fSimVar);
   double sim_const = 1.;
   GaussFunction gauss_sim(fSimMean, sim_sig, sim_const);

   for (unsigned int i = 0; i < n; i++) {

      // x-position, from -5sigma < mean < +5sigma
      double position = fSimMean - 5. * sim_sig + double(i) * 10. * sim_sig / double(n);
      fPositions.push_back(position);

      // y-position (function Value)
      double epsilon = rand_mvar();
      fMeasurements.push_back(gauss_sim(position) + epsilon);
      fVariances.push_back(mvariance);
   }
}

} // namespace Minuit2

} // namespace ROOT
