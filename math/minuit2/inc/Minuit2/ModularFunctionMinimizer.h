// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_ModularFunctionMinimizer
#define ROOT_Minuit2_ModularFunctionMinimizer

#include "Minuit2/MnConfig.h"

#include "Minuit2/FunctionMinimizer.h"

#include <vector>

namespace ROOT {

namespace Minuit2 {

class MinimumSeedGenerator;
class MinimumBuilder;
class MinimumSeed;
class MnFcn;
class GradientCalculator;
class MnUserParameterState;
class MnUserParameters;
class MnUserCovariance;
class MnStrategy;
class FumiliFCNBase;

//_____________________________________________________________
/**
   Base common class providing the API for all the minimizer
   Various Minimize methods are provided varying on the type of
   FCN function passesd and on the objects used for the parameters
 */
class ModularFunctionMinimizer : public FunctionMinimizer {

public:
   virtual ~ModularFunctionMinimizer() {}

   // inherited interface
   virtual FunctionMinimum Minimize(const FCNBase &, const std::vector<double> &, const std::vector<double> &,
                                    unsigned int stra = 1, unsigned int maxfcn = 0, double toler = 0.1) const;

   virtual FunctionMinimum Minimize(const FCNGradientBase &, const std::vector<double> &, const std::vector<double> &,
                                    unsigned int stra = 1, unsigned int maxfcn = 0, double toler = 0.1) const;

   virtual FunctionMinimum Minimize(const FCNBase &, const std::vector<double> &, unsigned int,
                                    const std::vector<double> &, unsigned int stra = 1, unsigned int maxfcn = 0,
                                    double toler = 0.1) const;

   virtual FunctionMinimum Minimize(const FCNGradientBase &, const std::vector<double> &, unsigned int,
                                    const std::vector<double> &, unsigned int stra = 1, unsigned int maxfcn = 0,
                                    double toler = 0.1) const;

   // extension
   virtual FunctionMinimum Minimize(const FCNBase &, const MnUserParameters &, const MnStrategy &,
                                    unsigned int maxfcn = 0, double toler = 0.1) const;

   virtual FunctionMinimum Minimize(const FCNGradientBase &, const MnUserParameters &, const MnStrategy &,
                                    unsigned int maxfcn = 0, double toler = 0.1) const;

   virtual FunctionMinimum Minimize(const FCNBase &, const MnUserParameters &, const MnUserCovariance &,
                                    const MnStrategy &, unsigned int maxfcn = 0, double toler = 0.1) const;

   virtual FunctionMinimum Minimize(const FCNGradientBase &, const MnUserParameters &, const MnUserCovariance &,
                                    const MnStrategy &, unsigned int maxfcn = 0, double toler = 0.1) const;

   virtual FunctionMinimum Minimize(const FCNBase &, const MnUserParameterState &, const MnStrategy &,
                                    unsigned int maxfcn = 0, double toler = 0.1) const;

   virtual FunctionMinimum Minimize(const FCNGradientBase &, const MnUserParameterState &, const MnStrategy &,
                                    unsigned int maxfcn = 0, double toler = 0.1) const;

   // for Fumili

   //   virtual FunctionMinimum Minimize(const FumiliFCNBase&, const std::vector<double>&, const std::vector<double>&,
   //   unsigned int stra=1, unsigned int maxfcn = 0, double toler = 0.1) const;

   //   virtual FunctionMinimum Minimize(const FumiliFCNBase&, const MnUserParameters&, const MnStrategy&, unsigned int
   //   maxfcn = 0, double toler = 0.1) const;

   //   virtual FunctionMinimum Minimize(const FumiliFCNBase&, const MnUserParameters&, const MnUserCovariance&, const
   //   MnStrategy&, unsigned int maxfcn = 0, double toler = 0.1) const;

   //   virtual FunctionMinimum Minimize(const FumiliFCNBase&, const MnUserParameterState&, const MnStrategy&, unsigned
   //   int maxfcn = 0, double toler = 0.1) const;

   virtual const MinimumSeedGenerator &SeedGenerator() const = 0;
   virtual const MinimumBuilder &Builder() const = 0;
   virtual MinimumBuilder &Builder() = 0;

public:
   virtual FunctionMinimum Minimize(const MnFcn &, const GradientCalculator &, const MinimumSeed &, const MnStrategy &,
                                    unsigned int, double) const;
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_ModularFunctionMinimizer
