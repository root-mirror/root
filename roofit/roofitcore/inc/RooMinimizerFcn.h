/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id$
 * Authors:                                                                  *
 *   AL, Alfio Lazzaro,   INFN Milan,        alfio.lazzaro@mi.infn.it        *
 *                                                                           *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef __ROOFIT_NOROOMINIMIZER

#ifndef ROO_MINIMIZER_FCN
#define ROO_MINIMIZER_FCN

#include "Math/IFunction.h"
#include "Fit/ParameterSettings.h"
#include "Fit/FitResult.h"

#include "RooAbsReal.h"
#include "RooArgList.h"

#include <fstream>
#include <vector>

#include <RooAbsMinimizerFcn.h>

template<typename T> class TMatrixTSym;
using TMatrixDSym = TMatrixTSym<double>;

// forward declaration
class RooMinimizer;

class RooMinimizerFcn : public RooAbsMinimizerFcn, public ROOT::Math::IBaseFunctionMultiDim {

public:
   RooMinimizerFcn(RooAbsReal *funct, RooMinimizer *context, bool verbose = false);
   RooMinimizerFcn(const RooMinimizerFcn &other);
   virtual ~RooMinimizerFcn();

   ROOT::Math::IBaseFunctionMultiDim *Clone() const override;
   unsigned int NDim() const override { return get_nDim(); }

   std::string getFunctionName() const override;
   std::string getFunctionTitle() const override;

   void setOptimizeConst(Int_t flag) override;

private:
   double DoEval(const double *x) const override;
   void optimizeConstantTerms(bool constStatChange, bool constValChange) override;

   RooAbsReal *_funct;
};

#endif
#endif
