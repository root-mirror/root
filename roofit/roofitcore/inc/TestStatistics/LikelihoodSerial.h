/*
 * Project: RooFit
 * Authors:
 *   PB, Patrick Bos, Netherlands eScience Center, p.bos@esciencecenter.nl
 *
 * Copyright (c) 2016-2020, Netherlands eScience Center
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 */
#ifndef ROOT_ROOFIT_LikelihoodSerial
#define ROOT_ROOFIT_LikelihoodSerial

#include <map>

#include "Math/MinimizerOptions.h"
#include <TestStatistics/LikelihoodWrapper.h>

#include "RooArgList.h"

namespace RooFit {
namespace TestStatistics {

class LikelihoodSerial : public LikelihoodWrapper {
public:
   LikelihoodSerial(std::shared_ptr<RooAbsL> likelihood, std::shared_ptr<WrapperCalculationCleanFlags> calculation_is_clean/*, RooMinimizer *minimizer*/);
   LikelihoodSerial* clone() const override;

   void init_vars();

   // TODO: implement override if necessary
//   void synchronize_with_minimizer(const ROOT::Math::MinimizerOptions & options) override;

   void evaluate() override;
   double return_result() const override;

private:
   double result = 0;
   double carry = 0;

   RooArgList _vars;      // Variables
   RooArgList _saveVars;  // Copy of variables

   LikelihoodType likelihood_type;
};

}
}

#endif // ROOT_ROOFIT_LikelihoodSerial
