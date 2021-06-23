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

#include <TestStatistics/LikelihoodSerial.h>
#include <RooFit/MultiProcess/JobManager.h>
#include <RooFit/MultiProcess/ProcessManager.h>
#include <RooFit/MultiProcess/Queue.h>
#include <RooFit/MultiProcess/Job.h>
#include <TestStatistics/kahan_sum.h>
#include <TestStatistics/RooAbsL.h>
#include <TestStatistics/RooUnbinnedL.h>
#include <TestStatistics/RooBinnedL.h>
#include <TestStatistics/RooSubsidiaryL.h>
#include <TestStatistics/RooSumL.h>

#include "RooRealVar.h"
#include <ROOT/RMakeUnique.hxx>

namespace RooFit {
namespace TestStatistics {

LikelihoodSerial::LikelihoodSerial(std::shared_ptr<RooAbsL> likelihood, std::shared_ptr<WrapperCalculationCleanFlags> calculation_is_clean/*,
                                   RooMinimizer *minimizer*/)
   : LikelihoodWrapper(std::move(likelihood), std::move(calculation_is_clean)/*, minimizer*/)
{
   init_vars();
   // determine likelihood type
   if (dynamic_cast<RooUnbinnedL *>(likelihood_.get()) != nullptr) {
      likelihood_type = LikelihoodType::unbinned;
   } else if (dynamic_cast<RooBinnedL *>(likelihood_.get()) != nullptr) {
      likelihood_type = LikelihoodType::binned;
   } else if (dynamic_cast<RooSumL *>(likelihood_.get()) != nullptr) {
      likelihood_type = LikelihoodType::sum;
   } else if (dynamic_cast<RooSubsidiaryL *>(likelihood_.get()) != nullptr) {
      likelihood_type = LikelihoodType::subsidiary;
   } else {
      throw std::logic_error("in LikelihoodSerial constructor: _likelihood is not of a valid subclass!");
   }
   // Note to future maintainers: take care when storing the minimizer_fcn pointer. The
   // RooAbsMinimizerFcn subclasses may get cloned inside MINUIT, which means the pointer
   // should also somehow be updated in this class.
}

LikelihoodSerial *LikelihoodSerial::clone() const
{
   return new LikelihoodSerial(*this);
}

// This is a separate function (instead of just in ctor) for historical reasons.
// Its predecessor RooRealMPFE::initVars() was used from multiple ctors, but also
// from RooRealMPFE::constOptimizeTestStatistic at the end, which makes sense,
// because it might change the set of variables. We may at some point want to do
// this here as well.
void LikelihoodSerial::init_vars()
{
   // Empty current lists
   _vars.removeAll();
   _saveVars.removeAll();

   // Retrieve non-constant parameters
   auto vars = std::make_unique<RooArgSet>(
      *likelihood_->getParameters()); // TODO: make sure this is the right list of parameters, compare to original
                                     // implementation in RooRealMPFE.cxx

//   std::cout << "vars size: " << vars->getSize() << std::endl;
//   auto iter = vars->fwdIterator();
//   RooAbsArg* var;
//   int ix = 0;
//   while((var = iter.next())) {
//      printf("LikelihoodSerial::init_vars var %d = %s %p\n", ix, var->GetName(), var);
//      ++ix;
//   }

   RooArgList varList(*vars);

   // Save in lists
   _vars.add(varList);
   _saveVars.addClone(varList);
}

void LikelihoodSerial::evaluate() {
//   std::size_t N_events = likelihood_->numDataEntries();

   switch (likelihood_type) {
   case LikelihoodType::unbinned:
   case LikelihoodType::binned: {
      result = likelihood_->evaluate_partition({0, 1}, 0, 0);
      carry = likelihood_->get_carry();
      break;
   }
   case LikelihoodType::sum: {
      result = likelihood_->evaluate_partition({0, 1}, 0, likelihood_->get_N_components());
      carry = likelihood_->get_carry();
      // TODO: this normalization part below came from RooOptTestStatistic::evaluate, probably this just means you need to do the normalization on master only when doing parallel calculation. Make sure of this! In any case, it is currently not relevant, because the norm term is 1 by default and is only overridden for the RooDataWeightAverage class.
//      // Only apply global normalization if SimMaster doesn't have MP master
//      if (numSets() == 1) {
//         const Double_t norm = globalNormalization();
//         result /= norm;
//         carry /= norm;
//      }
      break;
   }
   default: {
      throw std::logic_error("in LikelihoodSerial::evaluate_task: likelihood types other than binned, unbinned and simultaneous not yet implemented!");
      break;
   }
   }

   apply_offsetting(result, carry);
}

double LikelihoodSerial::return_result() const
{
   return result;
}

} // namespace TestStatistics
} // namespace RooFit