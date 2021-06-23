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
#include <TestStatistics/LikelihoodJob.h>
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


LikelihoodJob::LikelihoodJob(std::shared_ptr<RooAbsL> likelihood, std::shared_ptr<WrapperCalculationCleanFlags> calculation_is_clean/*, RooMinimizer *minimizer*/)
  : LikelihoodWrapper(std::move(likelihood), std::move(calculation_is_clean)/*, minimizer*/)
{
   init_vars();
   // determine likelihood type
   if (dynamic_cast<RooUnbinnedL*>(likelihood_.get()) != nullptr) {
      likelihood_type = LikelihoodType::unbinned;
   } else if (dynamic_cast<RooBinnedL*>(likelihood_.get()) != nullptr) {
      likelihood_type = LikelihoodType::binned;
   } else if (dynamic_cast<RooSumL *>(likelihood_.get()) != nullptr) {
      likelihood_type = LikelihoodType::sum;
   } else if (dynamic_cast<RooSubsidiaryL*>(likelihood_.get()) != nullptr) {
      likelihood_type = LikelihoodType::subsidiary;
   } else {
      throw std::logic_error("in LikelihoodJob constructor: likelihood is not of a valid subclass!");
   }
   // Note to future maintainers: take care when storing the minimizer_fcn pointer. The
   // RooAbsMinimizerFcn subclasses may get cloned inside MINUIT, which means the pointer
   // should also somehow be updated in this class.
}

LikelihoodJob* LikelihoodJob::clone() const {
   return new LikelihoodJob(*this);
}


// This is a separate function (instead of just in ctor) for historical reasons.
// Its predecessor RooRealMPFE::initVars() was used from multiple ctors, but also
// from RooRealMPFE::constOptimizeTestStatistic at the end, which makes sense,
// because it might change the set of variables. We may at some point want to do
// this here as well.
void LikelihoodJob::init_vars() {
   // Empty current lists
   _vars.removeAll() ;
   _saveVars.removeAll() ;

   // Retrieve non-constant parameters
   auto vars = std::make_unique<RooArgSet>(*likelihood_->getParameters());  // TODO: make sure this is the right list of parameters, compare to original implementation in RooRealMPFE.cxx
   RooArgList varList(*vars);

   // Save in lists
   _vars.add(varList);
   _saveVars.addClone(varList);
}


void LikelihoodJob::update_real(std::size_t ix, double val, bool is_constant) {
   if (get_manager()->process_manager().is_master()) {
      auto msg = RooFit::MultiProcess::M2Q::update_real;
      get_manager()->messenger().send_from_master_to_queue(msg, id, ix, val, is_constant);
   } else if (get_manager()->process_manager().is_worker()) {
      RooRealVar *rvar = (RooRealVar *) _vars.at(ix);
      rvar->setVal(static_cast<Double_t>(val));
      if (rvar->isConstant() != is_constant) {
         rvar->setConstant(static_cast<Bool_t>(is_constant));
      }
   } else {
      throw std::logic_error("LikelihoodJob::update_real only implemented on master and worker processes.");
   }
}


void LikelihoodJob::update_bool(std::size_t ix, bool value) {
   if (get_manager()->process_manager().is_master()) {
      auto msg = RooFit::MultiProcess::M2Q::update_bool;
      get_manager()->messenger().send_from_queue_to_master(msg, ix, value);
   } else if (get_manager()->process_manager().is_worker()) {
      switch(ix) {
      case 0: {
         LikelihoodWrapper::enable_offsetting(value);
         break;
      }
      default: {
         throw std::logic_error("LikelihoodJob::update_bool only supports ix = 0!");
         break;
      }
      }
   } else {
      throw std::logic_error("LikelihoodJob::update_bool only implemented on worker processes.");
   }
}


void LikelihoodJob::update_parameters() {
   if (get_manager()->process_manager().is_master()) {
      for (std::size_t ix = 0u; ix < static_cast<std::size_t>(_vars.getSize()); ++ix) {
         bool valChanged = !_vars[ix].isIdentical(_saveVars[ix], kTRUE);
         bool constChanged = (_vars[ix].isConstant() != _saveVars[ix].isConstant());

         if (valChanged || constChanged) {
            if (constChanged) {
               ((RooRealVar *) &_saveVars[ix])->setConstant(_vars[ix].isConstant());
            }
            // TODO: Check with Wouter why he uses copyCache in MPFE; makes it very difficult to extend, because copyCache is protected (so must be friend). Moved setting value to if-block below.
            //          _saveVars[ix].copyCache(&_vars[ix]);

            // send message to queue (which will relay to workers)
            RooAbsReal * rar_val = dynamic_cast<RooAbsReal *>(&_vars[ix]);
            if (rar_val) {
               Double_t val = rar_val->getVal();
               dynamic_cast<RooRealVar *>(&_saveVars[ix])->setVal(val);
               Bool_t isC = _vars[ix].isConstant();
               update_real(ix, val, isC);
            }
         }
      }
   }
}


double LikelihoodJob::return_result() const {
   return result;
}

void LikelihoodJob::evaluate() {
   if (get_manager()->process_manager().is_master()) {
      // update parameters that changed since last calculation (or creation if first time)
      update_parameters();

      // master fills queue with tasks
      for (std::size_t ix = 0; ix < get_manager()->process_manager().N_workers(); ++ix) {
         get_manager()->queue().add({id, ix});
      }
      N_tasks_at_workers = get_manager()->process_manager().N_workers();

      // wait for task results back from workers to master
      gather_worker_results();

      // put the results in vectors for calling sum_of_kahan_sums (TODO: make a map-friendly sum_of_kahan_sums)
//      std::vector<double> results_vec, carrys_vec;
//      for (auto const &item : results) {
//         results_vec.emplace_back(item.second);
//         carrys_vec.emplace_back(carrys[item.first]);
//      }
//
//      // sum task results
//      std::tie(result, carry) = sum_of_kahan_sums(results_vec, carrys_vec);
      std::tie(result, carry) = sum_of_kahan_sums(results, carrys);
      apply_offsetting(result, carry);
   }
}

// --- RESULT LOGISTICS ---

//void LikelihoodJob::send_back_task_result_from_worker(std::size_t task) {
//   get_manager()->messenger().send_from_worker_to_queue(id, task, result, carry);
//}

void LikelihoodJob::send_back_task_result_from_worker(std::size_t /*task*/) {
   get_manager()->messenger().send_from_worker_to_master(result, carry);
}

void LikelihoodJob::receive_task_result_on_queue(std::size_t task, std::size_t worker_id) {
   result = get_manager()->messenger().receive_from_worker_on_queue<double>(worker_id);
   carry = get_manager()->messenger().receive_from_worker_on_queue<double>(worker_id);
   results[task] = result;
   carrys[task] = carry;
}

void LikelihoodJob::send_back_results_from_queue_to_master() {
   get_manager()->messenger().send_from_queue_to_master(results.size());
   for (auto const &item : results) {
      get_manager()->messenger().send_from_queue_to_master(item.first, item.second, carrys[item.first]);
   }
}

bool LikelihoodJob::receive_task_result_on_master(const zmq::message_t & /*message*/) {
   std::size_t task = get_manager()->messenger().receive_from_worker_on_master<std::size_t>();
   results[task] = get_manager()->messenger().receive_from_worker_on_master<double>();
   carrys[task] = get_manager()->messenger().receive_from_worker_on_master<double>();
   --N_tasks_at_workers;
   bool job_completed = (N_tasks_at_workers == 0);
   return job_completed;
}

void LikelihoodJob::clear_results() {
   // empty results caches
   results.clear();
   carrys.clear();
}

void LikelihoodJob::receive_results_on_master() {
   std::size_t N_job_tasks = get_manager()->messenger().receive_from_queue_on_master<std::size_t>();
   for (std::size_t task_ix = 0ul; task_ix < N_job_tasks; ++task_ix) {
      std::size_t task_id = get_manager()->messenger().receive_from_queue_on_master<std::size_t>();
      results[task_id] = get_manager()->messenger().receive_from_queue_on_master<double>();
      carrys[task_id] = get_manager()->messenger().receive_from_queue_on_master<double>();
   }
}

// --- END OF RESULT LOGISTICS ---

void LikelihoodJob::evaluate_task(std::size_t task) {
   assert(get_manager()->process_manager().is_worker());

   std::size_t N_events = likelihood_->numDataEntries();

   // used to have multiple modes here, but only kept "bulk" mode; dropped interleaved, single_event and all_events from old MultiProcess::NLLVar
   std::size_t first = N_events * task / get_manager()->process_manager().N_workers();
   std::size_t last  = N_events * (task + 1) / get_manager()->process_manager().N_workers();

   switch (likelihood_type) {
   case LikelihoodType::unbinned:
   case LikelihoodType::binned: {
      result = likelihood_->evaluate_partition({static_cast<double>(first)/N_events, static_cast<double>(last)/N_events}, 0, 0);
      break;
   }
   default: {
      throw std::logic_error("in LikelihoodJob::evaluate_task: likelihood types other than binned and unbinned not yet implemented!");
      break;
   }
   }
   carry = likelihood_->get_carry();
}

void LikelihoodJob::enable_offsetting(bool flag) {
   update_bool(0, flag);
   LikelihoodWrapper::enable_offsetting(flag);
}

}
}
