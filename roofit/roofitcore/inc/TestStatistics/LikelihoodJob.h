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
#ifndef ROOT_ROOFIT_TESTSTATISTICS_LikelihoodJob
#define ROOT_ROOFIT_TESTSTATISTICS_LikelihoodJob

#include <map>

#include "Math/MinimizerOptions.h"
#include <RooFit/MultiProcess/Job.h>
#include <RooFit/MultiProcess/types.h>
#include <TestStatistics/LikelihoodWrapper.h>

#include "RooArgList.h"

namespace RooFit {
namespace TestStatistics {

class LikelihoodJob : public MultiProcess::Job, public LikelihoodWrapper {
public:
   LikelihoodJob(std::shared_ptr<RooAbsL> _likelihood, std::shared_ptr<WrapperCalculationCleanFlags> calculation_is_clean/*, RooMinimizer *minimizer*/);
   LikelihoodJob* clone() const override;

   void init_vars();

   // TODO: implement override if necessary
//   void synchronize_with_minimizer(const ROOT::Math::MinimizerOptions & options) override;

   void evaluate() override;
   double return_result() const override;

   void update_parameters();  // helper for evaluate

   // Job overrides:
   void evaluate_task(std::size_t task) override;
   void update_real(std::size_t ix, double val, bool is_constant) override;
   void update_bool(std::size_t ix, bool value) override;
   // --- RESULT LOGISTICS ---
   void send_back_task_result_from_worker(std::size_t task) override;
   void receive_task_result_on_queue(std::size_t task, std::size_t worker_id) override;
   void send_back_results_from_queue_to_master() override;
   void clear_results() override;
   void receive_results_on_master() override;
   bool receive_task_result_on_master(const zmq::message_t & message) override;

   void enable_offsetting(bool flag) override;

private:
   double result = 0;
   double carry = 0;
   std::map<MultiProcess::Task, double> results;
   std::map<MultiProcess::Task, double> carrys;

   RooArgList _vars;      // Variables
   RooArgList _saveVars;  // Copy of variables

   LikelihoodType likelihood_type;
   std::size_t N_tasks_at_workers = 0;
};

}
}

#endif // ROOT_ROOFIT_LikelihoodJob
