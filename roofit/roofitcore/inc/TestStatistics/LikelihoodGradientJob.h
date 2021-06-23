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
#ifndef ROOT_ROOFIT_TESTSTATISTICS_LikelihoodGradientJob
#define ROOT_ROOFIT_TESTSTATISTICS_LikelihoodGradientJob

#include <vector>
#include "Math/MinimizerOptions.h"
#include <RooFit/MultiProcess/Job.h>
#include <TestStatistics/LikelihoodGradientWrapper.h>
#include <NumericalDerivatorMinuit2.h>
#include "Minuit2/MnMatrix.h"

namespace RooFit {
namespace TestStatistics {

class LikelihoodGradientJob : public MultiProcess::Job, public LikelihoodGradientWrapper {
public:
   LikelihoodGradientJob(std::shared_ptr<RooAbsL> likelihood, std::shared_ptr<WrapperCalculationCleanFlags> calculation_is_clean, std::size_t N_dim, RooMinimizer* minimizer);
   LikelihoodGradientJob* clone() const override;
   LikelihoodGradientJob(const LikelihoodGradientJob& other);

   void fill_gradient(double *grad) override;
   void fill_second_derivative(double *g2) override;
   void fill_step_size(double *gstep) override;

   void update_state() override;

   // ----- BEGIN PASTE UIT RooGradientFunction.h -----
   // ----- BEGIN PASTE UIT RooGradientFunction.h -----
   // ----- BEGIN PASTE UIT RooGradientFunction.h -----
   // ----- BEGIN PASTE UIT RooGradientFunction.h -----
   // ----- BEGIN PASTE UIT RooGradientFunction.h -----

public:
   enum class GradientCalculatorMode {
      ExactlyMinuit2, AlmostMinuit2
   };

private:
   // TODO: are mutables here still necessary?
   // mutables below are because ROOT::Math::IMultiGradFunction::DoDerivative is const

   mutable std::vector<RooFit::MinuitDerivatorElement> _grad;
   mutable RooFit::NumericalDerivatorMinuit2 _gradf;

   void run_derivator(unsigned int i_component) const;

   // ----- END PASTE UIT RooGradientFunction.h -----
   // ----- END PASTE UIT RooGradientFunction.h -----
   // ----- END PASTE UIT RooGradientFunction.h -----
   // ----- END PASTE UIT RooGradientFunction.h -----
   // ----- END PASTE UIT RooGradientFunction.h -----

   void synchronize_parameter_settings(ROOT::Math::IMultiGenFunction* function, const std::vector<ROOT::Fit::ParameterSettings> &parameter_settings) override;

   void synchronize_with_minimizer(const ROOT::Math::MinimizerOptions & options) override;
   void set_strategy(int istrat);
   void set_step_tolerance(double step_tolerance) const;
   void set_grad_tolerance(double grad_tolerance) const;
   void set_ncycles(unsigned int ncycles) const;
   void set_error_level(double error_level) const;

   void update_minuit_internal_parameter_values(const std::vector<double>& minuit_internal_x) override;

   bool uses_minuit_internal_values() override;

   // Job overrides:
   void evaluate_task(std::size_t task) override;
   void update_real(std::size_t ix, double val, bool is_constant) override;
   void update_bool(std::size_t ix, bool value) override;
   void send_back_task_result_from_worker(std::size_t task) override;
   void receive_task_result_on_queue(std::size_t task, std::size_t worker_id) override;
   void send_back_results_from_queue_to_master() override;
   void clear_results() override;
   void receive_results_on_master() override;

   struct task_result_t {
      std::size_t job_id;
      std::size_t task_id;
      MinuitDerivatorElement grad;
   };
   bool receive_task_result_on_master(const zmq::message_t & message) override;

   void update_workers_state();
   void calculate_all();

   // members
   std::size_t N_tasks = 0;
   std::size_t N_tasks_at_workers = 0;
   std::vector<std::size_t> completed_task_ids;
   std::vector<double> minuit_internal_x_;
};

}
}


#endif // ROOT_ROOFIT_LikelihoodGradientJob
