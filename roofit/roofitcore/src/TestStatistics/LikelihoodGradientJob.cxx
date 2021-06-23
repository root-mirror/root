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
#include <TestStatistics/LikelihoodGradientJob.h>
#include <Minuit2/MnStrategy.h>

#include <RooTimer.h>
#include "RooMsgService.h"
#include "RooFit/MultiProcess/JobManager.h"
#include "RooFit/MultiProcess/Messenger.h"
#include "RooFit/MultiProcess/Queue.h"
#include "RooMinimizer.h"

namespace RooFit {
namespace TestStatistics {

LikelihoodGradientJob::LikelihoodGradientJob(std::shared_ptr<RooAbsL> likelihood, std::shared_ptr<WrapperCalculationCleanFlags> calculation_is_clean, std::size_t N_dim,
                                             RooMinimizer *minimizer)
   : LikelihoodGradientWrapper(std::move(likelihood), std::move(calculation_is_clean), N_dim, minimizer), _grad(N_dim)
{
   // Note to future maintainers: take care when storing the minimizer_fcn pointer. The
   // RooAbsMinimizerFcn subclasses may get cloned inside MINUIT, which means the pointer
   // should also somehow be updated in this class.
//   N_tasks = minimizer_fcn->get_nDim();
   N_tasks = N_dim;
   completed_task_ids.reserve(N_tasks);
   minuit_internal_x_.reserve(N_dim);
   // TODO: make sure that the full gradients are sent back so that the
   // derivator will depart from correct state next step everywhere!
}

LikelihoodGradientJob::LikelihoodGradientJob(const LikelihoodGradientJob &other)
   : LikelihoodGradientWrapper(other), _grad(other._grad), _gradf(other._gradf), N_tasks(other.N_tasks),
     completed_task_ids(other.completed_task_ids), minuit_internal_x_(other.minuit_internal_x_)
{
}

LikelihoodGradientJob *LikelihoodGradientJob::clone() const
{
   return new LikelihoodGradientJob(*this);
}

void LikelihoodGradientJob::synchronize_parameter_settings(
   ROOT::Math::IMultiGenFunction *function, const std::vector<ROOT::Fit::ParameterSettings> &parameter_settings)
{
   _gradf.SetInitialGradient(function, parameter_settings, _grad);
}

void LikelihoodGradientJob::synchronize_with_minimizer(const ROOT::Math::MinimizerOptions &options)
{
   set_strategy(options.Strategy());
   set_error_level(options.ErrorDef());
}

void LikelihoodGradientJob::set_strategy(int istrat)
{
   assert(istrat >= 0);
   ROOT::Minuit2::MnStrategy strategy(static_cast<unsigned int>(istrat));

   set_step_tolerance(strategy.GradientStepTolerance());
   set_grad_tolerance(strategy.GradientTolerance());
   set_ncycles(strategy.GradientNCycles());
}

void LikelihoodGradientJob::set_step_tolerance(double step_tolerance) const
{
   _gradf.set_step_tolerance(step_tolerance);
}

void LikelihoodGradientJob::set_grad_tolerance(double grad_tolerance) const
{
   _gradf.set_grad_tolerance(grad_tolerance);
}

void LikelihoodGradientJob::set_ncycles(unsigned int ncycles) const
{
   _gradf.set_ncycles(ncycles);
}

void LikelihoodGradientJob::set_error_level(double error_level) const
{
   _gradf.set_error_level(error_level);
}

///////////////////////////////////////////////////////////////////////////////
/// Job overrides:

void LikelihoodGradientJob::evaluate_task(std::size_t task)
{
//   RooWallTimer timer;
//   RooCPUTimer ctimer;
   run_derivator(task);
//   ctimer.stop();
//   timer.stop();
//   oocxcoutD((TObject *)nullptr, Benchmarking1)
//      << "worker_id: " << get_manager()->process_manager().worker_id() << ", task: " << task
//      << ", partial derivative time: " << timer.timing_s() << "s -- cputime: " << ctimer.timing_s() << "s" << std::endl;
}

void LikelihoodGradientJob::update_real(std::size_t ix, double val, bool /*is_constant*/)
{
   if (get_manager()->process_manager().is_worker()) {
      // ix is defined in "flat" FunctionGradient space ix_dim * size + ix_component
//      std::cout << "on worker, ix = " << ix << ", ix / _minimizer->getNPar() = " << ix / _minimizer->getNPar() << ", ix % _minimizer->getNPar() = " << ix % _minimizer->getNPar() << std::endl;
      switch (ix / _minimizer->getNPar()) {
      case 0: {
         _grad[ix % _minimizer->getNPar()].derivative = val;
         break;
      }
      case 1: {
         _grad[ix % _minimizer->getNPar()].second_derivative = val;
         break;
      }
      case 2: {
         _grad[ix % _minimizer->getNPar()].step_size = val;
         break;
      }
      case 3: {
         std::size_t ix_component = ix % _minimizer->getNPar();
         minuit_internal_x_[ix_component] = val;
//         _minimizer->set_function_parameter_value(ix_component, val);  // if we want to update this, we should send over external values!
         break;
      }
      default:
         std::stringstream ss;
         ss << "ix = " << ix << " out of range in LikelihoodGradientJob::update_real! NPar = " << _minimizer->getNPar();
         throw std::runtime_error(ss.str());
      }
   }
}

void LikelihoodGradientJob::update_bool(std::size_t /*ix*/, bool /*value*/) {}

// SYNCHRONIZATION FROM WORKERS TO MASTER

//void LikelihoodGradientJob::send_back_task_result_from_worker(std::size_t task)
//{
////   std::cout << "worker " << get_manager()->process_manager().worker_id() << " sends task result id=" << id
////             << " task=" << task << " grad=" << _grad[task].derivative << " g2=" << _grad[task].second_derivative
////             << " gstep=" << _grad[task].step_size << std::endl;
//   get_manager()->messenger().send_from_worker_to_queue(_grad[task].derivative, _grad[task].second_derivative, _grad[task].step_size);
//}

void LikelihoodGradientJob::send_back_task_result_from_worker(std::size_t task)
{
   task_result_t task_result {id, task, _grad[task]};
   zmq::message_t message(sizeof(task_result_t));
   memcpy(message.data(), &task_result, sizeof(task_result_t));
   get_manager()->messenger().send_from_worker_to_master(std::move(message));
}

void LikelihoodGradientJob::receive_task_result_on_queue(std::size_t task, std::size_t worker_id)
{
   completed_task_ids.push_back(task);
   _grad[task].derivative = get_manager()->messenger().receive_from_worker_on_queue<double>(worker_id);
   _grad[task].second_derivative = get_manager()->messenger().receive_from_worker_on_queue<double>(worker_id);
   _grad[task].step_size = get_manager()->messenger().receive_from_worker_on_queue<double>(worker_id);
//   std::cout << "queue receives for job id " << id << " from worker " << worker_id << " task result task=" << task
//             << " grad=" << _grad[task].derivative << " g2=" << _grad[task].second_derivative << " gstep=" << _grad[task].step_size
//             << std::endl;
}

void LikelihoodGradientJob::send_back_results_from_queue_to_master()
{
   get_manager()->messenger().send_from_queue_to_master(completed_task_ids.size());
//   std::cout << "sending from queue to master " << completed_task_ids.size() << " results: ";
   for (auto task : completed_task_ids) {
      get_manager()->messenger().send_from_queue_to_master(task, _grad[task].derivative, _grad[task].second_derivative,
                                                           _grad[task].step_size);
//      std::cout << "\t"
//                << "task=" << task << " grad=" << _grad[task].derivative << " g2=" << _grad[task].second_derivative
//                << " gstep=" << _grad[task].step_size;
   }
//   std::cout << std::endl;
}

void LikelihoodGradientJob::clear_results()
{
   completed_task_ids.clear();
}

void LikelihoodGradientJob::receive_results_on_master()
{
   auto get_time = []() {
      return std::chrono::duration_cast<std::chrono::nanoseconds>(
         std::chrono::high_resolution_clock::now().time_since_epoch())
         .count();
   };
   decltype(get_time()) t1, t2;
   t1 = get_time();
   std::size_t N_completed_tasks = get_manager()->messenger().receive_from_queue_on_master<std::size_t>();
   for (unsigned int sync_ix = 0u; sync_ix < N_completed_tasks; ++sync_ix) {
      std::size_t task = get_manager()->messenger().receive_from_queue_on_master<std::size_t>();
      _grad[task].derivative = get_manager()->messenger().receive_from_queue_on_master<double>();
      _grad[task].second_derivative = get_manager()->messenger().receive_from_queue_on_master<double>();
      _grad[task].step_size = get_manager()->messenger().receive_from_queue_on_master<double>();
   }
   t2 = get_time();
   printf("timestamps LikelihoodGradientJob::receive_results_on_master: %lld %lld\n", t1, t2);
}

//bool LikelihoodGradientJob::receive_task_result_on_master()
//{
//   std::size_t task = get_manager()->messenger().receive_from_worker_on_master<std::size_t>();
//   _grad[task] = get_manager()->messenger().receive_from_worker_on_master<RooFit::MinuitDerivatorElement>();
//   --N_tasks_at_workers;
//   bool job_completed = (N_tasks_at_workers == 0);
//   return job_completed;
//}

bool LikelihoodGradientJob::receive_task_result_on_master(const zmq::message_t & message)
{
   auto result = message.data<task_result_t>();
   _grad[result->task_id] = result->grad;
   --N_tasks_at_workers;
   bool job_completed = (N_tasks_at_workers == 0);
   return job_completed;
}


// END SYNCHRONIZATION FROM WORKERS TO MASTER

///////////////////////////////////////////////////////////////////////////////
/// Calculation stuff (mostly duplicates of RooGradMinimizerFcn code):

void LikelihoodGradientJob::run_derivator(unsigned int i_component) const
{
   // Calculate the derivative etc for these parameters
//   auto parameter_values = _minimizer->get_function_parameter_values();
   _grad[i_component] =
      _gradf.fast_partial_derivative(_minimizer->getMultiGenFcn(),
                                     _minimizer->fitter()->Config().ParamsSettings(), i_component, _grad[i_component]);
}

///////////////////////////////////////////////////////////////////////////////
/// copy pasted and adapted from old MP::GradMinimizerFcn:

//void LikelihoodGradientJob::update_workers_state()
//{
//   auto get_time = []() {
//      return std::chrono::duration_cast<std::chrono::nanoseconds>(
//         std::chrono::high_resolution_clock::now().time_since_epoch())
//         .count();
//   };
//   decltype(get_time()) t1, t2;
//   t1 = get_time();
//   // TODO optimization: only send changed parameters (now sending all)
//   std::size_t ix = 0;
//   RooFit::MultiProcess::M2Q msg = RooFit::MultiProcess::M2Q::update_real;
//   for (ix = 0; ix < static_cast<std::size_t>(_minimizer->getNPar()); ++ix) {
//      get_manager()->messenger().send_from_master_to_queue(msg, id, ix, _grad[ix].derivative, false);
//   }
//   for (ix = 0; ix < static_cast<std::size_t>(_minimizer->getNPar()); ++ix) {
//      get_manager()->messenger().send_from_master_to_queue(msg, id, ix + 1 * _minimizer->getNPar(), _grad[ix].second_derivative,
//                                                           false);
//   }
//   for (ix = 0; ix < static_cast<std::size_t>(_minimizer->getNPar()); ++ix) {
//      get_manager()->messenger().send_from_master_to_queue(msg, id, ix + 2 * _minimizer->getNPar(), _grad[ix].step_size,
//                                                           false);
//   }
//
////   ix = 0;
////   auto parameter_values = _minimizer->get_function_parameter_values();
////   for (auto &parameter : parameter_values) {
//   for (ix = 0; ix < static_cast<std::size_t>(_minimizer->getNPar()); ++ix) {
//      get_manager()->messenger().send_from_master_to_queue(msg, id, ix + 3 * _minimizer->getNPar(), /*parameter*/ minuit_internal_x_[ix], false);
////      ++ix;
//   }
//   t2 = get_time();
//   printf("timestamps LikelihoodGradientJob::update_workers_state: %lld %lld\n", t1, t2);
//}

//void LikelihoodGradientJob::update_workers_state()
//{
//   // TODO optimization: only send changed parameters (now sending all)
//   std::size_t ix;
//   for (ix = 0; ix < static_cast<std::size_t>(_minimizer->getNPar()); ++ix) {
//      get_manager()->messenger().publish_from_master_to_workers(_grad[ix].derivative);
//   }
//   for (ix = 0; ix < static_cast<std::size_t>(_minimizer->getNPar()); ++ix) {
//      get_manager()->messenger().publish_from_master_to_workers(_grad[ix].second_derivative);
//   }
//   for (ix = 0; ix < static_cast<std::size_t>(_minimizer->getNPar()); ++ix) {
//      get_manager()->messenger().publish_from_master_to_workers(_grad[ix].step_size);
//   }
//   for (ix = 0; ix < static_cast<std::size_t>(_minimizer->getNPar()); ++ix) {
//      get_manager()->messenger().publish_from_master_to_workers(minuit_internal_x_[ix]);
//   }
//}

void LikelihoodGradientJob::update_workers_state()
{
   // TODO optimization: only send changed parameters (now sending all)
   get_manager()->messenger().publish_from_master_to_workers(id);
   zmq::message_t gradient_message(_grad.begin(), _grad.end());
   get_manager()->messenger().publish_from_master_to_workers(std::move(gradient_message));
   zmq::message_t minuit_internal_x_message(minuit_internal_x_.begin(), minuit_internal_x_.end());
   get_manager()->messenger().publish_from_master_to_workers(std::move(minuit_internal_x_message));
}


void LikelihoodGradientJob::calculate_all()
{
//   auto get_time = []() {
//      return std::chrono::duration_cast<std::chrono::nanoseconds>(
//                std::chrono::high_resolution_clock::now().time_since_epoch())
//         .count();
//   };
//   decltype(get_time()) t1, t2;

//   std::cout << "BABBELBOX" << std::endl;

   if (get_manager()->process_manager().is_master()) {
//      std::cout << "HAAAAAA" << std::endl;
      // do Grad, G2 and Gstep here and then just return results from the
      // separate functions below

      // update parameters and object states that changed since last calculation (or creation if first time)
//      t1 = get_time();
      update_workers_state();
//      t2 = get_time();

//      printf("wallclock [master] update_workers_state: %f\n", (t2 - t1) / 1.e9);

//      t1 = get_time();
      // master fills queue with tasks
      for (std::size_t ix = 0; ix < N_tasks; ++ix) {
         MultiProcess::JobTask job_task(id, ix);
         get_manager()->queue().add(job_task);
      }
      N_tasks_at_workers = N_tasks;
//      t2 = get_time();

//      printf("wallclock [master] put job tasks in queue: %f\n", (t2 - t1) / 1.e9);

      // wait for task results back from workers to master (put into _grad)
//      t1 = get_time();
      gather_worker_results();
//      t2 = get_time();

//      printf("wallclock [master] gather_worker_results: %f\n", (t2 - t1) / 1.e9);

      calculation_is_clean->gradient = true;
   }
}

void LikelihoodGradientJob::fill_gradient(double *grad)
{
   if (get_manager()->process_manager().is_master()) {
      if (!calculation_is_clean->gradient) {
         calculate_all();
      }

      // TODO: maybe make a flag to avoid this copy operation, but maybe not worth the effort
      // put the results from _grad into *grad
      for (Int_t ix = 0; ix < _minimizer->getNPar(); ++ix) {
         grad[ix] = _grad[ix].derivative;
      }
   }
}

void LikelihoodGradientJob::fill_second_derivative(double *g2)
{
   if (get_manager()->process_manager().is_master()) {
      if (!calculation_is_clean->gradient) {
         calculate_all();
      }

      // TODO: maybe make a flag to avoid this copy operation, but maybe not worth the effort
      // put the results from _grad into *grad
      for (Int_t ix = 0; ix < _minimizer->getNPar(); ++ix) {
         g2[ix] = _grad[ix].second_derivative;
      }
   }
}

void LikelihoodGradientJob::fill_step_size(double *gstep)
{
   if (get_manager()->process_manager().is_master()) {
      if (!calculation_is_clean->gradient) {
         calculate_all();
      }

      // TODO: maybe make a flag to avoid this copy operation, but maybe not worth the effort
      // put the results from _grad into *grad
      for (Int_t ix = 0; ix < _minimizer->getNPar(); ++ix) {
         gstep[ix] = _grad[ix].step_size;
      }
   }
}

void LikelihoodGradientJob::update_minuit_internal_parameter_values(const std::vector<double>& minuit_internal_x)
{
   minuit_internal_x_ = minuit_internal_x;
}

bool LikelihoodGradientJob::uses_minuit_internal_values()
{
   return true;
}

void LikelihoodGradientJob::update_state()
{
//   std::size_t ix;
//
//   for (ix = 0; ix < static_cast<std::size_t>(_minimizer->getNPar()); ++ix) {
//      _grad[ix].derivative = get_manager()->messenger().receive_from_master_on_worker<double>();
//   }
//   for (ix = 0; ix < static_cast<std::size_t>(_minimizer->getNPar()); ++ix) {
//      _grad[ix].second_derivative = get_manager()->messenger().receive_from_master_on_worker<double>();
//   }
//   for (ix = 0; ix < static_cast<std::size_t>(_minimizer->getNPar()); ++ix) {
//      _grad[ix].step_size = get_manager()->messenger().receive_from_master_on_worker<double>();
//   }
//   for (ix = 0; ix < static_cast<std::size_t>(_minimizer->getNPar()); ++ix) {
//      minuit_internal_x_[ix] = get_manager()->messenger().receive_from_master_on_worker<double>();
//   }

   auto gradient_message = get_manager()->messenger().receive_from_master_on_worker<zmq::message_t>();
   auto gradient_message_begin = gradient_message.data<RooFit::MinuitDerivatorElement>();
   auto gradient_message_end = gradient_message_begin + gradient_message.size()/sizeof(RooFit::MinuitDerivatorElement);
   std::copy(gradient_message_begin, gradient_message_end, _grad.begin());

   auto minuit_internal_x_message = get_manager()->messenger().receive_from_master_on_worker<zmq::message_t>();
   auto minuit_internal_x_message_begin = minuit_internal_x_message.data<double>();
   auto minuit_internal_x_message_end = minuit_internal_x_message_begin + minuit_internal_x_message.size()/sizeof(double);
   std::copy(minuit_internal_x_message_begin, minuit_internal_x_message_end, minuit_internal_x_.begin());

   _gradf.setup_differentiate(_minimizer->getMultiGenFcn(), minuit_internal_x_.data(), _minimizer->fitter()->Config().ParamsSettings());
}


} // namespace TestStatistics
} // namespace RooFit
