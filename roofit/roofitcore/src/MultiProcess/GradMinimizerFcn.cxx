/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id$
 * Authors:                                                                  *
 *   PB, Patrick Bos,     NL eScience Center, p.bos@esciencecenter.nl        *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#include <MultiProcess/GradMinimizer.h>
#include <MultiProcess/messages.h>
#include <stdexcept>

#include <RooTimer.h>

namespace RooFit {
  namespace MultiProcessV1 {
    GradMinimizerFcn::GradMinimizerFcn(RooAbsReal *funct, RooMinimizerGenericPtr context, std::size_t _N_workers,
                                       bool verbose) :
        RooFit::MultiProcessV1::Vector<RooGradMinimizerFcn>(_N_workers, funct, context, verbose) {
      N_tasks = NDim();
      completed_task_ids.reserve(N_tasks);
      // TODO: make sure that the full gradients are sent back so that the
      // derivator will depart from correct state next step everywhere!
    }

    // copy ctor (necessary for Clone)
    GradMinimizerFcn::GradMinimizerFcn(const GradMinimizerFcn& other) :
        RooFit::MultiProcessV1::Vector<RooGradMinimizerFcn>(other),
        N_tasks(other.N_tasks),
        completed_task_ids(other.completed_task_ids) {}

    ROOT::Math::IMultiGradFunction* GradMinimizerFcn::Clone() const {
      return new GradMinimizerFcn(*this) ;
    }

    // SYNCHRONIZATION FROM MASTER TO WORKERS

    void GradMinimizerFcn::update_state() {
      // TODO optimization: only send changed parameters (now sending all)
      RooFit::MultiProcessV1::M2Q msg = RooFit::MultiProcessV1::M2Q::update_real;
      for (std::size_t ix = 0; ix < NDim(); ++ix) {
        get_manager()->send_from_master_to_queue(msg, id, ix, _grad.Grad()(ix), false);
      }
      for (std::size_t ix = 0; ix < NDim(); ++ix) {
        get_manager()->send_from_master_to_queue(msg, id, ix + 1 * NDim(), _grad.G2()(ix), false);
      }
      for (std::size_t ix = 0; ix < NDim(); ++ix) {
        get_manager()->send_from_master_to_queue(msg, id, ix + 2 * NDim(), _grad.Gstep()(ix), false);
      }

      std::size_t ix = 0;
      for (auto parameter : _grad_params) {
        get_manager()->send_from_master_to_queue(msg, id, ix + 3 * NDim(), parameter, false);
        ++ix;
      }
    }

    void GradMinimizerFcn::update_real(std::size_t ix, double val, bool /*is_constant*/)
    {
       if (get_manager()->is_worker()) {
          // ix is defined in "flat" FunctionGradient space ix_dim * size + ix_component
          switch (ix / NDim()) {
          case 0: {
             _grad[ix % NDim()].derivative = val;
             break;
          }
          case 1: {
             _grad[ix % NDim()].second_derivative = val;
             break;
          }
          case 2: {
             _grad[ix % NDim()].step_size = val;
             break;
          }
          case 3: {
             sync_parameter(val, ix % NDim());
             break;
          }
          default: throw std::runtime_error("ix out of range in GradMinimizerFcn::update_real!");
          }
       }
    }

    // END SYNCHRONIZATION FROM MASTER TO WORKERS


    // SYNCHRONIZATION FROM WORKERS TO MASTER

    void GradMinimizerFcn::send_back_task_result_from_worker(std::size_t task) {
      get_manager()->send_from_worker_to_queue(id, task, _grad[task].derivative, _grad[task].second_derivative, _grad[task].step_size);
    }

    void GradMinimizerFcn::receive_task_result_on_queue(std::size_t task, std::size_t worker_id)
    {
       completed_task_ids.push_back(task);
       _grad[task].derivative = get_manager()->receive_from_worker_on_queue<double>(worker_id);
       _grad[task].second_derivative = get_manager()->receive_from_worker_on_queue<double>(worker_id);
       _grad[task].step_size = get_manager()->receive_from_worker_on_queue<double>(worker_id);
    }

    void GradMinimizerFcn::send_back_results_from_queue_to_master() {
      get_manager()->send_from_queue_to_master(completed_task_ids.size());
      for (auto task : completed_task_ids) {
        get_manager()->send_from_queue_to_master(task, _grad[task].derivative, _grad[task].second_derivative, _grad[task].step_size);
      }
    }

    void GradMinimizerFcn::clear_results() {
      completed_task_ids.clear();
    }

    void GradMinimizerFcn::receive_results_on_master()
    {
       std::size_t N_completed_tasks = get_manager()->receive_from_queue_on_master<std::size_t>();
       for (unsigned int sync_ix = 0u; sync_ix < N_completed_tasks; ++sync_ix) {
          std::size_t task = get_manager()->receive_from_queue_on_master<std::size_t>();
          _grad[task].derivative = get_manager()->receive_from_queue_on_master<double>();
          _grad[task].second_derivative = get_manager()->receive_from_queue_on_master<double>();
          _grad[task].step_size = get_manager()->receive_from_queue_on_master<double>();
       }
    }

    // END SYNCHRONIZATION FROM WORKERS TO MASTER


    // ACTUAL WORK

    void GradMinimizerFcn::evaluate_task(std::size_t task) {
      RooWallTimer timer;
      RooCPUTimer ctimer;
      run_derivator(task);
      ctimer.stop();
      timer.stop();
      oocxcoutD((TObject*)nullptr,Benchmarking1) << "worker_id: " << get_manager()->get_worker_id() << ", task: " << task << ", partial derivative time: " << timer.timing_s() << "s -- cputime: " << ctimer.timing_s() << "s" << std::endl;
    }

    double GradMinimizerFcn::get_task_result(std::size_t task) {
      // this is useless here
      return _grad.Grad()(task);
    }


    void GradMinimizerFcn::CalculateAll(const double *x) {
      auto get_time = [](){return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();};
      decltype(get_time()) t1, t2;

      if (get_manager()->is_master()) {
        // do Grad, G2 and Gstep here and then just return results from the
        // separate functions below
        bool was_not_synced = sync_parameters(x);
        if (was_not_synced) {
          // update parameters and object states that changed since last calculation (or creation if first time)
          t1 = get_time();
          update_state();
          t2 = get_time();

          RooWallTimer timer;

          // master fills queue with tasks
          for (std::size_t ix = 0; ix < N_tasks; ++ix) {
            JobTask job_task(id, ix);
            get_manager()->to_queue(job_task);
          }
          waiting_for_queued_tasks = true;

          // wait for task results back from workers to master (put into _grad)
          gather_worker_results();

          timer.stop();

          oocxcoutD((TObject*)nullptr,Benchmarking1) << "update_state: " << (t2 - t1)/1.e9 << "s (from " << t1 << " to " << t2 << "ns), gradient work: " << timer.timing_s() << "s" << std::endl;
        }
      }
    }

    void GradMinimizerFcn::Gradient(const double *x, double *grad) const {
      const_cast<GradMinimizerFcn *>(this)->mutable_Gradient(x, grad);
    }

    void GradMinimizerFcn::mutable_Gradient(const double *x, double *grad) {
      if (get_manager()->is_master()) {
        CalculateAll(x);

        // TODO: maybe make a flag to avoid this copy operation, but maybe not worth the effort
        // put the results from _grad into *grad
        for (std::size_t ix = 0; ix < NDim(); ++ix) {
          grad[ix] = _grad.Grad()(ix);
        }
      }
    }

    void GradMinimizerFcn::G2ndDerivative(const double *x, double *g2) const {
      const_cast<GradMinimizerFcn *>(this)->mutable_G2ndDerivative(x, g2);
    }

    void GradMinimizerFcn::mutable_G2ndDerivative(const double *x, double *g2) {
      if (get_manager()->is_master()) {
        CalculateAll(x);

        // TODO: maybe make a flag to avoid this copy operation, but maybe not worth the effort
        // put the results from _grad into *grad
        for (std::size_t ix = 0; ix < NDim(); ++ix) {
          g2[ix] = _grad.G2()(ix);
        }
      }
    }

    void GradMinimizerFcn::GStepSize(const double *x, double *gstep) const {
      const_cast<GradMinimizerFcn *>(this)->mutable_GStepSize(x, gstep);
    }

    void GradMinimizerFcn::mutable_GStepSize(const double *x, double *gstep) {
      if (get_manager()->is_master()) {
        CalculateAll(x);

        // TODO: maybe make a flag to avoid this copy operation, but maybe not worth the effort
        // put the results from _grad into *grad
        for (std::size_t ix = 0; ix < NDim(); ++ix) {
          gstep[ix] = _grad.Gstep()(ix);
        }
      }
    }

    // END ACTUAL WORK

//      RooGradMinimizerFcn(funct, context, verbose), N_workers() {}
  } // namespace MultiProcess
} // namespace RooFit
