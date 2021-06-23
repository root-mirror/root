/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id$
 * Authors:                                                                  *
 *   PB, Patrick Bos,     NL eScience Center, p.bos@esciencecenter.nl        *
 *   IP, Inti Pelupessy,  NL eScience Center, i.pelupessy@esciencecenter.nl  *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef ROOFIT_MULTIPROCESS_VECTOR_H
#define ROOFIT_MULTIPROCESS_VECTOR_H

#include <MultiProcess/Job.h>
#include <MultiProcess/TaskManager.h>
#include <RooListProxy.h>
#include <RooArgList.h>
#include <RooRealVar.h>

namespace RooFit {
  namespace MultiProcessV1 {

    // Vector defines an interface and communication machinery to build a
    // parallelized subclass of an existing non-concurrent numerical class that
    // can be expressed as a vector of independent sub-calculations.
    //
    // The way to use Vector is to define a new class that inherits from
    // Vector<ConcurrentClass>. Vector itself inherits from its template
    // parameter class, so the inheritance tree will become:
    //
    //   ParallelClass -> Vector<ConcurrentClass> -> ConcurrentClass
    //
    // where -> denotes "inherits from".
    //
    // While Vector is an abstract class, it provides some default method
    // implementations that define the result of each task as a single vector
    // element, i.e. a double precision floating point number. By overriding
    // these methods, the result of a Vector task can be redefined as, for
    // instance, a sum and a carry value for when the actual result of the task
    // is given by Kahan summation and the carry value needs to be propagated
    // back to the master process for calculating the total sum value. These
    // methods are:
    //
    // - receive_task_result_on_queue
    // - send_back_results_from_queue_to_master
    // - clear_results
    // - receive_results_on_master
    // - send_back_task_result_from_worker
    //
    // Note that the results vector<double> "defines" result type as well. In
    // some cases this may not be wanted, but we expect that in such a case the
    // Vector class itself is not suitable and a different Job subclass should
    // be defined. One clear exception is the case where a vector is wanted,
    // but not one of doubles but of another type. This use-case is
    // accommodated through the result_t template parameter.
    //
    // Classes inheriting from Vector must implement the pure virtual methods:
    // - void evaluate_task(std::size_t task)
    // - double get_task_result(std::size_t task)
    //
    // and can optionally override the virtual methods (in addition to the
    // result logistics methods mentioned above):
    // - void update_real(std::size_t ix, double val, bool is_constant)
    // - double call_double_const_method(std::string key)
    //
    template <typename Base, typename result_t = double>
    class Vector : public Base, public Job {
     public:
      template<typename... Targs>
      Vector(std::size_t _N_workers, Targs ...args) :
          Base(args...),
          Job(_N_workers)
      {
        id = TaskManager::add_job_object(this);
      }

      Vector(const Vector & other) :
          Base(other),
          Job(other),
          results(other.results),
          _vars(other._vars),
          _saveVars(other._saveVars),
          _forceCalc(other._forceCalc)
      {
        id = TaskManager::add_job_object(this);
      }

      ~Vector() {
        TaskManager::remove_job_object(id);
      }

      void update_real(std::size_t ix, double val, bool is_constant) override {
        if (get_manager()->is_worker()) {
          RooRealVar *rvar = (RooRealVar *) _vars.at(ix);
          rvar->setVal(static_cast<Double_t>(val));
          if (rvar->isConstant() != is_constant) {
            rvar->setConstant(static_cast<Bool_t>(is_constant));
          }
        }
      }

     protected:
      void gather_worker_results() {
        if (waiting_for_queued_tasks) {
          get_manager()->retrieve();
          waiting_for_queued_tasks = false;
        }
      }

      void receive_task_result_on_queue(std::size_t task, std::size_t worker_id) override {
        result_t result = get_manager()->template receive_from_worker_on_queue<result_t>(worker_id);
        results[task] = result;
      }

      void send_back_results_from_queue_to_master() override {
        get_manager()->send_from_queue_to_master(results.size());
        for (auto const &item : results) {
          get_manager()->send_from_queue_to_master(item.first, item.second);
        }
      }

      void clear_results() override {
        // empty results cache
        results.clear();
      }

      void receive_results_on_master() override {
        std::size_t N_job_tasks = get_manager()->template receive_from_queue_on_master<std::size_t>();
        for (std::size_t task_ix = 0ul; task_ix < N_job_tasks; ++task_ix) {
          std::size_t task_id = get_manager()->template receive_from_queue_on_master<std::size_t>();
          results[task_id] = get_manager()->template receive_from_queue_on_master<result_t>();
        }
      }

      // Here we define maps of functions that a subclass may want to call on
      // the worker process and have the result sent back to master. Each
      // function type needs custom implementation, so we only allow a selected
      // number of function pointer types. Templates could make the Vector
      // header slightly more compact, but the implementation would not change
      // much, so this explicit approach seems preferable.
      using double_const_method_t = double (Vector<Base>::*)() const;
      std::map<std::string, double_const_method_t> double_const_methods;
      double call_double_const_method(std::string key) override {
        return (this->*double_const_methods[key])();
      }
      // Another example would be:
      //   std::map<std::string, double (Vector<Base>::*)(double)> double_from_double_methods;
      // We leave this out for now, as we don't currently need it.
      //
      // Every type also needs a corresponding set of:
      // - messages from master to queue
      // - messages from queue to worker
      // - method to return the method pointer from the object to be able to
      //   call it from the static worker_loop.
      // The method could be implemented using templates, but the messages
      // cannot, further motivating the use of explicit implementation of
      // every specific case.
      //
      // Note that due to the relatively expensive nature of these calls,
      // they should be used only in non-work-mode.

      // -- members --
     protected:
      std::map<Task, result_t> results;
      
      // the following members are used for syncing, but might be replaced by
      // remote "mapped function calls" (see above):

      // used in NLLVar
      RooListProxy _vars;    // Variables
      RooArgList _saveVars;  // Copy of variables
      bool _forceCalc = false;

    };  // class Vector

  }  // namespace MultiProcessV1
}  // namespace RooFit

#endif //ROOFIT_MULTIPROCESS_VECTOR_H
