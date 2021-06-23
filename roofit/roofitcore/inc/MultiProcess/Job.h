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

#ifndef ROOFIT_MULTIPROCESS_JOB_H
#define ROOFIT_MULTIPROCESS_JOB_H

#include <string>
#include <memory>  // shared_ptr

namespace RooFit {
  namespace MultiProcessV1 {
    // forward declaration
    class TaskManager;

    /*
     * @brief interface class for defining the actual work that TaskManager must do
     *
     * Think of "job" as in "employment", e.g. the job of a baker, which
     * involves *tasks* like baking and selling bread. The Job must define the
     * tasks through its execution (evaluate_task) and returning its result
     * (get_task_result), based on a task index argument.
     *
     * Classes inheriting from Job must implement the pure virtual methods:
     * - void evaluate_task(std::size_t task)
     * - double get_task_result(std::size_t task)
     * - void update_real(std::size_t ix, double val, bool is_constant)
     * - void receive_task_result_on_queue(std::size_t task, std::size_t worker_id)
     * - void send_back_results_from_queue_to_master()
     * - void clear_results()
     * - void receive_results_on_master()
     *
     * and can optionally override the virtual methods:
     * - double call_double_const_method(std::string key)
     * - void send_back_task_result_from_worker(std::size_t task)
     *
     * Note that Job::call_double_const_method throws a logic_error if not
     * overridden.
     */
    class Job {
     public:
      explicit Job(std::size_t _N_workers);
      Job(const Job & other);

      virtual void evaluate_task(std::size_t task) = 0;
      // TODO: replace get_task_result return type (double) with something more flexible
      virtual double get_task_result(std::size_t task) = 0;
      virtual void update_real(std::size_t ix, double val, bool is_constant) = 0;

      virtual double call_double_const_method(std::string /*key*/);
      virtual void send_back_task_result_from_worker(std::size_t task);

      virtual void receive_task_result_on_queue(std::size_t task, std::size_t worker_id) = 0;
      // an example implementation that receives only one double and stores it in
      // a std::map<JobTask, double>:
      //      {
      //        double result;
      //        pipe >> job_object_id  >> result;
      //        pipe << Q2W::result_received << BidirMMapPipe::flush;
      //        JobTask job_task(job_object_id, task);
      //        results[job_task] = result;
      //      }

      virtual void send_back_results_from_queue_to_master() = 0;
      // after results have been retrieved, they should be cleared to ensure
      // they won't be retrieved the next time again
      virtual void clear_results() = 0;

      virtual void receive_results_on_master() = 0;

      TaskManager* get_manager();

      static void worker_loop();

     protected:
      std::size_t N_workers;
      std::size_t id;
      bool waiting_for_queued_tasks = false;

     private:
      // do not use _manager directly, it must first be initialized! use get_manager()
      TaskManager* _manager = nullptr;

      static bool worker_loop_running;
    };

  } // namespace MultiProcessV1
} // namespace RooFit
#endif //ROOFIT_MULTIPROCESS_JOB_H
