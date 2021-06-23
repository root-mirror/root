/*
 * Project: RooFit
 * Authors:
 *   PB, Patrick Bos, Netherlands eScience Center, p.bos@esciencecenter.nl
 *   IP, Inti Pelupessy, Netherlands eScience Center, i.pelupessy@esciencecenter.nl
 *
 * Copyright (c) 2016-2019, Netherlands eScience Center
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 */
#ifndef ROOT_ROOFIT_MultiProcess_Job_decl
#define ROOT_ROOFIT_MultiProcess_Job_decl

#include <string>
#include "RooFit_ZMQ/zmq.hxx"

namespace RooFit {
namespace MultiProcess {

// forward declaration
class JobManager;

/*
 * @brief interface class for defining the actual work that must be done
 *
 * Think of "job" as in "employment", e.g. the job of a baker, which
 * involves *tasks* like baking and selling bread. The Job must define the
 * tasks through its execution (evaluate_task), based on a task index argument.
 *
 * Classes inheriting from Job must implement the pure virtual methods:
 * - void evaluate_task(std::size_t task)
 * - void update_real(std::size_t ix, double val, bool is_constant)
 * - void update_bool(std::size_t ix, bool value)
 * - void receive_task_result_on_queue(std::size_t task, std::size_t worker_id)
 * - void send_back_results_from_queue_to_master()
 * - void clear_results()
 * - void receive_results_on_master()
 * - void send_back_task_result_from_worker(std::size_t task)
 *
 * The latter five deal with sending results back from workers to master.
 * update_real is used to sync parameters that changed on master to the workers
 * before the next set of tasks is queued; it must have a part that updates the
 * actual parameters on the worker process, but can optionally also be used to
 * send them from the master to the queue (Q2W is handled by the Queue::loop).
 * An example implementation:



 *
 * The type of result from each task is strongly dependent on the Job at hand
 * and so Job does not provide a default results vector or anything like this.
 * It is up to the inheriting class to implement this. We would have liked a
 * template parameter task_result_t, so that we could also provide a default
 * "boilerplate" calculate function to show a typical Job use-case of all the
 * above infrastructure. This is not trivial, because the JobManager has to
 * keep a list of Job pointers, so if there would be different template
 * instantiations of Jobs, this would complicate this list. As an example,
 * a boilerplate calculation could look something like this:
 *
return_type CertainJob::evaluate() {
   if (get_manager()->process_manager().is_master()) {
      // update parameters that changed since last calculation (or creation if first time)
      for (size_t ix = 0; ix < changed_parameters.size(); ++ix) {  // changed_parameters would be a member of CertainJob
         auto msg = RooFit::MultiProcess::M2Q::update_real;
         get_manager()->messenger().send_from_master_to_queue(msg, id, ix,
                                                              changed_parameters[ix].getVal(),
                                                              changed_parameters[ix].isConstant());
      }

      // master fills queue with tasks
      for (std::size_t ix = 0; ix < N_tasks; ++ix) {  // N_tasks would also be a member of CertainJob
         JobTask job_task(id, ix);
         get_manager()->queue().add(job_task);
      }

      // wait for task results back from workers to master
      gather_worker_results();

      // possibly do some final manipulation of the results that were gathered
   }
   return result;  // or not, it could be stored and accessed later; also here, result would be a member of CertainJob
}
 *
 * Child classes should refrain from direct access to the JobManager instance
 * (through JobManager::instance), but rather use the here provided
 * Job::get_manager(). This function starts the worker_loop on the worker when
 * first called, meaning that the workers will not
 */
class Job {
public:
   explicit Job();
   Job(const Job &other);

   ~Job();

   virtual void evaluate_task(std::size_t task) = 0;
   virtual void update_real(std::size_t ix, double val, bool is_constant) = 0;
   virtual void update_bool(std::size_t ix, bool value) = 0;
   virtual void update_state();

   virtual void send_back_task_result_from_worker(std::size_t task) = 0;
   virtual void receive_task_result_on_queue(std::size_t task, std::size_t worker_id) = 0;
   virtual void send_back_results_from_queue_to_master() = 0;
   // after results have been retrieved, they may need to be cleared to ensure
   // they won't be retrieved the next time again, e.g. when using a map to
   // collect results; if not needed it can just be left empty
   virtual void clear_results() = 0;
   virtual void receive_results_on_master() = 0;
   virtual bool receive_task_result_on_master(const zmq::message_t & message) = 0;

   void gather_worker_results();

protected:
   JobManager *get_manager();

   std::size_t id;

private:
   // do not use _manager directly, it must first be initialized! use get_manager()
   JobManager *_manager = nullptr;
};

} // namespace MultiProcess
} // namespace RooFit

#endif // ROOT_ROOFIT_MultiProcess_Job_decl
