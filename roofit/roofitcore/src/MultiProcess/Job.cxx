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

#include <unistd.h> // getpid

#include <RooMsgService.h>

#include <MultiProcess/Job.h>
#include <MultiProcess/TaskManager.h>
#include <MultiProcess/messages.h>

namespace RooFit {
  namespace MultiProcessV1 {
    Job::Job(std::size_t _N_workers) : N_workers(_N_workers) {}
    Job::Job(const Job & other) :
      N_workers(other.N_workers),
      waiting_for_queued_tasks(other.waiting_for_queued_tasks),
      _manager(other._manager)
    {}

    double Job::call_double_const_method(std::string /*key*/) {
      throw std::logic_error("call_double_const_method not implemented for this Job");
    }

    // This default sends back only one double as a result; can be overloaded
    // e.g. for RooAbsCategorys, for tuples, etc. The queue_loop and master
    // process must implement corresponding result receivers.
    void Job::send_back_task_result_from_worker(std::size_t task) {
      double result = get_task_result(task);
      get_manager()->send_from_worker_to_queue(id, task, result);
    }

    // static function
    void Job::worker_loop() {
      assert(TaskManager::instance()->is_worker());
      worker_loop_running = true;
      bool carry_on = true;
      Task task;
      std::size_t job_id;
      Q2W message_q2w;

      // use a flag to not ask twice
      bool dequeue_acknowledged = true;

      auto get_time = [](){return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();};

      while (carry_on) {
        decltype(get_time()) t1 = 0, t2 = 0, t3 = 0;

        // try to dequeue a task
        if (dequeue_acknowledged) {  // don't ask twice
          t1 = get_time();
          TaskManager::instance()->send_from_worker_to_queue(W2Q::dequeue);
          dequeue_acknowledged = false;
        }

        // receive handshake
        message_q2w = TaskManager::instance()->receive_from_queue_on_worker<Q2W>();

        switch (message_q2w) {
          case Q2W::terminate: {
            carry_on = false;
            break;
          }

          case Q2W::dequeue_rejected: {
            t2 = get_time();
            oocxcoutD((TObject*)nullptr,Benchmarking2) << "no work: worker " << TaskManager::instance()->get_worker_id() << " asked at " << t1 << " and got rejected at " << t2 << std::endl;

            dequeue_acknowledged = true;
            break;
          }
          case Q2W::dequeue_accepted: {
            dequeue_acknowledged = true;
            job_id = TaskManager::instance()->receive_from_queue_on_worker<std::size_t>();
            task = TaskManager::instance()->receive_from_queue_on_worker<Task>();

            t2 = get_time();
            TaskManager::get_job_object(job_id)->evaluate_task(task);

            t3 = get_time();
            oocxcoutD((TObject*)nullptr,Benchmarking2) << "job done: worker " << TaskManager::instance()->get_worker_id() << " asked at " << t1 << ", started at " << t2 << " and finished at " << t3 << std::endl;

            TaskManager::instance()->send_from_worker_to_queue(W2Q::send_result);
            TaskManager::get_job_object(job_id)->send_back_task_result_from_worker(task);

            message_q2w = TaskManager::instance()->receive_from_queue_on_worker<Q2W>();
            if (message_q2w != Q2W::result_received) {
              std::cerr << "worker " << getpid() << " sent result, but did not receive Q2W::result_received handshake! Got " << message_q2w << " instead." << std::endl;
              throw std::runtime_error("");
            }
            break;
          }

          // previously this was non-work mode

          case Q2W::update_real: {
            auto t1 = get_time();

            job_id = TaskManager::instance()->receive_from_queue_on_worker<std::size_t>();
            std::size_t ix = TaskManager::instance()->receive_from_queue_on_worker<std::size_t>();
            double val = TaskManager::instance()->receive_from_queue_on_worker<double>();
            bool is_constant = TaskManager::instance()->receive_from_queue_on_worker<bool>();
            TaskManager::get_job_object(job_id)->update_real(ix, val, is_constant);

            auto t2 = get_time();
            oocxcoutD((TObject*)nullptr,Benchmarking1) << "update_real on worker " << TaskManager::instance()->get_worker_id() << ": " << (t2 - t1)/1.e9 << "s (from " << t1 << " to " << t2 << "ns)" << std::endl;

            break;
          }

          case Q2W::call_double_const_method: {
            job_id = TaskManager::instance()->receive_from_queue_on_worker<std::size_t>();
            std::string key = TaskManager::instance()->receive_from_queue_on_worker<std::string>();
            Job * job = TaskManager::get_job_object(job_id);
//                double (* method)() = job->get_double_const_method(key);
            double result = job->call_double_const_method(key);
            TaskManager::instance()->send_from_worker_to_queue(result);
            break;
          }

          case Q2W::flush_ostreams: {
            TaskManager::instance()->flush_ostreams();
            break;
          }

          case Q2W::result_received: {
            std::cerr << "In worker_loop: " << message_q2w << " message received, but should only be received as handshake!" << std::endl;
            break;
          }

        }
      }
    }

    TaskManager* Job::get_manager() {
      if (!_manager) {
        _manager = TaskManager::instance(N_workers);

//        _manager->identify_processes();
//        sleep(10);
      }

      // N.B.: must check for activation here, otherwise get_manager is not callable
      //       from queue loop!
      if (!_manager->is_activated()) {
        _manager->activate();
      }

      if (!worker_loop_running && _manager->is_worker()) {
        Job::worker_loop();
        _manager->close_worker_connections();
        // flush remaining output
        _manager->flush_ostreams();
//        std::cout << "exiting worker process " << getpid() << std::endl;
        std::_Exit(0);
      }

      return _manager;
    }

    // initialize static members
    bool Job::worker_loop_running = false;

  } // namespace MultiProcess
} // namespace RooFit
