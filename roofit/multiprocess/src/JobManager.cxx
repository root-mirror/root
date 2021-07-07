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
#include <unistd.h>  // getpid

#include <thread>  // std::thread::hardware_concurrency()

#include "ROOT/RMakeUnique.hxx"  // make_unique in C++11

#include "RooFit/MultiProcess/ProcessManager.h"
#include "RooFit/MultiProcess/Messenger.h"
#include "RooFit/MultiProcess/JobManager.h"
#include "RooFit/MultiProcess/Job.h"
#include "RooFit/MultiProcess/Queue.h"  // complete type for JobManager::queue()
#include "RooFit/MultiProcess/worker.h"
#include "RooFit/MultiProcess/util.h"


namespace RooFit {
namespace MultiProcess {

// static function
JobManager* JobManager::instance() {
   if (!JobManager::is_instantiated()) {
      assert(default_N_workers != 0);
      _instance.reset(new JobManager(default_N_workers));  // can't use make_unique, because ctor is private
      _instance->messenger().test_connections(_instance->process_manager());
      // set send to non blocking on all processes after checking the connections are working:
      _instance->messenger().set_send_flag(ZMQ_DONTWAIT);
   }
   return _instance.get();
}


// static function
bool JobManager::is_instantiated() {
   return static_cast<bool>(_instance);
}

// (private) constructor
// Don't construct JobManager objects manually, use the static instance if
// you need to run multiple jobs.
JobManager::JobManager(std::size_t N_workers) {
   queue_ptr = std::make_unique<Queue>();
   process_manager_ptr = std::make_unique<ProcessManager>(N_workers);
   messenger_ptr = std::make_unique<Messenger>(*process_manager_ptr);
}

JobManager::~JobManager() {
   // The instance gets created by some Job. Once all Jobs are gone, the
   // JM will get destroyed. In this case, the job_objects map should have
   // been emptied. This check makes sure:
   assert(JobManager::job_objects.empty());
   messenger_ptr.reset(nullptr);
   process_manager_ptr.reset(nullptr);
   queue_ptr.reset(nullptr);
}


// static function
// returns job_id for added job_object
std::size_t JobManager::add_job_object(Job *job_object) {
   if (JobManager::is_instantiated()) {
      if (_instance->process_manager().is_initialized()) {
         std::stringstream ss;
         ss << "Cannot add Job to JobManager instantiation, forking has already taken place! Instance object at raw ptr " << _instance.get();
         throw std::logic_error("Cannot add Job to JobManager instantiation, forking has already taken place! Call terminate() on the instance before adding new Jobs.");
      }
   }
   std::size_t job_id = job_counter++;
   job_objects[job_id] = job_object;
   return job_id;
}

// static function
Job* JobManager::get_job_object(std::size_t job_object_id) {
   return job_objects[job_object_id];
}

// static function
bool JobManager::remove_job_object(std::size_t job_object_id) {
   bool removed_succesfully = job_objects.erase(job_object_id) == 1;
   if (job_objects.empty()) {
      _instance.reset(nullptr);
   }
   return removed_succesfully;
}


ProcessManager & JobManager::process_manager() const {
   return *process_manager_ptr;
}


Messenger & JobManager::messenger() const {
   return *messenger_ptr;
}


Queue & JobManager::queue() const {
   return *queue_ptr;
}


//void JobManager::retrieve() {
//   if (process_manager().is_master()) {
//      bool carry_on = true;
//      while (carry_on) {
//         messenger().send_from_master_to_queue(M2Q::retrieve);
//         try {
//            auto handshake = messenger().receive_from_queue_on_master<Q2M>();
//            switch (handshake) {
//            case Q2M::retrieve_accepted: {
//               carry_on = false;
//               auto N_jobs = messenger().receive_from_queue_on_master<std::size_t>();
//               for (std::size_t job_ix = 0; job_ix < N_jobs; ++job_ix) {
//                  auto job_object_id = messenger().receive_from_queue_on_master<std::size_t>();
//                  JobManager::get_job_object(job_object_id)->receive_results_on_master();
//               }
//            } break;
//            case Q2M::retrieve_later: {
//               carry_on = true;
//            } break;
//            case Q2M::retrieve_rejected: {
//               carry_on = false;
//               throw std::logic_error(
//                  "Master sent M2Q::retrieve, but queue had no tasks yet: Q2M::retrieve_rejected. Aborting!");
//            } break;
//            }
//         } catch (ZMQ::ppoll_error_t &e) {
//            zmq_ppoll_error_response response;
//            try {
//               response = handle_zmq_ppoll_error(e);
//            } catch (std::logic_error& e) {
//               printf("JobManager::retrieve got unhandleable ZMQ::ppoll_error_t\n");
//               throw;
//            }
//            if (response == zmq_ppoll_error_response::abort) {
//               throw std::logic_error("in JobManager::retrieve: master received a SIGTERM, aborting");
//            } else if (response == zmq_ppoll_error_response::unknown_eintr) {
//               printf("EINTR in JobManager::retrieve, continuing\n");
//               continue;
//            } else if (response == zmq_ppoll_error_response::retry) {
//               printf("EAGAIN from ppoll in JobManager::retrieve, continuing\n");
//               continue;
//            }
//         } catch (zmq::error_t& e) {
//            printf("unhandled zmq::error_t (not a ppoll_error_t) in JobManager::retrieve with errno %d: %s\n", e.num(), e.what());
//            throw;
//         }
//
//      }
//   }
//}


void JobManager::retrieve(std::size_t requesting_job_id) {
   if (process_manager().is_master()) {
//      auto get_time = []() {
//         return std::chrono::duration_cast<std::chrono::nanoseconds>(
//            std::chrono::high_resolution_clock::now().time_since_epoch())
//            .count();
//      };
//      decltype(get_time()) t1, t2;

      bool job_fully_retrieved = false;
      while (not job_fully_retrieved) {
         try {
//            t1 = get_time();
            auto task_result_message = messenger().receive_from_worker_on_master<zmq::message_t>();
            auto job_object_id = *reinterpret_cast<std::size_t *>(task_result_message.data());  // job_id must always be the first element of the result message!
            bool this_job_fully_retrieved = JobManager::get_job_object(job_object_id)->receive_task_result_on_master(task_result_message);
            if (requesting_job_id == job_object_id) {
               job_fully_retrieved = this_job_fully_retrieved;
            }
//            t2 = get_time();
//            printf("wallclock [master] retrieve got task %lu: %f\n", *(reinterpret_cast<std::size_t *>(task_result_message.data()) + 1), (t2 - t1) / 1.e9);
         } catch (ZMQ::ppoll_error_t &e) {
            zmq_ppoll_error_response response;
            try {
               response = handle_zmq_ppoll_error(e);
            } catch (std::logic_error& e) {
               printf("JobManager::retrieve got unhandleable ZMQ::ppoll_error_t\n");
               throw;
            }
            if (response == zmq_ppoll_error_response::abort) {
               throw std::logic_error("in JobManager::retrieve: master received a SIGTERM, aborting");
            } else if (response == zmq_ppoll_error_response::unknown_eintr) {
               printf("EINTR in JobManager::retrieve, continuing\n");
               continue;
            } else if (response == zmq_ppoll_error_response::retry) {
               printf("EAGAIN from ppoll in JobManager::retrieve, continuing\n");
               continue;
            }
         } catch (zmq::error_t& e) {
            printf("unhandled zmq::error_t (not a ppoll_error_t) in JobManager::retrieve with errno %d: %s\n", e.num(), e.what());
            throw;
         }
      }
   }
}


void JobManager::results_from_queue_to_master()
{
   assert(process_manager().is_queue());
   messenger().send_from_queue_to_master(Q2M::retrieve_accepted, JobManager::job_objects.size());
   for (auto job_tuple : JobManager::job_objects) {
      messenger().send_from_queue_to_master(job_tuple.first);                 // job id
      job_tuple.second->send_back_results_from_queue_to_master(); // N_job_tasks, task_ids and results
      job_tuple.second->clear_results();
   }
}


void JobManager::activate()
{
   // This function exists purely because activation from the constructor is
   // impossible; the constructor must return a constructed instance, which it
   // can't do if it's stuck in an infinite loop. This means the Job that first
   // creates the JobManager instance must also activate it (or any other user
   // of this class).
   // This should be called soon after creation of instance, because everything
   // between construction and activation gets executed both on the master
   // process and on the slaves.

   activated = true;

   if (process_manager().is_queue()) {
      queue().loop();
      std::_Exit(0);
   }

   if (!is_worker_loop_running() && process_manager().is_worker()) {
      RooFit::MultiProcess::worker_loop();
      std::_Exit(0);
   }
}


bool JobManager::is_activated() const
{
   return activated;
}

// initialize static members
std::map<std::size_t, Job *> JobManager::job_objects;
std::size_t JobManager::job_counter = 0;
std::unique_ptr<JobManager> JobManager::_instance {nullptr};
unsigned int JobManager::default_N_workers = std::thread::hardware_concurrency();

} // namespace MultiProcess
} // namespace RooFit
