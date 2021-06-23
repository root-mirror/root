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

#include <stdexcept>  // logic_error
#include <sstream>
#include <unistd.h> // getpid
// for cpu affinity
#if !defined(__APPLE__) && !defined(_WIN32)
#include <sched.h>
#endif

#include <RooMsgService.h>

#include <ROOT/RMakeUnique.hxx>  // make_unique in C++11
#include <MultiProcess/TaskManager.h>
#include <MultiProcess/Job.h>
#include <MultiProcess/messages.h>
#include <MultiProcess/util.h>
#include <RooFit_ZMQ/ZeroMQPoller.h>
#include <algorithm>  // for_each

namespace RooFit {
  namespace MultiProcessV1 {

    void set_socket_immediate(ZmqLingeringSocketPtr<> & socket) {
      int optval = 1;
      socket->setsockopt(ZMQ_IMMEDIATE, &optval, sizeof(optval));
    }

    // static function
    TaskManager* TaskManager::instance(std::size_t N_workers) {
      if (!TaskManager::is_instantiated()) {
        assert(N_workers != 0);
        _instance = std::make_unique<TaskManager>(N_workers);
      }
      // these sanity checks no longer make sense with the worker_pipes only being maintained on the queue process
//      } else {
//        // some sanity checks
//        if(_instance->is_master() && N_workers != _instance->worker_pipes.size()) {
//          std::cerr << "On PID " << getpid() << ": N_workers != tmp->worker_pipes.size())! N_workers = " << N_workers << ", tmp->worker_pipes.size() = " << _instance->worker_pipes.size() << std::endl;
//          throw std::logic_error("");
//        } else if (_instance->is_worker()) {
//          if (_instance->get_worker_id() + 1 != _instance->worker_pipes.size()) {
//            std::cerr << "On PID " << getpid() << ": tmp->get_worker_id() + 1 != tmp->worker_pipes.size())! tmp->get_worker_id() = " << _instance->get_worker_id() << ", tmp->worker_pipes.size() = " << _instance->worker_pipes.size() << std::endl;
//            throw std::logic_error("");
//          }
//        }
//      }
      return _instance.get();
    }

    // static function
    TaskManager* TaskManager::instance() {
      if (!TaskManager::is_instantiated()) {
        throw std::runtime_error("in TaskManager::instance(): no instance was created yet! Call TaskManager::instance(std::size_t N_workers) first.");
      }
      return _instance.get();
    }

    // static function
    bool TaskManager::is_instantiated() {
      return static_cast<bool>(_instance);
    }


    void TaskManager::identify_processes() {
      // identify yourselves (for debugging)
      if (!(_is_master || _is_queue)) {
        std::cout << "I'm a worker, PID " << getpid() << '\n';
      } else if (_is_master) {
        std::cout << "I'm master, PID " << getpid() << '\n';
      } else if (_is_queue) {
        std::cout << "I'm queue, PID " << getpid() << '\n';
      }
    }

    // constructor
    // Don't construct IPQM objects manually, use the static instance if
    // you need to run multiple jobs.
    TaskManager::TaskManager(std::size_t N_workers) : N_workers(N_workers) {
      // This class defines three types of processes:
      // 1. master: the initial main process. It defines and enqueues tasks
      //    and processes results.
      // 2. workers: a pool of processes that will try to take tasks from the
      //    queue. These are first forked from master.
      // 3. queue: communication between the other types (necessarily) goes
      //    through this process. This process runs the queue_loop and
      //    maintains the queue of tasks. It is forked last and initialized
      //    with third BidirMMapPipe parameter false, which makes it the
      //    process that manages all pipes, though the pool of pages remains
      //    on the master process.
      // The reason for using this layout is that we use BidirMMapPipe for
      // forking and communicating between processes, and BidirMMapPipe only
      // supports forking from one process, not from an already forked
      // process (if forked using BidirMMapPipe). The latter layout would
      // allow us to fork first the queue from the main process and then fork
      // the workers from the queue, which may feel more natural.

      initialize_processes();
    }

    void TaskManager::initialize_processes(bool cpu_pinning) {
      // Initialize processes;
      // ... first workers:
      worker_pids.resize(N_workers);
      pid_t child_pid {};
      for (std::size_t ix = 0; ix < N_workers; ++ix) {
        child_pid = fork();
        if (!child_pid) {  // we're on the worker
          worker_id = ix;
          break;
        } else {          // we're on master
          worker_pids[ix] = child_pid;
        }
      }

      // ... then queue:
      if (child_pid) {  // we're on master
        queue_pid = fork();
        if (!queue_pid) { // we're now on queue
          _is_queue = true;
        } else {
          _is_master = true;
        }
      }

      // after all forks, create zmq connections (zmq context is automatically created in the ZeroMQSvc class and maintained as singleton)
      try {
        if (is_master()) {
          mq_socket.reset(zmqSvc().socket_ptr(zmq::PAIR)); //REQ));
          set_socket_immediate(mq_socket);
          mq_socket->bind("ipc:///tmp/roofitMP_master_queue");
  //        mq_socket->bind("tcp://*:55555");
        } else if (is_queue()) {
          // first the queue-worker sockets
          qw_sockets.resize(N_workers); // do resize instead of reserve so that the unique_ptrs are initialized (to nullptr) so that we can do reset below, alternatively you can do push/emplace_back with move or something
          for (std::size_t ix = 0; ix < N_workers; ++ix) {
            qw_sockets[ix].reset(zmqSvc().socket_ptr(zmq::PAIR)); //REP));
            set_socket_immediate(qw_sockets[ix]);
            std::stringstream socket_name;
            socket_name << "ipc:///tmp/roofitMP_queue_worker_" << ix;
  //          socket_name << "tcp://*:" << 55556 + ix;
            qw_sockets[ix]->bind(socket_name.str());
          }
          // then the master-queue socket
          mq_socket.reset(zmqSvc().socket_ptr(zmq::PAIR)); //REP));
          set_socket_immediate(mq_socket);
          mq_socket->connect("ipc:///tmp/roofitMP_master_queue");
  //        mq_socket->connect("tcp://127.0.0.1:55555");
        } else if (is_worker()) {
          this_worker_qw_socket.reset(zmqSvc().socket_ptr(zmq::PAIR)); //REQ));
          set_socket_immediate(this_worker_qw_socket);
          std::stringstream socket_name;
          socket_name << "ipc:///tmp/roofitMP_queue_worker_" << worker_id;
  //        socket_name << "tcp://127.0.0.1:" << 55556 + worker_id;
          this_worker_qw_socket->connect(socket_name.str());
        } else {
          // should never get here
          throw std::runtime_error("TaskManager::initialize_processes: I'm neither master, nor queue, nor a worker");
        }
      } catch (zmq::error_t& e) {
        std::cerr << e.what() << " -- errnum: " << e.num() << std::endl;
        throw;
      };

      if (cpu_pinning) {
        #if defined(__APPLE__)
        static bool affinity_warned = false;
        if (is_master() & !affinity_warned) {
          ooccoutD(static_cast<RooAbsArg*>(nullptr),Eval) << "CPU affinity cannot be set on macOS";
        }
        #elif defined(_WIN32)
        if (is_master()) std::cerr << "WARNING: CPU affinity setting not implemented on Windows, continuing...\n";
        #else
        cpu_set_t mask;
        // zero all bits in mask
        CPU_ZERO(&mask);
        // set correct bit
        std::size_t set_cpu;
        if (is_master()) {
          set_cpu = N_workers + 1;
        } else if (is_queue()) {
          set_cpu = N_workers;
        } else {
          set_cpu = worker_id;
        }
        CPU_SET(set_cpu, &mask);
        /* sched_setaffinity returns 0 in success */

        if (sched_setaffinity(0, sizeof(mask), &mask) == -1) {
          std::cerr << "WARNING: Could not set CPU affinity, continuing...\n";
        } else {
          std::cerr << "CPU affinity set to cpu " << set_cpu << " in process " << getpid() << '\n';
        }
        #endif
      }

//      identify_processes();

      processes_initialized = true;
    }


    TaskManager::~TaskManager() {
      ooccoutD(static_cast<RooAbsArg*>(nullptr),Eval) << "destroying TaskManager on PID " << getpid() << (is_worker() ? " worker" : (is_queue()? " queue" : " master"));
      // The TM instance gets created by some Job. Once all Jobs are gone, the
      // TM will get destroyed. In this case, the job_objects map should have
      // been emptied. This check makes sure:
      assert(TaskManager::job_objects.empty());
      terminate();
    }


    // static function
    // returns job_id for added job_object
    std::size_t TaskManager::add_job_object(Job *job_object) {
      if (TaskManager::is_instantiated()) {
        if (_instance->is_activated()) {
          std::stringstream ss;
          ss << "Cannot add Job to activated TaskManager instantiation (forking has already taken place)! Instance object at raw ptr " << _instance.get();
          throw std::logic_error("Cannot add Job to activated TaskManager instantiation (forking has already taken place)! Call terminate() on the instance before adding new Jobs.");
        }
      }
      std::size_t job_id = job_counter++;
      job_objects[job_id] = job_object;
      return job_id;
    }

    // static function
    Job* TaskManager::get_job_object(std::size_t job_object_id) {
      return job_objects[job_object_id];
    }

    // static function
    bool TaskManager::remove_job_object(std::size_t job_object_id) {
      bool removed_succesfully = job_objects.erase(job_object_id) == 1;
      if (job_objects.empty()) {
        _instance.reset(nullptr);
      }
      return removed_succesfully;
    }


    void TaskManager::terminate() noexcept {
      try {
        if (is_master()) {
          send_from_master_to_queue(M2Q::terminate);
//          int period = 0;
//          mq_socket->setsockopt(ZMQ_LINGER, &period, sizeof(period));
          mq_socket.reset(nullptr);
          zmqSvc().close_context();
          queue_activated = false;
          shutdown_processes();
        }
//      } catch (const BidirMMapPipe::Exception& e) {
//        std::cerr << "WARNING: in TaskManager::terminate, something in BidirMMapPipe threw an exception! Message:\n\t" << e.what() << std::endl;
      } catch (const std::exception& e) {
        std::cerr << "WARNING: something in TaskManager::terminate threw an exception! Original exception message:\n" << e.what() << '\n';
      }
    }


    void TaskManager::shutdown_processes() {
      if (is_master()) {
//        send_from_master_to_queue(M2Q::terminate);
//        *queue_pipe << M2Q::terminate << BidirMMapPipe::flush;

        // first wait for the workers that will be terminated by the queue
        for (auto pid : worker_pids) {
          wait_for_child(pid, true, 5);
        }
        // then wait for the queue
        wait_for_child(queue_pid, true, 5);

//        // delete queue_pipe (not worker_pipes, only on queue process!)
//        // CAUTION: the following invalidates a possibly created PollVector
//        queue_pipe.reset(); // sets to nullptr

//        for (auto it = worker_pids.begin(); it != worker_pids.end(); ++it) {
//          BidirMMapPipe::wait_for_child(*it, true);
//        }
      }

      processes_initialized = false;
      _is_master = false;
    }

    void TaskManager::close_worker_connections() {
//      int period = 0;
      if (is_worker()) {
//        this_worker_qw_socket->setsockopt(ZMQ_LINGER, &period, sizeof(period));
        this_worker_qw_socket.reset(nullptr);
        zmqSvc().close_context();
      } else if (is_queue()) {
        for (std::size_t worker_ix = 0ul; worker_ix < N_workers; ++worker_ix) {
//          qw_sockets[worker_ix]->setsockopt(ZMQ_LINGER, &period, sizeof(period));
          qw_sockets[worker_ix].reset(nullptr);
        }
      }
    }

    void TaskManager::terminate_workers() {
      if (is_queue()) {
//        for (std::unique_ptr<BidirMMapPipe> &worker_pipe : worker_pipes) {
//          *worker_pipe << Q2W::terminate << BidirMMapPipe::flush;
//        for(auto& worker_socket : qw_sockets) {
//          zmqSvc().send(*worker_socket, Q2W::terminate);
        for(std::size_t worker_ix = 0; worker_ix < N_workers; ++worker_ix) {
          send_from_queue_to_worker(worker_ix, Q2W::terminate);
//          int retval = worker_pipe->close();
//          if (0 != retval) {
//            std::cerr << "error terminating worker_pipe for worker with PID " << worker_pipe->pidOtherEnd() << "; child return value is " << retval << std::endl;
//          }
        }
        close_worker_connections();
      }
    }


    // start message loops on child processes and quit processes afterwards
    void TaskManager::activate() {
//      std::cout << "activating" << std::endl;
      // should be called soon after creation of this object, because everything in
      // between construction and activate gets executed both on the master process
      // and on the slaves
      if (!processes_initialized) {
//        std::cout << "intializing" << std::endl;
        initialize_processes();
      }

      queue_activated = true; // set on all processes, master, queue and slaves

      if (is_queue()) {
        queue_loop();
        terminate_workers();
//        int period = 0;
//        mq_socket->setsockopt(ZMQ_LINGER, &period, sizeof(period));
        mq_socket.reset(nullptr);   // delete/close master-queue socket from queue side
        zmqSvc().close_context();
        std::_Exit(0);
      }
    }


    bool TaskManager::is_activated() {
      return queue_activated;
    }


    // CAUTION:
    // this function returns a vector of pointers that may get invalidated by
    // the terminate function!
//    BidirMMapPipe::PollVector TaskManager::get_poll_vector() {
//      BidirMMapPipe::PollVector poll_vector;
//      poll_vector.reserve(1 + worker_pipes.size());
//      poll_vector.emplace_back(queue_pipe.get(), BidirMMapPipe::Readable);
//      for (std::unique_ptr<BidirMMapPipe>& pipe : worker_pipes) {
//        poll_vector.emplace_back(pipe.get(), BidirMMapPipe::Readable);
//      }
//      return poll_vector;
//    }


    bool TaskManager::process_queue_pipe_message(M2Q message) {
      bool carry_on = true;

      switch (message) {
        case M2Q::terminate: {
          carry_on = false;
        }
        break;

        case M2Q::enqueue: {
          // enqueue task
          auto job_object_id = receive_from_master_on_queue<std::size_t>();
          auto task = receive_from_master_on_queue<Task>();
          JobTask job_task(job_object_id, task);
          to_queue(job_task);
          N_tasks++;
        }
        break;

        case M2Q::retrieve: {
          // retrieve task results after queue is empty and all
          // tasks have been completed
          if (queue.empty() && N_tasks_completed == N_tasks) {
            send_from_queue_to_master(Q2M::retrieve_accepted, job_objects.size());
            for (auto job_tuple : job_objects) {
              send_from_queue_to_master(job_tuple.first);  // job id
              job_tuple.second->send_back_results_from_queue_to_master();  // N_job_tasks, task_ids and results
              job_tuple.second->clear_results();
            }
            // reset number of received and completed tasks
            N_tasks = 0;
            N_tasks_completed = 0;
          } else {
            send_from_queue_to_master(Q2M::retrieve_rejected);  // handshake message: tasks not done yet, try again
          }
        }
        break;

        case M2Q::update_real: {
          auto get_time = [](){return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();};
          auto t1 = get_time();
          auto job_id = receive_from_master_on_queue<std::size_t>();
          auto ix = receive_from_master_on_queue<std::size_t>();
          auto val = receive_from_master_on_queue<double>();
          auto is_constant = receive_from_master_on_queue<bool>();
          for (std::size_t worker_ix = 0; worker_ix < N_workers; ++worker_ix) {
            send_from_queue_to_worker(worker_ix, Q2W::update_real, job_id, ix, val, is_constant);
          }
          auto t2 = get_time();
          oocxcoutD((TObject*)nullptr,Benchmarking1) << "update_real on queue: " << (t2 - t1)/1.e9 << "s (from " << t1 << " to " << t2 << "ns)" << std::endl;
        }
        break;

        case M2Q::call_double_const_method: {
          auto job_id = receive_from_master_on_queue<std::size_t>();
          auto worker_id_call = receive_from_master_on_queue<std::size_t>();
          auto key = receive_from_master_on_queue<std::string>();
          send_from_queue_to_worker(worker_id_call, Q2W::call_double_const_method, job_id, key);
          auto result = receive_from_worker_on_queue<double>(worker_id_call);
          send_from_queue_to_master(result);
        }
        break;

        case M2Q::flush_ostreams: {
          flush_ostreams();
        }
        break;
      }

      return carry_on;
    }


    void TaskManager::retrieve() {
      if (_is_master) {
        bool carry_on = true;
        while (carry_on) {
          send_from_master_to_queue(M2Q::retrieve);
          auto handshake = receive_from_queue_on_master<Q2M>();
          if (handshake == Q2M::retrieve_accepted) {
            carry_on = false;
            auto N_jobs = receive_from_queue_on_master<std::size_t>();
            for (std::size_t job_ix = 0; job_ix < N_jobs; ++job_ix) {
              auto job_object_id = receive_from_queue_on_master<std::size_t>();
              TaskManager::get_job_object(job_object_id)->receive_results_on_master();
            }
          }
        }
      }
    }


    double TaskManager::call_double_const_method(const std::string& method_key, std::size_t job_id, std::size_t worker_id_call) {
      send_from_master_to_queue(M2Q::call_double_const_method, job_id, worker_id_call, method_key);
      auto result = receive_from_queue_on_master<double>();
      return result;
    }

    // -- WORKER - QUEUE COMMUNICATION --

    void TaskManager::send_from_worker_to_queue() {
//      *this_worker_pipe << BidirMMapPipe::flush;
    }

    void TaskManager::send_from_queue_to_worker(std::size_t /*this_worker_id*/) {
//      *worker_pipes[this_worker_id] << BidirMMapPipe::flush;
    }

    // -- QUEUE - MASTER COMMUNICATION --

    void TaskManager::send_from_queue_to_master() {
//      *queue_pipe << BidirMMapPipe::flush;
    }

    void TaskManager::send_from_master_to_queue() {
      send_from_queue_to_master();
    }


    void TaskManager::process_worker_pipe_message(std::size_t this_worker_id, W2Q message) {
      switch (message) {
        case W2Q::dequeue: {
          // dequeue task
          JobTask job_task;
          if (from_queue(job_task)) {
            send_from_queue_to_worker(this_worker_id, Q2W::dequeue_accepted, job_task.first, job_task.second);
          } else {
            send_from_queue_to_worker(this_worker_id, Q2W::dequeue_rejected);
          }
          break;
        }

        case W2Q::send_result: {
          // receive back task result
          auto job_object_id = receive_from_worker_on_queue<std::size_t>(this_worker_id);
          auto task = receive_from_worker_on_queue<Task>(this_worker_id);
          TaskManager::get_job_object(job_object_id)->receive_task_result_on_queue(task, this_worker_id);
          send_from_queue_to_worker(this_worker_id, Q2W::result_received);
          N_tasks_completed++;
          break;
        }
      }
    }


    void TaskManager::queue_loop() {
      if (_is_queue) {
        bool carry_on = true;
        ZeroMQPoller poller;
        auto mq_index = poller.register_socket(*mq_socket, zmq::POLLIN);
        for(auto& s : qw_sockets) {
          poller.register_socket(*s, zmq::POLLIN);
        }

        while (carry_on) {
          // poll: wait until status change (-1: infinite timeout)
          auto poll_result = poller.poll(-1);
          // then process incoming messages from sockets
          for (auto readable_socket : poll_result) {
            // message comes from the master/queue socket (first element):
            if (readable_socket.first == mq_index) {
              auto message = receive_from_master_on_queue<M2Q>();
              carry_on = process_queue_pipe_message(message);
              // on terminate, also stop for-loop, no need to check other
              // sockets anymore:
              if (!carry_on) {
                break;
              }
            } else { // from a worker socket
              auto this_worker_id = readable_socket.first - 1;
              auto message = receive_from_worker_on_queue<W2Q>(this_worker_id);
              process_worker_pipe_message(this_worker_id, message);
            }
          }
        }
      }
    }


    // Have a worker ask for a task-message from the queue
    bool TaskManager::from_queue(JobTask &job_task) {
      if (queue.empty()) {
        return false;
      } else {
        job_task = queue.front();
        queue.pop();
        return true;
      }
    }


    // Enqueue a task
    void TaskManager::to_queue(JobTask job_task) {
      if (is_master()) {
        if (!queue_activated) {
          activate();
        }
        send_from_master_to_queue(M2Q::enqueue, job_task.first, job_task.second);
      } else if (is_queue()) {
        queue.push(job_task);
      } else {
        throw std::logic_error("calling Communicator::to_master_queue from slave process");
      }
    }


    bool TaskManager::is_master() {
      return _is_master;
    }

    bool TaskManager::is_queue() {
      return _is_queue;
    }

    bool TaskManager::is_worker() {
      return !(_is_master || _is_queue);
    }

    std::size_t TaskManager::get_worker_id() {
      return worker_id;
    }

    void TaskManager::flush_ostreams() {
      if (is_master()) {
        send_from_master_to_queue(M2Q::flush_ostreams);
      } else if (is_queue()) {
        // flush output
        std::cout << std::flush;
        std::clog << std::flush;
        for (std::size_t worker_ix = 0; worker_ix < N_workers; ++worker_ix) {
          send_from_queue_to_worker(worker_ix, Q2W::flush_ostreams);
        }
      } else if (is_worker()) {
        // flush output
        std::cout << std::flush;
        std::clog << std::flush;
      } else {
        throw std::logic_error("calling TaskManager::flush_ostreams from unknown type of process (not master, queue or worker)");
      }
    }


    // initialize static members
    std::map<std::size_t, Job *> TaskManager::job_objects;
    std::size_t TaskManager::job_counter = 0;
    std::unique_ptr<TaskManager> TaskManager::_instance {nullptr};

  } // namespace MultiProcess
} // namespace RooFit