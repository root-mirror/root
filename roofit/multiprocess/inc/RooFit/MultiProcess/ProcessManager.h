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
#ifndef ROOT_ROOFIT_MultiProcess_ProcessManager
#define ROOT_ROOFIT_MultiProcess_ProcessManager

#include <sys/types.h>  // pid_t
#include <csignal>  // sig_atomic_t and for sigterm handling on child processes (in ProcessManager.cxx)
#include <vector>

// forward declaration
class Queue;

namespace RooFit {
namespace MultiProcess {

class ProcessManager {
   friend Queue;

public:
   explicit ProcessManager(std::size_t N_workers);
   ~ProcessManager();

   bool is_initialized() const;

   void terminate() noexcept;
//   void terminate_workers();
   void wait_for_sigterm_then_exit();

   bool is_master() const;
   bool is_queue() const;
   bool is_worker() const;
   std::size_t worker_id() const;
   std::size_t N_workers() const;

   void identify_processes() const;

   static void handle_sigterm(int signum);
   static bool sigterm_received();

   // for debugging/testing:
   pid_t get_queue_pid() {return queue_pid;}
   std::vector<pid_t> get_worker_pids() {return worker_pids;}

private:
   void initialize_processes(bool cpu_pinning = true);
   void shutdown_processes();

   bool _is_master = false;
   bool _is_queue = false;
   bool _is_worker = false;
   std::size_t _worker_id;
   std::size_t _N_workers;

   // master must wait for workers after completion, for which it needs their PIDs
   std::vector<pid_t> worker_pids;
   pid_t queue_pid;

   bool initialized = false;

   static volatile sig_atomic_t _sigterm_received;
};

} // namespace MultiProcess
} // namespace RooFit

#endif // ROOT_ROOFIT_MultiProcess_ProcessManager
