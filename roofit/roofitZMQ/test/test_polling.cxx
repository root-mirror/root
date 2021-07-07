/*
 * Project: RooFit
 * Authors:
 *   PB, Patrick Bos, Netherlands eScience Center, p.bos@esciencecenter.nl
 *
 * Copyright (c) 2016-2019, Netherlands eScience Center
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 */

#include "gtest/gtest.h"

#include <unistd.h> // fork, usleep
#include <csignal>  // signal blocking
#include <cstring>  // strsignal()

#include <sstream>

#include "RooFit_ZMQ/ZeroMQSvc.h"
#include "RooFit_ZMQ/ZeroMQPoller.h"

static bool terminated = false;

void handle_sigterm(int signum)
{
   terminated = true;
   std::cout << "handled signal " << strsignal(signum) << " on PID " << getpid() << std::endl;
}

TEST(Polling, doublePoll)
{
   pid_t child_pid{0};
   do {
      child_pid = fork();
   } while (child_pid == -1); // retry if fork fails

   if (child_pid > 0) { // master
      sigset_t sigmask, sigmask_old;
      sigemptyset(&sigmask);
      sigaddset(&sigmask, SIGCHLD);
      sigprocmask(SIG_BLOCK, &sigmask, &sigmask_old);

      //      std::cout << "master PID: " << getpid() << std::endl;
      ZmqLingeringSocketPtr<> pusher, puller;
      pusher.reset(zmqSvc().socket_ptr(zmq::PUSH));
      pusher->bind("ipc:///tmp/ZMQ_test_fork_polling_M2C.ipc");
      puller.reset(zmqSvc().socket_ptr(zmq::PULL));
      puller->bind("ipc:///tmp/ZMQ_test_fork_polling_C2M.ipc");

      ZeroMQPoller poller1, poller2;
      poller1.register_socket(*puller, zmq::POLLIN);
      poller2.register_socket(*puller, zmq::POLLIN);

      // start test
      zmqSvc().send(*pusher, std::string("breaker breaker"));

      auto result1a = poller1.poll(-1);
      auto result1b = poller1.poll(-1);
      auto result2 = poller2.poll(-1);
      EXPECT_EQ(result1a.size(), result1b.size());
      EXPECT_EQ(result1a.size(), result2.size());

      auto receipt = zmqSvc().receive<int>(*puller, ZMQ_DONTWAIT);

      EXPECT_EQ(receipt, 1212);

      kill(child_pid, SIGTERM);

      //    socket->close(); // this gives exception of type zmq::error_t: Socket operation on non-socket
      pusher.reset(nullptr);
      puller.reset(nullptr);
      zmqSvc().close_context(); // if you don't close context in parent process as well, the next repeat will hang

      // wait for child
      int status = 0;
      pid_t pid;
      do {
         pid = waitpid(child_pid, &status, 0);
      } while (-1 == pid && EINTR == errno); // retry on interrupted system call

      if (0 != status) {
         if (WIFEXITED(status)) {
            printf("exited, status=%d\n", WEXITSTATUS(status));
         } else if (WIFSIGNALED(status)) {
            printf("killed by signal %d\n", WTERMSIG(status));
         } else if (WIFSTOPPED(status)) {
            printf("stopped by signal %d\n", WSTOPSIG(status));
         } else if (WIFCONTINUED(status)) {
            printf("continued\n");
         }
      }

      if (-1 == pid) {
         throw std::runtime_error(std::string("waitpid, errno ") + std::to_string(errno));
      }

      sigprocmask(SIG_SETMASK, &sigmask_old, nullptr);
   } else { // child
            //      std::cout << "child PID: " << getpid() << std::endl;

      sigset_t sigmask, sigmask_old;
      sigemptyset(&sigmask);
      sigaddset(&sigmask, SIGTERM);
      sigprocmask(SIG_BLOCK, &sigmask, &sigmask_old);

      struct sigaction sa;
      memset(&sa, '\0', sizeof(sa));
      sa.sa_handler = handle_sigterm;

      if (sigaction(SIGTERM, &sa, NULL) < 0) {
         std::perror("sigaction failed");
         std::exit(1);
      }

      ZmqLingeringSocketPtr<> puller, pusher;
      puller.reset(zmqSvc().socket_ptr(zmq::PULL));
      puller->connect("ipc:///tmp/ZMQ_test_fork_polling_M2C.ipc");
      pusher.reset(zmqSvc().socket_ptr(zmq::PUSH));
      pusher->connect("ipc:///tmp/ZMQ_test_fork_polling_C2M.ipc");

      ZeroMQPoller poller1, poller2;
      poller1.register_socket(*puller, zmq::POLLIN);
      poller2.register_socket(*puller, zmq::POLLIN);

      // start test
      auto result1a = poller1.poll(-1);
      auto result1b = poller1.poll(-1);
      auto result2 = poller2.poll(-1);
      EXPECT_EQ(result1a.size(), result1b.size());
      EXPECT_EQ(result1a.size(), result2.size());

      auto receipt = zmqSvc().receive<std::string>(*puller, ZMQ_DONTWAIT);
      if (receipt == "breaker breaker") {
         zmqSvc().send(*pusher, 1212);
      }
      // take care, don't just use _exit, it will not cleanly destroy context etc!
      // if you really need to, at least close and destroy everything properly
      //    socket->close(); // this gives exception of type zmq::error_t: Socket operation on non-socket

      sigprocmask(SIG_SETMASK, &sigmask_old, nullptr);

      while (!terminated) {
      }

      //      std::cout << "child terminated" << std::endl;

      puller.reset(nullptr);
      pusher.reset(nullptr);
      zmqSvc().close_context();
      _Exit(0);
   }
}
