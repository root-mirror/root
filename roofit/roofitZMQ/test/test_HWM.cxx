/*
 * Project: RooFit
 * Authors:
 *   PB, Patrick Bos, Netherlands eScience Center, p.bos@esciencecenter.nl
 *
 * Copyright (c) 2016-2021, Netherlands eScience Center
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 */

#include "gtest/gtest.h"

#include <unistd.h> // fork, usleep

#include "RooFit_ZMQ/ZeroMQSvc.h"

class HighWaterMarkTest : public ::testing::Test {
protected:
   void SetUp() override {
      do {
         child_pid = fork();
      } while (child_pid == -1);  // retry if fork fails

      if (child_pid > 0) {                // parent
         pusher.reset(zmqSvc().socket_ptr(zmq::PUSH));
         if (set_hwm) {
            auto rc = zmq_setsockopt(*pusher, ZMQ_SNDHWM, &hwm, sizeof hwm);
            assert(rc == 0);
         }
         pusher->bind("ipc:///tmp/ZMQ_test_fork_polling_P2C.ipc");
      } else {                            // child
         puller.reset(zmqSvc().socket_ptr(zmq::PULL));
         if (set_hwm) {
            auto rc = zmq_setsockopt(*puller, ZMQ_RCVHWM, &hwm, sizeof hwm);
            assert(rc == 0);
         }
         puller->connect("ipc:///tmp/ZMQ_test_fork_polling_P2C.ipc");
      }
   }

   void TearDown() override {
      if (child_pid > 0) {  // parent
         // wait for child
         int status = 0;
         pid_t pid;
         do {
            pid = waitpid(child_pid, &status, 0);
         } while (-1 == pid && EINTR == errno); // retry on interrupted system call
         if (0 != status) {
            if (WIFEXITED(status)) { printf("exited, status=%d\n", WEXITSTATUS(status)); } else if (WIFSIGNALED(status)) { printf("killed by signal %d\n", WTERMSIG(status)); } else if (WIFSTOPPED(status)) { printf("stopped by signal %d\n", WSTOPSIG(status)); } else if (WIFCONTINUED(status)) { printf("continued\n"); }
         }
         if (-1 == pid) {
            throw std::runtime_error(std::string("waitpid, errno ") + std::to_string(errno));
         }
         pusher.reset(nullptr);
         zmqSvc().close_context();
      } else {              // child
         puller.reset(nullptr);
         zmqSvc().close_context();
      }
   }

   void run_test() {
      std::size_t max_sends = 2000;

      if (child_pid > 0) {                // parent
         // start test
         for (std::size_t ix = 0; ix < max_sends; ++ix) {
            zmqSvc().send(*pusher, 0.1f);
            zmqSvc().send(*pusher, 1);
            zmqSvc().send(*pusher, true);
            if (ix % 100 == 0) {
               printf("parent at ix = %lu\n", ix);
            }
         }
      } else {                            // child
         // wait a few seconds before reading to allow the parent to overflow the HWM
         printf("child waiting for 2 seconds...\n");
         sleep(2);
         printf("child starts receiving\n");

         for (std::size_t ix = 0; ix < max_sends; ++ix) {
            zmqSvc().receive<float>(*puller);
            zmqSvc().receive<int>(*puller);
            zmqSvc().receive<bool>(*puller);
            if (ix % 100 == 0) {
               printf("child at ix = %lu\n", ix);
            }
         }
      }
   }

   pid_t child_pid {0};
   ZmqLingeringSocketPtr<> pusher, puller;
   bool set_hwm = false;
   int hwm = 0;
};


TEST_F(HighWaterMarkTest, demonstrateHittingDefaultHWM) {
   run_test();
}

class HighWaterMarkZeroTest : public HighWaterMarkTest {
   void SetUp() override {
      set_hwm = true;
      hwm = 0;
      HighWaterMarkTest::SetUp();
   }
};

TEST_F(HighWaterMarkZeroTest, zeroHWM) {
   run_test();
}

class HighWaterMark2kTest : public HighWaterMarkTest {
   void SetUp() override {
      set_hwm = true;
      hwm = 2000;
      HighWaterMarkTest::SetUp();
   }
};

TEST_F(HighWaterMark2kTest, HWM2k) {
   run_test();
}
