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

#include "RooFit/MultiProcess/Messenger.h"
#include <RooFit/MultiProcess/ProcessManager.h>  // ... JobManager::process_manager()

#include "gtest/gtest.h"

TEST(TestMPMessenger, Connections)
{
   RooFit::MultiProcess::ProcessManager pm(2);
   RooFit::MultiProcess::Messenger messenger(pm);
   messenger.test_connections(pm);
}

TEST(TestMPMessenger, ConnectionsManualExit)
{
   // the point of this test is to see whether clean-up of ZeroMQ resources is done properly without calling any
   // destructors (which is what happens when you call _Exit() instead of regularly ending the program by reaching the
   // end of main()).
   RooFit::MultiProcess::ProcessManager pm( 2);
   RooFit::MultiProcess::Messenger messenger(pm);
   messenger.test_connections(pm);
   if (!pm.is_master()) {
      std::_Exit(0);
   }
}


#include <sys/wait.h>

TEST(TestMPMessenger, SigStop)
{
   // For some reason, I cannot get my debugger (lldb) to play nicely with ZeroMQ; every time I pause a process, it
   // says it crashed because of SIGSTOP, but SIGSTOP should just pause it... let's see whether we can reproduce this
   // behavior here and, if possible, fix it.
   RooFit::MultiProcess::ProcessManager pm(2);
   RooFit::MultiProcess::Messenger messenger(pm);
//   messenger.test_connections(pm);

   // all this works fine without a sigmask, maybe the sigmask is the issue? let's try:


   if (pm.is_master()) {
      printf("first send message to queue...\n");
      messenger.send_from_master_to_queue(1);
      printf("sleep for 2 seconds...\n");
      sleep(2);
      printf("SIGSTOPping all children...\n");
      kill(pm.get_queue_pid(), SIGSTOP);
      kill(pm.get_worker_pids()[0], SIGSTOP);
      kill(pm.get_worker_pids()[1], SIGSTOP);
      printf("send another message to queue...\n");
      messenger.send_from_master_to_queue(2);
      printf("sleep for 2 seconds...\n");
      sleep(2);
      printf("SIGCONT queue and worker 1...\n");
      kill(pm.get_queue_pid(), SIGCONT);
      kill(pm.get_worker_pids()[0], SIGCONT);
      printf("sleep for 2 seconds...\n");
      sleep(2);
      printf("SIGCONT worker 2...\n");
      kill(pm.get_worker_pids()[1], SIGCONT);
      EXPECT_EQ(messenger.receive_from_queue_on_master<int>(), 3);
   }

   if (pm.is_queue()) {
      auto number = messenger.receive_from_master_on_queue<int>();
      printf("received %d on queue\n", number);
      number = messenger.receive_from_master_on_queue<int>();
      printf("received %d on queue\n", number);
      messenger.send_from_queue_to_master(3);
   }

   if (!pm.is_master()) {
      while (!pm.sigterm_received()) {}
   }
}

TEST(TestMPMessenger, StressPubSub)
{
   std::size_t N_workers = 64;
   RooFit::MultiProcess::ProcessManager pm(N_workers);
   RooFit::MultiProcess::Messenger messenger(pm);

   if (pm.is_master()) {
      messenger.publish_from_master_to_workers("bert");
      std::size_t ernies = 0;
      for (; ernies < N_workers; ++ernies) {
         auto receipt = messenger.receive_from_worker_on_master<std::string>();
         if (receipt != "ernie") {
            printf("whoops, got %s instead of ernie!\n", receipt.c_str());
            FAIL();
         }
      }
      EXPECT_EQ(ernies, N_workers);
   } else if (pm.is_worker()) {
      auto receipt = messenger.receive_from_master_on_worker<std::string>();
      if (receipt == "bert") {
         messenger.send_from_worker_to_master("ernie");
      } else {
         printf("no bert on worker %lu\n", pm.worker_id());
      }
   }
}
