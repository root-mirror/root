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

#include <RooFit/MultiProcess/ProcessManager.h>  // ... JobManager::process_manager()

#include "gtest/gtest.h"

TEST(TestMPProcessManager, birthAndDeath)
{
   auto pm = RooFit::MultiProcess::ProcessManager(2);
}

// disabled: it is not working anymore (16 March 2021), due to changes introduced to fix Messenger::test_connections
TEST(TestMPProcessManager, DISABLED_multiBirth)
{
   // This test doesn't actually represent a current usecase in RooFit, but let's
   // showcase the possibility anyway for future reference.

   // first create a regular pm like in the birthAndDeath test
   auto pm = RooFit::MultiProcess::ProcessManager(2);
   // then from each forked node, spin up another set of workers+queue!
   auto pm2 = RooFit::MultiProcess::ProcessManager(2);
}


TEST(TestMPProcessManager, checkState)
{
   std::size_t N_workers = 2;
   auto master_pid = getpid();
   auto pm = RooFit::MultiProcess::ProcessManager(N_workers);
   ASSERT_TRUE(pm.is_initialized());
   ASSERT_EQ(pm.N_workers(), N_workers);

   if (pm.is_master()) {
      EXPECT_EQ(getpid(), master_pid);
      EXPECT_FALSE(pm.is_queue());
      EXPECT_FALSE(pm.is_worker());
   } else if (pm.is_queue()) {
      EXPECT_FALSE(pm.is_master());
      EXPECT_FALSE(pm.is_worker());
   } else if (pm.is_worker()) {
      EXPECT_FALSE(pm.is_master());
      EXPECT_FALSE(pm.is_queue());
      EXPECT_TRUE(pm.worker_id() == 0 || pm.worker_id() == 1);
   }
}