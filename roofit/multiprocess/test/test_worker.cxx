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

#include "RooFit/MultiProcess/worker.h"

#include "gtest/gtest.h"

// TODO: implement something along these lines (the following is pseudocode)
//TEST(worker, loop)
//{
//   do_fork();
//   if (on_master()) {
//      wait();
//      terminate_worker();
//   } else if (on_worker()) {
//      worker_loop();
//   }
//}
