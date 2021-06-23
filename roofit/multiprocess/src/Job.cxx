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

#include "RooFit/MultiProcess/JobManager.h"
#include "RooFit/MultiProcess/Messenger.h"
#include "RooFit/MultiProcess/worker.h"
#include "RooFit/MultiProcess/Job.h"
#include "RooFit/MultiProcess/Queue.h"

namespace RooFit {
namespace MultiProcess {

Job::Job()
{
   id = JobManager::add_job_object(this);
}

Job::Job(const Job &other)
   : _manager(other._manager)
{
   id = JobManager::add_job_object(this);
}

Job::~Job()
{
   JobManager::remove_job_object(id);
}

// This function is necessary here, because the Job knows about the number
// of workers, so only from the Job can the JobManager be instantiated.
JobManager *Job::get_manager()
{
   if (!_manager) {
      _manager = JobManager::instance();
   }

   if (!_manager->is_activated()) {
      _manager->activate();
   }

   return _manager;
}

void Job::gather_worker_results()
{
   get_manager()->retrieve(id);
}

void Job::update_state() {}


} // namespace MultiProcess
} // namespace RooFit
