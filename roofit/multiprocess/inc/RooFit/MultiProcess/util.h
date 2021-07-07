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

#ifndef ROOT_ROOFIT_MultiProcess_util
#define ROOT_ROOFIT_MultiProcess_util

#include <unistd.h> // getpid, pid_t

#include <RooFit_ZMQ/ppoll.h> // for ZMQ::ppoll_error_t
#include <RooFit_ZMQ/ZeroMQPoller.h>

namespace RooFit {
namespace MultiProcess {

int wait_for_child(pid_t child_pid, bool may_throw, int retries_before_killing);

enum class zmq_ppoll_error_response { abort, unknown_eintr, retry };
zmq_ppoll_error_response handle_zmq_ppoll_error(ZMQ::ppoll_error_t &e);
std::tuple<std::vector<std::pair<size_t, int>>, bool> careful_ppoll(ZeroMQPoller &poller, const sigset_t &ppoll_sigmask, std::size_t max_tries = 2);

} // namespace MultiProcess
} // namespace RooFit
#endif // ROOT_ROOFIT_MultiProcess_util
