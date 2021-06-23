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
#ifndef ROOT_ROOFIT_ZMQ_ppoll
#define ROOT_ROOFIT_ZMQ_ppoll

#include <vector>
#include "zmq.hxx"

namespace ZMQ {

int zmq_ppoll (zmq_pollitem_t *items_, int nitems_, long timeout_, const sigset_t * sigmask_);
int ppoll(zmq_pollitem_t *items_, size_t nitems_, long timeout_, const sigset_t * sigmask_);
int ppoll(std::vector<zmq_pollitem_t> &items, long timeout_, const sigset_t * sigmask_);
class ppoll_error_t : public zmq::error_t {};

}

#endif // ROOT_ROOFIT_ZMQ_ppoll
