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
#ifndef ROOT_ROOFIT_MultiProcess_Messenger
#define ROOT_ROOFIT_MultiProcess_Messenger

#include "RooFit/MultiProcess/Messenger_decl.h"

namespace RooFit {
namespace MultiProcess {

// -- WORKER - QUEUE COMMUNICATION --

template <typename T, typename... Ts>
void Messenger::send_from_worker_to_queue(T item, Ts... items)
{
   #ifndef NDEBUG
   std::stringstream ss;
   ss << "PID " << getpid() << " sends W2Q " << item;
   debug_print(ss.str());
   #endif

   zmqSvc().send(*this_worker_qw_push, item, send_flag);
   //      if (sizeof...(items) > 0) {  // this will only work with if constexpr, c++17
   send_from_worker_to_queue(items...);
}

template <typename value_t>
value_t Messenger::receive_from_worker_on_queue(std::size_t this_worker_id)
{
   qw_pull_poller[this_worker_id].ppoll(-1, &ppoll_sigmask);
   auto value = zmqSvc().receive<value_t>(*qw_pull[this_worker_id], ZMQ_DONTWAIT);

   #ifndef NDEBUG
   std::stringstream ss;
   ss << "PID " << getpid() << " receives W(" << this_worker_id << ")2Q " << value;
   debug_print(ss.str());
   #endif

   return value;
}

template <typename T, typename... Ts>
void Messenger::send_from_queue_to_worker(std::size_t this_worker_id, T item, Ts... items)
{
#ifndef NDEBUG
   std::stringstream ss;
   ss << "PID " << getpid() << " sends Q2W(" << this_worker_id << ") " << item;
   debug_print(ss.str());
#endif

   zmqSvc().send(*qw_push[this_worker_id], item, send_flag);
   //      if (sizeof...(items) > 0) {  // this will only work with if constexpr, c++17
   send_from_queue_to_worker(this_worker_id, items...);
}

template <typename value_t>
value_t Messenger::receive_from_queue_on_worker()
{
   qw_pull_poller[0].ppoll(-1, &ppoll_sigmask);
   auto value = zmqSvc().receive<value_t>(*this_worker_qw_pull, ZMQ_DONTWAIT);

   #ifndef NDEBUG
   std::stringstream ss;
   ss << "PID " << getpid() << " receives Q2W " << value;
   debug_print(ss.str());
   #endif

   return value;
}


// -- QUEUE - MASTER COMMUNICATION --

template <typename T, typename... Ts>
void Messenger::send_from_queue_to_master(T item, Ts... items)
{
   #ifndef NDEBUG
   std::stringstream ss;
   ss << "PID " << getpid() << " sends Q2M " << item;
   debug_print(ss.str());
   #endif

   zmqSvc().send(*mq_push, item, send_flag);
   //      if (sizeof...(items) > 0) {  // this will only work with if constexpr, c++17
   send_from_queue_to_master(items...);
}

template <typename value_t>
value_t Messenger::receive_from_queue_on_master()
{
   mq_pull_poller.ppoll(-1, &ppoll_sigmask);
   auto value = zmqSvc().receive<value_t>(*mq_pull, ZMQ_DONTWAIT);

   #ifndef NDEBUG
   std::stringstream ss;
   ss << "PID " << getpid() << " receives Q2M " << value;
   debug_print(ss.str());
   #endif

   return value;
}

template <typename T, typename... Ts>
void Messenger::send_from_master_to_queue(T item, Ts... items)
{
   #ifndef NDEBUG
   std::stringstream ss;
   ss << "PID " << getpid() << " sends M2Q " << item;
   debug_print(ss.str());
   #endif

   zmqSvc().send(*mq_push, item, send_flag);
   //      if (sizeof...(items) > 0) {  // this will only work with if constexpr, c++17
   send_from_master_to_queue(items...);
}

template <typename value_t>
value_t Messenger::receive_from_master_on_queue()
{
   mq_pull_poller.ppoll(-1, &ppoll_sigmask);
   auto value = zmqSvc().receive<value_t>(*mq_pull, ZMQ_DONTWAIT);

   #ifndef NDEBUG
   std::stringstream ss;
   ss << "PID " << getpid() << " receives M2Q " << value;
   debug_print(ss.str());
   #endif

   return value;
}


// -- MASTER - WORKER COMMUNICATION --

template <typename T, typename... Ts>
void Messenger::publish_from_master_to_workers(T item, Ts... items)
{
#ifndef NDEBUG
   std::stringstream ss;
   ss << "PID " << getpid() << " sends M2W " << item;
   debug_print(ss.str());
#endif

   zmqSvc().send(*mw_pub, item, send_flag);
   //      if (sizeof...(items) > 0) {  // this will only work with if constexpr, c++17
   publish_from_master_to_workers(items...);
}

template <typename value_t>
value_t Messenger::receive_from_master_on_worker()
{
   mw_sub_poller.ppoll(-1, &ppoll_sigmask);
   auto value = zmqSvc().receive<value_t>(*mw_sub, ZMQ_DONTWAIT);

#ifndef NDEBUG
   std::stringstream ss;
   ss << "PID " << getpid() << " receives M2W " << value;
   debug_print(ss.str());
#endif

   return value;
}

template <typename T, typename... Ts>
void Messenger::send_from_worker_to_master(T item, Ts... items)
{
#ifndef NDEBUG
   std::stringstream ss;
   ss << "PID " << getpid() << " sends M2W " << item;
   debug_print(ss.str());
#endif

   zmqSvc().send(*wm_push, item, send_flag);
   //      if (sizeof...(items) > 0) {  // this will only work with if constexpr, c++17
   send_from_worker_to_master(items...);
}

template <typename value_t>
value_t Messenger::receive_from_worker_on_master()
{
   wm_pull_poller.ppoll(-1, &ppoll_sigmask);
   auto value = zmqSvc().receive<value_t>(*wm_pull, ZMQ_DONTWAIT);

#ifndef NDEBUG
   std::stringstream ss;
   ss << "PID " << getpid() << " receives M2W " << value;
   debug_print(ss.str());
#endif

   return value;
}

} // namespace MultiProcess
} // namespace RooFit

#endif // ROOT_ROOFIT_MultiProcess_Messenger
