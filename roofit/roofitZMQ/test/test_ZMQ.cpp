/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id$
 * Authors:                                                                  *
 *   PB, Patrick Bos,     NL eScience Center, p.bos@esciencecenter.nl        *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#include "gtest/gtest.h"

#include <unistd.h> // fork, usleep

#include <sstream>

#include "RooFit_ZMQ/ZeroMQSvc.h"


// N.B.: wait_for_child is identically defined in RooFit::MultiProcess, but we copy it here to reduce module
// interdependencies. It also requires an extra include:
#include <csignal>    // kill, SIGKILL

int wait_for_child(pid_t child_pid, bool may_throw, int retries_before_killing)
{
   int status = 0;
   int patience = retries_before_killing;
   pid_t tmp;
   do {
      if (patience-- < 1) {
         ::kill(child_pid, SIGKILL);
      }
      tmp = waitpid(child_pid, &status, WNOHANG);
   } while (tmp == 0                         // child has not yet changed state, try again
            || (-1 == tmp && EINTR == errno) // retry on interrupted system call
      );

   if (patience < 1) {
      std::cout << "Had to send PID " << child_pid << " " << (-patience + 1) << " SIGKILLs";
   }

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

   if (-1 == tmp && may_throw)
      throw std::runtime_error(std::string("waitpid, errno ") + std::to_string(errno));

   return status;
}

void elaborate_bind(const ZmqLingeringSocketPtr<>& socket, std::string name) {
  try {
    socket->bind(name);
  } catch (const zmq::error_t& e) {
    if (e.num() == EADDRINUSE) {
      std::cerr << "address already in use, retrying bind in 500ms\n";
      usleep(500000);
      try {
        socket->bind(name);
      } catch (const zmq::error_t& e2) {
        if (e2.num() == EADDRINUSE) {
          std::cerr
              << "again: address already in use, aborting; please check whether there are any remaining improperly exited processes (zombies) around or whether some other program is using port 6660\n";
        }
        throw e2;
      }
      // Sometimes, the socket from the previous test needs some time to close, so
      // we introduce a latency here. A more robust and fast approach might be to
      // do the following on the bind side:
      // 1. first try another port, e.g. increase by one
      // 2. if that doesn't work, do the latency and retry the original port
      // The connect side then also needs to change, because it doesn't know which
      // port the bind side will bind to. The connect side could try connecting to
      // both options asynchronously, and then in a loop check both for signs of
      // life. If one comes alive, transfer ownership of that pointer to the pointer
      // you want to eventually use (`socket`) and that's it.
    } else {
      throw e;
    }
  }
}

class AllSocketTypes : public ::testing::TestWithParam< std::tuple<int, std::pair<zmq::SocketTypes, zmq::SocketTypes>, std::pair<std::string, std::string> /* socket_names */> > {};

TEST_P(AllSocketTypes, forkHandshake) {
  auto socket_names = std::get<2>(GetParam());
  pid_t child_pid {0};
  do {
    child_pid = fork();
  } while (child_pid == -1);  // retry if fork fails

  if (child_pid > 0) {                // master
    ZmqLingeringSocketPtr<> socket;
    socket.reset(zmqSvc().socket_ptr(std::get<1>(GetParam()).first));
    elaborate_bind(socket, socket_names.second);
    // bind is on the master process to avoid zombie children to hold on to binds

    // start test
    zmqSvc().send(*socket, std::string("breaker breaker"));

    auto receipt = zmqSvc().receive<int>(*socket);

    EXPECT_EQ(receipt, 1212);

//    socket->close(); // this gives exception of type zmq::error_t: Socket operation on non-socket
    socket.reset(nullptr);
    zmqSvc().close_context(); // if you don't close context in parent process as well, the next repeat will hang

    wait_for_child(child_pid, true, 5);
  } else {                            // child
    ZmqLingeringSocketPtr<> socket;
    socket.reset(zmqSvc().socket_ptr(std::get<1>(GetParam()).second));
    socket->connect(socket_names.first);

    // start test
    auto receipt = zmqSvc().receive<std::string>(*socket);
    if (receipt == "breaker breaker") {
      zmqSvc().send(*socket, 1212);
    }
    // take care, don't just use _exit, it will not cleanly destroy context etc!
    // if you really need to, at least close and destroy everything properly
//    socket->close(); // this gives exception of type zmq::error_t: Socket operation on non-socket
    socket.reset(nullptr);
    zmqSvc().close_context();
    _Exit(0);
  }
}

std::string ipc {"ipc:///tmp/ZMQ_test_fork.ipc"};
std::string tcp_server {"tcp://127.0.0.1:6660"};
std::string tcp_client {"tcp://*:6660"};
auto socket_name_options = ::testing::Values(
    std::make_pair(tcp_server, tcp_client),
    std::make_pair(ipc, ipc)
    );


INSTANTIATE_TEST_SUITE_P(REQREP, AllSocketTypes,
                        ::testing::Combine(::testing::Range(0, 10), // repeat to probe connection stability
                                           ::testing::Values(std::make_pair(zmq::REQ, zmq::REP)),
                                           socket_name_options
                        ));
INSTANTIATE_TEST_SUITE_P(PAIRPAIR, AllSocketTypes,
                        ::testing::Combine(::testing::Range(0, 10), // repeat to probe connection stability
                                           ::testing::Values(std::make_pair(zmq::PAIR, zmq::PAIR)),
                                           socket_name_options
                        ));


class AsyncSocketTypes : public ::testing::TestWithParam< std::tuple<int /* repeat_nr */, std::pair<zmq::SocketTypes, zmq::SocketTypes>, std::pair<std::string, std::string> /* socket_names */, bool /* expect_throw */> > {};

TEST_P(AsyncSocketTypes, forkMultiSendReceive) {
  bool expect_throw = std::get<3>(GetParam());
  ZmqLingeringSocketPtr<> socket;
  auto socket_names = std::get<2>(GetParam());
  pid_t child_pid {0};
  do {
    child_pid = fork();
  } while (child_pid == -1);  // retry if fork fails

  if (child_pid > 0) {                // master
    socket.reset(zmqSvc().socket_ptr(std::get<1>(GetParam()).first));
    elaborate_bind(socket, socket_names.second);
    // bind is on the master process to avoid zombie children to hold on to binds

    // start test: send 2 things, receive 1, send 1 more, finish
    zmqSvc().send(*socket, std::string("breaker breaker"));

    if (expect_throw) {
      EXPECT_ANY_THROW(zmqSvc().send(*socket, std::string("anybody out there?")));
      // NOTE: also in case of a throw, be sure to properly close down the connection!
      // Otherwise, you may get zombies waiting for a reply.
      socket.reset(nullptr);
      zmqSvc().close_context(); // if you don't close context in parent process as well, the next repeat will hang
      wait_for_child(child_pid, true, 5);
      return;
    } else {
      EXPECT_NO_THROW(zmqSvc().send(*socket, std::string("anybody out there?")));
    }

    auto receipt = zmqSvc().receive<int>(*socket);

    EXPECT_EQ(receipt, 1212);

    zmqSvc().send(*socket, std::string("kthxbye"));

//    socket->close(); // this gives exception of type zmq::error_t: Socket operation on non-socket
    socket.reset(nullptr);
    zmqSvc().close_context(); // if you don't close context in parent process as well, the next repeat will hang

    wait_for_child(child_pid, true, 5);
  } else {                            // child
    socket.reset(zmqSvc().socket_ptr(std::get<1>(GetParam()).second));
    socket->connect(socket_names.first);

    // start test, receive something
    auto receipt1 = zmqSvc().receive<std::string>(*socket);
    auto receipt2 = zmqSvc().receive<std::string>(*socket);
    if (receipt1 == "breaker breaker" && receipt2 == "anybody out there?") {
      zmqSvc().send(*socket, 1212);
    }
    auto receipt3 = zmqSvc().receive<std::string>(*socket);
    if (receipt3 != "kthxbye") {
      std::cerr << "did not receive final reply correctly\n";
    }

    // take care, don't just use _exit, it will not cleanly destroy context etc!
    // if you really need to, at least close and destroy everything properly
//    socket->close(); // this gives exception of type zmq::error_t: Socket operation on non-socket
    socket.reset(nullptr);
    zmqSvc().close_context();
    _Exit(0);
  }
}

TEST_P(AsyncSocketTypes, forkIgnoreSomeMessages) {
  bool expect_throw = std::get<3>(GetParam());
  ZmqLingeringSocketPtr<> socket;
  auto socket_names = std::get<2>(GetParam());
  pid_t child_pid {0};
  do {
    child_pid = fork();
  } while (child_pid == -1);  // retry if fork fails

  if (child_pid > 0) {                // master
    socket.reset(zmqSvc().socket_ptr(std::get<1>(GetParam()).first));
    elaborate_bind(socket, socket_names.second);
    // bind is on the master process to avoid zombie children to hold on to binds

    // start test: send 2 things, receive 1, send 1 more, finish
    zmqSvc().send(*socket, std::string("breaker breaker"));

    if (expect_throw) {
      EXPECT_ANY_THROW(zmqSvc().send(*socket, std::string("anybody out there?")));
      // NOTE: also in case of a throw, be sure to properly close down the connection!
      // Otherwise, you may get zombies waiting for a reply.
      socket.reset(nullptr);
      zmqSvc().close_context(); // if you don't close context in parent process as well, the next repeat will hang
      wait_for_child(child_pid, true, 5);
      return;
    } else {
      EXPECT_NO_THROW(zmqSvc().send(*socket, std::string("anybody out there?")));
    }

    auto receipt = zmqSvc().receive<int>(*socket);

    EXPECT_EQ(receipt, 1212);

    zmqSvc().send(*socket, std::string("kthxbye"));

//    socket->close(); // this gives exception of type zmq::error_t: Socket operation on non-socket
    socket.reset(nullptr);
    zmqSvc().close_context(); // if you don't close context in parent process as well, the next repeat will hang

    wait_for_child(child_pid, true, 5);
  } else {                            // child
    socket.reset(zmqSvc().socket_ptr(std::get<1>(GetParam()).second));
    socket->connect(socket_names.first);

    // start test, receive first thing
    auto receipt = zmqSvc().receive<std::string>(*socket);
    if (receipt == "breaker breaker") {
      zmqSvc().send(*socket, 1212);
    }

    // ignore the rest of the sent messages, but give the other end a second to
    // actually send its stuff, instead of hanging in retry_send because the
    // connection has died; a better solution would be if retry_send (in
    // ZeroMQSvc::send) had a callback mechanism that could be used to break
    // out when a child has died, but ok
    sleep(1);

    // take care, don't just use _exit, it will not cleanly destroy context etc!
    // if you really need to, at least close and destroy everything properly
//    socket->close(); // this gives exception of type zmq::error_t: Socket operation on non-socket
    socket.reset(nullptr);
    zmqSvc().close_context();
    _Exit(0);
  }
}


INSTANTIATE_TEST_SUITE_P(PAIRPAIR, AsyncSocketTypes,
                        ::testing::Combine(::testing::Range(0, 10), // repeat to probe connection stability
                                           ::testing::Values(std::make_pair(zmq::PAIR, zmq::PAIR)),
                                           socket_name_options,
                                           ::testing::Values(false) // don't expect throw
                        ));

INSTANTIATE_TEST_SUITE_P(REQREP, AsyncSocketTypes,
                        ::testing::Combine(::testing::Values(0), // no repeats, we only care about the throw
                                           ::testing::Values(std::make_pair(zmq::REQ, zmq::REP)),
                                           socket_name_options,
                                           ::testing::Values(true) // expect throw
                        ));
