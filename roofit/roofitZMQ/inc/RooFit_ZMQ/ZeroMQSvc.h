#ifndef ZEROMQ_IZEROMQSVC_H
#define ZEROMQ_IZEROMQSVC_H 1

// Include files
// from STL
#include <type_traits>
#include <string>
#include <vector>
#include <sstream>
#include <ios>
#include <iostream>  // std::cerr

// ZeroMQ
#include <zmq.hpp>
#include "RooFit_ZMQ/Utility.h"
#include "RooFit_ZMQ/functions.h"

// debugging
#include <unistd.h> // getpid

namespace ZMQ {

struct TimeOutException : std::exception {
   TimeOutException() = default;
   TimeOutException(const TimeOutException&) = default;
   ~TimeOutException() = default;
   TimeOutException& operator=(const TimeOutException&) = default;
};

struct MoreException : std::exception {
   MoreException() = default;
   MoreException(const MoreException&) = default;
   ~MoreException() = default;
   MoreException& operator=(const MoreException&) = default;
};

}

template <int PERIOD = 0>
struct ZmqLingeringSocketPtrDeleter {
   void operator()(zmq::socket_t* socket) {
      auto period = PERIOD;

      int tries = 0;
      int max_tries = 3;
      while(true) {
       try {
         // the actual work this function should do, plus the delete socket below:
         if (socket) socket->setsockopt(ZMQ_LINGER, &period, sizeof(period));
         break;
       } catch (zmq::error_t& e) {
         if (++tries == max_tries
             || e.num() == EINVAL || e.num() == ETERM || e.num() == ENOTSOCK  // not recoverable from here
             ) {
           std::cerr << "ERROR in ZeroMQSvc::socket: " << e.what() << " (errno: " << e.num() << ")\n";
           throw;
         }
         std::cerr << "RETRY " << tries << "/" << (max_tries - 1) << " in ZmqLingeringSocketPtrDeleter: call interrupted (errno: " << e.num() << ")\n";
       }
      }

      delete socket;
   }
};

template <int PERIOD = 0>
using ZmqLingeringSocketPtr = std::unique_ptr<zmq::socket_t, ZmqLingeringSocketPtrDeleter<PERIOD>>;


// We retry send and receive only on EINTR, all other errors are either fatal, or can only
// be handled at the caller.
template <typename... args_t>
auto retry_send(zmq::socket_t& socket, int max_tries, args_t ...args) -> decltype(socket.send(args...)) {
  int tries = 0;
  while(true) {
    try {
      // the actual work this function should do, all the rest is error handling:
      return socket.send(args...);
    } catch (zmq::error_t& e) {
      if (++tries == max_tries
          || e.num() != EINTR  // only recoverable error
          ) {
//        std::cerr << "ERROR in ZeroMQSvc::send (retry_send) on pid " << getpid() << ": " << e.what() << " (errno: " << e.num() << ")\n";
        throw;
      }
      std::cerr << "RETRY " << tries << "/" << (max_tries - 1) << " in ZeroMQSvc::send (retry_send) on pid " << getpid() << ": " << e.what() << ")\n";
    }
  }
}

template <typename... args_t>
auto retry_recv(zmq::socket_t& socket, int max_tries, args_t ...args) -> decltype(socket.recv(args...)) {
  int tries = 0;
  while(true) {
    try {
      // the actual work this function should do, all the rest is error handling:
      return socket.recv(args...);
    } catch (zmq::error_t& e) {
      if (++tries == max_tries
          || e.num() != EINTR  // only recoverable error
          ) {
//        std::cerr << "ERROR in ZeroMQSvc::recv (retry_recv) on pid " << getpid() << ": " << e.what() << " (errno: " << e.num() << ")\n";
        throw;
      }
      std::cerr << "RETRY " << tries << "/" << (max_tries - 1) << " in ZeroMQSvc::recv (retry_recv) on pid " << getpid() << ": " << e.what() << ")\n";
    }
  }
}


/** @class IZeroMQSvc IZeroMQSvc.h ZeroMQ/IZeroMQSvc.h
 *
 *
 *  @author
 *  @date   2015-06-22
 */
class ZeroMQSvc {
  // Note on error handling:
  // Creating message_t can throw, but only when memory ran out (errno ENOMEM),
  // and that is something only the caller can fix, so we don't catch it here.

public:

   enum Encoding {
      Text = 0,
      Binary
   };

   Encoding encoding() const;
   void setEncoding(const Encoding& e);
   zmq::context_t& context() const;
   zmq::socket_t socket(int type) const;
   zmq::socket_t* socket_ptr(int type) const;
   void close_context() const;

  // decode message with ZMQ, POD version
   template <class T, typename std::enable_if<!std::is_pointer<T>::value
                                              && ZMQ::Detail::is_trivial<T>::value, T>::type* = nullptr>
   T decode(const zmq::message_t& msg) const {
      T object;
      memcpy(&object, msg.data(), msg.size());
      return object;
   }

   // decode ZMQ message, string version
   template <class T, typename std::enable_if<std::is_same<T, std::string>::value, T>::type* = nullptr>
   std::string decode(const zmq::message_t& msg) const {
      std::string r(msg.size() + 1, char{});
      r.assign(static_cast<const char*>(msg.data()), msg.size());
      return r;
   }

   // receive message with ZMQ, general version
   // FIXME: what to do with flags=0.... more is a pointer, that might prevent conversion
   template <class T, typename std::enable_if<!(std::is_same<zmq::message_t, T>::value), T>::type* = nullptr>
   T receive(zmq::socket_t& socket, int flags = 0, bool* more = nullptr) const {
      // receive message
      zmq::message_t msg;
      auto nbytes = retry_recv(socket, 2, &msg, flags);
      if (0 == nbytes) {
         throw ZMQ::TimeOutException{};
      }
      if (more) *more = msg.more();

      // decode message
      return decode<T>(msg);
   }

   // receive message with ZMQ
   template <class T, typename std::enable_if<std::is_same<zmq::message_t, T>::value, T>::type* = nullptr>
   T receive(zmq::socket_t& socket, int flags = 0, bool* more = nullptr) const {
      // receive message
      zmq::message_t msg;
      auto nbytes = retry_recv(socket, 2, &msg, flags);
      if (0 == nbytes) {
         throw ZMQ::TimeOutException{};
      }
      if (more) *more = msg.more();
      return msg;
   }

   // encode message to ZMQ
   template <class T, typename std::enable_if<!std::is_pointer<T>::value
                                              && ZMQ::Detail::is_trivial<T>::value, T>::type* = nullptr>
   zmq::message_t encode(const T& item, std::function<size_t(const T& t)> sizeFun = ZMQ::defaultSizeOf<T>) const {
      size_t s = sizeFun(item);
      zmq::message_t msg{s};
      memcpy((void *)msg.data(), &item, s);
      return msg;
   }

   zmq::message_t encode(const char* item) const;
   zmq::message_t encode(const std::string& item) const;

   // Send message with ZMQ
   template <class T, typename std::enable_if<!std::is_same<T, zmq::message_t>::value, T>::type* = nullptr>
      bool send(zmq::socket_t& socket, const T& item, int flags = 0) const {
      return retry_send(socket, 1, encode(item), flags);
   }

   bool send(zmq::socket_t& socket, const char* item, int flags = 0) const;
   bool send(zmq::socket_t& socket, zmq::message_t& msg, int flags = 0) const;
   bool send(zmq::socket_t& socket, zmq::message_t&& msg, int flags = 0) const;

private:

   Encoding m_enc = Text;
   mutable zmq::context_t* m_context = nullptr;

};

ZeroMQSvc& zmqSvc();

#endif // ZEROMQ_IZEROMQSVC_H
