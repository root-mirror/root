#include <zmq.hpp>
#include "RooFit_ZMQ/ZeroMQSvc.h"
#include "RooFit_ZMQ/functions.h"

namespace ZMQ {
   size_t stringLength(const char& cs) {
      return strlen(&cs);
   }
}
