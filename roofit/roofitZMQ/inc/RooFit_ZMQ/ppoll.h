#ifndef ROOT_ROOFIT_ZMQ_ppoll
#define ROOT_ROOFIT_ZMQ_ppoll

#include <vector>
#include <zmq.hpp>

namespace ZMQ {

int zmq_ppoll(zmq_pollitem_t *items_, int nitems_, long timeout_, const sigset_t *sigmask_);
int ppoll(zmq_pollitem_t *items_, size_t nitems_, long timeout_, const sigset_t *sigmask_);
int ppoll(std::vector<zmq_pollitem_t> &items, long timeout_, const sigset_t *sigmask_);
class ppoll_error_t : public zmq::error_t {
};

} // namespace ZMQ

#endif // ROOT_ROOFIT_ZMQ_ppoll
