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

#include <unistd.h>  // pselect
#include <time.h>  // timespec
#include <poll.h>

#include "likely.hpp"  // unlikely
#include "clock.hpp"   // zmq::clock_t
#include "err.hpp"  // zmq_assert
#include "polling_util.hpp"  // zmq::optimized_fd_set_t, zmq::compute_timeout, zmq::timeout_t
#include "fd.hpp"  // zmq::fd_t

#include <RooFit_ZMQ/ppoll.h>


namespace ZMQ {

int zmq_ppoll (zmq_pollitem_t *items_, int nitems_, long timeout_,
               const sigset_t * sigmask_)
{
#if defined ZMQ_HAVE_POLLER
   // if poller is present, use that if there is at least 1 thread-safe socket,
    // otherwise fall back to the previous implementation as it's faster.
    for (int i = 0; i != nitems_; i++) {
        if (items_[i].socket
            && as_socket_base_t (items_[i].socket)->is_thread_safe ()) {
            return zmq_poller_poll (items_, nitems_, timeout_);
        }
    }
#endif // ZMQ_HAVE_POLLER
#if defined ZMQ_POLL_BASED_ON_POLL || defined ZMQ_POLL_BASED_ON_SELECT
   if (unlikely (nitems_ < 0)) {
        errno = EINVAL;
        return -1;
    }
    if (unlikely (nitems_ == 0)) {
        if (timeout_ == 0)
            return 0;
#if defined ZMQ_HAVE_WINDOWS
        Sleep (timeout_ > 0 ? timeout_ : INFINITE);
        return 0;
#elif defined ZMQ_HAVE_VXWORKS
        struct timespec ns_;
        ns_.tv_sec = timeout_ / 1000;
        ns_.tv_nsec = timeout_ % 1000 * 1000000;
        return nanosleep (&ns_, 0);
#else
        return usleep (timeout_ * 1000);
#endif
    }
    if (!items_) {
        errno = EFAULT;
        return -1;
    }

    zmq::clock_t clock;
    uint64_t now = 0;
    uint64_t end = 0;
#if defined ZMQ_POLL_BASED_ON_POLL
    zmq::fast_vector_t<pollfd, ZMQ_POLLITEMS_DFLT> pollfds (nitems_);

    //  Build pollset for poll () system call.
    for (int i = 0; i != nitems_; i++) {
        //  If the poll item is a 0MQ socket, we poll on the file descriptor
        //  retrieved by the ZMQ_FD socket option.
        if (items_[i].socket) {
            size_t zmq_fd_size = sizeof (zmq::fd_t);
            if (zmq_getsockopt (items_[i].socket, ZMQ_FD, &pollfds[i].fd,
                                &zmq_fd_size)
                == -1) {
                return -1;
            }
            pollfds[i].events = items_[i].events ? POLLIN : 0;
        }
        //  Else, the poll item is a raw file descriptor. Just convert the
        //  events to normal POLLIN/POLLOUT for poll ().
        else {
            pollfds[i].fd = items_[i].fd;
            pollfds[i].events =
              (items_[i].events & ZMQ_POLLIN ? POLLIN : 0)
              | (items_[i].events & ZMQ_POLLOUT ? POLLOUT : 0)
              | (items_[i].events & ZMQ_POLLPRI ? POLLPRI : 0);
        }
    }
#else
    //  Ensure we do not attempt to select () on more than FD_SETSIZE
    //  file descriptors.
    //  TODO since this function is called by a client, we could return errno EINVAL/ENOMEM/... here
    zmq_assert (nitems_ <= FD_SETSIZE);

    zmq::optimized_fd_set_t pollset_in (nitems_);
    FD_ZERO (pollset_in.get ());
    zmq::optimized_fd_set_t pollset_out (nitems_);
    FD_ZERO (pollset_out.get ());
    zmq::optimized_fd_set_t pollset_err (nitems_);
    FD_ZERO (pollset_err.get ());

    zmq::fd_t maxfd = 0;

    //  Build the fd_sets for passing to select ().
    for (int i = 0; i != nitems_; i++) {
        //  If the poll item is a 0MQ socket we are interested in input on the
        //  notification file descriptor retrieved by the ZMQ_FD socket option.
        if (items_[i].socket) {
            size_t zmq_fd_size = sizeof (zmq::fd_t);
            zmq::fd_t notify_fd;
            if (zmq_getsockopt (items_[i].socket, ZMQ_FD, &notify_fd,
                                &zmq_fd_size)
                == -1)
                return -1;
            if (items_[i].events) {
                FD_SET (notify_fd, pollset_in.get ());
                if (maxfd < notify_fd)
                    maxfd = notify_fd;
            }
        }
        //  Else, the poll item is a raw file descriptor. Convert the poll item
        //  events to the appropriate fd_sets.
        else {
            if (items_[i].events & ZMQ_POLLIN)
                FD_SET (items_[i].fd, pollset_in.get ());
            if (items_[i].events & ZMQ_POLLOUT)
                FD_SET (items_[i].fd, pollset_out.get ());
            if (items_[i].events & ZMQ_POLLERR)
                FD_SET (items_[i].fd, pollset_err.get ());
            if (maxfd < items_[i].fd)
                maxfd = items_[i].fd;
        }
    }

    zmq::optimized_fd_set_t inset (nitems_);
    zmq::optimized_fd_set_t outset (nitems_);
    zmq::optimized_fd_set_t errset (nitems_);
#endif

    bool first_pass = true;
    int nevents = 0;

    while (true) {
#if defined ZMQ_POLL_BASED_ON_POLL

        //  Compute the timeout for the subsequent poll.
        zmq::timeout_t timeout =
          zmq::compute_timeout (first_pass, timeout_, now, end);

        //  Wait for events.
        {
            int rc = ppoll (&pollfds[0], nitems_, timeout);
            if (rc == -1 && errno == EINTR) {
                return -1;
            }
            errno_assert (rc >= 0);
        }
        //  Check for the events.
        for (int i = 0; i != nitems_; i++) {
            items_[i].revents = 0;

            //  The poll item is a 0MQ socket. Retrieve pending events
            //  using the ZMQ_EVENTS socket option.
            if (items_[i].socket) {
                size_t zmq_events_size = sizeof (uint32_t);
                uint32_t zmq_events;
                if (zmq_getsockopt (items_[i].socket, ZMQ_EVENTS, &zmq_events,
                                    &zmq_events_size)
                    == -1) {
                    return -1;
                }
                if ((items_[i].events & ZMQ_POLLOUT)
                    && (zmq_events & ZMQ_POLLOUT))
                    items_[i].revents |= ZMQ_POLLOUT;
                if ((items_[i].events & ZMQ_POLLIN)
                    && (zmq_events & ZMQ_POLLIN))
                    items_[i].revents |= ZMQ_POLLIN;
            }
            //  Else, the poll item is a raw file descriptor, simply convert
            //  the events to zmq_pollitem_t-style format.
            else {
                if (pollfds[i].revents & POLLIN)
                    items_[i].revents |= ZMQ_POLLIN;
                if (pollfds[i].revents & POLLOUT)
                    items_[i].revents |= ZMQ_POLLOUT;
                if (pollfds[i].revents & POLLPRI)
                    items_[i].revents |= ZMQ_POLLPRI;
                if (pollfds[i].revents & ~(POLLIN | POLLOUT | POLLPRI))
                    items_[i].revents |= ZMQ_POLLERR;
            }

            if (items_[i].revents)
                nevents++;
        }

#else

        //  Compute the timeout for the subsequent poll.
        struct timespec timeout;
        struct timespec *ptimeout;
        if (first_pass) {
            timeout.tv_sec = 0;
            timeout.tv_nsec = 0;
            ptimeout = &timeout;
        } else if (timeout_ < 0)
            ptimeout = NULL;
        else {
            timeout.tv_sec = static_cast<long> ((end - now) / 1000);
            timeout.tv_nsec = static_cast<long> ((end - now) % 1000 * 1000000);
            ptimeout = &timeout;
        }

        //  Wait for events. Ignore interrupts if there's infinite timeout.
        while (true) {
            memcpy (inset.get (), pollset_in.get (),
                    zmq::valid_pollset_bytes (*pollset_in.get ()));
            memcpy (outset.get (), pollset_out.get (),
                    zmq::valid_pollset_bytes (*pollset_out.get ()));
            memcpy (errset.get (), pollset_err.get (),
                    zmq::valid_pollset_bytes (*pollset_err.get ()));
#if defined ZMQ_HAVE_WINDOWS
            int rc =
              select (0, inset.get (), outset.get (), errset.get (), ptimeout);
            if (unlikely (rc == SOCKET_ERROR)) {
                errno = zmq::wsa_error_to_errno (WSAGetLastError ());
                wsa_assert (errno == ENOTSOCK);
                return -1;
            }
#else
            int rc = pselect (maxfd + 1, inset.get (), outset.get (),
                             errset.get (), ptimeout, sigmask_);
            if (unlikely (rc == -1)) {
                errno_assert (errno == EINTR || errno == EBADF);
                return -1;
            }
#endif
            break;
        }

        //  Check for the events.
        for (int i = 0; i != nitems_; i++) {
            items_[i].revents = 0;

            //  The poll item is a 0MQ socket. Retrieve pending events
            //  using the ZMQ_EVENTS socket option.
            if (items_[i].socket) {
                size_t zmq_events_size = sizeof (uint32_t);
                uint32_t zmq_events;
                if (zmq_getsockopt (items_[i].socket, ZMQ_EVENTS, &zmq_events,
                                    &zmq_events_size)
                    == -1)
                    return -1;
                if ((items_[i].events & ZMQ_POLLOUT)
                    && (zmq_events & ZMQ_POLLOUT))
                    items_[i].revents |= ZMQ_POLLOUT;
                if ((items_[i].events & ZMQ_POLLIN)
                    && (zmq_events & ZMQ_POLLIN))
                    items_[i].revents |= ZMQ_POLLIN;
            }
            //  Else, the poll item is a raw file descriptor, simply convert
            //  the events to zmq_pollitem_t-style format.
            else {
                if (FD_ISSET (items_[i].fd, inset.get ()))
                    items_[i].revents |= ZMQ_POLLIN;
                if (FD_ISSET (items_[i].fd, outset.get ()))
                    items_[i].revents |= ZMQ_POLLOUT;
                if (FD_ISSET (items_[i].fd, errset.get ()))
                    items_[i].revents |= ZMQ_POLLERR;
            }

            if (items_[i].revents)
                nevents++;
        }
#endif

        //  If timeout is zero, exit immediately whether there are events or not.
        if (timeout_ == 0)
            break;

        //  If there are events to return, we can exit immediately.
        if (nevents)
            break;

        //  At this point we are meant to wait for events but there are none.
        //  If timeout is infinite we can just loop until we get some events.
        if (timeout_ < 0) {
            if (first_pass)
                first_pass = false;
            continue;
        }

        //  The timeout is finite and there are no events. In the first pass
        //  we get a timestamp of when the polling have begun. (We assume that
        //  first pass have taken negligible time). We also compute the time
        //  when the polling should time out.
        if (first_pass) {
            now = clock.now_ms ();
            end = now + timeout_;
            if (now == end)
                break;
            first_pass = false;
            continue;
        }

        //  Find out whether timeout have expired.
        now = clock.now_ms ();
        if (now >= end)
            break;
    }

    return nevents;
#else
   //  Exotic platforms that support neither poll() nor select().
   errno = ENOTSUP;
   return -1;
#endif
}


// This function can throw, so wrap in try-catch!
int ppoll(zmq_pollitem_t *items_, size_t nitems_, long timeout_, const sigset_t * sigmask_)
{
   int rc = zmq_ppoll(items_, static_cast<int>(nitems_), timeout_, sigmask_);
   if (rc < 0)
      throw ppoll_error_t();
   return rc;
}


// This function can throw, so wrap in try-catch!
int ppoll(std::vector<zmq_pollitem_t> &items, long timeout_, const sigset_t * sigmask_)
{
   return ppoll(items.data(), items.size(), timeout_, sigmask_);
}

}

