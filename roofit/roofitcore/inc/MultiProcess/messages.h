/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id$
 * Authors:                                                                  *
 *   PB, Patrick Bos,     NL eScience Center, p.bos@esciencecenter.nl        *
 *   IP, Inti Pelupessy,  NL eScience Center, i.pelupessy@esciencecenter.nl  *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef ROOFIT_MULTIPROCESS_MESSAGES_H
#define ROOFIT_MULTIPROCESS_MESSAGES_H

#include <iostream>

namespace RooFit {
  namespace MultiProcessV1 {

    // Messages from master to queue
    enum class M2Q : int {
      terminate = 100,
      enqueue = 10,
      retrieve = 11,
      update_real = 12,
//      update_cat = 13,
      call_double_const_method = 15,
      flush_ostreams = 16
    };

    // Messages from queue to master
    enum class Q2M : int {
      retrieve_rejected = 20,
      retrieve_accepted = 21
    };

    // Messages from worker to queue
    enum class W2Q : int {
      dequeue = 30,
      send_result = 31
    };

    // Messages from queue to worker
    enum class Q2W : int {
      terminate = 400,
      dequeue_rejected = 40,
      dequeue_accepted = 41,
      result_received = 43,
      update_real = 44,
//      update_cat = 45
      call_double_const_method = 46,
      flush_ostreams = 47
    };

    // stream output operators for debugging
    std::ostream& operator<<(std::ostream& out, const M2Q value);
    std::ostream& operator<<(std::ostream& out, const Q2M value);
    std::ostream& operator<<(std::ostream& out, const Q2W value);
    std::ostream& operator<<(std::ostream& out, const W2Q value);

  } // namespace MultiProcessV1

//  // forward declaration:
//  class BidirMMapPipe;
//
//  // bipe stream operators for message enum classes
//  BidirMMapPipe& operator<<(BidirMMapPipe& bipe, const MultiProcess::M2Q& sent);
//  BidirMMapPipe& operator>>(BidirMMapPipe& bipe, MultiProcess::M2Q& received);
//  BidirMMapPipe& operator<<(BidirMMapPipe& bipe, const MultiProcess::Q2M& sent);
//  BidirMMapPipe& operator>>(BidirMMapPipe& bipe, MultiProcess::Q2M& received);
//  BidirMMapPipe& operator<<(BidirMMapPipe& bipe, const MultiProcess::Q2W& sent);
//  BidirMMapPipe& operator>>(BidirMMapPipe& bipe, MultiProcess::Q2W& received);
//  BidirMMapPipe& operator<<(BidirMMapPipe& bipe, const MultiProcess::W2Q& sent);
//  BidirMMapPipe& operator>>(BidirMMapPipe& bipe, MultiProcess::W2Q& received);

} // namespace RooFit

#endif //ROOFIT_MULTIPROCESS_MESSAGES_H
