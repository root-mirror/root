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

#ifndef ROOFIT_MULTIPROCESS_UTIL_H
#define ROOFIT_MULTIPROCESS_UTIL_H

#include <unistd.h> // getpid, pid_t
namespace RooFit {
  namespace MultiProcessV1 {
    int wait_for_child(pid_t child_pid, bool may_throw, int retries_before_killing);
  }
}
#endif //ROOFIT_MULTIPROCESS_UTIL_H
