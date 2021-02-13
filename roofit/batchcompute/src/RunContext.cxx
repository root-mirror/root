// Author: Stephan Hageboeck, CERN  Jul 2020

/*****************************************************************************
 * RooFit
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2020, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#include "RunContext.h"

#include <limits>

class RooAbsReal;

namespace RooBatchCompute {

/// Check if there is a span of data corresponding to the object passed as owner.
RooSpan<const double> RunContext::getBatch(const RooAbsReal* owner) const {
  const auto item = spans.find(owner);
  if (item != spans.end())
    return item->second;

  return {};
}


/// Check if there is a writable span of data corresponding to the object passed as owner.
/// The span can be used both for reading and writing.
RooSpan<double> RunContext::getWritableBatch(const RooAbsReal* owner) {
  auto item = ownedMemory.find(owner);
  if (item != ownedMemory.end()) {
    assert(spans.count(owner) > 0); // If we can write, the span must also be registered for reading
    return RooSpan<double>(item->second);
  }

  return {};
}


/// Create a writable batch. If the RunContext already owns memory for the object
/// `owner`, just resize the memory. If it doesn't exist yet, allocate it.
/// \warning The memory will be uninitialised, so every entry **must** be overwritten.
/// On first use, all values are initialised to `NaN` to help detect such errors.
///
/// A read-only reference to the memory will be stored in `spans`.
/// \param owner RooFit object whose value should be written into the memory.
/// \param size Requested size of the span.
/// \return A writeable RooSpan of the requested size, whose memory is owned by
/// the RunContext.
RooSpan<double> RunContext::makeBatch(const RooAbsReal* owner, std::size_t size) {
  auto item = ownedMemory.find(owner);
  if (item == ownedMemory.end() || item->second.size() != size) {
    std::vector<double>& data = ownedMemory[owner];
    data.resize(size, std::numeric_limits<double>::quiet_NaN());
#ifndef NDEBUG
    data.assign(size, std::numeric_limits<double>::quiet_NaN());
#endif
    spans[owner] = RooSpan<const double>(data);
    return {data};
  }

  spans[owner] = RooSpan<const double>(item->second);
  return {item->second};
}

}
