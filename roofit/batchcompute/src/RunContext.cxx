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

/**
 * \struct RooBatchCompute::RunContext
 *
 * This struct enables passing computation data around between elements of a computation graph.
 *
 * ### Separating data and computation graph
 * The RunContext stores read-only spans to data that has already been computed.
 * This can be data of the observables (which is constant during a fit)
 * or intermediate computation results from evaluating PDFs or formulae for every point in a dataset.
 * The latter may change as fit parameters change.
 *
 * Instead of storing information about these data *inside* nodes of the computation graph (requiring a change
 * of their state, possibly violating const-correctness), this information is stored in RunContext::spans using
 * the pointer of the element that produced those results as a key. In this way, one or multiple RunContext
 * instances can be passed around when computations are running, leaving the objects of the computation graph
 * invariant.
 *
 * ### Memory ownership model
 * The RunContext can provide memory for temporary data, that is, data that can vanish after a fit converges. Using
 * RunContext::makeBatch(), a suitable amount of memory is allocated to store computation results.
 * When intermediate data are cleared, this memory is *not freed*. In this way, temporary data can be invalidated
 * when fit parameters change, but the memory is only allocated once per fit.
 *
 * When a RunContext goes out of scope, the memory is freed. That means that in between fit cycles, a RunContext should
 * be cleared using clear(), or single results should be invalidated by removing these from RunContext::spans.
 * The RunContext object should be destroyed only *after* a fit completes.
 */


#include "RunContext.h"
#include "RooBatchCompute.h"

#include <limits>
#include <iostream>

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

void RunContext::clear() {
  spans.clear();
  //~  dispatch->clearCuda(*this);
  rangeName = nullptr;
}

}
