// Author: Stephan Hageboeck, CERN  25 Feb 2019

/*****************************************************************************
 * RooFit
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2019, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef ROOFIT_BATCHCOMPUTE_BRACKETADAPTER_H
#define ROOFIT_BATCHCOMPUTE_BRACKETADAPTER_H

#include "RooSpan.h"

#include <vector>

namespace RooBatchCompute {

///Little adapter that gives a bracket operator to types that don't
///have one. It completely ignores the index and returns a constant.
template <class T = double>
class BracketAdapter {
  public:

    constexpr BracketAdapter(T payload) noexcept :
    _payload{payload} { }

    constexpr BracketAdapter(RooSpan<const T> payload) noexcept :
    _payload{payload[0]} { }

    constexpr double operator[](std::size_t) const {
      return _payload;
    }

    constexpr operator double() const {
      return _payload;
    }

    constexpr bool isBatch() const noexcept {
      return false;
    }

  private:
    const T _payload;
};


class BracketAdapterWithMask {
  public:
    /// Construct adapter from a fallback value and a batch of values.
    /// - If `batch.size() == 0`, always return `payload`.
    /// - Else, return `batch[i]`.
    BracketAdapterWithMask(double payload, const RooSpan<const double>& batch) noexcept :
    _isBatch(!batch.empty()),
    _payload(payload),
    _pointer(batch.empty() ? &_payload : batch.data()),
    _mask(batch.size() > 1 ? ~static_cast<size_t>(0): 0)
    {
    }
    
    /// Construct adapter from a batch of values.
    /// - If `batch.size() == 1`, always return the value at `batch[0]`.
    /// - Else, return `batch[i]`.
    BracketAdapterWithMask(RooSpan<const double> batch) :
    _isBatch(batch.size() > 1),
    _payload(batch[0]),
    _pointer(batch.data()),
    _mask(batch.size() > 1 ? ~static_cast<size_t>(0): 0)
    {
      assert(batch.size() > 0);
    }

    BracketAdapterWithMask(const BracketAdapterWithMask& other) noexcept:
    _isBatch(other._isBatch),
    _payload(other._payload),
    _pointer(other._isBatch ? other._pointer : &_payload),
    _mask(other._mask)
    {
    }
    
    BracketAdapterWithMask& operator= (const BracketAdapterWithMask& other) = delete;

    inline double operator[](std::size_t i) const noexcept {
      return _pointer[ i & _mask];
    }
    
    inline bool isBatch() const noexcept {
      return _isBatch;
    }

  private:
    const bool _isBatch;
    const double _payload;
    const double* __restrict const _pointer;
    const size_t _mask;
};

}

#endif /* ROOFIT_BATCHCOMPUTE_BRACKETADAPTER_H */
