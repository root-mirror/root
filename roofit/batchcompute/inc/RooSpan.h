// Author: Stephan Hageboeck, CERN  7 Feb 2019

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

#ifndef ROOFIT_BATCHCOMPUTE_ROOSPAN_H
#define ROOFIT_BATCHCOMPUTE_ROOSPAN_H

#include "ROOT/RSpan.hxx"
#include <vector>
#include <cassert>

////////////////////////////////////////////////////////////////////////////
/// A simple container to hold a batch of data values.
/// It can operate in two modes:
/// * Span: It holds only a pointer to the storage held by another object
/// like a std::span does.
/// * Temp data: It holds its own data, and exposes the span.
/// This mode is necessary to ship data that are not available in
/// a contiguous storage like e.g. data from a TTree. This means, however, that
/// data have to be copied, and only live as long as the span.
template<class T>
class RooSpan {
public:
  using iterator = typename std::span<T>::iterator;
  using value_type = typename std::remove_cv<T>::type;

  constexpr RooSpan() :
  _span{} { }

  constexpr RooSpan(RooSpan&& other) :
  _span{other._span.data(), other._span.size()}
  { }

  constexpr RooSpan(const RooSpan& other) :
  _span{other._span}
  { }


  /// Conversion constructor from <T> to <const T>
  /// If the input span owns some memory, the const-version of the
  /// span will copy the shared_ptr.
  template<typename NON_CONST_T,
      typename = typename std::enable_if<std::is_same<const NON_CONST_T, T>::value>::type >
  constexpr RooSpan(const RooSpan<NON_CONST_T>& other) :
  _span{other.data(), other.size()}
  { }


  /// Construct from a range. Data held by foreign object.
  template < class InputIterator>
  constexpr RooSpan(InputIterator beginIn, InputIterator endIn) :
  _span{beginIn, endIn}
  { }


  /// Construct from start and end pointers.
  constexpr RooSpan(typename std::span<T>::pointer beginIn,
      typename std::span<T>::pointer endIn) :
    _span{beginIn, endIn}
  { }


  /// Construct from start pointer and size.
  constexpr RooSpan(typename std::span<T>::pointer beginIn,
      typename std::span<T>::index_type sizeIn) :
  _span{beginIn, sizeIn}
  { }


  constexpr RooSpan(const std::vector<typename std::remove_cv<T>::type>& vec) noexcept :
  _span{vec}
  { }

  constexpr RooSpan(std::vector<typename std::remove_cv<T>::type>& vec) noexcept :
  _span{vec}
  { }


  /// We cannot point to temporary data.
  constexpr RooSpan(std::vector<value_type>&& payload) = delete;


  RooSpan<T>& operator=(const RooSpan<T>& other) = default;


  constexpr typename std::span<T>::iterator begin() const {
    return _span.begin();
  }

  constexpr typename std::span<T>::iterator end() const {
    return _span.end();
  }

  constexpr typename std::span<T>::pointer data() const {
    return _span.data();
  }

#ifdef NDEBUG
  constexpr typename std::span<T>::reference operator[](typename std::span<T>::index_type i) const noexcept {
    return _span[i];
  }
#else
  typename std::span<T>::reference operator[](typename std::span<T>::index_type i) const noexcept {
    assert(i < _span.size());
    return _span[i];
  }
#endif

  constexpr typename std::span<T>::index_type size() const noexcept {
    return _span.size();
  }

  constexpr bool empty() const noexcept {
    return _span.empty();
  }
  
  constexpr bool isBatch() const noexcept {
    return true;
  }


  ///Test if the span overlaps with `other`.
  template <class Span_t>
  constexpr bool overlaps(const Span_t& other) const {
    return insideSpan(other.begin()) || insideSpan(other.end()-1)
        || other.insideSpan(begin()) || other.insideSpan(end()-1);
  }

  ///Test if the given pointer/iterator is inside the span.
  template <typename ptr_t>
  constexpr bool insideSpan(ptr_t ptr) const {
    return begin() <= ptr && ptr < end();
  }

private:

  std::span<T> _span;
};


#endif /* ROOFIT_BATCHCOMPUTE_ROOSPAN_H */
