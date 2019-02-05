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

#ifndef ROOFIT_ROOFITCORE_INC_ROOSPAN_H_
#define ROOFIT_ROOFITCORE_INC_ROOSPAN_H_

#include "ROOT/RSpan.hxx"

////////////////////////////////////////////////////////////////////////////
/// A simple container to hold a batch of data values.
/// It can operate in two modes:
/// * Span: It holds only references to the storage held by another object
/// like a std::span does.
/// * Temp data: It holds its own data, and exposes the span.
/// This mode is necessary to ship data that are not available in
/// a contiguous storage like e.g. data from a TTree. This means, however, that
/// data have to be copied, and follow the DataBatch.
template<class T>
class RooSpan {
public:
  using iterator = typename std::span<T>::iterator;
  using value_type = typename std::remove_cv<T>::type;

  constexpr RooSpan(RooSpan&& other) :
  _auxStorage{std::move(other._auxStorage)},
  _span{other._span.data(), other._span.size()}
  { }

  constexpr RooSpan(const RooSpan& other) :
  _auxStorage{other._auxStorage},
  _span{other._span}
  { }

  /// Construct from a range. Data held by foreign object.
  template < class InputIterator>
  constexpr RooSpan(InputIterator begin, InputIterator end) :
  _auxStorage{},
  _span{begin, end}
  { }


  /// Construct from start and end pointers.
  constexpr RooSpan(typename std::span<T>::pointer begin,
      typename std::span<T>::pointer end) :
    _auxStorage{},
    _span{begin, end}
  { }


  /// Construct from start pointer and size.
  constexpr RooSpan(typename std::span<T>::pointer begin,
      typename std::span<T>::size_type size) :
  _auxStorage{},
  _span{begin, size}
  { }


  constexpr RooSpan(const std::vector<value_type>& vec) noexcept :
  _auxStorage{},
  _span{vec}
  { }


  /// Hand data over to this span. This will mean that the data will get
  /// deleted when it goes out of scope. Try to avoid this because
  /// unnecessary copies will be made.
  constexpr RooSpan(std::vector<value_type>&& payload) :
  _auxStorage{new std::vector<value_type>(std::forward<std::vector<value_type>>(payload))},
  _span{_auxStorage->begin(), _auxStorage->end()}
  { }


  constexpr typename std::span<T>::iterator begin() const {
    return _span.begin();
  }

  constexpr typename std::span<T>::iterator end() const {
    return _span.end();
  }

  constexpr typename std::span<T>::reference operator[](typename std::span<T>::index_type i) const noexcept {
    return *(static_cast<typename std::span<T>::pointer __restrict__>(_span.data())+i);
  }

  constexpr typename std::span<T>::index_type size() const noexcept {
    return _span.size();
  }

  ///Test if the span overlaps with `other`.
  template <class Span_t>
  bool overlaps(const Span_t& other) const {
    return insideSpan(other.begin()) || insideSpan(other.end()-1)
        || other.insideSpan(begin()) || other.insideSpan(end()-1);
  }

  ///Test if the given pointer/iterator is inside the span.
  template <typename ptr_t>
  bool insideSpan(ptr_t ptr) const {
    return begin() <= ptr && ptr < end();
  }

private:

  /// If a class has does not own a contiguous block of memory, which
  /// could be used to create a span, the memory has to be kept alive
  /// until all referring spans are destroyed.
  std::shared_ptr<std::vector<value_type>> _auxStorage;
  std::span<T> _span;
};


#endif /* ROOFIT_ROOFITCORE_INC_ROOSPAN_H_ */
