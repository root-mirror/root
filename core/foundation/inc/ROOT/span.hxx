/// \file ROOT/span.hxx
/// \ingroup Base StdExt
/// \author Axel Naumann <axel@cern.ch>
/// \date 2015-09-06

/*************************************************************************
 * Copyright (C) 1995-2015, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_RHYSD_SPAN_H
#define ROOT_RHYSD_SPAN_H

// Necessary to compile in c++11 mode
#if __cplusplus >= 201402L
#define R__CONSTEXPR_IF_CXX14 constexpr
#else
#define R__CONSTEXPR_IF_CXX14
#endif

// From https://github.com/rhysd/array_view/blob/master/include/array_view.hpp

#include <cstddef>
#include <iterator>
#include <array>
#include <vector>
#include <stdexcept>
#include <memory>
#include <type_traits>
#include <initializer_list>

namespace ROOT {
namespace Detail {
using std::size_t;

// detail meta functions {{{
template<class Array>
struct is_array_class {
  static bool const value = false;
};
template<class T, size_t N>
struct is_array_class<std::array<T, N>> {
  static bool const value = true;
};
template<class T>
struct is_array_class<std::vector<T>> {
  static bool const value = true;
};
template<class T>
struct is_array_class<std::initializer_list<T>> {
  static bool const value = true;
};
// }}}

// index sequences {{{
template< size_t... Indices >
struct indices{
  static constexpr size_t value[sizeof...(Indices)] = {Indices...};
};

template<class IndicesType, size_t Next>
struct make_indices_next;

template<size_t... Indices, size_t Next>
struct make_indices_next<indices<Indices...>, Next> {
  typedef indices<Indices..., (Indices + Next)...> type;
};

template<class IndicesType, size_t Next, size_t Tail>
struct make_indices_next2;

template<size_t... Indices, size_t Next, size_t Tail>
struct make_indices_next2<indices<Indices...>, Next, Tail> {
  typedef indices<Indices..., (Indices + Next)..., Tail> type;
};

template<size_t First, size_t Step, size_t N, class = void>
struct make_indices_impl;

template<size_t First, size_t Step, size_t N>
struct make_indices_impl<
   First,
   Step,
   N,
   typename std::enable_if<(N == 0)>::type
> {
  typedef indices<> type;
};

template<size_t First, size_t Step, size_t N>
struct make_indices_impl<
   First,
   Step,
   N,
   typename std::enable_if<(N == 1)>::type
> {
  typedef indices<First> type;
};

template<size_t First, size_t Step, size_t N>
struct make_indices_impl<
   First,
   Step,
   N,
   typename std::enable_if<(N > 1 && N % 2 == 0)>::type
>
   : ROOT::Detail::make_indices_next<
      typename ROOT::Detail::make_indices_impl<First, Step, N / 2>::type,
      First + N / 2 * Step
   >
{};

template<size_t First, size_t Step, size_t N>
struct make_indices_impl<
   First,
   Step,
   N,
   typename std::enable_if<(N > 1 && N % 2 == 1)>::type
>
   : ROOT::Detail::make_indices_next2<
      typename ROOT::Detail::make_indices_impl<First, Step, N / 2>::type,
      First + N / 2 * Step,
      First + (N - 1) * Step
   >
{};

template<size_t First, size_t Last, size_t Step = 1>
struct make_indices_
   : ROOT::Detail::make_indices_impl<
      First,
      Step,
      ((Last - First) + (Step - 1)) / Step
   >
{};

template < size_t Start, size_t Last, size_t Step = 1 >
using make_indices = typename make_indices_< Start, Last, Step >::type;
// }}}
} // namespace Detail
}

namespace std {

inline namespace __ROOT {

// span {{{

struct check_bound_t {};
static constexpr check_bound_t check_bound{};

template<class T>
class span {
public:
  /*
   * types
   */
  typedef T element_type;
  typedef std::remove_cv<T>     value_type;
  typedef element_type *        pointer;
  typedef element_type const*   const_pointer;
  typedef element_type &        reference;
  typedef element_type const&   const_reference;
  typedef element_type *        iterator;
  typedef element_type const*   const_iterator;
  typedef ptrdiff_t             difference_type;
  typedef std::size_t           index_type;
  typedef std::reverse_iterator<iterator> reverse_iterator;
  typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

  /*
   * ctors and assign operators
   */
  constexpr span() noexcept
     : length_(0), data_(nullptr)
  {}

  constexpr span(span const&) noexcept = default;
  constexpr span(span &&) noexcept = default;

  // Note:
  // This constructor can't be constexpr because & operator can't be constexpr.
  template<index_type N>
  /*implicit*/ span(std::array<T, N> & a) noexcept
     : length_(N), data_(N > 0 ? a.data() : nullptr)
  {}

  // Note:
  // This constructor can't be constexpr because & operator can't be constexpr.
  template<index_type N>
  /*implicit*/ span(T(& a)[N]) noexcept
     : length_(N), data_(N > 0 ? std::addressof(a[0]) : nullptr)
  {
    static_assert(N > 0, "Zero-length array is not permitted in ISO C++.");
  }

  /*implicit*/ span(std::vector<typename std::remove_cv<T>::type> const& v) noexcept
     : length_(v.size()), data_(v.empty() ? nullptr : v.data())
  {}

  /*implicit*/ span(std::vector<typename std::remove_cv<T>::type> & v) noexcept
     : length_(v.size()), data_(v.empty() ? nullptr : v.data())
  {}

  /*implicit*/ constexpr span(pointer a, index_type const n) noexcept
     : length_(n), data_(a)
  {}

  template<
     class InputIterator,
     class = typename std::enable_if<
        std::is_same<
           typename std::remove_cv<T>::type,
           typename std::iterator_traits<InputIterator>::value_type
        >::value
     >::type
  >
  span(InputIterator start, InputIterator last)
     : length_(std::distance(start, last)), data_(&*start)
  {}

  span(std::initializer_list<T> const& l)
     : length_(l.size()), data_(std::begin(l))
  {}

  span& operator=(span const&) noexcept = default;
  span& operator=(span &&) noexcept = delete;

  /*
   * iterator interfaces
   */
  constexpr iterator begin() const noexcept
  {
    return data_;
  }
  constexpr iterator end() const noexcept
  {
    return data_ + length_;
  }
  constexpr const_iterator cbegin() const noexcept
  {
    return begin();
  }
  constexpr const_iterator cend() const noexcept
  {
    return end();
  }
  reverse_iterator rbegin() const
  {
    return {end()};
  }
  reverse_iterator rend() const
  {
    return {begin()};
  }
  const_reverse_iterator crbegin() const
  {
    return rbegin();
  }
  const_reverse_iterator crend() const
  {
    return rend();
  }

  /*
   * access
   */
  constexpr index_type size() const noexcept
  {
    return length_;
  }
  constexpr index_type length() const noexcept
  {
    return size();
  }
  constexpr index_type max_size() const noexcept
  {
    return size();
  }
  constexpr bool empty() const noexcept
  {
    return length_ == 0;
  }
  constexpr reference operator[](index_type const n) const noexcept
  {
    return *(data_ + n);
  }
  constexpr reference at(index_type const n) const
  {
    //Works only in C++14
    //if (n >= length_) throw std::out_of_range("span::at()");
    //return *(data_ + n);
    return n >= length_ ? throw std::out_of_range("span::at()") : *(data_ + n);
  }
  constexpr pointer data() const noexcept
  {
    return data_;
  }
  constexpr const_reference front() const noexcept
  {
    return *data_;
  }
  constexpr const_reference back() const noexcept
  {
    return *(data_ + length_ - 1);
  }

  /*
   * slices
   */
  // slice with indices {{{
  // check bound {{{
  constexpr span<T> slice(check_bound_t, index_type const pos, index_type const slicelen) const
  {
    //Works only in C++14
    //if (pos >= length_ || pos + slicelen >= length_) {
    //  throw std::out_of_range("span::slice()");
    //}
    //return span<T>{begin() + pos, begin() + pos + slicelen};
    return pos >= length_ || pos + slicelen >= length_ ? throw std::out_of_range("span::slice()") : span<T>{begin() + pos, begin() + pos + slicelen};
  }
  constexpr span<T> slice_before(check_bound_t, index_type const pos) const
  {
    //Works only in C++14
    //if (pos >= length_) {
    //  throw std::out_of_range("span::slice()");
    //}
    //return span<T>{begin(), begin() + pos};
    return pos >= length_ ? std::out_of_range("span::slice()") : span<T>{begin(), begin() + pos};
  }
  constexpr span<T> slice_after(check_bound_t, index_type const pos) const
  {
    //Works only in C++14
    //if (pos >= length_) {
    //  throw std::out_of_range("span::slice()");
    //}
    //return span<T>{begin() + pos, end()};
    return pos >= length_ ? std::out_of_range("span::slice()") : span<T>{begin() + pos, end()};
  }
  // }}}
  // not check bound {{{
  constexpr span<T> slice(index_type const pos, index_type const slicelen) const
  {
    return span<T>{begin() + pos, begin() + pos + slicelen};
  }
  constexpr span<T> slice_before(index_type const pos) const
  {
    return span<T>{begin(), begin() + pos};
  }
  constexpr span<T> slice_after(index_type const pos) const
  {
    return span<T>{begin() + pos, end()};
  }
  // }}}
  // }}}
  // slice with iterators {{{
  // check bound {{{
  constexpr span<T> slice(check_bound_t, iterator start, iterator last) const
  {
    //Works only in C++14
    //if ( start >= end() ||
    //     last > end() ||
    //     start > last ||
    //     static_cast<size_t>(std::distance(start, last > end() ? end() : last)) > length_ - std::distance(begin(), start) ) {
    //  throw std::out_of_range("span::slice()");
    //}
    //return span<T>{start, last > end() ? end() : last};
    return ( start >= end() ||
             last > end() ||
             start > last ||
             static_cast<size_t>(std::distance(start, last > end() ? end() : last)) > length_ - std::distance(begin(), start) ) ? throw std::out_of_range("span::slice()") : span<T>{start, last > end() ? end() : last};
  }
  constexpr span<T> slice_before(check_bound_t, iterator const pos) const
  {
    //Works only in C++14
    //if (pos < begin() || pos > end()) {
    //  throw std::out_of_range("span::slice()");
    //}
    //return span<T>{begin(), pos > end() ? end() : pos};
    return pos < begin() || pos > end() ? throw std::out_of_range("span::slice()") : span<T>{begin(), pos > end() ? end() : pos};
  }
  constexpr span<T> slice_after(check_bound_t, iterator const pos) const
  {
    //Works only in C++14
    // if (pos < begin() || pos > end()) {
    //  throw std::out_of_range("span::slice()");
    //}
    //return span<T>{pos < begin() ? begin() : pos, end()};
    return pos < begin() || pos > end() ? throw std::out_of_range("span::slice()") : span<T>{pos < begin() ? begin() : pos, end()};
  }
  // }}}
  // not check bound {{{
  constexpr span<T> slice(iterator start, iterator last) const
  {
    return span<T>{start, last};
  }
  constexpr span<T> slice_before(iterator const pos) const
  {
    return span<T>{begin(), pos};
  }
  constexpr span<T> slice_after(iterator const pos) const
  {
    return span<T>{pos, end()};
  }
  // }}}
  // }}}

  /*
   * others
   */
  template<class Allocator = std::allocator<typename std::remove_cv<T>::type>>
  auto to_vector(Allocator const& alloc = Allocator{}) const
  -> std::vector<typename std::remove_cv<T>::type, Allocator>
  {
    return {begin(), end(), alloc};
  }

  template<size_t N>
  auto to_array() const
  -> std::array<T, N>
  {
    return to_array_impl(ROOT::Detail::make_indices<0, N>{});
  }
private:
  template<size_t... I>
  auto to_array_impl(ROOT::Detail::indices<I...>) const
  -> std::array<T, sizeof...(I)>
  {
    return {{(I < length_ ? *(data_ + I) : T{} )...}};
  }

private:
  index_type length_;
  pointer data_;
};
// }}}
} // inline namespace __ROOT
} // namespace std

namespace ROOT {
// compare operators {{{
namespace Detail {

template< class ArrayL, class ArrayR >
inline R__CONSTEXPR_IF_CXX14
bool operator_equal_impl(ArrayL const& lhs, size_t const lhs_size, ArrayR const& rhs, size_t const rhs_size)
{
  if (lhs_size != rhs_size) {
    return false;
  }

  auto litr = std::begin(lhs);
  auto ritr = std::begin(rhs);
  for (; litr != std::end(lhs); ++litr, ++ritr) {
    if (!(*litr == *ritr)) {
      return false;
    }
  }

  return true;
}
} // namespace Detail
} // namespace ROOT

namespace std {
inline namespace __ROOT {

template<class T1, class T2>
inline constexpr
bool operator==(span<T1> const& lhs, span<T2> const& rhs)
{
  return ROOT::Detail::operator_equal_impl(lhs, lhs.length(), rhs, rhs.length());
}

template<
   class T,
   class Array,
   class = typename std::enable_if<
      ROOT::Detail::is_array_class<Array>::value
   >::type
>
inline constexpr
bool operator==(span<T> const& lhs, Array const& rhs)
{
  return ROOT::Detail::operator_equal_impl(lhs, lhs.length(), rhs, rhs.size());
}

template<class T1, class T2, size_t N>
inline constexpr
bool operator==(span<T1> const& lhs, T2 const (& rhs)[N])
{
  return ROOT::Detail::operator_equal_impl(lhs, lhs.length(), rhs, N);
}

template<
   class T,
   class Array,
   class = typename std::enable_if<
      is_array<Array>::value
   >::type
>
inline constexpr
bool operator!=(span<T> const& lhs, Array const& rhs)
{
  return !(lhs == rhs);
}

template<
   class Array,
   class T,
   class = typename std::enable_if<
      is_array<Array>::value
   >::type
>
inline constexpr
bool operator==(Array const& lhs, span<T> const& rhs)
{
  return rhs == lhs;
}

template<
   class Array,
   class T,
   class = typename std::enable_if<
      is_array<Array>::value,
      Array
   >::type
>
inline constexpr
bool operator!=(Array const& lhs, span<T> const& rhs)
{
  return !(rhs == lhs);
}
// }}}

// helpers to construct view {{{
template<
   class Array,
   class = typename std::enable_if<
      ROOT::Detail::is_array_class<Array>::value
   >::type
>
inline constexpr
auto make_view(Array const& a)
-> span<typename Array::value_type>
{
  return {a};
}

template< class T, size_t N>
inline constexpr
span<T> make_view(T const (&a)[N])
{
  return {a};
}

template<class T>
inline constexpr
span<T> make_view(T const* p, typename span<T>::index_type const n)
{
  return span<T>{p, n};
}

template<class InputIterator, class Result = span<typename std::iterator_traits<InputIterator>::value_type>>
inline constexpr
Result make_view(InputIterator begin, InputIterator end)
{
  return Result{begin, end};
}

template<class T>
inline constexpr
span<T> make_view(std::initializer_list<T> const& l)
{
  return {l};
}
// }}}

} // inline namespace __ROOT
} // namespace std





#if 0
// This stuff is too complex for our simple use case!

#include <cstddef>
#include <array>
#include <type_traits>

// See N3851

namespace std {

template<int Rank>
class index;

template<int Rank>
class bounds {
public:
  static constexpr int rank = Rank;
  using reference = ptrdiff_t &;
  using const_reference = const ptrdiff_t &;
  using size_type = size_t;
  using value_type = ptrdiff_t;

private:
  std::array<value_type, Rank> m_B;

public:
  constexpr bounds() noexcept;

  constexpr bounds(value_type b) noexcept: m_B{{b}} { };
  //constexpr bounds(const initializer_list<value_type>&) noexcept;
  //constexpr bounds(const bounds&) noexcept;
  //bounds& operator=(const bounds&) noexcept;

  reference operator[](size_type i) noexcept { return m_B[i]; }

  constexpr const_reference operator[](
     size_type i) const noexcept { return m_B[i]; };


  bool operator==(const bounds &rhs) const noexcept;

  bool operator!=(const bounds &rhs) const noexcept;

  bounds operator+(const index<rank> &rhs) const noexcept;

  bounds operator-(const index<rank> &rhs) const noexcept;

  bounds &operator+=(const index<rank> &rhs) noexcept;

  bounds &operator-=(const index<rank> &rhs) noexcept;

  constexpr size_type size() const noexcept;

  bool contains(const index<rank> &idx) const noexcept;
  //bounds_iterator<rank> begin() const noexcept;
  //bounds_iterator<rank> end() const noexcept;

};

//bounds operator+(const index<rank>& lhs, const bounds& rhs) noexcept;

template<int Rank>
class index {
public:
  static constexpr int rank = Rank;
  using reference = ptrdiff_t &;
  using const_reference = const ptrdiff_t &;
  using size_type = size_t;
  using value_type = ptrdiff_t;

// For index<rank>:
  constexpr index() noexcept;

  constexpr index(value_type) noexcept;

  constexpr index(const initializer_list<value_type> &) noexcept;

  constexpr index(const index &) noexcept;

  index &operator=(const index &) noexcept;

  reference operator[](size_type component_idx) noexcept;

  constexpr const_reference operator[](size_type component_idx) const noexcept;

  bool operator==(const index &rhs) const noexcept;

  bool operator!=(const index &rhs) const noexcept;

  index operator+(const index &rhs) const noexcept;

  index operator-(const index &rhs) const noexcept;

  index &operator+=(const index &rhs) noexcept;

  index &operator-=(const index &rhs) noexcept;

  index &operator++() noexcept;

  index operator++(int) noexcept;

  index &operator--() noexcept;

  index operator--(int) noexcept;

  index operator+() const noexcept;

  index operator-() const noexcept;
};

/// Mock-up of future atd::(experimental::)span.
/// Supports only what we need for THist, e.g. Rank := 1.
template<typename ValueType, int Rank = 1>
class span {
public:
  static constexpr int rank = Rank;
  using index_type = index<rank>;
  using bounds_type = bounds<rank>;
  using size_type = typename bounds_type::size_type;
  using value_type = ValueType;
  using pointer = typename std::add_pointer_t<value_type>;
  using reference = typename std::add_lvalue_reference_t<value_type>;

  constexpr span() noexcept;

  constexpr explicit span(std::vector<ValueType> &cont) noexcept;

  template<typename ArrayType>
  constexpr explicit span(ArrayType &data) noexcept;

  template<typename ViewValueType>
  constexpr span(const span<ViewValueType, rank> &rhs) noexcept;

  template<typename Container>
  constexpr span(bounds_type bounds, Container &cont) noexcept;

  constexpr span(bounds_type bounds, pointer data) noexcept;

  template<typename ViewValueType>
  span &operator=(const span<ViewValueType, rank> &rhs) noexcept;

  constexpr bounds_type bounds() const noexcept;
  constexpr size_type size() const noexcept;
  constexpr index_type stride() const noexcept;

  constexpr pointer data() const noexcept;
  constexpr reference operator[](const index_type& idx) const noexcept;
};

}
#endif // too complex!
#endif
