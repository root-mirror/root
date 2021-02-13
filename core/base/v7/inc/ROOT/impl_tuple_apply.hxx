/// \file ROOT/impl_tuple_apply.hxx
/// \ingroup Base StdExt ROOT7
/// \author Axel Naumann <axel@cern.ch>
/// \date 2015-07-09
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!

/*************************************************************************
 * Copyright (C) 1995-2015, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_Impl_Tuple_Apply
#define ROOT7_Impl_Tuple_Apply

#include "RConfigure.h"

#include <functional>

// std::experimental::apply, invoke until it's there...
// from http://en.cppreference.com/w/cpp/utility/functional/invoke

#ifndef R__HAS_STD_INVOKE
namespace ROOT {
namespace Detail {
template <class F, class... Args>
inline auto INVOKE(F&& f, Args&&... args) ->
decltype(std::forward<F>(f)(std::forward<Args>(args)...)) {
  return std::forward<F>(f)(std::forward<Args>(args)...);
}

template <class Base, class T, class Derived>
inline auto INVOKE(T Base::*pmd, Derived&& ref) ->
decltype(std::forward<Derived>(ref).*pmd) {
  return std::forward<Derived>(ref).*pmd;
}

template <class PMD, class Pointer>
inline auto INVOKE(PMD pmd, Pointer&& ptr) ->
decltype((*std::forward<Pointer>(ptr)).*pmd) {
  return (*std::forward<Pointer>(ptr)).*pmd;
}

template <class Base, class T, class Derived, class... Args>
inline auto INVOKE(T Base::*pmf, Derived&& ref, Args&&... args) ->
decltype((std::forward<Derived>(ref).*pmf)(std::forward<Args>(args)...)) {
  return (std::forward<Derived>(ref).*pmf)(std::forward<Args>(args)...);
}

template <class PMF, class Pointer, class... Args>
inline auto INVOKE(PMF pmf, Pointer&& ptr, Args&&... args) ->
decltype(((*std::forward<Pointer>(ptr)).*pmf)(std::forward<Args>(args)...)) {
  return ((*std::forward<Pointer>(ptr)).*pmf)(std::forward<Args>(args)...);
}
} // namespace Detail
} // namespace ROOT

namespace std {
inline namespace __ROOT {

template< class F, class... ArgTypes>
decltype(auto) invoke(F&& f, ArgTypes&&... args) {
  return ROOT::Detail::INVOKE(std::forward<F>(f), std::forward<ArgTypes>(args)...);
}

} // inline namespace __ROOT {
} // namespace std
#endif // ndef R__HAS_STD_INVOKE

#ifndef R__HAS_STD_APPLY
// From http://en.cppreference.com/w/cpp/experimental/apply
namespace ROOT {
namespace Detail {
template<class F, class Tuple, std::size_t... I>
constexpr decltype(auto) apply_impl(F &&f, Tuple &&t,
                                    std::index_sequence<I...>) {
  return std::invoke(std::forward<F>(f),
                     std::get<I>(std::forward<Tuple>(t))...);
  // Note: std::invoke is a C++17 feature
}
} // namespace Detail
} // namespace ROOT

namespace std {
inline namespace __ROOT {
template<class F, class Tuple>
constexpr decltype(auto) apply(F &&f, Tuple &&t) {
  return ROOT::Detail::apply_impl(std::forward<F>(f), std::forward<Tuple>(t),
                            std::make_index_sequence < std::tuple_size <
                            std::decay_t < Tuple >> {} > {});
}
} // inline namespace __ROOT
} // namespace std
#endif // ndef R__HAS_STD_APPLY

#endif //ROOT7_TUPLE_APPLY_H
