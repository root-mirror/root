// @(#)root/foundation:
// Author: Axel Naumann, Enrico Guiraud, June 2017

/*************************************************************************
 * Copyright (C) 1995-2017, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TypeTraits
#define ROOT_TypeTraits

#include <memory> // shared_ptr, unique_ptr for IsSmartOrDumbPtr
#include <type_traits>

namespace ROOT {

/// ROOT type_traits extensions
namespace TypeTraits {
/// Lightweight storage for a collection of types.
/// Differently from std::tuple, no instantiation of objects of stored types is performed
template <typename... Types>
struct TypeList {
   static constexpr std::size_t list_size = sizeof...(Types);
};
} // end ns TypeTraits

namespace Detail {
template <typename T> constexpr auto HasCallOp(int /*goodOverload*/) -> decltype(&T::operator(), true) { return true; }
template <typename T> constexpr bool HasCallOp(char /*badOverload*/) { return false; }

/// Extract types from the signature of a callable object. See CallableTraits.
template <typename T, bool HasCallOp = ROOT::Detail::HasCallOp<T>(0)>
struct CallableTraitsImpl {};

// Extract signature of operator() and delegate to the appropriate CallableTraitsImpl overloads
template <typename T>
struct CallableTraitsImpl<T, true> {
   using arg_types = typename CallableTraitsImpl<decltype(&T::operator())>::arg_types;
   using arg_types_nodecay = typename CallableTraitsImpl<decltype(&T::operator())>::arg_types_nodecay;
   using ret_type = typename CallableTraitsImpl<decltype(&T::operator())>::ret_type;
};

// lambdas, std::function, const member functions
template <typename R, typename T, typename... Args>
struct CallableTraitsImpl<R (T::*)(Args...) const, false> {
   using arg_types = ROOT::TypeTraits::TypeList<typename std::decay<Args>::type...>;
   using arg_types_nodecay = ROOT::TypeTraits::TypeList<Args...>;
   using ret_type = R;
};

// mutable lambdas and functor classes, non-const member functions
template <typename R, typename T, typename... Args>
struct CallableTraitsImpl<R (T::*)(Args...), false> {
   using arg_types = ROOT::TypeTraits::TypeList<typename std::decay<Args>::type...>;
   using arg_types_nodecay = ROOT::TypeTraits::TypeList<Args...>;
   using ret_type = R;
};

// function pointers
template <typename R, typename... Args>
struct CallableTraitsImpl<R (*)(Args...), false> {
   using arg_types = ROOT::TypeTraits::TypeList<typename std::decay<Args>::type...>;
   using arg_types_nodecay = ROOT::TypeTraits::TypeList<Args...>;
   using ret_type = R;
};

// free functions
template <typename R, typename... Args>
struct CallableTraitsImpl<R(Args...), false> {
   using arg_types = ROOT::TypeTraits::TypeList<typename std::decay<Args>::type...>;
   using arg_types_nodecay = ROOT::TypeTraits::TypeList<Args...>;
   using ret_type = R;
};
} // end ns Detail

namespace TypeTraits {

///\class ROOT::TypeTraits::
template <class T>
class IsSmartOrDumbPtr : public std::integral_constant<bool, std::is_pointer<T>::value> {
};

template <class P>
class IsSmartOrDumbPtr<std::shared_ptr<P>> : public std::true_type {
};

template <class P>
class IsSmartOrDumbPtr<std::unique_ptr<P>> : public std::true_type {
};

/// Checks for signed integers types that are not characters
template<class T>
struct IsSignedNumeral : std::integral_constant<bool,
   std::is_integral<T>::value &&
   std::is_signed<T>::value &&
   !std::is_same<T, char>::value
> {};

/// Checks for unsigned integer types that are not characters
template<class T>
struct IsUnsignedNumeral : std::integral_constant<bool,
   std::is_integral<T>::value &&
   !std::is_signed<T>::value &&
   !std::is_same<T, char>::value
> {};

/// Checks for floating point types (that are not characters)
template<class T>
using IsFloatNumeral = std::is_floating_point<T>;

/// Extract types from the signature of a callable object.
/// The `CallableTraits` struct contains three type aliases:
///   - arg_types: a `TypeList` of all types in the signature, decayed through std::decay
///   - arg_types_nodecay: a `TypeList` of all types in the signature, including cv-qualifiers
template<typename F>
using CallableTraits = ROOT::Detail::CallableTraitsImpl<F>;

// Return first of a variadic list of types.
template <typename T, typename... Rest>
struct TakeFirstType {
   using type = T;
};

template <typename... Types>
using TakeFirstType_t = typename TakeFirstType<Types...>::type;

// Remove first type from a variadic list of types, return a TypeList containing the rest.
// e.g. RemoveFirst_t<A,B,C> is TypeList<B,C>
template <typename T, typename... Rest>
struct RemoveFirst {
   using type = TypeList<Rest...>;
};

template <typename... Args>
using RemoveFirst_t = typename RemoveFirst<Args...>::type;

/// Return first of possibly many template parameters.
/// For non-template types, the result is void
/// e.g. TakeFirstParameter<U<A,B>> is A
///      TakeFirstParameter<T> is void
template <typename T>
struct TakeFirstParameter {
   using type = void;
};

template <template <typename...> class Template, typename T, typename... Rest>
struct TakeFirstParameter<Template<T, Rest...>> {
   using type = T;
};

template <typename T>
using TakeFirstParameter_t = typename TakeFirstParameter<T>::type;

/// Remove first of possibly many template parameters.
/// e.g. RemoveFirstParameter_t<U<A,B>> is U\<B\>
template <typename>
struct RemoveFirstParameter {
};

template <typename T, template <typename...> class U, typename... Rest>
struct RemoveFirstParameter<U<T, Rest...>> {
   using type = U<Rest...>;
};

template <typename T>
using RemoveFirstParameter_t = typename RemoveFirstParameter<T>::type;

template <typename T>
struct HasBeginAndEnd {

   template <typename V>
   using Begin_t = typename V::const_iterator (V::*)() const;

   template <typename V>
   using End_t = typename V::const_iterator (V::*)() const;

   template <typename V>
   static constexpr auto Check(int)
      -> decltype(static_cast<Begin_t<V>>(&V::begin), static_cast<End_t<V>>(&V::end), true)
   {
      return true;
   }

   template <typename V>
   static constexpr bool Check(...)
   {
      return false;
   }

   static constexpr bool const value = Check<T>(0);
};

} // ns TypeTraits
} // ns ROOT
#endif // ROOT_TTypeTraits
