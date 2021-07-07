#ifndef SERIALIZE_UTILITY_H
#define SERIALIZE_UTILITY_H 1

#if defined(__clang__) && (__clang_major__ > 3 || (__clang_major__ == 3 && __clang_minor__ >= 10))
#define HAVE_TRIVIALLY_COPYABLE 1
#elif defined(__GNUC__) && __GNUC__ >= 5
#define HAVE_TRIVIALLY_COPYABLE 1
#else
#undef HAVE_TRIVIALLY_COPYABLE
#endif

#include <type_traits>

namespace ZMQ {
namespace Detail {

#if defined(HAVE_TRIVIALLY_COPYABLE)
template <class T>
using simple_object = std::is_trivially_copyable<T>;
#else
template <class T>
using simple_object = std::is_pod<T>;
#endif

// is trivial
template <class T>
struct is_trivial
   : std::conditional<simple_object<typename std::decay<T>::type>::value, std::true_type, std::false_type>::type {
};

} // namespace Detail
} // namespace ZMQ

#endif // SERIALIZE_UTILITY_H
