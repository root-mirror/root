#include "ROOT/TVec.hxx"

#define _VECOPS_USE_EXTERN_TEMPLATES true

// We do not support extern templates on Win
#ifdef _WIN32
#undef _VECOPS_USE_EXTERN_TEMPLATES
#define _VECOPS_USE_EXTERN_TEMPLATES false
#endif // _WIN32

// We do not support extern templates on Linux if the compiler is old
#ifdef R_LINUX
#if (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__) <= 40800
#undef _VECOPS_USE_EXTERN_TEMPLATES
#define _VECOPS_USE_EXTERN_TEMPLATES false
#endif // GCC version
#ifdef __clang__
#undef _VECOPS_USE_EXTERN_TEMPLATES
#define _VECOPS_USE_EXTERN_TEMPLATES false
#endif
#endif // R__LINUX

#if(_VECOPS_USE_EXTERN_TEMPLATES)

namespace ROOT {
namespace Experimental {
namespace VecOps {

#define TVEC_DECLARE_UNARY_OPERATOR(T, OP) \
   template TVec<T> operator OP(const TVec<T> &);

#define TVEC_DECLARE_BINARY_OPERATOR(T, OP)                                              \
   template auto operator OP(const TVec<T> &v, const T &y) -> TVec<decltype(v[0] OP y)>; \
   template auto operator OP(const T &x, const TVec<T> &v) -> TVec<decltype(x OP v[0])>; \
   template auto operator OP(const TVec<T> &v0, const TVec<T> &v1) -> TVec<decltype(v0[0] OP v1[0])>;

#define TVEC_DECLARE_LOGICAL_OPERATOR(T, OP)                   \
   template TVec<int> operator OP(const TVec<T> &, const T &); \
   template TVec<int> operator OP(const T &, const TVec<T> &); \
   template TVec<int> operator OP(const TVec<T> &, const TVec<T> &);

#define TVEC_DECLARE_ASSIGN_OPERATOR(T, OP)             \
   template TVec<T> &operator OP(TVec<T> &, const T &); \
   template TVec<T> &operator OP(TVec<T> &, const TVec<T> &);

#define TVEC_DECLARE_FLOAT_TEMPLATE(T)  \
   template class TVec<T>;              \
   TVEC_DECLARE_UNARY_OPERATOR(T, +)    \
   TVEC_DECLARE_UNARY_OPERATOR(T, -)    \
   TVEC_DECLARE_UNARY_OPERATOR(T, !)    \
   TVEC_DECLARE_BINARY_OPERATOR(T, +)   \
   TVEC_DECLARE_BINARY_OPERATOR(T, -)   \
   TVEC_DECLARE_BINARY_OPERATOR(T, *)   \
   TVEC_DECLARE_BINARY_OPERATOR(T, /)   \
   TVEC_DECLARE_ASSIGN_OPERATOR(T, +=)  \
   TVEC_DECLARE_ASSIGN_OPERATOR(T, -=)  \
   TVEC_DECLARE_ASSIGN_OPERATOR(T, *=)  \
   TVEC_DECLARE_ASSIGN_OPERATOR(T, /=)  \
   TVEC_DECLARE_LOGICAL_OPERATOR(T, <)  \
   TVEC_DECLARE_LOGICAL_OPERATOR(T, >)  \
   TVEC_DECLARE_LOGICAL_OPERATOR(T, ==) \
   TVEC_DECLARE_LOGICAL_OPERATOR(T, !=) \
   TVEC_DECLARE_LOGICAL_OPERATOR(T, <=) \
   TVEC_DECLARE_LOGICAL_OPERATOR(T, >=) \
   TVEC_DECLARE_LOGICAL_OPERATOR(T, &&) \
   TVEC_DECLARE_LOGICAL_OPERATOR(T, ||)

#define TVEC_DECLARE_INTEGER_TEMPLATE(T) \
   template class TVec<T>;               \
   TVEC_DECLARE_UNARY_OPERATOR(T, +)     \
   TVEC_DECLARE_UNARY_OPERATOR(T, -)     \
   TVEC_DECLARE_UNARY_OPERATOR(T, ~)     \
   TVEC_DECLARE_UNARY_OPERATOR(T, !)     \
   TVEC_DECLARE_BINARY_OPERATOR(T, +)    \
   TVEC_DECLARE_BINARY_OPERATOR(T, -)    \
   TVEC_DECLARE_BINARY_OPERATOR(T, *)    \
   TVEC_DECLARE_BINARY_OPERATOR(T, /)    \
   TVEC_DECLARE_BINARY_OPERATOR(T, %)    \
   TVEC_DECLARE_BINARY_OPERATOR(T, &)    \
   TVEC_DECLARE_BINARY_OPERATOR(T, |)    \
   TVEC_DECLARE_BINARY_OPERATOR(T, ^)    \
   TVEC_DECLARE_ASSIGN_OPERATOR(T, +=)   \
   TVEC_DECLARE_ASSIGN_OPERATOR(T, -=)   \
   TVEC_DECLARE_ASSIGN_OPERATOR(T, *=)   \
   TVEC_DECLARE_ASSIGN_OPERATOR(T, /=)   \
   TVEC_DECLARE_ASSIGN_OPERATOR(T, %=)   \
   TVEC_DECLARE_ASSIGN_OPERATOR(T, &=)   \
   TVEC_DECLARE_ASSIGN_OPERATOR(T, |=)   \
   TVEC_DECLARE_ASSIGN_OPERATOR(T, ^=)   \
   TVEC_DECLARE_ASSIGN_OPERATOR(T, >>=)  \
   TVEC_DECLARE_ASSIGN_OPERATOR(T, <<=)  \
   TVEC_DECLARE_LOGICAL_OPERATOR(T, <)   \
   TVEC_DECLARE_LOGICAL_OPERATOR(T, >)   \
   TVEC_DECLARE_LOGICAL_OPERATOR(T, ==)  \
   TVEC_DECLARE_LOGICAL_OPERATOR(T, !=)  \
   TVEC_DECLARE_LOGICAL_OPERATOR(T, <=)  \
   TVEC_DECLARE_LOGICAL_OPERATOR(T, >=)  \
   TVEC_DECLARE_LOGICAL_OPERATOR(T, &&)  \
   TVEC_DECLARE_LOGICAL_OPERATOR(T, ||)

TVEC_DECLARE_INTEGER_TEMPLATE(char)
TVEC_DECLARE_INTEGER_TEMPLATE(short)
TVEC_DECLARE_INTEGER_TEMPLATE(int)
TVEC_DECLARE_INTEGER_TEMPLATE(long)
TVEC_DECLARE_INTEGER_TEMPLATE(long long)

TVEC_DECLARE_INTEGER_TEMPLATE(unsigned char)
TVEC_DECLARE_INTEGER_TEMPLATE(unsigned short)
TVEC_DECLARE_INTEGER_TEMPLATE(unsigned int)
TVEC_DECLARE_INTEGER_TEMPLATE(unsigned long)
TVEC_DECLARE_INTEGER_TEMPLATE(unsigned long long)

TVEC_DECLARE_FLOAT_TEMPLATE(float)
TVEC_DECLARE_FLOAT_TEMPLATE(double)

#define TVEC_DECLARE_UNARY_FUNCTION(T, NAME, FUNC) \
   template TVec<T> NAME(const TVec<T> &);

#define TVEC_DECLARE_STD_UNARY_FUNCTION(T, F) TVEC_DECLARE_UNARY_FUNCTION(T, F, ::std::F)

#define TVEC_DECLARE_BINARY_FUNCTION(T0, T1, NAME, FUNC) \
   template TVec<decltype((T0() + T1()))> NAME(const TVec<T0> &, const T1 &); \
   template TVec<decltype((T0() + T1()))> NAME(const T0 &, const TVec<T1> &); \
   template TVec<decltype((T0() + T1()))> NAME(const TVec<T0> &, const TVec<T1> &);

#define TVEC_DECLARE_STD_BINARY_FUNCTION(T, F) TVEC_DECLARE_BINARY_FUNCTION(T, T, F, ::std::F)

#define TVEC_DECLARE_STD_FUNCTIONS(T)             \
   TVEC_DECLARE_STD_UNARY_FUNCTION(T, abs)        \
   TVEC_DECLARE_STD_BINARY_FUNCTION(T, fdim)      \
   TVEC_DECLARE_STD_BINARY_FUNCTION(T, fmod)      \
   TVEC_DECLARE_STD_BINARY_FUNCTION(T, remainder) \
   TVEC_DECLARE_STD_UNARY_FUNCTION(T, exp)        \
   TVEC_DECLARE_STD_UNARY_FUNCTION(T, exp2)       \
   TVEC_DECLARE_STD_UNARY_FUNCTION(T, expm1)      \
   TVEC_DECLARE_STD_UNARY_FUNCTION(T, log)        \
   TVEC_DECLARE_STD_UNARY_FUNCTION(T, log10)      \
   TVEC_DECLARE_STD_UNARY_FUNCTION(T, log2)       \
   TVEC_DECLARE_STD_UNARY_FUNCTION(T, log1p)      \
   TVEC_DECLARE_STD_BINARY_FUNCTION(T, pow)       \
   TVEC_DECLARE_STD_UNARY_FUNCTION(T, sqrt)       \
   TVEC_DECLARE_STD_UNARY_FUNCTION(T, cbrt)       \
   TVEC_DECLARE_STD_BINARY_FUNCTION(T, hypot)     \
   TVEC_DECLARE_STD_UNARY_FUNCTION(T, sin)        \
   TVEC_DECLARE_STD_UNARY_FUNCTION(T, cos)        \
   TVEC_DECLARE_STD_UNARY_FUNCTION(T, tan)        \
   TVEC_DECLARE_STD_UNARY_FUNCTION(T, asin)       \
   TVEC_DECLARE_STD_UNARY_FUNCTION(T, acos)       \
   TVEC_DECLARE_STD_UNARY_FUNCTION(T, atan)       \
   TVEC_DECLARE_STD_BINARY_FUNCTION(T, atan2)     \
   TVEC_DECLARE_STD_UNARY_FUNCTION(T, sinh)       \
   TVEC_DECLARE_STD_UNARY_FUNCTION(T, cosh)       \
   TVEC_DECLARE_STD_UNARY_FUNCTION(T, tanh)       \
   TVEC_DECLARE_STD_UNARY_FUNCTION(T, asinh)      \
   TVEC_DECLARE_STD_UNARY_FUNCTION(T, acosh)      \
   TVEC_DECLARE_STD_UNARY_FUNCTION(T, atanh)      \
   TVEC_DECLARE_STD_UNARY_FUNCTION(T, floor)      \
   TVEC_DECLARE_STD_UNARY_FUNCTION(T, ceil)       \
   TVEC_DECLARE_STD_UNARY_FUNCTION(T, trunc)      \
   TVEC_DECLARE_STD_UNARY_FUNCTION(T, round)      \
   TVEC_DECLARE_STD_UNARY_FUNCTION(T, lround)     \
   TVEC_DECLARE_STD_UNARY_FUNCTION(T, llround)    \
   TVEC_DECLARE_STD_UNARY_FUNCTION(T, erf)        \
   TVEC_DECLARE_STD_UNARY_FUNCTION(T, erfc)       \
   TVEC_DECLARE_STD_UNARY_FUNCTION(T, lgamma)     \
   TVEC_DECLARE_STD_UNARY_FUNCTION(T, tgamma)     \

TVEC_DECLARE_STD_FUNCTIONS(float)
TVEC_DECLARE_STD_FUNCTIONS(double)
#undef TVEC_DECLARE_STD_UNARY_FUNCTION
#undef TVEC_DECLARE_STD_BINARY_FUNCTION
#undef TVEC_DECLARE_STD_UNARY_FUNCTIONS

#ifdef R__HAS_VDT

#define TVEC_DECLARE_VDT_UNARY_FUNCTION(T, F)    \
   TVEC_DECLARE_UNARY_FUNCTION(T, F, vdt::F)

TVEC_DECLARE_VDT_UNARY_FUNCTION(float, fast_expf)
TVEC_DECLARE_VDT_UNARY_FUNCTION(float, fast_logf)
TVEC_DECLARE_VDT_UNARY_FUNCTION(float, fast_sinf)
TVEC_DECLARE_VDT_UNARY_FUNCTION(float, fast_cosf)
TVEC_DECLARE_VDT_UNARY_FUNCTION(float, fast_tanf)
TVEC_DECLARE_VDT_UNARY_FUNCTION(float, fast_asinf)
TVEC_DECLARE_VDT_UNARY_FUNCTION(float, fast_acosf)
TVEC_DECLARE_VDT_UNARY_FUNCTION(float, fast_atanf)

TVEC_DECLARE_VDT_UNARY_FUNCTION(double, fast_exp)
TVEC_DECLARE_VDT_UNARY_FUNCTION(double, fast_log)
TVEC_DECLARE_VDT_UNARY_FUNCTION(double, fast_sin)
TVEC_DECLARE_VDT_UNARY_FUNCTION(double, fast_cos)
TVEC_DECLARE_VDT_UNARY_FUNCTION(double, fast_tan)
TVEC_DECLARE_VDT_UNARY_FUNCTION(double, fast_asin)
TVEC_DECLARE_VDT_UNARY_FUNCTION(double, fast_acos)
TVEC_DECLARE_VDT_UNARY_FUNCTION(double, fast_atan)

#endif // R__HAS_VDT

} // namespace VecOps
} // namespace Experimental
} // namespace ROOT

#endif // _VECOPS_USE_EXTERN_TEMPLATES
