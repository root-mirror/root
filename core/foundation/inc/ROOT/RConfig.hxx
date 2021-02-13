/* @(#)root/base:$Id$ */

/*************************************************************************
 * Copyright (C) 1995-2002, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_RConfig
#define ROOT_RConfig

/*************************************************************************
 *                                                                       *
 * RConfig                                                               *
 *                                                                       *
 * Defines used by ROOT.                                                 *
 *                                                                       *
 *************************************************************************/

#include "../RVersion.h"
#include "RConfigure.h"


/*---- new C++ features ------------------------------------------------------*/

#if defined __has_feature
# if __has_feature(modules)
#  define R__CXXMODULES
# endif
#endif

#define R__USE_SHADOW_CLASS

/* Now required, thus defined by default for backward compatibility */
#define R__ANSISTREAM      /* ANSI C++ Standard Library conformant */
#define R__SSTREAM         /* use sstream or strstream header */

#if defined(_MSC_VER)
# if (_MSC_VER < 1910)
#  error "ROOT requires Visual Studio 2017 or higher."
# else
#  define R__NULLPTR
# endif
#else
# if defined(__cplusplus) && (__cplusplus < 201103L)
#  error "ROOT requires support for C++11 or higher."
#  if defined(__GNUC__) || defined(__clang__)
#   error "Pass `-std=c++11` as compiler argument."
#  endif
# endif
#endif

/*---- machines --------------------------------------------------------------*/

#ifdef __hpux
    /* R__HPUX10 or R__HPUX11 is determined in the Makefile */
#   define R__HPUX
#   define R__UNIX
#   define ANSICPP
#   ifdef __LP64__
#      define R__B64
#   endif
#   ifdef R__HPUX10
#      define NEED_SNPRINTF
#   endif
#endif

#ifdef _AIX
#   define R__AIX
#   define R__UNIX
#   define ANSICPP
#   define R__SEEK64
#   define R__PLACEMENTINLINE /* placement new/delete is inline in <new> */
#   define NEED_STRCASECMP
#endif

#if defined(__linux) || defined(__linux__)
#   ifndef linux
#      define linux
#   endif
#endif

#if defined(__CYGWIN__) && defined(__GNUC__)
#   ifndef linux
#      define linux
#   endif
#   ifndef R__WINGCC
#      define R__WINGCC
#   endif
#endif

#if defined(__sun) && !(defined(linux) || defined(__FCC_VERSION))
#   ifdef __SVR4
#      define R__SOLARIS
#      define R__SEEK64
#      define ANSICPP
#      ifdef __i386
#         define R__BYTESWAP
#      endif
#      ifdef __x86_64
#         define R__B64
#         define R__BYTESWAP
#      endif
#   else
#      define R__SUN
#      include <stdlib.h>
#   endif
#   define R__UNIX
#   define NEED_STRING
#   define NEED_SIGJMP
#   if __SUNPRO_CC > 0x420
#      define R__SOLARIS_CC50
#      define R__PLACEMENTINLINE /* placement new/delete is inline in <new> */
#   endif
#   if __SUNPRO_CC >= 0x420
#      define R__SUNCCBUG        /* to work around a compiler bug */
#   endif
#   if __GNUC__ >= 3 || __GNUC_MINOR__ >= 90   /* modern egcs/gcc */
#      define R__SUNGCC3
#   endif
#endif

#if defined(__FCC_VERSION)    /* Solaris with Fujitsu compiler */
#   define R__SOLARIS
#   define R__SEEK64
#   define ANSICPP
#   define R__UNIX
#   define NEED_STRING
#   define NEED_SIGJMP
#endif

#if defined(linux)
#   ifndef _LARGEFILE64_SOURCE
#      define _LARGEFILE64_SOURCE
#   endif
#   include <features.h>
#   if __GNU_LIBRARY__ == 6
#      ifndef R__GLIBC
#         define R__GLIBC
#      endif
#   endif
#   if __GLIBC__ == 2 && __GLIBC_MINOR__ >= 2
#      define R__NONSCALARFPOS2
#      define R__USESTHROW
#      define R__SEEK64
#   endif
#endif

#if defined(linux) && defined(__i386__)
#   define R__LINUX
#   define R__UNIX
#   define R__BYTESWAP
#   ifndef __i486__
#      define __i486__       /* turn off if you really want to run on an i386 */
#   endif
#   define NEED_SIGJMP
#endif

#if defined(linux) && defined(__ia64__)
#   define R__LINUX
#   define R__UNIX
#   define R__BYTESWAP
#   define R__B64
#   define NEED_SIGJMP
#endif

#if defined(linux) && defined(__x86_64__)
#   define R__LINUX
#   define R__UNIX
#   define R__BYTESWAP
#   define R__B64
#   define NEED_SIGJMP
#endif

#if defined(linux) && defined(__arm__)
#   define R__LINUX
#   define R__UNIX
#   define R__BYTESWAP
#   define NEED_SIGJMP
#endif

#if defined(linux) && defined(__sparc__)
#   define R__LINUX
#   define R__UNIX
#   define NEED_SIGJMP
/*#   define R__B64 */     /* enable when 64 bit machine */
#endif

#if defined(linux) && defined(__hppa)
#   define R__LINUX
#   define R__UNIX
#   define NEED_SIGJMP
#endif

#if defined(linux) && defined(__powerpc__)
#   define R__LINUX
#   define R__UNIX
#   define NEED_SIGJMP
#   if defined(R__ppc64)
#      define R__B64
#   endif
#   if defined(_LITTLE_ENDIAN)
#      define R__BYTESWAP
#   endif
#endif

#if defined(linux) && defined(__aarch64__)
#   define R__LINUX
#   define R__UNIX
#   define R__BYTESWAP
#   define R__B64
#   define NEED_SIGJMP
#endif

#if defined(linux) && defined(__s390__)
#   define R__LINUX
#   define R__UNIX
#   define NEED_SIGJMP
#endif

#if defined(linux) && defined(__s390x__)
#   define R__LINUX
#   define R__UNIX
#   define R__B64
#   define NEED_SIGJMP
#endif

#if defined(__MACH__) && defined(__i386__) && !defined(__APPLE__)
#   define R__HURD
#   define f2cFortran   /* cfortran.h does not know HURD - sigh */
#   define R__UNIX
#   define R__BYTESWAP
#   define R__GLIBC     /* GNU/Hurd always use GLIBC 2.x :-) */
#   define NEED_SIGJMP
#endif

#if defined(__Lynx__) && defined(__powerpc__)
#   define R__LYNXOS
#   define R__UNIX
#   define ANSICPP
#   define NEED_SIGJMP
#   define NEED_STRCASECMP
#   define NEED_SNPRINTF
#endif

#if defined(__FreeBSD__)
#   define R__FBSD
#   define R__UNIX
#   define R__BYTESWAP
#   if defined(__i386__)
#      ifndef __i486__
#         define __i486__    /* turn off if you really want to run on an i386 */
#      endif
#   endif
#   if defined(__amd64__)
#      define R__B64
#   endif
#   define HAS_STRLCPY
#endif

#if defined(__OpenBSD__)
#   define R__OBSD
#   define R__UNIX
#   define R__BYTESWAP
#   if defined(__i386__)
#      ifndef __i486__
#         define __i486__    /* turn off if you really want to run on an i386 */
#      endif
#   endif
#   if defined(__amd64__)
#      define R__B64
#   endif
#   define HAS_STRLCPY
#endif

#if defined(__APPLE__)       /* MacOS X support, initially following FreeBSD */
#   include <AvailabilityMacros.h>
#   ifndef __CINT__
#   include <TargetConditionals.h>
#   endif
#   define R__MACOSX
#   define R__UNIX
#   if defined(__xlC__) || defined(__xlc__)
#      define ANSICPP
#      define R__PLACEMENTINLINE /* placement new/delete is inline in <new> */
#   endif
#   if defined(__ppc64__)
#      define R__B64      /* enable when 64 bit machine */
#   endif
#   if defined(__i386__)
#      define R__BYTESWAP
#   endif
#   if defined(__arm__) || defined (__arm64__)
#      define R__BYTESWAP
#   endif
#   if defined(__x86_64__)
#      define R__BYTESWAP
#      define R__B64      /* enable when 64 bit machine */
#   endif
#   define HAS_STRLCPY
#endif

#ifdef _HIUX_SOURCE
#   define R__HIUX
#   define R__UNIX
#   define NEED_SIGJMP
#   define NEED_SNPRINTF
#   define ANSICPP
#endif

#ifdef __GNUC__
#   define R__GNU
#   define ANSICPP
#   if __GNUC__ >= 3 || __GNUC_MINOR__ >= 90    /* egcs 1.0.3 */
#      define R__VECNEWDELETE    /* supports overloading of new[] and delete[] */
#      define R__PLACEMENTDELETE /* supports overloading placement delete */
#   endif
#   if __GNUC__ >= 3 || ( __GNUC__ == 2 && __GNUC_MINOR__ >= 95)
#         define R__PLACEMENTINLINE /* placement new/delete is inline in <new> */
#   endif
#   if defined(__ia64__) &&  __GNUC__ < 3       /* gcc 2.9x (MINOR is 9!) */
#      define R__VECNEWDELETE    /* supports overloading of new[] and delete[] */
#      define R__PLACEMENTDELETE /* supports overloading placement delete */
#   endif
#   if __GNUC__ > 4 || ( __GNUC__ == 4 && __GNUC_MINOR__ > 1)
#      define R__PRAGMA_DIAGNOSTIC
#   endif
#endif

#if __cplusplus >= 201402L
#   if defined(R__MACOSX) && !defined(MAC_OS_X_VERSION_10_12)
      // At least on 10.11, the compiler defines but the c++ library does not provide the size operator delete.
      // See for example https://llvm.org/bugs/show_bug.cgi?id=22951 or
      // https://github.com/gperftools/gperftools/issues/794.
#   elif !defined(__GNUC__)
#      define R__SIZEDDELETE
#   elif __GNUC__ > 4
#      define R__SIZEDDELETE
#   endif
#endif

/* allows symbols to be hidden from the shared library export symbol table */
/* use typically on file statics and private methods */
#if defined(__GNUC__) && (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 3))
#    define R__HIDDEN __attribute__((__visibility__("hidden")))
#else
#    define R__HIDDEN
#endif

#ifdef __INTEL_COMPILER
#   define R__INTEL_COMPILER
#   define R__VECNEWDELETE    /* supports overloading of new[] and delete[] */
#   define R__PLACEMENTDELETE /* supports overloading placement delete */
#   define R__PLACEMENTINLINE /* placement new/delete is inline in <new> */
#   define ANSICPP
#endif

#ifdef __HP_aCC
#   define R__ACC
#   define R__VECNEWDELETE    /* supports overloading of new[] and delete[] */
#   define R__PLACEMENTINLINE /* placement new/delete is inline in <new> */
#   if __HP_aCC <= 015000
#      define R__OLDHPACC
#      define R__TEMPLATE_OVERLOAD_BUG
#      define R__GLOBALSTL       /* STL in global name space */
#      error "ROOT requires proper support for C++11 or higher"
#   else
#      define R__PLACEMENTDELETE /* supports overloading placement delete */
#      define R__TMPLTSTREAM     /* std::iostream implemented with templates */
#   endif
#   ifndef _INCLUDE_LONGLONG
#      define _INCLUDE_LONGLONG
#   endif
#endif

#ifdef _WIN32
#   define R__WIN32
#   ifndef WIN32
#      define WIN32
#   endif
#   define R__BYTESWAP
#   define R__ACCESS_IN_SYMBOL
//#   define __attribute__(X)
//#   define thread_local static __declspec(thread)
#endif

#ifdef __SC__
#   define SC
#   define R__SC
#   if defined(WIN32)
#      define NEED_STRING
#      define NEED_STRCASECMP
#      define NEED_SNPRINTF
#      define ANSICPP
#   else
#      define MSDOS
#      define NEED_STRCASECMP
#      define R__BYTESWAP
#   endif
#endif

#ifdef _MSC_VER
#   define R__VISUAL_CPLUSPLUS
#   define NEED_STRING
#   define NEED_STRCASECMP
#   if _MSC_VER < 1900
#     define NEED_SNPRINTF
#   endif
#   define ANSICPP
#   define R__VECNEWDELETE    /* supports overloading of new[] and delete[] */
#   define R__PLACEMENTDELETE /* supports overloading placement delete */
#   define R__PLACEMENTINLINE /* placement new/delete is inline in <new> */
#   if _MSC_VER >= 1400
#     define DONTNEED_VSNPRINTF
#   endif
#   if _MSC_VER < 1310
#      define R__NO_CLASS_TEMPLATE_SPECIALIZATION
#   endif
#   if _MSC_VER <= 1800
#      define R__NO_ATOMIC_FUNCTION_POINTER
#   endif
#endif

/*--- memory and object statistics -------------------------------------------*/

/* #define R__NOSTATS */

/*--- cpp --------------------------------------------------------------------*/

#ifdef ANSICPP
    /* symbol concatenation operator */
#   define _NAME1_(name) name
#   define _NAME2_(name1,name2) name1##name2
#   define _NAME3_(name1,name2,name3) name1##name2##name3

    /* stringizing */
#   define _QUOTE_(name) #name

#else

#   define _NAME1_(name) name
#   define _NAME2_(name1,name2) _NAME1_(name1)name2
#   define _NAME3_(name1,name2,name3) _NAME2_(name1,name2)name3

#   define _QUOTE_(name) "name"

#endif

/* produce an identifier that is almost unique inside a file */
#ifndef __CINT__
#   define _R__JOIN_(X,Y) _NAME2_(X,Y)
#   define _R__JOIN3_(F,X,Y) _NAME3_(F,X,Y)
#   define _R__UNIQUE_DICT_(X) _R__JOIN3_(R__DICTIONARY_FILENAME,X,__LINE__)
#   define _R__UNIQUE_(X) _R__JOIN_(X,__LINE__)
#else
    /* Currently CINT does not really mind to have duplicates and     */
    /* does not work correctly as far as merging tokens is concerned. */
#   define _R__UNIQUE_(X) X
#endif

/*---- deprecation -----------------------------------------------------------*/
#if defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER)
# if (__GNUC__ == 5 && (__GNUC_MINOR__ == 1 || __GNUC_MINOR__ == 2)) || defined(R__NO_DEPRECATION)
/* GCC 5.1, 5.2: false positives due to https://gcc.gnu.org/bugzilla/show_bug.cgi?id=15269 
   or deprecation turned off */
#   define _R__DEPRECATED_LATER(REASON)
# else
#   define _R__DEPRECATED_LATER(REASON) __attribute__((deprecated(REASON)))
# endif
#elif defined(_MSC_VER) || !defined(R__NO_DEPRECATION)
#   define _R__DEPRECATED_LATER(REASON) __pragma(deprecated(REASON))
#else
/* Deprecation not supported for this compiler. */
#   define _R__DEPRECATED_LATER(REASON)
#endif

#ifdef R__WIN32
#define _R_DEPRECATED_REMOVE_NOW(REASON)
#else
#define _R_DEPRECATED_REMOVE_NOW(REASON) __attribute__((REMOVE_THIS_NOW))
#endif

/* To be removed by 6.20 */
#if ROOT_VERSION_CODE < ROOT_VERSION(6,19,0)
# define _R__DEPRECATED_620(REASON) _R__DEPRECATED_LATER(REASON)
#else
# define _R__DEPRECATED_620(REASON) _R_DEPRECATED_REMOVE_NOW(REASON)
#endif

/* To be removed by 6.22 */
#if ROOT_VERSION_CODE < ROOT_VERSION(6,21,0)
# define _R__DEPRECATED_622(REASON) _R__DEPRECATED_LATER(REASON)
#else
# define _R__DEPRECATED_622(REASON) _R_DEPRECATED_REMOVE_NOW(REASON)
#endif

/* To be removed by 6.24 */
#if ROOT_VERSION_CODE <= ROOT_VERSION(6,23,0)
# define _R__DEPRECATED_624(REASON) _R__DEPRECATED_LATER(REASON)
#else
# define _R__DEPRECATED_624(REASON) _R_DEPRECATED_REMOVE_NOW(REASON)
#endif

/* To be removed by 6.26 */
#if ROOT_VERSION_CODE <= ROOT_VERSION(6,25,0)
# define _R__DEPRECATED_626(REASON) _R__DEPRECATED_LATER(REASON)
#else
# define _R__DEPRECATED_626(REASON) _R_DEPRECATED_REMOVE_NOW(REASON)
#endif

/* To be removed by 7.00 */
#if ROOT_VERSION_CODE < ROOT_VERSION(6,99,0)
# define _R__DEPRECATED_700(REASON) _R__DEPRECATED_LATER(REASON)
#else
# define _R__DEPRECATED_700(REASON) _R_DEPRECATED_REMOVE_NOW(REASON)
#endif


/* Spell as R__DEPRECATED(6,04, "Not threadsafe; use TFoo::Bar().") */
#define R__DEPRECATED(MAJOR, MINOR, REASON) \
  _R__JOIN3_(_R__DEPRECATED_,MAJOR,MINOR)("will be removed in ROOT v" #MAJOR "." #MINOR ": " REASON)

/* Mechanisms to advise users to avoid legacy functions and classes that will not be removed */
#if defined R__SUGGEST_NEW_INTERFACE
#  define R__SUGGEST_ALTERNATIVE(ALTERNATIVE) \
      _R__DEPRECATED_LATER("There is a superior alternative: " ALTERNATIVE)
#else
#  define R__SUGGEST_ALTERNATIVE(ALTERNATIVE)
#endif

#define R__ALWAYS_SUGGEST_ALTERNATIVE(ALTERNATIVE) \
    _R__DEPRECATED_LATER("There is a superior alternative: " ALTERNATIVE)



/*---- misc ------------------------------------------------------------------*/

#ifdef R__GNU
#   define SafeDelete(p) { if (p) { delete p; p = 0; } }
#else
#   define SafeDelete(p) { delete p; p = 0; }
#endif

#ifdef __FAST_MATH__
#define R__FAST_MATH
#endif

#if (__GNUC__ >= 7)
#define R__DO_PRAGMA(x) _Pragma (#x)
# define R__INTENTIONALLY_UNINIT_BEGIN \
  R__DO_PRAGMA(GCC diagnostic push) \
  R__DO_PRAGMA(GCC diagnostic ignored "-Wmaybe-uninitialized") \
  R__DO_PRAGMA(GCC diagnostic ignored "-Wuninitialized")
# define R__INTENTIONALLY_UNINIT_END \
  R__DO_PRAGMA(GCC diagnostic pop)
#else
# define R__INTENTIONALLY_UNINIT_BEGIN
# define R__INTENTIONALLY_UNINIT_END

#endif

#ifdef R__HAS_ATTRIBUTE_ALWAYS_INLINE
#define R__ALWAYS_INLINE inline __attribute__((always_inline))
#else
#if defined(_MSC_VER)
#define R__ALWAYS_INLINE __forceinline
#else
#define R__ALWAYS_INLINE inline
#endif
#endif

// See also https://nemequ.github.io/hedley/api-reference.html#HEDLEY_NEVER_INLINE
// for other platforms.
#ifdef R__HAS_ATTRIBUTE_NOINLINE
#define R__NEVER_INLINE inline __attribute__((noinline))
#else
#if defined(_MSC_VER)
#define R__NEVER_INLINE inline  __declspec(noinline)
#else
#define R__NEVER_INLINE inline
#endif
#endif

/*---- unlikely / likely expressions -----------------------------------------*/
// These are meant to use in cases like:
//   if (R__unlikely(expression)) { ... }
// in performance-critical sections.  R__unlikely / R__likely provide hints to
// the compiler code generation to heavily optimize one side of a conditional,
// causing the other branch to have a heavy performance cost.
//
// It is best to use this for conditionals that test for rare error cases or
// backward compatibility code.

#if (__GNUC__ >= 3) || defined(__INTEL_COMPILER)
#if !defined(R__unlikely)
  #define R__unlikely(expr) __builtin_expect(!!(expr), 0)
#endif
#if !defined(R__likely)
  #define R__likely(expr) __builtin_expect(!!(expr), 1)
#endif
#else
  #define R__unlikely(expr) expr
  #define R__likely(expr) expr
#endif

// Setting this define causes ROOT to keep statistics about memory buffer allocation
// time within the TTree.  Given that this is a "hot-path", we provide a mechanism
// for enabling / disabling this at compile time by developers; default is disabled.
#ifndef R__TRACK_BASKET_ALLOC_TIME
//#define R__TRACK_BASKET_ALLOC_TIME 1
#endif

#endif
