// -*- C++ -*-
// Author: Philippe Canal, March 2015

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//===---------------------wrap string_view ----------------------------===//
// Wrapper header adapting the snapshot of sring_view.h to build without
// the libc++ infrastructure header files.

#ifndef RWrap_libcpp_string_view_h
#define RWrap_libcpp_string_view_h

#ifndef RStringView_H
#error "Do not use RWrap_libcpp_string_view.h directly. #include \"ROOT/RStringView.hxx\" instead."
#endif // RStringView_H

// In case we are connected with a libc++ which defines those, we need
// to include it first, so we avoid being silently over-ridden.

// To import a new version of the original file do:
//
/*
 svn co http://llvm.org/svn/llvm-project/libcxx/trunk libcxx

 cat original/string_view | \
 sed -e 's:_LIBCPP_BEGIN_NAMESPACE_LFTS:_ROOT_LIBCPP_BEGIN_NAMESPACE_LFTS:' \
     -e 's:_LIBCPP_END_NAMESPACE_LFTS:_ROOT_LIBCPP_END_NAMESPACE_LFTS:' \
     -e 's:#include <__debug>://#include <__debug>:' \
     -e 's:#include <experimental/__config>://#include <experimental/__config>:' \
     -e 's:__put_character_sequence:R__put_character_sequence:g' \
 > core/base/inc/libcpp_string_view.h
*/

#include <string>
#include <algorithm>
#include <iterator>
#include <ostream>
#include <iomanip>
#include <stdexcept>

#ifndef R__WIN32

#define _ROOT_LIBCPP_BEGIN_NAMESPACE_LFTS \
namespace std { \
namespace experimental { inline namespace __ROOT {
#define _ROOT_LIBCPP_END_NAMESPACE_LFTS } } }
#else

// Microsoft compiler does not support inline namespace yet.
#define _ROOT_LIBCPP_BEGIN_NAMESPACE_LFTS \
namespace std { \
namespace experimental { namespace __ROOT {
#define _ROOT_LIBCPP_END_NAMESPACE_LFTS } using namespace __ROOT; } }

#endif


#ifndef _LIBCPP_BEGIN_NAMESPACE_STD
#define _LOCAL_LIBCPP_BEGIN_NAMESPACE_STD
#define _LIBCPP_BEGIN_NAMESPACE_STD namespace std {
#define _LIBCPP_END_NAMESPACE_STD }
#endif

#ifndef _LIBCPP_CONSTEXPR
#define _LOCAL_LIBCPP_CONSTEXPR
#define _LIBCPP_CONSTEXPR constexpr
#endif

#ifndef _VSTD
#define _LOCAL_VSTD
#define _VSTD ::std
#endif

#ifndef _LIBCPP_INLINE_VISIBILITY
#define _LOCAL_LIBCPP_INLINE_VISIBILITY
#define _LIBCPP_INLINE_VISIBILITY inline
#endif

#ifndef _LIBCPP_EXPLICIT
#define _LOCAL_LIBCPP_EXPLICIT
#define _LIBCPP_EXPLICIT explicit
#endif

//#ifndef _LIBCPP_CONSTEXPR_AFTER_CXX11
//#define _LOCAL_LIBCPP_CONSTEXPR_AFTER_CXX11
//#define _LIBCPP_CONSTEXPR_AFTER_CXX11 constexpr
//#endif

#ifdef _LIBCPP_STD_VER
#define _LOCAL_LIBCPP_STD_VER
#define _LIBCPP_STD_VER 11
#endif

#ifndef _LIBCPP_TYPE_VIS_ONLY
#define _LOCAL_LIBCPP_TYPE_VIS_ONLY
#define _LIBCPP_TYPE_VIS_ONLY
#endif

#ifndef _LIBCPP_CONSTEXPR_AFTER_CXX11
#define _LOCAL_LIBCPP_CONSTEXPR_AFTER_CXX11
#define _LIBCPP_CONSTEXPR_AFTER_CXX11
#endif

#ifndef _NOEXCEPT
#define _LOCAL_NOEXCEPT
#define _NOEXCEPT
#endif

#ifndef _LIBCPP_ASSERT
#define _LOCAL_LIBCPP_ASSERT
#define _LIBCPP_ASSERT(X,Y) ((void)0)
#endif

/* Also used:
 _LIBCPP_TYPE_VIS_ONLY
 _LIBCPP_CONSTEXPR_AFTER_CXX11
 */

namespace std {
#ifdef _LOCAL_VSTD
inline namespace __ROOT {
#else
// libC++ wins.
inline namespace __1 {
#endif

   //   template <class _Traits>
   //   struct _LIBCPP_HIDDEN __traits_eq
   //   {
   //      typedef typename _Traits::char_type char_type;
   //      _LIBCPP_INLINE_VISIBILITY
   //      bool operator()(const char_type& __x, const char_type& __y) _NOEXCEPT
   //      {return _Traits::eq(__x, __y);}
   //   };
   template <class _CharT, class _OutputIterator>
   _OutputIterator
   __pad_and_output(_OutputIterator __s,
                    const _CharT* __ob, const _CharT* __op, const _CharT* __oe,
                    ios_base& __iob, _CharT __fl)
   {
       streamsize __sz = __oe - __ob;
       streamsize __ns = __iob.width();
       if (__ns > __sz)
           __ns -= __sz;
       else
           __ns = 0;
       for (;__ob < __op; ++__ob, ++__s)
           *__s = *__ob;
       for (; __ns; --__ns, ++__s)
           *__s = __fl;
       for (; __ob < __oe; ++__ob, ++__s)
           *__s = *__ob;
       __iob.width(0);
       return __s;
   }

   template<class _CharT, class _Traits>
   basic_ostream<_CharT, _Traits>&
   R__put_character_sequence(basic_ostream<_CharT, _Traits>& __os,
                            const _CharT* __str, size_t __len)
   {
#if 0
//#ifndef _LIBCPP_NO_EXCEPTIONS
      try
      {
#endif  // _LIBCPP_NO_EXCEPTIONS
         typename basic_ostream<_CharT, _Traits>::sentry __s(__os);
         if (__s)
         {
            typedef ostreambuf_iterator<_CharT, _Traits> _Ip;
            if (__pad_and_output(_Ip(__os),
                                 __str,
                                 (__os.flags() & ios_base::adjustfield) == ios_base::left ?
                                 __str + __len :
                                 __str,
                                 __str + __len,
                                 __os,
                                 __os.fill()).failed())
               __os.setstate(ios_base::badbit | ios_base::failbit);
         }
#if 0
//#ifndef _LIBCPP_NO_EXCEPTIONS
      }
      catch (...)
      {
         __os.__set_badbit_and_consider_rethrow();
      }
#endif  // _LIBCPP_NO_EXCEPTIONS
      return __os;
   }

#ifdef _LOCAL_VSTD
   // search

   template <class _BinaryPredicate, class _ForwardIterator1, class _ForwardIterator2>
   _ForwardIterator1
   __search(_ForwardIterator1 __first1, _ForwardIterator1 __last1,
            _ForwardIterator2 __first2, _ForwardIterator2 __last2, _BinaryPredicate __pred,
            forward_iterator_tag, forward_iterator_tag)
   {
      if (__first2 == __last2)
         return __first1;  // Everything matches an empty sequence
      while (true)
      {
         // Find first element in sequence 1 that matchs *__first2, with a mininum of loop checks
         while (true)
         {
            if (__first1 == __last1)  // return __last1 if no element matches *__first2
               return __last1;
            if (__pred(*__first1, *__first2))
               break;
            ++__first1;
         }
         // *__first1 matches *__first2, now match elements after here
         _ForwardIterator1 __m1 = __first1;
         _ForwardIterator2 __m2 = __first2;
         while (true)
         {
            if (++__m2 == __last2)  // If pattern exhausted, __first1 is the answer (works for 1 element pattern)
               return __first1;
            if (++__m1 == __last1)  // Otherwise if source exhaused, pattern not found
               return __last1;
            if (!__pred(*__m1, *__m2))  // if there is a mismatch, restart with a new __first1
            {
               ++__first1;
               break;
            }  // else there is a match, check next elements
         }
      }
   }

   template <class _BinaryPredicate, class _RandomAccessIterator1, class _RandomAccessIterator2>
   _LIBCPP_CONSTEXPR_AFTER_CXX11 _RandomAccessIterator1
   __search(_RandomAccessIterator1 __first1, _RandomAccessIterator1 __last1,
            _RandomAccessIterator2 __first2, _RandomAccessIterator2 __last2, _BinaryPredicate __pred,
            random_access_iterator_tag, random_access_iterator_tag)
   {
      typedef typename std::iterator_traits<_RandomAccessIterator1>::difference_type _D1;
      typedef typename std::iterator_traits<_RandomAccessIterator2>::difference_type _D2;
      // Take advantage of knowing source and pattern lengths.  Stop short when source is smaller than pattern
      _D2 __len2 = __last2 - __first2;
      if (__len2 == 0)
         return __first1;
      _D1 __len1 = __last1 - __first1;
      if (__len1 < __len2)
         return __last1;
      const _RandomAccessIterator1 __s = __last1 - (__len2 - 1);  // Start of pattern match can't go beyond here
      while (true)
      {
#if !_LIBCPP_UNROLL_LOOPS
         while (true)
         {
            if (__first1 == __s)
               return __last1;
            if (__pred(*__first1, *__first2))
               break;
            ++__first1;
         }
#else  // !_LIBCPP_UNROLL_LOOPS
         for (_D1 __loop_unroll = (__s - __first1) / 4; __loop_unroll > 0; --__loop_unroll)
         {
            if (__pred(*__first1, *__first2))
               goto __phase2;
            if (__pred(*++__first1, *__first2))
               goto __phase2;
            if (__pred(*++__first1, *__first2))
               goto __phase2;
            if (__pred(*++__first1, *__first2))
               goto __phase2;
            ++__first1;
         }
         switch (__s - __first1)
         {
            case 3:
               if (__pred(*__first1, *__first2))
                  break;
               ++__first1;
            case 2:
               if (__pred(*__first1, *__first2))
                  break;
               ++__first1;
            case 1:
               if (__pred(*__first1, *__first2))
                  break;
            case 0:
               return __last1;
         }
      __phase2:
#endif  // !_LIBCPP_UNROLL_LOOPS
         _RandomAccessIterator1 __m1 = __first1;
         _RandomAccessIterator2 __m2 = __first2;
#if !_LIBCPP_UNROLL_LOOPS
         while (true)
         {
            if (++__m2 == __last2)
               return __first1;
            ++__m1;          // no need to check range on __m1 because __s guarantees we have enough source
            if (!__pred(*__m1, *__m2))
            {
               ++__first1;
               break;
            }
         }
#else  // !_LIBCPP_UNROLL_LOOPS
         ++__m2;
         ++__m1;
         for (_D2 __loop_unroll = (__last2 - __m2) / 4; __loop_unroll > 0; --__loop_unroll)
         {
            if (!__pred(*__m1, *__m2))
               goto __continue;
            if (!__pred(*++__m1, *++__m2))
               goto __continue;
            if (!__pred(*++__m1, *++__m2))
               goto __continue;
            if (!__pred(*++__m1, *++__m2))
               goto __continue;
            ++__m1;
            ++__m2;
         }
         switch (__last2 - __m2)
         {
            case 3:
               if (!__pred(*__m1, *__m2))
                  break;
               ++__m1;
               ++__m2;
            case 2:
               if (!__pred(*__m1, *__m2))
                  break;
               ++__m1;
               ++__m2;
            case 1:
               if (!__pred(*__m1, *__m2))
                  break;
            case 0:
               return __first1;
         }
      __continue:
         ++__first1;
#endif  // !_LIBCPP_UNROLL_LOOPS
      }
   }

#endif // _LOCAL_VSTD for __search

   template <class _ForwardIterator1, class _ForwardIterator2, class _BinaryPredicate>
   _LIBCPP_CONSTEXPR_AFTER_CXX11 _ForwardIterator1
   __find_first_of_ce(_ForwardIterator1 __first1, _ForwardIterator1 __last1,
                      _ForwardIterator2 __first2, _ForwardIterator2 __last2, _BinaryPredicate __pred)
   {
      for (; __first1 != __last1; ++__first1)
         for (_ForwardIterator2 __j = __first2; __j != __last2; ++__j)
            if (__pred(*__first1, *__j))
               return __first1;
      return __last1;
   }

   // __str_find
   template<class _CharT, class _SizeT, class _Traits, _SizeT __npos>
   _SizeT _LIBCPP_CONSTEXPR_AFTER_CXX11 _LIBCPP_INLINE_VISIBILITY
   __str_find(const _CharT *__p, _SizeT __sz,
              _CharT __c, _SizeT __pos) _NOEXCEPT
   {
   if (__pos >= __sz)
      return __npos;
   const _CharT* __r = _Traits::find(__p + __pos, __sz - __pos, __c);
   if (__r == 0)
      return __npos;
   return static_cast<_SizeT>(__r - __p);
}

   template<class _CharT, class _SizeT, class _Traits, _SizeT __npos>
   _SizeT _LIBCPP_CONSTEXPR_AFTER_CXX11 _LIBCPP_INLINE_VISIBILITY
   __str_find(const _CharT *__p, _SizeT __sz,
              const _CharT* __s, _SizeT __pos, _SizeT __n)
   {
     if (__pos > __sz || __sz - __pos < __n)
        return __npos;
      if (__n == 0)
         return __pos;
      const _CharT* __r =
      _VSTD::__search(__p + __pos, __p + __sz,
                      __s, __s + __n, _Traits::eq,
                      random_access_iterator_tag(), random_access_iterator_tag());
      if (__r == __p + __sz)
         return __npos;
      return static_cast<_SizeT>(__r - __p);
   }


   // __str_rfind

   template<class _CharT, class _SizeT, class _Traits, _SizeT __npos>
   _SizeT _LIBCPP_CONSTEXPR_AFTER_CXX11 _LIBCPP_INLINE_VISIBILITY
   __str_rfind(const _CharT *__p, _SizeT __sz,
               _CharT __c, _SizeT __pos)
   {
      if (__sz < 1)
         return __npos;
      if (__pos < __sz)
         ++__pos;
      else
         __pos = __sz;
      for (const _CharT* __ps = __p + __pos; __ps != __p;)
      {
         if (_Traits::eq(*--__ps, __c))
            return static_cast<_SizeT>(__ps - __p);
      }
      return __npos;
   }

   template<class _CharT, class _SizeT, class _Traits, _SizeT __npos>
   _SizeT _LIBCPP_CONSTEXPR_AFTER_CXX11 _LIBCPP_INLINE_VISIBILITY
   __str_rfind(const _CharT *__p, _SizeT __sz,
               const _CharT* __s, _SizeT __pos, _SizeT __n)
   {
      __pos = _VSTD::min(__pos, __sz);
      if (__n < __sz - __pos)
         __pos += __n;
      else
         __pos = __sz;
      const _CharT* __r = _VSTD::find_end(
                                            __p, __p + __pos, __s, __s + __n, _Traits::eq,
                                            random_access_iterator_tag(), random_access_iterator_tag());
      if (__n > 0 && __r == __p + __pos)
         return __npos;
      return static_cast<_SizeT>(__r - __p);
   }

   // __str_find_first_of
   template<class _CharT, class _SizeT, class _Traits, _SizeT __npos>
   _SizeT _LIBCPP_CONSTEXPR_AFTER_CXX11 _LIBCPP_INLINE_VISIBILITY
   __str_find_first_of(const _CharT *__p, _SizeT __sz,
                       const _CharT* __s, _SizeT __pos, _SizeT __n)
   {
      if (__pos >= __sz || __n == 0)
         return __npos;
      const _CharT* __r = _VSTD::__find_first_of_ce
      (__p + __pos, __p + __sz, __s, __s + __n, _Traits::eq );
      if (__r == __p + __sz)
         return __npos;
      return static_cast<_SizeT>(__r - __p);
   }


   // __str_find_last_of
   template<class _CharT, class _SizeT, class _Traits, _SizeT __npos>
   _SizeT _LIBCPP_CONSTEXPR_AFTER_CXX11 _LIBCPP_INLINE_VISIBILITY
   __str_find_last_of(const _CharT *__p, _SizeT __sz,
                      const _CharT* __s, _SizeT __pos, _SizeT __n)
   {
      if (__n != 0)
      {
         if (__pos < __sz)
            ++__pos;
         else
            __pos = __sz;
         for (const _CharT* __ps = __p + __pos; __ps != __p;)
         {
            const _CharT* __r = _Traits::find(__s, __n, *--__ps);
            if (__r)
               return static_cast<_SizeT>(__ps - __p);
         }
      }
      return __npos;
   }


   // __str_find_first_not_of
   template<class _CharT, class _SizeT, class _Traits, _SizeT __npos>
   _SizeT _LIBCPP_CONSTEXPR_AFTER_CXX11 _LIBCPP_INLINE_VISIBILITY
   __str_find_first_not_of(const _CharT *__p, _SizeT __sz,
                           const _CharT* __s, _SizeT __pos, _SizeT __n)
   {
      if (__pos < __sz)
      {
         const _CharT* __pe = __p + __sz;
         for (const _CharT* __ps = __p + __pos; __ps != __pe; ++__ps)
            if (_Traits::find(__s, __n, *__ps) == 0)
               return static_cast<_SizeT>(__ps - __p);
      }
      return __npos;
   }


   template<class _CharT, class _SizeT, class _Traits, _SizeT __npos>
   _SizeT _LIBCPP_CONSTEXPR_AFTER_CXX11 _LIBCPP_INLINE_VISIBILITY
   __str_find_first_not_of(const _CharT *__p, _SizeT __sz,
                           _CharT __c, _SizeT __pos)
   {
      if (__pos < __sz)
      {
         const _CharT* __pe = __p + __sz;
         for (const _CharT* __ps = __p + __pos; __ps != __pe; ++__ps)
            if (!_Traits::eq(*__ps, __c))
               return static_cast<_SizeT>(__ps - __p);
      }
      return __npos;
   }


   // __str_find_last_not_of
   template<class _CharT, class _SizeT, class _Traits, _SizeT __npos>
   _SizeT _LIBCPP_CONSTEXPR_AFTER_CXX11 _LIBCPP_INLINE_VISIBILITY
   __str_find_last_not_of(const _CharT *__p, _SizeT __sz,
                          const _CharT* __s, _SizeT __pos, _SizeT __n)
   {
      if (__pos < __sz)
         ++__pos;
      else
         __pos = __sz;
      for (const _CharT* __ps = __p + __pos; __ps != __p;)
         if (_Traits::find(__s, __n, *--__ps) == 0)
            return static_cast<_SizeT>(__ps - __p);
      return __npos;
   }


   template<class _CharT, class _SizeT, class _Traits, _SizeT __npos>
   _SizeT _LIBCPP_CONSTEXPR_AFTER_CXX11 _LIBCPP_INLINE_VISIBILITY
   __str_find_last_not_of(const _CharT *__p, _SizeT __sz,
                          _CharT __c, _SizeT __pos)
   {
      if (__pos < __sz)
         ++__pos;
      else
         __pos = __sz;
      for (const _CharT* __ps = __p + __pos; __ps != __p;)
         if (!_Traits::eq(*--__ps, __c))
            return static_cast<_SizeT>(__ps - __p);
      return __npos;
   }

} // namespace __1 or __ROOT

} // namespace std

// Now include the main meat
#include "libcpp_string_view.h"


// And undo
#ifdef _LOCAL_LIBCPP_BEGIN_NAMESPACE_LFTS
#undef _LIBCPP_BEGIN_NAMESPACE_LFTS
#undef _LIBCPP_END_NAMESPACE_LFTS
#endif

#ifdef _LOCAL_LIBCPP_BEGIN_NAMESPACE_STD
#undef _LIBCPP_BEGIN_NAMESPACE_STD
#undef _LIBCPP_END_NAMESPACE_STD
#endif

#ifdef _LOCAL_LIBCPP_CONSTEXPR
#undef _LIBCPP_CONSTEXPR
#endif

#ifdef _LOCAL_VSTD
#undef _VSTD
#endif

#ifdef _LOCAL_LIBCPP_INLINE_VISIBILITY
#undef _LIBCPP_INLINE_VISIBILITY
#endif

#ifdef _LOCAL_LIBCPP_STD_VER
#undef _LIBCPP_STD_VER
#endif

#ifdef _LOCAL_LIBCPP_TYPE_VIS_ONLY
#undef _LIBCPP_TYPE_VIS_ONLY
#endif

#ifdef _LOCAL_LIBCPP_CONSTEXPR_AFTER_CXX11
#undef _LIBCPP_CONSTEXPR_AFTER_CXX11
#endif

#ifdef _LOCAL_NOEXCEPT
#undef _NOEXCEPT
#endif

#ifdef _LOCAL_LIBCPP_ASSERT
#undef _LIBCPP_ASSERT
#endif

#endif // RWrap_libcpp_string_view_h
