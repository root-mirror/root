/// \file ROOT/RIndexIter.hxx
/// \ingroup Base ROOT7
/// \author Axel Naumann <axel@cern.ch>
/// \date 2016-01-19
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2016, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RIndexIter
#define ROOT7_RIndexIter

#include <iterator>

namespace ROOT {
namespace Experimental {
namespace Internal {

/**
 \class RIndexIter
 Iterates over an index; the REFERENCE is defined by the REFERENCE template parameter.

 Derived classes are expected to implement `const REFERENCE& operator*()` and
 `const POINTER operator->()`.
 */

template <class REFERENCE,
          class POINTER = typename std::add_pointer<typename std::remove_reference<REFERENCE>::type>::type>
class RIndexIter: public std::iterator<std::random_access_iterator_tag, REFERENCE, POINTER> {
   size_t fIndex;

protected:
   static constexpr size_t fgEndIndex = (size_t)-1;

public:
   RIndexIter(size_t idx): fIndex(idx) {}

   /// Get the current index value.
   size_t GetIndex() const noexcept { return fIndex; }

   ///\{
   ///\name Index modifiers
   /// ++i
   RIndexIter &operator++() noexcept
   {
      ++fIndex;
      return *this;
   }

   /// --i
   RIndexIter &operator--() noexcept
   {
      if (fIndex != fgEndIndex)
         --fIndex;
      return *this;
   }

   /// i++
   RIndexIter operator++(int)noexcept
   {
      RIndexIter old(*this);
      ++(*this);
      return old;
   }

   // i--
   RIndexIter operator--(int)noexcept
   {
      RIndexIter old(*this);
      --(*this);
      return old;
   }

   RIndexIter &operator+=(int d) noexcept
   {
      fIndex += d;
      return *this;
   }

   RIndexIter &operator-=(int d) noexcept
   {
      if (d > fIndex) {
         fIndex = fgEndIndex;
      } else {
         fIndex -= d;
      }
      return *this;
   }

   RIndexIter operator+(int d) noexcept
   {
      RIndexIter ret(*this);
      ret += d;
      return ret;
   }

   RIndexIter operator-(int d) noexcept
   {
      RIndexIter ret(*this);
      ret -= d;
      return ret;
   }
   ///\}
};

///\{
///\name Relational operators.
template <class REFERENCE, class POINTER>
bool operator<(RIndexIter<REFERENCE, POINTER> lhs, RIndexIter<REFERENCE, POINTER> rhs) noexcept
{
   return lhs.GetIndex() < rhs.GetIndex();
}

template <class REFERENCE, class POINTER>
bool operator>(RIndexIter<REFERENCE, POINTER> lhs, RIndexIter<REFERENCE, POINTER> rhs) noexcept
{
   return lhs.GetIndex() > rhs.GetIndex();
}

template <class REFERENCE, class POINTER>
bool operator<=(RIndexIter<REFERENCE, POINTER> lhs, RIndexIter<REFERENCE, POINTER> rhs) noexcept
{
   return lhs.GetIndex() <= rhs.GetIndex();
}

template <class REFERENCE, class POINTER>
inline bool operator>=(RIndexIter<REFERENCE, POINTER> lhs, RIndexIter<REFERENCE, POINTER> rhs) noexcept
{
   return lhs.GetIndex() >= rhs.GetIndex();
}

template <class REFERENCE, class POINTER>
inline bool operator==(RIndexIter<REFERENCE, POINTER> lhs, RIndexIter<REFERENCE, POINTER> rhs) noexcept
{
   return lhs.GetIndex() == rhs.GetIndex();
}

template <class REFERENCE, class POINTER>
inline bool operator!=(RIndexIter<REFERENCE, POINTER> lhs, RIndexIter<REFERENCE, POINTER> rhs) noexcept
{
   return lhs.GetIndex() != rhs.GetIndex();
}
///\}

} // namespace Internal
} // namespace Experimental
} // namespace ROOT

#endif // ROOT7_RIndexIter
