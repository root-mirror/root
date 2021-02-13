/// \file ROOT/RHistUtils.hxx
/// \ingroup Hist ROOT7
/// \author Axel Naumann <axel@cern.ch>
/// \date 2016-06-01
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!

/*************************************************************************
 * Copyright (C) 1995-2016, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RHistUtils
#define ROOT7_RHistUtils

#include <array>
#include <type_traits>

namespace ROOT {
namespace Experimental {

class RLogChannel;
/// Log channel for Hist diagnostics.
RLogChannel &HistLog(); // implemented in RAxis.cxx

namespace Hist {

template <int DIMENSIONS>
struct RCoordArray: std::array<double, DIMENSIONS> {
   using Base_t = std::array<double, DIMENSIONS>;

   /// Default construction.
   RCoordArray() = default;

   /// Construction with one `double` per `DIMENSION`.
   template<class...ELEMENTS, class = typename std::enable_if<sizeof...(ELEMENTS) + 1 == DIMENSIONS>::type>
   RCoordArray(double x, ELEMENTS...el): Base_t{{x, el...}} {}

   /// Fallback constructor, invoked if the one above fails because of the wrong number of
   /// arguments / coordinates.
   template<class T, class...ELEMENTS, class = typename std::enable_if<sizeof...(ELEMENTS) + 1 != DIMENSIONS>::type>
   RCoordArray(T, ELEMENTS...) {
      static_assert(sizeof...(ELEMENTS) + 1 == DIMENSIONS, "Number of coordinates does not match DIMENSIONS");
   }

   /// Construction from a C-style array.
   RCoordArray(double (&arr)[DIMENSIONS]): Base_t(arr) {}

   /// Copy-construction from a C++-style array.
   /// (No need for a move-constructor, it isn't any better for doubles)
   RCoordArray(const std::array<double, DIMENSIONS>& arr): Base_t(arr) {}
};

template <int DIMENSIONS>
//using CoordArray_t = std::array<double, DIMENSIONS>;
using CoordArray_t = RCoordArray<DIMENSIONS>;


} // namespace Hist
} // namespace Experimental
} // namespace ROOT

#endif //ROOT7_THistUtils_h
