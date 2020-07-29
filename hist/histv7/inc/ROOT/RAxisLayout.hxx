/// \file ROOT/RAxisLayout.h
/// \ingroup Hist ROOT7
/// \author Axel Naumann <axel@cern.ch>
/// \date 2020-07-06
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RAxisLayout
#define ROOT7_RAxisLayout

#include <string>
#include <utility>
#include <vector>

namespace ROOT {
namespace Experimental {
namespace Internal {

/**
\class RAxisLayout
Calculuates bin index from a set of axes and a coordinate tuple.
*/
template <class AxisTuple>
class RAxisLayout {
private:
   ///\{
   /// Internal types to translate tuple<AXES...> => tuple<AXES::Coord_t...>
   template <class T>
   struct CoordTuple;

   template <class... Axes>
   struct CoordTuple<std::tuple<Axes...>> {
      using type = std::tuple<typename Axes::Coord_t...>;
   };
   ///\}

public:
   using AxisTuple_t = AxisTuple;
   using CoordTuple_t = typename CoordTuple<AxisTuple>::type;

private:
   AxisTuple fAxes;

public:
   RAxisLayout(const AxisTuple &axes): fAxes(axes) {}

   ssize_t GetBinIndex(const CoordTuple_t &x) const;
};

} // namespace Internal
} // namespace Experimental
} // namespace ROOT

#endif // ROOT7_RAxisLayout header guard