/// \file ROOT/RHistData.h
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
namespace Hist {

namespace Stat {
enum EStat { // not enum class - want int-casts!
   kUncertainty = 1, ///< Poisson uncertainty per bin
   kRuntime = 2, ///< Runtime statistics - but how to set them?
   k1stMoment = 4, ///< 1st moment
   k2ndMoment = 8, ///< 2nd moment
   // k3rdMoment = 16, ///< 3rd moment
   // k4thMoment = 32, ///< 4th moment
};
}

} // namespace Hist
} // namespace Experimental
} // namespace ROOT

#endif //ROOT7_THistUtils_h
