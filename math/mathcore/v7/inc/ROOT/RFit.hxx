/// \file ROOT/RFit.hxx
/// \ingroup MathCore ROOT7
/// \author Axel Naumann <axel@cern.ch>
/// \date 2015-09-06
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2015, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_RFit
#define ROOT_RFit

#include <array>
#include <functional>

#include "ROOT/RSpan.hxx"

#include "ROOT/RHist.hxx"

namespace ROOT {
namespace Experimental {

class RFitResult {
};

template <int DIMENSION>
class RFunction {
public:
   RFunction(std::function<double(const std::array<double, DIMENSION> &, const std::span<const double> par)> func) {}
};

template <int DIMENSIONS, class PRECISION, template <int D_, class P_> class... STAT>
RFitResult FitTo(const RHist<DIMENSIONS, PRECISION, STAT...> &hist, const RFunction<DIMENSIONS> &func,
                 std::span<const double> paramInit)
{
   return RFitResult();
}

} // namespace Experimental
} // namespace ROOT

#endif
