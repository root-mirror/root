// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_MnContours
#define ROOT_Minuit2_MnContours

#include "Minuit2/MnConfig.h"
#include "Minuit2/MnStrategy.h"

#include <vector>
#include <utility>

namespace ROOT {

namespace Minuit2 {

class FCNBase;
class FunctionMinimum;
class ContoursError;

//_____________________________________________________________
/**
   API class for Contours Error analysis (2-dim errors);
   minimization has to be done before and Minimum must be valid;
   possibility to ask only for the points or the points and associated Minos
   errors;
 */

class MnContours {

public:
   /// construct from FCN + Minimum
   MnContours(const FCNBase &fcn, const FunctionMinimum &min) : fFCN(fcn), fMinimum(min), fStrategy(MnStrategy(1)) {}

   /// construct from FCN + Minimum + strategy
   MnContours(const FCNBase &fcn, const FunctionMinimum &min, unsigned int stra)
      : fFCN(fcn), fMinimum(min), fStrategy(MnStrategy(stra))
   {
   }

   /// construct from FCN + Minimum + strategy
   MnContours(const FCNBase &fcn, const FunctionMinimum &min, const MnStrategy &stra)
      : fFCN(fcn), fMinimum(min), fStrategy(stra)
   {
   }

   ~MnContours() {}

   /// ask for one Contour (points only) from number of points (>=4) and parameter indeces
   std::vector<std::pair<double, double>> operator()(unsigned int, unsigned int, unsigned int npoints = 20) const;

   /// ask for one Contour ContoursError (MinosErrors + points)
   /// from number of points (>=4) and parameter indeces
   /// can be printed via std::cout
   ContoursError Contour(unsigned int, unsigned int, unsigned int npoints = 20) const;

   const MnStrategy &Strategy() const { return fStrategy; }

private:
   const FCNBase &fFCN;
   const FunctionMinimum &fMinimum;
   MnStrategy fStrategy;
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_MnContours
