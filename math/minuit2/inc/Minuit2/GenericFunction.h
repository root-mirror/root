// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_GenericFunction
#define ROOT_Minuit2_GenericFunction

#include "Minuit2/MnConfig.h"

#include <vector>

namespace ROOT {

namespace Minuit2 {

//_____________________________________________________________________
/**

Class from which all the other classes, representing functions,
inherit. That is why it defines only one method, the operator(),
which allows to call the function.

@author Andras Zsenei and Lorenzo Moneta, Creation date: 23 Sep 2004

@ingroup Minuit

 */

class GenericFunction {

public:
   virtual ~GenericFunction() {}

   /**

      Evaluates the function using the vector containing the input values.

      @param x vector of the coordinates (for example the x coordinate for a
      one-dimensional Gaussian)

      @return the result of the evaluation of the function.

   */

   virtual double operator()(const std::vector<double> &x) const = 0;
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_GenericFunction
