// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "Minuit2/MnTiny.h"

namespace ROOT {

namespace Minuit2 {

double MnTiny::One() const
{
   return fOne;
}

double MnTiny::operator()(volatile double epsp1) const
{
   // evaluate minimal diference between two floating points
   double result = epsp1 - One();
   return result;
}

} // namespace Minuit2

} // namespace ROOT
