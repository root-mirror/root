// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_MnMatrix
#define ROOT_Minuit2_MnMatrix

// add MnConfig file to define before everything compiler
// dependent macros

#include "Minuit2/MnConfig.h"

#include "Minuit2/LASymMatrix.h"
#include "Minuit2/LAVector.h"
#include "Minuit2/LaInverse.h"
#include "Minuit2/LaOuterProduct.h"

namespace ROOT {

namespace Minuit2 {

typedef LASymMatrix MnAlgebraicSymMatrix;
typedef LAVector MnAlgebraicVector;

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_MnMatrix
