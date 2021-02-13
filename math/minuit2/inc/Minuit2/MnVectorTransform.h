// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_MnVectorTransform
#define ROOT_Minuit2_MnVectorTransform

#include "Minuit2/MnMatrix.h"

#include <vector>

namespace ROOT {

namespace Minuit2 {

class MnVectorTransform {

public:
   MnVectorTransform() {}

   ~MnVectorTransform() {}

   std::vector<double> operator()(const MnAlgebraicVector &avec) const
   {

      std::vector<double> result;
      result.reserve(avec.size());

      for (unsigned int i = 0; i < avec.size(); i++)
         result.push_back(avec(i));

      return result;
   }
};

} // namespace Minuit2

} // namespace ROOT

#endif // ROOT_Minuit2_MnVectorTransform
