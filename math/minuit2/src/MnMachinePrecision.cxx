// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "Minuit2/MnMachinePrecision.h"
#include "Minuit2/MnTiny.h"

namespace ROOT {

   namespace Minuit2 {

MnStaticMachinePrecision::MnStaticMachinePrecision() :
   eps(4.0E-7),
   eps2(2.*sqrt(4.0E-7)) {

   //determine machine precision
   /*
       char e[] = {"e"};
       eps = 8.*dlamch_(e);
       eps2 = 2.*sqrt(eps);
   */

   //   std::cout<<"machine precision eps: "<<eps()<<std::endl;

   MnTiny mytiny;

   //calculate machine precision
   double epstry = 0.5;
   double epsbak = 0.;
   volatile double epsp1 = 0.; // allow to run this method with fast-math
   double one = 1.0;
   for(int i = 0; i < 100; i++) {
      epstry *= 0.5;
      epsp1 = one + epstry;
      epsbak = mytiny(epsp1);
      if(epsbak < epstry) {
         eps = 8.*epstry;
         eps2 = 2.*sqrt(eps);
         break;
      }
   }
}

   }  // namespace Minuit2

}  // namespace ROOT
