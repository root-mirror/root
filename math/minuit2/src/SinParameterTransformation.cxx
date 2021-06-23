// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei, E.G.P. Bos   2003-2017

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "Minuit2/SinParameterTransformation.h"
#include "Minuit2/MnMachinePrecision.h"

#include <cmath>

namespace ROOT {

  namespace Minuit2 {



    long double SinParameterTransformation::Int2ext(long double Value, long double Upper, long double Lower) const {
      // transformation from  to internal (unlimited) to external values (limited by Lower/Upper )
      return Lower + 0.5*(Upper - Lower)*(std::sin(Value) + 1.);
    }

    long double SinParameterTransformation::Ext2int(long double Value, long double Upper, long double Lower, const MnMachinePrecision& prec) const {
      // transformation from external (limited by Lower/Upper )  to internal (unlimited) values given the lower/upper limits

      long double piby2 = 2.*std::atan(1.);
      long double distnn = 8.*std::sqrt(prec.Eps2());
      long double vlimhi = piby2 - distnn;
      long double vlimlo = -piby2 + distnn;

      long double yy = 2.*(Value - Lower)/(Upper - Lower) - 1.;
      long double yy2 = yy*yy;
      if(yy2 > (1. - prec.Eps2())) {
        if(yy < 0.) {
          // Lower limit
          //       std::cout<<"SinParameterTransformation warning: is at its Lower allowed limit. "<<Value<<std::endl;
          return vlimlo;
        } else {
          // Upper limit
          //       std::cout<<"SinParameterTransformation warning: is at its Upper allowed limit."<<std::endl;
          return vlimhi;
        }

      } else {
        return std::asin(yy);
      }
    }

    long double SinParameterTransformation::DInt2Ext(long double Value, long double Upper, long double Lower) const {
      // return the derivative of the transformation d Ext/ d Int
      return 0.5*((Upper - Lower)*std::cos(Value));
    }

  }  // namespace Minuit2

} // namespace ROOT
