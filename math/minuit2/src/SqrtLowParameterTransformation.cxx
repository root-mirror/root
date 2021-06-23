// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei, E.G.P. Bos   2003-2017

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

// Project   : LCG
// Package   : Minuit
// Author    : Lorenzo.MONETA@cern.ch
// Created by: moneta  at Thu Apr  8 10:26:22 2004

#include "Minuit2/SqrtLowParameterTransformation.h"
#include "Minuit2/MnMachinePrecision.h"

namespace ROOT {

  namespace Minuit2 {



    long double SqrtLowParameterTransformation::Int2ext(long double value, long double lower) const {
      /// internal to external transformation
      long double val = lower - 1. + std::sqrt( value*value + 1.);
      return val;
    }


    long double SqrtLowParameterTransformation::Ext2int(long double value, long double lower, const MnMachinePrecision& ) const {
      // external to internal transformation
      long double yy = value - lower + 1.;
      long double yy2 = yy*yy;
      if (yy2 < 1. )
        return 0;
      else
        return std::sqrt( yy2 -1);
    }


    long double SqrtLowParameterTransformation::DInt2Ext(long double value, long double) const {
      // derivative of internal to external transofrmation   :  d (Int2Ext) / d Int
      long double val = value/( std::sqrt( value*value + 1.) );
      return val;
    }

  }  // namespace Minuit2

} // namespace ROOT
