// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

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

double SqrtLowParameterTransformation::Int2ext(double value, double lower) const
{
   /// internal to external transformation
   double val = lower - 1. + std::sqrt(value * value + 1.);
   return val;
}

double SqrtLowParameterTransformation::Ext2int(double value, double lower, const MnMachinePrecision &) const
{
   // external to internal transformation
   double yy = value - lower + 1.;
   double yy2 = yy * yy;
   if (yy2 < 1.)
      return 0;
   else
      return std::sqrt(yy2 - 1);
}

double SqrtLowParameterTransformation::DInt2Ext(double value, double) const
{
   // derivative of internal to external transofrmation   :  d (Int2Ext) / d Int
   double val = value / (std::sqrt(value * value + 1.));
   return val;
}

} // namespace Minuit2

} // namespace ROOT
