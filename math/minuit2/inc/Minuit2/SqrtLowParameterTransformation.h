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

#ifndef ROOT_Minuit2_SqrtLowParameterTransformation
#define ROOT_Minuit2_SqrtLowParameterTransformation

namespace ROOT {

namespace Minuit2 {

class MnMachinePrecision;

/**
 * Transformation from external to internal Parameter based on  sqrt(1 + x**2)
 *
 * This transformation applies for the case of single side Lower Parameter limits
 */

class SqrtLowParameterTransformation /* : public ParameterTransformation */ {

public:
   SqrtLowParameterTransformation() {}

   ~SqrtLowParameterTransformation() {}

   // transformation from internal to external
   double Int2ext(double Value, double Lower) const;

   // transformation from external to internal
   double Ext2int(double Value, double Lower, const MnMachinePrecision &) const;

   // derivative of transformation from internal to external
   double DInt2Ext(double Value, double Lower) const;

private:
};

} // namespace Minuit2

} // namespace ROOT

#endif
