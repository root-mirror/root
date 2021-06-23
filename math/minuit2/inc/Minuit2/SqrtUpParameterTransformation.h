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

#ifndef ROOT_Minuit2_SqrtUpParameterTransformation
#define ROOT_Minuit2_SqrtUpParameterTransformation

namespace ROOT {

namespace Minuit2 {

    class MnMachinePrecision;

/**
 * Transformation from external to internal Parameter based on  sqrt(1 + x**2)
 *
 * This transformation applies for the case of single side Upper Parameter limits
 */

    class SqrtUpParameterTransformation /* : public ParameterTransformation */ {

    public:

      // create with user defined precision
      SqrtUpParameterTransformation() {}

      ~SqrtUpParameterTransformation() {}

      // transformation from internal to external
      long double Int2ext(long double Value, long double Upper) const;

      // transformation from external to internal
      long double Ext2int(long double Value, long double Upper, const MnMachinePrecision&) const;

      // derivative of transformation from internal to external
      long double DInt2Ext(long double Value, long double Upper) const;

    private:
    };

} // namespace Minuit2

} // namespace ROOT

#endif
