// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Joerg Stelzer, Helge Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : RootFinder                                                            *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Implementation (see header file for description)                          *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
 *      Helge Voss      <Helge.Voss@cern.ch>     - MPI-K Heidelberg, Germany      *
 *      Kai Voss        <Kai.Voss@cern.ch>       - U. of Victoria, Canada         *
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         *
 *      U. of Victoria, Canada                                                    *
 *      MPI-K Heidelberg, Germany                                                 *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

/*! \class TMVA::RootFinder
\ingroup TMVA
Root finding using Brents algorithm (translated from CERNLIB function RZERO)
*/

#include "TMVA/RootFinder.h"

#include "TMVA/MethodBase.h"
#include "TMVA/MsgLogger.h"
#include "TMVA/Types.h"

#include "TMath.h"

ClassImp(TMVA::RootFinder);

////////////////////////////////////////////////////////////////////////////////
/// constructor

TMVA::RootFinder::RootFinder(TMVA::MethodBase *method ,
                              Double_t rootMin,
                              Double_t rootMax,
                              Int_t maxIterations,
                              Double_t absTolerance )
: fRootMin( rootMin ),
   fRootMax( rootMax ),
   fMaxIter( maxIterations ),
   fAbsTol ( absTolerance  ),
   fLogger ( new MsgLogger("RootFinder") )
{
  fMethod=method;
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TMVA::RootFinder::~RootFinder( void )
{
   delete fLogger;
}

////////////////////////////////////////////////////////////////////////////////
/// Root finding using Brents algorithm; taken from CERNLIB function RZERO

Double_t TMVA::RootFinder::Root( Double_t refValue  )
{
   Double_t a  = fRootMin, b = fRootMax;
   Double_t fa = fMethod->GetValueForRoot( a ) - refValue;
   Double_t fb = fMethod->GetValueForRoot( b ) - refValue;
   if (fb*fa > 0) {
      Log() << kWARNING << "<Root> initial interval w/o root: "
            << "(a=" << a << ", b=" << b << "),"
            << " (Eff_a=" << fMethod->GetValueForRoot( a )
            << ", Eff_b=" << fMethod->GetValueForRoot( b ) << "), "
            << "(fa=" << fa << ", fb=" << fb << "), "
            << "refValue = " << refValue << Endl;
      return 1;
   }

   Bool_t   ac_equal(kFALSE);
   Double_t fc = fb;
   Double_t c  = 0, d = 0, e = 0;
   for (Int_t iter= 0; iter <= fMaxIter; iter++) {
      if ((fb < 0 && fc < 0) || (fb > 0 && fc > 0)) {

         // Rename a,b,c and adjust bounding interval d
         ac_equal = kTRUE;
         c  = a; fc = fa;
         d  = b - a; e  = b - a;
      }

      if (TMath::Abs(fc) < TMath::Abs(fb)) {
         ac_equal = kTRUE;
         a  = b;  b  = c;  c  = a;
         fa = fb; fb = fc; fc = fa;
      }

      Double_t tol = 0.5 * 2.2204460492503131e-16 * TMath::Abs(b);
      Double_t m   = 0.5 * (c - b);
      if (fb == 0 || TMath::Abs(m) <= tol || TMath::Abs(fb) < fAbsTol) return b;

      // Bounds decreasing too slowly: use bisection
      if (TMath::Abs (e) < tol || TMath::Abs (fa) <= TMath::Abs (fb)) { d = m; e = m; }
      else {
         // Attempt inverse cubic interpolation
         Double_t p, q, r;
         Double_t s = fb / fa;

         if (ac_equal) { p = 2 * m * s; q = 1 - s; }
         else {
            q = fa / fc; r = fb / fc;
            p = s * (2 * m * q * (q - r) - (b - a) * (r - 1));
            q = (q - 1) * (r - 1) * (s - 1);
         }
         // Check whether we are in bounds
         if (p > 0) q = -q;
         else       p = -p;

         Double_t min1 = 3 * m * q - TMath::Abs (tol * q);
         Double_t min2 = TMath::Abs (e * q);
         if (2 * p < (min1 < min2 ? min1 : min2)) {
            // Accept the interpolation
            e = d;        d = p / q;
         }
         else { d = m; e = m; } // Interpolation failed: use bisection.
      }
      // Move last best guess to a
      a  = b; fa = fb;
      // Evaluate new trial root
      if (TMath::Abs(d) > tol) b += d;
      else                     b += (m > 0 ? +tol : -tol);

      fb = fMethod->GetValueForRoot( b ) - refValue;

   }

   // Return our best guess if we run out of iterations
   Log() << kWARNING << "<Root> maximum iterations (" << fMaxIter
         << ") reached before convergence" << Endl;

   return b;
}
