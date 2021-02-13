/* @(#)root/hist:$Id$ */

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_Foption
#define ROOT_Foption

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Foption                                                              //
//                                                                      //
// Histogram fit options structure.                                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "ROOT/EExecutionPolicy.hxx"

struct Foption_t {
//*-*      chopt may be the concatenation of the following options:
//*-*      =======================================================
//*-*
//*-*   The following structure members are set to 1 if the option is selected:
   int Quiet;       // "Q"  Quiet mode. No print
   int Verbose;     // "V"  Verbose mode. Print results after each iteration
   int Bound;       // "B"  When using pre-defined functions user parameter settings are used instead of default one
   int Chi2;        // "X"  For fitting THnsparse use chi2 method (default is likelihood)
   int PChi2;       // "P"  Use Pearson chi2 built with the expected error instead of the observed ones
   int Like;        // "L"  Use Log Likelihood. Default is chisquare method except fitting THnsparse
   int User;        // "U"  Use a User specified fitting algorithm (via SetFCN)
   int W1;          // "W"  Set all the weights to 1. Ignore error bars
   int Errors;      // "E"  Performs a better error evaluation, calling HESSE and MINOS
   int More;        // "M"  Improve fit results.
   int Range;       // "R"  Use the range stored in function
   int Gradient;    // "G"  Option to compute derivatives analytically
   int Nostore;     // "N"  If set, do not store the function graph
   int Nograph;     // "0"  If set, do not display the function graph
   int Plus;        // "+"  Add new function (default is replace)
   int Integral;    // "I"  Use function integral instead of function in center of bin
   int Nochisq;     // "C"  In case of linear fitting, don't calculate the chisquare
   int Minuit;      // "F"  If fitting a polN, switch to minuit fitter
   int NoErrX;      // "EX0" or "T" When fitting a TGraphErrors do not consider error in coordinates
   int Robust;      // "ROB" or "H":  For a TGraph use robust fitting
   int StoreResult; // "S": Stores the result in a TFitResult structure
   int BinVolume;   // "WIDTH": scale content by the bin width/volume
   double hRobust;  //  value of h parameter used in robust fitting
   ROOT::EExecutionPolicy ExecPolicy;  //  Choose the execution Policy: "SERIAL", "MULTITHREAD" or "MULTIPROCESS"

  Foption_t() :
      Quiet        (0),
      Verbose      (0),
      Bound        (0),
      Chi2         (0),
      PChi2        (0),
      Like         (0),
      User         (0),
      W1           (0),
      Errors       (0),
      More         (0),
      Range        (0),
      Gradient     (0),
      Nostore      (0),
      Nograph      (0),
      Plus         (0),
      Integral     (0),
      Nochisq      (0),
      Minuit       (0),
      NoErrX       (0),
      Robust       (0),
      StoreResult  (0),
      BinVolume    (0),
      hRobust      (0),
      ExecPolicy   (ROOT::EExecutionPolicy::kSequential)
   {}
};

#endif
