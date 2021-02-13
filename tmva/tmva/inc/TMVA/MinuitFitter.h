// @(#)root/tmva $Id$ 
// Author: Andreas Hoecker

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : MinuitFitter                                                          *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *       Fitter using MINUIT                                                      *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker  <Andreas.Hocker@cern.ch> - CERN, Switzerland             *
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         * 
 *      MPI-K Heidelberg, Germany                                                 * 
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

#ifndef ROOT_TMVA_MinuitFitter
#define ROOT_TMVA_MinuitFitter

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// MinuitFitter                                                         //
//                                                                      //
// Fitter using MINUIT                                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TMVA/FitterBase.h"
#include "TMVA/IFitterTarget.h"
#include <vector>

class TFitter;

namespace TMVA {

   class IFitterTarget;
   class Interval;
   class MinuitWrapper;

   class MinuitFitter : public FitterBase, public IFitterTarget {

   public:

      MinuitFitter( IFitterTarget& target, const TString& name, 
                    std::vector<TMVA::Interval*>& ranges, const TString& theOption );

      virtual ~MinuitFitter();
      
      void Init();

      Double_t Run( std::vector<Double_t>& pars );
      Double_t EstimatorFunction( std::vector<Double_t>& pars );

   protected:

      MinuitWrapper *fMinWrap; // holds a wrapper around TMinuit

   private:

      void DeclareOptions();

      Int_t      fErrorLevel;              // minuit error level
      Int_t      fPrintLevel;              // minuit printout level
      Int_t      fFitStrategy;             // minuit strategy level
      Bool_t     fPrintWarnings;           // minuit warnings level
      Bool_t     fUseImprove;              // flag for 'IMPROVE' use
      Bool_t     fUseMinos;                // flag for 'MINOS' use
      Bool_t     fBatch;                   // batch mode
      Int_t      fMaxCalls;                // (approximate) maximum number of function calls
      Double_t   fTolerance;               // tolerance to the function value at the minimum
      
      ClassDef(MinuitFitter,0); // Fitter using a Genetic Algorithm
   };

} // namespace TMVA

#endif


