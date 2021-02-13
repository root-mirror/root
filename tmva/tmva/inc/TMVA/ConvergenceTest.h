// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Joerg Stelzer, Helge Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : ConvergenceTest                                                             *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Contains all the data information                                         *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
 *      Peter Speckmayer <Peter.Speckmayer@cern.ch>  - CERN, Switzerland          *
 *      Joerg Stelzer   <Joerg.Stelzer@cern.ch>  - CERN, Switzerland              *
 *      Helge Voss      <Helge.Voss@cern.ch>     - MPI-K Heidelberg, Germany      *
 *                                                                                *
 * Copyright (c) 2006:                                                            *
 *      CERN, Switzerland                                                         *
 *      U. of Victoria, Canada                                                    *
 *      MPI-K Heidelberg, Germany                                                 *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

#ifndef ROOT_TMVA_ConvergenceTest
#define ROOT_TMVA_ConvergenceTest

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// ConvergenceTest                                                      //
//                                                                      //
// check for convergence                                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <deque>

#include "RtypesCore.h"

namespace TMVA {
   
   class ConvergenceTest {

   public:
      
      ConvergenceTest();
      ~ConvergenceTest();

      // setters
      void                       SetConvergenceParameters(  Int_t steps, Double_t improvement )  
      { fSteps = steps; fImprovement = improvement; }
      void                       SetCurrentValue(  Float_t value )  { fCurrentValue = value; }
      Float_t                    GetCurrentValue()                  { return fCurrentValue; }
      void                       ResetConvergenceCounter()  { fCounter = -1; fMaxCounter = 0; }

      // getters
      Bool_t                     HasConverged( Bool_t withinConvergenceBand = kFALSE );
      Float_t                    Progress();          // from 0 (just started) to 1 (finished)
      Float_t                    SpeedControl( UInt_t ofSteps );  


   protected:

      Float_t                    fCurrentValue;      //! current value

      Float_t                    fImprovement;       //! minimum improvement which counts as improvement
      Int_t                      fSteps;             //! number of steps without improvement required for convergence

   private:
      
      Int_t                      fCounter;           //! counts the number of steps without improvement 
      Float_t                    fConvValue;         //! the best "fitness" value
      Int_t                      fMaxCounter;        //! maximum value for the counter so far 

      // speed-control (gives back the learning speed = improvement-rate in the last N steps)
      // successList keeps track of the improvements to be able
      Float_t                    fBestResult;        // 
      Float_t                    fLastResult;        // 
      std::deque<Short_t>        fSuccessList;       // to calculate the improvement-speed

   };
}

#endif
