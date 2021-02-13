// @(#)root/tmva $Id$ 
// Author: Andreas Hoecker, Peter Speckmayer, Joerg Stelzer, Helge Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : MCFitter                                                              *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Fitter using Monte Carlo sampling of parameters                           *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker  <Andreas.Hocker@cern.ch> - CERN, Switzerland             *
 *      Peter Speckmayer <speckmay@mail.cern.ch>  - CERN, Switzerland             *
 *      Joerg Stelzer    <Joerg.Stelzer@cern.ch>  - CERN, Switzerland             *
 *      Helge Voss       <Helge.Voss@cern.ch>     - MPI-K Heidelberg, Germany     *
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         * 
 *      MPI-K Heidelberg, Germany                                                 * 
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

#ifndef ROOT_TMVA_MCFitter
#define ROOT_TMVA_MCFitter

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// MCFitter                                                             //
//                                                                      //
// Fitter using Monte Carlo sampling of parameters                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TMVA/FitterBase.h"
#include <vector>

namespace TMVA {

   class MCFitter : public FitterBase {

   public:

      MCFitter( IFitterTarget& target, const TString& name, 
                const std::vector<TMVA::Interval*>& ranges, const TString& theOption );

      virtual ~MCFitter() {}

      void SetParameters( Int_t cycles );

      Double_t Run( std::vector<Double_t>& pars );

   private:

      void DeclareOptions();

      Int_t    fSamples;     // number of MC samples
      Double_t fSigma;       // new samples are generated randomly with a gaussian probability with fSigma around the current best value
      UInt_t   fSeed;        // Seed for the random generator (0 takes random seeds)

      ClassDef(MCFitter,0); //  Fitter using Monte Carlo sampling of parameters 
   };

} // namespace TMVA

#endif


