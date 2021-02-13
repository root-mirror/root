// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Krzysztof Danielowski, Kamil Kraszewski, Maciej Kruk

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : SimulatedAnnealingFitter                                              *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *       Fitter using Simulated Annealing algorithm                               *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Krzysztof Danielowski <danielow@cern.ch>       - IFJ & AGH, Poland        *
 *      Andreas Hoecker       <Andreas.Hocker@cern.ch> - CERN, Switzerland        *
 *      Kamil Kraszewski      <kalq@cern.ch>           - IFJ & UJ, Poland         *
 *      Maciej Kruk           <mkruk@cern.ch>          - IFJ & AGH, Poland        *
 *                                                                                *
 * Copyright (c) 2008:                                                            *
 *      IFJ-Krakow, Poland                                                        *
 *      CERN, Switzerland                                                         *
 *      MPI-K Heidelberg, Germany                                                 *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

#ifndef ROOT_TMVA_SimulatedAnnealingFitter
#define ROOT_TMVA_SimulatedAnnealingFitter

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// SimulatedAnnealingFitter                                             //
//                                                                      //
// Fitter using a Simulated Annealing Algorithm                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TMVA/FitterBase.h"

#include <vector>

namespace TMVA {

   class IFitterTarget;
   class Interval;

   class SimulatedAnnealingFitter : public FitterBase {

   public:

      SimulatedAnnealingFitter( IFitterTarget& target, const TString& name,
                                const std::vector<TMVA::Interval*>& ranges, const TString& theOption );

      virtual ~SimulatedAnnealingFitter() {}

      void SetParameters( Int_t    fMaxCalls,
                          Double_t fInitialTemperature,
                          Double_t fMinTemperature,
                          Double_t fEps,
                          TString  fKernelTemperatureS,
                          Double_t fTemperatureScale,
                          Double_t fTemperatureAdaptiveStep,
                          Bool_t   fUseDefaultScale,
                          Bool_t   fUseDefaultTemperature );

      Double_t Run( std::vector<Double_t>& pars );

   private:

      void DeclareOptions();

      Int_t              fMaxCalls;                // max number of FCN calls
      Double_t           fInitialTemperature;      // initial temperature (depends on FCN)
      Double_t           fMinTemperature;          // minimum temperature before SA quit
      Double_t           fEps;                     // relative required FCN accuracy at minimum
      TString            fKernelTemperatureS;      // string just to set fKernelTemperature
      Double_t           fTemperatureScale;        // how fast temperature change
      Double_t           fAdaptiveSpeed;           // how fast temperature change in adaptive (in adaptive two variables describe
                                                   // the change of temperature, but fAdaptiveSpeed should be 1.0 and its not
                                                   // recommended to change it)
      Double_t           fTemperatureAdaptiveStep; // used to calculate InitialTemperature if fUseDefaultTemperature
      Bool_t             fUseDefaultScale;         // if TRUE, SA calculates its own TemperatureScale
      Bool_t             fUseDefaultTemperature;   // if TRUE, SA calculates its own InitialTemperature (MinTemperautre)

      ClassDef(SimulatedAnnealingFitter,0); // Fitter using a Simulated Annealing Algorithm
   };

} // namespace TMVA

#endif


