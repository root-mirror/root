// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Joerg Stelzer, Helge Voss, Krzysztof Danielowski, Kamil Kraszewski, Maciej Kruk

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : SimulatedAnnealing                                                    *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Implementation of simulated annealing fitting procedure                   *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Krzysztof Danielowski <danielow@cern.ch>       - IFJ & AGH, Poland        *
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

#ifndef ROOT_TMVA_SimulatedAnnealing
#define ROOT_TMVA_SimulatedAnnealing

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// SimulatedAnnealing                                                   //
//                                                                      //
// Base implementation of simulated annealing fitting procedure         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <vector>

#include "TMVA/Types.h"

class TRandom;

namespace TMVA {

   class IFitterTarget;
   class Interval;
   class MsgLogger;

   class SimulatedAnnealing {

   public:

      SimulatedAnnealing( IFitterTarget& target, const std::vector<TMVA::Interval*>& ranges );
      virtual ~SimulatedAnnealing();

      // returns FCN value at minimum
      Double_t Minimize( std::vector<Double_t>& parameters );

      // accessors
      void SetMaxCalls          ( Int_t    mc    ) { fMaxCalls = mc; }
      void SetInitTemp          ( Double_t it    ) { fInitialTemperature = it; }
      void SetMinTemp           ( Double_t min   ) { fMinTemperature = min; }
      void SetAccuracy          ( Double_t eps   ) { fEps = eps; }
      void SetTemperatureScale  ( Double_t scale ) { fTemperatureScale = scale; }
      void SetAdaptiveSpeed     ( Double_t speed ) { fAdaptiveSpeed = speed; }

      void SetOptions( Int_t maxCalls, Double_t initialTemperature, Double_t minTemperature, Double_t eps,
                       TString  kernelTemperatureS, Double_t temperatureScale, Double_t adaptiveSpeed,
                       Double_t temperatureAdaptiveStep, Bool_t useDefaultScale, Bool_t useDefaultTemperature );

      //setting up helper variables for JsMVA
      void SetIPythonInteractive(bool* ExitFromTraining, UInt_t *fIPyCurrentIter_){
        fExitFromTraining = ExitFromTraining;
        fIPyCurrentIter = fIPyCurrentIter_;
      }

   private:

      enum EKernelTemperature {
         kSqrt = 0,
         kIncreasingAdaptive,
         kDecreasingAdaptive,
         kLog,
         kHomo,
         kSin,
         kGeo
      } fKernelTemperature;

      void FillWithRandomValues( std::vector<Double_t>& parameters );
      void ReWriteParameters( std::vector<Double_t>& from, std::vector<Double_t>& to );
      void GenerateNewTemperature(Double_t& currentTemperature, Int_t Iter );
      void GenerateNeighbour( std::vector<Double_t>& parameters, std::vector<Double_t>& oldParameters, Double_t currentTemperature );
      Bool_t ShouldGoIn( Double_t currentFit, Double_t localFit, Double_t currentTemperature );
      void SetDefaultScale();
      Double_t GenerateMaxTemperature( std::vector<Double_t>& parameters );
      std::vector<Double_t> GenerateNeighbour( std::vector<Double_t>& parameters, Double_t currentTemperature );

      IFitterTarget&                fFitterTarget;           // the fitter target
      TRandom*                      fRandom;                 // random generator
      const std::vector<TMVA::Interval*>& fRanges;                 // parameter ranges

      // fitter setup
      Int_t                         fMaxCalls;               // maximum number of minimisation calls
      Double_t                      fInitialTemperature;     // initial temperature
      Double_t                      fMinTemperature;         // minimum temperature
      Double_t                      fEps;                    // epsilon
      Double_t                      fTemperatureScale;       // how fast temperature change
      Double_t                      fAdaptiveSpeed;          // how fast temperature change in adaptive (in adaptive two variables describe
                                                             // the change of temperature, but fAdaptiveSpeed should be 1.0 and its not
                                                             // recommended to change it)
      Double_t                      fTemperatureAdaptiveStep;// used to calculate InitialTemperature if fUseDefaultTemperature

      Bool_t                        fUseDefaultScale;        // if TRUE, SA calculates its own TemperatureScale
      Bool_t                        fUseDefaultTemperature;  // if TRUE, SA calculates its own InitialTemperature (MinTemperautre)

      mutable MsgLogger*            fLogger;   // message logger
      MsgLogger& Log() const { return *fLogger; }

      Double_t fProgress;

      // variables for JsMVA
      UInt_t *fIPyCurrentIter = nullptr;
      bool * fExitFromTraining = nullptr;

      ClassDef(SimulatedAnnealing,0);  // Base class for Simulated Annealing fitting
   };

} // namespace TMVA

#endif

