// @(#)root/tmva $Id$
// Author: Peter Speckmayer

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : GeneticFitter                                                         *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *       Fitter using a Genetic Algorithm                                         *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Peter Speckmayer <speckmay@mail.cern.ch>  - CERN, Switzerland             *
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         *
 *      MPI-K Heidelberg, Germany                                                 *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

#ifndef ROOT_TMVA_GeneticFitter
#define ROOT_TMVA_GeneticFitter

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// GeneticFitter                                                        //
//                                                                      //
// Fitter using a Genetic Algorithm                                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TMVA/FitterBase.h"
#include <vector>

namespace TMVA {

   class IFitterTarget;
   class Interval;

   class GeneticFitter : public FitterBase {

   public:

      GeneticFitter( IFitterTarget& target, const TString& name,
                     const std::vector<TMVA::Interval*>& ranges, const TString& theOption );

      virtual ~GeneticFitter() {}

      void SetParameters( Int_t cycles,
                          Int_t nsteps,
                          Int_t popSize,
                          Int_t SC_steps,
                          Int_t SC_rate,
                          Double_t SC_factor,
                          Double_t convCrit );

      Double_t Run( std::vector<Double_t>& pars );

      Double_t NewFitness( Double_t oldF, Double_t newF ) { return oldF + newF; }

   private:

      void DeclareOptions();

      Int_t fCycles;                    // number of (nearly) independent calculation cycles
      Int_t fNsteps;                    // convergence criteria: if no improvements > fConvCrit was achieved within the last fNsteps: cycle has "converged"
      Int_t fPopSize;                   // number of individuals to start with
      Int_t fSC_steps;                  // regulates how strong the mutations for the coordinates are: if within fSC_steps there were more than...
      Int_t fSC_rate;                   // ... fSC_rate improvements, than multiply the sigma of the gaussian which defines how the random numbers are generated ...
      Double_t fSC_factor;              // ... with fSC_factor; if there were less improvements: divide by that factor; if there were exactly fSC_rate improvements, dont change anything
      Double_t fConvCrit;               // improvements bigger than fConvCrit are counted as "improvement"
      Int_t fSaveBestFromGeneration;    // store the best individuals from one generation (these are included as "hints" in the last cycle of GA calculation)
      Int_t fSaveBestFromCycle;         // store the best individuals from one cycle (these are included as "hints" in the last cycle of GA calculation)
      Bool_t fTrim;                     // take care, that the number of individuals is less fPopSize (trimming is done after the fitness of the individuals is assessed)
      UInt_t fSeed;                     // Seed for the random generator (0 takes random seeds)

      ClassDef(GeneticFitter,0); // Fitter using a Genetic Algorithm
   };

} // namespace TMVA

#endif


