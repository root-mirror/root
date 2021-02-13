// @(#)root/tmva $Id$
// Author: Peter Speckmayer

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : GeneticFitter                                                         *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Implementation                                                            *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Peter Speckmayer <speckmay@mail.cern.ch> - CERN, Switzerland              *
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         *
 *      MPI-K Heidelberg, Germany                                                 *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

/*! \class TMVA::GeneticFitter
\ingroup TMVA

Fitter using a Genetic Algorithm.

*/

#include "TMVA/GeneticFitter.h"

#include "TMVA/Configurable.h"
#include "TMVA/GeneticAlgorithm.h"
#include "TMVA/Interval.h"
#include "TMVA/FitterBase.h"
#include "TMVA/MsgLogger.h"
#include "TMVA/Timer.h"
#include "TMVA/Types.h"

#include "Rtypes.h"
#include "TString.h"

ClassImp(TMVA::GeneticFitter);

////////////////////////////////////////////////////////////////////////////////
/// constructor

TMVA::GeneticFitter::GeneticFitter( IFitterTarget& target,
                                    const TString& name,
                                    const std::vector<TMVA::Interval*>& ranges,
                                    const TString& theOption )
: FitterBase( target, name, ranges, theOption )
{
   // default parameters settings for Genetic Algorithm
   DeclareOptions();
   ParseOptions();
}

////////////////////////////////////////////////////////////////////////////////
/// declare GA options

void TMVA::GeneticFitter::DeclareOptions()
{
   DeclareOptionRef( fPopSize=300,    "PopSize",   "Population size for GA" );
   DeclareOptionRef( fNsteps=40,      "Steps",     "Number of steps for convergence" );
   DeclareOptionRef( fCycles=3,       "Cycles",    "Independent cycles of GA fitting" );
   DeclareOptionRef( fSC_steps=10,    "SC_steps",  "Spread control, steps" );
   DeclareOptionRef( fSC_rate=5,      "SC_rate",   "Spread control, rate: factor is changed depending on the rate" );
   DeclareOptionRef( fSC_factor=0.95, "SC_factor", "Spread control, factor" );
   DeclareOptionRef( fConvCrit=0.001, "ConvCrit",   "Convergence criteria" );

   DeclareOptionRef( fSaveBestFromGeneration=1, "SaveBestGen",
                     "Saves the best n results from each generation. They are included in the last cycle" );
   DeclareOptionRef( fSaveBestFromCycle=10,     "SaveBestCycle",
                     "Saves the best n results from each cycle. They are included in the last cycle. The value should be set to at least 1.0" );

   DeclareOptionRef( fTrim=kFALSE, "Trim",
                     "Trim the population to PopSize after assessing the fitness of each individual" );
   DeclareOptionRef( fSeed=100, "Seed", "Set seed of random generator (0 gives random seeds)" );
}

////////////////////////////////////////////////////////////////////////////////
/// set GA configuration parameters

void TMVA::GeneticFitter::SetParameters(  Int_t cycles,
                                          Int_t nsteps,
                                          Int_t popSize,
                                          Int_t SC_steps,
                                          Int_t SC_rate,
                                          Double_t SC_factor,
                                          Double_t convCrit)
{
   fNsteps    = nsteps;
   fCycles    = cycles;
   fPopSize   = popSize;
   fSC_steps  = SC_steps;
   fSC_rate   = SC_rate;
   fSC_factor = SC_factor;
   fConvCrit  = convCrit;
}

////////////////////////////////////////////////////////////////////////////////
/// Execute fitting

Double_t TMVA::GeneticFitter::Run( std::vector<Double_t>& pars )
{
   Log() << kHEADER << "<GeneticFitter> Optimisation, please be patient "
         << "... (inaccurate progress timing for GA)" << Endl;

   GetFitterTarget().ProgressNotifier( "GA", "init" );

   GeneticAlgorithm gstore( GetFitterTarget(),  fPopSize, fRanges);
   //   gstore.SetMakeCopies(kTRUE);  // commented out, because it reduces speed

   // timing of GA
   Timer timer( 100*(fCycles), GetName() );
   if (fIPyMaxIter) *fIPyMaxIter = 100*(fCycles);
   timer.DrawProgressBar( 0 );

   Double_t progress = 0.;

   for (Int_t cycle = 0; cycle < fCycles; cycle++) {
     if (fIPyCurrentIter) *fIPyCurrentIter = 100*(cycle);
     if (fExitFromTraining && *fExitFromTraining) break;
      GetFitterTarget().ProgressNotifier( "GA", "cycle" );
      // ---- perform series of fits to achieve best convergence

      // "m_ga_spread" times the number of variables
      GeneticAlgorithm ga( GetFitterTarget(), fPopSize, fRanges, fSeed );
      //      ga.SetMakeCopies(kTRUE);  // commented out, because it reduces speed

      if ( pars.size() == fRanges.size() ){
         ga.GetGeneticPopulation().GiveHint( pars, 0.0 );
      }
      if (cycle==fCycles-1) {
         GetFitterTarget().ProgressNotifier( "GA", "last" );
         ga.GetGeneticPopulation().AddPopulation( gstore.GetGeneticPopulation() );
      }

      GetFitterTarget().ProgressNotifier( "GA", "iteration" );

      ga.CalculateFitness();
      ga.GetGeneticPopulation().TrimPopulation();

      Double_t n=0.;
      do {
         GetFitterTarget().ProgressNotifier( "GA", "iteration" );
         ga.Init();
         ga.CalculateFitness();
         if ( fTrim ) ga.GetGeneticPopulation().TrimPopulation();
         ga.SpreadControl( fSC_steps, fSC_rate, fSC_factor );

         // monitor progrss
         if (ga.fConvCounter > n) n = Double_t(ga.fConvCounter);
         progress = 100*((Double_t)cycle) + 100*(n/Double_t(fNsteps));

         timer.DrawProgressBar( (Int_t)progress );

         // Copy the best genes of the generation
         ga.GetGeneticPopulation().Sort();
         for ( Int_t i = 0; i<fSaveBestFromGeneration && i<fPopSize; i++ ) {
            gstore.GetGeneticPopulation().GiveHint( ga.GetGeneticPopulation().GetGenes(i)->GetFactors(),
                                                    ga.GetGeneticPopulation().GetGenes(i)->GetFitness() );
         }
      } while (!ga.HasConverged( fNsteps, fConvCrit ));

      timer.DrawProgressBar( 100*(cycle+1) );

      ga.GetGeneticPopulation().Sort();
      for ( Int_t i = 0; i<fSaveBestFromGeneration && i<fPopSize; i++ ) {
         gstore.GetGeneticPopulation().GiveHint( ga.GetGeneticPopulation().GetGenes(i)->GetFactors(),
                                                 ga.GetGeneticPopulation().GetGenes(i)->GetFitness() );
      }
   }

   // get elapsed time
   Log() << kINFO << "Elapsed time: " << timer.GetElapsedTime()
         << "                            " << Endl;

   Double_t fitness = gstore.CalculateFitness();
   gstore.GetGeneticPopulation().Sort();
   pars.swap( gstore.GetGeneticPopulation().GetGenes(0)->GetFactors() );

   GetFitterTarget().ProgressNotifier( "GA", "stop" );
   return fitness;
}
