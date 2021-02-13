// @(#)root/tmva $Id$
// Author: Peter Speckmayer

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : TMVA::GeneticAlgorithm                                                *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Implementation (see header for description)                               *
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

/*! \class TMVA::GeneticAlgorithm
\ingroup TMVA

Base definition for genetic algorithm

*/

#include <algorithm>
#include <cfloat>

#ifdef _GLIBCXX_PARALLEL
#include <omp.h>
#endif

#include "TMVA/GeneticAlgorithm.h"
#include "TMVA/Interval.h"
#include "TMVA/IFitterTarget.h"
#include "TMVA/MsgLogger.h"
#include "TMVA/Types.h"

#include "RtypesCore.h"
#include "Rtypes.h"
#include "TMath.h"

namespace TMVA {
   const Bool_t GeneticAlgorithm__DEBUG__ = kFALSE;
}

ClassImp(TMVA::GeneticAlgorithm);

////////////////////////////////////////////////////////////////////////////////
/// Constructor
///
/// Parameters:
///
///  - int populationSize : defines the number of "Individuals" which are created and tested
///                          within one Generation (Iteration of the Evolution)
///  - std::vector<TMVA::Interval*> ranges : Interval holds the information of an interval, where the GetMin
///                          gets the low and GetMax gets the high constraint of the variable
///                          the size of "ranges" is the number of coefficients which are optimised
/// Purpose:
///
///     Creates a random population with individuals of the size ranges.size()

TMVA::GeneticAlgorithm::GeneticAlgorithm( IFitterTarget& target, Int_t populationSize,
                                          const std::vector<Interval*>& ranges, UInt_t seed )
: fConvCounter(-1),
   fFitterTarget( target ),
   fConvValue(0.),
   fLastResult(DBL_MAX),
   fSpread(0.1),
   fMirror(kTRUE),
   fFirstTime(kTRUE),
   fMakeCopies(kFALSE),
   fPopulationSize(populationSize),
   fRanges( ranges ),
   fPopulation(ranges, populationSize, seed),
   fBestFitness(DBL_MAX),
   fLogger( new MsgLogger("GeneticAlgorithm") )
{
   fPopulation.SetRandomSeed( seed );
}

TMVA::GeneticAlgorithm::~GeneticAlgorithm()
{
   // destructor; deletes fLogger
   delete fLogger;
}


////////////////////////////////////////////////////////////////////////////////
/// calls evolution, but if it is not the first time.
/// If it's the first time, the random population created by the
/// constructor is still not evaluated, .. therefore we wait for the
/// second time init is called.

void TMVA::GeneticAlgorithm::Init()
{
   if ( fFirstTime ) fFirstTime = kFALSE;
   else {
      Evolution();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// if the "fitnessFunction" is called multiple times for one set of
/// factors (because i.e. each event of a TTree has to be assessed with
/// each set of Factors proposed by the Genetic Algorithm) the value
/// of the current calculation has to be added(? or else) to the value
/// obtained up to now.
/// example: some chi-square is calculated for every event,
/// after every event the new chi-square (newValue) has to be simply
/// added to the oldValue.
///
/// this function has to be overridden eventually
/// it might contain only the following return statement.
///        return oldValue + newValue;

Double_t TMVA::GeneticAlgorithm::NewFitness( Double_t /*oldValue*/, Double_t newValue )
{
   return newValue;
}

////////////////////////////////////////////////////////////////////////////////
/// starts the evaluation of the fitness of all different individuals of
/// the population.
///
/// this function calls implicitly (many times) the "fitnessFunction" which
/// has been overridden by the user.

Double_t TMVA::GeneticAlgorithm::CalculateFitness()
{
   fBestFitness = DBL_MAX;
#ifdef _GLIBCXX_PARALLEL

   const int nt = omp_get_num_threads();
   Double_t bests[nt];
   for ( int i =0; i < nt; ++i )
      bests[i] = fBestFitness;

#pragma omp parallel
   {
      int thread_number = omp_get_thread_num();
#pragma omp for
      for ( int index = 0; index < fPopulation.GetPopulationSize(); ++index )
         {
            GeneticGenes* genes = fPopulation.GetGenes(index);
            Double_t fitness = NewFitness( genes->GetFitness(),
                                           fFitterTarget.EstimatorFunction(genes->GetFactors()) );
            genes->SetFitness( fitness );

            if ( bests[thread_number] > fitness )
               bests[thread_number] = fitness;
         }
   }

   fBestFitness = *std::min_element(bests, bests+nt);

#else

   for ( int index = 0; index < fPopulation.GetPopulationSize(); ++index ) {
      GeneticGenes* genes = fPopulation.GetGenes(index);
      Double_t fitness = NewFitness( genes->GetFitness(),
                                     fFitterTarget.EstimatorFunction(genes->GetFactors()) );
      genes->SetFitness( fitness );

      if ( fBestFitness  > fitness )
         fBestFitness = fitness;

   }

#endif

   fPopulation.Sort();

   return fBestFitness;
}

////////////////////////////////////////////////////////////////////////////////
/// this function is called from "init" and controls the evolution of the
/// individuals.
///
/// The function can be overridden to change the parameters for mutation rate
/// sexual reproduction and so on.

void TMVA::GeneticAlgorithm::Evolution()
{
   if ( fMakeCopies )
      fPopulation.MakeCopies( 5 );
   fPopulation.MakeChildren();

   fPopulation.Mutate( 10, 3, kTRUE, fSpread, fMirror );
   fPopulation.Mutate( 40, fPopulation.GetPopulationSize()*3/4 );
}

////////////////////////////////////////////////////////////////////////////////
/// this function provides the ability to change the stepSize of a mutation according to
/// the success of the last generations.
///
/// Parameters:
///
///  - int ofSteps :  = if OF the number of STEPS given in this variable (ofSteps)
///  - int successSteps : >sucessSteps Generations could improve the result
///  - double factor : than multiply the stepSize ( spread ) by this factor
///
/// (if ofSteps == successSteps nothing is changed, if ofSteps < successSteps, the spread
/// is divided by the factor)
///
/// using this function one can increase the stepSize of the mutation when we have
/// good success (to pass fast through the easy phase-space) and reduce the stepSize
/// if we are in a difficult "territory" of the phase-space.

Double_t TMVA::GeneticAlgorithm::SpreadControl( Int_t ofSteps, Int_t successSteps, Double_t factor )
{
   // < is valid for "less" comparison
   if ( fBestFitness < fLastResult || fSuccessList.size() <=0 ) {
      fLastResult = fBestFitness;
      fSuccessList.push_front( 1 ); // it got better
   }
   else {
      fSuccessList.push_front( 0 ); // it stayed the same
   }
   Int_t n = 0;
   Int_t sum = 0;
   std::deque<Int_t>::iterator vec = fSuccessList.begin();
   for (; vec != fSuccessList.end() ; ++vec) {
      sum += *vec;
      n++;
   }

   if ( n >= ofSteps ) {
      fSuccessList.pop_back();
      if ( sum > successSteps ) { // too much success
         fSpread /= factor;
         if (GeneticAlgorithm__DEBUG__) Log() << kINFO << ">" << std::flush;
      }
      else if ( sum == successSteps ) { // on the optimal path
         if (GeneticAlgorithm__DEBUG__) Log() << "=" << std::flush;
      }
      else {        // not very successful
         fSpread *= factor;
         if (GeneticAlgorithm__DEBUG__) Log() << "<" << std::flush;
      }
   }

   return fSpread;
}

////////////////////////////////////////////////////////////////////////////////
/// gives back true if the last "steps" steps have lead to an improvement of the
/// "fitness" of the "individuals" of at least "improvement"
///
/// this gives a simple measure of if the fitness of the individuals is
/// converging and no major improvement is to be expected soon.

Bool_t TMVA::GeneticAlgorithm::HasConverged( Int_t steps, Double_t improvement )
{
   if (fConvCounter < 0) {
      fConvValue = fBestFitness;
   }
   if (TMath::Abs(fBestFitness - fConvValue) <= improvement || steps<0) {
      fConvCounter ++;
   }
   else {
      fConvCounter = 0;
      fConvValue = fBestFitness;
   }
   if (GeneticAlgorithm__DEBUG__) Log() << "." << std::flush;
   if (fConvCounter < steps) return kFALSE;
   return kTRUE;
}
