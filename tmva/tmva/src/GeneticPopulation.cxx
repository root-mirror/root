// @(#)root/tmva $Id$
// Author: Peter Speckmayer

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : TMVA::GeneticPopulation                                               *
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

/*! \class TMVA::GeneticPopulation
\ingroup TMVA

Population definition for genetic algorithm.

*/

#include <algorithm>

#include "Rstrstream.h"
#include "TRandom3.h"
#include "TH1.h"

#include "TMVA/GeneticPopulation.h"
#include "TMVA/GeneticGenes.h"
#include "TMVA/MsgLogger.h"

ClassImp(TMVA::GeneticPopulation);

using namespace std;

////////////////////////////////////////////////////////////////////////////////
/// Constructor

TMVA::GeneticPopulation::GeneticPopulation(const std::vector<Interval*>& ranges, Int_t size, UInt_t seed)
   : fGenePool(size),
     fRanges(ranges.size()),
     fLogger( new MsgLogger("GeneticPopulation") )
{
   // create a randomGenerator for this population and set a seed
   // create the genePools
   //
   fRandomGenerator = new TRandom3( 100 ); //please check
   fRandomGenerator->Uniform(0.,1.);
   fRandomGenerator->SetSeed( seed );

   for ( unsigned int i = 0; i < ranges.size(); ++i )
      fRanges[i] = new TMVA::GeneticRange( fRandomGenerator, ranges[i] );

   vector<Double_t> newEntry( fRanges.size() );
   for ( int i = 0; i < size; ++i )
      {
         for ( unsigned int rIt = 0; rIt < fRanges.size(); ++rIt )
            newEntry[rIt] = fRanges[rIt]->Random();
         fGenePool[i] = TMVA::GeneticGenes( newEntry);
      }

   fPopulationSizeLimit = size;
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TMVA::GeneticPopulation::~GeneticPopulation()
{
   if (fRandomGenerator != NULL) delete fRandomGenerator;

   std::vector<GeneticRange*>::iterator it = fRanges.begin();
   for (;it!=fRanges.end(); ++it) delete *it;

   delete fLogger;
}



////////////////////////////////////////////////////////////////////////////////
/// the random seed of the random generator

void TMVA::GeneticPopulation::SetRandomSeed( UInt_t seed )
{
   fRandomGenerator->SetSeed( seed );
}

////////////////////////////////////////////////////////////////////////////////
/// Produces offspring which is are copies of their parents.
///
/// Parameters:
///  - int number : the number of the last individual to be copied

void TMVA::GeneticPopulation::MakeCopies( int number )
{
   int i=0;
   for (std::vector<TMVA::GeneticGenes>::iterator it = fGenePool.begin();
        it != fGenePool.end() && i < number;
        ++it, ++i ) {
      GiveHint( it->GetFactors(), it->GetFitness() );
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Creates children out of members of the current generation.
///
/// Children have a combination of the coefficients of their parents

void TMVA::GeneticPopulation::MakeChildren()
{
#ifdef _GLIBCXX_PARALLEL
#pragma omp parallel
#pragma omp for
#endif
   for ( int it = 0; it < (int) (fGenePool.size() / 2); ++it )
      {
         Int_t pos = (Int_t)fRandomGenerator->Integer( fGenePool.size()/2 );
         fGenePool[(fGenePool.size() / 2) + it] = MakeSex( fGenePool[it], fGenePool[pos] );
      }
}

////////////////////////////////////////////////////////////////////////////////
/// this function takes two individuals and produces offspring by mixing
/// (recombining) their coefficients.

TMVA::GeneticGenes TMVA::GeneticPopulation::MakeSex( TMVA::GeneticGenes male,
                                                     TMVA::GeneticGenes female )
{
   vector< Double_t > child(fRanges.size());
   for (unsigned int i = 0; i < fRanges.size(); ++i) {
      if (fRandomGenerator->Integer( 2 ) == 0) {
         child[i] = male.GetFactors()[i];
      }else{
         child[i] = female.GetFactors()[i];
      }
   }
   return TMVA::GeneticGenes( child );
}

////////////////////////////////////////////////////////////////////////////////
/// Mutates the individuals in the genePool.
///
/// Parameters:
///
///  - double probability : gives the probability (in percent) of a mutation of a coefficient
///  - int startIndex : leaves unchanged (without mutation) the individuals which are better ranked
///                     than indicated by "startIndex". This means: if "startIndex==3", the first (and best)
///                     three individuals are not mutated. This allows to preserve the best result of the
///                     current Generation for the next generation.
///  - Bool_t near : if true, the mutation will produce a new coefficient which is "near" the old one
///                     (gaussian around the current value)
///  - double spread : if near==true, spread gives the sigma of the gaussian
///  - Bool_t mirror : if the new value obtained would be outside of the given constraints
///                    the value is mapped between the constraints again. This can be done either
///                    by a kind of periodic boundary conditions or mirrored at the boundary.
///                    (mirror = true seems more "natural")

void TMVA::GeneticPopulation::Mutate( Double_t probability , Int_t startIndex,
                                      Bool_t near, Double_t spread, Bool_t mirror )
{
   vector< Double_t>::iterator vec;
   vector< TMVA::GeneticRange* >::iterator vecRange;

   //#ifdef _GLIBCXX_PARALLEL
   // #pragma omp parallel
   // #pragma omp for
   //#endif
   // The range methods are not thread safe!
   for (int it = startIndex; it < (int) fGenePool.size(); ++it) {
      vecRange = fRanges.begin();
      for (vec = (fGenePool[it].GetFactors()).begin(); vec < (fGenePool[it].GetFactors()).end(); ++vec) {
         if (fRandomGenerator->Uniform( 100 ) <= probability) {
            (*vec) = (*vecRange)->Random( near, (*vec), spread, mirror );
         }
         ++vecRange;
      }
   }
}


////////////////////////////////////////////////////////////////////////////////
/// gives back the "Genes" of the population with the given index.

TMVA::GeneticGenes* TMVA::GeneticPopulation::GetGenes( Int_t index )
{
   return &(fGenePool[index]);
}

////////////////////////////////////////////////////////////////////////////////
/// make a little printout of the individuals up to index "untilIndex"
/// this means, .. write out the best "untilIndex" individuals.

void TMVA::GeneticPopulation::Print( Int_t untilIndex )
{
   for ( unsigned int it = 0; it < fGenePool.size(); ++it )
      {
         Int_t n=0;
         if (untilIndex >= -1 ) {
            if (untilIndex == -1 ) return;
            untilIndex--;
         }
         Log() << "fitness: " << fGenePool[it].GetFitness() << "    ";
         for (vector< Double_t >::iterator vec = fGenePool[it].GetFactors().begin();
              vec < fGenePool[it].GetFactors().end(); ++vec ) {
            Log() << "f_" << n++ << ": " << (*vec) << "     ";
         }
         Log() << Endl;
      }
}

////////////////////////////////////////////////////////////////////////////////
/// make a little printout to the stream "out" of the individuals up to index "untilIndex"
/// this means, .. write out the best "untilIndex" individuals.

void TMVA::GeneticPopulation::Print( ostream & out, Int_t untilIndex )
{
   for ( unsigned int it = 0; it < fGenePool.size(); ++it ) {
      Int_t n=0;
      if (untilIndex >= -1 ) {
         if (untilIndex == -1 ) return;
         untilIndex--;
      }
      out << "fitness: " << fGenePool[it].GetFitness() << "    ";
      for (vector< Double_t >::iterator vec = fGenePool[it].GetFactors().begin();
           vec < fGenePool[it].GetFactors().end(); ++vec ) {
         out << "f_" << n++ << ": " << (*vec) << "     ";
      }
      out << std::endl;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// give back a histogram with the distribution of the coefficients.
///
/// Parameters:
///
///  - int bins : number of bins of the histogram
///  - int min : histogram minimum
///  - int max : maximum value of the histogram

TH1F* TMVA::GeneticPopulation::VariableDistribution( Int_t varNumber, Int_t bins,
                                                     Int_t min, Int_t max )
{
   std::cout << "FAILED! TMVA::GeneticPopulation::VariableDistribution" << std::endl;

   std::stringstream histName;
   histName.clear();
   histName.str("v");
   histName << varNumber;
   TH1F *hist = new TH1F( histName.str().c_str(),histName.str().c_str(), bins,min,max );

   return hist;
}

////////////////////////////////////////////////////////////////////////////////
/// gives back all the values of coefficient "varNumber" of the current generation

vector<Double_t> TMVA::GeneticPopulation::VariableDistribution( Int_t /*varNumber*/ )
{
   std::cout << "FAILED! TMVA::GeneticPopulation::VariableDistribution" << std::endl;

   vector< Double_t > varDist;

   return varDist;
}

////////////////////////////////////////////////////////////////////////////////
/// add another population (strangers) to the one of this GeneticPopulation

void TMVA::GeneticPopulation::AddPopulation( GeneticPopulation *strangers )
{
   for (std::vector<TMVA::GeneticGenes>::iterator it = strangers->fGenePool.begin();
        it != strangers->fGenePool.end(); ++it ) {
      GiveHint( it->GetFactors(), it->GetFitness() );
   }
}

////////////////////////////////////////////////////////////////////////////////
/// add another population (strangers) to the one of this GeneticPopulation

void TMVA::GeneticPopulation::AddPopulation( GeneticPopulation &strangers )
{
   AddPopulation(&strangers);
}

////////////////////////////////////////////////////////////////////////////////
/// trim the population to the predefined size

void TMVA::GeneticPopulation::TrimPopulation()
{
   std::sort(fGenePool.begin(), fGenePool.end());
   while ( fGenePool.size() > (unsigned int) fPopulationSizeLimit )
      fGenePool.pop_back();
}

////////////////////////////////////////////////////////////////////////////////
/// add an individual (a set of variables) to the population
/// if there is a set of variables which is known to perform good, they can be given as a hint to the population

void TMVA::GeneticPopulation::GiveHint( std::vector< Double_t >& hint, Double_t fitness )
{
   TMVA::GeneticGenes g(hint);
   g.SetFitness(fitness);

   fGenePool.push_back( g );
}

////////////////////////////////////////////////////////////////////////////////
/// sort the genepool according to the fitness of the individuals

void TMVA::GeneticPopulation::Sort()
{
   std::sort(fGenePool.begin(), fGenePool.end());
}

