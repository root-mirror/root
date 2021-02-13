// @(#)root/tmva $Id$    
// Author: Peter Speckmayer

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : GeneticPopulation                                                     *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *    Population definition for genetic algorithm                                 *
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

#ifndef ROOT_TMVA_GeneticPopulation
#define ROOT_TMVA_GeneticPopulation

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// GeneticPopulation                                                    //
//                                                                      //
// Population definition for genetic algorithm                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <vector>

#include "TMVA/GeneticGenes.h"
#include "TMVA/Interval.h"
#include "TMVA/GeneticRange.h"

class TH1F;

namespace TMVA {

   class MsgLogger;

   class GeneticPopulation {

   public:

      GeneticPopulation(const std::vector<TMVA::Interval*>& ranges, Int_t size, UInt_t seed = 0);
      virtual ~GeneticPopulation();

      void SetRandomSeed( UInt_t seed = 0);

      void MakeChildren();
      void Mutate( Double_t probability = 20, Int_t startIndex = 0, Bool_t near = kFALSE, 
                   Double_t spread = 0.1, Bool_t mirror = kFALSE  );

      GeneticGenes* GetGenes( Int_t index );
      Int_t         GetPopulationSize() const { return fGenePool.size(); }
      Double_t      GetFitness() const { return fGenePool.size()>0? fGenePool[0].GetFitness() : 0; }

      const std::vector<TMVA::GeneticGenes>& GetGenePool() const { return fGenePool; }
      const std::vector<TMVA::GeneticRange*>& GetRanges() const { return fRanges; }

      std::vector<TMVA::GeneticGenes>&  GetGenePool() { return fGenePool; }
      std::vector<TMVA::GeneticRange*>& GetRanges()   { return fRanges; }

      void Print( Int_t untilIndex = -1 );
      void Print( std::ostream & out, Int_t utilIndex = -1 );

      TH1F* VariableDistribution( Int_t varNumber, Int_t bins, Int_t min, Int_t max  );
      std::vector< Double_t > VariableDistribution( Int_t varNumber );

      // To keep compatibility: These methods might be reimplemented
      // or just eliminated later on. They are used by the
      // GeneticFitter class.
     
      void MakeCopies( int number );
      void NextGeneration() {}
      void AddPopulation( GeneticPopulation *strangers );
      void AddPopulation( GeneticPopulation &strangers );
      void TrimPopulation();
      void GiveHint( std::vector< Double_t >& hint, Double_t fitness = 0 );
      void Sort();

   private:
      GeneticGenes MakeSex( GeneticGenes male, GeneticGenes female );
  
   private:

      std::vector<TMVA::GeneticGenes>  fGenePool;    // the "genePool" where the individuals of the current generation are stored
      std::vector<TMVA::GeneticRange*> fRanges;      // contains the ranges inbetween the values of the coefficients have to be

      TRandom3*fRandomGenerator;    // random Generator for this population

      mutable MsgLogger* fLogger;   // message logger
      MsgLogger& Log() const { return *fLogger; }    

      Int_t fPopulationSizeLimit;

      ClassDef(GeneticPopulation,0); //Population definition for genetic algorithm
   };

} // namespace TMVA

#endif
