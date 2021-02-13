// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Joerg Stelzer, Fredrik Tegenfeldt, Helge Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : RuleFit                                                               *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      A class implementing various fits of rule ensembles                       *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Fredrik Tegenfeldt <Fredrik.Tegenfeldt@cern.ch> - Iowa State U., USA      *
 *      Helge Voss         <Helge.Voss@cern.ch>         - MPI-KP Heidelberg, Ger. *
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         *
 *      Iowa State U.                                                             *
 *      MPI-K Heidelberg, Germany                                                 *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

#ifndef ROOT_TMVA_RuleFit
#define ROOT_TMVA_RuleFit

#include "TMVA/DecisionTree.h"
#include "TMVA/RuleEnsemble.h"
#include "TMVA/RuleFitParams.h"
#include "TMVA/Event.h"

#include <algorithm>
#include <random>
#include <vector>

namespace TMVA {


   class MethodBase;
   class MethodRuleFit;
   class MsgLogger;

   class RuleFit {

   public:

      // main constructor
      RuleFit( const TMVA::MethodBase *rfbase );

      // empty constructor
      RuleFit( void );

      virtual ~RuleFit( void );

      void InitNEveEff();
      void InitPtrs( const TMVA::MethodBase *rfbase );
      void Initialize(  const TMVA::MethodBase *rfbase );

      void SetMsgType( EMsgType t );

      void SetTrainingEvents( const std::vector<const TMVA::Event *> & el );

      void ReshuffleEvents()
      {
         std::shuffle(fTrainingEventsRndm.begin(), fTrainingEventsRndm.end(), fRNGEngine);
      }

      void SetMethodBase( const MethodBase *rfbase );

      // make the forest of trees for rule generation
      void MakeForest();

      // build a tree
      void BuildTree( TMVA::DecisionTree *dt );

      // save event weights
      void SaveEventWeights();

      // restore saved event weights
      void RestoreEventWeights();

      // boost events based on the given tree
      void Boost( TMVA::DecisionTree *dt );

      // calculate and print some statistics on the given forest
      void ForestStatistics();

      // calculate the discriminating variable for the given event
      Double_t EvalEvent( const Event& e );

      // calculate sum of 
      Double_t CalcWeightSum( const std::vector<const TMVA::Event *> *events, UInt_t neve=0 );

      // do the fitting of the coefficients
      void     FitCoefficients();

      // calculate variable and rule importance from a set of events
      void     CalcImportance();

      // set usage of linear term
      void     SetModelLinear()                      { fRuleEnsemble.SetModelLinear(); }
      // set usage of rules
      void     SetModelRules()                       { fRuleEnsemble.SetModelRules(); }
      // set usage of linear term
      void     SetModelFull()                        { fRuleEnsemble.SetModelFull(); }
      // set minimum importance allowed
      void     SetImportanceCut( Double_t minimp=0 ) { fRuleEnsemble.SetImportanceCut(minimp); }
      // set minimum rule distance - see RuleEnsemble
      void     SetRuleMinDist( Double_t d )          { fRuleEnsemble.SetRuleMinDist(d); }
      // set path related parameters
      void     SetGDTau( Double_t t=0.0 )       { fRuleFitParams.SetGDTau(t); }
      void     SetGDPathStep( Double_t s=0.01 ) { fRuleFitParams.SetGDPathStep(s); }
      void     SetGDNPathSteps( Int_t n=100 )   { fRuleFitParams.SetGDNPathSteps(n); }
      // make visualization histograms
      void     SetVisHistsUseImp( Bool_t f ) { fVisHistsUseImp = f; }
      void     UseImportanceVisHists()       { fVisHistsUseImp = kTRUE; }
      void     UseCoefficientsVisHists()     { fVisHistsUseImp = kFALSE; }
      void     MakeVisHists();
      void     FillVisHistCut(const Rule * rule, std::vector<TH2F *> & hlist);
      void     FillVisHistCorr(const Rule * rule, std::vector<TH2F *> & hlist);
      void     FillCut(TH2F* h2,const TMVA::Rule *rule,Int_t vind);
      void     FillLin(TH2F* h2,Int_t vind);
      void     FillCorr(TH2F* h2,const TMVA::Rule *rule,Int_t v1, Int_t v2);
      void     NormVisHists(std::vector<TH2F *> & hlist);
      void     MakeDebugHists();
      Bool_t   GetCorrVars(TString & title, TString & var1, TString & var2);
      // accessors
      UInt_t        GetNTreeSample()            const { return fNTreeSample; }
      Double_t      GetNEveEff()                const { return fNEveEffTrain; } // reweighted number of events = sum(wi)
      const Event*  GetTrainingEvent(UInt_t i)  const { return static_cast< const Event *>(fTrainingEvents[i]); }
      Double_t      GetTrainingEventWeight(UInt_t i)  const { return fTrainingEvents[i]->GetWeight(); }

      //      const Event*  GetTrainingEvent(UInt_t i, UInt_t isub)  const { return &(fTrainingEvents[fSubsampleEvents[isub]])[i]; }

      const std::vector< const TMVA::Event * > & GetTrainingEvents()  const { return fTrainingEvents; }
      //      const std::vector< Int_t >               & GetSubsampleEvents() const { return fSubsampleEvents; }

      //      void  GetSubsampleEvents(Int_t sub, UInt_t & ibeg, UInt_t & iend) const;
      void  GetRndmSampleEvents(std::vector< const TMVA::Event * > & evevec, UInt_t nevents);
      //
      const std::vector< const TMVA::DecisionTree *> & GetForest()     const { return fForest; }
      const RuleEnsemble                       & GetRuleEnsemble()     const { return fRuleEnsemble; }
      RuleEnsemble                       * GetRuleEnsemblePtr()        { return &fRuleEnsemble; }
      const RuleFitParams                      & GetRuleFitParams()    const { return fRuleFitParams; }
      RuleFitParams                      * GetRuleFitParamsPtr()       { return &fRuleFitParams; }
      const MethodRuleFit                      * GetMethodRuleFit()    const { return fMethodRuleFit; }
      const MethodBase                         * GetMethodBase()       const { return fMethodBase; }

   private:

      // copy constructor
      RuleFit( const RuleFit & other );

      // copy method
      void Copy( const RuleFit & other );

      std::vector<const TMVA::Event *>    fTrainingEvents;      // all training events
      std::vector<const TMVA::Event *>    fTrainingEventsRndm;  // idem, but randomly shuffled
      std::vector<Double_t>               fEventWeights;        // original weights of the events - follows fTrainingEvents
      UInt_t                              fNTreeSample;         // number of events in sub sample = frac*neve

      Double_t                            fNEveEffTrain;    // reweighted number of events = sum(wi)
      std::vector< const TMVA::DecisionTree *>  fForest;    // the input forest of decision trees
      RuleEnsemble                        fRuleEnsemble;    // the ensemble of rules
      RuleFitParams                       fRuleFitParams;   // fit rule parameters
      const MethodRuleFit                *fMethodRuleFit;   // pointer the method which initialized this RuleFit instance
      const MethodBase                   *fMethodBase;      // pointer the method base which initialized this RuleFit instance
      Bool_t                              fVisHistsUseImp;  // if true, use importance as weight; else coef in vis hists

      mutable MsgLogger*                  fLogger;   // message logger
      MsgLogger& Log() const { return *fLogger; }    

      static const Int_t randSEED = 0; // set to 1 for debugging purposes or to zero for random seeds
      std::default_random_engine fRNGEngine;

      ClassDef(RuleFit,0);  // Calculations for Friedman's RuleFit method
   };
}

#endif
