// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Joerg Stelzer, Helge Voss, Kai Voss, Jan Therhaag

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : MethodBDT  (Boosted Decision Trees)                                   *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Analysis of Boosted Decision Trees                                        *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
 *      Helge Voss      <Helge.Voss@cern.ch>     - MPI-K Heidelberg, Germany      *
 *      Kai Voss        <Kai.Voss@cern.ch>       - U. of Victoria, Canada         *
 *      Doug Schouten   <dschoute@sfu.ca>        - Simon Fraser U., Canada        *
 *      Jan Therhaag    <jan.therhaag@cern.ch>   - U. of Bonn, Germany            *
 *                                                                                *
 * Copyright (c) 2005-2011:                                                       *
 *      CERN, Switzerland                                                         *
 *      U. of Victoria, Canada                                                    *
 *      MPI-K Heidelberg, Germany                                                 *
 *      U. of Bonn, Germany                                                       *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

#ifndef ROOT_TMVA_MethodBDT
#define ROOT_TMVA_MethodBDT

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// MethodBDT                                                            //
//                                                                      //
// Analysis of Boosted Decision Trees                                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <vector>
#include <memory>
#include <map>

#include "TH2.h"
#include "TTree.h"
#include "TMVA/MethodBase.h"
#include "TMVA/DecisionTree.h"
#include "TMVA/Event.h"
#include "TMVA/LossFunction.h"

// Multithreading only if the compilation flag is turned on
#ifdef R__USE_IMT
#include <ROOT/TThreadExecutor.hxx>
#include "TSystem.h"
#endif

namespace TMVA {

   class SeparationBase;

   class MethodBDT : public MethodBase {

   public:

      // constructor for training and reading
      MethodBDT( const TString& jobName,
                 const TString& methodTitle,
                 DataSetInfo& theData,
                 const TString& theOption = "");

      // constructor for calculating BDT-MVA using previously generatad decision trees
      MethodBDT( DataSetInfo& theData,
                 const TString& theWeightFile);

      virtual ~MethodBDT( void );

      virtual Bool_t HasAnalysisType( Types::EAnalysisType type, UInt_t numberClasses, UInt_t numberTargets );


      // write all Events from the Tree into a vector of Events, that are
      // more easily manipulated
      void InitEventSample();

      // optimize tuning parameters
      virtual std::map<TString,Double_t> OptimizeTuningParameters(TString fomType="ROCIntegral", TString fitType="FitGA");
      virtual void SetTuneParameters(std::map<TString,Double_t> tuneParameters);

      // training method
      void Train( void );

      // revoke training
      void Reset( void );

      using MethodBase::ReadWeightsFromStream;

      // write weights to file
      void AddWeightsXMLTo( void* parent ) const;

      // read weights from file
      void ReadWeightsFromStream( std::istream& istr );
      void ReadWeightsFromXML(void* parent);

      // write method specific histos to target file
      void WriteMonitoringHistosToFile( void ) const;

      // calculate the MVA value
      Double_t GetMvaValue( Double_t* err = 0, Double_t* errUpper = 0);

      // get the actual forest size (might be less than fNTrees, the requested one, if boosting is stopped early
      UInt_t   GetNTrees() const {return fForest.size();}
   private:

      Double_t GetMvaValue( Double_t* err, Double_t* errUpper, UInt_t useNTrees );
      Double_t PrivateGetMvaValue( const TMVA::Event *ev, Double_t* err=0, Double_t* errUpper=0, UInt_t useNTrees=0 );
      void     BoostMonitor(Int_t iTree);

   public:
      const std::vector<Float_t>& GetMulticlassValues();

      // regression response
      const std::vector<Float_t>& GetRegressionValues();

      // apply the boost algorithm to a tree in the collection
      Double_t Boost( std::vector<const TMVA::Event*>&, DecisionTree *dt, UInt_t cls = 0);

      // ranking of input variables
      const Ranking* CreateRanking();

      // the option handling methods
      void DeclareOptions();
      void ProcessOptions();
      void SetMaxDepth(Int_t d){fMaxDepth = d;}
      void SetMinNodeSize(Double_t sizeInPercent);
      void SetMinNodeSize(TString sizeInPercent);

      void SetNTrees(Int_t d){fNTrees = d;}
      void SetAdaBoostBeta(Double_t b){fAdaBoostBeta = b;}
      void SetNodePurityLimit(Double_t l){fNodePurityLimit = l;} 
      void SetShrinkage(Double_t s){fShrinkage = s;}
      void SetUseNvars(Int_t n){fUseNvars = n;}
      void SetBaggedSampleFraction(Double_t f){fBaggedSampleFraction = f;}


      // get the forest
      inline const std::vector<TMVA::DecisionTree*> & GetForest() const;

      // get the forest
      inline const std::vector<const TMVA::Event*> & GetTrainingEvents() const;

      inline const std::vector<double> & GetBoostWeights() const;

      //return the individual relative variable importance
      std::vector<Double_t> GetVariableImportance();
      Double_t GetVariableImportance(UInt_t ivar);

      Double_t TestTreeQuality( DecisionTree *dt );

      // make ROOT-independent C++ class for classifier response (classifier-specific implementation)
      void MakeClassSpecific( std::ostream&, const TString& ) const;

      // header and auxiliary classes
      void MakeClassSpecificHeader( std::ostream&, const TString& ) const;

      void MakeClassInstantiateNode( DecisionTreeNode *n, std::ostream& fout,
                                     const TString& className ) const;

      void GetHelpMessage() const;

   protected:
      void DeclareCompatibilityOptions();

   private:
      // Init used in the various constructors
      void Init( void );

      void PreProcessNegativeEventWeights();

      // boosting algorithm (adaptive boosting)
      Double_t AdaBoost( std::vector<const TMVA::Event*>&, DecisionTree *dt );

      // boosting algorithm (adaptive boosting with cost matrix)
      Double_t AdaCost( std::vector<const TMVA::Event*>&, DecisionTree *dt );

      // boosting as a random re-weighting
      Double_t Bagging( );

      // boosting special for regression
      Double_t RegBoost( std::vector<const TMVA::Event*>&, DecisionTree *dt );

      // adaboost adapted to regression
      Double_t AdaBoostR2( std::vector<const TMVA::Event*>&, DecisionTree *dt );

      // binomial likelihood gradient boost for classification
      // (see Friedman: "Greedy Function Approximation: a Gradient Boosting Machine"
      // Technical report, Dept. of Statistics, Stanford University)
      Double_t GradBoost( std::vector<const TMVA::Event*>&, DecisionTree *dt, UInt_t cls = 0);
      Double_t GradBoostRegression(std::vector<const TMVA::Event*>&, DecisionTree *dt );
      void InitGradBoost( std::vector<const TMVA::Event*>&);
      void UpdateTargets( std::vector<const TMVA::Event*>&, UInt_t cls = 0);
      void UpdateTargetsRegression( std::vector<const TMVA::Event*>&,Bool_t first=kFALSE);
      Double_t GetGradBoostMVA(const TMVA::Event *e, UInt_t nTrees);
      void     GetBaggedSubSample(std::vector<const TMVA::Event*>&);

      std::vector<const TMVA::Event*>       fEventSample;     // the training events
      std::vector<const TMVA::Event*>       fValidationSample;// the Validation events
      std::vector<const TMVA::Event*>       fSubSample;       // subsample for bagged grad boost
      std::vector<const TMVA::Event*>      *fTrainSample;     // pointer to sample actually used in training (fEventSample or fSubSample) for example

      Int_t                           fNTrees;          // number of decision trees requested
      std::vector<DecisionTree*>      fForest;          // the collection of decision trees
      std::vector<double>             fBoostWeights;    // the weights applied in the individual boosts
      Double_t                        fSigToBkgFraction;// Signal to Background fraction assumed during training
      TString                         fBoostType;       // string specifying the boost type
      Double_t                        fAdaBoostBeta;    // beta parameter for AdaBoost algorithm
      TString                         fAdaBoostR2Loss;  // loss type used in AdaBoostR2 (Linear,Quadratic or Exponential)
      //Double_t                        fTransitionPoint; // break-down point for gradient regression
      Double_t                        fShrinkage;       // learning rate for gradient boost;
      Bool_t                          fBaggedBoost;     // turn bagging in combination with boost on/off
      Bool_t                          fBaggedGradBoost; // turn bagging in combination with grad boost on/off
      //Double_t                        fSumOfWeights;    // sum of all event weights
      //std::map< const TMVA::Event*, std::pair<Double_t, Double_t> >       fWeightedResiduals;  // weighted regression residuals
      std::map< const TMVA::Event*, LossFunctionEventInfo>                fLossFunctionEventInfo;  // map event to true value, predicted value, and weight
                                                                                                   // used by different loss functions for BDT regression
      std::map< const TMVA::Event*,std::vector<double> > fResiduals; // individual event residuals for gradient boost

      //options for the decision Tree
      SeparationBase                 *fSepType;         // the separation used in node splitting
      TString                         fSepTypeS;        // the separation (option string) used in node splitting
      Int_t                           fMinNodeEvents;   // min number of events in node
      Float_t                         fMinNodeSize;     // min percentage of training events in node
      TString                         fMinNodeSizeS;    // string containing min percentage of training events in node

      Int_t                           fNCuts;           // grid used in cut applied in node splitting
      Bool_t                          fUseFisherCuts;   // use multivariate splits using the Fisher criterium
      Double_t                        fMinLinCorrForFisher; // the minimum linear correlation between two variables demanded for use in fisher criterium in node splitting
      Bool_t                          fUseExclusiveVars; // individual variables already used in fisher criterium are not anymore analysed individually for node splitting
      Bool_t                          fUseYesNoLeaf;    // use sig or bkg classification in leave nodes or sig/bkg
      Double_t                        fNodePurityLimit; // purity limit for sig/bkg nodes
      UInt_t                          fNNodesMax;       // max # of nodes
      UInt_t                          fMaxDepth;        // max depth

      DecisionTree::EPruneMethod       fPruneMethod;     // method used for prunig
      TString                          fPruneMethodS;    // prune method option String
      Double_t                         fPruneStrength;   // a parameter to set the "amount" of pruning..needs to be adjusted
      Double_t                         fFValidationEvents;    // fraction of events to use for pruning
      Bool_t                           fAutomatic;       // use user given prune strength or automatically determined one using a validation sample
      Bool_t                           fRandomisedTrees; // choose a random subset of possible cut variables at each node during training
      UInt_t                           fUseNvars;        // the number of variables used in the randomised tree splitting
      Bool_t                           fUsePoissonNvars; // use "fUseNvars" not as fixed number but as mean of a possion distr. in each split
      UInt_t                           fUseNTrainEvents; // number of randomly picked training events used in randomised (and bagged) trees

      Double_t                         fBaggedSampleFraction;     // relative size of bagged event sample to original sample size
      TString                          fNegWeightTreatment;     // variable that holds the option of how to treat negative event weights in training
      Bool_t                           fNoNegWeightsInTraining; // ignore negative event weights in the training
      Bool_t                           fInverseBoostNegWeights; // boost ev. with neg. weights with 1/boostweight rathre than boostweight
      Bool_t                           fPairNegWeightsGlobal;   // pair ev. with neg. and pos. weights in traning sample and "annihilate" them 
      Bool_t                           fTrainWithNegWeights; // yes there are negative event weights and we don't ignore them
      Bool_t                           fDoBoostMonitor; //create control plot with ROC integral vs tree number


      //some histograms for monitoring
      TTree*                           fMonitorNtuple;   // monitoring ntuple
      Int_t                            fITree;           // ntuple var: ith tree
      Double_t                         fBoostWeight;     // ntuple var: boost weight
      Double_t                         fErrorFraction;   // ntuple var: misclassification error fraction

      Double_t                         fCss;             // Cost factor
      Double_t                         fCts_sb;          // Cost factor
      Double_t                         fCtb_ss;          // Cost factor
      Double_t                         fCbb;             // Cost factor
      
      Bool_t                           fDoPreselection;  // do or do not perform automatic pre-selection of 100% eff. cuts

      Bool_t                           fSkipNormalization; // true for skipping normalization at initialization of trees

      std::vector<Double_t>            fVariableImportance; // the relative importance of the different variables


      void                             DeterminePreselectionCuts(const std::vector<const TMVA::Event*>& eventSample);
      Double_t                         ApplyPreselectionCuts(const Event* ev);
      
      std::vector<Double_t> fLowSigCut;
      std::vector<Double_t> fLowBkgCut;
      std::vector<Double_t> fHighSigCut;
      std::vector<Double_t> fHighBkgCut;
      
      std::vector<Bool_t>  fIsLowSigCut;  
      std::vector<Bool_t>  fIsLowBkgCut;  
      std::vector<Bool_t>  fIsHighSigCut; 
      std::vector<Bool_t>  fIsHighBkgCut; 
      
      Bool_t fHistoricBool; //historic variable, only needed for "CompatibilityOptions" 

      TString                         fRegressionLossFunctionBDTGS;       // the option string determining the loss function for BDT regression
      Double_t                        fHuberQuantile;                     // the option string determining the quantile for the Huber Loss Function
                                                                          // in BDT regression.
      LossFunctionBDT* fRegressionLossFunctionBDTG;

      // debugging flags
      static const Int_t               fgDebugLevel;     // debug level determining some printout/control plots etc.

      // for backward compatibility
      ClassDef(MethodBDT,0);  // Analysis of Boosted Decision Trees
   };

} // namespace TMVA

const std::vector<TMVA::DecisionTree*>& TMVA::MethodBDT::GetForest()         const { return fForest; }
const std::vector<const TMVA::Event*> & TMVA::MethodBDT::GetTrainingEvents() const { return fEventSample; }
const std::vector<double>&              TMVA::MethodBDT::GetBoostWeights()   const { return fBoostWeights; }

#endif
