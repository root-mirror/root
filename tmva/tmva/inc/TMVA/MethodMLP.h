// @(#)root/tmva $Id$
// Author: Krzysztof Danielowski, Andreas Hoecker, Matt Jachowski, Kamil Kraszewski, Maciej Kruk, Peter Speckmayer, Joerg Stelzer, Eckhard von Toerne, Jan Therhaag, Jiahang Zhong

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : MethodMLP                                                             *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      ANN Multilayer Perceptron  class for the discrimination of signal         *
 *      from background.  BFGS implementation based on TMultiLayerPerceptron      *
 *      class from ROOT (http://root.cern.ch).                                    *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Krzysztof Danielowski <danielow@cern.ch>       - IFJ & AGH, Poland        *
 *      Andreas Hoecker       <Andreas.Hocker@cern.ch> - CERN, Switzerland        *
 *      Matt Jachowski        <jachowski@stanford.edu> - Stanford University, USA *
 *      Kamil Kraszewski      <kalq@cern.ch>           - IFJ & UJ, Poland         *
 *      Maciej Kruk           <mkruk@cern.ch>          - IFJ & AGH, Poland        *
 *      Peter Speckmayer      <peter.speckmayer@cern.ch> - CERN, Switzerland      *
 *      Joerg Stelzer         <stelzer@cern.ch>        - DESY, Germany            *
 *      Jan Therhaag          <Jan.Therhaag@cern.ch>   - U of Bonn, Germany       *
 *      Eckhard v. Toerne     <evt@uni-bonn.de>        - U of Bonn, Germany       *
 *      Jiahang Zhong         <Jiahang.Zhong@cern.ch>  - Academia Sinica, Taipei  *
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

#ifndef ROOT_TMVA_MethodMLP
#define ROOT_TMVA_MethodMLP

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// MethodMLP                                                            //
//                                                                      //
// Multilayer Perceptron built off of MethodANNBase                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <vector>
#include <utility>
#include "TString.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TH1F.h"
#include "TMatrixDfwd.h"

#include "TMVA/IFitterTarget.h"
#include "TMVA/MethodBase.h"
#include "TMVA/MethodANNBase.h"
#include "TMVA/TNeuron.h"
#include "TMVA/TActivation.h"
#include "TMVA/ConvergenceTest.h"

#define MethodMLP_UseMinuit__
#undef  MethodMLP_UseMinuit__

namespace TMVA {

   class MethodMLP : public MethodANNBase, public IFitterTarget, public ConvergenceTest {

   public:

      // standard constructors
      MethodMLP( const TString& jobName,
                 const TString&  methodTitle,
                 DataSetInfo& theData,
                 const TString& theOption );

      MethodMLP( DataSetInfo& theData,
                 const TString& theWeightFile );

      virtual ~MethodMLP();

      virtual Bool_t HasAnalysisType( Types::EAnalysisType type, UInt_t numberClasses, UInt_t numberTargets );

      void Train();
      // for GA
      Double_t ComputeEstimator ( std::vector<Double_t>& parameters );
      Double_t EstimatorFunction( std::vector<Double_t>& parameters );

      enum ETrainingMethod { kBP=0, kBFGS, kGA };
      enum EBPTrainingMode { kSequential=0, kBatch };

      bool     HasInverseHessian() { return fCalculateErrors; }
      Double_t GetMvaValue( Double_t* err=0, Double_t* errUpper=0 );

   protected:

      // make ROOT-independent C++ class for classifier response (classifier-specific implementation)
      void MakeClassSpecific( std::ostream&, const TString& ) const;

      // get help message text
      void GetHelpMessage() const;


   private:

      // the option handling methods
      void DeclareOptions();
      void ProcessOptions();

      // general helper functions
      void     Train( Int_t nEpochs );
      void     Init();
      void     InitializeLearningRates(); // although this is only needed by backprop

      // used as a measure of success in all minimization techniques
      Double_t CalculateEstimator( Types::ETreeType treeType = Types::kTraining, Int_t iEpoch = -1 );

      // BFGS functions
      void     BFGSMinimize( Int_t nEpochs );
      void     SetGammaDelta( TMatrixD &Gamma, TMatrixD &Delta, std::vector<Double_t> &Buffer );
      void     SteepestDir( TMatrixD &Dir );
      Bool_t   GetHessian( TMatrixD &Hessian, TMatrixD &Gamma, TMatrixD &Delta );
      void     SetDir( TMatrixD &Hessian, TMatrixD &Dir );
      Double_t DerivDir( TMatrixD &Dir );
      Bool_t   LineSearch( TMatrixD &Dir, std::vector<Double_t> &Buffer, Double_t* dError=0 ); //zjh
      void     ComputeDEDw();
      void     SimulateEvent( const Event* ev );
      void     SetDirWeights( std::vector<Double_t> &Origin, TMatrixD &Dir, Double_t alpha );
      Double_t GetError();
      Double_t GetMSEErr( const Event* ev, UInt_t index = 0 );   //zjh
      Double_t GetCEErr( const Event* ev, UInt_t index = 0 );   //zjh

      // backpropagation functions
      void     BackPropagationMinimize( Int_t nEpochs );
      void     TrainOneEpoch();
      void     Shuffle( Int_t* index, Int_t n );
      void     DecaySynapseWeights(Bool_t lateEpoch );
      void     TrainOneEvent( Int_t ievt);
      Double_t GetDesiredOutput( const Event* ev );
      void     UpdateNetwork( Double_t desired, Double_t eventWeight=1.0 );
      void     UpdateNetwork(const std::vector<Float_t>& desired, Double_t eventWeight=1.0);
      void     CalculateNeuronDeltas();
      void     UpdateSynapses();
      void     AdjustSynapseWeights();

      // faster backpropagation
      void     TrainOneEventFast( Int_t ievt, Float_t*& branchVar, Int_t& type );

      // genetic algorithm functions
      void GeneticMinimize();


#ifdef MethodMLP_UseMinuit__
      // minuit functions -- commented out because they rely on a static pointer
      void MinuitMinimize();
      static MethodMLP* GetThisPtr();
      static void IFCN( Int_t& npars, Double_t* grad, Double_t &f, Double_t* fitPars, Int_t ifl );
      void FCN( Int_t& npars, Double_t* grad, Double_t &f, Double_t* fitPars, Int_t ifl );
#endif

      // general
      bool               fUseRegulator;         // zjh
      bool               fCalculateErrors;      // compute inverse hessian matrix at the end of the training
      Double_t           fPrior;                // zjh
      std::vector<Double_t> fPriorDev;          // zjh
      void               GetApproxInvHessian ( TMatrixD& InvHessian, bool regulate=true );   //rank-1 approximation, neglect 2nd derivatives. //zjh
      void               UpdateRegulators();    // zjh
      void               UpdatePriors();        // zjh
      Int_t              fUpdateLimit;          // zjh

      ETrainingMethod fTrainingMethod; // method of training, BP or GA
      TString         fTrainMethodS;   // training method option param

      Float_t         fSamplingFraction;  // fraction of events which is sampled for training
      Float_t         fSamplingEpoch;     // fraction of epochs where sampling is used
      Float_t         fSamplingWeight;    // changing factor for event weights when sampling is turned on
      Bool_t          fSamplingTraining;  // The training sample is sampled
      Bool_t          fSamplingTesting;   // The testing sample is sampled

      // BFGS variables
      Double_t        fLastAlpha;      // line search variable
      Double_t        fTau;            // line search variable
      Int_t           fResetStep;      // reset time (how often we clear hessian matrix)

      // back propagation variable
      Double_t        fLearnRate;      // learning rate for synapse weight adjustments
      Double_t        fDecayRate;      // decay rate for above learning rate
      EBPTrainingMode fBPMode;         // backprop learning mode (sequential or batch)
      TString         fBpModeS;        // backprop learning mode option string (sequential or batch)
      Int_t           fBatchSize;      // batch size, only matters if in batch learning mode
      Int_t           fTestRate;       // test for overtraining performed at each #th epochs
      Bool_t          fEpochMon;       // create and fill epoch-wise monitoring histograms (makes outputfile big!)

      // genetic algorithm variables
      Int_t           fGA_nsteps;      // GA settings: number of steps
      Int_t           fGA_preCalc;     // GA settings: number of pre-calc steps
      Int_t           fGA_SC_steps;    // GA settings: SC_steps
      Int_t           fGA_SC_rate; // GA settings: SC_rate
      Double_t        fGA_SC_factor;   // GA settings: SC_factor

      // regression, storage of deviations
      std::vector<std::pair<Float_t,Float_t> >* fDeviationsFromTargets; // deviation from the targets, event weight

      Float_t         fWeightRange;    // suppress outliers for the estimator calculation

#ifdef MethodMLP_UseMinuit__
      // minuit variables -- commented out because they rely on a static pointer
      Int_t          fNumberOfWeights; // Minuit: number of weights
      static MethodMLP* fgThis;        // Minuit: this pointer
#endif

      // debugging flags
      static const Int_t  fgPRINT_ESTIMATOR_INC = 10;     // debug flags
      static const Bool_t fgPRINT_SEQ           = kFALSE; // debug flags
      static const Bool_t fgPRINT_BATCH         = kFALSE; // debug flags

      ClassDef(MethodMLP,0); // Multi-layer perceptron implemented specifically for TMVA
   };

} // namespace TMVA

#endif
