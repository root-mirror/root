// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Joerg Stelzer, Fredrik Tegenfeldt, Helge Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : RuleFitParams                                                         *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      A class doing the actual fitting of a linear model using rules as         *
 *      base functions.                                                           *
 *      Reference paper: 1.Gradient Directed Regularization                       *
 *                         Friedman, Popescu, 2004                                *
 *                       2.Predictive Learning with Rule Ensembles                *
 *                         Friedman, Popescu, 2005                                *
 *                                                                                *
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

#ifndef ROOT_TMVA_RuleFitParams
#define ROOT_TMVA_RuleFitParams

#include "TMathBase.h"

#include "TMVA/Event.h"

#include <vector>

class TTree;

namespace TMVA {

   class RuleEnsemble;
   class MsgLogger;
   class RuleFit;
   class RuleFitParams {

   public:

      RuleFitParams();
      virtual ~RuleFitParams();

      void Init();

      // set message type
      void SetMsgType( EMsgType t );

      // set RuleFit ptr
      void SetRuleFit( RuleFit *rf )    { fRuleFit = rf; }
      //
      // GD path: set N(path steps)
      void SetGDNPathSteps( Int_t np )  { fGDNPathSteps = np; }

      // GD path: set path step size
      void SetGDPathStep( Double_t s )  { fGDPathStep = s; }

      // GD path: set tau search range
      void SetGDTauRange( Double_t t0, Double_t t1 )
      {
         fGDTauMin = (t0>1.0 ? 1.0:(t0<0.0 ? 0.0:t0));
         fGDTauMax = (t1>1.0 ? 1.0:(t1<0.0 ? 0.0:t1));
         if (fGDTauMax<fGDTauMin) fGDTauMax = fGDTauMin;
      }

      // GD path: set number of steps in tau search range
      void SetGDTauScan( UInt_t n )        { fGDTauScan = n; }

      // GD path: set tau
      void SetGDTau( Double_t t ) { fGDTau = t; }


      void SetGDErrScale( Double_t s ) { fGDErrScale = s; }
      void SetGDTauPrec( Double_t p )  { fGDTauPrec=p; CalcGDNTau(); fGDTauVec.resize(fGDNTau); }

      // return type such that +1 = signal and -1 = background
      Int_t Type( const Event * e ) const; // return (fRuleFit->GetMethodRuleFit()->DataInfo().IsSignal(e) ? 1:-1); }
      //
      UInt_t                            GetPathIdx1() const { return fPathIdx1; }
      UInt_t                            GetPathIdx2() const { return fPathIdx2; }
      UInt_t                            GetPerfIdx1() const { return fPerfIdx1; }
      UInt_t                            GetPerfIdx2() const { return fPerfIdx2; }

      // Loss function; Huber loss eq 33
      Double_t LossFunction( const Event& e ) const;

      // same but using evt idx (faster)
      Double_t LossFunction( UInt_t evtidx ) const;
      Double_t LossFunction( UInt_t evtidx, UInt_t itau ) const;

      // Empirical risk
      Double_t Risk(UInt_t ind1, UInt_t ind2, Double_t neff) const;
      Double_t Risk(UInt_t ind1, UInt_t ind2, Double_t neff, UInt_t itau) const;

      // Risk evaluation for fPathIdx and fPerfInd
      Double_t RiskPath() const { return Risk(fPathIdx1,fPathIdx2,fNEveEffPath); }
      Double_t RiskPerf() const { return Risk(fPerfIdx1,fPerfIdx2,fNEveEffPerf); }
      Double_t RiskPerf( UInt_t itau ) const { return Risk(fPerfIdx1,fPerfIdx2,fNEveEffPerf,itau); }

      // Risk evaluation for all tau
      UInt_t RiskPerfTst();

      // Penalty function; Lasso function (eq 8)
      Double_t Penalty() const;

      // initialize GD path
      void InitGD();

      // find best tau and return the number of scan steps used
      Int_t FindGDTau();

      // make path for binary classification (squared-error ramp, sect 6 in ref 1)
      void MakeGDPath();

   protected:

      // typedef of an Event const iterator
      typedef std::vector<const TMVA::Event *>::const_iterator  EventItr;

      // init ntuple
      void InitNtuple();

      // calculate N(tau) in scan - limit to 100000.
      void CalcGDNTau()  { fGDNTau = static_cast<UInt_t>(1.0/fGDTauPrec)+1; if (fGDNTau>100000) fGDNTau=100000; }

      // fill ntuple with coefficient info
      void FillCoefficients();

      // estimate the optimum scoring function
      void CalcFStar();

      // estimate of binary error rate
      Double_t ErrorRateBin();

      // estimate of scale average error rate
      Double_t ErrorRateReg();

      // estimate 1-area under ROC
      Double_t ErrorRateRocRaw( std::vector<Double_t> & sFsig, std::vector<Double_t> & sFbkg );
      Double_t ErrorRateRoc();
      void     ErrorRateRocTst();

      // estimate optimism
      Double_t Optimism();

      // make gradient vector (eq 44 in ref 1)
      void MakeGradientVector();

      // Calculate the direction in parameter space (eq 25, ref 1) and update coeffs (eq 22, ref 1)
      void UpdateCoefficients();

      // calculate average of responses of F
      Double_t CalcAverageResponse();
      Double_t CalcAverageResponseOLD();

      // calculate average of true response (initial estimate of a0)
      Double_t CalcAverageTruth();

      // calculate the average of each variable over the range
      void EvaluateAverage(UInt_t ind1, UInt_t ind2,
                           std::vector<Double_t> &avsel,
                           std::vector<Double_t> &avrul);

      // evaluate using fPathIdx1,2
      void EvaluateAveragePath() { EvaluateAverage( fPathIdx1, fPathIdx2, fAverageSelectorPath, fAverageRulePath ); }

      // evaluate using fPerfIdx1,2
      void EvaluateAveragePerf() { EvaluateAverage( fPerfIdx1, fPerfIdx2, fAverageSelectorPerf, fAverageRulePerf ); }

      // the same as above but for the various tau
      void MakeTstGradientVector();
      void UpdateTstCoefficients();
      void CalcTstAverageResponse();


      RuleFit             * fRuleFit;      // rule fit
      RuleEnsemble        * fRuleEnsemble; // rule ensemble
      //
      UInt_t                fNRules;       // number of rules
      UInt_t                fNLinear;      // number of linear terms
      //
      // Event indices for path/validation - TODO: should let the user decide
      // Now it is just a simple one-fold cross validation.
      //
      UInt_t                fPathIdx1;       // first event index for path search
      UInt_t                fPathIdx2;       // last event index for path search
      UInt_t                fPerfIdx1;       // first event index for performance evaluation
      UInt_t                fPerfIdx2;       // last event index for performance evaluation
      Double_t              fNEveEffPath;    // sum of weights for Path events
      Double_t              fNEveEffPerf;    // idem for Perf events

      std::vector<Double_t> fAverageSelectorPath; // average of each variable over the range fPathIdx1,2
      std::vector<Double_t> fAverageRulePath;     // average of each rule, same range
      std::vector<Double_t> fAverageSelectorPerf; // average of each variable over the range fPerfIdx1,2
      std::vector<Double_t> fAverageRulePerf;     // average of each rule, same range

      std::vector<Double_t> fGradVec;        // gradient vector - dimension = number of rules in ensemble
      std::vector<Double_t> fGradVecLin;     // gradient vector - dimension = number of variables

      std::vector< std::vector<Double_t> > fGradVecTst;    // gradient vector - one per tau
      std::vector< std::vector<Double_t> > fGradVecLinTst; // gradient vector, linear terms - one per tau
      //
      std::vector<Double_t> fGDErrTst;     // error rates per tau
      std::vector<Char_t>   fGDErrTstOK;   // error rate is sufficiently low <--- stores boolean
      std::vector< std::vector<Double_t> > fGDCoefTst;    // rule coeffs - one per tau
      std::vector< std::vector<Double_t> > fGDCoefLinTst; // linear coeffs - one per tau
      std::vector<Double_t> fGDOfsTst;       // offset per tau
      std::vector< Double_t > fGDTauVec;     // the tau's
      UInt_t                fGDNTauTstOK;    // number of tau in the test-phase that are ok
      UInt_t                fGDNTau;         // number of tau-paths - calculated in SetGDTauPrec
      Double_t              fGDTauPrec;      // precision in tau
      UInt_t                fGDTauScan;      // number scan for tau-paths
      Double_t              fGDTauMin;       // min threshold parameter (tau in eq 26, ref 1)
      Double_t              fGDTauMax;       // max threshold parameter (tau in eq 26, ref 1)
      Double_t              fGDTau;          // selected threshold parameter (tau in eq 26, ref 1)
      Double_t              fGDPathStep;     // step size along path (delta nu in eq 22, ref 1)
      Int_t                 fGDNPathSteps;   // number of path steps
      Double_t              fGDErrScale;     // stop scan at error = scale*errmin
      //
      Double_t              fAverageTruth;   // average truth, ie sum(y)/N, y=+-1
      //
      std::vector<Double_t> fFstar;          // vector of F*() - filled in CalcFStar()
      Double_t              fFstarMedian;    // median value of F*() using
      //
      TTree                *fGDNtuple;       // Gradient path ntuple, contains params for each step along the path
      Double_t              fNTRisk;         // GD path: risk
      Double_t              fNTErrorRate;    // GD path: error rate (or performance)
      Double_t              fNTNuval;        // GD path: value of nu
      Double_t              fNTCoefRad;      // GD path: 'radius' of all rulecoeffs
      Double_t              fNTOffset;       // GD path: model offset
      Double_t             *fNTCoeff;        // GD path: rule coefficients
      Double_t             *fNTLinCoeff;     // GD path: linear coefficients

      Double_t              fsigave;         // Sigma of current signal score function F(sig)
      Double_t              fsigrms;         // Rms of F(sig)
      Double_t              fbkgave;         // Average of F(bkg)
      Double_t              fbkgrms;         // Rms of F(bkg)

   private:

      mutable MsgLogger*    fLogger;         //! message logger
      MsgLogger& Log() const { return *fLogger; }

   };

   // --------------------------------------------------------

   class AbsValue {

   public:

      Bool_t operator()( Double_t first, Double_t second ) const { return TMath::Abs(first) < TMath::Abs(second); }
   };
}


#endif
