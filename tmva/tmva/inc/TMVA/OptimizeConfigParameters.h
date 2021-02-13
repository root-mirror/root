/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : OptimizeConfigParameters                                              *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description: The OptimizeConfigParameters takes care of "scanning/fitting"     *
 *              different tuning parameters in order to find the best set of      *
 *              tuning paraemters which will be used in the end                   *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Helge Voss      <Helge.Voss@cern.ch>     - MPI-K Heidelberg, Germany      *
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         *
 *      MPI-K Heidelberg, Germany                                                 *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://ttmva.sourceforge.net/LICENSE)                                         *
 **********************************************************************************/

#ifndef ROOT_TMVA_OptimizeConfigParameters
#define ROOT_TMVA_OptimizeConfigParameters

#include "Rtypes.h"

#include "TString.h"

#include "TMVA/MethodBase.h"

#include "TMVA/Interval.h"

#include "TMVA/DataSet.h"

#include "IFitterTarget.h"

#include "TH1.h"

#include <map>
#include <vector>

class TestOptimizeConfigParameters;

namespace TMVA {

   class MethodBase;
   class MsgLogger;
   class OptimizeConfigParameters : public IFitterTarget  {

   public:
      friend TestOptimizeConfigParameters;

      //default constructor
      OptimizeConfigParameters(MethodBase * const method, std::map<TString,TMVA::Interval*> tuneParameters, TString fomType="Separation", TString optimizationType = "GA");

      // destructor
      virtual ~OptimizeConfigParameters();
      // could later be changed to be set via option string...
      // but for now it's simpler like this
      std::map<TString,Double_t> optimize();

   private:
      std::vector< int > GetScanIndices( int val, std::vector<int> base);
      void optimizeScan();
      void optimizeFit();

      Double_t EstimatorFunction( std::vector<Double_t> & );

      Double_t GetFOM();

      MethodBase* GetMethod(){return fMethod;}

      void GetMVADists();
      Double_t GetSeparation();
      Double_t GetROCIntegral();
      Double_t GetSigEffAtBkgEff( Double_t bkgEff = 0.1);
      Double_t GetBkgEffAtSigEff( Double_t sigEff = 0.5);
      Double_t GetBkgRejAtSigEff( Double_t sigEff = 0.5);


      MethodBase* const fMethod; // The MVA method to be evaluated
      std::vector<Float_t>             fFOMvsIter; // graph showing the development of the Figure Of Merit values during the fit
      std::map<TString,TMVA::Interval*> fTuneParameters; // parameters included in the tuning
      std::map<TString,Double_t>       fTunedParameters; // parameters included in the tuning
      std::map< std::vector<Double_t> , Double_t>  fAlreadyTrainedParCombination; // save parameters for which the FOM is already known (GA seems to evaluate the same parameters several times)
      TString           fFOMType;    // the FOM type (Separation, ROC integra.. whatever you implemented..
      TString           fOptimizationFitType; // which type of optimisation procedure to be used
      TH1D             *fMvaSig; // MVA distribution for signal events, used for spline fit
      TH1D             *fMvaBkg; // MVA distribution for bakgr. events, used for spline fit

      TH1D             *fMvaSigFineBin; // MVA distribution for signal events
      TH1D             *fMvaBkgFineBin; // MVA distribution for bakgr. events

      Bool_t           fNotDoneYet; // flat to indicate of Method Transformations have been obtained yet or not (normally done in MethodBase::TrainMethod)

      mutable MsgLogger*         fLogger;   // message logger
      MsgLogger& Log() const { return *fLogger; }

      ClassDef(OptimizeConfigParameters,0); // Interface to different separation criteria used in training algorithms
   };
} // namespace TMVA

#endif
