// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Matt Jachowski, Peter Speckmayer, Helge Voss, Kai Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : MethodCuts                                                            *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Multivariate optimisation of signal efficiency for given background       *
 *      efficiency, using rectangular minimum and maximum requirements on         *
 *      input variables                                                           *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker  <Andreas.Hocker@cern.ch> - CERN, Switzerland             *
 *      Matt Jachowski   <jachowski@stanford.edu> - Stanford University, USA      *
 *      Peter Speckmayer <speckmay@mail.cern.ch>  - CERN, Switzerland             *
 *      Helge Voss       <Helge.Voss@cern.ch>     - MPI-K Heidelberg, Germany     *
 *      Kai Voss         <Kai.Voss@cern.ch>       - U. of Victoria, Canada        *
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         *
 *      U. of Victoria, Canada                                                    *
 *      MPI-K Heidelberg, Germany                                                 *
 *      LAPP, Annecy, France                                                      *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

#ifndef ROOT_TMVA_MethodCuts
#define ROOT_TMVA_MethodCuts

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// MethodCuts                                                           //
//                                                                      //
// Multivariate optimisation of signal efficiency for given background  //
// efficiency, using rectangular minimum and maximum requirements on    //
// input variables                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <vector>


#include "TMVA/MethodBase.h"
#include "TMVA/BinarySearchTree.h"
#include "TMVA/PDF.h"
#include "TMatrixDfwd.h"
#include "IFitterTarget.h"

class TRandom;

namespace TMVA {

   class Interval;

   class MethodCuts : public MethodBase, public IFitterTarget {

   public:

      MethodCuts( const TString& jobName,
                  const TString& methodTitle, 
                  DataSetInfo& theData,
                  const TString& theOption = "MC:150:10000:");

      MethodCuts( DataSetInfo& theData,
                  const TString& theWeightFile);

      // this is a workaround which is necessary since CINT is not capable of handling dynamic casts
      static MethodCuts* DynamicCast( IMethod* method ) { return dynamic_cast<MethodCuts*>(method); }

      virtual ~MethodCuts( void );

      virtual Bool_t HasAnalysisType( Types::EAnalysisType type, UInt_t numberClasses, UInt_t numberTargets );

      // training method
      void Train( void );

      using MethodBase::ReadWeightsFromStream;

      void AddWeightsXMLTo      ( void* parent ) const;

      void ReadWeightsFromStream( std::istream & i );
      void ReadWeightsFromXML   ( void* wghtnode );

      // calculate the MVA value (for CUTs this is just a dummy)
      Double_t GetMvaValue( Double_t* err = 0, Double_t* errUpper = 0 );

      // write method specific histos to target file
      void WriteMonitoringHistosToFile( void ) const;

      // test the method
      void TestClassification();
     
      // also overwrite --> not computed for cuts
      Double_t GetSeparation  ( TH1*, TH1* ) const { return -1; }
      Double_t GetSeparation  ( PDF* = 0, PDF* = 0 ) const { return -1; }
      Double_t GetSignificance( void )       const { return -1; }
      Double_t GetmuTransform ( TTree *)           { return -1; }
      Double_t GetEfficiency  ( const TString&, Types::ETreeType, Double_t& );
      Double_t GetTrainingEfficiency(const TString& );

      // rarity distributions (signal or background (default) is uniform in [0,1])
      Double_t GetRarity( Double_t, Types::ESBType ) const { return 0; }

      // accessors for Minuit
      Double_t ComputeEstimator( std::vector<Double_t> & );
      
      Double_t EstimatorFunction( std::vector<Double_t> & );
      Double_t EstimatorFunction( Int_t ievt1, Int_t ievt2 );

      void     SetTestSignalEfficiency( Double_t effS ) { fTestSignalEff = effS; }

      // retrieve cut values for given signal efficiency
      void     PrintCuts( Double_t effS ) const;
      Double_t GetCuts  ( Double_t effS, std::vector<Double_t>& cutMin, std::vector<Double_t>& cutMax ) const;
      Double_t GetCuts  ( Double_t effS, Double_t* cutMin, Double_t* cutMax ) const;

      // ranking of input variables (not available for cuts)
      const Ranking* CreateRanking() { return 0; }

      void DeclareOptions();
      void ProcessOptions();

      // maximum |cut| value
      static const Double_t fgMaxAbsCutVal;

      // no check of options at this place
      void CheckSetup() {}

   protected:

      // make ROOT-independent C++ class for classifier response (classifier-specific implementation)
      void MakeClassSpecific( std::ostream&, const TString& ) const;

      // get help message text
      void GetHelpMessage() const;

   private:

      // optimisation method
      enum EFitMethodType { kUseMonteCarlo = 0,
                            kUseGeneticAlgorithm,
                            kUseSimulatedAnnealing,
                            kUseMinuit,
                            kUseEventScan,
                            kUseMonteCarloEvents };

      // efficiency calculation method
      // - kUseEventSelection: computes efficiencies from given data sample
      // - kUsePDFs          : creates smoothed PDFs from data samples, and 
      //                       uses this to compute efficiencies
      enum EEffMethod     { kUseEventSelection = 0,
                            kUsePDFs };

      // improve the Monte Carlo by providing some additional information
      enum EFitParameters { kNotEnforced = 0,
                            kForceMin,
                            kForceMax,
                            kForceSmart };

      // general
      TString                 fFitMethodS;         // chosen fit method (string)
      EFitMethodType          fFitMethod;          // chosen fit method
      TString                 fEffMethodS;         // chosen efficiency calculation method (string)
      EEffMethod              fEffMethod;          // chosen efficiency calculation method
      std::vector<EFitParameters>* fFitParams;     // vector for series of fit methods
      Double_t                fTestSignalEff;      // used to test optimized signal efficiency
      Double_t                fEffSMin;            // used to test optimized signal efficiency
      Double_t                fEffSMax;            // used to test optimized signal efficiency
      Double_t*               fCutRangeMin;        // minimum of allowed cut range
      Double_t*               fCutRangeMax;        // maximum of allowed cut range
      std::vector<Interval*>  fCutRange;           // allowed ranges for cut optimisation

      // for the use of the binary tree method
      BinarySearchTree*       fBinaryTreeS;
      BinarySearchTree*       fBinaryTreeB;

      // MC method
      Double_t**              fCutMin;             // minimum requirement
      Double_t**              fCutMax;             // maximum requirement
      Double_t*               fTmpCutMin;          // temporary minimum requirement
      Double_t*               fTmpCutMax;          // temporary maximum requirement
      TString*                fAllVarsI;           // what to do with variables

      // relevant for all methods
      Int_t                   fNpar;               // number of parameters in fit (default: 2*Nvar)
      Double_t                fEffRef;             // reference efficiency
      std::vector<Int_t>*     fRangeSign;          // used to match cuts to fit parameters (and vice versa)
      TRandom*                fRandom;             // random generator for MC optimisation method

      // basic statistics
      std::vector<Double_t>*  fMeanS;              // means of variables (signal)
      std::vector<Double_t>*  fMeanB;              // means of variables (background)
      std::vector<Double_t>*  fRmsS;               // RMSs of variables (signal)
      std::vector<Double_t>*  fRmsB;               // RMSs of variables (background)

      TH1*                    fEffBvsSLocal;       // intermediate eff. background versus eff signal histo

      // PDF section
      std::vector<TH1*>*      fVarHistS;           // reference histograms (signal)
      std::vector<TH1*>*      fVarHistB;           // reference histograms (background)
      std::vector<TH1*>*      fVarHistS_smooth;    // smoothed reference histograms (signal)        
      std::vector<TH1*>*      fVarHistB_smooth;    // smoothed reference histograms (background)
      std::vector<PDF*>*      fVarPdfS;            // reference PDFs (signal)
      std::vector<PDF*>*      fVarPdfB;            // reference PDFs (background)

      // negative efficiencies
      Bool_t                  fNegEffWarning;      // flag risen in case of negative efficiency warning


      // the definition of fit parameters can be different from the actual 
      // cut requirements; these functions provide the matching
      void     MatchParsToCuts( const std::vector<Double_t>&, Double_t*, Double_t* );
      void     MatchParsToCuts( Double_t*, Double_t*, Double_t* );

      void     MatchCutsToPars( std::vector<Double_t>&, Double_t*, Double_t* );
      void     MatchCutsToPars( std::vector<Double_t>&, Double_t**, Double_t**, Int_t ibin );

      // creates PDFs in case these are used to compute efficiencies 
      // (corresponds to: EffMethod == kUsePDFs)
      void     CreateVariablePDFs( void );

      // returns signal and background efficiencies for given cuts - using event counting
      void     GetEffsfromSelection( Double_t* cutMin, Double_t* cutMax,
                                     Double_t& effS, Double_t& effB );
      // returns signal and background efficiencies for given cuts - using PDFs
      void     GetEffsfromPDFs( Double_t* cutMin, Double_t* cutMax,
                                Double_t& effS, Double_t& effB );

      // default initialisation method called by all constructors
      void     Init( void );

      ClassDef(MethodCuts,0);  // Multivariate optimisation of signal efficiency
   };

} // namespace TMVA

#endif
