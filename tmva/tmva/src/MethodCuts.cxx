// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Matt Jachowski, Peter Speckmayer, Eckhard von Toerne, Helge Voss, Kai Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate Data analysis       *
 * Package: TMVA                                                                  *
 * Class  : TMVA::MethodCuts                                                      *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Implementation (see header for description)                               *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
 *      Matt Jachowski  <jachowski@stanford.edu> - Stanford University, USA       *
 *      Peter Speckmayer <speckmay@mail.cern.ch> - CERN, Switzerland              *
 *      Eckhard von Toerne <evt@physik.uni-bonn.de> - U. of Bonn, Germany         *
 *      Helge Voss      <Helge.Voss@cern.ch>     - MPI-K Heidelberg, Germany      *
 *      Kai Voss        <Kai.Voss@cern.ch>       - U. of Victoria, Canada         *
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         *
 *      U. of Victoria, Canada                                                    *
 *      MPI-K Heidelberg, Germany                                                 *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

/*! \class TMVA::MethodCuts
\ingroup TMVA

   Multivariate optimisation of signal efficiency for given background
   efficiency, applying rectangular minimum and maximum requirements.

   Also implemented is a "decorrelate/diagonalized cuts approach",
   which improves over the uncorrelated cuts approach by
   transforming linearly the input variables into a diagonal space,
   using the square-root of the covariance matrix.

   Other optimisation criteria, such as maximising the signal significance-
   squared, \f$ \frac{S^2}{(S+B)} \f$, with S and B being the signal and background yields,
   correspond to a particular point in the optimised background rejection
   versus signal efficiency curve. This working point requires the knowledge
   of the expected yields, which is not the case in general. Note also that
   for rare signals, Poissonian statistics should be used, which modifies
   the significance criterion.

   The rectangular cut of a volume in the variable space is performed using
   a binary tree to sort the training events. This provides a significant
   reduction in computing time (up to several orders of magnitudes, depending
   on the complexity of the problem at hand).

   Technically, optimisation is achieved in TMVA by two methods:

   1. Monte Carlo generation using uniform priors for the lower cut value,
   and the cut width, thrown within the variable ranges.

   2. A Genetic Algorithm (GA) searches for the optimal ("fittest") cut sample.
   The GA is configurable by many external settings through the option
   string. For difficult cases (such as many variables), some tuning
   may be necessary to achieve satisfying results

   Attempts to use Minuit fits (Simplex ot Migrad) instead have not shown
   superior results, and often failed due to convergence at local minima.

   The tests we have performed so far showed that in generic applications,
   the GA is superior to MC sampling, and hence GA is the default method.
   It is worthwhile trying both anyway.

   **Decorrelated (or "diagonalized") Cuts**

   See class description for Method Likelihood for a detailed explanation.
*/

#include "TMVA/MethodCuts.h"

#include "TMVA/BinarySearchTree.h"
#include "TMVA/ClassifierFactory.h"
#include "TMVA/Config.h"
#include "TMVA/Configurable.h"
#include "TMVA/DataSet.h"
#include "TMVA/DataSetInfo.h"
#include "TMVA/Event.h"
#include "TMVA/IFitterTarget.h"
#include "TMVA/IMethod.h"
#include "TMVA/GeneticFitter.h"
#include "TMVA/Interval.h"
#include "TMVA/FitterBase.h"
#include "TMVA/MCFitter.h"
#include "TMVA/MethodBase.h"
#include "TMVA/MethodFDA.h"
#include "TMVA/MinuitFitter.h"
#include "TMVA/MsgLogger.h"
#include "TMVA/PDF.h"
#include "TMVA/Results.h"
#include "TMVA/SimulatedAnnealingFitter.h"
#include "TMVA/Timer.h"
#include "TMVA/Tools.h"
#include "TMVA/TransformationHandler.h"
#include "TMVA/Types.h"
#include "TMVA/TSpline1.h"
#include "TMVA/VariableTransformBase.h"

#include "TH1F.h"
#include "TObjString.h"
#include "TDirectory.h"
#include "TMath.h"
#include "TGraph.h"
#include "TSpline.h"
#include "TRandom3.h"

#include <cstdlib>
#include <iostream>
#include <iomanip>

using std::atof;

REGISTER_METHOD(Cuts)

ClassImp(TMVA::MethodCuts);

   const Double_t TMVA::MethodCuts::fgMaxAbsCutVal = 1.0e30;

////////////////////////////////////////////////////////////////////////////////
/// standard constructor

TMVA::MethodCuts::MethodCuts( const TString& jobName,
                              const TString& methodTitle,
                              DataSetInfo& theData,
                              const TString& theOption ) :
   MethodBase( jobName, Types::kCuts, methodTitle, theData, theOption),
   fFitMethod  ( kUseGeneticAlgorithm ),
   fEffMethod  ( kUseEventSelection ),
   fFitParams (0),
   fTestSignalEff(0.7),
   fEffSMin    ( 0 ),
   fEffSMax    ( 0 ),
   fCutRangeMin( 0 ),
   fCutRangeMax( 0 ),
   fBinaryTreeS( 0 ),
   fBinaryTreeB( 0 ),
   fCutMin     ( 0 ),
   fCutMax     ( 0 ),
   fTmpCutMin  ( 0 ),
   fTmpCutMax  ( 0 ),
   fAllVarsI   ( 0 ),
   fNpar       ( 0 ),
   fEffRef     ( 0 ),
   fRangeSign  ( 0 ),
   fRandom     ( 0 ),
   fMeanS      ( 0 ),
   fMeanB      ( 0 ),
   fRmsS       ( 0 ),
   fRmsB       ( 0 ),
   fEffBvsSLocal( 0 ),
   fVarHistS   ( 0 ),
   fVarHistB   ( 0 ),
   fVarHistS_smooth( 0 ),
   fVarHistB_smooth( 0 ),
   fVarPdfS    ( 0 ),
   fVarPdfB    ( 0 ),
   fNegEffWarning( kFALSE )
{
}

////////////////////////////////////////////////////////////////////////////////
/// construction from weight file

TMVA::MethodCuts::MethodCuts( DataSetInfo& theData,
                              const TString& theWeightFile) :
   MethodBase( Types::kCuts, theData, theWeightFile),
   fFitMethod  ( kUseGeneticAlgorithm ),
   fEffMethod  ( kUseEventSelection ),
   fFitParams (0),
   fTestSignalEff(0.7),
   fEffSMin    ( 0 ),
   fEffSMax    ( 0 ),
   fCutRangeMin( 0 ),
   fCutRangeMax( 0 ),
   fBinaryTreeS( 0 ),
   fBinaryTreeB( 0 ),
   fCutMin     ( 0 ),
   fCutMax     ( 0 ),
   fTmpCutMin  ( 0 ),
   fTmpCutMax  ( 0 ),
   fAllVarsI   ( 0 ),
   fNpar       ( 0 ),
   fEffRef     ( 0 ),
   fRangeSign  ( 0 ),
   fRandom     ( 0 ),
   fMeanS      ( 0 ),
   fMeanB      ( 0 ),
   fRmsS       ( 0 ),
   fRmsB       ( 0 ),
   fEffBvsSLocal( 0 ),
   fVarHistS   ( 0 ),
   fVarHistB   ( 0 ),
   fVarHistS_smooth( 0 ),
   fVarHistB_smooth( 0 ),
   fVarPdfS    ( 0 ),
   fVarPdfB    ( 0 ),
   fNegEffWarning( kFALSE )
{
}

////////////////////////////////////////////////////////////////////////////////
/// Cuts can only handle classification with 2 classes

Bool_t TMVA::MethodCuts::HasAnalysisType( Types::EAnalysisType type, UInt_t numberClasses,
                                          UInt_t /*numberTargets*/ )
{
   return (type == Types::kClassification && numberClasses == 2);
}

////////////////////////////////////////////////////////////////////////////////
/// default initialisation called by all constructors

void TMVA::MethodCuts::Init( void )
{
   fVarHistS          = fVarHistB = 0;
   fVarHistS_smooth   = fVarHistB_smooth = 0;
   fVarPdfS           = fVarPdfB = 0;
   fFitParams         = 0;
   fBinaryTreeS       = fBinaryTreeB = 0;
   fEffSMin           = 0;
   fEffSMax           = 0;

   // vector with fit results
   fNpar      = 2*GetNvar();
   fRangeSign = new std::vector<Int_t>   ( GetNvar() );
   for (UInt_t ivar=0; ivar<GetNvar(); ivar++) (*fRangeSign)[ivar] = +1;

   fMeanS     = new std::vector<Double_t>( GetNvar() );
   fMeanB     = new std::vector<Double_t>( GetNvar() );
   fRmsS      = new std::vector<Double_t>( GetNvar() );
   fRmsB      = new std::vector<Double_t>( GetNvar() );

   // get the variable specific options, first initialize default
   fFitParams = new std::vector<EFitParameters>( GetNvar() );
   for (UInt_t ivar=0; ivar<GetNvar(); ivar++) (*fFitParams)[ivar] = kNotEnforced;

   fFitMethod = kUseMonteCarlo;
   fTestSignalEff = -1;

   // create LUT for cuts
   fCutMin = new Double_t*[GetNvar()];
   fCutMax = new Double_t*[GetNvar()];
   for (UInt_t i=0; i<GetNvar(); i++) {
      fCutMin[i] = new Double_t[fNbins];
      fCutMax[i] = new Double_t[fNbins];
   }

   // init
   for (UInt_t ivar=0; ivar<GetNvar(); ivar++) {
      for (Int_t ibin=0; ibin<fNbins; ibin++) {
         fCutMin[ivar][ibin] = 0;
         fCutMax[ivar][ibin] = 0;
      }
   }

   fTmpCutMin = new Double_t[GetNvar()];
   fTmpCutMax = new Double_t[GetNvar()];
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TMVA::MethodCuts::~MethodCuts( void )
{
   delete fRangeSign;
   delete fMeanS;
   delete fMeanB;
   delete fRmsS;
   delete fRmsB;
   delete fFitParams;
   delete fEffBvsSLocal;

   if (NULL != fCutRangeMin) delete [] fCutRangeMin;
   if (NULL != fCutRangeMax) delete [] fCutRangeMax;
   if (NULL != fAllVarsI)    delete [] fAllVarsI;

   for (UInt_t i=0;i<GetNvar();i++) {
      if (NULL != fCutMin[i]  ) delete [] fCutMin[i];
      if (NULL != fCutMax[i]  ) delete [] fCutMax[i];
      if (NULL != fCutRange[i]) delete fCutRange[i];
   }

   if (NULL != fCutMin) delete [] fCutMin;
   if (NULL != fCutMax) delete [] fCutMax;

   if (NULL != fTmpCutMin) delete [] fTmpCutMin;
   if (NULL != fTmpCutMax) delete [] fTmpCutMax;

   if (NULL != fBinaryTreeS) delete fBinaryTreeS;
   if (NULL != fBinaryTreeB) delete fBinaryTreeB;
}

////////////////////////////////////////////////////////////////////////////////
/// define the options (their key words) that can be set in the option string.
///
/// know options:
///  - Method <string> Minimisation method. Available values are:
///    - MC Monte Carlo <default>
///    - GA Genetic Algorithm
///    - SA Simulated annealing
///
///  - EffMethod <string> Efficiency selection method. Available values are:
///    - EffSel <default>
///    - EffPDF
///
///  - VarProp <string> Property of variable 1 for the MC method (taking precedence over the
///    globale setting. The same values as for the global option are available. Variables 1..10 can be
///    set this way
///
///  - CutRangeMin/Max <float>  user-defined ranges in which cuts are varied

void TMVA::MethodCuts::DeclareOptions()
{
   DeclareOptionRef(fFitMethodS = "GA", "FitMethod", "Minimisation Method (GA, SA, and MC are the primary methods to be used; the others have been introduced for testing purposes and are depreciated)");
   AddPreDefVal(TString("GA"));
   AddPreDefVal(TString("SA"));
   AddPreDefVal(TString("MC"));
   AddPreDefVal(TString("MCEvents"));
   AddPreDefVal(TString("MINUIT"));
   AddPreDefVal(TString("EventScan"));

   // selection type
   DeclareOptionRef(fEffMethodS = "EffSel", "EffMethod", "Selection Method");
   AddPreDefVal(TString("EffSel"));
   AddPreDefVal(TString("EffPDF"));

   // cut ranges
   fCutRange.resize(GetNvar());
   fCutRangeMin = new Double_t[GetNvar()];
   fCutRangeMax = new Double_t[GetNvar()];
   for (UInt_t ivar=0; ivar<GetNvar(); ivar++) {
      fCutRange[ivar] = 0;
      fCutRangeMin[ivar] = fCutRangeMax[ivar] = -1;
   }

   DeclareOptionRef( fCutRangeMin, GetNvar(), "CutRangeMin", "Minimum of allowed cut range (set per variable)" );
   DeclareOptionRef( fCutRangeMax, GetNvar(), "CutRangeMax", "Maximum of allowed cut range (set per variable)" );

   fAllVarsI = new TString[GetNvar()];

   for (UInt_t i=0; i<GetNvar(); i++) fAllVarsI[i] = "NotEnforced";

   DeclareOptionRef(fAllVarsI, GetNvar(), "VarProp", "Categorisation of cuts");
   AddPreDefVal(TString("NotEnforced"));
   AddPreDefVal(TString("FMax"));
   AddPreDefVal(TString("FMin"));
   AddPreDefVal(TString("FSmart"));
}

////////////////////////////////////////////////////////////////////////////////
/// process user options.
///
/// sanity check, do not allow the input variables to be normalised, because this
/// only creates problems when interpreting the cuts

void TMVA::MethodCuts::ProcessOptions()
{
   if (IsNormalised()) {
      Log() << kWARNING << "Normalisation of the input variables for cut optimisation is not" << Endl;
      Log() << kWARNING << "supported because this provides intransparent cut values, and no" << Endl;
      Log() << kWARNING << "improvement in the performance of the algorithm." << Endl;
      Log() << kWARNING << "Please remove \"Normalise\" option from booking option string" << Endl;
      Log() << kWARNING << "==> Will reset normalisation flag to \"False\"" << Endl;
      SetNormalised( kFALSE );
   }

   if (IgnoreEventsWithNegWeightsInTraining()) {
      Log() << kFATAL << "Mechanism to ignore events with negative weights in training not yet available for method: "
            << GetMethodTypeName()
            << " --> Please remove \"IgnoreNegWeightsInTraining\" option from booking string."
            << Endl;
   }

   if      (fFitMethodS == "MC"      ) fFitMethod = kUseMonteCarlo;
   else if (fFitMethodS == "MCEvents") fFitMethod = kUseMonteCarloEvents;
   else if (fFitMethodS == "GA"      ) fFitMethod = kUseGeneticAlgorithm;
   else if (fFitMethodS == "SA"      ) fFitMethod = kUseSimulatedAnnealing;
   else if (fFitMethodS == "MINUIT"  ) {
      fFitMethod = kUseMinuit;
      Log() << kWARNING << "poor performance of MINUIT in MethodCuts; preferred fit method: GA" << Endl;
   }
   else if (fFitMethodS == "EventScan" ) fFitMethod = kUseEventScan;
   else Log() << kFATAL << "unknown minimisation method: " << fFitMethodS << Endl;

   if      (fEffMethodS == "EFFSEL" ) fEffMethod = kUseEventSelection; // highly recommended
   else if (fEffMethodS == "EFFPDF" ) fEffMethod = kUsePDFs;
   else                               fEffMethod = kUseEventSelection;

   // options output
   Log() << kINFO << Form("Use optimization method: \"%s\"",
                          (fFitMethod == kUseMonteCarlo) ? "Monte Carlo" :
                          (fFitMethod == kUseMonteCarlo) ? "Monte-Carlo-Event sampling" :
                          (fFitMethod == kUseEventScan)  ? "Full Event Scan (slow)" :
                          (fFitMethod == kUseMinuit)     ? "MINUIT" : "Genetic Algorithm" ) << Endl;
   Log() << kINFO << Form("Use efficiency computation method: \"%s\"",
                          (fEffMethod == kUseEventSelection) ? "Event Selection" : "PDF" ) << Endl;

   // cut ranges
   for (UInt_t ivar=0; ivar<GetNvar(); ivar++) {
      fCutRange[ivar] = new Interval( fCutRangeMin[ivar], fCutRangeMax[ivar] );
   }

   // individual options
   for (UInt_t ivar=0; ivar<GetNvar(); ivar++) {
      EFitParameters theFitP = kNotEnforced;
      if (fAllVarsI[ivar] == "" || fAllVarsI[ivar] == "NotEnforced") theFitP = kNotEnforced;
      else if (fAllVarsI[ivar] == "FMax" )                           theFitP = kForceMax;
      else if (fAllVarsI[ivar] == "FMin" )                           theFitP = kForceMin;
      else if (fAllVarsI[ivar] == "FSmart" )                         theFitP = kForceSmart;
      else {
         Log() << kFATAL << "unknown value \'" << fAllVarsI[ivar]
               << "\' for fit parameter option " << Form("VarProp[%i]",ivar) << Endl;
      }
      (*fFitParams)[ivar] = theFitP;

      if (theFitP != kNotEnforced)
         Log() << kINFO << "Use \"" << fAllVarsI[ivar]
               << "\" cuts for variable: " << "'" << (*fInputVars)[ivar] << "'" << Endl;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// cut evaluation: returns 1.0 if event passed, 0.0 otherwise

Double_t TMVA::MethodCuts::GetMvaValue( Double_t* err, Double_t* errUpper )
{
   // cannot determine error
   NoErrorCalc(err, errUpper);

   // sanity check
   if (fCutMin == NULL || fCutMax == NULL || fNbins == 0) {
      Log() << kFATAL << "<Eval_Cuts> fCutMin/Max have zero pointer. "
            << "Did you book Cuts ?" << Endl;
   }

   const Event* ev = GetEvent();

   // sanity check
   if (fTestSignalEff > 0) {
      // get efficiency bin
      Int_t ibin = fEffBvsSLocal->FindBin( fTestSignalEff );
      if (ibin < 0      ) ibin = 0;
      else if (ibin >= fNbins) ibin = fNbins - 1;

      Bool_t passed = kTRUE;
      for (UInt_t ivar=0; ivar<GetNvar(); ivar++)
         passed &= ( (ev->GetValue(ivar) >  fCutMin[ivar][ibin]) &&
                     (ev->GetValue(ivar) <= fCutMax[ivar][ibin]) );

      return passed ? 1. : 0. ;
   }
   else return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// print cuts

void TMVA::MethodCuts::PrintCuts( Double_t effS ) const
{
   std::vector<Double_t> cutsMin;
   std::vector<Double_t> cutsMax;
   Int_t ibin = fEffBvsSLocal->FindBin( effS );

   Double_t trueEffS = GetCuts( effS, cutsMin, cutsMax );

   // retrieve variable expressions (could be transformations)
   std::vector<TString>* varVec = 0;
   if (GetTransformationHandler().GetNumOfTransformations() == 0) {
      // no transformation applied, replace by current variables
      varVec = new std::vector<TString>;
      for (UInt_t ivar=0; ivar<cutsMin.size(); ivar++) {
         varVec->push_back( DataInfo().GetVariableInfo(ivar).GetLabel() );
      }
   }
   else if (GetTransformationHandler().GetNumOfTransformations() == 1) {
      // get transformation string
      varVec = GetTransformationHandler().GetTransformationStringsOfLastTransform();
   }
   else {
      // replace transformation print by current variables and indicated incompleteness
      varVec = new std::vector<TString>;
      for (UInt_t ivar=0; ivar<cutsMin.size(); ivar++) {
         varVec->push_back( DataInfo().GetVariableInfo(ivar).GetLabel() + " [transformed]" );
      }
   }

   UInt_t maxL = 0;
   for (UInt_t ivar=0; ivar<cutsMin.size(); ivar++) {
      if ((UInt_t)(*varVec)[ivar].Length() > maxL) maxL = (*varVec)[ivar].Length();
   }
   UInt_t maxLine = 20+maxL+16;

   for (UInt_t i=0; i<maxLine; i++) Log() << "-";
   Log() << Endl;
   Log() << kHEADER << "Cut values for requested signal efficiency: " << trueEffS << Endl;
   Log() << kINFO << "Corresponding background efficiency       : " << fEffBvsSLocal->GetBinContent( ibin ) << Endl;
   if (GetTransformationHandler().GetNumOfTransformations() == 1) {
      Log() << kINFO << "Transformation applied to input variables : \""
            << GetTransformationHandler().GetNameOfLastTransform() << "\"" << Endl;
   }
   else if (GetTransformationHandler().GetNumOfTransformations() > 1) {
      Log() << kINFO << "[ More than one (=" << GetTransformationHandler().GetNumOfTransformations() << ") "
            << " transformations applied in transformation chain; cuts applied on transformed quantities ] " << Endl;
   }
   else {
      Log() << kINFO << "Transformation applied to input variables : None"  << Endl;
   }
   for (UInt_t i=0; i<maxLine; i++) Log() << "-";
   Log() << Endl;
   for (UInt_t ivar=0; ivar<cutsMin.size(); ivar++) {
      Log() << kINFO
            << "Cut[" << std::setw(2) << ivar << "]: "
            << std::setw(10) << cutsMin[ivar]
            << " < "
            << std::setw(maxL) << (*varVec)[ivar]
            << " <= "
            << std::setw(10) << cutsMax[ivar] << Endl;
   }
   for (UInt_t i=0; i<maxLine; i++) Log() << "-";
   Log() << Endl;

   delete varVec; // yes, ownership has been given to us
}

////////////////////////////////////////////////////////////////////////////////
/// retrieve cut values for given signal efficiency
/// assume vector of correct size !!

Double_t TMVA::MethodCuts::GetCuts( Double_t effS, Double_t* cutMin, Double_t* cutMax ) const
{
   std::vector<Double_t> cMin( GetNvar() );
   std::vector<Double_t> cMax( GetNvar() );
   Double_t trueEffS = GetCuts( effS, cMin, cMax );
   for (UInt_t ivar=0; ivar<GetNvar(); ivar++) {
      cutMin[ivar] = cMin[ivar];
      cutMax[ivar] = cMax[ivar];
   }
   return trueEffS;
}

////////////////////////////////////////////////////////////////////////////////
/// retrieve cut values for given signal efficiency

Double_t TMVA::MethodCuts::GetCuts( Double_t effS,
                                    std::vector<Double_t>& cutMin,
                                    std::vector<Double_t>& cutMax ) const
{
   // find corresponding bin
   Int_t ibin = fEffBvsSLocal->FindBin( effS );

   // get the true efficiency which is the one on the "left hand" side of the bin
   Double_t trueEffS = fEffBvsSLocal->GetBinLowEdge( ibin );

   ibin--; // the 'cut' vector has 0...fNbins indices
   if      (ibin < 0      ) ibin = 0;
   else if (ibin >= fNbins) ibin = fNbins - 1;

   cutMin.clear();
   cutMax.clear();
   for (UInt_t ivar=0; ivar<GetNvar(); ivar++) {
      cutMin.push_back( fCutMin[ivar][ibin]  );
      cutMax.push_back( fCutMax[ivar][ibin] );
   }

   return trueEffS;
}

////////////////////////////////////////////////////////////////////////////////
/// training method: here the cuts are optimised for the training sample

void  TMVA::MethodCuts::Train( void )
{
   if (fEffMethod == kUsePDFs) CreateVariablePDFs(); // create PDFs for variables

   // create binary trees (global member variables) for signal and background
   if (fBinaryTreeS != 0) { delete fBinaryTreeS; fBinaryTreeS = 0; }
   if (fBinaryTreeB != 0) { delete fBinaryTreeB; fBinaryTreeB = 0; }

   // the variables may be transformed by a transformation method: to coherently
   // treat signal and background one must decide which transformation type shall
   // be used: our default is signal-type

   fBinaryTreeS = new BinarySearchTree();
   fBinaryTreeS->Fill( GetEventCollection(Types::kTraining), fSignalClass );
   fBinaryTreeB = new BinarySearchTree();
   fBinaryTreeB->Fill( GetEventCollection(Types::kTraining), fBackgroundClass );

   for (UInt_t ivar =0; ivar < Data()->GetNVariables(); ivar++) {
      (*fMeanS)[ivar] = fBinaryTreeS->Mean(Types::kSignal, ivar);
      (*fRmsS)[ivar]  = fBinaryTreeS->RMS (Types::kSignal, ivar);
      (*fMeanB)[ivar] = fBinaryTreeB->Mean(Types::kBackground, ivar);
      (*fRmsB)[ivar]  = fBinaryTreeB->RMS (Types::kBackground, ivar);

      // update interval ?
      Double_t xmin = TMath::Min(fBinaryTreeS->Min(Types::kSignal,     ivar),
                                 fBinaryTreeB->Min(Types::kBackground, ivar));
      Double_t xmax = TMath::Max(fBinaryTreeS->Max(Types::kSignal,     ivar),
                                 fBinaryTreeB->Max(Types::kBackground, ivar));

      // redefine ranges to be slightly smaller and larger than xmin and xmax, respectively
      Double_t eps = 0.01*(xmax - xmin);
      xmin -= eps;
      xmax += eps;

      if (TMath::Abs(fCutRange[ivar]->GetMin() - fCutRange[ivar]->GetMax()) < 1.0e-300 ) {
         fCutRange[ivar]->SetMin( xmin );
         fCutRange[ivar]->SetMax( xmax );
      }
      else if (xmin > fCutRange[ivar]->GetMin()) fCutRange[ivar]->SetMin( xmin );
      else if (xmax < fCutRange[ivar]->GetMax()) fCutRange[ivar]->SetMax( xmax );
   }

   std::vector<TH1F*> signalDist, bkgDist;

   // this is important: reset the branch addresses of the training tree to the current event
   delete fEffBvsSLocal;
   fEffBvsSLocal = new TH1F( GetTestvarName() + "_effBvsSLocal",
                             TString(GetName()) + " efficiency of B vs S", fNbins, 0.0, 1.0 );
   fEffBvsSLocal->SetDirectory(0); // it's local

   // init
   for (Int_t ibin=1; ibin<=fNbins; ibin++) fEffBvsSLocal->SetBinContent( ibin, -0.1 );

   // --------------------------------------------------------------------------
   if (fFitMethod == kUseGeneticAlgorithm ||
       fFitMethod == kUseMonteCarlo       ||
       fFitMethod == kUseMinuit           ||
       fFitMethod == kUseSimulatedAnnealing) {

      // ranges
      std::vector<Interval*> ranges;

      for (UInt_t ivar=0; ivar<GetNvar(); ivar++) {

         Int_t nbins = 0;
         if (DataInfo().GetVariableInfo(ivar).GetVarType() == 'I') {
            nbins = Int_t(fCutRange[ivar]->GetMax() - fCutRange[ivar]->GetMin()) + 1;
         }

         if ((*fFitParams)[ivar] == kForceSmart) {
            if ((*fMeanS)[ivar] > (*fMeanB)[ivar]) (*fFitParams)[ivar] = kForceMax;
            else                                   (*fFitParams)[ivar] = kForceMin;
         }

         if ((*fFitParams)[ivar] == kForceMin) {
            ranges.push_back( new Interval( fCutRange[ivar]->GetMin(), fCutRange[ivar]->GetMin(), nbins ) );
            ranges.push_back( new Interval( 0, fCutRange[ivar]->GetMax() - fCutRange[ivar]->GetMin(), nbins ) );
         }
         else if ((*fFitParams)[ivar] == kForceMax) {
            ranges.push_back( new Interval( fCutRange[ivar]->GetMin(), fCutRange[ivar]->GetMax(), nbins ) );
            ranges.push_back( new Interval( fCutRange[ivar]->GetMax() - fCutRange[ivar]->GetMin(),
                                            fCutRange[ivar]->GetMax() - fCutRange[ivar]->GetMin(), nbins ) );
         }
         else {
            ranges.push_back( new Interval( fCutRange[ivar]->GetMin(), fCutRange[ivar]->GetMax(), nbins ) );
            ranges.push_back( new Interval( 0, fCutRange[ivar]->GetMax() - fCutRange[ivar]->GetMin(), nbins ) );
         }
      }

      // create the fitter
      FitterBase* fitter = NULL;

      switch (fFitMethod) {
      case kUseGeneticAlgorithm:
         fitter = new GeneticFitter( *this, Form("%sFitter_GA",     GetName()), ranges, GetOptions() );
         break;
      case kUseMonteCarlo:
         fitter = new MCFitter     ( *this, Form("%sFitter_MC",     GetName()), ranges, GetOptions() );
         break;
      case kUseMinuit:
         fitter = new MinuitFitter ( *this, Form("%sFitter_MINUIT", GetName()), ranges, GetOptions() );
         break;
      case kUseSimulatedAnnealing:
         fitter = new SimulatedAnnealingFitter( *this, Form("%sFitter_SA", GetName()), ranges, GetOptions() );
         break;
      default:
         Log() << kFATAL << "Wrong fit method: " << fFitMethod << Endl;
      }

      if (fInteractive) fitter->SetIPythonInteractive(&fExitFromTraining, &fIPyMaxIter, &fIPyCurrentIter);

      fitter->CheckForUnusedOptions();

      // perform the fit
      fitter->Run();

      // clean up
      for (UInt_t ivar=0; ivar<ranges.size(); ivar++) delete ranges[ivar];
      delete fitter;

   }
   // --------------------------------------------------------------------------
   else if (fFitMethod == kUseEventScan) {

      Int_t nevents = Data()->GetNEvents();
      Int_t ic = 0;

      // timing of MC
      Int_t nsamples = Int_t(0.5*nevents*(nevents - 1));
      Timer timer( nsamples, GetName() );
      fIPyMaxIter = nsamples;

      Log() << kINFO << "Running full event scan: " << Endl;
      for (Int_t ievt1=0; ievt1<nevents; ievt1++) {
         for (Int_t ievt2=ievt1+1; ievt2<nevents; ievt2++) {

           fIPyCurrentIter = ic;
           if (fExitFromTraining) break;
            EstimatorFunction( ievt1, ievt2 );

            // what's the time please?
            ic++;
            if ((nsamples<10000) || ic%10000 == 0) timer.DrawProgressBar( ic );
         }
      }
   }
   // --------------------------------------------------------------------------
   else if (fFitMethod == kUseMonteCarloEvents) {

      Int_t  nsamples = 200000;
      UInt_t seed     = 100;
      DeclareOptionRef( nsamples, "SampleSize", "Number of Monte-Carlo-Event samples" );
      DeclareOptionRef( seed,     "Seed",       "Seed for the random generator (0 takes random seeds)" );
      ParseOptions();

      Int_t nevents = Data()->GetNEvents();
      Int_t ic = 0;

      // timing of MC
      Timer timer( nsamples, GetName() );
      fIPyMaxIter = nsamples;

      // random generator
      TRandom3*rnd = new TRandom3( seed );

      Log() << kINFO << "Running Monte-Carlo-Event sampling over " << nsamples << " events" << Endl;
      std::vector<Double_t> pars( 2*GetNvar() );

      for (Int_t itoy=0; itoy<nsamples; itoy++) {
        fIPyCurrentIter = ic;
        if (fExitFromTraining) break;

         for (UInt_t ivar=0; ivar<GetNvar(); ivar++) {

            // generate minimum and delta cuts for this variable

            // retrieve signal events
            Bool_t isSignal = kFALSE;
            Int_t    ievt1, ievt2;
            Double_t evt1 = 0., evt2 = 0.;
            Int_t nbreak = 0;
            while (!isSignal) {
               ievt1 = Int_t(rnd->Uniform(0.,1.)*nevents);
               ievt2 = Int_t(rnd->Uniform(0.,1.)*nevents);

               const Event *ev1 = GetEvent(ievt1);
               isSignal = DataInfo().IsSignal(ev1);
               evt1 = ev1->GetValue( ivar );

               const Event *ev2 = GetEvent(ievt2);
               isSignal &= DataInfo().IsSignal(ev2);
               evt2 = ev2->GetValue( ivar );

               if (nbreak++ > 10000) {
                  Log() << kFATAL << "<MCEvents>: could not find signal events"
                                           << " after 10000 trials - do you have signal events in your sample ?"
                                           << Endl;
                  isSignal = 1;
               }
            }

            // sort
            if (evt1 > evt2) { Double_t z = evt1; evt1 = evt2; evt2 = z; }
            pars[2*ivar]   = evt1;
            pars[2*ivar+1] = evt2 - evt1;
         }

         // compute estimator
         EstimatorFunction( pars );

         // what's the time please?
         ic++;
         if ((nsamples<1000) || ic%1000 == 0) timer.DrawProgressBar( ic );
      }

      delete rnd;
   }
   // --------------------------------------------------------------------------
   else Log() << kFATAL << "Unknown minimisation method: " << fFitMethod << Endl;

   if (fBinaryTreeS != 0) { delete fBinaryTreeS; fBinaryTreeS = 0; }
   if (fBinaryTreeB != 0) { delete fBinaryTreeB; fBinaryTreeB = 0; }

   // force cut ranges within limits
   for (UInt_t ivar=0; ivar<GetNvar(); ivar++) {
      for (Int_t ibin=0; ibin<fNbins; ibin++) {

         if ((*fFitParams)[ivar] == kForceMin && fCutMin[ivar][ibin] > -fgMaxAbsCutVal) {
            fCutMin[ivar][ibin] = -fgMaxAbsCutVal;
         }
         if ((*fFitParams)[ivar] == kForceMax && fCutMax[ivar][ibin] < fgMaxAbsCutVal) {
            fCutMax[ivar][ibin] = fgMaxAbsCutVal;
         }
      }
   }

   // some output
   // the efficiency which is asked for has to be slightly higher than the bin-borders.
   // if not, then the wrong bin is taken in some cases.
   Double_t epsilon = 0.0001;
   for (Double_t eff=0.1; eff<0.95; eff += 0.1) PrintCuts( eff+epsilon );

   if (!fExitFromTraining) fIPyMaxIter = fIPyCurrentIter;
   ExitFromTraining();
}

////////////////////////////////////////////////////////////////////////////////
/// nothing to test

void TMVA::MethodCuts::TestClassification()
{
}

////////////////////////////////////////////////////////////////////////////////
/// for full event scan

Double_t TMVA::MethodCuts::EstimatorFunction( Int_t ievt1, Int_t ievt2 )
{
   const Event *ev1 = GetEvent(ievt1);
   if (!DataInfo().IsSignal(ev1)) return -1;

   const Event *ev2 = GetEvent(ievt2);
   if (!DataInfo().IsSignal(ev2)) return -1;

   const Int_t nvar = GetNvar();
   Double_t* evt1 = new Double_t[nvar];
   Double_t* evt2 = new Double_t[nvar];

   for (Int_t ivar=0; ivar<nvar; ivar++) {
      evt1[ivar] = ev1->GetValue( ivar );
      evt2[ivar] = ev2->GetValue( ivar );
   }

   // determine cuts
   std::vector<Double_t> pars;
   for (Int_t ivar=0; ivar<nvar; ivar++) {
      Double_t cutMin;
      Double_t cutMax;
      if (evt1[ivar] < evt2[ivar]) {
         cutMin = evt1[ivar];
         cutMax = evt2[ivar];
      }
      else {
         cutMin = evt2[ivar];
         cutMax = evt1[ivar];
      }

      pars.push_back( cutMin );
      pars.push_back( cutMax - cutMin );
   }

   delete [] evt1;
   delete [] evt2;

   return ComputeEstimator( pars );
}

////////////////////////////////////////////////////////////////////////////////
/// returns estimator for "cut fitness" used by GA

Double_t TMVA::MethodCuts::EstimatorFunction( std::vector<Double_t>& pars )
{
   return ComputeEstimator( pars );
}

////////////////////////////////////////////////////////////////////////////////
/// returns estimator for "cut fitness" used by GA.
///
/// there are two requirements:
///  1. the signal efficiency must be equal to the required one in the
///    efficiency scan
///  2. the background efficiency must be as small as possible
///
/// the requirement 1. has priority over 2.

Double_t TMVA::MethodCuts::ComputeEstimator( std::vector<Double_t>& pars )
{
   // caution: the npar gives the _free_ parameters
   // however: the "pars" array contains all parameters

   // determine cuts
   Double_t effS = 0, effB = 0;
   this->MatchParsToCuts( pars, &fTmpCutMin[0], &fTmpCutMax[0] );

   // retrieve signal and background efficiencies for given cut
   switch (fEffMethod) {
   case kUsePDFs:
      this->GetEffsfromPDFs      (&fTmpCutMin[0], &fTmpCutMax[0], effS, effB);
      break;
   case kUseEventSelection:
      this->GetEffsfromSelection (&fTmpCutMin[0], &fTmpCutMax[0], effS, effB);
      break;
   default:
      this->GetEffsfromSelection (&fTmpCutMin[0], &fTmpCutMax[0], effS, effB);
   }

   Double_t eta = 0;

   // test for a estimator function which optimizes on the whole background-rejection signal-efficiency plot

   // get the backg-reject. and sig-eff for the parameters given to this function
   // effS, effB

   // get best background rejection for given signal efficiency
   Int_t ibinS = fEffBvsSLocal->FindBin( effS );

   Double_t effBH       = fEffBvsSLocal->GetBinContent( ibinS );
   Double_t effBH_left  = (ibinS > 1     ) ? fEffBvsSLocal->GetBinContent( ibinS-1 ) : effBH;
   Double_t effBH_right = (ibinS < fNbins) ? fEffBvsSLocal->GetBinContent( ibinS+1 ) : effBH;

   Double_t average = 0.5*(effBH_left + effBH_right);
   if (effBH < effB) average = effBH;

   // if the average of the bin right and left is larger than this one, add the difference to
   // the current value of the estimator (because you can do at least so much better)
   eta = ( -TMath::Abs(effBH-average) + (1.0 - (effBH - effB))) / (1.0 + effS);
   // alternative idea
   //if (effBH<0) eta = (1.e-6+effB)/(1.0 + effS);
   //else eta =  (effB - effBH) * (1.0 + 10.* effS);

   // if a point is found which is better than an existing one, ... replace it.
   // preliminary best event -> backup
   if (effBH < 0 || effBH > effB) {
      fEffBvsSLocal->SetBinContent( ibinS, effB );
      for (UInt_t ivar=0; ivar<GetNvar(); ivar++) {
         fCutMin[ivar][ibinS-1] = fTmpCutMin[ivar]; // bin 1 stored in index 0
         fCutMax[ivar][ibinS-1] = fTmpCutMax[ivar];
      }
   }

   // caution (!) this value is not good for a decision for MC, .. it is designed for GA
   // but .. it doesn't matter, as MC samplings are independent from the former ones
   // and the replacement of the best variables by better ones is done about 10 lines above.
   // ( if (effBH < 0 || effBH > effB) { .... )

   if (ibinS<=1) {
      // add penalty for effS=0 bin
      // to avoid that the minimizer gets stuck in the zero-bin
      // force it towards higher efficiency
      Double_t penalty=0.,diff=0.;
      for (UInt_t ivar=0; ivar<GetNvar(); ivar++) {
         diff=(fCutRange[ivar]->GetMax()-fTmpCutMax[ivar])/(fCutRange[ivar]->GetMax()-fCutRange[ivar]->GetMin());
         penalty+=diff*diff;
         diff=(fCutRange[ivar]->GetMin()-fTmpCutMin[ivar])/(fCutRange[ivar]->GetMax()-fCutRange[ivar]->GetMin());
         penalty+=4.*diff*diff;
      }

      if (effS<1.e-4) return 10.0+penalty;
      else return 10.*(1.-10.*effS);
   }
   return eta;
}

////////////////////////////////////////////////////////////////////////////////
/// translates parameters into cuts

void TMVA::MethodCuts::MatchParsToCuts( const std::vector<Double_t> & pars,
                                        Double_t* cutMin, Double_t* cutMax )
{
   for (UInt_t ivar=0; ivar<GetNvar(); ivar++) {
      Int_t ipar = 2*ivar;
      cutMin[ivar] = ((*fRangeSign)[ivar] > 0) ? pars[ipar] : pars[ipar] - pars[ipar+1];
      cutMax[ivar] = ((*fRangeSign)[ivar] > 0) ? pars[ipar] + pars[ipar+1] : pars[ipar];
   }
}

////////////////////////////////////////////////////////////////////////////////
/// translate the cuts into parameters (obsolete function)

void TMVA::MethodCuts::MatchCutsToPars( std::vector<Double_t>& pars,
                                        Double_t** cutMinAll, Double_t** cutMaxAll, Int_t ibin )
{
   if (ibin < 1 || ibin > fNbins) Log() << kFATAL << "::MatchCutsToPars: bin error: "
                                        << ibin << Endl;

   const UInt_t nvar = GetNvar();
   Double_t *cutMin = new Double_t[nvar];
   Double_t *cutMax = new Double_t[nvar];
   for (UInt_t ivar=0; ivar<nvar; ivar++) {
      cutMin[ivar] = cutMinAll[ivar][ibin-1];
      cutMax[ivar] = cutMaxAll[ivar][ibin-1];
   }

   MatchCutsToPars( pars, cutMin, cutMax );
   delete [] cutMin;
   delete [] cutMax;
}

////////////////////////////////////////////////////////////////////////////////
/// translates cuts into parameters

void TMVA::MethodCuts::MatchCutsToPars( std::vector<Double_t>& pars,
                                        Double_t* cutMin, Double_t* cutMax )
{
   for (UInt_t ivar=0; ivar<GetNvar(); ivar++) {
      Int_t ipar = 2*ivar;
      pars[ipar]   = ((*fRangeSign)[ivar] > 0) ? cutMin[ivar] : cutMax[ivar];
      pars[ipar+1] = cutMax[ivar] - cutMin[ivar];
   }
}

////////////////////////////////////////////////////////////////////////////////
/// compute signal and background efficiencies from PDFs
/// for given cut sample

void TMVA::MethodCuts::GetEffsfromPDFs( Double_t* cutMin, Double_t* cutMax,
                                        Double_t& effS, Double_t& effB )
{
   effS = 1.0;
   effB = 1.0;
   for (UInt_t ivar=0; ivar<GetNvar(); ivar++) {
      effS *= (*fVarPdfS)[ivar]->GetIntegral( cutMin[ivar], cutMax[ivar] );
      effB *= (*fVarPdfB)[ivar]->GetIntegral( cutMin[ivar], cutMax[ivar] );
   }

   // quick fix to prevent from efficiencies < 0
   if( effS < 0.0 ) {
      effS = 0.0;
      if( !fNegEffWarning ) Log() << kWARNING << "Negative signal efficiency found and set to 0. This is probably due to many events with negative weights in a certain cut-region." << Endl;
      fNegEffWarning = kTRUE;
   }
   if( effB < 0.0 ) {
      effB = 0.0;
      if( !fNegEffWarning ) Log() << kWARNING << "Negative background efficiency found and set to 0. This is probably due to many events with negative weights in a certain cut-region." << Endl;
      fNegEffWarning = kTRUE;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// compute signal and background efficiencies from event counting
/// for given cut sample

void TMVA::MethodCuts::GetEffsfromSelection( Double_t* cutMin, Double_t* cutMax,
                                             Double_t& effS, Double_t& effB)
{
   Float_t nTotS = 0, nTotB = 0;
   Float_t nSelS = 0, nSelB = 0;

   Volume* volume = new Volume( cutMin, cutMax, GetNvar() );

   // search for all events lying in the volume, and add up their weights
   nSelS = fBinaryTreeS->SearchVolume( volume );
   nSelB = fBinaryTreeB->SearchVolume( volume );

   delete volume;

   // total number of "events" (sum of weights) as reference to compute efficiency
   nTotS = fBinaryTreeS->GetSumOfWeights();
   nTotB = fBinaryTreeB->GetSumOfWeights();

   // sanity check
   if (nTotS == 0 && nTotB == 0) {
      Log() << kFATAL << "<GetEffsfromSelection> fatal error in zero total number of events:"
            << " nTotS, nTotB: " << nTotS << " " << nTotB << " ***" << Endl;
   }

   // efficiencies
   if (nTotS == 0 ) {
      effS = 0;
      effB = nSelB/nTotB;
      Log() << kWARNING << "<ComputeEstimator> zero number of signal events" << Endl;
   }
   else if (nTotB == 0) {
      effB = 0;
      effS = nSelS/nTotS;
      Log() << kWARNING << "<ComputeEstimator> zero number of background events" << Endl;
   }
   else {
      effS = nSelS/nTotS;
      effB = nSelB/nTotB;
   }

   // quick fix to prevent from efficiencies < 0
   if( effS < 0.0 ) {
      effS = 0.0;
      if( !fNegEffWarning ) Log() << kWARNING << "Negative signal efficiency found and set to 0. This is probably due to many events with negative weights in a certain cut-region." << Endl;
      fNegEffWarning = kTRUE;
   }
   if( effB < 0.0 ) {
      effB = 0.0;
      if( !fNegEffWarning ) Log() << kWARNING << "Negative background efficiency found and set to 0. This is probably due to many events with negative weights in a certain cut-region." << Endl;
      fNegEffWarning = kTRUE;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// for PDF method: create efficiency reference histograms and PDFs

void TMVA::MethodCuts::CreateVariablePDFs( void )
{
   // create list of histograms and PDFs
   fVarHistS        = new std::vector<TH1*>( GetNvar() );
   fVarHistB        = new std::vector<TH1*>( GetNvar() );
   fVarHistS_smooth = new std::vector<TH1*>( GetNvar() );
   fVarHistB_smooth = new std::vector<TH1*>( GetNvar() );
   fVarPdfS         = new std::vector<PDF*>( GetNvar() );
   fVarPdfB         = new std::vector<PDF*>( GetNvar() );

   Int_t nsmooth = 0;

   // get min and max values of all events
   Double_t minVal = DBL_MAX;
   Double_t maxVal = -DBL_MAX;
   for( UInt_t ievt=0; ievt<Data()->GetNEvents(); ievt++ ){
      const Event *ev = GetEvent(ievt);
      Float_t val = ev->GetValue(ievt);
      if( val > minVal ) minVal = val;
      if( val < maxVal ) maxVal = val;
   }

   for (UInt_t ivar=0; ivar<GetNvar(); ivar++) {

      // ---- signal
      TString histTitle = (*fInputVars)[ivar] + " signal training";
      TString histName  = (*fInputVars)[ivar] + "_sig";
      //      TString drawOpt   = (*fInputVars)[ivar] + ">>h(";
      //      drawOpt += fNbins;
      //      drawOpt += ")";

      // selection
      //      Data().GetTrainingTree()->Draw( drawOpt, "type==1", "goff" );
      //      (*fVarHistS)[ivar] = (TH1F*)gDirectory->Get("h");
      //      (*fVarHistS)[ivar]->SetName(histName);
      //      (*fVarHistS)[ivar]->SetTitle(histTitle);

      (*fVarHistS)[ivar] = new TH1F(histName.Data(), histTitle.Data(), fNbins, minVal, maxVal );

      // ---- background
      histTitle = (*fInputVars)[ivar] + " background training";
      histName  = (*fInputVars)[ivar] + "_bgd";
      //      drawOpt   = (*fInputVars)[ivar] + ">>h(";
      //      drawOpt += fNbins;
      //      drawOpt += ")";

      //      Data().GetTrainingTree()->Draw( drawOpt, "type==0", "goff" );
      //      (*fVarHistB)[ivar] = (TH1F*)gDirectory->Get("h");
      //      (*fVarHistB)[ivar]->SetName(histName);
      //      (*fVarHistB)[ivar]->SetTitle(histTitle);


      (*fVarHistB)[ivar] = new TH1F(histName.Data(), histTitle.Data(), fNbins, minVal, maxVal );

      for( UInt_t ievt=0; ievt<Data()->GetNEvents(); ievt++ ){
         const Event *ev = GetEvent(ievt);
         Float_t val = ev->GetValue(ievt);
         if( DataInfo().IsSignal(ev) ){
            (*fVarHistS)[ivar]->Fill( val );
         }else{
            (*fVarHistB)[ivar]->Fill( val );
         }
      }



      // make copy for smoothed histos
      (*fVarHistS_smooth)[ivar] = (TH1F*)(*fVarHistS)[ivar]->Clone();
      histTitle =  (*fInputVars)[ivar] + " signal training  smoothed ";
      histTitle += nsmooth;
      histTitle +=" times";
      histName =  (*fInputVars)[ivar] + "_sig_smooth";
      (*fVarHistS_smooth)[ivar]->SetName(histName);
      (*fVarHistS_smooth)[ivar]->SetTitle(histTitle);

      // smooth
      (*fVarHistS_smooth)[ivar]->Smooth(nsmooth);

      // ---- background
      //      histTitle = (*fInputVars)[ivar] + " background training";
      //      histName  = (*fInputVars)[ivar] + "_bgd";
      //      drawOpt   = (*fInputVars)[ivar] + ">>h(";
      //      drawOpt += fNbins;
      //      drawOpt += ")";

      //      Data().GetTrainingTree()->Draw( drawOpt, "type==0", "goff" );
      //      (*fVarHistB)[ivar] = (TH1F*)gDirectory->Get("h");
      //      (*fVarHistB)[ivar]->SetName(histName);
      //      (*fVarHistB)[ivar]->SetTitle(histTitle);

      // make copy for smoothed histos
      (*fVarHistB_smooth)[ivar] = (TH1F*)(*fVarHistB)[ivar]->Clone();
      histTitle  = (*fInputVars)[ivar]+" background training  smoothed ";
      histTitle += nsmooth;
      histTitle +=" times";
      histName   = (*fInputVars)[ivar]+"_bgd_smooth";
      (*fVarHistB_smooth)[ivar]->SetName(histName);
      (*fVarHistB_smooth)[ivar]->SetTitle(histTitle);

      // smooth
      (*fVarHistB_smooth)[ivar]->Smooth(nsmooth);

      // create PDFs
      (*fVarPdfS)[ivar] = new PDF( TString(GetName()) + " PDF Var Sig " + GetInputVar( ivar ), (*fVarHistS_smooth)[ivar], PDF::kSpline2 );
      (*fVarPdfB)[ivar] = new PDF( TString(GetName()) + " PDF Var Bkg " + GetInputVar( ivar ), (*fVarHistB_smooth)[ivar], PDF::kSpline2 );
   }
}

////////////////////////////////////////////////////////////////////////////////
/// read the cuts from stream

void  TMVA::MethodCuts::ReadWeightsFromStream( std::istream& istr )
{
   TString dummy;
   UInt_t  dummyInt;

   // first the dimensions
   istr >> dummy >> dummy;
   // coverity[tainted_data_argument]
   istr >> dummy >> fNbins;

   // get rid of one read-in here because we read in once all ready to check for decorrelation
   istr >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummyInt >> dummy ;

   // sanity check
   if (dummyInt != Data()->GetNVariables()) {
      Log() << kFATAL << "<ReadWeightsFromStream> fatal error: mismatch "
            << "in number of variables: " << dummyInt << " != " << Data()->GetNVariables() << Endl;
   }
   //SetNvar(dummyInt);

   // print some information
   if (fFitMethod == kUseMonteCarlo) {
      Log() << kWARNING << "Read cuts optimised using sample of MC events" << Endl;
   }
   else if (fFitMethod == kUseMonteCarloEvents) {
      Log() << kWARNING << "Read cuts optimised using sample of MC events" << Endl;
   }
   else if (fFitMethod == kUseGeneticAlgorithm) {
      Log() << kINFO << "Read cuts optimised using Genetic Algorithm" << Endl;
   }
   else if (fFitMethod == kUseSimulatedAnnealing) {
      Log() << kINFO << "Read cuts optimised using Simulated Annealing algorithm" << Endl;
   }
   else if (fFitMethod == kUseEventScan) {
      Log() << kINFO << "Read cuts optimised using Full Event Scan" << Endl;
   }
   else {
      Log() << kWARNING << "unknown method: " << fFitMethod << Endl;
   }
   Log() << kINFO << "in " << fNbins << " signal efficiency bins and for " << GetNvar() << " variables" << Endl;

   // now read the cuts
   char buffer[200];
   istr.getline(buffer,200);
   istr.getline(buffer,200);

   Int_t   tmpbin;
   Float_t tmpeffS, tmpeffB;
   if (fEffBvsSLocal != 0) delete fEffBvsSLocal;
   fEffBvsSLocal = new TH1F( GetTestvarName() + "_effBvsSLocal",
                             TString(GetName()) + " efficiency of B vs S", fNbins, 0.0, 1.0 );
   fEffBvsSLocal->SetDirectory(0); // it's local

   for (Int_t ibin=0; ibin<fNbins; ibin++) {
      istr >> tmpbin >> tmpeffS >> tmpeffB;
      fEffBvsSLocal->SetBinContent( ibin+1, tmpeffB );

      for (UInt_t ivar=0; ivar<GetNvar(); ivar++) {
         istr >> fCutMin[ivar][ibin] >> fCutMax[ivar][ibin];
      }
   }

   fEffSMin = fEffBvsSLocal->GetBinCenter(1);
   fEffSMax = fEffBvsSLocal->GetBinCenter(fNbins);
}

////////////////////////////////////////////////////////////////////////////////
/// create XML description for LD classification and regression
/// (for arbitrary number of output classes/targets)

void TMVA::MethodCuts::AddWeightsXMLTo( void* parent ) const
{
   // write all necessary information to the stream
   std::vector<Double_t> cutsMin;
   std::vector<Double_t> cutsMax;

   void* wght = gTools().AddChild(parent, "Weights");
   gTools().AddAttr( wght, "OptimisationMethod", (Int_t)fEffMethod);
   gTools().AddAttr( wght, "FitMethod",          (Int_t)fFitMethod );
   gTools().AddAttr( wght, "nbins",              fNbins );
   gTools().AddComment( wght, Form( "Below are the optimised cuts for %i variables: Format: ibin(hist) effS effB cutMin[ivar=0] cutMax[ivar=0] ... cutMin[ivar=n-1] cutMax[ivar=n-1]", GetNvar() ) );

   // NOTE: The signal efficiency written out into
   //       the weight file does not correspond to the center of the bin within which the
   //       background rejection is maximised (as before) but to the lower left edge of it.
   //       This is because the cut optimisation algorithm determines the best background
   //       rejection for all signal efficiencies belonging into a bin. Since the best background
   //       rejection is in general obtained for the lowest possible signal efficiency, the
   //       reference signal efficiency is the lowest value in the bin.

   for (Int_t ibin=0; ibin<fNbins; ibin++) {
      Double_t effS     = fEffBvsSLocal->GetBinCenter ( ibin + 1 );
      Double_t trueEffS = GetCuts( effS, cutsMin, cutsMax );
      if (TMath::Abs(trueEffS) < 1e-10) trueEffS = 0;

      void* binxml = gTools().AddChild( wght, "Bin" );
      gTools().AddAttr( binxml, "ibin", ibin+1   );
      gTools().AddAttr( binxml, "effS", trueEffS );
      gTools().AddAttr( binxml, "effB", fEffBvsSLocal->GetBinContent( ibin + 1 ) );
      void* cutsxml = gTools().AddChild( binxml, "Cuts" );
      for (UInt_t ivar=0; ivar<GetNvar(); ivar++) {
         gTools().AddAttr( cutsxml, Form( "cutMin_%i", ivar ), cutsMin[ivar] );
         gTools().AddAttr( cutsxml, Form( "cutMax_%i", ivar ), cutsMax[ivar] );
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// read coefficients from xml weight file

void TMVA::MethodCuts::ReadWeightsFromXML( void* wghtnode )
{
   // delete old min and max
   for (UInt_t i=0; i<GetNvar(); i++) {
      if (fCutMin[i] != 0) delete [] fCutMin[i];
      if (fCutMax[i] != 0) delete [] fCutMax[i];
   }
   if (fCutMin != 0) delete [] fCutMin;
   if (fCutMax != 0) delete [] fCutMax;

   Int_t tmpEffMethod, tmpFitMethod;
   gTools().ReadAttr( wghtnode, "OptimisationMethod", tmpEffMethod );
   gTools().ReadAttr( wghtnode, "FitMethod",          tmpFitMethod );
   gTools().ReadAttr( wghtnode, "nbins",              fNbins       );

   fEffMethod = (EEffMethod)tmpEffMethod;
   fFitMethod = (EFitMethodType)tmpFitMethod;

   // print some information
   if (fFitMethod == kUseMonteCarlo) {
      Log() << kINFO << "Read cuts optimised using sample of MC events" << Endl;
   }
   else if (fFitMethod == kUseMonteCarloEvents) {
      Log() << kINFO << "Read cuts optimised using sample of MC-Event events" << Endl;
   }
   else if (fFitMethod == kUseGeneticAlgorithm) {
      Log() << kINFO << "Read cuts optimised using Genetic Algorithm" << Endl;
   }
   else if (fFitMethod == kUseSimulatedAnnealing) {
      Log() << kINFO << "Read cuts optimised using Simulated Annealing algorithm" << Endl;
   }
   else if (fFitMethod == kUseEventScan) {
      Log() << kINFO << "Read cuts optimised using Full Event Scan" << Endl;
   }
   else {
      Log() << kWARNING << "unknown method: " << fFitMethod << Endl;
   }
   Log() << kINFO << "Reading " << fNbins << " signal efficiency bins for " << GetNvar() << " variables" << Endl;

   delete fEffBvsSLocal;
   fEffBvsSLocal = new TH1F( GetTestvarName() + "_effBvsSLocal",
                             TString(GetName()) + " efficiency of B vs S", fNbins, 0.0, 1.0 );
   fEffBvsSLocal->SetDirectory(0); // it's local
   for (Int_t ibin=1; ibin<=fNbins; ibin++) fEffBvsSLocal->SetBinContent( ibin, -0.1 ); // Init

   fCutMin = new Double_t*[GetNvar()];
   fCutMax = new Double_t*[GetNvar()];
   for (UInt_t i=0;i<GetNvar();i++) {
      fCutMin[i] = new Double_t[fNbins];
      fCutMax[i] = new Double_t[fNbins];
   }

   // read efficiencies and cuts
   Int_t   tmpbin;
   Float_t tmpeffS, tmpeffB;
   void* ch = gTools().GetChild(wghtnode,"Bin");
   while (ch) {
      //       if (strcmp(gTools().GetName(ch),"Bin") !=0) {
      //          ch = gTools().GetNextChild(ch);
      //          continue;
      //       }

      gTools().ReadAttr( ch, "ibin", tmpbin  );
      gTools().ReadAttr( ch, "effS", tmpeffS );
      gTools().ReadAttr( ch, "effB", tmpeffB );

      // sanity check
      if (tmpbin-1 >= fNbins || tmpbin-1 < 0) {
         Log() << kFATAL << "Mismatch in bins: " << tmpbin-1 << " >= " << fNbins << Endl;
      }

      fEffBvsSLocal->SetBinContent( tmpbin, tmpeffB );
      void* ct = gTools().GetChild(ch);
      for (UInt_t ivar=0; ivar<GetNvar(); ivar++) {
         gTools().ReadAttr( ct, Form( "cutMin_%i", ivar ), fCutMin[ivar][tmpbin-1] );
         gTools().ReadAttr( ct, Form( "cutMax_%i", ivar ), fCutMax[ivar][tmpbin-1] );
      }
      ch = gTools().GetNextChild(ch, "Bin");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// write histograms and PDFs to file for monitoring purposes

void TMVA::MethodCuts::WriteMonitoringHistosToFile( void ) const
{
   Log() << kINFO << "Write monitoring histograms to file: " << BaseDir()->GetPath() << Endl;

   fEffBvsSLocal->Write();

   // save reference histograms to file
   if (fEffMethod == kUsePDFs) {
      for (UInt_t ivar=0; ivar<GetNvar(); ivar++) {
         (*fVarHistS)[ivar]->Write();
         (*fVarHistB)[ivar]->Write();
         (*fVarHistS_smooth)[ivar]->Write();
         (*fVarHistB_smooth)[ivar]->Write();
         (*fVarPdfS)[ivar]->GetPDFHist()->Write();
         (*fVarPdfB)[ivar]->GetPDFHist()->Write();
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Overloaded function to create background efficiency (rejection) versus
/// signal efficiency plot (first call of this function).
///
/// The function returns the signal efficiency at background efficiency
/// indicated in theString
///
/// "theString" must have two entries:
///  - `[0]`: "Efficiency"
///  - `[1]`: the value of background efficiency at which the signal efficiency
///           is to be returned

Double_t TMVA::MethodCuts::GetTrainingEfficiency(const TString& theString)
{
   // parse input string for required background efficiency
   TList* list  = gTools().ParseFormatLine( theString );
   // sanity check
   if (list->GetSize() != 2) {
      Log() << kFATAL << "<GetTrainingEfficiency> wrong number of arguments"
            << " in string: " << theString
            << " | required format, e.g., Efficiency:0.05" << Endl;
      return -1;
   }

   Results* results = Data()->GetResults(GetMethodName(), Types::kTesting, GetAnalysisType());

   // that will be the value of the efficiency retured (does not affect
   // the efficiency-vs-bkg plot which is done anyway.
   Float_t effBref  = atof( ((TObjString*)list->At(1))->GetString() );

   delete list;

   // first round ? --> create histograms
   if (results->GetHist("EFF_BVSS_TR")==0) {

      if (fBinaryTreeS != 0) { delete fBinaryTreeS; fBinaryTreeS = 0; }
      if (fBinaryTreeB != 0) { delete fBinaryTreeB; fBinaryTreeB = 0; }

      fBinaryTreeS = new BinarySearchTree();
      fBinaryTreeS->Fill( GetEventCollection(Types::kTraining), fSignalClass );
      fBinaryTreeB = new BinarySearchTree();
      fBinaryTreeB->Fill( GetEventCollection(Types::kTraining), fBackgroundClass );
      // there is no really good equivalent to the fEffS; fEffB (efficiency vs cutvalue)
      // for the "Cuts" method (unless we had only one cut). Maybe later I might add here
      // histograms for each of the cuts...but this would require also a change in the
      // base class, and it is not really necessary, as we get exactly THIS info from the
      // "evaluateAllVariables" anyway.

      // now create efficiency curve: background versus signal
      TH1* eff_bvss_tr = new TH1F( GetTestvarName() + "_trainingEffBvsS", GetTestvarName() + "", fNbins, 0, 1 );
      for (Int_t ibin=1; ibin<=fNbins; ibin++) eff_bvss_tr->SetBinContent( ibin, -0.1 ); // Init
      TH1* rej_bvss_tr = new TH1F( GetTestvarName() + "_trainingRejBvsS", GetTestvarName() + "", fNbins, 0, 1 );
      for (Int_t ibin=1; ibin<=fNbins; ibin++) rej_bvss_tr->SetBinContent( ibin, 0. ); // Init
      results->Store(eff_bvss_tr, "EFF_BVSS_TR");
      results->Store(rej_bvss_tr, "REJ_BVSS_TR");

      // use root finder

      // make the background-vs-signal efficiency plot
      Double_t* tmpCutMin = new Double_t[GetNvar()];
      Double_t* tmpCutMax = new Double_t[GetNvar()];
      Int_t nFailedBins=0;
      for (Int_t bini=1; bini<=fNbins; bini++) {
         for (UInt_t ivar=0; ivar <GetNvar(); ivar++){
            tmpCutMin[ivar] = fCutMin[ivar][bini-1];
            tmpCutMax[ivar] = fCutMax[ivar][bini-1];
         }
         // find cut value corresponding to a given signal efficiency
         Double_t effS, effB;
         this->GetEffsfromSelection( &tmpCutMin[0], &tmpCutMax[0], effS, effB);
         // check that effS matches bini
         Int_t effBin = eff_bvss_tr->GetXaxis()->FindBin(effS);
         if (effBin != bini){
            Log()<< kVERBOSE << "unable to fill efficiency bin " << bini<< " " << effBin <<Endl;
            nFailedBins++;
         }
         else{
            // and fill histograms
            eff_bvss_tr->SetBinContent( bini, effB     );
            rej_bvss_tr->SetBinContent( bini, 1.0-effB );
         }
      }
      if (nFailedBins>0) Log()<< kWARNING << " unable to fill "<< nFailedBins <<" efficiency bins " <<Endl;

      delete [] tmpCutMin;
      delete [] tmpCutMax;

      // create splines for histogram
      fSplTrainEffBvsS = new TSpline1( "trainEffBvsS", new TGraph( eff_bvss_tr ) );
   }

   // must exist...
   if (NULL == fSplTrainEffBvsS) return 0.0;

   // now find signal efficiency that corresponds to required background efficiency
   Double_t effS = 0., effB, effS_ = 0., effB_ = 0.;
   Int_t    nbins_ = 1000;

   // loop over efficiency bins until the background eff. matches the requirement
   for (Int_t bini=1; bini<=nbins_; bini++) {
      // get corresponding signal and background efficiencies
      effS = (bini - 0.5)/Float_t(nbins_);
      effB = fSplTrainEffBvsS->Eval( effS );

      // find signal efficiency that corresponds to required background efficiency
      if ((effB - effBref)*(effB_ - effBref) < 0) break;
      effS_ = effS;
      effB_ = effB;
   }

   return 0.5*(effS + effS_);
}

////////////////////////////////////////////////////////////////////////////////
/// Overloaded function to create background efficiency (rejection) versus
/// signal efficiency plot (first call of this function).
///
/// The function returns the signal efficiency at background efficiency
/// indicated in theString
///
/// "theString" must have two entries:
///  - `[0]`: "Efficiency"
///  - `[1]`: the value of background efficiency at which the signal efficiency
///           is to be returned

Double_t TMVA::MethodCuts::GetEfficiency( const TString& theString, Types::ETreeType type, Double_t& effSerr )
{
   Data()->SetCurrentType(type);

   Results* results = Data()->GetResults( GetMethodName(), Types::kTesting, GetAnalysisType() );

   // parse input string for required background efficiency
   TList* list  = gTools().ParseFormatLine( theString, ":" );

   if (list->GetSize() > 2) {
      delete list;
      Log() << kFATAL << "<GetEfficiency> wrong number of arguments"
            << " in string: " << theString
            << " | required format, e.g., Efficiency:0.05, or empty string" << Endl;
      return -1;
   }

   // sanity check
   Bool_t computeArea = (list->GetSize() < 2); // the area is computed

   // that will be the value of the efficiency retured (does not affect
   // the efficiency-vs-bkg plot which is done anyway.
   Float_t effBref = (computeArea?1.:atof( ((TObjString*)list->At(1))->GetString() ));

   delete list;


   // first round ? --> create histograms
   if (results->GetHist("MVA_EFF_BvsS")==0) {

      if (fBinaryTreeS!=0) { delete fBinaryTreeS; fBinaryTreeS = 0; }
      if (fBinaryTreeB!=0) { delete fBinaryTreeB; fBinaryTreeB = 0; }

      // the variables may be transformed by a transformation method: to coherently
      // treat signal and background one must decide which transformation type shall
      // be used: our default is signal-type
      fBinaryTreeS = new BinarySearchTree();
      fBinaryTreeS->Fill( GetEventCollection(Types::kTesting), fSignalClass );
      fBinaryTreeB = new BinarySearchTree();
      fBinaryTreeB->Fill( GetEventCollection(Types::kTesting), fBackgroundClass );

      // there is no really good equivalent to the fEffS; fEffB (efficiency vs cutvalue)
      // for the "Cuts" method (unless we had only one cut). Maybe later I might add here
      // histograms for each of the cuts...but this would require also a change in the
      // base class, and it is not really necessary, as we get exactly THIS info from the
      // "evaluateAllVariables" anyway.

      // now create efficiency curve: background versus signal
      TH1* eff_BvsS = new TH1F( GetTestvarName() + "_effBvsS", GetTestvarName() + "", fNbins, 0, 1 );
      for (Int_t ibin=1; ibin<=fNbins; ibin++) eff_BvsS->SetBinContent( ibin, -0.1 ); // Init
      TH1* rej_BvsS = new TH1F( GetTestvarName() + "_rejBvsS", GetTestvarName() + "", fNbins, 0, 1 );
      for (Int_t ibin=1; ibin<=fNbins; ibin++) rej_BvsS->SetBinContent( ibin, 0.0 ); // Init
      results->Store(eff_BvsS, "MVA_EFF_BvsS");
      results->Store(rej_BvsS);

      Double_t xmin = 0.;
      Double_t xmax = 1.000001;

      TH1* eff_s = new TH1F( GetTestvarName() + "_effS", GetTestvarName() + " (signal)",     fNbins, xmin, xmax);
      for (Int_t ibin=1; ibin<=fNbins; ibin++) eff_s->SetBinContent( ibin, -0.1 ); // Init
      TH1* eff_b = new TH1F( GetTestvarName() + "_effB", GetTestvarName() + " (background)", fNbins, xmin, xmax);
      for (Int_t ibin=1; ibin<=fNbins; ibin++) eff_b->SetBinContent( ibin, -0.1 ); // Init
      results->Store(eff_s, "MVA_S");
      results->Store(eff_b, "MVA_B");

      // use root finder

      // make the background-vs-signal efficiency plot
      Double_t* tmpCutMin = new Double_t[GetNvar()];
      Double_t* tmpCutMax = new Double_t[GetNvar()];
      TGraph* tmpBvsS = new TGraph(fNbins+1);
      tmpBvsS->SetPoint(0, 0., 0.);

      for (Int_t bini=1; bini<=fNbins; bini++) {
         for (UInt_t ivar=0; ivar <GetNvar(); ivar++) {
            tmpCutMin[ivar] = fCutMin[ivar][bini-1];
            tmpCutMax[ivar] = fCutMax[ivar][bini-1];
         }
         // find cut value corresponding to a given signal efficiency
         Double_t effS, effB;
         this->GetEffsfromSelection( &tmpCutMin[0], &tmpCutMax[0], effS, effB);
         tmpBvsS->SetPoint(bini, effS, effB);

         eff_s->SetBinContent(bini, effS);
         eff_b->SetBinContent(bini, effB);
      }
      tmpBvsS->SetPoint(fNbins+1, 1., 1.);

      delete [] tmpCutMin;
      delete [] tmpCutMax;

      // create splines for histogram
      fSpleffBvsS = new TSpline1( "effBvsS", tmpBvsS );
      for (Int_t bini=1; bini<=fNbins; bini++) {
         Double_t effS = (bini - 0.5)/Float_t(fNbins);
         Double_t effB = fSpleffBvsS->Eval( effS );
         eff_BvsS->SetBinContent( bini, effB     );
         rej_BvsS->SetBinContent( bini, 1.0-effB );
      }
   }

   // must exist...
   if (NULL == fSpleffBvsS) return 0.0;

   // now find signal efficiency that corresponds to required background efficiency
   Double_t effS = 0, effB = 0, effS_ = 0, effB_ = 0;
   Int_t    nbins_ = 1000;

   if (computeArea) {

      // compute area of rej-vs-eff plot
      Double_t integral = 0;
      for (Int_t bini=1; bini<=nbins_; bini++) {

         // get corresponding signal and background efficiencies
         effS = (bini - 0.5)/Float_t(nbins_);
         effB = fSpleffBvsS->Eval( effS );
         integral += (1.0 - effB);
      }
      integral /= nbins_;

      return integral;
   }
   else {

      // loop over efficiency bins until the background eff. matches the requirement
      for (Int_t bini=1; bini<=nbins_; bini++) {
         // get corresponding signal and background efficiencies
         effS = (bini - 0.5)/Float_t(nbins_);
         effB = fSpleffBvsS->Eval( effS );

         // find signal efficiency that corresponds to required background efficiency
         if ((effB - effBref)*(effB_ - effBref) < 0) break;
         effS_ = effS;
         effB_ = effB;
      }

      effS    = 0.5*(effS + effS_);
      effSerr = 0;
      if (Data()->GetNEvtSigTest() > 0)
         effSerr = TMath::Sqrt( effS*(1.0 - effS)/Double_t(Data()->GetNEvtSigTest()) );

      return effS;

   }

   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// write specific classifier response

void TMVA::MethodCuts::MakeClassSpecific( std::ostream& fout, const TString& className ) const
{
   fout << "   // not implemented for class: \"" << className << "\"" << std::endl;
   fout << "};" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// get help message text
///
/// typical length of text line:
///         "|--------------------------------------------------------------|"

void TMVA::MethodCuts::GetHelpMessage() const
{
   TString bold    = gConfig().WriteOptionsReference() ? "<b>" : "";
   TString resbold = gConfig().WriteOptionsReference() ? "</b>" : "";
   TString brk     = gConfig().WriteOptionsReference() ? "<br>" : "";

   Log() << Endl;
   Log() << gTools().Color("bold") << "--- Short description:" << gTools().Color("reset") << Endl;
   Log() << Endl;
   Log() << "The optimisation of rectangular cuts performed by TMVA maximises " << Endl;
   Log() << "the background rejection at given signal efficiency, and scans " << Endl;
   Log() << "over the full range of the latter quantity. Three optimisation" << Endl;
   Log() << "methods are optional: Monte Carlo sampling (MC), a Genetics" << Endl;
   Log() << "Algorithm (GA), and Simulated Annealing (SA). GA and SA are"  << Endl;
   Log() << "expected to perform best." << Endl;
   Log() << Endl;
   Log() << "The difficulty to find the optimal cuts strongly increases with" << Endl;
   Log() << "the dimensionality (number of input variables) of the problem." << Endl;
   Log() << "This behavior is due to the non-uniqueness of the solution space."<<  Endl;
   Log() << Endl;
   Log() << gTools().Color("bold") << "--- Performance optimisation:" << gTools().Color("reset") << Endl;
   Log() << Endl;
   Log() << "If the dimensionality exceeds, say, 4 input variables, it is " << Endl;
   Log() << "advisable to scrutinize the separation power of the variables," << Endl;
   Log() << "and to remove the weakest ones. If some among the input variables" << Endl;
   Log() << "can be described by a single cut (e.g., because signal tends to be" << Endl;
   Log() << "larger than background), this can be indicated to MethodCuts via" << Endl;
   Log() << "the \"Fsmart\" options (see option string). Choosing this option" << Endl;
   Log() << "reduces the number of requirements for the variable from 2 (min/max)" << Endl;
   Log() << "to a single one (TMVA finds out whether it is to be interpreted as" << Endl;
   Log() << "min or max)." << Endl;
   Log() << Endl;
   Log() << gTools().Color("bold") << "--- Performance tuning via configuration options:" << gTools().Color("reset") << Endl;
   Log() << "" << Endl;
   Log() << bold << "Monte Carlo sampling:" << resbold << Endl;
   Log() << "" << Endl;
   Log() << "Apart form the \"Fsmart\" option for the variables, the only way" << Endl;
   Log() << "to improve the MC sampling is to increase the sampling rate. This" << Endl;
   Log() << "is done via the configuration option \"MC_NRandCuts\". The execution" << Endl;
   Log() << "time scales linearly with the sampling rate." << Endl;
   Log() << "" << Endl;
   Log() << bold << "Genetic Algorithm:" << resbold << Endl;
   Log() << "" << Endl;
   Log() << "The algorithm terminates if no significant fitness increase has" << Endl;
   Log() << "been achieved within the last \"nsteps\" steps of the calculation." << Endl;
   Log() << "Wiggles in the ROC curve or constant background rejection of 1" << Endl;
   Log() << "indicate that the GA failed to always converge at the true maximum" << Endl;
   Log() << "fitness. In such a case, it is recommended to broaden the search " << Endl;
   Log() << "by increasing the population size (\"popSize\") and to give the GA " << Endl;
   Log() << "more time to find improvements by increasing the number of steps" << Endl;
   Log() << "(\"nsteps\")" << Endl;
   Log() << "  -> increase \"popSize\" (at least >10 * number of variables)" << Endl;
   Log() << "  -> increase \"nsteps\"" << Endl;
   Log() << "" << Endl;
   Log() << bold << "Simulated Annealing (SA) algorithm:" << resbold << Endl;
   Log() << "" << Endl;
   Log() << "\"Increasing Adaptive\" approach:" << Endl;
   Log() << "" << Endl;
   Log() << "The algorithm seeks local minima and explores their neighborhoods, while" << Endl;
   Log() << "changing the ambient temperature depending on the number of failures" << Endl;
   Log() << "in the previous steps. The performance can be improved by increasing" << Endl;
   Log() << "the number of iteration steps (\"MaxCalls\"), or by adjusting the" << Endl;
   Log() << "minimal temperature (\"MinTemperature\"). Manual adjustments of the" << Endl;
   Log() << "speed of the temperature increase (\"TemperatureScale\" and \"AdaptiveSpeed\")" << Endl;
   Log() << "to individual data sets should also help. Summary:" << brk << Endl;
   Log() << "  -> increase \"MaxCalls\"" << brk << Endl;
   Log() << "  -> adjust   \"MinTemperature\"" << brk << Endl;
   Log() << "  -> adjust   \"TemperatureScale\"" << brk << Endl;
   Log() << "  -> adjust   \"AdaptiveSpeed\"" << Endl;
   Log() << "" << Endl;
   Log() << "\"Decreasing Adaptive\" approach:" << Endl;
   Log() << "" << Endl;
   Log() << "The algorithm calculates the initial temperature (based on the effect-" << Endl;
   Log() << "iveness of large steps) and the multiplier that ensures to reach the" << Endl;
   Log() << "minimal temperature with the requested number of iteration steps." << Endl;
   Log() << "The performance can be improved by adjusting the minimal temperature" << Endl;
   Log() << " (\"MinTemperature\") and by increasing number of steps (\"MaxCalls\"):" << brk << Endl;
   Log() << "  -> increase \"MaxCalls\"" << brk << Endl;
   Log() << "  -> adjust   \"MinTemperature\"" << Endl;
   Log() << " " << Endl;
   Log() << "Other kernels:" << Endl;
   Log() << "" << Endl;
   Log() << "Alternative ways of counting the temperature change are implemented. " << Endl;
   Log() << "Each of them starts with the maximum temperature (\"MaxTemperature\")" << Endl;
   Log() << "and decreases while changing the temperature according to a given" << Endl;
   Log() << "prescription:" << brk << Endl;
   Log() << "CurrentTemperature =" << brk << Endl;
   Log() << "  - Sqrt: InitialTemperature / Sqrt(StepNumber+2) * TemperatureScale" << brk << Endl;
   Log() << "  - Log:  InitialTemperature / Log(StepNumber+2) * TemperatureScale" << brk << Endl;
   Log() << "  - Homo: InitialTemperature / (StepNumber+2) * TemperatureScale" << brk << Endl;
   Log() << "  - Sin:  (Sin(StepNumber / TemperatureScale) + 1) / (StepNumber + 1)*InitialTemperature + Eps" << brk << Endl;
   Log() << "  - Geo:  CurrentTemperature * TemperatureScale" << Endl;
   Log() << "" << Endl;
   Log() << "Their performance can be improved by adjusting initial temperature" << Endl;
   Log() << "(\"InitialTemperature\"), the number of iteration steps (\"MaxCalls\")," << Endl;
   Log() << "and the multiplier that scales the temperature decrease" << Endl;
   Log() << "(\"TemperatureScale\")" << brk << Endl;
   Log() << "  -> increase \"MaxCalls\"" << brk << Endl;
   Log() << "  -> adjust   \"InitialTemperature\"" << brk << Endl;
   Log() << "  -> adjust   \"TemperatureScale\"" << brk << Endl;
   Log() << "  -> adjust   \"KernelTemperature\"" << Endl;
}
