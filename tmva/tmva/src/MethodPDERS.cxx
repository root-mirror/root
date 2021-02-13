// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Yair Mahalalel, Joerg Stelzer, Helge Voss, Kai Voss

/***********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis        *
 * Package: TMVA                                                                   *
 * Class  : MethodPDERS                                                            *
 * Web    : http://tmva.sourceforge.net                                            *
 *                                                                                 *
 * Description:                                                                    *
 *      Implementation                                                             *
 *                                                                                 *
 * Authors (alphabetical):                                                         *
 *      Krzysztof Danielowski <danielow@cern.ch>       - IFJ PAN & AGH, Poland     *
 *      Andreas Hoecker       <Andreas.Hocker@cern.ch> - CERN, Switzerland         *
 *      Kamil Kraszewski      <kalq@cern.ch>           - IFJ PAN & UJ, Poland      *
 *      Maciej Kruk           <mkruk@cern.ch>          - IFJ PAN & AGH, Poland     *
 *      Yair Mahalalel        <Yair.Mahalalel@cern.ch> - CERN, Switzerland         *
 *      Helge Voss            <Helge.Voss@cern.ch>     - MPI-K Heidelberg, Germany *
 *      Kai Voss              <Kai.Voss@cern.ch>       - U. of Victoria, Canada    *
 *                                                                                 *
 * Copyright (c) 2005:                                                             *
 *      CERN, Switzerland                                                          *
 *      U. of Victoria, Canada                                                     *
 *      MPI-K Heidelberg, Germany                                                  *
 *                                                                                 *
 * Redistribution and use in source and binary forms, with or without              *
 * modification, are permitted according to the terms listed in LICENSE            *
 * (http://tmva.sourceforge.net/LICENSE)                                           *
 ***********************************************************************************/

/*! \class TMVA::MethodPDERS
\ingroup TMVA

This is a generalization of the above Likelihood methods to \f$ N_{var} \f$
dimensions, where \f$ N_{var} \f$ is the number of input variables
used in the MVA. If the multi-dimensional probability density functions
(PDFs) for signal and background were known, this method contains the entire
physical information, and is therefore optimal. Usually, kernel estimation
methods are used to approximate the PDFs using the events from the
training sample.

A very simple probability density estimator (PDE) has been suggested
in [hep-ex/0211019](http://arxiv.org/abs/hep-ex/0211019). The
PDE for a given test event is obtained from counting the (normalized)
number of signal and background (training) events that occur in the
"vicinity" of the test event. The volume that describes "vicinity" is
user-defined. A [search method based on binary-trees](http://arxiv.org/abs/hep-ex/0211019)
is used to effectively reduce the
selection time for the range search. Three different volume definitions
are optional:

  - *MinMax:* the volume is defined in each dimension with respect
      to the full variable range found in the training sample.
  - *RMS:* the volume is defined in each dimensions with respect
      to the RMS estimated from the training sample.
  - *Adaptive:* a volume element is defined in each dimensions with
      respect to the RMS estimated from the training sample. The overall
      scale of the volume element is then determined for each event so
      that the total number of events confined in the volume be within
      a user-defined range.

The adaptive range search is used by default.
*/

#include "TMVA/MethodPDERS.h"

#include "TMVA/BinaryTree.h"
#include "TMVA/BinarySearchTree.h"
#include "TMVA/Configurable.h"
#include "TMVA/ClassifierFactory.h"
#include "TMVA/Event.h"
#include "TMVA/IMethod.h"
#include "TMVA/MethodBase.h"
#include "TMVA/MsgLogger.h"
#include "TMVA/RootFinder.h"
#include "TMVA/Tools.h"
#include "TMVA/TransformationHandler.h"
#include "TMVA/Types.h"

#include "ThreadLocalStorage.h"
#include "TFile.h"
#include "TMath.h"

#include <cassert>
#include <algorithm>

namespace TMVA {
   const Bool_t MethodPDERS_UseFindRoot = kFALSE;
};


REGISTER_METHOD(PDERS)

ClassImp(TMVA::MethodPDERS);

////////////////////////////////////////////////////////////////////////////////
/// standard constructor for the PDERS method

   TMVA::MethodPDERS::MethodPDERS( const TString& jobName,
                                   const TString& methodTitle,
                                   DataSetInfo& theData,
                                   const TString& theOption) :
   MethodBase( jobName, Types::kPDERS, methodTitle, theData, theOption),
   fFcnCall(0),
   fVRangeMode(kAdaptive),
   fKernelEstimator(kBox),
   fDelta(0),
   fShift(0),
   fScaleS(0),
   fScaleB(0),
   fDeltaFrac(0),
   fGaussSigma(0),
   fGaussSigmaNorm(0),
   fNRegOut(0),
   fNEventsMin(0),
   fNEventsMax(0),
   fMaxVIterations(0),
   fInitialScale(0),
   fInitializedVolumeEle(0),
   fkNNMin(0),
   fkNNMax(0),
   fMax_distance(0),
   fPrinted(0),
   fNormTree(0)
{
      fHelpVolume = NULL;
      fBinaryTree = NULL;
}

////////////////////////////////////////////////////////////////////////////////
/// construct MethodPDERS through from file

TMVA::MethodPDERS::MethodPDERS( DataSetInfo& theData,
                                const TString& theWeightFile) :
   MethodBase( Types::kPDERS, theData, theWeightFile),
   fFcnCall(0),
   fVRangeMode(kAdaptive),
   fKernelEstimator(kBox),
   fDelta(0),
   fShift(0),
   fScaleS(0),
   fScaleB(0),
   fDeltaFrac(0),
   fGaussSigma(0),
   fGaussSigmaNorm(0),
   fNRegOut(0),
   fNEventsMin(0),
   fNEventsMax(0),
   fMaxVIterations(0),
   fInitialScale(0),
   fInitializedVolumeEle(0),
   fkNNMin(0),
   fkNNMax(0),
   fMax_distance(0),
   fPrinted(0),
   fNormTree(0)
{
      fHelpVolume = NULL;
      fBinaryTree = NULL;
}

////////////////////////////////////////////////////////////////////////////////
/// PDERS can handle classification with 2 classes and regression with one or more regression-targets

Bool_t TMVA::MethodPDERS::HasAnalysisType( Types::EAnalysisType type, UInt_t numberClasses, UInt_t /*numberTargets*/ )
{
   if (type == Types::kClassification && numberClasses == 2) return kTRUE;
   if (type == Types::kRegression) return kTRUE;
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// default initialisation routine called by all constructors

void TMVA::MethodPDERS::Init( void )
{
   fBinaryTree = NULL;

   UpdateThis();

   // default options
   fDeltaFrac       = 3.0;
   fVRangeMode      = kAdaptive;
   fKernelEstimator = kBox;

   // special options for Adaptive mode
   fNEventsMin      = 100;
   fNEventsMax      = 200;
   fMaxVIterations  = 150;
   fInitialScale    = 0.99;
   fGaussSigma      = 0.1;
   fNormTree        = kFALSE;

   fkNNMin      = Int_t(fNEventsMin);
   fkNNMax      = Int_t(fNEventsMax);

   fInitializedVolumeEle = kFALSE;
   fAverageRMS.clear();

   // the minimum requirement to declare an event signal-like
   SetSignalReferenceCut( 0.0 );
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TMVA::MethodPDERS::~MethodPDERS( void )
{
   if (fDelta) delete fDelta;
   if (fShift) delete fShift;

   if (NULL != fBinaryTree) delete fBinaryTree;
}

////////////////////////////////////////////////////////////////////////////////
/// define the options (their key words) that can be set in the option string.
///
/// know options:
///  - VolumeRangeMode   <string>  Method to determine volume range
///    available values are:
///    - MinMax
///    - Unscaled
///    - RMS
///    - kNN
///    - Adaptive <default>
///
///  - KernelEstimator   <string>  Kernel estimation function
///    available values are:
///    - Box <default>
///    - Sphere
///    - Teepee
///    - Gauss
///    - Sinc3
///    - Sinc5
///    - Sinc7
///    - Sinc9
///    - Sinc11
///    - Lanczos2
///    - Lanczos3
///    - Lanczos5
///    - Lanczos8
///    - Trim
///
///  - DeltaFrac         <float>   Ratio of #EventsMin/#EventsMax for MinMax and RMS volume range
///  - NEventsMin        <int>     Minimum number of events for adaptive volume range
///  - NEventsMax        <int>     Maximum number of events for adaptive volume range
///  - MaxVIterations    <int>     Maximum number of iterations for adaptive volume range
///  - InitialScale      <float>   Initial scale for adaptive volume range
///  - GaussSigma        <float>   Width with respect to the volume size of Gaussian kernel estimator

void TMVA::MethodPDERS::DeclareOptions()
{
   DeclareOptionRef(fVolumeRange="Adaptive", "VolumeRangeMode", "Method to determine volume size");
   AddPreDefVal(TString("Unscaled"));
   AddPreDefVal(TString("MinMax"));
   AddPreDefVal(TString("RMS"));
   AddPreDefVal(TString("Adaptive"));
   AddPreDefVal(TString("kNN"));

   DeclareOptionRef(fKernelString="Box", "KernelEstimator", "Kernel estimation function");
   AddPreDefVal(TString("Box"));
   AddPreDefVal(TString("Sphere"));
   AddPreDefVal(TString("Teepee"));
   AddPreDefVal(TString("Gauss"));
   AddPreDefVal(TString("Sinc3"));
   AddPreDefVal(TString("Sinc5"));
   AddPreDefVal(TString("Sinc7"));
   AddPreDefVal(TString("Sinc9"));
   AddPreDefVal(TString("Sinc11"));
   AddPreDefVal(TString("Lanczos2"));
   AddPreDefVal(TString("Lanczos3"));
   AddPreDefVal(TString("Lanczos5"));
   AddPreDefVal(TString("Lanczos8"));
   AddPreDefVal(TString("Trim"));

   DeclareOptionRef(fDeltaFrac     , "DeltaFrac",      "nEventsMin/Max for minmax and rms volume range");
   DeclareOptionRef(fNEventsMin    , "NEventsMin",     "nEventsMin for adaptive volume range");
   DeclareOptionRef(fNEventsMax    , "NEventsMax",     "nEventsMax for adaptive volume range");
   DeclareOptionRef(fMaxVIterations, "MaxVIterations", "MaxVIterations for adaptive volume range");
   DeclareOptionRef(fInitialScale  , "InitialScale",   "InitialScale for adaptive volume range");
   DeclareOptionRef(fGaussSigma    , "GaussSigma",     "Width (wrt volume size) of Gaussian kernel estimator");
   DeclareOptionRef(fNormTree      , "NormTree",       "Normalize binary search tree");
}

////////////////////////////////////////////////////////////////////////////////
/// process the options specified by the user

void TMVA::MethodPDERS::ProcessOptions()
{
   if (IgnoreEventsWithNegWeightsInTraining()) {
      Log() << kFATAL << "Mechanism to ignore events with negative weights in training not yet available for method: "
            << GetMethodTypeName()
            << " --> please remove \"IgnoreNegWeightsInTraining\" option from booking string."
            << Endl;
   }

   fGaussSigmaNorm = fGaussSigma; // * TMath::Sqrt( Double_t(GetNvar()) );

   fVRangeMode = MethodPDERS::kUnsupported;

   if      (fVolumeRange == "MinMax"    ) fVRangeMode = kMinMax;
   else if (fVolumeRange == "RMS"       ) fVRangeMode = kRMS;
   else if (fVolumeRange == "Adaptive"  ) fVRangeMode = kAdaptive;
   else if (fVolumeRange == "Unscaled"  ) fVRangeMode = kUnscaled;
   else if (fVolumeRange == "kNN"   ) fVRangeMode = kkNN;
   else {
      Log() << kFATAL << "VolumeRangeMode parameter '" << fVolumeRange << "' unknown" << Endl;
   }

   if      (fKernelString == "Box"      ) fKernelEstimator = kBox;
   else if (fKernelString == "Sphere"   ) fKernelEstimator = kSphere;
   else if (fKernelString == "Teepee"   ) fKernelEstimator = kTeepee;
   else if (fKernelString == "Gauss"    ) fKernelEstimator = kGauss;
   else if (fKernelString == "Sinc3"    ) fKernelEstimator = kSinc3;
   else if (fKernelString == "Sinc5"    ) fKernelEstimator = kSinc5;
   else if (fKernelString == "Sinc7"    ) fKernelEstimator = kSinc7;
   else if (fKernelString == "Sinc9"    ) fKernelEstimator = kSinc9;
   else if (fKernelString == "Sinc11"   ) fKernelEstimator = kSinc11;
   else if (fKernelString == "Lanczos2" ) fKernelEstimator = kLanczos2;
   else if (fKernelString == "Lanczos3" ) fKernelEstimator = kLanczos3;
   else if (fKernelString == "Lanczos5" ) fKernelEstimator = kLanczos5;
   else if (fKernelString == "Lanczos8" ) fKernelEstimator = kLanczos8;
   else if (fKernelString == "Trim"     ) fKernelEstimator = kTrim;
   else {
      Log() << kFATAL << "KernelEstimator parameter '" << fKernelString << "' unknown" << Endl;
   }

   // TODO: Add parameter validation

   Log() << kVERBOSE << "interpreted option string: vRangeMethod: '"
         << (const char*)((fVRangeMode == kMinMax) ? "MinMax" :
                          (fVRangeMode == kUnscaled) ? "Unscaled" :
                          (fVRangeMode == kRMS   ) ? "RMS" : "Adaptive") << "'" << Endl;
   if (fVRangeMode == kMinMax || fVRangeMode == kRMS)
      Log() << kVERBOSE << "deltaFrac: " << fDeltaFrac << Endl;
   else
      Log() << kVERBOSE << "nEventsMin/Max, maxVIterations, initialScale: "
            << fNEventsMin << "  " << fNEventsMax
            << "  " << fMaxVIterations << "  " << fInitialScale << Endl;
   Log() << kVERBOSE << "KernelEstimator = " << fKernelString << Endl;
}

////////////////////////////////////////////////////////////////////////////////
/// this is a dummy training: the preparation work to do is the construction
/// of the binary tree as a pointer chain. It is easier to directly save the
/// trainingTree in the weight file, and to rebuild the binary tree in the
/// test phase from scratch

void TMVA::MethodPDERS::Train( void )
{
   if (IsNormalised()) Log() << kFATAL << "\"Normalise\" option cannot be used with PDERS; "
                             << "please remove the option from the configuration string, or "
                             << "use \"!Normalise\""
                             << Endl;

   CreateBinarySearchTree( Types::kTraining );

   CalcAverages();
   SetVolumeElement();

   fInitializedVolumeEle = kTRUE;
   ExitFromTraining();
}

////////////////////////////////////////////////////////////////////////////////
/// init the size of a volume element using a defined fraction of the
/// volume containing the entire events

Double_t TMVA::MethodPDERS::GetMvaValue( Double_t* err, Double_t* errUpper )
{
   if (fInitializedVolumeEle == kFALSE) {
      fInitializedVolumeEle = kTRUE;

      // binary trees must exist
      assert( fBinaryTree );

      CalcAverages();
      SetVolumeElement();
   }

   // cannot determine error
   NoErrorCalc(err, errUpper);

   return this->CRScalc( *GetEvent() );
}

////////////////////////////////////////////////////////////////////////////////

const std::vector< Float_t >& TMVA::MethodPDERS::GetRegressionValues()
{
   if (fRegressionReturnVal == 0) fRegressionReturnVal = new std::vector<Float_t>;
   fRegressionReturnVal->clear();
   // init the size of a volume element using a defined fraction of the
   // volume containing the entire events
   if (fInitializedVolumeEle == kFALSE) {
      fInitializedVolumeEle = kTRUE;

      // binary trees must exist
      assert( fBinaryTree );

      CalcAverages();

      SetVolumeElement();
   }

   const Event* ev = GetEvent();
   this->RRScalc( *ev, fRegressionReturnVal );

   Event * evT = new Event(*ev);
   UInt_t ivar = 0;
   for (std::vector<Float_t>::iterator it = fRegressionReturnVal->begin(); it != fRegressionReturnVal->end(); ++it ) {
      evT->SetTarget(ivar,(*it));
      ivar++;
   }

   const Event* evT2 = GetTransformationHandler().InverseTransform( evT );
   fRegressionReturnVal->clear();

   for (ivar = 0; ivar<evT2->GetNTargets(); ivar++) {
      fRegressionReturnVal->push_back(evT2->GetTarget(ivar));
   }

   delete evT;


   return (*fRegressionReturnVal);
}

////////////////////////////////////////////////////////////////////////////////
/// compute also average RMS values required for adaptive Gaussian

void TMVA::MethodPDERS::CalcAverages()
{
   if (fVRangeMode == kAdaptive || fVRangeMode == kRMS || fVRangeMode == kkNN  ) {
      fAverageRMS.clear();
      fBinaryTree->CalcStatistics();

      for (UInt_t ivar=0; ivar<GetNvar(); ivar++) {
         if (!DoRegression()){ //why there are separate rms for signal and background?
            Float_t rmsS = fBinaryTree->RMS(Types::kSignal, ivar);
            Float_t rmsB = fBinaryTree->RMS(Types::kBackground, ivar);
            fAverageRMS.push_back( (rmsS + rmsB)*0.5 );
         } else {
            Float_t rms = fBinaryTree->RMS( ivar );
            fAverageRMS.push_back( rms );
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// create binary search trees for signal and background

void TMVA::MethodPDERS::CreateBinarySearchTree( Types::ETreeType type )
{
   if (NULL != fBinaryTree) delete fBinaryTree;
   fBinaryTree = new BinarySearchTree();
   if (fNormTree) {
      fBinaryTree->SetNormalize( kTRUE );
   }

   fBinaryTree->Fill( GetEventCollection(type) );

   if (fNormTree) {
      fBinaryTree->NormalizeTree();
   }

   if (!DoRegression()) {
      // these are the signal and background scales for the weights
      fScaleS = 1.0/fBinaryTree->GetSumOfWeights( Types::kSignal );
      fScaleB = 1.0/fBinaryTree->GetSumOfWeights( Types::kBackground );

      Log() << kVERBOSE << "Signal and background scales: " << fScaleS << " " << fScaleB << Endl;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// defines volume dimensions

void TMVA::MethodPDERS::SetVolumeElement( void ) {
   if (GetNvar()==0) {
      Log() << kFATAL << "GetNvar() == 0" << Endl;
      return;
   }

   // init relative scales
   fkNNMin      = Int_t(fNEventsMin);
   fkNNMax      = Int_t(fNEventsMax);

   if (fDelta) delete fDelta;
   if (fShift) delete fShift;
   fDelta = new std::vector<Float_t>( GetNvar() );
   fShift = new std::vector<Float_t>( GetNvar() );

   for (UInt_t ivar=0; ivar<GetNvar(); ivar++) {
      switch (fVRangeMode) {

      case kRMS:
      case kkNN:
      case kAdaptive:
         // sanity check
         if (fAverageRMS.size() != GetNvar())
            Log() << kFATAL << "<SetVolumeElement> RMS not computed: " << fAverageRMS.size() << Endl;
         (*fDelta)[ivar] = fAverageRMS[ivar]*fDeltaFrac;
         Log() << kVERBOSE << "delta of var[" << (*fInputVars)[ivar]
               << "\t]: " << fAverageRMS[ivar]
               << "\t  |  comp with |max - min|: " << (GetXmax( ivar ) - GetXmin( ivar ))
               << Endl;
         break;
      case kMinMax:
         (*fDelta)[ivar] = (GetXmax( ivar ) - GetXmin( ivar ))*fDeltaFrac;
         break;
      case kUnscaled:
         (*fDelta)[ivar] = fDeltaFrac;
         break;
      default:
         Log() << kFATAL << "<SetVolumeElement> unknown range-set mode: "
               << fVRangeMode << Endl;
      }
      (*fShift)[ivar] = 0.5; // volume is centered around test value
   }

}

////////////////////////////////////////////////////////////////////////////////
/// Interface to RootFinder

Double_t TMVA::MethodPDERS::IGetVolumeContentForRoot( Double_t scale )
{
   return ThisPDERS()->GetVolumeContentForRoot( scale );
}

////////////////////////////////////////////////////////////////////////////////
/// count number of events in rescaled volume

Double_t TMVA::MethodPDERS::GetVolumeContentForRoot( Double_t scale )
{
   Volume v( *fHelpVolume );
   v.ScaleInterval( scale );

   Double_t count = GetBinaryTree()->SearchVolume( &v );

   v.Delete();
   return count;
}

void TMVA::MethodPDERS::GetSample( const Event& e,
                                   std::vector<const BinarySearchTreeNode*>& events,
                                   Volume *volume )
{
   Float_t count = 0;

   // -------------------------------------------------------------------------
   //
   // ==== test of volume search =====
   //
   // #define TMVA::MethodPDERS__countByHand__Debug__

#ifdef  TMVA_MethodPDERS__countByHand__Debug__

   // starting values
   count = fBinaryTree->SearchVolume( volume );

   Int_t iS = 0, iB = 0;
   UInt_t nvar = GetNvar();
   for (UInt_t ievt_=0; ievt_<Data()->GetNTrainingEvents(); ievt_++) {
      const Event * ev = GetTrainingEvent(ievt_);
      Bool_t inV;
      for (Int_t ivar=0; ivar<nvar; ivar++) {
         Float_t x = ev->GetValue(ivar);
         inV = (x > (*volume->Lower)[ivar] && x <= (*volume->Upper)[ivar]);
         if (!inV) break;
      }
      if (inV) {
         in++;
      }
   }
   Log() << kVERBOSE << "debug: my test: " << in << Endl;// <- ***********tree
   Log() << kVERBOSE << "debug: binTree: " << count << Endl << Endl;// <- ***********tree

#endif

   // -------------------------------------------------------------------------

   if (fVRangeMode == kRMS || fVRangeMode == kMinMax || fVRangeMode == kUnscaled) { // Constant volume

      std::vector<Double_t> *lb = new std::vector<Double_t>( GetNvar() );
      for (UInt_t ivar=0; ivar<GetNvar(); ivar++) (*lb)[ivar] = e.GetValue(ivar);
      std::vector<Double_t> *ub = new std::vector<Double_t>( *lb );
      for (UInt_t ivar=0; ivar<GetNvar(); ivar++) {
         (*lb)[ivar] -= (*fDelta)[ivar]*(1.0 - (*fShift)[ivar]);
         (*ub)[ivar] += (*fDelta)[ivar]*(*fShift)[ivar];
      }
      Volume* svolume = new Volume( lb, ub );
      // starting values

      fBinaryTree->SearchVolume( svolume, &events );
   }
   else if (fVRangeMode == kAdaptive) {      // adaptive volume

      // -----------------------------------------------------------------------

      // TODO: optimize, perhaps multi stage with broadening limits,
      // or a different root finding method entirely,
      if (MethodPDERS_UseFindRoot) {

         // that won't need to search through large volume, where the bottle neck probably is

         fHelpVolume = volume;

         UpdateThis(); // necessary update of static pointer
         RootFinder rootFinder( this, 0.01, 50, 200, 10 );
         Double_t scale = rootFinder.Root( (fNEventsMin + fNEventsMax)/2.0 );

         volume->ScaleInterval( scale );

         fBinaryTree->SearchVolume( volume, &events );

         fHelpVolume = NULL;
      }
      // -----------------------------------------------------------------------
      else {

         // starting values
         count = fBinaryTree->SearchVolume( volume );

         Float_t nEventsO = count;
         Int_t i_=0;

         while (nEventsO < fNEventsMin) { // this isn't a sain start... try again
            volume->ScaleInterval( 1.15 );
            count = fBinaryTree->SearchVolume( volume );
            nEventsO = count;
            i_++;
         }
         if (i_ > 50) Log() << kWARNING << "warning in event: " << e
                            << ": adaptive volume pre-adjustment reached "
                            << ">50 iterations in while loop (" << i_ << ")" << Endl;

         Float_t nEventsN    = nEventsO;
         Float_t nEventsE    = 0.5*(fNEventsMin + fNEventsMax);
         Float_t scaleO      = 1.0;
         Float_t scaleN      = fInitialScale;
         Float_t scale       = scaleN;
         Float_t scaleBest   = scaleN;
         Float_t nEventsBest = nEventsN;

         for (Int_t ic=1; ic<fMaxVIterations; ic++) {
            if (nEventsN < fNEventsMin || nEventsN > fNEventsMax) {

               // search for events in rescaled volume
               Volume* v = new Volume( *volume );
               v->ScaleInterval( scale );
               nEventsN  = fBinaryTree->SearchVolume( v );

               // determine next iteration (linear approximation)
               if (nEventsN > 1 && nEventsN - nEventsO != 0)
                  if (scaleN - scaleO != 0)
                     scale += (scaleN - scaleO)/(nEventsN - nEventsO)*(nEventsE - nEventsN);
                  else
                     scale += (-0.01); // should not actually occur...
               else
                  scale += 0.5; // use much larger volume

               // save old scale
               scaleN   = scale;

               // take if better (don't accept it if too small number of events)
               if (TMath::Abs(nEventsN - nEventsE) < TMath::Abs(nEventsBest - nEventsE) &&
                   (nEventsN >= fNEventsMin || nEventsBest < nEventsN)) {
                  nEventsBest = nEventsN;
                  scaleBest   = scale;
               }

               v->Delete();
               delete v;
            }
            else break;
         }

         // last sanity check
         nEventsN = nEventsBest;
         // include "1" to cover float precision
         if (nEventsN < fNEventsMin-1 || nEventsN > fNEventsMax+1)
            Log() << kWARNING << "warning in event " << e
                  << ": adaptive volume adjustment reached "
                  << "max. #iterations (" << fMaxVIterations << ")"
                  << "[ nEvents: " << nEventsN << "  " << fNEventsMin << "  " << fNEventsMax << "]"
                  << Endl;

         volume->ScaleInterval( scaleBest );
         fBinaryTree->SearchVolume( volume, &events );
      }

      // end of adaptive method

   } else if (fVRangeMode == kkNN) {
      Volume v(*volume);

      events.clear();
      // check number of signals in beginning volume
      Int_t kNNcount = fBinaryTree->SearchVolumeWithMaxLimit( &v, &events, fkNNMax+1 );
      //if this number is too large return fkNNMax+1

      Int_t t_times = 0;  // number of iterations

      while ( !(kNNcount >= fkNNMin && kNNcount <= fkNNMax) ) {
         if (kNNcount < fkNNMin) {         //if we have too less points
            Float_t scale = 2;      //better scale needed
            volume->ScaleInterval( scale );
         }
         else if (kNNcount > fkNNMax) {    //if we have too many points
            Float_t scale = 0.1;      //better scale needed
            volume->ScaleInterval( scale );
         }
         events.clear();

         v = *volume ;

         kNNcount = fBinaryTree->SearchVolumeWithMaxLimit( &v, &events, fkNNMax+1 );  //new search

         t_times++;

         if (t_times == fMaxVIterations) {
            Log() << kWARNING << "warning in event" << e
                  << ": kNN volume adjustment reached "
                  << "max. #iterations (" << fMaxVIterations << ")"
                  << "[ kNN: " << fkNNMin << " " << fkNNMax << Endl;
            break;
         }
      }

      //vector to normalize distance in each dimension
      Double_t *dim_normalization = new Double_t[GetNvar()];
      for (UInt_t ivar=0; ivar<GetNvar(); ivar++) {
         dim_normalization [ivar] = 1.0 / ((*v.fUpper)[ivar] - (*v.fLower)[ivar]);
      }

      std::vector<const BinarySearchTreeNode*> tempVector;    // temporary vector for events

      if (kNNcount >= fkNNMin) {
         std::vector<Double_t> *distances = new std::vector<Double_t>( kNNcount );

         //counting the distance for each event
         for (Int_t j=0;j< Int_t(events.size()) ;j++)
            (*distances)[j] = GetNormalizedDistance ( e, *events[j], dim_normalization );

         //counting the fkNNMin-th element
         std::vector<Double_t>::iterator wsk = distances->begin();
         for (Int_t j=0;j<fkNNMin-1;++j) ++wsk;
         std::nth_element( distances->begin(), wsk, distances->end() );

         //getting all elements that are closer than fkNNMin-th element
         //signals
         for (Int_t j=0;j<Int_t(events.size());j++) {
            Double_t dist = GetNormalizedDistance( e, *events[j], dim_normalization );

            if (dist <= (*distances)[fkNNMin-1])
               tempVector.push_back( events[j] );
         }
         fMax_distance = (*distances)[fkNNMin-1];
         delete distances;
      }
      delete[] dim_normalization;
      events = tempVector;

   } else {

      // troubles ahead...
      Log() << kFATAL << "<GetSample> unknown RangeMode: " << fVRangeMode << Endl;
   }
   // -----------------------------------------------------------------------
}

////////////////////////////////////////////////////////////////////////////////

Double_t TMVA::MethodPDERS::CRScalc( const Event& e )
{
   std::vector<const BinarySearchTreeNode*> events;

   // computes event weight by counting number of signal and background
   // events (of reference sample) that are found within given volume
   // defined by the event
   std::vector<Double_t> *lb = new std::vector<Double_t>( GetNvar() );
   for (UInt_t ivar=0; ivar<GetNvar(); ivar++) (*lb)[ivar] = e.GetValue(ivar);

   std::vector<Double_t> *ub = new std::vector<Double_t>( *lb );
   for (UInt_t ivar=0; ivar<GetNvar(); ivar++) {
      (*lb)[ivar] -= (*fDelta)[ivar]*(1.0 - (*fShift)[ivar]);
      (*ub)[ivar] += (*fDelta)[ivar]*(*fShift)[ivar];
   }

   Volume *volume = new Volume( lb, ub );

   GetSample( e, events, volume );
   Double_t count = CKernelEstimate( e, events, *volume );
   delete volume;
   delete lb;
   delete ub;

   return count;
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::MethodPDERS::RRScalc( const Event& e, std::vector<Float_t>* count )
{
   std::vector<const BinarySearchTreeNode*> events;

   // computes event weight by counting number of signal and background
   // events (of reference sample) that are found within given volume
   // defined by the event
   std::vector<Double_t> *lb = new std::vector<Double_t>( GetNvar() );
   for (UInt_t ivar=0; ivar<GetNvar(); ivar++) (*lb)[ivar] = e.GetValue(ivar);

   std::vector<Double_t> *ub = new std::vector<Double_t>( *lb );
   for (UInt_t ivar=0; ivar<GetNvar(); ivar++) {
      (*lb)[ivar] -= (*fDelta)[ivar]*(1.0 - (*fShift)[ivar]);
      (*ub)[ivar] += (*fDelta)[ivar]*(*fShift)[ivar];
   }
   Volume *volume = new Volume( lb, ub );

   GetSample( e, events, volume );
   RKernelEstimate( e, events, *volume, count );

   delete volume;

   return;
}
////////////////////////////////////////////////////////////////////////////////
/// normalization factors so we can work with radius 1 hyperspheres

Double_t TMVA::MethodPDERS::CKernelEstimate( const Event & event,
                                             std::vector<const BinarySearchTreeNode*>& events, Volume& v )
{
   Double_t *dim_normalization = new Double_t[GetNvar()];
   for (UInt_t ivar=0; ivar<GetNvar(); ivar++)
      dim_normalization [ivar] = 2 / ((*v.fUpper)[ivar] - (*v.fLower)[ivar]);

   Double_t pdfSumS = 0;
   Double_t pdfSumB = 0;

   // Iteration over sample points
   for (std::vector<const BinarySearchTreeNode*>::iterator iev = events.begin(); iev != events.end(); ++iev) {

      // First switch to the one dimensional distance
      Double_t normalized_distance = GetNormalizedDistance (event, *(*iev), dim_normalization);

      // always working within the hyperelipsoid, except for when we don't
      // note that rejection ratio goes to 1 as nvar goes to infinity
      if (normalized_distance > 1 && fKernelEstimator != kBox) continue;

      if ( (*iev)->GetClass()==fSignalClass )
         pdfSumS += ApplyKernelFunction (normalized_distance) * (*iev)->GetWeight();
      else
         pdfSumB += ApplyKernelFunction (normalized_distance) * (*iev)->GetWeight();
   }
   pdfSumS = KernelNormalization( pdfSumS < 0. ? 0. : pdfSumS );
   pdfSumB = KernelNormalization( pdfSumB < 0. ? 0. : pdfSumB );

   delete[] dim_normalization;

   if (pdfSumS < 1e-20 && pdfSumB < 1e-20) return 0.5;
   if (pdfSumB < 1e-20) return 1.0;
   if (pdfSumS < 1e-20) return 0.0;

   Float_t r = pdfSumB*fScaleB/(pdfSumS*fScaleS);
   return 1.0/(r + 1.0);   // TODO: propagate errors from here
}

////////////////////////////////////////////////////////////////////////////////
/// normalization factors so we can work with radius 1 hyperspheres

void TMVA::MethodPDERS::RKernelEstimate( const Event & event,
                                         std::vector<const BinarySearchTreeNode*>& events, Volume& v,
                                         std::vector<Float_t>* pdfSum )
{
   Double_t *dim_normalization = new Double_t[GetNvar()];
   for (UInt_t ivar=0; ivar<GetNvar(); ivar++)
      dim_normalization [ivar] = 2 / ((*v.fUpper)[ivar] - (*v.fLower)[ivar]);

   //   std::vector<Float_t> pdfSum;
   pdfSum->clear();
   Float_t pdfDiv = 0;
   fNRegOut = 1; // for now, regression is just for 1 dimension

   for (Int_t ivar = 0; ivar < fNRegOut ; ivar++)
      pdfSum->push_back( 0 );

   // Iteration over sample points
   for (std::vector<const BinarySearchTreeNode*>::iterator iev = events.begin(); iev != events.end(); ++iev) {

      // First switch to the one dimensional distance
      Double_t normalized_distance = GetNormalizedDistance (event, *(*iev), dim_normalization);

      // always working within the hyperelipsoid, except for when we don't
      // note that rejection ratio goes to 1 as nvar goes to infinity
      if (normalized_distance > 1 && fKernelEstimator != kBox) continue;

      for (Int_t ivar = 0; ivar < fNRegOut ; ivar++) {
         pdfSum->at(ivar) += ApplyKernelFunction (normalized_distance) * (*iev)->GetWeight() * (*iev)->GetTargets()[ivar];
         pdfDiv           += ApplyKernelFunction (normalized_distance) * (*iev)->GetWeight();
      }
   }

   delete[] dim_normalization;

   if (pdfDiv == 0)
      return;

   for (Int_t ivar = 0; ivar < fNRegOut ; ivar++)
      pdfSum->at(ivar) /= pdfDiv;

   return;
}

////////////////////////////////////////////////////////////////////////////////
/// from the normalized euclidean distance calculate the distance
/// for a certain kernel

Double_t TMVA::MethodPDERS::ApplyKernelFunction (Double_t normalized_distance)
{
   switch (fKernelEstimator) {
   case kBox:
   case kSphere:
      return 1;
      break;
   case kTeepee:
      return (1 - normalized_distance);
      break;
   case kGauss:
      return TMath::Gaus( normalized_distance, 0, fGaussSigmaNorm, kFALSE);
      break;
   case kSinc3:
   case kSinc5:
   case kSinc7:
   case kSinc9:
   case kSinc11: {
      Double_t side_crossings = 2 + ((int) fKernelEstimator) - ((int) kSinc3);
      return NormSinc (side_crossings * normalized_distance);
   }
      break;
   case kLanczos2:
      return LanczosFilter (2, normalized_distance);
      break;
   case kLanczos3:
      return LanczosFilter (3, normalized_distance);
      break;
   case kLanczos5:
      return LanczosFilter (5, normalized_distance);
      break;
   case kLanczos8:
      return LanczosFilter (8, normalized_distance);
      break;
   case kTrim: {
      Double_t x = normalized_distance / fMax_distance;
      x = 1 - x*x*x;
      return x*x*x;
   }
      break;
   default:
      Log() << kFATAL << "Kernel estimation function unsupported. Enumerator is " << fKernelEstimator << Endl;
      break;
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Calculating the normalization factor only once (might need a reset at some point.
/// Can the method be restarted with different params?)

Double_t TMVA::MethodPDERS::KernelNormalization (Double_t pdf)
{
   // Caching jammed to disable function.
   // It's not really useful after all, badly implemented and untested :-)
   TTHREAD_TLS(Double_t) ret = 1.0;

   if (ret != 0.0) return ret*pdf;

   // We first normalize by the volume of the hypersphere.
   switch (fKernelEstimator) {
   case kBox:
   case kSphere:
      ret = 1.;
      break;
   case kTeepee:
      ret =   (GetNvar() * (GetNvar() + 1) * TMath::Gamma (((Double_t) GetNvar()) / 2.)) /
         ( TMath::Power (2., (Double_t) GetNvar() + 1) * TMath::Power (TMath::Pi(), ((Double_t) GetNvar()) / 2.));
      break;
   case kGauss:
      // We use full range integral here. Reasonable because of the fast function decay.
      ret = 1. / TMath::Power ( 2 * TMath::Pi() * fGaussSigmaNorm * fGaussSigmaNorm, ((Double_t) GetNvar()) / 2.);
      break;
   case kSinc3:
   case kSinc5:
   case kSinc7:
   case kSinc9:
   case kSinc11:
   case kLanczos2:
   case kLanczos3:
   case kLanczos5:
   case kLanczos8:
      // We use the full range integral here. Reasonable because the central lobe dominates it.
      ret = 1 / TMath::Power ( 2., (Double_t) GetNvar() );
      break;
   default:
      Log() << kFATAL << "Kernel estimation function unsupported. Enumerator is " << fKernelEstimator << Endl;
   }

   // Normalizing by the full volume
   ret *= ( TMath::Power (2., static_cast<Int_t>(GetNvar())) * TMath::Gamma (1 + (((Double_t) GetNvar()) / 2.)) ) /
      TMath::Power (TMath::Pi(), ((Double_t) GetNvar()) / 2.);

   return ret*pdf;
}

////////////////////////////////////////////////////////////////////////////////
/// We use Euclidian metric here. Might not be best or most efficient.

Double_t TMVA::MethodPDERS::GetNormalizedDistance ( const Event &base_event,
                                                    const BinarySearchTreeNode &sample_event,
                                                    Double_t *dim_normalization)
{
   Double_t ret=0;
   for (UInt_t ivar=0; ivar<GetNvar(); ivar++) {
      Double_t dist = dim_normalization[ivar] * (sample_event.GetEventV()[ivar] - base_event.GetValue(ivar));
      ret += dist*dist;
   }
   // DD 26.11.2008: bugfix: division by (GetNvar()) was missing for correct normalisation
   ret /= GetNvar();
   return TMath::Sqrt (ret);
}

////////////////////////////////////////////////////////////////////////////////
/// NormSinc

Double_t TMVA::MethodPDERS::NormSinc (Double_t x)
{
   if (x < 10e-10 && x > -10e-10) {
      return 1; // Poor man's l'Hopital
   }

   Double_t pix = TMath::Pi() * x;
   Double_t sinc = TMath::Sin(pix) / pix;
   Double_t ret;

   if (GetNvar() % 2)
      ret = TMath::Power (sinc, static_cast<Int_t>(GetNvar()));
   else
      ret = TMath::Abs (sinc) * TMath::Power (sinc, static_cast<Int_t>(GetNvar() - 1));

   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// Lanczos Filter

Double_t TMVA::MethodPDERS::LanczosFilter (Int_t level, Double_t x)
{
   if (x < 10e-10 && x > -10e-10) {
      return 1; // Poor man's l'Hopital
   }

   Double_t pix = TMath::Pi() * x;
   Double_t pixtimesn = pix * ((Double_t) level);
   Double_t lanczos = (TMath::Sin(pix) / pix) * (TMath::Sin(pixtimesn) / pixtimesn);
   Double_t ret;

   if (GetNvar() % 2) ret = TMath::Power (lanczos, static_cast<Int_t>(GetNvar()));
   else               ret = TMath::Abs (lanczos) * TMath::Power (lanczos, static_cast<Int_t>(GetNvar() - 1));

   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// statistical error estimate for RS estimator

Float_t TMVA::MethodPDERS::GetError( Float_t countS, Float_t countB,
                                     Float_t sumW2S, Float_t sumW2B ) const
{
   Float_t c = fScaleB/fScaleS;
   Float_t d = countS + c*countB; d *= d;

   if (d < 1e-10) return 1; // Error is zero because of B = S = 0

   Float_t f = c*c/d/d;
   Float_t err = f*countB*countB*sumW2S + f*countS*countS*sumW2B;

   if (err < 1e-10) return 1; // Error is zero because of B or S = 0

   return sqrt(err);
}

////////////////////////////////////////////////////////////////////////////////
/// write weights to xml file

void TMVA::MethodPDERS::AddWeightsXMLTo( void* parent ) const
{
   void* wght = gTools().AddChild(parent, "Weights");
   if (fBinaryTree)
      fBinaryTree->AddXMLTo(wght);
   else
      Log() << kFATAL << "Signal and background binary search tree not available" << Endl;
   //Log() << kFATAL << "Please implement writing of weights as XML" << Endl;
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::MethodPDERS::ReadWeightsFromXML( void* wghtnode)
{
   if (NULL != fBinaryTree) delete fBinaryTree;
   void* treenode = gTools().GetChild(wghtnode);
   fBinaryTree = TMVA::BinarySearchTree::CreateFromXML(treenode);
   if(!fBinaryTree)
      Log() << kFATAL << "Could not create BinarySearchTree from XML" << Endl;
   if(!fBinaryTree)
      Log() << kFATAL << "Could not create BinarySearchTree from XML" << Endl;
   fBinaryTree->SetPeriode( GetNvar() );
   fBinaryTree->CalcStatistics();
   fBinaryTree->CountNodes();
   if (fBinaryTree->GetSumOfWeights( Types::kSignal ) > 0)
      fScaleS = 1.0/fBinaryTree->GetSumOfWeights( Types::kSignal );
   else fScaleS = 1;
   if (fBinaryTree->GetSumOfWeights( Types::kBackground ) > 0)
      fScaleB = 1.0/fBinaryTree->GetSumOfWeights( Types::kBackground );
   else fScaleB = 1;
   Log() << kINFO << "signal and background scales: " << fScaleS << " " << fScaleB << Endl;
   CalcAverages();
   SetVolumeElement();
   fInitializedVolumeEle = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// read weight info from file

void TMVA::MethodPDERS::ReadWeightsFromStream( std::istream& istr)
{
   if (NULL != fBinaryTree) delete fBinaryTree;

   fBinaryTree = new BinarySearchTree();

   istr >> *fBinaryTree;

   fBinaryTree->SetPeriode( GetNvar() );

   fBinaryTree->CalcStatistics();

   fBinaryTree->CountNodes();

   // these are the signal and background scales for the weights
   fScaleS = 1.0/fBinaryTree->GetSumOfWeights( Types::kSignal );
   fScaleB = 1.0/fBinaryTree->GetSumOfWeights( Types::kBackground );

   Log() << kINFO << "signal and background scales: " << fScaleS << " " << fScaleB << Endl;

   CalcAverages();

   SetVolumeElement();

   fInitializedVolumeEle = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// write training sample (TTree) to file

void TMVA::MethodPDERS::WriteWeightsToStream( TFile& ) const
{
}

////////////////////////////////////////////////////////////////////////////////
/// read training sample from file

void TMVA::MethodPDERS::ReadWeightsFromStream( TFile& /*rf*/ )
{
}

////////////////////////////////////////////////////////////////////////////////
/// static pointer to this object

TMVA::MethodPDERS* TMVA::MethodPDERS::ThisPDERS( void )
{
   return GetMethodPDERSThreadLocal();
}
////////////////////////////////////////////////////////////////////////////////
/// update static this pointer

void TMVA::MethodPDERS::UpdateThis( void )
{
   GetMethodPDERSThreadLocal() = this;
}

////////////////////////////////////////////////////////////////////////////////
/// write specific classifier response

void TMVA::MethodPDERS::MakeClassSpecific( std::ostream& fout, const TString& className ) const
{
   fout << "   // not implemented for class: \"" << className << "\"" << std::endl;
   fout << "};" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// get help message text
///
/// typical length of text line:
///         "|--------------------------------------------------------------|"

void TMVA::MethodPDERS::GetHelpMessage() const
{
   Log() << Endl;
   Log() << gTools().Color("bold") << "--- Short description:" << gTools().Color("reset") << Endl;
   Log() << Endl;
   Log() << "PDERS is a generalization of the projective likelihood classifier " << Endl;
   Log() << "to N dimensions, where N is the number of input variables used." << Endl;
   Log() << "In its adaptive form it is mostly equivalent to k-Nearest-Neighbor" << Endl;
   Log() << "(k-NN) methods. If the multidimensional PDF for signal and background" << Endl;
   Log() << "were known, this classifier would exploit the full information" << Endl;
   Log() << "contained in the input variables, and would hence be optimal. In " << Endl;
   Log() << "practice however, huge training samples are necessary to sufficiently " << Endl;
   Log() << "populate the multidimensional phase space. " << Endl;
   Log() << Endl;
   Log() << "The simplest implementation of PDERS counts the number of signal" << Endl;
   Log() << "and background events in the vicinity of a test event, and returns" << Endl;
   Log() << "a weight according to the majority species of the neighboring events." << Endl;
   Log() << "A more involved version of PDERS (selected by the option \"KernelEstimator\")" << Endl;
   Log() << "uses Kernel estimation methods to approximate the shape of the PDF." << Endl;
   Log() << Endl;
   Log() << gTools().Color("bold") << "--- Performance optimisation:" << gTools().Color("reset") << Endl;
   Log() << Endl;
   Log() << "PDERS can be very powerful in case of strongly non-linear problems, " << Endl;
   Log() << "e.g., distinct islands of signal and background regions. Because of " << Endl;
   Log() << "the exponential growth of the phase space, it is important to restrict" << Endl;
   Log() << "the number of input variables (dimension) to the strictly necessary." << Endl;
   Log() << Endl;
   Log() << "Note that PDERS is a slowly responding classifier. Moreover, the necessity" << Endl;
   Log() << "to store the entire binary tree in memory, to avoid accessing virtual " << Endl;
   Log() << "memory, limits the number of training events that can effectively be " << Endl;
   Log() << "used to model the multidimensional PDF." << Endl;
   Log() << Endl;
   Log() << gTools().Color("bold") << "--- Performance tuning via configuration options:" << gTools().Color("reset") << Endl;
   Log() << Endl;
   Log() << "If the PDERS response is found too slow when using the adaptive volume " << Endl;
   Log() << "size (option \"VolumeRangeMode=Adaptive\"), it might be found beneficial" << Endl;
   Log() << "to reduce the number of events required in the volume, and/or to enlarge" << Endl;
   Log() << "the allowed range (\"NeventsMin/Max\"). PDERS is relatively insensitive" << Endl;
   Log() << "to the width (\"GaussSigma\") of the Gaussian kernel (if used)." << Endl;
}
