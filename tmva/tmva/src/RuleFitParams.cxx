// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Joerg Stelzer, Fredrik Tegenfeldt, Helge Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : RuleFitParams                                                         *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Implementation                                                            *
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

/*! \class TMVA::RuleFitParams
\ingroup TMVA
A class doing the actual fitting of a linear model using rules as base functions.
*/
#include "TMVA/RuleFitParams.h"

#include "TMVA/DataSetInfo.h"
#include "TMVA/MethodRuleFit.h"
#include "TMVA/MsgLogger.h"
#include "TMVA/RuleEnsemble.h"
#include "TMVA/RuleFit.h"
#include "TMVA/Timer.h"
#include "TMVA/Tools.h"
#include "TMVA/Types.h"

#include "TTree.h"
#include "TMath.h"

#include <iomanip>
#include <numeric>
#include <algorithm>
#include <functional>

Bool_t gFIRSTTST=kTRUE;
Bool_t gFIRSTORG=kTRUE;

Double_t gGDInit=0;
Double_t gGDPtr=0;
Double_t gGDNorm=0;
Double_t gGDEval=0;
Double_t gGDEvalRule=0;
Double_t gGDRuleLoop=0;
Double_t gGDLinLoop=0;

////////////////////////////////////////////////////////////////////////////////
/// constructor

TMVA::RuleFitParams::RuleFitParams()
   : fRuleFit ( 0 )
   , fRuleEnsemble ( 0 )
   , fNRules ( 0 )
   , fNLinear ( 0 )
   , fPathIdx1 ( 0 )
   , fPathIdx2 ( 0 )
   , fPerfIdx1 ( 0 )
   , fPerfIdx2 ( 0 )
   , fGDNTauTstOK( 0 )
   , fGDNTau     ( 51 )
   , fGDTauPrec  ( 0.02 )
   , fGDTauScan  ( 1000 )
   , fGDTauMin   ( 0.0 )
   , fGDTauMax   ( 1.0 )
   , fGDTau      ( -1.0 )
   , fGDPathStep ( 0.01 )
   , fGDNPathSteps ( 1000 )
   , fGDErrScale ( 1.1 )
   , fAverageTruth( 0 )
   , fFstarMedian ( 0 )
   , fGDNtuple    ( 0 )
   , fNTRisk      ( 0 )
   , fNTErrorRate ( 0 )
   , fNTNuval     ( 0 )
   , fNTCoefRad   ( 0 )
   , fNTOffset ( 0 )
   , fNTCoeff ( 0 )
   , fNTLinCoeff ( 0 )
   , fsigave( 0 )
   , fsigrms( 0 )
   , fbkgave( 0 )
   , fbkgrms( 0 )
   , fLogger( new MsgLogger("RuleFit") )
{
   Init();
}
////////////////////////////////////////////////////////////////////////////////
/// destructor

TMVA::RuleFitParams::~RuleFitParams()
{
   if (fNTCoeff)     { delete [] fNTCoeff; fNTCoeff = nullptr; }
   if (fNTLinCoeff)  { delete [] fNTLinCoeff; fNTLinCoeff = nullptr; }
   delete fLogger;
}

////////////////////////////////////////////////////////////////////////////////
/// Initializes all parameters using the RuleEnsemble and the training tree

void TMVA::RuleFitParams::Init()
{
   if (fRuleFit==0) return;
   if (fRuleFit->GetMethodRuleFit()==0) {
      Log() << kFATAL << "RuleFitParams::Init() - MethodRuleFit ptr is null" << Endl;
   }
   UInt_t neve = fRuleFit->GetTrainingEvents().size();
   //
   fRuleEnsemble   = fRuleFit->GetRuleEnsemblePtr();
   fNRules         = fRuleEnsemble->GetNRules();
   fNLinear        = fRuleEnsemble->GetNLinear();

   //
   // Fraction of events used for validation should be close of unity..
   // Always selection from the END
   //
   UInt_t ofs;
   fPerfIdx1 = 0;
   if (neve>1) {
      fPerfIdx2 = static_cast<UInt_t>((neve-1)*fRuleFit->GetMethodRuleFit()->GetGDValidEveFrac());
   }
   else {
      fPerfIdx2 = 0;
   }
   ofs = neve - fPerfIdx2 - 1;
   fPerfIdx1 += ofs;
   fPerfIdx2 += ofs;
   //
   // Fraction of events used for the path search can be allowed to be a smaller value, say 0.5
   // Alwas select events from the BEGINNING.
   // This means that the validation and search samples will not overlap if both fractions are <0.5.
   //
   fPathIdx1 = 0;
   if (neve>1) {
      fPathIdx2 = static_cast<UInt_t>((neve-1)*fRuleFit->GetMethodRuleFit()->GetGDPathEveFrac());
   }
   else {
      fPathIdx2 = 0;
   }
   //
   // summarize weights
   //
   fNEveEffPath = 0;;
   for (UInt_t ie=fPathIdx1; ie<fPathIdx2+1; ie++) {
      fNEveEffPath += fRuleFit->GetTrainingEventWeight(ie);
   }

   fNEveEffPerf=0;
   for (UInt_t ie=fPerfIdx1; ie<fPerfIdx2+1; ie++) {
      fNEveEffPerf += fRuleFit->GetTrainingEventWeight(ie);
   }
   //
   Log() << kVERBOSE << "Path constr. - event index range = [ " << fPathIdx1 << ", " << fPathIdx2 << " ]"
         << ", effective N(events) = " << fNEveEffPath << Endl;
   Log() << kVERBOSE << "Error estim. - event index range = [ " << fPerfIdx1 << ", " << fPerfIdx2 << " ]"
         << ", effective N(events) = " << fNEveEffPerf << Endl;
   //
   if (fRuleEnsemble->DoRules())
      Log() << kDEBUG << "Number of rules in ensemble: " << fNRules << Endl;
   else
      Log() << kDEBUG << "Rules are disabled " << Endl;

   if (fRuleEnsemble->DoLinear())
      Log() << kDEBUG << "Number of linear terms: " << fNLinear << Endl;
   else
      Log() << kDEBUG << "Linear terms are disabled " << Endl;
}

////////////////////////////////////////////////////////////////////////////////
/// initializes the ntuple

void TMVA::RuleFitParams::InitNtuple()
{
   fGDNtuple= new TTree("MonitorNtuple_RuleFitParams","RuleFit path search");
   fGDNtuple->Branch("risk",    &fNTRisk,     "risk/D");
   fGDNtuple->Branch("error",   &fNTErrorRate,"error/D");
   fGDNtuple->Branch("nuval",   &fNTNuval,    "nuval/D");
   fGDNtuple->Branch("coefrad", &fNTCoefRad,  "coefrad/D");
   fGDNtuple->Branch("offset",  &fNTOffset,   "offset/D");
   //
   fNTCoeff    = (fNRules >0 ? new Double_t[fNRules]  : 0);
   fNTLinCoeff = (fNLinear>0 ? new Double_t[fNLinear] : 0);

   for (UInt_t i=0; i<fNRules; i++) {
      fGDNtuple->Branch(Form("a%d",i+1),&fNTCoeff[i],Form("a%d/D",i+1));
   }
   for (UInt_t i=0; i<fNLinear; i++) {
      fGDNtuple->Branch(Form("b%d",i+1),&fNTLinCoeff[i],Form("b%d/D",i+1));
   }
}

////////////////////////////////////////////////////////////////////////////////
/// evaluate the average of each variable and f(x) in the given range

void TMVA::RuleFitParams::EvaluateAverage( UInt_t ind1, UInt_t ind2,
                                           std::vector<Double_t> &avsel,
                                           std::vector<Double_t> &avrul )
{
   UInt_t neve = ind2-ind1+1;
   if (neve<1) {
      Log() << kFATAL << "<EvaluateAverage> - no events selected for path search -> BUG!" << Endl;
   }
   //
   avsel.clear();
   avrul.clear();
   //
   if (fNLinear>0) avsel.resize(fNLinear,0);
   if (fNRules>0)  avrul.resize(fNRules,0);
   const std::vector<UInt_t> *eventRuleMap=0;
   Double_t ew;
   Double_t sumew=0;
   //
   // Loop over events and calculate average of linear terms (normalised) and rule response.
   //
   if (fRuleEnsemble->IsRuleMapOK()) { // MakeRuleMap() has been called
      for ( UInt_t i=ind1; i<ind2+1; i++) {
         ew = fRuleFit->GetTrainingEventWeight(i);
         sumew += ew;
         for ( UInt_t sel=0; sel<fNLinear; sel++ ) {
            avsel[sel] += ew*fRuleEnsemble->EvalLinEvent(i,sel);
         }
         // loop over rules
         UInt_t nrules=0;
         if (fRuleEnsemble->DoRules()) {
            eventRuleMap = &(fRuleEnsemble->GetEventRuleMap(i));
            nrules = (*eventRuleMap).size();
         }
         for (UInt_t r=0; r<nrules; r++) {
            avrul[(*eventRuleMap)[r]] += ew;
         }
      }
   }
   else { // MakeRuleMap() has not yet been called
      const std::vector<const Event *> *events = &(fRuleFit->GetTrainingEvents());
      for ( UInt_t i=ind1; i<ind2+1; i++) {
         ew = fRuleFit->GetTrainingEventWeight(i);
         sumew += ew;
         // first cache rule/lin response
         /* Double_t val = */ fRuleEnsemble->EvalLinEvent(*((*events)[i]));
         /* val = */ fRuleEnsemble->EvalEvent(*((*events)[i]));
         // loop over linear terms
         for ( UInt_t sel=0; sel<fNLinear; sel++ ) {
            avsel[sel] += ew*fRuleEnsemble->GetEventLinearValNorm(sel);
         }
         // loop over rules
         for (UInt_t r=0; r<fNRules; r++) {
            avrul[r] += ew*fRuleEnsemble->GetEventRuleVal(r);
         }
      }
   }
   // average variable
   for ( UInt_t sel=0; sel<fNLinear; sel++ ) {
      avsel[sel] = avsel[sel] / sumew;
   }
   // average rule response, excl coeff
   for (UInt_t r=0; r<fNRules; r++) {
      avrul[r] = avrul[r] / sumew;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Implementation of squared-error ramp loss function (eq 39,40 in ref 1)
/// This is used for binary Classifications where y = {+1,-1} for (sig,bkg)

Double_t TMVA::RuleFitParams::LossFunction( const Event& e ) const
{
   Double_t h = TMath::Max( -1.0, TMath::Min(1.0,fRuleEnsemble->EvalEvent( e )) );
   Double_t diff = (fRuleFit->GetMethodRuleFit()->DataInfo().IsSignal(&e)?1:-1) - h;
   //
   return diff*diff*e.GetWeight();
}

////////////////////////////////////////////////////////////////////////////////
/// Implementation of squared-error ramp loss function (eq 39,40 in ref 1)
/// This is used for binary Classifications where y = {+1,-1} for (sig,bkg)

Double_t TMVA::RuleFitParams::LossFunction( UInt_t evtidx ) const
{
   Double_t h = TMath::Max( -1.0, TMath::Min(1.0,fRuleEnsemble->EvalEvent( evtidx )) );
   Double_t diff = (fRuleFit->GetMethodRuleFit()->DataInfo().IsSignal(fRuleEnsemble->GetRuleMapEvent( evtidx ))?1:-1) - h;
   //
   return diff*diff*fRuleFit->GetTrainingEventWeight(evtidx);
}

////////////////////////////////////////////////////////////////////////////////
/// Implementation of squared-error ramp loss function (eq 39,40 in ref 1)
/// This is used for binary Classifications where y = {+1,-1} for (sig,bkg)

Double_t TMVA::RuleFitParams::LossFunction( UInt_t evtidx, UInt_t itau ) const
{
   Double_t e = fRuleEnsemble->EvalEvent( evtidx , fGDOfsTst[itau], fGDCoefTst[itau], fGDCoefLinTst[itau]);
   Double_t h = TMath::Max( -1.0, TMath::Min(1.0,e) );
   Double_t diff = (fRuleFit->GetMethodRuleFit()->DataInfo().IsSignal(fRuleEnsemble->GetRuleMapEvent( evtidx ))?1:-1) - h;
   //
   return diff*diff*fRuleFit->GetTrainingEventWeight(evtidx);
}

////////////////////////////////////////////////////////////////////////////////
/// risk assessment

Double_t TMVA::RuleFitParams::Risk(UInt_t ind1,UInt_t ind2, Double_t neff) const
{
   UInt_t neve = ind2-ind1+1;
   if (neve<1) {
      Log() << kFATAL << "<Risk> Invalid start/end indices! BUG!!!" << Endl;
   }
   Double_t rval=0;
   //
   //   const std::vector<Event *> *events = &(fRuleFit->GetTrainingEvents());
   for ( UInt_t i=ind1; i<ind2+1; i++) {
      rval += LossFunction(i);
   }
   rval  = rval/neff;

   return rval;
}

////////////////////////////////////////////////////////////////////////////////
/// risk assessment for tau model <itau>

Double_t TMVA::RuleFitParams::Risk(UInt_t ind1,UInt_t ind2, Double_t neff, UInt_t itau) const
{
   UInt_t neve = ind2-ind1+1;
   if (neve<1) {
      Log() << kFATAL << "<Risk> Invalid start/end indices! BUG!!!" << Endl;
   }
   Double_t rval=0;
   //
   //   const std::vector<Event *> *events = &(fRuleFit->GetTrainingEvents());
   for ( UInt_t i=ind1; i<ind2+1; i++) {
      rval += LossFunction(i,itau);
   }
   rval  = rval/neff;

   return rval;
}

////////////////////////////////////////////////////////////////////////////////
/// This is the "lasso" penalty
/// To be used for regression.
/// --- NOT USED ---

Double_t TMVA::RuleFitParams::Penalty() const
{
   Log() << kWARNING << "<Penalty> Using unverified code! Check!" << Endl;
   Double_t rval=0;
   const std::vector<Double_t> *lincoeff = & (fRuleEnsemble->GetLinCoefficients());
   for (UInt_t i=0; i<fNRules; i++) {
      rval += TMath::Abs(fRuleEnsemble->GetRules(i)->GetCoefficient());
   }
   for (UInt_t i=0; i<fNLinear; i++) {
      rval += TMath::Abs((*lincoeff)[i]);
   }
   return rval;
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize GD path search

void TMVA::RuleFitParams::InitGD()
{
   if (fGDNTau<2) {
      fGDNTau    = 1;
      fGDTauScan = 0;
   }
   if (fGDTau<0.0) {
      //      fGDNTau    = 50; already set in MethodRuleFit
      fGDTauScan = 1000;
      fGDTauMin  = 0.0;
      fGDTauMax  = 1.0;
   }
   else {
      fGDNTau    = 1;
      fGDTauScan = 0;
   }
   // set all taus
   fGDTauVec.clear();
   fGDTauVec.resize( fGDNTau );
   if (fGDNTau==1) {
      fGDTauVec[0] = fGDTau;
   }
   else {
      // set tau vector - TODO: make a smarter choice of range in tau
      Double_t dtau = (fGDTauMax - fGDTauMin)/static_cast<Double_t>(fGDNTau-1);
      for (UInt_t itau=0; itau<fGDNTau; itau++) {
         fGDTauVec[itau] = static_cast<Double_t>(itau)*dtau + fGDTauMin;
         if (fGDTauVec[itau]>1.0) fGDTauVec[itau]=1.0;
      }
   }
   // initialize path search vectors

   fGradVec.clear();
   fGradVecLin.clear();
   fGradVecTst.clear();
   fGradVecLinTst.clear();
   fGDErrTst.clear();
   fGDErrTstOK.clear();
   fGDOfsTst.clear();
   fGDCoefTst.clear();
   fGDCoefLinTst.clear();
   //
   // rules
   //
   fGDCoefTst.resize(fGDNTau);
   fGradVec.resize(fNRules,0);
   fGradVecTst.resize(fGDNTau);
   for (UInt_t i=0; i<fGDNTau; i++) {
      fGradVecTst[i].resize(fNRules,0);
      fGDCoefTst[i].resize(fNRules,0);
   }
   //
   // linear terms
   //
   fGDCoefLinTst.resize(fGDNTau);
   fGradVecLin.resize(fNLinear,0);
   fGradVecLinTst.resize(fGDNTau);
   for (UInt_t i=0; i<fGDNTau; i++) {
      fGradVecLinTst[i].resize(fNLinear,0);
      fGDCoefLinTst[i].resize(fNLinear,0);
   }
   //
   // error, coefs etc
   //
   fGDErrTst.resize(fGDNTau,0);
   fGDErrTstOK.resize(fGDNTau,kTRUE);
   fGDOfsTst.resize(fGDNTau,0);
   fGDNTauTstOK = fGDNTau;
   //
   // calculate average selectors and rule responses for the path sample size
   //
}

////////////////////////////////////////////////////////////////////////////////
/// This finds the cutoff parameter tau by scanning several different paths

Int_t TMVA::RuleFitParams::FindGDTau()
{
   if (fGDNTau<2) return 0;
   if (fGDTauScan==0) return 0;

   if (fGDOfsTst.size()<1)
      Log() << kFATAL << "BUG! FindGDTau() has been called BEFORE InitGD()." << Endl;
   //
   Log() << kINFO << "Estimating the cutoff parameter tau. The estimated time is a pessimistic maximum." << Endl;
   //
   // Find how many points to scan and how often to calculate the error
   UInt_t nscan = fGDTauScan; //std::min(static_cast<Int_t>(fGDTauScan),fGDNPathSteps);
   UInt_t netst = std::min(nscan,UInt_t(100));
   UInt_t nscanned=0;
   //
   //--------------------
   // loop over the paths
   //--------------------
   // The number of MAXIMUM loops is given by nscan.
   // At each loop, the paths being far away from the minimum
   // are rejected. Hence at each check (every netst events), the number
   // of paths searched will be reduced.
   // The maximum 'distance' from the minimum error rate is
   // 1 sigma. See RiskPerfTst() for details.
   //
   Bool_t doloop=kTRUE;
   UInt_t ip=0;
   UInt_t itauMin=0;
   Timer timer( nscan, "RuleFit" );
   while (doloop) {
      // make gradvec
      MakeTstGradientVector();
      // update coefs
      UpdateTstCoefficients();
      // estimate error and do the sum
      // do this at index=0, netst-1, 2*netst-1 ...
      nscanned++;
      if ( (ip==0) || ((ip+1)%netst==0) ) {
         //         ErrorRateRocTst( );
         itauMin = RiskPerfTst();
         Log() << kVERBOSE << Form("%4d",ip+1) << ". tau = " << Form("%4.4f",fGDTauVec[itauMin])
               << " => error rate = " << fGDErrTst[itauMin] << Endl;
      }
      ip++;
      doloop = ((ip<nscan) && (fGDNTauTstOK>3));
      gFIRSTTST=kFALSE;
      if (Log().GetMinType()>kVERBOSE)
         timer.DrawProgressBar(ip);
   }
   //
   // Set tau and coefs
   // Downscale tau slightly in order to avoid numerical problems
   //
   if (nscanned==0) {
      Log() << kERROR << "<FindGDTau> number of scanned loops is zero! Should NOT see this message." << Endl;
   }
   fGDTau = fGDTauVec[itauMin];
   fRuleEnsemble->SetCoefficients( fGDCoefTst[itauMin] );
   fRuleEnsemble->SetLinCoefficients( fGDCoefLinTst[itauMin] );
   fRuleEnsemble->SetOffset( fGDOfsTst[itauMin] );
   Log() << kINFO << "Best path found with tau = " << Form("%4.4f",fGDTau)
         << " after " << timer.GetElapsedTime() << "      " << Endl;

   return nscan;
}

////////////////////////////////////////////////////////////////////////////////
/// The following finds the gradient directed path in parameter space.
/// More work is needed... FT, 24/9/2006
///
/// The algorithm is currently as follows (if not otherwise stated, the sample
/// used below is [fPathIdx1,fPathIdx2]):
///
/// 1. Set offset to -average(y(true)) and all coefs=0 => average of F(x)==0
/// 2. FindGDTau() : start scanning using several paths defined by different tau
///                  choose the tau yielding the best path
/// 3. start the scanning the chosen path
/// 4. check error rate at a given frequency
///    data used for check: [fPerfIdx1,fPerfIdx2]
/// 5. stop when either of the following conditions are fullfilled:
///    1. loop index==fGDNPathSteps
///    2. error > fGDErrScale*errmin
///    3. only in DEBUG mode: risk is not monotonously decreasing
///
/// The algorithm will warn if:
///   1. the error rate was still decreasing when loop finished -> increase fGDNPathSteps!
///   2. minimum was found at an early stage -> decrease fGDPathStep
///   3. DEBUG: risk > previous risk -> entered chaotic region (regularization is too small)

void TMVA::RuleFitParams::MakeGDPath()
{
   Log() << kINFO << "GD path scan - the scan stops when the max num. of steps is reached or a min is found"
         << Endl;
   Log() << kVERBOSE << "Number of events used per path step = " << fPathIdx2-fPathIdx1+1 << Endl;
   Log() << kVERBOSE << "Number of events used for error estimation = " << fPerfIdx2-fPerfIdx1+1 << Endl;

   // check if debug mode
   const Bool_t isVerbose = (Log().GetMinType()<=kVERBOSE);
   const Bool_t isDebug   = (Log().GetMinType()<=kDEBUG);

   // init GD parameters and clear coeff vectors
   InitGD();

   // evaluate average response of rules/linear terms (with event weights)
   EvaluateAveragePath();
   EvaluateAveragePerf();

   // initial estimate; all other a(i) are zero
   Log() << kVERBOSE << "Creating GD path"  << Endl;
   Log() << kVERBOSE << "  N(steps)     = "   << fGDNPathSteps << Endl;
   Log() << kVERBOSE << "  step         = "   << fGDPathStep   << Endl;
   Log() << kVERBOSE << "  N(tau)       = "   << fGDNTau       << Endl;
   Log() << kVERBOSE << "  N(tau steps) = "   << fGDTauScan    << Endl;
   Log() << kVERBOSE << "  tau range    = [ " << fGDTauVec[0]  << " , " << fGDTauVec[fGDNTau-1] << " ]" << Endl;

   // init ntuple
   if (isDebug) InitNtuple();

   // DEBUG: risk scan
   Int_t    nbadrisk=0;                // number of points where risk(i+1)>risk(i)
   Double_t trisk=0;                   // time per risk evaluation
   Double_t strisk=0;                  // total time
   Double_t rprev=1e32;                // previous risk

   // parameters set at point with min error
   Double_t              errmin=1e32;  // min error
   // Double_t              riskMin=0;    // risk
   Int_t                 indMin=-1;    // index
   std::vector<Double_t> coefsMin;     // rule coefs
   std::vector<Double_t> lincoefsMin;  // linear coefs
   Double_t              offsetMin;    // offset


   // DEBUG: timing
   clock_t  t0=0;
   Double_t tgradvec;
   Double_t tupgrade;
   Double_t tperf;
   Double_t stgradvec=0;
   Double_t stupgrade=0;
   Double_t stperf=0;

   // linear regression to estimate slope of error rate evolution
   const UInt_t npreg=5;
   std::vector<Double_t> valx;
   std::vector<Double_t> valy;
   std::vector<Double_t> valxy;

   // loop related
   Bool_t  docheck;          // true if an error rate check is to be done
   Int_t   iloop=0;          // loop index
   Bool_t  found=kFALSE;     // true if minimum is found
   Bool_t  riskFlat=kFALSE;  // DEBUG: flag is true if risk evolution behaves badly
   Bool_t  done = kFALSE;    // flag that the scan is done

   // calculate how often to check error rate
   int imod = fGDNPathSteps/100;
   if (imod<100) imod = std::min(100,fGDNPathSteps);
   if (imod>100) imod=100;

   // reset coefficients
   fAverageTruth = -CalcAverageTruth();
   offsetMin     = fAverageTruth;
   fRuleEnsemble->SetOffset(offsetMin);
   fRuleEnsemble->ClearCoefficients(0);
   fRuleEnsemble->ClearLinCoefficients(0);
   for (UInt_t i=0; i<fGDOfsTst.size(); i++) {
      fGDOfsTst[i] = offsetMin;
   }
   Log() << kVERBOSE << "Obtained initial offset = " << offsetMin << Endl;

   // find the best tau - returns the number of steps performed in scan
   Int_t nprescan = FindGDTau();

   //
   //
   // calculate F*
   //
   //   CalcFStar(fPerfIdx1, fPerfIdx2);
   //

   // set some ntuple values
   fNTRisk = rprev;
   fNTCoefRad = -1.0;
   fNTErrorRate = 0;

   // a local flag indicating for what reason the search was stopped
   Int_t stopCondition=0;

   Log() << kINFO << "Fitting model..." << Endl;
   // start loop with timer
   Timer timer( fGDNPathSteps, "RuleFit" );
   Log() << kWARNING;
   while (!done) {
      // Make gradient vector (eq 44, ref 1)
      if (isVerbose) t0 = clock();
      MakeGradientVector();
      if (isVerbose) {
         tgradvec = Double_t(clock()-t0)/CLOCKS_PER_SEC;
         stgradvec += tgradvec;
      }

      // Calculate the direction in parameter space (eq 25, ref 1) and update coeffs (eq 22, ref 1)
      if (isVerbose) t0 = clock();
      UpdateCoefficients();
      if (isVerbose) {
         tupgrade = Double_t(clock()-t0)/CLOCKS_PER_SEC;
         stupgrade += tupgrade;
      }

      // don't check error rate every loop
      docheck = ((iloop==0) ||((iloop+1)%imod==0));

      if (docheck) {
         // draw progressbar only if not debug
         if (!isVerbose)
            timer.DrawProgressBar(iloop);
         fNTNuval = Double_t(iloop)*fGDPathStep;
         fNTRisk = 0.0;

         // check risk evolution

         if (isDebug) FillCoefficients();
         fNTCoefRad = fRuleEnsemble->CoefficientRadius();

         // calculate risk
         t0 = clock();
         fNTRisk = RiskPath();
         trisk =  Double_t(clock()-t0)/CLOCKS_PER_SEC;
         strisk += trisk;
         //
         // Check for an increase in risk.
         // Such an increase would imply that the regularization is too small.
         // Stop the iteration if this happens.
         //
         if (fNTRisk>=rprev) {
            if (fNTRisk>rprev) {
               nbadrisk++;
               Log() << "Risk(i+1)>=Risk(i) in path" << Endl;
               riskFlat=(nbadrisk>3);
               if (riskFlat) {
                  Log() << kWARNING << "Chaotic behaviour of risk evolution" << Endl;
                  Log() << "--- STOPPING MINIMISATION ---" << Endl;
                  Log() << "This may be OK if minimum is already found" << Endl;
               }
            }
         }
         rprev = fNTRisk;
         //
         // Estimate the error rate using cross validation
         // Well, not quite full cross validation since we only
         // use ONE model.
         //
         if (isVerbose) t0 = clock();
         fNTErrorRate = 0;

         // Check error rate
         Double_t errroc;//= ErrorRateRoc();
         Double_t riskPerf = RiskPerf();
         //         Double_t optimism = Optimism();
         //
         errroc = riskPerf;
         //
         fNTErrorRate = errroc;
         //
         if (isVerbose) {
            tperf = Double_t(clock()-t0)/CLOCKS_PER_SEC;
            stperf +=tperf;
         }
         //
         // Always take the last min.
         // For each step the risk is reduced.
         //
         if (fNTErrorRate<=errmin) {
            errmin  = fNTErrorRate;
            // riskMin = fNTRisk;
            indMin  = iloop;
            fRuleEnsemble->GetCoefficients(coefsMin);
            lincoefsMin = fRuleEnsemble->GetLinCoefficients();
            offsetMin   = fRuleEnsemble->GetOffset();
         }
         if ( fNTErrorRate > fGDErrScale*errmin) found = kTRUE;
         //
         // check slope of last couple of points
         //
         if (valx.size()==npreg) {
            valx.erase(valx.begin());
            valy.erase(valy.begin());
            valxy.erase(valxy.begin());
         }
         valx.push_back(fNTNuval);
         valy.push_back(fNTErrorRate);
         valxy.push_back(fNTErrorRate*fNTNuval);

         gFIRSTORG=kFALSE;

         //
         if (isDebug) fGDNtuple->Fill();
         if (isVerbose) {
            Log() << kVERBOSE << "ParamsIRE : "
                  << std::setw(10)
                  << Form("%8d",iloop+1) << " "
                  << Form("%4.4f",fNTRisk) << " "
                  << Form("%4.4f",riskPerf)  << " "
                  << Form("%4.4f",fNTRisk+riskPerf)  << " "
               //                     << Form("%4.4f",fsigave+fbkgave) << " "
               //                     << Form("%4.4f",fsigave) << " "
               //                     << Form("%4.4f",fsigrms) << " "
               //                     << Form("%4.4f",fbkgave) << " "
               //                     << Form("%4.4f",fbkgrms) << " "

               //                    << Form("%4.4f",fRuleEnsemble->CoefficientRadius())
                  << Endl;
         }
      }
      iloop++;
      // Stop iteration under various conditions
      // * The condition R(i+1)<R(i) is no longer true (when then implicit regularization is too weak)
      // * If the current error estimate is > factor*errmin (factor = 1.1)
      // * We have reach the last step...
      Bool_t endOfLoop = (iloop==fGDNPathSteps);
      if ( ((riskFlat) || (endOfLoop)) && (!found) ) {
         if (riskFlat) {
            stopCondition = 1;
         }
         else if (endOfLoop) {
            stopCondition = 2;
         }
         if (indMin<0) {
            Log() << kWARNING << "BUG TRAP: should not be here - still, this bug is harmless;)" << Endl;
            errmin  = fNTErrorRate;
            // riskMin = fNTRisk;
            indMin  = iloop;
            fRuleEnsemble->GetCoefficients(coefsMin);
            lincoefsMin = fRuleEnsemble->GetLinCoefficients();
            offsetMin   = fRuleEnsemble->GetOffset();
         }
         found = kTRUE;
      }
      done = (found);
   }
   Log() << Endl;
   Log() << kINFO << "Minimisation elapsed time : " << timer.GetElapsedTime() << "                      " << Endl;
   Log() << kINFO << "----------------------------------------------------------------"  << Endl;
   Log() << kINFO << "Found minimum at step " << indMin+1 << " with error = " << errmin << Endl;
   Log() << kINFO << "Reason for ending loop: ";
   switch (stopCondition) {
   case 0:
      Log() << kINFO << "clear minima found";
      break;
   case 1:
      Log() << kINFO << "chaotic behaviour of risk";
      break;
   case 2:
      Log() << kINFO << "end of loop reached";
      break;
   default:
      Log() << kINFO << "unknown!";
      break;
   }
   Log() << Endl;
   Log() << kINFO << "----------------------------------------------------------------"  << Endl;

   // check if early minima - might be an indication of too large stepsize
   if ( Double_t(indMin)/Double_t(nprescan+fGDNPathSteps) < 0.05 ) {
      Log() << kWARNING << "Reached minimum early in the search" << Endl;
      Log() << "Check results and maybe decrease GDStep size" << Endl;
   }
   //
   // quick check of the sign of the slope for the last npreg points
   //
   Double_t sumx  = std::accumulate( valx.begin(), valx.end(), Double_t() );
   Double_t sumxy = std::accumulate( valxy.begin(), valxy.end(), Double_t() );
   Double_t sumy  = std::accumulate( valy.begin(), valy.end(), Double_t() );
   Double_t slope = Double_t(valx.size())*sumxy - sumx*sumy;
   if (slope<0) {
      Log() << kINFO << "The error rate was still decreasing at the end of the path" << Endl;
      Log() << kINFO << "Increase number of steps (GDNSteps)." << Endl;
   }
   //
   // set coefficients
   //
   if (found) {
      fRuleEnsemble->SetCoefficients( coefsMin );
      fRuleEnsemble->SetLinCoefficients( lincoefsMin );
      fRuleEnsemble->SetOffset( offsetMin );
   }
   else {
      Log() << kFATAL << "BUG TRAP: minimum not found in MakeGDPath()" << Endl;
   }

   //
   // print timing info (VERBOSE mode)
   //
   if (isVerbose) {
      Double_t stloop  = strisk +stupgrade + stgradvec + stperf;
      Log() << kVERBOSE << "Timing per loop (ms):" << Endl;
      Log() << kVERBOSE << "   gradvec = " << 1000*stgradvec/iloop << Endl;
      Log() << kVERBOSE << "   upgrade = " << 1000*stupgrade/iloop << Endl;
      Log() << kVERBOSE << "   risk    = " << 1000*strisk/iloop    << Endl;
      Log() << kVERBOSE << "   perf    = " << 1000*stperf/iloop    << Endl;
      Log() << kVERBOSE << "   loop    = " << 1000*stloop/iloop    << Endl;
      //
      Log() << kVERBOSE << "   GDInit      = " << 1000*gGDInit/iloop    << Endl;
      Log() << kVERBOSE << "   GDPtr       = " << 1000*gGDPtr/iloop    << Endl;
      Log() << kVERBOSE << "   GDEval      = " << 1000*gGDEval/iloop    << Endl;
      Log() << kVERBOSE << "   GDEvalRule  = " << 1000*gGDEvalRule/iloop    << Endl;
      Log() << kVERBOSE << "   GDNorm      = " << 1000*gGDNorm/iloop    << Endl;
      Log() << kVERBOSE << "   GDRuleLoop  = " << 1000*gGDRuleLoop/iloop    << Endl;
      Log() << kVERBOSE << "   GDLinLoop   = " << 1000*gGDLinLoop/iloop    << Endl;
   }
   //
   // write ntuple (DEBUG)
   if (isDebug) fGDNtuple->Write();
}

////////////////////////////////////////////////////////////////////////////////
/// helper function to store the rule coefficients in local arrays

void TMVA::RuleFitParams::FillCoefficients()
{
   fNTOffset = fRuleEnsemble->GetOffset();
   //
   for (UInt_t i=0; i<fNRules; i++) {
      fNTCoeff[i] = fRuleEnsemble->GetRules(i)->GetCoefficient();
   }
   for (UInt_t i=0; i<fNLinear; i++) {
      fNTLinCoeff[i] = fRuleEnsemble->GetLinCoefficients(i);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Estimates F* (optimum scoring function) for all events for the given sets.
/// The result is used in ErrorRateReg().
/// --- NOT USED ---

void TMVA::RuleFitParams::CalcFStar()
{
   Log() << kWARNING << "<CalcFStar> Using unverified code! Check!" << Endl;
   UInt_t neve = fPerfIdx2-fPerfIdx1+1;
   if (neve<1) {
      Log() << kFATAL << "<CalcFStar> Invalid start/end indices!" << Endl;
      return;
   }
   //
   const std::vector<const Event *> *events = &(fRuleFit->GetTrainingEvents());
   //
   fFstar.clear();
   std::vector<Double_t> fstarSorted;
   Double_t fstarVal;
   // loop over all events and estimate F* for each event
   for (UInt_t i=fPerfIdx1; i<fPerfIdx2+1; i++) {
      const Event& e = *(*events)[i];
      fstarVal = fRuleEnsemble->FStar(e);
      fFstar.push_back(fstarVal);
      fstarSorted.push_back(fstarVal);
      if (TMath::IsNaN(fstarVal)) Log() << kFATAL << "F* is NAN!" << Endl;
   }
   // sort F* and find median
   std::sort( fstarSorted.begin(), fstarSorted.end() );
   UInt_t ind = neve/2;
   if (neve&1) { // odd number of events
      fFstarMedian = 0.5*(fstarSorted[ind]+fstarSorted[ind-1]);
   }
   else { // even
      fFstarMedian = fstarSorted[ind];
   }
}

////////////////////////////////////////////////////////////////////////////////
/// implementation of eq. 7.17 in Hastie,Tibshirani & Friedman book
/// this is the covariance between the estimated response yhat and the
/// true value y.
/// NOT REALLY SURE IF THIS IS CORRECT!
/// --- THIS IS NOT USED ---

Double_t TMVA::RuleFitParams::Optimism()
{
   Log() << kWARNING << "<Optimism> Using unverified code! Check!" << Endl;
   UInt_t neve = fPerfIdx2-fPerfIdx1+1;
   if (neve<1) {
      Log() << kFATAL << "<Optimism> Invalid start/end indices!" << Endl;
   }
   //
   const std::vector<const Event *> *events = &(fRuleFit->GetTrainingEvents());
   //
   Double_t sumy=0;
   Double_t sumyhat=0;
   Double_t sumyhaty=0;
   Double_t sumw2=0;
   Double_t sumw=0;
   Double_t yhat;
   Double_t y;
   Double_t w;
   //
   for (UInt_t i=fPerfIdx1; i<fPerfIdx2+1; i++) {
      const Event& e = *(*events)[i];
      yhat = fRuleEnsemble->EvalEvent(i);         // evaluated using the model
      y    = (fRuleFit->GetMethodRuleFit()->DataInfo().IsSignal(&e) ? 1.0:-1.0);           // the truth
      w    = fRuleFit->GetTrainingEventWeight(i)/fNEveEffPerf; // the weight, reweighted such that sum=1
      sumy     += w*y;
      sumyhat  += w*yhat;
      sumyhaty += w*yhat*y;
      sumw2    += w*w;
      sumw     += w;
   }
   Double_t div = 1.0-sumw2;
   Double_t cov = sumyhaty - sumyhat*sumy;
   return 2.0*cov/div;
}

////////////////////////////////////////////////////////////////////////////////
/// Estimates the error rate with the current set of parameters
/// This code is pretty messy at the moment.
/// Cleanup is needed.
/// -- NOT USED ---

Double_t TMVA::RuleFitParams::ErrorRateReg()
{
   Log() << kWARNING << "<ErrorRateReg> Using unverified code! Check!" << Endl;
   UInt_t neve = fPerfIdx2-fPerfIdx1+1;
   if (neve<1) {
      Log() << kFATAL << "<ErrorRateReg> Invalid start/end indices!" << Endl;
   }
   if (fFstar.size()!=neve) {
      Log() << kFATAL << "--- RuleFitParams::ErrorRateReg() - F* not initialized! BUG!!!"
            << " Fstar.size() = " << fFstar.size() << " , N(events) = " << neve << Endl;
   }
   //
   Double_t sF;
   //
   const std::vector<const Event *> *events = &(fRuleFit->GetTrainingEvents());
   //
   Double_t sumdf = 0;
   Double_t sumdfmed = 0;
   //
   // A bit messy here.
   // I believe the binary error classification is appropriate here.
   // The problem is stability.
   //
   for (UInt_t i=fPerfIdx1; i<fPerfIdx2+1; i++) {
      const Event& e = *(*events)[i];
      sF = fRuleEnsemble->EvalEvent( e );
      // scaled abs error, eq 20 in RuleFit paper
      sumdf += TMath::Abs(fFstar[i-fPerfIdx1] - sF);
      sumdfmed += TMath::Abs(fFstar[i-fPerfIdx1] - fFstarMedian);
   }
   // scaled abs error, eq 20
   // This error (df) is large - need to think on how to compensate...
   //
   return sumdf/sumdfmed;
}

////////////////////////////////////////////////////////////////////////////////
///
/// Estimates the error rate with the current set of parameters
/// It uses a binary estimate of (y-F*(x))
/// (y-F*(x)) = (Num of events where sign(F)!=sign(y))/Neve
/// y = {+1 if event is signal, -1 otherwise}
/// --- NOT USED ---

Double_t TMVA::RuleFitParams::ErrorRateBin()
{
   Log() << kWARNING << "<ErrorRateBin> Using unverified code! Check!" << Endl;
   UInt_t neve = fPerfIdx2-fPerfIdx1+1;
   if (neve<1) {
      Log() << kFATAL << "<ErrorRateBin> Invalid start/end indices!" << Endl;
   }
   //
   const std::vector<const Event *> *events = &(fRuleFit->GetTrainingEvents());
   //
   Double_t sumdfbin = 0;
   Double_t dneve = Double_t(neve);
   Int_t signF, signy;
   Double_t sF;
   //
   for (UInt_t i=fPerfIdx1; i<fPerfIdx2+1; i++) {
      const Event& e = *(*events)[i];
      sF     = fRuleEnsemble->EvalEvent( e );
      //      Double_t sFstar = fRuleEnsemble->FStar(e); // THIS CAN BE CALCULATED ONCE!
      signF = (sF>0 ? +1:-1);
      //      signy = (sFStar>0 ? +1:-1);
      signy = (fRuleFit->GetMethodRuleFit()->DataInfo().IsSignal(&e) ? +1:-1);
      sumdfbin += TMath::Abs(Double_t(signF-signy))*0.5;
   }
   Double_t f = sumdfbin/dneve;
   //   Double_t   df = f*TMath::Sqrt((1.0/sumdfbin) + (1.0/dneve));
   return f;
}

////////////////////////////////////////////////////////////////////////////////
/// Estimates the error rate with the current set of parameters.
/// It calculates the area under the bkg rejection vs signal efficiency curve.
/// The value returned is 1-area.

Double_t TMVA::RuleFitParams::ErrorRateRocRaw( std::vector<Double_t> & sFsig,
                                               std::vector<Double_t> & sFbkg )

{
   std::sort(sFsig.begin(), sFsig.end());
   std::sort(sFbkg.begin(), sFbkg.end());
   const Double_t minsig = sFsig.front();
   const Double_t minbkg = sFbkg.front();
   const Double_t maxsig = sFsig.back();
   const Double_t maxbkg = sFbkg.back();
   const Double_t minf   = std::min(minsig,minbkg);
   const Double_t maxf   = std::max(maxsig,maxbkg);
   const Int_t    nsig   = Int_t(sFsig.size());
   const Int_t    nbkg   = Int_t(sFbkg.size());
   const Int_t    np     = std::min((nsig+nbkg)/4,50);
   const Double_t df     = (maxf-minf)/(np-1);
   //
   // calculate area under rejection/efficiency curve
   //
   Double_t fcut;
   std::vector<Double_t>::const_iterator indit;
   Int_t nrbkg;
   Int_t nesig;
   Int_t pnesig=0;
   Double_t rejb=0;
   Double_t effs=1.0;
   Double_t prejb=0;
   Double_t peffs=1.0;
   // Double_t drejb;
   Double_t deffs;
   Double_t area=0;
   Int_t    npok=0;
   //
   // loop over range of F [minf,maxf]
   //
   for (Int_t i=0; i<np; i++) {
      fcut = minf + df*Double_t(i);
      indit = std::find_if(sFsig.begin(), sFsig.end(),
                           std::bind(std::greater_equal<Double_t>(), std::placeholders::_1, fcut));
      nesig = sFsig.end()-indit; // number of sig accepted with F>cut
      if (TMath::Abs(pnesig-nesig)>0) {
         npok++;
         indit = std::find_if(sFbkg.begin(), sFbkg.end(),
                              std::bind(std::greater_equal<Double_t>(), std::placeholders::_1, fcut));
         nrbkg = indit-sFbkg.begin(); // number of bkg rejected with F>cut
         rejb = Double_t(nrbkg)/Double_t(nbkg);
         effs = Double_t(nesig)/Double_t(nsig);
         //
         // drejb = rejb-prejb;
         deffs = effs-peffs;
         area += 0.5*TMath::Abs(deffs)*(rejb+prejb); // trapezoid
         prejb = rejb;
         peffs = effs;
      }
      pnesig = nesig;
   }
   area += 0.5*(1+rejb)*effs; // extrapolate to the end point

   return (1.0-area);
}

////////////////////////////////////////////////////////////////////////////////
/// Estimates the error rate with the current set of parameters.
/// It calculates the area under the bkg rejection vs signal efficiency curve.
/// The value returned is 1-area.
/// This works but is less efficient than calculating the Risk using RiskPerf().

Double_t TMVA::RuleFitParams::ErrorRateRoc()
{
   Log() << kWARNING << "<ErrorRateRoc> Should not be used in the current version! Check!" << Endl;
   UInt_t neve = fPerfIdx2-fPerfIdx1+1;
   if (neve<1) {
      Log() << kFATAL << "<ErrorRateRoc> Invalid start/end indices!" << Endl;
   }
   //
   const std::vector<const Event *> *events = &(fRuleFit->GetTrainingEvents());
   //
   Double_t sF;
   //
   std::vector<Double_t> sFsig;
   std::vector<Double_t> sFbkg;
   Double_t sumfsig=0;
   Double_t sumfbkg=0;
   Double_t sumf2sig=0;
   Double_t sumf2bkg=0;
   //
   for (UInt_t i=fPerfIdx1; i<fPerfIdx2+1; i++) {
      const Event& e = *(*events)[i];
      sF = fRuleEnsemble->EvalEvent(i);// * fRuleFit->GetTrainingEventWeight(i);
      if (fRuleFit->GetMethodRuleFit()->DataInfo().IsSignal(&e)) {
         sFsig.push_back(sF);
         sumfsig  +=sF;
         sumf2sig +=sF*sF;
      }
      else {
         sFbkg.push_back(sF);
         sumfbkg  +=sF;
         sumf2bkg +=sF*sF;
      }
   }
   fsigave = sumfsig/sFsig.size();
   fbkgave = sumfbkg/sFbkg.size();
   fsigrms = TMath::Sqrt(gTools().ComputeVariance(sumf2sig,sumfsig,sFsig.size()));
   fbkgrms = TMath::Sqrt(gTools().ComputeVariance(sumf2bkg,sumfbkg,sFbkg.size()));
   //
   return ErrorRateRocRaw( sFsig, sFbkg );
}

////////////////////////////////////////////////////////////////////////////////
/// Estimates the error rate with the current set of parameters.
/// It calculates the area under the bkg rejection vs signal efficiency curve.
/// The value returned is 1-area.
///
/// See comment under ErrorRateRoc().

void TMVA::RuleFitParams::ErrorRateRocTst()
{
   Log() << kWARNING << "<ErrorRateRocTst> Should not be used in the current version! Check!" << Endl;
   UInt_t neve = fPerfIdx2-fPerfIdx1+1;
   if (neve<1) {
      Log() << kFATAL << "<ErrorRateRocTst> Invalid start/end indices!" << Endl;
      return;
   }
   //
   const std::vector<const Event *> *events = &(fRuleFit->GetTrainingEvents());
   //
   //   std::vector<Double_t> sF;
   Double_t sF;
   std::vector< std::vector<Double_t> > sFsig;
   std::vector< std::vector<Double_t> > sFbkg;
   //
   sFsig.resize( fGDNTau );
   sFbkg.resize( fGDNTau );
   //   sF.resize( fGDNTau );

   for (UInt_t i=fPerfIdx1; i<fPerfIdx2+1; i++) {
      for (UInt_t itau=0; itau<fGDNTau; itau++) {
         //         if (itau==0) sF = fRuleEnsemble->EvalEvent( *(*events)[i], fGDOfsTst[itau], fGDCoefTst[itau], fGDCoefLinTst[itau] );
         //         else         sF = fRuleEnsemble->EvalEvent(                fGDOfsTst[itau], fGDCoefTst[itau], fGDCoefLinTst[itau] );
         sF = fRuleEnsemble->EvalEvent( i, fGDOfsTst[itau], fGDCoefTst[itau], fGDCoefLinTst[itau] );
         if (fRuleFit->GetMethodRuleFit()->DataInfo().IsSignal((*events)[i])) {
            sFsig[itau].push_back(sF);
         }
         else {
            sFbkg[itau].push_back(sF);
         }
      }
   }
   Double_t err;

   for (UInt_t itau=0; itau<fGDNTau; itau++) {
      err = ErrorRateRocRaw( sFsig[itau], sFbkg[itau] );
      fGDErrTst[itau] = err;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Estimates the error rate with the current set of parameters.
/// using the <Perf> subsample.
/// Return the tau index giving the lowest error

UInt_t TMVA::RuleFitParams::RiskPerfTst()
{
   UInt_t neve = fPerfIdx2-fPerfIdx1+1;
   if (neve<1) {
      Log() << kFATAL << "<ErrorRateRocTst> Invalid start/end indices!" << Endl;
      return 0;
   }
   //
   Double_t sumx    = 0;
   Double_t sumx2   = 0;
   Double_t maxx    = -100.0;
   Double_t minx    = 1e30;
   UInt_t   itaumin = 0;
   UInt_t   nok=0;
   for (UInt_t itau=0; itau<fGDNTau; itau++) {
      if (fGDErrTstOK[itau]) {
         nok++;
         fGDErrTst[itau] = RiskPerf(itau);
         sumx  += fGDErrTst[itau];
         sumx2 += fGDErrTst[itau]*fGDErrTst[itau];
         if (fGDErrTst[itau]>maxx) maxx=fGDErrTst[itau];
         if (fGDErrTst[itau]<minx) {
            minx=fGDErrTst[itau];
            itaumin = itau;
         }
      }
   }
   Double_t sigx = TMath::Sqrt(gTools().ComputeVariance( sumx2, sumx, nok ) );
   Double_t maxacc = minx+sigx;
   //
   if (nok>0) {
      nok = 0;
      for (UInt_t itau=0; itau<fGDNTau; itau++) {
         if (fGDErrTstOK[itau]) {
            if (fGDErrTst[itau] > maxacc) {
               fGDErrTstOK[itau] = kFALSE;
            }
            else {
               nok++;
            }
         }
      }
   }
   fGDNTauTstOK = nok;
   Log() << kVERBOSE << "TAU: "
         << itaumin << "   "
         << nok     << "   "
         << minx    << "   "
         << maxx    << "   "
         << sigx    << Endl;
   //
   return itaumin;
}

////////////////////////////////////////////////////////////////////////////////
/// make test gradient vector for all tau
/// same algorithm as MakeGradientVector()

void TMVA::RuleFitParams::MakeTstGradientVector()
{
   UInt_t neve = fPathIdx1-fPathIdx2+1;
   if (neve<1) {
      Log() << kFATAL << "<MakeTstGradientVector> Invalid start/end indices!" << Endl;
      return;
   }
   //
   Double_t norm   = 2.0/fNEveEffPath;
   //
   const std::vector<const Event *> *events = &(fRuleFit->GetTrainingEvents());

   // Clear gradient vectors
   for (UInt_t itau=0; itau<fGDNTau; itau++) {
      if (fGDErrTstOK[itau]) {
         for (UInt_t ir=0; ir<fNRules; ir++) {
            fGradVecTst[itau][ir]=0;
         }
         for (UInt_t il=0; il<fNLinear; il++) {
            fGradVecLinTst[itau][il]=0;
         }
      }
   }
   //
   //   Double_t val; // temp store
   Double_t sF;   // score function value
   Double_t r;   // eq 35, ref 1
   Double_t y;   // true score (+1 or -1)
   const std::vector<UInt_t> *eventRuleMap=0;
   UInt_t rind;
   //
   // Loop over all events
   //
   UInt_t nsfok=0;
   for (UInt_t i=fPathIdx1; i<fPathIdx2+1; i++) {
      const Event *e = (*events)[i];
      UInt_t nrules=0;
      if (fRuleEnsemble->DoRules()) {
         eventRuleMap = &(fRuleEnsemble->GetEventRuleMap(i));
         nrules = (*eventRuleMap).size();
      }
      for (UInt_t itau=0; itau<fGDNTau; itau++) { // loop over tau
         //         if (itau==0) sF = fRuleEnsemble->EvalEvent( *e, fGDOfsTst[itau], fGDCoefTst[itau], fGDCoefLinTst[itau] );
         //         else         sF = fRuleEnsemble->EvalEvent(     fGDOfsTst[itau], fGDCoefTst[itau], fGDCoefLinTst[itau] );
         if (fGDErrTstOK[itau]) {
            sF = fRuleEnsemble->EvalEvent( i, fGDOfsTst[itau], fGDCoefTst[itau], fGDCoefLinTst[itau] );
            if (TMath::Abs(sF)<1.0) {
               nsfok++;
               r = 0;
               y = (fRuleFit->GetMethodRuleFit()->DataInfo().IsSignal(e)?1.0:-1.0);
               r = norm*(y - sF) * fRuleFit->GetTrainingEventWeight(i);
               // rule gradient vector
               for (UInt_t ir=0; ir<nrules; ir++) {
                  rind = (*eventRuleMap)[ir];
                  fGradVecTst[itau][rind] += r;
               }
               // linear terms
               for (UInt_t il=0; il<fNLinear; il++) {
                  fGradVecLinTst[itau][il] += r*fRuleEnsemble->EvalLinEventRaw( il,i, kTRUE );
               }
            } // if (TMath::Abs(F)<xxx)
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Establish maximum gradient for rules, linear terms and the offset
/// for all taus TODO: do not need index range!

void TMVA::RuleFitParams::UpdateTstCoefficients()
{
   Double_t maxr, maxl, cthresh, val;
   // loop over all taus
   for (UInt_t itau=0; itau<fGDNTau; itau++) {
      if (fGDErrTstOK[itau]) {
         // find max gradient
         maxr = ( (fNRules>0 ?
                   TMath::Abs(*(std::max_element( fGradVecTst[itau].begin(), fGradVecTst[itau].end(), AbsValue()))):0) );
         maxl = ( (fNLinear>0 ?
                   TMath::Abs(*(std::max_element( fGradVecLinTst[itau].begin(), fGradVecLinTst[itau].end(), AbsValue()))):0) );

         // Use the maximum as a threshold
         Double_t maxv = (maxr>maxl ? maxr:maxl);
         cthresh = maxv * fGDTauVec[itau];

         // Add to offset, if gradient is large enough:
         // Loop over the gradient vector and move to next set of coefficients
         // size of GradVec (and GradVecLin) should be 0 if learner is disabled
         //
         // Step-size is divided by 10 when looking for the best path.
         //
         if (maxv>0) {
            const Double_t stepScale=1.0;
            for (UInt_t i=0; i<fNRules; i++) {
               val = fGradVecTst[itau][i];

               if (TMath::Abs(val)>=cthresh) {
                  fGDCoefTst[itau][i] += fGDPathStep*val*stepScale;
               }
            }
            // Loop over the gradient vector for the linear part and move to next set of coefficients
            for (UInt_t i=0; i<fNLinear; i++) {
               val = fGradVecLinTst[itau][i];
               if (TMath::Abs(val)>=cthresh) {
                  fGDCoefLinTst[itau][i] += fGDPathStep*val*stepScale/fRuleEnsemble->GetLinNorm(i);
               }
            }
         }
      }
   }
   // set the offset - should be outside the itau loop!
   CalcTstAverageResponse();
}

////////////////////////////////////////////////////////////////////////////////
/// make gradient vector

void TMVA::RuleFitParams::MakeGradientVector()
{
   clock_t t0;
   //   clock_t t10;
   t0 = clock();
   //
   UInt_t neve = fPathIdx2-fPathIdx1+1;
   if (neve<1) {
      Log() << kFATAL << "<MakeGradientVector> Invalid start/end indices!" << Endl;
      return;
   }
   //
   const Double_t norm   = 2.0/fNEveEffPath;
   //
   const std::vector<const Event *> *events = &(fRuleFit->GetTrainingEvents());

   // Clear gradient vectors
   for (UInt_t ir=0; ir<fNRules; ir++) {
      fGradVec[ir]=0;
   }
   for (UInt_t il=0; il<fNLinear; il++) {
      fGradVecLin[il]=0;
   }
   //
   //   Double_t val; // temp store
   Double_t sF;   // score function value
   Double_t r;   // eq 35, ref 1
   Double_t y;   // true score (+1 or -1)
   const std::vector<UInt_t> *eventRuleMap=0;
   UInt_t rind;
   //
   gGDInit += Double_t(clock()-t0)/CLOCKS_PER_SEC;

   for (UInt_t i=fPathIdx1; i<fPathIdx2+1; i++) {
      const Event *e = (*events)[i];

      //    t0 = clock(); //DEB
      sF = fRuleEnsemble->EvalEvent( i ); // should not contain the weight
      //    gGDEval += Double_t(clock()-t0)/CLOCKS_PER_SEC;
      if (TMath::Abs(sF)<1.0) {
         UInt_t nrules=0;
         if (fRuleEnsemble->DoRules()) {
            eventRuleMap = &(fRuleEnsemble->GetEventRuleMap(i));
            nrules = (*eventRuleMap).size();
         }
         y = (fRuleFit->GetMethodRuleFit()->DataInfo().IsSignal(e)?1.0:-1.0);
         r = norm*(y - sF) * fRuleFit->GetTrainingEventWeight(i);
         // rule gradient vector
         for (UInt_t ir=0; ir<nrules; ir++) {
            rind = (*eventRuleMap)[ir];
            fGradVec[rind] += r;
         }
         //       gGDRuleLoop += Double_t(clock()-t0)/CLOCKS_PER_SEC;
         // linear terms
         //       t0 = clock(); //DEB
         for (UInt_t il=0; il<fNLinear; il++) {
            fGradVecLin[il] += r*fRuleEnsemble->EvalLinEventRaw( il, i, kTRUE );
         }
         //       gGDLinLoop += Double_t(clock()-t0)/CLOCKS_PER_SEC;
      } // if (TMath::Abs(F)<xxx)
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Establish maximum gradient for rules, linear terms and the offset

void TMVA::RuleFitParams::UpdateCoefficients()
{
   Double_t maxr = ( (fRuleEnsemble->DoRules() ?
                      TMath::Abs(*(std::max_element( fGradVec.begin(), fGradVec.end(), AbsValue()))):0) );
   Double_t maxl = ( (fRuleEnsemble->DoLinear() ?
                      TMath::Abs(*(std::max_element( fGradVecLin.begin(), fGradVecLin.end(), AbsValue()))):0) );
   // Use the maximum as a threshold
   Double_t maxv = (maxr>maxl ? maxr:maxl);
   Double_t cthresh = maxv * fGDTau;

   Double_t useRThresh;
   Double_t useLThresh;
   //
   // Choose thresholds.
   //
   useRThresh = cthresh;
   useLThresh = cthresh;

   Double_t gval, lval, coef, lcoef;

   // Add to offset, if gradient is large enough:
   // Loop over the gradient vector and move to next set of coefficients
   // size of GradVec (and GradVecLin) should be 0 if learner is disabled
   if (maxv>0) {
      for (UInt_t i=0; i<fGradVec.size(); i++) {
         gval = fGradVec[i];
         if (TMath::Abs(gval)>=useRThresh) {
            coef = fRuleEnsemble->GetRulesConst(i)->GetCoefficient() + fGDPathStep*gval;
            fRuleEnsemble->GetRules(i)->SetCoefficient(coef);
         }
      }

      // Loop over the gradient vector for the linear part and move to next set of coefficients
      for (UInt_t i=0; i<fGradVecLin.size(); i++) {
         lval = fGradVecLin[i];
         if (TMath::Abs(lval)>=useLThresh) {
            lcoef = fRuleEnsemble->GetLinCoefficients(i) + (fGDPathStep*lval/fRuleEnsemble->GetLinNorm(i));
            fRuleEnsemble->SetLinCoefficient(i,lcoef);
         }
      }
      // Set the offset
      Double_t offset = CalcAverageResponse();
      fRuleEnsemble->SetOffset( offset );
   }
}

////////////////////////////////////////////////////////////////////////////////
/// calc average response for all test paths - TODO: see comment under CalcAverageResponse()
/// note that 0 offset is used

void TMVA::RuleFitParams::CalcTstAverageResponse()
{
   for (UInt_t itau=0; itau<fGDNTau; itau++) {
      if (fGDErrTstOK[itau]) {
         fGDOfsTst[itau] = 0;
         for (UInt_t s=0; s<fNLinear; s++) {
            fGDOfsTst[itau] -= fGDCoefLinTst[itau][s] * fAverageSelectorPath[s];
         }
         for (UInt_t r=0; r<fNRules; r++) {
            fGDOfsTst[itau] -= fGDCoefTst[itau][r] * fAverageRulePath[r];
         }
      }
   }
   //
}

////////////////////////////////////////////////////////////////////////////////
/// calculate the average response - TODO : rewrite bad dependancy on EvaluateAverage() !
///
/// note that 0 offset is used

Double_t TMVA::RuleFitParams::CalcAverageResponse()
{
   Double_t ofs = 0;
   for (UInt_t s=0; s<fNLinear; s++) {
      ofs -= fRuleEnsemble->GetLinCoefficients(s) * fAverageSelectorPath[s];
   }
   for (UInt_t r=0; r<fNRules; r++) {
      ofs -= fRuleEnsemble->GetRules(r)->GetCoefficient() * fAverageRulePath[r];
   }
   return ofs;
}

////////////////////////////////////////////////////////////////////////////////
/// calculate the average truth

Double_t TMVA::RuleFitParams::CalcAverageTruth()
{
   if (fPathIdx2<=fPathIdx1) {
      Log() << kFATAL << "<CalcAverageTruth> Invalid start/end indices!" << Endl;
      return 0;
   }
   Double_t sum=0;
   Double_t ensig=0;
   Double_t enbkg=0;
   const std::vector<const Event *> *events = &(fRuleFit->GetTrainingEvents());
   for (UInt_t i=fPathIdx1; i<fPathIdx2+1; i++) {
      Double_t ew = fRuleFit->GetTrainingEventWeight(i);
      if (fRuleFit->GetMethodRuleFit()->DataInfo().IsSignal((*events)[i])) ensig += ew;
      else                          enbkg += ew;
      sum += ew*(fRuleFit->GetMethodRuleFit()->DataInfo().IsSignal((*events)[i])?1.0:-1.0);
   }
   Log() << kVERBOSE << "Effective number of signal / background = " << ensig << " / " << enbkg << Endl;

   return sum/fNEveEffPath;
}

////////////////////////////////////////////////////////////////////////////////

Int_t  TMVA::RuleFitParams::Type( const Event * e ) const {
   return (fRuleFit->GetMethodRuleFit()->DataInfo().IsSignal(e) ? 1:-1);
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::RuleFitParams::SetMsgType( EMsgType t ) {
   fLogger->SetMinType(t);
}
