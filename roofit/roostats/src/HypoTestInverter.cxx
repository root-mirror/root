// @(#)root/roostats:$Id$
// Author: Kyle Cranmer, Lorenzo Moneta, Gregory Schott, Wouter Verkerke
// Contributions: Giovanni Petrucciani and Annapaola Decosa
/*************************************************************************
 * Copyright (C) 1995-2008, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class RooStats::HypoTestInverter
    \ingroup Roostats

A class for performing a hypothesis test inversion by scanning
the hypothesis test results of a HypoTestCalculator  for various values of the
parameter of interest. By looking at the confidence level curve of the result, an
upper limit can be derived by computing the intersection of the confidence level curve with the desired confidence level.
The class implements the RooStats::IntervalCalculator interface, and returns a
RooStats::HypoTestInverterResult. The result is a SimpleInterval, which
via the method UpperLimit() returns to the user the upper limit value.

## Scanning options
The  HypoTestInverter implements various options for performing the scan.
- HypoTestInverter::RunFixedScan will scan the parameter of interest using a fixed grid.
- HypoTestInverter::SetAutoScan will perform an automatic scan to find
optimally the curve. It will stop when the desired precision is obtained.
- HypoTestInverter::RunOnePoint computes the confidence level at a given point.

### CLs presciption
The class can scan the CLs+b values or alternatively CLs. For the latter,
call HypoTestInverter::UseCLs().
*/

// include other header files

#include "RooAbsData.h"
#
#include "TMath.h"

#include "RooStats/HybridResult.h"

#include "RooStats/HypoTestInverter.h"
#include <cassert>
#include <cmath>

#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooAbsPdf.h"
#include "RooRandom.h"
#include "RooConstVar.h"
#include "RooMsgService.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/HybridCalculator.h"
#include "RooStats/FrequentistCalculator.h"
#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/SimpleLikelihoodRatioTestStat.h"
#include "RooStats/RatioOfProfiledLikelihoodsTestStat.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/HypoTestPlot.h"
#include "RooStats/HypoTestInverterPlot.h"

#include "RooStats/ProofConfig.h"

ClassImp(RooStats::HypoTestInverter);

using namespace RooStats;
using namespace std;

// static variable definitions
double HypoTestInverter::fgCLAccuracy = 0.005;
unsigned int HypoTestInverter::fgNToys = 500;

double HypoTestInverter::fgAbsAccuracy = 0.05;
double HypoTestInverter::fgRelAccuracy = 0.05;
std::string HypoTestInverter::fgAlgo = "logSecant";

bool HypoTestInverter::fgCloseProof = false;

// helper class to wrap the functionality of the various HypoTestCalculators

template<class HypoTestType>
struct HypoTestWrapper {

   static void SetToys(HypoTestType * h, int toyNull, int toyAlt) { h->SetToys(toyNull,toyAlt); }

};

////////////////////////////////////////////////////////////////////////////////
/// set flag to close proof for every new run

void HypoTestInverter::SetCloseProof(Bool_t flag) {
   fgCloseProof  = flag;
}

////////////////////////////////////////////////////////////////////////////////
/// get  the variable to scan
/// try first with null model if not go to alternate model

RooRealVar * HypoTestInverter::GetVariableToScan(const HypoTestCalculatorGeneric &hc) {

   RooRealVar * varToScan = 0;
   const ModelConfig * mc = hc.GetNullModel();
   if (mc) {
      const RooArgSet * poi  = mc->GetParametersOfInterest();
      if (poi) varToScan = dynamic_cast<RooRealVar*> (poi->first() );
   }
   if (!varToScan) {
      mc = hc.GetAlternateModel();
      if (mc) {
         const RooArgSet * poi  = mc->GetParametersOfInterest();
         if (poi) varToScan = dynamic_cast<RooRealVar*> (poi->first() );
      }
   }
   return varToScan;
}

////////////////////////////////////////////////////////////////////////////////
/// check  the model given the given hypotestcalculator

void HypoTestInverter::CheckInputModels(const HypoTestCalculatorGeneric &hc,const RooRealVar & scanVariable) {
   const ModelConfig * modelSB = hc.GetNullModel();
   const ModelConfig * modelB = hc.GetAlternateModel();
   if (!modelSB || ! modelB)
      oocoutF((TObject*)0,InputArguments) << "HypoTestInverter - model are not existing" << std::endl;
   assert(modelSB && modelB);

   oocoutI((TObject*)0,InputArguments) << "HypoTestInverter ---- Input models: \n"
                                       << "\t\t using as S+B (null) model     : "
                                       << modelSB->GetName() << "\n"
                                       << "\t\t using as B (alternate) model  : "
                                       << modelB->GetName() << "\n" << std::endl;

   // check if scanVariable is included in B model pdf
   RooAbsPdf * bPdf = modelB->GetPdf();
   const RooArgSet * bObs = modelB->GetObservables();
   if (!bPdf || !bObs) {
      oocoutE((TObject*)0,InputArguments) << "HypoTestInverter - B model has no pdf or observables defined" <<  std::endl;
      return;
   }
   RooArgSet * bParams = bPdf->getParameters(*bObs);
   if (!bParams) {
      oocoutE((TObject*)0,InputArguments) << "HypoTestInverter - pdf of B model has no parameters" << std::endl;
      return;
   }
   if (bParams->find(scanVariable.GetName() ) ) {
      const RooArgSet * poiB  = modelB->GetSnapshot();
      if (!poiB ||  !poiB->find(scanVariable.GetName()) ||
          ( (RooRealVar*)  poiB->find(scanVariable.GetName()) )->getVal() != 0 )
         oocoutW((TObject*)0,InputArguments) << "HypoTestInverter - using a B model  with POI "
                                             <<    scanVariable.GetName()  << " not equal to zero "
                                             << " user must check input model configurations " << endl;
      if (poiB) delete poiB;
   }
   delete bParams;
}

////////////////////////////////////////////////////////////////////////////////
/// default constructor (doesn't do anything)

HypoTestInverter::HypoTestInverter( ) :
   fTotalToysRun(0),
   fMaxToys(0),
   fCalculator0(0),
   fScannedVariable(0),
   fResults(0),
   fUseCLs(false),
   fScanLog(false),
   fSize(0),
   fVerbose(0),
   fCalcType(kUndefined),
   fNBins(0), fXmin(1), fXmax(1),
   fNumErr(0)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor from a HypoTestCalculatorGeneric
/// The HypoTest calculator must be a FrequentistCalculator or HybridCalculator type
/// Other type of calculators are not supported.
/// The calculator must be created before by using the S+B model for the null and
/// the B model for the alt
/// If no variable to scan are given they are assumed to be the first variable
/// from the parameter of interests of the null model

HypoTestInverter::HypoTestInverter( HypoTestCalculatorGeneric& hc,
                                    RooRealVar* scannedVariable, double size ) :
   fTotalToysRun(0),
   fMaxToys(0),
   fCalculator0(0),
   fScannedVariable(scannedVariable),
   fResults(0),
   fUseCLs(false),
   fScanLog(false),
   fSize(size),
   fVerbose(0),
   fCalcType(kUndefined),
   fNBins(0), fXmin(1), fXmax(1),
   fNumErr(0)
{

   if (!fScannedVariable) {
      fScannedVariable = HypoTestInverter::GetVariableToScan(hc);
   }
   if (!fScannedVariable)
      oocoutE((TObject*)0,InputArguments) << "HypoTestInverter - Cannot guess the variable to scan " << std::endl;
   else
      CheckInputModels(hc,*fScannedVariable);

   HybridCalculator * hybCalc = dynamic_cast<HybridCalculator*>(&hc);
   if (hybCalc) {
      fCalcType = kHybrid;
      fCalculator0 = hybCalc;
      return;
   }
   FrequentistCalculator * freqCalc = dynamic_cast<FrequentistCalculator*>(&hc);
   if (freqCalc) {
      fCalcType = kFrequentist;
      fCalculator0 = freqCalc;
      return;
   }
   AsymptoticCalculator * asymCalc = dynamic_cast<AsymptoticCalculator*>(&hc);
   if (asymCalc) {
      fCalcType = kAsymptotic;
      fCalculator0 = asymCalc;
      return;
   }
   oocoutE((TObject*)0,InputArguments) << "HypoTestInverter - Type of hypotest calculator is not supported " <<std::endl;
   fCalculator0 = &hc;
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor from a reference to a HybridCalculator
/// The calculator must be created before by using the S+B model for the null and
/// the B model for the alt
/// If no variable to scan are given they are assumed to be the first variable
/// from the parameter of interests of the null model

HypoTestInverter::HypoTestInverter( HybridCalculator& hc,
                                          RooRealVar* scannedVariable, double size ) :
   fTotalToysRun(0),
   fMaxToys(0),
   fCalculator0(&hc),
   fScannedVariable(scannedVariable),
   fResults(0),
   fUseCLs(false),
   fScanLog(false),
   fSize(size),
   fVerbose(0),
   fCalcType(kHybrid),
   fNBins(0), fXmin(1), fXmax(1),
   fNumErr(0)
{

   if (!fScannedVariable) {
      fScannedVariable = GetVariableToScan(hc);
   }
   if (!fScannedVariable)
      oocoutE((TObject*)0,InputArguments) << "HypoTestInverter - Cannot guess the variable to scan " << std::endl;
   else
      CheckInputModels(hc,*fScannedVariable);

}

////////////////////////////////////////////////////////////////////////////////
/// Constructor from a reference to a FrequentistCalculator
/// The calculator must be created before by using the S+B model for the null and
/// the B model for the alt
/// If no variable to scan are given they are assumed to be the first variable
/// from the parameter of interests of the null model

HypoTestInverter::HypoTestInverter( FrequentistCalculator& hc,
                                    RooRealVar* scannedVariable, double size ) :
   fTotalToysRun(0),
   fMaxToys(0),
   fCalculator0(&hc),
   fScannedVariable(scannedVariable),
   fResults(0),
   fUseCLs(false),
   fScanLog(false),
   fSize(size),
   fVerbose(0),
   fCalcType(kFrequentist),
   fNBins(0), fXmin(1), fXmax(1),
   fNumErr(0)
{

   if (!fScannedVariable) {
      fScannedVariable = GetVariableToScan(hc);
   }
   if (!fScannedVariable)
      oocoutE((TObject*)0,InputArguments) << "HypoTestInverter - Cannot guess the variable to scan " << std::endl;
   else
      CheckInputModels(hc,*fScannedVariable);
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor from a reference to a AsymptoticCalculator
/// The calculator must be created before by using the S+B model for the null and
/// the B model for the alt
/// If no variable to scan are given they are assumed to be the first variable
/// from the parameter of interests of the null model

HypoTestInverter::HypoTestInverter( AsymptoticCalculator& hc,
                                          RooRealVar* scannedVariable, double size ) :
   fTotalToysRun(0),
   fMaxToys(0),
   fCalculator0(&hc),
   fScannedVariable(scannedVariable),
   fResults(0),
   fUseCLs(false),
   fScanLog(false),
   fSize(size),
   fVerbose(0),
   fCalcType(kAsymptotic),
   fNBins(0), fXmin(1), fXmax(1),
   fNumErr(0)
{

   if (!fScannedVariable) {
      fScannedVariable = GetVariableToScan(hc);
   }
   if (!fScannedVariable)
      oocoutE((TObject*)0,InputArguments) << "HypoTestInverter - Cannot guess the variable to scan " << std::endl;
   else
      CheckInputModels(hc,*fScannedVariable);

}

////////////////////////////////////////////////////////////////////////////////
/// Constructor from a model for B model and a model for S+B.
/// An HypoTestCalculator (Hybrid of Frequentis) will be created using the
/// S+B model as the null and the B model as the alternate
/// If no variable to scan are given they are assumed to be the first variable
/// from the parameter of interests of the null model

HypoTestInverter::HypoTestInverter( RooAbsData& data, ModelConfig &sbModel, ModelConfig &bModel,
                RooRealVar * scannedVariable,  ECalculatorType type, double size) :
   fTotalToysRun(0),
   fMaxToys(0),
   fCalculator0(0),
   fScannedVariable(scannedVariable),
   fResults(0),
   fUseCLs(false),
   fScanLog(false),
   fSize(size),
   fVerbose(0),
   fCalcType(type),
   fNBins(0), fXmin(1), fXmax(1),
   fNumErr(0)
{
   if(fCalcType==kFrequentist) fHC.reset(new FrequentistCalculator(data, bModel, sbModel));
   if(fCalcType==kHybrid) fHC.reset( new HybridCalculator(data, bModel, sbModel)) ;
   if(fCalcType==kAsymptotic) fHC.reset( new AsymptoticCalculator(data, bModel, sbModel));
   fCalculator0 = fHC.get();
   // get scanned variable
   if (!fScannedVariable) {
      fScannedVariable = GetVariableToScan(*fCalculator0);
   }
   if (!fScannedVariable)
      oocoutE((TObject*)0,InputArguments) << "HypoTestInverter - Cannot guess the variable to scan " << std::endl;
   else
      CheckInputModels(*fCalculator0,*fScannedVariable);

}

////////////////////////////////////////////////////////////////////////////////
/// copy-constructor
/// NOTE: this class does not copy the contained result and
/// the HypoTestCalculator, but only the pointers
/// It requires the original HTI to be alive

HypoTestInverter::HypoTestInverter(const HypoTestInverter & rhs) :
   IntervalCalculator(),
   fTotalToysRun(0),
   fCalculator0(0), fScannedVariable(0),  // add these for Coverity
   fResults(0)
{
   (*this) = rhs;
}

////////////////////////////////////////////////////////////////////////////////
/// assignment operator
/// NOTE: this class does not copy the contained result and
/// the HypoTestCalculator, but only the pointers
/// It requires the original HTI to be alive

HypoTestInverter & HypoTestInverter::operator= (const HypoTestInverter & rhs) {
   if (this == &rhs) return *this;
   fTotalToysRun = 0;
   fMaxToys = rhs.fMaxToys;
   fCalculator0 = rhs.fCalculator0;
   fScannedVariable = rhs.fScannedVariable;
   fUseCLs = rhs.fUseCLs;
   fScanLog = rhs.fScanLog;
   fSize = rhs.fSize;
   fVerbose = rhs.fVerbose;
   fCalcType = rhs.fCalcType;
   fNBins = rhs.fNBins;
   fXmin = rhs.fXmin;
   fXmax = rhs.fXmax;
   fNumErr = rhs.fNumErr;

   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// destructor (delete the HypoTestInverterResult)

HypoTestInverter::~HypoTestInverter()
{
   if (fResults) delete fResults;
   fCalculator0 = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// return the test statistic which is or will be used by the class

TestStatistic * HypoTestInverter::GetTestStatistic( ) const
{
   if(fCalculator0 &&  fCalculator0->GetTestStatSampler()){
      return fCalculator0->GetTestStatSampler()->GetTestStatistic();
   }
   else
      return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// set the test statistic to use

bool HypoTestInverter::SetTestStatistic(TestStatistic& stat)
{
   if(fCalculator0 &&  fCalculator0->GetTestStatSampler()){
      fCalculator0->GetTestStatSampler()->SetTestStatistic(&stat);
      return true;
   }
   else return false;
}

////////////////////////////////////////////////////////////////////////////////
/// delete contained result and graph

void  HypoTestInverter::Clear()  {
   if (fResults) delete fResults;
   fResults = 0;
   fLimitPlot.reset(nullptr);
}

////////////////////////////////////////////////////////////////////////////////
/// create a new HypoTestInverterResult to hold all computed results

void  HypoTestInverter::CreateResults() const {
   if (fResults == 0) {
      TString results_name = "result_";
      results_name += fScannedVariable->GetName();
      fResults = new HypoTestInverterResult(results_name,*fScannedVariable,ConfidenceLevel());
      TString title = "HypoTestInverter Result For ";
      title += fScannedVariable->GetName();
      fResults->SetTitle(title);
   }
   fResults->UseCLs(fUseCLs);
   fResults->SetConfidenceLevel(1.-fSize);
   // check if one or two sided scan
   if (fCalculator0) {
      // if asymptotic calculator
      AsymptoticCalculator * ac = dynamic_cast<AsymptoticCalculator*>(fCalculator0);
      if (ac)
         fResults->fIsTwoSided = ac->IsTwoSided();
      else {
         // in case of the other calculators
         TestStatSampler * sampler = fCalculator0->GetTestStatSampler();
         if (sampler) {
            ProfileLikelihoodTestStat * pl = dynamic_cast<ProfileLikelihoodTestStat*>(sampler->GetTestStatistic());
            if (pl && pl->IsTwoSided() ) fResults->fIsTwoSided = true;
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Run a fixed scan or the automatic scan depending on the configuration.
/// Return if needed a copy of the result object which will be managed by the user.

HypoTestInverterResult* HypoTestInverter::GetInterval() const {

   // if having a result with at least  one point return it
   if (fResults && fResults->ArraySize() >= 1) {
      oocoutI((TObject*)0,Eval) << "HypoTestInverter::GetInterval - return an already existing interval " << std::endl;
      return  (HypoTestInverterResult*)(fResults->Clone());
   }

   if (fNBins > 0) {
      oocoutI((TObject*)0,Eval) << "HypoTestInverter::GetInterval - run a fixed scan" << std::endl;
      bool ret = RunFixedScan(fNBins, fXmin, fXmax, fScanLog);
      if (!ret)
         oocoutE((TObject*)0,Eval) << "HypoTestInverter::GetInterval - error running a fixed scan " << std::endl;
   }
   else {
      oocoutI((TObject*)0,Eval) << "HypoTestInverter::GetInterval - run an automatic scan" << std::endl;
      double limit(0),err(0);
      bool ret = RunLimit(limit,err);
      if (!ret)
         oocoutE((TObject*)0,Eval) << "HypoTestInverter::GetInterval - error running an auto scan " << std::endl;
   }

   if (fgCloseProof) ProofConfig::CloseProof();

   return (HypoTestInverterResult*) (fResults->Clone());
}

////////////////////////////////////////////////////////////////////////////////
/// Run the Hypothesis test at a previous configured point
/// (internal function called by RunOnePoint)

HypoTestResult * HypoTestInverter::Eval(HypoTestCalculatorGeneric &hc, bool adaptive, double clsTarget) const {
   //for debug
   // std::cout << ">>>>>>>>>>> " << std::endl;
   // std::cout << "alternate model " << std::endl;
   // hc.GetAlternateModel()->GetNuisanceParameters()->Print("V");
   // hc.GetAlternateModel()->GetParametersOfInterest()->Print("V");
   // std::cout << "Null model " << std::endl;
   // hc.GetNullModel()->GetNuisanceParameters()->Print("V");
   // hc.GetNullModel()->GetParametersOfInterest()->Print("V");
   // std::cout << "<<<<<<<<<<<<<<< " << std::endl;

   // run the hypothesis test
   HypoTestResult *  hcResult = hc.GetHypoTest();
   if (hcResult == 0) {
      oocoutE((TObject*)0,Eval) << "HypoTestInverter::Eval - HypoTest failed" << std::endl;
      return hcResult;
   }

   // since the b model is the alt need to set the flag
   hcResult->SetBackgroundAsAlt(true);


   // bool flipPvalue = false;
   // if (flipPValues)
   //       hcResult->SetPValueIsRightTail(!hcResult->GetPValueIsRightTail());

   // adjust for some numerical error in discrete models and == is not anymore
   if (hcResult->GetPValueIsRightTail() )
      hcResult->SetTestStatisticData(hcResult->GetTestStatisticData()-fNumErr); // issue with < vs <= in discrete models
   else
      hcResult->SetTestStatisticData(hcResult->GetTestStatisticData()+fNumErr); // issue with < vs <= in discrete models

   double clsMid    = (fUseCLs ? hcResult->CLs()      : hcResult->CLsplusb());
   double clsMidErr = (fUseCLs ? hcResult->CLsError() : hcResult->CLsplusbError());

   //if (fVerbose) std::cout << (fUseCLs ? "\tCLs = " : "\tCLsplusb = ") << clsMid << " +/- " << clsMidErr << std::endl;

   if (adaptive) {

      if (fCalcType == kHybrid) HypoTestWrapper<HybridCalculator>::SetToys((HybridCalculator*)&hc, fUseCLs ? fgNToys : 1, 4*fgNToys);
      if (fCalcType == kFrequentist) HypoTestWrapper<FrequentistCalculator>::SetToys((FrequentistCalculator*)&hc, fUseCLs ? fgNToys : 1, 4*fgNToys);

   while (clsMidErr >= fgCLAccuracy && (clsTarget == -1 || fabs(clsMid-clsTarget) < 3*clsMidErr) ) {
      std::unique_ptr<HypoTestResult> more(hc.GetHypoTest());

      // if (flipPValues)
      //    more->SetPValueIsRightTail(!more->GetPValueIsRightTail());

      hcResult->Append(more.get());
      clsMid    = (fUseCLs ? hcResult->CLs()      : hcResult->CLsplusb());
      clsMidErr = (fUseCLs ? hcResult->CLsError() : hcResult->CLsplusbError());
      if (fVerbose) std::cout << (fUseCLs ? "\tCLs = " : "\tCLsplusb = ") << clsMid << " +/- " << clsMidErr << std::endl;
   }

   }
   if (fVerbose ) {
      oocoutP((TObject*)0,Eval) << "P values for  " << fScannedVariable->GetName()  << " =  " <<
         fScannedVariable->getVal() << "\n" <<
         "\tCLs      = " << hcResult->CLs()      << " +/- " << hcResult->CLsError()      << "\n" <<
         "\tCLb      = " << hcResult->CLb()      << " +/- " << hcResult->CLbError()      << "\n" <<
         "\tCLsplusb = " << hcResult->CLsplusb() << " +/- " << hcResult->CLsplusbError() << "\n" <<
         std::endl;
   }

   if (fCalcType == kFrequentist || fCalcType == kHybrid)  {
      fTotalToysRun += (hcResult->GetAltDistribution()->GetSize() + hcResult->GetNullDistribution()->GetSize());

      // set sampling distribution name
      TString nullDistName = TString::Format("%s_%s_%4.2f",hcResult->GetNullDistribution()->GetName(),
                                             fScannedVariable->GetName(), fScannedVariable->getVal() );
      TString altDistName = TString::Format("%s_%s_%4.2f",hcResult->GetAltDistribution()->GetName(),
                                            fScannedVariable->GetName(), fScannedVariable->getVal() );

      hcResult->GetNullDistribution()->SetName(nullDistName);
      hcResult->GetAltDistribution()->SetName(altDistName);
   }

   return hcResult;
}

////////////////////////////////////////////////////////////////////////////////
/// Run a Fixed scan in npoints between min and max

bool HypoTestInverter::RunFixedScan( int nBins, double xMin, double xMax, bool scanLog ) const
{

   CreateResults();
   // interpolate the limits
   fResults->fFittedLowerLimit = false;
   fResults->fFittedUpperLimit = false;

   // safety checks
   if ( nBins<=0 ) {
      oocoutE((TObject*)0,InputArguments) << "HypoTestInverter::RunFixedScan - Please provide nBins>0\n";
      return false;
   }
   if ( nBins==1 && xMin!=xMax ) {
      oocoutW((TObject*)0,InputArguments) << "HypoTestInverter::RunFixedScan - nBins==1 -> I will run for xMin (" << xMin << ")\n";
   }
   if ( xMin==xMax && nBins>1 ) {
      oocoutW((TObject*)0,InputArguments) << "HypoTestInverter::RunFixedScan - xMin==xMax -> I will enforce nBins==1\n";
      nBins = 1;
   }
   if ( xMin>xMax ) {
      oocoutE((TObject*)0,InputArguments) << "HypoTestInverter::RunFixedScan - Please provide xMin ("
                                          << xMin << ") smaller than xMax (" << xMax << ")\n";
      return false;
   }

   if (xMin < fScannedVariable->getMin()) {
      xMin = fScannedVariable->getMin();
      oocoutW((TObject*)0,InputArguments) << "HypoTestInverter::RunFixedScan - xMin < lower bound, using xmin = "
                                          << xMin << std::endl;
   }
   if (xMax > fScannedVariable->getMax()) {
      xMax = fScannedVariable->getMax();
      oocoutW((TObject*)0,InputArguments) << "HypoTestInverter::RunFixedScan - xMax > upper bound, using xmax = "
                                          << xMax << std::endl;
   }

   if (xMin <= 0. && scanLog) {
     oocoutE((TObject*)nullptr, InputArguments) << "HypoTestInverter::RunFixedScan - cannot go in log steps if xMin <= 0" << std::endl;
     return false;
   }

   double thisX = xMin;
   for (int i=0; i<nBins; i++) {

      if (i > 0) { // avoids case of nBins = 1
         if (scanLog)
            thisX = exp(  log(xMin) +  i*(log(xMax)-log(xMin))/(nBins-1)  );  // scan in log x
         else
            thisX = xMin + i*(xMax-xMin)/(nBins-1);          // linear scan in x
      }

      const bool status = RunOnePoint(thisX);

      // check if failed status
      if ( status==false ) {
        oocoutW((TObject*)0,Eval) << "HypoTestInverter::RunFixedScan - The hypo test for point " << thisX << " failed. Skipping." << std::endl;
      }
   }

   return true;
}

////////////////////////////////////////////////////////////////////////////////
/// run only one point at the given POI value

bool HypoTestInverter::RunOnePoint( double rVal, bool adaptive, double clTarget) const
{

   CreateResults();

   // check if rVal is in the range specified for fScannedVariable
   if ( rVal < fScannedVariable->getMin() ) {
      oocoutE((TObject*)0,InputArguments) << "HypoTestInverter::RunOnePoint - Out of range: using the lower bound "
                                          << fScannedVariable->getMin()
                                          << " on the scanned variable rather than " << rVal<< "\n";
     rVal = fScannedVariable->getMin();
   }
   if ( rVal > fScannedVariable->getMax() ) {
      // print a message when you have a significative difference since rval is computed
      if ( rVal > fScannedVariable->getMax()*(1.+1.E-12) )
         oocoutE((TObject*)0,InputArguments) << "HypoTestInverter::RunOnePoint - Out of range: using the upper bound "
                                             << fScannedVariable->getMax()
                                             << " on the scanned variable rather than " << rVal<< "\n";
     rVal = fScannedVariable->getMax();
   }

   // save old value
   double oldValue = fScannedVariable->getVal();

   // evaluate hybrid calculator at a single point
   fScannedVariable->setVal(rVal);
   // need to set value of rval in hybridcalculator
   // assume null model is S+B and alternate is B only
   const ModelConfig * sbModel = fCalculator0->GetNullModel();
   RooArgSet poi; poi.add(*sbModel->GetParametersOfInterest());
   // set poi to right values
   poi = RooArgSet(*fScannedVariable);
   const_cast<ModelConfig*>(sbModel)->SetSnapshot(poi);

   if (fVerbose > 0)
      oocoutP((TObject*)0,Eval) << "Running for " << fScannedVariable->GetName() << " = " << fScannedVariable->getVal() << endl;

   // compute the results
   HypoTestResult* result =   Eval(*fCalculator0,adaptive,clTarget);
   if (!result) {
      oocoutE((TObject*)0,Eval) << "HypoTestInverter - Error running point " << fScannedVariable->GetName() << " = " <<
   fScannedVariable->getVal() << endl;
      return false;
   }
   // in case of a dummy result
   const double nullPV = result->NullPValue();
   const double altPV = result->AlternatePValue();
   if (!std::isfinite(nullPV) || nullPV < 0. || nullPV > 1. || !std::isfinite(altPV) || altPV < 0. || altPV > 1.) {
      oocoutW((TObject*)0,Eval) << "HypoTestInverter - Skipping invalid result for  point " << fScannedVariable->GetName() << " = " <<
         fScannedVariable->getVal() << ". null p-value=" << nullPV << ", alternate p-value=" << altPV << endl;
      return false;
   }

   double lastXtested;
   if ( fResults->ArraySize()!=0 ) lastXtested = fResults->GetXValue(fResults->ArraySize()-1);
   else lastXtested = -999;

   if ( (std::abs(rVal) < 1 && TMath::AreEqualAbs(rVal, lastXtested,1.E-12) ) ||
        (std::abs(rVal) >= 1 && TMath::AreEqualRel(rVal, lastXtested,1.E-12) ) ) {

      oocoutI((TObject*)0,Eval) << "HypoTestInverter::RunOnePoint - Merge with previous result for "
                                << fScannedVariable->GetName() << " = " << rVal << std::endl;
      HypoTestResult* prevResult =  fResults->GetResult(fResults->ArraySize()-1);
      if (prevResult && prevResult->GetNullDistribution() && prevResult->GetAltDistribution()) {
         prevResult->Append(result);
         delete result;  // we can delete the result
      }
      else {
         // if it was empty we re-use it
         oocoutI((TObject*)0,Eval) << "HypoTestInverter::RunOnePoint - replace previous empty result\n";
         fResults->fYObjects.Remove( prevResult);
         fResults->fYObjects.Add(result);
      }

   } else {

     // fill the results in the HypoTestInverterResult array
     fResults->fXValues.push_back(rVal);
     fResults->fYObjects.Add(result);

   }

   fScannedVariable->setVal(oldValue);

   return true;
}

////////////////////////////////////////////////////////////////////////////////
/// Run an automatic scan until the desired accuracy is reached.
/// Start by default from the full interval (min,max) of the POI and then via bisection find the line crossing
/// the target line.
/// Optionally, a hint can be provided and the scan will be done closer to that value.
/// If by bisection the desired accuracy will not be reached, a fit to the points is performed.
/// \param[out] limit The limit.
/// \param[out] limitErr The error of the limit.
/// \param[in] absAccuracy Desired absolute accuracy.
/// \param[in] relAccuracy Desired relative accuracy.
/// \param[in] hint Hint to start from or nullptr for no hint.

bool HypoTestInverter::RunLimit(double &limit, double &limitErr, double absAccuracy, double relAccuracy, const double*hint) const {


   // routine from G. Petrucciani (from HiggsCombination CMS package)

   RooRealVar *r = fScannedVariable;

  if ((hint != 0) && (*hint > r->getMin())) {
     r->setMax(std::min<double>(3.0 * (*hint), r->getMax()));
     r->setMin(std::max<double>(0.3 * (*hint), r->getMin()));
     oocoutI((TObject*)0,InputArguments) << "HypoTestInverter::RunLimit - Use hint value " << *hint
                                         << " search in interval " << r->getMin() << " , " << r->getMax() << std::endl;
  }

  // if not specified use the default values for rel and absolute accuracy
  if (absAccuracy <= 0) absAccuracy = fgAbsAccuracy;
  if (relAccuracy <= 0) relAccuracy = fgRelAccuracy;

  typedef std::pair<double,double> CLs_t;
  double clsTarget = fSize;
  CLs_t clsMin(1,0), clsMax(0,0), clsMid(0,0);
  double rMin = r->getMin(), rMax = r->getMax();
  limit    = 0.5*(rMax + rMin);
  limitErr = 0.5*(rMax - rMin);
  bool done = false;

  TF1 expoFit("expoFit","[0]*exp([1]*(x-[2]))", rMin, rMax);

  fLimitPlot.reset(new TGraphErrors());

  if (fVerbose > 0) std::cout << "Search for upper limit to the limit" << std::endl;
  for (int tries = 0; tries < 6; ++tries) {
     if (! RunOnePoint(rMax) ) {
        oocoutE((TObject*)0,Eval) << "HypoTestInverter::RunLimit - Hypotest failed at upper limit of scan range: " << rMax << std::endl;
        rMax *= 0.95;
        continue;
     }
     clsMax = std::make_pair( fResults->GetLastYValue(), fResults->GetLastYError() );
     if (clsMax.first == 0 || clsMax.first + 3 * fabs(clsMax.second) < clsTarget ) break;
     rMax += rMax;
     if (tries == 5) {
        oocoutE((TObject*)0,Eval) << "HypoTestInverter::RunLimit - Cannot determine upper limit of scan range. At " << r->GetName()
                                  << " = " << rMax  << " still getting "
                                  << (fUseCLs ? "CLs" : "CLsplusb") << " = " << clsMax.first << std::endl;
        return false;
     }
  }
  if (fVerbose > 0) {
     oocoutI((TObject*)0,Eval) << "HypoTestInverter::RunLimit - Search for lower limit to the limit" << std::endl;
  }

  if ( fUseCLs && rMin == 0 ) {
     clsMin =  CLs_t(1,0);
  }
  else {
     if (! RunOnePoint(rMin) ) {
       oocoutE((TObject*)0,Eval) << "HypoTestInverter::RunLimit - Hypotest failed at lower limit of scan range: " << rMin << std::endl;
       return false;
     }
     clsMin = std::make_pair( fResults->GetLastYValue(), fResults->GetLastYError() );
  }
  if (clsMin.first != 1 && clsMin.first - 3 * fabs(clsMin.second) < clsTarget) {
     if (fUseCLs) {
        rMin = 0;
        clsMin = CLs_t(1,0); // this is always true for CLs
     } else {
        rMin = -rMax / 4;
        for (int tries = 0; tries < 6; ++tries) {
           if (! RunOnePoint(rMin) ) {
             oocoutE((TObject*)0,Eval) << "HypoTestInverter::RunLimit - Hypotest failed at lower limit of scan range: " << rMin << std::endl;
             rMin = rMin == 0. ? 0.1 : rMin * 1.1;
             continue;
           }
           clsMin = std::make_pair( fResults->GetLastYValue(), fResults->GetLastYError() );
           if (clsMin.first == 1 || clsMin.first - 3 * fabs(clsMin.second) > clsTarget) break;
           rMin += rMin;
           if (tries == 5) {
              oocoutE((TObject*)0,Eval) << "HypoTestInverter::RunLimit - Cannot determine lower limit of scan range. At " << r->GetName()
                                        << " = " << rMin << " still get " << (fUseCLs ? "CLs" : "CLsplusb")
                                        << " = " << clsMin.first << std::endl;
              return false;
           }
        }
     }
  }

  if (fVerbose > 0)
      oocoutI((TObject*)0,Eval) << "HypoTestInverter::RunLimit - Now doing proper bracketing & bisection" << std::endl;
  do {

     // break loop in case max toys is reached
     if (fMaxToys > 0 && fTotalToysRun > fMaxToys ) {
        oocoutW((TObject*)0,Eval) << "HypoTestInverter::RunLimit - maximum number of toys reached  " << std::endl;
        done = false; break;
     }


     // determine point by bisection or interpolation
     limit = 0.5*(rMin+rMax); limitErr = 0.5*(rMax-rMin);
     if (fgAlgo == "logSecant" && clsMax.first != 0) {
        double logMin = log(clsMin.first), logMax = log(clsMax.first), logTarget = log(clsTarget);
        limit = rMin + (rMax-rMin) * (logTarget - logMin)/(logMax - logMin);
        if (clsMax.second != 0 && clsMin.second != 0) {
           limitErr = hypot((logTarget-logMax) * (clsMin.second/clsMin.first), (logTarget-logMin) * (clsMax.second/clsMax.first));
           limitErr *= (rMax-rMin)/((logMax-logMin)*(logMax-logMin));
        }
     }
     r->setError(limitErr);

     // exit if reached accuracy on r
     if (limitErr < std::max(absAccuracy, relAccuracy * limit)) {
        if (fVerbose > 1)
            oocoutI((TObject*)0,Eval) << "HypoTestInverter::RunLimit - reached accuracy " << limitErr
                                      << " below " << std::max(absAccuracy, relAccuracy * limit)  << std::endl;
        done = true; break;
     }

     // evaluate point
     if (! RunOnePoint(limit, true, clsTarget) ) {
       oocoutE((TObject*)0,Eval) << "HypoTestInverter::RunLimit - Hypo test failed at x=" << limit << " when trying to find limit." << std::endl;
       return false;
     }
     clsMid = std::make_pair( fResults->GetLastYValue(), fResults->GetLastYError() );

     if (clsMid.second == -1) {
        std::cerr << "Hypotest failed" << std::endl;
        return false;
     }

     // if sufficiently far away, drop one of the points
     if (fabs(clsMid.first-clsTarget) >= 2*clsMid.second) {
       if ((clsMid.first>clsTarget) == (clsMax.first>clsTarget)) {
         rMax = limit; clsMax = clsMid;
       } else {
         rMin = limit; clsMin = clsMid;
       }
     } else {
       if (fVerbose > 0) std::cout << "Trying to move the interval edges closer" << std::endl;
       double rMinBound = rMin, rMaxBound = rMax;
       // try to reduce the size of the interval
       while (clsMin.second == 0 || fabs(rMin-limit) > std::max(absAccuracy, relAccuracy * limit)) {
         rMin = 0.5*(rMin+limit);
         if (!RunOnePoint(rMin,true, clsTarget) ) {
           oocoutE((TObject*)0,Eval) << "HypoTestInverter::RunLimit - Hypo test failed at x=" << rMin << " when trying to find limit from below." << std::endl;
           return false;
         }
         clsMin = std::make_pair( fResults->GetLastYValue(), fResults->GetLastYError() );
         if (fabs(clsMin.first-clsTarget) <= 2*clsMin.second) break;
         rMinBound = rMin;
       }
       while (clsMax.second == 0 || fabs(rMax-limit) > std::max(absAccuracy, relAccuracy * limit)) {
         rMax = 0.5*(rMax+limit);
         if (!RunOnePoint(rMax,true,clsTarget) ) {
           oocoutE((TObject*)0,Eval) << "HypoTestInverter::RunLimit - Hypo test failed at x=" << rMin << " when trying to find limit from above." << std::endl;
           return false;
         }
         clsMax = std::make_pair( fResults->GetLastYValue(), fResults->GetLastYError() );
         if (fabs(clsMax.first-clsTarget) <= 2*clsMax.second) break;
         rMaxBound = rMax;
       }
       expoFit.SetRange(rMinBound,rMaxBound);
       break;
     }
  } while (true);

  if (!done) { // didn't reach accuracy with scan, now do fit
      if (fVerbose) {
         oocoutI((TObject*)0,Eval) << "HypoTestInverter::RunLimit - Before fit   --- \n";
         std::cout << "Limit: " << r->GetName() << " < " << limit << " +/- " << limitErr << " [" << rMin << ", " << rMax << "]\n";
      }

      expoFit.FixParameter(0,clsTarget);
      expoFit.SetParameter(1,log(clsMax.first/clsMin.first)/(rMax-rMin));
      expoFit.SetParameter(2,limit);
      double rMinBound, rMaxBound; expoFit.GetRange(rMinBound, rMaxBound);
      limitErr = std::max(fabs(rMinBound-limit), fabs(rMaxBound-limit));
      int npoints = 0;

      HypoTestInverterPlot plot("plot","plot",fResults);
      fLimitPlot.reset(plot.MakePlot() );


      for (int j = 0; j < fLimitPlot->GetN(); ++j) {
         if (fLimitPlot->GetX()[j] >= rMinBound && fLimitPlot->GetX()[j] <= rMaxBound) npoints++;
      }
      for (int i = 0, imax = /*(readHybridResults_ ? 0 : */  8; i <= imax; ++i, ++npoints) {
          fLimitPlot->Sort();
          fLimitPlot->Fit(&expoFit,(fVerbose <= 1 ? "QNR EX0" : "NR EXO"));
          if (fVerbose) {
               oocoutI((TObject*)0,Eval) << "Fit to " << npoints << " points: " << expoFit.GetParameter(2) << " +/- " << expoFit.GetParError(2) << std::endl;
          }
          if ((rMin < expoFit.GetParameter(2))  && (expoFit.GetParameter(2) < rMax) && (expoFit.GetParError(2) < 0.5*(rMaxBound-rMinBound))) {
              // sanity check fit result
              limit = expoFit.GetParameter(2);
              limitErr = expoFit.GetParError(2);
              if (limitErr < std::max(absAccuracy, relAccuracy * limit)) break;
          }
          // add one point in the interval.
          double rTry = RooRandom::uniform()*(rMaxBound-rMinBound)+rMinBound;
          if (i != imax) {
             if (!RunOnePoint(rTry,true,clsTarget) ) return false;
             //eval(w, mc_s, mc_b, data, rTry, true, clsTarget);
          }

      }
  }

//if (!plot_.empty() && fLimitPlot.get()) {
  if (fLimitPlot.get() && fLimitPlot->GetN() > 0) {
       //new TCanvas("c1","c1");
      fLimitPlot->Sort();
      fLimitPlot->SetLineWidth(2);
      double xmin = r->getMin(), xmax = r->getMax();
      for (int j = 0; j < fLimitPlot->GetN(); ++j) {
        if (fLimitPlot->GetY()[j] > 1.4*clsTarget || fLimitPlot->GetY()[j] < 0.6*clsTarget) continue;
        xmin = std::min(fLimitPlot->GetX()[j], xmin);
        xmax = std::max(fLimitPlot->GetX()[j], xmax);
      }
      fLimitPlot->GetXaxis()->SetRangeUser(xmin,xmax);
      fLimitPlot->GetYaxis()->SetRangeUser(0.5*clsTarget, 1.5*clsTarget);
      fLimitPlot->Draw("AP");
      expoFit.Draw("SAME");
      TLine line(fLimitPlot->GetX()[0], clsTarget, fLimitPlot->GetX()[fLimitPlot->GetN()-1], clsTarget);
      line.SetLineColor(kRed); line.SetLineWidth(2); line.Draw();
      line.DrawLine(limit, 0, limit, fLimitPlot->GetY()[0]);
      line.SetLineWidth(1); line.SetLineStyle(2);
      line.DrawLine(limit-limitErr, 0, limit-limitErr, fLimitPlot->GetY()[0]);
      line.DrawLine(limit+limitErr, 0, limit+limitErr, fLimitPlot->GetY()[0]);
      //c1->Print(plot_.c_str());
  }

  oocoutI((TObject*)0,Eval) << "HypoTestInverter::RunLimit - Result:    \n"
                            << "\tLimit: " << r->GetName() << " < " << limit << " +/- " << limitErr << " @ " << (1-fSize) * 100 << "% CL\n";
  if (fVerbose > 1) oocoutI((TObject*)0,Eval) << "Total toys: " << fTotalToysRun << std::endl;

  // set value in results
  fResults->fUpperLimit = limit;
  fResults->fUpperLimitError = limitErr;
  fResults->fFittedUpperLimit = true;
  // lower limit are always min of p value
  fResults->fLowerLimit = fScannedVariable->getMin();
  fResults->fLowerLimitError = 0;
  fResults->fFittedLowerLimit = true;

  return true;
}

////////////////////////////////////////////////////////////////////////////////
/// get the distribution of lower limit
/// if rebuild = false (default) it will re-use the results of the scan done
/// for obtained the observed limit and no extra toys will be generated
/// if rebuild a new set of B toys will be done and the procedure will be repeated
/// for each toy

SamplingDistribution * HypoTestInverter::GetLowerLimitDistribution(bool rebuild, int nToys) {
   if (!rebuild) {
      if (!fResults) {
         oocoutE((TObject*)0,InputArguments) << "HypoTestInverter::GetLowerLimitDistribution(false) - result not existing\n";
         return 0;
      }
      return fResults->GetLowerLimitDistribution();
   }

   TList * clsDist = 0;
   TList * clsbDist = 0;
   if (fUseCLs) clsDist = &fResults->fExpPValues;
   else clsbDist = &fResults->fExpPValues;

   return RebuildDistributions(false, nToys,clsDist, clsbDist);

}

////////////////////////////////////////////////////////////////////////////////
/// get the distribution of lower limit
/// if rebuild = false (default) it will re-use the results of the scan done
/// for obtained the observed limit and no extra toys will be generated
/// if rebuild a new set of B toys will be done and the procedure will be repeated
/// for each toy
/// The nuisance parameter value used for rebuild is the current one in the model
/// so it is user responsibility to set to the desired value (nomi

SamplingDistribution * HypoTestInverter::GetUpperLimitDistribution(bool rebuild, int nToys) {
   if (!rebuild) {
      if (!fResults) {
         oocoutE((TObject*)0,InputArguments) << "HypoTestInverter::GetUpperLimitDistribution(false) - result not existing\n";
         return 0;
      }
      return fResults->GetUpperLimitDistribution();
   }

   TList * clsDist = 0;
   TList * clsbDist = 0;
   if (fUseCLs) clsDist = &fResults->fExpPValues;
   else clsbDist = &fResults->fExpPValues;

   return RebuildDistributions(true, nToys,clsDist, clsbDist);
}

////////////////////////////////////////////////////////////////////////////////

void HypoTestInverter::SetData(RooAbsData & data) {
   if (fCalculator0) fCalculator0->SetData(data);
}

////////////////////////////////////////////////////////////////////////////////
/// rebuild the sampling distributions by
/// generating some toys and find for each of them a new upper limit
/// Return the upper limit distribution and optionally also the pValue distributions for Cls, Clsb and Clbxs
/// as a TList for each scanned point
/// The method uses the present parameter value. It is user responsibility to give the current parameters to rebuild the distributions
/// It returns a upper or lower limit distribution depending on the isUpper flag, however it computes also the lower limit distribution and it is saved in the
/// output file as an histogram

SamplingDistribution * HypoTestInverter::RebuildDistributions(bool isUpper, int nToys, TList * clsDist, TList * clsbDist, TList * clbDist, const char *outputfile) {

   if (!fScannedVariable || !fCalculator0) return 0;
   // get first background snapshot
   const ModelConfig * bModel = fCalculator0->GetAlternateModel();
   const ModelConfig * sbModel = fCalculator0->GetNullModel();
   if (!bModel || ! sbModel) return 0;
   RooArgSet paramPoint;
   if (!sbModel->GetParametersOfInterest()) return 0;
   paramPoint.add(*sbModel->GetParametersOfInterest());

   const RooArgSet * poibkg = bModel->GetSnapshot();
   if (!poibkg) {
      oocoutW((TObject*)0,InputArguments) << "HypoTestInverter::RebuildDistribution - background snapshot not existing"
                                          << " assume is for POI = 0" << std::endl;
      fScannedVariable->setVal(0);
      paramPoint = RooArgSet(*fScannedVariable);
   }
   else
      paramPoint = *poibkg;
   // generate data at bkg parameter point

   ToyMCSampler * toymcSampler = dynamic_cast<ToyMCSampler *>(fCalculator0->GetTestStatSampler() );
   if (!toymcSampler) {
      oocoutE((TObject*)0,InputArguments) << "HypoTestInverter::RebuildDistribution - no toy MC sampler existing" << std::endl;
      return 0;
   }
   // set up test stat sampler in case of asymptotic calculator
   if (dynamic_cast<RooStats::AsymptoticCalculator*>(fCalculator0) ) {
      toymcSampler->SetObservables(*sbModel->GetObservables() );
      toymcSampler->SetParametersForTestStat(*sbModel->GetParametersOfInterest());
      toymcSampler->SetPdf(*sbModel->GetPdf());
      toymcSampler->SetNuisanceParameters(*sbModel->GetNuisanceParameters());
      if (sbModel->GetGlobalObservables() )  toymcSampler->SetGlobalObservables(*sbModel->GetGlobalObservables() );
      // set number of events
      if (!sbModel->GetPdf()->canBeExtended())
         toymcSampler->SetNEventsPerToy(1);
   }

   // loop on data to generate
   int nPoints = fNBins;

   bool storePValues = clsDist || clsbDist || clbDist;
   if (fNBins <=0  && storePValues) {
      oocoutW((TObject*)0,InputArguments) << "HypoTestInverter::RebuildDistribution - cannot return p values distribution with the auto scan" << std::endl;
      storePValues = false;
      nPoints = 0;
   }

   if (storePValues) {
      if (fResults) nPoints = fResults->ArraySize();
      if (nPoints <=0) {
         oocoutE((TObject*)0,InputArguments) << "HypoTestInverter - result is not existing and number of point to scan is not set"
                                             << std::endl;
         return 0;
      }
   }

   if (nToys <= 0) nToys = 100; // default value

   std::vector<std::vector<double> > CLs_values(nPoints);
   std::vector<std::vector<double> > CLsb_values(nPoints);
   std::vector<std::vector<double> > CLb_values(nPoints);

   if (storePValues) {
      for (int i = 0; i < nPoints; ++i) {
         CLs_values[i].reserve(nToys);
         CLb_values[i].reserve(nToys);
         CLsb_values[i].reserve(nToys);
      }
   }

   std::vector<double> limit_values; limit_values.reserve(nToys);

   oocoutI((TObject*)0,InputArguments) << "HypoTestInverter - rebuilding  the p value distributions by generating ntoys = "
                                       << nToys << std::endl;


   oocoutI((TObject*)0,InputArguments) << "Rebuilding using parameter of interest point:  ";
   RooStats::PrintListContent(paramPoint, oocoutI((TObject*)0,InputArguments) );
   if (sbModel->GetNuisanceParameters() ) {
      oocoutI((TObject*)0,InputArguments) << "And using nuisance parameters: ";
      RooStats::PrintListContent(*sbModel->GetNuisanceParameters(), oocoutI((TObject*)0,InputArguments) );
   }
   // save all parameters to restore them later
   assert(bModel->GetPdf() );
   assert(bModel->GetObservables() );
   RooArgSet * allParams = bModel->GetPdf()->getParameters( *bModel->GetObservables() );
   RooArgSet saveParams;
   allParams->snapshot(saveParams);

   TFile * fileOut = TFile::Open(outputfile,"RECREATE");
   if (!fileOut) {
      oocoutE((TObject*)0,InputArguments) << "HypoTestInverter - RebuildDistributions - Error opening file " << outputfile
                                          << " - the resulting limits will not be stored" << std::endl;
   }
   // create  temporary histograms to store the limit result
   TH1D * hL = new TH1D("lowerLimitDist","Rebuilt lower limit distribution",100,1.,0.);
   TH1D * hU = new TH1D("upperLimitDist","Rebuilt upper limit distribution",100,1.,0.);
   TH1D * hN = new TH1D("nObs","Observed events",100,1.,0.);
   hL->SetBuffer(2*nToys);
   hU->SetBuffer(2*nToys);
   std::vector<TH1*> hCLb;
   std::vector<TH1*> hCLsb;
   std::vector<TH1*> hCLs;
   if (storePValues) {
      for (int i = 0; i < nPoints; ++i) {
         hCLb.push_back(new TH1D(TString::Format("CLbDist_bin%d",i),"CLb distribution",100,1.,0.));
         hCLs.push_back(new TH1D(TString::Format("ClsDist_bin%d",i),"CLs distribution",100,1.,0.));
         hCLsb.push_back(new TH1D(TString::Format("CLsbDist_bin%d",i),"CLs+b distribution",100,1.,0.));
      }
   }


   // loop now on the toys
   for (int itoy = 0; itoy < nToys; ++itoy) {

      oocoutP((TObject*)0,Eval) << "\nHypoTestInverter - RebuildDistributions - running toy # " << itoy << " / "
                                       << nToys << std::endl;


      printf("\n\nshnapshot of s+b model \n");
      sbModel->GetSnapshot()->Print("v");

      // reset parameters to initial values to be sure in case they are not reset
      if (itoy> 0) *allParams = saveParams;

      // need to set th epdf to clear the cache in ToyMCSampler
      // pdf we must use is background pdf
      toymcSampler->SetPdf(*bModel->GetPdf() );


      RooAbsData * bkgdata = toymcSampler->GenerateToyData(paramPoint);

      double nObs = bkgdata->sumEntries();
      // for debugging in case of number counting models
      if (bkgdata->numEntries() ==1 && !bModel->GetPdf()->canBeExtended()) {
         oocoutP((TObject*)0,Generation) << "Generate observables are : ";
         RooArgList  genObs(*bkgdata->get(0));
         RooStats::PrintListContent(genObs, oocoutP((TObject*)0,Generation) );
         nObs = 0;
         for (int i = 0; i < genObs.getSize(); ++i) {
            RooRealVar * x = dynamic_cast<RooRealVar*>(&genObs[i]);
            if (x) nObs += x->getVal();
         }
      }
      hN->Fill(nObs);

      // by copying I will have the same min/max as previous ones
      HypoTestInverter inverter = *this;
      inverter.SetData(*bkgdata);

      // print global observables
      auto gobs = bModel->GetPdf()->getVariables()->selectCommon(* sbModel->GetGlobalObservables() );
      gobs->Print("v");

      HypoTestInverterResult * r  = inverter.GetInterval();

      if (r == 0) continue;

      double value = (isUpper) ? r->UpperLimit() : r->LowerLimit();
      limit_values.push_back( value );
      hU->Fill(r->UpperLimit() );
      hL->Fill(r->LowerLimit() );


      std::cout << "The computed upper limit for toy #" << itoy << " is " << value << std::endl;

      // write every 10 toys
      if (itoy%10 == 0 || itoy == nToys-1) {
         hU->Write("",TObject::kOverwrite);
         hL->Write("",TObject::kOverwrite);
         hN->Write("",TObject::kOverwrite);
      }

      if (!storePValues) continue;

      if (nPoints < r->ArraySize()) {
         oocoutW((TObject*)0,InputArguments) << "HypoTestInverter: skip extra points" << std::endl;
      }
      else if (nPoints > r->ArraySize()) {
         oocoutW((TObject*)0,InputArguments) << "HypoTestInverter: missing some points" << std::endl;
      }


      for (int ipoint = 0; ipoint < nPoints; ++ipoint) {
         HypoTestResult * hr = r->GetResult(ipoint);
         if (hr) {
            CLs_values[ipoint].push_back( hr->CLs() );
            CLsb_values[ipoint].push_back( hr->CLsplusb() );
            CLb_values[ipoint].push_back( hr->CLb() );
            hCLs[ipoint]->Fill(  hr->CLs() );
            hCLb[ipoint]->Fill(  hr->CLb() );
            hCLsb[ipoint]->Fill(  hr->CLsplusb() );
         }
         else {
            oocoutW((TObject*)0,InputArguments) << "HypoTestInverter: missing result for point: x = "
                                                << fResults->GetXValue(ipoint) <<  std::endl;
         }
      }
      // write every 10 toys
      if (itoy%10 == 0 || itoy == nToys-1) {
         for (int ipoint = 0; ipoint < nPoints; ++ipoint) {
            hCLs[ipoint]->Write("",TObject::kOverwrite);
            hCLb[ipoint]->Write("",TObject::kOverwrite);
            hCLsb[ipoint]->Write("",TObject::kOverwrite);
         }
      }


      delete r;
      delete bkgdata;
   }


   if (storePValues) {
      if (clsDist) clsDist->SetOwner(true);
      if (clbDist) clbDist->SetOwner(true);
      if (clsbDist) clsbDist->SetOwner(true);

      oocoutI((TObject*)0,InputArguments) << "HypoTestInverter: storing rebuilt p values  " << std::endl;

      for (int ipoint = 0; ipoint < nPoints; ++ipoint) {
         if (clsDist) {
            TString name = TString::Format("CLs_distrib_%d",ipoint);
            clsDist->Add( new SamplingDistribution(name,name,CLs_values[ipoint] ) );
         }
         if (clbDist) {
            TString name = TString::Format("CLb_distrib_%d",ipoint);
            clbDist->Add( new SamplingDistribution(name,name,CLb_values[ipoint] ) );
         }
         if (clsbDist) {
            TString name = TString::Format("CLsb_distrib_%d",ipoint);
            clsbDist->Add( new SamplingDistribution(name,name,CLsb_values[ipoint] ) );
         }
      }
   }

   if (fileOut) {
      fileOut->Close();
   }
   else {
      // delete all the histograms
      delete hL;
      delete hU;
      for (int i = 0; i < nPoints && storePValues; ++i) {
         delete hCLs[i];
         delete hCLb[i];
         delete hCLsb[i];
      }
   }

   const char * disName = (isUpper) ? "upperLimit_dist" : "lowerLimit_dist";
   return new SamplingDistribution(disName, disName, limit_values);
}
