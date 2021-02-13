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

/*! \class TMVA::OptimizeConfigParameters
\ingroup TMVA

*/

#include "TMVA/OptimizeConfigParameters.h"
#include "TMVA/Config.h"
#include "TMVA/DataSet.h"
#include "TMVA/DataSetInfo.h"
#include "TMVA/Event.h"
#include "TMVA/IFitterTarget.h"
#include "TMVA/FitterBase.h"
#include "TMVA/GeneticFitter.h"
#include "TMVA/IMethod.h"
#include "TMVA/Interval.h"
#include "TMVA/MethodBase.h"
#include "TMVA/MethodFDA.h"
#include "TMVA/MsgLogger.h"
#include "TMVA/MinuitFitter.h"
#include "TMVA/PDF.h"
#include "TMVA/Tools.h"
#include "TMVA/Types.h"

#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"

#include <cstdlib>
#include <limits>


ClassImp(TMVA::OptimizeConfigParameters);

////////////////////////////////////////////////////////////////////////////////
/// Constructor which sets either "Classification or Regression"

TMVA::OptimizeConfigParameters::OptimizeConfigParameters(MethodBase * const method, std::map<TString,TMVA::Interval*> tuneParameters, TString fomType, TString optimizationFitType)
:  fMethod(method),
   fTuneParameters(tuneParameters),
   fFOMType(fomType),
   fOptimizationFitType(optimizationFitType),
   fMvaSig(NULL),
   fMvaBkg(NULL),
   fMvaSigFineBin(NULL),
   fMvaBkgFineBin(NULL),
   fNotDoneYet(kFALSE)
{
   std::string name = "OptimizeConfigParameters_";
   name += std::string(GetMethod()->GetName());
   fLogger = new MsgLogger(name);
   if (fMethod->DoRegression()){
      Log() << kFATAL << " ERROR: Sorry, Regression is not yet implement for automatic parameter optimization"
            << " --> exit" << Endl;
   }

   Log() << kINFO << "Automatic optimisation of tuning parameters in "
         << GetMethod()->GetName() << " uses:" << Endl;

   std::map<TString,TMVA::Interval*>::iterator it;
   for (it=fTuneParameters.begin(); it!=fTuneParameters.end();++it) {
      Log() << kINFO << it->first
            << " in range from: " << it->second->GetMin()
            << " to: " << it->second->GetMax()
            << " in : " << it->second->GetNbins()  << " steps"
            << Endl;
   }
   Log() << kINFO << " using the options: " << fFOMType << " and " << fOptimizationFitType << Endl;
}

////////////////////////////////////////////////////////////////////////////////
/// the destructor (delete the OptimizeConfigParameters, store the graph and .. delete it)

TMVA::OptimizeConfigParameters::~OptimizeConfigParameters()
{
   if(!GetMethod()->IsSilentFile()) GetMethod()->BaseDir()->cd();
   Int_t n=Int_t(fFOMvsIter.size());
   Float_t *x = new Float_t[n];
   Float_t *y = new Float_t[n];
   Float_t  ymin=(Float_t)+999999999;
   Float_t  ymax=(Float_t)-999999999;

   for (Int_t i=0;i<n;i++){
      x[i] = Float_t(i);
      y[i] = fFOMvsIter[i];
      if (ymin>y[i]) ymin=y[i];
      if (ymax<y[i]) ymax=y[i];
   }

   TH2D   *h=new TH2D(TString(GetMethod()->GetName())+"_FOMvsIterFrame","",2,0,n,2,ymin*0.95,ymax*1.05);
   h->SetXTitle("#iteration "+fOptimizationFitType);
   h->SetYTitle(fFOMType);
   TGraph *gFOMvsIter = new TGraph(n,x,y);
   gFOMvsIter->SetName((TString(GetMethod()->GetName())+"_FOMvsIter").Data());
   if(!GetMethod()->IsSilentFile()) gFOMvsIter->Write();
   if(!GetMethod()->IsSilentFile()) h->Write();

   delete [] x;
   delete [] y;
   // delete fFOMvsIter;
}

////////////////////////////////////////////////////////////////////////////////

std::map<TString,Double_t> TMVA::OptimizeConfigParameters::optimize()
{
   if      (fOptimizationFitType == "Scan"    ) this->optimizeScan();
   else if (fOptimizationFitType == "FitGA" || fOptimizationFitType == "Minuit" ) this->optimizeFit();
   else {
      Log() << kFATAL << "You have chosen as optimization type " << fOptimizationFitType
            << " that is not (yet) coded --> exit()" << Endl;
   }

   Log() << kINFO << "For " << GetMethod()->GetName() << " the optimized Parameters are: " << Endl;
   std::map<TString,Double_t>::iterator it;
   for(it=fTunedParameters.begin(); it!= fTunedParameters.end(); ++it){
      Log() << kINFO << it->first << " = " << it->second << Endl;
   }
   return fTunedParameters;

}

////////////////////////////////////////////////////////////////////////////////
/// helper function to scan through the all the combinations in the
/// parameter space

std::vector< int > TMVA::OptimizeConfigParameters::GetScanIndices( int val, std::vector<int> base){
   std::vector < int > indices;
   for (UInt_t i=0; i< base.size(); i++){
      indices.push_back(val % base[i] );
      val = int( floor( float(val) / float(base[i]) ) );
   }
   return indices;
}

////////////////////////////////////////////////////////////////////////////////
/// do the actual optimization using a simple scan method,
/// i.e. calculate the FOM for
/// different tuning paraemters and remember which one is
/// gave the best FOM

void TMVA::OptimizeConfigParameters::optimizeScan()
{

   Double_t      bestFOM=-1000000, currentFOM;

   std::map<TString,Double_t> currentParameters;
   std::map<TString,TMVA::Interval*>::iterator it;

   // for the scan, start at the lower end of the interval and then "move upwards"
   // initialize all parameters in currentParameter
   currentParameters.clear();
   fTunedParameters.clear();

   for (it=fTuneParameters.begin(); it!=fTuneParameters.end(); ++it){
      currentParameters.insert(std::pair<TString,Double_t>(it->first,it->second->GetMin()));
      fTunedParameters.insert(std::pair<TString,Double_t>(it->first,it->second->GetMin()));
   }
   // now loop over all the parameters and get for each combination the figure of merit

   // in order to loop over all the parameters, I first create an "array" (tune parameters)
   // of arrays (the different values of the tune parameter)

   std::vector< std::vector <Double_t> > v;
   for (it=fTuneParameters.begin(); it!=fTuneParameters.end(); ++it){
      std::vector< Double_t > tmp;
      for (Int_t k=0; k<it->second->GetNbins(); k++){
         tmp.push_back(it->second->GetElement(k));
      }
      v.push_back(tmp);
   }
   Int_t Ntot = 1;
   std::vector< int > Nindividual;
   for (UInt_t i=0; i<v.size(); i++) {
      Ntot *= v[i].size();
      Nindividual.push_back(v[i].size());
   }
   //loop on the total number of different combinations

   for (int i=0; i<Ntot; i++){
      UInt_t index=0;
      std::vector<int> indices = GetScanIndices(i, Nindividual );
      for (it=fTuneParameters.begin(), index=0; index< indices.size(); ++index, ++it){
         currentParameters[it->first] = v[index][indices[index]];
      }
      Log() << kINFO << "--------------------------" << Endl;
      Log() << kINFO <<"Settings being evaluated:" << Endl;
      for (std::map<TString,Double_t>::iterator it_print=currentParameters.begin();
           it_print!=currentParameters.end(); ++it_print){
         Log() << kINFO << "  " << it_print->first  << " = " << it_print->second << Endl;
      }

      GetMethod()->Reset();
      GetMethod()->SetTuneParameters(currentParameters);
      // now do the training for the current parameters:
      if(!GetMethod()->IsSilentFile()) GetMethod()->BaseDir()->cd();
      if (i==0) GetMethod()->GetTransformationHandler().CalcTransformations(
                                                                            GetMethod()->Data()->GetEventCollection());
      Event::SetIsTraining(kTRUE);
      GetMethod()->Train();
      Event::SetIsTraining(kFALSE);
      currentFOM = GetFOM();
      Log() << kINFO << "FOM was found : " << currentFOM << "; current best is " << bestFOM << Endl;

      if (currentFOM > bestFOM) {
         bestFOM = currentFOM;
         for (std::map<TString,Double_t>::iterator iter=currentParameters.begin();
              iter != currentParameters.end(); ++iter){
            fTunedParameters[iter->first]=iter->second;
         }
      }
   }

   GetMethod()->Reset();
   GetMethod()->SetTuneParameters(fTunedParameters);
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::OptimizeConfigParameters::optimizeFit()
{
   // ranges (intervals) in which the fit varies the parameters
   std::vector<TMVA::Interval*> ranges; // intervals of the fit ranges
   std::map<TString, TMVA::Interval*>::iterator it;
   std::vector<Double_t> pars;    // current (starting) fit parameters

   for (it=fTuneParameters.begin(); it != fTuneParameters.end(); ++it){
      ranges.push_back(new TMVA::Interval(*(it->second)));
      pars.push_back( (it->second)->GetMean() );  // like this the order is "right". Always keep the
      // order in the vector "pars" the same as the iterator
      // iterates through the tuneParameters !!!!
   }

   // added to allow for transformation on input variables i.e. norm
   GetMethod()->GetTransformationHandler().CalcTransformations(GetMethod()->Data()->GetEventCollection());

   // create the fitter

   FitterBase* fitter = NULL;

   if ( fOptimizationFitType == "Minuit"  ) {
      TString opt="FitStrategy=0:UseImprove=False:UseMinos=False:Tolerance=100";
      if (!TMVA::gConfig().IsSilent() ) opt += TString(":PrintLevel=0");

      fitter = new MinuitFitter(  *this,
                                  "FitterMinuit_BDTOptimize",
                                  ranges, opt );
   }else if ( fOptimizationFitType == "FitGA"  ) {
      TString opt="PopSize=20:Steps=30:Cycles=3:ConvCrit=0.01:SaveBestCycle=5";
      fitter = new GeneticFitter( *this,
                                  "FitterGA_BDTOptimize",
                                  ranges, opt );
   } else {
      Log() << kWARNING << " you did not specify a valid OptimizationFitType "
            << " will use the default (FitGA) " << Endl;
      TString opt="PopSize=20:Steps=30:Cycles=3:ConvCrit=0.01:SaveBestCycle=5";
      fitter = new GeneticFitter( *this,
                                  "FitterGA_BDTOptimize",
                                  ranges, opt );
   }

   fitter->CheckForUnusedOptions();

   // perform the fit
   fitter->Run(pars);

   // clean up
   for (UInt_t ipar=0; ipar<ranges.size(); ipar++) delete ranges[ipar];

   GetMethod()->Reset();

   fTunedParameters.clear();
   Int_t jcount=0;
   for (it=fTuneParameters.begin(); it!=fTuneParameters.end(); ++it){
      fTunedParameters.insert(std::pair<TString,Double_t>(it->first,pars[jcount++]));
   }

   GetMethod()->SetTuneParameters(fTunedParameters);

}

////////////////////////////////////////////////////////////////////////////////
/// return the estimator (from current FOM) for the fitting interface

Double_t TMVA::OptimizeConfigParameters::EstimatorFunction( std::vector<Double_t> & pars)
{
   std::map< std::vector<Double_t> , Double_t>::const_iterator iter;
   iter = fAlreadyTrainedParCombination.find(pars);

   if (iter != fAlreadyTrainedParCombination.end()) {
      // std::cout << "I  had trained  Depth=" <<Int_t(pars[0])
      //           <<" MinEv=" <<Int_t(pars[1])
      //           <<" already --> FOM="<< iter->second <<std::endl;
      return iter->second;
   }else{
      std::map<TString,Double_t> currentParameters;
      Int_t icount =0; // map "pars" to the  map of Tuneparameter, make sure
                       // you never screw up this order!!
      std::map<TString, TMVA::Interval*>::iterator it;
      for (it=fTuneParameters.begin(); it!=fTuneParameters.end(); ++it){
         currentParameters[it->first] = pars[icount++];
      }
      GetMethod()->Reset();
      GetMethod()->SetTuneParameters(currentParameters);
      if(!GetMethod()->IsSilentFile()) GetMethod()->BaseDir()->cd();

      if (fNotDoneYet){
         GetMethod()->GetTransformationHandler().
            CalcTransformations(GetMethod()->Data()->GetEventCollection());
         fNotDoneYet=kFALSE;
      }
      Event::SetIsTraining(kTRUE);
      GetMethod()->Train();
      Event::SetIsTraining(kFALSE);


      Double_t currentFOM = GetFOM();

      fAlreadyTrainedParCombination.insert(std::make_pair(pars,-currentFOM));
      return  -currentFOM;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return the Figure of Merit (FOM) used in the parameter
///  optimization process

Double_t TMVA::OptimizeConfigParameters::GetFOM()
{
   auto parsePercent = [this](TString input) -> Double_t {
      // Expects input e.g. SigEffAtBkgEff0 (14 chars) followed by a fraction
      // either as e.g. 01 or .01 (meaning the same thing 1 %).
      TString percent = TString(input(14, input.Sizeof()));
      if (!percent.CountChar('.')) percent.Insert(1,".");

      if (percent.IsFloat()) {
         return percent.Atof();
      } else {
         Log() << kFATAL << " ERROR, " << percent << " in " << fFOMType
               << " is not a valid floating point number" << Endl;
         return 0; // Cannot happen
      }
   };

   Double_t fom = 0;
   if (fMethod->DoRegression()){
      std::cout << " ERROR: Sorry, Regression is not yet implement for automatic parameter optimisation"
                << " --> exit" << std::endl;
      std::exit(1);
   } else {
      if      (fFOMType == "Separation")  fom = GetSeparation();
      else if (fFOMType == "ROCIntegral") fom = GetROCIntegral();
      else if (fFOMType.BeginsWith("SigEffAtBkgEff0")) fom = GetSigEffAtBkgEff(parsePercent(fFOMType));
      else if (fFOMType.BeginsWith("BkgRejAtSigEff0")) fom = GetBkgRejAtSigEff(parsePercent(fFOMType));
      else if (fFOMType.BeginsWith("BkgEffAtSigEff0")) fom = GetBkgEffAtSigEff(parsePercent(fFOMType));
      else {
         Log()<< kFATAL << " ERROR, you've specified as Figure of Merit in the "
              << " parameter optimisation " << fFOMType << " which has not"
              << " been implemented yet!! ---> exit " << Endl;
      }
   }

   fFOMvsIter.push_back(fom);
   //   std::cout << "fom="<<fom<<std::endl; // should write that into a debug log (as option)
   return fom;
}

////////////////////////////////////////////////////////////////////////////////
/// fill the private histograms with the mva distributions for sig/bkg

void TMVA::OptimizeConfigParameters::GetMVADists()
{
   if (fMvaSig) fMvaSig->Delete();
   if (fMvaBkg) fMvaBkg->Delete();
   if (fMvaSigFineBin) fMvaSigFineBin->Delete();
   if (fMvaBkgFineBin) fMvaBkgFineBin->Delete();

   // maybe later on this should be done a bit more clever (time consuming) by
   // first determining proper ranges, removing outliers, as we do in the
   // MVA output calculation in MethodBase::TestClassifier...
   // --> then it might be possible also to use the splined PDF's which currently
   // doesn't seem to work

   fMvaSig        = new TH1D("fMvaSig","",100,-1.5,1.5); //used for spline fit
   fMvaBkg        = new TH1D("fMvaBkg","",100,-1.5,1.5); //used for spline fit
   fMvaSigFineBin = new TH1D("fMvaSigFineBin","",100000,-1.5,1.5);
   fMvaBkgFineBin = new TH1D("fMvaBkgFineBin","",100000,-1.5,1.5);

   const std::vector< Event*> events=fMethod->Data()->GetEventCollection(Types::kTesting);

   UInt_t signalClassNr = fMethod->DataInfo().GetClassInfo("Signal")->GetNumber();

   //   fMethod->GetTransformationHandler().CalcTransformations(fMethod->Data()->GetEventCollection(Types::kTesting));

   for (UInt_t iev=0; iev < events.size() ; iev++){
      //      std::cout << " GetMVADists event " << iev << std::endl;
      //      std::cout << " Class  = " << events[iev]->GetClass() << std::endl;
      //         std::cout << " MVA Value = " << fMethod->GetMvaValue(events[iev]) << std::endl;
      if (events[iev]->GetClass() == signalClassNr) {
         fMvaSig->Fill(fMethod->GetMvaValue(events[iev]),events[iev]->GetWeight());
         fMvaSigFineBin->Fill(fMethod->GetMvaValue(events[iev]),events[iev]->GetWeight());
      } else {
         fMvaBkg->Fill(fMethod->GetMvaValue(events[iev]),events[iev]->GetWeight());
         fMvaBkgFineBin->Fill(fMethod->GetMvaValue(events[iev]),events[iev]->GetWeight());
      }
   }
}
////////////////////////////////////////////////////////////////////////////////
/// return the separation between the signal and background
/// MVA ouput distribution

Double_t TMVA::OptimizeConfigParameters::GetSeparation()
{
   GetMVADists();
   if (1){
      PDF *splS = new PDF( " PDF Sig", fMvaSig, PDF::kSpline2 );
      PDF *splB = new PDF( " PDF Bkg", fMvaBkg, PDF::kSpline2 );
      return gTools().GetSeparation(*splS,*splB);
   }else{
      std::cout << "Separation calculation via histograms (not PDFs) seems to give still strange results!! Don't do that, check!!"<<std::endl;
      return gTools().GetSeparation(fMvaSigFineBin,fMvaBkgFineBin); // somehow still gives strange results!!!! Check!!!
   }
}

////////////////////////////////////////////////////////////////////////////////
/// calculate the area (integral) under the ROC curve as a
/// overall quality measure of the classification
///
/// making pdfs out of the MVA-output distributions doesn't work
/// reliably for cases where the MVA-output isn't a smooth distribution.
/// this happens "frequently" in BDTs for example when the number of
/// trees is small resulting in only some discrete possible MVA output values.
/// (I still leave the code here, but use this with care!!! The default
/// however is to use the distributions!!!

Double_t TMVA::OptimizeConfigParameters::GetROCIntegral()
{
   GetMVADists();

   Double_t integral = 0;
   if (0){
      PDF *pdfS = new PDF( " PDF Sig", fMvaSig, PDF::kSpline2 );
      PDF *pdfB = new PDF( " PDF Bkg", fMvaBkg, PDF::kSpline2 );

      Double_t xmin = TMath::Min(pdfS->GetXmin(), pdfB->GetXmin());
      Double_t xmax = TMath::Max(pdfS->GetXmax(), pdfB->GetXmax());

      UInt_t   nsteps = 1000;
      Double_t step = (xmax-xmin)/Double_t(nsteps);
      Double_t cut = xmin;
      for (UInt_t i=0; i<nsteps; i++){
         integral += (1-pdfB->GetIntegral(cut,xmax)) * pdfS->GetVal(cut);
         cut+=step;
      }
      integral*=step;
   }else{
      // sanity checks
      if ( (fMvaSigFineBin->GetXaxis()->GetXmin() !=  fMvaBkgFineBin->GetXaxis()->GetXmin()) ||
           (fMvaSigFineBin->GetNbinsX() !=  fMvaBkgFineBin->GetNbinsX()) ){
         std::cout << " Error in OptimizeConfigParameters GetROCIntegral, unequal histograms for sig and bkg.." << std::endl;
         std::exit(1);
      }else{

         Double_t *cumulator  = fMvaBkgFineBin->GetIntegral();
         Int_t    nbins       = fMvaSigFineBin->GetNbinsX();
         // get the true signal integral (ComputeIntegral just return 1 as they
         // automatically normalize. IN ADDITION, they do not account for variable
         // bin sizes (which you might perhaps use later on for the fMvaSig/Bkg histograms)
         Double_t sigIntegral = 0;
         for (Int_t ibin=1; ibin<=nbins; ibin++){
            sigIntegral += fMvaSigFineBin->GetBinContent(ibin) * fMvaSigFineBin->GetBinWidth(ibin);
         }
         //gTools().NormHist( fMvaSigFineBin  ); // also doesn't  use variable bin width. And calls TH1::Scale, which oddly enough does not change the SumOfWeights !!!

         for (Int_t ibin=1; ibin <= nbins; ibin++){ // don't include under- and overflow bin
            integral += (cumulator[ibin]) * fMvaSigFineBin->GetBinContent(ibin)/sigIntegral * fMvaSigFineBin->GetBinWidth(ibin) ;
         }
      }
   }

   return integral;
}

////////////////////////////////////////////////////////////////////////////////
/// calculate the signal efficiency for a given background efficiency

Double_t TMVA::OptimizeConfigParameters::GetSigEffAtBkgEff(Double_t bkgEff)
{
   GetMVADists();
   Double_t sigEff=0;

   // sanity checks
   if ( (fMvaSigFineBin->GetXaxis()->GetXmin() !=  fMvaBkgFineBin->GetXaxis()->GetXmin()) ||
        (fMvaSigFineBin->GetNbinsX() !=  fMvaBkgFineBin->GetNbinsX()) ){
      std::cout << " Error in OptimizeConfigParameters GetSigEffAt, unequal histograms for sig and bkg.." << std::endl;
      std::exit(1);
   }else{
      Double_t *bkgCumulator   = fMvaBkgFineBin->GetIntegral();
      Double_t *sigCumulator   = fMvaSigFineBin->GetIntegral();

      Int_t nbins=fMvaBkgFineBin->GetNbinsX();
      Int_t ibin=0;

      // std::cout << " bkgIntegral="<<bkgIntegral
      //           << " sigIntegral="<<sigIntegral
      //           << " bkgCumulator[nbins]="<<bkgCumulator[nbins]
      //           << " sigCumulator[nbins]="<<sigCumulator[nbins]
      //           << std::endl;

      while (bkgCumulator[nbins-ibin] > (1-bkgEff)) {
         sigEff = sigCumulator[nbins]-sigCumulator[nbins-ibin];
         ibin++;
      }
   }
   return sigEff;
}


////////////////////////////////////////////////////////////////////////////////
/// calculate the background efficiency for a given signal efficiency
///
/// adapted by marc-olivier.bettler@cern.ch

Double_t TMVA::OptimizeConfigParameters::GetBkgEffAtSigEff(Double_t sigEff)
{
   GetMVADists();
   Double_t bkgEff=0;

   // sanity checks
   if ( (fMvaSigFineBin->GetXaxis()->GetXmin() !=  fMvaBkgFineBin->GetXaxis()->GetXmin()) ||
        (fMvaSigFineBin->GetNbinsX() !=  fMvaBkgFineBin->GetNbinsX()) ){
      std::cout << " Error in OptimizeConfigParameters GetBkgEffAt, unequal histograms for sig and bkg.." << std::endl;
      std::exit(1);
   }else{

      Double_t *bkgCumulator   = fMvaBkgFineBin->GetIntegral();
      Double_t *sigCumulator   = fMvaSigFineBin->GetIntegral();

      Int_t nbins=fMvaBkgFineBin->GetNbinsX();
      Int_t ibin=0;

      // std::cout << " bkgIntegral="<<bkgIntegral
      //           << " sigIntegral="<<sigIntegral
      //           << " bkgCumulator[nbins]="<<bkgCumulator[nbins]
      //           << " sigCumulator[nbins]="<<sigCumulator[nbins]
      //           << std::endl;

      while ( sigCumulator[nbins]-sigCumulator[nbins-ibin] < sigEff) {
         bkgEff = bkgCumulator[nbins]-bkgCumulator[nbins-ibin];
         ibin++;
      }
   }
   return bkgEff;
}

////////////////////////////////////////////////////////////////////////////////
/// calculate the background rejection for a given signal efficiency
///
/// adapted by marc-olivier.bettler@cern.ch

Double_t TMVA::OptimizeConfigParameters::GetBkgRejAtSigEff(Double_t sigEff)
{
   GetMVADists();
   Double_t bkgRej=0;

   // sanity checks
   if ( (fMvaSigFineBin->GetXaxis()->GetXmin() !=  fMvaBkgFineBin->GetXaxis()->GetXmin()) ||
        (fMvaSigFineBin->GetNbinsX() !=  fMvaBkgFineBin->GetNbinsX()) ){
      std::cout << " Error in OptimizeConfigParameters GetBkgEffAt, unequal histograms for sig and bkg.." << std::endl;
      std::exit(1);
   }else{

      Double_t *bkgCumulator   = fMvaBkgFineBin->GetIntegral();
      Double_t *sigCumulator   = fMvaSigFineBin->GetIntegral();

      Int_t nbins=fMvaBkgFineBin->GetNbinsX();
      Int_t ibin=0;

      // std::cout << " bkgIntegral="<<bkgIntegral
      //           << " sigIntegral="<<sigIntegral
      //           << " bkgCumulator[nbins]="<<bkgCumulator[nbins]
      //           << " sigCumulator[nbins]="<<sigCumulator[nbins]
      //           << std::endl;

      while ( sigCumulator[nbins]-sigCumulator[nbins-ibin] < sigEff) {
         bkgRej = bkgCumulator[nbins-ibin];
         ibin++;
      }
   }
   return bkgRej;
}
