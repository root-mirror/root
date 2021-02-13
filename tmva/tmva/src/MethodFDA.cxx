// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Peter Speckmayer, Joerg Stelzer

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : MethodFDA                                                             *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Implementation                                                            *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker  <Andreas.Hocker@cern.ch> - CERN, Switzerland             *
 *      Peter Speckmayer <speckmay@mail.cern.ch>  - CERN, Switzerland             *
 *      Joerg Stelzer    <stelzer@cern.ch>        - DESY, Germany                 *
 *      Maciej Kruk      <mkruk@cern.ch>          - IFJ PAN & AGH, Poland         *
 *                                                                                *
 * Copyright (c) 2005-2006:                                                       *
 *      CERN, Switzerland                                                         *
 *      MPI-K Heidelberg, Germany                                                 *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

/*! \class TMVA::MethodFDA
\ingroup TMVA

Function discriminant analysis (FDA).

This simple classifier
fits any user-defined TFormula (via option configuration string) to
the training data by requiring a formula response of 1 (0) to signal
(background) events. The parameter fitting is done via the abstract
class FitterBase, featuring Monte Carlo sampling, Genetic
Algorithm, Simulated Annealing, MINUIT and combinations of these.

Can compute regression value for one dimensional output
*/

#include "TMVA/MethodFDA.h"

#include "TMVA/ClassifierFactory.h"
#include "TMVA/Config.h"
#include "TMVA/Configurable.h"
#include "TMVA/DataSetInfo.h"
#include "TMVA/FitterBase.h"
#include "TMVA/GeneticFitter.h"
#include "TMVA/Interval.h"
#include "TMVA/IFitterTarget.h"
#include "TMVA/IMethod.h"
#include "TMVA/MCFitter.h"
#include "TMVA/MethodBase.h"
#include "TMVA/MinuitFitter.h"
#include "TMVA/MsgLogger.h"
#include "TMVA/Timer.h"
#include "TMVA/Tools.h"
#include "TMVA/TransformationHandler.h"
#include "TMVA/Types.h"
#include "TMVA/SimulatedAnnealingFitter.h"

#include "TList.h"
#include "TFormula.h"
#include "TString.h"
#include "TObjString.h"
#include "TRandom3.h"
#include "TMath.h"

#include <iostream>
#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <sstream>

using std::stringstream;

REGISTER_METHOD(FDA)

ClassImp(TMVA::MethodFDA);

////////////////////////////////////////////////////////////////////////////////
/// standard constructor

   TMVA::MethodFDA::MethodFDA( const TString& jobName,
                               const TString& methodTitle,
                               DataSetInfo& theData,
                               const TString& theOption)
   : MethodBase( jobName, Types::kFDA, methodTitle, theData, theOption),
   IFitterTarget   (),
   fFormula        ( 0 ),
   fNPars          ( 0 ),
   fFitter         ( 0 ),
   fConvergerFitter( 0 ),
   fSumOfWeightsSig( 0 ),
   fSumOfWeightsBkg( 0 ),
   fSumOfWeights   ( 0 ),
   fOutputDimensions( 0 )
{
}

////////////////////////////////////////////////////////////////////////////////
/// constructor from weight file

TMVA::MethodFDA::MethodFDA( DataSetInfo& theData,
                            const TString& theWeightFile)
   : MethodBase( Types::kFDA, theData, theWeightFile),
     IFitterTarget   (),
     fFormula        ( 0 ),
     fNPars          ( 0 ),
     fFitter         ( 0 ),
     fConvergerFitter( 0 ),
     fSumOfWeightsSig( 0 ),
     fSumOfWeightsBkg( 0 ),
     fSumOfWeights   ( 0 ),
     fOutputDimensions( 0 )
{
}

////////////////////////////////////////////////////////////////////////////////
/// default initialisation

void TMVA::MethodFDA::Init( void )
{
   fNPars    = 0;

   fBestPars.clear();

   fSumOfWeights    = 0;
   fSumOfWeightsSig = 0;
   fSumOfWeightsBkg = 0;

   fFormulaStringP  = "";
   fParRangeStringP = "";
   fFormulaStringT  = "";
   fParRangeStringT = "";

   fFitMethod       = "";
   fConverger       = "";

   if( DoMulticlass() )
      if (fMulticlassReturnVal == NULL) fMulticlassReturnVal = new std::vector<Float_t>();

}

////////////////////////////////////////////////////////////////////////////////
/// define the options (their key words) that can be set in the option string
///
/// format of function string:
///
///      "x0*(0)+((1)/x1)**(2)..."
///
/// where "[i]" are the parameters, and "xi" the input variables
///
/// format of parameter string:
///
///      "(-1.2,3.4);(-2.3,4.55);..."
///
/// where the numbers in "(a,b)" correspond to the a=min, b=max parameter ranges;
/// each parameter defined in the function string must have a corresponding range

void TMVA::MethodFDA::DeclareOptions()
{
   DeclareOptionRef( fFormulaStringP  = "(0)", "Formula",   "The discrimination formula" );
   DeclareOptionRef( fParRangeStringP = "()", "ParRanges", "Parameter ranges" );

   // fitter
   DeclareOptionRef( fFitMethod = "MINUIT", "FitMethod", "Optimisation Method");
   AddPreDefVal(TString("MC"));
   AddPreDefVal(TString("GA"));
   AddPreDefVal(TString("SA"));
   AddPreDefVal(TString("MINUIT"));

   DeclareOptionRef( fConverger = "None", "Converger", "FitMethod uses Converger to improve result");
   AddPreDefVal(TString("None"));
   AddPreDefVal(TString("MINUIT"));
}

////////////////////////////////////////////////////////////////////////////////
/// translate formula string into TFormula, and parameter string into par ranges

void TMVA::MethodFDA::CreateFormula()
{
   // process transient strings
   fFormulaStringT  = fFormulaStringP;

   // interpret formula string

   // replace the parameters "(i)" by the TFormula style "[i]"
   for (UInt_t ipar=0; ipar<fNPars; ipar++) {
      fFormulaStringT.ReplaceAll( Form("(%i)",ipar), Form("[%i]",ipar) );
   }

   // sanity check, there should be no "(i)", with 'i' a number anymore
   for (Int_t ipar=fNPars; ipar<1000; ipar++) {
      if (fFormulaStringT.Contains( Form("(%i)",ipar) ))
         Log() << kFATAL
               << "<CreateFormula> Formula contains expression: \"" << Form("(%i)",ipar) << "\", "
               << "which cannot be attributed to a parameter; "
               << "it may be that the number of variable ranges given via \"ParRanges\" "
               << "does not match the number of parameters in the formula expression, please verify!"
               << Endl;
   }

   // write the variables "xi" as additional parameters "[npar+i]"
   for (Int_t ivar=GetNvar()-1; ivar >= 0; ivar--) {
      fFormulaStringT.ReplaceAll( Form("x%i",ivar), Form("[%i]",ivar+fNPars) );
   }

   // sanity check, there should be no "xi", with 'i' a number anymore
   for (UInt_t ivar=GetNvar(); ivar<1000; ivar++) {
      if (fFormulaStringT.Contains( Form("x%i",ivar) ))
         Log() << kFATAL
               << "<CreateFormula> Formula contains expression: \"" << Form("x%i",ivar) << "\", "
               << "which cannot be attributed to an input variable" << Endl;
   }

   Log() << "User-defined formula string       : \"" << fFormulaStringP << "\"" << Endl;
   Log() << "TFormula-compatible formula string: \"" << fFormulaStringT << "\"" << Endl;
   Log() << kDEBUG << "Creating and compiling formula" << Endl;

   // create TF1
   if (fFormula) delete fFormula;
   fFormula = new TFormula( "FDA_Formula", fFormulaStringT );

   // is formula correct ?
   if (!fFormula->IsValid())
      Log() << kFATAL << "<ProcessOptions> Formula expression could not be properly compiled" << Endl;

   // other sanity checks
   if (fFormula->GetNpar() > (Int_t)(fNPars + GetNvar()))
      Log() << kFATAL << "<ProcessOptions> Dubious number of parameters in formula expression: "
            << fFormula->GetNpar() << " - compared to maximum allowed: " << fNPars + GetNvar() << Endl;
}

////////////////////////////////////////////////////////////////////////////////
/// the option string is decoded, for available options see "DeclareOptions"

void TMVA::MethodFDA::ProcessOptions()
{
   // process transient strings
   fParRangeStringT = fParRangeStringP;

   // interpret parameter string
   fParRangeStringT.ReplaceAll( " ", "" );
   fNPars = fParRangeStringT.CountChar( ')' );

   TList* parList = gTools().ParseFormatLine( fParRangeStringT, ";" );
   if ((UInt_t)parList->GetSize() != fNPars) {
      Log() << kFATAL << "<ProcessOptions> Mismatch in parameter string: "
            << "the number of parameters: " << fNPars << " != ranges defined: "
            << parList->GetSize() << "; the format of the \"ParRanges\" string "
            << "must be: \"(-1.2,3.4);(-2.3,4.55);...\", "
            << "where the numbers in \"(a,b)\" correspond to the a=min, b=max parameter ranges; "
            << "each parameter defined in the function string must have a corresponding rang."
            << Endl;
   }

   fParRange.resize( fNPars );
   for (UInt_t ipar=0; ipar<fNPars; ipar++) fParRange[ipar] = 0;

   for (UInt_t ipar=0; ipar<fNPars; ipar++) {
      // parse (a,b)
      TString str = ((TObjString*)parList->At(ipar))->GetString();
      Ssiz_t istr = str.First( ',' );
      TString pminS(str(1,istr-1));
      TString pmaxS(str(istr+1,str.Length()-2-istr));

      stringstream stmin; Float_t pmin=0; stmin << pminS.Data(); stmin >> pmin;
      stringstream stmax; Float_t pmax=0; stmax << pmaxS.Data(); stmax >> pmax;

      // sanity check
      if (TMath::Abs(pmax-pmin) < 1.e-30) pmax = pmin;
      if (pmin > pmax) Log() << kFATAL << "<ProcessOptions> max > min in interval for parameter: ["
                             << ipar << "] : [" << pmin  << ", " << pmax << "] " << Endl;

      Log() << kINFO << "Create parameter interval for parameter " << ipar << " : [" << pmin << "," << pmax << "]" << Endl;
      fParRange[ipar] = new Interval( pmin, pmax );
   }
   delete parList;

   // create formula
   CreateFormula();


   // copy parameter ranges for each output dimension ==================
   fOutputDimensions = 1;
   if( DoRegression() )
      fOutputDimensions = DataInfo().GetNTargets();
   if( DoMulticlass() )
      fOutputDimensions = DataInfo().GetNClasses();

   for( Int_t dim = 1; dim < fOutputDimensions; ++dim ){
      for( UInt_t par = 0; par < fNPars; ++par ){
         fParRange.push_back( fParRange.at(par) );
      }
   }
   // ====================

   // create minimiser
   fConvergerFitter = (IFitterTarget*)this;
   if (fConverger == "MINUIT") {
      fConvergerFitter = new MinuitFitter( *this, Form("%s_Converger_Minuit", GetName()), fParRange, GetOptions() );
      SetOptions(dynamic_cast<Configurable*>(fConvergerFitter)->GetOptions());
   }

   if(fFitMethod == "MC")
      fFitter = new MCFitter( *fConvergerFitter, Form("%s_Fitter_MC", GetName()), fParRange, GetOptions() );
   else if (fFitMethod == "GA")
      fFitter = new GeneticFitter( *fConvergerFitter, Form("%s_Fitter_GA", GetName()), fParRange, GetOptions() );
   else if (fFitMethod == "SA")
      fFitter = new SimulatedAnnealingFitter( *fConvergerFitter, Form("%s_Fitter_SA", GetName()), fParRange, GetOptions() );
   else if (fFitMethod == "MINUIT")
      fFitter = new MinuitFitter( *fConvergerFitter, Form("%s_Fitter_Minuit", GetName()), fParRange, GetOptions() );
   else {
      Log() << kFATAL << "<Train> Do not understand fit method:" << fFitMethod << Endl;
   }

   fFitter->CheckForUnusedOptions();
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TMVA::MethodFDA::~MethodFDA( void )
{
   ClearAll();
}

////////////////////////////////////////////////////////////////////////////////
/// FDA can handle classification with 2 classes and regression with one regression-target

Bool_t TMVA::MethodFDA::HasAnalysisType( Types::EAnalysisType type, UInt_t numberClasses, UInt_t /*numberTargets*/ )
{
   if (type == Types::kClassification && numberClasses == 2) return kTRUE;
   if (type == Types::kMulticlass ) return kTRUE;
   if (type == Types::kRegression ) return kTRUE;
   return kFALSE;
}


////////////////////////////////////////////////////////////////////////////////
/// delete and clear all class members

void TMVA::MethodFDA::ClearAll( void )
{
   // if there is more than one output dimension, the paramater ranges are the same again (object has been copied).
   // hence, ... erase the copied pointers to assure, that they are deleted only once.
   //   fParRange.erase( fParRange.begin()+(fNPars), fParRange.end() );
   for (UInt_t ipar=0; ipar<fParRange.size() && ipar<fNPars; ipar++) {
      if (fParRange[ipar] != 0) { delete fParRange[ipar]; fParRange[ipar] = 0; }
   }
   fParRange.clear();

   if (fFormula  != 0) { delete fFormula; fFormula = 0; }
   fBestPars.clear();
}

////////////////////////////////////////////////////////////////////////////////
/// FDA training

void TMVA::MethodFDA::Train( void )
{
   // cache training events
   fSumOfWeights    = 0;
   fSumOfWeightsSig = 0;
   fSumOfWeightsBkg = 0;

   for (UInt_t ievt=0; ievt<GetNEvents(); ievt++) {

      // read the training event
      const Event* ev = GetEvent(ievt);

      // true event copy
      Float_t w  = ev->GetWeight();

      if (!DoRegression()) {
         if (DataInfo().IsSignal(ev)) { fSumOfWeightsSig += w; }
         else                { fSumOfWeightsBkg += w; }
      }
      fSumOfWeights += w;
   }

   // sanity check
   if (!DoRegression()) {
      if (fSumOfWeightsSig <= 0 || fSumOfWeightsBkg <= 0) {
         Log() << kFATAL << "<Train> Troubles in sum of weights: "
               << fSumOfWeightsSig << " (S) : " << fSumOfWeightsBkg << " (B)" << Endl;
      }
   }
   else if (fSumOfWeights <= 0) {
      Log() << kFATAL << "<Train> Troubles in sum of weights: "
            << fSumOfWeights << Endl;
   }

   // starting values (not used by all fitters)
   fBestPars.clear();
   for (std::vector<Interval*>::const_iterator parIt = fParRange.begin(); parIt != fParRange.end(); ++parIt) {
      fBestPars.push_back( (*parIt)->GetMean() );
   }

   // execute the fit
   Double_t estimator = fFitter->Run( fBestPars );

   // print results
   PrintResults( fFitMethod, fBestPars, estimator );

   delete fFitter; fFitter = 0;
   if (fConvergerFitter!=0 && fConvergerFitter!=(IFitterTarget*)this) {
      delete fConvergerFitter;
      fConvergerFitter = 0;
   }
   ExitFromTraining();
}

////////////////////////////////////////////////////////////////////////////////
/// display fit parameters
/// check maximum length of variable name

void TMVA::MethodFDA::PrintResults( const TString& fitter, std::vector<Double_t>& pars, const Double_t estimator ) const
{
   Log() << kINFO;
   Log() << kHEADER << "Results for parameter fit using \"" << fitter << "\" fitter:" << Endl;
   std::vector<TString>  parNames;
   for (UInt_t ipar=0; ipar<pars.size(); ipar++) parNames.push_back( Form("Par(%i)",ipar ) );
   gTools().FormattedOutput( pars, parNames, "Parameter" , "Fit result", Log(), "%g" );
   Log() << "Discriminator expression: \"" << fFormulaStringP << "\"" << Endl;
   Log() << "Value of estimator at minimum: " << estimator << Endl;
}

////////////////////////////////////////////////////////////////////////////////
/// compute estimator for given parameter set (to be minimised)

Double_t TMVA::MethodFDA::EstimatorFunction( std::vector<Double_t>& pars )
{
   const Double_t sumOfWeights[]                = { fSumOfWeightsBkg, fSumOfWeightsSig, fSumOfWeights };
   Double_t estimator[]                         = { 0, 0, 0 };

   Double_t result, deviation;
   Double_t desired = 0.0;

   // calculate the deviation from the desired value
   if( DoRegression() ){
      for (UInt_t ievt=0; ievt<GetNEvents(); ievt++) {
         // read the training event
         const TMVA::Event* ev = GetEvent(ievt);

         for( Int_t dim = 0; dim < fOutputDimensions; ++dim ){
            desired = ev->GetTarget( dim );
            result    = InterpretFormula( ev, pars.begin(), pars.end() );
            deviation = TMath::Power(result - desired, 2);
            estimator[2]  += deviation * ev->GetWeight();
         }
      }
      estimator[2] /= sumOfWeights[2];
      // return value is sum over normalised signal and background contributions
      return estimator[2];

   }else if( DoMulticlass() ){
      for (UInt_t ievt=0; ievt<GetNEvents(); ievt++) {
         // read the training event
         const TMVA::Event* ev = GetEvent(ievt);

         CalculateMulticlassValues( ev, pars, *fMulticlassReturnVal );

         Double_t crossEntropy = 0.0;
         for( Int_t dim = 0; dim < fOutputDimensions; ++dim ){
            Double_t y = fMulticlassReturnVal->at(dim);
            Double_t t = (ev->GetClass() == static_cast<UInt_t>(dim) ? 1.0 : 0.0 );
            crossEntropy += t*log(y);
         }
         estimator[2] += ev->GetWeight()*crossEntropy;
      }
      estimator[2] /= sumOfWeights[2];
      // return value is sum over normalised signal and background contributions
      return estimator[2];

   }else{
      for (UInt_t ievt=0; ievt<GetNEvents(); ievt++) {
         // read the training event
         const TMVA::Event* ev = GetEvent(ievt);

         desired = (DataInfo().IsSignal(ev) ? 1.0 : 0.0);
         result    = InterpretFormula( ev, pars.begin(), pars.end() );
         deviation = TMath::Power(result - desired, 2);
         estimator[Int_t(desired)] += deviation * ev->GetWeight();
      }
      estimator[0] /= sumOfWeights[0];
      estimator[1] /= sumOfWeights[1];
      // return value is sum over normalised signal and background contributions
      return estimator[0] + estimator[1];
   }
}

////////////////////////////////////////////////////////////////////////////////
/// formula interpretation

Double_t TMVA::MethodFDA::InterpretFormula( const Event* event, std::vector<Double_t>::iterator parBegin, std::vector<Double_t>::iterator parEnd )
{
   Int_t ipar = 0;
   //    std::cout << "pars ";
   for( std::vector<Double_t>::iterator it = parBegin; it != parEnd; ++it ){
      //       std::cout << " i" << ipar << " val" << (*it);
      fFormula->SetParameter( ipar, (*it) );
      ++ipar;
   }
   for (UInt_t ivar=0;  ivar<GetNvar();  ivar++) fFormula->SetParameter( ivar+ipar, event->GetValue(ivar) );

   Double_t result = fFormula->Eval( 0 );
   //    std::cout << "  result " << result << std::endl;
   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// returns MVA value for given event

Double_t TMVA::MethodFDA::GetMvaValue( Double_t* err, Double_t* errUpper )
{
   const Event* ev = GetEvent();

   // cannot determine error
   NoErrorCalc(err, errUpper);

   return InterpretFormula( ev, fBestPars.begin(), fBestPars.end() );
}

////////////////////////////////////////////////////////////////////////////////

const std::vector<Float_t>& TMVA::MethodFDA::GetRegressionValues()
{
   if (fRegressionReturnVal == NULL) fRegressionReturnVal = new std::vector<Float_t>();
   fRegressionReturnVal->clear();

   const Event* ev = GetEvent();

   Event* evT = new Event(*ev);

   for( Int_t dim = 0; dim < fOutputDimensions; ++dim ){
      Int_t offset = dim*fNPars;
      evT->SetTarget(dim,InterpretFormula( ev, fBestPars.begin()+offset, fBestPars.begin()+offset+fNPars ) );
   }
   const Event* evT2 = GetTransformationHandler().InverseTransform( evT );
   fRegressionReturnVal->push_back(evT2->GetTarget(0));

   delete evT;

   return (*fRegressionReturnVal);
}

////////////////////////////////////////////////////////////////////////////////

const std::vector<Float_t>& TMVA::MethodFDA::GetMulticlassValues()
{
   if (fMulticlassReturnVal == NULL) fMulticlassReturnVal = new std::vector<Float_t>();
   fMulticlassReturnVal->clear();
   std::vector<Float_t> temp;

   // returns MVA value for given event
   const TMVA::Event* evt = GetEvent();

   CalculateMulticlassValues( evt, fBestPars, temp );

   UInt_t nClasses = DataInfo().GetNClasses();
   for(UInt_t iClass=0; iClass<nClasses; iClass++){
      Double_t norm = 0.0;
      for(UInt_t j=0;j<nClasses;j++){
         if(iClass!=j)
            norm+=exp(temp[j]-temp[iClass]);
      }
      (*fMulticlassReturnVal).push_back(1.0/(1.0+norm));
   }

   return (*fMulticlassReturnVal);
}


////////////////////////////////////////////////////////////////////////////////
/// calculate the values for multiclass

void TMVA::MethodFDA::CalculateMulticlassValues( const TMVA::Event*& evt, std::vector<Double_t>& parameters, std::vector<Float_t>& values)
{
   values.clear();

   //    std::copy( parameters.begin(), parameters.end(), std::ostream_iterator<double>( std::cout, " " ) );
   //    std::cout << std::endl;

   //    char inp;
   //    std::cin >> inp;

   Double_t sum=0;
   for( Int_t dim = 0; dim < fOutputDimensions; ++dim ){ // check for all other dimensions (=classes)
      Int_t offset = dim*fNPars;
      Double_t value = InterpretFormula( evt, parameters.begin()+offset, parameters.begin()+offset+fNPars );
      //       std::cout << "dim : " << dim << " value " << value << "    offset " << offset << std::endl;
      values.push_back( value );
      sum += value;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// read back the training results from a file (stream)

void  TMVA::MethodFDA::ReadWeightsFromStream( std::istream& istr )
{
   // retrieve best function parameters
   // coverity[tainted_data_argument]
   istr >> fNPars;

   fBestPars.clear();
   fBestPars.resize( fNPars );
   for (UInt_t ipar=0; ipar<fNPars; ipar++) istr >> fBestPars[ipar];
}

////////////////////////////////////////////////////////////////////////////////
/// create XML description for LD classification and regression
/// (for arbitrary number of output classes/targets)

void TMVA::MethodFDA::AddWeightsXMLTo( void* parent ) const
{
   void* wght = gTools().AddChild(parent, "Weights");
   gTools().AddAttr( wght, "NPars",  fNPars );
   gTools().AddAttr( wght, "NDim",   fOutputDimensions );
   for (UInt_t ipar=0; ipar<fNPars*fOutputDimensions; ipar++) {
      void* coeffxml = gTools().AddChild( wght, "Parameter" );
      gTools().AddAttr( coeffxml, "Index", ipar   );
      gTools().AddAttr( coeffxml, "Value", fBestPars[ipar] );
   }

   // write formula
   gTools().AddAttr( wght, "Formula", fFormulaStringP );
}

////////////////////////////////////////////////////////////////////////////////
/// read coefficients from xml weight file

void TMVA::MethodFDA::ReadWeightsFromXML( void* wghtnode )
{
   gTools().ReadAttr( wghtnode, "NPars", fNPars );

   if(gTools().HasAttr( wghtnode, "NDim")) {
      gTools().ReadAttr( wghtnode, "NDim" , fOutputDimensions );
   } else {
      // older weight files don't have this attribute
      fOutputDimensions = 1;
   }

   fBestPars.clear();
   fBestPars.resize( fNPars*fOutputDimensions );

   void* ch = gTools().GetChild(wghtnode);
   Double_t par;
   UInt_t    ipar;
   while (ch) {
      gTools().ReadAttr( ch, "Index", ipar );
      gTools().ReadAttr( ch, "Value", par  );

      // sanity check
      if (ipar >= fNPars*fOutputDimensions) Log() << kFATAL << "<ReadWeightsFromXML> index out of range: "
                                                  << ipar << " >= " << fNPars << Endl;
      fBestPars[ipar] = par;

      ch = gTools().GetNextChild(ch);
   }

   // read formula
   gTools().ReadAttr( wghtnode, "Formula", fFormulaStringP );

   // create the TFormula
   CreateFormula();
}

////////////////////////////////////////////////////////////////////////////////
/// write FDA-specific classifier response

void TMVA::MethodFDA::MakeClassSpecific( std::ostream& fout, const TString& className ) const
{
   fout << "   double              fParameter[" << fNPars << "];" << std::endl;
   fout << "};" << std::endl;
   fout << "" << std::endl;
   fout << "inline void " << className << "::Initialize() " << std::endl;
   fout << "{" << std::endl;
   for(UInt_t ipar=0; ipar<fNPars; ipar++) {
      fout << "   fParameter[" << ipar << "] = " << fBestPars[ipar] << ";" << std::endl;
   }
   fout << "}" << std::endl;
   fout << std::endl;
   fout << "inline double " << className << "::GetMvaValue__( const std::vector<double>& inputValues ) const" << std::endl;
   fout << "{" << std::endl;
   fout << "   // interpret the formula" << std::endl;

   // replace parameters
   TString str = fFormulaStringT;
   for (UInt_t ipar=0; ipar<fNPars; ipar++) {
      str.ReplaceAll( Form("[%i]", ipar), Form("fParameter[%i]", ipar) );
   }

   // replace input variables
   for (UInt_t ivar=0; ivar<GetNvar(); ivar++) {
      str.ReplaceAll( Form("[%i]", ivar+fNPars), Form("inputValues[%i]", ivar) );
   }

   fout << "   double retval = " << str << ";" << std::endl;
   fout << std::endl;
   fout << "   return retval; " << std::endl;
   fout << "}" << std::endl;
   fout << std::endl;
   fout << "// Clean up" << std::endl;
   fout << "inline void " << className << "::Clear() " << std::endl;
   fout << "{" << std::endl;
   fout << "   // nothing to clear" << std::endl;
   fout << "}" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// get help message text
///
/// typical length of text line:
///         "|--------------------------------------------------------------|"

void TMVA::MethodFDA::GetHelpMessage() const
{
   Log() << Endl;
   Log() << gTools().Color("bold") << "--- Short description:" << gTools().Color("reset") << Endl;
   Log() << Endl;
   Log() << "The function discriminant analysis (FDA) is a classifier suitable " << Endl;
   Log() << "to solve linear or simple nonlinear discrimination problems." << Endl;
   Log() << Endl;
   Log() << "The user provides the desired function with adjustable parameters" << Endl;
   Log() << "via the configuration option string, and FDA fits the parameters to" << Endl;
   Log() << "it, requiring the signal (background) function value to be as close" << Endl;
   Log() << "as possible to 1 (0). Its advantage over the more involved and" << Endl;
   Log() << "automatic nonlinear discriminators is the simplicity and transparency " << Endl;
   Log() << "of the discrimination expression. A shortcoming is that FDA will" << Endl;
   Log() << "underperform for involved problems with complicated, phase space" << Endl;
   Log() << "dependent nonlinear correlations." << Endl;
   Log() << Endl;
   Log() << "Please consult the Users Guide for the format of the formula string" << Endl;
   Log() << "and the allowed parameter ranges:" << Endl;
   if (gConfig().WriteOptionsReference()) {
      Log() << "<a href=\"http://tmva.sourceforge.net/docu/TMVAUsersGuide.pdf\">"
            << "http://tmva.sourceforge.net/docu/TMVAUsersGuide.pdf</a>" << Endl;
   }
   else Log() << "http://tmva.sourceforge.net/docu/TMVAUsersGuide.pdf" << Endl;
   Log() << Endl;
   Log() << gTools().Color("bold") << "--- Performance optimisation:" << gTools().Color("reset") << Endl;
   Log() << Endl;
   Log() << "The FDA performance depends on the complexity and fidelity of the" << Endl;
   Log() << "user-defined discriminator function. As a general rule, it should" << Endl;
   Log() << "be able to reproduce the discrimination power of any linear" << Endl;
   Log() << "discriminant analysis. To reach into the nonlinear domain, it is" << Endl;
   Log() << "useful to inspect the correlation profiles of the input variables," << Endl;
   Log() << "and add quadratic and higher polynomial terms between variables as" << Endl;
   Log() << "necessary. Comparison with more involved nonlinear classifiers can" << Endl;
   Log() << "be used as a guide." << Endl;
   Log() << Endl;
   Log() << gTools().Color("bold") << "--- Performance tuning via configuration options:" << gTools().Color("reset") << Endl;
   Log() << Endl;
   Log() << "Depending on the function used, the choice of \"FitMethod\" is" << Endl;
   Log() << "crucial for getting valuable solutions with FDA. As a guideline it" << Endl;
   Log() << "is recommended to start with \"FitMethod=MINUIT\". When more complex" << Endl;
   Log() << "functions are used where MINUIT does not converge to reasonable" << Endl;
   Log() << "results, the user should switch to non-gradient FitMethods such" << Endl;
   Log() << "as GeneticAlgorithm (GA) or Monte Carlo (MC). It might prove to be" << Endl;
   Log() << "useful to combine GA (or MC) with MINUIT by setting the option" << Endl;
   Log() << "\"Converger=MINUIT\". GA (MC) will then set the starting parameters" << Endl;
   Log() << "for MINUIT such that the basic quality of GA (MC) of finding global" << Endl;
   Log() << "minima is combined with the efficacy of MINUIT of finding local" << Endl;
   Log() << "minima." << Endl;
}
