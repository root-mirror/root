// @(#)root/tmva $Id$
// Author: Marcin Wolter, Andrzej Zemla

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : MethodSVM                                                             *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Implementation                                                            *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Marcin Wolter  <Marcin.Wolter@cern.ch> - IFJ PAN, Krakow, Poland          *
 *      Andrzej Zemla  <azemla@cern.ch>          - IFJ PAN, Krakow, Poland        *
 *      (IFJ PAN: Henryk Niewodniczanski Inst. Nucl. Physics, Krakow, Poland)     *
 *                                                                                *
 * Introduction of regression by:                                                 *
 *      Krzysztof Danielowski <danielow@cern.ch> - IFJ PAN & AGH, Krakow, Poland  *
 *      Kamil Kraszewski      <kalq@cern.ch>     - IFJ PAN & UJ, Krakow, Poland   *
 *      Maciej Kruk           <mkruk@cern.ch>    - IFJ PAN & AGH, Krakow, Poland  *
 *                                                                                *
 * Introduction of kernel parameter optimisation                                  *
 *            and additional kernel functions by:                                 *
 *      Adrian Bevan          <adrian.bevan@cern.ch> -   Queen Mary               *
 *                                                       University of London, UK *
 *      Tom Stevenson <thomas.james.stevenson@cern.ch> - Queen Mary               *
 *                                                       University of London, UK *
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         *
 *      MPI-K Heidelberg, Germany                                                 *
 *      PAN, Krakow, Poland                                                       *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

/*! \class TMVA::MethodSVM
\ingroup TMVA
SMO Platt's SVM classifier with Keerthi & Shavade improvements
*/

#include "TMVA/MethodSVM.h"

#include "TMVA/Tools.h"
#include "TMVA/Timer.h"

#include "TMVA/SVWorkingSet.h"

#include "TMVA/SVEvent.h"

#include "TMVA/SVKernelFunction.h"

#include "TMVA/ClassifierFactory.h"
#include "TMVA/Configurable.h"
#include "TMVA/DataSet.h"
#include "TMVA/DataSetInfo.h"
#include "TMVA/Event.h"
#include "TMVA/IMethod.h"
#include "TMVA/MethodBase.h"
#include "TMVA/MsgLogger.h"
#include "TMVA/Types.h"
#include "TMVA/Interval.h"
#include "TMVA/OptimizeConfigParameters.h"
#include "TMVA/Results.h"
#include "TMVA/ResultsClassification.h"
#include "TMVA/VariableInfo.h"

#include "TFile.h"
#include "TVectorD.h"
#include "TMath.h"

#include <iostream>
#include <string>

using std::vector;
using std::string;
using std::stringstream;

//const Int_t basketsize__ = 1280000;
REGISTER_METHOD(SVM)

ClassImp(TMVA::MethodSVM);

////////////////////////////////////////////////////////////////////////////////
/// standard constructor

   TMVA::MethodSVM::MethodSVM( const TString& jobName, const TString& methodTitle, DataSetInfo& theData,
                               const TString& theOption )
   : MethodBase( jobName, Types::kSVM, methodTitle, theData, theOption)
   , fCost(0)
   , fTolerance(0)
   , fMaxIter(0)
   , fNSubSets(0)
   , fBparm(0)
   , fGamma(0)
   , fWgSet(0)
   , fInputData(0)
   , fSupportVectors(0)
   , fSVKernelFunction(0)
   , fMinVars(0)
   , fMaxVars(0)
   , fDoubleSigmaSquared(0)
   , fOrder(0)
   , fTheta(0)
   , fKappa(0)
   , fMult(0)
   ,fNumVars(0)
   , fGammas("")
   , fGammaList("")
   , fDataSize(0)
   , fLoss(0)
{
   fVarNames.clear();
   fNumVars = theData.GetVariableInfos().size();
   for( int i=0; i<fNumVars; i++){
      fVarNames.push_back(theData.GetVariableInfos().at(i).GetTitle());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// constructor from weight file

TMVA::MethodSVM::MethodSVM( DataSetInfo& theData, const TString& theWeightFile)
   : MethodBase( Types::kSVM, theData, theWeightFile)
   , fCost(0)
   , fTolerance(0)
   , fMaxIter(0)
   , fNSubSets(0)
   , fBparm(0)
   , fGamma(0)
   , fWgSet(0)
   , fInputData(0)
   , fSupportVectors(0)
   , fSVKernelFunction(0)
   , fMinVars(0)
   , fMaxVars(0)
   , fDoubleSigmaSquared(0)
   , fOrder(0)
   , fTheta(0)
   , fKappa(0)
   , fMult(0)
   , fNumVars(0)
   , fGammas("")
   , fGammaList("")
   , fDataSize(0)
   , fLoss(0)
{
   fVarNames.clear();
   fNumVars = theData.GetVariableInfos().size();
   for( int i=0;i<fNumVars; i++){
      fVarNames.push_back(theData.GetVariableInfos().at(i).GetTitle());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TMVA::MethodSVM::~MethodSVM()
{
   fSupportVectors->clear();
   for (UInt_t i=0; i<fInputData->size(); i++) {
      delete fInputData->at(i);
   }
   if (fWgSet !=0)           { delete fWgSet; fWgSet=0; }
   if (fSVKernelFunction !=0 ) { delete fSVKernelFunction; fSVKernelFunction = 0; }
}

////////////////////////////////////////////////////////////////////////////////
// reset the method, as if it had just been instantiated (forget all training etc.)

void TMVA::MethodSVM::Reset( void )
{
   // reset the method, as if it had just been instantiated (forget all training etc.)
   fSupportVectors->clear();
   for (UInt_t i=0; i<fInputData->size(); i++){
      delete fInputData->at(i);
      fInputData->at(i)=0;
   }
   fInputData->clear();
   if (fWgSet !=0)           { fWgSet=0; }
   if (fSVKernelFunction !=0 ) { fSVKernelFunction = 0; }
   if (Data()){
      Data()->DeleteResults(GetMethodName(), Types::kTraining, GetAnalysisType());
   }

   Log() << kDEBUG << " successfully(?) reset the method " << Endl;
}

////////////////////////////////////////////////////////////////////////////////
/// SVM can handle classification with 2 classes and regression with one regression-target

Bool_t TMVA::MethodSVM::HasAnalysisType( Types::EAnalysisType type, UInt_t numberClasses, UInt_t numberTargets )
{
   if (type == Types::kClassification && numberClasses == 2) return kTRUE;
   if (type == Types::kRegression     && numberTargets == 1) return kTRUE;
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// default initialisation

void TMVA::MethodSVM::Init()
{
   // SVM always uses normalised input variables
   SetNormalised( kTRUE );

   // Helge: do not book a event vector of given size but rather fill the vector
   //        later with pus_back. Anyway, this is NOT what is time consuming in
   //        SVM and it allows to skip totally events with weights == 0 ;)
   fInputData = new std::vector<TMVA::SVEvent*>(0);
   fSupportVectors = new std::vector<TMVA::SVEvent*>(0);
}

////////////////////////////////////////////////////////////////////////////////
/// declare options available for this method

void TMVA::MethodSVM::DeclareOptions()
{
   DeclareOptionRef( fTheKernel = "RBF", "Kernel", "Pick which kernel ( RBF or MultiGauss )");
   // for gaussian kernel parameter(s)
   DeclareOptionRef( fGamma = 1., "Gamma", "RBF kernel parameter: Gamma (size of the Kernel)");
   // for polynomial kernel parameter(s)
   DeclareOptionRef( fOrder = 3, "Order", "Polynomial Kernel parameter: polynomial order");
   DeclareOptionRef( fTheta = 1., "Theta", "Polynomial Kernel parameter: polynomial theta");
   // for multi-gaussian kernel parameter(s)
   DeclareOptionRef( fGammas = "", "GammaList", "MultiGauss parameters" );

   // for range and step number for kernel parameter optimisation
   DeclareOptionRef( fTune = "All", "Tune", "Tune Parameters");
   // for list of kernels to be used with product or sum kernel
   DeclareOptionRef( fMultiKernels = "None", "KernelList", "Sum or product of kernels");
   DeclareOptionRef( fLoss = "hinge", "Loss", "Loss function");

   DeclareOptionRef( fCost,   "C",        "Cost parameter" );
   if (DoRegression()) {
      fCost = 0.002;
   }else{
      fCost = 1.;
   }
   DeclareOptionRef( fTolerance = 0.01, "Tol",      "Tolerance parameter" );  //should be fixed
   DeclareOptionRef( fMaxIter   = 1000, "MaxIter",  "Maximum number of training loops" );

}

////////////////////////////////////////////////////////////////////////////////
/// options that are used ONLY for the READER to ensure backward compatibility

void TMVA::MethodSVM::DeclareCompatibilityOptions()
{
   MethodBase::DeclareCompatibilityOptions();
   DeclareOptionRef( fNSubSets  = 1,    "NSubSets", "Number of training subsets" );
   DeclareOptionRef( fTheKernel = "Gauss", "Kernel", "Uses kernel function");
   // for gaussian kernel parameter(s)
   DeclareOptionRef( fDoubleSigmaSquared = 2., "Sigma", "Kernel parameter: sigma");
   // for polynomial kernel parameter(s)
   DeclareOptionRef( fOrder = 3, "Order", "Polynomial Kernel parameter: polynomial order");
   // for sigmoid kernel parameters
   DeclareOptionRef( fTheta = 1., "Theta", "Sigmoid Kernel parameter: theta");
   DeclareOptionRef( fKappa = 1., "Kappa", "Sigmoid Kernel parameter: kappa");
}

////////////////////////////////////////////////////////////////////////////////
/// option post processing (if necessary)

void TMVA::MethodSVM::ProcessOptions()
{
   if (IgnoreEventsWithNegWeightsInTraining()) {
      Log() << kFATAL << "Mechanism to ignore events with negative weights in training not yet available for method: "
            << GetMethodTypeName()
            << " --> please remove \"IgnoreNegWeightsInTraining\" option from booking string."
            << Endl;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Train SVM

void TMVA::MethodSVM::Train()
{
   fIPyMaxIter = fMaxIter;
   Data()->SetCurrentType(Types::kTraining);

   Log() << kDEBUG << "Create event vector"<< Endl;

   fDataSize = Data()->GetNEvents();
   Int_t nSignal = Data()->GetNEvtSigTrain();
   Int_t nBackground = Data()->GetNEvtBkgdTrain();
   Double_t CSig;
   Double_t CBkg;

   // Use number of signal and background from above to weight the cost parameter
   // so that the training is not biased towards the larger dataset when the signal
   // and background samples are significantly different sizes.
   if(nSignal < nBackground){
      CSig = fCost;
      CBkg = CSig*((double)nSignal/nBackground);
   }
   else{
      CBkg = fCost;
      CSig = CBkg*((double)nSignal/nBackground);
   }

   // Loop over events and assign the correct cost parameter.
   for (Int_t ievnt=0; ievnt<Data()->GetNEvents(); ievnt++){
      if (GetEvent(ievnt)->GetWeight() != 0){
         if(DataInfo().IsSignal(GetEvent(ievnt))){
            fInputData->push_back(new SVEvent(GetEvent(ievnt), CSig, DataInfo().IsSignal\
                                              (GetEvent(ievnt))));
         }
         else{
            fInputData->push_back(new SVEvent(GetEvent(ievnt), CBkg, DataInfo().IsSignal\
                                              (GetEvent(ievnt))));
         }
      }
   }

   // Set the correct kernel function.
   // Here we only use valid Mercer kernels. In the literature some people have reported reasonable
   // results using Sigmoid kernel function however that is not a valid Mercer kernel and is not used here.
   if( fTheKernel == "RBF"){
      fSVKernelFunction = new SVKernelFunction( SVKernelFunction::kRBF, fGamma);
   }
   else if( fTheKernel == "MultiGauss" ){
      if(fGammas!=""){
         SetMGamma(fGammas);
         fGammaList=fGammas;
      }
      else{
         if(fmGamma.size()!=0){ GetMGamma(fmGamma); } // Set fGammas if empty to write to XML file
         else{
            for(Int_t ngammas=0; ngammas<fNumVars; ++ngammas){
               fmGamma.push_back(1.0);
            }
            GetMGamma(fmGamma);
         }
      }
      fSVKernelFunction = new SVKernelFunction(fmGamma);
   }
   else if( fTheKernel == "Polynomial" ){
      fSVKernelFunction = new SVKernelFunction( SVKernelFunction::kPolynomial, fOrder,fTheta);
   }
   else if( fTheKernel == "Prod" ){
      if(fGammas!=""){
         SetMGamma(fGammas);
         fGammaList=fGammas;
      }
      else{
         if(fmGamma.size()!=0){ GetMGamma(fmGamma); } // Set fGammas if empty to write to XML file
      }
      fSVKernelFunction = new SVKernelFunction( SVKernelFunction::kProd, MakeKernelList(fMultiKernels,fTheKernel), fmGamma, fGamma, fOrder, fTheta );
   }
   else if( fTheKernel == "Sum" ){
      if(fGammas!=""){
         SetMGamma(fGammas);
         fGammaList=fGammas;
      }
      else{
         if(fmGamma.size()!=0){ GetMGamma(fmGamma); } // Set fGammas if empty to write to XML file
      }
      fSVKernelFunction = new SVKernelFunction( SVKernelFunction::kSum, MakeKernelList(fMultiKernels,fTheKernel), fmGamma, fGamma, fOrder, fTheta );
   }
   else {
      Log() << kWARNING << fTheKernel << " is not a recognised kernel function." << Endl;
      exit(1);
   }

   Log()<< kINFO << "Building SVM Working Set...with "<<fInputData->size()<<" event instances"<< Endl;
   Timer bldwstime( GetName());
   fWgSet = new SVWorkingSet( fInputData, fSVKernelFunction,fTolerance, DoRegression() );
   Log() << kINFO <<"Elapsed time for Working Set build: "<< bldwstime.GetElapsedTime()<<Endl;

   // timing
   Timer timer( GetName() );
   Log() << kINFO << "Sorry, no computing time forecast available for SVM, please wait ..." << Endl;

   if (fInteractive) fWgSet->SetIPythonInteractive(&fExitFromTraining, &fIPyCurrentIter);

   fWgSet->Train(fMaxIter);

   Log() << kINFO << "Elapsed time: " << timer.GetElapsedTime()
         << "                                          " << Endl;

   fBparm          = fWgSet->GetBpar();
   fSupportVectors = fWgSet->GetSupportVectors();
   delete fWgSet;
   fWgSet=0;

    if (!fExitFromTraining) fIPyMaxIter = fIPyCurrentIter;
    ExitFromTraining();
}

////////////////////////////////////////////////////////////////////////////////
/// write configuration to xml file

void TMVA::MethodSVM::AddWeightsXMLTo( void* parent ) const
{
   void* wght = gTools().AddChild(parent, "Weights");
   gTools().AddAttr(wght,"fBparm",fBparm);
   gTools().AddAttr(wght,"fGamma",fGamma);
   gTools().AddAttr(wght,"fGammaList",fGammaList);
   gTools().AddAttr(wght,"fTheta",fTheta);
   gTools().AddAttr(wght,"fOrder",fOrder);
   gTools().AddAttr(wght,"NSupVec",fSupportVectors->size());

   for (std::vector<TMVA::SVEvent*>::iterator veciter=fSupportVectors->begin();
        veciter!=fSupportVectors->end() ; ++veciter ) {
      TVectorD temp(GetNvar()+4);
      temp[0] = (*veciter)->GetNs();
      temp[1] = (*veciter)->GetTypeFlag();
      temp[2] = (*veciter)->GetAlpha();
      temp[3] = (*veciter)->GetAlpha_p();
      for (UInt_t ivar = 0; ivar < GetNvar(); ivar++)
         temp[ivar+4] = (*(*veciter)->GetDataVector())[ivar];
      gTools().WriteTVectorDToXML(wght,"SupportVector",&temp);
   }
   // write max/min data values
   void* maxnode = gTools().AddChild(wght, "Maxima");
   for (UInt_t ivar = 0; ivar < GetNvar(); ivar++)
      gTools().AddAttr(maxnode, "Var"+gTools().StringFromInt(ivar), GetXmax(ivar));
   void* minnode = gTools().AddChild(wght, "Minima");
   for (UInt_t ivar = 0; ivar < GetNvar(); ivar++)
      gTools().AddAttr(minnode, "Var"+gTools().StringFromInt(ivar), GetXmin(ivar));
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::MethodSVM::ReadWeightsFromXML( void* wghtnode )
{
   gTools().ReadAttr( wghtnode, "fBparm",fBparm   );
   gTools().ReadAttr( wghtnode, "fGamma",fGamma);
   gTools().ReadAttr( wghtnode, "fGammaList",fGammaList);
   gTools().ReadAttr( wghtnode, "fOrder",fOrder);
   gTools().ReadAttr( wghtnode, "fTheta",fTheta);
   UInt_t fNsupv=0;
   gTools().ReadAttr( wghtnode, "NSupVec",fNsupv   );

   Float_t alpha=0.;
   Float_t alpha_p = 0.;

   Int_t typeFlag=-1;
   // UInt_t ns = 0;
   std::vector<Float_t>* svector = new std::vector<Float_t>(GetNvar());

   if (fMaxVars!=0) delete fMaxVars;
   fMaxVars = new TVectorD( GetNvar() );
   if (fMinVars!=0) delete fMinVars;
   fMinVars = new TVectorD( GetNvar() );
   if (fSupportVectors!=0) {
      for (vector< SVEvent* >::iterator it = fSupportVectors->begin(); it!=fSupportVectors->end(); ++it)
         delete *it;
      delete fSupportVectors;
   }
   fSupportVectors = new std::vector<TMVA::SVEvent*>(0);
   void* supportvectornode = gTools().GetChild(wghtnode);
   for (UInt_t ievt = 0; ievt < fNsupv; ievt++) {
      TVectorD temp(GetNvar()+4);
      gTools().ReadTVectorDFromXML(supportvectornode,"SupportVector",&temp);
      // ns=(UInt_t)temp[0];
      typeFlag=(int)temp[1];
      alpha=temp[2];
      alpha_p=temp[3];
      for (UInt_t ivar = 0; ivar < GetNvar(); ivar++) (*svector)[ivar]=temp[ivar+4];

      fSupportVectors->push_back(new SVEvent(svector,alpha,alpha_p,typeFlag));
      supportvectornode = gTools().GetNextChild(supportvectornode);
   }

   void* maxminnode = supportvectornode;
   for (UInt_t ivar = 0; ivar < GetNvar(); ivar++)
      gTools().ReadAttr( maxminnode,"Var"+gTools().StringFromInt(ivar),(*fMaxVars)[ivar]);
   maxminnode = gTools().GetNextChild(maxminnode);
   for (UInt_t ivar = 0; ivar < GetNvar(); ivar++)
      gTools().ReadAttr( maxminnode,"Var"+gTools().StringFromInt(ivar),(*fMinVars)[ivar]);
   if (fSVKernelFunction!=0) delete fSVKernelFunction;
   if( fTheKernel == "RBF" ){
      fSVKernelFunction = new SVKernelFunction(SVKernelFunction::kRBF, fGamma);
   }
   else if( fTheKernel == "MultiGauss" ){
      SetMGamma(fGammaList);
      fSVKernelFunction = new SVKernelFunction(fmGamma);
   }
   else if( fTheKernel == "Polynomial" ){
      fSVKernelFunction = new SVKernelFunction(SVKernelFunction::kPolynomial, fOrder, fTheta);
   }
   else if( fTheKernel == "Prod" ){
      SetMGamma(fGammaList);
      fSVKernelFunction = new SVKernelFunction(SVKernelFunction::kSum, MakeKernelList(fMultiKernels,fTheKernel), fmGamma, fGamma, fOrder, fTheta);
   }
   else if( fTheKernel == "Sum" ){
      SetMGamma(fGammaList);
      fSVKernelFunction = new SVKernelFunction(SVKernelFunction::kSum, MakeKernelList(fMultiKernels,fTheKernel), fmGamma, fGamma, fOrder, fTheta);
   }
   else {
      Log() << kWARNING << fTheKernel << " is not a recognised kernel function." << Endl;
      exit(1);
   }
   delete svector;
}

////////////////////////////////////////////////////////////////////////////////
///TODO write IT
/// write training sample (TTree) to file

void TMVA::MethodSVM::WriteWeightsToStream( TFile& ) const
{
}

////////////////////////////////////////////////////////////////////////////////

void  TMVA::MethodSVM::ReadWeightsFromStream( std::istream& istr )
{
   if (fSupportVectors !=0) { delete fSupportVectors; fSupportVectors = 0;}
   fSupportVectors = new std::vector<TMVA::SVEvent*>(0);

   // read configuration from input stream
   istr >> fBparm;

   UInt_t fNsupv;
   // coverity[tainted_data_argument]
   istr >> fNsupv;
   fSupportVectors->reserve(fNsupv);

   Float_t typeTalpha=0.;
   Float_t alpha=0.;
   Int_t typeFlag=-1;
   UInt_t ns = 0;
   std::vector<Float_t>* svector = new std::vector<Float_t>(GetNvar());

   fMaxVars = new TVectorD( GetNvar() );
   fMinVars = new TVectorD( GetNvar() );

   for (UInt_t ievt = 0; ievt < fNsupv; ievt++) {
      istr>>ns;
      istr>>typeTalpha;
      typeFlag = typeTalpha<0?-1:1;
      alpha = typeTalpha<0?-typeTalpha:typeTalpha;
      for (UInt_t ivar = 0; ivar < GetNvar(); ivar++) istr >> svector->at(ivar);

      fSupportVectors->push_back(new SVEvent(svector,alpha,typeFlag,ns));
   }

   for (UInt_t ivar = 0; ivar < GetNvar(); ivar++) istr >> (*fMaxVars)[ivar];

   for (UInt_t ivar = 0; ivar < GetNvar(); ivar++) istr >> (*fMinVars)[ivar];

   delete fSVKernelFunction;
   if (fTheKernel == "Gauss" ) {
      fSVKernelFunction = new SVKernelFunction(1/fDoubleSigmaSquared);
   }
   else {
      SVKernelFunction::EKernelType k = SVKernelFunction::kLinear;
      if(fTheKernel == "Linear")           k = SVKernelFunction::kLinear;
      else if (fTheKernel == "Polynomial") k = SVKernelFunction::kPolynomial;
      else if (fTheKernel == "Sigmoid"   ) k = SVKernelFunction::kSigmoidal;
      else {
         Log() << kFATAL <<"Unknown kernel function found in weight file!" << Endl;
      }
      fSVKernelFunction = new SVKernelFunction();
      fSVKernelFunction->setCompatibilityParams(k, fOrder, fTheta, fKappa);
   }
   delete svector;
}

////////////////////////////////////////////////////////////////////////////////
/// TODO write IT

void TMVA::MethodSVM::ReadWeightsFromStream( TFile& /* fFin */ )
{
}

////////////////////////////////////////////////////////////////////////////////
/// returns MVA value for given event

Double_t TMVA::MethodSVM::GetMvaValue( Double_t* err, Double_t* errUpper )
{
   Double_t myMVA = 0;

   // TODO: avoid creation of a new SVEvent every time (Joerg)
   SVEvent* ev = new SVEvent( GetEvent(), 0. ); // check for specificators

   for (UInt_t ievt = 0; ievt < fSupportVectors->size() ; ievt++) {
      myMVA += ( fSupportVectors->at(ievt)->GetAlpha()
                 * fSupportVectors->at(ievt)->GetTypeFlag()
                 * fSVKernelFunction->Evaluate( fSupportVectors->at(ievt), ev ) );
   }

   delete ev;

   myMVA -= fBparm;

   // cannot determine error
   NoErrorCalc(err, errUpper);

   // 08/12/09: changed sign here to make results agree with convention signal=1
   return 1.0/(1.0 + TMath::Exp(myMVA));
}
////////////////////////////////////////////////////////////////////////////////

const std::vector<Float_t>& TMVA::MethodSVM::GetRegressionValues()
{
   if( fRegressionReturnVal == NULL )
      fRegressionReturnVal = new std::vector<Float_t>();
   fRegressionReturnVal->clear();

   Double_t myMVA = 0;

   const Event *baseev = GetEvent();
   SVEvent* ev = new SVEvent( baseev,0. ); //check for specificators

   for (UInt_t ievt = 0; ievt < fSupportVectors->size() ; ievt++) {
      myMVA += ( fSupportVectors->at(ievt)->GetDeltaAlpha()
                 *fSVKernelFunction->Evaluate( fSupportVectors->at(ievt), ev ) );
   }
   myMVA += fBparm;
   Event * evT = new Event(*baseev);
   evT->SetTarget(0,myMVA);

   const Event* evT2 = GetTransformationHandler().InverseTransform( evT );

   fRegressionReturnVal->push_back(evT2->GetTarget(0));

   delete evT;

   delete ev;

   return *fRegressionReturnVal;
}

////////////////////////////////////////////////////////////////////////////////
/// write specific classifier response

void TMVA::MethodSVM::MakeClassSpecific( std::ostream& fout, const TString& className ) const
{
   const int fNsupv = fSupportVectors->size();
   fout << "   // not implemented for class: \"" << className << "\"" << std::endl;
   fout << "   float        fBparameter;" << std::endl;
   fout << "   int          fNOfSuppVec;" << std::endl;
   fout << "   static float fAllSuppVectors[][" << fNsupv << "];" << std::endl;
   fout << "   static float fAlphaTypeCoef[" << fNsupv << "];" << std::endl;
   fout << std::endl;
   fout << "   // Kernel parameter(s) " << std::endl;
   fout << "   float fGamma;"  << std::endl;
   fout << "};" << std::endl;
   fout << "" << std::endl;

   //Initialize function definition
   fout << "inline void " << className << "::Initialize() " << std::endl;
   fout << "{" << std::endl;
   fout << "   fBparameter = " << fBparm << ";" << std::endl;
   fout << "   fNOfSuppVec = " << fNsupv << ";" << std::endl;
   fout << "   fGamma = " << fGamma << ";" <<std::endl;
   fout << "}" << std::endl;
   fout << std::endl;

   // GetMvaValue__ function definition
   fout << "inline double " << className << "::GetMvaValue__(const std::vector<double>& inputValues ) const" << std::endl;
   fout << "{" << std::endl;
   fout << "   double mvaval = 0; " << std::endl;
   fout << "   double temp = 0; " << std::endl;
   fout << std::endl;
   fout << "   for (int ievt = 0; ievt < fNOfSuppVec; ievt++ ){" << std::endl;
   fout << "      temp = 0;" << std::endl;
   fout << "      for ( unsigned int ivar = 0; ivar < GetNvar(); ivar++ ) {" << std::endl;

   fout << "         temp += (fAllSuppVectors[ivar][ievt] - inputValues[ivar])  " << std::endl;
   fout << "               * (fAllSuppVectors[ivar][ievt] - inputValues[ivar]); " << std::endl;
   fout << "      }" << std::endl;
   fout << "      mvaval += fAlphaTypeCoef[ievt] * exp( -fGamma * temp ); " << std::endl;

   fout << "   }" << std::endl;
   fout << "   mvaval -= fBparameter;" << std::endl;
   fout << "   return 1./(1. + exp(mvaval));" << std::endl;
   fout << "}" << std::endl;
   fout << "// Clean up" << std::endl;
   fout << "inline void " << className << "::Clear() " << std::endl;
   fout << "{" << std::endl;
   fout << "   // nothing to clear " << std::endl;
   fout << "}" << std::endl;
   fout << "" << std::endl;

   // define support vectors
   fout << "float " << className << "::fAlphaTypeCoef[] =" << std::endl;
   fout << "{ ";
   for (Int_t isv = 0; isv < fNsupv; isv++) {
      fout << fSupportVectors->at(isv)->GetDeltaAlpha() * fSupportVectors->at(isv)->GetTypeFlag();
      if (isv < fNsupv-1) fout << ", ";
   }
   fout << " };" << std::endl << std::endl;

   fout << "float " << className << "::fAllSuppVectors[][" << fNsupv << "] =" << std::endl;
   fout << "{";
   for (UInt_t ivar = 0; ivar < GetNvar(); ivar++) {
      fout << std::endl;
      fout << "   { ";
      for (Int_t isv = 0; isv < fNsupv; isv++){
         fout << fSupportVectors->at(isv)->GetDataVector()->at(ivar);
         if (isv < fNsupv-1) fout << ", ";
      }
      fout << " }";
      if (ivar < GetNvar()-1) fout << ", " << std::endl;
      else                    fout << std::endl;
   }
   fout << "};" << std::endl<< std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// get help message text
///
/// typical length of text line:
///         "|--------------------------------------------------------------|"

void TMVA::MethodSVM::GetHelpMessage() const
{
   Log() << Endl;
   Log() << gTools().Color("bold") << "--- Short description:" << gTools().Color("reset") << Endl;
   Log() << Endl;
   Log() << "The Support Vector Machine (SVM) builds a hyperplane separating" << Endl;
   Log() << "signal and background events (vectors) using the minimal subset of " << Endl;
   Log() << "all vectors used for training (support vectors). The extension to" << Endl;
   Log() << "the non-linear case is performed by mapping input vectors into a " << Endl;
   Log() << "higher-dimensional feature space in which linear separation is " << Endl;
   Log() << "possible. The use of the kernel functions thereby eliminates the " << Endl;
   Log() << "explicit transformation to the feature space. The implemented SVM " << Endl;
   Log() << "algorithm performs the classification tasks using linear, polynomial, " << Endl;
   Log() << "Gaussian and sigmoidal kernel functions. The Gaussian kernel allows " << Endl;
   Log() << "to apply any discriminant shape in the input space." << Endl;
   Log() << Endl;
   Log() << gTools().Color("bold") << "--- Performance optimisation:" << gTools().Color("reset") << Endl;
   Log() << Endl;
   Log() << "SVM is a general purpose non-linear classification method, which " << Endl;
   Log() << "does not require data preprocessing like decorrelation or Principal " << Endl;
   Log() << "Component Analysis. It generalises quite well and can handle analyses " << Endl;
   Log() << "with large numbers of input variables." << Endl;
   Log() << Endl;
   Log() << gTools().Color("bold") << "--- Performance tuning via configuration options:" << gTools().Color("reset") << Endl;
   Log() << Endl;
   Log() << "Optimal performance requires primarily a proper choice of the kernel " << Endl;
   Log() << "parameters (the width \"Sigma\" in case of Gaussian kernel) and the" << Endl;
   Log() << "cost parameter \"C\". The user must optimise them empirically by running" << Endl;
   Log() << "SVM several times with different parameter sets. The time needed for " << Endl;
   Log() << "each evaluation scales like the square of the number of training " << Endl;
   Log() << "events so that a coarse preliminary tuning should be performed on " << Endl;
   Log() << "reduced data sets." << Endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Optimize Tuning Parameters
/// This is used to optimise the kernel function parameters and cost. All kernel parameters
/// are optimised by default with default ranges, however the parameters to be optimised can
/// be set when booking the method with the option Tune.
///
/// Example:
///
/// "Tune=Gamma[0.01;1.0;100]" would only tune the RBF Gamma between 0.01 and 1.0
/// with 100 steps.

std::map<TString,Double_t> TMVA::MethodSVM::OptimizeTuningParameters(TString fomType, TString fitType)
{
   // Call the Optimizer with the set of kernel parameters and ranges that are meant to be tuned.
   std::map< TString,std::vector<Double_t> > optVars;
   // Get parameters and options specified in booking of method.
   if(fTune != "All"){
      optVars= GetTuningOptions();
   }
   std::map< TString,std::vector<Double_t> >::iterator iter;
   // Fill all the tuning parameters that should be optimized into a map
   std::map<TString,TMVA::Interval*> tuneParameters;
   std::map<TString,Double_t> tunedParameters;
   // Note: the 3rd parameter in the interval is the "number of bins", NOT the stepsize!!
   // The actual values are always read from the middle of the bins.
   Log() << kINFO << "Using the " << fTheKernel << " kernel." << Endl;
   // Setup map of parameters based on the specified options or defaults.
   if( fTheKernel == "RBF" ){
      if(fTune == "All"){
         tuneParameters.insert(std::pair<TString,Interval*>("Gamma",new Interval(0.01,1.,100)));
         tuneParameters.insert(std::pair<TString,Interval*>("C",new Interval(0.01,1.,100)));
      }
      else{
         for(iter=optVars.begin(); iter!=optVars.end(); ++iter){
            if( iter->first == "Gamma" || iter->first == "C"){
               tuneParameters.insert(std::pair<TString,Interval*>(iter->first, new Interval(iter->second.at(0),iter->second.at(1),iter->second.at(2))));
            }
            else{
               Log() << kWARNING << iter->first << " is not a recognised tuneable parameter." << Endl;
               exit(1);
            }
         }
      }
   }
   else if( fTheKernel == "Polynomial" ){
      if (fTune == "All"){
         tuneParameters.insert(std::pair<TString,Interval*>("Order", new Interval(1,10,10)));
         tuneParameters.insert(std::pair<TString,Interval*>("Theta", new Interval(0.01,1.,100)));
         tuneParameters.insert(std::pair<TString,Interval*>("C", new Interval(0.01,1.,100)));
      }
      else{
         for(iter=optVars.begin(); iter!=optVars.end(); ++iter){
            if( iter->first == "Theta" || iter->first == "C"){
               tuneParameters.insert(std::pair<TString,Interval*>(iter->first, new Interval(iter->second.at(0),iter->second.at(1),iter->second.at(2))));
            }
            else if( iter->first == "Order"){
               tuneParameters.insert(std::pair<TString,Interval*>(iter->first, new Interval(iter->second.at(0),iter->second.at(1),iter->second.at(2))));
            }
            else{
               Log() << kWARNING << iter->first << " is not a recognised tuneable parameter." << Endl;
               exit(1);
            }
         }
      }
   }
   else if( fTheKernel == "MultiGauss" ){
      if (fTune == "All"){
         for(int i=0; i<fNumVars; i++){
            stringstream s;
            s << fVarNames.at(i);
            string str = "Gamma_" + s.str();
            tuneParameters.insert(std::pair<TString,Interval*>(str,new Interval(0.01,1.,100)));
         }
         tuneParameters.insert(std::pair<TString,Interval*>("C",new Interval(0.01,1.,100)));
      } else {
         for(iter=optVars.begin(); iter!=optVars.end(); ++iter){
            if( iter->first == "GammaList"){
               for(int j=0; j<fNumVars; j++){
                  stringstream s;
                  s << fVarNames.at(j);
                  string str = "Gamma_" + s.str();
                  tuneParameters.insert(std::pair<TString,Interval*>(str, new Interval(iter->second.at(0),iter->second.at(1),iter->second.at(2))));
               }
            }
            else if( iter->first == "C"){
               tuneParameters.insert(std::pair<TString,Interval*>(iter->first, new Interval(iter->second.at(0),iter->second.at(1),iter->second.at(2))));
            }
            else{
               Log() << kWARNING << iter->first << " is not a recognised tuneable parameter." << Endl;
               exit(1);
            }
         }
      }
   }
   else if( fTheKernel == "Prod" ){
      std::stringstream tempstring(fMultiKernels);
      std::string value;
      while (std::getline(tempstring,value,'*')){
         if(value == "RBF"){
            tuneParameters.insert(std::pair<TString,Interval*>("Gamma",new Interval(0.01,1.,100)));
         }
         else if(value == "MultiGauss"){
            for(int i=0; i<fNumVars; i++){
               stringstream s;
               s << fVarNames.at(i);
               string str = "Gamma_" + s.str();
               tuneParameters.insert(std::pair<TString,Interval*>(str,new Interval(0.01,1.,100)));
            }
         }
         else if(value == "Polynomial"){
            tuneParameters.insert(std::pair<TString,Interval*>("Order",new Interval(1,10,10)));
            tuneParameters.insert(std::pair<TString,Interval*>("Theta",new Interval(0.0,1.0,101)));
         }
         else {
            Log() << kWARNING << value << " is not a recognised kernel function." << Endl;
            exit(1);
         }
      }
      tuneParameters.insert(std::pair<TString,Interval*>("C",new Interval(0.01,1.,100)));
   }
   else if( fTheKernel == "Sum" ){
      std::stringstream tempstring(fMultiKernels);
      std::string value;
      while (std::getline(tempstring,value,'+')){
         if(value == "RBF"){
            tuneParameters.insert(std::pair<TString,Interval*>("Gamma",new Interval(0.01,1.,100)));
         }
         else if(value == "MultiGauss"){
            for(int i=0; i<fNumVars; i++){
               stringstream s;
               s << fVarNames.at(i);
               string str = "Gamma_" + s.str();
               tuneParameters.insert(std::pair<TString,Interval*>(str,new Interval(0.01,1.,100)));
            }
         }
         else if(value == "Polynomial"){
            tuneParameters.insert(std::pair<TString,Interval*>("Order",new Interval(1,10,10)));
            tuneParameters.insert(std::pair<TString,Interval*>("Theta",new Interval(0.0,1.0,101)));
         }
         else {
            Log() << kWARNING << value << " is not a recognised kernel function." << Endl;
            exit(1);
         }
      }
      tuneParameters.insert(std::pair<TString,Interval*>("C",new Interval(0.01,1.,100)));
   }
   else {
      Log() << kWARNING << fTheKernel << " is not a recognised kernel function." << Endl;
      exit(1);
   }
   Log() << kINFO << " the following SVM parameters will be tuned on the respective *grid*\n" << Endl;
   std::map<TString,TMVA::Interval*>::iterator it;
   for(it=tuneParameters.begin(); it!=tuneParameters.end(); ++it){
      Log() << kWARNING << it->first <<Endl;
      std::ostringstream oss;
      (it->second)->Print(oss);
      Log()<<oss.str();
      Log()<<Endl;
   }
   OptimizeConfigParameters optimize(this, tuneParameters, fomType, fitType);
   tunedParameters=optimize.optimize();

   return tunedParameters;

}

////////////////////////////////////////////////////////////////////////////////
/// Set the tuning parameters according to the argument
void TMVA::MethodSVM::SetTuneParameters(std::map<TString,Double_t> tuneParameters)
{
   std::map<TString,Double_t>::iterator it;
   if( fTheKernel == "RBF" ){
      for(it=tuneParameters.begin(); it!=tuneParameters.end(); ++it){
         Log() << kWARNING << it->first << " = " << it->second << Endl;
         if (it->first == "Gamma"){
            SetGamma (it->second);
         }
         else if(it->first == "C"){
            SetCost (it->second);
         }
         else {
            Log() << kFATAL << " SetParameter for " << it->first << " not implemented " << Endl;
         }
      }
   }
   else if( fTheKernel == "MultiGauss" ){
      fmGamma.clear();
      for(int i=0; i<fNumVars; i++){
         stringstream s;
         s << fVarNames.at(i);
         string str = "Gamma_" + s.str();
         Log() << kWARNING << tuneParameters.find(str)->first << " = " << tuneParameters.find(str)->second << Endl;
         fmGamma.push_back(tuneParameters.find(str)->second);
      }
      for(it=tuneParameters.begin(); it!=tuneParameters.end(); ++it){
         if (it->first == "C"){
            Log() << kWARNING << it->first << " = " << it->second << Endl;
            SetCost(it->second);
            break;
         }
      }
   }
   else if( fTheKernel == "Polynomial" ){
      for(it=tuneParameters.begin(); it!=tuneParameters.end(); ++it){
         Log() << kWARNING << it->first << " = " << it->second << Endl;
         if (it->first == "Order"){
            SetOrder(it->second);
         }
         else if (it->first == "Theta"){
            SetTheta(it->second);
         }
         else if(it->first == "C"){ SetCost (it->second);
         }
         else if(it->first == "Mult"){
            SetMult(it->second);
         }
         else{
            Log() << kFATAL << " SetParameter for " << it->first << " not implemented " << Endl;
         }
      }
   }
   else if( fTheKernel == "Prod" || fTheKernel == "Sum"){
      fmGamma.clear();
      for(it=tuneParameters.begin(); it!=tuneParameters.end(); ++it){
         bool foundParam = false;
         Log() << kWARNING << it->first << " = " << it->second << Endl;
         for(int i=0; i<fNumVars; i++){
            stringstream s;
            s << fVarNames.at(i);
            string str = "Gamma_" + s.str();
            if(it->first == str){
               fmGamma.push_back(it->second);
               foundParam = true;
            }
         }
         if (it->first == "Gamma"){
            SetGamma (it->second);
            foundParam = true;
         }
         else if (it->first == "Order"){
            SetOrder (it->second);
            foundParam = true;
         }
         else if (it->first == "Theta"){
            SetTheta (it->second);
            foundParam = true;
         }
         else if (it->first == "C"){ SetCost (it->second);
            SetCost (it->second);
            foundParam = true;
         }
         else{
            if(!foundParam){
               Log() << kFATAL << " SetParameter for " << it->first << " not implemented " << Endl;
            }
         }
      }
   }
   else {
      Log() << kWARNING << fTheKernel << " is not a recognised kernel function." << Endl;
      exit(1);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Takes as input a string of values for multigaussian gammas and splits it, filling the
/// gamma vector required by the SVKernelFunction. Example: "GammaList=0.1,0.2,0.3" would
/// make a vector with Gammas of 0.1,0.2 & 0.3 corresponding to input variables 1,2 & 3
/// respectively.
void TMVA::MethodSVM::SetMGamma(std::string & mg){
   std::stringstream tempstring(mg);
   Float_t value;
   while (tempstring >> value){
      fmGamma.push_back(value);

      if (tempstring.peek() == ','){
         tempstring.ignore();
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Produces GammaList string for multigaussian kernel to be written to xml file
void TMVA::MethodSVM::GetMGamma(const std::vector<float> & gammas){
   std::ostringstream tempstring;
   for(UInt_t i = 0; i<gammas.size(); ++i){
      tempstring << gammas.at(i);
      if(i!=(gammas.size()-1)){
         tempstring << ",";
      }
   }
   fGammaList= tempstring.str();
}

////////////////////////////////////////////////////////////////////////////////
/// MakeKernelList
/// Function providing string manipulation for product or sum of kernels functions
/// to take list of kernels specified in the booking of the method and provide a vector
/// of SV kernels to iterate over in SVKernelFunction.
///
/// Example:
///
/// "KernelList=RBF*Polynomial" would use a product of the RBF and Polynomial
/// kernels.

std::vector<TMVA::SVKernelFunction::EKernelType> TMVA::MethodSVM::MakeKernelList(std::string multiKernels, TString kernel)
{
   std::vector<TMVA::SVKernelFunction::EKernelType> kernelsList;
   std::stringstream tempstring(multiKernels);
   std::string value;
   if(kernel=="Prod"){
      while (std::getline(tempstring,value,'*')){
         if(value == "RBF"){ kernelsList.push_back(SVKernelFunction::kRBF);}
         else if(value == "MultiGauss"){
            kernelsList.push_back(SVKernelFunction::kMultiGauss);
            if(fGammas!=""){
               SetMGamma(fGammas);
            }
         }
         else if(value == "Polynomial"){ kernelsList.push_back(SVKernelFunction::kPolynomial);}
         else {
            Log() << kWARNING << value << " is not a recognised kernel function." << Endl;
            exit(1);
         }
      }
   }
   else if(kernel=="Sum"){
      while (std::getline(tempstring,value,'+')){
         if(value == "RBF"){ kernelsList.push_back(SVKernelFunction::kRBF);}
         else if(value == "MultiGauss"){
            kernelsList.push_back(SVKernelFunction::kMultiGauss);
            if(fGammas!=""){
               SetMGamma(fGammas);
            }
         }
         else if(value == "Polynomial"){ kernelsList.push_back(SVKernelFunction::kPolynomial);}
         else {
            Log() << kWARNING << value << " is not a recognised kernel function." << Endl;
            exit(1);
         }
      }
   }
   else {
      Log() << kWARNING << "Unable to split MultiKernels. Delimiters */+ required." << Endl;
      exit(1);
   }
   return kernelsList;
}

////////////////////////////////////////////////////////////////////////////////
/// GetTuningOptions
/// Function to allow for ranges and number of steps (for scan) when optimising kernel
/// function parameters. Specified when booking the method after the parameter to be
/// optimised between square brackets with each value separated by ;, the first value
/// is the lower limit, the second the upper limit and the third is the number of steps.
/// Example: "Tune=Gamma[0.01;1.0;100]" would only tune the RBF Gamma between 0.01 and
/// 100 steps.
std::map< TString,std::vector<Double_t> > TMVA::MethodSVM::GetTuningOptions()
{
   std::map< TString,std::vector<Double_t> > optVars;
   std::stringstream tempstring(fTune);
   std::string value;
   while (std::getline(tempstring,value,',')){
      unsigned first = value.find('[')+1;
      unsigned last = value.find_last_of(']');
      std::string optParam = value.substr(0,first-1);
      std::stringstream strNew (value.substr(first,last-first));
      Double_t optInterval;
      std::vector<Double_t> tempVec;
      UInt_t i = 0;
      while (strNew >> optInterval){
         tempVec.push_back(optInterval);
         if (strNew.peek() == ';'){
            strNew.ignore();
         }
         ++i;
      }
      if(i != 3 && i == tempVec.size()){
         if(optParam == "C" || optParam == "Gamma" || optParam == "GammaList" || optParam == "Theta"){
            switch(i){
            case 0:
               tempVec.push_back(0.01);
            case 1:
               tempVec.push_back(1.);
            case 2:
               tempVec.push_back(100);
            }
         }
         else if(optParam == "Order"){
            switch(i){
            case 0:
               tempVec.push_back(1);
            case 1:
               tempVec.push_back(10);
            case 2:
               tempVec.push_back(10);
            }
         }
         else{
            Log() << kWARNING << optParam << " is not a recognised tuneable parameter." << Endl;
            exit(1);
         }
      }
      optVars.insert(std::pair<TString,std::vector<Double_t> >(optParam,tempVec));
   }
   return optVars;
}

////////////////////////////////////////////////////////////////////////////////
/// getLoss
/// Calculates loss for testing dataset. The loss function can be specified when
/// booking the method, otherwise defaults to hinge loss. Currently not used however
/// is accesible if required.

Double_t TMVA::MethodSVM::getLoss(TString lossFunction){
   Double_t loss = 0.0;
   Double_t sumW = 0.0;
   Double_t temp = 0.0;
   Data()->SetCurrentType(Types::kTesting);
   ResultsClassification* mvaRes = dynamic_cast<ResultsClassification*> ( Data()->GetResults(GetMethodName(),Types::kTesting, Types::kClassification) );
   for (Long64_t ievt=0; ievt<GetNEvents(); ievt++) {
      const Event* ev = GetEvent(ievt);
      Float_t v = (*mvaRes)[ievt][0];
      Float_t w = ev->GetWeight();
      if(DataInfo().IsSignal(ev)){
         if(lossFunction == "hinge"){
            temp += w*(1-v);
         }
         else if(lossFunction  == "exp"){
            temp += w*TMath::Exp(-v);
         }
         else if(lossFunction == "binomial"){
            temp += w*TMath::Log(1+TMath::Exp(-2*v));
         }
         else{
            Log() << kWARNING << lossFunction << " is not a recognised loss function." << Endl;
            exit(1);
         }
      }
      else{
         if(lossFunction == "hinge"){
            temp += w*v;
         }
         else if(lossFunction == "exp"){
            temp += w*TMath::Exp(-(1-v));
         }
         else if(lossFunction == "binomial"){
            temp += w*TMath::Log(1+TMath::Exp(-2*(1-v)));
         }
         else{
            Log() << kWARNING << lossFunction << " is not a recognised loss function." << Endl;
            exit(1);
         }
      }
      sumW += w;
   }
   loss = temp/sumW;

   return loss;
}
