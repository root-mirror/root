// @(#)root/tmva $Id$
// Author: Andreas Hoecker, Joerg Stelzer, Helge Voss, Kai Voss,Or Cohen, Jan Therhaag, Eckhard von Toerne

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : MethodCompositeBase                                                   *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Virtual base class for all MVA method                                     *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker    <Andreas.Hocker@cern.ch>   - CERN, Switzerland         *
 *      Peter Speckmayer   <Peter.Speckmazer@cern.ch> - CERN, Switzerland         *
 *      Joerg Stelzer      <Joerg.Stelzer@cern.ch>    - CERN, Switzerland         *
 *      Helge Voss         <Helge.Voss@cern.ch>       - MPI-K Heidelberg, Germany *
 *      Jan Therhaag       <Jan.Therhaag@cern.ch>     - U of Bonn, Germany        *
 *      Eckhard v. Toerne  <evt@uni-bonn.de>          - U of Bonn, Germany        *
 *                                                                                *
 * Copyright (c) 2005-2011:                                                       *
 *      CERN, Switzerland                                                         *
 *      U. of Victoria, Canada                                                    *
 *      MPI-K Heidelberg, Germany                                                 *
 *      U. of Bonn, Germany                                                       *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

/*! \class TMVA::MethodBoost
\ingroup TMVA

Class for boosting a TMVA method

This class is meant to boost a single classifier. Boosting means
training the classifier a few times. Every time the weights of the
events are modified according to how well the classifier performed
on the test sample.

*/

#include "TMVA/MethodBoost.h"

#include "TMVA/ClassifierFactory.h"
#include "TMVA/Config.h"
#include "TMVA/Configurable.h"
#include "TMVA/DataSet.h"
#include "TMVA/DataSetInfo.h"
#include "TMVA/IMethod.h"
#include "TMVA/MethodBase.h"
#include "TMVA/MethodCategory.h"
#include "TMVA/MethodCompositeBase.h"
#include "TMVA/MethodDT.h"
#include "TMVA/MethodFisher.h"
#include "TMVA/PDF.h"
#include "TMVA/Results.h"
#include "TMVA/Timer.h"
#include "TMVA/Tools.h"
#include "TMVA/Types.h"

#include "TMVA/SeparationBase.h"
#include "TMVA/MisClassificationError.h"
#include "TMVA/GiniIndex.h"
#include "TMVA/CrossEntropy.h"
#include "TMVA/RegressionVariance.h"
#include "TMVA/QuickMVAProbEstimator.h"

#include "TRandom3.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TSpline.h"
#include "TDirectory.h"
#include "TTree.h"

#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>


REGISTER_METHOD(Boost)

ClassImp(TMVA::MethodBoost);

////////////////////////////////////////////////////////////////////////////////

   TMVA::MethodBoost::MethodBoost( const TString& jobName,
                                   const TString& methodTitle,
                                   DataSetInfo& theData,
                                   const TString& theOption ) :
   TMVA::MethodCompositeBase( jobName, Types::kBoost, methodTitle, theData, theOption)
   , fBoostNum(0)
   , fDetailedMonitoring(kFALSE)
   , fAdaBoostBeta(0)
   , fRandomSeed(0)
   , fBaggedSampleFraction(0)
   , fBoostedMethodTitle(methodTitle)
   , fBoostedMethodOptions(theOption)
   , fMonitorBoostedMethod(kFALSE)
   , fMonitorTree(0)
   , fBoostWeight(0)
   , fMethodError(0)
   , fROC_training(0.0)
   , fOverlap_integral(0.0)
   , fMVAvalues(0)
{
   fMVAvalues = new std::vector<Float_t>;
   fDataSetManager = NULL;
   fHistoricBoolOption = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////

TMVA::MethodBoost::MethodBoost( DataSetInfo& dsi,
                                const TString& theWeightFile)
   : TMVA::MethodCompositeBase( Types::kBoost, dsi, theWeightFile)
   , fBoostNum(0)
   , fDetailedMonitoring(kFALSE)
   , fAdaBoostBeta(0)
   , fRandomSeed(0)
   , fBaggedSampleFraction(0)
   , fBoostedMethodTitle("")
   , fBoostedMethodOptions("")
   , fMonitorBoostedMethod(kFALSE)
   , fMonitorTree(0)
   , fBoostWeight(0)
   , fMethodError(0)
   , fROC_training(0.0)
   , fOverlap_integral(0.0)
   , fMVAvalues(0)
{
   fMVAvalues = new std::vector<Float_t>;
   fDataSetManager = NULL;
   fHistoricBoolOption = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TMVA::MethodBoost::~MethodBoost( void )
{
   fMethodWeight.clear();

   // the histogram themselves are deleted when the file is closed

   fTrainSigMVAHist.clear();
   fTrainBgdMVAHist.clear();
   fBTrainSigMVAHist.clear();
   fBTrainBgdMVAHist.clear();
   fTestSigMVAHist.clear();
   fTestBgdMVAHist.clear();

   if (fMVAvalues) {
      delete fMVAvalues;
      fMVAvalues = 0;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Boost can handle classification with 2 classes and regression with one regression-target

Bool_t TMVA::MethodBoost::HasAnalysisType( Types::EAnalysisType type, UInt_t numberClasses, UInt_t /*numberTargets*/ )
{
   if (type == Types::kClassification && numberClasses == 2) return kTRUE;
   //   if (type == Types::kRegression && numberTargets == 1) return kTRUE;
   return kFALSE;
}


////////////////////////////////////////////////////////////////////////////////

void TMVA::MethodBoost::DeclareOptions()
{
   DeclareOptionRef( fBoostNum = 1, "Boost_Num",
                     "Number of times the classifier is boosted" );

   DeclareOptionRef( fMonitorBoostedMethod = kTRUE, "Boost_MonitorMethod",
                     "Write monitoring histograms for each boosted classifier" );

   DeclareOptionRef( fDetailedMonitoring = kFALSE, "Boost_DetailedMonitoring",
                     "Produce histograms for detailed boost  monitoring" );

   DeclareOptionRef( fBoostType  = "AdaBoost", "Boost_Type", "Boosting type for the classifiers" );
   AddPreDefVal(TString("RealAdaBoost"));
   AddPreDefVal(TString("AdaBoost"));
   AddPreDefVal(TString("Bagging"));

   DeclareOptionRef(fBaggedSampleFraction=.6,"Boost_BaggedSampleFraction","Relative size of bagged event sample to original size of the data sample (used whenever bagging is used)" );

   DeclareOptionRef( fAdaBoostBeta = 1.0, "Boost_AdaBoostBeta",
                     "The ADA boost parameter that sets the effect of every boost step on the events' weights" );

   DeclareOptionRef( fTransformString = "step", "Boost_Transform",
                     "Type of transform applied to every boosted method linear, log, step" );
   AddPreDefVal(TString("step"));
   AddPreDefVal(TString("linear"));
   AddPreDefVal(TString("log"));
   AddPreDefVal(TString("gauss"));

   DeclareOptionRef( fRandomSeed = 0, "Boost_RandomSeed",
                     "Seed for random number generator used for bagging" );

   TMVA::MethodCompositeBase::fMethods.reserve(fBoostNum);
}

////////////////////////////////////////////////////////////////////////////////
/// options that are used ONLY for the READER to ensure backward compatibility
///   they are hence without any effect (the reader is only reading the training
///   options that HAD been used at the training of the .xml weight file at hand

void TMVA::MethodBoost::DeclareCompatibilityOptions()
{

   MethodBase::DeclareCompatibilityOptions();

   DeclareOptionRef( fHistoricOption = "ByError", "Boost_MethodWeightType",
                     "How to set the final weight of the boosted classifiers" );
   AddPreDefVal(TString("ByError"));
   AddPreDefVal(TString("Average"));
   AddPreDefVal(TString("ByROC"));
   AddPreDefVal(TString("ByOverlap"));
   AddPreDefVal(TString("LastMethod"));

   DeclareOptionRef( fHistoricOption = "step", "Boost_Transform",
                     "Type of transform applied to every boosted method linear, log, step" );
   AddPreDefVal(TString("step"));
   AddPreDefVal(TString("linear"));
   AddPreDefVal(TString("log"));
   AddPreDefVal(TString("gauss"));

   // this option here
   //DeclareOptionRef( fBoostType  = "AdaBoost", "Boost_Type", "Boosting type for the classifiers" );
   // still exists, but these two possible values
   AddPreDefVal(TString("HighEdgeGauss"));
   AddPreDefVal(TString("HighEdgeCoPara"));
   // have been deleted .. hope that works :)

   DeclareOptionRef( fHistoricBoolOption, "Boost_RecalculateMVACut",
                     "Recalculate the classifier MVA Signallike cut at every boost iteration" );

}

////////////////////////////////////////////////////////////////////////////////
/// just registering the string from which the boosted classifier will be created

Bool_t TMVA::MethodBoost::BookMethod( Types::EMVA theMethod, TString methodTitle, TString theOption )
{
   fBoostedMethodName     = Types::Instance().GetMethodName( theMethod );
   fBoostedMethodTitle    = methodTitle;
   fBoostedMethodOptions  = theOption;
   TString opts=theOption;
   opts.ToLower();
   //    if (opts.Contains("vartransform")) Log() << kFATAL << "It is not possible to use boost in conjunction with variable transform. Please remove either Boost_Num or VarTransform from the option string"<< methodTitle<<Endl;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::MethodBoost::Init()
{
}

////////////////////////////////////////////////////////////////////////////////
/// initialisation routine

void TMVA::MethodBoost::InitHistos()
{

   Results* results = Data()->GetResults(GetMethodName(), Types::kTraining, GetAnalysisType());

   results->Store(new TH1F("MethodWeight","Normalized Classifier Weight",fBoostNum,0,fBoostNum),"ClassifierWeight");
   results->Store(new TH1F("BoostWeight","Boost Weight",fBoostNum,0,fBoostNum),"BoostWeight");
   results->Store(new TH1F("ErrFraction","Error Fraction (by boosted event weights)",fBoostNum,0,fBoostNum),"ErrorFraction");
   if (fDetailedMonitoring){
      results->Store(new TH1F("ROCIntegral_test","ROC integral of single classifier (testing sample)",fBoostNum,0,fBoostNum),"ROCIntegral_test");
      results->Store(new TH1F("ROCIntegralBoosted_test","ROC integral of boosted method (testing sample)",fBoostNum,0,fBoostNum),"ROCIntegralBoosted_test");
      results->Store(new TH1F("ROCIntegral_train","ROC integral of single classifier (training sample)",fBoostNum,0,fBoostNum),"ROCIntegral_train");
      results->Store(new TH1F("ROCIntegralBoosted_train","ROC integral of boosted method (training sample)",fBoostNum,0,fBoostNum),"ROCIntegralBoosted_train");
      results->Store(new TH1F("OverlapIntegal_train","Overlap integral (training sample)",fBoostNum,0,fBoostNum),"Overlap");
   }


   results->GetHist("ClassifierWeight")->GetXaxis()->SetTitle("Index of boosted classifier");
   results->GetHist("ClassifierWeight")->GetYaxis()->SetTitle("Classifier Weight");
   results->GetHist("BoostWeight")->GetXaxis()->SetTitle("Index of boosted classifier");
   results->GetHist("BoostWeight")->GetYaxis()->SetTitle("Boost Weight");
   results->GetHist("ErrorFraction")->GetXaxis()->SetTitle("Index of boosted classifier");
   results->GetHist("ErrorFraction")->GetYaxis()->SetTitle("Error Fraction");
   if (fDetailedMonitoring){
      results->GetHist("ROCIntegral_test")->GetXaxis()->SetTitle("Index of boosted classifier");
      results->GetHist("ROCIntegral_test")->GetYaxis()->SetTitle("ROC integral of single classifier");
      results->GetHist("ROCIntegralBoosted_test")->GetXaxis()->SetTitle("Number of boosts");
      results->GetHist("ROCIntegralBoosted_test")->GetYaxis()->SetTitle("ROC integral boosted");
      results->GetHist("ROCIntegral_train")->GetXaxis()->SetTitle("Index of boosted classifier");
      results->GetHist("ROCIntegral_train")->GetYaxis()->SetTitle("ROC integral of single classifier");
      results->GetHist("ROCIntegralBoosted_train")->GetXaxis()->SetTitle("Number of boosts");
      results->GetHist("ROCIntegralBoosted_train")->GetYaxis()->SetTitle("ROC integral boosted");
      results->GetHist("Overlap")->GetXaxis()->SetTitle("Index of boosted classifier");
      results->GetHist("Overlap")->GetYaxis()->SetTitle("Overlap integral");
   }

   results->Store(new TH1F("SoverBtotal","S/B in reweighted training sample",fBoostNum,0,fBoostNum),"SoverBtotal");
   results->GetHist("SoverBtotal")->GetYaxis()->SetTitle("S/B (boosted sample)");
   results->GetHist("SoverBtotal")->GetXaxis()->SetTitle("Index of boosted classifier");

   results->Store(new TH1F("SeparationGain","SeparationGain",fBoostNum,0,fBoostNum),"SeparationGain");
   results->GetHist("SeparationGain")->GetYaxis()->SetTitle("SeparationGain");
   results->GetHist("SeparationGain")->GetXaxis()->SetTitle("Index of boosted classifier");



   fMonitorTree= new TTree("MonitorBoost","Boost variables");
   fMonitorTree->Branch("iMethod",&fCurrentMethodIdx,"iMethod/I");
   fMonitorTree->Branch("boostWeight",&fBoostWeight,"boostWeight/D");
   fMonitorTree->Branch("errorFraction",&fMethodError,"errorFraction/D");
   fMonitorBoostedMethod = kTRUE;

}


////////////////////////////////////////////////////////////////////////////////

void TMVA::MethodBoost::CheckSetup()
{
   Log() << kDEBUG << "CheckSetup: fBoostType="<<fBoostType << Endl;
   Log() << kDEBUG << "CheckSetup: fAdaBoostBeta="<<fAdaBoostBeta<<Endl;
   Log() << kDEBUG << "CheckSetup: fBoostWeight="<<fBoostWeight<<Endl;
   Log() << kDEBUG << "CheckSetup: fMethodError="<<fMethodError<<Endl;
   Log() << kDEBUG << "CheckSetup: fBoostNum="<<fBoostNum << Endl;
   Log() << kDEBUG << "CheckSetup: fRandomSeed=" << fRandomSeed<< Endl;
   Log() << kDEBUG << "CheckSetup: fTrainSigMVAHist.size()="<<fTrainSigMVAHist.size()<<Endl;
   Log() << kDEBUG << "CheckSetup: fTestSigMVAHist.size()="<<fTestSigMVAHist.size()<<Endl;
   Log() << kDEBUG << "CheckSetup: fMonitorBoostedMethod=" << (fMonitorBoostedMethod? "true" : "false") << Endl;
   Log() << kDEBUG << "CheckSetup: MName=" << fBoostedMethodName << " Title="<< fBoostedMethodTitle<< Endl;
   Log() << kDEBUG << "CheckSetup: MOptions="<< fBoostedMethodOptions << Endl;
   Log() << kDEBUG << "CheckSetup: fMonitorTree=" << fMonitorTree <<Endl;
   Log() << kDEBUG << "CheckSetup: fCurrentMethodIdx=" <<fCurrentMethodIdx << Endl;
   if (fMethods.size()>0) Log() << kDEBUG << "CheckSetup: fMethods[0]" <<fMethods[0]<<Endl;
   Log() << kDEBUG << "CheckSetup: fMethodWeight.size()" << fMethodWeight.size() << Endl;
   if (fMethodWeight.size()>0) Log() << kDEBUG << "CheckSetup: fMethodWeight[0]="<<fMethodWeight[0]<<Endl;
   Log() << kDEBUG << "CheckSetup: trying to repair things" << Endl;

}
////////////////////////////////////////////////////////////////////////////////

void TMVA::MethodBoost::Train()
{
   TDirectory* methodDir( 0 );
   TString     dirName,dirTitle;
   Int_t       StopCounter=0;
   Results* results = Data()->GetResults(GetMethodName(), Types::kTraining, GetAnalysisType());


   InitHistos();

   if (Data()->GetNTrainingEvents()==0) Log() << kFATAL << "<Train> Data() has zero events" << Endl;
   Data()->SetCurrentType(Types::kTraining);

   if (fMethods.size() > 0) fMethods.clear();
   fMVAvalues->resize(Data()->GetNTrainingEvents(), 0.0);

   Log() << kINFO << "Training "<< fBoostNum << " " << fBoostedMethodName << " with title " << fBoostedMethodTitle << " Classifiers ... patience please" << Endl;
   Timer timer( fBoostNum, GetName() );

   ResetBoostWeights();

   // clean boosted method options
   CleanBoostOptions();


   // remove transformations for individual boosting steps
   // the transformation of the main method will be rerouted to each of the boost steps
   Ssiz_t varTrafoStart=fBoostedMethodOptions.Index("~VarTransform=");
   if (varTrafoStart >0) {
      Ssiz_t varTrafoEnd  =fBoostedMethodOptions.Index(":",varTrafoStart);
      if (varTrafoEnd<varTrafoStart)
         varTrafoEnd=fBoostedMethodOptions.Length();
      fBoostedMethodOptions.Remove(varTrafoStart,varTrafoEnd-varTrafoStart);
   }

   //
   // training and boosting the classifiers
   for (fCurrentMethodIdx=0;fCurrentMethodIdx<fBoostNum;fCurrentMethodIdx++) {
      // the first classifier shows the option string output, the rest not
      if (fCurrentMethodIdx>0) TMVA::MsgLogger::InhibitOutput();

      IMethod *method = ClassifierFactory::Instance().Create(
         fBoostedMethodName.Data(), GetJobName(), Form("%s_B%04i", fBoostedMethodTitle.Data(), fCurrentMethodIdx),
         DataInfo(), fBoostedMethodOptions);
      TMVA::MsgLogger::EnableOutput();

      // suppressing the rest of the classifier output the right way
      fCurrentMethod  = (dynamic_cast<MethodBase*>(method));

      if (fCurrentMethod==0) {
         Log() << kFATAL << "uups.. guess the booking of the " << fCurrentMethodIdx << "-th classifier somehow failed" << Endl;
         return; // hope that makes coverity happy (as if fears I might use the pointer later on, not knowing that FATAL exits
      }

      // set fDataSetManager if MethodCategory (to enable Category to create datasetinfo objects) // DSMTEST
      if (fCurrentMethod->GetMethodType() == Types::kCategory) { // DSMTEST
         MethodCategory *methCat = (dynamic_cast<MethodCategory*>(fCurrentMethod)); // DSMTEST
         if (!methCat) // DSMTEST
            Log() << kFATAL << "Method with type kCategory cannot be casted to MethodCategory. /MethodBoost" << Endl; // DSMTEST
         methCat->fDataSetManager = fDataSetManager; // DSMTEST
      } // DSMTEST

      fCurrentMethod->SetMsgType(kWARNING);
      fCurrentMethod->SetupMethod();
      fCurrentMethod->ParseOptions();
      // put SetAnalysisType here for the needs of MLP
      fCurrentMethod->SetAnalysisType( GetAnalysisType() );
      fCurrentMethod->ProcessSetup();
      fCurrentMethod->CheckSetup();


      // reroute transformationhandler
      fCurrentMethod->RerouteTransformationHandler (&(this->GetTransformationHandler()));


      // creating the directory of the classifier
      if(!IsSilentFile())
      {
        if (fMonitorBoostedMethod) {
            methodDir=GetFile()->GetDirectory(dirName=Form("%s_B%04i",fBoostedMethodName.Data(),fCurrentMethodIdx));
            if (methodDir==0) {
                methodDir=BaseDir()->mkdir(dirName,dirTitle=Form("Directory Boosted %s #%04i", fBoostedMethodName.Data(),fCurrentMethodIdx));
            }
            fCurrentMethod->SetMethodDir(methodDir);
            fCurrentMethod->BaseDir()->cd();
        }
      }

      // training
      TMVA::MethodCompositeBase::fMethods.push_back(method);
      timer.DrawProgressBar( fCurrentMethodIdx );
      if (fCurrentMethodIdx==0) MonitorBoost(Types::kBoostProcBegin,fCurrentMethodIdx);
      MonitorBoost(Types::kBeforeTraining,fCurrentMethodIdx);
      TMVA::MsgLogger::InhibitOutput(); //suppressing Logger outside the method
      if (fBoostType=="Bagging") Bagging();  // you want also to train the first classifier on a bagged sample
      SingleTrain();
      TMVA::MsgLogger::EnableOutput();
      if(!IsSilentFile())fCurrentMethod->WriteMonitoringHistosToFile();

      // calculate MVA values of current method for all events in training sample
      // (used later on to get 'misclassified events' etc for the boosting
      CalcMVAValues();

      if(!IsSilentFile()) if (fCurrentMethodIdx==0 && fMonitorBoostedMethod) CreateMVAHistorgrams();

      // get ROC integral and overlap integral for single method on
      // training sample if fMethodWeightType == "ByROC" or the user
      // wants detailed monitoring

      // boosting (reweight training sample)
      MonitorBoost(Types::kBeforeBoosting,fCurrentMethodIdx);
      SingleBoost(fCurrentMethod);

      MonitorBoost(Types::kAfterBoosting,fCurrentMethodIdx);
      results->GetHist("BoostWeight")->SetBinContent(fCurrentMethodIdx+1,fBoostWeight);
      results->GetHist("ErrorFraction")->SetBinContent(fCurrentMethodIdx+1,fMethodError);

      if (fDetailedMonitoring) {
         fROC_training = GetBoostROCIntegral(kTRUE, Types::kTraining, kTRUE);
         results->GetHist("ROCIntegral_test")->SetBinContent(fCurrentMethodIdx+1, GetBoostROCIntegral(kTRUE,  Types::kTesting));
         results->GetHist("ROCIntegralBoosted_test")->SetBinContent(fCurrentMethodIdx+1, GetBoostROCIntegral(kFALSE, Types::kTesting));
         results->GetHist("ROCIntegral_train")->SetBinContent(fCurrentMethodIdx+1, fROC_training);
         results->GetHist("ROCIntegralBoosted_train")->SetBinContent(fCurrentMethodIdx+1, GetBoostROCIntegral(kFALSE, Types::kTraining));
         results->GetHist("Overlap")->SetBinContent(fCurrentMethodIdx+1, fOverlap_integral);
      }



      fMonitorTree->Fill();

      // stop boosting if needed when error has reached 0.5
      // thought of counting a few steps, but it doesn't seem to be necessary
      Log() << kDEBUG << "AdaBoost (methodErr) err = " << fMethodError << Endl;
      if (fMethodError > 0.49999) StopCounter++;
      if (StopCounter > 0 && fBoostType != "Bagging") {
         timer.DrawProgressBar( fBoostNum );
         fBoostNum = fCurrentMethodIdx+1;
         Log() << kINFO << "Error rate has reached 0.5 ("<< fMethodError<<"), boosting process stopped at #" << fBoostNum << " classifier" << Endl;
         if (fBoostNum < 5)
            Log() << kINFO << "The classifier might be too strong to boost with Beta = " << fAdaBoostBeta << ", try reducing it." <<Endl;
         break;
      }
   }

   //as MethodBoost acts not on a private event sample (like MethodBDT does), we need to remember not
   // to leave "boosted" events to the next classifier in the factory

   ResetBoostWeights();

   Timer* timer1= new Timer( fBoostNum, GetName() );
   // normalizing the weights of the classifiers
   for (fCurrentMethodIdx=0;fCurrentMethodIdx<fBoostNum;fCurrentMethodIdx++) {
      // performing post-boosting actions

      timer1->DrawProgressBar( fCurrentMethodIdx );

      if (fCurrentMethodIdx==fBoostNum) {
         Log() << kINFO << "Elapsed time: " << timer1->GetElapsedTime()
               << "                              " << Endl;
      }

      TH1F* tmp = dynamic_cast<TH1F*>( results->GetHist("ClassifierWeight") );
      if (tmp) tmp->SetBinContent(fCurrentMethodIdx+1,fMethodWeight[fCurrentMethodIdx]);

   }

   // Ensure that in case of only 1 boost the method weight equals
   // 1.0.  This avoids unexpected behaviour in case of very bad
   // classifiers which have fBoostWeight=1 or fMethodError=0.5,
   // because their weight would be set to zero.  This behaviour is
   // not ok if one boosts just one time.
   if (fMethods.size()==1)  fMethodWeight[0] = 1.0;

   MonitorBoost(Types::kBoostProcEnd);

   delete timer1;
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::MethodBoost::CleanBoostOptions()
{
   fBoostedMethodOptions=GetOptions();
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::MethodBoost::CreateMVAHistorgrams()
{
   if (fBoostNum <=0) Log() << kFATAL << "CreateHistograms called before fBoostNum is initialized" << Endl;
   // calculating histograms boundaries and creating histograms..
   // nrms = number of rms around the average to use for outline (of the 0 classifier)
   Double_t meanS, meanB, rmsS, rmsB, xmin, xmax, nrms = 10;
   Int_t signalClass = 0;
   if (DataInfo().GetClassInfo("Signal") != 0) {
      signalClass = DataInfo().GetClassInfo("Signal")->GetNumber();
   }
   gTools().ComputeStat( GetEventCollection( Types::kMaxTreeType ), fMVAvalues,
                         meanS, meanB, rmsS, rmsB, xmin, xmax, signalClass );

   fNbins = gConfig().fVariablePlotting.fNbinsXOfROCCurve;
   xmin = TMath::Max( TMath::Min(meanS - nrms*rmsS, meanB - nrms*rmsB ), xmin );
   xmax = TMath::Min( TMath::Max(meanS + nrms*rmsS, meanB + nrms*rmsB ), xmax ) + 0.00001;

   // creating all the histograms
   for (UInt_t imtd=0; imtd<fBoostNum; imtd++) {
      fTrainSigMVAHist .push_back( new TH1F( Form("MVA_Train_S_%04i",imtd), "MVA_Train_S",        fNbins, xmin, xmax ) );
      fTrainBgdMVAHist .push_back( new TH1F( Form("MVA_Train_B%04i", imtd), "MVA_Train_B",        fNbins, xmin, xmax ) );
      fBTrainSigMVAHist.push_back( new TH1F( Form("MVA_BTrain_S%04i",imtd), "MVA_BoostedTrain_S", fNbins, xmin, xmax ) );
      fBTrainBgdMVAHist.push_back( new TH1F( Form("MVA_BTrain_B%04i",imtd), "MVA_BoostedTrain_B", fNbins, xmin, xmax ) );
      fTestSigMVAHist  .push_back( new TH1F( Form("MVA_Test_S%04i",  imtd), "MVA_Test_S",         fNbins, xmin, xmax ) );
      fTestBgdMVAHist  .push_back( new TH1F( Form("MVA_Test_B%04i",  imtd), "MVA_Test_B",         fNbins, xmin, xmax ) );
   }
}

////////////////////////////////////////////////////////////////////////////////
/// resetting back the boosted weights of the events to 1

void TMVA::MethodBoost::ResetBoostWeights()
{
   for (Long64_t ievt=0; ievt<GetNEvents(); ievt++) {
      const Event *ev = Data()->GetEvent(ievt);
      ev->SetBoostWeight( 1.0 );
   }
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::MethodBoost::WriteMonitoringHistosToFile( void ) const
{
   TDirectory* dir=0;
   if (fMonitorBoostedMethod) {
      for (UInt_t imtd=0;imtd<fBoostNum;imtd++) {

         //writing the histograms in the specific classifier's directory
         MethodBase* m = dynamic_cast<MethodBase*>(fMethods[imtd]);
         if (!m) continue;
         dir = m->BaseDir();
         dir->cd();
         fTrainSigMVAHist[imtd]->SetDirectory(dir);
         fTrainSigMVAHist[imtd]->Write();
         fTrainBgdMVAHist[imtd]->SetDirectory(dir);
         fTrainBgdMVAHist[imtd]->Write();
         fBTrainSigMVAHist[imtd]->SetDirectory(dir);
         fBTrainSigMVAHist[imtd]->Write();
         fBTrainBgdMVAHist[imtd]->SetDirectory(dir);
         fBTrainBgdMVAHist[imtd]->Write();
      }
   }

   // going back to the original folder
   BaseDir()->cd();

   fMonitorTree->Write();
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::MethodBoost::TestClassification()
{
   MethodBase::TestClassification();
   if (fMonitorBoostedMethod) {
      UInt_t nloop = fTestSigMVAHist.size();
      if (fMethods.size()<nloop) nloop = fMethods.size();
      //running over all the events and populating the test MVA histograms
      Data()->SetCurrentType(Types::kTesting);
      for (Long64_t ievt=0; ievt<GetNEvents(); ievt++) {
         const Event* ev = GetEvent(ievt);
         Float_t w = ev->GetWeight();
         if (DataInfo().IsSignal(ev)) {
            for (UInt_t imtd=0; imtd<nloop; imtd++) {
               fTestSigMVAHist[imtd]->Fill(fMethods[imtd]->GetMvaValue(),w);
            }
         }
         else {
            for (UInt_t imtd=0; imtd<nloop; imtd++) {
               fTestBgdMVAHist[imtd]->Fill(fMethods[imtd]->GetMvaValue(),w);
            }
         }
      }
      Data()->SetCurrentType(Types::kTraining);
   }
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::MethodBoost::WriteEvaluationHistosToFile(Types::ETreeType treetype)
{
   MethodBase::WriteEvaluationHistosToFile(treetype);
   if (treetype==Types::kTraining) return;
   UInt_t nloop = fTestSigMVAHist.size();
   if (fMethods.size()<nloop) nloop = fMethods.size();
   if (fMonitorBoostedMethod) {
      TDirectory* dir=0;
      for (UInt_t imtd=0;imtd<nloop;imtd++) {
         //writing the histograms in the specific classifier's directory
         MethodBase* mva = dynamic_cast<MethodBase*>(fMethods[imtd]);
         if (!mva) continue;
         dir = mva->BaseDir();
         if (dir==0) continue;
         dir->cd();
         fTestSigMVAHist[imtd]->SetDirectory(dir);
         fTestSigMVAHist[imtd]->Write();
         fTestBgdMVAHist[imtd]->SetDirectory(dir);
         fTestBgdMVAHist[imtd]->Write();
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// process user options

void TMVA::MethodBoost::ProcessOptions()
{
}

////////////////////////////////////////////////////////////////////////////////
/// initialization

void TMVA::MethodBoost::SingleTrain()
{
   Data()->SetCurrentType(Types::kTraining);
   MethodBase* meth = dynamic_cast<MethodBase*>(GetLastMethod());
   if (meth){
       meth->SetSilentFile(IsSilentFile());
       if(IsModelPersistence()){
           TString _fFileDir= DataInfo().GetName();
           _fFileDir+="/"+gConfig().GetIONames().fWeightFileDir;
           meth->SetWeightFileDir(_fFileDir);
       }
       meth->SetModelPersistence(IsModelPersistence());
       meth->TrainMethod();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// find the CUT on the individual MVA that defines an event as
/// correct or misclassified (to be used in the boosting process)

void TMVA::MethodBoost::FindMVACut(MethodBase *method)
{
   if (!method || method->GetMethodType() == Types::kDT ){ return;}

   // creating a fine histograms containing the error rate
   const Int_t nBins=10001;
   Double_t minMVA=150000;
   Double_t maxMVA=-150000;
   for (Long64_t ievt=0; ievt<Data()->GetNEvents(); ievt++) {
      GetEvent(ievt);
      Double_t val=method->GetMvaValue();
      //Helge .. I think one could very well use fMVAValues for that ... -->to do
      if (val>maxMVA) maxMVA=val;
      if (val<minMVA) minMVA=val;
   }
   maxMVA = maxMVA+(maxMVA-minMVA)/nBins;

   Double_t sum = 0.;

   TH1D *mvaS  = new TH1D(Form("MVAS_%d",fCurrentMethodIdx) ,"",nBins,minMVA,maxMVA);
   TH1D *mvaB  = new TH1D(Form("MVAB_%d",fCurrentMethodIdx) ,"",nBins,minMVA,maxMVA);
   TH1D *mvaSC = new TH1D(Form("MVASC_%d",fCurrentMethodIdx),"",nBins,minMVA,maxMVA);
   TH1D *mvaBC = new TH1D(Form("MVABC_%d",fCurrentMethodIdx),"",nBins,minMVA,maxMVA);


   Results* results = Data()->GetResults(GetMethodName(), Types::kTraining, GetAnalysisType());
   if (fDetailedMonitoring){
      results->Store(mvaS, Form("MVAS_%d",fCurrentMethodIdx));
      results->Store(mvaB, Form("MVAB_%d",fCurrentMethodIdx));
      results->Store(mvaSC,Form("MVASC_%d",fCurrentMethodIdx));
      results->Store(mvaBC,Form("MVABC_%d",fCurrentMethodIdx));
   }

   for (Long64_t ievt=0; ievt<Data()->GetNEvents(); ievt++) {

      Double_t weight = GetEvent(ievt)->GetWeight();
      Double_t mvaVal=method->GetMvaValue();
      sum +=weight;
      if (DataInfo().IsSignal(GetEvent(ievt))){
         mvaS->Fill(mvaVal,weight);
      }else {
         mvaB->Fill(mvaVal,weight);
      }
   }
   SeparationBase *sepGain;


   // Boosting should use Misclassification not Gini Index (changed, Helge 31.5.2013)
   // WARNING! It works with Misclassification only if you fix the signal to
   // background at every step. Strangely enough, there are better results
   // ( as seen with BDT ) if you use Gini Index, and accept that sometimes no
   // sensible cut is found - i.e. the cut is then outside the MVA value range,
   // all events are classified as background and then according to the Boost
   // algorithm something is renormed 'automatically' ... so that in the next
   // step again the result is something sensible.
   // Strange ... that THIS is supposed to be right?

   //   SeparationBase *sepGain2 = new MisClassificationError();
   //sepGain = new MisClassificationError();
   sepGain = new GiniIndex();
   //sepGain = new CrossEntropy();

   Double_t sTot = mvaS->GetSum();
   Double_t bTot = mvaB->GetSum();

   mvaSC->SetBinContent(1,mvaS->GetBinContent(1));
   mvaBC->SetBinContent(1,mvaB->GetBinContent(1));
   Double_t sSel=0;
   Double_t bSel=0;
   Double_t separationGain=sepGain->GetSeparationGain(sSel,bSel,sTot,bTot);
   Double_t mvaCut=mvaSC->GetBinLowEdge(1);
   Double_t sSelCut=sSel;
   Double_t bSelCut=bSel;
   //      std::cout << "minMVA =" << minMVA << " maxMVA = " << maxMVA << " width = " << mvaSC->GetBinWidth(1) <<  std::endl;

   //      for (Int_t ibin=1;ibin<=nBins;ibin++) std::cout << " cutvalues[" << ibin<<"]="<<mvaSC->GetBinLowEdge(ibin) << "  " << mvaSC->GetBinCenter(ibin) << std::endl;
   Double_t mvaCutOrientation=1; // 1 if mva > mvaCut --> Signal and -1 if mva < mvaCut (i.e. mva*-1 > mvaCut*-1) --> Signal
   for (Int_t ibin=1;ibin<=nBins;ibin++){
      mvaSC->SetBinContent(ibin,mvaS->GetBinContent(ibin)+mvaSC->GetBinContent(ibin-1));
      mvaBC->SetBinContent(ibin,mvaB->GetBinContent(ibin)+mvaBC->GetBinContent(ibin-1));

      sSel=mvaSC->GetBinContent(ibin);
      bSel=mvaBC->GetBinContent(ibin);

      // if (ibin==nBins){
      //    std::cout << "Last bin s="<< sSel <<" b="<<bSel << " s="<< sTot-sSel <<" b="<<bTot-bSel << endl;
      // }

      if (separationGain < sepGain->GetSeparationGain(sSel,bSel,sTot,bTot)
          //  &&           (mvaSC->GetBinCenter(ibin) >0 || (fCurrentMethodIdx+1)%2 )
          ){
         separationGain = sepGain->GetSeparationGain(sSel,bSel,sTot,bTot);
         //         mvaCut=mvaSC->GetBinCenter(ibin);
         mvaCut=mvaSC->GetBinLowEdge(ibin+1);
         //         if (sSel/bSel > (sTot-sSel)/(bTot-bSel)) mvaCutOrientation=-1;
         if (sSel*(bTot-bSel) > (sTot-sSel)*bSel) mvaCutOrientation=-1;
         else                                     mvaCutOrientation=1;
         sSelCut=sSel;
         bSelCut=bSel;
         //         std::cout << "new cut at " << mvaCut << "with s="<<sTot-sSel << " b="<<bTot-bSel << std::endl;
      }
      /*
        Double_t ori;
        if (sSel/bSel > (sTot-sSel)/(bTot-bSel)) ori=-1;
        else                                     ori=1;
        std::cout << ibin << " mvacut="<<mvaCut
        << " sTot=" << sTot
        << " bTot=" << bTot
        << " sSel=" << sSel
        << " bSel=" << bSel
        << " s/b(1)=" << sSel/bSel
        << " s/b(2)=" << (sTot-sSel)/(bTot-bSel)
        << " sepGain="<<sepGain->GetSeparationGain(sSel,bSel,sTot,bTot)
        << " sepGain2="<<sepGain2->GetSeparationGain(sSel,bSel,sTot,bTot)
        << "      " <<ori
        << std::endl;
      */

   }

   if (0){
      double parentIndex=sepGain->GetSeparationIndex(sTot,bTot);
      double leftIndex  =sepGain->GetSeparationIndex(sSelCut,bSelCut);
      double rightIndex  =sepGain->GetSeparationIndex(sTot-sSelCut,bTot-bSelCut);
      std::cout
         << " sTot=" << sTot
         << " bTot=" << bTot
         << " s="<<sSelCut
         << " b="<<bSelCut
         << " s2="<<(sTot-sSelCut)
         << " b2="<<(bTot-bSelCut)
         << " s/b(1)=" << sSelCut/bSelCut
         << " s/b(2)=" << (sTot-sSelCut)/(bTot-bSelCut)
         << " index before cut=" << parentIndex
         << " after: left=" << leftIndex
         << " after: right=" << rightIndex
         << " sepGain=" << parentIndex-( (sSelCut+bSelCut) * leftIndex + (sTot-sSelCut+bTot-bSelCut) * rightIndex )/(sTot+bTot)
         << " sepGain="<<separationGain
         << " sepGain="<<sepGain->GetSeparationGain(sSelCut,bSelCut,sTot,bTot)
         << " cut=" << mvaCut
         << " idx="<<fCurrentMethodIdx
         << " cutOrientation="<<mvaCutOrientation
         << std::endl;
   }
   method->SetSignalReferenceCut(mvaCut);
   method->SetSignalReferenceCutOrientation(mvaCutOrientation);

   results->GetHist("SeparationGain")->SetBinContent(fCurrentMethodIdx+1,separationGain);


   Log() << kDEBUG << "(old step) Setting method cut to " <<method->GetSignalReferenceCut()<< Endl;

   if(IsSilentFile())
   {
        mvaS ->Delete();
        mvaB ->Delete();
        mvaSC->Delete();
        mvaBC->Delete();
   }
}

////////////////////////////////////////////////////////////////////////////////

Double_t TMVA::MethodBoost::SingleBoost(MethodBase* method)
{
   Double_t returnVal=-1;


   if      (fBoostType=="AdaBoost")      returnVal = this->AdaBoost  (method,1);
   else if (fBoostType=="RealAdaBoost")  returnVal = this->AdaBoost  (method,0);
   else if (fBoostType=="Bagging")       returnVal = this->Bagging   ();
   else{
      Log() << kFATAL << "<Boost> unknown boost option " << fBoostType<< " called" << Endl;
   }
   fMethodWeight.push_back(returnVal);
   return returnVal;
}
////////////////////////////////////////////////////////////////////////////////
/// the standard (discrete or real) AdaBoost algorithm

Double_t TMVA::MethodBoost::AdaBoost(MethodBase* method, Bool_t discreteAdaBoost)
{
   if (!method) {
      Log() << kWARNING << " AdaBoost called without classifier reference - needed for calculating AdaBoost " << Endl;
      return 0;
   }

   Float_t w,v; Bool_t sig=kTRUE;
   Double_t sumAll=0, sumWrong=0;
   Bool_t* WrongDetection=new Bool_t[GetNEvents()];
   QuickMVAProbEstimator *MVAProb=NULL;

   if (discreteAdaBoost) {
      FindMVACut(method);
      Log() << kDEBUG  << " individual mva cut value = " << method->GetSignalReferenceCut() << Endl;
   } else {
      MVAProb=new TMVA::QuickMVAProbEstimator();
      // the RealAdaBoost does use a simple "yes (signal)" or "no (background)"
      // answer from your single MVA, but a "signal probability" instead (in the BDT case,
      // that would be the 'purity' in the leaf node. For some MLP parameter, the MVA output
      // can also interpreted as a probability, but here I try a general approach to get this
      // probability from the MVA distributions...

      for (Long64_t evt=0; evt<GetNEvents(); evt++) {
         const Event* ev =  Data()->GetEvent(evt);
         MVAProb->AddEvent(fMVAvalues->at(evt),ev->GetWeight(),ev->GetClass());
      }
   }


   for (Long64_t ievt=0; ievt<GetNEvents(); ievt++) WrongDetection[ievt]=kTRUE;

   // finding the wrong events and calculating their total weights
   for (Long64_t ievt=0; ievt<GetNEvents(); ievt++) {
      const Event* ev = GetEvent(ievt);
      sig=DataInfo().IsSignal(ev);
      v = fMVAvalues->at(ievt);
      w = ev->GetWeight();
      sumAll += w;
      if(!IsSilentFile())
      {
        if (fMonitorBoostedMethod) {
            if (sig) {
                fBTrainSigMVAHist[fCurrentMethodIdx]->Fill(v,w);
                fTrainSigMVAHist[fCurrentMethodIdx]->Fill(v,ev->GetOriginalWeight());
            }
            else {
                fBTrainBgdMVAHist[fCurrentMethodIdx]->Fill(v,w);
                fTrainBgdMVAHist[fCurrentMethodIdx]->Fill(v,ev->GetOriginalWeight());
            }
        }
      }

      if (discreteAdaBoost){
         if (sig  == method->IsSignalLike(fMVAvalues->at(ievt))){
            WrongDetection[ievt]=kFALSE;
         }else{
            WrongDetection[ievt]=kTRUE;
            sumWrong+=w;
         }
      }else{
         Double_t mvaProb = MVAProb->GetMVAProbAt((Float_t)fMVAvalues->at(ievt));
         mvaProb = 2*(mvaProb-0.5);
         Int_t    trueType;
         if (DataInfo().IsSignal(ev)) trueType = 1;
         else trueType = -1;
         sumWrong+= w*trueType*mvaProb;
      }
   }

   fMethodError=sumWrong/sumAll;

   // calculating the fMethodError and the boostWeight out of it uses the formula
   // w = ((1-err)/err)^beta

   Double_t boostWeight=0;

   if (fMethodError == 0) { //no misclassification made.. perfect, no boost ;)
      Log() << kWARNING << "Your classifier worked perfectly on the training sample --> serious overtraining expected and no boosting done " << Endl;
   }else{

      if (discreteAdaBoost)
         boostWeight = TMath::Log((1.-fMethodError)/fMethodError)*fAdaBoostBeta;
      else
         boostWeight = TMath::Log((1.+fMethodError)/(1-fMethodError))*fAdaBoostBeta;


      //   std::cout << "boostweight = " << boostWeight << std::endl;

      // ADA boosting, rescaling the weight of the wrong events according to the error level
      // over the entire test sample rescaling all the weights to have the same sum, but without
      // touching the original weights (changing only the boosted weight of all the events)
      // first reweight

      Double_t newSum=0., oldSum=0.;


      Double_t boostfactor = TMath::Exp(boostWeight);


      for (Long64_t ievt=0; ievt<GetNEvents(); ievt++) {
         const Event* ev =  Data()->GetEvent(ievt);
         oldSum += ev->GetWeight();
         if (discreteAdaBoost){
            // events are classified as Signal OR background .. right or wrong
            if (WrongDetection[ievt] && boostWeight != 0) {
               if (ev->GetWeight() > 0) ev->ScaleBoostWeight(boostfactor);
               else                     ev->ScaleBoostWeight(1./boostfactor);
            }
            //         if (ievt<30) std::cout<<ievt<<" var0="<<ev->GetValue(0)<<" var1="<<ev->GetValue(1)<<" weight="<<ev->GetWeight() << "  boostby:"<<boostfactor<<std::endl;

         }else{
            // events are classified by their probability of being signal or background
            // (eventually you should write this one - i.e. re-use the MVA value that were already
            // calculated and stored..   however ,for the moment ..
            Double_t mvaProb = MVAProb->GetMVAProbAt((Float_t)fMVAvalues->at(ievt));
            mvaProb = 2*(mvaProb-0.5);
            // mvaProb = (1-mvaProb);

            Int_t    trueType=1;
            if (DataInfo().IsSignal(ev)) trueType = 1;
            else trueType = -1;

            boostfactor = TMath::Exp(-1*boostWeight*trueType*mvaProb);
            if (ev->GetWeight() > 0) ev->ScaleBoostWeight(boostfactor);
            else                     ev->ScaleBoostWeight(1./boostfactor);

         }
         newSum += ev->GetWeight();
      }

      Double_t normWeight = oldSum/newSum;
      // next normalize the weights
      Double_t normSig=0, normBkg=0;
      for (Long64_t ievt=0; ievt<GetNEvents(); ievt++) {
         const Event* ev = Data()->GetEvent(ievt);
         ev->ScaleBoostWeight(normWeight);
         if (ev->GetClass()) normSig+=ev->GetWeight();
         else                normBkg+=ev->GetWeight();
      }

      Results* results = Data()->GetResults(GetMethodName(), Types::kTraining, GetAnalysisType());
      results->GetHist("SoverBtotal")->SetBinContent(fCurrentMethodIdx+1, normSig/normBkg);

      for (Long64_t ievt=0; ievt<GetNEvents(); ievt++) {
         const Event* ev = Data()->GetEvent(ievt);

         if (ev->GetClass()) ev->ScaleBoostWeight(oldSum/normSig/2);
         else                ev->ScaleBoostWeight(oldSum/normBkg/2);
      }
   }

   delete[] WrongDetection;
   if (MVAProb) delete MVAProb;

   fBoostWeight = boostWeight;  // used ONLY for the monitoring tree

   return boostWeight;
}


////////////////////////////////////////////////////////////////////////////////
/// Bagging or Bootstrap boosting, gives new random poisson weight for every event

Double_t TMVA::MethodBoost::Bagging()
{
   TRandom3  *trandom   = new TRandom3(fRandomSeed+fMethods.size());
   for (Long64_t ievt=0; ievt<GetNEvents(); ievt++) {
      const Event* ev = Data()->GetEvent(ievt);
      ev->SetBoostWeight(trandom->PoissonD(fBaggedSampleFraction));
   }
   fBoostWeight = 1; // used ONLY for the monitoring tree
   return 1.;
}


////////////////////////////////////////////////////////////////////////////////
/// Get help message text
///
/// typical length of text line:
///         "|--------------------------------------------------------------|"

void TMVA::MethodBoost::GetHelpMessage() const
{
   Log() << Endl;
   Log() << gTools().Color("bold") << "--- Short description:" << gTools().Color("reset") << Endl;
   Log() << Endl;
   Log() << "This method combines several classifier of one species in a "<<Endl;
   Log() << "single multivariate quantity via the boost algorithm." << Endl;
   Log() << "the output is a weighted sum over all individual classifiers" <<Endl;
   Log() << "By default, the AdaBoost method is employed, which gives " << Endl;
   Log() << "events that were misclassified in the previous tree a larger " << Endl;
   Log() << "weight in the training of the following classifier."<<Endl;
   Log() << "Optionally, Bagged boosting can also be applied." << Endl;
   Log() << Endl;
   Log() << gTools().Color("bold") << "--- Performance tuning via configuration options:" << gTools().Color("reset") << Endl;
   Log() << Endl;
   Log() << "The most important parameter in the configuration is the "<<Endl;
   Log() << "number of boosts applied (Boost_Num) and the choice of boosting"<<Endl;
   Log() << "(Boost_Type), which can be set to either AdaBoost or Bagging." << Endl;
   Log() << "AdaBoosting: The most important parameters in this configuration" <<Endl;
   Log() << "is the beta parameter (Boost_AdaBoostBeta)  " << Endl;
   Log() << "When boosting a linear classifier, it is sometimes advantageous"<<Endl;
   Log() << "to transform the MVA output non-linearly. The following options" <<Endl;
   Log() << "are available: step, log, and minmax, the default is no transform."<<Endl;
   Log() <<Endl;
   Log() << "Some classifiers are hard to boost and do not improve much in"<<Endl;
   Log() << "their performance by boosting them, some even slightly deteriorate"<< Endl;
   Log() << "due to the boosting." <<Endl;
   Log() << "The booking of the boost method is special since it requires"<<Endl;
   Log() << "the booing of the method to be boosted and the boost itself."<<Endl;
   Log() << "This is solved by booking the method to be boosted and to add"<<Endl;
   Log() << "all Boost parameters, which all begin with \"Boost_\" to the"<<Endl;
   Log() << "options string. The factory separates the options and initiates"<<Endl;
   Log() << "the boost process. The TMVA macro directory contains the example"<<Endl;
   Log() << "macro \"Boost.C\"" <<Endl;
}

////////////////////////////////////////////////////////////////////////////////

const TMVA::Ranking* TMVA::MethodBoost::CreateRanking()
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// return boosted MVA response

Double_t TMVA::MethodBoost::GetMvaValue( Double_t* err, Double_t* errUpper )
{
   Double_t mvaValue = 0;
   Double_t norm = 0;
   Double_t epsilon = TMath::Exp(-1.);
   //Double_t fact    = TMath::Exp(-1.)+TMath::Exp(1.);
   for (UInt_t i=0;i< fMethods.size(); i++){
      MethodBase* m = dynamic_cast<MethodBase*>(fMethods[i]);
      if (m==0) continue;
      Double_t val = fTmpEvent ? m->GetMvaValue(fTmpEvent) : m->GetMvaValue();
      Double_t sigcut = m->GetSignalReferenceCut();

      // default is no transform
      if (fTransformString == "linear"){

      }
      else if (fTransformString == "log"){
         if (val < sigcut) val = sigcut;

         val = TMath::Log((val-sigcut)+epsilon);
      }
      else if (fTransformString == "step" ){
         if (m->IsSignalLike(val)) val = 1.;
         else val = -1.;
      }
      else if (fTransformString == "gauss"){
         val = TMath::Gaus((val-sigcut),1);
      }
      else {
         Log() << kFATAL << "error unknown transformation " << fTransformString<<Endl;
      }
      mvaValue+=val*fMethodWeight[i];
      norm    +=fMethodWeight[i];
      //      std::cout << "mva("<<i<<") = "<<val<<" " << valx<< " " << mvaValue<<"  and sigcut="<<sigcut << std::endl;
   }
   mvaValue/=norm;
   // cannot determine error
   NoErrorCalc(err, errUpper);

   return mvaValue;
}

////////////////////////////////////////////////////////////////////////////////
/// Calculate the ROC integral of a single classifier or even the
/// whole boosted classifier.  The tree type (training or testing
/// sample) is specified by 'eTT'.
///
/// If tree type kTraining is set, the original training sample is
/// used to compute the ROC integral (original weights).
///
/// - singleMethod - if kTRUE, return ROC integral of single (last
///                  trained) classifier; if kFALSE, return ROC
///                  integral of full classifier
///
/// - eTT - tree type (Types::kTraining / Types::kTesting)
///
/// - CalcOverlapIntergral - if kTRUE, the overlap integral of the
///                          signal/background MVA distributions
///                          is calculated and stored in
///                          'fOverlap_integral'

Double_t TMVA::MethodBoost::GetBoostROCIntegral(Bool_t singleMethod, Types::ETreeType eTT, Bool_t CalcOverlapIntergral)
{
   // set data sample training / testing
   Data()->SetCurrentType(eTT);

   MethodBase* method = singleMethod ? dynamic_cast<MethodBase*>(fMethods.back()) : 0; // ToDo CoVerity flags this line as there is no protection against a zero-pointer delivered by dynamic_cast
   // to make CoVerity happy (although, OF COURSE, the last method in the committee
   // has to be also of type MethodBase as ANY method is... hence the dynamic_cast
   // will never by "zero" ...
   if (singleMethod && !method) {
      Log() << kFATAL << " What do you do? Your method:"
            << fMethods.back()->GetName()
            << " seems not to be a propper TMVA method"
            << Endl;
      std::exit(1);
   }
   Double_t err = 0.0;

   // temporary renormalize the method weights in case of evaluation
   // of full classifier.
   // save the old normalization of the methods
   std::vector<Double_t> OldMethodWeight(fMethodWeight);
   if (!singleMethod) {
      // calculate sum of weights of all methods
      Double_t AllMethodsWeight = 0;
      for (UInt_t i=0; i<=fCurrentMethodIdx; i++)
         AllMethodsWeight += fMethodWeight.at(i);
      // normalize the weights of the classifiers
      if (AllMethodsWeight != 0.0) {
         for (UInt_t i=0; i<=fCurrentMethodIdx; i++)
            fMethodWeight[i] /= AllMethodsWeight;
      }
   }

   // calculate MVA values
   Double_t meanS, meanB, rmsS, rmsB, xmin, xmax, nrms = 10;
   std::vector <Float_t>* mvaRes;
   if (singleMethod && eTT==Types::kTraining)
      mvaRes = fMVAvalues; // values already calculated
   else {
      mvaRes = new std::vector <Float_t>(GetNEvents());
      for (Long64_t ievt=0; ievt<GetNEvents(); ievt++) {
         GetEvent(ievt);
         (*mvaRes)[ievt] = singleMethod ? method->GetMvaValue(&err) : GetMvaValue(&err);
      }
   }

   // restore the method weights
   if (!singleMethod)
      fMethodWeight = OldMethodWeight;

   // now create histograms for calculation of the ROC integral
   Int_t signalClass = 0;
   if (DataInfo().GetClassInfo("Signal") != 0) {
      signalClass = DataInfo().GetClassInfo("Signal")->GetNumber();
   }
   gTools().ComputeStat( GetEventCollection(eTT), mvaRes,
                         meanS, meanB, rmsS, rmsB, xmin, xmax, signalClass );

   fNbins = gConfig().fVariablePlotting.fNbinsXOfROCCurve;
   xmin = TMath::Max( TMath::Min(meanS - nrms*rmsS, meanB - nrms*rmsB ), xmin );
   xmax = TMath::Min( TMath::Max(meanS + nrms*rmsS, meanB + nrms*rmsB ), xmax ) + 0.0001;

   // calculate ROC integral
   TH1* mva_s = new TH1F( "MVA_S", "MVA_S", fNbins, xmin, xmax );
   TH1* mva_b = new TH1F( "MVA_B", "MVA_B", fNbins, xmin, xmax );
   TH1 *mva_s_overlap=0, *mva_b_overlap=0;
   if (CalcOverlapIntergral) {
      mva_s_overlap = new TH1F( "MVA_S_OVERLAP", "MVA_S_OVERLAP", fNbins, xmin, xmax );
      mva_b_overlap = new TH1F( "MVA_B_OVERLAP", "MVA_B_OVERLAP", fNbins, xmin, xmax );
   }
   for (Long64_t ievt=0; ievt<GetNEvents(); ievt++) {
      const Event* ev = GetEvent(ievt);
      Float_t w = (eTT==Types::kTesting ? ev->GetWeight() : ev->GetOriginalWeight());
      if (DataInfo().IsSignal(ev))  mva_s->Fill( (*mvaRes)[ievt], w );
      else                          mva_b->Fill( (*mvaRes)[ievt], w );

      if (CalcOverlapIntergral) {
         Float_t w_ov = ev->GetWeight();
         if (DataInfo().IsSignal(ev))
            mva_s_overlap->Fill( (*mvaRes)[ievt], w_ov );
         else
            mva_b_overlap->Fill( (*mvaRes)[ievt], w_ov );
      }
   }
   gTools().NormHist( mva_s );
   gTools().NormHist( mva_b );
   PDF *fS = new PDF( "PDF Sig", mva_s, PDF::kSpline2 );
   PDF *fB = new PDF( "PDF Bkg", mva_b, PDF::kSpline2 );

   // calculate ROC integral from fS, fB
   Double_t ROC = MethodBase::GetROCIntegral(fS, fB);

   // calculate overlap integral
   if (CalcOverlapIntergral) {
      gTools().NormHist( mva_s_overlap );
      gTools().NormHist( mva_b_overlap );

      fOverlap_integral = 0.0;
      for (Int_t bin=1; bin<=mva_s_overlap->GetNbinsX(); bin++){
         Double_t bc_s = mva_s_overlap->GetBinContent(bin);
         Double_t bc_b = mva_b_overlap->GetBinContent(bin);
         if (bc_s > 0.0 && bc_b > 0.0)
            fOverlap_integral += TMath::Min(bc_s, bc_b);
      }

      delete mva_s_overlap;
      delete mva_b_overlap;
   }

   delete mva_s;
   delete mva_b;
   delete fS;
   delete fB;
   if (!(singleMethod && eTT==Types::kTraining))  delete mvaRes;

   Data()->SetCurrentType(Types::kTraining);

   return ROC;
}

void TMVA::MethodBoost::CalcMVAValues()
{
   // Calculate MVA values of current method fMethods.back() on
   // training sample

   Data()->SetCurrentType(Types::kTraining);
   MethodBase* method = dynamic_cast<MethodBase*>(fMethods.back());
   if (!method) {
      Log() << kFATAL << "dynamic cast to MethodBase* failed" <<Endl;
      return;
   }
   // calculate MVA values
   for (Long64_t ievt=0; ievt<GetNEvents(); ievt++) {
      GetEvent(ievt);
      fMVAvalues->at(ievt) = method->GetMvaValue();
   }

   // fill cumulative mva distribution


}

////////////////////////////////////////////////////////////////////////////////
/// fill various monitoring histograms from information of the individual classifiers that
/// have been boosted.
/// of course.... this depends very much on the individual classifiers, and so far, only for
/// Decision Trees, this monitoring is actually implemented

void TMVA::MethodBoost::MonitorBoost( Types::EBoostStage stage , UInt_t methodIndex )
{
   Results* results = Data()->GetResults(GetMethodName(), Types::kTraining, GetAnalysisType());

   if (GetCurrentMethod(methodIndex)->GetMethodType() == TMVA::Types::kDT) {
      TMVA::MethodDT* currentDT=dynamic_cast<TMVA::MethodDT*>(GetCurrentMethod(methodIndex));
      if (currentDT){
         if (stage == Types::kBoostProcBegin){
            results->Store(new TH1I("NodesBeforePruning","nodes before pruning",this->GetBoostNum(),0,this->GetBoostNum()),"NodesBeforePruning");
            results->Store(new TH1I("NodesAfterPruning","nodes after pruning",this->GetBoostNum(),0,this->GetBoostNum()),"NodesAfterPruning");
         }

         if (stage == Types::kBeforeTraining){
         }
         else if (stage == Types::kBeforeBoosting){
            results->GetHist("NodesBeforePruning")->SetBinContent(methodIndex+1,currentDT->GetNNodesBeforePruning());
            results->GetHist("NodesAfterPruning")->SetBinContent(methodIndex+1,currentDT->GetNNodes());
         }
         else if (stage == Types::kAfterBoosting){

         }
         else if (stage != Types::kBoostProcEnd){
            Log() << kINFO << "<Train> average number of nodes before/after pruning : "
                  <<   results->GetHist("NodesBeforePruning")->GetMean() << " / "
                  <<   results->GetHist("NodesAfterPruning")->GetMean()
                  << Endl;
         }
      }

   }else if (GetCurrentMethod(methodIndex)->GetMethodType() == TMVA::Types::kFisher) {
      if (stage == Types::kAfterBoosting){
         TMVA::MsgLogger::EnableOutput();
      }
   }else{
      if (methodIndex < 3){
         Log() << kDEBUG << "No detailed boost monitoring for "
               << GetCurrentMethod(methodIndex)->GetMethodName()
               << " yet available " << Endl;
      }
   }

   //boosting plots universal for all classifiers 'typically for debug purposes only as they are not general enough'

   if (stage == Types::kBeforeBoosting){
      // if you want to display the weighted events for 2D case at each boost step:
      if (fDetailedMonitoring){
         // the following code is useful only for 2D examples - mainly illustration for debug/educational purposes:
         if (DataInfo().GetNVariables() == 2) {
            results->Store(new TH2F(Form("EventDistSig_%d",methodIndex),Form("EventDistSig_%d",methodIndex),100,0,7,100,0,7));
            results->GetHist(Form("EventDistSig_%d",methodIndex))->SetMarkerColor(4);
            results->Store(new TH2F(Form("EventDistBkg_%d",methodIndex),Form("EventDistBkg_%d",methodIndex),100,0,7,100,0,7));
            results->GetHist(Form("EventDistBkg_%d",methodIndex))->SetMarkerColor(2);

            Data()->SetCurrentType(Types::kTraining);
            for (Long64_t ievt=0; ievt<GetNEvents(); ievt++) {
               const Event* ev = GetEvent(ievt);
               Float_t w = ev->GetWeight();
               Float_t v0= ev->GetValue(0);
               Float_t v1= ev->GetValue(1);
               //         if (ievt<3) std::cout<<ievt<<" var0="<<v0<<" var1="<<v1<<" weight="<<w<<std::endl;
               TH2* h;
               if (DataInfo().IsSignal(ev)) h=results->GetHist2D(Form("EventDistSig_%d",methodIndex));
               else                         h=results->GetHist2D(Form("EventDistBkg_%d",methodIndex));
               if (h) h->Fill(v0,v1,w);
            }
         }
      }
   }

   return;
}


