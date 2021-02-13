// @(#)Root/tmva $Id$
// Author: Andreas Hoecker, Peter Speckmayer, Joerg Stelzer, Helge Voss, Kai Voss, Eckhard von Toerne, Jan Therhaag
// Updated by: Omar Zapata, Kim Albertsson
/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : Factory                                                               *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Implementation (see header for description)                               *
 *                                                                                *
 * Authors :                                                                      *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
 *      Joerg Stelzer   <stelzer@cern.ch>        - DESY, Germany                  *
 *      Peter Speckmayer <peter.speckmayer@cern.ch> - CERN, Switzerland           *
 *      Jan Therhaag          <Jan.Therhaag@cern.ch>   - U of Bonn, Germany       *
 *      Eckhard v. Toerne     <evt@uni-bonn.de>        - U of Bonn, Germany       *
 *      Helge Voss      <Helge.Voss@cern.ch>     - MPI-K Heidelberg, Germany      *
 *      Kai Voss        <Kai.Voss@cern.ch>       - U. of Victoria, Canada         *
 *      Omar Zapata     <Omar.Zapata@cern.ch>    - UdeA/ITM Colombia              *
 *      Lorenzo Moneta  <Lorenzo.Moneta@cern.ch> - CERN, Switzerland              *
 *      Sergei Gleyzer  <Sergei.Gleyzer@cern.ch> - U of Florida & CERN            *
 *      Kim Albertsson  <kim.albertsson@cern.ch> - LTU & CERN                     *
 *                                                                                *
 * Copyright (c) 2005-2015:                                                       *
 *      CERN, Switzerland                                                         *
 *      U. of Victoria, Canada                                                    *
 *      MPI-K Heidelberg, Germany                                                 *
 *      U. of Bonn, Germany                                                       *
 *      UdeA/ITM, Colombia                                                        *
 *      U. of Florida, USA                                                        *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

/*! \class TMVA::Factory
\ingroup TMVA

This is the main MVA steering class.
It creates all MVA methods, and guides them through the training, testing and
evaluation phases.
*/

#include "TMVA/Factory.h"

#include "TMVA/ClassifierFactory.h"
#include "TMVA/Config.h"
#include "TMVA/Configurable.h"
#include "TMVA/Tools.h"
#include "TMVA/Ranking.h"
#include "TMVA/DataSet.h"
#include "TMVA/IMethod.h"
#include "TMVA/MethodBase.h"
#include "TMVA/DataInputHandler.h"
#include "TMVA/DataSetManager.h"
#include "TMVA/DataSetInfo.h"
#include "TMVA/DataLoader.h"
#include "TMVA/MethodBoost.h"
#include "TMVA/MethodCategory.h"
#include "TMVA/ROCCalc.h"
#include "TMVA/ROCCurve.h"
#include "TMVA/MsgLogger.h"

#include "TMVA/VariableInfo.h"
#include "TMVA/VariableTransform.h"

#include "TMVA/Results.h"
#include "TMVA/ResultsClassification.h"
#include "TMVA/ResultsRegression.h"
#include "TMVA/ResultsMulticlass.h"
#include <list>
#include <bitset>
#include <set>

#include "TMVA/Types.h"

#include "TROOT.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TEventList.h"
#include "TH2.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TMatrixF.h"
#include "TMatrixDSym.h"
#include "TMultiGraph.h"
#include "TPrincipal.h"
#include "TMath.h"
#include "TSystem.h"
#include "TCanvas.h"

const Int_t MinNoTrainingEvents = 10;
// const Int_t  MinNoTestEvents     = 1;

ClassImp(TMVA::Factory);

#define READXML kTRUE

// number of bits for bitset
#define VIBITS 32

////////////////////////////////////////////////////////////////////////////////
/// Standard constructor.
///
///  - jobname       : this name will appear in all weight file names produced by the MVAs
///  - theTargetFile : output ROOT file; the test tree and all evaluation plots
///                   will be stored here
///  - theOption     : option string; currently: "V" for verbose

TMVA::Factory::Factory(TString jobName, TFile *theTargetFile, TString theOption)
   : Configurable(theOption), fTransformations("I"), fVerbose(kFALSE), fVerboseLevel(kINFO), fCorrelations(kFALSE),
     fROC(kTRUE), fSilentFile(theTargetFile == nullptr), fJobName(jobName), fAnalysisType(Types::kClassification),
     fModelPersistence(kTRUE)
{
   fName = "Factory";
   fgTargetFile = theTargetFile;
   fLogger->SetSource(fName.Data());

   // render silent
   if (gTools().CheckForSilentOption(GetOptions()))
      Log().InhibitOutput(); // make sure is silent if wanted to

   // init configurable
   SetConfigDescription("Configuration options for Factory running");
   SetConfigName(GetName());

   // histograms are not automatically associated with the current
   // directory and hence don't go out of scope when closing the file
   // TH1::AddDirectory(kFALSE);
   Bool_t silent = kFALSE;
#ifdef WIN32
   // under Windows, switch progress bar and color off by default, as the typical windows shell doesn't handle these
   // (would need different sequences..)
   Bool_t color = kFALSE;
   Bool_t drawProgressBar = kFALSE;
#else
   Bool_t color = !gROOT->IsBatch();
   Bool_t drawProgressBar = kTRUE;
#endif
   DeclareOptionRef(fVerbose, "V", "Verbose flag");
   DeclareOptionRef(fVerboseLevel = TString("Info"), "VerboseLevel", "VerboseLevel (Debug/Verbose/Info)");
   AddPreDefVal(TString("Debug"));
   AddPreDefVal(TString("Verbose"));
   AddPreDefVal(TString("Info"));
   DeclareOptionRef(color, "Color", "Flag for coloured screen output (default: True, if in batch mode: False)");
   DeclareOptionRef(
      fTransformations, "Transformations",
      "List of transformations to test; formatting example: \"Transformations=I;D;P;U;G,D\", for identity, "
      "decorrelation, PCA, Uniform and Gaussianisation followed by decorrelation transformations");
   DeclareOptionRef(fCorrelations, "Correlations", "boolean to show correlation in output");
   DeclareOptionRef(fROC, "ROC", "boolean to show ROC in output");
   DeclareOptionRef(silent, "Silent",
                    "Batch mode: boolean silent flag inhibiting any output from TMVA after the creation of the factory "
                    "class object (default: False)");
   DeclareOptionRef(drawProgressBar, "DrawProgressBar",
                    "Draw progress bar to display training, testing and evaluation schedule (default: True)");
   DeclareOptionRef(fModelPersistence, "ModelPersistence",
                    "Option to save the trained model in xml file or using serialization");

   TString analysisType("Auto");
   DeclareOptionRef(analysisType, "AnalysisType",
                    "Set the analysis type (Classification, Regression, Multiclass, Auto) (default: Auto)");
   AddPreDefVal(TString("Classification"));
   AddPreDefVal(TString("Regression"));
   AddPreDefVal(TString("Multiclass"));
   AddPreDefVal(TString("Auto"));

   ParseOptions();
   CheckForUnusedOptions();

   if (Verbose())
      fLogger->SetMinType(kVERBOSE);
   if (fVerboseLevel.CompareTo("Debug") == 0)
      fLogger->SetMinType(kDEBUG);
   if (fVerboseLevel.CompareTo("Verbose") == 0)
      fLogger->SetMinType(kVERBOSE);
   if (fVerboseLevel.CompareTo("Info") == 0)
      fLogger->SetMinType(kINFO);

   // global settings
   gConfig().SetUseColor(color);
   gConfig().SetSilent(silent);
   gConfig().SetDrawProgressBar(drawProgressBar);

   analysisType.ToLower();
   if (analysisType == "classification")
      fAnalysisType = Types::kClassification;
   else if (analysisType == "regression")
      fAnalysisType = Types::kRegression;
   else if (analysisType == "multiclass")
      fAnalysisType = Types::kMulticlass;
   else if (analysisType == "auto")
      fAnalysisType = Types::kNoAnalysisType;

   //   Greetings();
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TMVA::Factory::Factory(TString jobName, TString theOption)
   : Configurable(theOption), fTransformations("I"), fVerbose(kFALSE), fCorrelations(kFALSE), fROC(kTRUE),
     fSilentFile(kTRUE), fJobName(jobName), fAnalysisType(Types::kClassification), fModelPersistence(kTRUE)
{
   fName = "Factory";
   fgTargetFile = nullptr;
   fLogger->SetSource(fName.Data());

   // render silent
   if (gTools().CheckForSilentOption(GetOptions()))
      Log().InhibitOutput(); // make sure is silent if wanted to

   // init configurable
   SetConfigDescription("Configuration options for Factory running");
   SetConfigName(GetName());

   // histograms are not automatically associated with the current
   // directory and hence don't go out of scope when closing the file
   TH1::AddDirectory(kFALSE);
   Bool_t silent = kFALSE;
#ifdef WIN32
   // under Windows, switch progress bar and color off by default, as the typical windows shell doesn't handle these
   // (would need different sequences..)
   Bool_t color = kFALSE;
   Bool_t drawProgressBar = kFALSE;
#else
   Bool_t color = !gROOT->IsBatch();
   Bool_t drawProgressBar = kTRUE;
#endif
   DeclareOptionRef(fVerbose, "V", "Verbose flag");
   DeclareOptionRef(fVerboseLevel = TString("Info"), "VerboseLevel", "VerboseLevel (Debug/Verbose/Info)");
   AddPreDefVal(TString("Debug"));
   AddPreDefVal(TString("Verbose"));
   AddPreDefVal(TString("Info"));
   DeclareOptionRef(color, "Color", "Flag for coloured screen output (default: True, if in batch mode: False)");
   DeclareOptionRef(
      fTransformations, "Transformations",
      "List of transformations to test; formatting example: \"Transformations=I;D;P;U;G,D\", for identity, "
      "decorrelation, PCA, Uniform and Gaussianisation followed by decorrelation transformations");
   DeclareOptionRef(fCorrelations, "Correlations", "boolean to show correlation in output");
   DeclareOptionRef(fROC, "ROC", "boolean to show ROC in output");
   DeclareOptionRef(silent, "Silent",
                    "Batch mode: boolean silent flag inhibiting any output from TMVA after the creation of the factory "
                    "class object (default: False)");
   DeclareOptionRef(drawProgressBar, "DrawProgressBar",
                    "Draw progress bar to display training, testing and evaluation schedule (default: True)");
   DeclareOptionRef(fModelPersistence, "ModelPersistence",
                    "Option to save the trained model in xml file or using serialization");

   TString analysisType("Auto");
   DeclareOptionRef(analysisType, "AnalysisType",
                    "Set the analysis type (Classification, Regression, Multiclass, Auto) (default: Auto)");
   AddPreDefVal(TString("Classification"));
   AddPreDefVal(TString("Regression"));
   AddPreDefVal(TString("Multiclass"));
   AddPreDefVal(TString("Auto"));

   ParseOptions();
   CheckForUnusedOptions();

   if (Verbose())
      fLogger->SetMinType(kVERBOSE);
   if (fVerboseLevel.CompareTo("Debug") == 0)
      fLogger->SetMinType(kDEBUG);
   if (fVerboseLevel.CompareTo("Verbose") == 0)
      fLogger->SetMinType(kVERBOSE);
   if (fVerboseLevel.CompareTo("Info") == 0)
      fLogger->SetMinType(kINFO);

   // global settings
   gConfig().SetUseColor(color);
   gConfig().SetSilent(silent);
   gConfig().SetDrawProgressBar(drawProgressBar);

   analysisType.ToLower();
   if (analysisType == "classification")
      fAnalysisType = Types::kClassification;
   else if (analysisType == "regression")
      fAnalysisType = Types::kRegression;
   else if (analysisType == "multiclass")
      fAnalysisType = Types::kMulticlass;
   else if (analysisType == "auto")
      fAnalysisType = Types::kNoAnalysisType;

   Greetings();
}

////////////////////////////////////////////////////////////////////////////////
/// Print welcome message.
/// Options are: kLogoWelcomeMsg, kIsometricWelcomeMsg, kLeanWelcomeMsg

void TMVA::Factory::Greetings()
{
   gTools().ROOTVersionMessage(Log());
   gTools().TMVAWelcomeMessage(Log(), gTools().kLogoWelcomeMsg);
   gTools().TMVAVersionMessage(Log());
   Log() << Endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TMVA::Factory::~Factory(void)
{
   std::vector<TMVA::VariableTransformBase *>::iterator trfIt = fDefaultTrfs.begin();
   for (; trfIt != fDefaultTrfs.end(); ++trfIt)
      delete (*trfIt);

   this->DeleteAllMethods();

   // problem with call of REGISTER_METHOD macro ...
   //   ClassifierFactory::DestroyInstance();
   //   Types::DestroyInstance();
   // Tools::DestroyInstance();
   // Config::DestroyInstance();
}

////////////////////////////////////////////////////////////////////////////////
/// Delete methods.

void TMVA::Factory::DeleteAllMethods(void)
{
   std::map<TString, MVector *>::iterator itrMap;

   for (itrMap = fMethodsMap.begin(); itrMap != fMethodsMap.end(); ++itrMap) {
      MVector *methods = itrMap->second;
      // delete methods
      MVector::iterator itrMethod = methods->begin();
      for (; itrMethod != methods->end(); ++itrMethod) {
         Log() << kDEBUG << "Delete method: " << (*itrMethod)->GetName() << Endl;
         delete (*itrMethod);
      }
      methods->clear();
      delete methods;
   }
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::Factory::SetVerbose(Bool_t v)
{
   fVerbose = v;
}

////////////////////////////////////////////////////////////////////////////////
/// Book a classifier or regression method.

TMVA::MethodBase *
TMVA::Factory::BookMethod(TMVA::DataLoader *loader, TString theMethodName, TString methodTitle, TString theOption)
{
   if (fModelPersistence)
      gSystem->MakeDirectory(loader->GetName()); // creating directory for DataLoader output

   TString datasetname = loader->GetName();

   if (fAnalysisType == Types::kNoAnalysisType) {
      if (loader->GetDataSetInfo().GetNClasses() == 2 && loader->GetDataSetInfo().GetClassInfo("Signal") != NULL &&
          loader->GetDataSetInfo().GetClassInfo("Background") != NULL) {
         fAnalysisType = Types::kClassification; // default is classification
      } else if (loader->GetDataSetInfo().GetNClasses() >= 2) {
         fAnalysisType = Types::kMulticlass; // if two classes, but not named "Signal" and "Background"
      } else
         Log() << kFATAL << "No analysis type for " << loader->GetDataSetInfo().GetNClasses() << " classes and "
               << loader->GetDataSetInfo().GetNTargets() << " regression targets." << Endl;
   }

   // booking via name; the names are translated into enums and the
   // corresponding overloaded BookMethod is called

   if (fMethodsMap.find(datasetname) != fMethodsMap.end()) {
      if (GetMethod(datasetname, methodTitle) != 0) {
         Log() << kFATAL << "Booking failed since method with title <" << methodTitle << "> already exists "
               << "in with DataSet Name <" << loader->GetName() << ">  " << Endl;
      }
   }

   Log() << kHEADER << "Booking method: " << gTools().Color("bold")
         << methodTitle
         // << gTools().Color("reset")<<" DataSet Name: "<<gTools().Color("bold")<<loader->GetName()
         << gTools().Color("reset") << Endl << Endl;

   // interpret option string with respect to a request for boosting (i.e., BostNum > 0)
   Int_t boostNum = 0;
   TMVA::Configurable *conf = new TMVA::Configurable(theOption);
   conf->DeclareOptionRef(boostNum = 0, "Boost_num", "Number of times the classifier will be boosted");
   conf->ParseOptions();
   delete conf;
   // this is name of weight file directory
   TString fileDir;
   if (fModelPersistence) {
      // find prefix in fWeightFileDir;
      TString prefix = gConfig().GetIONames().fWeightFileDirPrefix;
      fileDir = prefix;
      if (!prefix.IsNull())
         if (fileDir[fileDir.Length() - 1] != '/')
            fileDir += "/";
      fileDir += loader->GetName();
      fileDir += "/" + gConfig().GetIONames().fWeightFileDir;
   }
   // initialize methods
   IMethod *im;
   if (!boostNum) {
      im = ClassifierFactory::Instance().Create(theMethodName.Data(), fJobName, methodTitle, loader->GetDataSetInfo(),
                                                theOption);
   } else {
      // boosted classifier, requires a specific definition, making it transparent for the user
      Log() << kDEBUG << "Boost Number is " << boostNum << " > 0: train boosted classifier" << Endl;
      im = ClassifierFactory::Instance().Create("Boost", fJobName, methodTitle, loader->GetDataSetInfo(), theOption);
      MethodBoost *methBoost = dynamic_cast<MethodBoost *>(im); // DSMTEST divided into two lines
      if (!methBoost) {                                         // DSMTEST
         Log() << kFATAL << "Method with type kBoost cannot be casted to MethodCategory. /Factory" << Endl; // DSMTEST
         return nullptr;
      }
      if (fModelPersistence)
         methBoost->SetWeightFileDir(fileDir);
      methBoost->SetModelPersistence(fModelPersistence);
      methBoost->SetBoostedMethodName(theMethodName);                            // DSMTEST divided into two lines
      methBoost->fDataSetManager = loader->GetDataSetInfo().GetDataSetManager(); // DSMTEST
      methBoost->SetFile(fgTargetFile);
      methBoost->SetSilentFile(IsSilentFile());
   }

   MethodBase *method = dynamic_cast<MethodBase *>(im);
   if (method == 0)
      return 0; // could not create method

   // set fDataSetManager if MethodCategory (to enable Category to create datasetinfo objects) // DSMTEST
   if (method->GetMethodType() == Types::kCategory) {                 // DSMTEST
      MethodCategory *methCat = (dynamic_cast<MethodCategory *>(im)); // DSMTEST
      if (!methCat) {                                                 // DSMTEST
         Log() << kFATAL << "Method with type kCategory cannot be casted to MethodCategory. /Factory"
               << Endl; // DSMTEST
         return nullptr;
      }
      if (fModelPersistence)
         methCat->SetWeightFileDir(fileDir);
      methCat->SetModelPersistence(fModelPersistence);
      methCat->fDataSetManager = loader->GetDataSetInfo().GetDataSetManager(); // DSMTEST
      methCat->SetFile(fgTargetFile);
      methCat->SetSilentFile(IsSilentFile());
   } // DSMTEST

   if (!method->HasAnalysisType(fAnalysisType, loader->GetDataSetInfo().GetNClasses(),
                                loader->GetDataSetInfo().GetNTargets())) {
      Log() << kWARNING << "Method " << method->GetMethodTypeName() << " is not capable of handling ";
      if (fAnalysisType == Types::kRegression) {
         Log() << "regression with " << loader->GetDataSetInfo().GetNTargets() << " targets." << Endl;
      } else if (fAnalysisType == Types::kMulticlass) {
         Log() << "multiclass classification with " << loader->GetDataSetInfo().GetNClasses() << " classes." << Endl;
      } else {
         Log() << "classification with " << loader->GetDataSetInfo().GetNClasses() << " classes." << Endl;
      }
      return 0;
   }

   if (fModelPersistence)
      method->SetWeightFileDir(fileDir);
   method->SetModelPersistence(fModelPersistence);
   method->SetAnalysisType(fAnalysisType);
   method->SetupMethod();
   method->ParseOptions();
   method->ProcessSetup();
   method->SetFile(fgTargetFile);
   method->SetSilentFile(IsSilentFile());

   // check-for-unused-options is performed; may be overridden by derived classes
   method->CheckSetup();

   if (fMethodsMap.find(datasetname) == fMethodsMap.end()) {
      MVector *mvector = new MVector;
      fMethodsMap[datasetname] = mvector;
   }
   fMethodsMap[datasetname]->push_back(method);
   return method;
}

////////////////////////////////////////////////////////////////////////////////
/// Books MVA method. The option configuration string is custom for each MVA
/// the TString field "theNameAppendix" serves to define (and distinguish)
/// several instances of a given MVA, eg, when one wants to compare the
/// performance of various configurations

TMVA::MethodBase *
TMVA::Factory::BookMethod(TMVA::DataLoader *loader, Types::EMVA theMethod, TString methodTitle, TString theOption)
{
   return BookMethod(loader, Types::Instance().GetMethodName(theMethod), methodTitle, theOption);
}

////////////////////////////////////////////////////////////////////////////////
/// Adds an already constructed method to be managed by this factory.
///
/// \note Private.
/// \note Know what you are doing when using this method. The method that you
/// are loading could be trained already.
///

TMVA::MethodBase *
TMVA::Factory::BookMethodWeightfile(DataLoader *loader, TMVA::Types::EMVA methodType, const TString &weightfile)
{
   TString datasetname = loader->GetName();
   std::string methodTypeName = std::string(Types::Instance().GetMethodName(methodType).Data());
   DataSetInfo &dsi = loader->GetDataSetInfo();

   IMethod *im = ClassifierFactory::Instance().Create(methodTypeName, dsi, weightfile);
   MethodBase *method = (dynamic_cast<MethodBase *>(im));

   if (method == nullptr)
      return nullptr;

   if (method->GetMethodType() == Types::kCategory) {
      Log() << kERROR << "Cannot handle category methods for now." << Endl;
   }

   TString fileDir;
   if (fModelPersistence) {
      // find prefix in fWeightFileDir;
      TString prefix = gConfig().GetIONames().fWeightFileDirPrefix;
      fileDir = prefix;
      if (!prefix.IsNull())
         if (fileDir[fileDir.Length() - 1] != '/')
            fileDir += "/";
      fileDir = loader->GetName();
      fileDir += "/" + gConfig().GetIONames().fWeightFileDir;
   }

   if (fModelPersistence)
      method->SetWeightFileDir(fileDir);
   method->SetModelPersistence(fModelPersistence);
   method->SetAnalysisType(fAnalysisType);
   method->SetupMethod();
   method->SetFile(fgTargetFile);
   method->SetSilentFile(IsSilentFile());

   method->DeclareCompatibilityOptions();

   // read weight file
   method->ReadStateFromFile();

   // method->CheckSetup();

   TString methodTitle = method->GetName();
   if (HasMethod(datasetname, methodTitle) != 0) {
      Log() << kFATAL << "Booking failed since method with title <" << methodTitle << "> already exists "
            << "in with DataSet Name <" << loader->GetName() << ">  " << Endl;
   }

   Log() << kINFO << "Booked classifier \"" << method->GetMethodName() << "\" of type: \""
         << method->GetMethodTypeName() << "\"" << Endl;

   if (fMethodsMap.count(datasetname) == 0) {
      MVector *mvector = new MVector;
      fMethodsMap[datasetname] = mvector;
   }

   fMethodsMap[datasetname]->push_back(method);

   return method;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns pointer to MVA that corresponds to given method title.

TMVA::IMethod *TMVA::Factory::GetMethod(const TString &datasetname, const TString &methodTitle) const
{
   if (fMethodsMap.find(datasetname) == fMethodsMap.end())
      return 0;

   MVector *methods = fMethodsMap.find(datasetname)->second;

   MVector::const_iterator itrMethod;
   //
   for (itrMethod = methods->begin(); itrMethod != methods->end(); ++itrMethod) {
      MethodBase *mva = dynamic_cast<MethodBase *>(*itrMethod);
      if ((mva->GetMethodName()) == methodTitle)
         return mva;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Checks whether a given method name is defined for a given dataset.

Bool_t TMVA::Factory::HasMethod(const TString &datasetname, const TString &methodTitle) const
{
   if (fMethodsMap.find(datasetname) == fMethodsMap.end())
      return 0;

   std::string methodName = methodTitle.Data();
   auto isEqualToMethodName = [&methodName](TMVA::IMethod *m) { return (0 == methodName.compare(m->GetName())); };

   TMVA::Factory::MVector *methods = this->fMethodsMap.at(datasetname);
   Bool_t isMethodNameExisting = std::any_of(methods->begin(), methods->end(), isEqualToMethodName);

   return isMethodNameExisting;
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::Factory::WriteDataInformation(DataSetInfo &fDataSetInfo)
{
   RootBaseDir()->cd();

   if (!RootBaseDir()->GetDirectory(fDataSetInfo.GetName()))
      RootBaseDir()->mkdir(fDataSetInfo.GetName());
   else
      return; // loader is now in the output file, we dont need to save again

   RootBaseDir()->cd(fDataSetInfo.GetName());
   fDataSetInfo.GetDataSet(); // builds dataset (including calculation of correlation matrix)

   // correlation matrix of the default DS
   const TMatrixD *m(0);
   const TH2 *h(0);

   if (fAnalysisType == Types::kMulticlass) {
      for (UInt_t cls = 0; cls < fDataSetInfo.GetNClasses(); cls++) {
         m = fDataSetInfo.CorrelationMatrix(fDataSetInfo.GetClassInfo(cls)->GetName());
         h = fDataSetInfo.CreateCorrelationMatrixHist(
            m, TString("CorrelationMatrix") + fDataSetInfo.GetClassInfo(cls)->GetName(),
            TString("Correlation Matrix (") + fDataSetInfo.GetClassInfo(cls)->GetName() + TString(")"));
         if (h != 0) {
            h->Write();
            delete h;
         }
      }
   } else {
      m = fDataSetInfo.CorrelationMatrix("Signal");
      h = fDataSetInfo.CreateCorrelationMatrixHist(m, "CorrelationMatrixS", "Correlation Matrix (signal)");
      if (h != 0) {
         h->Write();
         delete h;
      }

      m = fDataSetInfo.CorrelationMatrix("Background");
      h = fDataSetInfo.CreateCorrelationMatrixHist(m, "CorrelationMatrixB", "Correlation Matrix (background)");
      if (h != 0) {
         h->Write();
         delete h;
      }

      m = fDataSetInfo.CorrelationMatrix("Regression");
      h = fDataSetInfo.CreateCorrelationMatrixHist(m, "CorrelationMatrix", "Correlation Matrix");
      if (h != 0) {
         h->Write();
         delete h;
      }
   }

   // some default transformations to evaluate
   // NOTE: all transformations are destroyed after this test
   TString processTrfs = "I"; //"I;N;D;P;U;G,D;"

   // plus some user defined transformations
   processTrfs = fTransformations;

   // remove any trace of identity transform - if given (avoid to apply it twice)
   std::vector<TMVA::TransformationHandler *> trfs;
   TransformationHandler *identityTrHandler = 0;

   std::vector<TString> trfsDef = gTools().SplitString(processTrfs, ';');
   std::vector<TString>::iterator trfsDefIt = trfsDef.begin();
   for (; trfsDefIt != trfsDef.end(); ++trfsDefIt) {
      trfs.push_back(new TMVA::TransformationHandler(fDataSetInfo, "Factory"));
      TString trfS = (*trfsDefIt);

      // Log() << kINFO << Endl;
      Log() << kDEBUG << "current transformation string: '" << trfS.Data() << "'" << Endl;
      TMVA::CreateVariableTransforms(trfS, fDataSetInfo, *(trfs.back()), Log());

      if (trfS.BeginsWith('I'))
         identityTrHandler = trfs.back();
   }

   const std::vector<Event *> &inputEvents = fDataSetInfo.GetDataSet()->GetEventCollection();

   // apply all transformations
   std::vector<TMVA::TransformationHandler *>::iterator trfIt = trfs.begin();

   for (; trfIt != trfs.end(); ++trfIt) {
      // setting a Root dir causes the variables distributions to be saved to the root file
      (*trfIt)->SetRootDir(RootBaseDir()->GetDirectory(fDataSetInfo.GetName())); // every dataloader have its own dir
      (*trfIt)->CalcTransformations(inputEvents);
   }
   if (identityTrHandler)
      identityTrHandler->PrintVariableRanking();

   // clean up
   for (trfIt = trfs.begin(); trfIt != trfs.end(); ++trfIt)
      delete *trfIt;
}

////////////////////////////////////////////////////////////////////////////////
/// Iterates through all booked methods and sees if they use parameter tuning and if so
/// does just that, i.e.\ calls "Method::Train()" for different parameter settings and
/// keeps in mind the "optimal one"...\ and that's the one that will later on be used
/// in the main training loop.

std::map<TString, Double_t> TMVA::Factory::OptimizeAllMethods(TString fomType, TString fitType)
{

   std::map<TString, MVector *>::iterator itrMap;
   std::map<TString, Double_t> TunedParameters;
   for (itrMap = fMethodsMap.begin(); itrMap != fMethodsMap.end(); ++itrMap) {
      MVector *methods = itrMap->second;

      MVector::iterator itrMethod;

      // iterate over methods and optimize
      for (itrMethod = methods->begin(); itrMethod != methods->end(); ++itrMethod) {
         Event::SetIsTraining(kTRUE);
         MethodBase *mva = dynamic_cast<MethodBase *>(*itrMethod);
         if (!mva) {
            Log() << kFATAL << "Dynamic cast to MethodBase failed" << Endl;
            return TunedParameters;
         }

         if (mva->Data()->GetNTrainingEvents() < MinNoTrainingEvents) {
            Log() << kWARNING << "Method " << mva->GetMethodName() << " not trained (training tree has less entries ["
                  << mva->Data()->GetNTrainingEvents() << "] than required [" << MinNoTrainingEvents << "]" << Endl;
            continue;
         }

         Log() << kINFO << "Optimize method: " << mva->GetMethodName() << " for "
               << (fAnalysisType == Types::kRegression
                      ? "Regression"
                      : (fAnalysisType == Types::kMulticlass ? "Multiclass classification" : "Classification"))
               << Endl;

         TunedParameters = mva->OptimizeTuningParameters(fomType, fitType);
         Log() << kINFO << "Optimization of tuning parameters finished for Method:" << mva->GetName() << Endl;
      }
   }

   return TunedParameters;
}

////////////////////////////////////////////////////////////////////////////////
/// Private method to generate a ROCCurve instance for a given method.
/// Handles the conversion from TMVA ResultSet to a format the ROCCurve class
/// understands.
///
/// \note You own the retured pointer.
///

TMVA::ROCCurve *
TMVA::Factory::GetROC(TMVA::DataLoader *loader, TString theMethodName, UInt_t iClass, Types::ETreeType type)
{
   return GetROC((TString)loader->GetName(), theMethodName, iClass, type);
}

////////////////////////////////////////////////////////////////////////////////
/// Private method to generate a ROCCurve instance for a given method.
/// Handles the conversion from TMVA ResultSet to a format the ROCCurve class
/// understands.
///
/// \note You own the retured pointer.
///

TMVA::ROCCurve *TMVA::Factory::GetROC(TString datasetname, TString theMethodName, UInt_t iClass, Types::ETreeType type)
{
   if (fMethodsMap.find(datasetname) == fMethodsMap.end()) {
      Log() << kERROR << Form("DataSet = %s not found in methods map.", datasetname.Data()) << Endl;
      return nullptr;
   }

   if (!this->HasMethod(datasetname, theMethodName)) {
      Log() << kERROR << Form("Method = %s not found with Dataset = %s ", theMethodName.Data(), datasetname.Data())
            << Endl;
      return nullptr;
   }

   std::set<Types::EAnalysisType> allowedAnalysisTypes = {Types::kClassification, Types::kMulticlass};
   if (allowedAnalysisTypes.count(this->fAnalysisType) == 0) {
      Log() << kERROR << Form("Can only generate ROC curves for analysis type kClassification and kMulticlass.")
            << Endl;
      return nullptr;
   }

   TMVA::MethodBase *method = dynamic_cast<TMVA::MethodBase *>(this->GetMethod(datasetname, theMethodName));
   TMVA::DataSet *dataset = method->Data();
   dataset->SetCurrentType(type);
   TMVA::Results *results = dataset->GetResults(theMethodName, type, this->fAnalysisType);

   UInt_t nClasses = method->DataInfo().GetNClasses();
   if (this->fAnalysisType == Types::kMulticlass && iClass >= nClasses) {
      Log() << kERROR
            << Form("Given class number (iClass = %i) does not exist. There are %i classes in dataset.", iClass,
                    nClasses)
            << Endl;
      return nullptr;
   }

   TMVA::ROCCurve *rocCurve = nullptr;
   if (this->fAnalysisType == Types::kClassification) {

      std::vector<Float_t> *mvaRes = dynamic_cast<ResultsClassification *>(results)->GetValueVector();
      std::vector<Bool_t> *mvaResTypes = dynamic_cast<ResultsClassification *>(results)->GetValueVectorTypes();
      std::vector<Float_t> mvaResWeights;

      auto eventCollection = dataset->GetEventCollection(type);
      mvaResWeights.reserve(eventCollection.size());
      for (auto ev : eventCollection) {
         mvaResWeights.push_back(ev->GetWeight());
      }

      rocCurve = new TMVA::ROCCurve(*mvaRes, *mvaResTypes, mvaResWeights);

   } else if (this->fAnalysisType == Types::kMulticlass) {
      std::vector<Float_t> mvaRes;
      std::vector<Bool_t> mvaResTypes;
      std::vector<Float_t> mvaResWeights;

      std::vector<std::vector<Float_t>> *rawMvaRes = dynamic_cast<ResultsMulticlass *>(results)->GetValueVector();

      // Vector transpose due to values being stored as
      //    [ [0, 1, 2], [0, 1, 2], ... ]
      // in ResultsMulticlass::GetValueVector.
      mvaRes.reserve(rawMvaRes->size());
      for (auto item : *rawMvaRes) {
         mvaRes.push_back(item[iClass]);
      }

      auto eventCollection = dataset->GetEventCollection(type);
      mvaResTypes.reserve(eventCollection.size());
      mvaResWeights.reserve(eventCollection.size());
      for (auto ev : eventCollection) {
         mvaResTypes.push_back(ev->GetClass() == iClass);
         mvaResWeights.push_back(ev->GetWeight());
      }

      rocCurve = new TMVA::ROCCurve(mvaRes, mvaResTypes, mvaResWeights);
   }

   return rocCurve;
}

////////////////////////////////////////////////////////////////////////////////
/// Calculate the integral of the ROC curve, also known as the area under curve
/// (AUC), for a given method.
///
/// Argument iClass specifies the class to generate the ROC curve in a
/// multiclass setting. It is ignored for binary classification.
///

Double_t
TMVA::Factory::GetROCIntegral(TMVA::DataLoader *loader, TString theMethodName, UInt_t iClass, Types::ETreeType type)
{
   return GetROCIntegral((TString)loader->GetName(), theMethodName, iClass, type);
}

////////////////////////////////////////////////////////////////////////////////
/// Calculate the integral of the ROC curve, also known as the area under curve
/// (AUC), for a given method.
///
/// Argument iClass specifies the class to generate the ROC curve in a
/// multiclass setting. It is ignored for binary classification.
///

Double_t TMVA::Factory::GetROCIntegral(TString datasetname, TString theMethodName, UInt_t iClass, Types::ETreeType type)
{
   if (fMethodsMap.find(datasetname) == fMethodsMap.end()) {
      Log() << kERROR << Form("DataSet = %s not found in methods map.", datasetname.Data()) << Endl;
      return 0;
   }

   if (!this->HasMethod(datasetname, theMethodName)) {
      Log() << kERROR << Form("Method = %s not found with Dataset = %s ", theMethodName.Data(), datasetname.Data())
            << Endl;
      return 0;
   }

   std::set<Types::EAnalysisType> allowedAnalysisTypes = {Types::kClassification, Types::kMulticlass};
   if (allowedAnalysisTypes.count(this->fAnalysisType) == 0) {
      Log() << kERROR << Form("Can only generate ROC integral for analysis type kClassification. and kMulticlass.")
            << Endl;
      return 0;
   }

   TMVA::ROCCurve *rocCurve = GetROC(datasetname, theMethodName, iClass, type);
   if (!rocCurve) {
      Log() << kFATAL
            << Form("ROCCurve object was not created in Method = %s not found with Dataset = %s ", theMethodName.Data(),
                    datasetname.Data())
            << Endl;
      return 0;
   }

   Int_t npoints = TMVA::gConfig().fVariablePlotting.fNbinsXOfROCCurve + 1;
   Double_t rocIntegral = rocCurve->GetROCIntegral(npoints);
   delete rocCurve;

   return rocIntegral;
}

////////////////////////////////////////////////////////////////////////////////
/// Argument iClass specifies the class to generate the ROC curve in a
/// multiclass setting. It is ignored for binary classification.
///
/// Returns a ROC graph for a given method, or nullptr on error.
///
/// Note: Evaluation of the given method must have been run prior to ROC
/// generation through Factory::EvaluateAllMetods.
///
/// NOTE: The ROC curve is 1 vs. all where the given class is considered signal
/// and the others considered background. This is ok in binary classification
/// but in in multi class classification, the ROC surface is an N dimensional
/// shape, where N is number of classes - 1.

TGraph *TMVA::Factory::GetROCCurve(DataLoader *loader, TString theMethodName, Bool_t setTitles, UInt_t iClass,
                                   Types::ETreeType type)
{
   return GetROCCurve((TString)loader->GetName(), theMethodName, setTitles, iClass, type);
}

////////////////////////////////////////////////////////////////////////////////
/// Argument iClass specifies the class to generate the ROC curve in a
/// multiclass setting. It is ignored for binary classification.
///
/// Returns a ROC graph for a given method, or nullptr on error.
///
/// Note: Evaluation of the given method must have been run prior to ROC
/// generation through Factory::EvaluateAllMetods.
///
/// NOTE: The ROC curve is 1 vs. all where the given class is considered signal
/// and the others considered background. This is ok in binary classification
/// but in in multi class classification, the ROC surface is an N dimensional
/// shape, where N is number of classes - 1.

TGraph *TMVA::Factory::GetROCCurve(TString datasetname, TString theMethodName, Bool_t setTitles, UInt_t iClass,
                                   Types::ETreeType type)
{
   if (fMethodsMap.find(datasetname) == fMethodsMap.end()) {
      Log() << kERROR << Form("DataSet = %s not found in methods map.", datasetname.Data()) << Endl;
      return nullptr;
   }

   if (!this->HasMethod(datasetname, theMethodName)) {
      Log() << kERROR << Form("Method = %s not found with Dataset = %s ", theMethodName.Data(), datasetname.Data())
            << Endl;
      return nullptr;
   }

   std::set<Types::EAnalysisType> allowedAnalysisTypes = {Types::kClassification, Types::kMulticlass};
   if (allowedAnalysisTypes.count(this->fAnalysisType) == 0) {
      Log() << kERROR << Form("Can only generate ROC curves for analysis type kClassification and kMulticlass.")
            << Endl;
      return nullptr;
   }

   TMVA::ROCCurve *rocCurve = GetROC(datasetname, theMethodName, iClass, type);
   TGraph *graph = nullptr;

   if (!rocCurve) {
      Log() << kFATAL
            << Form("ROCCurve object was not created in Method = %s not found with Dataset = %s ", theMethodName.Data(),
                    datasetname.Data())
            << Endl;
      return nullptr;
   }

   graph = (TGraph *)rocCurve->GetROCCurve()->Clone();
   delete rocCurve;

   if (setTitles) {
      graph->GetYaxis()->SetTitle("Background rejection (Specificity)");
      graph->GetXaxis()->SetTitle("Signal efficiency (Sensitivity)");
      graph->SetTitle(Form("Signal efficiency vs. Background rejection (%s)", theMethodName.Data()));
   }

   return graph;
}

////////////////////////////////////////////////////////////////////////////////
/// Generate a collection of graphs, for all methods for a given class. Suitable
/// for comparing method performance.
///
/// Argument iClass specifies the class to generate the ROC curve in a
/// multiclass setting. It is ignored for binary classification.
///
/// NOTE: The ROC curve is 1 vs. all where the given class is considered signal
/// and the others considered background. This is ok in binary classification
/// but in in multi class classification, the ROC surface is an N dimensional
/// shape, where N is number of classes - 1.

TMultiGraph *TMVA::Factory::GetROCCurveAsMultiGraph(DataLoader *loader, UInt_t iClass, Types::ETreeType type)
{
   return GetROCCurveAsMultiGraph((TString)loader->GetName(), iClass, type);
}

////////////////////////////////////////////////////////////////////////////////
/// Generate a collection of graphs, for all methods for a given class. Suitable
/// for comparing method performance.
///
/// Argument iClass specifies the class to generate the ROC curve in a
/// multiclass setting. It is ignored for binary classification.
///
/// NOTE: The ROC curve is 1 vs. all where the given class is considered signal
/// and the others considered background. This is ok in binary classification
/// but in in multi class classification, the ROC surface is an N dimensional
/// shape, where N is number of classes - 1.

TMultiGraph *TMVA::Factory::GetROCCurveAsMultiGraph(TString datasetname, UInt_t iClass, Types::ETreeType type)
{
   UInt_t line_color = 1;

   TMultiGraph *multigraph = new TMultiGraph();

   MVector *methods = fMethodsMap[datasetname.Data()];
   for (auto *method_raw : *methods) {
      TMVA::MethodBase *method = dynamic_cast<TMVA::MethodBase *>(method_raw);
      if (method == nullptr) {
         continue;
      }

      TString methodName = method->GetMethodName();
      UInt_t nClasses = method->DataInfo().GetNClasses();

      if (this->fAnalysisType == Types::kMulticlass && iClass >= nClasses) {
         Log() << kERROR
               << Form("Given class number (iClass = %i) does not exist. There are %i classes in dataset.", iClass,
                       nClasses)
               << Endl;
         continue;
      }

      TString className = method->DataInfo().GetClassInfo(iClass)->GetName();

      TGraph *graph = this->GetROCCurve(datasetname, methodName, false, iClass, type);
      graph->SetTitle(methodName);

      graph->SetLineWidth(2);
      graph->SetLineColor(line_color++);
      graph->SetFillColor(10);

      multigraph->Add(graph);
   }

   if (multigraph->GetListOfGraphs() == nullptr) {
      Log() << kERROR << Form("No metohds have class %i defined.", iClass) << Endl;
      return nullptr;
   }

   return multigraph;
}

////////////////////////////////////////////////////////////////////////////////
/// Draws ROC curves for all methods booked with the factory for a given class
/// onto a canvas.
///
/// Argument iClass specifies the class to generate the ROC curve in a
/// multiclass setting. It is ignored for binary classification.
///
/// NOTE: The ROC curve is 1 vs. all where the given class is considered signal
/// and the others considered background. This is ok in binary classification
/// but in in multi class classification, the ROC surface is an N dimensional
/// shape, where N is number of classes - 1.

TCanvas *TMVA::Factory::GetROCCurve(TMVA::DataLoader *loader, UInt_t iClass, Types::ETreeType type)
{
   return GetROCCurve((TString)loader->GetName(), iClass, type);
}

////////////////////////////////////////////////////////////////////////////////
/// Draws ROC curves for all methods booked with the factory for a given class.
///
/// Argument iClass specifies the class to generate the ROC curve in a
/// multiclass setting. It is ignored for binary classification.
///
/// NOTE: The ROC curve is 1 vs. all where the given class is considered signal
/// and the others considered background. This is ok in binary classification
/// but in in multi class classification, the ROC surface is an N dimensional
/// shape, where N is number of classes - 1.

TCanvas *TMVA::Factory::GetROCCurve(TString datasetname, UInt_t iClass, Types::ETreeType type)
{
   if (fMethodsMap.find(datasetname) == fMethodsMap.end()) {
      Log() << kERROR << Form("DataSet = %s not found in methods map.", datasetname.Data()) << Endl;
      return 0;
   }

   TString name = Form("ROCCurve %s class %i", datasetname.Data(), iClass);
   TCanvas *canvas = new TCanvas(name, "ROC Curve", 200, 10, 700, 500);
   canvas->SetGrid();

   TMultiGraph *multigraph = this->GetROCCurveAsMultiGraph(datasetname, iClass, type);

   if (multigraph) {
      multigraph->Draw("AL");

      multigraph->GetYaxis()->SetTitle("Background rejection (Specificity)");
      multigraph->GetXaxis()->SetTitle("Signal efficiency (Sensitivity)");

      TString titleString = Form("Signal efficiency vs. Background rejection");
      if (this->fAnalysisType == Types::kMulticlass) {
         titleString = Form("%s (Class=%i)", titleString.Data(), iClass);
      }

      // Workaround for TMultigraph not drawing title correctly.
      multigraph->GetHistogram()->SetTitle(titleString);
      multigraph->SetTitle(titleString);

      canvas->BuildLegend(0.15, 0.15, 0.35, 0.3, "MVA Method");
   }

   return canvas;
}

////////////////////////////////////////////////////////////////////////////////
/// Iterates through all booked methods and calls training

void TMVA::Factory::TrainAllMethods()
{
   Log() << kHEADER << gTools().Color("bold") << "Train all methods" << gTools().Color("reset") << Endl;
   // iterates over all MVAs that have been booked, and calls their training methods

   // don't do anything if no method booked
   if (fMethodsMap.empty()) {
      Log() << kINFO << "...nothing found to train" << Endl;
      return;
   }

   // here the training starts
   // Log() << kINFO << " " << Endl;
   Log() << kDEBUG << "Train all methods for "
         << (fAnalysisType == Types::kRegression
                ? "Regression"
                : (fAnalysisType == Types::kMulticlass ? "Multiclass" : "Classification"))
         << " ..." << Endl;

   std::map<TString, MVector *>::iterator itrMap;

   for (itrMap = fMethodsMap.begin(); itrMap != fMethodsMap.end(); ++itrMap) {
      MVector *methods = itrMap->second;
      MVector::iterator itrMethod;

      // iterate over methods and train
      for (itrMethod = methods->begin(); itrMethod != methods->end(); ++itrMethod) {
         Event::SetIsTraining(kTRUE);
         MethodBase *mva = dynamic_cast<MethodBase *>(*itrMethod);

         if (mva == 0)
            continue;

         if (mva->DataInfo().GetDataSetManager()->DataInput().GetEntries() <=
             1) { // 0 entries --> 0 events, 1 entry --> dynamical dataset (or one entry)
            Log() << kFATAL << "No input data for the training provided!" << Endl;
         }

         if (fAnalysisType == Types::kRegression && mva->DataInfo().GetNTargets() < 1)
            Log() << kFATAL << "You want to do regression training without specifying a target." << Endl;
         else if ((fAnalysisType == Types::kMulticlass || fAnalysisType == Types::kClassification) &&
                  mva->DataInfo().GetNClasses() < 2)
            Log() << kFATAL << "You want to do classification training, but specified less than two classes." << Endl;

         // first print some information about the default dataset
         if (!IsSilentFile())
            WriteDataInformation(mva->fDataSetInfo);

         if (mva->Data()->GetNTrainingEvents() < MinNoTrainingEvents) {
            Log() << kWARNING << "Method " << mva->GetMethodName() << " not trained (training tree has less entries ["
                  << mva->Data()->GetNTrainingEvents() << "] than required [" << MinNoTrainingEvents << "]" << Endl;
            continue;
         }

         Log() << kHEADER << "Train method: " << mva->GetMethodName() << " for "
               << (fAnalysisType == Types::kRegression
                      ? "Regression"
                      : (fAnalysisType == Types::kMulticlass ? "Multiclass classification" : "Classification"))
               << Endl << Endl;
         mva->TrainMethod();
         Log() << kHEADER << "Training finished" << Endl << Endl;
      }

      if (fAnalysisType != Types::kRegression) {

         // variable ranking
         // Log() << Endl;
         Log() << kINFO << "Ranking input variables (method specific)..." << Endl;
         for (itrMethod = methods->begin(); itrMethod != methods->end(); ++itrMethod) {
            MethodBase *mva = dynamic_cast<MethodBase *>(*itrMethod);
            if (mva && mva->Data()->GetNTrainingEvents() >= MinNoTrainingEvents) {

               // create and print ranking
               const Ranking *ranking = (*itrMethod)->CreateRanking();
               if (ranking != 0)
                  ranking->Print();
               else
                  Log() << kINFO << "No variable ranking supplied by classifier: "
                        << dynamic_cast<MethodBase *>(*itrMethod)->GetMethodName() << Endl;
            }
         }
      }

      // save training history in case we are not in the silent mode
      if (!IsSilentFile()) {
         for (UInt_t i = 0; i < methods->size(); i++) {
            MethodBase *m = dynamic_cast<MethodBase *>((*methods)[i]);
            if (m == 0)
               continue;
            m->BaseDir()->cd();
            m->fTrainHistory.SaveHistory(m->GetMethodName());
         }
      }

      // delete all methods and recreate them from weight file - this ensures that the application
      // of the methods (in TMVAClassificationApplication) is consistent with the results obtained
      // in the testing
      // Log() << Endl;
      if (fModelPersistence) {

         Log() << kHEADER << "=== Destroy and recreate all methods via weight files for testing ===" << Endl << Endl;

         if (!IsSilentFile())
            RootBaseDir()->cd();

         // iterate through all booked methods
         for (UInt_t i = 0; i < methods->size(); i++) {

            MethodBase *m = dynamic_cast<MethodBase *>((*methods)[i]);
            if (m == nullptr)
               continue;

            TMVA::Types::EMVA methodType = m->GetMethodType();
            TString weightfile = m->GetWeightFileName();

            // decide if .txt or .xml file should be read:
            if (READXML)
               weightfile.ReplaceAll(".txt", ".xml");

            DataSetInfo &dataSetInfo = m->DataInfo();
            TString testvarName = m->GetTestvarName();
            delete m; // itrMethod[i];

            // recreate
            m = dynamic_cast<MethodBase *>(ClassifierFactory::Instance().Create(
               Types::Instance().GetMethodName(methodType).Data(), dataSetInfo, weightfile));
            if (m->GetMethodType() == Types::kCategory) {
               MethodCategory *methCat = (dynamic_cast<MethodCategory *>(m));
               if (!methCat)
                  Log() << kFATAL << "Method with type kCategory cannot be casted to MethodCategory. /Factory" << Endl;
               else
                  methCat->fDataSetManager = m->DataInfo().GetDataSetManager();
            }
            // ToDo, Do we need to fill the DataSetManager of MethodBoost here too?

            TString wfileDir = m->DataInfo().GetName();
            wfileDir += "/" + gConfig().GetIONames().fWeightFileDir;
            m->SetWeightFileDir(wfileDir);
            m->SetModelPersistence(fModelPersistence);
            m->SetSilentFile(IsSilentFile());
            m->SetAnalysisType(fAnalysisType);
            m->SetupMethod();
            m->ReadStateFromFile();
            m->SetTestvarName(testvarName);

            // replace trained method by newly created one (from weight file) in methods vector
            (*methods)[i] = m;
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Evaluates all booked methods on the testing data and adds the output to the
/// Results in the corresponiding DataSet.
///

void TMVA::Factory::TestAllMethods()
{
   Log() << kHEADER << gTools().Color("bold") << "Test all methods" << gTools().Color("reset") << Endl;

   // don't do anything if no method booked
   if (fMethodsMap.empty()) {
      Log() << kINFO << "...nothing found to test" << Endl;
      return;
   }
   std::map<TString, MVector *>::iterator itrMap;

   for (itrMap = fMethodsMap.begin(); itrMap != fMethodsMap.end(); ++itrMap) {
      MVector *methods = itrMap->second;
      MVector::iterator itrMethod;

      // iterate over methods and test
      for (itrMethod = methods->begin(); itrMethod != methods->end(); ++itrMethod) {
         Event::SetIsTraining(kFALSE);
         MethodBase *mva = dynamic_cast<MethodBase *>(*itrMethod);
         if (mva == 0)
            continue;
         Types::EAnalysisType analysisType = mva->GetAnalysisType();
         Log() << kHEADER << "Test method: " << mva->GetMethodName() << " for "
               << (analysisType == Types::kRegression
                      ? "Regression"
                      : (analysisType == Types::kMulticlass ? "Multiclass classification" : "Classification"))
               << " performance" << Endl << Endl;
         mva->AddOutput(Types::kTesting, analysisType);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::Factory::MakeClass(const TString &datasetname, const TString &methodTitle) const
{
   if (methodTitle != "") {
      IMethod *method = GetMethod(datasetname, methodTitle);
      if (method)
         method->MakeClass();
      else {
         Log() << kWARNING << "<MakeClass> Could not find classifier \"" << methodTitle << "\" in list" << Endl;
      }
   } else {

      // no classifier specified, print all help messages
      MVector *methods = fMethodsMap.find(datasetname)->second;
      MVector::const_iterator itrMethod;
      for (itrMethod = methods->begin(); itrMethod != methods->end(); ++itrMethod) {
         MethodBase *method = dynamic_cast<MethodBase *>(*itrMethod);
         if (method == 0)
            continue;
         Log() << kINFO << "Make response class for classifier: " << method->GetMethodName() << Endl;
         method->MakeClass();
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Print predefined help message of classifier.
/// Iterate over methods and test.

void TMVA::Factory::PrintHelpMessage(const TString &datasetname, const TString &methodTitle) const
{
   if (methodTitle != "") {
      IMethod *method = GetMethod(datasetname, methodTitle);
      if (method)
         method->PrintHelpMessage();
      else {
         Log() << kWARNING << "<PrintHelpMessage> Could not find classifier \"" << methodTitle << "\" in list" << Endl;
      }
   } else {

      // no classifier specified, print all help messages
      MVector *methods = fMethodsMap.find(datasetname)->second;
      MVector::const_iterator itrMethod;
      for (itrMethod = methods->begin(); itrMethod != methods->end(); ++itrMethod) {
         MethodBase *method = dynamic_cast<MethodBase *>(*itrMethod);
         if (method == 0)
            continue;
         Log() << kINFO << "Print help message for classifier: " << method->GetMethodName() << Endl;
         method->PrintHelpMessage();
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Iterates over all MVA input variables and evaluates them.

void TMVA::Factory::EvaluateAllVariables(DataLoader *loader, TString options)
{
   Log() << kINFO << "Evaluating all variables..." << Endl;
   Event::SetIsTraining(kFALSE);

   for (UInt_t i = 0; i < loader->GetDataSetInfo().GetNVariables(); i++) {
      TString s = loader->GetDataSetInfo().GetVariableInfo(i).GetLabel();
      if (options.Contains("V"))
         s += ":V";
      this->BookMethod(loader, "Variable", s);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Iterates over all MVAs that have been booked, and calls their evaluation methods.

void TMVA::Factory::EvaluateAllMethods(void)
{
   Log() << kHEADER << gTools().Color("bold") << "Evaluate all methods" << gTools().Color("reset") << Endl;

   // don't do anything if no method booked
   if (fMethodsMap.empty()) {
      Log() << kINFO << "...nothing found to evaluate" << Endl;
      return;
   }
   std::map<TString, MVector *>::iterator itrMap;

   for (itrMap = fMethodsMap.begin(); itrMap != fMethodsMap.end(); ++itrMap) {
      MVector *methods = itrMap->second;

      // -----------------------------------------------------------------------
      // First part of evaluation process
      // --> compute efficiencies, and other separation estimators
      // -----------------------------------------------------------------------

      // although equal, we now want to separate the output for the variables
      // and the real methods
      Int_t isel;                   // will be 0 for a Method; 1 for a Variable
      Int_t nmeth_used[2] = {0, 0}; // 0 Method; 1 Variable

      std::vector<std::vector<TString>> mname(2);
      std::vector<std::vector<Double_t>> sig(2), sep(2), roc(2);
      std::vector<std::vector<Double_t>> eff01(2), eff10(2), eff30(2), effArea(2);
      std::vector<std::vector<Double_t>> eff01err(2), eff10err(2), eff30err(2);
      std::vector<std::vector<Double_t>> trainEff01(2), trainEff10(2), trainEff30(2);

      std::vector<std::vector<Float_t>> multiclass_testEff;
      std::vector<std::vector<Float_t>> multiclass_trainEff;
      std::vector<std::vector<Float_t>> multiclass_testPur;
      std::vector<std::vector<Float_t>> multiclass_trainPur;

      std::vector<std::vector<Float_t>> train_history;

      // Multiclass confusion matrices.
      std::vector<TMatrixD> multiclass_trainConfusionEffB01;
      std::vector<TMatrixD> multiclass_trainConfusionEffB10;
      std::vector<TMatrixD> multiclass_trainConfusionEffB30;
      std::vector<TMatrixD> multiclass_testConfusionEffB01;
      std::vector<TMatrixD> multiclass_testConfusionEffB10;
      std::vector<TMatrixD> multiclass_testConfusionEffB30;

      std::vector<std::vector<Double_t>> biastrain(1); // "bias" of the regression on the training data
      std::vector<std::vector<Double_t>> biastest(1);  // "bias" of the regression on test data
      std::vector<std::vector<Double_t>> devtrain(1);  // "dev" of the regression on the training data
      std::vector<std::vector<Double_t>> devtest(1);   // "dev" of the regression on test data
      std::vector<std::vector<Double_t>> rmstrain(1);  // "rms" of the regression on the training data
      std::vector<std::vector<Double_t>> rmstest(1);   // "rms" of the regression on test data
      std::vector<std::vector<Double_t>> minftrain(1); // "minf" of the regression on the training data
      std::vector<std::vector<Double_t>> minftest(1);  // "minf" of the regression on test data
      std::vector<std::vector<Double_t>> rhotrain(1);  // correlation of the regression on the training data
      std::vector<std::vector<Double_t>> rhotest(1);   // correlation of the regression on test data

      // same as above but for 'truncated' quantities (computed for events within 2sigma of RMS)
      std::vector<std::vector<Double_t>> biastrainT(1);
      std::vector<std::vector<Double_t>> biastestT(1);
      std::vector<std::vector<Double_t>> devtrainT(1);
      std::vector<std::vector<Double_t>> devtestT(1);
      std::vector<std::vector<Double_t>> rmstrainT(1);
      std::vector<std::vector<Double_t>> rmstestT(1);
      std::vector<std::vector<Double_t>> minftrainT(1);
      std::vector<std::vector<Double_t>> minftestT(1);

      // following vector contains all methods - with the exception of Cuts, which are special
      MVector methodsNoCuts;

      Bool_t doRegression = kFALSE;
      Bool_t doMulticlass = kFALSE;

      // iterate over methods and evaluate
      for (MVector::iterator itrMethod = methods->begin(); itrMethod != methods->end(); ++itrMethod) {
         Event::SetIsTraining(kFALSE);
         MethodBase *theMethod = dynamic_cast<MethodBase *>(*itrMethod);
         if (theMethod == 0)
            continue;
         theMethod->SetFile(fgTargetFile);
         theMethod->SetSilentFile(IsSilentFile());
         if (theMethod->GetMethodType() != Types::kCuts)
            methodsNoCuts.push_back(*itrMethod);

         if (theMethod->DoRegression()) {
            doRegression = kTRUE;

            Log() << kINFO << "Evaluate regression method: " << theMethod->GetMethodName() << Endl;
            Double_t bias, dev, rms, mInf;
            Double_t biasT, devT, rmsT, mInfT;
            Double_t rho;

            Log() << kINFO << "TestRegression (testing)" << Endl;
            theMethod->TestRegression(bias, biasT, dev, devT, rms, rmsT, mInf, mInfT, rho, TMVA::Types::kTesting);
            biastest[0].push_back(bias);
            devtest[0].push_back(dev);
            rmstest[0].push_back(rms);
            minftest[0].push_back(mInf);
            rhotest[0].push_back(rho);
            biastestT[0].push_back(biasT);
            devtestT[0].push_back(devT);
            rmstestT[0].push_back(rmsT);
            minftestT[0].push_back(mInfT);

            Log() << kINFO << "TestRegression (training)" << Endl;
            theMethod->TestRegression(bias, biasT, dev, devT, rms, rmsT, mInf, mInfT, rho, TMVA::Types::kTraining);
            biastrain[0].push_back(bias);
            devtrain[0].push_back(dev);
            rmstrain[0].push_back(rms);
            minftrain[0].push_back(mInf);
            rhotrain[0].push_back(rho);
            biastrainT[0].push_back(biasT);
            devtrainT[0].push_back(devT);
            rmstrainT[0].push_back(rmsT);
            minftrainT[0].push_back(mInfT);

            mname[0].push_back(theMethod->GetMethodName());
            nmeth_used[0]++;
            if (!IsSilentFile()) {
               Log() << kDEBUG << "\tWrite evaluation histograms to file" << Endl;
               theMethod->WriteEvaluationHistosToFile(Types::kTesting);
               theMethod->WriteEvaluationHistosToFile(Types::kTraining);
            }
         } else if (theMethod->DoMulticlass()) {
            // ====================================================================
            // === Multiclass evaluation
            // ====================================================================
            doMulticlass = kTRUE;
            Log() << kINFO << "Evaluate multiclass classification method: " << theMethod->GetMethodName() << Endl;

            // This part uses a genetic alg. to evaluate the optimal sig eff * sig pur.
            // This is why it is disabled for now.
            // Find approximate optimal working point w.r.t. signalEfficiency * signalPurity.
            // theMethod->TestMulticlass(); // This is where the actual GA calc is done
            // multiclass_testEff.push_back(theMethod->GetMulticlassEfficiency(multiclass_testPur));

            theMethod->TestMulticlass();

            // Confusion matrix at three background efficiency levels
            multiclass_trainConfusionEffB01.push_back(theMethod->GetMulticlassConfusionMatrix(0.01, Types::kTraining));
            multiclass_trainConfusionEffB10.push_back(theMethod->GetMulticlassConfusionMatrix(0.10, Types::kTraining));
            multiclass_trainConfusionEffB30.push_back(theMethod->GetMulticlassConfusionMatrix(0.30, Types::kTraining));

            multiclass_testConfusionEffB01.push_back(theMethod->GetMulticlassConfusionMatrix(0.01, Types::kTesting));
            multiclass_testConfusionEffB10.push_back(theMethod->GetMulticlassConfusionMatrix(0.10, Types::kTesting));
            multiclass_testConfusionEffB30.push_back(theMethod->GetMulticlassConfusionMatrix(0.30, Types::kTesting));

            if (!IsSilentFile()) {
               Log() << kDEBUG << "\tWrite evaluation histograms to file" << Endl;
               theMethod->WriteEvaluationHistosToFile(Types::kTesting);
               theMethod->WriteEvaluationHistosToFile(Types::kTraining);
            }

            nmeth_used[0]++;
            mname[0].push_back(theMethod->GetMethodName());
         } else {

            Log() << kHEADER << "Evaluate classifier: " << theMethod->GetMethodName() << Endl << Endl;
            isel = (theMethod->GetMethodTypeName().Contains("Variable")) ? 1 : 0;

            // perform the evaluation
            theMethod->TestClassification();

            // evaluate the classifier
            mname[isel].push_back(theMethod->GetMethodName());
            sig[isel].push_back(theMethod->GetSignificance());
            sep[isel].push_back(theMethod->GetSeparation());
            roc[isel].push_back(theMethod->GetROCIntegral());

            Double_t err;
            eff01[isel].push_back(theMethod->GetEfficiency("Efficiency:0.01", Types::kTesting, err));
            eff01err[isel].push_back(err);
            eff10[isel].push_back(theMethod->GetEfficiency("Efficiency:0.10", Types::kTesting, err));
            eff10err[isel].push_back(err);
            eff30[isel].push_back(theMethod->GetEfficiency("Efficiency:0.30", Types::kTesting, err));
            eff30err[isel].push_back(err);
            effArea[isel].push_back(theMethod->GetEfficiency("", Types::kTesting, err)); // computes the area (average)

            trainEff01[isel].push_back(
               theMethod->GetTrainingEfficiency("Efficiency:0.01")); // the first pass takes longer
            trainEff10[isel].push_back(theMethod->GetTrainingEfficiency("Efficiency:0.10"));
            trainEff30[isel].push_back(theMethod->GetTrainingEfficiency("Efficiency:0.30"));

            nmeth_used[isel]++;

            if (!IsSilentFile()) {
               Log() << kDEBUG << "\tWrite evaluation histograms to file" << Endl;
               theMethod->WriteEvaluationHistosToFile(Types::kTesting);
               theMethod->WriteEvaluationHistosToFile(Types::kTraining);
            }
         }
      }
      if (doRegression) {

         std::vector<TString> vtemps = mname[0];
         std::vector<std::vector<Double_t>> vtmp;
         vtmp.push_back(devtest[0]); // this is the vector that is ranked
         vtmp.push_back(devtrain[0]);
         vtmp.push_back(biastest[0]);
         vtmp.push_back(biastrain[0]);
         vtmp.push_back(rmstest[0]);
         vtmp.push_back(rmstrain[0]);
         vtmp.push_back(minftest[0]);
         vtmp.push_back(minftrain[0]);
         vtmp.push_back(rhotest[0]);
         vtmp.push_back(rhotrain[0]);
         vtmp.push_back(devtestT[0]); // this is the vector that is ranked
         vtmp.push_back(devtrainT[0]);
         vtmp.push_back(biastestT[0]);
         vtmp.push_back(biastrainT[0]);
         vtmp.push_back(rmstestT[0]);
         vtmp.push_back(rmstrainT[0]);
         vtmp.push_back(minftestT[0]);
         vtmp.push_back(minftrainT[0]);
         gTools().UsefulSortAscending(vtmp, &vtemps);
         mname[0] = vtemps;
         devtest[0] = vtmp[0];
         devtrain[0] = vtmp[1];
         biastest[0] = vtmp[2];
         biastrain[0] = vtmp[3];
         rmstest[0] = vtmp[4];
         rmstrain[0] = vtmp[5];
         minftest[0] = vtmp[6];
         minftrain[0] = vtmp[7];
         rhotest[0] = vtmp[8];
         rhotrain[0] = vtmp[9];
         devtestT[0] = vtmp[10];
         devtrainT[0] = vtmp[11];
         biastestT[0] = vtmp[12];
         biastrainT[0] = vtmp[13];
         rmstestT[0] = vtmp[14];
         rmstrainT[0] = vtmp[15];
         minftestT[0] = vtmp[16];
         minftrainT[0] = vtmp[17];
      } else if (doMulticlass) {
         // TODO: fill in something meaningful
         // If there is some ranking of methods to be done it should be done here.
         // However, this is not so easy to define for multiclass so it is left out for now.

      } else {
         // now sort the variables according to the best 'eff at Beff=0.10'
         for (Int_t k = 0; k < 2; k++) {
            std::vector<std::vector<Double_t>> vtemp;
            vtemp.push_back(effArea[k]); // this is the vector that is ranked
            vtemp.push_back(eff10[k]);
            vtemp.push_back(eff01[k]);
            vtemp.push_back(eff30[k]);
            vtemp.push_back(eff10err[k]);
            vtemp.push_back(eff01err[k]);
            vtemp.push_back(eff30err[k]);
            vtemp.push_back(trainEff10[k]);
            vtemp.push_back(trainEff01[k]);
            vtemp.push_back(trainEff30[k]);
            vtemp.push_back(sig[k]);
            vtemp.push_back(sep[k]);
            vtemp.push_back(roc[k]);
            std::vector<TString> vtemps = mname[k];
            gTools().UsefulSortDescending(vtemp, &vtemps);
            effArea[k] = vtemp[0];
            eff10[k] = vtemp[1];
            eff01[k] = vtemp[2];
            eff30[k] = vtemp[3];
            eff10err[k] = vtemp[4];
            eff01err[k] = vtemp[5];
            eff30err[k] = vtemp[6];
            trainEff10[k] = vtemp[7];
            trainEff01[k] = vtemp[8];
            trainEff30[k] = vtemp[9];
            sig[k] = vtemp[10];
            sep[k] = vtemp[11];
            roc[k] = vtemp[12];
            mname[k] = vtemps;
         }
      }

      // -----------------------------------------------------------------------
      // Second part of evaluation process
      // --> compute correlations among MVAs
      // --> compute correlations between input variables and MVA (determines importance)
      // --> count overlaps
      // -----------------------------------------------------------------------
      if (fCorrelations) {
         const Int_t nmeth = methodsNoCuts.size();
         MethodBase *method = dynamic_cast<MethodBase *>(methods[0][0]);
         const Int_t nvar = method->fDataSetInfo.GetNVariables();
         if (!doRegression && !doMulticlass) {

            if (nmeth > 0) {

               //              needed for correlations
               Double_t *dvec = new Double_t[nmeth + nvar];
               std::vector<Double_t> rvec;

               //              for correlations
               TPrincipal *tpSig = new TPrincipal(nmeth + nvar, "");
               TPrincipal *tpBkg = new TPrincipal(nmeth + nvar, "");

               //              set required tree branch references
               Int_t ivar = 0;
               std::vector<TString> *theVars = new std::vector<TString>;
               std::vector<ResultsClassification *> mvaRes;
               for (MVector::iterator itrMethod = methodsNoCuts.begin(); itrMethod != methodsNoCuts.end();
                    ++itrMethod, ++ivar) {
                  MethodBase *m = dynamic_cast<MethodBase *>(*itrMethod);
                  if (m == 0)
                     continue;
                  theVars->push_back(m->GetTestvarName());
                  rvec.push_back(m->GetSignalReferenceCut());
                  theVars->back().ReplaceAll("MVA_", "");
                  mvaRes.push_back(dynamic_cast<ResultsClassification *>(
                     m->Data()->GetResults(m->GetMethodName(), Types::kTesting, Types::kMaxAnalysisType)));
               }

               //              for overlap study
               TMatrixD *overlapS = new TMatrixD(nmeth, nmeth);
               TMatrixD *overlapB = new TMatrixD(nmeth, nmeth);
               (*overlapS) *= 0; // init...
               (*overlapB) *= 0; // init...

               //              loop over test tree
               DataSet *defDs = method->fDataSetInfo.GetDataSet();
               defDs->SetCurrentType(Types::kTesting);
               for (Int_t ievt = 0; ievt < defDs->GetNEvents(); ievt++) {
                  const Event *ev = defDs->GetEvent(ievt);

                  //                 for correlations
                  TMatrixD *theMat = 0;
                  for (Int_t im = 0; im < nmeth; im++) {
                     //                    check for NaN value
                     Double_t retval = (Double_t)(*mvaRes[im])[ievt][0];
                     if (TMath::IsNaN(retval)) {
                        Log() << kWARNING << "Found NaN return value in event: " << ievt << " for method \""
                              << methodsNoCuts[im]->GetName() << "\"" << Endl;
                        dvec[im] = 0;
                     } else
                        dvec[im] = retval;
                  }
                  for (Int_t iv = 0; iv < nvar; iv++)
                     dvec[iv + nmeth] = (Double_t)ev->GetValue(iv);
                  if (method->fDataSetInfo.IsSignal(ev)) {
                     tpSig->AddRow(dvec);
                     theMat = overlapS;
                  } else {
                     tpBkg->AddRow(dvec);
                     theMat = overlapB;
                  }

                  //                 count overlaps
                  for (Int_t im = 0; im < nmeth; im++) {
                     for (Int_t jm = im; jm < nmeth; jm++) {
                        if ((dvec[im] - rvec[im]) * (dvec[jm] - rvec[jm]) > 0) {
                           (*theMat)(im, jm)++;
                           if (im != jm)
                              (*theMat)(jm, im)++;
                        }
                     }
                  }
               }

               //              renormalise overlap matrix
               (*overlapS) *= (1.0 / defDs->GetNEvtSigTest());  // init...
               (*overlapB) *= (1.0 / defDs->GetNEvtBkgdTest()); // init...

               tpSig->MakePrincipals();
               tpBkg->MakePrincipals();

               const TMatrixD *covMatS = tpSig->GetCovarianceMatrix();
               const TMatrixD *covMatB = tpBkg->GetCovarianceMatrix();

               const TMatrixD *corrMatS = gTools().GetCorrelationMatrix(covMatS);
               const TMatrixD *corrMatB = gTools().GetCorrelationMatrix(covMatB);

               //              print correlation matrices
               if (corrMatS != 0 && corrMatB != 0) {

                  //                 extract MVA matrix
                  TMatrixD mvaMatS(nmeth, nmeth);
                  TMatrixD mvaMatB(nmeth, nmeth);
                  for (Int_t im = 0; im < nmeth; im++) {
                     for (Int_t jm = 0; jm < nmeth; jm++) {
                        mvaMatS(im, jm) = (*corrMatS)(im, jm);
                        mvaMatB(im, jm) = (*corrMatB)(im, jm);
                     }
                  }

                  //                 extract variables - to MVA matrix
                  std::vector<TString> theInputVars;
                  TMatrixD varmvaMatS(nvar, nmeth);
                  TMatrixD varmvaMatB(nvar, nmeth);
                  for (Int_t iv = 0; iv < nvar; iv++) {
                     theInputVars.push_back(method->fDataSetInfo.GetVariableInfo(iv).GetLabel());
                     for (Int_t jm = 0; jm < nmeth; jm++) {
                        varmvaMatS(iv, jm) = (*corrMatS)(nmeth + iv, jm);
                        varmvaMatB(iv, jm) = (*corrMatB)(nmeth + iv, jm);
                     }
                  }

                  if (nmeth > 1) {
                     Log() << kINFO << Endl;
                     Log() << kINFO << Form("Dataset[%s] : ", method->fDataSetInfo.GetName())
                           << "Inter-MVA correlation matrix (signal):" << Endl;
                     gTools().FormattedOutput(mvaMatS, *theVars, Log());
                     Log() << kINFO << Endl;

                     Log() << kINFO << Form("Dataset[%s] : ", method->fDataSetInfo.GetName())
                           << "Inter-MVA correlation matrix (background):" << Endl;
                     gTools().FormattedOutput(mvaMatB, *theVars, Log());
                     Log() << kINFO << Endl;
                  }

                  Log() << kINFO << Form("Dataset[%s] : ", method->fDataSetInfo.GetName())
                        << "Correlations between input variables and MVA response (signal):" << Endl;
                  gTools().FormattedOutput(varmvaMatS, theInputVars, *theVars, Log());
                  Log() << kINFO << Endl;

                  Log() << kINFO << Form("Dataset[%s] : ", method->fDataSetInfo.GetName())
                        << "Correlations between input variables and MVA response (background):" << Endl;
                  gTools().FormattedOutput(varmvaMatB, theInputVars, *theVars, Log());
                  Log() << kINFO << Endl;
               } else
                  Log() << kWARNING << Form("Dataset[%s] : ", method->fDataSetInfo.GetName())
                        << "<TestAllMethods> cannot compute correlation matrices" << Endl;

               //              print overlap matrices
               Log() << kINFO << Form("Dataset[%s] : ", method->fDataSetInfo.GetName())
                     << "The following \"overlap\" matrices contain the fraction of events for which " << Endl;
               Log() << kINFO << Form("Dataset[%s] : ", method->fDataSetInfo.GetName())
                     << "the MVAs 'i' and 'j' have returned conform answers about \"signal-likeness\"" << Endl;
               Log() << kINFO << Form("Dataset[%s] : ", method->fDataSetInfo.GetName())
                     << "An event is signal-like, if its MVA output exceeds the following value:" << Endl;
               gTools().FormattedOutput(rvec, *theVars, "Method", "Cut value", Log());
               Log() << kINFO << Form("Dataset[%s] : ", method->fDataSetInfo.GetName())
                     << "which correspond to the working point: eff(signal) = 1 - eff(background)" << Endl;

               //              give notice that cut method has been excluded from this test
               if (nmeth != (Int_t)methods->size())
                  Log() << kINFO << Form("Dataset[%s] : ", method->fDataSetInfo.GetName())
                        << "Note: no correlations and overlap with cut method are provided at present" << Endl;

               if (nmeth > 1) {
                  Log() << kINFO << Endl;
                  Log() << kINFO << Form("Dataset[%s] : ", method->fDataSetInfo.GetName())
                        << "Inter-MVA overlap matrix (signal):" << Endl;
                  gTools().FormattedOutput(*overlapS, *theVars, Log());
                  Log() << kINFO << Endl;

                  Log() << kINFO << Form("Dataset[%s] : ", method->fDataSetInfo.GetName())
                        << "Inter-MVA overlap matrix (background):" << Endl;
                  gTools().FormattedOutput(*overlapB, *theVars, Log());
               }

               //              cleanup
               delete tpSig;
               delete tpBkg;
               delete corrMatS;
               delete corrMatB;
               delete theVars;
               delete overlapS;
               delete overlapB;
               delete[] dvec;
            }
         }
      }
      // -----------------------------------------------------------------------
      // Third part of evaluation process
      // --> output
      // -----------------------------------------------------------------------

      if (doRegression) {

         Log() << kINFO << Endl;
         TString hLine =
            "--------------------------------------------------------------------------------------------------";
         Log() << kINFO << "Evaluation results ranked by smallest RMS on test sample:" << Endl;
         Log() << kINFO << "(\"Bias\" quotes the mean deviation of the regression from true target." << Endl;
         Log() << kINFO << " \"MutInf\" is the \"Mutual Information\" between regression and target." << Endl;
         Log() << kINFO << " Indicated by \"_T\" are the corresponding \"truncated\" quantities ob-" << Endl;
         Log() << kINFO << " tained when removing events deviating more than 2sigma from average.)" << Endl;
         Log() << kINFO << hLine << Endl;
         // Log() << kINFO << "DataSet Name:        MVA Method:        <Bias>   <Bias_T>    RMS    RMS_T  |  MutInf
         // MutInf_T" << Endl;
         Log() << kINFO << hLine << Endl;

         for (Int_t i = 0; i < nmeth_used[0]; i++) {
            MethodBase *theMethod = dynamic_cast<MethodBase *>((*methods)[i]);
            if (theMethod == 0)
               continue;

            Log() << kINFO
                  << Form("%-20s %-15s:%#9.3g%#9.3g%#9.3g%#9.3g  |  %#5.3f  %#5.3f", theMethod->fDataSetInfo.GetName(),
                          (const char *)mname[0][i], biastest[0][i], biastestT[0][i], rmstest[0][i], rmstestT[0][i],
                          minftest[0][i], minftestT[0][i])
                  << Endl;
         }
         Log() << kINFO << hLine << Endl;
         Log() << kINFO << Endl;
         Log() << kINFO << "Evaluation results ranked by smallest RMS on training sample:" << Endl;
         Log() << kINFO << "(overtraining check)" << Endl;
         Log() << kINFO << hLine << Endl;
         Log() << kINFO
               << "DataSet Name:         MVA Method:        <Bias>   <Bias_T>    RMS    RMS_T  |  MutInf MutInf_T"
               << Endl;
         Log() << kINFO << hLine << Endl;

         for (Int_t i = 0; i < nmeth_used[0]; i++) {
            MethodBase *theMethod = dynamic_cast<MethodBase *>((*methods)[i]);
            if (theMethod == 0)
               continue;
            Log() << kINFO
                  << Form("%-20s %-15s:%#9.3g%#9.3g%#9.3g%#9.3g  |  %#5.3f  %#5.3f", theMethod->fDataSetInfo.GetName(),
                          (const char *)mname[0][i], biastrain[0][i], biastrainT[0][i], rmstrain[0][i], rmstrainT[0][i],
                          minftrain[0][i], minftrainT[0][i])
                  << Endl;
         }
         Log() << kINFO << hLine << Endl;
         Log() << kINFO << Endl;
      } else if (doMulticlass) {
         // ====================================================================
         // === Multiclass Output
         // ====================================================================

         TString hLine =
            "-------------------------------------------------------------------------------------------------------";

         // This part uses a genetic alg. to evaluate the optimal sig eff * sig pur.
         // This is why it is disabled for now.
         //
         // // --- Acheivable signal efficiency * signal purity
         // // --------------------------------------------------------------------
         // Log() << kINFO << Endl;
         // Log() << kINFO << "Evaluation results ranked by best signal efficiency times signal purity " << Endl;
         // Log() << kINFO << hLine << Endl;

         // // iterate over methods and evaluate
         // for (MVector::iterator itrMethod = methods->begin(); itrMethod != methods->end(); itrMethod++) {
         //    MethodBase *theMethod = dynamic_cast<MethodBase *>(*itrMethod);
         //    if (theMethod == 0) {
         //       continue;
         //    }

         //    TString header = "DataSet Name     MVA Method     ";
         //    for (UInt_t icls = 0; icls < theMethod->fDataSetInfo.GetNClasses(); ++icls) {
         //       header += Form("%-12s ", theMethod->fDataSetInfo.GetClassInfo(icls)->GetName());
         //    }

         //    Log() << kINFO << header << Endl;
         //    Log() << kINFO << hLine << Endl;
         //    for (Int_t i = 0; i < nmeth_used[0]; i++) {
         //       TString res = Form("[%-14s] %-15s", theMethod->fDataSetInfo.GetName(), (const char *)mname[0][i]);
         //       for (UInt_t icls = 0; icls < theMethod->fDataSetInfo.GetNClasses(); ++icls) {
         //          res += Form("%#1.3f        ", (multiclass_testEff[i][icls]) * (multiclass_testPur[i][icls]));
         //       }
         //       Log() << kINFO << res << Endl;
         //    }

         //    Log() << kINFO << hLine << Endl;
         //    Log() << kINFO << Endl;
         // }

         // --- 1 vs Rest ROC AUC, signal efficiency @ given background efficiency
         // --------------------------------------------------------------------
         TString header1 = Form("%-15s%-15s%-15s%-15s%-15s%-15s", "Dataset", "MVA Method", "ROC AUC", "Sig eff@B=0.01",
                                "Sig eff@B=0.10", "Sig eff@B=0.30");
         TString header2 = Form("%-15s%-15s%-15s%-15s%-15s%-15s", "Name:", "/ Class:", "test  (train)", "test  (train)",
                                "test  (train)", "test  (train)");
         Log() << kINFO << Endl;
         Log() << kINFO << "1-vs-rest performance metrics per class" << Endl;
         Log() << kINFO << hLine << Endl;
         Log() << kINFO << Endl;
         Log() << kINFO << "Considers the listed class as signal and the other classes" << Endl;
         Log() << kINFO << "as background, reporting the resulting binary performance." << Endl;
         Log() << kINFO << "A score of 0.820 (0.850) means 0.820 was acheived on the" << Endl;
         Log() << kINFO << "test set and 0.850 on the training set." << Endl;

         Log() << kINFO << Endl;
         Log() << kINFO << header1 << Endl;
         Log() << kINFO << header2 << Endl;
         for (Int_t k = 0; k < 2; k++) {
            for (Int_t i = 0; i < nmeth_used[k]; i++) {
               if (k == 1) {
                  mname[k][i].ReplaceAll("Variable_", "");
               }

               const TString datasetName = itrMap->first;
               const TString mvaName = mname[k][i];

               MethodBase *theMethod = dynamic_cast<MethodBase *>(GetMethod(datasetName, mvaName));
               if (theMethod == 0) {
                  continue;
               }

               Log() << kINFO << Endl;
               TString row = Form("%-15s%-15s", datasetName.Data(), mvaName.Data());
               Log() << kINFO << row << Endl;
               Log() << kINFO << "------------------------------" << Endl;

               UInt_t numClasses = theMethod->fDataSetInfo.GetNClasses();
               for (UInt_t iClass = 0; iClass < numClasses; ++iClass) {

                  ROCCurve *rocCurveTrain = GetROC(datasetName, mvaName, iClass, Types::kTraining);
                  ROCCurve *rocCurveTest = GetROC(datasetName, mvaName, iClass, Types::kTesting);

                  const TString className = theMethod->DataInfo().GetClassInfo(iClass)->GetName();
                  const Double_t rocaucTrain = rocCurveTrain->GetROCIntegral();
                  const Double_t effB01Train = rocCurveTrain->GetEffSForEffB(0.01);
                  const Double_t effB10Train = rocCurveTrain->GetEffSForEffB(0.10);
                  const Double_t effB30Train = rocCurveTrain->GetEffSForEffB(0.30);
                  const Double_t rocaucTest = rocCurveTest->GetROCIntegral();
                  const Double_t effB01Test = rocCurveTest->GetEffSForEffB(0.01);
                  const Double_t effB10Test = rocCurveTest->GetEffSForEffB(0.10);
                  const Double_t effB30Test = rocCurveTest->GetEffSForEffB(0.30);
                  const TString rocaucCmp = Form("%5.3f (%5.3f)", rocaucTest, rocaucTrain);
                  const TString effB01Cmp = Form("%5.3f (%5.3f)", effB01Test, effB01Train);
                  const TString effB10Cmp = Form("%5.3f (%5.3f)", effB10Test, effB10Train);
                  const TString effB30Cmp = Form("%5.3f (%5.3f)", effB30Test, effB30Train);
                  row = Form("%-15s%-15s%-15s%-15s%-15s%-15s", "", className.Data(), rocaucCmp.Data(), effB01Cmp.Data(),
                             effB10Cmp.Data(), effB30Cmp.Data());
                  Log() << kINFO << row << Endl;

                  delete rocCurveTrain;
                  delete rocCurveTest;
               }
            }
         }
         Log() << kINFO << Endl;
         Log() << kINFO << hLine << Endl;
         Log() << kINFO << Endl;

         // --- Confusion matrices
         // --------------------------------------------------------------------
         auto printMatrix = [](TMatrixD const &matTraining, TMatrixD const &matTesting, std::vector<TString> classnames,
                               UInt_t numClasses, MsgLogger &stream) {
            // assert (classLabledWidth >= valueLabelWidth + 2)
            // if (...) {Log() << kWARN << "..." << Endl; }

            // TODO: Ensure matrices are same size.

            TString header = Form(" %-14s", " ");
            TString headerInfo = Form(" %-14s", " ");
            ;
            for (UInt_t iCol = 0; iCol < numClasses; ++iCol) {
               header += Form(" %-14s", classnames[iCol].Data());
               headerInfo += Form(" %-14s", " test (train)");
            }
            stream << kINFO << header << Endl;
            stream << kINFO << headerInfo << Endl;

            for (UInt_t iRow = 0; iRow < numClasses; ++iRow) {
               stream << kINFO << Form(" %-14s", classnames[iRow].Data());

               for (UInt_t iCol = 0; iCol < numClasses; ++iCol) {
                  if (iCol == iRow) {
                     stream << kINFO << Form(" %-14s", "-");
                  } else {
                     Double_t trainValue = matTraining[iRow][iCol];
                     Double_t testValue = matTesting[iRow][iCol];
                     TString entry = Form("%-5.3f (%-5.3f)", testValue, trainValue);
                     stream << kINFO << Form(" %-14s", entry.Data());
                  }
               }
               stream << kINFO << Endl;
            }
         };

         Log() << kINFO << Endl;
         Log() << kINFO << "Confusion matrices for all methods" << Endl;
         Log() << kINFO << hLine << Endl;
         Log() << kINFO << Endl;
         Log() << kINFO << "Does a binary comparison between the two classes given by a " << Endl;
         Log() << kINFO << "particular row-column combination. In each case, the class " << Endl;
         Log() << kINFO << "given by the row is considered signal while the class given " << Endl;
         Log() << kINFO << "by the column index is considered background." << Endl;
         Log() << kINFO << Endl;
         for (UInt_t iMethod = 0; iMethod < methods->size(); ++iMethod) {
            MethodBase *theMethod = dynamic_cast<MethodBase *>(methods->at(iMethod));
            if (theMethod == nullptr) {
               continue;
            }
            UInt_t numClasses = theMethod->fDataSetInfo.GetNClasses();

            std::vector<TString> classnames;
            for (UInt_t iCls = 0; iCls < numClasses; ++iCls) {
               classnames.push_back(theMethod->fDataSetInfo.GetClassInfo(iCls)->GetName());
            }
            Log() << kINFO
                  << "=== Showing confusion matrix for method : " << Form("%-15s", (const char *)mname[0][iMethod])
                  << Endl;
            Log() << kINFO << "(Signal Efficiency for Background Efficiency 0.01%)" << Endl;
            Log() << kINFO << "---------------------------------------------------" << Endl;
            printMatrix(multiclass_testConfusionEffB01[iMethod], multiclass_trainConfusionEffB01[iMethod], classnames,
                        numClasses, Log());
            Log() << kINFO << Endl;

            Log() << kINFO << "(Signal Efficiency for Background Efficiency 0.10%)" << Endl;
            Log() << kINFO << "---------------------------------------------------" << Endl;
            printMatrix(multiclass_testConfusionEffB10[iMethod], multiclass_trainConfusionEffB10[iMethod], classnames,
                        numClasses, Log());
            Log() << kINFO << Endl;

            Log() << kINFO << "(Signal Efficiency for Background Efficiency 0.30%)" << Endl;
            Log() << kINFO << "---------------------------------------------------" << Endl;
            printMatrix(multiclass_testConfusionEffB30[iMethod], multiclass_trainConfusionEffB30[iMethod], classnames,
                        numClasses, Log());
            Log() << kINFO << Endl;
         }
         Log() << kINFO << hLine << Endl;
         Log() << kINFO << Endl;

      } else {
         // Binary classification
         if (fROC) {
            Log().EnableOutput();
            gConfig().SetSilent(kFALSE);
            Log() << Endl;
            TString hLine = "------------------------------------------------------------------------------------------"
                            "-------------------------";
            Log() << kINFO << "Evaluation results ranked by best signal efficiency and purity (area)" << Endl;
            Log() << kINFO << hLine << Endl;
            Log() << kINFO << "DataSet       MVA                       " << Endl;
            Log() << kINFO << "Name:         Method:          ROC-integ" << Endl;

            //       Log() << kDEBUG << "DataSet              MVA              Signal efficiency at bkg eff.(error):
            //       | Sepa-    Signifi- "   << Endl; Log() << kDEBUG << "Name:                Method:          @B=0.01
            //       @B=0.10    @B=0.30    ROC-integ    ROCCurve| ration:  cance:   "   << Endl;
            Log() << kDEBUG << hLine << Endl;
            for (Int_t k = 0; k < 2; k++) {
               if (k == 1 && nmeth_used[k] > 0) {
                  Log() << kINFO << hLine << Endl;
                  Log() << kINFO << "Input Variables: " << Endl << hLine << Endl;
               }
               for (Int_t i = 0; i < nmeth_used[k]; i++) {
                  TString datasetName = itrMap->first;
                  TString methodName = mname[k][i];

                  if (k == 1) {
                     methodName.ReplaceAll("Variable_", "");
                  }

                  MethodBase *theMethod = dynamic_cast<MethodBase *>(GetMethod(datasetName, methodName));
                  if (theMethod == 0) {
                     continue;
                  }

                  TMVA::DataSet *dataset = theMethod->Data();
                  TMVA::Results *results = dataset->GetResults(methodName, Types::kTesting, this->fAnalysisType);
                  std::vector<Bool_t> *mvaResType =
                     dynamic_cast<ResultsClassification *>(results)->GetValueVectorTypes();

                  Double_t rocIntegral = 0.0;
                  if (mvaResType->size() != 0) {
                     rocIntegral = GetROCIntegral(datasetName, methodName);
                  }

                  if (sep[k][i] < 0 || sig[k][i] < 0) {
                     // cannot compute separation/significance -> no MVA (usually for Cuts)
                     Log() << kINFO << Form("%-13s %-15s: %#1.3f", datasetName.Data(), methodName.Data(), effArea[k][i])
                           << Endl;

                     //               Log() << kDEBUG << Form("%-20s %-15s: %#1.3f(%02i)  %#1.3f(%02i)  %#1.3f(%02i)
                     //               %#1.3f       %#1.3f | --       --",
                     //                       datasetName.Data(),
                     //                       methodName.Data(),
                     //                       eff01[k][i], Int_t(1000*eff01err[k][i]),
                     //                       eff10[k][i], Int_t(1000*eff10err[k][i]),
                     //                       eff30[k][i], Int_t(1000*eff30err[k][i]),
                     //                       effArea[k][i],rocIntegral) << Endl;
                  } else {
                     Log() << kINFO << Form("%-13s %-15s: %#1.3f", datasetName.Data(), methodName.Data(), rocIntegral)
                           << Endl;
                     //               Log() << kDEBUG << Form("%-20s %-15s: %#1.3f(%02i)  %#1.3f(%02i)  %#1.3f(%02i)
                     //               %#1.3f       %#1.3f | %#1.3f    %#1.3f",
                     //                       datasetName.Data(),
                     //                       methodName.Data(),
                     //                       eff01[k][i], Int_t(1000*eff01err[k][i]),
                     //                       eff10[k][i], Int_t(1000*eff10err[k][i]),
                     //                       eff30[k][i], Int_t(1000*eff30err[k][i]),
                     //                       effArea[k][i],rocIntegral,
                     //                       sep[k][i], sig[k][i]) << Endl;
                  }
               }
            }
            Log() << kINFO << hLine << Endl;
            Log() << kINFO << Endl;
            Log() << kINFO << "Testing efficiency compared to training efficiency (overtraining check)" << Endl;
            Log() << kINFO << hLine << Endl;
            Log() << kINFO
                  << "DataSet              MVA              Signal efficiency: from test sample (from training sample) "
                  << Endl;
            Log() << kINFO << "Name:                Method:          @B=0.01             @B=0.10            @B=0.30   "
                  << Endl;
            Log() << kINFO << hLine << Endl;
            for (Int_t k = 0; k < 2; k++) {
               if (k == 1 && nmeth_used[k] > 0) {
                  Log() << kINFO << hLine << Endl;
                  Log() << kINFO << "Input Variables: " << Endl << hLine << Endl;
               }
               for (Int_t i = 0; i < nmeth_used[k]; i++) {
                  if (k == 1)
                     mname[k][i].ReplaceAll("Variable_", "");
                  MethodBase *theMethod = dynamic_cast<MethodBase *>((*methods)[i]);
                  if (theMethod == 0)
                     continue;

                  Log() << kINFO
                        << Form("%-20s %-15s: %#1.3f (%#1.3f)       %#1.3f (%#1.3f)      %#1.3f (%#1.3f)",
                                theMethod->fDataSetInfo.GetName(), (const char *)mname[k][i], eff01[k][i],
                                trainEff01[k][i], eff10[k][i], trainEff10[k][i], eff30[k][i], trainEff30[k][i])
                        << Endl;
               }
            }
            Log() << kINFO << hLine << Endl;
            Log() << kINFO << Endl;

            if (gTools().CheckForSilentOption(GetOptions()))
               Log().InhibitOutput();
         } // end fROC
      }
      if (!IsSilentFile()) {
         std::list<TString> datasets;
         for (Int_t k = 0; k < 2; k++) {
            for (Int_t i = 0; i < nmeth_used[k]; i++) {
               MethodBase *theMethod = dynamic_cast<MethodBase *>((*methods)[i]);
               if (theMethod == 0)
                  continue;
               // write test/training trees
               RootBaseDir()->cd(theMethod->fDataSetInfo.GetName());
               if (std::find(datasets.begin(), datasets.end(), theMethod->fDataSetInfo.GetName()) == datasets.end()) {
                  theMethod->fDataSetInfo.GetDataSet()->GetTree(Types::kTesting)->Write("", TObject::kOverwrite);
                  theMethod->fDataSetInfo.GetDataSet()->GetTree(Types::kTraining)->Write("", TObject::kOverwrite);
                  datasets.push_back(theMethod->fDataSetInfo.GetName());
               }
            }
         }
      }
   } // end for MethodsMap
   // references for citation
   gTools().TMVACitation(Log(), Tools::kHtmlLink);
}

////////////////////////////////////////////////////////////////////////////////
/// Evaluate Variable Importance

TH1F *TMVA::Factory::EvaluateImportance(DataLoader *loader, VIType vitype, Types::EMVA theMethod, TString methodTitle,
                                        const char *theOption)
{
   fModelPersistence = kFALSE;
   fSilentFile = kTRUE; // we need silent file here because we need fast classification results

   // getting number of variables and variable names from loader
   const int nbits = loader->GetDataSetInfo().GetNVariables();
   if (vitype == VIType::kShort)
      return EvaluateImportanceShort(loader, theMethod, methodTitle, theOption);
   else if (vitype == VIType::kAll)
      return EvaluateImportanceAll(loader, theMethod, methodTitle, theOption);
   else if (vitype == VIType::kRandom) {
      if ( nbits > 10 && nbits < 30) {
         // limit nbits to less than 30 to avoid error converting from double to uint and also cannot deal with too many combinations
         return EvaluateImportanceRandom(loader, static_cast<UInt_t>( pow(2, nbits) ), theMethod, methodTitle, theOption);
      } else if (nbits < 10)  {
         Log() << kERROR << "Error in Variable Importance: Random mode require more that 10 variables in the dataset."
                << Endl;
      } else if (nbits > 30) {
         Log() << kERROR << "Error in Variable Importance: Number of variables is too large for Random mode"
                  << Endl;
      }
   }
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////

TH1F *TMVA::Factory::EvaluateImportanceAll(DataLoader *loader, Types::EMVA theMethod, TString methodTitle,
                                           const char *theOption)
{

   uint64_t x = 0;
   uint64_t y = 0;

   // getting number of variables and variable names from loader
   const int nbits = loader->GetDataSetInfo().GetNVariables();
   std::vector<TString> varNames = loader->GetDataSetInfo().GetListOfVariables();

   if (nbits > 60) {
      Log() << kERROR << "Number of combinations is too large , is 2^" << nbits << Endl;
      return nullptr;
   }
   if (nbits > 20) {
      Log() << kWARNING << "Number of combinations is very large , is 2^" << nbits << Endl;
   }
   uint64_t range = static_cast<uint64_t>(pow(2, nbits));


   // vector to save importances
   std::vector<Double_t> importances(nbits);
   // vector to save ROC
   std::vector<Double_t> ROC(range);
   ROC[0] = 0.5;
   for (int i = 0; i < nbits; i++)
      importances[i] = 0;

   Double_t SROC, SSROC; // computed ROC value
   for (x = 1; x < range; x++) {

      std::bitset<VIBITS> xbitset(x);
      if (x == 0)
         continue; // data loader need at least one variable

      // creating loader for seed
      TMVA::DataLoader *seedloader = new TMVA::DataLoader(xbitset.to_string());

      // adding variables from seed
      for (int index = 0; index < nbits; index++) {
         if (xbitset[index])
            seedloader->AddVariable(varNames[index], 'F');
      }

      DataLoaderCopy(seedloader, loader);
      seedloader->PrepareTrainingAndTestTree(loader->GetDataSetInfo().GetCut("Signal"),
                                             loader->GetDataSetInfo().GetCut("Background"),
                                             loader->GetDataSetInfo().GetSplitOptions());

      // Booking Seed
      BookMethod(seedloader, theMethod, methodTitle, theOption);

      // Train/Test/Evaluation
      TrainAllMethods();
      TestAllMethods();
      EvaluateAllMethods();

      // getting ROC
      ROC[x] = GetROCIntegral(xbitset.to_string(), methodTitle);

      // cleaning information to process sub-seeds
      TMVA::MethodBase *smethod = dynamic_cast<TMVA::MethodBase *>(fMethodsMap[xbitset.to_string().c_str()][0][0]);
      TMVA::ResultsClassification *sresults = (TMVA::ResultsClassification *)smethod->Data()->GetResults(
         smethod->GetMethodName(), Types::kTesting, Types::kClassification);
      delete sresults;
      delete seedloader;
      this->DeleteAllMethods();

      fMethodsMap.clear();
      // removing global result because it is requiring a lot of RAM for all seeds
   }

   for (x = 0; x < range; x++) {
      SROC = ROC[x];
      for (uint32_t i = 0; i < VIBITS; ++i) {
         if (x & (uint64_t(1) << i)) {
            y = x & ~(1 << i);
            std::bitset<VIBITS> ybitset(y);
            // need at least one variable
            // NOTE: if sub-seed is zero then is the special case
            // that count in xbitset is 1
            uint32_t ny = static_cast<uint32_t>( log(x - y) / 0.693147 ) ;
            if (y == 0) {
               importances[ny] = SROC - 0.5;
               continue;
            }

            // getting ROC
            SSROC = ROC[y];
            importances[ny] += SROC - SSROC;
            // cleaning information
         }
      }
   }
   std::cout << "--- Variable Importance Results (All)" << std::endl;
   return GetImportance(nbits, importances, varNames);
}

static uint64_t sum(uint64_t i)
{
   // add a limit for overflows
   if (i > 62) return 0;
   return static_cast<uint64_t>( std::pow(2, i + 1)) - 1;
   // uint64_t _sum = 0;
   // for (uint64_t n = 0; n < i; n++)
   //    _sum += pow(2, n);
   // return _sum;
}

////////////////////////////////////////////////////////////////////////////////

TH1F *TMVA::Factory::EvaluateImportanceShort(DataLoader *loader, Types::EMVA theMethod, TString methodTitle,
                                             const char *theOption)
{
   uint64_t x = 0;
   uint64_t y = 0;

   // getting number of variables and variable names from loader
   const int nbits = loader->GetDataSetInfo().GetNVariables();
   std::vector<TString> varNames = loader->GetDataSetInfo().GetListOfVariables();

   if (nbits > 60) {
      Log() << kERROR << "Number of combinations is too large , is 2^" << nbits << Endl;
      return nullptr;
   }
   long int range = sum(nbits);
   //   std::cout<<range<<std::endl;
   // vector to save importances
   std::vector<Double_t> importances(nbits);
   for (int i = 0; i < nbits; i++)
      importances[i] = 0;

   Double_t SROC, SSROC; // computed ROC value

   x = range;

   std::bitset<VIBITS> xbitset(x);
   if (x == 0)
      Log() << kFATAL << "Error: need at least one variable."; // data loader need at least one variable

   // creating loader for seed
   TMVA::DataLoader *seedloader = new TMVA::DataLoader(xbitset.to_string());

   // adding variables from seed
   for (int index = 0; index < nbits; index++) {
      if (xbitset[index])
         seedloader->AddVariable(varNames[index], 'F');
   }

   // Loading Dataset
   DataLoaderCopy(seedloader, loader);

   // Booking Seed
   BookMethod(seedloader, theMethod, methodTitle, theOption);

   // Train/Test/Evaluation
   TrainAllMethods();
   TestAllMethods();
   EvaluateAllMethods();

   // getting ROC
   SROC = GetROCIntegral(xbitset.to_string(), methodTitle);

   // cleaning information to process sub-seeds
   TMVA::MethodBase *smethod = dynamic_cast<TMVA::MethodBase *>(fMethodsMap[xbitset.to_string().c_str()][0][0]);
   TMVA::ResultsClassification *sresults = (TMVA::ResultsClassification *)smethod->Data()->GetResults(
      smethod->GetMethodName(), Types::kTesting, Types::kClassification);
   delete sresults;
   delete seedloader;
   this->DeleteAllMethods();
   fMethodsMap.clear();

   // removing global result because it is requiring a lot of RAM for all seeds

   for (uint32_t i = 0; i < VIBITS; ++i) {
      if (x & (1 << i)) {
         y = x & ~(uint64_t(1) << i);
         std::bitset<VIBITS> ybitset(y);
         // need at least one variable
         // NOTE: if sub-seed is zero then is the special case
         // that count in xbitset is 1
         uint32_t ny = static_cast<uint32_t>(log(x - y) / 0.693147);
         if (y == 0) {
            importances[ny] = SROC - 0.5;
            continue;
         }

         // creating loader for sub-seed
         TMVA::DataLoader *subseedloader = new TMVA::DataLoader(ybitset.to_string());
         // adding variables from sub-seed
         for (int index = 0; index < nbits; index++) {
            if (ybitset[index])
               subseedloader->AddVariable(varNames[index], 'F');
         }

         // Loading Dataset
         DataLoaderCopy(subseedloader, loader);

         // Booking SubSeed
         BookMethod(subseedloader, theMethod, methodTitle, theOption);

         // Train/Test/Evaluation
         TrainAllMethods();
         TestAllMethods();
         EvaluateAllMethods();

         // getting ROC
         SSROC = GetROCIntegral(ybitset.to_string(), methodTitle);
         importances[ny] += SROC - SSROC;

         // cleaning information
         TMVA::MethodBase *ssmethod = dynamic_cast<TMVA::MethodBase *>(fMethodsMap[ybitset.to_string().c_str()][0][0]);
         TMVA::ResultsClassification *ssresults = (TMVA::ResultsClassification *)ssmethod->Data()->GetResults(
            ssmethod->GetMethodName(), Types::kTesting, Types::kClassification);
         delete ssresults;
         delete subseedloader;
         this->DeleteAllMethods();
         fMethodsMap.clear();
      }
   }
   std::cout << "--- Variable Importance Results (Short)" << std::endl;
   return GetImportance(nbits, importances, varNames);
}

////////////////////////////////////////////////////////////////////////////////

TH1F *TMVA::Factory::EvaluateImportanceRandom(DataLoader *loader, UInt_t nseeds, Types::EMVA theMethod,
                                              TString methodTitle, const char *theOption)
{
   TRandom3 *rangen = new TRandom3(0); // Random Gen.

   uint64_t x = 0;
   uint64_t y = 0;

   // getting number of variables and variable names from loader
   const int nbits = loader->GetDataSetInfo().GetNVariables();
   std::vector<TString> varNames = loader->GetDataSetInfo().GetListOfVariables();

   long int range = pow(2, nbits);

   // vector to save importances
   std::vector<Double_t> importances(nbits);
   Double_t importances_norm = 0;
   for (int i = 0; i < nbits; i++)
      importances[i] = 0;

   Double_t SROC, SSROC; // computed ROC value
   for (UInt_t n = 0; n < nseeds; n++) {
      x = rangen->Integer(range);

      std::bitset<32> xbitset(x);
      if (x == 0)
         continue; // data loader need at least one variable

      // creating loader for seed
      TMVA::DataLoader *seedloader = new TMVA::DataLoader(xbitset.to_string());

      // adding variables from seed
      for (int index = 0; index < nbits; index++) {
         if (xbitset[index])
            seedloader->AddVariable(varNames[index], 'F');
      }

      // Loading Dataset
      DataLoaderCopy(seedloader, loader);

      // Booking Seed
      BookMethod(seedloader, theMethod, methodTitle, theOption);

      // Train/Test/Evaluation
      TrainAllMethods();
      TestAllMethods();
      EvaluateAllMethods();

      // getting ROC
      SROC = GetROCIntegral(xbitset.to_string(), methodTitle);
      //       std::cout << "Seed: n " << n << " x " << x << " xbitset:" << xbitset << "  ROC " << SROC << std::endl;

      // cleaning information to process sub-seeds
      TMVA::MethodBase *smethod = dynamic_cast<TMVA::MethodBase *>(fMethodsMap[xbitset.to_string().c_str()][0][0]);
      TMVA::ResultsClassification *sresults = (TMVA::ResultsClassification *)smethod->Data()->GetResults(
         smethod->GetMethodName(), Types::kTesting, Types::kClassification);
      delete sresults;
      delete seedloader;
      this->DeleteAllMethods();
      fMethodsMap.clear();

      // removing global result because it is requiring a lot of RAM for all seeds

      for (uint32_t i = 0; i < 32; ++i) {
         if (x & (uint64_t(1) << i)) {
            y = x & ~(1 << i);
            std::bitset<32> ybitset(y);
            // need at least one variable
            // NOTE: if sub-seed is zero then is the special case
            // that count in xbitset is 1
            Double_t ny = log(x - y) / 0.693147;
            if (y == 0) {
               importances[ny] = SROC - 0.5;
               importances_norm += importances[ny];
               //  std::cout << "SubSeed: " << y << " y:" << ybitset << "ROC " << 0.5 << std::endl;
               continue;
            }

            // creating loader for sub-seed
            TMVA::DataLoader *subseedloader = new TMVA::DataLoader(ybitset.to_string());
            // adding variables from sub-seed
            for (int index = 0; index < nbits; index++) {
               if (ybitset[index])
                  subseedloader->AddVariable(varNames[index], 'F');
            }

            // Loading Dataset
            DataLoaderCopy(subseedloader, loader);

            // Booking SubSeed
            BookMethod(subseedloader, theMethod, methodTitle, theOption);

            // Train/Test/Evaluation
            TrainAllMethods();
            TestAllMethods();
            EvaluateAllMethods();

            // getting ROC
            SSROC = GetROCIntegral(ybitset.to_string(), methodTitle);
            importances[ny] += SROC - SSROC;
            // std::cout << "SubSeed: " << y << " y:" << ybitset << " x-y " << x - y << " " << std::bitset<32>(x - y) <<
            // " ny " << ny << " SROC " << SROC << " SSROC " << SSROC << " Importance = " << importances[ny] <<
            // std::endl; cleaning information
            TMVA::MethodBase *ssmethod =
               dynamic_cast<TMVA::MethodBase *>(fMethodsMap[ybitset.to_string().c_str()][0][0]);
            TMVA::ResultsClassification *ssresults = (TMVA::ResultsClassification *)ssmethod->Data()->GetResults(
               ssmethod->GetMethodName(), Types::kTesting, Types::kClassification);
            delete ssresults;
            delete subseedloader;
            this->DeleteAllMethods();
            fMethodsMap.clear();
         }
      }
   }
   std::cout << "--- Variable Importance Results (Random)" << std::endl;
   return GetImportance(nbits, importances, varNames);
}

////////////////////////////////////////////////////////////////////////////////

TH1F *TMVA::Factory::GetImportance(const int nbits, std::vector<Double_t> importances, std::vector<TString> varNames)
{
   TH1F *vih1 = new TH1F("vih1", "", nbits, 0, nbits);

   gStyle->SetOptStat(000000);

   Float_t normalization = 0.0;
   for (int i = 0; i < nbits; i++) {
      normalization = normalization + importances[i];
   }

   Float_t roc = 0.0;

   gStyle->SetTitleXOffset(0.4);
   gStyle->SetTitleXOffset(1.2);

   std::vector<Double_t> x_ie(nbits), y_ie(nbits);
   for (Int_t i = 1; i < nbits + 1; i++) {
      x_ie[i - 1] = (i - 1) * 1.;
      roc = 100.0 * importances[i - 1] / normalization;
      y_ie[i - 1] = roc;
      std::cout << "--- " << varNames[i - 1] << " = " << roc << " %" << std::endl;
      vih1->GetXaxis()->SetBinLabel(i, varNames[i - 1].Data());
      vih1->SetBinContent(i, roc);
   }
   TGraph *g_ie = new TGraph(nbits + 2, &x_ie[0], &y_ie[0]);
   g_ie->SetTitle("");

   vih1->LabelsOption("v >", "X");
   vih1->SetBarWidth(0.97);
   Int_t ca = TColor::GetColor("#006600");
   vih1->SetFillColor(ca);
   // Int_t ci = TColor::GetColor("#990000");

   vih1->GetYaxis()->SetTitle("Importance (%)");
   vih1->GetYaxis()->SetTitleSize(0.045);
   vih1->GetYaxis()->CenterTitle();
   vih1->GetYaxis()->SetTitleOffset(1.24);

   vih1->GetYaxis()->SetRangeUser(-7, 50);
   vih1->SetDirectory(0);

   //   vih1->Draw("B");
   return vih1;
}
