// @(#)root/tmva $Id$
// Author: Peter Speckmayer

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : MethodDNN                                                             *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      A neural network implementation                                           *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Simon Pfreundschuh    <s.pfreundschuh@gmail.com> - CERN, Switzerland      *
 *      Peter Speckmayer      <peter.speckmayer@gmx.ch>  - CERN, Switzerland      *
 *                                                                                *
 * Copyright (c) 2005-2015:                                                       *
 *      CERN, Switzerland                                                         *
 *      U. of Victoria, Canada                                                    *
 *      MPI-K Heidelberg, Germany                                                 *
 *      U. of Bonn, Germany                                                       *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

/*! \class TMVA::MethodDNN
\ingroup TMVA
Deep Neural Network Implementation.
*/

#include "TMVA/MethodDNN.h"

#include "TString.h"
#include "TFormula.h"
#include "TObjString.h"

#include "TMVA/ClassifierFactory.h"
#include "TMVA/Configurable.h"
#include "TMVA/IMethod.h"
#include "TMVA/MsgLogger.h"
#include "TMVA/MethodBase.h"
#include "TMVA/Timer.h"
#include "TMVA/Types.h"
#include "TMVA/Tools.h"
#include "TMVA/Config.h"
#include "TMVA/Ranking.h"

#include "TMVA/DNN/Net.h"
#include "TMVA/DNN/Architectures/Reference.h"

#include "TMVA/NeuralNet.h"
#include "TMVA/Monitoring.h"

#ifdef R__HAS_TMVACPU
#include "TMVA/DNN/Architectures/Cpu.h"
#endif
#ifdef R__HAS_TMVAGPU
#include "TMVA/DNN/Architectures/Cuda.h"
#endif

#include <algorithm>
#include <iostream>
#include <string>
#include <iomanip>

REGISTER_METHOD(DNN)

ClassImp(TMVA::MethodDNN);

namespace TMVA
{
   using namespace DNN;

   ////////////////////////////////////////////////////////////////////////////////
   /// standard constructor

   TMVA::MethodDNN::MethodDNN(const TString &jobName, const TString &methodTitle, DataSetInfo &theData,
                              const TString &theOption)
      : MethodBase(jobName, Types::kDNN, methodTitle, theData, theOption), fWeightInitialization(), fOutputFunction(),
        fLayoutString(), fErrorStrategy(), fTrainingStrategyString(), fWeightInitializationString(),
        fArchitectureString(), fTrainingSettings(), fResume(false), fSettings()
   {
}

////////////////////////////////////////////////////////////////////////////////
/// constructor from a weight file

TMVA::MethodDNN::MethodDNN(DataSetInfo& theData,
                           const TString& theWeightFile)
    : MethodBase( Types::kDNN, theData, theWeightFile),
     fWeightInitialization(), fOutputFunction(), fLayoutString(), fErrorStrategy(),
     fTrainingStrategyString(), fWeightInitializationString(), fArchitectureString(),
     fTrainingSettings(), fResume(false), fSettings()
{
        fWeightInitialization = DNN::EInitialization::kGauss;
        fOutputFunction = DNN::EOutputFunction::kSigmoid;
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TMVA::MethodDNN::~MethodDNN()
{
        fWeightInitialization = DNN::EInitialization::kGauss;
        fOutputFunction = DNN::EOutputFunction::kSigmoid;
}

////////////////////////////////////////////////////////////////////////////////
/// MLP can handle classification with 2 classes and regression with
/// one regression-target

Bool_t TMVA::MethodDNN::HasAnalysisType(Types::EAnalysisType type,
                                        UInt_t numberClasses,
                                        UInt_t /*numberTargets*/ )
{
   if (type == Types::kClassification && numberClasses == 2 ) return kTRUE;
   if (type == Types::kMulticlass ) return kTRUE;
   if (type == Types::kRegression ) return kTRUE;

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// default initializations

void TMVA::MethodDNN::Init() {
      Log() << kWARNING
            << "MethodDNN is deprecated and it will be removed in future ROOT version. "
               "Please use MethodDL ( TMVA::kDL)"
            << Endl;

}

////////////////////////////////////////////////////////////////////////////////
/// Options to be set in the option string:
///
///  - LearningRate    <float>      DNN learning rate parameter.
///  - DecayRate       <float>      Decay rate for learning parameter.
///  - TestRate        <int>        Period of validation set error computation.
///  - BatchSize       <int>        Number of event per batch.
///
///  - ValidationSize  <string>     How many events to use for validation. "0.2"
///                                 or "20%" indicates that a fifth of the
///                                 training data should be used. "100"
///                                 indicates that 100 events should be used.

void TMVA::MethodDNN::DeclareOptions()
{

   DeclareOptionRef(fLayoutString="SOFTSIGN|(N+100)*2,LINEAR",
                                  "Layout",
                                  "Layout of the network.");

   DeclareOptionRef(fValidationSize = "20%", "ValidationSize",
                    "Part of the training data to use for "
                    "validation. Specify as 0.2 or 20% to use a "
                    "fifth of the data set as validation set. "
                    "Specify as 100 to use exactly 100 events. "
                    "(Default: 20%)");

   DeclareOptionRef(fErrorStrategy="CROSSENTROPY",
                    "ErrorStrategy",
                    "Loss function: Mean squared error (regression)"
                    " or cross entropy (binary classification).");
   AddPreDefVal(TString("CROSSENTROPY"));
   AddPreDefVal(TString("SUMOFSQUARES"));
   AddPreDefVal(TString("MUTUALEXCLUSIVE"));

   DeclareOptionRef(fWeightInitializationString="XAVIER",
                    "WeightInitialization",
                    "Weight initialization strategy");
   AddPreDefVal(TString("XAVIER"));
   AddPreDefVal(TString("XAVIERUNIFORM"));

   DeclareOptionRef(fArchitectureString = "CPU", "Architecture", "Which architecture to perform the training on.");
   AddPreDefVal(TString("STANDARD"));
   AddPreDefVal(TString("CPU"));
   AddPreDefVal(TString("GPU"));
   AddPreDefVal(TString("OPENCL"));

   DeclareOptionRef(
       fTrainingStrategyString = "LearningRate=1e-1,"
                                 "Momentum=0.3,"
                                 "Repetitions=3,"
                                 "ConvergenceSteps=50,"
                                 "BatchSize=30,"
                                 "TestRepetitions=7,"
                                 "WeightDecay=0.0,"
                                 "Renormalize=L2,"
                                 "DropConfig=0.0,"
                                 "DropRepetitions=5|LearningRate=1e-4,"
                                 "Momentum=0.3,"
                                 "Repetitions=3,"
                                 "ConvergenceSteps=50,"
                                 "BatchSize=20,"
                                 "TestRepetitions=7,"
                                 "WeightDecay=0.001,"
                                 "Renormalize=L2,"
                                 "DropConfig=0.0+0.5+0.5,"
                                 "DropRepetitions=5,"
                                 "Multithreading=True",
                                 "TrainingStrategy",
                                 "Defines the training strategies.");
}

////////////////////////////////////////////////////////////////////////////////
/// parse layout specification string and return a vector, each entry
/// containing the number of neurons to go in each successive layer

auto TMVA::MethodDNN::ParseLayoutString(TString layoutString)
    -> LayoutVector_t
{
   LayoutVector_t layout;
   const TString layerDelimiter(",");
   const TString subDelimiter("|");

   const size_t inputSize = GetNvar();

   TObjArray* layerStrings = layoutString.Tokenize(layerDelimiter);
   TIter       nextLayer (layerStrings);
   TObjString* layerString = (TObjString*)nextLayer ();

   for (; layerString != nullptr; layerString = (TObjString*) nextLayer()) {
      int numNodes = 0;
      EActivationFunction activationFunction = EActivationFunction::kTanh;

      TObjArray* subStrings = layerString->GetString().Tokenize(subDelimiter);
      TIter nextToken (subStrings);
      TObjString* token = (TObjString *) nextToken();
      int idxToken = 0;
      for (; token != nullptr; token = (TObjString *) nextToken()) {
         switch (idxToken)
         {
         case 0:
         {
            TString strActFnc (token->GetString ());
            if (strActFnc == "RELU") {
                activationFunction = DNN::EActivationFunction::kRelu;
            } else if (strActFnc == "TANH") {
                activationFunction = DNN::EActivationFunction::kTanh;
            } else if (strActFnc == "SYMMRELU") {
                activationFunction = DNN::EActivationFunction::kSymmRelu;
            } else if (strActFnc == "SOFTSIGN") {
                activationFunction = DNN::EActivationFunction::kSoftSign;
            } else if (strActFnc == "SIGMOID") {
                activationFunction = DNN::EActivationFunction::kSigmoid;
            } else if (strActFnc == "LINEAR") {
                activationFunction = DNN::EActivationFunction::kIdentity;
            } else if (strActFnc == "GAUSS") {
                activationFunction = DNN::EActivationFunction::kGauss;
            }
         }
         break;
         case 1: // number of nodes
         {
            TString strNumNodes (token->GetString ());
            TString strN ("x");
            strNumNodes.ReplaceAll ("N", strN);
            strNumNodes.ReplaceAll ("n", strN);
            TFormula fml ("tmp",strNumNodes);
            numNodes = fml.Eval (inputSize);
         }
         break;
         }
         ++idxToken;
      }
      layout.push_back(std::make_pair(numNodes, activationFunction));
      }
   return layout;
}

////////////////////////////////////////////////////////////////////////////////
/// parse key value pairs in blocks -> return vector of blocks with map of key value pairs

auto TMVA::MethodDNN::ParseKeyValueString(TString parseString,
                                          TString blockDelim,
                                          TString tokenDelim)
    -> KeyValueVector_t
{
   KeyValueVector_t blockKeyValues;
   const TString keyValueDelim ("=");

   TObjArray* blockStrings = parseString.Tokenize (blockDelim);
   TIter nextBlock (blockStrings);
   TObjString* blockString = (TObjString *) nextBlock();

   for (; blockString != nullptr; blockString = (TObjString *) nextBlock())
   {
      blockKeyValues.push_back (std::map<TString,TString>());
      std::map<TString,TString>& currentBlock = blockKeyValues.back ();

      TObjArray* subStrings = blockString->GetString ().Tokenize (tokenDelim);
      TIter nextToken (subStrings);
      TObjString* token = (TObjString*)nextToken ();

      for (; token != nullptr; token = (TObjString *)nextToken())
      {
         TString strKeyValue (token->GetString ());
         int delimPos = strKeyValue.First (keyValueDelim.Data ());
         if (delimPos <= 0)
             continue;

         TString strKey = TString (strKeyValue (0, delimPos));
         strKey.ToUpper();
         TString strValue = TString (strKeyValue (delimPos+1, strKeyValue.Length ()));

         strKey.Strip (TString::kBoth, ' ');
         strValue.Strip (TString::kBoth, ' ');

         currentBlock.insert (std::make_pair (strKey, strValue));
      }
   }
   return blockKeyValues;
}

////////////////////////////////////////////////////////////////////////////////

TString fetchValue (const std::map<TString, TString>& keyValueMap, TString key)
{
   key.ToUpper ();
   std::map<TString, TString>::const_iterator it = keyValueMap.find (key);
   if (it == keyValueMap.end()) {
      return TString ("");
   }
   return it->second;
}

////////////////////////////////////////////////////////////////////////////////

template <typename T>
T fetchValue(const std::map<TString,TString>& keyValueMap,
              TString key,
              T defaultValue);

////////////////////////////////////////////////////////////////////////////////

template <>
int fetchValue(const std::map<TString,TString>& keyValueMap,
               TString key,
               int defaultValue)
{
   TString value (fetchValue (keyValueMap, key));
   if (value == "") {
      return defaultValue;
   }
   return value.Atoi ();
}

////////////////////////////////////////////////////////////////////////////////

template <>
double fetchValue (const std::map<TString,TString>& keyValueMap,
                   TString key, double defaultValue)
{
   TString value (fetchValue (keyValueMap, key));
   if (value == "") {
      return defaultValue;
   }
   return value.Atof ();
}

////////////////////////////////////////////////////////////////////////////////

template <>
TString fetchValue (const std::map<TString,TString>& keyValueMap,
                    TString key, TString defaultValue)
{
   TString value (fetchValue (keyValueMap, key));
   if (value == "") {
      return defaultValue;
   }
   return value;
}

////////////////////////////////////////////////////////////////////////////////

template <>
bool fetchValue (const std::map<TString,TString>& keyValueMap,
                 TString key, bool defaultValue)
{
   TString value (fetchValue (keyValueMap, key));
   if (value == "") {
      return defaultValue;
   }
   value.ToUpper ();
   if (value == "TRUE" || value == "T" || value == "1") {
      return true;
   }
   return false;
}

////////////////////////////////////////////////////////////////////////////////

template <>
std::vector<double> fetchValue(const std::map<TString, TString> & keyValueMap,
                               TString key,
                               std::vector<double> defaultValue)
{
   TString parseString (fetchValue (keyValueMap, key));
   if (parseString == "") {
      return defaultValue;
   }
   parseString.ToUpper ();
   std::vector<double> values;

   const TString tokenDelim ("+");
   TObjArray* tokenStrings = parseString.Tokenize (tokenDelim);
   TIter nextToken (tokenStrings);
   TObjString* tokenString = (TObjString*)nextToken ();
   for (; tokenString != NULL; tokenString = (TObjString*)nextToken ()) {
      std::stringstream sstr;
      double currentValue;
      sstr << tokenString->GetString ().Data ();
      sstr >> currentValue;
      values.push_back (currentValue);
   }
   return values;
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::MethodDNN::ProcessOptions()
{
   if (IgnoreEventsWithNegWeightsInTraining()) {
      Log() << kINFO
            << "Will ignore negative events in training!"
            << Endl;
   }

   if (fArchitectureString == "STANDARD") {
      Log() << kERROR << "The STANDARD architecture has been deprecated. "
                         "Please use Architecture=CPU or Architecture=CPU."
                         "See the TMVA Users' Guide for instructions if you "
                         "encounter problems."
            << Endl;
      Log() << kFATAL << "The STANDARD architecture has been deprecated. "
                         "Please use Architecture=CPU or Architecture=CPU."
                         "See the TMVA Users' Guide for instructions if you "
                         "encounter problems."
            << Endl;
   }

   if (fArchitectureString == "OPENCL") {
      Log() << kERROR << "The OPENCL architecture has not been implemented yet. "
                         "Please use Architecture=CPU or Architecture=CPU for the "
                         "time being. See the TMVA Users' Guide for instructions "
                         "if you encounter problems."
            << Endl;
      Log() << kFATAL << "The OPENCL architecture has not been implemented yet. "
                         "Please use Architecture=CPU or Architecture=CPU for the "
                         "time being. See the TMVA Users' Guide for instructions "
                         "if you encounter problems."
            << Endl;
   }

   if (fArchitectureString == "GPU") {
#ifndef DNNCUDA // Included only if DNNCUDA flag is _not_ set.
      Log() << kERROR << "CUDA backend not enabled. Please make sure "
                         "you have CUDA installed and it was successfully "
                         "detected by CMAKE."
            << Endl;
      Log() << kFATAL << "CUDA backend not enabled. Please make sure "
                         "you have CUDA installed and it was successfully "
                         "detected by CMAKE."
            << Endl;
#endif // DNNCUDA
   }

   if (fArchitectureString == "CPU") {
#ifndef DNNCPU // Included only if DNNCPU flag is _not_ set.
      Log() << kERROR << "Multi-core CPU backend not enabled. Please make sure "
                         "you have a BLAS implementation and it was successfully "
                         "detected by CMake as well that the imt CMake flag is set."
            << Endl;
      Log() << kFATAL << "Multi-core CPU backend not enabled. Please make sure "
                         "you have a BLAS implementation and it was successfully "
                         "detected by CMake as well that the imt CMake flag is set."
            << Endl;
#endif // DNNCPU
   }

   //
   // Set network structure.
   //

   fLayout = TMVA::MethodDNN::ParseLayoutString (fLayoutString);
   size_t inputSize = GetNVariables ();
   size_t outputSize = 1;
   if (fAnalysisType == Types::kRegression && GetNTargets() != 0) {
      outputSize = GetNTargets();
   } else if (fAnalysisType == Types::kMulticlass && DataInfo().GetNClasses() >= 2) {
      outputSize = DataInfo().GetNClasses();
   }

   fNet.SetBatchSize(1);
   fNet.SetInputWidth(inputSize);

   auto itLayout    = std::begin (fLayout);
   auto itLayoutEnd = std::end (fLayout)-1;
   for ( ; itLayout != itLayoutEnd; ++itLayout) {
      fNet.AddLayer((*itLayout).first, (*itLayout).second);
   }
   fNet.AddLayer(outputSize, EActivationFunction::kIdentity);

   //
   // Loss function and output.
   //

   fOutputFunction = EOutputFunction::kSigmoid;
   if (fAnalysisType == Types::kClassification)
   {
      if (fErrorStrategy == "SUMOFSQUARES") {
         fNet.SetLossFunction(ELossFunction::kMeanSquaredError);
      }
      if (fErrorStrategy == "CROSSENTROPY") {
         fNet.SetLossFunction(ELossFunction::kCrossEntropy);
      }
      fOutputFunction = EOutputFunction::kSigmoid;
   } else if (fAnalysisType == Types::kRegression) {
      if (fErrorStrategy != "SUMOFSQUARES") {
         Log () << kWARNING << "For regression only SUMOFSQUARES is a valid "
                << " neural net error function. Setting error function to "
                << " SUMOFSQUARES now." << Endl;
      }
      fNet.SetLossFunction(ELossFunction::kMeanSquaredError);
      fOutputFunction = EOutputFunction::kIdentity;
   } else if (fAnalysisType == Types::kMulticlass) {
      if (fErrorStrategy == "SUMOFSQUARES") {
         fNet.SetLossFunction(ELossFunction::kMeanSquaredError);
      }
      if (fErrorStrategy == "CROSSENTROPY") {
         fNet.SetLossFunction(ELossFunction::kCrossEntropy);
      }
      if (fErrorStrategy == "MUTUALEXCLUSIVE") {
         fNet.SetLossFunction(ELossFunction::kSoftmaxCrossEntropy);
      }
      fOutputFunction = EOutputFunction::kSoftmax;
   }

   //
   // Initialization
   //

   if (fWeightInitializationString == "XAVIER") {
      fWeightInitialization = DNN::EInitialization::kGauss;
   }
   else if (fWeightInitializationString == "XAVIERUNIFORM") {
      fWeightInitialization = DNN::EInitialization::kUniform;
   }
   else {
      fWeightInitialization = DNN::EInitialization::kGauss;
   }

   //
   // Training settings.
   //

   // Force validation of the ValidationSize option
   GetNumValidationSamples();

   KeyValueVector_t strategyKeyValues = ParseKeyValueString(fTrainingStrategyString,
                                                            TString ("|"),
                                                            TString (","));

   std::cout << "Parsed Training DNN string " << fTrainingStrategyString << std::endl;
   std::cout << "STring has size " << strategyKeyValues.size() << std::endl;
   for (auto& block : strategyKeyValues) {
      TTrainingSettings settings;

      settings.convergenceSteps = fetchValue(block, "ConvergenceSteps", 100);
      settings.batchSize        = fetchValue(block, "BatchSize", 30);
      settings.testInterval     = fetchValue(block, "TestRepetitions", 7);
      settings.weightDecay      = fetchValue(block, "WeightDecay", 0.0);
      settings.learningRate         = fetchValue(block, "LearningRate", 1e-5);
      settings.momentum             = fetchValue(block, "Momentum", 0.3);
      settings.dropoutProbabilities = fetchValue(block, "DropConfig",
                                                 std::vector<Double_t>());

      TString regularization = fetchValue(block, "Regularization",
                                          TString ("NONE"));
      if (regularization == "L1") {
         settings.regularization = DNN::ERegularization::kL1;
      } else if (regularization == "L2") {
         settings.regularization = DNN::ERegularization::kL2;
      } else {
         settings.regularization = DNN::ERegularization::kNone;
      }

      TString strMultithreading = fetchValue(block, "Multithreading",
                                             TString ("True"));
      if (strMultithreading.BeginsWith ("T")) {
         settings.multithreading = true;
      } else {
         settings.multithreading = false;
      }

      fTrainingSettings.push_back(settings);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Validation of the ValidationSize option. Allowed formats are 20%, 0.2 and
/// 100 etc.
///    - 20% and 0.2 selects 20% of the training set as validation data.
///    - 100 selects 100 events as the validation data.
///
/// @return number of samples in validation set
///

UInt_t TMVA::MethodDNN::GetNumValidationSamples()
{
   Int_t nValidationSamples = 0;
   UInt_t trainingSetSize = GetEventCollection(Types::kTraining).size();

   // Parsing + Validation
   // --------------------
   if (fValidationSize.EndsWith("%")) {
      // Relative spec. format 20%
      TString intValStr = TString(fValidationSize.Strip(TString::kTrailing, '%'));

      if (intValStr.IsFloat()) {
         Double_t valSizeAsDouble = fValidationSize.Atof() / 100.0;
         nValidationSamples = GetEventCollection(Types::kTraining).size() * valSizeAsDouble;
      } else {
         Log() << kFATAL << "Cannot parse number \"" << fValidationSize
               << "\". Expected string like \"20%\" or \"20.0%\"." << Endl;
      }
   } else if (fValidationSize.IsFloat()) {
      Double_t valSizeAsDouble = fValidationSize.Atof();

      if (valSizeAsDouble < 1.0) {
         // Relative spec. format 0.2
         nValidationSamples = GetEventCollection(Types::kTraining).size() * valSizeAsDouble;
      } else {
         // Absolute spec format 100 or 100.0
         nValidationSamples = valSizeAsDouble;
      }
   } else {
      Log() << kFATAL << "Cannot parse number \"" << fValidationSize << "\". Expected string like \"0.2\" or \"100\"."
            << Endl;
   }

   // Value validation
   // ----------------
   if (nValidationSamples < 0) {
      Log() << kFATAL << "Validation size \"" << fValidationSize << "\" is negative." << Endl;
   }

   if (nValidationSamples == 0) {
      Log() << kFATAL << "Validation size \"" << fValidationSize << "\" is zero." << Endl;
   }

   if (nValidationSamples >= (Int_t)trainingSetSize) {
      Log() << kFATAL << "Validation size \"" << fValidationSize
            << "\" is larger than or equal in size to training set (size=\"" << trainingSetSize << "\")." << Endl;
   }

   return nValidationSamples;
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::MethodDNN::Train()
{
   if (fInteractive && fInteractive->NotInitialized()){
      std::vector<TString> titles = {"Error on training set", "Error on test set"};
      fInteractive->Init(titles);
      // JsMVA progress bar maximum (100%)
      fIPyMaxIter = 100;
   }

   for (TTrainingSettings & settings : fTrainingSettings) {
      size_t nValidationSamples = GetNumValidationSamples();
      size_t nTrainingSamples = GetEventCollection(Types::kTraining).size() - nValidationSamples;
      size_t nTestSamples = nValidationSamples;

      if (nTrainingSamples < settings.batchSize ||
          nValidationSamples < settings.batchSize ||
          nTestSamples < settings.batchSize) {
         Log() << kFATAL << "Number of samples in the datasets are train: "
                         << nTrainingSamples << " valid: " << nValidationSamples
                         << " test: " << nTestSamples << ". "
                         << "One of these is smaller than the batch size of "
                         << settings.batchSize << ". Please increase the batch"
                         << " size to be at least the same size as the smallest"
                         << " of these values." << Endl;
      }
  }

   if (fArchitectureString == "GPU") {
       TrainGpu();
       if (!fExitFromTraining) fIPyMaxIter = fIPyCurrentIter;
       ExitFromTraining();
       return;
   } else if (fArchitectureString == "OpenCL") {
      Log() << kFATAL << "OpenCL backend not yet supported." << Endl;
      return;
   } else if (fArchitectureString == "CPU") {
      TrainCpu();
      if (!fExitFromTraining) fIPyMaxIter = fIPyCurrentIter;
      ExitFromTraining();
      return;
   }

   Log() << kINFO << "Using Standard Implementation.";

   std::vector<Pattern> trainPattern;
   std::vector<Pattern> testPattern;

   size_t nValidationSamples = GetNumValidationSamples();
   size_t nTrainingSamples = GetEventCollection(Types::kTraining).size() - nValidationSamples;

   const std::vector<TMVA::Event *> &allData = GetEventCollection(Types::kTraining);
   const std::vector<TMVA::Event *> eventCollectionTraining{allData.begin(), allData.begin() + nTrainingSamples};
   const std::vector<TMVA::Event *> eventCollectionTesting{allData.begin() + nTrainingSamples, allData.end()};

   for (auto &event : eventCollectionTraining) {
      const std::vector<Float_t>& values = event->GetValues();
      if (fAnalysisType == Types::kClassification) {
         double outputValue = event->GetClass () == 0 ? 0.9 : 0.1;
         trainPattern.push_back(Pattern (values.begin(),
                                         values.end(),
                                         outputValue,
                                         event->GetWeight()));
         trainPattern.back().addInput(1.0);
      } else if (fAnalysisType == Types::kMulticlass) {
         std::vector<Float_t> oneHot(DataInfo().GetNClasses(), 0.0);
         oneHot[event->GetClass()] = 1.0;
         trainPattern.push_back(Pattern (values.begin(), values.end(),
                                        oneHot.cbegin(), oneHot.cend(),
                                        event->GetWeight()));
         trainPattern.back().addInput(1.0);
      } else {
         const std::vector<Float_t>& targets = event->GetTargets ();
         trainPattern.push_back(Pattern(values.begin(),
                                        values.end(),
                                        targets.begin(),
                                        targets.end(),
                                        event->GetWeight ()));
         trainPattern.back ().addInput (1.0); // bias node
      }
   }

   for (auto &event : eventCollectionTesting) {
      const std::vector<Float_t>& values = event->GetValues();
      if (fAnalysisType == Types::kClassification) {
         double outputValue = event->GetClass () == 0 ? 0.9 : 0.1;
         testPattern.push_back(Pattern (values.begin(),
                                         values.end(),
                                         outputValue,
                                         event->GetWeight()));
         testPattern.back().addInput(1.0);
      } else if (fAnalysisType == Types::kMulticlass) {
         std::vector<Float_t> oneHot(DataInfo().GetNClasses(), 0.0);
         oneHot[event->GetClass()] = 1.0;
         testPattern.push_back(Pattern (values.begin(), values.end(),
                                        oneHot.cbegin(), oneHot.cend(),
                                        event->GetWeight()));
         testPattern.back().addInput(1.0);
      } else {
         const std::vector<Float_t>& targets = event->GetTargets ();
         testPattern.push_back(Pattern(values.begin(),
                                        values.end(),
                                        targets.begin(),
                                        targets.end(),
                                        event->GetWeight ()));
         testPattern.back ().addInput (1.0); // bias node
      }
   }

   TMVA::DNN::Net      net;
   std::vector<double> weights;

   net.SetIpythonInteractive(fInteractive, &fExitFromTraining, &fIPyMaxIter, &fIPyCurrentIter);

   net.setInputSize(fNet.GetInputWidth() + 1);
   net.setOutputSize(fNet.GetOutputWidth() + 1);

   for (size_t i = 0; i < fNet.GetDepth(); i++) {
      EActivationFunction f = fNet.GetLayer(i).GetActivationFunction();
      EnumFunction        g = EnumFunction::LINEAR;
      switch(f) {
         case EActivationFunction::kIdentity: g = EnumFunction::LINEAR;   break;
         case EActivationFunction::kRelu:     g = EnumFunction::RELU;     break;
         case EActivationFunction::kSigmoid:  g = EnumFunction::SIGMOID;  break;
         case EActivationFunction::kTanh:     g = EnumFunction::TANH;     break;
         case EActivationFunction::kFastTanh: g = EnumFunction::TANH;     break;
         case EActivationFunction::kSymmRelu: g = EnumFunction::SYMMRELU; break;
         case EActivationFunction::kSoftSign: g = EnumFunction::SOFTSIGN; break;
         case EActivationFunction::kGauss:    g = EnumFunction::GAUSS;    break;
      }
      if (i < fNet.GetDepth() - 1) {
         net.addLayer(Layer(fNet.GetLayer(i).GetWidth(), g));
      } else {
         ModeOutputValues h = ModeOutputValues::DIRECT;
         switch(fOutputFunction) {
            case EOutputFunction::kIdentity: h = ModeOutputValues::DIRECT;  break;
            case EOutputFunction::kSigmoid:  h = ModeOutputValues::SIGMOID; break;
            case EOutputFunction::kSoftmax:  h = ModeOutputValues::SOFTMAX; break;
         }
         net.addLayer(Layer(fNet.GetLayer(i).GetWidth(), g, h));
      }
   }

   switch(fNet.GetLossFunction()) {
      case ELossFunction::kMeanSquaredError:
         net.setErrorFunction(ModeErrorFunction::SUMOFSQUARES);
         break;
      case ELossFunction::kCrossEntropy:
         net.setErrorFunction(ModeErrorFunction::CROSSENTROPY);
         break;
      case ELossFunction::kSoftmaxCrossEntropy:
         net.setErrorFunction(ModeErrorFunction::CROSSENTROPY_MUTUALEXCLUSIVE);
         break;
   }

   switch(fWeightInitialization) {
      case EInitialization::kGauss:
          net.initializeWeights(WeightInitializationStrategy::XAVIER,
                                std::back_inserter(weights));
          break;
      case EInitialization::kUniform:
          net.initializeWeights(WeightInitializationStrategy::XAVIERUNIFORM,
                                std::back_inserter(weights));
          break;
      default:
          net.initializeWeights(WeightInitializationStrategy::XAVIER,
                                std::back_inserter(weights));
          break;
   }

   int idxSetting = 0;
   for (auto s : fTrainingSettings) {

      EnumRegularization r = EnumRegularization::NONE;
      switch(s.regularization) {
         case ERegularization::kNone: r = EnumRegularization::NONE; break;
         case ERegularization::kL1:   r = EnumRegularization::L1;   break;
         case ERegularization::kL2:   r = EnumRegularization::L2;   break;
      }

      Settings * settings = new Settings(TString(), s.convergenceSteps, s.batchSize,
                                         s.testInterval, s.weightDecay, r,
                                         MinimizerType::fSteepest, s.learningRate,
                                         s.momentum, 1, s.multithreading);
      std::shared_ptr<Settings> ptrSettings(settings);
      ptrSettings->setMonitoring (0);
      Log() << kINFO
            << "Training with learning rate = " << ptrSettings->learningRate ()
            << ", momentum = " << ptrSettings->momentum ()
            << ", repetitions = " << ptrSettings->repetitions ()
            << Endl;

      ptrSettings->setProgressLimits ((idxSetting)*100.0/(fSettings.size ()),
                                      (idxSetting+1)*100.0/(fSettings.size ()));

      const std::vector<double>& dropConfig = ptrSettings->dropFractions ();
      if (!dropConfig.empty ()) {
         Log () << kINFO << "Drop configuration" << Endl
                << "    drop repetitions = " << ptrSettings->dropRepetitions()
                << Endl;
      }

      int idx = 0;
      for (auto f : dropConfig) {
         Log () << kINFO << "    Layer " << idx << " = " << f << Endl;
         ++idx;
      }
      Log () << kINFO << Endl;

      DNN::Steepest minimizer(ptrSettings->learningRate(),
                              ptrSettings->momentum(),
                              ptrSettings->repetitions());
      net.train(weights, trainPattern, testPattern, minimizer, *ptrSettings.get());
      ptrSettings.reset();
      Log () << kINFO << Endl;
      idxSetting++;
   }
   size_t weightIndex = 0;
   for (size_t l = 0; l < fNet.GetDepth(); l++) {
      auto & layerWeights = fNet.GetLayer(l).GetWeights();
      for (Int_t j = 0; j < layerWeights.GetNcols(); j++) {
         for (Int_t i = 0; i < layerWeights.GetNrows(); i++) {
            layerWeights(i,j) = weights[weightIndex];
            weightIndex++;
         }
      }
      auto & layerBiases = fNet.GetLayer(l).GetBiases();
      if (l == 0) {
         for (Int_t i = 0; i < layerBiases.GetNrows(); i++) {
            layerBiases(i,0) = weights[weightIndex];
            weightIndex++;
         }
      } else {
         for (Int_t i = 0; i < layerBiases.GetNrows(); i++) {
            layerBiases(i,0) = 0.0;
         }
      }
   }
   if (!fExitFromTraining) fIPyMaxIter = fIPyCurrentIter;
   ExitFromTraining();
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::MethodDNN::TrainGpu()
{

#ifdef DNNCUDA // Included only if DNNCUDA flag is set.
   Log() << kINFO << "Start of neural network training on GPU." << Endl << Endl;

   size_t nValidationSamples = GetNumValidationSamples();
   size_t nTrainingSamples = GetEventCollection(Types::kTraining).size() - nValidationSamples;
   size_t nTestSamples = nValidationSamples;

   Log() << kDEBUG << "Using " << nValidationSamples << " validation samples." << Endl;
   Log() << kDEBUG << "Using " << nTestSamples << " training samples." << Endl;

   size_t trainingPhase = 1;
   fNet.Initialize(fWeightInitialization);
   for (TTrainingSettings & settings : fTrainingSettings) {

      if (fInteractive){
         fInteractive->ClearGraphs();
      }

      TNet<TCuda<>> net(settings.batchSize, fNet);
      net.SetWeightDecay(settings.weightDecay);
      net.SetRegularization(settings.regularization);

      // Need to convert dropoutprobabilities to conventions used
      // by backend implementation.
      std::vector<Double_t> dropoutVector(settings.dropoutProbabilities);
      for (auto & p : dropoutVector) {
         p = 1.0 - p;
      }
      net.SetDropoutProbabilities(dropoutVector);

      net.InitializeGradients();
      auto testNet = net.CreateClone(settings.batchSize);

      Log() << kINFO << "Training phase " << trainingPhase << " of "
            << fTrainingSettings.size() << ":" << Endl;
      trainingPhase++;

      using DataLoader_t = TDataLoader<TMVAInput_t, TCuda<>>;

      // Split training data into training and validation set
      const std::vector<Event *> &allData = GetEventCollection(Types::kTraining);
      const std::vector<Event *> trainingInputData =
         std::vector<Event *>(allData.begin(), allData.begin() + nTrainingSamples);
      const std::vector<Event *> testInputData =
         std::vector<Event *>(allData.begin() + nTrainingSamples, allData.end());

      if (trainingInputData.size() != nTrainingSamples) {
         Log() << kFATAL << "Inconsistent training sample size" << Endl;
      }
      if (testInputData.size() != nTestSamples) {
         Log() << kFATAL << "Inconsistent test sample size" << Endl;
      }

      size_t nThreads = 1;
      TMVAInput_t trainingTuple = std::tie(trainingInputData, DataInfo());
      TMVAInput_t testTuple = std::tie(testInputData, DataInfo());
      DataLoader_t trainingData(trainingTuple, nTrainingSamples,
                                net.GetBatchSize(), net.GetInputWidth(),
                                net.GetOutputWidth(), nThreads);
      DataLoader_t testData(testTuple, nTestSamples, testNet.GetBatchSize(),
                            net.GetInputWidth(), net.GetOutputWidth(),
                            nThreads);
      DNN::TGradientDescent<TCuda<>> minimizer(settings.learningRate,
                                             settings.convergenceSteps,
                                             settings.testInterval);

      std::vector<TNet<TCuda<>>> nets{};
      std::vector<TBatch<TCuda<>>> batches{};
      nets.reserve(nThreads);
      for (size_t i = 0; i < nThreads; i++) {
         nets.push_back(net);
         for (size_t j = 0; j < net.GetDepth(); j++)
         {
            auto &masterLayer = net.GetLayer(j);
            auto &layer = nets.back().GetLayer(j);
            TCuda<>::Copy(layer.GetWeights(),
                          masterLayer.GetWeights());
            TCuda<>::Copy(layer.GetBiases(),
                          masterLayer.GetBiases());
         }
      }

      bool   converged = false;
      size_t stepCount = 0;
      size_t batchesInEpoch = nTrainingSamples / net.GetBatchSize();

      std::chrono::time_point<std::chrono::system_clock> start, end;
      start = std::chrono::system_clock::now();

      if (!fInteractive) {
         Log() << std::setw(10) << "Epoch" << " | "
               << std::setw(12) << "Train Err."
               << std::setw(12) << "Test  Err."
               << std::setw(12) << "GFLOP/s"
               << std::setw(12) << "Conv. Steps" << Endl;
         std::string separator(62, '-');
         Log() << separator << Endl;
      }

      while (!converged)
      {
         stepCount++;

         // Perform minimization steps for a full epoch.
         trainingData.Shuffle();
         for (size_t i = 0; i < batchesInEpoch; i += nThreads) {
             batches.clear();
             for (size_t j = 0; j < nThreads; j++) {
                 batches.reserve(nThreads);
                 batches.push_back(trainingData.GetBatch());
             }
             if (settings.momentum > 0.0) {
                 minimizer.StepMomentum(net, nets, batches, settings.momentum);
             } else {
                 minimizer.Step(net, nets, batches);
             }
         }

         if ((stepCount % minimizer.GetTestInterval()) == 0) {

            // Compute test error.
            Double_t testError = 0.0;
            for (auto batch : testData) {
               auto inputMatrix  = batch.GetInput();
               auto outputMatrix = batch.GetOutput();
               testError += testNet.Loss(inputMatrix, outputMatrix);
            }
            testError /= (Double_t) (nTestSamples / settings.batchSize);

            //Log the loss value
            fTrainHistory.AddValue("testError",stepCount,testError);

            end   = std::chrono::system_clock::now();

            // Compute training error.
            Double_t trainingError = 0.0;
            for (auto batch : trainingData) {
               auto inputMatrix  = batch.GetInput();
               auto outputMatrix = batch.GetOutput();
               trainingError += net.Loss(inputMatrix, outputMatrix);
            }
            trainingError /= (Double_t) (nTrainingSamples / settings.batchSize);
            //Log the loss value
            fTrainHistory.AddValue("trainingError",stepCount,trainingError);

            // Compute numerical throughput.
            std::chrono::duration<double> elapsed_seconds = end - start;
            double seconds = elapsed_seconds.count();
            double nFlops  = (double) (settings.testInterval * batchesInEpoch);
            nFlops *= net.GetNFlops() * 1e-9;

            converged = minimizer.HasConverged(testError);
            start = std::chrono::system_clock::now();

            if (fInteractive) {
               fInteractive->AddPoint(stepCount, trainingError, testError);
               fIPyCurrentIter = 100.0 * minimizer.GetConvergenceCount()
                                  / minimizer.GetConvergenceSteps ();
               if (fExitFromTraining) break;
            } else {
               Log() << std::setw(10) << stepCount << " | "
                     << std::setw(12) << trainingError
                     << std::setw(12) << testError
                     << std::setw(12) << nFlops / seconds
                     << std::setw(12) << minimizer.GetConvergenceCount() << Endl;
               if (converged) {
                  Log() << Endl;
               }
            }
         }
      }
      for (size_t l = 0; l < net.GetDepth(); l++) {
         fNet.GetLayer(l).GetWeights() = (TMatrixT<Scalar_t>) net.GetLayer(l).GetWeights();
         fNet.GetLayer(l).GetBiases()  = (TMatrixT<Scalar_t>) net.GetLayer(l).GetBiases();
      }
   }

#else // DNNCUDA flag not set.

   Log() << kFATAL << "CUDA backend not enabled. Please make sure "
                      "you have CUDA installed and it was successfully "
                      "detected by CMAKE." << Endl;
#endif // DNNCUDA
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::MethodDNN::TrainCpu()
{

#ifdef DNNCPU // Included only if DNNCPU flag is set.
   Log() << kINFO << "Start of neural network training on CPU." << Endl << Endl;

   size_t nValidationSamples = GetNumValidationSamples();
   size_t nTrainingSamples = GetEventCollection(Types::kTraining).size() - nValidationSamples;
   size_t nTestSamples = nValidationSamples;

   Log() << kDEBUG << "Using " << nValidationSamples << " validation samples." << Endl;
   Log() << kDEBUG << "Using " << nTestSamples << " training samples." << Endl;

   fNet.Initialize(fWeightInitialization);

   size_t trainingPhase = 1;
   for (TTrainingSettings & settings : fTrainingSettings) {

      if (fInteractive){
         fInteractive->ClearGraphs();
      }

      Log() << "Training phase " << trainingPhase << " of "
            << fTrainingSettings.size() << ":" << Endl;
      trainingPhase++;

      TNet<TCpu<>> net(settings.batchSize, fNet);
      net.SetWeightDecay(settings.weightDecay);
      net.SetRegularization(settings.regularization);
      // Need to convert dropoutprobabilities to conventions used
      // by backend implementation.
      std::vector<Double_t> dropoutVector(settings.dropoutProbabilities);
      for (auto & p : dropoutVector) {
         p = 1.0 - p;
      }
      net.SetDropoutProbabilities(dropoutVector);
      net.InitializeGradients();
      auto testNet = net.CreateClone(settings.batchSize);

      using DataLoader_t = TDataLoader<TMVAInput_t, TCpu<>>;

      // Split training data into training and validation set
      const std::vector<Event *> &allData = GetEventCollection(Types::kTraining);
      const std::vector<Event *> trainingInputData =
         std::vector<Event *>(allData.begin(), allData.begin() + nTrainingSamples);
      const std::vector<Event *> testInputData =
         std::vector<Event *>(allData.begin() + nTrainingSamples, allData.end());

      if (trainingInputData.size() != nTrainingSamples) {
         Log() << kFATAL << "Inconsistent training sample size" << Endl;
      }
      if (testInputData.size() != nTestSamples) {
         Log() << kFATAL << "Inconsistent test sample size" << Endl;
      }

      size_t nThreads = 1;
      TMVAInput_t trainingTuple = std::tie(trainingInputData, DataInfo());
      TMVAInput_t testTuple = std::tie(testInputData, DataInfo());
      DataLoader_t trainingData(trainingTuple, nTrainingSamples,
                                net.GetBatchSize(), net.GetInputWidth(),
                                net.GetOutputWidth(), nThreads);
      DataLoader_t testData(testTuple, nTestSamples, testNet.GetBatchSize(),
                            net.GetInputWidth(), net.GetOutputWidth(),
                            nThreads);
      DNN::TGradientDescent<TCpu<>> minimizer(settings.learningRate,
                                               settings.convergenceSteps,
                                               settings.testInterval);

      std::vector<TNet<TCpu<>>>   nets{};
      std::vector<TBatch<TCpu<>>> batches{};
      nets.reserve(nThreads);
      for (size_t i = 0; i < nThreads; i++) {
         nets.push_back(net);
         for (size_t j = 0; j < net.GetDepth(); j++)
         {
            auto &masterLayer = net.GetLayer(j);
            auto &layer = nets.back().GetLayer(j);
            TCpu<>::Copy(layer.GetWeights(),
                          masterLayer.GetWeights());
            TCpu<>::Copy(layer.GetBiases(),
                          masterLayer.GetBiases());
         }
      }

      bool   converged = false;
      size_t stepCount = 0;
      size_t batchesInEpoch = nTrainingSamples / net.GetBatchSize();

      std::chrono::time_point<std::chrono::system_clock> start, end;
      start = std::chrono::system_clock::now();

      if (!fInteractive) {
         Log() << std::setw(10) << "Epoch" << " | "
               << std::setw(12) << "Train Err."
               << std::setw(12) << "Test  Err."
               << std::setw(12) << "GFLOP/s"
               << std::setw(12) << "Conv. Steps" << Endl;
         std::string separator(62, '-');
         Log() << separator << Endl;
      }

      while (!converged)
      {
         stepCount++;
         // Perform minimization steps for a full epoch.
         trainingData.Shuffle();
         for (size_t i = 0; i < batchesInEpoch; i += nThreads) {
             batches.clear();
             for (size_t j = 0; j < nThreads; j++) {
                 batches.reserve(nThreads);
                 batches.push_back(trainingData.GetBatch());
             }
             if (settings.momentum > 0.0) {
                 minimizer.StepMomentum(net, nets, batches, settings.momentum);
             } else {
                 minimizer.Step(net, nets, batches);
             }
         }

         if ((stepCount % minimizer.GetTestInterval()) == 0) {

            // Compute test error.
            Double_t testError = 0.0;
            for (auto batch : testData) {
               auto inputMatrix  = batch.GetInput();
               auto outputMatrix = batch.GetOutput();
               auto weightMatrix = batch.GetWeights();
               testError += testNet.Loss(inputMatrix, outputMatrix, weightMatrix);
            }
            testError /= (Double_t) (nTestSamples / settings.batchSize);

            //Log the loss value
            fTrainHistory.AddValue("testError",stepCount,testError);

            end   = std::chrono::system_clock::now();

            // Compute training error.
            Double_t trainingError = 0.0;
            for (auto batch : trainingData) {
               auto inputMatrix  = batch.GetInput();
               auto outputMatrix = batch.GetOutput();
               auto weightMatrix = batch.GetWeights();
               trainingError += net.Loss(inputMatrix, outputMatrix, weightMatrix);
            }
            trainingError /= (Double_t) (nTrainingSamples / settings.batchSize);

            //Log the loss value
            fTrainHistory.AddValue("trainingError",stepCount,trainingError);

            if (fInteractive){
               fInteractive->AddPoint(stepCount, trainingError, testError);
               fIPyCurrentIter = 100*(double)minimizer.GetConvergenceCount() /(double)settings.convergenceSteps;
               if (fExitFromTraining) break;
            }

            // Compute numerical throughput.
            std::chrono::duration<double> elapsed_seconds = end - start;
            double seconds = elapsed_seconds.count();
            double nFlops  = (double) (settings.testInterval * batchesInEpoch);
            nFlops *= net.GetNFlops() * 1e-9;

            converged = minimizer.HasConverged(testError);
            start = std::chrono::system_clock::now();

            if (fInteractive) {
               fInteractive->AddPoint(stepCount, trainingError, testError);
               fIPyCurrentIter = 100.0 * minimizer.GetConvergenceCount()
                                  / minimizer.GetConvergenceSteps ();
               if (fExitFromTraining) break;
            } else {
               Log() << std::setw(10) << stepCount << " | "
                     << std::setw(12) << trainingError
                     << std::setw(12) << testError
                     << std::setw(12) << nFlops / seconds
                     << std::setw(12) << minimizer.GetConvergenceCount() << Endl;
               if (converged) {
                  Log() << Endl;
               }
            }
         }
      }


      for (size_t l = 0; l < net.GetDepth(); l++) {
         auto & layer = fNet.GetLayer(l);
         layer.GetWeights() = (TMatrixT<Scalar_t>) net.GetLayer(l).GetWeights();
         layer.GetBiases()  = (TMatrixT<Scalar_t>) net.GetLayer(l).GetBiases();
      }
   }

#else // DNNCPU flag not set.
   Log() << kFATAL << "Multi-core CPU backend not enabled. Please make sure "
                      "you have a BLAS implementation and it was successfully "
                      "detected by CMake as well that the imt CMake flag is set." << Endl;
#endif // DNNCPU
}

////////////////////////////////////////////////////////////////////////////////

Double_t TMVA::MethodDNN::GetMvaValue( Double_t* /*errLower*/, Double_t* /*errUpper*/ )
{
   size_t nVariables = GetEvent()->GetNVariables();
   Matrix_t X(1, nVariables);
   Matrix_t YHat(1, 1);

   const std::vector<Float_t>& inputValues = GetEvent()->GetValues();
   for (size_t i = 0; i < nVariables; i++) {
      X(0,i) = inputValues[i];
   }

   fNet.Prediction(YHat, X, fOutputFunction);
   return YHat(0,0);
}

////////////////////////////////////////////////////////////////////////////////

const std::vector<Float_t> & TMVA::MethodDNN::GetRegressionValues()
{
   size_t nVariables = GetEvent()->GetNVariables();
   Matrix_t X(1, nVariables);

   const Event *ev = GetEvent();
   const std::vector<Float_t>& inputValues = ev->GetValues();
   for (size_t i = 0; i < nVariables; i++) {
       X(0,i) = inputValues[i];
   }

   size_t nTargets = std::max(1u, ev->GetNTargets());
   Matrix_t YHat(1, nTargets);
   std::vector<Float_t> output(nTargets);
   auto net = fNet.CreateClone(1);
   net.Prediction(YHat, X, fOutputFunction);

   for (size_t i = 0; i < nTargets; i++)
       output[i] = YHat(0, i);

   if (fRegressionReturnVal == NULL) {
       fRegressionReturnVal = new std::vector<Float_t>();
   }
   fRegressionReturnVal->clear();

   Event * evT = new Event(*ev);
   for (size_t i = 0; i < nTargets; ++i) {
      evT->SetTarget(i, output[i]);
   }

   const Event* evT2 = GetTransformationHandler().InverseTransform(evT);
   for (size_t i = 0; i < nTargets; ++i) {
      fRegressionReturnVal->push_back(evT2->GetTarget(i));
   }
   delete evT;
   return *fRegressionReturnVal;
}

const std::vector<Float_t> & TMVA::MethodDNN::GetMulticlassValues()
{
   size_t nVariables = GetEvent()->GetNVariables();
   Matrix_t X(1, nVariables);
   Matrix_t YHat(1, DataInfo().GetNClasses());
   if (fMulticlassReturnVal == NULL) {
      fMulticlassReturnVal = new std::vector<Float_t>(DataInfo().GetNClasses());
   }

   const std::vector<Float_t>& inputValues = GetEvent()->GetValues();
   for (size_t i = 0; i < nVariables; i++) {
      X(0,i) = inputValues[i];
   }

   fNet.Prediction(YHat, X, fOutputFunction);
   for (size_t i = 0; i < (size_t) YHat.GetNcols(); i++) {
      (*fMulticlassReturnVal)[i] = YHat(0, i);
   }
   return *fMulticlassReturnVal;
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::MethodDNN::AddWeightsXMLTo( void* parent ) const
{
   void* nn = gTools().xmlengine().NewChild(parent, 0, "Weights");
   Int_t inputWidth = fNet.GetInputWidth();
   Int_t depth      = fNet.GetDepth();
   char  lossFunction = static_cast<char>(fNet.GetLossFunction());
   gTools().xmlengine().NewAttr(nn, 0, "InputWidth",
                                gTools().StringFromInt(inputWidth));
   gTools().xmlengine().NewAttr(nn, 0, "Depth", gTools().StringFromInt(depth));
   gTools().xmlengine().NewAttr(nn, 0, "LossFunction", TString(lossFunction));
   gTools().xmlengine().NewAttr(nn, 0, "OutputFunction",
                                TString(static_cast<char>(fOutputFunction)));

   for (Int_t i = 0; i < depth; i++) {
      const auto& layer = fNet.GetLayer(i);
      auto layerxml = gTools().xmlengine().NewChild(nn, 0, "Layer");
      int activationFunction = static_cast<int>(layer.GetActivationFunction());
      gTools().xmlengine().NewAttr(layerxml, 0, "ActivationFunction",
                                   TString::Itoa(activationFunction, 10));
      WriteMatrixXML(layerxml, "Weights", layer.GetWeights());
      WriteMatrixXML(layerxml, "Biases",  layer.GetBiases());
   }
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::MethodDNN::ReadWeightsFromXML(void* rootXML)
{
   auto netXML = gTools().GetChild(rootXML, "Weights");
   if (!netXML){
      netXML = rootXML;
   }

   fNet.Clear();
   fNet.SetBatchSize(1);

   size_t inputWidth, depth;
   gTools().ReadAttr(netXML, "InputWidth", inputWidth);
   gTools().ReadAttr(netXML, "Depth", depth);
   char lossFunctionChar;
   gTools().ReadAttr(netXML, "LossFunction", lossFunctionChar);
   char outputFunctionChar;
   gTools().ReadAttr(netXML, "OutputFunction", outputFunctionChar);

   fNet.SetInputWidth(inputWidth);
   fNet.SetLossFunction(static_cast<ELossFunction>(lossFunctionChar));
   fOutputFunction = static_cast<EOutputFunction>(outputFunctionChar);

   size_t previousWidth = inputWidth;
   auto layerXML = gTools().xmlengine().GetChild(netXML, "Layer");
   for (size_t i = 0; i < depth; i++) {
      TString fString;
      EActivationFunction f;

      // Read activation function.
      gTools().ReadAttr(layerXML, "ActivationFunction", fString);
      f = static_cast<EActivationFunction>(fString.Atoi());

      // Read number of neurons.
      size_t width;
      auto matrixXML = gTools().GetChild(layerXML, "Weights");
      gTools().ReadAttr(matrixXML, "rows", width);

      fNet.AddLayer(width, f);
      TMatrixT<Double_t> weights(width, previousWidth);
      TMatrixT<Double_t> biases(width, 1);
      ReadMatrixXML(layerXML, "Weights", weights);
      ReadMatrixXML(layerXML, "Biases",  biases);
      fNet.GetLayer(i).GetWeights() = weights;
      fNet.GetLayer(i).GetBiases()  = biases;

      layerXML = gTools().GetNextChild(layerXML);
      previousWidth = width;
   }
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::MethodDNN::ReadWeightsFromStream( std::istream & /*istr*/)
{
}

////////////////////////////////////////////////////////////////////////////////

const TMVA::Ranking* TMVA::MethodDNN::CreateRanking()
{
   fRanking = new Ranking( GetName(), "Importance" );
   for (UInt_t ivar=0; ivar<GetNvar(); ivar++) {
      fRanking->AddRank( Rank( GetInputLabel(ivar), 1.0));
   }
   return fRanking;
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::MethodDNN::MakeClassSpecific( std::ostream& /*fout*/,
                                         const TString& /*className*/ ) const
{
}

////////////////////////////////////////////////////////////////////////////////

void TMVA::MethodDNN::GetHelpMessage() const
{
   // get help message text
   //
   // typical length of text line:
   //         "|--------------------------------------------------------------|"
   TString col    = gConfig().WriteOptionsReference() ? TString() : gTools().Color("bold");
   TString colres = gConfig().WriteOptionsReference() ? TString() : gTools().Color("reset");

   Log() << Endl;
   Log() << col << "--- Short description:" << colres << Endl;
   Log() << Endl;
   Log() << "The DNN neural network is a feedforward" << Endl;
   Log() << "multilayer perceptron implementation. The DNN has a user-" << Endl;
   Log() << "defined hidden layer architecture, where the number of input (output)" << Endl;
   Log() << "nodes is determined by the input variables (output classes, i.e., " << Endl;
   Log() << "signal and one background, regression or multiclass). " << Endl;
   Log() << Endl;
   Log() << col << "--- Performance optimisation:" << colres << Endl;
   Log() << Endl;

   const char* txt = "The DNN supports various options to improve performance in terms of training speed and \n \
reduction of overfitting: \n \
\n \
      - different training settings can be stacked. Such that the initial training  \n\
        is done with a large learning rate and a large drop out fraction whilst \n \
        in a later stage learning rate and drop out can be reduced. \n \
      - drop out  \n \
        [recommended: \n \
         initial training stage: 0.0 for the first layer, 0.5 for later layers. \n \
         later training stage: 0.1 or 0.0 for all layers \n \
         final training stage: 0.0] \n \
        Drop out is a technique where a at each training cycle a fraction of arbitrary  \n \
        nodes is disabled. This reduces co-adaptation of weights and thus reduces overfitting. \n \
      - L1 and L2 regularization are available \n \
      - Minibatches  \n \
        [recommended 10 - 150] \n \
        Arbitrary mini-batch sizes can be chosen. \n \
      - Multithreading \n \
        [recommended: True] \n \
        Multithreading can be turned on. The minibatches are distributed to the available \n \
        cores. The algorithm is lock-free (\"Hogwild!\"-style) for each cycle. \n \
 \n \
      Options: \n \
      \"Layout\": \n \
          - example: \"TANH|(N+30)*2,TANH|(N+30),LINEAR\" \n \
          - meaning:  \n \
              . two hidden layers (separated by \",\") \n \
              . the activation function is TANH (other options: RELU, SOFTSIGN, LINEAR) \n \
              . the activation function for the output layer is LINEAR \n \
              . the first hidden layer has (N+30)*2 nodes where N is the number of input neurons \n \
              . the second hidden layer has N+30 nodes, where N is the number of input neurons \n \
              . the number of nodes in the output layer is determined by the number of output nodes \n \
                and can therefore not be chosen freely.  \n \
 \n \
       \"ErrorStrategy\": \n \
           - SUMOFSQUARES \n \
             The error of the neural net is determined by a sum-of-squares error function \n \
             For regression, this is the only possible choice.  \n \
           - CROSSENTROPY \n \
             The error of the neural net is determined by a cross entropy function. The \n \
             output values are automatically (internally) transformed into probabilities \n \
             using a sigmoid function. \n \
             For signal/background classification this is the default choice.  \n \
             For multiclass using cross entropy more than one or no output classes  \n \
             can be equally true or false (e.g. Event 0: A and B are true, Event 1:  \n \
             A and C is true, Event 2: C is true, ...) \n \
           - MUTUALEXCLUSIVE \n \
             In multiclass settings, exactly one of the output classes can be true (e.g. either A or B or C) \n \
 \n \
        \"WeightInitialization\" \n \
           - XAVIER \n \
             [recommended] \n \
             \"Xavier Glorot & Yoshua Bengio\"-style of initializing the weights. The weights are chosen randomly \n \
             such that the variance of the values of the nodes is preserved for each layer.  \n \
           - XAVIERUNIFORM \n \
             The same as XAVIER, but with uniformly distributed weights instead of gaussian weights \n \
           - LAYERSIZE \n \
             Random values scaled by the layer size \n \
 \n \
         \"TrainingStrategy\" \n \
           - example: \"LearningRate=1e-1,Momentum=0.3,ConvergenceSteps=50,BatchSize=30,TestRepetitions=7,WeightDecay=0.0,Renormalize=L2,DropConfig=0.0,DropRepetitions=5|LearningRate=1e-4,Momentum=0.3,ConvergenceSteps=50,BatchSize=20,TestRepetitions=7,WeightDecay=0.001,Renormalize=L2,DropFraction=0.0,DropRepetitions=5\" \n \
           - explanation: two stacked training settings separated by \"|\" \n \
             . first training setting: \"LearningRate=1e-1,Momentum=0.3,ConvergenceSteps=50,BatchSize=30,TestRepetitions=7,WeightDecay=0.0,Renormalize=L2,DropConfig=0.0,DropRepetitions=5\" \n \
             . second training setting : \"LearningRate=1e-4,Momentum=0.3,ConvergenceSteps=50,BatchSize=20,TestRepetitions=7,WeightDecay=0.001,Renormalize=L2,DropFractions=0.0,DropRepetitions=5\" \n \
             . LearningRate :  \n \
               - recommended for classification: 0.1 initially, 1e-4 later \n \
               - recommended for regression: 1e-4 and less \n \
             . Momentum : \n \
               preserve a fraction of the momentum for the next training batch [fraction = 0.0 - 1.0] \n \
             . Repetitions : \n \
               train \"Repetitions\" repetitions with the same minibatch before switching to the next one \n \
             . ConvergenceSteps :  \n \
               Assume that convergence is reached after \"ConvergenceSteps\" cycles where no improvement \n \
               of the error on the test samples has been found. (Mind that only at each \"TestRepetitions\"  \n \
               cycle the test samples are evaluated and thus the convergence is checked) \n \
             . BatchSize \n \
               Size of the mini-batches.  \n \
             . TestRepetitions \n \
               Perform testing the neural net on the test samples each \"TestRepetitions\" cycle \n \
             . WeightDecay \n \
               If \"Renormalize\" is set to L1 or L2, \"WeightDecay\" provides the renormalization factor \n \
             . Renormalize \n \
               NONE, L1 (|w|) or L2 (w^2) \n \
             . DropConfig \n \
               Drop a fraction of arbitrary nodes of each of the layers according to the values given \n \
               in the DropConfig.  \n \
               [example: DropConfig=0.0+0.5+0.3 \n \
                meaning: drop no nodes in layer 0 (input layer), half of the nodes in layer 1 and 30% of the nodes \n \
                in layer 2 \n \
                recommended: leave all the nodes turned on for the input layer (layer 0) \n \
                turn off half of the nodes in later layers for the initial training; leave all nodes \n \
                turned on (0.0) in later training stages] \n \
             . DropRepetitions \n \
               Each \"DropRepetitions\" cycle the configuration of which nodes are dropped is changed \n \
               [recommended : 1] \n \
             . Multithreading \n \
               turn on multithreading [recommended: True] \n \
               \n";
   Log () << txt << Endl;
}

} // namespace TMVA
