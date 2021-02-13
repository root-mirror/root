// @(#)root/tmva/tmva/cnn:$Id$Ndl
// Authors: Vladimir Ilievski, Lorenzo Moneta, Saurav Shekhar, Ravi Kiran
/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : MethodDL                                                              *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Deep Neural Network Method                                                *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Vladimir Ilievski  <ilievski.vladimir@live.com> - CERN, Switzerland       *
 *      Saurav Shekhar     <sauravshekhar01@gmail.com> - ETH Zurich, Switzerland  *
 *      Ravi Kiran S       <sravikiran0606@gmail.com> - CERN, Switzerland         *
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

#include "TFormula.h"
#include "TString.h"
#include "TMath.h"
#include "TObjString.h"

#include "TMVA/Tools.h"
#include "TMVA/Configurable.h"
#include "TMVA/IMethod.h"
#include "TMVA/ClassifierFactory.h"
#include "TMVA/MethodDL.h"
#include "TMVA/Types.h"
#include "TMVA/DNN/TensorDataLoader.h"
#include "TMVA/DNN/Functions.h"
#include "TMVA/DNN/DLMinimizers.h"
#include "TMVA/DNN/SGD.h"
#include "TMVA/DNN/Adam.h"
#include "TMVA/DNN/Adagrad.h"
#include "TMVA/DNN/RMSProp.h"
#include "TMVA/DNN/Adadelta.h"
#include "TMVA/Timer.h"

#ifdef R__HAS_TMVAGPU
#include "TMVA/DNN/Architectures/Cuda.h"
#ifdef R__HAS_CUDNN
#include "TMVA/DNN/Architectures/TCudnn.h"
#endif
#endif

#include <chrono>

REGISTER_METHOD(DL)
ClassImp(TMVA::MethodDL);

using namespace TMVA::DNN::CNN;
using namespace TMVA::DNN;

using TMVA::DNN::EActivationFunction;
using TMVA::DNN::ELossFunction;
using TMVA::DNN::EInitialization;
using TMVA::DNN::EOutputFunction;
using TMVA::DNN::EOptimizer;


namespace TMVA {


////////////////////////////////////////////////////////////////////////////////
TString fetchValueTmp(const std::map<TString, TString> &keyValueMap, TString key)
{
   key.ToUpper();
   std::map<TString, TString>::const_iterator it = keyValueMap.find(key);
   if (it == keyValueMap.end()) {
      return TString("");
   }
   return it->second;
}

////////////////////////////////////////////////////////////////////////////////
template <typename T>
T fetchValueTmp(const std::map<TString, TString> &keyValueMap, TString key, T defaultValue);

////////////////////////////////////////////////////////////////////////////////
template <>
int fetchValueTmp(const std::map<TString, TString> &keyValueMap, TString key, int defaultValue)
{
   TString value(fetchValueTmp(keyValueMap, key));
   if (value == "") {
      return defaultValue;
   }
   return value.Atoi();
}

////////////////////////////////////////////////////////////////////////////////
template <>
double fetchValueTmp(const std::map<TString, TString> &keyValueMap, TString key, double defaultValue)
{
   TString value(fetchValueTmp(keyValueMap, key));
   if (value == "") {
      return defaultValue;
   }
   return value.Atof();
}

////////////////////////////////////////////////////////////////////////////////
template <>
TString fetchValueTmp(const std::map<TString, TString> &keyValueMap, TString key, TString defaultValue)
{
   TString value(fetchValueTmp(keyValueMap, key));
   if (value == "") {
      return defaultValue;
   }
   return value;
}

////////////////////////////////////////////////////////////////////////////////
template <>
bool fetchValueTmp(const std::map<TString, TString> &keyValueMap, TString key, bool defaultValue)
{
   TString value(fetchValueTmp(keyValueMap, key));
   if (value == "") {
      return defaultValue;
   }

   value.ToUpper();
   if (value == "TRUE" || value == "T" || value == "1") {
      return true;
   }

   return false;
}

////////////////////////////////////////////////////////////////////////////////
template <>
std::vector<double> fetchValueTmp(const std::map<TString, TString> &keyValueMap, TString key,
                                  std::vector<double> defaultValue)
{
   TString parseString(fetchValueTmp(keyValueMap, key));
   if (parseString == "") {
      return defaultValue;
   }

   parseString.ToUpper();
   std::vector<double> values;

   const TString tokenDelim("+");
   TObjArray *tokenStrings = parseString.Tokenize(tokenDelim);
   TIter nextToken(tokenStrings);
   TObjString *tokenString = (TObjString *)nextToken();
   for (; tokenString != NULL; tokenString = (TObjString *)nextToken()) {
      std::stringstream sstr;
      double currentValue;
      sstr << tokenString->GetString().Data();
      sstr >> currentValue;
      values.push_back(currentValue);
   }
   return values;
}

////////////////////////////////////////////////////////////////////////////////
void MethodDL::DeclareOptions()
{
   // Set default values for all option strings

   DeclareOptionRef(fInputLayoutString = "0|0|0", "InputLayout", "The Layout of the input");

   DeclareOptionRef(fBatchLayoutString = "0|0|0", "BatchLayout", "The Layout of the batch");

   DeclareOptionRef(fLayoutString = "DENSE|(N+100)*2|SOFTSIGN,DENSE|0|LINEAR", "Layout", "Layout of the network.");

   DeclareOptionRef(fErrorStrategy = "CROSSENTROPY", "ErrorStrategy", "Loss function: Mean squared error (regression)"
                                                                      " or cross entropy (binary classification).");
   AddPreDefVal(TString("CROSSENTROPY"));
   AddPreDefVal(TString("SUMOFSQUARES"));
   AddPreDefVal(TString("MUTUALEXCLUSIVE"));

   DeclareOptionRef(fWeightInitializationString = "XAVIER", "WeightInitialization", "Weight initialization strategy");
   AddPreDefVal(TString("XAVIER"));
   AddPreDefVal(TString("XAVIERUNIFORM"));
   AddPreDefVal(TString("GAUSS"));
   AddPreDefVal(TString("UNIFORM"));
   AddPreDefVal(TString("IDENTITY"));
   AddPreDefVal(TString("ZERO"));

   DeclareOptionRef(fRandomSeed = 0, "RandomSeed", "Random seed used for weight initialization and batch shuffling");

   DeclareOptionRef(fNumValidationString = "20%", "ValidationSize", "Part of the training data to use for validation. "
                    "Specify as 0.2 or 20% to use a fifth of the data set as validation set. "
                    "Specify as 100 to use exactly 100 events. (Default: 20%)");

   DeclareOptionRef(fArchitectureString = "CPU", "Architecture", "Which architecture to perform the training on.");
   AddPreDefVal(TString("STANDARD"));   // deprecated and not supported anymore
   AddPreDefVal(TString("CPU"));
   AddPreDefVal(TString("GPU"));
   AddPreDefVal(TString("OPENCL"));    // not yet implemented
   AddPreDefVal(TString("CUDNN"));     // not needed (by default GPU is now CUDNN if available)

   // define training strategy separated by a separator "|"
   DeclareOptionRef(fTrainingStrategyString = "LearningRate=1e-3,"
                                              "Momentum=0.0,"
                                              "ConvergenceSteps=100,"
                                              "MaxEpochs=2000,"
                                              "Optimizer=ADAM,"
                                              "BatchSize=30,"
                                              "TestRepetitions=1,"
                                              "WeightDecay=0.0,"
                                              "Regularization=None,"
                                              "DropConfig=0.0",
                    "TrainingStrategy", "Defines the training strategies.");
}

////////////////////////////////////////////////////////////////////////////////
void MethodDL::ProcessOptions()
{

   if (IgnoreEventsWithNegWeightsInTraining()) {
      Log() << kINFO << "Will ignore negative events in training!" << Endl;
   }

   if (fArchitectureString == "STANDARD") {
      Log() << kWARNING << "The STANDARD architecture is not supported anymore. "
                         "Please use Architecture=CPU or Architecture=CPU."
                         "See the TMVA Users' Guide for instructions if you "
                         "encounter problems."
            << Endl;
      Log() << kINFO << "We will use instead the CPU architecture" << Endl;
      fArchitectureString = "CPU";
   }
   if (fArchitectureString == "OPENCL") {
      Log() << kERROR << "The OPENCL architecture has not been implemented yet. "
                         "Please use Architecture=CPU or Architecture=CPU for the "
                         "time being. See the TMVA Users' Guide for instructions "
                         "if you encounter problems."
            << Endl;
      // use instead GPU
      Log() << kINFO << "We will try using the GPU-CUDA architecture if available" << Endl;
      fArchitectureString = "GPU";
   }

   // the architecture can now be set at runtime as an option


   if (fArchitectureString == "GPU" || fArchitectureString == "CUDNN") {
#ifdef R__HAS_TMVAGPU
      Log() << kINFO << "Will now use the GPU architecture !" << Endl;
#else  // case TMVA does not support GPU
      Log() << kERROR << "CUDA backend not enabled. Please make sure "
         "you have CUDA installed and it was successfully "
         "detected by CMAKE by using -Dtmva-gpu=On  "
            << Endl;
      fArchitectureString = "CPU";
      Log() << kINFO << "Will now use instead the CPU architecture !" << Endl;
#endif
   }

   if (fArchitectureString == "CPU") {
#ifdef R__HAS_TMVACPU  // TMVA has CPU BLAS and IMT support
      Log() << kINFO << "Will now use the CPU architecture with BLAS and IMT support !" << Endl;
#else  // TMVA has no CPU BLAS or IMT support
      Log() << kINFO << "Multi-core CPU backend not enabled. For better performances, make sure "
                          "you have a BLAS implementation and it was successfully "
                         "detected by CMake as well that the imt CMake flag is set."
            << Endl;
      Log() << kINFO << "Will use anyway the CPU architecture but with slower performance" << Endl;
#endif
   }

   // Input Layout
   ParseInputLayout();
   ParseBatchLayout();

   // Loss function and output.
   fOutputFunction = EOutputFunction::kSigmoid;
   if (fAnalysisType == Types::kClassification) {
      if (fErrorStrategy == "SUMOFSQUARES") {
         fLossFunction = ELossFunction::kMeanSquaredError;
      }
      if (fErrorStrategy == "CROSSENTROPY") {
         fLossFunction = ELossFunction::kCrossEntropy;
      }
      fOutputFunction = EOutputFunction::kSigmoid;
   } else if (fAnalysisType == Types::kRegression) {
      if (fErrorStrategy != "SUMOFSQUARES") {
         Log() << kWARNING << "For regression only SUMOFSQUARES is a valid "
               << " neural net error function. Setting error function to "
               << " SUMOFSQUARES now." << Endl;
      }

      fLossFunction = ELossFunction::kMeanSquaredError;
      fOutputFunction = EOutputFunction::kIdentity;
   } else if (fAnalysisType == Types::kMulticlass) {
      if (fErrorStrategy == "SUMOFSQUARES") {
         fLossFunction = ELossFunction::kMeanSquaredError;
      }
      if (fErrorStrategy == "CROSSENTROPY") {
         fLossFunction = ELossFunction::kCrossEntropy;
      }
      if (fErrorStrategy == "MUTUALEXCLUSIVE") {
         fLossFunction = ELossFunction::kSoftmaxCrossEntropy;
      }
      fOutputFunction = EOutputFunction::kSoftmax;
   }

   // Initialization
   // the biases will be always initialized to zero
   if (fWeightInitializationString == "XAVIER") {
      fWeightInitialization = DNN::EInitialization::kGlorotNormal;
   } else if (fWeightInitializationString == "XAVIERUNIFORM") {
      fWeightInitialization = DNN::EInitialization::kGlorotUniform;
   } else if (fWeightInitializationString == "GAUSS") {
      fWeightInitialization = DNN::EInitialization::kGauss;
   } else if (fWeightInitializationString == "UNIFORM") {
      fWeightInitialization = DNN::EInitialization::kUniform;
   } else if (fWeightInitializationString == "ZERO") {
      fWeightInitialization = DNN::EInitialization::kZero;
   } else if (fWeightInitializationString == "IDENTITY") {
      fWeightInitialization = DNN::EInitialization::kIdentity;
   } else {
      fWeightInitialization = DNN::EInitialization::kGlorotUniform;
   }

   // Training settings.

   KeyValueVector_t strategyKeyValues = ParseKeyValueString(fTrainingStrategyString, TString("|"), TString(","));
   for (auto &block : strategyKeyValues) {
      TTrainingSettings settings;

      settings.convergenceSteps = fetchValueTmp(block, "ConvergenceSteps", 100);
      settings.batchSize = fetchValueTmp(block, "BatchSize", 30);
      settings.maxEpochs = fetchValueTmp(block, "MaxEpochs", 2000);
      settings.testInterval = fetchValueTmp(block, "TestRepetitions", 7);
      settings.weightDecay = fetchValueTmp(block, "WeightDecay", 0.0);
      settings.learningRate = fetchValueTmp(block, "LearningRate", 1e-5);
      settings.momentum = fetchValueTmp(block, "Momentum", 0.3);
      settings.dropoutProbabilities = fetchValueTmp(block, "DropConfig", std::vector<Double_t>());

      TString regularization = fetchValueTmp(block, "Regularization", TString("NONE"));
      if (regularization == "L1") {
         settings.regularization = DNN::ERegularization::kL1;
      } else if (regularization == "L2") {
         settings.regularization = DNN::ERegularization::kL2;
      } else {
         settings.regularization = DNN::ERegularization::kNone;
      }

      TString optimizer = fetchValueTmp(block, "Optimizer", TString("ADAM"));
      settings.optimizerName = optimizer;
      if (optimizer == "SGD") {
         settings.optimizer = DNN::EOptimizer::kSGD;
      } else if (optimizer == "ADAM") {
         settings.optimizer = DNN::EOptimizer::kAdam;
      } else if (optimizer == "ADAGRAD") {
         settings.optimizer = DNN::EOptimizer::kAdagrad;
      } else if (optimizer == "RMSPROP") {
         settings.optimizer = DNN::EOptimizer::kRMSProp;
      } else if (optimizer == "ADADELTA") {
         settings.optimizer = DNN::EOptimizer::kAdadelta;
      } else {
         // Make Adam as default choice if the input string is
         // incorrect.
         settings.optimizer = DNN::EOptimizer::kAdam;
         settings.optimizerName = "ADAM";
      }

      fTrainingSettings.push_back(settings);
   }

   // this set fInputShape[0] = batchSize
   this->SetBatchSize(fTrainingSettings.front().batchSize);

   // case inputlayout and batch layout was not given. Use default then
   // (1, batchsize, nvariables)
   // fInputShape[0] -> BatchSize
   // fInputShape[1] -> InputDepth
   // fInputShape[2] -> InputHeight
   // fInputShape[3] -> InputWidth
   if (fInputShape[3] == 0 && fInputShape[2] == 0 && fInputShape[1] == 0) {
      fInputShape[1] = 1;
      fInputShape[2] = 1;
      fInputShape[3] = GetNVariables();
   }
   // case when batch layout is not provided (all zero)
   // batch layout can be determined by the input layout + batch size
   // case DNN : { 1, B, W }
   // case CNN :  { B, C, H*W}
   // case RNN :  { B, T, H*W }

   if (fBatchWidth == 0 && fBatchHeight == 0 && fBatchDepth == 0) {
      // case first layer is DENSE
      if (fInputShape[2] == 1 && fInputShape[1] == 1) {
         // case of (1, batchsize, input features)
         fBatchDepth  = 1;
         fBatchHeight = fTrainingSettings.front().batchSize;
         fBatchWidth  = fInputShape[3];
      }
      else { // more general cases (e.g. for CNN)
         // case CONV or RNN
         fBatchDepth  = fTrainingSettings.front().batchSize;
         fBatchHeight = fInputShape[1];
         fBatchWidth  = fInputShape[3]*fInputShape[2];
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// default initializations
void MethodDL::Init()
{
   // Nothing to do here
}

////////////////////////////////////////////////////////////////////////////////
/// Parse the input layout
void MethodDL::ParseInputLayout()
{
   // Define the delimiter
   const TString delim("|");

   // Get the input layout string
   TString inputLayoutString = this->GetInputLayoutString();

   // Split the input layout string
   TObjArray *inputDimStrings = inputLayoutString.Tokenize(delim);
   TIter nextInputDim(inputDimStrings);
   TObjString *inputDimString = (TObjString *)nextInputDim();

   // Go through every token and save its absolute value in the shape array
   // The first token is the batch size for easy compatibility with cudnn
   int subDim = 1;
   std::vector<size_t> inputShape;
   inputShape.reserve(inputLayoutString.Length()/2 + 2);
   inputShape.push_back(0);    // Will be set later by Trainingsettings, use 0 value now
   for (; inputDimString != nullptr; inputDimString = (TObjString *)nextInputDim()) {
      // size_t is unsigned
      subDim = (size_t) abs(inputDimString->GetString().Atoi());
      // Size among unused dimensions should be set to 1 for cudnn
      //if (subDim == 0) subDim = 1;
      inputShape.push_back(subDim);
   }
   // it is expected that empty Shape has at least 4 dimensions. We pad the missing one's with 1
   // for example in case of dense layer input layouts
   // when we will support 3D convolutions we would need to add extra 1's
   if (inputShape.size() == 2) {
      // case of dense layer where only width is specified
      inputShape.insert(inputShape.begin() + 1, {1,1});
   }
   else if (inputShape.size() == 3) {
      //e.g. case of RNN T,W -> T,1,W
      inputShape.insert(inputShape.begin() + 2, 1);
   }

   this->SetInputShape(inputShape);
}

////////////////////////////////////////////////////////////////////////////////
/// Parse the input layout
void MethodDL::ParseBatchLayout()
{
   // Define the delimiter
   const TString delim("|");

   // Get the input layout string
   TString batchLayoutString = this->GetBatchLayoutString();

   size_t batchDepth = 0;
   size_t batchHeight = 0;
   size_t batchWidth = 0;

   // Split the input layout string
   TObjArray *batchDimStrings = batchLayoutString.Tokenize(delim);
   TIter nextBatchDim(batchDimStrings);
   TObjString *batchDimString = (TObjString *)nextBatchDim();
   int idxToken = 0;

   for (; batchDimString != nullptr; batchDimString = (TObjString *)nextBatchDim()) {
      switch (idxToken) {
      case 0: // input depth
      {
         TString strDepth(batchDimString->GetString());
         batchDepth = (size_t)strDepth.Atoi();
      } break;
      case 1: // input height
      {
         TString strHeight(batchDimString->GetString());
         batchHeight = (size_t)strHeight.Atoi();
      } break;
      case 2: // input width
      {
         TString strWidth(batchDimString->GetString());
         batchWidth = (size_t)strWidth.Atoi();
      } break;
      }
      ++idxToken;
   }

   this->SetBatchDepth(batchDepth);
   this->SetBatchHeight(batchHeight);
   this->SetBatchWidth(batchWidth);
}

////////////////////////////////////////////////////////////////////////////////
/// Create a deep net based on the layout string
template <typename Architecture_t, typename Layer_t>
void MethodDL::CreateDeepNet(DNN::TDeepNet<Architecture_t, Layer_t> &deepNet,
                             std::vector<DNN::TDeepNet<Architecture_t, Layer_t>> &nets)
{
   // Layer specification, layer details
   const TString layerDelimiter(",");
   const TString subDelimiter("|");

   TString layoutString = this->GetLayoutString();

   //std::cout << "Create Deepnet - layout string " << layoutString << "\t layers : " << deepNet.GetLayers().size() << std::endl;

   // Split layers
   TObjArray *layerStrings = layoutString.Tokenize(layerDelimiter);
   TIter nextLayer(layerStrings);
   TObjString *layerString = (TObjString *)nextLayer();


   for (; layerString != nullptr; layerString = (TObjString *)nextLayer()) {

      // Split layer details
      TObjArray *subStrings = layerString->GetString().Tokenize(subDelimiter);
      TIter nextToken(subStrings);
      TObjString *token = (TObjString *)nextToken();

      // Determine the type of the layer
      TString strLayerType = token->GetString();


      if (strLayerType == "DENSE") {
         ParseDenseLayer(deepNet, nets, layerString->GetString(), subDelimiter);
      } else if (strLayerType == "CONV") {
         ParseConvLayer(deepNet, nets, layerString->GetString(), subDelimiter);
      } else if (strLayerType == "MAXPOOL") {
         ParseMaxPoolLayer(deepNet, nets, layerString->GetString(), subDelimiter);
      } else if (strLayerType == "RESHAPE") {
         ParseReshapeLayer(deepNet, nets, layerString->GetString(), subDelimiter);
      } else if (strLayerType == "BNORM") {
         ParseBatchNormLayer(deepNet, nets, layerString->GetString(), subDelimiter);
      } else if (strLayerType == "RNN") {
         ParseRecurrentLayer(kLayerRNN, deepNet, nets, layerString->GetString(), subDelimiter);
      } else if (strLayerType == "LSTM") {
         ParseRecurrentLayer(kLayerLSTM, deepNet, nets, layerString->GetString(), subDelimiter);
      } else if (strLayerType == "GRU") {
         ParseRecurrentLayer(kLayerGRU, deepNet, nets, layerString->GetString(), subDelimiter);
      } else {
         // no type of layer specified - assume is dense layer as in old DNN interface
         ParseDenseLayer(deepNet, nets, layerString->GetString(), subDelimiter);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Pases the layer string and creates the appropriate dense layer
template <typename Architecture_t, typename Layer_t>
void MethodDL::ParseDenseLayer(DNN::TDeepNet<Architecture_t, Layer_t> &deepNet,
                               std::vector<DNN::TDeepNet<Architecture_t, Layer_t>> & /*nets*/, TString layerString,
                               TString delim)
{
   int width = 0;
   EActivationFunction activationFunction = EActivationFunction::kTanh;

   // this return number of input variables for the method
   // it can be used to deduce width of dense layer if specified as N+10
   // where N is the number of input variables
   const size_t inputSize = GetNvar();

   // Split layer details
   TObjArray *subStrings = layerString.Tokenize(delim);
   TIter nextToken(subStrings);
   TObjString *token = (TObjString *)nextToken();
   int idxToken = 0;

   // loop on the tokens
   // order of sepcifying width and activation function is not relevant
   // both  100|TANH and TANH|100 are valid cases
   for (; token != nullptr; token = (TObjString *)nextToken()) {
      idxToken++;
      // try a match with the activation function
      TString strActFnc(token->GetString());
      // if first token defines the layer type- skip it
      if (strActFnc =="DENSE") continue;

      if (strActFnc == "RELU") {
         activationFunction = DNN::EActivationFunction::kRelu;
      } else if (strActFnc == "TANH") {
         activationFunction = DNN::EActivationFunction::kTanh;
      } else if (strActFnc == "FTANH") {
         activationFunction = DNN::EActivationFunction::kFastTanh;
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
      } else if (width == 0) {
         // no match found try to parse as text showing the width
         // support for input a formula where the variable 'x' is 'N' in the string
         // use TFormula for the evaluation
         TString  strNumNodes = strActFnc;
         // number of nodes
         TString strN("x");
         strNumNodes.ReplaceAll("N", strN);
         strNumNodes.ReplaceAll("n", strN);
         TFormula fml("tmp", strNumNodes);
         width = fml.Eval(inputSize);
      }
   }
   // avoid zero width. assume is last layer and give width = output width
   // Determine the number of outputs
   size_t outputSize = 1;
   if (fAnalysisType == Types::kRegression && GetNTargets() != 0) {
      outputSize = GetNTargets();
   } else if (fAnalysisType == Types::kMulticlass && DataInfo().GetNClasses() >= 2) {
      outputSize = DataInfo().GetNClasses();
   }
   if (width == 0) width = outputSize;

   // Add the dense layer, initialize the weights and biases and copy
   TDenseLayer<Architecture_t> *denseLayer = deepNet.AddDenseLayer(width, activationFunction);
   denseLayer->Initialize();

   // add same layer to fNet
   if (fBuildNet) fNet->AddDenseLayer(width, activationFunction);

   //TDenseLayer<Architecture_t> *copyDenseLayer = new TDenseLayer<Architecture_t>(*denseLayer);

   // add the copy to all slave nets
   //for (size_t i = 0; i < nets.size(); i++) {
   //   nets[i].AddDenseLayer(copyDenseLayer);
   //}

   // check compatibility of added layer
   // for a dense layer input should be 1 x 1 x DxHxW
}

////////////////////////////////////////////////////////////////////////////////
/// Pases the layer string and creates the appropriate convolutional layer
template <typename Architecture_t, typename Layer_t>
void MethodDL::ParseConvLayer(DNN::TDeepNet<Architecture_t, Layer_t> &deepNet,
                              std::vector<DNN::TDeepNet<Architecture_t, Layer_t>> & /*nets*/, TString layerString,
                              TString delim)
{
   int depth = 0;
   int fltHeight = 0;
   int fltWidth = 0;
   int strideRows = 0;
   int strideCols = 0;
   int zeroPadHeight = 0;
   int zeroPadWidth = 0;
   EActivationFunction activationFunction = EActivationFunction::kTanh;

   // Split layer details
   TObjArray *subStrings = layerString.Tokenize(delim);
   TIter nextToken(subStrings);
   TObjString *token = (TObjString *)nextToken();
   int idxToken = 0;

   for (; token != nullptr; token = (TObjString *)nextToken()) {
      switch (idxToken) {
      case 1: // depth
      {
         TString strDepth(token->GetString());
         depth = strDepth.Atoi();
      } break;
      case 2: // filter height
      {
         TString strFltHeight(token->GetString());
         fltHeight = strFltHeight.Atoi();
      } break;
      case 3: // filter width
      {
         TString strFltWidth(token->GetString());
         fltWidth = strFltWidth.Atoi();
      } break;
      case 4: // stride in rows
      {
         TString strStrideRows(token->GetString());
         strideRows = strStrideRows.Atoi();
      } break;
      case 5: // stride in cols
      {
         TString strStrideCols(token->GetString());
         strideCols = strStrideCols.Atoi();
      } break;
      case 6: // zero padding height
      {
         TString strZeroPadHeight(token->GetString());
         zeroPadHeight = strZeroPadHeight.Atoi();
      } break;
      case 7: // zero padding width
      {
         TString strZeroPadWidth(token->GetString());
         zeroPadWidth = strZeroPadWidth.Atoi();
      } break;
      case 8: // activation function
      {
         TString strActFnc(token->GetString());
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
      } break;
      }
      ++idxToken;
   }

   // Add the convolutional layer, initialize the weights and biases and copy
   TConvLayer<Architecture_t> *convLayer = deepNet.AddConvLayer(depth, fltHeight, fltWidth, strideRows, strideCols,
                                                                zeroPadHeight, zeroPadWidth, activationFunction);
   convLayer->Initialize();

   // Add same layer to fNet
   if (fBuildNet) fNet->AddConvLayer(depth, fltHeight, fltWidth, strideRows, strideCols,
                      zeroPadHeight, zeroPadWidth, activationFunction);

   //TConvLayer<Architecture_t> *copyConvLayer = new TConvLayer<Architecture_t>(*convLayer);

   //// add the copy to all slave nets
   //for (size_t i = 0; i < nets.size(); i++) {
   //   nets[i].AddConvLayer(copyConvLayer);
   //}
}

////////////////////////////////////////////////////////////////////////////////
/// Pases the layer string and creates the appropriate max pool layer
template <typename Architecture_t, typename Layer_t>
void MethodDL::ParseMaxPoolLayer(DNN::TDeepNet<Architecture_t, Layer_t> &deepNet,
                                 std::vector<DNN::TDeepNet<Architecture_t, Layer_t>> & /*nets*/, TString layerString,
                                 TString delim)
{

   int filterHeight = 0;
   int filterWidth = 0;
   int strideRows = 0;
   int strideCols = 0;

   // Split layer details
   TObjArray *subStrings = layerString.Tokenize(delim);
   TIter nextToken(subStrings);
   TObjString *token = (TObjString *)nextToken();
   int idxToken = 0;

   for (; token != nullptr; token = (TObjString *)nextToken()) {
      switch (idxToken) {
      case 1: // filter height
      {
         TString strFrmHeight(token->GetString());
         filterHeight = strFrmHeight.Atoi();
      } break;
      case 2: // filter width
      {
         TString strFrmWidth(token->GetString());
         filterWidth = strFrmWidth.Atoi();
      } break;
      case 3: // stride in rows
      {
         TString strStrideRows(token->GetString());
         strideRows = strStrideRows.Atoi();
      } break;
      case 4: // stride in cols
      {
         TString strStrideCols(token->GetString());
         strideCols = strStrideCols.Atoi();
      } break;
      }
      ++idxToken;
   }

   // Add the Max pooling layer
   // TMaxPoolLayer<Architecture_t> *maxPoolLayer =
   deepNet.AddMaxPoolLayer(filterHeight, filterWidth, strideRows, strideCols);

   // Add the same layer to fNet
   if (fBuildNet) fNet->AddMaxPoolLayer(filterHeight, filterWidth, strideRows, strideCols);


   //TMaxPoolLayer<Architecture_t> *copyMaxPoolLayer = new TMaxPoolLayer<Architecture_t>(*maxPoolLayer);

   //// add the copy to all slave nets
   //for (size_t i = 0; i < nets.size(); i++) {
   //   nets[i].AddMaxPoolLayer(copyMaxPoolLayer);
   //}
}

////////////////////////////////////////////////////////////////////////////////
/// Pases the layer string and creates the appropriate reshape layer
template <typename Architecture_t, typename Layer_t>
void MethodDL::ParseReshapeLayer(DNN::TDeepNet<Architecture_t, Layer_t> &deepNet,
                                 std::vector<DNN::TDeepNet<Architecture_t, Layer_t>> & /*nets*/, TString layerString,
                                 TString delim)
{
   int depth = 0;
   int height = 0;
   int width = 0;
   bool flattening = false;

   // Split layer details
   TObjArray *subStrings = layerString.Tokenize(delim);
   TIter nextToken(subStrings);
   TObjString *token = (TObjString *)nextToken();
   int idxToken = 0;

   for (; token != nullptr; token = (TObjString *)nextToken()) {
      if (token->GetString() == "FLAT") idxToken=4;
      switch (idxToken) {
      case 1: {
         TString strDepth(token->GetString());
         depth = strDepth.Atoi();
      } break;
      case 2: // height
      {
         TString strHeight(token->GetString());
         height = strHeight.Atoi();
      } break;
      case 3: // width
      {
         TString strWidth(token->GetString());
         width = strWidth.Atoi();
      } break;
      case 4: // flattening
      {
         TString flat(token->GetString());
         if (flat == "FLAT") {
            flattening = true;
         }
      } break;
      }
      ++idxToken;
   }

   // Add the reshape layer
   // TReshapeLayer<Architecture_t> *reshapeLayer =
   deepNet.AddReshapeLayer(depth, height, width, flattening);

   // Add the same layer to fNet
   if (fBuildNet) fNet->AddReshapeLayer(depth, height, width, flattening);

   //TReshapeLayer<Architecture_t> *copyReshapeLayer = new TReshapeLayer<Architecture_t>(*reshapeLayer);

   //// add the copy to all slave nets
   //for (size_t i = 0; i < nets.size(); i++) {
   //   nets[i].AddReshapeLayer(copyReshapeLayer);
   //}
}

////////////////////////////////////////////////////////////////////////////////
/// Pases the layer string and creates the appropriate reshape layer
template <typename Architecture_t, typename Layer_t>
void MethodDL::ParseBatchNormLayer(DNN::TDeepNet<Architecture_t, Layer_t> &deepNet,
                                 std::vector<DNN::TDeepNet<Architecture_t, Layer_t>> & /*nets*/, TString layerString,
                                 TString delim)
{

   // default values
   double momentum = -1; //0.99;
   double epsilon = 0.0001;

   // Split layer details
   TObjArray *subStrings = layerString.Tokenize(delim);
   TIter nextToken(subStrings);
   TObjString *token = (TObjString *)nextToken();
   int idxToken = 0;

   for (; token != nullptr; token = (TObjString *)nextToken()) {
      switch (idxToken) {
      case 1: {
         momentum = std::atof(token->GetString().Data());
      } break;
      case 2: // height
      {
         epsilon = std::atof(token->GetString().Data());
      } break;
      }
      ++idxToken;
   }

   // Add the batch norm  layer
   //
   auto layer = deepNet.AddBatchNormLayer(momentum, epsilon);
   layer->Initialize();

   // Add the same layer to fNet
   if (fBuildNet) fNet->AddBatchNormLayer(momentum, epsilon);

}

////////////////////////////////////////////////////////////////////////////////
/// Pases the layer string and creates the appropriate rnn layer
template <typename Architecture_t, typename Layer_t>
void MethodDL::ParseRecurrentLayer(ERecurrentLayerType rnnType, DNN::TDeepNet<Architecture_t, Layer_t> & deepNet,
                             std::vector<DNN::TDeepNet<Architecture_t, Layer_t>> & /*nets */, TString layerString,
                             TString delim)
{
   //    int depth = 0;
   int stateSize = 0;
   int inputSize = 0;
   int timeSteps = 0;
   bool rememberState = false;
   bool returnSequence = false;
   bool resetGateAfter = false;

   // Split layer details
   TObjArray *subStrings = layerString.Tokenize(delim);
   TIter nextToken(subStrings);
   TObjString *token = (TObjString *)nextToken();
   int idxToken = 0;

   for (; token != nullptr; token = (TObjString *)nextToken()) {
      switch (idxToken) {
         case 1:  // state size
         {
            TString strstateSize(token->GetString());
            stateSize = strstateSize.Atoi();
            break;
         }
         case 2:  // input size
         {
            TString strinputSize(token->GetString());
            inputSize = strinputSize.Atoi();
            break;
         }
         case 3:  // time steps
         {
            TString strtimeSteps(token->GetString());
            timeSteps = strtimeSteps.Atoi();
            break;
         }
         case 4: // returnSequence (option stateful in Keras)
         {
            TString strrememberState(token->GetString());
            rememberState = (bool) strrememberState.Atoi();
            break;
         }
         case 5: // return full output sequence (1 or 0)
         {
            TString str(token->GetString());
            returnSequence = (bool)str.Atoi();
            break;
         }
         case 6: // resetGate after option (only for GRU)
         {
            TString str(token->GetString());
            resetGateAfter = (bool)str.Atoi();
         }
      }
      ++idxToken;
   }

   // Add the recurrent layer, initialize the weights and biases and copy
   if (rnnType == kLayerRNN) {
      auto  * recurrentLayer = deepNet.AddBasicRNNLayer(stateSize, inputSize, timeSteps, rememberState, returnSequence);
      recurrentLayer->Initialize();
      // Add same layer to fNet
      if (fBuildNet) fNet->AddBasicRNNLayer(stateSize, inputSize, timeSteps, rememberState, returnSequence);
   }
   else if (rnnType == kLayerLSTM ) {
      auto *recurrentLayer = deepNet.AddBasicLSTMLayer(stateSize, inputSize, timeSteps, rememberState, returnSequence);
      recurrentLayer->Initialize();
      // Add same layer to fNet
      if (fBuildNet)
         fNet->AddBasicLSTMLayer(stateSize, inputSize, timeSteps, rememberState, returnSequence);
   }
   else if (rnnType == kLayerGRU) {
      if (Architecture_t::IsCudnn()) resetGateAfter = true; // needed for Cudnn
      auto *recurrentLayer = deepNet.AddBasicGRULayer(stateSize, inputSize, timeSteps, rememberState, returnSequence, resetGateAfter);
      recurrentLayer->Initialize();
      // Add same layer to fNet
      if (fBuildNet)
         fNet->AddBasicGRULayer(stateSize, inputSize, timeSteps, rememberState, returnSequence, resetGateAfter);
   }
   else {
      Log() << kFATAL << "Invalid Recurrent layer type " << Endl;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Standard constructor.
MethodDL::MethodDL(const TString &jobName, const TString &methodTitle, DataSetInfo &theData, const TString &theOption)
   : MethodBase(jobName, Types::kDL, methodTitle, theData, theOption), fInputShape(4,0),
     fBatchHeight(), fBatchWidth(), fRandomSeed(0), fWeightInitialization(),
     fOutputFunction(), fLossFunction(), fInputLayoutString(), fBatchLayoutString(),
     fLayoutString(), fErrorStrategy(), fTrainingStrategyString(), fWeightInitializationString(),
     fArchitectureString(), fResume(false), fBuildNet(true), fTrainingSettings(),
     fXInput()
{
   // Nothing to do here
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor from a weight file.
MethodDL::MethodDL(DataSetInfo &theData, const TString &theWeightFile)
   : MethodBase(Types::kDL, theData, theWeightFile), fInputShape(4,0), fBatchHeight(),
     fBatchWidth(), fRandomSeed(0), fWeightInitialization(), fOutputFunction(),
     fLossFunction(), fInputLayoutString(), fBatchLayoutString(), fLayoutString(),
     fErrorStrategy(), fTrainingStrategyString(), fWeightInitializationString(),
     fArchitectureString(), fResume(false), fBuildNet(true), fTrainingSettings(),
     fXInput()
{
   // Nothing to do here
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.
MethodDL::~MethodDL()
{
   // Nothing to do here
}

////////////////////////////////////////////////////////////////////////////////
/// Parse key value pairs in blocks -> return vector of blocks with map of key value pairs.
auto MethodDL::ParseKeyValueString(TString parseString, TString blockDelim, TString tokenDelim) -> KeyValueVector_t
{
   // remove empty spaces
   parseString.ReplaceAll(" ","");
   KeyValueVector_t blockKeyValues;
   const TString keyValueDelim("=");

   TObjArray *blockStrings = parseString.Tokenize(blockDelim);
   TIter nextBlock(blockStrings);
   TObjString *blockString = (TObjString *)nextBlock();

   for (; blockString != nullptr; blockString = (TObjString *)nextBlock()) {
      blockKeyValues.push_back(std::map<TString, TString>());
      std::map<TString, TString> &currentBlock = blockKeyValues.back();

      TObjArray *subStrings = blockString->GetString().Tokenize(tokenDelim);
      TIter nextToken(subStrings);
      TObjString *token = (TObjString *)nextToken();

      for (; token != nullptr; token = (TObjString *)nextToken()) {
         TString strKeyValue(token->GetString());
         int delimPos = strKeyValue.First(keyValueDelim.Data());
         if (delimPos <= 0) continue;

         TString strKey = TString(strKeyValue(0, delimPos));
         strKey.ToUpper();
         TString strValue = TString(strKeyValue(delimPos + 1, strKeyValue.Length()));

         strKey.Strip(TString::kBoth, ' ');
         strValue.Strip(TString::kBoth, ' ');

         currentBlock.insert(std::make_pair(strKey, strValue));
      }
   }
   return blockKeyValues;
}

////////////////////////////////////////////////////////////////////////////////
/// What kind of analysis type can handle the CNN
Bool_t MethodDL::HasAnalysisType(Types::EAnalysisType type, UInt_t numberClasses, UInt_t /*numberTargets*/)
{
   if (type == Types::kClassification && numberClasses == 2) return kTRUE;
   if (type == Types::kMulticlass) return kTRUE;
   if (type == Types::kRegression) return kTRUE;

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Validation of the ValidationSize option. Allowed formats are 20%, 0.2 and
/// 100 etc.
///    - 20% and 0.2 selects 20% of the training set as validation data.
///    - 100 selects 100 events as the validation data.
///
/// @return number of samples in validation set
///
UInt_t TMVA::MethodDL::GetNumValidationSamples()
{
   Int_t nValidationSamples = 0;
   UInt_t trainingSetSize = GetEventCollection(Types::kTraining).size();

   // Parsing + Validation
   // --------------------
   if (fNumValidationString.EndsWith("%")) {
      // Relative spec. format 20%
      TString intValStr = TString(fNumValidationString.Strip(TString::kTrailing, '%'));

      if (intValStr.IsFloat()) {
         Double_t valSizeAsDouble = fNumValidationString.Atof() / 100.0;
         nValidationSamples = GetEventCollection(Types::kTraining).size() * valSizeAsDouble;
      } else {
         Log() << kFATAL << "Cannot parse number \"" << fNumValidationString
               << "\". Expected string like \"20%\" or \"20.0%\"." << Endl;
      }
   } else if (fNumValidationString.IsFloat()) {
      Double_t valSizeAsDouble = fNumValidationString.Atof();

      if (valSizeAsDouble < 1.0) {
         // Relative spec. format 0.2
         nValidationSamples = GetEventCollection(Types::kTraining).size() * valSizeAsDouble;
      } else {
         // Absolute spec format 100 or 100.0
         nValidationSamples = valSizeAsDouble;
      }
   } else {
      Log() << kFATAL << "Cannot parse number \"" << fNumValidationString << "\". Expected string like \"0.2\" or \"100\"."
            << Endl;
   }

   // Value validation
   // ----------------
   if (nValidationSamples < 0) {
      Log() << kFATAL << "Validation size \"" << fNumValidationString << "\" is negative." << Endl;
   }

   if (nValidationSamples == 0) {
      Log() << kFATAL << "Validation size \"" << fNumValidationString << "\" is zero." << Endl;
   }

   if (nValidationSamples >= (Int_t)trainingSetSize) {
      Log() << kFATAL << "Validation size \"" << fNumValidationString
            << "\" is larger than or equal in size to training set (size=\"" << trainingSetSize << "\")." << Endl;
   }

   return nValidationSamples;
}


////////////////////////////////////////////////////////////////////////////////
///  Implementation of architecture specific train method
///
template <typename Architecture_t>
void MethodDL::TrainDeepNet()
{

   using Scalar_t = typename Architecture_t::Scalar_t;
   using Layer_t = TMVA::DNN::VGeneralLayer<Architecture_t>;
   using DeepNet_t = TMVA::DNN::TDeepNet<Architecture_t, Layer_t>;
   using TensorDataLoader_t = TTensorDataLoader<TMVAInput_t, Architecture_t>;

   bool debug = Log().GetMinType() == kDEBUG;


   // set the random seed for weight initialization
   Architecture_t::SetRandomSeed(fRandomSeed);

   ///split training data in training and validation data
   // and determine the number of training and testing examples

   size_t nValidationSamples = GetNumValidationSamples();
   size_t nTrainingSamples = GetEventCollection(Types::kTraining).size() - nValidationSamples;

   const std::vector<TMVA::Event *> &allData = GetEventCollection(Types::kTraining);
   const std::vector<TMVA::Event *> eventCollectionTraining{allData.begin(), allData.begin() + nTrainingSamples};
   const std::vector<TMVA::Event *> eventCollectionValidation{allData.begin() + nTrainingSamples, allData.end()};

   size_t trainingPhase = 1;

   for (TTrainingSettings &settings : this->GetTrainingSettings()) {

      size_t nThreads = 1;       // FIXME threads are hard coded to 1, no use of slave threads or multi-threading


      // After the processing of the options, initialize the master deep net
      size_t batchSize = settings.batchSize;
      this->SetBatchSize(batchSize);
      // Should be replaced by actual implementation. No support for this now.
      size_t inputDepth  = this->GetInputDepth();
      size_t inputHeight = this->GetInputHeight();
      size_t inputWidth  = this->GetInputWidth();
      size_t batchDepth  = this->GetBatchDepth();
      size_t batchHeight = this->GetBatchHeight();
      size_t batchWidth  = this->GetBatchWidth();
      ELossFunction J    = this->GetLossFunction();
      EInitialization I  = this->GetWeightInitialization();
      ERegularization R    = settings.regularization;
      EOptimizer O         = settings.optimizer;
      Scalar_t weightDecay = settings.weightDecay;

      //Batch size should be included in batch layout as well. There are two possibilities:
      //  1.  Batch depth = batch size   one will input tensorsa as (batch_size x d1 x d2)
      //       This is case for example if first layer is a conv layer and d1 = image depth, d2 = image width x image height
      //  2.  Batch depth = 1, batch height = batch size  batxch width = dim of input features
      //        This should be case if first layer is a Dense 1 and input tensor must be ( 1 x batch_size x input_features )

      if (batchDepth != batchSize && batchDepth > 1) {
         Error("Train","Given batch depth of %zu (specified in BatchLayout)  should be equal to given batch size %zu",batchDepth,batchSize);
         return;
      }
      if (batchDepth == 1 && batchSize > 1 && batchSize != batchHeight ) {
         Error("Train","Given batch height of %zu (specified in BatchLayout)  should be equal to given batch size %zu",batchHeight,batchSize);
         return;
      }


      //check also that input layout compatible with batch layout
      bool badLayout = false;
      // case batch depth == batch size
      if (batchDepth == batchSize)
         badLayout = ( inputDepth * inputHeight * inputWidth != batchHeight * batchWidth ) ;
      // case batch Height is batch size
      if (batchHeight == batchSize && batchDepth == 1)
         badLayout |=  ( inputDepth * inputHeight * inputWidth !=  batchWidth);
      if (badLayout) {
         Error("Train","Given input layout %zu x %zu x %zu is not compatible with  batch layout %zu x %zu x  %zu ",
               inputDepth,inputHeight,inputWidth,batchDepth,batchHeight,batchWidth);
         return;
      }

      // check batch size is compatible with number of events
      if (nTrainingSamples < settings.batchSize || nValidationSamples < settings.batchSize) {
         Log() << kFATAL << "Number of samples in the datasets are train: ("
               << nTrainingSamples << ") test: (" << nValidationSamples
               << "). One of these is smaller than the batch size of "
               << settings.batchSize << ". Please increase the batch"
               << " size to be at least the same size as the smallest"
               << " of them." << Endl;
      }

      DeepNet_t deepNet(batchSize, inputDepth, inputHeight, inputWidth, batchDepth, batchHeight, batchWidth, J, I, R, weightDecay);

      // create a copy of DeepNet for evaluating but with batch size = 1
      // fNet is the saved network and will be with CPU or Referrence architecture
      if (trainingPhase == 1) {
         fNet = std::unique_ptr<DeepNetImpl_t>(new DeepNetImpl_t(1, inputDepth, inputHeight, inputWidth, batchDepth,
                                                                 batchHeight, batchWidth, J, I, R, weightDecay));
         fBuildNet = true;
      }
      else
         fBuildNet = false;

      // Initialize the vector of slave nets
      std::vector<DeepNet_t> nets{};
      nets.reserve(nThreads);
      for (size_t i = 0; i < nThreads; i++) {
         // create a copies of the master deep net
         nets.push_back(deepNet);
      }


      // Add all appropriate layers to deepNet and (if fBuildNet is true) also to fNet
      CreateDeepNet(deepNet, nets);


      // set droput probabilities
      // use convention to store in the layer 1.- dropout probabilities
      std::vector<Double_t> dropoutVector(settings.dropoutProbabilities);
      for (auto & p : dropoutVector) {
         p = 1.0 - p;
      }
      deepNet.SetDropoutProbabilities(dropoutVector);

      if (trainingPhase > 1) {
         // copy initial weights from fNet to deepnet
         for (size_t i = 0; i < deepNet.GetDepth(); ++i) {
            deepNet.GetLayerAt(i)->CopyParameters(*fNet->GetLayerAt(i));
         }
      }

      // when fNet is built create also input matrix that will be used to evaluate it
      if (fBuildNet) {
         //int n1 = batchHeight;
         //int n2 = batchWidth;
         // treat case where batchHeight is the batchSize in case of first Dense layers (then we need to set to fNet batch size)
         //if (batchDepth == 1 && GetInputHeight() == 1 && GetInputDepth() == 1) n1 = fNet->GetBatchSize();
         //fXInput = TensorImpl_t(1,n1,n2);
         fXInput = ArchitectureImpl_t::CreateTensor(fNet->GetBatchSize(), GetInputDepth(), GetInputHeight(), GetInputWidth() );
         if (batchDepth == 1 && GetInputHeight() == 1 && GetInputDepth() == 1)
            fXInput = TensorImpl_t( fNet->GetBatchSize(), GetInputWidth() );
         fXInputBuffer = HostBufferImpl_t( fXInput.GetSize() );


         // create pointer to output matrix used for the predictions
         fYHat = std::unique_ptr<MatrixImpl_t>(new MatrixImpl_t(fNet->GetBatchSize(),  fNet->GetOutputWidth() ) );

         // print the created network
         Log()  << "*****   Deep Learning Network *****" << Endl;
         if (Log().GetMinType() <= kINFO)
            deepNet.Print();
      }
      Log() << "Using " << nTrainingSamples << " events for training and " <<  nValidationSamples << " for testing" << Endl;

      // Loading the training and validation datasets
      TMVAInput_t trainingTuple = std::tie(eventCollectionTraining, DataInfo());
      TensorDataLoader_t trainingData(trainingTuple, nTrainingSamples, batchSize,
                                      {inputDepth, inputHeight, inputWidth},
                                     {deepNet.GetBatchDepth(), deepNet.GetBatchHeight(), deepNet.GetBatchWidth()} ,
                                      deepNet.GetOutputWidth(), nThreads);

      TMVAInput_t validationTuple = std::tie(eventCollectionValidation, DataInfo());
      TensorDataLoader_t validationData(validationTuple, nValidationSamples, batchSize,
                                       {inputDepth, inputHeight, inputWidth},
                                       { deepNet.GetBatchDepth(),deepNet.GetBatchHeight(), deepNet.GetBatchWidth()} ,
                                        deepNet.GetOutputWidth(), nThreads);



      // do an evaluation of the network to compute initial  minimum test error

      Bool_t includeRegularization = (R != DNN::ERegularization::kNone);

      Double_t minValError = 0.0;
      Log() << "Compute initial loss  on the validation data " << Endl;
      for (auto batch : validationData) {
         auto inputTensor = batch.GetInput();
         auto outputMatrix = batch.GetOutput();
         auto weights = batch.GetWeights();

         //std::cout << " input use count " << inputTensor.GetBufferUseCount() << std::endl;
         // should we apply droput to the loss ??
         minValError += deepNet.Loss(inputTensor, outputMatrix, weights, false, includeRegularization);
      }
      // add Regularization term
      Double_t regzTerm = (includeRegularization) ? deepNet.RegularizationTerm() : 0.0;
      minValError /= (Double_t)(nValidationSamples / settings.batchSize);
      minValError += regzTerm;


      // create a pointer to base class VOptimizer
      std::unique_ptr<DNN::VOptimizer<Architecture_t, Layer_t, DeepNet_t>> optimizer;

      // initialize the base class pointer with the corresponding derived class object.
      switch (O) {

      case EOptimizer::kSGD:
         optimizer = std::unique_ptr<DNN::TSGD<Architecture_t, Layer_t, DeepNet_t>>(
            new DNN::TSGD<Architecture_t, Layer_t, DeepNet_t>(settings.learningRate, deepNet, settings.momentum));
         break;

      case EOptimizer::kAdam:
         optimizer = std::unique_ptr<DNN::TAdam<Architecture_t, Layer_t, DeepNet_t>>(
            new DNN::TAdam<Architecture_t, Layer_t, DeepNet_t>(deepNet, settings.learningRate));
         break;

      case EOptimizer::kAdagrad:
         optimizer = std::unique_ptr<DNN::TAdagrad<Architecture_t, Layer_t, DeepNet_t>>(
            new DNN::TAdagrad<Architecture_t, Layer_t, DeepNet_t>(deepNet, settings.learningRate));
         break;

      case EOptimizer::kRMSProp:
         optimizer = std::unique_ptr<DNN::TRMSProp<Architecture_t, Layer_t, DeepNet_t>>(
            new DNN::TRMSProp<Architecture_t, Layer_t, DeepNet_t>(deepNet, settings.learningRate, settings.momentum));
         break;

      case EOptimizer::kAdadelta:
         optimizer = std::unique_ptr<DNN::TAdadelta<Architecture_t, Layer_t, DeepNet_t>>(
            new DNN::TAdadelta<Architecture_t, Layer_t, DeepNet_t>(deepNet, settings.learningRate));
         break;
      }


      // Initialize the vector of batches, one batch for one slave network
      std::vector<TTensorBatch<Architecture_t>> batches{};

      bool converged = false;
      size_t convergenceCount = 0;
      size_t batchesInEpoch = nTrainingSamples / deepNet.GetBatchSize();

      // start measuring
      std::chrono::time_point<std::chrono::system_clock> tstart, tend;
      tstart = std::chrono::system_clock::now();

      Log() << "Training phase " << trainingPhase << " of " << this->GetTrainingSettings().size() << ": "
            << " Optimizer " << settings.optimizerName
            << " Learning rate = " << settings.learningRate
            << " regularization " << (char) settings.regularization
            << " minimum error = " << minValError
            << Endl;
      if (!fInteractive) {
         std::string separator(62, '-');
         Log() << separator << Endl;
         Log() << std::setw(10) << "Epoch"
               << " | " << std::setw(12) << "Train Err." << std::setw(12) << "Val. Err."
               << std::setw(12) << "t(s)/epoch" << std::setw(12)  << "t(s)/Loss"
               << std::setw(12) << "nEvents/s"
               << std::setw(12) << "Conv. Steps" << Endl;
         Log() << separator << Endl;
      }

      // set up generator for shuffling the batches
      // if seed is zero we have always a different order in the batches
      size_t shuffleSeed = 0;
      if (fRandomSeed != 0) shuffleSeed = fRandomSeed + trainingPhase;
      RandomGenerator<TRandom3> rng(shuffleSeed);

      // print weights before
      if (fBuildNet && debug) {
         Log() << "Initial Deep Net Weights " << Endl;
         auto & weights_tensor = deepNet.GetLayerAt(0)->GetWeights();
         for (size_t l = 0; l < weights_tensor.size(); ++l)
            weights_tensor[l].Print();
         auto & bias_tensor = deepNet.GetLayerAt(0)->GetBiases();
         bias_tensor[0].Print();
      }

      Log() << "   Start epoch iteration ..." << Endl;
      bool debugFirstEpoch = false;
      bool computeLossInTraining = true;  // compute loss in training or at test time
      size_t nTrainEpochs = 0;
      while (!converged) {
         nTrainEpochs++;
         trainingData.Shuffle(rng);

         // execute all epochs
         //for (size_t i = 0; i < batchesInEpoch; i += nThreads) {

         Double_t trainingError = 0;
         for (size_t i = 0; i < batchesInEpoch; ++i ) {
            // Clean and load new batches, one batch for one slave net
            //batches.clear();
            //batches.reserve(nThreads);
            //for (size_t j = 0; j < nThreads; j++) {
            //   batches.push_back(trainingData.GetTensorBatch());
            //}
            if (debugFirstEpoch) std::cout << "\n\n----- batch # " << i << "\n\n";

            auto my_batch = trainingData.GetTensorBatch();

            if (debugFirstEpoch)
               std::cout << "got batch data - doing forward \n";

#ifdef DEBUG

            Architecture_t::PrintTensor(my_batch.GetInput(),"input tensor",true);
            typename Architecture_t::Tensor_t tOut(my_batch.GetOutput());
            typename Architecture_t::Tensor_t tW(my_batch.GetWeights());
            Architecture_t::PrintTensor(tOut,"label tensor",true)   ;
            Architecture_t::PrintTensor(tW,"weight tensor",true)  ;
#endif

            deepNet.Forward(my_batch.GetInput(), true);
            // compute also loss
            if (computeLossInTraining) {
               auto outputMatrix = my_batch.GetOutput();
               auto weights = my_batch.GetWeights();
               trainingError += deepNet.Loss(outputMatrix, weights, false);
            }

               if (debugFirstEpoch)
                  std::cout << "- doing backward \n";

#ifdef DEBUG
            size_t nlayers = deepNet.GetLayers().size();
            for (size_t l = 0; l < nlayers; ++l) {
               if (deepNet.GetLayerAt(l)->GetWeights().size() > 0)
                  Architecture_t::PrintTensor(deepNet.GetLayerAt(l)->GetWeightsAt(0),
                                              TString::Format("initial weights layer %d", l).Data());

               Architecture_t::PrintTensor(deepNet.GetLayerAt(l)->GetOutput(),
                                           TString::Format("output tensor layer %d", l).Data());
            }
#endif

            //Architecture_t::PrintTensor(deepNet.GetLayerAt(nlayers-1)->GetOutput(),"output tensor last layer" );

            deepNet.Backward(my_batch.GetInput(), my_batch.GetOutput(), my_batch.GetWeights());

            if (debugFirstEpoch)
               std::cout << "- doing optimizer update  \n";

            // increment optimizer step that is used in some algorithms (e.g. ADAM)
            optimizer->IncrementGlobalStep();
            optimizer->Step();

#ifdef DEBUG
            std::cout << "minmimizer step - momentum " << settings.momentum << " learning rate " << optimizer->GetLearningRate() << std::endl;
            for (size_t l = 0; l < nlayers; ++l) {
               if (deepNet.GetLayerAt(l)->GetWeights().size() > 0) {
                  Architecture_t::PrintTensor(deepNet.GetLayerAt(l)->GetWeightsAt(0),TString::Format("weights after step layer %d",l).Data());
                  Architecture_t::PrintTensor(deepNet.GetLayerAt(l)->GetWeightGradientsAt(0),"weight gradients");
               }
            }
#endif

         }

         if (debugFirstEpoch) std::cout << "\n End batch loop - compute validation loss   \n";
         //}
         debugFirstEpoch = false;
         if ((nTrainEpochs % settings.testInterval) == 0) {

            std::chrono::time_point<std::chrono::system_clock> t1,t2;

            t1 = std::chrono::system_clock::now();

            // Compute validation error.


            Double_t valError = 0.0;
            bool inTraining = false;
            for (auto batch : validationData) {
               auto inputTensor = batch.GetInput();
               auto outputMatrix = batch.GetOutput();
               auto weights = batch.GetWeights();
               // should we apply droput to the loss ??
               valError += deepNet.Loss(inputTensor, outputMatrix, weights, inTraining, includeRegularization);
            }
            // normalize loss to number of batches and add regularization term
            Double_t regTerm = (includeRegularization) ? deepNet.RegularizationTerm() : 0.0;
            valError /= (Double_t)(nValidationSamples / settings.batchSize);
            valError += regTerm;

            //Log the loss value
            fTrainHistory.AddValue("valError",nTrainEpochs,valError);

            t2 = std::chrono::system_clock::now();

            // checking for convergence
            if (valError < minValError) {
               convergenceCount = 0;
            } else {
               convergenceCount += settings.testInterval;
            }

            // copy configuration when reached a minimum error
            if (valError < minValError ) {
               // Copy weights from deepNet to fNet
               Log() << std::setw(10) << nTrainEpochs
                     << " Minimum Test error found - save the configuration " << Endl;
               for (size_t i = 0; i < deepNet.GetDepth(); ++i) {
                  fNet->GetLayerAt(i)->CopyParameters(*deepNet.GetLayerAt(i));
                  // if (i == 0 && deepNet.GetLayerAt(0)->GetWeights().size() > 1) {
                  //    Architecture_t::PrintTensor(deepNet.GetLayerAt(0)->GetWeightsAt(0), " input weights");
                  //    Architecture_t::PrintTensor(deepNet.GetLayerAt(0)->GetWeightsAt(1), " state weights");
                  // }
               }
               // Architecture_t::PrintTensor(deepNet.GetLayerAt(1)->GetWeightsAt(0), " cudnn weights");
               // ArchitectureImpl_t::PrintTensor(fNet->GetLayerAt(1)->GetWeightsAt(0), " cpu weights");

               minValError = valError;
            }
            else if ( minValError <= 0. )
               minValError = valError;

            if (!computeLossInTraining) {
               trainingError = 0.0;
               // Compute training error.
               for (auto batch : trainingData) {
                  auto inputTensor = batch.GetInput();
               auto outputMatrix = batch.GetOutput();
               auto weights = batch.GetWeights();
               trainingError += deepNet.Loss(inputTensor, outputMatrix, weights, false, false);
               }
            }
            // normalize loss to number of batches and add regularization term
            trainingError /= (Double_t)(nTrainingSamples / settings.batchSize);
            trainingError += regTerm;

            //Log the loss value
            fTrainHistory.AddValue("trainingError",nTrainEpochs,trainingError);

            // stop measuring
            tend = std::chrono::system_clock::now();

            // Compute numerical throughput.
            std::chrono::duration<double> elapsed_seconds = tend - tstart;
            std::chrono::duration<double> elapsed1 = t1-tstart;
            // std::chrono::duration<double> elapsed2 = t2-tstart;
            // time to compute training and test errors
            std::chrono::duration<double> elapsed_testing = tend-t1;

            double seconds = elapsed_seconds.count();
            // double nGFlops = (double)(settings.testInterval * batchesInEpoch * settings.batchSize)*1.E-9;
            // nGFlops *= deepnet.GetNFlops() * 1e-9;
            double eventTime = elapsed1.count()/( batchesInEpoch * settings.testInterval * settings.batchSize);

            converged =
               convergenceCount > settings.convergenceSteps || nTrainEpochs >= settings.maxEpochs;


            Log() << std::setw(10) << nTrainEpochs  << " | "
                  << std::setw(12) << trainingError
                  << std::setw(12) << valError
                  << std::setw(12) << seconds / settings.testInterval
                  << std::setw(12)  << elapsed_testing.count()
                  << std::setw(12) << 1. / eventTime
                  << std::setw(12) << convergenceCount
                  << Endl;

            if (converged) {
               Log() << Endl;
            }
            tstart = std::chrono::system_clock::now();
         }

         // if (stepCount % 10 == 0 || converged) {
         if (converged && debug) {
            Log() << "Final Deep Net Weights for phase  " << trainingPhase << " epoch " << nTrainEpochs
                  << Endl;
            auto & weights_tensor = deepNet.GetLayerAt(0)->GetWeights();
            auto & bias_tensor = deepNet.GetLayerAt(0)->GetBiases();
            for (size_t l = 0; l < weights_tensor.size(); ++l)
               weights_tensor[l].Print();
            bias_tensor[0].Print();
         }

      }

      trainingPhase++;
   }  // end loop on training Phase
}

////////////////////////////////////////////////////////////////////////////////
void MethodDL::Train()
{
   if (fInteractive) {
      Log() << kFATAL << "Not implemented yet" << Endl;
      return;
   }

   // using for training same scalar type defined for the prediction
   if (this->GetArchitectureString() == "GPU") {
#ifdef R__HAS_TMVAGPU
      Log() << kINFO << "Start of deep neural network training on GPU." << Endl << Endl;
#ifdef R__HAS_CUDNN
      TrainDeepNet<DNN::TCudnn<ScalarImpl_t> >();
#else
      TrainDeepNet<DNN::TCuda<ScalarImpl_t>>();
#endif
#else
      Log() << kFATAL << "CUDA backend not enabled. Please make sure "
         "you have CUDA installed and it was successfully "
         "detected by CMAKE."
             << Endl;
      return;
#endif
   } else if (this->GetArchitectureString() == "CPU") {
#ifdef R__HAS_TMVACPU
      // note that number of threads used for BLAS might be different
      // e.g use openblas_set_num_threads(num_threads) for OPENBLAS backend
      Log() << kINFO << "Start of deep neural network training on CPU using MT,  nthreads = "
            << gConfig().GetNCpu() << Endl << Endl;
#else
      Log() << kINFO << "Start of deep neural network training on single thread CPU (without ROOT-MT support) " << Endl
            << Endl;
#endif
      TrainDeepNet<DNN::TCpu<ScalarImpl_t> >();
      return;
   }
   else {
      Log() << kFATAL << this->GetArchitectureString() <<
                      " is not  a supported architecture for TMVA::MethodDL"
            << Endl;
   }

}


////////////////////////////////////////////////////////////////////////////////
Double_t MethodDL::GetMvaValue(Double_t * /*errLower*/, Double_t * /*errUpper*/)
{

   // note that fNet  should have been build with a batch size of  1

   if (!fNet || fNet->GetDepth() == 0) {
       Log() << kFATAL << "The network has not been trained and fNet is not built"
             << Endl;
   }

   // input  size must be equal to  1 which is the batch size of fNet
   R__ASSERT(fNet->GetBatchSize() == 1);

   // int batchWidth = fNet->GetBatchWidth();
   // int batchDepth = fNet->GetBatchDepth();
   // int batchHeight = fNet->GetBatchHeight();
//   int noutput = fNet->GetOutputWidth();


   // get current event
   const std::vector<Float_t> &inputValues = GetEvent()->GetValues();

   size_t nVariables = GetEvent()->GetNVariables();

   // for Columnlayout tensor memory layout is   HWC while for rowwise is CHW
   if (fXInput.GetLayout() == TMVA::Experimental::MemoryLayout::ColumnMajor) {
      R__ASSERT(fXInput.GetShape().size() < 4);
      size_t nc, nhw = 0;
      if (fXInput.GetShape().size() == 2) {
         nc =  fXInput.GetShape()[0];
         if (nc != 1 ) {
             ArchitectureImpl_t::PrintTensor(fXInput);
             Log() << kFATAL << "First tensor dimension should be equal to batch size, i.e. = 1"
                   << Endl;
         }
         nhw = fXInput.GetShape()[1];
      } else {
         nc = fXInput.GetCSize();
         nhw = fXInput.GetWSize();
      }
      if ( nVariables != nc * nhw)  {
          Log() << kFATAL << "Input Event variable dimensions are not compatible with the built network architecture"
            << " n-event variables " << nVariables << " expected input tensor " << nc << " x " << nhw
            << Endl;
      }
      for (size_t j = 0; j < nc; j++) {
         for (size_t k = 0; k < nhw; k++) {
            // note that in TMVA events images are stored as C H W while in the buffer we stored as H W C
            fXInputBuffer[ k * nc + j] = inputValues[j*nhw + k];  // for column layout !!!
         }
      }
   } else {
      // row-wise layout
      assert(fXInput.GetShape().size() >= 4);
      size_t nc = fXInput.GetCSize();
      size_t nh = fXInput.GetHSize();
      size_t nw = fXInput.GetWSize();
      size_t n = nc * nh * nw;
      if ( nVariables != n) {
         Log() << kFATAL << "Input Event variable dimensions are not compatible with the built network architecture"
            << " n-event variables " << nVariables << " expected input tensor " << nc << " x " << nh << " x " << nw
            << Endl;
      }
      for (size_t j = 0; j < n; j++) {
         // in this case TMVA event has same order as input tensor
         fXInputBuffer[ j ] = inputValues[j];  // for column layout !!!
      }
   }
   // copy buffer in input
   fXInput.GetDeviceBuffer().CopyFrom( fXInputBuffer);

   // perform the prediction
   fNet->Prediction(*fYHat, fXInput, fOutputFunction);

   // return value
   double mvaValue = (*fYHat)(0, 0);

   // for debugging
#ifdef DEBUG_MVAVALUE
   using Tensor_t = std::vector<MatrixImpl_t>;
    TMatrixF  xInput(n1,n2, inputValues.data() );
    std::cout << "Input data - class " << GetEvent()->GetClass() << std::endl;
    xInput.Print();
    std::cout << "Output of DeepNet " << mvaValue << std::endl;
    auto & deepnet = *fNet;
    std::cout << "Loop on layers " << std::endl;
    for (int l = 0; l < deepnet.GetDepth(); ++l) {
       std::cout << "Layer " << l;
       const auto *  layer = deepnet.GetLayerAt(l);
       const Tensor_t & layer_output = layer->GetOutput();
       layer->Print();
       std::cout << "DNN output " << layer_output.size() << std::endl;
       for (size_t i = 0; i < layer_output.size(); ++i) {
#ifdef R__HAS_TMVAGPU
          //TMatrixD m(layer_output[i].GetNrows(), layer_output[i].GetNcols() , layer_output[i].GetDataPointer()  );
          TMatrixD m = layer_output[i];
#else
          TMatrixD m(layer_output[i].GetNrows(), layer_output[i].GetNcols() , layer_output[i].GetRawDataPointer()  );
#endif
          m.Print();
       }
       const Tensor_t & layer_weights = layer->GetWeights();
       std::cout << "DNN weights " << layer_weights.size() << std::endl;
       if (layer_weights.size() > 0) {
          int i = 0;
#ifdef R__HAS_TMVAGPU
          TMatrixD m = layer_weights[i];
//          TMatrixD m(layer_weights[i].GetNrows(), layer_weights[i].GetNcols() , layer_weights[i].GetDataPointer()  );
#else
          TMatrixD m(layer_weights[i].GetNrows(), layer_weights[i].GetNcols() , layer_weights[i].GetRawDataPointer()  );
#endif
          m.Print();
       }
    }
#endif

   return (TMath::IsNaN(mvaValue)) ? -999. : mvaValue;
}
////////////////////////////////////////////////////////////////////////////////
/// Evaluate the DeepNet on a vector of input values stored in the TMVA Event class
////////////////////////////////////////////////////////////////////////////////
template <typename Architecture_t>
std::vector<Double_t> MethodDL::PredictDeepNet(Long64_t firstEvt, Long64_t lastEvt, size_t batchSize, Bool_t logProgress)
{

   // Check whether the model is setup
   if (!fNet || fNet->GetDepth() == 0) {
       Log() << kFATAL << "The network has not been trained and fNet is not built"
             << Endl;
   }

   // rebuild the networks
   this->SetBatchSize(batchSize);
   size_t inputDepth  = this->GetInputDepth();
   size_t inputHeight = this->GetInputHeight();
   size_t inputWidth  = this->GetInputWidth();
   size_t batchDepth  = this->GetBatchDepth();
   size_t batchHeight = this->GetBatchHeight();
   size_t batchWidth  = this->GetBatchWidth();
   ELossFunction J      = fNet->GetLossFunction();
   EInitialization I    = fNet->GetInitialization();
   ERegularization R    = fNet->GetRegularization();
   Double_t weightDecay = fNet->GetWeightDecay();

   using DeepNet_t          = TMVA::DNN::TDeepNet<Architecture_t>;
   using Matrix_t           = typename Architecture_t::Matrix_t;
   using TensorDataLoader_t = TTensorDataLoader<TMVAInput_t, Architecture_t>;

   // create the deep neural network
   DeepNet_t deepNet(batchSize, inputDepth, inputHeight, inputWidth, batchDepth, batchHeight, batchWidth, J, I, R, weightDecay);
   std::vector<DeepNet_t> nets{};
   fBuildNet = false;
   CreateDeepNet(deepNet,nets);

   // copy weights from the saved fNet to the built DeepNet
   for (size_t i = 0; i < deepNet.GetDepth(); ++i) {
      deepNet.GetLayerAt(i)->CopyParameters(*fNet->GetLayerAt(i));
      // if (i == 0 && deepNet.GetLayerAt(0)->GetWeights().size() > 1) {
      //    Architecture_t::PrintTensor(deepNet.GetLayerAt(0)->GetWeightsAt(0), "Inference: input weights");
      //    Architecture_t::PrintTensor(deepNet.GetLayerAt(0)->GetWeightsAt(1), "Inference: state weights");
      // }
   }

   size_t n1 = deepNet.GetBatchHeight();
   size_t n2 = deepNet.GetBatchWidth();
   size_t n0 = deepNet.GetBatchSize();
   // treat case where batchHeight is the batchSize in case of first Dense layers (then we need to set to fNet batch size)
   if (batchDepth == 1 && GetInputHeight() == 1 && GetInputDepth() == 1) {
      n1 = deepNet.GetBatchSize();
      n0 = 1;
   }
   //this->SetBatchDepth(n0);
   Long64_t nEvents = lastEvt - firstEvt;
   TMVAInput_t testTuple = std::tie(GetEventCollection(Data()->GetCurrentType()), DataInfo());
   TensorDataLoader_t testData(testTuple, nEvents, batchSize, {inputDepth, inputHeight, inputWidth}, {n0, n1, n2}, deepNet.GetOutputWidth(), 1);


   // Tensor_t xInput;
   // for (size_t i = 0; i < n0; ++i)
   //    xInput.emplace_back(Matrix_t(n1,n2));

   // create pointer to output matrix used for the predictions
   Matrix_t yHat(deepNet.GetBatchSize(), deepNet.GetOutputWidth() );

   // use timer
   Timer timer( nEvents, GetName(), kTRUE );

   if (logProgress)
      Log() << kHEADER << Form("[%s] : ",DataInfo().GetName())
            << "Evaluation of " << GetMethodName() << " on "
            << (Data()->GetCurrentType() == Types::kTraining ? "training" : "testing")
            << " sample (" << nEvents << " events)" << Endl;


   // eventg loop
   std::vector<double> mvaValues(nEvents);


   for ( Long64_t ievt = firstEvt;  ievt < lastEvt; ievt+=batchSize) {

      Long64_t ievt_end = ievt + batchSize;
      // case of batch prediction for
      if (ievt_end <=  lastEvt) {

         if (ievt == firstEvt) {
            Data()->SetCurrentEvent(ievt);
            size_t nVariables = GetEvent()->GetNVariables();

            if (n1 == batchSize && n0 == 1)  {
               if (n2 != nVariables) {
                  Log() << kFATAL << "Input Event variable dimensions are not compatible with the built network architecture"
                        << " n-event variables " << nVariables << " expected input matrix " << n1 << " x " << n2
                        << Endl;
               }
            } else {
               if (n1*n2 != nVariables || n0 != batchSize) {
                  Log() << kFATAL << "Input Event variable dimensions are not compatible with the built network architecture"
                        << " n-event variables " << nVariables << " expected input tensor " << n0 << " x " << n1 << " x " << n2
                        << Endl;
               }
            }
         }

         auto batch = testData.GetTensorBatch();
         auto inputTensor = batch.GetInput();

         auto xInput = batch.GetInput();
         // make the prediction
         deepNet.Prediction(yHat, xInput, fOutputFunction);
         for (size_t i = 0; i < batchSize; ++i) {
            double value =  yHat(i,0);
            mvaValues[ievt + i] =  (TMath::IsNaN(value)) ? -999. : value;
         }
      }
      else {
         // case of remaining events: compute prediction by single event !
         for (Long64_t i = ievt; i < lastEvt; ++i) {
            Data()->SetCurrentEvent(i);
            mvaValues[i] = GetMvaValue();
         }
      }
   }

   if (logProgress) {
      Log() << kINFO
            << "Elapsed time for evaluation of " << nEvents <<  " events: "
            << timer.GetElapsedTime() << "       " << Endl;
   }

   return mvaValues;
}

const std::vector<Float_t> & TMVA::MethodDL::GetRegressionValues()
{
   size_t nVariables = GetEvent()->GetNVariables();
   MatrixImpl_t X(1, nVariables);
   TensorImpl_t X_vec ( 1,  1, nVariables);  // needs to be really 1
   const Event *ev = GetEvent();
   const std::vector<Float_t>& inputValues = ev->GetValues();
   for (size_t i = 0; i < nVariables; i++) {
       X_vec(0,i,0) = inputValues[i];   // in case of column format !!
   }
   //X_vec.emplace_back(X);

   size_t nTargets = std::max(1u, ev->GetNTargets());
   MatrixImpl_t YHat(1, nTargets);
   std::vector<Float_t> output(nTargets);
   fNet->Prediction(YHat, X_vec, fOutputFunction);

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

const std::vector<Float_t> & TMVA::MethodDL::GetMulticlassValues()
{
   size_t nVariables = GetEvent()->GetNVariables();
   MatrixImpl_t X(1, nVariables);
   TensorImpl_t X_vec ( 1, 1, nVariables);
   MatrixImpl_t YHat(1, DataInfo().GetNClasses());
   if (fMulticlassReturnVal == NULL) {
      fMulticlassReturnVal = new std::vector<Float_t>(DataInfo().GetNClasses());
   }

   const std::vector<Float_t>& inputValues = GetEvent()->GetValues();
   for (size_t i = 0; i < nVariables; i++) {
      X_vec(0,i, 0) = inputValues[i];
   }
   //X_vec.emplace_back(X);
   fNet->Prediction(YHat, X_vec, fOutputFunction);
   for (size_t i = 0; i < (size_t) YHat.GetNcols(); i++) {
      (*fMulticlassReturnVal)[i] = YHat(0, i);
   }
   return *fMulticlassReturnVal;
}


////////////////////////////////////////////////////////////////////////////////
/// Evaluate the DeepNet on a vector of input values stored in the TMVA Event class
////////////////////////////////////////////////////////////////////////////////
std::vector<Double_t> MethodDL::GetMvaValues(Long64_t firstEvt, Long64_t lastEvt, Bool_t logProgress)
{
   // Long64_t nEvents = Data()->GetNEvents();
   // std::vector<Double_t>  v(nEvents);
   // for (Long64_t i = 0; i < nEvents; ++i) {
   //    Data()->SetCurrentEvent(i);
   //    v[i] = GetMvaValue();
   // }
   // return v;


   Long64_t nEvents = Data()->GetNEvents();
   if (firstEvt > lastEvt || lastEvt > nEvents) lastEvt = nEvents;
   if (firstEvt < 0) firstEvt = 0;
   nEvents = lastEvt-firstEvt;

   // use same batch size as for training (from first strategy)
   size_t defaultEvalBatchSize = (fXInput.GetSize() > 1000) ? 100 : 1000;
   size_t batchSize = (fTrainingSettings.empty()) ? defaultEvalBatchSize :  fTrainingSettings.front().batchSize;
   if  ( size_t(nEvents) < batchSize ) batchSize = nEvents;

   // using for training same scalar type defined for the prediction
   if (this->GetArchitectureString() == "GPU") {
#ifdef R__HAS_TMVAGPU
      Log() << kINFO << "Evaluate deep neural network on GPU using batches with size = " <<  batchSize << Endl << Endl;
#ifdef R__HAS_CUDNN
      return PredictDeepNet<DNN::TCudnn<ScalarImpl_t>>(firstEvt, lastEvt, batchSize, logProgress);
#else
      return PredictDeepNet<DNN::TCuda<ScalarImpl_t>>(firstEvt, lastEvt, batchSize, logProgress);
#endif

#endif
   }
   Log() << kINFO << "Evaluate deep neural network on CPU using batches with size = " << batchSize << Endl << Endl;
   return PredictDeepNet<DNN::TCpu<ScalarImpl_t> >(firstEvt, lastEvt, batchSize, logProgress);
}
////////////////////////////////////////////////////////////////////////////////
void MethodDL::AddWeightsXMLTo(void * parent) const
{
      // Create the parent XML node with name "Weights"
   auto & xmlEngine = gTools().xmlengine();
   void* nn = xmlEngine.NewChild(parent, 0, "Weights");

   /*! Get all necessary information, in order to be able to reconstruct the net
    *  if we read the same XML file. */

   // Deep Net specific info
   Int_t depth = fNet->GetDepth();

   Int_t inputDepth = fNet->GetInputDepth();
   Int_t inputHeight = fNet->GetInputHeight();
   Int_t inputWidth = fNet->GetInputWidth();

   Int_t batchSize = fNet->GetBatchSize();

   Int_t batchDepth = fNet->GetBatchDepth();
   Int_t batchHeight = fNet->GetBatchHeight();
   Int_t batchWidth = fNet->GetBatchWidth();

   char lossFunction = static_cast<char>(fNet->GetLossFunction());
   char initialization = static_cast<char>(fNet->GetInitialization());
   char regularization = static_cast<char>(fNet->GetRegularization());

   Double_t weightDecay = fNet->GetWeightDecay();

   // Method specific info (not sure these are needed)
   char outputFunction = static_cast<char>(this->GetOutputFunction());
   //char lossFunction = static_cast<char>(this->GetLossFunction());

   // Add attributes to the parent node
   xmlEngine.NewAttr(nn, 0, "NetDepth", gTools().StringFromInt(depth));

   xmlEngine.NewAttr(nn, 0, "InputDepth", gTools().StringFromInt(inputDepth));
   xmlEngine.NewAttr(nn, 0, "InputHeight", gTools().StringFromInt(inputHeight));
   xmlEngine.NewAttr(nn, 0, "InputWidth", gTools().StringFromInt(inputWidth));

   xmlEngine.NewAttr(nn, 0, "BatchSize", gTools().StringFromInt(batchSize));
   xmlEngine.NewAttr(nn, 0, "BatchDepth", gTools().StringFromInt(batchDepth));
   xmlEngine.NewAttr(nn, 0, "BatchHeight", gTools().StringFromInt(batchHeight));
   xmlEngine.NewAttr(nn, 0, "BatchWidth", gTools().StringFromInt(batchWidth));

   xmlEngine.NewAttr(nn, 0, "LossFunction", TString(lossFunction));
   xmlEngine.NewAttr(nn, 0, "Initialization", TString(initialization));
   xmlEngine.NewAttr(nn, 0, "Regularization", TString(regularization));
   xmlEngine.NewAttr(nn, 0, "OutputFunction", TString(outputFunction));

   gTools().AddAttr(nn, "WeightDecay", weightDecay);


   for (Int_t i = 0; i < depth; i++)
   {
      fNet->GetLayerAt(i) -> AddWeightsXMLTo(nn);
   }


}

////////////////////////////////////////////////////////////////////////////////
void MethodDL::ReadWeightsFromXML(void * rootXML)
{

   auto netXML = gTools().GetChild(rootXML, "Weights");
   if (!netXML){
      netXML = rootXML;
   }

   size_t netDepth;
   gTools().ReadAttr(netXML, "NetDepth", netDepth);

   size_t inputDepth, inputHeight, inputWidth;
   gTools().ReadAttr(netXML, "InputDepth", inputDepth);
   gTools().ReadAttr(netXML, "InputHeight", inputHeight);
   gTools().ReadAttr(netXML, "InputWidth", inputWidth);

   size_t batchSize, batchDepth, batchHeight, batchWidth;
   gTools().ReadAttr(netXML, "BatchSize", batchSize);
   // use always batchsize = 1
   //batchSize = 1;
   gTools().ReadAttr(netXML, "BatchDepth", batchDepth);
   gTools().ReadAttr(netXML, "BatchHeight", batchHeight);
   gTools().ReadAttr(netXML, "BatchWidth",  batchWidth);

   char lossFunctionChar;
   gTools().ReadAttr(netXML, "LossFunction", lossFunctionChar);
   char initializationChar;
   gTools().ReadAttr(netXML, "Initialization", initializationChar);
   char regularizationChar;
   gTools().ReadAttr(netXML, "Regularization", regularizationChar);
   char outputFunctionChar;
   gTools().ReadAttr(netXML, "OutputFunction", outputFunctionChar);
   double weightDecay;
   gTools().ReadAttr(netXML, "WeightDecay", weightDecay);

   // create the net

   // DeepNetCpu_t is defined in MethodDL.h
   this->SetInputDepth(inputDepth);
   this->SetInputHeight(inputHeight);
   this->SetInputWidth(inputWidth);
   this->SetBatchDepth(batchDepth);
   this->SetBatchHeight(batchHeight);
   this->SetBatchWidth(batchWidth);



   fNet = std::unique_ptr<DeepNetImpl_t>(new DeepNetImpl_t(batchSize, inputDepth, inputHeight, inputWidth, batchDepth,
                                                   batchHeight, batchWidth,
                                                   static_cast<ELossFunction>(lossFunctionChar),
                                                   static_cast<EInitialization>(initializationChar),
                                                   static_cast<ERegularization>(regularizationChar),
                                                   weightDecay));

   fOutputFunction = static_cast<EOutputFunction>(outputFunctionChar);


   //size_t previousWidth = inputWidth;
   auto layerXML = gTools().xmlengine().GetChild(netXML);

   // loop on the layer and add them to the network
   for (size_t i = 0; i < netDepth; i++) {

      TString layerName = gTools().xmlengine().GetNodeName(layerXML);

      // case of dense layer
      if (layerName == "DenseLayer") {

         // read width and activation function and then we can create the layer
         size_t width = 0;
         gTools().ReadAttr(layerXML, "Width", width);

         // Read activation function.
         TString funcString;
         gTools().ReadAttr(layerXML, "ActivationFunction", funcString);
         EActivationFunction func = static_cast<EActivationFunction>(funcString.Atoi());


         fNet->AddDenseLayer(width, func, 0.0); // no need to pass dropout probability

      }
      // Convolutional Layer
      else if (layerName == "ConvLayer") {

         // read width and activation function and then we can create the layer
         size_t depth = 0;
         gTools().ReadAttr(layerXML, "Depth", depth);
         size_t fltHeight, fltWidth = 0;
         size_t strideRows, strideCols = 0;
         size_t padHeight, padWidth = 0;
         gTools().ReadAttr(layerXML, "FilterHeight", fltHeight);
         gTools().ReadAttr(layerXML, "FilterWidth", fltWidth);
         gTools().ReadAttr(layerXML, "StrideRows", strideRows);
         gTools().ReadAttr(layerXML, "StrideCols", strideCols);
         gTools().ReadAttr(layerXML, "PaddingHeight", padHeight);
         gTools().ReadAttr(layerXML, "PaddingWidth", padWidth);

         // Read activation function.
         TString funcString;
         gTools().ReadAttr(layerXML, "ActivationFunction", funcString);
         EActivationFunction actFunction = static_cast<EActivationFunction>(funcString.Atoi());


         fNet->AddConvLayer(depth, fltHeight, fltWidth, strideRows, strideCols,
                            padHeight, padWidth, actFunction);

      }

      // MaxPool Layer
      else if (layerName == "MaxPoolLayer") {

         // read maxpool layer info
         size_t filterHeight, filterWidth = 0;
         size_t strideRows, strideCols = 0;
         gTools().ReadAttr(layerXML, "FilterHeight", filterHeight);
         gTools().ReadAttr(layerXML, "FilterWidth", filterWidth);
         gTools().ReadAttr(layerXML, "StrideRows", strideRows);
         gTools().ReadAttr(layerXML, "StrideCols", strideCols);

         fNet->AddMaxPoolLayer(filterHeight, filterWidth, strideRows, strideCols);
      }
      // Reshape Layer
      else if (layerName == "ReshapeLayer") {

         // read reshape layer info
         size_t depth, height, width = 0;
         gTools().ReadAttr(layerXML, "Depth", depth);
         gTools().ReadAttr(layerXML, "Height", height);
         gTools().ReadAttr(layerXML, "Width", width);
         int flattening = 0;
         gTools().ReadAttr(layerXML, "Flattening",flattening );

         fNet->AddReshapeLayer(depth, height, width, flattening);

      }
      // RNN Layer
      else if (layerName == "RNNLayer") {

         // read RNN layer info
         size_t  stateSize,inputSize, timeSteps = 0;
         int rememberState= 0;
         int returnSequence = 0;
         gTools().ReadAttr(layerXML, "StateSize", stateSize);
         gTools().ReadAttr(layerXML, "InputSize", inputSize);
         gTools().ReadAttr(layerXML, "TimeSteps", timeSteps);
         gTools().ReadAttr(layerXML, "RememberState", rememberState );
         gTools().ReadAttr(layerXML, "ReturnSequence", returnSequence);

         fNet->AddBasicRNNLayer(stateSize, inputSize, timeSteps, rememberState, returnSequence);

      }
      // LSTM Layer
      else if (layerName == "LSTMLayer") {

         // read RNN layer info
         size_t  stateSize,inputSize, timeSteps = 0;
         int rememberState, returnSequence = 0;
         gTools().ReadAttr(layerXML, "StateSize", stateSize);
         gTools().ReadAttr(layerXML, "InputSize", inputSize);
         gTools().ReadAttr(layerXML, "TimeSteps", timeSteps);
         gTools().ReadAttr(layerXML, "RememberState", rememberState );
         gTools().ReadAttr(layerXML, "ReturnSequence", returnSequence);

         fNet->AddBasicLSTMLayer(stateSize, inputSize, timeSteps, rememberState, returnSequence);

      }
      // GRU Layer
      else if (layerName == "GRULayer") {

         // read RNN layer info
         size_t  stateSize,inputSize, timeSteps = 0;
         int rememberState, returnSequence, resetGateAfter = 0;
         gTools().ReadAttr(layerXML, "StateSize", stateSize);
         gTools().ReadAttr(layerXML, "InputSize", inputSize);
         gTools().ReadAttr(layerXML, "TimeSteps", timeSteps);
         gTools().ReadAttr(layerXML, "RememberState", rememberState );
         gTools().ReadAttr(layerXML, "ReturnSequence", returnSequence);
         gTools().ReadAttr(layerXML, "ResetGateAfter", resetGateAfter);

         if (!resetGateAfter && ArchitectureImpl_t::IsCudnn())
            Warning("ReadWeightsFromXML",
                    "Cannot use a reset gate after to false with CudNN - use implementation with resetgate=true");

         fNet->AddBasicGRULayer(stateSize, inputSize, timeSteps, rememberState, returnSequence, resetGateAfter);
      }
      // BatchNorm Layer
      else if (layerName == "BatchNormLayer") {
         // use some dammy value which will be overwrittem in BatchNormLayer::ReadWeightsFromXML
         fNet->AddBatchNormLayer(0., 0.0);
      }
      // read weights and biases
      fNet->GetLayers().back()->ReadWeightsFromXML(layerXML);

      // read next layer
      layerXML = gTools().GetNextChild(layerXML);
   }

   fBuildNet = false;
   // create now the input and output matrices
   //int n1 = batchHeight;
   //int n2 = batchWidth;
   // treat case where batchHeight is the batchSize in case of first Dense layers (then we need to set to fNet batch size)
   //if (fXInput.size() > 0) fXInput.clear();
   //fXInput.emplace_back(MatrixImpl_t(n1,n2));
   fXInput = ArchitectureImpl_t::CreateTensor(fNet->GetBatchSize(), GetInputDepth(), GetInputHeight(), GetInputWidth() );
   if (batchDepth == 1 && GetInputHeight() == 1 && GetInputDepth() == 1)
      // make here a ColumnMajor tensor
      fXInput = TensorImpl_t( fNet->GetBatchSize(), GetInputWidth(),TMVA::Experimental::MemoryLayout::ColumnMajor );
   fXInputBuffer =  HostBufferImpl_t( fXInput.GetSize());

   // create pointer to output matrix used for the predictions
   fYHat = std::unique_ptr<MatrixImpl_t>(new MatrixImpl_t(fNet->GetBatchSize(),  fNet->GetOutputWidth() ) );


}


////////////////////////////////////////////////////////////////////////////////
void MethodDL::ReadWeightsFromStream(std::istream & /*istr*/)
{
}

////////////////////////////////////////////////////////////////////////////////
const Ranking *TMVA::MethodDL::CreateRanking()
{
   // TODO
   return NULL;
}

////////////////////////////////////////////////////////////////////////////////
void MethodDL::GetHelpMessage() const
{
   // TODO
}

} // namespace TMVA
