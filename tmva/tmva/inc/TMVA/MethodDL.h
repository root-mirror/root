// @(#)root/tmva/tmva/dnn:$Id$
// Author: Vladimir Ilievski, Saurav Shekhar

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

#ifndef ROOT_TMVA_MethodDL
#define ROOT_TMVA_MethodDL

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// MethodDL                                                             //
//                                                                      //
// Method class for all Deep Learning Networks                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TString.h"

#include "TMVA/MethodBase.h"
#include "TMVA/Types.h"

#include "TMVA/DNN/Architectures/Reference.h"

//#ifdef R__HAS_TMVACPU
#include "TMVA/DNN/Architectures/Cpu.h"
//#endif

#if 0
#ifdef R__HAS_TMVAGPU
#include "TMVA/DNN/Architectures/Cuda.h"
#ifdef R__HAS_CUDNN
#include "TMVA/DNN/Architectures/TCudnn.h"
#endif
#endif
#endif

#include "TMVA/DNN/Functions.h"
#include "TMVA/DNN/DeepNet.h"

#include <vector>
#include <map>

#ifdef R__HAS_TMVAGPU
//#define USE_GPU_INFERENCE
#endif

namespace TMVA {

/*! All of the options that can be specified in the training string */
struct TTrainingSettings {
   size_t batchSize;
   size_t testInterval;
   size_t convergenceSteps;
   size_t maxEpochs;
   DNN::ERegularization regularization;
   DNN::EOptimizer optimizer;
   TString optimizerName;
   Double_t learningRate;
   Double_t momentum;
   Double_t weightDecay;
   std::vector<Double_t> dropoutProbabilities;
   bool multithreading;
};


class MethodDL : public MethodBase {

private:
   // Key-Value vector type, contining the values for the training options
   using KeyValueVector_t = std::vector<std::map<TString, TString>>;

// #ifdef R__HAS_TMVAGPU
// #ifdef R__HAS_CUDNN
//    using ArchitectureImpl_t = TMVA::DNN::TCudnn<Float_t>;
// #else
//   using ArchitectureImpl_t = TMVA::DNN::TCuda<Float_t>;
// #endif
// #else
// do not use arch GPU for evaluation. It is too slow for batch size=1
   using ArchitectureImpl_t = TMVA::DNN::TCpu<Float_t>;
// #endif

   using DeepNetImpl_t = TMVA::DNN::TDeepNet<ArchitectureImpl_t>;
   using MatrixImpl_t =  typename ArchitectureImpl_t::Matrix_t;
   using TensorImpl_t =  typename ArchitectureImpl_t::Tensor_t;
   using ScalarImpl_t =  typename ArchitectureImpl_t::Scalar_t;
   using HostBufferImpl_t = typename ArchitectureImpl_t::HostBuffer_t;

   /*! The option handling methods */
   void DeclareOptions();
   void ProcessOptions();

   void Init();

   // Function to parse the layout of the input
   void ParseInputLayout();
   void ParseBatchLayout();

   /*! After calling the ProcesOptions(), all of the options are parsed,
    *  so using the parsed options, and given the architecture and the
    *  type of the layers, we build the Deep Network passed as
    *  a reference in the function. */
   template <typename Architecture_t, typename Layer_t>
   void CreateDeepNet(DNN::TDeepNet<Architecture_t, Layer_t> &deepNet,
                      std::vector<DNN::TDeepNet<Architecture_t, Layer_t>> &nets);

   template <typename Architecture_t, typename Layer_t>
   void ParseDenseLayer(DNN::TDeepNet<Architecture_t, Layer_t> &deepNet,
                        std::vector<DNN::TDeepNet<Architecture_t, Layer_t>> &nets, TString layerString, TString delim);

   template <typename Architecture_t, typename Layer_t>
   void ParseConvLayer(DNN::TDeepNet<Architecture_t, Layer_t> &deepNet,
                       std::vector<DNN::TDeepNet<Architecture_t, Layer_t>> &nets, TString layerString, TString delim);

   template <typename Architecture_t, typename Layer_t>
   void ParseMaxPoolLayer(DNN::TDeepNet<Architecture_t, Layer_t> &deepNet,
                          std::vector<DNN::TDeepNet<Architecture_t, Layer_t>> &nets, TString layerString,
                          TString delim);

   template <typename Architecture_t, typename Layer_t>
   void ParseReshapeLayer(DNN::TDeepNet<Architecture_t, Layer_t> &deepNet,
                          std::vector<DNN::TDeepNet<Architecture_t, Layer_t>> &nets, TString layerString,
                          TString delim);

   template <typename Architecture_t, typename Layer_t>
   void ParseBatchNormLayer(DNN::TDeepNet<Architecture_t, Layer_t> &deepNet,
                          std::vector<DNN::TDeepNet<Architecture_t, Layer_t>> &nets, TString layerString,
                          TString delim);

   enum ERecurrentLayerType { kLayerRNN = 0, kLayerLSTM = 1, kLayerGRU = 2 };
   template <typename Architecture_t, typename Layer_t>
   void ParseRecurrentLayer(ERecurrentLayerType type, DNN::TDeepNet<Architecture_t, Layer_t> &deepNet,
                       std::vector<DNN::TDeepNet<Architecture_t, Layer_t>> &nets, TString layerString, TString delim);


   /// train of deep neural network using the defined architecture
   template <typename Architecture_t>
   void TrainDeepNet();

   /// perform prediction of the deep neural network
   /// using batches (called by GetMvaValues)
   template <typename Architecture_t>
   std::vector<Double_t> PredictDeepNet(Long64_t firstEvt, Long64_t lastEvt, size_t batchSize, Bool_t logProgress);

   /// parce the validation string and return the number of event data used for validation
   UInt_t GetNumValidationSamples();

   // cudnn implementation needs this format
   /** Contains the batch size (no. of images in the batch), input depth (no. channels)
    *  and furhter input dimensios of the data (image height, width ...)*/
   std::vector<size_t> fInputShape;

   // The size of the batch, i.e. the number of images that are contained in the batch, is either set to be the depth
   // or the height of the batch
   size_t fBatchDepth;  ///< The depth of the batch used to train the deep net.
   size_t fBatchHeight; ///< The height of the batch used to train the deep net.
   size_t fBatchWidth;  ///< The width of the batch used to train the deep net.

   size_t fRandomSeed;  ///<The random seed used to initialize the weights and shuffling batches (default is zero)

   DNN::EInitialization fWeightInitialization; ///< The initialization method
   DNN::EOutputFunction fOutputFunction;       ///< The output function for making the predictions
   DNN::ELossFunction   fLossFunction;         ///< The loss function

   TString fInputLayoutString;          ///< The string defining the layout of the input
   TString fBatchLayoutString;          ///< The string defining the layout of the batch
   TString fLayoutString;               ///< The string defining the layout of the deep net
   TString fErrorStrategy;              ///< The string defining the error strategy for training
   TString fTrainingStrategyString;     ///< The string defining the training strategy
   TString fWeightInitializationString; ///< The string defining the weight initialization method
   TString fArchitectureString;         ///< The string defining the architecure: CPU or GPU
   TString fNumValidationString;        ///< The string defining the number (or percentage) of training data used for validation
   bool fResume;
   bool fBuildNet;                     ///< Flag to control whether to build fNet, the stored network used for the evaluation

   KeyValueVector_t fSettings;                       ///< Map for the training strategy
   std::vector<TTrainingSettings> fTrainingSettings; ///< The vector defining each training strategy

   TensorImpl_t fXInput;                 // input tensor used to evaluate fNet
   HostBufferImpl_t fXInputBuffer;        // input hist buffer corresponding to X (needed for GPU implementation)
   std::unique_ptr<MatrixImpl_t> fYHat;   // output prediction matrix of fNet
   std::unique_ptr<DeepNetImpl_t> fNet;


   ClassDef(MethodDL, 0);

protected:
   // provide a help message
   void GetHelpMessage() const;

   virtual std::vector<Double_t> GetMvaValues(Long64_t firstEvt, Long64_t lastEvt, Bool_t logProgress);


public:
   /*! Constructor */
   MethodDL(const TString &jobName, const TString &methodTitle, DataSetInfo &theData, const TString &theOption);

   /*! Constructor */
   MethodDL(DataSetInfo &theData, const TString &theWeightFile);

   /*! Virtual Destructor */
   virtual ~MethodDL();

   /*! Function for parsing the training settings, provided as a string
    *  in a key-value form.  */
   KeyValueVector_t ParseKeyValueString(TString parseString, TString blockDelim, TString tokenDelim);

   /*! Check the type of analysis the deep learning network can do */
   Bool_t HasAnalysisType(Types::EAnalysisType type, UInt_t numberClasses, UInt_t numberTargets);

   /*! Methods for training the deep learning network */
   void Train();

   Double_t GetMvaValue(Double_t *err = 0, Double_t *errUpper = 0);
   virtual const std::vector<Float_t>& GetRegressionValues();
   virtual const std::vector<Float_t>& GetMulticlassValues();

   /*! Methods for writing and reading weights */
   using MethodBase::ReadWeightsFromStream;
   void AddWeightsXMLTo(void *parent) const;
   void ReadWeightsFromXML(void *wghtnode);
   void ReadWeightsFromStream(std::istream &);

   /* Create ranking */
   const Ranking *CreateRanking();

   /* Getters */
   size_t GetInputDepth()  const { return fInputShape[1]; }   //< no. of channels for an image
   size_t GetInputHeight() const { return fInputShape[2]; }
   size_t GetInputWidth()  const { return fInputShape[3]; }
   size_t GetInputDim()    const { return fInputShape.size() - 2; }
   std::vector<size_t> GetInputShape() const { return fInputShape; }

   size_t GetBatchSize()   const { return fInputShape[0]; }
   size_t GetBatchDepth()  const { return fBatchDepth; }
   size_t GetBatchHeight() const { return fBatchHeight; }
   size_t GetBatchWidth()  const { return fBatchWidth; }

   const DeepNetImpl_t & GetDeepNet() const { return *fNet; }

   DNN::EInitialization GetWeightInitialization() const { return fWeightInitialization; }
   DNN::EOutputFunction GetOutputFunction()       const { return fOutputFunction; }
   DNN::ELossFunction GetLossFunction()           const { return fLossFunction; }

   TString GetInputLayoutString()          const { return fInputLayoutString; }
   TString GetBatchLayoutString()          const { return fBatchLayoutString; }
   TString GetLayoutString()               const { return fLayoutString; }
   TString GetErrorStrategyString()        const { return fErrorStrategy; }
   TString GetTrainingStrategyString()     const { return fTrainingStrategyString; }
   TString GetWeightInitializationString() const { return fWeightInitializationString; }
   TString GetArchitectureString()         const { return fArchitectureString; }

   const std::vector<TTrainingSettings> &GetTrainingSettings() const { return fTrainingSettings; }
   std::vector<TTrainingSettings>       &GetTrainingSettings()       { return fTrainingSettings; }
   const KeyValueVector_t               &GetKeyValueSettings() const { return fSettings; }
   KeyValueVector_t                     &GetKeyValueSettings()       { return fSettings; }

   /** Setters */
   void SetInputDepth (int inputDepth)  { fInputShape[1] = inputDepth; }
   void SetInputHeight(int inputHeight) { fInputShape[2] = inputHeight; }
   void SetInputWidth (int inputWidth)  { fInputShape[3] = inputWidth; }
   void SetInputShape (std::vector<size_t> inputShape) { fInputShape = std::move(inputShape); }

   void SetBatchSize  (size_t batchSize)   { fInputShape[0] = batchSize; }
   void SetBatchDepth (size_t batchDepth)  { fBatchDepth    = batchDepth; }
   void SetBatchHeight(size_t batchHeight) { fBatchHeight   = batchHeight; }
   void SetBatchWidth (size_t batchWidth)  { fBatchWidth    = batchWidth; }

   void SetWeightInitialization(DNN::EInitialization weightInitialization)
   {
      fWeightInitialization = weightInitialization;
   }
   void SetOutputFunction            (DNN::EOutputFunction outputFunction) { fOutputFunction = outputFunction; }
   void SetErrorStrategyString       (TString errorStrategy)               { fErrorStrategy = errorStrategy; }
   void SetTrainingStrategyString    (TString trainingStrategyString)      { fTrainingStrategyString = trainingStrategyString; }
   void SetWeightInitializationString(TString weightInitializationString)
   {
      fWeightInitializationString = weightInitializationString;
   }
   void SetArchitectureString        (TString architectureString)          { fArchitectureString = architectureString; }
   void SetLayoutString              (TString layoutString)                { fLayoutString = layoutString; }
};

} // namespace TMVA

#endif
