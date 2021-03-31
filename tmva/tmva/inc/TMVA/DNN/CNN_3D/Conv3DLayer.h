// @(#)root/tmva/tmva/dnn:$Id$
// Author: Surya S Dwivedi

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : TConv3DLayer                                                            *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      3D Convolutional Deep Neural Network Layer                                *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Surya S Dwivedi  <surya1997@utexas.edu> - Univ of Texas Austin            *
 *                                                                                *                                                  
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/


#ifndef TMVA_CNN_3D_CONV3DLAYER
#define TMVA_CNN_3D_CONV3DLAYER


#include "TMatrix.h"

#include "TMVA/DNN/GeneralLayer.h"
#include "TMVA/DNN/Functions.h"
#include "TMVA/DNN/CNN_3D/ContextHandles.h"

#include <vector>
#include <iostream>

namespace TMVA {
namespace DNN {
namespace CNN_3D {


typedef struct TConv3DParams {

public:
   size_t batchSize; ///< Batch size used for training and evaluation

   size_t inputDepth;  ///< The depth of the previous layer or input.
   size_t inputHeight; ///< The height of the previous layer or input.
   size_t inputWidth;  ///< The width of the previous layer or input.

   size_t input4D;     ///< The depth along the 4th dimension (or noOfmaps)


   size_t outputDepth;  ///< The depth of the previous layer or input ; GetDepth() will return this
   size_t outputHeight; ///< The height of the previous layer or input
   size_t outputWidth;  ///< The width of the previous layer or input
   size_t output4D;     ///< The depth along the 4th dimension, same as numFilters
   
   size_t numberFilters; ///< The number of the filters, which is equal to the output4D.
   size_t filterHeight;  ///< The height of the filter.
   size_t filterWidth;   ///< The width of the filter.
   size_t filterDepth;   ///< The depth of the filter

   size_t strideX;    ///< The number of row pixels to slid the filter each step.
   size_t strideY;    ///< The number of column pixels to slid the filter each step.
   size_t strideZ;

   size_t paddingHeight; ///< The number of zero layers added top and bottom of the input.
   size_t paddingWidth;  ///< The number of zero layers left and right of the input.
   size_t paddingDepth;  ///< The number of zero layers front and back of the input.


   size_t calculateDimension(size_t imgDim, size_t fltDim, size_t padding, size_t stride)
   {
      size_t temp = imgDim - fltDim + 2 * padding;
      if (temp % stride || temp + stride <= 0) {
        Fatal("calculateDimension", "Not compatible hyper parameters for layer - (imageDim, filterDim, padding, stride) "
        "%zu, %zu, %zu, %zu", imgDim, fltDim, padding, stride);
      }
      return temp / stride + 1;
   }


   TConv3DParams(size_t _batchSize, size_t _inputDepth, size_t _inputHeight, size_t _inputWidth, size_t _input4D, size_t _output4D,
               size_t _filterHeight, size_t _filterWidth, size_t _filterDepth, size_t _strideX, size_t _strideY, size_t _strideZ,
               size_t _paddingHeight, size_t _paddingWidth, size_t _paddingDepth)
             : batchSize(_batchSize), inputDepth(_inputDepth), inputHeight(_inputHeight), inputWidth(_inputWidth), input4D(_input4D),
               outputDepth(calculateDimension(_inputDepth, _filterDepth, _paddingDepth, _strideZ)), 
               outputHeight(calculateDimension(_inputHeight, _filterHeight, _paddingHeight, _strideX)), outputWidth(calculateDimension(_inputWidth, _filterWidth, _paddingWidth, _strideY)), 
               output4D(_output4D), numberFilters(_output4D), filterHeight(_filterHeight), filterWidth(_filterWidth), filterDepth(_filterDepth), 
               strideX(_strideX), strideY(_strideY), strideZ(_strideZ), paddingHeight(_paddingHeight),
               paddingWidth(_paddingWidth), paddingDepth(_paddingDepth)
   {}
} TConv3DParams;



template <typename Architecture_t>
class TConv3DLayer : public VGeneralLayer<Architecture_t> {
public:
   using Tensor_t = typename Architecture_t::Tensor_t;
   using Matrix_t = typename Architecture_t::Matrix_t;
   using Scalar_t = typename Architecture_t::Scalar_t;

   using LayerDescriptor_t   = typename Architecture_t::ConvolutionDescriptor_t;
   using WeightsDescriptor_t = typename Architecture_t::FilterDescriptor_t;
   using HelperDescriptor_t  = typename Architecture_t::ActivationDescriptor_t;

   using AlgorithmForward_t  = typename Architecture_t::AlgorithmForward_t;  // Forward layer operation
   using AlgorithmBackward_t = typename Architecture_t::AlgorithmBackward_t; // Backward layer operation
   using AlgorithmHelper_t   = typename Architecture_t::AlgorithmHelper_t;   // Used for weight grad backward pass
   using ReduceTensorDescriptor_t = typename Architecture_t::ReduceTensorDescriptor_t; // used for reduction of tensor(bias grad)

   // FIXME: Add other cudnn types (algorithm preference etc.)
   using AlgorithmDataType_t = typename Architecture_t::AlgorithmDataType_t;

   /* Calculate the output dimension of the convolutional layer */
   static size_t calculateDimension(size_t imgDim, size_t fltDim, size_t padding, size_t stride);

   /* Calculate the number of pixels in a single receptive field */
   static size_t inline calculateNLocalViewPixels(size_t depth, size_t height, size_t width, size_t depth4D) { return depth * height * width * depth4D; }

   /* Calculate the number of receptive fields in an image given the filter and image sizes */
   static size_t calculateNLocalViews(size_t inputHeight, size_t filterHeight, size_t paddingHeight, size_t strideX,
                               size_t inputWidth, size_t filterWidth, size_t paddingWidth, size_t strideY,
                               size_t inputDepth, size_t filterDepth, size_t paddingDepth, size_t strideZ);

protected:
   size_t finput4D;
   size_t foutput4D;
   size_t fFilterDepth;  ///< The depth of the filter.
   size_t fFilterHeight; ///< The height of the filter.
   size_t fFilterWidth;  ///< The width of the filter.

   size_t fStrideX;   ///< The number of row pixels to slid the filter each step.
   size_t fStrideY;   ///< The number of column pixels to slid the filter each step.
   size_t fStrideZ;   ///< The number of depth pixels to slid filter each step

   size_t fNLocalViewPixels;     ///< The number of pixels in one local image view.
   size_t fNLocalViews;          ///< The number of local views in one image.

   Scalar_t fDropoutProbability; ///< Probability that an input is active.

   TDescriptors * fDescriptors = nullptr;  ///< Keeps the convolution, activations and filter descriptors

   TWorkspace * fWorkspace = nullptr;
private:
   size_t fPaddingHeight;        ///< The number of zero layers added top and bottom of the input.
   size_t fPaddingWidth;         ///< The number of zero layers left and right of the input.

   size_t fPaddingDepth;          ///< The number of zero layers front and back of the input

   Tensor_t fInputActivation;        ///< First output of this layer after conv, before activation.

   std::vector<int> fBackwardIndices;  ///< Vector of indices used for a fast Im2Col in backward pass

   EActivationFunction fF;             ///< Activation function of the layer.
   ERegularization fReg;               ///< The regularization method.
   Scalar_t fWeightDecay;              ///< The weight decay.

   Tensor_t fForwardTensor;            ///< Cache tensor used for speeding-up the forward pass.

   void InitializeDescriptors();
   void ReleaseDescriptors();
   void InitializeWorkspace();
   void FreeWorkspace();

public:
   /*! Constructor. */
   TConv3DLayer(size_t BatchSize, size_t InputDepth, size_t InputHeight, size_t InputWidth, EInitialization Init, size_t input4D, size_t output4D,
              size_t FilterHeight, size_t FilterWidth, size_t FilterDepth, size_t strideX, size_t strideY, size_t strideZ, size_t PaddingHeight,
              size_t PaddingWidth, size_t PaddingDepth, Scalar_t DropoutProbability, EActivationFunction f, ERegularization Reg,
              Scalar_t WeightDecay);

   /*! Copy the conv layer provided as a pointer */
   TConv3DLayer(TConv3DLayer<Architecture_t> *layer);

   /*! Copy constructor. */
   TConv3DLayer(const TConv3DLayer &);

   /*! Destructor. */
   virtual ~TConv3DLayer();

   //virtual void Initialize();

   /*! Computes activation of the layer for the given input. The input
   * must be in 3D tensor form with the different matrices corresponding to
   * different events in the batch. Computes activations as well as
   * the first partial derivative of the activation function at those
   * activations. */
   void Forward(Tensor_t &input, bool applyDropout = false);

   /*! Compute weight, bias and activation gradients. Uses the precomputed
    *  first partial derviatives of the activation function computed during
    *  forward propagation and modifies them. Must only be called directly
    *  at the corresponding call to Forward(...). */
   void Backward(Tensor_t &gradients_backward, const Tensor_t &activations_backward);
   ////              Tensor_t &inp1, Tensor_t &inp2);

   /*! Prints the info about the layer. */
   void Print() const;

   /*! Writes the information and the weights about the layer in an XML node. */
   virtual void AddWeightsXMLTo(void *parent);

   /*! Read the information and the weights about the layer from XML node. */
   virtual void ReadWeightsFromXML(void *parent);

   /*! Getters */
   size_t GetInput4D() const { return finput4D; }
   size_t GetOutput4D() const { return foutput4D; }
   size_t GetFilterDepth() const { return fFilterDepth; }
   size_t GetFilterHeight() const { return fFilterHeight; }
   size_t GetFilterWidth() const { return fFilterWidth; }

   size_t GetstrideX() const { return fStrideX; }
   size_t GetstrideY() const { return fStrideY; }
   size_t GetstrideZ() const { return fStrideZ; }

   size_t GetPaddingHeight() const { return fPaddingHeight; }
   size_t GetPaddingWidth() const { return fPaddingWidth; }
   size_t GetPaddingDepth() const { return fPaddingDepth; }

   size_t GetNLocalViewPixels() const { return fNLocalViewPixels; }
   size_t GetNLocalViews() const { return fNLocalViews; }

   Scalar_t GetDropoutProbability() const { return fDropoutProbability; }

   const Tensor_t &GetInputActivation() const { return fInputActivation; }
   Tensor_t &GetInputActivation() { return fInputActivation; }

   Matrix_t &GetInputActivationAt(size_t i) { return fInputActivation[i]; }
   const Matrix_t &GetInputActivationAt(size_t i) const { return fInputActivation[i]; }

   const Tensor_t &GetForwardMatrices() const { return fForwardTensor; }
   Tensor_t &GetForwardMatrices() { return fForwardTensor; }

   EActivationFunction GetActivationFunction() const { return fF; }
   ERegularization GetRegularization() const { return fReg; }
   Scalar_t GetWeightDecay() const { return fWeightDecay; }

   // The following getters are used for testing
   TDescriptors * GetDescriptors() {return fDescriptors;}
   const TDescriptors * GetDescriptors() const {return fDescriptors;}

   TWorkspace * GetWorkspace() {return fWorkspace;}
   const TWorkspace * GetWorkspace() const {return fWorkspace;}
};


//
//
//  Conv Layer Class - Implementation
//______________________________________________________________________________
template <typename Architecture_t>
TConv3DLayer<Architecture_t>::TConv3DLayer(size_t batchSize, size_t inputDepth, size_t inputHeight, size_t inputWidth,
                                      EInitialization init, size_t input4D, size_t output4D, size_t filterHeight, size_t filterWidth, 
                                       size_t filterDepth, size_t strideX, size_t strideY, size_t strideZ, size_t paddingHeight, 
                                       size_t paddingWidth, size_t paddingDepth, Scalar_t dropoutProbability, EActivationFunction f, 
                                       ERegularization reg, Scalar_t weightDecay)
   : VGeneralLayer<Architecture_t>(batchSize, inputDepth, inputHeight, inputWidth, inputDepth,
                                   calculateDimension(inputHeight, filterHeight, paddingHeight, strideX),
                                   calculateDimension(inputWidth, filterWidth, paddingWidth, strideY),
                                   calculateDimension(inputDepth, filterDepth, paddingDepth, strideZ),
                                   output4D, calculateNLocalViewPixels(filterDepth, filterHeight, filterWidth, input4D),
                                   1, output4D, 1, batchSize, output4D,
                                   calculateNLocalViews(inputHeight, filterHeight, paddingHeight, strideX,
                                                        inputWidth, filterWidth, paddingWidth, strideY,
                                                        inputDepth, filterDepth, paddingDepth, strideZ),
                                   init),
     fFilterDepth(filterDepth), fFilterHeight(filterHeight), fFilterWidth(filterWidth), fStrideX(strideX),
     fStrideY(strideY), fStrideZ(strideZ), fNLocalViewPixels(calculateNLocalViewPixels(filterDepth, filterHeight, filterWidth, input4D)),
     fNLocalViews(calculateNLocalViews(inputHeight, filterHeight, paddingHeight, strideX,
                                       inputWidth, filterWidth, paddingWidth, strideY, 
                                       inputDepth, filterDepth, paddingDepth, strideZ)),
     fDropoutProbability(dropoutProbability), fPaddingHeight(paddingHeight), fPaddingWidth(paddingWidth), fPaddingDepth(paddingDepth),
     fInputActivation(), fF(f), fReg(reg), fWeightDecay(weightDecay)
{
   /** Each element in the vector is a `T_Matrix` representing an event, therefore `vec.size() == batchSize`.
    *  Cells in these matrices are distributed in the following manner:
    *  Each row represents a single feature map, therefore we have `nRows == depth`.
    *  Each column represents a single pixel in that feature map, therefore we have `nCols == nLocalViews`.
    **/

   fInputActivation = Tensor_t( batchSize, output4D, fNLocalViews);     // create tensor (shape is B x C x LV)
   fForwardTensor = Tensor_t ( batchSize, fNLocalViews, fNLocalViewPixels );
   finput4D = input4D;
   foutput4D = output4D;

   InitializeDescriptors();
   InitializeWorkspace();
}

//______________________________________________________________________________
template <typename Architecture_t>
TConv3DLayer<Architecture_t>::TConv3DLayer(TConv3DLayer<Architecture_t> *layer)
   : VGeneralLayer<Architecture_t>(layer), finput4D(layer->GetInput4D()), foutput4D(layer->GetOutput4D()), fFilterDepth(layer->GetFilterDepth()),
     fFilterHeight(layer->GetFilterHeight()), fFilterWidth(layer->GetFilterWidth()),
     fStrideX(layer->GetstrideX()), fStrideY(layer->GetstrideY()), fStrideZ(layer->GetstrideZ()),
     fNLocalViewPixels(layer->GetNLocalViewPixels()), fNLocalViews(layer->GetNLocalViews()),
     fDropoutProbability(layer->GetDropoutProbability()), fPaddingHeight(layer->GetPaddingHeight()),
     fPaddingWidth(layer->GetPaddingWidth()), fPaddingDepth(layer->GetPaddingDepth()),
     fInputActivation( layer->GetInputActivation().GetShape() ),
     fF(layer->GetActivationFunction()),
     fReg(layer->GetRegularization()), fWeightDecay(layer->GetWeightDecay()),
     fForwardTensor( layer->GetForwardMatrices().GetShape() )
{
   InitializeDescriptors();
   InitializeWorkspace();

}

//______________________________________________________________________________
template <typename Architecture_t>
TConv3DLayer<Architecture_t>::TConv3DLayer(const TConv3DLayer &conv3DLayer)
   :  VGeneralLayer<Architecture_t>(conv3DLayer), fFilterDepth(conv3DLayer.fFilterDepth),
      fFilterHeight(conv3DLayer.fFilterHeight), fFilterWidth(conv3DLayer.fFilterWidth), fStrideX(conv3DLayer.fStrideX),
      fStrideY(conv3DLayer.fStrideY), fStrideZ(conv3DLayer.fStrideZ), fNLocalViewPixels(conv3DLayer.fNLocalViewPixels),
      fNLocalViews(conv3DLayer.fNLocalViews), fDropoutProbability(conv3DLayer.fDropoutProbability),
      fPaddingHeight(conv3DLayer.fPaddingHeight), fPaddingWidth(conv3DLayer.fPaddingWidth), fPaddingDepth(conv3DLayer.fPaddingDepth),
      fInputActivation( conv3DLayer.GetInputActivation().GetShape() ),
      fF(conv3DLayer.fF),
      fReg(conv3DLayer.fReg), fWeightDecay(conv3DLayer.fWeightDecay),
      fForwardTensor( conv3DLayer.GetForwardMatrices().GetShape() )
{
   InitializeDescriptors();
   InitializeWorkspace();
}

//______________________________________________________________________________
//FIXME: Add function for cudaFree
template <typename Architecture_t>
TConv3DLayer<Architecture_t>::~TConv3DLayer()
{
   //std::cout << "!!!!Delete conv layer " << this->GetOutput().GetShape()[1] << "  " << this->GetOutput().GetShape()[2] << "  " << this->GetOutput().GetShape()[3] << std::endl;
   if (fDescriptors) {
      ReleaseDescriptors();
      delete fDescriptors;
   }

   if (fWorkspace) {
      FreeWorkspace();
      delete fWorkspace;
   }
}


//______________________________________________________________________________
template <typename Architecture_t>
auto TConv3DLayer<Architecture_t>::Forward(Tensor_t &input, bool /*applyDropout*/) -> void
{

   TConv3DParams params(this->GetBatchSize(), this->GetInputDepth(), this->GetInputHeight(), this->GetInputWidth(), this->GetInput4D(), this->GetOutput4D(),
                       this->GetFilterHeight(), this->GetFilterWidth(), this->GetFilterDepth(),
                      this->GetstrideX(), this->GetstrideY(), this->GetstrideZ(), this->GetPaddingHeight(), this->GetPaddingWidth(), this->GetPaddingDepth());

   //R__ASSERT( input.size() > 0);
   Architecture_t::Conv3DLayerForward(this->GetOutput(), this->GetInputActivation(), input, this->GetWeightsAt(0),
                                    this->GetBiasesAt(0), params, this->GetActivationFunction(),
                                    this->GetForwardMatrices(), (CNN::TCNNDescriptors<TMVA::DNN::CNN::TConvLayer<Architecture_t>> &) (*fDescriptors),
                                    (CNN::TCNNWorkspace<TMVA::DNN::CNN::TConvLayer<Architecture_t>> &) (*fWorkspace));
}

//______________________________________________________________________________
template <typename Architecture_t>
auto TConv3DLayer<Architecture_t>::Backward(Tensor_t &gradients_backward,
                                          const Tensor_t &activations_backward) -> void
//                                          Tensor_t & /*inp1*/, Tensor_t &
//                                          /*inp2*/) -> void
{

	TConv3DParams params(this->GetBatchSize(), this->GetInputDepth(), this->GetInputHeight(), this->GetInputWidth(), this->GetInput4D(), this->GetOutput4D(),
                       this->GetFilterHeight(), this->GetFilterWidth(), this->GetFilterDepth(),
                      this->GetstrideX(), this->GetstrideY(), this->GetstrideZ(), this->GetPaddingHeight(), this->GetPaddingWidth(), this->GetPaddingDepth());


   Architecture_t::Conv3DLayerBackward(
      gradients_backward, this->GetWeightGradientsAt(0), this->GetBiasGradientsAt(0), this->GetInputActivation(),
      this->GetActivationGradients(), this->GetWeightsAt(0), activations_backward, this->GetOutput(),
      this->GetActivationFunction(),
      (CNN::TCNNDescriptors<TMVA::DNN::CNN::TConvLayer<Architecture_t>> &) (*fDescriptors),
      (CNN::TCNNWorkspace<TMVA::DNN::CNN::TConvLayer<Architecture_t>> &) (*fWorkspace),
      params);

   addRegularizationGradients<Architecture_t>(this->GetWeightGradientsAt(0), this->GetWeightsAt(0),
                                              this->GetWeightDecay(), this->GetRegularization());
}

//______________________________________________________________________________
template <typename Architecture_t>
auto TConv3DLayer<Architecture_t>::Print() const -> void
{
   std::cout << " CONV LAYER: \t";
   std::cout << "( W = " << this->GetWidth() << " , ";
   std::cout << " H = " << this->GetHeight() << " , ";
   std::cout << " D = " << this->GetDepth() << " ) ";

   std::cout << "\t Filter ( W = " << this->GetFilterWidth() << " , ";
   std::cout << " H = " << this->GetFilterHeight() << " , ";
   std::cout << " D = " << this->GetFilterDepth() << " ) ";
   //std::cout << "\t Local Views = " << this->GetNLocalViews()  << " " ;
   if (this->GetOutput().GetSize() > 0) {
      std::cout << "\tOutput = ( " << this->GetOutput().GetFirstSize() << " , "
                << this->GetOutput().GetCSize() << " , " << this->GetOutput().GetHSize() << " , " << this->GetOutput().GetWSize()
                << " ) ";
   }
   std::vector<std::string> activationNames = { "Identity","Relu","Sigmoid","Tanh","SymmRelu","SoftSign","Gauss" };
   std::cout << "\t Activation Function = ";
   std::cout << activationNames[ static_cast<int>(fF) ] << std::endl;
}

//______________________________________________________________________________
template <typename Architecture_t>
void TConv3DLayer<Architecture_t>::AddWeightsXMLTo(void *parent)
{
   auto layerxml = gTools().xmlengine().NewChild(parent, 0, "Conv3DLayer");

   gTools().xmlengine().NewAttr(layerxml, 0, "Depth", gTools().StringFromInt(this->GetDepth()));
   gTools().xmlengine().NewAttr(layerxml, 0, "FilterHeight", gTools().StringFromInt(this->GetFilterHeight()));
   gTools().xmlengine().NewAttr(layerxml, 0, "FilterWidth", gTools().StringFromInt(this->GetFilterWidth()));
   gTools().xmlengine().NewAttr(layerxml, 0, "FilterDepth", gTools().StringFromInt(this->GetFilterDepth()));
   gTools().xmlengine().NewAttr(layerxml, 0, "strideX", gTools().StringFromInt(this->GetstrideX()));
   gTools().xmlengine().NewAttr(layerxml, 0, "strideY", gTools().StringFromInt(this->GetstrideY()));
   gTools().xmlengine().NewAttr(layerxml, 0, "strideZ", gTools().StringFromInt(this->GetstrideZ()));
   gTools().xmlengine().NewAttr(layerxml, 0, "PaddingHeight", gTools().StringFromInt(this->GetPaddingHeight()));
   gTools().xmlengine().NewAttr(layerxml, 0, "PaddingWidth", gTools().StringFromInt(this->GetPaddingWidth()));
   gTools().xmlengine().NewAttr(layerxml, 0, "PaddingDepth", gTools().StringFromInt(this->GetPaddingDepth()));


   int activationFunction = static_cast<int>(this -> GetActivationFunction());
   gTools().xmlengine().NewAttr(layerxml, 0, "ActivationFunction",
                                TString::Itoa(activationFunction, 10));

   // write weights and bias matrix
   this->WriteMatrixToXML(layerxml, "Weights", this -> GetWeightsAt(0));
   this->WriteMatrixToXML(layerxml, "Biases",  this -> GetBiasesAt(0));
}

//______________________________________________________________________________
template <typename Architecture_t>
void TConv3DLayer<Architecture_t>::ReadWeightsFromXML(void *parent)
{
   // read weights and biases
   // the meta information is read before because it is needed before creating the Conv layer
   this->ReadMatrixXML(parent,"Weights", this -> GetWeightsAt(0));
   this->ReadMatrixXML(parent,"Biases", this -> GetBiasesAt(0));
}

template <typename Architecture_t>
size_t TConv3DLayer<Architecture_t>::calculateDimension(size_t imgDim, size_t fltDim, size_t padding, size_t stride)
{
   size_t temp = imgDim - fltDim + 2 * padding;
   if (temp % stride || temp + stride <= 0) {
      Fatal("calculateDimension", "Not compatible hyper parameters for layer - (imageDim, filterDim, padding, stride) "
            "%zu, %zu, %zu, %zu", imgDim, fltDim, padding, stride);
   }
   return temp / stride + 1;
}

template <typename Architecture_t>

size_t TConv3DLayer<Architecture_t>::calculateNLocalViews(size_t inputHeight, size_t filterHeight, size_t paddingHeight, size_t strideX,
                               size_t inputWidth, size_t filterWidth, size_t paddingWidth, size_t strideY,
                               size_t inputDepth, size_t filterDepth, size_t paddingDepth, size_t strideZ)
{
    int height = calculateDimension(inputHeight, filterHeight, paddingHeight, strideX);
    int width = calculateDimension(inputWidth, filterWidth, paddingWidth, strideY);
    int depth = calculateDimension(inputDepth, filterDepth, paddingDepth, strideZ);

    return height * width * depth;
}

//______________________________________________________________________________
template <typename Architecture_t>
void TConv3DLayer<Architecture_t>::InitializeDescriptors() {
  // Architecture_t::InitializeConvDescriptors(fDescriptors, this);
}

template <typename Architecture_t>
void TConv3DLayer<Architecture_t>::ReleaseDescriptors() {
   //Architecture_t::ReleaseConvDescriptors(fDescriptors);
}

//______________________________________________________________________________
template <typename Architecture_t>
void TConv3DLayer<Architecture_t>::InitializeWorkspace() {


   //Architecture_t::InitializeConvWorkspace(fWorkspace, fDescriptors, params, this);
}

template <typename Architecture_t>
void TConv3DLayer<Architecture_t>::FreeWorkspace() {
   //Architecture_t::FreeConvWorkspace(fWorkspace);
}

//______________________________________________________________________________

} // namespace CNN
} // namespace DNN
} // namespace TMVA

#endif