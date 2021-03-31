// @(#)root/tmva/tmva/dnn:$Id$
// Author: Surya S Dwivedi

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : TMaxPool3DLayer                                                         *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      3D MaxPool Deep Neural Network Layer                                      *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Surya S Dwivedi  <surya1997@utexas.edu> - Univ of Texas Austin            *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/


#ifndef MAXPOOL3DLAYER_H_
#define MAXPOOL3DLAYER_H_

#include "TMatrix.h"

#include "TMVA/DNN/CNN_3D/Conv3DLayer.h"
#include "TMVA/DNN/CNN/MaxPoolLayer.h"

#include "TMVA/DNN/Functions.h"
#include "TMVA/DNN/CNN_3D/ContextHandles.h"

#include <iostream>

namespace TMVA {
namespace DNN {

namespace CNN_3D {



/** \class TMaxPool3DLayer

    Generic Max Pooling Layer class.

    This generic Max Pooling Layer Class represents a pooling layer of
    a CNN. It inherits all of the properties of the convolutional layer
    TConvLayer, but it overrides the propagation methods. In a sense, max pooling
    can be seen as non-linear convolution: a filter slides over the input and produces
    one element as a function of the the elements within the receptive field.
    In addition to that, it contains a matrix of winning units.

    The height and width of the weights and biases is set to 0, since this
    layer does not contain any weights.

 */
template <typename Architecture_t>
class TMaxPool3DLayer : public VGeneralLayer<Architecture_t>  {

public:
   using Tensor_t = typename Architecture_t::Tensor_t;
   using Matrix_t = typename Architecture_t::Matrix_t;
   using Scalar_t = typename Architecture_t::Scalar_t;

   using LayerDescriptor_t = typename Architecture_t::PoolingDescriptor_t;
   using WeightsDescriptor_t = typename Architecture_t::EmptyDescriptor_t ;
   using HelperDescriptor_t = typename Architecture_t::DropoutDescriptor_t;

   // do we need thse ???
   using AlgorithmForward_t = typename Architecture_t::AlgorithmForward_t;   // Forward layer operation
   using AlgorithmBackward_t = typename Architecture_t::AlgorithmBackward_t; // Backward layer operation
   using AlgorithmHelper_t = typename Architecture_t::AlgorithmHelper_t;     // Used for weight grad backward pass

   // FIXME: Add other cudnn types (algorithm preference etc.)
   using AlgorithmDataType_t = typename Architecture_t::AlgorithmDataType_t;

   using ReduceTensorDescriptor_t = typename Architecture_t::ReduceTensorDescriptor_t; // used for reduction of tensor(bias grad)

protected:
   size_t fFilterDepth;  ///< The depth of the filter.
   size_t fFilterHeight; ///< The height of the filter.
   size_t fFilterWidth;  ///< The width of the filter.

   size_t fStrideX; ///< The number of row pixels to slid the filter each step.
   size_t fStrideY; ///< The number of column pixels to slid the filter each step.
   size_t fstrideZ; ///< The number of depth pixels to slid the filter each step.

   size_t fNLocalViewPixels; ///< The number of pixels in one local image view.
   size_t fNLocalViews;      ///< The number of local views in one image.

   Scalar_t fDropoutProbability; ///< Probability that an input is active.

   TDescriptors *fDescriptors = nullptr; ///< Keeps the convolution, activations and filter descriptors

   TWorkspace *fWorkspace = nullptr;

private:
   Tensor_t fIndexTensor; ///< Matrix of indices for the backward pass.

   void InitializeDescriptors();
   void ReleaseDescriptors();
   void InitializeWorkspace();
   void FreeWorkspace();

public:
   /*! Constructor. */
   TMaxPool3DLayer(size_t BatchSize, size_t InputDepth, size_t InputHeight, size_t InputWidth, size_t FilterHeight,
                 size_t FilterWidth, size_t FilterDepth, size_t strideX, size_t strideY, size_t strideZ, Scalar_t DropoutProbability);

   /*! Copy the max pooling layer provided as a pointer */
   TMaxPool3DLayer(TMaxPool3DLayer<Architecture_t> *layer);

   /*! Copy constructor. */
   TMaxPool3DLayer(const TMaxPool3DLayer &);

   /*! Destructor. */
   virtual ~TMaxPool3DLayer();

   // virtual void Initialize();

   /*! Computes activation of the layer for the given input. The input
    *  must be in 3D tensor form with the different matrices corresponding to
    *  different events in the batch. It spatially downsamples the input
    *  matrices. */
   void Forward(Tensor_t &input, bool applyDropout = true);

   /*! Depending on the winning units determined during the Forward pass,
    *  it only forwards the derivatives to the right units in the previous
    *  layer. Must only be called directly at the corresponding call
    *  to Forward(...). */
   void Backward(Tensor_t &gradients_backward, const Tensor_t &activations_backward);
   //             Tensor_t &inp1, Tensor_t &inp2);

   /*! Writes the information and the weights about the layer in an XML node. */
   virtual void AddWeightsXMLTo(void *parent);

   /*! Read the information and the weights about the layer from XML node. */
   virtual void ReadWeightsFromXML(void *parent);

   /*! Prints the info about the layer. */
   void Print() const;

   /*! Getters */
   size_t GetFilterDepth() const { return fFilterDepth; }
   size_t GetFilterHeight() const { return fFilterHeight; }
   size_t GetFilterWidth() const { return fFilterWidth; }

   size_t GetStrideX() const { return fStrideX; }
   size_t GetStrideY() const { return fStrideY; }
   size_t GetStrideZ() const { return fstrideZ; }

   size_t GetNLocalViews() const { return fNLocalViews; }

   Scalar_t GetDropoutProbability() const { return fDropoutProbability; }

   const Tensor_t & GetIndexTensor() const { return fIndexTensor; }
   Tensor_t & GetIndexTensor() { return fIndexTensor; }

   // The following getters are used for testing
   TDescriptors *GetDescriptors() { return fDescriptors; }
   const TDescriptors *GetDescriptors() const { return fDescriptors; }

   TWorkspace *GetWorkspace() { return fWorkspace; }
   const TWorkspace *GetWorkspace() const { return fWorkspace; }
};

//______________________________________________________________________________
template <typename Architecture_t>
TMaxPool3DLayer<Architecture_t>::TMaxPool3DLayer(size_t batchSize, size_t inputDepth, size_t inputHeight, size_t inputWidth,
                                             size_t filterHeight, size_t filterWidth, size_t filterDepth, size_t strideX,
                                             size_t strideY, size_t strideZ, Scalar_t dropoutProbability)
   : VGeneralLayer<Architecture_t>(
        batchSize, inputDepth, inputHeight, inputWidth, inputDepth,
        TConv3DLayer<Architecture_t>::calculateDimension(inputHeight, filterHeight, 0, strideX),
        TConv3DLayer<Architecture_t>::calculateDimension(inputWidth, filterWidth, 0, strideY),
        TConv3DLayer<Architecture_t>::calculateDimension(inputDepth, filterDepth, 0, strideZ), 0, 0, 0, 0,
        0, // weights dimensions
        batchSize, inputDepth,
        TConv3DLayer<Architecture_t>::calculateNLocalViews(inputHeight, filterHeight, 0, strideX, inputWidth,
                                                         filterWidth, 0, strideY, inputDepth, filterDepth, 0, strideZ),
        EInitialization::kZero),
     fFilterDepth(inputDepth), fFilterHeight(filterHeight), fFilterWidth(filterWidth), fStrideX(strideX),
     fStrideY(strideY), fstrideZ(strideZ),
     fNLocalViews(TConv3DLayer<Architecture_t>::calculateNLocalViews(inputHeight, filterHeight, 0, strideX,
                                                                   inputWidth, filterWidth, 0, strideY,
                                                                   inputDepth, filterDepth, 0, strideZ)),
     fDropoutProbability(dropoutProbability), fIndexTensor(batchSize, inputDepth, fNLocalViews)
{
   InitializeDescriptors();
   InitializeWorkspace();
}

//______________________________________________________________________________
template <typename Architecture_t>
TMaxPool3DLayer<Architecture_t>::TMaxPool3DLayer(TMaxPool3DLayer<Architecture_t> *layer)
   : VGeneralLayer<Architecture_t>(layer), fFilterDepth(layer->GetFilterDepth()),
     fFilterHeight(layer->GetFilterHeight()), fFilterWidth(layer->GetFilterWidth()),
     fStrideX(layer->GetStrideX()), fStrideY(layer->GetStrideY()), fstrideZ(layer->GetStrideZ()), fNLocalViews(layer->GetNLocalViews()),
     fDropoutProbability(layer->GetDropoutProbability()), fIndexTensor(layer->GetIndexTensor().GetShape())
{
   InitializeDescriptors();
   InitializeWorkspace();
}

//______________________________________________________________________________
template <typename Architecture_t>
TMaxPool3DLayer<Architecture_t>::TMaxPool3DLayer(const TMaxPool3DLayer &layer)
   : VGeneralLayer<Architecture_t>(layer), fFilterDepth(layer.fFilterDepth), fFilterHeight(layer.fFilterHeight),
     fFilterWidth(layer.fFilterWidth), fStrideX(layer.fStrideX), fStrideY(layer.fStrideY), fstrideZ(layer.fstrideZ),
     fNLocalViews(layer.fNLocalViews), fDropoutProbability(layer.fDropoutProbability),
     fIndexTensor(layer.GetIndexTensor().GetShape())
{
   InitializeDescriptors();
   InitializeWorkspace();
}

//______________________________________________________________________________
template <typename Architecture_t>
TMaxPool3DLayer<Architecture_t>::~TMaxPool3DLayer()
{
   if (fDescriptors) {
      ReleaseDescriptors();
      delete fDescriptors;
      fDescriptors = nullptr;
   }

   if (fWorkspace) {
      FreeWorkspace();
      delete fWorkspace;
      fWorkspace = nullptr;
   }
}

//______________________________________________________________________________
template <typename Architecture_t>
auto TMaxPool3DLayer<Architecture_t>::Forward(Tensor_t &input, bool applyDropout) -> void
{
   if (applyDropout && (this->GetDropoutProbability() != 1.0)) {
      Architecture_t::DropoutForward(input, fDescriptors, fWorkspace, this->GetDropoutProbability());
   }

   Architecture_t::Downsample3D(
      this->GetOutput(), fIndexTensor, input, (CNN::TCNNDescriptors<TMVA::DNN::CNN::TMaxPoolLayer<Architecture_t>> &) *fDescriptors,
      (CNN::TCNNWorkspace<TMVA::DNN::CNN::TMaxPoolLayer<Architecture_t>> &) *fWorkspace, this->GetInputHeight(), this->GetInputWidth(), this->GetInputDepth(),
      this->GetFilterHeight(), this->GetFilterWidth(), this->GetFilterDepth(), this->GetStrideX(), this->GetStrideY(), this->GetStrideZ());
}

//______________________________________________________________________________
template <typename Architecture_t>
auto TMaxPool3DLayer<Architecture_t>::Backward(Tensor_t &gradients_backward, const Tensor_t &activations_backward) -> void
//                                             Tensor_t & /*inp1*/, Tensor_t &
{

   if (this->GetDropoutProbability() != 1.0) {
      Architecture_t::DropoutBackward(this->GetActivationGradients(), fDescriptors, fWorkspace);
   }
   Architecture_t::MaxPoolLayer3DBackward(
      gradients_backward, this->GetActivationGradients(), fIndexTensor, activations_backward, this->GetOutput(),
      (CNN::TCNNDescriptors<TMVA::DNN::CNN::TMaxPoolLayer<Architecture_t>> &) (*fDescriptors),
      (CNN::TCNNWorkspace<TMVA::DNN::CNN::TMaxPoolLayer<Architecture_t>> &) (*fWorkspace), this->GetInputHeight(), this->GetInputWidth(),
      this->GetFilterHeight(), this->GetFilterWidth(),  this->GetStrideX(), this->GetStrideY(),
      this->GetNLocalViews());
}

//______________________________________________________________________________
template <typename Architecture_t>
auto TMaxPool3DLayer<Architecture_t>::Print() const -> void
{
   std::cout << " POOL Layer: \t";
   std::cout << "( W = " << this->GetWidth() << " , ";
   std::cout << " H = " << this->GetHeight() << " , ";
   std::cout << " D = " << this->GetDepth() << " ) ";

   std::cout << "\t Filter ( W = " << this->GetFilterWidth() << " , ";
   std::cout << " H = " << this->GetFilterHeight() << " , ";
   std::cout << " D = " << this->GetFilterDepth() << " ) ";

   if (this->GetOutput().GetSize() > 0) {
      std::cout << "\tOutput = ( " << this->GetOutput().GetFirstSize() << " , " << this->GetOutput().GetCSize()
                << " , " << this->GetOutput().GetHSize() << " , " << this->GetOutput().GetWSize() << " ) ";
   }
   std::cout << std::endl;
}

//______________________________________________________________________________
template <typename Architecture_t>
void TMaxPool3DLayer<Architecture_t>::AddWeightsXMLTo(void * /*parent*/ )
{
   // auto layerxml = gTools().xmlengine().NewChild(parent, 0, "MaxPoolLayer");
   // TODO
   // write  maxpool layer info
   // gTools().xmlengine().NewAttr(layerxml, 0, "FilterHeight", gTools().StringFromInt(this->GetFilterHeight()));
   // gTools().xmlengine().NewAttr(layerxml, 0, "FilterWidth", gTools().StringFromInt(this->GetFilterWidth()));
   // gTools().xmlengine().NewAttr(layerxml, 0, "FilterDepth", gTools().StringFromInt(this->GetFilterDepth()));
   // gTools().xmlengine().NewAttr(layerxml, 0, "strideX", gTools().StringFromInt(this->GetStrideX()));
   // gTools().xmlengine().NewAttr(layerxml, 0, "strideY", gTools().StringFromInt(this->GetStrideY()));
   // gTools().xmlengine().NewAttr(layerxml, 0, "strideZ", gTools().StringFromInt(this->GetStrideZ()));


}

//______________________________________________________________________________
template <typename Architecture_t>
void TMaxPool3DLayer<Architecture_t>::ReadWeightsFromXML(void * /*parent */)
{
   // all info is read before - nothing to do
}

//______________________________________________________________________________
template <typename Architecture_t>
void TMaxPool3DLayer<Architecture_t>::InitializeDescriptors() {
   //Architecture_t::InitializePoolDescriptors(fDescriptors, this);
}

template <typename Architecture_t>
void TMaxPool3DLayer<Architecture_t>::ReleaseDescriptors() {
   //Architecture_t::ReleasePoolDescriptors(fDescriptors);
}

//______________________________________________________________________________
template <typename Architecture_t>
void TMaxPool3DLayer<Architecture_t>::InitializeWorkspace() {
// #if 0  // this is for dropout and needs to be fixed
//    TConvParams params(this->GetBatchSize(), this->GetInputDepth(), this->GetInputHeight(), this->GetInputWidth(),
//                       this->GetDepth(), this->GetFilterHeight(), this->GetFilterWidth(),
//                       this->GetStrideRows(), this->GetStrideCols(), this->GetPaddingHeight(), this->GetPaddingWidth());
//
//    Architecture_t::InitializePoolDropoutWorkspace(fWorkspace, fDescriptors, params, this);
// #endif
}

template <typename Architecture_t>
void TMaxPool3DLayer<Architecture_t>::FreeWorkspace() {
  // Architecture_t::FreePoolDropoutWorkspace(fWorkspace, this);
}

} // namespace CNN
} // namespace DNN
} // namespace TMVA

#endif
