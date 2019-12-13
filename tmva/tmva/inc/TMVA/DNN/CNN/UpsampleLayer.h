// @(#)root/tmva/tmva/dnn:$Id$
// Author: Ashish Kshirsagar

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : TUpsampleLayer                                                         *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Upsample Deep Neural Network Layer                                        *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Ashish Kshirsagar      <ashishkshirsagar10@gmail.com>                     *
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

#ifndef UPSAMPLELAYER_H_
#define UPSAMPLELAYER_H_

#include "TMatrix.h"

#include "TMVA/DNN/CNN/TransConvLayer.h"
#include "TMVA/DNN/Functions.h"

#include <iostream>

namespace TMVA {
namespace DNN {
namespace CNN {

/** \class TUpsampleLayer

    Generic Upsample Layer class.

    This generic Upsample Layer Class represents a upsampling layer of
    a CNN. It inherits all of the properties of the transpose convolutional layer
    TTransConvLayer, but it overrides the propagation methods. In a sense, upsample
    can be seen as a linear convolution.
    
    The height and width of the weights and biases is set to 0, since this
    layer does not contain any weights.

 */
template <typename Architecture_t>
class TUpsampleLayer : public TTransConvLayer<Architecture_t> {

public:
    using Matrix_t = typename Architecture_t::Matrix_t;
    using Scalar_t = typename Architecture_t::Scalar_t;

private:
   std::vector<Matrix_t> indexMatrix; ///< Matrix of indices for the backward pass.

public:
   /*! Constructor. */
   TUpsampleLayer(size_t BatchSize, size_t InputDepth, size_t InputHeight, size_t InputWidth, size_t FilterHeight,
                 size_t FilterWidth, size_t StrideRows, size_t StrideCols, Scalar_t DropoutProbability);

   /*! Copy the upsample layer provided as a pointer */
   TUpsampleLayer(TUpsampleLayer<Architecture_t> *layer);

   /*! Copy constructor. */
   TUpsampleLayer(const TUpsampleLayer &);

   /*! Destructor. */
   ~TUpsampleLayer();

   /*! Computes activation of the layer for the given input. The input
    *  must be in 3D tensor form with the different matrices corresponding to
    *  different events in the batch. It spatially upsamples the input
    *  matrices. */
   void Forward(std::vector<Matrix_t> &input, bool applyDropout = false);

   /*! Depending on the winning units determined during the Forward pass,
    *  it only forwards the derivatives to the right units in the previous
    *  layer. Must only be called directly at the corresponding call
    *  to Forward(...). */
   void Backward(std::vector<Matrix_t> &gradients_backward, const std::vector<Matrix_t> &activations_backward,
                 std::vector<Matrix_t> &inp1, std::vector<Matrix_t> &inp2);

   /*! Writes the information and the weights about the layer in an XML node. */
   virtual void AddWeightsXMLTo(void *parent);

   /*! Read the information and the weights about the layer from XML node. */
   virtual void ReadWeightsFromXML(void *parent);

   /*! Prints the info about the layer. */
   void Print() const;

   /*! Getters */
   const std::vector<Matrix_t> &GetIndexMatrix() const { return indexMatrix; }
   std::vector<Matrix_t> &GetIndexMatrix() { return indexMatrix; }

};

//______________________________________________________________________________
template <typename Architecture_t>
TUpsampleLayer<Architecture_t>::TUpsampleLayer(size_t batchSize, size_t inputDepth, size_t inputHeight, size_t inputWidth,
                                             size_t filterHeight, size_t filterWidth, size_t strideRows,
                                             size_t strideCols, Scalar_t dropoutProbability)

        : TTransConvLayer<Architecture_t>(batchSize, inputDepth, inputHeight, inputWidth, inputDepth, EInitialization::kZero,
                                     filterHeight, filterWidth, strideRows, strideCols, 0, 0, dropoutProbability,
                                     EActivationFunction::kIdentity, ERegularization::kNone, 0),
          indexMatrix()
{
   for (size_t i = 0; i < this->GetBatchSize(); i++) {
      indexMatrix.emplace_back(this->GetDepth(), this->GetNLocalViews());
   }
}

//______________________________________________________________________________
template <typename Architecture_t>
TUpsampleLayer<Architecture_t>::TUpsampleLayer(TUpsampleLayer<Architecture_t> *layer)
   : TTransConvLayer<Architecture_t>(layer), indexMatrix()
{
   for (size_t i = 0; i < layer->GetBatchSize(); i++) {
      indexMatrix.emplace_back(layer->GetDepth(), layer->GetNLocalViews());
   }
}

//______________________________________________________________________________
template <typename Architecture_t>
TUpsampleLayer<Architecture_t>::TUpsampleLayer(const TUpsampleLayer &layer)
   : TTransConvLayer<Architecture_t>(layer), indexMatrix()
{
   for (size_t i = 0; i < layer.fBatchSize; i++) {
      indexMatrix.emplace_back(layer.fDepth, layer.fNLocalViews);
   }
}

//______________________________________________________________________________
template <typename Architecture_t>
TUpsampleLayer<Architecture_t>::~TUpsampleLayer()
{
}

//______________________________________________________________________________
template <typename Architecture_t>
auto TUpsampleLayer<Architecture_t>::Forward(std::vector<Matrix_t> &input, bool applyDropout) -> void
{
   for (size_t i = 0; i < this->GetBatchSize(); i++) {

      if (applyDropout && (this->GetDropoutProbability() != 1.0)) {
         Architecture_t::Dropout(input[i], this->GetDropoutProbability());
      }

      Architecture_t::Upsample(this->GetOutputAt(i), indexMatrix[i], input[i], this->GetInputHeight(),
                                 this->GetInputWidth(), this->GetFilterHeight(), this->GetFilterWidth(),
                                 this->GetStrideRows(), this->GetStrideCols());
   }
}

//______________________________________________________________________________
template <typename Architecture_t>
auto TUpsampleLayer<Architecture_t>::Backward(std::vector<Matrix_t> &gradients_backward,
                                             const std::vector<Matrix_t> & /*activations_backward*/,
                                             std::vector<Matrix_t> & /*inp1*/, std::vector<Matrix_t> &
                                             /*inp2*/) -> void
{
   for (size_t i = 0; i < this->GetBatchSize(); i++) {
      Architecture_t::UpsampleLayerBackward(gradients_backward[i], this->GetActivationGradients()[i],
                                           this->GetIndexMatrix()[i],
                                           this->GetInputHeight(), this->GetInputWidth(),
                                           this->GetFilterHeight(), this->GetFilterWidth(),
                                           this->GetStrideRows(), this->GetStrideCols(), this->GetNLocalViews());
   }
}

//______________________________________________________________________________
template <typename Architecture_t>
auto TUpsampleLayer<Architecture_t>::Print() const -> void
{
   std::cout << " UPSAMPLE Layer: \t";
   std::cout << "( W = " << this->GetWidth() << " , ";
   std::cout << " H = " << this->GetHeight() << " , ";
   std::cout << " D = " << this->GetDepth() << " ) ";

   std::cout << "\t Filter ( W = " << this->GetFilterWidth() << " , ";
   std::cout << " H = " << this->GetFilterHeight() << " ) ";

   if (this->GetOutput().size() > 0) {
      std::cout << "\tOutput = ( " << this->GetOutput().size() << " , " << this->GetOutput()[0].GetNrows() << " , " << this->GetOutput()[0].GetNcols() << " ) ";
   }
   std::cout << std::endl;
}

//______________________________________________________________________________
template <typename Architecture_t>
void TUpsampleLayer<Architecture_t>::AddWeightsXMLTo(void *parent)
{
   auto layerxml = gTools().xmlengine().NewChild(parent, 0, "UpsampleLayer");

   // write  maxpool layer info
   gTools().xmlengine().NewAttr(layerxml, 0, "FilterHeight", gTools().StringFromInt(this->GetFilterHeight()));
   gTools().xmlengine().NewAttr(layerxml, 0, "FilterWidth", gTools().StringFromInt(this->GetFilterWidth()));
   gTools().xmlengine().NewAttr(layerxml, 0, "StrideRows", gTools().StringFromInt(this->GetStrideRows()));
   gTools().xmlengine().NewAttr(layerxml, 0, "StrideCols", gTools().StringFromInt(this->GetStrideCols()));

}

//______________________________________________________________________________
template <typename Architecture_t>
void TUpsampleLayer<Architecture_t>::ReadWeightsFromXML(void * /*parent */)
{
   // all info is read before - nothing to do 
}

} // namespace CNN
} // namespace DNN
} // namespace TMVA

#endif
