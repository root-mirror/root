// @(#)root/tmva/tmva/dnn:$Id$
// Author: Vladimir Ilievski

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : TMaxPoolLayer                                                         *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Max Pool Deep Neural Network Layer                                        *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Vladimir Ilievski      <ilievski.vladimir@live.com>  - CERN, Switzerland  *
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

#ifndef MAXPOOLLAYER_H_
#define MAXPOOLLAYER_H_

#include "TMatrix.h"

#include "TMVA/DNN/GeneralLayer.h"
#include "TMVA/DNN/Functions.h"

#include <iostream>

namespace TMVA {
namespace DNN {
namespace CNN {

/** \class TMaxPoolLayer

    Generic Max Pooling Layer class.

    This generic Max Pooling Layer Class represents a pooling layer of
    a CNN. It inherits all of the properties of the generic virtual base class
    VGeneralLayer. In addition to that, it contains a matrix of winning units.

    The height and width of the weights and biases is set to 0, since this
    layer does not contain any weights.

 */
template <typename Architecture_t>
class TMaxPoolLayer : public VGeneralLayer<Architecture_t> {
public:
   using Matrix_t = typename Architecture_t::Matrix_t;
   using Scalar_t = typename Architecture_t::Scalar_t;

private:
    bool inline isInteger(Scalar_t x) const { return x == floor(x); }

    /* Calculate the output dimension of the convolutional layer */
    size_t calculateDimension(int imgDim, int fltDim, int stride);

    /* Calculate the number of pixels in a single receptive field */
    size_t inline calculateNLocalViewPixels(int depth, int height, int width) { return depth * height * width; }

    /* Calculate the number of receptive fields in an image given the filter and image sizes */
    size_t calculateNLocalViews(int inputHeight, int filterHeight, int strideRows,
                                int inputWidth, int filterWidth, int strideCols);

private:
   std::vector<Matrix_t> indexMatrix; ///< Matrix of indices for the backward pass.

   size_t fFrameHeight; ///< The height of the frame.
   size_t fFrameWidth;  ///< The width of the frame.

   size_t fStrideRows; ///< The number of row pixels to slid the filter each step.
   size_t fStrideCols; ///< The number of column pixels to slid the filter each step.

   size_t fNLocalViewPixels; ///< The number of pixels in one local image view.
   size_t fNLocalViews;      ///< The number of local views in one image.

   Scalar_t fDropoutProbability; ///< Probability that an input is active.

public:
   /*! Constructor. */
   TMaxPoolLayer(size_t BatchSize, size_t InputDepth, size_t InputHeight, size_t InputWidth, size_t FrameHeight,
                 size_t FrameWidth, size_t StrideRows, size_t StrideCols, Scalar_t DropoutProbability);

   /*! Copy the max pooling layer provided as a pointer */
   TMaxPoolLayer(TMaxPoolLayer<Architecture_t> *layer);

   /*! Copy constructor. */
   TMaxPoolLayer(const TMaxPoolLayer &);

   /*! Destructor. */
   ~TMaxPoolLayer();

   /*! Computes activation of the layer for the given input. The input
    *  must be in 3D tensor form with the different matrices corresponding to
    *  different events in the batch. It spatially downsamples the input
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

   size_t GetFrameHeight() const { return fFrameHeight; }
   size_t GetFrameWidth() const { return fFrameWidth; }

   size_t GetStrideRows() const { return fStrideRows; }
   size_t GetStrideCols() const { return fStrideCols; }

   size_t GetNLocalViewPixels() const { return fNLocalViewPixels; }
   size_t GetNLocalViews() const { return fNLocalViews; }

   Scalar_t GetDropoutProbability() const { return fDropoutProbability; }
};

//______________________________________________________________________________
template <typename Architecture_t>
TMaxPoolLayer<Architecture_t>::TMaxPoolLayer(size_t batchSize, size_t inputDepth, size_t inputHeight, size_t inputWidth,
                                             size_t frameHeight, size_t frameWidth, size_t strideRows,
                                             size_t strideCols, Scalar_t dropoutProbability)
   : VGeneralLayer<Architecture_t>(batchSize, inputDepth, inputHeight, inputWidth, inputDepth,
                                   calculateDimension(inputHeight, frameHeight, strideRows),
                                   calculateDimension(inputWidth, frameWidth, strideCols),
                                   0, 0, 0, 0, 0, 0, batchSize, inputDepth /* I suspect this should be 1 instead */,
                                   calculateNLocalViews(inputHeight, frameHeight, strideRows,
                                                        inputWidth, frameWidth, strideCols),
                                   EInitialization::kZero),
     indexMatrix(), fFrameHeight(frameHeight), fFrameWidth(frameWidth), fStrideRows(strideRows),
     fStrideCols(strideCols), fNLocalViewPixels(calculateNLocalViewPixels(inputDepth, frameHeight, frameWidth)),
     fNLocalViews(calculateNLocalViews(inputHeight, frameHeight, strideRows, inputWidth, frameWidth, strideCols)),
     fDropoutProbability(dropoutProbability)
{
   for (size_t i = 0; i < this->GetBatchSize(); i++) {
      indexMatrix.emplace_back(this->GetDepth(), this->GetNLocalViews());
   }
}

//______________________________________________________________________________
template <typename Architecture_t>
TMaxPoolLayer<Architecture_t>::TMaxPoolLayer(TMaxPoolLayer<Architecture_t> *layer)
   : VGeneralLayer<Architecture_t>(layer), indexMatrix(), fFrameHeight(layer->GetFrameHeight()),
     fFrameWidth(layer->GetFrameWidth()), fStrideRows(layer->GetStrideRows()), fStrideCols(layer->GetStrideCols()),
     fNLocalViewPixels(layer->GetNLocalViewPixels()), fNLocalViews(layer->GetNLocalViews()),
     fDropoutProbability(layer->GetDropoutProbability())
{
   for (size_t i = 0; i < layer->GetBatchSize(); i++) {
      indexMatrix.emplace_back(layer->GetDepth(), layer->GetNLocalViews());
   }
}

//______________________________________________________________________________
template <typename Architecture_t>
TMaxPoolLayer<Architecture_t>::TMaxPoolLayer(const TMaxPoolLayer &layer)
   : VGeneralLayer<Architecture_t>(layer), indexMatrix(), fFrameHeight(layer.fFrameHeight),
     fFrameWidth(layer.fFrameWidth), fStrideRows(layer.fStrideRows), fStrideCols(layer.fStrideCols),
     fNLocalViewPixels(layer.fNLocalViewPixels), fNLocalViews(layer.fNLocalViews),
     fDropoutProbability(layer.fDropoutProbability)
{
   for (size_t i = 0; i < layer.fBatchSize; i++) {
      indexMatrix.emplace_back(layer.fDepth, layer.fNLocalViews);
   }
}

//______________________________________________________________________________
template <typename Architecture_t>
TMaxPoolLayer<Architecture_t>::~TMaxPoolLayer()
{
}

//______________________________________________________________________________
template <typename Architecture_t>
auto TMaxPoolLayer<Architecture_t>::Forward(std::vector<Matrix_t> &input, bool applyDropout) -> void
{
   for (size_t i = 0; i < this->GetBatchSize(); i++) {

      if (applyDropout && (this->GetDropoutProbability() != 1.0)) {
         Architecture_t::Dropout(input[i], this->GetDropoutProbability());
      }

      Architecture_t::Downsample(this->GetOutputAt(i), indexMatrix[i], input[i], this->GetInputHeight(),
                                 this->GetInputWidth(), this->GetFrameHeight(), this->GetFrameWidth(),
                                 this->GetStrideRows(), this->GetStrideCols());
   }
}

//______________________________________________________________________________
template <typename Architecture_t>
auto TMaxPoolLayer<Architecture_t>::Backward(std::vector<Matrix_t> &gradients_backward,
                                             const std::vector<Matrix_t> & /*activations_backward*/,
                                             std::vector<Matrix_t> & /*inp1*/, std::vector<Matrix_t> &
                                             /*inp2*/) -> void
{
   Architecture_t::MaxPoolLayerBackward(gradients_backward, this->GetActivationGradients(), indexMatrix,
                                        this->GetBatchSize(), this->GetDepth(), this->GetNLocalViews());
}

//______________________________________________________________________________
template <typename Architecture_t>
auto TMaxPoolLayer<Architecture_t>::Print() const -> void
{
   std::cout << " POOL Layer: \t";
   std::cout << "( W = " << this->GetWidth() << " , ";
   std::cout << " H = " << this->GetHeight() << " , ";
   std::cout << " D = " << this->GetDepth() << " ) ";

   std::cout << "\t Frame ( W = " << this->GetFrameWidth() << " , ";
   std::cout << " H = " << this->GetFrameHeight() << " ) ";

   if (this->GetOutput().size() > 0) {
      std::cout << "\tOutput = ( " << this->GetOutput().size() << " , " << this->GetOutput()[0].GetNrows() << " , " << this->GetOutput()[0].GetNcols() << " ) ";
   }
   std::cout << std::endl;
}

//______________________________________________________________________________
template <typename Architecture_t>
void TMaxPoolLayer<Architecture_t>::AddWeightsXMLTo(void *parent)
{
   auto layerxml = gTools().xmlengine().NewChild(parent, 0, "MaxPoolLayer");

   // write  maxpool layer info
   gTools().xmlengine().NewAttr(layerxml, 0, "FrameHeight", gTools().StringFromInt(this->GetFrameHeight()));
   gTools().xmlengine().NewAttr(layerxml, 0, "FrameWidth", gTools().StringFromInt(this->GetFrameWidth()));
   gTools().xmlengine().NewAttr(layerxml, 0, "StrideRows", gTools().StringFromInt(this->GetStrideRows()));
   gTools().xmlengine().NewAttr(layerxml, 0, "StrideCols", gTools().StringFromInt(this->GetStrideCols()));

}

//______________________________________________________________________________
template <typename Architecture_t>
void TMaxPoolLayer<Architecture_t>::ReadWeightsFromXML(void * /*parent */)
{
   // all info is read before - nothing to do 
}

template <typename Architecture_t>
size_t TMaxPoolLayer<Architecture_t>::calculateDimension(int imgDim, int fltDim, int stride)
{
   Scalar_t dimension = ((imgDim - fltDim) / stride) + 1;
   if (!isInteger(dimension) || dimension <= 0) {
      Fatal("calculateDimension", "Not compatible hyper parameters for layer - (imageDim, filterDim, stride) %d , %d , %d",
            imgDim, fltDim, stride);
   }

   return (size_t)dimension;
}

template <typename Architecture_t>
size_t TMaxPoolLayer<Architecture_t>::calculateNLocalViews(int inputHeight, int filterHeight, int strideRows,
                                                           int inputWidth, int filterWidth, int strideCols)
{
   int height = calculateDimension(inputHeight, filterHeight, strideRows);
   int width = calculateDimension(inputWidth, filterWidth, strideCols);

   return height * width;
}

} // namespace CNN
} // namespace DNN
} // namespace TMVA

#endif
