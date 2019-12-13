// @(#)root/tmva/tmva/dnn:$Id$
// Author: Ashish Kshirsagar

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : TTransConvLayer                                                       *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Transpose Convolutional Deep Neural Network Layer                         *
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

#ifndef TMVA_CNN_CONVLAYER
#define TMVA_CNN_CONVLAYER

#include "TMatrix.h"

#include "TMVA/DNN/GeneralLayer.h"
#include "TMVA/DNN/Functions.h"

#include <vector>
#include <iostream>

namespace TMVA {
namespace DNN {
namespace CNN {

template <typename Architecture_t>
class TTransConvLayer : public VGeneralLayer<Architecture_t> {
public:
   using Matrix_t = typename Architecture_t::Matrix_t;
   using Scalar_t = typename Architecture_t::Scalar_t;

private:
   /* Calculate the output dimension of the transpose convolutional layer */
   size_t calculateDimension(size_t imgDim, size_t fltDim, size_t padding, size_t stride);

   /* Calculate the number of pixels in a single receptive field */
   size_t inline calculateNLocalViewPixels(size_t depth, size_t height, size_t width) { return depth * height * width; }

   /* Calculate the number of receptive fields in an image given the filter and image sizes */
   size_t calculateNLocalViews(size_t inputHeight, size_t filterHeight, size_t paddingHeight, size_t strideRows,
                               size_t inputWidth, size_t filterWidth, size_t paddingWidth, size_t strideCols);

protected:
   size_t fFilterDepth;  ///< The depth of the filter.
   size_t fFilterHeight; ///< The height of the filter.
   size_t fFilterWidth;  ///< The width of the filter.

   size_t fStrideRows; ///< The number of row pixels to slid the filter each step.
   size_t fStrideCols; ///< The number of column pixels to slid the filter each step.

   size_t fNLocalViewPixels; ///< The number of pixels in one local image view.
   size_t fNLocalViews;      ///< The number of local views in one image.

   Scalar_t fDropoutProbability; ///< Probability that an input is active.

private:
   size_t fPaddingHeight; ///< The number of zero layers added top and bottom of the input.
   size_t fPaddingWidth;  ///< The number of zero layers left and right of the input.

   std::vector<Matrix_t> fDerivatives; ///< First fDerivatives of the activations of this layer.

   std::vector<int> fBackwardIndices;  ///< Vector of indices used for a fast Im2Col in backward pass

   EActivationFunction fF; ///< Activation function of the layer.
   ERegularization fReg;   ///< The regularization method.
   Scalar_t fWeightDecay;  ///< The weight decay.

   std::vector<Matrix_t> fForwardMatrices; ///< Vector of matrices used for speeding-up the forward pass.

public:
   /*! Constructor. */
   TTransConvLayer(size_t BatchSize, size_t InputDepth, size_t InputHeight, size_t InputWidth, size_t Depth, EInitialization Init,
              size_t FilterHeight, size_t FilterWidth, size_t StrideRows, size_t StrideCols, size_t PaddingHeight,
              size_t PaddingWidth, Scalar_t DropoutProbability, EActivationFunction f, ERegularization Reg,
              Scalar_t WeightDecay);

   /*! Copy the trans conv layer provided as a pointer */
   TTransConvLayer(TTransConvLayer<Architecture_t> *layer);

   /*! Copy constructor. */
   TTransConvLayer(const TTransConvLayer &);

   /*! Destructor. */
   ~TTransConvLayer();

   /*! Computes activation of the layer for the given input. The input
   * must be in 3D tensor form with the different matrices corresponding to
   * different events in the batch. Computes activations as well as
   * the first partial derivative of the activation function at those
   * activations. */
   void Forward(std::vector<Matrix_t> &input, bool applyDropout = false);

   /*! Compute weight, bias and activation gradients. Uses the precomputed
    *  first partial derviatives of the activation function computed during
    *  forward propagation and modifies them. Must only be called directly
    *  at the corresponding call to Forward(...). */
   void Backward(std::vector<Matrix_t> &gradients_backward, const std::vector<Matrix_t> &activations_backward,
                 std::vector<Matrix_t> &inp1, std::vector<Matrix_t> &inp2);

   /*! Prints the info about the layer. */
   void Print() const;

   /*! Writes the information and the weights about the layer in an XML node. */
   virtual void AddWeightsXMLTo(void *parent);

   /*! Read the information and the weights about the layer from XML node. */
   virtual void ReadWeightsFromXML(void *parent);

   /*! Getters */
   size_t GetFilterDepth() const { return fFilterDepth; }
   size_t GetFilterHeight() const { return fFilterHeight; }
   size_t GetFilterWidth() const { return fFilterWidth; }

   size_t GetStrideRows() const { return fStrideRows; }
   size_t GetStrideCols() const { return fStrideCols; }

   size_t GetPaddingHeight() const { return fPaddingHeight; }
   size_t GetPaddingWidth() const { return fPaddingWidth; }

   size_t GetNLocalViewPixels() const { return fNLocalViewPixels; }
   size_t GetNLocalViews() const { return fNLocalViews; }

   Scalar_t GetDropoutProbability() const { return fDropoutProbability; }

   const std::vector<Matrix_t> &GetDerivatives() const { return fDerivatives; }
   std::vector<Matrix_t> &GetDerivatives() { return fDerivatives; }

   Matrix_t &GetDerivativesAt(size_t i) { return fDerivatives[i]; }
   const Matrix_t &GetDerivativesAt(size_t i) const { return fDerivatives[i]; }

   const std::vector<Matrix_t> &GetForwardMatrices() const { return fForwardMatrices; }
   std::vector<Matrix_t> &GetForwardMatrices() { return fForwardMatrices; }

   EActivationFunction GetActivationFunction() const { return fF; }
   ERegularization GetRegularization() const { return fReg; }
   Scalar_t GetWeightDecay() const { return fWeightDecay; }
};

typedef struct TTransConvParams {

public:
   size_t batchSize; ///< Batch size used for training and evaluation

   size_t inputDepth;  ///< The depth of the previous layer or input.
   size_t inputHeight; ///< The height of the previous layer or input.
   size_t inputWidth;  ///< The width of the previous layer or input.

   size_t numberFilters; ///< The number of the filters, which is equal to the output's depth.
   size_t filterHeight; ///< The height of the filter.
   size_t filterWidth;  ///< The width of the filter.

   size_t strideRows; ///< The number of row pixels to slid the filter each step.
   size_t strideCols; ///< The number of column pixels to slid the filter each step.
   size_t paddingHeight; ///< The number of zero layers added top and bottom of the input.
   size_t paddingWidth;  ///< The number of zero layers left and right of the input.

   TTransConvParams(size_t _batchSize, size_t _inputDepth, size_t _inputHeight, size_t _inputWidth, size_t _numberFilters,
               size_t _filterHeight, size_t _filterWidth, size_t _strideRows, size_t _strideCols,
               size_t _paddingHeight, size_t _paddingWidth)
           : batchSize(_batchSize), inputDepth(_inputDepth), inputHeight(_inputHeight), inputWidth(_inputWidth),
             numberFilters(_numberFilters), filterHeight(_filterHeight), filterWidth(_filterWidth),
             strideRows(_strideRows), strideCols(_strideCols), paddingHeight(_paddingHeight),
             paddingWidth(_paddingWidth)
   {}
} TTransConvParams;


//
//
//  Trans Conv Layer Class - Implementation
//______________________________________________________________________________
template <typename Architecture_t>
TTransConvLayer<Architecture_t>::TTransConvLayer(size_t batchSize, size_t inputDepth, size_t inputHeight, size_t inputWidth,
                                       size_t depth, EInitialization init, size_t filterHeight, size_t filterWidth,
                                       size_t strideRows, size_t strideCols, size_t paddingHeight, size_t paddingWidth,
                                       Scalar_t dropoutProbability, EActivationFunction f, ERegularization reg,
                                       Scalar_t weightDecay)
   : VGeneralLayer<Architecture_t>(batchSize, inputDepth, inputHeight, inputWidth, depth,
                                   calculateDimension(inputHeight, filterHeight, paddingHeight, strideRows),
                                   calculateDimension(inputWidth, filterWidth, paddingWidth, strideCols),
                                   1, depth, calculateNLocalViewPixels(inputDepth, filterHeight, filterWidth),
                                   1, depth, 1, batchSize, depth,
                                   calculateNLocalViews(inputHeight, filterHeight, paddingHeight, strideRows,
                                                        inputWidth, filterWidth, paddingWidth, strideCols),
                                   init),
     fFilterDepth(inputDepth), fFilterHeight(filterHeight), fFilterWidth(filterWidth), fStrideRows(strideRows),
     fStrideCols(strideCols), fNLocalViewPixels(calculateNLocalViewPixels(inputDepth, filterHeight, filterWidth)),
     fNLocalViews(calculateNLocalViews(inputHeight, filterHeight, paddingHeight, strideRows,
                                       inputWidth, filterWidth, paddingWidth, strideCols)),
     fDropoutProbability(dropoutProbability), fPaddingHeight(paddingHeight), fPaddingWidth(paddingWidth),
     fDerivatives(), fF(f), fReg(reg), fWeightDecay(weightDecay)
{
   /** Each element in the vector is a `T_Matrix` representing an event, therefore `vec.size() == batchSize`.
    *  Cells in these matrices are distributed in the following manner:
    *  Each row represents a single feature map, therefore we have `nRows == depth`.
    *  Each column represents a single pixel in that feature map, therefore we have `nCols == nLocalViews`.
    **/
   for (size_t i = 0; i < batchSize; i++) {
      fDerivatives.emplace_back(depth, fNLocalViews);
      fForwardMatrices.emplace_back(fNLocalViews, fNLocalViewPixels);
   }
   Architecture_t::PrepareInternals(fForwardMatrices);
}

//______________________________________________________________________________
template <typename Architecture_t>
TTransConvLayer<Architecture_t>::TTransConvLayer(TTransConvLayer<Architecture_t> *layer)
   : VGeneralLayer<Architecture_t>(layer), fFilterDepth(layer->GetFilterDepth()),
     fFilterHeight(layer->GetFilterHeight()), fFilterWidth(layer->GetFilterWidth()),
     fStrideRows(layer->GetStrideRows()), fStrideCols(layer->GetStrideCols()),
     fNLocalViewPixels(layer->GetNLocalViewPixels()), fNLocalViews(layer->GetNLocalViews()),
     fDropoutProbability(layer->GetDropoutProbability()), fPaddingHeight(layer->GetPaddingHeight()),
     fPaddingWidth(layer->GetPaddingWidth()), fF(layer->GetActivationFunction()),
     fReg(layer->GetRegularization()), fWeightDecay(layer->GetWeightDecay())
{
   size_t outputNSlices = (layer->GetDerivatives()).size();
   size_t outputNRows = 0;
   size_t outputNCols = 0;

   for (size_t i = 0; i < outputNSlices; i++) {
      outputNRows = (layer->GetDerivativesAt(i)).GetNrows();
      outputNCols = (layer->GetDerivativesAt(i)).GetNcols();
      fDerivatives.emplace_back(outputNRows, outputNCols);
      fForwardMatrices.emplace_back(layer->GetNLocalViews(), layer->GetNLocalViewPixels());
   }
}

//______________________________________________________________________________
template <typename Architecture_t>
TTransConvLayer<Architecture_t>::TTransConvLayer(const TTransConvLayer &transConvLayer)
   : VGeneralLayer<Architecture_t>(transConvLayer), fFilterDepth(transConvLayer.fFilterDepth),
     fFilterHeight(transConvLayer.fFilterHeight), fFilterWidth(transConvLayer.fFilterWidth), fStrideRows(transConvLayer.fStrideRows),
     fStrideCols(transConvLayer.fStrideCols), fNLocalViewPixels(transConvLayer.fNLocalViewPixels),
     fNLocalViews(transConvLayer.fNLocalViews), fDropoutProbability(transConvLayer.fDropoutProbability),
     fPaddingHeight(transConvLayer.fPaddingHeight), fPaddingWidth(transConvLayer.fPaddingWidth),  fF(transConvLayer.fF),
     fReg(transConvLayer.fReg), fWeightDecay(transConvLayer.fWeightDecay)
{
   size_t outputNSlices = transConvLayer.fDerivatives.size();
   size_t outputNRows = transConvLayer.GetDerivativesAt(0).GetNrows();
   size_t outputNCols = transConvLayer.GetDerivativesAt(0).GetNcols();

   for (size_t i = 0; i < outputNSlices; i++) {
      fDerivatives.emplace_back(outputNRows, outputNCols);
      fForwardMatrices.emplace_back(transConvLayer.fNLocalViews, transConvLayer.fNLocalViewPixels);
   }
}

//______________________________________________________________________________
template <typename Architecture_t>
TTransConvLayer<Architecture_t>::~TTransConvLayer()
{
}

//______________________________________________________________________________
template <typename Architecture_t>
auto TTransConvLayer<Architecture_t>::Forward(std::vector<Matrix_t> &input, bool /*applyDropout*/) -> void
{
   TTransConvParams params(this->GetBatchSize(), this->GetInputDepth(), this->GetInputHeight(), this->GetInputWidth(),
                      this->GetDepth(), this->GetFilterHeight(), this->GetFilterWidth(),
                      this->GetStrideRows(), this->GetStrideCols(), this->GetPaddingHeight(), this->GetPaddingWidth());

   R__ASSERT( input.size() > 0);
   Architecture_t::TransConvLayerForward(this->GetOutput(), this->GetDerivatives(), input, this->GetWeightsAt(0),
                                    this->GetBiasesAt(0), params, this->GetActivationFunction(),
                                    this->GetForwardMatrices());

#if 0
   // in printciple I could make the indices data member of the class
   Matrix_t inputTr(this->GetNLocalViews(), this->GetNLocalViewPixels());
   //Matrix_t inputTr2(this->GetNLocalViews(), this->GetNLocalViewPixels());
   std::vector<int> vIndices(inputTr.GetNrows() * inputTr.GetNcols() );
   R__ASSERT( input.size() > 0);
   Architecture_t::Im2colIndices(vIndices, input[0], this->GetNLocalViews(), this->GetInputHeight(), this->GetInputWidth(), this->GetFilterHeight(),
                             this->GetFilterWidth(), this->GetStrideRows(), this->GetStrideCols(),
                             this->GetPaddingHeight(), this->GetPaddingWidth());
   // batch size loop
   for (size_t i = 0; i < this->GetBatchSize(); i++) {

      if (applyDropout && (this->GetDropoutProbability() != 1.0)) {
         Architecture_t::Dropout(input[i], this->GetDropoutProbability());
      }

      inputTr.Zero();
      //inputTr2.Zero();
      // Architecture_t::Im2col(inputTr2, input[i], this->GetInputHeight(), this->GetInputWidth(), this->GetFilterHeight(),
      //                         this->GetFilterWidth(), this->GetStrideRows(), this->GetStrideCols(),
      //                         this->GetPaddingHeight(), this->GetPaddingWidth());
      Architecture_t::Im2colFast(inputTr, input[i], vIndices);
      // bool diff = false;
      // for (int j = 0; j < inputTr.GetNrows(); ++j) {
      //    for (int k = 0; k < inputTr.GetNcols(); ++k) {
      //       if ( inputTr2(j,k) != inputTr(j,k) ) {
      //          diff = true;
      //          std::cout <<  "different im2col for " << j << " , " << k << "  " << inputTr(j,k) << "  shoud be " << inputTr2(j,k) << std::endl;
      //       }
      //    }
      // }
      // if (diff) {
      //    std::cout << "TransConvLayer:: Different Im2Col for batch " << i  << std::endl;
      //    printf("Layer parameters : %d x %d , filter %d x %d , stride %d %d , pad %d %d \n",this->GetInputHeight(), this->GetInputWidth(), this->GetFilterHeight(),
      //                         this->GetFilterWidth(), this->GetStrideRows(), this->GetStrideCols(),
      //           this->GetPaddingHeight(), this->GetPaddingWidth() );
      //    // TMVA_DNN_PrintTCpuMatrix(inputTr);
      //    // TMVA_DNN_PrintTCpuMatrix(inputTr2);
      // }
      // R__ASSERT(!diff);
      Architecture_t::MultiplyTranspose(this->GetOutputAt(i), this->GetWeightsAt(0), inputTr);
      Architecture_t::AddConvBiases(this->GetOutputAt(i), this->GetBiasesAt(0));

      evaluateDerivative<Architecture_t>(this->GetDerivativesAt(i), fF, this->GetOutputAt(i));
      evaluate<Architecture_t>(this->GetOutputAt(i), fF);
   }
#endif
}

//______________________________________________________________________________
template <typename Architecture_t>
auto TTransConvLayer<Architecture_t>::Backward(std::vector<Matrix_t> &gradients_backward,
                                          const std::vector<Matrix_t> &activations_backward,
                                          std::vector<Matrix_t> & /*inp1*/, std::vector<Matrix_t> &
                                          /*inp2*/) -> void
{
   Architecture_t::TransConvLayerBackward(
      gradients_backward, this->GetWeightGradientsAt(0), this->GetBiasGradientsAt(0), this->GetDerivatives(),
      this->GetActivationGradients(), this->GetWeightsAt(0), activations_backward, this->GetBatchSize(),
      this->GetInputHeight(), this->GetInputWidth(), this->GetDepth(), this->GetHeight(), this->GetWidth(),
      this->GetFilterDepth(), this->GetFilterHeight(), this->GetFilterWidth(), this->GetNLocalViews());

   addRegularizationGradients<Architecture_t>(this->GetWeightGradientsAt(0), this->GetWeightsAt(0),
                                              this->GetWeightDecay(), this->GetRegularization());
}

//______________________________________________________________________________
template <typename Architecture_t>
auto TTransConvLayer<Architecture_t>::Print() const -> void
{
   std::cout << " TRANS CONV LAYER: \t";
   std::cout << "( W = " << this->GetWidth() << " , ";
   std::cout << " H = " << this->GetHeight() << " , ";
   std::cout << " D = " << this->GetDepth() << " ) ";

   std::cout << "\t Filter ( W = " << this->GetFilterWidth() << " , ";
   std::cout << " H = " << this->GetFilterHeight() << " ) ";
   //std::cout << "\t Local Views = " << this->GetNLocalViews()  << " " ;
   if (this->GetOutput().size() > 0) {
      std::cout << "\tOutput = ( " << this->GetOutput().size() << " , " << this->GetOutput()[0].GetNrows() << " , " << this->GetOutput()[0].GetNcols() << " ) ";
   }
   std::vector<std::string> activationNames = { "Identity","Relu","Sigmoid","Tanh","SymmRelu","SoftSign","Gauss" };
   std::cout << "\t Activation Function = ";
   std::cout << activationNames[ static_cast<int>(fF) ] << std::endl;
}

//______________________________________________________________________________
template <typename Architecture_t>
void TTransConvLayer<Architecture_t>::AddWeightsXMLTo(void *parent)
{
   auto layerxml = gTools().xmlengine().NewChild(parent, 0, "TransConvLayer");

   gTools().xmlengine().NewAttr(layerxml, 0, "Depth", gTools().StringFromInt(this->GetDepth()));
   gTools().xmlengine().NewAttr(layerxml, 0, "FilterHeight", gTools().StringFromInt(this->GetFilterHeight()));
   gTools().xmlengine().NewAttr(layerxml, 0, "FilterWidth", gTools().StringFromInt(this->GetFilterWidth()));
   gTools().xmlengine().NewAttr(layerxml, 0, "StrideRows", gTools().StringFromInt(this->GetStrideRows()));
   gTools().xmlengine().NewAttr(layerxml, 0, "StrideCols", gTools().StringFromInt(this->GetStrideCols()));
   gTools().xmlengine().NewAttr(layerxml, 0, "PaddingHeight", gTools().StringFromInt(this->GetPaddingHeight()));
   gTools().xmlengine().NewAttr(layerxml, 0, "PaddingWidth", gTools().StringFromInt(this->GetPaddingWidth()));


   int activationFunction = static_cast<int>(this -> GetActivationFunction());
   gTools().xmlengine().NewAttr(layerxml, 0, "ActivationFunction",
                                TString::Itoa(activationFunction, 10));

   // write weights and bias matrix
   this->WriteMatrixToXML(layerxml, "Weights", this -> GetWeightsAt(0));
   this->WriteMatrixToXML(layerxml, "Biases",  this -> GetBiasesAt(0));

}

//______________________________________________________________________________
template <typename Architecture_t>
void TTransConvLayer<Architecture_t>::ReadWeightsFromXML(void *parent)
{
   // read weights and biases
   // the meta information is read before because it is needed before creating the Conv layer
   this->ReadMatrixXML(parent,"Weights", this -> GetWeightsAt(0));
   this->ReadMatrixXML(parent,"Biases", this -> GetBiasesAt(0));
}

template <typename Architecture_t>
size_t TTransConvLayer<Architecture_t>::calculateDimension(size_t imgDim, size_t fltDim, size_t padding, size_t stride)
{
   size_t temp = imgDim - fltDim + 2 * padding;
   if (temp % stride || temp + stride <= 0) {
      Fatal("calculateDimension", "Not compatible hyper parameters for layer - (imageDim, filterDim, padding, stride) "
            "%zu, %zu, %zu, %zu", imgDim, fltDim, padding, stride);
   }
   return temp / stride + 1;
}

template <typename Architecture_t>
size_t TTransConvLayer<Architecture_t>::calculateNLocalViews(size_t inputHeight, size_t filterHeight, size_t paddingHeight,
                                                        size_t strideRows, size_t inputWidth, size_t filterWidth,
                                                        size_t paddingWidth, size_t strideCols)
{
    int height = calculateDimension(inputHeight, filterHeight, paddingHeight, strideRows);
    int width = calculateDimension(inputWidth, filterWidth, paddingWidth, strideCols);

    return height * width;
}

} // namespace CNN
} // namespace DNN
} // namespace TMVA

#endif
