// @(#)root/tmva/tmva/dnn:$Id$
// Author: Vladimir Ilievski

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : TGeneralLayer                                                         *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      General Deep Neural Network Layer                                         *
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

#ifndef TMVA_DNN_GENERALLAYER
#define TMVA_DNN_GENERALLAYER

#include <sstream>
#include <limits>
#include <vector>
#include <string>

// for xml
#include "TMVA/Tools.h"
#include "TError.h"   // for R__ASSERT

#include "TMVA/DNN/Functions.h"

namespace TMVA {
namespace DNN {

/** \class VGeneralLayer
    Generic General Layer class.

    This class represents the general class for all layers in the Deep Learning
    Module.
 */
template <typename Architecture_t>
class VGeneralLayer {

   using Tensor_t = typename Architecture_t::Tensor_t;
   using Matrix_t = typename Architecture_t::Matrix_t;
   using Scalar_t = typename Architecture_t::Scalar_t;


protected:
   size_t fBatchSize; ///< Batch size used for training and evaluation

   size_t fInputDepth;  ///< The depth of the previous layer or input.
   size_t fInputHeight; ///< The height of the previous layer or input.
   size_t fInputWidth;  ///< The width of the previous layer or input.

   size_t fDepth;  ///< The depth of the layer.
   size_t fHeight; ///< The height of the layer.
   size_t fWidth;  ///< The width of this layer.

   bool fIsTraining; ///< Flag indicating the mode

   std::vector<Matrix_t> fWeights; ///< The weights associated to the layer.
   std::vector<Matrix_t> fBiases;  ///< The biases associated to the layer.

   std::vector<Matrix_t> fWeightGradients; ///< Gradients w.r.t. the weights of the layer.
   std::vector<Matrix_t> fBiasGradients;   ///< Gradients w.r.t. the bias values of the layer.

   Tensor_t fOutput;              ///< Activations of this layer.
   Tensor_t fActivationGradients; ///< Gradients w.r.t. the activations of this layer.

   EInitialization fInit; ///< The initialization method.

public:
   /*! Constructor */
   VGeneralLayer(size_t BatchSize, size_t InputDepth, size_t InputHeight, size_t InputWidth, size_t Depth,
                 size_t Height, size_t Width, size_t WeightsNSlices, size_t WeightsNRows, size_t WeightsNCols,
                 size_t BiasesNSlices, size_t BiasesNRows, size_t BiasesNCols, size_t OutputNSlices, size_t OutputNRows,
                 size_t OutputNCols, EInitialization Init);

   /*! General Constructor with different weights dimension */
   VGeneralLayer(size_t BatchSize, size_t InputDepth, size_t InputHeight, size_t InputWidth, size_t Depth,
                 size_t Height, size_t Width, size_t WeightsNSlices, std::vector<size_t> WeightsNRows,
                 std::vector<size_t> WeightsNCols, size_t BiasesNSlices, std::vector<size_t> BiasesNRows,
                 std::vector<size_t> BiasesNCols, size_t OutputNSlices, size_t OutputNRows, size_t OutputNCols,
                 EInitialization Init);

   /*! Copy the layer provided as a pointer */
   VGeneralLayer(VGeneralLayer<Architecture_t> *layer);

   /*! Copy Constructor */
   VGeneralLayer(const VGeneralLayer &);

   /*! Virtual Destructor. */
   virtual ~VGeneralLayer();

   /*! Initialize the weights and biases according to the given initialization method. */
   virtual void Initialize();

   /*! Computes activation of the layer for the given input. The input
    * must be in 3D tensor form with the different matrices corresponding to
    * different events in the batch.  */
   virtual void Forward(Tensor_t &input, bool applyDropout = false) = 0;

   /*! Backpropagates the error. Must only be called directly at the corresponding
    *  call to Forward(...). */
   virtual void Backward(Tensor_t &gradients_backward, const Tensor_t &activations_backward ) = 0;
   /////                      std::vector<Matrix_t> &inp1, std::vector<Matrix_t> &inp2) = 0;

    /*! Reset some training flags after a loop on all batches
       Some layer (e.g. batchnormalization) might need to implement the function in case some operations
       are needed after looping an all batches                                                 */
   virtual void ResetTraining() {}

   /*! Updates the weights and biases, given the learning rate */
   void Update(const Scalar_t learningRate);

   /*! Updates the weights, given the gradients and the learning rate, */
   void UpdateWeights(const std::vector<Matrix_t> &weightGradients, const Scalar_t learningRate);

   /*! Updates the biases, given the gradients and the learning rate. */
   void UpdateBiases(const std::vector<Matrix_t> &biasGradients, const Scalar_t learningRate);

   /*! Updates the weight gradients, given some other weight gradients and learning rate. */
   void UpdateWeightGradients(const std::vector<Matrix_t> &weightGradients, const Scalar_t learningRate);

   /*! Updates the bias gradients, given some other weight gradients and learning rate. */
   void UpdateBiasGradients(const std::vector<Matrix_t> &biasGradients, const Scalar_t learningRate);

   /*! Copies the weights provided as an input.  */
   void CopyWeights(const std::vector<Matrix_t> &otherWeights);

   /*! Copies the biases provided as an input. */
   void CopyBiases(const std::vector<Matrix_t> &otherBiases);

   /*! Copy all trainable weight and biases from another equivalent layer but with different architecture
       The function can copy also extra parameters in addition to weights and biases if they are return
       by the function GetExtraLayerParameters */
   template <typename Arch>
   void CopyParameters(const VGeneralLayer<Arch> &layer);

   /*! Prints the info about the layer. */
   virtual void Print() const = 0;

   /*! Writes the information and the weights about the layer in an XML node. */
   virtual void AddWeightsXMLTo(void *parent) = 0;

   /*! Read the information and the weights about the layer from XML node. */
   virtual void ReadWeightsFromXML(void *parent) = 0;

   /*! Set Dropout probability. Reimplemented for layesrs supporting droput */
   virtual void SetDropoutProbability(Scalar_t ) {}

   /*! Getters */
   size_t GetBatchSize() const { return fBatchSize; }
   size_t GetInputDepth() const { return fInputDepth; }
   size_t GetInputHeight() const { return fInputHeight; }
   size_t GetInputWidth() const { return fInputWidth; }
   size_t GetDepth() const { return fDepth; }
   size_t GetHeight() const { return fHeight; }
   size_t GetWidth() const { return fWidth; }
   bool IsTraining() const { return fIsTraining; }

   const std::vector<Matrix_t> &GetWeights() const { return fWeights; }
   std::vector<Matrix_t> &GetWeights() { return fWeights; }

   const Matrix_t &GetWeightsAt(size_t i) const { return fWeights[i]; }
   Matrix_t &GetWeightsAt(size_t i) { return fWeights[i]; }

   const std::vector<Matrix_t> &GetBiases() const { return fBiases; }
   std::vector<Matrix_t> &GetBiases() { return fBiases; }

   const Matrix_t &GetBiasesAt(size_t i) const { return fBiases[i]; }
   Matrix_t &GetBiasesAt(size_t i) { return fBiases[i]; }

   const std::vector<Matrix_t> &GetWeightGradients() const { return fWeightGradients; }
   std::vector<Matrix_t> &GetWeightGradients() { return fWeightGradients; }

   const Matrix_t &GetWeightGradientsAt(size_t i) const { return fWeightGradients[i]; }
   Matrix_t &GetWeightGradientsAt(size_t i) { return fWeightGradients[i]; }

   const std::vector<Matrix_t> &GetBiasGradients() const { return fBiasGradients; }
   std::vector<Matrix_t> &GetBiasGradients() { return fBiasGradients; }

   const Matrix_t &GetBiasGradientsAt(size_t i) const { return fBiasGradients[i]; }
   Matrix_t &GetBiasGradientsAt(size_t i) { return fBiasGradients[i]; }

   const Tensor_t &GetOutput() const { return fOutput; }
   Tensor_t &GetOutput() { return fOutput; }

   const Tensor_t &GetActivationGradients() const { return fActivationGradients; }
   Tensor_t &GetActivationGradients() { return fActivationGradients; }

   Matrix_t GetOutputAt(size_t i) { return fOutput.At(i).GetMatrix(); }
   const Matrix_t &GetOutputAt(size_t i) const { return fOutput.At(i).GetMatrix(); }

   Matrix_t GetActivationGradientsAt(size_t i) { return fActivationGradients.At(i).GetMatrix(); }
   const Matrix_t &GetActivationGradientsAt(size_t i) const { return fActivationGradients.At(i).GetMatrix(); }

   // function to retrieve additional layer parameters which are learned during training but they are not weights
   // an example are the mean and std of batch normalization layer
   virtual std::vector<Matrix_t> GetExtraLayerParameters() const { return std::vector<Matrix_t>(); }
   // same thing but to set these extra parameters
   virtual void  SetExtraLayerParameters(const std::vector<Matrix_t> & ) {}

   EInitialization GetInitialization() const { return fInit; }

   /*! Setters */
   void SetBatchSize(size_t batchSize) { fBatchSize = batchSize; }
   void SetInputDepth(size_t inputDepth) { fInputDepth = inputDepth; }
   void SetInputHeight(size_t inputHeight) { fInputHeight = inputHeight; }
   void SetInputWidth(size_t inputWidth) { fInputWidth = inputWidth; }
   void SetDepth(size_t depth) { fDepth = depth; }
   void SetHeight(size_t height) { fHeight = height; }
   void SetWidth(size_t width) { fWidth = width; }
   void SetIsTraining(bool isTraining) { fIsTraining = isTraining; }

   /// helper functions for XML
   void WriteTensorToXML( void * node, const char * name, const std::vector<Matrix_t> & tensor);
   void WriteMatrixToXML( void * node, const char * name, const Matrix_t & matrix);

   void ReadMatrixXML( void * node, const char * name, Matrix_t & matrix);

};

//
//
//  The General Layer Class - Implementation
//_________________________________________________________________________________________________
template <typename Architecture_t>
VGeneralLayer<Architecture_t>::VGeneralLayer(size_t batchSize, size_t inputDepth, size_t inputHeight, size_t inputWidth,
                                             size_t depth, size_t height, size_t width, size_t weightsNSlices,
                                             size_t weightsNRows, size_t weightsNCols, size_t biasesNSlices,
                                             size_t biasesNRows, size_t biasesNCols, size_t outputNSlices,
                                             size_t outputNRows, size_t outputNCols, EInitialization init)
   :  fBatchSize(batchSize), fInputDepth(inputDepth), fInputHeight(inputHeight), fInputWidth(inputWidth), fDepth(depth),
      fHeight(height), fWidth(width), fIsTraining(true), fWeights(), fBiases(), fWeightGradients(), fBiasGradients(),
      fOutput( outputNSlices, outputNRows, outputNCols ),
      fActivationGradients( outputNSlices, outputNRows, outputNCols ),
      fInit(init)
{

   for (size_t i = 0; i < weightsNSlices; i++) {
      fWeights.emplace_back(weightsNRows, weightsNCols);
      fWeightGradients.emplace_back(weightsNRows, weightsNCols);
   }

   for (size_t i = 0; i < biasesNSlices; i++) {
      fBiases.emplace_back(biasesNRows, biasesNCols);
      fBiasGradients.emplace_back(biasesNRows, biasesNCols);
   }
}

//_________________________________________________________________________________________________
template <typename Architecture_t>
VGeneralLayer<Architecture_t>::VGeneralLayer(size_t batchSize, size_t inputDepth, size_t inputHeight, size_t inputWidth,
                                             size_t depth, size_t height, size_t width, size_t weightsNSlices,
                                             std::vector<size_t> weightsNRows, std::vector<size_t> weightsNCols,
                                             size_t biasesNSlices, std::vector<size_t> biasesNRows,
                                             std::vector<size_t> biasesNCols, size_t outputNSlices, size_t outputNRows,
                                             size_t outputNCols, EInitialization init)
   :  fBatchSize(batchSize), fInputDepth(inputDepth), fInputHeight(inputHeight), fInputWidth(inputWidth), fDepth(depth),
      fHeight(height), fWidth(width), fIsTraining(true), fWeights(), fBiases(), fWeightGradients(), fBiasGradients(),
      fOutput( outputNSlices, outputNRows, outputNCols ),
      fActivationGradients( outputNSlices, outputNRows, outputNCols ),
      fInit(init)
{
   // add constructor for weights with different shapes (e.g. in recurrent layers)
   for (size_t i = 0; i < weightsNSlices; i++) {
      fWeights.emplace_back(weightsNRows[i], weightsNCols[i]);
      fWeightGradients.emplace_back(weightsNRows[i], weightsNCols[i]);
   }

   for (size_t i = 0; i < biasesNSlices; i++) {
      fBiases.emplace_back(biasesNRows[i], biasesNCols[i]);
      fBiasGradients.emplace_back(biasesNRows[i], biasesNCols[i]);
   }

   // for (size_t i = 0; i < outputNSlices; i++) {
   //    fOutput.emplace_back(outputNRows, outputNCols);
   //    fActivationGradients.emplace_back(outputNRows, outputNCols);
   // }
}

//_________________________________________________________________________________________________
template <typename Architecture_t>
VGeneralLayer<Architecture_t>::VGeneralLayer(VGeneralLayer<Architecture_t> *layer)
   :  fBatchSize(layer->GetBatchSize()), fInputDepth(layer->GetInputDepth()), fInputHeight(layer->GetInputHeight()),
      fInputWidth(layer->GetInputWidth()), fDepth(layer->GetDepth()), fHeight(layer->GetHeight()),
      fWidth(layer->GetWidth()), fIsTraining(layer->IsTraining()), fWeights(), fBiases(), fWeightGradients(),
      fBiasGradients(),
      fOutput( layer->GetOutput().GetShape() ),   // construct from shape of other tensor
      fActivationGradients( layer->GetActivationGradients().GetShape() ),
      fInit(layer->GetInitialization() )
{
   // Constructor from another layer pointer of a different architecture
   size_t weightsNSlices = (layer->GetWeights()).size();
   size_t weightsNRows = 0;
   size_t weightsNCols = 0;

   for (size_t i = 0; i < weightsNSlices; i++) {
      weightsNRows = (layer->GetWeightsAt(i)).GetNrows();
      weightsNCols = (layer->GetWeightsAt(i)).GetNcols();

      fWeights.emplace_back(weightsNRows, weightsNCols);
      fWeightGradients.emplace_back(weightsNRows, weightsNCols);

      Architecture_t::Copy(fWeights[i], layer->GetWeightsAt(i));
   }

   size_t biasesNSlices = (layer->GetBiases()).size();
   size_t biasesNRows = 0;
   size_t biasesNCols = 0;

   for (size_t i = 0; i < biasesNSlices; i++) {
      biasesNRows = (layer->GetBiasesAt(i)).GetNrows();
      biasesNCols = (layer->GetBiasesAt(i)).GetNcols();

      fBiases.emplace_back(biasesNRows, biasesNCols);
      fBiasGradients.emplace_back(biasesNRows, biasesNCols);

      Architecture_t::Copy(fBiases[i], layer->GetBiasesAt(i));
   }
}

//_________________________________________________________________________________________________
template <typename Architecture_t>
VGeneralLayer<Architecture_t>::VGeneralLayer(const VGeneralLayer &layer)
   :  fBatchSize(layer.fBatchSize), fInputDepth(layer.fInputDepth), fInputHeight(layer.fInputHeight),
      fInputWidth(layer.fInputWidth), fDepth(layer.fDepth), fHeight(layer.fHeight), fWidth(layer.fWidth),
      fIsTraining(layer.fIsTraining), fWeights(), fBiases(), fWeightGradients(), fBiasGradients(),
      fOutput( layer.GetOutput() ),
      fActivationGradients( layer.GetActivationGradients() ),
      fInit( layer.GetInitialization())
{
   // copy constructor
   size_t weightsNSlices = layer.fWeights.size();
   size_t weightsNRows = 0;
   size_t weightsNCols = 0;

   for (size_t i = 0; i < weightsNSlices; i++) {
      weightsNRows = (layer.fWeights[i]).GetNrows();
      weightsNCols = (layer.fWeights[i]).GetNcols();

      fWeights.emplace_back(weightsNRows, weightsNCols);
      fWeightGradients.emplace_back(weightsNRows, weightsNCols);

      Architecture_t::Copy(fWeights[i], layer.fWeights[i]);
   }

   size_t biasesNSlices = layer.fBiases.size();
   size_t biasesNRows = 0;
   size_t biasesNCols = 0;

   for (size_t i = 0; i < biasesNSlices; i++) {
      biasesNRows = (layer.fBiases[i]).GetNrows();
      biasesNCols = (layer.fBiases[i]).GetNcols();

      fBiases.emplace_back(biasesNRows, biasesNCols);
      fBiasGradients.emplace_back(biasesNRows, biasesNCols);

      Architecture_t::Copy(fBiases[i], layer.fBiases[i]);
   }

   size_t outputNSlices = layer.fOutput.size();
   size_t outputNRows = 0;
   size_t outputNCols = 0;

   for (size_t i = 0; i < outputNSlices; i++) {
      outputNRows = (layer.fOutput[i]).GetNrows();
      outputNCols = (layer.fOutput[i]).GetNcols();

      fOutput.emplace_back(outputNRows, outputNCols);
      fActivationGradients.emplace_back(outputNRows, outputNCols);
   }
}

//_________________________________________________________________________________________________
template <typename Architecture_t>
VGeneralLayer<Architecture_t>::~VGeneralLayer()
{
   // Nothing to do here.
}

//_________________________________________________________________________________________________
template <typename Architecture_t>
auto VGeneralLayer<Architecture_t>::Initialize() -> void
{
   for (size_t i = 0; i < fWeights.size(); i++) {
      initialize<Architecture_t>(fWeights[i], this->GetInitialization());
      initialize<Architecture_t>(fWeightGradients[i], EInitialization::kZero);
   }

   for (size_t i = 0; i < fBiases.size(); i++) {
      initialize<Architecture_t>(fBiases[i], EInitialization::kZero);
      initialize<Architecture_t>(fBiasGradients[i], EInitialization::kZero);
   }
}

//_________________________________________________________________________________________________
template <typename Architecture_t>
auto VGeneralLayer<Architecture_t>::Update(const Scalar_t learningRate) -> void
{
   this->UpdateWeights(fWeightGradients, learningRate);
   this->UpdateBiases(fBiasGradients, learningRate);
}

//_________________________________________________________________________________________________
template <typename Architecture_t>
auto VGeneralLayer<Architecture_t>::UpdateWeights(const std::vector<Matrix_t> &weightGradients,
                                                  const Scalar_t learningRate) -> void
{
   for (size_t i = 0; i < fWeights.size(); i++) {
      Architecture_t::ScaleAdd(fWeights[i], weightGradients[i], -learningRate);
   }
}

//_________________________________________________________________________________________________
template <typename Architecture_t>
auto VGeneralLayer<Architecture_t>::UpdateBiases(const std::vector<Matrix_t> &biasGradients,
                                                 const Scalar_t learningRate) -> void
{
   for (size_t i = 0; i < fBiases.size(); i++) {
      Architecture_t::ScaleAdd(fBiases[i], biasGradients[i], -learningRate);
   }
}

//_________________________________________________________________________________________________
template <typename Architecture_t>
auto VGeneralLayer<Architecture_t>::UpdateWeightGradients(const std::vector<Matrix_t> &weightGradients,
                                                          const Scalar_t learningRate) -> void
{
   for (size_t i = 0; i < fWeightGradients.size(); i++) {
      Architecture_t::ScaleAdd(fWeightGradients[i], weightGradients[i], -learningRate);
   }
}

//_________________________________________________________________________________________________
template <typename Architecture_t>
auto VGeneralLayer<Architecture_t>::UpdateBiasGradients(const std::vector<Matrix_t> &biasGradients,
                                                        const Scalar_t learningRate) -> void
{
   for (size_t i = 0; i < fBiasGradients.size(); i++) {
      Architecture_t::ScaleAdd(fBiasGradients[i], biasGradients[i], -learningRate);
   }
}

//_________________________________________________________________________________________________
template <typename Architecture_t>
auto VGeneralLayer<Architecture_t>::CopyWeights(const std::vector<Matrix_t> &otherWeights) -> void
{

   for (size_t i = 0; i < fWeights.size(); i++) {
      Architecture_t::Copy(fWeights[i], otherWeights[i]);
   }
}

//_________________________________________________________________________________________________
template <typename Architecture_t>
auto VGeneralLayer<Architecture_t>::CopyBiases(const std::vector<Matrix_t> &otherBiases) -> void
{
   for (size_t i = 0; i < fBiases.size(); i++) {
      Architecture_t::Copy(fBiases[i], otherBiases[i]);
   }
}

//_________________________________________________________________________________________________
template <typename Architecture_t>
template <typename Arch>
void VGeneralLayer<Architecture_t>::CopyParameters(const VGeneralLayer<Arch> &layer)
{
   //assert(!std::is_same<Arch, Architecture_t>::value);
   // copy weights from a different arhcitecture- default generic implementation
   Architecture_t::CopyDiffArch(this->GetWeights(), layer.GetWeights());
   Architecture_t::CopyDiffArch(this->GetBiases(), layer.GetBiases());

   // copy also the additional layer parameters
   auto params = layer.GetExtraLayerParameters();
   if (params.size() > 0) {
      auto paramsToCopy = GetExtraLayerParameters();
      Architecture_t::CopyDiffArch(paramsToCopy, params );
      SetExtraLayerParameters(paramsToCopy);
   }
}

//_________________________________________________________________________________________________
template <typename Architecture_t>
auto VGeneralLayer<Architecture_t>::WriteTensorToXML(void * node, const char * name, const std::vector<Matrix_t> & tensor) -> void
{
   auto xmlengine = gTools().xmlengine();
   void* matnode = xmlengine.NewChild(node, 0, name);
   if (tensor.size() == 0) return;
   xmlengine.NewAttr(matnode,0,"Depth", gTools().StringFromInt(tensor.size()) );
   // assume same number of rows and columns for every matrix in std::vector
   xmlengine.NewAttr(matnode,0,"Rows", gTools().StringFromInt(tensor[0].GetNrows()) );
   xmlengine.NewAttr(matnode,0,"Columns", gTools().StringFromInt(tensor[0].GetNcols()) );
   std::stringstream s;
   for (size_t i = 0; i < tensor.size(); ++i) {
      auto & mat = tensor[i];
      for (Int_t row = 0; row < mat.GetNrows(); row++) {
         for (Int_t col = 0; col < mat.GetNcols(); col++) {
            // TString tmp = TString::Format( "%5.15e ", (mat)(row,col) );
            // s << tmp.Data();
            s << std::scientific << mat(row, col) << " ";
         }
      }
   }
   xmlengine.AddRawLine( matnode, s.str().c_str() );
}

//_________________________________________________________________________________________________
template <typename Architecture_t>
auto VGeneralLayer<Architecture_t>::WriteMatrixToXML(void * node, const char * name, const Matrix_t & matrix) -> void
{
   auto xmlengine = gTools().xmlengine();
   void* matnode = xmlengine.NewChild(node, 0, name);

   xmlengine.NewAttr(matnode,0,"Rows", gTools().StringFromInt(matrix.GetNrows()) );
   xmlengine.NewAttr(matnode,0,"Columns", gTools().StringFromInt(matrix.GetNcols()) );
   std::stringstream s;
   s.precision( std::numeric_limits<Scalar_t>::digits10 );
   size_t nrows = matrix.GetNrows();
   size_t ncols = matrix.GetNcols();
   for (size_t row = 0; row < nrows; row++) {
      for (size_t col = 0; col < ncols; col++) {
         //TString tmp = TString::Format( "%5.15e ", matrix(row,col) );
         s << std::scientific <<  matrix(row,col) << "  ";
      }
   }

   xmlengine.AddRawLine( matnode, s.str().c_str() );
}

//_________________________________________________________________________________________________
template <typename Architecture_t>
auto VGeneralLayer<Architecture_t>::ReadMatrixXML(void * node, const char * name, Matrix_t & matrix) -> void
{
   void *matrixXML = gTools().GetChild(node, name);
   size_t rows, cols;
   gTools().ReadAttr(matrixXML, "Rows", rows);
   gTools().ReadAttr(matrixXML, "Columns", cols);

   R__ASSERT((size_t) matrix.GetNrows() == rows);
   R__ASSERT((size_t) matrix.GetNcols() == cols);

   TMatrixT<Scalar_t> tmatrix(rows, cols);

   const char * matrixString = gTools().xmlengine().GetNodeContent(matrixXML);
   std::stringstream matrixStringStream(matrixString);

   for (size_t i = 0; i < rows; i++)
   {
      for (size_t j = 0; j < cols; j++)
      {
#ifndef R__HAS_TMVAGPU
         matrixStringStream >> tmatrix(i,j);
#else
         Scalar_t value;
         matrixStringStream >> value;
         tmatrix(i,j) = value;
#endif

      }
   }

   // copy from tmatrix to matrix
   Matrix_t tmp( tmatrix);
   Architecture_t::Copy(matrix, tmp);

}


template <typename Architecture>
auto debugTensor(const typename Architecture::Tensor_t & A, const std::string name = "tensor") -> void
{
   Architecture::PrintTensor(A,name);
}

} // namespace DNN
} // namespace TMVA

#endif
