
// Author: Vladimir Ilievski

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : TDenseLayer                                                           *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Dense Layer Class                                                         *
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

#ifndef TMVA_DNN_DENSELAYER
#define TMVA_DNN_DENSELAYER

#include "TMatrix.h"

#include "TMVA/DNN/GeneralLayer.h"
#include "TMVA/DNN/Functions.h"
#include "TMVA/DNN/CNN/ContextHandles.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

namespace TMVA {
namespace DNN {
/** \class TDenseLayer

Generic layer class.

This generic layer class represents a dense layer of a neural network with
a given width n and activation function f. The activation function of each
layer is given by \f$\mathbf{u} = \mathbf{W}\mathbf{x} + \boldsymbol{\theta}\f$.

In addition to the weight and bias matrices, each layer allocates memory
for its activations and the corresponding input tensor before evaluation  of
the activation function as well as the gradients of the weights and biases.

The layer provides member functions for the forward propagation of
activations through the given layer.
*/
template <typename Architecture_t>
class TDenseLayer : public VGeneralLayer<Architecture_t> {
public:

   using Scalar_t = typename Architecture_t::Scalar_t;
   using Matrix_t = typename Architecture_t::Matrix_t;
   using Tensor_t = typename Architecture_t::Tensor_t;

private:

   Tensor_t fInputActivation; /// output of GEMM and input to activation function
   Tensor_t fDerivatives; /// activation functgion gradient 

   Scalar_t fDropoutProbability; ///< Probability that an input is active.

   EActivationFunction fF; ///< Activation function of the layer.
   ERegularization fReg;   ///< The regularization method.
   Scalar_t fWeightDecay;  ///< The weight decay.

   typename Architecture_t::ActivationDescriptor_t fActivationDesc; // the descriptor for the activation function  

public:
   /*! Constructor */
   TDenseLayer(size_t BatchSize, size_t InputWidth, size_t Width, EInitialization init, Scalar_t DropoutProbability,
               EActivationFunction f, ERegularization reg, Scalar_t weightDecay);

   /*! Copy the dense layer provided as a pointer */
   TDenseLayer(TDenseLayer<Architecture_t> *layer);

   /*! Copy Constructor */
   TDenseLayer(const TDenseLayer &);

   /*! Destructor */
   ~TDenseLayer();

   /*! Compute activation of the layer for the given input. The input
    * must be in 3D tensor form with the different matrices corresponding to
    * different events in the batch. Computes activations as well as
    * the first partial derivative of the activation function at those
    * activations. */
   void Forward(Tensor_t &input, bool applyDropout = false);

   /*! Compute weight, bias and activation gradients. Uses the precomputed
    *  first partial derviatives of the activation function computed during
    *  forward propagation and modifies them. Must only be called directly
    *  a the corresponding call to Forward(...). */
   void Backward(Tensor_t &gradients_backward, const Tensor_t &activations_backward );
   ///              std::vector<Matrix_t> &inp1, std::vector<Matrix_t> &inp2);

   /*! Printing the layer info. */
   void Print() const;

   /*! Writes the information and the weights about the layer in an XML node. */
   virtual void AddWeightsXMLTo(void *parent);

   /*! Read the information and the weights about the layer from XML node. */
   virtual void ReadWeightsFromXML(void *parent);

   /*! Set dropout probabilities */
   virtual void SetDropoutProbability(Scalar_t dropoutProbability) { fDropoutProbability = dropoutProbability; }

   /*! Getters */
   Scalar_t GetDropoutProbability() const { return fDropoutProbability; }

   /* return output of Gemm before computing the activation function */ 
   const Tensor_t &GetInputActivation() const { return fInputActivation; }
   Tensor_t &GetInputActivation() { return fInputActivation; }
   
   EActivationFunction GetActivationFunction() const { return fF; }
   ERegularization GetRegularization() const { return fReg; }
   Scalar_t GetWeightDecay() const { return fWeightDecay; }
};

//
//
//  The Dense Layer Class - Implementation
//______________________________________________________________________________
template <typename Architecture_t>
TDenseLayer<Architecture_t>::TDenseLayer(size_t batchSize, size_t inputWidth, size_t width, EInitialization init,
                                         Scalar_t dropoutProbability, EActivationFunction f, ERegularization reg,
                                         Scalar_t weightDecay)
   :  VGeneralLayer<Architecture_t>(batchSize, 1, 1, inputWidth, 1, 1, width, 1, width, inputWidth, 1, width, 1, 1,
                                   batchSize, width, init),
      fInputActivation(), fDropoutProbability(dropoutProbability), fF(f), fReg(reg), fWeightDecay(weightDecay)
{
   // should be  {1, batchSize, width} but take from output
   fInputActivation = Tensor_t ( this->GetOutput().GetShape() );
   fDerivatives = Tensor_t ( this->GetOutput().GetShape() );

   Architecture_t::InitializeActivationDescriptor(fActivationDesc,fF);
}

//______________________________________________________________________________
template <typename Architecture_t>
TDenseLayer<Architecture_t>::TDenseLayer(TDenseLayer<Architecture_t> *layer) :
   VGeneralLayer<Architecture_t>(layer), 
   fInputActivation( layer->GetInputActivation().GetShape() ), 
   fDropoutProbability(layer->GetDropoutProbability()),
   fF(layer->GetActivationFunction()), fReg(layer->GetRegularization()), fWeightDecay(layer->GetWeightDecay())
{
   fDerivatives = Tensor_t ( this->GetOutput().GetShape() );
   Architecture_t::InitializeActivationDescriptor(fActivationDesc,fF);
}

//______________________________________________________________________________
template <typename Architecture_t>
TDenseLayer<Architecture_t>::TDenseLayer(const TDenseLayer &layer) :
   VGeneralLayer<Architecture_t>(layer), 
   fInputActivation( layer->GetInputActivation()), 
   fDropoutProbability(layer.fDropoutProbability), 
   fF(layer.fF), fReg(layer.fReg), fWeightDecay(layer.fWeightDecay)
{
   fDerivatives = Tensor_t ( this->GetOutput().GetShape() );
   Architecture_t::InitializeActivationDescriptor(fActivationDesc,fF);
}

//______________________________________________________________________________
template <typename Architecture_t>
TDenseLayer<Architecture_t>::~TDenseLayer()
{
   // release activation descriptor
   Architecture_t::ReleaseDescriptor(fActivationDesc);
}




//______________________________________________________________________________
template <typename Architecture_t>
auto TDenseLayer<Architecture_t>::Forward( Tensor_t &input, bool applyDropout) -> void
{
   if (applyDropout && (this->GetDropoutProbability() != 1.0)) {
      // 
      Architecture_t::DropoutForward(input, static_cast<TDescriptors *> (nullptr), 
                                     static_cast<TWorkspace *> (nullptr), 
                                     this->GetDropoutProbability());
   }
   Architecture_t::MultiplyTranspose(this->GetOutput() , input, this->GetWeightsAt(0));
   Architecture_t::AddRowWise(this->GetOutput(), this->GetBiasesAt(0));
  
   //evaluate<Architecture_t>(this->GetOutput(), this->GetActivationFunction());
   Architecture_t::Copy(this->GetInputActivation(),this->GetOutput());
   
   Architecture_t::ActivationFunctionForward(this->GetOutput(), this->GetActivationFunction(), fActivationDesc);
}

//______________________________________________________________________________
template <typename Architecture_t>
auto TDenseLayer<Architecture_t>::Backward(Tensor_t &gradients_backward, const Tensor_t &activations_backward) -> void
///                                           std::vector<Matrix_t> & /*inp1*/, std::vector<Matrix_t> &
////                                           /*inp2*/) -> void
{

   if (this->GetDropoutProbability() != 1.0) {
      Architecture_t::DropoutBackward(this->GetActivationGradients(), 
      static_cast<TDescriptors *> (nullptr), 
      static_cast<TWorkspace *> (nullptr));
   }

   Architecture_t::ActivationFunctionBackward(fDerivatives, this->GetOutput(), 
                                              this->GetActivationGradients(), this->GetInputActivation(),
                                              this->GetActivationFunction(), fActivationDesc);

   Architecture_t::Backward(gradients_backward, this->GetWeightGradientsAt(0), this->GetBiasGradientsAt(0),
                            fDerivatives, this->GetActivationGradients(), this->GetWeightsAt(0),
                            activations_backward);

   addRegularizationGradients<Architecture_t>(this->GetWeightGradientsAt(0), this->GetWeightsAt(0),
                                              this->GetWeightDecay(), this->GetRegularization());
}

//______________________________________________________________________________
template <typename Architecture_t>
void TDenseLayer<Architecture_t>::Print() const
{
   std::cout << " DENSE Layer: \t";
   std::cout << " ( Input =" << std::setw(6) << this->GetWeightsAt(0).GetNcols();  // input size 
   std::cout << " , Width =" << std::setw(6) << this->GetWeightsAt(0).GetNrows() << " ) ";  // layer width
  
   std::cout << "\tOutput = ( " << std::setw(2) << this->GetOutput().GetFirstSize() << " ," << std::setw(6) << this->GetOutput().GetShape()[0] << " ," << std::setw(6) << this->GetOutput().GetShape()[1] << " ) ";
   
   std::vector<std::string> activationNames = { "Identity","Relu","Sigmoid","Tanh","SymmRelu","SoftSign","Gauss" };
   std::cout << "\t Activation Function = ";
   std::cout << activationNames[ static_cast<int>(fF) ];
   if (fDropoutProbability != 1.) std::cout << "\t Dropout prob. = " << fDropoutProbability;
   std::cout << std::endl;
}

//______________________________________________________________________________

template <typename Architecture_t>
void TDenseLayer<Architecture_t>::AddWeightsXMLTo(void *parent)
{
  // write layer width activation function + weigbht and bias matrices

   auto layerxml = gTools().xmlengine().NewChild(parent, 0, "DenseLayer");

   gTools().xmlengine().NewAttr(layerxml, 0, "Width", gTools().StringFromInt(this->GetWidth()));

   int activationFunction = static_cast<int>(this -> GetActivationFunction());
   gTools().xmlengine().NewAttr(layerxml, 0, "ActivationFunction",
                                TString::Itoa(activationFunction, 10));
   // write weights and bias matrix 
   this->WriteMatrixToXML(layerxml, "Weights", this -> GetWeightsAt(0));
   this->WriteMatrixToXML(layerxml, "Biases",  this -> GetBiasesAt(0));
}

//______________________________________________________________________________
template <typename Architecture_t>
void TDenseLayer<Architecture_t>::ReadWeightsFromXML(void *parent)
{
   // Read layer weights and biases from XML
   this->ReadMatrixXML(parent,"Weights", this -> GetWeightsAt(0));
   this->ReadMatrixXML(parent,"Biases", this -> GetBiasesAt(0));
   
}


} // namespace DNN
} // namespace TMVA

#endif
