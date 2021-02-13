// @(#)root/tmva/tmva/dnn:$Id$
// Author: Simon Pfreundschuh 20/06/16

/*************************************************************************
 * Copyright (C) 2016, Simon Pfreundschuh                                *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

///////////////////////////////////////////////////////////////////////
// Declaration of the TReference architecture, which provides a      //
// reference implementation of the low-level interface for the DNN   //
// implementation based on ROOT's TMatrixT matrix type.              //
///////////////////////////////////////////////////////////////////////

#ifndef TMVA_DNN_ARCHITECTURES_REFERENCE
#define TMVA_DNN_ARCHITECTURES_REFERENCE

#include "TMatrix.h"
//#include "TMVA/RTensor.hxx"
#include "TMVA/DNN/Functions.h"
#include "TMVA/DNN/CNN/ConvLayer.h"
#include "TMVA/DNN/Architectures/Reference/DataLoader.h"
#include "TMVA/DNN/Architectures/Reference/TensorDataLoader.h"
#include <vector>


class TRandom;

namespace TMVA
{
namespace DNN
{
// struct TDescriptors {
// };
// struct TWorkspace {
// };

/*! The reference architecture class.
*
* Class template that contains the reference implementation of the low-level
* interface for the DNN implementation. The reference implementation uses the
* TMatrixT class template to represent matrices.
*
* \tparam AReal The floating point type used to represent scalars.
*/


template<typename AReal>
class TReference
{
private:
   static TRandom * fgRandomGen;
public:
   using Scalar_t     = AReal;
   using Matrix_t     = TMatrixT<AReal>;
   using Tensor_t     = TMatrixT<AReal>;
   //using Tensor_t     = TMVA::Experimental::RTensor<AReal>;

   //____________________________________________________________________________
   //
   // Propagation
   //____________________________________________________________________________

   /** @name Forward Propagation
    * Low-level functions required for the forward propagation of activations
    * through the network.
    */
   ///@{
   /** Matrix-multiply \p input with the transpose of \pweights and
    *  write the results into \p output. */

   static void MultiplyTranspose(TMatrixT<Scalar_t> &output,
                                 const TMatrixT<Scalar_t> &input,
                                 const TMatrixT<Scalar_t> &weights);


   /** Add the vectors biases row-wise to the matrix output */
   static void AddRowWise(TMatrixT<Scalar_t> &output,
                          const TMatrixT<Scalar_t> &biases);
   ///@}

   /** @name Backward Propagation
    * Low-level functions required for the forward propagation of activations
    * through the network.
    */
   ///@{
   /** Perform the complete backward propagation step. If the provided
    *  \p activationGradientsBackward matrix is not empty, compute the
    *  gradients of the objective function with respect to the activations
    *  of the previous layer (backward direction).
    *  Also compute the weight and the bias gradients. Modifies the values
    *  in \p df and thus produces only a valid result, if it is applied the
    *  first time after the corresponding forward propagation has been per-
    *  formed. */
   static void Backward(TMatrixT<Scalar_t> & activationGradientsBackward,
                        TMatrixT<Scalar_t> & weightGradients,
                        TMatrixT<Scalar_t> & biasGradients,
                        TMatrixT<Scalar_t> & df,
                        const TMatrixT<Scalar_t> & activationGradients,
                        const TMatrixT<Scalar_t> & weights,
                        const TMatrixT<Scalar_t> & activationBackward);
   /** Backpropagation step for a Recurrent Neural Network */
   static Matrix_t & RecurrentLayerBackward(TMatrixT<Scalar_t> & state_gradients_backward, // BxH
                                            TMatrixT<Scalar_t> & input_weight_gradients,
                                            TMatrixT<Scalar_t> & state_weight_gradients,
                                            TMatrixT<Scalar_t> & bias_gradients,
                                            TMatrixT<Scalar_t> & df, //DxH
                                            const TMatrixT<Scalar_t> & state, // BxH
                                            const TMatrixT<Scalar_t> & weights_input, // HxD
                                            const TMatrixT<Scalar_t> & weights_state, // HxH
                                            const TMatrixT<Scalar_t> & input,  // BxD
                                            TMatrixT<Scalar_t> & input_gradient);



   /** Backward pass for LSTM Network */
   static Matrix_t & LSTMLayerBackward(TMatrixT<Scalar_t> & state_gradients_backward,
			                              TMatrixT<Scalar_t> & cell_gradients_backward,
			                              TMatrixT<Scalar_t> & input_weight_gradients,
                                       TMatrixT<Scalar_t> & forget_weight_gradients,
                                       TMatrixT<Scalar_t> & candidate_weight_gradients,
                                       TMatrixT<Scalar_t> & output_weight_gradients,
                                       TMatrixT<Scalar_t> & input_state_weight_gradients,
                                       TMatrixT<Scalar_t> & forget_state_weight_gradients,
                                       TMatrixT<Scalar_t> & candidate_state_weight_gradients,
                                       TMatrixT<Scalar_t> & output_state_weight_gradients,
                                       TMatrixT<Scalar_t> & input_bias_gradients,
                                       TMatrixT<Scalar_t> & forget_bias_gradients,
                                       TMatrixT<Scalar_t> & candidate_bias_gradients,
                                       TMatrixT<Scalar_t> & output_bias_gradients,
                                       TMatrixT<Scalar_t> & di,
                                       TMatrixT<Scalar_t> & df,
                                       TMatrixT<Scalar_t> & dc,
                                       TMatrixT<Scalar_t> & dout,
                                       const TMatrixT<Scalar_t> & precStateActivations,
                                       const TMatrixT<Scalar_t> & precCellActivations,
                                       const TMatrixT<Scalar_t> & fInput,
                                       const TMatrixT<Scalar_t> & fForget,
                                       const TMatrixT<Scalar_t> & fCandidate,
                                       const TMatrixT<Scalar_t> & fOutput,
                                       const TMatrixT<Scalar_t> & weights_input,
                                       const TMatrixT<Scalar_t> & weights_forget,
                                       const TMatrixT<Scalar_t> & weights_candidate,
                                       const TMatrixT<Scalar_t> & weights_output,
                                       const TMatrixT<Scalar_t> & weights_input_state,
                                       const TMatrixT<Scalar_t> & weights_forget_state,
                                       const TMatrixT<Scalar_t> & weights_candidate_state,
                                       const TMatrixT<Scalar_t> & weights_output_state,
                                       const TMatrixT<Scalar_t> & input,
                                       TMatrixT<Scalar_t> & input_gradient,
                                       TMatrixT<Scalar_t> & cell_gradient,
                                       TMatrixT<Scalar_t> & cell_tanh);


 /** Backward pass for GRU Network */
   static Matrix_t & GRULayerBackward(TMatrixT<Scalar_t> & state_gradients_backward,
                                      TMatrixT<Scalar_t> & reset_weight_gradients,
                                      TMatrixT<Scalar_t> & update_weight_gradients,
                                      TMatrixT<Scalar_t> & candidate_weight_gradients,
                                      TMatrixT<Scalar_t> & reset_state_weight_gradients,
                                      TMatrixT<Scalar_t> & update_state_weight_gradients,
                                      TMatrixT<Scalar_t> & candidate_state_weight_gradients,
                                      TMatrixT<Scalar_t> & reset_bias_gradients,
                                      TMatrixT<Scalar_t> & update_bias_gradients,
                                      TMatrixT<Scalar_t> & candidate_bias_gradients,
                                      TMatrixT<Scalar_t> & dr,
                                      TMatrixT<Scalar_t> & du,
                                      TMatrixT<Scalar_t> & dc,
                                      const TMatrixT<Scalar_t> & precStateActivations,
                                      const TMatrixT<Scalar_t> & fReset,
                                      const TMatrixT<Scalar_t> & fUpdate,
                                      const TMatrixT<Scalar_t> & fCandidate,
                                      const TMatrixT<Scalar_t> & weights_reset,
                                      const TMatrixT<Scalar_t> & weights_update,
                                      const TMatrixT<Scalar_t> & weights_candidate,
                                      const TMatrixT<Scalar_t> & weights_reset_state,
                                      const TMatrixT<Scalar_t> & weights_update_state,
                                      const TMatrixT<Scalar_t> & weights_candidate_state,
                                      const TMatrixT<Scalar_t> & input,
                                      TMatrixT<Scalar_t> & input_gradient);

   /** Adds a the elements in matrix B scaled by c to the elements in
    *  the matrix A. This is required for the weight update in the gradient
    *  descent step.*/
   static void ScaleAdd(TMatrixT<Scalar_t> & A,
                        const TMatrixT<Scalar_t> & B,
                        Scalar_t beta = 1.0);

   static void Copy(TMatrixT<Scalar_t> & A,
                    const TMatrixT<Scalar_t> & B);

   // copy from another type of matrix
   template<typename AMatrix_t>
   static void CopyDiffArch(TMatrixT<Scalar_t> & A, const AMatrix_t & B);


   /** Above functions extended to vectors */
   static void ScaleAdd(std::vector<TMatrixT<Scalar_t>> & A,
                        const std::vector<TMatrixT<Scalar_t>> & B,
                        Scalar_t beta = 1.0);

   static void Copy(std::vector<TMatrixT<Scalar_t>> & A, const std::vector<TMatrixT<Scalar_t>> & B);

   // copy from another architecture
   template<typename AMatrix_t>
   static void CopyDiffArch(std::vector<TMatrixT<Scalar_t> > & A, const std::vector<AMatrix_t> & B);


   ///@}

   //____________________________________________________________________________
   //
   // Activation Functions
   //____________________________________________________________________________

   /** @name Activation Functions
    * For each activation function, the low-level interface contains two routines.
    * One that applies the acitvation function to a matrix and one that evaluate
    * the derivatives of the activation function at the elements of a given matrix
    * and writes the results into the result matrix.
    */
   ///@{
   static void Identity(TMatrixT<AReal> & B);
   static void IdentityDerivative(TMatrixT<AReal> & B,
                                  const TMatrixT<AReal> & A);

   static void Relu(TMatrixT<AReal> & B);
   static void ReluDerivative(TMatrixT<AReal> & B,
                              const TMatrixT<AReal> & A);

   static void Sigmoid(TMatrixT<AReal> & B);
   static void SigmoidDerivative(TMatrixT<AReal> & B,
                                 const TMatrixT<AReal> & A);

   static void Tanh(TMatrixT<AReal> & B);
   static void TanhDerivative(TMatrixT<AReal> & B,
                              const TMatrixT<AReal> & A);

   static void FastTanh(Tensor_t &B) { return Tanh(B); }
   static void FastTanhDerivative(Tensor_t &B, const Tensor_t &A) { return TanhDerivative(B, A); }

   static void SymmetricRelu(TMatrixT<AReal> & B);
   static void SymmetricReluDerivative(TMatrixT<AReal> & B,
                                       const TMatrixT<AReal> & A);

   static void SoftSign(TMatrixT<AReal> & B);
   static void SoftSignDerivative(TMatrixT<AReal> & B,
                                  const TMatrixT<AReal> & A);

   static void Gauss(TMatrixT<AReal> & B);
   static void GaussDerivative(TMatrixT<AReal> & B,
                               const TMatrixT<AReal> & A);


   ///@}

   //____________________________________________________________________________
   //
   // Loss Functions
   //____________________________________________________________________________

   /** @name Loss Functions
    * Loss functions compute a scalar value given the \p output of the network
    * for a given training input and the expected network prediction \p Y that
    * quantifies the quality of the prediction. For each function also a routing
    * that computes the gradients (suffixed by Gradients) must be provided for
    * the starting of the backpropagation algorithm.
    */
   ///@{

   static AReal MeanSquaredError(const TMatrixT<AReal> &Y, const TMatrixT<AReal> &output,
                                 const TMatrixT<AReal> &weights);
   static void MeanSquaredErrorGradients(TMatrixT<AReal> &dY, const TMatrixT<AReal> &Y, const TMatrixT<AReal> &output,
                                         const TMatrixT<AReal> &weights);

   /** Sigmoid transformation is implicitly applied, thus \p output should
    *  hold the linear activations of the last layer in the net. */
   static AReal CrossEntropy(const TMatrixT<AReal> &Y, const TMatrixT<AReal> &output, const TMatrixT<AReal> &weights);

   static void CrossEntropyGradients(TMatrixT<AReal> &dY, const TMatrixT<AReal> &Y, const TMatrixT<AReal> &output,
                                     const TMatrixT<AReal> &weights);

   /** Softmax transformation is implicitly applied, thus \p output should
    *  hold the linear activations of the last layer in the net. */
   static AReal SoftmaxCrossEntropy(const TMatrixT<AReal> &Y, const TMatrixT<AReal> &output,
                                    const TMatrixT<AReal> &weights);
   static void SoftmaxCrossEntropyGradients(TMatrixT<AReal> &dY, const TMatrixT<AReal> &Y,
                                            const TMatrixT<AReal> &output, const TMatrixT<AReal> &weights);
   ///@}

   //____________________________________________________________________________
   //
   // Output Functions
   //____________________________________________________________________________

   /** @name Output Functions
    * Output functions transform the activations \p output of the
    * output layer in the network to a valid prediction \p YHat for
    * the desired usage of the network, e.g.  the identity function
    * for regression or the sigmoid transformation for two-class
    * classification.
    */
   ///@{
   static void Sigmoid(TMatrixT<AReal> &YHat,
                       const TMatrixT<AReal> & );
   static void Softmax(TMatrixT<AReal> &YHat,
                       const TMatrixT<AReal> & );
   ///@}

   //____________________________________________________________________________
   //
   // Regularization
   //____________________________________________________________________________

   /** @name Regularization
    * For each regularization type two functions are required, one named
    * <tt><Type>Regularization</tt> that evaluates the corresponding
    * regularization functional for a given weight matrix and the
    * <tt>Add<Type>RegularizationGradients</tt>, that adds the regularization
    * component in the gradients to the provided matrix.
    */
   ///@{

   static AReal L1Regularization(const TMatrixT<AReal> & W);
   static void AddL1RegularizationGradients(TMatrixT<AReal> & A,
                                            const TMatrixT<AReal> & W,
                                            AReal weightDecay);

   static AReal L2Regularization(const TMatrixT<AReal> & W);
   static void AddL2RegularizationGradients(TMatrixT<AReal> & A,
                                            const TMatrixT<AReal> & W,
                                            AReal weightDecay);
   ///@}

   //____________________________________________________________________________
   //
   // Initialization
   //____________________________________________________________________________

   /** @name Initialization
    * For each initialization method, one function in the low-level interface
    * is provided. The naming scheme is <p>Initialize<Type></p> for a given
    * initialization method Type.
    */
   ///@{

   static void InitializeGauss(TMatrixT<AReal> & A);

   static void InitializeUniform(TMatrixT<AReal> & A);

   static void InitializeIdentity(TMatrixT<AReal> & A);

   static void InitializeZero(TMatrixT<AReal> & A);

   static void InitializeGlorotUniform(TMatrixT<AReal> & A);

   static void InitializeGlorotNormal(TMatrixT<AReal> & A);

   // return static instance of random generator used for initialization
   // if generator does not exist it is created the first time with a random seed (e.g. seed = 0)
   static TRandom & GetRandomGenerator();
   // set random seed for the static geenrator
   // if the static geneerator does not exists it is created
   static void SetRandomSeed(size_t seed);


   ///@}

   //____________________________________________________________________________
   //
   // Dropout
   //____________________________________________________________________________

   /** @name Dropout
    */
   ///@{

   /** Apply dropout with activation probability \p p to the given
    *  matrix \p A and scale the result by reciprocal of \p p. */
   //static void Dropout(TMatrixT<AReal> & A, AReal dropoutProbability);
   static void DropoutForward(Tensor_t &A, TDescriptors *descriptors, TWorkspace *workspace, Scalar_t p);
   static void DropoutForward(Matrix_t &A, Scalar_t p)
   {
      Tensor_t & tA = A;  // Tensor and matrix are same types
      DropoutForward(tA, static_cast<TDescriptors *>(nullptr), static_cast<TWorkspace *>(nullptr), p);
   }

   ///@}


   //____________________________________________________________________________
   //
   //  Convolutional Layer Propagation
   //____________________________________________________________________________

   /** @name Forward Propagation in Convolutional Layer
    */
   ///@{

   /** Transform the matrix \p B in local view format, suitable for
    *  convolution, and store it in matrix \p A. */
   static void Im2col(TMatrixT<AReal> &A,
                      const TMatrixT<AReal> &B,
                      size_t imgHeight,
                      size_t imgWidth,
                      size_t fltHeight,
                      size_t fltWidth,
                      size_t strideRows,
                      size_t strideCols,
                      size_t zeroPaddingHeight,
                      size_t zeroPaddingWidth);

   static void Im2colIndices(std::vector<int> &, const TMatrixT<AReal> &, size_t, size_t, size_t, size_t ,
                      size_t , size_t , size_t , size_t ,size_t ) {
      Fatal("Im2ColIndices","This function is not implemented for ref architectures");
   }
   static void Im2colFast(TMatrixT<AReal> &, const TMatrixT<AReal> &, const std::vector<int> & ) {
       Fatal("Im2ColFast","This function is not implemented for ref architectures");
   }

   /** Rotates the matrix \p B, which is representing a weights,
    *  and stores them in the matrix \p A. */
   static void RotateWeights(TMatrixT<AReal> &A, const TMatrixT<AReal> &B, size_t filterDepth, size_t filterHeight,
                             size_t filterWidth, size_t numFilters);

   /** Add the biases in the Convolutional Layer.  */
   static void AddConvBiases(TMatrixT<AReal> &output, const TMatrixT<AReal> &biases);
   ///@}

   /** Dummy placeholder - preparation is currently only required for the CUDA architecture. */
   static void PrepareInternals(std::vector<TMatrixT<AReal>> &) {}

   /** Forward propagation in the Convolutional layer */
   static void ConvLayerForward(std::vector<TMatrixT<AReal>> & /*output*/,
                                std::vector<TMatrixT<AReal>> & /*derivatives*/,
                                const std::vector<TMatrixT<AReal>> & /*input*/,
                                const TMatrixT<AReal> & /*weights*/, const TMatrixT<AReal> & /*biases*/,
                                const DNN::CNN::TConvParams & /*params*/, EActivationFunction /*activFunc*/,
                                std::vector<TMatrixT<AReal>> & /*inputPrime*/) {
      Fatal("ConvLayerForward","This function is not implemented for ref architectures");
   }


   /** @name Backward Propagation in Convolutional Layer
    */
   ///@{

   /** Perform the complete backward propagation step in a Convolutional Layer.
    *  If the provided \p activationGradientsBackward matrix is not empty, compute the
    *  gradients of the objective function with respect to the activations
    *  of the previous layer (backward direction).
    *  Also compute the weight and the bias gradients. Modifies the values
    *  in \p df and thus produces only a valid result, if it is applied the
    *  first time after the corresponding forward propagation has been per-
    *  formed. */
   static void ConvLayerBackward(std::vector<TMatrixT<AReal>> &,
                                 TMatrixT<AReal> &, TMatrixT<AReal> &,
                                 std::vector<TMatrixT<AReal>> &,
                                 const std::vector<TMatrixT<AReal>> &,
                                 const TMatrixT<AReal> &, const std::vector<TMatrixT<AReal>> &,
                                 size_t , size_t , size_t , size_t , size_t,
                                 size_t , size_t , size_t , size_t , size_t) {
      Fatal("ConvLayerBackward","This function is not implemented for ref architectures");

   }

#ifdef HAVE_CNN_REFERENCE
   /** Utility function for calculating the activation gradients of the layer
    *  before the convolutional layer. */
   static void CalculateConvActivationGradients(std::vector<TMatrixT<AReal>> &activationGradientsBackward,
                                                const std::vector<TMatrixT<AReal>> &df, const TMatrixT<AReal> &weights,
                                                size_t batchSize, size_t inputHeight, size_t inputWidth, size_t depth,
                                                size_t height, size_t width, size_t filterDepth, size_t filterHeight,
                                                size_t filterWidth);

   /** Utility function for calculating the weight gradients of the convolutional
    *  layer. */
   static void CalculateConvWeightGradients(TMatrixT<AReal> &weightGradients, const std::vector<TMatrixT<AReal>> &df,
                                            const std::vector<TMatrixT<AReal>> &activationBackward, size_t batchSize,
                                            size_t inputHeight, size_t inputWidth, size_t depth, size_t height,
                                            size_t width, size_t filterDepth, size_t filterHeight, size_t filterWidth,
                                            size_t nLocalViews);

   /** Utility function for calculating the bias gradients of the convolutional
    *  layer. */
   static void CalculateConvBiasGradients(TMatrixT<AReal> &biasGradients, const std::vector<TMatrixT<AReal>> &df,
                                          size_t batchSize, size_t depth, size_t nLocalViews);
   ///@}

#endif

   //____________________________________________________________________________
   //
   //  Max Pooling Layer Propagation
   //____________________________________________________________________________
   /** @name Forward Propagation in Max Pooling Layer
    */
   ///@{

  /** Downsample the matrix \p C to the matrix \p A, using max
    *  operation, such that the winning indices are stored in matrix
    *  \p B. */
   static void Downsample(TMatrixT<AReal> &A, TMatrixT<AReal> &B, const TMatrixT<AReal> &C, size_t imgHeight,
                          size_t imgWidth, size_t fltHeight, size_t fltWidth, size_t strideRows, size_t strideCols);

   ///@}

   /** @name Backward Propagation in Max Pooling Layer
    */
   ///@{

   /** Perform the complete backward propagation step in a Max Pooling Layer. Based on the
    *  winning idices stored in the index matrix, it just forwards the actiovation
    *  gradients to the previous layer. */
   static void MaxPoolLayerBackward(TMatrixT<AReal> &activationGradientsBackward,
                                    const TMatrixT<AReal> &activationGradients,
                                    const TMatrixT<AReal> &indexMatrix,
                                    size_t imgHeight,
                                    size_t imgWidth,
                                    size_t fltHeight,
                                    size_t fltWidth,
                                    size_t strideRows,
                                    size_t strideCol,
                                    size_t nLocalViews);
   ///@}
   //____________________________________________________________________________
   //
   //  Reshape Layer Propagation
   //____________________________________________________________________________
   /** @name Forward and Backward Propagation in Reshape Layer
    */
   ///@{

   /** Transform the matrix \p B to a matrix with different dimensions \p A */
   static void Reshape(TMatrixT<AReal> &A, const TMatrixT<AReal> &B);

   /** Flattens the tensor \p B, such that each matrix, is stretched in one row, resulting with a matrix \p A. */
   static void Flatten(TMatrixT<AReal> &A, const std::vector<TMatrixT<AReal>> &B, size_t size, size_t nRows,
                       size_t nCols);

   /** Transforms each row of \p B to a matrix and stores it in the tensor \p B. */
   static void Deflatten(std::vector<TMatrixT<AReal>> &A, const TMatrixT<Scalar_t> &B, size_t index, size_t nRows,
                         size_t nCols);
   /** Rearrage data accoring to time fill B x T x D out with T x B x D matrix in*/
   static void Rearrange(std::vector<TMatrixT<AReal>> &out, const std::vector<TMatrixT<AReal>> &in);

   ///@}

   //____________________________________________________________________________
   //
   // Additional Arithmetic Functions
   //____________________________________________________________________________

   /** Sum columns of (m x n) matrixx \p A and write the results into the first
    * m elements in \p A.
    */
   static void SumColumns(TMatrixT<AReal> &B, const TMatrixT<AReal> &A);

   /** In-place Hadamard (element-wise) product of matrices \p A and \p B
    *  with the result being written into \p A.
    */
   static void Hadamard(TMatrixT<AReal> &A, const TMatrixT<AReal> &B);

   /** Add the constant \p beta to all the elements of matrix \p A and write the
    * result into \p A.
    */
   static void ConstAdd(TMatrixT<AReal> &A, AReal beta);

   /** Multiply the constant \p beta to all the elements of matrix \p A and write the
    * result into \p A.
    */
   static void ConstMult(TMatrixT<AReal> &A, AReal beta);

   /** Reciprocal each element of the matrix \p A and write the result into
    * \p A
    */
   static void ReciprocalElementWise(TMatrixT<AReal> &A);

   /** Square each element of the matrix \p A and write the result into
    * \p A
    */
   static void SquareElementWise(TMatrixT<AReal> &A);

   /** Square root each element of the matrix \p A and write the result into
    * \p A
    */
   static void SqrtElementWise(TMatrixT<AReal> &A);

   // optimizer update functions

   /// Update functions for ADAM optimizer
   static void AdamUpdate(TMatrixT<AReal> & A, const TMatrixT<AReal> & M, const TMatrixT<AReal> & V, AReal alpha, AReal eps);
   static void AdamUpdateFirstMom(TMatrixT<AReal> & A, const TMatrixT<AReal> & B, AReal beta);
   static void AdamUpdateSecondMom(TMatrixT<AReal> & A, const TMatrixT<AReal> & B, AReal beta);



   //____________________________________________________________________________
   //
   // AutoEncoder Propagation
   //____________________________________________________________________________

   // Add Biases to the output
   static void AddBiases(TMatrixT<AReal> &A,
                         const TMatrixT<AReal> &biases);

   // Updating parameters after every backward pass. Weights and biases are
   // updated.
   static void
   UpdateParams(TMatrixT<AReal> &x, TMatrixT<AReal> &tildeX, TMatrixT<AReal> &y,
                TMatrixT<AReal> &z, TMatrixT<AReal> &fVBiases,
                TMatrixT<AReal> &fHBiases, TMatrixT<AReal> &fWeights,
                TMatrixT<AReal> &VBiasError, TMatrixT<AReal> &HBiasError,
                AReal learningRate, size_t fBatchSize);

   // Softmax functions redifined
   static void SoftmaxAE(TMatrixT<AReal> & A);


   // Corrupt the input values randomly on corruption Level.
   //Basically inputs are masked currently.
   static void CorruptInput(TMatrixT<AReal> & input,
                            TMatrixT<AReal> & corruptedInput,
                            AReal corruptionLevel);

   //Encodes the input Values in the compressed form.
   static void EncodeInput(TMatrixT<AReal> &input,
                           TMatrixT<AReal> &compressedInput,
                           TMatrixT<AReal> &Weights);

   // reconstructs the input. The reconstructed Input has same dimensions as that
   // of the input.
   static void ReconstructInput(TMatrixT<AReal> & compressedInput,
                                TMatrixT<AReal> & reconstructedInput,
                                TMatrixT<AReal> &fWeights);


   static void ForwardLogReg(TMatrixT<AReal> &input,
                             TMatrixT<AReal> &p,
                             TMatrixT<AReal> &fWeights);

   static void UpdateParamsLogReg(TMatrixT<AReal> &input,
                                  TMatrixT<AReal> &output,
                                  TMatrixT<AReal> &difference,
                                  TMatrixT<AReal> &p,
                                  TMatrixT<AReal> &fWeights,
                                  TMatrixT<AReal> &fBiases,
                                  AReal learningRate,
                                  size_t fBatchSize);

};


// implement the templated member functions
template <typename AReal>
template <typename AMatrix_t>
void TReference<AReal>::CopyDiffArch(TMatrixT<AReal> &A, const AMatrix_t &B)
{
   TMatrixT<AReal> tmp = B;
   A = tmp;
}

template <typename AReal>
template <typename AMatrix_t>
void TReference<AReal>::CopyDiffArch(std::vector<TMatrixT<AReal>> &A, const std::vector<AMatrix_t> &B)
{
   for (size_t i = 0; i < A.size(); ++i) {
      CopyDiffArch(A[i], B[i]);
   }
}



} // namespace DNN
} // namespace TMVA

#endif
