// @(#)root/tmva/tmva/dnn:$Id$ // Author: Simon Pfreundschuh 10/07/16

/*************************************************************************
 * Copyright (C) 2016, Simon Pfreundschuh                                *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/////////////////////////////////////////////////////////////////////
// Implementation of the functions required for the forward and    //
// backward propagation of activations through a neural network in //
// the reference implementation.                                   //
/////////////////////////////////////////////////////////////////////

#include "TMVA/DNN/Architectures/Reference.h"

namespace TMVA {
namespace DNN {

template <typename AReal>
void TReference<AReal>::MultiplyTranspose(TMatrixT<AReal> &output, const TMatrixT<AReal> &input,
                                          const TMatrixT<AReal> &weights)
{
   output.MultT(input, weights);
}

template <typename AReal>
void TReference<AReal>::AddRowWise(TMatrixT<AReal> &output, const TMatrixT<AReal> &biases)
{
   for (size_t i = 0; i < (size_t)output.GetNrows(); i++) {
      for (size_t j = 0; j < (size_t)output.GetNcols(); j++) {
         output(i, j) += biases(j, 0);
      }
   }
}

template <typename AReal>
void TReference<AReal>::Backward(TMatrixT<AReal> &activation_gradients_backward, TMatrixT<AReal> &weight_gradients,
                                 TMatrixT<AReal> &bias_gradients, TMatrixT<AReal> &df,
                                 const TMatrixT<AReal> &activation_gradients, const TMatrixT<AReal> &weights,
                                 const TMatrixT<AReal> &activations_backward)
{

   // Compute element-wise product.
   for (size_t i = 0; i < (size_t)df.GetNrows(); i++) {
      for (size_t j = 0; j < (size_t)df.GetNcols(); j++) {
         df(i, j) *= activation_gradients(i, j);
      }
   }

   // Activation gradients.
   if (activation_gradients_backward.GetNoElements() > 0) {
      activation_gradients_backward.Mult(df, weights);
   }

   // Weights gradients.
   if (weight_gradients.GetNoElements() > 0) {
      weight_gradients.TMult(df, activations_backward);
   }

   // Bias gradients.
   if (bias_gradients.GetNoElements() > 0) {
      for (size_t j = 0; j < (size_t)df.GetNcols(); j++) {
         AReal sum = 0.0;
         for (size_t i = 0; i < (size_t)df.GetNrows(); i++) {
            sum += df(i, j);
         }
         bias_gradients(j, 0) = sum;
      }
   }
}

template <typename AReal>
void TReference<AReal>::ScaleAdd(TMatrixT<AReal> &A, const TMatrixT<AReal> &B, AReal beta)
{
   for (size_t i = 0; i < (size_t)A.GetNrows(); i++) {
      for (size_t j = 0; j < (size_t)A.GetNcols(); j++) {
         A(i, j) += beta * B(i, j);
      }
   }
}

template <typename AReal>
void TReference<AReal>::Copy(TMatrixT<AReal> &A, const TMatrixT<AReal> &B)
{
   A = B;
}

template <typename AReal>
void TReference<AReal>::ScaleAdd(std::vector<TMatrixT<AReal>> &A, const std::vector<TMatrixT<AReal>> &B, AReal beta)
{
   for (size_t i = 0; i < A.size(); ++i) {
      ScaleAdd(A[i], B[i], beta);
   }
}

template <typename AReal>
void TReference<AReal>::Copy(std::vector<TMatrixT<AReal>> &A, const std::vector<TMatrixT<AReal>> &B)
{
   for (size_t i = 0; i < A.size(); ++i) {
      Copy(A[i], B[i]);
   }
}

//______________________________________________________________________________
template <typename AReal>
void TReference<AReal>::Im2col(TMatrixT<AReal> &A, const TMatrixT<AReal> &B, size_t imgHeight, size_t imgWidth,
                               size_t fltHeight, size_t fltWidth, size_t strideRows, size_t strideCols,
                               size_t zeroPaddingHeight, size_t zeroPaddingWidth)
{
   // image boudaries
   int imgHeightBound = imgHeight + zeroPaddingHeight - (fltHeight - 1) / 2 - 1;
   int imgWidthBound = imgWidth + zeroPaddingWidth - (fltWidth - 1) / 2 - 1;
   size_t currLocalView = 0;

   // convolution centers
   for (int i = -1*int(zeroPaddingHeight) + fltHeight / 2; i <= imgHeightBound; i += strideRows) {
      for (int j = -1*int(zeroPaddingWidth) + fltWidth / 2; j <= imgWidthBound; j += strideCols) {
         size_t currLocalViewPixel = 0;

         // within the local view
         for (int m = 0; m < B.GetNrows(); m++) {
            for (Int_t k = i - Int_t(fltHeight) / 2; k <= i + (Int_t(fltHeight) - 1) / 2; k++) {
               for (Int_t l = j - Int_t(fltWidth) / 2; l <= j + (Int_t(fltWidth) - 1) / 2; l++) {

                  // Check the boundaries
                  if (k < 0 || k >= Int_t(imgHeight) || l < 0 || l >= Int_t(imgWidth))
                     A(currLocalView, currLocalViewPixel++) = 0;
                  else
                     A(currLocalView, currLocalViewPixel++) = B(m, k * imgWidth + l);
               }
            }
         }

         currLocalView++;
      }
   }
}

//______________________________________________________________________________
template <typename AReal>
void TReference<AReal>::RotateWeights(TMatrixT<AReal> &A, const TMatrixT<AReal> &B, size_t filterDepth,
                                      size_t filterHeight, size_t filterWidth, size_t numFilters)
{
   size_t jump = filterHeight * filterWidth;
   for (size_t j = 0; j < filterDepth; j++) {
      for (size_t k = 0; k < numFilters; k++) {
         for (size_t i = 0; i < jump; i++) {
            A(j, k * jump + i) = B(k, ((j + 1) * jump - 1) - i);
         }
      }
   }
}

//______________________________________________________________________________
template <typename AReal>
void TReference<AReal>::AddConvBiases(TMatrixT<AReal> &output, const TMatrixT<AReal> &biases)
{
   for (size_t i = 0; i < (size_t)output.GetNrows(); i++) {
      for (size_t j = 0; j < (size_t)output.GetNcols(); j++) {
         output(i, j) += biases(i, 0);
      }
   }
}

#ifdef HAVE_CNN_REFERENCE
//______________________________________________________________________________
template <typename AReal>
void TReference<AReal>::ConvLayerBackward(std::vector<TMatrixT<AReal>> &activation_gradients_backward,
                                          TMatrixT<AReal> &weight_gradients, TMatrixT<AReal> &bias_gradients,
                                          std::vector<TMatrixT<AReal>> &df,
                                          const std::vector<TMatrixT<AReal>> &activation_gradients,
                                          const TMatrixT<AReal> &weights,
                                          const std::vector<TMatrixT<AReal>> &activations_backward, size_t batchSize,
                                          size_t inputHeight, size_t inputWidth, size_t depth, size_t height,
                                          size_t width, size_t filterDepth, size_t filterHeight, size_t filterWidth,
                                          size_t nLocalViews)
{

   // Update derivatives
   size_t m, n;
   m = activation_gradients[0].GetNrows();
   n = activation_gradients[0].GetNcols();

   for (size_t i = 0; i < batchSize; i++) {
      for (size_t j = 0; j < (size_t)m; j++) {
         for (size_t k = 0; k < (size_t)n; k++) {
            df[i](j, k) *= activation_gradients[i](j, k);
         }
      }
   }

   // Calculate the activation gradients of the previous layer
   CalculateConvActivationGradients(activation_gradients_backward, df, weights, batchSize, inputHeight, inputWidth,
                                    depth, height, width, filterDepth, filterHeight, filterWidth);

   // Calculate the weight gradients
   CalculateConvWeightGradients(weight_gradients, df, activations_backward, batchSize, inputHeight, inputWidth, depth,
                                height, width, filterDepth, filterHeight, filterWidth, nLocalViews);

   // Calculate the bias gradients
   CalculateConvBiasGradients(bias_gradients, df, batchSize, depth, nLocalViews);
}

//______________________________________________________________________________
template <typename AReal>
void TReference<AReal>::CalculateConvActivationGradients(std::vector<TMatrixT<AReal>> &activation_gradients_backward,
                                                         const std::vector<TMatrixT<AReal>> &df,
                                                         const TMatrixT<AReal> &weights, size_t batchSize,
                                                         size_t inputHeight, size_t inputWidth, size_t depth,
                                                         size_t height, size_t width, size_t filterDepth,
                                                         size_t filterHeight, size_t filterWidth)
{

   if (activation_gradients_backward.size() == 0) return;
   // need to implement
   // Transform the weights
   TMatrixT<AReal> rotWeights(filterDepth, depth * filterHeight * filterWidth);
   RotateWeights(rotWeights, weights, filterDepth, filterHeight, filterWidth, weights.GetNrows());

   // Calculate the zero paddings
   size_t tempZeroPaddingHeight = (size_t)(floor((inputHeight - height + filterHeight - 1) / 2));
   size_t tempZeroPaddingWidth = (size_t)(floor((inputWidth - width + filterWidth - 1) / 2));

   // Calculate the number of local views and the number of pixles in each view
   size_t tempNLocalViews = inputHeight * inputWidth;
   size_t tempNLocalViewPixels = depth * filterHeight * filterWidth;

   size_t tempStrideRows = 1;
   size_t tempStrideCols = 1;

   // An entire convolution follows
   for (size_t i = 0; i < batchSize; i++) {
      TMatrixT<AReal> dfTr(tempNLocalViews, tempNLocalViewPixels);
      Im2col(dfTr, df[i], inputHeight, inputWidth, filterHeight, filterWidth, tempStrideRows, tempStrideCols,
             tempZeroPaddingHeight, tempZeroPaddingWidth);

      activation_gradients_backward[i].MultT(rotWeights, dfTr);
   }

   return ;
}

//______________________________________________________________________________
template <typename AReal>
void TReference<AReal>::CalculateConvWeightGradients(TMatrixT<AReal> &weight_gradients,
                                                     const std::vector<TMatrixT<AReal>> &df,
                                                     const std::vector<TMatrixT<AReal>> &activations_backward,
                                                     size_t batchSize, size_t inputHeight, size_t inputWidth,
                                                     size_t depth, size_t height, size_t width, size_t filterDepth,
                                                     size_t filterHeight, size_t filterWidth, size_t nLocalViews)
{

   // reinitialize the weight gradients to 0
   for (Int_t i = 0; i < weight_gradients.GetNrows(); i++) {
      for (Int_t j = 0; j < weight_gradients.GetNcols(); j++) {
         weight_gradients(i, j) = 0;
      }
   }
   for (size_t i = 0; i < batchSize; i++) {
      // Calculate the zero paddings
      size_t tempZeroPaddingHeight = (filterHeight - height + inputHeight - 1) / 2;
      size_t tempZeroPaddingWidth = (filterWidth - width + inputWidth - 1) / 2;

      size_t tempNLocalViews = filterHeight * filterWidth;
      size_t tempNLocalViewPixels = inputHeight * inputWidth;

      size_t tempStrideRows = 1;
      size_t tempStrideCols = 1;

      for (size_t j = 0; j < depth; j++) {

         // row matrix
         TMatrixT<AReal> rowDelta(1, nLocalViews);
         for (size_t k = 0; k < nLocalViews; k++) {
            rowDelta(0, k) = df[i](j, k);
         }

         // convolution
         TMatrixT<AReal> res(filterDepth, filterHeight * filterWidth);

         TMatrixT<AReal> rowDeltaTr(tempNLocalViews, tempNLocalViewPixels);
         Im2col(rowDeltaTr, rowDelta, height, width, inputHeight, inputWidth, tempStrideRows, tempStrideCols,
                tempZeroPaddingHeight, tempZeroPaddingWidth);

         res.MultT(activations_backward[i], rowDeltaTr);

         for (size_t k = 0; k < filterDepth; k++) {
            for (size_t l = 0; l < filterHeight * filterWidth; l++) {
               weight_gradients(j, k * (filterHeight * filterWidth) + l) += res(k, (tempNLocalViews - 1) - l);
            }
         }
      }
   }
#if 0
   // to remove warning
   (void)weight_gradients;
   (void)df;
   (void)activations_backward;
   (void) batchSize;
   (void) inputHeight;
   (void)inputWidth;
   (void)depth;
   (void)height;
   (void) width;
   (void)filterDepth;
   (void)filterHeight;
   (void)filterWidth;
   (void)nLocalViews;
#endif
}

//______________________________________________________________________________
template <typename AReal>
void TReference<AReal>::CalculateConvBiasGradients(TMatrixT<AReal> &bias_gradients, const std::vector<TMatrixT<AReal>> &df,
                                                   size_t batchSize, size_t depth, size_t nLocalViews)
{
   for (size_t i = 0; i < depth; i++) {
      AReal sum = 0;
      for (size_t j = 0; j < nLocalViews; j++) {
         for (size_t k = 0; k < batchSize; k++) {
            sum += df[k](i, j);
         }
      }
      bias_gradients(i, 0) = sum;
   }
}
#endif

//______________________________________________________________________________
template <typename AReal>
void TReference<AReal>::Downsample(TMatrixT<AReal> &A, TMatrixT<AReal> &B, const TMatrixT<AReal> &C, size_t imgHeight,
                                   size_t imgWidth, size_t fltHeight, size_t fltWidth, size_t strideRows,
                                   size_t strideCols)
{
   // image boudaries
   int imgHeightBound = imgHeight - (fltHeight - 1) / 2 - 1;
   int imgWidthBound = imgWidth - (fltWidth - 1) / 2 - 1;
   size_t currLocalView = 0;

   // centers
   for (int i = fltHeight / 2; i <= imgHeightBound; i += strideRows) {
      for (int j = fltWidth / 2; j <= imgWidthBound; j += strideCols) {
         // within local views
         for (int m = 0; m < C.GetNrows(); m++) {
            AReal value = -std::numeric_limits<AReal>::max();

            for (int k = i - Int_t(fltHeight) / 2; k <= i + (Int_t(fltHeight) - 1) / 2; k++) {
               for (int l = j - Int_t(fltWidth) / 2; l <= j + (Int_t(fltWidth) - 1) / 2; l++) {
                  if (C(m, k * imgWidth + l) > value) {
                     value = C(m, k * imgWidth + l);
                     B(m, currLocalView) = k * imgWidth + l;
                  }
               }
            }
            A(m, currLocalView) = value;
         }
         currLocalView++;
      }
   }
}

//______________________________________________________________________________
template <typename AReal>
void TReference<AReal>::MaxPoolLayerBackward(TMatrixT<AReal> &activationGradientsBackward,
                                             const TMatrixT<AReal> &activationGradients,
                                             const TMatrixT<AReal> &indexMatrix,
                                             size_t /* imgHeight */, size_t /* imgWidth */, size_t /* fltHeight */,
                                             size_t /* fltWidth */, size_t /* strideRows */, size_t /* strideCols */,
                                             size_t nLocalViews)
{
    size_t depth = activationGradientsBackward.GetNrows();

   for (size_t j = 0; j < depth; j++) {
      // initialize to zeros
      for (size_t t = 0; t < (size_t)activationGradientsBackward.GetNcols(); t++) {
         activationGradientsBackward[j][t] = 0;
      }

      // set values
      for (size_t k = 0; k < nLocalViews; k++) {
         AReal grad = activationGradients[j][k];
         size_t winningIdx = indexMatrix[j][k];
         activationGradientsBackward[j][winningIdx] += grad;
      }
   }
}

//______________________________________________________________________________
template <typename AReal>
void TReference<AReal>::Reshape(TMatrixT<AReal> &A, const TMatrixT<AReal> &B)
{
   auto nColsA = A.GetNcols();
   auto nColsB = B.GetNcols();

   for (Int_t i = 0; i < A.GetNrows(); i++) {
      for (Int_t j = 0; j < A.GetNcols(); j++) {
         auto nElem = i * nColsA + j;
         A(i, j) = B(nElem / nColsB, nElem % nColsB);
      }
   }
}

//______________________________________________________________________________
template <typename AReal>
void TReference<AReal>::Flatten(TMatrixT<AReal> &A, const std::vector<TMatrixT<AReal>> &B, size_t size, size_t nRows,
                                size_t nCols)
{
   for (size_t i = 0; i < (size_t)size; i++) {
      for (size_t j = 0; j < (size_t)nRows; j++) {
         for (size_t k = 0; k < (size_t)nCols; k++) {
            A(i, j * nCols + k) = B[i](j, k);
         }
      }
   }
}

//______________________________________________________________________________
template <typename AReal>
void TReference<AReal>::Deflatten(std::vector<TMatrixT<AReal>> &A, const TMatrixT<AReal> &B, size_t size, size_t nRows,
                                  size_t nCols)
{
   for (size_t i = 0; i < (size_t)size; i++) {
      for (size_t j = 0; j < (size_t)nRows; j++) {
         for (size_t k = 0; k < (size_t)nCols; k++) {
            A[i](j, k) = B(i, j * nCols + k);
         }
      }
   }
}

//______________________________________________________________________________
template <typename AReal>
void TReference<AReal>::Rearrange(std::vector<TMatrixT<AReal>> &out, const std::vector<TMatrixT<AReal>> &in)
{
   // B x T x D out --- T x B x D in*/
   auto B = out.size();
   auto T = out[0].GetNrows();
   auto D = out[0].GetNcols();
   if ((T != (Int_t)in.size()) || (Int_t(B) != in[0].GetNrows()) || (D != in[0].GetNcols())) {
      std::cout << "Incompatible Dimensions\n"
                << in.size() << "x" << in[0].GetNrows() << "x" << in[0].GetNcols() << " --> " << B << "x" << T << "x"
                << D << "\n";
      return;
   }
   for (size_t i = 0; i < B; ++i) {
      for (Int_t j = 0; j < T; ++j) {
         for (Int_t k = 0; k < D; ++k) {
            out[i](j, k) = in[j](i, k);
         }
      }
   }
   return;
}

} // namespace DNN
} // namespace TMVA
