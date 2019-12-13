// @(#)root/tmva/tmva/cnn:$Id$
// Author: Vladimir Ilievski

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  :                                                                       *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Testing Conv Net Features                                                 *
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

#ifndef TMVA_TEST_DNN_TEST_CNN_TEST_CONV_NET_H
#define TMVA_TEST_DNN_TEST_CNN_TEST_CONV_NET_H

#include "../Utility.h"

#include "TMVA/DNN/Functions.h"
#include "TMVA/DNN/DeepNet.h"

using namespace TMVA::DNN;
using namespace TMVA::DNN::CNN;

enum ETestType { kLinearNet, kRndmActNet };

/** Testing the image to column function. Check wheether the matrix A after
 *  the transformation will be equal to the matrix B. */
//______________________________________________________________________________
template <typename Architecture>
auto testIm2col(typename Architecture::Matrix_t &A, typename Architecture::Matrix_t &B, size_t imgHeight,
                size_t imgWidth, size_t fltHeight, size_t fltWidth, size_t strideRows, size_t strideCols,
                size_t zeroPaddingHeight, size_t zeroPaddingWidth) -> bool
{

   size_t m, n;
   m = B.GetNrows();
   n = B.GetNcols();

   typename Architecture::Matrix_t ATr(m, n);
   Architecture::Im2col(ATr, A, imgHeight, imgWidth, fltHeight, fltWidth, strideRows, strideCols, zeroPaddingHeight,
                        zeroPaddingWidth);

   for (size_t i = 0; i < m; i++) {
      for (size_t j = 0; j < n; j++) {
         if (ATr(i, j) != B(i, j)) {
            return false;
         }
      }
   }

   return true;
}

/** Testing the rotation of the weights function. Check whether the rotated
 *  weight matrix A, will be equal to the matrix B. */
//______________________________________________________________________________
template <typename Architecture>
auto testRotateWeights(typename Architecture::Matrix_t &A, typename Architecture::Matrix_t &B, size_t filterDepth,
                       size_t filterHeight, size_t filterWidth, size_t numFilters) -> bool
{

   size_t m, n;
   m = B.GetNrows();
   n = B.GetNcols();

   typename Architecture::Matrix_t ARot(m, n);
   Architecture::RotateWeights(ARot, A, filterDepth, filterHeight, filterWidth, numFilters);

   for (size_t i = 0; i < m; i++) {
      for (size_t j = 0; j < n; j++) {
         if (ARot(i, j) != B(i, j)) {
            return false;
         }
      }
   }

   return true;
}

/** Downsample the matrix A and check whether the downsampled version
 *  is equal to B, and if the winning indices are equal to the matrix ind. */
//______________________________________________________________________________
template <typename Architecture>
auto testDownsample(const typename Architecture::Matrix_t &A, const typename Architecture::Matrix_t &ind,
                    const typename Architecture::Matrix_t &B, size_t imgHeight, size_t imgWidth, size_t fltHeight,
                    size_t fltWidth, size_t strideRows, size_t strideCols) -> bool
{

   size_t m1, n1;
   m1 = B.GetNrows();
   n1 = B.GetNcols();

   typename Architecture::Matrix_t ADown(m1, n1);

   size_t m2, n2;
   m2 = ind.GetNrows();
   n2 = ind.GetNcols();

   typename Architecture::Matrix_t AInd(m2, n2);

   Architecture::Downsample(ADown, AInd, A, imgHeight, imgWidth, fltHeight, fltWidth, strideRows, strideCols);

   for (size_t i = 0; i < m1; i++) {
      for (size_t j = 0; j < n1; j++) {
         if (ADown(i, j) != B(i, j)) {
            return false;
         }
      }
   }

   for (size_t i = 0; i < m2; i++) {
      for (size_t j = 0; j < n2; j++) {
         if (AInd(i, j) != ind(i, j)) {
            return false;
         }
      }
   }

   return true;
}

/** Back propagate the activation gradients through the max-pooling layer and check whether the
+ * computed gradients are equal to the matrix A. */
//______________________________________________________________________________
template <typename Architecture>
auto testPoolingBackward(const typename Architecture::Matrix_t &input, const typename Architecture::Matrix_t &output,
                         const typename Architecture::Matrix_t &indexMatrix, size_t imgHeight, size_t imgWidth,
                         size_t fltHeight, size_t fltWidth, size_t strideRows, size_t strideCols, size_t nLocalViews,
                         double epsilon = 0.01) -> bool {
    size_t depth = output.GetNrows();

    typename Architecture::Matrix_t ABack(output.GetNrows(), output.GetNcols());
    Architecture::MaxPoolLayerBackward(ABack, input, indexMatrix, imgHeight, imgWidth, fltHeight, fltWidth,
                                       strideRows, strideCols, nLocalViews);

    /* Needed to support double (almost) equality */
    auto almostEqual = [epsilon](double a, double b) {
        // Using a magic EPSILON value (makes sense for the existing tests).
        return fabs(a - b) < epsilon;
    };

    for (size_t d = 0; d < depth; d++) {
        for (size_t i = 0; i < nLocalViews; i++) {
            if (!almostEqual(ABack(d, i), output(d, i))) return false;
        }
    }
    return true;
}

/** Upsample the matrix A and check whether the upsampled version
+ *  is equal to B, and if the winning indices are equal to the matrix ind. */
//______________________________________________________________________________
template <typename Architecture>
auto testUpsample(const std::vector<typename Architecture::Matrix_t> &input, const std::vector<typename Architecture::Matrix_t> &output) -> bool 
{

  for(size_t v_i = 0; v_i < input.size(); v_i++){
    size_t m1,n1;
    m1 = input[v_i].GetNrows();
    n1 = input[v_i].GetNcols();

    size_t m2,n2;
    m2 = output[v_i].GetNrows();
    n2 = output[v_i].GetNcols();

    if ( m1 <= m2 && n1 <= n2){
      typename Architecture::Matrix_t AUp(m2,n2);

      std::vector< typename Architecture::Matrix_t > AUpVector;
      AUpVector.emplace_back(AUp);

      std::vector< typename Architecture::Matrix_t > AOutput;

      Architecture::Upsample(AUpVector,input[v_i]);

      for ( size_t i = 0; i < m2; i++){
        for ( size_t j = 0; j < n2; j++){
          if ( AUpVector[0](i,j) != output[v_i](i,j) ){
            return false;
          }
        }
      }
    }
    else{
      return false;
    }

  }
  return true;
}

/** Back propagate the activation gradients through the upsample layer and check whether the
+ * computed gradients are equal to the matrix A. */
//______________________________________________________________________________
template <typename Architecture>
auto testUpsampleBackward(const std::vector<typename Architecture::Matrix_t> &input, const std::vector<typename Architecture::Matrix_t> &output) -> bool 
{

  for(size_t v_i = 0; v_i < input.size(); v_i++){

    size_t m1,n1;
    m1 = input[v_i].GetNrows();
    n1 = input[v_i].GetNcols();

    size_t m2,n2;
    m2 = output[v_i].GetNrows();
    n2 = output[v_i].GetNcols();

    if ( m1 >= m2 && n1 >= n2){
      typename Architecture::Matrix_t AUp(m2,n2);

      std::vector< typename Architecture::Matrix_t > AVector;
      AVector.emplace_back(AUp);

      std::vector< typename Architecture::Matrix_t > AOutput;

      Architecture::UpsampleLayerBackward(AVector,input[v_i]);

      for ( size_t i = 0; i < m2; i++){
        for ( size_t j = 0; j < n2; j++){
          if ( AVector[0](i,j) != output[v_i](i,j) ){
            return false;
          }
        }
      }
    }
    else{
      return false;
    }
  }
  return true;
}

/** Reshape the matrix A using the Reshape function and compare it to
 *  the result in matrix B. */
//______________________________________________________________________________
template <typename Architecture_t>
auto testReshape(const typename Architecture_t::Matrix_t &A, const typename Architecture_t::Matrix_t &B) -> bool
{

    size_t m, n;
    m = B.GetNrows();
    n = B.GetNcols();

    typename Architecture_t::Matrix_t AReshaped(m, n);
    Architecture_t::Reshape(AReshaped, A);

    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            if (AReshaped(i, j) != B(i, j)) {
                return false;
            }
        }
    }
    return true;
}

/** Flatten the 3D tensor A using the Flatten function and compare it to
 *  the result in the flat matrix B. */
//______________________________________________________________________________
template <typename Architecture_t>
auto testFlatten(std::vector<typename Architecture_t::Matrix_t> &A, const typename Architecture_t::Matrix_t &B, size_t size,
                 size_t nRows, size_t nCols) -> bool
{

   size_t m, n;
   m = B.GetNrows();
   n = B.GetNcols();

   typename Architecture_t::Matrix_t AFlat(m, n);
   Architecture_t::Flatten(AFlat, A, size, nRows, nCols);

   for (size_t i = 0; i < m; i++) {
      for (size_t j = 0; j < n; j++) {
         if (AFlat(i, j) != B(i, j)) {
            return false;
         }
      }
   }

   return true;
}

template <typename Architecture>
auto testConvLayerForward(const std::vector<typename Architecture::Matrix_t> &input,
                          const std::vector<typename Architecture::Matrix_t> &expectedOutput,
                          const typename Architecture::Matrix_t &weights, const typename Architecture::Matrix_t &biases,
                          size_t inputHeight, size_t inputWidth, size_t inputDepth, size_t fltHeight,
                          size_t fltWidth, size_t numberFilters, size_t strideRows, size_t strideCols,
                          size_t zeroPaddingHeight, size_t zeroPaddingWidth) -> bool
{
    size_t nRows = expectedOutput[0].GetNrows();
    size_t nCols = expectedOutput[0].GetNcols();
    // batchSize == 1.
    std::vector<typename Architecture::Matrix_t> computedOutput;
    computedOutput.emplace_back(nRows, nCols);

    std::vector<typename Architecture::Matrix_t> computedDerivatives;
    computedDerivatives.emplace_back(nRows, nCols);

    TConvParams params(1, inputDepth, inputHeight, inputWidth, numberFilters, fltHeight, fltWidth, strideRows,
                       strideCols, zeroPaddingHeight, zeroPaddingWidth);


    size_t height = (inputHeight - fltHeight + 2 * zeroPaddingHeight) / strideRows + 1;
    size_t width = (inputWidth - fltWidth + 2 * zeroPaddingWidth) / strideCols + 1;
    size_t nLocalViews = height * width;
    size_t nLocalViewPixels = inputDepth * fltHeight * fltWidth;

    std::vector<typename Architecture::Matrix_t> forwardMatrices;
    forwardMatrices.emplace_back(nLocalViews, nLocalViewPixels);

    Architecture::ConvLayerForward(computedOutput, computedDerivatives, input, weights, biases, params,
                                   EActivationFunction::kIdentity, forwardMatrices);

    for (size_t slice = 0; slice < nRows; slice++) {
        for (size_t localView = 0; localView < nCols; localView++) {
            if (expectedOutput[0](slice, localView) != computedOutput[0](slice, localView)) return false;
        }
    }
    return true;
}

template <typename Architecture>
auto testTransConvLayerForward(const std::vector<typename Architecture::Matrix_t> &input,
                          const std::vector<typename Architecture::Matrix_t> &expectedOutput,
                          const typename Architecture::Matrix_t &weights, const typename Architecture::Matrix_t &biases,
                          size_t inputHeight, size_t inputWidth, size_t inputDepth, size_t fltHeight,
                          size_t fltWidth, size_t numberFilters, size_t strideRows, size_t strideCols,
                          size_t zeroPaddingHeight, size_t zeroPaddingWidth) -> bool
{
   for(size_t v_i = 0; v_i < expectedOutput.size(); v_i++)
   {
    size_t nRows = expectedOutput[v_i].GetNrows();
    size_t nCols = expectedOutput[v_i].GetNcols();
    // batchSize == 1.
    std::vector<typename Architecture::Matrix_t> computedOutput;
    computedOutput.emplace_back(nRows, nCols);

    std::vector<typename Architecture::Matrix_t> computedDerivatives;
    computedDerivatives.emplace_back(nRows, nCols);

    TConvParams params(1, inputDepth, inputHeight, inputWidth, numberFilters, fltHeight, fltWidth, strideRows,
                       strideCols, zeroPaddingHeight, zeroPaddingWidth);


    size_t height = ( ( inputHeight - 1)*strideRows )- 2* zeroPaddingHeight +  fltHeight;
    size_t width = ( ( inputWidth - 1)*strideCols )- 2* zeroPaddingWidth +  fltWidth;
    size_t nLocalViews = height * width;
    size_t nLocalViewPixels = inputDepth * fltHeight * fltWidth;

    std::vector<typename Architecture::Matrix_t> forwardMatrices;
    forwardMatrices.emplace_back(nLocalViews, nLocalViewPixels);

    Architecture::TransConvLayerForward(computedOutput, computedDerivatives, input, weights, biases, params,
                                   EActivationFunction::kIdentity, forwardMatrices);

    for (size_t slice = 0; slice < nRows; slice++) {
        for (size_t localView = 0; localView < nCols; localView++) {
            if (expectedOutput[v_i](slice, localView) != computedOutput[v_i](slice, localView)) return false;
        }
    }
  }
    return true;
}

template <typename Architecture>
auto testTransConvLayerBackward(const std::vector<typename Architecture::Matrix_t> &input,
                          const std::vector<typename Architecture::Matrix_t> &expectedOutput,
                          const typename Architecture::Matrix_t &weights, const typename Architecture::Matrix_t &biases,
                          size_t inputHeight, size_t inputWidth, size_t inputDepth, size_t fltHeight,
                          size_t fltWidth, size_t numberFilters, size_t strideRows, size_t strideCols,
                          size_t zeroPaddingHeight, size_t zeroPaddingWidth) -> bool
{
   for(size_t v_i = 0 ; v_i < expectedOutput.size(); v_i++)
   {
    size_t nRows = expectedOutput[v_i].GetNrows();
    size_t nCols = expectedOutput[v_i].GetNcols();
    // batchSize == 1.
    std::vector<typename Architecture::Matrix_t> computedOutput;
    computedOutput.emplace_back(nRows, nCols);

    std::vector<typename Architecture::Matrix_t> computedDerivatives;
    computedDerivatives.emplace_back(nRows, nCols);

    TConvParams params(1, inputDepth, inputHeight, inputWidth, numberFilters, fltHeight, fltWidth, strideRows,
                       strideCols, zeroPaddingHeight, zeroPaddingWidth);


    size_t height = ( ( inputHeight - 1)*strideRows )- 2* zeroPaddingHeight +  fltHeight;
    size_t width = ( ( inputWidth - 1)*strideCols )- 2* zeroPaddingWidth +  fltWidth;
    size_t nLocalViews = height * width;
    size_t nLocalViewPixels = inputDepth * fltHeight * fltWidth;

    std::vector<typename Architecture::Matrix_t> forwardMatrices;
    forwardMatrices.emplace_back(nLocalViews, nLocalViewPixels);

    Architecture::TransConvLayerBackward(computedOutput, computedDerivatives, input, weights, biases, params,
                                   EActivationFunction::kIdentity, forwardMatrices);

    for (size_t slice = 0; slice < nRows; slice++) {
        for (size_t localView = 0; localView < nCols; localView++) {
            if (expectedOutput[v_i](slice, localView) != computedOutput[v_i](slice, localView)) return false;
        }
    }

  }
    return true;
}

/** Deflatten the 2D tensor A using the Deflatten function and compare it to
 *  the result in the 3D tensor B. */
//______________________________________________________________________________
template <typename Architecture_t>
auto testDeflatten(const typename Architecture_t::Matrix_t &A, const std::vector<typename Architecture_t::Matrix_t> &B,
                   size_t size, size_t nRows, size_t nCols) -> bool
{
    std::vector<typename Architecture_t::Matrix_t> AComputed;
    for (size_t i = 0; i < size; i++) {
        AComputed.emplace_back(nRows, nCols);
    }

    Architecture_t::Deflatten(AComputed, A, size, nRows, nCols);

    for (size_t i = 0; i < size; i++) {
        for (size_t j = 0; j < nRows; j++) {
            for (size_t k = 0; k < nCols; k++) {
                if (AComputed[i](j, k) != B[i](j, k)) return false;
            }
        }
    }

    return true;
}

/*! Generate a conv net, perform forward pass */
//______________________________________________________________________________
template <typename Architecture>
auto testConvForwardPass(size_t batchSize, size_t imgDepth, size_t imgHeight, size_t imgWidth, size_t batchDepth,
                         size_t batchHeight, size_t batchWidth) -> void
{
   using Matrix_t = typename Architecture::Matrix_t;
   using Net_t = TDeepNet<Architecture>;

   Net_t convNet(batchSize, imgDepth, imgHeight, imgWidth, batchDepth, batchHeight, batchWidth,
                 ELossFunction::kMeanSquaredError, EInitialization::kGauss);
   constructConvNet(convNet);
   convNet.Initialize();

   std::vector<Matrix_t> X;
   for (size_t i = 0; i < batchSize; i++) {
      X.emplace_back(imgDepth, imgHeight * imgWidth);
      randomMatrix(X[i]);
   }

   convNet.Forward(X);
}

/*! Generate a conv net, get the loss. */
//______________________________________________________________________________
template <typename Architecture>
auto testConvLossFunction(size_t batchSize, size_t imgDepth, size_t imgHeight, size_t imgWidth, size_t batchDepth,
                          size_t batchHeight, size_t batchWidth) -> void
{
   using Scalar_t = typename Architecture::Scalar_t;
   using Matrix_t = typename Architecture::Matrix_t;
   using Net_t = TDeepNet<Architecture>;

   Net_t convNet(batchSize, imgDepth, imgHeight, imgWidth, batchDepth, batchHeight, batchWidth,
                 ELossFunction::kMeanSquaredError, EInitialization::kGauss);
   constructConvNet(convNet);
   convNet.Initialize();

   std::vector<Matrix_t> X;
   for (size_t i = 0; i < batchSize; i++) {
      X.emplace_back(imgDepth, imgHeight * imgWidth);
      randomMatrix(X[i]);
   }

   Matrix_t Y(batchSize, convNet.GetOutputWidth());
   Matrix_t W(batchSize, 1);
   randomMatrix(Y);
   randomMatrix(W);

   Scalar_t loss = convNet.Loss(X, Y, W);
   std::cout << "The loss is: " << loss << std::endl;
}

/*! Generate a conv net, get the predictios */
//______________________________________________________________________________
template <typename Architecture>
auto testConvPrediction(size_t batchSize, size_t imgDepth, size_t imgHeight, size_t imgWidth, size_t batchDepth,
                        size_t batchHeight, size_t batchWidth, EOutputFunction f) -> void
{
   using Matrix_t = typename Architecture::Matrix_t;
   using Net_t = TDeepNet<Architecture>;

   Net_t convNet(batchSize, imgDepth, imgHeight, imgWidth, batchDepth, batchHeight, batchWidth,
                 ELossFunction::kMeanSquaredError, EInitialization::kGauss);
   constructConvNet(convNet);
   convNet.Initialize();

   std::vector<Matrix_t> X;
   for (size_t i = 0; i < batchSize; i++) {
      X.emplace_back(imgDepth, imgHeight * imgWidth);
      randomMatrix(X[i]);
   }

   Matrix_t Predictions(batchSize, convNet.GetOutputWidth());
   convNet.Prediction(Predictions, X, f);

   for (size_t i = 0; i < batchSize; i++) {
      for (size_t j = 0; j < convNet.GetOutputWidth(); j++) {
         std::cout << Predictions(i, j) << " ";
      }
      std::cout << "" << std::endl;
   }
}

/*! Generate a conv net, test the backward pass, always with stride 1. */
//______________________________________________________________________________
template <typename Architecture>
auto testConvBackwardPassOnly(size_t batchSize, size_t imgDepth, size_t imgHeight, size_t imgWidth, size_t batchDepth,
                          size_t batchHeight, size_t batchWidth) -> void
{
   using Matrix_t = typename Architecture::Matrix_t;
   using Net_t = TDeepNet<Architecture>;

   Net_t convNet(batchSize, imgDepth, imgHeight, imgWidth, batchDepth, batchHeight, batchWidth,
                 ELossFunction::kMeanSquaredError, EInitialization::kGauss);
   constructConvNet(convNet);
   convNet.Initialize();

   std::vector<Matrix_t> X;
   for (size_t i = 0; i < batchSize; i++) {
      X.emplace_back(imgDepth, imgHeight * imgWidth);
      randomMatrix(X[i]);
   }

   Matrix_t Y(batchSize, convNet.GetOutputWidth());
   Matrix_t W(batchSize, 1);
   randomMatrix(Y);
   randomMatrix(W);

   convNet.Forward(X);
   convNet.Backward(X, Y, W);
}


/*! Compute the loss of the net as a function of the weight at index (i,j) in
 *  layer l. dx is added as an offset to the current value of the weight. */
//______________________________________________________________________________
template <typename Architecture>
auto evaluate_net_weight(TDeepNet<Architecture> &net, std::vector<typename Architecture::Matrix_t> &X,
                         const typename Architecture::Matrix_t &Y, const typename Architecture::Matrix_t &W, size_t l,
                         size_t i, size_t j, size_t k, typename Architecture::Scalar_t dx) -> typename Architecture::Scalar_t
{
    using Scalar_t = typename Architecture::Scalar_t;
    //using Matrix_t = typename Architecture::Matrix_t;

    // shift the weight value and compute the Loss
    auto & netW = net.GetLayerAt(l)->GetWeights();
    netW[k](i,j) += dx;
    Scalar_t res = net.Loss(X, Y, W);
    // rest weight to original value
    netW[k](i,j) -= dx;
    //std::cout << "loss(w+dx = " << res << " loss(w) " << net.Loss(X,Y,W) << std::endl;
    return res;
}

/*! Compute the loss of the net as a function of the weight at index i in
 *  layer l. dx is added as an offset to the current value of the weight. */
//______________________________________________________________________________
template <typename Architecture>
auto evaluate_net_bias(TDeepNet<Architecture> &net, std::vector<typename Architecture::Matrix_t> &X,
                       const typename Architecture::Matrix_t &Y, const typename Architecture::Matrix_t &W, size_t l,
                       size_t i, size_t k, typename Architecture::Scalar_t dx) -> typename Architecture::Scalar_t
{
    using Scalar_t = typename Architecture::Scalar_t;
    //using Matrix_t = typename Architecture::Matrix_t;

    auto & netB = net.GetLayerAt(l)->GetBiases();
    netB[k](i,0) += dx;
    Scalar_t res = net.Loss(X, Y, W);
    netB[k](i,0) -= dx;
    return res;
}



template <typename Architecture>
auto testConvBackwardPass(size_t batchSize, size_t imgDepth, size_t imgHeight, size_t imgWidth, size_t batchDepth,
                          size_t batchHeight, size_t batchWidth, typename Architecture::Scalar_t dx, ETestType testType) -> bool
{
   using Matrix_t = typename Architecture::Matrix_t;
   using Net_t = TDeepNet<Architecture>;
   using Scalar_t = typename Architecture::Scalar_t;

   Net_t convNet(batchSize, imgDepth, imgHeight, imgWidth, batchDepth, batchHeight, batchWidth,
                 ELossFunction::kMeanSquaredError, EInitialization::kGauss);
   // tyoe of network
   if (testType == kLinearNet) 
      constructLinearConvNet(convNet);
   else
      constructConvNet(convNet);

   convNet.Initialize();

//    auto & w0 = convNet.GetLayerAt(0)->GetWeights();
// #ifdef DEBUG  
//    std::cout << "Netwrok weights for Layer 0  " << std::endl; 
//    std::cout << "weight depth = " << w0.size() << std::endl;
//    for (size_t i = 0; i < w0.size();   ++i)
//       TMVA_DNN_PrintTCpuMatrix(w0[i],"weight-layer0");
// #endif  
   
   std::vector<Matrix_t> X;
   for (size_t i = 0; i < batchSize; i++) {
      X.emplace_back(batchHeight , batchWidth);
      randomMatrix(X[i]);

       // print input
#ifdef DEBUGH  
      std::cout << "INPUT - batch " << i << std::endl;
      TMVA_DNN_PrintTCpuMatrix(X[i],"input");
#endif  
   }

   Matrix_t Y(batchSize, convNet.GetOutputWidth());
   Matrix_t W(batchSize, 1);   // this are the data weights 
   randomMatrix(Y);
   //randomMatrix(W);
   // for the moment use weights equal to 1
   fillMatrix(W, 1.0);

   std::cout << "Do Forward Pass " << std::endl;
   convNet.Forward(X);

   // print layer derivatives
#ifdef DEBUG
   using ConvLayer_t = TConvLayer<Architecture>;
   auto convLayer = dynamic_cast<ConvLayer_t*>(convNet.GetLayerAt(0) );
   if (convLayer) { 
      auto & df = convLayer->GetDerivatives();
      std::cout << "Derivatives - size " << df.size() << std::endl;
      for (size_t ii=0; ii< df.size(); ++ii)
         TMVA_DNN_PrintTCpuMatrix(df[ii],"Derivatives");
   }
#endif

   //if (testType == kRndmActNet)  return true; 

   std::cout << "Do Backward Pass " << std::endl;
   convNet.Backward(X, Y, W);


   // now compare derivatives using finite differences and compare the result
    Scalar_t maximum_error = 0.0;

    for (size_t l = 0; l < convNet.GetDepth(); l++) {
      std::cout << "\rTesting weight gradients:      layer: " << l << " / " << convNet.GetDepth();
      std::cout << std::flush;
      auto & layer = *(convNet.GetLayerAt(l));
      std::vector<Matrix_t> &gw = layer.GetWeightGradients();

      std::cout << std::endl;

      if (gw.size() > 0) { 
         std::cout << "Weight gradient from back-propagation - vector size is " << gw.size()  << std::endl;

         if (gw[0].GetNoElements() < 100 ) {
            gw[0].Print();
         }
         else
            std::cout << "BP Weight Gradient ( " << gw[0].GetNrows() << " x " << gw[0].GetNcols() << " ) , ...... skip printing (too many elements ) " << std::endl;  
      }
      else {
         std::cout << "Layer " << l << " has no weights " << std::endl;
      }
      auto & actGrad = layer.GetActivationGradients();
      if (actGrad.size() > 0)  {
         std::cout << "Activation gradient from back-propagation  - vector size is " << actGrad.size() << std::endl;
         if (actGrad[0].GetNoElements() < 100 ) { 
            for (size_t ii = 0; ii < actGrad.size(); ++ii) 
               actGrad[ii].Print();
         } else
            std::cout << "Activation Gradient ( " << actGrad[0].GetNrows() << " x " << actGrad[0].GetNcols() << " ) , ...... skip printing (too many elements ) " << std::endl;
      }

      std::cout << "Layer " << l << " :  output  D x H x W " << layer.GetDepth() << "  " << layer.GetHeight() << "  " << layer.GetWidth()
                << "\t input D x H x W " << layer.GetInputDepth() << "  " << layer.GetInputHeight() << "  " << layer.GetInputWidth() << std::endl;
      

      // print output
      auto & outL = layer.GetOutput();
      std::cout << "layer output size " << outL.size() << std::endl;
      if (outL.size() > 0) {
         if (outL[0].GetNoElements() < 100 ) { 
            outL[0].Print();
         } else
            std::cout << "Layer Output ( " << outL[0].GetNrows() << " x " << outL[0].GetNcols() << " ) , ...... skip printing (too many elements ) " << std::endl;
      }
      
      std::cout << "Evaluate the Derivatives with Finite difference and compare with BP for Layer " << l << std::endl;
      int nerrors = 0;
      int ngood = 0;
#ifdef DEBUG
      int ngoodPrint = 10000;
#else
      int ngoodPrint = 3;
#endif 
      for (size_t k = 0; k <  gw.size() ; ++k) { 
         //for (size_t i = 0; i < layer.GetWidth(); i++) {
         for (size_t i = 0; i < gw[k].GetNrows(); i++) {
            //for (size_t j = 0; j < layer.GetInputWidth(); j++) {
            for (size_t j = 0; j < gw[k].GetNcols(); j++) {
               auto f = [&convNet, &X, &Y, &W, l, i, j, k](Scalar_t x) {
                  return evaluate_net_weight(convNet, X, Y, W, l, i, j, k, x);
               };
               Scalar_t dy = finiteDifference(f, dx) / (2.0 * dx);
               Scalar_t dy_ref = gw[k](i, j);
               // Compute the relative error if dy != 0.
               Scalar_t error;
               if (std::fabs(dy_ref) > 1e-15) {
                  error = std::fabs((dy - dy_ref) / dy_ref);
               } else {
                  error = std::fabs(dy - dy_ref);
               }
               maximum_error = std::max(error, maximum_error);
               if (error > 1.E-3) {
                  std::cout << k << " - " <<  i << " , " << j << " : " << dy << " from BP " << dy_ref << "   " << error << " ERROR " << std::endl;
                  nerrors ++; 
               }
               else {
                  if (ngood < ngoodPrint) std::cout << k << " - " <<  i << " , " << j << " : " << dy << " from BP " << dy_ref << "   " << error << std::endl;
                  ngood++;
               }
               if (nerrors > 10) {
                  std::cout << "Reached error limit skip..." << std::endl;
                  break;
               }
            }
            if (nerrors > 10) break; 
         }
         if (nerrors > 10) break; 
      }
    }
   std::cout << "\rTesting weight gradients:      ";
   std::cout << "maximum relative error: " << print_error(maximum_error) << std::endl;
   if (maximum_error > 1.E-3) {
      std::cout << "ERROR - BackPropagation test failed in computing  weight Derivatives " << std::endl;
      return false; 
   }
   //return maximum_error;
   return true; 
}

#endif
