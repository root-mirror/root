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

#include <string>
#include "TMath.h"
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

   typename Architecture::Tensor_t tDown(ADown,3); // convert to tensors of dims 3
   typename Architecture::Tensor_t tInd(AInd,3);
   typename Architecture::Tensor_t tA(A,3);

   std::cout << "Testing downsample with size = " << fltHeight << " , " << fltWidth
   << " stride " << strideRows << " . " << strideCols << std::endl;

   TDescriptors * poolDescriptors = nullptr;
   TWorkspace   * poolWorkspace   = nullptr;
   TConvParams params(1, imgHeight, imgWidth, imgWidth, 1, fltHeight, fltWidth, strideRows,
                      strideCols, 0, 0);

   TMaxPoolLayer<Architecture> *layer = nullptr;
   Architecture::InitializePoolDescriptors(poolDescriptors, layer);
   Architecture::InitializePoolDropoutWorkspace(poolWorkspace, poolDescriptors, params, layer);

   Architecture::Downsample(tDown, tInd, tA,
                            (typename Architecture::PoolingDescriptors_t &) *poolDescriptors,
                            (typename Architecture::PoolingWorkspace_t &) *poolWorkspace,
                            imgHeight, imgWidth, fltHeight, fltWidth, strideRows, strideCols);

   for (size_t i = 0; i < m1; i++) {
      for (size_t j = 0; j < n1; j++) {
         if (ADown(i, j) != B(i, j)) {
            std::cout << "Error - downsample failed for " << i << "," << j << std::endl;
            Architecture::PrintTensor(tDown,"downsample tensor");
            typename Architecture::Tensor_t tB(B,3);
            Architecture::PrintTensor(tB,"expected tensor");
            Architecture::PrintTensor(tA,"input tensor");
            return false;
         }
      }
   }

   for (size_t i = 0; i < m2; i++) {
      for (size_t j = 0; j < n2; j++) {
         if (AInd(i, j) != ind(i, j)) {
            std::cout << "Error - downsample index failed for " << i << "," << j << std::endl;
            Architecture::PrintTensor(tInd,"index tensor");
            typename Architecture::Tensor_t texInd(ind,3);
            Architecture::PrintTensor(texInd,"expected tensor");
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

    typename Architecture::Tensor_t ABack(1,output.GetNrows(), output.GetNcols());
    typename Architecture::Tensor_t tInput( input, 3);
    typename Architecture::Tensor_t tIndexMatrix( indexMatrix, 3);
    // only needed in cuDNN backward pass
    typename Architecture::Tensor_t inputActivation;
    typename Architecture::Tensor_t outputActivation;

    TDescriptors * poolDescriptors = nullptr;
    TWorkspace   * poolWorkspace   = nullptr;
    TConvParams params(1, imgHeight, imgWidth, imgWidth, 1, fltHeight, fltWidth, strideRows,
                      strideCols, 0, 0);

    TMaxPoolLayer<Architecture> *layer = nullptr;
    Architecture::InitializePoolDescriptors(poolDescriptors, layer);
    Architecture::InitializePoolDropoutWorkspace(poolWorkspace, poolDescriptors, params, layer);

    Architecture::MaxPoolLayerBackward(ABack, tInput, tIndexMatrix, inputActivation, outputActivation,
                                       (typename Architecture::PoolingDescriptors_t &) *poolDescriptors,
                                       (typename Architecture::PoolingWorkspace_t &) *poolWorkspace,
                                       imgHeight, imgWidth, fltHeight, fltWidth,
                                       strideRows, strideCols, nLocalViews);

    /* Needed to support double (almost) equality */
    auto almostEqual = [epsilon](double a, double b) {
        // Using a magic EPSILON value (makes sense for the existing tests).
        return fabs(a - b) < epsilon;
    };

    for (size_t d = 0; d < depth; d++) {
        for (size_t i = 0; i < nLocalViews; i++) {
            if (!almostEqual(ABack(0, d, i), output(d, i))) return false;
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
auto testFlatten(typename Architecture_t::Tensor_t &A, const typename Architecture_t::Tensor_t &B ) -> bool
{

   size_t m, n;
   m = B.GetHSize();
   n = B.GetWSize();

   typename Architecture_t::Tensor_t AFlat( B.GetShape() );

   Architecture_t::Flatten(AFlat, A );

   for (size_t i = 0; i < m; i++) {
      for (size_t j = 0; j < n; j++) {
         if (AFlat(0, i, j) != B(0, i, j)) {
            std::cout << "Error - test flatten failed for element " << i << "  " << j << std::endl;
            Architecture_t::PrintTensor(AFlat,"Flatten tensor");
            Architecture_t::PrintTensor(A,"Input tensor");
            return false;
         }
      }
   }

   return true;
}

template <typename Architecture>
auto testConvLayerForward(const typename Architecture::Tensor_t &input,
                          const typename Architecture::Tensor_t &expectedOutput,
                          const typename Architecture::Matrix_t &weights, const typename Architecture::Matrix_t &biases,
                          size_t inputHeight, size_t inputWidth, size_t inputDepth, size_t fltHeight,
                          size_t fltWidth, size_t numberFilters, size_t strideRows, size_t strideCols,
                          size_t zeroPaddingHeight, size_t zeroPaddingWidth) -> bool
{
    size_t nRows = expectedOutput.GetHSize();
    size_t nCols = expectedOutput.GetWSize();
    size_t batchSize = 1;

    typename Architecture::Tensor_t computedOutput( batchSize, nRows, nCols);
    typename Architecture::Tensor_t computedDerivatives(batchSize, nRows, nCols);

    TConvParams params(1, inputDepth, inputHeight, inputWidth, numberFilters, fltHeight, fltWidth, strideRows,
                       strideCols, zeroPaddingHeight, zeroPaddingWidth);


    size_t height = (inputHeight - fltHeight + 2 * zeroPaddingHeight) / strideRows + 1;
    size_t width = (inputWidth - fltWidth + 2 * zeroPaddingWidth) / strideCols + 1;
    size_t nLocalViews = height * width;
    size_t nLocalViewPixels = inputDepth * fltHeight * fltWidth;

    typename Architecture::Tensor_t forwardMatrices (1 , nLocalViews, nLocalViewPixels);

    TDescriptors * convDescriptors = nullptr;
    TWorkspace   * convWorkspace   = nullptr;

    TConvLayer<Architecture> *layer = nullptr;
    Architecture::InitializeConvDescriptors(convDescriptors, layer);
    Architecture::InitializeConvWorkspace(convWorkspace, convDescriptors, params, layer);

    Architecture::ConvLayerForward(computedOutput, computedDerivatives, input, weights, biases, params,
                                   EActivationFunction::kIdentity, forwardMatrices,
                                   (typename Architecture::ConvDescriptors_t &) *convDescriptors,
                                   (typename Architecture::ConvWorkspace_t &) *convWorkspace);

    for (size_t slice = 0; slice < nRows; slice++) {
        for (size_t localView = 0; localView < nCols; localView++) {
            if (expectedOutput(0, slice, localView) != computedOutput(0, slice, localView)) {
               std::cout << "Error - computed output different than expected for " << slice
               << " , " << localView << std::endl;
               Architecture::PrintTensor(computedOutput,"computed output tensor");
               Architecture::PrintTensor(expectedOutput,"expected output tensor");
               return false;
            }
        }
    }
    return true;
}

/** Deflatten the 2D tensor A using the Deflatten function and compare it to
 *  the result in the 3D tensor B. */
//______________________________________________________________________________
template <typename Architecture_t>
auto testDeflatten(const typename Architecture_t::Tensor_t &A, const typename Architecture_t::Tensor_t &B ) -> bool
{
   typename Architecture_t::Tensor_t AComputed( B.GetShape() );
   size_t size = B.GetFirstSize();
   size_t nRows = B.GetHSize();
   size_t nCols = B.GetWSize();

   Architecture_t::Deflatten(AComputed, A );

    for (size_t i = 0; i < size; i++) {
        for (size_t j = 0; j < nRows; j++) {
            for (size_t k = 0; k < nCols; k++) {
                if (AComputed(i, j, k) != B(i, j, k)) return false;
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
   //using Matrix_t = typename Architecture::Matrix_t;
   using Net_t = TDeepNet<Architecture>;

   Net_t convNet(batchSize, imgDepth, imgHeight, imgWidth, batchDepth, batchHeight, batchWidth,
                 ELossFunction::kMeanSquaredError, EInitialization::kGauss);
   constructConvNet(convNet);
   convNet.Initialize();

   //typename Architecture::Tensor_t X (batchSize, imgDepth, imgHeight * imgWidth);

   typename Architecture::Tensor_t X =  Architecture::CreateTensor(batchSize, imgDepth, imgHeight , imgWidth);
   randomBatch(X);


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

   typename Architecture::Tensor_t X(batchSize, imgDepth, imgHeight * imgWidth);

   randomBatch(X);


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

   typename Architecture::Tensor_t X(batchSize, imgDepth, imgHeight * imgWidth);
   randomBatch(X);

   Matrix_t Predictions(batchSize, convNet.GetOutputWidth());
   convNet.Prediction(Predictions, X, f);

   for (size_t i = 0; i < batchSize; i++) {
      for (size_t j = 0; j < convNet.GetOutputWidth(); j++) {
         std::cout << Predictions(i, j) << " ";
      }
      std::cout << "" << std::endl;
   }
}

/*! Generate a conv net, perform p */
//______________________________________________________________________________
template <typename Architecture1, typename Architecture2>
auto testMixedConvForwardPass(size_t batchSize, size_t imgDepth, size_t imgHeight, size_t imgWidth, size_t batchDepth,
                         size_t batchHeight, size_t batchWidth, bool debug) -> double
{
   // using Matrix_t = typename Architecture::Matrix_t;
   using Net1_t = TDeepNet<Architecture1>;
   using Net2_t = TDeepNet<Architecture2>;

   Net1_t convNet1(batchSize, imgDepth, imgHeight, imgWidth, batchDepth, batchHeight, batchWidth,
                 ELossFunction::kMeanSquaredError, EInitialization::kGauss);
   constructLinearConvNet(convNet1);
   convNet1.Initialize();

   Net2_t convNet2(batchSize, imgDepth, imgHeight, imgWidth, batchDepth, batchHeight, batchWidth,
                   ELossFunction::kMeanSquaredError, EInitialization::kGauss);
   constructLinearConvNet(convNet2);
   convNet2.Initialize();





   typename Architecture1::Tensor_t X =  Architecture1::CreateTensor(batchSize, imgDepth, imgHeight , imgWidth);

   randomBatch(X);

   typename Architecture2::Tensor_t X2 =  Architecture2::CreateTensor(batchSize, imgDepth, imgHeight , imgWidth);
   Architecture2::CopyDiffArch(X2, X);

   if (debug) {
      convNet1.Print();
      Architecture1::PrintTensor(X,"X input arch 1");
      Architecture2::PrintTensor(X2,"X input arch 2");
      // std::cout << "buffer of tensor2  - size " << X2.GetSize() << std::endl;
      // for (size_t i = 0; i < X2.GetSize(); ++i)
      //    std::cout << X2.GetDeviceBuffer()[i] << ", ";
      // std::cout << std::endl;

      std::cout << "Forward net1 " << std::endl;
   }

   convNet1.Forward(X);

   if (debug) {
      std::cout << "look at net2 " << std::endl;
      convNet2.Print();
   }

   // copy weights from net1 to net2
   for (size_t i = 0; i < convNet1.GetDepth(); ++i) {
      const auto &layer1 = convNet1.GetLayerAt(i);
      const auto &layer2 = convNet2.GetLayerAt(i);
      Architecture2::CopyDiffArch(layer2->GetWeights(), layer1->GetWeights());
      Architecture2::CopyDiffArch(layer2->GetBiases(), layer1->GetBiases());
   }
   if (debug) {
      Architecture1::PrintTensor(convNet1.GetLayerAt(0)->GetWeightsAt(0), " con layer 1 weights arch 1");
      Architecture2::PrintTensor(convNet2.GetLayerAt(0)->GetWeightsAt(0), " conv layer 1 weights arch 2 ");

      Architecture1::PrintTensor(convNet1.GetLayerAt(0)->GetBiasesAt(0), " conv layer 1 bias arch 1");
      Architecture2::PrintTensor(convNet2.GetLayerAt(0)->GetBiasesAt(0), " conv layer 1 bias arch 2 ");

      // typename Architecture2::Tensor_t t(convNet2.GetLayerAt(1)->GetWeightsAt(0));
      // std::cout << "buffer of tensor2  - size " << t.GetSize() << std::endl;
      // for (size_t i = 0; i < t.GetSize(); ++i)
      //    std::cout << t.GetDeviceBuffer()[i] << ", ";
      // std::cout << std::endl;
   }





   // examine the output
   convNet2.Forward(X2);

   if (debug) {
      for (size_t i = 0; i < convNet1.GetDepth()-1; ++i) {
         std::cout << "layer  " << i << std::endl;
         Architecture1::PrintTensor(convNet1.GetLayerAt(i)->GetOutput(), "output layer-i arch 1");
         Architecture2::PrintTensor(convNet2.GetLayerAt(i)->GetOutput(), "output layer-i arch 2");
      }
   }

   Architecture1::PrintTensor(convNet1.GetLayers().back()->GetOutput(), "output layer architecture 1");
   Architecture2::PrintTensor(convNet2.GetLayers().back()->GetOutput(), "output layer architecture 2");

   TMatrixT<typename Architecture1::Scalar_t> out1 = convNet1.GetLayers().back()->GetOutput();
   TMatrixT<typename Architecture2::Scalar_t> out2 = convNet2.GetLayers().back()->GetOutput();

   Double_t error = maximumRelativeError(out1, out2);
   return error;

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

   typename Architecture::Tensor_t X(batchSize, imgDepth, imgHeight * imgWidth);
   randomBatch(X);


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
auto evaluate_net_weight(TDeepNet<Architecture> &net, typename Architecture::Tensor_t &X,
                         const typename Architecture::Matrix_t &Y, const typename Architecture::Matrix_t &W, size_t l,
                         size_t i, size_t j, size_t k, typename Architecture::Scalar_t dx) -> typename Architecture::Scalar_t
{
    using Scalar_t = typename Architecture::Scalar_t;
    //using Matrix_t = typename Architecture::Matrix_t;

    // shift the weight value and compute the Loss
    auto & netW = net.GetLayerAt(l)->GetWeightsAt(k);

    netW(i,j) += dx;

    Scalar_t res = net.Loss(X, Y, W);
    // rest weight to original value
    netW(i,j) -= dx;
    //std::cout << "loss(w+dx = " << res << " loss(w) " << net.Loss(X,Y,W) << std::endl;
    return res;
}

/*! Compute the loss of the net as a function of the weight at index i in
 *  layer l. dx is added as an offset to the current value of the weight. */
//______________________________________________________________________________
template <typename Architecture>
auto evaluate_net_bias(TDeepNet<Architecture> &net, typename Architecture::Tensor_t &X,
                       const typename Architecture::Matrix_t &Y, const typename Architecture::Matrix_t &W, size_t l,
                       size_t i, size_t k, typename Architecture::Scalar_t dx) -> typename Architecture::Scalar_t
{
    using Scalar_t = typename Architecture::Scalar_t;
    //using Matrix_t = typename Architecture::Matrix_t;

    auto & netB = net.GetLayerAt(l)->GetBiasesAt(k);

    netB(i,0) += dx;

    //std::cout << "------ net biases layer " << l << " k  " << k <<  " element " << i << " dx = " <<  dx << std::endl;
    //Architecture::PrintTensor( netB, "netB bias tensor");
    //typename Architecture::Tensor_t temp(netB);


    Scalar_t res = net.Loss(X, Y, W);

    netB(i,0) -= dx;


    return res;
}



template <typename Architecture>
auto testConvBackwardPass(size_t batchSize, size_t imgDepth, size_t imgHeight, size_t imgWidth, size_t batchDepth,
                          size_t batchHeight, size_t batchWidth, typename Architecture::Scalar_t dx, ETestType testType) -> bool
{
   using Matrix_t = typename Architecture::Matrix_t;
   using Tensor_t = typename Architecture::Tensor_t;
   using Net_t = TDeepNet<Architecture>;
   using Scalar_t = typename Architecture::Scalar_t;

   Net_t convNet(batchSize, imgDepth, imgHeight, imgWidth, batchDepth, batchHeight, batchWidth,
                 ELossFunction::kMeanSquaredError, EInitialization:: kGlorotUniform);
   // tyoe of network
   if (testType == kLinearNet)
      constructLinearConvNet(convNet);
   else
      constructConvNet(convNet);

   Architecture::SetRandomSeed(111); // use fixed seed
   convNet.Initialize();

   std::cout << "test backward on this network " << std::endl;
   convNet.Print();

//    auto & w0 = convNet.GetLayerAt(0)->GetWeights();
// #ifdef DEBUG
//    std::cout << "Netwrok weights for Layer 0  " << std::endl;
//    std::cout << "weight depth = " << w0.size() << std::endl;
//    for (size_t i = 0; i < w0.size();   ++i)
//       TMVA_DNN_PrintTCpuMatrix(w0[i],"weight-layer0");
// #endif

   //typename Architecture::Tensor_t X(batchSize, imgDepth, imgHeight * imgWidth, Architecture::GetTensorLayout() );
   auto X =  Architecture::CreateTensor( batchSize, imgDepth, imgHeight , imgWidth);
   randomBatch(X);

   Architecture::PrintTensor( X, "X input tensor");

   Matrix_t Y(batchSize, convNet.GetOutputWidth());
   Matrix_t W(batchSize, 1);   // this are the data weights
   randomMatrix(Y);
   //randomMatrix(W);
   // for the moment use weights equal to 1
   fillMatrix(W, 1.0);

   for (size_t l = 0; l < convNet.GetLayers().size(); l++) {
      auto & theLayer = *(convNet.GetLayers()[l]);
      auto & vW = theLayer.GetWeights();
      if (vW.size() > 0) {
         TString tname = TString::Format("weight-tensor-layer-%d",(int)l);
         Tensor_t tW( vW[0] );
         Architecture::PrintTensor( tW ,std::string(tname));
      }
      auto & vB = theLayer.GetBiases();
      if (vB.size() > 0) {
         TString tname = TString::Format("bias-tensor-layer-%d",(int)l);
         Architecture::InitializeIdentity(vB[0]);
         Tensor_t tB( vB[0] );
         Architecture::PrintTensor( tB ,std::string(tname));
      }
   }

   auto & lLayer = *(convNet.GetLayers().back());

   Architecture::PrintTensor(Tensor_t(lLayer.GetWeights()[0]),"weights-before-fw");
   std::cout << "Do Forward Pass " << std::endl;
   convNet.Forward(X);

   Architecture::PrintTensor(Tensor_t(lLayer.GetWeights()[0]),"weights-after-fw");


   // print layer derivatives
#ifdef DEBUG
   using ConvLayer_t = TConvLayer<Architecture>;
   auto convLayer = dynamic_cast<ConvLayer_t*>(convNet.GetLayerAt(0) );
   if (convLayer) {
      auto & df = convLayer->GetDerivatives();
      std::cout << "Derivatives - size " << df.GetFirstSize() << std::endl;
      for (size_t ii=0; ii< df.GetFirstSize(); ++ii)
         TMVA_DNN_PrintTCpuMatrix(df[ii],"Derivatives");
   }
#endif

   //if (testType == kRndmActNet)  return true;

   std::cout << "Do Backward Pass " << std::endl;
   Architecture::PrintTensor(X,"input");
   convNet.Backward(X, Y, W);

   //Architecture::PrintTensor(Tensor_t(lLayer.GetWeights()[0]),"weights-after-bw");


   // now compare derivatives using finite differences and compare the result
   Scalar_t maximum_error = 0.0;

   for (size_t l = convNet.GetDepth()-1; (int) l >= 0; l--) {
      std::cout << "\n\n************************************* \n";
      std::cout << "\tTesting weight gradients:      layer: " << l << " / " << convNet.GetDepth();
      std::cout << std::flush;
      auto & layer = *(convNet.GetLayerAt(l));
      std::cout << std::endl;
      layer.Print();
      std::cout << "************************************* \n\n";

      //Architecture::PrintTensor(layer.GetOutput(),"output tensor");
      //Architecture::PrintTensor(layer.GetActivationGradients(),"activation gradient");

      auto &gw = layer.GetWeightGradients();
      auto &gb = layer.GetBiasGradients();

      if (gw.size() > 0) {
#ifdef DEBUG
         std::cout << "Weight gradient from back-propagation - vector size is " << gw.size()  << std::endl;
         Architecture::PrintTensor(Tensor_t(layer.GetWeights()[0]),"weights");
         Architecture::PrintTensor(Tensor_t(layer.GetWeightGradients()[0]),"weight gradients");
         Architecture::PrintTensor(Tensor_t(layer.GetBiasGradients()[0]),"bias gradients");
         std::cout << std::endl;
#endif
         // if (gw[0].GetNoElements() < 100 ) {
            //    gw[0].Print();
         //    // }
         // else
         //    std::cout << "BP Weight Gradient ( " << gw[0].GetNrows() << " x " << gw[0].GetNcols() << " ) , ...... skip printing (too many elements ) " << std::endl;
         // }
      }
      else {
         std::cout << "Layer " << l << " has no weights " << std::endl;
         continue;
      }
      std::cout << "Layer " << l << " :  output  D x H x W " << layer.GetDepth() << "  " << layer.GetHeight() << "  "
                << layer.GetWidth() << "\t input D x H x W " << layer.GetInputDepth() << "  " << layer.GetInputHeight()
                << "  " << layer.GetInputWidth() << std::endl;

      // print output and activation gradients

#ifdef DEBUG
      auto &outL = layer.GetOutput();
      auto & actGrad = layer.GetActivationGradients();
      if (Architecture::GetTensorLayout() == TMVA::Experimental::MemoryLayout::ColumnMajor) {
         std::cout << "layer output size " << outL.GetFirstSize() << std::endl;
         if (outL.GetFirstSize() > 0) {
            if (outL.GetSize() < 100) {
               outL.At(0).GetMatrix().Print();
            } else
               std::cout << "Layer Output ( " << outL.GetHSize() << " x " << outL.GetWSize()
                         << " ) , ...... skip printing (too many elements ) " << std::endl;
         }
         if (actGrad.GetFirstSize() > 0) {
            std::cout << "Activation gradient from back-propagation  - vector size is " << actGrad.GetFirstSize()
                      << std::endl;
            if (actGrad.GetSize() < 100) {
               for (size_t ii = 0; ii < actGrad.GetFirstSize(); ++ii)
                  actGrad.At(ii).GetMatrix().Print();
            } else
               std::cout << "Activation Gradient ( " << actGrad.GetHSize() << " x " << actGrad.GetWSize()
                         << " ) , ...... skip printing (too many elements ) " << std::endl;
         }
      } else {
         Architecture::PrintTensor(outL, "layer output");
         Architecture::PrintTensor(actGrad,"activation gradients");
      }
#endif

      std::cout << "Evaluate the Derivatives with Finite difference and compare with BP for Layer " << l << std::endl;
      int nerrorsW = 0;
      int nerrorsB = 0;
      int ngoodW = 0;
      int ngoodB = 0;
#ifdef DEBUG
      int ngoodPrint = 10000;
#else
      int ngoodPrint = 10;
#endif
      //  conv layer weights is a matrix
      // for (size_t k = 0; k <  gw.GetFirstSize() ; ++k) {
      // for (size_t i = 0; i < layer.GetWidth(); i++) {
      size_t k = 0;
      Matrix_t & gwm = gw[k];
      Matrix_t & gbm = gb[k];
      // size_t n1 = (gwm.GetNDim() == 2 ) ? gwm.GetNrows() : gwm.GetCSize();
      // size_t n1b =  (gwb.GetNDim() == 2 ) ? gwb.GetNrows() : gwb.GetCSize();
      // assert(n1 == n1b);
      // size_t n2 = (gwm.GetNDim() == 2 ) ? gwm.GetNcols() : gwm.GetHSize()*gwm.GetWSize();
      size_t n1 =  gwm.GetNrows();
      size_t n2 = gwm.GetNcols();
      size_t n1b =  gbm.GetNrows();
      //assert(n1b <= n1);
      //Matrix_t & wm = layer.GetWeightsAt(k);
      // make matrices in case they are not
      for (size_t i = 0; i < n1; i++) {
         // for (size_t j = 0; j < layer.GetInputWidth(); j++) {
         for (size_t j = 0; j < n2; j++) {
            auto f = [&convNet, &X, &Y, &W, l, i, j, k](Scalar_t x) {
               return evaluate_net_weight(convNet, X, Y, W, l, i, j, k, x);
            };
            Scalar_t dy = finiteDifference(f, dx) / (2.0 * dx);
            Scalar_t dy_ref = gwm(i, j);
            // Compute the relative error if dy != 0.
            Scalar_t error;
            if (std::fabs(dy_ref) > 1e-15) {
               error = std::fabs((dy - dy_ref) / dy_ref);
            } else {
               error = std::fabs(dy - dy_ref);
            }
            maximum_error = std::max(error, maximum_error);
            if (error > 1.E-3) {
               std::cout << k << " weight - " << i << " , " << j << " : " << dy << " from BP " << dy_ref << "   " << error
                         << " ERROR " << std::endl;
               nerrorsW++;
            } else {
               if (ngoodW < ngoodPrint)
                  std::cout << k << " weight - " << i << " , " << j << " : " << dy << " from BP " << dy_ref << "   " << error
                            << std::endl;
               ngoodW++;
            }
            if (nerrorsW > 10) {
               std::cout << "Reached error limit skip..." << std::endl;
               break;
            }
         }
      }

      for (size_t i = 0; i < n1b; i++)
      {
         // evaluate the bias gradients
         //std::cout << " i " << i << "  " << n1b << "  " << n1 << std::endl;

         auto fb = [&convNet, &X, &Y, &W, l, i, k](Scalar_t x) {
            return evaluate_net_bias(convNet, X, Y, W, l, i, k, x);
         };
         Scalar_t dy = finiteDifference(fb, dx) / (2.0 * dx);
         Scalar_t dy_ref = gbm(i, 0);
         // Compute the relative error if dy != 0.
         Scalar_t error;
         if (std::fabs(dy_ref) > 1e-15)
         {
            error = std::fabs((dy - dy_ref) / dy_ref);
         }
         else
         {
            error = std::fabs(dy - dy_ref);
         }
         maximum_error = std::max(error, maximum_error);
         if (error > 1.E-3)
         {
            std::cout << k << " bias - " << i << " : " << dy << " from BP " << dy_ref << "   " << error
                      << " ERROR " << std::endl;
            nerrorsB++;
         }
         else
         {
            if (ngoodB < ngoodPrint)
               std::cout << k << " bias - " << i << " : " << dy << " from BP " << dy_ref << "   " << error
                         << std::endl;
            ngoodB++;
         }
         if (nerrorsB > 10)
         {
            std::cout << "Reached error limit skip..." << std::endl;
            break;
         }


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
