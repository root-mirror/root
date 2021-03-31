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
 *      Surya S Dwivedi        <surya2191997@gmail.com>      - UT, Austin         *
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
using namespace TMVA::DNN::CNN_3D;

enum ETestType { kLinearNet, kRndmActNet };

/** Testing the image to column function. Check wheether the matrix A after
 *  the transformation will be equal to the matrix B. */




//______________________________________________________________________________
template <typename Architecture>
auto testIm2col3D(typename Architecture::Matrix_t &A, typename Architecture::Matrix_t &B, size_t imgHeight,
                size_t imgWidth, size_t imgDepth, size_t fltHeight, size_t fltWidth, size_t fltDepth, size_t strideRows, 
                size_t strideCols, size_t strideDepth, size_t zeroPaddingHeight, size_t zeroPaddingWidth, size_t zeroPaddingDepth) -> bool
{

   size_t m, n;
   m = B.GetNrows();
   n = B.GetNcols();

   typename Architecture::Matrix_t ATr(m, n);
   Architecture::Im2col3D(ATr, A, imgHeight, imgWidth, imgDepth,  fltHeight, fltWidth, fltDepth, strideRows, strideCols, strideDepth, zeroPaddingHeight,
                        zeroPaddingWidth, zeroPaddingDepth);

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
                       size_t filterHeight, size_t filterWidth, size_t numFilters, size_t input4D) -> bool
{

   size_t m, n;
   m = B.GetNrows();
   n = B.GetNcols();

   typename Architecture::Matrix_t ARot(m, n);
   Architecture::RotateWeights3D(ARot, A, filterDepth, filterHeight, filterWidth, numFilters, input4D);

   for (size_t i = 0; i < m; i++) {
      for (size_t j = 0; j < n; j++) {
         if (ARot(i, j) != B(i, j)) {
            return false;
         }
      }
   }

   return true;
}


template <typename Architecture>
auto testMaxPoolLayerForward(const typename Architecture::Tensor_t &input,
                          const typename Architecture::Tensor_t &expectedOutput,
                          typename Architecture::Tensor_t &index,
                          size_t inputHeight, size_t inputWidth, size_t inputDepth, size_t fltHeight,
                          size_t fltWidth, size_t fltDepth, size_t strideRows, size_t strideCols, size_t strideDepth) -> bool
{
    size_t nRows = expectedOutput.GetHSize();
    size_t nCols = expectedOutput.GetWSize();
    size_t batchSize = 1;
    size_t input4D  = input.GetHSize();
    // size_t output4D = input4D;


    typename Architecture::Tensor_t computedOutput( batchSize, nRows, nCols);
    typename Architecture::Tensor_t computedDerivatives(batchSize, nRows, nCols);


    size_t height = (inputHeight - fltHeight ) / strideRows + 1;
    size_t width =  (inputWidth - fltWidth ) / strideCols + 1;
    size_t depth =  (inputDepth - fltDepth ) / strideDepth + 1;
    size_t nLocalViews = height * width * depth;
    size_t nLocalViewPixels = input4D * fltHeight * fltWidth * fltDepth;

    typename Architecture::Tensor_t forwardMatrices(1 , nLocalViews, nLocalViewPixels);

    TDescriptors * convDescriptors = nullptr;
    TWorkspace   * convWorkspace   = nullptr;

    // TConvLayer<Architecture> *layer = nullptr;
    // Architecture::InitializeConvDescriptors(convDescriptors, layer);
    // Architecture::InitializeConvWorkspace(convWorkspace, convDescriptors, params, layer);

    Architecture::Downsample3D(computedOutput, index, input,
                              (typename Architecture::PoolingDescriptors_t &) *convDescriptors,
                              (typename Architecture::PoolingWorkspace_t &) *convWorkspace,
                              inputHeight, inputWidth, inputDepth, fltHeight, fltWidth, fltDepth, 
                              strideRows, strideCols, strideDepth);



    Architecture::PrintTensor(computedOutput,"computed output tensor");
    Architecture::PrintTensor(expectedOutput,"expected output tensor");
    
    bool val = true;
    for (size_t slice = 0; slice < nRows; slice++) {
        for (size_t localView = 0; localView < nCols; localView++) {
            if (expectedOutput(0, slice, localView) != computedOutput(0, slice, localView)) {
               val = false;             
            }
        }
    }
    if (!val) 
      return false;
    return true;
}


template <typename Architecture>
auto testConvLayerForward(const typename Architecture::Tensor_t &input,
                          const typename Architecture::Tensor_t &expectedOutput,
                          const typename Architecture::Matrix_t &weights, const typename Architecture::Matrix_t &biases,
                          size_t inputHeight, size_t inputWidth, size_t inputDepth, size_t fltHeight,
                          size_t fltWidth, size_t fltDepth, size_t strideRows, size_t strideCols, size_t strideDepth,
                          size_t zeroPaddingHeight, size_t zeroPaddingWidth, size_t zeroPaddingDepth) -> bool
{
    size_t nRows = expectedOutput.GetHSize();
    size_t nCols = expectedOutput.GetWSize();
    size_t batchSize = 1;
    size_t input4D  = input.GetHSize();
    size_t output4D = weights.GetNrows();


    typename Architecture::Tensor_t computedOutput( batchSize, nRows, nCols);
    typename Architecture::Tensor_t computedDerivatives(batchSize, nRows, nCols);


   TConv3DParams params(1, inputDepth, inputHeight, inputWidth, input4D, output4D,
                       fltHeight, fltWidth, fltDepth, strideRows, strideCols, strideCols, 
                       zeroPaddingHeight, zeroPaddingWidth, zeroPaddingDepth);


    size_t height = (inputHeight - fltHeight + 2 * zeroPaddingHeight) / strideRows + 1;
    size_t width =  (inputWidth - fltWidth + 2 * zeroPaddingWidth) / strideCols + 1;
    size_t depth =  (inputDepth - fltDepth + 2 * zeroPaddingDepth) / strideDepth + 1;
    size_t nLocalViews = height * width * depth;
    size_t nLocalViewPixels = input4D * fltHeight * fltWidth * fltDepth;

    typename Architecture::Tensor_t forwardMatrices (1 , nLocalViews, nLocalViewPixels);

    TDescriptors * convDescriptors = nullptr;
    TWorkspace   * convWorkspace   = nullptr;

    // TConvLayer<Architecture> *layer = nullptr;
    // Architecture::InitializeConvDescriptors(convDescriptors, layer);
    // Architecture::InitializeConvWorkspace(convWorkspace, convDescriptors, params, layer);

    Architecture::Conv3DLayerForward(computedOutput, computedDerivatives, input, weights, biases, params,
                                   EActivationFunction::kIdentity, forwardMatrices,
                                   (typename Architecture::ConvDescriptors_t &) *convDescriptors,
                                   (typename Architecture::ConvWorkspace_t &) *convWorkspace);



    Architecture::PrintTensor(computedOutput,"computed output tensor");
    Architecture::PrintTensor(expectedOutput,"expected output tensor");
    

    bool val = true;
    for (size_t slice = 0; slice < nRows; slice++) {
        for (size_t localView = 0; localView < nCols; localView++) {
            if (expectedOutput(0, slice, localView) != computedOutput(0, slice, localView)) {
               val = false;
            }
        }
    }
    
    if (!val) 
      return false;
    return true;
}


template<typename Architecture>
void test_func()
{
   using Matrix_t = typename Architecture::Matrix_t;
   using Tensor_t = typename Architecture::Tensor_t;


   // instantiating net
   double img[][27] = {
           {1, 2 , 3, 4, 5, 6, 7, 8, 9,
            1, 2, 3, 4, 5, 6, 7, 8, 9,
             1, 2, 3, 4, 5, 6, 7, 8 ,9},

           { 0, 0, 0, 0, 0, 0, 0, 0 , 0,
            0, 0, 0, 0, 0, 0, 0, 0 , 0,
             0, 0, 0, 0, 0, 0, 0, 0 , 0}
   };


   // num_rows = num_filter, each of row of size, filter_size*input4D = 8*2 = 16
   double weights[][16] = {
           {1, 1, 1, 1,
            1,  1, 1, 1,
            // repeat same thing
             1, 1, 1, 1,
            1,  1, 1, 1,},

           {0, 0, 0, 0,
             // repeat same thing
           0, 0, 0, 0}
   };

   double biases[][1] = {
           {0},

           {0}
   };

   double expected_cnn[][8] = {

           {24, 32, 48, 56, 24, 32, 48, 56},

           {0, 0, 0, 0, 0, 0, 0, 0}
    };


   size_t input4D = 2;
   size_t imgDepth = 3;
   size_t imgHeight = 3;
   size_t imgWidth = 3;
   size_t numberFilters = 2;
   size_t output4D = 2;
   size_t fltHeight = 2;
   size_t fltWidth = 2;
   size_t fltDepth = 2;
   size_t strideRows = 1;
   size_t strideCols = 1;
   size_t strideDepth = 1;
   size_t zeroPaddingHeight = 0;
   size_t zeroPaddingWidth = 0;
   size_t zeroPaddingDepth = 0;

   Matrix_t inputEvent(input4D, imgHeight * imgWidth * imgDepth);

   for (size_t i = 0; i < input4D; i++) {
      for (size_t j = 0; j < imgHeight * imgWidth * imgDepth; j++) {
         inputEvent(i, j) = img[i][j];
      }
   }
   Tensor_t input(inputEvent, 3);

   Matrix_t weightsMatrix(numberFilters, fltHeight * fltWidth * fltDepth * input4D);
   Matrix_t biasesMatrix(numberFilters, 1);
   for (size_t i = 0; i < numberFilters; i++) {
       for (size_t j = 0; j < fltHeight * fltWidth * fltDepth * input4D; j++){
           weightsMatrix(i, j) = weights[i][j];
       }
       biasesMatrix(i, 0) = biases[i][0];
   }

    size_t height = (imgHeight - fltHeight + 2 * zeroPaddingHeight) / strideRows + 1;
    size_t width =  (imgWidth - fltWidth + 2 * zeroPaddingWidth) / strideCols + 1;
    size_t depth =  (imgDepth - fltDepth + 2 * zeroPaddingDepth) / strideDepth + 1;

   Matrix_t output_cnn(numberFilters, height * width * depth);

   for (size_t i = 0; i < numberFilters; i++) {
      for (size_t j = 0; j < height * width * depth; j++) {
         output_cnn(i, j) = expected_cnn[i][j];
      }
   }
   Tensor_t expectedOutput_cnn (output_cnn, 3);
   

   // conv3d forward pass
   bool status = false; 
   std::cout << "************Testing Conv 3D forward pass************" << std::endl;
   status = testConvLayerForward<Architecture>(input, expectedOutput_cnn, weightsMatrix, biasesMatrix, imgHeight,
                                                    imgWidth, imgDepth, fltHeight, fltWidth, fltDepth, strideRows,
                                                    strideCols, strideDepth, zeroPaddingHeight, zeroPaddingWidth, zeroPaddingDepth);
   if(status)
   		std::cout << "Forward test passed" << std::endl;
   else
   		std::cout << "Forward test failed" << std::endl ;

    std::cout<< std::endl << std::endl;




   // max_pool 3d forward pass
   double expected_maxpool[][8] = {

           {5, 6, 8, 9, 5, 6, 8, 9},

           {0, 0, 0, 0, 0, 0, 0, 0}

    };

    Matrix_t output_maxpool(input4D, height * width * depth);

   for (size_t i = 0; i < input4D; i++) {
      for (size_t j = 0; j < height * width * depth; j++) {
         output_maxpool(i, j) = expected_maxpool[i][j];
      }
   }
  
  Tensor_t expectedOutput_maxpool (output_maxpool, 3);

  Tensor_t index(1, input4D, height * width * depth);  
  std::cout << "************Testing MaxPool 3D forward pass************" << std::endl;
  status = testMaxPoolLayerForward<Architecture>(input, expectedOutput_maxpool, index, imgHeight,
                                                    imgWidth, imgDepth, fltHeight, fltWidth, fltDepth, strideRows,
                                                    strideCols, strideDepth);

  if(status)
    std::cout << "Forward test passed" << std::endl;
  else
    std::cout << "Forward pass test failed" << std::endl ;

  std::cout << std::endl << std::endl;




  // CNN 3D backward pass

  std::cout << "************Testing Conv 3D backward pass************" << std::endl;


  // wrt to input, dl/dI
  Tensor_t activationGradientsBackward(input);
  Matrix_t weightGradients(weightsMatrix);
  Matrix_t biasGradients(biasesMatrix);

  // o/p without activation, basically expected_o/p
  Tensor_t activationsBackward(expectedOutput_cnn);

  Matrix_t tmp(input4D, height * width * depth);
  for (size_t i = 0; i < output4D; i++) {
    for (size_t j = 0; j < height * width * depth; j++) {
         tmp(i, j) = 1;
    }
  }

  for (size_t i = 0; i < input4D; i++) {
    for (size_t j = 0; j < imgHeight * imgWidth * imgDepth; j++) {
         activationGradientsBackward(0, i, j) = 0;
    }
  }
  
  TDescriptors * convDescriptors = nullptr;
  TWorkspace   * convWorkspace   = nullptr;

  // need to provide this, all 1's for now
  Tensor_t activationGradients(tmp, 3);

  TConv3DParams params(1, imgHeight, imgWidth, imgDepth, input4D, output4D,
                       fltHeight, fltWidth, fltDepth, strideRows, strideCols, strideCols, 
                       zeroPaddingHeight, zeroPaddingWidth, zeroPaddingDepth);

  Architecture::Conv3DLayerBackward(activationGradientsBackward,
                                    weightGradients, biasGradients,
                                    input,
                                    activationGradients,
                                    weightsMatrix,
                                    activationsBackward,
                                    expectedOutput_cnn,
                                    EActivationFunction::kIdentity,
                                    (typename Architecture::ConvDescriptors_t &) *convDescriptors,
                                    (typename Architecture::ConvWorkspace_t &) *convWorkspace,
                                    params);

    Architecture::PrintTensor(activationGradientsBackward,"activationGradientsBackward");
    Architecture::PrintTensor(weightGradients,"weightGradients");
    Architecture::PrintTensor(biasGradients,"biasGradients");
    std::cout << std::endl << std::endl;





    // MaxPool 3D backward pass

    std::cout << "************Testing MaxPool3D 3D backward pass************" << std::endl;

    for (size_t i = 0; i < input4D; i++) {
      for (size_t j = 0; j < imgHeight * imgWidth * imgDepth; j++) {
           activationGradientsBackward(0, i, j) = 0;
      }
    }

    height = (imgHeight - fltHeight ) / strideRows + 1;
    width =  (imgWidth  - fltWidth ) / strideCols + 1;
    depth =  (imgDepth  - fltDepth ) / strideDepth + 1;

    Architecture::MaxPoolLayer3DBackward(activationGradientsBackward, activationGradients, index,
                              input, expectedOutput_maxpool, 
                             (typename Architecture::PoolingDescriptors_t &) *convDescriptors,
                              (typename Architecture::PoolingWorkspace_t &) *convWorkspace,
                              imgHeight, imgWidth, fltHeight, fltWidth,  strideRows,
                              strideCols, height*width*depth);


    Architecture::PrintTensor(activationGradientsBackward,"activationGradientsBackward");


  // MaxPool3D Backward Pass test
}


#endif

