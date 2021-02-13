// @(#)root/tmva/tmva/cnn:$Id$
// Author: Manos Stergiadis

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  :                                                                       *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Testing Downsample method on a CPU architecture                           *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Manos Stergiadis       <em.stergiadis@gmail.com>  - CERN, Switzerland     *
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

////////////////////////////////////////////////////////////////////
// Testing the Convolutional Layer                                //
////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include <vector>

#include "TestConvNet.h"

using namespace TMVA::DNN;
using namespace TMVA::DNN::CNN;
using namespace TMVA::Experimental;

inline bool isInteger(double x)
{
   return x == floor(x);
}

size_t calculateDimension(size_t imgDim, size_t fltDim, size_t padding, size_t stride)
{
   double dimension = ((imgDim - fltDim + 2 * padding) / stride) + 1;
   if (!isInteger(dimension)) {
      std::cout << "Not compatible hyper parameters" << std::endl;
      std::exit(EXIT_FAILURE);
   }

   return (size_t)dimension;
}

template<typename AFloat>
bool almostEqual(AFloat expected, AFloat computed, double epsilon = 0.0001) {
    return abs(computed - expected) < epsilon;
}

/*************************************************************************
 * Test 1: Forward Propagation
 *  batch size = 1
 *  image depth = 2, image height = 4, image width = 4,
 *  num frames = 3, filter height = 2, filter width = 2,
 *  stride rows = 2, stride cols = 2,
 *  zero-padding height = 0, zero-padding width = 0,
 *************************************************************************/
template<typename Architecture>
bool testForward1()
{
   using Matrix_t = typename Architecture::Matrix_t;
   using Tensor_t = typename Architecture::Tensor_t;

   double img[][16] = {
           {166, 212, 213, 150,
            114, 119, 109, 115,
             88, 144, 227, 208,
            208, 235,  57,  58},

           { 57,  255, 184, 162,
            204,  220,  11, 192,
            183,  174,   2, 153,
            184,  175,  10,  55}
   };

   double weights[][8] = {
           {2.0,  3.0,  0.5, -1.5,
            1.0,  1.5, -2.0, -3.0},

           {-0.5,  1.0,  2.5, -1.0,
             2.0,  1.5, -0.5,  1.0},

           {-1.0, -2.0, 1.5, 0.5,
             2.0, -1.5, 0.5, 1.0}
   };

   double biases[][1] = {
           {45},

           {60},

           {12}
   };

   double expected[][9] = {

           {263.0, 1062.0,  632.0,
            104.0,  224.0,  245.5,
            -44.5,  843.0, 1111.0},

           { 969.5, 1042.5, 1058.5,
            1018.5,  614.0,  942.0,
            1155.0, 1019.0,  522.5},

           {-294.0, -38.0,   42.5,
             207.5, 517.0,    5.5,
             437.5, 237.5, -682.0}
    };



   size_t imgDepth = 2;
   size_t imgHeight = 4;
   size_t imgWidth = 4;
   size_t numberFilters = 3;
   size_t fltHeight = 2;
   size_t fltWidth = 2;
   size_t strideRows = 1;
   size_t strideCols = 1;
   size_t zeroPaddingHeight = 0;
   size_t zeroPaddingWidth = 0;

   Matrix_t inputEvent(imgDepth, imgHeight * imgWidth);

   for (size_t i = 0; i < imgDepth; i++) {
      for (size_t j = 0; j < imgHeight * imgWidth; j++) {
         inputEvent(i, j) = img[i][j];
      }
   }
   Tensor_t input(inputEvent, 4);

   Matrix_t weightsMatrix(numberFilters, fltHeight * fltWidth * imgDepth);
   Matrix_t biasesMatrix(numberFilters, 1);
   for (size_t i = 0; i < numberFilters; i++) {
       for (size_t j = 0; j < fltHeight * fltWidth * imgDepth; j++){
           weightsMatrix(i, j) = weights[i][j];
       }
       biasesMatrix(i, 0) = biases[i][0];
   }

   size_t height = calculateDimension(imgHeight, fltHeight, zeroPaddingHeight, strideRows);
   size_t width = calculateDimension(imgWidth, fltWidth, zeroPaddingWidth, strideCols);

   Matrix_t outputEvent(numberFilters, height * width);

   for (size_t i = 0; i < numberFilters; i++) {
      for (size_t j = 0; j < height * width; j++) {
         outputEvent(i, j) = expected[i][j];
      }
   }
   Tensor_t expectedOutput (outputEvent, 3);
   std::cout<< "Calling TestNet.h" << std::endl;

   bool status = testConvLayerForward<Architecture>(input, expectedOutput, weightsMatrix, biasesMatrix, imgHeight,
                                                    imgWidth, imgDepth, fltHeight, fltWidth, numberFilters, strideRows,
                                                    strideCols, zeroPaddingHeight, zeroPaddingWidth);

   return status;
}


template<typename Architecture>
bool testForward1_cudnn()
{
   double img[][16] = {
           {166, 212, 213, 150,
            114, 119, 109, 115,
             88, 144, 227, 208,
            208, 235,  57,  58},

           { 57,  255, 184, 162,
            204,  220,  11, 192,
            183,  174,   2, 153,
            184,  175,  10,  55}
   };

   double weights[][8] = {
           {2.0,  3.0,  0.5, -1.5,
            1.0,  1.5, -2.0, -3.0},

           {-0.5,  1.0,  2.5, -1.0,
             2.0,  1.5, -0.5,  1.0},

           {-1.0, -2.0, 1.5, 0.5,
             2.0, -1.5, 0.5, 1.0}
   };

   double biases[][9] = {
           {45, 45, 45,
            45, 45, 45,
            45, 45, 45},

           {60, 60, 60,
            60, 60, 60,
            60, 60, 60},

           {12, 12, 12,
            12, 12, 12,
            12, 12, 12}
   };

   double expected[][9] = {

           {263.0, 1062.0,  632.0,
            104.0,  224.0,  245.5,
            -44.5,  843.0, 1111.0},

           { 969.5, 1042.5, 1058.5,
            1018.5,  614.0,  942.0,
            1155.0, 1019.0,  522.5},

           {-294.0, -38.0,   42.5,
             207.5, 517.0,    5.5,
             437.5, 237.5, -682.0}
    };


   size_t imgDepth = 2;
   size_t imgHeight = 4;
   size_t imgWidth = 4;
   size_t numberFilters = 3;
   size_t fltHeight = 2;
   size_t fltWidth = 2;
   size_t strideRows = 1;
   size_t strideCols = 1;
   size_t zeroPaddingHeight = 0;
   size_t zeroPaddingWidth = 0;

   using Matrix_t = typename Architecture::Matrix_t;
   using Tensor_t = typename Architecture::Tensor_t;
   using HostBuffer_t = typename Architecture::HostBuffer_t;

   std::vector<size_t> inputShape {1, imgDepth, imgHeight, imgWidth};
   HostBuffer_t    input_hostbuffer(imgDepth * imgHeight * imgWidth);
   for (size_t i = 0; i < imgDepth; i++) {
      for (size_t j = 0; j < imgHeight * imgWidth; j++) {
         input_hostbuffer[i*imgHeight * imgWidth + j] = img[i][j];
      }
   }
   Tensor_t input(input_hostbuffer, inputShape, MemoryLayout::RowMajor, 0, 0);

   std::vector<size_t> weightShape {numberFilters, imgDepth, fltHeight, fltWidth};
   HostBuffer_t weight_hostbuffer(numberFilters * fltHeight * fltWidth * imgDepth);
   for (size_t i = 0; i < numberFilters; i++) {
       for (size_t j = 0; j < fltHeight * fltWidth * imgDepth; j++){
           weight_hostbuffer[i*fltHeight * fltWidth * imgDepth + j] = weights[i][j];
       }
   }

   Matrix_t weightsTensor(weight_hostbuffer, weightShape, MemoryLayout::RowMajor, 0, 0);

   size_t height = calculateDimension(imgHeight, fltHeight, zeroPaddingHeight, strideRows);
   size_t width = calculateDimension(imgWidth, fltWidth, zeroPaddingWidth, strideCols);

   std::vector<size_t> biasesShape {1, numberFilters, height, width};
   HostBuffer_t biases_hostbuffer(numberFilters * height * width);

   for (size_t i = 0; i < numberFilters; i++) {
      for (size_t j = 0; j < height * width; j++) {
         biases_hostbuffer[i * height * width + j] = biases[i][j];
      }
   }
   Tensor_t biasesTensor(biases_hostbuffer, biasesShape, MemoryLayout::RowMajor, 0, 0);

   Tensor_t computedDerivatives(1, 3, height, width, MemoryLayout::RowMajor,0,0);
   Tensor_t computedOutput (1, 3, height, width, MemoryLayout::RowMajor,0,0);

   TConvParams params(1, imgDepth, imgHeight, imgWidth, numberFilters, fltHeight, fltWidth, strideRows,
                      strideCols, zeroPaddingHeight, zeroPaddingWidth);

   Tensor_t forwardMatrix;

   EActivationFunction AFunct = EActivationFunction::kIdentity;
   //EActivationFunction AFunct = EActivationFunction::kRelu;
   TConvLayer<Architecture> convLayer (1, imgDepth, imgHeight, imgWidth, numberFilters,
                                           EInitialization::kIdentity, fltHeight, fltWidth,
                                           strideRows, strideCols, zeroPaddingHeight, zeroPaddingWidth,
                                           0.0, AFunct, ERegularization::kNone, 0.0);

   auto& convDescriptors = static_cast<typename Architecture::ConvDescriptors_t &> (*convLayer.GetDescriptors());
   auto& convWorkspace   = static_cast<typename Architecture::ConvWorkspace_t &> (*convLayer.GetWorkspace());

   Architecture::ConvLayerForward(computedOutput, computedDerivatives, input, weightsTensor, biasesTensor, params,
                                       AFunct, forwardMatrix, convDescriptors, convWorkspace);

   HostBuffer_t expectedOutput_buffer(numberFilters * height * width);
   for (size_t i = 0; i < numberFilters; i++) {
      for (size_t j = 0; j < height * width; j++) {
         expectedOutput_buffer[i * height * width + j] = expected[i][j];
      }
   }
   Architecture::PrintTensor(computedOutput ,"Convolution output: ");

   return computedOutput.isEqual(expectedOutput_buffer, expectedOutput_buffer.GetSize());
}

/*************************************************************************
* Test 1: Backward Propagation
*  batch size = 1
*  image depth = 2, image height = 5, image width = 5,
*  num frames = 2, filter height = 3, filter width = 3,
*  stride rows = 1, stride cols = 1,
*  zero-padding height = 0, zero-padding width = 0,
*************************************************************************/
template<typename Architecture>
bool testBackward1()
{
    using Matrix_t = typename Architecture::Matrix_t;
    using Tensor_t = typename Architecture::Tensor_t;

    size_t imgDepth = 2;
    size_t imgHeight = 5;
    size_t imgWidth = 5;
    size_t numberFilters = 2;
    size_t fltHeight = 3;
    size_t fltWidth = 3;
    size_t strideRows = 1;
    size_t strideCols = 1;
    size_t zeroPaddingHeight = 0;
    size_t zeroPaddingWidth = 0;

    size_t height = calculateDimension(imgHeight, fltHeight, zeroPaddingHeight, strideRows);
    size_t width = calculateDimension(imgWidth, fltWidth, zeroPaddingWidth, strideCols);
    size_t nLocalViews =  height * width;
    size_t batchSize = 1;

    double grad[][9] = {
            {0, 1.37, 0, 0, 0, 0, 0, -0.90, 0},
            {0, -0.37, 0, 0, 0, -0.25, 0.26, 0, 0}
    };

    Matrix_t gradEvent(numberFilters, height * width);

    for (size_t i = 0; i < numberFilters; i++) {
        for (size_t j = 0; j < height * width; j++) {
            gradEvent(i, j) = grad[i][j];
        }
    }
    Tensor_t activationGradients(gradEvent, 3);

    double derivatives[][9] = {
            {1, 1, 1, 1 , 1, 1, 1, 1, 1},
            {1, 1, 1, 1 , 1, 1, 1, 1, 1}
    };

    Matrix_t dfEvent(numberFilters, height * width);

    for (size_t i = 0; i < numberFilters; i++) {
        for (size_t j = 0; j < height * width; j++) {
            dfEvent(i, j) = derivatives[i][j];
        }
    }
    Tensor_t df (dfEvent, 3);

    double W[][18] = {
            {1, 0.31, -0.35, -0.33, 0.40, 0.26, -0.30, 0.29, -0.31,
             0.21, 0.44, -0.36, 0.03, -0.27, -0.53, 0.24, 0.22, -0.35},

            {-0.06, -0.35, 0.10, -0.49, -0.88, 0.35, 0.03, 0, -0.19,
             -0.09, 0.29, 0.10, -0.15, -0.11, 0.02, 0.08, 0.17, 0.35}
    };

    Matrix_t weights(numberFilters, imgDepth * fltHeight * fltWidth);

    for (size_t i = 0; i < numberFilters; i++) {
        for (size_t j = 0; j < imgDepth * fltHeight * fltWidth; j++) {
            weights(i, j) = W[i][j];
        }
    }

    double activationsPreviousLayer[][25] = {
            {-0.33, -0.05, -0.34, -0.15, -0.13, 0.44, 0.51, 0.26, 1.21,
             -0.11, 0.86, 0, -0.56, -1.11, 1.22, -0.02, 0.29, 2.32,
             1.87, -0.31, -0.08, -0.20, -1.62, -1.40, -0.80},

            {0.42, -1.12, 1.10, -0.73, 0.38, 0.66, -2.36, 1.13, -1.19, 0.39,
             0.28, 0.04, 1.01, 0.31, -1.87, -0.09, 0.21, -2.51, 0.11, 0.84,
             -0.04, 0.56, 0.48, 0.09, -0.08}
    };

    Matrix_t activationsBackwardEvent(imgDepth, imgHeight * imgWidth);

    for (size_t i = 0; i < imgDepth; i++) {
        for (size_t j = 0; j < imgHeight * imgWidth; j++) {
            activationsBackwardEvent(i, j) = activationsPreviousLayer[i][j];
        }
    }
    Tensor_t activationsBackward ( activationsBackwardEvent, 3);


    /////////////////////// Fill the expected output //////////////////////////
    double expectedActivationGradsBackward[][25] = {
            {0, 1.39, 0.55, -0.51, 0, 0, -0.27, 0.88,  0.31,  -0.02, -0.01, -1.41, 0.26,  0.18,  -0.08, -0.12, 0.06,
             -0.27, -0.23, 0.04,  0,    0.27,  -0.31, 0.27, 0},

            {0, 0.32, 0.49, -0.53, 0, 0, 0.09,  -0.30, -0.80, -0.02, -0.02, 0.18,  -0.09, -0.25, 0,     -0.03, -0.05,
             0.22,  0.43,  -0.08, 0.02, -0.17, -0.10, 0.31, 0}
    };

    Matrix_t expectedActivationGradientsBackwardEvent(imgDepth, imgHeight * imgWidth);

    for (size_t i = 0; i < imgDepth; i++) {
        for (size_t j = 0; j < imgHeight * imgWidth; j++) {
            expectedActivationGradientsBackwardEvent(i, j) = expectedActivationGradsBackward[i][j];
        }
    }

    Tensor_t computedActivationGradientsBackward (1,imgDepth, imgHeight * imgWidth);

    /////////////////////// Fill the expected weights gradients //////////////////////////
    double expectedWeightGrads[][18] = {
            {-0.06, 0.03, 0.79, 0.43, -1.73, -0.02, 0.18, 0.69, -0.26, -1.57, 0.59, -1.27, -3.42, 3.80,
             -1.72, -0.44, 0.95, 0.34},

            {0.17, -0.17, -0.06, -0.05, 0.25, -0.14, -0.60, -0.31, 0.06, 0.20, -0.09, 0.43, 0.59, -0.44,
             0.25, 0.60, -0.25, -0.19}
    };

    Matrix_t expectedWeightGradients(numberFilters, imgDepth * fltHeight * fltWidth);
    for (size_t i = 0; i < numberFilters; i++) {
        for (size_t j = 0; j < imgDepth * fltHeight * fltWidth; j++) {
            expectedWeightGradients(i, j) = expectedWeightGrads[i][j];
        }
    }

    /////////////////////// Fill the expected bias gradients //////////////////////////
    double expectedBiasGrads[][1] = {
            {0.47},
            {-0.36}
    };

    Matrix_t expectedBiasGradients(imgDepth, 1);
    for (size_t i = 0; i < imgDepth; i++) {
        expectedBiasGradients(i, 0) = expectedBiasGrads[i][0];
    }


    // Init outputs - these should be filled by the computation.
    Matrix_t computedWeightGradients(numberFilters, imgDepth * fltHeight * fltWidth);
    Matrix_t computedBiasGradients(numberFilters, 1);

    TDescriptors * convDescriptors = nullptr;
    TWorkspace   * convWorkspace   = nullptr;

    TConvParams params(1, imgDepth, imgHeight, imgWidth,
                       numberFilters, fltHeight, fltWidth,
                       strideRows, strideCols, zeroPaddingHeight, zeroPaddingWidth);

    TConvLayer<Architecture> *layer = nullptr;
    Architecture::InitializeConvDescriptors(convDescriptors, layer);
    Architecture::InitializeConvWorkspace(convWorkspace, convDescriptors, params, layer);

    Tensor_t output = df;

    Architecture::ConvLayerBackward(computedActivationGradientsBackward, computedWeightGradients, computedBiasGradients,
                                    df, activationGradients, weights, activationsBackward, output,
                                    EActivationFunction::kIdentity,
                                    (typename Architecture::ConvDescriptors_t &) * convDescriptors,
                                    (typename Architecture::ConvWorkspace_t &) * convWorkspace,
                                    batchSize, imgHeight, imgWidth, numberFilters, height,
                                    width, imgDepth, fltHeight, fltWidth, nLocalViews);


    // Check correctness.
    bool status = true;
    status &= Architecture::AlmostEquals(expectedActivationGradientsBackwardEvent, computedActivationGradientsBackward.At(0).GetMatrix());
    status &= Architecture::AlmostEquals(expectedWeightGradients, computedWeightGradients);
    status &= Architecture::AlmostEquals(expectedWeightGradients, computedWeightGradients);
    return status;
}

template<typename Architecture>
bool testBackward1_cudnn()
{
    using Matrix_t = typename Architecture::Matrix_t;
    using Tensor_t = typename Architecture::Tensor_t;
    using HostBuffer_t = typename Architecture::HostBuffer_t;

    size_t imgDepth = 2;
    size_t imgHeight = 5;
    size_t imgWidth = 5;
    size_t numberFilters = 2;
    size_t fltHeight = 3;
    size_t fltWidth = 3;
    size_t strideRows = 1;
    size_t strideCols = 1;
    size_t zeroPaddingHeight = 0;
    size_t zeroPaddingWidth = 0;

    size_t height = calculateDimension(imgHeight, fltHeight, zeroPaddingHeight, strideRows);
    size_t width = calculateDimension(imgWidth, fltWidth, zeroPaddingWidth, strideCols);
    size_t nLocalViews =  height * width;
    size_t batchSize = 1;

    double grad[][9] = {
            {0, 1.37, 0, 0, 0, 0, 0, -0.90, 0},
            {0, -0.37, 0, 0, 0, -0.25, 0.26, 0, 0}
    };

    std::vector<size_t> gradShape {1, numberFilters, height, width};
    HostBuffer_t grad_hostbuffer(numberFilters * height * width);
    for (size_t i = 0; i < numberFilters; i++) {
       for (size_t j = 0; j < height * width; j++) {
          grad_hostbuffer[i * height * width + j] = grad[i][j];
       }
    }
    Matrix_t activationGradients(grad_hostbuffer, gradShape, MemoryLayout::RowMajor, 0, 0);
    //activationGradients.Print();
    double derivatives[][9] = {
            {1, 1, 1, 1 , 1, 1, 1, 1, 1},
            {1, 1, 1, 1 , 1, 1, 1, 1, 1}
    };

    std::vector<size_t> derivShape {1, numberFilters, height, width};
    HostBuffer_t deriv_hostbuffer(numberFilters * height * width);
    for (size_t i = 0; i < numberFilters; i++) {
       for (size_t j = 0; j < height * width; j++) {
          deriv_hostbuffer[i * height * width + j] = derivatives[i][j];
       }
    }
    Matrix_t df(deriv_hostbuffer, derivShape, MemoryLayout::RowMajor, 0, 0);
    //df.Print();
    double W[][18] = {
            {1, 0.31, -0.35, -0.33, 0.40, 0.26, -0.30, 0.29, -0.31,
             0.21, 0.44, -0.36, 0.03, -0.27, -0.53, 0.24, 0.22, -0.35},

            {-0.06, -0.35, 0.10, -0.49, -0.88, 0.35, 0.03, 0, -0.19,
             -0.09, 0.29, 0.10, -0.15, -0.11, 0.02, 0.08, 0.17, 0.35}
    };


    std::vector<size_t> weightsShape {numberFilters, imgDepth, fltHeight, fltWidth};
    HostBuffer_t weights_hostbuffer(numberFilters * imgDepth * fltHeight * fltWidth);
    for (size_t i = 0; i < numberFilters; i++) {
        for (size_t j = 0; j < imgDepth * fltHeight * fltWidth; j++) {
            weights_hostbuffer[i * imgDepth * fltHeight * fltWidth + j] = W[i][j];
        }
    }
    Tensor_t weights(weights_hostbuffer, weightsShape, MemoryLayout::RowMajor, 0, 0);
    //weights.Print();


    /// input x
    double activationsPreviousLayer[][25] = {
            {-0.33, -0.05, -0.34, -0.15, -0.13, 0.44, 0.51, 0.26, 1.21,
             -0.11, 0.86, 0, -0.56, -1.11, 1.22, -0.02, 0.29, 2.32,
             1.87, -0.31, -0.08, -0.20, -1.62, -1.40, -0.80},

            {0.42, -1.12, 1.10, -0.73, 0.38, 0.66, -2.36, 1.13, -1.19, 0.39,
             0.28, 0.04, 1.01, 0.31, -1.87, -0.09, 0.21, -2.51, 0.11, 0.84,
             -0.04, 0.56, 0.48, 0.09, -0.08}
    };

    std::vector<size_t> inputActvShape {1, imgDepth, imgHeight, imgWidth};
    HostBuffer_t inputActv_hostbuffer(imgDepth * imgHeight * imgWidth);
    for (size_t i = 0; i < imgDepth; i++) {
        for (size_t j = 0; j < imgHeight * imgWidth; j++) {
            inputActv_hostbuffer[i * imgHeight * imgWidth + j] = activationsPreviousLayer[i][j];
        }
    }
    Tensor_t input(inputActv_hostbuffer, inputActvShape, MemoryLayout::RowMajor, 0, 0);
    //activationsBackward.Print();
    /////////////////////// Fill the expected output //////////////////////////
    // dx
    double expectedActivationGradsBackward[][25] = {
            {0, 1.39, 0.55, -0.51, 0, 0, -0.27, 0.88,  0.31,  -0.02, -0.01, -1.41, 0.26,  0.18,  -0.08, -0.12, 0.06,
             -0.27, -0.23, 0.04,  0,    0.27,  -0.31, 0.27, 0},

            {0, 0.32, 0.49, -0.53, 0, 0, 0.09,  -0.30, -0.80, -0.02, -0.02, 0.18,  -0.09, -0.25, 0,     -0.03, -0.05,
             0.22,  0.43,  -0.08, 0.02, -0.17, -0.10, 0.31, 0}
    };

    std::vector<size_t> expcActvShape {1, imgDepth, imgHeight, imgWidth};
    HostBuffer_t expcActv_hostbuffer(imgDepth * imgHeight * imgWidth);
    for (size_t i = 0; i < imgDepth; i++) {
        for (size_t j = 0; j < imgHeight * imgWidth; j++) {
            expcActv_hostbuffer[i * imgHeight * imgWidth + j] = expectedActivationGradsBackward[i][j];
        }
    }
    Tensor_t expectedActivationGradientsBackward(expcActv_hostbuffer, expcActvShape, MemoryLayout::RowMajor, 0, 0);
    //computedActivationGradientsBackward.Print();

    /////////////////////// Fill the expected weights gradients //////////////////////////
    double expectedWeightGrads[][18] = {
            {-0.06, 0.03, 0.79, 0.43, -1.73, -0.02, 0.18, 0.69, -0.26, -1.57, 0.59, -1.27, -3.42, 3.80,
             -1.72, -0.44, 0.95, 0.34},

            {0.17, -0.17, -0.06, -0.05, 0.25, -0.14, -0.60, -0.31, 0.06, 0.20, -0.09, 0.43, 0.59, -0.44,
             0.25, 0.60, -0.25, -0.19}
    };

    std::vector<size_t> expcWeightsShape {numberFilters, imgDepth, fltHeight, fltWidth};
    HostBuffer_t expc_weights_hostbuffer(numberFilters * imgDepth * fltHeight * fltWidth);
    for (size_t i = 0; i < numberFilters; i++) {
        for (size_t j = 0; j < imgDepth * fltHeight * fltWidth; j++) {
            expc_weights_hostbuffer[i * imgDepth * fltHeight * fltWidth + j] = expectedWeightGrads[i][j];
        }
    }
    Tensor_t expectedWeightGradients(expc_weights_hostbuffer, expcWeightsShape, MemoryLayout::RowMajor, 0, 0);
    //expectedWeightGradients.Print();
    /////////////////////// Fill the expected bias gradients //////////////////////////

    std::vector<size_t> biasShape = {1,2, 1, 1 };
    Tensor_t biasesTensor( biasShape, MemoryLayout::RowMajor, 0, 0);
    biasesTensor.Zero();


    double expectedBiasGrads[][1] = {
            {0.47},
            {-0.36}
    };

    Matrix_t expectedBiasGradients( biasShape,MemoryLayout::RowMajor, 0, 0);
    for (size_t i = 0; i < imgDepth; i++) {
         expectedBiasGradients(0, i, 0, 0) = expectedBiasGrads[i][0];
    }

    // Init outputs - these should be filled by the computation.
    Matrix_t computedWeightGradients( weightsShape,  MemoryLayout::RowMajor);
    Matrix_t computedBiasGradients(  biasShape,MemoryLayout::RowMajor, 0, 0);

    Tensor_t computedOutput( { 1,2,3,3}, MemoryLayout::RowMajor, 0, 0);
    Tensor_t computedInputActivFunc( { 1,2,3,3}, MemoryLayout::RowMajor, 0, 0);

    // Make a forward pass in preparation
    EActivationFunction AFunct = EActivationFunction::kIdentity;
    //EActivationFunction AFunct = EActivationFunction::kRelu;
    TConvLayer<Architecture> convLayer (1, imgDepth, imgHeight, imgWidth, numberFilters,
                                            EInitialization::kIdentity, fltHeight, fltWidth,
                                            strideRows, strideCols, zeroPaddingHeight, zeroPaddingWidth,
                                            0.0, AFunct, ERegularization::kNone, 0.0);

    auto& convDescriptors = static_cast<typename Architecture::ConvDescriptors_t &> (*convLayer.GetDescriptors());
    auto& convWorkspace = static_cast<typename Architecture::ConvWorkspace_t &> (*convLayer.GetWorkspace());

    TConvParams params(1, imgDepth, imgHeight, imgWidth, numberFilters, fltHeight, fltWidth, strideRows,
                      strideCols, zeroPaddingHeight, zeroPaddingWidth);
    Tensor_t dummy;
    Architecture::ConvLayerForward(computedOutput, computedInputActivFunc, input, weights, biasesTensor, params,
                                       AFunct, dummy, convDescriptors, convWorkspace);


    // Backward pass
    HostBuffer_t comp_actvGrad_hostbuffer(numberFilters * imgDepth * fltHeight * fltWidth);
    Tensor_t computedActivationGradientsBackward(comp_actvGrad_hostbuffer, expcActvShape, MemoryLayout::RowMajor);

    Architecture::ConvLayerBackward(computedActivationGradientsBackward, computedWeightGradients,
                                    computedBiasGradients,
                                    computedInputActivFunc , activationGradients, weights, input, computedOutput,
                                    EActivationFunction::kIdentity,  // this is not used in cudnn
                                    convDescriptors, convWorkspace, batchSize, imgHeight, imgWidth, numberFilters, height,
                                      width, imgDepth, fltHeight, fltWidth, nLocalViews);
    // Check correctness.
    bool status = true;

    Architecture::PrintTensor( expectedActivationGradientsBackward, "expected dx");
    Architecture::PrintTensor( computedActivationGradientsBackward, "computed dx");

    Architecture::PrintTensor( expectedWeightGradients, "expected dw");
    Architecture::PrintTensor( computedWeightGradients, "computed dw");

    Architecture::PrintTensor( expectedBiasGradients, "expected db");
    Architecture::PrintTensor( computedBiasGradients, "computed db");

    // FIXME: Implement Almost Equals for CudaTensor
    /*status &= Architecture::AlmostEquals(expectedActivationGradientsBackwardEvent,    computedActivationGradientsBackward.At(0).GetMatrix());
    status &= Architecture::AlmostEquals(expectedWeightGradients, computedWeightGradients);
    status &= Architecture::AlmostEquals(expectedWeightGradients, computedWeightGradients);*/
    return status;
};
