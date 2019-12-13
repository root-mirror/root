// @(#)root/tmva/tmva/cnn:$Id$
// Author: Ashish Kshirsagar

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  :                                                                       *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Testing Transpose convolution method on a CPU architecture                *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Ashish Kshirsagar       <ashishkshirsagar10@gmail.com>                    *
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
// Testing the Transpose Convolutional Layer                                //
////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>

#include "TestConvNet.h"

using namespace TMVA::DNN;
using namespace TMVA::DNN::CNN;

/*************************************************************************
 * Test 1: Forward Propagation
 *  batch size = 1
 *  image depth = 1, image height = 2, image width = 2,
 *  num frames = 1, filter height = 3, filter width = 3,
 *  stride rows = 1, stride cols = 1,
 *  zero-padding height = 0, zero-padding width = 0,
 *************************************************************************/
template<typename Architecture>
bool testForward1()
{
   using Matrix_t = typename Architecture::Matrix_t;

   double expected[][16] = {

    {
      2, 9, 6, 3, 
      9, 14, 29, 29,
      13, 21, 28, 24,
      24, 16, 4, 0
    }

   };

   double weights[][9] = {

    {
      1, 4, 1, 
      1, 4, 3,
      3, 3, 1
    }

   };

   double biases[][1] = {

      {
        0
      }

   };

   double img[][4] = {

      {
        2, 1, 4, 4
      }
    };
   
   size_t imgDepth = 1;
   size_t imgHeight = 2;
   size_t imgWidth = 2;
   size_t numberFilters = 1;
   size_t fltHeight = 3;
   size_t fltWidth = 3;
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
   std::vector<Matrix_t> input;
   input.push_back(inputEvent);

   Matrix_t weightsMatrix(numberFilters, fltHeight * fltWidth * imgDepth);
   Matrix_t biasesMatrix(numberFilters, 1);
   for (size_t i = 0; i < numberFilters; i++) {
       for (size_t j = 0; j < fltHeight * fltWidth * imgDepth; j++){
           weightsMatrix(i, j) = weights[i][j];
       }
       biasesMatrix(i, 0) = biases[i][0];
   }

   size_t height = 4;//calculateDimension(imgHeight, fltHeight, zeroPaddingHeight, strideRows);
   size_t width = 4;//calculateDimension(imgWidth, fltWidth, zeroPaddingWidth, strideCols);

   Matrix_t outputEvent(numberFilters, height * width);

   for (size_t i = 0; i < numberFilters; i++) {
      for (size_t j = 0; j < height * width; j++) {
         outputEvent(i, j) = expected[i][j];
      }
   }
   
   std::vector<Matrix_t> expectedOutput;
   expectedOutput.push_back(outputEvent);

   bool status = testTransConvLayerForward<Architecture>(input, expectedOutput, weightsMatrix, biasesMatrix, imgHeight,
                                                    imgWidth, imgDepth, fltHeight, fltWidth, numberFilters, strideRows,
                                                    strideCols, zeroPaddingHeight, zeroPaddingWidth);

   return status;
}

/*************************************************************************
 * Test 1: Backward Propagation
 *  batch size = 1
 *  image depth = 1, image height = 4, image width = 4,
 *  output depth = 1, output height = 2, output width = 2,
 *  zero-padding height = 0, zero-padding width = 0,
 *************************************************************************/
template<typename Architecture>
bool testBackward1()
{
   using Matrix_t = typename Architecture::Matrix_t;

   double img[][16] = {

    {
      2, 9, 6, 3, 
      9, 14, 29, 29,
      13, 21, 28, 24,
      24, 16, 4, 0
    }

   };

   double weights[][9] = {

    {
      1, 4, 1, 
      1, 4, 3,
      3, 3, 1
    }

   };

   double biases[][1] = {

      {
        0
      }

   };

   double expected[][4] = {

      {
        312, 335, 500, 487
      }
    };
   
   size_t imgDepth = 1;
   size_t imgHeight = 4;
   size_t imgWidth = 4;
   size_t numberFilters = 1;
   size_t fltHeight = 3;
   size_t fltWidth = 3;
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
   std::vector<Matrix_t> input;
   input.push_back(inputEvent);

   Matrix_t weightsMatrix(numberFilters, fltHeight * fltWidth * imgDepth);
   Matrix_t biasesMatrix(numberFilters, 1);
   for (size_t i = 0; i < numberFilters; i++) {
       for (size_t j = 0; j < fltHeight * fltWidth * imgDepth; j++){
           weightsMatrix(i, j) = weights[i][j];
       }
       biasesMatrix(i, 0) = biases[i][0];
   }

   size_t height = 2;//calculateDimension(imgHeight, fltHeight, zeroPaddingHeight, strideRows);
   size_t width = 2;//calculateDimension(imgWidth, fltWidth, zeroPaddingWidth, strideCols);

   Matrix_t outputEvent(numberFilters, height * width);

   for (size_t i = 0; i < numberFilters; i++) {
      for (size_t j = 0; j < height * width; j++) {
         outputEvent(i, j) = expected[i][j];
      }
   }
   
   std::vector<Matrix_t> expectedOutput;
   expectedOutput.push_back(outputEvent);

   bool status = testTransConvLayerBackward<Architecture>(input, expectedOutput, weightsMatrix, biasesMatrix, imgHeight,
                                                    imgWidth, imgDepth, fltHeight, fltWidth, numberFilters, strideRows,
                                                    strideCols, zeroPaddingHeight, zeroPaddingWidth);

   return status;
}