// @(#)root/tmva $Id$
// Author: Manos Stergiadis

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  :                                                                       *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Testing Flatten function for every architecture using templates           *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Manos Stergiadis      <em.stergiadis@gmail.com - CERN, Switzerland        *
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


#ifndef ROOT_TESTFLATTEN_H
#define ROOT_TESTFLATTEN_H

#include "TestConvNet.h"

using namespace TMVA::DNN;
using namespace TMVA::DNN::CNN;

/////////////////////////////////////////////////////////////////////////
/// Test Reshape:
/// Input Shape: (3, 5)
/// Output Shape: (5, 3)
/////////////////////////////////////////////////////////////////////////
template<typename Architecture_t>
bool testReshape()
{
   using Matrix_t = typename Architecture_t::Matrix_t;

   double input[][5] = {{158, 157, 22,  166, 179},
                        { 68, 179, 233, 110, 163},
                        {168, 216,  76,   8, 102}};

   size_t nRowsA = 3;
   size_t nColsA = 5;
   Matrix_t A(nRowsA, nColsA);
   for (size_t i = 0; i < nRowsA; i++) {
      for (size_t j = 0; j < nColsA; j++) {
         A(i, j) = input[i][j];
      }
   }

   double expected[][3] = {{158, 157,  22},
                           {166, 179,  68},
                           {179, 233, 110},
                           {163, 168, 216},
                           { 76,   8, 102}};

   size_t nRowsB = 5;
   size_t nColsB = 3;
   Matrix_t B(nRowsB, nColsB);
   for (size_t i = 0; i < nRowsB; i++) {
      for (size_t j = 0; j < nColsB; j++) {
         B(i, j) = expected[i][j];
      }
   }

   return testReshape<Architecture_t>(A, B);
}

/*************************************************************************
 * Test 1:
 * depth = 3, width = 5, height = 5
 *************************************************************************/
template<typename Architecture_t>
bool testFlatten()
{
    //using Matrix_t = typename Architecture_t::Matrix_t;
    using Tensor_t = typename Architecture_t::Tensor_t;

    double input[][5][5] = {{{158, 157, 22, 166, 179},
                                {68, 179, 233, 110, 163},
                                {168, 216, 76, 8, 102},
                                {159, 163, 25, 78, 119},
                                {116, 50, 206, 102, 247},
                               },

                               {{187, 166, 121, 112, 136},
                                {237, 30, 180, 7, 248},
                                {52, 172, 146, 130, 92},
                                {124, 244, 214, 175, 9},
                                {80, 232, 139, 224, 237}},

                               {{53, 147, 103, 53, 110},
                                {112, 222, 19, 156, 232},
                                {81, 19, 188, 224, 220},
                                {255, 190, 76, 219, 95},
                                {245, 4, 217, 22, 22}}};

    double expected[][25] = {{158, 157, 22,  166, 179, 68, 179, 233, 110, 163, 168, 216, 76,
                                   8, 102, 159, 163, 25,  78, 119, 116, 50,  206, 102, 247},

                                {187, 166, 121, 112, 136, 237, 30, 180, 7,   248, 52,  172, 146,
                                 130, 92,  124, 244, 214, 175, 9,  80,  232, 139, 224, 237},

                                { 53,  147, 103, 53,  110, 112, 222, 19,  156, 232, 81, 19, 188,
                                 224, 220, 255, 190, 76,  219, 95,  245, 4,   217, 22, 22}};

    size_t size = 3;
    size_t nRows = 5;
    size_t nCols = 5;

    Tensor_t A(size, nRows, nCols);
    for (size_t i = 0; i < size; i++) {
        for (size_t j = 0; j < nRows; j++) {
            for (size_t k = 0; k < nCols; k++) {
                A(i, j, k) = input[i][j][k];
            }
        }
    }

    Tensor_t B(1, size, nRows * nCols);
    for (size_t i = 0; i < size; i++) {
        for (size_t j = 0; j < nRows * nCols; j++) {
            B(0, i, j) = expected[i][j];
        }
    }

    return testFlatten<Architecture_t>(A, B);
}

/*************************************************************************
 * Test 1:
 * depth = 3, width = 5, height = 5
 *************************************************************************/
template<typename Architecture_t>
bool testDeflatten()
{
    //using Matrix_t = typename Architecture_t::Matrix_t;
    using Tensor_t = typename Architecture_t::Tensor_t;

    double input[][25] = {{158, 157, 22,  166, 179, 68, 179, 233, 110, 163, 168, 216, 76,
                           8, 102, 159, 163, 25,  78, 119, 116, 50,  206, 102, 247},

                           {187, 166, 121, 112, 136, 237, 30, 180, 7,   248, 52,  172, 146,
                            130, 92,  124, 244, 214, 175, 9,  80,  232, 139, 224, 237},

                           { 53,  147, 103, 53,  110, 112, 222, 19,  156, 232, 81, 19, 188,
                            224, 220, 255, 190, 76,  219, 95,  245, 4,   217, 22, 22}};

    double expected[][5][5] = {{{158, 157, 22, 166, 179},
                                {68, 179, 233, 110, 163},
                                {168, 216, 76, 8, 102},
                                {159, 163, 25, 78, 119},
                                {116, 50, 206, 102, 247},
                               },

                               {{187, 166, 121, 112, 136},
                                {237, 30, 180, 7, 248},
                                {52, 172, 146, 130, 92},
                                {124, 244, 214, 175, 9},
                                {80, 232, 139, 224, 237}},

                               {{53, 147, 103, 53, 110},
                                {112, 222, 19, 156, 232},
                                {81, 19, 188, 224, 220},
                                {255, 190, 76, 219, 95},
                                {245, 4, 217, 22, 22}}};

    size_t size = 3;
    size_t nRows = 5;
    size_t nCols = 5;

    Tensor_t A(1, size, nRows * nCols);
    for (size_t i = 0; i < size; i++) {
        for (size_t j = 0; j < nRows * nCols; j++) {
            A(0, i, j) = input[i][j];
        }
    }

    Tensor_t B(size, nRows, nCols);
    for (size_t i = 0; i < size; i++) {
        for (size_t j = 0; j < nRows; j++) {
            for (size_t k = 0; k < nCols; k++) {
                B(i, j, k) = expected[i][j][k];
            }
        }
    }

    return testDeflatten<Architecture_t>(A, B);

}

#endif //ROOT_TESTFLATTEN_H
