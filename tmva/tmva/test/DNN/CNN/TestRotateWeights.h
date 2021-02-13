// @(#)root/tmva/tmva/cnn:$Id$
// Author: Manos Stergiadis

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  :                                                                       *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Testing RotateWeights method on a CPU architecture                        *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Manos Stergiadis    <em.stergiadis@gmail.com>   - CERN, Switzerland       *
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
// Testing the Rotate Weights function                            //
////////////////////////////////////////////////////////////////////

#include <cmath>

#include "TestConvNet.h"

using namespace TMVA::DNN;
using namespace TMVA::DNN::CNN;

/*************************************************************************
 * Test 1:
 *  filter depth = 3, filter height = 2, filter width = 2, num. filters = 4
 *************************************************************************/
template<typename Architecture>
bool test1()
{
    using Matrix_t = typename Architecture::Matrix_t;

    double weights[][12] = {{252, 116, 155, 246, 170, 149, 227, 113, 166, 227, 119, 57},
                            {92, 103, 151, 37, 110, 46, 70, 8, 88, 182, 43, 236},
                            {153, 246, 216, 102, 179, 248, 187, 227, 66, 102, 180, 169},
                            {5, 215, 115, 103, 35, 138, 193, 28, 213, 93, 117, 208}};

    double answer[][16] = {{246, 155, 116, 252, 37, 151, 103, 92, 102, 216, 246, 153, 103, 115, 215, 5},
                           {113, 227, 149, 170, 8, 70, 46, 110, 227, 187, 248, 179, 28, 193, 138, 35},
                           {57, 119, 227, 166, 236, 43, 182, 88, 169, 180, 102, 66, 208, 117, 93, 213}};

    size_t filterDepth = 3;
    size_t filterHeight = 2;
    size_t filterWidth = 2;
    size_t numFilters = 4;

    Matrix_t A(numFilters, filterDepth * filterHeight * filterWidth);

    for (size_t i = 0; i < (size_t)A.GetNrows(); i++) {
        for (size_t j = 0; j < (size_t)A.GetNcols(); j++) {
            A(i, j) = weights[i][j];
        }
    }

    Matrix_t B(filterDepth, numFilters * filterHeight * filterWidth);

    for (size_t i = 0; i < (size_t)B.GetNrows(); i++) {
        for (size_t j = 0; j < (size_t)B.GetNcols(); j++) {
            B(i, j) = answer[i][j];
        }
    }

    return testRotateWeights<Architecture>(A, B, filterDepth, filterHeight, filterWidth, numFilters);
}

/*************************************************************************
 * Test 2:
 *  filter depth = 2, filter height = 2, filter width = 3, num. filters = 4
 *************************************************************************/
template<typename Architecture>
bool test2()
{
    using Matrix_t = typename Architecture::Matrix_t;

    double input[][12] = {{252, 116, 155, 246, 170, 149, 227, 113, 166, 227, 119, 57},
                          {92, 103, 151, 37, 110, 46, 70, 8, 88, 182, 43, 236},
                          {153, 246, 216, 102, 179, 248, 187, 227, 66, 102, 180, 169},
                          {5, 215, 115, 103, 35, 138, 193, 28, 213, 93, 117, 208}};

    double output[][24] = {{149,170,246,155,116,252,46,110,37,151,103,92,248,179,102,216,246,153,138,35,103,115,215,5},
                           {57,119,227,166,113,227,236,43,182,88,8,70,169,180,102,66,227,187,208,117,93,213,28,193}};

    size_t filterDepth = 2;
    size_t filterHeight = 2;
    size_t filterWidth = 3;
    size_t numFilters = 4;

    Matrix_t A(numFilters, filterDepth * filterHeight * filterWidth);

    for (size_t i = 0; i < (size_t)A.GetNrows(); i++) {
        for (size_t j = 0; j < (size_t)A.GetNcols(); j++) {
            A(i, j) = input[i][j];
        }
    }

    Matrix_t B(filterDepth, numFilters * filterHeight * filterWidth);

    for (size_t i = 0; i < (size_t)B.GetNrows(); i++) {
        for (size_t j = 0; j < (size_t)B.GetNcols(); j++) {
            B(i, j) = output[i][j];
        }
    }

    return testRotateWeights<Architecture>(A, B, filterDepth, filterHeight, filterWidth, numFilters);
}
