// @(#)root/tmva/tmva/cnn:$Id$
// Author: Ashish Kshirsagar

#ifndef ROOT_TESTUPSAMPLELAYER_H
#define ROOT_TESTUPSAMPLELAYER_H

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  :                                                                       *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Testing the Upsample layer in an architecture agnostic manner.             *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Ashish Kshirsagar      <ashishkshirsagar10@gmail.com>                     *
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
// Testing the Upsampling Layer                                  //
////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>

#include "TestConvNet.h"

using namespace TMVA::DNN;
using namespace TMVA::DNN::CNN;


/*************************************************************************
 * Test 1:
 *  batch size = 2
 *  input depth = 2, input image height = 3, input image width = 3,
 *  output depth = 4, output image height = 6, output image width = 6,
 *  input depth = 2, input image height = 4, input image width = 2,
 *  output depth = 4, output image height = 8, output image width = 2,
 *  zero-padding height = 0, zero-padding width = 0,
 *************************************************************************/
template<typename Architecture>
bool testUpsample1()
{
    using Matrix_t = typename Architecture::Matrix_t;

     double imgTest1[][9] =
     {
       {
         10, 4, 22, 2, 18, 7, 9, 14, 25
       },
       {
         122, 31, 23, 54, 85, 14, 87, 67, 92
       }
     };
   
     double answerTest1[][36] = 
     {
       {
         10, 10, 10, 10, 4, 4, 4, 4, 22, 22, 22, 22, 2, 2, 2, 2, 18, 18, 18, 18, 7, 7, 7, 7, 9, 9, 9, 9, 14, 14, 14, 14, 25, 25, 25, 25
       },
       {
         10, 10, 10, 10, 4, 4, 4, 4, 22, 22, 22, 22, 2, 2, 2, 2, 18, 18, 18, 18, 7, 7, 7, 7, 9, 9, 9, 9, 14, 14, 14, 14, 25, 25, 25, 25
       }, 
       {
         122, 122, 122, 122, 31, 31, 31, 31, 23, 23, 23, 23, 54, 54, 54, 54, 85, 85, 85, 85, 14, 14, 14, 14, 87, 87, 87, 87, 67, 67, 67, 67, 92, 92, 92, 92
        },
       {
         122, 122, 122, 122, 31, 31, 31, 31, 23, 23, 23, 23, 54, 54, 54, 54, 85, 85, 85, 85, 14, 14, 14, 14, 87, 87, 87, 87, 67, 67, 67, 67, 92, 92, 92, 92
       }
     };

    double imgTest2[][8] =
      {
        {
          212,  213,  213,  150,
         227,  250,  250,  235
       },
        
        {
          255,  255,  192,  204,
         153,  246,  246,  175
       }
      };

    double answerTest2[][16] = 
    {
      {
          212, 212, 213, 213, 213, 213, 150, 150, 
          227, 227, 250, 250, 250, 250, 235, 235 
      },
      {
          212, 212, 213, 213, 213, 213, 150, 150, 
          227, 227, 250, 250, 250, 250, 235, 235 
      },
      {
          255, 255, 255, 255, 192, 192, 204, 204, 
          153, 153, 246, 246, 246, 246, 175, 175 
      },
      {
          255, 255, 255, 255, 192, 192, 204, 204, 
          153, 153, 246, 246, 246, 246, 175, 175
      }
    };

    size_t imgDepthTest1 = 2;
    size_t imgHeightTest1 = 3;
    size_t imgWidthTest1 = 3;
    size_t answerDepthTest1 = 4;
    size_t answerHeightTest1 = 6;
    size_t answerWidthTest1 = 6;
    size_t imgDepthTest2 = 2;
    size_t imgHeightTest2 = 4;
    size_t imgWidthTest2 = 2;
    size_t answerDepthTest2 = 4;
    size_t answerHeightTest2 = 8;
    size_t answerWidthTest2 = 2;


    Matrix_t A1(imgDepthTest1, imgHeightTest1 * imgWidthTest1);

    for (size_t i = 0; i < (size_t)A1.GetNrows(); i++) {
        for (size_t j = 0; j < (size_t)A1.GetNcols(); j++) {
          A1(i, j) = imgTest1[i][j];
        }
    }

    Matrix_t A2(imgDepthTest2, imgHeightTest2 * imgWidthTest2);

    for (size_t i = 0; i < (size_t)A2.GetNrows(); i++) {
        for (size_t j = 0; j < (size_t)A2.GetNcols(); j++) {
            A2(i, j) = imgTest2[i][j];
        }
    }

    std::vector<Matrix_t> input;
    input.emplace_back(A1);
    input.emplace_back(A2);

    Matrix_t B1(answerDepthTest1, answerHeightTest1 * answerWidthTest1);

    for (size_t i = 0; i < (size_t)B1.GetNrows(); i++) {
      for (size_t j = 0; j < (size_t)B1.GetNcols(); j++) {
          B1(i, j) = answerTest1[i][j];
      }
    }

    Matrix_t B2(answerDepthTest2, answerHeightTest2 * answerWidthTest2);

    for (size_t i = 0; i < (size_t)B2.GetNrows(); i++) {
        for (size_t j = 0; j < (size_t)B2.GetNcols(); j++) {
            B2(i, j) = answerTest2[i][j];
        }
    }

    std::vector<Matrix_t> output;
    output.emplace_back(B1);
    output.emplace_back(B2);

    bool status = testUpsample<Architecture>(input, output);

    return status;
}


/*************************************************************************
 * Test 1:
 *  batch size = 2
 *  input depth = 4, input image height = 8, input image width = 2,
 *  output depth = 2, output image height = 4, output image width = 2,
 *  input depth = 4, input image height = 6, input image width = 6,
 *  output depth = 2, output image height = 3, output image width = 3,
 *  zero-padding height = 0, zero-padding width = 0,
 *************************************************************************/
template<typename Architecture>
bool testBackward1()
{
    using Matrix_t = typename Architecture::Matrix_t;

    double answerTest1[][8] =
    {
      {
        212,  213,  213,  150,
       227,  250,  250,  235
     },
      
      {
        255,  255,  192,  204,
       153,  246,  246,  175
     }
    };



  double imgTest1[][16] = 
  {
    {
        212, 212, 213, 213, 213, 213, 150, 150, 
        227, 227, 250, 250, 250, 250, 235, 235 
    },
    {
        212, 212, 213, 213, 213, 213, 150, 150, 
        227, 227, 250, 250, 250, 250, 235, 235 
    },
    {
        255, 255, 255, 255, 192, 192, 204, 204, 
        153, 153, 246, 246, 246, 246, 175, 175 
    },
    {
        255, 255, 255, 255, 192, 192, 204, 204, 
        153, 153, 246, 246, 246, 246, 175, 175
    }
  };

    double imgTest2[][36] = 
    {
      {
        10, 10, 10, 10, 4, 4, 4, 4, 22, 22, 22, 22, 2, 2, 2, 2, 18, 18, 18, 18, 7, 7, 7, 7, 9, 9, 9, 9, 14, 14, 14, 14, 25, 25, 25, 25
      },
      {
        10, 10, 10, 10, 4, 4, 4, 4, 22, 22, 22, 22, 2, 2, 2, 2, 18, 18, 18, 18, 7, 7, 7, 7, 9, 9, 9, 9, 14, 14, 14, 14, 25, 25, 25, 25
      }, 
      {
        122, 122, 122, 122, 31, 31, 31, 31, 23, 23, 23, 23, 54, 54, 54, 54, 85, 85, 85, 85, 14, 14, 14, 14, 87, 87, 87, 87, 67, 67, 67, 67, 92, 92, 92, 92
      },
      {
        122, 122, 122, 122, 31, 31, 31, 31, 23, 23, 23, 23, 54, 54, 54, 54, 85, 85, 85, 85, 14, 14, 14, 14, 87, 87, 87, 87, 67, 67, 67, 67, 92, 92, 92, 92
      }
    };

    double answerTest2[][9] =
    {
      {
        10, 4, 22, 2, 18, 7, 9, 14, 25
      },
      {
        122, 31, 23, 54, 85, 14, 87, 67, 92
      }
    };

    size_t imgDepthTest1 = 4;
    size_t imgHeightTest1 = 8;
    size_t imgWidthTest1 = 2;
    size_t answerDepthTest1 = 2;
    size_t answerHeightTest1 = 4;
    size_t answerWidthTest1 = 2;
    size_t imgDepthTest2 = 4;
    size_t imgHeightTest2 = 6;
    size_t imgWidthTest2 = 6;
    size_t answerDepthTest2 = 2;
    size_t answerHeightTest2 = 3;
    size_t answerWidthTest2 = 3;


    Matrix_t A1(imgDepthTest1, imgHeightTest1 * imgWidthTest1);

    for (size_t i = 0; i < (size_t)A1.GetNrows(); i++) {
        for (size_t j = 0; j < (size_t)A1.GetNcols(); j++) {
            A1(i, j) = imgTest1[i][j];
        }
    }

    Matrix_t A2(imgDepthTest2, imgHeightTest2 * imgWidthTest2);

    for (size_t i = 0; i < (size_t)A2.GetNrows(); i++) {
        for (size_t j = 0; j < (size_t)A2.GetNcols(); j++) {
            A2(i, j) = imgTest2[i][j];
        }
    }

    std::vector<Matrix_t> input;
    input.emplace_back(A1);
    input.emplace_back(A2);
    
    Matrix_t B1(answerDepthTest1, answerHeightTest1 * answerWidthTest1);

    for (size_t i = 0; i < (size_t)B1.GetNrows(); i++) {
        for (size_t j = 0; j < (size_t)B1.GetNcols(); j++) {
            B1(i, j) = answerTest1[i][j];
        }
    }

    Matrix_t B2(answerDepthTest2, answerHeightTest2 * answerWidthTest2);

    for (size_t i = 0; i < (size_t)B2.GetNrows(); i++) {
        for (size_t j = 0; j < (size_t)B2.GetNcols(); j++) {
            B2(i, j) = answerTest2[i][j];
        }
    }
    
    std::vector<Matrix_t> output;
    output.emplace_back(B1);
    output.emplace_back(B2);

    bool status = testUpsampleBackward<Architecture>(input,output);

    return status;
}

#endif //ROOT_TESTPOOLINGLAYER_H
