// @(#)root/tmva $Id$    
// Author: Andrzej Zemla

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : SVKernelFunction                                                      *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Kernel for Support Vector Machine                                         *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Marcin Wolter  <Marcin.Wolter@cern.ch> - IFJ PAN, Krakow, Poland          *
 *      Andrzej Zemla  <azemla@cern.ch>        - IFJ PAN, Krakow, Poland          *
 *      (IFJ PAN: Henryk Niewodniczanski Inst. Nucl. Physics, Krakow, Poland)     *   
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         * 
 *      MPI-K Heidelberg, Germany                                                 * 
 *      PAN, Krakow, Poland                                                       *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

#ifndef ROOT_TMVA_SVKernelFunction
#define ROOT_TMVA_SVKernelFunction

#include "RtypesCore.h"
#include <vector>

namespace TMVA {

   class SVEvent;
   class SVKernelFunction {

   public:

      enum EKernelType { kLinear , kRBF, kPolynomial, kSigmoidal, kMultiGauss, kProd, kSum};

      SVKernelFunction();
      SVKernelFunction( Float_t );
      SVKernelFunction( EKernelType, Float_t, Float_t=0);
      SVKernelFunction( std::vector<float> params );
      SVKernelFunction(EKernelType k, std::vector<EKernelType> kernels, std::vector<Float_t> gammas, Float_t gamma, Float_t order, Float_t theta);
      ~SVKernelFunction();
      
      Float_t Evaluate( SVEvent* ev1, SVEvent* ev2 );

      void setCompatibilityParams(EKernelType k, UInt_t order, Float_t theta, Float_t kappa);
         
   private:

      Float_t fGamma;   // documentation

      // vector of gammas for multidimensional gaussian
      std::vector<Float_t> fmGamma;

      // kernel, order, theta, and kappa are for backward compatibility
      EKernelType fKernel;
      UInt_t      fOrder;
      Float_t     fTheta;
      Float_t     fKappa;

      std::vector<EKernelType> fKernelsList;
   };
}

#endif
