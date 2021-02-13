// @(#)root/tmva $Id$
// Author: Dominik Dannheim, Alexander Voigt

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Classes: PDEFoamKernelGauss                                                    *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      PDEFoam kernel, which weights all cell values by a gauss function.        *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      S. Jadach        - Institute of Nuclear Physics, Cracow, Poland           *
 *      Tancredi Carli   - CERN, Switzerland                                      *
 *      Dominik Dannheim - CERN, Switzerland                                      *
 *      Alexander Voigt  - TU Dresden, Germany                                    *
 *                                                                                *
 * Copyright (c) 2008, 2010:                                                      *
 *      CERN, Switzerland                                                         *
 *      MPI-K Heidelberg, Germany                                                 *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

#ifndef ROOT_TMVA_PDEFoamKernelGauss
#define ROOT_TMVA_PDEFoamKernelGauss

#include "TMVA/PDEFoam.h"
#include "TMVA/PDEFoamKernelBase.h"
#include <vector>

namespace TMVA
{

   class PDEFoamKernelGauss : public PDEFoamKernelBase
   {

   protected:
      Float_t fSigma;          // width of gauss curve

      // Square function (fastest implementation)
      template<typename T> T Sqr(T x) const { return x * x; }

      // calculate gaussian weight
      Float_t WeightGaus(PDEFoam*, PDEFoamCell*, std::vector<Float_t>&);

      // estimate the cell value by its neighbors
      Float_t GetAverageNeighborsValue(PDEFoam*, std::vector<Float_t>&, ECellValue);

   public:
      PDEFoamKernelGauss(Float_t sigma);    // Constructor
      PDEFoamKernelGauss(const PDEFoamKernelGauss&); // Copy Constructor
      virtual ~PDEFoamKernelGauss() {}      // Destructor

      // kernel estimator
      virtual Float_t Estimate(PDEFoam*, std::vector<Float_t>&, ECellValue);

      ClassDef(PDEFoamKernelGauss, 1) // Gaussian PDEFoam kernel estimator
         }; // end of PDEFoamKernelGauss
}  // namespace TMVA

#endif
