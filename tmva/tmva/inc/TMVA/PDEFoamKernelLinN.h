// @(#)root/tmva $Id$
// Author: Dominik Dannheim, Alexander Voigt

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Classes: PDEFoamKernelLinN                                                     *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      PDEFoam kernel, which linear weights with the neighbor cells.             *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      S. Jadach        - Institute of Nuclear Physics, Cracow, Poland           *
 *      Tancredi Carli   - CERN, Switzerland                                      *
 *      Dominik Dannheim - CERN, Switzerland                                      *
 *      Alexander Voigt  - TU Dresden, Germany                                    *
 *                                                                                *
 * Copyright (c) 2010:                                                            *
 *      CERN, Switzerland                                                         *
 *      MPI-K Heidelberg, Germany                                                 *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

#ifndef ROOT_TMVA_PDEFoamKernelLinN
#define ROOT_TMVA_PDEFoamKernelLinN

#include "TMVA/PDEFoam.h"
#include "TMVA/PDEFoamKernelBase.h"

#include <vector>

namespace TMVA
{

   class PDEFoamKernelLinN : public PDEFoamKernelBase
   {

   protected:
      // helper function, which linear weights with the neighbor cells
      Float_t WeightLinNeighbors(PDEFoam*, std::vector<Float_t>&, ECellValue, Bool_t);

      // helper function for WeightLinNeighbors()
      Float_t GetAverageNeighborsValue(PDEFoam*, std::vector<Float_t>&, ECellValue);

   public:
      PDEFoamKernelLinN();                 // Constructor
      PDEFoamKernelLinN(const PDEFoamKernelLinN&); // Copy Constructor
      virtual ~PDEFoamKernelLinN() {}      // Destructor

      // kernel estimator
      virtual Float_t Estimate(PDEFoam*, std::vector<Float_t>&, ECellValue);

      ClassDef(PDEFoamKernelLinN, 1) // next neighbor PDEFoam kernel estimator
         }; // end of PDEFoamKernelLinN
}  // namespace TMVA

#endif
