// @(#)root/tmva $Id$
// Author: Dominik Dannheim, Alexander Voigt

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Classes: PDEFoamKernelBase                                                     *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      PDEFoam kernel interface                                                  *
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

#ifndef ROOT_TMVA_PDEFoamKernelBase
#define ROOT_TMVA_PDEFoamKernelBase

#include "TObject.h"

#include "TMVA/PDEFoam.h"

#include <vector>

namespace TMVA
{

   class PDEFoamKernelBase : public TObject
   {

   protected:
      mutable MsgLogger* fLogger;  //! message logger

   public:
      PDEFoamKernelBase();                 // Constructor
      PDEFoamKernelBase(const PDEFoamKernelBase&); // Copy constructor
      virtual ~PDEFoamKernelBase();        // Destructor

      // kernel estimator
      virtual Float_t Estimate(PDEFoam*, std::vector<Float_t>&, ECellValue) = 0;

      // Message logger
      MsgLogger& Log() const { return *fLogger; }

      ClassDef(PDEFoamKernelBase, 1) // PDEFoam kernel interface
         }; // end of PDEFoamKernelBase
}  // namespace TMVA

#endif
