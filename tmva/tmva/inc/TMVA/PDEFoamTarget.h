// @(#)root/tmva $Id$
// Author: Tancredi Carli, Dominik Dannheim, Alexander Voigt

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Classes: PDEFoamTarget                                                         *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Concrete PDEFoam sub-class.  This foam stores the first target            *
 *      (index 0) with every cell, as well as the statistical error on            *
 *      the target.                                                               *
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

#ifndef ROOT_TMVA_PDEFoamTarget
#define ROOT_TMVA_PDEFoamTarget

#include "TMVA/PDEFoam.h"

#include <vector>

namespace TMVA
{

   class PDEFoamTarget : public PDEFoam
   {

   protected:
      UInt_t fTarget; // the target to fill the cells with

      // specific function used during evaluation; determines, whether
      // a cell value is undefined
      Bool_t CellValueIsUndefined(PDEFoamCell* cell);

      // calculate the average of the neighbor cell values
      Float_t GetAverageNeighborsValue(std::vector<Float_t>&, ECellValue);

      PDEFoamTarget(const PDEFoamTarget&); // Copy Constructor  NOT USED

      // ---------- Public functions ----------------------------------
   public:
      PDEFoamTarget();                  // Default constructor (used only by ROOT streamer)
      PDEFoamTarget(const TString&, UInt_t); // Principal user-defined constructor
      virtual ~PDEFoamTarget() {}       // Default destructor

      // function to fill created cell with given value
      virtual void FillFoamCells(const Event* ev, Float_t wt);

      // function to call after foam is grown
      virtual void Finalize();

      virtual Float_t GetCellValue(const std::vector<Float_t> &xvec, ECellValue cv, PDEFoamKernelBase*);
      using PDEFoam::GetCellValue;

      // ---------- ROOT class definition
      ClassDef(PDEFoamTarget, 1) // Tree of PDEFoamCells
         }; // end of PDEFoamTarget

}  // namespace TMVA

// ---------- Inline functions

#endif
