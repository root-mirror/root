// @(#)root/tmva $Id$
// Author: Tancredi Carli, Dominik Dannheim, Alexander Voigt

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Classes: PDEFoamDiscriminantDensity                                            *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Class PDEFoamDiscriminantDensity is a class representing                  *
 *      n-dimensional real positive integrand function                            *
 *      The main function is Density() which provides the event density at a      *
 *      given point during the foam build-up (sampling).                          *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Tancredi Carli   - CERN, Switzerland                                      *
 *      Dominik Dannheim - CERN, Switzerland                                      *
 *      S. Jadach        - Institute of Nuclear Physics, Cracow, Poland           *
 *      Alexander Voigt  - TU Dresden, Germany                                    *
 *      Peter Speckmayer - CERN, Switzerland                                      *
 *                                                                                *
 * Copyright (c) 2008, 2010:                                                      *
 *      CERN, Switzerland                                                         *
 *      MPI-K Heidelberg, Germany                                                 *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

#ifndef ROOT_TMVA_PDEFoamDiscriminantDensity
#define ROOT_TMVA_PDEFoamDiscriminantDensity

#include "TMVA/PDEFoamDensityBase.h"
#include <vector>

namespace TMVA
{

   // class definition of underlying event density
   class PDEFoamDiscriminantDensity : public PDEFoamDensityBase
   {

   protected:
      UInt_t fClass; // signal class

   public:
      PDEFoamDiscriminantDensity();
      PDEFoamDiscriminantDensity(std::vector<Double_t> box, UInt_t cls);
      PDEFoamDiscriminantDensity(const PDEFoamDiscriminantDensity&);
      virtual ~PDEFoamDiscriminantDensity() {}

      // main function used by PDEFoam
      // returns discriminant density N_class/N_total at a given point
      // by range searching in BST
      virtual Double_t Density(std::vector<Double_t> &Xarg, Double_t &event_density);

      ClassDef(PDEFoamDiscriminantDensity, 1) //Class for Discriminant density
         };  //end of PDEFoamDiscriminantDensity

}  // namespace TMVA

#endif
