// @(#)root/tmva $Id$
// Author: Tancredi Carli, Dominik Dannheim, Alexander Voigt

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Classes: PDEFoamDensityBase                                                    *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Class PDEFoamDensityBase is an Abstract class representing                *
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

#ifndef ROOT_TMVA_PDEFoamDensityBase
#define ROOT_TMVA_PDEFoamDensityBase

#include "TObject.h"
#include <vector>

#include "TMVA/BinarySearchTree.h"
#include "TMVA/Event.h"
#include "TMVA/MsgLogger.h"

namespace TMVA
{

   // class definition of underlying density
   class PDEFoamDensityBase : public ::TObject
      {
      private:
         std::vector<Double_t> fBox; // range-searching box
         Double_t fBoxVolume;        // volume of range searching box
         Bool_t fBoxHasChanged;      // range searching box has changed

      protected:
         BinarySearchTree *fBst;     // Binary tree to find events within a volume
         mutable MsgLogger *fLogger; //! message logger

         MsgLogger& Log() const { return *fLogger; }

         // calculate volume of fBox
         Double_t GetBoxVolume();

      public:
         PDEFoamDensityBase();
         PDEFoamDensityBase(std::vector<Double_t> box);
         PDEFoamDensityBase(const PDEFoamDensityBase&);
         virtual ~PDEFoamDensityBase();

         // fill event into binary search tree
         void FillBinarySearchTree(const Event* ev);

         // set the range-searching box
         void SetBox(std::vector<Double_t> box) { fBox = box; fBoxHasChanged = kTRUE; }

         // get the range-searching box
         const std::vector<Double_t>& GetBox() const { return fBox; }

         // main function used by PDEFoam
         // returns density at a given point by range searching in BST
         virtual Double_t Density(std::vector<Double_t> &Xarg, Double_t &event_density) = 0;

         ClassDef(PDEFoamDensityBase, 1) // PDEFoam event density interface
            };  //end of PDEFoamDensityBase

}  // namespace TMVA

#endif
