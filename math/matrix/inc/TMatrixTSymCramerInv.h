// @(#)root/base:$Id$
// Authors: Fons Rademakers, Eddy Offermann  Oct 2004

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TMatrixTSymCramerInv
#define ROOT_TMatrixTSymCramerInv

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMatrixTSymCramerInv                                                 //
//                                                                      //
// Encapsulate templates of Cramer Inversion routines.                  //
//                                                                      //
// The 4x4, 5x5 and 6x6 are adapted from routines written by            //
// Mark Fischler and Steven Haywood as part of the CLHEP package        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "RtypesCore.h"

#include "TMatrixTSym.h"

namespace TMatrixTSymCramerInv {

   template<class Element> Bool_t Inv2x2(TMatrixTSym<Element> &m,Double_t *determ);
   template<class Element> Bool_t Inv3x3(TMatrixTSym<Element> &m,Double_t *determ);
   template<class Element> Bool_t Inv4x4(TMatrixTSym<Element> &m,Double_t *determ);
   template<class Element> Bool_t Inv5x5(TMatrixTSym<Element> &m,Double_t *determ);
   template<class Element> Bool_t Inv6x6(TMatrixTSym<Element> &m,Double_t *determ);

}

#endif
