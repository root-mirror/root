/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TMatrixDfwd
#define ROOT_TMatrixDfwd

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMatrixD                                                             //
//                                                                      //
//  Forward declaration of TMatrixT<Double_t>                           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "RtypesCore.h"

template<class Element> class TMatrixT;
typedef TMatrixT<Double_t> TMatrixD;

#endif
