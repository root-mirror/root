// @(#)root/matrix:$Id$
// Authors: Fons Rademakers, Eddy Offermann   Nov 2003

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TMatrixDUtilsfwd
#define ROOT_TMatrixDUtilsfwd

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Matrix utility classes.                                              //
//                                                                      //
//  Forward declaration of                                              //
//   TMatrixTRow_const       <Double_t>  TMatrixTRow       <Double_t>   //
//   TMatrixTColumn_const    <Double_t>  TMatrixTColumn    <Double_t>   //
//   TMatrixTDiag_const      <Double_t>  TMatrixTDiag      <Double_t>   //
//   TMatrixTFlat_const      <Double_t>  TMatrixTFlat      <Double_t>   //
//   TMatrixTSub_const       <Double_t>  TMatrixTSub       <Double_t>   //
//   TMatrixTSparseRow_const <Double_t>  TMatrixTSparseRow <Double_t>   //
//   TMatrixTSparseDiag_const<Double_t>  TMatrixTSparseDiag<Double_t>   //
//                                                                      //
//   TElementActionT   <Double_t>                                       //
//   TElementPosActionT<Double_t>                                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "RtypesCore.h"

template<class Element> class TMatrixTRow_const;
template<class Element> class TMatrixTColumn_const;
template<class Element> class TMatrixTDiag_const;
template<class Element> class TMatrixTFlat_const;
template<class Element> class TMatrixTSub_const;
template<class Element> class TMatrixTSparseRow_const;
template<class Element> class TMatrixTSparseDiag_const;

template<class Element> class TMatrixTRow;
template<class Element> class TMatrixTColumn;
template<class Element> class TMatrixTDiag;
template<class Element> class TMatrixTFlat;
template<class Element> class TMatrixTSub;
template<class Element> class TMatrixTSparseRow;
template<class Element> class TMatrixTSparseDiag;

template<class Element> class TElementActionT;
template<class Element> class TElementPosActionT;

typedef TMatrixTRow_const       <Double_t> TMatrixDRow_const;
typedef TMatrixTColumn_const    <Double_t> TMatrixDColumn_const;
typedef TMatrixTDiag_const      <Double_t> TMatrixDDiag_const;
typedef TMatrixTFlat_const      <Double_t> TMatrixDFlat_const;
typedef TMatrixTSub_const       <Double_t> TMatrixDSub_const;
typedef TMatrixTSparseRow_const <Double_t> TMatrixDSparseRow_const;
typedef TMatrixTSparseDiag_const<Double_t> TMatrixDSparseDiag_const;

typedef TMatrixTRow             <Double_t> TMatrixDRow;
typedef TMatrixTColumn          <Double_t> TMatrixDColumn;
typedef TMatrixTDiag            <Double_t> TMatrixDDiag;
typedef TMatrixTFlat            <Double_t> TMatrixDFlat;
typedef TMatrixTSub             <Double_t> TMatrixDSub;
typedef TMatrixTSparseRow       <Double_t> TMatrixDSparseRow;
typedef TMatrixTSparseDiag      <Double_t> TMatrixDSparseDiag;

typedef TElementActionT         <Double_t> TElementActionD;
typedef TElementPosActionT      <Double_t> TElementPosActionD;

#endif
