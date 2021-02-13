// @(#)root/matrix:$Id$
// Authors: Fons Rademakers, Eddy Offermann   Nov 2003

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TMatrixFUtilsfwd
#define ROOT_TMatrixFUtilsfwd

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Matrix utility classes.                                              //
//                                                                      //
//  Forward declaration of                                              //
//   TMatrixTRow_const       <Float_t>  TMatrixTRow       <Float_t>     //
//   TMatrixTColumn_const    <Float_t>  TMatrixTColumn    <Float_t>     //
//   TMatrixTDiag_const      <Float_t>  TMatrixTDiag      <Float_t>     //
//   TMatrixTFlat_const      <Float_t>  TMatrixTFlat      <Float_t>     //
//   TMatrixTSub_const       <Float_t>  TMatrixTSub       <Float_t>     //
//   TMatrixTSparseRow_const <Float_t>  TMatrixTSparseRow <Float_t>     //
//   TMatrixTSparseDiag_const<Float_t>  TMatrixTSparseDiag<Float_t>     //
//                                                                      //
//   TElementActionT   <Float_t>                                        //
//   TElementPosActionT<Float_t>                                        //
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

typedef TMatrixTRow_const       <Float_t> TMatrixFRow_const;
typedef TMatrixTColumn_const    <Float_t> TMatrixFColumn_const;
typedef TMatrixTDiag_const      <Float_t> TMatrixFDiag_const;
typedef TMatrixTFlat_const      <Float_t> TMatrixFFlat_const;
typedef TMatrixTSub_const       <Float_t> TMatrixFSub_const;
typedef TMatrixTSparseRow_const <Float_t> TMatrixFSparseRow_const;
typedef TMatrixTSparseDiag_const<Float_t> TMatrixFSparseDiag_const;

typedef TMatrixTRow             <Float_t> TMatrixFRow;
typedef TMatrixTColumn          <Float_t> TMatrixFColumn;
typedef TMatrixTDiag            <Float_t> TMatrixFDiag;
typedef TMatrixTFlat            <Float_t> TMatrixFFlat;
typedef TMatrixTSub             <Float_t> TMatrixFSub;
typedef TMatrixTSparseRow       <Float_t> TMatrixFSparseRow;
typedef TMatrixTSparseDiag      <Float_t> TMatrixFSparseDiag;

typedef TElementActionT         <Float_t> TElementActionF;
typedef TElementPosActionT      <Float_t> TElementPosActionF;

#endif
