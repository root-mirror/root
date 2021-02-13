// @(#)root/quadp:$Id$
// Author: Eddy Offermann   May 2004

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/*************************************************************************
 * Parts of this file are copied from the OOQP distribution and          *
 * are subject to the following license:                                 *
 *                                                                       *
 * COPYRIGHT 2001 UNIVERSITY OF CHICAGO                                  *
 *                                                                       *
 * The copyright holder hereby grants you royalty-free rights to use,    *
 * reproduce, prepare derivative works, and to redistribute this software*
 * to others, provided that any changes are clearly documented. This     *
 * software was authored by:                                             *
 *                                                                       *
 *   E. MICHAEL GERTZ      gertz@mcs.anl.gov                             *
 *   Mathematics and Computer Science Division                           *
 *   Argonne National Laboratory                                         *
 *   9700 S. Cass Avenue                                                 *
 *   Argonne, IL 60439-4844                                              *
 *                                                                       *
 *   STEPHEN J. WRIGHT     swright@cs.wisc.edu                           *
 *   Computer Sciences Department                                        *
 *   University of Wisconsin                                             *
 *   1210 West Dayton Street                                             *
 *   Madison, WI 53706   FAX: (608)262-9777                              *
 *                                                                       *
 * Any questions or comments may be directed to one of the authors.      *
 *                                                                       *
 * ARGONNE NATIONAL LABORATORY (ANL), WITH FACILITIES IN THE STATES OF   *
 * ILLINOIS AND IDAHO, IS OWNED BY THE UNITED STATES GOVERNMENT, AND     *
 * OPERATED BY THE UNIVERSITY OF CHICAGO UNDER PROVISION OF A CONTRACT   *
 * WITH THE DEPARTMENT OF ENERGY.                                        *
 *************************************************************************/

#include "TQpDataSparse.h"

////////////////////////////////////////////////////////////////////////////////
///
/// \class TQpDataSparse
///
/// Data for the sparse QP formulation
///
////////////////////////////////////////////////////////////////////////////////

ClassImp(TQpDataSparse);

////////////////////////////////////////////////////////////////////////////////
/// Constructor

TQpDataSparse::TQpDataSparse(Int_t nx,Int_t my,Int_t mz)
: TQpDataBase(nx,my,mz)
{
   fQ.ResizeTo(fNx,fNx);
   fA.ResizeTo(fMy,fNx);
   fC.ResizeTo(fMz,fNx);
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor

TQpDataSparse::TQpDataSparse(TVectorD       &c_in,   TMatrixDSparse &Q_in,
                             TVectorD       &xlow_in,TVectorD       &ixlow_in,
                             TVectorD       &xupp_in,TVectorD       &ixupp_in,
                             TMatrixDSparse &A_in,   TVectorD       &bA_in,
                             TMatrixDSparse &C_in,
                             TVectorD       &clow_in,TVectorD       &iclow_in,
                             TVectorD       &cupp_in,TVectorD       &icupp_in)
{
   fG       .ResizeTo(c_in)    ; fG        = c_in;
   fBa      .ResizeTo(bA_in)   ; fBa       = bA_in;
   fXloBound.ResizeTo(xlow_in) ; fXloBound = xlow_in;
   fXloIndex.ResizeTo(ixlow_in); fXloIndex = ixlow_in;
   fXupBound.ResizeTo(xupp_in) ; fXupBound = xupp_in;
   fXupIndex.ResizeTo(ixupp_in); fXupIndex = ixupp_in;
   fCloBound.ResizeTo(clow_in) ; fCloBound = clow_in;
   fCloIndex.ResizeTo(iclow_in); fCloIndex = iclow_in;
   fCupBound.ResizeTo(cupp_in) ; fCupBound = cupp_in;
   fCupIndex.ResizeTo(icupp_in); fCupIndex = icupp_in;

   fNx = fG.GetNrows();
   fQ.Use(Q_in);

   if (A_in.GetNrows() > 0) {
      fA.Use(A_in);
      fMy = fA.GetNrows();
   } else
   fMy = 0;

   if (C_in.GetNrows()) {
      fC.Use(C_in);
      fMz = fC.GetNrows();
   } else
   fMz = 0;
   // fQ.Print();
   // fA.Print();
   // fC.Print();
   // printf("fNx: %d\n",fNx);
   // printf("fMy: %d\n",fMy);
   // printf("fMz: %d\n",fMz);
}


////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

TQpDataSparse::TQpDataSparse(const TQpDataSparse &another) : TQpDataBase(another)
{
   *this = another;
}


////////////////////////////////////////////////////////////////////////////////
/// Allocate space for the appropriate number of non-zeros in the matrices

void TQpDataSparse::SetNonZeros(Int_t nnzQ,Int_t nnzA,Int_t nnzC)
{
   fQ.SetSparseIndex(nnzQ);
   fA.SetSparseIndex(nnzA);
   fC.SetSparseIndex(nnzC);
}


////////////////////////////////////////////////////////////////////////////////
/// calculate y = beta*y + alpha*(fQ*x)

void TQpDataSparse::Qmult(Double_t beta,TVectorD &y,Double_t alpha,const TVectorD &x )
{
   y *= beta;
   if (fQ.GetNoElements() > 0)
      y += alpha*(fQ*x);
}


////////////////////////////////////////////////////////////////////////////////
/// calculate y = beta*y + alpha*(fA*x)

void TQpDataSparse::Amult(Double_t beta,TVectorD &y,Double_t alpha,const TVectorD &x)
{
   y *= beta;
   if (fA.GetNoElements() > 0)
      y += alpha*(fA*x);
}


////////////////////////////////////////////////////////////////////////////////
/// calculate y = beta*y + alpha*(fC*x)

void TQpDataSparse::Cmult(Double_t beta,TVectorD &y,Double_t alpha,const TVectorD &x)
{
   y *= beta;
   if (fC.GetNoElements() > 0)
      y += alpha*(fC*x);
}


////////////////////////////////////////////////////////////////////////////////
/// calculate y = beta*y + alpha*(fA^T*x)

void TQpDataSparse::ATransmult(Double_t beta,TVectorD &y,Double_t alpha,const TVectorD &x)
{
   y *= beta;
   if (fA.GetNoElements() > 0)
      y += alpha*(TMatrixDSparse(TMatrixDSparse::kTransposed,fA)*x);
}


////////////////////////////////////////////////////////////////////////////////
/// calculate y = beta*y + alpha*(fC^T*x)

void TQpDataSparse::CTransmult(Double_t beta,TVectorD &y,Double_t alpha,const TVectorD &x)
{
   y *= beta;
   if (fC.GetNoElements() > 0)
      y += alpha*(TMatrixDSparse(TMatrixDSparse::kTransposed,fC)*x);
}


////////////////////////////////////////////////////////////////////////////////
/// Return the largest component of several vectors in the data class

Double_t TQpDataSparse::DataNorm()
{
   Double_t norm = 0.0;

   Double_t componentNorm = fG.NormInf();
   if (componentNorm > norm) norm = componentNorm;

   TMatrixDSparse fQ_abs(fQ);
   componentNorm = (fQ_abs.Abs()).Max();
   if (componentNorm > norm) norm = componentNorm;

   componentNorm = fBa.NormInf();
   if (componentNorm > norm) norm = componentNorm;

   TMatrixDSparse fA_abs(fA);
   componentNorm = (fA_abs.Abs()).Max();
   if (componentNorm > norm) norm = componentNorm;

   TMatrixDSparse fC_abs(fC);
   componentNorm = (fC_abs.Abs()).Max();
   if (componentNorm > norm) norm = componentNorm;

   R__ASSERT(fXloBound.MatchesNonZeroPattern(fXloIndex));
   componentNorm = fXloBound.NormInf();
   if (componentNorm > norm) norm = componentNorm;

   R__ASSERT(fXupBound.MatchesNonZeroPattern(fXupIndex));
   componentNorm = fXupBound.NormInf();
   if (componentNorm > norm) norm = componentNorm;

   R__ASSERT(fCloBound.MatchesNonZeroPattern(fCloIndex));
   componentNorm = fCloBound.NormInf();
   if (componentNorm > norm) norm = componentNorm;

   R__ASSERT(fCupBound.MatchesNonZeroPattern(fCupIndex));
   componentNorm = fCupBound.NormInf();
   if (componentNorm > norm) norm = componentNorm;

   return norm;
}


////////////////////////////////////////////////////////////////////////////////
/// Print class members

void TQpDataSparse::Print(Option_t * /*opt*/) const
{
   fQ.Print("Q");
   fG.Print("c");

   fXloBound.Print("xlow");
   fXloIndex.Print("ixlow");

   fXupBound.Print("xupp");
   fXupIndex.Print("ixupp");

   fA.Print("A");
   fBa.Print("b");
   fC.Print("C");

   fCloBound.Print("clow");
   fCloIndex.Print("iclow");

   fCupBound.Print("cupp");
   fCupIndex.Print("icupp");
}


////////////////////////////////////////////////////////////////////////////////
/// Insert the Hessian Q into the matrix M at index (row,col) for the fundamental
/// linear system

void TQpDataSparse::PutQIntoAt(TMatrixDBase &m,Int_t row,Int_t col)
{
   m.SetSub(row,col,fQ);
}


////////////////////////////////////////////////////////////////////////////////
/// Insert the constraint matrix A into the matrix M at index (row,col) for the fundamental
/// linear system

void TQpDataSparse::PutAIntoAt(TMatrixDBase &m,Int_t row,Int_t col)
{
   m.SetSub(row,col,fA);
}


////////////////////////////////////////////////////////////////////////////////
/// Insert the constraint matrix C into the matrix M at index (row,col) for the fundamental
/// linear system

void TQpDataSparse::PutCIntoAt(TMatrixDBase &m,Int_t row,Int_t col)
{
   m.SetSub(row,col,fC);
}


////////////////////////////////////////////////////////////////////////////////
/// Return in vector dq the diagonal of matrix fQ

void TQpDataSparse::GetDiagonalOfQ(TVectorD &dq)
{
   const Int_t n = TMath::Min(fQ.GetNrows(),fQ.GetNcols());
   dq.ResizeTo(n);
   dq = TMatrixDSparseDiag(fQ);
}


////////////////////////////////////////////////////////////////////////////////
/// Return value of the objective function

Double_t TQpDataSparse::ObjectiveValue(TQpVar *vars)
{
   TVectorD tmp(fG);
   this->Qmult(1.0,tmp,0.5,vars->fX);

   return tmp*vars->fX;
}


////////////////////////////////////////////////////////////////////////////////
/// Choose randomly a QP problem

void TQpDataSparse::DataRandom(TVectorD &x,TVectorD &y,TVectorD &z,TVectorD &s)
{
   Double_t ix = 3074.20374;

   TVectorD xdual(fNx);
   this->RandomlyChooseBoundedVariables(x,xdual,fXloBound,fXloIndex,fXupBound,fXupIndex,ix,.25,.25,.25);

   TVectorD sprime(fMz);
   this->RandomlyChooseBoundedVariables(sprime,z,fCloBound,fCloIndex,fCupBound,fCupIndex,ix,.25,.25,.5);

   fQ.RandomizePD(0.0,1.0,ix);
   fA.Randomize(-10.0,10.0,ix);
   fC.Randomize(-10.0,10.0,ix);
   y .Randomize(-10.0,10.0,ix);

   // fG = - fQ x + A^T y + C^T z + xdual

   fG = xdual;
   fG -= fQ*x;

   fG += TMatrixDSparse(TMatrixDSparse::kTransposed,fA)*y;
   fG += TMatrixDSparse(TMatrixDSparse::kTransposed,fC)*z;

   fBa = fA*x;
   s   = fC*x;

   // Now compute the real q = s-sprime
   const TVectorD q = s-sprime;

   // Adjust fCloBound and fCupBound appropriately
   Add(fCloBound,1.0,q);
   Add(fCupBound,1.0,q);

   fCloBound.SelectNonZeros(fCloIndex);
   fCupBound.SelectNonZeros(fCupIndex);
}


////////////////////////////////////////////////////////////////////////////////
/// Assignment operator

TQpDataSparse &TQpDataSparse::operator=(const TQpDataSparse &source)
{
   if (this != &source) {
      TQpDataBase::operator=(source);
      fQ.ResizeTo(source.fQ); fQ = source.fQ;
      fA.ResizeTo(source.fA); fA = source.fA;
      fC.ResizeTo(source.fC); fC = source.fC;
   }
   return *this;
}
