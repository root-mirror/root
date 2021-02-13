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

////////////////////////////////////////////////////////////////////////////////
///
/// \class TQpLinSolverDens
///
/// Implements the aspects of the solvers for dense general QP
/// formulation that are specific to the dense case.
///
////////////////////////////////////////////////////////////////////////////////

#include "TQpLinSolverDens.h"

ClassImp(TQpLinSolverDens);

////////////////////////////////////////////////////////////////////////////////
/// Constructor

TQpLinSolverDens::TQpLinSolverDens(TQpProbDens *factory,TQpDataDens *data) :
                  TQpLinSolverBase(factory,data)
{
   const Int_t n = factory->fNx+factory->fMy+factory->fMz;
   fKkt.ResizeTo(n,n);

   data->PutQIntoAt(fKkt,0,      0);
   if (fMy > 0) data->PutAIntoAt(fKkt,fNx,    0);
   if (fMz > 0) data->PutCIntoAt(fKkt,fNx+fMy,0);
   for (Int_t ix = fNx; ix < fNx+fMy+fMz; ix++) {
      for (Int_t iy = fNx; iy < fNx+fMy+fMz; iy++)
         fKkt(ix,iy) = 0.0;
   }

   fSolveLU = TDecompLU(n);
}


////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

TQpLinSolverDens::TQpLinSolverDens(const TQpLinSolverDens &another) : TQpLinSolverBase(another)
{
   *this = another;
}


////////////////////////////////////////////////////////////////////////////////
/// Sets up the matrix for the main linear system in "augmented system" form.

void TQpLinSolverDens::Factor(TQpDataBase *prob,TQpVar *vars)
{
   TQpLinSolverBase::Factor(prob,vars);
   fSolveLU.SetMatrix(fKkt);
}


////////////////////////////////////////////////////////////////////////////////
/// Places the diagonal resulting from the bounds on x into the augmented system matrix

void TQpLinSolverDens::PutXDiagonal(TVectorD &xdiag)
{
   TMatrixDDiag diag(fKkt);
   for (Int_t i = 0; i < xdiag.GetNrows(); i++)
      diag[i] = xdiag[i];
}


////////////////////////////////////////////////////////////////////////////////
/// Places the diagonal resulting from the bounds on Cx into the augmented system matrix

void TQpLinSolverDens::PutZDiagonal(TVectorD &zdiag)
{
   TMatrixDDiag diag(fKkt);
   for (Int_t i = 0; i < zdiag.GetNrows(); i++)
      diag[i+fNx+fMy] = zdiag[i];
}


////////////////////////////////////////////////////////////////////////////////
/// Perform the actual solve using the factors produced in factor.
/// rhs on input contains the aggregated right-hand side of the augmented system;
///  on output contains the solution in aggregated form .

void TQpLinSolverDens::SolveCompressed(TVectorD &compressedRhs)
{
   fSolveLU.Solve(compressedRhs);
}


////////////////////////////////////////////////////////////////////////////////
/// Assignment operator

TQpLinSolverDens &TQpLinSolverDens::operator=(const TQpLinSolverDens &source)
{
   if (this != &source) {
      TQpLinSolverBase::operator=(source);
      fKkt.ResizeTo(source.fKkt); fKkt = source.fKkt;
      fSolveLU = source.fSolveLU;
   }
   return *this;
}
