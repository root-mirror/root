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
/// \class TQpResidual
///
/// The Residuals class calculates and stores the quantities that appear
/// on the right-hand side of the linear systems that arise at each
/// interior-point iteration. These residuals can be partitioned into
/// two fundamental categories: the components arising from the linear
/// equations in the KKT conditions, and the components arising from the
/// complementarity conditions.
///
////////////////////////////////////////////////////////////////////////////////

#include "TQpResidual.h"

ClassImp(TQpResidual);

////////////////////////////////////////////////////////////////////////////////
/// Constructor

TQpResidual::TQpResidual()
{
   fNx   = 0;
   fMy   = 0;
   fMz   = 0;

   fNxup = 0.0;
   fNxlo = 0.0;
   fMcup = 0.0;
   fMclo = 0.0;
   fResidualNorm = 0.0;
   fDualityGap = 0.0;
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor

TQpResidual::TQpResidual(Int_t nx,Int_t my,Int_t mz,TVectorD &ixlo,TVectorD &ixup,
                         TVectorD &iclo,TVectorD &icup)
{
   fNx = nx;
   fMy = my;
   fMz = mz;

   if (ixlo.GetNrows() > 0) fXloIndex.Use(ixlo.GetNrows(),ixlo.GetMatrixArray());
   if (ixup.GetNrows() > 0) fXupIndex.Use(ixup.GetNrows(),ixup.GetMatrixArray());
   if (iclo.GetNrows() > 0) fCloIndex.Use(iclo.GetNrows(),iclo.GetMatrixArray());
   if (icup.GetNrows() > 0) fCupIndex.Use(icup.GetNrows(),icup.GetMatrixArray());
   fNxlo = ixlo.NonZeros();
   fNxup = ixup.NonZeros();
   fMclo = iclo.NonZeros();
   fMcup = icup.NonZeros();

   fRQ.ResizeTo(fNx);
   fRA.ResizeTo(fMy);
   fRC.ResizeTo(fMz);

   fRz.ResizeTo(fMz);
   if (fMclo > 0) {
      fRt.ResizeTo(fMz);
      fRlambda.ResizeTo(fMz);
   }
   if (fMcup > 0) {
      fRu.ResizeTo(fMz);
      fRpi.ResizeTo(fMz);
   }
   if (fNxlo > 0) {
      fRv.ResizeTo(fNx);
      fRgamma.ResizeTo(fNx);
   }
   if (fNxup > 0) {
      fRw.ResizeTo(fNx);
      fRphi.ResizeTo(fNx);
   }

   fResidualNorm = 0.0;
   fDualityGap = 0.0;
}


////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

TQpResidual::TQpResidual(const TQpResidual &another) : TObject(another)
{
   *this = another;
}


////////////////////////////////////////////////////////////////////////////////
/// Calculate residuals, their norms, and duality complementarity gap,
/// given a problem and variable set.

void TQpResidual::CalcResids(TQpDataBase *prob_in,TQpVar *vars)
{
   TQpDataDens *prob = (TQpDataDens *) prob_in;

   fRQ.ResizeTo(prob->fG); fRQ = prob->fG;
   prob->Qmult(1.0,fRQ,1.0,vars->fX);

   // calculate x^T (g+Qx) - contribution to the duality gap
   Double_t gap = fRQ*vars->fX;

   prob->ATransmult(1.0,fRQ,-1.0,vars->fY);
   prob->CTransmult(1.0,fRQ,-1.0,vars->fZ);
   if (fNxlo > 0) Add(fRQ,-1.0,vars->fGamma);
   if (fNxup > 0) Add(fRQ, 1.0,vars->fPhi);

   Double_t norm = 0.0;
   Double_t componentNorm = fRQ.NormInf();
   if (componentNorm > norm) norm = componentNorm;

   fRA.ResizeTo(prob->fBa); fRA = prob->fBa;
   prob->Amult(-1.0,fRA,1.0,vars->fX);

   // contribution -d^T y to duality gap
   gap -= prob->fBa*vars->fY;

   componentNorm = fRA.NormInf();
   if( componentNorm > norm ) norm = componentNorm;

   fRC.ResizeTo(vars->fS); fRC = vars->fS;
   prob->Cmult(-1.0,fRC,1.0,vars->fX);

   componentNorm = fRC.NormInf();
   if( componentNorm > norm ) norm = componentNorm;

   fRz.ResizeTo(vars->fZ); fRz = vars->fZ;

   if (fMclo > 0) {
      Add(fRz,-1.0,vars->fLambda);

      fRt.ResizeTo(vars->fS); fRt = vars->fS;
      Add(fRt,-1.0,prob->GetSlowerBound());
      fRt.SelectNonZeros(fCloIndex);
      Add(fRt,-1.0,vars->fT);

      gap -= prob->fCloBound*vars->fLambda;

      componentNorm = fRt.NormInf();
      if( componentNorm > norm ) norm = componentNorm;
   }

   if (fMcup > 0) {
      Add(fRz,1.0,vars->fPi);

      fRu.ResizeTo(vars->fS); fRu = vars->fS;
      Add(fRu,-1.0,prob->GetSupperBound() );
      fRu.SelectNonZeros(fCupIndex);
      Add(fRu,1.0,vars->fU);

      gap += prob->fCupBound*vars->fPi;

      componentNorm = fRu.NormInf();
      if( componentNorm > norm ) norm = componentNorm;
   }

   componentNorm = fRz.NormInf();
   if( componentNorm > norm ) norm = componentNorm;

   if (fNxlo > 0) {
      fRv.ResizeTo(vars->fX); fRv = vars->fX;
      Add(fRv,-1.0,prob->GetXlowerBound());
      fRv.SelectNonZeros(fXloIndex);
      Add(fRv,-1.0,vars->fV);

      gap -= prob->fXloBound*vars->fGamma;

      componentNorm = fRv.NormInf();
      if( componentNorm > norm ) norm = componentNorm;
   }

   if (fNxup > 0) {
      fRw.ResizeTo(vars->fX); fRw = vars->fX;
      Add(fRw,-1.0,prob->GetXupperBound());
      fRw.SelectNonZeros(fXupIndex);
      Add(fRw,1.0,vars->fW);

      gap += prob->fXupBound*vars->fPhi;

      componentNorm = fRw.NormInf();
      if (componentNorm > norm) norm = componentNorm;
   }

   fDualityGap   = gap;
   fResidualNorm = norm;
}


////////////////////////////////////////////////////////////////////////////////
/// Modify the "complementarity" component of the residuals, by adding the pairwise
/// products of the complementary variables plus a constant alpha to this term.

void TQpResidual::Add_r3_xz_alpha(TQpVar *vars,Double_t alpha)
{
   if (fMclo > 0) AddElemMult(fRlambda,1.0,vars->fT,vars->fLambda);
   if (fMcup > 0) AddElemMult(fRpi    ,1.0,vars->fU,vars->fPi);
   if (fNxlo > 0) AddElemMult(fRgamma ,1.0,vars->fV,vars->fGamma);
   if (fNxup > 0) AddElemMult(fRphi   ,1.0,vars->fW,vars->fPhi);

   if (alpha != 0.0) {
      if (fMclo > 0) fRlambda.AddSomeConstant(alpha,fCloIndex);
      if (fMcup > 0) fRpi    .AddSomeConstant(alpha,fCupIndex);
      if (fNxlo > 0) fRgamma .AddSomeConstant(alpha,fXloIndex);
      if (fNxup > 0) fRphi   .AddSomeConstant(alpha,fXupIndex);
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Set the "complementarity" component of the residuals to the pairwise products of
/// the complementary variables plus a constant alpha .

void TQpResidual::Set_r3_xz_alpha(TQpVar *vars,Double_t alpha)
{
   this->Clear_r3();
   this->Add_r3_xz_alpha(vars,alpha);
}


////////////////////////////////////////////////////////////////////////////////
/// set the complementarity component of the residuals to 0.

void TQpResidual::Clear_r3()
{
   if (fMclo > 0) fRlambda.Zero();
   if (fMcup > 0) fRpi    .Zero();
   if (fNxlo > 0) fRgamma .Zero();
   if (fNxup > 0) fRphi   .Zero();
}


////////////////////////////////////////////////////////////////////////////////
/// set the noncomplementarity components of the residual (the terms arising from
/// the linear equalities in the KKT conditions) to 0.

void TQpResidual::Clear_r1r2()
{
   fRQ.Zero();
   fRA.Zero();
   fRC.Zero();
   fRz.Zero();
   if (fNxlo > 0) fRv.Zero();
   if (fNxup > 0) fRw.Zero();
   if (fMclo > 0) fRt.Zero();
   if (fMcup > 0) fRu.Zero();
}


////////////////////////////////////////////////////////////////////////////////
/// Perform the projection operation required by Gondzio algorithm: replace each
/// component r3_i of the complementarity component of the residuals by r3p_i-r3_i,
/// where r3p_i is the projection of r3_i onto the box [rmin, rmax]. Then if the
/// resulting value is less than -rmax, replace it by -rmax.

void TQpResidual::Project_r3(Double_t rmin,Double_t rmax)
{
   if (fMclo > 0) {
      GondzioProjection(fRlambda,rmin,rmax);
      fRlambda.SelectNonZeros(fCloIndex);
   }
   if (fMcup > 0) {
      GondzioProjection(fRpi,rmin,rmax);
      fRpi.SelectNonZeros(fCupIndex);
   }
   if (fNxlo > 0) {
      GondzioProjection(fRgamma,rmin,rmax);
      fRgamma.SelectNonZeros(fXloIndex);
   }
   if (fNxup > 0) {
      GondzioProjection(fRphi,rmin,rmax);
      fRphi.SelectNonZeros(fXupIndex);
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Check if vector elements as selected through array indices are non-zero

Bool_t TQpResidual::ValidNonZeroPattern()
{
   if (fNxlo > 0 &&
      (!fRv    .MatchesNonZeroPattern(fXloIndex) ||
       !fRgamma.MatchesNonZeroPattern(fXloIndex)) ) {
      return kFALSE;
   }

   if (fNxup > 0 &&
      (!fRw  .MatchesNonZeroPattern(fXupIndex) ||
       !fRphi.MatchesNonZeroPattern(fXupIndex)) ) {
      return kFALSE;
   }
   if (fMclo > 0 &&
      (!fRt     .MatchesNonZeroPattern(fCloIndex) ||
       !fRlambda.MatchesNonZeroPattern(fCloIndex)) ) {
      return kFALSE;
   }

   if (fMcup > 0 &&
      (!fRu .MatchesNonZeroPattern(fCupIndex) ||
       !fRpi.MatchesNonZeroPattern(fCupIndex)) ) {
      return kFALSE;
   }

   return kTRUE;
}


////////////////////////////////////////////////////////////////////////////////
/// Replace each component r3_i of the complementarity component of the residuals
/// by r3p_i-r3_i, where r3p_i is the projection of r3_i onto the box [rmin, rmax].
/// Then if the resulting value is less than -rmax, replace it by -rmax.

void TQpResidual::GondzioProjection(TVectorD &v,Double_t rmin,Double_t rmax)
{
         Double_t *        ep = v.GetMatrixArray();
   const Double_t * const fep = ep+v.GetNrows();

   while (ep < fep) {
      if (*ep < rmin)
         *ep = rmin - *ep;
      else if (*ep > rmax)
         *ep = rmax - *ep;
      else
         *ep = 0.0;
      if (*ep < -rmax) *ep = -rmax;
      ep++;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Assignment operator

TQpResidual &TQpResidual::operator=(const TQpResidual &source)
{
   if (this != &source) {
      TObject::operator=(source);

      fNx   = source.fNx;
      fMy   = source.fMy;
      fMz   = source.fMz;

      fNxup = source.fNxup;
      fNxlo = source.fNxlo;
      fMcup = source.fMcup;
      fMclo = source.fMclo;

      fXupIndex.ResizeTo(source.fXupIndex); fXupIndex = source.fXupIndex;
      fXloIndex.ResizeTo(source.fXloIndex); fXloIndex = source.fXloIndex;
      fCupIndex.ResizeTo(source.fCupIndex); fCupIndex = source.fCupIndex;
      fCloIndex.ResizeTo(source.fCloIndex); fCupIndex = source.fCupIndex;

      fRQ     .ResizeTo(source.fRQ);      fRQ      = source.fRQ;
      fRA     .ResizeTo(source.fRA);      fRA      = source.fRA;
      fRC     .ResizeTo(source.fRC);      fRC      = source.fRC;
      fRz     .ResizeTo(source.fRz);      fRz      = source.fRz;
      fRv     .ResizeTo(source.fRv);      fRv      = source.fRv;
      fRw     .ResizeTo(source.fRw);      fRw      = source.fRw;
      fRt     .ResizeTo(source.fRt);      fRt      = source.fRt;
      fRu     .ResizeTo(source.fRu);      fRu      = source.fRu;
      fRgamma .ResizeTo(source.fRgamma);  fRgamma  = source.fRgamma;
      fRphi   .ResizeTo(source.fRphi);    fRphi    = source.fRphi;
      fRlambda.ResizeTo(source.fRlambda); fRlambda = source.fRlambda;
      fRpi    .ResizeTo(source.fRpi);     fRpi     = source.fRpi;

      // LM: copy also these data members
      fResidualNorm = source.fResidualNorm;
      fDualityGap = source.fDualityGap;
   }
   return *this;
}
