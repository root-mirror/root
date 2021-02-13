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

/** \class TQpVar

Class containing the variables for the general QP formulation

In terms of in our abstract problem formulation, these variables are
the vectors x, y, z and s.
*/

#include <iostream>
#include "TQpVar.h"

ClassImp(TQpVar);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TQpVar::TQpVar()
{
   fNx   = 0;
   fMy   = 0;
   fMz   = 0;
   fNxup = 0;
   fNxlo = 0;
   fMcup = 0;
   fMclo = 0;
   fNComplementaryVariables = 0;
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor

TQpVar::TQpVar(TVectorD &x_in,TVectorD &s_in,TVectorD &y_in,TVectorD &z_in,
               TVectorD &v_in,TVectorD &gamma_in,TVectorD &w_in,TVectorD &phi_in,
               TVectorD &t_in,TVectorD &lambda_in,TVectorD &u_in,TVectorD &pi_in,
               TVectorD &ixlow_in,TVectorD &ixupp_in,TVectorD &iclow_in,TVectorD &icupp_in)
{
   if (x_in     .GetNrows() > 0) fX.       Use(x_in     .GetNrows(),x_in     .GetMatrixArray());
   if (s_in     .GetNrows() > 0) fS.       Use(s_in     .GetNrows(),s_in     .GetMatrixArray());
   if (y_in     .GetNrows() > 0) fY.       Use(y_in     .GetNrows(),y_in     .GetMatrixArray());
   if (z_in     .GetNrows() > 0) fZ.       Use(z_in     .GetNrows(),z_in     .GetMatrixArray());
   if (v_in     .GetNrows() > 0) fV.       Use(v_in     .GetNrows(),v_in     .GetMatrixArray());
   if (phi_in   .GetNrows() > 0) fPhi.     Use(phi_in   .GetNrows(),phi_in   .GetMatrixArray());
   if (w_in     .GetNrows() > 0) fW.       Use(w_in     .GetNrows(),w_in     .GetMatrixArray());
   if (gamma_in .GetNrows() > 0) fGamma.   Use(gamma_in .GetNrows(),gamma_in .GetMatrixArray());
   if (t_in     .GetNrows() > 0) fT.       Use(t_in     .GetNrows(),t_in     .GetMatrixArray());
   if (lambda_in.GetNrows() > 0) fLambda.  Use(lambda_in.GetNrows(),lambda_in.GetMatrixArray());
   if (u_in     .GetNrows() > 0) fU.       Use(u_in     .GetNrows(),u_in     .GetMatrixArray());
   if (pi_in    .GetNrows() > 0) fPi.      Use(pi_in    .GetNrows(),pi_in    .GetMatrixArray());
   if (ixlow_in .GetNrows() > 0) fXloIndex.Use(ixlow_in .GetNrows(),ixlow_in .GetMatrixArray());
   if (ixupp_in .GetNrows() > 0) fXupIndex.Use(ixupp_in .GetNrows(),ixupp_in .GetMatrixArray());
   if (iclow_in .GetNrows() > 0) fCloIndex.Use(iclow_in .GetNrows(),iclow_in .GetMatrixArray());
   if (icupp_in .GetNrows() > 0) fCupIndex.Use(icupp_in .GetNrows(),icupp_in .GetMatrixArray());

   fNx = fX.GetNrows();
   fMy = fY.GetNrows();
   fMz = fZ.GetNrows();

   R__ASSERT(fNx == fXloIndex.GetNrows() || 0 == fXloIndex.GetNrows());
   R__ASSERT(fNx == fXloIndex.GetNrows() || 0 == fXloIndex.GetNrows());
   R__ASSERT(fMz == fCloIndex.GetNrows() || 0 == fCloIndex.GetNrows());
   R__ASSERT(fMz == fCupIndex.GetNrows() || 0 == fCupIndex.GetNrows());

   fNxlo = fXloIndex.NonZeros();
   fNxup = fXupIndex.NonZeros();
   fMclo = fCloIndex.NonZeros();
   fMcup = fCupIndex.NonZeros();
   fNComplementaryVariables = fMclo+fMcup+fNxlo+fNxup;

   R__ASSERT(fMz == fS.GetNrows());
   R__ASSERT(fNx == fV     .GetNrows() || (0 == fV     .GetNrows() && fNxlo == 0));
   R__ASSERT(fNx == fGamma .GetNrows() || (0 == fGamma .GetNrows() && fNxlo == 0));

   R__ASSERT(fNx == fW     .GetNrows() || (0 == fW     .GetNrows() && fNxup == 0));
   R__ASSERT(fNx == fPhi   .GetNrows() || (0 == fPhi   .GetNrows() && fNxup == 0));

   R__ASSERT(fMz == fT     .GetNrows() || (0 == fT     .GetNrows() && fMclo == 0));
   R__ASSERT(fMz == fLambda.GetNrows() || (0 == fLambda.GetNrows() && fMclo == 0));

   R__ASSERT(fMz == fU     .GetNrows() || (0 == fU     .GetNrows() && fMcup == 0));
   R__ASSERT(fMz == fPi    .GetNrows() || (0 == fPi    .GetNrows() && fMcup == 0));
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor

TQpVar::TQpVar(Int_t nx,Int_t my,Int_t mz,TVectorD &ixlow,TVectorD &ixupp,
               TVectorD &iclow,TVectorD &icupp)
{
   R__ASSERT(nx == ixlow.GetNrows() || 0 == ixlow.GetNrows());
   R__ASSERT(nx == ixlow.GetNrows() || 0 == ixlow.GetNrows());
   R__ASSERT(mz == iclow.GetNrows() || 0 == iclow.GetNrows());
   R__ASSERT(mz == icupp.GetNrows() || 0 == icupp.GetNrows());

   fNxlo = ixlow.NonZeros();
   fNxup = ixupp.NonZeros();
   fMclo = iclow.NonZeros();
   fMcup = icupp.NonZeros();

   if (ixlow.GetNrows() > 0) fXloIndex.Use(ixlow.GetNrows(),ixlow.GetMatrixArray());
   if (ixupp.GetNrows() > 0) fXupIndex.Use(ixupp.GetNrows(),ixupp.GetMatrixArray());
   if (iclow.GetNrows() > 0) fCloIndex.Use(iclow.GetNrows(),iclow.GetMatrixArray());
   if (icupp.GetNrows() > 0) fCupIndex.Use(icupp.GetNrows(),icupp.GetMatrixArray());

   fNx = nx;
   fMy = my;
   fMz = mz;

   if (fMclo > 0) {
      fT.ResizeTo(fMz);
      fLambda.ResizeTo(fMz);
   }
   if (fMcup > 0) {
      fU.ResizeTo(fMz);
      fPi.ResizeTo(fMz);
   }
   if (fNxlo > 0) {
      fV.ResizeTo(fNx);
      fGamma.ResizeTo(fNx);
   }

   if (fNxup > 0) {
      fW.ResizeTo(fNx);
      fPhi.ResizeTo(fNx);
   }

   fS.ResizeTo(fMz);
   fX.ResizeTo(fNx);
   fY.ResizeTo(fMy);
   fZ.ResizeTo(fMz);
   fNComplementaryVariables = fMclo+fMcup+fNxlo+fNxup;
}


////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

TQpVar::TQpVar(const TQpVar &another) : TObject(another)
{
   *this = another;
}


////////////////////////////////////////////////////////////////////////////////
/// compute complementarity gap, obtained by taking the inner product of the
/// complementary vectors and dividing by the total number of components
/// computes mu = (t'lambda +u'pi + v'gamma + w'phi)/(mclow+mcupp+nxlow+nxupp)

Double_t TQpVar::GetMu()
{
   Double_t mu = 0.0;
   if (fNComplementaryVariables > 0 ) {
      if (fMclo > 0) mu += fT*fLambda;
      if (fMcup > 0) mu += fU*fPi;
      if (fNxlo > 0) mu += fV*fGamma;
      if (fNxup > 0) mu += fW*fPhi;

      mu /= fNComplementaryVariables;
   }
   return mu;
}


////////////////////////////////////////////////////////////////////////////////
/// Compute the complementarity gap resulting from a step of length "alpha" along
/// direction "step"

Double_t TQpVar::MuStep(TQpVar *step,Double_t alpha)
{
   Double_t mu = 0.0;
   if (fNComplementaryVariables > 0) {
      if (fMclo > 0)
         mu += (fT+alpha*step->fT)*(fLambda+alpha*step->fLambda);
      if (fMcup > 0)
         mu += (fU+alpha*step->fU)*(fPi+alpha*step->fPi);
      if (fNxlo > 0)
         mu += (fV+alpha*step->fV)*(fGamma+alpha*step->fGamma);
      if (fNxup > 0)
         mu += (fW+alpha*step->fW)*(fPhi+alpha*step->fPhi);
      mu /= fNComplementaryVariables;
   }
   return mu;
}


////////////////////////////////////////////////////////////////////////////////
/// Perform a "saxpy" operation on all data vectors : x += alpha*y

void TQpVar::Saxpy(TQpVar *b,Double_t alpha)
{
   Add(fX,alpha,b->fX);
   Add(fY,alpha,b->fY);
   Add(fZ,alpha,b->fZ);
   Add(fS,alpha,b->fS);
   if (fMclo > 0) {
      R__ASSERT((b->fT)     .MatchesNonZeroPattern(fCloIndex) &&
         (b->fLambda).MatchesNonZeroPattern(fCloIndex));

      Add(fT     ,alpha,b->fT);
      Add(fLambda,alpha,b->fLambda);
   }
   if (fMcup > 0) {
      R__ASSERT((b->fU) .MatchesNonZeroPattern(fCupIndex) &&
         (b->fPi).MatchesNonZeroPattern(fCupIndex));

      Add(fU ,alpha,b->fU);
      Add(fPi,alpha,b->fPi);
   }
   if (fNxlo > 0) {
      R__ASSERT((b->fV)    .MatchesNonZeroPattern(fXloIndex) &&
         (b->fGamma).MatchesNonZeroPattern(fXloIndex));

      Add(fV    ,alpha,b->fV);
      Add(fGamma,alpha,b->fGamma);
   }
   if (fNxup > 0) {
      R__ASSERT((b->fW)  .MatchesNonZeroPattern(fXupIndex) &&
         (b->fPhi).MatchesNonZeroPattern(fXupIndex));

      Add(fW  ,alpha,b->fW);
      Add(fPhi,alpha,b->fPhi);
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Perform a "negate" operation on all data vectors : x =  -x

void TQpVar::Negate()
{
   fS *= -1.;
   fX *= -1.;
   fY *= -1.;
   fZ *= -1.;
   if (fMclo > 0) {
      fT      *= -1.;
      fLambda *= -1.;
   }
   if (fMcup > 0) {
      fU  *= -1.;
      fPi *= -1.;
   }
   if (fNxlo > 0) {
      fV     *= -1.;
      fGamma *= -1.;
   }
   if (fNxup > 0) {
      fW   *= -1.;
      fPhi *= -1.;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// calculate the largest alpha in (0,1] such that the/ nonnegative variables stay
/// nonnegative in the given search direction. In the general QP problem formulation
/// this is the largest value of alpha such that
///     (t,u,v,w,lambda,pi,phi,gamma) + alpha * (b->t,b->u,b->v,b->w,b->lambda,b->pi,
///                                                b->phi,b->gamma) >= 0.

Double_t TQpVar::StepBound(TQpVar *b)
{
   Double_t maxStep = 1.0;

   if (fMclo > 0 ) {
      R__ASSERT(fT     .SomePositive(fCloIndex));
      R__ASSERT(fLambda.SomePositive(fCloIndex));

      maxStep = this->StepBound(fT,     b->fT,     maxStep);
      maxStep = this->StepBound(fLambda,b->fLambda,maxStep);
   }

   if (fMcup > 0 ) {
      R__ASSERT(fU .SomePositive(fCupIndex));
      R__ASSERT(fPi.SomePositive(fCupIndex));

      maxStep = this->StepBound(fU, b->fU, maxStep);
      maxStep = this->StepBound(fPi,b->fPi,maxStep);
   }

   if (fNxlo > 0 ) {
      R__ASSERT(fV    .SomePositive(fXloIndex));
      R__ASSERT(fGamma.SomePositive(fXloIndex));

      maxStep = this->StepBound(fV,    b->fV,    maxStep);
      maxStep = this->StepBound(fGamma,b->fGamma,maxStep);
   }

   if (fNxup > 0 ) {
      R__ASSERT(fW  .SomePositive(fXupIndex));
      R__ASSERT(fPhi.SomePositive(fXupIndex));

      maxStep = this->StepBound(fW,  b->fW,  maxStep);
      maxStep = this->StepBound(fPhi,b->fPhi,maxStep);
   }

   return maxStep;
}


////////////////////////////////////////////////////////////////////////////////
/// Find the maximum stepsize of v in direction dir
/// before violating the nonnegativity constraints

Double_t TQpVar::StepBound(TVectorD &v,TVectorD &dir,Double_t maxStep)
{
   if (!AreCompatible(v,dir)) {
      ::Error("StepBound(TVectorD &,TVectorD &,Double_t)","vector's not compatible");
      return kFALSE;
   }

   const Int_t n = v.GetNrows();
   const Double_t * const pD = dir.GetMatrixArray();
   const Double_t * const pV = v.GetMatrixArray();

   Double_t bound = maxStep;
   for (Int_t i = 0; i < n; i++) {
      Double_t tmp = pD[i];
      if ( pV[i] >= 0 && tmp < 0 ) {
         tmp = -pV[i]/tmp;
         if (tmp < bound)
            bound = tmp;
      }
   }
   return bound;
}


////////////////////////////////////////////////////////////////////////////////
/// Is the current position an interior point  ?

Bool_t TQpVar::IsInteriorPoint()
{
   Bool_t interior = kTRUE;
   if (fMclo > 0)
      interior = interior &&
         fT.SomePositive(fCloIndex) && fLambda.SomePositive(fCloIndex);

   if (fMcup > 0)
      interior = interior &&
         fU.SomePositive(fCupIndex) && fPi.SomePositive(fCupIndex);

   if (fNxlo > 0)
      interior = interior &&
         fV.SomePositive(fXloIndex) && fGamma.SomePositive(fXloIndex);

   if (fNxup > 0)
      interior = interior &&
         fW.SomePositive(fXupIndex) && fPhi.SomePositive(fXupIndex);

   return interior;
}


////////////////////////////////////////////////////////////////////////////////
/// Performs the same function as StepBound, and supplies additional information about
/// which component of the nonnegative variables is responsible for restricting alpha.
/// In terms of the abstract formulation, the components have the following meanings :
///
///  primalValue   : the value of the blocking component of the primal variables (u,t,v,w).
///  primalStep    : the corresponding value of the blocking component of the primal step
///                  variables (b->u,b->t,b->v,b->w)
///  dualValue     : the value of the blocking component of the dual variables/
///                  (lambda,pi,phi,gamma).
///  dualStep      : the corresponding value of the blocking component of the dual step
///                   variables (b->lambda,b->pi,b->phi,b->gamma)
///  firstOrSecond : 1 if the primal step is blocking,
///                  2 if the dual step is block,
///                  0 if no step is blocking.

Double_t TQpVar::FindBlocking(TQpVar   *step,
                              Double_t &primalValue,
                              Double_t &primalStep,
                              Double_t &dualValue,
                              Double_t &dualStep,
                              Int_t    &fIrstOrSecond)
{
   fIrstOrSecond = 0;
   Double_t alpha = 1.0;
   if (fMclo > 0)
      alpha = FindBlocking(fT,step->fT,fLambda,step->fLambda,alpha,
         primalValue,primalStep,dualValue,dualStep,fIrstOrSecond);

   if (fMcup > 0)
      alpha = FindBlocking(fU,step->fU,fPi,step->fPi,alpha,
         primalValue,primalStep,dualValue,dualStep,fIrstOrSecond);

   if (fNxlo > 0)
      alpha = FindBlocking(fV,step->fV,fGamma,step->fGamma,alpha,
         primalValue,primalStep,dualValue,dualStep,fIrstOrSecond);

   if (fNxup > 0)
      alpha = FindBlocking(fW,step->fW,fPhi,step->fPhi,alpha,
         primalValue,primalStep,dualValue,dualStep,fIrstOrSecond);

   return alpha;
}


////////////////////////////////////////////////////////////////////////////////
/// See other FindBlocking function

Double_t TQpVar::FindBlocking(TVectorD &w,TVectorD &wstep,TVectorD &u,TVectorD &ustep,
                              Double_t maxStep,Double_t &w_elt,Double_t &wstep_elt,Double_t &u_elt,
                              Double_t &ustep_elt,int& fIrst_or_second)
{
   return FindBlockingSub(w.GetNrows(),
      w.GetMatrixArray(),    1,
      wstep.GetMatrixArray(),1,
      u.GetMatrixArray(),    1,
      ustep.GetMatrixArray(),1,
      maxStep,
      w_elt,wstep_elt,
      u_elt,ustep_elt,
      fIrst_or_second);
}


////////////////////////////////////////////////////////////////////////////////
/// See FindBlocking function

Double_t TQpVar::FindBlockingSub(Int_t n,
                                 Double_t *w,    Int_t incw,
                                 Double_t *wstep,Int_t incwstep,
                                 Double_t *u,    Int_t incu,
                                 Double_t *ustep,Int_t incustep,
                                 Double_t maxStep,
                                 Double_t &w_elt,Double_t &wstep_elt,
                                 Double_t &u_elt,Double_t &ustep_elt,
                                 Int_t &fIrst_or_second)
{
   Double_t bound = maxStep;

   Int_t i = n-1;
   Int_t lastBlocking = -1;

   // Search backward so that we fInd the blocking constraint of lowest
   // index. We do this to make things consistent with MPI's MPI_MINLOC,
   // which returns the processor with smallest rank where a min occurs.
   //
   // Still, going backward is ugly!
   Double_t *pw     = w    +(n-1)*incw;
   Double_t *pwstep = wstep+(n-1)*incwstep;
   Double_t *pu     = u    +(n-1)*incu;
   Double_t *pustep = ustep+(n-1)*incustep;

   while (i >= 0) {
      Double_t tmp = *pwstep;
      if (*pw > 0 && tmp < 0) {
         tmp = -*pw/tmp;
         if( tmp <= bound ) {
            bound = tmp;
            lastBlocking = i;
            fIrst_or_second = 1;
         }
      }
      tmp = *pustep;
      if (*pu > 0 && tmp < 0) {
         tmp = -*pu/tmp;
         if( tmp <= bound ) {
            bound = tmp;
            lastBlocking = i;
            fIrst_or_second = 2;
         }
      }

      i--;
      if (i >= 0) {
         // It is safe to decrement the pointers
         pw     -= incw;
         pwstep -= incwstep;
         pu     -= incu;
         pustep -= incustep;
      }
   }

   if (lastBlocking > -1) {
      // fIll out the elements
      w_elt     = w[lastBlocking];
      wstep_elt = wstep[lastBlocking];
      u_elt     = u[lastBlocking];
      ustep_elt = ustep[lastBlocking];
   }
   return bound;
}


////////////////////////////////////////////////////////////////////////////////
/// Sets components of (u,t,v,w) to alpha and of (lambda,pi,phi,gamma) to beta

void TQpVar::InteriorPoint(Double_t alpha,Double_t beta)
{
   fS.Zero();
   fX.Zero();
   fY.Zero();
   fZ.Zero();

   if (fNxlo > 0) {
      fV = alpha;
      fV.SelectNonZeros(fXloIndex);
      fGamma = beta;
      fGamma.SelectNonZeros(fXloIndex);
   }

   if (fNxup > 0) {
      fW = alpha;
      fW.SelectNonZeros(fXupIndex);
      fPhi = beta;
      fPhi.SelectNonZeros(fXupIndex);
   }

   if (fMclo > 0 ) {
      fT = alpha;
      fT.SelectNonZeros(fCloIndex);
      fLambda = beta;
      fLambda.SelectNonZeros(fCloIndex);
   }

   if (fMcup > 0) {
      fU = alpha;
      fU.SelectNonZeros(fCupIndex);
      fPi = beta;
      fPi.SelectNonZeros(fCupIndex);
   }
}


////////////////////////////////////////////////////////////////////////////////
/// The amount by which the current variables violate the  non-negativity constraints.

Double_t TQpVar::Violation()
{
   Double_t viol = 0.0;
   Double_t cmin;

   if (fNxlo > 0) {
      cmin = fV.Min();
      if (cmin < viol) viol = cmin;

      cmin = fGamma.Min();
      if (cmin < viol) viol = cmin;
   }
   if (fNxup > 0) {
      cmin = fW.Min();
      if (cmin < viol) viol = cmin;

      cmin = fPhi.Min();
      if (cmin < viol) viol = cmin;
   }
   if (fMclo > 0) {
      cmin = fT.Min();
      if (cmin < viol) viol = cmin;

      cmin = fLambda.Min();
      if (cmin < viol) viol = cmin;
   }
   if (fMcup > 0) {
      cmin = fU.Min();
      if (cmin < viol) viol = cmin;

      cmin = fPi.Min();
      if (cmin < viol) viol = cmin;
   }

   return -viol;
}


////////////////////////////////////////////////////////////////////////////////
/// Add alpha to components of (u,t,v,w) and beta to components of (lambda,pi,phi,gamma)

void TQpVar::ShiftBoundVariables(Double_t alpha,Double_t beta)
{
   if (fNxlo > 0) {
      fV    .AddSomeConstant(alpha,fXloIndex);
      fGamma.AddSomeConstant(beta, fXloIndex);
   }
   if (fNxup > 0) {
      fW  .AddSomeConstant(alpha,fXupIndex);
      fPhi.AddSomeConstant(beta, fXupIndex);
   }
   if (fMclo > 0) {
      fT     .AddSomeConstant(alpha,fCloIndex);
      fLambda.AddSomeConstant(beta, fCloIndex);
   }
   if (fMcup > 0) {
      fU .AddSomeConstant(alpha,fCupIndex);
      fPi.AddSomeConstant(beta, fCupIndex);
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Print class members

void TQpVar::Print(Option_t * /*option*/) const
{
   std::cout << "fNx  : " << fNx   << std::endl;
   std::cout << "fMy  : " << fMy   << std::endl;
   std::cout << "fMz  : " << fMz   << std::endl;
   std::cout << "fNxup: " << fNxup << std::endl;
   std::cout << "fNxlo: " << fNxlo << std::endl;
   std::cout << "fMcup: " << fMcup << std::endl;
   std::cout << "fMclo: " << fMclo << std::endl;

   fXloIndex.Print("fXloIndex");
   fXupIndex.Print("fXupIndex");
   fCupIndex.Print("fCupIndex");
   fCloIndex.Print("fCloIndex");

   fX.Print("fX");
   fS.Print("fS");
   fY.Print("fY");
   fZ.Print("fZ");

   fV.Print("fV");
   fPhi.Print("fPhi");

   fW.Print("fW");
   fGamma.Print("fGamma");

   fT.Print("fT");
   fLambda.Print("fLambda");

   fU.Print("fU");
   fPi.Print("fPi");
}


////////////////////////////////////////////////////////////////////////////////
/// Return the sum of the vector-norm1's

Double_t TQpVar::Norm1()
{
   Double_t norm = 0.0;
   norm += fX.Norm1();
   norm += fS.Norm1();
   norm += fY.Norm1();
   norm += fZ.Norm1();

   norm += fV.Norm1();
   norm += fPhi.Norm1();
   norm += fW.Norm1();
   norm += fGamma.Norm1();
   norm += fT.Norm1();
   norm += fLambda.Norm1();
   norm += fU.Norm1();
   norm += fPi.Norm1();

   return norm;
}


////////////////////////////////////////////////////////////////////////////////
/// Return the sum of the vector-normInf's

Double_t TQpVar::NormInf()
{
   Double_t norm = 0.0;

   Double_t tmp = fX.NormInf();
   if (tmp > norm) norm = tmp;
   tmp = fS.NormInf();
   if (tmp > norm) norm = tmp;
   tmp = fY.NormInf();
   if (tmp > norm) norm = tmp;
   tmp = fZ.NormInf();
   if (tmp > norm) norm = tmp;

   tmp = fV.NormInf();
   if (tmp > norm) norm = tmp;
   tmp = fPhi.NormInf();
   if (tmp > norm) norm = tmp;

   tmp = fW.NormInf();
   if (tmp > norm) norm = tmp;
   tmp = fGamma.NormInf();
   if (tmp > norm) norm = tmp;

   tmp = fT.NormInf();
   if (tmp > norm) norm = tmp;
   tmp = fLambda.NormInf();
   if (tmp > norm) norm = tmp;

   tmp = fU.NormInf();
   if (tmp > norm) norm = tmp;
   tmp = fPi.NormInf();
   if (tmp > norm) norm = tmp;

   return norm;
}


////////////////////////////////////////////////////////////////////////////////
/// Check that the variables conform to the non-zero indices

Bool_t TQpVar::ValidNonZeroPattern()
{
   if (fNxlo > 0 &&
      ( !fV    .MatchesNonZeroPattern(fXloIndex) ||
        !fGamma.MatchesNonZeroPattern(fXloIndex) ) ) {
      return kFALSE;
   }

   if (fNxup > 0 &&
      ( !fW  .MatchesNonZeroPattern(fXupIndex) ||
        !fPhi.MatchesNonZeroPattern(fXupIndex) ) ) {
      return kFALSE;
   }
   if (fMclo > 0 &&
      ( !fT     .MatchesNonZeroPattern(fCloIndex) ||
        !fLambda.MatchesNonZeroPattern(fCloIndex) ) ) {
      return kFALSE;
   }

   if (fMcup > 0 &&
      ( !fU .MatchesNonZeroPattern(fCupIndex) ||
        !fPi.MatchesNonZeroPattern(fCupIndex) ) ) {
      return kFALSE;
   }

   return kTRUE;
}


////////////////////////////////////////////////////////////////////////////////
/// Assignment operator

TQpVar &TQpVar::operator=(const TQpVar &source)
{
   if (this != &source) {
      TObject::operator=(source);
      fNx       = source.fNx;
      fMy       = source.fMy;
      fMz       = source.fMz;
      fNxup     = source.fNxup;
      fNxlo     = source.fNxlo;
      fMcup     = source.fMcup;
      fMclo     = source.fMclo;

      fXloIndex = source.fXloIndex;
      fXupIndex = source.fXupIndex;
      fCupIndex = source.fCupIndex;
      fCloIndex = source.fCloIndex;

      fX     .ResizeTo(source.fX);      fX      = source.fX;
      fS     .ResizeTo(source.fS);      fS      = source.fS;
      fY     .ResizeTo(source.fY);      fY      = source.fY;
      fZ     .ResizeTo(source.fZ);      fZ      = source.fZ;

      fV     .ResizeTo(source.fV);      fV      = source.fV;
      fPhi   .ResizeTo(source.fPhi);    fPhi    = source.fPhi;

      fW     .ResizeTo(source.fW);      fW      = source.fW;
      fGamma .ResizeTo(source.fGamma) ; fGamma  = source.fGamma;

      fT     .ResizeTo(source.fT);      fT      = source.fT;
      fLambda.ResizeTo(source.fLambda); fLambda = source.fLambda;

      fU     .ResizeTo(source.fU);      fU      = source.fU;
      fPi    .ResizeTo(source.fPi);     fPi     = source.fPi;

      // LM: copy also this data member
      fNComplementaryVariables = source.fNComplementaryVariables;
   }
   return *this;
}
