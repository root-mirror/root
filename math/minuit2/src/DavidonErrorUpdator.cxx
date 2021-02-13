// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "Minuit2/DavidonErrorUpdator.h"
#include "Minuit2/MinimumState.h"
#include "Minuit2/LaSum.h"
#include "Minuit2/LaProd.h"
#include "Minuit2/MnPrint.h"

namespace ROOT {

namespace Minuit2 {

double inner_product(const LAVector &, const LAVector &);
double similarity(const LAVector &, const LASymMatrix &);
double sum_of_elements(const LASymMatrix &);

MinimumError
DavidonErrorUpdator::Update(const MinimumState &s0, const MinimumParameters &p1, const FunctionGradient &g1) const
{

   // update of the covarianze matrix (Davidon formula, see Tutorial, par. 4.8 pag 26)
   // in case of delgam > gvg (PHI > 1) use rank one formula
   // see  par 4.10 pag 30

   MnPrint print("DavidonErrorUpdator");

   const MnAlgebraicSymMatrix &v0 = s0.Error().InvHessian();
   MnAlgebraicVector dx = p1.Vec() - s0.Vec();
   MnAlgebraicVector dg = g1.Vec() - s0.Gradient().Vec();

   double delgam = inner_product(dx, dg);
   double gvg = similarity(dg, v0);

   print.Debug("\ndx", dx, "\ndg", dg, "\ndelgam", delgam, "gvg", gvg);

   if (delgam == 0) {
      print.Warn("delgam = 0 : cannot update - return same matrix");
      return s0.Error();
   }

   if (delgam < 0) {
      print.Warn("delgam < 0 : first derivatives increasing along search line");
   }

   if (gvg <= 0) {
      // since v0 is pos def this gvg can be only = 0 if  dg = 0 - should never be here
      print.Warn("gvg <= 0 : cannot update - return same matrix");
      return s0.Error();
   }

   MnAlgebraicVector vg = v0 * dg;

   MnAlgebraicSymMatrix vUpd = Outer_product(dx) / delgam - Outer_product(vg) / gvg;

   if (delgam > gvg) {
      // use rank 1 formula
      vUpd += gvg * Outer_product(MnAlgebraicVector(dx / delgam - vg / gvg));
   }

   double sum_upd = sum_of_elements(vUpd);
   vUpd += v0;

   double dcov = 0.5 * (s0.Error().Dcovar() + sum_upd / sum_of_elements(vUpd));

   return MinimumError(vUpd, dcov);
}

/*
MinimumError DavidonErrorUpdator::Update(const MinimumState& s0,
                const MinimumParameters& p1,
                const FunctionGradient& g1) const {

  const MnAlgebraicSymMatrix& v0 = s0.Error().InvHessian();
  MnAlgebraicVector dx = p1.Vec() - s0.Vec();
  MnAlgebraicVector dg = g1.Vec() - s0.Gradient().Vec();

  double delgam = inner_product(dx, dg);
  double gvg = similarity(dg, v0);

//   std::cout<<"delgam= "<<delgam<<" gvg= "<<gvg<<std::endl;
  MnAlgebraicVector vg = v0*dg;
//   MnAlgebraicSymMatrix vUpd(v0.Nrow());

//   MnAlgebraicSymMatrix dd = ( 1./delgam )*outer_product(dx);
//   dd *= ( 1./delgam );
//   MnAlgebraicSymMatrix VggV = ( 1./gvg )*outer_product(vg);
//   VggV *= ( 1./gvg );
//   vUpd = dd - VggV;
//   MnAlgebraicSymMatrix vUpd = ( 1./delgam )*outer_product(dx) - ( 1./gvg )*outer_product(vg);
  MnAlgebraicSymMatrix vUpd = Outer_product(dx)/delgam - Outer_product(vg)/gvg;

  if(delgam > gvg) {
//     dx *= ( 1./delgam );
//     vg *= ( 1./gvg );
//     MnAlgebraicVector flnu = dx - vg;
//     MnAlgebraicSymMatrix tmp = Outer_product(flnu);
//     tmp *= gvg;
//     vUpd = vUpd + tmp;
    vUpd += gvg*outer_product(dx/delgam - vg/gvg);
  }

//
//     MnAlgebraicSymMatrix dd = Outer_product(dx);
//     dd *= ( 1./delgam );
//     MnAlgebraicSymMatrix VggV = Outer_product(vg);
//     VggV *= ( 1./gvg );
//     vUpd = dd - VggV;
//
//
//   double phi = delgam/(delgam - gvg);

//   MnAlgebraicSymMatrix vUpd(v0.Nrow());
//   if(phi < 0) {
//     // rank-2 Update
//     MnAlgebraicSymMatrix dd = Outer_product(dx);
//     dd *= ( 1./delgam );
//     MnAlgebraicSymMatrix VggV = Outer_product(vg);
//     VggV *= ( 1./gvg );
//     vUpd = dd - VggV;
//   }
//   if(phi > 1) {
//     // rank-1 Update
//     MnAlgebraicVector tmp = dx - vg;
//     vUpd = Outer_product(tmp);
//     vUpd *= ( 1./(delgam - gvg) );
//   }
//

//
//   if(delgam > gvg) {
//     // rank-1 Update
//     MnAlgebraicVector tmp = dx - vg;
//     vUpd = Outer_product(tmp);
//     vUpd *= ( 1./(delgam - gvg) );
//   } else {
//     // rank-2 Update
//     MnAlgebraicSymMatrix dd = Outer_product(dx);
//     dd *= ( 1./delgam );
//     MnAlgebraicSymMatrix VggV = Outer_productn(vg);
//     VggV *= ( 1./gvg );
//     vUpd = dd - VggV;
//   }
//

  double sum_upd = sum_of_elements(vUpd);
  vUpd += v0;

//   MnAlgebraicSymMatrix V1 = v0 + vUpd;

  double dcov =
    0.5*(s0.Error().Dcovar() + sum_upd/sum_of_elements(vUpd));

  return MinimumError(vUpd, dcov);
}
*/

} // namespace Minuit2

} // namespace ROOT
