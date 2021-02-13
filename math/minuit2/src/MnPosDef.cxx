// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "Minuit2/MnPosDef.h"
#include "Minuit2/MinimumState.h"
#include "Minuit2/MnMachinePrecision.h"
#include "Minuit2/MnPrint.h"

#include <algorithm>

namespace ROOT {

namespace Minuit2 {

LAVector eigenvalues(const LASymMatrix &);

MinimumState MnPosDef::operator()(const MinimumState &st, const MnMachinePrecision &prec) const
{
   // interface from minimum state
   MinimumError err = (*this)(st.Error(), prec);
   return MinimumState(st.Parameters(), err, st.Gradient(), st.Edm(), st.NFcn());
}

MinimumError MnPosDef::operator()(const MinimumError &e, const MnMachinePrecision &prec) const
{
   MnPrint print("MnPosDef");

   // make error matrix positive defined returning a new corrected minimum error state

   MnAlgebraicSymMatrix err(e.InvHessian());
   if (err.size() == 1 && err(0, 0) < prec.Eps()) {
      err(0, 0) = 1.;
      return MinimumError(err, MinimumError::MnMadePosDef());
   }
   if (err.size() == 1 && err(0, 0) > prec.Eps()) {
      return e;
   }
   //   std::cout<<"MnPosDef init matrix= "<<err<<std::endl;

   double epspdf = std::max(1.e-6, prec.Eps2());
   double dgmin = err(0, 0);

   for (unsigned int i = 0; i < err.Nrow(); i++) {
      if (err(i, i) <= 0 /* prec.Eps2() */)
         print.Warn("non-positive diagonal element in covariance matrix[", i, "] =", err(i, i));

      if (err(i, i) < dgmin)
         dgmin = err(i, i);
   }
   double dg = 0.;
   if (dgmin <= 0) {
      // dg = 1. + epspdf - dgmin;
      dg = 0.5 + epspdf - dgmin;
      //     dg = 0.5*(1. + epspdf - dgmin);

      print.Warn("Added to diagonal of Error matrix a value", dg);

      // std::cout << "Error matrix " << err << std::endl;
   }

   MnAlgebraicVector s(err.Nrow());
   MnAlgebraicSymMatrix p(err.Nrow());
   for (unsigned int i = 0; i < err.Nrow(); i++) {
      err(i, i) += dg;
      if (err(i, i) < 0.)
         err(i, i) = 1.;
      s(i) = 1. / std::sqrt(err(i, i));
      for (unsigned int j = 0; j <= i; j++) {
         p(i, j) = err(i, j) * s(i) * s(j);
      }
   }

   // std::cout<<"MnPosDef p: "<<p<<std::endl;
   MnAlgebraicVector eval = eigenvalues(p);
   double pmin = eval(0);
   double pmax = eval(eval.size() - 1);
   // std::cout<<"pmin= "<<pmin<<" pmax= "<<pmax<<std::endl;
   pmax = std::max(std::fabs(pmax), 1.);
   if (pmin > epspdf * pmax)
      return MinimumError(err, e.Dcovar());

   double padd = 0.001 * pmax - pmin;

   for (unsigned int i = 0; i < err.Nrow(); i++)
      err(i, i) *= (1. + padd);

   print.Debug([&](std::ostream &os) {
      os << "Eigenvalues:";
      for (unsigned i = 0; i < err.Nrow(); ++i)
         os << "\n  " << eval(i);
   });

   //   std::cout<<"MnPosDef final matrix: "<<err<<std::endl;

   print.Warn("Matrix forced pos-def by adding to diagonal", padd);

   return MinimumError(err, MinimumError::MnMadePosDef());
}

} // namespace Minuit2

} // namespace ROOT
