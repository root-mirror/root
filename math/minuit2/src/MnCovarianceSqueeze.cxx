// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "Minuit2/MnCovarianceSqueeze.h"
#include "Minuit2/MnUserCovariance.h"
#include "Minuit2/MinimumError.h"
#include "Minuit2/MnPrint.h"

namespace ROOT {

namespace Minuit2 {

MnUserCovariance MnCovarianceSqueeze::operator()(const MnUserCovariance &cov, unsigned int n) const
{
   // squeeze MnUserCovariance class
   // MnUserCovariance contasins the error matrix. Need to invert first to get the hessian, then
   // after having squuezed the hessian, need to invert again to get the new error matrix
   assert(cov.Nrow() > 0);
   assert(n < cov.Nrow());

   MnPrint print("MnCovarianceSqueeze");

   MnAlgebraicSymMatrix hess(cov.Nrow());
   for (unsigned int i = 0; i < cov.Nrow(); i++) {
      for (unsigned int j = i; j < cov.Nrow(); j++) {
         hess(i, j) = cov(i, j);
      }
   }

   int ifail = Invert(hess);

   if (ifail != 0) {
      print.Warn("inversion failed; return diagonal matrix;");
      MnUserCovariance result(cov.Nrow() - 1);
      for (unsigned int i = 0, j = 0; i < cov.Nrow(); i++) {
         if (i == n)
            continue;
         result(j, j) = cov(i, i);
         j++;
      }
      return result;
   }

   MnAlgebraicSymMatrix squeezed = (*this)(hess, n);

   ifail = Invert(squeezed);
   if (ifail != 0) {
      print.Warn("back-inversion failed; return diagonal matrix;");
      MnUserCovariance result(squeezed.Nrow());
      for (unsigned int i = 0; i < squeezed.Nrow(); i++) {
         result(i, i) = 1. / squeezed(i, i);
      }
      return result;
   }

   return MnUserCovariance(std::vector<double>(squeezed.Data(), squeezed.Data() + squeezed.size()), squeezed.Nrow());
}

MinimumError MnCovarianceSqueeze::operator()(const MinimumError &err, unsigned int n) const
{

   MnPrint print("MnCovarianceSqueeze");

   // squueze the minimum error class
   // Remove index-row on the Hessian matrix and the get the new correct error matrix
   // (inverse of new Hessian)
   MnAlgebraicSymMatrix hess = err.Hessian();
   MnAlgebraicSymMatrix squeezed = (*this)(hess, n);
   int ifail = Invert(squeezed);
   if (ifail != 0) {
      print.Warn("MinimumError inversion fails; return diagonal matrix.");

      MnAlgebraicSymMatrix tmp(squeezed.Nrow());
      for (unsigned int i = 0; i < squeezed.Nrow(); i++) {
         tmp(i, i) = 1. / squeezed(i, i);
      }
      return MinimumError(tmp, MinimumError::MnInvertFailed());
   }

   return MinimumError(squeezed, err.Dcovar());
}

MnAlgebraicSymMatrix MnCovarianceSqueeze::operator()(const MnAlgebraicSymMatrix &hess, unsigned int n) const
{
   // squueze a symmetrix matrix (remove entire row and column n)
   assert(hess.Nrow() > 0);
   assert(n < hess.Nrow());

   MnAlgebraicSymMatrix hs(hess.Nrow() - 1);
   for (unsigned int i = 0, j = 0; i < hess.Nrow(); i++) {
      if (i == n)
         continue;
      for (unsigned int k = i, l = j; k < hess.Nrow(); k++) {
         if (k == n)
            continue;
         hs(j, l) = hess(i, k);
         l++;
      }
      j++;
   }

   return hs;
}

} // namespace Minuit2

} // namespace ROOT
