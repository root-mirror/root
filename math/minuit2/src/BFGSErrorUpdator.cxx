// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei   2003-2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 *                                                                    *
 **********************************************************************/

#include "Minuit2/BFGSErrorUpdator.h"
#include "Minuit2/MinimumState.h"
#include "Minuit2/LaSum.h"
#include "Minuit2/LaProd.h"
#include "Minuit2/MnPrint.h"

#include <vector>

namespace ROOT {

namespace Minuit2 {

double inner_product(const LAVector &, const LAVector &);
double similarity(const LAVector &, const LASymMatrix &);
double sum_of_elements(const LASymMatrix &);

// define here a square matrix that it is needed for computingthe BFGS update
//  define just the class, no need for defining operatipons as dane for the Symmetric matrices
// since the square matrix will be converted in a symmetric one afterwards

class LASquareMatrix {
public:
   LASquareMatrix(unsigned int n) : fNRow(n), fData(std::vector<double>(n * n)) {}

   double operator()(unsigned int row, unsigned int col) const
   {
      assert(row < fNRow && col < fNRow);
      return fData[col + row * fNRow];
   }

   double &operator()(unsigned int row, unsigned int col)
   {
      assert(row < fNRow && col < fNRow);
      return fData[col + row * fNRow];
   }

   unsigned int Nrow() const { return fNRow; }

private:
   unsigned int fNRow;
   std::vector<double> fData;
};

// compute outer product of two vector of same size to return a square matrix
LASquareMatrix OuterProduct(const LAVector &v1, const LAVector &v2)
{
   assert(v1.size() == v2.size());
   LASquareMatrix a(v1.size());
   for (unsigned int i = 0; i < v1.size(); ++i) {
      for (unsigned int j = 0; j < v2.size(); ++j) {
         a(i, j) = v1[i] * v2[j];
      }
   }
   return a;
}
// compute product of symmetric matrix with square matrix

LASquareMatrix MatrixProduct(const LASymMatrix &m1, const LASquareMatrix &m2)
{
   unsigned int n = m1.Nrow();
   assert(n == m2.Nrow());
   LASquareMatrix a(n);
   for (unsigned int i = 0; i < n; ++i) {
      for (unsigned int j = 0; j < n; ++j) {
         a(i, j) = 0;
         for (unsigned int k = 0; k < n; ++k) {
            a(i, j) += m1(i, k) * m2(k, j);
         }
      }
   }
   return a;
}

MinimumError
BFGSErrorUpdator::Update(const MinimumState &s0, const MinimumParameters &p1, const FunctionGradient &g1) const
{

   // update of the covarianze matrix (BFGS formula, see Tutorial, par. 4.8 pag 26)
   // in case of delgam > gvg (PHI > 1) use rank one formula
   // see  par 4.10 pag 30

   const MnAlgebraicSymMatrix &v0 = s0.Error().InvHessian();
   MnAlgebraicVector dx = p1.Vec() - s0.Vec();
   MnAlgebraicVector dg = g1.Vec() - s0.Gradient().Vec();

   double delgam = inner_product(dx, dg); // this is s^T y  using wikipedia conventions
   double gvg = similarity(dg, v0);       // this is y^T B^-1 y

   MnPrint print("BFGSErrorUpdator");
   print.Debug("dx", dx, "dg", dg, "delgam", delgam, "gvg", gvg);

   if (delgam == 0) {
      print.Warn("delgam = 0 : cannot update - return same matrix");
      return s0.Error();
   }

   if (delgam < 0) {
      print.Warn("delgam < 0 : first derivatives increasing along search line");
   }

   if (gvg <= 0) {
      // since v0 is pos def this gvg can be only = 0 if  dg = 0 - should never be here
      print.Warn("gvg <= 0");
      // return s0.Error();
   }

   // compute update formula for BFGS
   // see wikipedia  https://en.wikipedia.org/wiki/Broyden–Fletcher–Goldfarb–Shanno_algorithm
   // need to compute outer product dg . dx and it is not symmetric anymore

   LASquareMatrix a = OuterProduct(dg, dx);
   LASquareMatrix b = MatrixProduct(v0, a);

   unsigned int n = v0.Nrow();
   MnAlgebraicSymMatrix v2(n);
   for (unsigned int i = 0; i < n; ++i) {
      for (unsigned int j = i; j < n; ++j) {
         v2(i, j) = (b(i, j) + b(j, i)) / (delgam);
      }
   }

   MnAlgebraicSymMatrix vUpd = (delgam + gvg) * Outer_product(dx) / (delgam * delgam);
   vUpd -= v2;

   double sum_upd = sum_of_elements(vUpd);
   vUpd += v0;

   double dcov = 0.5 * (s0.Error().Dcovar() + sum_upd / sum_of_elements(vUpd));

   print.Debug("dcov", dcov);

   return MinimumError(vUpd, dcov);
}

} // namespace Minuit2

} // namespace ROOT
