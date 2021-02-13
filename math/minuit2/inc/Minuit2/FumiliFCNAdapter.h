// @(#)root/minuit2:$Id$
// Author: L. Moneta    10/2006

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006 ROOT Foundation,  CERN/PH-SFT                   *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_FumiliFCNAdapter
#define ROOT_Minuit2_FumiliFCNAdapter

#include "Minuit2/FumiliFCNBase.h"

#include "Math/FitMethodFunction.h"

#include "Minuit2/MnPrint.h"

// #ifndef ROOT_Math_Util
// #include "Math/Util.h"
// #endif

#include <cmath>
#include <cassert>
#include <vector>

namespace ROOT {

namespace Minuit2 {

/**


template wrapped class for adapting to FumiliFCNBase signature

@author Lorenzo Moneta

@ingroup Minuit

*/

template <class Function>
class FumiliFCNAdapter : public FumiliFCNBase {

public:
   //   typedef ROOT::Math::FitMethodFunction Function;
   typedef typename Function::Type_t Type_t;

   FumiliFCNAdapter(const Function &f, unsigned int ndim, double up = 1.) : FumiliFCNBase(ndim), fFunc(f), fUp(up) {}

   ~FumiliFCNAdapter() {}

   double operator()(const std::vector<double> &v) const { return fFunc.operator()(&v[0]); }
   double operator()(const double *v) const { return fFunc.operator()(v); }
   double Up() const { return fUp; }

   void SetErrorDef(double up) { fUp = up; }

   // virtual std::vector<double> Gradient(const std::vector<double>&) const;

   // forward interface
   // virtual double operator()(int npar, double* params,int iflag = 4) const;

   /**
       evaluate gradient hessian and function value needed by fumili
     */
   void EvaluateAll(const std::vector<double> &v);

private:
   // data member

   const Function &fFunc;
   double fUp;
};

template <class Function>
void FumiliFCNAdapter<Function>::EvaluateAll(const std::vector<double> &v)
{
   MnPrint print("FumiliFCNAdaptor");

   // typedef FumiliFCNAdapter::Function Function;

   // evaluate all elements
   unsigned int npar = Dimension();
   if (npar != v.size())
      print.Error("npar", npar, "v.size()", v.size());
   assert(npar == v.size());
   // must distinguish case of likelihood or LS

   std::vector<double> &grad = Gradient();
   std::vector<double> &hess = Hessian();
   // reset
   assert(grad.size() == npar);
   grad.assign(npar, 0.0);
   hess.assign(hess.size(), 0.0);

   double sum = 0;
   unsigned int ndata = fFunc.NPoints();

   std::vector<double> gf(npar);

   // loop on the data points

   // assume for now least-square
   if (fFunc.Type() == Function::kLeastSquare) {

      for (unsigned int i = 0; i < ndata; ++i) {
         // calculate data element and gradient
         double fval = fFunc.DataElement(&v.front(), i, &gf[0]);

         // t.b.d should protect for bad  values of fval
         sum += fval * fval;

         for (unsigned int j = 0; j < npar; ++j) {
            grad[j] += 2. * fval * gf[j];
            for (unsigned int k = j; k < npar; ++k) {
               int idx = j + k * (k + 1) / 2;
               hess[idx] += 2.0 * gf[j] * gf[k];
            }
         }
      }
   } else if (fFunc.Type() == Function::kLogLikelihood) {

      for (unsigned int i = 0; i < ndata; ++i) {

         // calculate data element and gradient
         // return value is log of pdf and derivative of the log(Pdf)
         double fval = fFunc.DataElement(&v.front(), i, &gf[0]);

         sum -= fval;

         for (unsigned int j = 0; j < npar; ++j) {
            double gfj = gf[j];
            grad[j] -= gfj;
            for (unsigned int k = j; k < npar; ++k) {
               int idx = j + k * (k + 1) / 2;
               hess[idx] += gfj * gf[k];
            }
         }
      }
   } else {
      print.Error("Type of fit method is not supported, it must be chi2 or log-likelihood");
   }
}

} // end namespace Minuit2

} // end namespace ROOT

#endif // ROOT_Minuit2_FCNAdapter
