// @(#)root/minuit2:$Id$
// Author: L. Moneta    10/2006

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006 ROOT Foundation,  CERN/PH-SFT                   *
 *                                                                    *
 **********************************************************************/

#ifndef ROOT_Minuit2_FCNGradAdapter
#define ROOT_Minuit2_FCNGradAdapter

#include "Minuit2/FCNGradientBase.h"
#include "Minuit2/MnPrint.h"

#include <vector>

namespace ROOT {

namespace Minuit2 {

/**


template wrapped class for adapting to FCNBase signature a IGradFunction

@author Lorenzo Moneta

@ingroup Minuit

*/

template <class Function>
class FCNGradAdapter : public FCNGradientBase {

public:
   FCNGradAdapter(const Function &f, double up = 1.) : fFunc(f), fUp(up), fGrad(std::vector<double>(fFunc.NDim())) {}

   ~FCNGradAdapter() {}

   double operator()(const std::vector<double> &v) const override { return fFunc.operator()(&v[0]); }
   double operator()(const double *v) const { return fFunc.operator()(v); }

   double Up() const override { return fUp; }

   std::vector<double> Gradient(const std::vector<double> &v) const override
   {
      fFunc.Gradient(&v[0], &fGrad[0]);

      MnPrint("FCNGradAdapter").Debug([&](std::ostream &os) {
         os << "gradient in FCNAdapter = {";
         for (unsigned int i = 0; i < fGrad.size(); ++i)
            os << fGrad[i] << (i == fGrad.size() - 1 ? '}' : '\t');
      });
      return fGrad;
   }
   std::vector<double> Gradient(const std::vector<double> &v, double *previous_grad, double *previous_g2,
                                double *previous_gstep) const override
   {
      fFunc.Gradient(&v[0], &fGrad[0], previous_grad, previous_g2, previous_gstep);

      MnPrint("FCNGradAdapter").Debug([&](std::ostream &os) {
         os << "gradient in FCNAdapter = {";
         for (unsigned int i = 0; i < fGrad.size(); ++i)
            os << fGrad[i] << (i == fGrad.size() - 1 ? '}' : '\t');
      });
      return fGrad;
   }
   // forward interface
   // virtual double operator()(int npar, double* params,int iflag = 4) const;
   bool CheckGradient() const override { return false; }

   GradientParameterSpace gradParameterSpace() const override {
      if (fFunc.returnsInMinuit2ParameterSpace()) {
         return GradientParameterSpace::Internal;
      } else {
         return GradientParameterSpace::External;
      }
   }

private:
   const Function &fFunc;
   double fUp;
   mutable std::vector<double> fGrad;
};

} // end namespace Minuit2

} // end namespace ROOT

#endif // ROOT_Minuit2_FCNGradAdapter
