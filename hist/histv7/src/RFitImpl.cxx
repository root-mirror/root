/// \file ROOT/RFitImpl.cxx
/// \ingroup Hist ROOT7
/// \author Claire Guyot
/// \date 2020-07
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "ROOT/RFitImpl.hxx"

namespace ROOT {
namespace Experimental {

bool RFit::AdjustError(const ROOT::Fit::DataOptions & option, double & error, double value) {
   if (error <= 0) {
      if (option.fUseEmpty || (option.fErrors1 && std::abs(value) > 0 ) )
         error = 1.;
      else
         return false;
   } else if (option.fErrors1)
      error = 1;
   return true;
}

int RFit::CheckFitFunction(const TF1 * f1, int dim)
{
   // Check validity of fitted function
   if (!f1) {
      Error("Fit", "function may not be null pointer");
      return -1;
   }
   if (f1->IsZombie()) {
      Error("Fit", "function is zombie");
      return -2;
   }

   int npar = f1->GetNpar();
   if (npar <= 0) {
      Error("Fit", "function %s has illegal number of parameters = %d", f1->GetName(), npar);
      return -3;
   }

   // Check that function has same dimension as histogram
   if (f1->GetNdim() > dim) {
      Error("Fit","function %s dimension, %d, is greater than fit object dimension, %d",
            f1->GetName(), f1->GetNdim(), dim);
      return -4;
   }
   if (f1->GetNdim() < dim-1) {
      Error("Fit","function %s dimension, %d, is smaller than fit object dimension -1, %d",
            f1->GetName(), f1->GetNdim(), dim);
      return -5;
   }

   return 0;

}

void RFit::GetFunctionRange(const TF1 & f1, ROOT::Fit::DataRange & range)
{
   // get the range form the function and fill and return the DataRange object
   Double_t fxmin, fymin, fzmin, fxmax, fymax, fzmax;
   f1.GetRange(fxmin, fymin, fzmin, fxmax, fymax, fzmax);
   // support only one range - so add only if was not set before
   if (range.Size(0) == 0) range.AddRange(0,fxmin,fxmax);
   if (range.Size(1) == 0) range.AddRange(1,fymin,fymax);
   if (range.Size(2) == 0) range.AddRange(2,fzmin,fzmax);
   return;
}

}// namespace Experimental
}// namespace ROOT
