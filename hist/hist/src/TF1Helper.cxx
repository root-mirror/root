// @(#)root/hist:$Id$
// Author: Lorenzo Moneta 12/06/07

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006  LCG ROOT Math Team, CERN/PH-SFT                *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// helper functions and classes used internally by TF1

#include "TF1Helper.h"
#include "TError.h"
#include "TMath.h"
#include <vector>
#include <cmath>
#include <cassert>

#include "TBackCompFitter.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TROOT.h"

#include "Math/IntegratorOptions.h"

int gDefMaxPts = ROOT::Math::IntegratorMultiDimOptions::DefaultNCalls();

namespace ROOT {



   namespace TF1Helper{





      double IntegralError(TF1 * func, Int_t ndim, const double * a, const double * b, const double * params, const double * covmat, double epsilon) {

   // calculate the eror on an integral from a to b of a parametetric function f when the parameters
   // are estimated from a fit and have an error represented by the covariance matrix of the fit.
   // The latest fit result is used

   // need to create the gradient functions w.r.t to the parameters


   // loop on all parameters
   bool onedim = ndim == 1;
   int npar = func->GetNpar();
   if (npar == 0) {
      Error("TF1Helper","Function has no parameters");
      return 0;
   }

   std::vector<double> oldParams;
   if (params) {
      // when using an external set of parameters
      oldParams.resize(npar);
      std::copy(func->GetParameters(), func->GetParameters()+npar, oldParams.begin());
      func->SetParameters(params);
   }


   TMatrixDSym covMatrix(npar);
   if (covmat == 0) {
      // with ROOT implicit MT there is no global TVirtualFitter
      if (ROOT::IsImplicitMTEnabled()) {
         Error("TF1Helper::IntegralError", "ROOT has enabled implicit MT. There is no existing lobal fitter, as shown in the documentation a pointer to the covariance matrix"
         "from the TFitResult must be passed to TF1::IntegralError");
         return 0;
      }
      // use matrix from last fit (needs to be a TBackCompFitter)
      TVirtualFitter * vfitter = TVirtualFitter::GetFitter();
      TBackCompFitter * fitter = dynamic_cast<TBackCompFitter*> (vfitter);
      if (fitter == 0) {
         Error("TF1Helper::IntegralError","No existing fitter can be used for computing the integral error");
         return 0;
      }
      // check that fitter and function are in sync
      if (fitter->GetNumberTotalParameters() != npar) {
         Error("TF1Helper::IntegralError","Last used fitter is not compatible with the current TF1");
         return 0;
      }
      // check that errors are provided
      if (int(fitter->GetFitResult().Errors().size()) != npar) {
         Warning("TF1Helper::INtegralError","Last used fitter does no provide parameter errors and a covariance matrix");
         return 0;
      }

      // check also the parameter values
      for (int i = 0; i < npar; ++i) {
         if (fitter->GetParameter(i) != func->GetParameter(i) ) {
            Error("TF1Helper::IntegralError","Last used Fitter has different parameter values");
            return 0;
         }
      }

      // fill the covariance matrix
      fitter->GetFitResult().GetCovarianceMatrix(covMatrix);
   }
   else {
      covMatrix.Use(npar,covmat);
   }



   // loop on the parameter and calculate the errors
   TVectorD ig(npar);

   if (epsilon <= 0) epsilon = 1.e7*ROOT::Math::IntegratorMultiDimOptions::DefaultRelTolerance();

   double numError2 = 0;
   for (int i=0; i < npar; ++i) {

      // use tolerance factor of 10 smaller than parameter errors
      double epsrel = TMath::Min(0.1, epsilon / std::abs(func->GetParameter(i))  );
      double epsabs = epsrel;

      // check that parameter error is not zero - otherwise skip it
      // should check the limits
      double integral  = 0;
      double error = 0;
      if (covMatrix(i,i) > 0 ) {
         TF1 gradFunc("gradFunc",TGradientParFunction(i,func),0,0,0);
         if (onedim) {
            integral = gradFunc.IntegralOneDim(*a,*b,epsrel,epsabs, error);
         }
         else {
            double relerr;
            int ifail = 0;
            int nfnevl = 0;
            int maxpts = ROOT::Math::IntegratorMultiDimOptions::DefaultNCalls();
            if (maxpts == gDefMaxPts) maxpts = 10*gDefMaxPts;   // increate points
            integral = gradFunc.IntegralMultiple(ndim,a,b,maxpts,epsrel,epsabs,relerr,nfnevl,ifail);
            //if (ifail) Warning("TF1Helper::IntegralError","n-dim integration failed code=%d I = %g, relerr =%g, ncall = %d, maxpts = %d, epsrel = %g, epsabs = %g, ",ifail,integral,relerr,nfnevl,maxpts,epsrel,epsabs);
            error = relerr*std::abs(integral);
         }
      }
      ig[i] = integral;
      // std::cout << " ipar  " << i << "  sigma " << sqrt(covMatrix(i,i)) << " rel " << sqrt(covMatrix(i,i))/std::abs(func->GetParameter(i)) << " integral " << integral << "  +/- " <<  error << "  " <<
      // error/std::abs(integral) << std::endl;

      // estimate numerical error (neglect correlations)
      numError2 += covMatrix(i,i)*covMatrix(i,i) * integral * integral * error * error;
   }
   double err2 =  covMatrix.Similarity(ig);
   double result = sqrt(err2);
   double numError = sqrt( numError2/sqrt(err2) );

   //std::cout << "integral error is " << result << " num error is " << numError << std::endl;
   if ( result > epsilon && numError > 2*epsilon*result )
      Warning("TF1Helper::IntegralError","numerical error from integration is too large. Integral error = %g +/- %g  - eps = %g",result,numError,epsilon);

   // restore old parameters in TF1
   if (!oldParams.empty()) {
      func->SetParameters(&oldParams.front());
   }



   return std::sqrt(err2);

}


} // end namespace TF1Helper


} // end namespace ROOT
