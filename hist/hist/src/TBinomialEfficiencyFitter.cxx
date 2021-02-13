// @(#)root/hist:$Id$
// Author: Frank Filthaut, Rene Brun   30/05/2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TBinomialEfficiencyFitter
    \ingroup Hist
    \brief Binomial fitter for the division of two histograms.

Use when you need to calculate a selection's efficiency from two histograms,
one containing all entries, and one containing the subset of these entries
that pass the selection, and when you have a parametrization available for
the efficiency as a function of the variable(s) under consideration.

A very common problem when estimating efficiencies is that of error estimation:
when no other information is available than the total number of events N and
the selected number n, the best estimate for the selection efficiency p is n/N.
Standard binomial statistics dictates that the uncertainty (this presupposes
sufficiently high statistics that an approximation by a normal distribution is
reasonable) on p, given N, is
\f[
    \sqrt{\frac{p(1-p)}{N}}
\f]
However, when p is estimated as n/N, fluctuations from the true p to its
estimate become important, especially for low numbers of events, and giving
rise to biased results.

When fitting a parametrized efficiency, these problems can largely be overcome,
as a hypothesized true efficiency is available by construction. Even so, simply
using the corresponding uncertainty still presupposes that Gaussian errors
yields a reasonable approximation. When using, instead of binned efficiency
histograms, the original numerator and denominator histograms, a binned maximum
likelihood can be constructed as the product of bin-by-bin binomial probabilities
to select n out of N events. Assuming that a correct parametrization of the
efficiency is provided, this construction in general yields less biased results
(and is much less sensitive to binning details).

A generic use of this method is given below (note that the method works for 2D
and 3D histograms as well):

~~~ {.cpp}
    {
        TH1* denominator;              // denominator histogram
        TH1* numerator;                // corresponding numerator histogram
        TF1* eff;                      // efficiency parametrization
        ....                           // set step sizes and initial parameter
        ....                           //   values for the fit function
        ....                           // possibly also set ranges, see TF1::SetRange()
        TBinomialEfficiencyFitter* f = new TBinomialEfficiencyFitter(
                                     numerator, denominator);
        Int_t status = f->Fit(eff, "I");
        if (status == 0) {
          // if the fit was successful, display bin-by-bin efficiencies
          // as well as the result of the fit
          numerator->Sumw2();
          TH1* hEff = dynamic_cast<TH1*>(numerator->Clone("heff"));
          hEff->Divide(hEff, denominator, 1.0, 1.0, "B");
          hEff->Draw("E");
          eff->Draw("same");
        }
    }
~~~

Note that this method cannot be expected to yield reliable results when using
weighted histograms (because the likelihood computation will be incorrect).

*/

#include "TBinomialEfficiencyFitter.h"

#include "TMath.h"
#include "TPluginManager.h"
#include "TROOT.h"
#include "TH1.h"
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "Fit/FitConfig.h"
#include "Fit/Fitter.h"
#include "TFitResult.h"
#include "Math/Functor.h"
#include "Math/WrappedMultiTF1.h"

#include <limits>


const Double_t kDefaultEpsilon = 1E-12;

ClassImp(TBinomialEfficiencyFitter);


////////////////////////////////////////////////////////////////////////////////
/// default constructor

TBinomialEfficiencyFitter::TBinomialEfficiencyFitter() {
   fNumerator   = 0;
   fDenominator = 0;
   fFunction    = 0;
   fFitDone     = kFALSE;
   fAverage     = kFALSE;
   fRange       = kFALSE;
   fEpsilon     = kDefaultEpsilon;
   fFitter      = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.
///
/// Note that no objects are copied, so it is up to the user to ensure that the
/// histogram pointers remain valid.
///
/// Both histograms need to be "consistent". This is not checked here, but in
/// TBinomialEfficiencyFitter::Fit().

TBinomialEfficiencyFitter::TBinomialEfficiencyFitter(const TH1 *numerator, const TH1 *denominator) {
   fEpsilon  = kDefaultEpsilon;
   fFunction = 0;
   fFitter   = 0;
   Set(numerator,denominator);
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TBinomialEfficiencyFitter::~TBinomialEfficiencyFitter() {
   if (fFitter) delete fFitter;
   fFitter = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize with a new set of inputs.

void TBinomialEfficiencyFitter::Set(const TH1 *numerator, const TH1 *denominator)
{
   fNumerator   = (TH1*)numerator;
   fDenominator = (TH1*)denominator;

   fFitDone     = kFALSE;
   fAverage     = kFALSE;
   fRange       = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set the required integration precision, see TF1::Integral()

void TBinomialEfficiencyFitter::SetPrecision(Double_t epsilon)
{
   fEpsilon = epsilon;
}

////////////////////////////////////////////////////////////////////////////////
/// Provide access to the underlying fitter object.
/// This may be useful e.g. for the retrieval of additional information (such
/// as the output covariance matrix of the fit).

ROOT::Fit::Fitter* TBinomialEfficiencyFitter::GetFitter()
{
   if (!fFitter)  fFitter = new ROOT::Fit::Fitter();
   return fFitter;

}

////////////////////////////////////////////////////////////////////////////////
/// Carry out the fit of the given function to the given histograms.
///
/// If option "I" is used, the fit function will be averaged over the
/// bin (the default is to evaluate it simply at the bin center).
///
/// If option "R" is used, the fit range will be taken from the fit
/// function (the default is to use the entire histogram).
///
/// If option "S" a TFitResult object is returned and it can be used to obtain
///  additional fit information, like covariance or correlation matrix.
///
/// Note that all parameter values, limits, and step sizes are copied
/// from the input fit function f1 (so they should be set before calling
/// this method. This is particularly relevant for the step sizes, taken
/// to be the "error" set on input, as a null step size usually fixes the
/// corresponding parameter. That is protected against, but in such cases
/// an arbitrary starting step size will be used, and the reliability of
/// the fit should be questioned). If parameters are to be fixed, this
/// should be done by specifying non-null parameter limits, with lower
/// limits larger than upper limits.
///
/// On output, f1 contains the fitted parameters and errors, as well as
/// the number of degrees of freedom, and the goodness-of-fit estimator
/// as given by S. Baker and R. Cousins, Nucl. Instr. Meth. A221 (1984) 437.

TFitResultPtr TBinomialEfficiencyFitter::Fit(TF1 *f1, Option_t* option)
{
   TString opt = option;
   opt.ToUpper();
   fAverage  = opt.Contains("I");
   fRange    = opt.Contains("R");
   Bool_t verbose    = opt.Contains("V");
   Bool_t quiet      = opt.Contains("Q");
   Bool_t saveResult  = opt.Contains("S");
   if (!f1) return -1;
   fFunction = (TF1*)f1;
   Int_t i, npar;
   npar = f1->GetNpar();
   if (npar <= 0) {
      Error("Fit", "function %s has illegal number of parameters = %d",
            f1->GetName(), npar);
      return -3;
   }

   // Check that function has same dimension as histogram
   if (!fNumerator || !fDenominator) {
      Error("Fit","No numerator or denominator histograms set");
      return -5;
   }
   if (f1->GetNdim() != fNumerator->GetDimension()) {
      Error("Fit","function %s dimension, %d, does not match histogram dimension, %d",
            f1->GetName(), f1->GetNdim(), fNumerator->GetDimension());
      return -4;
   }
   // Check that the numbers of bins for the histograms match
   if (fNumerator->GetNbinsX() != fDenominator->GetNbinsX() ||
       (f1->GetNdim() > 1 && fNumerator->GetNbinsY() != fDenominator->GetNbinsY()) ||
       (f1->GetNdim() > 2 && fNumerator->GetNbinsZ() != fDenominator->GetNbinsZ())) {
      Error("Fit", "numerator and denominator histograms do not have identical numbers of bins");
      return -6;
   }

   // initialize the fitter

   if (!fFitter) {
      fFitter = new ROOT::Fit::Fitter();
   }


   std::vector<ROOT::Fit::ParameterSettings> & parameters = fFitter->Config().ParamsSettings();
   parameters.reserve(npar);
   for (i = 0; i < npar; i++) {

      // assign an ARBITRARY starting error to ensure the parameter won't be fixed!
      Double_t we = f1->GetParError(i);
      if (we <= 0) we = 0.3*TMath::Abs(f1->GetParameter(i));
      if (we == 0) we = 0.01;

      parameters.push_back(ROOT::Fit::ParameterSettings(f1->GetParName(i), f1->GetParameter(i), we) );

      Double_t plow, pup;
      f1->GetParLimits(i,plow,pup);
      // when a parameter is fixed is having plow and pup equal to the value (if this is not zero)
      // we handle special case when fixed parameter has zero value (in that case plow=1 and pup =1 )
      if (plow >= pup && (plow==f1->GetParameter(i) || pup==f1->GetParameter(i) ||
      ( f1->GetParameter(i) == 0 && plow==1. && pup == 1.) ) ) {
         parameters.back().Fix();
         Info("Fit", "Fixing parameter %s to value %f", f1->GetParName(i), f1->GetParameter(i));
      } else if (plow < pup) {
         parameters.back().SetLimits(plow,pup);
         Info("Fit", "Setting limits for parameter %s to [%f,%f]", f1->GetParName(i), plow,pup);
      }
   }

   // fcn must be set after setting the parameters
   ROOT::Math::Functor fcnFunction(this, &TBinomialEfficiencyFitter::EvaluateFCN, npar);

   // set also model function in fitter to have it in FitResult
   // in this way one can compute for example the confidence intervals
   fFitter->SetFCN(static_cast<ROOT::Math::IMultiGenFunction &>(fcnFunction), ROOT::Math::WrappedMultiTF1(*f1));

   // in case default value of 1.0 is used
   if (fFitter->Config().MinimizerOptions().ErrorDef() == 1.0 ) {
      fFitter->Config().MinimizerOptions().SetErrorDef(0.5);
   }

   if (verbose)   {
      fFitter->Config().MinimizerOptions().SetPrintLevel(3);
   }
   else if (quiet) {
      fFitter->Config().MinimizerOptions().SetPrintLevel(0);
   }

   // perform the actual fit

   fFitDone = kTRUE;
   Bool_t status = fFitter->FitFCN();
   if ( !status && !quiet)
      Warning("Fit","Abnormal termination of minimization.");


   //Store fit results in fitFunction
   const ROOT::Fit::FitResult & fitResult = fFitter->Result();
   if (!fitResult.IsEmpty() ) {
      // set in f1 the result of the fit
      f1->SetNDF(fitResult.Ndf() );

      //f1->SetNumberFitPoints(...);  // this is set in ComputeFCN

      f1->SetParameters( &(fitResult.Parameters().front()) );
      if ( int( fitResult.Errors().size()) >= f1->GetNpar() )
         f1->SetParErrors( &(fitResult.Errors().front()) );

      f1->SetChisquare(2.*fitResult.MinFcnValue());    // store goodness of fit (Baker&Cousins)
      f1->SetNDF(f1->GetNumberFitPoints()- fitResult.NFreeParameters());
      Info("result"," chi2 %f ndf %d ",2.*fitResult.MinFcnValue(), fitResult.Ndf() );

   }
   // create a new result class if needed
   if (saveResult) {
      TFitResult* fr = new TFitResult(fitResult);
      TString name = TString::Format("TBinomialEfficiencyFitter_result_of_%s",f1->GetName() );
      fr->SetName(name); fr->SetTitle(name);
      return TFitResultPtr(fr);
   }
   else {
      return TFitResultPtr(fitResult.Status() );
   }

}

////////////////////////////////////////////////////////////////////////////////
/// Compute the likelihood.

void TBinomialEfficiencyFitter::ComputeFCN(Double_t& f, const Double_t* par)
{
   int nDim = fDenominator->GetDimension();

   int xlowbin  = fDenominator->GetXaxis()->GetFirst();
   int xhighbin = fDenominator->GetXaxis()->GetLast();
   int ylowbin = 0, yhighbin = 0, zlowbin = 0, zhighbin = 0;
   if (nDim > 1) {
      ylowbin  = fDenominator->GetYaxis()->GetFirst();
      yhighbin = fDenominator->GetYaxis()->GetLast();
      if (nDim > 2) {
         zlowbin  = fDenominator->GetZaxis()->GetFirst();
         zhighbin = fDenominator->GetZaxis()->GetLast();
      }
   }

   fFunction->SetParameters(par);

   if (fRange) {
      double xmin, xmax, ymin, ymax, zmin, zmax;

      // This way to ensure that a minimum range chosen exactly at a
      // bin boundary is far from elegant, but is hopefully adequate.

      if (nDim == 1) {
         fFunction->GetRange(xmin, xmax);
         xlowbin  = fDenominator->GetXaxis()->FindBin(xmin);
         xhighbin = fDenominator->GetXaxis()->FindBin(xmax);
      } else if (nDim == 2) {
         fFunction->GetRange(xmin, ymin, xmax, ymax);
         xlowbin  = fDenominator->GetXaxis()->FindBin(xmin);
         xhighbin = fDenominator->GetXaxis()->FindBin(xmax);
         ylowbin  = fDenominator->GetYaxis()->FindBin(ymin);
         yhighbin = fDenominator->GetYaxis()->FindBin(ymax);
      } else if (nDim == 3) {
         fFunction->GetRange(xmin, ymin, zmin, xmax, ymax, zmax);
         xlowbin  = fDenominator->GetXaxis()->FindBin(xmin);
         xhighbin = fDenominator->GetXaxis()->FindBin(xmax);
         ylowbin  = fDenominator->GetYaxis()->FindBin(ymin);
         yhighbin = fDenominator->GetYaxis()->FindBin(ymax);
         zlowbin  = fDenominator->GetZaxis()->FindBin(zmin);
         zhighbin = fDenominator->GetZaxis()->FindBin(zmax);
      }
   }

   // The coding below is perhaps somewhat awkward -- but it is done
   // so that 1D, 2D, and 3D cases can be covered in the same loops.

   f = 0.;

   Int_t npoints = 0;
   Double_t nmax = 0;
   for (int xbin = xlowbin; xbin <= xhighbin; ++xbin) {

      // compute the bin edges
      Double_t xlow = fDenominator->GetXaxis()->GetBinLowEdge(xbin);
      Double_t xup  = fDenominator->GetXaxis()->GetBinLowEdge(xbin+1);

      for (int ybin = ylowbin; ybin <= yhighbin; ++ybin) {

         // compute the bin edges (if applicable)
         Double_t ylow  = (nDim > 1) ? fDenominator->GetYaxis()->GetBinLowEdge(ybin) : 0;
         Double_t yup   = (nDim > 1) ? fDenominator->GetYaxis()->GetBinLowEdge(ybin+1) : 0;

         for (int zbin = zlowbin; zbin <= zhighbin; ++zbin) {

            // compute the bin edges (if applicable)
            Double_t zlow  = (nDim > 2) ? fDenominator->GetZaxis()->GetBinLowEdge(zbin) : 0;
            Double_t zup   = (nDim > 2) ? fDenominator->GetZaxis()->GetBinLowEdge(zbin+1) : 0;

            int bin = fDenominator->GetBin(xbin, ybin, zbin);
            Double_t nDen = fDenominator->GetBinContent(bin);
            Double_t nNum = fNumerator->GetBinContent(bin);

            // count maximum value to use in the likelihood for inf
            // i.e. a number much larger than the other terms
            if (nDen> nmax) nmax = nDen;
            if (nDen <= 0.) continue;
            npoints++;

            // mu is the average of the function over the bin OR
            // the function evaluated at the bin centre
            // As yet, there is nothing to prevent mu from being
            // outside the range <0,1> !!

            Double_t mu = 0;
            switch (nDim) {
               case 1:
                  mu = (fAverage) ?
                     fFunction->Integral(xlow, xup, fEpsilon)
                        / (xup-xlow) :
                     fFunction->Eval(fDenominator->GetBinCenter(bin));
                  break;
               case 2:
                  {
                     mu = (fAverage) ?
                        ((TF2*)fFunction)->Integral(xlow, xup, ylow, yup, fEpsilon)
                        / ((xup-xlow)*(yup-ylow)) :
                     fFunction->Eval(fDenominator->GetXaxis()->GetBinCenter(xbin),
                     fDenominator->GetYaxis()->GetBinCenter(ybin));
                  }
                  break;
               case 3:
                  {
                     mu = (fAverage) ?
                        ((TF3*)fFunction)->Integral(xlow, xup, ylow, yup, zlow, zup, fEpsilon)
                           / ((xup-xlow)*(yup-ylow)*(zup-zlow)) :
                        fFunction->Eval(fDenominator->GetXaxis()->GetBinCenter(xbin),
                                 fDenominator->GetYaxis()->GetBinCenter(ybin),
                                 fDenominator->GetZaxis()->GetBinCenter(zbin));
                  }
            }

            // binomial formula (forgetting about the factorials)
            if (nNum != 0.) {
               if (mu > 0.)
                  f -= nNum * TMath::Log(mu*nDen/nNum);
               else
                  f -= nmax * -1E30; // crossing our fingers
            }
            if (nDen - nNum != 0.) {
               if (1. - mu > 0.)
                  f -= (nDen - nNum) * TMath::Log((1. - mu)*nDen/(nDen-nNum));
               else
                  f -= nmax * -1E30; // crossing our fingers
            }
         }
      }
   }

   fFunction->SetNumberFitPoints(npoints);
}
