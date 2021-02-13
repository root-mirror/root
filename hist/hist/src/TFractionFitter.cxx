// @(#)root/hist:$Id$
// Author: Frank Filthaut F.Filthaut@science.ru.nl  20/05/2002
// with additions by Bram Wijngaarden <dwijngaa@hef.kun.nl>

/** \class TFractionFitter
Fits MC fractions to data histogram. A la HMCMLL, see R. Barlow and C. Beeston,
Comp. Phys. Comm. 77 (1993) 219-228, and http://www.hep.man.ac.uk/~roger/hfrac.f

The virtue of this fit is that it takes into account both data and Monte Carlo
statistical uncertainties. The way in which this is done is through a standard
likelihood fit using Poisson statistics; however, the template (MC) predictions
are also varied within statistics, leading to additional contributions to the
overall likelihood. This leads to many more fit parameters (one per bin per
template), but the minimisation with respect to these additional parameters is
done analytically rather than introducing them as formal fit parameters. Some
special care needs to be taken in the case of bins with zero content. For more
details please see the original publication cited above.

An example application of this fit is given below. For a TH1* histogram
("data") fitted as the sum of three Monte Carlo sources ("mc"):

~~~{.cpp}
{
  TH1F *data;                              //data histogram
  TH1F *mc0;                               // first MC histogram
  TH1F *mc1;                               // second MC histogram
  TH1F *mc2;                               // third MC histogram
  ....                                     // retrieve histograms
   TObjArray *mc = new TObjArray(3);        // MC histograms are put in this array
   mc->Add(mc0);
   mc->Add(mc1);
   mc->Add(mc2);
   TFractionFitter* fit = new TFractionFitter(data, mc); // initialise
   fit->Constrain(1,0.0,1.0);               // constrain fraction 1 to be between 0 and 1
   fit->SetRangeX(1,15);                    // use only the first 15 bins in the fit
   Int_t status = fit->Fit();               // perform the fit
   std::cout << "fit status: " << status << std::endl;
   if (status == 0) {                       // check on fit status
     TH1F* result = (TH1F*) fit->GetPlot();
     data->Draw("Ep");
     result->Draw("same");
   }
}
~~~

## Assumptions
A few assumptions need to be made for the fit procedure to be carried out:
 1 The total number of events in each template is not too small
    (so that its Poisson uncertainty can be neglected).
 2 The number of events in each bin is much smaller than the total
    number of events in each template (so that multinomial
    uncertainties can be replaced with Poisson uncertainties).

Biased fit uncertainties may result if these conditions are not fulfilled
(see e.g. arXiv:0803.2711).

## Instantiation
A fit object is instantiated through
    TFractionFitter* fit = new TFractionFitter(data, mc);
A number of basic checks (intended to ensure that the template histograms
represent the same "kind" of distribution as the data one) are carried out.
The TVirtualFitter object is then addressed and all fit parameters (the
template fractions) declared (initially unbounded).

## Applying constraints
Fit parameters can be constrained through

 fit->Constrain(parameter #, lower bound, upper bound);

Setting lower bound = upper bound = 0 removes the constraint (a la Minuit);
however, a function

    fit->Unconstrain(parameter #)

is also provided to simplify this.

## Setting parameter values
The function

    ROOT::Fit::Fitter* fitter = fit->GetFitter();

is provided for direct access to the ROOT::Fit::Fitter object. This allows to
set and fix parameter values, limits and set step sizes directly via

    fitter->Config().ParSettings(parameter #).Set(const std::string &name, double value, double step, double lower, double upper);

## Restricting the fit range
The fit range can be restricted through

    fit->SetRangeX(first bin #, last bin #);
and freed using

    fit->ReleaseRangeX();
For 2D histograms the Y range can be similarly restricted using

    fit->SetRangeY(first bin #, last bin #);
    fit->ReleaseRangeY();
and for 3D histograms also

    fit->SetRangeZ(first bin #, last bin #);
    fit->ReleaseRangeZ();
It is also possible to exclude individual bins from the fit through

    fit->ExcludeBin(bin #);
where the given bin number is assumed to follow the TH1::GetBin() numbering.
Any bins excluded in this way can be included again using the corresponding

    fit->IncludeBin(bin #);

## Weights histograms
Weights histograms (for a motivation see the above publication) can be specified
for the individual MC sources through

    fit->SetWeight(parameter #, pointer to weights histogram);
and unset by specifying a null pointer.

## Obtaining fit results
The fit is carried out through

    Int_t status = fit->Fit();
where  status  is the code returned from the "MINIMIZE" command. For fits
that converged, parameter values and errors can be obtained through

    fit->GetResult(parameter #, value, error);
and the histogram corresponding to the total Monte Carlo prediction (which
is not the same as a simple weighted sum of the input Monte Carlo distributions)
can be obtained by

    TH1* result = fit->GetPlot();
## Using different histograms
It is possible to change the histogram being fitted through

    fit->SetData(TH1* data);
and to change the template histogram for a given parameter number through

    fit->SetMC(parameter #, TH1* MC);
This can speed up code in case of multiple data or template histograms;
however, it should be done with care as any settings are taken over from
the previous fit. In addition, neither the dimensionality nor the numbers of
bins of the histograms should change (in that case it is better to instantiate
a new TFractionFitter object).

## Errors
Any serious inconsistency results in an error.
*/

#include "TH1.h"
#include "TMath.h"
#include "TClass.h"

#include "Fit/FitConfig.h"
#include "Fit/Fitter.h"
#include "TFitResult.h"
#include "Math/Functor.h"
#include "TFractionFitter.h"

ClassImp(TFractionFitter);


////////////////////////////////////////////////////////////////////////////////
/// TFractionFitter default constructor.

TFractionFitter::TFractionFitter() :
fFitDone(kFALSE),
fLowLimitX(0), fHighLimitX(0),
fLowLimitY(0), fHighLimitY(0),
fLowLimitZ(0), fHighLimitZ(0),
fData(0), fIntegralData(0),
fPlot(0)
{
   fFractionFitter = 0;
   fIntegralMCs   = 0;
   fFractions     = 0;

   fNpfits        = 0;
   fNDF           = 0;
   fChisquare     = 0;
   fNpar          = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// TFractionFitter constructor. Does a complete initialisation (including
/// consistency checks, default fit range as the whole histogram but without
/// under- and overflows, and declaration of the fit parameters). Note that
/// the histograms are not copied, only references are used.
/// \param[in] data histogram to be fitted
/// \param[in] MCs  array of TH1* corresponding template distributions
/// \param[in] option   can be used to control the print level of the minimization algorithm
///            - option = "Q"  : quite - no message is printed
///            - option = "V"  : verbose - max print out
///            - option = ""   : default: print initial fraction values and result

TFractionFitter::TFractionFitter(TH1* data, TObjArray  *MCs, Option_t *option) :
fFitDone(kFALSE), fChisquare(0), fPlot(0)  {
   fData = data;
   // Default: include all of the histogram (but without under- and overflows)
   fLowLimitX = 1;
   fHighLimitX = fData->GetNbinsX();
   if (fData->GetDimension() > 1) {
      fLowLimitY = 1;
      fHighLimitY = fData->GetNbinsY();
      if (fData->GetDimension() > 2) {
         fLowLimitZ = 1;
         fHighLimitZ = fData->GetNbinsZ();
      }
   }
   fNpar = MCs->GetEntries();
   Int_t par;
   for (par = 0; par < fNpar; ++par) {
      fMCs.Add(MCs->At(par));
      // Histogram containing template prediction
      TString s = Form("Prediction for MC sample %i",par);
      TH1* pred = (TH1*) ((TH1*)MCs->At(par))->Clone(s);
      // TFractionFitter manages these histograms
      pred->SetDirectory(0);
      pred->SetTitle(s);
      fAji.Add(pred);
   }
   fIntegralMCs = new Double_t[fNpar];
   fFractions = new Double_t[fNpar];

   CheckConsistency();
   fWeights.Expand(fNpar);

   fFractionFitter = new ROOT::Fit::Fitter();

   // set print level
   TString opt(option);
   opt.ToUpper();
   if (opt.Contains("Q") ) {
      fFractionFitter->Config().MinimizerOptions().SetPrintLevel(0);
   }
   else if (opt.Contains("V") ) {
      fFractionFitter->Config().MinimizerOptions().SetPrintLevel(2);
   }
   else
      fFractionFitter->Config().MinimizerOptions().SetPrintLevel(1);

   Double_t defaultFraction = 1.0/((Double_t)fNpar);
   Double_t defaultStep = 0.01;
   // set the parameters
   std::vector<ROOT::Fit::ParameterSettings> & parameters = fFractionFitter->Config().ParamsSettings();
   parameters.reserve(fNpar);
   for (par = 0; par < fNpar; ++par) {
      TString name("frac"); name += par;
      parameters.push_back(ROOT::Fit::ParameterSettings(name.Data(), defaultFraction, defaultStep) );
   }

   if (fFractionFitter->Config().MinimizerOptions().ErrorDef() == 1.0 )
      fFractionFitter->Config().MinimizerOptions().SetErrorDef(0.5);

}

////////////////////////////////////////////////////////////////////////////////
/// TFractionFitter default destructor

TFractionFitter::~TFractionFitter() {
   if (fFractionFitter) delete fFractionFitter;
   delete[] fIntegralMCs;
   delete[] fFractions;
   if (fPlot) delete fPlot;
   fAji.Delete();
}

////////////////////////////////////////////////////////////////////////////////
/// Change the histogram to be fitted to. Notes:
/// - Parameter constraints and settings are retained from a possible previous fit.
/// - Modifying the dimension or number of bins results in an error (in this case
///   rather instantiate a new TFractionFitter object)

void TFractionFitter::SetData(TH1* data) {
   fData = data;
   fFitDone = kFALSE;
   CheckConsistency();
}

////////////////////////////////////////////////////////////////////////////////
/// Change the histogram for template number <parm>. Notes:
/// - Parameter constraints and settings are retained from a possible previous fit.
/// - Modifying the dimension or number of bins results in an error (in this case
///   rather instantiate a new TFractionFitter object)

void TFractionFitter::SetMC(Int_t parm, TH1* MC) {
   CheckParNo(parm);
   fMCs.RemoveAt(parm);
   fMCs.AddAt(MC,parm);
   fFitDone = kFALSE;
   CheckConsistency();
}

////////////////////////////////////////////////////////////////////////////////
/// Set bin by bin weights for template number <parm> (the parameter numbering
/// follows that of the input template vector).
/// Weights can be "unset" by passing a null pointer.
/// Consistency of the weights histogram with the data histogram is checked at
/// this point, and an error in case of problems.

void TFractionFitter::SetWeight(Int_t parm, TH1* weight) {
   CheckParNo(parm);
   if (fWeights[parm]) {
      fWeights.RemoveAt(parm);
   }
   if (weight) {
      if (weight->GetNbinsX() != fData->GetNbinsX() ||
          (fData->GetDimension() > 1 && weight->GetNbinsY() != fData->GetNbinsY()) ||
          (fData->GetDimension() > 2 && weight->GetNbinsZ() != fData->GetNbinsZ())) {
         Error("SetWeight","Inconsistent weights histogram for source %d", parm);
         return;
      }
      TString ts = "weight hist: "; ts += weight->GetName();
      fWeights.AddAt(weight,parm);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Give direct access to the underlying fitter class. This can be
/// used e.g. to modify parameter values or step sizes.

ROOT::Fit::Fitter* TFractionFitter::GetFitter() const {
   return fFractionFitter;
}

////////////////////////////////////////////////////////////////////////////////
/// Function for internal use, checking parameter validity
/// An invalid parameter results in an error.

void TFractionFitter::CheckParNo(Int_t parm) const {
   if (parm < 0 || parm > fNpar) {
      Error("CheckParNo","Invalid parameter number %d",parm);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the X range of the histogram to be used in the fit.
/// Use ReleaseRangeX() to go back to fitting the full histogram.
/// The consistency check ensures that no empty fit range occurs (and also
/// recomputes the bin content integrals).
/// \param[in] low lower X bin number
/// \param[in] high upper X bin number

void TFractionFitter::SetRangeX(Int_t low, Int_t high) {
   fLowLimitX = (low > 0 ) ? low : 1;
   fHighLimitX = ( high > 0 && high <= fData->GetNbinsX()) ? high : fData->GetNbinsX();
   CheckConsistency();
}

////////////////////////////////////////////////////////////////////////////////
/// Release restrictions on the X range of the histogram to be used in the fit.

void TFractionFitter::ReleaseRangeX() {
   fLowLimitX = 1;
   fHighLimitX = fData->GetNbinsX();
   CheckConsistency();
}

////////////////////////////////////////////////////////////////////////////////
/// Set the Y range of the histogram to be used in the fit (2D or 3D histograms only).
/// Use ReleaseRangeY() to go back to fitting the full histogram.
/// The consistency check ensures that no empty fit range occurs (and also
/// recomputes the bin content integrals).
/// \param[in] low lower X bin number
/// \param[in] high upper X bin number

void TFractionFitter::SetRangeY(Int_t low, Int_t high) {
   if (fData->GetDimension() < 2) {
      Error("SetRangeY","Y range cannot be set for 1D histogram");
      return;
   }

   fLowLimitY = (low > 0) ? low : 1;
   fHighLimitY = (high > 0 && high <= fData->GetNbinsY()) ? high : fData->GetNbinsY();
   CheckConsistency();
}

////////////////////////////////////////////////////////////////////////////////
/// Release restrictions on the Y range of the histogram to be used in the fit.

void TFractionFitter::ReleaseRangeY() {
   fLowLimitY = 1;
   fHighLimitY = fData->GetNbinsY();
   CheckConsistency();
}


////////////////////////////////////////////////////////////////////////////////
/// Set the Z range of the histogram to be used in the fit (3D histograms only).
/// Use ReleaseRangeY() to go back to fitting the full histogram.
/// The consistency check ensures that no empty fit range occurs (and also
/// recomputes the bin content integrals).
/// \param[in] low lower X bin number
/// \param[in] high upper X bin number

void TFractionFitter::SetRangeZ(Int_t low, Int_t high) {
   if (fData->GetDimension() < 3) {
      Error("SetRangeZ","Z range cannot be set for 1D or 2D histogram");
      return;
   }


   fLowLimitZ = (low > 0) ? low : 1;
   fHighLimitZ = (high > 0 && high <= fData->GetNbinsZ()) ? high : fData->GetNbinsZ();
   CheckConsistency();
}

////////////////////////////////////////////////////////////////////////////////
/// Release restrictions on the Z range of the histogram to be used in the fit.

void TFractionFitter::ReleaseRangeZ() {
   fLowLimitZ = 1;
   fHighLimitZ = fData->GetNbinsZ();
   CheckConsistency();
}

////////////////////////////////////////////////////////////////////////////////
/// Exclude the given bin from the fit. The bin numbering to be used is that
/// of TH1::GetBin().

void TFractionFitter::ExcludeBin(Int_t bin) {
   int excluded = fExcludedBins.size();
   for (int b = 0; b < excluded; ++b) {
      if (fExcludedBins[b] == bin) {
         Error("ExcludeBin", "bin %d already excluded", bin);
         return;
      }
   }
   fExcludedBins.push_back(bin);
   // This call serves to properly (re)determine the number of degrees of freeom
   CheckConsistency();
}

////////////////////////////////////////////////////////////////////////////////
/// Include the given bin in the fit, if it was excluded before using ExcludeBin().
/// The bin numbering to be used is that of TH1::GetBin().

void TFractionFitter::IncludeBin(Int_t bin) {
   for (std::vector<Int_t>::iterator it = fExcludedBins.begin();
        it != fExcludedBins.end(); ++it) {
      if (*it == bin) {
         fExcludedBins.erase(it);
         // This call serves to properly (re)determine the number of degrees of freeom
         CheckConsistency();
         return;
      }
   }
   Error("IncludeBin", "bin %d was not excluded", bin);
}

////////////////////////////////////////////////////////////////////////////////
/// Function for internal use, checking whether the given bin is
/// excluded from the fit or not.

bool TFractionFitter::IsExcluded(Int_t bin) const {
   for (unsigned int b = 0; b < fExcludedBins.size(); ++b)
      if (fExcludedBins[b] == bin) return true;
   return false;
}

////////////////////////////////////////////////////////////////////////////////
/// Constrain the values of parameter number <parm> (the parameter numbering
/// follows that of the input template vector).
/// Use UnConstrain() to remove this constraint.

void TFractionFitter::Constrain(Int_t parm, Double_t low, Double_t high) {
   CheckParNo(parm);
   assert( parm >= 0 && parm <    (int) fFractionFitter->Config().ParamsSettings().size() );
   fFractionFitter->Config().ParSettings(parm).SetLimits(low,high);
}

////////////////////////////////////////////////////////////////////////////////
/// Remove the constraints on the possible values of parameter <parm>.

void TFractionFitter::UnConstrain(Int_t parm) {
   CheckParNo(parm);
   fFractionFitter->Config().ParSettings(parm).RemoveLimits();
}

////////////////////////////////////////////////////////////////////////////////
/// Function used internally to check the consistency between the
/// various histograms. Checks are performed on nonexistent or empty
/// histograms, the precise histogram class, and the number of bins.
/// In addition, integrals over the "allowed" bin ranges are computed.
/// Any inconsistency results in a error.

void TFractionFitter::CheckConsistency() {
   if (! fData) {
      Error("CheckConsistency","Nonexistent data histogram");
      return;
   }
   Int_t minX, maxX, minY, maxY, minZ, maxZ;
   Int_t x,y,z,par;
   GetRanges(minX, maxX, minY, maxY, minZ, maxZ);
   fIntegralData = 0;
   fNpfits = 0;
   for (z = minZ; z <= maxZ; ++z) {
      for (y = minY; y <= maxY; ++y) {
         for (x = minX; x <= maxX; ++x) {
            if (IsExcluded(fData->GetBin(x, y, z))) continue;
            fNpfits++;
            fIntegralData += fData->GetBinContent(x, y, z);
         }
      }
   }
   if (fIntegralData <= 0) {
      Error("CheckConsistency","Empty data histogram");
      return;
   }
   TClass* cl = fData->Class();

   fNDF = fNpfits - fNpar;

   if (fNpar < 2) {
      Error("CheckConsistency","Need at least two MC histograms");
      return;
   }

   for (par = 0; par < fNpar; ++par) {
      TH1 *h = (TH1*)fMCs.At(par);
      if (! h) {
         Error("CheckConsistency","Nonexistent MC histogram for source #%d",par);
         return;
      }
      if ((! h->Class()->InheritsFrom(cl)) || h->GetNbinsX() != fData->GetNbinsX() ||
          (fData->GetDimension() > 1 && h->GetNbinsY() != fData->GetNbinsY()) ||
          (fData->GetDimension() > 2 && h->GetNbinsZ() != fData->GetNbinsZ())) {
         Error("CheckConsistency","Histogram inconsistency for source #%d",par);
         return;
      }
      fIntegralMCs[par] = 0;
      for (z = minZ; z <= maxZ; ++z) {
         for (y = minY; y <= maxY; ++y) {
            for (x = minX; x <= maxX; ++x) {
               Int_t bin = fData->GetBin(x, y, z);
               if (IsExcluded(bin)) continue;
               Double_t MCEvents = h->GetBinContent(bin);
               if (MCEvents < 0) {
                  Error("CheckConsistency", "Number of MC events (bin = %d, par = %d) cannot be negative: "
                        " their distribution is binomial (see paper)", bin, par);
               }
               fIntegralMCs[par] += MCEvents;
            }
         }
      }
      if (fIntegralMCs[par] <= 0) {
         Error("CheckConsistency","Empty MC histogram #%d",par);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Perform the fit with the default UP value.
/// The value returned is the minimisation status.

TFitResultPtr TFractionFitter::Fit() {

   // remove any existing output histogram
   if (fPlot) {
      delete fPlot; fPlot = 0;
   }

   // Make sure the correct likelihood computation is used
   ROOT::Math::Functor fcnFunction(this, &TFractionFitter::EvaluateFCN, fNpar);
   fFractionFitter->SetFCN(static_cast<ROOT::Math::IMultiGenFunction&>(fcnFunction));

   // fit
   Bool_t status = fFractionFitter->FitFCN();
   if (!status)   Warning("Fit","Abnormal termination of minimization.");

   fFitDone = kTRUE;

   // determine goodness of fit
   ComputeChisquareLambda();

   // create a new result class
   TFitResult* fr = new TFitResult(fFractionFitter->Result());
   TString name = TString::Format("TFractionFitter_result_of_%s",fData->GetName() );
   fr->SetName(name); fr->SetTitle(name);
   return TFitResultPtr(fr);
}

////////////////////////////////////////////////////////////////////////////////
/// Set UP to the given value (see class TMinuit), and perform a MINOS minimisation.

void TFractionFitter::ErrorAnalysis(Double_t UP) {
   if (! fFitDone) {
      Error("ErrorAnalysis","Fit not yet performed");
      return;
   }


   Double_t up = UP > 0 ? UP : 0.5;
   fFractionFitter->Config().MinimizerOptions().SetErrorDef(up);
   Bool_t status = fFractionFitter->CalculateMinosErrors();
   if (!status) {
      Error("ErrorAnalysis","Error return from MINOS: %d",fFractionFitter->Result().Status());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Obtain the fit result for parameter <parm> (the parameter numbering
/// follows that of the input template vector).

void TFractionFitter::GetResult(Int_t parm, Double_t& value, Double_t& error) const {
   CheckParNo(parm);
   if (! fFitDone) {
      Error("GetResult","Fit not yet performed");
      return;
   }
   value = fFractionFitter->Result().Parameter(parm);
   error = fFractionFitter->Result().Error(parm);
}

////////////////////////////////////////////////////////////////////////////////
/// Return the "template prediction" corresponding to the fit result (this is not
/// the same as the weighted sum of template distributions, as template statistical
/// uncertainties are taken into account).
/// Note that the name of this histogram will simply be the same as that of the
/// "data" histogram, prefixed with the string "Fraction fit to hist: ".
/// Note also that the histogram is managed by the TFractionFitter class, so the returned pointer will be invalid if
/// the class is deleted

TH1* TFractionFitter::GetPlot() {
   if (! fFitDone) {
      Error("GetPlot","Fit not yet performed");
      return 0;
   }
   if (! fPlot) {
      Double_t f = 0;
      const Double_t * par = fFractionFitter->Result().GetParams();
      assert(par);
      ComputeFCN(f, par, 3);
   }
   return fPlot;
}

////////////////////////////////////////////////////////////////////////////////
/// Used internally to obtain the bin ranges according to the dimensionality of
/// the histogram and the limits set by hand.

void TFractionFitter::GetRanges(Int_t& minX, Int_t& maxX, Int_t& minY, Int_t& maxY,
                                Int_t& minZ, Int_t& maxZ) const {
   if (fData->GetDimension() < 2) {
      minY = maxY = minZ = maxZ = 0;
      minX = fLowLimitX;
      maxX = fHighLimitX;
   } else if (fData->GetDimension() < 3) {
      minZ = maxZ = 0;
      minX = fLowLimitX;
      maxX = fHighLimitX;
      minY = fLowLimitY;
      maxY = fHighLimitY;
   } else {
      minX = fLowLimitX;
      maxX = fHighLimitX;
      minY = fLowLimitY;
      maxY = fHighLimitY;
      minZ = fLowLimitZ;
      maxZ = fHighLimitZ;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Used internally to compute the likelihood value.

void TFractionFitter::ComputeFCN(Double_t& f, const Double_t* xx, Int_t flag)
{
   // normalise the fit parameters
   Int_t bin, mc;
   Int_t minX, maxX, minY, maxY, minZ, maxZ;
   Int_t x,y,z;
   GetRanges(minX, maxX, minY, maxY, minZ, maxZ);
   for (mc = 0; mc < fNpar; ++mc) {
      Double_t tot;
      TH1 *h  = (TH1*)fMCs[mc];
      TH1 *hw = (TH1*)fWeights[mc];
      if (hw) {
         tot = 0;
         for (z = minZ; z <= maxZ; ++z) {
            for (y = minY; y <= maxY; ++y) {
               for (x = minX; x <= maxX; ++x) {
                  if (IsExcluded(fData->GetBin(x, y, z))) continue;
                  Double_t weight = hw->GetBinContent(x, y, z);
                  if (weight <= 0) {
                     Error("ComputeFCN","Invalid weight encountered for MC source %d",mc);
                     return;
                  }
                  tot += weight * h->GetBinContent(x, y, z);
               }
            }
         }
      } else tot = fIntegralMCs[mc];
      fFractions[mc] = xx[mc] * fIntegralData / tot;
   }

   if (flag == 3) {
      TString ts = "Fraction fit to hist: "; ts += fData->GetName();
      fPlot = (TH1*) fData->Clone(ts.Data());
      // plot histogram is managed by TFractionFitter
      fPlot->SetDirectory(0);
      fPlot->Reset();
   }
   // likelihood computation
   Double_t result = 0;
   for (z = minZ; z <= maxZ; ++z) {
      for (y = minY; y <= maxY; ++y) {
         for (x = minX; x <= maxX; ++x) {
            bin = fData->GetBin(x, y, z);
            if (IsExcluded(bin)) continue;

            // Solve for the "predictions"
            int k0 = 0;
            Double_t ti = 0.0; Double_t aki = 0.0;
            FindPrediction(bin, ti, k0, aki);

            Double_t prediction = 0;
            for (mc = 0; mc < fNpar; ++mc) {
               TH1 *h  = (TH1*)fMCs[mc];
               TH1 *hw = (TH1*)fWeights[mc];
               Double_t binPrediction;
               Double_t binContent = h->GetBinContent(bin);
               Double_t weight = hw ? hw->GetBinContent(bin) : 1;
               if (k0 >= 0 && fFractions[mc] == fFractions[k0]) {
                  binPrediction = aki;
               } else {
                  binPrediction = binContent > 0 ? binContent / (1+weight*fFractions[mc]*ti) : 0;
               }

               prediction += fFractions[mc]*weight*binPrediction;
               result -= binPrediction;
               if (binContent > 0 && binPrediction > 0)
                  result += binContent*TMath::Log(binPrediction);

               if (flag == 3) {
                  ((TH1*)fAji.At(mc))->SetBinContent(bin, binPrediction);
               }
            }

            if (flag == 3) {
               fPlot->SetBinContent(bin, prediction);
            }

            result -= prediction;
            Double_t found = fData->GetBinContent(bin);
            if (found > 0 && prediction > 0)
               result += found*TMath::Log(prediction);
         }
      }
   }

   f = -result;
}

////////////////////////////////////////////////////////////////////////////////
/// Function used internally to obtain the template prediction in the individual bins
/// 'bin' <=> 'i' (paper)
/// 'par' <=> 'j' (paper)

void TFractionFitter::FindPrediction(int bin, Double_t &t_i, int& k_0, Double_t &A_ki) const {
   std::vector<Double_t> wgtFrac(fNpar); // weighted fractions (strengths of the sources)
   std::vector<Double_t> a_ji(fNpar); // number of observed MC events for bin 'i' and par (source) 'j'
   Double_t d_i = fData->GetBinContent(bin); // number of events in the real data for bin 'i'

   // Cache the weighted fractions and the number of observed MC events
   // Sanity check: none of the fractions should be == 0
   for (Int_t par = 0; par < fNpar; ++par) {
      a_ji[par] = ((TH1*)fMCs.At(par))->GetBinContent(bin);
      TH1* hw = (TH1*)fWeights.At(par);
      wgtFrac[par] = hw ? hw->GetBinContent(bin) * fFractions[par] : fFractions[par];
      if (wgtFrac[par] == 0) {
         Error("FindPrediction", "Fraction[%d] = 0!", par);
         return;
      }
   }

   // Case data = 0
   if (TMath::Nint(d_i) == 0) {
      t_i = 1;
      k_0 = -1;
      A_ki = 0;
      return;
   }

   // Case one or more of the MC bin contents == 0 -> find largest fraction
   // k_0 stores the source index of the largest fraction
   k_0 = 0;
   Double_t maxWgtFrac = wgtFrac[0];
   for (Int_t par = 1; par < fNpar; ++par) {
      if (wgtFrac[par] > maxWgtFrac) {
         k_0 = par;
         maxWgtFrac = wgtFrac[par];
      }
   }
   Double_t t_min = -1 / maxWgtFrac; // t_i cannot be smaller than this value (see paper, par 5)

   // Determine if there are more sources which have the same maximum contribution (fraction)
   Int_t nMax = 1; Double_t contentsMax = a_ji[k_0];
   for (Int_t par = 0; par < fNpar; ++par) {
      if (par == k_0) continue;
      if (wgtFrac[par] == maxWgtFrac) {
         nMax++;
         contentsMax += a_ji[par];
      }
   }

   // special action if there is a zero in the number of entries for the MC source with
   // the largest strength (fraction) -> See Paper, Paragraph 5
   if (contentsMax == 0) {
      A_ki = d_i / (1.0 + maxWgtFrac);
      for (Int_t par = 0; par < fNpar; ++par) {
         if (par == k_0 || wgtFrac[par] == maxWgtFrac) continue;
         A_ki -= a_ji[par] * wgtFrac[par] / (maxWgtFrac - wgtFrac[par]);
      }
      if (A_ki > 0) {
         A_ki /= nMax;
         t_i = t_min;
         return;
      }
   }
   k_0 = -1;

   // Case of nonzero histogram contents: solve for t_i using Newton's method
   // The equation that needs to be solved:
   //    func(t_i) = \sum\limits_j{\frac{ p_j a_{ji} }{1 + p_j t_i}} - \frac{d_i}{1 - t_i} = 0
   t_i = 0; Double_t step = 0.2;
   Int_t maxIter = 100000; // maximum number of iterations
   for(Int_t i = 0; i < maxIter; ++i) {
      if (t_i >= 1 || t_i < t_min) {
         step /= 10;
         t_i = 0;
      }
      Double_t func   = - d_i / (1.0 - t_i);
      Double_t deriv = func / (1.0 - t_i);
      for (Int_t par = 0; par < fNpar; ++par) {
         Double_t r = 1.0 / (t_i + 1.0 / wgtFrac[par]);
         func   += a_ji[par] * r;
         deriv -= a_ji[par] * r * r;
      }
      if (TMath::Abs(func) < 1e-12) return; // solution found
      Double_t delta = - func / deriv; // update delta
      if (TMath::Abs(delta) > step)
         delta = (delta > 0) ? step : -step; // correct delta if it becomes too large
      t_i += delta;
      if (TMath::Abs(delta) < 1e-13) return; // solution found
   } // the loop breaks when the solution is found, or when the maximum number of iterations is exhausted

   Warning("FindPrediction", "Did not find solution for t_i in %d iterations", maxIter);

   return;
}

#ifdef OLD
////////////////////////////////////////////////////////////////////////////////
/// Function called by the minimisation package. The actual functionality is passed
/// on to the TFractionFitter::ComputeFCN member function.

void TFractionFitFCN(Int_t& npar, Double_t* gin, Double_t& f, Double_t* par, Int_t flag) {
   TFractionFitter* fitter = dynamic_cast<TFractionFitter*>(fFractionFitter->GetObjectFit());
   if (!fitter) {
      Error("TFractionFitFCN","Invalid fit object encountered!");
      return;
   }
   fitter->ComputeFCN(npar, gin, f, par, flag);
}
#endif

////////////////////////////////////////////////////////////////////////////////
/// Return the likelihood ratio Chi-squared (chi2) for the fit.
/// The value is computed when the fit is executed successfully.
/// Chi2 calculation is based on the "likelihood ratio" lambda,
/// lambda = L(y;n) / L(m;n),
/// where L(y;n) is the likelihood of the fit result <y> describing the data <n>
/// and L(m;n) is the likelihood of an unknown "true" underlying distribution
/// <m> describing the data <n>. Since <m> is unknown, the data distribution is
/// used instead,
/// lambda = L(y;n) / L(n;n).
/// Note that this ratio is 1 if the fit is perfect. The chi2 value is then
/// computed according to
/// chi2 = -2*ln(lambda).
/// This parameter can be shown to follow a Chi-square distribution. See for
/// example S. Baker and R. Cousins, "Clarification of the use of chi-square
/// and likelihood functions in fits to histograms", Nucl. Instr. Meth. A221,
/// pp. 437-442 (1984)

Double_t TFractionFitter::GetChisquare() const
{
   return fChisquare;
}

////////////////////////////////////////////////////////////////////////////////
/// return the number of degrees of freedom in the fit
/// the fNDF parameter has been previously computed during a fit.
/// The number of degrees of freedom corresponds to the number of points
/// used in the fit minus the number of templates.

Int_t TFractionFitter::GetNDF() const
{
   if (fNDF == 0) return fNpfits-fNpar;
   return fNDF;
}

////////////////////////////////////////////////////////////////////////////////
/// return the fit probability

Double_t TFractionFitter::GetProb() const
{
   Int_t ndf = fNpfits - fNpar;
   if (ndf <= 0) return 0;
   return TMath::Prob(fChisquare,ndf);
}

////////////////////////////////////////////////////////////////////////////////
/// Method used internally to compute the likelihood ratio chi2
/// See the function GetChisquare() for details

void TFractionFitter::ComputeChisquareLambda()
{
   if ( !fFitDone ) {
      Error("ComputeChisquareLambda","Fit not yet (successfully) performed");
      fChisquare = 0;
      return;
   }

   // fPlot must be initialized and filled. Leave this to the GetPlot() method.
   if (! fPlot)
      GetPlot();

   Int_t minX, maxX, minY, maxY, minZ, maxZ;
   GetRanges(minX, maxX, minY, maxY, minZ, maxZ);

   Double_t logLyn = 0; // likelihood of prediction
   Double_t logLmn = 0; // likelihood of data ("true" distribution)
   for(Int_t x = minX; x <= maxX; x++) {
      for(Int_t y = minY; y <= maxY; y++) {
         for(Int_t z = minZ; z <= maxZ; z++) {
            if (IsExcluded(fData->GetBin(x, y, z))) continue;
            Double_t di = fData->GetBinContent(x, y, z);
            Double_t fi = fPlot->GetBinContent(x, y, z);
            if(fi != 0) logLyn += di * TMath::Log(fi) - fi;
            if(di != 0) logLmn += di * TMath::Log(di) - di;
            for(Int_t j = 0; j < fNpar; j++) {
               Double_t aji = ((TH1*)fMCs.At(j))->GetBinContent(x, y, z);
               Double_t bji = ((TH1*)fAji.At(j))->GetBinContent(x, y, z);
               if(bji != 0) logLyn += aji * TMath::Log(bji) - bji;
               if(aji != 0) logLmn += aji * TMath::Log(aji) - aji;
            }
         }
      }
   }

   fChisquare = -2*logLyn + 2*logLmn;

   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the adjusted MC template (Aji) for template (parm).
/// Note that the (Aji) times fractions only sum to the total prediction
/// of the fit if all weights are 1.
/// Note also that the histogram is managed by the TFractionFitter class, so the returned pointer will be invalid if
/// the class is deleted

TH1* TFractionFitter::GetMCPrediction(Int_t parm) const
{
   CheckParNo(parm);
   if ( !fFitDone ) {
      Error("GetMCPrediction","Fit not yet performed");
      return 0;
   }
   return (TH1*) fAji.At(parm);
}
