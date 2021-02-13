/// \file
/// \ingroup tutorial_roofit
/// \notebook -js
/// Data and categories: using weights in unbinned datasets
///
/// \macro_image
/// \macro_output
/// \macro_code
///
/// \date July 2008
/// \author Wouter Verkerke

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooFormulaVar.h"
#include "RooGenericPdf.h"
#include "RooPolynomial.h"
#include "RooChi2Var.h"
#include "RooMinimizer.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooFitResult.h"
using namespace RooFit;

void rf403_weightedevts()
{
   // C r e a t e   o b s e r v a b l e   a n d   u n w e i g h t e d   d a t a s e t
   // -------------------------------------------------------------------------------

   // Declare observable
   RooRealVar x("x", "x", -10, 10);
   x.setBins(40);

   // Construction a uniform pdf
   RooPolynomial p0("px", "px", x);

   // Sample 1000 events from pdf
   RooDataSet *data = p0.generate(x, 1000);

   // C a l c u l a t e   w e i g h t   a n d   m a k e   d a t a s e t   w e i g h t e d
   // -----------------------------------------------------------------------------------

   // Construct formula to calculate (fake) weight for events
   RooFormulaVar wFunc("w", "event weight", "(x*x+10)", x);

   // Add column with variable w to previously generated dataset
   RooRealVar *w = (RooRealVar *)data->addColumn(wFunc);

   // Dataset d is now a dataset with two observable (x,w) with 1000 entries
   data->Print();

   // Instruct dataset wdata in interpret w as event weight rather than as observable
   RooDataSet wdata(data->GetName(), data->GetTitle(), data, *data->get(), 0, w->GetName());

   // Dataset d is now a dataset with one observable (x) with 1000 entries and a sum of weights of ~430K
   wdata.Print();

   // U n b i n n e d   M L   f i t   t o   w e i g h t e d   d a t a
   // ---------------------------------------------------------------

   // Construction quadratic polynomial pdf for fitting
   RooRealVar a0("a0", "a0", 1);
   RooRealVar a1("a1", "a1", 0, -1, 1);
   RooRealVar a2("a2", "a2", 1, 0, 10);
   RooPolynomial p2("p2", "p2", x, RooArgList(a0, a1, a2), 0);

   // Fit quadratic polynomial to weighted data

   // NOTE: A plain Maximum likelihood fit to weighted data does in general
   //       NOT result in correct error estimates, unless individual
   //       event weights represent Poisson statistics themselves.
   //
   // Fit with 'wrong' errors
   RooFitResult *r_ml_wgt = p2.fitTo(wdata, Save());

   // A first order correction to estimated parameter errors in an
   // (unbinned) ML fit can be obtained by calculating the
   // covariance matrix as
   //
   //    V' = V C-1 V
   //
   // where V is the covariance matrix calculated from a fit
   // to -logL = - sum [ w_i log f(x_i) ] and C is the covariance
   // matrix calculated from -logL' = -sum [ w_i^2 log f(x_i) ]
   // (i.e. the weights are applied squared)
   //
   // A fit in this mode can be performed as follows:

   RooFitResult *r_ml_wgt_corr = p2.fitTo(wdata, Save(), SumW2Error(kTRUE));

   // P l o t   w e i g h e d   d a t a   a n d   f i t   r e s u l t
   // ---------------------------------------------------------------

   // Construct plot frame
   RooPlot *frame = x.frame(Title("Unbinned ML fit, binned chi^2 fit to weighted data"));

   // Plot data using sum-of-weights-squared error rather than Poisson errors
   wdata.plotOn(frame, DataError(RooAbsData::SumW2));

   // Overlay result of 2nd order polynomial fit to weighted data
   p2.plotOn(frame);

   // ML Fit of pdf to equivalent unweighted dataset
   // -----------------------------------------------------------------------------------------

   // Construct a pdf with the same shape as p0 after weighting
   RooGenericPdf genPdf("genPdf", "x*x+10", x);

   // Sample a dataset with the same number of events as data
   RooDataSet *data2 = genPdf.generate(x, 1000);

   // Sample a dataset with the same number of weights as data
   RooDataSet *data3 = genPdf.generate(x, 43000);

   // Fit the 2nd order polynomial to both unweighted datasets and save the results for comparison
   RooFitResult *r_ml_unw10 = p2.fitTo(*data2, Save());
   RooFitResult *r_ml_unw43 = p2.fitTo(*data3, Save());

   // C h i 2   f i t   o f   p d f   t o   b i n n e d   w e i g h t e d   d a t a s e t
   // ------------------------------------------------------------------------------------

   // Construct binned clone of unbinned weighted dataset
   RooDataHist *binnedData = wdata.binnedClone();
   binnedData->Print("v");

   // Perform chi2 fit to binned weighted dataset using sum-of-weights errors
   //
   // NB: Within the usual approximations of a chi2 fit, a chi2 fit to weighted
   // data using sum-of-weights-squared errors does give correct error
   // estimates
   RooChi2Var chi2("chi2", "chi2", p2, *binnedData, DataError(RooAbsData::SumW2));
   RooMinimizer m(chi2);
   m.migrad();
   m.hesse();

   // Plot chi^2 fit result on frame as well
   RooFitResult *r_chi2_wgt = m.save();
   p2.plotOn(frame, LineStyle(kDashed), LineColor(kRed));

   // C o m p a r e   f i t   r e s u l t s   o f   c h i 2 , M L   f i t s   t o   ( u n ) w e i g h t e d   d a t a
   // ---------------------------------------------------------------------------------------------------------------

   // Note that ML fit on 1Kevt of weighted data is closer to result of ML fit on 43Kevt of unweighted data
   // than to 1Kevt of unweighted data, whereas the reference chi^2 fit with SumW2 error gives a result closer to
   // that of an unbinned ML fit to 1Kevt of unweighted data.

   cout << "==> ML Fit results on 1K unweighted events" << endl;
   r_ml_unw10->Print();
   cout << "==> ML Fit results on 43K unweighted events" << endl;
   r_ml_unw43->Print();
   cout << "==> ML Fit results on 1K weighted events with a summed weight of 43K" << endl;
   r_ml_wgt->Print();
   cout << "==> Corrected ML Fit results on 1K weighted events with a summed weight of 43K" << endl;
   r_ml_wgt_corr->Print();
   cout << "==> Chi2 Fit results on 1K weighted events with a summed weight of 43K" << endl;
   r_chi2_wgt->Print();

   new TCanvas("rf403_weightedevts", "rf403_weightedevts", 600, 600);
   gPad->SetLeftMargin(0.15);
   frame->GetYaxis()->SetTitleOffset(1.8);
   frame->Draw();
}
