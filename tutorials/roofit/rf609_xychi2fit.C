/// \file
/// \ingroup tutorial_roofit
/// \notebook -js
/// Likelihood and minimization: setting up a chi^2 fit to an unbinned dataset with X,Y,err(Y)
/// values (and optionally err(X) values)
///
/// \macro_image
/// \macro_output
/// \macro_code
///
/// \date July 2008
/// \author Wouter Verkerke

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooPolyVar.h"
#include "RooConstVar.h"
#include "RooChi2Var.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TRandom.h"

using namespace RooFit;

void rf609_xychi2fit()
{
   // C r e a t e   d a t a s e t   w i t h   X   a n d   Y   v a l u e s
   // -------------------------------------------------------------------

   // Make weighted XY dataset with asymmetric errors stored
   // The StoreError() argument is essential as it makes
   // the dataset store the error in addition to the values
   // of the observables. If errors on one or more observables
   // are asymmetric, one can store the asymmetric error
   // using the StoreAsymError() argument

   RooRealVar x("x", "x", -11, 11);
   RooRealVar y("y", "y", -10, 200);
   RooDataSet dxy("dxy", "dxy", RooArgSet(x, y), StoreError(RooArgSet(x, y)));

   // Fill an example dataset with X,err(X),Y,err(Y) values
   for (int i = 0; i <= 10; i++) {

      // Set X value and error
      x = -10 + 2 * i;
      x.setError(i < 5 ? 0.5 / 1. : 1.0 / 1.);

      // Set Y value and error
      y = x.getVal() * x.getVal() + 4 * fabs(gRandom->Gaus());
      y.setError(sqrt(y.getVal()));

      dxy.add(RooArgSet(x, y));
   }

   // P e r f o r m   c h i 2   f i t   t o   X + / - d x   a n d   Y + / - d Y   v a l u e s
   // ---------------------------------------------------------------------------------------

   // Make fit function
   RooRealVar a("a", "a", 0.0, -10, 10);
   RooRealVar b("b", "b", 0.0, -100, 100);
   RooPolyVar f("f", "f", x, RooArgList(b, a, RooConst(1)));

   // Plot dataset in X-Y interpretation
   RooPlot *frame = x.frame(Title("Chi^2 fit of function set of (X#pmdX,Y#pmdY) values"));
   dxy.plotOnXY(frame, YVar(y));

   // Fit chi^2 using X and Y errors
   f.chi2FitTo(dxy, YVar(y));

   // Overlay fitted function
   f.plotOn(frame);

   // Alternative: fit chi^2 integrating f(x) over ranges defined by X errors, rather
   // than taking point at center of bin
   f.chi2FitTo(dxy, YVar(y), Integrate(kTRUE));

   // Overlay alternate fit result
   f.plotOn(frame, LineStyle(kDashed), LineColor(kRed));

   // Draw the plot on a canvas
   new TCanvas("rf609_xychi2fit", "rf609_xychi2fit", 600, 600);
   gPad->SetLeftMargin(0.15);
   frame->GetYaxis()->SetTitleOffset(1.4);
   frame->Draw();
}
