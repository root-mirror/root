/// \file
/// \ingroup tutorial_roofit
/// \notebook
/// Multidimensional models: marginizalization of multi-dimensional pdfs through integration
///
/// \macro_image
/// \macro_output
/// \macro_code
///
/// \date July 2008
/// \author Wouter Verkerke

#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooPolyVar.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooNumIntConfig.h"
#include "RooConstVar.h"
using namespace RooFit;

void rf315_projectpdf()
{
   // C r e a t e   p d f   m ( x , y )  =  g x ( x | y ) * g ( y )
   // --------------------------------------------------------------

   // Increase default precision of numeric integration
   // as this exercise has high sensitivity to numeric integration precision
   RooAbsPdf::defaultIntegratorConfig()->setEpsRel(1e-8);
   RooAbsPdf::defaultIntegratorConfig()->setEpsAbs(1e-8);

   // Create observables
   RooRealVar x("x", "x", -5, 5);
   RooRealVar y("y", "y", -2, 2);

   // Create function f(y) = a0 + a1*y
   RooRealVar a0("a0", "a0", 0);
   RooRealVar a1("a1", "a1", -1.5, -3, 1);
   RooPolyVar fy("fy", "fy", y, RooArgSet(a0, a1));

   // Create gaussx(x,f(y),sx)
   RooRealVar sigmax("sigmax", "width of gaussian", 0.5);
   RooGaussian gaussx("gaussx", "Gaussian in x with shifting mean in y", x, fy, sigmax);

   // Create gaussy(y,0,2)
   RooGaussian gaussy("gaussy", "Gaussian in y", y, RooConst(0), RooConst(2));

   // Create gaussx(x,sx|y) * gaussy(y)
   RooProdPdf model("model", "gaussx(x|y)*gaussy(y)", gaussy, Conditional(gaussx, x));

   // M a r g i n a l i z e   m ( x , y )   t o   m ( x )
   // ----------------------------------------------------

   // modelx(x) = Int model(x,y) dy
   RooAbsPdf *modelx = model.createProjection(y);

   // U s e   m a r g i n a l i z e d   p . d . f .   a s   r e g u l a r   1 - D   p . d . f .
   // ------------------------------------------------------------------------------------------

   // Sample 1000 events from modelx
   RooAbsData *data = modelx->generateBinned(x, 1000);

   // Fit modelx to toy data
   modelx->fitTo(*data, Verbose());

   // Plot modelx over data
   RooPlot *frame = x.frame(40);
   data->plotOn(frame);
   modelx->plotOn(frame);

   // Make 2D histogram of model(x,y)
   TH1 *hh = model.createHistogram("x,y");
   hh->SetLineColor(kBlue);

   TCanvas *c = new TCanvas("rf315_projectpdf", "rf315_projectpdf", 800, 400);
   c->Divide(2);
   c->cd(1);
   gPad->SetLeftMargin(0.15);
   frame->GetYaxis()->SetTitleOffset(1.4);
   frame->Draw();
   c->cd(2);
   gPad->SetLeftMargin(0.20);
   hh->GetZaxis()->SetTitleOffset(2.5);
   hh->Draw("surf");
}
