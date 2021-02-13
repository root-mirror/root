/// \file
/// \ingroup tutorial_roofit
/// \notebook
/// Special pdf's: linear interpolation between pdf shapes using the 'Alex Read' algorithm
///
/// \macro_image
/// \macro_output
/// \macro_code
///
/// \date July 2008
/// \author Wouter Verkerke

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooPolynomial.h"
#include "RooIntegralMorph.h"
#include "RooNLLVar.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TH1.h"
using namespace RooFit;

void rf705_linearmorph()
{
   // C r e a t e   e n d   p o i n t   p d f   s h a p e s
   // ------------------------------------------------------

   // Observable
   RooRealVar x("x", "x", -20, 20);

   // Lower end point shape: a Gaussian
   RooRealVar g1mean("g1mean", "g1mean", -10);
   RooGaussian g1("g1", "g1", x, g1mean, RooConst(2));

   // Upper end point shape: a Polynomial
   RooPolynomial g2("g2", "g2", x, RooArgSet(RooConst(-0.03), RooConst(-0.001)));

   // C r e a t e   i n t e r p o l a t i n g   p d f
   // -----------------------------------------------

   // Create interpolation variable
   RooRealVar alpha("alpha", "alpha", 0, 1.0);

   // Specify sampling density on observable and interpolation variable
   x.setBins(1000, "cache");
   alpha.setBins(50, "cache");

   // Construct interpolating pdf in (x,a) represent g1(x) at a=a_min
   // and g2(x) at a=a_max
   RooIntegralMorph lmorph("lmorph", "lmorph", g1, g2, x, alpha);

   // P l o t   i n t e r p o l a t i n g   p d f   a t   v a r i o u s   a l p h a
   // -----------------------------------------------------------------------------

   // Show end points as blue curves
   RooPlot *frame1 = x.frame();
   g1.plotOn(frame1);
   g2.plotOn(frame1);

   // Show interpolated shapes in red
   alpha.setVal(0.125);
   lmorph.plotOn(frame1, LineColor(kRed));
   alpha.setVal(0.25);
   lmorph.plotOn(frame1, LineColor(kRed));
   alpha.setVal(0.375);
   lmorph.plotOn(frame1, LineColor(kRed));
   alpha.setVal(0.50);
   lmorph.plotOn(frame1, LineColor(kRed));
   alpha.setVal(0.625);
   lmorph.plotOn(frame1, LineColor(kRed));
   alpha.setVal(0.75);
   lmorph.plotOn(frame1, LineColor(kRed));
   alpha.setVal(0.875);
   lmorph.plotOn(frame1, LineColor(kRed));
   alpha.setVal(0.95);
   lmorph.plotOn(frame1, LineColor(kRed));

   // S h o w   2 D   d i s t r i b u t i o n   o f   p d f ( x , a l p h a )
   // -----------------------------------------------------------------------

   // Create 2D histogram
   TH1 *hh = lmorph.createHistogram("hh", x, Binning(40), YVar(alpha, Binning(40)));
   hh->SetLineColor(kBlue);

   // F i t   p d f   t o   d a t a s e t   w i t h   a l p h a = 0 . 8
   // -----------------------------------------------------------------

   // Generate a toy dataset at alpha = 0.8
   alpha = 0.8;
   RooDataSet *data = lmorph.generate(x, 1000);

   // Fit pdf to toy data
   lmorph.setCacheAlpha(kTRUE);
   lmorph.fitTo(*data, Verbose(kTRUE));

   // Plot fitted pdf and data overlaid
   RooPlot *frame2 = x.frame(Bins(100));
   data->plotOn(frame2);
   lmorph.plotOn(frame2);

   // S c a n   - l o g ( L )   v s   a l p h a
   // -----------------------------------------

   // Show scan -log(L) of dataset w.r.t alpha
   RooPlot *frame3 = alpha.frame(Bins(100), Range(0.1, 0.9));

   // Make 2D pdf of histogram
   RooNLLVar nll("nll", "nll", lmorph, *data);
   nll.plotOn(frame3, ShiftToZero());

   lmorph.setCacheAlpha(kFALSE);

   TCanvas *c = new TCanvas("rf705_linearmorph", "rf705_linearmorph", 800, 800);
   c->Divide(2, 2);
   c->cd(1);
   gPad->SetLeftMargin(0.15);
   frame1->GetYaxis()->SetTitleOffset(1.6);
   frame1->Draw();
   c->cd(2);
   gPad->SetLeftMargin(0.20);
   hh->GetZaxis()->SetTitleOffset(2.5);
   hh->Draw("surf");
   c->cd(3);
   gPad->SetLeftMargin(0.15);
   frame3->GetYaxis()->SetTitleOffset(1.4);
   frame3->Draw();
   c->cd(4);
   gPad->SetLeftMargin(0.15);
   frame2->GetYaxis()->SetTitleOffset(1.4);
   frame2->Draw();

   return;
}
