/// \file
/// \ingroup tutorial_roofit
/// \notebook -js
/// Special pdf's: unbinned maximum likelihood fit of an efficiency eff(x) function
///
/// to a dataset D(x,cut), where cut is a category encoding a selection, of which the
/// efficiency as function of x should be described by eff(x)
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
#include "RooFormulaVar.h"
#include "RooProdPdf.h"
#include "RooEfficiency.h"
#include "RooPolynomial.h"
#include "RooCategory.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
using namespace RooFit;

void rf701_efficiencyfit()
{
   // C o n s t r u c t   e f f i c i e n c y   f u n c t i o n   e ( x )
   // -------------------------------------------------------------------

   // Declare variables x,mean,sigma with associated name, title, initial value and allowed range
   RooRealVar x("x", "x", -10, 10);

   // Efficiency function eff(x;a,b)
   RooRealVar a("a", "a", 0.4, 0, 1);
   RooRealVar b("b", "b", 5);
   RooRealVar c("c", "c", -1, -10, 10);
   RooFormulaVar effFunc("effFunc", "(1-a)+a*cos((x-c)/b)", RooArgList(a, b, c, x));

   // C o n s t r u c t   c o n d i t i o n a l    e f f i c i e n c y   p d f   E ( c u t | x )
   // ------------------------------------------------------------------------------------------

   // Acceptance state cut (1 or 0)
   RooCategory cut("cut", "cutr", { {"accept", 1}, {"reject", 0} });

   // Construct efficiency pdf eff(cut|x)
   RooEfficiency effPdf("effPdf", "effPdf", effFunc, cut, "accept");

   // G e n e r a t e   d a t a   ( x ,   c u t )   f r o m   a   t o y   m o d e l
   // -----------------------------------------------------------------------------

   // Construct global shape pdf shape(x) and product model(x,cut) = eff(cut|x)*shape(x)
   // (These are _only_ needed to generate some toy MC here to be used later)
   RooPolynomial shapePdf("shapePdf", "shapePdf", x, RooConst(-0.095));
   RooProdPdf model("model", "model", shapePdf, Conditional(effPdf, cut));

   // Generate some toy data from model
   RooDataSet *data = model.generate(RooArgSet(x, cut), 10000);

   // F i t   c o n d i t i o n a l   e f f i c i e n c y   p d f   t o   d a t a
   // --------------------------------------------------------------------------

   // Fit conditional efficiency pdf to data
   effPdf.fitTo(*data, ConditionalObservables(x));

   // P l o t   f i t t e d ,   d a t a   e f f i c i e n c y
   // --------------------------------------------------------

   // Plot distribution of all events and accepted fraction of events on frame
   RooPlot *frame1 = x.frame(Bins(20), Title("Data (all, accepted)"));
   data->plotOn(frame1);
   data->plotOn(frame1, Cut("cut==cut::accept"), MarkerColor(kRed), LineColor(kRed));

   // Plot accept/reject efficiency on data overlay fitted efficiency curve
   RooPlot *frame2 = x.frame(Bins(20), Title("Fitted efficiency"));
   data->plotOn(frame2, Efficiency(cut)); // needs ROOT version >= 5.21
   effFunc.plotOn(frame2, LineColor(kRed));

   // Draw all frames on a canvas
   TCanvas *ca = new TCanvas("rf701_efficiency", "rf701_efficiency", 800, 400);
   ca->Divide(2);
   ca->cd(1);
   gPad->SetLeftMargin(0.15);
   frame1->GetYaxis()->SetTitleOffset(1.6);
   frame1->Draw();
   ca->cd(2);
   gPad->SetLeftMargin(0.15);
   frame2->GetYaxis()->SetTitleOffset(1.4);
   frame2->Draw();
}
