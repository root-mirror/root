/// \file
/// \ingroup tutorial_roofit
/// \notebook
/// Special pdf's: unbinned maximum likelihood fit of an efficiency eff(x) function
/// to a dataset D(x,cut), cut is a category encoding a selection whose efficiency as function
/// of x should be described by eff(x)
///
/// \macro_image
/// \macro_output
/// \macro_code
///
/// \date February 2018
/// \authors Clemens Lange, Wouter Verkerke (C++ version)

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooCategory.h"
#include "RooEfficiency.h"
#include "RooPolynomial.h"
#include "RooProdPdf.h"
#include "RooFormulaVar.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
#include "RooPlot.h"
using namespace RooFit;

void rf702_efficiencyfit_2D(Bool_t flat = kFALSE)
{
   // C o n s t r u c t   e f f i c i e n c y   f u n c t i o n   e ( x , y )
   // -----------------------------------------------------------------------

   // Declare variables x,mean,sigma with associated name, title, initial value and allowed range
   RooRealVar x("x", "x", -10, 10);
   RooRealVar y("y", "y", -10, 10);

   // Efficiency function eff(x;a,b)
   RooRealVar ax("ax", "ay", 0.6, 0, 1);
   RooRealVar bx("bx", "by", 5);
   RooRealVar cx("cx", "cy", -1, -10, 10);

   RooRealVar ay("ay", "ay", 0.2, 0, 1);
   RooRealVar by("by", "by", 5);
   RooRealVar cy("cy", "cy", -1, -10, 10);

   RooFormulaVar effFunc("effFunc", "((1-ax)+ax*cos((x-cx)/bx))*((1-ay)+ay*cos((y-cy)/by))",
                         RooArgList(ax, bx, cx, x, ay, by, cy, y));

   // Acceptance state cut (1 or 0)
   RooCategory cut("cut", "cutr", { {"accept", 1}, {"reject", 0} });

   // C o n s t r u c t   c o n d i t i o n a l    e f f i c i e n c y   p d f   E ( c u t | x , y )
   // ---------------------------------------------------------------------------------------------

   // Construct efficiency pdf eff(cut|x)
   RooEfficiency effPdf("effPdf", "effPdf", effFunc, cut, "accept");

   // G e n e r a t e   d a t a   ( x , y , c u t )   f r o m   a   t o y   m o d e l
   // -------------------------------------------------------------------------------

   // Construct global shape pdf shape(x) and product model(x,cut) = eff(cut|x)*shape(x)
   // (These are _only_ needed to generate some toy MC here to be used later)
   RooPolynomial shapePdfX("shapePdfX", "shapePdfX", x, RooConst(flat ? 0 : -0.095));
   RooPolynomial shapePdfY("shapePdfY", "shapePdfY", y, RooConst(flat ? 0 : +0.095));
   RooProdPdf shapePdf("shapePdf", "shapePdf", RooArgSet(shapePdfX, shapePdfY));
   RooProdPdf model("model", "model", shapePdf, Conditional(effPdf, cut));

   // Generate some toy data from model
   RooDataSet *data = model.generate(RooArgSet(x, y, cut), 10000);

   // F i t   c o n d i t i o n a l   e f f i c i e n c y   p d f   t o   d a t a
   // --------------------------------------------------------------------------

   // Fit conditional efficiency pdf to data
   effPdf.fitTo(*data, ConditionalObservables(RooArgSet(x, y)));

   // P l o t   f i t t e d ,   d a t a   e f f i c i e n c y
   // --------------------------------------------------------

   // Make 2D histograms of all data, selected data and efficiency function
   TH1 *hh_data_all = data->createHistogram("hh_data_all", x, Binning(8), YVar(y, Binning(8)));
   TH1 *hh_data_sel = data->createHistogram("hh_data_sel", x, Binning(8), YVar(y, Binning(8)), Cut("cut==cut::accept"));
   TH1 *hh_eff = effFunc.createHistogram("hh_eff", x, Binning(50), YVar(y, Binning(50)));

   // Some adjustment for good visualization
   hh_data_all->SetMinimum(0);
   hh_data_sel->SetMinimum(0);
   hh_eff->SetMinimum(0);
   hh_eff->SetLineColor(kBlue);

   // Draw all frames on a canvas
   TCanvas *ca = new TCanvas("rf702_efficiency_2D", "rf702_efficiency_2D", 1200, 400);
   ca->Divide(3);
   ca->cd(1);
   gPad->SetLeftMargin(0.15);
   hh_data_all->GetZaxis()->SetTitleOffset(1.8);
   hh_data_all->Draw("lego");
   ca->cd(2);
   gPad->SetLeftMargin(0.15);
   hh_data_sel->GetZaxis()->SetTitleOffset(1.8);
   hh_data_sel->Draw("lego");
   ca->cd(3);
   gPad->SetLeftMargin(0.15);
   hh_eff->GetZaxis()->SetTitleOffset(1.8);
   hh_eff->Draw("surf");

   return;
}
