/// \file
/// \ingroup tutorial_roofit
/// \notebook -js
/// Validation and MC studies: using RooMCStudy on models with constrains
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
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooMCStudy.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
using namespace RooFit;

void rf804_mcstudy_constr()
{
   // C r e a t e   m o d e l   w i t h   p a r a m e t e r   c o n s t r a i n t
   // ---------------------------------------------------------------------------

   // Observable
   RooRealVar x("x", "x", -10, 10);

   // Signal component
   RooRealVar m("m", "m", 0, -10, 10);
   RooRealVar s("s", "s", 2, 0.1, 10);
   RooGaussian g("g", "g", x, m, s);

   // Background component
   RooPolynomial p("p", "p", x);

   // Composite model
   RooRealVar f("f", "f", 0.4, 0., 1.);
   RooAddPdf sum("sum", "sum", RooArgSet(g, p), f);

   // Construct constraint on parameter f
   RooGaussian fconstraint("fconstraint", "fconstraint", f, RooConst(0.7), RooConst(0.1));

   // Multiply constraint with pdf
   RooProdPdf sumc("sumc", "sum with constraint", RooArgSet(sum, fconstraint));

   // S e t u p   t o y   s t u d y   w i t h   m o d e l
   // ---------------------------------------------------

   // Perform toy study with internal constraint on f
   RooMCStudy mcs(sumc, x, Constrain(f), Silence(), Binned(), FitOptions(PrintLevel(-1)));

   // Run 500 toys of 2000 events.
   // Before each toy is generated, a value for the f is sampled from the constraint pdf and
   // that value is used for the generation of that toy.
   mcs.generateAndFit(500, 2000);

   // Make plot of distribution of generated value of f parameter
   TH1 *h_f_gen = mcs.fitParDataSet().createHistogram("f_gen", -40);

   // Make plot of distribution of fitted value of f parameter
   RooPlot *frame1 = mcs.plotParam(f, Bins(40));
   frame1->SetTitle("Distribution of fitted f values");

   // Make plot of pull distribution on f
   RooPlot *frame2 = mcs.plotPull(f, Bins(40), FitGauss());
   frame1->SetTitle("Distribution of f pull values");

   TCanvas *c = new TCanvas("rf804_mcstudy_constr", "rf804_mcstudy_constr", 1200, 400);
   c->Divide(3);
   c->cd(1);
   gPad->SetLeftMargin(0.15);
   h_f_gen->GetYaxis()->SetTitleOffset(1.4);
   h_f_gen->Draw();
   c->cd(2);
   gPad->SetLeftMargin(0.15);
   frame1->GetYaxis()->SetTitleOffset(1.4);
   frame1->Draw();
   c->cd(3);
   gPad->SetLeftMargin(0.15);
   frame2->GetYaxis()->SetTitleOffset(1.4);
   frame2->Draw();
}
