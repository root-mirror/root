/// \file
/// \ingroup tutorial_roostats
/// \notebook
/// Bayesian calculator: basic example
///
/// \macro_image
/// \macro_output
/// \macro_code
///
/// \author Gregory Schott

#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooMsgService.h"

#include "RooStats/BayesianCalculator.h"
#include "RooStats/SimpleInterval.h"
#include "TCanvas.h"

using namespace RooFit;
using namespace RooStats;

void rs701_BayesianCalculator(bool useBkg = true, double confLevel = 0.90)
{

   RooWorkspace *w = new RooWorkspace("w", true);
   w->factory("SUM::pdf(s[0.001,15]*Uniform(x[0,1]),b[1,0,2]*Uniform(x))");
   w->factory("Gaussian::prior_b(b,1,1)");
   w->factory("PROD::model(pdf,prior_b)");
   RooAbsPdf *model = w->pdf("model"); // pdf*priorNuisance
   RooArgSet nuisanceParameters(*(w->var("b")));

   RooAbsRealLValue *POI = w->var("s");
   RooAbsPdf *priorPOI = (RooAbsPdf *)w->factory("Uniform::priorPOI(s)");
   RooAbsPdf *priorPOI2 = (RooAbsPdf *)w->factory("GenericPdf::priorPOI2('1/sqrt(@0)',s)");

   w->factory("n[3]"); // observed number of events
   // create a data set with n observed events
   RooDataSet data("data", "", RooArgSet(*(w->var("x")), *(w->var("n"))), "n");
   data.add(RooArgSet(*(w->var("x"))), w->var("n")->getVal());

   // to suppress messages when pdf goes to zero
   RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

   RooArgSet *nuisPar = 0;
   if (useBkg)
      nuisPar = &nuisanceParameters;
   // if (!useBkg) ((RooRealVar *)w->var("b"))->setVal(0);

   double size = 1. - confLevel;
   std::cout << "\nBayesian Result using a Flat prior " << std::endl;
   BayesianCalculator bcalc(data, *model, RooArgSet(*POI), *priorPOI, nuisPar);
   bcalc.SetTestSize(size);
   SimpleInterval *interval = bcalc.GetInterval();
   double cl = bcalc.ConfidenceLevel();
   std::cout << cl << "% CL central interval: [ " << interval->LowerLimit() << " - " << interval->UpperLimit()
             << " ] or " << cl + (1. - cl) / 2 << "% CL limits\n";
   RooPlot *plot = bcalc.GetPosteriorPlot();
   TCanvas *c1 = new TCanvas("c1", "Bayesian Calculator Result");
   c1->Divide(1, 2);
   c1->cd(1);
   plot->Draw();
   c1->Update();

   std::cout << "\nBayesian Result using a 1/sqrt(s) prior  " << std::endl;
   BayesianCalculator bcalc2(data, *model, RooArgSet(*POI), *priorPOI2, nuisPar);
   bcalc2.SetTestSize(size);
   SimpleInterval *interval2 = bcalc2.GetInterval();
   cl = bcalc2.ConfidenceLevel();
   std::cout << cl << "% CL central interval: [ " << interval2->LowerLimit() << " - " << interval2->UpperLimit()
             << " ] or " << cl + (1. - cl) / 2 << "% CL limits\n";

   RooPlot *plot2 = bcalc2.GetPosteriorPlot();
   c1->cd(2);
   plot2->Draw();
   gPad->SetLogy();
   c1->Update();

   // observe one event while expecting one background event -> the 95% CL upper limit on s is 4.10
   // observe one event while expecting zero background event -> the 95% CL upper limit on s is 4.74
}
