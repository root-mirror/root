/// \file
/// \ingroup tutorial_roofit
/// \notebook -js
/// Likelihood and minimization: setting up a multi-core parallelized unbinned maximum likelihood fit
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
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
using namespace RooFit;

void rf603_multicpu()
{

   // C r e a t e   3 D   p d f   a n d   d a t a
   // -------------------------------------------

   // Create observables
   RooRealVar x("x", "x", -5, 5);
   RooRealVar y("y", "y", -5, 5);
   RooRealVar z("z", "z", -5, 5);

   // Create signal pdf gauss(x)*gauss(y)*gauss(z)
   RooGaussian gx("gx", "gx", x, RooConst(0), RooConst(1));
   RooGaussian gy("gy", "gy", y, RooConst(0), RooConst(1));
   RooGaussian gz("gz", "gz", z, RooConst(0), RooConst(1));
   RooProdPdf sig("sig", "sig", RooArgSet(gx, gy, gz));

   // Create background pdf poly(x)*poly(y)*poly(z)
   RooPolynomial px("px", "px", x, RooArgSet(RooConst(-0.1), RooConst(0.004)));
   RooPolynomial py("py", "py", y, RooArgSet(RooConst(0.1), RooConst(-0.004)));
   RooPolynomial pz("pz", "pz", z);
   RooProdPdf bkg("bkg", "bkg", RooArgSet(px, py, pz));

   // Create composite pdf sig+bkg
   RooRealVar fsig("fsig", "signal fraction", 0.1, 0., 1.);
   RooAddPdf model("model", "model", RooArgList(sig, bkg), fsig);

   // Generate large dataset
   RooDataSet *data = model.generate(RooArgSet(x, y, z), 200000);

   // P a r a l l e l   f i t t i n g
   // -------------------------------

   // In parallel mode the likelihood calculation is split in N pieces,
   // that are calculated in parallel and added a posteriori before passing
   // it back to MINUIT.

   // Use four processes and time results both in wall time and CPU time
   model.fitTo(*data, NumCPU(4), Timer(kTRUE));

   // P a r a l l e l   M C   p r o j e c t i o n s
   // ----------------------------------------------

   // Construct signal, total likelihood projection on (y,z) observables and likelihood ratio
   RooAbsPdf *sigyz = sig.createProjection(x);
   RooAbsPdf *totyz = model.createProjection(x);
   RooFormulaVar llratio_func("llratio", "log10(@0)-log10(@1)", RooArgList(*sigyz, *totyz));

   // Calculate likelihood ratio for each event, define subset of events with high signal likelihood
   data->addColumn(llratio_func);
   RooDataSet *dataSel = (RooDataSet *)data->reduce(Cut("llratio>0.7"));

   // Make plot frame and plot data
   RooPlot *frame = x.frame(Title("Projection on X with LLratio(y,z)>0.7"), Bins(40));
   dataSel->plotOn(frame);

   // Perform parallel projection using MC integration of pdf using given input dataSet.
   // In this mode the data-weighted average of the pdf is calculated by splitting the
   // input dataset in N equal pieces and calculating in parallel the weighted average
   // one each subset. The N results of those calculations are then weighted into the
   // final result

   // Use four processes
   model.plotOn(frame, ProjWData(*dataSel), NumCPU(4));

   new TCanvas("rf603_multicpu", "rf603_multicpu", 600, 600);
   gPad->SetLeftMargin(0.15);
   frame->GetYaxis()->SetTitleOffset(1.6);
   frame->Draw();
}
