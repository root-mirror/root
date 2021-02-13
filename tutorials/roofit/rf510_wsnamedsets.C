/// \file
/// \ingroup tutorial_roofit
/// \notebook -js
/// Organization and simultaneous fits: working with named parameter sets and parameter
/// snapshots in workspaces
///
/// \macro_image
/// \macro_output
/// \macro_code
///
/// \date April 2009
/// \author Wouter Verkerke

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooWorkspace.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TFile.h"
#include "TH1.h"

using namespace RooFit;

void fillWorkspace(RooWorkspace &w);

void rf510_wsnamedsets()
{
   // C r e a t e   m o d e l   a n d   d a t a s e t
   // -----------------------------------------------

   RooWorkspace *w = new RooWorkspace("w");
   fillWorkspace(*w);

   // Exploit convention encoded in named set "parameters" and "observables"
   // to use workspace contents w/o need for introspected
   RooAbsPdf *model = w->pdf("model");

   // Generate data from pdf in given observables
   RooDataSet *data = model->generate(*w->set("observables"), 1000);

   // Fit model to data
   model->fitTo(*data);

   // Plot fitted model and data on frame of first (only) observable
   RooPlot *frame = ((RooRealVar *)w->set("observables")->first())->frame();
   data->plotOn(frame);
   model->plotOn(frame);

   // Overlay plot with model with reference parameters as stored in snapshots
   w->loadSnapshot("reference_fit");
   model->plotOn(frame, LineColor(kRed));
   w->loadSnapshot("reference_fit_bkgonly");
   model->plotOn(frame, LineColor(kRed), LineStyle(kDashed));

   // Draw the frame on the canvas
   new TCanvas("rf510_wsnamedsets", "rf503_wsnamedsets", 600, 600);
   gPad->SetLeftMargin(0.15);
   frame->GetYaxis()->SetTitleOffset(1.4);
   frame->Draw();

   // Print workspace contents
   w->Print();

   // Workspace will remain in memory after macro finishes
   gDirectory->Add(w);
}

void fillWorkspace(RooWorkspace &w)
{
   // C r e a t e   m o d e l
   // -----------------------

   // Declare observable x
   RooRealVar x("x", "x", 0, 10);

   // Create two Gaussian PDFs g1(x,mean1,sigma) anf g2(x,mean2,sigma) and their parameters
   RooRealVar mean("mean", "mean of gaussians", 5, 0, 10);
   RooRealVar sigma1("sigma1", "width of gaussians", 0.5);
   RooRealVar sigma2("sigma2", "width of gaussians", 1);

   RooGaussian sig1("sig1", "Signal component 1", x, mean, sigma1);
   RooGaussian sig2("sig2", "Signal component 2", x, mean, sigma2);

   // Build Chebychev polynomial pdf
   RooRealVar a0("a0", "a0", 0.5, 0., 1.);
   RooRealVar a1("a1", "a1", 0.2, 0., 1.);
   RooChebychev bkg("bkg", "Background", x, RooArgSet(a0, a1));

   // Sum the signal components into a composite signal pdf
   RooRealVar sig1frac("sig1frac", "fraction of component 1 in signal", 0.8, 0., 1.);
   RooAddPdf sig("sig", "Signal", RooArgList(sig1, sig2), sig1frac);

   // Sum the composite signal and background
   RooRealVar bkgfrac("bkgfrac", "fraction of background", 0.5, 0., 1.);
   RooAddPdf model("model", "g1+g2+a", RooArgList(bkg, sig), bkgfrac);

   // Import model into pdf
   w.import(model);

   // E n c o d e   d e f i n i t i o n   o f   p a r a m e t e r s   i n   w o r k s p a c e
   // ---------------------------------------------------------------------------------------

   // Define named sets "parameters" and "observables", which list which variables should be considered
   // parameters and observables by the users convention
   //
   // Variables appearing in sets _must_ live in the workspace already, or the autoImport flag
   // of defineSet must be set to import them on the fly. Named sets contain only references
   // to the original variables, therefore the value of observables in named sets already
   // reflect their 'current' value
   RooArgSet *params = (RooArgSet *)model.getParameters(x);
   w.defineSet("parameters", *params);
   w.defineSet("observables", x);

   // E n c o d e   r e f e r e n c e   v a l u e   f o r   p a r a m e t e r s   i n   w o r k s p a c e
   // ---------------------------------------------------------------------------------------------------

   // Define a parameter 'snapshot' in the pdf
   // Unlike a named set, a parameter snapshot stores an independent set of values for
   // a given set of variables in the workspace. The values can be stored and reloaded
   // into the workspace variable objects using the loadSnapshot() and saveSnapshot()
   // methods. A snapshot saves the value of each variable, any errors that are stored
   // with it as well as the 'Constant' flag that is used in fits to determine if a
   // parameter is kept fixed or not.

   // Do a dummy fit to a (supposedly) reference dataset here and store the results
   // of that fit into a snapshot
   RooDataSet *refData = model.generate(x, 10000);
   model.fitTo(*refData, PrintLevel(-1));

   // The kTRUE flag imports the values of the objects in (*params) into the workspace
   // If not set, the present values of the workspace parameters objects are stored
   w.saveSnapshot("reference_fit", *params, kTRUE);

   // Make another fit with the signal component forced to zero
   // and save those parameters too

   bkgfrac.setVal(1);
   bkgfrac.setConstant(kTRUE);
   bkgfrac.removeError();
   model.fitTo(*refData, PrintLevel(-1));

   w.saveSnapshot("reference_fit_bkgonly", *params, kTRUE);
}
