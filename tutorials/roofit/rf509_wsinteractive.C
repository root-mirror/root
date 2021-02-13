/// \file
/// \ingroup tutorial_roofit
/// \notebook -js
/// Organization and simultaneous fits: easy interactive access to workspace contents - CINT
/// to CLING code migration
///
/// \macro_image
/// \macro_output
/// \macro_code
///
/// \date April 2009
/// \author Wouter Verkerke

using namespace RooFit;

void fillWorkspace(RooWorkspace &w);

void rf509_wsinteractive()
{
   // C r e a t e  a n d   f i l l   w o r k s p a c e
   // ------------------------------------------------

   // Create a workspace named 'w'
   // With CINT w could exports its contents to
   // a same-name C++ namespace in CINT 'namespace w'.
   // but this does not work anymore in CLING.
   // so this tutorial is an example on how to
   // change the code
   RooWorkspace *w1 = new RooWorkspace("w", kTRUE);

   // Fill workspace with pdf and data in a separate function
   fillWorkspace(*w1);

   // Print workspace contents
   w1->Print();

   // this does not work anymore with CLING
   // use normal workspace functionality

   // U s e   w o r k s p a c e   c o n t e n t s
   // ----------------------------------------------

   // Old syntax to use the name space prefix operator to access the workspace contents
   //
   // RooDataSet* d = w::model.generate(w::x,1000) ;
   // RooFitResult* r = w::model.fitTo(*d) ;

   // use normal workspace methods
   RooAbsPdf *model = w1->pdf("model");
   RooRealVar *x = w1->var("x");

   RooDataSet *d = model->generate(*x, 1000);
   RooFitResult *r = model->fitTo(*d);

   // old syntax to access the variable x
   // RooPlot* frame = w::x.frame() ;

   RooPlot *frame = x->frame();
   d->plotOn(frame);

   // OLD syntax to omit x::
   // NB: The 'w::' prefix can be omitted if namespace w is imported in local namespace
   // in the usual C++ way
   //
   // using namespace w;
   // model.plotOn(frame) ;
   // model.plotOn(frame,Components(bkg),LineStyle(kDashed)) ;

   // new correct syntax
   RooAbsPdf *bkg = w1->pdf("bkg");
   model->plotOn(frame);
   model->plotOn(frame, Components(*bkg), LineStyle(kDashed));

   // Draw the frame on the canvas
   new TCanvas("rf509_wsinteractive", "rf509_wsinteractive", 600, 600);
   gPad->SetLeftMargin(0.15);
   frame->GetYaxis()->SetTitleOffset(1.4);
   frame->Draw();
}

void fillWorkspace(RooWorkspace &w)
{
   // C r e a t e  p d f   a n d   f i l l   w o r k s p a c e
   // --------------------------------------------------------

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

   w.import(model);
}
