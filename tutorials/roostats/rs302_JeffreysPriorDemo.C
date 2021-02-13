/// \file
/// \ingroup tutorial_roostats
/// \notebook -js
/// tutorial demonstrating and validates the RooJeffreysPrior class
///
/// Jeffreys's prior is an 'objective prior' based on formal rules.
/// It is calculated from the Fisher information matrix.
///
/// Read more:
/// http://en.wikipedia.org/wiki/Jeffreys_prior
///
/// The analytic form is not known for most PDFs, but it is for
/// simple cases like the Poisson mean, Gaussian mean, Gaussian sigma.
///
/// This class uses numerical tricks to calculate the Fisher Information Matrix
/// efficiently.  In particular, it takes advantage of a property of the
/// 'Asimov data' as described in
/// Asymptotic formulae for likelihood-based tests of new physics
/// Glen Cowan, Kyle Cranmer, Eilam Gross, Ofer Vitells
/// http://arxiv.org/abs/arXiv:1007.1727
///
/// This Demo has four parts:
///  1. TestJeffreysPriorDemo -- validates Poisson mean case 1/sqrt(mu)
///  2. TestJeffreysGaussMean -- validates Gaussian mean case
///  3. TestJeffreysGaussSigma -- validates Gaussian sigma case 1/sigma
///  4. TestJeffreysGaussMeanAndSigma -- demonstrates 2-d example
///
/// \macro_image
/// \macro_output
/// \macro_code
///
/// \author Kyle Cranmer

#include "RooJeffreysPrior.h"

#include "RooWorkspace.h"
#include "RooDataHist.h"
#include "RooGenericPdf.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooNumIntConfig.h"

#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMatrixDSym.h"

using namespace RooFit;

void rs302_JeffreysPriorDemo()
{
   RooWorkspace w("w");
   w.factory("Uniform::u(x[0,1])");
   w.factory("mu[100,1,200]");
   w.factory("ExtendPdf::p(u,mu)");

   RooDataHist *asimov = w.pdf("p")->generateBinned(*w.var("x"), ExpectedData());

   RooFitResult *res = w.pdf("p")->fitTo(*asimov, Save(), SumW2Error(kTRUE));

   asimov->Print();
   res->Print();
   TMatrixDSym cov = res->covarianceMatrix();
   cout << "variance = " << (cov.Determinant()) << endl;
   cout << "stdev = " << sqrt(cov.Determinant()) << endl;
   cov.Invert();
   cout << "jeffreys = " << sqrt(cov.Determinant()) << endl;

   w.defineSet("poi", "mu");
   w.defineSet("obs", "x");

   RooJeffreysPrior pi("jeffreys", "jeffreys", *w.pdf("p"), *w.set("poi"), *w.set("obs"));

   RooGenericPdf *test = new RooGenericPdf("Expected", "Expected = 1/#sqrt{#mu}", "1./sqrt(mu)",
      *w.set("poi"));

   TCanvas *c1 = new TCanvas;
   RooPlot *plot = w.var("mu")->frame();
   pi.plotOn(plot);
   test->plotOn(plot, LineColor(kRed), LineStyle(kDashDotted));
   plot->Draw();

   auto legend = plot->BuildLegend();
   legend->DrawClone();
}

//_________________________________________________
void TestJeffreysGaussMean()
{
   RooWorkspace w("w");
   w.factory("Gaussian::g(x[0,-20,20],mu[0,-5.,5],sigma[1,0,10])");
   w.factory("n[10,.1,200]");
   w.factory("ExtendPdf::p(g,n)");
   w.var("sigma")->setConstant();
   w.var("n")->setConstant();

   RooDataHist *asimov = w.pdf("p")->generateBinned(*w.var("x"), ExpectedData());

   RooFitResult *res = w.pdf("p")->fitTo(*asimov, Save(), SumW2Error(kTRUE));

   asimov->Print();
   res->Print();
   TMatrixDSym cov = res->covarianceMatrix();
   cout << "variance = " << (cov.Determinant()) << endl;
   cout << "stdev = " << sqrt(cov.Determinant()) << endl;
   cov.Invert();
   cout << "jeffreys = " << sqrt(cov.Determinant()) << endl;

   w.defineSet("poi", "mu");
   w.defineSet("obs", "x");

   RooJeffreysPrior pi("jeffreys", "jeffreys", *w.pdf("p"), *w.set("poi"), *w.set("obs"));

   const RooArgSet *temp = w.set("poi");
   pi.getParameters(*temp)->Print();

   //  return;
   RooGenericPdf *test = new RooGenericPdf("constant", "Expected = constant", "1", *w.set("poi"));

   TCanvas *c1 = new TCanvas;
   RooPlot *plot = w.var("mu")->frame();
   pi.plotOn(plot);
   test->plotOn(plot, LineColor(kRed), LineStyle(kDashDotted));
   plot->Draw();

   auto legend = plot->BuildLegend();
   legend->DrawClone();
}

//_________________________________________________
void TestJeffreysGaussSigma()
{
   // this one is VERY sensitive
   // if the Gaussian is narrow ~ range(x)/nbins(x) then the peak isn't resolved
   //   and you get really bizarre shapes
   // if the Gaussian is too wide range(x) ~ sigma then PDF gets renormalized
   //   and the PDF falls off too fast at high sigma
   RooWorkspace w("w");
   w.factory("Gaussian::g(x[0,-20,20],mu[0,-5,5],sigma[1,1,5])");
   w.factory("n[100,.1,2000]");
   w.factory("ExtendPdf::p(g,n)");
   //  w.var("sigma")->setConstant();
   w.var("mu")->setConstant();
   w.var("n")->setConstant();
   w.var("x")->setBins(301);

   RooDataHist *asimov = w.pdf("p")->generateBinned(*w.var("x"), ExpectedData());

   RooFitResult *res = w.pdf("p")->fitTo(*asimov, Save(), SumW2Error(kTRUE));

   asimov->Print();
   res->Print();
   TMatrixDSym cov = res->covarianceMatrix();
   cout << "variance = " << (cov.Determinant()) << endl;
   cout << "stdev = " << sqrt(cov.Determinant()) << endl;
   cov.Invert();
   cout << "jeffreys = " << sqrt(cov.Determinant()) << endl;

   w.defineSet("poi", "sigma");
   w.defineSet("obs", "x");

   RooJeffreysPrior pi("jeffreys", "jeffreys", *w.pdf("p"), *w.set("poi"), *w.set("obs"));

   const RooArgSet *temp = w.set("poi");
   pi.getParameters(*temp)->Print();

   RooGenericPdf *test = new RooGenericPdf("test", "Expected = #sqrt{2}/#sigma", "sqrt(2.)/sigma", *w.set("poi"));

   TCanvas *c1 = new TCanvas;
   RooPlot *plot = w.var("sigma")->frame();
   pi.plotOn(plot);
   test->plotOn(plot, LineColor(kRed), LineStyle(kDashDotted));
   plot->Draw();

   auto legend = plot->BuildLegend();
   legend->DrawClone();
}

//_________________________________________________
void TestJeffreysGaussMeanAndSigma()
{
   // this one is VERY sensitive
   // if the Gaussian is narrow ~ range(x)/nbins(x) then the peak isn't resolved
   //   and you get really bizarre shapes
   // if the Gaussian is too wide range(x) ~ sigma then PDF gets renormalized
   //   and the PDF falls off too fast at high sigma
   RooWorkspace w("w");
   w.factory("Gaussian::g(x[0,-20,20],mu[0,-5,5],sigma[1,1.,5.])");
   w.factory("n[100,.1,2000]");
   w.factory("ExtendPdf::p(g,n)");

   w.var("n")->setConstant();
   w.var("x")->setBins(301);

   RooDataHist *asimov = w.pdf("p")->generateBinned(*w.var("x"), ExpectedData());

   RooFitResult *res = w.pdf("p")->fitTo(*asimov, Save(), SumW2Error(kTRUE));

   asimov->Print();
   res->Print();
   TMatrixDSym cov = res->covarianceMatrix();
   cout << "variance = " << (cov.Determinant()) << endl;
   cout << "stdev = " << sqrt(cov.Determinant()) << endl;
   cov.Invert();
   cout << "jeffreys = " << sqrt(cov.Determinant()) << endl;

   w.defineSet("poi", "mu,sigma");
   w.defineSet("obs", "x");

   RooJeffreysPrior pi("jeffreys", "jeffreys", *w.pdf("p"), *w.set("poi"), *w.set("obs"));

   const RooArgSet *temp = w.set("poi");
   pi.getParameters(*temp)->Print();
   //  return;

   TCanvas *c1 = new TCanvas;
   TH1 *Jeff2d = pi.createHistogram("2dJeffreys", *w.var("mu"), Binning(10, -5., 5), YVar(*w.var("sigma"), Binning(10, 1., 5.)));
   Jeff2d->Draw("surf");
}
