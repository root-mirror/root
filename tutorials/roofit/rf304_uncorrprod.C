/// \file
/// \ingroup tutorial_roofit
/// \notebook -js
/// Multidimensional models: simple uncorrelated multi-dimensional pdfs
///
/// `pdf = gauss(x,mx,sx) * gauss(y,my,sy)`
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
#include "RooProdPdf.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
using namespace RooFit;

void rf304_uncorrprod()
{

   // C r e a t e   c o m p o n e n t   p d f s   i n   x   a n d   y
   // ----------------------------------------------------------------

   // Create two pdfs gaussx(x,meanx,sigmax) gaussy(y,meany,sigmay) and its variables
   RooRealVar x("x", "x", -5, 5);
   RooRealVar y("y", "y", -5, 5);

   RooRealVar meanx("mean1", "mean of gaussian x", 2);
   RooRealVar meany("mean2", "mean of gaussian y", -2);
   RooRealVar sigmax("sigmax", "width of gaussian x", 1);
   RooRealVar sigmay("sigmay", "width of gaussian y", 5);

   RooGaussian gaussx("gaussx", "gaussian PDF", x, meanx, sigmax);
   RooGaussian gaussy("gaussy", "gaussian PDF", y, meany, sigmay);

   // C o n s t r u c t   u n c o r r e l a t e d   p r o d u c t   p d f
   // -------------------------------------------------------------------

   // Multiply gaussx and gaussy into a two-dimensional pdf gaussxy
   RooProdPdf gaussxy("gaussxy", "gaussx*gaussy", RooArgList(gaussx, gaussy));

   // S a m p l e   p d f ,   p l o t   p r o j e c t i o n   o n   x   a n d   y
   // ---------------------------------------------------------------------------

   // Generate 10000 events in x and y from gaussxy
   RooDataSet *data = gaussxy.generate(RooArgSet(x, y), 10000);

   // Plot x distribution of data and projection of gaussxy on x = Int(dy) gaussxy(x,y)
   RooPlot *xframe = x.frame(Title("X projection of gauss(x)*gauss(y)"));
   data->plotOn(xframe);
   gaussxy.plotOn(xframe);

   // Plot x distribution of data and projection of gaussxy on y = Int(dx) gaussxy(x,y)
   RooPlot *yframe = y.frame(Title("Y projection of gauss(x)*gauss(y)"));
   data->plotOn(yframe);
   gaussxy.plotOn(yframe);

   // Make canvas and draw RooPlots
   TCanvas *c = new TCanvas("rf304_uncorrprod", "rf304_uncorrprod", 800, 400);
   c->Divide(2);
   c->cd(1);
   gPad->SetLeftMargin(0.15);
   xframe->GetYaxis()->SetTitleOffset(1.4);
   xframe->Draw();
   c->cd(2);
   gPad->SetLeftMargin(0.15);
   yframe->GetYaxis()->SetTitleOffset(1.4);
   yframe->Draw();
}
