/// \file
/// \ingroup tutorial_roofit
/// \notebook -js
/// Multidimensional models: use of tailored pdf as conditional pdfs.s
///
/// `pdf = gauss(x,f(y),sx | y )` with `f(y) = a0 + a1*y`
///
/// \macro_image
/// \macro_output
/// \macro_code
///
/// \date July 2008
/// \author Wouter Verkerke

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooPolyVar.h"
#include "RooProdPdf.h"
#include "RooPlot.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"

using namespace RooFit;

RooDataSet *makeFakeDataXY();

void rf303_conditional()
{
   // S e t u p   c o m p o s e d   m o d e l   g a u s s ( x , m ( y ) , s )
   // -----------------------------------------------------------------------

   // Create observables
   RooRealVar x("x", "x", -10, 10);
   RooRealVar y("y", "y", -10, 10);

   // Create function f(y) = a0 + a1*y
   RooRealVar a0("a0", "a0", -0.5, -5, 5);
   RooRealVar a1("a1", "a1", -0.5, -1, 1);
   RooPolyVar fy("fy", "fy", y, RooArgSet(a0, a1));

   // Create gauss(x,f(y),s)
   RooRealVar sigma("sigma", "width of gaussian", 0.5, 0.1, 2.0);
   RooGaussian model("model", "Gaussian with shifting mean", x, fy, sigma);

   // Obtain fake external experimental dataset with values for x and y
   RooDataSet *expDataXY = makeFakeDataXY();

   // G e n e r a t e   d a t a   f r o m   c o n d i t i o n a l   p . d . f   m o d e l ( x | y )
   // ---------------------------------------------------------------------------------------------

   // Make subset of experimental data with only y values
   RooDataSet *expDataY = (RooDataSet *)expDataXY->reduce(y);

   // Generate 10000 events in x obtained from _conditional_ model(x|y) with y values taken from experimental data
   RooDataSet *data = model.generate(x, ProtoData(*expDataY));
   data->Print();

   // F i t   c o n d i t i o n a l   p . d . f   m o d e l ( x | y )   t o   d a t a
   // ---------------------------------------------------------------------------------------------

   model.fitTo(*expDataXY, ConditionalObservables(y));

   // P r o j e c t   c o n d i t i o n a l   p . d . f   o n   x   a n d   y   d i m e n s i o n s
   // ---------------------------------------------------------------------------------------------

   // Plot x distribution of data and projection of model on x = 1/Ndata sum(data(y_i)) model(x;y_i)
   RooPlot *xframe = x.frame();
   expDataXY->plotOn(xframe);
   model.plotOn(xframe, ProjWData(*expDataY));

   // Speed up (and approximate) projection by using binned clone of data for projection
   RooAbsData *binnedDataY = expDataY->binnedClone();
   model.plotOn(xframe, ProjWData(*binnedDataY), LineColor(kCyan), LineStyle(kDotted));

   // Show effect of projection with too coarse binning
   ((RooRealVar *)expDataY->get()->find("y"))->setBins(5);
   RooAbsData *binnedDataY2 = expDataY->binnedClone();
   model.plotOn(xframe, ProjWData(*binnedDataY2), LineColor(kRed));

   // Make canvas and draw RooPlots
   new TCanvas("rf303_conditional", "rf303_conditional", 600, 460);
   gPad->SetLeftMargin(0.15);
   xframe->GetYaxis()->SetTitleOffset(1.2);
   xframe->Draw();
}

RooDataSet *makeFakeDataXY()
{
   RooRealVar x("x", "x", -10, 10);
   RooRealVar y("y", "y", -10, 10);
   RooArgSet coord(x, y);

   RooDataSet *d = new RooDataSet("d", "d", RooArgSet(x, y));

   for (int i = 0; i < 10000; i++) {
      Double_t tmpy = gRandom->Gaus(0, 10);
      Double_t tmpx = gRandom->Gaus(0.5 * tmpy, 1);
      if (fabs(tmpy) < 10 && fabs(tmpx) < 10) {
         x = tmpx;
         y = tmpy;
         d->add(coord);
      }
   }

   return d;
}
