/// \file
/// \ingroup tutorial_roofit
/// \notebook -js
/// Basic functionality: fitting, plotting, toy data generation on one-dimensional PDFs.
///
///  pdf = gauss(x,m,s)
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
#include "TCanvas.h"
#include "RooPlot.h"
#include "TAxis.h"
using namespace RooFit;

void rf101_basics()
{
   // S e t u p   m o d e l
   // ---------------------

   // Declare variables x,mean,sigma with associated name, title, initial value and allowed range
   RooRealVar x("x", "x", -10, 10);
   RooRealVar mean("mean", "mean of gaussian", 1, -10, 10);
   RooRealVar sigma("sigma", "width of gaussian", 1, 0.1, 10);

   // Build gaussian pdf in terms of x,mean and sigma
   RooGaussian gauss("gauss", "gaussian PDF", x, mean, sigma);

   // Construct plot frame in 'x'
   RooPlot *xframe = x.frame(Title("Gaussian pdf."));

   // P l o t   m o d e l   a n d   c h a n g e   p a r a m e t e r   v a l u e s
   // ---------------------------------------------------------------------------

   // Plot gauss in frame (i.e. in x)
   gauss.plotOn(xframe);

   // Change the value of sigma to 3
   sigma.setVal(3);

   // Plot gauss in frame (i.e. in x) and draw frame on canvas
   gauss.plotOn(xframe, LineColor(kRed));

   // G e n e r a t e   e v e n t s
   // -----------------------------

   // Generate a dataset of 1000 events in x from gauss
   RooDataSet *data = gauss.generate(x, 10000);

   // Make a second plot frame in x and draw both the
   // data and the pdf in the frame
   RooPlot *xframe2 = x.frame(Title("Gaussian pdf with data"));
   data->plotOn(xframe2);
   gauss.plotOn(xframe2);

   // F i t   m o d e l   t o   d a t a
   // -----------------------------

   // Fit pdf to data
   gauss.fitTo(*data);

   // Print values of mean and sigma (that now reflect fitted values and errors)
   mean.Print();
   sigma.Print();

   // Draw all frames on a canvas
   TCanvas *c = new TCanvas("rf101_basics", "rf101_basics", 800, 400);
   c->Divide(2);
   c->cd(1);
   gPad->SetLeftMargin(0.15);
   xframe->GetYaxis()->SetTitleOffset(1.6);
   xframe->Draw();
   c->cd(2);
   gPad->SetLeftMargin(0.15);
   xframe2->GetYaxis()->SetTitleOffset(1.6);
   xframe2->Draw();
}
