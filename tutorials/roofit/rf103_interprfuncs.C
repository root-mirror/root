/// \file
/// \ingroup tutorial_roofit
/// \notebook -js
/// Basic functionality: interpreted functions and PDFs.
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
#include "TAxis.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooGenericPdf.h"
#include "RooConstVar.h"

using namespace RooFit;

void rf103_interprfuncs()
{
   // ----------------------------------------------------
   // G e n e r i c   i n t e r p r e t e d   p . d . f .
   // ====================================================

   // Declare observable x
   RooRealVar x("x", "x", -20, 20);

   // C o n s t r u c t   g e n e r i c   p d f   f r o m   i n t e r p r e t e d   e x p r e s s i o n
   // -------------------------------------------------------------------------------------------------

   // To construct a proper pdf, the formula expression is explicitly normalized internally by dividing
   // it by a numeric integral of the expression over x in the range [-20,20]
   //
   RooRealVar alpha("alpha", "alpha", 5, 0.1, 10);
   RooGenericPdf genpdf("genpdf", "genpdf", "(1+0.1*abs(x)+sin(sqrt(abs(x*alpha+0.1))))", RooArgSet(x, alpha));

   // S a m p l e ,   f i t   a n d   p l o t   g e n e r i c   p d f
   // ---------------------------------------------------------------

   // Generate a toy dataset from the interpreted pdf
   RooDataSet *data = genpdf.generate(x, 10000);

   // Fit the interpreted pdf to the generated data
   genpdf.fitTo(*data);

   // Make a plot of the data and the pdf overlaid
   RooPlot *xframe = x.frame(Title("Interpreted expression pdf"));
   data->plotOn(xframe);
   genpdf.plotOn(xframe);

   // -----------------------------------------------------------------------------------------------------------
   // S t a n d a r d   p . d . f   a d j u s t   w i t h   i n t e r p r e t e d   h e l p e r   f u n c t i o n
   // ==========================================================================================================
   // Make a gauss(x,sqrt(mean2),sigma) from a standard RooGaussian

   // C o n s t r u c t   s t a n d a r d   p d f  w i t h   f o r m u l a   r e p l a c i n g   p a r a m e t e r
   // ------------------------------------------------------------------------------------------------------------

   // Construct parameter mean2 and sigma
   RooRealVar mean2("mean2", "mean^2", 10, 0, 200);
   RooRealVar sigma("sigma", "sigma", 3, 0.1, 10);

   // Construct interpreted function mean = sqrt(mean^2)
   RooFormulaVar mean("mean", "mean", "sqrt(mean2)", mean2);

   // Construct a gaussian g2(x,sqrt(mean2),sigma) ;
   RooGaussian g2("g2", "h2", x, mean, sigma);

   // G e n e r a t e   t o y   d a t a
   // ---------------------------------

   // Construct a separate gaussian g1(x,10,3) to generate a toy Gaussian dataset with mean 10 and width 3
   RooGaussian g1("g1", "g1", x, RooConst(10), RooConst(3));
   RooDataSet *data2 = g1.generate(x, 1000);

   // F i t   a n d   p l o t   t a i l o r e d   s t a n d a r d   p d f
   // -------------------------------------------------------------------

   // Fit g2 to data from g1
   RooFitResult *r = g2.fitTo(*data2, Save());
   r->Print();

   // Plot data on frame and overlay projection of g2
   RooPlot *xframe2 = x.frame(Title("Tailored Gaussian pdf"));
   data2->plotOn(xframe2);
   g2.plotOn(xframe2);

   // Draw all frames on a canvas
   TCanvas *c = new TCanvas("rf103_interprfuncs", "rf103_interprfuncs", 800, 400);
   c->Divide(2);
   c->cd(1);
   gPad->SetLeftMargin(0.15);
   xframe->GetYaxis()->SetTitleOffset(1.4);
   xframe->Draw();
   c->cd(2);
   gPad->SetLeftMargin(0.15);
   xframe2->GetYaxis()->SetTitleOffset(1.4);
   xframe2->Draw();
}
