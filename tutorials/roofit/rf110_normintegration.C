/// \file
/// \ingroup tutorial_roofit
/// \notebook -js
/// Basic functionality: normalization and integration of pdfs, construction of cumulative distribution
/// monodimensional functions
///
/// \macro_image
/// \macro_output
/// \macro_code
///
/// \date July 2008
/// \author Wouter Verkerke

#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooAbsReal.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
using namespace RooFit;

void rf110_normintegration()
{
   // S e t u p   m o d e l
   // ---------------------

   // Create observables x,y
   RooRealVar x("x", "x", -10, 10);

   // Create pdf gaussx(x,-2,3)
   RooGaussian gx("gx", "gx", x, RooConst(-2), RooConst(3));

   // R e t r i e v e   r a w  &   n o r m a l i z e d   v a l u e s   o f   R o o F i t   p . d . f . s
   // --------------------------------------------------------------------------------------------------

   // Return 'raw' unnormalized value of gx
   cout << "gx = " << gx.getVal() << endl;

   // Return value of gx normalized over x in range [-10,10]
   RooArgSet nset(x);
   cout << "gx_Norm[x] = " << gx.getVal(&nset) << endl;

   // Create object representing integral over gx
   // which is used to calculate  gx_Norm[x] == gx / gx_Int[x]
   RooAbsReal *igx = gx.createIntegral(x);
   cout << "gx_Int[x] = " << igx->getVal() << endl;

   // I n t e g r a t e   n o r m a l i z e d   p d f   o v e r   s u b r a n g e
   // ----------------------------------------------------------------------------

   // Define a range named "signal" in x from -5,5
   x.setRange("signal", -5, 5);

   // Create an integral of gx_Norm[x] over x in range "signal"
   // This is the fraction of of pdf gx_Norm[x] which is in the
   // range named "signal"
   RooAbsReal *igx_sig = gx.createIntegral(x, NormSet(x), Range("signal"));
   cout << "gx_Int[x|signal]_Norm[x] = " << igx_sig->getVal() << endl;

   // C o n s t r u c t   c u m u l a t i v e   d i s t r i b u t i o n   f u n c t i o n   f r o m   p d f
   // -----------------------------------------------------------------------------------------------------

   // Create the cumulative distribution function of gx
   // i.e. calculate Int[-10,x] gx(x') dx'
   RooAbsReal *gx_cdf = gx.createCdf(x);

   // Plot cdf of gx versus x
   RooPlot *frame = x.frame(Title("cdf of Gaussian pdf"));
   gx_cdf->plotOn(frame);

   // Draw plot on canvas
   new TCanvas("rf110_normintegration", "rf110_normintegration", 600, 600);
   gPad->SetLeftMargin(0.15);
   frame->GetYaxis()->SetTitleOffset(1.6);
   frame->Draw();
}
