/// \file
/// \ingroup tutorial_roofit
/// \notebook -js
/// Basic functionality: adding boxes with parameters, statistics to RooPlots, decorating with arrows, text etc...
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
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
using namespace RooFit;

void rf106_plotdecoration()
{
   // S e t u p   m o d e l
   // ---------------------

   // Create observables
   RooRealVar x("x", "x", -10, 10);

   // Create Gaussian
   RooRealVar sigma("sigma", "sigma", 1, 0.1, 10);
   RooRealVar mean("mean", "mean", -3, -10, 10);
   RooGaussian gauss("gauss", "gauss", x, mean, sigma);

   // Generate a sample of 1000 events with sigma=3
   RooDataSet *data = gauss.generate(x, 1000);

   // Fit pdf to data
   gauss.fitTo(*data);

   // P l o t   p . d . f   a n d   d a t a
   // -------------------------------------

   // Overlay projection of gauss on data
   RooPlot *frame = x.frame(Name("xframe"), Title("RooPlot with decorations"), Bins(40));
   data->plotOn(frame);
   gauss.plotOn(frame);

   // A d d   b o x   w i t h   p d f   p a r a m e t e r s
   // -----------------------------------------------------

   // Left edge of box starts at 55% of Xaxis)
   gauss.paramOn(frame, Layout(0.55));

   // A d d   b o x   w i t h   d a t a   s t a t i s t i c s
   // -------------------------------------------------------

   // X size of box is from 55% to 99% of Xaxis range, top of box is at 80% of Yaxis range)
   data->statOn(frame, Layout(0.55, 0.99, 0.8));

   // A d d   t e x t   a n d   a r r o w
   // -----------------------------------

   // Add text to frame
   TText *txt = new TText(2, 100, "Signal");
   txt->SetTextSize(0.04);
   txt->SetTextColor(kRed);
   frame->addObject(txt);

   // Add arrow to frame
   TArrow *arrow = new TArrow(2, 100, -1, 50, 0.01, "|>");
   arrow->SetLineColor(kRed);
   arrow->SetFillColor(kRed);
   arrow->SetLineWidth(3);
   frame->addObject(arrow);

   // P e r s i s t   f r a m e   w i t h   a l l   d e c o r a t i o n s   i n   R O O T   f i l e
   // ---------------------------------------------------------------------------------------------

   TFile f("rf106_plotdecoration.root", "RECREATE");
   frame->Write();
   f.Close();

   // To read back and plot frame with all decorations in clean root session do
   // root> TFile f("rf106_plotdecoration.root") ;
   // root>  xframe->Draw() ;

   new TCanvas("rf106_plotdecoration", "rf106_plotdecoration", 600, 600);
   gPad->SetLeftMargin(0.15);
   frame->GetYaxis()->SetTitleOffset(1.6);
   frame->Draw();
}
