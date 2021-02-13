/// \file
/// \ingroup tutorial_roofit
/// \notebook -js
/// Special pdf's: using a product of an (acceptance) efficiency and a pdf as pdf
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
#include "RooExponential.h"
#include "RooEffProd.h"
#include "RooFormulaVar.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
using namespace RooFit;

void rf703_effpdfprod()
{
   // D e f i n e   o b s e r v a b l e s   a n d   d e c a y   p d f
   // ---------------------------------------------------------------

   // Declare observables
   RooRealVar t("t", "t", 0, 5);

   // Make pdf
   RooRealVar tau("tau", "tau", -1.54, -4, -0.1);
   RooExponential model("model", "model", t, tau);

   // D e f i n e   e f f i c i e n c y   f u n c t i o n
   // ---------------------------------------------------

   // Use error function to simulate turn-on slope
   RooFormulaVar eff("eff", "0.5*(TMath::Erf((t-1)/0.5)+1)", t);

   // D e f i n e   d e c a y   p d f   w i t h   e f f i c i e n c y
   // ---------------------------------------------------------------

   // Multiply pdf(t) with efficiency in t
   RooEffProd modelEff("modelEff", "model with efficiency", model, eff);

   // P l o t   e f f i c i e n c y ,   p d f
   // ----------------------------------------

   RooPlot *frame1 = t.frame(Title("Efficiency"));
   eff.plotOn(frame1, LineColor(kRed));

   RooPlot *frame2 = t.frame(Title("Pdf with and without efficiency"));

   model.plotOn(frame2, LineStyle(kDashed));
   modelEff.plotOn(frame2);

   // G e n e r a t e   t o y   d a t a ,   f i t   m o d e l E f f   t o   d a t a
   // ------------------------------------------------------------------------------

   // Generate events. If the input pdf has an internal generator, the internal generator
   // is used and an accept/reject sampling on the efficiency is applied.
   RooDataSet *data = modelEff.generate(t, 10000);

   // Fit pdf. The normalization integral is calculated numerically.
   modelEff.fitTo(*data);

   // Plot generated data and overlay fitted pdf
   RooPlot *frame3 = t.frame(Title("Fitted pdf with efficiency"));
   data->plotOn(frame3);
   modelEff.plotOn(frame3);

   TCanvas *c = new TCanvas("rf703_effpdfprod", "rf703_effpdfprod", 1200, 400);
   c->Divide(3);
   c->cd(1);
   gPad->SetLeftMargin(0.15);
   frame1->GetYaxis()->SetTitleOffset(1.4);
   frame1->Draw();
   c->cd(2);
   gPad->SetLeftMargin(0.15);
   frame2->GetYaxis()->SetTitleOffset(1.6);
   frame2->Draw();
   c->cd(3);
   gPad->SetLeftMargin(0.15);
   frame3->GetYaxis()->SetTitleOffset(1.6);
   frame3->Draw();
}
