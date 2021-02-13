/// \file
/// \ingroup tutorial_roofit
/// \notebook -js
/// Addition and convolution: convolution in cyclical angular observables theta
///
/// and construction of pdf in terms of transformed angular coordinates, e.g. cos(theta),
/// where the convolution is performed in theta rather than cos(theta)
///
/// ```
///   pdf(theta)    = T(theta)          (x) gauss(theta)
///   pdf(cosTheta) = T(acos(cosTheta)) (x) gauss(acos(cosTheta))
/// ```
///
/// This tutorial requires FFT3 to be enabled.
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
#include "RooGenericPdf.h"
#include "RooFormulaVar.h"
#include "RooFFTConvPdf.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
using namespace RooFit;

void rf210_angularconv()
{
   // S e t u p   c o m p o n e n t   p d f s
   // ---------------------------------------

   // Define angle psi
   RooRealVar psi("psi", "psi", 0, 3.14159268);

   // Define physics pdf T(psi)
   RooGenericPdf Tpsi("Tpsi", "1+sin(2*@0)", psi);

   // Define resolution R(psi)
   RooRealVar gbias("gbias", "gbias", 0.2, 0., 1);
   RooRealVar greso("greso", "greso", 0.3, 0.1, 1.0);
   RooGaussian Rpsi("Rpsi", "Rpsi", psi, gbias, greso);

   // Define cos(psi) and function psif that calculates psi from cos(psi)
   RooRealVar cpsi("cpsi", "cos(psi)", -1, 1);
   RooFormulaVar psif("psif", "acos(cpsi)", cpsi);

   // Define physics pdf also as function of cos(psi): T(psif(cpsi)) = T(cpsi) ;
   RooGenericPdf Tcpsi("T", "1+sin(2*@0)", psif);

   // C o n s t r u c t   c o n v o l u t i o n   p d f  i n   p s i
   // --------------------------------------------------------------

   // Define convoluted pdf as function of psi: M=[T(x)R](psi) = M(psi)
   RooFFTConvPdf Mpsi("Mf", "Mf", psi, Tpsi, Rpsi);

   // Set the buffer fraction to zero to obtain a true cyclical convolution
   Mpsi.setBufferFraction(0);

   // S a m p l e ,   f i t   a n d   p l o t   c o n v o l u t e d   p d f  ( p s i )
   // --------------------------------------------------------------------------------

   // Generate some events in observable psi
   RooDataSet *data_psi = Mpsi.generate(psi, 10000);

   // Fit convoluted model as function of angle psi
   Mpsi.fitTo(*data_psi);

   // Plot cos(psi) frame with Mf(cpsi)
   RooPlot *frame1 = psi.frame(Title("Cyclical convolution in angle psi"));
   data_psi->plotOn(frame1);
   Mpsi.plotOn(frame1);

   // Overlay comparison to unsmeared physics pdf T(psi)
   Tpsi.plotOn(frame1, LineColor(kRed));

   // C o n s t r u c t   c o n v o l u t i o n   p d f   i n   c o s ( p s i )
   // --------------------------------------------------------------------------

   // Define convoluted pdf as function of cos(psi): M=[T(x)R](psif(cpsi)) = M(cpsi)
   //
   // Need to give both observable psi here (for definition of convolution)
   // and function psif here (for definition of observables, ultimately in cpsi)
   RooFFTConvPdf Mcpsi("Mf", "Mf", psif, psi, Tpsi, Rpsi);

   // Set the buffer fraction to zero to obtain a true cyclical convolution
   Mcpsi.setBufferFraction(0);

   // S a m p l e ,   f i t   a n d   p l o t   c o n v o l u t e d   p d f  ( c o s p s i )
   // --------------------------------------------------------------------------------

   // Generate some events
   RooDataSet *data_cpsi = Mcpsi.generate(cpsi, 10000);

   // set psi constant to exclude to be a parameter of the fit
   psi.setConstant(true);

   // Fit convoluted model as function of cos(psi)
   Mcpsi.fitTo(*data_cpsi);

   // Plot cos(psi) frame with Mf(cpsi)
   RooPlot *frame2 = cpsi.frame(Title("Same convolution in psi, expressed in cos(psi)"));
   data_cpsi->plotOn(frame2);
   Mcpsi.plotOn(frame2);

   // Overlay comparison to unsmeared physics pdf Tf(cpsi)
   Tcpsi.plotOn(frame2, LineColor(kRed));

   // Draw frame on canvas
   TCanvas *c = new TCanvas("rf210_angularconv", "rf210_angularconv", 800, 400);
   c->Divide(2);
   c->cd(1);
   gPad->SetLeftMargin(0.15);
   frame1->GetYaxis()->SetTitleOffset(1.4);
   frame1->Draw();
   c->cd(2);
   gPad->SetLeftMargin(0.15);
   frame2->GetYaxis()->SetTitleOffset(1.4);
   frame2->Draw();
}
