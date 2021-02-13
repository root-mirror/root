/// \file
/// \ingroup tutorial_roofit
/// \notebook
/// Multidimensional models: conditional pdf with per-event errors
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
#include "RooGaussModel.h"
#include "RooConstVar.h"
#include "RooDecay.h"
#include "RooLandau.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
using namespace RooFit;

void rf306_condpereventerrors()
{
   // B - p h y s i c s   p d f   w i t h   p e r - e v e n t  G a u s s i a n   r e s o l u t i o n
   // ----------------------------------------------------------------------------------------------

   // Observables
   RooRealVar dt("dt", "dt", -10, 10);
   RooRealVar dterr("dterr", "per-event error on dt", 0.01, 10);

   // Build a gaussian resolution model scaled by the per-event error = gauss(dt,bias,sigma*dterr)
   RooRealVar bias("bias", "bias", 0, -10, 10);
   RooRealVar sigma("sigma", "per-event error scale factor", 1, 0.1, 10);
   RooGaussModel gm("gm1", "gauss model scaled bt per-event error", dt, bias, sigma, dterr);

   // Construct decay(dt) (x) gauss1(dt|dterr)
   RooRealVar tau("tau", "tau", 1.548);
   RooDecay decay_gm("decay_gm", "decay", dt, tau, gm, RooDecay::DoubleSided);

   // C o n s t r u c t   f a k e   ' e x t e r n a l '   d a t a    w i t h   p e r - e v e n t   e r r o r
   // ------------------------------------------------------------------------------------------------------

   // Use landau pdf to get somewhat realistic distribution with long tail
   RooLandau pdfDtErr("pdfDtErr", "pdfDtErr", dterr, RooConst(1), RooConst(0.25));
   RooDataSet *expDataDterr = pdfDtErr.generate(dterr, 10000);

   // S a m p l e   d a t a   f r o m   c o n d i t i o n a l   d e c a y _ g m ( d t | d t e r r )
   // ---------------------------------------------------------------------------------------------

   // Specify external dataset with dterr values to use decay_dm as conditional pdf
   RooDataSet *data = decay_gm.generate(dt, ProtoData(*expDataDterr));

   // F i t   c o n d i t i o n a l   d e c a y _ d m ( d t | d t e r r )
   // ---------------------------------------------------------------------

   // Specify dterr as conditional observable
   decay_gm.fitTo(*data, ConditionalObservables(dterr));

   // P l o t   c o n d i t i o n a l   d e c a y _ d m ( d t | d t e r r )
   // ---------------------------------------------------------------------

   // Make two-dimensional plot of conditional pdf in (dt,dterr)
   TH1 *hh_decay = decay_gm.createHistogram("hh_decay", dt, Binning(50), YVar(dterr, Binning(50)));
   hh_decay->SetLineColor(kBlue);

   // Plot decay_gm(dt|dterr) at various values of dterr
   RooPlot *frame = dt.frame(Title("Slices of decay(dt|dterr) at various dterr"));
   for (Int_t ibin = 0; ibin < 100; ibin += 20) {
      dterr.setBin(ibin);
      decay_gm.plotOn(frame, Normalization(5.));
   }

   // Make projection of data an dt
   RooPlot *frame2 = dt.frame(Title("Projection of decay(dt|dterr) on dt"));
   data->plotOn(frame2);

   // Make projection of decay(dt|dterr) on dt.
   //
   // Instead of integrating out dterr, make a weighted average of curves
   // at values dterr_i as given in the external dataset.
   // (The kTRUE argument bins the data before projection to speed up the process)
   decay_gm.plotOn(frame2, ProjWData(*expDataDterr, kTRUE));

   // Draw all frames on canvas
   TCanvas *c = new TCanvas("rf306_condpereventerrors", "rf306_condperventerrors", 1200, 400);
   c->Divide(3);
   c->cd(1);
   gPad->SetLeftMargin(0.20);
   hh_decay->GetZaxis()->SetTitleOffset(2.5);
   hh_decay->Draw("surf");
   c->cd(2);
   gPad->SetLeftMargin(0.15);
   frame->GetYaxis()->SetTitleOffset(1.6);
   frame->Draw();
   c->cd(3);
   gPad->SetLeftMargin(0.15);
   frame2->GetYaxis()->SetTitleOffset(1.6);
   frame2->Draw();
}
