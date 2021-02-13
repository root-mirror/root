/// \file
/// \ingroup tutorial_roofit
/// \notebook
/// Multidimensional models: full pdf with per-event errors
///
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
#include "RooProdPdf.h"
#include "RooHistPdf.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
using namespace RooFit;

void rf307_fullpereventerrors()
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

   // C o n s t r u c t   e m p i r i c a l   p d f   f o r   p e r - e v e n t   e r r o r
   // -----------------------------------------------------------------

   // Use landau pdf to get empirical distribution with long tail
   RooLandau pdfDtErr("pdfDtErr", "pdfDtErr", dterr, RooConst(1), RooConst(0.25));
   RooDataSet *expDataDterr = pdfDtErr.generate(dterr, 10000);

   // Construct a histogram pdf to describe the shape of the dtErr distribution
   RooDataHist *expHistDterr = expDataDterr->binnedClone();
   RooHistPdf pdfErr("pdfErr", "pdfErr", dterr, *expHistDterr);

   // C o n s t r u c t   c o n d i t i o n a l   p r o d u c t   d e c a y _ d m ( d t | d t e r r ) * p d f ( d t e r
   // r )
   // ----------------------------------------------------------------------------------------------------------------------

   // Construct production of conditional decay_dm(dt|dterr) with empirical pdfErr(dterr)
   RooProdPdf model("model", "model", pdfErr, Conditional(decay_gm, dt));

   // (Alternatively you could also use the landau shape pdfDtErr)
   // RooProdPdf model("model","model",pdfDtErr,Conditional(decay_gm,dt)) ;

   // S a m p l e,   f i t   a n d   p l o t   p r o d u c t   m o d e l
   // ------------------------------------------------------------------

   // Specify external dataset with dterr values to use model_dm as conditional pdf
   RooDataSet *data = model.generate(RooArgSet(dt, dterr), 10000);

   // F i t   c o n d i t i o n a l   d e c a y _ d m ( d t | d t e r r )
   // ---------------------------------------------------------------------

   // Specify dterr as conditional observable
   model.fitTo(*data);

   // P l o t   c o n d i t i o n a l   d e c a y _ d m ( d t | d t e r r )
   // ---------------------------------------------------------------------

   // Make two-dimensional plot of conditional pdf in (dt,dterr)
   TH1 *hh_model = model.createHistogram("hh_model", dt, Binning(50), YVar(dterr, Binning(50)));
   hh_model->SetLineColor(kBlue);

   // Make projection of data an dt
   RooPlot *frame = dt.frame(Title("Projection of model(dt|dterr) on dt"));
   data->plotOn(frame);
   model.plotOn(frame);

   // Draw all frames on canvas
   TCanvas *c = new TCanvas("rf307_fullpereventerrors", "rf307_fullperventerrors", 800, 400);
   c->Divide(2);
   c->cd(1);
   gPad->SetLeftMargin(0.20);
   hh_model->GetZaxis()->SetTitleOffset(2.5);
   hh_model->Draw("surf");
   c->cd(2);
   gPad->SetLeftMargin(0.15);
   frame->GetYaxis()->SetTitleOffset(1.6);
   frame->Draw();
}
