/// \file
/// \ingroup tutorial_roofit
/// \notebook -js
/// Fitting and plotting in sub ranges.
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
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
using namespace RooFit;

void rf203_ranges()
{
   // S e t u p   m o d e l
   // ---------------------

   // Construct observables x
   RooRealVar x("x", "x", -10, 10);

   // Construct gaussx(x,mx,1)
   RooRealVar mx("mx", "mx", 0, -10, 10);
   RooGaussian gx("gx", "gx", x, mx, RooConst(1));

   // Construct px = 1 (flat in x)
   RooPolynomial px("px", "px", x);

   // Construct model = f*gx + (1-f)px
   RooRealVar f("f", "f", 0., 1.);
   RooAddPdf model("model", "model", RooArgList(gx, px), f);

   // Generated 10000 events in (x,y) from pdf model
   RooDataSet *modelData = model.generate(x, 10000);

   // F i t   f u l l   r a n g e
   // ---------------------------

   // Fit pdf to all data
   RooFitResult *r_full = model.fitTo(*modelData, Save(kTRUE));

   // F i t   p a r t i a l   r a n g e
   // ----------------------------------

   // Define "signal" range in x as [-3,3]
   x.setRange("signal", -3, 3);

   // Fit pdf only to data in "signal" range
   RooFitResult *r_sig = model.fitTo(*modelData, Save(kTRUE), Range("signal"));

   // P l o t   /   p r i n t   r e s u l t s
   // ---------------------------------------

   // Make plot frame in x and add data and fitted model
   RooPlot *frame = x.frame(Title("Fitting a sub range"));
   modelData->plotOn(frame);
   model.plotOn(frame, Range("Full"), LineStyle(kDashed), LineColor(kRed)); // Add shape in full ranged dashed
   model.plotOn(frame);                                                     // By default only fitted range is shown

   // Print fit results
   cout << "result of fit on all data " << endl;
   r_full->Print();
   cout << "result of fit in in signal region (note increased error on signal fraction)" << endl;
   r_sig->Print();

   // Draw frame on canvas
   new TCanvas("rf203_ranges", "rf203_ranges", 600, 600);
   gPad->SetLeftMargin(0.15);
   frame->GetYaxis()->SetTitleOffset(1.4);
   frame->Draw();

   return;
}
