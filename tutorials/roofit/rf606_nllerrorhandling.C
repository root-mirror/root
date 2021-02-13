/// \file
/// \ingroup tutorial_roofit
/// \notebook -js
/// Likelihood and minimization: understanding and customizing error handling in likelihood evaluations
///
/// \macro_image
/// \macro_output
/// \macro_code
///
/// \date July 2008
/// \author Wouter Verkerke

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooArgusBG.h"
#include "RooNLLVar.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
using namespace RooFit;

void rf606_nllerrorhandling()
{
   // C r e a t e   m o d e l  a n d   d a t a s e t
   // ----------------------------------------------

   // Observable
   RooRealVar m("m", "m", 5.20, 5.30);

   // Parameters
   RooRealVar m0("m0", "m0", 5.291, 5.20, 5.30);
   RooRealVar k("k", "k", -30, -50, -10);

   // Pdf
   RooArgusBG argus("argus", "argus", m, m0, k);

   // Sample 1000 events in m from argus
   RooDataSet *data = argus.generate(m, 1000);

   // P l o t   m o d e l   a n d   d a t a
   // --------------------------------------

   RooPlot *frame1 = m.frame(Bins(40), Title("Argus model and data"));
   data->plotOn(frame1);
   argus.plotOn(frame1);

   // F i t   m o d e l   t o   d a t a
   // ---------------------------------

   // The ARGUS background shape has a sharp kinematic cutoff at m=m0
   // and is prone to evaluation errors if the cutoff parameter m0
   // is floated: when the pdf cutoff value is lower than that in data
   // events with m>m0 will have zero probability

   // Perform unbinned ML fit. Print detailed error messages for up to
   // 10 events per likelihood evaluation. The default error handling strategy
   // is to return a very high value of the likelihood to MINUIT if errors occur,
   // which will force MINUIT to retreat from the problematic area

   argus.fitTo(*data, PrintEvalErrors(10));

   // Perform another fit. In this configuration only the number of errors per
   // likelihood evaluation is shown, if it is greater than zero. The
   // EvalErrorWall(kFALSE) arguments disables the default error handling strategy
   // and will cause the actual (problematic) value of the likelihood to be passed
   // to MINUIT.
   //
   // NB: Use of this option is NOT recommended as default strategy as broken -log(L) values
   // can often be lower than 'good' ones because offending events are removed.
   // This may effectively create a false minimum in problem areas. This is clearly
   // illustrated in the second plot

   m0.setError(0.1);
   argus.fitTo(*data, PrintEvalErrors(0), EvalErrorWall(kFALSE));

   // P l o t   l i k e l i h o o d   a s   f u n c t i o n   o f   m 0
   // ------------------------------------------------------------------

   // Construct likelihood function of model and data
   RooNLLVar nll("nll", "nll", argus, *data);

   // Plot likelihood in m0 in range that includes problematic values
   // In this configuration no messages are printed for likelihood evaluation errors,
   // but if an likelihood value evaluates with error, the corresponding value
   // on the curve will be set to the value given in EvalErrorValue().

   RooPlot *frame2 = m0.frame(Range(5.288, 5.293), Title("-log(L) scan vs m0, problematic regions masked"));
   nll.plotOn(frame2, PrintEvalErrors(-1), ShiftToZero(), EvalErrorValue(nll.getVal() + 10), LineColor(kRed));
   frame2->SetMaximum(15);
   frame2->SetMinimum(0);

   TCanvas *c = new TCanvas("rf606_nllerrorhandling", "rf606_nllerrorhandling", 1200, 400);
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
