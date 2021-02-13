/// \file
/// \ingroup tutorial_roofit
/// \notebook -nodraw
/// Multidimensional models: performing fits in multiple (disjoint) ranges in one or more dimensions
///
/// \macro_output
/// \macro_code
///
/// \date July 2008
/// \author Wouter Verkerke

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooProdPdf.h"
#include "RooAddPdf.h"
#include "RooPolynomial.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooFitResult.h"
using namespace RooFit;

void rf312_multirangefit()
{

   // C r e a t e   2 D   p d f   a n d   d a t a
   // -------------------------------------------

   // Define observables x,y
   RooRealVar x("x", "x", -10, 10);
   RooRealVar y("y", "y", -10, 10);

   // Construct the signal pdf gauss(x)*gauss(y)
   RooRealVar mx("mx", "mx", 1, -10, 10);
   RooRealVar my("my", "my", 1, -10, 10);

   RooGaussian gx("gx", "gx", x, mx, RooConst(1));
   RooGaussian gy("gy", "gy", y, my, RooConst(1));

   RooProdPdf sig("sig", "sig", gx, gy);

   // Construct the background pdf (flat in x,y)
   RooPolynomial px("px", "px", x);
   RooPolynomial py("py", "py", y);
   RooProdPdf bkg("bkg", "bkg", px, py);

   // Construct the composite model sig+bkg
   RooRealVar f("f", "f", 0., 1.);
   RooAddPdf model("model", "model", RooArgList(sig, bkg), f);

   // Sample 10000 events in (x,y) from the model
   RooDataSet *modelData = model.generate(RooArgSet(x, y), 10000);

   // D e f i n e   s i g n a l   a n d   s i d e b a n d   r e g i o n s
   // -------------------------------------------------------------------

   // Construct the SideBand1,SideBand2,Signal regions
   //
   //                    |
   //      +-------------+-----------+
   //      |             |           |
   //      |    Side     |   Sig     |
   //      |    Band1    |   nal     |
   //      |             |           |
   //    --+-------------+-----------+--
   //      |                         |
   //      |           Side          |
   //      |           Band2         |
   //      |                         |
   //      +-------------+-----------+
   //                    |

   x.setRange("SB1", -10, +10);
   y.setRange("SB1", -10, 0);

   x.setRange("SB2", -10, 0);
   y.setRange("SB2", 0, +10);

   x.setRange("SIG", 0, +10);
   y.setRange("SIG", 0, +10);

   x.setRange("FULL", -10, +10);
   y.setRange("FULL", -10, +10);

   // P e r f o r m   f i t s   i n   i n d i v i d u a l   s i d e b a n d   r e g i o n s
   // -------------------------------------------------------------------------------------

   // Perform fit in SideBand1 region (RooAddPdf coefficients will be interpreted in full range)
   RooFitResult *r_sb1 = model.fitTo(*modelData, Range("SB1"), Save());

   // Perform fit in SideBand2 region (RooAddPdf coefficients will be interpreted in full range)
   RooFitResult *r_sb2 = model.fitTo(*modelData, Range("SB2"), Save());

   // P e r f o r m   f i t s   i n   j o i n t    s i d e b a n d   r e g i o n s
   // -----------------------------------------------------------------------------

   // Now perform fit to joint 'L-shaped' sideband region 'SB1|SB2'
   // (RooAddPdf coefficients will be interpreted in full range)
   RooFitResult *r_sb12 = model.fitTo(*modelData, Range("SB1,SB2"), Save());

   // Print results for comparison
   r_sb1->Print();
   r_sb2->Print();
   r_sb12->Print();
}
