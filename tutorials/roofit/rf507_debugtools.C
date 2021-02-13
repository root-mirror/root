/// \file
/// \ingroup tutorial_roofit
/// \notebook -nodraw
/// Organization and simultaneous fits: RooFit memory tracing debug tool
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
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooTrace.h"
using namespace RooFit;

void rf507_debugtools()
{
   // Activate RooFit memory tracing
   RooTrace::active(kTRUE);

   // Construct gauss(x,m,s)
   RooRealVar x("x", "x", -10, 10);
   RooRealVar m("m", "m", 0, -10, 10);
   RooRealVar s("s", "s", 1, -10, 10);
   RooGaussian gauss("g", "g", x, m, s);

   // Show dump of all RooFit object in memory
   RooTrace::dump();

   // Activate verbose mode
   RooTrace::verbose(kTRUE);

   // Construct poly(x,p0)
   RooRealVar p0("p0", "p0", 0.01, 0., 1.);
   RooPolynomial poly("p", "p", x, p0);

   // Put marker in trace list for future reference
   RooTrace::mark();

   // Construct model = f*gauss(x) + (1-f)*poly(x)
   RooRealVar f("f", "f", 0.5, 0., 1.);
   RooAddPdf model("model", "model", RooArgSet(gauss, poly), f);

   // Show object added to memory since marker
   RooTrace::printObjectCounts();

   // Since verbose mode is still on, you will see messages
   // pertaining to destructor calls of all RooFit objects
   // made in this macro
   //
   // A call to RooTrace::dump() at the end of this macro
   // should show that there a no RooFit object left in memory
}
