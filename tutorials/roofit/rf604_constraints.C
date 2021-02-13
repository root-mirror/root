/// \file
/// \ingroup tutorial_roofit
/// \notebook -nodraw
/// Likelihood and minimization: fitting with constraints
///
/// \macro_output
/// \macro_code
///
/// \date 07/2008
/// \author Wouter Verkerke

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
using namespace RooFit;

void rf604_constraints()
{

   // C r e a t e   m o d e l  a n d   d a t a s e t
   // ----------------------------------------------

   // Construct a Gaussian pdf
   RooRealVar x("x", "x", -10, 10);

   RooRealVar m("m", "m", 0, -10, 10);
   RooRealVar s("s", "s", 2, 0.1, 10);
   RooGaussian gauss("gauss", "gauss(x,m,s)", x, m, s);

   // Construct a flat pdf (polynomial of 0th order)
   RooPolynomial poly("poly", "poly(x)", x);

   // Construct model = f*gauss + (1-f)*poly
   RooRealVar f("f", "f", 0.5, 0., 1.);
   RooAddPdf model("model", "model", RooArgSet(gauss, poly), f);

   // Generate small dataset for use in fitting below
   RooDataSet *d = model.generate(x, 50);

   // C r e a t e   c o n s t r a i n t   p d f
   // -----------------------------------------

   // Construct Gaussian constraint pdf on parameter f at 0.8 with resolution of 0.1
   RooGaussian fconstraint("fconstraint", "fconstraint", f, RooConst(0.8), RooConst(0.2));

   // M E T H O D   1   -   A d d   i n t e r n a l   c o n s t r a i n t   t o   m o d e l
   // -------------------------------------------------------------------------------------

   // Multiply constraint term with regular pdf using RooProdPdf
   // Specify in fitTo() that internal constraints on parameter f should be used

   // Multiply constraint with pdf
   RooProdPdf modelc("modelc", "model with constraint", RooArgSet(model, fconstraint));

   // Fit model (without use of constraint term)
   RooFitResult *r1 = model.fitTo(*d, Save());

   // Fit modelc with constraint term on parameter f
   RooFitResult *r2 = modelc.fitTo(*d, Constrain(f), Save());

   // M E T H O D   2   -     S p e c i f y   e x t e r n a l   c o n s t r a i n t   w h e n   f i t t i n g
   // -------------------------------------------------------------------------------------------------------

   // Construct another Gaussian constraint pdf on parameter f at 0.2 with resolution of 0.1
   RooGaussian fconstext("fconstext", "fconstext", f, RooConst(0.2), RooConst(0.1));

   // Fit with external constraint
   RooFitResult *r3 = model.fitTo(*d, ExternalConstraints(fconstext), Save());

   // Print the fit results
   cout << "fit result without constraint (data generated at f=0.5)" << endl;
   r1->Print("v");
   cout << "fit result with internal constraint (data generated at f=0.5, constraint is f=0.8+/-0.2)" << endl;
   r2->Print("v");
   cout << "fit result with (another) external constraint (data generated at f=0.5, constraint is f=0.2+/-0.1)" << endl;
   r3->Print("v");
}
