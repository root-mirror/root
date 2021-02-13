/// \file
/// \ingroup tutorial_roofit
/// \notebook -nodraw
/// Organisation and simultaneous fits: reading and writing ASCII configuration files
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
using namespace RooFit;

void rf505_asciicfg()
{
   // C r e a t e  p d f
   // ------------------

   // Construct gauss(x,m,s)
   RooRealVar x("x", "x", -10, 10);
   RooRealVar m("m", "m", 0, -10, 10);
   RooRealVar s("s", "s", 1, -10, 10);
   RooGaussian gauss("g", "g", x, m, s);

   // Construct poly(x,p0)
   RooRealVar p0("p0", "p0", 0.01, 0., 1.);
   RooPolynomial poly("p", "p", x, p0);

   // Construct model = f*gauss(x) + (1-f)*poly(x)
   RooRealVar f("f", "f", 0.5, 0., 1.);
   RooAddPdf model("model", "model", RooArgSet(gauss, poly), f);

   // F i t   m o d e l   t o   t o y   d a t a
   // -----------------------------------------

   RooDataSet *d = model.generate(x, 1000);
   model.fitTo(*d);

   // W r i t e   p a r a m e t e r s   t o   a s c i i   f i l e
   // -----------------------------------------------------------

   // Obtain set of parameters
   RooArgSet *params = model.getParameters(x);

   // Write parameters to file
   params->writeToFile("rf505_asciicfg_example.txt");

   TString dir1 = gROOT->GetTutorialDir() ;
   dir1.Append("/roofit/rf505_asciicfg.txt") ;
   TString dir2 = "rf505_asciicfg_example.txt";

   // R e a d    p a r a m e t e r s   f r o m    a s c i i   f i l e
   // ----------------------------------------------------------------

   // Read parameters from file
   params->readFromFile(dir2);
   params->Print("v");

   // Read parameters from section 'Section2' of file
   params->readFromFile(dir1, 0, "Section2");
   params->Print("v");

   // Read parameters from section 'Section3' of file. Mark all
   // variables that were processed with the "READ" attribute
   params->readFromFile(dir1, "READ", "Section3");

   // Print the list of parameters that were not read from Section3
   cout << "The following parameters of the were _not_ read from Section3: "
        << (*params->selectByAttrib("READ", kFALSE)) << endl;

   // Read parameters from section 'Section4' of file, which contains
   // 'include file' statement of rf505_asciicfg_example.txt
   // so that we effective read the same

   params->readFromFile(dir1, 0, "Section4");
   params->Print("v");
}
