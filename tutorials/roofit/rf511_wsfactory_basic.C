/// \file
/// \ingroup tutorial_roofit
/// \notebook -nodraw
/// Organization and simultaneous fits: basic use of the 'object factory' associated with
/// a workspace to rapidly build pdfs functions and their parameter components
///
/// \macro_output
/// \macro_code
///
/// \date July 2009
/// \author Wouter Verkerke

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooWorkspace.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
using namespace RooFit;

void rf511_wsfactory_basic(Bool_t compact = kFALSE)
{
   RooWorkspace *w = new RooWorkspace("w");

   // C r e a t i n g   a n d   a d d i n g   b a s i c  p . d . f . s
   // ----------------------------------------------------------------

   // Remake example pdf of tutorial rs502_wspacewrite.C:
   //
   // Basic pdf construction: ClassName::ObjectName(constructor arguments)
   // Variable construction    : VarName[x,xlo,xhi], VarName[xlo,xhi], VarName[x]
   // P.d.f. addition          : SUM::ObjectName(coef1*pdf1,...coefM*pdfM,pdfN)
   //

   if (!compact) {

      // Use object factory to build pdf of tutorial rs502_wspacewrite
      w->factory("Gaussian::sig1(x[-10,10],mean[5,0,10],0.5)");
      w->factory("Gaussian::sig2(x,mean,1)");
      w->factory("Chebychev::bkg(x,{a0[0.5,0.,1],a1[0.2,0.,1.]})");
      w->factory("SUM::sig(sig1frac[0.8,0.,1.]*sig1,sig2)");
      w->factory("SUM::model(bkgfrac[0.5,0.,1.]*bkg,sig)");

   } else {

      // Use object factory to build pdf of tutorial rs502_wspacewrite but
      //  - Contracted to a single line recursive expression,
      //  - Omitting explicit names for components that are not referred to explicitly later

      w->factory("SUM::model(bkgfrac[0.5,0.,1.]*Chebychev::bkg(x[-10,10],{a0[0.5,0.,1],a1[0.2,0.,1.]}),"
                 "SUM(sig1frac[0.8,0.,1.]*Gaussian(x,mean[5,0,10],0.5), Gaussian(x,mean,1)))");
   }

   // A d v a n c e d   p . d . f .  c o n s t r u c t o r   a r g u m e n t s
   // ----------------------------------------------------------------
   //
   // P.d.f. constructor arguments may by any type of RooAbsArg, but also
   //
   // Double_t --> converted to RooConst(...)
   // {a,b,c} --> converted to RooArgSet() or RooArgList() depending on required ctor arg
   // dataset name --> converted to RooAbsData reference for any dataset residing in the workspace
   // enum --> Any enum label that belongs to an enum defined in the (base) class

   // Make a dummy dataset pdf 'model' and import it in the workspace
   RooDataSet *data = w->pdf("model")->generate(*w->var("x"), 1000);
   w->import(*data, Rename("data"));

   // Construct a KEYS pdf passing a dataset name and an enum type defining the
   // mirroring strategy
   w->factory("KeysPdf::k(x,data,NoMirror,0.2)");

   // Print workspace contents
   w->Print();
}
