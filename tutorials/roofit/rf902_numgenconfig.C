/// \file
/// \ingroup tutorial_roofit
/// \notebook -nodraw
/// Numeric algorithm tuning: configuration and customization of how MC sampling algorithms on specific pdfs are
/// executed
///
/// \macro_output
/// \macro_code
///
/// \date July 2008
/// \author Wouter Verkerke

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooConstVar.h"
#include "RooChebychev.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooNumGenConfig.h"
#include "RooArgSet.h"
#include <iomanip>
using namespace RooFit;

void rf902_numgenconfig()
{

   // A d j u s t   g l o b a l   MC   s a m p l i n g   s t r a t e g y
   // ------------------------------------------------------------------

   // Example pdf for use below
   RooRealVar x("x", "x", 0, 10);
   RooChebychev model("model", "model", x, RooArgList(RooConst(0), RooConst(0.5), RooConst(-0.1)));

   // Change global strategy for 1D sampling problems without conditional observable
   // (1st kFALSE) and without discrete observable (2nd kFALSE) from RooFoamGenerator,
   // ( an interface to the TFoam MC generator with adaptive subdivisioning strategy ) to RooAcceptReject,
   // a plain accept/reject sampling algorithm [ RooFit default before ROOT 5.23/04 ]
   RooAbsPdf::defaultGeneratorConfig()->method1D(kFALSE, kFALSE).setLabel("RooAcceptReject");

   // Generate 10Kevt using RooAcceptReject
   RooDataSet *data_ar = model.generate(x, 10000, Verbose(kTRUE));
   data_ar->Print();

   // A d j u s t i n g   d e f a u l t   c o n f i g   f o r   a   s p e c i f i c   p d f
   // -------------------------------------------------------------------------------------

   // Another possibility: associate custom MC sampling configuration as default for object 'model'
   // The kTRUE argument will install a clone of the default configuration as specialized configuration
   // for this model if none existed so far
   model.specialGeneratorConfig(kTRUE)->method1D(kFALSE, kFALSE).setLabel("RooFoamGenerator");

   // A d j u s t i n g   p a r a m e t e r s   o f   a   s p e c i f i c   t e c h n i q u e
   // ---------------------------------------------------------------------------------------

   // Adjust maximum number of steps of RooIntegrator1D in the global default configuration
   RooAbsPdf::defaultGeneratorConfig()->getConfigSection("RooAcceptReject").setRealValue("nTrial1D", 2000);

   // Example of how to change the parameters of a numeric integrator
   // (Each config section is a RooArgSet with RooRealVars holding real-valued parameters
   //  and RooCategories holding parameters with a finite set of options)
   model.specialGeneratorConfig()->getConfigSection("RooFoamGenerator").setRealValue("chatLevel", 1);

   // Generate 10Kevt using RooFoamGenerator (FOAM verbosity increased with above chatLevel adjustment for illustration
   // purposes)
   RooDataSet *data_foam = model.generate(x, 10000, Verbose());
   data_foam->Print();
}
