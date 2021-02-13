/// \file
/// \ingroup tutorial_roofit
/// \notebook -nodraw
/// Numeric algorithm tuning: configuration and customization of how numeric (partial) integrals are executed
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
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooNumIntConfig.h"
#include "RooLandau.h"
#include "RooArgSet.h"
#include <iomanip>
using namespace RooFit;

void rf901_numintconfig()
{

   // A d j u s t   g l o b a l   1 D   i n t e g r a t i o n   p r e c i s i o n
   // ----------------------------------------------------------------------------

   // Print current global default configuration for numeric integration strategies
   RooAbsReal::defaultIntegratorConfig()->Print("v");

   // Example: Change global precision for 1D integrals from 1e-7 to 1e-6
   //
   // The relative epsilon (change as fraction of current best integral estimate) and
   // absolute epsilon (absolute change w.r.t last best integral estimate) can be specified
   // separately. For most pdf integrals the relative change criterium is the most important,
   // however for certain non-pdf functions that integrate out to zero a separate absolute
   // change criterium is necessary to declare convergence of the integral
   //
   // NB: This change is for illustration only. In general the precision should be at least 1e-7
   // for normalization integrals for MINUIT to succeed.
   //
   RooAbsReal::defaultIntegratorConfig()->setEpsAbs(1e-6);
   RooAbsReal::defaultIntegratorConfig()->setEpsRel(1e-6);

   // N u m e r i c   i n t e g r a t i o n   o f   l a n d a u   p d f
   // ------------------------------------------------------------------

   // Construct pdf without support for analytical integrator for demonstration purposes
   RooRealVar x("x", "x", -10, 10);
   RooLandau landau("landau", "landau", x, RooConst(0), RooConst(0.1));

   // Activate debug-level messages for topic integration to be able to follow actions below
   RooMsgService::instance().addStream(DEBUG, Topic(Integration));

   // Calculate integral over landau with default choice of numeric integrator
   RooAbsReal *intLandau = landau.createIntegral(x);
   Double_t val = intLandau->getVal();
   cout << " [1] int_dx landau(x) = " << setprecision(15) << val << endl;

   // S a m e   w i t h   c u s t o m   c o n f i g u r a t i o n
   // -----------------------------------------------------------

   // Construct a custom configuration which uses the adaptive Gauss-Kronrod technique
   // for closed 1D integrals
   RooNumIntConfig customConfig(*RooAbsReal::defaultIntegratorConfig());
#ifdef R__HAS_MATHMORE
   customConfig.method1D().setLabel("RooAdaptiveGaussKronrodIntegrator1D");
#else
   Warning("rf901_numintconfig","ROOT is built without Mathmore (GSL) support. Cannot use RooAdaptiveGaussKronrodIntegrator1D");
#endif

   // Calculate integral over landau with custom integral specification
   RooAbsReal *intLandau2 = landau.createIntegral(x, NumIntConfig(customConfig));
   Double_t val2 = intLandau2->getVal();
   cout << " [2] int_dx landau(x) = " << val2 << endl;

   // A d j u s t i n g   d e f a u l t   c o n f i g   f o r   a   s p e c i f i c   p d f
   // -------------------------------------------------------------------------------------

   // Another possibility: associate custom numeric integration configuration as default for object 'landau'
   landau.setIntegratorConfig(customConfig);

   // Calculate integral over landau custom numeric integrator specified as object default
   RooAbsReal *intLandau3 = landau.createIntegral(x);
   Double_t val3 = intLandau3->getVal();
   cout << " [3] int_dx landau(x) = " << val3 << endl;

   // Another possibility: Change global default for 1D numeric integration strategy on finite domains
#ifdef R__HAS_MATHMORE
   RooAbsReal::defaultIntegratorConfig()->method1D().setLabel("RooAdaptiveGaussKronrodIntegrator1D");

   // A d j u s t i n g   p a r a m e t e r s   o f   a   s p e c i f i c   t e c h n i q u e
   // ---------------------------------------------------------------------------------------

   // Adjust maximum number of steps of RooIntegrator1D in the global default configuration
   RooAbsReal::defaultIntegratorConfig()->getConfigSection("RooIntegrator1D").setRealValue("maxSteps", 30);

   // Example of how to change the parameters of a numeric integrator
   // (Each config section is a RooArgSet with RooRealVars holding real-valued parameters
   //  and RooCategories holding parameters with a finite set of options)
   customConfig.getConfigSection("RooAdaptiveGaussKronrodIntegrator1D").setRealValue("maxSeg", 50);
   customConfig.getConfigSection("RooAdaptiveGaussKronrodIntegrator1D").setCatLabel("method", "15Points");

   // Example of how to print set of possible values for "method" category
   customConfig.getConfigSection("RooAdaptiveGaussKronrodIntegrator1D").find("method")->Print("v");
#endif

}
