/// \file
/// \ingroup tutorial_roofit
/// \notebook -nodraw
/// Organisation and simultaneous fits: tuning and customizing the RooFit message logging facility
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
#include "RooMsgService.h"

using namespace RooFit;

void rf506_msgservice()
{
   // C r e a t e   p d f
   // --------------------

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

   RooDataSet *data = model.generate(x, 10);

   // P r i n t   c o n f i g u r a t i o n   o f   m e s s a g e   s e r v i c e
   // ---------------------------------------------------------------------------

   // Print streams configuration
   RooMsgService::instance().Print();
   cout << endl;

   // A d d i n g   I n t e g r a t i o n   t o p i c   t o   e x i s t i n g   I N F O   s t r e a m
   // -----------------------------------------------------------------------------------------------

   // Print streams configuration
   RooMsgService::instance().Print();
   cout << endl;

   // Add Integration topic to existing INFO stream
   RooMsgService::instance().getStream(1).addTopic(Integration);

   // Construct integral over gauss to demonstrate new message stream
   RooAbsReal *igauss = gauss.createIntegral(x);
   igauss->Print();

   // Print streams configuration in verbose, which also shows inactive streams
   cout << endl;
   RooMsgService::instance().Print();
   cout << endl;

   // Remove stream
   RooMsgService::instance().getStream(1).removeTopic(Integration);

   // E x a m p l e s   o f   p d f   v a l u e   t r a c i n g   s t r e a m
   // -----------------------------------------------------------------------

   // Show DEBUG level message on function tracing, trace RooGaussian only
   RooMsgService::instance().addStream(DEBUG, Topic(Tracing), ClassName("RooGaussian"));

   // Perform a fit to generate some tracing messages
   model.fitTo(*data, Verbose(kTRUE));

   // Reset message service to default stream configuration
   RooMsgService::instance().reset();

   // Show DEBUG level message on function tracing on all objects, redirect output to file
   RooMsgService::instance().addStream(DEBUG, Topic(Tracing), OutputFile("rf506_debug.log"));

   // Perform a fit to generate some tracing messages
   model.fitTo(*data, Verbose(kTRUE));

   // Reset message service to default stream configuration
   RooMsgService::instance().reset();

   // E x a m p l e   o f   a n o t h e r   d e b u g g i n g   s t r e a m
   // ---------------------------------------------------------------------

   // Show DEBUG level messages on client/server link state management
   RooMsgService::instance().addStream(DEBUG, Topic(LinkStateMgmt));
   RooMsgService::instance().Print("v");

   // Clone composite pdf g to trigger some link state management activity
   RooAbsArg *gprime = gauss.cloneTree();
   gprime->Print();

   // Reset message service to default stream configuration
   RooMsgService::instance().reset();
}
