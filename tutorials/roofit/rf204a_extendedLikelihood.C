/// \file
/// \ingroup tutorial_roofit
/// \notebook -js
/// Extended maximum likelihood fit in multiple ranges.
///
///  When an extended pdf and multiple ranges are used, the
///  RooExtendPdf cannot correctly interpret the coefficients
///  used for extension.
///  This can be solved by using a RooAddPdf for extending the model.
///
/// \macro_output
/// \macro_code
///
/// \date December 2018
/// \author Stephan Hageboeck


#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooFitResult.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
using namespace RooFit ;


void rf204a_extendedLikelihood()
{


   // S e t u p   c o m p o n e n t   p d f s
   // ---------------------------------------

   // Declare observable x
   RooRealVar x("x","x",0,11) ;

   // Create two Gaussian PDFs g1(x,mean1,sigma) anf g2(x,mean2,sigma) and their parameters
   RooRealVar mean("mean","mean of gaussians",5) ;
   RooRealVar sigma1("sigma1","width of gaussians",0.5) ;
   RooRealVar sigma2("sigma2","width of gaussians",1) ;

   RooGaussian sig1("sig1","Signal component 1",x,mean,sigma1) ;
   RooGaussian sig2("sig2","Signal component 2",x,mean,sigma2) ;

   // Build Chebychev polynomial pdf
   RooRealVar a0("a0","a0",0.5,0.,1.) ;
   RooRealVar a1("a1","a1",0.2,0.,1.) ;
   RooChebychev bkg("bkg","Background",x,RooArgSet(a0,a1)) ;

   // Sum the signal components into a composite signal pdf
   RooRealVar sig1frac("sig1frac","fraction of component 1 in signal",0.8,0.,1.) ;
   RooAddPdf sig("sig","Signal",RooArgList(sig1,sig2),sig1frac) ;


   // E x t e n d   t h e   p d f s
   // -----------------------------


   // Define signal range in which events counts are to be defined
   x.setRange("signalRange",4,6) ;

   // Associated nsig/nbkg as expected number of events with sig/bkg _in_the_range_ "signalRange"
   RooRealVar nsig("nsig","number of signal events in signalRange",500,0.,10000) ;
   RooRealVar nbkg("nbkg","number of background events in signalRange",500,0,10000) ;

   // Use AddPdf to extend the model. Giving as many coefficients as pdfs switches
   // on extension.
   RooAddPdf  model("model","(g1+g2)+a", RooArgList(bkg,sig), RooArgList(nbkg,nsig)) ;


   // S a m p l e   d a t a ,   f i t   m o d e l
   // -------------------------------------------

   // Generate 1000 events from model so that nsig,nbkg come out to numbers <<500 in fit
   RooDataSet *data = model.generate(x,1000) ;



   auto canv = new TCanvas("Canvas", "Canvas", 1500, 600);
   canv->Divide(3,1);

   // Fit full range
   // -------------------------------------------

   canv->cd(1);

   // Perform unbinned ML fit to data, full range

   // IMPORTANT:
   // The model needs to be copied when fitting with different ranges because
   // the interpretation of the coefficients is tied to the fit range
   // that's used in the first fit
   RooAddPdf model1(model);
   RooFitResult* r = model1.fitTo(*data,Save()) ;
   r->Print() ;

   RooPlot * frame = x.frame(Title("Full range fitted"));
   data->plotOn(frame);
   model1.plotOn(frame, VisualizeError(*r));
   model1.plotOn(frame);
   model1.paramOn(frame);
   frame->Draw();


   // Fit in two regions
   // -------------------------------------------

   canv->cd(2);
   x.setRange("left",  0., 4.);
   x.setRange("right", 6., 10.);

   RooAddPdf model2(model);
   RooFitResult* r2 = model2.fitTo(*data,
      Range("left,right"),
      Save()) ;
   r2->Print();


   RooPlot * frame2 = x.frame(Title("Fit in left/right sideband"));
   data->plotOn(frame2);
   model2.plotOn(frame2, VisualizeError(*r2));
   model2.plotOn(frame2);
   model2.paramOn(frame2);
   frame2->Draw();


   // Fit in one region
   // -------------------------------------------
   // Note how restricting the region to only the left tail increases
   // the fit uncertainty

   canv->cd(3);
   x.setRange("leftToMiddle",  0., 5.);

   RooAddPdf model3(model);
   RooFitResult* r3 = model3.fitTo(*data,
      Range("leftToMiddle"),
      Save()) ;
   r3->Print();


   RooPlot * frame3 = x.frame(Title("Fit from left to middle"));
   data->plotOn(frame3);
   model3.plotOn(frame3, VisualizeError(*r3));
   model3.plotOn(frame3);
   model3.paramOn(frame3);
   frame3->Draw();

   canv->Draw();
}
