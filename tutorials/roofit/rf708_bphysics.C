/// \file
/// \ingroup tutorial_roofit
/// \notebook -js
/// Special pdf's: special decay pdf for B physics with mixing and/or CP violation
///
/// \macro_image
/// \macro_output
/// \macro_code
///
/// \date July 2008
/// \author Wouter Verkerke

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooCategory.h"
#include "RooBMixDecay.h"
#include "RooBCPEffDecay.h"
#include "RooBDecay.h"
#include "RooFormulaVar.h"
#include "RooTruthModel.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
using namespace RooFit;

void rf708_bphysics()
{
   // -------------------------------------
   // B - D e c a y   w i t h   m i x i n g
   // =====================================

   // C o n s t r u c t   p d f
   // -------------------------

   // Observable
   RooRealVar dt("dt", "dt", -10, 10);
   dt.setBins(40);

   // Parameters
   RooRealVar dm("dm", "delta m(B0)", 0.472);
   RooRealVar tau("tau", "tau (B0)", 1.547);
   RooRealVar w("w", "flavour mistag rate", 0.1);
   RooRealVar dw("dw", "delta mistag rate for B0/B0bar", 0.1);

   RooCategory mixState("mixState", "B0/B0bar mixing state");
   mixState.defineType("mixed", -1);
   mixState.defineType("unmixed", 1);

   RooCategory tagFlav("tagFlav", "Flavour of the tagged B0");
   tagFlav.defineType("B0", 1);
   tagFlav.defineType("B0bar", -1);

   // Use delta function resolution model
   RooTruthModel truthModel("tm", "truth model", dt);

   // Construct Bdecay with mixing
   RooBMixDecay bmix("bmix", "decay", dt, mixState, tagFlav, tau, dm, w, dw, truthModel, RooBMixDecay::DoubleSided);

   // P l o t   p d f   i n   v a r i o u s   s l i c e s
   // ---------------------------------------------------

   // Generate some data
   RooDataSet *data = bmix.generate(RooArgSet(dt, mixState, tagFlav), 10000);

   // Plot B0 and B0bar tagged data separately
   // For all plots below B0 and B0 tagged data will look somewhat differently
   // if the flavor tagging mistag rate for B0 and B0 is different (i.e. dw!=0)
   RooPlot *frame1 = dt.frame(Title("B decay distribution with mixing (B0/B0bar)"));

   data->plotOn(frame1, Cut("tagFlav==tagFlav::B0"));
   bmix.plotOn(frame1, Slice(tagFlav, "B0"));

   data->plotOn(frame1, Cut("tagFlav==tagFlav::B0bar"), MarkerColor(kCyan));
   bmix.plotOn(frame1, Slice(tagFlav, "B0bar"), LineColor(kCyan));

   // Plot mixed slice for B0 and B0bar tagged data separately
   RooPlot *frame2 = dt.frame(Title("B decay distribution of mixed events (B0/B0bar)"));

   data->plotOn(frame2, Cut("mixState==mixState::mixed&&tagFlav==tagFlav::B0"));
   bmix.plotOn(frame2, Slice(tagFlav, "B0"), Slice(mixState, "mixed"));

   data->plotOn(frame2, Cut("mixState==mixState::mixed&&tagFlav==tagFlav::B0bar"), MarkerColor(kCyan));
   bmix.plotOn(frame2, Slice(tagFlav, "B0bar"), Slice(mixState, "mixed"), LineColor(kCyan));

   // Plot unmixed slice for B0 and B0bar tagged data separately
   RooPlot *frame3 = dt.frame(Title("B decay distribution of unmixed events (B0/B0bar)"));

   data->plotOn(frame3, Cut("mixState==mixState::unmixed&&tagFlav==tagFlav::B0"));
   bmix.plotOn(frame3, Slice(tagFlav, "B0"), Slice(mixState, "unmixed"));

   data->plotOn(frame3, Cut("mixState==mixState::unmixed&&tagFlav==tagFlav::B0bar"), MarkerColor(kCyan));
   bmix.plotOn(frame3, Slice(tagFlav, "B0bar"), Slice(mixState, "unmixed"), LineColor(kCyan));

   // -------------------------------------------------
   // B - D e c a y   w i t h   C P   v i o l a t i o n
   // =================================================

   // C o n s t r u c t   p d f
   // -------------------------

   // Additional parameters needed for B decay with CPV
   RooRealVar CPeigen("CPeigen", "CP eigen value", -1);
   RooRealVar absLambda("absLambda", "|lambda|", 1, 0, 2);
   RooRealVar argLambda("absLambda", "|lambda|", 0.7, -1, 1);
   RooRealVar effR("effR", "B0/B0bar reco efficiency ratio", 1);

   // Construct Bdecay with CP violation
   RooBCPEffDecay bcp("bcp", "bcp", dt, tagFlav, tau, dm, w, CPeigen, absLambda, argLambda, effR, dw, truthModel,
                      RooBCPEffDecay::DoubleSided);

   // P l o t   s c e n a r i o   1   -   s i n ( 2 b )   =   0 . 7 ,   | l | = 1
   // ---------------------------------------------------------------------------

   // Generate some data
   RooDataSet *data2 = bcp.generate(RooArgSet(dt, tagFlav), 10000);

   // Plot B0 and B0bar tagged data separately
   RooPlot *frame4 = dt.frame(Title("B decay distribution with CPV(|l|=1,Im(l)=0.7) (B0/B0bar)"));

   data2->plotOn(frame4, Cut("tagFlav==tagFlav::B0"));
   bcp.plotOn(frame4, Slice(tagFlav, "B0"));

   data2->plotOn(frame4, Cut("tagFlav==tagFlav::B0bar"), MarkerColor(kCyan));
   bcp.plotOn(frame4, Slice(tagFlav, "B0bar"), LineColor(kCyan));

   // P l o t   s c e n a r i o   2   -   s i n ( 2 b )   =   0 . 7 ,   | l | = 0 . 7
   // -------------------------------------------------------------------------------

   absLambda = 0.7;

   // Generate some data
   RooDataSet *data3 = bcp.generate(RooArgSet(dt, tagFlav), 10000);

   // Plot B0 and B0bar tagged data separately (sin2b = 0.7 plus direct CPV |l|=0.5)
   RooPlot *frame5 = dt.frame(Title("B decay distribution with CPV(|l|=0.7,Im(l)=0.7) (B0/B0bar)"));

   data3->plotOn(frame5, Cut("tagFlav==tagFlav::B0"));
   bcp.plotOn(frame5, Slice(tagFlav, "B0"));

   data3->plotOn(frame5, Cut("tagFlav==tagFlav::B0bar"), MarkerColor(kCyan));
   bcp.plotOn(frame5, Slice(tagFlav, "B0bar"), LineColor(kCyan));

   // ---------------------------------------------------------------------------
   // G e n e r i c   B   d e c a y  w i t h    u s e r   c o e f f i c i e n t s
   // ===========================================================================

   // C o n s t r u c t   p d f
   // -------------------------

   // Model parameters
   RooRealVar DGbG("DGbG", "DGamma/GammaAvg", 0.5, -1, 1);
   RooRealVar Adir("Adir", "-[1-abs(l)**2]/[1+abs(l)**2]", 0);
   RooRealVar Amix("Amix", "2Im(l)/[1+abs(l)**2]", 0.7);
   RooRealVar Adel("Adel", "2Re(l)/[1+abs(l)**2]", 0.7);

   // Derived input parameters for pdf
   RooFormulaVar DG("DG", "Delta Gamma", "@1/@0", RooArgList(tau, DGbG));

   // Construct coefficient functions for sin,cos,sinh modulations of decay distribution
   RooFormulaVar fsin("fsin", "fsin", "@0*@1*(1-2*@2)", RooArgList(Amix, tagFlav, w));
   RooFormulaVar fcos("fcos", "fcos", "@0*@1*(1-2*@2)", RooArgList(Adir, tagFlav, w));
   RooFormulaVar fsinh("fsinh", "fsinh", "@0", RooArgList(Adel));

   // Construct generic B decay pdf using above user coefficients
   RooBDecay bcpg("bcpg", "bcpg", dt, tau, DG, RooConst(1), fsinh, fcos, fsin, dm, truthModel, RooBDecay::DoubleSided);

   // P l o t   -   I m ( l ) = 0 . 7 ,   R e ( l ) = 0 . 7   | l | = 1,   d G / G = 0 . 5
   // -------------------------------------------------------------------------------------

   // Generate some data
   RooDataSet *data4 = bcpg.generate(RooArgSet(dt, tagFlav), 10000);

   // Plot B0 and B0bar tagged data separately
   RooPlot *frame6 = dt.frame(Title("B decay distribution with CPV(Im(l)=0.7,Re(l)=0.7,|l|=1,dG/G=0.5) (B0/B0bar)"));

   data4->plotOn(frame6, Cut("tagFlav==tagFlav::B0"));
   bcpg.plotOn(frame6, Slice(tagFlav, "B0"));

   data4->plotOn(frame6, Cut("tagFlav==tagFlav::B0bar"), MarkerColor(kCyan));
   bcpg.plotOn(frame6, Slice(tagFlav, "B0bar"), LineColor(kCyan));

   TCanvas *c = new TCanvas("rf708_bphysics", "rf708_bphysics", 1200, 800);
   c->Divide(3, 2);
   c->cd(1);
   gPad->SetLeftMargin(0.15);
   frame1->GetYaxis()->SetTitleOffset(1.6);
   frame1->Draw();
   c->cd(2);
   gPad->SetLeftMargin(0.15);
   frame2->GetYaxis()->SetTitleOffset(1.6);
   frame2->Draw();
   c->cd(3);
   gPad->SetLeftMargin(0.15);
   frame3->GetYaxis()->SetTitleOffset(1.6);
   frame3->Draw();
   c->cd(4);
   gPad->SetLeftMargin(0.15);
   frame4->GetYaxis()->SetTitleOffset(1.6);
   frame4->Draw();
   c->cd(5);
   gPad->SetLeftMargin(0.15);
   frame5->GetYaxis()->SetTitleOffset(1.6);
   frame5->Draw();
   c->cd(6);
   gPad->SetLeftMargin(0.15);
   frame6->GetYaxis()->SetTitleOffset(1.6);
   frame6->Draw();
}
