/// \file
/// \ingroup tutorial_roofit
/// \notebook -nodraw
/// Organization and simultaneous fits: RooCustomizer and RooSimWSTool interface in factory
/// workspace tool in a complex standalone B physics example
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

void rf513_wsfactory_tools()
{
   RooWorkspace *w = new RooWorkspace("w");

   // B u i l d   a   c o m p l e x   e x a m p l e   p . d . f .
   // -----------------------------------------------------------

   // Make signal model for CPV: A bmixing decay function in t (convoluted with a triple Gaussian resolution model)
   //                            times a Gaussian function the reconstructed mass
   w->factory("PROD::sig(  BMixDecay::sig_t( dt[-20,20], mixState[mixed=1,unmix=-1], tagFlav[B0=1,B0bar=-1], "
              "tau[1.54], dm[0.472], w[0.05], dw[0],"
              "AddModel::gm({GaussModel(dt,biasC[-10,10],sigmaC[0.1,3],dterr[0.01,0.2]),"
              "GaussModel(dt,0,sigmaT[3,10]),"
              "GaussModel(dt,0,20)},{fracC[0,1],fracT[0,1]}),"
              "DoubleSided ),"
              "Gaussian::sig_m( mes[5.20,5.30], mB0[5.20,5.30], sigmB0[0.01,0.05] ))");

   // Make background component: A plain decay function in t times an Argus function in the reconstructed mass
   w->factory("PROD::bkg(  Decay::bkg_t( dt, tau, gm, DoubleSided),"
              "ArgusBG::bkg_m( mes, 5.291, k[-100,-10]))");

   // Make composite model from the signal and background component
   w->factory("SUM::model( Nsig[5000,0,10000]*sig, NBkg[500,0,10000]*bkg )");

   // E x a m p l e   o f   R o o S i m W S T o o l   i n t e r f a c e
   // ------------------------------------------------------------------

   // Introduce a flavour tagging category tagCat as observable with 4 states corresponding
   // to 4 flavour tagging techniques with different performance that require different
   // parameterizations of the fit model
   //
   // RooSimWSTool operation:
   //     - Make 4 clones of model (for each tagCat) state, that will gain an individual
   //       copy of parameters w,dw and biasC. The other parameters remain common
   //     - Make a simultaneous pdf of the 4 clones assigning each to the appropriate
   //       state of the tagCat index category

   // RooSimWSTool is interfaced as meta-type SIMCLONE in the factory. The $SplitParam()
   // argument maps to the SplitParam() named argument in the RooSimWSTool constructor
   w->factory("SIMCLONE::model_sim( model, $SplitParam({w,dw,biasC},tagCat[Lep,Kao,NT1,NT2]))");

   // E x a m p l e   o f   R o o C u s t o m i z e r   i n t e r f a c e
   // -------------------------------------------------------------------
   //
   // Class RooCustomizer makes clones of existing pdfs with certain prescribed
   // modifications (branch of leaf node replacements)
   //
   // Here we take our model (the original before RooSimWSTool modifications)
   // and request that the parameter w (the mistag rate) is replaced with
   // an expression-based function that calculates w in terms of the Dilution
   // parameter D that is defined as D = 1-2*w

   // Make a clone model_D of original 'model' replacing 'w' with 'expr('0.5-D/2',D[0,1])'
   w->factory("EDIT::model_D(model, w=expr('0.5-D/2',D[0,1]) )");

   // Print workspace contents
   w->Print();

   // Make workspace visible on command line
   gDirectory->Add(w);
}
