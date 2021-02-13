/// \file
/// \ingroup tutorial_roostats
/// \notebook
/// This example is a generalization of the on/off problem.
///
///  This example is a generalization of the on/off problem.
/// It's a common setup for SUSY searches.  Imagine that one has two
/// variables "x" and "y" (eg. missing ET and SumET), see figure.
/// The signal region has high values of both of these variables (top right).
/// One can see low values of "x" or "y" acting as side-bands.  If we
/// just used "y" as a sideband, we would have the on/off problem.
///  - In the signal region we observe non events and expect s+b events.
///  - In the region with low values of "y" (bottom right)
///    we observe noff events and expect tau*b events.
/// Note the significance of tau.  In the background only case:
///
/// ~~~{.cpp}
///    tau ~ <expectation off> / <expectation on>
/// ~~~
///
/// If tau is known, this model is sufficient, but often tau is not known exactly.
/// So one can use low values of "x" as an additional constraint for tau.
/// Note that this technique critically depends on the notion that the
/// joint distribution for "x" and "y" can be factorized.
/// Generally, these regions have many events, so it the ratio can be
/// measured very precisely there.  So we extend the model to describe the
/// left two boxes... denoted with "bar".
///   - In the upper left we observe nonbar events and expect bbar events
///   - In the bottom left we observe noffbar events and expect tau bbar events
/// Note again we have:
///
/// ~~~{.cpp}
///    tau ~ <expectation off bar> / <expectation on bar>
/// ~~~
///
/// One can further expand the model to account for the systematic associated
/// to assuming the distribution of "x" and "y" factorizes (eg. that
/// tau is the same for off/on and offbar/onbar). This can be done in several
/// ways, but here we introduce an additional parameter rho, which so that
/// one set of models will use tau and the other tau*rho. The choice is arbitrary,
/// but it has consequences on the numerical stability of the algorithms.
/// The "bar" measurements typically have more events (& smaller relative errors).
/// If we choose
///
/// ~~~{.cpp}
/// <expectation noffbar> = tau * rho * <expectation noonbar>
/// ~~~
///
/// the product tau*rho will be known very precisely (~1/sqrt(bbar)) and the contour
/// in those parameters will be narrow and have a non-trivial tau~1/rho shape.
/// However, if we choose to put rho on the non/noff measurements (where the
/// product will have an error `~1/sqrt(b))`, the contours will be more amenable
/// to numerical techniques.  Thus, here we choose to define:
///
/// ~~~{.cpp}
///    tau := <expectation off bar> / (<expectation on bar>)
///    rho := <expectation off> / (<expectation on> * tau)
///
/// ^ y
/// |
/// |---------------------------+
/// |               |           |
/// |     nonbar    |    non    |
/// |      bbar     |    s+b    |
/// |               |           |
/// |---------------+-----------|
/// |               |           |
/// |    noffbar    |    noff   |
/// |    tau bbar   | tau b rho |
/// |               |           |
/// +-----------------------------> x
/// ~~~
///
/// Left in this way, the problem is under-constrained.  However, one may
/// have some auxiliary measurement (usually based on Monte Carlo) to
/// constrain rho.  Let us call this auxiliary measurement that gives
/// the nominal value of rho "rhonom".  Thus, there is a 'constraint' term in
/// the full model: P(rhonom | rho).  In this case, we consider a Gaussian
/// constraint with standard deviation sigma.
///
/// In the example, the initial values of the parameters are:
///
/// ~~~{.cpp}
///   - s    = 40
///   - b    = 100
///   - tau  = 5
///   - bbar = 1000
///   - rho  = 1
///   (sigma for rho = 20%)
/// ~~~
///
/// and in the toy dataset:
///
/// ~~~{.cpp}
///    - non = 139
///    - noff = 528
///    - nonbar = 993
///    - noffbar = 4906
///    - rhonom = 1.27824
/// ~~~
///
/// Note, the covariance matrix of the parameters has large off-diagonal terms.
/// Clearly s,b are anti-correlated.  Similarly, since noffbar >> nonbar, one would
/// expect bbar,tau to be anti-correlated.
///
/// This can be seen below.
///
/// ~~~{.cpp}
///             GLOBAL      b    bbar   rho      s     tau
///         b  0.96820   1.000  0.191 -0.942 -0.762 -0.209
///      bbar  0.91191   0.191  1.000  0.000 -0.146 -0.912
///       rho  0.96348  -0.942  0.000  1.000  0.718 -0.000
///         s  0.76250  -0.762 -0.146  0.718  1.000  0.160
///       tau  0.92084  -0.209 -0.912 -0.000  0.160  1.000
/// ~~~
///
/// Similarly, since tau*rho appears as a product, we expect rho,tau
/// to be anti-correlated. When the error on rho is significantly
/// larger than 1/sqrt(bbar), tau is essentially known and the
/// correlation is minimal (tau mainly cares about bbar, and rho about b,s).
/// In the alternate parametrization (bbar* tau * rho) the correlation coefficient
/// for rho,tau is large (and negative).
///
/// The code below uses best-practices for RooFit & RooStats as of June 2010.
///
/// It proceeds as follows:
///  - create a workspace to hold the model
///  - use workspace factory to quickly create the terms of the model
///  - use workspace factory to define total model (a prod pdf)
///  - create a RooStats ModelConfig to specify observables, parameters of interest
///  - add to the ModelConfig a prior on the parameters for Bayesian techniques
///    note, the pdf it is factorized for parameters of interest & nuisance params
///  - visualize the model
///  - write the workspace to a file
///  - use several of RooStats IntervalCalculators & compare results
///
/// \macro_image
/// \macro_output
/// \macro_code
///
/// \authors  authors: Kyle Cranmer, Tanja Rommerskirchen

#include "TStopwatch.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "RooPlot.h"
#include "RooAbsPdf.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooGlobalFunc.h"
#include "RooFitResult.h"
#include "RooRandom.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/BayesianCalculator.h"
#include "RooStats/MCMCCalculator.h"
#include "RooStats/MCMCInterval.h"
#include "RooStats/MCMCIntervalPlot.h"
#include "RooStats/ProposalHelper.h"
#include "RooStats/SimpleInterval.h"
#include "RooStats/FeldmanCousins.h"
#include "RooStats/PointSetInterval.h"

using namespace RooFit;
using namespace RooStats;

void FourBinInstructional(bool doBayesian = false, bool doFeldmanCousins = false, bool doMCMC = false)
{

   // let's time this challenging example
   TStopwatch t;
   t.Start();

   // set RooFit random seed for reproducible results
   RooRandom::randomGenerator()->SetSeed(4357);

   // make model
   RooWorkspace *wspace = new RooWorkspace("wspace");
   wspace->factory("Poisson::on(non[0,1000], sum::splusb(s[40,0,100],b[100,0,300]))");
   wspace->factory("Poisson::off(noff[0,5000], prod::taub(b,tau[5,3,7],rho[1,0,2]))");
   wspace->factory("Poisson::onbar(nonbar[0,10000], bbar[1000,500,2000])");
   wspace->factory("Poisson::offbar(noffbar[0,1000000], prod::lambdaoffbar(bbar, tau))");
   wspace->factory("Gaussian::mcCons(rhonom[1.,0,2], rho, sigma[.2])");
   wspace->factory("PROD::model(on,off,onbar,offbar,mcCons)");
   wspace->defineSet("obs", "non,noff,nonbar,noffbar,rhonom");

   wspace->factory("Uniform::prior_poi({s})");
   wspace->factory("Uniform::prior_nuis({b,bbar,tau, rho})");
   wspace->factory("PROD::prior(prior_poi,prior_nuis)");

   // ----------------------------------
   // Control some interesting variations
   // define parameers of interest
   // for 1-d plots
   wspace->defineSet("poi", "s");
   wspace->defineSet("nuis", "b,tau,rho,bbar");
   // for 2-d plots to inspect correlations:
   //  wspace->defineSet("poi","s,rho");

   // test simpler cases where parameters are known.
   //  wspace->var("tau")->setConstant();
   //  wspace->var("rho")->setConstant();
   //  wspace->var("b")->setConstant();
   //  wspace->var("bbar")->setConstant();

   // inspect workspace
   //  wspace->Print();

   // ----------------------------------------------------------
   // Generate toy data
   // generate toy data assuming current value of the parameters
   // import into workspace.
   // add Verbose() to see how it's being generated
   RooDataSet *data = wspace->pdf("model")->generate(*wspace->set("obs"), 1);
   //  data->Print("v");
   wspace->import(*data);

   // ----------------------------------
   // Now the statistical tests
   // model config
   ModelConfig *modelConfig = new ModelConfig("FourBins");
   modelConfig->SetWorkspace(*wspace);
   modelConfig->SetPdf(*wspace->pdf("model"));
   modelConfig->SetPriorPdf(*wspace->pdf("prior"));
   modelConfig->SetParametersOfInterest(*wspace->set("poi"));
   modelConfig->SetNuisanceParameters(*wspace->set("nuis"));
   wspace->import(*modelConfig);
   wspace->writeToFile("FourBin.root");

   // -------------------------------------------------
   // If you want to see the covariance matrix uncomment
   //  wspace->pdf("model")->fitTo(*data);

   // use ProfileLikelihood
   ProfileLikelihoodCalculator plc(*data, *modelConfig);
   plc.SetConfidenceLevel(0.95);
   LikelihoodInterval *plInt = plc.GetInterval();
   RooFit::MsgLevel msglevel = RooMsgService::instance().globalKillBelow();
   RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
   plInt->LowerLimit(*wspace->var("s")); // get ugly print out of the way. Fix.
   RooMsgService::instance().setGlobalKillBelow(msglevel);

   // use FeldmaCousins (takes ~20 min)
   FeldmanCousins fc(*data, *modelConfig);
   fc.SetConfidenceLevel(0.95);
   // number counting: dataset always has 1 entry with N events observed
   fc.FluctuateNumDataEntries(false);
   fc.UseAdaptiveSampling(true);
   fc.SetNBins(40);
   PointSetInterval *fcInt = NULL;
   if (doFeldmanCousins) {                          // takes 7 minutes
      fcInt = (PointSetInterval *)fc.GetInterval(); // fix cast
   }

   // use BayesianCalculator (only 1-d parameter of interest, slow for this problem)
   BayesianCalculator bc(*data, *modelConfig);
   bc.SetConfidenceLevel(0.95);
   SimpleInterval *bInt = NULL;
   if (doBayesian && wspace->set("poi")->getSize() == 1) {
      bInt = bc.GetInterval();
   } else {
      cout << "Bayesian Calc. only supports on parameter of interest" << endl;
   }

   // use MCMCCalculator  (takes about 1 min)
   // Want an efficient proposal function, so derive it from covariance
   // matrix of fit
   RooFitResult *fit = wspace->pdf("model")->fitTo(*data, Save());
   ProposalHelper ph;
   ph.SetVariables((RooArgSet &)fit->floatParsFinal());
   ph.SetCovMatrix(fit->covarianceMatrix());
   ph.SetUpdateProposalParameters(kTRUE); // auto-create mean vars and add mappings
   ph.SetCacheSize(100);
   ProposalFunction *pf = ph.GetProposalFunction();

   MCMCCalculator mc(*data, *modelConfig);
   mc.SetConfidenceLevel(0.95);
   mc.SetProposalFunction(*pf);
   mc.SetNumBurnInSteps(500); // first N steps to be ignored as burn-in
   mc.SetNumIters(50000);
   mc.SetLeftSideTailFraction(0.5); // make a central interval
   MCMCInterval *mcInt = NULL;
   if (doMCMC)
      mcInt = mc.GetInterval();

   // ----------------------------------
   // Make some  plots
   TCanvas *c1 = (TCanvas *)gROOT->Get("c1");
   if (!c1)
      c1 = new TCanvas("c1");

   if (doBayesian && doMCMC) {
      c1->Divide(3);
      c1->cd(1);
   } else if (doBayesian || doMCMC) {
      c1->Divide(2);
      c1->cd(1);
   }

   LikelihoodIntervalPlot *lrplot = new LikelihoodIntervalPlot(plInt);
   lrplot->Draw();

   if (doBayesian && wspace->set("poi")->getSize() == 1) {
      c1->cd(2);
      // the plot takes a long time and print lots of error
      // using a scan it is better
      bc.SetScanOfPosterior(20);
      RooPlot *bplot = bc.GetPosteriorPlot();
      bplot->Draw();
   }

   if (doMCMC) {
      if (doBayesian && wspace->set("poi")->getSize() == 1)
         c1->cd(3);
      else
         c1->cd(2);
      MCMCIntervalPlot mcPlot(*mcInt);
      mcPlot.Draw();
   }

   // ----------------------------------
   // querry intervals
   cout << "Profile Likelihood interval on s = [" << plInt->LowerLimit(*wspace->var("s")) << ", "
        << plInt->UpperLimit(*wspace->var("s")) << "]" << endl;
   // Profile Likelihood interval on s = [12.1902, 88.6871]

   if (doBayesian && wspace->set("poi")->getSize() == 1) {
      cout << "Bayesian interval on s = [" << bInt->LowerLimit() << ", " << bInt->UpperLimit() << "]" << endl;
   }

   if (doFeldmanCousins) {
      cout << "Feldman Cousins interval on s = [" << fcInt->LowerLimit(*wspace->var("s")) << ", "
           << fcInt->UpperLimit(*wspace->var("s")) << "]" << endl;
      // Feldman Cousins interval on s = [18.75 +/- 2.45, 83.75 +/- 2.45]
   }

   if (doMCMC) {
      cout << "MCMC interval on s = [" << mcInt->LowerLimit(*wspace->var("s")) << ", "
           << mcInt->UpperLimit(*wspace->var("s")) << "]" << endl;
      // MCMC interval on s = [15.7628, 84.7266]
   }

   t.Print();
}
