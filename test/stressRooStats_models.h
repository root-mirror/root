// RooFit headers
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooDataSet.h"

// RooStats headers
#include "RooStats/ModelConfig.h"

using namespace RooFit;
using namespace RooStats;


////////////////////////////////////////////////////////////////////////////////
/// Build model

void buildSimultaneousModel(RooWorkspace *w)
{
   w->factory("sig[8,0,50]");
   w->factory("Uniform::u1(x1[0,1])");
   w->factory("Uniform::u2(x2[0,1])");
   w->factory("Gaussian::constr1(gbkg1[50,0,100], bkg1[50,0,100], 3)");
   w->factory("Gaussian::constr2(gbkg2[50,0,100], bkg2[50,0,100], 2)");

   w->factory("ExtendPdf::ext_pdf1(PROD::p1(u1,constr1), expr::n1('sig+bkg1', sig, bkg1))");
   w->factory("ExtendPdf::ext_pdf2(PROD::p2(u2,constr2), expr::n2('sig+bkg2', sig, bkg2))");
   w->factory("SIMUL::sim_pdf(index[cat1,cat2],cat1=ext_pdf1,cat2=ext_pdf2)");

   // create combined signal + background model configuration
   ModelConfig *sbModel = new ModelConfig("S+B", w);
   sbModel->SetObservables("x1,x2,index");
   sbModel->SetParametersOfInterest("sig");
   sbModel->SetGlobalObservables("gbkg1,gbkg2");
   sbModel->SetNuisanceParameters("bkg1,bkg2");
   sbModel->SetPdf("sim_pdf");
   w->import(*sbModel);

   // create combined background model configuration
   ModelConfig *bModel = new ModelConfig(*sbModel);
   bModel->SetName("B");
   w->import(*bModel);

   // define data set
   RooRandom::randomGenerator()->Rndm();  //wast a number to get a better dataset (not too high significance) and closer to expected
   RooDataSet *data = w->pdf("sim_pdf")->generate(*sbModel->GetObservables(), Extended(), Name("data"));
   w->import(*data);
}

////////////////////////////////////////////////////////////////////////////////
/// Build product model

void buildPoissonProductModel(RooWorkspace *w)
{
   w->factory("expr::comp_sig('2*sig*pow(1.2, beta)', sig[0,20], beta[-3,3])");
   w->factory("Poisson::poiss1(x[0,40], sum::splusb1(sig, bkg1[0,10]))");
   w->factory("Poisson::poiss2(y[0,120], sum::splusb2(comp_sig, bkg2[0,10]))");
   w->factory("Poisson::constr1(gbkg1[5,0,10], bkg1)");
   w->factory("Poisson::constr2(gbkg2[5,0,10], bkg2)");
   w->factory("Gaussian::constr3(beta0[0,-3,3], beta, 1)");
   w->factory("PROD::pdf(poiss1, poiss2, constr1, constr2, constr3)");

   // POI prior Pdf (for BayesianCalculator and other Bayesian methods)
   w->factory("Uniform::prior(sig)");
   // Nuisance parameters Pdf (for HybridCalculator)
   w->factory("PROD::prior_nuis(constr1,constr2,constr3)");

   // create signal + background model configuration
   ModelConfig *sbModel = new ModelConfig("S+B", w);
   sbModel->SetObservables("x,y");
   sbModel->SetGlobalObservables("beta0,gbkg1,gbkg2");
   sbModel->SetParametersOfInterest("sig");
   sbModel->SetNuisanceParameters("beta,bkg1,bkg2");
   sbModel->SetPdf("pdf");
   w->import(*sbModel);

   // create background model configuration
   ModelConfig *bModel = new ModelConfig(*sbModel);
   bModel->SetName("B");
   w->import(*bModel);

   // define data set
   RooDataSet *data = new RooDataSet("data", "data", *sbModel->GetObservables());
   w->import(*data);
}


//__________________________________________________________________________________
// Insightful comments on model courtesy of Kyle Cranmer, Wouter Verkerke, Sven Kreiss
//    from $ROOTSYS/tutorials/roostats/HybridInstructional.C
void buildOnOffModel(RooWorkspace *w)
{
   // Build model for prototype on/off problem
   // Poiss(x | s+b) * Poiss(y | tau b )
   w->factory("Poisson::on_pdf(n_on[0,300],sum::splusb(sig[0,100],bkg[0,200]))");
   w->factory("Poisson::off_pdf(n_off[0,1100],prod::taub(tau[0.1,5.0],bkg))");
   w->factory("PROD::prod_pdf(on_pdf, off_pdf)");

   // construct the Bayesian-averaged model (eg. a projection pdf)
   // p'(x|s) = \int db p(x|s+b) * [ p(y|b) * prior(b) ]
   w->factory("Uniform::prior(bkg)");
   w->factory("PROJ::averagedModel(PROD::foo(on_pdf|bkg,off_pdf,prior),bkg)") ;

   // create signal + background model configuration
   ModelConfig *sbModel = new ModelConfig("S+B", w);
   sbModel->SetPdf("prod_pdf");
   sbModel->SetObservables("n_on,n_off");
   sbModel->SetParametersOfInterest("sig");
   sbModel->SetNuisanceParameters("bkg");
   w->import(*sbModel);

   // create background model configuration
   ModelConfig *bModel = new ModelConfig(*sbModel);
   bModel->SetName("B");
   w->import(*bModel);

   // alternate priors
   w->factory("Gaussian::gauss_prior(bkg, n_off, expr::sqrty('sqrt(n_off)', n_off))");
   w->factory("Lognormal::lognorm_prior(bkg, n_off, expr::kappa('1+1./sqrt(n_off)',n_off))");

   // define data set
   RooDataSet *data = new RooDataSet("data", "data", *sbModel->GetObservables());
   w->import(*data);
}


void buildPoissonEfficiencyModel(RooWorkspace *w)
{

   // build models
   w->factory("Gaussian::constrb(b0[-5,5], b1[-5,5], 1)");
   w->factory("Gaussian::constre(e0[-5,5], e1[-5,5], 1)");
   w->factory("expr::bkg('5 * pow(1.3, b1)', b1)"); // background
   w->factory("expr::eff('0.5 * pow(1.2, e1)', e1)"); // efficiency
   w->factory("expr::esb('eff * sig + bkg', eff, bkg, sig[0,50])");
   w->factory("Poisson::poiss(x[0,50], esb)");
   w->factory("PROD::pdf(poiss, constrb, constre)");

   // create model configuration
   ModelConfig *sbModel = new ModelConfig("S+B", w);
   sbModel->SetObservables("x");
   sbModel->SetParametersOfInterest("sig");
   sbModel->SetNuisanceParameters("b1,e1");
   sbModel->SetGlobalObservables("b0,e0");
   sbModel->SetPdf("pdf");
   w->import(*sbModel);

   ModelConfig *bModel = new ModelConfig(*sbModel);
   bModel->SetName("B");
   w->import(*bModel);

   // define data set and import it into workspace
   RooDataSet *data = new RooDataSet("data", "data", *sbModel->GetObservables());
   w->import(*data);
}


