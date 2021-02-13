/// \file
/// \ingroup tutorial_roostats
/// \notebook -js
/// SPlot tutorial
///
/// This tutorial shows an example of using SPlot to unfold two distributions.
/// The physics context for the example is that we want to know
/// the isolation distribution for real electrons from Z events
/// and fake electrons from QCD.  Isolation is our 'control' variable.
/// To unfold them, we need a model for an uncorrelated variable that
/// discriminates between Z and QCD.  To do this, we use the invariant
/// mass of two electrons.  We model the Z with a Gaussian and the QCD
/// with a falling exponential.
///
/// Note, since we don't have real data in this tutorial, we need to generate
/// toy data.  To do that we need a model for the isolation variable for
/// both Z and QCD.  This is only used to generate the toy data, and would
/// not be needed if we had real data.
///
/// \macro_image
/// \macro_code
/// \macro_output
///
/// \author Kyle Cranmer

#include "RooRealVar.h"
#include "RooStats/SPlot.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "TCanvas.h"
#include "RooAbsPdf.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooConstVar.h"
#include <iomanip>

// use this order for safety on library loading
using namespace RooFit;
using namespace RooStats;

// see below for implementation
void AddModel(RooWorkspace *);
void AddData(RooWorkspace *);
void DoSPlot(RooWorkspace *);
void MakePlots(RooWorkspace *);

void rs301_splot()
{

   // Create a new workspace to manage the project.
   RooWorkspace *wspace = new RooWorkspace("myWS");

   // add the signal and background models to the workspace.
   // Inside this function you will find a description of our model.
   AddModel(wspace);

   // add some toy data to the workspace
   AddData(wspace);

   // inspect the workspace if you wish
   //  wspace->Print();

   // do sPlot.
   // This will make a new dataset with sWeights added for every event.
   DoSPlot(wspace);

   // Make some plots showing the discriminating variable and
   // the control variable after unfolding.
   MakePlots(wspace);

   // cleanup
   delete wspace;
}

//____________________________________
void AddModel(RooWorkspace *ws)
{

   // Make models for signal (Higgs) and background (Z+jets and QCD)
   // In real life, this part requires an intelligent modeling
   // of signal and background -- this is only an example.

   // set range of observable
   Double_t lowRange = 0., highRange = 200.;

   // make a RooRealVar for the observables
   RooRealVar invMass("invMass", "M_{inv}", lowRange, highRange, "GeV");
   RooRealVar isolation("isolation", "isolation", 0., 20., "GeV");

   // --------------------------------------
   // make 2-d model for Z including the invariant mass
   // distribution  and an isolation distribution which we want to
   // unfold from QCD.
   std::cout << "make z model" << std::endl;
   // mass model for Z
   RooRealVar mZ("mZ", "Z Mass", 91.2, lowRange, highRange);
   RooRealVar sigmaZ("sigmaZ", "Width of Gaussian", 2, 0, 10, "GeV");
   RooGaussian mZModel("mZModel", "Z+jets Model", invMass, mZ, sigmaZ);
   // we know Z mass
   mZ.setConstant();
   // we leave the width of the Z free during the fit in this example.

   // isolation model for Z.  Only used to generate toy MC.
   // the exponential is of the form exp(c*x).  If we want
   // the isolation to decay an e-fold every R GeV, we use
   // c = -1/R.
   RooConstVar zIsolDecayConst("zIsolDecayConst", "z isolation decay  constant", -1);
   RooExponential zIsolationModel("zIsolationModel", "z isolation model", isolation, zIsolDecayConst);

   // make the combined Z model
   RooProdPdf zModel("zModel", "2-d model for Z", RooArgSet(mZModel, zIsolationModel));

   // --------------------------------------
   // make QCD model

   std::cout << "make qcd model" << std::endl;
   // mass model for QCD.
   // the exponential is of the form exp(c*x).  If we want
   // the mass to decay an e-fold every R GeV, we use
   // c = -1/R.
   // We can leave this parameter free during the fit.
   RooRealVar qcdMassDecayConst("qcdMassDecayConst", "Decay const for QCD mass spectrum", -0.01, -100, 100, "1/GeV");
   RooExponential qcdMassModel("qcdMassModel", "qcd Mass Model", invMass, qcdMassDecayConst);

   // isolation model for QCD.  Only used to generate toy MC
   // the exponential is of the form exp(c*x).  If we want
   // the isolation to decay an e-fold every R GeV, we use
   // c = -1/R.
   RooConstVar qcdIsolDecayConst("qcdIsolDecayConst", "Et resolution constant", -.1);
   RooExponential qcdIsolationModel("qcdIsolationModel", "QCD isolation model", isolation, qcdIsolDecayConst);

   // make the 2-d model
   RooProdPdf qcdModel("qcdModel", "2-d model for QCD", RooArgSet(qcdMassModel, qcdIsolationModel));

   // --------------------------------------
   // combined model

   // These variables represent the number of Z or QCD events
   // They will be fitted.
   RooRealVar zYield("zYield", "fitted yield for Z", 50, 0., 1000);
   RooRealVar qcdYield("qcdYield", "fitted yield for QCD", 100, 0., 1000);

   // now make the combined model
   std::cout << "make full model" << std::endl;
   RooAddPdf model("model", "z+qcd background models", RooArgList(zModel, qcdModel), RooArgList(zYield, qcdYield));

   // interesting for debugging and visualizing the model
   model.graphVizTree("fullModel.dot");

   std::cout << "import model" << std::endl;

   ws->import(model);
}

//____________________________________
void AddData(RooWorkspace *ws)
{
   // Add a toy dataset

   // how many events do we want?
   Int_t nEvents = 1000;

   // get what we need out of the workspace to make toy data
   RooAbsPdf *model = ws->pdf("model");
   RooRealVar *invMass = ws->var("invMass");
   RooRealVar *isolation = ws->var("isolation");

   // make the toy data
   std::cout << "make data set and import to workspace" << std::endl;
   RooDataSet *data = model->generate(RooArgSet(*invMass, *isolation), nEvents);

   // import data into workspace
   ws->import(*data, Rename("data"));
}

//____________________________________
void DoSPlot(RooWorkspace *ws)
{
   std::cout << "Calculate sWeights" << std::endl;

   // get what we need out of the workspace to do the fit
   RooAbsPdf *model = ws->pdf("model");
   RooRealVar *zYield = ws->var("zYield");
   RooRealVar *qcdYield = ws->var("qcdYield");
   RooDataSet *data = (RooDataSet *)ws->data("data");

   // fit the model to the data.
   model->fitTo(*data, Extended());

   // The sPlot technique requires that we fix the parameters
   // of the model that are not yields after doing the fit.
   //
   // This *could* be done with the lines below, however this is taken care of
   // by the RooStats::SPlot constructor (or more precisely the AddSWeight
   // method).
   //
   // RooRealVar* sigmaZ = ws->var("sigmaZ");
   // RooRealVar* qcdMassDecayConst = ws->var("qcdMassDecayConst");
   // sigmaZ->setConstant();
   // qcdMassDecayConst->setConstant();

   RooMsgService::instance().setSilentMode(true);

   std::cout << "\n\n------------------------------------------\nThe dataset before creating sWeights:\n";
   data->Print();

   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

   // Now we use the SPlot class to add SWeights to our data set
   // based on our model and our yield variables
   RooStats::SPlot *sData = new RooStats::SPlot("sData", "An SPlot", *data, model, RooArgList(*zYield, *qcdYield));

   std::cout << "\n\nThe dataset after creating sWeights:\n";
   data->Print();

   // Check that our weights have the desired properties

   std::cout << "\n\n------------------------------------------\n\nCheck SWeights:" << std::endl;

   std::cout << std::endl
             << "Yield of Z is\t" << zYield->getVal() << ".  From sWeights it is "
             << sData->GetYieldFromSWeight("zYield") << std::endl;

   std::cout << "Yield of QCD is\t" << qcdYield->getVal() << ".  From sWeights it is "
             << sData->GetYieldFromSWeight("qcdYield") << std::endl
             << std::endl;

   for (Int_t i = 0; i < 10; i++) {
      std::cout << "z Weight for event " << i << std::right << std::setw(12) << sData->GetSWeight(i, "zYield") << "  qcd Weight"
                << std::setw(12) << sData->GetSWeight(i, "qcdYield") << "  Total Weight" << std::setw(12) << sData->GetSumOfEventSWeight(i)
                << std::endl;
   }

   std::cout << std::endl;

   // import this new dataset with sWeights
   std::cout << "import new dataset with sWeights" << std::endl;
   ws->import(*data, Rename("dataWithSWeights"));

   RooMsgService::instance().setGlobalKillBelow(RooFit::INFO);
}

void MakePlots(RooWorkspace *ws)
{

   // Here we make plots of the discriminating variable (invMass) after the fit
   // and of the control variable (isolation) after unfolding with sPlot.
   std::cout << "make plots" << std::endl;

   // make our canvas
   TCanvas *cdata = new TCanvas("sPlot", "sPlot demo", 400, 600);
   cdata->Divide(1, 3);

   // get what we need out of the workspace
   RooAbsPdf *model = ws->pdf("model");
   RooAbsPdf *zModel = ws->pdf("zModel");
   RooAbsPdf *qcdModel = ws->pdf("qcdModel");

   RooRealVar *isolation = ws->var("isolation");
   RooRealVar *invMass = ws->var("invMass");

   // note, we get the dataset with sWeights
   RooDataSet *data = (RooDataSet *)ws->data("dataWithSWeights");

   // this shouldn't be necessary, need to fix something with workspace
   // do this to set parameters back to their fitted values.
//   model->fitTo(*data, Extended());

   // plot invMass for data with full model and individual components overlaid
   //  TCanvas* cdata = new TCanvas();
   cdata->cd(1);
   RooPlot *frame = invMass->frame();
   data->plotOn(frame);
   model->plotOn(frame, Name("FullModel"));
   model->plotOn(frame, Components(*zModel), LineStyle(kDashed), LineColor(kRed), Name("ZModel"));
   model->plotOn(frame, Components(*qcdModel), LineStyle(kDashed), LineColor(kGreen), Name("QCDModel"));

   TLegend leg(0.11, 0.5, 0.5, 0.8);
   leg.AddEntry(frame->findObject("FullModel"), "Full model", "L");
   leg.AddEntry(frame->findObject("ZModel"), "Z model", "L");
   leg.AddEntry(frame->findObject("QCDModel"), "QCD model", "L");
   leg.SetBorderSize(0);
   leg.SetFillStyle(0);

   frame->SetTitle("Fit of model to discriminating variable");
   frame->Draw();
   leg.DrawClone();

   // Now use the sWeights to show isolation distribution for Z and QCD.
   // The SPlot class can make this easier, but here we demonstrate in more
   // detail how the sWeights are used.  The SPlot class should make this
   // very easy and needs some more development.

   // Plot isolation for Z component.
   // Do this by plotting all events weighted by the sWeight for the Z component.
   // The SPlot class adds a new variable that has the name of the corresponding
   // yield + "_sw".
   cdata->cd(2);

   // create weighted data set
   RooDataSet *dataw_z = new RooDataSet(data->GetName(), data->GetTitle(), data, *data->get(), 0, "zYield_sw");

   RooPlot *frame2 = isolation->frame();
   // Since the data are weighted, we use SumW2 to compute the errors.
   dataw_z->plotOn(frame2, DataError(RooAbsData::SumW2));

   frame2->SetTitle("Isolation distribution with s weights to project out Z");
   frame2->Draw();

   // Plot isolation for QCD component.
   // Eg. plot all events weighted by the sWeight for the QCD component.
   // The SPlot class adds a new variable that has the name of the corresponding
   // yield + "_sw".
   cdata->cd(3);
   RooDataSet *dataw_qcd = new RooDataSet(data->GetName(), data->GetTitle(), data, *data->get(), 0, "qcdYield_sw");
   RooPlot *frame3 = isolation->frame();
   dataw_qcd->plotOn(frame3, DataError(RooAbsData::SumW2));

   frame3->SetTitle("Isolation distribution with s weights to project out QCD");
   frame3->Draw();

   //  cdata->SaveAs("SPlot.gif");
}
