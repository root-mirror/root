/// \file
/// \ingroup tutorial_roostats
/// \notebook -js
/// TwoSidedFrequentistUpperLimitWithBands
///
///
/// This is a standard demo that can be used with any ROOT file
/// prepared in the standard way.  You specify:
///  - name for input ROOT file
///  - name of workspace inside ROOT file that holds model and data
///  - name of ModelConfig that specifies details for calculator tools
///  - name of dataset
///
/// With default parameters the macro will attempt to run the
/// standard hist2workspace example and read the ROOT file
/// that it produces.
///
/// You may want to control:
/// ~~~{.cpp}
///   double confidenceLevel=0.95;
///   double additionalToysFac = 1.;
///   int nPointsToScan = 12;
///   int nToyMC = 200;
/// ~~~
///
/// This uses a modified version of the profile likelihood ratio as
/// a test statistic for upper limits (eg. test stat = 0 if muhat>mu).
///
/// Based on the observed data, one defines a set of parameter points
/// to be tested based on the value of the parameter of interest
/// and the conditional MLE (eg. profiled) values of the nuisance parameters.
///
/// At each parameter point, pseudo-experiments are generated using this
/// fixed reference model and then the test statistic is evaluated.
/// The auxiliary measurements (global observables) associated with the
/// constraint terms in nuisance parameters are also fluctuated in the
/// process of generating the pseudo-experiments in a frequentist manner
/// forming an 'unconditional ensemble'.  One could form a 'conditional'
/// ensemble in which these auxiliary measurements are fixed.  Note that the
/// nuisance parameters are not randomized, which is a Bayesian procedure.
/// Note, the nuisance parameters are floating in the fits.  For each point,
/// the threshold that defines the 95% acceptance region is found.  This
/// forms a "Confidence Belt".
///
/// After constructing the confidence belt, one can find the confidence
/// interval for any particular dataset by finding the intersection
/// of the observed test statistic and the confidence belt.  First
/// this is done on the observed data to get an observed 1-sided upper limt.
///
/// Finally, there expected limit and bands (from background-only) are
/// formed by generating background-only data and finding the upper limit.
/// The background-only is defined as such that the nuisance parameters are
/// fixed to their best fit value based on the data with the signal rate fixed to 0.
/// The bands are done by hand for now, will later be part of the RooStats tools.
///
/// On a technical note, this technique IS the generalization of Feldman-Cousins
/// with nuisance parameters.
///
/// Building the confidence belt can be computationally expensive.
/// Once it is built, one could save it to a file and use it in a separate step.
///
/// We can use PROOF to speed things along in parallel, however,
/// the test statistic has to be installed on the workers
/// so either turn off PROOF or include the modified test statistic
/// in your $ROOTSYS/roofit/roostats/inc directory,
/// add the additional line to the LinkDef.h file,
/// and recompile root.
///
/// Note, if you have a boundary on the parameter of interest (eg. cross-section)
/// the threshold on the two-sided test statistic starts off at moderate values and plateaus.
///
/// [#0] PROGRESS:Generation -- generated toys: 500 / 999
/// NeymanConstruction: Prog: 12/50 total MC = 39 this test stat = 0
///  SigXsecOverSM=0.69 alpha_syst1=0.136515 alpha_syst3=0.425415 beta_syst2=1.08496 [-1e+30, 0.011215]  in interval = 1
///
/// this tells you the values of the parameters being used to generate the pseudo-experiments
/// and the threshold in this case is 0.011215.  One would expect for 95% that the threshold
/// would be ~1.35 once the cross-section is far enough away from 0 that it is essentially
/// unaffected by the boundary.  As one reaches the last points in the scan, the
/// theshold starts to get artificially high.  This is because the range of the parameter in
/// the fit is the same as the range in the scan.  In the future, these should be independently
/// controlled, but they are not now.  As a result the ~50% of pseudo-experiments that have an
/// upward fluctuation end up with muhat = muMax.  Because of this, the upper range of the
/// parameter should be well above the expected upper limit... but not too high or one will
/// need a very large value of nPointsToScan to resolve the relevant region.  This can be
/// improved, but this is the first version of this script.
///
/// Important note: when the model includes external constraint terms, like a Gaussian
/// constraint to a nuisance parameter centered around some nominal value there is
/// a subtlety.  The asymptotic results are all based on the assumption that all the
/// measurements fluctuate... including the nominal values from auxiliary measurements.
/// If these do not fluctuate, this corresponds to an "conditional ensemble".  The
/// result is that the distribution of the test statistic can become very non-chi^2.
/// This results in thresholds that become very large.
///
/// \macro_image
/// \macro_output
/// \macro_code
///
/// \authors Kyle Cranmer,Contributions from Aaron Armbruster, Haoshuang Ji, Haichen Wang and Daniel Whiteson

#include "TFile.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TSystem.h"
#include <iostream>

#include "RooWorkspace.h"
#include "RooSimultaneous.h"
#include "RooAbsData.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/FeldmanCousins.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/PointSetInterval.h"
#include "RooStats/ConfidenceBelt.h"

#include "RooStats/RooStatsUtils.h"
#include "RooStats/ProfileLikelihoodTestStat.h"

using namespace RooFit;
using namespace RooStats;
using namespace std;

bool useProof = false; // flag to control whether to use Proof
int nworkers = 0;      // number of workers (default use all available cores)

// -------------------------------------------------------

void TwoSidedFrequentistUpperLimitWithBands(const char *infile = "", const char *workspaceName = "combined",
                                            const char *modelConfigName = "ModelConfig",
                                            const char *dataName = "obsData")
{

   double confidenceLevel = 0.95;
   // degrade/improve number of pseudo-experiments used to define the confidence belt.
   // value of 1 corresponds to default number of toys in the tail, which is 50/(1-confidenceLevel)
   double additionalToysFac = 0.5;
   int nPointsToScan = 20; // number of steps in the parameter of interest
   int nToyMC = 200;       // number of toys used to define the expected limit and band

   // -------------------------------------------------------
   // First part is just to access a user-defined file
   // or create the standard example file if it doesn't exist
   const char *filename = "";
   if (!strcmp(infile, "")) {
      filename = "results/example_combined_GaussExample_model.root";
      bool fileExist = !gSystem->AccessPathName(filename); // note opposite return code
      // if file does not exists generate with histfactory
      if (!fileExist) {
#ifdef _WIN32
         cout << "HistFactory file cannot be generated on Windows - exit" << endl;
         return;
#endif
         // Normally this would be run on the command line
         cout << "will run standard hist2workspace example" << endl;
         gROOT->ProcessLine(".! prepareHistFactory .");
         gROOT->ProcessLine(".! hist2workspace config/example.xml");
         cout << "\n\n---------------------" << endl;
         cout << "Done creating example input" << endl;
         cout << "---------------------\n\n" << endl;
      }

   } else
      filename = infile;

   // Try to open the file
   TFile *inputFile = TFile::Open(filename);

   // if input file was specified byt not found, quit
   if (!inputFile) {
      cout << "StandardRooStatsDemoMacro: Input file " << filename << " is not found" << endl;
      return;
   }

   // -------------------------------------------------------
   // Now get the data and workspace

   // get the workspace out of the file
   RooWorkspace *w = (RooWorkspace *)inputFile->Get(workspaceName);
   if (!w) {
      cout << "workspace not found" << endl;
      return;
   }

   // get the modelConfig out of the file
   ModelConfig *mc = (ModelConfig *)w->obj(modelConfigName);

   // get the modelConfig out of the file
   RooAbsData *data = w->data(dataName);

   // make sure ingredients are found
   if (!data || !mc) {
      w->Print();
      cout << "data or ModelConfig was not found" << endl;
      return;
   }

   cout << "Found data and ModelConfig:" << endl;
   mc->Print();

   // -------------------------------------------------------
   // Now get the POI for convenience
   // you may want to adjust the range of your POI
   RooRealVar *firstPOI = (RooRealVar *)mc->GetParametersOfInterest()->first();
   /*  firstPOI->setMin(0);*/
   /*  firstPOI->setMax(10);*/

   // -------------------------------------------------------
   // create and use the FeldmanCousins tool
   // to find and plot the 95% confidence interval
   // on the parameter of interest as specified
   // in the model config
   // REMEMBER, we will change the test statistic
   // so this is NOT a Feldman-Cousins interval
   FeldmanCousins fc(*data, *mc);
   fc.SetConfidenceLevel(confidenceLevel);
   fc.AdditionalNToysFactor(additionalToysFac); // improve sampling that defines confidence belt
   //  fc.UseAdaptiveSampling(true); // speed it up a bit, but don't use for expected limits
   fc.SetNBins(nPointsToScan); // set how many points per parameter of interest to scan
   fc.CreateConfBelt(true);    // save the information in the belt for plotting

   // -------------------------------------------------------
   // Feldman-Cousins is a unified limit by definition
   // but the tool takes care of a few things for us like which values
   // of the nuisance parameters should be used to generate toys.
   // so let's just change the test statistic and realize this is
   // no longer "Feldman-Cousins" but is a fully frequentist Neyman-Construction.
   //  fc.GetTestStatSampler()->SetTestStatistic(&onesided);
   // ((ToyMCSampler*) fc.GetTestStatSampler())->SetGenerateBinned(true);
   ToyMCSampler *toymcsampler = (ToyMCSampler *)fc.GetTestStatSampler();
   ProfileLikelihoodTestStat *testStat = dynamic_cast<ProfileLikelihoodTestStat *>(toymcsampler->GetTestStatistic());

   // Since this tool needs to throw toy MC the PDF needs to be
   // extended or the tool needs to know how many entries in a dataset
   // per pseudo experiment.
   // In the 'number counting form' where the entries in the dataset
   // are counts, and not values of discriminating variables, the
   // datasets typically only have one entry and the PDF is not
   // extended.
   if (!mc->GetPdf()->canBeExtended()) {
      if (data->numEntries() == 1)
         fc.FluctuateNumDataEntries(false);
      else
         cout << "Not sure what to do about this model" << endl;
   }

   // We can use PROOF to speed things along in parallel
   // However, the test statistic has to be installed on the workers
   // so either turn off PROOF or include the modified test statistic
   // in your $ROOTSYS/roofit/roostats/inc directory,
   // add the additional line to the LinkDef.h file,
   // and recompile root.
   if (useProof) {
      ProofConfig pc(*w, nworkers, "", false);
      toymcsampler->SetProofConfig(&pc); // enable proof
   }

   if (mc->GetGlobalObservables()) {
      cout << "will use global observables for unconditional ensemble" << endl;
      mc->GetGlobalObservables()->Print();
      toymcsampler->SetGlobalObservables(*mc->GetGlobalObservables());
   }

   // Now get the interval
   PointSetInterval *interval = fc.GetInterval();
   ConfidenceBelt *belt = fc.GetConfidenceBelt();

   // print out the interval on the first Parameter of Interest
   cout << "\n95% interval on " << firstPOI->GetName() << " is : [" << interval->LowerLimit(*firstPOI) << ", "
        << interval->UpperLimit(*firstPOI) << "] " << endl;

   // get observed UL and value of test statistic evaluated there
   RooArgSet tmpPOI(*firstPOI);
   double observedUL = interval->UpperLimit(*firstPOI);
   firstPOI->setVal(observedUL);
   double obsTSatObsUL = fc.GetTestStatSampler()->EvaluateTestStatistic(*data, tmpPOI);

   // Ask the calculator which points were scanned
   RooDataSet *parameterScan = (RooDataSet *)fc.GetPointsToScan();
   RooArgSet *tmpPoint;

   // make a histogram of parameter vs. threshold
   TH1F *histOfThresholds =
      new TH1F("histOfThresholds", "", parameterScan->numEntries(), firstPOI->getMin(), firstPOI->getMax());
   histOfThresholds->GetXaxis()->SetTitle(firstPOI->GetName());
   histOfThresholds->GetYaxis()->SetTitle("Threshold");

   // loop through the points that were tested and ask confidence belt
   // what the upper/lower thresholds were.
   // For FeldmanCousins, the lower cut off is always 0
   for (Int_t i = 0; i < parameterScan->numEntries(); ++i) {
      tmpPoint = (RooArgSet *)parameterScan->get(i)->clone("temp");
      // cout <<"get threshold"<<endl;
      double arMax = belt->GetAcceptanceRegionMax(*tmpPoint);
      double poiVal = tmpPoint->getRealValue(firstPOI->GetName());
      histOfThresholds->Fill(poiVal, arMax);
   }
   TCanvas *c1 = new TCanvas();
   c1->Divide(2);
   c1->cd(1);
   histOfThresholds->SetMinimum(0);
   histOfThresholds->Draw();
   c1->cd(2);

   // -------------------------------------------------------
   // Now we generate the expected bands and power-constraint

   // First: find parameter point for mu=0, with conditional MLEs for nuisance parameters
   RooAbsReal *nll = mc->GetPdf()->createNLL(*data);
   RooAbsReal *profile = nll->createProfile(*mc->GetParametersOfInterest());
   firstPOI->setVal(0.);
   profile->getVal(); // this will do fit and set nuisance parameters to profiled values
   RooArgSet *poiAndNuisance = new RooArgSet();
   if (mc->GetNuisanceParameters())
      poiAndNuisance->add(*mc->GetNuisanceParameters());
   poiAndNuisance->add(*mc->GetParametersOfInterest());
   w->saveSnapshot("paramsToGenerateData", *poiAndNuisance);
   RooArgSet *paramsToGenerateData = (RooArgSet *)poiAndNuisance->snapshot();
   cout << "\nWill use these parameter points to generate pseudo data for bkg only" << endl;
   paramsToGenerateData->Print("v");

   RooArgSet unconditionalObs;
   unconditionalObs.add(*mc->GetObservables());
   unconditionalObs.add(*mc->GetGlobalObservables()); // comment this out for the original conditional ensemble

   double CLb = 0;
   double CLbinclusive = 0;

   // Now we generate background only and find distribution of upper limits
   TH1F *histOfUL = new TH1F("histOfUL", "", 100, 0, firstPOI->getMax());
   histOfUL->GetXaxis()->SetTitle("Upper Limit (background only)");
   histOfUL->GetYaxis()->SetTitle("Entries");
   for (int imc = 0; imc < nToyMC; ++imc) {

      // set parameters back to values for generating pseudo data
      //    cout << "\n get current nuis, set vals, print again" << endl;
      w->loadSnapshot("paramsToGenerateData");
      //    poiAndNuisance->Print("v");

      RooDataSet *toyData = 0;
      // now generate a toy dataset for the main measurement
      if (!mc->GetPdf()->canBeExtended()) {
         if (data->numEntries() == 1)
            toyData = mc->GetPdf()->generate(*mc->GetObservables(), 1);
         else
            cout << "Not sure what to do about this model" << endl;
      } else {
         //      cout << "generating extended dataset"<<endl;
         toyData = mc->GetPdf()->generate(*mc->GetObservables(), Extended());
      }

      // generate global observables
      // need to be careful for simpdf.
      // In ROOT 5.28 there is a problem with generating global observables
      // with a simultaneous PDF.  In 5.29 there is a solution with
      // RooSimultaneous::generateSimGlobal, but this may change to
      // the standard generate interface in 5.30.

      RooSimultaneous *simPdf = dynamic_cast<RooSimultaneous *>(mc->GetPdf());
      if (!simPdf) {
         RooDataSet *one = mc->GetPdf()->generate(*mc->GetGlobalObservables(), 1);
         const RooArgSet *values = one->get();
         RooArgSet *allVars = mc->GetPdf()->getVariables();
         *allVars = *values;
         delete allVars;
         delete one;
      } else {
         RooDataSet *one = simPdf->generateSimGlobal(*mc->GetGlobalObservables(), 1);
         const RooArgSet *values = one->get();
         RooArgSet *allVars = mc->GetPdf()->getVariables();
         *allVars = *values;
         delete allVars;
         delete one;
      }

      // get test stat at observed UL in observed data
      firstPOI->setVal(observedUL);
      double toyTSatObsUL = fc.GetTestStatSampler()->EvaluateTestStatistic(*toyData, tmpPOI);
      //    toyData->get()->Print("v");
      //    cout <<"obsTSatObsUL " <<obsTSatObsUL << "toyTS " << toyTSatObsUL << endl;
      if (obsTSatObsUL < toyTSatObsUL) // not sure about <= part yet
         CLb += (1.) / nToyMC;
      if (obsTSatObsUL <= toyTSatObsUL) // not sure about <= part yet
         CLbinclusive += (1.) / nToyMC;

      // loop over points in belt to find upper limit for this toy data
      double thisUL = 0;
      for (Int_t i = 0; i < parameterScan->numEntries(); ++i) {
         tmpPoint = (RooArgSet *)parameterScan->get(i)->clone("temp");
         double arMax = belt->GetAcceptanceRegionMax(*tmpPoint);
         firstPOI->setVal(tmpPoint->getRealValue(firstPOI->GetName()));
         //   double thisTS = profile->getVal();
         double thisTS = fc.GetTestStatSampler()->EvaluateTestStatistic(*toyData, tmpPOI);

         //   cout << "poi = " << firstPOI->getVal()
         // << " max is " << arMax << " this profile = " << thisTS << endl;
         //      cout << "thisTS = " << thisTS<<endl;
         if (thisTS <= arMax) {
            thisUL = firstPOI->getVal();
         } else {
            break;
         }
      }

      histOfUL->Fill(thisUL);

      // for few events, data is often the same, and UL is often the same
      //    cout << "thisUL = " << thisUL<<endl;

      delete toyData;
   }
   histOfUL->Draw();
   c1->SaveAs("two-sided_upper_limit_output.pdf");

   // if you want to see a plot of the sampling distribution for a particular scan point:
   /*
   SamplingDistPlot sampPlot;
   int indexInScan = 0;
   tmpPoint = (RooArgSet*) parameterScan->get(indexInScan)->clone("temp");
   firstPOI->setVal( tmpPoint->getRealValue(firstPOI->GetName()) );
   toymcsampler->SetParametersForTestStat(tmpPOI);
   SamplingDistribution* samp = toymcsampler->GetSamplingDistribution(*tmpPoint);
   sampPlot.AddSamplingDistribution(samp);
   sampPlot.Draw();
      */

   // Now find bands and power constraint
   Double_t *bins = histOfUL->GetIntegral();
   TH1F *cumulative = (TH1F *)histOfUL->Clone("cumulative");
   cumulative->SetContent(bins);
   double band2sigDown = 0, band1sigDown = 0, bandMedian = 0, band1sigUp = 0, band2sigUp = 0;
   for (int i = 1; i <= cumulative->GetNbinsX(); ++i) {
      if (bins[i] < RooStats::SignificanceToPValue(2))
         band2sigDown = cumulative->GetBinCenter(i);
      if (bins[i] < RooStats::SignificanceToPValue(1))
         band1sigDown = cumulative->GetBinCenter(i);
      if (bins[i] < 0.5)
         bandMedian = cumulative->GetBinCenter(i);
      if (bins[i] < RooStats::SignificanceToPValue(-1))
         band1sigUp = cumulative->GetBinCenter(i);
      if (bins[i] < RooStats::SignificanceToPValue(-2))
         band2sigUp = cumulative->GetBinCenter(i);
   }
   cout << "-2 sigma  band " << band2sigDown << endl;
   cout << "-1 sigma  band " << band1sigDown << " [Power Constraint)]" << endl;
   cout << "median of band " << bandMedian << endl;
   cout << "+1 sigma  band " << band1sigUp << endl;
   cout << "+2 sigma  band " << band2sigUp << endl;

   // print out the interval on the first Parameter of Interest
   cout << "\nobserved 95% upper-limit " << interval->UpperLimit(*firstPOI) << endl;
   cout << "CLb strict [P(toy>obs|0)] for observed 95% upper-limit " << CLb << endl;
   cout << "CLb inclusive [P(toy>=obs|0)] for observed 95% upper-limit " << CLbinclusive << endl;

   delete profile;
   delete nll;
}
