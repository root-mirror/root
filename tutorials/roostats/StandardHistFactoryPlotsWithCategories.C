/// \file
/// \ingroup tutorial_roostats
/// \notebook -js
/// StandardHistFactoryPlotsWithCategories
///
///  This is a standard demo that can be used with any ROOT file
///  prepared in the standard way.  You specify:
///  - name for input ROOT file
///  - name of workspace inside ROOT file that holds model and data
///  - name of ModelConfig that specifies details for calculator tools
///  - name of dataset
///
///  With default parameters the macro will attempt to run the
///  standard hist2workspace example and read the ROOT file
///  that it produces.
///
///  The macro will scan through all the categories in a simPdf find the corresponding
///  observable.  For each category, it will loop through each of the nuisance parameters
///  and plot
///  - the data
///  - the nominal model (blue)
///  - the +Nsigma (red)
///  - the -Nsigma (green)
///
///  You can specify how many sigma to vary by changing nSigmaToVary.
///  You can also change the signal rate by changing muVal.
///
///  The script produces a lot plots, you can merge them by doing:
/// ~~~{.cpp}
///  gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=merged.pdf `ls *pdf`
/// ~~~
///
/// \macro_image
/// \macro_output
/// \macro_code
///
/// \author Kyle Cranmer

#include "TFile.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TList.h"
#include "TMath.h"
#include "TSystem.h"
#include "RooWorkspace.h"
#include "RooAbsData.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileInspector.h"

using namespace RooFit;
using namespace RooStats;
using namespace std;

void StandardHistFactoryPlotsWithCategories(const char *infile = "", const char *workspaceName = "combined",
                                            const char *modelConfigName = "ModelConfig",
                                            const char *dataName = "obsData")
{

   double nSigmaToVary = 5.;
   double muVal = 0;
   bool doFit = false;

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
   TFile *file = TFile::Open(filename);

   // if input file was specified byt not found, quit
   if (!file) {
      cout << "StandardRooStatsDemoMacro: Input file " << filename << " is not found" << endl;
      return;
   }

   // -------------------------------------------------------
   // Tutorial starts here
   // -------------------------------------------------------

   // get the workspace out of the file
   RooWorkspace *w = (RooWorkspace *)file->Get(workspaceName);
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

   // -------------------------------------------------------
   // now use the profile inspector

   RooRealVar *obs = (RooRealVar *)mc->GetObservables()->first();
   TList *list = new TList();

   RooRealVar *firstPOI = dynamic_cast<RooRealVar *>(mc->GetParametersOfInterest()->first());

   firstPOI->setVal(muVal);
   //  firstPOI->setConstant();
   if (doFit) {
      mc->GetPdf()->fitTo(*data);
   }

   // -------------------------------------------------------

   mc->GetNuisanceParameters()->Print("v");
   int nPlotsMax = 1000;
   cout << " check expectedData by category" << endl;
   RooDataSet *simData = NULL;
   RooSimultaneous *simPdf = NULL;
   if (strcmp(mc->GetPdf()->ClassName(), "RooSimultaneous") == 0) {
      cout << "Is a simultaneous PDF" << endl;
      simPdf = (RooSimultaneous *)(mc->GetPdf());
   } else {
      cout << "Is not a simultaneous PDF" << endl;
   }

   if (doFit) {
      RooCategory *channelCat = (RooCategory *)(&simPdf->indexCat());
      TIterator *iter = channelCat->typeIterator();
      RooCatType *tt = NULL;
      tt = (RooCatType *)iter->Next();
      RooAbsPdf *pdftmp = ((RooSimultaneous *)mc->GetPdf())->getPdf(tt->GetName());
      RooArgSet *obstmp = pdftmp->getObservables(*mc->GetObservables());
      obs = ((RooRealVar *)obstmp->first());
      RooPlot *frame = obs->frame();
      cout << Form("%s==%s::%s", channelCat->GetName(), channelCat->GetName(), tt->GetName()) << endl;
      cout << tt->GetName() << " " << channelCat->getLabel() << endl;
      data->plotOn(frame, MarkerSize(1),
                   Cut(Form("%s==%s::%s", channelCat->GetName(), channelCat->GetName(), tt->GetName())),
                   DataError(RooAbsData::None));

      Double_t normCount =
         data->sumEntries(Form("%s==%s::%s", channelCat->GetName(), channelCat->GetName(), tt->GetName()));

      pdftmp->plotOn(frame, LineWidth(2.), Normalization(normCount, RooAbsReal::NumEvent));
      frame->Draw();
      cout << "expected events = " << mc->GetPdf()->expectedEvents(*data->get()) << endl;
      return;
   }

   int nPlots = 0;
   if (!simPdf) {

      TIterator *it = mc->GetNuisanceParameters()->createIterator();
      RooRealVar *var = NULL;
      while ((var = (RooRealVar *)it->Next()) != NULL) {
         RooPlot *frame = obs->frame();
         frame->SetYTitle(var->GetName());
         data->plotOn(frame, MarkerSize(1));
         var->setVal(0);
         mc->GetPdf()->plotOn(frame, LineWidth(1.));
         var->setVal(1);
         mc->GetPdf()->plotOn(frame, LineColor(kRed), LineStyle(kDashed), LineWidth(1));
         var->setVal(-1);
         mc->GetPdf()->plotOn(frame, LineColor(kGreen), LineStyle(kDashed), LineWidth(1));
         list->Add(frame);
         var->setVal(0);
      }

   } else {
      RooCategory *channelCat = (RooCategory *)(&simPdf->indexCat());
      //    TIterator* iter = simPdf->indexCat().typeIterator() ;
      TIterator *iter = channelCat->typeIterator();
      RooCatType *tt = NULL;
      while (nPlots < nPlotsMax && (tt = (RooCatType *)iter->Next())) {

         cout << "on type " << tt->GetName() << " " << endl;
         // Get pdf associated with state from simpdf
         RooAbsPdf *pdftmp = simPdf->getPdf(tt->GetName());

         // Generate observables defined by the pdf associated with this state
         RooArgSet *obstmp = pdftmp->getObservables(*mc->GetObservables());
         //      obstmp->Print();

         obs = ((RooRealVar *)obstmp->first());

         TIterator *it = mc->GetNuisanceParameters()->createIterator();
         RooRealVar *var = NULL;
         while (nPlots < nPlotsMax && (var = (RooRealVar *)it->Next())) {
            TCanvas *c2 = new TCanvas("c2");
            RooPlot *frame = obs->frame();
            frame->SetName(Form("frame%d", nPlots));
            frame->SetYTitle(var->GetName());

            cout << Form("%s==%s::%s", channelCat->GetName(), channelCat->GetName(), tt->GetName()) << endl;
            cout << tt->GetName() << " " << channelCat->getLabel() << endl;
            data->plotOn(frame, MarkerSize(1),
                         Cut(Form("%s==%s::%s", channelCat->GetName(), channelCat->GetName(), tt->GetName())),
                         DataError(RooAbsData::None));

            Double_t normCount =
               data->sumEntries(Form("%s==%s::%s", channelCat->GetName(), channelCat->GetName(), tt->GetName()));

            if (strcmp(var->GetName(), "Lumi") == 0) {
               cout << "working on lumi" << endl;
               var->setVal(w->var("nominalLumi")->getVal());
               var->Print();
            } else {
               var->setVal(0);
            }
            // w->allVars().Print("v");
            // mc->GetNuisanceParameters()->Print("v");
            // pdftmp->plotOn(frame,LineWidth(2.));
            // mc->GetPdf()->plotOn(frame,LineWidth(2.),Slice(*channelCat,tt->GetName()),ProjWData(*data));
            // pdftmp->plotOn(frame,LineWidth(2.),Slice(*channelCat,tt->GetName()),ProjWData(*data));
            normCount = pdftmp->expectedEvents(*obs);
            pdftmp->plotOn(frame, LineWidth(2.), Normalization(normCount, RooAbsReal::NumEvent));

            if (strcmp(var->GetName(), "Lumi") == 0) {
               cout << "working on lumi" << endl;
               var->setVal(w->var("nominalLumi")->getVal() + 0.05);
               var->Print();
            } else {
               var->setVal(nSigmaToVary);
            }
            // pdftmp->plotOn(frame,LineColor(kRed),LineStyle(kDashed),LineWidth(2));
            // mc->GetPdf()->plotOn(frame,LineColor(kRed),LineStyle(kDashed),LineWidth(2.),Slice(*channelCat,tt->GetName()),ProjWData(*data));
            // pdftmp->plotOn(frame,LineColor(kRed),LineStyle(kDashed),LineWidth(2.),Slice(*channelCat,tt->GetName()),ProjWData(*data));
            normCount = pdftmp->expectedEvents(*obs);
            pdftmp->plotOn(frame, LineWidth(2.), LineColor(kRed), LineStyle(kDashed),
                           Normalization(normCount, RooAbsReal::NumEvent));

            if (strcmp(var->GetName(), "Lumi") == 0) {
               cout << "working on lumi" << endl;
               var->setVal(w->var("nominalLumi")->getVal() - 0.05);
               var->Print();
            } else {
               var->setVal(-nSigmaToVary);
            }
            // pdftmp->plotOn(frame,LineColor(kGreen),LineStyle(kDashed),LineWidth(2));
            // mc->GetPdf()->plotOn(frame,LineColor(kGreen),LineStyle(kDashed),LineWidth(2),Slice(*channelCat,tt->GetName()),ProjWData(*data));
            // pdftmp->plotOn(frame,LineColor(kGreen),LineStyle(kDashed),LineWidth(2),Slice(*channelCat,tt->GetName()),ProjWData(*data));
            normCount = pdftmp->expectedEvents(*obs);
            pdftmp->plotOn(frame, LineWidth(2.), LineColor(kGreen), LineStyle(kDashed),
                           Normalization(normCount, RooAbsReal::NumEvent));

            // set them back to normal
            if (strcmp(var->GetName(), "Lumi") == 0) {
               cout << "working on lumi" << endl;
               var->setVal(w->var("nominalLumi")->getVal());
               var->Print();
            } else {
               var->setVal(0);
            }

            list->Add(frame);

            // quit making plots
            ++nPlots;

            frame->Draw();
            c2->SaveAs(Form("%s_%s_%s.pdf", tt->GetName(), obs->GetName(), var->GetName()));
            delete c2;
         }
      }
   }

   // -------------------------------------------------------

   // now make plots
   TCanvas *c1 = new TCanvas("c1", "ProfileInspectorDemo", 800, 200);
   if (list->GetSize() > 4) {
      double n = list->GetSize();
      int nx = (int)sqrt(n);
      int ny = TMath::CeilNint(n / nx);
      nx = TMath::CeilNint(sqrt(n));
      c1->Divide(ny, nx);
   } else
      c1->Divide(list->GetSize());
   for (int i = 0; i < list->GetSize(); ++i) {
      c1->cd(i + 1);
      list->At(i)->Draw();
   }
}
