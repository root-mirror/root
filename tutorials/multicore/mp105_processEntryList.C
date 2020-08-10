/// \file
/// \ingroup tutorial_multicore
/// \notebook -nodraw
/// Illustrate the usage of the multiproc to process TEntryList with the H1 analysis
/// example.
///
/// \macro_code
///
/// \author Gerardo Ganis

#include "TString.h"
#include "TROOT.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TEntryList.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
#include "TSystem.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TF1.h"
#include "TLine.h"
#include "TPaveStats.h"
#include "TStopwatch.h"
#include "ROOT/TTreeProcessorMP.hxx"

static std::string tutname = "mp105_processEntryList: ";
static std::string logfile = "mp105_processEntryList.log";
static RedirectHandle_t gRH;

std::vector<std::string> files = {"http://root.cern/files/h1/dstarmb.root",
                                  "http://root.cern/files/h1/dstarp1a.root",
                                  "http://root.cern/files/h1/dstarp1b.root",
                                  "http://root.cern/files/h1/dstarp2.root"};

int mp105_processEntryList()
{

   // MacOSX may generate connection to WindowServer errors
   gROOT->SetBatch(kTRUE);

   TStopwatch stp;

#include "mp_H1_lambdas.C"

   ROOT::TTreeProcessorMP pool(3);

   std::cout << tutname << "creating the entry list \n";

   auto sumElist = pool.Process(files, doH1fillList, "h42");

   // Print the entry list
   if (sumElist) {
      sumElist->Print();
   } else {
      std::cout << tutname << " ERROR creating the entry list \n";
      return -1;
   }

   // Time taken
   stp.Print();
   stp.Start();

   // Let's analyse H1 with the list
   std::cout << tutname << "processing the entry list with a lambda \n";

   // Run the analysis
   auto hListFun = pool.Process(files, doH1useList, *sumElist, "h42");

   // Check the output
   if (checkH1(hListFun) < 0)
      return -1;

   // Do the fit
   if (doFit(hListFun, logfile.c_str()) < 0)
      return -1;

   stp.Print();
   stp.Start();

   // Run the analysis with a selector
   TString selectorPath = gROOT->GetTutorialDir();
   selectorPath += "/tree/h1analysisTreeReader.C+";
   std::cout << tutname << "processing the entry list with selector '" << selectorPath << "'\n";
   auto sel = TSelector::GetSelector(selectorPath);

   // In a second run we use sel
   sel->SetOption("useList");
   gSystem->RedirectOutput(logfile.c_str(), "w", &gRH);
   auto hListSel = pool.Process(files, *sel, *sumElist, "h42");
   gSystem->RedirectOutput(0, 0, &gRH);

   // Check the output
   if (checkH1(hListSel) < 0)
      return -1;

   // Do the fit
   if (doFit(hListSel, logfile.c_str()) < 0)
      return -1;

   stp.Print();
   stp.Start();

   return 0;
}
