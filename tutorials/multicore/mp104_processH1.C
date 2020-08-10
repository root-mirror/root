/// \file
/// \ingroup tutorial_multicore
/// \notebook -nodraw
/// Illustrate the usage of the multiproc to process the H1 analysis
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

static std::string tutname = "mp104_processH1: ";
static std::string logfile = "mp104_processH1.log";
static RedirectHandle_t gRH;

std::vector<std::string> files {"http://root.cern/files/h1/dstarmb.root",
                                "http://root.cern/files/h1/dstarp1a.root",
                                "http://root.cern/files/h1/dstarp1b.root",
                                "http://root.cern/files/h1/dstarp2.root"};

int mp104_processH1()
{

   // MacOSX may generate connection to WindowServer errors
   gROOT->SetBatch(kTRUE);

   TStopwatch stp;

   // Check and fit lambdas
   #include "mp_H1_lambdas.C"

   ROOT::TTreeProcessorMP pool(3);

   std::cout << tutname << "processing the H1 dataset with a lambda \n";

   auto hListFun = pool.Process(files, doH1, "h42");

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
   std::cout << tutname << "processing the H1 dataset with selector '" << selectorPath << "'\n";
   auto sel = TSelector::GetSelector(selectorPath);

   // In a second run we use sel
   gSystem->RedirectOutput(logfile.c_str(), "w", &gRH);
   auto hListSel = pool.Process(files, *sel, "h42");
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
