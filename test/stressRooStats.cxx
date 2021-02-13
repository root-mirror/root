// @(#)root/roofitcore:$name:  $:$id$
// Authors: Wouter Verkerke  November 2007

// C/C++ headers
#include <string>
#include <list>
#include <iostream>
#include <iomanip>
#include <cmath>

// ROOT headers
#include "TWebFile.h"
#include "TSystem.h"
#include "TString.h"
#include "TStopwatch.h"
#include "TROOT.h"
#include "TLine.h"
#include "TFile.h"
#include "TClass.h"
#include "TF1.h"
#include "TBenchmark.h"

// RooFit headers
#include "RooGlobalFunc.h"
#include "RooNumIntConfig.h"
#include "RooMsgService.h"
#include "RooResolutionModel.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooDouble.h"
#include "RooWorkspace.h"
#include "Roo1DTable.h"
#include "RooCurve.h"
#include "RooHist.h"
#include "RooRandom.h"
#include "RooTrace.h"
#include "RooMath.h"

// Tests file
#include "stressRooStats_tests.h"

using namespace std ;
using namespace RooFit ;


//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*//
//                                                                           //
// RooStats Unit Test S.T.R.E.S.S. Suite                                     //
// Authors: Ioan Gabriel Bucur, Lorenzo Moneta, Wouter Verkerke              //
//                                                                           //
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*//


////////////////////////////////////////////////////////////////////////////////
/// Print test program number and its title

void StatusPrint(const Int_t id, const TString &title, const Int_t status, const Int_t lineWidth)
{
   TString header = TString::Format("Test %d : %s ", id, title.Data());
   cout << left << setw(lineWidth) << setfill('.') << header << " " << (status > 0 ? "OK" : (status < 0 ? "SKIPPED" : "FAILED")) << endl;
}

////////////////////////////////////////////////////////////////////////////////
/// width of lines when printing test results

Int_t stressRooStats(const char* refFile, Bool_t writeRef, Int_t verbose, Bool_t allTests, Bool_t oneTest, Int_t testNumber, Bool_t dryRun, Bool_t doDump, Bool_t doTreeStore)
{
   const Int_t lineWidth = 120;

   // Save memory directory location
   auto memDir = gDirectory;
   RooUnitTest::setMemDir(gDirectory) ;

   if (doTreeStore) {
      RooAbsData::setDefaultStorageType(RooAbsData::Tree) ;
   }

   TFile* fref = 0 ;
   if (!dryRun) {
      if (TString(refFile).Contains("http:")) {
         if (writeRef) {
            cout << "stressRooStats ERROR: reference file must be local file in writing mode" << endl ;
            return -1 ;
         }
         fref = new TWebFile(refFile) ;
      } else {
         fref = new TFile(refFile, writeRef ? "RECREATE" : "") ;
      }
      if (fref->IsZombie()) {
         cout << "stressRooStats ERROR: cannot open reference file " << refFile << endl ;
         return -1 ;
      }
   }

   if (dryRun) {
      // Preload singletons here so they don't show up in trace accounting
      RooNumIntConfig::defaultConfig() ;
      RooResolutionModel::identity() ;

      RooTrace::active(1) ;
   }

   // Add dedicated logging stream for errors that will remain active in silent mode
   RooMsgService::instance().addStream(RooFit::ERROR) ;

   cout << left << setw(lineWidth) << setfill('*') << "" << endl;
   cout << "*" << setw(lineWidth - 2) << setfill(' ') << " RooStats S.T.R.E.S.S. suite " << "*" << endl;
   cout << setw(lineWidth) << setfill('*') << "" << endl;
   cout << setw(lineWidth) << setfill('*') << "" << endl;


   TStopwatch timer;
   timer.Start();

   list<RooUnitTest*> testList;


   // 1-5 TEST PLC CONFINT SIMPLE GAUSSIAN : Confidence Level range is (0,1)
   testList.push_back(new TestProfileLikelihoodCalculator1(fref, writeRef, verbose, 0.99999)); // boundary case CL -> 1
   testList.push_back(new TestProfileLikelihoodCalculator1(fref, writeRef, verbose, 2 * ROOT::Math::normal_cdf(3) - 1)); // 3 sigma
   testList.push_back(new TestProfileLikelihoodCalculator1(fref, writeRef, verbose, 2 * ROOT::Math::normal_cdf(2) - 1)); // 2 sigma
   testList.push_back(new TestProfileLikelihoodCalculator1(fref, writeRef, verbose, 2 * ROOT::Math::normal_cdf(1) - 1)); // 1 sigma
   testList.push_back(new TestProfileLikelihoodCalculator1(fref, writeRef, verbose, 0.00001)); // boundary case CL -> 0

   // 6-10 TEST PLC CONFINT SIMPLE POISSON : Observed value range is [0,1000]
   testList.push_back(new TestProfileLikelihoodCalculator2(fref, writeRef, verbose, 0)); // boundary Poisson value (0)
   testList.push_back(new TestProfileLikelihoodCalculator2(fref, writeRef, verbose, 1));
   testList.push_back(new TestProfileLikelihoodCalculator2(fref, writeRef, verbose, 5));
   testList.push_back(new TestProfileLikelihoodCalculator2(fref, writeRef, verbose, 100));
   testList.push_back(new TestProfileLikelihoodCalculator2(fref, writeRef, verbose, 800)); // boundary Poisson value

   // 11-13 TEST PLC CONFINT PRODUCT POISSON : Observed value range is [0,30] for x=s+b and [0,80] for y=2*s*1.2^beta
   testList.push_back(new TestProfileLikelihoodCalculator3(fref, writeRef, verbose, 10, 30));
   testList.push_back(new TestProfileLikelihoodCalculator3(fref, writeRef, verbose, 20, 25));
   testList.push_back(new TestProfileLikelihoodCalculator3(fref, writeRef, verbose, 15, 20, 2 * ROOT::Math::normal_cdf(2) - 1));

   // 14 TEST PLC HYPOTEST ON/OFF MODEL
   testList.push_back(new TestProfileLikelihoodCalculator4(fref, writeRef, verbose));

   // 15-18 TEST BC CONFINT CENTRAL SIMPLE POISSON : Observed value range is [0,100]
   testList.push_back(new TestBayesianCalculator1(fref, writeRef, verbose, 1));
   testList.push_back(new TestBayesianCalculator1(fref, writeRef, verbose, 3));
   testList.push_back(new TestBayesianCalculator1(fref, writeRef, verbose, 10));
   testList.push_back(new TestBayesianCalculator1(fref, writeRef, verbose, 50));

   // 19 TEST BC CONFINT SHORTEST SIMPLE POISSON
   testList.push_back(new TestBayesianCalculator2(fref, writeRef, verbose));

   // 20-22 TEST BC CONFINT CENTRAL PRODUCT POISSON : Observed value range is [0,30] for x=s+b and [0,80] for y=2*s*1.2^beta
   testList.push_back(new TestBayesianCalculator3(fref, writeRef, verbose, 10, 30));
   testList.push_back(new TestBayesianCalculator3(fref, writeRef, verbose, 20, 25));
   testList.push_back(new TestBayesianCalculator3(fref, writeRef, verbose, 15, 20, 2 * ROOT::Math::normal_cdf(2) - 1));

   // 23-25 TEST MCMCC CONFINT PRODUCT POISSON : Observed value range is [0,30] for x=s+b and [0,80] for y=2*s*1.2^beta
   testList.push_back(new TestMCMCCalculator(fref, writeRef, verbose, 10, 30));
   testList.push_back(new TestMCMCCalculator(fref, writeRef, verbose, 20, 25));
   testList.push_back(new TestMCMCCalculator(fref, writeRef, verbose, 15, 20, 2 * ROOT::Math::normal_cdf(2) - 1));

   // 26 TEST ZBI SIGNIFICANCE
   testList.push_back(new TestZBi(fref, writeRef, verbose));

   // 27-31 TEST PLC VS AC SIGNIFICANCE : Observed value range is [0,300] for on source and [0,1100] for off-source; tau has the range [0.1,5.0]
   testList.push_back(new TestHypoTestCalculator1(fref, writeRef, verbose, 150, 100, 1.0));
   testList.push_back(new TestHypoTestCalculator1(fref, writeRef, verbose, 200, 100, 1.0));
   testList.push_back(new TestHypoTestCalculator1(fref, writeRef, verbose, 105, 100, 1.0));
   testList.push_back(new TestHypoTestCalculator1(fref, writeRef, verbose, 150, 10, 0.1));
   testList.push_back(new TestHypoTestCalculator1(fref, writeRef, verbose, 150, 400, 4.0));

   // 32-36 TEST HTC SIGNIFICANCE
   testList.push_back(new TestHypoTestCalculator2(fref, writeRef, verbose, kAsymptotic));
   testList.push_back(new TestHypoTestCalculator2(fref, writeRef, verbose, kFrequentist, kSimpleLR));
   testList.push_back(new TestHypoTestCalculator2(fref, writeRef, verbose, kFrequentist, kRatioLR));
   testList.push_back(new TestHypoTestCalculator2(fref, writeRef, verbose, kFrequentist, kProfileLROneSidedDiscovery));
   testList.push_back(new TestHypoTestCalculator2(fref, writeRef, verbose, kHybrid, kProfileLROneSidedDiscovery));

   // 37-43 TEST HTI PRODUCT POISSON : Observed value range is [0,30] for x=s+b and [0,80] for y=2*s*1.2^beta
   testList.push_back(new TestHypoTestInverter1(fref, writeRef, verbose, kAsymptotic, kProfileLR, 10, 30));
   testList.push_back(new TestHypoTestInverter1(fref, writeRef, verbose, kAsymptotic, kProfileLR, 20, 25));
   testList.push_back(new TestHypoTestInverter1(fref, writeRef, verbose, kAsymptotic, kProfileLR, 15, 20));
   testList.push_back(new TestHypoTestInverter1(fref, writeRef, verbose, kFrequentist, kProfileLR, 10, 30));
   testList.push_back(new TestHypoTestInverter1(fref, writeRef, verbose, kFrequentist, kProfileLR, 20, 25));
   testList.push_back(new TestHypoTestInverter1(fref, writeRef, verbose, kFrequentist, kProfileLR, 15, 20));
   testList.push_back(new TestHypoTestInverter1(fref, writeRef, verbose, kHybrid, kProfileLR, 10, 30));

   // 44-48 TEST HTI S+B+E POISSON : Observed value range is [0,50] for x = e*s+b
   testList.push_back(new TestHypoTestInverter2(fref, writeRef, verbose, kAsymptotic, kProfileLROneSided, 10, 0.95));
   testList.push_back(new TestHypoTestInverter2(fref, writeRef, verbose, kAsymptotic, kProfileLROneSided, 20));
//    testList.push_back(new TestHypoTestInverter2(fref, writeRef, verbose, kFrequentist, kSimpleLR, 10));
//    testList.push_back(new TestHypoTestInverter2(fref, writeRef, verbose, kFrequentist, kSimpleLR, 20));
   testList.push_back(new TestHypoTestInverter2(fref, writeRef, verbose, kFrequentist, kRatioLR, 10, 0.95));
   testList.push_back(new TestHypoTestInverter2(fref, writeRef, verbose, kFrequentist, kProfileLROneSided, 10, 0.95));
   testList.push_back(new TestHypoTestInverter2(fref, writeRef, verbose, kHybrid, kSimpleLR, 10, 0.95));


   TString suiteType = TString::Format(" Starting S.T.R.E.S.S. %s",
                                       allTests ? "full suite" : (oneTest ? TString::Format("test %d", testNumber).Data() : "basic suite")
                                      );

   cout << "*" << setw(lineWidth - 3) << setfill(' ') << suiteType << " *" << endl;
   cout << setw(lineWidth) << setfill('*') << "" << endl;

   if (doDump) {
      TFile fdbg("stressRooStats_DEBUG.root", "RECREATE") ;
   }

   gBenchmark->Start("stressRooStats");

   int nFailed = 0;
   {
      Int_t i;
      list<RooUnitTest*>::iterator iter;

      if (oneTest && (testNumber <= 0 || (UInt_t) testNumber > testList.size())) {
         cout << "Tests are numbered from 1 to " << testList.size() << endl;
      } else {
         for (iter = testList.begin(), i = 1; iter != testList.end(); iter++, i++) {
            if (!oneTest || testNumber == i) {
               if (doDump) {
                  (*iter)->setDebug(kTRUE);
               }
               int status =  (*iter)->isTestAvailable() ? (*iter)->runTest() : -1;
               StatusPrint(i, (*iter)->GetName(), status , lineWidth);
               if (!status) nFailed++;  // do not count the skipped tests
            }
            delete *iter;
         }
      }
   }

   if (dryRun) {
      RooTrace::dump();
   }

   gBenchmark->Stop("stressRooStats");


   //Print table with results
   Bool_t UNIX = strcmp(gSystem->GetName(), "Unix") == 0;
   cout << setw(lineWidth) << setfill('*') << "" << endl;
   if (UNIX) {
      TString sp = gSystem->GetFromPipe("uname -a");
      cout << "* SYS: " << sp << endl;
      if (strstr(gSystem->GetBuildNode(), "Linux")) {
         sp = gSystem->GetFromPipe("lsb_release -d -s");
         cout << "* SYS: " << sp << endl;
      }
      if (strstr(gSystem->GetBuildNode(), "Darwin")) {
         sp  = gSystem->GetFromPipe("sw_vers -productVersion");
         sp += " Mac OS X ";
         cout << "* SYS: " << sp << endl;
      }
   } else {
      const Char_t *os = gSystem->Getenv("OS");
      if (!os) cout << "*  SYS: Windows 95" << endl;
      else     cout << "*  SYS: " << os << " " << gSystem->Getenv("PROCESSOR_IDENTIFIER") << endl;
   }

   cout << setw(lineWidth) << setfill('*') << "" << endl;
   gBenchmark->Print("stressRooStats");
#ifdef __CINT__
   Double_t reftime = 186.34; //pcbrun4 interpreted
#else
   Double_t reftime = 93.59; //pcbrun4 compiled
#endif
   const Double_t rootmarks = 860 * reftime / gBenchmark->GetCpuTime("stressRooStats");

   cout << setw(lineWidth) << setfill('*') << "" << endl;
   cout << TString::Format("*  ROOTMARKS = %6.1f  *  Root %-8s %d/%d", rootmarks, gROOT->GetVersion(),
                           gROOT->GetVersionDate(), gROOT->GetVersionTime()) << endl;
   cout << setw(lineWidth) << setfill('*') << "" << endl;

   // NOTE: The function TStopwatch::CpuTime() calls Tstopwatch::Stop(), so you do not need to stop the timer separately.
   cout << "Time at the end of job = " << timer.CpuTime() << " seconds" << endl;

   if (fref) {
      fref->Close() ;
      delete fref ;
   }

   delete gBenchmark ;
   gBenchmark = 0 ;

   // Some of the object are multiple times in the list, let's make sure they
   // are not deleted twice.
   // The addition of memDir to the list of Cleanups is not needed if it already
   // there, for example if memDir is gROOT.
   bool needCleanupAdd = nullptr == gROOT->GetListOfCleanups()->FindObject(memDir->GetList());
   if (needCleanupAdd) gROOT->GetListOfCleanups()->Add(memDir->GetList());

   memDir->GetList()->Delete("slow");

   if (needCleanupAdd) gROOT->GetListOfCleanups()->Remove(memDir->GetList());

   return nFailed;
}

//_____________________________batch only_____________________
#ifndef __CINT__

int main(int argc, const char *argv[])
{
   Bool_t doWrite     = kFALSE;
   Int_t  verbose     =      0;
   Bool_t allTests    = kFALSE;
   Bool_t oneTest     = kFALSE;
   Int_t testNumber   =      0;
   Bool_t dryRun      = kFALSE;
   Bool_t doDump      = kFALSE;
   Bool_t doTreeStore = kFALSE;

   //string refFileName = "http://root.cern.ch/files/stressRooStats_v534_ref.root" ;
   string refFileName = "$ROOTSYS/test/stressRooStats_ref.root" ;
   string minimizerName = "Minuit";

   // Parse command line arguments
   for (Int_t i = 1 ;  i < argc ; i++) {
      string arg = argv[i] ;

      if (arg == "-f") {
         cout << "stressRooStats: using reference file " << argv[i + 1] << endl ;
         refFileName = argv[++i] ;
      } else if (arg == "-w") {
         cout << "stressRooStats: running in writing mode to update reference file" << endl ;
         doWrite = kTRUE ;
      } else if (arg == "-mc") {
         cout << "stressRooStats: running in memcheck mode, no regression tests are performed" << endl;
         dryRun = kTRUE;
      } else if (arg == "-min" || arg == "-minim") {
         cout << "stressRooStats: running using minimizer " << argv[i +1]  << endl;
         minimizerName = argv[++i] ;
      } else if (arg == "-ts") {
         cout << "stressRooStats: setting tree-based storage for datasets" << endl;
         doTreeStore = kTRUE;
      } else if (arg == "-v") {
         cout << "stressRooStats: running in verbose mode" << endl;
         verbose = 1;
      } else if (arg == "-vv") {
         cout << "stressRooStats: running in very verbose mode" << endl;
         verbose = 2;
      } else if (arg == "-vvv") {
         cout << "stressRooStats: running in very very verbose mode" << endl;
         verbose = 3;
      } else if (arg == "-a") {
         cout << "stressRooStats: deploying full suite of tests" << endl;
         allTests = kTRUE;
      } else if (arg == "-n") {
         cout << "stressRooStats: running single test" << endl;
         oneTest = kTRUE;
         testNumber = atoi(argv[++i]);
      } else if (arg == "-d") {
         cout << "stressRooStats: setting gDebug to " << argv[i + 1] << endl;
         gDebug = atoi(argv[++i]);
      } else if (arg == "-c") {
         cout << "stressRooStats: dumping comparison file for failed tests " << endl;
         doDump = kTRUE;
      } else if (arg == "-h") {
         cout << "usage: stressRooStats [ options ] " << endl;
         cout << "" << endl;
         cout << "       -f <file> : use given reference file instead of default (" << refFileName << ")" << endl;
         cout << "       -w        : write reference file, instead of reading file and running comparison tests" << endl;
         cout << "       -n N      : only run test with sequential number N" << endl;
         cout << "       -a        : run full suite of tests (default is basic suite); this overrides the -n single test option" << endl;
         cout << "       -c        : dump file stressRooStats_DEBUG.root to which results of both current result and reference for each failed test are written" << endl;
         cout << "       -mc       : memory check mode, no regression test are performed. Set this flag when running with valgrind" << endl;
         cout << "     -min <name> : minimizer name (default is Minuit2_ if available" << endl;
         cout << "       -vs       : use vector-based storage for all datasets (default is tree-based storage)" << endl;
         cout << "       -v/-vv    : set verbose mode (show result of each regression test) or very verbose mode (show all roofit output as well)" << endl;
         cout << "       -d N      : set ROOT gDebug flag to N" << endl ;
         cout << " " << endl ;
         return 0 ;
      }

   }

//    if (doWrite && refFileName.find("http:") == 0) {

//       // Locate file name part in URL and update refFileName accordingly
//       char* buf = new char[refFileName.size() + 1];
//       strcpy(buf, refFileName.c_str());
//       char *ptr = strrchr(buf, '/');
//       if (!ptr) ptr = strrchr(buf, ':');
//       refFileName = ptr + 1;
//       delete[] buf;

//       cout << "stressRooStats: WARNING running in write mode, but reference file is web file, writing local file instead: "
//            << refFileName << endl;
//    }

   // set minimizer
    // use Minut2 if available
   std::string minimizerType = minimizerName;
   // check in case of Minuit2 and set Minuit in case we cannot use it
   if (minimizerType == "Minuit2") {
      int prec = gErrorIgnoreLevel;
      gErrorIgnoreLevel = kFatal;
      if (gSystem->Load("libMinuit2") < 0) minimizerType = "Minuit";
      gErrorIgnoreLevel=prec;
   }
   ROOT::Math::MinimizerOptions::SetDefaultMinimizer(minimizerType.c_str());


   gBenchmark = new TBenchmark();
   return stressRooStats(refFileName.c_str(), doWrite, verbose, allTests, oneTest, testNumber, dryRun, doDump, doTreeStore);
}

////////////////////////////////////////////////////////////////////////////////

Int_t stressRooStats()
{
   Bool_t doWrite     = kFALSE;
   Int_t  verbose     =      0;
   Bool_t allTests    = kFALSE;
   Bool_t oneTest     = kFALSE;
   Int_t testNumber   =      0;
   Bool_t dryRun      = kFALSE;
   Bool_t doDump      = kFALSE;
   Bool_t doTreeStore = kFALSE;
   string refFileName = "$ROOTSYS/test/stressRooStats_ref.root";
   return stressRooStats(refFileName.c_str(), doWrite, verbose, allTests, oneTest, testNumber, dryRun, doDump, doTreeStore);
}

#endif
