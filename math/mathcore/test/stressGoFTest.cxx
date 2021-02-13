#include "TBenchmark.h"
#include "Math/GoFTest.h"
#include "Math/Functor.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStopwatch.h"

#ifdef R__HAS_MATHMORE
#include "Math/GSLRndmEngines.h"
#include "Math/Random.h"
#endif

#include "TRandom3.h"

#include <iostream>

#include <cassert>

/*N.B.: The tests' expected values (expectedDn and expectedA2) were computed on Pcphsft54.cern.ch i386 GNU/Linux computer (slc4_ia32_gcc34)

  LM. (16/9/14)  Expected values for AD2 test have been computed with R kSamples package
*/
struct GoFTStress {

   static enum EDebugLevelTypes {
      kNoDebug,       // no test results printing
      kBasicDebug,    // prints either "OK" or "FAIL" test results
      kStandardDebug, // prints summarized test results
   } fgDebugLevel;

   Int_t RunTests() {
      Int_t result = 0;
      result += UnitTest1();
      result += UnitTest2();
      result += UnitTest3();
      result += UnitTest4();
      result += UnitTest5();
      result += UnitTest6();
      result += UnitTest7();
      return result;
   }

   void PrintBenchmark() {
      if (fgDebugLevel == kNoDebug) {
         return;
      }

      //Print table with results
      Bool_t UNIX = strcmp(gSystem->GetName(), "Unix") == 0;
      printf("******************************************************************\n");
      if (UNIX) {
         TString sp = gSystem->GetFromPipe("uname -a");
         sp.Resize(60);
         printf("*  SYS: %s\n",sp.Data());
         if (strstr(gSystem->GetBuildNode(),"Linux")) {
            sp = gSystem->GetFromPipe("lsb_release -d -s");
            printf("*  SYS: %s\n",sp.Data());
         }
         if (strstr(gSystem->GetBuildNode(),"Darwin")) {
            sp  = gSystem->GetFromPipe("sw_vers -productVersion");
            sp += " Mac OS X ";
            printf("*  SYS: %s\n",sp.Data());
         }
      } else {
         const Char_t *os = gSystem->Getenv("OS");
         if (!os) {
            printf("*  SYS: Windows 95\n");
         } else {
            printf("*  SYS: %s %s \n",os,gSystem->Getenv("PROCESSOR_IDENTIFIER"));
         }
      }

      printf("*****************************************************************************************\n");
      gBenchmark->Print("GoFTestStress");
#ifdef __CINT__
      Double_t reftime = 0.02 ; // slc4rabacal interpreted (CPU time taken to run complete tests with ACLiC)
#else
      Double_t reftime = 0.04; // slc4rabacal compiled (CPU time taken to run complete tests)
#endif

      Double_t rootmarks = 800. * reftime / gBenchmark->GetCpuTime("GoFTestStress");

      printf("*****************************************************************************************\n");
      printf("*  ROOTMARKS = %6.1f   *  Root%-8s  %d/%d\n", rootmarks, gROOT->GetVersion(),
            gROOT->GetVersionDate(),gROOT->GetVersionTime());
      printf("*****************************************************************************************\n");
   }

   private:

   Int_t UnitTest1() {
      std::cout << "UNIT TEST 1" << std::endl;

      UInt_t nsmps = 2;
      const UInt_t smpSize1 = 71;
      const UInt_t smpSize2 = 135;
   /*
      Data set adapted from the paper (1)
      "K-Sample Anderson-Darling Tests of Fit for continuous and discrete cases" by Scholz and Stephens
      values of expected A2 taken by running R kSamples code
   */
      const Double_t smp1[smpSize1] = {194, 15, 41, 29, 33, 181, 413, 14, 58, 37, 100, 65, 9, 169, 447, 18, 4, 36, 201, 118, 34, 31, 18, 18, 67, 57, 62, 7, 22, 34, 90, 10, 60, 186, 61, 49, 14, 24, 56, 20, 79, 84, 44, 59, 29, 118, 25, 156, 310, 76, 26, 44, 23, 62, 130, 208, 70, 101, 208, 74, 57, 48, 29, 502, 12, 70, 21, 29, 386, 59, 27};

      const Double_t smp2[smpSize2] = {55, 320, 56, 104, 220, 239, 47, 246, 176, 182, 33, 23, 261, 87, 7, 120, 14, 62, 47, 225, 71, 246, 21, 42, 20, 5, 12, 120, 11, 3, 14, 71, 11, 14, 11, 16, 90, 1, 16, 52, 95, 97, 51, 11, 4, 141, 18, 142, 68, 77, 80, 1, 16, 106, 206, 82, 54, 31, 216, 46, 111, 39, 63, 18, 191, 18, 163, 24, 50, 44, 102, 72, 22, 39, 3, 15, 197, 188, 79, 88, 46, 5, 5, 36, 22, 139, 210, 97, 30, 23, 13, 14, 359, 9, 12, 270, 603, 3, 104, 2, 438, 50, 254, 5, 283, 35, 12, 487, 18, 100, 7, 98, 5, 85, 91, 43, 230, 3, 130, 102, 209, 14, 57, 54, 32, 67, 59, 134, 152, 27, 14, 230, 66, 61, 34};

      ROOT::Math::GoFTest* goft = new ROOT::Math::GoFTest( smpSize1, smp1, smpSize2, smp2 );

      Double_t A2 = goft->AndersonDarling2SamplesTest("t"); // standartized A2_akN
   //     Double_t pvalueAD = goft->AndersonDarling2SamplesTest();
      Double_t pvalueAD = (*goft)(ROOT::Math::GoFTest::kAD2s);

      Double_t expectedA2_akN = 1.5686 ; // A2_akN in (1)

      Double_t sigmaN = 0.754539;    // sigmaN in (1)

      Double_t zScore = 0.75360;    // zScore in (1)

      Int_t result = PrintResultAD2Samples(nsmps, A2, expectedA2_akN, sigmaN, zScore, pvalueAD);

   //     Double_t Dn = goft->KolmogorovSmirnov2SamplesTest("t"); // standartized A2_akN
      Double_t Dn = (*goft)(ROOT::Math::GoFTest::kKS2s, "t");
      Double_t pvalueKS = goft->KolmogorovSmirnov2SamplesTest();

      Double_t expectedDn = 0.139176;

      result += PrintResultKS(nsmps, Dn, expectedDn, pvalueKS);

      delete goft;
      return result;
   }

   Int_t UnitTest2() {
      std::cout << "UNIT TEST 2" << std::endl;

      const UInt_t nsmps = 2;
      const UInt_t smpSize = 16;
   /*
      Data sets adapted from the paper (1)
      "K-Sample Anderson-Darling Tests of Fit for continuous and discrete cases" by Scholz and Stephens
   */
      const Double_t smp1[smpSize] = {38.7, 41.5, 43.8, 44.5, 45.5, 46.0, 47.7, 58.0, 39.2, 39.3, 39.7, 41.4, 41.8, 42.9, 43.3, 45.8};
      const Double_t smp2[smpSize] = {34.0, 35.0, 39.0, 40.0, 43.0, 43.0, 44.0, 45.0, 34.0, 34.8, 34.8, 35.4, 37.2, 37.8, 41.2, 42.8};

      ROOT::Math::GoFTest* goft = new ROOT::Math::GoFTest( smpSize, smp1, smpSize, smp2 );

   //     Double_t A2 = goft->AndersonDarling2SamplesTest("t"); // standartized A2_akN
      Double_t A2 = (*goft)(ROOT::Math::GoFTest::kAD2s, "t");
      Double_t pvalueAD = goft->AndersonDarling2SamplesTest();

      Double_t expectedA2_akN = 4.5516; // unstandardized A2_akN in (1) (values verified and obtained with R kSamples )

      Double_t sigmaN = 0.71939;   // sigmaN in (1)

      Double_t zScore = 4.9369;    // zScore in (1) (version 1)

      Int_t result = PrintResultAD2Samples(nsmps, A2, expectedA2_akN, sigmaN, zScore, pvalueAD);

      Double_t Dn = goft->KolmogorovSmirnov2SamplesTest("t"); // standartized A2_akN
   //     Double_t pvalueKS = goft->KolmogorovSmirnov2SamplesTest();
      Double_t pvalueKS = (*goft)(ROOT::Math::GoFTest::kKS2s);

      Double_t expectedDn = 0.5;

      result += PrintResultKS(nsmps, Dn, expectedDn, pvalueKS);

      delete goft;
      return result;
   }

   Int_t UnitTest3() {
      std::cout << "UNIT TEST 3" << std::endl;

#ifndef R__HAS_MATHMORE
      std::cout << "SKIPPED (Mathmore is not present) " << std::endl;
      return 0;
#else

      UInt_t nEvents = 1000;
      UInt_t nsmps = 1;

      ROOT::Math::Random<ROOT::Math::GSLRngMT> r;

      Double_t* sample = new Double_t[nEvents];

      for (UInt_t i = 0; i < nEvents; ++i) {
         Double_t data = r.LogNormal(5.0, 2.0);
         sample[i] = data;
         assert(sample[i] == data);
      }

      ROOT::Math::GoFTest* goft = new ROOT::Math::GoFTest(nEvents, sample, ROOT::Math::GoFTest::kLogNormal);

      if (GoFTStress::fgDebugLevel == GoFTStress::kStandardDebug)
         std::cout << "LogNormal fitting" << std::endl;

   //     Double_t A2 = goft->AndersonDarlingTest("t");
      Double_t A2 = goft->operator()(ROOT::Math::GoFTest::kAD, "t");
      Double_t pvalueAD = goft->AndersonDarlingTest();

      Double_t expectedA2 = 0.422771;

      Int_t result = PrintResultAD1Sample(A2, expectedA2, pvalueAD);

      Double_t Dn = goft->KolmogorovSmirnovTest("t");
   //     Double_t pvalueKS = goft->KolmogorovSmirnovTest();
      Double_t pvalueKS = goft->operator()(ROOT::Math::GoFTest::kKS);

      Double_t expectedDn = 0.0204916;

      result += PrintResultKS(nsmps, Dn, expectedDn, pvalueKS);

      delete goft;
      delete [] sample;

      return result;
#endif
   }

   Int_t UnitTest4() {
      std::cout << "UNIT TEST 4" << std::endl;

      UInt_t nEvents = 1000;
      UInt_t nsmps = 1;

      TRandom3 r;

      Double_t* sample = new Double_t[nEvents];

      for (UInt_t i = 0; i < nEvents; ++i) {
         Double_t data = r.Exp(1.54);
         sample[i] = data;
         assert(sample[i] == data);
      }

      ROOT::Math::GoFTest* goft = new ROOT::Math::GoFTest(nEvents, sample, ROOT::Math::GoFTest::kExponential);

      if (GoFTStress::fgDebugLevel == GoFTStress::kStandardDebug)
         std::cout << "**Exponential fitting**" << std::endl;

      Double_t A2 = goft->AndersonDarlingTest("t");
   //     Double_t pvalueAD = goft->AndersonDarlingTest();
      Double_t pvalueAD = goft->operator()();

      Double_t expectedA2 = 0.521153;

      Int_t result = PrintResultAD1Sample(A2, expectedA2, pvalueAD);

   //     Double_t Dn = goft->KolmogorovSmirnovTest("t");
      Double_t Dn = goft->operator()(ROOT::Math::GoFTest::kKS, "t");
      Double_t pvalueKS = goft->KolmogorovSmirnovTest();

      Double_t expectedDn = 0.0218148;

      result += PrintResultKS(nsmps, Dn, expectedDn, pvalueKS);

      delete goft;
      delete [] sample;

      return result;
   }

   Int_t UnitTest5() {
      std::cout << "UNIT TEST 5" << std::endl;

      UInt_t nEvents = 1000;
      UInt_t nsmps = 1;

      TRandom3 r;

      Double_t* sample = new Double_t[nEvents];

      for (UInt_t i = 0; i < nEvents; ++i) {
         Double_t data = r.Gaus(300, 50);
         sample[i] = data;
         assert(sample[i] == data);
      }

      ROOT::Math::GoFTest* goft = new ROOT::Math::GoFTest(nEvents, sample);
      goft->SetDistribution(ROOT::Math::GoFTest::kGaussian);

      if (GoFTStress::fgDebugLevel == GoFTStress::kStandardDebug)
         std::cout << "**Gaussian fitting**" << std::endl;

      Double_t A2 =  goft->AndersonDarlingTest("t");
      Double_t pvalueAD = goft->AndersonDarlingTest();

      Double_t expectedA2 = 0.441755;

      Int_t result = PrintResultAD1Sample(A2, expectedA2, pvalueAD);

      Double_t Dn = goft->KolmogorovSmirnovTest("t");
      Double_t pvalueKS = goft->KolmogorovSmirnovTest();

      Double_t expectedDn = 0.0282508;

      result += PrintResultKS(nsmps, Dn, expectedDn, pvalueKS);

      delete goft;

      delete [] sample;

      return result;
   }

   Int_t UnitTest6() {
      std::cout << "UNIT TEST 6" << std::endl;
      UInt_t nEvents = 1000;
      UInt_t nsmps = 1;

      TRandom3 r;

      Double_t* sample = new Double_t[nEvents];

      for (UInt_t i = 0; i < nEvents; ++i) {
         Double_t data = r.Landau();
         sample[i] = data;
         assert(sample[i] == data);
      }

      ROOT::Math::Functor1D userCdf(&TMath::LandauI);
      ROOT::Math::GoFTest* goft = new ROOT::Math::GoFTest(nEvents, sample, userCdf, ROOT::Math::GoFTest::kCDF);

      if (GoFTStress::fgDebugLevel == GoFTStress::kStandardDebug)
         std::cout << "**Landau fitting**" << std::endl;

      Double_t A2 =  goft->AndersonDarlingTest("t");
      Double_t pvalueAD = goft->AndersonDarlingTest();

      Double_t expectedA2 = 0.544658;

      Int_t result = PrintResultAD1Sample(A2, expectedA2, pvalueAD);

      Double_t Dn = goft->KolmogorovSmirnovTest("t");
      Double_t pvalueKS = goft->KolmogorovSmirnovTest();

      Double_t expectedDn = 0.0203432;

      result += PrintResultKS(nsmps, Dn, expectedDn, pvalueKS);

      delete goft;
      delete [] sample;

      return result;
   }

   Int_t UnitTest7() {
      std::cout << "UNIT TEST 7" << std::endl;

      UInt_t nEvents = 1000;
      UInt_t nsmps = 1;

      TRandom3 r;

      Double_t* sample = new Double_t[nEvents];

      for (UInt_t i = 0; i < nEvents; ++i) {
         Double_t data = r.Landau();
         sample[i] = data;
      }

      // need to specify min and max otherwise pdf does not converge
      ROOT::Math::GoFTest* goft = new ROOT::Math::GoFTest(nEvents, sample);

      ROOT::Math::Functor1D userPdf([](double x){ return TMath::Landau(x);});
      // need to use a reasanble range for the Landau
      // but must be bigger than xmin and xmax
      double xmin = 3*TMath::MinElement(nEvents, sample);
      double xmax = 3*TMath::MaxElement(nEvents, sample);

      if (GoFTStress::fgDebugLevel == GoFTStress::kStandardDebug)
         std::cout << "**Landau fitting**" << " in [ " << xmin << " , " << xmax << " ]" << std::endl;

      goft->SetUserDistribution(userPdf,ROOT::Math::GoFTest::kPDF, xmin, xmax);

      Double_t A2 =  goft->AndersonDarlingTest("t");
      Double_t pvalueAD = goft->AndersonDarlingTest();

      Double_t expectedA2 =  0.544658;

      // use larger tolerance due to truncation error of Landau
      Int_t result = PrintResultAD1Sample(A2, expectedA2, pvalueAD,0.002);

      Double_t Dn = goft->KolmogorovSmirnovTest("t");
      Double_t pvalueKS = goft->KolmogorovSmirnovTest();

      Double_t expectedDn = 0.0203432;

      result += PrintResultKS(nsmps, Dn, expectedDn, pvalueKS,0.001);

      delete goft;
      delete [] sample;

      return result;
   }

   Int_t PrintResultAD2Samples(UInt_t nsmps, Double_t A2, Double_t expectedA2_akN, Double_t sigmaN, Double_t zScore, Double_t pvalue, Double_t limit = 0.0001) {

      Double_t A2_akN = A2  * sigmaN + (nsmps - 1);

      if (GoFTStress::fgDebugLevel == GoFTStress::kStandardDebug) {
         std::cout << "Anderson-Darling 2-Samples Test" << std::endl;

         std::cout << "pvalue = " << pvalue << std::endl;

         std::cout << "zScore = " << zScore << std::endl;

         std::cout << "Expected A2 value = " << A2 << std::endl;

         std::cout << "A2_akN = " << A2_akN << std::endl;

         std::cout << "Expected A2_akN value = " << expectedA2_akN << std::endl;
      }
      if (TMath::Abs(A2_akN - expectedA2_akN) < limit * expectedA2_akN || TMath::Abs(A2 - zScore) < limit * zScore) {
         if (GoFTStress::fgDebugLevel >= GoFTStress::kBasicDebug)
         std::cout << "+++++++ TEST OK +++++++ \n" << std::endl;
         return EXIT_SUCCESS;
      } else {
         if (GoFTStress::fgDebugLevel >= GoFTStress::kBasicDebug)
         std::cout << "------- TEST FAIL ------- \n" << std::endl;
         return EXIT_FAILURE;
      }
   }

   Int_t PrintResultAD1Sample(Double_t A2, Double_t expectedA2, Double_t pvalue, Double_t limit = 0.0001) {

      if (GoFTStress::fgDebugLevel == GoFTStress::kStandardDebug) {
         std::cout << "Anderson-Darling 1-Sample Test" << std::endl;

         std::cout << "pvalue = " << pvalue << std::endl;

         std::cout << "A2 value = " << A2 << std::endl;

         std::cout << "Expected A2 value = " << expectedA2 << std::endl;
      }
      if (TMath::Abs(A2 - expectedA2) < limit * expectedA2) {
         if (GoFTStress::fgDebugLevel >= GoFTStress::kBasicDebug)
         std::cout << "+++++++ TEST OK +++++++ \n" << std::endl;
         return EXIT_SUCCESS;
      } else {
         if (GoFTStress::fgDebugLevel >= GoFTStress::kBasicDebug)
         std::cout << "------- TEST FAIL ------- \n" << std::endl;
         return EXIT_FAILURE;
      }
   }

   Int_t PrintResultKS(UInt_t nsmps, Double_t Dn, Double_t expectedDn, Double_t pvalue, Double_t limit = 0.0001) {

      if (GoFTStress::fgDebugLevel == GoFTStress::kStandardDebug) {
         std::cout << "Kolmogorov-Smirnov "<< nsmps << "-Samples Test" << std::endl;

         std::cout << "pvalue = " << pvalue << std::endl;

         std::cout << "Dn value = " << Dn << std::endl;

         std::cout << "Expected Dn value = " << expectedDn << std::endl;
      }
      if (TMath::Abs(Dn - expectedDn) < limit * expectedDn) {
         if (GoFTStress::fgDebugLevel >= GoFTStress::kBasicDebug)
         std::cout << "+++++++ TEST OK +++++++ \n" << std::endl;
         return EXIT_SUCCESS;
      } else {
         if (GoFTStress::fgDebugLevel >= GoFTStress::kBasicDebug)
         std::cout << "------- TEST FAIL ------- \n" << std::endl;
         return EXIT_FAILURE;
      }
   }
};

#ifdef __MAKECINT__
#pragma link C++ class GoFTStress-;
#endif

GoFTStress::EDebugLevelTypes GoFTStress::fgDebugLevel;

Int_t RunTests(Int_t argc, Char_t* argv[]) {
   Int_t result = 0;
   Bool_t validOption(kFALSE);

   if (argc == 1) {
      GoFTStress::fgDebugLevel = GoFTStress::kBasicDebug;
      validOption = kTRUE;
   } else if (argc == 3 && strcmp(argv[1], "--o") == 0) {
      if (strcmp(argv[2], "0") == 0 || strcmp(argv[2], "no-debug") == 0) {
         GoFTStress::fgDebugLevel = GoFTStress::kNoDebug;
         validOption = kTRUE;
      } else if (strcmp(argv[2], "1") == 0 || strcmp(argv[2], "basic-debug") == 0) {
         GoFTStress::fgDebugLevel = GoFTStress::kBasicDebug;
         validOption = kTRUE;
         std::cout << "Debug level: basic." << std::endl;
      } else if (strcmp(argv[2], "2") == 0 || strcmp(argv[2], "standard-debug") == 0) {
         GoFTStress::fgDebugLevel = GoFTStress::kStandardDebug;
         validOption = kTRUE;
         std::cout << "Debug level: standard." << std::endl;
      } else if (strcmp(argv[2], "help") == 0) {
         std::cout << "Showing help: " << std::endl;
      }
   }

   if (!validOption) {
      std::cout << "Please type just one of the following debug levels in either formats preceded by \"--o\":" << std::endl;
      std::cout << "a) no-debug(0) [no test results printing]" << std::endl;
      std::cout << "b) basic-debug(1) [prints either \"OK\" or \"FAIL\" test results]" << std::endl;
      std::cout << "c) standard-debug(2) [prints summarized test results]" << std::endl;
      return 0;
   }

   if (GoFTStress::fgDebugLevel >= GoFTStress::kBasicDebug) {
      std::cout << "*****************************************************************************************" << std::endl;
      std::cout << "*                           Goodness of Fit Test STRESS suite                           *" <<std::endl;
      std::cout << "*****************************************************************************************" << std::endl;
}

   gBenchmark = new TBenchmark();

   TStopwatch timer;
   timer.Start();

   gBenchmark->Start("GoFTestStress");

   GoFTStress* goftStressTest = new GoFTStress;

   result = goftStressTest->RunTests();

   gBenchmark->Stop("GoFTestStress");

   goftStressTest->PrintBenchmark();

   delete goftStressTest;

   delete gBenchmark;

   return result;
}

Int_t stressGoFTest(Int_t argc = 1 , Char_t* argv[] = 0) {
   return RunTests(argc, argv);
}

#if !defined(__CINT__) && !defined(__MAKECINT__)
Int_t main(Int_t argc, Char_t* argv[]) {
   return RunTests(argc, argv);
}
#endif
