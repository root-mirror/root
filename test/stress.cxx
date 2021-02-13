// @(#)root/test:$Id$
// Author: Rene Brun   05/11/98

/////////////////////////////////////////////////////////////////
//
//    R O O T   T E S T  S U I T E  and  B E N C H M A R K S
//    ======================================================
//
// The suite of programs below test the essential parts of Root.
// In particular, there is an extensive test of the I/O and Trees.
// The test can be run in batch or with the interpreter.
// You must run
//   gmake  on Unix systems
//   nmake  on Windows
// To run in batch, do
//   stress -b 1000  : with no parameters, run standard test with 1000 events
//   stress -b 30:     run test with 30 events only
//
// To run interactively, do
// root -b
//  Root > .L stress.cxx
//  Root > stress(1000)  run standard test with 1000 events
//  Root > stress(30)    run with 30 events only
//
// The standard test with 1000 events will create several files.
// The size of all files is around 100 Mbytes.
// The test with 30 events only require around  20 Mbytes
// NB: The test must be run with more than 10 events
//
// The tests runs sequentially 16 tests. Each test will produce
// one line (Test OK or Test failed) with some result parameters.
// At the end of the test a table is printed showing the global results
// with the amount of I/O, Real Time and Cpu Time.
// One single number (ROOTMARKS) is also calculated showing the relative
// performance of your machine compared to a reference machine
// a Pentium IV 2.4 Ghz) with 512 MBytes of memory
// and 120 GBytes IDE disk.
//
// An example of output when all the tests run OK is shown below:
// ******************************************************************
// *  Starting  R O O T - S T R E S S test suite with 1000 events
// ******************************************************************
// Test  1 : Functions, Random Numbers, Histogram Fits............. OK
// Test  2 : Check size & compression factor of a Root file........ OK
// Test  3 : Purge, Reuse of gaps in TFile......................... OK
// Test  4 : Test of 2-d histograms, functions, 2-d fits........... OK
// Test  5 : Test graphics & Postscript............................ OK
// Test  6 : Test subdirectories in a Root file.................... OK
// Test  7 : TNtuple, selections, TCut, TCutG, TEventList.......... OK
// Test  8 : Trees split and compression modes..................... OK
// Test  9 : Analyze Event.root file of stress 8................... OK
// Test 10 : Create 10 files starting from Event.root.............. OK
// Test 11 : Test chains of Trees using the 10 files............... OK
// Test 12 : Compare histograms of test 9 and 11................... OK
// Test 13 : Test merging files of a chain......................... OK
// Test 14 : Check correct rebuilt of Event.root in test 13........ OK
// Test 15 : Divert Tree branches to separate files................ OK
// Test 16 : CINT test (3 nested loops) with LHCb trigger.......... OK
// Test 17 : Test mkdir............................................ OK
// ******************************************************************
//*  Linux pcbrun.cern.ch 2.4.20 #1 Thu Jan 9 12:21:02 MET 2003
//******************************************************************
//stress    : Total I/O =  703.7 Mbytes, I =  535.2, O = 168.5
//stress    : Compr I/O =  557.0 Mbytes, I =  425.1, O = 131.9
//stress    : Real Time =  64.84 seconds Cpu Time =  61.00 seconds
//******************************************************************
//*  ROOTMARKS = 600.1   *  Root4.02/00   20041217/1146
//******************************************************************
//
//_____________________________batch only_____________________
#ifndef __CINT__

#include <cstdlib>
#include <TROOT.h>
#include <TSystem.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TMath.h>
#include <TF1.h>
#include <TF2.h>
#include <TProfile.h>
#include <TKey.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <TPostScript.h>
#include <TNtuple.h>
#include <TBranch.h>
#include <TTreeCache.h>
#include <TChain.h>
#include <TCut.h>
#include <TCutG.h>
#include <TEventList.h>
#include <TBenchmark.h>
#include <TSystem.h>
#include <TApplication.h>
#include <TClassTable.h>
#include <Compression.h>
#include <snprintf.h>
#include "Event.h"

void stress(Int_t nevent, Int_t style, Int_t printSubBenchmark, UInt_t portion );
void stress1();
void stress2();
void stress3();
void stress4();
void stress5();
void stress6();
void stress7();
void stress8(Int_t nevent);
void stress9tree(TTree *tree, Int_t realTestNum);
void stress9();
void stress10();
void stress11();
void stress12(Int_t testid);
void stress13();
void stress14();
void stress15();
void stress16();
void stress17();
void cleanup();


int main(int argc, char **argv)
{
   std::string inclRootSys = ("-I" + TROOT::GetRootSys() + "/test").Data();
   TROOT::AddExtraInterpreterArgs({inclRootSys});

   gROOT->SetBatch();
   TApplication theApp("App", &argc, argv);
   gBenchmark = new TBenchmark();
   Int_t nevent = 1000;      // by default create 1000 events
   if (argc > 1)  nevent = atoi(argv[1]);
   Int_t style  = 1;        // by default the new branch style
   if (argc > 2) style  = atoi(argv[2]);
   Int_t printSubBench = kFALSE;
   if (argc > 3) printSubBench = atoi(argv[3]);
   Int_t portion = 131071;
   if (argc > 4) portion  = atoi(argv[4]);
   stress(nevent, style, printSubBench, portion);
   return 0;
}

#endif

class TH1;
class TTree;

int gPrintSubBench = 0;

//_______________________common part_________________________

Double_t ntotin=0, ntotout=0;

void stress(Int_t nevent, Int_t style = 1,
            Int_t printSubBenchmark = kFALSE, UInt_t portion = 131071)
{
   //Main control function invoking all test programs

   gPrintSubBench = printSubBenchmark;

   if (nevent < 11) nevent = 11; // must have at least 10 events
   //Delete all possible objects in memory (to execute stress several times)
   gROOT->GetListOfFunctions()->Delete();
   gROOT->GetList()->Delete();

   printf("******************************************************************\n");
   printf("*  Starting  R O O T - S T R E S S test suite with %d events\n",nevent);
   printf("******************************************************************\n");
   // select the branch style
   TTree::SetBranchStyle(style);

   //Run the standard test suite
   gBenchmark->Start("stress");
   if (portion&1) stress1();
   if (portion&2) stress2();
   if (portion&4) stress3();
   if (portion&8) stress4();
   if (portion&16) stress5();
   if (portion&32) stress6();
   if (portion&64) stress7();
   if (portion&128) stress8(nevent);
   if (portion&256) stress9();
   if (portion&512) stress10();
   if (portion&1024) stress11();
   if (portion&2048) stress12(12);
   if (portion&4096) stress13();
   if (portion&8192) stress14();
   if (portion&16384) stress15();
   if (portion&32768) stress16();
   if (portion&65536) stress17();
   gBenchmark->Stop("stress");

   cleanup();

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
      const char *os = gSystem->Getenv("OS");
      if (!os) printf("*  SYS: Windows 95\n");
      else     printf("*  SYS: %s %s \n",os,gSystem->Getenv("PROCESSOR_IDENTIFIER"));
   }

   printf("******************************************************************\n");
   Float_t mbtot = (Float_t)(ntotin+ntotout)/1000000.;
   Float_t mbin  = (Float_t)ntotin/1000000.;
   Float_t mbout = (Float_t)ntotout/1000000.;
   printf("stress    : Total I/O =%7.1f Mbytes, I =%7.1f, O =%6.1f\n",mbtot,mbin,mbout);
   Float_t mbin1  = (Float_t)(TFile::GetFileBytesRead()/1000000.);
   Float_t mbout1 = (Float_t)(TFile::GetFileBytesWritten()/1000000.);
   Float_t mbtot1 = mbin1+mbout1;
   printf("stress    : Compr I/O =%7.1f Mbytes, I =%7.1f, O =%6.1f\n",mbtot1,mbin1,mbout1);
   gBenchmark->Print("stress");
#ifndef __CINT__
   Float_t cp_brun_30   = 12.73;
   Float_t cp_brun_1000 = 61.88;
#else
   Float_t cp_brun_30   = 31.03;  //The difference is essentially coming from stress16
   Float_t cp_brun_1000 = 84.30;
#endif
   Float_t cp_brun = cp_brun_1000 - (cp_brun_1000 - cp_brun_30)*(1000-nevent)/(1000-30);
   Float_t ct = gBenchmark->GetCpuTime("stress");
   Float_t rootmarks = 600*cp_brun/ct;
   printf("******************************************************************\n");
   printf("*  ROOTMARKS =%6.1f   *  Root%-8s  %d/%d\n",rootmarks,gROOT->GetVersion(),gROOT->GetVersionDate(),gROOT->GetVersionTime());
   printf("******************************************************************\n");

   delete gBenchmark;
}

////////////////////////////////////////////////////////////////////////////////
///Compute a function sum of 3 gaussians

Double_t f1int(Double_t *x, Double_t *p)
{
   Double_t e1 = (x[0]-p[1])/p[2];
   Double_t e2 = (x[0]-p[4])/p[5];
   Double_t e3 = (x[0]-p[7])/p[8];
   Double_t f  = p[0]*TMath::Exp(-0.5*e1*e1)
                +p[3]*TMath::Exp(-0.5*e2*e2)
                +p[6]*TMath::Exp(-0.5*e3*e3);
   return f;
}

////////////////////////////////////////////////////////////////////////////////
/// Print test program number and its title

void Bprint(Int_t id, const char *title)
{
   const Int_t kMAX = 65;
   char header[80];
   snprintf(header,80,"Test %2d : %s",id,title);
   Int_t nch = strlen(header);
   for (Int_t i=nch;i<kMAX;i++) header[i] = '.';
   header[kMAX] = 0;
   header[kMAX-1] = ' ';
   printf("%s",header);
}

////////////////////////////////////////////////////////////////////////////////
///Generate two functions supposed to produce the same result
///One function "f1form" will be computed by the TFormula class
///The second function "f1int" will be
///   - compiled when running in batch mode
///   - interpreted by CINT when running in interactive mode

void stress1()
{
   Bprint(1,"Functions, Random Numbers, Histogram Fits");

   //Start with a function inline expression (managed by TFormula)
   Double_t f1params[9] = {100,-3,3,60,0,0.5,40,4,0.7};
   TF1 *f1form = new TF1("f1form","gaus(0)+gaus(3)+gaus(6)",-10,10);
   f1form->SetParameters(f1params);

   //Create an histogram and fill it randomly with f1form
   gRandom->SetSeed(65539);
   TH1F *h1form = new TH1F("h1form","distribution from f1form",100,-10,10);
   TH1F *h1diff = (TH1F*)h1form->Clone();
   h1diff->SetName("h1diff");
   h1form->FillRandom("f1form",10000);

   //Fit h1form with original function f1form
   h1form->Fit("f1form","q0");

   // std::cout << "done formula" << std::endl;
   // f1form->Print("v");

   //same operation with an interpreted function f1int
   TF1 *f1 = new TF1("f1int",f1int,-10,10,9);
   f1->SetParameters(f1params);

   //Create an histogram and fill it randomly with f1int
   gRandom->SetSeed(65539); //make sure we start with the same random numbers
   TH1F *h1int = new TH1F("h1int","distribution from f1int",100,-10,10);
   h1int->FillRandom("f1int",10000);

   //Fit h1int with original function f1int
   h1int->Fit("f1int","q0");

   //The difference between the two histograms must be null
   h1diff->Add(h1form, h1int, 1, -1);
   Double_t hdiff = h1diff->Integral(0,101);

   //Compare fitted parameters and value of integral of f1form in [-8,6]
   Int_t npar = f1form->GetNpar();
   Double_t pdiff, pdifftot = 0;
   for (Int_t i=0;i<npar;i++) {
      pdiff = (f1form->GetParameter(i) - f1->GetParameter(i))/f1form->GetParameter(i);
      pdifftot += TMath::Abs(pdiff);
   }
   // The integral in the range [-8,6] must be = 1923.74578
   Double_t rint = TMath::Abs(f1form->Integral(-8,6) - 1923.74578);

   //Some slight differences are authorized to take into account
   //different math libraries used by the compiler, CINT and TFormula
   Bool_t OK = kTRUE;
   if (hdiff > 0.1 || pdifftot > 2.e-3 || rint > 10) OK = kFALSE;
   if (OK) printf("OK\n");
   else    {
      printf("FAILED\n");
      printf("%-8s hdiff=%g, pdifftot=%g, rint=%g\n"," ",hdiff,pdifftot,rint);
   }
   if (gPrintSubBench) { printf("Test  1 : "); gBenchmark->Show("stress");gBenchmark->Start("stress"); }
   //Save all objects in a Root file (will be checked by stress2)
   TFile local("stress.root","recreate");
   f1form->Write();
   f1->Write();
   h1form->Write();
   h1int->Write();
   ntotout += local.GetBytesWritten();
   //do not close the file. should be done by the destructor automatically
   // delete h1int;
   // delete h1form;
   // delete h1diff;
}

////////////////////////////////////////////////////////////////////////////////
///check length and compression factor in stress.root

void stress2()
{
   Bprint(2,"Check size & compression factor of a Root file");
   TFile f("stress.root");
   Long64_t last = f.GetEND();
   Float_t comp = f.GetCompressionFactor();

   Bool_t OK = kTRUE;
   //Long64_t lastgood = 12383; //9428;
   //Long64_t lastgood = 9789;  // changes for new TFormula
   //Long64_t lastgood = 9797;  // changes for TH1 v8 ROOT-9173 on 32-bits
#ifdef R__HAS_DEFAULT_LZ4
      Long64_t lastgood = 10733;
      if (last < lastgood - 200 || last > lastgood + 200 || comp < 1.5 || comp > 2.1)
         OK = kFALSE;
#else
#ifdef R__HAS_CLOUDFLARE_ZLIB
      Long64_t lastgood = 9813;
#else
      Long64_t lastgood = 10034;  // changes in TFormula (v12)
#endif
      if (last < lastgood - 200 || last > lastgood + 200 || comp < 2.0 || comp > 2.4)
         OK = kFALSE;
#endif
   if (OK) printf("OK\n");
   else    {
      printf("FAILED\n");
      printf("%-8s last =%lld, comp=%f\n"," ",last,comp);
   }
   if (gPrintSubBench) { printf("Test  2 : "); gBenchmark->Show("stress");gBenchmark->Start("stress"); }
}

////////////////////////////////////////////////////////////////////////////////
///Open stress.root, read all objects, save 10 times and purge
///This function tests the generation and reuse of gaps in files

void stress3()
{
   Bprint(3,"Purge, Reuse of gaps in TFile");
   TFile f("stress.root","update");
   f.ReadAll();
   for (Int_t i=0;i<10;i++) {
      f.Write();
   }
   f.Purge();
   f.Write();

   //check length and compression level in stress.root
   ntotin  += f.GetBytesRead();
   ntotout += f.GetBytesWritten();
   Long64_t last = f.GetEND();
   Float_t comp = f.GetCompressionFactor();
   Bool_t OK = kTRUE;
#ifdef R__HAS_CLOUDFLARE_ZLIB
   constexpr Long64_t lastgood = 52027;
#else
   constexpr Long64_t lastgood = 51886;
#endif
   constexpr Long64_t tolerance = 150;
#ifdef R__HAS_DEFAULT_LZ4
      constexpr Long64_t difflastgoodlz4 = 5500;
      if (last < lastgood - tolerance || last > lastgood + difflastgoodlz4 + tolerance || comp < 1.5 || comp > 2.1)
         OK = kFALSE;
#else
      if (last < lastgood - tolerance || last > lastgood + tolerance || comp < 1.8 || comp > 2.4)
         OK = kFALSE;
#endif
   if (OK) printf("OK\n");
   else    {
      printf("FAILED\n");
#ifdef R__HAS_DEFAULT_LZ4
      printf("%-8s LZ4 file size= %lld (expected %lld +/- %lld)\n"
             "%-8s Comp Fact=  %3.2f (expected 1.8 +/- 0.3)\n",
             " ", last, lastgood + difflastgoodlz4, tolerance, " ", comp);
#else
       printf("%-8s File size= %lld (expected %lld +/- %lld)\n"
             "%-8s Comp Fact=  %3.2f (expected 2.1 +/- 0.3)\n",
             " ", last, lastgood, tolerance, " ", comp);
#endif
   }
   if (gPrintSubBench) { printf("Test  3 : "); gBenchmark->Show("stress");gBenchmark->Start("stress"); }
}

////////////////////////////////////////////////////////////////////////////////
/// Test of 2-d histograms, functions, 2-d fits

void stress4()
{
   Bprint(4,"Test of 2-d histograms, functions, 2-d fits");

   Double_t f2params[15] = {100,-3,3,-3,3,160,0,0.8,0,0.9,40,4,0.7,4,0.7};
   TF2 *f2form = new TF2("f2form","xygaus(0)+xygaus(5)+xygaus(10)",-10,10,-10,10);
   f2form->SetParameters(f2params);

   //Create an histogram and fill it randomly with f2form
   gRandom->SetSeed(65539);
   TH2F *h2form = new TH2F("h2form","distribution from f2form",40,-10,10,40,-10,10);
   Int_t nentries = 100000;
   h2form->FillRandom("f2form",nentries);
   //Fit h2form with original function f2form
   Float_t ratio = 4*nentries/100000;
   f2params[ 0] *= ratio;
   f2params[ 5] *= ratio;
   f2params[10] *= ratio;
   f2form->SetParameters(f2params);
   h2form->Fit("f2form","q0");
   //Update stress.root
   TFile f("stress.root","update");
   h2form->Write();
   f2form->Write();

   ntotin  += f.GetBytesRead();
   ntotout += f.GetBytesWritten();

   //Compare results of fit with expected parameters
   Bool_t OK = kTRUE;
   for (int k = 0; k < 3; ++k) {
      for (int  l = 1; l < 5; ++l) {
         int idx = k*5+l;
         Double_t dp0  = TMath::Abs((f2form->GetParameter(idx) -f2params[idx]));
         if (f2params[idx] != 0.) dp0 /=  f2params[idx];
         bool testok =  (dp0 < 5.e-2);
         if (!testok) {
            printf("\nFAILED:   ipar=%d delta=%g, par=%g, nom=%g",idx,dp0,f2form->GetParameter(idx),f2params[idx]);
         }
         OK &= testok;
      }
   }
   if (OK) printf("OK\n");
   else    printf("\ntest FAILED !\n");
   if (gPrintSubBench) { printf("Test  4 : "); gBenchmark->Show("stress");gBenchmark->Start("stress"); }
}

////////////////////////////////////////////////////////////////////////////////
/// Test of Postscript.
/// Make a complex picture. Verify number of lines on ps file
/// Testing automatically the graphics package is a complex problem.
/// The best way we have found is to generate a Postscript image
/// of a complex canvas containing many objects.
/// The number of lines in the ps file is compared with a reference run.
/// A few lines (up to 2 or 3) of difference may be expected because
/// Postscript works with floats. The date and time of the run are also
/// different.
/// You can also inspect visually the ps file with a ps viewer.

void stress5()
{
   Bprint(5,"Test graphics & Postscript");

   TCanvas *c1 = new TCanvas("c1","stress canvas",800,600);
   gROOT->LoadClass("TPostScript","Postscript");
   TPostScript ps("stress.ps",112);

   //Get objects generated in previous test
   TFile f("stress.root");
   TF1  *f1form = (TF1*)f.Get("f1form");
   TF2  *f2form = (TF2*)f.Get("f2form");
   TH1F *h1form = (TH1F*)f.Get("h1form");
   TH2F *h2form = (TH2F*)f.Get("h2form");

   //Divide the canvas in subpads. Plot with different options
   c1->Divide(2,2);
   c1->cd(1);
   f1form->Draw();
   c1->cd(2);
   h1form->Draw();
   c1->cd(3);
   h2form->Draw("box");
   f2form->Draw("cont1same");
   c1->cd(4);
   f2form->Draw("surf");

   ps.Close();

   //count number of lines in ps file
   FILE *fp = fopen("stress.ps","r");
   char line[260];
   Int_t nlines = 0;
   Int_t nlinesGood = 632;
   while (fgets(line,255,fp)) {
      nlines++;
   }
   fclose(fp);
   ntotin  += f.GetBytesRead();
   ntotout += f.GetBytesWritten();
   Bool_t OK = kTRUE;
   if (nlines < nlinesGood-110 || nlines > nlinesGood+110) OK = kFALSE;
   if (OK) printf("OK\n");
   else    {
      printf("FAILED\n");
      printf("%-8s nlines in stress.ps file = %d\n"," ",nlines);
   }
   delete c1;
   if (gPrintSubBench) { printf("Test  5 : "); gBenchmark->Show("stress");gBenchmark->Start("stress"); }

}

////////////////////////////////////////////////////////////////////////////////
/// Test subdirectories in a Root file
/// Create many TH1S histograms, make operations between them

void stress6()
{
   Bprint(6,"Test subdirectories in a Root file");

   TFile f("stress.root","update");
   // create a new subdirectory for each plane
   gRandom->SetSeed(65539);
   const Int_t nplanes = 10;
   const Int_t ncounters = 100;
   char dirname[50];
   char hname[20];
   char htitle[80];
   TH1S *hn[ncounters];
   TH1S *hs[ncounters];
   Int_t i,j,k,id;
   TH1F *hsumPlanes = new TH1F("hsumPlanes","Sum of all planes",100,0,100);
   //Create a subdirectory per detector plane
   for (i=0;i<nplanes;i++) {
      snprintf(dirname,50,"plane%d",i);
      TDirectory *cdplane = f.mkdir(dirname);
      if (cdplane == 0) continue;
      cdplane->cd();
      // create counter histograms
      for (j=0;j<ncounters;j++) {
         snprintf(hname,20,"h%d_%dN",i,j);
         snprintf(htitle,80,"hist for counter:%d in plane:%d North",j,i);
         hn[j] = new TH1S(hname,htitle,100,0,100);
         snprintf(hname,20,"h%d_%dS",i,j);
         snprintf(htitle,80,"hist for counter:%d in plane:%d South",j,i);
         hs[j] = new TH1S(hname,htitle,100,0,100);
      }
      // fill counter histograms randomly
      for (k=0;k<10000;k++) {
         id = Int_t(ncounters*gRandom->Rndm());
         hn[id]->Fill(gRandom->Gaus(60,10));
         hs[id]->Fill(gRandom->Gaus(40,5));
      }
      // Write all objects in directory in memory to disk
      cdplane->Write();
      // Delete all objects from memory
      cdplane->GetList()->Delete();
      f.cd();
   }
   // Now read back all objects from all subdirectories
   // Add North and south histograms in hsumPlanes
   for (i=0;i<nplanes;i++) {
      snprintf(dirname,50,"plane%d",i);
      f.cd(dirname);
      for (j=0;j<ncounters;j++) {
         snprintf(hname,20,"h%d_%dN",i,j);
         TH1S *hnorth; gDirectory->GetObject(hname,hnorth);
         snprintf(hname,20,"h%d_%dS",i,j);
         TH1S *hsouth; gDirectory->GetObject(hname,hsouth);
         if (hnorth == 0 || hsouth == 0) continue;
         hsumPlanes->Add(hnorth);
         hsumPlanes->Add(hsouth);
         delete hnorth; delete hsouth;
      }
      f.cd();    // change current directory to top
   }
   // Verify number of entries, rms and mean value
   ntotin  += f.GetBytesRead();
   ntotout += f.GetBytesWritten();
   Int_t nentries = (Int_t)hsumPlanes->GetEntries();
   Double_t rms   = hsumPlanes->GetRMS();
   Double_t mean  = hsumPlanes->GetMean();
   Int_t nentriesGood = 200000;
   Double_t rmsGood  = 12.745;
   Double_t meanGood = 50.01;
   Double_t diffrms  = TMath::Abs(rmsGood -rms)/rmsGood;
   Double_t diffmean = TMath::Abs(meanGood -mean)/meanGood;
   Bool_t OK = kTRUE;
   if (nentriesGood != nentries || diffrms > 1.e-2 || diffmean > 1.e-2) OK = kFALSE;
   if (OK) printf("OK\n");
   else    {
      printf("FAILED\n");
      printf("%-8s nentries=%d, diffmean=%g, diffrms=%g\n"," ",nentries,diffmean,diffrms);
   }
   if (gPrintSubBench) { printf("Test  6 : "); gBenchmark->Show("stress");gBenchmark->Start("stress"); }
}

////////////////////////////////////////////////////////////////////////////////
/// Test TNtuple class with several selection mechanisms
/// Test expression cuts
/// Test graphical cuts
/// Test event lists and operations on event lists
/// Compare results of TTree::Draw with results of an explict loop

void stress7()
{
   Bprint(7,"TNtuple, selections, TCut, TCutG, TEventList");

   TFile f("stress.root","update");
   // Create and fill a TNtuple
   gRandom->SetSeed(65539);
   TNtuple *ntuple = new TNtuple("ntuple","Demo ntuple","px:py:pz:random:i");
   Float_t px, py, pz;
   Int_t nall = 50000;
   Int_t i;
   for (i = 0; i < nall; i++) {
      gRandom->Rannor(px,py);
      pz = px*px + py*py;
      Float_t random = gRandom->Rndm();
      ntuple->Fill(px,py,pz,random,i);
   }
   ntuple->Write();

   // Create a graphical cut. Select only events in cut
   TCutG *cutg = new TCutG("cutg",9);
   cutg->SetVarX("py");
   cutg->SetVarY("px");
   cutg->SetPoint(0,-1.75713,2.46193);
   cutg->SetPoint(1,-2.58656,-0.786802);
   cutg->SetPoint(2,-0.179195,-0.101523);
   cutg->SetPoint(3,2.12702,-1.49746);
   cutg->SetPoint(4,2.2484,1.95431);
   cutg->SetPoint(5,0.630004,0.583756);
   cutg->SetPoint(6,-0.381495,2.28426);
   cutg->SetPoint(7,-1.27161,1.01523);
   cutg->SetPoint(8,-1.75713,2.46193);
   TH2F *hpxpy = new TH2F("hpxpy","px vx py with cutg",40,-4,4,40,-4,4);
   ntuple->Draw("px:py>>hpxpy","cutg","goff");
   Int_t npxpy = (Int_t)hpxpy->GetEntries();
   Int_t npxpyGood = 27918;
   hpxpy->Write();
   cutg->Write();
   delete cutg;

   // Fill a TEventList using the standard cut
   ntuple->Draw(">>elist","py<0 && pz>4 && random<0.5","goff");
   TEventList *elist; gDirectory->GetObject("elist",elist);
   // Fill hist htemp using the standard cut
   ntuple->Draw("px>>htemp0","py<0 && pz>4 && random<0.5","goff");
   TH1 *htemp0;  gDirectory->GetObject("htemp0",htemp0);
   Double_t pxmean0 = htemp0->GetMean();
   Double_t pxrms0  = htemp0->GetRMS();

   // Fill hist hcut using a TCut = the standard cut
   TCut cut1 = "py<0 && pz>4 && random<0.5";
   TCut vcut = "px>>hcut";
   ntuple->Draw(vcut,cut1,"goff");
   // Fill hist helist looping on the eventlist in TTree::Draw
   ntuple->SetEventList(elist);
   ntuple->Draw("px>>helist","","goff");
   ntuple->SetEventList(0);
   TH1 *hcut;   gDirectory->GetObject("hcut",hcut);
   TH1 *helist; gDirectory->GetObject("helist",helist);
   Int_t n1 = (Int_t)hcut->GetEntries();
   Int_t n2 = (Int_t)helist->GetEntries();
   htemp0->Write();
   cut1.Write();
   helist->Write();
   hcut->Write();

   // now loop on eventlist explicitly and fill helist again
   Float_t pxr;
   ntuple->SetBranchAddress("px",&pxr);
   TH1 *helistc = (TH1*)helist->Clone();
   helistc->Reset();
   helistc->SetName("helistc");
   Int_t nlist = elist->GetN();
   for (i=0;i<nlist;i++) {
      Long64_t event = elist->GetEntry(i);
      ntuple->GetEntry(event);
      helistc->Fill(pxr);
   }
   Int_t n3 = (Int_t)helistc->GetEntries();
   Double_t pxmean2 = helistc->GetMean();
   Double_t pxrms2  = helistc->GetRMS();
   helistc->Write();
   elist->Write();

   // Generate several TEventlist objects + total and save them
   char elistname[20];
   char cutname[20];
   TEventList *el[10];
   TEventList *elistall = new TEventList("elistall","Sum of all cuts");
   for (i=0;i<10;i++) {
      snprintf(elistname,20,">>elist%d",i);
      snprintf(cutname,20,"i 10 == %d",i); cutname[1] ='%';
      ntuple->Draw(elistname,cutname,"goff");
      gDirectory->GetObject(&elistname[2],el[i]);
      el[i]->Write();
      elistall->Add(el[i]);
   }
   elistall->Write();

   // Read big list from file and check that the distribution with the list
   // correspond to all events (no cuts)
   delete ntuple;
   TNtuple *nt; gDirectory->GetObject("ntuple",nt);
   nt->SetBranchAddress("px",&pxr);
   TH1F *hpx = new TH1F("hpx","hpx",100,-3,3);
   nt->Draw("px>>hpx","","goff");
   TEventList *all; gDirectory->GetObject("elistall",all);
   nt->SetEstimate(nall); //must be done because the order in eventlist is different
   nt->SetEventList(all);
   TH1F *hall = (TH1F*)hpx->Clone();
   hall->SetName("hall");
   nt->Draw("px>>hall","","goff");
   // Take the difference between the two histograms. Must be empty
   //TH1F hcomp = (*hall) - (*hpx);
   //Double_t compsum = hcomp.GetSum();
   hall->Add(hpx,-1);
   Double_t compsum = hall->GetSum();
   ntotin  += f.GetBytesRead();
   ntotout += f.GetBytesWritten();

   // We can compare entries, means and rms
   Bool_t OK = kTRUE;
   if (n1 != n2 || n1 != n3 || n3 != nlist || nall !=elistall->GetN()
                || npxpy != npxpyGood
                || compsum != 0
                || TMath::Abs(pxmean0-pxmean2) > 0.1
                || TMath::Abs(pxrms0-pxrms2) > 0.01) OK = kFALSE;
   if (OK) printf("OK\n");
   else    {
      printf("FAILED\n");
      printf("%-8s n1=%d, n2=%d, n3=%d, elistallN=%d\n"," ",n1,n2,n3,elistall->GetN());
      printf("%-8s pxmean0=%g, pxmean2=%g, pxrms0=%g\n"," ",pxmean0,pxmean2,pxrms0);
      printf("%-8s pxrms2=%g, compsum=%g, npxpy=%d\n"," ",pxrms2,compsum,npxpy);
   }
   if (gPrintSubBench) { printf("Test  7 : "); gBenchmark->Show("stress");gBenchmark->Start("stress"); }
}

////////////////////////////////////////////////////////////////////////////////
///  Read the event file
///  Loop on all events in the file (reading everything).
///  Count number of bytes read

Int_t stress8read(Int_t nevent)
{
   TFile *hfile = new TFile("Event.root");
   TTree *tree; hfile->GetObject("T",tree);
   Event *event = 0;
   tree->SetBranchAddress("event",&event);
   Int_t nentries = (Int_t)tree->GetEntries();
   Int_t nev = TMath::Max(nevent,nentries);
   //activate the treeCache
   Int_t cachesize = 10000000; //this is the default value: 10 MBytes
   tree->SetCacheSize(cachesize);
   TTreeCache::SetLearnEntries(1); //one entry is sufficient to learn
   TTreeCache *tc = (TTreeCache*)hfile->GetCacheRead();
   tc->SetEntryRange(0,nevent);
   Int_t nb = 0;
   for (Int_t ev = 0; ev < nev; ev++) {
      nb += tree->GetEntry(ev);        //read complete event in memory
   }
   ntotin  += hfile->GetBytesRead();

   delete event;
   delete hfile;
   return nb;
}


////////////////////////////////////////////////////////////////////////////////
///  Create the Event file in various modes
/// comp = compression level
/// split = 1 split mode, 0 = no split

Int_t stress8write(Int_t nevent, Int_t comp, Int_t split)
{
   // Create the Event file, the Tree and the branches
   TFile *hfile = new TFile("Event.root","RECREATE","TTree benchmark ROOT file");
   hfile->SetCompressionLevel(comp);

   // Create one event
   Event *event = new Event();

   // Create a ROOT Tree and one superbranch
   TTree *tree = new TTree("T","An example of a ROOT tree");
   tree->SetAutoSave(100000000);  // autosave when 100 Mbytes written
   Int_t bufsize = 64000;
   if (split)  bufsize /= 4;
   tree->Branch("event", &event, bufsize,split);

   //Fill the Tree
   Int_t ev, nb=0, meanTracks=600;
   Float_t ptmin = 1;
   for (ev = 0; ev < nevent; ev++) {
      event->Build(ev,meanTracks,ptmin);

      nb += tree->Fill();  //fill the tree
   }
   hfile->Write();
   ntotout += hfile->GetBytesWritten();
   delete event;
   delete hfile;
   return nb;
}


////////////////////////////////////////////////////////////////////////////////
///  Run the $ROOTSYS/test/Event program in several configurations.

void stress8(Int_t nevent)
{
   Bprint(8,"Trees split and compression modes");

  // First step: make sure the Event shared library exists
  // This test dynamic linking when running in interpreted mode
   if (!TClassTable::GetDict("Event")) {
      Int_t st1 = -1;
      if (gSystem->DynamicPathName("$ROOTSYS/test/libEvent",kTRUE)) {
         st1 = gSystem->Load("$(ROOTSYS)/test/libEvent");
      }
      if (st1 == -1) {
         if (gSystem->DynamicPathName("test/libEvent",kTRUE)) {
            st1 = gSystem->Load("test/libEvent");
         }
         if (st1 == -1) {
            printf("===>stress8 will try to build the libEvent library\n");
            Bool_t UNIX = strcmp(gSystem->GetName(), "Unix") == 0;
            if (UNIX) gSystem->Exec("(cd $ROOTSYS/test; make Event)");
            else      gSystem->Exec("(cd %ROOTSYS%\\test && nmake libEvent.dll)");
            st1 = gSystem->Load("$(ROOTSYS)/test/libEvent");
         }
      }
   }

   // Create the file not compressed, in no-split mode and read it back
   gRandom->SetSeed(65539);
   Int_t nbw0 = stress8write(100,0,0);
   Int_t nbr0 = stress8read(0);
   Event::Reset();

   // Create the file compressed, in no-split mode and read it back
   gRandom->SetSeed(65539);
   Int_t nbw1 = stress8write(100,1,0);
   Int_t nbr1 = stress8read(0);
   Event::Reset();

   // Create the file compressed, in split mode and read it back
   gRandom->SetSeed(65539);
   Int_t nbw2 = stress8write(nevent,1,9);
   Int_t nbr2 = stress8read(0);
   Event::Reset();

   Bool_t OK = kTRUE;
   if (nbw0 != nbr0 || nbw1 != nbr1 || nbw2 != nbr2) OK = kFALSE;
   if (nbw0 != nbw1) OK = kFALSE;
   if (OK) printf("OK\n");
   else    {
      printf("FAILED\n");
      printf("%-8s nbw0=%d, nbr0=%d, nbw1=%d\n"," ",nbw0,nbr0,nbw1);
      printf("%-8s nbr1=%d, nbw2=%d, nbr2=%d\n"," ",nbr1,nbw2,nbr2);
   }
   if (gPrintSubBench) { printf("Test  8 : "); gBenchmark->Show("stress");gBenchmark->Start("stress"); }
}

////////////////////////////////////////////////////////////////////////////////
/// Compare histograms h1 and h2
/// Check number of entries, mean and rms
/// if means differ by more than 1/1000 of the range return -1
/// if rms differs in percent by more than 1/1000 return -2
/// Otherwise return difference of number of entries

Int_t HistCompare(TH1 *h1, TH1 *h2)
{
   Int_t n1       = (Int_t)h1->GetEntries();
   Double_t mean1 = h1->GetMean();
   Double_t rms1  = h1->GetRMS();
   Int_t n2       = (Int_t)h2->GetEntries();
   Double_t mean2 = h2->GetMean();
   Double_t rms2  = h2->GetRMS();
   Float_t xrange = h1->GetXaxis()->GetXmax() - h1->GetXaxis()->GetXmin();
   if (TMath::Abs((mean1-mean2)/xrange) > 0.001*xrange) return -1;
   if (rms1 && TMath::Abs((rms1-rms2)/rms1) > 0.001)    return -2;
   return n1-n2;
}

////////////////////////////////////////////////////////////////////////////////
/// Test selections via TreeFormula
/// tree is a TTree when called by stress9
/// tree is a TChain when called from stress11
/// This is a quite complex test checking the results of TTree::Draw
/// or TChain::Draw with an explicit loop on events.
/// Also a good test for the interpreter

void stress9tree(TTree *tree, Int_t realTestNum)
{
   Event *event = 0;
   tree->SetBranchAddress("event",&event);
   gROOT->cd();
   TDirectory *hfile = gDirectory;
   Double_t nrsave = TFile::GetFileBytesRead();

   // Each tree->Draw generates an histogram
   tree->Draw("fNtrack>>hNtrack",    "","goff");
   tree->Draw("fNseg>>hNseg",        "","goff");
   tree->Draw("fTemperature>>hTemp", "","goff");
   tree->Draw("fH.GetMean()>>hHmean","","goff");
   tree->Draw("fTracks.fPx>>hPx","fEvtHdr.fEvtNum%10 == 0","goff");
   tree->Draw("fTracks.fPy>>hPy","fEvtHdr.fEvtNum%10 == 0","goff");
   tree->Draw("fTracks.fPz>>hPz","fEvtHdr.fEvtNum%10 == 0","goff");
   tree->Draw("fRandom>>hRandom","fEvtHdr.fEvtNum%10 == 1","goff");
   tree->Draw("fMass2>>hMass2",  "fEvtHdr.fEvtNum%10 == 1","goff");
   tree->Draw("fBx>>hBx",        "fEvtHdr.fEvtNum%10 == 1","goff");
   tree->Draw("fBy>>hBy",        "fEvtHdr.fEvtNum%10 == 1","goff");
   tree->Draw("fXfirst>>hXfirst","fEvtHdr.fEvtNum%10 == 2","goff");
   tree->Draw("fYfirst>>hYfirst","fEvtHdr.fEvtNum%10 == 2","goff");
   tree->Draw("fZfirst>>hZfirst","fEvtHdr.fEvtNum%10 == 2","goff");
   tree->Draw("fXlast>>hXlast",  "fEvtHdr.fEvtNum%10 == 3","goff");
   tree->Draw("fYlast>>hYlast",  "fEvtHdr.fEvtNum%10 == 3","goff");
   tree->Draw("fZlast>>hZlast",  "fEvtHdr.fEvtNum%10 == 3","goff");
   tree->Draw("fCharge>>hCharge","fPx < 0","goff");
   tree->Draw("fNpoint>>hNpoint","fPx < 0","goff");
   tree->Draw("fValid>>hValid",  "fPx < 0","goff");

   tree->Draw("fMatrix>>hFullMatrix","","goff");
   tree->Draw("fMatrix[][0]>>hColMatrix","","goff");
   tree->Draw("fMatrix[1][]>>hRowMatrix","","goff");
   tree->Draw("fMatrix[2][2]>>hCellMatrix","","goff");

   tree->Draw("fMatrix - fVertex>>hFullOper","","goff");
   tree->Draw("fMatrix[2][1] - fVertex[5][1]>>hCellOper","","goff");
   tree->Draw("fMatrix[][1]  - fVertex[5][1]>>hColOper","","goff");
   tree->Draw("fMatrix[2][]  - fVertex[5][2]>>hRowOper","","goff");
   tree->Draw("fMatrix[2][]  - fVertex[5][]>>hMatchRowOper","","goff");
   tree->Draw("fMatrix[][2]  - fVertex[][1]>>hMatchColOper","","goff");
   tree->Draw("fMatrix[][2]  - fVertex[][]>>hRowMatOper","","goff");
   tree->Draw("fMatrix[][2]  - fVertex[5][]>>hMatchDiffOper","","goff");
   tree->Draw("fMatrix[][]   - fVertex[][]>>hFullOper2","","goff");

   if (gPrintSubBench) { printf("\n"); printf("Test %2dD: ",realTestNum); gBenchmark->Show("stress");gBenchmark->Start("stress"); }

   ntotin  += TFile::GetFileBytesRead() -nrsave;

   //Get pointers to the histograms generated above
   TH1 *hNtrack = (TH1*)hfile->Get("hNtrack");
   TH1 *hNseg   = (TH1*)hfile->Get("hNseg");
   TH1 *hTemp   = (TH1*)hfile->Get("hTemp");
   TH1 *hHmean  = (TH1*)hfile->Get("hHmean");
   TH1 *hPx     = (TH1*)hfile->Get("hPx");
   TH1 *hPy     = (TH1*)hfile->Get("hPy");
   TH1 *hPz     = (TH1*)hfile->Get("hPz");
   TH1 *hRandom = (TH1*)hfile->Get("hRandom");
   TH1 *hMass2  = (TH1*)hfile->Get("hMass2");
   TH1 *hBx     = (TH1*)hfile->Get("hBx");
   TH1 *hBy     = (TH1*)hfile->Get("hBy");
   TH1 *hXfirst = (TH1*)hfile->Get("hXfirst");
   TH1 *hYfirst = (TH1*)hfile->Get("hYfirst");
   TH1 *hZfirst = (TH1*)hfile->Get("hZfirst");
   TH1 *hXlast  = (TH1*)hfile->Get("hXlast");
   TH1 *hYlast  = (TH1*)hfile->Get("hYlast");
   TH1 *hZlast  = (TH1*)hfile->Get("hZlast");
   TH1 *hCharge = (TH1*)hfile->Get("hCharge");
   TH1 *hNpoint = (TH1*)hfile->Get("hNpoint");
   TH1 *hValid  = (TH1*)hfile->Get("hValid");

   TH1F *hFullMatrix    = (TH1F*)hfile->Get("hFullMatrix");
   TH1F *hColMatrix     = (TH1F*)hfile->Get("hColMatrix");
   TH1F *hRowMatrix     = (TH1F*)hfile->Get("hRowMatrix");
   TH1F *hCellMatrix    = (TH1F*)hfile->Get("hCellMatrix");
   TH1F *hFullOper      = (TH1F*)hfile->Get("hFullOper");
   TH1F *hCellOper      = (TH1F*)hfile->Get("hCellOper");
   TH1F *hColOper       = (TH1F*)hfile->Get("hColOper");
   TH1F *hRowOper       = (TH1F*)hfile->Get("hRowOper");
   TH1F *hMatchRowOper  = (TH1F*)hfile->Get("hMatchRowOper");
   TH1F *hMatchColOper  = (TH1F*)hfile->Get("hMatchColOper");
   TH1F *hRowMatOper    = (TH1F*)hfile->Get("hRowMatOper");
   TH1F *hMatchDiffOper = (TH1F*)hfile->Get("hMatchDiffOper");
   TH1F *hFullOper2     = (TH1F*)hfile->Get("hFullOper2");

   //We make clones of the generated histograms
   //We set new names and reset the clones.
   //We want to have identical histogram limits
   TH1F *bNtrack = (TH1F*)hNtrack->Clone(); bNtrack->SetName("bNtrack"); bNtrack->Reset();
   TH1F *bNseg   = (TH1F*)hNseg->Clone();   bNseg->SetName("bNseg");     bNseg->Reset();
   TH1F *bTemp   = (TH1F*)hTemp->Clone();   bTemp->SetName("bTemp");     bTemp->Reset();
   TH1F *bHmean  = (TH1F*)hHmean->Clone();  bHmean->SetName("bHmean");   bHmean->Reset();
   TH1F *bPx     = (TH1F*)hPx->Clone();     bPx->SetName("bPx");         bPx->Reset();
   TH1F *bPy     = (TH1F*)hPy->Clone();     bPy->SetName("bPy");         bPy->Reset();
   TH1F *bPz     = (TH1F*)hPz->Clone();     bPz->SetName("bPz");         bPz->Reset();
   TH1F *bRandom = (TH1F*)hRandom->Clone(); bRandom->SetName("bRandom"); bRandom->Reset();
   TH1F *bMass2  = (TH1F*)hMass2->Clone();  bMass2->SetName("bMass2");   bMass2->Reset();
   TH1F *bBx     = (TH1F*)hBx->Clone();     bBx->SetName("bBx");         bBx->Reset();
   TH1F *bBy     = (TH1F*)hBy->Clone();     bBy->SetName("bBy");         bBy->Reset();
   TH1F *bXfirst = (TH1F*)hXfirst->Clone(); bXfirst->SetName("bXfirst"); bXfirst->Reset();
   TH1F *bYfirst = (TH1F*)hYfirst->Clone(); bYfirst->SetName("bYfirst"); bYfirst->Reset();
   TH1F *bZfirst = (TH1F*)hZfirst->Clone(); bZfirst->SetName("bZfirst"); bZfirst->Reset();
   TH1F *bXlast  = (TH1F*)hXlast->Clone();  bXlast->SetName("bXlast");   bXlast->Reset();
   TH1F *bYlast  = (TH1F*)hYlast->Clone();  bYlast->SetName("bYlast");   bYlast->Reset();
   TH1F *bZlast  = (TH1F*)hZlast->Clone();  bZlast->SetName("bZlast");   bZlast->Reset();
   TH1F *bCharge = (TH1F*)hCharge->Clone(); bCharge->SetName("bCharge"); bCharge->Reset();
   TH1F *bNpoint = (TH1F*)hNpoint->Clone(); bNpoint->SetName("bNpoint"); bNpoint->Reset();
   TH1F *bValid  = (TH1F*)hValid->Clone();  bValid->SetName("bValid");   bValid->Reset();

   TH1F *bFullMatrix    =(TH1F*)hFullMatrix->Clone();    bFullMatrix->SetName("bFullMatrix");       bFullMatrix->Reset();
   TH1F *bColMatrix    = (TH1F*)hColMatrix->Clone();     bColMatrix->SetName("bColMatrix");         bColMatrix->Reset();
   TH1F *bRowMatrix    = (TH1F*)hRowMatrix->Clone();     bRowMatrix->SetName("bRowMatrix");         bRowMatrix->Reset();
   TH1F *bCellMatrix   = (TH1F*)hCellMatrix->Clone();    bCellMatrix->SetName("bCellMatrix");       bCellMatrix->Reset();
   TH1F *bFullOper     = (TH1F*)hFullOper->Clone();      bFullOper->SetName("bFullOper");           bFullOper->Reset();
   TH1F *bCellOper     = (TH1F*)hCellOper->Clone();      bCellOper->SetName("bCellOper");           bCellOper->Reset();
   TH1F *bColOper      = (TH1F*)hColOper->Clone();       bColOper->SetName("bColOper");             bColOper->Reset();
   TH1F *bRowOper      = (TH1F*)hRowOper->Clone();       bRowOper->SetName("bRowOper");             bRowOper->Reset();
   TH1F *bMatchRowOper = (TH1F*)hMatchRowOper->Clone();  bMatchRowOper->SetName("bMatchRowOper");   bMatchRowOper->Reset();
   TH1F *bMatchColOper = (TH1F*)hMatchColOper->Clone();  bMatchColOper->SetName("bMatchColOper");   bMatchColOper->Reset();
   TH1F *bRowMatOper   = (TH1F*)hRowMatOper->Clone();    bRowMatOper->SetName("bRowMatOper");       bRowMatOper->Reset();
   TH1F *bMatchDiffOper= (TH1F*)hMatchDiffOper->Clone(); bMatchDiffOper->SetName("bMatchDiffOper"); bMatchDiffOper->Reset();
   TH1F *bFullOper2    = (TH1F*)hFullOper2->Clone();     bFullOper2->SetName("bFullOper2");         bFullOper2->Reset();

   // Loop with user code on all events and fill the b histograms
   // The code below should produce identical results to the tree->Draw above

   TClonesArray *tracks = event->GetTracks();
   Int_t nev = (Int_t)tree->GetEntries();
   Int_t i, ntracks, evmod,i0,i1;
   Track *t;
   EventHeader *head;
   Int_t nbin = 0;
   for (Int_t ev=0;ev<nev;ev++) {
      nbin += tree->GetEntry(ev);
      head = event->GetHeader();
      evmod = head->GetEvtNum()%10;
      bNtrack->Fill(event->GetNtrack());
      bNseg->Fill(event->GetNseg());
      bTemp->Fill(event->GetTemperature());
      bHmean->Fill(event->GetHistogram()->GetMean());
      ntracks = event->GetNtrack();
      for(i0=0;i0<4;i0++) {
         for(i1=0;i1<4;i1++) {
            bFullMatrix->Fill(event->GetMatrix(i0,i1));
         }
         bColMatrix->Fill(event->GetMatrix(i0,0));
         bRowMatrix->Fill(event->GetMatrix(1,i0)); // done here because the matrix is square!
      }
      bCellMatrix->Fill(event->GetMatrix(2,2));
      if ( 5 < ntracks ) {
         t = (Track*)tracks->UncheckedAt(5);
         for(i0=0;i0<4;i0++) {
            for(i1=0;i1<4;i1++) {
            }
            bColOper->Fill( event->GetMatrix(i0,1) - t->GetVertex(1) );
            bRowOper->Fill( event->GetMatrix(2,i0) - t->GetVertex(2) );
         }
         for(i0=0;i0<3;i0++) {
            bMatchRowOper->Fill( event->GetMatrix(2,i0) - t->GetVertex(i0) );
            bMatchDiffOper->Fill( event->GetMatrix(i0,2) - t->GetVertex(i0) );
         }
         bCellOper->Fill( event->GetMatrix(2,1) - t->GetVertex(1) );
      }
      for (i=0;i<ntracks;i++) {
         t = (Track*)tracks->UncheckedAt(i);
         if (evmod == 0) bPx->Fill(t->GetPx());
         if (evmod == 0) bPy->Fill(t->GetPy());
         if (evmod == 0) bPz->Fill(t->GetPz());
         if (evmod == 1) bRandom->Fill(t->GetRandom());
         if (evmod == 1) bMass2->Fill(t->GetMass2());
         if (evmod == 1) bBx->Fill(t->GetBx());
         if (evmod == 1) bBy->Fill(t->GetBy());
         if (evmod == 2) bXfirst->Fill(t->GetXfirst());
         if (evmod == 2) bYfirst->Fill(t->GetYfirst());
         if (evmod == 2) bZfirst->Fill(t->GetZfirst());
         if (evmod == 3) bXlast->Fill(t->GetXlast());
         if (evmod == 3) bYlast->Fill(t->GetYlast());
         if (evmod == 3) bZlast->Fill(t->GetZlast());
         if (t->GetPx() < 0) {
            bCharge->Fill(t->GetCharge());
            bNpoint->Fill(t->GetNpoint());
            bValid->Fill(t->GetValid());
         }
         if (i<4) {
            for(i1=0;i1<3;i1++) { // 3 is the min of the 2nd dim of Matrix and Vertex
               bFullOper ->Fill( event->GetMatrix(i,i1) - t->GetVertex(i1) );
               bFullOper2->Fill( event->GetMatrix(i,i1) - t->GetVertex(i1) );
               bRowMatOper->Fill( event->GetMatrix(i,2) - t->GetVertex(i1) );
            }
            bMatchColOper->Fill( event->GetMatrix(i,2) - t->GetVertex(1) );
         }
      }
   }

   // Compare h and b histograms
   Int_t cNtrack = HistCompare(hNtrack,bNtrack);
   Int_t cNseg   = HistCompare(hNseg,bNseg);
   Int_t cTemp   = HistCompare(hTemp,bTemp);
   Int_t cHmean  = HistCompare(hHmean,bHmean);
   Int_t cPx     = HistCompare(hPx,bPx);
   Int_t cPy     = HistCompare(hPy,bPy);
   Int_t cPz     = HistCompare(hPz,bPz);
   Int_t cRandom = HistCompare(hRandom,bRandom);
   Int_t cMass2  = HistCompare(hMass2,bMass2);
   Int_t cBx     = HistCompare(hBx,bBx);
   Int_t cBy     = HistCompare(hBy,bBy);
   Int_t cXfirst = HistCompare(hXfirst,bXfirst);
   Int_t cYfirst = HistCompare(hYfirst,bYfirst);
   Int_t cZfirst = HistCompare(hZfirst,bZfirst);
   Int_t cXlast  = HistCompare(hXlast,bXlast);
   Int_t cYlast  = HistCompare(hYlast,bYlast);
   Int_t cZlast  = HistCompare(hZlast,bZlast);
   Int_t cCharge = HistCompare(hCharge,bCharge);
   Int_t cNpoint = HistCompare(hNpoint,bNpoint);
   Int_t cValid  = HistCompare(hValid,bValid);

   Int_t cFullMatrix   = HistCompare(hFullMatrix,bFullMatrix);
   Int_t cColMatrix    = HistCompare(hColMatrix,bColMatrix);
   Int_t cRowMatrix    = HistCompare(hRowMatrix,bRowMatrix);
   Int_t cCellMatrix   = HistCompare(hCellMatrix,bCellMatrix);
   Int_t cFullOper     = HistCompare(hFullOper,bFullOper);
   Int_t cCellOper     = HistCompare(hCellOper,bCellOper);
   Int_t cColOper      = HistCompare(hColOper,bColOper);
   Int_t cRowOper      = HistCompare(hRowOper,bRowOper);
   Int_t cMatchRowOper = HistCompare(hMatchRowOper,bMatchRowOper);
   Int_t cMatchColOper = HistCompare(hMatchColOper,bMatchColOper);
   Int_t cRowMatOper   = HistCompare(hRowMatOper,bRowMatOper);
   Int_t cMatchDiffOper= HistCompare(hMatchDiffOper,bMatchDiffOper);
   Int_t cFullOper2    = HistCompare(hFullOper2,bFullOper2);

   delete event;
   Event::Reset();
   ntotin += nbin;

   if (gPrintSubBench) {
      printf("Test %2dC: ",realTestNum);
      gBenchmark->Show("stress");gBenchmark->Start("stress");
      // Since we disturbed the flow (due to the double benchmark printing),
      // let's repeat the header!
      printf("Test %2d : ",realTestNum);
   }

   Bool_t OK = kTRUE;
   if (cNtrack || cNseg   || cTemp  || cHmean || cPx    || cPy     || cPz) OK = kFALSE;
   if (cRandom || cMass2  || cBx    || cBy    || cXfirst|| cYfirst || cZfirst) OK = kFALSE;
   if (cXlast  || cYlast  || cZlast || cCharge|| cNpoint|| cValid) OK = kFALSE;
   if (cFullMatrix || cColMatrix || cRowMatrix || cCellMatrix || cFullOper ) OK = kFALSE;
   if (cCellOper || cColOper || cRowOper || cMatchRowOper || cMatchColOper ) OK = kFALSE;
   if (cRowMatOper || cMatchDiffOper || cFullOper2 ) OK = kFALSE;
   if (OK) printf("OK\n");
   else    {
      printf("FAILED\n");
      printf("%-8s cNtrak =%d, cNseg  =%d, cTemp  =%d, cHmean =%d\n"," ",cNtrack,cNseg,cTemp,cHmean);
      printf("%-8s cPx    =%d, cPy    =%d, cPz    =%d, cRandom=%d\n"," ",cPx,cPy,cPz,cRandom);
      printf("%-8s cMass2 =%d, cbx    =%d, cBy    =%d, cXfirst=%d\n"," ",cMass2,cBx,cBy,cXfirst);
      printf("%-8s cYfirst=%d, cZfirst=%d, cXlast =%d, cYlast =%d\n"," ",cYfirst,cZfirst,cXlast,cYlast);
      printf("%-8s cZlast =%d, cCharge=%d, cNpoint=%d, cValid =%d\n"," ",cZlast,cCharge,cNpoint,cValid);
      printf("%-8s cFullMatrix=%d, cColMatrix=%d, cRowMatrix=%d, cCellMatrix=%d\n"," ",cFullMatrix,cColMatrix,cRowMatrix,cCellMatrix);
      printf("%-8s cFullOper=%d, cCellOper=%d, cColOper=%d, cRowOper=%d\n"," ",cFullOper,cCellOper,cColOper,cRowOper);
      printf("%-8s cMatchRowOper=%d, cMatchColOper=%d, cRowMatOper=%d, cMatchDiffOper=%d\n"," ",cMatchRowOper,cMatchColOper,cRowMatOper,cMatchDiffOper);
      printf("%-8s cFullOper2=%d\n"," ",cFullOper2);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Analyse the file Event.root generated in the last part of test8

void stress9()
{
   Bprint(9,"Analyze Event.root file of stress 8");

   gROOT->GetList()->Delete();
   TFile *hfile = new TFile("Event.root");
   TTree *tree; hfile->GetObject("T",tree);

   stress9tree(tree,9);

   // Save test9 histograms
   TFile f("stress_test9.root","recreate");
   gROOT->GetList()->Write();
   gROOT->GetList()->Delete();
   ntotout += f.GetBytesWritten();


   delete hfile;
}

////////////////////////////////////////////////////////////////////////////////
/// Make 10 Trees starting from the Event.root tree.
/// Events for which event_number%10 == 0 go to Event_0.root
/// Events for which event_number%10 == 1 go to Event_1.root
///...
/// Events for which event_number%10 == 9 go to Event_9.root

void stress10()
{
   Bprint(10,"Create 10 files starting from Event.root");

   TFile *hfile = new TFile("Event.root");
   if (hfile==0 || hfile->IsZombie()) {
      delete hfile;
      printf("FAILED\n");
      return;
   }
   TTree *tree; hfile->GetObject("T",tree);

   Event *event = 0;
   tree->SetBranchAddress("event",&event);

   // Create 10 clones of this tree
   char filename[20];
   TTree *chTree[10];
   TFile *chfile[10];
   Int_t file;
   for (file=0;file<10;file++) {
      snprintf(filename,20,"Event_%d.root",file);
      chfile[file] = new TFile(filename,"recreate");
      if (file>=5) {
         chfile[file]->SetCompressionAlgorithm(ROOT::kLZMA);
      }
      chTree[file] = (TTree*)tree->CloneTree(0);
   }

   // Fill the small trees
   Int_t nev = (Int_t)tree->GetEntries();
   Int_t evmod, nbin=0, nbout=0;
   EventHeader *head;
   for (Int_t ev=0;ev<nev;ev++) {
      nbin += tree->GetEntry(ev);
      head = event->GetHeader();
      evmod = head->GetEvtNum()%10;
      nbout += chTree[evmod]->Fill();
      event->Clear();
   }
   // save headers
   Int_t ntot = 0;
   for (file=0;file<10;file++) {
      ntot += (Int_t)chTree[file]->GetEntries();
      chfile[file]->Write();
      delete chfile[file];
   }
   delete event;
   delete hfile;
   Event::Reset();
   ntotin  += nbin;
   ntotout += nbout;

   //We compare the number of bytes read from the big file
   //with the total number of bytes written in the 10 small files
   Bool_t OK = kTRUE;
   if (nbin != nbout || nev != ntot) OK = kFALSE;
   if (OK) printf("OK\n");
   else    {
      printf("FAILED\n");
      printf("%-8s nbin=%d, nbout=%d, nev=%d, ntot=%d\n"," ",nbin,nbout,nev,ntot);
   }
   if (gPrintSubBench) { printf("Test 10 : "); gBenchmark->Show("stress");gBenchmark->Start("stress"); }
}

////////////////////////////////////////////////////////////////////////////////
/// Test chains of Trees
/// We make a TChain using the 10 files generated in test10
/// We expect the same results when analyzing the chain than
/// in the analysis of the original big file Event.root in test9.
/// Because TChain derives from TTree, we can use the same
/// analysis procedure "stress9tree"

void stress11()
{
   Bprint(11,"Test chains of Trees using the 10 files");

   gROOT->GetList()->Delete();
   TChain *chain = new TChain("T");
   char filename[20];
   Int_t file;
   for (file=0;file<10;file++) {
      snprintf(filename,20,"Event_%d.root",file);
      chain->Add(filename);
   }

   stress9tree(chain,11);

   // Save test11 histograms
   delete chain;
   TFile f("stress_test11.root","recreate");
   gROOT->GetList()->Write();
   gROOT->GetList()->Delete();
   ntotout += f.GetBytesWritten();
}

////////////////////////////////////////////////////////////////////////////////
/// Compare histograms of stress9 with stress11

void stress12(Int_t testid)
{
   if (testid == 12) Bprint(12,"Compare histograms of test 9 and 11");

   TFile f9("stress_test9.root");
   TFile f11("stress_test11.root");
   //Let's loop on all keys of second file
   //We expect to find the same keys in the original stress9 file
   TIter next(f11.GetListOfKeys());
   TKey *key;
   TH1F *h9, *h11;
   Int_t comp, ngood = 0;
   while ((key=(TKey*)next())) {
      if (strcmp(key->GetClassName(),"TH1F")) continue; //may be a TList of TStreamerInfo
      h9  = (TH1F*)f9.Get(key->GetName());
      h11 = (TH1F*)f11.Get(key->GetName());
      if (h9 == 0 || h11 == 0) continue;
      comp = HistCompare(h9,h11);
      if (comp == 0) ngood++;
   }
   ntotin += f9.GetBytesRead();
   ntotin += f11.GetBytesRead();
   Bool_t OK = kTRUE;
   if (ngood < 40) OK = kFALSE;
   if (OK) printf("OK\n");
   else    {
      printf("FAILED\n");
      printf("%-8s ngood=%d\n"," ",ngood);
   }
   if (gPrintSubBench) { printf("Test 12 : "); gBenchmark->Show("stress");gBenchmark->Start("stress"); }
}

////////////////////////////////////////////////////////////////////////////////
/// test of TChain::Merge
/// The 10 small Tree files generated in stress10 are again merged
/// into one single file.
/// Should be the same as the file generated in stress8, except
/// that events will be in a different order.
/// But global analysis histograms should be identical (checked by stress14)

void stress13()
{
   Bprint(13,"Test merging files of a chain");

   gROOT->GetList()->Delete();
   TChain *chain = new TChain("T");
   char filename[20];
   Int_t file;
   for (file=0;file<10;file++) {
      snprintf(filename,20,"Event_%d.root",file);
      chain->Add(filename);
   }

   chain->Merge("Event.root");

   Double_t chentries = chain->GetEntries();
   delete chain;

   Event::Reset();
   gROOT->GetList()->Delete();

   TFile f("Event.root");
   TTree *tree = (TTree*)f.Get("T");
   ntotin  += (Double_t)f.GetEND();
   ntotout += (Double_t)f.GetEND();

   Bool_t OK = kTRUE;
   if (chentries != tree->GetEntries()) OK = kFALSE;
   if (OK) printf("OK\n");
   else    {
      printf("FAILED\n");
   }
   if (gPrintSubBench) { printf("Test 13 : "); gBenchmark->Show("stress");gBenchmark->Start("stress"); }
}

////////////////////////////////////////////////////////////////////////////////
/// Verify that stress13 has correctly rebuild the original Event.root

void stress14()
{
   Bprint(14,"Check correct rebuilt of Event.root in test 13");

   stress12(14);
}

////////////////////////////////////////////////////////////////////////////////
/// Divert some branches to separate files

void stress15()
{
   Bprint(15,"Divert Tree branches to separate files");

   //Get old file, old tree and set top branch address
   //We want to copy only a few branches.
   TFile *oldfile = new TFile("Event.root");
   if (oldfile->IsZombie()) {
      printf("FAILED\n");
      return;
   }
   TTree *oldtree; oldfile->GetObject("T",oldtree);
   Event *event   = 0;
   oldtree->SetBranchAddress("event",&event);
   oldtree->SetBranchStatus("*",0);
   oldtree->SetBranchStatus("event",1);
   oldtree->SetBranchStatus("fNtrack",1);
   oldtree->SetBranchStatus("fNseg",1);
   oldtree->SetBranchStatus("fH",1);


   //Create a new file + a clone of old tree header. Do not copy events
   TFile *newfile = new TFile("stress_small.root","recreate");
   TTree *newtree = oldtree->CloneTree(0);

   //Divert branch fH to a separate file and copy all events
   newtree->GetBranch("fH")->SetFile("stress_fH.root");
   newtree->CopyEntries(oldtree);

   newfile->Write();
   ntotin  += oldfile->GetBytesRead();
   ntotout += newfile->GetBytesWritten();
   delete event;
   delete newfile;
   delete oldfile;
   Event::Reset();
   gROOT->GetList()->Delete();

   // Open small file, histogram fNtrack and fH
   newfile = new TFile("stress_small.root");
   newfile->GetObject("T", newtree);
   newtree->Draw("fNtrack>>hNtrack","","goff");
   newtree->Draw("fH.GetMean()>>hHmean","","goff");
   TH1 *hNtrack; newfile->GetObject("hNtrack",hNtrack);
   TH1 *hHmean; newfile->GetObject("hHmean",hHmean);
   ntotin  += newfile->GetBytesRead();

   // Open old reference file of stress9
   oldfile = new TFile("stress_test9.root");
   if (oldfile->IsZombie()) {
      printf("FAILED\n");
      return;
   }
   TH1 *bNtrack; oldfile->GetObject("bNtrack",bNtrack);
   TH1 *bHmean;  oldfile->GetObject("bHmean",bHmean);
   Int_t cNtrack = HistCompare(hNtrack,bNtrack);
   Int_t cHmean  = HistCompare(hHmean, bHmean);
   delete newfile;
   delete oldfile;
   Event::Reset();
   gROOT->GetList()->Delete();

   Bool_t OK = kTRUE;
   if (cNtrack || cHmean) OK = kFALSE;
   if (OK) printf("OK\n");
   else    {
      printf("FAILED\n");
      printf("%-8s cNtrack=%d, cHmean=%d\n"," ",cNtrack,cHmean);
   }
   if (gPrintSubBench) { printf("Test 15 : "); gBenchmark->Show("stress");gBenchmark->Start("stress"); }
}

void stress16()
{
// Prototype trigger simulation for the LHCb experiment
// This test nested loops with the interpreter.
// Expected to run fast with the compiler, slow with the interpreter.
// This code is extracted from an original macro by Hans Dijkstra (LHCb)
// The program generates histograms and profile histograms.
// A canvas with subpads containing the results is sent to Postscript.
// We check graphics results by counting the number of lines in the ps file.

   Bprint(16,"CINT test (3 nested loops) with LHCb trigger");

   const int nbuf    = 153;    // buffer size
   const int nlev    = 4;      // number of trigger levels
   const int nstep   = 50000;  // number of steps
   const int itt[4]  = { 1000, 4000, 40000, 400000 }; // time needed per trigger
   const float a[4]  = { 0.25, 0.04, 0.25, 0 };       // acceptance/trigger (last always 0)

   int   i, il, istep, itim[192], itrig[192], it, im, ipass;
   float dead, sum[10];

   // create histogram and array of profile histograms
   TCanvas *c = new TCanvas("laten","latency simulation",700,600);
   gROOT->LoadClass("TPostScript","Postscript");
   TPostScript ps("stress_lhcb.ps",112);
   gRandom->SetSeed(65539);
   TFile f("stress_lhcb.root", "recreate");
   TH1F *pipe = new TH1F("pipe", "free in pipeline", nbuf+1, -0.5, nbuf+0.5);
   pipe->SetLineColor(2);
   pipe->SetFillColor(2);

   TProfile *hp[nlev+1];
   TProfile::Approximate();
   for (i = 0; i <= nlev; i++) {
      char s[64];
      snprintf(s,64, "buf%d", i);
      hp[i] = new TProfile(s, "in buffers", 1000, 0,nstep, -1., 1000.);
      hp[i]->SetLineColor(2);
   }

   dead   = 0;
   sum[0] = nbuf;
   for (i = 1; i <= nlev; i++) sum[i] = 0;
   for (i = 0; i < nbuf; i++) { itrig[i] = 0; itim[i] = 0; }

   for (istep = 0; istep < nstep; istep++) {
      // evaluate status of buffer
      pipe->Fill(sum[0]);
      if ((istep+1)%10 == 0) {
         for (i = 0; i <= nlev; i++)
            hp[i]->Fill((float)istep, sum[i], 1.);
      }

      ipass = 0;
      for (i = 0; i < nbuf; i++) {
         it = itrig[i];
         if (it >= 1) {
            // add 25 ns to all times
            itim[i] += 25;
            im = itim[i];
            // level decisions
            for (il = 0; il < nlev; il++) {
               if (it == il+1 && im > itt[il]) {
                  if (gRandom->Rndm() > a[il]) {
                     itrig[i] = -1;
                     sum[0]++;
                     sum[il+1]--;
                  } else {
                     itrig[i]++;
                     sum[il+1]--;
                     sum[il+2]++;
                  }
               }
            }
         } else if (ipass == 0) {
            itrig[i] = 1;
            itim[i]  = 25;
            sum[0]--;
            sum[1]++;
            ipass++;
         }
      }
      if (ipass == 0) dead++;
   }
//   Float_t deadTime = 100.*dead/nstep;

   // View results in the canvas and make the Postscript file

   c->Divide(2,3);
   c->cd(1); pipe->Draw();
   c->cd(2); hp[0]->Draw();
   c->cd(3); hp[1]->Draw();
   c->cd(4); hp[2]->Draw();
   c->cd(5); hp[3]->Draw();
   c->cd(6); hp[4]->Draw();
   ps.Close();

   f.Write();
   ntotout += f.GetBytesWritten();

   // Check length of Postscript file
   FILE *fp = fopen("stress_lhcb.ps","r");
   char line[260];
   Int_t nlines = 0;
   Int_t nlinesGood = 2121;
   Bool_t counting = kFALSE;
   while (fgets(line,255,fp)) {
      if (counting) nlines++;
      if (strstr(line,"%%EndProlog")) counting = kTRUE;
   }
   fclose(fp);
   delete c;
   Bool_t OK = kTRUE;
   if (nlines < nlinesGood-100 || nlines > nlinesGood+100) OK = kFALSE;
   if (OK) printf("OK\n");
   else    {
      printf("FAILED\n");
      printf("%-8s nlines in stress_lhcb.ps file = %d\n"," ",nlines);
   }
   if (gPrintSubBench) { printf("Test 16 : "); gBenchmark->Show("stress");gBenchmark->Start("stress"); }
}

////////////////////////////////////////////////////////////////////////////////
/// Test mkdir returnExistingDirectory outside and within a Root file
/// Create some directories, ensure they point to the expected places

void stress17()
{
   Bprint(17,"Test mkdir");

   // check mkdir functions as expected in TDirectory
   TDirectory *free_motherdir = new TDirectory("free_motherdir", "free_motherdir");
   TDirectory *free_daughterdir = free_motherdir->mkdir("free_daughterdir");
   TDirectory *free_daughter2 = free_motherdir->mkdir("free_daughterdir");
   TDirectory *free_daughtersame = free_motherdir->mkdir("free_daughterdir", "", true);

   // check mkdir functions as expected inside a file
   TFile f("stress.root","update");
   TDirectory *motherdir = f.mkdir("motherdir");
   TDirectory *daughterdir = motherdir->mkdir("daughterdir");
   TDirectory *daughtersame = motherdir->mkdir("daughterdir", "", true);

   Bool_t OK = kTRUE;
   if (daughtersame != daughterdir || free_daughter2 == free_daughterdir || free_daughtersame != free_daughterdir) OK = kFALSE;
   if (OK) printf("OK\n");
   else    {
      printf("FAILED\n");
      printf("%-8s free_daughterdir=%p, free_daughter2=%p, free_daughtersame=%p, daughterdir=%p, daughtersame=%p \n"," ",free_daughterdir,free_daughter2,free_daughtersame,daughterdir,daughtersame);
   }
   if (gPrintSubBench) { printf("Test 17 : "); gBenchmark->Show("stress");gBenchmark->Start("stress"); }
}

void cleanup()
{
   gSystem->Unlink("Event.root");
   gSystem->Unlink("Event_0.root");
   gSystem->Unlink("Event_1.root");
   gSystem->Unlink("Event_2.root");
   gSystem->Unlink("Event_3.root");
   gSystem->Unlink("Event_4.root");
   gSystem->Unlink("Event_5.root");
   gSystem->Unlink("Event_6.root");
   gSystem->Unlink("Event_7.root");
   gSystem->Unlink("Event_8.root");
   gSystem->Unlink("Event_9.root");
   gSystem->Unlink("stress.ps");
   gSystem->Unlink("stress.root");
   gSystem->Unlink("stress_fH.root");
   gSystem->Unlink("stress_lhcb.ps");
   gSystem->Unlink("stress_lhcb.root");
   gSystem->Unlink("stress_small.root");
   gSystem->Unlink("stress_test9.root");
   gSystem->Unlink("stress_test11.root");
}
