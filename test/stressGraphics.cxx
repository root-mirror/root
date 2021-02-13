// @(#)root/test:$name:  $:$id: stressGraphics.cxx,v 1.0 exp $
// Author: O.Couet

//
//    ROOT  Graphics test suite and benchmarks.
//
// The suite of programs below tests many elements of the graphics classes
//
// The test can be run as a standalone program or with the interpreter.
//
// To run as a standalone program:
//
//    make stressGraphics
//    stressGraphics
//
// To get a short help:
//
//    stressGraphics -h
//
// To run interactively, do:
//
//    root -b
//    root [0] .L stressGraphics.cxx
//    root [1] stressGraphics()
//
// or:
//
//    root -b
//    root [0] .x  stressGraphics.cxx

#include <stdlib.h>
#include <Riostream.h>
#include <time.h>
#include <TString.h>
#include <TROOT.h>
#include <TError.h>
#include <TRandom.h>
#include <TBenchmark.h>
#include <TSystem.h>
#include <TApplication.h>
#include <TDatime.h>
#include <TFile.h>
#include <TF1.h>
#include "TF2.h"
#include <TF3.h>
#include <TH2.h>
#include <TH2Poly.h>
#include <TNtuple.h>
#include <TKey.h>
#include <TProfile.h>
#include "TString.h"

#include <TStyle.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TFrame.h>
#include <TPostScript.h>
#include <TPDF.h>
#include <TLine.h>
#include <TMarker.h>
#include <TPolyLine.h>
#include <TLatex.h>
#include <TMathText.h>
#include <TLegend.h>
#include <TEllipse.h>
#include <TCurlyArc.h>
#include <TArc.h>
#include <TPaveText.h>
#include <TPaveStats.h>
#include <TPaveLabel.h>
#include <TRatioPlot.h>
#include <TGaxis.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TGraphBentErrors.h>
#include <TMultiGraph.h>
#include <TGraph2D.h>
#include <TParallelCoord.h>
#include <TImage.h>
#include <TMath.h>


void     stressGraphics (Int_t verbose);
Int_t    StatusPrint    (TString &filename, Int_t id, const TString &title, Int_t res, Int_t ref, Int_t err);
Int_t    AnalysePS      (const TString &filename);
Int_t    FileSize       (char *filename);
TCanvas *StartTest      (Int_t w, Int_t h);
void     TestReport1    (TCanvas *C, const TString &title, Int_t IPS=0);
void     TestReport2    (Int_t IPS=0);
void     DoCcode        (TCanvas *C);
TString  stime(time_t* t, bool utc = false, bool display_time_zone = true);


// Tests functions.
void     clonepad       ();
void     earth          ();
void     feynman        ();
void     hbars          ();
void     itbf           ();
void     kerning        ();
void     labels1        ();
void     ntuple1        ();
void     options2d1     ();
void     options2d2     ();
void     options2d3     ();
void     options2d4     ();
void     options2d5     ();
void     parallelcoord  ();
void     patterns       ();
void     quarks         ();
void     ratioplot      ();
void     statfitparam   ();
void     tellipse       ();
void     tgaxis1        ();
void     tgaxis2        ();
void     tgaxis3        ();
void     tgaxis4        ();
void     tgaxis5        ();
void     tgraph1        ();
void     tgraph2        ();
void     tgraph3        ();
void     tgraph4        ();
void     tgraph2d1      ();
void     tgraph2d2      ();
void     tgraph2d3      ();
void     th2poly        ();
void     timage         ();
void     tlatex1        ();
void     tlatex2        ();
void     tlatex3        ();
void     tlatex4        ();
void     tlatex5        ();
void     tline          ();
void     tmarker        ();
void     tmathtext      ();
void     tmultigraph1   ();
void     tmultigraph2   ();
void     tpolyline      ();
void     transpad       ();
void     ttext1         ();
void     ttext2         ();
void     waves          ();
void     zoomfit        ();
void     zoomtf1        ();

// Auxiliary functions
void     patterns_box   (Int_t pat, Double_t x1, Double_t y1, Double_t x2, Double_t  y2);
Double_t interference   (Double_t *x, Double_t *par);
Double_t result         (Double_t *x, Double_t *par);
void     cleanup        ();

// Global variables.
Int_t     gVerbose;
Int_t     gTestNum;
Int_t     gTestsFailed;
Int_t     gPS1RefNb[60];
Int_t     gPS1ErrNb[60];
Int_t     gPDFRefNb[60];
Int_t     gPDFErrNb[60];
Int_t     gGIFRefNb[60];
Int_t     gGIFErrNb[60];
Int_t     gJPGRefNb[60];
Int_t     gJPGErrNb[60];
Int_t     gPNGRefNb[60];
Int_t     gPNGErrNb[60];
Int_t     gPS2RefNb[60];
Int_t     gPS2ErrNb[60];
Bool_t    gOptionR;
Bool_t    gOptionK;
TH2F     *gH2;
TFile    *gHsimple;
TFile    *gCernstaff;
char      gCfile[16];
char      outfile[16];
char      gLine[80];


#ifndef __CLING__
////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
   TString opt;
   Int_t verbose = 0;
   if (argc > 1) verbose = atoi(argv[1]);
   opt = argv[1];

   if (argc > 2) {
      for (int i = 2; i<argc; i++) opt += argv[i];
   }

   if (opt.Contains("-h")) {
      printf("Usage: stressGraphics [-h] [-r] [-k]\n");
      printf("Options:\n");
      printf("  -r : Generate de reference output.\n");
      printf("       Redirect the output in the file \"stressGraphics.ref\"\n");
      printf("       to redefine the reference file.\n");
      printf("\n");
      printf("  -k : Keep the output files even for passed tests.\n");
      printf("       By default output files for passed tests are deleted.\n");
      printf("\n");
      printf("  -h : Print usage\n");
      printf("\n");
      printf("  Any other option is ignored.\n");
      return 0;
   }

   gROOT->SetBatch();
   TApplication theApp("App", &argc, argv);
   gBenchmark = new TBenchmark();

   if (opt.Contains("-r")) {
      gOptionR = kTRUE;
   } else {
      gOptionR = kFALSE;
   }

   if (opt.Contains("-k")) {
      gOptionK = kTRUE;
   } else {
      gOptionK = kFALSE;
   }

   stressGraphics(verbose);

   cleanup();
   return 0;
}
#else
// Unlinking causes the file system to re-use inodes. This confuses clang:
// "Ummm I know the content of the file with that inode number already",
// even though the deleted file has been replaced by a new one (which happens
// to re-use the old inode number).
#define ClingWorkAroundDeletedSourceFile 1
#endif


////////////////////////////////////////////////////////////////////////////////
/// Run all graphics stress tests.

void stressGraphics(Int_t verbose = 0)
{
   gErrorIgnoreLevel = 9999;
   gROOT->SetStyle("Classic");

   // Check if $ROOTSYS/tutorials/hsimple.root exists
   gHsimple = new TFile("$(ROOTSYS)/tutorials/hsimple.root");
   if (gHsimple->IsZombie()) {
      delete gHsimple;
      gHsimple = new TFile("hsimple.root");
      if (gHsimple->IsZombie()) {
         delete gHsimple;
         printf("Create $(ROOTSYS)/tutorials/hsimple.root\n");
         gROOT->Macro("$(ROOTSYS)/tutorials/hsimple.C");
         gHsimple = new TFile("$(ROOTSYS)/tutorials/hsimple.root");
         if (gHsimple->IsZombie()) {
            delete gHsimple;
            printf("Could not create $(ROOTSYS)/tutorials/hsimple.root\n");
            return;
         }
      }
   }

   // Check if cernstaff.root exists
   gCernstaff = new TFile("cernstaff.root");
   if (gCernstaff->IsZombie()) {
      delete gCernstaff;
      gCernstaff = new TFile("cernstaff.root");
      if (gCernstaff->IsZombie()) {
         delete gCernstaff;
         printf("Create ./cernstaff.root\n");
         gROOT->Macro("$(ROOTSYS)/tutorials/tree/cernbuild.C(0,0)");
         gCernstaff = new TFile("cernstaff.root");
         if (gCernstaff->IsZombie()) {
            delete gCernstaff;
            printf("Could not create ./cernstaff.root\n");
            return;
         }
      }
   }

   gErrorIgnoreLevel = 0;

   // Read the reference file "stressGraphics.ref"
#ifdef R__HAS_CLOUDFLARE_ZLIB
   FILE *sg = fopen("stressGraphics_builtinzlib.ref","r");
#else
   FILE *sg = fopen("stressGraphics.ref","r");
#endif
   if (!sg) {
      printf("Could not open stressGraphics.ref/stressGraphics_builtinzlib.ref\n");
      return;
   }
   char line[160];
   Int_t i = -1;
   while (fgets(line,160,sg)) {
      if (i>=0) {
         sscanf(&line[7]  ,"%d",&gPS1RefNb[i]);
         sscanf(&line[18] ,"%d",&gPS1ErrNb[i]);
         sscanf(&line[28] ,"%d",&gPDFRefNb[i]);
         sscanf(&line[38] ,"%d",&gPDFErrNb[i]);
         sscanf(&line[48] ,"%d",&gGIFRefNb[i]);
         sscanf(&line[58] ,"%d",&gGIFErrNb[i]);
         sscanf(&line[68] ,"%d",&gJPGRefNb[i]);
         sscanf(&line[78] ,"%d",&gJPGErrNb[i]);
         sscanf(&line[88] ,"%d",&gPNGRefNb[i]);
         sscanf(&line[98] ,"%d",&gPNGErrNb[i]);
         sscanf(&line[107],"%d",&gPS2RefNb[i]);
         sscanf(&line[118],"%d",&gPS2ErrNb[i]);
      }
      i++;
   }
   fclose(sg);

   gRandom->SetSeed(65539);

   if (gOptionR) {
      std::cout << "Test#   PS1Ref#   PS1Err#   PDFRef#   PDFErr#   GIFRef#   GIFErr#   JPGRef#   JPGErr#   PNGRef#   PNGErr#   PS2Ref#   PS2Err#" <<std::endl;
   } else {
      std::cout << "**********************************************************************" <<std::endl;
      std::cout << "*  Starting  Graphics - S T R E S S suite                            *" <<std::endl;
      std::cout << "**********************************************************************" <<std::endl;
   }

   gVerbose     = verbose;
   gTestNum     = 0;
   gTestsFailed = 0;

   gBenchmark->Start("stressGraphics");

   if (!gOptionR) {
      std::cout << "*  Starting Basic Graphics - S T R E S S                             *" <<std::endl;
      std::cout << "**********************************************************************" <<std::endl;
   }
   tline         ();
   tmarker       ();
   tpolyline     ();
   patterns      ();
   ttext1        ();
   ttext2        ();
   tlatex1       ();
   tlatex2       ();
   tlatex3       ();
   tlatex4       ();
   tlatex5       ();
   kerning       ();
   itbf          ();
   tmathtext     ();
   transpad      ();
   statfitparam  ();
   if (!gOptionR) {
      std::cout << "**********************************************************************" <<std::endl;
      std::cout << "*  Starting High Level 2D Primitives - S T R E S S                   *" <<std::endl;
      std::cout << "**********************************************************************" <<std::endl;
   }
   tgaxis1       ();
   tgaxis2       ();
   tgaxis3       ();
   tgaxis4       ();
   tgaxis5       ();
   labels1       ();
   tellipse      ();
   feynman       ();
   ratioplot     ();
   tgraph1       ();
   tgraph2       ();
   tgraph3       ();
   tgraph4       ();
   tmultigraph1  ();
   tmultigraph2  ();
   waves         ();
   if (!gOptionR) {
      std::cout << "**********************************************************************" <<std::endl;
      std::cout << "*  Starting High Level 3D Primitives - S T R E S S                   *" <<std::endl;
      std::cout << "**********************************************************************" <<std::endl;
   }
   options2d1    ();
   options2d2    ();
   options2d3    ();
   options2d4    ();
   options2d5    ();
   earth         ();
   tgraph2d1     ();
   tgraph2d2     ();
   tgraph2d3     ();
   if (!gOptionR) {
      std::cout << "**********************************************************************" <<std::endl;
      std::cout << "*  Starting complex drawing and TPad - S T R E S S                   *" <<std::endl;
      std::cout << "**********************************************************************" <<std::endl;
   }
   ntuple1       ();
   quarks        ();
   timage        ();
   zoomtf1       ();
   zoomfit       ();
   parallelcoord ();
   clonepad      ();
   hbars         ();
   th2poly       ();
   if (!gOptionR) {
      std::cout << "**********************************************************************" <<std::endl;
      if (!gTestsFailed) {
         std::cout << "*  All the tests passed. :-)" <<std::endl;
      } else {
         if (gTestsFailed>1) std::cout << "*  " << gTestsFailed <<" tests failed. :-(" <<std::endl;
         else                std::cout << "*  " << gTestsFailed <<" test failed. :-(" <<std::endl;
      }
      std::cout << "**********************************************************************" <<std::endl;

      gBenchmark->Stop("stressGraphics");

      //Print table with results
      Bool_t UNIX = strcmp(gSystem->GetName(), "Unix") == 0;
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
         if (!os) printf("*  SYS: Windows 95\n");
         else     printf("*  SYS: %s %s \n",os,gSystem->Getenv("PROCESSOR_IDENTIFIER"));
      }

      printf("**********************************************************************\n");
      printf("* ");
      gBenchmark->Print("stressGraphics");

      Double_t ct = gBenchmark->GetCpuTime("stressGraphics");
      //normalize at 860 rootmarks on pcbrun4
      const Double_t rootmarks = 860*(47.12/ct);

      printf("**********************************************************************\n");
      printf("*  ROOTMARKS =%6.1f   *  Root %-8s  %d/%d\n",rootmarks,gROOT->GetVersion(),
             gROOT->GetVersionDate(),gROOT->GetVersionTime());
      printf("**********************************************************************\n");
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Print test program number and its title

Int_t StatusPrint(TString &filename, Int_t id, const TString &title,
                  Int_t res, Int_t ref, Int_t err)
{
   if (!gOptionR) {
      if (id>0) {
         sprintf(gLine,"Test %2d: %s",id,title.Data());
      } else {
         sprintf(gLine,"       %s",title.Data());
      }

      const Int_t nch = strlen(gLine);
      if (TMath::Abs(res-ref)<=err) {
         std::cout << gLine;
         for (Int_t i = nch; i < 67; i++) std::cout << ".";
         std::cout << " OK" << std::endl;
#ifndef ClingWorkAroundDeletedSourceFile
         if (!gOptionK) gSystem->Unlink(filename.Data());
#endif
      } else {
         std::cout << gLine;
         Int_t ndots = 60;
         Int_t w = 3;
         if (gTestNum < 10) { ndots++; w--;}
         for (Int_t i = nch; i < ndots; i++) std::cout << ".";
         std::cout << std::setw(w) << gTestNum << " FAILED" << std::endl;
         std::cout << "         Result    = "  << res << std::endl;
         std::cout << "         Reference = "  << ref << std::endl;
         std::cout << "         Error     = "  << TMath::Abs(res-ref)
                                          << " (was " << err << ")"<< std::endl;
         gTestsFailed++;
         return 1;
      }
   } else {
      if (id>0)  printf("%5d%10d%10d",id,res,err);
      if (id==0) printf("%10d%10d",res,err);
      if (id<0)  printf("%10d%10d\n",res,err);
   }
   return 0;
}


////////////////////////////////////////////////////////////////////////////////
/// Return the size of filename

Int_t FileSize (char *filename)
{
   FileStat_t fs;
   if (!gSystem->GetPathInfo(filename, fs)) {
      return (Int_t)fs.fSize;
   } else {
      return 0;
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Analyse the PS file "filename" and return the number of character in the
/// meaningful part of the file. The variable part (date etc..) are not
/// counted.

Int_t AnalysePS(const TString &filename)
{
   Bool_t counting = kFALSE;
   Int_t count = 0;

   FILE *fp;
   Int_t status;
   if ((fp=fopen(filename.Data(), "r"))==NULL) {
      printf("ERROR1 : File can not open !..\n");
      return 0;
   }

   char *line = new char[251];
   TString l;
   while((status=fscanf(fp, "%s", line)) != EOF) {
      l = line;
      if (l.Contains("%!PS-Adobe"))  counting = kFALSE;
      if (l.Contains("%%EndProlog")) counting = kTRUE;
      if (counting) count = count+l.Length();
   }
   if (gVerbose==1) printf(">>>>>>>>> Number of characters found in %s: %d\n",filename.Data(),count);
   fclose(fp);
   delete [] line;
   return count;
}


////////////////////////////////////////////////////////////////////////////////
/// Start Test:
/// Open the TCanvas C and set the acceptable error (number of characters)

TCanvas *StartTest(Int_t w, Int_t h)
{
   gTestNum++;
   gStyle->Reset();
   TCanvas *old = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("C");
   if (old) {
      if (old->IsOnHeap()) delete old;
   }
   TCanvas *C = new TCanvas("C","C",0,0,w,h);
   return C;
}


////////////////////////////////////////////////////////////////////////////////
/// Report 1:
/// Draw the canvas generate as PostScript, count the number of characters in
/// the PS file and compare the result with the reference value.

void TestReport1(TCanvas *C, const TString &title, Int_t IPS)
{
   gErrorIgnoreLevel = 9999;
   sprintf(outfile,"sg1_%2.2d.ps",gTestNum);

   TPostScript *ps1 = new TPostScript(outfile, 111);
   C->Draw();
   ps1->Close();
   TString psfile = outfile;
   if (IPS) {
      StatusPrint(psfile,  gTestNum, title, FileSize(outfile) ,
                                            gPS1RefNb[gTestNum-1],
                                            gPS1ErrNb[gTestNum-1]);
   } else {
      StatusPrint(psfile,  gTestNum, title, AnalysePS(outfile) ,
                                            gPS1RefNb[gTestNum-1],
                                            gPS1ErrNb[gTestNum-1]);
   }

   sprintf(outfile,"sg%2.2d.pdf",gTestNum);
   C->cd(0);
   TPDF *pdf = new TPDF(outfile,111);
   C->Draw();
   pdf->Close();
   TString pdffile = outfile;
   StatusPrint(pdffile, 0, "  PDF output", FileSize(outfile),
                                           gPDFRefNb[gTestNum-1],
                                           gPDFErrNb[gTestNum-1]);

   sprintf(outfile,"sg%2.2d.gif",gTestNum);
   C->cd(0);
   C->SaveAs(outfile);
   TString giffile = outfile;
   StatusPrint(giffile, 0, "  GIF output", FileSize(outfile),
                                           gGIFRefNb[gTestNum-1],
                                           gGIFErrNb[gTestNum-1]);

   sprintf(outfile,"sg%2.2d.jpg",gTestNum);
   C->cd(0);
   C->SaveAs(outfile);
   TString jpgfile = outfile;
   StatusPrint(jpgfile, 0, "  JPG output", FileSize(outfile),
                                           gJPGRefNb[gTestNum-1],
                                           gJPGErrNb[gTestNum-1]);

   sprintf(outfile,"sg%2.2d.png",gTestNum);
   C->cd(0);
   C->SaveAs(outfile);
   TString pngfile = outfile;
   StatusPrint(pngfile, 0, "  PNG output", FileSize(outfile),
                                           gPNGRefNb[gTestNum-1],
                                           gPNGErrNb[gTestNum-1]);

   gErrorIgnoreLevel = 0;

   return;
}


////////////////////////////////////////////////////////////////////////////////
/// Generate the C code corresponding to the canvas C.

void DoCcode(TCanvas *C)
{
   gErrorIgnoreLevel = 9999;

   sprintf(gCfile,"sg%2.2d.C",gTestNum);

   if (C) {
      C->SaveAs(gCfile);
      delete C;
      C = 0;
   }

   gErrorIgnoreLevel = 0;
   return;
}


////////////////////////////////////////////////////////////////////////////////
/// Report 2:
/// Draw the canvas generate as .C, generate the corresponding PostScript
/// file (using gPad), count the number of characters in it and compare the
/// result with the reference value.

void TestReport2(Int_t IPS)
{
   sprintf(outfile,"sg2_%2.2d.ps",gTestNum);

   gErrorIgnoreLevel = 9999;
   sprintf(gCfile,".x sg%2.2d.C",gTestNum);
   gROOT->ProcessLine(gCfile);
   gPad->SaveAs(outfile);
   gErrorIgnoreLevel = 0;
   Int_t i;

   TString psfile = outfile;
   if (IPS) {
      i = StatusPrint(psfile,-1, "  C file result", FileSize(outfile),
                                                    gPS2RefNb[gTestNum-1],
                                                    gPS2ErrNb[gTestNum-1]);
   } else {
      i = StatusPrint(psfile,-1, "  C file result", AnalysePS(outfile),
                                                    gPS2RefNb[gTestNum-1],
                                                    gPS2ErrNb[gTestNum-1]);
   }

   sprintf(gCfile,"sg%2.2d.C",gTestNum);
#ifndef ClingWorkAroundDeletedSourceFile
   if (!gOptionK && !i) gSystem->Unlink(gCfile);
#endif
   return;
}


////////////////////////////////////////////////////////////////////////////////
/// Test TLine.

void tline()
{
   TCanvas *C = StartTest(800,800);

   TLine *l1 = new TLine(0.1,0.1,0.9,0.1);
   l1->SetLineColor(1); l1->SetLineWidth(1) ; l1->SetLineStyle(1) ; l1->Draw();
   TLine *l2 = new TLine(0.1,0.2,0.9,0.2);
   l2->SetLineColor(2); l2->SetLineWidth(2) ; l2->SetLineStyle(2) ; l2->Draw();
   TLine *l3 = new TLine(0.1,0.3,0.9,0.3);
   l3->SetLineColor(3); l3->SetLineWidth(3) ; l3->SetLineStyle(3) ; l3->Draw();
   TLine *l4 = new TLine(0.1,0.4,0.9,0.4);
   l4->SetLineColor(4); l4->SetLineWidth(4) ; l4->SetLineStyle(4) ; l4->Draw();
   TLine *l5 = new TLine(0.1,0.5,0.9,0.5);
   l5->SetLineColor(5); l5->SetLineWidth(5) ; l5->SetLineStyle(5) ; l5->Draw();
   TLine *l6 = new TLine(0.1,0.6,0.9,0.6);
   l6->SetLineColor(6); l6->SetLineWidth(6) ; l6->SetLineStyle(6) ; l6->Draw();
   TLine *l7 = new TLine(0.1,0.7,0.9,0.7);
   l7->SetLineColor(7); l7->SetLineWidth(7) ; l7->SetLineStyle(7) ; l7->Draw();
   TLine *l8 = new TLine(0.1,0.8,0.9,0.8);
   l8->SetLineColor(8); l8->SetLineWidth(8) ; l8->SetLineStyle(8) ; l8->Draw();
   TLine *l9 = new TLine(0.1,0.9,0.9,0.9);
   l9->SetLineColor(9); l9->SetLineWidth(9) ; l9->SetLineStyle(9) ; l9->Draw();

   TestReport1(C, "TLine");
   DoCcode(C);
   TestReport2();
}


////////////////////////////////////////////////////////////////////////////////
/// Test TMarker

void tmarker()
{
   TCanvas *C = StartTest(500,200);

   TMarker m;
   m.DisplayMarkerTypes();

   TestReport1(C, "TMarker");
   DoCcode(C);
   TestReport2();
}


////////////////////////////////////////////////////////////////////////////////
/// Test TPolyLine

void tpolyline()
{
   TCanvas *C = StartTest(700,500);

   C->Range(0,30,11,650);
   Double_t x[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
   Double_t y[10] = {200, 300, 300, 200, 200, 100, 10, 580, 10, 600};
   TPolyLine *p = new TPolyLine(10,x,y);
   p->SetLineWidth(3);
   p->SetLineColor(2);
   p->Draw("F");
   p->Draw("");

   TestReport1(C, "TPolyLine");
   DoCcode(C);
   TestReport2();
}


////////////////////////////////////////////////////////////////////////////////
/// Test Patterns

void patterns()
{
   TCanvas *C = StartTest(700,900);

   C->Range(0,0,1,1);
   C->SetBorderSize(2);
   C->SetFrameFillColor(0);
   Double_t bh = 0.059;
   Double_t db = 0.01;
   Double_t y  = 0.995;
   Int_t i,j=3001;
   for (i=1; i<=5; i++) {
      patterns_box(j++, 0.01, y-bh, 0.19, y);
      patterns_box(j++, 0.21, y-bh, 0.39, y);
      patterns_box(j++, 0.41, y-bh, 0.59, y);
      patterns_box(j++, 0.61, y-bh, 0.79, y);
      patterns_box(j++, 0.81, y-bh, 0.99, y);
      y = y-bh-db;
   }
   y = y-3*db;
   gStyle->SetHatchesSpacing(2.0);
   gStyle->SetHatchesLineWidth(3);
   Int_t j1 = 3144;
   Int_t j2 = 3305;
   Int_t j3 = 3350;
   Int_t j4 = 3490;
   Int_t j5 = 3609;
   for (i=1; i<=9; i++) {
      if (i==6) {j2 += 10; j3 += 1; j4 += 1; j5 += 10;}
      if (i==5) {j4 -= 10; j5 -= 1;}
      patterns_box(j1, 0.01, y-bh, 0.19, y);
      patterns_box(j2, 0.21, y-bh, 0.39, y);
      patterns_box(j3, 0.41, y-bh, 0.59, y);
      patterns_box(j4, 0.61, y-bh, 0.79, y);
      patterns_box(j5, 0.81, y-bh, 0.99, y);
      j1 += 100;
      j2 += 10;
      j3 += 1;
      j4 -= 9;
      j5 += 9;
      y = y-bh-db;
   }

   TestReport1(C, "Fill patterns");
   DoCcode(C);
   TestReport2();
}


////////////////////////////////////////////////////////////////////////////////
/// Auxiliary function used by "patterns"

void patterns_box(Int_t pat, Double_t x1, Double_t y1, Double_t x2, Double_t  y2)
{
   TBox b;
   b.SetFillColor(1);
   b.SetFillStyle(pat);
   b.DrawBox(x1,y1,x2,y2);
   b.SetFillStyle(0);
   b.DrawBox(x1,y1,x2,y2);
   b.SetFillColor(0);
   b.SetFillStyle(1000);
   Double_t dx = (x2-x1)/3;
   Double_t dy = (y2-y1)/3;
   Double_t h  = (y2-y1)/3;
   b.DrawBox(x1+dx, y1+dy, x2-dx, y2-dy);
   b.SetFillStyle(0);
   b.DrawBox(x1+dx, y1+dy, x2-dx, y2-dy);
   TLatex l;
   l.SetTextAlign(22);
   l.SetTextSize(h);
   l.DrawLatex((x1+x2)/2, (y1+y2)/2, Form("%d",pat));
}


////////////////////////////////////////////////////////////////////////////////
/// 1st TText test.

void ttext1()
{
   TCanvas *C = StartTest(900,500);

   C->Range(0,0,1,1);
   C->SetBorderSize(2);
   C->SetFrameFillColor(0);
   TLine *lv = new TLine(0.5,0.0,0.5,1.0);
   lv->Draw();
   for (float s=0.1; s<1.0 ; s+=0.1) {
      TLine *lh = new TLine(0.,s,1.,s);
      lh->Draw();
   }
   TText *tex1b = new TText(0.02,0.4,"jgabcdefhiklmnopqrstuvwxyz_{}");
   // 161 is not a valid text font. This tests if the protection against
   // invalid text font is working.
   tex1b->SetTextFont(161);
   tex1b->SetTextColor(2);
   tex1b->SetTextAngle(0);
   tex1b->SetTextAlign(11);
   tex1b->SetTextSize(0.1);
   tex1b->Draw();
   TText *tex1 = new TText(0.5,0.1,"j0al {&`ag}_:^)Jj");
   tex1->SetTextFont(41);
   tex1->SetTextColor(2);
   tex1->SetTextAngle( 0);
   tex1->SetTextAlign(21);
   tex1->SetTextSize(0.15);
   tex1->Draw();
   TText *tex2 = new TText(0.5,0.5,"j0Al {&`ag}_:^)Jj");
   tex2->SetTextColor(3);
   tex2->SetTextFont(21);
   tex2->SetTextAlign(21);
   tex2->SetTextSize(0.1);
   tex2->Draw();
   TText *tex3 = new TText(0.5,0.3,"j0Al {&`ag}_:^)Jj");
   tex3->SetTextColor(4);
   tex3->SetTextFont(31);
   tex3->SetTextAlign(31);
   tex3->SetTextSize(0.1);
   tex3->Draw();
   TText *tex4 = new TText(0.5,0.8,"j0Al {&`ag}_:^)Jj");
   tex4->SetTextColor(5);
   tex4->SetTextFont(71);
   tex4->SetTextAlign(22);
   tex4->SetTextSize(0.07);
   tex4->Draw();
   TText *tex5 = new TText(0.5,0.7,"13 j0Al {&`ag}_:^)Jj");
   tex5->SetTextColor(6);
   tex5->SetTextFont(51);
   tex5->SetTextAlign(13);
   tex5->SetTextSize(0.1);
   tex5->Draw();

   TestReport1(C, "TText 1 (Text attributes)");
   DoCcode(C);
   TestReport2();
}


////////////////////////////////////////////////////////////////////////////////
/// 2nd TText test. A very long text string.

void ttext2()
{
   TCanvas *C = StartTest(600,600);

   TText t(0.001,0.5,"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789");
   t.SetTextFont(42); t.SetTextSize(0.02);
   t.Draw();

   TestReport1(C, "TText 2 (A very long text string)");
   DoCcode(C);
   TestReport2();
}


////////////////////////////////////////////////////////////////////////////////
/// 1st TLatex test.

void tlatex1()
{
   TCanvas *C = StartTest(600,700);

   TLatex l;
   l.SetTextAlign(12);
   l.SetTextSize(0.04);
   l.DrawLatex(0.1,0.9,"1) C(x) = d #sqrt{#frac{2}{#lambdaD}}  #int^{x}_{0}cos(#frac{#pi}{2}t^{2})dt");
   l.DrawLatex(0.1,0.7,"2) C(x) = d #sqrt{#frac{2}{#lambdaD}}  #int^{x}cos(#frac{#pi}{2}t^{2})dt");
   l.DrawLatex(0.1,0.5,"3) R = |A|^{2} = #frac{1}{2}#left(#[]{#frac{1}{2}+C(V)}^{2}+#[]{#frac{1}{2}+S(V)}^{2}#right)");
   l.DrawLatex(0.1,0.3,"4) F(t) = #sum_{i=-#infty}^{#infty}A(i)cos#[]{#frac{i}{t+i}}");
   l.DrawLatex(0.1,0.1,"5) {}_{3}^{7}Li");

   TestReport1(C, "TLatex 1");
   DoCcode(C);
   TestReport2();
}


////////////////////////////////////////////////////////////////////////////////
/// 2nd TLatex test.

void tlatex2()
{
   TCanvas *C = StartTest(700,500);

   TLatex l;
   l.SetTextAlign(23);
   l.SetTextSize(0.1);
   l.DrawLatex(0.5,0.95,"e^{+}e^{-}#rightarrowZ^{0}#rightarrowI#bar{I}, q#bar{q}");
   l.DrawLatex(0.5,0.75,"|#vec{a}#bullet#vec{b}|=#Sigmaa^{i}_{jk}+b^{bj}_{i}");
   l.DrawLatex(0.5,0.5,"i(#partial_{#mu}#bar{#psi}#gamma^{#mu}+m#bar{#psi})=0#Leftrightarrow(#Box+m^{2})#psi=0");
   l.DrawLatex(0.5,0.3,"L_{em}=eJ^{#mu}_{em}A_{#mu} , ^{}J^{#mu}_{em}=#bar{I}#gamma_{#mu}I , M^{j}_{i}=#SigmaA_{#alpha}#tau^{#alphaj}_{i}");

   TestReport1(C, "TLatex 2");
   DoCcode(C);
   TestReport2();
}


////////////////////////////////////////////////////////////////////////////////
/// 3rd TLatex test.

void tlatex3()
{
   TCanvas *C = StartTest(700,500);

   TPaveText pt(.05,.1,.95,.8);
   pt.AddText("#frac{2s}{#pi#alpha^{2}}  #frac{d#sigma}{dcos#theta} (e^{+}e^{-}#rightarrowf#bar{f} ) = \
#left| #frac{1}{1 - #Delta#alpha} #right|^{2} (1+cos^{2}#theta)");
   pt.AddText("+ 4 Re #left{ #frac{2}{1 - #Delta#alpha} #chi(s) #[]{#hat{g}_{#nu}^{e}#hat{g}_{#nu}^{f} \
(1 + cos^{2}#theta) + 2 #hat{g}_{a}^{e}#hat{g}_{a}^{f} cos#theta) } #right}");
   pt.AddText("+ 16#left|#chi(s)#right|^{2}\
#left[(#hat{g}_{a}^{e}^{2}+#hat{g}_{v}^{e}^{2})\
(#hat{g}_{a}^{f}^{2} + #hat{g}_{v}^{f}^{2})(1+cos^{2}#theta) \
+ 8 #hat{g}_{a}^{e} #hat{g}_{a}^{f} #hat{g}_{v}^{e} #hat{g}_{v}^{f}cos#theta#right] ");
   pt.SetLabel("Born equation");
   pt.Draw();

   TestReport1(C, "TLatex 3 (TLatex in TPaveText)");
   DoCcode(C);
   TestReport2();
}


////////////////////////////////////////////////////////////////////////////////
/// 4th TLatex test.

void tlatex4()
{
   TCanvas *C = StartTest(600,700);

   TLatex l;
   l.SetTextSize(0.03);
   l.SetTextAlign(22);
   l.DrawLatex(0.165, 0.95, "Lower case");
   l.DrawLatex(0.495, 0.95, "Upper case");
   l.DrawLatex(0.825, 0.95, "Variations");
   l.SetTextAlign(12);
   float y, x1, x2;
   y = 0.90; x1 = 0.07; x2 = x1+0.2;
                 l.DrawLatex(x1, y, "alpha : ")   ; l.DrawLatex(x2, y, "#alpha");
   y -= 0.0375 ; l.DrawLatex(x1, y, "beta : ")    ; l.DrawLatex(x2, y, "#beta");
   y -= 0.0375 ; l.DrawLatex(x1, y, "gamma : ")   ; l.DrawLatex(x2, y, "#gamma");
   y -= 0.0375 ; l.DrawLatex(x1, y, "delta : ")   ; l.DrawLatex(x2, y, "#delta");
   y -= 0.0375 ; l.DrawLatex(x1, y, "epsilon : ") ; l.DrawLatex(x2, y, "#epsilon");
   y -= 0.0375 ; l.DrawLatex(x1, y, "zeta : ")    ; l.DrawLatex(x2, y, "#zeta");
   y -= 0.0375 ; l.DrawLatex(x1, y, "eta : ")     ; l.DrawLatex(x2, y, "#eta");
   y -= 0.0375 ; l.DrawLatex(x1, y, "theta : ")   ; l.DrawLatex(x2, y, "#theta");
   y -= 0.0375 ; l.DrawLatex(x1, y, "iota : ")    ; l.DrawLatex(x2, y, "#iota");
   y -= 0.0375 ; l.DrawLatex(x1, y, "kappa : ")   ; l.DrawLatex(x2, y, "#kappa");
   y -= 0.0375 ; l.DrawLatex(x1, y, "lambda : ")  ; l.DrawLatex(x2, y, "#lambda");
   y -= 0.0375 ; l.DrawLatex(x1, y, "mu : ")      ; l.DrawLatex(x2, y, "#mu");
   y -= 0.0375 ; l.DrawLatex(x1, y, "nu : ")      ; l.DrawLatex(x2, y, "#nu");
   y -= 0.0375 ; l.DrawLatex(x1, y, "xi : ")      ; l.DrawLatex(x2, y, "#xi");
   y -= 0.0375 ; l.DrawLatex(x1, y, "omicron : ") ; l.DrawLatex(x2, y, "#omicron");
   y -= 0.0375 ; l.DrawLatex(x1, y, "pi : ")      ; l.DrawLatex(x2, y, "#pi");
   y -= 0.0375 ; l.DrawLatex(x1, y, "rho : ")     ; l.DrawLatex(x2, y, "#rho");
   y -= 0.0375 ; l.DrawLatex(x1, y, "sigma : ")   ; l.DrawLatex(x2, y, "#sigma");
   y -= 0.0375 ; l.DrawLatex(x1, y, "tau : ")     ; l.DrawLatex(x2, y, "#tau");
   y -= 0.0375 ; l.DrawLatex(x1, y, "upsilon : ") ; l.DrawLatex(x2, y, "#upsilon");
   y -= 0.0375 ; l.DrawLatex(x1, y, "phi : ")     ; l.DrawLatex(x2, y, "#phi");
   y -= 0.0375 ; l.DrawLatex(x1, y, "chi : ")     ; l.DrawLatex(x2, y, "#chi");
   y -= 0.0375 ; l.DrawLatex(x1, y, "psi : ")     ; l.DrawLatex(x2, y, "#psi");
   y -= 0.0375 ; l.DrawLatex(x1, y, "omega : ")   ; l.DrawLatex(x2, y, "#omega");
   y = 0.90; x1 = 0.40; x2 = x1+0.2;
                 l.DrawLatex(x1, y, "Alpha : ")   ; l.DrawLatex(x2, y, "#Alpha");
   y -= 0.0375 ; l.DrawLatex(x1, y, "Beta : ")    ; l.DrawLatex(x2, y, "#Beta");
   y -= 0.0375 ; l.DrawLatex(x1, y, "Gamma : ")   ; l.DrawLatex(x2, y, "#Gamma");
   y -= 0.0375 ; l.DrawLatex(x1, y, "Delta : ")   ; l.DrawLatex(x2, y, "#Delta");
   y -= 0.0375 ; l.DrawLatex(x1, y, "Epsilon : ") ; l.DrawLatex(x2, y, "#Epsilon");
   y -= 0.0375 ; l.DrawLatex(x1, y, "Zeta : ")    ; l.DrawLatex(x2, y, "#Zeta");
   y -= 0.0375 ; l.DrawLatex(x1, y, "Eta : ")     ; l.DrawLatex(x2, y, "#Eta");
   y -= 0.0375 ; l.DrawLatex(x1, y, "Theta : ")   ; l.DrawLatex(x2, y, "#Theta");
   y -= 0.0375 ; l.DrawLatex(x1, y, "Iota : ")    ; l.DrawLatex(x2, y, "#Iota");
   y -= 0.0375 ; l.DrawLatex(x1, y, "Kappa : ")   ; l.DrawLatex(x2, y, "#Kappa");
   y -= 0.0375 ; l.DrawLatex(x1, y, "Lambda : ")  ; l.DrawLatex(x2, y, "#Lambda");
   y -= 0.0375 ; l.DrawLatex(x1, y, "Mu : ")      ; l.DrawLatex(x2, y, "#Mu");
   y -= 0.0375 ; l.DrawLatex(x1, y, "Nu : ")      ; l.DrawLatex(x2, y, "#Nu");
   y -= 0.0375 ; l.DrawLatex(x1, y, "Xi : ")      ; l.DrawLatex(x2, y, "#Xi");
   y -= 0.0375 ; l.DrawLatex(x1, y, "Omicron : ") ; l.DrawLatex(x2, y, "#Omicron");
   y -= 0.0375 ; l.DrawLatex(x1, y, "Pi : ")      ; l.DrawLatex(x2, y, "#Pi");
   y -= 0.0375 ; l.DrawLatex(x1, y, "Rho : ")     ; l.DrawLatex(x2, y, "#Rho");
   y -= 0.0375 ; l.DrawLatex(x1, y, "Sigma : ")   ; l.DrawLatex(x2, y, "#Sigma");
   y -= 0.0375 ; l.DrawLatex(x1, y, "Tau : ")     ; l.DrawLatex(x2, y, "#Tau");
   y -= 0.0375 ; l.DrawLatex(x1, y, "Upsilon : ") ; l.DrawLatex(x2, y, "#Upsilon");
   y -= 0.0375 ; l.DrawLatex(x1, y, "Phi : ")     ; l.DrawLatex(x2, y, "#Phi");
   y -= 0.0375 ; l.DrawLatex(x1, y, "Chi : ")     ; l.DrawLatex(x2, y, "#Chi");
   y -= 0.0375 ; l.DrawLatex(x1, y, "Psi : ")     ; l.DrawLatex(x2, y, "#Psi");
   y -= 0.0375 ; l.DrawLatex(x1, y, "Omega : ")   ; l.DrawLatex(x2, y, "#Omega");
   x1 = 0.73; x2 = x1+0.2;
   y = 0.7500 ; l.DrawLatex(x1, y, "varepsilon : ") ; l.DrawLatex(x2, y, "#varepsilon");
   y = 0.6375 ; l.DrawLatex(x1, y, "vartheta : ")   ; l.DrawLatex(x2, y, "#vartheta");
   y = 0.2625 ; l.DrawLatex(x1, y, "varsigma : ")   ; l.DrawLatex(x2, y, "#varsigma");
   y = 0.1875 ; l.DrawLatex(x1, y, "varUpsilon : ") ; l.DrawLatex(x2, y, "#varUpsilon");
   y = 0.1500 ; l.DrawLatex(x1, y, "varphi : ")     ; l.DrawLatex(x2, y, "#varphi");
   y = 0.0375 ; l.DrawLatex(x1, y, "varomega : ")   ; l.DrawLatex(x2, y, "#varomega");

   TestReport1(C, "TLatex 4 (Greek letters)");
   DoCcode(C);
   TestReport2();
}


////////////////////////////////////////////////////////////////////////////////
/// 5th TLatex test.

void tlatex5()
{
   TCanvas *C = StartTest(600,600);

   TLatex l;
   l.SetTextSize(0.03);
   l.SetTextAlign(12);
   float y, step, x1, x2;
   y = 0.96; step = 0.0465; x1 = 0.02; x2 = x1+0.04;
               l.DrawLatex(x1, y, "#club")           ; l.DrawText(x2, y, "#club");
   y -= step ; l.DrawLatex(x1, y, "#voidn")          ; l.DrawText(x2, y, "#voidn");
   y -= step ; l.DrawLatex(x1, y, "#leq")            ; l.DrawText(x2, y, "#leq");
   y -= step ; l.DrawLatex(x1, y, "#approx")         ; l.DrawText(x2, y, "#approx");
   y -= step ; l.DrawLatex(x1, y, "#in")             ; l.DrawText(x2, y, "#in");
   y -= step ; l.DrawLatex(x1, y, "#supset")         ; l.DrawText(x2, y, "#supset");
   y -= step ; l.DrawLatex(x1, y, "#cap")            ; l.DrawText(x2, y, "#cap");
   y -= step ; l.DrawLatex(x1, y, "#ocopyright")     ; l.DrawText(x2, y, "#ocopyright");
   y -= step ; l.DrawLatex(x1, y, "#trademark")      ; l.DrawText(x2, y, "#trademark");
   y -= step ; l.DrawLatex(x1, y, "#times")          ; l.DrawText(x2, y, "#times");
   y -= step ; l.DrawLatex(x1, y, "#bullet")         ; l.DrawText(x2, y, "#bullet");
   y -= step ; l.DrawLatex(x1, y, "#voidb")          ; l.DrawText(x2, y, "#voidb");
   y -= step ; l.DrawLatex(x1, y, "#doublequote")    ; l.DrawText(x2, y, "#doublequote");
   y -= step ; l.DrawLatex(x1, y, "#lbar")           ; l.DrawText(x2, y, "#lbar");
   y -= step ; l.DrawLatex(x1, y, "#arcbottom")      ; l.DrawText(x2, y, "#arcbottom");
   y -= step ; l.DrawLatex(x1, y, "#downarrow")      ; l.DrawText(x2, y, "#downarrow");
   y -= step ; l.DrawLatex(x1, y, "#leftrightarrow") ; l.DrawText(x2, y, "#leftrightarrow");
   y -= step ; l.DrawLatex(x1, y, "#Downarrow")      ; l.DrawText(x2, y, "#Downarrow");
   y -= step ; l.DrawLatex(x1, y, "#Leftrightarrow") ; l.DrawText(x2, y, "#Leftrightarrow");
   y -= step ; l.DrawLatex(x1, y, "#void8")          ; l.DrawText(x2, y, "#void8");
   y -= step ; l.DrawLatex(x1, y, "#hbar")           ; l.DrawText(x2, y, "#hbar");
   y = 0.96; step = 0.0465; x1 = 0.27; x2 = x1+0.04;
               l.DrawLatex(x1, y, "#diamond")        ; l.DrawText(x2, y, "#diamond");
   y -= step ; l.DrawLatex(x1, y, "#aleph")          ; l.DrawText(x2, y, "#aleph");
   y -= step ; l.DrawLatex(x1, y, "#geq")            ; l.DrawText(x2, y, "#geq");
   y -= step ; l.DrawLatex(x1, y, "#neq")            ; l.DrawText(x2, y, "#neq");
   y -= step ; l.DrawLatex(x1, y, "#notin")          ; l.DrawText(x2, y, "#notin");
   y -= step ; l.DrawLatex(x1, y, "#subseteq")       ; l.DrawText(x2, y, "#subseteq");
   y -= step ; l.DrawLatex(x1, y, "#cup")            ; l.DrawText(x2, y, "#cup");
   y -= step ; l.DrawLatex(x1, y, "#copyright")      ; l.DrawText(x2, y, "#copyright");
   y -= step ; l.DrawLatex(x1, y, "#void3")          ; l.DrawText(x2, y, "#void3");
   y -= step ; l.DrawLatex(x1, y, "#divide")         ; l.DrawText(x2, y, "#divide");
   y -= step ; l.DrawLatex(x1, y, "#circ")           ; l.DrawText(x2, y, "#circ");
   y -= step ; l.DrawLatex(x1, y, "#infty")          ; l.DrawText(x2, y, "#infty");
   y -= step ; l.DrawLatex(x1, y, "#angle")          ; l.DrawText(x2, y, "#angle");
   y -= step ; l.DrawLatex(x1, y, "#cbar")           ; l.DrawText(x2, y, "#cbar");
   y -= step ; l.DrawLatex(x1, y, "#arctop")         ; l.DrawText(x2, y, "#arctop");
   y -= step ; l.DrawLatex(x1, y, "#leftarrow")      ; l.DrawText(x2, y, "#leftarrow");
   y -= step ; l.DrawLatex(x1, y, "#otimes")         ; l.DrawText(x2, y, "#otimes");
   y -= step ; l.DrawLatex(x1, y, "#Leftarrow")      ; l.DrawText(x2, y, "#Leftarrow");
   y -= step ; l.DrawLatex(x1, y, "#prod")           ; l.DrawText(x2, y, "#prod");
   y -= step ; l.DrawLatex(x1, y, "#Box")            ; l.DrawText(x2, y, "#Box");
   y -= step ; l.DrawLatex(x1, y, "#parallel")       ; l.DrawText(x2, y, "#parallel");
   y = 0.96; step = 0.0465; x1 = 0.52; x2 = x1+0.04;
               l.DrawLatex(x1, y, "#heart")          ; l.DrawText(x2, y, "#heart");
   y -= step ; l.DrawLatex(x1, y, "#Jgothic")        ; l.DrawText(x2, y, "#Jgothic");
   y -= step ; l.DrawLatex(x1, y, "#LT")             ; l.DrawText(x2, y, "#LT");
   y -= step ; l.DrawLatex(x1, y, "#equiv")          ; l.DrawText(x2, y, "#equiv");
   y -= step ; l.DrawLatex(x1, y, "#subset")         ; l.DrawText(x2, y, "#subset");
   y -= step ; l.DrawLatex(x1, y, "#supseteq")       ; l.DrawText(x2, y, "#supseteq");
   y -= step ; l.DrawLatex(x1, y, "#wedge")          ; l.DrawText(x2, y, "#wedge");
   y -= step ; l.DrawLatex(x1, y, "#oright")         ; l.DrawText(x2, y, "#oright");
   y -= step ; l.DrawLatex(x1, y, "#AA")             ; l.DrawText(x2, y, "#AA");
   y -= step ; l.DrawLatex(x1, y, "#pm")             ; l.DrawText(x2, y, "#pm");
   y -= step ; l.DrawLatex(x1, y, "#3dots")          ; l.DrawText(x2, y, "#3dots");
   y -= step ; l.DrawLatex(x1, y, "#nabla")          ; l.DrawText(x2, y, "#nabla");
   y -= step ; l.DrawLatex(x1, y, "#downleftarrow")  ; l.DrawText(x2, y, "#downleftarrow");
   y -= step ; l.DrawLatex(x1, y, "#topbar")         ; l.DrawText(x2, y, "#topbar");
   y -= step ; l.DrawLatex(x1, y, "#arcbar")         ; l.DrawText(x2, y, "#arcbar");
   y -= step ; l.DrawLatex(x1, y, "#uparrow")        ; l.DrawText(x2, y, "#uparrow");
   y -= step ; l.DrawLatex(x1, y, "#oplus")          ; l.DrawText(x2, y, "#oplus");
   y -= step ; l.DrawLatex(x1, y, "#Uparrow")        ; l.DrawText(x2, y, "#Uparrow");
   y -= step ; l.DrawLatex(x1, y-0.01, "#sum")       ; l.DrawText(x2, y, "#sum");
   y -= step ; l.DrawLatex(x1, y, "#perp")           ; l.DrawText(x2, y, "#perp");
   y = 0.96; step = 0.0465; x1 = 0.77; x2 = x1+0.04;
               l.DrawLatex(x1, y, "#spade")          ; l.DrawText(x2, y, "#spade");
   y -= step ; l.DrawLatex(x1, y, "#Rgothic")        ; l.DrawText(x2, y, "#Rgothic");
   y -= step ; l.DrawLatex(x1, y, "#GT")             ; l.DrawText(x2, y, "#GT");
   y -= step ; l.DrawLatex(x1, y, "#propto")         ; l.DrawText(x2, y, "#propto");
   y -= step ; l.DrawLatex(x1, y, "#notsubset")      ; l.DrawText(x2, y, "#notsubset");
   y -= step ; l.DrawLatex(x1, y, "#oslash")         ; l.DrawText(x2, y, "#oslash");
   y -= step ; l.DrawLatex(x1, y, "#vee")            ; l.DrawText(x2, y, "#vee");
   y -= step ; l.DrawLatex(x1, y, "#void1")          ; l.DrawText(x2, y, "#void1");
   y -= step ; l.DrawLatex(x1, y, "#aa")             ; l.DrawText(x2, y, "#aa");
   y -= step ; l.DrawLatex(x1, y, "#/")              ; l.DrawText(x2, y, "#/");
   y -= step ; l.DrawLatex(x1, y, "#upoint")         ; l.DrawText(x2, y, "#upoint");
   y -= step ; l.DrawLatex(x1, y, "#partial")        ; l.DrawText(x2, y, "#partial");
   y -= step ; l.DrawLatex(x1, y, "#corner")         ; l.DrawText(x2, y, "#corner");
   y -= step ; l.DrawLatex(x1, y, "#ltbar")          ; l.DrawText(x2, y, "#ltbar");
   y -= step ; l.DrawLatex(x1, y, "#bottombar")      ; l.DrawText(x2, y, "#bottombar");
   y -= step ; l.DrawLatex(x1, y, "#rightarrow")     ; l.DrawText(x2, y, "#rightarrow");
   y -= step ; l.DrawLatex(x1, y, "#surd")           ; l.DrawText(x2, y, "#surd");
   y -= step ; l.DrawLatex(x1, y, "#Rightarrow")     ; l.DrawText(x2, y, "#Rightarrow");
   y -= step ; l.DrawLatex(x1, y-0.015, "#int")      ; l.DrawText(x2, y, "#int");
   y -= step ; l.DrawLatex(x1, y, "#odot")           ; l.DrawText(x2, y, "#odot");

   TestReport1(C, "TLatex 5 (Mathematical Symbols)");
   DoCcode(C);
   TestReport2();
}


////////////////////////////////////////////////////////////////////////////////
/// Text kerning.

void kerning()
{
   TCanvas *C = StartTest(1000, 700);

   for (Int_t i = 0;i < 25;i++) {
      TLine *ln = new TLine(0, 0.04 * (i - 0.2), 1, 0.04 * (i - 0.2));
      ln->Draw();
      Float_t sz = 0.0016 * i;
      TLatex *l = new TLatex(0.10, 0.04 * i, "AVAVAVAVAVAVAVAVAVAVAVAVAVAVAVAVAVAVA#color[2]{X}");
      l->SetTextSize(sz);
      l->Draw();
      TLatex *l1 = new TLatex(0.05, 0.04 * i, Form("%g", sz));
      l1->SetTextSize(0.02);
      l1->Draw();
   }

   TestReport1(C, "Text kerning");
   DoCcode(C);
   TestReport2();
}


////////////////////////////////////////////////////////////////////////////////
/// TLatex commands #kern, #lower, #it and #bf

void itbf()
{
   TCanvas *C = StartTest(700, 500);

   gStyle->SetTextFont(132);

   (new TLatex(0.01, 0.9, "Positive k#kern[0.3]{e}#kern[0.3]{r}#kern[0.3]{n}#kern[0.3]{i}#kern[0.3]{n}#kern[0.3]{g} with #^{}kern[0.3]"))->Draw();
   (new TLatex(0.01, 0.7, "Negative k#kern[-0.3]{e}#kern[-0.3]{r}#kern[-0.3]{n}#kern[-0.3]{i}#kern[-0.3]{n}#kern[-0.3]{g} with #^{}kern[-0.3]"))->Draw();
   (new TLatex(0.01, 0.5, "Vertical a#lower[0.2]{d}#lower[0.4]{j}#lower[0.1]{u}#lower[-0.1]{s}#lower[-0.3]{t}#lower[-0.4]{m}#lower[-0.2]{e}#lower[0.1]{n}t with #^{}lower[-0.4...+0.4]"))->Draw();
   (new TLatex(0.01, 0.3, "Font styles: #^{}bf{#bf{bold}}, #^{}it{#it{italic}}, #^{}bf{#^{}it{#bf{#it{bold italic}}}}, #^{}bf{#^{}bf{#bf{#bf{unbold}}}}"))->Draw();
   (new TLatex(0.01, 0.1, "Font styles: abc#alpha#beta#gamma, #^{}it{#it{abc#alpha#beta#gamma}}, #^{}it{#^{}it{#it{#it{abc#alpha#beta#gamma}}}}"))->Draw();

   TestReport1(C, "TLatex commands #kern, #lower, #it and #bf");
   DoCcode(C);
   TestReport2();
}


////////////////////////////////////////////////////////////////////////////////
/// TMatText test

void tmathtext()
{
   TCanvas *C = StartTest(700, 500);

   TMathText l;
   l.SetTextAlign(23);
   l.SetTextSize(0.06);
   l.DrawMathText(0.50, 1.000, "\\prod_{j\\ge0} \\left(\\sum_{k\\ge0} a_{jk}z^k\\right) = \\sum_{n\\ge0} z^n \\left(\\sum_{k_0,k_1,\\ldots\\ge0\\atop k_0+k_1+\\cdots=n} a_{0k_0}a_{1k_1} \\cdots \\right)");
   l.DrawMathText(0.50, 0.800, "W_{\\delta_1\\rho_1\\sigma_2}^{3\\beta} = U_{\\delta_1\\rho_1\\sigma_2}^{3\\beta} + {1\\over 8\\pi^2} \\int_{\\alpha_1}^{\\alpha_2} d\\alpha_2^\\prime \\left[ {U_{\\delta_1\\rho_1}^{2\\beta} - \\alpha_2^\\prime U_{\\rho_1\\sigma_2}^{1\\beta} \\over U_{\\rho_1\\sigma_2}^{0\\beta}} \\right]");
   l.DrawMathText(0.50, 0.600, "d\\Gamma = {1\\over 2m_A} \\left( \\prod_f {d^3p_f\\over (2\\pi)^3} {1\\over 2E_f} \\right) \\left| \\mathscr{M} \\left(m_A - \\left\\{p_f\\right\\} \\right) \\right|^2 (2\\pi)^4 \\delta^{(4)} \\left(p_A - \\sum p_f \\right)");
   l.DrawMathText(0.50, 0.425, "4\\mathrm{Re}\\left\\{{2\\over 1-\\Delta\\alpha} \\chi(s) \\left[ \\^{g}_\\nu^e \\^{g}_\\nu^f (1 + \\cos^2\\theta) + \\^{g}_a^e \\^{g}_a^f \\cos\\theta \\right] \\right\\}");
   l.DrawMathText(0.50, 0.330, "p(n) = {1\\over\\pi\\sqrt{2}} \\sum_{k = 1}^\\infty \\sqrt{k} A_k(n) {d\\over dn} {\\sinh \\left\\{ {\\pi\\over k} \\sqrt{2\\over 3} \\sqrt{n - {1\\over 24}} \\right\\} \\over \\sqrt{n - {1\\over 24}}}");
   l.DrawMathText(0.13, 0.150, "{(\\ell+1)C_{\\ell}^{TE} \\over 2\\pi}");
   l.DrawMathText(0.27, 0.110, "\\mathbb{N} \\subset \\mathbb{R}");
   l.DrawMathText(0.63, 0.100, "\\hbox{RHIC スピン物理 Нью-Йорк}");

   TestReport1(C, "TMathText",1);
   DoCcode(C);
   TestReport2(1);
}


////////////////////////////////////////////////////////////////////////////////
/// Transparent pad.

void transpad()
{
   TCanvas *C = StartTest(700,500);

   TPad *pad1 = new TPad("pad1","",0,0,1,1);
   TPad *pad2 = new TPad("pad2","",0,0,1,1);
   pad2->SetFillStyle(4000); //will be transparent
   pad1->Draw();
   pad1->cd();

   TH1F *ht1 = new TH1F("ht1","ht1",100,-3,3);
   TH1F *ht2 = new TH1F("ht2","ht2",100,-3,3);
   TRandom r;
   for (Int_t i=0;i<100000;i++) {
      Double_t x1 = r.Gaus(-1,0.5);
      Double_t x2 = r.Gaus(1,1.5);
      if (i <1000) ht1->Fill(x1);
      ht2->Fill(x2);
   }
   ht1->Draw();
   pad1->Update(); //this will force the generation of the "stats" box
   TPaveStats *ps1 = (TPaveStats*)ht1->GetListOfFunctions()->FindObject("stats");
   ps1->SetX1NDC(0.4); ps1->SetX2NDC(0.6);
   pad1->Modified();
   C->cd();

   //compute the pad range with suitable margins
   Double_t ymin = 0;
   Double_t ymax = 2000;
   Double_t dy = (ymax-ymin)/0.8; //10 per cent margins top and bottom
   Double_t xmin = -3;
   Double_t xmax = 3;
   Double_t dx = (xmax-xmin)/0.8; //10 per cent margins left and right
   pad2->Range(xmin-0.1*dx,ymin-0.1*dy,xmax+0.1*dx,ymax+0.1*dy);
   pad2->Draw();
   pad2->cd();
   ht2->SetLineColor(kRed);
   ht2->Draw("][sames");
   pad2->Update();
   TPaveStats *ps2 = (TPaveStats*)ht2->GetListOfFunctions()->FindObject("stats");
   ps2->SetX1NDC(0.65); ps2->SetX2NDC(0.85);
   ps2->SetTextColor(kRed);

   // draw axis on the right side of the pad
   TGaxis *axis = new TGaxis(xmax,ymin,xmax,ymax,ymin,ymax,50510,"+L");
   axis->SetLabelColor(kRed);
   axis->Draw();

   TestReport1(C, "Transparent pad");
   DoCcode(C);
   TestReport2();
}


////////////////////////////////////////////////////////////////////////////////
/// Stat and fit parameters with errors.

void statfitparam ()
{
   TCanvas *C = StartTest(800,500);

   C->Divide(3,2);
   gStyle->SetOptFit(1111);
   gStyle->SetOptStat(111111);
   gStyle->SetStatW(0.43);
   gStyle->SetStatH(0.35);

   TH1 *hsf1 = new TH1F("hsf1","hsf1", 2,0.,1.);
   TH1 *hsf2 = new TH1F("hsf2","hsf2", 2,0.,1.);
   TH1 *hsf3 = new TH1F("hsf3","hsf3", 2,0.,1.);
   TH1 *hsf4 = new TH1F("hsf4","hsf4", 2,0.,1.);
   TH1 *hsf5 = new TH1F("hsf5","hsf5", 2,0.,1.);

   C->cd(1);
   hsf1->SetBinContent (1, 5.3E5); hsf1->SetBinError (1, 0.9);
   hsf1->SetBinContent (2, 5.3E5); hsf1->SetBinError (2, 0.1);
   hsf1->Fit("pol0","Q");

   C->cd(2);
   hsf2->SetBinContent (1, 5.0E15); hsf2->SetBinError (1, 4.9E15);
   hsf2->SetBinContent (2, 5.0E15); hsf2->SetBinError (2, 4.9E11);
   hsf2->Fit("pol0","Q");

   C->cd(3);
   hsf3->SetBinContent (1, 5.0E-15); hsf3->SetBinError (1, 4.9E-15);
   hsf3->SetBinContent (2, 5.0E-15); hsf3->SetBinError (2, 4.9E-11);
   hsf3->Fit("pol0","Q");

   C->cd(4);
   hsf4->SetBinContent (1, 5); hsf4->SetBinError (1, 3);
   hsf4->SetBinContent (2, 5); hsf4->SetBinError (2, 1);
   hsf4->Fit("pol0","Q");

   C->cd(5);
   hsf5->SetBinContent (1, 5.3); hsf5->SetBinError (1, 0.9);
   hsf5->SetBinContent (2, 5.3); hsf5->SetBinError (2, 0.1);
   hsf5->Fit("pol0","Q");

   C->cd(6);
   TPaveText *pt = new TPaveText(0.02,0.2,0.98,0.8,"brNDC");
   pt->SetFillColor(18);
   pt->SetTextAlign(12);
   pt->AddText("This example test all the possible cases");
   pt->AddText("handled by THistPainter::GetBestFormat.");
   pt->AddText("This method returns the best format to");
   pt->AddText("paint the fit parameters errors.");
   pt->Draw();

   TestReport1(C, "Stat and fit parameters with errors");
   DoCcode(C);
   TestReport2();
}

////////////////////////////////////////////////////////////////////////////////
/// 1st TGaxis test.

void tgaxis1()
{
   TCanvas *C = StartTest(700,500);

   C->Range(-10,-1,10,1);
   TGaxis *axis1 = new TGaxis(-4.5,-0.2,5.5,-0.2,-6,8,510,"");
   axis1->SetName("axis1");
   axis1->Draw();
   TGaxis *axis2 = new TGaxis(-4.5,0.2,5.5,0.2,0.001,10000,510,"G");
   axis2->SetName("axis2");
   axis2->Draw();
   TGaxis *axis3 = new TGaxis(-9,-0.8,-9,0.8,-8,8,50510,"");
   axis3->SetName("axis3");
   axis3->SetTitle("axis3");
   axis3->SetTitleOffset(0.5);
   axis3->Draw();
   TGaxis *axis4 = new TGaxis(-7,-0.8,-7,0.8,1,10000,50510,"G");
   axis4->SetName("axis4");
   axis4->SetTitle("axis4");
   axis4->Draw();
   TGaxis *axis5 = new TGaxis(-4.5,-0.6,5.5,-0.6,1.2,1.32,80506,"-+");
   axis5->SetName("axis5");
   axis5->SetLabelSize(0.03);
   axis5->SetTextFont(72);
   axis5->Draw();
   TGaxis *axis6 = new TGaxis(-4.5,0.5,5.5,0.5,100,900,50510,"-");
   axis6->SetName("axis6");
   axis6->Draw();
   TGaxis *axis6a = new TGaxis(-5.5,0.85,5.5,0.85,0,4.3e-6,510,"");
   axis6a->SetName("axis6a");
   axis6a->Draw();
   TGaxis *axis7 = new TGaxis(8,-0.8,8,0.8,0,9000,50510,"+L");
   axis7->SetName("axis7");
   axis7->Draw();
   TGaxis *axis8 = new TGaxis(6.5,0.8,6.499,-0.8,0,90,50510,"-");
   axis8->SetName("axis8");
   axis8->Draw();

   TestReport1(C, "TGaxis 1");
   DoCcode(C);
   TestReport2();
}


////////////////////////////////////////////////////////////////////////////////
/// 2nd TGaxis test.

void tgaxis2()
{
   TCanvas *C = StartTest(600,700);

   C->Range(-10,-1,10,1);
   TGaxis *axis1 = new TGaxis(-5,-0.2,6,-0.2,-6,8,510,"");
   axis1->SetName("axis1");
   axis1->Draw();
   TGaxis *axis2 = new TGaxis(-5,0.2,6,0.2,0.001,10000,510,"G");
   axis2->SetName("axis2");
   axis2->Draw();
   TGaxis *axis3 = new TGaxis(-9,-0.8,-9,0.8,-8,8,50510,"");
   axis3->SetName("axis3");
   axis3->Draw();
   TGaxis *axis4 = new TGaxis(-7,-0.8,-7,0.8,1,10000,50510,"G");
   axis4->SetName("axis4");
   axis4->Draw();
   TGaxis *axis5 = new TGaxis(-5,-0.6,6,-0.6,1.2,1.32,80506,"-+");
   axis5->SetName("axis5");
   axis5->SetLabelSize(0.03);
   axis5->SetTextFont(72);
   axis5->SetLabelOffset(0.025);
   axis5->Draw();
   TGaxis *axis6 = new TGaxis(-5,0.6,6,0.6,100,900,50510,"-");
   axis6->SetName("axis6");
   axis6->Draw();
   TGaxis *axis7 = new TGaxis(8,-0.8,8,0.8,0,9000,50510,"+L");
   axis7->SetName("axis7");
   axis7->SetLabelOffset(0.01);
   axis7->Draw();

   TestReport1(C, "TGaxis 2");
   DoCcode(C);
   TestReport2();
}


////////////////////////////////////////////////////////////////////////////////
/// 3rd TGaxis test.

void tgaxis3()
{
   TCanvas *C = StartTest(700,900);

   time_t script_time;
   script_time = time(0);
   script_time = 3600*(int)(script_time/3600);
   gStyle->SetTimeOffset(script_time);
   C->Divide(1,3);
   C->SetFillColor(28);
   int i;
   gStyle->SetTitleH(0.08);
   float noise;
   TH1F *ht = new TH1F("ht","Love at first sight",3000,0.,2000.);
   for (i=1;i<3000;i++) {
      noise = gRandom->Gaus(0,120);
      if (i>700) {
         noise += 1000*sin((i-700)*6.28/30)*exp((double)(700-i)/300);
      }
      ht->SetBinContent(i,noise);
   }
   C->cd(1);
   ht->SetLineColor(2);
   ht->GetXaxis()->SetLabelSize(0.05);
   ht->Draw();
   ht->GetXaxis()->SetTimeDisplay(1);
   float x[100], t[100];
   for (i=0;i<100;i++) {
      x[i] = sin(i*4*3.1415926/50)*exp(-(double)i/20);
      t[i] = 6000+(double)i/20;
   }
   TGraph *gt = new TGraph(100,t,x);
   gt->SetTitle("Politics");
   C->cd(2);
   gt->SetFillColor(19);
   gt->SetLineColor(5);
   gt->SetLineWidth(2);
   gt->Draw("AL");
   gt->GetXaxis()->SetLabelSize(0.05);
   gt->GetXaxis()->SetTimeDisplay(1);
   gPad->Modified();
   float x2[10], t2[10];
   for (i=0;i<10;i++) {
      x2[i] = gRandom->Gaus(500,100)*i;
      t2[i] = i*365*86400;
   }
   TGraph *gt2 = new TGraph(10,t2,x2);
   gt2->SetTitle("Number of monkeys on the moon");
   C->cd(3);
   gt2->SetFillColor(19);
   gt2->SetMarkerColor(4);
   gt2->SetMarkerStyle(29);
   gt2->SetMarkerSize(1.3);
   gt2->Draw("AP");
   gt2->GetXaxis()->SetLabelSize(0.05);
   gt2->GetXaxis()->SetTimeDisplay(1);
   gt2->GetXaxis()->SetTimeFormat("y. %Y");

   TestReport1(C, "TGaxis 3 (Time on axis)");
   DoCcode(C);
   TestReport2();
}


////////////////////////////////////////////////////////////////////////////////
/// 4th TGaxis test.

void tgaxis4()
{
   TCanvas *C = StartTest(600,700);

   TDatime T0(2003,1,1,0,0,0);
   int X0 = T0.Convert();
   gStyle->SetTimeOffset(X0);
   TDatime T1(2002,9,23,0,0,0);
   int X1 = T1.Convert()-X0;
   TDatime T2(2003,3,7,0,0,0);
   int X2 = T2.Convert()-X0;
   TH1F * h1 = new TH1F("h1","test",100,X1,X2);
   TRandom r;
   for (Int_t i=0;i<30000;i++) {
      Double_t noise = r.Gaus(0.5*(X1+X2),0.1*(X2-X1));
      h1->Fill(noise);
   }
   h1->GetXaxis()->SetTimeDisplay(1);
   h1->GetXaxis()->SetLabelSize(0.03);
   h1->GetXaxis()->SetTimeFormat("%Y:%m:%d");
   h1->Draw();

   TestReport1(C, "TGaxis 4 (Time on axis)");
   DoCcode(C);
   TestReport2();
   delete h1;
}


////////////////////////////////////////////////////////////////////////////////
/// 5th TGaxis test.

void tgaxis5()
{
   TCanvas *C = StartTest(800,570);

   double f = 1.8;

   TLatex* tex1 = new TLatex;
   tex1->SetNDC();
   tex1->SetTextFont(102);
   tex1->SetTextSize(0.07*f);

   TLatex* tex3 = new TLatex;
   tex3->SetNDC();
   tex3->SetTextFont(102);
   tex3->SetTextSize(0.07*f);
   tex3->SetTextColor(kBlue+2);

   TLatex* tex2 = new TLatex;
   tex2->SetNDC();
   tex2->SetTextFont(102);
   tex2->SetTextSize(0.07*f);
   tex2->SetTextColor(kOrange+3);

   time_t offset[] = {0,                   0, 1325376000, 1341100800};
   time_t t[]      = {1331150400, 1336417200,          0, 36000};

   C->SetTopMargin(0);  C->SetBottomMargin(0);
   C->SetLeftMargin(0); C->SetRightMargin(0);
   C->Divide(2, 4, -1, -1);
   TLine l;
   l.DrawLine(0.5, 0, 0.5, 1.);

   for(int i = 0; i < 4; ++i){
      for(int gmt = 0; gmt < 2; ++gmt){
         const char* opt = (gmt ? "gmt" : "local");
         TVirtualPad* p = C->cd(2*i + gmt + 1);
         p->SetTopMargin(0); p->SetBottomMargin(0);
         p->SetLeftMargin(0); p->SetRightMargin(0);
         p->SetFillStyle(4000);

         TGaxis* ga = new TGaxis (.4, .25, 5., .25, t[i], t[i] + 1,  1, "t");
         ga->SetTimeFormat("TGaxis label: #color[2]{%Y-%m-%d %H:%M:%S}");
         ga->SetLabelFont(102);
         ga->SetLabelColor(kBlue+2);

         ga->SetTimeOffset(offset[i], opt);
         ga->SetLabelOffset(0.04*f);
         ga->SetLabelSize(0.07*f);
         ga->SetLineColor(0);
         ga->Draw();

         // Get offset string of axis time format: there is not accessor
         // to time format in TGaxis.
         // Assumes TAxis use the same format.
         TAxis a(10, 0, 1600000000);
         a.SetTimeOffset(offset[i], opt);
         const char* offsettimeformat = a.GetTimeFormat();

         char buf[256];
         if (offset[i] < t[i]) {
            sprintf(buf, "#splitline{%s, %s}{offset: %ld, option %s}",
                    stime(t+i).Data(), stime(t+i, true).Data(), (long) offset[i], opt);
         } else {
            int h = t[i] / 3600;
            int m = (t[i] - 3600 * h) / 60 ;
            int s = (t[i] - h * 3600 - m * 60);
            sprintf(buf, "#splitline{%d h %d m %d s}{offset: %s, option %s}",
                    h, m, s, stime(offset + i, gmt).Data(), opt);
         }
         tex1->DrawLatex(.01, .75, buf);
         tex2->DrawLatex(.01, .50, offsettimeformat);
         time_t t_ = t[i] + offset[i];
         sprintf(buf, "Expecting:    #color[2]{%s}", stime(&t_, gmt, false).Data());
         tex3->DrawLatex(.01, .24, buf);
         if(i > 0) l.DrawLine(0, 0.95, 1, 0.95);
      }
   }

   TestReport1(C, "TGaxis 5 (Time on axis: reference test)");
   DoCcode(C);
   TestReport2();
}


////////////////////////////////////////////////////////////////////////////////
/// Function used by tgaxis5

TString stime(time_t* t, bool utc, bool display_time_zone)
{
   struct tm* tt;
   if (utc) tt = gmtime(t);
   else     tt = localtime(t);
   char buf[256];
   if (display_time_zone) strftime(buf, sizeof(buf), "%H:%M:%S %Z", tt);
   else                   strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", tt);
   return TString(buf);
}


////////////////////////////////////////////////////////////////////////////////
/// Alphanumeric labels in a 1-d histogram

void labels1()
{
   TCanvas *C = StartTest(900,500);

   const Int_t nx = 20;

   C->SetGrid();
   C->SetBottomMargin(0.15);
   TH1F *hlab1 = new TH1F("hlab1","hlab1",nx,0,nx);
   hlab1->SetFillColor(38);
   for (Int_t i=0;i<5000;i++) {
      hlab1->Fill(gRandom->Gaus(0.5*nx,0.2*nx));
   }
   hlab1->SetStats(0);
   TAxis *xa = hlab1->GetXaxis();
   xa->SetBinLabel( 1, "Jean");
   xa->SetBinLabel( 2, "Pierre");
   xa->SetBinLabel( 3, "Marie");
   xa->SetBinLabel( 4, "Odile");
   xa->SetBinLabel( 5, "Sebastien");
   xa->SetBinLabel( 6, "Fons");
   xa->SetBinLabel( 7, "Rene");
   xa->SetBinLabel( 8, "Nicolas");
   xa->SetBinLabel( 9, "Xavier");
   xa->SetBinLabel(10, "Greg");
   xa->SetBinLabel(11, "Bjarne");
   xa->SetBinLabel(12, "Anton");
   xa->SetBinLabel(13, "Otto");
   xa->SetBinLabel(14, "Eddy");
   xa->SetBinLabel(15, "Peter");
   xa->SetBinLabel(16, "Pasha");
   xa->SetBinLabel(17, "Philippe");
   xa->SetBinLabel(18, "Suzanne");
   xa->SetBinLabel(19, "Jeff");
   xa->SetBinLabel(20, "Valery");
   hlab1->Draw();
   TPaveText *pt = new TPaveText(0.6,0.7,0.98,0.98,"brNDC");
   pt->SetFillColor(18);
   pt->SetTextAlign(12);
   pt->AddText("Use the axis Context Menu LabelsOption");
   pt->AddText(" \"a\"   to sort by alphabetic order");
   pt->AddText(" \">\"   to sort by decreasing values");
   pt->AddText(" \"<\"   to sort by increasing values");
   pt->Draw();

   TestReport1(C, "Alphanumeric labels in a 1-d histogram");
   DoCcode(C);
   TestReport2();
   delete hlab1;
}


////////////////////////////////////////////////////////////////////////////////
/// TEllipse test.

void tellipse()
{
   TCanvas *C = StartTest(700,800);

   C->Range(0,0,1,1);
   TPaveLabel pel(0.1,0.8,0.9,0.95,"Examples of Ellipses");
   pel.SetFillColor(42);
   pel.Draw();
   TEllipse el1(0.25,0.25,.1,.2);
   el1.Draw();
   el1.SetFillStyle(0);
   TEllipse el2(0.25,0.6,.2,.1);
   el2.SetFillColor(6);
   el2.SetFillStyle(3008);
   el2.Draw();
   TEllipse el3(0.75,0.6,.2,.1,45,315);
   el3.SetFillColor(2);
   el3.SetFillStyle(1001);
   el3.SetLineColor(4);
   el3.Draw();
   TEllipse el4(0.75,0.25,.2,.15,45,315,62);
   el4.SetFillColor(5);
   el4.SetFillStyle(1001);
   el4.SetLineColor(4);
   el4.SetLineWidth(6);
   el4.Draw();

   TestReport1(C, "TEllipse");
   DoCcode(C);
   TestReport2();
}


////////////////////////////////////////////////////////////////////////////////
/// Feynman diagrams test.

void feynman()
{
   TCanvas *C = StartTest(600,300);

   C->Range(0, 0, 140, 60);
   Int_t linsav = gStyle->GetLineWidth();
   gStyle->SetLineWidth(3);
   TLatex t;
   t.SetTextAlign(22);
   t.SetTextSize(0.1);
   TLine * l;
   l = new TLine(10, 10, 30, 30); l->Draw();
   l = new TLine(10, 50, 30, 30); l->Draw();
   TCurlyArc *ginit = new TCurlyArc(30, 30, 12.5*TMath::Sqrt(2), 135, 225);
   ginit->SetWavy();
   ginit->Draw();
   t.DrawLatex(7,6,"e^{-}");
   t.DrawLatex(7,55,"e^{+}");
   t.DrawLatex(7,30,"#gamma");
   TCurlyLine *gamma = new TCurlyLine(30, 30, 55, 30);
   gamma->SetWavy();
   gamma->Draw();
   t.DrawLatex(42.5,37.7,"#gamma");
   TArc *a = new TArc(70, 30, 15);
   a->Draw();
   a->SetFillStyle(0);
   t.DrawLatex(55, 45,"#bar{q}");
   t.DrawLatex(85, 15,"q");
   TCurlyLine *gluon = new TCurlyLine(70, 45, 70, 15);
   gluon->Draw();
   t.DrawLatex(77.5,30,"g");
   TCurlyLine *z0 = new TCurlyLine(85, 30, 110, 30);
   z0->SetWavy();
   z0->Draw();
   t.DrawLatex(100, 37.5,"Z^{0}");
   l = new TLine(110, 30, 130, 10); l->Draw();
   l = new TLine(110, 30, 130, 50); l->Draw();
   TCurlyArc *gluon1 = new TCurlyArc(110, 30, 12.5*TMath::Sqrt(2), 315, 45);
   gluon1->Draw();
   t.DrawLatex(135,6,"#bar{q}");
   t.DrawLatex(135,55,"q");
   t.DrawLatex(135,30,"g");
   C->Update();
   gStyle->SetLineWidth(linsav);

   TestReport1(C, "Feynman diagrams");
   DoCcode(C);
   TestReport2();
}


////////////////////////////////////////////////////////////////////////////////
/// Ratio plot test.

void ratioplot()
{
   TCanvas *C = StartTest(600,300);

   gStyle->SetOptStat(0);
   TH1D* rh1 = new TH1D("rh1", "rh1", 50, 0, 10);
   TH1D* rh2 = new TH1D("rh2", "rh2", 50, 0, 10);
   TF1*  rf1 = new TF1("rf1", "exp(- x/[0] )");
   rf1->SetParameter(0, 3);
   rh1->FillRandom("rf1", 1900);
   rh2->FillRandom("rf1", 2000);
   rh1->Sumw2();
   rh2->Scale(1.9 / 2.,"nosw2");
   TRatioPlot *rp = new TRatioPlot(rh1, rh2);
   C->SetTicks(0, 1);
   rp->Draw();

   TestReport1(C, "Ratio plot");
   DoCcode(C);
   TestReport2();
}


////////////////////////////////////////////////////////////////////////////////
/// 1st TGraph test.

void tgraph1()
{
   TCanvas *C = StartTest(700,500);

   C->SetFillColor(42);
   C->SetGrid();
   const Int_t n = 20;
   Double_t x[n], y[n];
   for (Int_t i=0;i<n;i++) {
     x[i] = i*0.1;
     y[i] = 10*sin(x[i]+0.2);
   }
   TGraph *gr = new TGraph(n,x,y);
   gr->SetLineColor(2);
   gr->SetLineWidth(4);
   gr->SetMarkerColor(4);
   gr->SetMarkerStyle(21);
   gr->SetTitle("a simple graph");
   gr->GetXaxis()->SetTitle("X title");
   gr->GetYaxis()->SetTitle("Y title");
   gr->Draw("ACP");
   C->Update();
   C->GetFrame()->SetFillColor(21);
   C->GetFrame()->SetBorderSize(12);

   TestReport1(C, "TGraph 1");
   DoCcode(C);
   TestReport2();
}


////////////////////////////////////////////////////////////////////////////////
/// 2nd TGraph test.

void tgraph2()
{
   TCanvas *C = StartTest(700,500);

   C->SetGrid();
   TMultiGraph *mg = new TMultiGraph();
   mg->SetTitle("Exclusion graphs");
   const Int_t n = 35;
   Double_t x1[n], x2[n], x3[n], y1[n], y2[n], y3[n];
   for (Int_t i=0;i<n;i++) {
     x1[i]  = i*0.1;
     x2[i]  = x1[i];
     x3[i]  = x1[i]+.5;
     y1[i] = 10*sin(x1[i]);
     y2[i] = 10*cos(x1[i]);
     y3[i] = 10*sin(x1[i])-2;
   }
   TGraph *gr1 = new TGraph(n,x1,y1);
   gr1->SetLineColor(2);
   gr1->SetLineWidth(1504);
   gr1->SetFillStyle(3005);
   TGraph *gr2 = new TGraph(n,x2,y2);
   gr2->SetLineColor(4);
   gr2->SetLineWidth(-2002);
   gr2->SetFillStyle(3004);
   gr2->SetFillColor(9);
   TGraph *gr3 = new TGraph(n,x3,y3);
   gr3->SetLineColor(5);
   gr3->SetLineWidth(-802);
   gr3->SetFillStyle(3002);
   gr3->SetFillColor(2);
   mg->Add(gr1);
   mg->Add(gr2);
   mg->Add(gr3);
   mg->Draw("AC");

   TestReport1(C, "TGraph 2 (Exclusion Zone)");
   DoCcode(C);
   TestReport2();
}


////////////////////////////////////////////////////////////////////////////////
/// 3rd TGraph test.

void tgraph3()
{
   TCanvas *C = StartTest(800,400);

   C->Divide(2,1);

   TGraph *g1 = new TGraph();
   g1->SetPoint(0, 1e-4, 1);
   g1->SetPoint(1, 1e-2, 2);
   g1->SetPoint(2, 1e-1, 3);
   g1->SetPoint(3, 1, 4);
   g1->SetPoint(4, 1e1, 5);
   g1->SetPoint(5, 1e2, 5);
   g1->SetPoint(6, 1e3, 4);
   g1->SetPoint(7, 1e4, 3);
   g1->SetPoint(8, 1e5, 2);
   g1->SetPoint(9, 1e6, 1);
   g1->SetTitle("10 blue circles should be visible");

   g1->SetMarkerStyle(kFullCircle);
   g1->SetMarkerSize(1.0);
   g1->SetMarkerColor(kBlue);
   g1->SetLineColor(kBlue);

   C->cd(1);
   g1->Fit("gaus","Q");
   g1->Draw("AP");
   gPad->SetLogx();

   C->cd(2);
   gPad->SetLogx();
   gPad->SetLogy();
   TGraph* g2 = new TGraph();
   for (int i = 0; i < 10; i++) g2->SetPoint(i, i + 1, i + 1);
   g2->SetTitle("2 log scales from 1e-2 to 1e2;x;y");
   g2->GetXaxis()->SetLimits(1e-2, 1e2);
   g2->GetHistogram()->SetMinimum(1e-2);
   g2->GetHistogram()->SetMaximum(1e2);
   g2->GetXaxis()->CenterTitle();
   g2->GetYaxis()->CenterTitle();
   g2->Draw("a*");

   TestReport1(C, "TGraph 3 (Fitting and log scales)");
   DoCcode(C);
   TestReport2();
}


////////////////////////////////////////////////////////////////////////////////
/// 4th TGraph test.

void tgraph4()
{
   TCanvas *C = StartTest(800,400);

   C->DivideSquare(2);

   const int numPoints = 10;
   double xValues[numPoints], yValues[numPoints];
   for (int i=0;i<numPoints;i++) {
      xValues[i] = i+2;
      yValues[i] = pow(10,i-3);
   }
   TGraph *g1 = new TGraph(numPoints, xValues, yValues);
   g1->SetTitle("These two plots should be the same");
   TGraph *g2 = new TGraph(numPoints, xValues, yValues);
   g2->SetTitle("");

   // Log x first
   C->cd(1);
   g1->Draw();
   gPad->SetLogx();
   C->Update();
   gPad->SetLogy();
   C->Update();

   // Log y first
   C->cd(2);
   g2->Draw();
   gPad->SetLogy();
   C->Update();
   gPad->SetLogx();

   TestReport1(C, "TGraph 4 (Log scales setting order)");
   DoCcode(C);
   TestReport2();
}


////////////////////////////////////////////////////////////////////////////////
/// TH2Poly test.

void th2poly()
{
   TCanvas *C = StartTest(800,400);

   Int_t i;
   const Int_t nx = 48;
   const char *states [nx] = {
      "alabama",      "arizona",        "arkansas",       "california",
      "colorado",     "connecticut",    "delaware",       "florida",
      "georgia",      "idaho",          "illinois",       "indiana",
      "iowa",         "kansas",         "kentucky",       "louisiana",
      "maine",        "maryland",       "massachusetts",  "michigan",
      "minnesota",    "mississippi",    "missouri",       "montana",
      "nebraska",     "nevada",         "new_hampshire",  "new_jersey",
      "new_mexico",   "new_york",       "north_carolina", "north_dakota",
      "ohio",         "oklahoma",       "oregon",         "pennsylvania",
      "rhode_island", "south_carolina", "south_dakota",   "tennessee",
      "texas",        "utah",           "vermont",        "virginia",
      "washington",   "west_virginia",  "wisconsin",      "wyoming"
   };
   Double_t pop[nx] = {
    4708708, 6595778,  2889450, 36961664, 5024748,  3518288,  885122, 18537969,
    9829211, 1545801, 12910409,  6423113, 3007856,  2818747, 4314113,  4492076,
    1318301, 5699478,  6593587,  9969727, 5266214,  2951996, 5987580,   974989,
    1796619, 2643085,  1324575,  8707739, 2009671, 19541453, 9380884,   646844,
   11542645, 3687050,  3825657, 12604767, 1053209,  4561242,  812383,  6296254,
   24782302, 2784572,   621760,  7882590, 6664195,  1819777, 5654774,   544270
   };

   Double_t lon1 = -130;
   Double_t lon2 = -65;
   Double_t lat1 = 24;
   Double_t lat2 = 50;
   TH2Poly *p = new TH2Poly("USA","USA Population",lon1,lon2,lat1,lat2);
   gErrorIgnoreLevel = 9999;
   TFile::SetCacheFileDir(".");
   TFile *f = TFile::Open("http://root.cern.ch/files/usa.root", "CACHEREAD");

   if (!f) {
      printf("Cannot access usa.root. Is internet working ?\n");
      return;
   }

   // Define the TH2Poly bins.
   TMultiGraph *mg;
   TKey *key;
   TIter nextkey(gDirectory->GetListOfKeys());
   while ((key = (TKey*)nextkey())) {
      TObject *obj = key->ReadObj();
      if (obj->InheritsFrom("TMultiGraph")) {
         mg = (TMultiGraph*)obj;
         p->AddBin(mg);
      }
   }

   // Fill TH2Poly.
   for (i=0; i<nx; i++) p->Fill(states[i], pop[i]);

   gStyle->SetOptStat(11);
   gStyle->SetPalette(kBird);
   p->DrawClone("COL");

   TestReport1(C, "TH2Poly.(DrawClone() and remote file access)");
   DoCcode(C);
   TestReport2();
}


////////////////////////////////////////////////////////////////////////////////
/// TMultigraph and TGraphErrors test

void tmultigraph1()
{
   TCanvas *C = StartTest(700,500);

   gStyle->SetOptFit();
   C->SetGrid();
   TMultiGraph *mg = new TMultiGraph();
   Int_t n1 = 10;
   Double_t x1[]  = {-0.1, 0.05, 0.25, 0.35, 0.5, 0.61,0.7,0.85,0.89,0.95};
   Double_t y1[]  = {-1,2.9,5.6,7.4,9,9.6,8.7,6.3,4.5,1};
   Double_t ex1[] = {.05,.1,.07,.07,.04,.05,.06,.07,.08,.05};
   Double_t ey1[] = {.8,.7,.6,.5,.4,.4,.5,.6,.7,.8};
   TGraphErrors *gr1 = new TGraphErrors(n1,x1,y1,ex1,ey1);
   gr1->SetMarkerColor(kBlue);
   gr1->SetMarkerStyle(21);
   gr1->Fit("pol6","q ex0");
   mg->Add(gr1);
   Int_t n2 = 10;
   Float_t x2[]  = {-0.28, 0.005, 0.19, 0.29, 0.45, 0.56,0.65,0.80,0.90,1.01};
   Float_t y2[]  = {2.1,3.86,7,9,10,10.55,9.64,7.26,5.42,2};
   Float_t ex2[] = {.04,.12,.08,.06,.05,.04,.07,.06,.08,.04};
   Float_t ey2[] = {.6,.8,.7,.4,.3,.3,.4,.5,.6,.7};
   TGraphErrors *gr2 = new TGraphErrors(n2,x2,y2,ex2,ey2);
   gr2->SetMarkerColor(kRed);
   gr2->SetMarkerStyle(20);
   gr2->Fit("pol5","q ex0");
   mg->Add(gr2);
   mg->Draw("ap");
   C->Update();
   TPaveStats *stats1 = (TPaveStats*)gr1->GetListOfFunctions()->FindObject("stats");
   TPaveStats *stats2 = (TPaveStats*)gr2->GetListOfFunctions()->FindObject("stats");
   stats1->SetTextColor(kBlue);
   stats2->SetTextColor(kRed);
   stats1->SetX1NDC(0.12); stats1->SetX2NDC(0.32); stats1->SetY1NDC(0.75);
   stats2->SetX1NDC(0.72); stats2->SetX2NDC(0.92); stats2->SetY1NDC(0.78);
   C->Modified();

   TestReport1(C, "TMultigraph and TGraphErrors");
   DoCcode(C);
   TestReport2();
}


////////////////////////////////////////////////////////////////////////////////
/// All Kind of TMultigraph test

void tmultigraph2()
{
   TCanvas *C = StartTest(800,800);

   gStyle->SetOptFit();

   C->SetGrid();
   C->Divide(2,2);

   // Draw a frame to define the range
   TMultiGraph *mg1 = new TMultiGraph();
   TMultiGraph *mg2 = new TMultiGraph();
   TMultiGraph *mg3 = new TMultiGraph();
   TMultiGraph *mg4 = new TMultiGraph();

   // Vectors used to build the graphs
   Int_t n1 = 10;
   Double_t x1[]    = {-0.1, 0.05, 0.25, 0.35, 0.5, 0.61,0.7,0.85,0.89,0.95};
   Double_t y1[]    = {-1,2.9,5.6,7.4,9,9.6,8.7,6.3,4.5,1};
   Double_t exl1[]  = {.05,.1,.07,.07,.04,.05,.06,.07,.08,.05};
   Double_t eyl1[]  = {.8,.7,.6,.5,.4,.4,.5,.6,.7,.8};
   Double_t exh1[]  = {.02,.08,.05,.05,.03,.03,.04,.05,.06,.03};
   Double_t eyh1[]  = {.6,.5,.4,.3,.2,.2,.3,.4,.5,.6};
   Double_t exld1[] = {.0,.0,.0,.0,.0,.0,.0,.0,.0,.0};
   Double_t eyld1[] = {.0,.0,.05,.0,.0,.0,.0,.0,.0,.0};
   Double_t exhd1[] = {.0,.0,.0,.0,.0,.0,.0,.0,.0,.0};
   Double_t eyhd1[] = {.0,.0,.0,.0,.0,.0,.0,.0,.05,.0};

   Int_t n2 = 10;
   Float_t  x2[]    = {-0.28, 0.005, 0.19, 0.29, 0.45, 0.56,0.65,0.80,0.90,1.01};
   Float_t  y2[]    = {2.1,3.86,7,9,10,10.55,9.64,7.26,5.42,2};
   Float_t  exl2[]  = {.04,.12,.08,.06,.05,.04,.07,.06,.08,.04};
   Float_t  eyl2[]  = {.6,.8,.7,.4,.3,.3,.4,.5,.6,.7};
   Float_t  exh2[]  = {.02,.08,.05,.05,.03,.03,.04,.05,.06,.03};
   Float_t  eyh2[]  = {.6,.5,.4,.3,.2,.2,.3,.4,.5,.6};
   Float_t  exld2[] = {.0,.0,.0,.0,.0,.0,.0,.0,.0,.0};
   Float_t  eyld2[] = {.0,.0,.05,.0,.0,.0,.0,.0,.0,.0};
   Float_t  exhd2[] = {.0,.0,.0,.0,.0,.0,.0,.0,.0,.0};
   Float_t  eyhd2[] = {.0,.0,.0,.0,.0,.0,.0,.0,.05,.0};

   // Create 1st multigraph
   C->cd(1);
   TGraph *gr11 = new TGraph(n1,x1,y1);
   gr11->SetMarkerColor(kBlue);
   gr11->SetMarkerStyle(20);
   TGraph *gr12 = new TGraph(n2,x2,y2);
   gr12->SetMarkerColor(kRed);
   gr12->SetMarkerStyle(23);
   mg1->Add(gr11,"pc");
   mg1->Add(gr12);
   mg1->Draw("ap");
   TLatex *tex1 = new TLatex(-0.3,10.0,"TGraph");
   tex1->Draw();

   // Create 2nd multigraph
   C->cd(2);
   TGraphErrors *gr21 = new TGraphErrors(n1,x1,y1,exl1,eyl1);
   gr21->SetMarkerColor(kBlue+1);
   gr21->SetMarkerStyle(21);
   TGraphErrors *gr22 = new TGraphErrors(n2,x2,y2,exl2,eyl2);
   gr22->SetMarkerColor(kRed+1);
   gr22->SetMarkerStyle(20);
   gr22->Fit("pol3","q ex0");
   mg2->Add(gr21,"pl");
   mg2->Add(gr22);
   mg2->Draw("ap");
   TLatex *tex2 = new TLatex(-0.3,10.0,"TGraphErrors");
   tex2->Draw();

   // Create 3rd multigraph
   C->cd(3);
   TGraphAsymmErrors *gr31 = new TGraphAsymmErrors(n1,x1,y1,exl1,exh1,eyl1,eyh1);
   gr31->SetMarkerColor(kBlue-1);
   gr31->SetMarkerStyle(21);
   TGraphAsymmErrors *gr32 = new TGraphAsymmErrors(n2,x2,y2,exl2,exh2,eyl2,eyh2);
   gr32->SetMarkerColor(kRed-1);
   gr32->SetMarkerStyle(20);
   gr32->Fit("pol4","q ex0");
   mg3->Add(gr31,"pl");
   mg3->Add(gr32);
   mg3->Draw("ap");
   TLatex *tex3 = new TLatex(-0.3,10.0,"TGraphAsymmErrors");
   tex3->Draw();

   // Create 4th multigraph
   C->cd(4);
   TGraphBentErrors *gr41 = new TGraphBentErrors(n1,x1,y1,exl1,exh1,eyl1,eyh1,exld1,exhd1,eyld1,eyhd1);
   gr41->SetMarkerColor(kGreen);
   gr41->SetMarkerStyle(21);
   TGraphBentErrors *gr42 = new TGraphBentErrors(n2,x2,y2,exl2,exh2,eyl2,eyh2,exld2,exhd2,eyld2,eyhd2);
   gr42->SetMarkerColor(kViolet);
   gr42->SetMarkerStyle(20);
   gr42->Fit("pol5","q ex0");
   mg4->Add(gr41,"pc");
   mg4->Add(gr42);
   mg4->Draw("ap");
   TLatex *tex4 = new TLatex(-0.3,10.0,"TGraphBentErrors");
   tex4->Draw();

   C->Modified();

   TestReport1(C, "All Kind of TMultigraph");
   DoCcode(C);
   TestReport2();
}


////////////////////////////////////////////////////////////////////////////////
/// 1st 2D options Test

void options2d1()
{
   TCanvas *C = StartTest(800,600);

   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetCanvasColor(33);
   gStyle->SetFrameFillColor(18);
   TF2 *f2 = new TF2("f2","xygaus + xygaus(5) + xylandau(10)",-4,4,-4,4);
   Double_t params[] = {130,-1.4,1.8,1.5,1, 150,2,0.5,-2,0.5, 3600,-2,0.7,-3,0.3};
   f2->SetParameters(params);
   gH2 = new TH2F("h2","xygaus + xygaus(5) + xylandau(10)",20,-4,4,20,-4,4);
   gH2->SetFillColor(46);
   gH2->FillRandom("f2",40000);

   TPaveLabel pl1;
   Float_t x1=0.67, y1=0.875, x2=0.85, y2=0.95;
   C->Divide(2,2);
   C->SetFillColor(17);
   C->cd(1);
   gH2->Draw();       pl1.DrawPaveLabel(x1,y1,x2,y2,"SCAT","brNDC");
   C->cd(2);
   gH2->Draw("box");  pl1.DrawPaveLabel(x1,y1,x2,y2,"BOX","brNDC");
   C->cd(3);
   gH2->Draw("arr");  pl1.DrawPaveLabel(x1,y1,x2,y2,"ARR","brNDC");
   C->cd(4);
   gH2->Draw("colz"); pl1.DrawPaveLabel(x1,y1,x2,y2,"COLZ","brNDC");

   TestReport1(C, "Basic 2D options");
   DoCcode(C);
   TestReport2();
}


////////////////////////////////////////////////////////////////////////////////
/// 2nd 2D options Test

void options2d2()
{
   TCanvas *C = StartTest(800,600);

   TPaveLabel pl2;
   Float_t x1=0.67, y1=0.875, x2=0.85, y2=0.95;
   gPad->SetGrid();
   C->SetFillColor(17);
   C->SetGrid();
   gH2->Draw("text"); pl2.DrawPaveLabel(x1,y1,x2,y2,"TEXT","brNDC");

   TestReport1(C, "Text option");
   DoCcode(C);
   TestReport2();
}


////////////////////////////////////////////////////////////////////////////////
/// 3rd 2D options Test

void options2d3()
{
   TCanvas *C = StartTest(800,600);

   TPaveLabel pl3;
   Float_t x1=0.67, y1=0.875, x2=0.85, y2=0.95;
   C->Divide(2,2);
   gPad->SetGrid();
   C->SetFillColor(17);
   gStyle->SetPalette(kBird);
   C->cd(1);
   gH2->Draw("cont4z"); pl3.DrawPaveLabel(x1,y1,x2,y2,"CONT4Z","brNDC");
   C->cd(2);
   gPad->SetGrid();
   gH2->Draw("cont1"); pl3.DrawPaveLabel(x1,y1,x2,y2,"CONT1","brNDC");
   C->cd(3);
   gPad->SetGrid();
   gH2->Draw("cont2"); pl3.DrawPaveLabel(x1,y1,x2,y2,"CONT2","brNDC");
   C->cd(4);
   gPad->SetGrid();
   gH2->Draw("cont3"); pl3.DrawPaveLabel(x1,y1,x2,y2,"CONT3","brNDC");

   TestReport1(C, "Contour options");
   DoCcode(C);
   TestReport2();
}


////////////////////////////////////////////////////////////////////////////////
/// 4th 2D options Test

void options2d4()
{
   TCanvas *C = StartTest(800,600);

   TPaveLabel pl4;
   Float_t x1=0.67, y1=0.875, x2=0.85, y2=0.95;
   C->Divide(2,2);
   C->SetFillColor(17);
   C->cd(1);
   gStyle->SetPalette(kBird);
   gH2->Draw("lego");     pl4.DrawPaveLabel(x1,y1,x2,y2,"LEGO","brNDC");
   C->cd(2);
   gH2->Draw("lego1");    pl4.DrawPaveLabel(x1,y1,x2,y2,"LEGO1","brNDC");
   C->cd(3);
   gPad->SetTheta(61); gPad->SetPhi(-82);
   gH2->Draw("surf1pol"); pl4.DrawPaveLabel(x1,y1,x2+0.05,y2,"SURF1POL","brNDC");
   C->cd(4);
   gPad->SetTheta(21); gPad->SetPhi(-90);
   gH2->Draw("surf1cyl"); pl4.DrawPaveLabel(x1,y1,x2+0.05,y2,"SURF1CYL","brNDC");

   TestReport1(C, "Lego options");
   DoCcode(C);
   TestReport2();
}


////////////////////////////////////////////////////////////////////////////////
/// 5th 2D options Test

void options2d5()
{
   TCanvas *C = StartTest(800,600);

   TPaveLabel pl5;
   Float_t x1=0.67, y1=0.875, x2=0.85, y2=0.95;
   C->Divide(2,2);
   C->SetFillColor(17);
   gStyle->SetPalette(kBird);
   C->cd(1);
   gH2->Draw("surf1");   pl5.DrawPaveLabel(x1,y1,x2,y2,"SURF1","brNDC");
   C->cd(2);
   gH2->Draw("surf2z");  pl5.DrawPaveLabel(x1,y1,x2,y2,"SURF2Z","brNDC");
   C->cd(3);
   gH2->Draw("surf3");   pl5.DrawPaveLabel(x1,y1,x2,y2,"SURF3","brNDC");
   C->cd(4);
   gH2->Draw("surf4");   pl5.DrawPaveLabel(x1,y1,x2,y2,"SURF4","brNDC");

   TestReport1(C, "Surface options");
   DoCcode(C);
   TestReport2();
   delete gH2;
}


////////////////////////////////////////////////////////////////////////////////
/// 6th 2D options Test

void earth()
{
   TCanvas *C = StartTest(1000,800);

   gStyle->SetPalette(kBird);
   gStyle->SetOptTitle(1);
   gStyle->SetOptStat(0);
   C->Divide(2,2);
   TH2F *h1 = new TH2F("h01","Aitoff",    50, -180, 180, 50, -89.5, 89.5);
   TH2F *h2 = new TH2F("h02","Mercator",  50, -180, 180, 50, -80.5, 80.5);
   TH2F *h3 = new TH2F("h03","Sinusoidal",50, -180, 180, 50, -90.5, 90.5);
   TH2F *h4 = new TH2F("h04","Parabolic", 50, -180, 180, 50, -90.5, 90.5);
   std::ifstream in;
   in.open("../tutorials/graphics/earth.dat");
   if (!in) {
      in.clear();
      in.open("earth.dat");
   }
   if (!in)
      printf("Cannot find earth.dat!\n");
   Float_t x,y;
   while (1) {
     in >> x >> y;
     if (!in.good()) break;
     h1->Fill(x,y, 1);
     h2->Fill(x,y, 1);
     h3->Fill(x,y, 1);
     h4->Fill(x,y, 1);
   }
   in.close();
   C->cd(1); h1->Draw("z aitoff");
   C->cd(2); h2->Draw("z mercator");
   C->cd(3); h3->Draw("z sinusoidal");
   C->cd(4); h4->Draw("z parabolic");

   TestReport1(C, "Special contour options (AITOFF etc.)");
   DoCcode(C);
   TestReport2();
   delete h1;
   delete h2;
   delete h3;
   delete h4;
}


////////////////////////////////////////////////////////////////////////////////
/// 1st TGraph2D Test

void tgraph2d1()
{
   TCanvas *C = StartTest(600,600);

   Double_t P = 5.;
   Int_t npx  = 20 ;
   Int_t npy  = 20 ;
   Double_t x = -P;
   Double_t y = -P;
   Double_t z;
   Int_t k = 0;
   Double_t dx = (2*P)/npx;
   Double_t dy = (2*P)/npy;
   TGraph2D *dt = new TGraph2D(npx*npy);
   dt->SetName("Graph2DA");
   dt->SetNpy(41);
   dt->SetNpx(40);
   for (Int_t i=0; i<npx; i++) {
      for (Int_t j=0; j<npy; j++) {
         z = sin(sqrt(x*x+y*y))+1;
         dt->SetPoint(k,x,y,z);
         k++;
         y = y+dy;
      }
       x = x+dx;
       y = -P;
   }
   gStyle->SetPalette(kBird);
   dt->SetFillColor(0);
   dt->SetLineColor(1);
   dt->SetMarkerSize(1);
   dt->Draw("tri2p0Z  ");

   TestReport1(C, "TGraph2D 1 (TRI2 and P0)");

   DoCcode(C);

   TObject *old = (TObject*)gDirectory->GetList()->FindObject(dt->GetName());
   if (old) gDirectory->GetList()->Remove(old);

   TestReport2();
   delete dt;
}


////////////////////////////////////////////////////////////////////////////////
/// 2nd TGraph2D Test

void tgraph2d2()
{
   TCanvas *C = StartTest(600,600);

   gStyle->SetPadBorderMode(0);
   gStyle->SetFrameBorderMode(0);
   gStyle->SetCanvasBorderMode(0);
   Double_t Px = 6.;
   Double_t Py = 6.;
   Int_t np    = 1000;
   Double_t *rx=0, *ry=0, *rz=0;
   rx = new Double_t[np];
   ry = new Double_t[np];
   rz = new Double_t[np];
   TRandom *r = new TRandom();
   for (Int_t N=0; N<np; N++) {
      rx[N]=2*Px*(r->Rndm())-Px;
      ry[N]=2*Py*(r->Rndm())-Py;
      rz[N]=sin(sqrt(rx[N]*rx[N]+ry[N]*ry[N]))+1;
   }
   gStyle->SetPalette(kBird);
   TGraph2D *dt = new TGraph2D( np, rx, ry, rz);
   dt->SetName("Graph2DA");
   dt->SetFillColor(0);
   dt->SetMarkerStyle(20);
   dt->Draw("PCOL");

   TestReport1(C, "TGraph2D 2 (COL and P)");
   DoCcode(C);

   TObject *old = (TObject*)gDirectory->GetList()->FindObject(dt->GetName());
   if (old) gDirectory->GetList()->Remove(old);

   TestReport2();
   delete dt;
   delete [] rx;
   delete [] ry;
   delete [] rz;
}


////////////////////////////////////////////////////////////////////////////////
/// 3rd TGraph2D Test

void tgraph2d3()
{
   TCanvas *C = StartTest(600,600);

   gStyle->SetPadBorderMode(0);
   gStyle->SetFrameBorderMode(0);
   gStyle->SetCanvasBorderMode(0);
   Double_t Px = 6.;
   Double_t Py = 6.;
   Int_t np    = 200;
   Double_t *rx=0, *ry=0, *rz=0;
   rx = new Double_t[np];
   ry = new Double_t[np];
   rz = new Double_t[np];
   TRandom *r = new TRandom();
   for (Int_t N=0; N<np; N++) {
      rx[N]=2*Px*(r->Rndm())-Px;
      ry[N]=2*Py*(r->Rndm())-Py;
      rz[N]=sin(sqrt(rx[N]*rx[N]+ry[N]*ry[N]))+1;
   }
   gStyle->SetPalette(kBird);
   TGraph2D *dt = new TGraph2D( np, rx, ry, rz);
   dt->SetName("Graph2DA");
   dt->SetFillColor(0);
   dt->Draw("CONT5  ");

   TestReport1(C, "TGraph2D 3 (CONT5)");
   DoCcode(C);

   TObject *old = (TObject*)gDirectory->GetList()->FindObject(dt->GetName());
   if (old) gDirectory->GetList()->Remove(old);

   TestReport2();
   delete dt;
   delete [] rx;
   delete [] ry;
   delete [] rz;
}


////////////////////////////////////////////////////////////////////////////////
/// 1st complex drawing and TPad test

void ntuple1()
{
   TCanvas *C = StartTest(700,780);

   TPad *pad1 = new TPad("pad1","This is pad1",0.02,0.52,0.48,0.98,21);
   TPad *pad2 = new TPad("pad2","This is pad2",0.52,0.52,0.98,0.98,21);
   TPad *pad3 = new TPad("pad3","This is pad3",0.02,0.02,0.48,0.48,21);
   TPad *pad4 = new TPad("pad4","This is pad4",0.52,0.02,0.98,0.48,1);
   pad1->Draw();
   pad2->Draw();
   pad3->Draw();
   pad4->Draw();
   gStyle->SetStatW(0.30);
   gStyle->SetStatH(0.20);
   gStyle->SetStatColor(42);
   pad1->cd();
   pad1->SetGrid();
   pad1->SetLogy();
   pad1->GetFrame()->SetFillColor(15);
   TNtuple *ntuple = (TNtuple*)gHsimple->Get("ntuple");
   ntuple->SetLineColor(1);
   ntuple->SetFillStyle(1001);
   ntuple->SetFillColor(45);
   ntuple->Draw("3*px+2","px**2+py**2>1");
   ntuple->SetFillColor(38);
   ntuple->Draw("2*px+2","pz>2","same");
   ntuple->SetFillColor(5);
   ntuple->Draw("1.3*px+2","(px^2+py^2>4) && py>0","same");
   pad1->RedrawAxis();
   pad2->cd();
   pad2->SetGrid();
   pad2->GetFrame()->SetFillColor(32);
   ntuple->Draw("pz:px>>hprofs","","goffprofs");
   TProfile *hprofs = (TProfile*)gDirectory->Get("hprofs");
   hprofs->SetMarkerColor(5);
   hprofs->SetMarkerSize(0.7);
   hprofs->SetMarkerStyle(21);
   hprofs->Fit("pol2","q");
   TF1 *fpol2 = hprofs->GetFunction("pol2");
   fpol2->SetLineWidth(4);
   fpol2->SetLineColor(2);
   pad3->cd();
   pad3->GetFrame()->SetFillColor(38);
   pad3->GetFrame()->SetBorderSize(8);
   ntuple->SetMarkerColor(1);
   ntuple->Draw("py:px","pz>1");
   ntuple->SetMarkerColor(2);
   ntuple->Draw("py:px","pz<1","same");
   pad4->cd();
   ntuple->Draw("pz:py:px","(pz<10 && pz>6)+(pz<4 && pz>3)");
   ntuple->SetMarkerColor(4);
   ntuple->Draw("pz:py:px","pz<6 && pz>4","same");
   ntuple->SetMarkerColor(5);
   ntuple->Draw("pz:py:px","pz<4 && pz>3","same");
   TPaveText *l4 = new TPaveText(-0.9,0.5,0.9,0.95);
   l4->SetFillColor(42);
   l4->SetTextAlign(12);
   l4->AddText("You can interactively rotate this view in 2 ways:");
   l4->AddText("  - With the RotateCube in clicking in this pad");
   l4->AddText("  - Selecting View with x3d in the View menu");
   l4->Draw();
   gStyle->SetStatColor(19);

   TestReport1(C, "Ntuple drawing and TPad");
   DoCcode(C);
   TestReport2();
}


////////////////////////////////////////////////////////////////////////////////
/// 2nd complex drawing and TPad test

void quarks()
{
   TCanvas *C = StartTest(630,760);

   C->SetFillColor(kBlack);
   Int_t quarkColor  = 50;
   Int_t leptonColor = 16;
   Int_t forceColor  = 38;
   Int_t titleColor  = kYellow;
   Int_t border = 8;
   TLatex *texf = new TLatex(0.90,0.455,"Force Carriers");
   texf->SetTextColor(forceColor);
   texf->SetTextAlign(22); texf->SetTextSize(0.07); texf->SetTextAngle(90);
   texf->Draw();
   TLatex *texl = new TLatex(0.11,0.288,"Leptons");
   texl->SetTextColor(leptonColor);
   texl->SetTextAlign(22); texl->SetTextSize(0.07); texl->SetTextAngle(90);
   texl->Draw();
   TLatex *texq = new TLatex(0.11,0.624,"Quarks");
   texq->SetTextColor(quarkColor);
   texq->SetTextAlign(22); texq->SetTextSize(0.07); texq->SetTextAngle(90);
   texq->Draw();
   TLatex tex1(0.5,0.5,"u");
   tex1.SetTextColor(titleColor); tex1.SetTextFont(32); tex1.SetTextAlign(22);
   tex1.SetTextSize(0.14); tex1.DrawLatex(0.5,0.93,"Elementary");
   tex1.SetTextSize(0.12); tex1.DrawLatex(0.5,0.84,"Particles");
   tex1.SetTextSize(0.05); tex1.DrawLatex(0.5,0.067,"Three Generations of Matter");
   tex1.SetTextColor(kBlack); tex1.SetTextSize(0.8);
   TPad *pad = new TPad("pad", "pad",0.15,0.11,0.85,0.79);
   pad->Draw();
   pad->cd();
   pad->Divide(4,4,0.0003,0.0003);
   pad->cd(1); gPad->SetFillColor(quarkColor);   gPad->SetBorderSize(border);
   tex1.DrawLatex(.5,.5,"u");
   pad->cd(2); gPad->SetFillColor(quarkColor);   gPad->SetBorderSize(border);
   tex1.DrawLatex(.5,.5,"c");
   pad->cd(3); gPad->SetFillColor(quarkColor);   gPad->SetBorderSize(border);
   tex1.DrawLatex(.5,.5,"t");
   pad->cd(4); gPad->SetFillColor(forceColor);   gPad->SetBorderSize(border);
   tex1.DrawLatex(.5,.55,"#gamma");
   pad->cd(5); gPad->SetFillColor(quarkColor);   gPad->SetBorderSize(border);
   tex1.DrawLatex(.5,.5,"d");
   pad->cd(6); gPad->SetFillColor(quarkColor);   gPad->SetBorderSize(border);
   tex1.DrawLatex(.5,.5,"s");
   pad->cd(7); gPad->SetFillColor(quarkColor);   gPad->SetBorderSize(border);
   tex1.DrawLatex(.5,.5,"b");
   pad->cd(8); gPad->SetFillColor(forceColor);   gPad->SetBorderSize(border);
   tex1.DrawLatex(.5,.55,"g");
   pad->cd(9); gPad->SetFillColor(leptonColor);  gPad->SetBorderSize(border);
   tex1.DrawLatex(.5,.5,"#nu_{e}");
   pad->cd(10); gPad->SetFillColor(leptonColor); gPad->SetBorderSize(border);
   tex1.DrawLatex(.5,.5,"#nu_{#mu}");
   pad->cd(11); gPad->SetFillColor(leptonColor); gPad->SetBorderSize(border);
   tex1.DrawLatex(.5,.5,"#nu_{#tau}");
   pad->cd(12); gPad->SetFillColor(forceColor);  gPad->SetBorderSize(border);
   tex1.DrawLatex(.5,.5,"Z");
   pad->cd(13); gPad->SetFillColor(leptonColor); gPad->SetBorderSize(border);
   tex1.DrawLatex(.5,.5,"e");
   pad->cd(14); gPad->SetFillColor(leptonColor); gPad->SetBorderSize(border);
   tex1.DrawLatex(.5,.56,"#mu");
   pad->cd(15); gPad->SetFillColor(leptonColor); gPad->SetBorderSize(border);
   tex1.DrawLatex(.5,.5,"#tau");
   pad->cd(16); gPad->SetFillColor(forceColor);  gPad->SetBorderSize(border);
   tex1.DrawLatex(.5,.5,"W");
   C->cd();

   TestReport1(C, "Divided pads and TLatex");
   DoCcode(C);
   TestReport2();
}


////////////////////////////////////////////////////////////////////////////////
/// TImage test

void timage()
{
   TCanvas *C = StartTest(800,800);

   TImage *img = TImage::Open("$(ROOTSYS)/tutorials/image/rose512.jpg");
   if (!img) {
      printf("Could not create an image... exit\n");
      return;
   }
   TImage *i1 = TImage::Open("$(ROOTSYS)/tutorials/image/rose512.jpg");
   i1->SetConstRatio(kFALSE);
   i1->Flip(90);
   TImage *i2 = TImage::Open("$(ROOTSYS)/tutorials/image/rose512.jpg");
   i2->SetConstRatio(kFALSE);
   i2->Flip(180);
   TImage *i3 = TImage::Open("$(ROOTSYS)/tutorials/image/rose512.jpg");
   i3->SetConstRatio(kFALSE);
   i3->Flip(270);
   TImage *i4 = TImage::Open("$(ROOTSYS)/tutorials/image/rose512.jpg");
   i4->SetConstRatio(kFALSE);
   i4->Mirror(kTRUE);
   float d = 0.40;
   TPad *p1 = new TPad("i1", "i1", 0.05, 0.55, 0.05+d*i1->GetWidth()/i1->GetHeight(), 0.95);
   TPad *p2 = new TPad("i2", "i2", 0.55, 0.55, 0.95, 0.55+d*i2->GetHeight()/i2->GetWidth());
   TPad *p3 = new TPad("i3", "i3", 0.55, 0.05, 0.55+d*i3->GetWidth()/i3->GetHeight(), 0.45);
   TPad *p4 = new TPad("i4", "i4", 0.05, 0.05, 0.45, 0.05+d*i4->GetHeight()/i4->GetWidth());
   p1->Draw();
   p1->cd();
   i1->Draw();
   C->cd();
   p2->Draw();
   p2->cd();
   i2->Draw();
   C->cd();
   p3->Draw();
   p3->cd();
   i3->Draw();
   C->cd();
   p4->Draw();
   p4->cd();
   i4->Draw();
   C->cd();

   TestReport1(C, "TImage");
   DoCcode(C);
   TestReport2();
}


////////////////////////////////////////////////////////////////////////////////
/// Zoom/UnZoom a collection of TF1

double fg(double *x, double *p) {return sin((*p)*(*x));}
void zoomtf1()
{
   TCanvas *C = StartTest(800,800);

   TF1* f[6];

   for (int i=0;i<6;++i) {
      f[i] = new TF1(Form("f%d",i),fg, 0,2, 1);
      f[i]->SetParameter(0,i+1);
      f[i]->SetLineColor(i+1);
      if (i==0) {
         f[i]->GetYaxis()->SetTitle("Y axis title on the right side of the plot");
         f[i]->GetXaxis()->SetTitle("X axis");
         f[i]->GetYaxis()->SetTitleOffset(-13.);
         f[i]->SetTitle("Zoom/UnZoom a collection of TF1 and change axis attributes");
      }
      f[i]->Draw(i?"same":"");
   }
   f[0]->GetXaxis()->SetRangeUser(.1,.3);
   gPad->Update();
   f[0]->GetXaxis()->UnZoom();
   gPad->Modified();

   TestReport1(C, "Zoom/UnZoom a collection of TF1");
   if (gOptionR) printf("%10d%10d\n",0,0);
}


////////////////////////////////////////////////////////////////////////////////
/// Zoom/UnZoom a fitted histogram

void zoomfit()
{
   TCanvas *C = StartTest(800,800);

   TH1 *hpx = (TH1*)gHsimple->Get("hpx");
   hpx->Fit("gaus","q");
   hpx->GetXaxis()->SetRangeUser(.1,.3);
   gPad->Modified();
   gPad->Update();
   hpx->GetXaxis()->UnZoom();
   gPad->Modified();
   gPad->Update();

   TestReport1(C, "Zoom/UnZoom a fitted histogram");
   DoCcode(C);
   TestReport2();
}


////////////////////////////////////////////////////////////////////////////////
/// Ntuple drawing with alphanumeric variables

void hbars()
{
   TCanvas *C = StartTest(700,800);

   TTree *T = (TTree*)gCernstaff->Get("T");
   T->SetFillColor(45);
   C->SetFillColor(42);
   C->Divide(1,2);

   //horizontal bar chart
   C->cd(1); gPad->SetGrid(); gPad->SetLogx(); gPad->SetFrameFillColor(33);
   T->Draw("Nation","","hbar2");

   //vertical bar chart
   C->cd(2); gPad->SetGrid(); gPad->SetFrameFillColor(33);
   T->Draw("Division>>hDiv","","goff");
   TH1F *hDiv   = (TH1F*)gDirectory->Get("hDiv");
   hDiv->SetStats(0);
   TH1F *hDivFR = (TH1F*)hDiv->Clone("hDivFR");
   T->Draw("Division>>hDivFR","Nation==\"FR\"","goff");
   hDiv->SetBarWidth(0.45);
   hDiv->SetBarOffset(0.1);
   hDiv->SetFillColor(49);
   TH1 *h1 = hDiv->DrawCopy("bar2");
   hDivFR->SetBarWidth(0.4);
   hDivFR->SetBarOffset(0.55);
   hDivFR->SetFillColor(50);
   TH1 *h2 = hDivFR->DrawCopy("bar2,same");

   TLegend *legend = new TLegend(0.55,0.65,0.76,0.82);
   legend->AddEntry(h1,"All nations","f");
   legend->AddEntry(h2,"French only","f");
   legend->Draw();

   gPad->Modified();
   gPad->Update();

   TestReport1(C, "Ntuple drawing with alphanumeric variables");
   DoCcode(C);
   TestReport2();
}


////////////////////////////////////////////////////////////////////////////////
/// Parallel Coordinates

void parallelcoord()
{
   TCanvas *C = StartTest(800,700);

   TNtuple *ntuple = (TNtuple*)gHsimple->Get("ntuple");

   C->Divide(1,2);

   C->cd(1);
   ntuple->Draw("px:py:pz:random:px*py*pz","","para");
   TParallelCoord* para = (TParallelCoord*)gPad->GetListOfPrimitives()->FindObject("ParaCoord");
   para->SetLineColor(25);
   TColor *col25 = gROOT->GetColor(25);
   if (col25) col25->SetAlpha(0.05);
   C->cd(2);
   ntuple->Draw("px:py:pz:random:px*py*pz","","candle");

   TestReport1(C, "Parallel Coordinates");
   if (gOptionR) printf("%10d%10d\n",0,0);
}


////////////////////////////////////////////////////////////////////////////////
/// Draw a pad and clone it

void clonepad()
{
   TCanvas *C = StartTest(700,500);

   TH1 *hpxpy = (TH1*)gHsimple->Get("hpxpy");
   hpxpy->Draw();
   TCanvas *C2 = (TCanvas*)C->DrawClone();

   TestReport1(C2, "Draw a pad and clone it");
   DoCcode(C2);
   TestReport2();
}


////////////////////////////////////////////////////////////////////////////////
/// Needed for the "waves" test

Double_t interference( Double_t *x, Double_t *par)
{
   Double_t x_p2 = x[0] * x[0];
   Double_t d_2 = 0.5 * par[2];
   Double_t ym_p2 = (x[1] - d_2) * (x[1] - d_2);
   Double_t yp_p2 = (x[1] + d_2) * (x[1] + d_2);
   Double_t  tpi_l = TMath::Pi() /  par[1];
   Double_t amplitude = par[0] * (cos(tpi_l  * sqrt(x_p2 + ym_p2))
                         + par[3] * cos(tpi_l  * sqrt(x_p2 + yp_p2)));
   return amplitude * amplitude;
}


////////////////////////////////////////////////////////////////////////////////
/// Needed for the "waves" test

Double_t result( Double_t *x, Double_t *par)
{
   Double_t xint[2];
   Double_t  maxintens = 0, xcur = 14;
   Double_t dlambda = 0.1 * par[1];
   for(Int_t i=0; i<10; i++){
      xint[0] = xcur;
      xint[1] = x[1];
      Double_t  intens = interference(xint, par);
      if(intens > maxintens) maxintens = intens;
      xcur -= dlambda;
   }
   return maxintens;
}


////////////////////////////////////////////////////////////////////////////////
/// TGraph, TArc, TPalette and TColor

void waves()
{
   TF2 * finter;
   Double_t d = 3;
   Double_t lambda = 1;
   Double_t amp = 10;

   TCanvas *C = StartTest(1004, 759);

   C->Range(0, -10,  30, 10);
   C->SetFillColor(0);
   TPad *pad = new TPad("pr","pr",  0.5, 0 , 1., 1);
   pad->Range(0, -10,  15, 10);
   pad->Draw();

   const Int_t colNum = 30;
   Int_t palette[colNum];
   Int_t color_offset = 2001;
   for (Int_t i=0;i<colNum;i++) {
      new TColor(color_offset+i
      ,    pow(i/((colNum)*1.0),0.3)
      ,    pow(i/((colNum)*1.0),0.3)
      ,0.5*(i/((colNum)*1.0)),"");
      palette[i] = color_offset+i;
   }
   gStyle->SetPalette(colNum,palette);
   C->cd();
   TF2 * f0 = new TF2("ray_source",interference, 0.02, 15, -8, 8, 4);

   f0->SetParameters(amp, lambda, 0, 0);
   f0->SetNpx(200);
   f0->SetNpy(200);
   f0->SetContour(colNum-2);
   f0->Draw("samecolz");

   TLatex title;
   title.DrawLatex(1.6, 8.5, "A double slit experiment");

   TGraph *graph = new TGraph(4);
   graph->SetFillColor(0);
   graph->SetFillStyle(1001);
   graph->SetLineWidth(0);
   graph->SetPoint(0, 0., 0.1);
   graph->SetPoint(1, 14.8, 8);
   graph->SetPoint(2, 0, 8);
   graph->SetPoint(3, 0, 0.1);
   graph->Draw("F");

   graph = new TGraph(4);
   graph->SetFillColor(0);
   graph->SetFillStyle(1001);
   graph->SetLineWidth(0);
   graph->SetPoint(0, 0, -0.1);
   graph->SetPoint(1, 14.8, -8);
   graph->SetPoint(2, 0, -8);
   graph->SetPoint(3, 0, -0.1);
   graph->Draw("F");

   TLine * line;
   line = new TLine(15,-10, 15, 0 - 0.5*d -0.2);
   line->SetLineWidth(10); line->Draw();
   line = new TLine(15, 0 - 0.5*d +0.2 ,15, 0 + 0.5*d -0.2);
   line->SetLineWidth(10); line->Draw();

   line = new TLine(15,0 + 0.5*d + 0.2,15, 10);
   line->SetLineWidth(10); line->Draw();

   pad ->cd();
   finter = new TF2("interference",interference, 0.01, 14, -10, 10, 4);

   finter->SetParameters(amp, lambda, d, 1);
   finter->SetNpx(200);
   finter->SetNpy(200);
   finter->SetContour(colNum-2);
   finter->Draw("samecolorz");

   TArc *arc = new TArc();;
   arc->SetFillStyle(0);
   arc->SetLineWidth(2);
   arc->SetLineColor(5);
   Float_t r = 0.5 * lambda, dr = lambda;
      for (Int_t i = 0; i < 15; i++) {
      arc->DrawArc(0,  0.5*d, r, 0., 360., "only");
      arc->DrawArc(0, -0.5*d, r, 0., 360., "only");
      r += dr;
   }

   pad ->cd();
   TF2 * fresult = new TF2("result",result, 14, 15, -10, 10, 4);

   fresult->SetParameters(amp, lambda, d, 1);
   fresult->SetNpx(300);
   fresult->SetNpy(300);
   fresult->SetContour(colNum-2);
   fresult->Draw("samecolor");
   line = new TLine(13.8,-10, 14, 10);
   line->SetLineWidth(10); line->SetLineColor(0); line->Draw();

   TestReport1(C, "TGraph, TArc, TPalette and TColor");
   if (gOptionR) printf("%10d%10d\n",0,0);
}


////////////////////////////////////////////////////////////////////////////////

void cleanup()
{
}
