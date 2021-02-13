#include <iostream>
#include <fstream>
#include <vector>

#include <cmath>

#include <TMath.h>
#include <Math/SpecFuncMathCore.h>

// #include "SpecFuncCephes.h"

#include <TApplication.h>

#include <TCanvas.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TLegend.h>

const double ERRORLIMIT = 1E-8;
const double MIN = -2.5;
const double MAX = +2.5;
const double INCREMENT = 0.01;
const int ARRAYSIZE = (int) (( MAX - MIN ) / INCREMENT) + 1;

bool showGraphics = false;
bool verbose = false;

using namespace std;

TGraph* drawPoints(Double_t x[], Double_t y[], int color, int style = 1)
{
   TGraph* g = new TGraph(ARRAYSIZE, x, y);
   g->SetLineColor(color);
   g->SetLineStyle(style);
   g->SetLineWidth(3);
   g->Draw("SAME");

   return g;
}

int testSpecFuncErf()
{
   vector<Double_t> x( ARRAYSIZE );
   vector<Double_t> yerf( ARRAYSIZE );
   vector<Double_t> ymerf( ARRAYSIZE );
   vector<Double_t> yerfc( ARRAYSIZE );
   vector<Double_t> ymerfc( ARRAYSIZE );
   vector<Double_t> yierf( ARRAYSIZE );
   vector<Double_t> yierfc( ARRAYSIZE );
//    vector<Double_t> yndtri( ARRAYSIZE );

   int status = 0;

//    ofstream outputFile ("values.txt");

   unsigned int index = 0;
   for ( double i = MIN; i < MAX; i += INCREMENT )
   {
//       outputFile << "i:"; outputFile.width(5); outputFile << i
//            << " index: "; outputFile.width(5); outputFile << index
//            << " TMath::Erf(x): "; outputFile.width(10); outputFile << TMath::Erf(i)
//            << " ROOT::Math::erf(x): "; outputFile.width(10); outputFile << ROOT::Math::erf(i)
//            << " TMath::Erfc(x): "; outputFile.width(10); outputFile << TMath::Erfc(i)
//            << " ROOT::Math::erfc(x): "; outputFile.width(10); outputFile << ROOT::Math::erfc(i)
//            << " TMath::ErfInverse(x): "; outputFile.width(10); outputFile << TMath::ErfInverse(i)
//            << " TMath::ErfcInverse(x): "; outputFile.width(10); outputFile << TMath::ErfcInverse(i)
// //            << " ROOT::Math::Cephes::ndtri(x): "; outputFile.width(10); outputFile << ROOT::Math::Cephes::ndtri(i)
//            << endl;

      x[index] = i;

      yerf[index] = TMath::Erf(i);
      ymerf[index] = ROOT::Math::erf(i);
      if ( std::fabs( yerf[index] - ymerf[index] ) > ERRORLIMIT )
      {
         cout << "i " << i
              << " yerf[index] " << yerf[index]
              << " ymerf[index] " << ymerf[index]
              << " " << std::fabs( yerf[index] - ymerf[index] )
              << endl;
         status += 1;
      }

      yerfc[index] = TMath::Erfc(i);
      ymerfc[index] = ROOT::Math::erfc(i);
      if ( std::fabs( yerfc[index] - ymerfc[index] ) > ERRORLIMIT )
      {
         cout << "i " << i
              << " yerfc[index] " << yerfc[index]
              << " ymerfc[index] " << ymerfc[index]
              << " " << std::fabs( yerfc[index] - ymerfc[index] )
              << endl;
         status += 1;
      }

      yierf[index] = TMath::ErfInverse(yerf[index]);
      if ( std::fabs( yierf[index] - i ) > ERRORLIMIT )
      {
         cout << "i " << i
              << " yierf[index] " << yierf[index]
              << " " << std::fabs( yierf[index] - i )
              << endl;
         status += 1;
      }

      yierfc[index] = TMath::ErfcInverse(yerfc[index]);
      if ( std::fabs( yierfc[index] - i ) > ERRORLIMIT )
      {
         cout << "i " << i
              << " yierfc[index] " << yierfc[index]
              << " " << std::fabs( yierfc[index] - i )
              << endl;
         status += 1;
      }

//       yndtri[index] = ROOT::Math::Cephes::ndtri(i);

      index += 1;
   }

   if ( showGraphics )
   {

      TCanvas* c1 = new TCanvas("c1", "Two Graphs", 600, 400);
      TH2F* hpx = new TH2F("hpx", "Two Graphs(hpx)", ARRAYSIZE, MIN, MAX, ARRAYSIZE, -1,2);
      hpx->SetStats(kFALSE);
      hpx->Draw();

      TGraph* gerf   = drawPoints(&x[0], &yerf[0], 14);
      TGraph* gmerf  = drawPoints(&x[0], &ymerf[0], 5, 7);
      TGraph* gerfc  = drawPoints(&x[0], &yerfc[0], 2);
      TGraph* gmerfc = drawPoints(&x[0], &ymerfc[0], 3, 7);
//   drawPoints(&x[0], &yierf[0], 21);
//   drawPoints(&x[0], &yierfc[0], 28);
//   drawPoints(&x[0], &yndtri[0], 9);

      TLegend* legend = new TLegend(0.61,0.62,0.86,0.86);
      legend->AddEntry(gerf,   "TMath::Erf()");
      legend->AddEntry(gmerf,  "ROOT:Math::erf()");
      legend->AddEntry(gerfc,  "TMath::Erfc()");
      legend->AddEntry(gmerfc, "ROOT::Math::erfInverse()");
      legend->Draw();

      c1->Show();
   }

   cout << "Test Done!" << endl;

   return status;
}

int main(int argc, char **argv)
{

  // Parse command line arguments
  for (Int_t i=1 ;  i<argc ; i++) {
     std::string arg = argv[i] ;
     if (arg == "-g") {
      showGraphics = true;
     }
     if (arg == "-v") {
      showGraphics = true;
      //verbose = true;
     }
     if (arg == "-h") {
        cerr << "Usage: " << argv[0] << " [-g] [-v]\n";
        cerr << "  where:\n";
        cerr << "     -g : graphics mode\n";
        cerr << "     -v : verbose  mode";
        cerr << endl;
        return -1;
     }
   }

   TApplication* theApp = 0;
   if ( showGraphics )
      theApp = new TApplication("App",&argc,argv);

   int status = testSpecFuncErf();

   if ( showGraphics )
   {
      theApp->Run();
      delete theApp;
      theApp = 0;
   }

   return status;
}
