#include <iostream>
#include <vector>

#include <cmath>

#include <TMath.h>
#include <Math/SpecFuncMathCore.h>

#include <TApplication.h>

#include <TCanvas.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TLegend.h>

const double ERRORLIMIT = 1E-8;
const double MIN = 0;
const double MAX = 1;
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

int testSpecFuncBetaI()
{
   vector<Double_t> x( ARRAYSIZE );
   vector<Double_t> yb( ARRAYSIZE );
   vector<Double_t> ymb( ARRAYSIZE );

   int status = 0;

   double b = 0.2, a= 0.9;
   cout << "** b = " << b << " **" << endl;
   unsigned int index = 0;
   for ( double i = MIN; i < MAX; i += INCREMENT )
   {
      if (verbose) {
         cout << "i:"; cout.width(5); cout << i
              << " index: "; cout.width(5); cout << index
              << " TMath::BetaIncomplete(x,a,b): "; cout.width(10); cout << TMath::BetaIncomplete(i,a,b)
              << " ROOT::Math::inc_beta(a,a,b): "; cout.width(10); cout << ROOT::Math::inc_beta(i,a,b)
              << endl;
      }

      x[index] = i;

      yb[index] = TMath::BetaIncomplete(i,a,b);
      ymb[index] = ROOT::Math::inc_beta(i,a,b);
      if ( std::fabs( yb[index] - ymb[index] ) > ERRORLIMIT )
      {
         cout << "i " << i
              << " yb[index] " << yb[index]
              << " ymb[index] " << ymb[index]
              << " " << std::fabs( yb[index] - ymb[index] )
              << endl;
         status += 1;
      }
      index += 1;
   }

   if ( showGraphics )
   {

      TCanvas* c1 = new TCanvas("c1", "Two Graphs", 600, 400);
      TH2F* hpx = new TH2F("hpx", "Two Graphs(hpx)", ARRAYSIZE, MIN, MAX, ARRAYSIZE, 0, 5);
      hpx->SetStats(kFALSE);
      hpx->Draw();

      int color = 2;

      TGraph *gb, *gmb;
      gb = drawPoints(&x[0], &yb[0], color++);
      gmb = drawPoints(&x[0], &ymb[0], color++, 7);

      TLegend* legend = new TLegend(0.61,0.72,0.86,0.86);
      legend->AddEntry(gb, "TMath::BetaIncomplete()");
      legend->AddEntry(gmb, "ROOT::Math::inc_beta()");
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
      verbose = true;
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

   int status = testSpecFuncBetaI();

   if ( showGraphics )
   {
      theApp->Run();
      delete theApp;
      theApp = 0;
   }

   return status;
}
