// Compare the AdaptiveIntegratorMultiDim and
// TF1::IntegralMultiple performance and results

// Compares time performance
// for different dimensions
// draws a graph
//
// Author: David Gonzalez Maline
//

#include "TMath.h"
#include "TStopwatch.h"
#include <cmath>
#include <iostream>

#include "Math/Integrator.h"
#include "Math/Functor.h"
#include "Math/IFunction.h"
#include "Math/WrappedParamFunction.h"
#include "Math/AdaptiveIntegratorMultiDim.h"
#include "Math/IFunctionfwd.h"
#include "TF1.h"

// for graphical comparison of performance
#include "TAxis.h"
#include "TVirtualPad.h"
#include "TApplication.h"
#include "TPaveLabel.h"
#include "TLegend.h"
#include "TH1.h"

bool showGraphics = false;
bool verbose = false;

using namespace std;

int NMAX = 6; //maximum dimension

Double_t Sum( const double* x, const double *p)
{
  double sum = 0.;
  for(int i = 0; i < p[0]; i++)
    sum += x[i];
  return  sum;

}

//multidim function to integrate
Double_t SimpleFun( const double* x, const double *p)
{
  double prod = 1.;
  for(int i = 0; i < p[0]; i++)
    prod *= TMath::Power(TMath::E(), -x[i]*x[i]);
  return  prod*Sum(x, p)*TMath::Sin(x[0]);

}

  // ################################################################
  //
  //      testing IntegratorMultiDim class
  //
  // ################################################################

double integral_num(unsigned int dim, double* a, double* b, double* p)
{

  if (verbose) std::cout << "\nTesting IntegratorMultiDim class.." << std::endl;
  //std::cout << p[0] << std::endl;dimensionality

  TStopwatch timer;
  timer.Start();
  ROOT::Math::WrappedParamFunction<> funptr1(&SimpleFun, dim, p, p+1);
  unsigned int nmax = (unsigned int) 1.E7; // apply cut off to avoid routine to explode
  ROOT::Math::AdaptiveIntegratorMultiDim ig1(funptr1, 1.E-5, 1.E-5, nmax);

  ig1.SetFunction(funptr1);
  ig1.Integral(a, b);
  timer.Stop();

  if (verbose) {
     std::cout.precision(12);
     std::cout << "result:  \t";
     std::cout << ig1.Result() << "\t" << "error: \t" << ig1.Error() << std::endl;
     std::cout << "Number of function evaluations: " << ig1.NEval() << std::endl;
     std::cout << "Time using IntegratorMultiDim: \t" << timer.RealTime() << std::endl;
     std::cout << "------------------------------------" << std::endl;
  }
  return timer.RealTime();
}

  // ################################################################
  //
  //      testing TF1::IntegralMultiple class
  //
  // ################################################################
double integral_TF1(unsigned int dim, double* a, double* b, double* p)
{

  double timeTF1;
  if (verbose) std::cout << "\ntesting TF1::IntegralMultiple.." << std::endl;

  TStopwatch timer;
  //timer.Start();
  //ROOT::Math::WrappedParamFunction<> funptr(&SimpleFun, dim, p, p+1);

  TF1 function("function", &SimpleFun, 0, 0, 1);
  function.SetParameters(p);
  double error, result;
  int iter, fail;
  result = function.IntegralMultiple(dim, a, b, 0, (int) 1.E7, 1.E-5, error, iter, fail);

//   ROOT::Math::GSLMCIntegrator ig1;
//   ig1.SetType(ROOT::Math::MCIntegration::VEGAS);
//   ig1.SetFunction(funptr);
//   ig1.Integral(a, b);

  timer.Stop();
//  std::cout << "result: \t";
//  std::cout << ig1.Result() << "\t" << "error: \t" << ig1.Error() << std::endl;
//   std::cout << "sigma: \t" << ig1.Sigma();
//   std::cout << "\t" << "chi2: \t" << ig1.ChiSqr() << std::endl;
//   std::cout << "\nTime using TF1::IntegralMultiple :\t" << timer.RealTime() << std::endl;
//   std::cout << "\n------------------------------------" << std::endl;
//   std::cout << "\t MISER.. \n" << std::endl;

  if (verbose) {
     std::cout.precision(12);
     std::cout << "result:  \t";
     std::cout << result << "\t" << "error: \t" << error << std::endl;
     std::cout << "Number of function evaluations: " << iter << std::endl;
     std::cout << "Time using TF1::IntegralMultiple : \t" << timer.RealTime() << std::endl;
     std::cout << "------------------------------------" << std::endl;
  }

  timeTF1 = timer.RealTime();

  return timeTF1;
}

void performance()
{
  //dimensionality
  unsigned int Nmax = NMAX;
  unsigned int size = Nmax-1;
  TH1D *num_performance = new TH1D("cubature", "", size, 1.5, Nmax+.5);
  TH1D *TF1_performance = new TH1D("montecarlo", "", size, 1.5, Nmax+.5);

   num_performance->SetBinContent(1, 0.0);
   TF1_performance->SetBinContent(1,0.0);
   std::cout << "Testing multidim integration performances for various dimensions\n";
   for(unsigned int N = 2; N <=Nmax; N++)//dim
   {
      if (verbose)  {
         std::cout<< "*********************************************" << std::endl;
         std::cout<< "Number of dimensions: "<< N << std::endl;
      }
      else {
         std::cout << "\tdim="<< N << std::endl;
      }
      //integration limits
      double * a = new double[N];
      double * b = new double[N];
      double p[1];
      p[0] = N;
      for (unsigned int i=0; i < N; i++)
      {
         a[i] = -1.;//-TMath::Pi();
         b[i] = 1;//TMath::Pi();
      }
      num_performance->SetBinContent(N-1, integral_num(N, a, b, p));
      TF1_performance->SetBinContent(N-1,integral_TF1(N, a, b, p));

      delete [] a;
      delete [] b;
   }

   if (showGraphics) {
      num_performance->SetBarWidth(0.45);
      num_performance->SetBarOffset(0.05);
      num_performance->SetFillColor(49);
      num_performance->SetStats(0);
      //num_performance->GetXaxis()->SetLimits(1.5, Nmax+0.5);
      num_performance->GetXaxis()->SetTitle("number of dimensions");
      num_performance->GetYaxis()->SetTitle("time [s]");
      num_performance->SetTitle("comparison of performance");
      TH1 *h1 = num_performance->DrawCopy("bar2");
      TF1_performance->SetBarWidth(0.40);
      TF1_performance->SetBarOffset(0.5);
      TF1_performance->SetFillColor(kRed);
      TH1 *h2 =  TF1_performance->DrawCopy("bar2,same");

      TLegend *legend = new TLegend(0.25,0.65,0.55,0.82);
      legend->AddEntry(h1,"AdaptiveIntegratorMultiDim","f");
      legend->AddEntry(h2,"TF1::IntegralMultiple()","f");
      legend->Draw();

      gPad->Update();
   }
   std::cout << "Test Timing results\n";
   std::cout << "N dim \t Adaptive time \t\t TF1 time \n";
   for (unsigned int i=1; i<=size; i++)
      std::cout << i+1 << " " << num_performance->GetBinContent(i) << "\t" << TF1_performance->GetBinContent(i)<<std::endl;

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

   performance();

   if ( showGraphics )
   {
      theApp->Run();
      delete theApp;
      theApp = 0;
   }

   return 0;

}
