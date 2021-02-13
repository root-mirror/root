// test of multidimentional Integration
// Calculates an intergal of a function
// in 2,3,..., 8 dimensions
// by using adaptive Genz Malik cubature
// and MonteCarlo methods:
// --PLAIN
// --VEGAS
// --MISER
//
// from
// IntegratorMultiDim class
// and GSLMCIntegrator class
//
// Compares time performance
// for different dimensions
// draws a graph
//
// Author: Magdalena Slawinska
//

#include "TMath.h"
#include "TStopwatch.h"

#include <cmath>
#include <iostream>
#include <iomanip>

#include "Math/Integrator.h"
#include "Math/Functor.h"
#include "Math/IFunction.h"
#include "Math/WrappedParamFunction.h"
#include "Math/AdaptiveIntegratorMultiDim.h"
#include "Math/IFunctionfwd.h"
#include "Math/GSLMCIntegrator.h"
#include "Math/Random.h"
#include "Math/GSLRndmEngines.h"


// for graphical comparison of performance
#include "TAxis.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TPaveLabel.h"
#include "TLegend.h"
#include "TH1.h"

bool showGraphics = false;
bool verbose = false;
int NMAX = 6; // maximum dimension used


using std::cout;
using std::endl;
using std::cerr;

//const int n = 3; //default dimensionality

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

//singularity at (0,0, ..., 0)
Double_t SingularFun( const double* x, const double *p)
{

  double prod = 1.;
  for(int i = 0; i < p[0]; i++)
    prod *=    TMath::Cos(x[i]);
  return 1./((1-prod)*8*TMath::Power(TMath::Pi(), 3));

}


  // ################################################################
  //
  //      testing IntegratorMultiDim class
  //
  // ################################################################

double integral_num(unsigned int dim, double* a, double* b, double* p, double & value, double & error)
{
  if (verbose) std::cout << "\nTesting IntegratorMultiDim class.." << std::endl;
  //std::cout << p[0] << std::endl;dimensionality

  TStopwatch timer;
  timer.Start();
  ROOT::Math::WrappedParamFunction<> funptr1(&SimpleFun, dim, p, p+1);
  unsigned int nmax = (unsigned int) 1.E8; // apply cut off to avoid routine to explode
  //unsigned int size = (unsigned int) 1.E8; // apply cut off to avoid routine to explode
  unsigned int size = 0; // use default  value defined by nmax
  ROOT::Math::AdaptiveIntegratorMultiDim ig1(funptr1, 1.E-5, 1.E-5, nmax,size);
  //  std::cout << "1. integral= " << std::endl;
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
  else { std::cout << " . "; }
  value = ig1.Result();
  error = ig1.Error();
  return timer.RealTime();
}

  // ################################################################
  //
  //      testing MCIntegrator class
  //
  // ################################################################
std::vector<double> integral_MC(unsigned int dim, double* a, double* b, double* p, double * value, double * error)
{

  double timeVegas;
  double timeMiser;
  double timePlain;
  if (verbose) {
     std::cout << "\nTesting GSLMCIntegrator class.." << std::endl;
     std::cout << "\t VEGAS.. " << std::endl;
     std::cout << "" << std::endl;
  }
  else { std::cout << "."; }

  TStopwatch timer;
  //timer.Start();
  ROOT::Math::WrappedParamFunction<> funptr(&SimpleFun, dim, p, p+1);

  ROOT::Math::GSLMCIntegrator ig1(ROOT::Math::MCIntegration::kVEGAS);
  ig1.SetFunction(funptr);
  //ig1.SetMode(ROOT::Math::MCIntegration::kIMPORTANCE_ONLY);

  /*
  VegasParameters param;
  param.iterations = 2;
  ig1.SetParameters(param);
  */

  ig1.Integral(a, b);
  timer.Stop();
  timeVegas = timer.RealTime();
  value[0] = ig1.Result();
  error[0] = ig1.Error();

  if (verbose) {
     std::cout << "result: \t";
     std::cout << ig1.Result() << "\t" << "error: \t" << ig1.Error() << std::endl;
     std::cout << "sigma: \t" << ig1.Sigma();
     std::cout << "\t" << "chi2: \t" << ig1.ChiSqr() << std::endl;
     std::cout << std::endl;
     std::cout << "Time using GSLMCIntegrator::VEGAS :\t" << timer.RealTime() << std::endl;
     //std::cout << "Number of function evaluations: " << ig1.Eval() << std::endl;
     //ig2_param.iterations = 1000;
     std::cout << "" <<std::endl;
     std::cout << "------------------------------------" << std::endl;
     std::cout << "\t MISER.. " << std::endl;
     std::cout << "" << std::endl;
  }
  else { std::cout << "."; }


  timer.Start();
  ROOT::Math::GSLMCIntegrator ig2(ROOT::Math::MCIntegration::kMISER);

  ig2.SetFunction(funptr);

  // test using a different generator
  ROOT::Math::Random<ROOT::Math::GSLRngCMRG> r;
  ig2.SetGenerator(r.Rng() );


  //par.min_calls = 4*dim;
  //par.min_calls_per_bisection = 8*par.min_calls;


  //MiserParameters par(dim);
  //ig2.SetParameters(par);
  ig2.Integral(a, b);
  timer.Stop();
  timeMiser = timer.RealTime();
  value[1] = ig2.Result();
  error[1] = ig2.Error();
  if (verbose) {
     std::cout << "result: \t";
     std::cout << ig2.Result() << "\t" << "error: \t" << ig2.Error() << std::endl;

     std::cout << "Time using GSLMCIntegrator::MISER :\t" << timer.RealTime() << std::endl;
     std::cout << "" << std::endl;
     std::cout << "------------------------------------" << std::endl;
     std::cout << "\t PLAIN.. " << std::endl;
  }
  else { std::cout << "."; }

  timer.Start();

  ROOT::Math::GSLMCIntegrator ig3(ROOT::Math::MCIntegration::kPLAIN);
  ig3.SetFunction(funptr);
  ig3.Integral(a, b);
  timer.Stop();
  timePlain = timer.RealTime();
  value[2] = ig3.Result();
  error[2] = ig3.Error();

  if (verbose) {
     std::cout << "" << std::endl;
     std::cout << "result: \t";
     std::cout << ig3.Result() << "\t" << "error: \t" << ig3.Error() << std::endl;
     std::cout << "Time using GSLMCIntegrator::PLAIN :\t" << timer.RealTime() << std::endl;
     std::cout << "" << std::endl;
  }
  std::vector<double> result(3);
  result[0] = timeVegas; result[1] = timeMiser; result[2] = timePlain;
  return result;
}

bool performance()
{
  //dimensionality
  unsigned int Nmax = NMAX;
  unsigned int size = Nmax-1;
  bool ok = true;

  TH1D *num_performance = new TH1D("cubature", "", size, 1.5, Nmax+.5);
  TH1D *Vegas_performance = new TH1D("VegasMC", "", size, 1.5, Nmax+.5);
  TH1D *Miser_performance = new TH1D("MiserMC", "", size, 1.5, Nmax+.5);
  TH1D *Plain_performance = new TH1D("PlainMC", "", size, 1.5, Nmax+.5);

   num_performance->SetBinContent(1, 0.0);
   Vegas_performance->SetBinContent(1,0.0);
   for(unsigned int N = 2; N <=Nmax; N++)//dim
  {
     if (verbose) {
        std::cout<< "*********************************************" << std::endl;
        std::cout<< "Number of dimensions: "<< N << std::endl;
     }
     else {
        std::cout << "\n\tdim="<< N << " : ";
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
    //x[N] = N;
    double val0, err0;
    double valMC[3], errMC[3];
    double timeNumInt = integral_num(N, a, b, p, val0, err0);
    std::vector<double> timeMCInt = integral_MC(N, a, b, p, valMC, errMC);

    // set the histograms
    num_performance->SetBinContent(N-1, timeNumInt);
    Vegas_performance->SetBinContent(N-1, timeMCInt[0]);
    Miser_performance->SetBinContent(N-1, timeMCInt[1]);
    Plain_performance->SetBinContent(N-1, timeMCInt[2]);

    // test the values
    for (int j = 0; j < 3; ++j) {
       if (TMath::Abs(val0-valMC[j] ) > 5 * std::sqrt( err0*err0 + errMC[j]*errMC[j] ) ) {
          Error("testMCIntegration","Result is not consistent for dim %d between adaptive and MC %d ",N,j);
          ok = false;
       }
    }

    delete [] a;
    delete [] b;

   }




   if ( showGraphics )
   {

      TCanvas * c1 = new TCanvas();
      c1->SetFillColor(kYellow-10);

      num_performance->SetBarWidth(0.23);
      num_performance->SetBarOffset(0.04);
      num_performance->SetFillColor(kRed+3);
      num_performance->SetStats(0);
      //num_performance->GetXaxis()->SetLimits(1.5, Nmax+0.5);
      num_performance->GetXaxis()->SetTitle("number of dimensions");
      num_performance->GetYaxis()->SetTitle("time [s]");
      num_performance->SetTitle("comparison of performance");
      TH1 *h1 = num_performance->DrawCopy("bar");
      Vegas_performance->SetBarWidth(0.23);
      Vegas_performance->SetBarOffset(0.27);
      Vegas_performance->SetFillColor(kRed);
      TH1 *h2 =  Vegas_performance->DrawCopy("bar,same");
      Miser_performance->SetBarWidth(0.23);
      Miser_performance->SetBarOffset(0.50);
      Miser_performance->SetFillColor(kRed-7);
      TH1 *h3 =  Miser_performance->DrawCopy("bar,same");
      Plain_performance->SetBarWidth(0.23);
      Plain_performance->SetBarOffset(0.73);
      Plain_performance->SetFillColor(kRed-10);
      TH1 *h4 =  Plain_performance->DrawCopy("bar,same");

      TLegend *legend = new TLegend(0.25,0.65,0.55,0.82);
      legend->AddEntry(h1,"Cubature","f");
      legend->AddEntry(h2,"MC Vegas","f");
      legend->AddEntry(h3,"MC Miser","f");
      legend->AddEntry(h4,"MC Plain","f");
      legend->Draw();
      gPad->SetLogy(true);
      gPad->Update();
   }

   std::cout << "\nTest Timing results\n";
   std::cout << "   N dim   \t     Adaptive    MC Vegas    MC Miser    MC Plain \n";
   for (unsigned int i=1; i<=size; i++) {
      std::cout.width(8);
      std::cout.precision(6);
      std::cout << i+1;
      std::cout << "\t " << std::setw(12) << num_performance->GetBinContent(i) << std::setw(12)       << Vegas_performance->GetBinContent(i)
                << std::setw(12) << Miser_performance->GetBinContent(i)   << std::setw(12) << Plain_performance->GetBinContent(i) << std::endl;
   }
   return ok;
}


int main(int argc, char **argv)
{
   int status = 0;

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
        cout << "Usage: " << argv[0] << " [-g] [-v]\n";
        cout << "  where:\n";
        cout << "     -g : graphics mode\n";
        cout << "     -v : verbose  mode";
        cout << endl;
        return -1;
     }
   }


   TApplication* theApp = 0;

   if ( showGraphics )
      theApp = new TApplication("App",&argc,argv);

   status = performance() ? 0 : 1;

   if ( showGraphics )
   {
      theApp->Run();
      delete theApp;
      theApp = 0;
   }

   return status;
}
