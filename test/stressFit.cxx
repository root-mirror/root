// @(#)root/test:$name:  $:$id: stressFit.cxx,v 1.15 2002/10/25 10:47:51 rdm exp $
// Authors: Rene Brun, Eddy Offermann  April 2006

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*//
//                                                                           //
// Function Minimization Examples, Fred James                                //
//                                                                           //
// from the                                                                  //
//   Proceedings of the 1972 CERN Computing and Data Processing School       //
//   Pertisau, Austria, 10-24 September, 1972 (CERN 72-21)                   //
//                                                                           //
// Here a collection of test problems is assembled which were found to be    //
// useful in verifying and comparing minimization routines. Many of these    //
// are standard functions upon which it has become conventional to try all   //
// new methods, quoting the performance in the publication of the algorithm  //
//                                                                           //
// Each test will produce one line (Test OK or Test FAILED) . At the end of  //
// the test a table is printed showing the global results Real Time and      //
// Cpu Time. One single number (ROOTMARKS) is also calculated showing the    //
// relative performance of your machine compared to a reference machine      //
// a Pentium IV 2.4 Ghz) with 512 MBytes of memory and 120 GBytes IDE disk.  //
//                                                                           //
// In the main routine the fitter can be chosen through TVirtualFitter :     //
//   - Minuit                                                                //
//   - Minuit2                                                               //
//   - Fumili                                                                //
//
//  To run the test, do, eg
// root -b -q stressFit.cxx
// root -b -q "stressFit.cxx(\"Minuit2\")"
// root -b -q "stressFit.cxx+(\"Minuit2\")"
//                                                                           //
// The verbosity can be set through the global parameter gVerbose :          //
//   -1: off  1: on                                                          //
// The tolerance on the parameter deviation from the minimum can be set      //
// through gAbsTolerance .                                                   //
//                                                                           //
// An example of output when all the tests run OK is shown below:            //
// *******************************************************************       //
// *  Minimization - S T R E S S suite                               *       //
// *******************************************************************       //
// *******************************************************************       //
// *  Starting  S T R E S S                                          *       //
// *******************************************************************       //
// Test  1 : Wood.................................................. OK       //
// Test  2 : RosenBrock............................................ OK       //
// Test  3 : Powell................................................ OK       //
// Test  4 : Fletcher.............................................. OK       //
// Test  5 : GoldStein1............................................ OK       //
// Test  6 : GoldStein2............................................ OK       //
// Test  7 : TrigoFletcher......................................... OK       //
// *******************************************************************       //
//                                                                           //
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*_*//

#include <cstdlib>
#include "TSystem.h"
#include "TROOT.h"
#include "TBenchmark.h"
#include "TMath.h"
#include "TStopwatch.h"
#include "Riostream.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "snprintf.h"

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/IFunction.h"
#include "Math/MinimizerOptions.h"
#include "Math/Minimizer.h"

Int_t stressFit(const char *type = "Minuit", const char *algo = "Migrad", Int_t N = 2000);
Int_t    gVerbose      = -1;
Double_t gToleranceMult = 1.e-3;

////////////////////////////////////////////////////////////////////////////////
/// Print test program number and its title

void StatusPrint(Int_t id,const TString &title, Int_t nsuccess, Int_t nattempts)
{
  const Int_t kMAX = 65;
  Char_t number[4];
  snprintf(number,4,"%2d",id);
  TString header = TString("Test ")+number+" : "+title;
  const Int_t nch = header.Length();
  for (Int_t i = nch; i < kMAX; i++) header += '.';
  std::cout << header << " " << nsuccess << " out of " << nattempts << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

Double_t RosenBrock(const Double_t *par)
{
  const Double_t x = par[0];
  const Double_t y = par[1];
  const Double_t tmp1 = y-x*x;
  const Double_t tmp2 = 1-x;
  return 100*tmp1*tmp1+tmp2*tmp2;
}

////////////////////////////////////////////////////////////////////////////////
///
/// F(x,y) = 100 (y-x^2)^2 + (1-x)^2
///
///   start point: F(-1.2,1.0) = 24.20
///   minimum    : F(1.0,1.0)  = 0.
///
/// This narrow, parabolic valley is probably the best known of all test cases. The floor
/// of the valley follows approximately the parabola y = x^2+1/200 .
/// There is a region where the covariance matrix is not positive-definite and even a path
/// where it is singular . Stepping methods tend to perform at least as well as gradient
///  method for this function .
/// [Reference: Comput. J. 3,175 (1960).]

Bool_t RunRosenBrock()
{
  Bool_t ok = kTRUE;
  const int nvars = 2;

  ROOT::Math::Minimizer* min =
    ROOT::Math::Factory::CreateMinimizer(ROOT::Math::MinimizerOptions::DefaultMinimizerType(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo());
  min->SetPrintLevel(gVerbose);

  ROOT::Math::Functor f(&RosenBrock, nvars);
  min->SetFunction(f);
  double step = 0.01;
  double xmin = -5.0;
  double xmax =  5.0;
  double tolerance = 1.e-3;

  min->SetVariable(0, "x", -1.2, step);
  min->SetVariable(1, "y",  1.0, step);
  for (int ivar = 0; ivar < nvars; ivar++)
    min->SetVariableInitialRange(ivar, xmin, xmax);

  min->SetMaxFunctionCalls(100000);
  min->SetTolerance(tolerance * gToleranceMult);

  min->Minimize();

  if (min->MinValue() > tolerance)
    ok = kFALSE;

  delete min;

  return ok;
}

////////////////////////////////////////////////////////////////////////////////

Double_t Wood4(const Double_t *par)
{
  const Double_t w = par[0];
  const Double_t x = par[1];
  const Double_t y = par[2];
  const Double_t z = par[3];

  const Double_t w1 = w-1;
  const Double_t x1 = x-1;
  const Double_t y1 = y-1;
  const Double_t z1 = z-1;
  const Double_t tmp1 = x-w*w;
  const Double_t tmp2 = z-y*y;

  return 100*tmp1*tmp1+w1*w1+90*tmp2*tmp2+y1*y1+10.1*(x1*x1+z1*z1)+19.8*x1*z1;
}

////////////////////////////////////////////////////////////////////////////////
///
/// F(w,x,y,z) = 100 (y-w^2)^2 + (w-1)^2 + 90 (z-y^2)^2
///              + (1-y)^2 + 10.1 [(x-1)^2 + (z-1)^2]
///              + 19.8 (x-1)(z-1)
///
///   start point: F(-3,-1,-3,-1) = 19192
///   minimum    : F(1,1,1,1)  =   0.
///
/// This is a fourth-degree polynomial which is reasonably well-behaved near the minimum,
/// but in order to get there one must cross a rather flat, four-dimensional "plateau"
/// which often causes minimization algorithm to get "stuck" far from the minimum. As
/// such it is a particularly good test of convergence criteria and simulates quite well a
/// feature of many physical problems in many variables where no good starting
/// approximation is known .
/// [Reference: Unpublished. See IBM Technical Report No. 320-2949.]

Bool_t RunWood4()
{

  Bool_t ok = kTRUE;
  const int nvars = 4;
  ROOT::Math::Minimizer* min =
    ROOT::Math::Factory::CreateMinimizer(ROOT::Math::MinimizerOptions::DefaultMinimizerType(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo());
  if (!min)
    {
      std::cerr << "RunWood4(): failed to create ROOT::Math::Minimizer" << std::endl;
      return 0;
    }
  min->SetPrintLevel(gVerbose);

  ROOT::Math::Functor f(&Wood4, nvars);
  min->SetFunction(f);
  double step = 0.01;
  double xmin = -5.0;
  double xmax =  5.0;
  double tolerance = 1.e-3;

  min->SetVariable(0, "w", -3.0, step);
  min->SetVariable(1, "x", -1.0, step);
  min->SetVariable(2, "y", -3.0, step);
  min->SetVariable(3, "z", -1.0, step);
  for (int ivar = 0; ivar < nvars; ivar++)
    min->SetVariableInitialRange(ivar, xmin, xmax);

  min->SetMaxFunctionCalls(100000);
  min->SetTolerance(tolerance * gToleranceMult);

  min->Minimize();

  if (min->MinValue() > tolerance)
    ok = kFALSE;
  delete min;

  return ok;
}

////////////////////////////////////////////////////////////////////////////////

Double_t Powell(const Double_t *par)
{
  const Double_t w = par[0];
  const Double_t x = par[1];
  const Double_t y = par[2];
  const Double_t z = par[3];

  const Double_t tmp1 = w+10*x;
  const Double_t tmp2 = y-z;
  const Double_t tmp3 = x-2*y;
  const Double_t tmp4 = w-z;

  return tmp1*tmp1+5*tmp2*tmp2+tmp3*tmp3*tmp3*tmp3+10*tmp4*tmp4*tmp4*tmp4;
}

////////////////////////////////////////////////////////////////////////////////
///
/// F(w,x,y,z) = (w+10x)^2 + 5(y-z)^2 + (x-2y)^4+ 10 (w-z)^4
///
///   start point: F(-3,-1,0,1) = 215
///   minimum    : F(0,0,0,0)  =   0.
///
/// This function is difficult because its matrix of second derivatives becomes singular
///  at the minimum. Near the minimum the function is given by (w + 10x)^2 + 5 (y-5)^2
/// which does not determine the minimum uniquely.
/// [Reference: Comput. J. 5, 147 (1962).]

Bool_t RunPowell()
{
  Bool_t ok = kTRUE;
  const int nvars = 4;

  ROOT::Math::Minimizer* min =
    ROOT::Math::Factory::CreateMinimizer(ROOT::Math::MinimizerOptions::DefaultMinimizerType(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo());
  min->SetPrintLevel(gVerbose);

  ROOT::Math::Functor f(&Powell, nvars);
  min->SetFunction(f);
  double step = 0.01;
  double xmin = -5.0;
  double xmax =  5.0;
  double tolerance = 1.e-3;

  min->SetVariable(0, "w", +3.0, step);
  min->SetVariable(1, "x", -1.0, step);
  min->SetVariable(2, "y",  0.0, step);
  min->SetVariable(3, "z", +1.0, step);
  for (int ivar = 0; ivar < nvars; ivar++)
    min->SetVariableInitialRange(ivar, xmin, xmax);

  min->SetMaxFunctionCalls(100000);
  min->SetTolerance(tolerance * gToleranceMult);

  min->Minimize();

  if (min->MinValue() > tolerance)
    ok = kFALSE;

  delete min;

  return ok;
}

////////////////////////////////////////////////////////////////////////////////

Double_t Fletcher(const Double_t *par)
{
  const Double_t x = par[0];
  const Double_t y = par[1];
  const Double_t z = par[2];

  Double_t psi;
  if (x > 0)
    psi = TMath::ATan(y/x)/2/TMath::Pi();
  else if (x < 0)
    psi = 0.5+TMath::ATan(y/x)/2/TMath::Pi();
  else
    psi = 0.0;

  const Double_t tmp1 = z-10*psi;
  const Double_t tmp2 = TMath::Sqrt(x*x+y*y)-1;

  return 100*(tmp1*tmp1+tmp2*tmp2)+z*z;
}

////////////////////////////////////////////////////////////////////////////////
///
/// F(x,y,z) = 100 {[z - 10 G(x,y)]^2 + ( (x^2+y^2)^1/2 - 1 )^2} + z^2
///
///                     | arctan(y/x)        for x > 0
/// where 2 pi G(x,y) = |
///                     | pi + arctan(y/x)   for x < 0
///
///   start point: F(-1,0,0) = 2500
///   minimum    : F(1,0,0)  =   0.
///
/// F is defined only for -0.25 < G(x,y) < 0.75
///
/// This is a curved valley problem, similar to Rosenbrock's, but in three dimensions .
/// [Reference: Comput. J. 6, 163 (1963).]

Bool_t RunFletcher()
{
  Bool_t ok = kTRUE;
  const int nvars = 3;

  ROOT::Math::Minimizer* min =
    ROOT::Math::Factory::CreateMinimizer(ROOT::Math::MinimizerOptions::DefaultMinimizerType(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo());
  min->SetPrintLevel(gVerbose);

  ROOT::Math::Functor f(&Fletcher, nvars);
  min->SetFunction(f);
  double step = 0.01;
  double xmin = -5.0;
  double xmax =  5.0;
  double tolerance = 1.e-3;

  min->SetVariable(0, "x", -1.0, step);
  min->SetVariable(1, "y",  0.0, step);
  min->SetVariable(2, "z",  0.0, step);
  for (int ivar = 0; ivar < nvars; ivar++)
    min->SetVariableInitialRange(ivar, xmin, xmax);

  min->SetMaxFunctionCalls(100000);
  min->SetTolerance(tolerance * gToleranceMult);

  min->Minimize();

  if (min->MinValue() > tolerance)
    ok = kFALSE;

  delete min;

  return ok;
}

////////////////////////////////////////////////////////////////////////////////

Double_t GoldStein1(const Double_t *par)
{
  const Double_t x = par[0];
  const Double_t y = par[1];

  const Double_t tmp1 = x+y+1;
  const Double_t tmp2 = 19-14*x+3*x*x-14*y+6*x*y+3*y*y;
  const Double_t tmp3 = 2*x-3*y;
  const Double_t tmp4 = 18-32*x+12*x*x+48*y-36*x*y+27*y*y;

  return (1+tmp1*tmp1*tmp2)*(30+tmp3*tmp3*tmp4);
}

////////////////////////////////////////////////////////////////////////////////
///
/// F(x,y) = (1 + (x+y+1)^2 * (19-14x+3x^2-14y+6xy+3y^2))
///           * (30 + (2x-3y)^2 * (18-32x+12x^2+48y-36xy+27y^2))
///
///   start point     : F(-0.4,-0,6) = 35
///   local  minima   : F(1.2,0.8)   = 840
///                     F(1.8,0.2)   = 84
///                     F(-0.6,-0.4) = 30
///   global minimum  : F(0.0,-1.0)  = 3
///
/// This is an eighth-order polynomial in two variables which is well behaved near each
/// minimum, but has four local minima and is of course non-positive-definite in many
/// regions. The saddle point between the two lowest minima occurs at F(-0.4,-0.6)=35
/// making this an interesting start point .
/// [Reference: Math. Comp. 25, 571 (1971).]

Bool_t RunGoldStein1()
{
  Bool_t ok = kTRUE;
  const int nvars = 2;

  ROOT::Math::Minimizer* min =
    ROOT::Math::Factory::CreateMinimizer(ROOT::Math::MinimizerOptions::DefaultMinimizerType(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo());
  min->SetPrintLevel(gVerbose);

  ROOT::Math::Functor f(&GoldStein1, nvars);
  min->SetFunction(f);
  double step = 0.01;
  double xmin = -2.0;
  double xmax =  2.0;
  double ymin = 3.0;
  double tolerance = 1.e-3;

  min->SetLimitedVariable(0, "x", -0.3999, step, xmin, xmax);
  min->SetLimitedVariable(1, "y", -0.6,    step, xmin, xmax);

  min->SetMaxFunctionCalls(100000);
  min->SetTolerance(tolerance * gToleranceMult);

  min->Minimize();

  if (min->MinValue() > ymin + tolerance)
    ok = kFALSE;

  delete min;

  return ok;
}

////////////////////////////////////////////////////////////////////////////////

Double_t GoldStein2(const Double_t *par)
{
  const Double_t x = par[0];
  const Double_t y = par[1];

  const Double_t tmp1 = x*x+y*y-25;
  const Double_t tmp2 = TMath::Sin(4*x-3*y);
  const Double_t tmp3 = 2*x+y-10;

  return TMath::Exp(0.5*tmp1*tmp1)+tmp2*tmp2*tmp2*tmp2+0.5*tmp3*tmp3;
}

////////////////////////////////////////////////////////////////////////////////
///
/// F(x,y) = (1 + (x+y+1)^2 * (19-14x+3x^2-14y+6xy+3y^2))
///           * (30 + (2x-3y)^2 * (18-32x+12x^2+48y-36xy+27y^2))
///
///   start point     : F(1.6,3.4) =
///   global minimum  : F(3,4)     = 1
///
/// This function has many local minima .
/// [Reference: Math. Comp. 25, 571 (1971).]

Bool_t RunGoldStein2()
{
  Bool_t ok = kTRUE;
  const int nvars = 2;

  ROOT::Math::Minimizer* min =
    ROOT::Math::Factory::CreateMinimizer(ROOT::Math::MinimizerOptions::DefaultMinimizerType(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo());
  min->SetPrintLevel(gVerbose);

  ROOT::Math::Functor f(&GoldStein2, nvars);
  min->SetFunction(f);
  double step = 0.01;
  double xmin = -5.0;
  double xmax =  5.0;
  double ymin = 1.0;
  double tolerance = 1.e-2;

  min->SetLimitedVariable(0, "x", +1.0, step, xmin, xmax);
  min->SetLimitedVariable(1, "y", +3.2, step, xmin, xmax);

  min->SetMaxFunctionCalls(100000);
  min->SetTolerance(tolerance * gToleranceMult);

  min->Minimize();

  if (min->MinValue() > ymin + tolerance)
    ok = kFALSE;

  delete min;

  return ok;
}

Double_t seed = 3;
Int_t  nf;
TMatrixD A;
TMatrixD B;
TVectorD xx0;
TVectorD sx0;
TVectorD cx0;
TVectorD sx;
TVectorD cx;
TVectorD vv0;
TVectorD vv;
TVectorD rr;

////////////////////////////////////////////////////////////////////////////////

Double_t TrigoFletcher(const Double_t *par)
{
  Int_t i;
  for (i = 0; i < nf ; i++) {
    cx0[i] = TMath::Cos(xx0[i]);
    sx0[i] = TMath::Sin(xx0[i]);
    cx [i] = TMath::Cos(par[i]);
    sx [i] = TMath::Sin(par[i]);
  }

  vv0 = A*sx0+B*cx0;
  vv  = A*sx +B*cx;
  rr  = vv0-vv;

  return rr * rr;
}

////////////////////////////////////////////////////////////////////////////////
///
/// F(\vec{x}) = \sum_{i=1}^n ( E_i - \sum_{j=1}^n (A_{ij} \sin x_j + B_{ij} \cos x_j) )^2
///
///   where E_i = \sum_{j=1}^n ( A_{ij} \sin x_{0j} + B_{ij} \cos x_{0j} )
///
///   B_{ij} and A_{ij} are random matrices composed of integers between -100 and 100;
///   for j = 1,...,n: x_{0j} are any random numbers, -\pi < x_{0j} < \pi;
///
///   start point : x_j = x_{0j} + 0.1 \delta_j,  -\pi < \delta_j < \pi
///   minimum     : F(\vec{x} = \vec{x}_0) = 0
///
/// This is a set of functions of any number of variables n, where the minimum is always
/// known in advance, but where the problem can be changed by choosing different
/// (random) values of the constants A_{ij}, B_{ij}, and x_{0j} . The difficulty can be
/// varied by choosing larger starting deviations \delta_j . In practice, most methods
/// find the "right" minimum, corresponding to \vec{x} = \vec{x}_0, but there are usually
/// many subsidiary minima.
/// [Reference: Comput. J. 6 163 (1963).]

Bool_t RunTrigoFletcher()
{

  const Double_t pi = TMath::Pi();
  Bool_t ok = kTRUE;
  Double_t delta = 0.1;
  Double_t tolerance = 1.e-2;

  for (nf = 5; nf<32;nf +=5) {
     ROOT::Math::Minimizer* min =
       ROOT::Math::Factory::CreateMinimizer(ROOT::Math::MinimizerOptions::DefaultMinimizerType(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo());
     min->SetPrintLevel(gVerbose);

     ROOT::Math::Functor f(&TrigoFletcher, nf);
     min->SetFunction(f);

     A.ResizeTo(nf,nf);
     B.ResizeTo(nf,nf);
     xx0.ResizeTo(nf);
     sx0.ResizeTo(nf);
     cx0.ResizeTo(nf);
     sx.ResizeTo(nf);
     cx.ResizeTo(nf);
     vv0.ResizeTo(nf);
     vv.ResizeTo(nf);
     rr.ResizeTo(nf);
     A.Randomize(-100.,100,seed);
     B.Randomize(-100.,100,seed);
     for (Int_t i = 0; i < nf; i++) {
       for (Int_t j = 0; j < nf; j++) {
         A(i,j) = Int_t(A(i,j));
         B(i,j) = Int_t(B(i,j));
       }
     }

     xx0.Randomize(-pi,pi,seed);
     TVectorD x1(nf); x1.Randomize(-delta*pi,delta*pi,seed);
     x1+= xx0;

     for (Int_t i = 0; i < nf; i++)
       min->SetLimitedVariable(i, Form("x_%d",i), x1[i], 0.01, -pi*(1+delta), +pi*(1+delta));

     min->SetMaxFunctionCalls(100000);
     min->SetTolerance(tolerance * gToleranceMult);

     min->Minimize();

     if (min->MinValue() > tolerance)
       ok = kFALSE;

     delete min;
  }

  return ok;
}

////////////////////////////////////////////////////////////////////////////////

Int_t stressFit(const char *type, const char *algo, Int_t N)
{
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer(type, algo);

  std::cout << "******************************************************************" <<std::endl;
  std::cout << "*  Minimization - S T R E S S suite                              *" <<std::endl;
  std::cout << "******************************************************************" <<std::endl;
  std::cout << "******************************************************************" <<std::endl;

   TStopwatch timer;
   timer.Start();

  std::cout << "*  Starting  S T R E S S  with fitter : "
            << ROOT::Math::MinimizerOptions::DefaultMinimizerType() << " / "
            << ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo() << std::endl;
  std::cout << "******************************************************************" << std::endl;

  gBenchmark->Start("stressFit");

  int okRosenBrock    = 0;
  int okWood          = 0;
  int okPowell        = 0;
  int okFletcher      = 0;
  int okGoldStein1    = 0;
  int okGoldStein2    = 0;
  int okTrigoFletcher = 0;
  Int_t i;
  for (i = 0; i < N; i++) if (RunWood4()) okWood++;
  StatusPrint(1, "Wood", okWood, N);
  for (i = 0; i < N; i++) if (RunRosenBrock()) okRosenBrock++;
  StatusPrint(2, "RosenBrock", okRosenBrock, N);
  for (i = 0; i < N; i++) if (RunPowell()) okPowell++;
  StatusPrint(3, "Powell", okPowell, N);
  for (i = 0; i < N; i++) if (RunFletcher()) okFletcher++;
  StatusPrint(4, "Fletcher", okFletcher, N);
  for (i = 0; i < N; i++) if (RunGoldStein1()) okGoldStein1++;
  StatusPrint(5, "GoldStein1", okGoldStein1, N);
  for (i = 0; i < N; i++) if (RunGoldStein2()) okGoldStein2++;
  StatusPrint(6, "GoldStein2", okGoldStein2, N);
  if (RunTrigoFletcher()) okTrigoFletcher++;
  StatusPrint(7, "TrigoFletcher", okTrigoFletcher, 1);

  gBenchmark->Stop("stressFit");


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
    if (!os) printf("*  SYS: Windows 95\n");
    else     printf("*  SYS: %s %s \n",os,gSystem->Getenv("PROCESSOR_IDENTIFIER"));
  }

  printf("******************************************************************\n");
  gBenchmark->Print("stressFit");
#ifdef __CINT__
  Double_t reftime = 86.34; //macbrun interpreted
#else
  Double_t reftime = 12.07; //macbrun compiled
#endif
  const Double_t rootmarks = 800.*reftime/gBenchmark->GetCpuTime("stressFit");

  printf("******************************************************************\n");
  printf("*  ROOTMARKS =%6.1f   *  Root%-8s  %d/%d\n",rootmarks,gROOT->GetVersion(),
         gROOT->GetVersionDate(),gROOT->GetVersionTime());
  printf("******************************************************************\n");

   return 0;
}

//_____________________________batch only_____________________
#ifndef __CINT__

int main(int argc,const char *argv[])
{
  gBenchmark = new TBenchmark();
  const char *fittertype = "Minuit";
  const char *fitteralgo = "Migrad";
  if (argc > 1)  fittertype = argv[1];
  if (argc > 2)  fitteralgo = argv[2];
  if (strcmp(fittertype, "Minuit") && strcmp(fittertype, "Minuit2") && strcmp(fittertype, "Fumili")) {
     printf("stressFit illegal option %s, using Minuit instead\n", fittertype);
     fittertype = "Minuit";
     fitteralgo = "Migrad";
  }
  Int_t N = 2000;
  if (argc > 3) N = atoi(argv[3]);
  stressFit(fittertype, fitteralgo, N);  //default is Minuit
  return 0;
}

#endif
