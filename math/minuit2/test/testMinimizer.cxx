// test of minimization using new minimizer classes

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include "TVirtualFitter.h"

#include "Math/IFunction.h"
#include "Math/Util.h"
#include <cmath>
#include <cassert>

#include <string>
#include <iostream>

#include "TStopwatch.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TRandom3.h"
#include "TMath.h"

int gNCall = 0;
int gNCall2 = 0;
int gNmin = 1000;
int gVerbose = 0;
bool useGradient = false;

bool minos = true;

double gAbsTolerance = 0.005;

// Rosenbrok function to be minimize

typedef void (*FCN)(Int_t &npar, Double_t *gin, Double_t &f, Double_t *u, Int_t flag);

// ROSENBROCK function
////////////////////////////////////////////////////////////////////////////////

void RosenBrock(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t /*iflag*/)
{
   gNCall++;
   const Double_t x = par[0];
   const Double_t y = par[1];
   const Double_t tmp1 = y - x * x;
   const Double_t tmp2 = 1 - x;
   f = 100 * tmp1 * tmp1 + tmp2 * tmp2;
}

class RosenBrockFunction : public ROOT::Math::IMultiGenFunction {

public:
   virtual ~RosenBrockFunction() {}

   unsigned int NDim() const { return 2; }

   ROOT::Math::IMultiGenFunction *Clone() const { return new RosenBrockFunction(); }

   const double *TrueMinimum() const
   {
      fTrueMin[0] = 1;
      fTrueMin[1] = 1;
      return fTrueMin;
   }

private:
   inline double DoEval(const double *x) const
   {
#ifdef USE_FREE_FUNC
      double f = 0;
      int ierr = 0;
      int i = 0;
      RosenBrock(i, 0, f, const_cast<double *>(x), ierr);
      return f;
#else
      gNCall++;
      const Double_t xx = x[0];
      const Double_t yy = x[1];
      const Double_t tmp1 = yy - xx * xx;
      const Double_t tmp2 = 1 - xx;
      return 100 * tmp1 * tmp1 + tmp2 * tmp2;
#endif
   }

   mutable double fTrueMin[2];
};

// TRIGONOMETRIC FLETCHER FUNCTION

class TrigoFletcherFunction : public ROOT::Math::IMultiGradFunction {

public:
   TrigoFletcherFunction(unsigned int dim) : fDim(dim)
   {
      double seed = 3;
      A.ResizeTo(dim, dim);
      B.ResizeTo(dim, dim);
      x0.ResizeTo(dim);
      sx0.ResizeTo(dim);
      cx0.ResizeTo(dim);
      sx.ResizeTo(dim);
      cx.ResizeTo(dim);
      v0.ResizeTo(dim);
      v.ResizeTo(dim);
      r.ResizeTo(dim);
      A.Randomize(-100., 100, seed);
      B.Randomize(-100., 100, seed);
      for (unsigned int i = 0; i < dim; i++) {
         for (unsigned int j = 0; j < dim; j++) {
            A(i, j) = int(A(i, j));
            B(i, j) = int(B(i, j));
         }
      }
      x0.Randomize(-TMath::Pi(), TMath::Pi(), seed);
      // calculate vector Ei
      for (unsigned int i = 0; i < fDim; ++i) {
         cx0[i] = std::cos(x0[i]);
         sx0[i] = std::sin(x0[i]);
      }
      v0 = A * sx0 + B * cx0;
   }

   unsigned int NDim() const { return fDim; }

   ROOT::Math::IMultiGenFunction *Clone() const
   {
      TrigoFletcherFunction *f = new TrigoFletcherFunction(*this);
      //       std::cerr <<"cannot clone this function" << std::endl;
      //       assert(0);
      return f;
   }

   void StartPoints(double *x, double *s)
   {
      TRandom3 rndm;
      const double stepSize = 0.01;
      const double deltaAmp = 0.1;
      const double pi = TMath::Pi();
      for (unsigned int i = 0; i < fDim; ++i) {
         double delta = rndm.Uniform(-deltaAmp * pi, deltaAmp * pi);
         x[i] = x0(i) + 0.1 * delta;
         if (x[i] <= -pi)
            x[i] += 2. * pi;
         if (x[i] > pi)
            x[i] -= 2. * pi;
         s[i] = stepSize;
      }
   }

   const double *TrueMinimum() const { return x0.GetMatrixArray(); }

   void Gradient(const double *x, double *g) const
   {
      gNCall2++;

      for (unsigned int i = 0; i < fDim; ++i) {
         cx[i] = std::cos(x[i]);
         sx[i] = std::sin(x[i]);
      }

      v = A * sx + B * cx;
      r = v0 - v;

      // calculate the grad components
      for (unsigned int i = 0; i < fDim; ++i) {
         g[i] = 0;
         for (unsigned int k = 0; k < fDim; ++k) {
            g[i] += 2. * r(k) * (-A(k, i) * cx(i) + B(k, i) * sx(i));
         }
      }
   }

#ifdef USE_FDF
   void FdF(const double *x, double &f, double *g) const
   {
      gNCall++;

      for (unsigned int i = 0; i < fDim; ++i) {
         cx[i] = std::cos(x[i]);
         sx[i] = std::sin(x[i]);
      }

      v = A * sx + B * cx;
      r = v0 - v;

      f = r * r;

      // calculate the grad components
      for (unsigned int i = 0; i < fDim; ++i) {
         g[i] = 0;
         for (unsigned int k = 0; k < fDim; ++k) {
            g[i] += 2. * r(k) * (-A(k, i) * cx(i) + B(k, i) * sx(i));
         }
      }
   }
#endif

private:
   //    TrigoFletcherFunction(const TrigoFletcherFunction & ) {}
   //    TrigoFletcherFunction & operator=(const TrigoFletcherFunction &) { return *this; }

   double DoEval(const double *x) const
   {
      gNCall++;

      for (unsigned int i = 0; i < fDim; ++i) {
         cx[i] = std::cos(x[i]);
         sx[i] = std::sin(x[i]);
      }

      v = A * sx + B * cx;
      r = v0 - v;

      return r * r;
   }

   double DoDerivative(const double *x, unsigned int i) const
   {
      std::vector<double> g(fDim);
      Gradient(x, &g[0]);
      return g[i];
   }

private:
   unsigned int fDim;

   TMatrixD A;
   TMatrixD B;
   TVectorD x0;
   mutable TVectorD sx0;
   mutable TVectorD cx0;
   mutable TVectorD sx;
   mutable TVectorD cx;
   mutable TVectorD v0;
   mutable TVectorD v;
   mutable TVectorD r;
};

// CHEBYQUAD FUNCTION

class ChebyQuadFunction : public ROOT::Math::IMultiGradFunction {

public:
   ChebyQuadFunction(unsigned int n) : fDim(n), fvec(std::vector<double>(n)), fTrueMin(std::vector<double>(n)) {}

   unsigned int NDim() const { return fDim; }

   ROOT::Math::IMultiGenFunction *Clone() const { return new ChebyQuadFunction(*this); }

   const double *TrueMinimum() const { return &fTrueMin[0]; }

   // use equally spaced points
   void StartPoints(double *x, double *s)
   {
      for (unsigned int i = 0; i < fDim; ++i) {
         s[i] = 0.01;
         x[i] = double(i) / (double(fDim) + 1.0);
      }
   }

   // compute gradient

   void Gradient(const double *x, double *g) const
   {
      gNCall2++;
      unsigned int n = fDim;
      // estimate first the fvec
      DoCalculatefi(x);

      for (unsigned int j = 0; j < n; ++j) {
         g[j] = 0.0;
         double t1 = 1.0;
         double t2 = 2.0 * x[j] - 1.0;
         double t = 2.0 * t2;
         double s1 = 0.0;
         double s2 = 2.0;
         for (unsigned int i = 0; i < n; ++i) {
            g[j] += fvec[i] * s2;
            double th = 4.0 * t2 + t * s2 - s1;
            s1 = s2;
            s2 = th;
            th = t * t2 - t1;
            t1 = t2;
            t2 = th;
         }
         g[j] = 2. * g[j] / double(n);
      }
   }

private:
   double DoEval(const double *x) const
   {

      gNCall++;
      DoCalculatefi(x);
      double f = 0;
      for (unsigned int i = 0; i < fDim; ++i)
         f += fvec[i] * fvec[i];

      return f;
   }

   double DoDerivative(const double *x, unsigned int i) const
   {
      std::vector<double> g(fDim);
      Gradient(x, &g[0]);
      return g[i];
   }

   void DoCalculatefi(const double *x) const
   {
      // calculate the i- element ; F(X) = Sum {fi]
      unsigned int n = fDim;
      for (unsigned int i = 0; i < n; ++i)
         fvec[i] = 0;

      for (unsigned int j = 0; j < n; ++j) {
         double t1 = 1.0;
         double t2 = 2.0 * x[j] - 1.0;
         double t = 2.0 * t2;
         for (unsigned int i = 0; i < n; ++i) {
            fvec[i] += t2;
            double th = t * t2 - t1;
            t1 = t2;
            t2 = th;
         }
      }

      // sum with the integral (integral is zero for odd Cheb polynomial and = 1/(i**2 -1) for the even ones
      for (unsigned int i = 1; i <= n; ++i) {
         int l = i - 1;
         fvec[l] /= double(n);
         if ((i % 2) == 0) {
            fvec[l] += 1.0 / (double(i * i) - 1.0);
         }
      }
   }

   unsigned int fDim;
   mutable std::vector<double> fvec;
   mutable std::vector<double> fTrueMin;
};

// WOOD function (4 dim function)

double WoodFunction(const double *par)
{
   gNCall++;

   const double w = par[0];
   const double x = par[1];
   const double y = par[2];
   const double z = par[3];

   const double w1 = w - 1;
   const double x1 = x - 1;
   const double y1 = y - 1;
   const double z1 = z - 1;
   const double tmp1 = x - w * w;
   const double tmp2 = z - y * y;

   double f = 100 * tmp1 * tmp1 + w1 * w1 + 90 * tmp2 * tmp2 + y1 * y1 + 10.1 * (x1 * x1 + z1 * z1) + 19.8 * x1 * z1;
   return f;
}

// Powell Function (4 dim function)
double PowellFunction(const double *par)
{
   gNCall++;

   const double w = par[0];
   const double x = par[1];
   const double y = par[2];
   const double z = par[3];

   const double tmp1 = w + 10 * x;
   const double tmp2 = y - z;
   const double tmp3 = x - 2 * y;
   const double tmp4 = w - z;

   double f = tmp1 * tmp1 + 5 * tmp2 * tmp2 + tmp3 * tmp3 * tmp3 * tmp3 + 10 * tmp4 * tmp4 * tmp4 * tmp4;
   return f;
}

double SimpleQuadFunction(const double *par)
{
   gNCall++;
   double x = par[0];
   double y = par[1];
   double f = x * x + 2 * y * y - x * y;
   std::cout << "Quadfunc " << gNCall << "\t" << x << " , " << y << " f = " << f << std::endl;
   return f;
}

const double *TrueMinimum(const ROOT::Math::IMultiGenFunction &func)
{

   const RosenBrockFunction *fRB = dynamic_cast<const RosenBrockFunction *>(&func);
   if (fRB != 0)
      return fRB->TrueMinimum();
   const TrigoFletcherFunction *fTF = dynamic_cast<const TrigoFletcherFunction *>(&func);
   if (fTF != 0)
      return fTF->TrueMinimum();
   //    const ChebyQuadFunction * fCQ = dynamic_cast< const ChebyQuadFunction *> (&func);
   //    if (fCQ != 0) return fCQ->TrueMinimum();
   return 0;
}

void printMinimum(const std::vector<double> &x)
{
   std::cout << "Minimum X values\n";
   std::cout << "\t";
   int pr = std::cout.precision(12);
   unsigned int n = x.size();
   for (unsigned int i = 0; i < n; ++i) {
      std::cout << x[i];
      if (i != n - 1)
         std::cout << " , ";
      if (i > 0 && i % 6 == 0)
         std::cout << "\n\t";
   }
   std::cout << std::endl;
   std::cout.precision(pr);
}

int DoNewMinimization(const ROOT::Math::IMultiGenFunction &func, const double *x0, const double *s0,
                      ROOT::Math::Minimizer *min, double &minval, double &edm, double *minx)
{

   int iret = 0;

   min->SetMaxFunctionCalls(1000000);
   min->SetMaxIterations(100000);
   min->SetTolerance(gAbsTolerance);
   if (func.NDim() >= 10) {
      min->SetTolerance(0.01);
   }

   min->SetPrintLevel(gVerbose);
   // check if func provides gradient
   const ROOT::Math::IMultiGradFunction *gfunc = dynamic_cast<const ROOT::Math::IMultiGradFunction *>(&func);
   if (gfunc != 0 && useGradient)
      min->SetFunction(*gfunc);
   else
      min->SetFunction(func);

   for (unsigned int i = 0; i < func.NDim(); ++i) {
#ifdef DEBUG
      std::cout << "set variable " << i << " to value " << x0[i] << std::endl;
#endif
      min->SetVariable(i, "x" + ROOT::Math::Util::ToString(i), x0[i], s0[i]);
   }

   bool ret = min->Minimize();
   minval = min->MinValue();
   edm = min->Edm();

   if (!ret) {
      return -1;
   }

   const double *xmin = min->X();

   bool ok = true;
   const double *trueMin = TrueMinimum(func);
   if (trueMin != 0) {
      for (unsigned int i = 0; i < func.NDim(); ++i)
         ok &= (std::fabs(xmin[i] - trueMin[i]) < gAbsTolerance);
   }

   if (!ok)
      iret = -2;
   int ncall_migrad = gNCall;

   // test Minos (use the default up of 1)
   if (minos) {

      double el, eu;
      for (unsigned int i = 0; i < func.NDim(); ++i) {
         ret = min->GetMinosError(i, el, eu);
         std::cout << "ncalls " << gNCall << "\t";
         if (ret)
            std::cout << "MINOS error  for " << i << " = " << el << "   " << eu << std::endl;
         else
            std::cout << "MINOS failed for " << i << std::endl;
      }
   }

   std::cout << "\nMigrad Ncalls:\t " << ncall_migrad << std::endl;
   std::cout << "Minos Ncalls :\t " << gNCall - ncall_migrad << std::endl;

   //   std::cout << "function at the minimum " << func(xmin) << std::endl;
   std::copy(xmin, xmin + func.NDim(), minx);
   min->Clear();

   return iret;
}

int DoOldMinimization(FCN func, TVirtualFitter *min, double &minval, double &edm)
{

   int iret = 0;

   assert(min != 0);
   min->SetFCN(func);

   Double_t arglist[100];
   arglist[0] = gVerbose - 1;
   min->ExecuteCommand("SET PRINT", arglist, 1);

   min->SetParameter(0, "x", -1.2, 0.01, 0, 0);
   min->SetParameter(1, "y", 1.0, 0.01, 0, 0);

   arglist[0] = 0;
   min->ExecuteCommand("SET NOW", arglist, 0);
   arglist[0] = 1000;          // number of function calls
   arglist[1] = gAbsTolerance; // tolerance
   // min->ExecuteCommand("MIGRAD",arglist,0);
   min->ExecuteCommand("MIGRAD", arglist, 2);

   if (minos)
      min->ExecuteCommand("MINOS", arglist, 0);

   Double_t parx, pary;
   Double_t we, al, bl;
   Char_t parName[32];
   min->GetParameter(0, parName, parx, we, al, bl);
   min->GetParameter(1, parName, pary, we, al, bl);

   bool ok = (TMath::Abs(parx - 1.) < gAbsTolerance && TMath::Abs(pary - 1.) < gAbsTolerance);

   double errdef = 0;
   int nvpar, nparx;
   min->GetStats(minval, edm, errdef, nvpar, nparx);
   if (!ok)
      iret = -2;

   min->Clear(); // for further use
   return iret;
}

int testNewMinimizer(const ROOT::Math::IMultiGenFunction &func, const double *x0, const double *s0,
                     const std::string &minimizer, const std::string &algoType)
{

   std::cout << "\n************************************************************\n";
   std::cout << "\tTest new ROOT::Math::Minimizer\n";
   std::cout << "\tMinimizer is " << minimizer << "  " << algoType << std::endl;

   int iret = 0;
   double minval = 0., edm = 0.;
   std::vector<double> xmin(func.NDim());

   TStopwatch w;
   w.Start();

   ROOT::Math::Minimizer *min = ROOT::Math::Factory::CreateMinimizer(minimizer, algoType);
   if (min == 0) {
      std::cout << "Error using minimizer " << minimizer << "  " << algoType << std::endl;
      return -1;
   }

   for (int i = 0; i < gNmin; ++i) {
      gNCall = 0;
      gNCall2 = 0;
      iret |= DoNewMinimization(func, x0, s0, min, minval, edm, &xmin[0]);
   }

   w.Stop();
   if (iret != 0)
      std::cout << "\n****** ERROR:   Minimization FAILED ! \n";
   int pr = std::cout.precision(18);
   std::cout << "\nNCalls: \t" << gNCall << " , " << gNCall2 << "\tMinValue: \t" << minval << "\tEdm: \t" << edm;
   std::cout.precision(pr);
   std::cout << "\nTime:   \t" << w.RealTime() << " , " << w.CpuTime() << std::endl;
   printMinimum(xmin);
   std::cout << "\n************************************************************\n";

#ifdef CHECK_WITHMINUIT
   // do Minuit after BFGS
   if (minimizer == "GSL_BFGS") {
      std::cout << "DO Minuit2 from last point\n";
      gNCall = 0;
      iret |= DoNewMinimization(func, &xmin.front(), s0, "Minuit2", "", minval, edm, &xmin[0]);
      int pr = std::cout.precision(18);
      std::cout << "\nNCalls: \t" << gNCall << "\tMinValue: \t" << minval << "\tEdm: \t" << edm;
      std::cout.precision(pr);
      std::cout << std::endl;
   }
#endif

   delete min;

   return iret;
}

int testOldMinimizer(FCN func, const std::string &fitter, int n = 25)
{

   std::cout << "\n************************************************************\n";
   std::cout << "\tTest using TVirtualFitter\n";
   std::cout << "\tFitter is " << fitter << std::endl;

   int iret = 0;
   double minval = 0., edm = 0.;

   TStopwatch w;
   w.Start();

   TVirtualFitter::SetDefaultFitter(fitter.c_str());

   TVirtualFitter *min = TVirtualFitter::Fitter(0, n);

   // min->Dump();

   for (int i = 0; i < gNmin; ++i) {
      gNCall = 0;
      iret |= DoOldMinimization(func, min, minval, edm);
   }

   w.Stop();
   if (iret != 0)
      std::cout << "\n****** ERROR:   Minimization FAILED ! \n";
   int pr = std::cout.precision(18);
   std::cout << "\nNCalls: \t" << gNCall << "\tMinValue: \t" << minval << "\tEdm: \t" << edm;
   std::cout.precision(pr);
   std::cout << "\nTime: \t" << w.RealTime() << " , " << w.CpuTime() << std::endl;
   std::cout << "\n************************************************************\n";

   return iret;
}

int testRosenBrock()
{

   int iret = 0;

   std::cout << "\n************************************************************\n";
   std::cout << "\tROSENBROCK function test\n\n";

   double s0[2] = {0.01, 0.01};

   // minimize using Rosenbrock Function

   gNmin = 1;

   RosenBrockFunction fRB;
   double xRB[2] = {-1., 1.2};

   iret |= testNewMinimizer(fRB, xRB, s0, "Minuit", "");
   iret |= testNewMinimizer(fRB, xRB, s0, "Minuit2", "");

   //    iret |= testNewMinimizer(fRB,xRB,s0,"GSLMultiMin","ConjugateFR");
   //    iret |= testNewMinimizer(fRB,xRB,s0,"GSLMultiMin","ConjugatePR");
   //    iret |= testNewMinimizer(fRB,xRB,s0,"GSLMultiMin","BFGS");
   //    iret |= testNewMinimizer(fRB,xRB,s0,"GSLMultiMin","BFGS2");

   //    iret |= testOldMinimizer(RosenBrock,"Minuit",2);
   //    iret |= testOldMinimizer(RosenBrock,"Minuit2",2);

   return iret;
}

int testChebyQuad()
{

   int iret = 0;

   // test with ChebyQuad function

   const int n = 8;
   ChebyQuadFunction func(n);

   gNmin = 1;

   double s0[n];
   double x0[n];
   func.StartPoints(x0, s0);

   std::cout << "\n************************************************************\n";
   std::cout << "\tCHEBYQUAD function test , n = " << n << "\n\n";

   //    double x[8] = {0.043153E+00, 0.193091E+00, 0.266329E+00, 0.500000E+00,
   //                   0.500000E+00, 0.733671E+00, 0.806910E+00, 0.956847E+00 };
   //    double x[2] = {0.5, 0.5001};
   //    std::cout << "FUNC " << func(x) << std::endl;
   double x1[100] = {
      0.00712780070646, 0.0123441993113, 0.0195428378255, 0.0283679084192, 0.0385291131289, 0.0492202424892,
      0.0591277976178,  0.0689433195252, 0.0791293590525, 0.088794974369,  0.0988949579193, 0.108607151294,
      0.118571075831,   0.128605446238,  0.137918291068,  0.149177761352,  0.156665324587,  0.170851061982,
      0.174688134016,   0.192838903364,  0.193078270803,  0.209255377225,  0.217740096779,  0.225888518345,
      0.241031047421,   0.244253844041,  0.257830449676,  0.269467652526,  0.274286498012,  0.288877029988,
      0.297549406597,   0.304950954529,  0.319230811642,  0.326387092206,  0.335229058731,  0.349178359226,
      0.355905988048,   0.365197862755,  0.379068092603,  0.385826036925,  0.394978252826,  0.408974425717,
      0.415968185065,   0.424621041584,  0.438837361714,  0.446214149031,  0.454242324351,  0.468614308013,
      0.476506553416,   0.483916944941,  0.498229247409,  0.506794629616,  0.513736742474,  0.527712475478,
      0.537073277673,   0.543731917673,  0.557187513963,  0.567346279639,  0.57379846397,   0.586691058785,
      0.597561941009,   0.60382873461,   0.616316037506,  0.627719652101,  0.633760038662,  0.646175283836,
      0.657809344891,   0.663569004722,  0.676314563639,  0.687674566849,  0.69332205923,   0.706839545953,
      0.716907408637,   0.723407327715,  0.738019389561,  0.744806584048,  0.754657613362,  0.769181875619,
      0.772250323489,   0.787104833193,  0.795856360905,  0.804099304478,  0.82142178741,   0.819589601284,
      0.839024540481,   0.842457233039,  0.857393475964,  0.86408033345,   0.876329840525,  0.884867318008,
      0.895744532071,   0.905113958163,  0.915445338697,  0.925148068352,  0.935344457785,  0.945127838313,
      0.955272197168,   0.965687518559,  0.975129521484,  0.982662007764};

   std::cout << "FUNC " << func(x1) << std::endl;

   iret |= testNewMinimizer(func, x0, s0, "Minuit", "");
   iret |= testNewMinimizer(func, x0, s0, "Minuit2", "");

   //    iret |= testNewMinimizer(func,x0,s0, "GSLMultiMin","ConjugateFR");
   //    iret |= testNewMinimizer(func,x0,s0, "GSLMultiMin","ConjugatePR");
   //    iret |= testNewMinimizer(func,x0,s0, "GSLMultiMin","BFGS");

   return iret;
}

int testTrigoFletcher()
{

   int iret = 0;

   // test with fletcher trigonometric function

   gNmin = 1;

   const int nT = 50;
   TrigoFletcherFunction fTrigo(nT);
   double sTrigo[nT];
   double xTrigo[nT];
   fTrigo.StartPoints(xTrigo, sTrigo);

   std::cout << "\n************************************************************\n";
   std::cout << "\tTRIGONOMETRIC Fletcher function test , n = " << nT << "\n\n";

   iret |= testNewMinimizer(fTrigo, xTrigo, sTrigo, "Minuit2", "");
   iret |= testNewMinimizer(fTrigo, xTrigo, sTrigo, "Minuit", "");

   //    iret |= testNewMinimizer(fTrigo,xTrigo,sTrigo,"GSLMultiMin","ConjugateFR");
   //    iret |= testNewMinimizer(fTrigo,xTrigo,sTrigo,"GSLMultiMin","ConjugatePR");
   //    iret |= testNewMinimizer(fTrigo,xTrigo,sTrigo,"GSLMultiMin","BFGS");

   return iret;
}

int testWood()
{

   int iret = 0;

   // test with Wood function (4d)
   //   minimum    : F(1,1,1,1)  =   0.

   gNmin = 1;

   ROOT::Math::Functor f(&WoodFunction, 4);

   double x0[4] = {-3, -1, -3, -1};
   double s0[4] = {0.1, 0.1, 0.1, 0.1};

   std::cout << "\n************************************************************\n";
   std::cout << "\tWOOD 4 function test  \n\n";

   iret |= testNewMinimizer(f, x0, s0, "Minuit", "");
   iret |= testNewMinimizer(f, x0, s0, "Minuit2", "");

   return iret;
}

int testPowell()
{

   int iret = 0;

   // test with Powell function (4d)
   // minimum is at F(0,0,0,0) = 0

   gNmin = 1;

   ROOT::Math::Functor f(&PowellFunction, 4);

   double x0[4] = {-3, -1, 0, 1};
   double s0[4] = {0.1, 0.1, 0.1, 0.1};

   std::cout << "\n************************************************************\n";
   std::cout << "\tPOWELL function test   \n\n";

   iret |= testNewMinimizer(f, x0, s0, "Minuit", "");
   iret |= testNewMinimizer(f, x0, s0, "Minuit2", "");

   return iret;
}

int testQuadFunc()
{

   int iret = 0;

   // test with a simple quadratic function 2d
   // minimum is at F(0,0) = 0

   gNmin = 1;

   ROOT::Math::Functor f(&SimpleQuadFunction, 2);

   double x0[4] = {-3, -3};
   double s0[4] = {0.1, 0.1};

   std::cout << "\n************************************************************\n";
   std::cout << "\tSIMPLE QUAD  function test   \n\n";

   iret |= testNewMinimizer(f, x0, s0, "Minuit", "");
   iret |= testNewMinimizer(f, x0, s0, "Minuit2", "");

   return iret;
}

int main()
{

   int iret = 0;

#ifdef DEBUG
   gVerbose = 4;
   gNmin = 1;
#endif

   iret |= testRosenBrock();
   iret |= testChebyQuad();
   iret |= testTrigoFletcher();

   iret |= testWood();
   iret |= testPowell();

   iret |= testQuadFunc();

   if (iret != 0)
      std::cerr << "testMinim :\t FAILED " << std::endl;
   else
      std::cerr << "testMinim :\t OK " << std::endl;
   return iret;
}
