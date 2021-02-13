// @(#)root/test:$Id$
// Author: Lorenzo Moneta   06/2005
///////////////////////////////////////////////////////////////////////////////////
//
//  MathMore Benchmark test suite
//  ==============================
//
//  This program performs tests :
//     - numerical integration, derivation and root finders
//    - it compares for various values of the gamma and beta distribution)
//          - the numerical calculated integral of pdf  with cdf function,
//           - the calculated derivative of cdf with pdf
//           - the inverse (using root finder) of cdf with quantile
//
//     to run the program outside ROOT do:
//        > make stressMathMore
//        > ./stressMathMore
//
//     to run the program in ROOT
//       root> gSystem->Load("libMathMore")
//       root> .x stressMathMore.cxx+
//

#include "Math/DistFunc.h"
#include "Math/IParamFunction.h"
#include "Math/Integrator.h"
#include "Math/Derivator.h"
#include "Math/Functor.h"
#include "Math/RootFinderAlgorithms.h"
#include "Math/RootFinder.h"

#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <cstdlib>
#include "TBenchmark.h"
#include "TROOT.h"
#include "TRandom3.h"

using namespace ROOT::Math;




#ifdef __CINT__
#define INF 1.7E308
#else
#define INF std::numeric_limits<double>::infinity()
#endif

//#define DEBUG

bool debug = true;  // print out reason of test failures
bool removeFiles = false; // remove Output root files

void PrintTest(std::string name) {
   std::cout << std::left << std::setw(40) << name;
}

void PrintStatus(int iret) {
   if (iret == 0)
      std::cout <<"\t\t................ OK" << std::endl;
   else
      std::cout <<"\t\t............ FAILED " << std::endl;
}


int compare( std::string name, double v1, double v2, double scale = 2.0) {
  //  ntest = ntest + 1;

   //std::cout << std::setw(50) << std::left << name << ":\t";

  // numerical double limit for epsilon
   double eps = scale* std::numeric_limits<double>::epsilon();
   int iret = 0;
   double delta = v2 - v1;
   double d = 0;
   if (delta < 0 ) delta = - delta;
   if (v1 == 0 || v2 == 0) {
      if  (delta > eps ) {
         iret = 1;
      }
   }
   // skip case v1 or v2 is infinity
   else {
      d = v1;

      if ( v1 < 0) d = -d;
      // add also case when delta is small by default
      if ( delta/d  > eps && delta > eps )
         iret =  1;
   }

   if (iret) {
      if (debug) {
         int pr = std::cout.precision (18);
         std::cout << "\nDiscrepancy in " << name.c_str() << "() :\n  " << v1 << " != " << v2 << " discr = " << int(delta/d/eps)
                   << "   (Allowed discrepancy is " << eps  << ")\n\n";
         std::cout.precision (pr);
      //nfail = nfail + 1;
      }
   }
   //else
      //  std::cout <<".";

   return iret;
}

// typedef for a free function like gamma(double x, double a, double b)
// (dont have blank spaces between for not confusing CINT parser)
typedef double (*FreeFunc3)(double, double, double );
typedef double (*FreeFunc4)(double, double, double, double );

//implement simple functor
struct Func {
   virtual ~Func() {}
   virtual double operator() (double , double, double) const = 0;

};
struct Func3 : public Func {
   Func3(FreeFunc3 f) : fFunc(f) {};
   double operator() (double x, double a, double b) const {
      return fFunc(x,a,b);
   }
   FreeFunc3 fFunc;
};

struct Func4 : public Func {
   Func4(FreeFunc4 f) : fFunc(f) {};
   double operator() (double x, double a, double b) const {
      return fFunc(x,a,b,0.);
   }
   FreeFunc4 fFunc;
};



// statistical function class
const int NPAR = 2;

class StatFunction : public ROOT::Math::IParamFunction {

public:

   StatFunction(Func & pdf, Func & cdf, Func & quant,
                double x1 = -INF,
                double x2 = INF ) :
      fPdf(&pdf), fCdf(&cdf), fQuant(&quant),
      xmin(0.), xmax(0.),
      xlow(x1), xup(x2),
      fHasLowRange(false), fHasUpRange(false), fStartRoot(0.)
   {
      fScaleIg = 10; //scale for integral test
      fScaleDer = 1;  //scale for der test
      fScaleInv = 100;  //scale for inverse test
      for(int i = 0; i< NPAR; ++i) fParams[i]=0;
      NFuncTest = 100;
      if (xlow > -INF) fHasLowRange = true;
      if (xup < INF) fHasUpRange = true;
   }



   unsigned int NPar() const { return NPAR; }
   const double * Parameters() const { return fParams; }
   ROOT::Math::IGenFunction * Clone() const { return new StatFunction(*fPdf,*fCdf,*fQuant); }

   void SetParameters(const double * p) { std::copy(p,p+NPAR,fParams); }

   void SetParameters(double p0, double p1) { *fParams = p0; *(fParams+1) = p1; }

   void SetTestRange(double x1, double x2) { xmin = x1; xmax = x2; }
   void SetNTest(int n) { NFuncTest = n; }
   void SetStartRoot(double x) { fStartRoot =x; }


   double Pdf(double x) const {
      return (*this)(x);
   }

   double Cdf(double x) const {
      return  (*fCdf) ( x, fParams[0], fParams[1] );
   }

   double Quantile(double x) const {
      return (*fQuant)( x, fParams[0], fParams[1]  );
   }


   // test integral with cdf function
   int TestIntegral(IntegrationOneDim::Type algotype);

   // test derivative from cdf to pdf function
   int TestDerivative();

   // test root finding algorithm for finding inverse of cdf
   int TestInverse1(RootFinder::EType algotype);

   // test root finding algorithm for finding inverse of cdf using drivatives
   int TestInverse2(RootFinder::EType algotype);


   void SetScaleIg(double s) { fScaleIg = s; }
   void SetScaleDer(double s) { fScaleDer = s; }
   void SetScaleInv(double s) { fScaleInv = s; }


private:

   double DoEvalPar(double x, const double * ) const {
      // use esplicity cached param values
      return (*fPdf)(x, *fParams, *(fParams+1));
   }

   Func * fPdf;
   Func *  fCdf;
   Func *  fQuant;
   double fParams[NPAR];
   double fScaleIg;
   double fScaleDer;
   double fScaleInv;
   int NFuncTest;
   double xmin;
   double xmax;
   double xlow;
   double xup;
   bool fHasLowRange;
   bool fHasUpRange;
   double fStartRoot;
};

// test integral of function

int StatFunction::TestIntegral(IntegrationOneDim::Type algoType = IntegrationOneDim::kADAPTIVESINGULAR) {

   int iret = 0;

   // scan all values from xmin to xmax
   double dx = (xmax-xmin)/NFuncTest;

   // create Integrator
   Integrator ig(algoType, 1.E-12,1.E-12,100000);
   ig.SetFunction(*this);

   for (int i = 0; i < NFuncTest; ++i) {
      double v1 = xmin + dx*i;  // value used  for testing
      double q1 = Cdf(v1);
      //std::cout << "v1 " << v1 << " pdf " << (*this)(v1) << " cdf " << q1 << " quantile " << Quantile(q1) << std::endl;
      // calculate integral of pdf
      double q2 = 0;

      // lower integral (cdf)
      if (!fHasLowRange)
         q2 = ig.IntegralLow(v1);
      else
         q2 = ig.Integral(xlow,v1);

      int r  = ig.Status();
      // use a larger scale (integral error is 10-9)
      double err = ig.Error();
      //std::cout << "integral result is = " << q2 << " error is " << err << std::endl;
      // Gauss integral sometimes returns an error of 0
      err = std::max(err,  std::numeric_limits<double>::epsilon() );
      double scale = std::max( fScaleIg * err / std::numeric_limits<double>::epsilon(), 1.);
      r |= compare("test integral", q1, q2, scale );
      if (r && debug)  {
         std::cout << "Failed test for x = " << v1 << " q1= " << q1 << " q2= " << q2 << " p = ";
         for (int j = 0; j < NPAR; ++j) std::cout << fParams[j] << "\t";
         std::cout << "ig error is " << err << " status " << ig.Status() << std::endl;
      }
      iret |= r;

   }
   return iret;

}

int StatFunction::TestDerivative() {

   int iret = 0;

   // scan all values from xmin to xmax
   double dx = (xmax-xmin)/NFuncTest;
   // create CDF function
   Functor1D func(this, &StatFunction::Cdf);
   Derivator d(func);

   for (int i = 0; i < NFuncTest; ++i) {
      double v1 = xmin + dx*i;  // value used  for testing
      double q1 = Pdf(v1);
      //std::cout << "v1 " << v1 << " pdf " << (*this)(v1) << " cdf " << q1 << " quantile " << Quantile(q1) << std::endl;
      // calculate derivative of cdf
      double q2 = 0;
      if (fHasLowRange  && v1 == xlow)
         q2 = d.EvalForward(v1);
      else if (fHasUpRange && v1 == xup)
         q2 = d.EvalBackward(v1);
      else
         q2 = d.Eval(v1);

      int r = d.Status();
      double err = d.Error();

      double scale = std::max(1.,fScaleDer * err / std::numeric_limits<double>::epsilon() );


      r |= compare("test Derivative", q1, q2, scale );
      if (r && debug)  {
         std::cout << "Failed test for x = " << v1 << " p = ";
         for (int j = 0; j < NPAR; ++j) std::cout << fParams[j] << "\t";
         std::cout << "der error is " << err << std::endl;
         std::cout << d.Eval(v1) << "\t" << d.EvalForward(v1) << std::endl;
      }
      iret |= r;

   }
   return iret;

}

// function to be used in ROOT finding algorithm
struct InvFunc {
   InvFunc(const StatFunction * f, double y) : fFunc(f), fY(y)  {}
   double operator() (double x) {
      return fFunc->Cdf(x) - fY;
   }
   const StatFunction * fFunc;
   double fY;
};


int StatFunction::TestInverse1(RootFinder::EType algoType) {

   int iret = 0;
   int maxitr = 2000;
   double abstol = 1.E-15;
   double reltol = 1.E-15;
   //NFuncTest = 4;


   // scan all values from 0.05 to 0.95  to avoid problem at the border of definitions
   double x1 = 0.05; double x2 = 0.95;
   double dx = (x2-x1)/NFuncTest;
   double vmin = Quantile(dx/2);
   double vmax = Quantile(1.-dx/2);

   // test ROOT finder algorithm function without derivative
   RootFinder  rf1(algoType);
   for (int i = 1; i < NFuncTest; ++i) {
      double v1 = x1 + dx*i;  // value used  for testing
      InvFunc finv(this,v1);
      Functor1D func(finv);
      rf1.SetFunction(func, vmin, vmax);
      //std::cout << "\nfun values for :" << v1 << " f:  " << func(0.0) << "  " << func(1.0) << std::endl;
      int ret = ! rf1.Solve(maxitr,abstol,reltol);
      if (ret && debug) {
         std::cout << "\nError in solving for inverse, niter = " << rf1.Iterations() << std::endl;
      }
      double q1 = rf1.Root();
      // test that quantile value correspond:
      double q2 = Quantile(v1);

      ret |= compare("test Inverse1", q1, q2, fScaleInv );
      if (ret && debug)  {
         std::cout << "\nFailed test for x = " << v1 << " p = ";
         for (int j = 0; j < NPAR; ++j) std::cout << fParams[j] << "\t";
         std::cout << std::endl;
      }
      iret |= ret;

   }
   return iret;

}

int StatFunction::TestInverse2(RootFinder::EType algoType) {

   int iret = 0;
   int maxitr = 2000;
   // put lower tolerance
   double abstol = 1.E-12;
   double reltol = 1.E-12;
   //NFuncTest = 10;

   // scan all values from 0.05 to 0.95  to avoid problem at the border of definitions
   double x1 = 0.05; double x2 = 0.95;
   double dx = (x2-x1)/NFuncTest;
   // starting root is always on the left to avoid to go negative
   // it is very sensible at the starting point
   double vstart = fStartRoot; //depends on function shape
   // test ROOT finder algorithm function with derivative
   RootFinder rf1(algoType);
   //RootFinder<Roots::Secant> rf1;

   for (int i = 1; i < NFuncTest; ++i) {
      double v1 = x1 + dx*i;  // value used  for testing

      InvFunc finv(this,v1);
      //make a gradient function using inv function and derivative (which is pdf)
      GradFunctor1D func(finv,*this);
      // use as estimate the quantile at 0.5
      //std::cout << "\nvstart : " << vstart << " fun/der values" << func(vstart) << "  " << func.Derivative(vstart) << std::endl;
      rf1.SetFunction(func,vstart );
      int ret = !rf1.Solve(maxitr,abstol,reltol);
      if (ret && debug) {
         std::cout << "\nError in solving for inverse using derivatives,  niter = " << rf1.Iterations() << std::endl;
      }
      double q1 = rf1.Root();
      // test that quantile value correspond:
      double q2 = Quantile(v1);

      ret |= compare("test InverseDeriv", q1, q2, fScaleInv );
      if (ret && debug)  {
         std::cout << "Failed test for x = " << v1 << " p = ";
         for (int j = 0; j < NPAR; ++j) std::cout << fParams[j] << "\t";
         std::cout << std::endl;
      }
      iret |= ret;
   }
   return iret;
}

// test intergal. derivative and inverse(Rootfinder)
int testGammaFunction(int n = 100) {

   int iret = 0;

   Func4 pdf(gamma_pdf);
   Func4 cdf(gamma_cdf);
   Func3 quantile(gamma_quantile);
   StatFunction dist(pdf, cdf, quantile, 0.);
   dist.SetNTest(n);
   dist.SetTestRange(0.,10.);
   dist.SetScaleDer(10); // few tests fail here
   dist.SetScaleInv(10000); // few tests fail here
   // vary shape of gamma parameter
   for (int i =1; i <= 5; ++i) {
      double k = std::pow(2.,double(i-1));
      double theta = 2./double(i);
      dist.SetParameters(k,theta);
      if (k <=1 )
         dist.SetStartRoot(0.1);
      else
         dist.SetStartRoot(k*theta-1.);

      std::string name = "Gamma("+Util::ToString(int(k))+","+Util::ToString(theta)+") ";
      std::cout << "\nTest " << name << " distribution\n";
      int ret = 0;

      PrintTest("\t test integral GSL adaptive");
      ret = dist.TestIntegral(IntegrationOneDim::kADAPTIVESINGULAR);
      PrintStatus(ret);
      iret |= ret;

      PrintTest("\t test integral Gauss");
      dist.SetScaleIg(100); // relax for Gauss integral
      ret = dist.TestIntegral(IntegrationOneDim::kGAUSS);
      PrintStatus(ret);
      iret |= ret;

      PrintTest("\t test derivative");
      ret = dist.TestDerivative();
      PrintStatus(ret);
      iret |= ret;

      PrintTest("\t test inverse with GSL Brent method");
      ret = dist.TestInverse1(RootFinder::kGSL_BRENT);
      PrintStatus(ret);
      iret |= ret;

      PrintTest("\t test inverse with Steffenson algo");
      ret = dist.TestInverse2(RootFinder::kGSL_STEFFENSON);
      PrintStatus(ret);
      iret |= ret;

      PrintTest("\t test inverse with Brent method");
      ret = dist.TestInverse1(RootFinder::kBRENT);
      PrintStatus(ret);
      iret |= ret;

   }

   return iret;
}

// test intergal. derivative and inverse(Rootfinder)
int testBetaFunction(int n = 100) {

   int iret = 0;

   Func3 pdf(beta_pdf);
   Func3 cdf(beta_cdf);
   Func3 quantile(beta_quantile);
   StatFunction dist(pdf, cdf, quantile, 0.,1.);

   dist.SetNTest(n);
   dist.SetTestRange(0.,1.);
   // vary shape of beta function parameters
   for (int i = 0; i < 5; ++i) {
      // avoid case alpha or beta = 1
      double alpha = i+2;
      double beta = 6-i;
      dist.SetParameters(alpha,beta);
      dist.SetStartRoot(alpha/(alpha+beta)); // use mean value

      std::string name = "Beta("+Util::ToString(int(alpha))+","+Util::ToString(beta)+") ";
      std::cout << "\nTest " << name << " distribution\n";
      int ret = 0;

      PrintTest("\t test integral GSL adaptive");
      ret = dist.TestIntegral(IntegrationOneDim::kADAPTIVESINGULAR);
      PrintStatus(ret);
      iret |= ret;

      PrintTest("\t test integral Gauss");
      dist.SetScaleIg(100); // relax for Gauss integral
      ret = dist.TestIntegral(IntegrationOneDim::kGAUSS);
      PrintStatus(ret);
      iret |= ret;

      PrintTest("\t test derivative");
      ret = dist.TestDerivative();
      PrintStatus(ret);
      iret |= ret;

      PrintTest("\t test inverse with Brent method");
      ret = dist.TestInverse1(RootFinder::kBRENT);
      PrintStatus(ret);
      iret |= ret;

      PrintTest("\t test inverse with GSL Brent method");
      ret = dist.TestInverse1(RootFinder::kGSL_BRENT);
      PrintStatus(ret);
      iret |= ret;

      if (i < 5) {  // test failed for k=5
         PrintTest("\t test inverse with Steffenson algo");
         ret = dist.TestInverse2(RootFinder::kGSL_STEFFENSON);
         PrintStatus(ret);
         iret |= ret;
      }
   }

   return iret;
}


int stressMathMore(double nscale = 1) {

   int iret = 0;

#ifdef __CINT__
   std::cout << "Test must be run in compile mode - please use ACLIC !!" << std::endl;
   return 0;
#endif

   TBenchmark bm;
   bm.Start("stressMathMore");

   const int ntest = 10000;
   int n = int(nscale*ntest);
   //std::cout << "StressMathMore: test number  n = " << n << std::endl;

   iret |= testGammaFunction(n);
   iret |= testBetaFunction(n);

   bm.Stop("stressMathMore");
   std::cout <<"******************************************************************************\n";
   bm.Print("stressMathMore");
   const double reftime = 7.24; //to be updated  // ref time on  pcbrun4
   double rootmarks = 860 * reftime / bm.GetCpuTime("stressMathMore");
   std::cout << " ROOTMARKS = " << rootmarks << " ROOT version: " << gROOT->GetVersion() << "\t"
             << gROOT->GetGitBranch() << "@" << gROOT->GetGitCommit() << std::endl;
   std::cout <<"*******************************************************************************\n";


   if (iret !=0) std::cerr << "stressMathMore Test Failed !!" << std::endl;
   return iret;
}


int main(int argc,const char *argv[]) {
   double nscale = 1;
   if (argc > 1) {
      nscale = atof(argv[1]);
      //nscale = std::pow(10.0,double(scale));
   }
   return stressMathMore(nscale);
}
