#include "TF1.h"
#include "Math/Functor.h"
#include "TStopwatch.h"
#include "RConfigure.h"

#include "Math/RootFinder.h"
#include "Math/DistFunc.h"

#include <iostream>
#include <cmath>

#include "gtest/gtest.h"

/**
   Test of ROOT finder for various function

   case = 0  simple function (polynomial)
   case = 1  function which fails for a bug in BrentMethod::MinimBrent fixed with r = 32544 (5.27.01)
 */

const double ERRORLIMIT = 1E-8;
int iterTest = 100;
int myfuncCalls = 0;

double Y0_P2 = 5.0;   // Y0 for test 1 (parabola)
// these are value which gave problems in 5.26 for gamma_cdf
const double Y0_GAMMA = 0.32;  // Y0 for test 2 (gamma cdf)
const double ALPHA_GAMMA = 16.; // alpha of gamma cdf
const double THETA_GAMMA = 0.4;  // theta of gamma cdf

bool debug =  false; 

int gTestCase = 0;

double myfunc ( double x ) {
   myfuncCalls ++;
   if (gTestCase == 0) // polynomial
      return x*x - Y0_P2;
   if (gTestCase == 1) // use function in logx
      return log(x)*log(x) - Y0_P2; 
   if (gTestCase == 2) // function showing bug in BrentMethod::
      return ROOT::Math::gamma_cdf(x,ALPHA_GAMMA,THETA_GAMMA)-Y0_GAMMA;
   return 0;
}

double ExactResult(double y0 = 0, int type = 1) {
   if (gTestCase == 0) {
      return type * std::sqrt(Y0_P2 + y0);
   }
   if (gTestCase == 1) {
      return std::exp( type * std::sqrt(Y0_P2 + y0) );
   }
   if (gTestCase == 2) {
      if (y0 == 0)  return 5.55680381022934800;
#ifdef R__HAS_MATHMORE
      return ROOT::Math::gamma_quantile(y0 + Y0_GAMMA,ALPHA_GAMMA,THETA_GAMMA);
#endif
   }
   return 0;
}


double myfunc_p (double *x, double *) { return myfunc(x[0]); }

void printStats(TStopwatch& timer, double y0, double root) {

   //std::cout << "Return code:  " << status << std::endl;
   if (debug) { 
      double difference = root - ExactResult(y0);
      int pr = std::cout.precision(16);
      //std::cout << "y0 :       " << y0 << std::endl;
      std::cout << "A Result :       " << root << std::endl;
      std::cout << "Exact result: " << ExactResult(y0);
      std::cout.precision(pr);
      std::cout << " difference: " << difference << std::endl;
   }
   std::cout << "Average Time: " << timer.RealTime()/(2. * iterTest) << std::endl;
   std::cout << "Average Number of calls to function: " << double(myfuncCalls)/(2. * iterTest) << std::endl;

   return; // difference > ERRORLIMIT;
}

void runTestTF1(int testcase = 0) {

   double xmin = -10;
   double xmax = 10;
   bool logx = false;
   Y0_P2 = 5; 

   std::cout << "*************************************************************\n";
   gTestCase = testcase;
   if (gTestCase == 0)
      std::cout << "Test for parabola function f(x) = x^2 - 5 " << std::endl;
   if (gTestCase == 1) {
      std::cout << "Test for parabola log function f(x) = (logx)^2 - 5" << std::endl;
      xmin = 1.E-6;
      xmax = 1.E6;
      logx = true;
   }
   if (gTestCase == 2) {
      std::cout << "Test for  f(x) = gamma_cdf  " << std::endl;
      xmin = 3.955687382047723;
      xmax = 9.3423159494328623;
      iterTest = 1; 
   }
   std::cout << "\t TF1::GetX()" << std::endl;
   std::cout << "*************************************************************\n";

   TStopwatch timer;
   double root = 0.;
   //int status = 0;

   double tol = 1.E-14;
   int maxiter = 100;

   TF1  f1("f1", myfunc_p, xmin, xmax);
   timer.Reset(); timer.Start(); myfuncCalls = 0;
   double y0 = 0;
   double delta = 0;
   if (gTestCase == 0) 
      delta = (xmax*xmax)/double(iterTest)/10.;
   else if (gTestCase == 1) 
      delta = (log(xmax)*log(xmax))/double(iterTest)/10.;


   for (int i = 0; i < iterTest; ++i)
   {
            
      if (gTestCase == 0) {
         y0 = (xmax*xmax)/double(iterTest)*(i+0.5) - Y0_P2; 
      }
      else if (gTestCase == 1) {
         y0 = (log(xmax)*log(xmax))/double(iterTest)*(i+0.5) - Y0_P2; 
      }
      
      double root1 = f1.GetX(y0, xmin, xmax,tol,maxiter,logx);
      EXPECT_NEAR(root1,  ExactResult(y0,-1), ERRORLIMIT );
      root = root1; 
      // for the parabola test cases find also second root
      if (gTestCase < 2) { 
         double root2 = f1.GetX(y0, root1+ delta, xmax,tol,maxiter, logx);
         EXPECT_NEAR(root2, ExactResult(y0, 1), ERRORLIMIT);
         root = root2;
         
         if (debug) std::cout << "tested #" << i << " y0=" << y0 << " in ["<< xmin << "," << xmax << "] and ["
                              << root1+delta << "," << xmax << "]  x : f(x)=y0 is " << root1  << " and " << root << std::endl;
      }
      else {
         if (debug) std::cout << "tested #" << i << " y0=" << y0 << " in ["<< xmin << "," << xmax
                              <<  "]  x : f(x)=y0 is " << root1  << std::endl; 
      }

   }
   timer.Stop();

   printStats(timer, y0, root);

   return;
}


void runTestBrent(int testcase = 0, ROOT::Math::RootFinder::EType rf_type = ROOT::Math::RootFinder::kBRENT) {

   double xmin = -10;
   double xmax = 10;

   std::cout << "*************************************************************\n";
   gTestCase = testcase;
   if (gTestCase == 0)
      std::cout << "Test for parabola function f(x) = x^2 + C " << std::endl;
   if (gTestCase == 1) { 
      std::cout << "Test for parabola log function f(x) = (logx)^2 - 5" << std::endl;
      xmin = 1.E-6;
      xmax = 1.E6;
   }
   if (gTestCase == 2) {
      std::cout << "Test for  f(x) = gamma_cdf " << std::endl;
      xmin = 3.955687382047723;
      xmax = 9.3423159494328623;
   }
   // Only BrentRootFinder will accept xmax1 = xmax for th eparabola 
   double xmax1 = xmax;
   if (rf_type != ROOT::Math::RootFinder::kBRENT) {
      if (gTestCase == 0) xmax1 = 0; 
      if (gTestCase == 1) xmax1 = 1; 
      if (gTestCase == 2) xmax1 = xmax; 
   }

   TStopwatch timer;
   double root = 0.;
   //int status = 0;

   double tol = 1.E-14;
   int maxiter = 100;

   ROOT::Math::Functor1D  func(&myfunc);

   double y0 = 0;
   double delta = 0;
   if (gTestCase == 0) 
      delta = (xmax*xmax)/double(iterTest)/10.;
   else if (gTestCase == 1) 
      delta = (log(xmax)*log(xmax))/double(iterTest)/10.;

   ROOT::Math::RootFinder brf(rf_type);

   std::cout << "\t RootFinder : " << brf.Name() << std::endl;
   std::cout << "*************************************************************\n";
   
   timer.Reset(); timer.Start(); myfuncCalls = 0;
   for (int i = 0; i < iterTest; ++i)
   {
      if (gTestCase == 0) {
         y0 = (xmax*xmax)/double(iterTest)*(i+0.5) - 5; 
         Y0_P2 = 5.0 + y0; 
      }
      else if (gTestCase == 1) {
         y0 = (log(xmax)*log(xmax))/double(iterTest)*(i+0.5) - 5.; 
         Y0_P2 = 5.0 + y0; 
      }
      brf.SetFunction( func, xmin, xmax1 );
      bool ret1 = brf.Solve(maxiter,tol,tol);
      if (!ret1 ) std::cout << "Error returned from RootFinder::Solve BRENT - interval [ "
                           << xmin << " , " << xmax1 << " ]  i = " << i << std::endl;
      ASSERT_EQ(ret1,true);
      double root1 = brf.Root();
      EXPECT_NEAR(root1,  ExactResult(0, -1), ERRORLIMIT );
      root = root1; 
      if (gTestCase < 2) { 
         brf.SetFunction( func, root1+delta, xmax );
         bool ret2 = brf.Solve(maxiter,tol,tol);
         if (!ret2) std::cout << "Error returned from RootFinder::Solve BRENT - interval [ "
                              << root1+delta << " , " << xmax << " ]  i = " << i << std::endl;
         ASSERT_EQ(ret2,true);
         double root2 = brf.Root();
         EXPECT_NEAR(root2, ExactResult(0, 1), ERRORLIMIT);
         root = root2;
         if (debug) std::cout << "tested #" << i << " y0=" << y0 << " in ["<< xmin << "," << xmax1 << "] and ["
                              << root1+delta << "," << xmax << "]  x : f(x)=y0 is " << root1  << " and " << root << std::endl; 
      }
      else {
         if (debug) std::cout << "tested #" << i << " y0=" << y0 << " in ["<< xmin << "," << xmax1
                              <<  "]  x : f(x)=y0 is " << root1  << std::endl; 
      }
   }
   timer.Stop();
   printStats(timer, 0, root);

   return;

}


/*
int testRootFinder() {

   int status = 0;
   status |= runTest(0);  // test pol function
   if (status) std::cerr << "Test pol function  FAILED" << std::endl;

   status |= runTest(1);  // test gamma_cdf
   if (status) std::cerr << "Test gamma function  FAILED" << std::endl;

   std::cerr << "*************************************************************\n";
   std::cerr << "\nTest RootFinder :\t";
   if (status == 0)
      std::cerr << "OK " << std::endl;
  else
     std::cerr << "Failed !" << std::endl;

   return status;
}
*/

//_________________________________________
// Parabola tests

TEST(Parabola,TF1GetX)
{
   //test parabola using TF1 
   runTestTF1(0);      
}

TEST(Parabola,BrentRootFinder)
{
   // test parabola using Brent Root Finder 
   runTestBrent(0);      
}

#ifdef R__HAS_MATHMORE
TEST(Parabola, GSLBrent)
{
   //test parabola using GSL Brent Root Finder 
   runTestBrent(0, ROOT::Math::RootFinder::kGSL_BRENT );      
}
TEST(Parabola, GSLBisection)
{
   //test parabola using GSL Brent Root Finder 
   runTestBrent(0, ROOT::Math::RootFinder::kGSL_BISECTION );      
}
TEST(Parabola, GSLFalsePos)
{
   //test parabola using GSL Brent Root Finder 
   runTestBrent(0, ROOT::Math::RootFinder::kGSL_FALSE_POS );      
}
// need to add test for GSL algorithms with derivatives
#endif

//_________________________________________
// Log-Parabola tests

TEST(LogParabola,TF1GetX)
{

   // test log-parabola using TF1
   runTestTF1(1);      
}

TEST(LogParabola,BrentRootFinder)
{
   // test parabola using Brent Root Finder 
   runTestBrent(0);      
}

#ifdef R__HAS_MATHMORE
TEST(LogParabola,GSLBrent)
{
   // test gamma cdf using GSL Brent
   runTestBrent(1, ROOT::Math::RootFinder::kGSL_BRENT);      
}
#endif

//_________________________________________
// GammaCDF tests

TEST(GammaCDF,TF1GetX)
{
   // test gamma cdf using TF1 
   runTestTF1(2);      
}

TEST(GammaCDF,BrentRootFinder)
{
   // test gamma cdf using Brent
   runTestBrent(2);      
}

#ifdef R__HAS_MATHMORE
TEST(GammaCDF,GSLBrent)
{
   // test gamma cdf using GSL Brent
   runTestBrent(2, ROOT::Math::RootFinder::kGSL_BRENT);      
}
TEST(GammaCDF, GSLBisection)
{
   //test parabola using GSL Brent Root Finder 
   runTestBrent(2, ROOT::Math::RootFinder::kGSL_BISECTION );      
}
TEST(GammaCDF, GSLFalsePos)
{
   //test parabola using GSL Brent Root Finder 
   runTestBrent(2, ROOT::Math::RootFinder::kGSL_FALSE_POS );      
}
#endif

int main(int argc, char **argv) {
   // Parse command line arguments
   for (Int_t i = 1 ;  i < argc ; i++) {
      std::string arg = argv[i] ;
      if (arg == "-v" || arg == "-vv") {
         std::cout << "running in verbose mode" << std::endl;
         debug = true;
      }
   }
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
