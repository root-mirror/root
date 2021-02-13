#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "Math/PdfFuncMathCore.h"

#include "gtest/gtest.h"

#include <iostream>

using namespace std;

bool debug = false;

void testTF1Moments(TF1 &f1, double eps = 1.E-6)
{

   double par[3] = {1, 1.5, 1.2};

   //f1.SetParameters(par);

   double mx1 = f1.Mean(-10, 10, par, eps);

   double mx2 = f1.Moment(2, -10, 10, par, eps);

   double varx = f1.Variance(-10, 10, par, eps);

   double delta = 100 * eps;

   EXPECT_NEAR(mx1 , par[1], delta);

   EXPECT_NEAR(sqrt(varx ), par[2], delta);
   EXPECT_NEAR(sqrt(mx2 - mx1 * mx1 ), par[2], delta);
}

// Test the moments for a TF2 function - a bi-variate gaussian with correlation

void testTF2Moments(TF2 & f1, double eps = 1.E-6) {

   double par[6] = {1, 2., 0.8, -1., 2., 0.3};

   f1.SetParameters(par);

   if (debug) {
      f1.Print("v");
      gDebug = 1;
      std::cout << f1.Eval(0, 0) << std::endl;
      return;
   }


   double mx1 = f1.Mean2X(-10, 10, -10, 10, eps);
   double my1 = f1.Mean2Y(-10, 10, -10, 10, eps);

   double mx2 = f1.Moment2(2, -10, 10, 0, -10, 10., eps);
   double my2 = f1.Moment2(0, -10, 10, 2, -10, 10., eps);

   double varx = f1.Variance2X(-10, 10, -10, 10, eps);
   double vary = f1.Variance2Y(-10, 10, -10, 10, eps);


   double covxy = f1.Covariance2XY(-10, 10, -10, 10, eps);
   double corrxy = covxy / sqrt(varx * vary );

   gDebug = 0;


   if (debug) {
      std::cout << "Results for function " << f1.GetName() << std::endl;


      std::cout << mx1 << " , " << my1 << std::endl;
      std::cout << sqrt(mx2 - mx1 * mx1 ) << " , " << sqrt(my2 - my1 * my1 ) << std::endl;
      std::cout << sqrt(varx ) << " , " << sqrt(vary ) << std::endl;
      std::cout << corrxy << std::endl;
   }

   double delta = 100*eps;
   EXPECT_NEAR(mx1 , par[1], delta);
   EXPECT_NEAR(my1 , par[3], delta);

   EXPECT_NEAR(sqrt(varx ), par[2], delta);
   EXPECT_NEAR(sqrt(vary ), par[4], delta);

   EXPECT_NEAR(sqrt(mx2 - mx1 * mx1 ), par[2], delta);
   EXPECT_NEAR(sqrt(my2 - my1 * my1 ), par[4], delta);

   EXPECT_NEAR(corrxy , par[5], delta);

}

void testTF3Moments(TF3 &f1, double eps = 1.E-6)
{

   double par[9] = {1, 2., 0.5, -1., 2., -0.4, 1., 1.5, 1.3};

   f1.SetParameters(par);

   // increase number of points for integral
   f1.SetNpx(100);
   f1.SetNpy(100);
   f1.SetNpz(100);

   if (debug) {
      f1.Print("v");
      gDebug = 1;
      return;
   }

   double mx1 = f1.Mean3X(-10, 10, -10, 10, -10, 10, eps);
   double my1 = f1.Mean3Y(-10, 10, -10, 10, -10, 10, eps);
   double mz1 = f1.Mean3Z(-10, 10, -10, 10, -10, 10, eps);

   double mx2 = f1.Moment3(2, -10, 10, 0, -10, 10., 0., -10, 10, eps);
   double my2 = f1.Moment3(0, -10, 10, 2, -10, 10., 0., -10., 10, eps);
   double mz2 = f1.Moment3(0, -10, 10, 0, -10, 10., 2, -10, 10, eps);

   double varx = f1.Variance3X(-10, 10, -10, 10, -10, 10, eps);
   double vary = f1.Variance3Y(-10, 10, -10, 10, -10, 10, eps);
   double varz = f1.Variance3Z(-10, 10, -10, 10, -10, 10, eps);

   double covxy = f1.Covariance3XY(-10, 10, -10, 10, -10, 10, eps);
   double covxz = f1.Covariance3XZ(-10, 10, -10, 10, -10, 10, eps);
   double covyz = f1.Covariance3YZ(-10, 10, -10, 10, -10, 10, eps);

   double corrxy = covxy / sqrt(varx * vary );

   gDebug = 0;

   if (debug) {
      std::cout << "Results for function " << f1.GetName() << std::endl;

      std::cout << mx1 << " , " << my1 << "  " << mz1 << std::endl;
      std::cout << sqrt(mx2 - mx1 * mx1 ) << " , " << sqrt(my2 - my1 * my1 ) << " , "
                << sqrt(mz2 - mz1 * mz1 ) << std::endl;
      std::cout << sqrt(varx ) << " , " << sqrt(vary ) << " , " << sqrt(varz ) << std::endl;
      std::cout << corrxy << " , " << covxz << " , " << covyz << std::endl;
   }

   double delta = 100*eps;
   EXPECT_NEAR(mx1 , par[1], delta);
   EXPECT_NEAR(my1 , par[3], delta);
   EXPECT_NEAR(mz1 , par[7], delta);

   EXPECT_NEAR(sqrt(varx ), par[2], delta);
   EXPECT_NEAR(sqrt(vary ), par[4], delta);
   EXPECT_NEAR(sqrt(varz ), par[8], delta);

   EXPECT_NEAR(sqrt(mx2 - mx1 * mx1 ), par[2], delta);
   EXPECT_NEAR(sqrt(my2 - my1 * my1 ), par[4], delta);
   EXPECT_NEAR(sqrt(mz2 - mz1 * mz1 ), par[8], delta);

   EXPECT_NEAR(corrxy , par[5], delta);
   EXPECT_NEAR(covxz , 0, delta);
   EXPECT_NEAR(covyz , 0, delta);
}

TEST(TF1Moments, FormulaTF1)
{
   TF1 f1("gaus_formula", "gaus", -10, 10);
   testTF1Moments(f1);
}

TEST(TF2Moments, FormulaTF2)
{
   TF2 f1("bigaus_formula", "bigaus", -10, 10, -10, 10);
   testTF2Moments(f1);
}

TEST(TF2Moments, CompiledTF2)
{
   auto bgaus = [](double *x, double *p) {
      return p[0] * ROOT::Math::bigaussian_pdf(x[0], x[1], p[2], p[4], p[5], p[1], p[3]);
   };

   TF2 f2("bigaus_lambda", bgaus, -10, 10, -10, 10, 6);
   testTF2Moments(f2);
}

TEST(TF3Moments, FormulaTF3)
{
   TF3 f1("trigaus_formula", "bigaus*gausn(z,[6],[7],[8])", -10, 10, -10, 10, -10, 10);
   testTF3Moments(f1);
}

TEST(TF3Moments, CompiledTF3)
{
   auto tgaus = [](double *x, double *p) {
      return p[0] * ROOT::Math::bigaussian_pdf(x[0], x[1], p[2], p[4], p[5], p[1], p[3]) * p[6] *
             ROOT::Math::normal_pdf(x[2], p[8], p[7]);
   };

   TF3 f2("f2", tgaus, -10, 10, -10, 10, -10, 10, 9);
   testTF3Moments(f2);
}
