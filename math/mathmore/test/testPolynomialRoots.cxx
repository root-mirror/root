// test finding Roots of polynomials

#include "gtest/gtest.h"

#include "Math/Polynomial.h"
#include <vector>
#include <complex>

using namespace ROOT::Math;
using namespace std;

class QuarticPolynomial : public ::testing::Test {

protected:
   enum EModeType { kAnalytical = 1, kNumerical = 2 };
   vector<double> coeff;
   vector<complex<double>> expectedResult;
   double a, b, c, d, e = 0;
   double tolerance = 0;
   bool debug = false;
   EModeType algoType;

   void checkResult(vector<complex<double>> &result)
   {
      // sort first result
      auto smallerComplex = [&](complex<double> c1, complex<double> c2) {
         double diff = std::abs( c1.real()-c2.real() );
         // in case of numerical roots difference can be small
         //if (c1.real() != c2.real())
         if (diff > std::max(tolerance, 4 * std::numeric_limits<double>::epsilon() ) )
            return c1.real() < c2.real();
         else
            return c1.imag() < c2.imag();
      };
      std::sort(result.begin(), result.end(), smallerComplex);

      ASSERT_EQ(result.size(), expectedResult.size());
      ASSERT_EQ(result.size(), 4 );

      // unrolll loop to have test printing component when failed
      if (tolerance <= 0) {
         // in this case check with 4 ulp
         EXPECT_DOUBLE_EQ(result[0].real(), expectedResult[0].real());
         EXPECT_DOUBLE_EQ(result[0].imag(), expectedResult[0].imag());

         EXPECT_DOUBLE_EQ(result[1].real(), expectedResult[1].real());
         EXPECT_DOUBLE_EQ(result[1].imag(), expectedResult[1].imag());

         EXPECT_DOUBLE_EQ(result[2].real(), expectedResult[2].real());
         EXPECT_DOUBLE_EQ(result[2].imag(), expectedResult[2].imag());

         EXPECT_DOUBLE_EQ(result[3].real(), expectedResult[3].real());
         EXPECT_DOUBLE_EQ(result[3].imag(), expectedResult[3].imag());
      }
      else {
         // check within a given absolute tolerance
         EXPECT_NEAR(result[0].real(), expectedResult[0].real(), tolerance);
         EXPECT_NEAR(result[0].imag(), expectedResult[0].imag(), tolerance);

         EXPECT_NEAR(result[1].real(), expectedResult[1].real(), tolerance);
         EXPECT_NEAR(result[1].imag(), expectedResult[1].imag(), tolerance);

         EXPECT_NEAR(result[2].real(), expectedResult[2].real(), tolerance);
         EXPECT_NEAR(result[2].imag(), expectedResult[2].imag(), tolerance);

         EXPECT_NEAR(result[3].real(), expectedResult[3].real(),tolerance);
         EXPECT_NEAR(result[3].imag(), expectedResult[3].imag(),tolerance);
      }
   }

   void print(const vector<complex<double>> &result)
   {
      std::string algoName = (algoType == kNumerical) ? "Numerical" : "Analytical";
      std::cout << algoName << " solution for " << a << "*x^4 + " << b << "*x^3 + " << c << "*x^2 + " << d << "*x + " << e << " = 0"
                << std::endl;
      for (unsigned int i = 0; i < result.size(); i++)
      {
         std::cout << " root " << i << " : " << result[i] << std::endl;
      }
   }

      void runTest(EModeType type, double tol = 0)
   {
      a = coeff[0];
      b = coeff[1];
      c = coeff[2];
      d = coeff[3];
      e = coeff[4];

      tolerance = tol;
      algoType = type;

      Polynomial pol(a, b, c, d, e);

      auto result = (type != kNumerical) ? pol.FindRoots() : pol.FindNumRoots();

      checkResult(result);

      if (debug)
         print(result);
   }
};

// test x^4 - 16 = 0
TEST_F(QuarticPolynomial, FindRoots_DegPol1)
{
   coeff = {1, 0, 0, 0, -16};

   expectedResult = {complex<double>(-2, 0),
                     complex<double>( 0,-2),
                     complex<double>( 0, 2),
                     complex<double>( 2, 0)};

   runTest(kAnalytical);
   runTest(kNumerical, 1.E-13);
}

// test (x+1)^4 = 0
TEST_F(QuarticPolynomial, FindRoots_DegPol2)
{
   coeff = {1, 4, 6, 4, 1};

   expectedResult = {complex<double>(-1, 0),
                     complex<double>(-1, 0),
                     complex<double>(-1, 0),
                     complex<double>(-1, 0)};

   runTest(kAnalytical);
   // bad tolerance for this test
   debug = true;
   // numerical method has a large error for this case
   runTest(kNumerical, 1.E-3);
}

// test x4 + 5x^2 + 4 = 0
//4 imaginary roots (x-i)(x+i)(x-2i)(x+2i)=0
TEST_F(QuarticPolynomial, FindRoots_4ImagRoots)
{
   coeff = {1, 0, 5, 0, 4};

   expectedResult = {complex<double>(0, -2),
                     complex<double>(0, -1),
                     complex<double>(0, 1),
                     complex<double>(0, 2)};

   //debug = true;
   runTest(kAnalytical, 0);
   runTest(kNumerical, 1.E-12);
}

//
// four full complex roots, (x-1-i)(x-1+i)(x-1-2i)(x-1+2i)=0")
TEST_F(QuarticPolynomial, FindRoots_4CompRoots)
{
   coeff = {1.0, -4.0, 11.0, -14.0, 10.0};

   expectedResult = {complex<double>(1, -2),
                     complex<double>(1, -1),
                     complex<double>(1, 1),
                     complex<double>(1, 2)};

   runTest(kAnalytical, 0);
   runTest(kNumerical, 1.E-12);
}

// four real roots (x+1)(x+2)(x+3)(x+4)=0
TEST_F(QuarticPolynomial, FindRoots_4RealRoots)
{

   coeff = {1.0, 10.0, 35.0, 50.0, 24.0};

   expectedResult = {complex<double>(-4, 0),
                     complex<double>(-3, 0),
                     complex<double>(-2, 0),
                     complex<double>(-1, 0)};

   runTest(kAnalytical, 0);
   runTest(kNumerical, 1.E-12);
}

// test   x^4-8x^3+12x^2+16x+4=0
// 4 real roots where 2 are degeenrates

TEST_F(QuarticPolynomial, FindRoots_4RealDegRoots)
{

   coeff = {1.0, -8.0, 12.0, 16.0, 4.0};

   expectedResult = {complex<double>(-0.44948974278317788, 0),
                     complex<double>(-0.44948974278317788, 0),
                     complex<double>(4.44948974278317788, 0),  // this is 4 - r0
                     complex<double>(4.44948974278317788, 0)};

   runTest(kAnalytical, 0);
   debug = true;
   runTest(kNumerical, 1.E-7);
}

// test 5x^4 + 4x^3 + 3x^2 + 2x + 1

TEST_F(QuarticPolynomial, FindRoots_54321)
{
   coeff = {5, 4, 3, 2, 1};

   expectedResult = {complex<double>(-0.5378322749029899, -0.358284686345128),
                     complex<double>(-0.5378322749029899, +0.358284686345128),
                     complex<double>(+0.13783227490298988, -0.6781543891053364),
                     complex<double>(+0.13783227490298988, +0.6781543891053364)};

   runTest(kAnalytical);
   runTest(kNumerical, 1.E-13);
}
// special case reported in issue #6900 by S. Binet
TEST_F(QuarticPolynomial, FindRoots_4RealDegRootsR0)
{
   coeff = {2.2206846808021337, 7.643281053997895, 8.831759446092846, 3.880673545129404, 0.5724551380144077};

   expectedResult = {complex<double>(-1.3429253, 0.),
                     complex<double>(-1.3427327, 0 ),
                     complex<double>(-0.3781962, 0 ),
                     complex<double>(-0.3780036, 0) };

   debug = true;
   runTest(kAnalytical, 1.E-6);
   runTest(kNumerical, 1.E-4);
}
