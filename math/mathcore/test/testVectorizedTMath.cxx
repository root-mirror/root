#include <cmath>
#include <random>
#include "VectorizedTMath.h"
#include <gtest/gtest.h>

#define N 16384

Double_t uniform_random(Double_t a, Double_t b)
{
   return a + (b - a) * drand48();
}

class VectorizedTMathTest : public ::testing::Test {
protected:
   VectorizedTMathTest() {}

   size_t kVS = vecCore::VectorSize<ROOT::Double_v>();
   Double_t input_array1[N] __attribute__((aligned(VECCORE_SIMD_ALIGN)));
   Double_t input_array2[N] __attribute__((aligned(VECCORE_SIMD_ALIGN)));

   Double_t output_array[N] __attribute__((aligned(VECCORE_SIMD_ALIGN)));
};

#define TEST_VECTORIZED_TMATH_FUNCTION(tmathfunc, a, b)                                                             \
   TEST_F(VectorizedTMathTest, tmathfunc)                                                                           \
   {                                                                                                                \
      int trials = N;                                                                                           \
      for (int i = 0; i < trials; i++)                                                                              \
         input_array1[i] = uniform_random(a, b);                                                                    \
                                                                                                                    \
      ROOT::Double_v x, y;                                                                                          \
      for (int j = 0; j < trials; j += kVS) {                                                                       \
         vecCore::Load<ROOT::Double_v>(x, &input_array1[j]);                                                        \
         y = TMath::tmathfunc(x);                                                                                   \
         vecCore::Store<ROOT::Double_v>(y, &output_array[j]);                                                       \
      }                                                                                                             \
      for (int j = 0; j < trials; j++) {                                                                            \
         Double_t scalar_output = TMath::tmathfunc(input_array1[j]);                                                \
         Double_t vec_output = output_array[j];                                                                     \
         Double_t re =                                                                                              \
            (scalar_output == vec_output && scalar_output == 0) ? 0 : (vec_output - scalar_output) / scalar_output; \
         EXPECT_NEAR(0, re, 1e9 * std::numeric_limits<double>::epsilon());                              \
      }                                                                                                             \
   }

#define TEST_VECTORIZED_TMATH_FUNCTION2(tmathfunc, a, b, c, d)                                                      \
   TEST_F(VectorizedTMathTest, tmathfunc)                                                                           \
   {                                                                                                                \
      int trials = N;                                                                                           \
      for (int i = 0; i < trials; i++) {                                                                            \
         input_array1[i] = uniform_random(a, b);                                                                    \
         input_array2[i] = uniform_random(c, d);                                                                    \
      }                                                                                                             \
      ROOT::Double_v x1, x2, y;                                                                                     \
      for (int j = 0; j < trials; j += kVS) {                                                                       \
         vecCore::Load<ROOT::Double_v>(x1, &input_array1[j]);                                                       \
         vecCore::Load<ROOT::Double_v>(x2, &input_array2[j]);                                                       \
         y = TMath::tmathfunc(x1, x2);                                                                              \
         vecCore::Store<ROOT::Double_v>(y, &output_array[j]);                                                       \
      }                                                                                                             \
      for (int j = 0; j < trials; j++) {                                                                            \
         Double_t scalar_output = TMath::tmathfunc(input_array1[j], input_array2[j]);                               \
         Double_t vec_output = output_array[j];                                                                     \
         Double_t re =                                                                                              \
            (scalar_output == vec_output && scalar_output == 0) ? 0 : (vec_output - scalar_output) / scalar_output; \
         EXPECT_NEAR(0, re, 1e10 * std::numeric_limits<double>::epsilon());                             \
      }                                                                                                             \
   }

#define TEST_VECTORIZED_TMATH_FUNCTION_FLT_RANGE(tmathfunc) TEST_VECTORIZED_TMATH_FUNCTION(tmathfunc, -FLT_MAX, FLT_MAX)

#define TEST_VECTORIZED_TMATH_FUNCTION_FLT_POS_RANGE(tmathfunc) \
   TEST_VECTORIZED_TMATH_FUNCTION(tmathfunc, FLT_MIN, FLT_MAX)

#define TEST_VECTORIZED_TMATH_FUNCTION_FLT_LOG(tmathfunc) TEST_VECTORIZED_TMATH_FUNCTION(tmathfunc, 1, FLT_MAX)

#define TEST_VECTORIZED_TMATH_FUNCTION_FLT_EXP(tmathfunc) \
   TEST_VECTORIZED_TMATH_FUNCTION(tmathfunc, FLT_MIN_10_EXP, FLT_MAX_10_EXP)

#define TEST_VECTORIZED_TMATH_FUNCTION_FLT_EXP_POS(tmathfunc) \
   TEST_VECTORIZED_TMATH_FUNCTION(tmathfunc, FLT_MIN, FLT_MAX_10_EXP)

#define TEST_VECTORIZED_TMATH_FUNCTION_FLT_PROB(tmathfunc) TEST_VECTORIZED_TMATH_FUNCTION(tmathfunc, 0, 1)

#define TEST_VECTORIZED_TMATH_FUNCTION_FLT_BESSEL(tmathfunc) TEST_VECTORIZED_TMATH_FUNCTION(tmathfunc, 0, 1e5)

TEST_VECTORIZED_TMATH_FUNCTION_FLT_LOG(Log2);
TEST_VECTORIZED_TMATH_FUNCTION_FLT_RANGE(BreitWigner);
TEST_VECTORIZED_TMATH_FUNCTION_FLT_RANGE(Gaus);
TEST_VECTORIZED_TMATH_FUNCTION_FLT_EXP(LaplaceDist);
TEST_VECTORIZED_TMATH_FUNCTION_FLT_EXP(LaplaceDistI);

// Freq function has exp(x^2). Hence range -6 to 6.
TEST_VECTORIZED_TMATH_FUNCTION(Freq, -6, 6);

TEST_VECTORIZED_TMATH_FUNCTION_FLT_EXP(BesselI0);
TEST_VECTORIZED_TMATH_FUNCTION(BesselI1, FLT_MIN, FLT_MAX_10_EXP);

TEST_VECTORIZED_TMATH_FUNCTION_FLT_BESSEL(BesselJ0);
TEST_VECTORIZED_TMATH_FUNCTION_FLT_BESSEL(BesselJ1);

int main(int argc, char *argv[])
{
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
