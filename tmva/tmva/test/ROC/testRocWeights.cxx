// @(#)root/tmva $Id$
// Author: Kim Albertsson

/*************************************************************************
 * Copyright (C) 2017, Kim Albertsson                                    *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////
//  Collection of tests to verify that roc curves and integrals are //
//  calculated correctly.                                           //
//////////////////////////////////////////////////////////////////////

#include "gtest/gtest.h"

#include "TMVA/ROCCurve.h"

#include <algorithm>
#include <random>
#include <vector>
#include <functional>

namespace TMVA {

std::default_random_engine generator;

class ROCCurveTest : public ::testing::Test {
   using fvec_t = std::vector<Float_t>;
   using gen_t = void (*)(fvec_t &, fvec_t &, fvec_t &, fvec_t &, size_t);

protected:
   ROCCurveTest() {}
   virtual ~ROCCurveTest() {}

   /**
    * Generates two variables, the signal class (A) will be uniformly
    * distributed while the background class (B) will be triangular.
    *     __________
    *     |\_      |
    *     |  \_ A  |
    *     | B  \_  |
    *     |______\_|
    *
    * The analytical roc curve has in this case an area of 2/3.
    */
   static void gen_ut(fvec_t &a, fvec_t &b, fvec_t &aw, fvec_t &bw, size_t num_samples)
   {
      std::uniform_real_distribution<Float_t> distribution(0., 1.);

      a.reserve(num_samples);
      b.reserve(num_samples);
      aw.reserve(num_samples);
      bw.reserve(num_samples);

      for (size_t i = 0; i < num_samples; ++i) {
         a.push_back(distribution(generator));
         b.push_back(distribution(generator));
         aw.push_back(1);
         bw.push_back(distribution(generator));
      }

      std::sort(std::begin(b), std::end(b));
      std::sort(std::begin(bw), std::end(bw), std::greater<Float_t>());
   }

   /**
    * Generates two uniformly distributed variables. When doing a
    * classification, you can do no better than random -> AUC = 0.5.
    */
   static void gen_uu(fvec_t &a, fvec_t &b, fvec_t &aw, fvec_t &bw, size_t num_samples)
   {
      std::uniform_real_distribution<Float_t> distribution(0., 1.);

      a.reserve(num_samples);
      b.reserve(num_samples);
      aw.reserve(num_samples);
      bw.reserve(num_samples);

      for (size_t i = 0; i < num_samples; ++i) {
         a.push_back(distribution(generator));
         b.push_back(distribution(generator));
         aw.push_back(1);
         bw.push_back(1);
      }
   }

   /**
    * Generates four fixed samples regardless of num_samples input.
    */
   static void gen_4samples(fvec_t &a, fvec_t &b, fvec_t &aw, fvec_t &bw, size_t)
   {
      a.push_back(0.5);
      a.push_back(1.0);
      b.push_back(0.0);
      b.push_back(0.5);
      aw.push_back(0.5);
      aw.push_back(1.0);
      bw.push_back(1.0);
      bw.push_back(4.0);
   }

   /**
    * Generates random data according to datagen_function and calculates the
    * resluting AUC score.
    */
   Float_t singleAuc(size_t num_samples, gen_t datagen_function)
   {
      fvec_t a;
      fvec_t b;
      fvec_t aw;
      fvec_t bw;

      datagen_function(a, b, aw, bw, num_samples);

      TMVA::ROCCurve roc = TMVA::ROCCurve(a, b, aw, bw);
      return roc.GetROCIntegral();
   }

   /**
    * Averages the AUC score from several runs.
    */
   Float_t avgAuc(size_t num_samples, size_t N, gen_t datagen_function)
   {
      Float_t sum = 0.;
      for (size_t i = 0; i < N; ++i) {
         const Float_t integral = this->singleAuc(num_samples, datagen_function);
         // std::cout << "AUC: " << integral << std::endl;
         sum += integral;
      }
      return sum / (Float_t)N;
   }
};

TEST_F(ROCCurveTest, aucSimple)
{
   // Simple sanity check to make sure weights are respected.
   // Simple in the sense it only uses 4 non-random datapoints.
   EXPECT_NEAR(singleAuc(0, gen_4samples), 0.866666666667, 0.005);
}

TEST_F(ROCCurveTest, aucRandom)
{
   // Larger sanity check. Two uniform distributions
   // with equal weights should yield AUC of 0.5.
   EXPECT_NEAR(avgAuc(10000, 10, gen_uu), 1. / 2., 0.005);
}

TEST_F(ROCCurveTest, aucRandomWithWeights)
{
   // Larger sanity check. Using uniform dist as signal
   // and triangular as background should yield AUC of 2/3.
   //
   // Using 0.002 as the limit should yield a probability
   // of an error of 0.0001%. (Estimated numerically).
   // If the error is triggered, consider changing the seed of "generator" and/or
   // increase either the error limit, the number of datapoints or the number of
   // averagings.
   EXPECT_NEAR(avgAuc(10000, 10, gen_ut), 2. / 3., 0.005);
}

int main(int argc, char *argv[])
{
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}

}