#include <gtest/gtest.h>

#include <TMVA/RStandardScaler.hxx>
#include <TMVA/RTensorUtils.hxx>

#include <cmath>

using namespace TMVA::Experimental;

static const std::string filename_ = "http://root.cern/files/tmva_class_example.root";
static const std::vector<std::string> variables = {"var1", "var2", "var3", "var4"};

TEST(RStandardScaler, TensorOutputShape)
{
   ROOT::RDataFrame df("TreeS", filename_);
   auto x = AsTensor<float>(df, variables);

   RStandardScaler<float> scaler;
   scaler.Fit(x);
   auto y = scaler.Compute(x);

   EXPECT_EQ(y.GetShape().size(), 2ul);
   EXPECT_EQ(x.GetShape()[0], y.GetShape()[0]);
   EXPECT_EQ(x.GetShape()[1], y.GetShape()[1]);
}

TEST(RStandardScaler, VectorOutputShape)
{
   ROOT::RDataFrame df("TreeS", filename_);
   auto x = AsTensor<float>(df, variables);

   RStandardScaler<float> scaler;
   scaler.Fit(x);
   auto y = scaler.Compute({1.0, 2.0, 3.0, 4.0});

   EXPECT_EQ(x.GetShape()[1], y.size());
}

TEST(RStandardScaler, GetMeans)
{
   ROOT::RDataFrame df("TreeS", filename_);
   auto x = AsTensor<float>(df, variables);

   RStandardScaler<float> scaler;
   scaler.Fit(x);
   const auto y = scaler.GetMeans();

   for (std::size_t j = 0; j < x.GetShape()[1]; j++) {
      float m = 0.0;
      for (std::size_t i = 0; i < x.GetShape()[0]; i++) {
         m += x(i, j);
      }
      m /= x.GetShape()[0];
      ASSERT_NEAR(y[j], m, 1e-4);
   }
}

TEST(RStandardScaler, ComputeMeans)
{
   ROOT::RDataFrame df("TreeS", filename_);
   auto x = AsTensor<float>(df, variables);

   RStandardScaler<float> scaler;
   scaler.Fit(x);
   const auto y = scaler.Compute(x);

   for (std::size_t j = 0; j < x.GetShape()[1]; j++) {
      float m = 0.0;
      for (std::size_t i = 0; i < x.GetShape()[0]; i++) {
         m += y(i, j);
      }
      m /= x.GetShape()[0];
      ASSERT_NEAR(0.0, m, 1e-4);
   }
}

TEST(RStandardScaler, ComputeStds)
{
   ROOT::RDataFrame df("TreeS", filename_);
   auto x = AsTensor<float>(df, variables);

   RStandardScaler<float> scaler;
   scaler.Fit(x);
   const auto y = scaler.Compute(x);

   for (std::size_t j = 0; j < x.GetShape()[1]; j++) {
      float s = 0.0;
      for (std::size_t i = 0; i < x.GetShape()[0]; i++) {
         s += y(i, j) * y(i, j);
      }
      s = std::sqrt(s / (x.GetShape()[0] - 1));
      ASSERT_NEAR(1.0, s, 1e-4);
   }
}

TEST(RStandardScaler, CompareVectorTensorOutput)
{
   ROOT::RDataFrame df("TreeS", filename_);
   auto x = AsTensor<float>(df, variables);

   RStandardScaler<float> scaler;
   scaler.Fit(x);
   const auto y = scaler.Compute(x);

   for (std::size_t i = 0; i < x.GetShape()[0]; i++) {
      const auto y2 = scaler.Compute({x(i, 0), x(i, 1), x(i, 2), x(i, 3)});
      for (std::size_t j = 0; j < x.GetShape()[1]; j++) {
          EXPECT_EQ(y2[j], y(i, j));
      }
   }
}

TEST(RStandardScaler, SaveLoad)
{
   ROOT::RDataFrame df("TreeS", filename_);
   auto x = AsTensor<float>(df, variables);

   RStandardScaler<float> scaler;
   scaler.Fit(x);
   scaler.Save("foo", "RStandardScalerSaveLoad.root");

   RStandardScaler<float> scaler2("foo", "RStandardScalerSaveLoad.root");
   auto y1 = scaler.Compute({1.0, 2.0, 3.0, 4.0});
   auto y2 = scaler2.Compute({1.0, 2.0, 3.0, 4.0});
   for(std::size_t i = 0; i < 4; i++) {
      EXPECT_EQ(y1[i], y2[i]);
   }
}
