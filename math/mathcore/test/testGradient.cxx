// @(#)root/test:$Id$
// Author: Alejandro García Montoro 07/2017

#include "Fit/BinData.h"
#include "ROOT/EExecutionPolicy.hxx"
#include "Fit/FitUtil.h"
#include "Fit/UnBinData.h"
#include "Math/WrappedMultiTF1.h"

#include "Fit/Chi2FCN.h"
#include "Fit/Fitter.h"
#include "Fit/LogLikelihoodFCN.h"
#include "Fit/PoissonLikelihoodFCN.h"
#include "HFitInterface.h"
#include "TF1.h"
#include "TH1.h"
#include "TF2.h"
#include "TH2.h"
#include "TROOT.h"

#include "TRandom3.h"

#include "gtest/gtest.h"

#include <chrono>
#include <iostream>
#include <string>

int printLevel = 0;

// Class to encapsulate the types that define how the gradient test is
// performed; it also stores information strings about the types.
//    DataType defines how to instantiate the gradient evaluation: Double_t,
//    Double_v.
//    ExecutionPolicyType defines the execution policy: kSerial, kMultithread...
//    DataInfoStr points to a human-readable string describing DataType (e.g.,
//    "Scalar", "Vectorial")
//    PolicyInfoStr points to a human-readable string describing
//    ExecutionPolicyType (e.g., "Serial", "Multithread")
template <typename U, ROOT::EExecutionPolicy V, int W, const char *dataInfoStr, const char *policyInfoStr>
struct GradientTestTraits {
   using DataType = U;
   static constexpr ROOT::EExecutionPolicy ExecutionPolicyType() { return V; };
   static constexpr int Dimensions() { return W; };

   static void PrintTypeInfo(const std::string &fittingInfo)
   {
      std::cout << "---------------- TEST INFO ----------------" << std::endl;
      std::cout << "- Fitting type:     " << fittingInfo << std::endl;
      std::cout << "- Data type:        " << dataInfoStr << std::endl;
      std::cout << "- Execution policy: " << policyInfoStr << std::endl;
      std::cout << "- Dimensions:       " << Dimensions() << std::endl;
      std::cout << "-------------------------------------------" << std::endl;
   }
};

// Info strings describing data types
char scalarStr[] = "Scalar";
char vectorStr[] = "Vectorial";

// Info strings describing execution policies
char serialStr[] = "Serial";
char mthreadStr[] = "Multithread";

// Typedefs of GradientTestTraits for scalar (serial and multithreaded)
// scenarios
using ScalarSerial1D = GradientTestTraits<Double_t, ROOT::EExecutionPolicy::kSequential, 1, scalarStr, serialStr>;
using ScalarMultithread1D =
   GradientTestTraits<Double_t, ROOT::EExecutionPolicy::kMultiThread, 1, scalarStr, mthreadStr>;
using ScalarSerial2D = GradientTestTraits<Double_t, ROOT::EExecutionPolicy::kSequential, 2, scalarStr, serialStr>;
using ScalarMultithread2D =
   GradientTestTraits<Double_t, ROOT::EExecutionPolicy::kMultiThread, 2, scalarStr, mthreadStr>;

#ifdef R__HAS_VECCORE

// Typedefs of GradientTestTraits for vectorial (serial and multithreaded)
// scenarios
using VectorialSerial1D =
   GradientTestTraits<ROOT::Double_v, ROOT::EExecutionPolicy::kSequential, 1, vectorStr, serialStr>;
using VectorialMultithread1D =
   GradientTestTraits<ROOT::Double_v, ROOT::EExecutionPolicy::kMultiThread, 1, vectorStr, mthreadStr>;
using VectorialSerial2D =
   GradientTestTraits<ROOT::Double_v, ROOT::EExecutionPolicy::kSequential, 2, vectorStr, serialStr>;
using VectorialMultithread2D =
   GradientTestTraits<ROOT::Double_v, ROOT::EExecutionPolicy::kMultiThread, 2, vectorStr, mthreadStr>;

#endif

// Interface abstracting the model function and its related data (number of parameters, parameters, the model function
// itself...)
template <class T>
struct Model {
   virtual void FillModelData(ROOT::Fit::BinData *&data) = 0;
   virtual void FillModelData(ROOT::Fit::UnBinData *&data) = 0;

   TF1 *fModelFunction;
   unsigned fNumParams;
   const Double_t *fParams;
};

// Class storing a 1D model function, its parameters and a histogram with data from the function
template <class T>
struct Model1D : public Model<T> {
   Model1D()
   {
      Model<T>::fNumParams = fNumParams;
      Model<T>::fParams = fParams;

      // Create TF1 from model function and initialize the fit function
      std::stringstream streamTF1;
      streamTF1 << "f" << this;
      std::string nameTF1 = streamTF1.str();

      Model<T>::fModelFunction = new TF1(nameTF1.c_str(), Function, 100, 200, 4);
      Model<T>::fModelFunction->SetParameters(Model<T>::fParams);

      // Assure the to-be-created histogram does not replace an old one
      std::stringstream streamTH1;
      streamTH1 << "h" << this;
      std::string nameTH1 = streamTH1.str();

      auto oldTH1 = gROOT->FindObject(nameTH1.c_str());
      if (oldTH1)
         delete oldTH1;

      // Create TH1 and fill it with values from model function
      fNumPoints = 12801;
      //fNumPoints = 11;
      fHistogram = new TH1D(nameTH1.c_str(), "Test random numbers", fNumPoints, 100, 200);
      gRandom->SetSeed(1);
      fHistogram->FillRandom(nameTF1.c_str(), 1000000);
   }

   static T Function(const T *data, const double *params)
   {
      return params[0] * exp(-(*data + (-130.)) * (*data + (-130.)) / 2) +
             params[1] * exp(-(params[2] * (*data * (0.01)) - params[3] * ((*data) * (0.01)) * ((*data) * (0.01))));
   }

   // Fill binned data
   void FillModelData(ROOT::Fit::BinData *&data) override
   {
      data = new ROOT::Fit::BinData(fNumPoints, 1);
      ROOT::Fit::FillData(*data, fHistogram, Model<T>::fModelFunction);
   }

   // Fill unbinned data
   void FillModelData(ROOT::Fit::UnBinData *&data) override
   {
      data = new ROOT::Fit::UnBinData(fNumPoints, 1);

      TAxis *coords = fHistogram->GetXaxis();
      for (unsigned i = 0; i < fNumPoints; i++)
         data->Add(coords->GetBinCenter(i));
   }

   TH1D *fHistogram;
   static const unsigned fNumParams = 4;
   const Double_t fParams[fNumParams] = {1, 1000, 7.5, 1.5};
   unsigned fNumPoints;
};

// Class storing a DD model function, its parameters and a histogram with data from the function
template <class T>
struct Model2D : public Model<T> {
   Model2D()
   {
      Model<T>::fNumParams = fNumParams;
      Model<T>::fParams = fParams;

      // Create TF1 from model function and initialize the fit function
      std::stringstream streamTF2;
      streamTF2 << "f" << this;
      std::string nameTF2 = streamTF2.str();

      Model<T>::fModelFunction = new TF2(nameTF2.c_str(), Function, -5, 5, -5, 5, fNumParams);
      Model<T>::fModelFunction->SetParameters(Model<T>::fParams);

      // Assure the to-be-created histogram does not replace an old one
      std::stringstream streamTH2;
      streamTH2 << "h" << this;
      std::string nameTH2 = streamTH2.str();

      auto oldTH2 = gROOT->FindObject(nameTH2.c_str());
      if (oldTH2)
         delete oldTH2;

      // Create TH1 and fill it with values from model function
      fNumPoints = 801;
      fHistogram = new TH2D(nameTH2.c_str(), "Test random numbers", fNumPoints, -5, 5, fNumPoints, -5, 5);
      gRandom->SetSeed(1);
      fHistogram->FillRandom(nameTF2.c_str(), 1000);
   }

#ifdef R__HAS_VECCORE
   static T TemplatedGaus(T x, Double_t mean, Double_t sigma, Bool_t norm = false)
   {
      if (sigma == 0)
         return 1.e30;

      T arg = (x - mean) / sigma;

      // for |arg| > 39  result is zero in double precision
      vecCore::Mask_v<T> mask = !(arg < -39.0 || arg > 39.0);

      // Initialize the result to 0.0
      T res(0.0);

      // Compute the function only when the arg meets the criteria, using the mask computed before
      vecCore::MaskedAssign<T>(res, mask, vecCore::math::Exp(-0.5 * arg * arg));

      if (!norm)
         return res;

      return res / (2.50662827463100024 * sigma); // sqrt(2*Pi)=2.50662827463100024
   }

   static T Function(const T *data, const double *params)
   {
      return params[0] * TemplatedGaus(data[0], params[1], params[2]) * TemplatedGaus(data[1], params[3], params[4]);
   }

#else

   static T Function(const T *data, const double *params)
   {
      return params[0] * TMath::Gaus(data[0], params[1], params[2]) * TMath::Gaus(data[1], params[3], params[4]);
   }
#endif

   // Fill binned data
   void FillModelData(ROOT::Fit::BinData *&data) override
   {
      data = new ROOT::Fit::BinData(fNumPoints, 2);
      ROOT::Fit::FillData(*data, fHistogram, Model<T>::fModelFunction);
   }

   // Fill unbinned data
   void FillModelData(ROOT::Fit::UnBinData *&data) override
   {
      data = new ROOT::Fit::UnBinData(fNumPoints, 2);

      TAxis *xCoords = fHistogram->GetXaxis();
      TAxis *yCoords = fHistogram->GetYaxis();
      for (unsigned i = 0; i < fNumPoints; i++)
         data->Add(xCoords->GetBinCenter(i), yCoords->GetBinCenter(i));
   }

   TH2D *fHistogram;
   static const unsigned fNumParams = 5;
   const Double_t fParams[fNumParams] = {300, 0., 2., 0., 3.};
   unsigned fNumPoints;
};

// Helper class used to encapsulate the calls to the gradient interfaces,
// templated with a GradientTestTraits type
template <class T, class U>
struct GradientTestEvaluation {
   // Types to instantiate Chi2FCN: the gradient function interface handles the
   // parameters, so its base typedef
   // has to be always a double; the parametric function interface is templated
   // to test both serial and vectorial
   // cases.
   using GradFunctionType = ROOT::Math::IGradientFunctionMultiDimTempl<double>;
   using BaseFunctionType = ROOT::Math::IParamMultiFunctionTempl<typename T::DataType>;

   // Basic type to compare against
   using ScalarSerial =
      GradientTestTraits<double, ROOT::EExecutionPolicy::kSequential, T::Dimensions(), scalarStr, serialStr>;

   GradientTestEvaluation()
   {
      if (T::Dimensions() == 1)
         fModel = new Model1D<typename T::DataType>();
      else if (T::Dimensions() == 2)
         fModel = new Model2D<typename T::DataType>();

      fNumParams = fModel->fNumParams;
      fFitFunction = new ROOT::Math::WrappedMultiTF1Templ<typename T::DataType>(*(fModel->fModelFunction),
                                                                                fModel->fModelFunction->GetNdim());

      // Fill (binned or unbinned) data
      fModel->FillModelData(fData);
   }

   virtual void SetFitter() = 0;

   Double_t BenchmarkSolution(Double_t *solution)
   {
      std::chrono::time_point<std::chrono::system_clock> start, end;

      start = std::chrono::system_clock::now();
      for (int i = 0; i < fNumRepetitions; i++)
         fFitter->Gradient(fModel->fParams, solution);
      end = std::chrono::system_clock::now();

      std::chrono::duration<Double_t> timeElapsed = end - start;

      if (printLevel > 0) {
         std::cout << "Gradient is : " << fFitter->NDim() << "  ";
         for (unsigned int i = 0; i < fNumParams; ++i)
            std::cout << "  " << solution[i];
         std::cout << std::endl;
         std::cout << "elapsed time is " << timeElapsed.count() << std::endl;
      }

      return timeElapsed.count() / fNumRepetitions;

   }

   static const int fNumRepetitions = 2;

   unsigned fNumParams;
   Model<typename T::DataType> *fModel;
   ROOT::Math::WrappedMultiTF1Templ<typename T::DataType> *fFitFunction;
   U *fData;
   ROOT::Fit::BasicFCN<GradFunctionType, BaseFunctionType, U> *fFitter;
};

template <class T>
struct Chi2GradientTestEvaluation : public GradientTestEvaluation<T, ROOT::Fit::BinData> {
   using typename GradientTestEvaluation<T, ROOT::Fit::BinData>::GradFunctionType;
   using typename GradientTestEvaluation<T, ROOT::Fit::BinData>::BaseFunctionType;
   using typename GradientTestEvaluation<T, ROOT::Fit::BinData>::ScalarSerial;

   Chi2GradientTestEvaluation() { SetFitter(); }

   virtual void SetFitter() override
   {
      this->fFitter = new ROOT::Fit::Chi2FCN<GradFunctionType, BaseFunctionType>(*(this->fData), *(this->fFitFunction),
                                                                                 T::ExecutionPolicyType());
   }
};

template <class T>
struct PoissonLikelihoodGradientTestEvaluation : public GradientTestEvaluation<T, ROOT::Fit::BinData> {
   using typename GradientTestEvaluation<T, ROOT::Fit::BinData>::GradFunctionType;
   using typename GradientTestEvaluation<T, ROOT::Fit::BinData>::BaseFunctionType;
   using typename GradientTestEvaluation<T, ROOT::Fit::BinData>::ScalarSerial;

   PoissonLikelihoodGradientTestEvaluation() { SetFitter(); }

   virtual void SetFitter() override
   {
      this->fFitter = new ROOT::Fit::PoissonLikelihoodFCN<GradFunctionType, BaseFunctionType>(
         *(this->fData), *(this->fFitFunction), 0, true, T::ExecutionPolicyType());
   }
};

template <class T>
struct LogLikelihoodGradientTestEvaluation : public GradientTestEvaluation<T, ROOT::Fit::UnBinData> {
   using typename GradientTestEvaluation<T, ROOT::Fit::UnBinData>::GradFunctionType;
   using typename GradientTestEvaluation<T, ROOT::Fit::UnBinData>::BaseFunctionType;
   using typename GradientTestEvaluation<T, ROOT::Fit::UnBinData>::ScalarSerial;

   LogLikelihoodGradientTestEvaluation() { SetFitter(); }

   virtual void SetFitter() override
   {
      this->fFitter = new ROOT::Fit::LogLikelihoodFCN<GradFunctionType, BaseFunctionType>(
         *(this->fData), *(this->fFitFunction), 0, false, T::ExecutionPolicyType());
   }
};

// Test class: creates a reference solution (computing the gradient with scalar
// values in a serial scenario), and
// compares its values and its performance against the evaluation of the
// gradient specified by the GradientTestTraits
// type.
template <class T>
class Chi2GradientTest : public ::testing::Test, public Chi2GradientTestEvaluation<T> {
   using typename Chi2GradientTestEvaluation<T>::ScalarSerial;

protected:
   virtual void SetUp()
   {
      T::PrintTypeInfo("Chi2FCN");

      // Set up the scalar serial reference fitter, to benchmark times.
      Chi2GradientTestEvaluation<ScalarSerial> reference;

      fReferenceSolution.resize(this->fNumParams);
      fReferenceTime = reference.BenchmarkSolution(fReferenceSolution.data());
   }

   Double_t fReferenceTime;
   std::vector<Double_t> fReferenceSolution;
};

// Test class: creates a reference solution (computing the gradient with scalar
// values in a serial scenario), and
// compares its values and its performance against the evaluation of the
// gradient specified by the GradientTestTraits
// type.
template <class T>
class PoissonLikelihoodGradientTest : public ::testing::Test, public PoissonLikelihoodGradientTestEvaluation<T> {
   using typename PoissonLikelihoodGradientTestEvaluation<T>::ScalarSerial;

protected:
   virtual void SetUp()
   {
      T::PrintTypeInfo("PoissonLikelihoodFCN");

      // Set up the scalar serial reference fitter, to benchmark times.
      PoissonLikelihoodGradientTestEvaluation<ScalarSerial> reference;

      fReferenceSolution.resize(this->fNumParams);
      fReferenceTime = reference.BenchmarkSolution(fReferenceSolution.data());
   }

   Double_t fReferenceTime;
   std::vector<Double_t> fReferenceSolution;
};

// Test class: creates a reference solution (computing the gradient with scalar
// values in a serial scenario), and
// compares its values and its performance against the evaluation of the
// gradient specified by the GradientTestTraits
// type.
template <class T>
class LogLikelihoodGradientTest : public ::testing::Test, public LogLikelihoodGradientTestEvaluation<T> {
   using typename LogLikelihoodGradientTestEvaluation<T>::ScalarSerial;

protected:
   virtual void SetUp()
   {
      T::PrintTypeInfo("LogLikelihoodFCN");

      // Set up the scalar serial reference fitter, to benchmark times.
      LogLikelihoodGradientTestEvaluation<ScalarSerial> reference;

      fReferenceSolution.resize(this->fNumParams);
      fReferenceTime = reference.BenchmarkSolution(fReferenceSolution.data());
   }

   Double_t fReferenceTime;
   std::vector<Double_t> fReferenceSolution;
};

// Types used by Google Test to instantiate the tests.
#ifdef R__HAS_VECCORE
#  ifdef R__USE_IMT
typedef ::testing::Types<ScalarMultithread1D, VectorialSerial1D, VectorialMultithread1D, ScalarMultithread2D,
                         VectorialSerial2D, VectorialMultithread2D>
   TestTypes;
#  else
typedef ::testing::Types<VectorialSerial1D, VectorialSerial2D> TestTypes;
#  endif
#else
#  ifdef R__USE_IMT
typedef ::testing::Types<ScalarMultithread1D, ScalarMultithread2D> TestTypes;
#  else
typedef ::testing::Types<ScalarSerial1D,ScalarSerial2D> TestTypes;
#  endif
#endif


TYPED_TEST_SUITE(LogLikelihoodGradientTest, TestTypes);

// Test EvalChi2Gradient and outputs its speedup against the scalar serial case.
TYPED_TEST(LogLikelihoodGradientTest, LogLikelihoodGradient)
{
   std::vector<Double_t> solution(TestFixture::fNumParams);

   Double_t benchmarkTime = TestFixture::BenchmarkSolution(solution.data());

   std::cout << std::fixed << std::setprecision(4);
   std::cout << "Speed-up with respect to scalar serial case: " << TestFixture::fReferenceTime / benchmarkTime;
   std::cout << std::endl;

   for (unsigned i = 0; i < TestFixture::fNumParams; i++) {
      EXPECT_NEAR(solution[i], TestFixture::fReferenceSolution[i], 1e-6);
   }
}


TYPED_TEST_SUITE(Chi2GradientTest, TestTypes);

// Test EvalChi2Gradient and outputs its speedup against the scalar serial case.
TYPED_TEST(Chi2GradientTest, Chi2Gradient)
{
   std::vector<Double_t> solution(TestFixture::fNumParams);

   Double_t benchmarkTime = TestFixture::BenchmarkSolution(solution.data());

   std::cout << std::fixed << std::setprecision(4);
   std::cout << "Speed-up with respect to scalar serial case: " << TestFixture::fReferenceTime / benchmarkTime;
   std::cout << std::endl;

   for (unsigned i = 0; i < TestFixture::fNumParams; i++) {
      EXPECT_NEAR(solution[i], TestFixture::fReferenceSolution[i], 1e-6);
   }
}

TYPED_TEST_SUITE(PoissonLikelihoodGradientTest, TestTypes);

// Test EvalChi2Gradient and outputs its speedup against the scalar serial case.
TYPED_TEST(PoissonLikelihoodGradientTest, PoissonLikelihoodGradient)
{
   std::vector<Double_t> solution(TestFixture::fNumParams);

   Double_t benchmarkTime = TestFixture::BenchmarkSolution(solution.data());

   std::cout << std::fixed << std::setprecision(4);
   std::cout << "Speed-up with respect to scalar serial case: " << TestFixture::fReferenceTime / benchmarkTime;
   std::cout << std::endl;

   for (unsigned i = 0; i < TestFixture::fNumParams; i++) {
      EXPECT_NEAR(solution[i], TestFixture::fReferenceSolution[i], 1e-6);
   }
}

// add main() to avoid a linking error
int main(int argc, char **argv)
{

   // Parse command line arguments
   for (Int_t i = 1; i < argc; i++) {
      std::string arg = argv[i];
      if (arg == "-v") {
         std::cout << "---running in verbose mode" << std::endl;
         printLevel = 1;
      } else if (arg == "-vv") {
         std::cout << "---running in very verbose mode" << std::endl;
         printLevel = 2;
      } else if (arg == "-vvv") {
         std::cout << "---running in very very verbose mode" << std::endl;
         printLevel = 3;
      }
   }

   // This allows the user to override the flag on the command line.
   ::testing::InitGoogleTest(&argc, argv);

   return RUN_ALL_TESTS();
}
