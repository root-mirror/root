#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooArgusBG.h"
#include "RooAddPdf.h"
#include "RooNLLVar.h"
#include "RooNLLVarNew.h"
#include "RooFitDriver.h"
#include "RooGaussian.h"
#include "RooMinuit.h"
#include "RooMinimizer.h"
#include "RooFitResult.h"

#include <utility>
#include <chrono>

#include "gtest/gtest.h"

struct FitOutput {
    int evalCount;
    long int elapsedTime;
    std::unique_ptr<RooFitResult> result;
};

TEST(testRooFitDriver, SimpleLikelihoodFit)
{
   RooRealVar x("x","x",-10,10);

   RooRealVar mean("mean", "mean",-1, -10, 10);
   RooRealVar width("width", "width", 1., 0.01, 10);
   RooGaussian model("model","model",x,mean,width);

   std::size_t nEvents = 2000;

   RooDataSet *data = model.generate(x, nEvents);

   auto resetGraph = [&](){
       // reset initial fit values so both fits have the same starting condition
       mean.setVal(1);
       width.setVal(2);

       // reset the error such that they are not used as the step size in the fit
       mean.setError(0.0);
       width.setError(0.0);
    };

   auto doFit = [](RooAbsReal& absReal){

      std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

      RooMinimizer minimizer(absReal);
      minimizer.setPrintLevel(-1);
      minimizer.minimize("Minuit", "minuit");
      std::unique_ptr<RooFitResult> result(minimizer.save());

      std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

      auto elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();

      return FitOutput{minimizer.evalCounter(),
                       elapsedTime, std::move(result)};
   };

   resetGraph();

   // do likelihood fitting the old way...
   auto nllScalar = model.createNLL(*data, RooFit::BatchMode(false));
   auto resultScalar = doFit(*nllScalar);
   std::cout << "- scalar mode fit took " << resultScalar.elapsedTime << " ms" << std::endl;

   resetGraph();

   // ...and now the new way with RooFitDriver
   RooNLLVarNew nll("nll", "nll", model);
   RooFitDriver driver("driver", "driver", nll, *data, x);
   auto resultBatchNew = doFit(driver);
   std::cout << "-  batch mode fit took " << resultBatchNew.elapsedTime << " ms" << std::endl;

   // make sure the fit result is the same and that there were the same number
   // of evaluations
   EXPECT_EQ(resultScalar.evalCount, resultBatchNew.evalCount);
   ASSERT_TRUE(resultBatchNew.result->isIdentical(*resultScalar.result));
}
