// Tests for the RooSimultaneous
// Authors: Jonas Rembser, CERN  06/2021

#include "RooAddPdf.h"
#include "RooConstVar.h"
#include "RooCategory.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooRealVar.h"
#include "RooSimultaneous.h"
#include "RooProdPdf.h"

#include "gtest/gtest.h"

#include <memory>

/// GitHub issue #8307.
/// A likelihood with a model wrapped in a RooSimultaneous ith one category
/// should give the same results as the likelihood with the model directly.
TEST(RooSimultaneous, ImportFromTreeWithCut)
{
   using namespace RooFit;

   RooRealVar x("x", "x", 0, 10);
   RooRealVar mean("mean", "mean", 1., 0, 10);
   RooRealVar width("width", "width", 1, 0.1, 10);
   RooRealVar nsig("nsig", "nsig", 500, 100, 1000);

   RooGaussian gauss1("gauss1", "gauss1", x, mean, width);
   RooGaussian fconstraint("fconstraint", "fconstraint", mean, RooConst(2.0), RooConst(0.2));

   RooAddPdf model("model", "model", RooArgList(gauss1), RooArgList(nsig));
   RooProdPdf modelConstrained("modelConstrained", "modelConstrained", RooArgSet(model, fconstraint));

   RooCategory cat("cat", "cat");
   cat.defineType("physics");

   RooSimultaneous modelSim("modelSim", "modelSim", RooArgList{modelConstrained}, cat);

   std::unique_ptr<RooDataSet> data{model.generate(x)};
   RooDataSet combData("combData", "combData", x, Index(cat), Import("physics", *data));

   RooArgSet constraints{fconstraint};

   std::unique_ptr<RooAbsReal> nllDirect{modelConstrained.createNLL(combData, Constrain(constraints))};
   std::unique_ptr<RooAbsReal> nllSimWrapped{modelSim.createNLL(combData, Constrain(constraints))};

   EXPECT_FLOAT_EQ(nllDirect->getVal(), nllSimWrapped->getVal());
}
