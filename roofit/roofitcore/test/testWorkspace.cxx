// Tests for the RooWorkspace
// Authors: Stephan Hageboeck, CERN  01/2019
#include "RooWorkspace.h"
#include "RooGlobalFunc.h"
#include "RooHelpers.h"
#include "RooGaussian.h"
#include "RooArgList.h"
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooStats/ModelConfig.h"

#include "TFile.h"
#include "TSystem.h"

#include "gtest/gtest.h"

using namespace RooStats;

/// ROOT-9777, cloning a RooWorkspace. The ModelConfig did not get updated
/// when a workspace was cloned, and was hence pointing to a non-existing workspace.
///
TEST(RooWorkspace, CloneModelConfig_ROOT_9777)
{
   const char* filename = "ROOT-9777.root";

   RooRealVar x("x", "x", 1, 0, 10);
   RooRealVar mu("mu", "mu", 1, 0, 10);
   RooRealVar sigma("sigma", "sigma", 1, 0, 10);

   RooGaussian pdf("Gauss", "Gauss", x, mu, sigma);

   {
      TFile outfile(filename, "RECREATE");
      
      // now create the model config for this problem
      RooWorkspace* w = new RooWorkspace("ws");
      ModelConfig modelConfig("ModelConfig", w);
      modelConfig.SetPdf(pdf);
      modelConfig.SetParametersOfInterest(RooArgSet(sigma));
      modelConfig.SetGlobalObservables(RooArgSet(mu));
      w->import(modelConfig);

      outfile.WriteObject(w, "ws");
      delete w;
   }
   
   RooWorkspace *w2;
   {
      TFile infile(filename, "READ");
      RooWorkspace *w;
      infile.GetObject("ws", w);
      ASSERT_TRUE(w) << "Workspace not read from file.";

      w2 = new RooWorkspace(*w);
      delete w;
   }
   
   w2->Print();

   ModelConfig *mc = dynamic_cast<ModelConfig*>(w2->genobj("ModelConfig"));
   ASSERT_TRUE(mc) << "ModelConfig not retrieved.";
   mc->Print();

   ASSERT_TRUE(mc->GetGlobalObservables()) << "GlobalObsevables in mc broken.";
   mc->GetGlobalObservables()->Print();

   ASSERT_TRUE(mc->GetParametersOfInterest()) << "ParametersOfInterest in mc broken.";
   mc->GetParametersOfInterest()->Print();

   gSystem->Unlink(filename);
}



/// Set up a simple workspace for later tests.
class TestRooWorkspaceWithGaussian : public ::testing::Test {
protected:
  TestRooWorkspaceWithGaussian() :
  Test()
  {
    RooRealVar x("x", "x", 1, 0, 10);
    RooRealVar mu("mu", "mu", 1, 0, 10);
    RooRealVar sigma("sigma", "sigma", 1, 0, 10);

    RooGaussian pdf("Gauss", "Gauss", x, mu, sigma);

    TFile outfile(_filename, "RECREATE");

    // now create the model config for this problem
    RooWorkspace w("ws");
    RooStats::ModelConfig modelConfig("ModelConfig", &w);
    modelConfig.SetPdf(pdf);
    modelConfig.SetParametersOfInterest(RooArgSet(sigma));
    modelConfig.SetGlobalObservables(RooArgSet(mu));
    w.import(modelConfig);

    outfile.WriteObject(&w, "ws");
  }

  ~TestRooWorkspaceWithGaussian() {
    gSystem->Unlink(_filename);
  }

  const char* _filename = "ROOT-9777.root";
};


/// Test the string tokeniser that does all the string splitting for the RooWorkspace
/// implementation.
TEST(RooHelpers, Tokeniser)
{
  std::vector<std::string> tok = RooHelpers::tokenise("abc, def, ghi", ", ");
  EXPECT_EQ(tok.size(), 3U);
  EXPECT_EQ(tok[0], "abc");
  EXPECT_EQ(tok[1], "def");
  EXPECT_EQ(tok[2], "ghi");

  std::vector<std::string> tok2 = RooHelpers::tokenise("abc, def", ":");
  EXPECT_EQ(tok2.size(), 1U);
  EXPECT_EQ(tok2[0], "abc, def");

  std::vector<std::string> tok3 = RooHelpers::tokenise(",  ,abc, def,", ", ");
  EXPECT_EQ(tok3.size(), 2U);
  EXPECT_EQ(tok3[0], "abc");
  EXPECT_EQ(tok3[1], "def");

  std::vector<std::string> tok4 = RooHelpers::tokenise(",  ,abc, def,", ",");
  EXPECT_EQ(tok4.size(), 3U);
  EXPECT_EQ(tok4[0], "  ");
  EXPECT_EQ(tok4[1], "abc");
  EXPECT_EQ(tok4[2], " def");


}

/// Test proper string handling when importing an object from a workspace
/// in a different file.
TEST_F(TestRooWorkspaceWithGaussian, ImportFromFile)
{
  std::ostringstream spec;
  spec << _filename << ":" << "ws:Gauss";

  RooWorkspace w("ws");

  //Expect successful import:
  EXPECT_FALSE(w.import(spec.str().c_str()));

#if !defined(_MSC_VER) || defined(R__ENABLE_BROKEN_WIN_TESTS)
  //Expect import failures:
  RooHelpers::HijackMessageStream hijack(RooFit::ERROR, RooFit::InputArguments, "ws");
  EXPECT_TRUE(w.import("bogus:abc"));
  EXPECT_FALSE(hijack.str().empty());

  hijack.stream().str("");
  ASSERT_TRUE(hijack.str().empty());
  EXPECT_TRUE(w.import( (spec.str()+"bogus").c_str()));
  EXPECT_FALSE(hijack.str().empty());
#endif
}

/// [ROOT-7921] When using EDIT, cannot build PDFs from edit PDF.
TEST_F(TestRooWorkspaceWithGaussian, RooCustomiserInterface) {
  TFile file(_filename, "READ");
  RooWorkspace* ws;
  file.GetObject("ws", ws);
  ASSERT_NE(ws, nullptr);

  // Prepare
  ASSERT_NE(ws->factory("SUM:sum(a[0.5,0,1]*Gauss,Gauss)"), nullptr);
  ASSERT_NE(ws->factory("expr:sig2(\"1 + @0 * @1\", {sigma_alpha[0.1], theta_alpha[0, -5, 5]})"), nullptr);
  ASSERT_NE(ws->factory("EDIT::editPdf(sum, sigma=sig2)"), nullptr);
  ASSERT_NE(ws->factory("Gaussian::constraint_alpha(global_alpha[0], theta_alpha, 1)"), nullptr);

  // Build a product using the edited pdf. This failed because of ROOT-7921
  // Problem was in RooCustomizer::CustIFace::create
  EXPECT_NE(ws->factory("PROD::model_constrained(editPdf, constraint_alpha)"), nullptr);

  // Test the other code path in RooCustomizer::CustIFace::create.
  // Edit the top-level pdf in-place, replacing all existing conflicting nodes in the workspace by <node>_orig
  ASSERT_NE(ws->factory("EDIT::model_constrained(model_constrained, mu=mu2[-1,-10,10])"), nullptr);

  // Test that the new model_constrained has been altered
  auto model_constrained = ws->pdf("model_constrained");
  ASSERT_NE(model_constrained, nullptr);
  EXPECT_TRUE(model_constrained->dependsOn(*ws->var("mu2")));
  EXPECT_FALSE(model_constrained->dependsOn(*ws->var("mu")));

  // Test that the old model still exists suffixed with _orig
  auto model_constrained_orig = ws->pdf("model_constrained_orig");
  ASSERT_NE(model_constrained_orig, nullptr);
  EXPECT_TRUE(model_constrained_orig->dependsOn(*ws->var("mu")));
  EXPECT_FALSE(model_constrained_orig->dependsOn(*ws->var("mu2")));
  EXPECT_NE(ws->pdf("Gauss_editPdf_orig"), nullptr);
}


/// Test that things still work when hash lookup for elements
/// is performed.
TEST_F(TestRooWorkspaceWithGaussian, HashLookupInWorkspace) {
  TFile file(_filename, "READ");
  RooWorkspace* ws;
  file.GetObject("ws", ws);
  ASSERT_NE(ws, nullptr);

  ws->useFindsWithHashLookup(true);

  // Prepare
  ASSERT_NE(ws->factory("SUM:sum(a[0.5,0,1]*Gauss,Gauss)"), nullptr);
  ASSERT_NE(ws->factory("expr:sig2(\"1 + @0 * @1\", {sigma_alpha[0.1], theta_alpha[0, -5, 5]})"), nullptr);
  ASSERT_NE(ws->factory("EDIT::editPdf(sum, sigma=sig2)"), nullptr);
  ASSERT_NE(ws->factory("Gaussian::constraint_alpha(global_alpha[0], theta_alpha, 1)"), nullptr);

  // Build a product using the edited pdf. This failed because of ROOT-7921
  // Problem was in RooCustomizer::CustIFace::create
  EXPECT_NE(ws->factory("PROD::model_constrained(editPdf, constraint_alpha)"), nullptr);

  // Test the other code path in RooCustomizer::CustIFace::create.
  // Edit the top-level pdf in-place, replacing all existing conflicting nodes in the workspace by <node>_orig
  ASSERT_NE(ws->factory("EDIT::model_constrained(model_constrained, mu=mu2[-1,-10,10])"), nullptr);

  // Test that the new model_constrained has been altered
  auto model_constrained = ws->pdf("model_constrained");
  ASSERT_NE(model_constrained, nullptr);
  EXPECT_TRUE(model_constrained->dependsOn(*ws->var("mu2")));
  EXPECT_FALSE(model_constrained->dependsOn(*ws->var("mu")));

  // Test that the old model still exists suffixed with _orig
  auto model_constrained_orig = ws->pdf("model_constrained_orig");
  ASSERT_NE(model_constrained_orig, nullptr);
  EXPECT_TRUE(model_constrained_orig->dependsOn(*ws->var("mu")));
  EXPECT_FALSE(model_constrained_orig->dependsOn(*ws->var("mu2")));
  EXPECT_NE(ws->pdf("Gauss_editPdf_orig"), nullptr);
}
