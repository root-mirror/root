// Tests for the HistFactory
// Authors: Stephan Hageboeck, CERN  01/2019

#include "RooStats/HistFactory/Sample.h"
#include "RooStats/ModelConfig.h"
#include "RooWorkspace.h"
#include "RooArgSet.h"

#include "TROOT.h"
#include "TFile.h"
#include "gtest/gtest.h"

using namespace RooStats;
using namespace RooStats::HistFactory;

TEST(Sample, CopyAssignment)
{
  Sample s("s");
  {
    Sample s1("s1");
    auto hist1 = new TH1D("hist1", "hist1", 10, 0, 10);
    s1.SetHisto(hist1);
    s = s1;
    //Now go out of scope. Should delete hist1, that's owned by s1.
  }
  
  auto hist = s.GetHisto();
  ASSERT_EQ(hist->GetNbinsX(), 10);
}


TEST(HistFactory, Read_ROOT6_16_Model) {
  std::string filename = "./ref_6.16_example_UsingC_channel1_meas_model.root";
  std::unique_ptr<TFile> file(TFile::Open(filename.c_str()));
  if (!file || !file->IsOpen()) {
    filename = TROOT::GetRootSys() + "/roofit/histfactory/test/" + filename;
    file.reset(TFile::Open(filename.c_str()));
  }

  ASSERT_TRUE(file && file->IsOpen());
  RooWorkspace* ws;
  file->GetObject("channel1", ws);
  ASSERT_NE(ws, nullptr);

  auto mc = dynamic_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));
  ASSERT_NE(mc, nullptr);

  RooAbsPdf* pdf = mc->GetPdf();
  ASSERT_NE(pdf, nullptr);

  const RooArgSet* obs = mc->GetObservables();
  ASSERT_NE(obs, nullptr);

  EXPECT_NEAR(pdf->getVal(), 0.17488817, 1.E-8);
  EXPECT_NEAR(pdf->getVal(*obs), 0.95652174, 1.E-8);
}


TEST(HistFactory, Read_ROOT6_16_Combined_Model) {
  std::string filename = "./ref_6.16_example_UsingC_combined_meas_model.root";
  std::unique_ptr<TFile> file(TFile::Open(filename.c_str()));
  if (!file || !file->IsOpen()) {
    filename = TROOT::GetRootSys() + "/roofit/histfactory/test/" + filename;
    file.reset(TFile::Open(filename.c_str()));
  }

  ASSERT_TRUE(file && file->IsOpen());
  RooWorkspace* ws;
  file->GetObject("combined", ws);
  ASSERT_NE(ws, nullptr);

  auto mc = dynamic_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));
  ASSERT_NE(mc, nullptr);

  RooAbsPdf* pdf = mc->GetPdf();
  ASSERT_NE(pdf, nullptr);

  const RooArgSet* obs = mc->GetObservables();
  ASSERT_NE(obs, nullptr);

  EXPECT_NEAR(pdf->getVal(), 0.17488817, 1.E-8);
  EXPECT_NEAR(pdf->getVal(*obs), 0.95652174, 1.E-8);
}
