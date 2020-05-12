// Tests for the RooDataSet
// Authors: Stephan Hageboeck, CERN  04/2020

#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooRealVar.h"
#include "RooHelpers.h"
#include "TTree.h"

#include <TRandom3.h>
#include <TH1F.h>
#include <TCut.h>

#include "gtest/gtest.h"

/// ROOT-10676
/// The RooDataSet warns that it's not using all variables if the selection string doesn't
/// make use of all variables. Although true, the user has no way to suppress this.
TEST(RooDataSet, ImportFromTreeWithCut)
{
  RooHelpers::HijackMessageStream hijack(RooFit::INFO, RooFit::InputArguments);

  TTree tree("tree", "tree");
  double thex, they;
  tree.Branch("x", &thex);
  tree.Branch("y", &they);
  tree.Branch("z", &they);
  thex = -0.337;
  they = 1.;
  tree.Fill();

  thex = 0.337;
  they = 1.;
  tree.Fill();

  thex = 1.337;
  they = 1.;
  tree.Fill();

  RooRealVar x("x", "x", 0);
  RooRealVar y("y", "y", 0);
  RooRealVar z("z", "z", 0);
  RooDataSet data("data", "data", &tree, RooArgSet(x, y, z), "x>y");

  EXPECT_TRUE(hijack.str().empty()) << "Messages issued were: " << hijack.str();
  EXPECT_EQ(data.numEntries(), 1);

  RooRealVar* theX = dynamic_cast<RooRealVar*>(data.get(0)->find("x"));
  ASSERT_NE(theX, nullptr);
  EXPECT_FLOAT_EQ(theX->getVal(), 1.337);
}


/// ROOT-9528 Branch names are capped after a certain number of characters
TEST(RooDataSet, ImportLongBranchNames) {

  TTree tree("theTree", "theTree");
  double doub = 0.;
  tree.Branch("HLT_mu6_mu4_bBmumux_BsmumuPhi_delayed_L1BPH_2M8_MU6MU4_BPH_0DR15_MU6MU4", &doub);
  doub = 2.;
  tree.Fill();
  doub = 4.;
  tree.Fill();

  RooRealVar *v = new RooRealVar("HLT_mu6_mu4_bBmumux_BsmumuPhi_delayed_L1BPH_2M8_MU6MU4_BPH_0DR15_MU6MU4", "HLT_mu6_mu4_bBmumux_BsmumuPhi_delayed_L1BPH_2M8_MU6MU4_BPH_0DR15_MU6MU4", 0., -100000., 100000.);

  RooDataSet ds("ds", "ds", RooArgSet(*v), RooFit::Import(tree));
  EXPECT_EQ(static_cast<RooRealVar*>(ds.get(0)->find(*v))->getVal(), 2.);
  EXPECT_EQ(static_cast<RooRealVar*>(ds.get(1)->find(*v))->getVal(), 4.);
  EXPECT_EQ(ds.numEntries(), 2);
  EXPECT_DOUBLE_EQ(ds.sumEntries("HLT_mu6_mu4_bBmumux_BsmumuPhi_delayed_L1BPH_2M8_MU6MU4_BPH_0DR15_MU6MU4 > 3."), 1.);
}



/// ROOT-4580, possibly solved by ROOT-10517
TEST(RooDataSet, ReducingData) {
  //Test Data hist and such.
  TTree mytree("tree","tree") ;
  double mass_x, track0_chi2_x, track1_chi2_x;

  mytree.Branch("track0_chi2",&track0_chi2_x,"track0_chi2/D");
  mytree.Branch("track1_chi2",&track1_chi2_x,"track1_chi2/D");
  mytree.Branch("mass",&mass_x,"mass/D");
  for (int i=0; i < 50; i++) {
    track0_chi2_x = gRandom->Landau(1,0.5);
    track1_chi2_x = gRandom->Landau(1,0.5);
    mass_x = gRandom->Gaus(20, 0.5);
    mytree.Fill() ;
  }

  Double_t chi2cutval = 1.0;
  constexpr Double_t massmin = 0;
  constexpr Double_t massmax = 40;

  //Now use roofit
  //observables from ttree
  RooRealVar mymass("mass", "mass", massmin, massmax);
  RooRealVar track0_chi2("track0_chi2", "track0_chi2", -10., 90);
  RooRealVar track1_chi2("track1_chi2", "track1_chi2", -10., 90);

  //get the datasets
  RooDataSet *data_unbinned = new RooDataSet("mass_example","mass example", &mytree, RooArgSet(mymass,track0_chi2,track1_chi2));
  std::unique_ptr<RooDataHist> data( data_unbinned->binnedClone("data") );

  for (int i=0; i<3; ++i){
    // Check with root:
    TH1F test_hist(Form("h%i", i), "histo", 10, massmin, massmax);
    chi2cutval+=0.5;

    TCut chi2_test_cut = Form("max(track0_chi2,track1_chi2)<%f",chi2cutval);

    Long64_t drawnEvents = mytree.Draw(Form("mass>>h%i", i), chi2_test_cut /*&& mass_cut*/);
    ASSERT_NE(drawnEvents, 0l);
    ASSERT_EQ(test_hist.Integral(), drawnEvents);

    // For unbinned data, reducing should be equivalent to the tree.
    std::unique_ptr<RooDataSet> data_unbinned_reduced ( static_cast<RooDataSet*>(data_unbinned->reduce(RooFit::Cut(chi2_test_cut))) );
    EXPECT_DOUBLE_EQ(data_unbinned_reduced->sumEntries(), test_hist.Integral());
    EXPECT_EQ(data_unbinned_reduced->numEntries(), test_hist.Integral());

    // When using binned data, reducing and expecting the ame number of entries as in the unbinned case is not possible,
    // since information is lost if entries to the left and right of the cut end up in the same bin.
    // Therefore, can only test <=
    std::unique_ptr<RooDataHist> reduced_binned_data ( static_cast<RooDataHist*>(data->reduce(RooFit::Cut(chi2_test_cut))) );
    if (floor(chi2cutval) == chi2cutval)
      EXPECT_FLOAT_EQ(reduced_binned_data->sumEntries(), test_hist.Integral());
    else
      EXPECT_LE(reduced_binned_data->sumEntries(), test_hist.Integral());
  }
}

