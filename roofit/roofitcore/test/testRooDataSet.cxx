// Tests for the RooDataSet
// Authors: Stephan Hageboeck, CERN  04/2020

#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooRealVar.h"
#include "RooHelpers.h"
#include "RooCategory.h"

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <RooDataHist.h>
#include <TRandom3.h>
#include <TH1F.h>
#include <TCut.h>
#include <TSystem.h>
#include <ROOT/RMakeUnique.hxx>

#include <TRandom3.h>
#include <TH1F.h>
#include <TCut.h>

#include <fstream>

#include "gtest/gtest.h"

#if !defined(_MSC_VER) || defined(R__ENABLE_BROKEN_WIN_TESTS)
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
#endif

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


/// ROOT-3579 Binned clone seemed to create problems with chains.
/// Code adapted from example in JIRA.
TEST(RooDataSet, BinnedClone) {
  const char* filename[2] = {"RooDataSet_BinnedClone1.root", "RooDataSet_BinnedClone2.root"};
  double sumW = 0;

  for (unsigned int i=0; i < 2; ++i) {
    TFile file(filename[i], "RECREATE");
    TTree tree("cand", "cand");
    double Mes,weight;
    tree.Branch("Mes", &Mes);
    tree.Branch("weight", &weight);

    for (unsigned int j=0; j<20; ++j) {
      Mes = 5.24 + j*0.05/20. + i*0.0003;
      weight = 1.3 + j + i;
      sumW += weight;
      tree.Fill();
    }
    file.WriteObject(&tree, "cand");
    file.Close();
  }

  TChain chain("cand");
  chain.Add(filename[0]);
  chain.Add(filename[1]);
  RooRealVar mes("Mes","Mes",5.28,5.24,5.29);
  mes.setBin(40);
  RooRealVar weight("weight","weight",1,0,100);
  RooDataSet* data = new RooDataSet("dataset", "dataset", &chain, RooArgSet(mes,weight), 0, weight.GetName());
  RooDataHist* hist = data->binnedClone();

  EXPECT_DOUBLE_EQ(hist->sumEntries(), sumW);

  delete hist;//invalid read here
  delete data;//and here too

  gSystem->Unlink(filename[0]);
  gSystem->Unlink(filename[1]);
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

/// ROOT-10845 IsOnHeap() always returned false.
TEST(RooDataSet, IsOnHeap) {
  auto setp = new RooDataSet();
  EXPECT_TRUE(setp->IsOnHeap());

  RooDataSet setStack;
  EXPECT_FALSE(setStack.IsOnHeap());
}

/// ROOT-10935. Cannot read a category from a text file if states are given as index instead of label.
TEST(RooDataSet, ReadCategory) {
  RooCategory cat("cat", "cat", {{"One", 1}, {"Two", 2}, {"Three", 3}});
  cat.setRange("OneTwo", "One,Two");

  RooRealVar x("x", "x", 0., 10.);

  constexpr auto filename = "datasetWithCategory.txt";
  std::ofstream file(filename);
  file << "1. One\n" << "2. Two\n" << "3. 3" << std::endl;

  auto dataset = RooDataSet::read(filename, RooArgList(x,cat));
  EXPECT_EQ(dataset->numEntries(), 3);
  for (int i=1; i < 4; ++i) {
    EXPECT_EQ(static_cast<RooRealVar*>(dataset->get(i-1)->find("x"))->getVal(), i);
    EXPECT_EQ(static_cast<RooCategory*>(dataset->get(i-1)->find("cat"))->getIndex(), i);
  }

  gSystem->Unlink(filename);
}


/// ROOT-8173. Reading negative exponents from file goes wrong.
TEST(RooDataSet, ReadNegativeExponent) {
  RooRealVar x("x", "x", 0., 10.);

  constexpr auto filename = "datasetWithCategory.txt";
  std::ofstream file(filename);
  file << "2.E-1\n" << "2.E-2\n" << "2.E-3" << std::endl;

  auto dataset = RooDataSet::read(filename, RooArgList(x));
  ASSERT_EQ(dataset->numEntries(), 3);
  const double solutions[] = { 2.E-1, 2.E-2, 2.E-3 };
  for (int i=0; i < 3; ++i) {
    EXPECT_EQ(static_cast<RooRealVar*>(dataset->get(i)->find("x"))->getVal(), solutions[i]);
  }

  gSystem->Unlink(filename);
}

/// root-project/root#6408: Importing from tree deletes the TFile with the original
TEST(RooDataSet, CrashAfterImportFromTree) {
  TTree* tree = new TTree("tree", "tree");
  double var = 1;
  tree->Branch("var", &var, "var/D");
  var = 1;
  tree->Fill();
  var = 2;
  tree->Fill();

  auto roovar = std::make_unique<RooRealVar>("var", "var", 0, 10);
  auto output_file = std::make_unique<TFile>("test.root", "RECREATE", "output_file");

  ASSERT_TRUE(output_file->IsOpen());
  auto data_set = std::make_unique<RooDataSet>("data_set", "data_set", tree, RooArgSet(*roovar));

  // Would crash, since the TFile would be deleted by importing:
  ASSERT_TRUE(output_file->IsOpen());

  EXPECT_EQ(data_set->sumEntries(), 2.);
  EXPECT_EQ(data_set->numEntries(), 2);
  EXPECT_EQ(static_cast<RooRealVar*>(data_set->get(0)->find("var"))->getVal(), 1.);
  EXPECT_EQ(static_cast<RooRealVar*>(data_set->get(1)->find("var"))->getVal(), 2.);

}
