#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TBasket.h"
#include "TRandom.h"
#include "gtest/gtest.h"


class TTreeTest : public ::testing::Test {
protected:
    virtual void SetUp()
    {
        auto random = new TRandom(836);
        auto file = new TFile("TTreeCompression.root", "RECREATE");
        auto tree = new TTree("tree", "A test tree");
        Double_t data = 0;
        auto branch = tree->Branch("branch", &data);
        for (Int_t ev = 0; ev < 1000; ev++) {
            data = random->Gaus(100, 7);
            tree->Fill();
        }
        file->Write();
        delete random;
        delete tree;
        delete file;
    }
};

TEST_F(TTreeTest, testDefaultCompression)
{
    auto file = new TFile("TTreeCompression.root");
    auto tree = (TTree*)file->Get("tree");
    auto branch = (TBranch*)tree->Get("branch");

    auto compress;
    compress = file->GetCompressionSettings();    
    ASSERT_EQ(compress, 101);

    compress = branch->GetCompressionSettings();
    ASSERT_EQ(compress, 101);

    delete file;
}