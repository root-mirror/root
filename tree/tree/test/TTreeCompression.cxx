#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include "gtest/gtest.h"


class TTreeTest : public ::testing::Test {
protected:
    virtual void SetUp()
    {
        auto random = new TRandom(836);
        auto file = new TFile("TTree.root", "RECREATE");
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
    auto file = new TFile("TTree.root");
    auto compress = file->GetCompressionSettings();
    
    ASSERT_EQ(compress, 101);

    delete file;
}