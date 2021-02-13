#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"

#include "gtest/gtest.h"

TEST(TTreeTruncatedDatatypes, float16double32leaves)
{
   const auto ofileName = "float16double32leaves.root";
   TFile f(ofileName, "RECREATE");
   auto t = new TTree("t", "t");

   Float16_t  f1 = 0., f2 = 0.;
   Double32_t d1 = 0., d2 = 0.;
   t->Branch("f1", &f1, "f1/f");
   t->Branch("f2", &f2, "f2/f[0,1000000,20]");
   t->Branch("d1", &d1, "d1/d");
   t->Branch("d2", &d2, "d2/d[0,1000000,20]");
   t->Fill();

   f1 = f2 = 1.;
   d1 = d2 = 1.;
   t->Fill();

   f1 = f2 = 999999.;
   d1 = d2 = 999999.;
   t->Fill();

   t->Write();
   delete t;

   t = (TTree*) f.Get("t");
   t->SetBranchAddress("f1", &f1);
   t->SetBranchAddress("f2", &f2);
   t->SetBranchAddress("d1", &d1);
   t->SetBranchAddress("d2", &d2);

   t->GetEntry(0);
   EXPECT_FLOAT_EQ (0., f1);
   EXPECT_FLOAT_EQ (0., f2);
   EXPECT_DOUBLE_EQ(0., d1);
   EXPECT_DOUBLE_EQ(0., d2);

   t->GetEntry(1);
   EXPECT_FLOAT_EQ (1.,               f1);
   EXPECT_FLOAT_EQ (0.95367432,       f2); // Here the precision kicks in
   EXPECT_DOUBLE_EQ(1.,               d1);
   EXPECT_DOUBLE_EQ(0.95367431640625, d2); // Here the precision kicks in

   t->GetEntry(2);
   EXPECT_FLOAT_EQ (999936.,            f1);  // Here the precision kicks in
   EXPECT_FLOAT_EQ (999999.,            f2);
   EXPECT_DOUBLE_EQ(999999.,            d1);
   EXPECT_DOUBLE_EQ(999999.04632568359, d2);

   delete t;
   f.Close();
   gSystem->Unlink(ofileName);
}
