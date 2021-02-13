#include "gtest/gtest.h"

#include "THn.h"
#include "TH1.h"
#include "TH2.h"

// Filling THn
TEST(THn, Fill) {
   Int_t bins[2] = {2, 3};
   Double_t xmin[2] = {0., -3.};
   Double_t xmax[2] = {10., 3.};
   THnD hn("hn", "hn", 2, bins, xmin, xmax);

   {
      // Non-overflow
      Double_t x0[2]{4., -0.01};
      EXPECT_EQ(7, hn.GetBin(x0));
      EXPECT_DOUBLE_EQ(0.0, hn.GetBinContent(7));

      EXPECT_EQ(7, hn.Fill(x0, 0.42));
      EXPECT_DOUBLE_EQ(0.42, hn.GetBinContent(7));
   }

   {
      Double_t x0[2]{-0.49, -2.90};
      EXPECT_EQ(1, hn.GetBin(x0));
      EXPECT_DOUBLE_EQ(0., hn.GetBinContent(10));

      EXPECT_EQ(1, hn.Fill(x0, 0.17));
      EXPECT_EQ(1, hn.GetBin(x0));
      EXPECT_DOUBLE_EQ(0.17, hn.GetBinContent(1));
   }

}


TEST(THn, Projection) {
   Int_t bins[2] = {2, 3};
   Double_t xmin[2] = {0., -3.};
   Double_t xmax[2] = {10., 3.};
   THnD hn("hn", "hn", 2, bins, xmin, xmax);

   {
      // Non-overflow
      Double_t x0[2]{4., -0.01};
      hn.Fill(x0, 0.42);
   }

   {
      // Overflow x
      Double_t x0[2]{-0.49, -2.90};
      hn.Fill(x0, 0.17);
   }

   {
      TH1D* hProj = hn.Projection(0);
      EXPECT_EQ(2, hProj->GetNbinsX());
      EXPECT_DOUBLE_EQ(0.17, hProj->GetBinContent(0));
      EXPECT_DOUBLE_EQ(0.42, hProj->GetBinContent(1));
      EXPECT_DOUBLE_EQ(0., hProj->GetBinContent(2));
      EXPECT_DOUBLE_EQ(0., hProj->GetBinContent(3));
      delete hProj;
   }

   {
      TH1D* hProj = hn.Projection(1);
      EXPECT_EQ(3, hProj->GetNbinsX());
      EXPECT_DOUBLE_EQ(0.0, hProj->GetBinContent(0));
      EXPECT_DOUBLE_EQ(0.17, hProj->GetBinContent(1));
      EXPECT_DOUBLE_EQ(0.42, hProj->GetBinContent(2));
      EXPECT_DOUBLE_EQ(0., hProj->GetBinContent(3));
      EXPECT_DOUBLE_EQ(0., hProj->GetBinContent(4));
      delete hProj;
   }

   {
      TH2D* hProj = hn.Projection(0, 1);
      // Yes, x|y are intentionally inverted...
      EXPECT_EQ(3, hProj->GetNbinsX());
      EXPECT_EQ(2, hProj->GetNbinsY());
      EXPECT_DOUBLE_EQ(0., hProj->GetBinContent(0));
      EXPECT_DOUBLE_EQ(.17, hProj->GetBinContent(1));
      EXPECT_DOUBLE_EQ(0., hProj->GetBinContent(2));
      EXPECT_DOUBLE_EQ(0., hProj->GetBinContent(3));
      EXPECT_DOUBLE_EQ(0., hProj->GetBinContent(4));
      EXPECT_DOUBLE_EQ(0., hProj->GetBinContent(5));
      EXPECT_DOUBLE_EQ(0., hProj->GetBinContent(6));
      EXPECT_DOUBLE_EQ(.42, hProj->GetBinContent(7));
      EXPECT_DOUBLE_EQ(0., hProj->GetBinContent(8));
      EXPECT_DOUBLE_EQ(0., hProj->GetBinContent(9));
      EXPECT_DOUBLE_EQ(0., hProj->GetBinContent(10));
      EXPECT_DOUBLE_EQ(0., hProj->GetBinContent(11));
      delete hProj;
   }

   {
      // Overflow x and y
      Double_t x0[2]{-0.49, 3.90};
      hn.Fill(x0, 43.);

      // With axis range, excluding above point; see ROOT-8457
      hn.GetAxis(1)->SetRange(1,2);
      TH1D* hProj = hn.Projection(1);
      EXPECT_EQ(2, hProj->GetNbinsX());
      EXPECT_DOUBLE_EQ(0.0, hProj->GetBinContent(0));
      EXPECT_DOUBLE_EQ(0.17, hProj->GetBinContent(1));
      EXPECT_DOUBLE_EQ(0.42, hProj->GetBinContent(2));
      EXPECT_DOUBLE_EQ(0., hProj->GetBinContent(3));
      EXPECT_DOUBLE_EQ(0., hProj->GetBinContent(4));
      delete hProj;
   }

}
