// @(#)root/test:$Id$
// Authors: David Gonzalez Maline November 2008
//          Martin Storø Nyfløtt  June 2017

#include <sstream>

#include "TH2.h"
#include "TH3.h"
#include "THn.h"
#include "THnSparse.h"

#include "TF1.h"

#include "HFitInterface.h"

#include "TRandom2.h"

#include "gtest/gtest.h"

#include "StressHistogramGlobal.h"

using namespace std;

TEST(StressHistorgram, TestClone1D)
{
   // Tests the clone method for 1D Histograms

   TH1D *h1 = new TH1D("cl1D-h1", "h1-Title", numberOfBins, minRange, maxRange);

   h1->Sumw2();

   for (Int_t e = 0; e < nEvents; ++e) {
      Double_t value = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      h1->Fill(value, 1.0);
   }

   TH1D *h2 = static_cast<TH1D *>(h1->Clone());

   EXPECT_TRUE(HistogramsEquals(h1, h2, cmpOptStats));
   delete h1;
   delete h2;
}

TEST(StressHistorgram, TestClone2D)
{
   // Tests the clone method for 2D Histograms

   TH2D *h1 = new TH2D("cl2D-h1", "h1-Title", numberOfBins, minRange, maxRange, numberOfBins + 2, minRange, maxRange);

   h1->Sumw2();

   for (Int_t e = 0; e < nEvents * nEvents; ++e) {
      Double_t x = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t y = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      h1->Fill(x, y, 1.0);
   }

   TH2D *h2 = static_cast<TH2D *>(h1->Clone());

   EXPECT_TRUE(HistogramsEquals(h1, h2, cmpOptStats));
   delete h1;
   delete h2;
}

TEST(StressHistorgram, TestClone3D)
{
   // Tests the clone method for 3D Histograms

   TH3D *h1 = new TH3D("cl3D-h1", "h1-Title", numberOfBins, minRange, maxRange, numberOfBins + 1, minRange, maxRange,
                       numberOfBins + 2, minRange, maxRange);

   h1->Sumw2();

   for (Int_t e = 0; e < nEvents * nEvents; ++e) {
      Double_t x = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t y = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t z = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      h1->Fill(x, y, z, 1.0);
   }

   TH3D *h2 = static_cast<TH3D *>(h1->Clone());

   EXPECT_TRUE(HistogramsEquals(h1, h2, cmpOptStats));
   delete h1;
   delete h2;
}

TEST(StressHistorgram, TestCloneProfile1D)
{
   // Tests the clone method for 1D Profiles

   TProfile *p1 = new TProfile("cl1D-p1", "p1-Title", numberOfBins, minRange, maxRange);

   for (Int_t e = 0; e < nEvents; ++e) {
      Double_t x = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t y = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      p1->Fill(x, y, 1.0);
   }

   TProfile *p2 = static_cast<TProfile *>(p1->Clone());

   EXPECT_TRUE(HistogramsEquals(p1, p2, cmpOptStats));
   delete p1;
   delete p2;
}

TEST(StressHistorgram, TestCloneProfile2D)
{
   // Tests the clone method for 2D Profiles

   TProfile2D *p1 =
      new TProfile2D("cl2D-p1", "p1-Title", numberOfBins, minRange, maxRange, numberOfBins + 2, minRange, maxRange);

   for (Int_t e = 0; e < nEvents * nEvents; ++e) {
      Double_t x = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t y = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t z = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      p1->Fill(x, y, z, 1.0);
   }

   TProfile2D *p2 = static_cast<TProfile2D *>(p1->Clone());

   EXPECT_TRUE(HistogramsEquals(p1, p2, cmpOptStats));
   delete p1;
   delete p2;
}

TEST(StressHistorgram, TestCloneProfile3D)
{
   // Tests the clone method for 3D Profiles

   TProfile3D *p1 = new TProfile3D("cl3D-p1", "p1-Title", numberOfBins, minRange, maxRange, numberOfBins + 1, minRange,
                                   maxRange, numberOfBins + 2, minRange, maxRange);

   for (Int_t e = 0; e < nEvents * nEvents; ++e) {
      Double_t x = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t y = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t z = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t t = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      p1->Fill(x, y, z, t, 1.0);
   }

   TProfile3D *p2 = static_cast<TProfile3D *>(p1->Clone());

   EXPECT_TRUE(HistogramsEquals(p1, p2));
   delete p1;
   delete p2;
}

TEST(StressHistorgram, TestCloneProfileVar1D)
{
   // Tests the clone method for 1D Profiles with variable bin size

   Double_t v[numberOfBins + 1];
   FillVariableRange(v);

   TProfile *p1 = new TProfile("cl1D-p1", "p1-Title", numberOfBins, v);

   for (Int_t e = 0; e < nEvents; ++e) {
      Double_t x = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t y = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      p1->Fill(x, y, 1.0);
   }

   TProfile *p2 = static_cast<TProfile *>(p1->Clone());

   EXPECT_TRUE(HistogramsEquals(p1, p2, cmpOptStats));
   delete p1;
   delete p2;
}

TEST(StressHistorgram, TestCloneVar1D)
{
   // Tests the clone method for 1D Histograms with variable bin size

   Double_t v[numberOfBins + 1];
   FillVariableRange(v);

   TH1D *h1 = new TH1D("cl1D-h1", "h1-Title", numberOfBins, v);

   h1->Sumw2();

   for (Int_t e = 0; e < nEvents; ++e) {
      Double_t value = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      h1->Fill(value, 1.0);
   }

   TH1D *h2 = static_cast<TH1D *>(h1->Clone());

   EXPECT_TRUE(HistogramsEquals(h1, h2, cmpOptStats));
   delete h1;
   delete h2;
}
