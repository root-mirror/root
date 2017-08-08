// @(#)root/test:$Id$
// Authors: David Gonzalez Maline November 2008
//          Martin Storø Nyfløtt  June 2017

#include <sstream>

#include "TH2.h"
#include "TH3.h"
#include "THn.h"
#include "THnSparse.h"

#include "TProfile3D.h"

#include "TF1.h"

#include "HFitInterface.h"

#include "TRandom2.h"

#include "gtest/gtest.h"

#include "../StressHistogramGlobal.h"

using namespace std;

TEST(StressHistogram, TestAdd3DProfile1)
{
   TH1::SetDefaultSumw2();

   // Tests the second Add method for 3D Profiles

   Double_t c1 = r.Rndm();
   Double_t c2 = r.Rndm();

   TProfile3D p1("t3D1-p1", "p1", numberOfBins, minRange, maxRange, numberOfBins + 1, minRange, maxRange,
                 numberOfBins + 2, minRange, maxRange);
   TProfile3D p2("t3D1-p2", "p2", numberOfBins, minRange, maxRange, numberOfBins + 1, minRange, maxRange,
                 numberOfBins + 2, minRange, maxRange);
   TProfile3D p3("t3D1-p3", "p3=c1*p1+c2*p2", numberOfBins, minRange, maxRange, numberOfBins + 1, minRange, maxRange,
                 numberOfBins + 2, minRange, maxRange);

   for (Int_t e = 0; e < nEvents * nEvents; ++e) {
      Double_t x = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t y = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t z = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t t = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      p1.Fill(x, y, z, t, 1.0);
      p3.Fill(x, y, z, t, c1);
   }

   for (Int_t e = 0; e < nEvents * nEvents; ++e) {
      Double_t x = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t y = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t z = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t t = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      p2.Fill(x, y, z, t, 1.0);
      p3.Fill(x, y, z, t, c2);
   }

   TProfile3D p4("t3D1-p4", "p4=c1*p1+c2*p2", numberOfBins, minRange, maxRange, numberOfBins + 1, minRange, maxRange,
                 numberOfBins + 2, minRange, maxRange);
   p4.Add(&p1, &p2, c1, c2);
   EXPECT_TRUE(HistogramsEquals(p3, p4, cmpOptStats, 1E-10));
}

TEST(StressHistogram, TestAdd3DProfile2)
{
   TH1::SetDefaultSumw2();
   // Tests the second Add method for 3D Profiles

   Double_t c2 = r.Rndm();

   TProfile3D p1("t3D2-p1", "p1", numberOfBins, minRange, maxRange, numberOfBins + 1, minRange, maxRange,
                 numberOfBins + 2, minRange, maxRange);
   TProfile3D p2("t3D2-p2", "p2", numberOfBins, minRange, maxRange, numberOfBins + 1, minRange, maxRange,
                 numberOfBins + 2, minRange, maxRange);
   TProfile3D p3("t3D2-p3", "p3=p1+c2*p2", numberOfBins, minRange, maxRange, numberOfBins + 1, minRange, maxRange,
                 numberOfBins + 2, minRange, maxRange);

   for (Int_t e = 0; e < nEvents * nEvents; ++e) {
      Double_t x = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t y = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t z = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t t = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      p1.Fill(x, y, z, t, 1.0);
      p3.Fill(x, y, z, t, 1.0);
   }

   for (Int_t e = 0; e < nEvents * nEvents; ++e) {
      Double_t x = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t y = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t z = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t t = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      p2.Fill(x, y, z, t, 1.0);
      p3.Fill(x, y, z, t, c2);
   }

   p1.Add(&p2, c2);
   EXPECT_TRUE(HistogramsEquals(p3, p1, cmpOptStats, 1E-10));
}
