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

#include "TApplication.h"
#include "TRandom2.h"
#include "TFile.h"
#include "TClass.h"

#include "gtest/gtest.h"

#include "../StressHistogramGlobal.h"

using namespace std;

TEST(StressHistorgram, TestMergeProf3DLabelAll)
{
   // Tests the merge method with fully equally labelled 3D Profiles

   TProfile3D *p1 = new TProfile3D("merge3DLabelAll-p1", "p1-Title", numberOfBins, minRange, maxRange, numberOfBins + 1,
                                   minRange, maxRange, numberOfBins + 2, minRange, maxRange);
   TProfile3D *p2 = new TProfile3D("merge3DLabelAll-p2", "p2-Title", numberOfBins, minRange, maxRange, numberOfBins + 1,
                                   minRange, maxRange, numberOfBins + 2, minRange, maxRange);
   TProfile3D *p3 = new TProfile3D("merge3DLabelAll-p3", "p3-Title", numberOfBins, minRange, maxRange, numberOfBins + 1,
                                   minRange, maxRange, numberOfBins + 2, minRange, maxRange);
   TProfile3D *p4 = new TProfile3D("merge3DLabelAll-p4", "p4-Title", numberOfBins, minRange, maxRange, numberOfBins + 1,
                                   minRange, maxRange, numberOfBins + 2, minRange, maxRange);

   for (Int_t e = 0; e < nEvents * nEvents; ++e) {
      Double_t x = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t y = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t z = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t t = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      p1->Fill(x, y, z, t, 1.0);
      p4->Fill(x, y, z, t, 1.0);
   }

   for (Int_t e = 0; e < nEvents * nEvents; ++e) {
      Double_t x = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t y = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t z = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t t = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      p2->Fill(x, y, z, t, 1.0);
      p4->Fill(x, y, z, t, 1.0);
   }

   for (Int_t e = 0; e < nEvents * nEvents; ++e) {
      Double_t x = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t y = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t z = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t t = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      p3->Fill(x, y, z, t, 1.0);
      p4->Fill(x, y, z, t, 1.0);
   }

   for (Int_t i = 1; i <= numberOfBins; ++i) {
      ostringstream name;
      name << (char)((int)'a' + i - 1);
      p1->GetXaxis()->SetBinLabel(i, name.str().c_str());
      p2->GetXaxis()->SetBinLabel(i, name.str().c_str());
      p3->GetXaxis()->SetBinLabel(i, name.str().c_str());
      p4->GetXaxis()->SetBinLabel(i, name.str().c_str());
   }

   TList *list = new TList;
   list->Add(p2);
   list->Add(p3);

   p1->Merge(list);

   EXPECT_TRUE(HistogramsEquals(p1, p4, cmpOptStats, 1E-10));
   delete p1;
   delete p2;
   delete p3;
   delete p4;
}

TEST(StressHistorgram, TestMergeProf3DLabelAllDiff)
{
   // Tests the merge method with fully differently labelled 3D Profiles

   TProfile3D *p1 = new TProfile3D("merge3DLabelAllDiff-p1", "p1-Title", numberOfBins, minRange, maxRange,
                                   numberOfBins + 1, minRange, maxRange, numberOfBins + 2, minRange, maxRange);
   TProfile3D *p2 = new TProfile3D("merge3DLabelAllDiff-p2", "p2-Title", numberOfBins, minRange, maxRange,
                                   numberOfBins + 1, minRange, maxRange, numberOfBins + 2, minRange, maxRange);
   TProfile3D *p3 = new TProfile3D("merge3DLabelAllDiff-p3", "p3-Title", numberOfBins, minRange, maxRange,
                                   numberOfBins + 1, minRange, maxRange, numberOfBins + 2, minRange, maxRange);
   TProfile3D *p4 = new TProfile3D("merge3DLabelAllDiff-p4", "p4-Title", numberOfBins, minRange, maxRange,
                                   numberOfBins + 1, minRange, maxRange, numberOfBins + 2, minRange, maxRange);

   for (Int_t e = 0; e < nEvents * nEvents; ++e) {
      Double_t x = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t y = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t z = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t t = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      p1->Fill(x, y, z, t, 1.0);
      p4->Fill(x, y, z, t, 1.0);
   }

   for (Int_t e = 0; e < nEvents * nEvents; ++e) {
      Double_t x = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t y = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t z = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t t = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      p2->Fill(x, y, z, t, 1.0);
      p4->Fill(x, y, z, t, 1.0);
   }

   for (Int_t e = 0; e < nEvents * nEvents; ++e) {
      Double_t x = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t y = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t z = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t t = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      p3->Fill(x, y, z, t, 1.0);
      p4->Fill(x, y, z, t, 1.0);
   }

   // It does not work properly! Look, the bins with the same labels
   // are different ones and still the tests passes! This is not
   // consistent with TH1::Merge()
   for (Int_t i = 1; i <= numberOfBins; ++i) {
      ostringstream name;
      name << (char)((int)'a' + i - 1);
      p1->GetXaxis()->SetBinLabel(i, name.str().c_str());
      p1->GetYaxis()->SetBinLabel(i, name.str().c_str());
      p1->GetZaxis()->SetBinLabel(i, name.str().c_str());
      name << 1;
      p2->GetXaxis()->SetBinLabel(i, name.str().c_str());
      p2->GetYaxis()->SetBinLabel(i, name.str().c_str());
      p2->GetZaxis()->SetBinLabel(i, name.str().c_str());
      name << 2;
      p3->GetXaxis()->SetBinLabel(i, name.str().c_str());
      p3->GetYaxis()->SetBinLabel(i, name.str().c_str());
      p3->GetZaxis()->SetBinLabel(i, name.str().c_str());
      name << 3;
      p4->GetXaxis()->SetBinLabel(i, name.str().c_str());
      p4->GetYaxis()->SetBinLabel(i, name.str().c_str());
      p4->GetZaxis()->SetBinLabel(i, name.str().c_str());
   }

   TList *list = new TList;
   list->Add(p2);
   list->Add(p3);

   p1->Merge(list);

   EXPECT_TRUE(HistogramsEquals(p1, p4, cmpOptStats, 1E-10));
   delete p1;
   delete p2;
   delete p3;
   delete p4;
}

TEST(StressHistorgram, TestMergeProf3DLabelDiff)
{
   // Tests the merge with some different labels method for 3D Profiles

   TProfile3D *p1 = new TProfile3D("merge3DLabelDiff-p1", "p1-Title", numberOfBins, minRange, maxRange,
                                   numberOfBins + 1, minRange, maxRange, numberOfBins + 2, minRange, maxRange);
   TProfile3D *p2 = new TProfile3D("merge3DLabelDiff-p2", "p2-Title", numberOfBins, minRange, maxRange,
                                   numberOfBins + 1, minRange, maxRange, numberOfBins + 2, minRange, maxRange);
   TProfile3D *p3 = new TProfile3D("merge3DLabelDiff-p3", "p3-Title", numberOfBins, minRange, maxRange,
                                   numberOfBins + 1, minRange, maxRange, numberOfBins + 2, minRange, maxRange);
   TProfile3D *p4 = new TProfile3D("merge3DLabelDiff-p4", "p4-Title", numberOfBins, minRange, maxRange,
                                   numberOfBins + 1, minRange, maxRange, numberOfBins + 2, minRange, maxRange);

   // It does not work properly! Look, the bins with the same labels
   // are different ones and still the tests passes! This is not
   // consistent with TH1::Merge()
   p1->GetXaxis()->SetBinLabel(2, "gamma");
   p2->GetXaxis()->SetBinLabel(6, "beta");
   p3->GetXaxis()->SetBinLabel(4, "alpha");
   p4->GetXaxis()->SetBinLabel(4, "alpha");

   for (Int_t e = 0; e < nEvents * nEvents; ++e) {
      Double_t x = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t y = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t z = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t t = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      p1->Fill(x, y, z, t, 1.0);
      p4->Fill(x, y, z, t, 1.0);
   }

   for (Int_t e = 0; e < nEvents * nEvents; ++e) {
      Double_t x = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t y = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t z = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t t = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      p2->Fill(x, y, z, t, 1.0);
      p4->Fill(x, y, z, t, 1.0);
   }

   for (Int_t e = 0; e < nEvents * nEvents; ++e) {
      Double_t x = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t y = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t z = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t t = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      p3->Fill(x, y, z, t, 1.0);
      p4->Fill(x, y, z, t, 1.0);
   }

   TList *list = new TList;
   list->Add(p2);
   list->Add(p3);

   p1->Merge(list);

   EXPECT_TRUE(HistogramsEquals(p1, p4, cmpOptStats, 1E-10));
   delete p1;
   delete p2;
   delete p3;
   delete p4;
}

TEST(StressHistorgram, TestMergeProf3DLabelSame)
{
   // Tests the merge with some equal labels method for 3D Profiles

   TProfile3D *p1 = new TProfile3D("merge3DLabelSame-p1", "p1-Title", numberOfBins, minRange, maxRange,
                                   numberOfBins + 1, minRange, maxRange, numberOfBins + 2, minRange, maxRange);
   TProfile3D *p2 = new TProfile3D("merge3DLabelSame-p2", "p2-Title", numberOfBins, minRange, maxRange,
                                   numberOfBins + 1, minRange, maxRange, numberOfBins + 2, minRange, maxRange);
   TProfile3D *p3 = new TProfile3D("merge3DLabelSame-p3", "p3-Title", numberOfBins, minRange, maxRange,
                                   numberOfBins + 1, minRange, maxRange, numberOfBins + 2, minRange, maxRange);
   TProfile3D *p4 = new TProfile3D("merge3DLabelSame-p4", "p4-Title", numberOfBins, minRange, maxRange,
                                   numberOfBins + 1, minRange, maxRange, numberOfBins + 2, minRange, maxRange);

   // It does not work properly! Look, the bins with the same labels
   // are different ones and still the tests passes! This is not
   // consistent with TH1::Merge()
   p1->GetXaxis()->SetBinLabel(4, "alpha");
   p2->GetXaxis()->SetBinLabel(6, "alpha");
   p3->GetXaxis()->SetBinLabel(8, "alpha");
   p4->GetXaxis()->SetBinLabel(4, "alpha");

   for (Int_t e = 0; e < nEvents * nEvents; ++e) {
      Double_t x = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t y = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t z = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t t = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      p1->Fill(x, y, z, t, 1.0);
      p4->Fill(x, y, z, t, 1.0);
   }

   for (Int_t e = 0; e < nEvents * nEvents; ++e) {
      Double_t x = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t y = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t z = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t t = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      p2->Fill(x, y, z, t, 1.0);
      p4->Fill(x, y, z, t, 1.0);
   }

   for (Int_t e = 0; e < nEvents * nEvents; ++e) {
      Double_t x = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t y = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t z = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      Double_t t = r.Uniform(0.9 * minRange, 1.1 * maxRange);
      p3->Fill(x, y, z, t, 1.0);
      p4->Fill(x, y, z, t, 1.0);
   }

   TList *list = new TList;
   list->Add(p2);
   list->Add(p3);

   p1->Merge(list);

   EXPECT_TRUE(HistogramsEquals(p1, p4, cmpOptStats, 1E-10));
   delete p1;
   delete p2;
   delete p3;
   delete p4;
}
