#include "ROOT/RDataFrame.hxx"
#include "ROOT/RForest.hxx"
#include "ROOT/RForestDS.hxx"
#include "ROOT/RForestModel.hxx"
#include "ROOT/RPageStorage.hxx"
#include "ROOT/RPageStorageRoot.hxx"

#include <TClass.h>
#include <TFile.h>
#include <TRandom3.h>

#include "gtest/gtest.h"

#include <exception>
#include <memory>
#include <string>
#include <utility>

using RInputForest = ROOT::Experimental::RInputForest;
using ROutputForest = ROOT::Experimental::ROutputForest;
using RForestModel = ROOT::Experimental::RForestModel;
using RPageSource = ROOT::Experimental::Detail::RPageSource;
using RPageSinkRoot = ROOT::Experimental::Detail::RPageSinkRoot;
using RPageSourceRoot = ROOT::Experimental::Detail::RPageSourceRoot;
using RFieldBase = ROOT::Experimental::Detail::RFieldBase;

TEST(RForest, Basics)
{
   auto model = RForestModel::Create();
   auto fieldPt = model->MakeField<float>("pt");

   //RInputTree tree(model, std::make_unique<RPageSource>("T"));
   //RInputTree tree2(std::make_unique<RPageSource>("T"));
}

TEST(RForest, ReconstructModel)
{
   auto model = RForestModel::Create();
   auto fieldPt = model->MakeField<float>("pt", 42.0);
   auto fieldNnlo = model->MakeField<std::vector<std::vector<float>>>("nnlo");
   auto fieldKlass = model->MakeField<ROOT::Experimental::RForestTest>("klass");
   {
      RPageSinkRoot sinkRoot("myTree", "test.root");
      sinkRoot.Create(model.get());
      sinkRoot.CommitDataset();
   }

   RPageSourceRoot sourceRoot("myTree", "test.root");
   sourceRoot.Attach();

   auto modelReconstructed = sourceRoot.GenerateModel();
   EXPECT_EQ(nullptr, modelReconstructed->GetDefaultEntry()->Get<float>("xyz"));
   auto vecPtr = modelReconstructed->GetDefaultEntry()->Get<std::vector<std::vector<float>>>("nnlo");
   EXPECT_TRUE(vecPtr != nullptr);
   // Don't crash
   vecPtr->push_back(std::vector<float>{1.0});
}

TEST(RForest, StorageRoot)
{
   TFile *file = TFile::Open("test.root", "RECREATE");
   RPageSinkRoot::RSettings settingsWrite;
   settingsWrite.fFile = file;
   RPageSinkRoot sinkRoot("myTree", settingsWrite);

   auto model = RForestModel::Create();
   auto fieldPt = model->MakeField<float>("pt", 42.0);
   auto fieldX = model->MakeField<float>("energy");
   auto fieldStr = model->MakeField<std::string>("string", "abc");

   //auto fieldFail = model->AddField<int>("jets");
   auto fieldJet = model->MakeField<std::vector<float>>("jets" /* TODO(jblomer), {1.0, 2.0}*/);
   auto nnlo = model->MakeField<std::vector<std::vector<float>>>("nnlo");

   sinkRoot.Create(model.get());
   sinkRoot.CommitDataset();
   file->Close();

   file = TFile::Open("test.root", "READ");
   RPageSourceRoot::RSettings settingsRead;
   settingsRead.fFile = file;
   RPageSourceRoot sourceRoot("myTree", settingsRead);
   sourceRoot.Attach();
   file->Close();
}


TEST(RForest, WriteRead)
{
   auto modelWrite = RForestModel::Create();
   auto wrPt = modelWrite->MakeField<float>("pt", 42.0);
   auto wrEnergy = modelWrite->MakeField<float>("energy", 7.0);
   auto wrTag = modelWrite->MakeField<std::string>("tag", "xyz");
   auto wrJets = modelWrite->MakeField<std::vector<float>>("jets");
   wrJets->push_back(1.0);
   wrJets->push_back(2.0);
   auto wrNnlo = modelWrite->MakeField<std::vector<std::vector<float>>>("nnlo");
   wrNnlo->push_back(std::vector<float>());
   wrNnlo->push_back(std::vector<float>{1.0});
   wrNnlo->push_back(std::vector<float>{1.0, 2.0, 4.0, 8.0});
   auto wrKlass = modelWrite->MakeField<ROOT::Experimental::RForestTest>("klass");
   wrKlass->s = "abc";

   auto modelRead = std::unique_ptr<RForestModel>(modelWrite->Clone());

   {
      ROutputForest forest(std::move(modelWrite), std::make_unique<RPageSinkRoot>("f", "test.root"));
      forest.Fill();
   }

   auto rdPt = modelRead->Get<float>("pt");
   auto rdEnergy = modelRead->Get<float>("energy");
   auto rdTag = modelRead->Get<std::string>("tag");
   auto rdJets = modelRead->Get<std::vector<float>>("jets");
   auto rdNnlo = modelRead->Get<std::vector<std::vector<float>>>("nnlo");
   auto rdKlass = modelRead->Get<ROOT::Experimental::RForestTest>("klass");

   RInputForest forest(std::move(modelRead), std::make_unique<RPageSourceRoot>("f", "test.root"));
   EXPECT_EQ(1U, forest.GetNEntries());
   forest.GetEntry(0);

   EXPECT_EQ(42.0, *rdPt);
   EXPECT_EQ(7.0, *rdEnergy);
   EXPECT_STREQ("xyz", rdTag->c_str());

   EXPECT_EQ(2U, rdJets->size());
   EXPECT_EQ(1.0, (*rdJets)[0]);
   EXPECT_EQ(2.0, (*rdJets)[1]);

   EXPECT_EQ(3U, rdNnlo->size());
   EXPECT_EQ(0U, (*rdNnlo)[0].size());
   EXPECT_EQ(1U, (*rdNnlo)[1].size());
   EXPECT_EQ(4U, (*rdNnlo)[2].size());
   EXPECT_EQ(1.0, (*rdNnlo)[1][0]);
   EXPECT_EQ(1.0, (*rdNnlo)[2][0]);
   EXPECT_EQ(2.0, (*rdNnlo)[2][1]);
   EXPECT_EQ(4.0, (*rdNnlo)[2][2]);
   EXPECT_EQ(8.0, (*rdNnlo)[2][3]);

   EXPECT_STREQ("abc", rdKlass->s.c_str());
}

TEST(RForest, View)
{
   auto model = RForestModel::Create();
   auto fieldPt = model->MakeField<float>("pt", 42.0);
   auto fieldTag = model->MakeField<std::string>("tag", "xyz");
   auto fieldJets = model->MakeField<std::vector<float>>("jets");
   fieldJets->push_back(1.0);
   fieldJets->push_back(2.0);

   {
      ROutputForest forest(std::move(model), std::make_unique<RPageSinkRoot>("f", "test.root"));
      forest.Fill();
   }

   RInputForest forest(std::make_unique<RPageSourceRoot>("f", "test.root"));
   auto viewPt = forest.GetView<float>("pt");
   int n = 0;
   for (auto i : forest.GetViewRange()) {
      EXPECT_EQ(42.0, viewPt(i));
      n++;
   }
   EXPECT_EQ(1, n);

   auto viewJets = forest.GetView<std::vector<float>>("jets");
   n = 0;
   for (auto i : forest.GetViewRange()) {
      EXPECT_EQ(2U, viewJets(i).size());
      EXPECT_EQ(1.0, viewJets(i)[0]);
      EXPECT_EQ(2.0, viewJets(i)[1]);
      n++;
   }
   EXPECT_EQ(1, n);
}

TEST(RForest, Capture) {
   auto model = RForestModel::Create();
   float pt;
   model->AddField("pt", &pt);
}

TEST(RForest, Composable)
{
   auto eventModel = RForestModel::Create();
   auto fldPt = eventModel->MakeField<float>("pt", 0.0);

   auto hitModel = RForestModel::Create();
   auto fldHitX = hitModel->MakeField<float>("x", 0.0);
   auto fldHitY = hitModel->MakeField<float>("y", 0.0);

   auto trackModel = RForestModel::Create();
   auto fldTrackEnergy = trackModel->MakeField<float>("energy", 0.0);

   auto fldHits = trackModel->MakeCollection("hits", std::move(hitModel));
   auto fldTracks = eventModel->MakeCollection("tracks", std::move(trackModel));

   {
      auto forest = ROutputForest::Create(std::move(eventModel), "f", "test.root");

      for (unsigned i = 0; i < 8; ++i) {
         for (unsigned t = 0; t < 3; ++t) {
            for (unsigned h = 0; h < 2; ++h) {
               *fldHitX = 4.0;
               *fldHitY = 8.0;
               fldHits->Fill();
            }
            *fldTrackEnergy = 1.0;
            fldTracks->Fill();
         }
         *fldPt = float(i);
         forest->Fill();
      }
   }

   RInputForest forest(std::make_unique<RPageSourceRoot>("f", "test.root"));
   auto viewPt = forest.GetView<float>("pt");
   auto viewTracks = forest.GetViewCollection("tracks");
   auto viewTrackEnergy = viewTracks.GetView<float>("energy");
   auto viewHits = viewTracks.GetViewCollection("hits");
   auto viewHitX = viewHits.GetView<float>("x");
   auto viewHitY = viewHits.GetView<float>("y");

   int nEv = 0;
   for (auto e : forest.GetViewRange()) {
      EXPECT_EQ(float(nEv), viewPt(e));
      EXPECT_EQ(3U, viewTracks(e));

      int nTr = 0;
      for (auto t : viewTracks.GetViewRange(e)) {
         nTr++;
         EXPECT_EQ(1.0, viewTrackEnergy(t));

         EXPECT_EQ(2.0, viewHits(t));
         for (auto h : viewHits.GetViewRange(t)) {
            EXPECT_EQ(4.0, viewHitX(h));
            EXPECT_EQ(8.0, viewHitY(h));
         }
      }
      EXPECT_EQ(3, nTr);

      nEv++;
   }
   EXPECT_EQ(8, nEv);
}

TEST(RForest, TypeName) {
   EXPECT_STREQ("float", ROOT::Experimental::RField<float>::MyTypeName().c_str());
   EXPECT_STREQ("std::vector<std::string>",
                ROOT::Experimental::RField<std::vector<std::string>>::MyTypeName().c_str());
   EXPECT_STREQ("ROOT::Experimental::RForestTest",
                ROOT::Experimental::RField<ROOT::Experimental::RForestTest>::MyTypeName().c_str());
}

namespace {
class RNoDictionary {};
} // namespace

TEST(RForest, TClass) {
   auto modelFail = RForestModel::Create();
   EXPECT_THROW(modelFail->MakeField<RNoDictionary>("nodict"), std::runtime_error);

   auto model = RForestModel::Create();
   auto ptrKlass = model->MakeField<ROOT::Experimental::RForestTest>("klass");

   ROutputForest forest(std::move(model), std::make_unique<RPageSinkRoot>("f", "test.root"));
}


TEST(RForest, RealWorld1)
{
   // See https://github.com/olifre/root-io-bench/blob/master/benchmark.cpp
   auto modelWrite = RForestModel::Create();
   auto& wrEvent   = *modelWrite->MakeField<std::uint32_t>("event");
   auto& wrEnergy  = *modelWrite->MakeField<double>("energy");
   auto& wrTimes   = *modelWrite->MakeField<std::vector<double>>("times");
   auto& wrIndices = *modelWrite->MakeField<std::vector<std::uint32_t>>("indices");

   TRandom3 rnd(42);
   double chksumWrite = 0.0;
   {
      auto forest = ROutputForest::Create(std::move(modelWrite), "f", "test.root");
      constexpr unsigned int nEvents = 60000;
      for (unsigned int i = 0; i < nEvents; ++i) {
         wrEvent = i;
         wrEnergy = rnd.Rndm() * 1000.;

         chksumWrite += double(wrEvent);
         chksumWrite += wrEnergy;

         auto nTimes = 1 + floor(rnd.Rndm() * 1000.);
         wrTimes.resize(nTimes);
         for (unsigned int n = 0; n < nTimes; ++n) {
            wrTimes[n] = 1 + rnd.Rndm()*1000. - 500.;
            chksumWrite += wrTimes[n];
         }

         auto nIndices = 1 + floor(rnd.Rndm() * 1000.);
         wrIndices.resize(nIndices);
         for (unsigned int n = 0; n < nIndices; ++n) {
            wrIndices[n] = 1 + floor(rnd.Rndm() * 1000.);
            chksumWrite += double(wrIndices[n]);
         }

         forest->Fill();
      }
   }

   auto modelRead  = RForestModel::Create();
   auto& rdEvent   = *modelRead->MakeField<std::uint32_t>("event");
   auto& rdEnergy  = *modelRead->MakeField<double>("energy");
   auto& rdTimes   = *modelRead->MakeField<std::vector<double>>("times");
   auto& rdIndices = *modelRead->MakeField<std::vector<std::uint32_t>>("indices");

   double chksumRead = 0.0;
   auto forest = RInputForest::Create(std::move(modelRead), "f", "test.root");
   for (unsigned int i = 0; i < forest->GetNEntries(); ++i) {
      forest->GetEntry(i);
      chksumRead += double(rdEvent) + rdEnergy;
      for (auto t : rdTimes) chksumRead += t;
      for (auto ind : rdIndices) chksumRead += double(ind);
   }

   EXPECT_EQ(chksumRead, chksumWrite);
}


TEST(RForest, RDF)
{
   RInputForest *forest = nullptr;
   auto rdf = std::make_unique<ROOT::RDataFrame>(std::make_unique<ROOT::RDF::RForestDS>(forest));
}
