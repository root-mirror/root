#include <ROOT/RDataFrame.hxx>
#include <ROOT/RNTupleDS.hxx>

#include <ROOT/RNTuple.hxx>
#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RPageStorage.hxx>

#include <gtest/gtest.h>

using ROOT::Experimental::RNTupleDS;
using ROOT::Experimental::RNTupleWriter;
using ROOT::Experimental::RNTupleModel;
using ROOT::Experimental::Detail::RPageSource;

class RNTupleDSTest : public ::testing::Test {
protected:
   std::string fFileName = "RNTupleDS_test.root";
   std::string fNtplName = "ntuple";
   std::unique_ptr<RPageSource> fPageSource;

   void SetUp() override {
      auto modelWrite = RNTupleModel::Create();
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
      {
         auto ntuple = RNTupleWriter::Recreate(std::move(modelWrite), fNtplName, fFileName);
         ntuple->Fill();
      }
      fPageSource = RPageSource::Create(fNtplName, fFileName);
   }

   void TearDown() override {
      std::remove(fFileName.c_str());
   }
};

TEST_F(RNTupleDSTest, ColTypeNames)
{
   RNTupleDS tds(std::move(fPageSource));

   auto colNames = tds.GetColumnNames();
   ASSERT_EQ(colNames.size(), 5);

   EXPECT_TRUE(tds.HasColumn("pt"));
   EXPECT_TRUE(tds.HasColumn("energy"));
   EXPECT_FALSE(tds.HasColumn("Address"));

   EXPECT_STREQ("std::string", tds.GetTypeName("tag").c_str());
   EXPECT_STREQ("float", tds.GetTypeName("energy").c_str());
}


void ReadTest(const std::string &name, const std::string &fname) {
   auto df = ROOT::Experimental::MakeNTupleDataFrame(name, fname);

   auto count = df.Count();
   auto sumpt = df.Sum<float>("pt");
   auto tag = df.Take<std::string>("tag");
   auto sumjets = df.Sum<std::vector<float>>("jets");
   auto sumvec = [](float red, const std::vector<std::vector<float>> &nnlo) {
      auto sum = 0.f;
      for (auto &v : nnlo)
         for (auto e : v)
            sum += e;
      return red + sum;
   };
   auto sumnnlo = df.Aggregate(sumvec, std::plus<float>{}, "nnlo", 0.f);

   EXPECT_EQ(count.GetValue(), 1ull);
   EXPECT_DOUBLE_EQ(sumpt.GetValue(), 42.f);
   EXPECT_EQ(tag.GetValue().size(), 1ull);
   EXPECT_EQ(tag.GetValue()[0], "xyz");
   EXPECT_EQ(sumjets.GetValue(), 3.f);
   EXPECT_EQ(sumnnlo.GetValue(), 16.f);
}

TEST_F(RNTupleDSTest, Read)
{
   ReadTest(fNtplName, fFileName);
}

struct IMTRAII {
   IMTRAII() { ROOT::EnableImplicitMT(); }
   ~IMTRAII() { ROOT::DisableImplicitMT(); }
};

TEST_F(RNTupleDSTest, ReadMT)
{
   IMTRAII _;

   ReadTest(fNtplName, fFileName);
}
