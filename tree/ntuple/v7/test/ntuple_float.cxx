#include <ROOT/RNTuple.hxx>
#include <ROOT/RNTupleModel.hxx>
#include <TFile.h>

#include "gtest/gtest.h"

#include <cmath>
#include <limits>
#include <string>
#include <vector>

using RNTupleReader = ROOT::Experimental::RNTupleReader;
using RNTupleWriter = ROOT::Experimental::RNTupleWriter;
using RNTupleModel = ROOT::Experimental::RNTupleModel;
using RFieldBase = ROOT::Experimental::Detail::RFieldBase;
using float24_t = ROOT::Experimental::float24_t;
using float16_t = ROOT::Experimental::float16_t;
using float8_t = ROOT::Experimental::float8_t;
using ENTupleInfo = ROOT::Experimental::ENTupleInfo;

template <class T>
using RField = ROOT::Experimental::RField<T>;

namespace {

/**
 * An RAII wrapper around an open temporary file on disk. It cleans up the guarded file when the wrapper object
 * goes out of scope.
 */
class FileRaii {
private:
   std::string fPath;
public:
   explicit FileRaii(const std::string &path) : fPath(path) { }
   FileRaii(const FileRaii&) = delete;
   FileRaii& operator=(const FileRaii&) = delete;
   ~FileRaii() { std::remove(fPath.c_str()); }
   std::string GetPath() const { return fPath; }
};

} // anonymous namespace

TEST(RNTupleFloat, float24_t)
{
   FileRaii fileGuard("test_ntuple_float_float24_t.root");
   const std::string_view ntupleName{"float24_NTuple"};
   {
      auto model = RNTupleModel::Create();
      auto fieldPt = model->MakeField<float24_t>("ft24");
      auto fieldVec = model->MakeField<std::vector<float24_t>>("ft24vec");
      auto ntuple = RNTupleWriter::Recreate(std::move(model), ntupleName, fileGuard.GetPath());
      for (int i = 0; i < 100; ++i) {
         *fieldPt = static_cast<float>(i + 0.1f * i);
         *fieldVec = { *fieldPt, *fieldPt };
         ntuple->Fill();
      }
      for (int i = 0; i < 100; ++i) {
         *fieldPt = static_cast<float>(i * 1200);
         *fieldVec = { *fieldPt, *fieldPt, *fieldPt };
         ntuple->Fill();
      }
      *fieldPt = std::numeric_limits<float>::infinity() * (-1.0f); // negative infinity
      *fieldVec = { *fieldPt, *fieldPt };
      ntuple->Fill();
      *fieldPt = std::numeric_limits<float>::infinity(); // positive infinity
      *fieldVec = { *fieldPt, *fieldPt };
      ntuple->Fill();
      *fieldPt = static_cast<float>(std::nanf("1")); // NaN (Not A Number)
      *fieldVec = { *fieldPt };
      ntuple->Fill();
   }
   auto model = RNTupleModel::Create();
   auto fieldPt = model->MakeField<float24_t>("ft24");
   auto fieldVec = model->MakeField<std::vector<float24_t>>("ft24vec");
   auto ntuple = RNTupleReader::Open(std::move(model), ntupleName, fileGuard.GetPath());
   auto float24View = ntuple->GetView<float24_t>("ft24");
   auto vec24View = ntuple->GetView<std::vector<float24_t>>("ft24vec");
   for (int i = 0; i < 100; ++i) {
      EXPECT_NEAR(static_cast<float>(i + 0.1*i), float24View(i), 0.001); // EXPECT_FLOAT_EQ will fail here, so use EXPECT_NEAR instead.
      EXPECT_NEAR(static_cast<float>(i + 0.1*i), vec24View(i).at(1), 0.001);
      EXPECT_EQ(static_cast<float>(i*1200.0f), float24View(i+100));
      EXPECT_EQ(static_cast<float>(i*1200.0f), vec24View(i+100).at(2));
      ntuple->LoadEntry(i);
      EXPECT_NEAR(static_cast<float>(i + 0.1*i), *fieldPt, 0.001);
      EXPECT_NEAR(static_cast<float>(i + 0.1*i), (*fieldVec).at(1), 0.001);
      ntuple->LoadEntry(i+100);
      EXPECT_EQ(static_cast<float>(i*1200.0f), *fieldPt);
      EXPECT_EQ(static_cast<float>(i*1200.0f), (*fieldVec).at(2));
   }
   EXPECT_EQ(true, std::isinf(float24View(200)));
   EXPECT_EQ(true, std::isinf(vec24View(200).at(1)));
   EXPECT_EQ(true, std::isinf(float24View(201)));
   EXPECT_EQ(true, std::isinf(vec24View(201).at(1)));
   EXPECT_EQ(true, std::isnan(float24View(202)));
   EXPECT_EQ(true, std::isnan(vec24View(202).at(0)));
   std::ostringstream os;
   // Only check if calling PrintInfo does not lead to error.
   ntuple->PrintInfo(ENTupleInfo::kSummary, os);
   ntuple->PrintInfo(ENTupleInfo::kStorageDetails, os);
   auto autoGenerateModel = RNTupleReader::Open(ntupleName, fileGuard.GetPath());
}

TEST(RNTupleFloat, float16_t)
{
   FileRaii fileGuard("test_ntuple_float_float16_t.root");
   const std::string_view ntupleName{"float16_t NTuple"};
   {
      auto model = RNTupleModel::Create();
      auto fieldPt = model->MakeField<float16_t>("ft16");
      auto fieldVec = model->MakeField<std::array<float16_t, 2>>("ft16vec");
      auto ntuple = RNTupleWriter::Recreate(std::move(model), ntupleName, fileGuard.GetPath());
      for (int i = 0; i < 100; ++i) {
         *fieldPt = static_cast<float>(0.1 * i);
         *fieldVec = { *fieldPt, *fieldPt };
         ntuple->Fill();
      }
      for (int i = 0; i < 100; ++i) {
         *fieldPt = static_cast<float>(i * 12);
         *fieldVec = { *fieldPt, *fieldPt };
         ntuple->Fill();
      }
      *fieldPt = std::numeric_limits<float>::infinity() * (-1.0f); // negative infinity
      *fieldVec = { *fieldPt, *fieldPt };
      ntuple->Fill();
      *fieldPt = std::numeric_limits<float>::infinity(); // positive infinity
      *fieldVec = { *fieldPt, *fieldPt };
      ntuple->Fill();
      *fieldPt = static_cast<float>(std::nanf("1")); // NaN (Not A Number)
      *fieldVec = { *fieldPt, *fieldPt };
      ntuple->Fill();
   }
   auto model = RNTupleModel::Create();
   auto fieldPt = model->MakeField<float16_t>("ft16");
   auto fieldVec = model->MakeField<std::array<float16_t, 2>>("ft16vec");
   auto ntuple = RNTupleReader::Open(std::move(model), ntupleName, fileGuard.GetPath());
   auto float16View = ntuple->GetView<float16_t>("ft16");
   auto vec16View = ntuple->GetView<std::array<float16_t, 2>>("ft16vec");
   for (int i = 0; i < 100; ++i) {
      EXPECT_NEAR(0.1*i, float16View(i), 0.01); // EXPECT_FLOAT_EQ will fail here, so use EXPECT_NEAR instead.
      EXPECT_NEAR(static_cast<float>(0.1*i), vec16View(i).at(1), 0.01);
      EXPECT_EQ(static_cast<float>(i*12.0f), float16View(i+100));
      EXPECT_EQ(static_cast<float>(i*12.0f), vec16View(i+100).at(1));
      ntuple->LoadEntry(i);
      EXPECT_NEAR(static_cast<float>(0.1*i), *fieldPt, 0.01);
      EXPECT_NEAR(static_cast<float>(0.1*i), (*fieldVec).at(1), 0.01);
      ntuple->LoadEntry(i+100);
      EXPECT_EQ(static_cast<float>(i*12.0f), *fieldPt);
      EXPECT_EQ(static_cast<float>(i*12.0f), (*fieldVec).at(1));
   }
   EXPECT_EQ(true, std::isinf(float16View(200)));
   EXPECT_EQ(true, std::isinf(vec16View(200).at(1)));
   EXPECT_EQ(true, std::isinf(float16View(201)));
   EXPECT_EQ(true, std::isinf(vec16View(201).at(1)));
   EXPECT_EQ(true, std::isnan(float16View(202)));
   EXPECT_EQ(true, std::isnan(vec16View(202).at(1)));
   std::ostringstream os;
   // Only check if calling PrintInfo does not lead to error.
   ntuple->PrintInfo(ENTupleInfo::kSummary, os);
   ntuple->PrintInfo(ENTupleInfo::kStorageDetails, os);
}

TEST(RNTupleFloat, float8_t)
{
   FileRaii fileGuard("test_ntuple_float_float8_t.root");
   const std::string_view ntupleName{"float8_t NTuple"};
   {
      auto model = RNTupleModel::Create();
      auto fieldPt = model->MakeField<float8_t>("ft8");
      auto fieldVec = model->MakeField<std::vector<float8_t>>("ft8vec");
      auto ntuple = RNTupleWriter::Recreate(std::move(model), ntupleName, fileGuard.GetPath());
      for (int i = 0; i < 32; ++i) {
         *fieldPt = static_cast<float>(0.5f * i);
         *fieldVec = { *fieldPt, *fieldPt, *fieldPt };
         ntuple->Fill();
      }
      *fieldPt = std::numeric_limits<float>::infinity() * (-1.0f); // negative infinity
      *fieldVec = { *fieldPt };
      ntuple->Fill();
      *fieldPt = std::numeric_limits<float>::infinity(); // positive infinity
      *fieldVec = { *fieldPt, *fieldPt, *fieldPt, *fieldPt };
      ntuple->Fill();
      *fieldPt = static_cast<float>(std::nanf("1")); // NaN (Not A Number)
      *fieldVec = { *fieldPt };
      ntuple->Fill();
   }
   auto model = RNTupleModel::Create();
   auto fieldPt = model->MakeField<float8_t>("ft8");
   auto fieldVec = model->MakeField<std::vector<float8_t>>("ft8vec");
   auto ntuple = RNTupleReader::Open(std::move(model), ntupleName, fileGuard.GetPath());
   auto float8View = ntuple->GetView<float8_t>("ft8");
   auto vec8View = ntuple->GetView<std::vector<float8_t>>("ft8vec");
   for (int i = 0; i < 32; ++i) {
      EXPECT_EQ(static_cast<float>(0.5f*i), float8View(i)); // EXPECT_FLOAT_EQ will fail here, so use EXPECT_NEAR instead.
      EXPECT_EQ(static_cast<float>(0.5f*i), vec8View(i).at(2));
      ntuple->LoadEntry(i);
      EXPECT_EQ(static_cast<float>(0.5f*i), *fieldPt);
      EXPECT_EQ(static_cast<float>(0.5f*i), (*fieldVec).at(2));
   }
   EXPECT_EQ(true, std::isinf(float8View(32)));
   EXPECT_EQ(true, std::isinf(vec8View(32).at(0)));
   EXPECT_EQ(true, std::isinf(float8View(33)));
   EXPECT_EQ(true, std::isinf(vec8View(33).at(3)));
   EXPECT_EQ(true, std::isnan(float8View(34)));
   EXPECT_EQ(true, std::isnan(vec8View(34).at(0)));
   std::ostringstream os;
   // Only check if calling PrintInfo does not lead to error.
   ntuple->PrintInfo(ENTupleInfo::kSummary, os);
   ntuple->PrintInfo(ENTupleInfo::kStorageDetails, os);
}

// Test the case where number of bits is divisible by 8.
TEST(RNTupleMinMaxDefinedFloat, 24bitDouble)
{
   FileRaii fileGuard("test_ntuple_float_MinMaxDefinedFloat_24bit.root");
   const std::string_view ntupleName{"CustomFloatNTuple24bit"};
   {
      auto model = RNTupleModel::Create();
      auto fieldPt = model->MakeField<double, 24, 0, 1>("ftcustom");
      auto fieldVec = model->MakeField<std::vector<double>, 24, 0, 1>("ftcustom2");
      auto ntuple = RNTupleWriter::Recreate(std::move(model), ntupleName, fileGuard.GetPath());
      for (int i = 0; i < 100; ++i) {
         *fieldPt = static_cast<double>(i*0.01);
         *fieldVec = { static_cast<double>(i*0.005), static_cast<double>(i*0.01) };
         ntuple->Fill();
      }
      *fieldPt = std::numeric_limits<double>::infinity() * (-1.0); // negative infinity
      *fieldVec = { *fieldPt };
      ntuple->Fill();
      *fieldPt = std::numeric_limits<double>::infinity(); // positive infinity
      *fieldVec = { *fieldPt, *fieldPt, *fieldPt, *fieldPt };
      ntuple->Fill();
      *fieldPt = static_cast<double>(std::nan("1")); // NaN (Not A Number)
      *fieldVec = { *fieldPt };
      ntuple->Fill();
   } // flush contents to .root file.

   auto model = RNTupleModel::Create();
   auto fieldPt = model->MakeField<double, 24, 0, 1>("ftcustom");
   auto fieldVec = model->MakeField<std::vector<double>, 24, 0, 1>("ftcustom2");
   auto ntuple = RNTupleReader::Open(std::move(model), ntupleName, fileGuard.GetPath());
   auto view = ntuple->GetView<double, 24, 0, 1>("ftcustom");
   for (int i = 0; i < 100; ++i) {
      EXPECT_NEAR(static_cast<double>(i*0.01), view(i), 0.000001);
      ntuple->LoadEntry(i);
      EXPECT_NEAR(static_cast<double>(i*0.01), *fieldPt, 0.000001);
      EXPECT_NEAR(static_cast<double>(i*0.01), (*fieldVec).at(1), 0.000001);
   }
   EXPECT_EQ(true, std::isinf(view(100)));
   ntuple->LoadEntry(100);
   EXPECT_EQ(true, std::isinf(*fieldPt));
   EXPECT_EQ(true, std::isinf((*fieldVec).at(0)));
   EXPECT_EQ(true, std::isinf(view(101)));
   ntuple->LoadEntry(101);
   EXPECT_EQ(true, std::isinf((*fieldVec).at(3)));
   EXPECT_EQ(true, std::isnan(view(102)));
   ntuple->LoadEntry(102);
   EXPECT_EQ(true, std::isnan((*fieldVec).at(0)));
}

TEST(RNTupleMinMaxDefinedFloat, 47bitDouble)
{
   FileRaii fileGuard("test_ntuple_float_MinMaxDefinedFloat_47bit.root");
   const std::string_view ntupleName{"CustomFloatNTuple47bit"};
   {
      auto model = RNTupleModel::Create();
      auto fieldPt = model->MakeField<double, 47, -50, 100>("ftcustom");
      auto fieldVec = model->MakeField<std::vector<double>, 47, -50, 100>("ftcustom2");
      auto ntuple = RNTupleWriter::Recreate(std::move(model), ntupleName, fileGuard.GetPath());
      for (int i = -500; i <= 1000; ++i) {
         *fieldPt = static_cast<double>(i*0.1);
         *fieldVec = { static_cast<double>(4.0), static_cast<double>(i*0.1) };
         ntuple->Fill();
      }
      *fieldPt = std::numeric_limits<double>::infinity() * (-1.0); // negative infinity
      *fieldVec = { *fieldPt };
      ntuple->Fill();
      *fieldPt = std::numeric_limits<double>::infinity(); // positive infinity
      *fieldVec = { *fieldPt, *fieldPt, *fieldPt, *fieldPt };
      ntuple->Fill();
      *fieldPt = static_cast<double>(std::nan("1")); // NaN (Not A Number)
      *fieldVec = { *fieldPt };
      ntuple->Fill();
   } // flush contents to .root file.

   auto model = RNTupleModel::Create();
   auto fieldPt = model->MakeField<double, 47, -50, 100>("ftcustom");
   auto fieldVec = model->MakeField<std::vector<double>, 47, -50, 100>("ftcustom2");
   auto ntuple = RNTupleReader::Open(std::move(model), ntupleName, fileGuard.GetPath());
   auto view = ntuple->GetView<double, 47, -50, 100>("ftcustom");
   auto viewVec = ntuple->GetView<std::vector<double>, 47, -50, 100>("ftcustom2");
   for (int i = 0; i < 1501; ++i) {
      EXPECT_NEAR(static_cast<double>(((int)i-500)*0.1), view(i), 0.000001);
      EXPECT_NEAR(static_cast<double>(((int)i-500)*0.1), viewVec(i).at(1), 0.000001);
      ntuple->LoadEntry(i);
      EXPECT_NEAR(static_cast<double>(((int)i-500)*0.1), *fieldPt, 0.000001);
      EXPECT_NEAR(static_cast<double>(((int)i-500)*0.1), (*fieldVec).at(1), 0.000001);
   }
   EXPECT_EQ(true, std::isinf(view(1501)));
   EXPECT_EQ(true, std::isinf(viewVec(1501).at(0)));
   ntuple->LoadEntry(1501);
   EXPECT_EQ(true, std::isinf(*fieldPt));
   EXPECT_EQ(true, std::isinf((*fieldVec).at(0)));
   EXPECT_EQ(true, std::isinf(view(1502)));
   EXPECT_EQ(true, std::isinf(viewVec(1502).at(3)));
   ntuple->LoadEntry(1502);
   EXPECT_EQ(true, std::isinf(*fieldPt));
   EXPECT_EQ(true, std::isinf((*fieldVec).at(2)));
   EXPECT_EQ(true, std::isnan(view(1503)));
   EXPECT_EQ(true, std::isnan(viewVec(1503).at(0)));
   ntuple->LoadEntry(1503);
   EXPECT_EQ(true, std::isnan(*fieldPt));
   EXPECT_EQ(true, std::isnan((*fieldVec).at(0)));
   
   // Only checks if calling these functions leads to error, output is not checked.
   auto ntuple2 = RNTupleReader::Open(ntupleName, fileGuard.GetPath());
   std::ostringstream os;
   ntuple->PrintInfo(ENTupleInfo::kSummary, os);
   ntuple->PrintInfo(ENTupleInfo::kStorageDetails, os);
}

TEST(RNTupleMinMaxDefinedFloat, 7bitDouble)
{
   FileRaii fileGuard("test_ntuple_float_MinMaxDefinedFloat_7bit.root");
   const std::string_view ntupleName{"CustomFloatNTuple7bit"};
   {
      auto model = RNTupleModel::Create();
      auto fieldPt = model->MakeField<double, 7, -101, -80>("ftcustom");
      auto fieldArray = model->MakeField<std::array<double, 2>, 7, -101, -80>("ftcustom2");
      auto ntuple = RNTupleWriter::Recreate(std::move(model), ntupleName, fileGuard.GetPath());
      for (int i = -100; i <= -80; ++i) {
         *fieldPt = i;
         *fieldArray = {static_cast<double>(i), static_cast<double>(i)};
         ntuple->Fill();
      }
      *fieldPt = std::numeric_limits<double>::infinity() * (-1.0); // negative infinity
      *fieldArray = { *fieldPt , 0};
      ntuple->Fill();
      *fieldPt = std::numeric_limits<double>::infinity(); // positive infinity
      *fieldArray = { *fieldPt, *fieldPt };
      ntuple->Fill();
      *fieldPt = static_cast<double>(std::nan("1")); // NaN (Not A Number)
      *fieldArray = { *fieldPt, 0 };
      ntuple->Fill();
   } // flush contents to .root file.

   // generate model from descriptor
   auto ntuple = RNTupleReader::Open( ntupleName, fileGuard.GetPath());
   auto view = ntuple->GetView<double, 7, -101, -80>("ftcustom");
   auto arrayView = ntuple->GetView<std::array<double, 2>, 7, -101, -80>("ftcustom2");
   for (int i = -100; i <= -80; ++i) {
      EXPECT_NEAR(static_cast<double>(i), view(i+100), 0.1);
      EXPECT_NEAR(static_cast<double>(i), arrayView(i+100).at(1), 0.1);
   }
   EXPECT_EQ(true, std::isinf(view(21)));
   EXPECT_EQ(true, std::isinf(arrayView(21).at(0)));
   EXPECT_EQ(true, std::isinf(view(22)));
   EXPECT_EQ(true, std::isinf(arrayView(22).at(1)));
   EXPECT_EQ(true, std::isnan(view(23)));
   EXPECT_EQ(true, std::isnan(arrayView(23).at(0)));
}

TEST(RNTupleMinMaxDefinedFloat, 40bitFloat)
{
   FileRaii fileGuard("test_ntuple_float_MinMaxDefinedFloat_40bitFloat.root");
   const std::string_view ntupleName{"CustomFloatNTuple40bitFloat"};
   {
      auto model = RNTupleModel::Create();
      auto fieldPt = model->MakeField<float, 40, 0, 1>("ftcustom");
      auto arrayField = model->MakeField<std::array<float, 3>, 40, 0, 1>("ftcustom2");
      auto ntuple = RNTupleWriter::Recreate(std::move(model), ntupleName, fileGuard.GetPath());
      for (int i = 0; i < 100; ++i) {
         *fieldPt = static_cast<float>(i*0.01f);
         *arrayField = {static_cast<float>(i*0.01f), static_cast<float>(i*0.005f), static_cast<float>(0)};
         ntuple->Fill();
      }
      *fieldPt = std::numeric_limits<float>::infinity() * (-1.0f); // negative infinity
      *arrayField = { *fieldPt, 0, 0 };
      ntuple->Fill();
      *fieldPt = std::numeric_limits<float>::infinity(); // positive infinity
      *arrayField = { *fieldPt, *fieldPt, *fieldPt };
      ntuple->Fill();
      *fieldPt = static_cast<float>(std::nanf("1")); // NaN (Not A Number)
      *arrayField = { *fieldPt, 0, *fieldPt };
      ntuple->Fill();
   } // flush contents to .root file.

   auto model = RNTupleModel::Create();
   auto fieldPt = model->MakeField<float, 40, 0, 1>("ftcustom");
   auto arrayField = model->MakeField<std::array<float, 3>, 40, 0, 1>("ftcustom2");
   auto ntuple = RNTupleReader::Open(std::move(model), ntupleName, fileGuard.GetPath());
   auto view = ntuple->GetView<float, 40, 0, 1>("ftcustom");
   auto arrayView = ntuple->GetView<std::array<float, 3>, 40, 0, 1>("ftcustom2");
   for (int i = 0; i < 100; ++i) {
      EXPECT_NEAR(static_cast<float>(i*0.01), view(i), 0.000001);
      EXPECT_NEAR(static_cast<float>(i*0.01), arrayView(i).at(0), 0.000001);
      ntuple->LoadEntry(i);
      EXPECT_NEAR(static_cast<float>(i*0.01), *fieldPt, 0.000001);
      EXPECT_NEAR(static_cast<float>(i*0.01), (*arrayField).at(0), 0.000001);
   }
   EXPECT_EQ(true, std::isinf(view(100)));
   EXPECT_EQ(true, std::isinf(arrayView(100).at(0)));
   ntuple->LoadEntry(100);
   EXPECT_EQ(true, std::isinf(*fieldPt));
   EXPECT_EQ(true, std::isinf((*arrayField).at(0)));
   EXPECT_EQ(true, std::isinf(view(101)));
   EXPECT_EQ(true, std::isinf(arrayView(101).at(2)));
   ntuple->LoadEntry(101);
   EXPECT_EQ(true, std::isinf(*fieldPt));
   EXPECT_EQ(true, std::isinf((*arrayField).at(2)));
   EXPECT_EQ(true, std::isnan(view(102)));
   EXPECT_EQ(true, std::isnan(arrayView(102).at(2)));
   ntuple->LoadEntry(102);
   EXPECT_EQ(true, std::isnan(*fieldPt));
   EXPECT_EQ(true, std::isnan((*arrayField).at(2)));
}

TEST(RNTupleMinMaxDefinedFloat, 23bitFloat)
{
   FileRaii fileGuard("test_ntuple_float_MinMaxDefinedFloat_23bitFloat.root");
   const std::string_view ntupleName{"CustomFloatNTuple23bitFloat"};
   {
      auto model = RNTupleModel::Create();
      auto fieldPt = model->MakeField<float, 23, -50, 100>("ftcustom");
      auto ntuple = RNTupleWriter::Recreate(std::move(model), ntupleName, fileGuard.GetPath());
      for (int i = -500; i <= 1000; ++i) {
         *fieldPt = static_cast<float>(i*0.1);
         ntuple->Fill();
      }
      *fieldPt = std::numeric_limits<float>::infinity() * (-1.0f); // negative infinity
      ntuple->Fill();
      *fieldPt = std::numeric_limits<float>::infinity(); // positive infinity
      ntuple->Fill();
      *fieldPt = static_cast<float>(std::nanf("1")); // NaN (Not A Number)
      ntuple->Fill();
   } // flush contents to .root file.

   auto model = RNTupleModel::Create();
   auto fieldPt = model->MakeField<float, 23, -50, 100>("ftcustom");
   auto ntuple = RNTupleReader::Open(std::move(model), ntupleName, fileGuard.GetPath());
   auto view = ntuple->GetView<float, 23, -50, 100>("ftcustom");
   for (int i = 0; i < 1501; ++i) {
      EXPECT_NEAR(static_cast<float>(((int)i-500)*0.1), view(i), 0.0001);
      ntuple->LoadEntry(i);
      EXPECT_NEAR(static_cast<float>(((int)i-500)*0.1), *fieldPt, 0.0001);
   }
   EXPECT_EQ(true, std::isinf(view(1501)));
   ntuple->LoadEntry(1501);
   EXPECT_EQ(true, std::isinf(*fieldPt));
   EXPECT_EQ(true, std::isinf(view(1502)));
   ntuple->LoadEntry(1502);
   EXPECT_EQ(true, std::isinf(*fieldPt));
   EXPECT_EQ(true, std::isnan(view(1503)));
   ntuple->LoadEntry(1503);
   EXPECT_EQ(true, std::isnan(*fieldPt));
}

TEST(RNTupleMinMaxDefinedFloat, 5bitFloat)
{
   FileRaii fileGuard("test_ntuple_float_MinMaxDefinedFloat_5bitFloat.root");
   const std::string_view ntupleName{"CustomFloatNTuple5bitFloat"};
   {
      auto model = RNTupleModel::Create();
      auto fieldPt = model->MakeField<float, 5, -101, -80>("ftcustom");
      auto vecField = model->MakeField<std::vector<float>, 5, -101, -80>("ftcustom2");
      auto ntuple = RNTupleWriter::Recreate(std::move(model), ntupleName, fileGuard.GetPath());
      for (int i = -100; i <= -80; ++i) { // eig. bis 80
         *fieldPt = static_cast<float>(i);
         *vecField = { static_cast<float>(i), -100.0f };
         ntuple->Fill();
      }
      *fieldPt = std::numeric_limits<float>::infinity() * (-1.0f); // negative infinity
      *vecField = { *fieldPt };
      ntuple->Fill();
      *fieldPt = std::numeric_limits<float>::infinity(); // positive infinity
      *vecField = { *fieldPt, *fieldPt, *fieldPt, *fieldPt };
      ntuple->Fill();
      *fieldPt = static_cast<float>(std::nanf("1")); // NaN (Not A Number)
      *vecField = { *fieldPt };
      ntuple->Fill();
   } // flush contents to .root file.
   auto model = RNTupleModel::Create();
   auto fieldPt = model->MakeField<float, 5, -101, -80>("ftcustom");
   auto vecField = model->MakeField<std::vector<float>, 5, -101, -80>("ftcustom2");
   auto ntuple = RNTupleReader::Open(std::move(model), ntupleName, fileGuard.GetPath());
   auto view = ntuple->GetView<float, 5, -101, -80>("ftcustom");
   auto vecView = ntuple->GetView<std::vector<float>, 5, -101, -80>("ftcustom2");
   for (int i = -100; i <= -80; ++i) {
      EXPECT_NEAR(static_cast<float>(i), view(i+100), 0.25);
      EXPECT_NEAR(static_cast<float>(i), vecView(i+100).at(0), 0.25);
      ntuple->LoadEntry(i+100);
      EXPECT_NEAR(static_cast<float>(i), *fieldPt, 0.25);
      EXPECT_NEAR(static_cast<float>(i), (*vecField).at(0), 0.25);
   }
   EXPECT_EQ(true, std::isinf(view(21)));
   EXPECT_EQ(true, std::isinf(vecView(21).at(0)));
   ntuple->LoadEntry(21);
   EXPECT_EQ(true, std::isinf(*fieldPt));
   EXPECT_EQ(true, std::isinf((*vecField).at(0)));
   EXPECT_EQ(true, std::isinf(view(22)));
   EXPECT_EQ(true, std::isinf(vecView(22).at(3)));
   ntuple->LoadEntry(22);
   EXPECT_EQ(true, std::isinf(*fieldPt));
   EXPECT_EQ(true, std::isinf((*vecField).at(3)));
   EXPECT_EQ(true, std::isnan(view(23)));
   EXPECT_EQ(true, std::isnan(vecView(23).at(0)));
   ntuple->LoadEntry(23);
   EXPECT_EQ(true, std::isnan(*fieldPt));
   EXPECT_EQ(true, std::isnan((*vecField).at(0)));

   // Only checks if calling these functions leads to error, output is not checked.
   auto ntuple2 = RNTupleReader::Open(ntupleName, fileGuard.GetPath());
   std::ostringstream os;
   ntuple->PrintInfo(ENTupleInfo::kSummary, os);
   ntuple->PrintInfo(ENTupleInfo::kStorageDetails, os);
}
