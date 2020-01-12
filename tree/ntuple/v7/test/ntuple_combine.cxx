// IMPORTANT: Do not shuffle the order of tests. Later tests assume that files created in previous tests already exist!
// TODO (lesimon): Check if RNTupleReader::Show() leads to desired output.
#include <ROOT/RField.hxx>
#include <ROOT/RNTuple.hxx>
#include <ROOT/RNTupleDescriptor.hxx>
#include <ROOT/RNTupleModel.hxx>

#include <TFile.h>

#include "gtest/gtest.h"

#include <iostream>
#include <memory>
#include <string>
#include <vector>

// Uncomment line below to test raw-files instead of .root-files
// #define ROOT7_RAWFILE

// Uncomment line below if the function RNTupleReader::GetDescriptor() exists. This allows a more extensive test.
// #define ROOT7_GETDESCRIPTOR_IN_RNTUPLEREADER

using RFieldBase = ROOT::Experimental::Detail::RFieldBase;
using RNTupleReader = ROOT::Experimental::RNTupleReader;
using RNTupleWriter = ROOT::Experimental::RNTupleWriter;
using RNTupleModel = ROOT::Experimental::RNTupleModel;
using RNTupleDescriptor = ROOT::Experimental::RNTupleDescriptor;
using EFileOpeningOptions = ROOT::Experimental::EFileOpeningOptions;


namespace {
   /**
    * An RAII wrapper around an open temporary file on disk. It cleans up the guarded file when the wrapper object
    * goes out of scope.
    */
   class FileRaii {
   private:
      std::string fPath;
      TFile *file;
   public:
      FileRaii(const std::string &path) :
         fPath(path),
         file(TFile::Open(fPath.c_str(), "RECREATE"))
         { }
      FileRaii(const FileRaii&) = delete;
      FileRaii& operator=(const FileRaii&) = delete;
      TFile* GetFile() const { return file; }
      std::string GetPath() const { return fPath; }
      ~FileRaii() {
         file->Close();
         std::remove(fPath.c_str());
      }
   };
} // anonymous namespace
#ifndef ROOT7_RAWFILE
static const std::string FileEnding(".root");
#endif
#ifdef ROOT7_RAWFILE
static const std::string FileEnding(".ntuple");
#endif
// Define static const fileGuarads here, a single file will be used for multiple tests.
static const FileRaii fileGuardChain1("test_ntuple_chain_File1"+FileEnding); // 10 entries
static const FileRaii fileGuardChain2("test_ntuple_chain_File2"+FileEnding); // 75000 entries
static const FileRaii fileGuardChain3("test_ntuple_chain_File3"+FileEnding); // 30000 entries
static const FileRaii fileGuardChain4("test_ntuple_chain_File4"+FileEnding); // 30000 entries
static const FileRaii fileGuardChain1twiceMerged("test_ntuple_chain_File1twice"+FileEnding);
static const FileRaii fileGuardChain1234Merged("test_ntuple_chain_File1234"+FileEnding);
static const FileRaii fileGuardChain134Merged("test_ntuple_chain_File134"+FileEnding);
static const FileRaii fileGuardChainWrongNumFields("test_ntuple_chain_WrongNumFields"+FileEnding);
static const FileRaii fileGuardChainWrongFieldName("test_ntuple_chain_WrongFieldName"+FileEnding);

static const FileRaii fileGuardFriend1("test_ntuple_friend_File1"+FileEnding); // 11000 entries
static const FileRaii fileGuardFriend2("test_ntuple_friend_File2"+FileEnding); // 11000 entries
static const FileRaii fileGuardFriend3("test_ntuple_friend_File3"+FileEnding); // 11000 entries
static const FileRaii fileGuardFriend4("test_ntuple_friend_File4"+FileEnding); // 11000 entries
static const FileRaii fileGuardFriend1234Merged("test_ntuple_friend_File1234"+FileEnding);
static const FileRaii fileGuardFriendWrongClusterData("test_ntuple_friend_WrongClusterData"+FileEnding); // 11000 entries

static const FileRaii &fileGuardMixed1(fileGuardChain3);
static const FileRaii &fileGuardMixed2(fileGuardChain4);
static const FileRaii fileGuardMixed3("test_ntuple_mixed_File3"+FileEnding);
static const FileRaii fileGuardMixed4("test_ntuple_mixed_File4"+FileEnding);
static const FileRaii fileGuardMixed1234("test_ntuple_mixed_File1234"+FileEnding);

// Note(lesimon): No separate names for files which will get merged together. Files which get merged together require the same ntupleName.
static const std::string_view ntupleNameChain("chainNTuple");
static const std::string_view ntupleNameChain1twiceMerged("chainNTuple1twice");
static const std::string_view ntupleNameChain1234Merged("chainNTuple1234");
static const std::string_view ntupleNameChain134Merged("chainNTuple134");

static const std::string_view ntupleNameFriend("friendNTuple");
static const std::string_view ntupleNameFriend1234Merged("friendNTuple1234Merged");

static const std::string_view ntupleNameMixed(ntupleNameChain);
static const std::string_view ntupleNameMixed1234Merged("mixedNTuple1234");

// ------------------------------------------- Tests ------------------------------------------

TEST(RNTupleChain, noFiles)
{
   const std::string_view ntupleName = "noFilesNTuple";
   auto ntuple = RNTupleReader::Open(ntupleName,  std::vector<std::string>{ });
   EXPECT_EQ(ntuple, nullptr);
}


TEST(RNTupleChain, oneFile)
{
   // Create: test_ntuple_chain_File1.root
   {
      auto model = RNTupleModel::Create();
      auto dbField = model->MakeField<double>("db");
      auto stField = model->MakeField<std::string>("st");
      auto itField = model->MakeField<std::vector<std::int32_t>>("it");
      auto ntuple = RNTupleWriter::Recreate(std::move(model), ntupleNameChain, fileGuardChain1.GetPath());

      for (int i = 1; i < 11; ++i) {
         *dbField = 5.0+i;
         *stField = "foo";
         std::vector<std::int32_t> intVec(i);
         for (int j = 0; j < i; ++j) {
            intVec.at(j) = j;
         }
         *itField = intVec;
         ntuple->Fill();
      }
   } // flush contents to test_ntuple_chain_File1.root

   auto ntupleNotChain = RNTupleReader::Open(ntupleNameChain, fileGuardChain1.GetPath());
   auto ntupleChain = RNTupleReader::Open(ntupleNameChain, std::vector<std::string>{ fileGuardChain1.GetPath() });
   auto dbView = ntupleChain->GetView<double>("db");
   auto stView = ntupleChain->GetView<std::string>("st");
   auto itView = ntupleChain->GetView<std::vector<std::int32_t>>("it");
   auto dbView2 = ntupleNotChain->GetView<double>("db");
   auto stView2 = ntupleNotChain->GetView<std::string>("st");
   auto itView2 = ntupleNotChain->GetView<std::vector<std::int32_t>>("it");
   
   EXPECT_NE(nullptr, ntupleChain);
   for (auto i : ntupleChain->GetViewRange()) {
      EXPECT_DOUBLE_EQ(dbView2(i), dbView(i));
      EXPECT_EQ(stView2(i), stView(i));
      EXPECT_EQ(itView2(i), itView(i));
   }
}


TEST(RNTupleChain, twiceSameFile)
{
   // Create: test_ntuple_chain_File1twice.root
   {
      auto model = RNTupleModel::Create();
      auto dbField = model->MakeField<double>("db");
      auto stField = model->MakeField<std::string>("st");
      auto itField = model->MakeField<std::vector<std::int32_t>>("it");
      auto ntuple = RNTupleWriter::Recreate(std::move(model), ntupleNameChain1twiceMerged, fileGuardChain1twiceMerged.GetPath());

      for (int j = 0; j < 2; ++j) {
         for (int i = 1; i < 11; ++i) {
            *dbField = 5.0+i;
            *stField = "foo";
            std::vector<std::int32_t> intVec(i);
            for (int k = 0; k < i; ++k) {
               intVec.at(k) = k;
            }
            *itField = intVec;
            ntuple->Fill();
         }
         if (j == 0)
            ntuple->CommitCluster(); // required for descriptorTest
      }
   } // flush contents to test_ntuple_chain_File1twice.root

   auto ntupleChain = RNTupleReader::Open(ntupleNameChain, { fileGuardChain1.GetPath(), fileGuardChain1.GetPath() } );
   auto ntupleNotChain = RNTupleReader::Open(ntupleNameChain1twiceMerged, fileGuardChain1twiceMerged.GetPath() );

   auto dbView = ntupleChain->GetView<double>("db");
   auto stView = ntupleChain->GetView<std::string>("st");
   auto itView = ntupleChain->GetView<std::vector<std::int32_t>>("it");
   EXPECT_DOUBLE_EQ(15.0, dbView(19));
   EXPECT_EQ("foo", stView(14));
   EXPECT_EQ(9, itView(19).at(9));
}


TEST(RNTupleChain, fourFiles)
{
   // Create: test_ntuple_chain_File2.root
   {
      auto model = RNTupleModel::Create();
      auto dbField = model->MakeField<double>("db");
      auto stField = model->MakeField<std::string>("st");
      auto itField = model->MakeField<std::vector<std::int32_t>>("it");
      auto ntuple = RNTupleWriter::Recreate(std::move(model), ntupleNameChain, fileGuardChain2.GetPath());
      for (int j = 0; j < 10; ++j) {
         for (int i = 0; i < 7500; ++i) {
            *dbField = 5.0+i;
            *stField = "foo" + std::to_string(i);
            std::vector<std::int32_t> intVec(j);
            for (int k = 0; k < j; ++k) {
               intVec.at(k) = k;
            }
            *itField = intVec;
            ntuple->Fill();
         }
      }
   } // flush contents to test_ntuple_chain_File2.root

   // Create: test_ntuple_chain_File3.root
   {
      auto model = RNTupleModel::Create();
      auto dbField = model->MakeField<double>("db");
      auto stField = model->MakeField<std::string>("st");
      auto itField = model->MakeField<std::vector<std::int32_t>>("it");
      auto ntuple = RNTupleWriter::Recreate(std::move(model), ntupleNameChain, fileGuardChain3.GetPath());
      for (int j = 0; j < 30000; ++j) {
         *dbField = 5.0+j;
         *stField = "foo" + std::to_string(j);
         std::vector<std::int32_t> intVec(3);
         for (int k = 0; k < 3; ++k) {
            intVec.at(k) = k;
         }
         *itField = intVec;
         ntuple->Fill();
      }
   } // flush contents to test_ntuple_chain_File3.root
   
   // Create: test_ntuple_chain_File4.root
   {
      auto model = RNTupleModel::Create();
      auto dbField = model->MakeField<double>("db");
      auto stField = model->MakeField<std::string>("st");
      auto itField = model->MakeField<std::vector<std::int32_t>>("it");
      auto ntuple = RNTupleWriter::Recreate(std::move(model), ntupleNameChain, fileGuardChain4.GetPath());
      for (int j = 0; j < 30000; ++j) {
         *dbField = 10.11+j;
         *stField = "goo" + std::to_string(j);
         std::vector<std::int32_t> intVec(2);
         for (int k = 0; k < 2; ++k) {
            intVec.at(k) = k;
         }
         *itField = intVec;
         ntuple->Fill();
      }
   } // flush contents to test_ntuple_chain_File4.root



   // Create: test_ntuple_chain_File1234.root
   {
      auto model = RNTupleModel::Create();
      auto dbField = model->MakeField<double>("db");
      auto stField = model->MakeField<std::string>("st");
      auto itField = model->MakeField<std::vector<std::int32_t>>("it");
      auto ntuple = RNTupleWriter::Recreate(std::move(model), ntupleNameChain1234Merged, fileGuardChain1234Merged.GetPath());
      for (int i = 1; i < 11; ++i) {
         *dbField = 5.0+i;
         *stField = "foo";
         std::vector<std::int32_t> intVec(i);
         for (int j = 0; j < i; ++j) {
            intVec.at(j) = j;
         }
         *itField = intVec;
         ntuple->Fill();
      }
      for (int j = 0; j < 10; ++j) {
         for (int i = 0; i < 7500; ++i) {
            *dbField = 5.0+i;
            *stField = "foo" + std::to_string(i);
            std::vector<std::int32_t> intVec(j);
            for (int k = 0; k < j; ++k) {
               intVec.at(k) = k;
            }
            *itField = intVec;
            ntuple->Fill();
         }
      }
      for (int j = 0; j < 30000; ++j) {
         *dbField = 5.0+j;
         *stField = "foo" + std::to_string(j);
         std::vector<std::int32_t> intVec(3);
         for (int k = 0; k < 3; ++k) {
            intVec.at(k) = k;
         }
         *itField = intVec;
         ntuple->Fill();
      }
      for (int j = 0; j < 30000; ++j) {
         *dbField = 10.11+j;
         *stField = "goo" + std::to_string(j);
         std::vector<std::int32_t> intVec(2);
         for (int k = 0; k < 2; ++k) {
            intVec.at(k) = k;
         }
         *itField = intVec;
         ntuple->Fill();
      }
   } // flush contents to test_ntuple_chain_File1234.root

   auto ntupleChain = RNTupleReader::Open(ntupleNameChain, { fileGuardChain1.GetPath(), fileGuardChain2.GetPath(), fileGuardChain3.GetPath(), fileGuardChain4.GetPath() });
   auto ntupleNotChain = RNTupleReader::Open(ntupleNameChain1234Merged, fileGuardChain1234Merged.GetPath());
   EXPECT_EQ(ntupleNotChain->GetNEntries(), ntupleChain->GetNEntries());
   
   auto dbView = ntupleChain->GetView<double>("db");
   auto stView = ntupleChain->GetView<std::string>("st");
   auto itView = ntupleChain->GetView<std::vector<std::int32_t>>("it");
   auto dbViewNotChain = ntupleNotChain->GetView<double>("db");
   auto stViewNotChain = ntupleNotChain->GetView<std::string>("st");
   auto itViewNotChain = ntupleNotChain->GetView<std::vector<std::int32_t>>("it");
   
   for (auto i : ntupleChain->GetViewRange()) {
      EXPECT_DOUBLE_EQ(dbViewNotChain(i), dbView(i));
      EXPECT_EQ(stViewNotChain(i), stView(i));
      EXPECT_EQ(itViewNotChain(i), itView(i));
   }
}


#ifdef ROOT7_GETDESCRIPTOR_IN_RNTUPLEREADER
// -------------------------- helper functions ---------------------------------

bool CompareFields (const RNTupleDescriptor &desc1, const RNTupleDescriptor &desc2) {
   if (desc1.GetNFields() != desc2.GetNFields())
      return false;
   for (std::size_t i = 0; i < desc1.GetNFields(); ++i) {
      if (!(desc1.GetFieldDescriptor(i) == desc2.GetFieldDescriptor(i)))
         return false;
   }
   return true;
}

bool CompareColumns (const RNTupleDescriptor &desc1, const RNTupleDescriptor &desc2) {
   if(desc1.GetNColumns() != desc2.GetNColumns())
      return false;
   for (std::size_t i = 0; i < desc1.GetNColumns(); ++i) {
      if (!(desc1.GetColumnDescriptor(i) == desc2.GetColumnDescriptor(i)))
         return false;
   }
   return true;
}

bool CompareClusters (const RNTupleDescriptor &desc1, const RNTupleDescriptor &desc2) {
   if(desc1.GetNClusters() != desc2.GetNClusters()) {
      return false; }
   for (std::size_t i = 0; i < desc1.GetNClusters(); ++i) {
      if (!(desc1.GetClusterDescriptor(i).BeFriendAble(desc2.GetClusterDescriptor(i)))) {
         return false;
      }
      for (int j = desc1.GetNColumns()-1; j >= 0; --j) {
         if (!(desc1.GetClusterDescriptor(i).GetColumnRange(j) == desc2.GetClusterDescriptor(i).GetColumnRange(j)))
            return false;
         if (!(desc1.GetClusterDescriptor(i).GetPageRange(j).fColumnId ==
               desc2.GetClusterDescriptor(i).GetPageRange(j).fColumnId))
            return false;
         if (!(desc1.GetClusterDescriptor(i).GetPageRange(j).fPageInfos.size() ==
               desc2.GetClusterDescriptor(i).GetPageRange(j).fPageInfos.size()))
         return false;
         for (int k = desc1.GetClusterDescriptor(i).GetPageRange(j).fPageInfos.size()-1; k >= 0; --k) {
            if (!(desc1.GetClusterDescriptor(i).GetPageRange(j).fPageInfos.at(k).fNElements ==
                  desc2.GetClusterDescriptor(i).GetPageRange(j).fPageInfos.at(k).fNElements))
               return false;
         }
      }
   }
   return true;
}

bool DescriptorsAreSame(const RNTupleDescriptor &desc1, const RNTupleDescriptor &desc2) {
   return ( CompareFields(desc1, desc2) && CompareColumns(desc1, desc2) && CompareClusters(desc1, desc2) );
}

// -------------------------- Test ---------------------------------

TEST(RNTupleChain, DescriptorTest)
{
   auto ntupleNotChain1234 = RNTupleReader::Open(ntupleNameChain1234Merged, fileGuardChain1234Merged.GetPath());
   // 60010 is the total number of entries in test_ntuple_chain_File1.root, ...3.root and ...4.root combined.
   // If the maximal number of elements a cluster can hold is below that number, the clusterDescriptors will not match.
   ASSERT_TRUE((ntupleNotChain1234->GetDescriptor().GetClusterDescriptor(0).GetNEntries() > 60010) && ("For this test to succeed a cluster should be able to hold more than 60010 entries!\n"));
   
   
   auto ntupleChainFile2 = RNTupleReader::Open(ntupleNameChain, { fileGuardChain2.GetPath() });
   auto ntupleNotChainFile2 = RNTupleReader::Open(ntupleNameChain, fileGuardChain2.GetPath());
   ASSERT_NE(nullptr, ntupleChainFile2);
   EXPECT_TRUE(DescriptorsAreSame(ntupleNotChainFile2->GetDescriptor(), ntupleChainFile2->GetDescriptor()));
   
   auto ntupleChainFile1Twice = RNTupleReader::Open(ntupleNameChain, { fileGuardChain1.GetPath(), fileGuardChain1.GetPath() });
   auto ntupleNotChainFile1Twice = RNTupleReader::Open(ntupleNameChain1twiceMerged, fileGuardChain1twiceMerged.GetPath() );
   ASSERT_NE(nullptr, ntupleChainFile1Twice);
   EXPECT_TRUE(DescriptorsAreSame(ntupleNotChainFile1Twice->GetDescriptor(), ntupleChainFile1Twice->GetDescriptor()));
   
   // Create: test_ntuple_chain_File134.root
   {
      auto model = RNTupleModel::Create();
      auto dbField = model->MakeField<double>("db");
      auto stField = model->MakeField<std::string>("st");
      auto itField = model->MakeField<std::vector<std::int32_t>>("it");
      auto ntuple = RNTupleWriter::Recreate(std::move(model), ntupleNameChain134Merged, fileGuardChain134Merged.GetPath());
      for (int i = 1; i < 11; ++i) {
         *dbField = 5.0+i;
         *stField = "foo";
         std::vector<std::int32_t> intVec(i);
         for (int j = 0; j < i; ++j) {
            intVec.at(j) = j;
         }
         *itField = intVec;
         ntuple->Fill();
      }
      ntuple->CommitCluster();
      for (int j = 0; j < 30000; ++j) {
         *dbField = 5.0+j;
         *stField = "foo" + std::to_string(j);
         std::vector<std::int32_t> intVec(3);
         for (int k = 0; k < 3; ++k) {
            intVec.at(k) = k;
         }
         *itField = intVec;
         ntuple->Fill();
      }
      ntuple->CommitCluster();
      for (int j = 0; j < 30000; ++j) {
         *dbField = 10.11+j;
         *stField = "goo" + std::to_string(j);
         std::vector<std::int32_t> intVec(2);
         for (int k = 0; k < 2; ++k) {
            intVec.at(k) = k;
         }
         *itField = intVec;
         ntuple->Fill();
      }
   } // flush contents to test_ntuple_chain_File134.root
   auto ntupleChainFiles134 = RNTupleReader::Open(ntupleNameChain, { fileGuardChain1.GetPath(), fileGuardChain3.GetPath(), fileGuardChain4.GetPath() });
   auto ntupleNotChainFiles134 = RNTupleReader::Open(ntupleNameChain134Merged, fileGuardChain134Merged.GetPath());
   ASSERT_NE(nullptr, ntupleChainFiles134);
   EXPECT_TRUE(DescriptorsAreSame(ntupleNotChainFiles134->GetDescriptor(), ntupleChainFiles134->GetDescriptor()));
}
#endif


TEST(RNTupleChain, WrongNumberOfFields)
{
   // Create: test_ntuple_chain_WrongNumFields.root
   {
      auto model = RNTupleModel::Create();
      auto dbField = model->MakeField<double>("db");
      auto stField = model->MakeField<std::string>("st");
      auto itField = model->MakeField<std::vector<std::int32_t>>("it");
      auto blField = model->MakeField<bool>("bl");
      auto ntuple = RNTupleWriter::Recreate(std::move(model), ntupleNameChain, fileGuardChainWrongNumFields.GetPath());

      for (int i = 1; i < 11; ++i) {
         *dbField = 5.0+i;
         *stField = "foo";
         std::vector<std::int32_t> intVec(i);
         for (int j = 0; j < i; ++j) {
            intVec.at(j) = j;
         }
         *itField = intVec;
         *blField = true;
         ntuple->Fill();
      }
   } // flush contents to test_ntuple_chain_WrongNumFields.root
   auto ntuple = RNTupleReader::Open(ntupleNameChain, { fileGuardChain1.GetPath(), fileGuardChainWrongNumFields.GetPath() });
   EXPECT_EQ(nullptr, ntuple);
}

TEST(RNTupleChain, WrongFieldName)
{
   // Create: test_ntuple_chain_WrongFieldName.root
   {
      auto model = RNTupleModel::Create();
      auto dbField = model->MakeField<double>("db");
      auto stField = model->MakeField<std::string>("st");
      auto itField = model->MakeField<std::vector<std::int32_t>>("wrongFieldName");
      auto ntuple = RNTupleWriter::Recreate(std::move(model), ntupleNameChain, fileGuardChainWrongFieldName.GetPath());

      for (int i = 1; i < 11; ++i) {
         *dbField = 5.0+i;
         *stField = "foo";
         std::vector<std::int32_t> intVec(i);
         for (int j = 0; j < i; ++j) {
            intVec.at(j) = j;
         }
         *itField = intVec;
         ntuple->Fill();
      }
   } // flush contents to test_ntuple_chain_WrongFieldName.root
   auto ntuple = RNTupleReader::Open(ntupleNameChain, { fileGuardChain1.GetPath(), fileGuardChainWrongFieldName.GetPath() });
   EXPECT_EQ(nullptr, ntuple);
}


// When a cached page is returned in RPageSourceRoot::PopulatePage(), 2 pages with the same buffer will exist.
// This test should check if RNTupleChain can deal with such a situation.
TEST(RNTupleChain, CachedPage)
{
   auto ntuple = RNTupleReader::Open(ntupleNameChain, { fileGuardChain1.GetPath() });
   
   auto dbView = ntuple->GetView<double>("db");
   auto stView = ntuple->GetView<std::string>("st");
   auto itView = ntuple->GetView<std::vector<std::int32_t>>("it");
   auto dbView2 = ntuple->GetView<double>("db");
   auto stView2 = ntuple->GetView<std::string>("st");
   auto itView2 = ntuple->GetView<std::vector<std::int32_t>>("it");
   // 2 viewers populate the same page while none of the pages are being released.
   for (int i = 0; i < 2; ++i) {
      EXPECT_DOUBLE_EQ(dbView(5), dbView2(5));
      EXPECT_EQ(stView(5), stView2(5));
      EXPECT_EQ(itView(5), itView2(5));
   }
}


TEST(RNTupleChain, ChainOfChainOfChain)
{
   auto ntuple1 = RNTupleReader::Open(ntupleNameChain, fileGuardChain1.GetPath());
   auto ntuple2 = RNTupleReader::Open(ntupleNameChain, fileGuardChain2.GetPath());
   auto ntuple3 = RNTupleReader::Open(ntupleNameChain, fileGuardChain3.GetPath());
   auto ntuple4 = RNTupleReader::Open(ntupleNameChain, fileGuardChain4.GetPath());
   auto ntupleChain = RNTupleReader::ChainReader(ntupleNameChain, ntuple1, ntuple2);
   auto ntupleChainOfChain = RNTupleReader::ChainReader(ntupleNameChain, ntupleChain, ntuple3);
   auto ntupleChainOfChainOfChain = RNTupleReader::ChainReader(ntupleNameChain, ntupleChainOfChain, ntuple4);
   auto ntuple1234merged = RNTupleReader::Open(ntupleNameChain1234Merged, fileGuardChain1234Merged.GetPath());
   EXPECT_EQ(ntuple1234merged->GetNEntries(), ntupleChainOfChainOfChain->GetNEntries());
   
   auto dbView = ntupleChainOfChainOfChain->GetView<double>("db");
   auto stView = ntupleChainOfChainOfChain->GetView<std::string>("st");
   auto itView = ntupleChainOfChainOfChain->GetView<std::vector<std::int32_t>>("it");
   auto dbViewNotChain = ntuple1234merged->GetView<double>("db");
   auto stViewNotChain = ntuple1234merged->GetView<std::string>("st");
   auto itViewNotChain = ntuple1234merged->GetView<std::vector<std::int32_t>>("it");
   
   for (auto i : ntupleChainOfChainOfChain->GetViewRange()) {
      EXPECT_DOUBLE_EQ(dbViewNotChain(i), dbView(i));
      EXPECT_EQ(stViewNotChain(i), stView(i));
      EXPECT_EQ(itViewNotChain(i), itView(i));
   }
}


TEST(RNTupleChain, ChainOfChainOfChainWithStdMove)
{
   auto ntuple1 = RNTupleReader::Open(ntupleNameChain, fileGuardChain1.GetPath());
   auto ntuple2 = RNTupleReader::Open(ntupleNameChain, fileGuardChain2.GetPath());
   auto ntuple3 = RNTupleReader::Open(ntupleNameChain, fileGuardChain3.GetPath());
   auto ntuple4 = RNTupleReader::Open(ntupleNameChain, fileGuardChain4.GetPath());
   auto ntupleChain = RNTupleReader::ChainReader(ntupleNameChain, std::move(ntuple1), std::move(ntuple2));
   auto ntupleChainOfChain = RNTupleReader::ChainReader(ntupleNameChain, std::move(ntupleChain), std::move(ntuple3));
   auto ntupleChainOfChainOfChain = RNTupleReader::ChainReader(ntupleNameChain, std::move(ntupleChainOfChain), std::move(ntuple4));
   auto ntuple1234merged = RNTupleReader::Open(ntupleNameChain1234Merged, fileGuardChain1234Merged.GetPath());
   EXPECT_EQ(ntuple1234merged->GetNEntries(), ntupleChainOfChainOfChain->GetNEntries());
   
   auto dbView = ntupleChainOfChainOfChain->GetView<double>("db");
   auto stView = ntupleChainOfChainOfChain->GetView<std::string>("st");
   auto itView = ntupleChainOfChainOfChain->GetView<std::vector<std::int32_t>>("it");
   auto dbViewNotChain = ntuple1234merged->GetView<double>("db");
   auto stViewNotChain = ntuple1234merged->GetView<std::string>("st");
   auto itViewNotChain = ntuple1234merged->GetView<std::vector<std::int32_t>>("it");
   
   for (auto i : ntupleChainOfChainOfChain->GetViewRange()) {
      EXPECT_DOUBLE_EQ(dbViewNotChain(i), dbView(i));
      EXPECT_EQ(stViewNotChain(i), stView(i));
      EXPECT_EQ(itViewNotChain(i), itView(i));
   }
}


TEST(RNTupleFriend, OneFile)
{
   // Create: test_ntuple_friend_File1.root
   {
      auto model = RNTupleModel::Create();
      auto boField = model->MakeField<bool>("bo");
      auto u64Field = model->MakeField<std::uint64_t>("u64");
      auto ntuple = RNTupleWriter::Recreate(std::move(model), ntupleNameFriend, fileGuardFriend1.GetPath());

      for (std::size_t i = 0; i < 11000; ++i) {
         *boField = i%2;
         *u64Field = i;
         ntuple->Fill();
         if (i == 6000)
            ntuple->CommitCluster();
      }
   } // flush contents to test_ntuple_friend_File1.root

   auto ntupleFriend = RNTupleReader::Open(ntupleNameFriend, { fileGuardFriend1.GetPath() }, EFileOpeningOptions::kFriend);
   auto ntupleNotFriend = RNTupleReader::Open(ntupleNameFriend, fileGuardFriend1.GetPath());
   EXPECT_NE(nullptr, ntupleFriend);
   
   auto boView = ntupleFriend->GetView<bool>("bo");
   auto u64View = ntupleFriend->GetView<std::uint64_t>("u64");
   auto boView2 = ntupleNotFriend->GetView<bool>("bo");
   auto u64View2 = ntupleNotFriend->GetView<std::uint64_t>("u64");

   for (auto i : ntupleFriend->GetViewRange()) {
      EXPECT_EQ(boView2(i), boView(i));
      EXPECT_EQ(u64View2(i), u64View(i));
   }
}


TEST(RNTupleFriend, TwiceSameFile)
{
   auto ntuple = RNTupleReader::Open(ntupleNameFriend, { fileGuardFriend1.GetPath(), fileGuardFriend1.GetPath() }, EFileOpeningOptions::kFriend);
   EXPECT_EQ(nullptr, ntuple);
}


TEST(RNTupleFriend, fourFiles)
{
   // Create: test_ntuple_friend_File2.ntuple
   {
      auto model = RNTupleModel::Create();
      auto stField = model->MakeField<std::string>("st");
      auto ntuple = RNTupleWriter::Recreate(std::move(model), ntupleNameFriend, fileGuardFriend2.GetPath());

      for (std::size_t i = 0; i < 11000; ++i) {
         *stField = "foo" + std::to_string(i);
         ntuple->Fill();
         if (i == 6000)
            ntuple->CommitCluster();
      }
   } // flush contents to test_ntuple_friend_File2.root
   
   // Create: test_ntuple_friend_File3.root
   {
      auto model = RNTupleModel::Create();
      auto ftField = model->MakeField<std::vector<float>>("ftVec");
      auto arField = model->MakeField<std::array<std::int32_t, 4>>("ar32");
      auto ntuple = RNTupleWriter::Recreate(std::move(model), ntupleNameFriend, fileGuardFriend3.GetPath());

      for (std::size_t i = 0; i < 11000; ++i) {
         std::vector<float> floatVec(i%3+1);
         for (int j = i%3; j >= 0; --j) {
            floatVec.at(j) = 4.0f;
         }
         *ftField = floatVec;
         *arField = {1, 2, 3, 4};
         ntuple->Fill();
         if (i == 6000)
            ntuple->CommitCluster();
      }
   } // flush contents to test_ntuple_friend_File3.root
   
   // Create: test_ntuple_friend_File4.root
   {
      auto model = RNTupleModel::Create();
      auto boField = model->MakeField<bool>("bo2");
      auto u64Field = model->MakeField<std::uint64_t>("u642");
      auto ntuple = RNTupleWriter::Recreate(std::move(model), ntupleNameFriend, fileGuardFriend4.GetPath());

      for (std::size_t i = 0; i < 11000; ++i) {
         *boField = i%2;
         *u64Field = i;
         ntuple->Fill();
         if (i == 6000)
            ntuple->CommitCluster();
      }
   } // flush contents to test_ntuple_friend_File4.root
   
   // Create: test_ntuple_friend_File1234.ntuple
   {
      auto model = RNTupleModel::Create();
      auto boField = model->MakeField<bool>("bo");
      auto u64Field = model->MakeField<std::uint64_t>("u64");
      auto stField = model->MakeField<std::string>("st");
      auto ftField = model->MakeField<std::vector<float>>("ftVec");
      auto arField = model->MakeField<std::array<std::int32_t, 4>>("ar32");
      auto boField2 = model->MakeField<bool>("bo2");
      auto u64Field2 = model->MakeField<std::uint64_t>("u642");
      auto ntuple = RNTupleWriter::Recreate(std::move(model), ntupleNameFriend1234Merged, fileGuardFriend1234Merged.GetPath());

      for (std::size_t i = 0; i < 11000; ++i) {
         *boField = i%2;
         *u64Field = i;
         *stField = "foo" + std::to_string(i);
         std::vector<float> floatVec(i%3+1);
         for (int j = i%3; j >= 0; --j) {
            floatVec.at(j) = 4.0f;
         }
         *ftField = floatVec;
         *arField = {1, 2, 3, 4};
         *boField2 = i%2;
         *u64Field2 = i;
         ntuple->Fill();
         if (i == 6000)
            ntuple->CommitCluster();
      }
   } // flush contents to test_ntuple_friend_File1234.ntuple

   auto ntupleFriend = RNTupleReader::Open(ntupleNameFriend, { fileGuardFriend1.GetPath(), fileGuardFriend2.GetPath(), fileGuardFriend3.GetPath(), fileGuardFriend4.GetPath() }, EFileOpeningOptions::kFriend);
   auto ntupleNotFriend = RNTupleReader::Open(ntupleNameFriend1234Merged, fileGuardFriend1234Merged.GetPath());
   EXPECT_NE(nullptr, ntupleFriend);

   auto boView = ntupleFriend->GetView<bool>("bo");
   auto u64View = ntupleFriend->GetView<std::uint64_t>("u64");
   auto stView = ntupleFriend->GetView<std::string>("st");
   auto ftView = ntupleFriend->GetView<std::vector<float>>("ftVec");
   auto arView = ntupleFriend->GetView<std::array<std::int32_t, 4>>("ar32");
   auto bo2View = ntupleFriend->GetView<bool>("bo2");
   auto u642View = ntupleFriend->GetView<std::uint64_t>("u642");
   
   auto boView2 = ntupleNotFriend->GetView<bool>("bo");
   auto u64View2 = ntupleNotFriend->GetView<std::uint64_t>("u64");
   auto stView2 = ntupleNotFriend->GetView<std::string>("st");
   auto ftView2 = ntupleNotFriend->GetView<std::vector<float>>("ftVec");
   auto arView2 = ntupleNotFriend->GetView<std::array<std::int32_t, 4>>("ar32");
   auto bo2View2 = ntupleNotFriend->GetView<bool>("bo2");
   auto u642View2 = ntupleNotFriend->GetView<std::uint64_t>("u642");

   for (auto i : ntupleFriend->GetViewRange()) {
      EXPECT_EQ(boView2(i), boView(i));
      EXPECT_EQ(u64View2(i), u64View(i));
      EXPECT_EQ(stView2(i), stView(i));
      EXPECT_EQ(ftView2(i), ftView(i));
      EXPECT_EQ(arView2(i), arView(i));
      EXPECT_EQ(bo2View2(i), bo2View(i));
      EXPECT_EQ(u642View2(i), u642View(i));
   }
}


#ifdef ROOT7_GETDESCRIPTOR_IN_RNTUPLEREADER
TEST(RNTupleFriend, CompareDescriptor)
{
   auto ntupleFriend = RNTupleReader::Open(ntupleNameFriend, { fileGuardFriend1.GetPath(), fileGuardFriend2.GetPath(), fileGuardFriend3.GetPath(), fileGuardFriend4.GetPath() }, EFileOpeningOptions::kFriend);
   auto ntupleNotFriend = RNTupleReader::Open(ntupleNameFriend1234Merged, fileGuardFriend1234Merged.GetPath());
   EXPECT_NE(nullptr, ntupleFriend);
   EXPECT_TRUE(DescriptorsAreSame(ntupleNotFriend->GetDescriptor(), ntupleFriend->GetDescriptor()));
}
#endif


TEST(RNTupleFriend, DifferentNumberOfEntryInCluster)
{
   // Create: test_ntuple_friend_WrongClusterData.root
   {
      auto model = RNTupleModel::Create();
      auto ftField = model->MakeField<float>("ft");
      auto ntuple = RNTupleWriter::Recreate(std::move(model), ntupleNameFriend, fileGuardFriendWrongClusterData.GetPath());

      for (std::size_t i = 0; i < 11000; ++i) {
         *ftField = 5.0f + i;
         ntuple->Fill();
         if (i == 5000)
            ntuple->CommitCluster();
      }
   } // flush contents to test_ntuple_friend_WrongClusterData.root
   
   auto ntuple = RNTupleReader::Open(ntupleNameFriend, { fileGuardFriend1.GetPath(), fileGuardFriendWrongClusterData.GetPath() }, EFileOpeningOptions::kFriend);
   EXPECT_EQ(nullptr, ntuple);
}


TEST(RNTupleFriend, CachedPage)
{
   auto ntuple = RNTupleReader::Open(ntupleNameFriend, { fileGuardFriend1.GetPath(), fileGuardFriend2.GetPath() }, EFileOpeningOptions::kFriend);
   auto boView = ntuple->GetView<bool>("bo");
   auto boView2 = ntuple->GetView<bool>("bo");
   auto boView3 = ntuple->GetView<bool>("bo");
   // all 3 viewers populate the same page while none are being released.
   EXPECT_EQ(true, boView(1));
   EXPECT_EQ(false, boView2(2));
   EXPECT_EQ(true, boView3(3));

   auto stView = ntuple->GetView<std::string>("st");
   auto stView2 = ntuple->GetView<std::string>("st");
   auto stView3 = ntuple->GetView<std::string>("st");
   // all 3 viewers populate the same page while none are being released.
   EXPECT_EQ("foo1", stView(1));
   EXPECT_EQ("foo2", stView2(2));
   EXPECT_EQ("foo3", stView3(3));
}


TEST(RNTupleFriend, FriendOfFriendOfFriend)
{
   auto ntuple = RNTupleReader::Open(ntupleNameFriend, fileGuardFriend1.GetPath());
   auto ntuple2 = RNTupleReader::Open(ntupleNameFriend, fileGuardFriend2.GetPath());
   auto ntuple3 = RNTupleReader::Open(ntupleNameFriend, fileGuardFriend3.GetPath());
   auto ntuple4 = RNTupleReader::Open(ntupleNameFriend, fileGuardFriend4.GetPath());
   
   auto ntupleFriend = RNTupleReader::ChainReader(ntupleNameFriend, ntuple, ntuple2, EFileOpeningOptions::kFriend);
   auto ntupleFriendOfFriend = RNTupleReader::ChainReader(ntupleNameFriend, ntupleFriend, ntuple3, EFileOpeningOptions::kFriend);
   auto ntupleFriendOfFriendOfFriend = RNTupleReader::ChainReader(ntupleNameFriend, ntupleFriendOfFriend, ntuple4, EFileOpeningOptions::kFriend);
   EXPECT_EQ(ntuple->GetNEntries(), ntupleFriendOfFriendOfFriend->GetNEntries());
   
   auto ntupleNotFriend = RNTupleReader::Open(ntupleNameFriend1234Merged, fileGuardFriend1234Merged.GetPath());
   
   auto boView = ntupleFriendOfFriendOfFriend->GetView<bool>("bo");
   auto u64View = ntupleFriendOfFriendOfFriend->GetView<std::uint64_t>("u64");
   auto stView = ntupleFriendOfFriendOfFriend->GetView<std::string>("st");
   auto ftView = ntupleFriendOfFriendOfFriend->GetView<std::vector<float>>("ftVec");
   auto arView = ntupleFriendOfFriendOfFriend->GetView<std::array<std::int32_t, 4>>("ar32");
   auto bo2View = ntupleFriendOfFriendOfFriend->GetView<bool>("bo2");
   auto u642View = ntupleFriendOfFriendOfFriend->GetView<std::uint64_t>("u642");
   
   auto boView2 = ntupleNotFriend->GetView<bool>("bo");
   auto u64View2 = ntupleNotFriend->GetView<std::uint64_t>("u64");
   auto stView2 = ntupleNotFriend->GetView<std::string>("st");
   auto ftView2 = ntupleNotFriend->GetView<std::vector<float>>("ftVec");
   auto arView2 = ntupleNotFriend->GetView<std::array<std::int32_t, 4>>("ar32");
   auto bo2View2 = ntupleNotFriend->GetView<bool>("bo2");
   auto u642View2 = ntupleNotFriend->GetView<std::uint64_t>("u642");

   for (auto i : ntupleFriendOfFriendOfFriend->GetViewRange()) {
      EXPECT_EQ(boView2(i), boView(i));
      EXPECT_EQ(u64View2(i), u64View(i));
      EXPECT_EQ(stView2(i), stView(i));
      EXPECT_EQ(ftView2(i), ftView(i));
      EXPECT_EQ(arView2(i), arView(i));
      EXPECT_EQ(bo2View2(i), bo2View(i));
      EXPECT_EQ(u642View2(i), u642View(i));
   }
}


TEST(RNTupleFriend, FriendOfFriendOfFriendWithStdMove)
{
   auto ntuple = RNTupleReader::Open(ntupleNameFriend, fileGuardFriend1.GetPath());
   auto ntuple2 = RNTupleReader::Open(ntupleNameFriend, fileGuardFriend2.GetPath());
   auto ntuple3 = RNTupleReader::Open(ntupleNameFriend, fileGuardFriend3.GetPath());
   auto ntuple4 = RNTupleReader::Open(ntupleNameFriend, fileGuardFriend4.GetPath());
   
   auto ntupleFriend = RNTupleReader::ChainReader(ntupleNameFriend, std::move(ntuple), std::move(ntuple2), EFileOpeningOptions::kFriend);
   auto ntupleFriendOfFriend = RNTupleReader::ChainReader(ntupleNameFriend, std::move(ntupleFriend), std::move(ntuple3), EFileOpeningOptions::kFriend);
   auto ntupleFriendOfFriendOfFriend = RNTupleReader::ChainReader(ntupleNameFriend, std::move(ntupleFriendOfFriend), std::move(ntuple4), EFileOpeningOptions::kFriend);
   EXPECT_EQ(ntuple->GetNEntries(), ntupleFriendOfFriendOfFriend->GetNEntries());
   
   auto ntupleNotFriend = RNTupleReader::Open(ntupleNameFriend1234Merged, fileGuardFriend1234Merged.GetPath());
   
   auto boView = ntupleFriendOfFriendOfFriend->GetView<bool>("bo");
   auto u64View = ntupleFriendOfFriendOfFriend->GetView<std::uint64_t>("u64");
   auto stView = ntupleFriendOfFriendOfFriend->GetView<std::string>("st");
   auto ftView = ntupleFriendOfFriendOfFriend->GetView<std::vector<float>>("ftVec");
   auto arView = ntupleFriendOfFriendOfFriend->GetView<std::array<std::int32_t, 4>>("ar32");
   auto bo2View = ntupleFriendOfFriendOfFriend->GetView<bool>("bo2");
   auto u642View = ntupleFriendOfFriendOfFriend->GetView<std::uint64_t>("u642");
   
   auto boView2 = ntupleNotFriend->GetView<bool>("bo");
   auto u64View2 = ntupleNotFriend->GetView<std::uint64_t>("u64");
   auto stView2 = ntupleNotFriend->GetView<std::string>("st");
   auto ftView2 = ntupleNotFriend->GetView<std::vector<float>>("ftVec");
   auto arView2 = ntupleNotFriend->GetView<std::array<std::int32_t, 4>>("ar32");
   auto bo2View2 = ntupleNotFriend->GetView<bool>("bo2");
   auto u642View2 = ntupleNotFriend->GetView<std::uint64_t>("u642");

   for (auto i : ntupleFriendOfFriendOfFriend->GetViewRange()) {
      EXPECT_EQ(boView2(i), boView(i));
      EXPECT_EQ(u64View2(i), u64View(i));
      EXPECT_EQ(stView2(i), stView(i));
      EXPECT_EQ(ftView2(i), ftView(i));
      EXPECT_EQ(arView2(i), arView(i));
      EXPECT_EQ(bo2View2(i), bo2View(i));
      EXPECT_EQ(u642View2(i), u642View(i));
   }
}


TEST(RNTupleCombined, FriendOfChain)
{
   // Create: test_ntuple_mixed_File3.root
   {
      auto model = RNTupleModel::Create();
      auto blField = model->MakeField<bool>("bl");
      auto ntuple = RNTupleWriter::Recreate(std::move(model), ntupleNameMixed, fileGuardMixed3.GetPath());
      for (int j = 0; j < 30000; ++j) {
         *blField = true;
         ntuple->Fill();
      }
   } // flush contents to test_ntuple_mixed_File4.root
   
   // Create: test_ntuple_mixed_File3.root
   {
      auto model = RNTupleModel::Create();
      auto blField = model->MakeField<bool>("bl");
      auto ntuple = RNTupleWriter::Recreate(std::move(model), ntupleNameMixed, fileGuardMixed4.GetPath());
      for (int j = 0; j < 30000; ++j) {
         *blField = false;
         ntuple->Fill();
      }
   } // flush contents to test_ntuple_mixed_File4.root
   
   // Create: test_ntuple_mixed_File1234.root
   {
      auto model = RNTupleModel::Create();
      auto dbField = model->MakeField<double>("db");
      auto stField = model->MakeField<std::string>("st");
      auto itField = model->MakeField<std::vector<std::int32_t>>("it");
      auto blField = model->MakeField<bool>("bl");
      auto ntuple = RNTupleWriter::Recreate(std::move(model), ntupleNameMixed1234Merged, fileGuardMixed1234.GetPath());
      for (int j = 0; j < 30000; ++j) {
         *dbField = 5.0+j;
         *stField = "foo" + std::to_string(j);
         std::vector<std::int32_t> intVec(3);
         for (int k = 0; k < 3; ++k) {
            intVec.at(k) = k;
         }
         *itField = intVec;
         *blField = true;
         ntuple->Fill();
      }
      ntuple->CommitCluster();
      for (int j = 0; j < 30000; ++j) {
         *dbField = 10.11+j;
         *stField = "goo" + std::to_string(j);
         std::vector<std::int32_t> intVec(2);
         for (int k = 0; k < 2; ++k) {
            intVec.at(k) = k;
         }
         *itField = intVec;
         *blField = false;
         ntuple->Fill();
      }
   } // flush contents to test_ntuple_mixed_File1234.root
   
   auto ntupleChain1 = RNTupleReader::Open(ntupleNameMixed, { fileGuardMixed1.GetPath(), fileGuardMixed2.GetPath() });
   auto ntupleChain2 = RNTupleReader::Open(ntupleNameMixed, { fileGuardMixed3.GetPath(), fileGuardMixed4.GetPath() });
   auto ntupleFriendOfChain = RNTupleReader::ChainReader(ntupleNameMixed, ntupleChain1, ntupleChain2, EFileOpeningOptions::kFriend);
   auto ntupleNotChained = RNTupleReader::Open(ntupleNameMixed1234Merged, fileGuardMixed1234.GetPath());
   
   auto dbView = ntupleFriendOfChain->GetView<double>("db");
   auto stView = ntupleFriendOfChain->GetView<std::string>("st");
   auto itView = ntupleFriendOfChain->GetView<std::vector<std::int32_t>>("it");
   auto blView = ntupleFriendOfChain->GetView<bool>("bl");
   
   auto dbView2 = ntupleNotChained->GetView<double>("db");
   auto stView2 = ntupleNotChained->GetView<std::string>("st");
   auto itView2 = ntupleNotChained->GetView<std::vector<std::int32_t>>("it");
   auto blView2 = ntupleNotChained->GetView<bool>("bl");
   
   EXPECT_EQ(ntupleNotChained->GetNEntries(), ntupleFriendOfChain->GetNEntries());
   for (auto i : ntupleFriendOfChain->GetViewRange()) {
      EXPECT_DOUBLE_EQ(dbView2(i), dbView(i));
      EXPECT_EQ(stView2(i), stView(i));
      EXPECT_EQ(itView2(i), itView(i));
      EXPECT_EQ(blView2(i), blView(i));
   }
}


TEST(RNTupleCombined, ChainOfFriend)
{
   auto ntupleFriend1 = RNTupleReader::Open(ntupleNameMixed, { fileGuardMixed1.GetPath(), fileGuardMixed3.GetPath() }, EFileOpeningOptions::kFriend);
   auto ntupleFriend2 = RNTupleReader::Open(ntupleNameMixed, { fileGuardMixed2.GetPath(), fileGuardMixed4.GetPath() }, EFileOpeningOptions::kFriend);
   auto ntupleChainOfFriend = RNTupleReader::ChainReader(ntupleNameMixed, ntupleFriend1, ntupleFriend2);
   auto ntupleNotChained = RNTupleReader::Open(ntupleNameMixed1234Merged, fileGuardMixed1234.GetPath());
   
   auto dbView = ntupleChainOfFriend->GetView<double>("db");
   auto stView = ntupleChainOfFriend->GetView<std::string>("st");
   auto itView = ntupleChainOfFriend->GetView<std::vector<std::int32_t>>("it");
   auto blView = ntupleChainOfFriend->GetView<bool>("bl");
   
   auto dbView2 = ntupleNotChained->GetView<double>("db");
   auto stView2 = ntupleNotChained->GetView<std::string>("st");
   auto itView2 = ntupleNotChained->GetView<std::vector<std::int32_t>>("it");
   auto blView2 = ntupleNotChained->GetView<bool>("bl");
   
   EXPECT_EQ(ntupleNotChained->GetNEntries(), ntupleChainOfFriend->GetNEntries());
   for (auto i : ntupleChainOfFriend->GetViewRange()) {
      EXPECT_DOUBLE_EQ(dbView2(i), dbView(i));
      EXPECT_EQ(stView2(i), stView(i));
      EXPECT_EQ(itView2(i), itView(i));
      EXPECT_EQ(blView2(i), blView(i));
   }
}


#ifdef ROOT7_GETDESCRIPTOR_IN_RNTUPLEREADER
TEST(RNTupleCombined, CompareDescriptor)
{
   auto ntupleChain1 = RNTupleReader::Open(ntupleNameMixed, { fileGuardMixed1.GetPath(), fileGuardMixed2.GetPath() });
   auto ntupleChain2 = RNTupleReader::Open(ntupleNameMixed, { fileGuardMixed3.GetPath(), fileGuardMixed4.GetPath() });
   auto ntupleFriendOfChain = RNTupleReader::ChainReader(ntupleNameMixed, ntupleChain1, ntupleChain2, EFileOpeningOptions::kFriend);
   auto ntupleFriend1 = RNTupleReader::Open(ntupleNameMixed, { fileGuardMixed1.GetPath(), fileGuardMixed3.GetPath() }, EFileOpeningOptions::kFriend);
   auto ntupleFriend2 = RNTupleReader::Open(ntupleNameMixed, { fileGuardMixed2.GetPath(), fileGuardMixed4.GetPath() }, EFileOpeningOptions::kFriend);
   auto ntupleChainOfFriend = RNTupleReader::ChainReader(ntupleNameMixed, ntupleFriend1, ntupleFriend2);
   auto ntupleNotChained = RNTupleReader::Open(ntupleNameMixed1234Merged, fileGuardMixed1234.GetPath());
   
   EXPECT_TRUE(DescriptorsAreSame(ntupleNotChained->GetDescriptor(), ntupleFriendOfChain->GetDescriptor()));
   EXPECT_TRUE(DescriptorsAreSame(ntupleNotChained->GetDescriptor(), ntupleChainOfFriend->GetDescriptor()));
}
#endif
