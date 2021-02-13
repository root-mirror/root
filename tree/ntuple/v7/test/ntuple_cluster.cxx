#include "gtest/gtest.h"

#include <ROOT/RCluster.hxx>
#include <ROOT/RClusterPool.hxx>
#include <ROOT/RColumn.hxx>
#include <ROOT/RColumnModel.hxx>
#include <ROOT/RNTuple.hxx>
#include <ROOT/RNTupleDescriptor.hxx>
#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RNTupleOptions.hxx>
#include <ROOT/RNTupleUtil.hxx>
#include <ROOT/RPage.hxx>
#include <ROOT/RPageStorage.hxx>
#include <ROOT/RPageStorageFile.hxx>
#include <ROOT/RStringView.hxx>

#include <memory>
#include <utility>
#include <vector>

using ClusterSize_t = ROOT::Experimental::ClusterSize_t;
using RCluster = ROOT::Experimental::Detail::RCluster;
using RClusterPool = ROOT::Experimental::Detail::RClusterPool;
using RNTupleDescriptor = ROOT::Experimental::RNTupleDescriptor;
using RNTupleVersion = ROOT::Experimental::RNTupleVersion;
using ROnDiskPage = ROOT::Experimental::Detail::ROnDiskPage;
using RPage = ROOT::Experimental::Detail::RPage;
using RPageSource = ROOT::Experimental::Detail::RPageSource;

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

/**
 * Used to track LoadCluster calls triggered by ClusterPool::GetCluster
 */
class RPageSourceMock : public RPageSource {
protected:
   RNTupleDescriptor AttachImpl() final { return RNTupleDescriptor(); }

public:
   /// Records the cluster IDs requests by LoadCluster() calls
   std::vector<ROOT::Experimental::DescriptorId_t> fReqsClusterIds;
   std::vector<ROOT::Experimental::Detail::RPageSource::ColumnSet_t> fReqsColumns;

   RPageSourceMock() : RPageSource("test", ROOT::Experimental::RNTupleReadOptions()) {
      ROOT::Experimental::RNTupleDescriptorBuilder descBuilder;
      descBuilder.AddCluster(0, RNTupleVersion(), 0, ClusterSize_t(1));
      descBuilder.AddCluster(1, RNTupleVersion(), 1, ClusterSize_t(1));
      descBuilder.AddCluster(2, RNTupleVersion(), 2, ClusterSize_t(1));
      descBuilder.AddCluster(3, RNTupleVersion(), 3, ClusterSize_t(1));
      descBuilder.AddCluster(4, RNTupleVersion(), 4, ClusterSize_t(1));
      fDescriptor = descBuilder.MoveDescriptor();
   }
   std::unique_ptr<RPageSource> Clone() const final { return nullptr; }
   RPage PopulatePage(ColumnHandle_t, ROOT::Experimental::NTupleSize_t) final { return RPage(); }
   RPage PopulatePage(ColumnHandle_t, const ROOT::Experimental::RClusterIndex &) final { return RPage(); }
   void ReleasePage(RPage &) final {}
   std::unique_ptr<RCluster> LoadCluster(
      ROOT::Experimental::DescriptorId_t clusterId,
      const ROOT::Experimental::Detail::RPageSource::ColumnSet_t &columns) final
   {
      fReqsClusterIds.emplace_back(clusterId);
      fReqsColumns.emplace_back(columns);
      auto cluster = std::make_unique<RCluster>(clusterId);
      auto pageMap = std::make_unique<ROOT::Experimental::Detail::ROnDiskPageMap>();
      for (auto colId : columns) {
         pageMap->Register(ROnDiskPage::Key(colId, 0), ROnDiskPage(nullptr, 0));
         cluster->SetColumnAvailable(colId);
      }
      cluster->Adopt(std::move(pageMap));
      return cluster;
   }
};

} // anonymous namespace


TEST(Cluster, Allocate)
{
   auto cluster = new ROOT::Experimental::Detail::ROnDiskPageMapHeap(nullptr);
   delete cluster;

   cluster = new ROOT::Experimental::Detail::ROnDiskPageMapHeap(std::make_unique<unsigned char []>(1));
   delete cluster;
}


TEST(Cluster, Basics)
{
   auto memory = new unsigned char[3];
   auto pageMap = std::make_unique<ROOT::Experimental::Detail::ROnDiskPageMapHeap>(
      std::unique_ptr<unsigned char []>(memory));
   pageMap->Register(ROnDiskPage::Key(5, 0), ROnDiskPage(&memory[0], 1));
   pageMap->Register(ROnDiskPage::Key(5, 1), ROnDiskPage(&memory[1], 2));
   auto cluster = std::make_unique<RCluster>(0);
   cluster->Adopt(std::move(pageMap));
   cluster->SetColumnAvailable(5);

   EXPECT_EQ(nullptr, cluster->GetOnDiskPage(ROnDiskPage::Key(5, 2)));
   EXPECT_EQ(nullptr, cluster->GetOnDiskPage(ROnDiskPage::Key(4, 0)));
   auto onDiskPage = cluster->GetOnDiskPage(ROnDiskPage::Key(5, 0));
   EXPECT_EQ(&memory[0], onDiskPage->GetAddress());
   EXPECT_EQ(1U, onDiskPage->GetSize());
   onDiskPage = cluster->GetOnDiskPage(ROnDiskPage::Key(5, 1));
   EXPECT_EQ(&memory[1], onDiskPage->GetAddress());
   EXPECT_EQ(2U, onDiskPage->GetSize());
}

TEST(Cluster, AdoptPageMaps)
{
   auto mem1 = new unsigned char[3];
   auto pageMap1 = std::make_unique<ROOT::Experimental::Detail::ROnDiskPageMapHeap>(
      std::unique_ptr<unsigned char []>(mem1));
   pageMap1->Register(ROnDiskPage::Key(5, 0), ROnDiskPage(&mem1[0], 1));
   pageMap1->Register(ROnDiskPage::Key(5, 1), ROnDiskPage(&mem1[1], 2));
   // Column 5 is in both mem1 and mem2 but that should not hurt
   auto mem2 = new unsigned char[4];
   auto pageMap2 = std::make_unique<ROOT::Experimental::Detail::ROnDiskPageMapHeap>(
      std::unique_ptr<unsigned char []>(mem2));
   pageMap2->Register(ROnDiskPage::Key(5, 0), ROnDiskPage(&mem2[0], 1));
   pageMap2->Register(ROnDiskPage::Key(5, 1), ROnDiskPage(&mem2[1], 2));
   pageMap2->Register(ROnDiskPage::Key(6, 0), ROnDiskPage(&mem2[3], 1));

   auto cluster = std::make_unique<RCluster>(0);
   cluster->Adopt(std::move(pageMap1));
   cluster->Adopt(std::move(pageMap2));
   cluster->SetColumnAvailable(5);
   cluster->SetColumnAvailable(6);

   EXPECT_EQ(3U, cluster->GetNOnDiskPages());
   EXPECT_TRUE(cluster->ContainsColumn(5));
   EXPECT_TRUE(cluster->ContainsColumn(6));
   EXPECT_FALSE(cluster->ContainsColumn(7));

   auto onDiskPage = cluster->GetOnDiskPage(ROnDiskPage::Key(6, 0));
   EXPECT_EQ(&mem2[3], onDiskPage->GetAddress());
   EXPECT_EQ(1U, onDiskPage->GetSize());
   onDiskPage = cluster->GetOnDiskPage(ROnDiskPage::Key(5, 0));
   EXPECT_TRUE((onDiskPage->GetAddress() == &mem1[0]) || (onDiskPage->GetAddress() == &mem2[0]));
   EXPECT_EQ(1U, onDiskPage->GetSize());
   onDiskPage = cluster->GetOnDiskPage(ROnDiskPage::Key(5, 1));
   EXPECT_TRUE((onDiskPage->GetAddress() == &mem1[1]) || (onDiskPage->GetAddress() == &mem2[1]));
   EXPECT_EQ(2U, onDiskPage->GetSize());
}


TEST(Cluster, AdoptClusters)
{
   auto mem1 = new unsigned char[3];
   auto pageMap1 = std::make_unique<ROOT::Experimental::Detail::ROnDiskPageMapHeap>(
      std::unique_ptr<unsigned char []>(mem1));
   pageMap1->Register(ROnDiskPage::Key(5, 0), ROnDiskPage(&mem1[0], 1));
   pageMap1->Register(ROnDiskPage::Key(5, 1), ROnDiskPage(&mem1[1], 2));
   auto cluster1 = std::make_unique<RCluster>(0);
   cluster1->Adopt(std::move(pageMap1));
   cluster1->SetColumnAvailable(5);

   // Column 5 is in both clusters but that should not hurt
   auto mem2 = new unsigned char[4];
   auto pageMap2 = std::make_unique<ROOT::Experimental::Detail::ROnDiskPageMapHeap>(
      std::unique_ptr<unsigned char []>(mem2));
   pageMap2->Register(ROnDiskPage::Key(5, 0), ROnDiskPage(&mem2[0], 1));
   pageMap2->Register(ROnDiskPage::Key(5, 1), ROnDiskPage(&mem2[1], 2));
   pageMap2->Register(ROnDiskPage::Key(6, 0), ROnDiskPage(&mem2[3], 1));
   auto cluster2 = std::make_unique<RCluster>(0);
   cluster2->Adopt(std::move(pageMap2));
   cluster2->SetColumnAvailable(5);
   cluster2->SetColumnAvailable(6);

   cluster2->Adopt(std::move(*cluster1));

   EXPECT_EQ(3U, cluster2->GetNOnDiskPages());
   EXPECT_TRUE(cluster2->ContainsColumn(5));
   EXPECT_TRUE(cluster2->ContainsColumn(6));
   EXPECT_FALSE(cluster2->ContainsColumn(7));

   auto onDiskPage = cluster2->GetOnDiskPage(ROnDiskPage::Key(6, 0));
   EXPECT_EQ(&mem2[3], onDiskPage->GetAddress());
   EXPECT_EQ(1U, onDiskPage->GetSize());
   onDiskPage = cluster2->GetOnDiskPage(ROnDiskPage::Key(5, 0));
   EXPECT_TRUE((onDiskPage->GetAddress() == &mem1[0]) || (onDiskPage->GetAddress() == &mem2[0]));
   EXPECT_EQ(1U, onDiskPage->GetSize());
   onDiskPage = cluster2->GetOnDiskPage(ROnDiskPage::Key(5, 1));
   EXPECT_TRUE((onDiskPage->GetAddress() == &mem1[1]) || (onDiskPage->GetAddress() == &mem2[1]));
   EXPECT_EQ(2U, onDiskPage->GetSize());
}


TEST(ClusterPool, Windows)
{
   RPageSourceMock ps;

   EXPECT_DEATH(RClusterPool(ps, 0), ".*");
   RClusterPool c1(ps, 1);
   EXPECT_EQ(0U, c1.GetWindowPre());
   EXPECT_EQ(1U, c1.GetWindowPost());
   RClusterPool c2(ps, 2);
   EXPECT_EQ(0U, c2.GetWindowPre());
   EXPECT_EQ(2U, c2.GetWindowPost());
   RClusterPool c3(ps, 3);
   EXPECT_EQ(1U, c3.GetWindowPre());
   EXPECT_EQ(2U, c3.GetWindowPost());
   RClusterPool c5(ps, 5);
   EXPECT_EQ(1U, c5.GetWindowPre());
   EXPECT_EQ(4U, c5.GetWindowPost());
   RClusterPool c6(ps, 6);
   EXPECT_EQ(2U, c6.GetWindowPre());
   EXPECT_EQ(4U, c6.GetWindowPost());
   RClusterPool c9(ps, 9);
   EXPECT_EQ(2U, c9.GetWindowPre());
   EXPECT_EQ(7U, c9.GetWindowPost());
   RClusterPool c10(ps, 10);
   EXPECT_EQ(3U, c10.GetWindowPre());
   EXPECT_EQ(7U, c10.GetWindowPost());
   RClusterPool c15(ps, 15);
   EXPECT_EQ(3U,  c15.GetWindowPre());
   EXPECT_EQ(12U, c15.GetWindowPost());
   RClusterPool c16(ps, 16);
   EXPECT_EQ(4U,  c16.GetWindowPre());
   EXPECT_EQ(12U, c16.GetWindowPost());
}

TEST(ClusterPool, GetClusterBasics)
{
   RPageSourceMock p1;
   RClusterPool c1(p1, 1);
   c1.GetCluster(3, {0});
   ASSERT_EQ(1U, p1.fReqsClusterIds.size());
   EXPECT_EQ(3U, p1.fReqsClusterIds[0]);
   EXPECT_EQ(RPageSource::ColumnSet_t({0}), p1.fReqsColumns[0]);

   RPageSourceMock p2;
   {
      RClusterPool c2(p2, 2);
      c2.GetCluster(0, {0});
   }
   ASSERT_EQ(2U, p2.fReqsClusterIds.size());
   EXPECT_EQ(0U, p2.fReqsClusterIds[0]);
   EXPECT_EQ(1U, p2.fReqsClusterIds[1]);
   EXPECT_EQ(RPageSource::ColumnSet_t({0}), p2.fReqsColumns[0]);
   EXPECT_EQ(RPageSource::ColumnSet_t({0}), p2.fReqsColumns[1]);

   RPageSourceMock p3;
   {
      RClusterPool c3(p3, 4);
      c3.GetCluster(2, {0});
   }
   ASSERT_EQ(3U, p3.fReqsClusterIds.size());
   EXPECT_EQ(2U, p3.fReqsClusterIds[0]);
   EXPECT_EQ(3U, p3.fReqsClusterIds[1]);
   EXPECT_EQ(4U, p3.fReqsClusterIds[2]);
   EXPECT_EQ(RPageSource::ColumnSet_t({0}), p3.fReqsColumns[0]);
   EXPECT_EQ(RPageSource::ColumnSet_t({0}), p3.fReqsColumns[1]);
   EXPECT_EQ(RPageSource::ColumnSet_t({0}), p3.fReqsColumns[2]);
}


TEST(ClusterPool, GetClusterIncrementally)
{
   RPageSourceMock p1;
   RClusterPool c1(p1, 1);
   c1.GetCluster(3, {0});
   ASSERT_EQ(1U, p1.fReqsClusterIds.size());
   EXPECT_EQ(3U, p1.fReqsClusterIds[0]);
   EXPECT_EQ(RPageSource::ColumnSet_t({0}), p1.fReqsColumns[0]);

   c1.GetCluster(3, {1});
   ASSERT_EQ(2U, p1.fReqsClusterIds.size());
   EXPECT_EQ(3U, p1.fReqsClusterIds[1]);
   EXPECT_EQ(RPageSource::ColumnSet_t({1}), p1.fReqsColumns[1]);
}


TEST(PageStorageFile, LoadCluster)
{
   FileRaii fileGuard("test_ntuple_clusters.root");

   auto modelWrite = ROOT::Experimental::RNTupleModel::Create();
   auto wrPt = modelWrite->MakeField<float>("pt", 42.0);
   auto wrTag = modelWrite->MakeField<int32_t>("tag", 0);

   {
      ROOT::Experimental::RNTupleWriter ntuple(
         std::move(modelWrite), std::make_unique<ROOT::Experimental::Detail::RPageSinkFile>(
            "myNTuple", fileGuard.GetPath(), ROOT::Experimental::RNTupleWriteOptions()));
      ntuple.Fill();
      ntuple.CommitCluster();
      *wrPt = 24.0;
      *wrTag = 1;
      ntuple.Fill();
   }

   ROOT::Experimental::Detail::RPageSourceFile source(
      "myNTuple", fileGuard.GetPath(), ROOT::Experimental::RNTupleReadOptions());
   source.Attach();

   auto ptId = source.GetDescriptor().FindFieldId("pt");
   EXPECT_NE(ROOT::Experimental::kInvalidDescriptorId, ptId);
   auto colId = source.GetDescriptor().FindColumnId(ptId, 0);
   EXPECT_NE(ROOT::Experimental::kInvalidDescriptorId, colId);

   auto cluster = source.LoadCluster(0, {});
   EXPECT_EQ(0U, cluster->GetId());
   EXPECT_EQ(0U, cluster->GetNOnDiskPages());

   auto column = std::unique_ptr<ROOT::Experimental::Detail::RColumn>(
      ROOT::Experimental::Detail::RColumn::Create<float, ROOT::Experimental::EColumnType::kReal32>(
         ROOT::Experimental::RColumnModel(ROOT::Experimental::EColumnType::kReal32, false), 0));
   column->Connect(ptId, &source);
   cluster = source.LoadCluster(1, {colId});
   EXPECT_EQ(1U, cluster->GetId());
   EXPECT_EQ(1U, cluster->GetNOnDiskPages());

   ROnDiskPage::Key key(colId, 0);
   EXPECT_NE(nullptr, cluster->GetOnDiskPage(key));
}
