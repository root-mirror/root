#include "ntuple_test.hxx"

TEST(Pages, Allocation)
{
   RPageAllocatorHeap allocator;

   auto page = allocator.NewPage(42, 4, 16);
   EXPECT_FALSE(page.IsNull());
   EXPECT_EQ(64U, page.GetCapacity());
   EXPECT_EQ(0U, page.GetNElements());
   EXPECT_EQ(0U, page.GetSize());
   allocator.DeletePage(page);
}

TEST(Pages, Pool)
{
   RPagePool pool;

   auto page = pool.GetPage(0, 0);
   EXPECT_TRUE(page.IsNull());
   pool.ReturnPage(page); // should not crash

   RPage::RClusterInfo clusterInfo(2, 40);
   page = RPage(1, &page, 10, 1);
   EXPECT_NE(nullptr, page.TryGrow(10));
   page.SetWindow(50, clusterInfo);
   EXPECT_FALSE(page.IsNull());
   unsigned int nCallDeleter = 0;
   pool.RegisterPage(page, RPageDeleter([&nCallDeleter](const RPage & /*page*/, void * /*userData*/) {
      nCallDeleter++;
   }));

   page = pool.GetPage(0, 0);
   EXPECT_TRUE(page.IsNull());
   page = pool.GetPage(0, 55);
   EXPECT_TRUE(page.IsNull());
   page = pool.GetPage(1, 55);
   EXPECT_FALSE(page.IsNull());
   EXPECT_EQ(50U, page.GetGlobalRangeFirst());
   EXPECT_EQ(59U, page.GetGlobalRangeLast());
   EXPECT_EQ(10U, page.GetClusterRangeFirst());
   EXPECT_EQ(19U, page.GetClusterRangeLast());

   page = pool.GetPage(1, ROOT::Experimental::RClusterIndex(0, 15));
   EXPECT_TRUE(page.IsNull());
   page = pool.GetPage(1, ROOT::Experimental::RClusterIndex(2, 15));
   EXPECT_FALSE(page.IsNull());

   pool.ReturnPage(page);
   EXPECT_EQ(0U, nCallDeleter);
   pool.ReturnPage(page);
   EXPECT_EQ(0U, nCallDeleter);
   pool.ReturnPage(page);
   EXPECT_EQ(1U, nCallDeleter);
   page = pool.GetPage(1, 55);
   EXPECT_TRUE(page.IsNull());
}
