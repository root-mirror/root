#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "ROOT/RMakeUnique.hxx"

#include <array>
#include <memory>

#include "gtest/gtest.h"

std::unique_ptr<TTree> MakeTestTree()
{
   auto tree = std::make_unique<TTree>("TTreeReaderArrayTree", "In-memory test tree");
   std::array<int, 6> values{{0, 1, 2, 3, 4, 5}};
   int n = 6;
   tree->Branch("n", &n);
   tree->Branch("arr", &values, "arr[n]/I");
   tree->Fill();
   n = 0;
   tree->Fill();
   tree->ResetBranchAddresses();
   return tree;
}

struct NonConstTag {};
struct ConstTag {};

template <typename ConstnessTag>
void TestReaderArray()
{
   static_assert(std::is_same<ConstnessTag, ConstTag>::value || std::is_same<ConstnessTag, NonConstTag>::value, "");
   using ReaderArray_t = typename std::conditional<std::is_same<ConstnessTag, NonConstTag>::value,
                                                   TTreeReaderArray<int>, const TTreeReaderArray<int>>::type;

   // create TTreeReaderArray
   auto tree = MakeTestTree();
   TTreeReader r(tree.get());
   ReaderArray_t arr(r, "arr");

   // check iterator type (iterator vs const_iterator)
   using ExpectedIterator_t =
      typename std::conditional<std::is_same<ConstnessTag, NonConstTag>::value,
                                TTreeReaderArray<int>::Iterator_t<TTreeReaderArray<int>>,
                                TTreeReaderArray<int>::Iterator_t<const TTreeReaderArray<int>>>::type;
   ::testing::StaticAssertTypeEq<decltype(arr.begin()), ExpectedIterator_t>();

   // load full entry
   r.Next();

   // begin and end, comparisons
   EXPECT_EQ(arr.begin(), arr.begin());
   EXPECT_EQ(arr.end(), arr.end());
   EXPECT_NE(arr.begin(), arr.end());

   // copy/move ctors
   auto it = arr.begin();
   const auto it2(it);
   (void)it2;
   auto it3(std::move(it));

   // increments, decrements, copy/move assignments
   EXPECT_EQ(++it, arr.begin() + 1);
   it = it3++;
   EXPECT_EQ(it += 1, it3);

   EXPECT_EQ((std::ptrdiff_t(2) + arr.begin()) - 1, --(arr.begin() + 2));
   it = std::move(it3);
   auto it4 = it--;
   EXPECT_EQ(it4 -= 1, it);

   // comparisons
   auto e = arr.end();
   auto b = arr.begin();
   EXPECT_GT(e, b);
   EXPECT_LT(b, b + 2);
   EXPECT_LE(b + 1, b + 1);
   EXPECT_GE(b + 1, b + 1);

   // difference
   const std::ptrdiff_t length = arr.end() - arr.begin();
   const auto size = arr.GetSize();
   EXPECT_EQ(static_cast<decltype(size)>(length), size);
   EXPECT_EQ(arr.begin() + length, arr.end());
   const auto it5 = arr.begin() + 2;
   EXPECT_EQ(arr.begin() + (it5 - arr.begin()), arr.begin() + 2);

   // operator[]
   EXPECT_EQ(arr.begin()[2u], *(arr.begin() + 2));

   // ranged for loop
   int truth = 0;
   for (auto v : arr)
      EXPECT_DOUBLE_EQ(v, truth++);

   // load empty entry
   r.Next();
   EXPECT_EQ(arr.begin(), arr.end());

   // end() < end()
   EXPECT_EQ(arr.end(), arr.end());
   EXPECT_FALSE(arr.end() < arr.end() || arr.end() > arr.end());
   EXPECT_FALSE(++arr.end() > arr.end());
}

TEST(ReaderArrayIterator, NonConst)
{
   TestReaderArray<NonConstTag>();
}

TEST(ReaderArrayIterator, ModifyContent)
{
   auto tree = MakeTestTree();
   TTreeReader r(tree.get());
   TTreeReaderArray<int> arr(r, "arr");
   r.Next();

   for (auto it = arr.begin(); it != arr.end(); ++it)
      *it = 42;

   auto it = arr.cbegin();
   const auto size = arr.GetSize();
   for (auto i = 0u; i < size; ++i)
      EXPECT_EQ(it[i], 42);
}

TEST(ReaderArrayIterator, Const)
{
   TestReaderArray<ConstTag>();
}