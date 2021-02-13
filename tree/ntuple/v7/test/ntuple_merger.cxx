#include "ntuple_test.hxx"

TEST(RFieldMerger, Merge)
{
   auto mergeResult = RFieldMerger::Merge(RFieldDescriptor(), RFieldDescriptor());
   EXPECT_FALSE(mergeResult);
}
