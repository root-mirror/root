#include "ROOTUnitTestSupport.h"

#include "TFileMerger.h"

#include "TMemFile.h"
#include "TTree.h"

static void CreateATuple(TMemFile &file, const char *name, double value)
{
   auto mytree = new TTree(name, "A tree");
   // FIXME: We inherit EnableImplicitIMT from TBufferMerger tests (we are sharing the same executable) where we call
   // EnableThreadSafety(). Here, we hit a race condition in TBranch::FlushBaskets. Once we get that fixed we probably
   // should re-enable implicit MT.
   //
   // In general, we should probably have a way to conditionally enable/disable thread safety.
   mytree->SetImplicitMT(false);

   mytree->SetDirectory(&file);
   mytree->Branch(name, &value);
   mytree->Fill();
   file.Write();
}

static void CheckTree(TMemFile &file, const char *name, double expectedValue)
{
   auto t = static_cast<TTree *>(file.Get(name));
   ASSERT_TRUE(t != nullptr);

   double d;
   t->SetBranchAddress(name, &d);
   t->GetEntry(0);
   EXPECT_EQ(expectedValue, d);
   // Setting branch address to a stack address requires to either call ResetBranchAddresses or delete the tree.
   t->ResetBranchAddresses();
}

TEST(TFileMerger, CreateWithTFilePointer)
{
   TMemFile a("a.root", "RECREATE");
   CreateATuple(a, "a_tree", 1.);

   TMemFile b("b.root", "RECREATE");
   CreateATuple(b, "b_tree", 2.);

   TFileMerger merger;
   auto output = std::unique_ptr<TMemFile>(new TMemFile("output.root", "CREATE"));
   bool success = merger.OutputFile(std::move(output));

   ASSERT_TRUE(success);

   merger.AddFile(&a, false);
   merger.AddFile(&b, false);
   // FIXME: Calling merger.Merge() will call Close() and *delete* output.
   merger.PartialMerge();

   auto &result = *static_cast<TMemFile *>(merger.GetOutputFile());
   CheckTree(result, "a_tree", 1);
   CheckTree(result, "b_tree", 2);
}

TEST(TFileMerger, CreateWithUnwritableTFilePointer)
{
   TFileMerger merger;
   auto output = std::unique_ptr<TMemFile>(new TMemFile("output.root", "RECREATE"));
   // FIXME: The ctor of TMemFile sets the 'zombie' flag to all TMemFiles whose options are different than CREATE and
   // RECREATE. We should probably fix the API but until then work around it.
   output->SetWritable(false);
   ROOT_EXPECT_ERROR(merger.OutputFile(std::move(output)), "TFileMerger::OutputFile",
                     "output file output.root is not writable");
}
