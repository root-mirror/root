#include <TChain.h>
#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>
 
#include "gtest/gtest.h"

class TTreeCache;

// ROOT-10672
TEST(TChain, GetReadCacheBug)
{
   const auto treename = "tree";
   const auto filename = "tchain_getreadcachebug.root";
   {
      TFile f(filename, "recreate");
      ASSERT_FALSE(f.IsZombie());
      TTree t(treename, treename);
      t.Fill();
      t.Write();
      f.Close();
   }

   TChain chain(treename);
   chain.Add(filename);
   chain.GetEntry(0);
   TFile *currentfile = chain.GetCurrentFile();
   TTreeCache *treecache = chain.GetReadCache(currentfile, true);
   EXPECT_NE(treecache, nullptr);

   gSystem->Unlink(filename);
}
