#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TH1D.h>

#include "gtest/gtest.h"

// ROOT-10023
TEST(TTreeReaderBasic, BranchObject)
{
   TTree t("t", "t");
   TObject *o = nullptr;
   t.Branch("o", &o, 32000, 0); // must be unsplit to generate a TBranchObject

   // Fill branch with different concrete types
   TNamed name("name", "title");
   TList list;
   list.Add(&name);
   o = &list;
   t.Fill();
   TH1D h("h", "h", 100, 0, 100);
   h.Fill(42);
   o = &h;
   t.Fill();
   o = nullptr;

   TTreeReader r(&t);
   TTreeReaderValue<TObject> rv(r, "o");
   EXPECT_TRUE(r.Next());
   EXPECT_STREQ(rv->ClassName(), "TList");
   auto *asList = static_cast<TList *>(rv.Get());
   EXPECT_EQ(asList->GetEntries(), 1);
   EXPECT_STREQ(asList->At(0)->GetTitle(), "title");
   EXPECT_STREQ(asList->At(0)->GetName(), "name");

   EXPECT_TRUE(r.Next());
   EXPECT_STREQ(rv->ClassName(), "TH1D");
   EXPECT_DOUBLE_EQ(static_cast<TH1D *>(rv.Get())->GetMean(), 42.);
}
