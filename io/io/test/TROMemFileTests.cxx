#include "TError.h"
#include "TMemFile.h"
#include "TROOT.h"
#include "TTree.h"

#include <cstring>

#include "gtest/gtest.h"

static std::shared_ptr<const std::vector<char>> CreateBuffer(const char* title) {
   TNamed n("name", title);
   TMemFile memFile("a.root", "RECREATE", "TMemFile shared data test file", 0 /*no compression*/);
   memFile.WriteTObject(&n);
   memFile.Write();

   std::vector<char> data;
   data.resize(memFile.GetSize());
   memFile.CopyTo(&data.front(), data.size());
   return std::make_shared<const std::vector<char>>(std::move(data));
}

TEST(TROMemFile, Basics)
{
   constexpr const char title[] = "This is a title for TMemFile shared data test Basics";
   std::shared_ptr<const std::vector<char>> dataPtr(CreateBuffer(title));

   TMemFile rosmf("romemfile.root", dataPtr);
   TObject *readN = rosmf.Get("name");
   ASSERT_NE(nullptr, readN);
   EXPECT_STREQ(title, readN->GetTitle());
}

TEST(TROMemFile, NoWriting)
{
   constexpr const char title[] = "This is a title for TMemFile shared data test NoWriting";
   std::shared_ptr<const std::vector<char>> dataPtr(CreateBuffer(title));

   TMemFile rosmf("romemfile.root", dataPtr);
   auto oldIgnoreLevel = gErrorIgnoreLevel;
   gErrorIgnoreLevel = kBreak;
   TNamed doNotWrite("doNotWrite", "doNotWrite Title");
   EXPECT_EQ(0, rosmf.WriteTObject(&doNotWrite));
   gErrorIgnoreLevel = oldIgnoreLevel;
}

/// Check that TMemFile uses the original buffer, not a copy
TEST(TROMemFile, NoMemCopy)
{
   // Create a TNamed with this original title, and open the TMemFile with it.
   constexpr const char title1[] = "This is a title for TMemFile shared data test NoMemCopy";
   std::shared_ptr<const std::vector<char>> dataPtr1(CreateBuffer(title1));
   TMemFile rosmf("romemfile.root", dataPtr1);

   // Swap rosmf's data buffer against another one, with a different title for the TNamed.
   constexpr const char title2[] = "Fish is a title for TMemFile shared data test NoMemCopy";
   std::shared_ptr<const std::vector<char>> dataPtr2(CreateBuffer(title2));
   std::vector<char> &dataVec1 = const_cast<std::vector<char>&>(*dataPtr1);
   std::vector<char> &dataVec2 = const_cast<std::vector<char>&>(*dataPtr2);
   std::copy(dataVec2.begin(), dataVec2.end(), dataVec1.begin());

   /// Make sure rosmf sees the changed buffer, because it doesn't copy the buffer:
   TObject *readN = rosmf.Get("name");
   ASSERT_NE(nullptr, readN);
   EXPECT_STREQ(title2, readN->GetTitle());
}

/// Check that the TMemFile is not registered in the list of files
TEST(TROMemFile, ListOfFilesRegistration)
{
   auto &lof = *gROOT->GetListOfFiles();

   TMemFile registered("registered.root","RECREATE");
   EXPECT_EQ(lof.GetSize(), 1);

   TMemFile notregistered("notregistered.root","RECREATEUNMANAGED");
   EXPECT_EQ(lof.GetSize(), 1);
}
