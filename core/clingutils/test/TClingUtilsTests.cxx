/// \file TClingUtilsTests.cxx
///
/// \brief The file contain unit tests which test the TClingUtils.h
///
/// \author Vassil Vassilev <vvasilev@cern.ch>
///
/// \date Aug, 2019
///
/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <TClingUtils.h>

#include <ROOT/FoundationUtils.hxx>

#include "gtest/gtest.h"

#include <fstream>

TEST(TClingUtilsTests, GetRealPath)
{
#ifndef R__WIN32
   // Create a file in the current folder
   std::string cwd = ROOT::FoundationUtils::GetCurrentDir();

   ASSERT_FALSE(cwd.empty());

   std::ofstream file;
   file.open("./realfile1", std::ios::out);
   file.close();
   file.open("./realfile2", std::ios::out);
   file.close();
   // Create symlinks
   ASSERT_GE(symlink("./realfile1", "symlink_realfile1"), 0);
   ASSERT_GE(symlink((cwd + "realfile2").c_str(), "symlink_realfile2"), 0);

   // Recursive symlinks
   ASSERT_GE(symlink("./symlink_realfile1", "symlink1_symlink_realfile1"), 0);
   ASSERT_GE(symlink((cwd + "symlink_realfile2").c_str(), "symlink1_symlink_realfile2"), 0);

   ASSERT_GE(symlink("./symlink_realfile2", "symlink2_symlink_realfile2"), 0);
   ASSERT_GE(symlink((cwd + "symlink_realfile1").c_str(), "symlink2_symlink_realfile1"), 0);

   ASSERT_GE(symlink((cwd + ".//.///symlink_realfile1").c_str(), "symlink3_symlink_realfile1"), 0);

   using namespace ROOT::TMetaUtils;
   std::string realfile1 = GetRealPath("./realfile1");
   std::string realfile2 = GetRealPath("./realfile2");
   ASSERT_EQ(realfile1, GetRealPath("./symlink_realfile1"));
   ASSERT_EQ(realfile1, GetRealPath("./symlink_realfile1"));
   ASSERT_EQ(realfile1, GetRealPath("./symlink1_symlink_realfile1"));
   ASSERT_EQ(realfile1, GetRealPath("./symlink2_symlink_realfile1"));
   ASSERT_EQ(realfile1, GetRealPath("./symlink3_symlink_realfile1"));

   ASSERT_EQ(realfile2, GetRealPath("./realfile2"));
   ASSERT_EQ(realfile2, GetRealPath("./symlink_realfile2"));
   ASSERT_EQ(realfile2, GetRealPath("./symlink1_symlink_realfile2"));
   ASSERT_EQ(realfile2, GetRealPath("./symlink2_symlink_realfile2"));

   std::remove("./symlink3_symlink_realfile1");
   std::remove("./symlink2_symlink_realfile2");
   std::remove("./symlink2_symlink_realfile1");
   std::remove("./symlink1_symlink_realfile2");
   std::remove("./symlink1_symlink_realfile1");
   std::remove("./symlink_realfile1");
   std::remove("./symlink_realfile2");
   std::remove("./realfile1");
   std::remove("./realfile2");
#endif // not R__WIN32
}
