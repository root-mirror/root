#include "TChain.h"
#include "gtest/gtest.h"

TEST(TChainParsing, RemoteAdd)
{
   TChain c;
   c.Add("root://some.domain/path/to/file.root/treename");
   c.Add("root://some.domain//path/to/file.root/treename");
   c.Add("root://some.domain/path/to/foo.something/file.root/treename");
   c.Add("root://some.domain/path/to/foo.root/file.root/treename"); // ROOT-9344
   c.Add("root://some.domain/path//to/file.root/treename"); // ROOT-10494
   c.Add("root://some.domain//path//to//file.root//treename");
   const auto files = c.GetListOfFiles();

   EXPECT_STREQ(files->At(0)->GetTitle(), "root://some.domain/path/to/file.root");
   EXPECT_STREQ(files->At(0)->GetName(), "treename");

   EXPECT_STREQ(files->At(1)->GetTitle(), "root://some.domain//path/to/file.root");
   EXPECT_STREQ(files->At(1)->GetName(), "treename");

   EXPECT_STREQ(files->At(2)->GetTitle(), "root://some.domain/path/to/foo.something/file.root");
   EXPECT_STREQ(files->At(2)->GetName(), "treename");

   EXPECT_STREQ(files->At(3)->GetTitle(), "root://some.domain/path/to/foo.root/file.root");
   EXPECT_STREQ(files->At(3)->GetName(), "treename");

   EXPECT_STREQ(files->At(4)->GetTitle(), "root://some.domain/path/to/file.root");
   EXPECT_STREQ(files->At(4)->GetName(), "treename");

   EXPECT_STREQ(files->At(5)->GetTitle(), "root://some.domain//path/to/file.root");
   EXPECT_STREQ(files->At(5)->GetName(), "treename");
}

TEST(TChainParsing, LocalAdd)
{
   TChain c;
   c.Add("/path/to/file.root");
   c.Add("/path/to/file.root/foo");
   c.Add("/path/to/file.root/foo/bar");
   c.Add("/path/to/file.root/foo.bar/treename");
   c.Add("/path/to/file.root/foo.root/treename");
   c.Add("/path/to/file.root/root/treename");
   c.Add("path/to/file.root/treename");
   c.Add("/path/to/file.root//treename");
   const auto files = c.GetListOfFiles();

   EXPECT_STREQ(files->At(0)->GetTitle(), "/path/to/file.root");
   EXPECT_STREQ(files->At(0)->GetName(), "");

   EXPECT_STREQ(files->At(1)->GetTitle(), "/path/to/file.root");
   EXPECT_STREQ(files->At(1)->GetName(), "foo");

   EXPECT_STREQ(files->At(2)->GetTitle(), "/path/to/file.root");
   EXPECT_STREQ(files->At(2)->GetName(), "foo/bar");

   EXPECT_STREQ(files->At(3)->GetTitle(), "/path/to/file.root");
   EXPECT_STREQ(files->At(3)->GetName(), "foo.bar/treename");

   EXPECT_STREQ(files->At(4)->GetTitle(), "/path/to/file.root/foo.root");
   EXPECT_STREQ(files->At(4)->GetName(), "treename");

   EXPECT_STREQ(files->At(5)->GetTitle(), "/path/to/file.root");
   EXPECT_STREQ(files->At(5)->GetName(), "root/treename");

   EXPECT_STREQ(files->At(6)->GetTitle(), "path/to/file.root");
   EXPECT_STREQ(files->At(6)->GetName(), "treename");

   EXPECT_STREQ(files->At(7)->GetTitle(), "/path/to/file.root");
   EXPECT_STREQ(files->At(7)->GetName(), "/treename");
}
