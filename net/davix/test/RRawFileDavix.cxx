#include "RConfigure.h"
#include "ROOT/RRawFileDavix.hxx"

#include <string>

#include "gtest/gtest.h"

using RRawFile = ROOT::Internal::RRawFile;
using RRawFileDavix = ROOT::Internal::RRawFileDavix;

TEST(RRawFileDavix, Idle)
{
   // Test construction and destruction if the URL is never opened
   RRawFile::ROptions options;
   RRawFileDavix("http://root.cern.ch/files/davix.test", options);
}

TEST(RRawFileDavix, Basics)
{
   std::string line;
   RRawFile::ROptions options;
   std::unique_ptr<RRawFileDavix> f(new RRawFileDavix("http://root.cern.ch/files/davix.test", options));
   f->Readln(line);
   EXPECT_STREQ("Hello, World", line.c_str());
   EXPECT_EQ(13u, f->GetSize());

   std::unique_ptr<RRawFileDavix> f2(new RRawFileDavix("http://root.cern.ch/files/davix.test.404", options));
   EXPECT_THROW(f2->Readln(line), std::runtime_error);
}


TEST(RRawFileDavix, Eof)
{
   char tail[4];
   tail[3] = '\0';
   RRawFile::ROptions options;
   std::unique_ptr<RRawFileDavix> f(new RRawFileDavix("http://root.cern.ch/files/davix.test", options));
   auto nbytes = f->ReadAt(tail, 10, f->GetSize() - 3);
   EXPECT_EQ(3u, nbytes);
   EXPECT_STREQ("ld\n", tail);
}


TEST(RRawFileDavix, ReadV)
{
   RRawFile::ROptions options;
   options.fBlockSize = 0;
   std::unique_ptr<RRawFileDavix> f(new RRawFileDavix("http://root.cern.ch/files/davix.test", options));

   char buffer[2];
   buffer[0] = buffer[1] = 0;
   RRawFile::RIOVec iovec[2];
   iovec[0].fBuffer = &buffer[0];
   iovec[0].fOffset = 0;
   iovec[0].fSize = 1;
   iovec[1].fBuffer = &buffer[1];
   iovec[1].fOffset = 11;
   iovec[1].fSize = 1;
   f->ReadV(iovec, 2);

   EXPECT_EQ(1U, iovec[0].fOutBytes);
   EXPECT_EQ(1U, iovec[1].fOutBytes);
   EXPECT_EQ('H', buffer[0]);
   EXPECT_EQ('d', buffer[1]);
}
