// @(#)root/io:$Id$
// Author: Jakob Blomer

/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_RRawFileWin
#define ROOT_RRawFileWin

#include <ROOT/RRawFile.hxx>

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <string>

namespace ROOT {
namespace Detail {

// TODO(jblomer): use Win32 API instead of C I/O to avoid double buffering
class RRawFileWin : public RRawFile {
private:
   FILE *fileptr;
   void EnsureOpen();
   void Seek(long offset, int whence);

protected:
   size_t DoPread(void *buffer, size_t nbytes, std::uint64_t offset) final;
   std::uint64_t DoGetSize() final;

public:
   RRawFileWin(const std::string &url, RRawFile::ROptions options);
   ~RRawFileWin();
};

} // namespace Detail
} // namespace ROOT

#endif
