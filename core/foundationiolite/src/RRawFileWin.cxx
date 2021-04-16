// @(#)root/io:$Id$
// Author: Jakob Blomer

/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "ROOT/RRawFileWin.hxx"
#include "ROOT/RMakeUnique.hxx"

#include "TError.h"

#include <cerrno>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>

namespace {
constexpr int kDefaultBlockSize = 4096; // Read files in 4k pages unless told otherwise
} // anonymous namespace

ROOT::Internal::RRawFileWin::RRawFileWin(std::string_view url, ROptions options)
   : RRawFile(url, options), fFilePtr(nullptr)
{
}

ROOT::Internal::RRawFileWin::~RRawFileWin()
{
   if (fFilePtr != nullptr)
      fclose(fFilePtr);
}

std::unique_ptr<ROOT::Internal::RRawFile> ROOT::Internal::RRawFileWin::Clone() const
{
   return std::make_unique<RRawFileWin>(fUrl, fOptions);
}

std::uint64_t ROOT::Internal::RRawFileWin::GetSizeImpl()
{
   Seek(0L, SEEK_END);
   long size = ftell(fFilePtr);
   if (size < 0)
      throw std::runtime_error("Cannot tell position counter in '" + fUrl +
                               "', error: " + std::string(strerror(errno)));
   Seek(fFilePos, SEEK_SET);
   return size;
}

void ROOT::Internal::RRawFileWin::OpenImpl()
{
   fFilePtr = fopen(GetLocation(fUrl).c_str(), "rb");
   if (fFilePtr == nullptr)
      throw std::runtime_error("Cannot open '" + fUrl + "', error: " + std::string(strerror(errno)));
   // Prevent double buffering
   int res = setvbuf(fFilePtr, nullptr, _IONBF, 0);
   R__ASSERT(res == 0);
   if (fOptions.fBlockSize < 0)
      fOptions.fBlockSize = kDefaultBlockSize;
}

size_t ROOT::Internal::RRawFileWin::ReadAtImpl(void *buffer, size_t nbytes, std::uint64_t offset)
{
   Seek(offset, SEEK_SET);
   size_t res = fread(buffer, 1, nbytes, fFilePtr);
   if ((res < nbytes) && (ferror(fFilePtr) != 0)) {
      clearerr(fFilePtr);
      throw std::runtime_error("Cannot read from '" + fUrl + "', error: " + std::string(strerror(errno)));
   }
   return res;
}

void ROOT::Internal::RRawFileWin::Seek(long offset, int whence)
{
   int res = fseek(fFilePtr, offset, whence);
   if (res != 0)
      throw std::runtime_error("Cannot seek in '" + fUrl + "', error: " + std::string(strerror(errno)));
}
