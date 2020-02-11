// @(#)root/io:$Id$
// Author: Jakob Blomer

/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "ROOT/RRawFileUnix.hxx"
#include "ROOT/RMakeUnique.hxx"

#include "TError.h"

#include <cerrno>
#include <cstring>
#include <stdexcept>
#include <string>
#include <utility>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

namespace {
constexpr int kDefaultBlockSize = 4096; // If fstat() does not provide a block size hint, use this value instead
} // anonymous namespace

ROOT::Internal::RRawFileUnix::RRawFileUnix(std::string_view url, ROptions options)
   : RRawFile(url, options), fFileDes(-1)
{
}

ROOT::Internal::RRawFileUnix::~RRawFileUnix()
{
   if (fFileDes >= 0)
      close(fFileDes);
}

std::unique_ptr<ROOT::Internal::RRawFile> ROOT::Internal::RRawFileUnix::Clone() const
{
   return std::make_unique<RRawFileUnix>(fUrl, fOptions);
}

std::uint64_t ROOT::Internal::RRawFileUnix::GetSizeImpl()
{
   struct stat info;
   int res = fstat(fFileDes, &info);
   if (res != 0)
      throw std::runtime_error("Cannot call fstat on '" + fUrl + "', error: " + std::string(strerror(errno)));
   return info.st_size;
}

void *ROOT::Internal::RRawFileUnix::MapImpl(size_t nbytes, std::uint64_t offset, std::uint64_t &mapdOffset)
{
   static std::uint64_t szPageBitmap = sysconf(_SC_PAGESIZE) - 1;
   mapdOffset = offset & ~szPageBitmap;
   nbytes += offset & szPageBitmap;

   void *result = mmap(nullptr, nbytes, PROT_READ, MAP_PRIVATE, fFileDes, mapdOffset);
   if (result == MAP_FAILED)
      throw std::runtime_error(std::string("Cannot perform memory mapping: ") + strerror(errno));
   return result;
}

void ROOT::Internal::RRawFileUnix::OpenImpl()
{
   fFileDes = open(GetLocation(fUrl).c_str(), O_RDONLY);
   if (fFileDes < 0) {
      throw std::runtime_error("Cannot open '" + fUrl + "', error: " + std::string(strerror(errno)));
   }

   if (fOptions.fBlockSize >= 0)
      return;

   struct stat info;
   int res = fstat(fFileDes, &info);
   if (res != 0) {
      throw std::runtime_error("Cannot call fstat on '" + fUrl + "', error: " + std::string(strerror(errno)));
   }
   if (info.st_blksize > 0) {
      fOptions.fBlockSize = info.st_blksize;
   } else {
      fOptions.fBlockSize = kDefaultBlockSize;
   }
}

size_t ROOT::Internal::RRawFileUnix::ReadAtImpl(void *buffer, size_t nbytes, std::uint64_t offset)
{
   size_t total_bytes = 0;
   while (nbytes) {
      ssize_t res = pread(fFileDes, buffer, nbytes, offset);
      if (res < 0) {
         if (errno == EINTR)
            continue;
         throw std::runtime_error("Cannot read from '" + fUrl + "', error: " + std::string(strerror(errno)));
      } else if (res == 0) {
         return total_bytes;
      }
      R__ASSERT(static_cast<size_t>(res) <= nbytes);
      buffer = reinterpret_cast<unsigned char *>(buffer) + res;
      nbytes -= res;
      total_bytes += res;
      offset += res;
   }
   return total_bytes;
}

void ROOT::Internal::RRawFileUnix::UnmapImpl(void *region, size_t nbytes)
{
   int rv = munmap(region, nbytes);
   if (rv != 0)
      throw std::runtime_error(std::string("Cannot remove memory mapping: ") + strerror(errno));
}
