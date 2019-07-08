/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <string>

#ifndef WIN32
#include <sys/types.h>
#include <sys/mman.h>
#endif

#include "TFile.h"
#include "TError.h"
#include "TSystem.h"

#include "TFileBufferRead.h"


// Download files in 128MB chunks.
#define CHUNK_SIZE (128*1024*1024)


TFileBufferRead::TFileBufferRead (TFile *file) :
   fFile(file)
{
   fSize = fFile->GetSize();
   fTotal = (fSize + CHUNK_SIZE - 1) / CHUNK_SIZE;
   fPresent.resize(fTotal, 0);

   if (TmpFile()) {
      fInvalid = false;
   }
#ifdef WIN32
   fInvalid = true;
#endif
}


TFileBufferRead::~TFileBufferRead() {
   if (fFd >= 0) {close(fFd);}
}

///
// Create a temporary file on disk and unlink it; use it as a buffer
// for the source TFile.
//
// - tmpdir: Directory where the temporary file should be created.
//
// Returns false on failure.
//
Bool_t TFileBufferRead::TmpFile(const std::string &tmpdir) {
   TString base = "root-shadow-";
   FILE *fp = gSystem->TempFileName(base, tmpdir.empty() ? nullptr : tmpdir.c_str());
   if (!fp) {
      Warning("TFileBufferRead", "Cannot create temporary file %s: %s (errno=%d)",
                                 base.Data(), strerror(errno), errno);
      return false;
   }
   int fd = fileno(fp);
   if (-1 == unlink(base.Data())) {
      Warning("TFileBufferRead", "Cannot unlink temporary file %s: %s (errno=%d)",
                                 base.Data(), strerror(errno), errno);
      return false;
   }
   if (-1 == ftruncate(fd, fSize)) {
      Warning("TFileBufferRead", "Cannot resize temporary file %s: %s (errno=%d)",
                                 base.Data(), strerror(errno), errno);
      // We ignore this error - it shouldn't be fatal.
   }
   fFd = fd;
   return true;
}


///
// Positional read of the file, using the underlying buffer file.
//
// If the buffer file is invalid, then this reads from the underlying file.
//
// On error, -1 is returned and errno is set appropriately.
Long_t TFileBufferRead::Pread(char *into, UInt_t n, Long64_t pos) {
   if (fInvalid) {
      // Note the order of arguments is different between POSIX read and
      // TFile's ReadBuffer.
      TFileCacheRead *origCache = fFile->GetCacheRead();
      fFile->SetCacheRead(nullptr, nullptr, TFile::kDoNotDisconnect);
      bool failed = fFile->ReadBuffer(into, pos, n);
      fFile->SetCacheRead(origCache, nullptr, TFile::kDoNotDisconnect);
      if (failed) {
        if (!errno) {errno = EIO;}
         return -1;
      }
   }

   if (!Cache(pos, pos + n)) {
      if (!errno) {errno = EIO;}
      return -1;
   }
   ssize_t retval;
   if (-1 == (retval = ::pread(fFd, into, n, pos))) {
      Warning("TFileBufferRead", "Failed to read from local buffer file: %s (errno=%d)",
                                 strerror(errno), errno);
      return -1;
   }
   return retval;
}


///
// Make sure the specified byte range is cached in the local buffer.
//
// Returns false on failure.
//
Bool_t TFileBufferRead::Cache(Long64_t start, Long64_t end) {
#ifdef WIN32
   return false;
#else
   if (fInvalid) return false;
   start = (start / CHUNK_SIZE) * CHUNK_SIZE;
   end = std::min(end, fSize);

   size_t index = start / CHUNK_SIZE;

   while (start < end) {
      size_t len = std::min(static_cast<off_t>(fSize - start), static_cast<off_t>(CHUNK_SIZE));
      if (!fPresent[index]) {

         void *window = mmap(0, len, PROT_READ | PROT_WRITE, MAP_SHARED, fFd, start);
         if (window == MAP_FAILED) {
            fInvalid = true;
            Warning("TFileBufferRead", "Unable to map a window of local buffer file: %s (errno=%d)",
                                       strerror(errno), errno);
            fInvalid = true;
            return false;
         }

         TFileCacheRead *origCache = fFile->GetCacheRead();
         fFile->SetCacheRead(nullptr, nullptr, TFile::kDoNotDisconnect);
         bool failed = fFile->ReadBuffer(static_cast<char *>(window), start, len);
         fFile->SetCacheRead(origCache, nullptr, TFile::kDoNotDisconnect);

         munmap(window, len);

         if (failed) {
            Warning("TFileBufferRead", "Failed to read into the buffer file.");
            fInvalid = true;
            return false;
         }

         fPresent[index] = 1;
         ++fCount;
         // TODO: If we provided a more complete wrapper around the source file,
         // we could close it once the file has been fully downloaded to local disk.
      }

      start += len;
      ++index;
   }
   return true;
#endif
}
