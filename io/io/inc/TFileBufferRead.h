/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TFileBufferRead
#define ROOT_TFileBufferRead

#include <string>
#include <vector>

#include "RtypesCore.h"

/**
 * Add a buffering layer in front of a TFile.
 *
 * This divides the file into 128MB chunks; the first time a byte from a chunk
 * is accessed, the whole chunk is downloaded to a temporary file on disk.  A
 * file chunk is only downloaded if accessed.
 *
 * This is *not* meant to be a cache: as soon as the TFile is closed, the buffer
 * is removed.  This means repeated access to the same file will still cause it
 * to be downloaded multiple times.
 *
 * Rather, the focus of this buffer is to hide the effects of network latency.
 *
 * Implementation is based on a similar approach from CMSSW:
 * https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/Utilities/StorageFactory/src/LocalCacheFile.cc
 */
class TFile;

class TFileBufferRead {
public:
   TFileBufferRead(TFile *file);
   ~TFileBufferRead();

   Long_t Pread(char *into, UInt_t n, Long64_t pos);

   ///
   // Return the number of file chunks downloaded to the local buffer.
   UInt_t GetCount() const {return fCount;}

private:
   Bool_t TmpFile(const std::string &tmpdir="");
   Bool_t Cache(Long64_t start, Long64_t end);

   std::vector<char> fPresent;       // A mask of all the currently present file chunks in the buffer.
   TFile            *fFile{nullptr}; // A copy of the TFile we are buffering.
   Long64_t          fSize{-1};      // Size of the source TFile.
   UInt_t            fCount{0};      // Number of file chunks we have buffered.
   Int_t             fTotal{-1};     // Total number of chunks in source TFile.
   Int_t             fFd{-1};        // File descriptor pointing to the local file.
   Bool_t            fInvalid{true}; // Set to true if this buffer is in an invalid state
};

#endif  // ROOT_TFileBufferRead
