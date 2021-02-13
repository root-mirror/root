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
#include <ROOT/RStringView.hxx>

#include <cstddef>
#include <cstdint>
#include <cstdio>

namespace ROOT {
namespace Internal {

/**
 * \class RRawFileWin RRawFileWin.hxx
 * \ingroup IO
 *
 * The RRawFileWin class uses portable C I/O calls to read from a drive. The standard C I/O buffering is turned off
 * for the buffering of RRawFile base class.
 */
class RRawFileWin : public RRawFile {
private:
   FILE *fFilePtr;
   void Seek(long offset, int whence);

protected:
   void OpenImpl() final;
   size_t ReadAtImpl(void *buffer, size_t nbytes, std::uint64_t offset) final;
   std::uint64_t GetSizeImpl() final;

public:
   RRawFileWin(std::string_view url, RRawFile::ROptions options);
   ~RRawFileWin();
   std::unique_ptr<RRawFile> Clone() const final;
   int GetFeatures() const final { return kFeatureHasSize; }
};

} // namespace Internal
} // namespace ROOT

#endif
