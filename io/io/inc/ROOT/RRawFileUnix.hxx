// @(#)root/io:$Id$
// Author: Jakob Blomer

/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_RRawFileUnix
#define ROOT_RRawFileUnix

#include <ROOT/RRawFile.hxx>
#include <ROOT/RStringView.hxx>

#include <cstddef>
#include <cstdint>

namespace ROOT {
namespace Experimental {
namespace Detail {

/**
 * \class RRawFileUnix RRawFileUnix.hxx
 * \ingroup IO
 *
 * The RRawFileUnix class uses POSIX calls to read from a mounted file system. Thus the path name can refer,
 * for instance, to a named pipe instead of a regular file.
 */
class RRawFileUnix : public RRawFile {
private:
   int fFileDes;

protected:
   void DoOpen() final;
   size_t DoReadAt(void *buffer, size_t nbytes, std::uint64_t offset) final;
   std::uint64_t DoGetSize() final;
   void *DoMap(size_t nbytes, std::uint64_t offset, std::uint64_t &mapdOffset) final;
   void DoUnmap(void *region, size_t nbytes) final;

public:
   RRawFileUnix(std::string_view url, RRawFile::ROptions options);
   ~RRawFileUnix();
   std::unique_ptr<RRawFile> Clone() const final;
   int GetFeatures() const final { return kFeatureHasSize | kFeatureHasMmap; }
};

} // namespace Detail
} // namespace Experimental
} // namespace ROOT

#endif
