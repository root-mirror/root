// @(#)root/io:$Id$
// Author: Jakob Blomer

/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_RRawFileDavix
#define ROOT_RRawFileDavix

#include <ROOT/RRawFile.hxx>
#include <ROOT/RStringView.hxx>

#include <cstddef>
#include <cstdint>
#include <memory>

namespace ROOT {
namespace Internal {

struct RDavixFileDes;

/**
 * \class RRawFileDavix RRawFileDavix.hxx
 *
 * The RRawFileDavix class provides read-only access to remote non-ROOT files.  It uses the Davix library for
 * the transport layer.  It instructs the RRawFile base class to buffer in larger chunks than the default for
 * local files, assuming that remote file access has high(er) latency.
 */
class RRawFileDavix : public RRawFile {
private:
   std::unique_ptr<Internal::RDavixFileDes> fFileDes;

protected:
   void OpenImpl() final;
   size_t ReadAtImpl(void *buffer, size_t nbytes, std::uint64_t offset) final;
   void ReadVImpl(RIOVec *ioVec, unsigned int nReq) final;
   std::uint64_t GetSizeImpl() final;

public:
   RRawFileDavix(std::string_view url, RRawFile::ROptions options);
   ~RRawFileDavix();
   std::unique_ptr<RRawFile> Clone() const final;
   int GetFeatures() const final { return kFeatureHasSize; }
};

} // namespace Internal
} // namespace ROOT

#endif
