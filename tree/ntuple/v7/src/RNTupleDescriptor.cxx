/// \file RNTupleDescriptor.cxx
/// \ingroup NTuple ROOT7
/// \author Jakob Blomer <jblomer@cern.ch>
/// \date 2018-10-04
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/RError.hxx>
#include <ROOT/RField.hxx>
#include <ROOT/RNTupleDescriptor.hxx>
#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RNTupleUtil.hxx>
#include <ROOT/RStringView.hxx>

#include <RZip.h>
#include <TError.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <utility>

namespace {

/// The machine-independent serialization of meta-data wraps the header and footer as well as sub structures in
/// frames.  The frame layout is
///
/// -----------------------------------------------------------
/// |  TYPE           | DESCRIPTION                           |
/// |----------------------------------------------------------
/// | std::uint16_t   | Version used to write the frame       |
/// | std::uint16_t   | Minimum version for reading the frame |
/// | std::uint32_t   | Length of the frame incl. preamble    |
/// -----------------------------------------------------------
///
/// In addition, the header and footer store a 4 byte CRC32 checksum of the frame immediately after the frame.
/// The footer also repeats the frame size just before the CRC32 checksum.  That means, one can read the last 8 bytes
/// to determine the footer length, and the first 8 bytes to determine the header length.
///
/// Within the frames, integers of different lengths are stored in a machine-independent representation. Strings and
/// vectors store the number of items followed by the items. Time stamps are stored in number of seconds since the
/// UNIX epoch.

std::uint32_t SerializeInt64(std::int64_t val, void *buffer)
{
   if (buffer != nullptr) {
      auto bytes = reinterpret_cast<unsigned char *>(buffer);
      bytes[0] = (val & 0x00000000000000FF);
      bytes[1] = (val & 0x000000000000FF00) >> 8;
      bytes[2] = (val & 0x0000000000FF0000) >> 16;
      bytes[3] = (val & 0x00000000FF000000) >> 24;
      bytes[4] = (val & 0x000000FF00000000) >> 32;
      bytes[5] = (val & 0x0000FF0000000000) >> 40;
      bytes[6] = (val & 0x00FF000000000000) >> 48;
      bytes[7] = (val & 0xFF00000000000000) >> 56;
   }
   return 8;
}

std::uint32_t SerializeUInt64(std::uint64_t val, void *buffer)
{
   return SerializeInt64(val, buffer);
}

std::uint32_t DeserializeInt64(const void *buffer, std::int64_t *val)
{
   auto bytes = reinterpret_cast<const unsigned char *>(buffer);
   *val = std::int64_t(bytes[0]) + (std::int64_t(bytes[1]) << 8) +
          (std::int64_t(bytes[2]) << 16) + (std::int64_t(bytes[3]) << 24) +
          (std::int64_t(bytes[4]) << 32) + (std::int64_t(bytes[5]) << 40) +
          (std::int64_t(bytes[6]) << 48) + (std::int64_t(bytes[7]) << 56);
   return 8;
}

std::uint32_t DeserializeUInt64(const void *buffer, std::uint64_t *val)
{
   return DeserializeInt64(buffer, reinterpret_cast<std::int64_t *>(val));
}

std::uint32_t SerializeInt32(std::int32_t val, void *buffer)
{
   if (buffer != nullptr) {
      auto bytes = reinterpret_cast<unsigned char *>(buffer);
      bytes[0] = (val & 0x000000FF);
      bytes[1] = (val & 0x0000FF00) >> 8;
      bytes[2] = (val & 0x00FF0000) >> 16;
      bytes[3] = (val & 0xFF000000) >> 24;
   }
   return 4;
}

std::uint32_t SerializeUInt32(std::uint32_t val, void *buffer)
{
   return SerializeInt32(val, buffer);
}

std::uint32_t DeserializeInt32(const void *buffer, std::int32_t *val)
{
   auto bytes = reinterpret_cast<const unsigned char *>(buffer);
   *val = std::int32_t(bytes[0]) + (std::int32_t(bytes[1]) << 8) +
          (std::int32_t(bytes[2]) << 16) + (std::int32_t(bytes[3]) << 24);
   return 4;
}

std::uint32_t DeserializeUInt32(const void *buffer, std::uint32_t *val)
{
   return DeserializeInt32(buffer, reinterpret_cast<std::int32_t *>(val));
}

std::uint32_t SerializeInt16(std::int16_t val, void *buffer)
{
   if (buffer != nullptr) {
      auto bytes = reinterpret_cast<unsigned char *>(buffer);
      bytes[0] = (val & 0x00FF);
      bytes[1] = (val & 0xFF00) >> 8;
   }
   return 2;
}

std::uint32_t SerializeUInt16(std::uint16_t val, void *buffer)
{
   return SerializeInt16(val, buffer);
}

std::uint32_t DeserializeInt16(const void *buffer, std::int16_t *val)
{
   auto bytes = reinterpret_cast<const unsigned char *>(buffer);
   *val = std::int16_t(bytes[0]) + (std::int16_t(bytes[1]) << 8);
   return 2;
}

std::uint32_t DeserializeUInt16(const void *buffer, std::uint16_t *val)
{
   return DeserializeInt16(buffer, reinterpret_cast<std::int16_t *>(val));
}

std::uint32_t SerializeClusterSize(ROOT::Experimental::ClusterSize_t val, void *buffer)
{
   return SerializeUInt32(val, buffer);
}

std::uint32_t DeserializeClusterSize(const void *buffer, ROOT::Experimental::ClusterSize_t *val)
{
   std::uint32_t size;
   auto nbytes = DeserializeUInt32(buffer, &size);
   *val = size;
   return nbytes;
}

std::uint32_t SerializeString(const std::string &val, void *buffer)
{
   if (buffer != nullptr) {
      auto pos = reinterpret_cast<unsigned char *>(buffer);
      pos += SerializeUInt32(val.length(), pos);
      memcpy(pos, val.data(), val.length());
   }
   return SerializeUInt32(val.length(), nullptr) + val.length();
}

std::uint32_t DeserializeString(const void *buffer, std::string *val)
{
   auto base = reinterpret_cast<const unsigned char *>(buffer);
   auto bytes = base;
   std::uint32_t length;
   bytes += DeserializeUInt32(buffer, &length);
   val->resize(length);
   memcpy(&(*val)[0], bytes, length);
   return bytes + length - base;
}

std::uint32_t SerializeLocator(const ROOT::Experimental::RClusterDescriptor::RLocator &val, void *buffer)
{
   // In order to keep the meta-data small, we don't wrap the locator in a frame
   if (buffer != nullptr) {
      auto pos = reinterpret_cast<unsigned char *>(buffer);
      pos += SerializeInt64(val.fPosition, pos);
      pos += SerializeUInt32(val.fBytesOnStorage, pos);
      pos += SerializeString(val.fUrl, pos);
   }
   return SerializeString(val.fUrl, nullptr) + 12;
}

std::uint32_t DeserializeLocator(const void *buffer, ROOT::Experimental::RClusterDescriptor::RLocator *val)
{
   auto bytes = reinterpret_cast<const unsigned char *>(buffer);
   bytes += DeserializeInt64(bytes, &val->fPosition);
   bytes += DeserializeUInt32(bytes, &val->fBytesOnStorage);
   bytes += DeserializeString(bytes, &val->fUrl);
   return SerializeString(val->fUrl, nullptr) + 12;
}

std::uint32_t SerializeFrame(std::uint16_t protocolVersionCurrent, std::uint16_t protocolVersionMin, void *buffer,
   void **ptrSize)
{
   if (buffer != nullptr) {
      auto pos = reinterpret_cast<unsigned char *>(buffer);
      pos += SerializeUInt16(protocolVersionCurrent, pos); // The protocol version used to write the structure
      pos += SerializeUInt16(protocolVersionMin, pos); // The minimum protocol version required to read the data
      *ptrSize = pos;
      pos += SerializeUInt32(0, pos); // placeholder for the size of the frame
   }
   return 8;
}

std::uint32_t DeserializeFrame(std::uint16_t protocolVersion, const void *buffer, std::uint32_t *size)
{
   auto bytes = reinterpret_cast<const unsigned char *>(buffer);
   std::uint16_t protocolVersionAtWrite;
   std::uint16_t protocolVersionMinRequired;
   bytes += DeserializeUInt16(bytes, &protocolVersionAtWrite);
   bytes += DeserializeUInt16(bytes, &protocolVersionMinRequired);
   R__ASSERT(protocolVersionAtWrite >= protocolVersionMinRequired);
   R__ASSERT(protocolVersion >= protocolVersionMinRequired);
   bytes += DeserializeUInt32(bytes, size);
   return 8;
}

std::uint32_t SerializeVersion(const ROOT::Experimental::RNTupleVersion &val, void *buffer)
{
   auto base = reinterpret_cast<unsigned char *>((buffer != nullptr) ? buffer : 0);
   auto pos = base;
   void** where = (buffer == nullptr) ? &buffer : reinterpret_cast<void**>(&pos);

   void *ptrSize = nullptr;
   pos += SerializeFrame(0, 0, *where, &ptrSize);

   pos += SerializeUInt32(val.GetVersionUse(), *where);
   pos += SerializeUInt32(val.GetVersionMin(), *where);
   pos += SerializeUInt64(val.GetFlags(), *where);

   auto size = pos - base;
   SerializeUInt32(size, ptrSize);
   return size;
}

std::uint32_t DeserializeVersion(const void *buffer, ROOT::Experimental::RNTupleVersion *version)
{
   auto bytes = reinterpret_cast<const unsigned char *>(buffer);
   std::uint32_t frameSize;
   bytes += DeserializeFrame(0, bytes, &frameSize);

   std::uint32_t versionUse;
   std::uint32_t versionMin;
   std::uint64_t flags;
   bytes += DeserializeUInt32(bytes, &versionUse);
   bytes += DeserializeUInt32(bytes, &versionMin);
   bytes += DeserializeUInt64(bytes, &flags);
   *version = ROOT::Experimental::RNTupleVersion(versionUse, versionMin, flags);

   return frameSize;
}

std::uint32_t SerializeUuid(const ROOT::Experimental::RNTupleUuid &val, void *buffer)
{
   auto base = reinterpret_cast<unsigned char *>((buffer != nullptr) ? buffer : 0);
   auto pos = base;
   void** where = (buffer == nullptr) ? &buffer : reinterpret_cast<void**>(&pos);

   void *ptrSize = nullptr;
   pos += SerializeFrame(0, 0, *where, &ptrSize);

   pos += SerializeString(val, *where);

   auto size = pos - base;
   SerializeUInt32(size, ptrSize);
   return size;
}

std::uint32_t DeserializeUuid(const void *buffer, ROOT::Experimental::RNTupleUuid *uuid)
{
   auto bytes = reinterpret_cast<const unsigned char *>(buffer);
   std::uint32_t frameSize;
   bytes += DeserializeFrame(0, bytes, &frameSize);

   bytes += DeserializeString(bytes, uuid);

   return frameSize;
}

std::uint32_t SerializeColumnModel(const ROOT::Experimental::RColumnModel &val, void *buffer)
{
   auto base = reinterpret_cast<unsigned char *>((buffer != nullptr) ? buffer : 0);
   auto pos = base;
   void** where = (buffer == nullptr) ? &buffer : reinterpret_cast<void**>(&pos);

   void *ptrSize = nullptr;
   pos += SerializeFrame(0, 0, *where, &ptrSize);

   pos += SerializeInt32(static_cast<int>(val.GetType()), *where);
   pos += SerializeInt32(static_cast<int>(val.GetIsSorted()), *where);

   auto size = pos - base;
   SerializeUInt32(size, ptrSize);
   return size;
}

std::uint32_t DeserializeColumnModel(const void *buffer, ROOT::Experimental::RColumnModel *columnModel)
{
   auto bytes = reinterpret_cast<const unsigned char *>(buffer);
   std::uint32_t frameSize;
   bytes += DeserializeFrame(0, bytes, &frameSize);

   std::int32_t type;
   std::int32_t isSorted;
   bytes += DeserializeInt32(bytes, &type);
   bytes += DeserializeInt32(bytes, &isSorted);
   *columnModel = ROOT::Experimental::RColumnModel(static_cast<ROOT::Experimental::EColumnType>(type), isSorted);

   return frameSize;
}

std::uint32_t SerializeTimeStamp(const std::chrono::system_clock::time_point &val, void *buffer)
{
   return SerializeInt64(std::chrono::system_clock::to_time_t(val), buffer);
}

std::uint32_t DeserializeTimeStamp(const void *buffer, std::chrono::system_clock::time_point *timeStamp)
{
   std::int64_t secSinceUnixEpoch;
   auto size = DeserializeInt64(buffer, &secSinceUnixEpoch);
   *timeStamp = std::chrono::system_clock::from_time_t(secSinceUnixEpoch);
   return size;
}

std::uint32_t SerializeColumnRange(const ROOT::Experimental::RClusterDescriptor::RColumnRange &val, void *buffer)
{
   // To keep the cluster footers small, we don't put a frame around individual column ranges.
   if (buffer != nullptr) {
      auto pos = reinterpret_cast<unsigned char *>(buffer);
      // The column id is stored in SerializeFooter() for the column range and the page range altogether
      pos += SerializeUInt64(val.fFirstElementIndex, pos);
      pos += SerializeClusterSize(val.fNElements, pos);
      pos += SerializeInt64(val.fCompressionSettings, pos);
   }
   return 20;
}

std::uint32_t DeserializeColumnRange(const void *buffer,
   ROOT::Experimental::RClusterDescriptor::RColumnRange *columnRange)
{
   auto bytes = reinterpret_cast<const unsigned char *>(buffer);
   // The column id is set elsewhere (see AddClustersFromFooter())
   bytes += DeserializeUInt64(bytes, &columnRange->fFirstElementIndex);
   bytes += DeserializeClusterSize(bytes, &columnRange->fNElements);
   bytes += DeserializeInt64(bytes, &columnRange->fCompressionSettings);
   return 20;
}

std::uint32_t SerializePageInfo(const ROOT::Experimental::RClusterDescriptor::RPageRange::RPageInfo &val, void *buffer)
{
   // To keep the cluster footers small, we don't put a frame around individual page infos.
   if (buffer != nullptr) {
      auto pos = reinterpret_cast<unsigned char *>(buffer);
      // The column id is stored in SerializeFooter() for the column range and the page range altogether
      pos += SerializeClusterSize(val.fNElements, pos);
      pos += SerializeLocator(val.fLocator, pos);
   }
   return 4 + SerializeLocator(val.fLocator, nullptr);
}

std::uint32_t DeserializePageInfo(const void *buffer,
   ROOT::Experimental::RClusterDescriptor::RPageRange::RPageInfo *pageInfo)
{
   auto base = reinterpret_cast<const unsigned char *>(buffer);
   auto bytes = base;
   // The column id is set elsewhere (see AddClustersFromFooter())
   bytes += DeserializeClusterSize(bytes, &pageInfo->fNElements);
   bytes += DeserializeLocator(bytes, &pageInfo->fLocator);
   return bytes - base;
}

std::uint32_t SerializeCrc32(const unsigned char *data, std::uint32_t length, void *buffer)
{
   auto checksum = R__crc32(0, nullptr, 0);
   if (buffer != nullptr) {
      checksum = R__crc32(checksum, data, length);
      SerializeUInt32(checksum, buffer);
   }
   return 4;
}

void VerifyCrc32(const unsigned char *data, std::uint32_t length)
{
   auto checksumReal = R__crc32(0, nullptr, 0);
   checksumReal = R__crc32(checksumReal, data, length);
   std::uint32_t checksumFound;
   DeserializeUInt32(data + length, &checksumFound);
   R__ASSERT(checksumFound == checksumReal);
}

std::uint32_t SerializeField(const ROOT::Experimental::RFieldDescriptor &val, void *buffer)
{
   auto base = reinterpret_cast<unsigned char *>((buffer != nullptr) ? buffer : 0);
   auto pos = base;
   void** where = (buffer == nullptr) ? &buffer : reinterpret_cast<void**>(&pos);

   void *ptrSize = nullptr;
   pos += SerializeFrame(ROOT::Experimental::RFieldDescriptor::kFrameVersionCurrent,
      ROOT::Experimental::RFieldDescriptor::kFrameVersionMin, *where, &ptrSize);

   pos += SerializeUInt64(val.GetId(), *where);
   pos += SerializeVersion(val.GetFieldVersion(), *where);
   pos += SerializeVersion(val.GetTypeVersion(), *where);
   pos += SerializeString(val.GetFieldName(), *where);
   pos += SerializeString(val.GetFieldDescription(), *where);
   pos += SerializeString(val.GetTypeName(), *where);
   pos += SerializeUInt64(val.GetNRepetitions(), *where);
   pos += SerializeUInt32(static_cast<int>(val.GetStructure()), *where);
   pos += SerializeUInt64(val.GetParentId(), *where);
   pos += SerializeUInt32(val.GetLinkIds().size(), *where);
   for (const auto& l : val.GetLinkIds())
      pos += SerializeUInt64(l, *where);

   auto size = pos - base;
   SerializeUInt32(size, ptrSize);
   return size;
}

std::uint32_t SerializeColumn(const ROOT::Experimental::RColumnDescriptor &val, void *buffer)
{
   auto base = reinterpret_cast<unsigned char *>((buffer != nullptr) ? buffer : 0);
   auto pos = base;
   void** where = (buffer == nullptr) ? &buffer : reinterpret_cast<void**>(&pos);

   void *ptrSize = nullptr;
   pos += SerializeFrame(ROOT::Experimental::RColumnDescriptor::kFrameVersionCurrent,
      ROOT::Experimental::RColumnDescriptor::kFrameVersionMin, *where, &ptrSize);

   pos += SerializeUInt64(val.GetId(), *where);
   pos += SerializeVersion(val.GetVersion(), *where);
   pos += SerializeColumnModel(val.GetModel(), *where);
   pos += SerializeUInt64(val.GetFieldId(), *where);
   pos += SerializeUInt32(val.GetIndex(), *where);

   auto size = pos - base;
   SerializeUInt32(size, ptrSize);
   return size;
}

std::uint32_t SerializeClusterSummary(const ROOT::Experimental::RClusterDescriptor &val, void *buffer)
{
   auto base = reinterpret_cast<unsigned char *>((buffer != nullptr) ? buffer : 0);
   auto pos = base;
   void** where = (buffer == nullptr) ? &buffer : reinterpret_cast<void**>(&pos);

   void *ptrSize = nullptr;
   pos += SerializeFrame(ROOT::Experimental::RClusterDescriptor::kFrameVersionCurrent,
      ROOT::Experimental::RClusterDescriptor::kFrameVersionMin, *where, &ptrSize);

   pos += SerializeUInt64(val.GetId(), *where);
   pos += SerializeVersion(val.GetVersion(), *where);
   pos += SerializeUInt64(val.GetFirstEntryIndex(), *where);
   pos += SerializeUInt64(val.GetNEntries(), *where);
   pos += SerializeLocator(val.GetLocator(), *where);

   auto size = pos - base;
   SerializeUInt32(size, ptrSize);
   return size;
}

} // anonymous namespace


////////////////////////////////////////////////////////////////////////////////


bool ROOT::Experimental::RFieldDescriptor::operator==(const RFieldDescriptor &other) const {
   return fFieldId == other.fFieldId &&
          fFieldVersion == other.fFieldVersion &&
          fTypeVersion == other.fTypeVersion &&
          fFieldName == other.fFieldName &&
          fFieldDescription == other.fFieldDescription &&
          fTypeName == other.fTypeName &&
          fNRepetitions == other.fNRepetitions &&
          fStructure == other.fStructure &&
          fParentId == other.fParentId &&
          fLinkIds == other.fLinkIds;
}

ROOT::Experimental::RFieldDescriptor
ROOT::Experimental::RFieldDescriptor::Clone() const {
   RFieldDescriptor clone;
   clone.fFieldId = fFieldId;
   clone.fFieldVersion = fFieldVersion;
   clone.fTypeVersion = fTypeVersion;
   clone.fFieldName = fFieldName;
   clone.fFieldDescription = fFieldDescription;
   clone.fTypeName = fTypeName;
   clone.fNRepetitions = fNRepetitions;
   clone.fStructure = fStructure;
   clone.fParentId = fParentId;
   clone.fLinkIds = fLinkIds;
   return clone;
}


////////////////////////////////////////////////////////////////////////////////


bool ROOT::Experimental::RColumnDescriptor::operator==(const RColumnDescriptor &other) const {
   return fColumnId == other.fColumnId &&
          fVersion == other.fVersion &&
          fModel == other.fModel &&
          fFieldId == other.fFieldId &&
          fIndex == other.fIndex;
}


////////////////////////////////////////////////////////////////////////////////


bool ROOT::Experimental::RClusterDescriptor::operator==(const RClusterDescriptor &other) const {
   return fClusterId == other.fClusterId &&
          fVersion == other.fVersion &&
          fFirstEntryIndex == other.fFirstEntryIndex &&
          fNEntries == other.fNEntries &&
          fLocator == other.fLocator &&
          fColumnRanges == other.fColumnRanges &&
          fPageRanges == other.fPageRanges;
}


////////////////////////////////////////////////////////////////////////////////


bool ROOT::Experimental::RNTupleDescriptor::operator==(const RNTupleDescriptor &other) const {
   return fName == other.fName &&
          fDescription == other.fDescription &&
          fAuthor == other.fAuthor &&
          fCustodian == other.fCustodian &&
          fTimeStampData == other.fTimeStampData &&
          fTimeStampWritten == other.fTimeStampWritten &&
          fVersion == other.fVersion &&
          fOwnUuid == other.fOwnUuid &&
          fGroupUuid == other.fGroupUuid &&
          fFieldDescriptors == other.fFieldDescriptors &&
          fColumnDescriptors == other.fColumnDescriptors &&
          fClusterDescriptors == other.fClusterDescriptors;
}


std::uint32_t ROOT::Experimental::RNTupleDescriptor::SerializeHeader(void* buffer) const
{
   auto base = reinterpret_cast<unsigned char *>((buffer != nullptr) ? buffer : 0);
   auto pos = base;
   void** where = (buffer == nullptr) ? &buffer : reinterpret_cast<void**>(&pos);

   void *ptrSize = nullptr;
   pos += SerializeFrame(
      RNTupleDescriptor::kFrameVersionCurrent, RNTupleDescriptor::kFrameVersionMin, *where, &ptrSize);
   pos += SerializeUInt64(0, *where); // reserved; can be at some point used, e.g., for compression flags

   pos += SerializeString(fName, *where);
   pos += SerializeString(fDescription, *where);
   pos += SerializeString(fAuthor, *where);
   pos += SerializeString(fCustodian, *where);
   pos += SerializeTimeStamp(fTimeStampData, *where);
   pos += SerializeTimeStamp(fTimeStampWritten, *where);
   pos += SerializeVersion(fVersion, *where);
   pos += SerializeUuid(fOwnUuid, *where);
   pos += SerializeUuid(fGroupUuid, *where);
   pos += SerializeUInt32(fFieldDescriptors.size(), *where);
   for (const auto& f : fFieldDescriptors) {
      pos += SerializeField(f.second, *where);
   }
   pos += SerializeUInt32(fColumnDescriptors.size(), *where);
   for (const auto& c : fColumnDescriptors) {
      pos += SerializeColumn(c.second, *where);
   }

   std::uint32_t size = pos - base;
   SerializeUInt32(size, ptrSize);
   size += SerializeCrc32(base, size, *where);

   return size;
}

std::uint32_t ROOT::Experimental::RNTupleDescriptor::SerializeFooter(void* buffer) const
{
   auto base = reinterpret_cast<unsigned char *>((buffer != nullptr) ? buffer : 0);
   auto pos = base;
   void** where = (buffer == nullptr) ? &buffer : reinterpret_cast<void**>(&pos);

   void *ptrSize = nullptr;
   pos += SerializeFrame(
      RNTupleDescriptor::kFrameVersionCurrent, RNTupleDescriptor::kFrameVersionMin, *where, &ptrSize);
   pos += SerializeUInt64(0, *where); // reserved; can be at some point used, e.g., for compression flags

   pos += SerializeUInt64(fClusterDescriptors.size(), *where);
   for (const auto& cluster : fClusterDescriptors) {
      pos += SerializeUuid(fOwnUuid, *where); // in order to verify that header and footer belong together
      pos += SerializeClusterSummary(cluster.second, *where);

      pos += SerializeUInt32(fColumnDescriptors.size(), *where);
      for (const auto& column : fColumnDescriptors) {
         auto columnId = column.first;
         pos += SerializeUInt64(columnId, *where);

         const auto &columnRange = cluster.second.GetColumnRange(columnId);
         R__ASSERT(columnRange.fColumnId == columnId);
         pos += SerializeColumnRange(columnRange, *where);

         const auto &pageRange = cluster.second.GetPageRange(columnId);
         R__ASSERT(pageRange.fColumnId == columnId);
         auto nPages = pageRange.fPageInfos.size();
         pos += SerializeUInt32(nPages, *where);
         for (unsigned int i = 0; i < nPages; ++i) {
            pos += SerializePageInfo(pageRange.fPageInfos[i], *where);
         }
      }
   }

   // The next 16 bytes make the ntuple's postscript
   pos += SerializeUInt16(kFrameVersionCurrent, *where);
   pos += SerializeUInt16(kFrameVersionMin, *where);
   // Add the CRC32 bytes to the header and footer sizes
   pos += SerializeUInt32(GetHeaderSize(), *where);
   std::uint32_t size = pos - base + 4;
   pos += SerializeUInt32(size + 4, *where);
   size += SerializeCrc32(base, size, *where);

   return size;
}


void ROOT::Experimental::RNTupleDescriptor::LocateMetadata(
   const void *postscript, std::uint32_t &szHeader, std::uint32_t &szFooter)
{
   auto pos = reinterpret_cast<const unsigned char *>(postscript);
   std::uint16_t dummy;
   pos += DeserializeUInt16(pos, &dummy);
   pos += DeserializeUInt16(pos, &dummy);
   pos += DeserializeUInt32(pos, &szHeader);
   pos += DeserializeUInt32(pos, &szFooter);
}


ROOT::Experimental::NTupleSize_t ROOT::Experimental::RNTupleDescriptor::GetNEntries() const
{
   NTupleSize_t result = 0;
   for (const auto &cd : fClusterDescriptors) {
      result = std::max(result, cd.second.GetFirstEntryIndex() + cd.second.GetNEntries());
   }
   return result;
}

ROOT::Experimental::NTupleSize_t ROOT::Experimental::RNTupleDescriptor::GetNElements(DescriptorId_t columnId) const
{
   NTupleSize_t result = 0;
   for (const auto &cd : fClusterDescriptors) {
      auto columnRange = cd.second.GetColumnRange(columnId);
      result = std::max(result, columnRange.fFirstElementIndex + columnRange.fNElements);
   }
   return result;
}


ROOT::Experimental::DescriptorId_t
ROOT::Experimental::RNTupleDescriptor::FindFieldId(std::string_view fieldName, DescriptorId_t parentId) const
{
   std::string leafName(fieldName);
   auto posDot = leafName.find_last_of('.');
   if (posDot != std::string::npos) {
      auto parentName = leafName.substr(0, posDot);
      leafName = leafName.substr(posDot + 1);
      parentId = FindFieldId(parentName, parentId);
   }
   for (const auto &fd : fFieldDescriptors) {
      if (fd.second.GetParentId() == parentId && fd.second.GetFieldName() == leafName)
         return fd.second.GetId();
   }
   return kInvalidDescriptorId;
}


std::string ROOT::Experimental::RNTupleDescriptor::GetQualifiedFieldName(DescriptorId_t fieldId) const
{
   if (fieldId == kInvalidDescriptorId)
      return "";

   const auto &fieldDescriptor = fFieldDescriptors.at(fieldId);
   auto prefix = GetQualifiedFieldName(fieldDescriptor.GetParentId());
   if (prefix.empty())
      return fieldDescriptor.GetFieldName();
   return prefix + "." + fieldDescriptor.GetFieldName();
}


ROOT::Experimental::DescriptorId_t
ROOT::Experimental::RNTupleDescriptor::GetFieldZeroId() const
{
   return FindFieldId("", kInvalidDescriptorId);
}


ROOT::Experimental::DescriptorId_t
ROOT::Experimental::RNTupleDescriptor::FindFieldId(std::string_view fieldName) const
{
   return FindFieldId(fieldName, GetFieldZeroId());
}


ROOT::Experimental::DescriptorId_t
ROOT::Experimental::RNTupleDescriptor::FindColumnId(DescriptorId_t fieldId, std::uint32_t columnIndex) const
{
   for (const auto &cd : fColumnDescriptors) {
      if (cd.second.GetFieldId() == fieldId && cd.second.GetIndex() == columnIndex)
        return cd.second.GetId();
   }
   return kInvalidDescriptorId;
}


ROOT::Experimental::DescriptorId_t
ROOT::Experimental::RNTupleDescriptor::FindClusterId(DescriptorId_t columnId, NTupleSize_t index) const
{
   // TODO(jblomer): binary search?
   for (const auto &cd : fClusterDescriptors) {
      auto columnRange = cd.second.GetColumnRange(columnId);
      if (columnRange.Contains(index))
         return cd.second.GetId();
   }
   return kInvalidDescriptorId;
}


ROOT::Experimental::DescriptorId_t
ROOT::Experimental::RNTupleDescriptor::FindNextClusterId(DescriptorId_t clusterId) const
{
   const auto &clusterDesc = GetClusterDescriptor(clusterId);
   auto firstEntryInNextCluster = clusterDesc.GetFirstEntryIndex() + clusterDesc.GetNEntries();
   // TODO(jblomer): binary search?
   for (const auto &cd : fClusterDescriptors) {
      if (cd.second.GetFirstEntryIndex() == firstEntryInNextCluster)
         return cd.second.GetId();
   }
   return kInvalidDescriptorId;
}


ROOT::Experimental::DescriptorId_t
ROOT::Experimental::RNTupleDescriptor::FindPrevClusterId(DescriptorId_t clusterId) const
{
   const auto &clusterDesc = GetClusterDescriptor(clusterId);
   // TODO(jblomer): binary search?
   for (const auto &cd : fClusterDescriptors) {
      if (cd.second.GetFirstEntryIndex() + cd.second.GetNEntries() == clusterDesc.GetFirstEntryIndex())
         return cd.second.GetId();
   }
   return kInvalidDescriptorId;
}


std::unique_ptr<ROOT::Experimental::RNTupleModel> ROOT::Experimental::RNTupleDescriptor::GenerateModel() const
{
   auto model = std::make_unique<RNTupleModel>();
   for (const auto &topDesc : GetTopLevelFields()) {
      auto field = Detail::RFieldBase::Create(topDesc.GetFieldName(), topDesc.GetTypeName());
      model->AddField(field.Unwrap());
   }
   return model;
}


////////////////////////////////////////////////////////////////////////////////


ROOT::Experimental::RResult<void>
ROOT::Experimental::RNTupleDescriptorBuilder::EnsureValidDescriptor() const {
   // Reuse field name validity check
   auto validName = Detail::RFieldBase::EnsureValidFieldName(fDescriptor.GetName());
   if (!validName) {
      return R__FORWARD_ERROR(validName);
   }
   // open-ended list of invariant checks
   for (const auto& key_val: fDescriptor.fFieldDescriptors) {
      const auto& id = key_val.first;
      const auto& desc = key_val.second;
      // parent not properly set
      if (id != DescriptorId_t(0) && desc.GetParentId() == kInvalidDescriptorId) {
         return R__FAIL("field with id '" + std::to_string(id) + "' has an invalid parent id");
      }
   }
   return RResult<void>::Success();
}

ROOT::Experimental::RNTupleDescriptor ROOT::Experimental::RNTupleDescriptorBuilder::MoveDescriptor()
{
   RNTupleDescriptor result;
   std::swap(result, fDescriptor);
   return result;
}

void ROOT::Experimental::RNTupleDescriptorBuilder::SetFromHeader(void* headerBuffer)
{
   auto pos = reinterpret_cast<unsigned char *>(headerBuffer);
   auto base = pos;

   std::uint32_t frameSize;
   pos += DeserializeFrame(RNTupleDescriptor::kFrameVersionCurrent, base, &frameSize);
   VerifyCrc32(base, frameSize);
   std::uint64_t reserved;
   pos += DeserializeUInt64(pos, &reserved);

   pos += DeserializeString(pos, &fDescriptor.fName);
   pos += DeserializeString(pos, &fDescriptor.fDescription);
   pos += DeserializeString(pos, &fDescriptor.fAuthor);
   pos += DeserializeString(pos, &fDescriptor.fCustodian);
   pos += DeserializeTimeStamp(pos, &fDescriptor.fTimeStampData);
   pos += DeserializeTimeStamp(pos, &fDescriptor.fTimeStampWritten);
   pos += DeserializeVersion(pos, &fDescriptor.fVersion);
   pos += DeserializeUuid(pos, &fDescriptor.fOwnUuid);
   pos += DeserializeUuid(pos, &fDescriptor.fGroupUuid);

   std::uint32_t nFields;
   pos += DeserializeUInt32(pos, &nFields);
   for (std::uint32_t i = 0; i < nFields; ++i) {
      auto fieldBase = pos;
      pos += DeserializeFrame(RFieldDescriptor::kFrameVersionCurrent, fieldBase, &frameSize);

      RFieldDescriptor f;
      pos += DeserializeUInt64(pos, &f.fFieldId);
      pos += DeserializeVersion(pos, &f.fFieldVersion);
      pos += DeserializeVersion(pos, &f.fTypeVersion);
      pos += DeserializeString(pos, &f.fFieldName);
      pos += DeserializeString(pos, &f.fFieldDescription);
      pos += DeserializeString(pos, &f.fTypeName);
      pos += DeserializeUInt64(pos, &f.fNRepetitions);
      std::int32_t structure;
      pos += DeserializeInt32(pos, &structure);
      f.fStructure = static_cast<ENTupleStructure>(structure);
      pos += DeserializeUInt64(pos, &f.fParentId);

      std::uint32_t nLinks;
      pos += DeserializeUInt32(pos, &nLinks);
      f.fLinkIds.resize(nLinks);
      for (std::uint32_t j = 0; j < nLinks; ++j) {
         pos += DeserializeUInt64(pos, &f.fLinkIds[j]);
      }

      pos = fieldBase + frameSize;
      fDescriptor.fFieldDescriptors.emplace(f.fFieldId, std::move(f));
   }

   std::uint32_t nColumns;
   pos += DeserializeUInt32(pos, &nColumns);
   for (std::uint32_t i = 0; i < nColumns; ++i) {
      auto columnBase = pos;
      pos += DeserializeFrame(RColumnDescriptor::kFrameVersionCurrent, columnBase, &frameSize);

      RColumnDescriptor c;
      pos += DeserializeUInt64(pos, &c.fColumnId);
      pos += DeserializeVersion(pos, &c.fVersion);
      pos += DeserializeColumnModel(pos, &c.fModel);
      pos += DeserializeUInt64(pos, &c.fFieldId);
      pos += DeserializeUInt32(pos, &c.fIndex);

      pos = columnBase + frameSize;
      fDescriptor.fColumnDescriptors.emplace(c.fColumnId, std::move(c));
   }
}

void ROOT::Experimental::RNTupleDescriptorBuilder::AddClustersFromFooter(void* footerBuffer) {
   auto pos = reinterpret_cast<unsigned char *>(footerBuffer);
   auto base = pos;

   std::uint32_t frameSize;
   pos += DeserializeFrame(RNTupleDescriptor::kFrameVersionCurrent, pos, &frameSize);
   VerifyCrc32(base, frameSize);
   std::uint64_t reserved;
   pos += DeserializeUInt64(pos, &reserved);

   std::uint64_t nClusters;
   pos += DeserializeUInt64(pos, &nClusters);
   for (std::uint64_t i = 0; i < nClusters; ++i) {
      RNTupleUuid uuid;
      pos += DeserializeUuid(pos, &uuid);
      R__ASSERT(uuid == fDescriptor.fOwnUuid);
      auto clusterBase = pos;
      pos += DeserializeFrame(RClusterDescriptor::kFrameVersionCurrent, clusterBase, &frameSize);

      std::uint64_t clusterId;
      RNTupleVersion version;
      std::uint64_t firstEntry;
      std::uint64_t nEntries;
      pos += DeserializeUInt64(pos, &clusterId);
      pos += DeserializeVersion(pos, &version);
      pos += DeserializeUInt64(pos, &firstEntry);
      pos += DeserializeUInt64(pos, &nEntries);
      AddCluster(clusterId, version, firstEntry, ROOT::Experimental::ClusterSize_t(nEntries));
      RClusterDescriptor::RLocator locator;
      pos += DeserializeLocator(pos, &locator);
      SetClusterLocator(clusterId, locator);

      pos = clusterBase + frameSize;

      std::uint32_t nColumns;
      pos += DeserializeUInt32(pos, &nColumns);
      for (std::uint32_t j = 0; j < nColumns; ++j) {
         uint64_t columnId;
         pos += DeserializeUInt64(pos, &columnId);

         RClusterDescriptor::RColumnRange columnRange;
         columnRange.fColumnId = columnId;
         pos += DeserializeColumnRange(pos, &columnRange);
         AddClusterColumnRange(clusterId, columnRange);

         RClusterDescriptor::RPageRange pageRange;
         pageRange.fColumnId = columnId;
         uint32_t nPages;
         pos += DeserializeUInt32(pos, &nPages);
         for (unsigned int k = 0; k < nPages; ++k) {
            RClusterDescriptor::RPageRange::RPageInfo pageInfo;
            pos += DeserializePageInfo(pos, &pageInfo);
            pageRange.fPageInfos.emplace_back(pageInfo);
         }
         AddClusterPageRange(clusterId, std::move(pageRange));
      }
   }
}

void ROOT::Experimental::RNTupleDescriptorBuilder::SetNTuple(
   const std::string_view name, const std::string_view description, const std::string_view author,
   const RNTupleVersion &version, const RNTupleUuid &uuid)
{
   fDescriptor.fName = std::string(name);
   fDescriptor.fDescription = std::string(description);
   fDescriptor.fAuthor = std::string(author);
   fDescriptor.fVersion = version;
   fDescriptor.fOwnUuid = uuid;
   fDescriptor.fGroupUuid = uuid;
}

ROOT::Experimental::RDanglingFieldDescriptor::RDanglingFieldDescriptor(
   const RFieldDescriptor& fieldDesc) : fField(fieldDesc.Clone())
{
   fField.fParentId = kInvalidDescriptorId;
   fField.fLinkIds = {};
}

ROOT::Experimental::RDanglingFieldDescriptor
ROOT::Experimental::RDanglingFieldDescriptor::FromField(const Detail::RFieldBase& field) {
   RDanglingFieldDescriptor fieldDesc;
   fieldDesc.FieldVersion(field.GetFieldVersion())
      .TypeVersion(field.GetTypeVersion())
      .FieldName(field.GetName())
      .TypeName(field.GetType())
      .Structure(field.GetStructure())
      .NRepetitions(field.GetNRepetitions());
   return fieldDesc;
}

ROOT::Experimental::RResult<ROOT::Experimental::RFieldDescriptor>
ROOT::Experimental::RDanglingFieldDescriptor::MakeDescriptor() const {
   if (fField.GetId() == kInvalidDescriptorId) {
      return R__FAIL("invalid field id");
   }
   if (fField.GetStructure() == ENTupleStructure::kInvalid) {
      return R__FAIL("invalid field structure");
   }
   // FieldZero is usually named "" and would be a false positive here
   if (fField.GetId() != DescriptorId_t(0)) {
      auto validName = Detail::RFieldBase::EnsureValidFieldName(fField.GetFieldName());
      if (!validName) {
         return R__FORWARD_ERROR(validName);
      }
   }
   return fField.Clone();
}

void ROOT::Experimental::RNTupleDescriptorBuilder::AddField(const RFieldDescriptor& fieldDesc) {
   fDescriptor.fFieldDescriptors.emplace(fieldDesc.GetId(), fieldDesc.Clone());
}

ROOT::Experimental::RResult<void>
ROOT::Experimental::RNTupleDescriptorBuilder::AddFieldLink(DescriptorId_t fieldId, DescriptorId_t linkId)
{
   if (linkId == DescriptorId_t(0)) {
      return R__FAIL("cannot make FieldZero a child field");
   }
   if (fDescriptor.fFieldDescriptors.count(fieldId) == 0) {
      return R__FAIL("field with id '" + std::to_string(fieldId) + "' doesn't exist in NTuple");
   }
   if (fDescriptor.fFieldDescriptors.count(linkId) == 0) {
      return R__FAIL("child field with id '" + std::to_string(linkId) + "' doesn't exist in NTuple");
   }
   // fail if field already has a valid parent
   auto parentId = fDescriptor.fFieldDescriptors.at(linkId).GetParentId();
   if (parentId != kInvalidDescriptorId) {
      return R__FAIL("field '" + std::to_string(linkId) + "' already has a parent ('" +
         std::to_string(parentId) + ")");
   }
   if (fieldId == linkId) {
      return R__FAIL("cannot make field '" + std::to_string(fieldId) + "' a child of itself");
   }
   fDescriptor.fFieldDescriptors.at(linkId).fParentId = fieldId;
   fDescriptor.fFieldDescriptors.at(fieldId).fLinkIds.push_back(linkId);
   return RResult<void>::Success();
}

void ROOT::Experimental::RNTupleDescriptorBuilder::AddColumn(
   DescriptorId_t columnId, DescriptorId_t fieldId, const RNTupleVersion &version, const RColumnModel &model,
   std::uint32_t index)
{
   RColumnDescriptor c;
   c.fColumnId = columnId;
   c.fFieldId = fieldId;
   c.fVersion = version;
   c.fModel = model;
   c.fIndex = index;
   fDescriptor.fColumnDescriptors.emplace(columnId, std::move(c));
}

void ROOT::Experimental::RNTupleDescriptorBuilder::AddCluster(
   DescriptorId_t clusterId, RNTupleVersion version, NTupleSize_t firstEntryIndex, ClusterSize_t nEntries)
{
   RClusterDescriptor c;
   c.fClusterId = clusterId;
   c.fVersion = version;
   c.fFirstEntryIndex = firstEntryIndex;
   c.fNEntries = nEntries;
   fDescriptor.fClusterDescriptors.emplace(clusterId, std::move(c));
}

void ROOT::Experimental::RNTupleDescriptorBuilder::SetClusterLocator(DescriptorId_t clusterId,
   RClusterDescriptor::RLocator locator)
{
   fDescriptor.fClusterDescriptors[clusterId].fLocator = locator;
}

void ROOT::Experimental::RNTupleDescriptorBuilder::AddClusterColumnRange(
   DescriptorId_t clusterId, const RClusterDescriptor::RColumnRange &columnRange)
{
   fDescriptor.fClusterDescriptors[clusterId].fColumnRanges[columnRange.fColumnId] = columnRange;
}

void ROOT::Experimental::RNTupleDescriptorBuilder::AddClusterPageRange(
   DescriptorId_t clusterId, RClusterDescriptor::RPageRange &&pageRange)
{
   fDescriptor.fClusterDescriptors[clusterId].fPageRanges.emplace(pageRange.fColumnId, std::move(pageRange));
}
