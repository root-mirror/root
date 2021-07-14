/// \file RMiniFile.cxx
/// \ingroup NTuple ROOT7
/// \author Jakob Blomer <jblomer@cern.ch>
/// \date 2019-12-22
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/RConfig.hxx>

#include "ROOT/RMiniFile.hxx"

#include <ROOT/RRawFile.hxx>
#include <ROOT/RNTupleZip.hxx>

#include <TError.h>
#include <TFile.h>
#include <TKey.h>

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <string>
#include <utility>
#include <chrono>

namespace {

// The following types are used to read and write the TFile binary format

/// Big-endian 16-bit unsigned integer
class RUInt16BE {
private:
   std::uint16_t fValBE = 0;
   static std::uint16_t Swap(std::uint16_t val) {
      return (val & 0x00FF) << 8 | (val & 0xFF00) >> 8;
   }
public:
   RUInt16BE() = default;
   explicit RUInt16BE(const std::uint16_t val) : fValBE(Swap(val)) {}
   operator std::uint16_t() const {
      return Swap(fValBE);
   }
   RUInt16BE& operator =(const std::uint16_t val) {
      fValBE = Swap(val);
      return *this;
   }
};

/// Big-endian 32-bit unsigned integer
class RUInt32BE {
private:
   std::uint32_t fValBE = 0;
   static std::uint32_t Swap(std::uint32_t val) {
      auto x = (val & 0x0000FFFF) << 16 | (val & 0xFFFF0000) >> 16;
      return (x & 0x00FF00FF) << 8 | (x & 0xFF00FF00) >> 8;
   }
public:
   RUInt32BE() = default;
   explicit RUInt32BE(const std::uint32_t val) : fValBE(Swap(val)) {}
   operator std::uint32_t() const {
      return Swap(fValBE);
   }
   RUInt32BE& operator =(const std::uint32_t val) {
      fValBE = Swap(val);
      return *this;
   }
};

/// Big-endian 32-bit signed integer
class RInt32BE {
private:
   std::int32_t fValBE = 0;
   static std::int32_t Swap(std::int32_t val) {
      auto x = (val & 0x0000FFFF) << 16 | (val & 0xFFFF0000) >> 16;
      return (x & 0x00FF00FF) << 8 | (x & 0xFF00FF00) >> 8;
   }
public:
   RInt32BE() = default;
   explicit RInt32BE(const std::int32_t val) : fValBE(Swap(val)) {}
   operator std::int32_t() const {
      return Swap(fValBE);
   }
   RInt32BE& operator =(const std::int32_t val) {
      fValBE = Swap(val);
      return *this;
   }
};

/// Big-endian 64-bit unsigned integer
class RUInt64BE {
private:
   std::uint64_t fValBE = 0;
   static std::uint64_t Swap(std::uint64_t val) {
      auto x = (val & 0x00000000FFFFFFFF) << 32 | (val & 0xFFFFFFFF00000000) >> 32;
      x = (x & 0x0000FFFF0000FFFF) << 16 | (x & 0xFFFF0000FFFF0000) >> 16;
      return (x & 0x00FF00FF00FF00FF) << 8  | (x & 0xFF00FF00FF00FF00) >> 8;
   }
public:
   RUInt64BE() = default;
   explicit RUInt64BE(const std::uint64_t val) : fValBE(Swap(val)) {}
   operator std::uint64_t() const {
      return Swap(fValBE);
   }
   RUInt64BE& operator =(const std::uint64_t val) {
      fValBE = Swap(val);
      return *this;
   }
};

/// Composition of class RNTuple as being interpreted by TFile
constexpr std::int32_t ChecksumRNTupleClass() {
   const char ident[] = "ROOT::Experimental::RNTuple"
      "fVersion"
      "unsigned int"
      "fSize"
      "unsigned int"
      "fSeekHeader"
      "unsigned long"
      "fNBytesHeader"
      "unsigned int"
      "fLenHeader"
      "unsigned int"
      "fSeekFooter"
      "unsigned long"
      "fNBytesFooter"
      "unsigned int"
      "fLenFooter"
      "unsigned int"
      "fReserved"
      "unsigned long";
   std::int32_t id = 0;
   for (unsigned i = 0; i < (sizeof(ident) - 1); i++)
      id = static_cast<std::int32_t>(static_cast<std::int64_t>(id) * 3 + ident[i]);
   return id;
}


#pragma pack(push, 1)
/// A name (type, identifies, ...) in the TFile binary format
struct RTFString {
   char fLName{0};
   char fData[255];
   RTFString() = default;
   RTFString(const std::string &str) {
      R__ASSERT(str.length() < 256);
      fLName = str.length();
      memcpy(fData, str.data(), fLName);
   }
   char GetSize() const { return 1 + fLName; }
};

/// The timestamp format used in TFile; the default constructor initializes with the current time
struct RTFDatetime {
   RUInt32BE fDatetime;
   RTFDatetime() {
      auto now = std::chrono::system_clock::now();
      auto tt = std::chrono::system_clock::to_time_t(now);
      auto tm = *localtime(&tt);
      fDatetime = (tm.tm_year + 1900 - 1995) << 26 | (tm.tm_mon + 1) << 22 | tm.tm_mday << 17 |
                  tm.tm_hour << 12 | tm.tm_min << 6 | tm.tm_sec;
   }
   explicit RTFDatetime(RUInt32BE val) : fDatetime(val) {}
};

/// The key part of a TFile record excluding the class, object, and title names
struct RTFKey {
   RInt32BE fNbytes{0};
   RUInt16BE fVersion{4};
   RUInt32BE fObjLen{0};
   RTFDatetime fDatetime;
   RUInt16BE fKeyLen{0};
   RUInt16BE fCycle{1};
   union {
      struct {
         RUInt32BE fSeekKey{0};
         RUInt32BE fSeekPdir{0};
      } fInfoShort;
      struct {
         RUInt64BE fSeekKey{0};
         RUInt64BE fSeekPdir{0};
      } fInfoLong;
   };

   std::uint32_t fKeyHeaderSize{18 + sizeof(fInfoShort)};  // not part of serialization

   RTFKey() : fInfoShort() {}
   RTFKey(std::uint64_t seekKey, std::uint64_t seekPdir,
          const RTFString &clName, const RTFString &objName, const RTFString &titleName,
          std::uint32_t szObjInMem, std::uint32_t szObjOnDisk = 0)
   {
      fObjLen = szObjInMem;
      if ((seekKey > static_cast<unsigned int>(std::numeric_limits<std::int32_t>::max())) ||
          (seekPdir > static_cast<unsigned int>(std::numeric_limits<std::int32_t>::max())))
      {
         fKeyHeaderSize = 18 + sizeof(fInfoLong);
         fKeyLen = fKeyHeaderSize + clName.GetSize() + objName.GetSize() + titleName.GetSize();
         fInfoLong.fSeekKey = seekKey;
         fInfoLong.fSeekPdir = seekPdir;
         fVersion = fVersion + 1000;
      } else {
         fKeyHeaderSize = 18 + sizeof(fInfoShort);
         fKeyLen = fKeyHeaderSize + clName.GetSize() + objName.GetSize() + titleName.GetSize();
         fInfoShort.fSeekKey = seekKey;
         fInfoShort.fSeekPdir = seekPdir;
      }
      fNbytes = fKeyLen + ((szObjOnDisk == 0) ? szObjInMem : szObjOnDisk);
   }

   std::uint32_t GetSize() const {
      // Negative size indicates a gap in the file
      if (fNbytes < 0)
         return -fNbytes;
      return fNbytes;
   }

   std::uint32_t GetHeaderSize() const {
      if (fVersion >= 1000)
         return 18 + sizeof(fInfoLong);
      return 18 + sizeof(fInfoShort);
   }

   std::uint64_t GetSeekKey() const {
      if (fVersion >= 1000)
         return fInfoLong.fSeekKey;
      return fInfoShort.fSeekKey;
   }
};

/// The TFile global header
struct RTFHeader {
   char fMagic[4]{ 'r', 'o', 'o', 't' };
   RUInt32BE fVersion{(ROOT_VERSION_CODE >> 16)*10000 +
                      ((ROOT_VERSION_CODE & 0xFF00) >> 8) * 100 +
                      (ROOT_VERSION_CODE & 0xFF)};
   RUInt32BE fBEGIN{100};
   union {
      struct {
         RUInt32BE fEND{0};
         RUInt32BE fSeekFree{0};
         RUInt32BE fNbytesFree{0};
         RUInt32BE fNfree{1};
         RUInt32BE fNbytesName{0};
         unsigned char fUnits{4};
         RUInt32BE fCompress{0};
         RUInt32BE fSeekInfo{0};
         RUInt32BE fNbytesInfo{0};
      } fInfoShort;
      struct {
         RUInt64BE fEND{0};
         RUInt64BE fSeekFree{0};
         RUInt32BE fNbytesFree{0};
         RUInt32BE fNfree{1};
         RUInt32BE fNbytesName{0};
         unsigned char fUnits{8};
         RUInt32BE fCompress{0};
         RUInt64BE fSeekInfo{0};
         RUInt32BE fNbytesInfo{0};
      } fInfoLong;
   };

   RTFHeader() : fInfoShort() {}
   RTFHeader(int compression) : fInfoShort() {
      fInfoShort.fCompress = compression;
   }

   void SetBigFile() {
      if (fVersion >= 1000000)
         return;

      std::uint32_t end        = fInfoShort.fEND;
      std::uint32_t seekFree   = fInfoShort.fSeekFree;
      std::uint32_t nbytesFree = fInfoShort.fNbytesFree;
      std::uint32_t nFree      = fInfoShort.fNfree;
      std::uint32_t nbytesName = fInfoShort.fNbytesName;
      std::uint32_t compress   = fInfoShort.fCompress;
      std::uint32_t seekInfo   = fInfoShort.fSeekInfo;
      std::uint32_t nbytesInfo = fInfoShort.fNbytesInfo;
      fInfoLong.fEND        = end;
      fInfoLong.fSeekFree   = seekFree;
      fInfoLong.fNbytesFree = nbytesFree;
      fInfoLong.fNfree      = nFree;
      fInfoLong.fNbytesName = nbytesName;
      fInfoLong.fUnits      = 8;
      fInfoLong.fCompress   = compress;
      fInfoLong.fSeekInfo   = seekInfo;
      fInfoLong.fNbytesInfo = nbytesInfo;
      fVersion = fVersion + 1000000;
   }

   bool IsBigFile(std::uint64_t offset = 0) const {
      return (fVersion >= 1000000) || (offset > static_cast<unsigned int>(std::numeric_limits<std::int32_t>::max()));
   }

   std::uint32_t GetSize() const {
      std::uint32_t sizeHead = 4 + sizeof(fVersion) + sizeof(fBEGIN);
      if (IsBigFile()) return sizeHead + sizeof(fInfoLong);
      return sizeHead + sizeof(fInfoShort);
   }

   std::uint64_t GetEnd() const {
      if (IsBigFile()) return fInfoLong.fEND;
      return fInfoShort.fEND;
   }

   void SetEnd(std::uint64_t value) {
      if (IsBigFile(value)) {
         SetBigFile();
         fInfoLong.fEND = value;
      } else {
         fInfoShort.fEND = value;
      }
   }

   std::uint64_t GetSeekFree() const {
      if (IsBigFile()) return fInfoLong.fSeekFree;
      return fInfoShort.fSeekFree;
   }

   void SetSeekFree(std::uint64_t value) {
      if (IsBigFile(value)) {
         SetBigFile();
         fInfoLong.fSeekFree = value;
      } else {
         fInfoShort.fSeekFree = value;
      }
   }

   void SetNbytesFree(std::uint32_t value) {
      if (IsBigFile()) {
         fInfoLong.fNbytesFree = value;
      } else {
         fInfoShort.fNbytesFree = value;
      }
   }

   void SetNbytesName(std::uint32_t value) {
      if (IsBigFile()) {
         fInfoLong.fNbytesName = value;
      } else {
         fInfoShort.fNbytesName = value;
      }
   }

   std::uint64_t GetSeekInfo() const {
      if (IsBigFile()) return fInfoLong.fSeekInfo;
      return fInfoShort.fSeekInfo;
   }

   void SetSeekInfo(std::uint64_t value) {
      if (IsBigFile(value)) {
         SetBigFile();
         fInfoLong.fSeekInfo = value;
      } else {
         fInfoShort.fSeekInfo = value;
      }
   }

   void SetNbytesInfo(std::uint32_t value) {
      if (IsBigFile()) {
         fInfoLong.fNbytesInfo = value;
      } else {
         fInfoShort.fNbytesInfo = value;
      }
   }

   void SetCompression(std::uint32_t value) {
      if (IsBigFile()) {
         fInfoLong.fCompress = value;
      } else {
         fInfoShort.fCompress = value;
      }
   }
};


/// A reference to an unused byte-range in a TFile
struct RTFFreeEntry {
   RUInt16BE fVersion{1};
   union {
      struct {
         RUInt32BE fFirst{0};
         RUInt32BE fLast{0};
      } fInfoShort;
      struct {
         RUInt64BE fFirst{0};
         RUInt64BE fLast{0};
      } fInfoLong;
   };

   RTFFreeEntry() : fInfoShort() {}
   void Set(std::uint64_t first, std::uint64_t last) {
      if (last > static_cast<unsigned int>(std::numeric_limits<std::int32_t>::max())) {
         fVersion = fVersion + 1000;
         fInfoLong.fFirst = first;
         fInfoLong.fLast = last;
      } else {
         fInfoShort.fFirst = first;
         fInfoShort.fLast = last;
      }
   }
   std::uint32_t GetSize() { return (fVersion >= 1000) ? 18 : 10; }
};

/// Streamer info for TObject
struct RTFObject {
   RUInt16BE fVersion{1};
   RUInt32BE fUniqueID{0};  // unused
   RUInt32BE fBits;
   explicit RTFObject(std::uint32_t bits) : fBits(bits) {}
};

/// Streamer info for data member RNTuple::fVersion
struct RTFStreamerElementVersion {
   RUInt32BE fByteCount{0x40000000 | (sizeof(RTFStreamerElementVersion) - sizeof(RUInt32BE))};
   RUInt16BE fVersion{4};

   RUInt32BE fByteCountNamed{0x40000000 | (sizeof(RUInt16BE) + sizeof(RTFObject) + 10)};
   RUInt16BE fVersionNamed{1};
   RTFObject fObjectNamed{0x02000000 | 0x01000000};
   char fLName = 8;
   char fName[8]{ 'f', 'V', 'e', 'r', 's', 'i', 'o', 'n' };
   char fLTitle = 0;

   RUInt32BE fType{13};
   RUInt32BE fSize{4};
   RUInt32BE fArrLength{0};
   RUInt32BE fArrDim{0};
   char fMaxIndex[20]{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   char fLTypeName = 12;
   char fTypeName[12]{ 'u', 'n', 's', 'i', 'g', 'n', 'e', 'd', ' ', 'i', 'n', 't' };
};

/// Streamer info for data member RNTuple::fSize
struct RTFStreamerElementSize {
   RUInt32BE fByteCount{0x40000000 | (sizeof(RTFStreamerElementSize) - sizeof(RUInt32BE))};
   RUInt16BE fVersion{4};

   RUInt32BE fByteCountNamed{0x40000000 | (sizeof(RUInt16BE) + sizeof(RTFObject) + 7)};
   RUInt16BE fVersionNamed{1};
   RTFObject fObjectNamed{0x02000000 | 0x01000000};
   char fLName = 5;
   char fName[5]{ 'f', 'S', 'i', 'z', 'e' };
   char fLTitle = 0;

   RUInt32BE fType{13};
   RUInt32BE fSize{4};
   RUInt32BE fArrLength{0};
   RUInt32BE fArrDim{0};
   char fMaxIndex[20]{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   char fLTypeName = 12;
   char fTypeName[12]{ 'u', 'n', 's', 'i', 'g', 'n', 'e', 'd', ' ', 'i', 'n', 't' };
};

/// Streamer info for data member RNTuple::fSeekHeader
struct RTFStreamerElementSeekHeader {
   RUInt32BE fByteCount{0x40000000 | (sizeof(RTFStreamerElementSeekHeader) - sizeof(RUInt32BE))};
   RUInt16BE fVersion{4};

   RUInt32BE fByteCountNamed{0x40000000 |
      (sizeof(RUInt16BE) + sizeof(RTFObject) + 13)};
   RUInt16BE fVersionNamed{1};
   RTFObject fObjectNamed{0x02000000 | 0x01000000};
   char fLName = 11;
   char fName[11]{ 'f', 'S', 'e', 'e', 'k', 'H', 'e', 'a', 'd', 'e', 'r' };
   char fLTitle = 0;

   RUInt32BE fType{14};
   RUInt32BE fSize{8};
   RUInt32BE fArrLength{0};
   RUInt32BE fArrDim{0};
   char fMaxIndex[20]{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   char fLTypeName = 13;
   char fTypeName[13]{ 'u', 'n', 's', 'i', 'g', 'n', 'e', 'd', ' ', 'l', 'o', 'n', 'g' };
};

/// Streamer info for data member RNTuple::fNbytesHeader
struct RTFStreamerElementNBytesHeader {
   RUInt32BE fByteCount{0x40000000 | (sizeof(RTFStreamerElementNBytesHeader) - sizeof(RUInt32BE))};
   RUInt16BE fVersion{4};

   RUInt32BE fByteCountNamed{0x40000000 |
      (sizeof(RUInt16BE) + sizeof(RTFObject) + 15)};
   RUInt16BE fVersionNamed{1};
   RTFObject fObjectNamed{0x02000000 | 0x01000000};
   char fLName = 13;
   char fName[13]{ 'f', 'N', 'B', 'y', 't', 'e', 's', 'H', 'e', 'a', 'd', 'e', 'r' };
   char fLTitle = 0;

   RUInt32BE fType{13};
   RUInt32BE fSize{4};
   RUInt32BE fArrLength{0};
   RUInt32BE fArrDim{0};
   char fMaxIndex[20]{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   char fLTypeName = 12;
   char fTypeName[12]{ 'u', 'n', 's', 'i', 'g', 'n', 'e', 'd', ' ', 'i', 'n', 't' };
};

/// Streamer info for data member RNTuple::fLenHeader
struct RTFStreamerElementLenHeader {
   RUInt32BE fByteCount{0x40000000 | (sizeof(RTFStreamerElementLenHeader) - sizeof(RUInt32BE))};
   RUInt16BE fVersion{4};

   RUInt32BE fByteCountNamed{0x40000000 |
      (sizeof(RUInt16BE) + sizeof(RTFObject) + 12)};
   RUInt16BE fVersionNamed{1};
   RTFObject fObjectNamed{0x02000000 | 0x01000000};
   char fLName = 10;
   char fName[10]{ 'f', 'L', 'e', 'n', 'H', 'e', 'a', 'd', 'e', 'r' };
   char fLTitle = 0;

   RUInt32BE fType{13};
   RUInt32BE fSize{4};
   RUInt32BE fArrLength{0};
   RUInt32BE fArrDim{0};
   char fMaxIndex[20]{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   char fLTypeName = 12;
   char fTypeName[12]{ 'u', 'n', 's', 'i', 'g', 'n', 'e', 'd', ' ', 'i', 'n', 't' };
};

/// Streamer info for data member RNTuple::fSeekFooter
struct RTFStreamerElementSeekFooter {
   RUInt32BE fByteCount{0x40000000 | (sizeof(RTFStreamerElementSeekFooter) - sizeof(RUInt32BE))};
   RUInt16BE fVersion{4};

   RUInt32BE fByteCountNamed{0x40000000 |
      (sizeof(RUInt16BE) + sizeof(RTFObject) + 13)};
   RUInt16BE fVersionNamed{1};
   RTFObject fObjectNamed{0x02000000 | 0x01000000};
   char fLName = 11;
   char fName[11]{ 'f', 'S', 'e', 'e', 'k', 'F', 'o', 'o', 't', 'e', 'r' };
   char fLTitle = 0;

   RUInt32BE fType{14};
   RUInt32BE fSize{8};
   RUInt32BE fArrLength{0};
   RUInt32BE fArrDim{0};
   char fMaxIndex[20]{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   char fLTypeName = 13;
   char fTypeName[13]{ 'u', 'n', 's', 'i', 'g', 'n', 'e', 'd', ' ', 'l', 'o', 'n', 'g' };
};

/// Streamer info for data member RNTuple::fNbytesFooter
struct RTFStreamerElementNBytesFooter {
   RUInt32BE fByteCount{0x40000000 | (sizeof(RTFStreamerElementNBytesFooter) - sizeof(RUInt32BE))};
   RUInt16BE fVersion{4};

   RUInt32BE fByteCountNamed{0x40000000 | (sizeof(RUInt16BE) + sizeof(RTFObject) + 15)};
   RUInt16BE fVersionNamed{1};
   RTFObject fObjectNamed{0x02000000 | 0x01000000};
   char fLName = 13;
   char fName[13]{ 'f', 'N', 'B', 'y', 't', 'e', 's', 'F', 'o', 'o', 't', 'e', 'r' };
   char fLTitle = 0;

   RUInt32BE fType{13};
   RUInt32BE fSize{4};
   RUInt32BE fArrLength{0};
   RUInt32BE fArrDim{0};
   char fMaxIndex[20]{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   char fLTypeName = 12;
   char fTypeName[12]{ 'u', 'n', 's', 'i', 'g', 'n', 'e', 'd', ' ', 'i', 'n', 't' };
};

/// Streamer info for data member RNTuple::fLenFooter
struct RTFStreamerElementLenFooter {
   RUInt32BE fByteCount{0x40000000 | (sizeof(RTFStreamerElementLenFooter) - sizeof(RUInt32BE))};
   RUInt16BE fVersion{4};

   RUInt32BE fByteCountNamed{0x40000000 | (sizeof(RUInt16BE) + sizeof(RTFObject) + 12)};
   RUInt16BE fVersionNamed{1};
   RTFObject fObjectNamed{0x02000000 | 0x01000000};
   char fLName = 10;
   char fName[10]{ 'f', 'L', 'e', 'n', 'F', 'o', 'o', 't', 'e', 'r' };
   char fLTitle = 0;

   RUInt32BE fType{13};
   RUInt32BE fSize{4};
   RUInt32BE fArrLength{0};
   RUInt32BE fArrDim{0};
   char fMaxIndex[20]{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   char fLTypeName = 12;
   char fTypeName[12]{ 'u', 'n', 's', 'i', 'g', 'n', 'e', 'd', ' ', 'i', 'n', 't' };
};

/// Streamer info for data member RNTuple::fReserved
struct RTFStreamerElementReserved {
   RUInt32BE fByteCount{0x40000000 | (sizeof(RTFStreamerElementReserved) - sizeof(RUInt32BE))};
   RUInt16BE fVersion{4};

   RUInt32BE fByteCountNamed{0x40000000 |
      (sizeof(RUInt16BE) + sizeof(RTFObject) + 11)};
   RUInt16BE fVersionNamed{1};
   RTFObject fObjectNamed{0x02000000 | 0x01000000};
   char fLName = 9;
   char fName[9]{ 'f', 'R', 'e', 's', 'e', 'r', 'v', 'e', 'd' };
   char fLTitle = 0;

   RUInt32BE fType{14};
   RUInt32BE fSize{8};
   RUInt32BE fArrLength{0};
   RUInt32BE fArrDim{0};
   char fMaxIndex[20]{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
   char fLTypeName = 13;
   char fTypeName[13]{ 'u', 'n', 's', 'i', 'g', 'n', 'e', 'd', ' ', 'l', 'o', 'n', 'g' };
};

/// Streamer info frame for data member RNTuple::fVersion
struct RTFStreamerVersion {
   RUInt32BE fByteCount{0x40000000 | (sizeof(RTFStreamerVersion) - sizeof(RUInt32BE))};
   RUInt32BE fNewClassTag{0xffffffff};
   char fClassName[19]{'T', 'S', 't', 'r', 'e', 'a', 'm', 'e', 'r', 'B', 'a', 's', 'i', 'c', 'T', 'y', 'p', 'e', '\0'};
   RUInt32BE fByteCountRemaining{0x40000000 |
      (sizeof(RTFStreamerVersion) - 2 * sizeof(RUInt32BE) - 19 /* strlen(fClassName) + 1 */ - sizeof(RUInt32BE))};
   RUInt16BE fVersion{2};
   RTFStreamerElementVersion fStreamerElementVersion;
};

/// Streamer info frame for data member RNTuple::fSize
struct RTFStreamerSize {
   RUInt32BE fByteCount{0x40000000 | (sizeof(RTFStreamerSize) - sizeof(RUInt32BE))};
   RUInt32BE fClassTag{0x80000000};  // Fix-up after construction, or'd with 0x80000000
   RUInt32BE fByteCountRemaining{0x40000000 | (sizeof(RTFStreamerSize) - 3 * sizeof(RUInt32BE))};
   RUInt16BE fVersion{2};
   RTFStreamerElementSize fStreamerElementSize;
};

/// Streamer info frame for data member RNTuple::fSeekHeader
struct RTFStreamerSeekHeader {
   RUInt32BE fByteCount{0x40000000 | (sizeof(RTFStreamerSeekHeader) - sizeof(RUInt32BE))};
   RUInt32BE fClassTag{0x80000000};  // Fix-up after construction, or'd with 0x80000000
   RUInt32BE fByteCountRemaining{0x40000000 | (sizeof(RTFStreamerSeekHeader) - 3 * sizeof(RUInt32BE))};
   RUInt16BE fVersion{2};
   RTFStreamerElementSeekHeader fStreamerElementSeekHeader;
};

/// Streamer info frame for data member RNTuple::fNbytesHeader
struct RTFStreamerNBytesHeader {
   RUInt32BE fByteCount{0x40000000 | (sizeof(RTFStreamerNBytesHeader) - sizeof(RUInt32BE))};
   RUInt32BE fClassTag{0x80000000};  // Fix-up after construction, or'd with 0x80000000
   RUInt32BE fByteCountRemaining{0x40000000 | (sizeof(RTFStreamerNBytesHeader) - 3 * sizeof(RUInt32BE))};
   RUInt16BE fVersion{2};
   RTFStreamerElementNBytesHeader fStreamerElementNBytesHeader;
};

/// Streamer info frame for data member RNTuple::fLenHeader
struct RTFStreamerLenHeader {
   RUInt32BE fByteCount{0x40000000 | (sizeof(RTFStreamerLenHeader) - sizeof(RUInt32BE))};
   RUInt32BE fClassTag{0x80000000};  // Fix-up after construction, or'd with 0x80000000
   RUInt32BE fByteCountRemaining{0x40000000 | (sizeof(RTFStreamerLenHeader) - 3 * sizeof(RUInt32BE))};
   RUInt16BE fVersion{2};
   RTFStreamerElementLenHeader fStreamerElementLenHeader;
};

/// Streamer info frame for data member RNTuple::fSeekFooter
struct RTFStreamerSeekFooter {
   RUInt32BE fByteCount{0x40000000 | (sizeof(RTFStreamerSeekFooter) - sizeof(RUInt32BE))};
   RUInt32BE fClassTag{0x80000000};  // Fix-up after construction, or'd with 0x80000000
   RUInt32BE fByteCountRemaining{0x40000000 | (sizeof(RTFStreamerSeekFooter) - 3 * sizeof(RUInt32BE))};
   RUInt16BE fVersion{2};
   RTFStreamerElementSeekFooter fStreamerElementSeekFooter;
};

/// Streamer info frame for data member RNTuple::fNBytesFooter
struct RTFStreamerNBytesFooter {
   RUInt32BE fByteCount{0x40000000 | (sizeof(RTFStreamerNBytesFooter) - sizeof(RUInt32BE))};
   RUInt32BE fClassTag{0x80000000};  // Fix-up after construction, or'd with 0x80000000
   RUInt32BE fByteCountRemaining{0x40000000 | (sizeof(RTFStreamerNBytesFooter) - 3 * sizeof(RUInt32BE))};
   RUInt16BE fVersion{2};
   RTFStreamerElementNBytesFooter fStreamerElementNBytesFooter;
};

/// Streamer info frame for data member RNTuple::fLenFooter
struct RTFStreamerLenFooter {
   RUInt32BE fByteCount{0x40000000 | (sizeof(RTFStreamerLenFooter) - sizeof(RUInt32BE))};
   RUInt32BE fClassTag{0x80000000};  // Fix-up after construction, or'd with 0x80000000
   RUInt32BE fByteCountRemaining{0x40000000 | (sizeof(RTFStreamerLenFooter) - 3 * sizeof(RUInt32BE))};
   RUInt16BE fVersion{2};
   RTFStreamerElementLenFooter fStreamerElementLenFooter;
};

/// Streamer info frame for data member RNTuple::fReserved
struct RTFStreamerReserved {
   RUInt32BE fByteCount{0x40000000 | (sizeof(RTFStreamerReserved) - sizeof(RUInt32BE))};
   RUInt32BE fClassTag{0x80000000};  // Fix-up after construction, or'd with 0x80000000
   RUInt32BE fByteCountRemaining{0x40000000 | (sizeof(RTFStreamerReserved) - 3 * sizeof(RUInt32BE))};
   RUInt16BE fVersion{2};
   RTFStreamerElementReserved fStreamerElementReserved;
};

/// Streamer info for class RNTuple
struct RTFStreamerInfoObject {
   RUInt32BE fByteCount{0x40000000 | (sizeof(RTFStreamerInfoObject) - sizeof(fByteCount))};
   RUInt32BE fNewClassTag{0xffffffff};
   char fClassName[14]{ 'T', 'S', 't', 'r', 'e', 'a', 'm', 'e', 'r', 'I', 'n', 'f', 'o', '\0' };
   RUInt32BE fByteCountRemaining{0x40000000 |
      (sizeof(RTFStreamerInfoObject) - 2 * sizeof(RUInt32BE) - 14 - sizeof(RUInt32BE))};
   RUInt16BE fVersion{9};

   RUInt32BE fByteCountNamed{0x40000000 |
      (sizeof(RUInt16BE) + sizeof(RTFObject) + 29 /* strlen("ROOT::Experimental::RNTuple") + 2 */)};
   RUInt16BE fVersionNamed{1};
   RTFObject fObjectNamed{0x02000000 | 0x01000000 | 0x00010000};
   char fLName = 27;
   char fName[27]{ 'R', 'O', 'O', 'T', ':', ':',
      'E', 'x', 'p', 'e', 'r', 'i', 'm', 'e', 'n', 't', 'a', 'l', ':', ':',
      'R', 'N', 'T', 'u', 'p', 'l', 'e'};
   char fLTitle = 0;

   RInt32BE fChecksum{ChecksumRNTupleClass()};
   RUInt32BE fVersionRNTuple{1};

   RUInt32BE fByteCountObjArr{0x40000000 |
      (sizeof(RUInt32BE) + 10 /* strlen(TObjArray) + 1 */ + sizeof(RUInt32BE) +
       sizeof(RUInt16BE) + sizeof(RTFObject) + 1 + 2*sizeof(RUInt32BE) +
       sizeof(fStreamers))};
   RUInt32BE fNewClassTagObjArray{0xffffffff};
   char fClassNameObjArray[10]{'T', 'O', 'b', 'j', 'A', 'r', 'r', 'a', 'y', '\0'};
   RUInt32BE fByteCountObjArrRemaining{0x40000000 |
      (sizeof(RUInt16BE) + sizeof(RTFObject) + 1 + 2*sizeof(RUInt32BE) +
       sizeof(fStreamers))};
   RUInt16BE fVersionObjArr{3};
   RTFObject fObjectObjArr{0x02000000};
   char fNameObjArr{0};

   RUInt32BE fNObjects{9};
   RUInt32BE fLowerBound{0};

   struct {
      RTFStreamerVersion fStreamerVersion;
      RTFStreamerSize fStreamerSize;
      RTFStreamerSeekHeader fStreamerSeekHeader;
      RTFStreamerNBytesHeader fStreamerNBytesHeader;
      RTFStreamerLenHeader fStreamerLenHeader;
      RTFStreamerSeekFooter fStreamerSeekFooter;
      RTFStreamerNBytesFooter fStreamerNBytesFooter;
      RTFStreamerLenFooter fStreamerLenFooter;
      RTFStreamerReserved fStreamerReserved;
   } fStreamers;
};

/// The list of streamer info objects, for a new ntuple contains only the RNTuple class
struct RTFStreamerInfoList {
   RUInt32BE fByteCount{0x40000000 | (sizeof(RTFStreamerInfoList) - sizeof(fByteCount))};
   RUInt16BE fVersion{5};
   RTFObject fObject{0x02000000};
   char fName{0};
   RUInt32BE fNObjects{1};
   RTFStreamerInfoObject fStreamerInfo;
   char fEnd{0};

   std::uint32_t GetSize() const { return sizeof(RTFStreamerInfoList); }
};

/// The header of the directory key index
struct RTFKeyList {
   RUInt32BE fNKeys;
   std::uint32_t GetSize() const { return sizeof(RTFKeyList); }
   explicit RTFKeyList(std::uint32_t nKeys) : fNKeys(nKeys) {}
};

/// A streamed TFile object
struct RTFFile {
   char fModified{0};
   char fWritable{1};
   RTFDatetime fDateC;
   RTFDatetime fDateM;
   RUInt32BE fNBytesKeys{0};
   RUInt32BE fNBytesName{0};
   // The version of the key has to tell whether offsets are 32bit or 64bit long
   union {
      struct {
         RUInt32BE fSeekDir{100};
         RUInt32BE fSeekParent{0};
         RUInt32BE fSeekKeys{0};
      } fInfoShort;
      struct {
         RUInt64BE fSeekDir{100};
         RUInt64BE fSeekParent{0};
         RUInt64BE fSeekKeys{0};
      } fInfoLong;
   };

   RTFFile() : fInfoShort() {}

   std::uint32_t GetSize(std::uint32_t versionKey = 0) const {
      if (versionKey >= 1000)
         return 18 + sizeof(fInfoLong);
      return 18 + sizeof(fInfoShort);
   }

   std::uint64_t GetSeekKeys(std::uint32_t versionKey = 0) const {
      if (versionKey >= 1000)
         return fInfoLong.fSeekKeys;
      return fInfoShort.fSeekKeys;
   }
};

/// A streamed RNTuple class
struct RTFNTuple {
   RUInt32BE fByteCount{0x40000000 | (sizeof(RTFNTuple) - sizeof(fByteCount))};
   RUInt16BE fVersionClass{0};
   RInt32BE fChecksum{ChecksumRNTupleClass()};
   RUInt32BE fVersionInternal{0};
   RUInt32BE fSize{sizeof(ROOT::Experimental::RNTuple)};
   RUInt64BE fSeekHeader{0};
   RUInt32BE fNBytesHeader{0};
   RUInt32BE fLenHeader{0};
   RUInt64BE fSeekFooter{0};
   RUInt32BE fNBytesFooter{0};
   RUInt32BE fLenFooter{0};
   RUInt64BE fReserved{0};

   RTFNTuple() = default;
   explicit RTFNTuple(const ROOT::Experimental::RNTuple &inMemoryNTuple) {
      fVersionInternal = inMemoryNTuple.fVersion;
      fSize            = inMemoryNTuple.fSize;
      fSeekHeader      = inMemoryNTuple.fSeekHeader;
      fNBytesHeader    = inMemoryNTuple.fNBytesHeader;
      fLenHeader       = inMemoryNTuple.fLenHeader;
      fSeekFooter      = inMemoryNTuple.fSeekFooter;
      fNBytesFooter    = inMemoryNTuple.fNBytesFooter;
      fLenFooter       = inMemoryNTuple.fLenFooter;
      fReserved        = inMemoryNTuple.fReserved;
   }
   std::uint32_t GetSize() const { return sizeof(RTFNTuple); }
   ROOT::Experimental::RNTuple ToRNTuple() const {
      ROOT::Experimental::RNTuple ntuple;
      ntuple.fVersion      = fVersionInternal;
      ntuple.fSize         = fSize;
      ntuple.fSeekHeader   = fSeekHeader;
      ntuple.fNBytesHeader = fNBytesHeader;
      ntuple.fLenHeader    = fLenHeader;
      ntuple.fSeekFooter   = fSeekFooter;
      ntuple.fNBytesFooter = fNBytesFooter;
      ntuple.fLenFooter    = fLenFooter;
      ntuple.fReserved     = fReserved;
      return ntuple;
   }
   void Print(std::ostream& output) {
      output << "RTFNTuple {\n";
      output << "    fByteCount: " << fByteCount << ",\n";
      output << "    fVersionClass: " << fVersionClass << ",\n";
      output << "    fChecksum: " << fChecksum << ",\n";
      output << "    fVersionInternal: " << fVersionInternal << ",\n";
      output << "    fSize: " << fSize << ",\n";
      output << "    fSeekHeader: " << fSeekHeader << ",\n";
      output << "    fNBytesHeader: " << fNBytesHeader << ",\n";
      output << "    fLenHeader: " << fLenHeader << ",\n";
      output << "    fSeekFooter: " << fSeekFooter << ",\n";
      output << "    fNBytesFooter: " << fNBytesFooter << ",\n";
      output << "    fLenFooter: " << fLenFooter << ",\n";
      output << "    fReserved: " << fReserved << ",\n";
      output << "}";
   }
};

/// The bare file global header
struct RBareFileHeader {
   char fMagic[7]{ 'r', 'n', 't', 'u', 'p', 'l', 'e' };
   RUInt32BE fRootVersion{(ROOT_VERSION_CODE >> 16) * 10000 +
                          ((ROOT_VERSION_CODE & 0xFF00) >> 8) * 100 +
                          (ROOT_VERSION_CODE & 0xFF)};
   RUInt32BE fFormatVersion{1};
   RUInt32BE fCompress{0};
   RTFNTuple fNTuple;
   // followed by the ntuple name
};
#pragma pack(pop)

/// The artifical class name shown for opaque RNTuple keys (see TBasket)
static constexpr char const *kBlobClassName = "RBlob";
/// The class name of the RNTuple anchor
static constexpr char const *kNTupleClassName = "ROOT::Experimental::RNTuple";

/// The RKeyBlob writes an invsisble key into a TFile.  That is, a key that is not indexed in the list of keys,
/// like a TBasket.
class RKeyBlob : public TKey {
public:
   explicit RKeyBlob(TFile *file) : TKey(file) {
      fClassName = kBlobClassName;
      fKeylen += strlen(kBlobClassName);
   }

   /// Register a new key for a data record of size nbytes
   void Reserve(size_t nbytes, std::uint64_t *seekKey)
   {
      Create(nbytes);
      *seekKey = fSeekKey;
   }
};

} // anonymous namespace


namespace ROOT {
namespace Experimental {
namespace Internal {
/// If a TFile container is written by a C stream (simple file), on dataset commit, the file header needs to be updated
/// and the RNTuple object needs to be written at the offset that was reserved on file creation
struct RTFileControlBlock {
   RTFHeader fHeader;
   std::uint32_t fSeekNTuple{0};
};
} // namespace ROOT
} // namespace Experimental
} // namespace Internal


ROOT::Experimental::Internal::RMiniFileReader::RMiniFileReader(ROOT::Internal::RRawFile *rawFile)
   : fRawFile(rawFile)
{
}

ROOT::Experimental::RResult<ROOT::Experimental::RNTuple>
ROOT::Experimental::Internal::RMiniFileReader::GetNTuple(std::string_view ntupleName)
{
   char ident[4];
   ReadBuffer(ident, 4, 0);
   if (std::string(ident, 4) == "root")
      return GetNTupleProper(ntupleName);
   fIsBare = true;
   return GetNTupleBare(ntupleName);
}


ROOT::Experimental::RResult<ROOT::Experimental::RNTuple>
ROOT::Experimental::Internal::RMiniFileReader::GetNTupleProper(std::string_view ntupleName)
{
   RTFHeader fileHeader;
   ReadBuffer(&fileHeader, sizeof(fileHeader), 0);

   RTFKey key;
   RTFString name;
   ReadBuffer(&key, sizeof(key), fileHeader.fBEGIN);
   auto offset = fileHeader.fBEGIN + key.fKeyLen;
   ReadBuffer(&name, 1, offset);
   offset += name.GetSize();
   ReadBuffer(&name, 1, offset);
   offset += name.GetSize();
   RTFFile file;
   ReadBuffer(&file, sizeof(file), offset);

   RUInt32BE nKeys;
   offset = file.GetSeekKeys(key.fVersion);
   ReadBuffer(&key, sizeof(key), offset);
   offset += key.fKeyLen;
   ReadBuffer(&nKeys, sizeof(nKeys), offset);
   offset += sizeof(nKeys);
   bool found = false;
   for (unsigned int i = 0; i < nKeys; ++i) {
      ReadBuffer(&key, sizeof(key), offset);
      auto offsetNextKey = offset + key.fKeyLen;

      offset += key.GetHeaderSize();
      ReadBuffer(&name, 1, offset);
      ReadBuffer(&name, name.GetSize(), offset);
      if (std::string_view(name.fData, name.fLName) != kNTupleClassName) {
         offset = offsetNextKey;
         continue;
      }
      offset += name.GetSize();
      ReadBuffer(&name, 1, offset);
      ReadBuffer(&name, name.GetSize(), offset);
      if (std::string_view(name.fData, name.fLName) == ntupleName) {
         found = true;
         break;
      }
      offset = offsetNextKey;
   }
   if (!found) {
      return R__FAIL("no RNTuple named '" + std::string(ntupleName)
         + "' in file '" + fRawFile->GetUrl() + "'");
   }

   ReadBuffer(&key, sizeof(key), key.GetSeekKey());
   offset = key.GetSeekKey() + key.fKeyLen;
   RTFNTuple ntuple;
   ReadBuffer(&ntuple, sizeof(ntuple), offset);
   return ntuple.ToRNTuple();
}

ROOT::Experimental::RResult<ROOT::Experimental::RNTuple>
ROOT::Experimental::Internal::RMiniFileReader::GetNTupleBare(std::string_view ntupleName)
{
   RBareFileHeader fileHeader;
   ReadBuffer(&fileHeader, sizeof(fileHeader), 0);
   RTFString name;
   auto offset = sizeof(fileHeader);
   ReadBuffer(&name, 1, offset);
   ReadBuffer(&name, name.GetSize(), offset);
   std::string_view foundName(name.fData, name.fLName);
   if (foundName != ntupleName) {
      return R__FAIL("expected RNTuple named '" + std::string(ntupleName)
         + "' but instead found '" + std::string(foundName)
         + "' in file '" + fRawFile->GetUrl() + "'");
   }
   offset += name.GetSize();

   RTFNTuple ntuple;
   ReadBuffer(&ntuple, sizeof(ntuple), offset);
   return ntuple.ToRNTuple();
}


void ROOT::Experimental::Internal::RMiniFileReader::ReadBuffer(void *buffer, size_t nbytes, std::uint64_t offset)
{
   auto nread = fRawFile->ReadAt(buffer, nbytes, offset);
   R__ASSERT(nread == nbytes);
}


////////////////////////////////////////////////////////////////////////////////


ROOT::Experimental::Internal::RNTupleFileWriter::RFileSimple::~RFileSimple()
{
   if (fFile)
      fclose(fFile);
}


void ROOT::Experimental::Internal::RNTupleFileWriter::RFileSimple::Write(
   const void *buffer, size_t nbytes, std::int64_t offset)
{
   R__ASSERT(fFile);
   size_t retval;
   if ((offset >= 0) && (static_cast<std::uint64_t>(offset) != fFilePos)) {
#ifdef R__SEEK64
      retval = fseeko64(fFile, offset, SEEK_SET);
#else
      retval = fseek(fFile, offset, SEEK_SET);
#endif
      R__ASSERT(retval == 0);
      fFilePos = offset;
   }
   retval = fwrite(buffer, 1, nbytes, fFile);
   R__ASSERT(retval == nbytes);
   fFilePos += nbytes;
}


std::uint64_t ROOT::Experimental::Internal::RNTupleFileWriter::RFileSimple::WriteKey(
   const void *buffer, std::size_t nbytes, std::size_t len, std::int64_t offset,
   std::uint64_t directoryOffset,
   const std::string &className,
   const std::string &objectName,
   const std::string &title)
{
   if (offset < 0)
      offset = fFilePos;
   RTFString strClass{className};
   RTFString strObject{objectName};
   RTFString strTitle{title};

   RTFKey key(offset, directoryOffset, strClass, strObject, strTitle, len, nbytes);
   Write(&key, key.fKeyHeaderSize, offset);
   Write(&strClass, strClass.GetSize());
   Write(&strObject, strObject.GetSize());
   Write(&strTitle, strTitle.GetSize());
   auto offsetData = fFilePos;
   if (buffer)
      Write(buffer, nbytes);

   return offsetData;
}


////////////////////////////////////////////////////////////////////////////////


void ROOT::Experimental::Internal::RNTupleFileWriter::RFileProper::Write(
   const void *buffer, size_t nbytes, std::int64_t offset)
{
   R__ASSERT(fFile);
   fFile->Seek(offset);
   bool rv = fFile->WriteBuffer((char *)(buffer), nbytes);
   R__ASSERT(!rv);
}


std::uint64_t ROOT::Experimental::Internal::RNTupleFileWriter::RFileProper::WriteKey(
   const void *buffer, size_t nbytes, size_t len)
{
   std::uint64_t offsetKey;
   RKeyBlob keyBlob(fFile);
   keyBlob.Reserve(nbytes, &offsetKey);

   auto offset = offsetKey;
   RTFString strClass{kBlobClassName};
   RTFString strObject;
   RTFString strTitle;
   RTFKey keyHeader(offset, offset, strClass, strObject, strTitle, len, nbytes);

   Write(&keyHeader, keyHeader.fKeyHeaderSize, offset);
   offset += keyHeader.fKeyHeaderSize;
   Write(&strClass, strClass.GetSize(), offset);
   offset += strClass.GetSize();
   Write(&strObject, strObject.GetSize(), offset);
   offset += strObject.GetSize();
   Write(&strTitle, strTitle.GetSize(), offset);
   offset += strTitle.GetSize();
   auto offsetData = offset;
   Write(buffer, nbytes, offset);

   return offsetData;
}


////////////////////////////////////////////////////////////////////////////////


ROOT::Experimental::Internal::RNTupleFileWriter::RNTupleFileWriter(std::string_view name)
   : fNTupleName(name)
{
   fFileSimple.fControlBlock = std::make_unique<ROOT::Experimental::Internal::RTFileControlBlock>();
}


ROOT::Experimental::Internal::RNTupleFileWriter::~RNTupleFileWriter()
{
}


ROOT::Experimental::Internal::RNTupleFileWriter *ROOT::Experimental::Internal::RNTupleFileWriter::Recreate(
   std::string_view ntupleName, std::string_view path, int defaultCompression, ENTupleContainerFormat containerFormat)
{
   std::string fileName(path);
   size_t idxDirSep = fileName.find_last_of("\\/");
   if (idxDirSep != std::string::npos) {
      fileName.erase(0, idxDirSep + 1);
   }
#ifdef R__SEEK64
   FILE *fileStream = fopen64(std::string(path.data(), path.size()).c_str(), "wb");
#else
   FILE *fileStream = fopen(std::string(path.data(), path.size()).c_str(), "wb");
#endif
   R__ASSERT(fileStream);

   auto writer = new RNTupleFileWriter(ntupleName);
   writer->fFileSimple.fFile = fileStream;
   writer->fFileName = fileName;

   switch (containerFormat) {
   case ENTupleContainerFormat::kTFile:
      writer->WriteTFileSkeleton(defaultCompression);
      break;
   case ENTupleContainerFormat::kBare:
      writer->fIsBare = true;
      writer->WriteBareFileSkeleton(defaultCompression);
      break;
   default:
      R__ASSERT(false && "Internal error: unhandled container format");
   }

   return writer;
}


ROOT::Experimental::Internal::RNTupleFileWriter *ROOT::Experimental::Internal::RNTupleFileWriter::Recreate(
   std::string_view ntupleName, std::string_view path, std::unique_ptr<TFile> &file)
{
   file = std::unique_ptr<TFile>(TFile::Open(std::string(path.data(), path.size()).c_str(), "RECREATE"));
   R__ASSERT(file && !file->IsZombie());

   auto writer = new RNTupleFileWriter(ntupleName);
   writer->fFileProper.fFile = file.get();
   return writer;
}


ROOT::Experimental::Internal::RNTupleFileWriter *ROOT::Experimental::Internal::RNTupleFileWriter::Append(
   std::string_view ntupleName, TFile &file)
{
   auto writer = new RNTupleFileWriter(ntupleName);
   writer->fFileProper.fFile = &file;
   return writer;
}


void ROOT::Experimental::Internal::RNTupleFileWriter::Commit()
{
   if (fFileProper) {
      // Easy case, the ROOT file header and the RNTuple streaming is taken care of by TFile
      fFileProper.fFile->WriteObject(&fNTupleAnchor, fNTupleName.c_str());
      fFileProper.fFile->Write();
      return;
   }

   // Writing by C file stream: prepare the container format header and stream the RNTuple anchor object
   R__ASSERT(fFileSimple);
   RTFNTuple ntupleOnDisk(fNTupleAnchor);

   if (fIsBare) {
      fFileSimple.Write(&ntupleOnDisk, ntupleOnDisk.GetSize(), fFileSimple.fControlBlock->fSeekNTuple);
      fflush(fFileSimple.fFile);
      return;
   }

   fFileSimple.fControlBlock->fHeader.SetSeekFree(fFileSimple.fFilePos);
   RTFString strEmpty;
   RTFFreeEntry freeEntry;
   RTFKey keyFreeList(fFileSimple.fControlBlock->fHeader.GetSeekFree(), 100,
                      strEmpty, strEmpty, strEmpty, freeEntry.GetSize());
   std::uint64_t firstFree = fFileSimple.fControlBlock->fHeader.GetSeekFree() + keyFreeList.GetSize();
   freeEntry.Set(firstFree, std::max(2000000000ULL, ((firstFree / 1000000000ULL) + 1) * 1000000000ULL));
   fFileSimple.WriteKey(
      &freeEntry, freeEntry.GetSize(), freeEntry.GetSize(), fFileSimple.fControlBlock->fHeader.GetSeekFree(),
      100, "", "", "");
   fFileSimple.fControlBlock->fHeader.SetNbytesFree(fFileSimple.fFilePos -
                                                    fFileSimple.fControlBlock->fHeader.GetSeekFree());
   fFileSimple.fControlBlock->fHeader.SetEnd(fFileSimple.fFilePos);

   auto szNTuple = ntupleOnDisk.GetSize();
   fFileSimple.WriteKey(&ntupleOnDisk, szNTuple, szNTuple, fFileSimple.fControlBlock->fSeekNTuple, 100,
                        kNTupleClassName, fNTupleName);

   fFileSimple.Write(&fFileSimple.fControlBlock->fHeader, fFileSimple.fControlBlock->fHeader.GetSize(), 0);
   fflush(fFileSimple.fFile);
}

int ROOT::Experimental::Internal::RNTupleFileWriter::GetTFileCompression() const {
   R__ASSERT(fFileProper && "Internal error: called GetTFileCompression on a bare file");
   return fFileProper.fFile->GetCompressionSettings();
}

std::uint64_t ROOT::Experimental::Internal::RNTupleFileWriter::WriteBlob(const void *data, size_t nbytes, size_t len)
{
   std::uint64_t offset;
   if (fFileSimple) {
      if (fIsBare) {
         offset = fFileSimple.fFilePos;
         fFileSimple.Write(data, nbytes);
      } else {
         offset = fFileSimple.WriteKey(data, nbytes, len, -1, 100, kBlobClassName);
      }
   } else {
      offset = fFileProper.WriteKey(data, nbytes, len);
   }
   return offset;
}


std::uint64_t ROOT::Experimental::Internal::RNTupleFileWriter::WriteNTupleHeader(
   const void *data, size_t nbytes, size_t lenHeader)
{
   auto offset = WriteBlob(data, nbytes, lenHeader);
   fNTupleAnchor.fLenHeader = lenHeader;
   fNTupleAnchor.fNBytesHeader = nbytes;
   fNTupleAnchor.fSeekHeader = offset;
   return offset;
}


std::uint64_t ROOT::Experimental::Internal::RNTupleFileWriter::WriteNTupleFooter(
   const void *data, size_t nbytes, size_t lenFooter)
{
   auto offset = WriteBlob(data, nbytes, lenFooter);
   fNTupleAnchor.fLenFooter = lenFooter;
   fNTupleAnchor.fNBytesFooter = nbytes;
   fNTupleAnchor.fSeekFooter = offset;
   return offset;
}


void ROOT::Experimental::Internal::RNTupleFileWriter::WriteBareFileSkeleton(int defaultCompression)
{
   RBareFileHeader bareHeader;
   bareHeader.fCompress = defaultCompression;
   fFileSimple.Write(&bareHeader, sizeof(bareHeader), 0);
   RTFString ntupleName{fNTupleName};
   fFileSimple.Write(&ntupleName, ntupleName.GetSize());

   // Write zero-initialized ntuple to reserve the space; will be overwritten on commit
   RTFNTuple ntupleOnDisk;
   fFileSimple.fControlBlock->fSeekNTuple = fFileSimple.fFilePos;
   fFileSimple.Write(&ntupleOnDisk, ntupleOnDisk.GetSize());
}


void ROOT::Experimental::Internal::RNTupleFileWriter::WriteTFileSkeleton(int defaultCompression)
{
   RTFString strTFile{"TFile"};
   RTFString strFileName{fFileName};
   RTFString strTList{"TList"};
   RTFString strStreamerInfo{"StreamerInfo"};
   RTFString strStreamerTitle{"Doubly linked list"};
   RTFString strRNTupleClass{"ROOT::Experimental::RNTuple"};
   RTFString strRNTupleName{fNTupleName};
   RTFString strEmpty;

   fFileSimple.fControlBlock->fHeader = RTFHeader(defaultCompression);

   // First record of the file: the TFile object at offset 100
   RTFFile fileRoot;
   RTFKey keyRoot(100, 0, strTFile, strFileName, strEmpty,
                  fileRoot.GetSize() + strFileName.GetSize() + strEmpty.GetSize());
   std::uint32_t nbytesName = keyRoot.fKeyLen + strFileName.GetSize() + 1;
   fileRoot.fNBytesName = nbytesName;
   fFileSimple.fControlBlock->fHeader.SetNbytesName(nbytesName);

   // Second record: the compressed StreamerInfo with the description of the RNTuple class.
   // This record usually comes at the end, but in this context we create a file with a single RNTuple only,
   // so we know beforehand all the types.
   fFileSimple.fControlBlock->fHeader.SetSeekInfo(100 + keyRoot.GetSize());
   RTFKey keyStreamerInfo(
      fFileSimple.fControlBlock->fHeader.GetSeekInfo(), 100, strTList, strStreamerInfo, strStreamerTitle, 0);
   RTFStreamerInfoList streamerInfo;
   auto classTagOffset = keyStreamerInfo.fKeyLen +
      offsetof(struct RTFStreamerInfoList, fStreamerInfo) +
      offsetof(struct RTFStreamerInfoObject, fStreamers) +
      offsetof(struct RTFStreamerVersion, fNewClassTag) + 2;
   streamerInfo.fStreamerInfo.fStreamers.fStreamerSize.fClassTag = 0x80000000 | classTagOffset;
   streamerInfo.fStreamerInfo.fStreamers.fStreamerSeekHeader.fClassTag = 0x80000000 | classTagOffset;
   streamerInfo.fStreamerInfo.fStreamers.fStreamerNBytesHeader.fClassTag = 0x80000000 | classTagOffset;
   streamerInfo.fStreamerInfo.fStreamers.fStreamerLenHeader.fClassTag = 0x80000000 | classTagOffset;
   streamerInfo.fStreamerInfo.fStreamers.fStreamerSeekFooter.fClassTag = 0x80000000 | classTagOffset;
   streamerInfo.fStreamerInfo.fStreamers.fStreamerNBytesFooter.fClassTag = 0x80000000 | classTagOffset;
   streamerInfo.fStreamerInfo.fStreamers.fStreamerLenFooter.fClassTag = 0x80000000 | classTagOffset;
   streamerInfo.fStreamerInfo.fStreamers.fStreamerReserved.fClassTag = 0x80000000 | classTagOffset;
   Detail::RNTupleCompressor compressor;
   auto szStreamerInfo = compressor.Zip(&streamerInfo, streamerInfo.GetSize(), 1);
   fFileSimple.WriteKey(compressor.GetZipBuffer(), szStreamerInfo, streamerInfo.GetSize(),
                        fFileSimple.fControlBlock->fHeader.GetSeekInfo(), 100,
                        "TList", "StreamerInfo", "Doubly linked list");
   fFileSimple.fControlBlock->fHeader.SetNbytesInfo(
      fFileSimple.fFilePos - fFileSimple.fControlBlock->fHeader.GetSeekInfo());

   // Reserve the space for the RNTuple record, which will be written on commit
   RTFNTuple ntupleOnDisk;
   fFileSimple.fControlBlock->fSeekNTuple = fFileSimple.fFilePos;
   RTFKey keyRNTuple(fFileSimple.fControlBlock->fSeekNTuple, 100, strRNTupleClass, strRNTupleName, strEmpty,
                     ntupleOnDisk.GetSize());
   fFileSimple.WriteKey(&ntupleOnDisk, ntupleOnDisk.GetSize(), ntupleOnDisk.GetSize(),
                        fFileSimple.fControlBlock->fSeekNTuple, 100, "ROOT::Experimental::RNTUple", fNTupleName, "");

   // The key index of the root TFile object, containing for the time being only the RNTuple key
   fileRoot.fInfoShort.fSeekKeys = fFileSimple.fFilePos;
   RTFKeyList keyList{1};
   RTFKey keyKeyList(fileRoot.GetSeekKeys(), 100, strEmpty, strEmpty, strEmpty, keyList.GetSize() + keyRNTuple.fKeyLen);
   fFileSimple.Write(&keyKeyList, keyKeyList.fKeyHeaderSize, fileRoot.GetSeekKeys());
   fFileSimple.Write(&strEmpty, strEmpty.GetSize());
   fFileSimple.Write(&strEmpty, strEmpty.GetSize());
   fFileSimple.Write(&strEmpty, strEmpty.GetSize());
   fFileSimple.Write(&keyList, keyList.GetSize());
   fFileSimple.Write(&keyRNTuple, keyRNTuple.fKeyHeaderSize);
   fFileSimple.Write(&strRNTupleClass, strRNTupleClass.GetSize());
   fFileSimple.Write(&strRNTupleName, strRNTupleName.GetSize());
   fFileSimple.Write(&strEmpty, strEmpty.GetSize());
   fileRoot.fNBytesKeys = fFileSimple.fFilePos - fileRoot.GetSeekKeys();

   auto tail = fFileSimple.fFilePos;
   fFileSimple.Write(&keyRoot, keyRoot.fKeyHeaderSize, 100);
   fFileSimple.Write(&strTFile, strTFile.GetSize());
   fFileSimple.Write(&strFileName, strFileName.GetSize());
   fFileSimple.Write(&strEmpty, strEmpty.GetSize());
   fFileSimple.Write(&strFileName, strFileName.GetSize());
   fFileSimple.Write(&strEmpty, strEmpty.GetSize());
   fFileSimple.Write(&fileRoot, fileRoot.GetSize());
   fFileSimple.fFilePos = tail;
#ifdef R__SEEK64
   auto retval = fseeko64(fFileSimple.fFile, tail, SEEK_SET);
#else
   auto retval = fseek(fFileSimple.fFile, tail, SEEK_SET);
#endif
   R__ASSERT(retval == 0);
   fFileSimple.fFilePos = tail;
}
