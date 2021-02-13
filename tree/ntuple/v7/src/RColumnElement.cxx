/// \file RColumnElement.cxx
/// \ingroup NTuple ROOT7
/// \author Jakob Blomer <jblomer@cern.ch>
/// \date 2019-08-11
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/RColumnElement.hxx>

#include <TError.h>

#include <algorithm>
#include <bitset>
#include <cstdint>
#include <memory>
#include <utility>

std::unique_ptr<ROOT::Experimental::Detail::RColumnElementBase>
ROOT::Experimental::Detail::RColumnElementBase::Generate(EColumnType type) {
   switch (type) {
   case EColumnType::kReal32:
      return std::make_unique<RColumnElement<float, EColumnType::kReal32>>(nullptr);
   case EColumnType::kReal64:
      return std::make_unique<RColumnElement<double, EColumnType::kReal64>>(nullptr);
   case EColumnType::kByte:
      return std::make_unique<RColumnElement<std::uint8_t, EColumnType::kByte>>(nullptr);
   case EColumnType::kInt32:
      return std::make_unique<RColumnElement<std::int32_t, EColumnType::kInt32>>(nullptr);
   case EColumnType::kInt64:
      return std::make_unique<RColumnElement<std::int64_t, EColumnType::kInt64>>(nullptr);
   case EColumnType::kBit:
      return std::make_unique<RColumnElement<bool, EColumnType::kBit>>(nullptr);
   case EColumnType::kIndex:
      return std::make_unique<RColumnElement<ClusterSize_t, EColumnType::kIndex>>(nullptr);
   case EColumnType::kSwitch:
      return std::make_unique<RColumnElement<RColumnSwitch, EColumnType::kSwitch>>(nullptr);
   default:
      R__ASSERT(false);
   }
   // never here
   return nullptr;
}

std::size_t ROOT::Experimental::Detail::RColumnElementBase::GetBitsOnStorage(EColumnType type) {
   switch (type) {
   case EColumnType::kReal32:
      return 32;
   case EColumnType::kReal64:
      return 64;
   case EColumnType::kByte:
      return 8;
   case EColumnType::kInt32:
      return 32;
   case EColumnType::kInt64:
      return 64;
   case EColumnType::kBit:
      return 1;
   case EColumnType::kIndex:
      return 32;
   case EColumnType::kSwitch:
      return 64;
   default:
      R__ASSERT(false);
   }
   // never here
   return 0;
}

void ROOT::Experimental::Detail::RColumnElement<bool, ROOT::Experimental::EColumnType::kBit>::Pack(
  void *dst, void *src, std::size_t count) const
{
   bool *boolArray = reinterpret_cast<bool *>(src);
   char *charArray = reinterpret_cast<char *>(dst);
   std::bitset<8> bitSet;
   std::size_t i = 0;
   for (; i < count; ++i) {
      bitSet.set(i % 8, boolArray[i]);
      if (i % 8 == 7) {
         char packed = bitSet.to_ulong();
         charArray[i / 8] = packed;
      }
   }
   if (i % 8 != 0) {
      char packed = bitSet.to_ulong();
      charArray[i / 8] = packed;
   }
}

void ROOT::Experimental::Detail::RColumnElement<bool, ROOT::Experimental::EColumnType::kBit>::Unpack(
  void *dst, void *src, std::size_t count) const
{
   bool *boolArray = reinterpret_cast<bool *>(dst);
   char *charArray = reinterpret_cast<char *>(src);
   std::bitset<8> bitSet;
   for (std::size_t i = 0; i < count; i += 8) {
      bitSet = charArray[i / 8];
      for (std::size_t j = i; j < std::min(count, i + 8); ++j) {
         boolArray[j] = bitSet[j % 8];
      }
   }
}
