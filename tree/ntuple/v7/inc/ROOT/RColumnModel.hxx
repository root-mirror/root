/// \file ROOT/RColumnModel.hxx
/// \ingroup NTuple ROOT7
/// \author Jakob Blomer <jblomer@cern.ch>
/// \date 2018-10-09
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RColumnModel
#define ROOT7_RColumnModel

#include <ROOT/RNTupleUtil.hxx>
#include <ROOT/RStringView.hxx>

#include <array>
#include <string>

namespace ROOT {
namespace Experimental {

// clang-format off
/**
\class ROOT::Experimental::EColumnType
\ingroup NTuple
\brief The available trivial, native content types of a column

More complex types, such as classes, get translated into columns of such simple types by the RField.
*/
// clang-format on
enum class EColumnType {
   kUnknown = 0,
   // type for root columns of (nested) collections; 32bit integers that count relative to the current cluster
   kIndex,
   // 64 bit column that uses the lower 32bits as kIndex and the higher 32bits as a dispatch tag; used, e.g.,
   // in order to serialize std::variant
   kSwitch,
   kByte,
   kBit,
   kReal64,
   kReal32,
   // kReal24, to uncomment after implementing custom-sized float
   kReal16,
   kReal8,
   kInt64,
   kInt32,
   kInt16,
   // kCustomDouble, to uncomment after implementing custom-sized float
   // kCustomFloat, to uncomment after implementing custom-sized float
};

// clang-format off
/**
\class ROOT::Experimental::RColumnTypeIdentifier
\ingroup NTuple
\brief Holds static arrays with EColumnType MetaData

Contains static arrays to obtain information about a specific columnType.
*/
// clang-format on
class RColumnTypeIdentifier {
public:
   static constexpr std::array<const char *, 12 /*15*/> fColumnTypeNames{
      "Unknown", "Index", "Switch", "Byte", "Bit", "Real64", "Real32", /*"Real24", */ "Real16",
      "Real8",   "Int64", "Int32",  "Int16" /*, "CustomDouble", "CustomFloat"*/};
   static constexpr std::array<ClusterSize_t::ValueType, 12 /*15*/> fColumnBitSizeOnDisk{
      0,
      sizeof(ClusterSize_t) * 8,
      sizeof(ROOT::Experimental::RColumnSwitch) * 8,
      sizeof(char) * 8,
      sizeof(bool) * 8,
      sizeof(double) * 8,
      sizeof(float) * 8,/*, 24*/
      16,
      8,
      64,
      32,
      16 /*, sizeof(double)*8, sizeof(float)*8*/};
   static const char *GetColumnTypeNames(std::uint32_t index)
   {
      assert(index < fColumnTypeNames.size());
      return fColumnTypeNames[index];
   }
   static ClusterSize_t::ValueType GetColumnBitSizeOnDisk(std::uint32_t index)
   {
      assert(index < fColumnBitSizeOnDisk.size());
      return fColumnBitSizeOnDisk[index];
   }
};

// clang-format off
/**
\class ROOT::Experimental::RColumnModel
\ingroup NTuple
\brief Holds the static meta-data of a column in a tree
*/
// clang-format on
class RColumnModel {
private:
   EColumnType fType;
   bool fIsSorted;

public:
   RColumnModel() : fType(EColumnType::kUnknown), fIsSorted(false) {}
   RColumnModel(EColumnType type, bool isSorted) : fType(type), fIsSorted(isSorted) {}

   EColumnType GetType() const { return fType; }
   bool GetIsSorted() const { return fIsSorted; }

   bool operator==(const RColumnModel &other) const { return (fType == other.fType) && (fIsSorted == other.fIsSorted); }
};

} // namespace Experimental
} // namespace ROOT

#endif
