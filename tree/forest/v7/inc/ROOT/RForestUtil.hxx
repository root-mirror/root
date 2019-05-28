/// \file ROOT/RForestUtil.hxx
/// \ingroup Forest ROOT7
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

#ifndef ROOT7_RForestUtil
#define ROOT7_RForestUtil

#include <cstdint>

#include <string>
#include <vector>

namespace ROOT {
namespace Experimental {


/**
 * The fields in the data model tree can carry different structural information about the type system.
 * Leaf fields contain just data, collection fields resolve to offset columns, record root fields have no
 * materialization on the primitive column layer.
 */
enum EForestStructure {
  kLeaf,
  kCollection,
  kRecord,
  // unimplemented so far
  kReference,
  kOptional,
  kVariant,
};

/// Integer type long enough to hold the maximum number of entries in a column
using ForestSize_t = std::uint64_t;
constexpr ForestSize_t kInvalidForestIndex = std::uint64_t(-1);
/// Wrap the 32bit integer in a struct in order to avoid template specialization clash with std::uint32_t
struct RClusterSize {
   RClusterSize() : fValue(0) {}
   explicit constexpr RClusterSize(std::uint32_t value) : fValue(value) {}
   RClusterSize& operator =(const std::uint32_t value) { fValue = value; return *this; }
   RClusterSize& operator +=(const std::uint32_t value) { fValue += value; return *this; }
   RClusterSize operator++(int) { auto result = *this; fValue++; return result; }
   operator std::uint32_t() const { return fValue; }

   std::uint32_t fValue;
};
using ClusterSize_t = RClusterSize;
constexpr ClusterSize_t kInvalidClusterIndex(std::uint32_t(-1));

/// Uniquely identifies a physical column within the scope of the current process, used to tag pages
using ColumnId_t = std::int64_t;
constexpr ColumnId_t kInvalidColumnId = -1;

/// Distriniguishes elements of the same type within a descriptor, e.g. different fields
using DescriptorId_t = std::uint64_t;
constexpr DescriptorId_t kInvalidDescriptorId = std::uint64_t(-1);


/// 64 possible flags to apply to all versioned entities (so far unused).
using ForestFlags_t = std::uint64_t;
/// For forward and backward compatibility, attach version information to
/// the consitituents of the file format (column, field, cluster, forest).
class RForestVersion {
private:
   /// The version used to write an entity
   std::uint32_t fVersionUse = 0;
   /// The minimum required version necessary to read an entity
   std::uint32_t fVersionMin = 0;
   ForestFlags_t fFlags = 0;

public:
   RForestVersion() = default;
   RForestVersion(std::uint32_t versionUse, std::uint32_t versionMin)
     : fVersionUse(versionUse), fVersionMin(versionMin)
   {}
   RForestVersion(std::uint32_t versionUse, std::uint32_t versionMin, ForestFlags_t flags)
     : fVersionUse(versionUse), fVersionMin(versionMin), fFlags(flags)
   {}

   std::uint32_t GetVersionUse() const { return fVersionUse; }
   std::uint32_t GetVersionMin() const { return fVersionMin; }
   ForestFlags_t GetFlags() const { return fFlags; }
};

} // namespace Experimental
} // namespace ROOT

#endif
