/// \file ROOT/RNTupleDescriptor.hxx
/// \ingroup NTuple ROOT7
/// \author Jakob Blomer <jblomer@cern.ch>
/// \date 2018-07-19
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RNTupleDescriptor
#define ROOT7_RNTupleDescriptor

#include <ROOT/RColumnModel.hxx>
#include <ROOT/RNTupleUtil.hxx>
#include <ROOT/RStringView.hxx>

#include <chrono>
#include <memory>
#include <ostream>
#include <vector>
#include <string>
#include <unordered_map>

#include <iostream>
#include <iomanip>
namespace ROOT {
namespace Experimental {

class RNTupleDescriptorBuilder;
class RNTupleModel;

// clang-format off
/**
\class ROOT::Experimental::RFieldDescriptor
\ingroup NTuple
\brief Meta-data stored for every field of an ntuple
*/
// clang-format on
class RFieldDescriptor {
   friend class RNTupleDescriptorBuilder;

private:
   DescriptorId_t fFieldId = kInvalidDescriptorId;
   /// The version of the C++-type-to-column translation mechanics
   RNTupleVersion fFieldVersion;
   /// The version of the C++ type itself
   RNTupleVersion fTypeVersion;
   /// The leaf name, not including parent fields
   std::string fFieldName;
   /// Free text set by the user
   std::string fFieldDescription;
   /// The C++ type that was used when writing the field
   std::string fTypeName;
   /// The number of elements per entry for fixed-size arrays
   std::uint64_t fNRepetitions;
   /// The structural information carried by this field in the data model tree
   ENTupleStructure fStructure;
   /// Establishes sub field relationships, such as classes and collections
   DescriptorId_t fParentId = kInvalidDescriptorId;
   /// The pointers in the other direction from parent to children. They are serialized, too, to keep the
   /// order of sub fields.
   std::vector<DescriptorId_t> fLinkIds;

public:
   /// In order to handle changes to the serialization routine in future ntuple versions
   static constexpr std::uint16_t kFrameVersionCurrent = 0;
   static constexpr std::uint16_t kFrameVersionMin = 0;

   RFieldDescriptor() = default;
   RFieldDescriptor(const RFieldDescriptor &other) = delete;
   RFieldDescriptor &operator =(const RFieldDescriptor &other) = delete;
   RFieldDescriptor(RFieldDescriptor &&other) = default;
   RFieldDescriptor &operator =(RFieldDescriptor &&other) = default;

   bool operator==(const RFieldDescriptor &other) const;

   DescriptorId_t GetId() const { return fFieldId; }
   RNTupleVersion GetFieldVersion() const { return fFieldVersion; }
   RNTupleVersion GetTypeVersion() const { return fTypeVersion; }
   std::string GetFieldName() const { return fFieldName; }
   std::string GetFieldDescription() const { return fFieldDescription; }
   std::string GetTypeName() const { return fTypeName; }
   std::uint64_t GetNRepetitions() const { return fNRepetitions; }
   ENTupleStructure GetStructure() const { return fStructure; }
   DescriptorId_t GetParentId() const { return fParentId; }
   const std::vector<DescriptorId_t> &GetLinkIds() const { return fLinkIds; }
};


// clang-format off
/**
\class ROOT::Experimental::RColumnDescriptor
\ingroup NTuple
\brief Meta-data stored for every column of an ntuple
*/
// clang-format on
class RColumnDescriptor {
   friend class RNTupleDescriptorBuilder;

private:
   DescriptorId_t fColumnId = kInvalidDescriptorId;
   /// Versions can change, e.g., when new column types are added
   RNTupleVersion fVersion;
   /// Contains the column type and whether it is sorted
   RColumnModel fModel;
   /// Every column belongs to one and only one field
   DescriptorId_t fFieldId = kInvalidDescriptorId;
   /// A field can be serialized into several columns, which are numbered from zero to $n$
   std::uint32_t fIndex;

public:
   /// In order to handle changes to the serialization routine in future ntuple versions
   static constexpr std::uint16_t kFrameVersionCurrent = 0;
   static constexpr std::uint16_t kFrameVersionMin = 0;

   RColumnDescriptor() = default;
   RColumnDescriptor(const RColumnDescriptor &other) = delete;
   RColumnDescriptor &operator =(const RColumnDescriptor &other) = delete;
   RColumnDescriptor(RColumnDescriptor &&other) = default;
   RColumnDescriptor &operator =(RColumnDescriptor &&other) = default;

   bool operator==(const RColumnDescriptor &other) const;

   DescriptorId_t GetId() const { return fColumnId; }
   RNTupleVersion GetVersion() const { return fVersion; }
   RColumnModel GetModel() const { return fModel; }
   std::uint32_t GetIndex() const { return fIndex; }
   DescriptorId_t GetFieldId() const { return fFieldId; }
};


// clang-format off
/**
\class ROOT::Experimental::RClusterDescriptor
\ingroup NTuple
\brief Meta-data for a set of ntuple clusters

The cluster descriptor might carry information of only a subset of available clusters, for instance if multiple
files are chained and not all of them have been processed yet.
*/
// clang-format on
class RClusterDescriptor {
   friend class RNTupleDescriptorBuilder;

public:
   /// Generic information about the physical location of data. Values depend on the concrete storage type.  E.g.,
   /// for a local file fUrl might be unsused and fPosition might be a file offset. Objects on storage can be compressed
   /// and therefore we need to store their actual size.
   struct RLocator {
      std::int64_t fPosition = 0;
      std::uint32_t fBytesOnStorage = 0;
      std::string fUrl;

      bool operator==(const RLocator &other) const {
         return fPosition == other.fPosition && fBytesOnStorage == other.fBytesOnStorage && fUrl == other.fUrl;
      }
   };

   /// The window of element indexes of a particular column in a particular cluster
   struct RColumnRange {
      DescriptorId_t fColumnId = kInvalidDescriptorId;
      /// A 64bit element index
      NTupleSize_t fFirstElementIndex = kInvalidNTupleIndex;
      /// A 32bit value for the number of column elements in the cluster
      ClusterSize_t fNElements = kInvalidClusterIndex;
      /// The usual format for ROOT compression settings (see Compression.h).
      /// The pages of a particular column in a particular cluster are all compressed with the same settings.
      std::int64_t fCompressionSettings = 0;

      // TODO(jblomer): we perhaps want to store summary information, such as average, min/max, etc.
      // Should this be done on the field level?

      bool operator==(const RColumnRange &other) const {
         return fColumnId == other.fColumnId && fFirstElementIndex == other.fFirstElementIndex &&
                fNElements == other.fNElements && fCompressionSettings == other.fCompressionSettings;
      }

      bool Contains(NTupleSize_t index) const {
         return (fFirstElementIndex <= index && (fFirstElementIndex + fNElements) > index);
      }
   };

   /// Records the parition of data into pages for a particular column in a particular cluster
   struct RPageRange {
      /// We do not need to store the element size / uncompressed page size because we know to which column
      /// the page belongs
      struct RPageInfo {
         /// The sum of the elements of all the pages must match the corresponding fNElements field in fColumnRanges
         ClusterSize_t fNElements = kInvalidClusterIndex;
         /// The meaning of fLocator depends on the storage backend.
         RLocator fLocator;

         bool operator==(const RPageInfo &other) const {
            return fNElements == other.fNElements && fLocator == other.fLocator;
         }
      };

      RPageRange() = default;
      RPageRange(const RPageRange &other) = delete;
      RPageRange &operator =(const RPageRange &other) = delete;
      RPageRange(RPageRange &&other) = default;
      RPageRange &operator =(RPageRange &&other) = default;

      DescriptorId_t fColumnId = kInvalidDescriptorId;
      std::vector<RPageInfo> fPageInfos;

      bool operator==(const RPageRange &other) const {
         return fColumnId == other.fColumnId && fPageInfos == other.fPageInfos;
      }
   };

private:
   DescriptorId_t fClusterId = kInvalidDescriptorId;
   /// Future versions of the cluster descriptor might add more meta-data, e.g. a semantic checksum
   RNTupleVersion fVersion;
   /// Clusters can be swapped by adjusting the entry offsets
   NTupleSize_t fFirstEntryIndex = kInvalidNTupleIndex;
   ClusterSize_t fNEntries = kInvalidClusterIndex;
   /// For pre-fetching / caching an entire contiguous cluster
   RLocator fLocator;

   std::unordered_map<DescriptorId_t, RColumnRange> fColumnRanges;
   std::unordered_map<DescriptorId_t, RPageRange> fPageRanges;

public:
   /// In order to handle changes to the serialization routine in future ntuple versions
   static constexpr std::uint16_t kFrameVersionCurrent = 0;
   static constexpr std::uint16_t kFrameVersionMin = 0;

   RClusterDescriptor() = default;
   RClusterDescriptor(const RClusterDescriptor &other) = delete;
   RClusterDescriptor &operator =(const RClusterDescriptor &other) = delete;
   RClusterDescriptor(RClusterDescriptor &&other) = default;
   RClusterDescriptor &operator =(RClusterDescriptor &&other) = default;

   bool operator==(const RClusterDescriptor &other) const;
   /// Tells if 2 RPageSources can be befriended as in RPageSourceFriend
   bool BeFriendAble(const RClusterDescriptor &other) const;

   DescriptorId_t GetId() const { return fClusterId; }
   RNTupleVersion GetVersion() const { return fVersion; }
   NTupleSize_t GetFirstEntryIndex() const { return fFirstEntryIndex; }
   ClusterSize_t GetNEntries() const { return fNEntries; }
   RLocator GetLocator() const { return fLocator; }
   const RColumnRange &GetColumnRange(DescriptorId_t columnId) const { return fColumnRanges.at(columnId); }
   const RPageRange &GetPageRange(DescriptorId_t columnId) const { return fPageRanges.at(columnId); }
};


// clang-format off
/**
\class ROOT::Experimental::RNTupleDescriptor
\ingroup NTuple
\brief The on-storage meta-data of an ntuple

Represents the on-disk (on storage) information about an ntuple. The meta-data consists of a header and one or
several footers. The header carries the ntuple schema, i.e. the fields and the associated columns and their
relationships. The footer(s) carry information about one or several clusters. For every cluster, a footer stores
its location and size, and for every column the range of element indexes as well as a list of pages and page
locations.

The descriptor provide machine-independent (de-)serialization of headers and footers, and it provides lookup routines
for ntuple objects (pages, clusters, ...).  It is supposed to be usable by all RPageStorage implementations.

The serialization does not use standard ROOT streamers in order to not let it depend on libCore. The serialization uses
the concept of frames: header, footer, and substructures have a preamble with version numbers and the size of the
writte struct. This allows for forward and backward compatibility when the meta-data evolves.
*/
// clang-format on
class RNTupleDescriptor {
   friend class RNTupleDescriptorBuilder;

private:
   /// The ntuple name needs to be unique in a given storage location (file)
   std::string fName;
   /// Free text from the user
   std::string fDescription;
   /// The origin of the data
   std::string fAuthor;
   /// The current responsible for storing the data
   std::string fCustodian;
   /// The time stamp of the ntuple data (immutable)
   std::chrono::system_clock::time_point fTimeStampData;
   /// The time stamp of writing the data to storage, which gets updated when re-written
   std::chrono::system_clock::time_point fTimeStampWritten;
   /// The version evolves with the ntuple summary meta-data
   RNTupleVersion fVersion;
   /// Every NTuple gets a unique identifier
   RNTupleUuid fOwnUuid;
   /// Column sets that are created as derived sets from existing NTuples share the same group id.
   /// NTuples in the same group have the same number of entries and are supposed to contain associated data.
   RNTupleUuid fGroupUuid;

   std::unordered_map<DescriptorId_t, RFieldDescriptor> fFieldDescriptors;
   std::unordered_map<DescriptorId_t, RColumnDescriptor> fColumnDescriptors;
   /// May contain only a subset of all the available clusters, e.g. the clusters of the current file
   /// from a chain of files
   std::unordered_map<DescriptorId_t, RClusterDescriptor> fClusterDescriptors;

public:
   /// In order to handle changes to the serialization routine in future ntuple versions
   static constexpr std::uint16_t kFrameVersionCurrent = 0;
   static constexpr std::uint16_t kFrameVersionMin = 0;
   /// The preamble is sufficient to get the length of the header
   static constexpr unsigned int kNBytesPreamble = 8;
   /// The last few bytes after the footer store the length of footer and header
   static constexpr unsigned int kNBytesPostscript = 16;

   RNTupleDescriptor() = default;
   RNTupleDescriptor(const RNTupleDescriptor &other) = delete;
   RNTupleDescriptor &operator=(const RNTupleDescriptor &other) = delete;
   RNTupleDescriptor(RNTupleDescriptor &&other) = default;
   RNTupleDescriptor &operator=(RNTupleDescriptor &&other) = default;

   bool operator ==(const RNTupleDescriptor &other) const;

   /// We deliberately do not use ROOT's built-in serialization in order to allow for use of RNTuple's without libCore
   /// Serializes the global ntuple information as well as the column and field schemata
   /// Returns the number of bytes and fills buffer if it is not nullptr.
   std::uint32_t SerializeHeader(void* buffer) const;
   /// Serializes cluster meta data. Returns the number of bytes and fills buffer if it is not nullptr.
   std::uint32_t SerializeFooter(void* buffer) const;
   /// Given kNBytesPostscript bytes, extract the header and footer lengths in bytes
   static void LocateMetadata(const void *postscript, std::uint32_t &szHeader, std::uint32_t &szFooter);

   const RFieldDescriptor& GetFieldDescriptor(DescriptorId_t fieldId) const { return fFieldDescriptors.at(fieldId); }
   const RColumnDescriptor& GetColumnDescriptor(DescriptorId_t columnId) const {
      return fColumnDescriptors.at(columnId);
   }
   const RClusterDescriptor& GetClusterDescriptor(DescriptorId_t clusterId) const {
      return fClusterDescriptors.at(clusterId);
   }
   std::string GetName() const { return fName; }
   std::string GetDescription() const { return fDescription; }
   std::string GetAuthor() const { return fAuthor; }
   std::string GetCustodian() const { return fCustodian; }
   std::chrono::system_clock::time_point GetTimeStampData() const { return fTimeStampData; }
   std::chrono::system_clock::time_point GetTimeStampWritten() const { return fTimeStampWritten; }
   RNTupleVersion GetVersion() const { return fVersion; }
   RNTupleUuid GetOwnUuid() const { return fOwnUuid; }
   RNTupleUuid GetGroupUuid() const { return fGroupUuid; }

   std::size_t GetNFields() const { return fFieldDescriptors.size(); }
   std::size_t GetNColumns() const { return fColumnDescriptors.size(); }
   std::size_t GetNClusters() const { return fClusterDescriptors.size(); }

   // The number of entries as seen with the currently loaded cluster meta-data; there might be more
   NTupleSize_t GetNEntries() const;
   NTupleSize_t GetNElements(DescriptorId_t columnId) const;

   DescriptorId_t FindFieldId(std::string_view fieldName, DescriptorId_t parentId) const;
   /// Searches for a top-level field
   DescriptorId_t FindFieldId(std::string_view fieldName) const;
   DescriptorId_t FindColumnId(DescriptorId_t fieldId, std::uint32_t columnIndex) const;
   DescriptorId_t FindClusterId(DescriptorId_t columnId, NTupleSize_t index) const;

   /// Re-create the C++ model from the stored meta-data
   std::unique_ptr<RNTupleModel> GenerateModel() const;
   bool CommonFieldNameExists(const RNTupleDescriptor &desc) const;
   void PrintInfo(std::ostream &output) const;
};


// clang-format off
/**
\class ROOT::Experimental::RNTupleDescriptorBuilder
\ingroup NTuple
\brief A helper class for piece-wise construction of an RNTupleDescriptor

Used by RPageStorage implementations in order to construct the RNTupleDescriptor from the various header parts.
*/
// clang-format on
class RNTupleDescriptorBuilder {
private:
   RNTupleDescriptor fDescriptor;

public:
   bool IsValid() const { return true; /* TODO(jblomer) */}
   const RNTupleDescriptor& GetDescriptor() const { return fDescriptor; }
   RNTupleDescriptor MoveDescriptor();

   void SetNTuple(const std::string_view name, const std::string_view description, const std::string_view author,
                  const RNTupleVersion &version, const RNTupleUuid &uuid);

   void AddField(DescriptorId_t fieldId, const RNTupleVersion &fieldVersion, const RNTupleVersion &typeVersion,
                 std::string_view fieldName, std::string_view typeName, std::uint64_t nRepetitions,
                 ENTupleStructure structure);
   void AddFieldLink(DescriptorId_t fieldId, DescriptorId_t linkId);

   void AddColumn(DescriptorId_t columnId, DescriptorId_t fieldId,
                  const RNTupleVersion &version, const RColumnModel &model, std::uint32_t index);

   void SetFromHeader(void* headerBuffer);

   void AddCluster(DescriptorId_t clusterId, RNTupleVersion version,
                   NTupleSize_t firstEntryIndex, ClusterSize_t nEntries);
   void SetClusterLocator(DescriptorId_t clusterId, RClusterDescriptor::RLocator locator);
   void AddClusterColumnRange(DescriptorId_t clusterId, const RClusterDescriptor::RColumnRange &columnRange);
   void AddClusterPageRange(DescriptorId_t clusterId, RClusterDescriptor::RPageRange &&pageRange);

   void AddClustersFromFooter(void* footerBuffer);
   
   void SetNTuple(const RNTupleDescriptor& desc);
   void AddFieldsAndColumnsFromDescriptor(const RNTupleDescriptor& desc);
   void AddClustersFromDescriptor(const RNTupleDescriptor& desc);
   // specific for RPageSourceFriend
   void AddFriendInfo(const RNTupleDescriptor& desc);
};

} // namespace Experimental
} // namespace ROOT

#endif
