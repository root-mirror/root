/// \file RNTupleDescriptorFmt.cxx
/// \ingroup NTuple ROOT7
/// \author Jakob Blomer <jblomer@cern.ch>
/// \date 2019-08-25
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
#include <ROOT/RColumnModel.hxx>
#include <ROOT/RNTupleDescriptor.hxx>
#include <ROOT/RNTupleUtil.hxx>

#include <algorithm>
#include <iomanip>
#include <ostream>
#include <unordered_map>
#include <vector>

namespace {

struct ClusterInfo {
   std::uint64_t fFirstEntry = 0;
   std::uint32_t fNPages = 0;
   std::uint32_t fNEntries = 0;
   std::uint32_t fBytesOnStorage = 0;
   std::uint32_t fBytesInMemory = 0;

   bool operator ==(const ClusterInfo &other) const {
      return fFirstEntry == other.fFirstEntry;
   }

   bool operator <(const ClusterInfo &other) const {
      return fFirstEntry < other.fFirstEntry;
   }
};

struct ColumnInfo {
   ROOT::Experimental::DescriptorId_t fColumnId = 0;
   ROOT::Experimental::DescriptorId_t fFieldId = 0;
   std::uint64_t fLocalOrder = 0;
   std::uint64_t fNElements = 0;
   std::uint64_t fNPages = 0;
   std::uint64_t fBytesOnStorage = 0;
   std::uint32_t fElementSize = 0;
   ROOT::Experimental::EColumnType fType;
   std::string fFieldName;

   bool operator <(const ColumnInfo &other) const {
      if (fFieldName == other.fFieldName)
         return fLocalOrder < other.fLocalOrder;
      return fFieldName < other.fFieldName;
   }
};

static std::string GetFieldName(ROOT::Experimental::DescriptorId_t fieldId,
   const ROOT::Experimental::RNTupleDescriptor &ntupleDesc)
{
   const auto &fieldDesc = ntupleDesc.GetFieldDescriptor(fieldId);
   if (fieldDesc.GetParentId() == ROOT::Experimental::kInvalidDescriptorId)
      return fieldDesc.GetFieldName();
   return GetFieldName(fieldDesc.GetParentId(), ntupleDesc) + "." + fieldDesc.GetFieldName();
}

static std::string GetColumnTypeName(ROOT::Experimental::EColumnType type)
{
   switch (type) {
   case ROOT::Experimental::EColumnType::kBit:
      return "Bit";
   case ROOT::Experimental::EColumnType::kByte:
      return "Byte";
   case ROOT::Experimental::EColumnType::kInt32:
      return "Int32";
   case ROOT::Experimental::EColumnType::kInt64:
      return "Int64";
   case ROOT::Experimental::EColumnType::kReal32:
      return "Real32";
   case ROOT::Experimental::EColumnType::kReal64:
      return "Real64";
   case ROOT::Experimental::EColumnType::kIndex:
      return "Index";
   case ROOT::Experimental::EColumnType::kSwitch:
      return "Switch";
   default:
      return "UNKNOWN";
   }
}

} // anonymous namespace

void ROOT::Experimental::RNTupleDescriptor::PrintInfo(std::ostream &output) const
{
   std::vector<ColumnInfo> columns;
   std::vector<ClusterInfo> clusters;
   std::unordered_map<DescriptorId_t, unsigned int> cluster2Idx;
   for (const auto &cluster : fClusterDescriptors) {
      ClusterInfo info;
      info.fFirstEntry = cluster.second.GetFirstEntryIndex();
      info.fNEntries = cluster.second.GetNEntries();
      cluster2Idx[cluster.first] = clusters.size();
      clusters.emplace_back(info);
   }

   std::uint64_t bytesOnStorage = 0;
   std::uint64_t bytesInMemory = 0;
   std::uint64_t nPages = 0;
   int compression = -1;
   for (const auto &column : fColumnDescriptors) {
      // We generate the default memory representation for the given column type in order
      // to report the size _in memory_ of column elements
      auto elementSize = Detail::RColumnElementBase::Generate(column.second.GetModel().GetType())->GetSize();

      ColumnInfo info;
      info.fColumnId = column.second.GetId();
      info.fFieldId = column.second.GetFieldId();
      info.fLocalOrder = column.second.GetIndex();
      info.fElementSize = elementSize;
      info.fType = column.second.GetModel().GetType();

      for (const auto &cluster : fClusterDescriptors) {
         auto columnRange = cluster.second.GetColumnRange(column.first);
         info.fNElements += columnRange.fNElements;
         if (compression == -1) {
            compression = columnRange.fCompressionSettings;
         }
         const auto &pageRange = cluster.second.GetPageRange(column.first);
         auto idx = cluster2Idx[cluster.first];
         for (const auto &page : pageRange.fPageInfos) {
            bytesOnStorage += page.fLocator.fBytesOnStorage;
            bytesInMemory += page.fNElements * elementSize;
            clusters[idx].fBytesOnStorage += page.fLocator.fBytesOnStorage;
            clusters[idx].fBytesInMemory += page.fNElements * elementSize;
            ++clusters[idx].fNPages;
            info.fBytesOnStorage += page.fLocator.fBytesOnStorage;
            ++info.fNPages;
            ++nPages;
         }
      }
      columns.emplace_back(info);
   }
   auto headerSize = GetHeaderSize();
   auto footerSize = GetFooterSize();
   output << "============================================================" << std::endl;
   output << "NTUPLE:      " << GetName() << std::endl;
   output << "Compression: " << compression << std::endl;
   output << "------------------------------------------------------------" << std::endl;
   output << "  # Entries:        " << GetNEntries() << std::endl;
   output << "  # Fields:         " << GetNFields() << std::endl;
   output << "  # Columns:        " << GetNColumns() << std::endl;
   output << "  # Pages:          " << nPages << std::endl;
   output << "  # Clusters:       " << GetNClusters() << std::endl;
   output << "  Size on storage:  " << bytesOnStorage << " B" << std::endl;
   output << "  Compression rate: " << std::fixed << std::setprecision(2)
                                    << float(bytesInMemory) / float(bytesOnStorage) << std::endl;
   output << "  Header size:      " << headerSize << " B" << std::endl;
   output << "  Footer size:      " << footerSize << " B" << std::endl;
   output << "  Meta-data / data: " << std::fixed << std::setprecision(3)
                                    << float(headerSize + footerSize) / float(bytesOnStorage) << std::endl;
   output << "------------------------------------------------------------" << std::endl;
   output << "CLUSTER DETAILS" << std::endl;
   output << "------------------------------------------------------------" << std::endl;

   std::sort(clusters.begin(), clusters.end());
   for (unsigned int i = 0; i < clusters.size(); ++i) {
      output << "  # " << std::setw(5) << i
             << "   Entry range:     [" << clusters[i].fFirstEntry << ".."
             << clusters[i].fFirstEntry + clusters[i].fNEntries - 1 << "]  --  " << clusters[i].fNEntries << std::endl;
      output << "         "
             << "   # Pages:         " << clusters[i].fNPages << std::endl;
      output << "         "
             << "   Size on storage: " << clusters[i].fBytesOnStorage << " B" << std::endl;
      output << "         "
             << "   Compression:     " << std::fixed << std::setprecision(2)
             << float(clusters[i].fBytesInMemory) / float(float(clusters[i].fBytesOnStorage)) << std::endl;
   }

   output << "------------------------------------------------------------" << std::endl;
   output << "COLUMN DETAILS" << std::endl;
   output << "------------------------------------------------------------" << std::endl;
   for (auto &col : columns)
      col.fFieldName = GetFieldName(col.fFieldId, *this).substr(1);
   std::sort(columns.begin(), columns.end());
   for (const auto &col : columns) {
      auto avgPageSize = (col.fNPages == 0) ? 0 : (col.fBytesOnStorage / col.fNPages);
      auto avgElementsPerPage = (col.fNPages == 0) ? 0 : (col.fNElements / col.fNPages);
      std::string nameAndType = std::string("  ") + col.fFieldName + " [#" + std::to_string(col.fLocalOrder) + "]"
         + "  --  " + GetColumnTypeName(col.fType);
      std::string id = std::string("{id:") + std::to_string(col.fColumnId) + "}";
      output << nameAndType << std::setw(60 - nameAndType.length()) << id << std::endl;
      output << "    # Elements:          " << col.fNElements << std::endl;
      output << "    # Pages:             " << col.fNPages << std::endl;
      output << "    Avg elements / page: " << avgElementsPerPage << std::endl;
      output << "    Avg page size:       " << avgPageSize << " B" << std::endl;
      output << "    Size on storage:     " << col.fBytesOnStorage << " B" << std::endl;
      output << "    Compression:         " << std::fixed << std::setprecision(2)
             << float(col.fElementSize * col.fNElements) / float(col.fBytesOnStorage) << std::endl;
      output << "............................................................" << std::endl;
   }
}
