/// \file RPageStorage.cxx
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

#include <ROOT/RField.hxx>
#include <ROOT/RNTupleDescriptor.hxx>
#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RPage.hxx>
#include <ROOT/RPageAllocator.hxx>
#include <ROOT/RPagePool.hxx>
#include <ROOT/RPageStorageRoot.hxx>
#include <ROOT/RLogger.hxx>

#include <TKey.h>

#include <cstdlib>
#include <iostream>
#include <utility>

namespace {

static constexpr const char* kKeyNTupleFooter = "NTPLF";
static constexpr const char* kKeyNTupleHeader = "NTPLH";

}

ROOT::Experimental::Detail::RPageSinkRoot::RPageSinkRoot(std::string_view ntupleName, RSettings settings)
   : RPageSink(ntupleName)
   , fPageAllocator(std::make_unique<RPageAllocatorHeap>())
   , fDirectory(nullptr)
   , fSettings(settings)
{
   R__WARNING_HERE("NTuple") << "The RNTuple file format will change. " <<
      "Do not store real data with this version of RNTuple!";
}

ROOT::Experimental::Detail::RPageSinkRoot::RPageSinkRoot(std::string_view ntupleName, std::string_view path)
   : RPageSink(ntupleName)
   , fPageAllocator(std::make_unique<RPageAllocatorHeap>())
   , fDirectory(nullptr)
{
   R__WARNING_HERE("NTuple") << "The RNTuple file format will change. " <<
      "Do not store real data with this version of RNTuple!";
   TFile *file = TFile::Open(std::string(path).c_str(), "UPDATE");
   fSettings.fFile = file;
   fSettings.fTakeOwnership = true;
}

ROOT::Experimental::Detail::RPageSinkRoot::~RPageSinkRoot()
{
   if (fSettings.fTakeOwnership) {
      fSettings.fFile->Close();
      delete fSettings.fFile;
   }
}

ROOT::Experimental::Detail::RPageStorage::ColumnHandle_t
ROOT::Experimental::Detail::RPageSinkRoot::AddColumn(DescriptorId_t fieldId, const RColumn &column)
{
   ROOT::Experimental::Internal::RColumnHeader columnHeader;
   columnHeader.fName = column.GetModel().GetName();
   columnHeader.fType = column.GetModel().GetType();
   columnHeader.fIsSorted = column.GetModel().GetIsSorted();
   if (column.GetOffsetColumn() != nullptr) {
      columnHeader.fOffsetColumn = column.GetOffsetColumn()->GetModel().GetName();
   }
   fNTupleHeader.fColumns.emplace_back(columnHeader);

   auto columnId = fLastColumnId++;
   fDescriptorBuilder.AddColumn(columnId, fieldId, column.GetVersion(), column.GetModel(), column.GetIndex());
   if (column.GetOffsetColumn() != nullptr) {
      auto parentId = column.GetOffsetColumn()->GetHandleSink().fId;
      fDescriptorBuilder.SetColumnOffset(columnId, parentId);
      fDescriptorBuilder.AddColumnLink(parentId, columnId);
   }

   //printf("Added column %s type %d\n", columnHeader.fName.c_str(), (int)columnHeader.fType);
   return ColumnHandle_t(columnId, &column);
}


void ROOT::Experimental::Detail::RPageSinkRoot::Create(RNTupleModel &model)
{
   fDirectory = fSettings.fFile->mkdir(fNTupleName.c_str());
   // In TBrowser, use RNTupleBrowser(TDirectory *directory) in order to show the ntuple contents
   fDirectory->SetBit(TDirectoryFile::kCustomBrowse);
   fDirectory->SetTitle("ROOT::Experimental::Detail::RNTupleBrowser");

   fDescriptorBuilder.SetNTuple(fNTupleName, model.GetDescription(), model.GetVersion(), model.GetUuid());

   std::unordered_map<const RFieldBase *, DescriptorId_t> fFieldPtr2Id; // necessary to find parent field ids
   for (auto& f : *model.GetRootField()) {
      ROOT::Experimental::Internal::RFieldHeader fieldHeader;
      fieldHeader.fName = f.GetName();
      fieldHeader.fType = f.GetType();
      //printf("Added field %s type [%s]\n", f.GetName().c_str(), f.GetType().c_str());
      if (f.GetParent()) fieldHeader.fParentName = f.GetParent()->GetName();
      fNTupleHeader.fFields.emplace_back(fieldHeader);
      fDescriptorBuilder.AddField(fLastFieldId, f.GetFieldVersion(), f.GetTypeVersion(),
                                  RFieldBase::GetLeafName(f.GetName()), f.GetType(), f.GetStructure());
      if (f.GetParent() != model.GetRootField()) {
         fDescriptorBuilder.SetFieldParent(fLastFieldId, fFieldPtr2Id[f.GetParent()]);
      }

      Detail::RFieldFuse::Connect(fLastFieldId, *this, f); // issues in turn one or several calls to AddColumn()
      fFieldPtr2Id[&f] = fLastFieldId++;
   }

   auto nColumns = fLastColumnId;
   fCurrentCluster.fPagesPerColumn.resize(nColumns);
   fNTupleFooter.fNElementsPerColumn.resize(nColumns, 0);
   for (DescriptorId_t i = 0; i < nColumns; ++i) {
      RClusterDescriptor::RColumnRange columnRange;
      columnRange.fColumnId = i;
      columnRange.fFirstElementIndex = 0;
      columnRange.fNElements = 0;
      fOpenColumnRanges.emplace_back(columnRange);
      RClusterDescriptor::RPageRange pageRange;
      pageRange.fColumnId = i;
      fOpenPageRanges.emplace_back(pageRange);
   }

   fDirectory->WriteObject(&fNTupleHeader, RMapper::kKeyNTupleHeader);
   const auto &descriptor = fDescriptorBuilder.GetDescriptor();
   auto szFooter = descriptor.SerializeHeader(nullptr);
   auto buffer = new unsigned char[szFooter];
   descriptor.SerializeHeader(buffer);
   ROOT::Experimental::Internal::RNTupleBlob blob(szFooter, buffer);
   fDirectory->WriteObject(&blob, kKeyNTupleHeader);
   delete[] buffer;
}

void ROOT::Experimental::Detail::RPageSinkRoot::CommitPage(ColumnHandle_t columnHandle, const RPage &page)
{
   ROOT::Experimental::Internal::RNTupleBlob pagePayload(
      page.GetSize(), static_cast<unsigned char *>(page.GetBuffer()));
   std::string keyName = std::string(RMapper::kKeyPagePayload) +
      std::to_string(fLastClusterId) + RMapper::kKeySeparator +
      std::to_string(fLastPageIdx);
   fDirectory->WriteObject(&pagePayload, keyName.c_str());

   auto columnId = columnHandle.fId;
   fCurrentCluster.fPagesPerColumn[columnId].fRangeStarts.push_back(page.GetRangeFirst());
   fNTupleFooter.fNElementsPerColumn[columnId] += page.GetNElements();

   fOpenColumnRanges[columnId].fNElements += page.GetNElements();
   RClusterDescriptor::RPageRange::RPageInfo pageInfo;
   pageInfo.fNElements = page.GetNElements();
   pageInfo.fLocator = fLastPageIdx++;
   fOpenPageRanges[columnId].fPageInfos.emplace_back(pageInfo);
}

void ROOT::Experimental::Detail::RPageSinkRoot::CommitCluster(ROOT::Experimental::NTupleSize_t nEntries)
{
   R__ASSERT((nEntries - fPrevClusterNEntries) < ClusterSize_t(-1));
   fDescriptorBuilder.AddCluster(fLastClusterId, RNTupleVersion(), fPrevClusterNEntries,
                                 ClusterSize_t(nEntries - fPrevClusterNEntries));
   for (auto &range : fOpenColumnRanges) {
      fDescriptorBuilder.AddClusterColumnRange(fLastClusterId, range);
      range.fFirstElementIndex += range.fNElements;
   }
   for (auto &range : fOpenPageRanges) {
      fDescriptorBuilder.AddClusterPageRange(fLastClusterId, range);
      range.fPageInfos.clear();
   }
   for (auto &range : fOpenColumnRanges) {
      range.fNElements = 0;
   }
   ++fLastClusterId;

   fCurrentCluster.fNEntries = nEntries - fPrevClusterNEntries;
   fPrevClusterNEntries = nEntries;
   std::string key = RMapper::kKeyClusterFooter + std::to_string(fNTupleFooter.fNClusters);
   fDirectory->WriteObject(&fCurrentCluster, key.c_str());
   fNTupleFooter.fNClusters++;
   fNTupleFooter.fNEntries = nEntries;

   for (auto& pageInfo : fCurrentCluster.fPagesPerColumn) {
      pageInfo.fRangeStarts.clear();
   }
   fCurrentCluster.fEntryRangeStart = fNTupleFooter.fNEntries;
   fLastPageIdx = 0;
}

void ROOT::Experimental::Detail::RPageSinkRoot::CommitDataset()
{
   if (!fDirectory)
      return;

   const auto &descriptor = fDescriptorBuilder.GetDescriptor();
   auto szFooter = descriptor.SerializeFooter(nullptr);
   auto buffer = new unsigned char[szFooter];
   descriptor.SerializeFooter(buffer);
   ROOT::Experimental::Internal::RNTupleBlob footerBlob(szFooter, buffer);
   fDirectory->WriteObject(&footerBlob, kKeyNTupleFooter);
   delete[] buffer;

   fDirectory->WriteObject(&fNTupleFooter, RMapper::kKeyNTupleFooter);
}

ROOT::Experimental::Detail::RPage
ROOT::Experimental::Detail::RPageSinkRoot::ReservePage(ColumnHandle_t columnHandle, std::size_t nElements)
{
   if (nElements == 0)
      nElements = kDefaultElementsPerPage;
   auto elementSize = columnHandle.fColumn->GetModel().GetElementSize();
   return fPageAllocator->NewPage(columnHandle.fId, elementSize, nElements);
}

void ROOT::Experimental::Detail::RPageSinkRoot::ReleasePage(RPage &page)
{
   fPageAllocator->DeletePage(page);
}


////////////////////////////////////////////////////////////////////////////////


ROOT::Experimental::Detail::RPage ROOT::Experimental::Detail::RPageAllocatorKey::NewPage(
   ColumnId_t columnId, void *mem, std::size_t elementSize, std::size_t nElements)
{
   RPage newPage(columnId, mem, elementSize * nElements, elementSize);
   newPage.TryGrow(nElements);
   return newPage;
}

void ROOT::Experimental::Detail::RPageAllocatorKey::DeletePage(
   const RPage& page, ROOT::Experimental::Internal::RNTupleBlob *payload)
{
   if (page.IsNull())
      return;
   R__ASSERT(page.GetBuffer() == payload->fContent);
   free(payload->fContent);
   delete payload;
}


////////////////////////////////////////////////////////////////////////////////


ROOT::Experimental::Detail::RPageSourceRoot::RPageSourceRoot(std::string_view ntupleName, RSettings settings)
   : RPageSource(ntupleName)
   , fPageAllocator(std::make_unique<RPageAllocatorKey>())
   , fPagePool(std::make_shared<RPagePool>())
   , fDirectory(nullptr)
   , fSettings(settings)
{
}

ROOT::Experimental::Detail::RPageSourceRoot::RPageSourceRoot(std::string_view ntupleName, std::string_view path)
   : RPageSource(ntupleName)
   , fPageAllocator(std::make_unique<RPageAllocatorKey>())
   , fPagePool(std::make_shared<RPagePool>())
   , fDirectory(nullptr)
{
   TFile *file = TFile::Open(std::string(path).c_str(), "READ");
   fSettings.fFile = file;
   fSettings.fTakeOwnership = true;
}


ROOT::Experimental::Detail::RPageSourceRoot::~RPageSourceRoot()
{
   if (fSettings.fTakeOwnership) {
      fSettings.fFile->Close();
      delete fSettings.fFile;
   }
}


ROOT::Experimental::Detail::RPageStorage::ColumnHandle_t
ROOT::Experimental::Detail::RPageSourceRoot::AddColumn(DescriptorId_t fieldId, const RColumn &column)
{
   R__ASSERT(fieldId != kInvalidDescriptorId);
   auto& model = column.GetModel();
   auto columnId = fDescriptor.FindColumnId(fieldId, column.GetIndex());
   R__ASSERT(columnId != kInvalidDescriptorId);
   R__ASSERT(model == *fMapper.fId2ColumnModel[columnId]);
   //printf("Attaching column %s id %d type %d length %lu\n",
   //   column->GetModel().GetName().c_str(), columnId, (int)(column->GetModel().GetType()),
   //   fMapper.fColumnIndex[columnId].fNElements);
   return ColumnHandle_t(columnId, &column);
}


void ROOT::Experimental::Detail::RPageSourceRoot::Attach()
{
   fDirectory = fSettings.fFile->GetDirectory(fNTupleName.c_str());
   RNTupleDescriptorBuilder descBuilder;

   auto keyRawNTupleHeader = fDirectory->GetKey(kKeyNTupleHeader);
   auto ntupleRawHeader = keyRawNTupleHeader->ReadObject<ROOT::Experimental::Internal::RNTupleBlob>();
   descBuilder.SetFromHeader(ntupleRawHeader->fContent);
   free(ntupleRawHeader->fContent);
   delete ntupleRawHeader;

   auto keyRawNTupleFooter = fDirectory->GetKey(kKeyNTupleFooter);
   auto ntupleRawFooter = keyRawNTupleFooter->ReadObject<ROOT::Experimental::Internal::RNTupleBlob>();
   descBuilder.AddClustersFromFooter(ntupleRawFooter->fContent);
   free(ntupleRawFooter->fContent);
   delete ntupleRawFooter;

   fDescriptor = descBuilder.GetDescriptor();

   auto keyNTupleHeader = fDirectory->GetKey(RMapper::kKeyNTupleHeader);
   auto ntupleHeader = keyNTupleHeader->ReadObject<ROOT::Experimental::Internal::RNTupleHeader>();
   //printf("Number of fields %lu, of columns %lu\n", ntupleHeader->fFields.size(), ntupleHeader->fColumns.size());

   for (auto &fieldHeader : ntupleHeader->fFields) {
      if (fieldHeader.fParentName.empty()) {
         fMapper.fRootFields.push_back(RMapper::RFieldDescriptor(fieldHeader.fName, fieldHeader.fType));
      }
   }

   auto nColumns = ntupleHeader->fColumns.size();
   //fPagePool = std::make_unique<RPagePool>();
   fMapper.fColumnIndex.resize(nColumns);

   std::int32_t columnId = 0;
   for (auto &columnHeader : ntupleHeader->fColumns) {
      auto columnModel = std::make_unique<RColumnModel>(
         columnHeader.fName, columnHeader.fType, columnHeader.fIsSorted);
      fMapper.fId2ColumnModel[columnId] = std::move(columnModel);
      fMapper.fColumnName2Id[columnHeader.fName] = columnId;
      columnId++;
   }

   /// Determine column dependencies (offset - pointee relationships)
   for (auto &columnHeader : ntupleHeader->fColumns) {
      if (columnHeader.fOffsetColumn.empty()) continue;
      fMapper.fColumn2Pointee[fMapper.fColumnName2Id[columnHeader.fOffsetColumn]] =
        fMapper.fColumnName2Id[columnHeader.fName];
   }

   auto keyNTupleFooter = fDirectory->GetKey(RMapper::kKeyNTupleFooter);
   auto ntupleFooter = keyNTupleFooter->ReadObject<ROOT::Experimental::Internal::RNTupleFooter>();
   //printf("Number of clusters: %d, entries %ld\n", ntupleFooter->fNClusters, ntupleFooter->fNEntries);

   for (std::int32_t iCluster = 0; iCluster < ntupleFooter->fNClusters; ++iCluster) {
      auto keyClusterFooter = fDirectory->GetKey((RMapper::kKeyClusterFooter + std::to_string(iCluster)).c_str());
      auto clusterFooter = keyClusterFooter->ReadObject<ROOT::Experimental::Internal::RClusterFooter>();
      R__ASSERT(clusterFooter->fPagesPerColumn.size() == nColumns);
      for (unsigned iColumn = 0; iColumn < nColumns; ++iColumn) {
         if (clusterFooter->fPagesPerColumn[iColumn].fRangeStarts.empty())
            continue;
         NTupleSize_t selfClusterOffset = clusterFooter->fPagesPerColumn[iColumn].fRangeStarts[0];
         NTupleSize_t pointeeClusterOffset = kInvalidNTupleIndex;
         auto itrPointee = fMapper.fColumn2Pointee.find(iColumn);
         if (itrPointee != fMapper.fColumn2Pointee.end()) {
            //printf("COLUMN %s wants to know pointee offset of column %s\n",
            //  fMapper.fId2ColumnModel[iColumn]->GetName().c_str(),
            //  fMapper.fId2ColumnModel[itrPointee->second]->GetName().c_str());
            /// The pointee might not have any pages in this cluster (e.g. all empty collections)
            if (!clusterFooter->fPagesPerColumn[itrPointee->second].fRangeStarts.empty())
               pointeeClusterOffset = clusterFooter->fPagesPerColumn[itrPointee->second].fRangeStarts[0];
         }
         NTupleSize_t pageInCluster = 0;
         for (auto rangeStart : clusterFooter->fPagesPerColumn[iColumn].fRangeStarts) {
            fMapper.fColumnIndex[iColumn].fRangeStarts.push_back(rangeStart);
            fMapper.fColumnIndex[iColumn].fClusterId.push_back(iCluster);
            fMapper.fColumnIndex[iColumn].fPageInCluster.push_back(pageInCluster);
            fMapper.fColumnIndex[iColumn].fSelfClusterOffset.push_back(selfClusterOffset);
            fMapper.fColumnIndex[iColumn].fPointeeClusterOffset.push_back(pointeeClusterOffset);
            pageInCluster++;
         }
      }
      delete clusterFooter;
   }

   for (unsigned iColumn = 0; iColumn < nColumns; ++iColumn) {
      fMapper.fColumnIndex[iColumn].fNElements = ntupleFooter->fNElementsPerColumn[iColumn];
   }
   fMapper.fNEntries = ntupleFooter->fNEntries;

   delete ntupleFooter;
   delete ntupleHeader;
}


ROOT::Experimental::Detail::RPage ROOT::Experimental::Detail::RPageSourceRoot::PopulatePage(
   ColumnHandle_t columnHandle, NTupleSize_t index)
{
   auto columnId = columnHandle.fId;
   auto cachedPage = fPagePool->GetPage(columnId, index);
   if (!cachedPage.IsNull())
      return cachedPage;

   auto columnDescriptor = fDescriptor.GetColumnDescriptor(columnId);

   // TODO(jblomer): save last used cluster id and check it first
   auto clusterId = fDescriptor.FindClusterId(columnId, index);
   R__ASSERT(clusterId != kInvalidDescriptorId);
   auto clusterDescriptor = fDescriptor.GetClusterDescriptor(clusterId);
   auto pageRange = clusterDescriptor.GetPageRange(columnId);
   auto selfOffset = clusterDescriptor.GetColumnRange(columnId).fFirstElementIndex;

   // TODO(jblomer): binary search
   RClusterDescriptor::RPageRange::RPageInfo pageInfo;
   auto firstInPage = selfOffset;
   for (const auto &pi : pageRange.fPageInfos) {
      if (firstInPage + pi.fNElements > index) {
         pageInfo = pi;
         break;
      }
      firstInPage += pi.fNElements;
   }
   R__ASSERT(firstInPage <= index);
   R__ASSERT((firstInPage + pageInfo.fNElements) > index);

   NTupleSize_t pointeeOffset = 0;
   // TODO(jblomer): deal with multiple linked columns
   if (!columnDescriptor.GetLinkIds().empty())
      pointeeOffset = clusterDescriptor.GetColumnRange(columnDescriptor.GetLinkIds()[0]).fFirstElementIndex;

   //printf("Populating page %lu/%lu [%lu] for column %d starting at %lu\n", clusterId, pageInCluster, pageIdx, columnId, firstInPage);

   std::string keyName = std::string(RMapper::kKeyPagePayload) +
      std::to_string(clusterId) + RMapper::kKeySeparator +
      std::to_string(pageInfo.fLocator);
   auto pageKey = fDirectory->GetKey(keyName.c_str());
   auto pagePayload = pageKey->ReadObject<ROOT::Experimental::Internal::RNTupleBlob>();
   auto elementSize = pagePayload->fSize / pageInfo.fNElements;
   R__ASSERT(pagePayload->fSize % pageInfo.fNElements == 0);
   auto newPage = fPageAllocator->NewPage(columnId, pagePayload->fContent, elementSize, pageInfo.fNElements);
   newPage.SetWindow(firstInPage, RPage::RClusterInfo(clusterId, selfOffset, pointeeOffset));
   fPagePool->RegisterPage(newPage,
      RPageDeleter([](const RPage &page, void *userData)
      {
         RPageAllocatorKey::DeletePage(page, reinterpret_cast<ROOT::Experimental::Internal::RNTupleBlob *>(userData));
      }, pagePayload));
   return newPage;
}

void ROOT::Experimental::Detail::RPageSourceRoot::ReleasePage(RPage &page)
{
   fPagePool->ReturnPage(page);
}

ROOT::Experimental::NTupleSize_t ROOT::Experimental::Detail::RPageSourceRoot::GetNEntries()
{
   return fMapper.fNEntries;
}

ROOT::Experimental::NTupleSize_t ROOT::Experimental::Detail::RPageSourceRoot::GetNElements(ColumnHandle_t columnHandle)
{
   return fMapper.fColumnIndex[columnHandle.fId].fNElements;
}

ROOT::Experimental::ColumnId_t ROOT::Experimental::Detail::RPageSourceRoot::GetColumnId(ColumnHandle_t columnHandle)
{
   // TODO(jblomer) distinguish trees
   return columnHandle.fId;
}
