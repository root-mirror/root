/// \file RPageSourceFriends.cxx
/// \ingroup NTuple ROOT7
/// \author Jakob Blomer <jblomer@cern.ch>
/// \date 2019-08-10
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/RCluster.hxx>
#include <ROOT/RError.hxx>
#include <ROOT/RLogger.hxx>
#include <ROOT/RNTupleOptions.hxx>
#include <ROOT/RPageSourceFriends.hxx>

#include <utility>

ROOT::Experimental::Detail::RPageSourceFriends::RPageSourceFriends(
   std::string_view ntupleName, std::span<std::unique_ptr<RPageSource>> sources)
   : RPageSource(ntupleName, RNTupleReadOptions())
   , fMetrics(std::string(ntupleName))
{
   for (auto &s : sources) {
      fSources.emplace_back(std::move(s));
      fMetrics.ObserveMetrics(fSources.back()->GetMetrics());
   }
}

ROOT::Experimental::Detail::RPageSourceFriends::~RPageSourceFriends() = default;


void ROOT::Experimental::Detail::RPageSourceFriends::AddVirtualField(
   const RNTupleDescriptor &originDesc,
   std::size_t originIdx,
   const RFieldDescriptor &originField,
   DescriptorId_t virtualParent,
   const std::string &virtualName)
{
   auto virtualFieldId = fNextId++;
   auto virtualField = RDanglingFieldDescriptor(originField)
      .FieldId(virtualFieldId)
      .FieldName(virtualName)
      .MakeDescriptor().Unwrap();
   fBuilder.AddField(virtualField);
   fBuilder.AddFieldLink(virtualParent, virtualFieldId);

   for (const auto &f : originDesc.GetFieldRange(originField))
      AddVirtualField(originDesc, originIdx, f, virtualFieldId, f.GetFieldName());

   for (const auto &c: originDesc.GetColumnRange(originField)) {
      fBuilder.AddColumn(fNextId, virtualFieldId, c.GetVersion(), c.GetModel(), c.GetIndex());
      fIdBiMap.Insert({originIdx, c.GetId()}, fNextId);
      fNextId++;
   }
}


ROOT::Experimental::RNTupleDescriptor ROOT::Experimental::Detail::RPageSourceFriends::AttachImpl()
{
   fBuilder.SetNTuple(fNTupleName, "", "", RNTupleVersion(), RNTupleUuid());
   fBuilder.AddField(RDanglingFieldDescriptor()
      .FieldId(0)
      .Structure(ENTupleStructure::kRecord)
      .MakeDescriptor()
      .Unwrap());

   for (std::size_t i = 0; i < fSources.size(); ++i) {
      fSources[i]->Attach();
      if (fSources[i]->GetNEntries() != fSources[0]->GetNEntries()) {
         fNextId = 1;
         fIdBiMap.Clear();
         fBuilder.Reset();
         throw RException(R__FAIL("mismatch in the number of entries of friend RNTuples"));
      }
      const auto &desc = fSources[i]->GetDescriptor();
      AddVirtualField(desc, i, desc.GetFieldZero(), 0, desc.GetName());

      for (const auto &c : desc.GetClusterRange()) {
         fBuilder.AddCluster(fNextId, c.GetVersion(), c.GetFirstEntryIndex(), c.GetNEntries());
         for (auto originColumnId : c.GetColumnIds()) {
            DescriptorId_t virtualColumnId = fIdBiMap.GetVirtualId({i, originColumnId});
            fBuilder.AddClusterPageRange(virtualColumnId, c.GetPageRange(originColumnId).Clone());
            fBuilder.AddClusterColumnRange(virtualColumnId, c.GetColumnRange(originColumnId));
         }
         fIdBiMap.Insert({i, c.GetId()}, fNextId);
         fNextId++;
      }
   }

   fBuilder.EnsureValidDescriptor();
   return fBuilder.MoveDescriptor();
}


std::unique_ptr<ROOT::Experimental::Detail::RPageSource>
ROOT::Experimental::Detail::RPageSourceFriends::Clone() const
{
   // TODO
   return nullptr;
}


ROOT::Experimental::Detail::RPage
ROOT::Experimental::Detail::RPageSourceFriends::PopulatePage(
   ColumnHandle_t columnHandle, NTupleSize_t globalIndex)
{
   auto virtualColumnId = columnHandle.fId;
   auto originColumnId = fIdBiMap.GetOriginId(virtualColumnId);
   columnHandle.fId = originColumnId.fId;

   auto page = fSources[originColumnId.fSourceIdx]->PopulatePage(columnHandle, globalIndex);

   fPage2SourceIdx[page.GetBuffer()] = originColumnId.fSourceIdx;
   auto virtualClusterId = fIdBiMap.GetVirtualId({originColumnId.fSourceIdx, page.GetClusterInfo().GetId()});
   page.ChangeIds(virtualColumnId, virtualClusterId);

   return page;
}


ROOT::Experimental::Detail::RPage
ROOT::Experimental::Detail::RPageSourceFriends::PopulatePage(
   ColumnHandle_t columnHandle, const RClusterIndex &clusterIndex)
{
   auto virtualColumnId = columnHandle.fId;
   auto originColumnId = fIdBiMap.GetOriginId(virtualColumnId);
   RClusterIndex originClusterIndex(
      fIdBiMap.GetOriginId(clusterIndex.GetClusterId()).fId,
      clusterIndex.GetIndex());
   columnHandle.fId = originColumnId.fId;

   auto page = fSources[originColumnId.fSourceIdx]->PopulatePage(columnHandle, originClusterIndex);

   fPage2SourceIdx[page.GetBuffer()] = originColumnId.fSourceIdx;
   page.ChangeIds(virtualColumnId, clusterIndex.GetClusterId());

   return RPage();
}

void ROOT::Experimental::Detail::RPageSourceFriends::ReleasePage(RPage &page)
{
   fSources[fPage2SourceIdx[page.GetBuffer()]]->ReleasePage(page);
}


std::unique_ptr<ROOT::Experimental::Detail::RCluster>
ROOT::Experimental::Detail::RPageSourceFriends::LoadCluster(DescriptorId_t /* clusterId */,
                                                            const ColumnSet_t & /* columns */)
{
   // The virtual friends page source does not pre-load any clusters itself. However, the underlying page sources
   // that are combined may well do it.
   return nullptr;
}
