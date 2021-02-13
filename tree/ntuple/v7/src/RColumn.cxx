/// \file RColumn.cxx
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

#include <ROOT/RColumn.hxx>
#include <ROOT/RColumnModel.hxx>
#include <ROOT/RPageStorage.hxx>

#include <TError.h>

#include <iostream>

ROOT::Experimental::Detail::RColumn::RColumn(const RColumnModel& model, std::uint32_t index)
   : fModel(model), fIndex(index), fPageSink(nullptr), fPageSource(nullptr), fHeadPage(), fNElements(0),
     fCurrentPage(),
     fColumnIdSource(kInvalidColumnId)
{
}

ROOT::Experimental::Detail::RColumn::~RColumn()
{
   if (!fHeadPage.IsNull())
      fPageSink->ReleasePage(fHeadPage);
   if (!fCurrentPage.IsNull())
      fPageSource->ReleasePage(fCurrentPage);
   if (fHandleSink)
      fPageSink->DropColumn(fHandleSink);
   if (fHandleSource)
      fPageSource->DropColumn(fHandleSource);
}

void ROOT::Experimental::Detail::RColumn::Connect(DescriptorId_t fieldId, RPageStorage *pageStorage)
{
   switch (pageStorage->GetType()) {
   case EPageStorageType::kSink:
      fPageSink = static_cast<RPageSink*>(pageStorage); // the page sink initializes fHeadPage on AddColumn
      fHandleSink = fPageSink->AddColumn(fieldId, *this);
      fHeadPage = fPageSink->ReservePage(fHandleSink);
      break;
   case EPageStorageType::kSource:
      fPageSource = static_cast<RPageSource*>(pageStorage);
      fHandleSource = fPageSource->AddColumn(fieldId, *this);
      fNElements = fPageSource->GetNElements(fHandleSource);
      fColumnIdSource = fPageSource->GetColumnId(fHandleSource);
      break;
   default:
      R__ASSERT(false);
   }
}

void ROOT::Experimental::Detail::RColumn::Flush()
{
   if (fHeadPage.GetSize() == 0) return;

   fPageSink->CommitPage(fHandleSink, fHeadPage);
   fHeadPage.Reset(fNElements);
}

void ROOT::Experimental::Detail::RColumn::MapPage(const NTupleSize_t index)
{
   fPageSource->ReleasePage(fCurrentPage);
   fCurrentPage = fPageSource->PopulatePage(fHandleSource, index);
}

void ROOT::Experimental::Detail::RColumn::MapPage(const RClusterIndex &clusterIndex)
{
   fPageSource->ReleasePage(fCurrentPage);
   fCurrentPage = fPageSource->PopulatePage(fHandleSource, clusterIndex);
}
