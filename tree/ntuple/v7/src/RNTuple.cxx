/// \file RNTuple.cxx
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

#include <ROOT/RNTuple.hxx>

#include <ROOT/RFieldVisitor.hxx>
#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RPageStorage.hxx>
#include "ROOT/RPageStorageFile.hxx"
#ifdef R__USE_IMT
#include <ROOT/TTaskGroup.hxx>
#endif

#include <TError.h>
#include <TROOT.h> // for IsImplicitMTEnabled()

#include <algorithm>
#include <exception>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>


void ROOT::Experimental::RNTupleImtTaskScheduler::Reset()
{
   fTaskGroup = std::make_unique<TTaskGroup>();
}


void ROOT::Experimental::RNTupleImtTaskScheduler::AddTask(const std::function<void(void)> &taskFunc)
{
   fTaskGroup->Run(taskFunc);
}


void ROOT::Experimental::RNTupleImtTaskScheduler::Wait()
{
   fTaskGroup->Wait();
}


//------------------------------------------------------------------------------


void ROOT::Experimental::RNTupleReader::ConnectModel(const RNTupleModel &model) {
   std::unordered_map<const Detail::RFieldBase *, DescriptorId_t> fieldPtr2Id;
   fieldPtr2Id[model.GetFieldZero()] = fSource->GetDescriptor().GetFieldZeroId();
   for (auto &field : *model.GetFieldZero()) {
      auto parentId = fieldPtr2Id[field.GetParent()];
      auto fieldId = fSource->GetDescriptor().FindFieldId(field.GetName(), parentId);
      R__ASSERT(fieldId != kInvalidDescriptorId);
      fieldPtr2Id[&field] = fieldId;
      Detail::RFieldFuse::Connect(fieldId, *fSource, field);
   }
}

void ROOT::Experimental::RNTupleReader::InitPageSource()
{
#ifdef R__USE_IMT
   if (IsImplicitMTEnabled()) {
      fSource->SetTaskScheduler(&fUnzipTasks);
   }
#endif
   fSource->Attach();
   fMetrics.ObserveMetrics(fSource->GetMetrics());
}

ROOT::Experimental::RNTupleReader::RNTupleReader(
   std::unique_ptr<ROOT::Experimental::RNTupleModel> model,
   std::unique_ptr<ROOT::Experimental::Detail::RPageSource> source)
   : fSource(std::move(source))
   , fModel(std::move(model))
   , fMetrics("RNTupleReader")
{
   InitPageSource();
   ConnectModel(*fModel);
}

ROOT::Experimental::RNTupleReader::RNTupleReader(std::unique_ptr<ROOT::Experimental::Detail::RPageSource> source)
   : fSource(std::move(source))
   , fModel(nullptr)
   , fMetrics("RNTupleReader")
{
   InitPageSource();
}

ROOT::Experimental::RNTupleReader::~RNTupleReader()
{
#ifdef R__USE_IMT
   fSource->SetTaskScheduler(nullptr);
#endif
}

std::unique_ptr<ROOT::Experimental::RNTupleReader> ROOT::Experimental::RNTupleReader::Open(
   std::unique_ptr<RNTupleModel> model,
   std::string_view ntupleName,
   std::string_view storage,
   const RNTupleReadOptions &options)
{
   return std::make_unique<RNTupleReader>(std::move(model), Detail::RPageSource::Create(ntupleName, storage, options));
}

std::unique_ptr<ROOT::Experimental::RNTupleReader> ROOT::Experimental::RNTupleReader::Open(
   std::string_view ntupleName,
   std::string_view storage,
   const RNTupleReadOptions &options)
{
   return std::make_unique<RNTupleReader>(Detail::RPageSource::Create(ntupleName, storage, options));
}

ROOT::Experimental::RNTupleModel *ROOT::Experimental::RNTupleReader::GetModel()
{
   if (!fModel) {
      fModel = fSource->GetDescriptor().GenerateModel();
      ConnectModel(*fModel);
   }
   return fModel.get();
}

void ROOT::Experimental::RNTupleReader::PrintInfo(const ENTupleInfo what, std::ostream &output)
{
   // TODO(lesimon): In a later version, these variables may be defined by the user or the ideal width may be read out from the terminal.
   char frameSymbol = '*';
   int width = 80;
   /*
   if (width < 30) {
      output << "The width is too small! Should be at least 30." << std::endl;
      return;
   }
   */
   std::string name = fSource->GetDescriptor().GetName();
   switch (what) {
   case ENTupleInfo::kSummary: {
      for (int i = 0; i < (width/2 + width%2 - 4); ++i)
            output << frameSymbol;
      output << " NTUPLE ";
      for (int i = 0; i < (width/2 - 4); ++i)
         output << frameSymbol;
      output << std::endl;
      // FitString defined in RFieldVisitor.cxx
      output << frameSymbol << " N-Tuple : " << RNTupleFormatter::FitString(name, width-13) << frameSymbol << std::endl; // prints line with name of ntuple
      output << frameSymbol << " Entries : " << RNTupleFormatter::FitString(std::to_string(GetNEntries()), width - 13) << frameSymbol << std::endl;  // prints line with number of entries

      // Traverses through all fields to gather information needed for printing.
      RPrepareVisitor prepVisitor;
      // Traverses through all fields to do the actual printing.
      RPrintSchemaVisitor printVisitor(output);

      // Note that we do not need to connect the model, we are only looking at its tree of fields
      auto fullModel = fSource->GetDescriptor().GenerateModel();
      fullModel->GetFieldZero()->AcceptVisitor(prepVisitor);

      printVisitor.SetFrameSymbol(frameSymbol);
      printVisitor.SetWidth(width);
      printVisitor.SetDeepestLevel(prepVisitor.GetDeepestLevel());
      printVisitor.SetNumFields(prepVisitor.GetNumFields());

      for (int i = 0; i < width; ++i)
         output << frameSymbol;
      output << std::endl;
      fullModel->GetFieldZero()->AcceptVisitor(printVisitor);
      for (int i = 0; i < width; ++i)
         output << frameSymbol;
      output << std::endl;
      break;
   }
   case ENTupleInfo::kStorageDetails:
      fSource->GetDescriptor().PrintInfo(output);
      break;
   case ENTupleInfo::kMetrics:
      fMetrics.Print(output);
      break;
   default:
      // Unhandled case, internal error
      R__ASSERT(false);
   }
}


ROOT::Experimental::RNTupleReader *ROOT::Experimental::RNTupleReader::GetDisplayReader()
{
   if (!fDisplayReader)
      fDisplayReader = Clone();
   return fDisplayReader.get();
}


void ROOT::Experimental::RNTupleReader::Show(NTupleSize_t index, const ENTupleShowFormat format, std::ostream &output)
{
   RNTupleReader *reader = this;
   REntry *entry = nullptr;
   // Don't accidentally trigger loading of the entire model
   if (fModel)
      entry = fModel->GetDefaultEntry();

   switch(format) {
   case ENTupleShowFormat::kCompleteJSON:
      reader = GetDisplayReader();
      entry = reader->GetModel()->GetDefaultEntry();
      // Fall through
   case ENTupleShowFormat::kCurrentModelJSON:
      if (!entry) {
         output << "{}" << std::endl;
         break;
      }

      reader->LoadEntry(index);
      output << "{";
      for (auto iValue = entry->begin(); iValue != entry->end(); ) {
         output << std::endl;
         RPrintValueVisitor visitor(*iValue, output, 1 /* level */);
         iValue->GetField()->AcceptVisitor(visitor);

         if (++iValue == entry->end()) {
            output << std::endl;
            break;
         } else {
            output << ",";
         }
      }
      output << "}" << std::endl;
      break;
   default:
      // Unhandled case, internal error
      R__ASSERT(false);
   }
}


//------------------------------------------------------------------------------


ROOT::Experimental::RNTupleWriter::RNTupleWriter(
   std::unique_ptr<ROOT::Experimental::RNTupleModel> model,
   std::unique_ptr<ROOT::Experimental::Detail::RPageSink> sink)
   : fSink(std::move(sink))
   , fModel(std::move(model))
   , fClusterSizeEntries(kDefaultClusterSizeEntries)
   , fLastCommitted(0)
   , fNEntries(0)
{
   fSink->Create(*fModel.get());
}

ROOT::Experimental::RNTupleWriter::~RNTupleWriter()
{
   CommitCluster();
   fSink->CommitDataset();
}

std::unique_ptr<ROOT::Experimental::RNTupleWriter> ROOT::Experimental::RNTupleWriter::Recreate(
   std::unique_ptr<RNTupleModel> model,
   std::string_view ntupleName,
   std::string_view storage,
   const RNTupleWriteOptions &options)
{
   return std::make_unique<RNTupleWriter>(std::move(model), Detail::RPageSink::Create(ntupleName, storage, options));
}

std::unique_ptr<ROOT::Experimental::RNTupleWriter> ROOT::Experimental::RNTupleWriter::Append(
   std::unique_ptr<RNTupleModel> model,
   std::string_view ntupleName,
   TFile &file,
   const RNTupleWriteOptions &options)
{
   auto sink = std::make_unique<Detail::RPageSinkFile>(ntupleName, file, options);
   return std::make_unique<RNTupleWriter>(std::move(model), std::move(sink));
}


void ROOT::Experimental::RNTupleWriter::CommitCluster()
{
   if (fNEntries == fLastCommitted) return;
   for (auto& field : *fModel->GetFieldZero()) {
      field.Flush();
      field.CommitCluster();
   }
   fSink->CommitCluster(fNEntries);
   fLastCommitted = fNEntries;
}


//------------------------------------------------------------------------------


ROOT::Experimental::RCollectionNTuple::RCollectionNTuple(std::unique_ptr<REntry> defaultEntry)
   : fOffset(0), fDefaultEntry(std::move(defaultEntry))
{
}
