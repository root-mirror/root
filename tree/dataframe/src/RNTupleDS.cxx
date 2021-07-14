/// \file RNTupleDS.cxx
/// \ingroup NTuple ROOT7
/// \author Jakob Blomer <jblomer@cern.ch>
/// \author Enrico Guiraud <enrico.guiraud@cern.ch>
/// \date 2018-10-04
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/RDF/RColumnReaderBase.hxx>
#include <ROOT/RField.hxx>
#include <ROOT/RFieldValue.hxx>
#include <ROOT/RNTupleDescriptor.hxx>
#include <ROOT/RNTupleDS.hxx>
#include <ROOT/RNTupleUtil.hxx>
#include <ROOT/RPageStorage.hxx>
#include <ROOT/RStringView.hxx>

#include <TError.h>

#include <string>
#include <vector>
#include <typeinfo>
#include <utility>

namespace ROOT {
namespace Experimental {
namespace Internal {

/// An artificial field that transforms an RNTuple column that contains the offset of collections into
/// collection sizes. It is used to provide the "number of" RDF columns for collections, e.g.
/// `__rdf_sizeof_jets` for a collection named `jets`.
///
/// This field owns the collection offset field but instead of exposing the collection offsets it exposes
/// the collection sizes (offset(N+1) - offset(N)).  For the time being, we offer this functionality only in RDataFrame.
/// TODO(jblomer): consider providing a general set of useful virtual fields as part of RNTuple.
class RRDFCardinalityField : public ROOT::Experimental::Detail::RFieldBase {
protected:
   std::unique_ptr<ROOT::Experimental::Detail::RFieldBase> CloneImpl(std::string_view /* newName */) const final
   {
      return std::make_unique<RRDFCardinalityField>();
   }

public:
   static std::string TypeName() { return "std::size_t"; }
   RRDFCardinalityField()
      : ROOT::Experimental::Detail::RFieldBase("", TypeName(), ENTupleStructure::kLeaf, false /* isSimple */) {}
   RRDFCardinalityField(RRDFCardinalityField &&other) = default;
   RRDFCardinalityField &operator=(RRDFCardinalityField &&other) = default;
   ~RRDFCardinalityField() = default;

   // Field is only used for reading
   void GenerateColumnsImpl() final { R__DEBUG_ASSERT(false && "Cardinality fields must only be used for reading"); }

   void GenerateColumnsImpl(const RNTupleDescriptor &) final
   {
      RColumnModel model(EColumnType::kIndex, true /* isSorted*/);
      fColumns.emplace_back(std::unique_ptr<ROOT::Experimental::Detail::RColumn>(
         ROOT::Experimental::Detail::RColumn::Create<ClusterSize_t, EColumnType::kIndex>(model, 0)));
      fPrincipalColumn = fColumns[0].get();
   }

   ROOT::Experimental::Detail::RFieldValue GenerateValue(void *where) final
   {
      return ROOT::Experimental::Detail::RFieldValue(this, static_cast<std::size_t *>(where));
   }
   ROOT::Experimental::Detail::RFieldValue CaptureValue(void *where) final
   {
      return ROOT::Experimental::Detail::RFieldValue(true /* captureFlag */, this, where);
   }
   size_t GetValueSize() const final { return sizeof(std::size_t); }

   /// Get the number of elements of the collection identified by globalIndex
   void
   ReadGlobalImpl(ROOT::Experimental::NTupleSize_t globalIndex, ROOT::Experimental::Detail::RFieldValue *value) final
   {
      RClusterIndex collectionStart;
      ClusterSize_t size;
      fPrincipalColumn->GetCollectionInfo(globalIndex, &collectionStart, &size);
      *value->Get<std::size_t>() = size;
   }

   /// Get the number of elements of the collection identified by clusterIndex
   void ReadInClusterImpl(const ROOT::Experimental::RClusterIndex &clusterIndex,
                          ROOT::Experimental::Detail::RFieldValue *value) final
   {
      RClusterIndex collectionStart;
      ClusterSize_t size;
      fPrincipalColumn->GetCollectionInfo(clusterIndex, &collectionStart, &size);
      *value->Get<std::size_t>() = size;
   }
};

/// Every RDF column is represented by exactly one RNTuple field
class RNTupleColumnReader : public ROOT::Detail::RDF::RColumnReaderBase {
   using RFieldBase = ROOT::Experimental::Detail::RFieldBase;
   using RFieldValue = ROOT::Experimental::Detail::RFieldValue;
   using RPageSource = ROOT::Experimental::Detail::RPageSource;

   std::unique_ptr<RFieldBase> fField; ///< The field backing the RDF column
   RFieldValue fValue;                 ///< The memory location used to read from fField
   Long64_t fLastEntry;                ///< Last entry number that was read

public:
   RNTupleColumnReader(std::unique_ptr<RFieldBase> f)
      : fField(std::move(f)), fValue(fField->GenerateValue()), fLastEntry(-1)
   {
   }
   virtual ~RNTupleColumnReader() { fField->DestroyValue(fValue); }

   /// Column readers are created as prototype and then cloned for every slot
   std::unique_ptr<RNTupleColumnReader> Clone()
   {
      return std::make_unique<RNTupleColumnReader>(fField->Clone(fField->GetName()));
   }

   /// Connect the field and its subfields to the page source
   void Connect(RPageSource &source)
   {
      fField->ConnectPageSource(source);
      for (auto &f : *fField)
         f.ConnectPageSource(source);
   }

   void *GetImpl(Long64_t entry) final
   {
      if (entry != fLastEntry) {
         fField->Read(entry, &fValue);
         fLastEntry = entry;
      }
      return fValue.GetRawPtr();
   }
};

} // namespace Internal

RNTupleDS::~RNTupleDS() = default;

void RNTupleDS::AddField(const RNTupleDescriptor &desc, std::string_view colName, DescriptorId_t fieldId,
                         std::vector<DescriptorId_t> skeinIDs)
{
   // As an example for the mapping of RNTuple fields to RDF columns, let's consider an RNTuple
   // using the following types and with a top-level field named "event" of type Event:
   //
   // struct Event {
   //    int id;
   //    std::vector<Track> tracks;
   // };
   // struct Track {
   //    std::vector<Hit> hits;
   // };
   // struct Hit {
   //    float x;
   //    float y;
   // };
   //
   // AddField() will be called from the constructor with the RNTuple root field (ENTupleStructure::kRecord).
   // From there, we recurse into the "event" sub field (also ENTupleStructure::kRecord) and further down the
   // tree of sub fields and expose the following RDF columns:
   // TODO(jblomer): Collections should be exposed as RVec<T> instead of std::vector<T>
   //
   // "event"                             [Event]
   // "event.id"                          [int]
   // "event.tracks"                      [std::vector<Track>]
   // "__rdf_sizeof_event.tracks"         [unsigned int]
   // "event.tracks.hits"                 [std::vector<std::vector<Hit>>]
   // "__rdf_sizeof_event.tracks.hits"    [std::vector<unsigned int>]
   // "event.tracks.hits.x"               [std::vector<std::vector<float>>]
   // "__rdf_sizeof_event.tracks.hits.x"  [std::vector<unsigned int>]
   // "event.tracks.hits.y"               [std::vector<std::vector<float>>]
   // "__rdf_sizeof_event.tracks.hits.y"  [std::vector<unsigned int>]

   const auto &fieldDesc = desc.GetFieldDescriptor(fieldId);
   if (fieldDesc.GetStructure() == ENTupleStructure::kCollection) {
      // Inner fields of collections are provided as projected collections of only that inner field,
      // E.g. we provide a projected collection vector<vector<float>> for "event.tracks.hits.x" in the example
      // above.

      // We open a new collection scope with fieldID being the inner most collection. E.g. for "event.tracks.hits",
      // skeinIDs would already contain the fieldID of "event.tracks"
      skeinIDs.emplace_back(fieldId);
      // There should only be one sub field but it's easiest to access via the sub field range
      for (const auto &f : desc.GetFieldIterable(fieldDesc.GetId())) {
         AddField(desc, colName, f.GetId(), skeinIDs);
      }
      // Note that at the end of the recursion, we handled the inner sub collections as well as the
      // collection as whole, so we are done.
      return;
   } else if (fieldDesc.GetStructure() == ENTupleStructure::kRecord) {
      // Inner fields of records are provided as individual RDF columns, e.g. "event.id"
      for (const auto &f : desc.GetFieldIterable(fieldDesc.GetId())) {
         auto innerName = colName.empty() ? f.GetFieldName() : (std::string(colName) + "." + f.GetFieldName());
         AddField(desc, innerName, f.GetId(), skeinIDs);
      }
   }

   // The fieldID could be the root field or the class of fieldId might not be loaded.
   // In these cases, only the inner fields are exposed as RDF columns.
   auto fieldOrException = Detail::RFieldBase::Create("", fieldDesc.GetTypeName());
   if (!fieldOrException)
      return;
   auto valueField = fieldOrException.Unwrap();
   valueField->SetOnDiskId(fieldId);
   std::unique_ptr<Detail::RFieldBase> cardinalityField;
   // Collections get the additional "number of" RDF column (e.g. "__rdf_sizeof_tracks")
   if (!skeinIDs.empty()) {
      cardinalityField = std::make_unique<ROOT::Experimental::Internal::RRDFCardinalityField>();
      cardinalityField->SetOnDiskId(skeinIDs.back());
   }

   std::string typeName;
   for (auto i = skeinIDs.rbegin(); i != skeinIDs.rend(); ++i) {
      valueField = std::make_unique<ROOT::Experimental::RVectorField>("", std::move(valueField));
      valueField->SetOnDiskId(*i);
      // Skip the inner-most collection level to construct the cardinality column
      if (i != skeinIDs.rbegin()) {
         cardinalityField = std::make_unique<ROOT::Experimental::RVectorField>("", std::move(cardinalityField));
         cardinalityField->SetOnDiskId(*i);
      }
   }

   if (cardinalityField) {
      fColumnNames.emplace_back("__rdf_sizeof_" + std::string(colName));
      fColumnTypes.emplace_back(cardinalityField->GetType());
      auto cardColReader = std::make_unique<ROOT::Experimental::Internal::RNTupleColumnReader>(
         std::move(cardinalityField));
      fColumnReaderPrototypes.emplace_back(std::move(cardColReader));
   }

   skeinIDs.emplace_back(fieldId);
   fColumnNames.emplace_back(colName);
   fColumnTypes.emplace_back(valueField->GetType());
   auto valColReader = std::make_unique<ROOT::Experimental::Internal::RNTupleColumnReader>(std::move(valueField));
   fColumnReaderPrototypes.emplace_back(std::move(valColReader));
}

RNTupleDS::RNTupleDS(std::unique_ptr<Detail::RPageSource> pageSource)
{
   pageSource->Attach();
   const auto &descriptor = pageSource->GetDescriptor();
   fSources.emplace_back(std::move(pageSource));

   AddField(descriptor, "", descriptor.GetFieldZeroId(), std::vector<DescriptorId_t>());
}

RDF::RDataSource::Record_t RNTupleDS::GetColumnReadersImpl(std::string_view /* name */, const std::type_info & /* ti */)
{
   // This datasource uses the GetColumnReaders2 API instead (better name in the works)
   return {};
}

std::unique_ptr<ROOT::Detail::RDF::RColumnReaderBase>
RNTupleDS::GetColumnReaders(unsigned int slot, std::string_view name, const std::type_info & /*tid*/)
{
   // at this point we can assume that `name` will be found in fColumnNames, RDF is in charge validation
   // TODO(jblomer): check incoming type
   const auto index = std::distance(fColumnNames.begin(), std::find(fColumnNames.begin(), fColumnNames.end(), name));
   auto clone = fColumnReaderPrototypes[index]->Clone();
   clone->Connect(*fSources[slot]);
   return clone;
}

bool RNTupleDS::SetEntry(unsigned int, ULong64_t)
{
   return true;
}

std::vector<std::pair<ULong64_t, ULong64_t>> RNTupleDS::GetEntryRanges()
{
   // TODO(jblomer): use cluster boundaries for the entry ranges
   std::vector<std::pair<ULong64_t, ULong64_t>> ranges;
   if (fHasSeenAllRanges)
      return ranges;

   auto nEntries = fSources[0]->GetNEntries();
   const auto chunkSize = nEntries / fNSlots;
   const auto reminder = 1U == fNSlots ? 0 : nEntries % fNSlots;
   auto start = 0UL;
   auto end = 0UL;
   for (auto i : ROOT::TSeqU(fNSlots)) {
      start = end;
      end += chunkSize;
      ranges.emplace_back(start, end);
      (void)i;
   }
   ranges.back().second += reminder;
   fHasSeenAllRanges = true;
   return ranges;
}

std::string RNTupleDS::GetTypeName(std::string_view colName) const
{
   const auto index = std::distance(fColumnNames.begin(), std::find(fColumnNames.begin(), fColumnNames.end(), colName));
   return fColumnTypes[index];
}

bool RNTupleDS::HasColumn(std::string_view colName) const
{
   return std::find(fColumnNames.begin(), fColumnNames.end(), colName) != fColumnNames.end();
}

void RNTupleDS::Initialise()
{
   fHasSeenAllRanges = false;
}

void RNTupleDS::Finalise() {}

void RNTupleDS::SetNSlots(unsigned int nSlots)
{
   R__ASSERT(fNSlots == 0);
   R__ASSERT(nSlots > 0);
   fNSlots = nSlots;

   for (unsigned int i = 1; i < fNSlots; ++i) {
      fSources.emplace_back(fSources[0]->Clone());
      R__DEBUG_ASSERT(i == (fSources.size() - 1));
      fSources[i]->Attach();
   }
}
} // namespace Experimental
} // namespace ROOT

ROOT::RDataFrame ROOT::Experimental::MakeNTupleDataFrame(std::string_view ntupleName, std::string_view fileName)
{
   auto pageSource = ROOT::Experimental::Detail::RPageSource::Create(ntupleName, fileName);
   ROOT::RDataFrame rdf(std::make_unique<RNTupleDS>(std::move(pageSource)));
   return rdf;
}
