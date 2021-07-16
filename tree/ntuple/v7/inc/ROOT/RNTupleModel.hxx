/// \file ROOT/RNTupleModel.hxx
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

#ifndef ROOT7_RNTupleModel
#define ROOT7_RNTupleModel

#include <ROOT/REntry.hxx>
#include <ROOT/RField.hxx>
#include <ROOT/RNTupleUtil.hxx>
#include <ROOT/RFieldValue.hxx>
#include <ROOT/RStringView.hxx>

#include <memory>
#include <unordered_set>
#include <utility>

namespace ROOT {
namespace Experimental {

class RCollectionNTupleWriter;

// clang-format off
/**
\class ROOT::Experimental::RNTupleModel
\ingroup NTuple
\brief The RNTupleModel encapulates the schema of an ntuple.

The ntuple model comprises a collection of hierarchically organized fields. From a model, "entries"
can be extracted. For convenience, the model provides a default entry. Models have a unique model identifier
that faciliates checking whether entries are compatible with it (i.e.: have been extracted from that model).
*/
// clang-format on
class RNTupleModel {
   /// Hierarchy of fields consisting of simple types and collections (sub trees)
   std::unique_ptr<RFieldZero> fFieldZero;
   /// Contains field values corresponding to the created top-level fields
   std::unique_ptr<REntry> fDefaultEntry;
   /// Keeps track of which field names are taken.
   std::unordered_set<std::string> fFieldNames;

   /// Checks that user-provided field names are valid in the context
   /// of this NTuple model. Throws an RException for invalid names.
   void EnsureValidFieldName(std::string_view fieldName);

   /// Free text set by the user
   std::string fDescription;

public:
   RNTupleModel();
   RNTupleModel(const RNTupleModel&) = delete;
   RNTupleModel& operator =(const RNTupleModel&) = delete;
   ~RNTupleModel() = default;

   std::unique_ptr<RNTupleModel> Clone() const;
   static std::unique_ptr<RNTupleModel> Create() { return std::make_unique<RNTupleModel>(); }

   /// Creates a new field and a corresponding tree value that is managed by a shared pointer.
   ///
   /// **Example: create some fields and fill an %RNTuple**
   /// ~~~ {.cpp}
   /// #include <ROOT/RNTuple.hxx>
   /// using ROOT::Experimental::RNTupleModel;
   /// using ROOT::Experimental::RNTupleWriter;
   ///
   /// #include <vector>
   ///
   /// auto model = RNTupleModel::Create();
   /// auto pt = model->MakeField<float>("pt");
   /// auto vec = model->MakeField<std::vector<int>>("vec");
   ///
   /// // The RNTuple is written to disk when the RNTupleWriter goes out of scope
   /// {
   ///    auto ntuple = RNTupleWriter::Recreate(std::move(model), "myNTuple", "myFile.root");
   ///    for (int i = 0; i < 100; i++) {
   ///       *pt = static_cast<float>(i);
   ///       *vec = {i, i+1, i+2};
   ///       ntuple->Fill();
   ///    }
   /// }
   /// ~~~
   /// **Example: create a field with an initial value**
   /// ~~~ {.cpp}
   /// #include <ROOT/RNTuple.hxx>
   /// using ROOT::Experimental::RNTupleModel;
   ///
   /// auto model = RNTupleModel::Create();
   /// // pt's initial value is 42.0
   /// auto pt = model->MakeField<float>("pt", 42.0);
   /// ~~~
   template <typename T, typename... ArgsT>
   std::shared_ptr<T> MakeField(std::string_view fieldName, ArgsT&&... args) {
      return MakeField<T>({fieldName, ""}, std::forward<ArgsT>(args)...);
   }

   /// Creates a new field given a `{name, description}` pair and a corresponding tree value that
   /// is managed by a shared pointer.
   ///
   /// **Example: create a field with a description**
   /// ~~~ {.cpp}
   /// #include <ROOT/RNTuple.hxx>
   /// using ROOT::Experimental::RNTupleModel;
   ///
   /// auto model = RNTupleModel::Create();
   /// auto hadronFlavour = model->MakeField<float>({
   ///    "hadronFlavour", "flavour from hadron ghost clustering"
   /// });
   /// ~~~
   template <typename T, typename... ArgsT>
   std::shared_ptr<T> MakeField(std::pair<std::string_view, std::string_view> fieldNameDesc,
      ArgsT&&... args)
   {
      EnsureValidFieldName(fieldNameDesc.first);
      auto field = std::make_unique<RField<T>>(fieldNameDesc.first);
      field->SetDescription(fieldNameDesc.second);
      auto ptr = fDefaultEntry->AddValue<T>(field.get(), std::forward<ArgsT>(args)...);
      fFieldZero->Attach(std::move(field));
      return ptr;
   }

   /// Adds a field whose type is not known at compile time.  Thus there is no shared pointer returned.
   ///
   /// Throws an exception if the field is null.
   void AddField(std::unique_ptr<Detail::RFieldBase> field);

   /// Throws an exception if fromWhere is null.
   template <typename T>
   void AddField(std::string_view fieldName, T* fromWhere) {
      AddField<T>({fieldName, ""}, fromWhere);
   }

   /// Throws an exception if fromWhere is null.
   template <typename T>
   void AddField(std::pair<std::string_view, std::string_view> fieldNameDesc, T* fromWhere) {
      EnsureValidFieldName(fieldNameDesc.first);
      if (!fromWhere) {
         throw RException(R__FAIL("null field fromWhere"));
      }
      auto field = std::make_unique<RField<T>>(fieldNameDesc.first);
      field->SetDescription(fieldNameDesc.second);
      fDefaultEntry->CaptureValue(field->CaptureValue(fromWhere));
      fFieldZero->Attach(std::move(field));
   }

   template <typename T>
   T* Get(std::string_view fieldName) {
      return fDefaultEntry->Get<T>(fieldName);
   }

   /// Ingests a model for a sub collection and attaches it to the current model
   ///
   /// Throws an exception if collectionModel is null.
   std::shared_ptr<RCollectionNTupleWriter> MakeCollection(
      std::string_view fieldName,
      std::unique_ptr<RNTupleModel> collectionModel);

   RFieldZero *GetFieldZero() const { return fFieldZero.get(); }
   Detail::RFieldBase *GetField(std::string_view fieldName);
   REntry *GetDefaultEntry() { return fDefaultEntry.get(); }
   std::unique_ptr<REntry> CreateEntry();
   RNTupleVersion GetVersion() const { return RNTupleVersion(); }
   std::string GetDescription() const { return fDescription; }
   void SetDescription(std::string_view description) { fDescription = std::string(description); }
   RNTupleUuid GetUuid() const { return RNTupleUuid(); /* TODO */ }
};

} // namespace Experimental
} // namespace ROOT

#endif
