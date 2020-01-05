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

#include <TError.h>

#include <memory>
#include <utility>

namespace ROOT {
namespace Experimental {

class RCollectionNTuple;

// clang-format off
/**
\class ROOT::Experimental::RNTupleModel
\ingroup NTuple
\brief The RNTupleModel encapulates the schema of an ntuple.

The ntuple model comprises a collection of hierarchically organized fields. From a frozen model, "entries"
can be extracted. For convenience, the model provides a default entry. Models have a unique model identifier
that faciliates checking whether entries are compatible with it (i.e.: have been extracted from that model).
A model needs to be frozen before it can be used to create a live ntuple.
*/
// clang-format on
class RNTupleModel {
   /// Hierarchy of fields consisting of simple types and collections (sub trees)
   std::unique_ptr<RFieldRoot> fRootField;
   /// Contains field values corresponding to the created top-level fields
   std::unique_ptr<REntry> fDefaultEntry;

public:
   RNTupleModel();
   RNTupleModel(const RNTupleModel&) = delete;
   RNTupleModel& operator =(const RNTupleModel&) = delete;
   ~RNTupleModel() = default;

   RNTupleModel* Clone();
   static std::unique_ptr<RNTupleModel> Create() { return std::make_unique<RNTupleModel>(); }

   /// Creates a new field and a corresponding tree value that is managed by a shared pointer.
   template <typename T, std::size_t nBits, std::int64_t min, std::int64_t max, typename... ArgsT>
   std::shared_ptr<T> MakeField(std::string_view fieldName, ArgsT&&... args) {
      if (nBits > 64) {
         std::cout << "The number of bits should be 64 or lower. Nothing was attached to the model and a nullptr was returned." << std::endl;
         return nullptr;
      } else if (nBits < 3) {
         std::cout << "The number of bits should be 3 or higher. Nothing was attached to the model and a nullptr was returned." << std::endl;
         return nullptr;
      }
      auto field = std::make_unique<RField<T, RCustomSizedFloat> >(fieldName, nBits, min, max);
      auto ptr = fDefaultEntry->AddValueEncoder<T>(field.get(), std::forward<ArgsT>(args)...);
      fRootField->Attach(std::move(field));
      return ptr;
   }

   /// Creates a new field and a corresponding tree value that is managed by a shared pointer.
   template <typename T, typename... ArgsT>
   std::shared_ptr<T> MakeField(std::string_view fieldName, ArgsT&&... args) {
      auto field = std::make_unique<RField<T>>(fieldName);
      auto ptr = fDefaultEntry->AddValue<T>(field.get(), std::forward<ArgsT>(args)...);
      fRootField->Attach(std::move(field));
      return ptr;
   }

   /// Adds a field whose type is not known at compile time.  Thus there is no shared pointer returned.
   void AddField(std::unique_ptr<Detail::RFieldBase> field);

   template <typename T>
   void AddField(std::string_view fieldName, T* fromWhere) {
      auto field = std::make_unique<RField<T>>(fieldName);
      fDefaultEntry->CaptureValue(field->CaptureValue(fromWhere));
      fRootField->Attach(std::move(field));
   }
   template <typename T, std::size_t nBits, std::int64_t min, std::int64_t max>
   void AddField(std::string_view fieldName, T* fromWhere) {
      auto field = std::make_unique<RField<T, RCustomSizedFloat>>(fieldName, nBits, min, max);
      fDefaultEntry->CaptureValue(field->CaptureValue(fromWhere));
      fRootField->Attach(std::move(field));
   }

   template <typename T>
   T* Get(std::string_view fieldName) {
      return fDefaultEntry->Get<T>(fieldName);
   }

   /// Ingests a model for a sub collection and attaches it to the current model
   std::shared_ptr<RCollectionNTuple> MakeCollection(
      std::string_view fieldName,
      std::unique_ptr<RNTupleModel> collectionModel);

   RFieldRoot *GetRootField() const { return fRootField.get(); }
   REntry* GetDefaultEntry() { return fDefaultEntry.get(); }
   std::unique_ptr<REntry> CreateEntry();
   RNTupleVersion GetVersion() const { return RNTupleVersion(); }
   std::string GetDescription() const { return ""; /* TODO */ }
   RNTupleUuid GetUuid() const { return RNTupleUuid(); /* TODO */ }
};

} // namespace Experimental
} // namespace ROOT

#endif
