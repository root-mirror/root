/// \file RNTupleModel.cxx
/// \ingroup NTuple ROOT7
/// \author Jakob Blomer <jblomer@cern.ch>
/// \date 2018-10-15
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/RError.hxx>
#include <ROOT/RField.hxx>
#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RNTuple.hxx>

#include <cstdlib>
#include <memory>
#include <utility>


void ROOT::Experimental::RNTupleModel::EnsureValidFieldName(std::string_view fieldName)
{
   RResult<void> nameValid = Detail::RFieldBase::EnsureValidFieldName(fieldName);
   if (!nameValid) {
      nameValid.Throw();
   }
   auto fieldNameStr = std::string(fieldName);
   if (fFieldNames.insert(fieldNameStr).second == false) {
      throw RException(R__FAIL("field name '" + fieldNameStr + "' already exists in NTuple model"));
   }
}

ROOT::Experimental::RNTupleModel::RNTupleModel()
  : fFieldZero(std::make_unique<RFieldZero>())
  , fDefaultEntry(std::make_unique<REntry>())
{}

std::unique_ptr<ROOT::Experimental::RNTupleModel> ROOT::Experimental::RNTupleModel::Clone() const
{
   auto cloneModel = std::make_unique<RNTupleModel>();
   auto cloneFieldZero = fFieldZero->Clone("");
   cloneModel->fFieldZero = std::unique_ptr<RFieldZero>(static_cast<RFieldZero *>(cloneFieldZero.release()));
   cloneModel->fDefaultEntry = cloneModel->fFieldZero->GenerateEntry();
   cloneModel->fFieldNames = fFieldNames;
   cloneModel->fDescription = fDescription;
   return cloneModel;
}


void ROOT::Experimental::RNTupleModel::AddField(std::unique_ptr<Detail::RFieldBase> field)
{
   if (!field) {
      throw RException(R__FAIL("null field"));
   }
   EnsureValidFieldName(field->GetName());
   fDefaultEntry->AddValue(field->GenerateValue());
   fFieldZero->Attach(std::move(field));
}


std::shared_ptr<ROOT::Experimental::RCollectionNTupleWriter> ROOT::Experimental::RNTupleModel::MakeCollection(
   std::string_view fieldName, std::unique_ptr<RNTupleModel> collectionModel)
{
   EnsureValidFieldName(fieldName);
   if (!collectionModel) {
      throw RException(R__FAIL("null collectionModel"));
   }
   auto collectionNTuple = std::make_shared<RCollectionNTupleWriter>(std::move(collectionModel->fDefaultEntry));
   auto field = std::make_unique<RCollectionField>(fieldName, collectionNTuple, std::move(collectionModel));
   fDefaultEntry->CaptureValue(field->CaptureValue(collectionNTuple->GetOffsetPtr()));
   fFieldZero->Attach(std::move(field));
   return collectionNTuple;
}

ROOT::Experimental::Detail::RFieldBase *ROOT::Experimental::RNTupleModel::GetField(std::string_view fieldName)
{
   auto split_on_dot = [](std::string_view in) {
      std::vector<std::string_view> splits;
      for (std::size_t pos = in.find('.'); pos != std::string_view::npos; pos = in.find('.')) {
         splits.push_back(in.substr(0, pos));
         in.remove_prefix(pos + 1);
      }
      if (in.size() > 0)
         splits.emplace_back(in);

      return splits;
   };

   auto *field = static_cast<ROOT::Experimental::Detail::RFieldBase *>(fFieldZero.get());
   for (auto subfieldName : split_on_dot(fieldName)) {
      const auto it =
         std::find_if(field->begin(), field->end(), [&](const auto &f) { return f.GetName() == subfieldName; });
      if (it != field->end()) {
         field = std::addressof(*it);
      } else {
         field = nullptr;
         break;
      }
   }

   return field;
}

std::unique_ptr<ROOT::Experimental::REntry> ROOT::Experimental::RNTupleModel::CreateEntry()
{
   auto entry = std::make_unique<REntry>();
   for (auto& f : *fFieldZero) {
      if (f.GetParent() != GetFieldZero())
         continue;
      entry->AddValue(f.GenerateValue());
   }
   return entry;
}
