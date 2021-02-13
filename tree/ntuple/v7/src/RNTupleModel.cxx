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
   return cloneModel;
}


void ROOT::Experimental::RNTupleModel::AddField(std::unique_ptr<Detail::RFieldBase> field)
{
   EnsureValidFieldName(field->GetName());
   fDefaultEntry->AddValue(field->GenerateValue());
   fFieldZero->Attach(std::move(field));
}


std::shared_ptr<ROOT::Experimental::RCollectionNTuple> ROOT::Experimental::RNTupleModel::MakeCollection(
   std::string_view fieldName, std::unique_ptr<RNTupleModel> collectionModel)
{
   EnsureValidFieldName(fieldName);
   auto collectionNTuple = std::make_shared<RCollectionNTuple>(std::move(collectionModel->fDefaultEntry));
   auto field = std::make_unique<RCollectionField>(fieldName, collectionNTuple, std::move(collectionModel));
   fDefaultEntry->CaptureValue(field->CaptureValue(collectionNTuple->GetOffsetPtr()));
   fFieldZero->Attach(std::move(field));
   return collectionNTuple;
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
