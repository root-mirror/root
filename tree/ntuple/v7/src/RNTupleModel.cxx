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
#include <ROOT/RFieldVisitor.hxx>
#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RNTuple.hxx>

#include <cstdlib>
#include <memory>
#include <sstream>
#include <utility>

namespace {

class RHasFieldVisitor : public ROOT::Experimental::Detail::RFieldVisitor {
private:
   std::string fRest;
   bool fHasField = false;
   bool StopRecurse() {
      return fRest.empty();
   }
   std::string PopTopLevelField() {
      std::istringstream iss(fRest);
      std::string topLevelField;
      std::getline(iss, topLevelField, '.');
      std::string rest;
      std::getline(iss, rest);
      fRest = rest;
      return topLevelField;
   }
public:
   RHasFieldVisitor(std::string_view fieldName) : fRest(std::string(fieldName)) {}
   bool HasField() const {
      return fHasField;
   }
   void VisitField(const ROOT::Experimental::Detail::RFieldBase &field) final {
      auto fieldName = PopTopLevelField();
      for (auto f: field.GetSubFields()) {
         if (fieldName == f->GetName()) {
            if (StopRecurse()) {
               fHasField = true;
               return;
            }
            VisitField(*f);
         }
      }
   }
};

} // anonymous namespace

void ROOT::Experimental::RNTupleModel::EnsureValidFieldName(std::string_view fieldName)
{
   RResult<void> nameValid = Detail::RFieldBase::EnsureValidFieldName(fieldName);
   if (!nameValid) {
      nameValid.Throw();
   }
   RHasFieldVisitor vis(fieldName);
   fFieldZero->AcceptVisitor(vis);
   if (vis.HasField()) {
      throw RException(R__FAIL(
         "field name '" + std::string(fieldName) + "' already exists in NTuple model"
      ));
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
   if (!field) {
      throw RException(R__FAIL("null field"));
   }
   EnsureValidFieldName(field->GetName());
   fDefaultEntry->AddValue(field->GenerateValue());
   fFieldZero->Attach(std::move(field));
}

bool ROOT::Experimental::RNTupleModel::HasField(std::string_view fieldName) const {
   RHasFieldVisitor vis(fieldName);
   fFieldZero->AcceptVisitor(vis);
   return vis.HasField();
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
