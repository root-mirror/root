/// \file RField.cxx
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

#include <ROOT/RColumn.hxx>
#include <ROOT/RColumnModel.hxx>
#include <ROOT/REntry.hxx>
#include <ROOT/RError.hxx>
#include <ROOT/RField.hxx>
#include <ROOT/RFieldValue.hxx>
#include <ROOT/RFieldVisitor.hxx>
#include <ROOT/RLogger.hxx>
#include <ROOT/RNTuple.hxx>
#include <ROOT/RNTupleModel.hxx>

#include <TClass.h>
#include <TCollection.h>
#include <TDataMember.h>
#include <TError.h>
#include <TList.h>

#include <algorithm>
#include <cctype> // for isspace
#include <cstdlib> // for malloc, free
#include <cstring> // for memset
#include <exception>
#include <iostream>
#include <type_traits>
#include <unordered_map>

namespace {

/// Used in CreateField() in order to get the comma-separated list of template types
/// E.g., gets {"int", "std::variant<double,int>"} from "int,std::variant<double,int>"
std::vector<std::string> TokenizeTypeList(std::string templateType) {
   std::vector<std::string> result;
   if (templateType.empty())
      return result;

   const char *eol = templateType.data() + templateType.length();
   const char *typeBegin = templateType.data();
   const char *typeCursor = templateType.data();
   unsigned int nestingLevel = 0;
   while (typeCursor != eol) {
      switch (*typeCursor) {
      case '<':
         ++nestingLevel;
         break;
      case '>':
         --nestingLevel;
         break;
      case ',':
         if (nestingLevel == 0) {
            result.push_back(std::string(typeBegin, typeCursor - typeBegin));
            typeBegin = typeCursor + 1;
         }
         break;
      }
      typeCursor++;
   }
   result.push_back(std::string(typeBegin, typeCursor - typeBegin));
   return result;
}

/// Remove leading and trailing white spaces
std::string Trim(const std::string &raw) {
  if (raw.empty()) return "";

  unsigned start_pos = 0;
  for (; (start_pos < raw.length()) && (raw[start_pos] == ' ' || raw[start_pos] == '\t'); ++start_pos) { }

  unsigned end_pos = raw.length() - 1;  // at least one character in raw
  for (; (end_pos >= start_pos) && (raw[end_pos] == ' ' || raw[end_pos] == '\t'); --end_pos) { }

  return raw.substr(start_pos, end_pos - start_pos + 1);
}

std::string GetNormalizedType(const std::string &typeName) {
   std::string normalizedType(Trim(typeName));

   // TODO(jblomer): use a type translation map
   if (normalizedType == "Bool_t") normalizedType = "bool";
   if (normalizedType == "Float_t") normalizedType = "float";
   if (normalizedType == "Double_t") normalizedType = "double";
   if (normalizedType == "UChar_t") normalizedType = "std::uint8_t";
   if (normalizedType == "unsigned char") normalizedType = "std::uint8_t";
   if (normalizedType == "uint8_t") normalizedType = "std::uint8_t";
   if (normalizedType == "Int_t") normalizedType = "std::int32_t";
   if (normalizedType == "int") normalizedType = "std::int32_t";
   if (normalizedType == "int32_t") normalizedType = "std::int32_t";
   if (normalizedType == "unsigned") normalizedType = "std::uint32_t";
   if (normalizedType == "unsigned int") normalizedType = "std::uint32_t";
   if (normalizedType == "UInt_t") normalizedType = "std::uint32_t";
   if (normalizedType == "uint32_t") normalizedType = "std::uint32_t";
   if (normalizedType == "ULong64_t") normalizedType = "std::uint64_t";
   if (normalizedType == "uint64_t") normalizedType = "std::uint64_t";
   if (normalizedType == "string") normalizedType = "std::string";
   if (normalizedType.substr(0, 7) == "vector<") normalizedType = "std::" + normalizedType;
   if (normalizedType.substr(0, 6) == "array<") normalizedType = "std::" + normalizedType;
   if (normalizedType.substr(0, 8) == "variant<") normalizedType = "std::" + normalizedType;

   return normalizedType;
}

} // anonymous namespace

void ROOT::Experimental::Detail::RFieldFuse::Connect(DescriptorId_t fieldId, RPageStorage &pageStorage, RFieldBase &field)
{
   if (field.fColumns.empty())
      field.GenerateColumnsImpl();
   for (auto& column : field.fColumns)
      column->Connect(fieldId, &pageStorage);
}


void ROOT::Experimental::Detail::RFieldFuse::ConnectRecursively(
   DescriptorId_t fieldId, RPageSource &pageSource, RFieldBase &field)
{
   Connect(fieldId, pageSource, field);
   std::unordered_map<const RFieldBase *, DescriptorId_t> field2Id;
   field2Id[&field] = fieldId;
   for (auto &f : field) {
      auto subFieldId = pageSource.GetDescriptor().FindFieldId(f.GetName(), field2Id[f.GetParent()]);
      Detail::RFieldFuse::Connect(subFieldId, pageSource, f);
      field2Id[&f] = subFieldId;
   }
}


//------------------------------------------------------------------------------


ROOT::Experimental::Detail::RFieldBase::RFieldBase(
   std::string_view name, std::string_view type, ENTupleStructure structure, bool isSimple, std::size_t nRepetitions)
   : fName(name), fType(type), fStructure(structure), fNRepetitions(nRepetitions), fIsSimple(isSimple),
     fParent(nullptr), fPrincipalColumn(nullptr)
{
}

ROOT::Experimental::Detail::RFieldBase::~RFieldBase()
{
}

ROOT::Experimental::RResult<std::unique_ptr<ROOT::Experimental::Detail::RFieldBase>>
ROOT::Experimental::Detail::RFieldBase::Create(const std::string &fieldName, const std::string &typeName)
{
   std::string normalizedType(GetNormalizedType(typeName));
   if (normalizedType.empty())
      return R__FAIL("no type name specified for Field " + fieldName);

   std::unique_ptr<ROOT::Experimental::Detail::RFieldBase> result;

   if (normalizedType == "ROOT::Experimental::ClusterSize_t") {
      result = std::make_unique<RField<ClusterSize_t>>(fieldName);
   } else if (normalizedType == "bool") {
      result = std::make_unique<RField<bool>>(fieldName);
   } else if (normalizedType == "std::uint8_t") {
      result = std::make_unique<RField<std::uint8_t>>(fieldName);
   } else if (normalizedType == "std::int32_t") {
      result = std::make_unique<RField<std::int32_t>>(fieldName);
   } else if (normalizedType == "std::uint32_t") {
      result = std::make_unique<RField<std::uint32_t>>(fieldName);
   } else if (normalizedType == "std::uint64_t") {
      result = std::make_unique<RField<std::uint64_t>>(fieldName);
   } else if (normalizedType == "float") {
      result = std::make_unique<RField<float>>(fieldName);
   } else if (normalizedType == "double") {
      result = std::make_unique<RField<double>>(fieldName);
   } else if (normalizedType == "std::string") {
      result = std::make_unique<RField<std::string>>(fieldName);
   } else if (normalizedType == "std::vector<bool>") {
      result = std::make_unique<RField<std::vector<bool>>>(fieldName);
   } else if (normalizedType.substr(0, 12) == "std::vector<") {
      std::string itemTypeName = normalizedType.substr(12, normalizedType.length() - 13);
      auto itemField = Create(GetNormalizedType(itemTypeName), itemTypeName);
      result = std::make_unique<RVectorField>(fieldName, itemField.Unwrap());
   } else if (normalizedType == "ROOT::VecOps::RVec<bool>") {
      result = std::make_unique<RField<ROOT::VecOps::RVec<bool>>>(fieldName);
   } else if (normalizedType.substr(0, 19) == "ROOT::VecOps::RVec<") {
      // For the time being, we silently read RVec fields as std::vector
      std::string itemTypeName = normalizedType.substr(19, normalizedType.length() - 20);
      auto itemField = Create(GetNormalizedType(itemTypeName), itemTypeName);
      result = std::make_unique<RVectorField>(fieldName, itemField.Unwrap());
   } else if (normalizedType.substr(0, 11) == "std::array<") {
      auto arrayDef = TokenizeTypeList(normalizedType.substr(11, normalizedType.length() - 12));
      R__ASSERT(arrayDef.size() == 2);
      auto arrayLength = std::stoi(arrayDef[1]);
      auto itemField = Create(GetNormalizedType(arrayDef[0]), arrayDef[0]);
      result = std::make_unique<RArrayField>(fieldName, itemField.Unwrap(), arrayLength);
   }
#if __cplusplus >= 201703L
   if (normalizedType.substr(0, 13) == "std::variant<") {
      auto innerTypes = TokenizeTypeList(normalizedType.substr(13, normalizedType.length() - 14));
      std::vector<RFieldBase *> items;
      for (unsigned int i = 0; i < innerTypes.size(); ++i) {
         items.emplace_back(Create("variant" + std::to_string(i), innerTypes[i]).Unwrap().release());
      }
      result = std::make_unique<RVariantField>(fieldName, items);
   }
#endif
   // TODO: create an RCollectionField?
   if (normalizedType == ":Collection:")
     result = std::make_unique<RField<ClusterSize_t>>(fieldName);

   if (!result) {
      auto cl = TClass::GetClass(normalizedType.c_str());
      if (cl != nullptr) {
         result = std::make_unique<RClassField>(fieldName, normalizedType);
      }
   }

   if (result)
      return result;
   return R__FAIL(std::string("Field ") + fieldName + " has unknown type " + normalizedType);
}

ROOT::Experimental::RResult<void>
ROOT::Experimental::Detail::RFieldBase::EnsureValidFieldName(std::string_view fieldName)
{
   if (fieldName == "") {
      return R__FAIL("name cannot be empty string \"\"");
   } else if (fieldName.find(".") != std::string::npos) {
      return R__FAIL("name '" + std::string(fieldName) + "' cannot contain dot characters '.'");
   }
   return RResult<void>::Success();
}

void ROOT::Experimental::Detail::RFieldBase::AppendImpl(const ROOT::Experimental::Detail::RFieldValue& /*value*/) {
   R__ASSERT(false);
}

void ROOT::Experimental::Detail::RFieldBase::ReadGlobalImpl(
   ROOT::Experimental::NTupleSize_t /*index*/,
   RFieldValue* /*value*/)
{
   R__ASSERT(false);
}

ROOT::Experimental::Detail::RFieldValue ROOT::Experimental::Detail::RFieldBase::GenerateValue()
{
   void *where = malloc(GetValueSize());
   R__ASSERT(where != nullptr);
   return GenerateValue(where);
}

void ROOT::Experimental::Detail::RFieldBase::DestroyValue(const RFieldValue &value, bool dtorOnly)
{
   if (!dtorOnly)
      free(value.GetRawPtr());
}

std::vector<ROOT::Experimental::Detail::RFieldValue>
ROOT::Experimental::Detail::RFieldBase::SplitValue(const RFieldValue & /*value*/) const
{
   return std::vector<RFieldValue>();
}

void ROOT::Experimental::Detail::RFieldBase::Attach(
   std::unique_ptr<ROOT::Experimental::Detail::RFieldBase> child)
{
   child->fParent = this;
   fSubFields.emplace_back(std::move(child));
}


std::vector<const ROOT::Experimental::Detail::RFieldBase *>
ROOT::Experimental::Detail::RFieldBase::GetSubFields() const
{
   std::vector<const RFieldBase *> result;
   for (const auto &f : fSubFields) {
      result.emplace_back(f.get());
   }
   return result;
}


void ROOT::Experimental::Detail::RFieldBase::Flush() const
{
   for (auto& column : fColumns) {
      column->Flush();
   }
}


void ROOT::Experimental::Detail::RFieldBase::AcceptVisitor(Detail::RFieldVisitor &visitor) const
{
   visitor.VisitField(*this);
}


ROOT::Experimental::Detail::RFieldBase::RSchemaIterator ROOT::Experimental::Detail::RFieldBase::begin()
{
   if (fSubFields.empty()) return RSchemaIterator(this, -1);
   return RSchemaIterator(this->fSubFields[0].get(), 0);
}


ROOT::Experimental::Detail::RFieldBase::RSchemaIterator ROOT::Experimental::Detail::RFieldBase::end()
{
   return RSchemaIterator(this, -1);
}


//-----------------------------------------------------------------------------


void ROOT::Experimental::Detail::RFieldBase::RSchemaIterator::Advance()
{
   auto itr = fStack.rbegin();
   if (!itr->fFieldPtr->fSubFields.empty()) {
      fStack.emplace_back(Position(itr->fFieldPtr->fSubFields[0].get(), 0));
      return;
   }

   unsigned int nextIdxInParent = ++(itr->fIdxInParent);
   while (nextIdxInParent >= itr->fFieldPtr->fParent->fSubFields.size()) {
      if (fStack.size() == 1) {
         itr->fFieldPtr = itr->fFieldPtr->fParent;
         itr->fIdxInParent = -1;
         return;
      }
      fStack.pop_back();
      itr = fStack.rbegin();
      nextIdxInParent = ++(itr->fIdxInParent);
   }
   itr->fFieldPtr = itr->fFieldPtr->fParent->fSubFields[nextIdxInParent].get();
}


//------------------------------------------------------------------------------


std::unique_ptr<ROOT::Experimental::Detail::RFieldBase>
ROOT::Experimental::RFieldZero::Clone(std::string_view /*newName*/) const
{
   auto result = std::make_unique<RFieldZero>();
   for (auto &f : fSubFields)
      result->Attach(f->Clone(f->GetName()));
   return result;
}


std::unique_ptr<ROOT::Experimental::REntry> ROOT::Experimental::RFieldZero::GenerateEntry() const
{
   auto entry = std::make_unique<REntry>();
   for (auto& f : fSubFields) {
      entry->AddValue(f->GenerateValue());
   }
   return entry;
}

void ROOT::Experimental::RFieldZero::AcceptVisitor(Detail::RFieldVisitor &visitor) const
{
   visitor.VisitFieldZero(*this);
}


//------------------------------------------------------------------------------


void ROOT::Experimental::RField<ROOT::Experimental::ClusterSize_t>::GenerateColumnsImpl()
{
   RColumnModel model(EColumnType::kIndex, true /* isSorted*/);
   fColumns.emplace_back(std::unique_ptr<Detail::RColumn>(
      Detail::RColumn::Create<ClusterSize_t, EColumnType::kIndex>(model, 0)));
   fPrincipalColumn = fColumns[0].get();
}

void ROOT::Experimental::RField<ROOT::Experimental::ClusterSize_t>::AcceptVisitor(Detail::RFieldVisitor &visitor) const
{
   visitor.VisitClusterSizeField(*this);
}

//------------------------------------------------------------------------------

void ROOT::Experimental::RField<std::uint8_t>::GenerateColumnsImpl()
{
   RColumnModel model(EColumnType::kByte, false /* isSorted*/);
   fColumns.emplace_back(std::unique_ptr<Detail::RColumn>(Detail::RColumn::Create<
      std::uint8_t, EColumnType::kByte>(model, 0)));
   fPrincipalColumn = fColumns[0].get();
}

void ROOT::Experimental::RField<std::uint8_t>::AcceptVisitor(Detail::RFieldVisitor &visitor) const
{
   visitor.VisitUInt8Field(*this);
}

//------------------------------------------------------------------------------


void ROOT::Experimental::RField<bool>::GenerateColumnsImpl()
{
   RColumnModel model(EColumnType::kBit, false /* isSorted*/);
   fColumns.emplace_back(std::unique_ptr<Detail::RColumn>(
      Detail::RColumn::Create<bool, EColumnType::kBit>(model, 0)));
   fPrincipalColumn = fColumns[0].get();
}

void ROOT::Experimental::RField<bool>::AcceptVisitor(Detail::RFieldVisitor &visitor) const
{
   visitor.VisitBoolField(*this);
}

//------------------------------------------------------------------------------


void ROOT::Experimental::RField<float>::GenerateColumnsImpl()
{
   RColumnModel model(EColumnType::kReal32, false /* isSorted*/);
   fColumns.emplace_back(std::unique_ptr<Detail::RColumn>(
      Detail::RColumn::Create<float, EColumnType::kReal32>(model, 0)));
   fPrincipalColumn = fColumns[0].get();
}

void ROOT::Experimental::RField<float>::AcceptVisitor(Detail::RFieldVisitor &visitor) const
{
   visitor.VisitFloatField(*this);
}


//------------------------------------------------------------------------------

void ROOT::Experimental::RField<double>::GenerateColumnsImpl()
{
   RColumnModel model(EColumnType::kReal64, false /* isSorted*/);
   fColumns.emplace_back(std::unique_ptr<Detail::RColumn>(
      Detail::RColumn::Create<double, EColumnType::kReal64>(model, 0)));
   fPrincipalColumn = fColumns[0].get();
}

void ROOT::Experimental::RField<double>::AcceptVisitor(Detail::RFieldVisitor &visitor) const
{
   visitor.VisitDoubleField(*this);
}

//------------------------------------------------------------------------------

void ROOT::Experimental::RField<std::int32_t>::GenerateColumnsImpl()
{
   RColumnModel model(EColumnType::kInt32, false /* isSorted*/);
   fColumns.emplace_back(std::unique_ptr<Detail::RColumn>(Detail::RColumn::Create<
      std::int32_t, EColumnType::kInt32>(model, 0)));
   fPrincipalColumn = fColumns[0].get();
}

void ROOT::Experimental::RField<std::int32_t>::AcceptVisitor(Detail::RFieldVisitor &visitor) const
{
   visitor.VisitIntField(*this);
}

//------------------------------------------------------------------------------

void ROOT::Experimental::RField<std::uint32_t>::GenerateColumnsImpl()
{
   RColumnModel model(EColumnType::kInt32, false /* isSorted*/);
   fColumns.emplace_back(std::unique_ptr<Detail::RColumn>(
      Detail::RColumn::Create<std::uint32_t, EColumnType::kInt32>(model, 0)));
   fPrincipalColumn = fColumns[0].get();
}

void ROOT::Experimental::RField<std::uint32_t>::AcceptVisitor(Detail::RFieldVisitor &visitor) const
{
   visitor.VisitUInt32Field(*this);
}

//------------------------------------------------------------------------------

void ROOT::Experimental::RField<std::uint64_t>::GenerateColumnsImpl()
{
   RColumnModel model(EColumnType::kInt64, false /* isSorted*/);
   fColumns.emplace_back(std::unique_ptr<Detail::RColumn>(
      Detail::RColumn::Create<std::uint64_t, EColumnType::kInt64>(model, 0)));
   fPrincipalColumn = fColumns[0].get();
}

void ROOT::Experimental::RField<std::uint64_t>::AcceptVisitor(Detail::RFieldVisitor &visitor) const
{
   visitor.VisitUInt64Field(*this);
}

//------------------------------------------------------------------------------


void ROOT::Experimental::RField<std::string>::GenerateColumnsImpl()
{
   RColumnModel modelIndex(EColumnType::kIndex, true /* isSorted*/);
   fColumns.emplace_back(std::unique_ptr<Detail::RColumn>(
      Detail::RColumn::Create<ClusterSize_t, EColumnType::kIndex>(modelIndex, 0)));

   RColumnModel modelChars(EColumnType::kByte, false /* isSorted*/);
   fColumns.emplace_back(std::unique_ptr<Detail::RColumn>(
      Detail::RColumn::Create<char, EColumnType::kByte>(modelChars, 1)));
   fPrincipalColumn = fColumns[0].get();
}

void ROOT::Experimental::RField<std::string>::AppendImpl(const ROOT::Experimental::Detail::RFieldValue& value)
{
   auto typedValue = value.Get<std::string>();
   auto length = typedValue->length();
   Detail::RColumnElement<char, EColumnType::kByte> elemChars(const_cast<char*>(typedValue->data()));
   fColumns[1]->AppendV(elemChars, length);
   fIndex += length;
   fColumns[0]->Append(fElemIndex);
}

void ROOT::Experimental::RField<std::string>::ReadGlobalImpl(
   ROOT::Experimental::NTupleSize_t globalIndex, ROOT::Experimental::Detail::RFieldValue *value)
{
   auto typedValue = value->Get<std::string>();
   RClusterIndex collectionStart;
   ClusterSize_t nChars;
   fPrincipalColumn->GetCollectionInfo(globalIndex, &collectionStart, &nChars);
   if (nChars == 0) {
      typedValue->clear();
   } else {
      typedValue->resize(nChars);
      Detail::RColumnElement<char, EColumnType::kByte> elemChars(const_cast<char*>(typedValue->data()));
      fColumns[1]->ReadV(collectionStart, nChars, &elemChars);
   }
}

void ROOT::Experimental::RField<std::string>::CommitCluster()
{
   fIndex = 0;
}

void ROOT::Experimental::RField<std::string>::AcceptVisitor(Detail::RFieldVisitor &visitor) const
{
   visitor.VisitStringField(*this);
}

//------------------------------------------------------------------------------


ROOT::Experimental::RClassField::RClassField(std::string_view fieldName, std::string_view className)
   : ROOT::Experimental::Detail::RFieldBase(fieldName, className, ENTupleStructure::kRecord, false /* isSimple */)
   , fClass(TClass::GetClass(std::string(className).c_str()))
{
   if (fClass == nullptr) {
      throw std::runtime_error("RField: no I/O support for type " + std::string(className));
   }
   TIter next(fClass->GetListOfDataMembers());
   while (auto dataMember = static_cast<TDataMember *>(next())) {
      //printf("Now looking at %s %s\n", dataMember->GetName(), dataMember->GetFullTypeName());
      auto subField = Detail::RFieldBase::Create(dataMember->GetName(), dataMember->GetFullTypeName()).Unwrap();
      fMaxAlignment = std::max(fMaxAlignment, subField->GetAlignment());
      Attach(std::move(subField));
   }
}

std::unique_ptr<ROOT::Experimental::Detail::RFieldBase>
ROOT::Experimental::RClassField::Clone(std::string_view newName) const
{
   return std::make_unique<RClassField>(newName, GetType());
}

void ROOT::Experimental::RClassField::AppendImpl(const Detail::RFieldValue& value) {
   TIter next(fClass->GetListOfDataMembers());
   unsigned i = 0;
   while (auto dataMember = static_cast<TDataMember *>(next())) {
      auto memberValue = fSubFields[i]->CaptureValue(value.Get<unsigned char>() + dataMember->GetOffset());
      fSubFields[i]->Append(memberValue);
      i++;
   }
}

void ROOT::Experimental::RClassField::ReadGlobalImpl(NTupleSize_t globalIndex, Detail::RFieldValue *value)
{
   TIter next(fClass->GetListOfDataMembers());
   unsigned i = 0;
   while (auto dataMember = static_cast<TDataMember *>(next())) {
      auto memberValue = fSubFields[i]->GenerateValue(value->Get<unsigned char>() + dataMember->GetOffset());
      fSubFields[i]->Read(globalIndex, &memberValue);
      i++;
   }
}

void ROOT::Experimental::RClassField::ReadInClusterImpl(const RClusterIndex &clusterIndex, Detail::RFieldValue *value)
{
   TIter next(fClass->GetListOfDataMembers());
   unsigned i = 0;
   while (auto dataMember = static_cast<TDataMember *>(next())) {
      auto memberValue = fSubFields[i]->GenerateValue(value->Get<unsigned char>() + dataMember->GetOffset());
      fSubFields[i]->Read(clusterIndex, &memberValue);
      i++;
   }
}

void ROOT::Experimental::RClassField::GenerateColumnsImpl()
{
}

ROOT::Experimental::Detail::RFieldValue ROOT::Experimental::RClassField::GenerateValue(void* where)
{
   return Detail::RFieldValue(true /* captureFlag */, this, fClass->New(where));
}

void ROOT::Experimental::RClassField::DestroyValue(const Detail::RFieldValue& value, bool dtorOnly)
{
   fClass->Destructor(value.GetRawPtr(), true /* dtorOnly */);
   if (!dtorOnly)
      free(value.GetRawPtr());
}

ROOT::Experimental::Detail::RFieldValue ROOT::Experimental::RClassField::CaptureValue(void* where)
{
   return Detail::RFieldValue(true /* captureFlat */, this, where);
}


std::vector<ROOT::Experimental::Detail::RFieldValue>
ROOT::Experimental::RClassField::SplitValue(const Detail::RFieldValue &value) const
{
   TIter next(fClass->GetListOfDataMembers());
   unsigned i = 0;
   std::vector<Detail::RFieldValue> result;
   while (auto dataMember = static_cast<TDataMember *>(next())) {
      auto memberValue = fSubFields[i]->CaptureValue(value.Get<unsigned char>() + dataMember->GetOffset());
      result.emplace_back(memberValue);
      i++;
   }
   return result;
}


size_t ROOT::Experimental::RClassField::GetValueSize() const
{
   return fClass->GetClassSize();
}

void ROOT::Experimental::RClassField::AcceptVisitor(Detail::RFieldVisitor &visitor) const
{
   visitor.VisitClassField(*this);
}

//------------------------------------------------------------------------------


ROOT::Experimental::RVectorField::RVectorField(
   std::string_view fieldName, std::unique_ptr<Detail::RFieldBase> itemField)
   : ROOT::Experimental::Detail::RFieldBase(
      fieldName, "std::vector<" + itemField->GetType() + ">", ENTupleStructure::kCollection, false /* isSimple */)
   , fItemSize(itemField->GetValueSize()), fNWritten(0)
{
   Attach(std::move(itemField));
}

std::unique_ptr<ROOT::Experimental::Detail::RFieldBase>
ROOT::Experimental::RVectorField::Clone(std::string_view newName) const
{
   auto newItemField = fSubFields[0]->Clone(fSubFields[0]->GetName());
   return std::make_unique<RVectorField>(newName, std::move(newItemField));
}

void ROOT::Experimental::RVectorField::AppendImpl(const Detail::RFieldValue& value) {
   auto typedValue = value.Get<std::vector<char>>();
   R__ASSERT((typedValue->size() % fItemSize) == 0);
   auto count = typedValue->size() / fItemSize;
   for (unsigned i = 0; i < count; ++i) {
      auto itemValue = fSubFields[0]->CaptureValue(typedValue->data() + (i * fItemSize));
      fSubFields[0]->Append(itemValue);
   }
   Detail::RColumnElement<ClusterSize_t, EColumnType::kIndex> elemIndex(&fNWritten);
   fNWritten += count;
   fColumns[0]->Append(elemIndex);
}

void ROOT::Experimental::RVectorField::ReadGlobalImpl(NTupleSize_t globalIndex, Detail::RFieldValue *value)
{
   auto typedValue = value->Get<std::vector<char>>();

   ClusterSize_t nItems;
   RClusterIndex collectionStart;
   fPrincipalColumn->GetCollectionInfo(globalIndex, &collectionStart, &nItems);

   typedValue->resize(nItems * fItemSize);
   for (unsigned i = 0; i < nItems; ++i) {
      auto itemValue = fSubFields[0]->GenerateValue(typedValue->data() + (i * fItemSize));
      fSubFields[0]->Read(collectionStart + i, &itemValue);
   }
}

void ROOT::Experimental::RVectorField::GenerateColumnsImpl()
{
   RColumnModel modelIndex(EColumnType::kIndex, true /* isSorted*/);
   fColumns.emplace_back(std::unique_ptr<Detail::RColumn>(
      Detail::RColumn::Create<ClusterSize_t, EColumnType::kIndex>(modelIndex, 0)));
   fPrincipalColumn = fColumns[0].get();
}

ROOT::Experimental::Detail::RFieldValue ROOT::Experimental::RVectorField::GenerateValue(void* where)
{
   return Detail::RFieldValue(this, reinterpret_cast<std::vector<char>*>(where));
}

void ROOT::Experimental::RVectorField::DestroyValue(const Detail::RFieldValue& value, bool dtorOnly)
{
   auto vec = static_cast<std::vector<char>*>(value.GetRawPtr());
   R__ASSERT((vec->size() % fItemSize) == 0);
   auto nItems = vec->size() / fItemSize;
   for (unsigned i = 0; i < nItems; ++i) {
      auto itemValue = fSubFields[0]->CaptureValue(vec->data() + (i * fItemSize));
      fSubFields[0]->DestroyValue(itemValue, true /* dtorOnly */);
   }
   vec->~vector();
   if (!dtorOnly)
      free(vec);
}

ROOT::Experimental::Detail::RFieldValue ROOT::Experimental::RVectorField::CaptureValue(void* where)
{
   return Detail::RFieldValue(true /* captureFlag */, this, where);
}

std::vector<ROOT::Experimental::Detail::RFieldValue>
ROOT::Experimental::RVectorField::SplitValue(const Detail::RFieldValue &value) const
{
   auto vec = static_cast<std::vector<char>*>(value.GetRawPtr());
   R__ASSERT((vec->size() % fItemSize) == 0);
   auto nItems = vec->size() / fItemSize;
   std::vector<Detail::RFieldValue> result;
   for (unsigned i = 0; i < nItems; ++i) {
      result.emplace_back(fSubFields[0]->CaptureValue(vec->data() + (i * fItemSize)));
   }
   return result;
}

void ROOT::Experimental::RVectorField::CommitCluster()
{
   fNWritten = 0;
}

void ROOT::Experimental::RVectorField::AcceptVisitor(Detail::RFieldVisitor &visitor) const
{
   visitor.VisitVectorField(*this);
}


//------------------------------------------------------------------------------


ROOT::Experimental::RField<std::vector<bool>>::RField(std::string_view name)
   : ROOT::Experimental::Detail::RFieldBase(name, "std::vector<bool>", ENTupleStructure::kCollection,
                                            false /* isSimple */)
{
   Attach(std::make_unique<RField<bool>>("bool"));
}

void ROOT::Experimental::RField<std::vector<bool>>::AppendImpl(const Detail::RFieldValue& value) {
   auto typedValue = value.Get<std::vector<bool>>();
   auto count = typedValue->size();
   for (unsigned i = 0; i < count; ++i) {
      bool bval = (*typedValue)[i];
      auto itemValue = fSubFields[0]->CaptureValue(&bval);
      fSubFields[0]->Append(itemValue);
   }
   Detail::RColumnElement<ClusterSize_t, EColumnType::kIndex> elemIndex(&fNWritten);
   fNWritten += count;
   fColumns[0]->Append(elemIndex);
}

void ROOT::Experimental::RField<std::vector<bool>>::ReadGlobalImpl(NTupleSize_t globalIndex, Detail::RFieldValue* value)
{
   auto typedValue = value->Get<std::vector<bool>>();

   ClusterSize_t nItems;
   RClusterIndex collectionStart;
   fPrincipalColumn->GetCollectionInfo(globalIndex, &collectionStart, &nItems);

   typedValue->resize(nItems);
   for (unsigned i = 0; i < nItems; ++i) {
      bool bval;
      auto itemValue = fSubFields[0]->GenerateValue(&bval);
      fSubFields[0]->Read(collectionStart + i, &itemValue);
      (*typedValue)[i] = bval;
   }
}

void ROOT::Experimental::RField<std::vector<bool>>::GenerateColumnsImpl()
{
   RColumnModel modelIndex(EColumnType::kIndex, true /* isSorted*/);
   fColumns.emplace_back(std::unique_ptr<Detail::RColumn>(
      Detail::RColumn::Create<ClusterSize_t, EColumnType::kIndex>(modelIndex, 0)));
   fPrincipalColumn = fColumns[0].get();
}

std::vector<ROOT::Experimental::Detail::RFieldValue>
ROOT::Experimental::RField<std::vector<bool>>::SplitValue(const Detail::RFieldValue& value) const
{
   const static bool trueValue = true;
   const static bool falseValue = false;

   auto typedValue = value.Get<std::vector<bool>>();
   auto count = typedValue->size();
   std::vector<Detail::RFieldValue> result;
   for (unsigned i = 0; i < count; ++i) {
      if ((*typedValue)[i])
         result.emplace_back(fSubFields[0]->CaptureValue(const_cast<bool *>(&trueValue)));
      else
         result.emplace_back(fSubFields[0]->CaptureValue(const_cast<bool *>(&falseValue)));
   }
   return result;
}


void ROOT::Experimental::RField<std::vector<bool>>::DestroyValue(const Detail::RFieldValue& value, bool dtorOnly)
{
   auto vec = static_cast<std::vector<bool>*>(value.GetRawPtr());
   vec->~vector();
   if (!dtorOnly)
      free(vec);
}

void ROOT::Experimental::RField<std::vector<bool>>::AcceptVisitor(Detail::RFieldVisitor &visitor) const
{
   visitor.VisitVectorBoolField(*this);
}


//------------------------------------------------------------------------------


ROOT::Experimental::RArrayField::RArrayField(
   std::string_view fieldName, std::unique_ptr<Detail::RFieldBase> itemField, std::size_t arrayLength)
   : ROOT::Experimental::Detail::RFieldBase(
      fieldName, "std::array<" + itemField->GetType() + "," + std::to_string(arrayLength) + ">",
      ENTupleStructure::kLeaf, false /* isSimple */, arrayLength)
   , fItemSize(itemField->GetValueSize()), fArrayLength(arrayLength)
{
   Attach(std::move(itemField));
}

std::unique_ptr<ROOT::Experimental::Detail::RFieldBase>
ROOT::Experimental::RArrayField::Clone(std::string_view newName) const
{
   auto newItemField = fSubFields[0]->Clone(fSubFields[0]->GetName());
   return std::make_unique<RArrayField>(newName, std::move(newItemField), fArrayLength);
}

void ROOT::Experimental::RArrayField::AppendImpl(const Detail::RFieldValue& value) {
   auto arrayPtr = value.Get<unsigned char>();
   for (unsigned i = 0; i < fArrayLength; ++i) {
      auto itemValue = fSubFields[0]->CaptureValue(arrayPtr + (i * fItemSize));
      fSubFields[0]->Append(itemValue);
   }
}

void ROOT::Experimental::RArrayField::ReadGlobalImpl(NTupleSize_t globalIndex, Detail::RFieldValue *value)
{
   auto arrayPtr = value->Get<unsigned char>();
   for (unsigned i = 0; i < fArrayLength; ++i) {
      auto itemValue = fSubFields[0]->GenerateValue(arrayPtr + (i * fItemSize));
      fSubFields[0]->Read(globalIndex * fArrayLength + i, &itemValue);
   }
}

void ROOT::Experimental::RArrayField::ReadInClusterImpl(const RClusterIndex &clusterIndex, Detail::RFieldValue *value)
{
   auto arrayPtr = value->Get<unsigned char>();
   for (unsigned i = 0; i < fArrayLength; ++i) {
      auto itemValue = fSubFields[0]->GenerateValue(arrayPtr + (i * fItemSize));
      fSubFields[0]->Read(RClusterIndex(clusterIndex.GetClusterId(), clusterIndex.GetIndex() * fArrayLength + i),
                          &itemValue);
   }
}

void ROOT::Experimental::RArrayField::GenerateColumnsImpl()
{
}

ROOT::Experimental::Detail::RFieldValue ROOT::Experimental::RArrayField::GenerateValue(void *where)
{
   auto arrayPtr = reinterpret_cast<unsigned char *>(where);
   for (unsigned i = 0; i < fArrayLength; ++i) {
      fSubFields[0]->GenerateValue(arrayPtr + (i * fItemSize));
   }
   return Detail::RFieldValue(true /* captureFlag */, this, where);
}

void ROOT::Experimental::RArrayField::DestroyValue(const Detail::RFieldValue& value, bool dtorOnly)
{
   auto arrayPtr = value.Get<unsigned char>();
   for (unsigned i = 0; i < fArrayLength; ++i) {
      auto itemValue = fSubFields[0]->CaptureValue(arrayPtr + (i * fItemSize));
      fSubFields[0]->DestroyValue(itemValue, true /* dtorOnly */);
   }
   if (!dtorOnly)
      free(arrayPtr);
}

ROOT::Experimental::Detail::RFieldValue ROOT::Experimental::RArrayField::CaptureValue(void *where)
{
   return Detail::RFieldValue(true /* captureFlag */, this, where);
}

std::vector<ROOT::Experimental::Detail::RFieldValue>
ROOT::Experimental::RArrayField::SplitValue(const Detail::RFieldValue &value) const
{
   auto arrayPtr = value.Get<unsigned char>();
   std::vector<Detail::RFieldValue> result;
   for (unsigned i = 0; i < fArrayLength; ++i) {
      auto itemValue = fSubFields[0]->CaptureValue(arrayPtr + (i * fItemSize));
      result.emplace_back(itemValue);
   }
   return result;
}

void ROOT::Experimental::RArrayField::AcceptVisitor(Detail::RFieldVisitor &visitor) const
{
   visitor.VisitArrayField(*this);
}

//------------------------------------------------------------------------------

#if __cplusplus >= 201703L
std::string ROOT::Experimental::RVariantField::GetTypeList(const std::vector<Detail::RFieldBase *> &itemFields)
{
   std::string result;
   for (size_t i = 0; i < itemFields.size(); ++i) {
      result += itemFields[i]->GetType() + ",";
   }
   R__ASSERT(!result.empty()); // there is always at least one variant
   result.pop_back(); // remove trailing comma
   return result;
}

ROOT::Experimental::RVariantField::RVariantField(
   std::string_view fieldName, const std::vector<Detail::RFieldBase *> &itemFields)
   : ROOT::Experimental::Detail::RFieldBase(fieldName,
      "std::variant<" + GetTypeList(itemFields) + ">", ENTupleStructure::kVariant, false /* isSimple */)
{
   auto nFields = itemFields.size();
   R__ASSERT(nFields > 0);
   fNWritten.resize(nFields, 0);
   for (unsigned int i = 0; i < nFields; ++i) {
      fMaxItemSize = std::max(fMaxItemSize, itemFields[i]->GetValueSize());
      fMaxAlignment = std::max(fMaxAlignment, itemFields[i]->GetAlignment());
      Attach(std::unique_ptr<Detail::RFieldBase>(itemFields[i]));
   }
   fTagOffset = (fMaxItemSize < fMaxAlignment) ? fMaxAlignment : fMaxItemSize;
}

std::unique_ptr<ROOT::Experimental::Detail::RFieldBase>
ROOT::Experimental::RVariantField::Clone(std::string_view newName) const
{
   auto nFields = fSubFields.size();
   std::vector<Detail::RFieldBase *> itemFields;
   for (unsigned i = 0; i < nFields; ++i) {
      // TODO(jblomer): use unique_ptr in RVariantField constructor
      itemFields.emplace_back(fSubFields[i]->Clone(fSubFields[i]->GetName()).release());
   }
   return std::make_unique<RVariantField>(newName, itemFields);
}

std::uint32_t ROOT::Experimental::RVariantField::GetTag(void *variantPtr) const
{
   auto index = *(reinterpret_cast<char *>(variantPtr) + fTagOffset);
   return (index < 0) ? 0 : index + 1;
}

void ROOT::Experimental::RVariantField::SetTag(void *variantPtr, std::uint32_t tag) const
{
   auto index = reinterpret_cast<char *>(variantPtr) + fTagOffset;
   *index = static_cast<char>(tag - 1);
}

void ROOT::Experimental::RVariantField::AppendImpl(const Detail::RFieldValue& value)
{
   auto tag = GetTag(value.GetRawPtr());
   auto index = 0;
   if (tag > 0) {
      auto itemValue = fSubFields[tag - 1]->CaptureValue(value.GetRawPtr());
      fSubFields[tag - 1]->Append(itemValue);
      index = fNWritten[tag - 1]++;
   }
   RColumnSwitch varSwitch(ClusterSize_t(index), tag);
   Detail::RColumnElement<RColumnSwitch, EColumnType::kSwitch> elemSwitch(&varSwitch);
   fColumns[0]->Append(elemSwitch);
}

void ROOT::Experimental::RVariantField::ReadGlobalImpl(NTupleSize_t globalIndex, Detail::RFieldValue *value)
{
   RClusterIndex variantIndex;
   std::uint32_t tag;
   fPrincipalColumn->GetSwitchInfo(globalIndex, &variantIndex, &tag);
   R__ASSERT(tag > 0); // TODO(jblomer): deal with invalid variants

   auto itemValue = fSubFields[tag - 1]->GenerateValue(value->GetRawPtr());
   fSubFields[tag - 1]->Read(variantIndex, &itemValue);
   SetTag(value->GetRawPtr(), tag);
}

void ROOT::Experimental::RVariantField::GenerateColumnsImpl()
{
   RColumnModel modelSwitch(EColumnType::kSwitch, false);
   fColumns.emplace_back(std::unique_ptr<Detail::RColumn>(
      Detail::RColumn::Create<RColumnSwitch, EColumnType::kSwitch>(modelSwitch, 0)));
   fPrincipalColumn = fColumns[0].get();
}

ROOT::Experimental::Detail::RFieldValue ROOT::Experimental::RVariantField::GenerateValue(void *where)
{
   memset(where, 0, GetValueSize());
   fSubFields[0]->GenerateValue(where);
   SetTag(where, 1);
   return Detail::RFieldValue(this, reinterpret_cast<unsigned char *>(where));
}

void ROOT::Experimental::RVariantField::DestroyValue(const Detail::RFieldValue& value, bool dtorOnly)
{
   auto variantPtr = value.GetRawPtr();
   auto tag = GetTag(variantPtr);
   if (tag > 0) {
      auto itemValue = fSubFields[tag - 1]->CaptureValue(variantPtr);
      fSubFields[tag - 1]->DestroyValue(itemValue, true /* dtorOnly */);
   }
   if (!dtorOnly)
      free(variantPtr);
}

ROOT::Experimental::Detail::RFieldValue ROOT::Experimental::RVariantField::CaptureValue(void *where)
{
   return Detail::RFieldValue(true /* captureFlag */, this, where);
}

size_t ROOT::Experimental::RVariantField::GetValueSize() const
{
   return fMaxItemSize + fMaxAlignment;  // TODO: fix for more than 255 items
}

void ROOT::Experimental::RVariantField::CommitCluster()
{
   std::fill(fNWritten.begin(), fNWritten.end(), 0);
}
#endif


//------------------------------------------------------------------------------


ROOT::Experimental::RCollectionField::RCollectionField(
   std::string_view name,
   std::shared_ptr<RCollectionNTuple> collectionNTuple,
   std::unique_ptr<RNTupleModel> collectionModel)
   : RFieldBase(name, ":Collection:", ENTupleStructure::kCollection, true /* isSimple */)
   , fCollectionNTuple(collectionNTuple)
{
   for (unsigned i = 0; i < collectionModel->GetFieldZero()->fSubFields.size(); ++i) {
      auto& subField = collectionModel->GetFieldZero()->fSubFields[i];
      Attach(std::move(subField));
   }
}


void ROOT::Experimental::RCollectionField::GenerateColumnsImpl()
{
   RColumnModel modelIndex(EColumnType::kIndex, true /* isSorted*/);
   fColumns.emplace_back(std::unique_ptr<Detail::RColumn>(
      Detail::RColumn::Create<ClusterSize_t, EColumnType::kIndex>(modelIndex, 0)));
   fPrincipalColumn = fColumns[0].get();
}


std::unique_ptr<ROOT::Experimental::Detail::RFieldBase>
ROOT::Experimental::RCollectionField::Clone(std::string_view /*newName*/) const
{
   // TODO(jblomer)
   return nullptr;
   //auto result = new RCollectionField(newName, fCollectionNTuple, RNTupleModel::Create());
   //for (auto& f : fSubFields) {
   //   // switch the name prefix for the new parent name
   //   std::string cloneName = std::string(newName) + f->GetName().substr(GetName().length());
   //   auto clone = f->Clone(cloneName);
   //   result->Attach(std::unique_ptr<RFieldBase>(clone));
   //}
   //return result;
}

void ROOT::Experimental::RCollectionField::CommitCluster() {
   *fCollectionNTuple->GetOffsetPtr() = 0;
}

