/// \file ROOT/RField.hxx
/// \ingroup NTuple ROOT7
/// \author Jakob Blomer <jblomer@cern.ch>
/// \date 2018-10-09
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RField
#define ROOT7_RField

#include <ROOT/RColumn.hxx>
#include <ROOT/RError.hxx>
#include <ROOT/RColumnElement.hxx>
#include <ROOT/RFieldValue.hxx>
#include <ROOT/RNTupleUtil.hxx>
#include <ROOT/RSpan.hxx>
#include <ROOT/RStringView.hxx>
#include <ROOT/RVec.hxx>
#include <ROOT/TypeTraits.hxx>

#include <TGenericClassInfo.h>

#include <algorithm>
#include <array>
#include <iostream>
#include <iterator>
#include <memory>
#include <string>
#include <type_traits>
#include <typeinfo>
#if __cplusplus >= 201703L
#include <variant>
#endif
#include <vector>
#include <utility>

class TClass;

namespace ROOT {
namespace Experimental {

class RCollectionField;
class RCollectionNTuple;
class REntry;
class RNTupleModel;

namespace Detail {

class RFieldFuse;
class RFieldVisitor;
class RPageStorage;

// clang-format off
/**
\class ROOT::Experimental::RFieldBase
\ingroup NTuple
\brief A field translates read and write calls from/to underlying columns to/from tree values

A field is a serializable C++ type or a container for a collection of sub fields. The RFieldBase and its
type-safe descendants provide the object to column mapper. They map C++ objects to primitive columns.  The
mapping is trivial for simple types such as 'double'. Complex types resolve to multiple primitive columns.
The field knows based on its type and the field name the type(s) and name(s) of the columns.
*/
// clang-format on
class RFieldBase {
   friend class ROOT::Experimental::Detail::RFieldFuse; // to connect the columns to a page storage
   friend class ROOT::Experimental::RCollectionField; // to change the field names when collections are attached

private:
   /// The field name relative to its parent field
   std::string fName;
   /// The C++ type captured by this field
   std::string fType;
   /// The role of this field in the data model structure
   ENTupleStructure fStructure;
   /// For fixed sized arrays, the array length
   std::size_t fNRepetitions;
   /// A field on a trivial type that maps as-is to a single column
   bool fIsSimple;

protected:
   /// Collections and classes own sub fields
   std::vector<std::unique_ptr<RFieldBase>> fSubFields;
   /// Sub fields point to their mother field
   RFieldBase* fParent;
   /// Points into fColumns.  All fields that have columns have a distinct main column. For simple fields
   /// (float, int, ...), the principal column corresponds to the field type. For collection fields expect std::array,
   /// the main column is the offset field.  Class fields have no column of their own.
   RColumn* fPrincipalColumn;
   /// The columns are connected either to a sink or to a source (not to both); they are owned by the field.
   std::vector<std::unique_ptr<RColumn>> fColumns;

   /// Creates the backing columns corresponsing to the field type and name
   virtual void GenerateColumnsImpl() = 0;

   /// Operations on values of complex types, e.g. ones that involve multiple columns or for which no direct
   /// column type exists.
   virtual void AppendImpl(const RFieldValue &value);
   virtual void ReadGlobalImpl(NTupleSize_t globalIndex, RFieldValue *value);
   virtual void ReadInClusterImpl(const RClusterIndex &clusterIndex, RFieldValue *value) {
      ReadGlobalImpl(fPrincipalColumn->GetGlobalIndex(clusterIndex), value);
   }

public:
   /// Iterates over the sub tree of fields in depth-first search order
   class RSchemaIterator {
   private:
      struct Position {
         Position() : fFieldPtr(nullptr), fIdxInParent(-1) { }
         Position(RFieldBase *fieldPtr, int idxInParent) : fFieldPtr(fieldPtr), fIdxInParent(idxInParent) { }
         RFieldBase *fFieldPtr;
         int fIdxInParent;
      };
      /// The stack of nodes visited when walking down the tree of fields
      std::vector<Position> fStack;
   public:
      using iterator = RSchemaIterator;
      using iterator_category = std::forward_iterator_tag;
      using value_type = RFieldBase;
      using difference_type = std::ptrdiff_t;
      using pointer = RFieldBase*;
      using reference = RFieldBase&;

      RSchemaIterator() { fStack.emplace_back(Position()); }
      RSchemaIterator(pointer val, int idxInParent) { fStack.emplace_back(Position(val, idxInParent)); }
      ~RSchemaIterator() {}
      /// Given that the iterator points to a valid field which is not the end iterator, go to the next field
      /// in depth-first search order
      void Advance();

      iterator  operator++(int) /* postfix */        { auto r = *this; Advance(); return r; }
      iterator& operator++()    /* prefix */         { Advance(); return *this; }
      reference operator* () const                   { return *fStack.back().fFieldPtr; }
      pointer   operator->() const                   { return fStack.back().fFieldPtr; }
      bool      operator==(const iterator& rh) const { return fStack.back().fFieldPtr == rh.fStack.back().fFieldPtr; }
      bool      operator!=(const iterator& rh) const { return fStack.back().fFieldPtr != rh.fStack.back().fFieldPtr; }
   };

   /// The constructor creates the underlying column objects and connects them to either a sink or a source.
   RFieldBase(std::string_view name, std::string_view type, ENTupleStructure structure, bool isSimple,
              std::size_t nRepetitions = 0);
   RFieldBase(const RFieldBase&) = delete;
   RFieldBase(RFieldBase&&) = default;
   RFieldBase& operator =(const RFieldBase&) = delete;
   RFieldBase& operator =(RFieldBase&&) = default;
   virtual ~RFieldBase();

   ///// Copies the field and its sub fields using a possibly new name and a new, unconnected set of columns
   virtual std::unique_ptr<RFieldBase> Clone(std::string_view newName) const = 0;

   /// Factory method to resurrect a field from the stored on-disk type information
   static RResult<std::unique_ptr<RFieldBase>> Create(const std::string &fieldName, const std::string &typeName);
   /// Check whether a given string is a valid field name
   static RResult<void> EnsureValidFieldName(std::string_view fieldName);

   /// Generates an object of the field type and allocates new initialized memory according to the type.
   RFieldValue GenerateValue();
   /// Generates a tree value in a given location of size at least GetValueSize(). Assumes that where has been
   /// allocated by malloc().
   virtual RFieldValue GenerateValue(void *where) = 0;
   /// Releases the resources acquired during GenerateValue (memory and constructor)
   /// This implementation works for simple types but needs to be overwritten for complex ones
   virtual void DestroyValue(const RFieldValue &value, bool dtorOnly = false);
   /// Creates a value from a memory location with an already constructed object
   virtual RFieldValue CaptureValue(void *where) = 0;
   /// Creates the list of direct child values given a value for this field.  E.g. a single value for the
   /// correct variant or all the elements of a collection.  The default implementation assumes no sub values
   /// and returns an empty vector.
   virtual std::vector<RFieldValue> SplitValue(const RFieldValue &value) const;
   /// The number of bytes taken by a value of the appropriate type
   virtual size_t GetValueSize() const = 0;
   /// For many types, the alignment requirement is equal to the size; otherwise override.
   virtual size_t GetAlignment() const { return GetValueSize(); }

   /// Write the given value into columns. The value object has to be of the same type as the field.
   void Append(const RFieldValue& value) {
      if (!fIsSimple) {
         AppendImpl(value);
         return;
      }
      fPrincipalColumn->Append(value.fMappedElement);
   }

   /// Populate a single value with data from the tree, which needs to be of the fitting type.
   /// Reading copies data into the memory wrapped by the ntuple value.
   void Read(NTupleSize_t globalIndex, RFieldValue *value) {
      if (!fIsSimple) {
         ReadGlobalImpl(globalIndex, value);
         return;
      }
      fPrincipalColumn->Read(globalIndex, &value->fMappedElement);
   }

   void Read(const RClusterIndex &clusterIndex, RFieldValue *value) {
      if (!fIsSimple) {
         ReadInClusterImpl(clusterIndex, value);
         return;
      }
      fPrincipalColumn->Read(clusterIndex, &value->fMappedElement);
   }

   /// Ensure that all received items are written from page buffers to the storage.
   void Flush() const;
   /// Perform housekeeping tasks for global to cluster-local index translation
   virtual void CommitCluster() {}

   /// Add a new subfield to the list of nested fields
   void Attach(std::unique_ptr<Detail::RFieldBase> child);

   std::string GetName() const { return fName; }
   std::string GetType() const { return fType; }
   ENTupleStructure GetStructure() const { return fStructure; }
   std::size_t GetNRepetitions() const { return fNRepetitions; }
   NTupleSize_t GetNElements() const { return fPrincipalColumn->GetNElements(); }
   const RFieldBase *GetParent() const { return fParent; }
   std::vector<const RFieldBase *> GetSubFields() const;
   bool IsSimple() const { return fIsSimple; }

   /// Indicates an evolution of the mapping scheme from C++ type to columns
   virtual RNTupleVersion GetFieldVersion() const { return RNTupleVersion(); }
   /// Indicates an evolution of the C++ type itself
   virtual RNTupleVersion GetTypeVersion() const { return RNTupleVersion(); }

   RSchemaIterator begin();
   RSchemaIterator end();

   virtual void AcceptVisitor(RFieldVisitor &visitor) const;
};

// clang-format off
/**
\class ROOT::Experimental::RFieldFuse
\ingroup NTuple
\brief A friend of RFieldBase responsible for connecting a field's columns to the physical page storage

Fields and their columns live in the void until connected to a physical page storage.  Only once connected, data
can be read or written.
*/
// clang-format on
class RFieldFuse {
public:
   static void Connect(DescriptorId_t fieldId, RPageStorage &pageStorage, RFieldBase &field);
   /// Connect the field columns and all sub field columns
   static void ConnectRecursively(DescriptorId_t fieldId, RPageSource &pageSource, RFieldBase &field);
};

} // namespace Detail



/// The container field for an ntuple model, which itself has no physical representation
class RFieldZero : public Detail::RFieldBase {
public:
   RFieldZero() : Detail::RFieldBase("", "", ENTupleStructure::kRecord, false /* isSimple */) { }
   std::unique_ptr<Detail::RFieldBase> Clone(std::string_view newName) const;

   void GenerateColumnsImpl() final {}
   using Detail::RFieldBase::GenerateValue;
   Detail::RFieldValue GenerateValue(void*) { return Detail::RFieldValue(); }
   Detail::RFieldValue CaptureValue(void*) final { return Detail::RFieldValue(); }
   size_t GetValueSize() const final { return 0; }

   /// Generates managed values for the top-level sub fields
   std::unique_ptr<REntry> GenerateEntry() const;
   void AcceptVisitor(Detail::RFieldVisitor &visitor) const final;
};

/// The field for a class with dictionary
class RClassField : public Detail::RFieldBase {
private:
   TClass* fClass;
   std::size_t fMaxAlignment = 1;

protected:
   void AppendImpl(const Detail::RFieldValue& value) final;
   void ReadGlobalImpl(NTupleSize_t globalIndex, Detail::RFieldValue *value) final;
   void ReadInClusterImpl(const RClusterIndex &clusterIndex, Detail::RFieldValue *value) final;

public:
   RClassField(std::string_view fieldName, std::string_view className);
   RClassField(RClassField&& other) = default;
   RClassField& operator =(RClassField&& other) = default;
   ~RClassField() = default;
   std::unique_ptr<Detail::RFieldBase> Clone(std::string_view newName) const final;

   void GenerateColumnsImpl() final;
   using Detail::RFieldBase::GenerateValue;
   Detail::RFieldValue GenerateValue(void* where) override;
   void DestroyValue(const Detail::RFieldValue& value, bool dtorOnly = false) final;
   Detail::RFieldValue CaptureValue(void *where) final;
   std::vector<Detail::RFieldValue> SplitValue(const Detail::RFieldValue &value) const final;
   size_t GetValueSize() const override;
   size_t GetAlignment() const final { return fMaxAlignment; }
   void AcceptVisitor(Detail::RFieldVisitor &visitor) const override;
};

/// The generic field for a (nested) std::vector<Type> except for std::vector<bool>
class RVectorField : public Detail::RFieldBase {
private:
   std::size_t fItemSize;
   ClusterSize_t fNWritten;

protected:
   void AppendImpl(const Detail::RFieldValue& value) final;
   void ReadGlobalImpl(NTupleSize_t globalIndex, Detail::RFieldValue *value) final;

public:
   RVectorField(std::string_view fieldName, std::unique_ptr<Detail::RFieldBase> itemField);
   RVectorField(RVectorField&& other) = default;
   RVectorField& operator =(RVectorField&& other) = default;
   ~RVectorField() = default;
   std::unique_ptr<Detail::RFieldBase> Clone(std::string_view newName) const final;

   void GenerateColumnsImpl() final;
   using Detail::RFieldBase::GenerateValue;
   Detail::RFieldValue GenerateValue(void* where) override;
   void DestroyValue(const Detail::RFieldValue& value, bool dtorOnly = false) final;
   Detail::RFieldValue CaptureValue(void *where) override;
   std::vector<Detail::RFieldValue> SplitValue(const Detail::RFieldValue &value) const final;
   size_t GetValueSize() const override { return sizeof(std::vector<char>); }
   size_t GetAlignment() const final { return std::alignment_of<std::vector<char>>(); }
   void CommitCluster() final;
   void AcceptVisitor(Detail::RFieldVisitor &visitor) const final;
   void GetCollectionInfo(NTupleSize_t globalIndex, RClusterIndex *collectionStart, ClusterSize_t *size) const {
      fPrincipalColumn->GetCollectionInfo(globalIndex, collectionStart, size);
   }
   void GetCollectionInfo(const RClusterIndex &clusterIndex, RClusterIndex *collectionStart, ClusterSize_t *size) const {
      fPrincipalColumn->GetCollectionInfo(clusterIndex, collectionStart, size);
   }
};


/// The generic field for fixed size arrays, which do not need an offset column
class RArrayField : public Detail::RFieldBase {
private:
   std::size_t fItemSize;
   std::size_t fArrayLength;

protected:
   void AppendImpl(const Detail::RFieldValue& value) final;
   void ReadGlobalImpl(NTupleSize_t globalIndex, Detail::RFieldValue *value) final;
   void ReadInClusterImpl(const RClusterIndex &clusterIndex, Detail::RFieldValue *value) final;

public:
   RArrayField(std::string_view fieldName, std::unique_ptr<Detail::RFieldBase> itemField, std::size_t arrayLength);
   RArrayField(RArrayField &&other) = default;
   RArrayField& operator =(RArrayField &&other) = default;
   ~RArrayField() = default;
   std::unique_ptr<Detail::RFieldBase> Clone(std::string_view newName) const final;

   void GenerateColumnsImpl() final;
   using Detail::RFieldBase::GenerateValue;
   Detail::RFieldValue GenerateValue(void *where) override;
   void DestroyValue(const Detail::RFieldValue &value, bool dtorOnly = false) final;
   Detail::RFieldValue CaptureValue(void *where) final;
   std::vector<Detail::RFieldValue> SplitValue(const Detail::RFieldValue &value) const final;
   size_t GetLength() const { return fArrayLength; }
   size_t GetValueSize() const final { return fItemSize * fArrayLength; }
   size_t GetAlignment() const final { return fSubFields[0]->GetAlignment(); }
   void AcceptVisitor(Detail::RFieldVisitor &visitor) const final;
};

#if __cplusplus >= 201703L
/// The generic field for std::variant types
class RVariantField : public Detail::RFieldBase {
private:
   size_t fMaxItemSize = 0;
   size_t fMaxAlignment = 1;
   /// In the std::variant memory layout, at which byte number is the index stored
   size_t fTagOffset = 0;
   std::vector<ClusterSize_t::ValueType> fNWritten;

   static std::string GetTypeList(const std::vector<Detail::RFieldBase *> &itemFields);
   /// Extracts the index from an std::variant and transforms it into the 1-based index used for the switch column
   std::uint32_t GetTag(void *variantPtr) const;
   void SetTag(void *variantPtr, std::uint32_t tag) const;

protected:
   void AppendImpl(const Detail::RFieldValue& value) final;
   void ReadGlobalImpl(NTupleSize_t globalIndex, Detail::RFieldValue *value) final;

public:
   // TODO(jblomer): use std::span in signature
   RVariantField(std::string_view fieldName, const std::vector<Detail::RFieldBase *> &itemFields);
   RVariantField(RVariantField &&other) = default;
   RVariantField& operator =(RVariantField &&other) = default;
   ~RVariantField() = default;
   std::unique_ptr<Detail::RFieldBase> Clone(std::string_view newName) const final;

   void GenerateColumnsImpl() final;
   using Detail::RFieldBase::GenerateValue;
   Detail::RFieldValue GenerateValue(void *where) override;
   void DestroyValue(const Detail::RFieldValue &value, bool dtorOnly = false) final;
   Detail::RFieldValue CaptureValue(void *where) final;
   size_t GetValueSize() const final;
   size_t GetAlignment() const final { return fMaxAlignment; }
   void CommitCluster() final;
};
#endif


/// Classes with dictionaries that can be inspected by TClass
template <typename T, typename=void>
class RField : public RClassField {
public:
   static std::string TypeName() { return ROOT::Internal::GetDemangledTypeName(typeid(T)); }
   RField(std::string_view name) : RClassField(name, TypeName()) {
      static_assert(std::is_class<T>::value, "no I/O support for this basic C++ type");
   }
   RField(RField&& other) = default;
   RField& operator =(RField&& other) = default;
   ~RField() = default;

   using Detail::RFieldBase::GenerateValue;
   template <typename... ArgsT>
   ROOT::Experimental::Detail::RFieldValue GenerateValue(void* where, ArgsT&&... args)
   {
      return Detail::RFieldValue(this, static_cast<T*>(where), std::forward<ArgsT>(args)...);
   }
   ROOT::Experimental::Detail::RFieldValue GenerateValue(void* where) final { return GenerateValue(where, T()); }
};


class RCollectionField : public ROOT::Experimental::Detail::RFieldBase {
private:
   /// Save the link to the collection ntuple in order to reset the offset counter when committing the cluster
   std::shared_ptr<RCollectionNTuple> fCollectionNTuple;
public:
   static std::string TypeName() { return ":RCollectionField:"; }
   RCollectionField(std::string_view name,
                    std::shared_ptr<RCollectionNTuple> collectionNTuple,
                    std::unique_ptr<RNTupleModel> collectionModel);
   RCollectionField(RCollectionField&& other) = default;
   RCollectionField& operator =(RCollectionField&& other) = default;
   ~RCollectionField() = default;
   std::unique_ptr<Detail::RFieldBase> Clone(std::string_view newName) const final;

   void GenerateColumnsImpl() final;

   using Detail::RFieldBase::GenerateValue;
   ROOT::Experimental::Detail::RFieldValue GenerateValue(void* where) final {
      return Detail::RFieldValue(
         Detail::RColumnElement<ClusterSize_t, EColumnType::kIndex>(static_cast<ClusterSize_t*>(where)),
         this, static_cast<ClusterSize_t*>(where));
   }
   Detail::RFieldValue CaptureValue(void* where) final {
      return Detail::RFieldValue(true /* captureFlag */,
         Detail::RColumnElement<ClusterSize_t, EColumnType::kIndex>(static_cast<ClusterSize_t*>(where)), this, where);
   }
   size_t GetValueSize() const final { return 0; }
   void CommitCluster() final;
};


/// Template specializations for concrete C++ types


template <>
class RField<ClusterSize_t> : public Detail::RFieldBase {
public:
   static std::string TypeName() { return "ROOT::Experimental::ClusterSize_t"; }
   explicit RField(std::string_view name)
     : Detail::RFieldBase(name, TypeName(), ENTupleStructure::kLeaf, true /* isSimple */) {}
   RField(RField&& other) = default;
   RField& operator =(RField&& other) = default;
   ~RField() = default;
   std::unique_ptr<Detail::RFieldBase> Clone(std::string_view newName) const final {
      return std::make_unique<RField>(newName);
   }

   void GenerateColumnsImpl() final;

   ClusterSize_t *Map(NTupleSize_t globalIndex) {
      return fPrincipalColumn->Map<ClusterSize_t, EColumnType::kIndex>(globalIndex);
   }
   ClusterSize_t *Map(const RClusterIndex &clusterIndex) {
      return fPrincipalColumn->Map<ClusterSize_t, EColumnType::kIndex>(clusterIndex);
   }

   using Detail::RFieldBase::GenerateValue;
   template <typename... ArgsT>
   ROOT::Experimental::Detail::RFieldValue GenerateValue(void* where, ArgsT&&... args)
   {
      return Detail::RFieldValue(
         Detail::RColumnElement<ClusterSize_t, EColumnType::kIndex>(static_cast<ClusterSize_t*>(where)),
         this, static_cast<ClusterSize_t*>(where), std::forward<ArgsT>(args)...);
   }
   ROOT::Experimental::Detail::RFieldValue GenerateValue(void* where) final { return GenerateValue(where, 0); }
   Detail::RFieldValue CaptureValue(void *where) final {
      return Detail::RFieldValue(true /* captureFlag */,
         Detail::RColumnElement<ClusterSize_t, EColumnType::kIndex>(static_cast<ClusterSize_t*>(where)), this, where);
   }
   size_t GetValueSize() const final { return sizeof(ClusterSize_t); }

   /// Special help for offset fields
   void GetCollectionInfo(NTupleSize_t globalIndex, RClusterIndex *collectionStart, ClusterSize_t *size) {
      fPrincipalColumn->GetCollectionInfo(globalIndex, collectionStart, size);
   }
   void GetCollectionInfo(const RClusterIndex &clusterIndex, RClusterIndex *collectionStart, ClusterSize_t *size) {
      fPrincipalColumn->GetCollectionInfo(clusterIndex, collectionStart, size);
   }
   void AcceptVisitor(Detail::RFieldVisitor &visitor) const final;
};


template <>
class RField<bool> : public Detail::RFieldBase {
public:
   static std::string TypeName() { return "bool"; }
   explicit RField(std::string_view name)
     : Detail::RFieldBase(name, TypeName(), ENTupleStructure::kLeaf, true /* isSimple */) {}
   RField(RField&& other) = default;
   RField& operator =(RField&& other) = default;
   ~RField() = default;
   std::unique_ptr<Detail::RFieldBase> Clone(std::string_view newName) const final {
      return std::make_unique<RField>(newName);
   }

   void GenerateColumnsImpl() final;

   bool *Map(NTupleSize_t globalIndex) {
      return fPrincipalColumn->Map<bool, EColumnType::kBit>(globalIndex);
   }
   bool *Map(const RClusterIndex &clusterIndex) {
      return fPrincipalColumn->Map<bool, EColumnType::kBit>(clusterIndex);
   }

   using Detail::RFieldBase::GenerateValue;
   template <typename... ArgsT>
   ROOT::Experimental::Detail::RFieldValue GenerateValue(void* where, ArgsT&&... args)
   {
      return Detail::RFieldValue(
         Detail::RColumnElement<bool, EColumnType::kBit>(static_cast<bool*>(where)),
         this, static_cast<bool*>(where), std::forward<ArgsT>(args)...);
   }
   ROOT::Experimental::Detail::RFieldValue GenerateValue(void* where) final { return GenerateValue(where, false); }
   Detail::RFieldValue CaptureValue(void *where) final {
      return Detail::RFieldValue(true /* captureFlag */,
         Detail::RColumnElement<bool, EColumnType::kBit>(static_cast<bool*>(where)), this, where);
   }
   size_t GetValueSize() const final { return sizeof(bool); }
   void AcceptVisitor(Detail::RFieldVisitor &visitor) const final;
};

template <>
class RField<float> : public Detail::RFieldBase {
public:
   static std::string TypeName() { return "float"; }
   explicit RField(std::string_view name)
     : Detail::RFieldBase(name, TypeName(), ENTupleStructure::kLeaf, true /* isSimple */) {}
   RField(RField&& other) = default;
   RField& operator =(RField&& other) = default;
   ~RField() = default;
   std::unique_ptr<Detail::RFieldBase> Clone(std::string_view newName) const final {
      return std::make_unique<RField>(newName);
   }

   void GenerateColumnsImpl() final;

   float *Map(NTupleSize_t globalIndex) {
      return fPrincipalColumn->Map<float, EColumnType::kReal32>(globalIndex);
   }
   float *Map(const RClusterIndex &clusterIndex) {
      return fPrincipalColumn->Map<float, EColumnType::kReal32>(clusterIndex);
   }

   using Detail::RFieldBase::GenerateValue;
   template <typename... ArgsT>
   ROOT::Experimental::Detail::RFieldValue GenerateValue(void* where, ArgsT&&... args)
   {
      return Detail::RFieldValue(
         Detail::RColumnElement<float, EColumnType::kReal32>(static_cast<float*>(where)),
         this, static_cast<float*>(where), std::forward<ArgsT>(args)...);
   }
   ROOT::Experimental::Detail::RFieldValue GenerateValue(void* where) final { return GenerateValue(where, 0.0); }
   Detail::RFieldValue CaptureValue(void *where) final {
      return Detail::RFieldValue(true /* captureFlag */,
         Detail::RColumnElement<float, EColumnType::kReal32>(static_cast<float*>(where)), this, where);
   }
   size_t GetValueSize() const final { return sizeof(float); }
   void AcceptVisitor(Detail::RFieldVisitor &visitor) const final;
};


template <>
class RField<double> : public Detail::RFieldBase {
public:
   static std::string TypeName() { return "double"; }
   explicit RField(std::string_view name)
     : Detail::RFieldBase(name, TypeName(), ENTupleStructure::kLeaf, true /* isSimple */) {}
   RField(RField&& other) = default;
   RField& operator =(RField&& other) = default;
   ~RField() = default;
   std::unique_ptr<Detail::RFieldBase> Clone(std::string_view newName) const final {
      return std::make_unique<RField>(newName);
   }

   void GenerateColumnsImpl() final;

   double *Map(NTupleSize_t globalIndex) {
      return fPrincipalColumn->Map<double, EColumnType::kReal64>(globalIndex);
   }
   double *Map(const RClusterIndex &clusterIndex) {
      return fPrincipalColumn->Map<double, EColumnType::kReal64>(clusterIndex);
   }

   using Detail::RFieldBase::GenerateValue;
   template <typename... ArgsT>
   ROOT::Experimental::Detail::RFieldValue GenerateValue(void* where, ArgsT&&... args)
   {
      return Detail::RFieldValue(
         Detail::RColumnElement<double, EColumnType::kReal64>(static_cast<double*>(where)),
         this, static_cast<double*>(where), std::forward<ArgsT>(args)...);
   }
   ROOT::Experimental::Detail::RFieldValue GenerateValue(void* where) final { return GenerateValue(where, 0.0); }
   Detail::RFieldValue CaptureValue(void *where) final {
      return Detail::RFieldValue(true /* captureFlag */,
         Detail::RColumnElement<double, EColumnType::kReal64>(static_cast<double*>(where)), this, where);
   }
   size_t GetValueSize() const final { return sizeof(double); }
   void AcceptVisitor(Detail::RFieldVisitor &visitor) const final;
};

template <>
class RField<std::uint8_t> : public Detail::RFieldBase {
public:
   static std::string TypeName() { return "std::uint8_t"; }
   explicit RField(std::string_view name)
     : Detail::RFieldBase(name, TypeName(), ENTupleStructure::kLeaf, true /* isSimple */) {}
   RField(RField&& other) = default;
   RField& operator =(RField&& other) = default;
   ~RField() = default;
   std::unique_ptr<Detail::RFieldBase> Clone(std::string_view newName) const final {
      return std::make_unique<RField>(newName);
   }

   void GenerateColumnsImpl() final;

   std::uint8_t *Map(NTupleSize_t globalIndex) {
      return fPrincipalColumn->Map<std::uint8_t, EColumnType::kByte>(globalIndex);
   }
   std::uint8_t *Map(const RClusterIndex &clusterIndex) {
      return fPrincipalColumn->Map<std::uint8_t, EColumnType::kByte>(clusterIndex);
   }

   using Detail::RFieldBase::GenerateValue;
   template <typename... ArgsT>
   ROOT::Experimental::Detail::RFieldValue GenerateValue(void *where, ArgsT&&... args)
   {
      return Detail::RFieldValue(
         Detail::RColumnElement<std::uint8_t, EColumnType::kByte>(static_cast<std::uint8_t*>(where)),
         this, static_cast<std::uint8_t*>(where), std::forward<ArgsT>(args)...);
   }
   ROOT::Experimental::Detail::RFieldValue GenerateValue(void *where) final { return GenerateValue(where, 0); }
   Detail::RFieldValue CaptureValue(void *where) final {
      return Detail::RFieldValue(true /* captureFlag */,
         Detail::RColumnElement<std::uint8_t, EColumnType::kByte>(static_cast<std::uint8_t*>(where)), this, where);
   }
   size_t GetValueSize() const final { return sizeof(std::uint8_t); }
   void AcceptVisitor(Detail::RFieldVisitor &visitor) const final;
};

template <>
class RField<std::int32_t> : public Detail::RFieldBase {
public:
   static std::string TypeName() { return "std::int32_t"; }
   explicit RField(std::string_view name)
     : Detail::RFieldBase(name, TypeName(), ENTupleStructure::kLeaf, true /* isSimple */) {}
   RField(RField&& other) = default;
   RField& operator =(RField&& other) = default;
   ~RField() = default;
   std::unique_ptr<Detail::RFieldBase> Clone(std::string_view newName) const final {
      return std::make_unique<RField>(newName);
   }

   void GenerateColumnsImpl() final;

   std::int32_t *Map(NTupleSize_t globalIndex) {
      return fPrincipalColumn->Map<std::int32_t, EColumnType::kInt32>(globalIndex);
   }
   std::int32_t *Map(const RClusterIndex &clusterIndex) {
      return fPrincipalColumn->Map<std::int32_t, EColumnType::kInt32>(clusterIndex);
   }

   using Detail::RFieldBase::GenerateValue;
   template <typename... ArgsT>
   ROOT::Experimental::Detail::RFieldValue GenerateValue(void* where, ArgsT&&... args)
   {
      return Detail::RFieldValue(
         Detail::RColumnElement<std::int32_t, EColumnType::kInt32>(static_cast<std::int32_t*>(where)),
         this, static_cast<std::int32_t*>(where), std::forward<ArgsT>(args)...);
   }
   ROOT::Experimental::Detail::RFieldValue GenerateValue(void* where) final { return GenerateValue(where, 0); }
   Detail::RFieldValue CaptureValue(void *where) final {
      return Detail::RFieldValue(true /* captureFlag */,
         Detail::RColumnElement<std::int32_t, EColumnType::kInt32>(static_cast<std::int32_t*>(where)), this, where);
   }
   size_t GetValueSize() const final { return sizeof(std::int32_t); }
   void AcceptVisitor(Detail::RFieldVisitor &visitor) const final;
};

template <>
class RField<std::uint32_t> : public Detail::RFieldBase {
public:
   static std::string TypeName() { return "std::uint32_t"; }
   explicit RField(std::string_view name)
     : Detail::RFieldBase(name, TypeName(), ENTupleStructure::kLeaf, true /* isSimple */) {}
   RField(RField&& other) = default;
   RField& operator =(RField&& other) = default;
   ~RField() = default;
   std::unique_ptr<Detail::RFieldBase> Clone(std::string_view newName) const final {
      return std::make_unique<RField>(newName);
   }

   void GenerateColumnsImpl() final;

   std::uint32_t *Map(NTupleSize_t globalIndex) {
      return fPrincipalColumn->Map<std::uint32_t, EColumnType::kInt32>(globalIndex);
   }
   std::uint32_t *Map(const RClusterIndex clusterIndex) {
      return fPrincipalColumn->Map<std::uint32_t, EColumnType::kInt32>(clusterIndex);
   }

   using Detail::RFieldBase::GenerateValue;
   template <typename... ArgsT>
   ROOT::Experimental::Detail::RFieldValue GenerateValue(void* where, ArgsT&&... args)
   {
      return Detail::RFieldValue(
         Detail::RColumnElement<std::uint32_t, EColumnType::kInt32>(static_cast<std::uint32_t*>(where)),
         this, static_cast<std::uint32_t*>(where), std::forward<ArgsT>(args)...);
   }
   ROOT::Experimental::Detail::RFieldValue GenerateValue(void* where) final { return GenerateValue(where, 0); }
   Detail::RFieldValue CaptureValue(void *where) final {
      return Detail::RFieldValue(true /* captureFlag */,
         Detail::RColumnElement<std::uint32_t, EColumnType::kInt32>(static_cast<std::uint32_t*>(where)), this, where);
   }
   size_t GetValueSize() const final { return sizeof(std::uint32_t); }
   void AcceptVisitor(Detail::RFieldVisitor &visitor) const final;
};

template <>
class RField<std::uint64_t> : public Detail::RFieldBase {
public:
   static std::string TypeName() { return "std::uint64_t"; }
   explicit RField(std::string_view name)
     : Detail::RFieldBase(name, TypeName(), ENTupleStructure::kLeaf, true /* isSimple */) {}
   RField(RField&& other) = default;
   RField& operator =(RField&& other) = default;
   ~RField() = default;
   std::unique_ptr<Detail::RFieldBase> Clone(std::string_view newName) const final {
      return std::make_unique<RField>(newName);
   }

   void GenerateColumnsImpl() final;

   std::uint64_t *Map(NTupleSize_t globalIndex) {
      return fPrincipalColumn->Map<std::uint64_t, EColumnType::kInt64>(globalIndex);
   }
   std::uint64_t *Map(const RClusterIndex &clusterIndex) {
      return fPrincipalColumn->Map<std::uint64_t, EColumnType::kInt64>(clusterIndex);
   }

   using Detail::RFieldBase::GenerateValue;
   template <typename... ArgsT>
   ROOT::Experimental::Detail::RFieldValue GenerateValue(void* where, ArgsT&&... args)
   {
      return Detail::RFieldValue(
         Detail::RColumnElement<std::uint64_t, EColumnType::kInt64>(static_cast<std::uint64_t*>(where)),
         this, static_cast<std::uint64_t*>(where), std::forward<ArgsT>(args)...);
   }
   ROOT::Experimental::Detail::RFieldValue GenerateValue(void* where) final { return GenerateValue(where, 0); }
   Detail::RFieldValue CaptureValue(void *where) final {
      return Detail::RFieldValue(true /* captureFlag */,
         Detail::RColumnElement<std::uint64_t, EColumnType::kInt64>(static_cast<std::uint64_t*>(where)), this, where);
   }
   size_t GetValueSize() const final { return sizeof(std::uint64_t); }
   void AcceptVisitor(Detail::RFieldVisitor &visitor) const final;
};


template <>
class RField<std::string> : public Detail::RFieldBase {
private:
   ClusterSize_t fIndex;
   Detail::RColumnElement<ClusterSize_t, EColumnType::kIndex> fElemIndex;

   void AppendImpl(const ROOT::Experimental::Detail::RFieldValue& value) final;
   void ReadGlobalImpl(ROOT::Experimental::NTupleSize_t globalIndex,
                       ROOT::Experimental::Detail::RFieldValue *value) final;

public:
   static std::string TypeName() { return "std::string"; }
   explicit RField(std::string_view name)
      : Detail::RFieldBase(name, TypeName(), ENTupleStructure::kLeaf, false /* isSimple */)
      , fIndex(0), fElemIndex(&fIndex) {}
   RField(RField&& other) = default;
   RField& operator =(RField&& other) = default;
   ~RField() = default;
   std::unique_ptr<Detail::RFieldBase> Clone(std::string_view newName) const final {
      return std::make_unique<RField>(newName);
   }

   void GenerateColumnsImpl() final;

   using Detail::RFieldBase::GenerateValue;
   template <typename... ArgsT>
   ROOT::Experimental::Detail::RFieldValue GenerateValue(void* where, ArgsT&&... args)
   {
      return Detail::RFieldValue(this, static_cast<std::string*>(where), std::forward<ArgsT>(args)...);
   }
   ROOT::Experimental::Detail::RFieldValue GenerateValue(void* where) final { return GenerateValue(where, ""); }
   void DestroyValue(const Detail::RFieldValue& value, bool dtorOnly = false) {
      auto str = value.Get<std::string>();
      str->~basic_string(); // TODO(jblomer) C++17 std::destroy_at
      if (!dtorOnly)
         free(str);
   }
   Detail::RFieldValue CaptureValue(void *where) {
      return Detail::RFieldValue(true /* captureFlag */, this, where);
   }
   size_t GetValueSize() const final { return sizeof(std::string); }
   size_t GetAlignment() const final { return std::alignment_of<std::string>(); }
   void CommitCluster() final;
   void AcceptVisitor(Detail::RFieldVisitor &visitor) const final;
};


template <typename ItemT, std::size_t N>
class RField<std::array<ItemT, N>> : public RArrayField {
   using ContainerT = typename std::array<ItemT, N>;
public:
   static std::string TypeName() {
      return "std::array<" + RField<ItemT>::TypeName() + "," + std::to_string(N) + ">";
   }
   explicit RField(std::string_view name)
      : RArrayField(name, std::make_unique<RField<ItemT>>(RField<ItemT>::TypeName()), N)
   {}
   RField(RField&& other) = default;
   RField& operator =(RField&& other) = default;
   ~RField() = default;

   using Detail::RFieldBase::GenerateValue;
   template <typename... ArgsT>
   ROOT::Experimental::Detail::RFieldValue GenerateValue(void *where, ArgsT&&... args)
   {
      return Detail::RFieldValue(this, static_cast<ContainerT*>(where), std::forward<ArgsT>(args)...);
   }
   ROOT::Experimental::Detail::RFieldValue GenerateValue(void *where) final {
      return GenerateValue(where, ContainerT());
   }
};


#if __cplusplus >= 201703L
template <typename... ItemTs>
class RField<std::variant<ItemTs...>> : public RVariantField {
   using ContainerT = typename std::variant<ItemTs...>;
private:
   template <typename HeadT, typename... TailTs>
   static std::string BuildItemTypes()
   {
      std::string result = RField<HeadT>::TypeName();
      if constexpr(sizeof...(TailTs) > 0)
         result += "," + BuildItemTypes<TailTs...>();
      return result;
   }

   template <typename HeadT, typename... TailTs>
   static std::vector<Detail::RFieldBase *> BuildItemFields(unsigned int index = 0)
   {
      std::vector<Detail::RFieldBase *> result;
      result.emplace_back(new RField<HeadT>("variant" + std::to_string(index)));
      if constexpr(sizeof...(TailTs) > 0) {
         auto tailFields = BuildItemFields<TailTs...>(index + 1);
         result.insert(result.end(), tailFields.begin(), tailFields.end());
      }
      return result;
   }

public:
   static std::string TypeName() { return "std::variant<" + BuildItemTypes<ItemTs...>() + ">"; }
   explicit RField(std::string_view name) : RVariantField(name, BuildItemFields<ItemTs...>()) {}
   RField(RField&& other) = default;
   RField& operator =(RField&& other) = default;
   ~RField() = default;

   using Detail::RFieldBase::GenerateValue;
   template <typename... ArgsT>
   ROOT::Experimental::Detail::RFieldValue GenerateValue(void *where, ArgsT&&... args)
   {
      return Detail::RFieldValue(this, static_cast<ContainerT*>(where), std::forward<ArgsT>(args)...);
   }
   ROOT::Experimental::Detail::RFieldValue GenerateValue(void *where) final {
      return GenerateValue(where, ContainerT());
   }
};
#endif

template <typename ItemT>
class RField<std::vector<ItemT>> : public RVectorField {
   using ContainerT = typename std::vector<ItemT>;
public:
   static std::string TypeName() { return "std::vector<" + RField<ItemT>::TypeName() + ">"; }
   explicit RField(std::string_view name)
      : RVectorField(name, std::make_unique<RField<ItemT>>(RField<ItemT>::TypeName()))
   {}
   RField(RField&& other) = default;
   RField& operator =(RField&& other) = default;
   ~RField() = default;

   using Detail::RFieldBase::GenerateValue;
   template <typename... ArgsT>
   ROOT::Experimental::Detail::RFieldValue GenerateValue(void* where, ArgsT&&... args)
   {
      return Detail::RFieldValue(this, static_cast<ContainerT*>(where), std::forward<ArgsT>(args)...);
   }
   ROOT::Experimental::Detail::RFieldValue GenerateValue(void* where) final {
      return GenerateValue(where, ContainerT());
   }
   Detail::RFieldValue CaptureValue(void *where) final {
      return Detail::RFieldValue(true /* captureFlag */, this, where);
   }
   size_t GetValueSize() const final { return sizeof(ContainerT); }
};

// std::vector<bool> is a template specialization and needs special treatment
template <>
class RField<std::vector<bool>> : public Detail::RFieldBase {
private:
   ClusterSize_t fNWritten{0};

protected:
   void AppendImpl(const Detail::RFieldValue& value) final;
   void ReadGlobalImpl(NTupleSize_t globalIndex, Detail::RFieldValue *value) final;
   void GenerateColumnsImpl() final;

public:
   static std::string TypeName() { return "std::vector<bool>"; }
   explicit RField(std::string_view name);
   RField(RField&& other) = default;
   RField& operator =(RField&& other) = default;
   ~RField() = default;
   std::unique_ptr<Detail::RFieldBase> Clone(std::string_view newName) const final {
      return std::make_unique<RField>(newName);
   }

   using Detail::RFieldBase::GenerateValue;
   template <typename... ArgsT>
   ROOT::Experimental::Detail::RFieldValue GenerateValue(void* where, ArgsT&&... args)
   {
      return Detail::RFieldValue(this, static_cast<std::vector<bool>*>(where), std::forward<ArgsT>(args)...);
   }
   ROOT::Experimental::Detail::RFieldValue GenerateValue(void* where) final {
      return GenerateValue(where, std::vector<bool>());
   }
   Detail::RFieldValue CaptureValue(void *where) final {
      return Detail::RFieldValue(true /* captureFlag */, this, where);
   }
   std::vector<Detail::RFieldValue> SplitValue(const Detail::RFieldValue &value) const final;
   void DestroyValue(const Detail::RFieldValue& value, bool dtorOnly = false) final;

   size_t GetValueSize() const final { return sizeof(std::vector<bool>); }
   size_t GetAlignment() const final { return std::alignment_of<std::vector<bool>>(); }
   void CommitCluster() final { fNWritten = 0; }
   void AcceptVisitor(Detail::RFieldVisitor &visitor) const final;
   void GetCollectionInfo(NTupleSize_t globalIndex, RClusterIndex *collectionStart, ClusterSize_t *size) const {
      fPrincipalColumn->GetCollectionInfo(globalIndex, collectionStart, size);
   }
   void GetCollectionInfo(const RClusterIndex &clusterIndex, RClusterIndex *collectionStart, ClusterSize_t *size) const {
      fPrincipalColumn->GetCollectionInfo(clusterIndex, collectionStart, size);
   }
};


/**
 * The RVec type has different layouts depending on the item type, therefore we cannot go with a generic
 * RVec implementation as we can with std::vector
 */
template <typename ItemT>
class RField<ROOT::VecOps::RVec<ItemT>> : public Detail::RFieldBase {
   using ContainerT = typename ROOT::VecOps::RVec<ItemT>;
private:
   size_t fItemSize;
   ClusterSize_t fNWritten;

protected:
   void AppendImpl(const Detail::RFieldValue& value) final {
      auto typedValue = value.Get<ContainerT>();
      auto count = typedValue->size();
      for (unsigned i = 0; i < count; ++i) {
         auto itemValue = fSubFields[0]->CaptureValue(&typedValue->data()[i]);
         fSubFields[0]->Append(itemValue);
      }
      Detail::RColumnElement<ClusterSize_t, EColumnType::kIndex> elemIndex(&fNWritten);
      fNWritten += count;
      fColumns[0]->Append(elemIndex);
   }
   void ReadGlobalImpl(NTupleSize_t globalIndex, Detail::RFieldValue *value) final {
      auto typedValue = value->Get<ContainerT>();
      ClusterSize_t nItems;
      RClusterIndex collectionStart;
      fPrincipalColumn->GetCollectionInfo(globalIndex, &collectionStart, &nItems);
      typedValue->resize(nItems);
      for (unsigned i = 0; i < nItems; ++i) {
         auto itemValue = fSubFields[0]->GenerateValue(&typedValue->data()[i]);
         fSubFields[0]->Read(collectionStart + i, &itemValue);
      }
   }

public:
   RField(std::string_view fieldName, std::unique_ptr<Detail::RFieldBase> itemField)
      : ROOT::Experimental::Detail::RFieldBase(
           fieldName, "ROOT::VecOps::RVec<" + itemField->GetType() + ">", ENTupleStructure::kCollection, false)
      , fItemSize(itemField->GetValueSize()), fNWritten(0)
   {
      Attach(std::move(itemField));
   }
   explicit RField(std::string_view name)
      : RField(name, std::make_unique<RField<ItemT>>(RField<ItemT>::TypeName()))
   {
   }
   RField(RField&& other) = default;
   RField& operator =(RField&& other) = default;
   ~RField() = default;
   std::unique_ptr<Detail::RFieldBase> Clone(std::string_view newName) const final {
      auto newItemField = fSubFields[0]->Clone(fSubFields[0]->GetName());
      return std::make_unique<RField<ROOT::VecOps::RVec<ItemT>>>(newName, std::move(newItemField));
   }

   void GenerateColumnsImpl() final {
      RColumnModel modelIndex(EColumnType::kIndex, true /* isSorted*/);
      fColumns.emplace_back(std::unique_ptr<Detail::RColumn>(
         Detail::RColumn::Create<ClusterSize_t, EColumnType::kIndex>(modelIndex, 0)));
      fPrincipalColumn = fColumns[0].get();
   }
   void DestroyValue(const Detail::RFieldValue& value, bool dtorOnly = false) final {
      auto vec = reinterpret_cast<ContainerT*>(value.GetRawPtr());
      auto nItems = vec->size();
      for (unsigned i = 0; i < nItems; ++i) {
         auto itemValue = fSubFields[0]->CaptureValue(vec->data() + (i * fItemSize));
         fSubFields[0]->DestroyValue(itemValue, true /* dtorOnly */);
      }
      vec->~RVec();
      if (!dtorOnly)
         free(vec);
   }
   void CommitCluster() final { fNWritten = 0; }

   static std::string TypeName() { return "ROOT::VecOps::RVec<" + RField<ItemT>::TypeName() + ">"; }

   using Detail::RFieldBase::GenerateValue;
   template <typename... ArgsT>
   ROOT::Experimental::Detail::RFieldValue GenerateValue(void* where, ArgsT&&... args)
   {
      return Detail::RFieldValue(this, static_cast<ContainerT*>(where), std::forward<ArgsT>(args)...);
   }
   ROOT::Experimental::Detail::RFieldValue GenerateValue(void* where) final {
      return GenerateValue(where, ContainerT());
   }
   Detail::RFieldValue CaptureValue(void *where) final {
      return Detail::RFieldValue(true /* captureFlag */, this, static_cast<ContainerT*>(where));
   }
   size_t GetValueSize() const final { return sizeof(ContainerT); }
   size_t GetAlignment() const final { return std::alignment_of<ContainerT>(); }
};

/**
 * RVec<bool> needs special treatment due to std::vector<bool> sepcialization
 */
template <>
class RField<ROOT::VecOps::RVec<bool>> : public Detail::RFieldBase {
   using ContainerT = typename ROOT::VecOps::RVec<bool>;
private:
   ClusterSize_t fNWritten{0};

protected:
   void AppendImpl(const Detail::RFieldValue& value) final {
      auto typedValue = value.Get<ContainerT>();
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
   void ReadGlobalImpl(NTupleSize_t globalIndex, Detail::RFieldValue *value) final {
      auto typedValue = value->Get<ContainerT>();
      ClusterSize_t nItems;
      RClusterIndex collectionStart;
      fPrincipalColumn->GetCollectionInfo(globalIndex, &collectionStart, &nItems);
      typedValue->resize(nItems);
      for (unsigned i = 0; i < nItems; ++i) {
         bool bval = (*typedValue)[i];
         auto itemValue = fSubFields[0]->GenerateValue(&bval);
         fSubFields[0]->Read(collectionStart + i, &itemValue);
         (*typedValue)[i] = bval;
      }
   }

public:
   RField(std::string_view name)
      : ROOT::Experimental::Detail::RFieldBase(name, "ROOT::VecOps::RVec<bool>", ENTupleStructure::kCollection, false)
   {
      Attach(std::make_unique<RField<bool>>("bool"));
   }
   RField(RField&& other) = default;
   RField& operator =(RField&& other) = default;
   ~RField() = default;
   std::unique_ptr<Detail::RFieldBase> Clone(std::string_view newName) const final {
      return std::make_unique<RField<ROOT::VecOps::RVec<bool>>>(newName);
   }

   void GenerateColumnsImpl() final {
      RColumnModel modelIndex(EColumnType::kIndex, true /* isSorted*/);
      fColumns.emplace_back(std::unique_ptr<Detail::RColumn>(
         Detail::RColumn::Create<ClusterSize_t, EColumnType::kIndex>(modelIndex, 0)));
      fPrincipalColumn = fColumns[0].get();
   }
   void DestroyValue(const Detail::RFieldValue& value, bool dtorOnly = false) final {
      auto vec = reinterpret_cast<ContainerT*>(value.GetRawPtr());
      vec->~RVec();
      if (!dtorOnly)
         free(vec);
   }
   void CommitCluster() final { fNWritten = 0; }

   static std::string TypeName() { return "ROOT::VecOps::RVec<bool>"; }

   using Detail::RFieldBase::GenerateValue;
   template <typename... ArgsT>
   ROOT::Experimental::Detail::RFieldValue GenerateValue(void *where, ArgsT&&... args)
   {
      return Detail::RFieldValue(this, static_cast<ContainerT*>(where), std::forward<ArgsT>(args)...);
   }
   ROOT::Experimental::Detail::RFieldValue GenerateValue(void *where) final {
      return GenerateValue(where, ContainerT());
   }
   Detail::RFieldValue CaptureValue(void *where) final {
      return Detail::RFieldValue(true /* captureFlag */, this, static_cast<ContainerT*>(where));
   }
   size_t GetValueSize() const final { return sizeof(ContainerT); }
   size_t GetAlignment() const final { return std::alignment_of<ContainerT>(); }
};

} // namespace Experimental
} // namespace ROOT

#endif
