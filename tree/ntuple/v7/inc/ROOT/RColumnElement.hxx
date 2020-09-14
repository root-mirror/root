/// \file ROOT/RColumnElement.hxx
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

#ifndef ROOT7_RColumnElement
#define ROOT7_RColumnElement

#include <ROOT/RColumnModel.hxx>
#include <ROOT/RNTupleUtil.hxx>

#include <TError.h>

#include <cstring> // for memcpy
#include <cstdint>
#include <memory>
#include <type_traits>

namespace ROOT {
namespace Experimental {

namespace Detail {

// clang-format off
/**
\class ROOT::Experimental::Detail::RColumnElement
\ingroup NTuple
\brief A column element points either to the content of an RFieldValue or into a memory mapped page.

The content pointed to by fRawContent can be a single element or the first element of an array.
Usually the on-disk element should map bitwise to the in-memory element. Sometimes that's not the case
though, for instance on big endian platforms and for exotic physical columns like 8 bit float.

This class does not provide protection around the raw pointer, fRawContent has to be managed correctly
by the user of this class.
*/
// clang-format on
class RColumnElementBase {
protected:
   /// Points to valid C++ data, either a single value or an array of values
   void* fRawContent;
   /// Size of the C++ value pointed to by fRawContent (not necessarily equal to the on-disk element size)
   std::size_t fSize;

public:
   RColumnElementBase()
     : fRawContent(nullptr)
     , fSize(0)
   {}
   RColumnElementBase(void *rawContent, std::size_t size) : fRawContent(rawContent), fSize(size)
   {}
   RColumnElementBase(const RColumnElementBase &elemArray, std::size_t at)
     : fRawContent(static_cast<unsigned char *>(elemArray.fRawContent) + elemArray.fSize * at)
     , fSize(elemArray.fSize)
   {}
   RColumnElementBase(const RColumnElementBase& other) = default;
   RColumnElementBase(RColumnElementBase&& other) = default;
   RColumnElementBase& operator =(const RColumnElementBase& other) = delete;
   RColumnElementBase& operator =(RColumnElementBase&& other) = default;
   virtual ~RColumnElementBase() = default;

   static std::unique_ptr<RColumnElementBase> Generate(EColumnType type);

   /// Write one or multiple column elements into destination
   void WriteTo(void *destination, std::size_t count) const {
      std::memcpy(destination, fRawContent, fSize * count);
   }

   /// Set the column element or an array of elements from the memory location source
   void ReadFrom(void *source, std::size_t count) {
      std::memcpy(fRawContent, source, fSize * count);
   }

   /// Derived, typed classes tell whether the on-storage layout is bitwise identical to the memory layout
   virtual bool IsMappable() const { R__ASSERT(false); return false; }
   virtual std::size_t GetBitsOnStorage() const { R__ASSERT(false); return 0; }

   /// If the on-storage layout and the in-memory layout differ, packing creates an on-disk page from an in-memory page
   virtual void Pack(void *destination, void *source, std::size_t count) const
   {
      std::memcpy(destination, source, count);
   }

   /// If the on-storage layout and the in-memory layout differ, unpacking creates a memory page from an on-storage page
   virtual void Unpack(void *destination, void *source, std::size_t count) const
   {
      std::memcpy(destination, source, count);
   }

   void *GetRawContent() const { return fRawContent; }
   std::size_t GetSize() const { return fSize; }
};

/**
 * Pairs of C++ type and column type, like float and EColumnType::kReal32
 */
template <typename CppT, EColumnType ColumnT>
class RColumnElement : public RColumnElementBase {
public:
   explicit RColumnElement(CppT* value) : RColumnElementBase(value, sizeof(CppT))
   {
      // Do not allow this template to be instantiated unless there is a specialization. The assert needs to depend
      // on the template type or else the static_assert will always fire.
      static_assert(sizeof(CppT) != sizeof(CppT), "No column mapping for this C++ type");
   }
};


template <>
class RColumnElement<float, EColumnType::kReal32> : public RColumnElementBase {
public:
   static constexpr bool kIsMappable = true;
   static constexpr std::size_t kSize = sizeof(float);
   static constexpr std::size_t kBitsOnStorage = kSize * 8;
   explicit RColumnElement(float *value) : RColumnElementBase(value, kSize) {}
   bool IsMappable() const final { return kIsMappable; }
   std::size_t GetBitsOnStorage() const final { return kBitsOnStorage; }
};

template <>
class RColumnElement<double, EColumnType::kReal64> : public RColumnElementBase {
public:
   static constexpr bool kIsMappable = true;
   static constexpr std::size_t kSize = sizeof(double);
   static constexpr std::size_t kBitsOnStorage = kSize * 8;
   explicit RColumnElement(double *value) : RColumnElementBase(value, kSize) {}
   bool IsMappable() const final { return kIsMappable; }
   std::size_t GetBitsOnStorage() const final { return kBitsOnStorage; }
};

template <>
class RColumnElement<std::uint8_t, EColumnType::kByte> : public RColumnElementBase {
public:
   static constexpr bool kIsMappable = true;
   static constexpr std::size_t kSize = sizeof(std::uint8_t);
   static constexpr std::size_t kBitsOnStorage = kSize * 8;
   explicit RColumnElement(std::uint8_t *value) : RColumnElementBase(value, kSize) {}
   bool IsMappable() const final { return kIsMappable; }
   std::size_t GetBitsOnStorage() const final { return kBitsOnStorage; }
};

template <>
class RColumnElement<std::int32_t, EColumnType::kInt32> : public RColumnElementBase {
public:
   static constexpr bool kIsMappable = true;
   static constexpr std::size_t kSize = sizeof(std::int32_t);
   static constexpr std::size_t kBitsOnStorage = kSize * 8;
   explicit RColumnElement(std::int32_t *value) : RColumnElementBase(value, kSize) {}
   bool IsMappable() const final { return kIsMappable; }
   std::size_t GetBitsOnStorage() const final { return kBitsOnStorage; }
};

template <>
class RColumnElement<std::uint32_t, EColumnType::kInt32> : public RColumnElementBase {
public:
   static constexpr bool kIsMappable = true;
   static constexpr std::size_t kSize = sizeof(std::uint32_t);
   static constexpr std::size_t kBitsOnStorage = kSize * 8;
   explicit RColumnElement(std::uint32_t *value) : RColumnElementBase(value, kSize) {}
   bool IsMappable() const final { return kIsMappable; }
   std::size_t GetBitsOnStorage() const final { return kBitsOnStorage; }
};

template <>
class RColumnElement<std::int64_t, EColumnType::kInt64> : public RColumnElementBase {
public:
   static constexpr bool kIsMappable = true;
   static constexpr std::size_t kSize = sizeof(std::int64_t);
   static constexpr std::size_t kBitsOnStorage = kSize * 8;
   explicit RColumnElement(std::int64_t *value) : RColumnElementBase(value, kSize) {}
   bool IsMappable() const final { return kIsMappable; }
   std::size_t GetBitsOnStorage() const final { return kBitsOnStorage; }
};

template <>
class RColumnElement<std::uint64_t, EColumnType::kInt64> : public RColumnElementBase {
public:
   static constexpr bool kIsMappable = true;
   static constexpr std::size_t kSize = sizeof(std::uint64_t);
   static constexpr std::size_t kBitsOnStorage = kSize * 8;
   explicit RColumnElement(std::uint64_t *value) : RColumnElementBase(value, kSize) {}
   bool IsMappable() const final { return kIsMappable; }
   std::size_t GetBitsOnStorage() const final { return kBitsOnStorage; }
};

template <>
class RColumnElement<ClusterSize_t, EColumnType::kIndex> : public RColumnElementBase {
public:
   static constexpr bool kIsMappable = true;
   static constexpr std::size_t kSize = sizeof(ROOT::Experimental::ClusterSize_t);
   static constexpr std::size_t kBitsOnStorage = kSize * 8;
   explicit RColumnElement(ClusterSize_t *value) : RColumnElementBase(value, kSize) {}
   bool IsMappable() const final { return kIsMappable; }
   std::size_t GetBitsOnStorage() const final { return kBitsOnStorage; }
};

template <>
class RColumnElement<RColumnSwitch, EColumnType::kSwitch> : public RColumnElementBase {
public:
   static constexpr bool kIsMappable = true;
   static constexpr std::size_t kSize = sizeof(ROOT::Experimental::RColumnSwitch);
   static constexpr std::size_t kBitsOnStorage = kSize * 8;
   explicit RColumnElement(RColumnSwitch *value) : RColumnElementBase(value, kSize) {}
   bool IsMappable() const final { return kIsMappable; }
   std::size_t GetBitsOnStorage() const final { return kBitsOnStorage; }
};

template <>
class RColumnElement<char, EColumnType::kByte> : public RColumnElementBase {
public:
   static constexpr bool kIsMappable = true;
   static constexpr std::size_t kSize = sizeof(char);
   static constexpr std::size_t kBitsOnStorage = kSize * 8;
   explicit RColumnElement(char *value) : RColumnElementBase(value, kSize) {}
   bool IsMappable() const final { return kIsMappable; }
   std::size_t GetBitsOnStorage() const final { return kBitsOnStorage; }
};

template <>
class RColumnElement<bool, EColumnType::kBit> : public RColumnElementBase {
public:
   static constexpr bool kIsMappable = false;
   static constexpr std::size_t kSize = sizeof(bool);
   static constexpr std::size_t kBitsOnStorage = 1;
   explicit RColumnElement(bool *value) : RColumnElementBase(value, kSize) {}
   bool IsMappable() const final { return kIsMappable; }
   std::size_t GetBitsOnStorage() const final { return kBitsOnStorage; }

   void Pack(void *dst, void *src, std::size_t count) const final;
   void Unpack(void *dst, void *src, std::size_t count) const final;
};

} // namespace Detail
} // namespace Experimental
} // namespace ROOT

#endif
