/// \file ROOT/RColumnElement.hxx
/// \ingroup Forest ROOT7
/// \author Jakob Blomer <jblomer@cern.ch>
/// \date 2018-10-09
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2015, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RColumnElement
#define ROOT7_RColumnElement

#include <ROOT/RColumnModel.hxx>
#include <ROOT/RForestUtil.hxx>

#include <cstring> // for memcpy
#include <type_traits>

namespace ROOT {
namespace Experimental {

namespace Detail {

// clang-format off
/**
\class ROOT::Experimental::Detail::RColumnElement
\ingroup Forest
\brief A column element points either to the content of an RTreeValue or into a memory mapped page.

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
   const unsigned int fSize;

   /// Indicates that *fRawContent is bitwise identical to the physical column element
   const bool fIsMappable;

   virtual void DoSerialize(void* /* destination */, std::size_t /*count*/) const { }
   virtual void DoDeserialize(void* /* source */, std::size_t /*count*/) const { }

public:
   RColumnElementBase()
     : fRawContent(nullptr)
     , fSize(0)
     , fIsMappable(false)
   {}
   RColumnElementBase(void* rawContent, unsigned int size, bool isMappable)
     : fRawContent(rawContent)
     , fSize(size)
     , fIsMappable(isMappable)
   {}
   RColumnElementBase(const RColumnElementBase &elemArray, std::size_t at)
     : fRawContent(static_cast<unsigned char *>(elemArray.fRawContent) + elemArray.fSize * at)
     , fSize(elemArray.fSize)
     , fIsMappable(elemArray.fIsMappable)
   {}
   RColumnElementBase(const RColumnElementBase& other) = default;
   RColumnElementBase(RColumnElementBase&& other) = default;
   RColumnElementBase& operator =(const RColumnElementBase& other) = delete;
   RColumnElementBase& operator =(RColumnElementBase&& other) = default;
   virtual ~RColumnElementBase() = default;

   void Serialize(void *destination, std::size_t count) const {
     if (!fIsMappable) {
       DoSerialize(destination, count);
       return;
     }
     std::memcpy(destination, fRawContent, fSize * count);
   }

   void Deserialize(void *source, std::size_t count) {
     if (!fIsMappable) {
       DoDeserialize(source, count);
       return;
     }
     std::memcpy(fRawContent, source, fSize * count);
   }

   void SetRawContent(void* content) { fRawContent = content; }
   void* GetRawContent() const { return fRawContent; }
   unsigned int GetSize() const { return fSize; }
};

/**
 * Pairs of C++ type and column type, like float and EColumnType::kReal32
 */
template <typename CppT, EColumnType ColumnT>
class RColumnElement : public RColumnElementBase {
public:
   static constexpr bool kIsMappable = false;
   explicit RColumnElement(CppT* value) : RColumnElementBase(value, sizeof(CppT), kIsMappable)
   {
      static_assert(sizeof(CppT) != sizeof(CppT), "No column mapping for this C++ type");
   }
};

} // namespace Detail

} // namespace Experimental
} // namespace ROOT

template <>
class ROOT::Experimental::Detail::RColumnElement<float, ROOT::Experimental::EColumnType::kReal32>
   : public ROOT::Experimental::Detail::RColumnElementBase {
public:
   static constexpr bool kIsMappable = true;
   explicit RColumnElement(float* value) : RColumnElementBase(value, sizeof(float), kIsMappable) {}
};

template <>
class ROOT::Experimental::Detail::RColumnElement<
   ROOT::Experimental::TreeIndex_t, ROOT::Experimental::EColumnType::kIndex>
   : public ROOT::Experimental::Detail::RColumnElementBase {
public:
   static constexpr bool kIsMappable = true;
   explicit RColumnElement(ROOT::Experimental::TreeIndex_t* value)
      : RColumnElementBase(value, sizeof(ROOT::Experimental::TreeIndex_t), kIsMappable) {}
};

template <>
class ROOT::Experimental::Detail::RColumnElement<char, ROOT::Experimental::EColumnType::kByte>
   : public ROOT::Experimental::Detail::RColumnElementBase {
public:
   static constexpr bool kIsMappable = true;
   explicit RColumnElement(char* value) : RColumnElementBase(value, sizeof(char), kIsMappable) {}
};

#endif
