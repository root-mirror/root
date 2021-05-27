// Author: Enrico Guiraud, Enric Tejedor, Danilo Piparo CERN  04/2021

/*************************************************************************
 * Copyright (C) 1995-2021, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/*
 * The implementation of RVec is based on LLVM's SmallVector.
 *
 * Compared to LLVM's SmallVectors, the ability to adopt an existing memory buffer has been added,
 * as well as other minor patches for backward compatibility with
 * the previous implementation of RVec (which did not have a small buffer
 * optimization and was based on std::vector with a custom allocator) and
 * for compatibility with cppyy.
 * RVec also adds a vectorized operator[] and many "numpy-like" helper functions.
 */

#ifndef ROOT_RVEC
#define ROOT_RVEC

#if __cplusplus > 201402L
#define R__RVEC_NODISCARD [[nodiscard]]
#else
#define R__RVEC_NODISCARD
#endif

#ifdef _WIN32
   #ifndef M_PI
      #ifndef _USE_MATH_DEFINES
         #define _USE_MATH_DEFINES
      #endif
      #include <math.h>
      #undef _USE_MATH_DEFINES
   #endif
   #define _VECOPS_USE_EXTERN_TEMPLATES false
#else
   #define _VECOPS_USE_EXTERN_TEMPLATES true
#endif

#include <ROOT/RAdoptAllocator.hxx>
#include <ROOT/RStringView.hxx>
#include <TError.h> // R__ASSERT
#include <ROOT/TypeTraits.hxx>

#include <algorithm>
#include <cmath>
#include <cassert>
#include <cstring>
#include <numeric> // for inner_product
#include <new>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <vector>
#include <utility>
#include <tuple>

#ifdef R__HAS_VDT
#include <vdt/vdtMath.h>
#endif


namespace ROOT {

namespace VecOps {
template<typename T>
class RVec;
}

namespace Detail {
namespace VecOps {

template<typename T>
using RVec = ROOT::VecOps::RVec<T>;

template <typename... T>
std::size_t GetVectorsSize(std::string_view id, const RVec<T> &... vs)
{
   constexpr const auto nArgs = sizeof...(T);
   const std::size_t sizes[] = {vs.size()...};
   if (nArgs > 1) {
      for (auto i = 1UL; i < nArgs; i++) {
         if (sizes[0] == sizes[i])
            continue;
         std::string msg(id);
         msg += ": input RVec instances have different lengths!";
         throw std::runtime_error(msg);
      }
   }
   return sizes[0];
}

template <typename F, typename... T>
auto MapImpl(F &&f, const RVec<T> &... vs) -> RVec<decltype(f(vs[0]...))>
{
   const auto size = GetVectorsSize("Map", vs...);
   RVec<decltype(f(vs[0]...))> ret(size);

   for (auto i = 0UL; i < size; i++)
      ret[i] = f(vs[i]...);

   return ret;
}

template <typename Tuple_t, std::size_t... Is>
auto MapFromTuple(Tuple_t &&t, std::index_sequence<Is...>)
   -> decltype(MapImpl(std::get<std::tuple_size<Tuple_t>::value - 1>(t), std::get<Is>(t)...))
{
   constexpr const auto tupleSizeM1 = std::tuple_size<Tuple_t>::value - 1;
   return MapImpl(std::get<tupleSizeM1>(t), std::get<Is>(t)...);
}

} // namespace VecOps
} // namespace Detail

namespace Internal {
namespace VecOps {
/// Return the next power of two (in 64-bits) that is strictly greater than A.
/// Return zero on overflow.
inline uint64_t NextPowerOf2(uint64_t A)
{
   A |= (A >> 1);
   A |= (A >> 2);
   A |= (A >> 4);
   A |= (A >> 8);
   A |= (A >> 16);
   A |= (A >> 32);
   return A + 1;
}

/// This is all the stuff common to all SmallVectors.
///
/// The template parameter specifies the type which should be used to hold the
/// Size and Capacity of the SmallVector, so it can be adjusted.
/// Using 32 bit size is desirable to shrink the size of the SmallVector.
/// Using 64 bit size is desirable for cases like SmallVector<char>, where a
/// 32 bit size would limit the vector to ~4GB. SmallVectors are used for
/// buffering bitcode output - which can exceed 4GB.
template <class Size_T>
class SmallVectorBase {
protected:
   void *BeginX;
   Size_T fSize = 0, Capacity;
   bool Owns = true;

   /// The maximum value of the Size_T used.
   static constexpr size_t SizeTypeMax() { return std::numeric_limits<Size_T>::max(); }

   SmallVectorBase() = delete;
   SmallVectorBase(void *FirstEl, size_t TotalCapacity) : BeginX(FirstEl), Capacity(TotalCapacity) {}

   /// This is an implementation of the grow() method which only works
   /// on POD-like data types and is out of line to reduce code duplication.
   /// This function will report a fatal error if it cannot increase capacity.
   void grow_pod(void *FirstEl, size_t MinSize, size_t TSize);

   /// Report that MinSize doesn't fit into this vector's size type. Throws
   /// std::length_error or calls report_fatal_error.
   static void report_size_overflow(size_t MinSize);
   /// Report that this vector is already at maximum capacity. Throws
   /// std::length_error or calls report_fatal_error.
   static void report_at_maximum_capacity();

public:
   size_t size() const { return fSize; }
   size_t capacity() const noexcept { return Capacity; }

   R__RVEC_NODISCARD bool empty() const { return !fSize; }

   /// Set the array size to \p N, which the current array must have enough
   /// capacity for.
   ///
   /// This does not construct or destroy any elements in the vector.
   ///
   /// Clients can use this in conjunction with capacity() to write past the end
   /// of the buffer when they know that more elements are available, and only
   /// update the size later. This avoids the cost of value initializing elements
   /// which will only be overwritten.
   void set_size(size_t N)
   {
      assert(N <= capacity());
      fSize = N;
   }

   // LLVM SmallVector does not have a shrink_to_fit method
   // it's technically ok to do nothing, but assuming no one uses this method for RVec anyway, I'd rather deprecate it
   R__DEPRECATED(6, 28, "This method will be removed.")
   void shrink_to_fit() { }
};

template <class T>
using SmallVectorSizeType = typename std::conditional<sizeof(T) < 4 && sizeof(void *) >= 8, uint64_t, uint32_t>::type;

/// Figure out the offset of the first element.
template <class T, typename = void>
struct SmallVectorAlignmentAndSize {
   alignas(SmallVectorBase<SmallVectorSizeType<T>>) char Base[sizeof(SmallVectorBase<SmallVectorSizeType<T>>)];
   alignas(T) char FirstEl[sizeof(T)];
};

/// This is the part of SmallVectorTemplateBase which does not depend on whether
/// the type T is a POD. The extra dummy template argument is used by ArrayRef
/// to avoid unnecessarily requiring T to be complete.
template <typename T, typename = void>
class SmallVectorTemplateCommon : public SmallVectorBase<SmallVectorSizeType<T>> {
   using Base = SmallVectorBase<SmallVectorSizeType<T>>;

   /// Find the address of the first element.  For this pointer math to be valid
   /// with small-size of 0 for T with lots of alignment, it's important that
   /// SmallVectorStorage is properly-aligned even for small-size of 0.
   void *getFirstEl() const
   {
      return const_cast<void *>(reinterpret_cast<const void *>(reinterpret_cast<const char *>(this) +
                                                               offsetof(SmallVectorAlignmentAndSize<T>, FirstEl)));
   }
   // Space after 'FirstEl' is clobbered, do not add any instance vars after it.

protected:
   SmallVectorTemplateCommon(size_t Size) : Base(getFirstEl(), Size) {}

   void grow_pod(size_t MinSize, size_t TSize) { Base::grow_pod(getFirstEl(), MinSize, TSize); }

   /// Return true if this is a smallvector which has not had dynamic
   /// memory allocated for it.
   bool isSmall() const { return this->BeginX == getFirstEl(); }

   /// Put this vector in a state of being small.
   void resetToSmall()
   {
      this->BeginX = getFirstEl();
      // from the original LLVM implementation:
      // FIXME: Setting Capacity to 0 is suspect.
      this->fSize = this->Capacity = 0;
   }

public:
   using size_type = size_t;
   using difference_type = ptrdiff_t;
   using value_type = T;
   using iterator = T *;
   using const_iterator = const T *;

   using const_reverse_iterator = std::reverse_iterator<const_iterator>;
   using reverse_iterator = std::reverse_iterator<iterator>;

   using reference = T &;
   using const_reference = const T &;
   using pointer = T *;
   using const_pointer = const T *;

   using Base::capacity;
   using Base::empty;
   using Base::size;

   // forward iterator creation methods.
   iterator begin() noexcept { return (iterator)this->BeginX; }
   const_iterator begin() const noexcept { return (const_iterator)this->BeginX; }
   const_iterator cbegin() const noexcept { return (const_iterator)this->BeginX; }
   iterator end() noexcept { return begin() + size(); }
   const_iterator end() const noexcept { return begin() + size(); }
   const_iterator cend() const noexcept { return begin() + size(); }

   // reverse iterator creation methods.
   reverse_iterator rbegin() noexcept { return reverse_iterator(end()); }
   const_reverse_iterator rbegin() const noexcept { return const_reverse_iterator(end()); }
   const_reverse_iterator crbegin() const noexcept { return const_reverse_iterator(end()); }
   reverse_iterator rend() noexcept { return reverse_iterator(begin()); }
   const_reverse_iterator rend() const noexcept { return const_reverse_iterator(begin()); }
   const_reverse_iterator crend() const noexcept { return const_reverse_iterator(begin()); }

   size_type size_in_bytes() const { return size() * sizeof(T); }
   size_type max_size() const noexcept { return std::min(this->SizeTypeMax(), size_type(-1) / sizeof(T)); }

   size_t capacity_in_bytes() const { return capacity() * sizeof(T); }

   /// Return a pointer to the vector's buffer, even if empty().
   pointer data() noexcept { return pointer(begin()); }
   /// Return a pointer to the vector's buffer, even if empty().
   const_pointer data() const noexcept { return const_pointer(begin()); }

   reference front()
   {
      assert(!empty());
      return begin()[0];
   }
   const_reference front() const
   {
      assert(!empty());
      return begin()[0];
   }

   reference back()
   {
      assert(!empty());
      return end()[-1];
   }
   const_reference back() const
   {
      assert(!empty());
      return end()[-1];
   }
};

/// SmallVectorTemplateBase<TriviallyCopyable = false> - This is where we put
/// method implementations that are designed to work with non-trivial T's.
///
/// We approximate is_trivially_copyable with trivial move/copy construction and
/// trivial destruction. While the standard doesn't specify that you're allowed
/// copy these types with memcpy, there is no way for the type to observe this.
/// This catches the important case of std::pair<POD, POD>, which is not
/// trivially assignable.
template <typename T, bool = (std::is_trivially_copy_constructible<T>::value) &&
                             (std::is_trivially_move_constructible<T>::value) &&
                             std::is_trivially_destructible<T>::value>
class SmallVectorTemplateBase : public SmallVectorTemplateCommon<T> {
protected:
   SmallVectorTemplateBase(size_t Size) : SmallVectorTemplateCommon<T>(Size) {}

   static void destroy_range(T *S, T *E)
   {
      while (S != E) {
         --E;
         E->~T();
      }
   }

   /// Move the range [I, E) into the uninitialized memory starting with "Dest",
   /// constructing elements as needed.
   template <typename It1, typename It2>
   static void uninitialized_move(It1 I, It1 E, It2 Dest)
   {
      std::uninitialized_copy(std::make_move_iterator(I), std::make_move_iterator(E), Dest);
   }

   /// Copy the range [I, E) onto the uninitialized memory starting with "Dest",
   /// constructing elements as needed.
   template <typename It1, typename It2>
   static void uninitialized_copy(It1 I, It1 E, It2 Dest)
   {
      std::uninitialized_copy(I, E, Dest);
   }

   /// Grow the allocated memory (without initializing new elements), doubling
   /// the size of the allocated memory. Guarantees space for at least one more
   /// element, or MinSize more elements if specified.
   void grow(size_t MinSize = 0);

public:
   void push_back(const T &Elt)
   {
      if (R__unlikely(this->size() >= this->capacity()))
         this->grow();
      ::new ((void *)this->end()) T(Elt);
      this->set_size(this->size() + 1);
   }

   void push_back(T &&Elt)
   {
      if (R__unlikely(this->size() >= this->capacity()))
         this->grow();
      ::new ((void *)this->end()) T(::std::move(Elt));
      this->set_size(this->size() + 1);
   }

   void pop_back()
   {
      this->set_size(this->size() - 1);
      this->end()->~T();
   }
};

// Define this out-of-line to dissuade the C++ compiler from inlining it.
template <typename T, bool TriviallyCopyable>
void SmallVectorTemplateBase<T, TriviallyCopyable>::grow(size_t MinSize)
{
   // Ensure we can fit the new capacity.
   // This is only going to be applicable when the capacity is 32 bit.
   if (MinSize > this->SizeTypeMax())
      this->report_size_overflow(MinSize);

   // Ensure we can meet the guarantee of space for at least one more element.
   // The above check alone will not catch the case where grow is called with a
   // default MinSize of 0, but the current capacity cannot be increased.
   // This is only going to be applicable when the capacity is 32 bit.
   if (this->capacity() == this->SizeTypeMax())
      this->report_at_maximum_capacity();

   // Always grow, even from zero.
   size_t NewCapacity = size_t(NextPowerOf2(this->capacity() + 2));
   NewCapacity = std::min(std::max(NewCapacity, MinSize), this->SizeTypeMax());
   T *NewElts = static_cast<T *>(malloc(NewCapacity * sizeof(T)));

   // Move the elements over.
   this->uninitialized_move(this->begin(), this->end(), NewElts);

   if (this->Owns) {
      // Destroy the original elements.
      destroy_range(this->begin(), this->end());

      // If this wasn't grown from the inline copy, deallocate the old space.
      if (!this->isSmall())
         free(this->begin());
   }

   this->BeginX = NewElts;
   this->Capacity = NewCapacity;
   this->Owns = true;
}

/// SmallVectorTemplateBase<TriviallyCopyable = true> - This is where we put
/// method implementations that are designed to work with trivially copyable
/// T's. This allows using memcpy in place of copy/move construction and
/// skipping destruction.
template <typename T>
class SmallVectorTemplateBase<T, true> : public SmallVectorTemplateCommon<T> {
   using SuperClass = SmallVectorTemplateCommon<T>;

protected:
   SmallVectorTemplateBase(size_t Size) : SmallVectorTemplateCommon<T>(Size) {}

   // No need to do a destroy loop for POD's.
   static void destroy_range(T *, T *) {}

   /// Move the range [I, E) onto the uninitialized memory
   /// starting with "Dest", constructing elements into it as needed.
   template <typename It1, typename It2>
   static void uninitialized_move(It1 I, It1 E, It2 Dest)
   {
      // Just do a copy.
      uninitialized_copy(I, E, Dest);
   }

   /// Copy the range [I, E) onto the uninitialized memory
   /// starting with "Dest", constructing elements into it as needed.
   template <typename It1, typename It2>
   static void uninitialized_copy(It1 I, It1 E, It2 Dest)
   {
      // Arbitrary iterator types; just use the basic implementation.
      std::uninitialized_copy(I, E, Dest);
   }

   /// Copy the range [I, E) onto the uninitialized memory
   /// starting with "Dest", constructing elements into it as needed.
   template <typename T1, typename T2>
   static void uninitialized_copy(
      T1 *I, T1 *E, T2 *Dest,
      typename std::enable_if<std::is_same<typename std::remove_const<T1>::type, T2>::value>::type * = nullptr)
   {
      // Use memcpy for PODs iterated by pointers (which includes SmallVector
      // iterators): std::uninitialized_copy optimizes to memmove, but we can
      // use memcpy here. Note that I and E are iterators and thus might be
      // invalid for memcpy if they are equal.
      if (I != E)
         memcpy(reinterpret_cast<void *>(Dest), I, (E - I) * sizeof(T));
   }

   /// Double the size of the allocated memory, guaranteeing space for at
   /// least one more element or MinSize if specified.
   void grow(size_t MinSize = 0)
   {
      this->grow_pod(MinSize, sizeof(T));
      this->Owns = true;
   }

public:
   using iterator = typename SuperClass::iterator;
   using const_iterator = typename SuperClass::const_iterator;
   using reference = typename SuperClass::reference;
   using size_type = typename SuperClass::size_type;

   void push_back(const T &Elt)
   {
      if (R__unlikely(this->size() >= this->capacity()))
         this->grow();
      memcpy(reinterpret_cast<void *>(this->end()), &Elt, sizeof(T));
      this->set_size(this->size() + 1);
   }

   void pop_back() { this->set_size(this->size() - 1); }
};

/// This class consists of common code factored out of the SmallVector class to
/// reduce code duplication based on the SmallVector 'N' template parameter.
template <typename T>
class SmallVectorImpl : public SmallVectorTemplateBase<T> {
   using SuperClass = SmallVectorTemplateBase<T>;

public:
   using iterator = typename SuperClass::iterator;
   using const_iterator = typename SuperClass::const_iterator;
   using reference = typename SuperClass::reference;
   using size_type = typename SuperClass::size_type;

protected:
   // Default ctor - Initialize to empty.
   explicit SmallVectorImpl(unsigned N) : SmallVectorTemplateBase<T>(N) {}

public:
   SmallVectorImpl(const SmallVectorImpl &) = delete;

   ~SmallVectorImpl()
   {
      // Subclass has already destructed this vector's elements.
      // If this wasn't grown from the inline copy, deallocate the old space.
      if (!this->isSmall() && this->Owns)
         free(this->begin());
   }

   // also give up adopted memory if applicable
   void clear()
   {
      if (this->Owns) {
         this->destroy_range(this->begin(), this->end());
         this->fSize = 0;
      } else {
         this->resetToSmall();
      }
   }

   void resize(size_type N)
   {
      if (N < this->size()) {
         if (this->Owns)
            this->destroy_range(this->begin() + N, this->end());
         this->set_size(N);
      } else if (N > this->size()) {
         if (this->capacity() < N)
            this->grow(N);
         for (auto I = this->end(), E = this->begin() + N; I != E; ++I)
            new (&*I) T();
         this->set_size(N);
      }
   }

   void resize(size_type N, const T &NV)
   {
      if (N < this->size()) {
         if (this->Owns)
            this->destroy_range(this->begin() + N, this->end());
         this->set_size(N);
      } else if (N > this->size()) {
         if (this->capacity() < N)
            this->grow(N);
         std::uninitialized_fill(this->end(), this->begin() + N, NV);
         this->set_size(N);
      }
   }

   void reserve(size_type N)
   {
      if (this->capacity() < N)
         this->grow(N);
   }

   void pop_back_n(size_type NumItems)
   {
      assert(this->size() >= NumItems);
      if (this->Owns)
         this->destroy_range(this->end() - NumItems, this->end());
      this->set_size(this->size() - NumItems);
   }

   R__RVEC_NODISCARD T pop_back_val()
   {
      T Result = ::std::move(this->back());
      this->pop_back();
      return Result;
   }

   void swap(SmallVectorImpl &RHS);

   /// Add the specified range to the end of the SmallVector.
   template <typename in_iter,
             typename = typename std::enable_if<std::is_convertible<
                typename std::iterator_traits<in_iter>::iterator_category, std::input_iterator_tag>::value>::type>
   void append(in_iter in_start, in_iter in_end)
   {
      size_type NumInputs = std::distance(in_start, in_end);
      if (NumInputs > this->capacity() - this->size())
         this->grow(this->size() + NumInputs);

      this->uninitialized_copy(in_start, in_end, this->end());
      this->set_size(this->size() + NumInputs);
   }

   /// Append \p NumInputs copies of \p Elt to the end.
   void append(size_type NumInputs, const T &Elt)
   {
      if (NumInputs > this->capacity() - this->size())
         this->grow(this->size() + NumInputs);

      std::uninitialized_fill_n(this->end(), NumInputs, Elt);
      this->set_size(this->size() + NumInputs);
   }

   void append(std::initializer_list<T> IL) { append(IL.begin(), IL.end()); }

   // from the original LLVM implementation:
   // FIXME: Consider assigning over existing elements, rather than clearing &
   // re-initializing them - for all assign(...) variants.

   void assign(size_type NumElts, const T &Elt)
   {
      clear();
      if (this->capacity() < NumElts)
         this->grow(NumElts);
      this->set_size(NumElts);
      std::uninitialized_fill(this->begin(), this->end(), Elt);
   }

   template <typename in_iter,
             typename = typename std::enable_if<std::is_convertible<
                typename std::iterator_traits<in_iter>::iterator_category, std::input_iterator_tag>::value>::type>
   void assign(in_iter in_start, in_iter in_end)
   {
      clear();
      append(in_start, in_end);
   }

   void assign(std::initializer_list<T> IL)
   {
      clear();
      append(IL);
   }

   iterator erase(const_iterator CI)
   {
      // Just cast away constness because this is a non-const member function.
      iterator I = const_cast<iterator>(CI);

      assert(I >= this->begin() && "Iterator to erase is out of bounds.");
      assert(I < this->end() && "Erasing at past-the-end iterator.");

      iterator N = I;
      // Shift all elts down one.
      std::move(I + 1, this->end(), I);
      // Drop the last elt.
      this->pop_back();
      return (N);
   }

   iterator erase(const_iterator CS, const_iterator CE)
   {
      // Just cast away constness because this is a non-const member function.
      iterator S = const_cast<iterator>(CS);
      iterator E = const_cast<iterator>(CE);

      assert(S >= this->begin() && "Range to erase is out of bounds.");
      assert(S <= E && "Trying to erase invalid range.");
      assert(E <= this->end() && "Trying to erase past the end.");

      iterator N = S;
      // Shift all elts down.
      iterator I = std::move(E, this->end(), S);
      // Drop the last elts.
      if (this->Owns)
         this->destroy_range(I, this->end());
      this->set_size(I - this->begin());
      return (N);
   }

   iterator insert(iterator I, T &&Elt)
   {
      if (I == this->end()) { // Important special case for empty vector.
         this->push_back(::std::move(Elt));
         return this->end() - 1;
      }

      assert(I >= this->begin() && "Insertion iterator is out of bounds.");
      assert(I <= this->end() && "Inserting past the end of the vector.");

      if (this->size() >= this->capacity()) {
         size_t EltNo = I - this->begin();
         this->grow();
         I = this->begin() + EltNo;
      }

      ::new ((void *)this->end()) T(::std::move(this->back()));
      // Push everything else over.
      std::move_backward(I, this->end() - 1, this->end());
      this->set_size(this->size() + 1);

      // If we just moved the element we're inserting, be sure to update
      // the reference.
      T *EltPtr = &Elt;
      if (I <= EltPtr && EltPtr < this->end())
         ++EltPtr;

      *I = ::std::move(*EltPtr);
      return I;
   }

   iterator insert(iterator I, const T &Elt)
   {
      if (I == this->end()) { // Important special case for empty vector.
         this->push_back(Elt);
         return this->end() - 1;
      }

      assert(I >= this->begin() && "Insertion iterator is out of bounds.");
      assert(I <= this->end() && "Inserting past the end of the vector.");

      if (this->size() >= this->capacity()) {
         size_t EltNo = I - this->begin();
         this->grow();
         I = this->begin() + EltNo;
      }
      ::new ((void *)this->end()) T(std::move(this->back()));
      // Push everything else over.
      std::move_backward(I, this->end() - 1, this->end());
      this->set_size(this->size() + 1);

      // If we just moved the element we're inserting, be sure to update
      // the reference.
      const T *EltPtr = &Elt;
      if (I <= EltPtr && EltPtr < this->end())
         ++EltPtr;

      *I = *EltPtr;
      return I;
   }

   iterator insert(iterator I, size_type NumToInsert, const T &Elt)
   {
      // Convert iterator to elt# to avoid invalidating iterator when we reserve()
      size_t InsertElt = I - this->begin();

      if (I == this->end()) { // Important special case for empty vector.
         append(NumToInsert, Elt);
         return this->begin() + InsertElt;
      }

      assert(I >= this->begin() && "Insertion iterator is out of bounds.");
      assert(I <= this->end() && "Inserting past the end of the vector.");

      // Ensure there is enough space.
      reserve(this->size() + NumToInsert);

      // Uninvalidate the iterator.
      I = this->begin() + InsertElt;

      // If there are more elements between the insertion point and the end of the
      // range than there are being inserted, we can use a simple approach to
      // insertion.  Since we already reserved space, we know that this won't
      // reallocate the vector.
      if (size_t(this->end() - I) >= NumToInsert) {
         T *OldEnd = this->end();
         append(std::move_iterator<iterator>(this->end() - NumToInsert), std::move_iterator<iterator>(this->end()));

         // Copy the existing elements that get replaced.
         std::move_backward(I, OldEnd - NumToInsert, OldEnd);

         std::fill_n(I, NumToInsert, Elt);
         return I;
      }

      // Otherwise, we're inserting more elements than exist already, and we're
      // not inserting at the end.

      // Move over the elements that we're about to overwrite.
      T *OldEnd = this->end();
      this->set_size(this->size() + NumToInsert);
      size_t NumOverwritten = OldEnd - I;
      this->uninitialized_move(I, OldEnd, this->end() - NumOverwritten);

      // Replace the overwritten part.
      std::fill_n(I, NumOverwritten, Elt);

      // Insert the non-overwritten middle part.
      std::uninitialized_fill_n(OldEnd, NumToInsert - NumOverwritten, Elt);
      return I;
   }

   template <typename ItTy,
             typename = typename std::enable_if<std::is_convertible<
                typename std::iterator_traits<ItTy>::iterator_category, std::input_iterator_tag>::value>::type>
   iterator insert(iterator I, ItTy From, ItTy To)
   {
      // Convert iterator to elt# to avoid invalidating iterator when we reserve()
      size_t InsertElt = I - this->begin();

      if (I == this->end()) { // Important special case for empty vector.
         append(From, To);
         return this->begin() + InsertElt;
      }

      assert(I >= this->begin() && "Insertion iterator is out of bounds.");
      assert(I <= this->end() && "Inserting past the end of the vector.");

      size_t NumToInsert = std::distance(From, To);

      // Ensure there is enough space.
      reserve(this->size() + NumToInsert);

      // Uninvalidate the iterator.
      I = this->begin() + InsertElt;

      // If there are more elements between the insertion point and the end of the
      // range than there are being inserted, we can use a simple approach to
      // insertion.  Since we already reserved space, we know that this won't
      // reallocate the vector.
      if (size_t(this->end() - I) >= NumToInsert) {
         T *OldEnd = this->end();
         append(std::move_iterator<iterator>(this->end() - NumToInsert), std::move_iterator<iterator>(this->end()));

         // Copy the existing elements that get replaced.
         std::move_backward(I, OldEnd - NumToInsert, OldEnd);

         std::copy(From, To, I);
         return I;
      }

      // Otherwise, we're inserting more elements than exist already, and we're
      // not inserting at the end.

      // Move over the elements that we're about to overwrite.
      T *OldEnd = this->end();
      this->set_size(this->size() + NumToInsert);
      size_t NumOverwritten = OldEnd - I;
      this->uninitialized_move(I, OldEnd, this->end() - NumOverwritten);

      // Replace the overwritten part.
      for (T *J = I; NumOverwritten > 0; --NumOverwritten) {
         *J = *From;
         ++J;
         ++From;
      }

      // Insert the non-overwritten middle part.
      this->uninitialized_copy(From, To, OldEnd);
      return I;
   }

   void insert(iterator I, std::initializer_list<T> IL) { insert(I, IL.begin(), IL.end()); }

   template <typename... ArgTypes>
   R__DEPRECATED(6, 28, "Please use RVec::insert instead.")
   reference emplace(iterator I, ArgTypes &&...Args)
   {
      // these are not the exact semantics of emplace, of course, hence the deprecation.
      return insert(I, T(std::forward<Args>...));
   }

   template <typename... ArgTypes>
   reference emplace_back(ArgTypes &&...Args)
   {
      if (R__unlikely(this->size() >= this->capacity()))
         this->grow();
      ::new ((void *)this->end()) T(std::forward<ArgTypes>(Args)...);
      this->set_size(this->size() + 1);
      return this->back();
   }

   SmallVectorImpl &operator=(const SmallVectorImpl &RHS);

   SmallVectorImpl &operator=(SmallVectorImpl &&RHS);

   bool operator==(const SmallVectorImpl &RHS) const
   {
      if (this->size() != RHS.size())
         return false;
      return std::equal(this->begin(), this->end(), RHS.begin());
   }
   bool operator!=(const SmallVectorImpl &RHS) const { return !(*this == RHS); }

   bool operator<(const SmallVectorImpl &RHS) const
   {
      return std::lexicographical_compare(this->begin(), this->end(), RHS.begin(), RHS.end());
   }
};

template <typename T>
void SmallVectorImpl<T>::swap(SmallVectorImpl<T> &RHS)
{
   if (this == &RHS)
      return;

   // We can only avoid copying elements if neither vector is small.
   if (!this->isSmall() && !RHS.isSmall()) {
      std::swap(this->BeginX, RHS.BeginX);
      std::swap(this->fSize, RHS.fSize);
      std::swap(this->Capacity, RHS.Capacity);
      std::swap(this->Owns, RHS.Owns);
      return;
   }
   if (RHS.size() > this->capacity())
      this->grow(RHS.size());
   if (this->size() > RHS.capacity())
      RHS.grow(this->size());

   // Swap the shared elements.
   size_t NumShared = this->size();
   if (NumShared > RHS.size())
      NumShared = RHS.size();
   for (size_type i = 0; i != NumShared; ++i)
      std::swap((*this)[i], RHS[i]);

   // Copy over the extra elts.
   if (this->size() > RHS.size()) {
      size_t EltDiff = this->size() - RHS.size();
      this->uninitialized_copy(this->begin() + NumShared, this->end(), RHS.end());
      RHS.set_size(RHS.size() + EltDiff);
      if (this->Owns)
         this->destroy_range(this->begin() + NumShared, this->end());
      this->set_size(NumShared);
   } else if (RHS.size() > this->size()) {
      size_t EltDiff = RHS.size() - this->size();
      this->uninitialized_copy(RHS.begin() + NumShared, RHS.end(), this->end());
      this->set_size(this->size() + EltDiff);
      if (RHS.Owns)
         this->destroy_range(RHS.begin() + NumShared, RHS.end());
      RHS.set_size(NumShared);
   }
}

template <typename T>
SmallVectorImpl<T> &SmallVectorImpl<T>::operator=(const SmallVectorImpl<T> &RHS)
{
   // Avoid self-assignment.
   if (this == &RHS)
      return *this;

   // If we already have sufficient space, assign the common elements, then
   // destroy any excess.
   size_t RHSSize = RHS.size();
   size_t CurSize = this->size();
   if (CurSize >= RHSSize) {
      // Assign common elements.
      iterator NewEnd;
      if (RHSSize)
         NewEnd = std::copy(RHS.begin(), RHS.begin() + RHSSize, this->begin());
      else
         NewEnd = this->begin();

      // Destroy excess elements.
      if (this->Owns)
         this->destroy_range(NewEnd, this->end());

      // Trim.
      this->set_size(RHSSize);
      return *this;
   }

   // If we have to grow to have enough elements, destroy the current elements.
   // This allows us to avoid copying them during the grow.
   // From the original LLVM implementation:
   // FIXME: don't do this if they're efficiently moveable.
   if (this->capacity() < RHSSize) {
      if (this->Owns) {
         // Destroy current elements.
         this->destroy_range(this->begin(), this->end());
      }
      this->set_size(0);
      CurSize = 0;
      this->grow(RHSSize);
   } else if (CurSize) {
      // Otherwise, use assignment for the already-constructed elements.
      std::copy(RHS.begin(), RHS.begin() + CurSize, this->begin());
   }

   // Copy construct the new elements in place.
   this->uninitialized_copy(RHS.begin() + CurSize, RHS.end(), this->begin() + CurSize);

   // Set end.
   this->set_size(RHSSize);
   return *this;
}

template <typename T>
SmallVectorImpl<T> &SmallVectorImpl<T>::operator=(SmallVectorImpl<T> &&RHS)
{
   // Avoid self-assignment.
   if (this == &RHS)
      return *this;

   // If the RHS isn't small, clear this vector and then steal its buffer.
   if (!RHS.isSmall()) {
      if (this->Owns) {
         this->destroy_range(this->begin(), this->end());
         if (!this->isSmall())
            free(this->begin());
      }
      this->BeginX = RHS.BeginX;
      this->fSize = RHS.fSize;
      this->Capacity = RHS.Capacity;
      this->Owns = RHS.Owns;
      RHS.resetToSmall();
      return *this;
   }

   // If we already have sufficient space, assign the common elements, then
   // destroy any excess.
   size_t RHSSize = RHS.size();
   size_t CurSize = this->size();
   if (CurSize >= RHSSize) {
      // Assign common elements.
      iterator NewEnd = this->begin();
      if (RHSSize)
         NewEnd = std::move(RHS.begin(), RHS.end(), NewEnd);

      // Destroy excess elements and trim the bounds.
      if (this->Owns)
         this->destroy_range(NewEnd, this->end());
      this->set_size(RHSSize);

      // Clear the RHS.
      RHS.clear();

      return *this;
   }

   // If we have to grow to have enough elements, destroy the current elements.
   // This allows us to avoid copying them during the grow.
   // From the original LLVM implementation:
   // FIXME: this may not actually make any sense if we can efficiently move
   // elements.
   if (this->capacity() < RHSSize) {
      if (this->Owns) {
         // Destroy current elements.
         this->destroy_range(this->begin(), this->end());
      }
      this->set_size(0);
      CurSize = 0;
      this->grow(RHSSize);
   } else if (CurSize) {
      // Otherwise, use assignment for the already-constructed elements.
      std::move(RHS.begin(), RHS.begin() + CurSize, this->begin());
   }

   // Move-construct the new elements in place.
   this->uninitialized_move(RHS.begin() + CurSize, RHS.end(), this->begin() + CurSize);

   // Set end.
   this->set_size(RHSSize);

   RHS.clear();
   return *this;
}

/// Storage for the SmallVector elements.  This is specialized for the N=0 case
/// to avoid allocating unnecessary storage.
template <typename T, unsigned N>
struct SmallVectorStorage {
   alignas(T) char InlineElts[N * sizeof(T)];
};

/// We need the storage to be properly aligned even for small-size of 0 so that
/// the pointer math in \a SmallVectorTemplateCommon::getFirstEl() is
/// well-defined.
template <typename T>
struct alignas(T) SmallVectorStorage<T, 0> {
};
} // namespace VecOps
} // namespace Internal

namespace VecOps {
// Note that we open here with @{ the Doxygen group vecops and it is
// closed again at the end of the C++ namespace VecOps
/**
  * \defgroup vecops VecOps
  * A "std::vector"-like collection of values implementing handy operation to analyse them
  * @{
*/

// clang-format off
/**
\class ROOT::VecOps::RVec
\brief A "std::vector"-like collection of values implementing handy operation to analyse them
\tparam T The type of the contained objects

A RVec is a container designed to make analysis of values' collections fast and easy.
Its storage is contiguous in memory and its interface is designed such to resemble to the one
of the stl vector. In addition the interface features methods and external functions to ease
the manipulation and analysis of the data in the RVec.

\htmlonly
<a href="https://doi.org/10.5281/zenodo.1253756"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.1253756.svg" alt="DOI"></a>
\endhtmlonly

## Table of Contents
- [Example](\ref example)
- [Owning and adopting memory](\ref owningandadoptingmemory)
- [Sorting and manipulation of indices](\ref sorting)
- [Usage in combination with RDataFrame](\ref usagetdataframe)
- [Reference for the RVec class](\ref RVecdoxyref)

Also see the [reference for RVec helper functions](https://root.cern/doc/master/namespaceROOT_1_1VecOps.html).

\anchor example
## Example
Suppose to have an event featuring a collection of muons with a certain pseudorapidity,
momentum and charge, e.g.:
~~~{.cpp}
std::vector<short> mu_charge {1, 1, -1, -1, -1, 1, 1, -1};
std::vector<float> mu_pt {56, 45, 32, 24, 12, 8, 7, 6.2};
std::vector<float> mu_eta {3.1, -.2, -1.1, 1, 4.1, 1.6, 2.4, -.5};
~~~
Suppose you want to extract the transverse momenta of the muons satisfying certain
criteria, for example consider only negatively charged muons with a pseudorapidity
smaller or equal to 2 and with a transverse momentum greater than 10 GeV.
Such a selection would require, among the other things, the management of an explicit
loop, for example:
~~~{.cpp}
std::vector<float> goodMuons_pt;
const auto size = mu_charge.size();
for (size_t i=0; i < size; ++i) {
   if (mu_pt[i] > 10 && abs(mu_eta[i]) <= 2. &&  mu_charge[i] == -1) {
      goodMuons_pt.emplace_back(mu_pt[i]);
   }
}
~~~
These operations become straightforward with RVec - we just need to *write what
we mean*:
~~~{.cpp}
auto goodMuons_pt = mu_pt[ (mu_pt > 10.f && abs(mu_eta) <= 2.f && mu_charge == -1) ]
~~~
Now the clean collection of transverse momenta can be used within the rest of the data analysis, for
example to fill a histogram.

\anchor owningandadoptingmemory
## Owning and adopting memory
RVec has contiguous memory associated to it. It can own it or simply adopt it. In the latter case,
it can be constructed with the address of the memory associated to it and its length. For example:
~~~{.cpp}
std::vector<int> myStlVec {1,2,3};
RVec<int> myRVec(myStlVec.data(), myStlVec.size());
~~~
In this case, the memory associated to myStlVec and myRVec is the same, myRVec simply "adopted it".
If any method which implies a re-allocation is called, e.g. *emplace_back* or *resize*, the adopted
memory is released and new one is allocated. The previous content is copied in the new memory and
preserved.

\anchor sorting
## Sorting and manipulation of indices

### Sorting
RVec complies to the STL interfaces when it comes to iterations. As a result, standard algorithms
can be used, for example sorting:
~~~{.cpp}
RVec<double> v{6., 4., 5.};
std::sort(v.begin(), v.end());
~~~

For convinience, helpers are provided too:
~~~{.cpp}
auto sorted_v = Sort(v);
auto reversed_v = Reverse(v);
~~~

### Manipulation of indices

It is also possible to manipulated the RVecs acting on their indices. For example,
the following syntax
~~~{.cpp}
RVec<double> v0 {9., 7., 8.};
auto v1 = Take(v0, {1, 2, 0});
~~~
will yield a new RVec<double> the content of which is the first, second and zeroth element of
v0, i.e. `{7., 8., 9.}`.

The `Argsort` helper extracts the indices which order the content of a `RVec`. For example,
this snippet accomplish in a more expressive way what we just achieved:
~~~{.cpp}
auto v1_indices = Argsort(v0); // The content of v1_indices is {1, 2, 0}.
v1 = Take(v0, v1_indices);
~~~

The `Take` utility allows to extract portions of the `RVec`. The content to be *taken*
can be specified with an `RVec` of indices or an integer. If the integer is negative,
elements will be picked starting from the end of the container:
~~~{.cpp}
RVec<float> vf {1.f, 2.f, 3.f, 4.f};
auto vf_1 = Take(vf, {1, 3}); // The content is {2.f, 4.f}
auto vf_2 = Take(vf, 2); // The content is {1.f, 2.f}
auto vf_3 = Take(vf, -3); // The content is {2.f, 3.f, 4.f}
~~~

\anchor usagetdataframe
## Usage in combination with RDataFrame
RDataFrame leverages internally RVecs. Suppose to have a dataset stored in a
TTree which holds these columns (here we choose C arrays to represent the
collections, they could be as well std::vector instances):
~~~{.bash}
  nPart            "nPart/I"            An integer representing the number of particles
  px               "px[nPart]/D"        The C array of the particles' x component of the momentum
  py               "py[nPart]/D"        The C array of the particles' y component of the momentum
  E                "E[nPart]/D"         The C array of the particles' Energy
~~~
Suppose you'd like to plot in a histogram the transverse momenta of all particles
for which the energy is greater than 200 MeV.
The code required would just be:
~~~{.cpp}
RDataFrame d("mytree", "myfile.root");
using doubles = RVec<double>;
auto cutPt = [](doubles &pxs, doubles &pys, doubles &Es) {
   auto all_pts = sqrt(pxs * pxs + pys * pys);
   auto good_pts = all_pts[Es > 200.];
   return good_pts;
   };

auto hpt = d.Define("pt", cutPt, {"px", "py", "E"})
            .Histo1D("pt");
hpt->Draw();
~~~
And if you'd like to express your selection as a string:
~~~{.cpp}
RDataFrame d("mytree", "myfile.root");
auto hpt = d.Define("pt", "sqrt(pxs * pxs + pys * pys)[E>200]")
            .Histo1D("pt");
hpt->Draw();
~~~
\anchor RVecdoxyref
**/
// clang-format on



/// This is a 'vector' (really, a variable-sized array), optimized
/// for the case when the array is small.  It contains some number of elements
/// in-place, which allows it to avoid heap allocation when the actual number of
/// elements is below that threshold.  This allows normal "small" cases to be
/// fast without losing generality for large inputs.
///
/// Note that this does not attempt to be exception safe.
///
template <typename T>
class RVec : public Internal::VecOps::SmallVectorImpl<T>, Internal::VecOps::SmallVectorStorage<T, 8> {
public:
   static constexpr unsigned N = 8;
   RVec() : Internal::VecOps::SmallVectorImpl<T>(N) {}

   ~RVec()
   {
      if (this->Owns) {
         // Destroy the constructed elements in the vector.
         this->destroy_range(this->begin(), this->end());
      }
   }

   explicit RVec(size_t Size, const T &Value = T()) : Internal::VecOps::SmallVectorImpl<T>(N)
   {
      this->assign(Size, Value);
   }

   template <typename ItTy,
             typename = typename std::enable_if<std::is_convertible<
                typename std::iterator_traits<ItTy>::iterator_category, std::input_iterator_tag>::value>::type>
   RVec(ItTy S, ItTy E) : Internal::VecOps::SmallVectorImpl<T>(N)
   {
      this->append(S, E);
   }

   RVec(std::initializer_list<T> IL) : Internal::VecOps::SmallVectorImpl<T>(N) { this->assign(IL); }

   RVec(const RVec &RHS) : Internal::VecOps::SmallVectorImpl<T>(N)
   {
      if (!RHS.empty())
         Internal::VecOps::SmallVectorImpl<T>::operator=(RHS);
   }

   RVec &operator=(const RVec &RHS)
   {
      Internal::VecOps::SmallVectorImpl<T>::operator=(RHS);
      return *this;
   }

   RVec(RVec &&RHS) : Internal::VecOps::SmallVectorImpl<T>(N)
   {
      if (!RHS.empty())
         Internal::VecOps::SmallVectorImpl<T>::operator=(::std::move(RHS));
   }

   RVec(Internal::VecOps::SmallVectorImpl<T> &&RHS) : Internal::VecOps::SmallVectorImpl<T>(N)
   {
      if (!RHS.empty())
         Internal::VecOps::SmallVectorImpl<T>::operator=(::std::move(RHS));
   }

   RVec(const std::vector<T> &RHS) : RVec(RHS.begin(), RHS.end()) {}

   RVec &operator=(RVec &&RHS)
   {
      Internal::VecOps::SmallVectorImpl<T>::operator=(::std::move(RHS));
      return *this;
   }

   RVec(T* p, size_t n) : Internal::VecOps::SmallVectorImpl<T>(N)
   {
      this->BeginX = p;
      this->fSize = n;
      this->Capacity = n;
      this->Owns = false;
   }

   RVec &operator=(Internal::VecOps::SmallVectorImpl<T> &&RHS)
   {
      Internal::VecOps::SmallVectorImpl<T>::operator=(::std::move(RHS));
      return *this;
   }

   RVec &operator=(std::initializer_list<T> IL)
   {
      this->assign(IL);
      return *this;
   }

   using reference = typename Internal::VecOps::SmallVectorTemplateCommon<T>::reference;
   using const_reference = typename Internal::VecOps::SmallVectorTemplateCommon<T>::const_reference;
   using size_type = typename Internal::VecOps::SmallVectorTemplateCommon<T>::size_type;
   using Internal::VecOps::SmallVectorTemplateCommon<T>::begin;
   using Internal::VecOps::SmallVectorTemplateCommon<T>::size;

   reference operator[](size_type idx)
   {
      assert(idx < size());
      return begin()[idx];
   }
   const_reference operator[](size_type idx) const
   {
      assert(idx < size());
      return begin()[idx];
   }

   template <typename V, typename = std::enable_if<std::is_convertible<V, bool>::value>>
   RVec operator[](const RVec<V> &conds) const
   {
      const auto n = conds.size();

      if (n != this->size())
         throw std::runtime_error("Cannot index RVec with condition vector of different size");

      RVec<T> ret;
      ret.reserve(n);
      for (auto i = 0u; i < n; ++i)
         if (conds[i])
            ret.emplace_back(this->operator[](i));
      return ret;
   }

   // conversion
   template <typename U, typename = std::enable_if<std::is_convertible<T, U>::value>>
   operator RVec<U>() const
   {
      return RVec<U>(this->begin(), this->end());
   }

   // from the original LLVM implementation:
   // FIXME find a less verbose way
   typename Internal::VecOps::SmallVectorTemplateCommon<T>::reference at(typename Internal::VecOps::SmallVectorTemplateCommon<T>::size_type pos)
   {
      if (pos >= this->fSize)
         throw std::out_of_range("RVec");
      return this->operator[](pos);
   }
   typename Internal::VecOps::SmallVectorTemplateCommon<T>::const_reference at(typename Internal::VecOps::SmallVectorTemplateCommon<T>::size_type pos) const
   {
      if (pos >= this->fSize)
         throw std::out_of_range("RVec");
      return this->operator[](pos);
   }
   /// No exception thrown. The user specifies the desired value in case the RVec is shorter than `pos`.
   typename Internal::VecOps::SmallVectorTemplateCommon<T>::value_type
   at(typename Internal::VecOps::SmallVectorTemplateCommon<T>::size_type pos, typename Internal::VecOps::SmallVectorTemplateCommon<T>::value_type fallback)
   {
      if (pos >= this->fSize)
         return fallback;
      return this->operator[](pos);
   }
   /// No exception thrown. The user specifies the desired value in case the RVec is shorter than `pos`.
   typename Internal::VecOps::SmallVectorTemplateCommon<T>::value_type
   at(typename Internal::VecOps::SmallVectorTemplateCommon<T>::size_type pos, typename Internal::VecOps::SmallVectorTemplateCommon<T>::value_type fallback) const
   {
      if (pos >= this->fSize)
         return fallback;
      return this->operator[](pos);
   }
};

template <typename T>
inline size_t capacity_in_bytes(const RVec<T> &X)
{
   return X.capacity_in_bytes();
}

///@name RVec Unary Arithmetic Operators
///@{

#define RVEC_UNARY_OPERATOR(OP)                                                \
template <typename T>                                                          \
RVec<T> operator OP(const RVec<T> &v)                                          \
{                                                                              \
   RVec<T> ret(v);                                                             \
   for (auto &x : ret)                                                         \
      x = OP x;                                                                \
return ret;                                                                    \
}                                                                              \

RVEC_UNARY_OPERATOR(+)
RVEC_UNARY_OPERATOR(-)
RVEC_UNARY_OPERATOR(~)
RVEC_UNARY_OPERATOR(!)
#undef RVEC_UNARY_OPERATOR

///@}
///@name RVec Binary Arithmetic Operators
///@{

#define ERROR_MESSAGE(OP) \
 "Cannot call operator " #OP " on vectors of different sizes."

#define RVEC_BINARY_OPERATOR(OP)                                               \
template <typename T0, typename T1>                                            \
auto operator OP(const RVec<T0> &v, const T1 &y)                               \
  -> RVec<decltype(v[0] OP y)>                                                 \
{                                                                              \
   RVec<decltype(v[0] OP y)> ret(v.size());                                    \
   auto op = [&y](const T0 &x) { return x OP y; };                             \
   std::transform(v.begin(), v.end(), ret.begin(), op);                        \
   return ret;                                                                 \
}                                                                              \
                                                                               \
template <typename T0, typename T1>                                            \
auto operator OP(const T0 &x, const RVec<T1> &v)                               \
  -> RVec<decltype(x OP v[0])>                                                 \
{                                                                              \
   RVec<decltype(x OP v[0])> ret(v.size());                                    \
   auto op = [&x](const T1 &y) { return x OP y; };                             \
   std::transform(v.begin(), v.end(), ret.begin(), op);                        \
   return ret;                                                                 \
}                                                                              \
                                                                               \
template <typename T0, typename T1>                                            \
auto operator OP(const RVec<T0> &v0, const RVec<T1> &v1)                       \
  -> RVec<decltype(v0[0] OP v1[0])>                                            \
{                                                                              \
   if (v0.size() != v1.size())                                                 \
      throw std::runtime_error(ERROR_MESSAGE(OP));                             \
                                                                               \
   RVec<decltype(v0[0] OP v1[0])> ret(v0.size());                              \
   auto op = [](const T0 &x, const T1 &y) { return x OP y; };                  \
   std::transform(v0.begin(), v0.end(), v1.begin(), ret.begin(), op);          \
   return ret;                                                                 \
}                                                                              \

RVEC_BINARY_OPERATOR(+)
RVEC_BINARY_OPERATOR(-)
RVEC_BINARY_OPERATOR(*)
RVEC_BINARY_OPERATOR(/)
RVEC_BINARY_OPERATOR(%)
RVEC_BINARY_OPERATOR(^)
RVEC_BINARY_OPERATOR(|)
RVEC_BINARY_OPERATOR(&)
#undef RVEC_BINARY_OPERATOR

///@}
///@name RVec Assignment Arithmetic Operators
///@{

#define RVEC_ASSIGNMENT_OPERATOR(OP)                                           \
template <typename T0, typename T1>                                            \
RVec<T0>& operator OP(RVec<T0> &v, const T1 &y)                                \
{                                                                              \
   auto op = [&y](T0 &x) { return x OP y; };                                   \
   std::transform(v.begin(), v.end(), v.begin(), op);                          \
   return v;                                                                   \
}                                                                              \
                                                                               \
template <typename T0, typename T1>                                            \
RVec<T0>& operator OP(RVec<T0> &v0, const RVec<T1> &v1)                        \
{                                                                              \
   if (v0.size() != v1.size())                                                 \
      throw std::runtime_error(ERROR_MESSAGE(OP));                             \
                                                                               \
   auto op = [](T0 &x, const T1 &y) { return x OP y; };                        \
   std::transform(v0.begin(), v0.end(), v1.begin(), v0.begin(), op);           \
   return v0;                                                                  \
}                                                                              \

RVEC_ASSIGNMENT_OPERATOR(+=)
RVEC_ASSIGNMENT_OPERATOR(-=)
RVEC_ASSIGNMENT_OPERATOR(*=)
RVEC_ASSIGNMENT_OPERATOR(/=)
RVEC_ASSIGNMENT_OPERATOR(%=)
RVEC_ASSIGNMENT_OPERATOR(^=)
RVEC_ASSIGNMENT_OPERATOR(|=)
RVEC_ASSIGNMENT_OPERATOR(&=)
RVEC_ASSIGNMENT_OPERATOR(>>=)
RVEC_ASSIGNMENT_OPERATOR(<<=)
#undef RVEC_ASSIGNMENT_OPERATOR

///@}
///@name RVec Comparison and Logical Operators
///@{

#define RVEC_LOGICAL_OPERATOR(OP)                                              \
template <typename T0, typename T1>                                            \
auto operator OP(const RVec<T0> &v, const T1 &y)                               \
  -> RVec<int> /* avoid std::vector<bool> */                                   \
{                                                                              \
   RVec<int> ret(v.size());                                                    \
   auto op = [y](const T0 &x) -> int { return x OP y; };                       \
   std::transform(v.begin(), v.end(), ret.begin(), op);                        \
   return ret;                                                                 \
}                                                                              \
                                                                               \
template <typename T0, typename T1>                                            \
auto operator OP(const T0 &x, const RVec<T1> &v)                               \
  -> RVec<int> /* avoid std::vector<bool> */                                   \
{                                                                              \
   RVec<int> ret(v.size());                                                    \
   auto op = [x](const T1 &y) -> int { return x OP y; };                       \
   std::transform(v.begin(), v.end(), ret.begin(), op);                        \
   return ret;                                                                 \
}                                                                              \
                                                                               \
template <typename T0, typename T1>                                            \
auto operator OP(const RVec<T0> &v0, const RVec<T1> &v1)                       \
  -> RVec<int> /* avoid std::vector<bool> */                                   \
{                                                                              \
   if (v0.size() != v1.size())                                                 \
      throw std::runtime_error(ERROR_MESSAGE(OP));                             \
                                                                               \
   RVec<int> ret(v0.size());                                                   \
   auto op = [](const T0 &x, const T1 &y) -> int { return x OP y; };           \
   std::transform(v0.begin(), v0.end(), v1.begin(), ret.begin(), op);          \
   return ret;                                                                 \
}                                                                              \

RVEC_LOGICAL_OPERATOR(<)
RVEC_LOGICAL_OPERATOR(>)
RVEC_LOGICAL_OPERATOR(==)
RVEC_LOGICAL_OPERATOR(!=)
RVEC_LOGICAL_OPERATOR(<=)
RVEC_LOGICAL_OPERATOR(>=)
RVEC_LOGICAL_OPERATOR(&&)
RVEC_LOGICAL_OPERATOR(||)
#undef RVEC_LOGICAL_OPERATOR

///@}
///@name RVec Standard Mathematical Functions
///@{

/// \cond
template <typename T> struct PromoteTypeImpl;

template <> struct PromoteTypeImpl<float>       { using Type = float;       };
template <> struct PromoteTypeImpl<double>      { using Type = double;      };
template <> struct PromoteTypeImpl<long double> { using Type = long double; };

template <typename T> struct PromoteTypeImpl { using Type = double; };

template <typename T>
using PromoteType = typename PromoteTypeImpl<T>::Type;

template <typename U, typename V>
using PromoteTypes = decltype(PromoteType<U>() + PromoteType<V>());

/// \endcond

#define RVEC_UNARY_FUNCTION(NAME, FUNC)                                        \
   template <typename T>                                                       \
   RVec<PromoteType<T>> NAME(const RVec<T> &v)                                 \
   {                                                                           \
      RVec<PromoteType<T>> ret(v.size());                                      \
      auto f = [](const T &x) { return FUNC(x); };                             \
      std::transform(v.begin(), v.end(), ret.begin(), f);                      \
      return ret;                                                              \
   }

#define RVEC_BINARY_FUNCTION(NAME, FUNC)                                       \
   template <typename T0, typename T1>                                         \
   RVec<PromoteTypes<T0, T1>> NAME(const T0 &x, const RVec<T1> &v)             \
   {                                                                           \
      RVec<PromoteTypes<T0, T1>> ret(v.size());                                \
      auto f = [&x](const T1 &y) { return FUNC(x, y); };                       \
      std::transform(v.begin(), v.end(), ret.begin(), f);                      \
      return ret;                                                              \
   }                                                                           \
                                                                               \
   template <typename T0, typename T1>                                         \
   RVec<PromoteTypes<T0, T1>> NAME(const RVec<T0> &v, const T1 &y)             \
   {                                                                           \
      RVec<PromoteTypes<T0, T1>> ret(v.size());                                \
      auto f = [&y](const T1 &x) { return FUNC(x, y); };                       \
      std::transform(v.begin(), v.end(), ret.begin(), f);                      \
      return ret;                                                              \
   }                                                                           \
                                                                               \
   template <typename T0, typename T1>                                         \
   RVec<PromoteTypes<T0, T1>> NAME(const RVec<T0> &v0, const RVec<T1> &v1)     \
   {                                                                           \
      if (v0.size() != v1.size())                                              \
         throw std::runtime_error(ERROR_MESSAGE(NAME));                        \
                                                                               \
      RVec<PromoteTypes<T0, T1>> ret(v0.size());                               \
      auto f = [](const T0 &x, const T1 &y) { return FUNC(x, y); };            \
      std::transform(v0.begin(), v0.end(), v1.begin(), ret.begin(), f);        \
      return ret;                                                              \
   }                                                                           \

#define RVEC_STD_UNARY_FUNCTION(F) RVEC_UNARY_FUNCTION(F, std::F)
#define RVEC_STD_BINARY_FUNCTION(F) RVEC_BINARY_FUNCTION(F, std::F)

RVEC_STD_UNARY_FUNCTION(abs)
RVEC_STD_BINARY_FUNCTION(fdim)
RVEC_STD_BINARY_FUNCTION(fmod)
RVEC_STD_BINARY_FUNCTION(remainder)

RVEC_STD_UNARY_FUNCTION(exp)
RVEC_STD_UNARY_FUNCTION(exp2)
RVEC_STD_UNARY_FUNCTION(expm1)

RVEC_STD_UNARY_FUNCTION(log)
RVEC_STD_UNARY_FUNCTION(log10)
RVEC_STD_UNARY_FUNCTION(log2)
RVEC_STD_UNARY_FUNCTION(log1p)

RVEC_STD_BINARY_FUNCTION(pow)
RVEC_STD_UNARY_FUNCTION(sqrt)
RVEC_STD_UNARY_FUNCTION(cbrt)
RVEC_STD_BINARY_FUNCTION(hypot)

RVEC_STD_UNARY_FUNCTION(sin)
RVEC_STD_UNARY_FUNCTION(cos)
RVEC_STD_UNARY_FUNCTION(tan)
RVEC_STD_UNARY_FUNCTION(asin)
RVEC_STD_UNARY_FUNCTION(acos)
RVEC_STD_UNARY_FUNCTION(atan)
RVEC_STD_BINARY_FUNCTION(atan2)

RVEC_STD_UNARY_FUNCTION(sinh)
RVEC_STD_UNARY_FUNCTION(cosh)
RVEC_STD_UNARY_FUNCTION(tanh)
RVEC_STD_UNARY_FUNCTION(asinh)
RVEC_STD_UNARY_FUNCTION(acosh)
RVEC_STD_UNARY_FUNCTION(atanh)

RVEC_STD_UNARY_FUNCTION(floor)
RVEC_STD_UNARY_FUNCTION(ceil)
RVEC_STD_UNARY_FUNCTION(trunc)
RVEC_STD_UNARY_FUNCTION(round)
RVEC_STD_UNARY_FUNCTION(lround)
RVEC_STD_UNARY_FUNCTION(llround)

RVEC_STD_UNARY_FUNCTION(erf)
RVEC_STD_UNARY_FUNCTION(erfc)
RVEC_STD_UNARY_FUNCTION(lgamma)
RVEC_STD_UNARY_FUNCTION(tgamma)
#undef RVEC_STD_UNARY_FUNCTION

///@}
///@name RVec Fast Mathematical Functions with Vdt
///@{

#ifdef R__HAS_VDT
#define RVEC_VDT_UNARY_FUNCTION(F) RVEC_UNARY_FUNCTION(F, vdt::F)

RVEC_VDT_UNARY_FUNCTION(fast_expf)
RVEC_VDT_UNARY_FUNCTION(fast_logf)
RVEC_VDT_UNARY_FUNCTION(fast_sinf)
RVEC_VDT_UNARY_FUNCTION(fast_cosf)
RVEC_VDT_UNARY_FUNCTION(fast_tanf)
RVEC_VDT_UNARY_FUNCTION(fast_asinf)
RVEC_VDT_UNARY_FUNCTION(fast_acosf)
RVEC_VDT_UNARY_FUNCTION(fast_atanf)

RVEC_VDT_UNARY_FUNCTION(fast_exp)
RVEC_VDT_UNARY_FUNCTION(fast_log)
RVEC_VDT_UNARY_FUNCTION(fast_sin)
RVEC_VDT_UNARY_FUNCTION(fast_cos)
RVEC_VDT_UNARY_FUNCTION(fast_tan)
RVEC_VDT_UNARY_FUNCTION(fast_asin)
RVEC_VDT_UNARY_FUNCTION(fast_acos)
RVEC_VDT_UNARY_FUNCTION(fast_atan)
#undef RVEC_VDT_UNARY_FUNCTION

#endif // R__HAS_VDT

#undef RVEC_UNARY_FUNCTION

///@}

/// Inner product
///
/// Example code, at the ROOT prompt:
/// ~~~{.cpp}
/// using namespace ROOT::VecOps;
/// RVec<float> v1 {1., 2., 3.};
/// RVec<float> v2 {4., 5., 6.};
/// auto v1_dot_v2 = Dot(v1, v2);
/// v1_dot_v2
/// // (float) 32.f
/// ~~~
template <typename T, typename V>
auto Dot(const RVec<T> &v0, const RVec<V> &v1) -> decltype(v0[0] * v1[0])
{
   if (v0.size() != v1.size())
      throw std::runtime_error("Cannot compute inner product of vectors of different sizes");
   return std::inner_product(v0.begin(), v0.end(), v1.begin(), decltype(v0[0] * v1[0])(0));
}

/// Sum elements of an RVec
///
/// Example code, at the ROOT prompt:
/// ~~~{.cpp}
/// using namespace ROOT::VecOps;
/// RVec<float> v {1.f, 2.f, 3.f};
/// auto v_sum = Sum(v);
/// v_sum
/// // (float) 6.f
/// ~~~
template <typename T>
T Sum(const RVec<T> &v)
{
   return std::accumulate(v.begin(), v.end(), T(0));
}

/// Get the mean of the elements of an RVec
///
/// The return type is a double precision floating point number.
/// Example code, at the ROOT prompt:
/// ~~~{.cpp}
/// using namespace ROOT::VecOps;
/// RVec<float> v {1.f, 2.f, 4.f};
/// auto v_mean = Mean(v);
/// v_mean
/// // (double) 2.3333333
/// ~~~
template <typename T>
double Mean(const RVec<T> &v)
{
   if (v.empty()) return 0.;
   return double(Sum(v)) / v.size();
}

/// Get the greatest element of an RVec
///
/// Example code, at the ROOT prompt:
/// ~~~~{.cpp}
/// using namespace ROOT::VecOps;
/// RVec<float> v {1.f, 2.f, 4.f};
/// auto v_max = Max(v)
/// v_max
/// (float) 4.f
/// ~~~~
template <typename T>
T Max(const RVec<T> &v)
{
   return *std::max_element(v.begin(), v.end());
}

/// Get the smallest element of an RVec
///
/// Example code, at the ROOT prompt:
/// ~~~~{.cpp}
/// using namespace ROOT::VecOps;
/// RVec<float> v {1.f, 2.f, 4.f};
/// auto v_min = Min(v)
/// v_min
/// (float) 1.f
/// ~~~~
template <typename T>
T Min(const RVec<T> &v)
{
   return *std::min_element(v.begin(), v.end());
}

/// Get the index of the greatest element of an RVec
/// In case of multiple occurrences of the maximum values,
/// the index corresponding to the first occurrence is returned.
///
/// Example code, at the ROOT prompt:
/// ~~~~{.cpp}
/// using namespace ROOT::VecOps;
/// RVec<float> v {1.f, 2.f, 4.f};
/// auto v_argmax = ArgMax(v);
/// v_argmax
/// // (int) 2
/// ~~~~
template <typename T>
std::size_t ArgMax(const RVec<T> &v)
{
   return std::distance(v.begin(), std::max_element(v.begin(), v.end()));
}

/// Get the index of the smallest element of an RVec
/// In case of multiple occurrences of the minimum values,
/// the index corresponding to the first occurrence is returned.
///
/// Example code, at the ROOT prompt:
/// ~~~~{.cpp}
/// using namespace ROOT::VecOps;
/// RVec<float> v {1.f, 2.f, 4.f};
/// auto v_argmin = ArgMin(v);
/// v_argmin
/// // (int) 0
/// ~~~~
template <typename T>
std::size_t ArgMin(const RVec<T> &v)
{
   return std::distance(v.begin(), std::min_element(v.begin(), v.end()));
}

/// Get the variance of the elements of an RVec
///
/// The return type is a double precision floating point number.
/// Example code, at the ROOT prompt:
/// ~~~{.cpp}
/// using namespace ROOT::VecOps;
/// RVec<float> v {1.f, 2.f, 4.f};
/// auto v_var = Var(v);
/// v_var
/// // (double) 2.3333333
/// ~~~
template <typename T>
double Var(const RVec<T> &v)
{
   const std::size_t size = v.size();
   if (size < std::size_t(2)) return 0.;
   T sum_squares(0), squared_sum(0);
   auto pred = [&sum_squares, &squared_sum](const T& x) {sum_squares+=x*x; squared_sum+=x;};
   std::for_each(v.begin(), v.end(), pred);
   squared_sum *= squared_sum;
   const auto dsize = (double) size;
   return 1. / (dsize - 1.) * (sum_squares - squared_sum / dsize );
}

/// Get the standard deviation of the elements of an RVec
///
/// The return type is a double precision floating point number.
/// Example code, at the ROOT prompt:
/// ~~~{.cpp}
/// using namespace ROOT::VecOps;
/// RVec<float> v {1.f, 2.f, 4.f};
/// auto v_sd = StdDev(v);
/// v_sd
/// // (double) 1.5275252
/// ~~~
template <typename T>
double StdDev(const RVec<T> &v)
{
   return std::sqrt(Var(v));
}

/// Create new collection applying a callable to the elements of the input collection
///
/// Example code, at the ROOT prompt:
/// ~~~{.cpp}
/// using namespace ROOT::VecOps;
/// RVec<float> v {1.f, 2.f, 4.f};
/// auto v_square = Map(v, [](float f){return f* 2.f;});
/// v_square
/// // (ROOT::VecOps::RVec<float> &) { 2.00000f, 4.00000f, 8.00000f }
///
/// RVec<float> x({1.f, 2.f, 3.f});
/// RVec<float> y({4.f, 5.f, 6.f});
/// RVec<float> z({7.f, 8.f, 9.f});
/// auto mod = [](float x, float y, float z) { return sqrt(x * x + y * y + z * z); };
/// auto v_mod = Map(x, y, z, mod);
/// v_mod
/// // (ROOT::VecOps::RVec<float> &) { 8.12404f, 9.64365f, 11.2250f }
/// ~~~
template <typename... Args>
auto Map(Args &&... args)
   -> decltype(ROOT::Detail::VecOps::MapFromTuple(std::forward_as_tuple(args...),
                                                  std::make_index_sequence<sizeof...(args) - 1>()))
{
   /*
   Here the strategy in order to generalise the previous implementation of Map, i.e.
   `RVec Map(RVec, F)`, here we need to move the last parameter of the pack in first
   position in order to be able to invoke the Map function with automatic type deduction.
   This is achieved in two steps:
   1. Forward as tuple the pack to MapFromTuple
   2. Invoke the MapImpl helper which has the signature `template<...T, F> RVec MapImpl(F &&f, RVec<T>...)`
   NOTA BENE: the signature is very heavy but it is one of the lightest ways to manage in C++11
   to build the return type based on the template args.
   */
   return ROOT::Detail::VecOps::MapFromTuple(std::forward_as_tuple(args...),
                                             std::make_index_sequence<sizeof...(args) - 1>());
}

/// Create a new collection with the elements passing the filter expressed by the predicate
///
/// Example code, at the ROOT prompt:
/// ~~~{.cpp}
/// using namespace ROOT::VecOps;
/// RVec<int> v {1, 2, 4};
/// auto v_even = Filter(v, [](int i){return 0 == i%2;});
/// v_even
/// // (ROOT::VecOps::RVec<int> &) { 2, 4 }
/// ~~~
template <typename T, typename F>
RVec<T> Filter(const RVec<T> &v, F &&f)
{
   const auto thisSize = v.size();
   RVec<T> w;
   w.reserve(thisSize);
   for (auto &&val : v) {
      if (f(val))
         w.emplace_back(val);
   }
   return w;
}

/// Return true if any of the elements equates to true, return false otherwise.
///
/// Example code, at the ROOT prompt:
/// ~~~{.cpp}
/// using namespace ROOT::VecOps;
/// RVec<int> v {0, 1, 0};
/// auto anyTrue = Any(v);
/// anyTrue
/// // (bool) true
/// ~~~
template <typename T>
auto Any(const RVec<T> &v) -> decltype(v[0] == true)
{
   for (auto &&e : v)
      if (static_cast<bool>(e) == true)
         return true;
   return false;
}

/// Return true if all of the elements equate to true, return false otherwise.
///
/// Example code, at the ROOT prompt:
/// ~~~{.cpp}
/// using namespace ROOT::VecOps;
/// RVec<int> v {0, 1, 0};
/// auto allTrue = All(v);
/// allTrue
/// // (bool) false
/// ~~~
template <typename T>
auto All(const RVec<T> &v) -> decltype(v[0] == false)
{
   for (auto &&e : v)
      if (static_cast<bool>(e) == false)
         return false;
   return true;
}

template <typename T>
void swap(RVec<T> &lhs, RVec<T> &rhs)
{
   lhs.swap(rhs);
}

/// Return an RVec of indices that sort the input RVec
///
/// Example code, at the ROOT prompt:
/// ~~~{.cpp}
/// using namespace ROOT::VecOps;
/// RVec<double> v {2., 3., 1.};
/// auto sortIndices = Argsort(v);
/// // (ROOT::VecOps::RVec<unsigned long> &) { 2, 0, 1 }
/// auto values = Take(v, sortIndices)
/// // (ROOT::VecOps::RVec<double> &) { 1., 2., 3. }
/// ~~~
template <typename T>
RVec<typename RVec<T>::size_type> Argsort(const RVec<T> &v)
{
   using size_type = typename RVec<T>::size_type;
   RVec<size_type> i(v.size());
   std::iota(i.begin(), i.end(), 0);
   std::sort(i.begin(), i.end(), [&v](size_type i1, size_type i2) { return v[i1] < v[i2]; });
   return i;
}

/// Return an RVec of indices that sort the input RVec based on a comparison function.
///
/// Example code, at the ROOT prompt:
/// ~~~{.cpp}
/// using namespace ROOT::VecOps;
/// RVec<double> v {2., 3., 1.};
/// auto sortIndices = Argsort(v, [](double x, double y) {return x > y;})
/// // (ROOT::VecOps::RVec<unsigned long> &) { 1, 0, 2 }
/// auto values = Take(v, sortIndices)
/// // (ROOT::VecOps::RVec<double> &) { 3., 2., 1. }
/// ~~~
template <typename T, typename Compare>
RVec<typename RVec<T>::size_type> Argsort(const RVec<T> &v, Compare &&c)
{
   using size_type = typename RVec<T>::size_type;
   RVec<size_type> i(v.size());
   std::iota(i.begin(), i.end(), 0);
   std::sort(i.begin(), i.end(),
             [&v, &c](size_type i1, size_type i2) { return c(v[i1], v[i2]); });
   return i;
}

/// Return elements of a vector at given indices
///
/// Example code, at the ROOT prompt:
/// ~~~{.cpp}
/// using namespace ROOT::VecOps;
/// RVec<double> v {2., 3., 1.};
/// auto vTaken = Take(v, {0,2});
/// vTaken
/// // (ROOT::VecOps::RVec<double>) { 2.0000000, 1.0000000 }
/// ~~~
template <typename T>
RVec<T> Take(const RVec<T> &v, const RVec<typename RVec<T>::size_type> &i)
{
   using size_type = typename RVec<T>::size_type;
   const size_type isize = i.size();
   RVec<T> r(isize);
   for (size_type k = 0; k < isize; k++)
      r[k] = v[i[k]];
   return r;
}

/// Return first or last `n` elements of an RVec
///
/// if `n > 0` and last elements if `n < 0`.
///
/// Example code, at the ROOT prompt:
/// ~~~{.cpp}
/// using namespace ROOT::VecOps;
/// RVec<double> v {2., 3., 1.};
/// auto firstTwo = Take(v, 2);
/// firstTwo
/// // (ROOT::VecOps::RVec<double>) { 2.0000000, 3.0000000 }
/// auto lastOne = Take(v, -1);
/// lastOne
/// // (ROOT::VecOps::RVec<double>) { 1.0000000 }
/// ~~~
template <typename T>
RVec<T> Take(const RVec<T> &v, const int n)
{
   using size_type = typename RVec<T>::size_type;
   const size_type size = v.size();
   const size_type absn = std::abs(n);
   if (absn > size) {
      std::stringstream ss;
      ss << "Try to take " << absn << " elements but vector has only size " << size << ".";
      throw std::runtime_error(ss.str());
   }
   RVec<T> r(absn);
   if (n < 0) {
      for (size_type k = 0; k < absn; k++)
         r[k] = v[size - absn + k];
   } else {
      for (size_type k = 0; k < absn; k++)
         r[k] = v[k];
   }
   return r;
}

/// Return copy of reversed vector
///
/// Example code, at the ROOT prompt:
/// ~~~{.cpp}
/// using namespace ROOT::VecOps;
/// RVec<double> v {2., 3., 1.};
/// auto v_reverse = Reverse(v);
/// v_reverse
/// // (ROOT::VecOps::RVec<double>) { 1.0000000, 3.0000000, 2.0000000 }
/// ~~~
template <typename T>
RVec<T> Reverse(const RVec<T> &v)
{
   RVec<T> r(v);
   std::reverse(r.begin(), r.end());
   return r;
}

/// Return copy of RVec with elements sorted in ascending order
///
/// This helper is different from ArgSort since it does not return an RVec of indices,
/// but an RVec of values.
///
/// Example code, at the ROOT prompt:
/// ~~~{.cpp}
/// using namespace ROOT::VecOps;
/// RVec<double> v {2., 3., 1.};
/// auto v_sorted = Sort(v);
/// v_sorted
/// // (ROOT::VecOps::RVec<double>) { 1.0000000, 2.0000000, 3.0000000 }
/// ~~~
template <typename T>
RVec<T> Sort(const RVec<T> &v)
{
   RVec<T> r(v);
   std::sort(r.begin(), r.end());
   return r;
}

/// Return copy of RVec with elements sorted based on a comparison operator
///
/// The comparison operator has to fullfill the same requirements of the
/// predicate of by std::sort.
///
///
/// This helper is different from ArgSort since it does not return an RVec of indices,
/// but an RVec of values.
///
/// Example code, at the ROOT prompt:
/// ~~~{.cpp}
/// using namespace ROOT::VecOps;
/// RVec<double> v {2., 3., 1.};
/// auto v_sorted = Sort(v, [](double x, double y) {return 1/x < 1/y;});
/// v_sorted
/// // (ROOT::VecOps::RVec<double>) { 3.0000000, 2.0000000, 1.0000000 }
/// ~~~
template <typename T, typename Compare>
RVec<T> Sort(const RVec<T> &v, Compare &&c)
{
   RVec<T> r(v);
   std::sort(r.begin(), r.end(), std::forward<Compare>(c));
   return r;
}

/// Return the indices that represent all combinations of the elements of two
/// RVecs.
///
/// The type of the return value is an RVec of two RVecs containing indices.
///
/// Example code, at the ROOT prompt:
/// ~~~{.cpp}
/// using namespace ROOT::VecOps;
/// auto comb_idx = Combinations(3, 2);
/// comb_idx
/// // (ROOT::VecOps::RVec<ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>::size_type> >) { { 0, 0, 1, 1, 2, 2 }, { 0, 1, 0, 1, 0, 1 } }
/// ~~~
inline RVec<RVec<std::size_t>> Combinations(const std::size_t size1, const std::size_t size2)
{
   using size_type = std::size_t;
   RVec<RVec<size_type>> r(2);
   r[0].resize(size1*size2);
   r[1].resize(size1*size2);
   size_type c = 0;
   for(size_type i=0; i<size1; i++) {
      for(size_type j=0; j<size2; j++) {
         r[0][c] = i;
         r[1][c] = j;
         c++;
      }
   }
   return r;
}

/// Return the indices that represent all combinations of the elements of two
/// RVecs.
///
/// The type of the return value is an RVec of two RVecs containing indices.
///
/// Example code, at the ROOT prompt:
/// ~~~{.cpp}
/// using namespace ROOT::VecOps;
/// RVec<double> v1 {1., 2., 3.};
/// RVec<double> v2 {-4., -5.};
/// auto comb_idx = Combinations(v1, v2);
/// comb_idx
/// // (ROOT::VecOps::RVec<ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>::size_type> >) { { 0, 0, 1, 1, 2, 2 }, { 0, 1,
/// 0, 1, 0, 1 } }
/// ~~~
template <typename T1, typename T2>
RVec<RVec<typename RVec<T1>::size_type>> Combinations(const RVec<T1> &v1, const RVec<T2> &v2)
{
   return Combinations(v1.size(), v2.size());
}

/// Return the indices that represent all unique combinations of the
/// elements of a given RVec.
///
/// ~~~{.cpp}
/// using namespace ROOT::VecOps;
/// RVec<double> v {1., 2., 3., 4.};
/// auto v_1 = Combinations(v, 1);
/// v_1
/// (ROOT::VecOps::RVec<ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>::size_type> >) { { 0, 1, 2, 3 } }
/// auto v_2 = Combinations(v, 2);
/// auto v_2
/// (ROOT::VecOps::RVec<ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>::size_type> >) { { 0, 0, 0, 1, 1, 2 }, { 1, 2, 3, 2, 3, 3 } }
/// auto v_3 = Combinations(v, 3);
/// v_3
/// (ROOT::VecOps::RVec<ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>::size_type> >) { { 0, 0, 0, 1 }, { 1, 1, 2, 2 }, { 2, 3, 3, 3 } }
/// auto v_4 = Combinations(v, 4);
/// v_4
/// (ROOT::VecOps::RVec<ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>::size_type> >) { { 0 }, { 1 }, { 2 }, { 3 } }
/// ~~~
template <typename T>
RVec<RVec<typename RVec<T>::size_type>> Combinations(const RVec<T>& v, const typename RVec<T>::size_type n)
{
   using size_type = typename RVec<T>::size_type;
   const size_type s = v.size();
   if (n > s) {
      std::stringstream ss;
      ss << "Cannot make unique combinations of size " << n << " from vector of size " << s << ".";
      throw std::runtime_error(ss.str());
   }
   RVec<size_type> indices(s);
   for(size_type k=0; k<s; k++)
      indices[k] = k;
   RVec<RVec<size_type>> c(n);
   for(size_type k=0; k<n; k++)
      c[k].emplace_back(indices[k]);
   while (true) {
      bool run_through = true;
      long i = n - 1;
      for (; i>=0; i--) {
         if (indices[i] != i + s - n){
            run_through = false;
            break;
         }
      }
      if (run_through) {
         return c;
      }
      indices[i]++;
      for (long j=i+1; j<(long)n; j++)
         indices[j] = indices[j-1] + 1;
      for(size_type k=0; k<n; k++)
         c[k].emplace_back(indices[k]);
   }
}

/// Return the indices of the elements which are not zero
///
/// Example code, at the ROOT prompt:
/// ~~~{.cpp}
/// using namespace ROOT::VecOps;
/// RVec<double> v {2., 0., 3., 0., 1.};
/// auto nonzero_idx = Nonzero(v);
/// nonzero_idx
/// // (ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>::size_type>) { 0, 2, 4 }
/// ~~~
template <typename T>
RVec<typename RVec<T>::size_type> Nonzero(const RVec<T> &v)
{
   using size_type = typename RVec<T>::size_type;
   RVec<size_type> r;
   const auto size = v.size();
   r.reserve(size);
   for(size_type i=0; i<size; i++) {
      if(v[i] != 0) {
         r.emplace_back(i);
      }
   }
   return r;
}

/// Return the intersection of elements of two RVecs.
///
/// Each element of v1 is looked up in v2 and added to the returned vector if
/// found. Following, the order of v1 is preserved. If v2 is already sorted, the
/// optional argument v2_is_sorted can be used to toggle of the internal sorting
/// step, therewith optimising runtime.
///
/// Example code, at the ROOT prompt:
/// ~~~{.cpp}
/// using namespace ROOT::VecOps;
/// RVec<double> v1 {1., 2., 3.};
/// RVec<double> v2 {-4., -5., 2., 1.};
/// auto v1_intersect_v2 = Intersect(v1, v2);
/// v1_intersect_v2
/// // (ROOT::VecOps::RVec<double>) { 1.0000000, 2.0000000 }
/// ~~~
template <typename T>
RVec<T> Intersect(const RVec<T>& v1, const RVec<T>& v2, bool v2_is_sorted = false)
{
   RVec<T> v2_sorted;
   if (!v2_is_sorted) v2_sorted = Sort(v2);
   const auto v2_begin = v2_is_sorted ? v2.begin() : v2_sorted.begin();
   const auto v2_end = v2_is_sorted ? v2.end() : v2_sorted.end();
   RVec<T> r;
   const auto size = v1.size();
   r.reserve(size);
   using size_type = typename RVec<T>::size_type;
   for(size_type i=0; i<size; i++) {
      if (std::binary_search(v2_begin, v2_end, v1[i])) {
         r.emplace_back(v1[i]);
      }
   }
   return r;
}

/// Return the elements of v1 if the condition c is true and v2 if the
/// condition c is false.
///
/// Example code, at the ROOT prompt:
/// ~~~{.cpp}
/// using namespace ROOT::VecOps;
/// RVec<double> v1 {1., 2., 3.};
/// RVec<double> v2 {-1., -2., -3.};
/// auto c = v1 > 1;
/// c
/// // (ROOT::VecOps::RVec<int> &) { 0, 1, 1 }
/// auto if_c_v1_else_v2 = Where(c, v1, v2);
/// if_c_v1_else_v2
/// // (ROOT::VecOps::RVec<double> &) { -1.0000000, 2.0000000, 3.0000000 }
/// ~~~
template <typename T>
RVec<T> Where(const RVec<int>& c, const RVec<T>& v1, const RVec<T>& v2)
{
   using size_type = typename RVec<T>::size_type;
   const size_type size = c.size();
   RVec<T> r;
   r.reserve(size);
   for (size_type i=0; i<size; i++) {
      r.emplace_back(c[i] != 0 ? v1[i] : v2[i]);
   }
   return r;
}

/// Return the elements of v1 if the condition c is true and sets the value v2
/// if the condition c is false.
///
/// Example code, at the ROOT prompt:
/// ~~~{.cpp}
/// using namespace ROOT::VecOps;
/// RVec<double> v1 {1., 2., 3.};
/// double v2 = 4.;
/// auto c = v1 > 1;
/// c
/// // (ROOT::VecOps::RVec<int> &) { 0, 1, 1 }
/// auto if_c_v1_else_v2 = Where(c, v1, v2);
/// if_c_v1_else_v2
/// // (ROOT::VecOps::RVec<double>) { 4.0000000, 2.0000000, 3.0000000 }
/// ~~~
template <typename T>
RVec<T> Where(const RVec<int>& c, const RVec<T>& v1, T v2)
{
   using size_type = typename RVec<T>::size_type;
   const size_type size = c.size();
   RVec<T> r;
   r.reserve(size);
   for (size_type i=0; i<size; i++) {
      r.emplace_back(c[i] != 0 ? v1[i] : v2);
   }
   return r;
}

/// Return the elements of v2 if the condition c is false and sets the value v1
/// if the condition c is true.
///
/// Example code, at the ROOT prompt:
/// ~~~{.cpp}
/// using namespace ROOT::VecOps;
/// double v1 = 4.;
/// RVec<double> v2 {1., 2., 3.};
/// auto c = v2 > 1;
/// c
/// // (ROOT::VecOps::RVec<int> &) { 0, 1, 1 }
/// auto if_c_v1_else_v2 = Where(c, v1, v2);
/// if_c_v1_else_v2
/// // (ROOT::VecOps::RVec<double>) { 1.0000000, 4.0000000, 4.0000000 }
/// ~~~
template <typename T>
RVec<T> Where(const RVec<int>& c, T v1, const RVec<T>& v2)
{
   using size_type = typename RVec<T>::size_type;
   const size_type size = c.size();
   RVec<T> r;
   r.reserve(size);
   for (size_type i=0; i<size; i++) {
      r.emplace_back(c[i] != 0 ? v1 : v2[i]);
   }
   return r;
}

/// Return a vector with the value v2 if the condition c is false and sets the
/// value v1 if the condition c is true.
///
/// Example code, at the ROOT prompt:
/// ~~~{.cpp}
/// using namespace ROOT::VecOps;
/// double v1 = 4.;
/// double v2 = 2.;
/// RVec<int> c {0, 1, 1};
/// auto if_c_v1_else_v2 = Where(c, v1, v2);
/// if_c_v1_else_v2
/// // (ROOT::VecOps::RVec<double>) { 2.0000000, 4.0000000, 4.0000000 }
/// ~~~
template <typename T>
RVec<T> Where(const RVec<int>& c, T v1, T v2)
{
   using size_type = typename RVec<T>::size_type;
   const size_type size = c.size();
   RVec<T> r;
   r.reserve(size);
   for (size_type i=0; i<size; i++) {
      r.emplace_back(c[i] != 0 ? v1 : v2);
   }
   return r;
}

/// Return the concatenation of two RVecs.
///
/// Example code, at the ROOT prompt:
/// ~~~{.cpp}
/// using namespace ROOT::VecOps;
/// RVec<float> rvf {0.f, 1.f, 2.f};
/// RVec<int> rvi {7, 8, 9};
/// Concatenate(rvf, rvi);
/// // (ROOT::VecOps::RVec<float>) { 2.0000000, 4.0000000, 4.0000000 }
/// ~~~
template <typename T0, typename T1, typename Common_t = typename std::common_type<T0, T1>::type>
RVec<Common_t> Concatenate(const RVec<T0> &v0, const RVec<T1> &v1)
{
   RVec<Common_t> res;
   res.reserve(v0.size() + v1.size());
   std::copy(v0.begin(), v0.end(), std::back_inserter(res));
   std::copy(v1.begin(), v1.end(), std::back_inserter(res));
   return res;
}

/// Return the angle difference \f$\Delta \phi\f$ of two scalars.
///
/// The function computes the closest angle from v1 to v2 with sign and is
/// therefore in the range \f$[-\pi, \pi]\f$.
/// The computation is done per default in radians \f$c = \pi\f$ but can be switched
/// to degrees \f$c = 180\f$.
template <typename T>
T DeltaPhi(T v1, T v2, const T c = M_PI)
{
   static_assert(std::is_floating_point<T>::value,
                 "DeltaPhi must be called with floating point values.");
   auto r = std::fmod(v2 - v1, 2.0 * c);
   if (r < -c) {
      r += 2.0 * c;
   }
   else if (r > c) {
      r -= 2.0 * c;
   }
   return r;
}

/// Return the angle difference \f$\Delta \phi\f$ in radians of two vectors.
///
/// The function computes the closest angle from v1 to v2 with sign and is
/// therefore in the range \f$[-\pi, \pi]\f$.
/// The computation is done per default in radians \f$c = \pi\f$ but can be switched
/// to degrees \f$c = 180\f$.
template <typename T>
RVec<T> DeltaPhi(const RVec<T>& v1, const RVec<T>& v2, const T c = M_PI)
{
   using size_type = typename RVec<T>::size_type;
   const size_type size = v1.size();
   auto r = RVec<T>(size);
   for (size_type i = 0; i < size; i++) {
      r[i] = DeltaPhi(v1[i], v2[i], c);
   }
   return r;
}

/// Return the angle difference \f$\Delta \phi\f$ in radians of a vector and a scalar.
///
/// The function computes the closest angle from v1 to v2 with sign and is
/// therefore in the range \f$[-\pi, \pi]\f$.
/// The computation is done per default in radians \f$c = \pi\f$ but can be switched
/// to degrees \f$c = 180\f$.
template <typename T>
RVec<T> DeltaPhi(const RVec<T>& v1, T v2, const T c = M_PI)
{
   using size_type = typename RVec<T>::size_type;
   const size_type size = v1.size();
   auto r = RVec<T>(size);
   for (size_type i = 0; i < size; i++) {
      r[i] = DeltaPhi(v1[i], v2, c);
   }
   return r;
}

/// Return the angle difference \f$\Delta \phi\f$ in radians of a scalar and a vector.
///
/// The function computes the closest angle from v1 to v2 with sign and is
/// therefore in the range \f$[-\pi, \pi]\f$.
/// The computation is done per default in radians \f$c = \pi\f$ but can be switched
/// to degrees \f$c = 180\f$.
template <typename T>
RVec<T> DeltaPhi(T v1, const RVec<T>& v2, const T c = M_PI)
{
   using size_type = typename RVec<T>::size_type;
   const size_type size = v2.size();
   auto r = RVec<T>(size);
   for (size_type i = 0; i < size; i++) {
      r[i] = DeltaPhi(v1, v2[i], c);
   }
   return r;
}

/// Return the square of the distance on the \f$\eta\f$-\f$\phi\f$ plane (\f$\Delta R\f$) from
/// the collections eta1, eta2, phi1 and phi2.
///
/// The function computes \f$\Delta R^2 = (\eta_1 - \eta_2)^2 + (\phi_1 - \phi_2)^2\f$
/// of the given collections eta1, eta2, phi1 and phi2. The angle \f$\phi\f$ can
/// be set to radian or degrees using the optional argument c, see the documentation
/// of the DeltaPhi helper.
template <typename T>
RVec<T> DeltaR2(const RVec<T>& eta1, const RVec<T>& eta2, const RVec<T>& phi1, const RVec<T>& phi2, const T c = M_PI)
{
   const auto dphi = DeltaPhi(phi1, phi2, c);
   return (eta1 - eta2) * (eta1 - eta2) + dphi * dphi;
}

/// Return the distance on the \f$\eta\f$-\f$\phi\f$ plane (\f$\Delta R\f$) from
/// the collections eta1, eta2, phi1 and phi2.
///
/// The function computes \f$\Delta R = \sqrt{(\eta_1 - \eta_2)^2 + (\phi_1 - \phi_2)^2}\f$
/// of the given collections eta1, eta2, phi1 and phi2. The angle \f$\phi\f$ can
/// be set to radian or degrees using the optional argument c, see the documentation
/// of the DeltaPhi helper.
template <typename T>
RVec<T> DeltaR(const RVec<T>& eta1, const RVec<T>& eta2, const RVec<T>& phi1, const RVec<T>& phi2, const T c = M_PI)
{
   return sqrt(DeltaR2(eta1, eta2, phi1, phi2, c));
}

/// Return the distance on the \f$\eta\f$-\f$\phi\f$ plane (\f$\Delta R\f$) from
/// the scalars eta1, eta2, phi1 and phi2.
///
/// The function computes \f$\Delta R = \sqrt{(\eta_1 - \eta_2)^2 + (\phi_1 - \phi_2)^2}\f$
/// of the given scalars eta1, eta2, phi1 and phi2. The angle \f$\phi\f$ can
/// be set to radian or degrees using the optional argument c, see the documentation
/// of the DeltaPhi helper.
template <typename T>
T DeltaR(T eta1, T eta2, T phi1, T phi2, const T c = M_PI)
{
   const auto dphi = DeltaPhi(phi1, phi2, c);
   return std::sqrt((eta1 - eta2) * (eta1 - eta2) + dphi * dphi);
}

/// Return the invariant mass of two particles given the collections of the quantities
/// transverse momentum (pt), rapidity (eta), azimuth (phi) and mass.
///
/// The function computes the invariant mass of two particles with the four-vectors
/// (pt1, eta2, phi1, mass1) and (pt2, eta2, phi2, mass2).
template <typename T>
RVec<T> InvariantMasses(
        const RVec<T>& pt1, const RVec<T>& eta1, const RVec<T>& phi1, const RVec<T>& mass1,
        const RVec<T>& pt2, const RVec<T>& eta2, const RVec<T>& phi2, const RVec<T>& mass2)
{
   std::size_t size = pt1.size();

   R__ASSERT(eta1.size() == size && phi1.size() == size && mass1.size() == size);
   R__ASSERT(pt2.size() == size && phi2.size() == size && mass2.size() == size);

   RVec<T> inv_masses(size);

   for (std::size_t i = 0u; i < size; ++i) {
      // Conversion from (pt, eta, phi, mass) to (x, y, z, e) coordinate system
      const auto x1 = pt1[i] * std::cos(phi1[i]);
      const auto y1 = pt1[i] * std::sin(phi1[i]);
      const auto z1 = pt1[i] * std::sinh(eta1[i]);
      const auto e1 = std::sqrt(x1 * x1 + y1 * y1 + z1 * z1 + mass1[i] * mass1[i]);

      const auto x2 = pt2[i] * std::cos(phi2[i]);
      const auto y2 = pt2[i] * std::sin(phi2[i]);
      const auto z2 = pt2[i] * std::sinh(eta2[i]);
      const auto e2 = std::sqrt(x2 * x2 + y2 * y2 + z2 * z2 + mass2[i] * mass2[i]);

      // Addition of particle four-vector elements
      const auto e = e1 + e2;
      const auto x = x1 + x2;
      const auto y = y1 + y2;
      const auto z = z1 + z2;

      inv_masses[i] = std::sqrt(e * e - x * x - y * y - z * z);
   }

   // Return invariant mass with (+, -, -, -) metric
   return inv_masses;
}

/// Return the invariant mass of multiple particles given the collections of the
/// quantities transverse momentum (pt), rapidity (eta), azimuth (phi) and mass.
///
/// The function computes the invariant mass of multiple particles with the
/// four-vectors (pt, eta, phi, mass).
template <typename T>
T InvariantMass(const RVec<T>& pt, const RVec<T>& eta, const RVec<T>& phi, const RVec<T>& mass)
{
   const std::size_t size = pt.size();

   R__ASSERT(eta.size() == size && phi.size() == size && mass.size() == size);

   T x_sum = 0.;
   T y_sum = 0.;
   T z_sum = 0.;
   T e_sum = 0.;

   for (std::size_t i = 0u; i < size; ++ i) {
      // Convert to (e, x, y, z) coordinate system and update sums
      const auto x = pt[i] * std::cos(phi[i]);
      x_sum += x;
      const auto y = pt[i] * std::sin(phi[i]);
      y_sum += y;
      const auto z = pt[i] * std::sinh(eta[i]);
      z_sum += z;
      const auto e = std::sqrt(x * x + y * y + z * z + mass[i] * mass[i]);
      e_sum += e;
   }

   // Return invariant mass with (+, -, -, -) metric
   return std::sqrt(e_sum * e_sum - x_sum * x_sum - y_sum * y_sum - z_sum * z_sum);
}

////////////////////////////////////////////////////////////////////////////
/// \brief Build an RVec of objects starting from RVecs of input to their constructors.
/// \tparam T Type of the objects contained in the created RVec.
/// \tparam Args_t Pack of types templating the input RVecs.
/// \param[in] args The RVecs containing the values used to initialise the output objects.
/// \return The RVec of objects initialised with the input parameters.
///
/// Example code, at the ROOT prompt:
/// ~~~{.cpp}
/// using namespace ROOT::VecOps;
/// RVec<float> pts = {15.5, 34.32, 12.95};
/// RVec<float> etas = {0.3, 2.2, 1.32};
/// RVec<float> phis = {0.1, 3.02, 2.2};
/// RVec<float> masses = {105.65, 105.65, 105.65};
/// auto fourVecs = Construct<ROOT::Math::PtEtaPhiMVector>(pts, etas, phis, masses);
/// cout << fourVecs << endl;
/// // { (15.5,0.3,0.1,105.65), (34.32,2.2,3.02,105.65), (12.95,1.32,2.2,105.65) }
/// ~~~
template <typename T, typename... Args_t>
RVec<T> Construct(const RVec<Args_t> &... args)
{
   const auto size = ::ROOT::Detail::VecOps::GetVectorsSize("Construct", args...);
   RVec<T> ret;
   ret.reserve(size);
   for (auto i = 0UL; i < size; ++i) {
      ret.emplace_back(args[i]...);
   }
   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// Print a RVec at the prompt:
template <class T>
std::ostream &operator<<(std::ostream &os, const RVec<T> &v)
{
   // In order to print properly, convert to 64 bit int if this is a char
   constexpr bool mustConvert = std::is_same<char, T>::value || std::is_same<signed char, T>::value ||
                                std::is_same<unsigned char, T>::value || std::is_same<wchar_t, T>::value ||
                                std::is_same<char16_t, T>::value || std::is_same<char32_t, T>::value;
   using Print_t = typename std::conditional<mustConvert, long long int, T>::type;
   os << "{ ";
   auto size = v.size();
   if (size) {
      for (std::size_t i = 0; i < size - 1; ++i) {
         os << (Print_t)v[i] << ", ";
      }
      os << (Print_t)v[size - 1];
   }
   os << " }";
   return os;
}

#if (_VECOPS_USE_EXTERN_TEMPLATES)

#define RVEC_EXTERN_UNARY_OPERATOR(T, OP) \
   extern template RVec<T> operator OP<T>(const RVec<T> &);

#define RVEC_EXTERN_BINARY_OPERATOR(T, OP)                                     \
   extern template auto operator OP<T, T>(const T &x, const RVec<T> &v)        \
      -> RVec<decltype(x OP v[0])>;                                            \
   extern template auto operator OP<T, T>(const RVec<T> &v, const T &y)        \
      -> RVec<decltype(v[0] OP y)>;                                            \
   extern template auto operator OP<T, T>(const RVec<T> &v0, const RVec<T> &v1)\
      -> RVec<decltype(v0[0] OP v1[0])>;

#define RVEC_EXTERN_ASSIGN_OPERATOR(T, OP)                           \
   extern template RVec<T> &operator OP<T, T>(RVec<T> &, const T &); \
   extern template RVec<T> &operator OP<T, T>(RVec<T> &, const RVec<T> &);

#define RVEC_EXTERN_LOGICAL_OPERATOR(T, OP)                                 \
   extern template RVec<int> operator OP<T, T>(const RVec<T> &, const T &); \
   extern template RVec<int> operator OP<T, T>(const T &, const RVec<T> &); \
   extern template RVec<int> operator OP<T, T>(const RVec<T> &, const RVec<T> &);

#define RVEC_EXTERN_FLOAT_TEMPLATE(T)   \
   extern template class RVec<T>;       \
   RVEC_EXTERN_UNARY_OPERATOR(T, +)     \
   RVEC_EXTERN_UNARY_OPERATOR(T, -)     \
   RVEC_EXTERN_UNARY_OPERATOR(T, !)     \
   RVEC_EXTERN_BINARY_OPERATOR(T, +)    \
   RVEC_EXTERN_BINARY_OPERATOR(T, -)    \
   RVEC_EXTERN_BINARY_OPERATOR(T, *)    \
   RVEC_EXTERN_BINARY_OPERATOR(T, /)    \
   RVEC_EXTERN_ASSIGN_OPERATOR(T, +=)   \
   RVEC_EXTERN_ASSIGN_OPERATOR(T, -=)   \
   RVEC_EXTERN_ASSIGN_OPERATOR(T, *=)   \
   RVEC_EXTERN_ASSIGN_OPERATOR(T, /=)   \
   RVEC_EXTERN_LOGICAL_OPERATOR(T, <)   \
   RVEC_EXTERN_LOGICAL_OPERATOR(T, >)   \
   RVEC_EXTERN_LOGICAL_OPERATOR(T, ==)  \
   RVEC_EXTERN_LOGICAL_OPERATOR(T, !=)  \
   RVEC_EXTERN_LOGICAL_OPERATOR(T, <=)  \
   RVEC_EXTERN_LOGICAL_OPERATOR(T, >=)  \
   RVEC_EXTERN_LOGICAL_OPERATOR(T, &&)  \
   RVEC_EXTERN_LOGICAL_OPERATOR(T, ||)

#define RVEC_EXTERN_INTEGER_TEMPLATE(T) \
   extern template class RVec<T>;       \
   RVEC_EXTERN_UNARY_OPERATOR(T, +)     \
   RVEC_EXTERN_UNARY_OPERATOR(T, -)     \
   RVEC_EXTERN_UNARY_OPERATOR(T, ~)     \
   RVEC_EXTERN_UNARY_OPERATOR(T, !)     \
   RVEC_EXTERN_BINARY_OPERATOR(T, +)    \
   RVEC_EXTERN_BINARY_OPERATOR(T, -)    \
   RVEC_EXTERN_BINARY_OPERATOR(T, *)    \
   RVEC_EXTERN_BINARY_OPERATOR(T, /)    \
   RVEC_EXTERN_BINARY_OPERATOR(T, %)    \
   RVEC_EXTERN_BINARY_OPERATOR(T, &)    \
   RVEC_EXTERN_BINARY_OPERATOR(T, |)    \
   RVEC_EXTERN_BINARY_OPERATOR(T, ^)    \
   RVEC_EXTERN_ASSIGN_OPERATOR(T, +=)   \
   RVEC_EXTERN_ASSIGN_OPERATOR(T, -=)   \
   RVEC_EXTERN_ASSIGN_OPERATOR(T, *=)   \
   RVEC_EXTERN_ASSIGN_OPERATOR(T, /=)   \
   RVEC_EXTERN_ASSIGN_OPERATOR(T, %=)   \
   RVEC_EXTERN_ASSIGN_OPERATOR(T, &=)   \
   RVEC_EXTERN_ASSIGN_OPERATOR(T, |=)   \
   RVEC_EXTERN_ASSIGN_OPERATOR(T, ^=)   \
   RVEC_EXTERN_ASSIGN_OPERATOR(T, >>=)  \
   RVEC_EXTERN_ASSIGN_OPERATOR(T, <<=)  \
   RVEC_EXTERN_LOGICAL_OPERATOR(T, <)   \
   RVEC_EXTERN_LOGICAL_OPERATOR(T, >)   \
   RVEC_EXTERN_LOGICAL_OPERATOR(T, ==)  \
   RVEC_EXTERN_LOGICAL_OPERATOR(T, !=)  \
   RVEC_EXTERN_LOGICAL_OPERATOR(T, <=)  \
   RVEC_EXTERN_LOGICAL_OPERATOR(T, >=)  \
   RVEC_EXTERN_LOGICAL_OPERATOR(T, &&)  \
   RVEC_EXTERN_LOGICAL_OPERATOR(T, ||)

RVEC_EXTERN_INTEGER_TEMPLATE(char)
RVEC_EXTERN_INTEGER_TEMPLATE(short)
RVEC_EXTERN_INTEGER_TEMPLATE(int)
RVEC_EXTERN_INTEGER_TEMPLATE(long)
//RVEC_EXTERN_INTEGER_TEMPLATE(long long)

RVEC_EXTERN_INTEGER_TEMPLATE(unsigned char)
RVEC_EXTERN_INTEGER_TEMPLATE(unsigned short)
RVEC_EXTERN_INTEGER_TEMPLATE(unsigned int)
RVEC_EXTERN_INTEGER_TEMPLATE(unsigned long)
//RVEC_EXTERN_INTEGER_TEMPLATE(unsigned long long)

RVEC_EXTERN_FLOAT_TEMPLATE(float)
RVEC_EXTERN_FLOAT_TEMPLATE(double)

#undef RVEC_EXTERN_UNARY_OPERATOR
#undef RVEC_EXTERN_BINARY_OPERATOR
#undef RVEC_EXTERN_ASSIGN_OPERATOR
#undef RVEC_EXTERN_LOGICAL_OPERATOR
#undef RVEC_EXTERN_INTEGER_TEMPLATE
#undef RVEC_EXTERN_FLOAT_TEMPLATE

#define RVEC_EXTERN_UNARY_FUNCTION(T, NAME, FUNC) \
   extern template RVec<PromoteType<T>> NAME(const RVec<T> &);

#define RVEC_EXTERN_STD_UNARY_FUNCTION(T, F) RVEC_EXTERN_UNARY_FUNCTION(T, F, std::F)

#define RVEC_EXTERN_BINARY_FUNCTION(T0, T1, NAME, FUNC)                            \
   extern template RVec<PromoteTypes<T0, T1>> NAME(const RVec<T0> &, const T1 &); \
   extern template RVec<PromoteTypes<T0, T1>> NAME(const T0 &, const RVec<T1> &); \
   extern template RVec<PromoteTypes<T0, T1>> NAME(const RVec<T0> &, const RVec<T1> &);

#define RVEC_EXTERN_STD_BINARY_FUNCTION(T, F) RVEC_EXTERN_BINARY_FUNCTION(T, T, F, std::F)

#define RVEC_EXTERN_STD_FUNCTIONS(T)             \
   RVEC_EXTERN_STD_UNARY_FUNCTION(T, abs)        \
   RVEC_EXTERN_STD_BINARY_FUNCTION(T, fdim)      \
   RVEC_EXTERN_STD_BINARY_FUNCTION(T, fmod)      \
   RVEC_EXTERN_STD_BINARY_FUNCTION(T, remainder) \
   RVEC_EXTERN_STD_UNARY_FUNCTION(T, exp)        \
   RVEC_EXTERN_STD_UNARY_FUNCTION(T, exp2)       \
   RVEC_EXTERN_STD_UNARY_FUNCTION(T, expm1)      \
   RVEC_EXTERN_STD_UNARY_FUNCTION(T, log)        \
   RVEC_EXTERN_STD_UNARY_FUNCTION(T, log10)      \
   RVEC_EXTERN_STD_UNARY_FUNCTION(T, log2)       \
   RVEC_EXTERN_STD_UNARY_FUNCTION(T, log1p)      \
   RVEC_EXTERN_STD_BINARY_FUNCTION(T, pow)       \
   RVEC_EXTERN_STD_UNARY_FUNCTION(T, sqrt)       \
   RVEC_EXTERN_STD_UNARY_FUNCTION(T, cbrt)       \
   RVEC_EXTERN_STD_BINARY_FUNCTION(T, hypot)     \
   RVEC_EXTERN_STD_UNARY_FUNCTION(T, sin)        \
   RVEC_EXTERN_STD_UNARY_FUNCTION(T, cos)        \
   RVEC_EXTERN_STD_UNARY_FUNCTION(T, tan)        \
   RVEC_EXTERN_STD_UNARY_FUNCTION(T, asin)       \
   RVEC_EXTERN_STD_UNARY_FUNCTION(T, acos)       \
   RVEC_EXTERN_STD_UNARY_FUNCTION(T, atan)       \
   RVEC_EXTERN_STD_BINARY_FUNCTION(T, atan2)     \
   RVEC_EXTERN_STD_UNARY_FUNCTION(T, sinh)       \
   RVEC_EXTERN_STD_UNARY_FUNCTION(T, cosh)       \
   RVEC_EXTERN_STD_UNARY_FUNCTION(T, tanh)       \
   RVEC_EXTERN_STD_UNARY_FUNCTION(T, asinh)      \
   RVEC_EXTERN_STD_UNARY_FUNCTION(T, acosh)      \
   RVEC_EXTERN_STD_UNARY_FUNCTION(T, atanh)      \
   RVEC_EXTERN_STD_UNARY_FUNCTION(T, floor)      \
   RVEC_EXTERN_STD_UNARY_FUNCTION(T, ceil)       \
   RVEC_EXTERN_STD_UNARY_FUNCTION(T, trunc)      \
   RVEC_EXTERN_STD_UNARY_FUNCTION(T, round)      \
   RVEC_EXTERN_STD_UNARY_FUNCTION(T, erf)        \
   RVEC_EXTERN_STD_UNARY_FUNCTION(T, erfc)       \
   RVEC_EXTERN_STD_UNARY_FUNCTION(T, lgamma)     \
   RVEC_EXTERN_STD_UNARY_FUNCTION(T, tgamma)     \

RVEC_EXTERN_STD_FUNCTIONS(float)
RVEC_EXTERN_STD_FUNCTIONS(double)
#undef RVEC_EXTERN_STD_UNARY_FUNCTION
#undef RVEC_EXTERN_STD_BINARY_FUNCTION
#undef RVEC_EXTERN_STD_UNARY_FUNCTIONS

#ifdef R__HAS_VDT

#define RVEC_EXTERN_VDT_UNARY_FUNCTION(T, F) RVEC_EXTERN_UNARY_FUNCTION(T, F, vdt::F)

RVEC_EXTERN_VDT_UNARY_FUNCTION(float, fast_expf)
RVEC_EXTERN_VDT_UNARY_FUNCTION(float, fast_logf)
RVEC_EXTERN_VDT_UNARY_FUNCTION(float, fast_sinf)
RVEC_EXTERN_VDT_UNARY_FUNCTION(float, fast_cosf)
RVEC_EXTERN_VDT_UNARY_FUNCTION(float, fast_tanf)
RVEC_EXTERN_VDT_UNARY_FUNCTION(float, fast_asinf)
RVEC_EXTERN_VDT_UNARY_FUNCTION(float, fast_acosf)
RVEC_EXTERN_VDT_UNARY_FUNCTION(float, fast_atanf)

RVEC_EXTERN_VDT_UNARY_FUNCTION(double, fast_exp)
RVEC_EXTERN_VDT_UNARY_FUNCTION(double, fast_log)
RVEC_EXTERN_VDT_UNARY_FUNCTION(double, fast_sin)
RVEC_EXTERN_VDT_UNARY_FUNCTION(double, fast_cos)
RVEC_EXTERN_VDT_UNARY_FUNCTION(double, fast_tan)
RVEC_EXTERN_VDT_UNARY_FUNCTION(double, fast_asin)
RVEC_EXTERN_VDT_UNARY_FUNCTION(double, fast_acos)
RVEC_EXTERN_VDT_UNARY_FUNCTION(double, fast_atan)

#endif // R__HAS_VDT

#endif // _VECOPS_USE_EXTERN_TEMPLATES

/** @} */ // end of Doxygen group vecops

} // End of VecOps NS

// Allow to use RVec as ROOT::RVec
using ROOT::VecOps::RVec;

} // End of ROOT NS

#endif // ROOT_RVEC
