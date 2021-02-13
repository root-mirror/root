#ifndef TMVA_RTENSOR
#define TMVA_RTENSOR

#include <vector>
#include <cstddef>     // std::size_t
#include <stdexcept>   // std::runtime_error
#include <sstream>     // std::stringstream
#include <memory>      // std::shared_ptr
#include <type_traits> // std::is_convertible
#include <algorithm>   // std::reverse
#include <iterator>    // std::random_access_iterator_tag

namespace TMVA {
namespace Experimental {

/// Memory layout type
enum class MemoryLayout : uint8_t {
   RowMajor = 0x01,
   ColumnMajor = 0x02
};

namespace Internal {

/// \brief Get size of tensor from shape vector
/// \param[in] shape Shape vector
/// \return Size of contiguous memory
template <typename T>
inline std::size_t GetSizeFromShape(const T &shape)
{
   if (shape.size() == 0)
      return 0;
   std::size_t size = 1;
   for (auto &s : shape)
      size *= s;
   return size;
}

/// \brief Compute strides from shape vector.
/// \param[in] shape Shape vector
/// \param[in] layout Memory layout
/// \return Size of contiguous memory
///
/// This information is needed for the multi-dimensional indexing. See here:
/// https://en.wikipedia.org/wiki/Row-_and_column-major_order
/// https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.strides.html
template <typename T>
inline std::vector<std::size_t> ComputeStridesFromShape(const T &shape, MemoryLayout layout)
{
   const auto size = shape.size();
   T strides(size);
   if (layout == MemoryLayout::RowMajor) {
      for (std::size_t i = 0; i < size; i++) {
         if (i == 0) {
            strides[size - 1 - i] = 1;
         } else {
            strides[size - 1 - i] = strides[size - 1 - i + 1] * shape[size - 1 - i + 1];
         }
      }
   } else if (layout == MemoryLayout::ColumnMajor) {
      for (std::size_t i = 0; i < size; i++) {
         if (i == 0) {
            strides[i] = 1;
         } else {
            strides[i] = strides[i - 1] * shape[i - 1];
         }
      }
   } else {
      std::stringstream ss;
      ss << "Memory layout type is not valid for calculating strides.";
      throw std::runtime_error(ss.str());
   }
   return strides;
}

/// \brief Compute indices from global index
/// \param[in] shape Shape vector
/// \param[in] idx Global index
/// \param[in] layout Memory layout
/// \return Indice vector
template <typename T>
inline T ComputeIndicesFromGlobalIndex(const T& shape, MemoryLayout layout, const typename T::value_type idx)
{
    const auto size = shape.size();
    auto strides = ComputeStridesFromShape(shape, layout);
    T indices(size);
    auto r = idx;
    for (std::size_t i = 0; i < size; i++) {
        indices[i] = int(r / strides[i]);
        r = r % strides[i];
    }
    return indices;
}

/// \brief Compute global index from indices
/// \param[in] strides Strides vector
/// \param[in] idx Indice vector
/// \return Global index
template <typename U, typename V>
inline std::size_t ComputeGlobalIndex(const U& strides, const V& idx)
{
   std::size_t globalIndex = 0;
   const auto size = idx.size();
   for (std::size_t i = 0; i < size; i++) {
      globalIndex += strides[size - 1 - i] * idx[size - 1 - i];
   }
   return globalIndex;
}

/// \brief Type checking for all types of a parameter pack, e.g., used in combination with std::is_convertible
template <class... Ts>
struct and_types : std::true_type {
};

template <class T0, class... Ts>
struct and_types<T0, Ts...> : std::integral_constant<bool, T0() && and_types<Ts...>()> {
};

/// \brief Copy slice of a tensor recursively from here to there
/// \param[in] here Source tensor
/// \param[in] there Target tensor (slice of source tensor)
/// \param[in] mins Minimum of indices for each dimension
/// \param[in] maxs Maximum of indices for each dimension
/// \param[in] idx Current indices
/// \param[in] active Active index needed to stop the recursion
///
/// Copy the content of a slice of a tensor from source to target. This is done
/// by recursively iterating over the ranges of the slice for each dimension.
template <typename T>
void RecursiveCopy(const T &here, T &there,
                   const std::vector<std::size_t> &mins, const std::vector<std::size_t> &maxs,
                   std::vector<std::size_t> idx, std::size_t active)
{
   const auto size = idx.size();
   for (std::size_t i = mins[active]; i < maxs[active]; i++) {
      idx[active] = i;
      if (active == size - 1) {
         auto idxThere = idx;
         for (std::size_t j = 0; j < size; j++) {
            idxThere[j] -= mins[j];
         }
         there(idxThere) = here(idx);
      } else {
         Internal::RecursiveCopy(here, there, mins, maxs, idx, active + 1);
      }
   }
}

} // namespace TMVA::Experimental::Internal

/// \class TMVA::Experimental::RTensor
/// \brief RTensor is a container with contiguous memory and shape information.
/// \tparam T Data-type of the tensor
///
/// An RTensor is a vector-like container, which has additional shape information.
/// The elements of the multi-dimensional container can be accessed by their
/// indices in a coherent way without taking care about the one-dimensional memory
/// layout of the contiguous storage. This also allows to manipulate the shape
/// of the container without moving the actual elements in memory. Another feature
/// is that an RTensor can own the underlying contiguous memory but can also represent
/// only a view on existing data without owning it.
template <typename V, typename C = std::vector<V>>
class RTensor {
public:
   // Typedefs
   using Value_t = V;
   using Shape_t = std::vector<std::size_t>;
   using Index_t = Shape_t;
   using Slice_t = std::vector<Shape_t>;
   using Container_t = C;

private:
   Shape_t fShape;
   Shape_t fStrides;
   std::size_t fSize;
   MemoryLayout fLayout;
   Value_t *fData;
   std::shared_ptr<Container_t> fContainer;

protected:
   void ReshapeInplace(const Shape_t &shape);

public:
   // Constructors

   /// \brief Construct a tensor as view on data
   /// \param[in] data Pointer to data contiguous in memory
   /// \param[in] shape Shape vector
   /// \param[in] layout Memory layout
   RTensor(Value_t *data, Shape_t shape, MemoryLayout layout = MemoryLayout::RowMajor)
      : fShape(shape), fLayout(layout), fData(data), fContainer(NULL)
   {
      fSize = Internal::GetSizeFromShape(shape);
      fStrides = Internal::ComputeStridesFromShape(shape, layout);
   }

   /// \brief Construct a tensor as view on data
   /// \param[in] data Pointer to data contiguous in memory
   /// \param[in] shape Shape vector
   /// \param[in] strides Strides vector
   /// \param[in] layout Memory layout
   RTensor(Value_t *data, Shape_t shape, Shape_t strides, MemoryLayout layout = MemoryLayout::RowMajor)
      : fShape(shape), fStrides(strides), fLayout(layout), fData(data), fContainer(NULL)
   {
      fSize = Internal::GetSizeFromShape(shape);
   }

   /// \brief Construct a tensor owning externally provided data
   /// \param[in] container Shared pointer to data container
   /// \param[in] shape Shape vector
   /// \param[in] layout Memory layout
   RTensor(std::shared_ptr<Container_t> container, Shape_t shape,
           MemoryLayout layout = MemoryLayout::RowMajor)
      : fShape(shape), fLayout(layout), fContainer(container)
   {
      fSize = Internal::GetSizeFromShape(shape);
      fStrides = Internal::ComputeStridesFromShape(shape, layout);
      fData = &(*container->begin());
   }

   /// \brief Construct a tensor owning data initialized with new container
   /// \param[in] shape Shape vector
   /// \param[in] layout Memory layout
   RTensor(Shape_t shape, MemoryLayout layout = MemoryLayout::RowMajor)
      : fShape(shape), fLayout(layout)
   {
      // TODO: Document how data pointer is determined using STL iterator interface.
      // TODO: Sanitize given container type with type traits
      fSize = Internal::GetSizeFromShape(shape);
      fStrides = Internal::ComputeStridesFromShape(shape, layout);
      fContainer = std::make_shared<Container_t>(fSize);
      fData = &(*fContainer->begin());
   }

   // Access elements
   Value_t &operator()(const Index_t &idx);
   const Value_t &operator() (const Index_t &idx) const;
   template <typename... Idx> Value_t &operator()(Idx... idx);
   template <typename... Idx> const Value_t &operator() (Idx... idx) const;

   // Access properties
   std::size_t GetSize() const { return fSize; }
   const Shape_t &GetShape() const { return fShape; }
   const Shape_t &GetStrides() const { return fStrides; }
   Value_t *GetData() { return fData; }
   const Value_t *GetData() const { return fData; }
   std::shared_ptr<Container_t> GetContainer() { return fContainer; }
   const std::shared_ptr<Container_t> GetContainer() const { return fContainer; }
   MemoryLayout GetMemoryLayout() const { return fLayout; }
   bool IsView() const { return fContainer == NULL; }
   bool IsOwner() const { return !IsView(); }

   // Copy
   RTensor<Value_t, Container_t> Copy(MemoryLayout layout = MemoryLayout::RowMajor) const;

   // Transformations
   RTensor<Value_t, Container_t> Transpose() const;
   RTensor<Value_t, Container_t> Squeeze() const;
   RTensor<Value_t, Container_t> ExpandDims(int idx) const;
   RTensor<Value_t, Container_t> Reshape(const Shape_t &shape) const;
   RTensor<Value_t, Container_t> Slice(const Slice_t &slice); 

   // Iterator class
   class Iterator : public std::iterator<std::random_access_iterator_tag, Value_t> {
   private:
      RTensor<Value_t, Container_t>& fTensor;
      Index_t::value_type fGlobalIndex;
   public:
      using difference_type = typename std::iterator<std::random_access_iterator_tag, Value_t>::difference_type;

      Iterator(RTensor<Value_t, Container_t>& x, typename Index_t::value_type idx) : fTensor(x), fGlobalIndex(idx) {}
      Iterator& operator++() { fGlobalIndex++; return *this; }
      Iterator operator++(int) { auto tmp = *this; operator++(); return tmp; }
      Iterator& operator--() { fGlobalIndex--; return *this; }
      Iterator operator--(int) { auto tmp = *this; operator--(); return tmp; }
      Iterator operator+(difference_type rhs) const { return Iterator(fTensor, fGlobalIndex + rhs); }
      Iterator operator-(difference_type rhs) const { return Iterator(fTensor, fGlobalIndex - rhs); }
      difference_type operator-(const Iterator& rhs) { return fGlobalIndex - rhs.GetGlobalIndex(); }
      Iterator& operator+=(difference_type rhs) { fGlobalIndex += rhs; return *this; }
      Iterator& operator-=(difference_type rhs) { fGlobalIndex -= rhs; return *this; }
      Value_t& operator*()
      {
         auto idx = Internal::ComputeIndicesFromGlobalIndex(fTensor.GetShape(), fTensor.GetMemoryLayout(), fGlobalIndex);
         return fTensor(idx);
      }
      bool operator==(const Iterator& rhs) const
      {
         if (fGlobalIndex == rhs.GetGlobalIndex()) return true;
         return false;
      }
      bool operator!=(const Iterator& rhs) const { return !operator==(rhs); };
      bool operator>(const Iterator& rhs) const { return fGlobalIndex > rhs.GetGlobalIndex(); }
      bool operator<(const Iterator& rhs) const { return fGlobalIndex < rhs.GetGlobalIndex(); }
      bool operator>=(const Iterator& rhs) const { return fGlobalIndex >= rhs.GetGlobalIndex(); }
      bool operator<=(const Iterator& rhs) const { return fGlobalIndex <= rhs.GetGlobalIndex(); }
      typename Index_t::value_type GetGlobalIndex() const { return fGlobalIndex; };
   };

   // Iterator interface
   // TODO: Document that the iterator always iterates following the physical memory layout.
   Iterator begin() noexcept {
      return Iterator(*this, 0);
   }
   Iterator end() noexcept {
      return Iterator(*this, fSize);
   }
};

/// \brief Reshape tensor in place
/// \param[in] shape Shape vector
/// Reshape tensor without changing the overall size
template <typename Value_t, typename Container_t>
inline void RTensor<Value_t, Container_t>::ReshapeInplace(const Shape_t &shape)
{
   const auto size = Internal::GetSizeFromShape(shape);
   if (size != fSize) {
      std::stringstream ss;
      ss << "Cannot reshape tensor with size " << fSize << " into shape { ";
      for (std::size_t i = 0; i < shape.size(); i++) {
         if (i != shape.size() - 1) {
            ss << shape[i] << ", ";
         } else {
            ss << shape[i] << " }.";
         }
      }
      throw std::runtime_error(ss.str());
   }

   // Compute new strides from shape
   auto strides = Internal::ComputeStridesFromShape(shape, fLayout);
   fShape = shape;
   fStrides = strides;
}


/// \brief Access elements
/// \param[in] idx Index vector
/// \return Reference to element
template <typename Value_t, typename Container_t>
inline Value_t &RTensor<Value_t, Container_t>::operator()(const Index_t &idx)
{
   const auto globalIndex = Internal::ComputeGlobalIndex(fStrides, idx);
   return fData[globalIndex];
}

/// \brief Access elements
/// \param[in] idx Index vector
/// \return Reference to element
template <typename Value_t, typename Container_t>
inline const Value_t &RTensor<Value_t, Container_t>::operator() (const Index_t &idx) const
{
   const auto globalIndex = Internal::ComputeGlobalIndex(fStrides, idx);
   return fData[globalIndex];
}

/// \brief Access elements
/// \param[in] idx Indices
/// \return Reference to element
template <typename Value_t, typename Container_t>
template <typename... Idx>
Value_t &RTensor<Value_t, Container_t>::operator()(Idx... idx)
{
   static_assert(Internal::and_types<std::is_convertible<Idx, std::size_t>...>{},
                 "Indices are not convertible to std::size_t.");
   return operator()({static_cast<std::size_t>(idx)...});
}

/// \brief Access elements
/// \param[in] idx Indices
/// \return Reference to element
template <typename Value_t, typename Container_t>
template <typename... Idx>
const Value_t &RTensor<Value_t, Container_t>::operator() (Idx... idx) const
{
   static_assert(Internal::and_types<std::is_convertible<Idx, std::size_t>...>{},
                 "Indices are not convertible to std::size_t.");
   return operator()({static_cast<std::size_t>(idx)...});
}

/// \brief Transpose
/// \returns New RTensor
/// The tensor is transposed by inverting the associated memory layout from row-
/// major to column-major and vice versa. Therefore, the underlying data is not
/// touched.
template <typename Value_t, typename Container_t>
inline RTensor<Value_t, Container_t> RTensor<Value_t, Container_t>::Transpose() const
{
   MemoryLayout layout;
   // Transpose by inverting memory layout
   if (fLayout == MemoryLayout::RowMajor) {
      layout = MemoryLayout::ColumnMajor;
   } else if (fLayout == MemoryLayout::ColumnMajor) {
      layout = MemoryLayout::RowMajor;
   } else {
      throw std::runtime_error("Memory layout is not known.");
   }

   // Create copy of container
   RTensor<Value_t, Container_t> x(fData, fShape, fStrides, layout);

   // Reverse shape
   std::reverse(x.fShape.begin(), x.fShape.end());

   // Reverse strides
   std::reverse(x.fStrides.begin(), x.fStrides.end());

   return x;
}

/// \brief Squeeze dimensions
/// \returns New RTensor
/// Squeeze removes the dimensions of size one from the shape.
template <typename Value_t, typename Container_t>
inline RTensor<Value_t, Container_t> RTensor<Value_t, Container_t>::Squeeze() const
{
   // Remove dimensions of one and associated strides
   Shape_t shape;
   Shape_t strides;
   for (std::size_t i = 0; i < fShape.size(); i++) {
      if (fShape[i] != 1) {
         shape.emplace_back(fShape[i]);
         strides.emplace_back(fStrides[i]);
      }
   }

   // If all dimensions are 1, we need to keep one.
   // This does not apply if the inital shape is already empty. Then, return
   // the empty shape.
   if (shape.size() == 0 && fShape.size() != 0) {
      shape.emplace_back(1);
      strides.emplace_back(1);
   }

   // Create copy, attach new shape and strides and return
   RTensor<Value_t, Container_t> x(*this);
   x.fShape = shape;
   x.fStrides = strides;
   return x;
}

/// \brief Expand dimensions
/// \param[in] idx Index in shape vector where dimension is added
/// \returns New RTensor
/// Inserts a dimension of one into the shape.
template <typename Value_t, typename Container_t>
inline RTensor<Value_t, Container_t> RTensor<Value_t, Container_t>::ExpandDims(int idx) const
{
   // Compose shape vector with additional dimensions and adjust strides
   const int len = fShape.size();
   auto shape = fShape;
   auto strides = fStrides;
   if (idx < 0) {
      if (len + idx + 1 < 0) {
         throw std::runtime_error("Given negative index is invalid.");
      }
      shape.insert(shape.end() + 1 + idx, 1);
      strides.insert(strides.begin() + 1 + idx, 1);
   } else {
      if (idx > len) {
         throw std::runtime_error("Given index is invalid.");
      }
      shape.insert(shape.begin() + idx, 1);
      strides.insert(strides.begin() + idx, 1);
   }

   // Create copy, attach new shape and strides and return
   RTensor<Value_t, Container_t> x(*this);
   x.fShape = shape;
   x.fStrides = strides;
   return x;
}

/// \brief Reshape tensor
/// \param[in] shape Shape vector
/// \returns New RTensor
/// Reshape tensor without changing the overall size
template <typename Value_t, typename Container_t>
inline RTensor<Value_t, Container_t> RTensor<Value_t, Container_t>::Reshape(const Shape_t &shape) const
{
   // Create copy, replace and return
   RTensor<Value_t, Container_t> x(*this);
   x.ReshapeInplace(shape);
   return x;
}

/// \brief Create a slice of the tensor
/// \param[in] slice Slice vector
/// \returns New RTensor
/// A slice is a subset of the tensor defined by a vector of pairs of indices.
template <typename Value_t, typename Container_t>
inline RTensor<Value_t, Container_t> RTensor<Value_t, Container_t>::Slice(const Slice_t &slice)
{
   // Sanitize size of slice
   const auto sliceSize = slice.size();
   const auto shapeSize = fShape.size();
   if (sliceSize != shapeSize) {
      std::stringstream ss;
      ss << "Size of slice (" << sliceSize << ") is unequal number of dimensions (" << shapeSize << ").";
      throw std::runtime_error(ss.str());
   }

   // Sanitize slice indices
   // TODO: Sanitize slice indices
   /*
   for (std::size_t i = 0; i < sliceSize; i++) {
   }
   */

   // Convert -1 in slice to proper pair of indices
   // TODO

   // Recompute shape and size
   Shape_t shape(sliceSize);
   for (std::size_t i = 0; i < sliceSize; i++) {
      shape[i] = slice[i][1] - slice[i][0];
   }
   auto size = Internal::GetSizeFromShape(shape);

   // Determine first element contributing to the slice and get the data pointer
   Value_t *data;
   Shape_t idx(sliceSize);
   for (std::size_t i = 0; i < sliceSize; i++) {
      idx[i] = slice[i][0];
   }
   data = &operator()(idx);

   // Create copy and modify properties
   RTensor<Value_t, Container_t> x(*this);
   x.fData = data;
   x.fShape = shape;
   x.fSize = size;

   // Squeeze tensor and return
   return x.Squeeze();
}

/// Copy RTensor to new object
/// \param[in] layout Memory layout of the new RTensor
/// \returns New RTensor
/// The operation copies all elements of the current RTensor to a new RTensor
/// with the given layout contiguous in memory. Note that this copies by default
/// to a row major memory layout.
template <typename Value_t, typename Container_t>
inline RTensor<Value_t, Container_t> RTensor<Value_t, Container_t>::Copy(MemoryLayout layout) const
{
   // Create new tensor with zeros owning the memory
   RTensor<Value_t, Container_t> r(fShape, layout);

   // Copy over the elements from this tensor
   const auto mins = Shape_t(fShape.size());
   const auto maxs = fShape;
   auto idx = mins;
   Internal::RecursiveCopy(*this, r, mins, maxs, idx, 0);

   return r;
}

/// \brief Pretty printing
/// \param[in] os Output stream
/// \param[in] x RTensor
/// \return Modified output stream
template <typename T>
std::ostream &operator<<(std::ostream &os, RTensor<T> &x)
{
   const auto shapeSize = x.GetShape().size();
   if (shapeSize == 1) {
      os << "{ ";
      const auto size = x.GetSize();
      for (std::size_t i = 0; i < size; i++) {
         os << x({i});
         if (i != size - 1)
            os << ", ";
      }
      os << " }";
   } else if (shapeSize == 2) {
      os << "{";
      const auto shape = x.GetShape();
      for (std::size_t i = 0; i < shape[0]; i++) {
         os << " { ";
         for (std::size_t j = 0; j < shape[1]; j++) {
            os << x({i, j});
            if (j < shape[1] - 1) {
               os << ", ";
            } else {
               os << " ";
            }
         }
         os << "}";
      }
      os << " }";
   } else {
      os << "{ printing not yet implemented for this rank }";
   }
   return os;
}

} // namespace TMVA::Experimental
} // namespace TMVA

namespace cling {
template <typename T>
std::string printValue(TMVA::Experimental::RTensor<T> *x)
{
   std::stringstream ss;
   ss << *x;
   return ss.str();
}
} // namespace cling

#endif // TMVA_RTENSOR
