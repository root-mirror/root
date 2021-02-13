// @(#)root/tmva/tmva/dnn:$Id$
// Authors: Sitong An, Lorenzo Moneta 10/2019

/*************************************************************************
 * Copyright (C) 2019, ROOT                                              *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////
// Definition of the CpuTensor class used to represent  //
// tensor data  in deep neural nets (CNN, RNN, etc..)   //
//////////////////////////////////////////////////////////

#ifndef TMVA_DNN_ARCHITECTURES_CPU_CPUTENSOR
#define TMVA_DNN_ARCHITECTURES_CPU_CPUTENSOR

#include <cstddef>


#include "TMatrix.h"
#include "TMVA/Config.h"
#include "CpuBuffer.h"
#include "CpuMatrix.h"
#include <TMVA/RTensor.hxx>

namespace TMVA {
namespace DNN {

// CPU Tensor Class
// It is a simple wrapper for TMVA RTensor based on
// memory owned by CPU Buffer
// We need to keep a pointer for CPUBuffer for fast conversion
// without copying to TCpuMatrix
// also provides compatibility with old interface

template <typename AFloat>
class TCpuTensor : public TMVA::Experimental::RTensor<AFloat, TCpuBuffer<AFloat>> {

private:
   //TCpuTensor will have no extra private members than RTensor
public:
   friend class TCpuMatrix<AFloat>;

   using Shape_t = typename TMVA::Experimental::RTensor<AFloat>::Shape_t;
   using MemoryLayout = TMVA::Experimental::MemoryLayout;
   using Matrix_t = TCpuMatrix<AFloat>;
   using Scalar_t = AFloat;

   // default constructor
   TCpuTensor(): TMVA::Experimental::RTensor<AFloat, TCpuBuffer<AFloat>>(std::make_shared<TCpuBuffer<AFloat>>(0), {0})
   {}

   /** constructors from n m */
   TCpuTensor(size_t n, size_t m, MemoryLayout memlayout = MemoryLayout::ColumnMajor)
      : TMVA::Experimental::RTensor<AFloat, TCpuBuffer<AFloat>>(std::make_shared<TCpuBuffer<AFloat>>(n * m), {n, m}, memlayout)
   {}

   /** constructors from batch size, depth, height*width */
   TCpuTensor(size_t bsize, size_t depth, size_t hw, MemoryLayout memlayout = MemoryLayout::ColumnMajor)
      : TMVA::Experimental::RTensor<AFloat, TCpuBuffer<AFloat>>(std::make_shared<TCpuBuffer<AFloat>>(bsize * depth * hw), {depth, hw, bsize}, memlayout)
   {
      if (memlayout == MemoryLayout::RowMajor)
         this->ReshapeInplace({bsize, depth, hw});
   }

   /** constructors from batch size, depth, height, width */
   TCpuTensor(size_t bsize, size_t depth, size_t height, size_t width,
              MemoryLayout memlayout = MemoryLayout::ColumnMajor)
      : TMVA::Experimental::RTensor<AFloat, TCpuBuffer<AFloat>>(std::make_shared<TCpuBuffer<AFloat>>(bsize * depth * height * width),
      {depth, height, width, bsize}, memlayout)
   {
      if (memlayout == MemoryLayout::RowMajor)
         this->ReshapeInplace({bsize, depth, height, width});
   }

   /** constructors from a shape.*/
   TCpuTensor(Shape_t shape, MemoryLayout memlayout = MemoryLayout::ColumnMajor)
      : TMVA::Experimental::RTensor<AFloat, TCpuBuffer<AFloat>>(std::make_shared<TCpuBuffer<AFloat>>(TMVA::Experimental::Internal::GetSizeFromShape(shape)),
      shape, memlayout)
   {}

    /* constructors from a AFloat pointer  and a shape. This is a copy */

   TCpuTensor(AFloat *data, const Shape_t &shape,
              MemoryLayout memlayout = MemoryLayout::ColumnMajor)
      : TMVA::Experimental::RTensor<AFloat, TCpuBuffer<AFloat>>(std::make_shared<TCpuBuffer<AFloat>>(TMVA::Experimental::Internal::GetSizeFromShape(shape)), shape, memlayout)
   {
      auto& container = *(this->GetContainer());
      for (size_t i = 0; i <  this->GetSize(); ++i) container[i] = data[i];
   }



   /** constructors from a TCpuBuffer and a shape */
   //unsafe method for backwards compatibility, const not promised. A view.
   TCpuTensor(const TCpuBuffer<AFloat>& buffer, Shape_t shape, MemoryLayout memlayout = MemoryLayout::ColumnMajor)
      : TMVA::Experimental::RTensor<AFloat, TCpuBuffer<AFloat>>(std::make_shared<TCpuBuffer<AFloat>>(buffer), shape, memlayout) {
         R__ASSERT(this->GetSize() <= this->GetContainer()->GetSize());
      }



   /** constructors from a TCpuMatrix. Memory layout is forced to be same as matrix (i.e. columnlayout) */
   //unsafe method for backwards compatibility, const not promised. A view of underlying data.
   TCpuTensor(const TCpuMatrix<AFloat> &matrix, size_t dim = 3, MemoryLayout memlayout = MemoryLayout::ColumnMajor)
      : TMVA::Experimental::RTensor<AFloat, TCpuBuffer<AFloat>>(std::make_shared<TCpuBuffer<AFloat>>(matrix.GetBuffer()),{matrix.GetNrows(), matrix.GetNcols()}, memlayout)
   {

      if (dim >  2) {
         Shape_t shape = this->GetShape();

         if (this->GetLayout() == MemoryLayout::ColumnMajor) {
            shape.insert(shape.end(),dim-2, 1);
         } else {
            shape.insert(shape.begin(), dim - 2, 1);
         }
         this->ReshapeInplace(shape);
      }
   }


   /** Convert to a TMatrixT<AFloat_t> object. Performs a deep copy of the matrix
    *  elements. */

   operator TMatrixT<AFloat>() const {
      // this should work only for size 2 or 4 tensors
      if (this->GetShape().size() == 2 || (this->GetShape().size() == 3 && GetFirstSize() == 1)) {
         TCpuMatrix<AFloat> temp = GetMatrix();
         return temp;
      }
      // convert as a flat vector
      return TMatrixT<AFloat>(1, this->GetSize(), this->GetData());
   }


   /** Return raw pointer to the elements stored contiguously in column-major
    *  order. */
   AFloat *GetRawDataPointer() { return *(this->GetContainer()); }
   const AFloat *GetRawDataPointer() const { return *(this->GetContainer()); }

   // for same API as CudaTensor (device buffer is the CpuBuffer)
   const TCpuBuffer<AFloat> & GetDeviceBuffer()     const {return *(this->GetContainer());}
   TCpuBuffer<AFloat>       & GetDeviceBuffer()           {return *(this->GetContainer());}


   size_t GetNoElements() const { return this->GetSize(); }

   // return the size of the first dimension (if in row order) or last dimension if in column order
   // Tensor is  F x H x W x...for row order layout FHWC
   // or      H x W x ... x F  for column order layout CHWF
   // logic copied from TCudaTensor
   size_t GetFirstSize() const
   {
      auto& shape = this->GetShape();
      return (this->GetMemoryLayout() == MemoryLayout::ColumnMajor) ? shape.back() : shape.front();
   }

   size_t GetCSize() const
   {
      auto& shape = this->GetShape();
      if (shape.size() == 2)  return 1;
      return (this->GetMemoryLayout() == MemoryLayout::ColumnMajor) ? shape.front() : shape[1]; // assume NHWC
   }
   //
   size_t GetHSize() const
   {
      auto& shape = this->GetShape();
      if (shape.size() == 2)  return shape[0];
      if (shape.size() == 3)  return (this->GetMemoryLayout() == MemoryLayout::ColumnMajor) ? shape[0] : shape[1] ;
      if (shape.size() >= 4)  return shape[2] ;
      return 0;

   }
   size_t GetWSize() const
   {
      auto& shape = this->GetShape();
      if (shape.size() == 2)  return shape[1];
      if (shape.size() == 3)  return (this->GetMemoryLayout() == MemoryLayout::ColumnMajor) ? shape[1] : shape[2] ;
      if (shape.size() >= 4)  return shape[3] ;
      return 0;

   }

   // for backward compatibility (assume column-major 
   // for backward compatibility : for CM tensor (n1,n2,n3,n4) -> ( n1*n2*n3, n4)
   //                              for RM tensor (n1,n2,n3,n4) -> ( n2*n3*n4, n1 ) ???
   size_t GetNrows() const { return (GetLayout() == MemoryLayout::ColumnMajor ) ? this->GetStrides().back() : this->GetShape().front();}
   size_t GetNcols() const { return (GetLayout() == MemoryLayout::ColumnMajor ) ? this->GetShape().back() : this->GetStrides().front(); }


   MemoryLayout GetLayout() const { return this->GetMemoryLayout(); }

   //this will be an unsafe view. Method exists for backwards compatibility only
   TCpuMatrix<AFloat> GetMatrix() const
   {
      size_t ndims = 0;
      auto& shape = this->GetShape();
      //check if squeezable but do not actually squeeze
      for (auto& shape_i : shape){
         if (shape_i != 1) {
            ndims++;
         }
      }
      assert(ndims <= 2 && shape.size() > 1);  // to support shape cases {n,1}
      return TCpuMatrix<AFloat>(*(this->GetContainer()), GetHSize(), GetWSize());
   }

   // Create copy, replace and return
   TCpuTensor<AFloat> Reshape(Shape_t shape) const
   {
      TCpuTensor<AFloat> x(*this);
      x.ReshapeInplace(shape);
      return x;
   }

      // return a view of slices in the first dimension (if row wise) or last dimension if colun wise
      // so single event slices
      TCpuTensor<AFloat> At(size_t i)
      {
         auto &shape = this->GetShape();
         auto layout = this->GetMemoryLayout();
         Shape_t sliced_shape = (layout == MemoryLayout::RowMajor) ? Shape_t(shape.begin() + 1, shape.end())
                                                                   : Shape_t(shape.begin(), shape.end() - 1);

         size_t buffsize = (layout == MemoryLayout::RowMajor) ? this->GetStrides().front() : this->GetStrides().back();
         size_t offset = i * buffsize;

         return TCpuTensor<AFloat>(this->GetContainer()->GetSubBuffer(offset, buffsize), sliced_shape, layout);
      }

      TCpuTensor<AFloat> At(size_t i) const { return (const_cast<TCpuTensor<AFloat> &>(*this)).At(i); }

      // for compatibility with old tensor (std::vector<matrix>)
      TCpuMatrix<AFloat> operator[](size_t i) const {
         assert(this->GetMemoryLayout() == MemoryLayout::ColumnMajor );  
         return At(i).GetMatrix(); 
      }

      // set all the tensor contents to zero
      void Zero()
      {
         AFloat *data = *(this->GetContainer());
         for (size_t i = 0; i < this->GetSize(); ++i)
            data[i] = 0;
      }

      // access single element - assume tensor dim is 2
      AFloat &operator()(size_t i, size_t j)
      {
         auto &shape = this->GetShape();
         assert(shape.size() == 2);
         return (this->GetMemoryLayout() == MemoryLayout::RowMajor) ? (*(this->GetContainer()))[i * shape[1] + j]
                                                                    : (*(this->GetContainer()))[j * shape[0] + i];
      }

      // access single element - assume tensor dim is 3. First index i is always the major  indipendent of row-major or
      // column major row- major  I - J - K    . Column- major  is  J - K - I
      AFloat &operator()(size_t i, size_t j, size_t k)
      {
         auto &shape = this->GetShape();
         assert(shape.size() == 3);

         return (this->GetMemoryLayout() == MemoryLayout::RowMajor)
                   ? (*(this->GetContainer()))[i * shape[1] * shape[2] + j * shape[2] + k]
                   : (*(this->GetContainer()))[i * shape[0] * shape[1] + k * shape[0] + j]; // note that is J-K-I
      }

      // access single element - assume tensor dim is 2
      AFloat operator()(size_t i, size_t j) const
      {
         auto &shape = this->GetShape();
         assert(shape.size() == 2);
         return (this->GetMemoryLayout() == MemoryLayout::RowMajor) ? (this->GetData())[i * shape[1] + j]
                                                                    : (this->GetData())[j * shape[0] + i];
      }

      AFloat operator()(size_t i, size_t j, size_t k) const
      {
         auto &shape = this->GetShape();
         assert(shape.size() == 3);

         return (this->GetMemoryLayout() == MemoryLayout::RowMajor)
                   ? (this->GetData())[i * shape[1] * shape[2] + j * shape[2] + k]
                   : (this->GetData())[i * shape[0] * shape[1] + k * shape[0] + j]; // note that is J-K-I
      }

      /** Map the given function over the matrix elements. Executed in parallel
       *  using TThreadExecutor. */
      template <typename Function_t>
      void Map(Function_t & f);

      /** Same as maps but takes the input values from the tensor \p A and writes
       *  the results in this tensor. */
      template <typename Function_t>
      void MapFrom(Function_t & f, const TCpuTensor<AFloat> &A);

      size_t GetBufferUseCount() const { return this->GetContainer()->GetUseCount(); }

      void Print(const char *name = "Tensor") const
      {
         PrintShape(name);

         for (size_t i = 0; i < this->GetSize(); i++)
            std::cout << (this->GetData())[i] << "  ";
         std::cout << std::endl;
      }
      void PrintShape(const char *name = "Tensor") const
      {
         std::string memlayout = (GetLayout() == MemoryLayout::RowMajor) ? "RowMajor" : "ColMajor";
         std::cout << name << " shape : { ";
         auto &shape = this->GetShape();
         for (size_t i = 0; i < shape.size() - 1; ++i)
            std::cout << shape[i] << " , ";
         std::cout << shape.back() << " } "
                   << " Layout : " << memlayout << std::endl;
      }
};

//______________________________________________________________________________
template <typename AFloat>
template <typename Function_t>
inline void TCpuTensor<AFloat>::Map(Function_t &f)
{
   AFloat *data = GetRawDataPointer();
   size_t nelements = GetNoElements();
   size_t nsteps = TCpuMatrix<AFloat>::GetNWorkItems(nelements);

   auto ff = [data, &nsteps, &nelements, &f](UInt_t workerID) {
      size_t jMax = std::min(workerID + nsteps, nelements);
      for (size_t j = workerID; j < jMax; ++j) {
         data[j] = f(data[j]);
      }
      return 0;
   };

   if (nsteps < nelements) {
      TMVA::Config::Instance().GetThreadExecutor().Foreach(ff, ROOT::TSeqI(0, nelements, nsteps));

      // for (size_t i = 0;  i < nelements; i+=nsteps)
      //    ff(i);

   } else {
      R__ASSERT(nelements == nsteps);
      ff(0);
   }
}

//______________________________________________________________________________
template <typename AFloat>
template <typename Function_t>
inline void TCpuTensor<AFloat>::MapFrom(Function_t &f, const TCpuTensor<AFloat> &A)
{
   AFloat *dataB = GetRawDataPointer();
   const AFloat *dataA = A.GetRawDataPointer();

   size_t nelements = GetNoElements();
   R__ASSERT(nelements == A.GetNoElements());
   size_t nsteps = TCpuMatrix<AFloat>::GetNWorkItems(nelements);

   auto ff = [&dataB, &dataA, &nsteps, &nelements, &f](UInt_t workerID) {
      size_t jMax = std::min(workerID + nsteps, nelements);
      for (size_t j = workerID; j < jMax; ++j) {
         dataB[j] = f(dataA[j]);
      }
      return 0;
   };
   if (nsteps < nelements) {
      TMVA::Config::Instance().GetThreadExecutor().Foreach(ff, ROOT::TSeqI(0, nelements, nsteps));
      // for (size_t i = 0;  i < nelements; i+=nsteps)
      //    ff(i);

   } else {
      R__ASSERT(nelements == nsteps);
      ff(0);
   }
}


} // namespace DNN
} // namespace TMVA

#endif
