// @(#)root/smatrix:$Id$
// Author: L. Moneta, J. Palacios    2006

#ifndef ROOT_Math_MatrixRepresentationsStatic
#define ROOT_Math_MatrixRepresentationsStatic 1

// Include files

/**
\defgroup MatRep SMatrix Storage Representation
\ingroup SMatrixGroup

Classes MatRepStd and MatRepSym for generic and symmetric matrix
data storage and manipulation. Define data storage and access, plus
operators =, +=, -=, ==.

\author Juan Palacios
\date   2006-01-15
 */

#include "Math/StaticCheck.h"

#include <cstddef>
#include <utility>
#include <type_traits>
#include <array>

namespace ROOT {

namespace Math {

/**
\defgroup MatRepStd Standard Matrix representation
\ingroup MatRep

Standard Matrix representation for a general D1 x D2 matrix.
This class is itself a template on the contained type T, the number of rows and the number of columns.
Its data member is an array T[nrows*ncols] containing the matrix data.
The data are stored in the row-major C convention.
For example, for a matrix, M, of size 3x3, the data \f$ \left[a_0,a_1,a_2,.......,a_7,a_8 \right] \f$d
are stored in the following order:

\f[
M = \left( \begin{array}{ccc}
a_0 & a_1 & a_2  \\
a_3 & a_4  & a_5  \\
a_6 & a_7  & a_8   \end{array} \right)
\f]

*/


   template <class T, unsigned int D1, unsigned int D2=D1>
   class MatRepStd {

   public:

      typedef T  value_type;

      inline const T& operator()(unsigned int i, unsigned int j) const {
         return fArray[i*D2+j];
      }
      inline T& operator()(unsigned int i, unsigned int j) {
         return fArray[i*D2+j];
      }
      inline T& operator[](unsigned int i) { return fArray[i]; }

      inline const T& operator[](unsigned int i) const { return fArray[i]; }

      inline T apply(unsigned int i) const { return fArray[i]; }

      inline T* Array() { return fArray; }

      inline const T* Array() const { return fArray; }

      template <class R>
      inline MatRepStd<T, D1, D2>& operator+=(const R& rhs) {
         for(unsigned int i=0; i<kSize; ++i) fArray[i] += rhs[i];
         return *this;
      }

      template <class R>
      inline MatRepStd<T, D1, D2>& operator-=(const R& rhs) {
         for(unsigned int i=0; i<kSize; ++i) fArray[i] -= rhs[i];
         return *this;
      }

      template <class R>
      inline MatRepStd<T, D1, D2>& operator=(const R& rhs) {
         for(unsigned int i=0; i<kSize; ++i) fArray[i] = rhs[i];
         return *this;
      }

      template <class R>
      inline bool operator==(const R& rhs) const {
         bool rc = true;
         for(unsigned int i=0; i<kSize; ++i) {
            rc = rc && (fArray[i] == rhs[i]);
         }
         return rc;
      }

      enum {
         /// return no. of matrix rows
         kRows = D1,
         /// return no. of matrix columns
         kCols = D2,
         /// return no of elements: rows*columns
         kSize = D1*D2
      };

   private:
      //T __attribute__ ((aligned (16))) fArray[kSize];
      T  fArray[kSize];
   };


//     template<unigned int D>
//     struct Creator {
//       static const RowOffsets<D> & Offsets() {
//          static RowOffsets<D> off;
//           return off;
//       }

   /**
      Static structure to keep the conversion from (i,j) to offsets in the storage data for a
      symmetric matrix
   */

   template<unsigned int D>
   struct RowOffsets {
      inline RowOffsets() {
         int v[D];
         v[0]=0;
         for (unsigned int i=1; i<D; ++i)
            v[i]=v[i-1]+i;
         for (unsigned int i=0; i<D; ++i) {
            for (unsigned int j=0; j<=i; ++j)
               fOff[i*D+j] = v[i]+j;
            for (unsigned int j=i+1; j<D; ++j)
               fOff[i*D+j] = v[j]+i ;
         }
      }
      inline int operator()(unsigned int i, unsigned int j) const { return fOff[i*D+j]; }
      inline int apply(unsigned int i) const { return fOff[i]; }
      int fOff[D*D];
   };

  namespace rowOffsetsUtils {

    ///////////
    // Some meta template stuff
    template<int...> struct indices{};

    template<int I, class IndexTuple, int N>
    struct make_indices_impl;

    template<int I, int... Indices, int N>
    struct make_indices_impl<I, indices<Indices...>, N>
    {
      typedef typename make_indices_impl<I + 1, indices<Indices..., I>,
					 N>::type type;
    };

    template<int N, int... Indices>
    struct make_indices_impl<N, indices<Indices...>, N> {
      typedef indices<Indices...> type;
    };

    template<int N>
    struct make_indices : make_indices_impl<0, indices<>, N> {};
    // end of stuff



    template<int I0, class F, int... I>
    constexpr std::array<decltype(std::declval<F>()(std::declval<int>())), sizeof...(I)>
    do_make(F f, indices<I...>)
    {
      return  std::array<decltype(std::declval<F>()(std::declval<int>())),
			 sizeof...(I)>{{ f(I0 + I)... }};
    }

    template<int N, int I0 = 0, class F>
    constexpr std::array<decltype(std::declval<F>()(std::declval<int>())), N>
    make(F f) {
      return do_make<I0>(f, typename make_indices<N>::type());
    }

  } // namespace rowOffsetsUtils


//_________________________________________________________________________________
   /**
      MatRepSym
      Matrix storage representation for a symmetric matrix of dimension NxN
      This class is a template on the contained type and on the symmetric matrix size, N.
      It has as data member an array of type T of size N*(N+1)/2,
      containing the lower diagonal block of the matrix.
      The order follows the lower diagonal block, still in a row-major convention.
      For example for a symmetric 3x3 matrix the order of the 6 elements
      \f$ \left[a_0,a_1.....a_5 \right]\f$ is:
      \f[
      M = \left( \begin{array}{ccc}
      a_0 & a_1  & a_3  \\
      a_1 & a_2  & a_4  \\
      a_3 & a_4 & a_5   \end{array} \right)
      \f]

      @ingroup MatRep
   */
   template <class T, unsigned int D>
   class MatRepSym {

   public:

    /* constexpr */ inline MatRepSym(){}

    typedef T  value_type;


    inline T & operator()(unsigned int i, unsigned int j)
     { return fArray[offset(i, j)]; }

     inline /* constexpr */ T const & operator()(unsigned int i, unsigned int j) const
     { return fArray[offset(i, j)]; }

     inline T& operator[](unsigned int i) {
       return fArray[off(i)];
     }

     inline /* constexpr */ T const & operator[](unsigned int i) const {
       return fArray[off(i)];
     }

     inline /* constexpr */ T apply(unsigned int i) const {
       return fArray[off(i)];
     }

     inline T* Array() { return fArray; }

     inline const T* Array() const { return fArray; }

      /**
         assignment : only symmetric to symmetric allowed
       */
      template <class R>
      inline MatRepSym<T, D>& operator=(const R&) {
         STATIC_CHECK(0==1,
                      Cannot_assign_general_to_symmetric_matrix_representation);
         return *this;
      }
      inline MatRepSym<T, D>& operator=(const MatRepSym& rhs) {
         for(unsigned int i=0; i<kSize; ++i) fArray[i] = rhs.Array()[i];
         return *this;
      }

      /**
         self addition : only symmetric to symmetric allowed
       */
      template <class R>
      inline MatRepSym<T, D>& operator+=(const R&) {
         STATIC_CHECK(0==1,
                      Cannot_add_general_to_symmetric_matrix_representation);
         return *this;
      }
      inline MatRepSym<T, D>& operator+=(const MatRepSym& rhs) {
         for(unsigned int i=0; i<kSize; ++i) fArray[i] += rhs.Array()[i];
         return *this;
      }

      /**
         self subtraction : only symmetric to symmetric allowed
       */
      template <class R>
      inline MatRepSym<T, D>& operator-=(const R&) {
         STATIC_CHECK(0==1,
                      Cannot_substract_general_to_symmetric_matrix_representation);
         return *this;
      }
      inline MatRepSym<T, D>& operator-=(const MatRepSym& rhs) {
         for(unsigned int i=0; i<kSize; ++i) fArray[i] -= rhs.Array()[i];
         return *this;
      }
      template <class R>
      inline bool operator==(const R& rhs) const {
         bool rc = true;
         for(unsigned int i=0; i<D*D; ++i) {
            rc = rc && (operator[](i) == rhs[i]);
         }
         return rc;
      }

      enum {
         /// return no. of matrix rows
         kRows = D,
         /// return no. of matrix columns
         kCols = D,
         /// return no of elements: rows*columns
         kSize = D*(D+1)/2
      };

     static constexpr int off0(int i) { return i==0 ? 0 : off0(i-1)+i;}
     static constexpr int off2(int i, int j) { return j<i ? off0(i)+j : off0(j)+i; }
     static constexpr int off1(int i) { return off2(i/D, i%D);}

     static int off(int i) {
       static constexpr auto v = rowOffsetsUtils::make<D*D>(off1);
       return v[i];
     }

     static inline constexpr unsigned int
     offset(unsigned int i, unsigned int j)
     {
       //if (j > i) std::swap(i, j);
       return off(i*D+j);
       // return (i>j) ? (i * (i+1) / 2) + j :  (j * (j+1) / 2) + i;
     }

   private:
      //T __attribute__ ((aligned (16))) fArray[kSize];
      T fArray[kSize];
   };



} // namespace Math
} // namespace ROOT


#endif // MATH_MATRIXREPRESENTATIONSSTATIC_H
