// @(#)root/smatrix:$Id$
// Authors: T. Glebe, L. Moneta    2005

#ifndef ROOT_Math_MatrixFunctions
#define ROOT_Math_MatrixFunctions
// ********************************************************************
//
// source:
//
// type:      source code
//
// created:   20. Mar 2001
//
// author:    Thorsten Glebe
//            HERA-B Collaboration
//            Max-Planck-Institut fuer Kernphysik
//            Saupfercheckweg 1
//            69117 Heidelberg
//            Germany
//            E-mail: T.Glebe@mpi-hd.mpg.de
//
// Description: Functions/Operators special to Matrix
//
// changes:
// 20 Mar 2001 (TG) creation
// 20 Mar 2001 (TG) Added Matrix * Vector multiplication
// 21 Mar 2001 (TG) Added transpose, product
// 11 Apr 2001 (TG) transpose() speed improvment by removing rows(), cols()
//                  through static members of Matrix and Expr class
//
// ********************************************************************


/**
\defgroup MatrixFunctions Matrix Template Functions
\ingroup SMatrixGroup

These function apply to matrices (and also Matrix expression) and can return a
matrix expression of a particular defined type, like in the matrix multiplication or
a vector, like in the matrix-vector product or a scalar like in the Similarity
vector-matrix product.

*/

#include "Math/BinaryOpPolicy.h"
#include "Math/Expression.h"
#include "Math/HelperOps.h"
#include "Math/CholeskyDecomp.h"

namespace ROOT {

  namespace Math {

    template <class T, unsigned int D>
    class SVector;

#ifdef XXX
//==============================================================================
// SMatrix * SVector
//==============================================================================
template <class T, unsigned int D1, unsigned int D2, class R>
SVector<T,D1> operator*(const SMatrix<T,D1,D2,R>& rhs, const SVector<T,D2>& lhs)
{
  SVector<T,D1> tmp;
  for(unsigned int i=0; i<D1; ++i) {
    const unsigned int rpos = i*D2;
    for(unsigned int j=0; j<D2; ++j) {
      tmp[i] += rhs.apply(rpos+j) * lhs.apply(j);
    }
  }
  return tmp;
}
#endif


// matrix-vector product:
// use apply(i) function for matrices. Tested  (11/05/06) with using (i,j) but
// performances are slightly worse (not clear  why)

//==============================================================================
// meta_row_dot
//==============================================================================
template <unsigned int I>
struct meta_row_dot {
  template <class A, class B>
  static inline typename A::value_type f(const A& lhs, const B& rhs,
                                         const unsigned int offset) {
    return lhs.apply(offset+I) * rhs.apply(I) + meta_row_dot<I-1>::f(lhs,rhs,offset);
  }
};


//==============================================================================
// meta_row_dot<0>
//==============================================================================
template <>
struct meta_row_dot<0> {
  template <class A, class B>
  static inline typename A::value_type f(const A& lhs, const B& rhs,
                                         const unsigned int offset) {
    return lhs.apply(offset) * rhs.apply(0);
  }
};

//==============================================================================
// VectorMatrixRowOp
//==============================================================================
template <class Matrix, class Vector, unsigned int D2>
class VectorMatrixRowOp {
public:

  typedef typename Vector::value_type T;

  ///
  VectorMatrixRowOp(const Matrix& lhs, const Vector& rhs) :
    lhs_(lhs), rhs_(rhs) {}

  ///
  ~VectorMatrixRowOp() {}

  /// calc \f$ \sum_{j} a_{ij} * v_j \f$
  inline typename Matrix::value_type apply(unsigned int i) const {
    return meta_row_dot<D2-1>::f(lhs_, rhs_, i*D2);
  }

   // check if passed pointer is in use
   // check only the vector since this is a vector expression
  inline bool IsInUse (const T * p) const {
    return rhs_.IsInUse(p);
  }


protected:
  const Matrix& lhs_;
  const Vector& rhs_;
};


//==============================================================================
// meta_col_dot
//==============================================================================
template <unsigned int I>
struct meta_col_dot {
  template <class Matrix, class Vector>
  static inline typename Matrix::value_type f(const Matrix& lhs, const Vector& rhs,
                                              const unsigned int offset) {
    return lhs.apply(Matrix::kCols*I+offset) * rhs.apply(I) +
           meta_col_dot<I-1>::f(lhs,rhs,offset);
  }
};


//==============================================================================
// meta_col_dot<0>
//==============================================================================
template <>
struct meta_col_dot<0> {
  template <class Matrix, class Vector>
  static inline typename Matrix::value_type f(const Matrix& lhs, const Vector& rhs,
                                              const unsigned int offset) {
    return lhs.apply(offset) * rhs.apply(0);
  }
};

//==============================================================================
// VectorMatrixColOp
//==============================================================================
/**
   Class for Vector-Matrix multiplication

   @ingroup Expression
 */
template <class Vector, class Matrix, unsigned int D1>
class VectorMatrixColOp {
public:

  typedef typename Vector::value_type T;
  ///
  VectorMatrixColOp(const Vector& lhs, const Matrix& rhs) :
    lhs_(lhs), rhs_(rhs) {}

  ///
  ~VectorMatrixColOp() {}

  /// calc \f$ \sum_{j} a_{ij} * v_j \f$
  inline typename Matrix::value_type apply(unsigned int i) const {
    return meta_col_dot<D1-1>::f(rhs_, lhs_, i);
  }

   // check if passed pointer is in use
   // check only the vector since this is a vector expression
   inline bool IsInUse (const T * p) const {
    return lhs_.IsInUse(p);
  }


protected:
  const Vector&    lhs_;
  const Matrix&    rhs_;
};

/**
   Matrix *  Vector multiplication   \f$ a(i) = \sum_{j} M(i,j) * b(j) \f$
   returning a vector expression

   @ingroup MatrixFunctions
 */
//==============================================================================
// operator*: SMatrix * SVector
//==============================================================================
template <class T, unsigned int D1, unsigned int D2, class R>
inline VecExpr<VectorMatrixRowOp<SMatrix<T,D1,D2,R>,SVector<T,D2>, D2>, T, D1>
 operator*(const SMatrix<T,D1,D2,R>& lhs, const SVector<T,D2>& rhs) {
  typedef VectorMatrixRowOp<SMatrix<T,D1,D2,R>,SVector<T,D2>, D2> VMOp;
  return VecExpr<VMOp, T, D1>(VMOp(lhs,rhs));
}

//==============================================================================
// operator*: SMatrix * Expr<A,T,D2>
//==============================================================================
template <class A, class T, unsigned int D1, unsigned int D2, class R>
inline VecExpr<VectorMatrixRowOp<SMatrix<T,D1,D2,R>, VecExpr<A,T,D2>, D2>, T, D1>
 operator*(const SMatrix<T,D1,D2,R>& lhs, const VecExpr<A,T,D2>& rhs) {
  typedef VectorMatrixRowOp<SMatrix<T,D1,D2,R>,VecExpr<A,T,D2>, D2> VMOp;
  return VecExpr<VMOp, T, D1>(VMOp(lhs,rhs));
}

//==============================================================================
// operator*: Expr<A,T,D1,D2> * SVector
//==============================================================================
template <class A, class T, unsigned int D1, unsigned int D2, class R>
inline VecExpr<VectorMatrixRowOp<Expr<A,T,D1,D2,R>, SVector<T,D2>, D2>, T, D1>
 operator*(const Expr<A,T,D1,D2,R>& lhs, const SVector<T,D2>& rhs) {
  typedef VectorMatrixRowOp<Expr<A,T,D1,D2,R>,SVector<T,D2>, D2> VMOp;
  return VecExpr<VMOp, T, D1>(VMOp(lhs,rhs));
}

//==============================================================================
// operator*: Expr<A,T,D1,D2> * VecExpr<B,T,D2>
//==============================================================================
template <class A, class B, class T, unsigned int D1, unsigned int D2, class R>
inline VecExpr<VectorMatrixRowOp<Expr<A,T,D1,D2,R>, VecExpr<B,T,D2>, D2>, T, D1>
 operator*(const Expr<A,T,D1,D2,R>& lhs, const VecExpr<B,T,D2>& rhs) {
  typedef VectorMatrixRowOp<Expr<A,T,D1,D2,R>,VecExpr<B,T,D2>, D2> VMOp;
  return VecExpr<VMOp, T, D1>(VMOp(lhs,rhs));
}

//==============================================================================
// operator*: SVector * SMatrix
//==============================================================================
template <class T, unsigned int D1, unsigned int D2, class R>
inline VecExpr<VectorMatrixColOp<SVector<T,D1>, SMatrix<T,D1,D2,R>, D1>, T, D2>
 operator*(const SVector<T,D1>& lhs, const SMatrix<T,D1,D2,R>& rhs) {
  typedef VectorMatrixColOp<SVector<T,D1>, SMatrix<T,D1,D2,R>, D1> VMOp;
  return VecExpr<VMOp, T, D2>(VMOp(lhs,rhs));
}

//==============================================================================
// operator*: SVector * Expr<A,T,D1,D2>
//==============================================================================
template <class A, class T, unsigned int D1, unsigned int D2, class R>
inline VecExpr<VectorMatrixColOp<SVector<T,D1>, Expr<A,T,D1,D2,R>, D1>, T, D2>
 operator*(const SVector<T,D1>& lhs, const Expr<A,T,D1,D2,R>& rhs) {
  typedef VectorMatrixColOp<SVector<T,D1>, Expr<A,T,D1,D2,R>, D1> VMOp;
  return VecExpr<VMOp, T, D2>(VMOp(lhs,rhs));
}

//==============================================================================
// operator*: VecExpr<A,T,D1> * SMatrix
//==============================================================================
template <class A, class T, unsigned int D1, unsigned int D2, class R>
inline VecExpr<VectorMatrixColOp<VecExpr<A,T,D1>, SMatrix<T,D1,D2,R>, D1>, T, D2>
 operator*(const VecExpr<A,T,D1>& lhs, const SMatrix<T,D1,D2,R>& rhs) {
  typedef VectorMatrixColOp<VecExpr<A,T,D1>, SMatrix<T,D1,D2,R>, D1> VMOp;
  return VecExpr<VMOp, T, D2>(VMOp(lhs,rhs));
}

//==============================================================================
// operator*: VecExpr<A,T,D1> * Expr<B,T,D1,D2>
//==============================================================================
template <class A, class B, class T, unsigned int D1, unsigned int D2, class R>
inline VecExpr<VectorMatrixColOp<VecExpr<A,T,D1>, Expr<B,T,D1,D2,R>, D1>, T, D2>
 operator*(const VecExpr<A,T,D1>& lhs, const Expr<B,T,D1,D2,R>& rhs) {
  typedef VectorMatrixColOp<VecExpr<A,T,D1>, Expr<B,T,D1,D2,R>, D1> VMOp;
  return VecExpr<VMOp, T, D2>(VMOp(lhs,rhs));
}

//==============================================================================
// meta_matrix_dot
//==============================================================================
template <unsigned int I>
struct meta_matrix_dot {

  template <class MatrixA, class MatrixB>
  static inline typename MatrixA::value_type f(const MatrixA& lhs,
                                               const MatrixB& rhs,
                                               const unsigned int offset) {
    return lhs.apply(offset/MatrixB::kCols*MatrixA::kCols + I) *
           rhs.apply(MatrixB::kCols*I + offset%MatrixB::kCols) +
           meta_matrix_dot<I-1>::f(lhs,rhs,offset);
  }

  // multiplication using i and j indeces
  template <class MatrixA, class MatrixB>
  static inline typename MatrixA::value_type g(const MatrixA& lhs,
                                               const MatrixB& rhs,
                                               unsigned int i,
                                               unsigned int j) {
    return lhs(i, I) * rhs(I , j) +
           meta_matrix_dot<I-1>::g(lhs,rhs,i,j);
  }
};


//==============================================================================
// meta_matrix_dot<0>
//==============================================================================
template <>
struct meta_matrix_dot<0> {

  template <class MatrixA, class MatrixB>
  static inline typename MatrixA::value_type f(const MatrixA& lhs,
                                               const MatrixB& rhs,
                                               const unsigned int offset) {
    return lhs.apply(offset/MatrixB::kCols*MatrixA::kCols) *
           rhs.apply(offset%MatrixB::kCols);
  }

  // multiplication using i and j
  template <class MatrixA, class MatrixB>
  static inline typename MatrixA::value_type g(const MatrixA& lhs,
                                               const MatrixB& rhs,
                                               unsigned int i, unsigned int j) {
    return lhs(i,0) * rhs(0,j);
  }

};

//==============================================================================
// MatrixMulOp
//==============================================================================
/**
   Class for Matrix-Matrix multiplication

   @ingroup Expression
 */
template <class MatrixA, class MatrixB, class T, unsigned int D>
class MatrixMulOp {
public:
  ///
  MatrixMulOp(const MatrixA& lhs, const MatrixB& rhs) :
    lhs_(lhs), rhs_(rhs) {}

  ///
  ~MatrixMulOp() {}

  /// calc \f$\sum_{j} a_{ik} * b_{kj}\f$
  inline T apply(unsigned int i) const {
    return meta_matrix_dot<D-1>::f(lhs_, rhs_, i);
  }

  inline T operator() (unsigned int i, unsigned j) const {
    return meta_matrix_dot<D-1>::g(lhs_, rhs_, i, j);
  }

  inline bool IsInUse (const T * p) const {
    return lhs_.IsInUse(p) || rhs_.IsInUse(p);
  }


protected:
  const MatrixA&    lhs_;
  const MatrixB&    rhs_;
};


/**
   Matrix *  Matrix multiplication , \f$ C(i,j) = \sum_{k} A(i,k) * B(k,j)\f$
   returning a matrix expression

   @ingroup MatrixFunctions
 */
//==============================================================================
// operator* (SMatrix * SMatrix, binary)
//==============================================================================
template <  class T, unsigned int D1, unsigned int D, unsigned int D2, class R1, class R2>
inline Expr<MatrixMulOp<SMatrix<T,D1,D,R1>, SMatrix<T,D,D2,R2>,T,D>, T, D1, D2, typename MultPolicy<T,R1,R2>::RepType>
 operator*(const SMatrix<T,D1,D,R1>& lhs, const SMatrix<T,D,D2,R2>& rhs) {
  typedef MatrixMulOp<SMatrix<T,D1,D,R1>, SMatrix<T,D,D2,R2>, T,D> MatMulOp;
  return Expr<MatMulOp,T,D1,D2,
    typename MultPolicy<T,R1,R2>::RepType>(MatMulOp(lhs,rhs));
}

//==============================================================================
// operator* (SMatrix * Expr, binary)
//==============================================================================
template <class A, class T, unsigned int D1, unsigned int D, unsigned int D2, class R1, class R2>
inline Expr<MatrixMulOp<SMatrix<T,D1,D,R1>, Expr<A,T,D,D2,R2>,T,D>, T, D1, D2, typename MultPolicy<T,R1,R2>::RepType>
 operator*(const SMatrix<T,D1,D,R1>& lhs, const Expr<A,T,D,D2,R2>& rhs) {
  typedef MatrixMulOp<SMatrix<T,D1,D,R1>, Expr<A,T,D,D2,R2>,T,D> MatMulOp;
  return Expr<MatMulOp,T,D1,D2,
    typename MultPolicy<T,R1,R2>::RepType>(MatMulOp(lhs,rhs));
}

//==============================================================================
// operator* (Expr * SMatrix, binary)
//==============================================================================
template <class A, class T, unsigned int D1, unsigned int D, unsigned int D2, class R1, class R2>
inline Expr<MatrixMulOp<Expr<A,T,D1,D,R1>, SMatrix<T,D,D2,R2>,T,D>, T, D1, D2, typename MultPolicy<T,R1,R2>::RepType>
 operator*(const Expr<A,T,D1,D,R1>& lhs, const SMatrix<T,D,D2,R2>& rhs) {
  typedef MatrixMulOp<Expr<A,T,D1,D,R1>, SMatrix<T,D,D2,R2>,T,D> MatMulOp;
  return Expr<MatMulOp,T,D1,D2,
    typename MultPolicy<T,R1,R2>::RepType>(MatMulOp(lhs,rhs));
}

//==============================================================================
// operator* (Expr * Expr, binary)
//==============================================================================
template <class A, class B, class T, unsigned int D1, unsigned int D, unsigned int D2, class R1, class R2>
inline Expr<MatrixMulOp<Expr<A,T,D1,D,R1>, Expr<B,T,D,D2,R2>,T,D>, T, D1, D2, typename MultPolicy<T,R1,R2>::RepType>
 operator*(const Expr<A,T,D1,D,R1>& lhs, const Expr<B,T,D,D2,R2>& rhs) {
  typedef MatrixMulOp<Expr<A,T,D1,D,R1>, Expr<B,T,D,D2,R2>, T,D> MatMulOp;
  return Expr<MatMulOp,T,D1,D2,typename MultPolicy<T,R1,R2>::RepType>(MatMulOp(lhs,rhs));
}



#ifdef XXX
//==============================================================================
// MatrixMulOp
//==============================================================================
template <class MatrixA, class MatrixB, unsigned int D>
class MatrixMulOp {
public:
  ///
  MatrixMulOp(const MatrixA& lhs, const MatrixB& rhs) :
    lhs_(lhs), rhs_(rhs) {}

  ///
  ~MatrixMulOp() {}

  /// calc \f$\sum_{j} a_{ik} * b_{kj}\f$
  inline typename MatrixA::value_type apply(unsigned int i) const {
    return meta_matrix_dot<D-1>::f(lhs_, rhs_, i);
  }

protected:
  const MatrixA&    lhs_;
  const MatrixB&    rhs_;
};


//==============================================================================
// operator* (SMatrix * SMatrix, binary)
//==============================================================================
template <  class T, unsigned int D1, unsigned int D, unsigned int D2, class R1, class R2>
inline Expr<MatrixMulOp<SMatrix<T,D1,D,R1>, SMatrix<T,D,D2,R2>, D>, T, D1, D2, typename MultPolicy<T,R1,R2>::RepType>
 operator*(const SMatrix<T,D1,D,R1>& lhs, const SMatrix<T,D,D2,R2>& rhs) {
  typedef MatrixMulOp<SMatrix<T,D1,D,R1>, SMatrix<T,D,D2,R2>, D> MatMulOp;
  return Expr<MatMulOp,T,D1,D2,typename MultPolicy<T,R1,R2>::RepType>(MatMulOp(lhs,rhs));
}

//==============================================================================
// operator* (SMatrix * Expr, binary)
//==============================================================================
template <class A, class T, unsigned int D1, unsigned int D, unsigned int D2, class R1, class R2>
inline Expr<MatrixMulOp<SMatrix<T,D1,D,R1>, Expr<A,T,D,D2,R2>, D>, T, D1, D2, typename MultPolicy<T,R1,R2>::RepType>
 operator*(const SMatrix<T,D1,D,R1>& lhs, const Expr<A,T,D,D2,R2>& rhs) {
  typedef MatrixMulOp<SMatrix<T,D1,D,R1>, Expr<A,T,D,D2,R2>, D> MatMulOp;
  return Expr<MatMulOp,T,D1,D2,typename MultPolicy<T,R1,R2>::RepType>(MatMulOp(lhs,rhs));
}

//==============================================================================
// operator* (Expr * SMatrix, binary)
//==============================================================================
template <class A, class T, unsigned int D1, unsigned int D, unsigned int D2, class R1, class R2>
inline Expr<MatrixMulOp<Expr<A,T,D1,D,R1>, SMatrix<T,D,D2,R2>, D>, T, D1, D2, typename MultPolicy<T,R1,R2>::RepType>
 operator*(const Expr<A,T,D1,D,R1>& lhs, const SMatrix<T,D,D2,R2>& rhs) {
  typedef MatrixMulOp<Expr<A,T,D1,D,R1>, SMatrix<T,D,D2,R2>, D> MatMulOp;
  return Expr<MatMulOp,T,D1,D2,typename MultPolicy<T,R1,R2>::RepType>(MatMulOp(lhs,rhs));
}

//=============================================================================
// operator* (Expr * Expr, binary)
//=============================================================================
template <class A, class B, class T, unsigned int D1, unsigned int D, unsigned int D2, class R1, class R2>
inline Expr<MatrixMulOp<Expr<A,T,D1,D,R1>, Expr<B,T,D,D2,R2>, D>, T, D1, D2, typename MultPolicy<T,R1,R2>::RepType>
 operator*(const Expr<A,T,D1,D,R1>& lhs, const Expr<B,T,D,D2,R2>& rhs) {
  typedef MatrixMulOp<Expr<A,T,D1,D,R1>, Expr<B,T,D,D2,R2>, D> MatMulOp;
  return Expr<MatMulOp,T,D1,D2,typename MultPolicy<T,R1,R2>::RepType>(MatMulOp(lhs,rhs));
}
#endif

//==============================================================================
// TransposeOp
//==============================================================================
/**
   Class for Transpose Operations

   @ingroup Expression
 */
template <class Matrix, class T, unsigned int D1, unsigned int D2=D1>
class TransposeOp {
public:
  ///
  TransposeOp( const Matrix& rhs) :
    rhs_(rhs) {}

  ///
  ~TransposeOp() {}

  ///
  inline T apply(unsigned int i) const {
    return rhs_.apply( (i%D1)*D2 + i/D1);
  }
  inline T operator() (unsigned int i, unsigned j) const {
    return rhs_( j, i);
  }

  inline bool IsInUse (const T * p) const {
    return rhs_.IsInUse(p);
  }

protected:
  const Matrix& rhs_;
};


/**
   Matrix Transpose   B(i,j) = A(j,i)
   returning a matrix expression

   @ingroup MatrixFunctions
 */
//==============================================================================
// transpose
//==============================================================================
template <class T, unsigned int D1, unsigned int D2, class R>
inline Expr<TransposeOp<SMatrix<T,D1,D2,R>,T,D1,D2>, T, D2, D1, typename TranspPolicy<T,D1,D2,R>::RepType>
 Transpose(const SMatrix<T,D1,D2, R>& rhs) {
  typedef TransposeOp<SMatrix<T,D1,D2,R>,T,D1,D2> MatTrOp;

  return Expr<MatTrOp, T, D2, D1, typename TranspPolicy<T,D1,D2,R>::RepType>(MatTrOp(rhs));
}

//==============================================================================
// transpose
//==============================================================================
template <class A, class T, unsigned int D1, unsigned int D2, class R>
inline Expr<TransposeOp<Expr<A,T,D1,D2,R>,T,D1,D2>, T, D2, D1, typename TranspPolicy<T,D1,D2,R>::RepType>
 Transpose(const Expr<A,T,D1,D2,R>& rhs) {
  typedef TransposeOp<Expr<A,T,D1,D2,R>,T,D1,D2> MatTrOp;

  return Expr<MatTrOp, T, D2, D1, typename TranspPolicy<T,D1,D2,R>::RepType>(MatTrOp(rhs));
}


#ifdef ENABLE_TEMPORARIES_TRANSPOSE
// sometimes is faster to create a temp, not clear why

//==============================================================================
// transpose
//==============================================================================
template <class T, unsigned int D1, unsigned int D2, class R>
inline SMatrix< T, D2, D1, typename TranspPolicy<T,D1,D2,R>::RepType>
 Transpose(const SMatrix<T,D1,D2, R>& rhs) {
  typedef TransposeOp<SMatrix<T,D1,D2,R>,T,D1,D2> MatTrOp;

  return SMatrix< T, D2, D1, typename TranspPolicy<T,D1,D2,R>::RepType>
    ( Expr<MatTrOp, T, D2, D1, typename TranspPolicy<T,D1,D2,R>::RepType>(MatTrOp(rhs)) );
}

//==============================================================================
// transpose
//==============================================================================
template <class A, class T, unsigned int D1, unsigned int D2, class R>
inline  SMatrix< T, D2, D1, typename TranspPolicy<T,D1,D2,R>::RepType>
 Transpose(const Expr<A,T,D1,D2,R>& rhs) {
  typedef TransposeOp<Expr<A,T,D1,D2,R>,T,D1,D2> MatTrOp;

  return SMatrix< T, D2, D1, typename TranspPolicy<T,D1,D2,R>::RepType>
    ( Expr<MatTrOp, T, D2, D1, typename TranspPolicy<T,D1,D2,R>::RepType>(MatTrOp(rhs)) );
}

#endif


#ifdef OLD
//==============================================================================
// product: SMatrix/SVector calculate v^T * A * v
//==============================================================================
template <class T, unsigned int D, class R>
inline T Product(const SMatrix<T,D,D,R>& lhs, const SVector<T,D>& rhs) {
  return Dot(rhs, lhs * rhs);
}

//==============================================================================
// product: SVector/SMatrix calculate v^T * A * v
//==============================================================================
template <class T, unsigned int D, class R>
inline T Product(const SVector<T,D>& lhs, const SMatrix<T,D,D,R>& rhs) {
  return Dot(lhs, rhs * lhs);
}

//==============================================================================
// product: SMatrix/Expr calculate v^T * A * v
//==============================================================================
template <class A, class T, unsigned int D, class R>
inline T Product(const SMatrix<T,D,D,R>& lhs, const VecExpr<A,T,D>& rhs) {
  return Dot(rhs, lhs * rhs);
}

//==============================================================================
// product: Expr/SMatrix calculate v^T * A * v
//==============================================================================
template <class A, class T, unsigned int D, class R>
inline T Product(const VecExpr<A,T,D>& lhs, const SMatrix<T,D,D,R>& rhs) {
  return Dot(lhs, rhs * lhs);
}

//==============================================================================
// product: SVector/Expr calculate v^T * A * v
//==============================================================================
template <class A, class T, unsigned int D, class R>
inline T Product(const SVector<T,D>& lhs, const Expr<A,T,D,D,R>& rhs) {
  return Dot(lhs, rhs * lhs);
}

//==============================================================================
// product: Expr/SVector calculate v^T * A * v
//==============================================================================
template <class A, class T, unsigned int D, class R>
inline T Product(const Expr<A,T,D,D,R>& lhs, const SVector<T,D>& rhs) {
  return Dot(rhs, lhs * rhs);
}

//==============================================================================
// product: Expr/Expr calculate v^T * A * v
//==============================================================================
template <class A, class B, class T, unsigned int D, class R>
inline T Product(const Expr<A,T,D,D,R>& lhs, const VecExpr<B,T,D>& rhs) {
  return Dot(rhs, lhs * rhs);
}

//==============================================================================
// product: Expr/Expr calculate v^T * A * v
//==============================================================================
template <class A, class B, class T, unsigned int D, class R>
inline T Product(const VecExpr<A,T,D>& lhs, const Expr<B,T,D,D,R>& rhs) {
  return Dot(lhs, rhs * lhs);
}
#endif

/**
   Similarity Vector - Matrix Product:  v^T * A * v
   returning a scalar value of type T   \f$ s = \sum_{i,j} v(i) * A(i,j) * v(j)\f$

   @ingroup MatrixFunctions
 */

//==============================================================================
// product: SMatrix/SVector calculate v^T * A * v
//==============================================================================
template <class T, unsigned int D, class R>
inline T Similarity(const SMatrix<T,D,D,R>& lhs, const SVector<T,D>& rhs) {
  return Dot(rhs, lhs * rhs);
}

//==============================================================================
// product: SVector/SMatrix calculate v^T * A * v
//==============================================================================
template <class T, unsigned int D, class R>
inline T Similarity(const SVector<T,D>& lhs, const SMatrix<T,D,D,R>& rhs) {
  return Dot(lhs, rhs * lhs);
}

//==============================================================================
// product: SMatrix/Expr calculate v^T * A * v
//==============================================================================
template <class A, class T, unsigned int D, class R>
inline T Similarity(const SMatrix<T,D,D,R>& lhs, const VecExpr<A,T,D>& rhs) {
  return Dot(rhs, lhs * rhs);
}

//==============================================================================
// product: Expr/SMatrix calculate v^T * A * v
//==============================================================================
template <class A, class T, unsigned int D, class R>
inline T Similarity(const VecExpr<A,T,D>& lhs, const SMatrix<T,D,D,R>& rhs) {
  return Dot(lhs, rhs * lhs);
}

//==============================================================================
// product: SVector/Expr calculate v^T * A * v
//==============================================================================
template <class A, class T, unsigned int D, class R>
inline T Similarity(const SVector<T,D>& lhs, const Expr<A,T,D,D,R>& rhs) {
  return Dot(lhs, rhs * lhs);
}

//==============================================================================
// product: Expr/SVector calculate v^T * A * v
//==============================================================================
template <class A, class T, unsigned int D, class R>
inline T Similarity(const Expr<A,T,D,D,R>& lhs, const SVector<T,D>& rhs) {
  return Dot(rhs, lhs * rhs);
}

//==============================================================================
// product: Expr/Expr calculate v^T * A * v
//==============================================================================
template <class A, class B, class T, unsigned int D, class R>
inline T Similarity(const Expr<A,T,D,D,R>& lhs, const VecExpr<B,T,D>& rhs) {
  return Dot(rhs, lhs * rhs);
}

//==============================================================================
// product: Expr/Expr calculate v^T * A * v
//==============================================================================
template <class A, class B, class T, unsigned int D, class R>
inline T Similarity(const VecExpr<A,T,D>& lhs, const Expr<B,T,D,D,R>& rhs) {
  return Dot(lhs, rhs * lhs);
}


/**
   Similarity Matrix Product :  B = U * A * U^T for A symmetric
   returning a symmetric matrix expression:
   \f$ B(i,j) = \sum_{k,l} U(i,k) * A(k,l) * U(j,l) \f$

   @ingroup MatrixFunctions
 */
//==============================================================================
// product: SMatrix/SMatrix calculate M * A * M^T where A is a symmetric matrix
// return matrix will be nrows M x nrows M
//==============================================================================
template <class T, unsigned int D1, unsigned int D2, class R>
inline SMatrix<T,D1,D1,MatRepSym<T,D1> > Similarity(const SMatrix<T,D1,D2,R>& lhs, const SMatrix<T,D2,D2,MatRepSym<T,D2> >& rhs) {
  SMatrix<T,D1,D2, MatRepStd<T,D1,D2> > tmp = lhs * rhs;
  typedef  SMatrix<T,D1,D1,MatRepSym<T,D1> > SMatrixSym;
  SMatrixSym mret;
  AssignSym::Evaluate(mret,  tmp * Transpose(lhs)  );
  return mret;
}

//==============================================================================
// product: SMatrix/SMatrix calculate M * A * M^T where A is a symmetric matrix
// return matrix will be nrowsM x nrows M
// M is a matrix expression
//==============================================================================
template <class A, class T, unsigned int D1, unsigned int D2, class R>
inline SMatrix<T,D1,D1,MatRepSym<T,D1> > Similarity(const Expr<A,T,D1,D2,R>& lhs, const SMatrix<T,D2,D2,MatRepSym<T,D2> >& rhs) {
  SMatrix<T,D1,D2,MatRepStd<T,D1,D2> > tmp = lhs * rhs;
  typedef  SMatrix<T,D1,D1,MatRepSym<T,D1> > SMatrixSym;
  SMatrixSym mret;
  AssignSym::Evaluate(mret,  tmp * Transpose(lhs)  );
  return mret;
}

#ifdef XXX
    // not needed (
//==============================================================================
// product: SMatrix/SMatrix calculate M * A * M where A and M are symmetric matrices
// return matrix will be nrows M x nrows M
//==============================================================================
template <class T, unsigned int D1>
inline SMatrix<T,D1,D1,MatRepSym<T,D1> > Similarity(const SMatrix<T,D1,D1,MatRepSym<T,D1> >& lhs, const SMatrix<T,D1,D1,MatRepSym<T,D1> >& rhs) {
  SMatrix<T,D1,D1, MatRepStd<T,D1,D1> > tmp = lhs * rhs;
  typedef  SMatrix<T,D1,D1,MatRepSym<T,D1> > SMatrixSym;
  SMatrixSym mret;
  AssignSym::Evaluate(mret,  tmp * lhs  );
  return mret;
}
#endif


/**
   Transpose Similarity Matrix Product :  B = U^T * A * U for A symmetric
   returning a symmetric matrix expression: \f$ B(i,j) = \sum_{k,l} U(k,i) * A(k,l) * U(l,j) \f$

   @ingroup MatrixFunctions
 */
//==============================================================================
// product: SMatrix/SMatrix calculate M^T * A * M where A is a symmetric matrix
// return matrix will be ncolsM x ncols M
//==============================================================================
template <class T, unsigned int D1, unsigned int D2, class R>
inline SMatrix<T,D2,D2,MatRepSym<T,D2> > SimilarityT(const SMatrix<T,D1,D2,R>& lhs, const SMatrix<T,D1,D1,MatRepSym<T,D1> >& rhs) {
  SMatrix<T,D1,D2,MatRepStd<T,D1,D2> > tmp = rhs * lhs;
  typedef  SMatrix<T,D2,D2,MatRepSym<T,D2> > SMatrixSym;
  SMatrixSym mret;
  AssignSym::Evaluate(mret,  Transpose(lhs) * tmp );
  return mret;
}

//==============================================================================
// product: SMatrix/SMatrix calculate M^T * A * M where A is a symmetric matrix
// return matrix will be ncolsM x ncols M
// M is a matrix expression
//==============================================================================
template <class A, class T, unsigned int D1, unsigned int D2, class R>
inline SMatrix<T,D2,D2,MatRepSym<T,D2> > SimilarityT(const Expr<A,T,D1,D2,R>& lhs, const SMatrix<T,D1,D1,MatRepSym<T,D1> >& rhs) {
  SMatrix<T,D1,D2,MatRepStd<T,D1,D2> > tmp = rhs * lhs;
  typedef  SMatrix<T,D2,D2,MatRepSym<T,D2> > SMatrixSym;
  SMatrixSym mret;
  AssignSym::Evaluate(mret,  Transpose(lhs) * tmp );
  return mret;
}





// //==============================================================================
// // Mult * (Expr * Expr, binary) with a symmetric result
// // the operation is done only for half
// //==============================================================================
// template <class A, class B, class T, unsigned int D1, unsigned int D, unsigned int D2, class R1, class R2>
// inline Expr<MatrixMulOp<Expr<A,T,D,D,MatRepSym<T,D> >, Expr<B,T,D,D2,R2>,T,D>, T, D1, D2, typename MultPolicy<T,R1,R2>::RepType>
//  operator*(const Expr<A,T,D1,D,R1>& lhs, const Expr<B,T,D,D2,R2>& rhs) {
//   typedef MatrixMulOp<Expr<A,T,D1,D,R1>, Expr<B,T,D,D2,R2>, T,D> MatMulOp;
//   return Expr<MatMulOp,T,D1,D2,typename MultPolicy<T,R1,R2>::RepType>(MatMulOp(lhs,rhs));
// }



//==============================================================================
// TensorMulOp
//==============================================================================
/**
   Class for Tensor Multiplication (outer product) of two vectors
   giving a matrix

   @ingroup Expression
 */
template <class Vector1, class Vector2>
class TensorMulOp {
public:
  ///
  TensorMulOp( const Vector1 & lhs, const Vector2 & rhs) :
    lhs_(lhs),
    rhs_(rhs) {}

  ///
  ~TensorMulOp() {}

  /// Vector2::kSize is the number of columns in the resulting matrix
  inline typename Vector1::value_type apply(unsigned int i) const {
    return lhs_.apply( i/ Vector2::kSize) * rhs_.apply( i % Vector2::kSize );
  }
  inline typename Vector1::value_type operator() (unsigned int i, unsigned j) const {
    return lhs_.apply(i) * rhs_.apply(j);
  }

  inline bool IsInUse (const typename Vector1::value_type * ) const {
    return false;
  }


protected:

  const Vector1 & lhs_;
  const Vector2 & rhs_;

};



/**
   Tensor Vector Product : M(i,j) = v(i) * v(j)
   returning a matrix expression

   @ingroup VectFunction
 */

#ifndef _WIN32

    // Tensor Prod (use default MatRepStd for the returned expression
    // cannot make a symmetric matrix
//==============================================================================
// TensorProd (SVector x SVector)
//==============================================================================
template <class T, unsigned int D1, unsigned int D2>
inline Expr<TensorMulOp<SVector<T,D1>, SVector<T,D2>  >, T, D1, D2 >
 TensorProd(const SVector<T,D1>& lhs, const SVector<T,D2>& rhs) {
  typedef TensorMulOp<SVector<T,D1>, SVector<T,D2> > TVMulOp;
  return Expr<TVMulOp,T,D1,D2>(TVMulOp(lhs,rhs));
}

//==============================================================================
// TensorProd (VecExpr x SVector)
//==============================================================================
 template <class T, unsigned int D1, unsigned int D2, class A>
inline Expr<TensorMulOp<VecExpr<A,T,D1>, SVector<T,D2>  >, T, D1, D2 >
 TensorProd(const VecExpr<A,T,D1>& lhs, const SVector<T,D2>& rhs) {
  typedef TensorMulOp<VecExpr<A,T,D1>, SVector<T,D2> > TVMulOp;
  return Expr<TVMulOp,T,D1,D2>(TVMulOp(lhs,rhs));
}

//==============================================================================
// TensorProd (SVector x VecExpr)
//==============================================================================
 template <class T, unsigned int D1, unsigned int D2, class A>
inline Expr<TensorMulOp<SVector<T,D1>, VecExpr<A,T,D2>  >, T, D1, D2 >
 TensorProd(const SVector<T,D1>& lhs, const VecExpr<A,T,D2>& rhs) {
  typedef TensorMulOp<SVector<T,D1>, VecExpr<A,T,D2> > TVMulOp;
  return Expr<TVMulOp,T,D1,D2>(TVMulOp(lhs,rhs));
}


//==============================================================================
// TensorProd (VecExpr x VecExpr)
//==============================================================================
 template <class T, unsigned int D1, unsigned int D2, class A, class B>
inline Expr<TensorMulOp<VecExpr<A,T,D1>, VecExpr<B,T,D2>  >, T, D1, D2 >
 TensorProd(const VecExpr<A,T,D1>& lhs, const VecExpr<B,T,D2>& rhs) {
  typedef TensorMulOp<VecExpr<A,T,D1>, VecExpr<B,T,D2> > TVMulOp;
  return Expr<TVMulOp,T,D1,D2>(TVMulOp(lhs,rhs));
}

#endif
#ifdef _WIN32
/// case of WINDOWS - problem using Expression (  C1001: INTERNAL COMPILER ERROR )

//==============================================================================
// TensorProd (SVector x SVector)
//==============================================================================
template <class T, unsigned int D1, unsigned int D2>
inline SMatrix<T,D1,D2,MatRepStd<T, D1, D2>> TensorProd(const SVector<T,D1>& lhs, const SVector<T,D2>& rhs) {
  SMatrix<T,D1,D2,MatRepStd<T, D1, D2>> tmp;
  for (unsigned int i=0; i< D1; ++i)
    for (unsigned int j=0; j< D2; ++j) {
      tmp(i,j) = lhs[i]*rhs[j];
    }

  return tmp;
}
//==============================================================================
// TensorProd (VecExpr x SVector)
//==============================================================================
template <class T, unsigned int D1, unsigned int D2, class A>
inline SMatrix<T,D1,D2,MatRepStd<T, D1, D2>>  TensorProd(const VecExpr<A,T,D1>& lhs, const SVector<T,D2>& rhs) {
  SMatrix<T,D1,D2,MatRepStd<T, D1, D2>> tmp;
  for (unsigned int i=0; i< D1; ++i)
    for (unsigned int j=0; j< D2; ++j)
      tmp(i,j) = lhs.apply(i) * rhs.apply(j);

  return tmp;
}
//==============================================================================
// TensorProd (SVector x VecExpr)
//==============================================================================
template <class T, unsigned int D1, unsigned int D2, class A>
inline SMatrix<T,D1,D2, MatRepStd<T, D1, D2>> TensorProd(const SVector<T,D1>& lhs, const VecExpr<A,T,D2>& rhs) {
  SMatrix<T,D1,D2,MatRepStd<T, D1, D2>> tmp;
  for (unsigned int i=0; i< D1; ++i)
    for (unsigned int j=0; j< D2; ++j)
      tmp(i,j) = lhs.apply(i) * rhs.apply(j);

  return tmp;
}

//==============================================================================
// TensorProd (VecExpr x VecExpr)
//==============================================================================

template <class T, unsigned int D1, unsigned int D2, class A, class B>
inline SMatrix<T,D1,D2,MatRepStd<T, D1, D2>> TensorProd(const VecExpr<A,T,D1>& lhs, const VecExpr<B,T,D2>& rhs) {
  SMatrix<T,D1,D2,MatRepStd<T, D1, D2>> tmp;
  for (unsigned int i=0; i< D1; ++i)
    for (unsigned int j=0; j< D2; ++j)
      tmp(i,j) = lhs.apply(i) * rhs.apply(j);

  return tmp;
}


#endif

// solving a positive defined symmetric linear system using Choleski decompositions
// matrix will be decomposed and the returned vector will be overwritten in vec
// If the user wants to pass const objects need to copy the matrices
// It will work only for symmetric matrices
template <class T, unsigned int D>
bool SolveChol( SMatrix<T, D, D, MatRepSym<T, D>  > & mat,  SVector<T, D> & vec ) {
   CholeskyDecomp<T, D> decomp(mat);
   return decomp.Solve(vec);
}

/// same function as before but not overwriting the matrix and returning a copy of the vector
/// (this is the slow version)
template <class T, unsigned int D>
SVector<T,D> SolveChol( const SMatrix<T, D, D, MatRepSym<T, D>  > & mat,  const SVector<T, D> & vec, int & ifail  ) {
   SMatrix<T, D, D, MatRepSym<T, D> > atmp(mat);
   SVector<T,D> vret(vec);
   bool ok = SolveChol( atmp, vret);
   ifail =  (ok) ? 0 : -1;
   return vret;
}



  }  // namespace Math

}  // namespace ROOT


#endif  /* ROOT_Math_MatrixFunctions */
