// @(#)root/smatrix:$Id$
// Authors: T. Glebe, L. Moneta    2005

#ifndef ROOT_Math_Functions
#define ROOT_Math_Functions
// ********************************************************************
//
// source:
//
// type:      source code
//
// created:   16. Mar 2001
//
// author:    Thorsten Glebe
//            HERA-B Collaboration
//            Max-Planck-Institut fuer Kernphysik
//            Saupfercheckweg 1
//            69117 Heidelberg
//            Germany
//            E-mail: T.Glebe@mpi-hd.mpg.de
//
// Description: Functions which are not applied like a unary operator
//
// changes:
// 16 Mar 2001 (TG) creation
// 03 Apr 2001 (TG) minimum added, doc++ comments added
// 07 Apr 2001 (TG) Lmag2, Lmag added
// 19 Apr 2001 (TG) added #include <cmath>
// 24 Apr 2001 (TG) added sign()
// 26 Jul 2001 (TG) round() added
// 27 Sep 2001 (TG) added Expr declaration
//
// ********************************************************************
#include <cmath>

#include "Math/Expression.h"

/**
\defgroup TempFunction Generic Template Functions
\ingroup SMatrixGroup

   These functions apply for any type T, such as a scalar, a vector or a matrix.
*/

/**
\defgroup VectFunction Vector Template Functions
\ingroup SMatrixGroup

These functions apply to SVector types (and also to Vector expressions) and can
return a vector expression or
a scalar, like in the Dot product, or a matrix, like in the Tensor product
*/


namespace ROOT {

  namespace Math {



template <class T, unsigned int D> class SVector;


/** square
    Template function to compute \f$x\cdot x \f$, for any type T returning a type T

    @ingroup TempFunction
    @author T. Glebe
*/
//==============================================================================
// square: x*x
//==============================================================================
template <class T>
inline const T Square(const T& x) { return x*x; }

/** maximum.
    Template to find max(a,b) where a,b are of type T

    @ingroup TempFunction
    @author T. Glebe
*/
//==============================================================================
// maximum
//==============================================================================
template <class T>
inline const T Maximum(const T& lhs, const T& rhs) {
  return (lhs > rhs) ? lhs : rhs;
}

/** minimum.
    Template to find min(a,b) where a,b are of type T

    @ingroup TempFunction
    @author T. Glebe
*/
//==============================================================================
// minimum
//==============================================================================
template <class T>
inline const T Minimum(const T& lhs, const T& rhs) {
  return (lhs < rhs) ? lhs : rhs;
}

/** round.
    Template to compute nearest integer value for any type T
    @ingroup TempFunction
    @author T. Glebe
*/
//==============================================================================
// round
//==============================================================================
template <class T>
inline int Round(const T& x) {
  return (x-static_cast<int>(x) < 0.5) ? static_cast<int>(x) : static_cast<int>(x+1);
}


/** sign.
    Template to compute the sign of a number

    @ingroup TempFunction
    @author T. Glebe
*/
//==============================================================================
// sign
//==============================================================================
template <class T>
inline int Sign(const T& x) { return (x==0)? 0 : (x<0)? -1 : 1; }

//==============================================================================
// meta_dot
//==============================================================================
template <unsigned int I>
struct meta_dot {
  template <class A, class B, class T>
  static inline T f(const A& lhs, const B& rhs, const T& x) {
    return lhs.apply(I) * rhs.apply(I) + meta_dot<I-1>::f(lhs,rhs,x);
  }
};


//==============================================================================
// meta_dot<0>
//==============================================================================
template <>
struct meta_dot<0> {
  template <class A, class B, class T>
  static inline T f(const A& lhs, const B& rhs, const T& /*x */) {
    return lhs.apply(0) * rhs.apply(0);
  }
};


/**
    Vector dot product.
    Template to compute \f$\vec{a}\cdot\vec{b} = \sum_i a_i\cdot b_i \f$.

    @ingroup VectFunction
    @author T. Glebe
*/
//==============================================================================
// dot
//==============================================================================
template <class T, unsigned int D>
inline T Dot(const SVector<T,D>& lhs, const SVector<T,D>& rhs) {
  return meta_dot<D-1>::f(lhs,rhs, T());
}

//==============================================================================
// dot
//==============================================================================
template <class A, class T, unsigned int D>
inline T Dot(const SVector<T,D>& lhs, const VecExpr<A,T,D>& rhs) {
  return meta_dot<D-1>::f(lhs,rhs, T());
}

//==============================================================================
// dot
//==============================================================================
template <class A, class T, unsigned int D>
inline T Dot(const VecExpr<A,T,D>& lhs, const SVector<T,D>& rhs) {
  return meta_dot<D-1>::f(lhs,rhs, T());
}


//==============================================================================
// dot
//==============================================================================
template <class A, class B, class T, unsigned int D>
inline T Dot(const VecExpr<A,T,D>& lhs, const VecExpr<B,T,D>& rhs) {
  return meta_dot<D-1>::f(rhs,lhs, T());
}


//==============================================================================
// meta_mag
//==============================================================================
template <unsigned int I>
struct meta_mag {
  template <class A, class T>
  static inline T f(const A& rhs, const T& x) {
    return Square(rhs.apply(I)) + meta_mag<I-1>::f(rhs, x);
  }
};


//==============================================================================
// meta_mag<0>
//==============================================================================
template <>
struct meta_mag<0> {
  template <class A, class T>
  static inline T f(const A& rhs, const T& ) {
    return Square(rhs.apply(0));
  }
};


/**
    Vector magnitude square
    Template to compute \f$|\vec{v}|^2 = \sum_iv_i^2 \f$.

    @ingroup VectFunction
    @author T. Glebe
*/
//==============================================================================
// mag2
//==============================================================================
template <class T, unsigned int D>
inline T Mag2(const SVector<T,D>& rhs) {
  return meta_mag<D-1>::f(rhs, T());
}

//==============================================================================
// mag2
//==============================================================================
template <class A, class T, unsigned int D>
inline T Mag2(const VecExpr<A,T,D>& rhs) {
  return meta_mag<D-1>::f(rhs, T());
}

/**
    Vector magnitude (Euclidian norm)
    Compute : \f$ |\vec{v}| = \sqrt{\sum_iv_i^2} \f$.

    @ingroup VectFunction
    @author T. Glebe
*/
//==============================================================================
// mag
//==============================================================================
template <class T, unsigned int D>
inline T Mag(const SVector<T,D>& rhs) {
  return std::sqrt(Mag2(rhs));
}

//==============================================================================
// mag
//==============================================================================
template <class A, class T, unsigned int D>
inline T Mag(const VecExpr<A,T,D>& rhs) {
  return std::sqrt(Mag2(rhs));
}


/** Lmag2: Square of Minkowski Lorentz-Vector norm (only for 4D Vectors)
    Template to compute \f$ |\vec{v}|^2 = v_0^2 - v_1^2 - v_2^2 -v_3^2 \f$.

    @ingroup VectFunction
    @author T. Glebe
*/
//==============================================================================
// Lmag2
//==============================================================================
template <class T>
inline T Lmag2(const SVector<T,4>& rhs) {
  return Square(rhs[0]) - Square(rhs[1]) - Square(rhs[2]) - Square(rhs[3]);
}

//==============================================================================
// Lmag2
//==============================================================================
template <class A, class T>
inline T Lmag2(const VecExpr<A,T,4>& rhs) {
  return Square(rhs.apply(0))
    - Square(rhs.apply(1)) - Square(rhs.apply(2)) - Square(rhs.apply(3));
}

/** Lmag: Minkowski Lorentz-Vector norm (only for 4-dim vectors)
    Length of a vector Lorentz-Vector:
    \f$ |\vec{v}| = \sqrt{v_0^2 - v_1^2 - v_2^2 -v_3^2} \f$.

    @ingroup VectFunction
    @author T. Glebe
*/
//==============================================================================
// Lmag
//==============================================================================
template <class T>
inline T Lmag(const SVector<T,4>& rhs) {
  return std::sqrt(Lmag2(rhs));
}

//==============================================================================
// Lmag
//==============================================================================
template <class A, class T>
inline T Lmag(const VecExpr<A,T,4>& rhs) {
  return std::sqrt(Lmag2(rhs));
}


/** Vector Cross Product (only for 3-dim vectors)
    \f$ \vec{c} = \vec{a}\times\vec{b} \f$.

    @ingroup VectFunction
    @author T. Glebe
*/
//==============================================================================
// cross product
//==============================================================================
template <class T>
inline SVector<T,3> Cross(const SVector<T,3>& lhs, const SVector<T,3>& rhs) {
  return SVector<T,3>(lhs.apply(1)*rhs.apply(2) -
                      lhs.apply(2)*rhs.apply(1),
                      lhs.apply(2)*rhs.apply(0) -
                      lhs.apply(0)*rhs.apply(2),
                      lhs.apply(0)*rhs.apply(1) -
                      lhs.apply(1)*rhs.apply(0));
}

//==============================================================================
// cross product
//==============================================================================
template <class A, class T>
inline SVector<T,3> Cross(const VecExpr<A,T,3>& lhs, const SVector<T,3>& rhs) {
  return SVector<T,3>(lhs.apply(1)*rhs.apply(2) -
                      lhs.apply(2)*rhs.apply(1),
                      lhs.apply(2)*rhs.apply(0) -
                      lhs.apply(0)*rhs.apply(2),
                      lhs.apply(0)*rhs.apply(1) -
                      lhs.apply(1)*rhs.apply(0));
}

//==============================================================================
// cross product
//==============================================================================
template <class T, class A>
inline SVector<T,3> Cross(const SVector<T,3>& lhs, const VecExpr<A,T,3>& rhs) {
  return SVector<T,3>(lhs.apply(1)*rhs.apply(2) -
                      lhs.apply(2)*rhs.apply(1),
                      lhs.apply(2)*rhs.apply(0) -
                      lhs.apply(0)*rhs.apply(2),
                      lhs.apply(0)*rhs.apply(1) -
                      lhs.apply(1)*rhs.apply(0));
}

//==============================================================================
// cross product
//==============================================================================
template <class A, class B, class T>
inline SVector<T,3> Cross(const VecExpr<A,T,3>& lhs, const VecExpr<B,T,3>& rhs) {
  return SVector<T,3>(lhs.apply(1)*rhs.apply(2) -
                      lhs.apply(2)*rhs.apply(1),
                      lhs.apply(2)*rhs.apply(0) -
                      lhs.apply(0)*rhs.apply(2),
                      lhs.apply(0)*rhs.apply(1) -
                      lhs.apply(1)*rhs.apply(0));
}


/** Unit.
    Return a vector of unit length: \f$ \vec{e}_v = \vec{v}/|\vec{v}| \f$.

    @ingroup VectFunction
    @author T. Glebe
*/
//==============================================================================
// unit: returns a unit vector
//==============================================================================
template <class T, unsigned int D>
inline SVector<T,D> Unit(const SVector<T,D>& rhs) {
  return SVector<T,D>(rhs).Unit();
}

//==============================================================================
// unit: returns a unit vector
//==============================================================================
template <class A, class T, unsigned int D>
inline SVector<T,D> Unit(const VecExpr<A,T,D>& rhs) {
  return SVector<T,D>(rhs).Unit();
}

#ifdef XXX
//==============================================================================
// unit: returns a unit vector (worse performance)
//==============================================================================
template <class T, unsigned int D>
inline VecExpr<BinaryOp<DivOp<T>, SVector<T,D>, Constant<T>, T>, T, D>
 unit(const SVector<T,D>& lhs) {
  typedef BinaryOp<DivOp<T>, SVector<T,D>, Constant<T>, T> DivOpBinOp;
  return VecExpr<DivOpBinOp,T,D>(DivOpBinOp(DivOp<T>(),lhs,Constant<T>(mag(lhs))));
}

//==============================================================================
// unit: returns a unit vector (worse performance)
//==============================================================================
template <class A, class T, unsigned int D>
inline VecExpr<BinaryOp<DivOp<T>, VecExpr<A,T,D>, Constant<T>, T>, T, D>
 unit(const VecExpr<A,T,D>& lhs) {
  typedef BinaryOp<DivOp<T>, VecExpr<A,T,D>, Constant<T>, T> DivOpBinOp;
  return VecExpr<DivOpBinOp,T,D>(DivOpBinOp(DivOp<T>(),lhs,Constant<T>(mag(lhs))));
}
#endif


  }  // namespace Math

}  // namespace ROOT



#endif   /* ROOT_Math_Functions */
