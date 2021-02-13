// @(#)root/smatrix:$Id$
// Authors: T. Glebe, L. Moneta    2005

#ifndef ROOT_Math_Expression
#define ROOT_Math_Expression
// ********************************************************************
//
// source:
//
// type:      source code
//
// created:   19. Mar 2001
//
// author:    Thorsten Glebe
//            HERA-B Collaboration
//            Max-Planck-Institut fuer Kernphysik
//            Saupfercheckweg 1
//            69117 Heidelberg
//            Germany
//            E-mail: T.Glebe@mpi-hd.mpg.de
//
// Description: Expression Template Elements for SVector
//
// changes:
// 19 Mar 2001 (TG) creation
// 20 Mar 2001 (TG) added rows(), cols() to Expr
// 21 Mar 2001 (TG) added Expr::value_type
// 11 Apr 2001 (TG) rows(), cols() replaced by rows, cols
// 10 Okt 2001 (TG) added print() and operator<<() for Expr class
//
// ********************************************************************

/**
\defgroup Expression Expression Template Classes
\ingroup SMatrixGroup
*/

//==============================================================================
// Expr: class representing SVector expressions
//=============================================================================

// modified BinaryOp with two extension BinaryOpCopyL and BinaryOpCopyR to store the
// object in BinaryOp by value and not reference. When used with constant BinaryOp reference give problems
// on some compilers (like Windows) where a temporary Constant object is ccreated and then destructed


#include <iomanip>
#include <iostream>

namespace ROOT {

  namespace Math {



//    template <class T, unsigned int D, unsigned int D2> class MatRepStd;

/**
    Expression wrapper class for Vector objects

    @ingroup Expression
*/
template <class ExprType, class T, unsigned int D >
class VecExpr {

public:
  typedef T  value_type;

  ///
  VecExpr(const ExprType& rhs) :
    rhs_(rhs) {}

  ///
  ~VecExpr() {}

   ///
  inline T apply(unsigned int i) const {
    return rhs_.apply(i);
  }

  inline T operator() (unsigned int i) const {
    return rhs_.apply(i);
  }


#ifdef OLD_IMPL
  ///
  static const unsigned int rows = D;
  ///
  ///static const unsigned int cols = D2;
#else
  // use enumerations
  enum {

    kRows = D

  };
#endif

  /**
      function to  determine if any use operand
      is being used (has same memory adress)
   */
  inline bool IsInUse (const T * p) const {
    return rhs_.IsInUse(p);
  }


  /// used by operator<<()
  std::ostream& print(std::ostream& os) const {
    os.setf(std::ios::right,std::ios::adjustfield);
    unsigned int i=0;
    os << "[ ";
    for(; i<D-1; ++i) {
      os << apply(i) << ", ";
    }
    os << apply(i);
    os << " ]";

    return os;
  }

private:
  ExprType rhs_; // cannot be a reference!
};


/**
    Expression wrapper class for Matrix objects

    @ingroup Expression
*/

template <class T, unsigned int D, unsigned int D2> class MatRepStd;

template <class ExprType, class T, unsigned int D, unsigned int D2 = 1,
     class R1=MatRepStd<T,D,D2> >
class Expr {
public:
  typedef T  value_type;

  ///
  Expr(const ExprType& rhs) :
    rhs_(rhs) {}

  ///
  ~Expr() {}

  ///
  inline T apply(unsigned int i) const {
    return rhs_.apply(i);
  }
  inline T operator() (unsigned int i, unsigned j) const {
    return rhs_(i,j);
  }

  /**
      function to  determine if any use operand
      is being used (has same memory adress)
   */
  inline bool IsInUse (const T * p) const {
    return rhs_.IsInUse(p);
  }



#ifdef OLD_IMPL
  ///
  static const unsigned int rows = D;
  ///
  static const unsigned int cols = D2;
#else
  // use enumerations
  enum {
    ///
    kRows = D,
  ///
    kCols = D2
  };
#endif

  /// used by operator<<()
  /// simplify to use apply(i,j)
  std::ostream& print(std::ostream& os) const {
    os.setf(std::ios::right,std::ios::adjustfield);
      os << "[ ";
      for (unsigned int i=0; i < D; ++i) {
         unsigned int d2 = D2; // to avoid some annoying warnings in case of vectors (D2 = 0)
         for (unsigned int j=0; j < D2; ++j) {
            os << std::setw(12) << this->operator() (i,j);
            if ((!((j+1)%12)) && (j < d2-1))
            os << std::endl << "         ...";
         }
         if (i != D - 1)
         os << std::endl  << "  ";
      }
     os << " ]";

     return os;
  }

private:
  ExprType rhs_; // cannot be a reference!
};

//==============================================================================
// operator<<
//==============================================================================
template <class A, class T, unsigned int D>
inline std::ostream& operator<<(std::ostream& os, const VecExpr<A,T,D>& rhs) {
  return rhs.print(os);
}

template <class A, class T, unsigned int D1, unsigned int D2, class R1>
inline std::ostream& operator<<(std::ostream& os, const Expr<A,T,D1,D2,R1>& rhs) {
  return rhs.print(os);
}

/**
    BinaryOperation class
    A class representing binary operators in the parse tree.
    This is the default case where objects are kept by reference

    @ingroup  Expression
    @author T. Glebe
*/



//==============================================================================
// BinaryOp
//==============================================================================
template <class Operator, class LHS, class RHS, class T>
class BinaryOp {
public:
  ///
  BinaryOp( Operator /* op */, const LHS& lhs, const RHS& rhs) :
    lhs_(lhs), rhs_(rhs) {}

  ///
  ~BinaryOp() {}

  ///
  inline T apply(unsigned int i) const {
    return Operator::apply(lhs_.apply(i), rhs_.apply(i));
  }
  inline T operator() (unsigned int i, unsigned int j) const {
    return Operator::apply(lhs_(i,j), rhs_(i,j) );
  }

  inline bool IsInUse (const T * p) const {
    return lhs_.IsInUse(p) || rhs_.IsInUse(p);
  }

protected:

  const LHS& lhs_;
  const RHS& rhs_;

};

//LM :: add specialization of BinaryOP when first or second argument needs to be copied
// (maybe it can be doen with a template specialization, but it is not worth, easier to have a separate class

//==============================================================================
/**
   Binary Operation class with value storage for the left argument.
   Special case of BinaryOp where for the left argument the passed object
   is copied and stored by value instead of a reference.
   This is used in the case of operations involving a constant, where we cannot store a
   reference to the constant (we get a temporary object) and we need to copy it.

   @ingroup  Expression
*/
//==============================================================================
template <class Operator, class LHS, class RHS, class T>
class BinaryOpCopyL {
public:
  ///
  BinaryOpCopyL( Operator /* op */, const LHS& lhs, const RHS& rhs) :
    lhs_(lhs), rhs_(rhs) {}

  ///
  ~BinaryOpCopyL() {}

  ///
  inline T apply(unsigned int i) const {
    return Operator::apply(lhs_.apply(i), rhs_.apply(i));
  }
  inline T operator() (unsigned int i, unsigned int j) const {
    return Operator::apply(lhs_(i,j), rhs_(i,j) );
  }

  inline bool IsInUse (const T * p) const {
    // no need to check left since we copy it
    return rhs_.IsInUse(p);
  }

protected:

  const LHS  lhs_;
  const RHS& rhs_;

};


//==============================================================================
/**
   Binary Operation class with value storage for the right argument.
   Special case of BinaryOp where for the wight argument a copy is stored instead of a reference
   This is use in the case for example of constant where we cannot store by reference
   but need to copy since Constant is a temporary object

   @ingroup  Expression
*/
//==============================================================================
template <class Operator, class LHS, class RHS, class T>
class BinaryOpCopyR {
public:
  ///
  BinaryOpCopyR( Operator /* op */, const LHS& lhs, const RHS& rhs) :
    lhs_(lhs), rhs_(rhs) {}

  ///
  ~BinaryOpCopyR() {}

  ///
  inline T apply(unsigned int i) const {
    return Operator::apply(lhs_.apply(i), rhs_.apply(i));
  }
  inline T operator() (unsigned int i, unsigned int j) const {
    return Operator::apply(lhs_(i,j), rhs_(i,j) );
  }

  inline bool IsInUse (const T * p) const {
    // no need for right since we copied
    return lhs_.IsInUse(p);
  }

protected:

  const LHS&  lhs_;
  const RHS rhs_;

};



/**
    UnaryOperation class
    A class representing unary operators in the parse tree.
    The objects are stored by reference

    @ingroup  Expression
    @author T. Glebe
*/
//==============================================================================
// UnaryOp
//==============================================================================
template <class Operator, class RHS, class T>
class UnaryOp {
public:
  ///
  UnaryOp( Operator /* op */ , const RHS& rhs) :
    rhs_(rhs) {}

  ///
  ~UnaryOp() {}

  ///
  inline T apply(unsigned int i) const {
    return Operator::apply(rhs_.apply(i));
  }
  inline T operator() (unsigned int i, unsigned int j) const {
    return Operator::apply(rhs_(i,j));
  }

  inline bool IsInUse (const T * p) const {
    return rhs_.IsInUse(p);
  }

protected:

  const RHS& rhs_;

};


/**
    Constant expression class
    A class representing constant expressions (literals) in the parse tree.

    @ingroup Expression
    @author T. Glebe
*/
//==============================================================================
// Constant
//==============================================================================
template <class T>
class Constant {
public:
  ///
  Constant( const T& rhs ) :
    rhs_(rhs) {}

  ///
  ~Constant() {}

  ///
  inline T apply(unsigned int /*i */ ) const { return rhs_; }

  inline T operator() (unsigned int /*i */, unsigned int /*j */ ) const { return rhs_; }

  //inline bool IsInUse (const T * ) const { return false; }

protected:

  const T rhs_;  // no need for reference. It is  a fundamental type normally


};



  }  // namespace Math

}  // namespace ROOT



#endif  /* ROOT_Math_Expression */
