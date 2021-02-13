// @(#)root/mathcore:$Id$
// Authors: W. Brown, M. Fischler, L. Moneta    2005

 /**********************************************************************
  *                                                                    *
  * Copyright (c) 2005 ROOT FNAL MathLib Team                          *
  *                                                                    *
  *                                                                    *
  **********************************************************************/

// Header file for Boost
//
// Created by: Mark Fischler  Mon Nov 1  2005
//
// Last update: $Id$
//
#ifndef ROOT_Math_GenVector_Boost
#define ROOT_Math_GenVector_Boost 1

#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/PxPyPzE4D.h"
#include "Math/GenVector/DisplacementVector3D.h"
#include "Math/GenVector/Cartesian3D.h"

#include "Math/GenVector/BoostX.h"
#include "Math/GenVector/BoostY.h"
#include "Math/GenVector/BoostZ.h"

namespace ROOT {

  namespace Math {

//__________________________________________________________________________________________
  /**
     Lorentz boost class with the (4D) transformation represented internally
     by a 4x4 orthosymplectic matrix.
     See also BoostX, BoostY and BoostZ for classes representing
     specialized Lorentz boosts.
     Also, the 3-D rotation classes can be considered to be special Lorentz
     transformations which do not mix space and time components.

     @ingroup GenVector

  */

class Boost {

public:

  typedef double Scalar;

  enum ELorentzRotationMatrixIndex {
      kLXX =  0, kLXY =  1, kLXZ =  2, kLXT =  3
    , kLYX =  4, kLYY =  5, kLYZ =  6, kLYT =  7
    , kLZX =  8, kLZY =  9, kLZZ = 10, kLZT = 11
    , kLTX = 12, kLTY = 13, kLTZ = 14, kLTT = 15
  };

  enum EBoostMatrixIndex {
      kXX =  0, kXY =  1, kXZ =  2, kXT =  3
     , kYY =  4, kYZ =  5, kYT =  6
     , kZZ =  7, kZT =  8
     , kTT =  9
  };

  // ========== Constructors and Assignment =====================

  /**
      Default constructor (identity transformation)
  */
  Boost() { SetIdentity(); }

  /**
     Construct given a three Scalars beta_x, beta_y, and beta_z
   */
  Boost(Scalar beta_x, Scalar beta_y, Scalar beta_z)
   { SetComponents(beta_x, beta_y, beta_z); }

  /**
     Construct given a beta vector (which must have methods x(), y(), z())
   */
  template <class Avector>
  explicit
  Boost(const Avector & beta) { SetComponents(beta); }

  /**
     Construct given a pair of pointers or iterators defining the
     beginning and end of an array of three Scalars to use as beta_x, _y, and _z
   */
  template<class IT>
  Boost(IT begin, IT end) { SetComponents(begin,end); }

   /**
      copy constructor
   */
   Boost(Boost const & b) {
      *this = b;
   }

  /**
     Construct from an axial boost
  */

  explicit Boost( BoostX const & bx ) {SetComponents(bx.BetaVector());}
  explicit Boost( BoostY const & by ) {SetComponents(by.BetaVector());}
  explicit Boost( BoostZ const & bz ) {SetComponents(bz.BetaVector());}

  // The compiler-generated copy ctor, copy assignment, and dtor are OK.

   /**
      Assignment operator
    */
   Boost &
   operator=(Boost const & rhs ) {
    for (unsigned int i=0; i < 10; ++i) {
       fM[i] = rhs.fM[i];
    }
    return *this;
   }

  /**
     Assign from an axial pure boost
  */
  Boost &
  operator=( BoostX const & bx ) { return operator=(Boost(bx)); }
  Boost &
  operator=( BoostY const & by ) { return operator=(Boost(by)); }
  Boost &
  operator=( BoostZ const & bz ) { return operator=(Boost(bz)); }

  /**
     Re-adjust components to eliminate small deviations from a perfect
     orthosyplectic matrix.
   */
  void Rectify();

  // ======== Components ==============

  /**
     Set components from beta_x, beta_y, and beta_z
  */
  void
  SetComponents (Scalar beta_x, Scalar beta_y, Scalar beta_z);

  /**
     Get components into beta_x, beta_y, and beta_z
  */
  void
  GetComponents (Scalar& beta_x, Scalar& beta_y, Scalar& beta_z) const;

  /**
     Set components from a beta vector
  */
  template <class Avector>
  void
  SetComponents (const Avector & beta)
   { SetComponents(beta.x(), beta.y(), beta.z()); }

  /**
     Set given a pair of pointers or iterators defining the beginning and end of
     an array of three Scalars to use as beta_x,beta _y, and beta_z
   */
  template<class IT>
  void SetComponents(IT begin, IT end) {
    IT a = begin; IT b = ++begin; IT c = ++begin;
    (void)end;
    assert (++begin==end);
    SetComponents (*a, *b, *c);
  }

  /**
     Get given a pair of pointers or iterators defining the beginning and end of
     an array of three Scalars into which to place beta_x, beta_y, and beta_z
   */
  template<class IT>
  void GetComponents(IT begin, IT end) const {
    IT a = begin; IT b = ++begin; IT c = ++begin;
    (void)end;
    assert (++begin==end);
    GetComponents (*a, *b, *c);
  }

  /**
     Get given a pointer or an iterator defining the beginning of
     an array into which to place beta_x, beta_y, and beta_z
   */
  template<class IT>
  void GetComponents(IT begin ) const {
     double bx,by,bz = 0;
     GetComponents (bx,by,bz);
     *begin++ = bx;
     *begin++ = by;
     *begin = bz;
  }

  /**
     The beta vector for this boost
   */
  typedef  DisplacementVector3D<Cartesian3D<double>, DefaultCoordinateSystemTag > XYZVector;
  XYZVector BetaVector() const;

  /**
     Get elements of internal 4x4 symmetric representation, into a data
     array suitable for direct use as the components of a LorentzRotation
     Note -- 16 Scalars will be written into the array; if the array is not
     that large, then this will lead to undefined behavior.
  */
  void
  GetLorentzRotation (Scalar r[]) const;

  // =========== operations ==============

  /**
     Lorentz transformation operation on a Minkowski ('Cartesian')
     LorentzVector
  */
  LorentzVector< ROOT::Math::PxPyPzE4D<double> >
  operator() (const LorentzVector< ROOT::Math::PxPyPzE4D<double> > & v) const;

  /**
     Lorentz transformation operation on a LorentzVector in any
     coordinate system
   */
  template <class CoordSystem>
  LorentzVector<CoordSystem>
  operator() (const LorentzVector<CoordSystem> & v) const {
    LorentzVector< PxPyPzE4D<double> > xyzt(v);
    LorentzVector< PxPyPzE4D<double> > r_xyzt = operator()(xyzt);
    return LorentzVector<CoordSystem> ( r_xyzt );
  }

  /**
     Lorentz transformation operation on an arbitrary 4-vector v.
     Preconditions:  v must implement methods x(), y(), z(), and t()
     and the arbitrary vector type must have a constructor taking (x,y,z,t)
   */
  template <class Foreign4Vector>
  Foreign4Vector
  operator() (const Foreign4Vector & v) const {
    LorentzVector< PxPyPzE4D<double> > xyzt(v);
    LorentzVector< PxPyPzE4D<double> > r_xyzt = operator()(xyzt);
    return Foreign4Vector ( r_xyzt.X(), r_xyzt.Y(), r_xyzt.Z(), r_xyzt.T() );
  }

  /**
     Overload operator * for boost on a vector
   */
  template <class A4Vector>
  inline
  A4Vector operator* (const A4Vector & v) const
  {
    return operator()(v);
  }

  /**
      Invert a Boost in place
   */
  void Invert();

  /**
      Return inverse of  a boost
   */
  Boost Inverse() const;

  /**
     Equality/inequality operators
   */
  bool operator == (const Boost & rhs) const {
    for (unsigned int i=0; i < 10; ++i) {
      if( fM[i] != rhs.fM[i] )  return false;
    }
    return true;
  }
  bool operator != (const Boost & rhs) const {
    return ! operator==(rhs);
  }

protected:

  void SetIdentity();

private:

  Scalar fM[10];

};  // Boost

// ============ Class Boost ends here ============

/**
   Stream Output and Input
 */
  // TODO - I/O should be put in the manipulator form

std::ostream & operator<< (std::ostream & os, const Boost & b);


} //namespace Math
} //namespace ROOT







#endif /* ROOT_Math_GenVector_Boost  */
