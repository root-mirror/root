// @(#)root/mathcore:$Id$
// Authors: W. Brown, M. Fischler, L. Moneta    2005

 /**********************************************************************
  *                                                                    *
  * Copyright (c) 2005 , LCG ROOT FNAL MathLib Team                    *
  *                                                                    *
  *                                                                    *
  **********************************************************************/

// Header file for class Rotation in 3 dimensions, represented by 3x3 matrix
//
// Created by: Mark Fischler Thurs June 9  2005
//
// Last update: $Id$
//
#ifndef ROOT_Math_GenVector_Rotation3D
#define ROOT_Math_GenVector_Rotation3D  1


#include "Math/GenVector/Cartesian3D.h"
#include "Math/GenVector/DisplacementVector3D.h"
#include "Math/GenVector/PositionVector3D.h"
#include "Math/GenVector/PxPyPzE4D.h"
#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/3DConversions.h"
#include "Math/GenVector/3DDistances.h"

#include "Math/GenVector/Rotation3Dfwd.h"
#include "Math/GenVector/AxisAnglefwd.h"
#include "Math/GenVector/EulerAnglesfwd.h"
#include "Math/GenVector/Quaternionfwd.h"
#include "Math/GenVector/RotationXfwd.h"
#include "Math/GenVector/RotationYfwd.h"
#include "Math/GenVector/RotationZfwd.h"


#include <algorithm>
#include <cassert>
#include <iostream>


namespace ROOT {
namespace Math {


//__________________________________________________________________________________________
  /**
     Rotation class with the (3D) rotation represented by
     a 3x3 orthogonal matrix.
     This is the optimal representation for application to vectors.
     See also ROOT::Math::AxisAngle, ROOT::Math::EulerAngles, and ROOT::Math::Quaternion for
     classes which have conversion operators to Rotation3D.

     All Rotations types (not only Rotation3D) can be applied to all 3D Vector classes
     (like ROOT::Math::DisplacementVector3D and ROOT::Math::PositionVector3D)
     and also to the 4D Vectors (ROOT::Math::LorentzVector classes), acting on the 3D components.
     A rotation operation is applied by using the operator() or the operator *.
     With the operator * is possible also to combine rotations.
     Note that the operator is NOT commutative, the order how the rotations are applied is relevant.

     @ingroup GenVector
  */

class Rotation3D {

public:

   typedef double Scalar;

   enum ERotation3DMatrixIndex {
      kXX = 0, kXY = 1, kXZ = 2
      , kYX = 3, kYY = 4, kYZ = 5
      , kZX = 6, kZY = 7, kZZ = 8
   };

   // ========== Constructors and Assignment =====================

   /**
      Default constructor (identity rotation)
   */
   Rotation3D();

   /**
      Construct given a pair of pointers or iterators defining the
      beginning and end of an array of nine Scalars
   */
   template<class IT>
   Rotation3D(IT begin, IT end) { SetComponents(begin,end); }

   /**
      copy constructor
   */
   Rotation3D ( Rotation3D const   & r ) {
      *this = r;
   }

   /**
      Construct from an AxisAngle
   */
   explicit Rotation3D( AxisAngle const   & a ) { gv_detail::convert(a, *this); }

   /**
      Construct from EulerAngles
   */
   explicit Rotation3D( EulerAngles const & e ) { gv_detail::convert(e, *this); }

   /**
      Construct from RotationZYX
   */
   explicit Rotation3D( RotationZYX const & e ) { gv_detail::convert(e, *this); }

   /**
      Construct from a Quaternion
   */
   explicit Rotation3D( Quaternion const  & q ) { gv_detail::convert(q, *this); }

   /**
      Construct from an axial rotation
   */
   explicit Rotation3D( RotationZ const & r ) { gv_detail::convert(r, *this); }
   explicit Rotation3D( RotationY const & r ) { gv_detail::convert(r, *this); }
   explicit Rotation3D( RotationX const & r ) { gv_detail::convert(r, *this); }

   /**
      Construct from a linear algebra matrix of size at least 3x3,
      which must support operator()(i,j) to obtain elements (0,0) thru (2,2).
      Precondition:  The matrix is assumed to be orthonormal.  No checking
      or re-adjusting is performed.
   */
   template<class ForeignMatrix>
   explicit Rotation3D(const ForeignMatrix & m) { SetComponents(m); }

   /**
      Construct from three orthonormal vectors (which must have methods
      x(), y() and z()) which will be used as the columns of the rotation
      matrix.  The orthonormality will be checked, and values adjusted
      so that the result will always be a good rotation matrix.
   */
   template<class ForeignVector>
   Rotation3D(const ForeignVector& v1,
              const ForeignVector& v2,
              const ForeignVector& v3 ) { SetComponents(v1, v2, v3); }

   // compiler generated destruuctor is ok

   /**
      Raw constructor from nine Scalar components (without any checking)
   */
   Rotation3D(Scalar  xx, Scalar  xy, Scalar  xz,
              Scalar  yx, Scalar  yy, Scalar  yz,
              Scalar  zx, Scalar  zy, Scalar  zz)
   {
      SetComponents (xx, xy, xz, yx, yy, yz, zx, zy, zz);
   }

   // need to implement assignment operator to avoid using the templated one

   /**
      Assignment operator
   */
   Rotation3D &
   operator=( Rotation3D const   & rhs ) {
      SetComponents( rhs.fM[0], rhs.fM[1], rhs.fM[2],
                     rhs.fM[3], rhs.fM[4], rhs.fM[5],
                     rhs.fM[6], rhs.fM[7], rhs.fM[8] );
      return *this;
   }

   /**
      Assign from an AxisAngle
   */
   Rotation3D &
   operator=( AxisAngle const   & a ) { return operator=(Rotation3D(a)); }

   /**
      Assign from EulerAngles
   */
   Rotation3D &
   operator=( EulerAngles const & e ) { return operator=(Rotation3D(e)); }

   /**
      Assign from RotationZYX
   */
   Rotation3D &
   operator=( RotationZYX const & r ) { return operator=(Rotation3D(r)); }

   /**
      Assign from a Quaternion
   */
   Rotation3D &
   operator=( Quaternion const  & q ) {return operator=(Rotation3D(q)); }

   /**
      Assign from an axial rotation
   */
   Rotation3D &
   operator=( RotationZ const & r ) { return operator=(Rotation3D(r)); }
   Rotation3D &
   operator=( RotationY const & r ) { return operator=(Rotation3D(r)); }
   Rotation3D &
   operator=( RotationX const & r ) { return operator=(Rotation3D(r)); }

   /**
      Assign from an orthonormal linear algebra matrix of size 3x3,
      which must support operator()(i,j) to obtain elements (0,0) thru (2,2).
   */
   template<class ForeignMatrix>
   Rotation3D &
   operator=(const ForeignMatrix & m) {
      SetComponents( m(0,0), m(0,1), m(0,2),
                     m(1,0), m(1,1), m(1,2),
                     m(2,0), m(2,1), m(2,2) );
      return *this;
   }

   /**
      Re-adjust components to eliminate small deviations from perfect
      orthonormality.
   */
   void Rectify();

   // ======== Components ==============

   /**
      Set components from three orthonormal vectors (which must have methods
      x(), y() and z()) which will be used as the columns of the rotation
      matrix.  The orthonormality will be checked, and values adjusted
      so that the result will always be a good rotation matrix.
   */
   template<class ForeignVector>
   void
   SetComponents (const ForeignVector& v1,
                  const ForeignVector& v2,
                  const ForeignVector& v3 ) {
      fM[kXX]=v1.x();  fM[kXY]=v2.x();  fM[kXZ]=v3.x();
      fM[kYX]=v1.y();  fM[kYY]=v2.y();  fM[kYZ]=v3.y();
      fM[kZX]=v1.z();  fM[kZY]=v2.z();  fM[kZZ]=v3.z();
      Rectify();
   }

   /**
      Get components into three vectors which will be the (orthonormal)
      columns of the rotation matrix.  (The vector class must have a
      constructor from 3 Scalars.)
   */
   template<class ForeignVector>
   void
   GetComponents ( ForeignVector& v1,
                   ForeignVector& v2,
                   ForeignVector& v3 ) const {
      v1 = ForeignVector ( fM[kXX], fM[kYX], fM[kZX] );
      v2 = ForeignVector ( fM[kXY], fM[kYY], fM[kZY] );
      v3 = ForeignVector ( fM[kXZ], fM[kYZ], fM[kZZ] );
   }

   /**
      Set the 9 matrix components given an iterator to the start of
      the desired data, and another to the end (9 past start).
   */
   template<class IT>
   void SetComponents(IT begin, IT end) {
      for (int i = 0; i <9; ++i) {
         fM[i] = *begin;
         ++begin;
      }
      (void)end;
      assert (end==begin);
   }

   /**
      Get the 9 matrix components into data specified by an iterator begin
      and another to the end of the desired data (9 past start).
   */
   template<class IT>
   void GetComponents(IT begin, IT end) const {
      for (int i = 0; i <9; ++i) {
         *begin = fM[i];
         ++begin;
      }
      (void)end;
      assert (end==begin);
   }

   /**
      Get the 9 matrix components into data specified by an iterator begin
   */
   template<class IT>
   void GetComponents(IT begin) const {
      std::copy ( fM, fM+9, begin );
   }

   /**
      Set components from a linear algebra matrix of size at least 3x3,
      which must support operator()(i,j) to obtain elements (0,0) thru (2,2).
      Precondition:  The matrix is assumed to be orthonormal.  NO checking
      or re-adjusting is performed.
   */
   template<class ForeignMatrix>
   void
   SetRotationMatrix (const ForeignMatrix & m) {
      fM[kXX]=m(0,0);  fM[kXY]=m(0,1);  fM[kXZ]=m(0,2);
      fM[kYX]=m(1,0);  fM[kYY]=m(1,1);  fM[kYZ]=m(1,2);
      fM[kZX]=m(2,0);  fM[kZY]=m(2,1);  fM[kZZ]=m(2,2);
   }

   /**
      Get components into a linear algebra matrix of size at least 3x3,
      which must support operator()(i,j) for write access to elements
      (0,0) thru (2,2).
   */
   template<class ForeignMatrix>
   void
   GetRotationMatrix (ForeignMatrix & m) const {
      m(0,0)=fM[kXX];  m(0,1)=fM[kXY];  m(0,2)=fM[kXZ];
      m(1,0)=fM[kYX];  m(1,1)=fM[kYY];  m(1,2)=fM[kYZ];
      m(2,0)=fM[kZX];  m(2,1)=fM[kZY];  m(2,2)=fM[kZZ];
   }

   /**
      Set the components from nine scalars -- UNCHECKED for orthonormaility
   */
   void
   SetComponents (Scalar  xx, Scalar  xy, Scalar  xz,
                  Scalar  yx, Scalar  yy, Scalar  yz,
                  Scalar  zx, Scalar  zy, Scalar  zz) {
      fM[kXX]=xx;  fM[kXY]=xy;  fM[kXZ]=xz;
      fM[kYX]=yx;  fM[kYY]=yy;  fM[kYZ]=yz;
      fM[kZX]=zx;  fM[kZY]=zy;  fM[kZZ]=zz;
   }

   /**
      Get the nine components into nine scalars
   */
   void
   GetComponents (Scalar &xx, Scalar &xy, Scalar &xz,
                  Scalar &yx, Scalar &yy, Scalar &yz,
                  Scalar &zx, Scalar &zy, Scalar &zz) const {
      xx=fM[kXX];  xy=fM[kXY];  xz=fM[kXZ];
      yx=fM[kYX];  yy=fM[kYY];  yz=fM[kYZ];
      zx=fM[kZX];  zy=fM[kZY];  zz=fM[kZZ];
   }

   // =========== operations ==============


   /**
      Rotation operation on a displacement vector in any coordinate system
   */
   template <class CoordSystem, class U>
   DisplacementVector3D<CoordSystem,U>
   operator() (const DisplacementVector3D<CoordSystem,U> & v) const {
      DisplacementVector3D< Cartesian3D<double>,U > xyz;
      xyz.SetXYZ( fM[kXX] * v.X() + fM[kXY] * v.Y() + fM[kXZ] * v.Z() ,
                  fM[kYX] * v.X() + fM[kYY] * v.Y() + fM[kYZ] * v.Z() ,
                  fM[kZX] * v.X() + fM[kZY] * v.Y() + fM[kZZ] * v.Z() );
      return  DisplacementVector3D<CoordSystem,U>( xyz );
   }

   /**
      Rotation operation on a position vector in any coordinate system
   */
   template <class CoordSystem, class U>
   PositionVector3D<CoordSystem,U>
   operator() (const PositionVector3D<CoordSystem,U> & v) const {
      DisplacementVector3D< Cartesian3D<double>,U > xyz(v);
      DisplacementVector3D< Cartesian3D<double>,U > rxyz = operator()(xyz);
      return PositionVector3D<CoordSystem,U> ( rxyz );
   }

   /**
      Rotation operation on a Lorentz vector in any spatial coordinate system
   */
   template <class CoordSystem>
   LorentzVector<CoordSystem>
   operator() (const LorentzVector<CoordSystem> & v) const {
      DisplacementVector3D< Cartesian3D<double> > xyz(v.Vect());
      xyz = operator()(xyz);
      LorentzVector< PxPyPzE4D<double> > xyzt (xyz.X(), xyz.Y(), xyz.Z(), v.E());
      return LorentzVector<CoordSystem> ( xyzt );
   }

   /**
      Rotation operation on an arbitrary vector v.
      Preconditions:  v must implement methods x(), y(), and z()
      and the arbitrary vector type must have a constructor taking (x,y,z)
   */
   template <class ForeignVector>
   ForeignVector
   operator() (const  ForeignVector & v) const {
      DisplacementVector3D< Cartesian3D<double> > xyz(v);
      DisplacementVector3D< Cartesian3D<double> > rxyz = operator()(xyz);
      return ForeignVector ( rxyz.X(), rxyz.Y(), rxyz.Z() );
   }

   /**
      Overload operator * for rotation on a vector
   */
   template <class AVector>
   inline
   AVector operator* (const AVector & v) const
   {
      return operator()(v);
   }

   /**
      Invert a rotation in place
   */
   void Invert();

   /**
      Return inverse of  a rotation
   */
   Rotation3D Inverse() const { Rotation3D t(*this); t.Invert(); return t; }

   // ========= Multi-Rotation Operations ===============

   /**
      Multiply (combine) two rotations
   */
   Rotation3D operator * (const Rotation3D  & r) const {
   return Rotation3D
   (  fM[kXX]*r.fM[kXX] + fM[kXY]*r.fM[kYX] + fM[kXZ]*r.fM[kZX]
    , fM[kXX]*r.fM[kXY] + fM[kXY]*r.fM[kYY] + fM[kXZ]*r.fM[kZY]
    , fM[kXX]*r.fM[kXZ] + fM[kXY]*r.fM[kYZ] + fM[kXZ]*r.fM[kZZ]

    , fM[kYX]*r.fM[kXX] + fM[kYY]*r.fM[kYX] + fM[kYZ]*r.fM[kZX]
    , fM[kYX]*r.fM[kXY] + fM[kYY]*r.fM[kYY] + fM[kYZ]*r.fM[kZY]
    , fM[kYX]*r.fM[kXZ] + fM[kYY]*r.fM[kYZ] + fM[kYZ]*r.fM[kZZ]

    , fM[kZX]*r.fM[kXX] + fM[kZY]*r.fM[kYX] + fM[kZZ]*r.fM[kZX]
    , fM[kZX]*r.fM[kXY] + fM[kZY]*r.fM[kYY] + fM[kZZ]*r.fM[kZY]
    , fM[kZX]*r.fM[kXZ] + fM[kZY]*r.fM[kYZ] + fM[kZZ]*r.fM[kZZ]   );

   }


   /**
      Multiplication with arbitrary rotations
    */
    // note: cannot have a  template method since it is ambigous with the operator * on vectors

   Rotation3D operator * (const AxisAngle   & a) const;
   Rotation3D operator * (const EulerAngles & e) const;
   Rotation3D operator * (const Quaternion  & q) const;
   Rotation3D operator * (const RotationZYX & r) const;
   Rotation3D operator * (const RotationX  & rx) const;
   Rotation3D operator * (const RotationY  & ry) const;
   Rotation3D operator * (const RotationZ  & rz) const;

   /**
      Post-Multiply (on right) by another rotation :  T = T*R
   */
   template <class R>
   Rotation3D & operator *= (const R & r) { return *this = (*this)*r; }

   /**
                    Equality/inequality operators
   */
   bool operator == (const Rotation3D & rhs) const {
      if( fM[0] != rhs.fM[0] )  return false;
      if( fM[1] != rhs.fM[1] )  return false;
      if( fM[2] != rhs.fM[2] )  return false;
      if( fM[3] != rhs.fM[3] )  return false;
      if( fM[4] != rhs.fM[4] )  return false;
      if( fM[5] != rhs.fM[5] )  return false;
      if( fM[6] != rhs.fM[6] )  return false;
      if( fM[7] != rhs.fM[7] )  return false;
      if( fM[8] != rhs.fM[8] )  return false;
      return true;
   }
   bool operator != (const Rotation3D & rhs) const {
      return ! operator==(rhs);
   }

private:

   Scalar fM[9];  // 9 elements (3x3 matrix) representing the rotation

};  // Rotation3D

// ============ Class Rotation3D ends here ============

/**
   Distance between two rotations
 */
template <class R>
inline
typename Rotation3D::Scalar
Distance ( const Rotation3D& r1, const R & r2) {return gv_detail::dist(r1,r2);}

/**
   Multiplication of an axial rotation by a Rotation3D
 */
Rotation3D operator* (RotationX const & r1, Rotation3D const & r2);
Rotation3D operator* (RotationY const & r1, Rotation3D const & r2);
Rotation3D operator* (RotationZ const & r1, Rotation3D const & r2);

/**
   Multiplication of an axial rotation by another axial Rotation
 */
Rotation3D operator* (RotationX const & r1, RotationY const & r2);
Rotation3D operator* (RotationX const & r1, RotationZ const & r2);

Rotation3D operator* (RotationY const & r1, RotationX const & r2);
Rotation3D operator* (RotationY const & r1, RotationZ const & r2);

Rotation3D operator* (RotationZ const & r1, RotationX const & r2);
Rotation3D operator* (RotationZ const & r1, RotationY const & r2);

/**
   Stream Output and Input
 */
  // TODO - I/O should be put in the manipulator form

std::ostream & operator<< (std::ostream & os, const Rotation3D & r);

} // namespace Math
} // namespace ROOT

#endif // ROOT_Math_GenVector_Rotation3D
