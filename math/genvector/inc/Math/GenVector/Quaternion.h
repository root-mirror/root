// @(#)root/mathcore:$Id$
// Authors: W. Brown, M. Fischler, L. Moneta    2005

 /**********************************************************************
  *                                                                    *
  * Copyright (c) 2005 , LCG ROOT FNAL MathLib Team                    *
  *                                                                    *
  *                                                                    *
  **********************************************************************/

// Header file for rotation in 3 dimensions, represented by a quaternion
// Created by: Mark Fischler Thurs June 9  2005
//
// Last update: $Id$
//
#ifndef ROOT_Math_GenVector_Quaternion
#define ROOT_Math_GenVector_Quaternion  1


#include "Math/GenVector/Cartesian3D.h"
#include "Math/GenVector/DisplacementVector3D.h"
#include "Math/GenVector/PositionVector3D.h"
#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/3DConversions.h"
#include "Math/GenVector/3DDistances.h"

#include <algorithm>
#include <cassert>


namespace ROOT {
namespace Math {


//__________________________________________________________________________________________
   /**
      Rotation class with the (3D) rotation represented by
      a unit quaternion (u, i, j, k).
      This is the optimal representation for multiplication of multiple
      rotations, and for computation of group-manifold-invariant distance
      between two rotations.
      See also ROOT::Math::AxisAngle, ROOT::Math::EulerAngles, and ROOT::Math::Rotation3D.

      @ingroup GenVector
   */

class Quaternion {

public:

  typedef double Scalar;

  // ========== Constructors and Assignment =====================

  /**
      Default constructor (identity rotation)
  */
   Quaternion()
      : fU(1.0)
      , fI(0.0)
      , fJ(0.0)
      , fK(0.0)
   { }

   /**
      Construct given a pair of pointers or iterators defining the
      beginning and end of an array of four Scalars
   */
   template<class IT>
   Quaternion(IT begin, IT end) { SetComponents(begin,end); }

   // ======== Construction From other Rotation Forms ==================

   /**
      Construct from another supported rotation type (see gv_detail::convert )
   */
   template <class OtherRotation>
   explicit Quaternion(const OtherRotation & r) {gv_detail::convert(r,*this);}


   /**
      Construct from four Scalars representing the coefficients of u, i, j, k
   */
   Quaternion(Scalar u, Scalar i, Scalar j, Scalar k) :
      fU(u), fI(i), fJ(j), fK(k) { }

   // The compiler-generated copy ctor, copy assignment, and dtor are OK.

   /**
      Re-adjust components to eliminate small deviations from |Q| = 1
      orthonormality.
   */
   void Rectify();

   /**
      Assign from another supported rotation type (see gv_detail::convert )
   */
   template <class OtherRotation>
   Quaternion & operator=( OtherRotation const  & r ) {
      gv_detail::convert(r,*this);
      return *this;
   }

   // ======== Components ==============

   /**
      Set the four components given an iterator to the start of
      the desired data, and another to the end (4 past start).
   */
   template<class IT>
   void SetComponents(IT begin, IT end) {
      fU = *begin++;
      fI = *begin++;
      fJ = *begin++;
      fK = *begin++;
      (void)end;
      assert (end==begin);
   }

   /**
      Get the components into data specified by an iterator begin
      and another to the end of the desired data (4 past start).
   */
   template<class IT>
   void GetComponents(IT begin, IT end) const {
      *begin++ = fU;
      *begin++ = fI;
      *begin++ = fJ;
      *begin++ = fK;
      (void)end;
      assert (end==begin);
   }

   /**
      Get the components into data specified by an iterator begin
   */
   template<class IT>
   void GetComponents(IT begin ) const {
      *begin++ = fU;
      *begin++ = fI;
      *begin++ = fJ;
      *begin   = fK;
   }

   /**
      Set the components based on four Scalars.  The sum of the squares of
      these Scalars should be 1; no checking is done.
   */
   void SetComponents(Scalar u, Scalar i, Scalar j, Scalar k) {
      fU=u; fI=i; fJ=j; fK=k;
   }

   /**
      Get the components into four Scalars.
   */
   void GetComponents(Scalar & u, Scalar & i, Scalar & j, Scalar & k) const {
      u=fU; i=fI; j=fJ; k=fK;
   }

   /**
      Access to the four quaternion components:
      U() is the coefficient of the identity Pauli matrix,
      I(), J() and K() are the coefficients of sigma_x, sigma_y, sigma_z
   */
   Scalar U() const { return fU; }
   Scalar I() const { return fI; }
   Scalar J() const { return fJ; }
   Scalar K() const { return fK; }

   // =========== operations ==============

   /**
      Rotation operation on a cartesian vector
   */
   typedef  DisplacementVector3D<Cartesian3D<double>, DefaultCoordinateSystemTag > XYZVector;
   XYZVector operator() (const XYZVector & v) const {

      const Scalar alpha = fU*fU - fI*fI - fJ*fJ - fK*fK;
      const Scalar twoQv = 2*(fI*v.X() + fJ*v.Y() + fK*v.Z());
      const Scalar twoU  = 2 * fU;
      return XYZVector  (  alpha * v.X() + twoU * (fJ*v.Z() - fK*v.Y()) + twoQv * fI ,
                           alpha * v.Y() + twoU * (fK*v.X() - fI*v.Z()) + twoQv * fJ ,
                           alpha * v.Z() + twoU * (fI*v.Y() - fJ*v.X()) + twoQv * fK );
   }

   /**
      Rotation operation on a displacement vector in any coordinate system
   */
   template <class CoordSystem,class Tag>
   DisplacementVector3D<CoordSystem,Tag>
   operator() (const DisplacementVector3D<CoordSystem,Tag> & v) const {
      DisplacementVector3D< Cartesian3D<double> > xyz(v.X(), v.Y(), v.Z());
      DisplacementVector3D< Cartesian3D<double> > rxyz = operator()(xyz);
      DisplacementVector3D< CoordSystem,Tag > vNew;
      vNew.SetXYZ( rxyz.X(), rxyz.Y(), rxyz.Z() );
      return vNew;
   }

   /**
      Rotation operation on a position vector in any coordinate system
   */
   template <class CoordSystem, class Tag>
   PositionVector3D<CoordSystem,Tag>
   operator() (const PositionVector3D<CoordSystem,Tag> & p) const {
      DisplacementVector3D< Cartesian3D<double>,Tag > xyz(p);
      DisplacementVector3D< Cartesian3D<double>,Tag > rxyz = operator()(xyz);
      return PositionVector3D<CoordSystem,Tag> ( rxyz );
   }

   /**
      Rotation operation on a Lorentz vector in any 4D coordinate system
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
   void Invert() { fI = -fI; fJ = -fJ; fK = -fK; }

   /**
      Return inverse of a rotation
   */
   Quaternion Inverse() const { return Quaternion(fU, -fI, -fJ, -fK); }

   // ========= Multi-Rotation Operations ===============

   /**
      Multiply (combine) two rotations
   */
   /**
      Multiply (combine) two rotations
   */
   Quaternion operator * (const Quaternion  & q) const {
      return Quaternion  (   fU*q.fU - fI*q.fI - fJ*q.fJ - fK*q.fK ,
                             fU*q.fI + fI*q.fU + fJ*q.fK - fK*q.fJ ,
                             fU*q.fJ - fI*q.fK + fJ*q.fU + fK*q.fI ,
                             fU*q.fK + fI*q.fJ - fJ*q.fI + fK*q.fU  );
   }

   Quaternion operator * (const Rotation3D  & r) const;
   Quaternion operator * (const AxisAngle   & a) const;
   Quaternion operator * (const EulerAngles & e) const;
   Quaternion operator * (const RotationZYX & r) const;
   Quaternion operator * (const RotationX  & rx) const;
   Quaternion operator * (const RotationY  & ry) const;
   Quaternion operator * (const RotationZ  & rz) const;

   /**
      Post-Multiply (on right) by another rotation :  T = T*R
   */
   template <class R>
   Quaternion & operator *= (const R & r) { return *this = (*this)*r; }


   /**
      Distance between two rotations in Quaternion form
      Note:  The rotation group is isomorphic to a 3-sphere
      with diametrically opposite points identified.
      The (rotation group-invariant) is the smaller
      of the two possible angles between the images of
      the two totations on that sphere.  Thus the distance
      is never greater than pi/2.
   */

   Scalar Distance(const Quaternion & q) const ;

   /**
      Equality/inequality operators
   */
   bool operator == (const Quaternion & rhs) const {
      if( fU != rhs.fU )  return false;
      if( fI != rhs.fI )  return false;
      if( fJ != rhs.fJ )  return false;
      if( fK != rhs.fK )  return false;
      return true;
   }
   bool operator != (const Quaternion & rhs) const {
      return ! operator==(rhs);
   }

private:

   Scalar fU;
   Scalar fI;
   Scalar fJ;
   Scalar fK;

};  // Quaternion

// ============ Class Quaternion ends here ============

/**
   Distance between two rotations
 */
template <class R>
inline
typename Quaternion::Scalar
Distance ( const Quaternion& r1, const R & r2) {return gv_detail::dist(r1,r2);}

/**
   Multiplication of an axial rotation by an AxisAngle
 */
Quaternion operator* (RotationX const & r1, Quaternion const & r2);
Quaternion operator* (RotationY const & r1, Quaternion const & r2);
Quaternion operator* (RotationZ const & r1, Quaternion const & r2);

/**
   Stream Output and Input
 */
  // TODO - I/O should be put in the manipulator form

std::ostream & operator<< (std::ostream & os, const Quaternion & q);


}  // namespace Math
}  // namespace ROOT

#endif // ROOT_Math_GenVector_Quaternion
