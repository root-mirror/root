// @(#)root/mathcore:$Id$
// Authors: W. Brown, M. Fischler, L. Moneta    2005

/**********************************************************************
  *                                                                    *
  * Copyright (c) 2005 , LCG ROOT MathLib Team                         *
  *                                                                    *
  *                                                                    *
  **********************************************************************/

// Header file for class AxisAngle
//
// Created by: Lorenzo Moneta  at Wed May 11 10:37:10 2005
//
// Last update: Wed May 11 10:37:10 2005
//
#ifndef ROOT_Math_GenVector_AxisAngle
#define ROOT_Math_GenVector_AxisAngle  1

#include "Math/GenVector/Rotation3D.h"
#include "Math/GenVector/DisplacementVector3D.h"
#include "Math/GenVector/PositionVector3D.h"
#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/3DConversions.h"
#include <algorithm>
#include <cassert>


namespace ROOT {
namespace Math {


//__________________________________________________________________________________________
   /**
      AxisAngle class describing rotation represented with direction axis (3D Vector) and an
      angle of rotation around that axis.

      @ingroup GenVector

   */
class AxisAngle {

public:

   typedef double Scalar;

   /**
      definition of vector axis
   */
   typedef DisplacementVector3D<Cartesian3D<Scalar> > AxisVector;


   /**
      Default constructor (axis is z and angle is zero)
   */
   AxisAngle() : fAxis(0,0,1), fAngle(0) { }

   /**
      Construct from a non-zero vector (x,y,z) and an angle.
      Precondition:  the Vector needs to implement x(), y(), z(), and unit()
   */
   template<class AnyVector>
   AxisAngle(const AnyVector & v, Scalar angle) :
      fAxis(v.unit()), fAngle(angle) { }

   /**
      Construct given a pair of pointers or iterators defining the
      beginning and end of an array of four Scalars, to be treated as
      the x, y, and z components of a unit axis vector, and the angle
      of rotation.
      Precondition:  The first three components are assumed to represent
      the rotation axis vector and the 4-th the rotation angle.
      The angle is assumed to be in the range (-pi,pi].
      The axis vector is automatically normalized to be a unit vector
   */
   template<class IT>
   AxisAngle(IT begin, IT end) { SetComponents(begin,end); }

   // The compiler-generated copy ctor, copy assignment, and dtor are OK.

   /**
      Re-adjust components to eliminate small deviations from the axis
      being a unit vector and angles out of the canonical range (-pi,pi]
   */
   void Rectify();

   // ======== Construction From other Rotation Forms ==================

   /**
      Construct from another supported rotation type (see gv_detail::convert )
   */
   template <class OtherRotation>
   explicit AxisAngle(const OtherRotation & r) {gv_detail::convert(r,*this);}


   /**
      Assign from another supported rotation type (see gv_detail::convert )
   */
   template <class OtherRotation>
   AxisAngle & operator=( OtherRotation const  & r ) {
      gv_detail::convert(r,*this);
      return *this;
   }

   // ======== Components ==============

   /**
      Set the axis and then the angle given a pair of pointers or iterators
      defining the beginning and end of an array of four Scalars.
      Precondition:  The first three components are assumed to represent
      the rotation axis vector and the 4-th the rotation angle.
      The angle is assumed to be in the range (-pi,pi].
      The axis vector is automatically normalized to be a unit vector
   */
   template<class IT>
   void SetComponents(IT begin, IT end) {
      IT a = begin; IT b = ++begin; IT c = ++begin;
      fAxis.SetCoordinates(*a,*b,*c);
      fAngle = *(++begin);
      (void)end;
      assert (++begin==end);
      // re-normalize the vector
      double tot = fAxis.R();
      if (tot >  0) fAxis /= tot;
   }

   /**
      Get the axis and then the angle into data specified by an iterator begin
      and another to the end of the desired data (4 past start).
   */
   template<class IT>
   void GetComponents(IT begin, IT end) const {
      IT a = begin; IT b = ++begin; IT c = ++begin;
      fAxis.GetCoordinates(*a,*b,*c);
      *(++begin) = fAngle;
      (void)end;
      assert (++begin==end);
   }

   /**
      Get the axis and then the angle into data specified by an iterator begin
   */
   template<class IT>
   void GetComponents(IT begin) const {
      double ax,ay,az = 0;
      fAxis.GetCoordinates(ax,ay,az);
      *begin++ = ax;
      *begin++ = ay;
      *begin++ = az;
      *begin = fAngle;
   }

   /**
      Set components from a non-zero vector (x,y,z) and an angle.
      Precondition:  the Vector needs to implement x(), y(), z(), and unit()
   */
   template<class AnyVector>
   void SetComponents(const AnyVector & v, Scalar angle) {
      fAxis=v.unit();
      fAngle=angle;
   }

   /**
      Set components into a non-zero vector (x,y,z) and an angle.
      The vector is intended to be a cartesian dispalcement vector
      but any vector class assignable from one will work.
   */
   template<class AnyVector>
   void GetComponents(AnyVector & axis, Scalar & angle) const {
      axis  = fAxis;
      angle = fAngle;
   }

   /**
      accesss to rotation axis
   */
   AxisVector Axis() const { return fAxis; }

   /**
      access to rotation angle
   */
   Scalar Angle() const { return fAngle; }

   // =========== operations ==============

   /**
      Rotation operation on a cartesian vector
   */
   typedef  DisplacementVector3D<Cartesian3D<double>, DefaultCoordinateSystemTag > XYZVector;
   XYZVector operator() (const XYZVector & v) const;

   /**
      Rotation operation on a displacement vector in any coordinate system
   */
   template <class CoordSystem, class Tag>
   DisplacementVector3D<CoordSystem, Tag>
   operator() (const DisplacementVector3D<CoordSystem, Tag> & v) const {
      DisplacementVector3D< Cartesian3D<double> > xyz(v.X(), v.Y(), v.Z());
      DisplacementVector3D< Cartesian3D<double> > rxyz = operator()(xyz);
      DisplacementVector3D< CoordSystem, Tag > vNew;
      vNew.SetXYZ( rxyz.X(), rxyz.Y(), rxyz.Z() );
      return vNew;
   }

   /**
      Rotation operation on a position vector in any coordinate system
   */
   template <class CoordSystem, class Tag>
   PositionVector3D<CoordSystem, Tag>
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
      Invert an AxisAngle rotation in place
   */
   void Invert() { fAngle = -fAngle; }

   /**
      Return inverse of an AxisAngle rotation
   */
   AxisAngle Inverse() const { AxisAngle result(*this); result.Invert(); return result; }

   // ========= Multi-Rotation Operations ===============

   /**
      Multiply (combine) two rotations
   */
   AxisAngle operator * (const Rotation3D  & r) const;
   AxisAngle operator * (const AxisAngle   & a) const;
   AxisAngle operator * (const EulerAngles & e) const;
   AxisAngle operator * (const Quaternion  & q) const;
   AxisAngle operator * (const RotationZYX & r) const;
   AxisAngle operator * (const RotationX  & rx) const;
   AxisAngle operator * (const RotationY  & ry) const;
   AxisAngle operator * (const RotationZ  & rz) const;

   /**
      Post-Multiply (on right) by another rotation :  T = T*R
   */
   template <class R>
   AxisAngle & operator *= (const R & r) { return *this = (*this)*r; }


   /**
      Distance between two rotations
   */
   template <class R>
   Scalar Distance ( const R & r ) const {return gv_detail::dist(*this,r);}

   /**
      Equality/inequality operators
   */
   bool operator == (const AxisAngle & rhs) const {
      if( fAxis  != rhs.fAxis  )  return false;
      if( fAngle != rhs.fAngle )  return false;
      return true;
   }
   bool operator != (const AxisAngle & rhs) const {
      return ! operator==(rhs);
   }

private:

   AxisVector  fAxis;      // rotation axis (3D vector)
   Scalar      fAngle;     // rotation angle

   void RectifyAngle();

   static double Pi() { return 3.14159265358979323; }

};  // AxisAngle

// ============ Class AxisAngle ends here ============

/**
   Distance between two rotations
 */
template <class R>
inline
typename AxisAngle::Scalar
Distance ( const AxisAngle& r1, const R & r2) {return gv_detail::dist(r1,r2);}

/**
   Multiplication of an axial rotation by an AxisAngle
 */
AxisAngle operator* (RotationX const & r1, AxisAngle const & r2);
AxisAngle operator* (RotationY const & r1, AxisAngle const & r2);
AxisAngle operator* (RotationZ const & r1, AxisAngle const & r2);

/**
   Stream Output and Input
 */
  // TODO - I/O should be put in the manipulator form

std::ostream & operator<< (std::ostream & os, const AxisAngle & a);

} // namespace Math
} // namespace ROOT


#endif /* ROOT_Math_GenVector_AxisAngle  */
