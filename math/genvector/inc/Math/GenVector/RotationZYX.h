// @(#)root/mathcore:$Id$
// Authors: J. Palacios, L. Moneta    2007

 /**********************************************************************
  *                                                                    *
  * Copyright (c) 2007 , LCG ROOT MathLib Team                         *
  *                                                                    *
  *                                                                    *
  **********************************************************************/

// Header file for class Rotation in 3 dimensions, described by 3 Z-Y-X  Euler angles
// representing a rotation along Z, Y and X
//
// Created by: Lorenzo Moneta, Wed. May 22, 2007
//
// Last update: $Id$
//
#ifndef ROOT_Math_GenVector_RotationZYX
#define ROOT_Math_GenVector_RotationZYX  1

#include "Math/Math.h"

#include "Math/GenVector/Rotation3D.h"


#include "Math/GenVector/DisplacementVector3D.h"

#include "Math/GenVector/PositionVector3D.h"

#include "Math/GenVector/LorentzVector.h"

#include "Math/GenVector/3DConversions.h"


#include <algorithm>
#include <cassert>
#include <iostream>


namespace ROOT {
namespace Math {


//__________________________________________________________________________________________
  /**
     Rotation class with the (3D) rotation represented by
     angles describing first a rotation of
     an angle phi (yaw) about the  Z axis,
     followed by a rotation of an angle theta (pitch) about the Y axis,
     followed by a third rotation of an angle psi (roll) about the X axis.
     Note that the rotations are extrinsic rotations happening around a fixed coordinate system. 
     This is  different than the convention of the ROOT::Math::EulerAngles class, where the rotation are intrinsic. 
     Also it has not to be confused with the typical Goldstein definition of the Euler Angles
     (Z-X-Z or 313 sequence) which is used by the ROOT::Math::EulerAngles class, while the sequence here is Z-Y-X or 321.
     Applying a RotationZYX(phi, theta, psi)  to a vector is then equal to applying RotationX(psi) * RotationY(theta) * RotationZ(phi) to the same vector. 


     @ingroup GenVector
  */

class RotationZYX {

public:

   typedef double Scalar;


   // ========== Constructors and Assignment =====================

   /**
      Default constructor
   */
   RotationZYX() : fPhi(0.0), fTheta(0.0), fPsi(0.0) { }

   /**
      Constructor from phi, theta and psi
   */
   RotationZYX( Scalar phi, Scalar theta, Scalar psi ) :
      fPhi(phi), fTheta(theta), fPsi(psi)
   {Rectify();}  // Added 27 Jan. 06   JMM

   /**
      Construct given a pair of pointers or iterators defining the
      beginning and end of an array of three Scalars, to be treated as
      the angles phi, theta and psi.
   */
   template<class IT>
   RotationZYX(IT begin, IT end) { SetComponents(begin,end); }

   // The compiler-generated copy ctor, copy assignment, and dtor are OK.

   /**
      Re-adjust components place angles in canonical ranges
   */
   void Rectify();


   // ======== Construction and Assignment From other Rotation Forms ==================

   /**
      Construct from another supported rotation type (see gv_detail::convert )
   */
   template <class OtherRotation>
   explicit RotationZYX(const OtherRotation & r) {gv_detail::convert(r,*this);}


   /**
      Assign from another supported rotation type (see gv_detail::convert )
   */
   template <class OtherRotation>
   RotationZYX & operator=( OtherRotation const  & r ) {
      gv_detail::convert(r,*this);
      return *this;
   }


   // ======== Components ==============

   /**
      Set the three Euler angles given a pair of pointers or iterators
      defining the beginning and end of an array of three Scalars.
   */
   template<class IT>
   void SetComponents(IT begin, IT end) {
      fPhi   = *begin++;
      fTheta = *begin++;
      fPsi   = *begin++;
      (void)end;
      assert(begin == end);
      Rectify();
   }

   /**
      Get the axis and then the angle into data specified by an iterator begin
      and another to the end of the desired data (4 past start).
   */
   template<class IT>
   void GetComponents(IT begin, IT end) const {
      *begin++ = fPhi;
      *begin++ = fTheta;
      *begin++ = fPsi;
      (void)end;
      assert(begin == end);
   }

   /**
      Get the axis and then the angle into data specified by an iterator begin
   */
   template<class IT>
   void GetComponents(IT begin) const {
      *begin++ = fPhi;
      *begin++ = fTheta;
      *begin   = fPsi;
   }

   /**
      Set the components phi, theta, psi based on three Scalars.
   */
   void SetComponents(Scalar phi, Scalar theta, Scalar psi) {
      fPhi=phi; fTheta=theta; fPsi=psi;
      Rectify();
   }

   /**
      Get the components phi, theta, psi into three Scalars.
   */
   void GetComponents(Scalar & phi, Scalar & theta, Scalar & psi) const {
      phi=fPhi; theta=fTheta; psi=fPsi;
   }

   /**
      Set Phi angle (Z rotation angle)
   */
   void SetPhi(Scalar phi) { fPhi=phi; Rectify(); }

   /**
      Return Phi angle (Z rotation angle)
   */
   Scalar Phi() const { return fPhi; }

   /**
      Set Theta angle (Y' rotation angle)
   */
   void SetTheta(Scalar theta) { fTheta=theta; Rectify(); }

   /**
      Return Theta angle (Y' rotation angle)
   */
   Scalar Theta() const { return fTheta; }

   /**
      Set Psi angle (X'' rotation angle)
   */
   void SetPsi(Scalar psi) { fPsi=psi; Rectify(); }

   /**
      Return Psi angle (X'' rotation angle)
   */
   Scalar Psi() const { return fPsi; }

   // =========== operations ==============


   /**
      Rotation operation on a displacement vector in any coordinate system and tag
   */
   template <class CoordSystem, class U>
   DisplacementVector3D<CoordSystem,U>
   operator() (const DisplacementVector3D<CoordSystem,U> & v) const {
      return Rotation3D(*this) ( v );
   }

   /**
      Rotation operation on a position vector in any coordinate system
   */
   template <class CoordSystem, class U>
   PositionVector3D<CoordSystem, U>
   operator() (const PositionVector3D<CoordSystem,U> & v) const {
      DisplacementVector3D< Cartesian3D<double>,U > xyz(v);
      DisplacementVector3D< Cartesian3D<double>,U > rxyz = operator()(xyz);
      return PositionVector3D<CoordSystem,U> ( rxyz );
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
   void Invert();

   /**
      Return inverse of a rotation
   */
   RotationZYX Inverse() const {
      RotationZYX r(*this);
      r.Invert();
      return r;
   }


   // ========= Multi-Rotation Operations ===============

   /**
      Multiply (combine) two rotations
   */
   RotationZYX operator * (const RotationZYX & e) const;
   RotationZYX operator * (const Rotation3D  & r) const;
   RotationZYX operator * (const AxisAngle   & a) const;
   RotationZYX operator * (const Quaternion  & q) const;
   RotationZYX operator * (const EulerAngles & q) const;
   RotationZYX operator * (const RotationX  & rx) const;
   RotationZYX operator * (const RotationY  & ry) const;
   RotationZYX operator * (const RotationZ  & rz) const;

   /**
      Post-Multiply (on right) by another rotation :  T = T*R
   */
   template <class R>
   RotationZYX & operator *= (const R & r) { return *this = (*this)*r; }

   /**
      Distance between two rotations
   */
   template <class R>
   Scalar Distance ( const R & r ) const {return gv_detail::dist(*this,r);}

   /**
      Equality/inequality operators
   */
   bool operator == (const RotationZYX & rhs) const {
      if( fPhi   != rhs.fPhi   ) return false;
      if( fTheta != rhs.fTheta ) return false;
      if( fPsi   != rhs.fPsi   ) return false;
      return true;
   }
   bool operator != (const RotationZYX & rhs) const {
      return ! operator==(rhs);
   }

private:

   double fPhi;      // Z rotation angle (yaw)    defined in (-PI,PI]
   double fTheta;    // Y' rotation angle (pitch) defined in [-PI/2,PI/2]
   double fPsi;      // X'' rotation angle (roll) defined in (-PI,PI]

   static double Pi() { return M_PI; }

};  // RotationZYX

/**
   Distance between two rotations
 */
template <class R>
inline
typename RotationZYX::Scalar
Distance ( const RotationZYX& r1, const R & r2) {return gv_detail::dist(r1,r2);}

/**
   Multiplication of an axial rotation by an AxisAngle
 */
RotationZYX operator* (RotationX const & r1, RotationZYX const & r2);
RotationZYX operator* (RotationY const & r1, RotationZYX const & r2);
RotationZYX operator* (RotationZ const & r1, RotationZYX const & r2);

/**
   Stream Output and Input
 */
  // TODO - I/O should be put in the manipulator form

std::ostream & operator<< (std::ostream & os, const RotationZYX & e);


} // namespace Math
} // namespace ROOT

#endif // ROOT_Math_GenVector_RotationZYX
