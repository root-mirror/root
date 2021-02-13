// @(#)root/mathcore:$Id$
// Authors: W. Brown, M. Fischler, L. Moneta    2005

 /**********************************************************************
  *                                                                    *
  * Copyright (c) 2005 ROOT MathLib Team                               *
  *                                                                    *
  *                                                                    *
  **********************************************************************/

// Header file for LorentzRotation
//
// Created by: Mark Fischler  Mon Aug 8  2005
//
// Last update: $Id$
//
#ifndef ROOT_Math_GenVector_LorentzRotation
#define ROOT_Math_GenVector_LorentzRotation  1

#include "Math/GenVector/LorentzRotationfwd.h"

#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/PxPyPzE4D.h"

#include "Math/GenVector/Rotation3Dfwd.h"
#include "Math/GenVector/AxisAnglefwd.h"
#include "Math/GenVector/EulerAnglesfwd.h"
#include "Math/GenVector/Quaternionfwd.h"
#include "Math/GenVector/RotationXfwd.h"
#include "Math/GenVector/RotationYfwd.h"
#include "Math/GenVector/RotationZfwd.h"
#include "Math/GenVector/Boost.h"
#include "Math/GenVector/BoostX.h"
#include "Math/GenVector/BoostY.h"
#include "Math/GenVector/BoostZ.h"

namespace ROOT {

  namespace Math {

//__________________________________________________________________________________________
  /**
     Lorentz transformation class with the (4D) transformation represented by
     a 4x4 orthosymplectic matrix.
     See also Boost, BoostX, BoostY and BoostZ for classes representing
     specialized Lorentz transformations.
     Also, the 3-D rotation classes can be considered to be special Lorentz
     transformations which do not mix space and time components.

     @ingroup GenVector

  */

class LorentzRotation {

public:

  typedef double Scalar;

  enum ELorentzRotationMatrixIndex {
      kXX =  0, kXY =  1, kXZ =  2, kXT =  3
    , kYX =  4, kYY =  5, kYZ =  6, kYT =  7
    , kZX =  8, kZY =  9, kZZ = 10, kZT = 11
    , kTX = 12, kTY = 13, kTZ = 14, kTT = 15
  };

  // ========== Constructors and Assignment =====================

  /**
      Default constructor (identity transformation)
  */
  LorentzRotation();

  /**
     Construct given a pair of pointers or iterators defining the
     beginning and end of an array of sixteen Scalars
   */
  template<class IT>
  LorentzRotation(IT begin, IT end) { SetComponents(begin,end); }

  // The compiler-generated and dtor are OK but we have implementwd the copy-ctor and
  // assignment operators since we have a template assignment

  /**
     Copy constructor
   */
   LorentzRotation( LorentzRotation const & r ) {
      *this = r;
   }

  /**
     Construct from a pure boost
  */
  explicit LorentzRotation( Boost  const & b  ) {  b.GetLorentzRotation( fM+0 ); }
  explicit LorentzRotation( BoostX const & bx ) { bx.GetLorentzRotation( fM+0 ); }
  explicit LorentzRotation( BoostY const & by ) { by.GetLorentzRotation( fM+0 ); }
  explicit LorentzRotation( BoostZ const & bz ) { bz.GetLorentzRotation( fM+0 ); }

  /**
     Construct from a 3-D rotation (no space-time mixing)
  */
  explicit LorentzRotation( Rotation3D  const & r );
  explicit LorentzRotation( AxisAngle   const & a );
  explicit LorentzRotation( EulerAngles const & e );
  explicit LorentzRotation( Quaternion  const & q );
  explicit LorentzRotation( RotationX   const & r );
  explicit LorentzRotation( RotationY   const & r );
  explicit LorentzRotation( RotationZ   const & r );

  /**
     Construct from a linear algebra matrix of size at least 4x4,
     which must support operator()(i,j) to obtain elements (0,3) thru (3,3).
     Precondition:  The matrix is assumed to be orthosymplectic.  NO checking
     or re-adjusting is performed.
     Note:  (0,0) refers to the XX component; (3,3) refers to the TT component.
  */
  template<class ForeignMatrix>
  explicit LorentzRotation(const ForeignMatrix & m) { SetComponents(m); }

  /**
     Construct from four orthosymplectic vectors (which must have methods
     x(), y(), z() and t()) which will be used as the columns of the Lorentz
     rotation matrix.  The orthosymplectic conditions will be checked, and
     values adjusted so that the result will always be a good Lorentz rotation
     matrix.
  */
  template<class Foreign4Vector>
  LorentzRotation(const Foreign4Vector& v1,
                  const Foreign4Vector& v2,
                  const Foreign4Vector& v3,
                  const Foreign4Vector& v4 ) { SetComponents(v1, v2, v3, v4); }


  /**
     Raw constructor from sixteen Scalar components (without any checking)
  */
  LorentzRotation(Scalar  xx, Scalar  xy, Scalar  xz, Scalar xt,
                  Scalar  yx, Scalar  yy, Scalar  yz, Scalar yt,
                  Scalar  zx, Scalar  zy, Scalar  zz, Scalar zt,
                  Scalar  tx, Scalar  ty, Scalar  tz, Scalar tt)
 {
    SetComponents (xx, xy, xz, xt,
                   yx, yy, yz, yt,
                   zx, zy, zz, zt,
                   tx, ty, tz, tt);
 }

  /**
      Assign from another LorentzRotation
  */
   LorentzRotation &
  operator=( LorentzRotation  const & rhs ) {
      SetComponents( rhs.fM[0],  rhs.fM[1],  rhs.fM[2],  rhs.fM[3],
                     rhs.fM[4],  rhs.fM[5],  rhs.fM[6],  rhs.fM[7],
                     rhs.fM[8],  rhs.fM[9],  rhs.fM[10], rhs.fM[11],
                     rhs.fM[12], rhs.fM[13], rhs.fM[14], rhs.fM[15] );
      return *this;
   }

  /**
     Assign from a pure boost
  */
  LorentzRotation &
  operator=( Boost  const & b ) { return operator=(LorentzRotation(b)); }
  LorentzRotation &
  operator=( BoostX const & b ) { return operator=(LorentzRotation(b)); }
  LorentzRotation &
  operator=( BoostY const & b ) { return operator=(LorentzRotation(b)); }
  LorentzRotation &
  operator=( BoostZ const & b ) { return operator=(LorentzRotation(b)); }

  /**
     Assign from a 3-D rotation
  */
  LorentzRotation &
  operator=( Rotation3D  const & r ) { return operator=(LorentzRotation(r)); }
  LorentzRotation &
  operator=( AxisAngle   const & a ) { return operator=(LorentzRotation(a)); }
  LorentzRotation &
  operator=( EulerAngles const & e ) { return operator=(LorentzRotation(e)); }
  LorentzRotation &
  operator=( Quaternion  const & q ) { return operator=(LorentzRotation(q)); }
  LorentzRotation &
  operator=( RotationZ   const & r ) { return operator=(LorentzRotation(r)); }
  LorentzRotation &
  operator=( RotationY   const & r ) { return operator=(LorentzRotation(r)); }
  LorentzRotation &
  operator=( RotationX   const & r ) { return operator=(LorentzRotation(r)); }

  /**
     Assign from a linear algebra matrix of size at least 4x4,
     which must support operator()(i,j) to obtain elements (0,3) thru (3,3).
     Precondition:  The matrix is assumed to be orthosymplectic.  NO checking
     or re-adjusting is performed.
  */
  template<class ForeignMatrix>
  LorentzRotation &
  operator=(const ForeignMatrix & m) {
     SetComponents( m(0,0), m(0,1), m(0,2), m(0,3),
                    m(1,0), m(1,1), m(1,2), m(1,3),
                    m(2,0), m(2,1), m(2,2), m(2,3),
                    m(3,0), m(3,1), m(3,2), m(3,3) );
     return *this;
  }

  /**
     Re-adjust components to eliminate small deviations from a perfect
     orthosyplectic matrix.
   */
  void Rectify();

  // ======== Components ==============

  /**
     Set components from four orthosymplectic vectors (which must have methods
     x(), y(), z(), and t()) which will be used as the columns of the
     Lorentz rotation matrix.  The values will be adjusted
     so that the result will always be a good Lorentz rotation matrix.
  */
  template<class Foreign4Vector>
  void
  SetComponents (const Foreign4Vector& v1,
                 const Foreign4Vector& v2,
                 const Foreign4Vector& v3,
                 const Foreign4Vector& v4 ) {
    fM[kXX]=v1.x();  fM[kXY]=v2.x();  fM[kXZ]=v3.x();  fM[kXT]=v4.x();
    fM[kYX]=v1.y();  fM[kYY]=v2.y();  fM[kYZ]=v3.y();  fM[kYT]=v4.y();
    fM[kZX]=v1.z();  fM[kZY]=v2.z();  fM[kZZ]=v3.z();  fM[kZT]=v4.z();
    fM[kTX]=v1.t();  fM[kTY]=v2.t();  fM[kTZ]=v3.t();  fM[kTT]=v4.t();
    Rectify();
  }

  /**
     Get components into four 4-vectors which will be the (orthosymplectic)
     columns of the rotation matrix.  (The 4-vector class must have a
     constructor from 4 Scalars used as x, y, z, t)
  */
  template<class Foreign4Vector>
  void
  GetComponents ( Foreign4Vector& v1,
                  Foreign4Vector& v2,
                  Foreign4Vector& v3,
                  Foreign4Vector& v4 ) const {
    v1 = Foreign4Vector ( fM[kXX], fM[kYX], fM[kZX], fM[kTX] );
    v2 = Foreign4Vector ( fM[kXY], fM[kYY], fM[kZY], fM[kTY] );
    v3 = Foreign4Vector ( fM[kXZ], fM[kYZ], fM[kZZ], fM[kTZ] );
    v4 = Foreign4Vector ( fM[kXT], fM[kYT], fM[kZT], fM[kTT] );
  }

  /**
     Set the 16 matrix components given an iterator to the start of
     the desired data, and another to the end (16 past start).
   */
  template<class IT>
  void SetComponents(IT begin, IT end) {
     for (int i = 0; i <16; ++i) {
        fM[i] = *begin;
        ++begin;
     }
     (void)end;
     assert (end==begin);
  }

  /**
     Get the 16 matrix components into data specified by an iterator begin
     and another to the end of the desired data (16 past start).
   */
  template<class IT>
  void GetComponents(IT begin, IT end) const {
     for (int i = 0; i <16; ++i) {
        *begin = fM[i];
        ++begin;
     }
     (void)end;
     assert (end==begin);
  }

  /**
     Get the 16 matrix components into data specified by an iterator begin
   */
  template<class IT>
  void GetComponents(IT begin) const {
    std::copy ( fM+0, fM+16, begin );
  }

  /**
     Set components from a linear algebra matrix of size at least 4x4,
     which must support operator()(i,j) to obtain elements (0,0) thru (3,3).
     Precondition:  The matrix is assumed to be orthosymplectic.  NO checking
     or re-adjusting is performed.
  */
  template<class ForeignMatrix>
  void
  SetRotationMatrix (const ForeignMatrix & m) {
    fM[kXX]=m(0,0);  fM[kXY]=m(0,1);  fM[kXZ]=m(0,2);  fM[kXT]=m(0,3);
    fM[kYX]=m(1,0);  fM[kYY]=m(1,1);  fM[kYZ]=m(1,2);  fM[kYT]=m(1,3);
    fM[kZX]=m(2,0);  fM[kZY]=m(2,1);  fM[kZZ]=m(2,2);  fM[kZT]=m(2,3);
    fM[kTX]=m(3,0);  fM[kTY]=m(3,1);  fM[kTZ]=m(3,2);  fM[kTT]=m(3,3);
  }

  /**
     Get components into a linear algebra matrix of size at least 4x4,
     which must support operator()(i,j) for write access to elements
     (0,0) thru (3,3).
  */
  template<class ForeignMatrix>
  void
  GetRotationMatrix (ForeignMatrix & m) const {
    m(0,0)=fM[kXX];  m(0,1)=fM[kXY];  m(0,2)=fM[kXZ]; m(0,3)=fM[kXT];
    m(1,0)=fM[kYX];  m(1,1)=fM[kYY];  m(1,2)=fM[kYZ]; m(1,3)=fM[kYT];
    m(2,0)=fM[kZX];  m(2,1)=fM[kZY];  m(2,2)=fM[kZZ]; m(2,3)=fM[kZT];
    m(3,0)=fM[kTX];  m(3,1)=fM[kTY];  m(3,2)=fM[kTZ]; m(3,3)=fM[kTT];
  }

  /**
     Set the components from sixteen scalars -- UNCHECKED for orthosymplectic
   */
  void
  SetComponents (Scalar  xx, Scalar  xy, Scalar  xz, Scalar  xt,
                 Scalar  yx, Scalar  yy, Scalar  yz, Scalar  yt,
                 Scalar  zx, Scalar  zy, Scalar  zz, Scalar  zt,
                 Scalar  tx, Scalar  ty, Scalar  tz, Scalar  tt) {
                 fM[kXX]=xx;  fM[kXY]=xy;  fM[kXZ]=xz;  fM[kXT]=xt;
                 fM[kYX]=yx;  fM[kYY]=yy;  fM[kYZ]=yz;  fM[kYT]=yt;
                 fM[kZX]=zx;  fM[kZY]=zy;  fM[kZZ]=zz;  fM[kZT]=zt;
                 fM[kTX]=tx;  fM[kTY]=ty;  fM[kTZ]=tz;  fM[kTT]=tt;
  }

  /**
     Get the sixteen components into sixteen scalars
   */
  void
  GetComponents (Scalar &xx, Scalar &xy, Scalar &xz, Scalar &xt,
                 Scalar &yx, Scalar &yy, Scalar &yz, Scalar &yt,
                 Scalar &zx, Scalar &zy, Scalar &zz, Scalar &zt,
                 Scalar &tx, Scalar &ty, Scalar &tz, Scalar &tt) const {
                 xx=fM[kXX];  xy=fM[kXY];  xz=fM[kXZ];  xt=fM[kXT];
                 yx=fM[kYX];  yy=fM[kYY];  yz=fM[kYZ];  yt=fM[kYT];
                 zx=fM[kZX];  zy=fM[kZY];  zz=fM[kZZ];  zt=fM[kZT];
                 tx=fM[kTX];  ty=fM[kTY];  tz=fM[kTZ];  tt=fM[kTT];
  }

  // =========== operations ==============

  /**
     Lorentz transformation operation on a Minkowski ('Cartesian')
     LorentzVector
  */
  LorentzVector< ROOT::Math::PxPyPzE4D<double> >
  operator() (const LorentzVector< ROOT::Math::PxPyPzE4D<double> > & v) const {
        Scalar x = v.Px();
        Scalar y = v.Py();
        Scalar z = v.Pz();
        Scalar t = v.E();
        return LorentzVector< PxPyPzE4D<double> >
           ( fM[kXX]*x + fM[kXY]*y + fM[kXZ]*z + fM[kXT]*t
             , fM[kYX]*x + fM[kYY]*y + fM[kYZ]*z + fM[kYT]*t
             , fM[kZX]*x + fM[kZY]*y + fM[kZZ]*z + fM[kZT]*t
             , fM[kTX]*x + fM[kTY]*y + fM[kTZ]*z + fM[kTT]*t );
  }

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
     Overload operator * for rotation on a vector
   */
  template <class A4Vector>
  inline
  A4Vector operator* (const A4Vector & v) const
  {
    return operator()(v);
  }

  /**
      Invert a Lorentz rotation in place
   */
  void Invert();

  /**
      Return inverse of  a rotation
   */
  LorentzRotation Inverse() const;

  // ========= Multi-Rotation Operations ===============

  /**
     Multiply (combine) this Lorentz rotation by another LorentzRotation
   */
  LorentzRotation operator * (const LorentzRotation & r) const;

  //#ifdef TODO_LATER
  /**
     Multiply (combine) this Lorentz rotation by a pure Lorentz boost
   */
  //TODO: implement directly in a more efficient way. Now are implemented
  // going through another LorentzRotation
  LorentzRotation operator * (const Boost  & b) const  { LorentzRotation tmp(b); return (*this)*tmp; }
  LorentzRotation operator * (const BoostX & b) const  { LorentzRotation tmp(b); return (*this)*tmp; }
  LorentzRotation operator * (const BoostY & b) const  { LorentzRotation tmp(b); return (*this)*tmp; }
  LorentzRotation operator * (const BoostZ & b) const  { LorentzRotation tmp(b); return (*this)*tmp; }

  /**
     Multiply (combine) this Lorentz rotation by a 3-D Rotation
   */
  LorentzRotation operator * (const Rotation3D  & r) const { LorentzRotation tmp(r); return (*this)*tmp; }
  LorentzRotation operator * (const AxisAngle   & a) const { LorentzRotation tmp(a); return (*this)*tmp; }
  LorentzRotation operator * (const EulerAngles & e) const { LorentzRotation tmp(e); return (*this)*tmp; }
  LorentzRotation operator * (const Quaternion  & q) const { LorentzRotation tmp(q); return (*this)*tmp; }
  LorentzRotation operator * (const RotationX  & rx) const { LorentzRotation tmp(rx); return (*this)*tmp; }
  LorentzRotation operator * (const RotationY  & ry) const { LorentzRotation tmp(ry); return (*this)*tmp; }
  LorentzRotation operator * (const RotationZ  & rz) const { LorentzRotation tmp(rz); return (*this)*tmp; }
  //#endif

  /**
     Post-Multiply (on right) by another LorentzRotation, Boost, or
     rotation :  T = T*R
   */
  template <class R>
  LorentzRotation & operator *= (const R & r) { return *this = (*this)*r; }

  /**
     Equality/inequality operators
   */
  bool operator == (const LorentzRotation & rhs) const {
    for (unsigned int i=0; i < 16; ++i) {
      if( fM[i] != rhs.fM[i] )  return false;
    }
    return true;
  }
  bool operator != (const LorentzRotation & rhs) const {
    return ! operator==(rhs);
  }

private:

  Scalar fM[16];

};  // LorentzRotation

// ============ Class LorentzRotation ends here ============


/**
   Stream Output and Input
 */
  // TODO - I/O should be put in the manipulator form

std::ostream & operator<< (std::ostream & os, const LorentzRotation & r);

// ============================================ vetted to here  ============

#ifdef NOTYET
/**
   Distance between two Lorentz rotations
 */
template <class R>
inline
typename Rotation3D::Scalar
Distance ( const Rotation3D& r1, const R & r2) {return gv_detail::dist(r1,r2);}
#endif

} //namespace Math
} //namespace ROOT







#endif /* ROOT_Math_GenVector_LorentzRotation  */
