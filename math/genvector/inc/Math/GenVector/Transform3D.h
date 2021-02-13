// @(#)root/mathcore:$Id$
// Authors: W. Brown, M. Fischler, L. Moneta    2005

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 , LCG ROOT MathLib Team                         *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Header file for class Transform3D
//
// Created by: Lorenzo Moneta  October 21 2005
//
//
#ifndef ROOT_Math_GenVector_Transform3D
#define ROOT_Math_GenVector_Transform3D  1



#include "Math/GenVector/DisplacementVector3D.h"

#include "Math/GenVector/PositionVector3D.h"

#include "Math/GenVector/Rotation3D.h"

#include "Math/GenVector/Translation3D.h"


#include "Math/GenVector/AxisAnglefwd.h"
#include "Math/GenVector/EulerAnglesfwd.h"
#include "Math/GenVector/Quaternionfwd.h"
#include "Math/GenVector/RotationZYXfwd.h"
#include "Math/GenVector/RotationXfwd.h"
#include "Math/GenVector/RotationYfwd.h"
#include "Math/GenVector/RotationZfwd.h"

#include <iostream>
#include <type_traits>
#include <cmath>

//#include "Math/Vector3Dfwd.h"



namespace ROOT {

namespace Math {

namespace Impl {

//_________________________________________________________________________________________
/**
    Basic 3D Transformation class describing  a rotation and then a translation
    The internal data are a 3D rotation data (represented as a 3x3 matrix) and a 3D vector data.
    They are represented and held in this class like a 3x4 matrix (a simple array of 12 numbers).

    The class can be constructed from any 3D rotation object
    (ROOT::Math::Rotation3D, ROOT::Math::AxisAngle, ROOT::Math::Quaternion, etc...) and/or
    a 3D Vector (ROOT::Math::DislacementVector3D or via ROOT::Math::Translation ) representing a Translation.
    The Transformation is defined by applying first the rotation and then the translation.
    A transformation defined by applying first a translation and then a rotation is equivalent to the
    transformation obtained applying first the rotation and then a translation equivalent to the rotated vector.
    The operator * can be used to obtain directly such transformations, in addition to combine various
    transformations.
    Keep in mind that the operator * (like in the case of rotations ) is not commutative.
    The operator * is used (in addition to operator() ) to apply a transformations on the vector
    (DisplacementVector3D and LorentzVector classes) and point (PositionVector3D)  classes.
    In the case of Vector objects the transformation only rotates them and does not translate them.
    Only Point objects are able to be both rotated and translated.


    @ingroup GenVector

*/

template <typename T = double>
class Transform3D {

public:
   typedef T Scalar;

   typedef DisplacementVector3D<Cartesian3D<T>, DefaultCoordinateSystemTag> Vector;
   typedef PositionVector3D<Cartesian3D<T>, DefaultCoordinateSystemTag>     Point;

   enum ETransform3DMatrixIndex {
      kXX = 0, kXY = 1, kXZ = 2, kDX = 3,
      kYX = 4, kYY = 5, kYZ = 6, kDY = 7,
      kZX = 8, kZY = 9, kZZ =10, kDZ = 11
   };



   /**
       Default constructor (identy rotation) + zero translation
   */
   Transform3D()
   {
      SetIdentity();
   }

   /**
      Construct given a pair of pointers or iterators defining the
      beginning and end of an array of 12 Scalars
   */
   template<class IT>
   Transform3D(IT begin, IT end)
   {
      SetComponents(begin,end);
   }

   /**
      Construct from a rotation and then a translation described by a Vector
   */
   Transform3D( const Rotation3D & r, const Vector & v)
   {
      AssignFrom( r, v );
   }
   /**
      Construct from a rotation and then a translation described by a Translation3D class
   */
   Transform3D(const Rotation3D &r, const Translation3D<T> &t) { AssignFrom(r, t.Vect()); }

   /**
      Construct from a rotation (any rotation object)  and then a translation
      (represented by any DisplacementVector)
      The requirements on the rotation and vector objects are that they can be transformed in a
      Rotation3D class and in a Cartesian3D Vector
   */
   template <class ARotation, class CoordSystem, class Tag>
   Transform3D( const ARotation & r, const DisplacementVector3D<CoordSystem,Tag> & v)
   {
      AssignFrom( Rotation3D(r), Vector (v.X(),v.Y(),v.Z()) );
   }

   /**
      Construct from a rotation (any rotation object)  and then a translation
      represented by a Translation3D class
      The requirements on the rotation is that it can be transformed in a
      Rotation3D class
   */
   template <class ARotation>
   Transform3D(const ARotation &r, const Translation3D<T> &t)
   {
      AssignFrom( Rotation3D(r), t.Vect() );
   }


#ifdef OLD_VERSION
   /**
      Construct from a translation and then a rotation (inverse assignment)
   */
   Transform3D( const Vector & v, const Rotation3D & r)
   {
      // is equivalent from having first the rotation and then the translation vector rotated
      AssignFrom( r, r(v) );
   }
#endif

   /**
      Construct from a 3D Rotation only with zero translation
   */
   explicit Transform3D( const Rotation3D & r) {
      AssignFrom(r);
   }

   // convenience methods for constructing a Transform3D from all the 3D rotations classes
   // (cannot use templates for conflict with LA)

   explicit Transform3D( const AxisAngle & r) {
      AssignFrom(Rotation3D(r));
   }
   explicit Transform3D( const EulerAngles & r) {
      AssignFrom(Rotation3D(r));
   }
   explicit Transform3D( const Quaternion & r) {
      AssignFrom(Rotation3D(r));
   }
   explicit Transform3D( const RotationZYX & r) {
      AssignFrom(Rotation3D(r));
   }

   // Constructors from axial rotations
   // TO DO: implement direct methods for axial rotations without going through Rotation3D
   explicit Transform3D( const RotationX & r) {
      AssignFrom(Rotation3D(r));
   }
   explicit Transform3D( const RotationY & r) {
      AssignFrom(Rotation3D(r));
   }
   explicit Transform3D( const RotationZ & r) {
      AssignFrom(Rotation3D(r));
   }

   /**
      Construct from a translation only, represented by any DisplacementVector3D
      and with an identity rotation
   */
   template<class CoordSystem, class Tag>
   explicit Transform3D( const DisplacementVector3D<CoordSystem,Tag> & v) {
      AssignFrom(Vector(v.X(),v.Y(),v.Z()));
   }
   /**
      Construct from a translation only, represented by a Cartesian 3D Vector,
      and with an identity rotation
   */
   explicit Transform3D( const Vector & v) {
      AssignFrom(v);
   }
   /**
      Construct from a translation only, represented by a Translation3D class
      and with an identity rotation
   */
   explicit Transform3D(const Translation3D<T> &t) { AssignFrom(t.Vect()); }

   //#if !defined(__MAKECINT__) && !defined(G__DICTIONARY)  // this is ambigous with double * , double *


#ifdef OLD_VERSION
   /**
      Construct from a translation (using any type of DisplacementVector )
      and then a rotation (any rotation object).
      Requirement on the rotation and vector objects are that they can be transformed in a
      Rotation3D class and in a Vector
   */
   template <class ARotation, class CoordSystem, class Tag>
   Transform3D(const DisplacementVector3D<CoordSystem,Tag> & v , const ARotation & r)
   {
      // is equivalent from having first the rotation and then the translation vector rotated
      Rotation3D r3d(r);
      AssignFrom( r3d, r3d( Vector(v.X(),v.Y(),v.Z()) ) );
   }
#endif

public:
   /**
      Construct transformation from one coordinate system defined by three
      points (origin + two axis) to
      a new coordinate system defined by other three points (origin + axis)
      Scalar version.
      @param fr0  point defining origin of original reference system
      @param fr1  point defining first axis of original reference system
      @param fr2  point defining second axis of original reference system
      @param to0  point defining origin of transformed reference system
      @param to1  point defining first axis transformed reference system
      @param to2  point defining second axis transformed reference system
   */
   template <typename SCALAR = T, typename std::enable_if<std::is_arithmetic<SCALAR>::value>::type * = nullptr>
   Transform3D(const Point &fr0, const Point &fr1, const Point &fr2, const Point &to0, const Point &to1,
               const Point &to2)
   {
      // takes impl. from CLHEP ( E.Chernyaev). To be checked

      Vector x1 = (fr1 - fr0).Unit();
      Vector y1 = (fr2 - fr0).Unit();
      Vector x2 = (to1 - to0).Unit();
      Vector y2 = (to2 - to0).Unit();

      //   C H E C K   A N G L E S

      const T cos1 = x1.Dot(y1);
      const T cos2 = x2.Dot(y2);

      if (std::fabs(T(1) - cos1) <= T(0.000001) || std::fabs(T(1) - cos2) <= T(0.000001)) {
         std::cerr << "Transform3D: Error : zero angle between axes" << std::endl;
         SetIdentity();
      } else {
         if (std::fabs(cos1 - cos2) > T(0.000001)) {
            std::cerr << "Transform3D: Warning: angles between axes are not equal" << std::endl;
         }

         //   F I N D   R O T A T I O N   M A T R I X

         Vector z1 = (x1.Cross(y1)).Unit();
         y1        = z1.Cross(x1);

         Vector z2 = (x2.Cross(y2)).Unit();
         y2        = z2.Cross(x2);

         T x1x = x1.x();
         T x1y = x1.y();
         T x1z = x1.z();
         T y1x = y1.x();
         T y1y = y1.y();
         T y1z = y1.z();
         T z1x = z1.x();
         T z1y = z1.y();
         T z1z = z1.z();

         T x2x = x2.x();
         T x2y = x2.y();
         T x2z = x2.z();
         T y2x = y2.x();
         T y2y = y2.y();
         T y2z = y2.z();
         T z2x = z2.x();
         T z2y = z2.y();
         T z2z = z2.z();

         T detxx = (y1y * z1z - z1y * y1z);
         T detxy = -(y1x * z1z - z1x * y1z);
         T detxz = (y1x * z1y - z1x * y1y);
         T detyx = -(x1y * z1z - z1y * x1z);
         T detyy = (x1x * z1z - z1x * x1z);
         T detyz = -(x1x * z1y - z1x * x1y);
         T detzx = (x1y * y1z - y1y * x1z);
         T detzy = -(x1x * y1z - y1x * x1z);
         T detzz = (x1x * y1y - y1x * x1y);

         T txx = x2x * detxx + y2x * detyx + z2x * detzx;
         T txy = x2x * detxy + y2x * detyy + z2x * detzy;
         T txz = x2x * detxz + y2x * detyz + z2x * detzz;
         T tyx = x2y * detxx + y2y * detyx + z2y * detzx;
         T tyy = x2y * detxy + y2y * detyy + z2y * detzy;
         T tyz = x2y * detxz + y2y * detyz + z2y * detzz;
         T tzx = x2z * detxx + y2z * detyx + z2z * detzx;
         T tzy = x2z * detxy + y2z * detyy + z2z * detzy;
         T tzz = x2z * detxz + y2z * detyz + z2z * detzz;

         //   S E T    T R A N S F O R M A T I O N

         T dx1 = fr0.x(), dy1 = fr0.y(), dz1 = fr0.z();
         T dx2 = to0.x(), dy2 = to0.y(), dz2 = to0.z();

         SetComponents(txx, txy, txz, dx2 - txx * dx1 - txy * dy1 - txz * dz1, tyx, tyy, tyz,
                       dy2 - tyx * dx1 - tyy * dy1 - tyz * dz1, tzx, tzy, tzz, dz2 - tzx * dx1 - tzy * dy1 - tzz * dz1);
      }
   }

   /**
      Construct transformation from one coordinate system defined by three
      points (origin + two axis) to
      a new coordinate system defined by other three points (origin + axis)
      Vectorised version.
      @param fr0  point defining origin of original reference system
      @param fr1  point defining first axis of original reference system
      @param fr2  point defining second axis of original reference system
      @param to0  point defining origin of transformed reference system
      @param to1  point defining first axis transformed reference system
      @param to2  point defining second axis transformed reference system
   */
   template <typename SCALAR = T, typename std::enable_if<!std::is_arithmetic<SCALAR>::value>::type * = nullptr>
   Transform3D(const Point &fr0, const Point &fr1, const Point &fr2, const Point &to0, const Point &to1,
               const Point &to2)
   {
      // takes impl. from CLHEP ( E.Chernyaev). To be checked

      Vector x1 = (fr1 - fr0).Unit();
      Vector y1 = (fr2 - fr0).Unit();
      Vector x2 = (to1 - to0).Unit();
      Vector y2 = (to2 - to0).Unit();

      //   C H E C K   A N G L E S

      const T cos1 = x1.Dot(y1);
      const T cos2 = x2.Dot(y2);

      const auto m1 = (abs(T(1) - cos1) <= T(0.000001) || abs(T(1) - cos2) <= T(0.000001));

      const auto m2 = (abs(cos1 - cos2) > T(0.000001));
      if (any_of(m2)) {
         std::cerr << "Transform3D: Warning: angles between axes are not equal" << std::endl;
      }

      //   F I N D   R O T A T I O N   M A T R I X

      Vector z1 = (x1.Cross(y1)).Unit();
      y1        = z1.Cross(x1);

      Vector z2 = (x2.Cross(y2)).Unit();
      y2        = z2.Cross(x2);

      T x1x = x1.x();
      T x1y = x1.y();
      T x1z = x1.z();
      T y1x = y1.x();
      T y1y = y1.y();
      T y1z = y1.z();
      T z1x = z1.x();
      T z1y = z1.y();
      T z1z = z1.z();

      T x2x = x2.x();
      T x2y = x2.y();
      T x2z = x2.z();
      T y2x = y2.x();
      T y2y = y2.y();
      T y2z = y2.z();
      T z2x = z2.x();
      T z2y = z2.y();
      T z2z = z2.z();

      T detxx = (y1y * z1z - z1y * y1z);
      T detxy = -(y1x * z1z - z1x * y1z);
      T detxz = (y1x * z1y - z1x * y1y);
      T detyx = -(x1y * z1z - z1y * x1z);
      T detyy = (x1x * z1z - z1x * x1z);
      T detyz = -(x1x * z1y - z1x * x1y);
      T detzx = (x1y * y1z - y1y * x1z);
      T detzy = -(x1x * y1z - y1x * x1z);
      T detzz = (x1x * y1y - y1x * x1y);

      T txx = x2x * detxx + y2x * detyx + z2x * detzx;
      T txy = x2x * detxy + y2x * detyy + z2x * detzy;
      T txz = x2x * detxz + y2x * detyz + z2x * detzz;
      T tyx = x2y * detxx + y2y * detyx + z2y * detzx;
      T tyy = x2y * detxy + y2y * detyy + z2y * detzy;
      T tyz = x2y * detxz + y2y * detyz + z2y * detzz;
      T tzx = x2z * detxx + y2z * detyx + z2z * detzx;
      T tzy = x2z * detxy + y2z * detyy + z2z * detzy;
      T tzz = x2z * detxz + y2z * detyz + z2z * detzz;

      //   S E T    T R A N S F O R M A T I O N

      T dx1 = fr0.x(), dy1 = fr0.y(), dz1 = fr0.z();
      T dx2 = to0.x(), dy2 = to0.y(), dz2 = to0.z();

      SetComponents(txx, txy, txz, dx2 - txx * dx1 - txy * dy1 - txz * dz1, tyx, tyy, tyz,
                    dy2 - tyx * dx1 - tyy * dy1 - tyz * dz1, tzx, tzy, tzz, dz2 - tzx * dx1 - tzy * dy1 - tzz * dz1);

      if (any_of(m1)) {
         std::cerr << "Transform3D: Error : zero angle between axes" << std::endl;
         SetIdentity(m1);
      }
   }

   // use compiler generated copy ctor, copy assignmet and dtor

   /**
      Construct from a linear algebra matrix of size at least 3x4,
      which must support operator()(i,j) to obtain elements (0,0) thru (2,3).
      The 3x3 sub-block is assumed to be the rotation part and the translations vector
      are described by the 4-th column
   */
   template<class ForeignMatrix>
   explicit Transform3D(const ForeignMatrix & m) {
      SetComponents(m);
   }

   /**
      Raw constructor from 12 Scalar components
   */
   Transform3D(T xx, T xy, T xz, T dx, T yx, T yy, T yz, T dy, T zx, T zy, T zz, T dz)
   {
      SetComponents (xx, xy, xz, dx, yx, yy, yz, dy, zx, zy, zz, dz);
   }


   /**
      Construct from a linear algebra matrix of size at least 3x4,
      which must support operator()(i,j) to obtain elements (0,0) thru (2,3).
      The 3x3 sub-block is assumed to be the rotation part and the translations vector
      are described by the 4-th column
   */
   template <class ForeignMatrix>
   Transform3D<T> &operator=(const ForeignMatrix &m)
   {
      SetComponents(m);
      return *this;
   }


   // ======== Components ==============


   /**
      Set the 12 matrix components given an iterator to the start of
      the desired data, and another to the end (12 past start).
   */
   template<class IT>
   void SetComponents(IT begin, IT end) {
      for (int i = 0; i <12; ++i) {
         fM[i] = *begin;
         ++begin;
      }
      (void)end;
      assert (end==begin);
   }

   /**
      Get the 12 matrix components into data specified by an iterator begin
      and another to the end of the desired data (12 past start).
   */
   template<class IT>
   void GetComponents(IT begin, IT end) const {
      for (int i = 0; i <12; ++i) {
         *begin = fM[i];
         ++begin;
      }
      (void)end;
      assert (end==begin);
   }

   /**
      Get the 12 matrix components into data specified by an iterator begin
   */
   template<class IT>
   void GetComponents(IT begin) const {
      std::copy(fM, fM + 12, begin);
   }

   /**
      Set components from a linear algebra matrix of size at least 3x4,
      which must support operator()(i,j) to obtain elements (0,0) thru (2,3).
      The 3x3 sub-block is assumed to be the rotation part and the translations vector
      are described by the 4-th column
   */
   template<class ForeignMatrix>
   void
   SetTransformMatrix (const ForeignMatrix & m) {
      fM[kXX]=m(0,0);  fM[kXY]=m(0,1);  fM[kXZ]=m(0,2); fM[kDX]=m(0,3);
      fM[kYX]=m(1,0);  fM[kYY]=m(1,1);  fM[kYZ]=m(1,2); fM[kDY]=m(1,3);
      fM[kZX]=m(2,0);  fM[kZY]=m(2,1);  fM[kZZ]=m(2,2); fM[kDZ]=m(2,3);
   }

   /**
      Get components into a linear algebra matrix of size at least 3x4,
      which must support operator()(i,j) for write access to elements
      (0,0) thru (2,3).
   */
   template<class ForeignMatrix>
   void
   GetTransformMatrix (ForeignMatrix & m) const {
      m(0,0)=fM[kXX];  m(0,1)=fM[kXY];  m(0,2)=fM[kXZ];  m(0,3)=fM[kDX];
      m(1,0)=fM[kYX];  m(1,1)=fM[kYY];  m(1,2)=fM[kYZ];  m(1,3)=fM[kDY];
      m(2,0)=fM[kZX];  m(2,1)=fM[kZY];  m(2,2)=fM[kZZ];  m(2,3)=fM[kDZ];
   }


   /**
      Set the components from 12 scalars
   */
   void SetComponents(T xx, T xy, T xz, T dx, T yx, T yy, T yz, T dy, T zx, T zy, T zz, T dz)
   {
      fM[kXX]=xx;  fM[kXY]=xy;  fM[kXZ]=xz;  fM[kDX]=dx;
      fM[kYX]=yx;  fM[kYY]=yy;  fM[kYZ]=yz;  fM[kDY]=dy;
      fM[kZX]=zx;  fM[kZY]=zy;  fM[kZZ]=zz;  fM[kDZ]=dz;
   }

   /**
      Get the components into 12 scalars
   */
   void GetComponents(T &xx, T &xy, T &xz, T &dx, T &yx, T &yy, T &yz, T &dy, T &zx, T &zy, T &zz, T &dz) const
   {
      xx=fM[kXX];  xy=fM[kXY];  xz=fM[kXZ];  dx=fM[kDX];
      yx=fM[kYX];  yy=fM[kYY];  yz=fM[kYZ];  dy=fM[kDY];
      zx=fM[kZX];  zy=fM[kZY];  zz=fM[kZZ];  dz=fM[kDZ];
   }


   /**
      Get the rotation and translation vector representing the 3D transformation
      in any rotation and any vector (the Translation class could also be used)
   */
   template<class AnyRotation, class V>
   void GetDecomposition(AnyRotation &r, V &v) const {
      GetRotation(r);
      GetTranslation(v);
   }


   /**
      Get the rotation and translation vector representing the 3D transformation
   */
   void GetDecomposition(Rotation3D &r, Vector &v) const {
      GetRotation(r);
      GetTranslation(v);
   }

   /**
      Get the 3D rotation representing the 3D transformation
   */
   Rotation3D Rotation() const {
      return Rotation3D( fM[kXX], fM[kXY], fM[kXZ],
                         fM[kYX], fM[kYY], fM[kYZ],
                         fM[kZX], fM[kZY], fM[kZZ] );
   }

   /**
      Get the rotation representing the 3D transformation
   */
   template <class AnyRotation>
   AnyRotation Rotation() const {
      return AnyRotation(Rotation3D(fM[kXX], fM[kXY], fM[kXZ], fM[kYX], fM[kYY], fM[kYZ], fM[kZX], fM[kZY], fM[kZZ]));
   }

   /**
      Get the  rotation (any type) representing the 3D transformation
   */
   template <class AnyRotation>
   void GetRotation(AnyRotation &r) const {
      r = Rotation();
   }

   /**
      Get the translation representing the 3D transformation in a Cartesian vector
   */
   Translation3D<T> Translation() const { return Translation3D<T>(fM[kDX], fM[kDY], fM[kDZ]); }

   /**
      Get the translation representing the 3D transformation in any vector
      which implements the SetXYZ method
   */
   template <class AnyVector>
   void GetTranslation(AnyVector &v) const {
      v.SetXYZ(fM[kDX], fM[kDY], fM[kDZ]);
   }



   // operations on points and vectors

   /**
      Transformation operation for Position Vector in Cartesian coordinate
      For a Position Vector first a rotation and then a translation is applied
   */
   Point operator() (const Point & p) const {
      return Point ( fM[kXX]*p.X() + fM[kXY]*p.Y() + fM[kXZ]*p.Z() + fM[kDX],
                     fM[kYX]*p.X() + fM[kYY]*p.Y() + fM[kYZ]*p.Z() + fM[kDY],
                     fM[kZX]*p.X() + fM[kZY]*p.Y() + fM[kZZ]*p.Z() + fM[kDZ] );
   }


   /**
      Transformation operation for Displacement Vectors in Cartesian coordinate
      For the Displacement Vectors only the rotation applies - no translations
   */
   Vector operator() (const Vector & v) const {
      return Vector( fM[kXX]*v.X() + fM[kXY]*v.Y() + fM[kXZ]*v.Z() ,
                     fM[kYX]*v.X() + fM[kYY]*v.Y() + fM[kYZ]*v.Z() ,
                     fM[kZX]*v.X() + fM[kZY]*v.Y() + fM[kZZ]*v.Z()  );
   }


   /**
      Transformation operation for Position Vector in any coordinate system
   */
   template <class CoordSystem>
   PositionVector3D<CoordSystem> operator()(const PositionVector3D<CoordSystem> &p) const
   {
      return PositionVector3D<CoordSystem>(operator()(Point(p)));
   }
   /**
      Transformation operation for Position Vector in any coordinate system
   */
   template <class CoordSystem>
   PositionVector3D<CoordSystem> operator*(const PositionVector3D<CoordSystem> &v) const
   {
      return operator()(v);
   }

   /**
      Transformation operation for Displacement Vector in any coordinate system
   */
   template<class CoordSystem >
   DisplacementVector3D<CoordSystem> operator() (const DisplacementVector3D <CoordSystem> & v) const {
      return DisplacementVector3D<CoordSystem>(operator()(Vector(v)));
   }
   /**
      Transformation operation for Displacement Vector in any coordinate system
   */
   template <class CoordSystem>
   DisplacementVector3D<CoordSystem> operator*(const DisplacementVector3D<CoordSystem> &v) const
   {
      return operator()(v);
   }

   /**
      Directly apply the inverse affine transformation on vectors.
      Avoids having to calculate the inverse as an intermediate result.
      This is possible since the inverse of a rotation is its transpose.
   */
   Vector ApplyInverse(const Vector &v) const
   {
      return Vector(fM[kXX] * v.X() + fM[kYX] * v.Y() + fM[kZX] * v.Z(),
                    fM[kXY] * v.X() + fM[kYY] * v.Y() + fM[kZY] * v.Z(),
                    fM[kXZ] * v.X() + fM[kYZ] * v.Y() + fM[kZZ] * v.Z());
   }

   /**
      Directly apply the inverse affine transformation on points
      (first inverse translation then inverse rotation).
      Avoids having to calculate the inverse as an intermediate result.
      This is possible since the inverse of a rotation is its transpose.
   */
   Point ApplyInverse(const Point &p) const
   {
      Point tmp(p.X() - fM[kDX], p.Y() - fM[kDY], p.Z() - fM[kDZ]);
      return Point(fM[kXX] * tmp.X() + fM[kYX] * tmp.Y() + fM[kZX] * tmp.Z(),
                   fM[kXY] * tmp.X() + fM[kYY] * tmp.Y() + fM[kZY] * tmp.Z(),
                   fM[kXZ] * tmp.X() + fM[kYZ] * tmp.Y() + fM[kZZ] * tmp.Z());
   }

   /**
      Directly apply the inverse affine transformation on an arbitrary
      coordinate-system point.
      Involves casting to Point(p) type.
   */
   template <class CoordSystem>
   PositionVector3D<CoordSystem> ApplyInverse(const PositionVector3D<CoordSystem> &p) const
   {
      return PositionVector3D<CoordSystem>(ApplyInverse(Point(p)));
   }

   /**
      Directly apply the inverse affine transformation on an arbitrary
      coordinate-system vector.
      Involves casting to Vector(p) type.
   */
   template <class CoordSystem>
   DisplacementVector3D<CoordSystem> ApplyInverse(const DisplacementVector3D<CoordSystem> &p) const
   {
      return DisplacementVector3D<CoordSystem>(ApplyInverse(Vector(p)));
   }

   /**
      Transformation operation for points between different coordinate system tags
   */
   template <class CoordSystem, class Tag1, class Tag2>
   void Transform(const PositionVector3D<CoordSystem, Tag1> &p1, PositionVector3D<CoordSystem, Tag2> &p2) const
   {
      const Point xyzNew = operator()(Point(p1.X(), p1.Y(), p1.Z()));
      p2.SetXYZ( xyzNew.X(), xyzNew.Y(), xyzNew.Z() );
   }


   /**
      Transformation operation for Displacement Vector of different coordinate systems
   */
   template <class CoordSystem, class Tag1, class Tag2>
   void Transform(const DisplacementVector3D<CoordSystem, Tag1> &v1, DisplacementVector3D<CoordSystem, Tag2> &v2) const
   {
      const Vector xyzNew = operator()(Vector(v1.X(), v1.Y(), v1.Z()));
      v2.SetXYZ( xyzNew.X(), xyzNew.Y(), xyzNew.Z() );
   }

   /**
      Transformation operation for a Lorentz Vector in any  coordinate system
   */
   template <class CoordSystem >
   LorentzVector<CoordSystem> operator() (const LorentzVector<CoordSystem> & q) const {
      const Vector xyzNew = operator()(Vector(q.Vect()));
      return LorentzVector<CoordSystem>(xyzNew.X(), xyzNew.Y(), xyzNew.Z(), q.E());
   }
   /**
      Transformation operation for a Lorentz Vector in any  coordinate system
   */
   template <class CoordSystem>
   LorentzVector<CoordSystem> operator*(const LorentzVector<CoordSystem> &q) const
   {
      return operator()(q);
   }

   /**
      Transformation on a 3D plane
   */
   template <typename TYPE>
   Plane3D<TYPE> operator()(const Plane3D<TYPE> &plane) const
   {
      // transformations on a 3D plane
      const auto n = plane.Normal();
      // take a point on the plane. Use origin projection on the plane
      // ( -ad, -bd, -cd) if (a**2 + b**2 + c**2 ) = 1
      const auto d = plane.HesseDistance();
      Point p(-d * n.X(), -d * n.Y(), -d * n.Z());
      return Plane3D<TYPE>(operator()(n), operator()(p));
   }

   /// Multiplication operator for 3D plane
   template <typename TYPE>
   Plane3D<TYPE> operator*(const Plane3D<TYPE> &plane) const
   {
      return operator()(plane);
   }

   // skip transformation for arbitrary vectors - not really defined if point or displacement vectors

   /**
      multiply (combine) with another transformation in place
   */
   inline Transform3D<T> &operator*=(const Transform3D<T> &t);

   /**
      multiply (combine) two transformations
   */
   inline Transform3D<T> operator*(const Transform3D<T> &t) const;

   /**
       Invert the transformation in place (scalar)
   */
   template <typename SCALAR = T, typename std::enable_if<std::is_arithmetic<SCALAR>::value>::type * = nullptr>
   void Invert()
   {
      //
      // Name: Transform3D::inverse                     Date:    24.09.96
      // Author: E.Chernyaev (IHEP/Protvino)            Revised:
      //
      // Function: Find inverse affine transformation.

      T detxx = fM[kYY] * fM[kZZ] - fM[kYZ] * fM[kZY];
      T detxy = fM[kYX] * fM[kZZ] - fM[kYZ] * fM[kZX];
      T detxz = fM[kYX] * fM[kZY] - fM[kYY] * fM[kZX];
      T det   = fM[kXX] * detxx - fM[kXY] * detxy + fM[kXZ] * detxz;
      if (det == T(0)) {
         std::cerr << "Transform3D::inverse error: zero determinant" << std::endl;
         return;
      }
      det = T(1) / det;
      detxx *= det;
      detxy *= det;
      detxz *= det;
      T detyx = (fM[kXY] * fM[kZZ] - fM[kXZ] * fM[kZY]) * det;
      T detyy = (fM[kXX] * fM[kZZ] - fM[kXZ] * fM[kZX]) * det;
      T detyz = (fM[kXX] * fM[kZY] - fM[kXY] * fM[kZX]) * det;
      T detzx = (fM[kXY] * fM[kYZ] - fM[kXZ] * fM[kYY]) * det;
      T detzy = (fM[kXX] * fM[kYZ] - fM[kXZ] * fM[kYX]) * det;
      T detzz = (fM[kXX] * fM[kYY] - fM[kXY] * fM[kYX]) * det;
      SetComponents(detxx, -detyx, detzx, -detxx * fM[kDX] + detyx * fM[kDY] - detzx * fM[kDZ], -detxy, detyy, -detzy,
                    detxy * fM[kDX] - detyy * fM[kDY] + detzy * fM[kDZ], detxz, -detyz, detzz,
                    -detxz * fM[kDX] + detyz * fM[kDY] - detzz * fM[kDZ]);
   }

   /**
       Invert the transformation in place (vectorised)
   */
   template <typename SCALAR = T, typename std::enable_if<!std::is_arithmetic<SCALAR>::value>::type * = nullptr>
   void Invert()
   {
      //
      // Name: Transform3D::inverse                     Date:    24.09.96
      // Author: E.Chernyaev (IHEP/Protvino)            Revised:
      //
      // Function: Find inverse affine transformation.

      T          detxx    = fM[kYY] * fM[kZZ] - fM[kYZ] * fM[kZY];
      T          detxy    = fM[kYX] * fM[kZZ] - fM[kYZ] * fM[kZX];
      T          detxz    = fM[kYX] * fM[kZY] - fM[kYY] * fM[kZX];
      T          det      = fM[kXX] * detxx - fM[kXY] * detxy + fM[kXZ] * detxz;
      const auto detZmask = (det == T(0));
      if (any_of(detZmask)) {
         std::cerr << "Transform3D::inverse error: zero determinant" << std::endl;
         det(detZmask) = T(1);
      }
      det = T(1) / det;
      detxx *= det;
      detxy *= det;
      detxz *= det;
      T detyx = (fM[kXY] * fM[kZZ] - fM[kXZ] * fM[kZY]) * det;
      T detyy = (fM[kXX] * fM[kZZ] - fM[kXZ] * fM[kZX]) * det;
      T detyz = (fM[kXX] * fM[kZY] - fM[kXY] * fM[kZX]) * det;
      T detzx = (fM[kXY] * fM[kYZ] - fM[kXZ] * fM[kYY]) * det;
      T detzy = (fM[kXX] * fM[kYZ] - fM[kXZ] * fM[kYX]) * det;
      T detzz = (fM[kXX] * fM[kYY] - fM[kXY] * fM[kYX]) * det;
      // Set det=0 cases to 0
      if (any_of(detZmask)) {
         detxx(detZmask) = T(0);
         detxy(detZmask) = T(0);
         detxz(detZmask) = T(0);
         detyx(detZmask) = T(0);
         detyy(detZmask) = T(0);
         detyz(detZmask) = T(0);
         detzx(detZmask) = T(0);
         detzy(detZmask) = T(0);
         detzz(detZmask) = T(0);
      }
      // set final components
      SetComponents(detxx, -detyx, detzx, -detxx * fM[kDX] + detyx * fM[kDY] - detzx * fM[kDZ], -detxy, detyy, -detzy,
                    detxy * fM[kDX] - detyy * fM[kDY] + detzy * fM[kDZ], detxz, -detyz, detzz,
                    -detxz * fM[kDX] + detyz * fM[kDY] - detzz * fM[kDZ]);
   }

   /**
      Return the inverse of the transformation.
   */
   Transform3D<T> Inverse() const
   {
      Transform3D<T> t(*this);
      t.Invert();
      return t;
   }

   /**
      Equality operator. Check equality for each element
      To do: use T tolerance
   */
   bool operator==(const Transform3D<T> &rhs) const
   {
      return (fM[0] == rhs.fM[0] && fM[1] == rhs.fM[1] && fM[2] == rhs.fM[2] && fM[3] == rhs.fM[3] &&
              fM[4] == rhs.fM[4] && fM[5] == rhs.fM[5] && fM[6] == rhs.fM[6] && fM[7] == rhs.fM[7] &&
              fM[8] == rhs.fM[8] && fM[9] == rhs.fM[9] && fM[10] == rhs.fM[10] && fM[11] == rhs.fM[11]);
   }

   /**
      Inequality operator. Check equality for each element
      To do: use T tolerance
   */
   bool operator!=(const Transform3D<T> &rhs) const { return !operator==(rhs); }

protected:

   /**
      make transformation from first a rotation then a translation
   */
   void AssignFrom(const Rotation3D &r, const Vector &v)
   {
      // assignment  from rotation + translation

      T rotData[9];
      r.GetComponents(rotData, rotData + 9);
      // first raw
      for (int i = 0; i < 3; ++i) fM[i] = rotData[i];
      // second raw
      for (int i = 0; i < 3; ++i) fM[kYX + i] = rotData[3 + i];
      // third raw
      for (int i = 0; i < 3; ++i) fM[kZX + i] = rotData[6 + i];

      // translation data
      T vecData[3];
      v.GetCoordinates(vecData, vecData + 3);
      fM[kDX] = vecData[0];
      fM[kDY] = vecData[1];
      fM[kDZ] = vecData[2];
   }

   /**
      make transformation from only rotations (zero translation)
   */
   void AssignFrom(const Rotation3D &r)
   {
      // assign from only a rotation  (null translation)
      T rotData[9];
      r.GetComponents(rotData, rotData + 9);
      for (int i = 0; i < 3; ++i) {
         for (int j = 0; j < 3; ++j) fM[4 * i + j] = rotData[3 * i + j];
         // empty vector data
         fM[4 * i + 3] = T(0);
      }
   }

   /**
      make transformation from only translation (identity rotations)
   */
   void AssignFrom(const Vector &v)
   {
      // assign from a translation only (identity rotations)
      fM[kXX] = T(1);
      fM[kXY] = T(0);
      fM[kXZ] = T(0);
      fM[kDX] = v.X();
      fM[kYX] = T(0);
      fM[kYY] = T(1);
      fM[kYZ] = T(0);
      fM[kDY] = v.Y();
      fM[kZX] = T(0);
      fM[kZY] = T(0);
      fM[kZZ] = T(1);
      fM[kDZ] = v.Z();
   }

   /**
      Set identity transformation (identity rotation , zero translation)
   */
   void SetIdentity()
   {
      // set identity ( identity rotation and zero translation)
      fM[kXX] = T(1);
      fM[kXY] = T(0);
      fM[kXZ] = T(0);
      fM[kDX] = T(0);
      fM[kYX] = T(0);
      fM[kYY] = T(1);
      fM[kYZ] = T(0);
      fM[kDY] = T(0);
      fM[kZX] = T(0);
      fM[kZY] = T(0);
      fM[kZZ] = T(1);
      fM[kDZ] = T(0);
   }

   /**
      Set identity transformation (identity rotation , zero translation)
      vectorised version that sets using a mask
   */
   template <typename SCALAR = T, typename std::enable_if<!std::is_arithmetic<SCALAR>::value>::type * = nullptr>
   void SetIdentity(const typename SCALAR::mask_type m)
   {
      // set identity ( identity rotation and zero translation)
      fM[kXX](m) = T(1);
      fM[kXY](m) = T(0);
      fM[kXZ](m) = T(0);
      fM[kDX](m) = T(0);
      fM[kYX](m) = T(0);
      fM[kYY](m) = T(1);
      fM[kYZ](m) = T(0);
      fM[kDY](m) = T(0);
      fM[kZX](m) = T(0);
      fM[kZY](m) = T(0);
      fM[kZZ](m) = T(1);
      fM[kDZ](m) = T(0);
   }

private:
   T fM[12]; // transformation elements (3x4 matrix)
};




// inline functions (combination of transformations)

template <class T>
inline Transform3D<T> &Transform3D<T>::operator*=(const Transform3D<T> &t)
{
   // combination of transformations

   SetComponents(fM[kXX]*t.fM[kXX]+fM[kXY]*t.fM[kYX]+fM[kXZ]*t.fM[kZX],
                 fM[kXX]*t.fM[kXY]+fM[kXY]*t.fM[kYY]+fM[kXZ]*t.fM[kZY],
                 fM[kXX]*t.fM[kXZ]+fM[kXY]*t.fM[kYZ]+fM[kXZ]*t.fM[kZZ],
                 fM[kXX]*t.fM[kDX]+fM[kXY]*t.fM[kDY]+fM[kXZ]*t.fM[kDZ]+fM[kDX],

                 fM[kYX]*t.fM[kXX]+fM[kYY]*t.fM[kYX]+fM[kYZ]*t.fM[kZX],
                 fM[kYX]*t.fM[kXY]+fM[kYY]*t.fM[kYY]+fM[kYZ]*t.fM[kZY],
                 fM[kYX]*t.fM[kXZ]+fM[kYY]*t.fM[kYZ]+fM[kYZ]*t.fM[kZZ],
                 fM[kYX]*t.fM[kDX]+fM[kYY]*t.fM[kDY]+fM[kYZ]*t.fM[kDZ]+fM[kDY],

                 fM[kZX]*t.fM[kXX]+fM[kZY]*t.fM[kYX]+fM[kZZ]*t.fM[kZX],
                 fM[kZX]*t.fM[kXY]+fM[kZY]*t.fM[kYY]+fM[kZZ]*t.fM[kZY],
                 fM[kZX]*t.fM[kXZ]+fM[kZY]*t.fM[kYZ]+fM[kZZ]*t.fM[kZZ],
                 fM[kZX]*t.fM[kDX]+fM[kZY]*t.fM[kDY]+fM[kZZ]*t.fM[kDZ]+fM[kDZ]);

   return *this;
}

template <class T>
inline Transform3D<T> Transform3D<T>::operator*(const Transform3D<T> &t) const
{
   // combination of transformations

   return Transform3D<T>(fM[kXX] * t.fM[kXX] + fM[kXY] * t.fM[kYX] + fM[kXZ] * t.fM[kZX],
                         fM[kXX] * t.fM[kXY] + fM[kXY] * t.fM[kYY] + fM[kXZ] * t.fM[kZY],
                         fM[kXX] * t.fM[kXZ] + fM[kXY] * t.fM[kYZ] + fM[kXZ] * t.fM[kZZ],
                         fM[kXX] * t.fM[kDX] + fM[kXY] * t.fM[kDY] + fM[kXZ] * t.fM[kDZ] + fM[kDX],

                         fM[kYX] * t.fM[kXX] + fM[kYY] * t.fM[kYX] + fM[kYZ] * t.fM[kZX],
                         fM[kYX] * t.fM[kXY] + fM[kYY] * t.fM[kYY] + fM[kYZ] * t.fM[kZY],
                         fM[kYX] * t.fM[kXZ] + fM[kYY] * t.fM[kYZ] + fM[kYZ] * t.fM[kZZ],
                         fM[kYX] * t.fM[kDX] + fM[kYY] * t.fM[kDY] + fM[kYZ] * t.fM[kDZ] + fM[kDY],

                         fM[kZX] * t.fM[kXX] + fM[kZY] * t.fM[kYX] + fM[kZZ] * t.fM[kZX],
                         fM[kZX] * t.fM[kXY] + fM[kZY] * t.fM[kYY] + fM[kZZ] * t.fM[kZY],
                         fM[kZX] * t.fM[kXZ] + fM[kZY] * t.fM[kYZ] + fM[kZZ] * t.fM[kZZ],
                         fM[kZX] * t.fM[kDX] + fM[kZY] * t.fM[kDY] + fM[kZZ] * t.fM[kDZ] + fM[kDZ]);
}




//--- global functions resulting in Transform3D -------


// ------ combination of a  translation (first)  and a rotation ------


/**
   combine a translation and a rotation to give a transform3d
   First the translation then the rotation
 */
template <class T>
inline Transform3D<T> operator*(const Rotation3D &r, const Translation3D<T> &t)
{
   return Transform3D<T>(r, r(t.Vect()));
}
template <class T>
inline Transform3D<T> operator*(const RotationX &r, const Translation3D<T> &t)
{
   Rotation3D r3(r);
   return Transform3D<T>(r3, r3(t.Vect()));
}
template <class T>
inline Transform3D<T> operator*(const RotationY &r, const Translation3D<T> &t)
{
   Rotation3D r3(r);
   return Transform3D<T>(r3, r3(t.Vect()));
}
template <class T>
inline Transform3D<T> operator*(const RotationZ &r, const Translation3D<T> &t)
{
   Rotation3D r3(r);
   return Transform3D<T>(r3, r3(t.Vect()));
}
template <class T>
inline Transform3D<T> operator*(const RotationZYX &r, const Translation3D<T> &t)
{
   Rotation3D r3(r);
   return Transform3D<T>(r3, r3(t.Vect()));
}
template <class T>
inline Transform3D<T> operator*(const AxisAngle &r, const Translation3D<T> &t)
{
   Rotation3D r3(r);
   return Transform3D<T>(r3, r3(t.Vect()));
}
template <class T>
inline Transform3D<T> operator*(const EulerAngles &r, const Translation3D<T> &t)
{
   Rotation3D r3(r);
   return Transform3D<T>(r3, r3(t.Vect()));
}
template <class T>
inline Transform3D<T> operator*(const Quaternion &r, const Translation3D<T> &t)
{
   Rotation3D r3(r);
   return Transform3D<T>(r3, r3(t.Vect()));
}

// ------ combination of a  rotation (first)  and then a translation ------

/**
   combine a rotation and a translation to give a transform3d
   First a rotation then the translation
 */
template <class T>
inline Transform3D<T> operator*(const Translation3D<T> &t, const Rotation3D &r)
{
   return Transform3D<T>(r, t.Vect());
}
template <class T>
inline Transform3D<T> operator*(const Translation3D<T> &t, const RotationX &r)
{
   return Transform3D<T>(Rotation3D(r), t.Vect());
}
template <class T>
inline Transform3D<T> operator*(const Translation3D<T> &t, const RotationY &r)
{
   return Transform3D<T>(Rotation3D(r), t.Vect());
}
template <class T>
inline Transform3D<T> operator*(const Translation3D<T> &t, const RotationZ &r)
{
   return Transform3D<T>(Rotation3D(r), t.Vect());
}
template <class T>
inline Transform3D<T> operator*(const Translation3D<T> &t, const RotationZYX &r)
{
   return Transform3D<T>(Rotation3D(r), t.Vect());
}
template <class T>
inline Transform3D<T> operator*(const Translation3D<T> &t, const EulerAngles &r)
{
   return Transform3D<T>(Rotation3D(r), t.Vect());
}
template <class T>
inline Transform3D<T> operator*(const Translation3D<T> &t, const Quaternion &r)
{
   return Transform3D<T>(Rotation3D(r), t.Vect());
}
template <class T>
inline Transform3D<T> operator*(const Translation3D<T> &t, const AxisAngle &r)
{
   return Transform3D<T>(Rotation3D(r), t.Vect());
}

// ------ combination of a Transform3D and a pure translation------

/**
   combine a transformation and a translation to give a transform3d
   First the translation then the transform3D
 */
template <class T>
inline Transform3D<T> operator*(const Transform3D<T> &t, const Translation3D<T> &d)
{
   Rotation3D r = t.Rotation();
   return Transform3D<T>(r, r(d.Vect()) + t.Translation().Vect());
}

/**
   combine a translation and a transformation to give a transform3d
   First the transformation then the translation
 */
template <class T>
inline Transform3D<T> operator*(const Translation3D<T> &d, const Transform3D<T> &t)
{
   return Transform3D<T>(t.Rotation(), t.Translation().Vect() + d.Vect());
}

// ------ combination of a Transform3D and any rotation------


/**
   combine a transformation and a rotation to give a transform3d
   First the rotation then the transform3D
 */
template <class T>
inline Transform3D<T> operator*(const Transform3D<T> &t, const Rotation3D &r)
{
   return Transform3D<T>(t.Rotation() * r, t.Translation());
}
template <class T>
inline Transform3D<T> operator*(const Transform3D<T> &t, const RotationX &r)
{
   return Transform3D<T>(t.Rotation() * r, t.Translation());
}
template <class T>
inline Transform3D<T> operator*(const Transform3D<T> &t, const RotationY &r)
{
   return Transform3D<T>(t.Rotation() * r, t.Translation());
}
template <class T>
inline Transform3D<T> operator*(const Transform3D<T> &t, const RotationZ &r)
{
   return Transform3D<T>(t.Rotation() * r, t.Translation());
}
template <class T>
inline Transform3D<T> operator*(const Transform3D<T> &t, const RotationZYX &r)
{
   return Transform3D<T>(t.Rotation() * r, t.Translation());
}
template <class T>
inline Transform3D<T> operator*(const Transform3D<T> &t, const EulerAngles &r)
{
   return Transform3D<T>(t.Rotation() * r, t.Translation());
}
template <class T>
inline Transform3D<T> operator*(const Transform3D<T> &t, const AxisAngle &r)
{
   return Transform3D<T>(t.Rotation() * r, t.Translation());
}
template <class T>
inline Transform3D<T> operator*(const Transform3D<T> &t, const Quaternion &r)
{
   return Transform3D<T>(t.Rotation() * r, t.Translation());
}



/**
   combine a rotation and a transformation to give a transform3d
   First the transformation then the rotation
 */
template <class T>
inline Transform3D<T> operator*(const Rotation3D &r, const Transform3D<T> &t)
{
   return Transform3D<T>(r * t.Rotation(), r * t.Translation().Vect());
}
template <class T>
inline Transform3D<T> operator*(const RotationX &r, const Transform3D<T> &t)
{
   Rotation3D r3d(r);
   return Transform3D<T>(r3d * t.Rotation(), r3d * t.Translation().Vect());
}
template <class T>
inline Transform3D<T> operator*(const RotationY &r, const Transform3D<T> &t)
{
   Rotation3D r3d(r);
   return Transform3D<T>(r3d * t.Rotation(), r3d * t.Translation().Vect());
}
template <class T>
inline Transform3D<T> operator*(const RotationZ &r, const Transform3D<T> &t)
{
   Rotation3D r3d(r);
   return Transform3D<T>(r3d * t.Rotation(), r3d * t.Translation().Vect());
}
template <class T>
inline Transform3D<T> operator*(const RotationZYX &r, const Transform3D<T> &t)
{
   Rotation3D r3d(r);
   return Transform3D<T>(r3d * t.Rotation(), r3d * t.Translation().Vect());
}
template <class T>
inline Transform3D<T> operator*(const EulerAngles &r, const Transform3D<T> &t)
{
   Rotation3D r3d(r);
   return Transform3D<T>(r3d * t.Rotation(), r3d * t.Translation().Vect());
}
template <class T>
inline Transform3D<T> operator*(const AxisAngle &r, const Transform3D<T> &t)
{
   Rotation3D r3d(r);
   return Transform3D<T>(r3d * t.Rotation(), r3d * t.Translation().Vect());
}
template <class T>
inline Transform3D<T> operator*(const Quaternion &r, const Transform3D<T> &t)
{
   Rotation3D r3d(r);
   return Transform3D<T>(r3d * t.Rotation(), r3d * t.Translation().Vect());
}


//---I/O functions
// TODO - I/O should be put in the manipulator form

/**
   print the 12 components of the Transform3D
 */
template <class T>
std::ostream &operator<<(std::ostream &os, const Transform3D<T> &t)
{
   // TODO - this will need changing for machine-readable issues
   //        and even the human readable form needs formatting improvements

   T m[12];
   t.GetComponents(m, m + 12);
   os << "\n" << m[0] << "  " << m[1] << "  " << m[2] << "  " << m[3];
   os << "\n" << m[4] << "  " << m[5] << "  " << m[6] << "  " << m[7];
   os << "\n" << m[8] << "  " << m[9] << "  " << m[10] << "  " << m[11] << "\n";
   return os;
}

} // end namespace Impl

// typedefs for double and float versions
typedef Impl::Transform3D<double> Transform3D;
typedef Impl::Transform3D<float>  Transform3DF;

} // end namespace Math

} // end namespace ROOT


#endif /* ROOT_Math_GenVector_Transform3D */
