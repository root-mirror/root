// @(#)root/gl:$Id$
// Author:  Richard Maunder  25/05/2005

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGLUtil
#define ROOT_TGLUtil

#include "Rtypes.h"
#include "TError.h"

#include <vector>
#include <cmath>
#include <cassert>
#include <utility>

class TString;
class TGLBoundingBox;
class TGLCamera;

class TAttMarker;
class TAttLine;

class GLUtesselator;

namespace Rgl
{
   enum EOverlap
   {
      kInside = 0,
      kPartial,
      kOutside
   };
}

enum EGLCoordType
{
   kGLCartesian,
   kGLPolar,
   kGLCylindrical,
   kGLSpherical
};

enum EGLPlotType
{
   kGLLegoPlot,
   kGLSurfacePlot,
   kGLBoxPlot,
   kGLTF3Plot,
   kGLStackPlot,
   kGLParametricPlot,
   kGLIsoPlot,
   kGL5D,
   kGLTH3Composition,
   kGLVoxel,
   kGLDefaultPlot
};


// TODO: Split these into own h/cxx files - too long now!

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGLVertex3                                                           //
//                                                                      //
// 3 component (x/y/z) vertex class                                     //
//                                                                      //
// This is part of collection of utility classes for GL in TGLUtil.h/cxx//
// These provide const and non-const accessors Arr() / CArr() to a GL   //
// compatible internal field - so can be used directly with OpenGL C API//
// calls. They are not intended to be fully featured just provide       //
// minimum required.                                                    //
//////////////////////////////////////////////////////////////////////////

class TGLVector3; // Forward declare for Shift()

class TGLVertex3
{
protected:
   // Fields
   Bool_t ValidIndex(UInt_t index) const { return (index < 3); }
   Double_t fVals[3];

public:
   TGLVertex3();
   TGLVertex3(Double_t x, Double_t y, Double_t z);
   TGLVertex3(Double_t* v);
   TGLVertex3(const TGLVertex3 & other);
   ~TGLVertex3();

   Bool_t       operator == (const TGLVertex3 & rhs) const;
   TGLVertex3 & operator =  (const TGLVertex3 & rhs);
   TGLVertex3 & operator *= (Double_t f);
   TGLVertex3   operator -  () const;
   const TGLVertex3 & operator -= (const TGLVector3 & val);
   const TGLVertex3 & operator += (const TGLVector3 & val);

   // Manipulators
   void Fill(Double_t val);
   void Set(Double_t x, Double_t y, Double_t z);
   void Set(const Double_t* xyz);
   void Set(const TGLVertex3 & other);
   void Shift(TGLVector3 & shift);
   void Shift(Double_t xDelta, Double_t yDelta, Double_t zDelta);
   void Negate();

   void Minimum(const TGLVertex3 & other);
   void Maximum(const TGLVertex3 & other);

   // Accessors
         Double_t & operator [] (Int_t index);
   const Double_t & operator [] (Int_t index) const;
   Double_t   X() const { return fVals[0]; }
   Double_t & X()       { return fVals[0]; }
   Double_t   Y() const { return fVals[1]; }
   Double_t & Y()       { return fVals[1]; }
   Double_t   Z() const { return fVals[2]; }
   Double_t & Z()       { return fVals[2]; }

   const Double_t * CArr() const { return fVals; }
   Double_t *       Arr()        { return fVals; }

   void Dump() const;

   ClassDefNV(TGLVertex3,1); // GL 3 component vertex helper/wrapper class
};

//______________________________________________________________________________
inline TGLVertex3 operator*(Double_t f, const TGLVertex3& v)
{
   return TGLVertex3(f*v.X(), f*v.Y(), f*v.Z());
}

//______________________________________________________________________________
inline void TGLVertex3::Negate()
{
   fVals[0] = -fVals[0];
   fVals[1] = -fVals[1];
   fVals[2] = -fVals[2];
}

//______________________________________________________________________________
inline Bool_t TGLVertex3::operator == (const TGLVertex3 & rhs) const
{
   return (fVals[0] == rhs.fVals[0] && fVals[1] == rhs.fVals[1] && fVals[2] == rhs.fVals[2]);
}

//______________________________________________________________________________
inline TGLVertex3 & TGLVertex3::operator = (const TGLVertex3 & rhs)
{
   // Check for self-assignment
   if (this != &rhs) {
      Set(rhs);
   }
   return *this;
}

// operator -= & operator += inline needs to be defered until full TGLVector3 definition

//______________________________________________________________________________
inline TGLVertex3 TGLVertex3::operator - () const
{
   return TGLVertex3(-fVals[0], -fVals[1], -fVals[2]);
}

//______________________________________________________________________________
inline TGLVertex3& TGLVertex3::operator *= (Double_t f)
{
   fVals[0] *= f;
   fVals[1] *= f;
   fVals[2] *= f;
   return *this;
}

//______________________________________________________________________________
inline Double_t & TGLVertex3::operator [] (Int_t index)
{
   /*if (!ValidIndex(index)) {
      assert(kFALSE);
      return fVals[0];
   } else {*/
      return fVals[index];
   //}
}

//______________________________________________________________________________
inline const Double_t& TGLVertex3::operator [] (Int_t index) const
{
   /*if (!ValidIndex(index)) {
      assert(kFALSE);
      return fVals[0];
   } else {*/
      return fVals[index];
   //}
}

//______________________________________________________________________________
inline void TGLVertex3::Fill(Double_t val)
{
   Set(val,val,val);
}

//______________________________________________________________________________
inline void TGLVertex3::Set(Double_t x, Double_t y, Double_t z)
{
   fVals[0]=x;
   fVals[1]=y;
   fVals[2]=z;
}

//______________________________________________________________________________
inline void TGLVertex3::Set(const Double_t* xyz)
{
   fVals[0]=xyz[0];
   fVals[1]=xyz[1];
   fVals[2]=xyz[2];
}

//______________________________________________________________________________
inline void TGLVertex3::Set(const TGLVertex3 & other)
{
   fVals[0]=other.fVals[0];
   fVals[1]=other.fVals[1];
   fVals[2]=other.fVals[2];
}


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGLVector3                                                           //
//                                                                      //
// 3 component (x/y/z) vector class                                     //
//                                                                      //
// This is part of collection of utility classes for GL in TGLUtil.h/cxx//
// These provide const and non-const accessors Arr() / CArr() to a GL   //
// compatible internal field - so can be used directly with OpenGL C API//
// calls. They are not intended to be fully featured just provide       //
// minimum required.                                                    //
//////////////////////////////////////////////////////////////////////////

class TGLVector3 : public TGLVertex3
{
public:
   TGLVector3() = default;
   TGLVector3(Double_t x, Double_t y, Double_t z);
   TGLVector3(const Double_t *src);

   TGLVector3& operator = (const TGLVertex3& v)
   { fVals[0] = v[0]; fVals[1] = v[1]; fVals[2] = v[2]; return *this; }

   TGLVector3 & operator /= (Double_t val);
   TGLVector3   operator -  () const;

   Double_t Mag() const;
   void     Normalise();

   ClassDefNV(TGLVector3,1); // GL 3 component vector helper/wrapper class
};

// Inline for TGLVertex3 requiring full TGLVector definition
//______________________________________________________________________________
inline const TGLVertex3 & TGLVertex3::operator -= (const TGLVector3 & vec)
{
   fVals[0] -= vec[0]; fVals[1] -= vec[1]; fVals[2] -= vec[2];
   return *this;
}

// Inline for TGLVertex3 requiring full TGLVector definition
//______________________________________________________________________________
inline const TGLVertex3 & TGLVertex3::operator += (const TGLVector3 & vec)
{
   fVals[0] += vec[0]; fVals[1] += vec[1]; fVals[2] += vec[2];
   return *this;
}

//______________________________________________________________________________
inline TGLVector3 & TGLVector3::operator /= (Double_t val)
{
   fVals[0] /= val;
   fVals[1] /= val;
   fVals[2] /= val;
   return *this;
}

//______________________________________________________________________________
inline TGLVector3 TGLVector3::operator - () const
{
   return TGLVector3(-fVals[0], -fVals[1], -fVals[2]);
}

//______________________________________________________________________________
inline Double_t TGLVector3::Mag() const
{
   return std::sqrt(fVals[0]*fVals[0] + fVals[1]*fVals[1] + fVals[2]*fVals[2]);
}

//______________________________________________________________________________
inline void TGLVector3::Normalise()
{
   Double_t mag = Mag();
   if ( mag == 0.0 ) {
      Error("TGLVector3::Normalise", "vector has zero magnitude");
      return;
   }
   fVals[0] /= mag;
   fVals[1] /= mag;
   fVals[2] /= mag;
}

//______________________________________________________________________________
inline Double_t Dot(const TGLVector3 & v1, const TGLVector3 & v2)
{
   return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

//______________________________________________________________________________
inline TGLVector3 Cross(const TGLVector3 & v1, const TGLVector3 & v2)
{
    return TGLVector3(v1[1]*v2[2] - v2[1]*v1[2],
                      v1[2]*v2[0] - v2[2]*v1[0],
                      v1[0]*v2[1] - v2[0]*v1[1]);
}

//______________________________________________________________________________
inline const TGLVector3 operator / (const TGLVector3 & vec, Double_t val)
{
   return TGLVector3(vec[0] / val, vec[1] / val, vec[2] / val);
}

//______________________________________________________________________________
inline const TGLVector3 operator * (const TGLVector3 & vec, Double_t val)
{
   return TGLVector3(vec[0] * val, vec[1] * val, vec[2] * val);
}

//______________________________________________________________________________
// Vertex + Vector => Vertex
inline TGLVertex3 operator + (const TGLVertex3 & vertex1, const TGLVector3 & vertex2)
{
   return TGLVertex3(vertex1[0] + vertex2[0], vertex1[1] + vertex2[1], vertex1[2] + vertex2[2]);
}

//______________________________________________________________________________
// Vertex - Vertex => Vector
inline TGLVector3 operator - (const TGLVertex3 & vertex1, const TGLVertex3 & vertex2)
{
   return TGLVector3(vertex1[0] - vertex2[0], vertex1[1] - vertex2[1], vertex1[2] - vertex2[2]);
}

//______________________________________________________________________________
// Vector + Vector => Vector
inline TGLVector3 operator + (const TGLVector3 & vector1, const TGLVector3 & vector2)
{
   return TGLVector3(vector1[0] + vector2[0], vector1[1] + vector2[1], vector1[2] + vector2[2]);
}

//______________________________________________________________________________
// Vector - Vector => Vector
inline TGLVector3 operator - (const TGLVector3 & vector1, const TGLVector3 & vector2)
{
   return TGLVector3(vector1[0] - vector2[0], vector1[1] - vector2[1], vector1[2] - vector2[2]);
}

//______________________________________________________________________________
// Dot-product
inline Double_t operator * (const TGLVector3 & a, const TGLVector3 & b)
{
   return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGLLine3                                                             //
//                                                                      //
// 3D space, fixed length, line class, with direction / length 'vector',//
// passing through point 'vertex'. Just wraps a TGLVector3 / TGLVertex3 //
// pair.                                                                //
//////////////////////////////////////////////////////////////////////////

class TGLLine3
{
private:
   // Fields
   TGLVertex3 fVertex; //! Start vertex of line
   TGLVector3 fVector; //! Vector of line from fVertex

public:
   TGLLine3(const TGLVertex3 & start, const TGLVertex3 & end);
   TGLLine3(const TGLVertex3 & start, const TGLVector3 & vector);
   ~TGLLine3() = default;

   void Set(const TGLVertex3 & start, const TGLVertex3 & end);
   void Set(const TGLVertex3 & start, const TGLVector3 & vector);

   // Bitwise copy constructor and = operator are fine

   // Accessors
   const TGLVertex3 & Start()  const { return fVertex; }
   const TGLVertex3   End()    const { return fVertex + fVector; }
   const TGLVector3 & Vector() const { return fVector; }

   // Debug
   void Draw() const;

   ClassDefNV(TGLLine3,0); // GL line wrapper class
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGLRect                                                              //
//                                                                      //
// Viewport (pixel base) 2D rectangle class                             //
//////////////////////////////////////////////////////////////////////////

class TGLRect
{
private:
   // Fields
   Int_t    fX, fY;           //! Corner
   Int_t    fWidth, fHeight;  //! Positive width/height

public:
   TGLRect();
   TGLRect(Int_t x, Int_t y, Int_t width, Int_t height);
   TGLRect(Int_t x, Int_t y, UInt_t width, UInt_t height);
   virtual ~TGLRect();

   // Bitwise copy const & =op are ok at present

   // Manipulators
   void Set(Int_t x, Int_t y, Int_t width, Int_t height);
   void SetCorner(Int_t x, Int_t y);
   void Offset(Int_t dX, Int_t dY);
   void Expand(Int_t x, Int_t y);

   // Accessors
   const Int_t* CArr() const { return &fX; }
         Int_t* CArr()       { return &fX; }

   Int_t    X()       const { return fX; }
   Int_t &  X()             { return fX; }
   Int_t    Y()       const { return fY; }
   Int_t &  Y()             { return fY; }
   Int_t    Width()   const { return fWidth; }
   Int_t &  Width()         { return fWidth; }
   Int_t    Height()  const { return fHeight; }
   Int_t &  Height()        { return fHeight; }
   Int_t    CenterX() const { return fX + fWidth/2; }
   Int_t    CenterY() const { return fY + fHeight/2; }
   Int_t    Left()    const { return fX; }
   Int_t    Right()   const { return fX + fWidth; }
   Int_t    Top()     const { return fY; }
   Int_t    Bottom()  const { return fY + fHeight; }

   Int_t Diagonal() const;
   Int_t Longest() const;

   Double_t Aspect() const;
   Rgl::EOverlap Overlap(const TGLRect & other) const;

   ClassDef(TGLRect,0); // GL rect helper/wrapper class
};

//______________________________________________________________________________
inline void TGLRect::Set(Int_t x, Int_t y, Int_t width, Int_t height)
{
   fX = x;
   fY = y;
   fWidth = width;
   fHeight = height;
}

//______________________________________________________________________________
inline void TGLRect::SetCorner(Int_t x, Int_t y)
{
   fX = x;
   fY = y;
}

//______________________________________________________________________________
inline void TGLRect::Offset(Int_t dX, Int_t dY)
{
   fX += dX;
   fY += dY;
}

//______________________________________________________________________________
inline Int_t TGLRect::Longest() const
{
   return fWidth > fHeight ? fWidth : fHeight;
}

//______________________________________________________________________________
inline Double_t TGLRect::Aspect() const
{
   // Return aspect ratio (width/height)
   if (fHeight == 0) {
      return 0.0;
   } else {
      return static_cast<Double_t>(fWidth) / static_cast<Double_t>(fHeight);
   }
}

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGLPlane                                                             //
//                                                                      //
// 3D plane class - of format Ax + By + Cz + D = 0                      //
//                                                                      //
// This is part of collection of simple utility classes for GL only in  //
// TGLUtil.h/cxx. These provide const and non-const accessors Arr() &   //
// CArr() to a GL compatible internal field - so can be used directly   //
// with OpenGL C API calls - which TVector3 etc cannot (easily).        //
// They are not intended to be fully featured just provide minimum      //
// required.                                                            //
//////////////////////////////////////////////////////////////////////////

class TGLPlane
{
private:
   // Fields
   Double_t fVals[4];

   // Methods
   void Normalise();

public:
   TGLPlane();
   TGLPlane(const TGLPlane & other);
   TGLPlane(Double_t a, Double_t b, Double_t c, Double_t d);
   TGLPlane(Double_t eq[4]);
   TGLPlane(const TGLVector3 & norm, const TGLVertex3 & point);
   TGLPlane(const TGLVertex3 & p1, const TGLVertex3 & p2, const TGLVertex3 & p3);
   ~TGLPlane() = default;

   TGLPlane &operator=(const TGLPlane &src);

   // Manipulators
   void Set(const TGLPlane & other);
   void Set(Double_t a, Double_t b, Double_t c, Double_t d);
   void Set(Double_t eq[4]);
   void Set(const TGLVector3 & norm, const TGLVertex3 & point);
   void Set(const TGLVertex3 & p1, const TGLVertex3 & p2, const TGLVertex3 & p3);
   void Negate();

   // Accessors
   Double_t A() const { return fVals[0]; }
   Double_t B() const { return fVals[1]; }
   Double_t C() const { return fVals[2]; }
   Double_t D() const { return fVals[3]; }

   TGLVector3 Norm() const { return TGLVector3( fVals[0], fVals[1], fVals[2]); }
   Double_t DistanceTo(const TGLVertex3 & vertex) const;
   TGLVertex3 NearestOn(const TGLVertex3 & point) const;

   // Internal data accessors - for GL API
   const Double_t * CArr() const { return fVals; }
   Double_t * Arr() { return fVals; }

   void Dump() const;

   ClassDefNV(TGLPlane,0); // GL plane helper/wrapper class
};

typedef std::vector<TGLPlane>                 TGLPlaneSet_t;
typedef std::vector<TGLPlane>::iterator       TGLPlaneSet_i;
typedef std::vector<TGLPlane>::const_iterator TGLPlaneSet_ci;

// Some free functions for planes
std::pair<Bool_t, TGLLine3>   Intersection(const TGLPlane & p1, const TGLPlane & p2);
std::pair<Bool_t, TGLVertex3> Intersection(const TGLPlane & p1, const TGLPlane & p2, const TGLPlane & p3);
std::pair<Bool_t, TGLVertex3> Intersection(const TGLPlane & plane, const TGLLine3 & line, Bool_t extend);


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGLMatrix                                                            //
//                                                                      //
// 16 component (4x4) transform matrix - column MAJOR as per GL.        //
// Provides limited support for adjusting the translation, scale and    //
// rotation components.                                                 //
//                                                                      //
// This is part of collection of simple utility classes for GL only in  //
// TGLUtil.h/cxx. These provide const and non-const accessors Arr() &   //
// CArr() to a GL compatible internal field - so can be used directly   //
// with OpenGL C API calls - which TVector3 etc cannot (easily).        //
// They are not intended to be fully featured just provide minimum      //
// required.                                                            //
//////////////////////////////////////////////////////////////////////////

class TGLMatrix
{
private:
   // Fields
   Double_t fVals[16]; // Column MAJOR as per OGL

   // Methods
   Bool_t ValidIndex(UInt_t index) const { return (index < 16); }

public:
   TGLMatrix();
   TGLMatrix(Double_t x, Double_t y, Double_t z);
   TGLMatrix(const TGLVertex3 & translation);
   TGLMatrix(const TGLVertex3 & origin, const TGLVector3 & zAxis, const TGLVector3 & xAxis);
   TGLMatrix(const TGLVertex3 & origin, const TGLVector3 & zAxis);
   TGLMatrix(const Double_t vals[16]);
   TGLMatrix(const TGLMatrix & other);
   virtual ~TGLMatrix();

   // Operators
   TGLMatrix & operator =(const TGLMatrix & rhs);
   Double_t  & operator [] (Int_t index);
   Double_t    operator [] (Int_t index) const;

   void MultRight(const TGLMatrix & rhs);
   void MultLeft (const TGLMatrix & lhs);
   TGLMatrix & operator*=(const TGLMatrix & rhs) { MultRight(rhs); return *this; }

   // Manipulators
   void Set(const TGLVertex3 & origin, const TGLVector3 & zAxis, const TGLVector3 & xAxis = 0);
   void Set(const Double_t vals[16]);
   void SetIdentity();

   void SetTranslation(Double_t x, Double_t y, Double_t z);
   void SetTranslation(const TGLVertex3 & translation);

   void Translate(const TGLVector3 & vect);
   void MoveLF(Int_t ai, Double_t amount);
   void Move3LF(Double_t x, Double_t y, Double_t z);

   void Scale(const TGLVector3 & scale);
   void Rotate(const TGLVertex3 & pivot, const TGLVector3 & axis, Double_t angle);
   void RotateLF(Int_t i1, Int_t i2, Double_t amount);
   void RotatePF(Int_t i1, Int_t i2, Double_t amount);
   void TransformVertex(TGLVertex3 & vertex) const;
   void Transpose3x3();
   Double_t Invert();

   // Accesors
   TGLVector3  GetTranslation() const;
   TGLVector3  GetScale() const;
   Bool_t      IsScalingForRender() const;

   void SetBaseVec(Int_t b, Double_t x, Double_t y, Double_t z);
   void SetBaseVec(Int_t b, const TGLVector3& v);
   void SetBaseVec(Int_t b, Double_t* x);

   TGLVector3 GetBaseVec(Int_t b) const;
   void       GetBaseVec(Int_t b, TGLVector3& v) const;
   void       GetBaseVec(Int_t b, Double_t* x) const;

   TGLVector3 Multiply(const TGLVector3& v, Double_t w=1) const;
   TGLVector3 Rotate(const TGLVector3& v) const;
   void       MultiplyIP(TGLVector3& v, Double_t w=1) const;
   void       RotateIP(TGLVector3& v) const;

   // Internal data accessors - for GL API
   const Double_t * CArr() const { return fVals; }
   Double_t * Arr() { return fVals; }

   void Dump() const;

   ClassDef(TGLMatrix,1); // GL matrix helper/wrapper class
};

//______________________________________________________________________________
inline TGLMatrix & TGLMatrix::operator =(const TGLMatrix & rhs)
{
   // Check for self-assignment
   if (this != &rhs) {
      Set(rhs.fVals);
   }
   return *this;
}

//______________________________________________________________________________
inline Double_t & TGLMatrix::operator [] (Int_t index)
{
   /*if (!ValidIndex(index)) {
      assert(kFALSE);
      return fVals[0];
   } else {*/
      return fVals[index];
   //}
}

//______________________________________________________________________________
inline Double_t TGLMatrix::operator [] (Int_t index) const
{
   /*if (!ValidIndex(index)) {
      assert(kFALSE);
      return fVals[0];
   } else {*/
      return fVals[index];
   //}
}

//______________________________________________________________________________
inline TGLMatrix operator * (const TGLMatrix & lhs, const TGLMatrix & rhs)
{
   TGLMatrix res;

   res[ 0] = rhs[ 0] * lhs[ 0] + rhs[ 1] * lhs[ 4] + rhs[ 2] * lhs[ 8] + rhs[ 3] * lhs[12];
   res[ 1] = rhs[ 0] * lhs[ 1] + rhs[ 1] * lhs[ 5] + rhs[ 2] * lhs[ 9] + rhs[ 3] * lhs[13];
   res[ 2] = rhs[ 0] * lhs[ 2] + rhs[ 1] * lhs[ 6] + rhs[ 2] * lhs[10] + rhs[ 3] * lhs[14];
   res[ 3] = rhs[ 0] * lhs[ 3] + rhs[ 1] * lhs[ 7] + rhs[ 2] * lhs[11] + rhs[ 3] * lhs[15];

   res[ 4] = rhs[ 4] * lhs[ 0] + rhs[ 5] * lhs[ 4] + rhs[ 6] * lhs[ 8] + rhs[ 7] * lhs[12];
   res[ 5] = rhs[ 4] * lhs[ 1] + rhs[ 5] * lhs[ 5] + rhs[ 6] * lhs[ 9] + rhs[ 7] * lhs[13];
   res[ 6] = rhs[ 4] * lhs[ 2] + rhs[ 5] * lhs[ 6] + rhs[ 6] * lhs[10] + rhs[ 7] * lhs[14];
   res[ 7] = rhs[ 4] * lhs[ 3] + rhs[ 5] * lhs[ 7] + rhs[ 6] * lhs[11] + rhs[ 7] * lhs[15];

   res[ 8] = rhs[ 8] * lhs[ 0] + rhs[ 9] * lhs[ 4] + rhs[10] * lhs[ 8] + rhs[11] * lhs[12];
   res[ 9] = rhs[ 8] * lhs[ 1] + rhs[ 9] * lhs[ 5] + rhs[10] * lhs[ 9] + rhs[11] * lhs[13];
   res[10] = rhs[ 8] * lhs[ 2] + rhs[ 9] * lhs[ 6] + rhs[10] * lhs[10] + rhs[11] * lhs[14];
   res[11] = rhs[ 8] * lhs[ 3] + rhs[ 9] * lhs[ 7] + rhs[10] * lhs[11] + rhs[11] * lhs[15];

   res[12] = rhs[12] * lhs[ 0] + rhs[13] * lhs[ 4] + rhs[14] * lhs[ 8] + rhs[15] * lhs[12];
   res[13] = rhs[12] * lhs[ 1] + rhs[13] * lhs[ 5] + rhs[14] * lhs[ 9] + rhs[15] * lhs[13];
   res[14] = rhs[12] * lhs[ 2] + rhs[13] * lhs[ 6] + rhs[14] * lhs[10] + rhs[15] * lhs[14];
   res[15] = rhs[12] * lhs[ 3] + rhs[13] * lhs[ 7] + rhs[14] * lhs[11] + rhs[15] * lhs[15];

   return res;
}

//______________________________________________________________________________
inline void TGLMatrix::SetBaseVec(Int_t b, Double_t x, Double_t y, Double_t z)
{
   Double_t* C = fVals + 4*--b;
   C[0] = x; C[1] = y; C[2] = z;
}

//______________________________________________________________________________
inline void TGLMatrix::SetBaseVec(Int_t b, const TGLVector3& v)
{
   Double_t* C = fVals + 4*--b;
   C[0] = v[0]; C[1] = v[1]; C[2] = v[2];
}

//______________________________________________________________________________
inline void TGLMatrix::SetBaseVec(Int_t b, Double_t* x)
{
   Double_t* C = fVals + 4*--b;
   C[0] = x[0]; C[1] = x[1]; C[2] = x[2];
}

//______________________________________________________________________________
inline TGLVector3 TGLMatrix::GetBaseVec(Int_t b) const
{
   return TGLVector3(&fVals[4*--b]);
}

//______________________________________________________________________________
inline void TGLMatrix::GetBaseVec(Int_t b, TGLVector3& v) const
{
   const Double_t* C = fVals + 4*--b;
   v[0] = C[0]; v[1] = C[1]; v[2] = C[2];
}

//______________________________________________________________________________
inline void TGLMatrix::GetBaseVec(Int_t b, Double_t* x) const
{
   const Double_t* C = fVals + 4*--b;
   x[0] = C[0], x[1] = C[1], x[2] = C[2];
}


//////////////////////////////////////////////////////////////////////////
//
// TGLColor
//
// Encapsulate color in preferred GL format - UChar_t RGBA array.
// Color index is also cached for easier interfacing with the
// traditional ROOT graphics.
//
//////////////////////////////////////////////////////////////////////////

class TGLColor
{
protected:
   UChar_t         fRGBA[4];
   mutable Short_t fIndex;

public:
   TGLColor();
   TGLColor(Int_t r, Int_t g, Int_t b, Int_t a=255);
   TGLColor(Float_t r, Float_t g, Float_t b, Float_t a=1);
   TGLColor(Color_t color_index, Char_t transparency=0);
   TGLColor(const TGLColor& c);

   TGLColor& operator=(const TGLColor& c);

   UChar_t*        Arr()       { return fRGBA; }
   const UChar_t* CArr() const { return fRGBA; }

   UChar_t GetRed()   const { return fRGBA[0]; }
   UChar_t GetGreen() const { return fRGBA[1]; }
   UChar_t GetBlue()  const { return fRGBA[2]; }
   UChar_t GetAlpha() const { return fRGBA[3]; }

   Color_t GetColorIndex()   const;
   Char_t  GetTransparency() const;

   void SetRed(Int_t v)   { fRGBA[0] = v; }
   void SetGreen(Int_t v) { fRGBA[1] = v; }
   void SetBlue(Int_t v)  { fRGBA[2] = v; }
   void SetAlpha(Int_t v) { fRGBA[3] = v; }

   void SetColor(Int_t r, Int_t g, Int_t b, Int_t a=255);
   void SetColor(Float_t r, Float_t g, Float_t b, Float_t a=1);
   void SetColor(Color_t color_index);
   void SetColor(Color_t color_index, Char_t transparency);
   void SetTransparency(Char_t transparency);

   TString AsString() const;

   ClassDefNV(TGLColor, 0); // Color in preferred GL format - RGBA.
};


//////////////////////////////////////////////////////////////////////////
//
// TGLColorSet
//
// A collection of colors used for OpenGL rendering.
//
//////////////////////////////////////////////////////////////////////////

class TGLColorSet
{
protected:
   TGLColor        fBackground;
   TGLColor        fForeground;
   TGLColor        fOutline;
   TGLColor        fMarkup;
   TGLColor        fSelection[5];   // Colors for shape-selection-levels

public:
   TGLColorSet();
   TGLColorSet(const TGLColorSet& s);
   ~TGLColorSet() = default;

   TGLColorSet& operator=(const TGLColorSet& s);

   TGLColor& Background()       { return fBackground; }
   TGLColor& Foreground()       { return fForeground; }
   TGLColor& Outline()          { return fOutline; }
   TGLColor& Markup()           { return fMarkup;  }
   TGLColor& Selection(Int_t i) { return fSelection[i]; }

   const TGLColor& Background()       const { return fBackground; }
   const TGLColor& Foreground()       const { return fForeground; }
   const TGLColor& Outline()          const { return fOutline; }
   const TGLColor& Markup()           const { return fMarkup;  }
   const TGLColor& Selection(Int_t i) const { return fSelection[i]; }

   void StdDarkBackground();
   void StdLightBackground();

   ClassDefNV(TGLColorSet, 0); // Collection of colors used for GL rendering.
};

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGLUtil                                                              //
//                                                                      //
// Wrapper class for various misc static functions - error checking,    //
// draw helpers etc.                                                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

class TGLUtil
{
public:
   class TColorLocker
   {
   public:
      TColorLocker()          { LockColor();   }
      virtual ~TColorLocker() { UnlockColor(); }

      ClassDef(TColorLocker,0); // Lock/unlock color in constructor/destructor.
   };

   class TDrawQualityModifier
   {
      Int_t fOldQuality;
   public:
      TDrawQualityModifier(Int_t dq) :
      fOldQuality(GetDrawQuality()) {SetDrawQuality(dq); }

      virtual ~TDrawQualityModifier()
      { SetDrawQuality(fOldQuality); }

      ClassDef(TDrawQualityModifier,0); // Set/restore draw quality in constructor/destructor.
   };

   class TDrawQualityScaler
   {
      Int_t fOldQuality;
   public:
      TDrawQualityScaler(Float_t fac) :
      fOldQuality(GetDrawQuality()) {SetDrawQuality((Int_t)(fac*fOldQuality)); }

      virtual ~TDrawQualityScaler()
      { SetDrawQuality(fOldQuality); }

      ClassDef(TDrawQualityScaler,0); // Multiply/restore draw quality in constructor/destructor.
   };

private:
   static UInt_t fgDefaultDrawQuality;
   static UInt_t fgDrawQuality;

   static UInt_t fgColorLockCount;

   static Float_t fgPointSize;
   static Float_t fgLineWidth;
   static Float_t fgPointSizeScale;
   static Float_t fgLineWidthScale;

   static Float_t fgScreenScalingFactor;
   static Float_t fgPointLineScalingFactor;
   static Int_t   fgPickingRadius;

   TGLUtil(const TGLUtil&);            // Not implemented.
   TGLUtil& operator=(const TGLUtil&); // Not implemented.

public:
   virtual ~TGLUtil() {}
   static void InitializeIfNeeded();

   // Error checking
   static Int_t  CheckError(const char * loc);

   // Polygon tesselator for direct drawing
   static GLUtesselator* GetDrawTesselator3fv();
   static GLUtesselator* GetDrawTesselator4fv();
   static GLUtesselator* GetDrawTesselator3dv();
   static GLUtesselator* GetDrawTesselator4dv();

   // Some simple shape drawing utils
   enum ELineHeadShape { kLineHeadNone, kLineHeadArrow, kLineHeadBox };
   enum EAxesType      { kAxesNone, kAxesEdge, kAxesOrigin };

   static UInt_t GetDrawQuality();
   static void   SetDrawQuality(UInt_t dq);
   static void   ResetDrawQuality();
   static UInt_t GetDefaultDrawQuality();
   static void   SetDefaultDrawQuality(UInt_t dq);

   static UInt_t LockColor();
   static UInt_t UnlockColor();
   static Bool_t IsColorLocked();

   static void Color(const TGLColor& color);
   static void ColorAlpha(const TGLColor& color, UChar_t alpha);
   static void ColorAlpha(const TGLColor& color, Float_t alpha);
   static void ColorAlpha(Color_t color_index, Float_t alpha=1);
   static void ColorTransparency(Color_t color_index, Char_t transparency=0);
   static void Color3ub(UChar_t r, UChar_t g, UChar_t b);
   static void Color4ub(UChar_t r, UChar_t g, UChar_t b, UChar_t a);
   static void Color3ubv(const UChar_t* rgb);
   static void Color4ubv(const UChar_t* rgba);
   static void Color3f(Float_t r, Float_t g, Float_t b);
   static void Color4f(Float_t r, Float_t g, Float_t b, Float_t a);
   static void Color3fv(const Float_t* rgb);
   static void Color4fv(const Float_t* rgba);

   // Coordinate conversion and extra scaling (needed for osx retina)
   static void    PointToViewport(Int_t& x, Int_t& y);
   static void    PointToViewport(Int_t& x, Int_t& y, Int_t& w, Int_t& h);
   static Float_t GetScreenScalingFactor();
   static Float_t GetPointLineScalingFactor();
   static Int_t   GetPickingRadius();

   static Float_t GetPointSizeScale();
   static void    SetPointSizeScale(Float_t scale);
   static Float_t GetLineWidthScale();
   static void    SetLineWidthScale(Float_t scale);

   static void    PointSize(Float_t point_size);
   static void    LineWidth(Float_t line_width);

   static Float_t PointSize();
   static Float_t LineWidth();

   static void BeginExtendPickRegion(Float_t scale);
   static void EndExtendPickRegion();

   static void RenderPolyMarkers(const TAttMarker& marker, Char_t transp,
                                 Float_t* p, Int_t n,
                                 Int_t pick_radius=0, Bool_t selection=kFALSE,
                                 Bool_t sec_selection=kFALSE);

   static void RenderPolyMarkers(const TAttMarker &marker, const std::vector<Double_t> &points,
                                 Double_t dX, Double_t dY, Double_t dZ);

   static void RenderPoints(const TAttMarker& marker,
                            Float_t* p, Int_t n,
                            Int_t pick_radius=0, Bool_t selection=kFALSE,
                            Bool_t sec_selection=kFALSE);

   static void RenderPoints(const TAttMarker& marker,
                            const std::vector<Double_t> &points);

   static void RenderCrosses(const TAttMarker& marker,
                             Float_t* p, Int_t n,
                             Bool_t sec_selection=kFALSE);

   static void RenderCrosses(const TAttMarker& marker,
                             const std::vector<Double_t> &points,
                             Double_t dX, Double_t dY, Double_t dZ);

   static void RenderPolyLine(const TAttLine& aline, Char_t transp,
                              Float_t* p, Int_t n,
                              Int_t pick_radius=0, Bool_t selection=kFALSE);

   static void BeginAttLine(const TAttLine& aline, Char_t transp,
                            Int_t pick_radius=0, Bool_t selection=kFALSE);
   static void EndAttLine(Int_t pick_radius=0, Bool_t selection=kFALSE);

   // TODO: These draw routines should take LOD hints
   static void SetDrawColors(const UChar_t rgba[4]);
   static void DrawSphere(const TGLVertex3 & position, Double_t radius, const UChar_t rgba[4]);
   static void DrawLine(const TGLLine3 & line, ELineHeadShape head, Double_t size, const UChar_t rgba[4]);
   static void DrawLine(const TGLVertex3 & start, const TGLVector3 & vector, ELineHeadShape head,
                        Double_t size, const UChar_t rgba[4]);
   static void DrawRing(const TGLVertex3 & center, const TGLVector3 & normal,
                        Double_t radius, const UChar_t* rgba);

   static void DrawReferenceMarker(const TGLCamera  & camera,
                                   const TGLVertex3 & pos,
                                         Float_t      radius = 3,
                                   const UChar_t    * rgba   = 0);
   static void DrawSimpleAxes(const TGLCamera      & camera,
                              const TGLBoundingBox & bbox,
                                    Int_t            axesType);
   static void DrawNumber(const TString    & num,
                          const TGLVertex3 & pos,
                                Bool_t       center = kFALSE);

   // Frequently used colors.
   static const UChar_t fgRed[4];
   static const UChar_t fgGreen[4];
   static const UChar_t fgBlue[4];
   static const UChar_t fgYellow[4];
   static const UChar_t fgWhite[4];
   static const UChar_t fgGrey[4];

   ClassDef(TGLUtil,0); // Wrapper class for misc GL pieces
};

/**************************************************************************/

class TGLCapabilitySwitch
{
private:
   TGLCapabilitySwitch(const TGLCapabilitySwitch &);
   TGLCapabilitySwitch &operator = (const TGLCapabilitySwitch &);

   Int_t    fWhat;
   Bool_t   fState;
   Bool_t   fFlip;

   void SetState(Bool_t s);

public:
   TGLCapabilitySwitch(Int_t what, Bool_t state);
   ~TGLCapabilitySwitch();
};

class TGLCapabilityEnabler
{
private:
   TGLCapabilityEnabler(const TGLCapabilityEnabler &);
   TGLCapabilityEnabler &operator = (const TGLCapabilityEnabler &);

   Int_t    fWhat;
   Bool_t   fFlip;

public:
   TGLCapabilityEnabler(Int_t what, Bool_t state);
   ~TGLCapabilityEnabler();
};

class TGLFloatHolder
{
   TGLFloatHolder(const TGLFloatHolder&);            // Not implemented
   TGLFloatHolder& operator=(const TGLFloatHolder&); // Not implemented

   Int_t    fWhat;
   Float_t  fState;
   Bool_t   fFlip;
   void   (*fFoo)(Float_t);

public:
   TGLFloatHolder(Int_t what, Float_t state, void (*foo)(Float_t));
   ~TGLFloatHolder();
};

class TGLEnableGuard {
private:
   Int_t fCap;

public:
   TGLEnableGuard(Int_t cap);
   ~TGLEnableGuard();

private:
   TGLEnableGuard(const TGLEnableGuard &);
   TGLEnableGuard &operator = (const TGLEnableGuard &);
};

class TGLDisableGuard {
private:
   Int_t fCap;

public:
   TGLDisableGuard(Int_t cap);
   ~TGLDisableGuard();

private:
   TGLDisableGuard(const TGLDisableGuard &);
   TGLDisableGuard &operator = (const TGLDisableGuard &);
};

class TGLSelectionBuffer {
   std::vector<UChar_t> fBuffer;
   Int_t                fWidth;
   Int_t                fHeight;

public:
   TGLSelectionBuffer();
   virtual ~TGLSelectionBuffer();

   void           ReadColorBuffer(Int_t width, Int_t height);
   void           ReadColorBuffer(Int_t x, Int_t y, Int_t width, Int_t height);
   const UChar_t *GetPixelColor(Int_t px, Int_t py)const;

private:
   TGLSelectionBuffer(const TGLSelectionBuffer &);
   TGLSelectionBuffer &operator = (const TGLSelectionBuffer &);

   ClassDef(TGLSelectionBuffer, 0); //Holds color buffer content for selection
};

template<class T>
class TGL2DArray : public std::vector<T> {
private:
   Int_t fRowLen;
   Int_t fMaxRow;
   typedef typename std::vector<T>::size_type size_type;

public:
   TGL2DArray() : fRowLen(0), fMaxRow(0){}
   void SetMaxRow(Int_t max)
   {
      fMaxRow = max;
   }
   void SetRowLen(Int_t len)
   {
      fRowLen = len;
   }
   const T *operator [] (size_type ind)const
   {
      return &std::vector<T>::operator [](ind * fRowLen);
   }
   T *operator [] (size_type ind)
   {
      return &std::vector<T>::operator [] (ind * fRowLen);
   }
};

class TGLPlotCoordinates;
class TGLQuadric;
class TAxis;

namespace Rgl {

extern const Float_t gRedEmission[];
extern const Float_t gGreenEmission[];
extern const Float_t gBlueEmission[];
extern const Float_t gOrangeEmission[];
extern const Float_t gWhiteEmission[];
extern const Float_t gGrayEmission[];
extern const Float_t gNullEmission[];

typedef std::pair<Int_t, Int_t> BinRange_t;
typedef std::pair<Double_t, Double_t> Range_t;

void ObjectIDToColor(Int_t objectID, Bool_t highColor);
Int_t ColorToObjectID(const UChar_t *color, Bool_t highColor);
void DrawQuadOutline(const TGLVertex3 &v1, const TGLVertex3 &v2,
                     const TGLVertex3 &v3, const TGLVertex3 &v4);
void DrawQuadFilled(const TGLVertex3 &v0, const TGLVertex3 &v1,
                    const TGLVertex3 &v2, const TGLVertex3 &v3,
                    const TGLVector3 &normal);
void DrawQuadFilled(const Double_t *v0, const Double_t *v1,
                    const Double_t *v2, const Double_t *v3,
                    const Double_t *normal);

void DrawSmoothFace(const TGLVertex3 &v1, const TGLVertex3 &v2,
                     const TGLVertex3 &v3, const TGLVector3 &norm1,
                     const TGLVector3 &norm2, const TGLVector3 &norm3);
void DrawBoxFront(Double_t xMin, Double_t xMax, Double_t yMin, Double_t yMax,
                  Double_t zMin, Double_t zMax, Int_t fp);

void DrawTransparentBox(Double_t xMin, Double_t xMax, Double_t yMin, Double_t yMax,
                        Double_t zMin, Double_t zMax, Int_t fp);


void DrawBoxFrontTextured(Double_t xMin, Double_t xMax, Double_t yMin,
                          Double_t yMax, Double_t zMin, Double_t zMax,
                          Double_t tMin, Double_t tMax, Int_t front);

void DrawBoxWithGradientFill(Double_t y1, Double_t y2, Double_t x1, Double_t x2,
                             const Double_t *rgba1, const Double_t *rgba2);

void DrawQuadStripWithRadialGradientFill(unsigned nPoints, const Double_t *inner, const Double_t *innerRGBA,
                                        const Double_t *outer, const Double_t *outerRGBA);

#ifndef __CINT__
void DrawTrapezoidTextured(const Double_t ver[][2], Double_t zMin, Double_t zMax,
                           Double_t tMin, Double_t tMax);
void DrawTrapezoidTextured(const Double_t ver[][3], Double_t texMin, Double_t texMax);
void DrawTrapezoidTextured2(const Double_t ver[][2], Double_t zMin, Double_t zMax,
                            Double_t tMin, Double_t tMax);
#endif

void DrawCylinder(TGLQuadric *quadric, Double_t xMin, Double_t xMax, Double_t yMin,
                  Double_t yMax, Double_t zMin, Double_t zMax);
void DrawSphere(TGLQuadric *quadric, Double_t xMin, Double_t xMax, Double_t yMin,
                Double_t yMax, Double_t zMin, Double_t zMax);
void DrawError(Double_t xMin, Double_t xMax, Double_t yMin,
               Double_t yMax, Double_t zMin, Double_t zMax);

#ifndef __CINT__
void DrawTrapezoid(const Double_t ver[][2], Double_t zMin, Double_t zMax, Bool_t color = kTRUE);
void DrawTrapezoid(const Double_t ver[][3]);
#endif

void DrawAxes(Int_t frontPoint, const Int_t *viewport, const TGLVertex3 *box2D,
              const TGLPlotCoordinates *plotCoord, TAxis *xAxis, TAxis *yAxis,
              TAxis *zAxis);
void SetZLevels(TAxis *zAxis, Double_t zMin, Double_t zMax,
                Double_t zScale, std::vector<Double_t> &zLevels);

void DrawFaceTextured(const TGLVertex3 &v1, const TGLVertex3 &v2, const TGLVertex3 &v3,
                      Double_t t1, Double_t t2, Double_t t3, const TGLVector3 &norm1,
                      const TGLVector3 &norm2, const TGLVector3 &norm3);
void DrawFaceTextured(const TGLVertex3 &v1, const TGLVertex3 &v2, const TGLVertex3 &v3,
                      Double_t t1, Double_t t2, Double_t t3, Double_t z, const TGLVector3 &planeNormal);
void GetColor(Float_t v, Float_t vmin, Float_t vmax, Int_t type, Float_t *rgba);

class TGuardBase {
private:
   mutable Bool_t fActive;

   TGuardBase &operator = (const TGuardBase &rhs);
protected:
   TGuardBase()
      : fActive(kTRUE)
   {
   }
   TGuardBase(const TGuardBase &rhs)
      : fActive(kTRUE)
   {
      rhs.fActive = kFALSE;
   }

   Bool_t IsActive()const
   {
      return fActive;
   }

public:
   void Stop()const
   {
      fActive = kFALSE;
   }
};

template<class Func, class Arg>
class TOneArgGuard : public TGuardBase {
private:
   Func fFunc;
   Arg  fArg;
public:
   TOneArgGuard(Func f, Arg a)
      : fFunc(f), fArg(a)
   {
   }
   ~TOneArgGuard()
   {
      if (IsActive())
         fFunc(fArg);
   }
};

template<class Func, class Arg1, class Arg2>
class TTwoArgsGuard : public TGuardBase {
private:
   Func fFunc;
   Arg1 fArg1;
   Arg2 fArg2;

public:
   TTwoArgsGuard(Func f, Arg1 a1, Arg2 a2)
      : fFunc(f), fArg1(a1), fArg2(a2)
   {
   }
   ~TTwoArgsGuard()
   {
      if (IsActive())
         fFunc(fArg1, fArg2);
   }
};

template<class Func, class Arg>
TOneArgGuard<Func, Arg> make_guard(Func f, Arg a)
{
   return TOneArgGuard<Func, Arg>(f, a);
}

template<class Func, class Arg1, class Arg2>
TTwoArgsGuard<Func, Arg1, Arg2> make_guard(Func f, Arg1 a1, Arg2 a2)
{
   return TTwoArgsGuard<Func, Arg1, Arg2>(f, a1, a2);
}

}//namespace Rgl.

class TGLLevelPalette {
private:
   std::vector<UChar_t>         fTexels;
   const std::vector<Double_t> *fContours;
   UInt_t                       fPaletteSize;
   mutable UInt_t               fTexture;
   Int_t                        fMaxPaletteSize;
   Rgl::Range_t                 fZRange;

   TGLLevelPalette(const TGLLevelPalette&);    // Not implemented
   TGLLevelPalette& operator=(const TGLLevelPalette&);  // Not implemented

public:
   TGLLevelPalette();

   Bool_t GeneratePalette(UInt_t paletteSize, const Rgl::Range_t &zRange, Bool_t checkSize = kTRUE);

   void   SetContours(const std::vector<Double_t> *contours);

   void   EnableTexture(Int_t mode)const;
   void   DisableTexture()const;

   Int_t  GetPaletteSize()const;

   Double_t       GetTexCoord(Double_t z)const;

   const UChar_t *GetColour(Double_t z)const;
   const UChar_t *GetColour(Int_t ind)const;
};

#endif // ROOT_TGLUtil
