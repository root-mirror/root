// @(#)root/eve:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TEveTrans.h"
#include "TEveUtil.h"

#include "TBuffer.h"
#include "TClass.h"
#include "TMath.h"

#include <iostream>
#include <cctype>

#define F00  0
#define F01  4
#define F02  8
#define F03 12

#define F10  1
#define F11  5
#define F12  9
#define F13 13

#define F20  2
#define F21  6
#define F22 10
#define F23 14

#define F30  3
#define F31  7
#define F32 11
#define F33 15

/** \class TEveTrans
\ingroup TEve
TEveTrans is a 4x4 transformation matrix for homogeneous coordinates
stored internally in a column-major order to allow direct usage by
GL. The element type is Double32_t as statically the floats would
be precise enough but continuous operations on the matrix must
retain precision of column vectors.

Cartan angles are stored in fA[1-3] (+z, -y, +x). They are
recalculated on demand.

Direct  element access (first two should be used with care):
  - operator[i]    direct access to elements,   i:0->15
  - CM(i,j)        element 4*j + i;           i,j:0->3    { CM ~ c-matrix }
  - operator(i,j)  element 4*(j-1) + i - 1    i,j:1->4

Column-vector access:
USet Get/SetBaseVec(), Get/SetPos() and Arr[XYZT]() methods.

For all methods taking the matrix indices:
1->X, 2->Y, 3->Z; 4->Position (if applicable). 0 reserved for time.

Shorthands in method-names:
LF ~ LocalFrame; PF ~ ParentFrame; IP ~ InPlace
*/

ClassImp(TEveTrans);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor.

TEveTrans::TEveTrans() :
   TObject(),
   fA1(0), fA2(0), fA3(0), fAsOK(kFALSE),
   fUseTrans (kTRUE),
   fEditTrans(kFALSE),
   fEditRotation(kTRUE),
   fEditScale(kTRUE)
{
   UnitTrans();
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TEveTrans::TEveTrans(const TEveTrans& t) :
   TObject(),
   fA1(t.fA1), fA2(t.fA2), fA3(t.fA3), fAsOK(t.fAsOK),
   fUseTrans (t.fUseTrans),
   fEditTrans(t.fEditTrans),
   fEditRotation(kTRUE),
   fEditScale(kTRUE)
{
   SetTrans(t, kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TEveTrans::TEveTrans(const Double_t arr[16]) :
   TObject(),
   fA1(0), fA2(0), fA3(0), fAsOK(kFALSE),
   fUseTrans (kTRUE),
   fEditTrans(kFALSE),
   fEditRotation(kTRUE),
   fEditScale(kTRUE)
{
   SetFromArray(arr);
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TEveTrans::TEveTrans(const Float_t arr[16]) :
   TObject(),
   fA1(0), fA2(0), fA3(0), fAsOK(kFALSE),
   fUseTrans (kTRUE),
   fEditTrans(kFALSE),
   fEditRotation(kTRUE),
   fEditScale(kTRUE)
{
   SetFromArray(arr);
}

////////////////////////////////////////////////////////////////////////////////
/// Reset matrix to unity.

void TEveTrans::UnitTrans()
{
   memset(fM, 0, 16*sizeof(Double_t));
   fM[F00] = fM[F11] = fM[F22] = fM[F33] = 1;
   fA1 = fA2 = fA3 = 0;
   fAsOK = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Reset matrix to zero, only the perspective scaling is set to w
/// (1 by default).

void TEveTrans::ZeroTrans(Double_t w)
{
   memset(fM, 0, 16*sizeof(Double_t));
   fM[F33] = w;
   fA1 = fA2 = fA3 = 0;
   fAsOK = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Reset rotation part of the matrix to unity.

void TEveTrans::UnitRot()
{
   memset(fM, 0, 12*sizeof(Double_t));
   fM[F00] = fM[F11] = fM[F22] = 1;
   fA1 = fA2 = fA3 = 0;
   fAsOK = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set matrix from another,

void TEveTrans::SetTrans(const TEveTrans& t, Bool_t copyAngles)
{
   memcpy(fM, t.fM, sizeof(fM));
   if (copyAngles && t.fAsOK) {
      fAsOK = kTRUE;
      fA1 = t.fA1; fA2 = t.fA2; fA3 = t.fA3;
   } else {
      fAsOK = kFALSE;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set matrix from Double_t array.

void TEveTrans::SetFromArray(const Double_t arr[16])
{
   for(Int_t i=0; i<16; ++i) fM[i] = arr[i];
   fAsOK = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set matrix from Float_t array.

void TEveTrans::SetFromArray(const Float_t arr[16])
{
   for(Int_t i=0; i<16; ++i) fM[i] = arr[i];
   fAsOK = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Setup the matrix as an elementary rotation.
/// Optimized versions of left/right multiplication with an elementary
/// rotation matrix are implemented in RotatePF/RotateLF.
/// Expects identity matrix.

void TEveTrans::SetupRotation(Int_t i, Int_t j, Double_t f)
{
   if(i == j) return;
   TEveTrans& t = *this;
   t(i,i) = t(j,j) = TMath::Cos(f);
   Double_t s = TMath::Sin(f);
   t(i,j) = -s; t(j,i) = s;
   fAsOK = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// A function for creating a rotation matrix that rotates a vector called
/// "from" into another vector called "to".
/// Input : from[3], to[3] which both must be *normalized* non-zero vectors
/// Output: mtx[3][3] -- a 3x3 matrix in column-major form
///
/// Authors: Tomas Möller, John Hughes
///          "Efficiently Building a Matrix to Rotate One Vector to Another"
///          Journal of Graphics Tools, 4(4):1-4, 1999

void TEveTrans::SetupFromToVec(const TEveVector& from, const TEveVector& to)
{
   static const float kFromToEpsilon = 0.000001f;

   ZeroTrans();

   Float_t e, f;
   e = from.Dot(to);
   f = (e < 0.0f) ? -e : e;

   if (f > 1.0f - kFromToEpsilon) /* "from" and "to"-vector almost parallel */
   {
      TEveVector u, v;       /* temporary storage vectors */
      TEveVector x;          /* vector most nearly orthogonal to "from" */
      Float_t    c1, c2, c3; /* coefficients for later use */

      x.fX = (from.fX > 0.0f) ? from.fX : -from.fX;
      x.fY = (from.fY > 0.0f) ? from.fY : -from.fY;
      x.fZ = (from.fZ > 0.0f) ? from.fZ : -from.fZ;

      if (x.fX < x.fY)
      {
         if (x.fX < x.fZ) {
            x.fX = 1.0f; x.fY = x.fZ = 0.0f;
         } else {
            x.fZ = 1.0f; x.fX = x.fY = 0.0f;
         }
      }
      else
      {
         if (x.fY < x.fZ) {
            x.fY = 1.0f; x.fX = x.fZ = 0.0f;
         } else {
            x.fZ = 1.0f; x.fX = x.fY = 0.0f;
         }
      }

      u.Sub(x, from);
      v.Sub(x, to);

      c1 = 2.0f / u.Mag2();
      c2 = 2.0f / v.Mag2();
      c3 = c1 * c2  * u.Dot(v);

      for (int i = 0; i < 3; i++) {
         for (int j = 0; j < 3; j++) {
            CM(i, j) =  - c1 * u[i] * u[j]
               - c2 * v[i] * v[j]
               + c3 * v[i] * u[j];
         }
         CM(i, i) += 1.0;
      }
   }
   else  /* the most common case, unless "from"="to", or "from"=-"to" */
   {
      TEveVector v = from.Cross(to);

      Float_t h, hvx, hvz, hvxy, hvxz, hvyz;
      h   = 1.0f/(1.0f + e);
      hvx = h * v.fX;
      hvz = h * v.fZ;
      hvxy = hvx * v.fY;
      hvxz = hvx * v.fZ;
      hvyz = hvz * v.fY;

      CM(0, 0) = e + hvx * v.fX;
      CM(0, 1) = hvxy - v.fZ;
      CM(0, 2) = hvxz + v.fY;

      CM(1, 0) = hvxy + v.fZ;
      CM(1, 1) = e + h * v.fY * v.fY;
      CM(1, 2) = hvyz - v.fX;

      CM(2, 0) = hvxz - v.fY;
      CM(2, 1) = hvyz + v.fX;
      CM(2, 2) = e + hvz * v.fZ;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Multiply from left: this = t * this.

void TEveTrans::MultLeft(const TEveTrans& t)
{
   Double_t  buf[4];
   Double_t* col = fM;
   for(int c=0; c<4; ++c, col+=4) {
      const Double_t* row = t.fM;
      for(int r=0; r<4; ++r, ++row)
         buf[r] = row[0]*col[0] + row[4]*col[1] + row[8]*col[2] + row[12]*col[3];
      col[0] = buf[0]; col[1] = buf[1]; col[2] = buf[2]; col[3] = buf[3];
   }
   fAsOK = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Multiply from right: this = this * t.

void TEveTrans::MultRight(const TEveTrans& t)
{
   Double_t  buf[4];
   Double_t* row = fM;
   for(int r=0; r<4; ++r, ++row) {
      const Double_t* col = t.fM;
      for(int c=0; c<4; ++c, col+=4)
         buf[c] = row[0]*col[0] + row[4]*col[1] + row[8]*col[2] + row[12]*col[3];
      row[0] = buf[0]; row[4] = buf[1]; row[8] = buf[2]; row[12] = buf[3];
   }
   fAsOK = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy, multiply from right and return product.
/// Avoid unless necessary.

TEveTrans TEveTrans::operator*(const TEveTrans& t)
{
   TEveTrans b(*this);
   b.MultRight(t);
   return b;
}

////////////////////////////////////////////////////////////////////////////////
/// Transpose 3x3 rotation sub-matrix.

void TEveTrans::TransposeRotationPart()
{
   Double_t x;
   x = fM[F01]; fM[F01] = fM[F10]; fM[F10] = x;
   x = fM[F02]; fM[F02] = fM[F20]; fM[F20] = x;
   x = fM[F12]; fM[F12] = fM[F21]; fM[F21] = x;
   fAsOK = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Move in local-frame along axis with index ai.

void TEveTrans::MoveLF(Int_t ai, Double_t amount)
{
   const Double_t *col = fM + 4*--ai;
   fM[F03] += amount*col[0]; fM[F13] += amount*col[1]; fM[F23] += amount*col[2];
}

////////////////////////////////////////////////////////////////////////////////
/// General move in local-frame.

void TEveTrans::Move3LF(Double_t x, Double_t y, Double_t z)
{
   fM[F03] += x*fM[0] + y*fM[4] + z*fM[8];
   fM[F13] += x*fM[1] + y*fM[5] + z*fM[9];
   fM[F23] += x*fM[2] + y*fM[6] + z*fM[10];
}

////////////////////////////////////////////////////////////////////////////////
/// Rotate in local frame. Does optimised version of MultRight.

void TEveTrans::RotateLF(Int_t i1, Int_t i2, Double_t amount)
{
   if(i1 == i2) return;
   // Algorithm: TEveTrans a; a.SetupRotation(i1, i2, amount); MultRight(a);
   // Optimized version:
   const Double_t cos = TMath::Cos(amount), sin = TMath::Sin(amount);
   Double_t  b1, b2;
   Double_t* row = fM;
   --i1 <<= 2; --i2 <<= 2; // column major
   for (int r=0; r<4; ++r, ++row) {
      b1 = cos*row[i1] + sin*row[i2];
      b2 = cos*row[i2] - sin*row[i1];
      row[i1] = b1; row[i2] = b2;
   }
   fAsOK = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Move in parent-frame along axis index ai.

void TEveTrans::MovePF(Int_t ai, Double_t amount)
{
   fM[F03 + --ai] += amount;
}

////////////////////////////////////////////////////////////////////////////////
/// General move in parent-frame.

void TEveTrans::Move3PF(Double_t x, Double_t y, Double_t z)
{
   fM[F03] += x;
   fM[F13] += y;
   fM[F23] += z;
}

////////////////////////////////////////////////////////////////////////////////
/// Rotate in parent frame. Does optimised version of MultLeft.

void TEveTrans::RotatePF(Int_t i1, Int_t i2, Double_t amount)
{
   if(i1 == i2) return;
   // Algorithm: TEveTrans a; a.SetupRotation(i1, i2, amount); MultLeft(a);

   // Optimized version:
   const Double_t cos = TMath::Cos(amount), sin = TMath::Sin(amount);
   Double_t  b1, b2;
   Double_t* col = fM;
   --i1; --i2;
   for(int c=0; c<4; ++c, col+=4) {
      b1 = cos*col[i1] - sin*col[i2];
      b2 = cos*col[i2] + sin*col[i1];
      col[i1] = b1; col[i2] = b2;
   }
   fAsOK = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Move in a's coord-system along axis-index ai.

void TEveTrans::Move(const TEveTrans& a, Int_t ai, Double_t amount)
{
   const Double_t* vec = a.fM + 4*--ai;
   fM[F03] += amount*vec[0];
   fM[F13] += amount*vec[1];
   fM[F23] += amount*vec[2];
}

////////////////////////////////////////////////////////////////////////////////
/// General move in a's coord-system.

void TEveTrans::Move3(const TEveTrans& a, Double_t x, Double_t y, Double_t z)
{
   const Double_t* m = a.fM;
   fM[F03] += x*m[F00] + y*m[F01] + z*m[F02];
   fM[F13] += x*m[F10] + y*m[F11] + z*m[F12];
   fM[F23] += x*m[F20] + y*m[F21] + z*m[F22];
}

////////////////////////////////////////////////////////////////////////////////
/// Rotate in a's coord-system, rotating base vector with index i1
/// into i2.

void TEveTrans::Rotate(const TEveTrans& a, Int_t i1, Int_t i2, Double_t amount)
{
   if(i1 == i2) return;
   TEveTrans x(a);
   x.Invert();
   MultLeft(x);
   RotatePF(i1, i2, amount);
   MultLeft(a);
   fAsOK = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set base-vector with index b.

void TEveTrans::SetBaseVec(Int_t b, Double_t x, Double_t y, Double_t z)
{
   Double_t* col = fM + 4*--b;
   col[0] = x; col[1] = y; col[2] = z;
   fAsOK = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set base-vector with index b.

void TEveTrans::SetBaseVec(Int_t b, const TVector3& v)
{
   Double_t* col = fM + 4*--b;
   v.GetXYZ(col);
   fAsOK = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Get base-vector with index b.

TVector3 TEveTrans::GetBaseVec(Int_t b) const
{
   return TVector3(&fM[4*--b]);
}

void TEveTrans::GetBaseVec(Int_t b, TVector3& v) const
{
   // Get base-vector with index b.

   const Double_t* col = fM + 4*--b;
   v.SetXYZ(col[0], col[1], col[2]);
}

////////////////////////////////////////////////////////////////////////////////
/// Set position (base-vec 4).

void TEveTrans::SetPos(Double_t x, Double_t y, Double_t z)
{
   fM[F03] = x; fM[F13] = y; fM[F23] = z;
}

void TEveTrans::SetPos(Double_t* x)
{
   // Set position (base-vec 4).
   fM[F03] = x[0]; fM[F13] = x[1]; fM[F23] = x[2];
}

void TEveTrans::SetPos(Float_t* x)
{
   // Set position (base-vec 4).
   fM[F03] = x[0]; fM[F13] = x[1]; fM[F23] = x[2];
}

void TEveTrans::SetPos(const TEveTrans& t)
{
   // Set position (base-vec 4).
   const Double_t* m = t.fM;
   fM[F03] = m[F03]; fM[F13] = m[F13]; fM[F23] = m[F23];
}

////////////////////////////////////////////////////////////////////////////////
/// Get position (base-vec 4).

void TEveTrans::GetPos(Double_t& x, Double_t& y, Double_t& z) const
{
   x = fM[F03]; y = fM[F13]; z = fM[F23];
}

void TEveTrans::GetPos(Double_t* x) const
{
   // Get position (base-vec 4).
   x[0] = fM[F03]; x[1] = fM[F13]; x[2] = fM[F23];
}

void TEveTrans::GetPos(Float_t* x) const
{
   // Get position (base-vec 4).
   x[0] = fM[F03]; x[1] = fM[F13]; x[2] = fM[F23];
}

void TEveTrans::GetPos(TVector3& v) const
{
   // Get position (base-vec 4).
   v.SetXYZ(fM[F03], fM[F13], fM[F23]);
}

TVector3 TEveTrans::GetPos() const
{
   // Get position (base-vec 4).
   return TVector3(fM[F03], fM[F13], fM[F23]);
}

namespace
{
inline void clamp_angle(Float_t& a)
{
   while(a < -TMath::TwoPi()) a += TMath::TwoPi();
   while(a >  TMath::TwoPi()) a -= TMath::TwoPi();
}
}

void TEveTrans::SetRotByAngles(Float_t a1, Float_t a2, Float_t a3)
{
   // Sets Rotation part as given by angles:
   // a1 around z, -a2 around y, a3 around x.

   clamp_angle(a1); clamp_angle(a2); clamp_angle(a3);

   Double_t a, b, c, d, e, f;
   a = TMath::Cos(a3); b = TMath::Sin(a3);
   c = TMath::Cos(a2); d = TMath::Sin(a2); // should be -sin(a2) for positive direction
   e = TMath::Cos(a1); f = TMath::Sin(a1);
   Double_t ad = a*d, bd = b*d;

   fM[F00] = c*e; fM[F01] = -bd*e - a*f; fM[F02] = -ad*e + b*f;
   fM[F10] = c*f; fM[F11] = -bd*f + a*e; fM[F12] = -ad*f - b*e;
   fM[F20] = d;   fM[F21] =  b*c;        fM[F22] =  a*c;

   fA1 = a1; fA2 = a2; fA3 = a3;
   fAsOK = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Sets Rotation part as given by angles a1, a1, a3 and pattern pat.
/// Pattern consists of "XxYyZz" characters.
/// eg: x means rotate about x axis, X means rotate in negative direction
/// xYz -> R_x(a3) * R_y(-a2) * R_z(a1); (standard Gled representation)
/// Note that angles and pattern elements have inverted order!
///
/// Implements Eulerian/Cardanian angles in a uniform way.

void TEveTrans::SetRotByAnyAngles(Float_t a1, Float_t a2, Float_t a3,
                                  const char* pat)
{
   int n = strspn(pat, "XxYyZz"); if(n > 3) n = 3;
   // Build Trans ... assign ...
   Float_t a[] = { a3, a2, a1 };
   UnitRot();
   for(int i=0; i<n; i++) {
      if(isupper(pat[i])) a[i] = -a[i];
      switch(pat[i]) {
         case 'x': case 'X': RotateLF(2, 3, a[i]); break;
         case 'y': case 'Y': RotateLF(3, 1, a[i]); break;
         case 'z': case 'Z': RotateLF(1, 2, a[i]); break;
      }
   }
   fAsOK = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Get Cardan rotation angles (pattern xYz above).

void TEveTrans::GetRotAngles(Float_t* x) const
{
   if(!fAsOK) {
      Double_t sx, sy, sz;
      GetScale(sx, sy, sz);
      Double_t d = fM[F20]/sx;
      if(d>1) d=1; else if(d<-1) d=-1; // Fix numerical errors
      fA2 = TMath::ASin(d);
      Double_t cos = TMath::Cos(fA2);
      if(TMath::Abs(cos) > 8.7e-6) {
         fA1 = TMath::ATan2(fM[F10], fM[F00]);
         fA3 = TMath::ATan2(fM[F21]/sy, fM[F22]/sz);
      } else {
         fA1 = TMath::ATan2(fM[F10]/sx, fM[F11]/sy);
         fA3 = 0;
      }
      fAsOK = kTRUE;
   }
   x[0] = fA1; x[1] = fA2; x[2] = fA3;
}

////////////////////////////////////////////////////////////////////////////////
/// Scale matrix. Translation part untouched.

void TEveTrans::Scale(Double_t sx, Double_t sy, Double_t sz)
{
   fM[F00] *= sx; fM[F10] *= sx; fM[F20] *= sx;
   fM[F01] *= sy; fM[F11] *= sy; fM[F21] *= sy;
   fM[F02] *= sz; fM[F12] *= sz; fM[F22] *= sz;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove scaling, make all base vectors of unit length.

Double_t TEveTrans::Unscale()
{
   Double_t sx, sy, sz;
   Unscale(sx, sy, sz);
   return (sx + sy + sz)/3;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove scaling, make all base vectors of unit length.

void TEveTrans::Unscale(Double_t& sx, Double_t& sy, Double_t& sz)
{
   GetScale(sx, sy, sz);
   fM[F00] /= sx; fM[F10] /= sx; fM[F20] /= sx;
   fM[F01] /= sy; fM[F11] /= sy; fM[F21] /= sy;
   fM[F02] /= sz; fM[F12] /= sz; fM[F22] /= sz;
}

////////////////////////////////////////////////////////////////////////////////
/// Deduce scales from sizes of base vectors.

void TEveTrans::GetScale(Double_t& sx, Double_t& sy, Double_t& sz) const
{
   sx = TMath::Sqrt( fM[F00]*fM[F00] + fM[F10]*fM[F10] + fM[F20]*fM[F20] );
   sy = TMath::Sqrt( fM[F01]*fM[F01] + fM[F11]*fM[F11] + fM[F21]*fM[F21] );
   sz = TMath::Sqrt( fM[F02]*fM[F02] + fM[F12]*fM[F12] + fM[F22]*fM[F22] );
}

////////////////////////////////////////////////////////////////////////////////
/// Set scaling.

void TEveTrans::SetScale(Double_t sx, Double_t sy, Double_t sz)
{
   sx /= TMath::Sqrt( fM[F00]*fM[F00] + fM[F10]*fM[F10] + fM[F20]*fM[F20] );
   sy /= TMath::Sqrt( fM[F01]*fM[F01] + fM[F11]*fM[F11] + fM[F21]*fM[F21] );
   sz /= TMath::Sqrt( fM[F02]*fM[F02] + fM[F12]*fM[F12] + fM[F22]*fM[F22] );

   fM[F00] *= sx; fM[F10] *= sx; fM[F20] *= sx;
   fM[F01] *= sy; fM[F11] *= sy; fM[F21] *= sy;
   fM[F02] *= sz; fM[F12] *= sz; fM[F22] *= sz;
}

////////////////////////////////////////////////////////////////////////////////
/// Change x scaling.

void TEveTrans::SetScaleX(Double_t sx)
{
   sx /= TMath::Sqrt( fM[F00]*fM[F00] + fM[F10]*fM[F10] + fM[F20]*fM[F20] );
   fM[F00] *= sx; fM[F10] *= sx; fM[F20] *= sx;
}

////////////////////////////////////////////////////////////////////////////////
/// Change y scaling.

void TEveTrans::SetScaleY(Double_t sy)
{
   sy /= TMath::Sqrt( fM[F01]*fM[F01] + fM[F11]*fM[F11] + fM[F21]*fM[F21] );
   fM[F01] *= sy; fM[F11] *= sy; fM[F21] *= sy;
}

////////////////////////////////////////////////////////////////////////////////
/// Change z scaling.

void TEveTrans::SetScaleZ(Double_t sz)
{
   sz /= TMath::Sqrt( fM[F02]*fM[F02] + fM[F12]*fM[F12] + fM[F22]*fM[F22] );
   fM[F02] *= sz; fM[F12] *= sz; fM[F22] *= sz;
}

////////////////////////////////////////////////////////////////////////////////
/// Multiply vector in-place.

void TEveTrans::MultiplyIP(TVector3& v, Double_t w) const
{
   v.SetXYZ(fM[F00]*v.x() + fM[F01]*v.y() + fM[F02]*v.z() + fM[F03]*w,
            fM[F10]*v.x() + fM[F11]*v.y() + fM[F12]*v.z() + fM[F13]*w,
            fM[F20]*v.x() + fM[F21]*v.y() + fM[F22]*v.z() + fM[F23]*w);
}

////////////////////////////////////////////////////////////////////////////////
/// Multiply vector in-place.

void TEveTrans::MultiplyIP(Double_t* v, Double_t w) const
{
   Double_t r[3] = { v[0], v[1], v[2] };
   v[0] = fM[F00]*r[0] + fM[F01]*r[1] + fM[F02]*r[2] + fM[F03]*w;
   v[1] = fM[F10]*r[0] + fM[F11]*r[1] + fM[F12]*r[2] + fM[F13]*w;
   v[2] = fM[F20]*r[0] + fM[F21]*r[1] + fM[F22]*r[2] + fM[F23]*w;
}

////////////////////////////////////////////////////////////////////////////////
/// Multiply vector in-place.

void TEveTrans::MultiplyIP(Float_t* v, Double_t w) const
{
   Double_t r[3] = { v[0], v[1], v[2] };
   v[0] = fM[F00]*r[0] + fM[F01]*r[1] + fM[F02]*r[2] + fM[F03]*w;
   v[1] = fM[F10]*r[0] + fM[F11]*r[1] + fM[F12]*r[2] + fM[F13]*w;
   v[2] = fM[F20]*r[0] + fM[F21]*r[1] + fM[F22]*r[2] + fM[F23]*w;
}

////////////////////////////////////////////////////////////////////////////////
/// Multiply vector and return it.

TVector3 TEveTrans::Multiply(const TVector3& v, Double_t w) const
{
   return TVector3(fM[F00]*v.x() + fM[F01]*v.y() + fM[F02]*v.z() + fM[F03]*w,
                   fM[F10]*v.x() + fM[F11]*v.y() + fM[F12]*v.z() + fM[F13]*w,
                   fM[F20]*v.x() + fM[F21]*v.y() + fM[F22]*v.z() + fM[F23]*w);
}

////////////////////////////////////////////////////////////////////////////////
/// Multiply vector and fill output array vout.

void TEveTrans::Multiply(const Double_t *vin, Double_t* vout, Double_t w) const
{
   vout[0] = fM[F00]*vin[0] + fM[F01]*vin[1] + fM[F02]*vin[2] + fM[F03]*w;
   vout[1] = fM[F10]*vin[0] + fM[F11]*vin[1] + fM[F12]*vin[2] + fM[F13]*w;
   vout[2] = fM[F20]*vin[0] + fM[F21]*vin[1] + fM[F22]*vin[2] + fM[F23]*w;
}

////////////////////////////////////////////////////////////////////////////////
/// Rotate vector in-place. Translation is NOT applied.

void TEveTrans::RotateIP(TVector3& v) const
{
   v.SetXYZ(fM[F00]*v.x() + fM[F01]*v.y() + fM[F02]*v.z(),
            fM[F10]*v.x() + fM[F11]*v.y() + fM[F12]*v.z(),
            fM[F20]*v.x() + fM[F21]*v.y() + fM[F22]*v.z());
}

////////////////////////////////////////////////////////////////////////////////
/// Rotate vector in-place. Translation is NOT applied.

void TEveTrans::RotateIP(Double_t* v) const
{
   Double_t t[3] = { v[0], v[1], v[2] };

   v[0] = fM[F00]*t[0] + fM[F01]*t[1] + fM[F02]*t[2];
   v[1] = fM[F10]*t[0] + fM[F11]*t[1] + fM[F12]*t[2];
   v[2] = fM[F20]*t[0] + fM[F21]*t[1] + fM[F22]*t[2];
}

////////////////////////////////////////////////////////////////////////////////
/// Rotate vector in-place. Translation is NOT applied.

void TEveTrans::RotateIP(Float_t* v) const
{
   Double_t t[3] = { v[0], v[1], v[2] };

   v[0] = fM[F00]*t[0] + fM[F01]*t[1] + fM[F02]*t[2];
   v[1] = fM[F10]*t[0] + fM[F11]*t[1] + fM[F12]*t[2];
   v[2] = fM[F20]*t[0] + fM[F21]*t[1] + fM[F22]*t[2];
}

////////////////////////////////////////////////////////////////////////////////
/// Rotate vector and return the rotated vector. Translation is NOT applied.

TVector3 TEveTrans::Rotate(const TVector3& v) const
{
   return TVector3(fM[F00]*v.x() + fM[F01]*v.y() + fM[F02]*v.z(),
                   fM[F10]*v.x() + fM[F11]*v.y() + fM[F12]*v.z(),
                   fM[F20]*v.x() + fM[F21]*v.y() + fM[F22]*v.z());
}

////////////////////////////////////////////////////////////////////////////////
/// Norm 3-vector in column col.

Double_t TEveTrans::Norm3Column(Int_t col)
{
   Double_t* c = fM + 4*--col;
   const Double_t  l = TMath::Sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
   c[0] /= l; c[1] /= l; c[2] /= l;
   return l;
}

////////////////////////////////////////////////////////////////////////////////
/// Orto-norm 3-vector in column col with respect to column ref.

Double_t TEveTrans::Orto3Column(Int_t col, Int_t ref)
{
   Double_t* c =  fM + 4*--col;
   Double_t* rc = fM + 4*--ref;
   const Double_t dp = c[0]*rc[0] + c[1]*rc[1] + c[2]*rc[2];
   c[0] -= rc[0]*dp; c[1] -= rc[1]*dp; c[2] -= rc[2]*dp;
   return dp;
}

////////////////////////////////////////////////////////////////////////////////
/// Orto-norm columns 1 to 3.

void TEveTrans::OrtoNorm3()
{
   Norm3Column(1);
   Orto3Column(2,1); Norm3Column(2);
   fM[F02] = fM[F10]*fM[F21] - fM[F11]*fM[F20];
   fM[F12] = fM[F20]*fM[F01] - fM[F21]*fM[F00];
   fM[F22] = fM[F00]*fM[F11] - fM[F01]*fM[F10];
   // Cross-product faster than the following.
   // Orto3Column(3,1); Orto3Column(3,2); Norm3Column(3);
}

////////////////////////////////////////////////////////////////////////////////
/// Invert matrix.
/// Copied from ROOT's TMatrixFCramerInv.

Double_t TEveTrans::Invert()
{
   static const TEveException eh("TEveTrans::Invert ");

   // Find all NECESSARY 2x2 dets:  (18 of them)
   const Double_t det2_12_01 = fM[F10]*fM[F21] - fM[F11]*fM[F20];
   const Double_t det2_12_02 = fM[F10]*fM[F22] - fM[F12]*fM[F20];
   const Double_t det2_12_03 = fM[F10]*fM[F23] - fM[F13]*fM[F20];
   const Double_t det2_12_13 = fM[F11]*fM[F23] - fM[F13]*fM[F21];
   const Double_t det2_12_23 = fM[F12]*fM[F23] - fM[F13]*fM[F22];
   const Double_t det2_12_12 = fM[F11]*fM[F22] - fM[F12]*fM[F21];
   const Double_t det2_13_01 = fM[F10]*fM[F31] - fM[F11]*fM[F30];
   const Double_t det2_13_02 = fM[F10]*fM[F32] - fM[F12]*fM[F30];
   const Double_t det2_13_03 = fM[F10]*fM[F33] - fM[F13]*fM[F30];
   const Double_t det2_13_12 = fM[F11]*fM[F32] - fM[F12]*fM[F31];
   const Double_t det2_13_13 = fM[F11]*fM[F33] - fM[F13]*fM[F31];
   const Double_t det2_13_23 = fM[F12]*fM[F33] - fM[F13]*fM[F32];
   const Double_t det2_23_01 = fM[F20]*fM[F31] - fM[F21]*fM[F30];
   const Double_t det2_23_02 = fM[F20]*fM[F32] - fM[F22]*fM[F30];
   const Double_t det2_23_03 = fM[F20]*fM[F33] - fM[F23]*fM[F30];
   const Double_t det2_23_12 = fM[F21]*fM[F32] - fM[F22]*fM[F31];
   const Double_t det2_23_13 = fM[F21]*fM[F33] - fM[F23]*fM[F31];
   const Double_t det2_23_23 = fM[F22]*fM[F33] - fM[F23]*fM[F32];

   // Find all NECESSARY 3x3 dets:   (16 of them)
   const Double_t det3_012_012 = fM[F00]*det2_12_12 - fM[F01]*det2_12_02 + fM[F02]*det2_12_01;
   const Double_t det3_012_013 = fM[F00]*det2_12_13 - fM[F01]*det2_12_03 + fM[F03]*det2_12_01;
   const Double_t det3_012_023 = fM[F00]*det2_12_23 - fM[F02]*det2_12_03 + fM[F03]*det2_12_02;
   const Double_t det3_012_123 = fM[F01]*det2_12_23 - fM[F02]*det2_12_13 + fM[F03]*det2_12_12;
   const Double_t det3_013_012 = fM[F00]*det2_13_12 - fM[F01]*det2_13_02 + fM[F02]*det2_13_01;
   const Double_t det3_013_013 = fM[F00]*det2_13_13 - fM[F01]*det2_13_03 + fM[F03]*det2_13_01;
   const Double_t det3_013_023 = fM[F00]*det2_13_23 - fM[F02]*det2_13_03 + fM[F03]*det2_13_02;
   const Double_t det3_013_123 = fM[F01]*det2_13_23 - fM[F02]*det2_13_13 + fM[F03]*det2_13_12;
   const Double_t det3_023_012 = fM[F00]*det2_23_12 - fM[F01]*det2_23_02 + fM[F02]*det2_23_01;
   const Double_t det3_023_013 = fM[F00]*det2_23_13 - fM[F01]*det2_23_03 + fM[F03]*det2_23_01;
   const Double_t det3_023_023 = fM[F00]*det2_23_23 - fM[F02]*det2_23_03 + fM[F03]*det2_23_02;
   const Double_t det3_023_123 = fM[F01]*det2_23_23 - fM[F02]*det2_23_13 + fM[F03]*det2_23_12;
   const Double_t det3_123_012 = fM[F10]*det2_23_12 - fM[F11]*det2_23_02 + fM[F12]*det2_23_01;
   const Double_t det3_123_013 = fM[F10]*det2_23_13 - fM[F11]*det2_23_03 + fM[F13]*det2_23_01;
   const Double_t det3_123_023 = fM[F10]*det2_23_23 - fM[F12]*det2_23_03 + fM[F13]*det2_23_02;
   const Double_t det3_123_123 = fM[F11]*det2_23_23 - fM[F12]*det2_23_13 + fM[F13]*det2_23_12;

   // Find the 4x4 det:
   const Double_t det = fM[F00]*det3_123_123 - fM[F01]*det3_123_023 +
      fM[F02]*det3_123_013 - fM[F03]*det3_123_012;

   if(det == 0) {
      throw(eh + "matrix is singular.");
   }

   const Double_t oneOverDet = 1.0/det;
   const Double_t mn1OverDet = - oneOverDet;

   fM[F00] = det3_123_123 * oneOverDet;
   fM[F01] = det3_023_123 * mn1OverDet;
   fM[F02] = det3_013_123 * oneOverDet;
   fM[F03] = det3_012_123 * mn1OverDet;

   fM[F10] = det3_123_023 * mn1OverDet;
   fM[F11] = det3_023_023 * oneOverDet;
   fM[F12] = det3_013_023 * mn1OverDet;
   fM[F13] = det3_012_023 * oneOverDet;

   fM[F20] = det3_123_013 * oneOverDet;
   fM[F21] = det3_023_013 * mn1OverDet;
   fM[F22] = det3_013_013 * oneOverDet;
   fM[F23] = det3_012_013 * mn1OverDet;

   fM[F30] = det3_123_012 * mn1OverDet;
   fM[F31] = det3_023_012 * oneOverDet;
   fM[F32] = det3_013_012 * mn1OverDet;
   fM[F33] = det3_012_012 * oneOverDet;

   fAsOK = kFALSE;
   return det;
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class TEveTrans.

void TEveTrans::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      TEveTrans::Class()->ReadBuffer(R__b, this);
      fAsOK = kFALSE;
   } else {
      TEveTrans::Class()->WriteBuffer(R__b, this);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Print in reasonable format.

void TEveTrans::Print(Option_t* /*option*/) const
{
   const Double_t* row = fM;
   for(Int_t i=0; i<4; ++i, ++row)
      printf("%8.3f %8.3f %8.3f | %8.3f\n", row[0], row[4], row[8], row[12]);
}

#include <iomanip>

////////////////////////////////////////////////////////////////////////////////
/// Print to std::ostream.

std::ostream& operator<<(std::ostream& s, const TEveTrans& t)
{
   s.setf(std::ios::fixed, std::ios::floatfield);
   s.precision(3);
   for(Int_t i=1; i<=4; i++)
      for(Int_t j=1; j<=4; j++)
         s << t(i,j) << ((j==4) ? "\n" : "\t");
   return s;
}

#include "TGeoMatrix.h"
#include "TBuffer3D.h"

void TEveTrans::SetFrom(Double_t* carr)
{
   // Initialize from array.

   fUseTrans = kTRUE;
   memcpy(fM, carr, 16*sizeof(Double_t));
   fAsOK = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize from TGeoMatrix.

void TEveTrans::SetFrom(const TGeoMatrix& mat)
{
   fUseTrans = kTRUE;
   const Double_t *r = mat.GetRotationMatrix();
   const Double_t *t = mat.GetTranslation();
   Double_t       *m = fM;
   if (mat.IsScale())
   {
      const Double_t *s = mat.GetScale();
      m[0]  = r[0]*s[0]; m[1]  = r[3]*s[0]; m[2]  = r[6]*s[0]; m[3]  = 0;
      m[4]  = r[1]*s[1]; m[5]  = r[4]*s[1]; m[6]  = r[7]*s[1]; m[7]  = 0;
      m[8]  = r[2]*s[2]; m[9]  = r[5]*s[2]; m[10] = r[8]*s[2]; m[11] = 0;
      m[12] = t[0];      m[13] = t[1];      m[14] = t[2];      m[15] = 1;
   }
   else
   {
      m[0]  = r[0];      m[1]  = r[3];      m[2]  = r[6];      m[3]  = 0;
      m[4]  = r[1];      m[5]  = r[4];      m[6]  = r[7];      m[7]  = 0;
      m[8]  = r[2];      m[9]  = r[5];      m[10] = r[8];      m[11] = 0;
      m[12] = t[0];      m[13] = t[1];      m[14] = t[2];      m[15] = 1;
   }
   fAsOK = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set TGeoHMatrix mat.

void TEveTrans::SetGeoHMatrix(TGeoHMatrix& mat)
{
   Double_t *r = mat.GetRotationMatrix();
   Double_t *t = mat.GetTranslation();
   Double_t *s = mat.GetScale();
   if (fUseTrans)
   {
      mat.SetBit(TGeoMatrix::kGeoGenTrans);
      Double_t *m = fM;
      GetScale(s[0], s[1], s[2]);
      r[0] = m[0]/s[0]; r[3] = m[1]/s[0]; r[6] = m[2]/s[0]; m += 4;
      r[1] = m[0]/s[1]; r[4] = m[1]/s[1]; r[7] = m[2]/s[1]; m += 4;
      r[2] = m[0]/s[2]; r[5] = m[1]/s[2]; r[8] = m[2]/s[2]; m += 4;
      t[0] = m[0];      t[1] = m[1];      t[2] = m[2];
   }
   else
   {
      mat.ResetBit(TGeoMatrix::kGeoGenTrans);
      r[0] = 1; r[3] = 0; r[6] = 0;
      r[1] = 0; r[4] = 1; r[7] = 0;
      r[2] = 0; r[5] = 0; r[8] = 1;
      s[0] = s[1] = s[2] = 1;
      t[0] = t[1] = t[2] = 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Fill transformation part TBuffer3D core section.

void TEveTrans::SetBuffer3D(TBuffer3D& buff)
{
   buff.fLocalFrame = fUseTrans;
   if (fUseTrans) {
      // In phys-shape ctor the rotation part is transposed, due to
      // TGeo's convention for rotation matrix. So we have to transpose
      // it here, also.
      Double_t *m = buff.fLocalMaster;
      m[0]  = fM[0];  m[1]  = fM[4];  m[2]  = fM[8];  m[3]  = fM[3];
      m[4]  = fM[1];  m[5]  = fM[5];  m[6]  = fM[9];  m[7]  = fM[7];
      m[8]  = fM[2];  m[9]  = fM[6];  m[10] = fM[10]; m[11] = fM[11];
      m[12] = fM[12]; m[13] = fM[13]; m[14] = fM[14]; m[15] = fM[15];
      // Otherwise this would do:
      // memcpy(buff.fLocalMaster, fM, 16*sizeof(Double_t));
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Test if the transformation is a scale.
/// To be used by ROOT TGLObject descendants that potentially need to
/// use GL_NORMALIZE.
/// The low/high limits are expected to be squares of actual limits.
///
/// Ideally this should be done by the TGLViewer [but is not].

Bool_t TEveTrans::IsScale(Double_t low, Double_t high) const
{
   if (!fUseTrans) return kFALSE;
   Double_t s;
   s = fM[F00]*fM[F00] + fM[F10]*fM[F10] + fM[F20]*fM[F20];
   if (s < low || s > high) return kTRUE;
   s = fM[F01]*fM[F01] + fM[F11]*fM[F11] + fM[F21]*fM[F21];
   if (s < low || s > high) return kTRUE;
   s = fM[F02]*fM[F02] + fM[F12]*fM[F12] + fM[F22]*fM[F22];
   if (s < low || s > high) return kTRUE;
   return kFALSE;
}
