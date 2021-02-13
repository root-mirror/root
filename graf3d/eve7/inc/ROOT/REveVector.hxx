// @(#)root/eve7:$Id$
// Author: Matevz Tadel 2007, 2018

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_REveVector
#define ROOT7_REveVector

#include "TMath.h"
#include <cstddef>

class TVector3;

namespace ROOT {
namespace Experimental {

////////////////////////////////////////////////////////////////////////////////
/// REveVectorT
/// A three-vector template without TObject inheritance and virtual functions.
////////////////////////////////////////////////////////////////////////////////

template <typename TT>
class REveVectorT {
public:
   TT fX{0}, fY{0}, fZ{0}; // Components of the vector.

   REveVectorT() = default;
   template <typename OO>
   REveVectorT(const REveVectorT<OO>& v) : fX(v.fX), fY(v.fY), fZ(v.fZ) {}
   REveVectorT(const Float_t*  v) : fX(v[0]), fY(v[1]), fZ(v[2]) {}
   REveVectorT(const Double_t* v) : fX(v[0]), fY(v[1]), fZ(v[2]) {}
   REveVectorT(TT x, TT y, TT  z) : fX(x), fY(y), fZ(z) {}

   void Dump() const;

#ifdef R__WIN32
   const TT *Arr() const
   {
      if (offsetof(REveVectorT, fZ) != offsetof(REveVectorT, fX) + 2 * sizeof(TT))
         Error("REveVectorT", "Subsequent members cannot be accessed as array!");
      return &fX;
   }
   TT *Arr()
   {
      if (offsetof(REveVectorT, fZ) != offsetof(REveVectorT, fX) + 2 * sizeof(TT))
         Error("REveVectorT", "Subsequent members cannot be accessed as array!");
      return &fX;
   }
#else
   const TT *Arr() const
   {
      static_assert(offsetof(REveVectorT, fZ) == offsetof(REveVectorT, fX) + 2 * sizeof(TT),
                    "Subsequent members cannot be accessed as array!");
      return &fX;
   }
   TT *Arr()
   {
      static_assert(offsetof(REveVectorT, fZ) == offsetof(REveVectorT, fX) + 2 * sizeof(TT),
                    "Subsequent members cannot be accessed as array!");
      return &fX;
   }
#endif

   operator const TT*() const { return Arr(); }
   operator       TT*()       { return Arr(); }

   TT  operator [] (Int_t idx) const { return Arr()[idx]; }
   TT& operator [] (Int_t idx)       { return Arr()[idx]; }

   REveVectorT& operator*=(TT s)                 { fX *= s;    fY *= s;    fZ *= s;    return *this; }
   REveVectorT& operator+=(const REveVectorT& v) { fX += v.fX; fY += v.fY; fZ += v.fZ; return *this; }
   REveVectorT& operator-=(const REveVectorT& v) { fX -= v.fX; fY -= v.fY; fZ -= v.fZ; return *this; }

   void Set(const Float_t*  v) { fX = v[0]; fY = v[1]; fZ = v[2]; }
   void Set(const Double_t* v) { fX = v[0]; fY = v[1]; fZ = v[2]; }
   void Set(TT x, TT  y, TT z) { fX = x; fY = y; fZ = z; }
   void Set(const TVector3& v);

   template <typename OO>
   void Set(const REveVectorT<OO>& v) { fX = v.fX;  fY = v.fY;  fZ = v.fZ; }

   void NegateXYZ() { fX = - fX; fY = -fY; fZ = -fZ; }
   TT   Normalize(TT length=1);

   TT   Phi()      const;
   TT   Theta()    const;
   TT   CosTheta() const;
   TT   Eta()      const;

   TT   Mag2()  const { return fX*fX + fY*fY + fZ*fZ; }
   TT   Mag()   const { return TMath::Sqrt(Mag2());   }

   TT   Perp2() const { return fX*fX + fY*fY;        }
   TT   Perp()  const { return TMath::Sqrt(Perp2()); }
   TT   R()     const { return Perp();               }

   TT   Distance(const REveVectorT& v) const;
   TT   SquareDistance(const REveVectorT& v) const;

   TT   Dot(const REveVectorT& a) const;

   REveVectorT  Cross(const REveVectorT& a) const;

   REveVectorT& Sub(const REveVectorT& a, const REveVectorT& b);
   REveVectorT& Mult(const REveVectorT& a, TT af);

   REveVectorT  Orthogonal() const;
   void         OrthoNormBase(REveVectorT& a, REveVectorT& b) const;

   Bool_t       IsZero() const { return fX == 0 && fY == 0 && fZ == 0; }
};

typedef REveVectorT<Float_t>  REveVector;
typedef REveVectorT<Float_t>  REveVectorF;
typedef REveVectorT<Double_t> REveVectorD;

//______________________________________________________________________________
template<typename TT>
inline TT REveVectorT<TT>::Phi() const
{
   return fX == 0 && fY == 0 ? 0 : TMath::ATan2(fY, fX);
}

//______________________________________________________________________________
template<typename TT>
inline TT REveVectorT<TT>::Theta() const
{
   return fX == 0 && fY == 0 && fZ == 0 ? 0 : TMath::ATan2(Perp(), fZ);
}

//______________________________________________________________________________
template<typename TT>
inline TT REveVectorT<TT>::CosTheta() const
{
   Float_t ptot = Mag(); return ptot == 0 ? 1 : fZ/ptot;
}

//______________________________________________________________________________
template<typename TT>
inline TT REveVectorT<TT>::Distance(const REveVectorT& b) const
{
   return TMath::Sqrt((fX - b.fX)*(fX - b.fX) +
                      (fY - b.fY)*(fY - b.fY) +
                      (fZ - b.fZ)*(fZ - b.fZ));
}

//______________________________________________________________________________
template<typename TT>
inline TT REveVectorT<TT>::SquareDistance(const REveVectorT& b) const
{
   return ((fX - b.fX) * (fX - b.fX) +
           (fY - b.fY) * (fY - b.fY) +
           (fZ - b.fZ) * (fZ - b.fZ));
}

//______________________________________________________________________________
template<typename TT>
inline TT REveVectorT<TT>::Dot(const REveVectorT& a) const
{
   return a.fX*fX + a.fY*fY + a.fZ*fZ;
}

//______________________________________________________________________________
template<typename TT>
inline REveVectorT<TT> REveVectorT<TT>::Cross(const REveVectorT<TT>& a) const
{
   REveVectorT<TT> r;
   r.fX = fY * a.fZ - fZ * a.fY;
   r.fY = fZ * a.fX - fX * a.fZ;
   r.fZ = fX * a.fY - fY * a.fX;
   return r;
}

//______________________________________________________________________________
template<typename TT>
inline REveVectorT<TT>& REveVectorT<TT>::Sub(const REveVectorT<TT>& a, const REveVectorT<TT>& b)
{
   fX = a.fX - b.fX;
   fY = a.fY - b.fY;
   fZ = a.fZ - b.fZ;
   return *this;
}

//______________________________________________________________________________
template<typename TT>
inline REveVectorT<TT>& REveVectorT<TT>::Mult(const REveVectorT<TT>& a, TT af)
{
   fX = a.fX * af;
   fY = a.fY * af;
   fZ = a.fZ * af;
   return *this;
}

//______________________________________________________________________________
template<typename TT>
inline REveVectorT<TT> operator+(const REveVectorT<TT>& a, const REveVectorT<TT>& b)
{
   REveVectorT<TT> r(a);
   return r += b;
}

//______________________________________________________________________________
template<typename TT>
inline REveVectorT<TT> operator-(const REveVectorT<TT>& a, const REveVectorT<TT>& b)
{
   REveVectorT<TT> r(a);
   return r -= b;
}

//______________________________________________________________________________
template<typename TT>
inline REveVectorT<TT> operator*(const REveVectorT<TT>& a, TT b)
{
   REveVectorT<TT> r(a);
   return r *= b;
}

//______________________________________________________________________________
template<typename TT>
inline REveVectorT<TT> operator*(TT b, const REveVectorT<TT>& a)
{
   REveVectorT<TT> r(a);
   return r *= b;
}

////////////////////////////////////////////////////////////////////////////////
/// REveVector4T
/// A four-vector template without TObject inheritance and virtual functions.
////////////////////////////////////////////////////////////////////////////////

template <typename TT>
class REveVector4T : public REveVectorT<TT>
{
   typedef REveVectorT<TT> TP;

public:
   TT fT;

   REveVector4T() : TP(),  fT(0) {}
   template <typename OO>
   REveVector4T(const REveVectorT<OO>& v) : TP(v.fX, v.fY, v.fZ), fT(0) {}
   template <typename OO>
   REveVector4T(const REveVectorT<OO>& v, Float_t t) : TP(v.fX, v.fY, v.fZ), fT(t) {}
   template <typename OO>
   REveVector4T(const REveVector4T<OO>& v) : TP(v.fX, v.fY, v.fZ), fT(v.fT) {}
   REveVector4T(const Float_t*  v) : TP(v), fT(v[3]) {}
   REveVector4T(const Double_t* v) : TP(v), fT(v[3]) {}
   REveVector4T(TT x, TT y, TT z, TT t=0) : TP(x, y, z), fT(t) {}

   void Dump() const;

   REveVector4T& operator*=(TT s)                  { TP::operator*=(s); fT *= s;    return *this; }
   REveVector4T& operator+=(const REveVector4T& v) { TP::operator+=(v); fT += v.fT; return *this; }
   REveVector4T& operator-=(const REveVector4T& v) { TP::operator-=(v); fT -= v.fT; return *this; }

   using TP::operator+=;
   using TP::operator-=;
};

typedef REveVector4T<Float_t>  REveVector4;
typedef REveVector4T<Float_t>  REveVector4F;
typedef REveVector4T<Double_t> REveVector4D;

//______________________________________________________________________________
template<typename TT>
inline REveVector4T<TT> operator+(const REveVector4T<TT>& a, const REveVector4T<TT>& b)
{
   return REveVector4T<TT>(a.fX + b.fX, a.fY + b.fY, a.fZ + b.fZ, a.fT + b.fT);
}

//______________________________________________________________________________
template<typename TT>
inline REveVector4T<TT> operator-(const REveVector4T<TT>& a, const REveVector4T<TT>& b)
{
   return REveVector4T<TT>(a.fX - b.fX, a.fY - b.fY, a.fZ - b.fZ, a.fT - b.fT);
}

//______________________________________________________________________________
template<typename TT>
inline REveVector4T<TT> operator*(const REveVector4T<TT>& a, TT b)
{
   return REveVector4T<TT>(a.fX*b, a.fY*b, a.fZ*b, a.fT*b);
}

//______________________________________________________________________________
template<typename TT>
inline REveVector4T<TT> operator*(TT b, const REveVector4T<TT>& a)
{
   return REveVector4T<TT>(a.fX*b, a.fY*b, a.fZ*b, a.fT*b);
}

////////////////////////////////////////////////////////////////////////////////
/// REveVector2T
/// A two-vector template without TObject inheritance and virtual functions.
////////////////////////////////////////////////////////////////////////////////

template <typename TT>
class REveVector2T
{
public:
   TT fX, fY; // Components of the point.

   REveVector2T() : fX(0), fY(0) {}
   template <typename OO>
   REveVector2T(const REveVector2T<OO>& v) : fX(v.fX), fY(v.fY) {}
   REveVector2T(const Float_t* v)  : fX(v[0]), fY(v[1]) {}
   REveVector2T(const Double_t* v) : fX(v[0]), fY(v[1]) {}
   REveVector2T(TT x, TT y) : fX(x), fY(y)    {}

   void Dump() const;

   operator const TT*() const { return &fX; }
   operator       TT*()       { return &fX; }

   REveVector2T& operator*=(TT s)                  { fX *= s;    fY *= s;    return *this; }
   REveVector2T& operator+=(const REveVector2T& v) { fX += v.fX; fY += v.fY; return *this; }
   REveVector2T& operator-=(const REveVector2T& v) { fX -= v.fX; fY -= v.fY; return *this; }

   TT& operator[](Int_t idx)       { return (&fX)[idx]; }
   TT  operator[](Int_t idx) const { return (&fX)[idx]; }

   const TT* Arr() const { return &fX; }
   TT* Arr()             { return &fX; }

   void Set(const Float_t*  v) { fX = v[0]; fY = v[1]; }
   void Set(const Double_t* v) { fX = v[0]; fY = v[1]; }
   void Set(TT x, TT y) { fX = x; fY = y; }

   template <typename OO>
   void Set(const REveVector2T<OO>& v) { fX = v.fX; fY = v.fY; }

   void NegateXY() { fX = - fX; fY = -fY; }
   void Normalize(TT length=1);

   TT Phi()  const;

   TT Mag2() const { return fX*fX + fY*fY;}
   TT Mag()  const { return TMath::Sqrt(Mag2());}

   TT Distance(const REveVector2T& v) const;
   TT SquareDistance(const REveVector2T& v) const;

   TT Dot(const REveVector2T& a) const;
   TT Cross(const REveVector2T& a) const;

   REveVector2T& Sub(const REveVector2T& p, const REveVector2T& q);

   REveVector2T& Mult(const REveVector2T& a, TT af);
};

typedef REveVector2T<Float_t>  REveVector2;
typedef REveVector2T<Float_t>  REveVector2F;
typedef REveVector2T<Double_t> REveVector2D;

//______________________________________________________________________________
template<typename TT>
inline TT REveVector2T<TT>::Phi() const
{
   return fX == 0.0 && fY == 0.0 ? 0.0 : TMath::ATan2(fY, fX);
}

//______________________________________________________________________________
template<typename TT>
inline TT REveVector2T<TT>::Distance( const REveVector2T<TT>& b) const
{
   return TMath::Sqrt((fX - b.fX)*(fX - b.fX) +
                      (fY - b.fY)*(fY - b.fY));
}

//______________________________________________________________________________
template<typename TT>
inline TT REveVector2T<TT>::SquareDistance(const REveVector2T<TT>& b) const
{
   return ((fX - b.fX) * (fX - b.fX) +
           (fY - b.fY) * (fY - b.fY));
}

//______________________________________________________________________________
template<typename TT>
inline TT REveVector2T<TT>::Dot(const REveVector2T<TT>& a) const
{
   return a.fX*fX + a.fY*fY;
}

//______________________________________________________________________________
template<typename TT>
inline TT REveVector2T<TT>::Cross(const REveVector2T<TT>& a) const
{
   return fX * a.fY - fY * a.fX;
}

//______________________________________________________________________________
template<typename TT>
inline REveVector2T<TT>& REveVector2T<TT>::Sub(const REveVector2T<TT>& p, const REveVector2T<TT>& q)
{
   fX = p.fX - q.fX;
   fY = p.fY - q.fY;
   return *this;
}

//______________________________________________________________________________
template<typename TT>
inline REveVector2T<TT>& REveVector2T<TT>::Mult(const REveVector2T<TT>& a, TT af)
{
   fX = a.fX * af;
   fY = a.fY * af;
   return *this;
}

//______________________________________________________________________________
template<typename TT>
inline REveVector2T<TT> operator+(const REveVector2T<TT>& a, const REveVector2T<TT>& b)
{
   REveVector2T<TT> r(a);
   return r += b;
}

//______________________________________________________________________________
template<typename TT>
inline REveVector2T<TT> operator-(const REveVector2T<TT>& a, const REveVector2T<TT>& b)
{
   REveVector2T<TT> r(a);
   return r -= b;
}

//______________________________________________________________________________
template<typename TT>
inline REveVector2T<TT> operator*(const REveVector2T<TT>& a, TT b)
{
   REveVector2T<TT> r(a);
   return r *= b;
}

//______________________________________________________________________________
template<typename TT>
inline REveVector2T<TT> operator*(TT b, const REveVector2T<TT>& a)
{
   REveVector2T<TT> r(a);
   return r *= b;
}

}}

#endif
