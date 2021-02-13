// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TEveVector
#define ROOT_TEveVector

#include "Rtypes.h"
#include "TMath.h"
#include <cstddef>

class TVector3;


//==============================================================================
// TEveVectorT
//==============================================================================

template <typename TT>
class TEveVectorT
{
public:
   TT fX, fY, fZ; // Components of the vector.

   TEveVectorT() : fX(0), fY(0), fZ(0) {}
   template <typename OO>
   TEveVectorT(const TEveVectorT<OO>& v) : fX(v.fX), fY(v.fY), fZ(v.fZ) {}
   TEveVectorT(const Float_t*  v) : fX(v[0]), fY(v[1]), fZ(v[2]) {}
   TEveVectorT(const Double_t* v) : fX(v[0]), fY(v[1]), fZ(v[2]) {}
   TEveVectorT(TT x, TT y, TT  z) : fX(x), fY(y), fZ(z) {}

   void Dump() const;

#ifdef R__WIN32
   // This fixes the following rootcling error when generating the dictionary:
   // error G34C21FBE: static_assert expression is not an integral constant expression
   // FIXME: check if the error is fixed when upgrading llvm/clang 
   const TT *Arr() const
   {
      if (offsetof(TEveVectorT, fZ) != offsetof(TEveVectorT, fX) + 2 * sizeof(TT))
         Error("TEveVectorT", "Subsequent nembers cannot be accessed as array!");
      return &fX;
   }
   TT *Arr()
   {
      if (offsetof(TEveVectorT, fZ) != offsetof(TEveVectorT, fX) + 2 * sizeof(TT))
         Error("TEveVectorT", "Subsequent nembers cannot be accessed as array!");
      return &fX;
   }
#else
   const TT *Arr() const
   {
      static_assert(offsetof(TEveVectorT, fZ) == offsetof(TEveVectorT, fX) + 2 * sizeof(TT),
                    "Subsequent nembers cannot be accessed as array!");
      return &fX;
   }
   TT *Arr()
   {
      static_assert(offsetof(TEveVectorT, fZ) == offsetof(TEveVectorT, fX) + 2 * sizeof(TT),
                    "Subsequent nembers cannot be accessed as array!");
      return &fX;
   }
#endif

   operator const TT*() const { return Arr(); }
   operator       TT*()       { return Arr(); }

   TT  operator [] (Int_t idx) const { return Arr()[idx]; }
   TT& operator [] (Int_t idx)       { return Arr()[idx]; }

   TEveVectorT& operator*=(TT s)                 { fX *= s;    fY *= s;    fZ *= s;    return *this; }
   TEveVectorT& operator+=(const TEveVectorT& v) { fX += v.fX; fY += v.fY; fZ += v.fZ; return *this; }
   TEveVectorT& operator-=(const TEveVectorT& v) { fX -= v.fX; fY -= v.fY; fZ -= v.fZ; return *this; }

   void Set(const Float_t*  v) { fX = v[0]; fY = v[1]; fZ = v[2]; }
   void Set(const Double_t* v) { fX = v[0]; fY = v[1]; fZ = v[2]; }
   void Set(TT x, TT  y, TT z) { fX = x; fY = y; fZ = z; }
   void Set(const TVector3& v);

   template <typename OO>
   void Set(const TEveVectorT<OO>& v) { fX = v.fX;  fY = v.fY;  fZ = v.fZ; }

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

   TT   Distance(const TEveVectorT& v) const;
   TT   SquareDistance(const TEveVectorT& v) const;

   TT   Dot(const TEveVectorT& a) const;

   TEveVectorT  Cross(const TEveVectorT& a) const;

   TEveVectorT& Sub(const TEveVectorT& a, const TEveVectorT& b);
   TEveVectorT& Mult(const TEveVectorT& a, TT af);

   TEveVectorT  Orthogonal() const;
   void         OrthoNormBase(TEveVectorT& a, TEveVectorT& b) const;

   Bool_t       IsZero() const { return fX == 0 && fY == 0 && fZ == 0; }

   ClassDefNV(TEveVectorT, 2); // A three-vector template without TObject inheritance and virtual functions.
};

typedef TEveVectorT<Float_t>  TEveVector;
typedef TEveVectorT<Float_t>  TEveVectorF;
typedef TEveVectorT<Double_t> TEveVectorD;

//______________________________________________________________________________
template<typename TT>
inline TT TEveVectorT<TT>::Phi() const
{
   return fX == 0 && fY == 0 ? 0 : TMath::ATan2(fY, fX);
}

//______________________________________________________________________________
template<typename TT>
inline TT TEveVectorT<TT>::Theta() const
{
   return fX == 0 && fY == 0 && fZ == 0 ? 0 : TMath::ATan2(Perp(), fZ);
}

//______________________________________________________________________________
template<typename TT>
inline TT TEveVectorT<TT>::CosTheta() const
{
   Float_t ptot = Mag(); return ptot == 0 ? 1 : fZ/ptot;
}

//______________________________________________________________________________
template<typename TT>
inline TT TEveVectorT<TT>::Distance(const TEveVectorT& b) const
{
   return TMath::Sqrt((fX - b.fX)*(fX - b.fX) +
                      (fY - b.fY)*(fY - b.fY) +
                      (fZ - b.fZ)*(fZ - b.fZ));
}

//______________________________________________________________________________
template<typename TT>
inline TT TEveVectorT<TT>::SquareDistance(const TEveVectorT& b) const
{
   return ((fX - b.fX) * (fX - b.fX) +
           (fY - b.fY) * (fY - b.fY) +
           (fZ - b.fZ) * (fZ - b.fZ));
}

//______________________________________________________________________________
template<typename TT>
inline TT TEveVectorT<TT>::Dot(const TEveVectorT& a) const
{
   return a.fX*fX + a.fY*fY + a.fZ*fZ;
}

//______________________________________________________________________________
template<typename TT>
inline TEveVectorT<TT> TEveVectorT<TT>::Cross(const TEveVectorT<TT>& a) const
{
   TEveVectorT<TT> r;
   r.fX = fY * a.fZ - fZ * a.fY;
   r.fY = fZ * a.fX - fX * a.fZ;
   r.fZ = fX * a.fY - fY * a.fX;
   return r;
}

//______________________________________________________________________________
template<typename TT>
inline TEveVectorT<TT>& TEveVectorT<TT>::Sub(const TEveVectorT<TT>& a, const TEveVectorT<TT>& b)
{
   fX = a.fX - b.fX;
   fY = a.fY - b.fY;
   fZ = a.fZ - b.fZ;
   return *this;
}

//______________________________________________________________________________
template<typename TT>
inline TEveVectorT<TT>& TEveVectorT<TT>::Mult(const TEveVectorT<TT>& a, TT af)
{
   fX = a.fX * af;
   fY = a.fY * af;
   fZ = a.fZ * af;
   return *this;
}

//______________________________________________________________________________
template<typename TT>
inline TEveVectorT<TT> operator+(const TEveVectorT<TT>& a, const TEveVectorT<TT>& b)
{
   TEveVectorT<TT> r(a);
   return r += b;
}

//______________________________________________________________________________
template<typename TT>
inline TEveVectorT<TT> operator-(const TEveVectorT<TT>& a, const TEveVectorT<TT>& b)
{
   TEveVectorT<TT> r(a);
   return r -= b;
}

//______________________________________________________________________________
template<typename TT>
inline TEveVectorT<TT> operator*(const TEveVectorT<TT>& a, TT b)
{
   TEveVectorT<TT> r(a);
   return r *= b;
}

//______________________________________________________________________________
template<typename TT>
inline TEveVectorT<TT> operator*(TT b, const TEveVectorT<TT>& a)
{
   TEveVectorT<TT> r(a);
   return r *= b;
}


//==============================================================================
// TEveVector4T
//==============================================================================

template <typename TT>
class TEveVector4T : public TEveVectorT<TT>
{
   typedef TEveVectorT<TT> TP;

public:
   TT fT;

   TEveVector4T() : TP(),  fT(0) {}
   template <typename OO>
   TEveVector4T(const TEveVectorT<OO>& v) : TP(v.fX, v.fY, v.fZ), fT(0) {}
   template <typename OO>
   TEveVector4T(const TEveVectorT<OO>& v, Float_t t) : TP(v.fX, v.fY, v.fZ), fT(t) {}
   template <typename OO>
   TEveVector4T(const TEveVector4T<OO>& v) : TP(v.fX, v.fY, v.fZ), fT(v.fT) {}
   TEveVector4T(const Float_t*  v) : TP(v), fT(v[3]) {}
   TEveVector4T(const Double_t* v) : TP(v), fT(v[3]) {}
   TEveVector4T(TT x, TT y, TT z, TT t=0) : TP(x, y, z), fT(t) {}

   void Dump() const;

   TEveVector4T& operator*=(TT s)                  { TP::operator*=(s); fT *= s;    return *this; }
   TEveVector4T& operator+=(const TEveVector4T& v) { TP::operator+=(v); fT += v.fT; return *this; }
   TEveVector4T& operator-=(const TEveVector4T& v) { TP::operator-=(v); fT -= v.fT; return *this; }

   using TP::operator+=;
   using TP::operator-=;

   ClassDefNV(TEveVector4T, 1); // A four-vector template without TObject inheritance and virtual functions.
};

typedef TEveVector4T<Float_t>  TEveVector4;
typedef TEveVector4T<Float_t>  TEveVector4F;
typedef TEveVector4T<Double_t> TEveVector4D;

//______________________________________________________________________________
template<typename TT>
inline TEveVector4T<TT> operator+(const TEveVector4T<TT>& a, const TEveVector4T<TT>& b)
{
   return TEveVector4T<TT>(a.fX + b.fX, a.fY + b.fY, a.fZ + b.fZ, a.fT + b.fT);
}

//______________________________________________________________________________
template<typename TT>
inline TEveVector4T<TT> operator-(const TEveVector4T<TT>& a, const TEveVector4T<TT>& b)
{
   return TEveVector4T<TT>(a.fX - b.fX, a.fY - b.fY, a.fZ - b.fZ, a.fT - b.fT);
}

//______________________________________________________________________________
template<typename TT>
inline TEveVector4T<TT> operator*(const TEveVector4T<TT>& a, TT b)
{
   return TEveVector4T<TT>(a.fX*b, a.fY*b, a.fZ*b, a.fT*b);
}

//______________________________________________________________________________
template<typename TT>
inline TEveVector4T<TT> operator*(TT b, const TEveVector4T<TT>& a)
{
   return TEveVector4T<TT>(a.fX*b, a.fY*b, a.fZ*b, a.fT*b);
}


//==============================================================================
// TEveVector2T
//==============================================================================

template <typename TT>
class TEveVector2T
{
public:
   TT fX, fY; // Components of the point.

   TEveVector2T() : fX(0), fY(0) {}
   template <typename OO>
   TEveVector2T(const TEveVector2T<OO>& v) : fX(v.fX), fY(v.fY) {}
   TEveVector2T(const Float_t* v)  : fX(v[0]), fY(v[1]) {}
   TEveVector2T(const Double_t* v) : fX(v[0]), fY(v[1]) {}
   TEveVector2T(TT x, TT y) : fX(x), fY(y)    {}

   void Dump() const;

   operator const TT*() const { return &fX; }
   operator       TT*()       { return &fX; }

   TEveVector2T& operator*=(TT s)                  { fX *= s;    fY *= s;    return *this; }
   TEveVector2T& operator+=(const TEveVector2T& v) { fX += v.fX; fY += v.fY; return *this; }
   TEveVector2T& operator-=(const TEveVector2T& v) { fX -= v.fX; fY -= v.fY; return *this; }

   TT& operator[](Int_t idx)       { return (&fX)[idx]; }
   TT  operator[](Int_t idx) const { return (&fX)[idx]; }

   const TT* Arr() const { return &fX; }
   TT* Arr()             { return &fX; }

   void Set(const Float_t*  v) { fX = v[0]; fY = v[1]; }
   void Set(const Double_t* v) { fX = v[0]; fY = v[1]; }
   void Set(TT x, TT y) { fX = x; fY = y; }

   template <typename OO>
   void Set(const TEveVector2T<OO>& v) { fX = v.fX; fY = v.fY; }

   void NegateXY() { fX = - fX; fY = -fY; }
   void Normalize(TT length=1);

   TT Phi()  const;

   TT Mag2() const { return fX*fX + fY*fY;}
   TT Mag()  const { return TMath::Sqrt(Mag2());}

   TT Distance(const TEveVector2T& v) const;
   TT SquareDistance(const TEveVector2T& v) const;

   TT Dot(const TEveVector2T& a) const;
   TT Cross(const TEveVector2T& a) const;

   TEveVector2T& Sub(const TEveVector2T& p, const TEveVector2T& q);

   TEveVector2T& Mult(const TEveVector2T& a, TT af);

   ClassDefNV(TEveVector2T, 1); // // A two-vector template without TObject inheritance and virtual functions.
};

typedef TEveVector2T<Float_t>  TEveVector2;
typedef TEveVector2T<Float_t>  TEveVector2F;
typedef TEveVector2T<Double_t> TEveVector2D;

//______________________________________________________________________________
template<typename TT>
inline TT TEveVector2T<TT>::Phi() const
{
   return fX == 0.0 && fY == 0.0 ? 0.0 : TMath::ATan2(fY, fX);
}

//______________________________________________________________________________
template<typename TT>
inline TT TEveVector2T<TT>::Distance( const TEveVector2T<TT>& b) const
{
   return TMath::Sqrt((fX - b.fX)*(fX - b.fX) +
                      (fY - b.fY)*(fY - b.fY));
}

//______________________________________________________________________________
template<typename TT>
inline TT TEveVector2T<TT>::SquareDistance(const TEveVector2T<TT>& b) const
{
   return ((fX - b.fX) * (fX - b.fX) +
           (fY - b.fY) * (fY - b.fY));
}

//______________________________________________________________________________
template<typename TT>
inline TT TEveVector2T<TT>::Dot(const TEveVector2T<TT>& a) const
{
   return a.fX*fX + a.fY*fY;
}

//______________________________________________________________________________
template<typename TT>
inline TT TEveVector2T<TT>::Cross(const TEveVector2T<TT>& a) const
{
   return fX * a.fY - fY * a.fX;
}

//______________________________________________________________________________
template<typename TT>
inline TEveVector2T<TT>& TEveVector2T<TT>::Sub(const TEveVector2T<TT>& p, const TEveVector2T<TT>& q)
{
   fX = p.fX - q.fX;
   fY = p.fY - q.fY;
   return *this;
}

//______________________________________________________________________________
template<typename TT>
inline TEveVector2T<TT>& TEveVector2T<TT>::Mult(const TEveVector2T<TT>& a, TT af)
{
   fX = a.fX * af;
   fY = a.fY * af;
   return *this;
}

//______________________________________________________________________________
template<typename TT>
inline TEveVector2T<TT> operator+(const TEveVector2T<TT>& a, const TEveVector2T<TT>& b)
{
   TEveVector2T<TT> r(a);
   return r += b;
}

//______________________________________________________________________________
template<typename TT>
inline TEveVector2T<TT> operator-(const TEveVector2T<TT>& a, const TEveVector2T<TT>& b)
{
   TEveVector2T<TT> r(a);
   return r -= b;
}

//______________________________________________________________________________
template<typename TT>
inline TEveVector2T<TT> operator*(const TEveVector2T<TT>& a, TT b)
{
   TEveVector2T<TT> r(a);
   return r *= b;
}

//______________________________________________________________________________
template<typename TT>
inline TEveVector2T<TT> operator*(TT b, const TEveVector2T<TT>& a)
{
   TEveVector2T<TT> r(a);
   return r *= b;
}

#endif
