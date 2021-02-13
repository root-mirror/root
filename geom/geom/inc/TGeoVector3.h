// @(#)root/geom:$Id$
// Author: Andrei Gheata   20/12/19

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGeoVector3
#define ROOT_TGeoVector3

#include <Riostream.h>
#include <TMath.h>

namespace ROOT {
namespace Geom {

struct Vertex_t {
   double fVec[3] = {0.};

   Vertex_t(const double a, const double b, const double c)
   {
      fVec[0] = a;
      fVec[1] = b;
      fVec[2] = c;
   }

   Vertex_t(const double a = 0.)
   {
      fVec[0] = a;
      fVec[1] = a;
      fVec[2] = a;
   }

   double &operator[](const int index) { return fVec[index]; }
   double const &operator[](const int index) const { return fVec[index]; }

   // Inplace binary operators

#define Vertex_t_INPLACE_BINARY_OP(OPERATOR)                 \
   inline Vertex_t &operator OPERATOR(const Vertex_t &other) \
   {                                                         \
      fVec[0] OPERATOR other.fVec[0];                        \
      fVec[1] OPERATOR other.fVec[1];                        \
      fVec[2] OPERATOR other.fVec[2];                        \
      return *this;                                          \
   }                                                         \
   inline Vertex_t &operator OPERATOR(const double &scalar)  \
   {                                                         \
      fVec[0] OPERATOR scalar;                               \
      fVec[1] OPERATOR scalar;                               \
      fVec[2] OPERATOR scalar;                               \
      return *this;                                          \
   }
   Vertex_t_INPLACE_BINARY_OP(+=) Vertex_t_INPLACE_BINARY_OP(-=) Vertex_t_INPLACE_BINARY_OP(*=)
      Vertex_t_INPLACE_BINARY_OP(/=)
#undef Vertex_t_INPLACE_BINARY_OP

         double &x()
   {
      return fVec[0];
   }
   double const &x() const { return fVec[0]; }

   double &y() { return fVec[1]; }
   double const &y() const { return fVec[1]; }

   double &z() { return fVec[2]; }
   double const &z() const { return fVec[2]; }

   inline void CopyTo(double *dest) const
   {
      dest[0] = fVec[0];
      dest[1] = fVec[1];
      dest[2] = fVec[2];
   }

   void Set(double const &a, double const &b, double const &c)
   {
      fVec[0] = a;
      fVec[1] = b;
      fVec[2] = c;
   }

   void Set(const double a) { Set(a, a, a); }

   /// \Return the length squared perpendicular to z direction
   double Perp2() const { return fVec[0] * fVec[0] + fVec[1] * fVec[1]; }

   /// \Return the length perpendicular to z direction
   double Perp() const { return TMath::Sqrt(Perp2()); }

   /// The dot product of two vector objects
   static double Dot(Vertex_t const &left, Vertex_t const &right)
   {
      return left[0] * right[0] + left[1] * right[1] + left[2] * right[2];
   }

   /// The dot product of two vector
   double Dot(Vertex_t const &right) const { return Dot(*this, right); }

   /// \return Squared magnitude of the vector.
   double Mag2() const { return Dot(*this, *this); }

   /// \return Magnitude of the vector.
   double Mag() const { return TMath::Sqrt(Mag2()); }

   double Length() const { return Mag(); }

   double Length2() const { return Mag2(); }

   /// Normalizes the vector by dividing each entry by the length.
   void Normalize() { *this *= (1. / Length()); }

   // Vertex_t Normalized() const { return Vertex_t(*this) * (1. / Length()); }

   // checks if vector is normalized
   bool IsNormalized() const
   {
      double norm = Mag2();
      constexpr double tolerance = 1.e-10;
      return 1. - tolerance < norm && norm < 1 + tolerance;
   }

   /// \return Azimuthal angle between -pi and pi.
   double Phi() const { return TMath::ATan2(fVec[1], fVec[0]); }

   /// \return Polar angle between 0 and pi.
   double Theta() const { return TMath::ACos(fVec[2] / Mag()); }

   /// The cross (vector) product of two Vector3D<T> objects
   static Vertex_t Cross(Vertex_t const &left, Vertex_t const &right)
   {
      return Vertex_t(left[1] * right[2] - left[2] * right[1], left[2] * right[0] - left[0] * right[2],
                      left[0] * right[1] - left[1] * right[0]);
   }

   Vertex_t Abs() const { return Vertex_t(TMath::Abs(fVec[0]), TMath::Abs(fVec[1]), TMath::Abs(fVec[2])); }

   double Min() const { return TMath::Min(TMath::Min(fVec[0], fVec[1]), fVec[2]); }

   double Max() const { return TMath::Max(TMath::Max(fVec[0], fVec[1]), fVec[2]); }

   Vertex_t Unit() const
   {
      constexpr double kMinimum = std::numeric_limits<double>::min();
      const double mag2 = Mag2();
      Vertex_t output(*this);
      output /= TMath::Sqrt(mag2 + kMinimum);
      return output;
   }
};

inline bool operator==(Vertex_t const &lhs, Vertex_t const &rhs)
{
   constexpr double kTolerance = 1.e-8;
   return TMath::Abs(lhs[0] - rhs[0]) < kTolerance && TMath::Abs(lhs[1] - rhs[1]) < kTolerance &&
          TMath::Abs(lhs[2] - rhs[2]) < kTolerance;
}

inline bool operator!=(Vertex_t const &lhs, Vertex_t const &rhs)
{
   return !(lhs == rhs);
}

#define Vertex_t_BINARY_OP(OPERATOR, INPLACE)                                  \
   inline Vertex_t operator OPERATOR(const Vertex_t &lhs, const Vertex_t &rhs) \
   {                                                                           \
      Vertex_t result(lhs);                                                    \
      result INPLACE rhs;                                                      \
      return result;                                                           \
   }                                                                           \
   inline Vertex_t operator OPERATOR(Vertex_t const &lhs, const double rhs)    \
   {                                                                           \
      Vertex_t result(lhs);                                                    \
      result INPLACE rhs;                                                      \
      return result;                                                           \
   }                                                                           \
   inline Vertex_t operator OPERATOR(const double lhs, Vertex_t const &rhs)    \
   {                                                                           \
      Vertex_t result(lhs);                                                    \
      result INPLACE rhs;                                                      \
      return result;                                                           \
   }
Vertex_t_BINARY_OP(+, +=) Vertex_t_BINARY_OP(-, -=) Vertex_t_BINARY_OP(*, *=) Vertex_t_BINARY_OP(/, /=)
#undef Vertex_t_BINARY_OP

} // namespace Geom
} // namespace ROOT

std::ostream &operator<<(std::ostream &os, ROOT::Geom::Vertex_t const &vec);

#endif
