/// \defgroup ROOT7PadCoordSystems ROOT7 RPad coordinate systems
/// \ingroup ROOT7Graphics
/// \brief The ROOT7 RPad coordinate systems.
///
/// \author Axel Naumann <axel@cern.ch>
/// \date 2017-07-06

/*************************************************************************
 * Copyright (C) 1995-2017, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RPadLength
#define ROOT7_RPadLength

#include <vector>
#include <string>
#include <cmath>

namespace ROOT {
namespace Experimental {

   /** \class RPadLength
       \ingroup ROOT7PadCoordSystems
       \brief A length in RPad.
    */

class RPadLength  {

protected:

   std::vector<double> fArr; ///< components [0] - normalized, [1] - pixel, [2] - user

public:
   template <class DERIVED>
   struct CoordSysBase {
      double  fVal{0.};      ///<! Coordinate value

      CoordSysBase() = default;
      CoordSysBase(double val): fVal(val) {}
      DERIVED &ToDerived() { return static_cast<DERIVED &>(*this); }

      DERIVED operator-() { return DERIVED(-fVal); }

      friend DERIVED operator+(DERIVED lhs, DERIVED rhs) { return DERIVED{lhs.fVal + rhs.fVal}; }
      friend DERIVED operator-(DERIVED lhs, DERIVED rhs) { return DERIVED{lhs.fVal - rhs.fVal}; }
      friend double operator/(DERIVED lhs, DERIVED rhs) { return lhs.fVal / rhs.fVal; }
      DERIVED &operator+=(const DERIVED &rhs)
      {
         fVal += rhs.fVal;
         return ToDerived();
      }
      DERIVED &operator-=(const DERIVED &rhs)
      {
         fVal -= rhs.fVal;
         return ToDerived();
      }
      DERIVED &operator*=(double scale)
      {
         fVal *= scale;
         return ToDerived();
      }
      friend DERIVED operator*(const DERIVED &lhs, double rhs) { return DERIVED(lhs.fVal * rhs); }
      friend DERIVED operator*(double lhs, const DERIVED &rhs) { return DERIVED(lhs * rhs.fVal); }
      friend DERIVED operator/(const DERIVED &lhs, double rhs) { return DERIVED(lhs.fVal * rhs); }
      friend bool operator<(const DERIVED &lhs, const DERIVED &rhs) { return lhs.fVal < rhs.fVal; }
      friend bool operator>(const DERIVED &lhs, const DERIVED &rhs) { return lhs.fVal > rhs.fVal; }
      friend bool operator<=(const DERIVED &lhs, const DERIVED &rhs) { return lhs.fVal <= rhs.fVal; }
      friend bool operator>=(const DERIVED &lhs, const DERIVED &rhs) { return lhs.fVal >= rhs.fVal; }
      // no ==, !=
   };

   /// \defgroup TypesafeCoordinates Typesafe Coordinates
   /// \ingroup ROOT7PadCoordSystems
   /// These define typesafe coordinates used by RPad to identify which coordinate system a coordinate is referring to
   /// The origin (0,0) is in the `RPad`'s bottom left corner for all of them.
   /// \{

   /** \struct Normal
       \brief A normalized coordinate.

     0 in the left, bottom corner, 1 in the top, right corner of the `RPad`.
     Resizing the pad will resize the objects with it.
    */
   struct Normal: CoordSysBase<Normal> {
      using CoordSysBase<Normal>::CoordSysBase;
   };

   /** \struct Pixel
       \brief A pixel coordinate.

     0 in the left, bottom corner, 1 in the top, right corner of the `RPad`.
     Resizing the pad will keep the pixel-position of the objects positioned in `Pixel` coordinates.
    */
   struct Pixel: CoordSysBase<Pixel> {
      using CoordSysBase<Pixel>::CoordSysBase;
   };

   /** \struct User
       \brief A user coordinate.

       as defined by the EUserCoordSystem parameter of the `RPad`.
    */
   struct User: CoordSysBase<User> {
      using CoordSysBase<User>::CoordSysBase;
   };
   /// \}

   RPadLength() {}

   /// Constructor from a `Normal` coordinate.
   RPadLength(Normal normal): RPadLength() { SetNormal(normal.fVal); }

   /// By default all numeric values are normal values
   RPadLength(double normal): RPadLength() { SetNormal(normal); }

   /// Constructor from a `Pixel` coordinate.
   RPadLength(Pixel px): RPadLength() { SetPixel(px.fVal); }

   /// Constructor from a `User` coordinate.
   RPadLength(User user) : RPadLength() { SetUser(user.fVal); }

   /// Constructor for normal and pixel coordinate.
   RPadLength(Normal normal, Pixel px): RPadLength() { SetPixel(px.fVal); SetNormal(normal.fVal);  }

   /// Constructor for normal, pixel and user coordinate.
   RPadLength(Normal normal, Pixel px, User user): RPadLength() { SetUser(user.fVal); SetPixel(px.fVal); SetNormal(normal.fVal);  }

   /// Constructor from string representation
   RPadLength(const std::string &csscode) : RPadLength() { if (!csscode.empty()) ParseString(csscode); }

   bool HasNormal() const { return fArr.size() > 0; }
   bool HasPixel() const { return fArr.size() > 1; }
   bool HasUser() const { return fArr.size() > 2; }

   RPadLength &SetNormal(double v)
   {
      if (fArr.size() < 1)
         fArr.resize(1);
      fArr[0] = v;
      return *this;
   }

   RPadLength &SetPixel(double v)
   {
      if (fArr.size() < 2)
         fArr.resize(2, 0.);
      fArr[1] = v;
      return *this;
   }

   RPadLength &SetUser(double v)
   {
      if (fArr.size() < 3)
         fArr.resize(3, 0.);
      fArr[2] = v;
      return *this;
   }

   double GetNormal() const { return fArr.size() > 0 ? fArr[0] : 0.; }
   double GetPixel() const { return fArr.size() > 1 ? fArr[1] : 0.; }
   double GetUser() const { return fArr.size() > 2 ? fArr[2] : 0.; }

   void ClearUser()
   {
      if (fArr.size() > 2)
         fArr.resize(2);
   }

   void ClearPixelAndUser()
   {
      if (fArr.size() > 1)
         fArr.resize(1);
   }

   void Clear() { fArr.clear(); }

   bool Empty() const { return fArr.size() == 0; }

   /// Add two `RPadLength`s.
   friend RPadLength operator+(RPadLength lhs, const RPadLength &rhs)
   {
      RPadLength res;
      if (lhs.HasUser() || rhs.HasUser())
         res.SetUser(lhs.GetUser() + rhs.GetUser());
      if (lhs.HasPixel() || rhs.HasPixel())
         res.SetPixel(lhs.GetPixel() + rhs.GetPixel());
      if (lhs.HasNormal() || rhs.HasNormal())
         res.SetNormal(lhs.GetNormal() + rhs.GetNormal());
      return res;
   }

   /// Subtract two `RPadLength`s.
   friend RPadLength operator-(RPadLength lhs, const RPadLength &rhs)
   {
      RPadLength res;
      if (lhs.HasUser() || rhs.HasUser())
         res.SetUser(lhs.GetUser() - rhs.GetUser());
      if (lhs.HasPixel() || rhs.HasPixel())
         res.SetPixel(lhs.GetPixel() - rhs.GetPixel());
      if (lhs.HasNormal() || rhs.HasNormal())
         res.SetNormal(lhs.GetNormal() - rhs.GetNormal());
      return res;
   }

   /// Unary -.
   RPadLength operator-()
   {
      RPadLength res;
      if (HasUser()) res.SetUser(-GetUser());
      if (HasPixel()) res.SetPixel(-GetPixel());
      if (HasNormal()) res.SetNormal(-GetNormal());
      return res;
   }

   /// Add a `RPadLength`.
   RPadLength &operator+=(const RPadLength &rhs)
   {
      if (HasUser() || rhs.HasUser())
         SetUser(GetUser() + rhs.GetUser());
      if (HasPixel() || rhs.HasPixel())
         SetPixel(GetPixel() + rhs.GetPixel());
      if (HasNormal() || rhs.HasNormal())
         SetNormal(GetNormal() + rhs.GetNormal());
      return *this;
   }

   /// Subtract a `RPadLength`.
   RPadLength &operator-=(const RPadLength &rhs)
   {
      if (HasUser() || rhs.HasUser())
         SetUser(GetUser() - rhs.GetUser());
      if (HasPixel() || rhs.HasPixel())
         SetPixel(GetPixel() - rhs.GetPixel());
      if (HasNormal() || rhs.HasNormal())
         SetNormal(GetNormal() - rhs.GetNormal());
      return *this;
   }

   /// Multiply a `RPadLength`.
   RPadLength &operator*=(double scale)
   {
      if (HasUser()) SetUser(scale*GetUser());
      if (HasPixel()) SetPixel(scale*GetPixel());
      if (HasNormal()) SetNormal(scale*GetNormal());
      return *this;
   }

   /// Compare a `RPadLength`.
   bool operator==(const RPadLength &rhs) const
   {
      if ((HasUser() != rhs.HasUser()) ||
          (HasUser() && std::fabs(GetUser() - rhs.GetUser()) > 1e-30)) return false;

      if(std::fabs(GetPixel() - rhs.GetPixel()) > 1e-4) return false;

      if (std::fabs(GetNormal() - rhs.GetNormal()) > 1e-6) return false;

      return true;
   };

   std::string AsString() const;

   bool ParseString(const std::string &val);

};

/// User-defined literal for `RPadLength::Normal`
///
/// Use as
/// ```
/// using namespace ROOT::Experimental;
/// RLine(0.1_normal, 0.0_normal, RLineExtent(0.2_normal, 0.5_normal));
/// ```
inline RPadLength::Normal operator"" _normal(long double val)
{
   return RPadLength::Normal{(double)val};
}
inline RPadLength::Normal operator"" _normal(unsigned long long int val)
{
   return RPadLength::Normal{(double)val};
}

/// User-defined literal for `RPadLength::Pixel`
///
/// Use as
/// ```
/// using namespace ROOT::Experimental;
/// RLine(100_px, 0_px, RLineExtent(20_px, 50_px));
/// ```
inline RPadLength::Pixel operator"" _px(long double val)
{
   return RPadLength::Pixel{(double)val};
}
inline RPadLength::Pixel operator"" _px(unsigned long long int val)
{
   return RPadLength::Pixel{(double)val};
}

/// User-defined literal for `RPadLength::User`
///
/// Use as
/// ```
/// using namespace ROOT::Experimental;
/// RLine(0.1_user, 0.0_user, RLineExtent(0.2_user, 0.5_user));
/// ```
inline RPadLength::User operator"" _user(long double val)
{
   return RPadLength::User{(double)val};
}
inline RPadLength::User operator"" _user(unsigned long long int val)
{
   return RPadLength::User{(double)val};
}

} // namespace Experimental
} // namespace ROOT

#endif
