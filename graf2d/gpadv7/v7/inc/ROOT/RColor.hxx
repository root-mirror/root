/// \file ROOT/RColorOld.hxx
/// \ingroup Gpad ROOT7
/// \author Axel Naumann <axel@cern.ch>
/// \date 2017-09-26
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2017, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RColorOld
#define ROOT7_RColorOld

#include <array>
#include <vector>
#include <string>

#include <ROOT/RDrawingAttr.hxx>

namespace ROOT {
namespace Experimental {

/** \class ROOT::Experimental::RColorOld
  A color: Red|Green|Blue|Alpha, or a position in a RPalette
  */
class RColorOld {
public:
   /** \class ROOT::Experimental::RColorOld::TAlpha
    The alpha value of a color: 0 is completely transparent, 1 is completely opaque.
    */
   struct Alpha {
      float fVal;
      explicit operator float() const { return fVal; }
   };
   /// An opaque color.
   static constexpr Alpha kOpaque{1.};
   /// A completely transparent color.
   static constexpr Alpha kTransparent{0.};

   enum class EKind {
      kRGBA, ///< The color is defined as specific RGBA values.
      kPalettePos, ///< The color is defined as a value in the `RFrame`'s `RPalette`.
      kAuto ///< The color will be set upon drawing the canvas choosing a `RPalette` color, see `RColorOld(Auto_t)`
   };

private:
   // TODO: use a `variant` here!
   /// The "R" in RGBA (0 <= R <= 1), or the palette pos if fKind is `kPalettePos`.
   float fRedOrPalettePos = 0.;

   /// The "G" in RGBA (0 <= G <= 1). Unused if `fKind != kRGBA`.
   float fGreen = 0.;

   /// The "B" in RGBA (0 <= B <= 1). Unused if `fKind != kRGBA`.
   float fBlue = 0.;

   /// The "A" in RGBA (0 <= A <= 1). Unused if `fKind != kRGBA`. `fAlpha == 0` means so transparent it's invisible,
   /// `fAlpha == 1` means completely opaque.
   float fAlpha = 1.;

   /// How the color is defined.
   EKind fKind = EKind::kRGBA;

   /// throw an exception if the color isn't specified as `kRGBA` or `kAuto`, the two cases where
   /// asking for RBGA members makes sense.
   bool AssertNotPalettePos() const;

public:
   using RGBA = std::array<float, 4>;

   // Default constructor: good old solid black.
   constexpr RColorOld() = default;

   /// Initialize a RColorOld with red, green, blue and alpha component.
   constexpr RColorOld(float r, float g, float b, float alpha): fRedOrPalettePos(r), fGreen(g), fBlue(b), fAlpha(alpha) {}

   /// Initialize a RColorOld with red, green, blue and alpha component.
   constexpr RColorOld(float r, float g, float b, Alpha alpha = kOpaque): RColorOld(r, g, b, alpha.fVal) {}

   /// Initialize a RColorOld with red, green, blue and alpha component as an array.
   constexpr RColorOld(const RGBA &rgba): RColorOld(rgba[0], rgba[1], rgba[2], rgba[3]) {}

   /// Initialize a `RColorOld` with a `RPalette` ordinal. The actual color is determined from the pad's
   /// (or rather its `RFrame`'s) `RPalette`
   constexpr RColorOld(float paletteOrdinal): fRedOrPalettePos(paletteOrdinal), fKind(EKind::kPalettePos) {}

   /**\class AutoTag
    Used to signal that this color shall be automatically chosen by the drawing routines, by picking a color
    from the `RPad`'s (or rather its `RFrame`'s) current `RPalette`.
   */
   class AutoTag {};

   /// Constructs an automatically assigned color. Call as `RColorOld col(RColorOld::kAuto)`.
   constexpr RColorOld(AutoTag): fKind(EKind::kAuto) {}

   /// Determine whether this RColorOld is storing RGBA (in contrast to an ordinal of a RPalette).
   bool IsRGBA() const { return fKind == EKind::kRGBA; }

   /// Determine whether this `RColorOld` is storing an ordinal of a RPalette (in contrast to RGBA).
   bool IsPaletteOrdinal() const { return fKind == EKind::kPalettePos; }

   /// Determine whether this `RColorOld` will be assigned a actual color upon drawing.
   bool IsAuto() const { return fKind == EKind::kAuto; }

   /// If this is an ordinal in a palette, resolve the
   float GetPaletteOrdinal() const;

   friend bool operator==(const RColorOld &lhs, const RColorOld &rhs)
   {
      if (lhs.fKind != rhs.fKind)
         return false;
      switch (lhs.fKind) {
      case EKind::kPalettePos:
         return lhs.fRedOrPalettePos == rhs.fRedOrPalettePos;
      case EKind::kRGBA:
         return lhs.fRedOrPalettePos == rhs.fRedOrPalettePos && lhs.fGreen == rhs.fGreen && lhs.fBlue == rhs.fBlue &&
             lhs.fAlpha == rhs.fAlpha;
      case EKind::kAuto:
         return true; // is that what we need?
      }
      return false;
   }

   /// For RGBA or auto colors, get the red component (0..1).
   float GetRed() const {
      if (AssertNotPalettePos())
         return fRedOrPalettePos;
      return 0.;
   }

   /// For RGBA or auto colors, get the green component (0..1).
   float GetGreen() const {
      if (AssertNotPalettePos())
         return fGreen;
      return 0.;
   }

   /// For RGBA or auto colors, get the blue component (0..1).
   float GetBlue() const {
      if (AssertNotPalettePos())
         return fBlue;
      return 0.;
   }

   /// For RGBA or auto colors, get the alpha component (0..1).
   float GetAlpha() const {
      if (AssertNotPalettePos())
         return fAlpha;
      return 0.;
   }

   /// For RGBA or auto colors, set the red component.
   void SetRed(float r) {
      if (AssertNotPalettePos())
         fRedOrPalettePos = r;
   }

   /// For RGBA or auto colors, set the green component.
   void SetGreen(float g) {
      if (AssertNotPalettePos())
         fGreen = g;
   }

   /// For RGBA or auto colors, set the blue component.
   void SetBlue(float b) {
      if (AssertNotPalettePos())
         fBlue = b;
   }

   /// For RGBA or auto colors, set the alpha component.
   void SetAlpha(float a) {
      if (AssertNotPalettePos())
         fAlpha = a;
   }

   /// For RGBA or auto colors, set the alpha component.
   void SetAlpha(Alpha a) {
      if (AssertNotPalettePos())
         fAlpha = (float)a;
   }

   /// Return the Hue, Light, Saturation (HLS) definition of this RColorOld
   void GetHLS(float &hue, float &light, float &satur) {
      hue = light = satur = 0.;
      if (AssertNotPalettePos()) {
         float rnorm, gnorm, bnorm, minval, maxval, msum, mdiff;
         minval = maxval =0 ;

         minval = fRedOrPalettePos;
         if (fGreen < minval) minval = fGreen;
         if (fBlue < minval)  minval = fBlue;
         maxval = fRedOrPalettePos;
         if (fGreen > maxval) maxval = fGreen;
         if (fBlue > maxval)  maxval = fBlue;

         rnorm = gnorm = bnorm = 0;
         mdiff = maxval - minval;
         msum  = maxval + minval;
         light = 0.5 * msum;
         if (maxval != minval) {
            rnorm = (maxval - fRedOrPalettePos)/mdiff;
            gnorm = (maxval - fGreen)/mdiff;
            bnorm = (maxval - fBlue)/mdiff;
         } else {
            satur = hue = 0;
            return;
         }

         if (light < 0.5) satur = mdiff/msum;
         else             satur = mdiff/(2.0 - msum);

         if      (fRedOrPalettePos == maxval) hue = 60.0 * (6.0 + bnorm - gnorm);
         else if (fGreen == maxval)           hue = 60.0 * (2.0 + rnorm - bnorm);
         else                                 hue = 60.0 * (4.0 + gnorm - rnorm);

         if (hue > 360) hue = hue - 360;
      }
   }

   /// Set the Red Green and Blue (RGB) values from the Hue, Light, Saturation (HLS).
   void SetRGBFromHLS(float hue, float light, float satur) {
      if (AssertNotPalettePos()) {
         float rh, rl, rs, rm1, rm2;
         rh = rl = rs = 0;
         if (hue   > 0) { rh = hue;   if (rh > 360) rh = 360; }
         if (light > 0) { rl = light; if (rl > 1)   rl = 1; }
         if (satur > 0) { rs = satur; if (rs > 1)   rs = 1; }

         if (rl <= 0.5) rm2 = rl*(1.0 + rs);
         else           rm2 = rl + rs - rl*rs;
         rm1 = 2.0*rl - rm2;

         if (!rs) { fRedOrPalettePos = fGreen = fBlue = rl; return; }

         auto toRGB = [rm1, rm2] (float h) {
            if (h > 360) h = h - 360;
            if (h < 0)   h = h + 360;
            if (h < 60 ) return rm1 + (rm2-rm1)*h/60;
            if (h < 180) return rm2;
            if (h < 240) return rm1 + (rm2-rm1)*(240-h)/60;
            return rm1;
         };

         fRedOrPalettePos = toRGB(rh+120);
         fGreen           = toRGB(rh);
         fBlue            = toRGB(rh-120);
      }
   }

   ///\{
   ///\name Default colors
   static constexpr RGBA kRed{{1., 0., 0., 1.}};
   static constexpr RGBA kGreen{{0., 1., 0., 1.}};
   static constexpr RGBA kBlue{{0., 0, 1., 1.}};
   static constexpr RGBA kWhite{{1., 1, 1., 1.}};
   static constexpr RGBA kBlack{{0., 0., 0., 1.}};
   static constexpr RGBA kInvisible{{0., 0., 0., 0.}};
   static constexpr AutoTag kAuto{};
   ///\}
};

// TODO: see also imagemagick's C++ interface for RColorOld operations!
// https://www.imagemagick.org/api/magick++-classes.php


class RColor : public RAttributesVisitor {

protected:
   const RDrawableAttributes::Map_t &GetDefaults() const override
   {
      static auto dflts = RDrawableAttributes::Map_t().AddString("rgb","0,0,0").AddDouble("a",1.);
      return dflts;
   }

public:

   using RGB_t = std::array<int, 3>;


   using RAttributesVisitor::RAttributesVisitor;

   RColor(int r, int g, int b) : RColor()
   {
      SetRGB(r,g,b);
   }

   RColor(int r, int g, int b, double alfa) : RColor()
   {
      SetRGB(r,g,b);
      SetAlfa(alfa);
   }

   RColor(const RGB_t &rgb) : RColor()
   {
      SetRGB(rgb[0],rgb[1],rgb[2]);
   }


   std::string GetRGB() const { return GetValue<std::string>("rgb"); }
   RColor &SetRGB(const std::string &_rgb) { SetValue("rgb", _rgb); return *this; }
   RColor &SetRGB(int r, int g, int b) { return SetRGB(std::to_string(r) + "," + std::to_string(g) + "," + std::to_string(b)); }

   double GetAlfa() const { return GetValue<double>("a"); }
   bool HasAlfa() const { return HasValue("a"); }
   RColor &SetAlfa(double _alfa) { SetValue("a", _alfa); return *this; }

   std::string AsSVG() const
   {
      auto rgb = GetRGB();
      if (HasAlfa())
         return std::string("rgba(") + rgb + "," + std::to_string(GetAlfa()) + ")";
       return std::string("rgb(") + rgb + ")";
   }

   static constexpr RGB_t kRed{{255, 0, 0}};
   static constexpr RGB_t kGreen{{0, 255, 0}};
   static constexpr RGB_t kBlue{{0, 0, 255}};
   static constexpr RGB_t kWhite{{255, 255, 255}};
   static constexpr RGB_t kBlack{{0, 0, 0}};

};


} // namespace Experimental
} // namespace ROOT

#endif
