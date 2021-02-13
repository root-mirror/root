/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RColor
#define ROOT7_RColor

#include <cstdint>
#include <vector>
#include <string>
#include <array>
#include <DllImport.h>

namespace ROOT {
namespace Experimental {

class RAttrColor;

// TODO: see also imagemagick's C++ interface for RColor operations!
// https://www.imagemagick.org/api/magick++-classes.php

/** \class RColor
\ingroup GpadROOT7
\brief The color class
\author Axel Naumann <axel@cern.ch>
\author Sergey Linev <S.Linev@gsi.de>
\date 2017-09-26
\warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
*/

class RColor {

   friend class RAttrColor;

   using RGB_t = std::array<uint8_t, 3>;

private:

   std::vector<uint8_t> fRGBA;  ///<  RGB + Alpha
   std::string fName;           ///<  name of color - if any

   static std::string toHex(uint8_t v);

   static bool ConvertToRGB(const std::string &name, std::vector<uint8_t> &rgba);

   bool SetRGBHex(const std::string &hex);
   bool SetAlphaHex(const std::string &hex);

public:

   RColor() = default;

   /** Construct color with provided r,g,b values */
   RColor(uint8_t r, uint8_t g, uint8_t b) { SetRGB(r, g, b); }

   /** Construct color with provided r,g,b and alpha values */
   RColor(uint8_t r, uint8_t g, uint8_t b, float alpha)
   {
      SetRGBA(r, g, b, alpha);
   }

   /** Construct color with provided RGB_t value */
   RColor(const RGB_t &rgb) { SetRGB(rgb[0], rgb[1], rgb[2]); };

   /** Construct color with provided name */
   RColor(const std::string &name) { SetName(name); };

   /** Set r/g/b components of color */
   void SetRGB(const RGB_t &rgb) { SetRGB(rgb[0], rgb[1], rgb[2]); }

   /** Set r/g/b components of color */
   void SetRGB(uint8_t r, uint8_t g, uint8_t b)
   {
      fName.clear();
      if (fRGBA.size() < 3)
         fRGBA.resize(3);
      fRGBA[0] = r;
      fRGBA[1] = g;
      fRGBA[2] = b;
   }

   /** Set r/g/b/a components of color */
   void SetRGBA(uint8_t r, uint8_t g, uint8_t b, uint8_t alpha)
   {
      fName.clear();
      fRGBA.resize(4);
      SetRGB(r,g,b);
      SetAlpha(alpha);
   }

   /** Set r/g/b/a components of color */
   void SetRGBA(uint8_t r, uint8_t g, uint8_t b, float alpha)
   {
      fName.clear();
      fRGBA.resize(4);
      SetRGB(r,g,b);
      SetAlphaFloat(alpha);
   }

   /** Set alpha as float value from range 0..1 */
   void SetAlphaFloat(float alpha)
   {
      if (alpha <= 0.)
         SetAlpha(0);
      else if (alpha >= 1.)
         SetAlpha(255);
      else
         SetAlpha((uint8_t)(alpha * 255));
   }

   /** Set alpha as value from range 0..255 */
   void SetAlpha(uint8_t alpha)
   {
      if (fRGBA.empty()) {
         ConvertToRGB(fName, fRGBA);
         fName.clear();
      }

      if (fRGBA.size() < 4)
         fRGBA.resize(4);

      fRGBA[3] = alpha;
   }

   /** Returns true if color alpha (opacity) was specified */
   bool HasAlpha() const { return fRGBA.size() == 4; }

   /** Returns true if no color is specified */
   bool IsEmpty() const { return fName.empty() && (fRGBA.size() == 0); }

   /** Returns color as RGBA array - when exists */
   const std::vector<uint8_t> &GetRGBA() const { return fRGBA; }

   /** Returns color as RGBA array, trying also convert color name into RGBA value */
   std::vector<uint8_t> AsRGBA() const;

   /** Returns red color component 0..255 */
   uint8_t GetRed() const
   {
      if (fRGBA.size() > 2)
         return fRGBA[0];

      std::vector<uint8_t> rgb;
      return ConvertToRGB(fName, rgb) ? rgb[0] : 0;
   }

   /** Returns green color component 0..255 */
   uint8_t GetGreen() const
   {
      if (fRGBA.size() > 2)
         return fRGBA[1];

      std::vector<uint8_t> rgb;
      return ConvertToRGB(fName, rgb) ? rgb[1] : 0;
   }

   /** Returns blue color component 0..255 */
   uint8_t GetBlue() const
   {
      if (fRGBA.size() > 2)
         return fRGBA[2];

      std::vector<uint8_t> rgb;
      return ConvertToRGB(fName, rgb) ? rgb[2] : 0;
   }

   /** Returns color alpha (opacity) as uint8_t 0..255 */
   uint8_t GetAlpha() const
   {
      if (fRGBA.size() > 0)
         return fRGBA.size() == 4 ? fRGBA[3] : 255;

      std::vector<uint8_t> rgba;
      if (ConvertToRGB(fName, rgba) && (rgba.size() == 3))
         return rgba[3];

      return 255;
   }

   /** Returns color alpha (opacity) as float from 0..1 */
   float GetAlphaFloat() const
   {
      return GetAlpha() / 255.;
   }

   /** Set color as plain SVG name like "white" or "lightblue". Clears RGB component before */
   RColor &SetName(const std::string &name)
   {
      fRGBA.clear();
      fName = name;
      return *this;
   }

   /** Returns color as plain SVG name like "white" or "lightblue" */
   const std::string &GetName() const { return fName; }

   /** Return the Hue, Light, Saturation (HLS) definition of this RColor */
   bool GetHLS(float &hue, float &light, float &satur) const;

   /** Set the Red Green and Blue (RGB) values from the Hue, Light, Saturation (HLS). */
   void SetHLS(float hue, float light, float satur);

   std::string AsHex(bool with_alpha = false) const;
   std::string AsSVG() const;

   void Clear()
   {
      fRGBA.clear();
      fName.clear();
   }

   R__DLLEXPORT static constexpr RGB_t kBlack{{0, 0, 0}};
   R__DLLEXPORT static constexpr RGB_t kGreen{{0, 0x80, 0}};
   R__DLLEXPORT static constexpr RGB_t kLime{{0, 0xFF, 0}};
   R__DLLEXPORT static constexpr RGB_t kAqua{{0, 0xFF, 0xFF}};
   R__DLLEXPORT static constexpr RGB_t kPurple{{0x80, 0, 0x80}};
   R__DLLEXPORT static constexpr RGB_t kGrey{{0x80, 0x80, 0x80}};
   R__DLLEXPORT static constexpr RGB_t kFuchsia{{0xFF, 0, 0xFF}};
   R__DLLEXPORT static constexpr RGB_t kNavy{{0, 0, 0x80}};
   R__DLLEXPORT static constexpr RGB_t kBlue{{0, 0, 0xff}};
   R__DLLEXPORT static constexpr RGB_t kTeal{{0, 0x80, 0x80}};
   R__DLLEXPORT static constexpr RGB_t kOlive{{0x80, 0x80, 0}};
   R__DLLEXPORT static constexpr RGB_t kSilver{{0xc0, 0xc0, 0xc0}};
   R__DLLEXPORT static constexpr RGB_t kMaroon{{0x80, 0, 0}};
   R__DLLEXPORT static constexpr RGB_t kRed{{0xff, 0, 0}};
   R__DLLEXPORT static constexpr RGB_t kYellow{{0xff, 0xff, 0}};
   R__DLLEXPORT static constexpr RGB_t kWhite{{0xff, 0xff, 0xff}};
   R__DLLEXPORT static constexpr float kTransparent{0.};
   R__DLLEXPORT static constexpr float kSemiTransparent{0.5};
   R__DLLEXPORT static constexpr float kOpaque{1.};

   friend bool operator==(const RColor &lhs, const RColor &rhs)
   {
      if ((lhs.fName == rhs.fName) && (lhs.fRGBA == rhs.fRGBA)) return true;

      auto l = lhs.AsRGBA();
      auto r = rhs.AsRGBA();

      return !l.empty() && (l==r);
   }
};

} // namespace Experimental
} // namespace ROOT

#endif
