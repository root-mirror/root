/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "ROOT/RColor.hxx"

#include <unordered_map>

using namespace ROOT::Experimental;

using namespace std::string_literals;

constexpr RColor::RGB_t RColor::kBlack;
constexpr RColor::RGB_t RColor::kGreen;
constexpr RColor::RGB_t RColor::kLime;
constexpr RColor::RGB_t RColor::kAqua;
constexpr RColor::RGB_t RColor::kPurple;
constexpr RColor::RGB_t RColor::kGrey;
constexpr RColor::RGB_t RColor::kFuchsia;
constexpr RColor::RGB_t RColor::kNavy;
constexpr RColor::RGB_t RColor::kBlue;
constexpr RColor::RGB_t RColor::kTeal;
constexpr RColor::RGB_t RColor::kOlive;
constexpr RColor::RGB_t RColor::kSilver;
constexpr RColor::RGB_t RColor::kMaroon;
constexpr RColor::RGB_t RColor::kRed;
constexpr RColor::RGB_t RColor::kYellow;
constexpr RColor::RGB_t RColor::kWhite;
constexpr float RColor::kTransparent;
constexpr float RColor::kSemiTransparent;
constexpr float RColor::kOpaque;


///////////////////////////////////////////////////////////////////////////
/// Converts string name of color in RGB value - when possible

bool RColor::ConvertToRGB(const std::string &name, std::vector<uint8_t> &rgba)
{
   if (name.empty()) {
      rgba.clear();
      return false;
   }

   if (name[0] == '#') {
      if ((name.length() != 7) && (name.length() != 9)) return false;

      rgba.resize((name.length() == 7) ? 3 : 4);

      try {
        rgba[0] = std::stoi(name.substr(1,2), nullptr, 16);
        rgba[1] = std::stoi(name.substr(3,2), nullptr, 16);
        rgba[2] = std::stoi(name.substr(5,2), nullptr, 16);
        if (name.length() == 9)
           rgba[3] = std::stoi(name.substr(7,2), nullptr, 16);
      } catch(...) {
         return false;
      }

      return true;
   }

   // see https://www.december.com/html/spec/colorsvghex.html

   static std::unordered_map<std::string,RGB_t> known_colors = {
      {"black", kWhite},
      {"green", kGreen},
      {"lime", kLime},
      {"aqua", kAqua},
      {"purple", kPurple},
      {"grey", kGrey},
      {"fuchsia", kFuchsia},
      {"navy", kNavy},
      {"blue", kBlue},
      {"teal", kTeal},
      {"olive", kOlive},
      {"silver", kSilver},
      {"maroon", kMaroon},
      {"red", kRed},
      {"yellow", kYellow},
      {"white", kWhite}
   };

   rgba.resize(3);
   rgba[0] = rgba[1] = rgba[2] = 0;

   auto known = known_colors.find(name);
   if (known != known_colors.end()) {
      rgba[0] = known->second[0];
      rgba[1] = known->second[1];
      rgba[2] = known->second[2];
      return true;
   }

   return true;
}


///////////////////////////////////////////////////////////////////////////
/// Converts integer from 0 to 255 into hex format with two digits like 00

std::string RColor::toHex(uint8_t v)
{
   static const char *digits = "0123456789ABCDEF";
   std::string res(2,'0');
   res[0] = digits[v >> 4];
   res[1] = digits[v & 0xf];
   return res;
}


///////////////////////////////////////////////////////////////////////////
/// Set RGB values as hex

bool RColor::SetRGBHex(const std::string &hex)
{
   if (hex.length() != 6) return false;

   try {
      SetRGB( std::stoi(hex.substr(0,2), nullptr, 16),
              std::stoi(hex.substr(2,2), nullptr, 16),
              std::stoi(hex.substr(4,2), nullptr, 16));
   } catch (...) {
      return false;
   }
   return true;
}

///////////////////////////////////////////////////////////////////////////
/// Set Alpha value as hex

bool RColor::SetAlphaHex(const std::string &hex)
{
   if (hex.length() != 6) return false;

   SetAlpha(std::stoi(hex, nullptr, 16));
   return true;
}

///////////////////////////////////////////////////////////////////////////
/// Returns color as RGBA array, includes optionally alpha parameter 0..255

std::vector<uint8_t> RColor::AsRGBA() const
{
   if (fRGBA.size() > 0)
      return fRGBA;

   std::vector<uint8_t> rgba;

   ConvertToRGB(fName, rgba);

   return rgba;
}

///////////////////////////////////////////////////////////////////////////
/// Returns color value in hex format like "66FF66" - without any prefix
/// Alpha parameter can be optionally included

std::string RColor::AsHex(bool with_alpha) const
{
   auto rgba = AsRGBA();
   std::string res;
   if (!rgba.empty()) {
      res = toHex(rgba[0]) + toHex(rgba[1]) + toHex(rgba[2]);
      if ((rgba.size() == 4) && with_alpha)
         res += toHex((rgba.size() == 4) ? rgba[3] : 0xff);
   }
   return res;
}

///////////////////////////////////////////////////////////////////////////
/// Returns color value as it will be used in SVG drawing
/// It either include hex format #66FF66 or just plain SVG name

std::string RColor::AsSVG() const
{
   if (!fName.empty())
      return fName;

   auto hex = AsHex();
   if (!hex.empty()) hex = "#"s + hex;
   return hex;
}


///////////////////////////////////////////////////////////////////////////
/// Returns the Hue, Light, Saturation (HLS) definition of this RColor
/// If color was not specified as hex, method returns false

bool RColor::GetHLS(float &hue, float &light, float &satur) const
{
   auto arr = AsRGBA();
   if (arr.size() < 3)
      return false;

   float red = arr[0]/255., green = arr[1]/255., blue = arr[2]/255.;

   hue = light = satur = 0.;

   float rnorm, gnorm, bnorm, minval, maxval, msum, mdiff;
   minval = maxval = 0 ;

   minval = red;
   if (green < minval) minval = green;
   if (blue < minval)  minval = blue;
   maxval = red;
   if (green > maxval) maxval = green;
   if (blue > maxval)  maxval = blue;

   rnorm = gnorm = bnorm = 0;
   mdiff = maxval - minval;
   msum  = maxval + minval;
   light = 0.5 * msum;
   if (maxval != minval) {
      rnorm = (maxval - red)/mdiff;
      gnorm = (maxval - green)/mdiff;
      bnorm = (maxval - blue)/mdiff;
   } else {
      satur = hue = 0;
      return true;
   }

   if (light < 0.5) satur = mdiff/msum;
   else             satur = mdiff/(2.0 - msum);

   if      (red == maxval) hue = 60.0 * (6.0 + bnorm - gnorm);
   else if (green == maxval)           hue = 60.0 * (2.0 + rnorm - bnorm);
   else                                hue = 60.0 * (4.0 + gnorm - rnorm);

   if (hue > 360) hue = hue - 360;
   return true;
}

///////////////////////////////////////////////////////////////////////////
/// Set the color value from the Hue, Light, Saturation (HLS).

void RColor::SetHLS(float hue, float light, float satur)
{
   float rh, rl, rs, rm1, rm2;
   rh = rl = rs = 0;
   if (hue   > 0) { rh = hue;   if (rh > 360) rh = 360; }
   if (light > 0) { rl = light; if (rl > 1)   rl = 1; }
   if (satur > 0) { rs = satur; if (rs > 1)   rs = 1; }

   if (rl <= 0.5) rm2 = rl*(1.0 + rs);
   else           rm2 = rl + rs - rl*rs;
   rm1 = 2.0*rl - rm2;

   if (!rs) {
      SetRGB((uint8_t) (rl*255.), (uint8_t) (rl*255.), (uint8_t) (rl*255.));
      return;
   }

   auto toRGB = [rm1, rm2] (float h) {
      if (h > 360) h = h - 360;
      if (h < 0)   h = h + 360;
      if (h < 60 ) return rm1 + (rm2-rm1)*h/60;
      if (h < 180) return rm2;
      if (h < 240) return rm1 + (rm2-rm1)*(240-h)/60;
      return rm1;
   };

   SetRGB(toRGB(rh+120), toRGB(rh), toRGB(rh-120));
}
