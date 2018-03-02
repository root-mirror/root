/// \file TColor.cxx
/// \ingroup Gpad ROOT7
/// \author Axel Naumann <axel@cern.ch>
/// \date 2017-09-27
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2017, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "ROOT/TColor.hxx"

#include <ROOT/TLogger.hxx>

#include <exception>

using namespace ROOT::Experimental;

// TColor constexpr values:
constexpr TColor::Alpha TColor::kOpaque;
constexpr TColor::Alpha TColor::kTransparent;
constexpr TColor::RGBA TColor::kRed;
constexpr TColor::RGBA TColor::kGreen;
constexpr TColor::RGBA TColor::kBlue;
constexpr TColor::RGBA TColor::kWhite;
constexpr TColor::RGBA TColor::kBlack;
constexpr TColor::RGBA TColor::kInvisible;
constexpr TColor::AutoTag TColor::kAuto;


float TColor::GetPaletteOrdinal() const
{
   if (fKind != EKind::kPalettePos)
      throw std::runtime_error("This color does not represent a palette ordinal!");
   return fRedOrPalettePos;
}

bool TColor::AssertNotPalettePos() const
{
   if (fKind == EKind::kPalettePos) {
      throw std::runtime_error("This color does not represent a palette ordinal!");
      return false;
   }
   return true;
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize an attribute `val` from a string value.
/// Colors can be specified as RGBA (red green blue alpha) or RRGGBBAA:
///     %fa7f %ffa07bff # hash introduces a comment!
/// For all predefined colors in TColor, colors can be specified as name without leading 'k', e.g. `red` for
/// `TColor::kRed`.
/// Prints an error and returns `TColor::kBlack` if the attribute string cannot be parsed or if the attribute has no
/// entry in `fAttrs`.
///
///\param[in] name - the attribute name (for diagnostic purposes).
///\param[in] strval - the attribute value as a string.
///\param[out] val - the value to be initialized.

void ROOT::Experimental::InitializeAttrFromString(const std::string &name, const std::string &strval, TColor& val)
{
   if (strval.empty())
      return;

   if (strval[0] == '#') {
      auto rgbalen = strval.length() - 1;
      if (rgbalen != 3 && rgbalen != 4 && rgbalen != 6 && rgbalen != 8) {
         R__ERROR_HERE("Graf2d") << "Invalid value for TColor default style " << name
            << " with value \"" << strval
            << "\": expect '#' followed by 3, 4, 6 or 8 hex digits (#rgb, #rgba, #rrggbbaa or #rrggbb).";
         return;
      }
      std::size_t pos;
      long long rgbaLL = std::stoll(strval.substr(1), &pos, /*base*/ 16);
      if (pos != 3 && pos != 4 && pos != 6 && pos != 8) {
         R__ERROR_HERE("Graf2d") << "Invalid value while parsing default style value for TColor " << name
            << " with value \"" << strval
            << "\": expect '#' followed by 3, 4, 6 or 8 hex digits (#rgb, #rgba, #rrggbbaa or #rrggbb).";
         return;
      }
      if (pos != rgbalen) {
         R__WARNING_HERE("Graf2d") << "Leftover characters while parsing default style value for TColor " << name
            << " with value \"" << strval << "\", remainder: \"" << strval.substr(pos - 1) << "\"";
         return;
      }
      float rgba[4] = {0};
      // #rrggbb[aa] has 8 bits per channel, #rgb[a] has 4.
      const int bitsPerChannel = (pos > 4) ? 8 : 4;
      const int bitMask = (1 << bitsPerChannel) - 1;
      const float bitMaskFloat = static_cast<float>(bitMask);
      for (auto& channel: rgba) {
         channel = (rgbaLL & bitMask) / bitMaskFloat;
         rgbaLL >>= bitsPerChannel;
      }
      if (pos == 3 || pos == 6) {
         // no alpha, set it to 255.
         rgba[3] = 1.;
      }
      val = TColor(rgba[0]);
   }
}
