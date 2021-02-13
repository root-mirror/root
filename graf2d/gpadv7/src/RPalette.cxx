/*************************************************************************
 * Copyright (C) 1995-2017, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "ROOT/RPalette.hxx"

#include "ROOT/RLogger.hxx"
#include "ROOT/RAttrBase.hxx" // for GPadLog()

#include <algorithm>
#include <cmath>
#include <exception>
#include <unordered_map>

using namespace ROOT::Experimental;

RPalette::RPalette(bool interpolate, bool knownNormalized, const std::vector<RPalette::OrdinalAndColor> &points)
   : fColors(points), fInterpolate(interpolate), fNormalized(knownNormalized)
{
   if (points.size() < 2)
      throw std::runtime_error("Must have at least two points to build a palette!");

   std::sort(fColors.begin(), fColors.end());

   if (!knownNormalized) {
      // Is this a normalized palette? I.e. are the first and last ordinals 0 and 1?
      double high = fColors.back().fOrdinal;
      double low = fColors.front().fOrdinal;
      double prec = (high - low) * 1e-8;

      auto reasonablyEqual = [&](double val, double expected) -> bool { return std::fabs(val - expected) < prec; };
      fNormalized = reasonablyEqual(low, 0.) && reasonablyEqual(high, 1.);
   }
}

namespace {
static std::vector<RPalette::OrdinalAndColor> AddOrdinals(const std::vector<RColor> &points)
{
   std::vector<RPalette::OrdinalAndColor> ret(points.size());
   auto addOneOrdinal = [&](const RColor &col) -> RPalette::OrdinalAndColor {
      return {1. / (points.size() - 1) * (&col - points.data()), col};
   };
   std::transform(points.begin(), points.end(), ret.begin(), addOneOrdinal);
   return ret;
}
} // unnamed namespace

RPalette::RPalette(bool interpolate, const std::vector<RColor> &points)
   : RPalette(interpolate, true, AddOrdinals(points))
{}

RColor RPalette::GetColor(double ordinal)
{
   if (fColors.size() == 0)
      return RColor();

   if (fColors.size() == 1)
      return fColors.front().fColor;

   constexpr float epsilon = 1e-8;
   if (ordinal < fColors.front().fOrdinal + epsilon)
      return fColors.front().fColor;

   if (ordinal > fColors.back().fOrdinal - epsilon)
      return fColors.back().fColor;

   auto iColor2 = std::lower_bound(fColors.begin(), fColors.end(), ordinal);
   auto iColor1 = iColor2 - 1;

   auto diff1 = ordinal - iColor1->fOrdinal;
   auto diff2 = iColor2->fOrdinal - ordinal;

   if ((diff1 < -epsilon) || (diff2 < -epsilon)) {
      R__LOG_ERROR(GPadLog()) << "Wrong palette settings";
      return fColors.back().fColor;
   }

   if (diff1 < epsilon)
      return iColor1->fColor;

   if (diff2 < epsilon)
      return iColor2->fColor;

   if (IsGradient()) {
      auto dist = diff1 + diff2;
      auto rgba1 = iColor1->fColor.AsRGBA();
      auto rgba2 = iColor2->fColor.AsRGBA();
      if ((dist > epsilon) && (rgba1.size() > 2) && (rgba2.size() > 2)) {
         if (rgba1.size() == 4)
            rgba2.resize(4, 0xff);
         else if (rgba2.size() == 4)
            rgba1.resize(4, 0xff);

         for (unsigned i = 0; i < rgba1.size(); ++i)
            rgba1[i] = (uint8_t) std::lround( (diff2*rgba1[i] + diff1*rgba2[i]) / dist);

         RColor res;
         res.SetRGB(rgba1[0], rgba1[1], rgba1[2]);
         if (rgba1.size() == 4)
            res.SetAlpha(rgba1[3]);

         return res;
      }

      R__LOG_ERROR(GPadLog()) << "Fail to interpolate color";
   }

   return (diff2 < diff1) ? iColor2->fColor : iColor1->fColor;
}

namespace {
using GlobalPalettes_t = std::unordered_map<std::string, RPalette>;
static GlobalPalettes_t CreateDefaultPalettes()
{
   GlobalPalettes_t ret;
   ret["default"] = RPalette({RColor::kRed, RColor::kBlue});
   ret["bw"] = RPalette({RColor::kBlack, RColor::kWhite});
   ret["bird"] = RPalette({RColor(53,42, 135), RColor(15,92,221),   RColor(20,129,214),
                           RColor(6,164,202),  RColor(46,183, 164), RColor(135,191,119),
                           RColor(209,187,89), RColor(254,200,50),  RColor(249,251,14)});
   ret["rainbow"] = RPalette({RColor(0,0,99), RColor(5,48,142), RColor(15,124,198),
                              RColor(35,192,201), RColor(102,206,90), RColor(196,226,22),
                              RColor(208,97,13), RColor(199,16,8), RColor(110,0,2)});
   return ret;
}

static GlobalPalettes_t &GetGlobalPalettes()
{
   static GlobalPalettes_t globalPalettes = CreateDefaultPalettes();
   return globalPalettes;
}
} // unnamed namespace

void RPalette::RegisterPalette(std::string_view name, const RPalette &palette)
{
   GetGlobalPalettes()[std::string(name)] = palette;
}

const RPalette &RPalette::GetPalette(std::string_view name)
{
   static const RPalette sNoPaletteWithThatName;
   if (name.empty()) name = "bird";
   auto iGlobalPalette = GetGlobalPalettes().find(std::string(name));
   if (iGlobalPalette == GetGlobalPalettes().end())
      return sNoPaletteWithThatName;
   return iGlobalPalette->second;
}
