/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RAttrColor
#define ROOT7_RAttrColor

#include <ROOT/RAttrBase.hxx>

#include <ROOT/RColor.hxx>

#include <array>

namespace ROOT {
namespace Experimental {

// TODO: see also imagemagick's C++ interface for RColor operations!
// https://www.imagemagick.org/api/magick++-classes.php

/** \class RAttrColor
\ingroup GpadROOT7
\brief Access RColor from drawable attributes
\author Sergey Linev <S.Linev@gsi.de>
\date 2020-03-09
\warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is
welcome!
*/

class RAttrColor : public RAttrBase {

   R__ATTR_CLASS(RAttrColor, "color");

protected:

   /** Provide possible attributes names and values */
   RAttrMap CollectDefaults() const override { return RAttrMap().AddString("rgb", "").AddString("a", "").AddString("name", "").AddBool("auto", false); }

   /** Set color as plain SVG name like "white" or "lightblue". Clears RGB component before */
   void SetName(const std::string &_name) { SetValue("name", _name); }

   /** Returns color as plain SVG name like "white" or "lightblue" */
   std::string GetName() const { return GetValue<std::string>("name"); }

   /** Clear color plain SVG name (if any) */
   void ClearName() { ClearValue("name"); }

   /** Returns true if color name was  specified */
   bool HasName() const { return HasValue("name"); }

   /** Set color as hex string like 00FF00 */
   void SetHex(const std::string &_hex) { SetValue("rgb", _hex); }

   /** Remove color hex value */
   void ClearHex() { ClearValue("rgb"); }

   /** Returns true if color hex value was  specified */
   bool HasHex() const { return HasValue("rgb"); }

   /** Returns color alpha (opacity) as hex string like FF. Default is empty */
   std::string GetAlphaHex() const { return GetValue<std::string>("a"); }

   /** Returns true if color alpha (opacity) was specified */
   bool HasAlpha() const { return HasValue("a"); }

   /** Set color alpha (opacity) value - from 0 to 1 */
   void SetAlpha(float alpha) { return SetAlphaHex(RColor::toHex((uint8_t)(alpha * 255))); }

   /** Set color alpha (opacity) value as hex */
   void SetAlphaHex(const std::string &val) { SetValue("a", val); }

   /** Remove color alpha value */
   void ClearAlpha() { ClearValue("a"); }

public:
   /** Set r/g/b components of color as hex code, default for the color */
   RAttrColor &SetColor(const RColor &col)
   {
      Clear();

      if (!col.GetName().empty()) {
         SetName(col.GetName());
      } else {
         auto rgba = col.GetRGBA();
         if (rgba.size() > 2) SetHex(RColor::toHex(rgba[0]) + RColor::toHex(rgba[1]) + RColor::toHex(rgba[2]));
         if (rgba.size() == 4) SetAlphaHex(RColor::toHex(rgba[3]));
      }

      return *this;
   }

   /** Extract RColor for given attribute */
   RColor GetColor() const
   {
      RColor res;
      if (HasName())
         res.SetName(GetName());
      else if (HasHex()) {
         res.SetRGBHex(GetHex());
         if (HasAlpha())
            res.SetAlphaHex(GetAlphaHex());
      }
      return res;
   }

   /** Remove all values which can correspond to RColor value */
   void Clear()
   {
      ClearHex();
      ClearAlpha();
      ClearName();
      ClearAuto();
   }

   /** Return color as hex string like 00FF00 */
   std::string GetHex() const { return GetValue<std::string>("rgb"); }


   /** Returns true if color should get auto value when primitive drawing is performed */
   bool IsAuto() const { return GetValue<bool>("auto"); }

   /** Set automatic mode for RAttrColor, will be assigned before primitive painted on the canvas */
   void SetAuto(bool on = true) { SetValue("auto", on); }

   /** Clear auto flag of the RAttrColor */
   void ClearAuto() { ClearValue("auto"); }


   friend bool operator==(const RAttrColor &lhs, const RAttrColor &rhs) { return lhs.GetColor() == rhs.GetColor(); }


   RAttrColor &operator=(const RColor &col)
   {
      SetColor(col);
      return *this;
   }

};

} // namespace Experimental
} // namespace ROOT

#endif
