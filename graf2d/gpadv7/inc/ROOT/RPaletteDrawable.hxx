/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RPaletteDrawable
#define ROOT7_RPaletteDrawable

#include <ROOT/RDrawable.hxx>
#include <ROOT/RAttrValue.hxx>
#include <ROOT/RPadPos.hxx>
#include <ROOT/RPalette.hxx>

#include <memory>
#include <string>

namespace ROOT {
namespace Experimental {

/** \class ROOT::Experimental::RPaletteDrawable
\ingroup GrafROOT7
\brief A color palette draw near the frame.
\author Sergey Linev <s.linev@gsi.de>
\date 2020-03-05
\warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
*/


class RPaletteDrawable final : public RDrawable {

   RPalette                fPalette;                              ///<  color palette to draw
   RAttrValue<bool>        fVisible{this, "visible", true};       ///<! visibility flag
   RAttrValue<RPadLength>  fMargin{this, "margin", 0.02_normal};  ///<! margin
   RAttrValue<RPadLength>  fSize{this, "size", 0.05_normal};      ///<! margin

protected:

   bool IsFrameRequired() const final { return true; }

   RPaletteDrawable() : RDrawable("palette") {}

public:

   RPaletteDrawable(const RPalette &palette) : RPaletteDrawable() { fPalette = palette; }
   RPaletteDrawable(const RPalette &palette, bool visible) : RPaletteDrawable() { fPalette = palette; SetVisible(visible); }
   const RPalette &GetPalette() const { return fPalette; }

   RPaletteDrawable &SetVisible(bool on = true) { fVisible = on; return *this; }
   bool GetVisible() const { return fVisible; }

   RPaletteDrawable &SetMargin(const RPadLength &pos) { fMargin = pos; return *this; }
   RPadLength GetMargin() const { return fMargin; }

   RPaletteDrawable &SetSize(const RPadLength &sz) { fSize = sz; return *this; }
   RPadLength GetSize() const { return fSize; }
};

//inline auto GetDrawable(const RPalette &palette)
//{
//   return std::make_shared<RPaletteDrawable>(palette);
//}


} // namespace Experimental
} // namespace ROOT

#endif
