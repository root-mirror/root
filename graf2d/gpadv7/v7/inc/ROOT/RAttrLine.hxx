/// \file ROOT/RAttrLine.hxx
/// \ingroup Graf ROOT7
/// \author Axel Naumann <axel@cern.ch>
/// \date 2018-10-12
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RAttrLine
#define ROOT7_RAttrLine

#include <ROOT/RColor.hxx>
#include <ROOT/RDrawingAttr.hxx>

namespace ROOT {
namespace Experimental {

/** class ROOT::Experimental::RAttrLine
 Drawing attributes for RLine.
 */
class RAttrLine: public RDrawingAttrBase {
public:
   using RDrawingAttrBase::RDrawingAttrBase;

   /// The color of the line.
   RAttrLine &SetColor(const RColor& col) { Set("color", col); return *this; }
   RColor GetColor() const { return Get<RColor>("color"); }

   ///The width of the line.
   RAttrLine &SetWidth(float width) { Set("width", width); return *this; }
   float GetWidth() const { return Get<float>("width"); }

   ///The style of the line.
   RAttrLine &SetStyle(int style) { Set("style", style); return *this; }
   int GetStyle() const { return Get<int>("style"); }
};

} // namespace Experimental
} // namespace ROOT

#endif
