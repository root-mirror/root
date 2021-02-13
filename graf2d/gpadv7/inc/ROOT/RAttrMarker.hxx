/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RAttrMarker
#define ROOT7_RAttrMarker

#include <ROOT/RAttrBase.hxx>
#include <ROOT/RAttrColor.hxx>
#include <ROOT/RAttrValue.hxx>

namespace ROOT {
namespace Experimental {

/** \class RAttrMarker
\ingroup GpadROOT7
\authors Axel Naumann <axel@cern.ch> Sergey Linev <s.linev@gsi.de>
\date 2018-10-12
\brief A marker attributes.
\warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
*/

class RAttrMarker : public RAttrBase {

   RAttrColor           fColor{this, "color"};      ///<! marker color
   RAttrValue<double>   fSize{this, "size", 1.};    ///<! marker size
   RAttrValue<int>      fStyle{this, "style", 1};   ///<! marker style

   R__ATTR_CLASS(RAttrMarker, "marker");

   RAttrMarker &SetColor(const RColor &color) { fColor = color; return *this; }
   RColor GetColor() const { return fColor.GetColor(); }
   RAttrColor &AttrColor() { return fColor; }

   /// The size of the marker.
   RAttrMarker &SetSize(double size) { fSize = size; return *this; }
   double GetSize() const { return fSize; }

   /// The style of the marker.
   RAttrMarker &SetStyle(int style) { fStyle = style; return *this; }
   int GetStyle() const { return fStyle; }

};

} // namespace Experimental
} // namespace ROOT

#endif
