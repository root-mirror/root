/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RAttrLine
#define ROOT7_RAttrLine

#include <ROOT/RAttrBase.hxx>
#include <ROOT/RAttrColor.hxx>
#include <ROOT/RAttrValue.hxx>

namespace ROOT {
namespace Experimental {

/** \class RAttrLine
\ingroup GpadROOT7
\author Axel Naumann <axel@cern.ch>
\date 2018-10-12
\brief Drawing line attributes for different objects.
\warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
*/

class RAttrLine : public RAttrBase {

   RAttrColor          fColor{this, "color"};      ///<! line color
   RAttrValue<double>  fWidth{this, "width", 1.};  ///<! line width
   RAttrValue<int>     fStyle{this, "style", 1};   ///<! line style

   R__ATTR_CLASS(RAttrLine, "line");

   ///The width of the line.
   RAttrLine &SetWidth(double width) { fWidth = width; return *this; }
   double GetWidth() const { return fWidth; }

   ///The style of the line.
   RAttrLine &SetStyle(int style) { fStyle = style; return *this; }
   int GetStyle() const { return fStyle; }

   ///The color of the line.
   RAttrLine &SetColor(const RColor &color) { fColor = color; return *this; }
   RColor GetColor() const { return fColor.GetColor(); }
   RAttrColor &AttrColor() { return fColor; }

};

} // namespace Experimental
} // namespace ROOT

#endif
