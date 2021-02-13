/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RAttrMargins
#define ROOT7_RAttrMargins

#include <ROOT/RAttrBase.hxx>
#include <ROOT/RAttrValue.hxx>
#include <ROOT/RPadLength.hxx>

namespace ROOT {
namespace Experimental {

/** \class RAttrMargins
\ingroup GpadROOT7
\author Sergey Linev <s.linev@gsi.de>
\date 2020-02-20
\brief A margins attributes. Only relative and pixel coordinates are allowed
\warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
*/

class RAttrMargins : public RAttrBase {

   RAttrValue<RPadLength>   fLeft{this, "left", 0._normal};
   RAttrValue<RPadLength>   fRight{this, "right", 0._normal};
   RAttrValue<RPadLength>   fTop{this, "top", 0._normal};
   RAttrValue<RPadLength>   fBottom{this, "bottom", 0._normal};

   R__ATTR_CLASS(RAttrMargins, "margin");

public:

   RAttrMargins &SetLeft(const RPadLength &len) { fLeft = len; return *this; }
   RPadLength GetLeft() const { return fLeft; }

   RAttrMargins &SetRight(const RPadLength &len) { fRight = len; return *this; }
   RPadLength GetRight() const { return fRight; }

   RAttrMargins &SetTop(const RPadLength &len) { fTop = len; return *this; }
   RPadLength GetTop() const { return fTop; }

   RAttrMargins &SetBottom(const RPadLength &len) { fBottom = len; return *this; }
   RPadLength GetBottom() const { return fBottom; }
};

} // namespace Experimental
} // namespace ROOT

#endif
