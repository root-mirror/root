/*************************************************************************
 * Copyright (C) 1995-2017, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RMarker
#define ROOT7_RMarker

#include <ROOT/RDrawable.hxx>
#include <ROOT/RAttrMarker.hxx>
#include <ROOT/RPadPos.hxx>

#include <initializer_list>

namespace ROOT {
namespace Experimental {

/** \class RMarker
\ingroup GrafROOT7
\brief A simple marker.
\author Olivier Couet <Olivier.Couet@cern.ch>
\date 2017-10-16
\warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
*/

class RMarker : public RDrawable {

   RPadPos fP;                              ///< position
   RAttrMarker fMarkerAttr{this, "marker"}; ///<! marker attributes

public:
   RMarker() : RDrawable("marker") {}

   RMarker(const RPadPos &p) : RMarker() { fP = p; }

   RMarker &SetP(const RPadPos &p)
   {
      fP = p;
      return *this;
   }
   const RPadPos &GetP() const { return fP; }

   const RAttrMarker &GetAttrMarker() const { return fMarkerAttr; }
   RMarker &SetAttrMarker(const RAttrMarker &attr)
   {
      fMarkerAttr = attr;
      return *this;
   }
   RAttrMarker &AttrMarker() { return fMarkerAttr; }
};

} // namespace Experimental
} // namespace ROOT

#endif
