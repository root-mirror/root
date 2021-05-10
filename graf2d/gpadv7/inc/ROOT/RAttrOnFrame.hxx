/*************************************************************************
 * Copyright (C) 1995-2021, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RAttrOnFrame
#define ROOT7_RAttrOnFrame

#include <ROOT/RAttrValue.hxx>

namespace ROOT {
namespace Experimental {

/** \class RAttrOnFrame
\ingroup GpadROOT7
\author Sergey Linev <s.linev@gsi.de>
\date 2021-05-06
\brief Class which add onframe property for drawable, used when drawable can be drawn on frame or outside frame
\warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!
*/

class RAttrOnFrame {

   RAttrValue<bool> fOnFrame;    ///<! is drawable drawn on the frame or not
   RAttrValue<bool> fClipping;   ///<! if drawable should be clipped by frame

public:

   RAttrOnFrame(RDrawable *drawable) : fOnFrame(drawable, "onframe", false), fClipping(drawable, "clipping", false) {}

   void SetOnFrame(bool on = true) { fOnFrame = on; }
   bool GetOnFrame() const { return fOnFrame; }

   void SetClipping(bool on = true) { fClipping = on; }
   bool GetClipping() const { return fClipping; }
};

} // namespace Experimental
} // namespace ROOT

#endif // ROOT7_RAttrOnFrame
