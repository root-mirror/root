// @(#)root/gl:$Id$
// Author:  Richard Maunder  16/09/2005

/*************************************************************************
 * Copyright (C) 1995-2005, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TGLManip.h"
#include "TGLUtil.h"
#include "TGLCamera.h"
#include "TGLPhysicalShape.h"
#include "TGLIncludes.h"

/** \class TGLManip
\ingroup opengl
Abstract base class for viewer manipulators, which allow direct in
viewer manipulation of a (TGlPhysicalShape) object - currently
translation, scaling and rotation along/round objects local axes.
See derived classes for these implementations.

This class provides binding to the zero or one manipulated physical,
hit testing (selection) for manipulator sub component (widget), and
some common mouse action handling/tracking.
*/

ClassImp(TGLManip);

////////////////////////////////////////////////////////////////////////////////
/// Construct a manipulator object, bound to supplied viewer, and no
/// physical shape.

TGLManip::TGLManip() :
   fShape(0),
   fSelectedWidget(0), fActive(kFALSE),
   fFirstMouse(0, 0),
   fLastMouse(0, 0)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Construct a manipulator object, bound to supplied physical shape.

TGLManip::TGLManip(TGLPhysicalShape* shape) :
   fShape(shape),
   fSelectedWidget(0), fActive(kFALSE),
   fFirstMouse(0, 0),
   fLastMouse(0, 0)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.

TGLManip::TGLManip(const TGLManip& gm) :
  TVirtualGLManip(gm),
  fShape(gm.fShape),
  fSelectedWidget(gm.fSelectedWidget),
  fActive(gm.fActive),
  fFirstMouse(gm.fFirstMouse),
  fLastMouse(gm.fLastMouse)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Assignment operator.

TGLManip& TGLManip::operator=(const TGLManip& gm)
{
   if(this!=&gm) {
      TVirtualGLManip::operator=(gm);
      fShape=gm.fShape;
      fSelectedWidget=gm.fSelectedWidget;
      fActive=gm.fActive;
      fFirstMouse=gm.fFirstMouse;
      fLastMouse=gm.fLastMouse;
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Destroy manipulator object.

TGLManip::~TGLManip()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Returns color to be used for given widget.

const UChar_t* TGLManip::ColorFor(UInt_t widget) const
{
   if (widget == fSelectedWidget)
   {
      return TGLUtil::fgYellow;
   }
   else
   {
      switch (widget)
      {
         case 1:  return TGLUtil::fgRed;
         case 2:  return TGLUtil::fgGreen;
         case 3:  return TGLUtil::fgBlue;
         default: return TGLUtil::fgGrey;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Handle a mouse button event - return kTRUE if processed, kFALSE otherwise

Bool_t TGLManip::HandleButton(const Event_t& event, const TGLCamera& /*camera*/)
{
   // Only interested in Left mouse button actions
   if (event.fCode != kButton1) {
      return kFALSE;
   }

   // Mouse down on selected widget?
   if (event.fType == kButtonPress && fSelectedWidget != 0) {
      fFirstMouse.SetX(event.fX);
      fFirstMouse.SetY(event.fY);
      fLastMouse.SetX(event.fX);
      fLastMouse.SetY(event.fY);
      fActive = kTRUE;
      return kTRUE;
   } else if (event.fType == kButtonRelease && fActive) {
      fActive = kFALSE;
      return kTRUE;
   } else {
      return kFALSE;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Handle a mouse button event - return kTRUE if widget selection change
/// kFALSE otherwise

Bool_t TGLManip::HandleMotion(const Event_t&   /*event*/,
                              const TGLCamera& /*camera*/)
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Calculates base and axis scale factor (in world units) for
/// drawing manipulators with reasonable size range in current
/// camera.

void TGLManip::CalcDrawScale(const TGLBoundingBox& box,
                             const TGLCamera&      camera,
                             Double_t&             base,
                             TGLVector3            axis[3]) const
{
   // Calculate a base scale
   base = box.Extents().Mag() / 100.0;

   // Clamp this base scale to a viewport pixel range
   // Allow some variation so zooming is noticeable
   TGLVector3 pixelInWorld = camera.ViewportDeltaToWorld(box.Center(), 1, 1);
   Double_t pixelScale = pixelInWorld.Mag();
   if (base < pixelScale * 3.0) {
      base = pixelScale * 3.0;
   } else if (base > pixelScale * 6.0) {
      base = pixelScale * 6.0;
   }

   // Calculate some axis scales
   for (UInt_t i = 0; i<3; i++) {
      if (box.IsEmpty()) {
         axis[i] = box.Axis(i, kTRUE)*base*-10.0;
      } else {
         axis[i] = box.Axis(i, kFALSE)*-0.51;
         if (axis[i].Mag() < base*10.0) {
            axis[i] = box.Axis(i, kTRUE)*base*-10.0;
         }
      }
   }
}
