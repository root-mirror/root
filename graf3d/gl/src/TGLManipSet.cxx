// @(#)root/gl:$Id$
// Author:  Matevz Tadel, Feb 2007

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TGLManipSet.h"

#include "TGLTransManip.h"
#include "TGLScaleManip.h"
#include "TGLRotateManip.h"

#include "TGLPhysicalShape.h"
#include "TGLRnrCtx.h"
#include "TGLSelectRecord.h"

#include "TGLIncludes.h"

#include <KeySymbols.h>

/** \class TGLManipSet
\ingroup opengl

Combine all available manipulators in a collection.

At first I wanted to merge them back into TGLManip (to have a
single class) but then it seemed somehow messy.
Maybe next time.
*/

ClassImp(TGLManipSet);

TGLManipSet::TGLManipSet() :
   TGLOverlayElement(kViewer),
   fType     (kTrans),
   fDrawBBox (kFALSE)
{
   // Constructor.

   fManip[kTrans]  = new TGLTransManip;
   fManip[kScale]  = new TGLScaleManip;
   fManip[kRotate] = new TGLRotateManip;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TGLManipSet::~TGLManipSet()
{
   for (Int_t i=kTrans; i<kEndType; ++i)
      delete fManip[i];
}

////////////////////////////////////////////////////////////////////////////////
/// Set phys-shape, override of virtual from TGLPShapeRef.
/// Forward to all managed manipulators.

void TGLManipSet::SetPShape(TGLPhysicalShape* shape)
{
   TGLPShapeRef::SetPShape(shape);
   for (Int_t i=kTrans; i<kEndType; ++i)
      fManip[i]->Attach(shape);
}

////////////////////////////////////////////////////////////////////////////////
/// Mouse has entered this element.
/// Always accept.

Bool_t TGLManipSet::MouseEnter(TGLOvlSelectRecord& /*selRec*/)
{
   TGLManip* manip = GetCurrentManip();
   manip->SetActive(kFALSE);
   manip->SetSelectedWidget(0);
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle overlay event.
/// Return TRUE if event was handled.

Bool_t TGLManipSet::Handle(TGLRnrCtx&          rnrCtx,
                           TGLOvlSelectRecord& selRec,
                           Event_t*            event)
{
   TGLManip* manip = GetCurrentManip();

   switch (event->fType)
   {
      case kButtonPress:
      {
         return manip->HandleButton(*event, rnrCtx.RefCamera());
      }
      case kButtonRelease:
      {
         manip->SetActive(kFALSE);
         return kTRUE;
      }
      case kMotionNotify:
      {
         if (manip->GetActive())
            return manip->HandleMotion(*event, rnrCtx.RefCamera());
         if (selRec.GetCurrItem() != manip->GetSelectedWidget())
         {
            manip->SetSelectedWidget(selRec.GetCurrItem());
            return kTRUE;
         }
         return kFALSE;
      }
      case kGKeyPress:
      {
         switch (rnrCtx.GetEventKeySym())
         {
            case kKey_V: case kKey_v:
               SetManipType(kTrans);
               return kTRUE;
            case kKey_C: case kKey_c:
               SetManipType(kRotate);
               return kTRUE;
            case kKey_X: case kKey_x:
               SetManipType(kScale);
               return kTRUE;
            default:
               return kFALSE;
         }
      }
      default:
      {
         return kFALSE;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Mouse has left the element.

void TGLManipSet::MouseLeave()
{
   TGLManip* manip = GetCurrentManip();
   manip->SetActive(kFALSE);
   manip->SetSelectedWidget(0);
}

////////////////////////////////////////////////////////////////////////////////
/// Render the manipulator and bounding-box.

void TGLManipSet::Render(TGLRnrCtx& rnrCtx)
{
   if (fPShape == 0)
      return;

   if (rnrCtx.Selection())
   {
      TGLUtil::SetDrawQuality(12);
      fManip[fType]->Draw(rnrCtx.RefCamera());
      TGLUtil::ResetDrawQuality();
   } else {
      fManip[fType]->Draw(rnrCtx.RefCamera());
   }

   if (fDrawBBox && ! rnrCtx.Selection())
   {
      // TODO: This must be replaced by some color in rnrCtx,
      // like def-overlay-color, background-color, foreground-color
      // Or at least bkgcol ... i can then find high contrast.
      TGLUtil::Color(rnrCtx.ColorSet().Markup());
      glDisable(GL_LIGHTING);
      fPShape->BoundingBox().Draw();
      glEnable(GL_LIGHTING);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set manipulator type, range checked.

void TGLManipSet::SetManipType(Int_t type)
{
   if (type < 0 || type >= kEndType)
      return;
   fType = (EManip) type;
}
