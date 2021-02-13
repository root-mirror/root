// @(#)root/gl:$Id$
// Author:  Matevz Tadel, Feb 2007

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TString.h"

#include "TGLRnrCtx.h"
#include "TGLSceneInfo.h"
#include "TGLSelectBuffer.h"
#include "TGLIncludes.h"
#include "TGLUtil.h"
#include "TGLCamera.h"
#include "TGLFontManager.h"
#include "TGLContext.h"

#include "TError.h"
#include "TMathBase.h"
#include "TMath.h"

#include <list>
#include <algorithm>
#include <cassert>

/** \class TGLRnrCtx
\ingroup opengl
The TGLRnrCtx class aggregates data for a given redering context as
needed by various parts of the ROOT's OpenGL infrastructure. It
serves as a connecting point between the steering part of the
infrastructure (viewer, scene) and concrete rendering classes
(logical, physical shape). It is just a data-holder, there is no
functionality in it.

Development notes:

One RnrCtx is created by each viewer and it is just an extension of
the viewer context that changes along the render
descend. Separating this also has some abstract benefit of hiding
viewer implementation from those that do not need to know it.

Current scene / scene-info part is always initialized by viewer,
scenes can assume they're ok.
*/

ClassImp(TGLRnrCtx);

////////////////////////////////////////////////////////////////////////////////

TGLRnrCtx::TGLRnrCtx(TGLViewerBase* viewer) :
   fViewer    (viewer),
   fCamera    (0),
   fSceneInfo (0),

   fViewerLOD    (kLODUndef),
   fSceneLOD     (kLODUndef),
   fCombiLOD     (kLODUndef),
   fShapeLOD     (kLODUndef),
   fShapePixSize (0),

   fViewerStyle  (kStyleUndef),
   fSceneStyle   (kStyleUndef),

   fViewerWFLineW (0),
   fSceneWFLineW  (0),
   fViewerOLLineW (0),
   fSceneOLLineW  (0),

   fViewerClip   (0),
   fSceneClip    (0),
   fClip         (0),
   fDrawPass     (kPassUndef),

   fStopwatch    (),
   fRenderTimeOut(0.0),
   fIsRunning    (kFALSE),
   fHasTimedOut  (kFALSE),

   fHighlight    (kFALSE),  fHighlightOutline (kFALSE),
   fSelection    (kFALSE),  fSecSelection     (kFALSE),
   fSelectTransparents (kIfNoOpaques),
   fPickRadius   (0),
   fPickRectangle(0),
   fSelectBuffer (0),

   fColorSetStack(0),
   fRenderScale  (1),

   fEventKeySym  (0),

   fDLCaptureOpen (kFALSE),
   fGLCtxIdentity (0),
   fQuadric       (0),

   fGrabImage     (kFALSE),
   fGrabBuffer    (-1),
   fGrabbedImage  (0)
{
   // Constructor.

   fColorSetStack = new lpTGLColorSet_t;
   fColorSetStack->push_back(0);

   fSelectBuffer = new TGLSelectBuffer;
   if (fViewer == 0)
   {
      // Assume external usage, initialize for highest quality.
      fViewerLOD = fSceneLOD = fCombiLOD = fShapeLOD = kLODHigh;
      fViewerStyle = fSceneStyle = kFill;
      fDrawPass = kPassFill;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TGLRnrCtx::~TGLRnrCtx()
{
   gluDeleteQuadric(fQuadric);
   delete fPickRectangle;
   delete fSelectBuffer;
   delete fColorSetStack;
}

////////////////////////////////////////////////////////////////////////////////
/// Return current scene (based on scene-info data).

TGLSceneBase * TGLRnrCtx::GetScene()
{
   return  fSceneInfo->GetScene();
}

////////////////////////////////////////////////////////////////////////////////
/// Return current scene (based on scene-info data).

TGLSceneBase & TGLRnrCtx::RefScene()
{
   return *fSceneInfo->GetScene();
}

/**************************************************************************/

////////////////////////////////////////////////////////////////////////////////
/// Returns true if current render-pass uses filled polygon style.

Bool_t TGLRnrCtx::IsDrawPassFilled() const
{
   return fDrawPass == kPassFill || fDrawPass == kPassOutlineFill;
}


/******************************************************************************/
// Stopwatch
/******************************************************************************/

////////////////////////////////////////////////////////////////////////////////
/// Start the stopwatch.

void TGLRnrCtx:: StartStopwatch()
{
   if (fIsRunning)
      return;

   fStopwatch.Start();
   fIsRunning   = kTRUE;
   fHasTimedOut = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Stop the stopwatch.

void TGLRnrCtx:: StopStopwatch()
{
   fHasTimedOut = fStopwatch.End() > fRenderTimeOut;
   fIsRunning = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Check if the stopwatch went beyond the render time limit.

Bool_t TGLRnrCtx::HasStopwatchTimedOut()
{
   if (fHasTimedOut) return kTRUE;
   if (fIsRunning && fStopwatch.Lap() > fRenderTimeOut)
      fHasTimedOut = kTRUE;
   return fHasTimedOut;
}


/******************************************************************************/
// Selection & picking
/******************************************************************************/

////////////////////////////////////////////////////////////////////////////////
/// Setup context for running selection.
/// x and y are in window coordinates.

void TGLRnrCtx::BeginSelection(Int_t x, Int_t y, Int_t r)
{
   fSelection    = kTRUE;
   fSecSelection = kFALSE;
   fPickRadius   = r;
   if (!fPickRectangle) fPickRectangle = new TGLRect;
   fPickRectangle->Set(x, y, r, r);

   glSelectBuffer(fSelectBuffer->GetBufSize(), fSelectBuffer->GetBuf());
}

////////////////////////////////////////////////////////////////////////////////
/// End selection.

void TGLRnrCtx::EndSelection(Int_t glResult)
{
   fSelection    = kFALSE;
   fSecSelection = kFALSE;
   fPickRadius   = 0;
   delete fPickRectangle; fPickRectangle = 0;

   if (glResult < 0)
   {
      if (fSelectBuffer->CanGrow())
      {
         Warning("TGLRnrCtx::EndSelection",
                 "Select buffer size (%d) insufficient, doubling it.",
                 fSelectBuffer->GetBufSize());
         fSelectBuffer->Grow();
      }
      else
      {
         Warning("TGLRnrCtx::EndSelection",
                 "Select buffer size (%d) insufficient. This is maximum.",
                 fSelectBuffer->GetBufSize());
      }
   }
   fSelectBuffer->ProcessResult(glResult);
}

////////////////////////////////////////////////////////////////////////////////
/// Return current pick rectangle. This is *zero* when
/// selection is not set.

TGLRect * TGLRnrCtx::GetPickRectangle()
{
   return fPickRectangle;
}

////////////////////////////////////////////////////////////////////////////////
/// Return pick radius. If selection is not active it returns 0.

Int_t TGLRnrCtx::GetPickRadius()
{
   return fPickRadius;
}


/**************************************************************************/
// ColorSet access & management
/******************************************************************************/

////////////////////////////////////////////////////////////////////////////////
/// Create copy of current color-set on the top of the stack.

void TGLRnrCtx::PushColorSet()
{
   fColorSetStack->push_back(new TGLColorSet(*fColorSetStack->back()));
}

////////////////////////////////////////////////////////////////////////////////
/// Return reference to current color-set (top of the stack).

TGLColorSet& TGLRnrCtx::ColorSet()
{
   return * fColorSetStack->back();
}

////////////////////////////////////////////////////////////////////////////////
/// Pops the top-most color-set.
/// If only one entry is available, error is printed and the entry remains.

void TGLRnrCtx::PopColorSet()
{
   if (fColorSetStack->size() >= 2)
   {
      delete fColorSetStack->back();
      fColorSetStack->pop_back();
   }
   else
   {
      Error("PopColorSet()", "Attempting to remove the last entry.");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Change the default/bottom color-set.
/// Returns the previous color-set.

TGLColorSet* TGLRnrCtx::ChangeBaseColorSet(TGLColorSet* set)
{
   TGLColorSet* old = fColorSetStack->front();
   fColorSetStack->front() = set;
   return old;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the current base color-set.

TGLColorSet* TGLRnrCtx::GetBaseColorSet()
{
   return fColorSetStack->front();
}

////////////////////////////////////////////////////////////////////////////////
/// Set col if it is different from background, otherwise use
/// current foreground color.

void TGLRnrCtx::ColorOrForeground(Color_t col)
{
   if (fColorSetStack->back()->Background().GetColorIndex() == col)
      TGLUtil::Color(fColorSetStack->back()->Foreground());
   else
      TGLUtil::Color(col);
}

/**************************************************************************/
// Display-list state
/******************************************************************************/

////////////////////////////////////////////////////////////////////////////////
/// Start display-list capture.

void TGLRnrCtx::OpenDLCapture()
{
   assert(fDLCaptureOpen == kFALSE);
   fDLCaptureOpen = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// End display list capture.

void TGLRnrCtx::CloseDLCapture()
{
   assert(fDLCaptureOpen == kTRUE);
   fDLCaptureOpen = kFALSE;
}

/******************************************************************************/
// TGLFont interface
/******************************************************************************/
////////////////////////////////////////////////////////////////////////////////
/// Release font in the GL rendering context.

void TGLRnrCtx::ReleaseFont(TGLFont& font)
{
   fGLCtxIdentity->GetFontManager()->ReleaseFont(font);
}

////////////////////////////////////////////////////////////////////////////////
/// Get font in the GL rendering context.

void TGLRnrCtx::RegisterFontNoScale(Int_t size, Int_t file, Int_t mode, TGLFont& out)
{
   fGLCtxIdentity->GetFontManager()->RegisterFont( size, file, (TGLFont::EMode)mode, out);
}

////////////////////////////////////////////////////////////////////////////////
/// Get font in the GL rendering context.

void TGLRnrCtx::RegisterFontNoScale(Int_t size, const char* name, Int_t mode, TGLFont& out)
{
   fGLCtxIdentity->GetFontManager()->RegisterFont(size, name, (TGLFont::EMode)mode, out);
}

////////////////////////////////////////////////////////////////////////////////
/// Get font in the GL rendering context.
/// The font is scaled relative to current render scale.

void TGLRnrCtx::RegisterFont(Int_t size, Int_t file, Int_t mode, TGLFont& out)
{
  RegisterFontNoScale(TMath::Nint(size*fRenderScale), file, mode, out);
}

////////////////////////////////////////////////////////////////////////////////
/// Get font in the GL rendering context.
/// The font is scaled relative to current render scale.

void TGLRnrCtx::RegisterFont(Int_t size, const char* name, Int_t mode, TGLFont& out)
{
  RegisterFontNoScale(TMath::Nint(size*fRenderScale), name, mode, out);
}

/******************************************************************************/
// fQuadric's initialization.
/******************************************************************************/

////////////////////////////////////////////////////////////////////////////////
/// Initialize fQuadric.

GLUquadric *TGLRnrCtx::GetGluQuadric()
{
   if (!fQuadric) {
      if ((fQuadric = gluNewQuadric())) {
         gluQuadricOrientation(fQuadric, (GLenum)GLU_OUTSIDE);
         gluQuadricNormals(fQuadric, (GLenum)GLU_SMOOTH);
      } else
         Error("TGLRnrCtx::GetGluQuadric", "gluNewQuadric failed");
   }

   return fQuadric;
}


/******************************************************************************/
// Matrix manipulation helpers
/******************************************************************************/

void TGLRnrCtx::ProjectionMatrixPushIdentity()
{
   glMatrixMode(GL_PROJECTION);
   glPushMatrix();
   glLoadIdentity();
   if (Selection())
   {
      TGLRect rect(*GetPickRectangle());
      GetCamera()->WindowToViewport(rect);
      gluPickMatrix(rect.X(), rect.Y(), rect.Width(), rect.Height(),
                    (Int_t*) GetCamera()->RefViewport().CArr());
   }
   glMatrixMode(GL_MODELVIEW);
}

void TGLRnrCtx::ProjectionMatrixPop()
{
   glMatrixMode(GL_PROJECTION);
   glPopMatrix();
   glMatrixMode(GL_MODELVIEW);
}


/**************************************************************************/
// Static helpers
/**************************************************************************/

////////////////////////////////////////////////////////////////////////////////
/// Return string describing the style.

const char* TGLRnrCtx::StyleName(Short_t style)
{
   switch (style)
   {
      case TGLRnrCtx::kFill:       return "Filled Polys";
      case TGLRnrCtx::kWireFrame:  return "Wireframe";
      case TGLRnrCtx::kOutline:    return "Outline";
      default:                     return "Oogaa-dooga style";
   }
}
