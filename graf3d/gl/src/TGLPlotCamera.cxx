// @(#)root/gl:$Id$
// Author: Timur Pocheptsov

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/
#include "TGLPlotCamera.h"
#include "TGLIncludes.h"
#include "TGLUtil.h"

/** \class TGLPlotCamera
\ingroup opengl
Camera for TGLPlotPainter and sub-classes.
*/

ClassImp(TGLPlotCamera);

////////////////////////////////////////////////////////////////////////////////
///Construct camera for plot painters.

TGLPlotCamera::TGLPlotCamera() :
   fZoom(1.), fShift(1.5), fCenter(),
   fVpChanged(kFALSE)
{
   fOrthoBox[0] = 1.;
   fOrthoBox[1] = 1.;
   fOrthoBox[2] = -100.;
   fOrthoBox[3] = 100.;
}

////////////////////////////////////////////////////////////////////////////////
///Setup viewport, if it was changed, plus reset arcball.

void TGLPlotCamera::SetViewport(const TGLRect &vp)
{
   if (vp.Width() != fViewport.Width() || vp.Height() != fViewport.Height() ||
       vp.X() != fViewport.X() || vp.Y() != fViewport.Y())
   {
      fVpChanged = kTRUE;
      fArcBall.SetBounds(vp.Width(), vp.Height());
      fViewport = vp;

   } else
      fVpChanged = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
///'box' is the TGLPlotPainter's back box's coordinates.

void TGLPlotCamera::SetViewVolume(const TGLVertex3* /* box */)
{
/*   fCenter[0] = (box[0].X() + box[1].X()) / 2;
   fCenter[1] = (box[0].Y() + box[2].Y()) / 2;
   fCenter[2] = (box[0].Z() + box[4].Z()) / 2;
   const Double_t maxDim = box[1].X() - box[0].X();
   fOrthoBox[0] = maxDim;
   fOrthoBox[1] = maxDim;
   fOrthoBox[2] = -100 * maxDim;//100?
   fOrthoBox[3] = 100 * maxDim;
   fShift = maxDim * 1.5;*/
}

////////////////////////////////////////////////////////////////////////////////
///User clicks somewhere (px, py).

void TGLPlotCamera::StartRotation(Int_t px, Int_t py)
{
   fArcBall.Click(TPoint(px, py));
}

////////////////////////////////////////////////////////////////////////////////
///Mouse movement.

void TGLPlotCamera::RotateCamera(Int_t px, Int_t py)
{
   fArcBall.Drag(TPoint(px, py));
}

////////////////////////////////////////////////////////////////////////////////
///User clicks somewhere (px, py).

void TGLPlotCamera::StartPan(Int_t px, Int_t py)
{
   fMousePos.fX = px;
   fMousePos.fY = fViewport.Height() - py;
}

////////////////////////////////////////////////////////////////////////////////
///Pan camera.

void TGLPlotCamera::Pan(Int_t px, Int_t py)
{
   py = fViewport.Height() - py;

   //Extract gl matrices.
   Double_t mv[16] = {0.};
   glGetDoublev(GL_MODELVIEW_MATRIX, mv);
   Double_t pr[16] = {0.};
   glGetDoublev(GL_PROJECTION_MATRIX, pr);
   Int_t vp[] = {0, 0, fViewport.Width(), fViewport.Height()};
   //Adjust pan vector.
   TGLVertex3 start, end;
   gluUnProject(fMousePos.fX, fMousePos.fY, 1., mv, pr, vp, &start.X(), &start.Y(), &start.Z());
   gluUnProject(px, py, 1., mv, pr, vp, &end.X(), &end.Y(), &end.Z());
   fTruck += (start - end) /= 2.;
   //
   fMousePos.fX = px;
   fMousePos.fY = py;
}

////////////////////////////////////////////////////////////////////////////////
///Viewport and projection.

void TGLPlotCamera::SetCamera()const
{
   glViewport(fViewport.X(), fViewport.Y(), fViewport.Width(), fViewport.Height());

   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   glOrtho(
           -fOrthoBox[0] * fZoom,
            fOrthoBox[0] * fZoom,
           -fOrthoBox[1] * fZoom,
            fOrthoBox[1] * fZoom,
            fOrthoBox[2],
            fOrthoBox[3]
          );

   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();
}

////////////////////////////////////////////////////////////////////////////////
///Applies rotations and translations before drawing

void TGLPlotCamera::Apply(Double_t phi, Double_t theta)const
{
   glTranslated(0., 0., -fShift);
   glMultMatrixd(fArcBall.GetRotMatrix());
   glRotated(theta - 90., 1., 0., 0.);
   glRotated(phi, 0., 0., 1.);
   glTranslated(-fTruck[0], -fTruck[1], -fTruck[2]);
//   glTranslated(-fCenter[0], -fCenter[1], -fCenter[2]);
}

////////////////////////////////////////////////////////////////////////////////
///viewport[0]

Int_t TGLPlotCamera::GetX()const
{
   return fViewport.X();
}

////////////////////////////////////////////////////////////////////////////////
///viewport[1]

Int_t TGLPlotCamera::GetY()const
{
   return fViewport.Y();
}

////////////////////////////////////////////////////////////////////////////////
///viewport[2]

Int_t TGLPlotCamera::GetWidth()const
{
   return Int_t(fViewport.Width());
}

////////////////////////////////////////////////////////////////////////////////
///viewport[3]

Int_t TGLPlotCamera::GetHeight()const
{
   return Int_t(fViewport.Height());
}

////////////////////////////////////////////////////////////////////////////////
///Zoom in.

void TGLPlotCamera::ZoomIn()
{
   fZoom /= 1.2;
}

////////////////////////////////////////////////////////////////////////////////
///Zoom out.

void TGLPlotCamera::ZoomOut()
{
   fZoom *= 1.2;
}
