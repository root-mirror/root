// @(#)root/gl:$Id$
// Author:  Matevz Tadel, Jun 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TGLParametricEquationGL.h"

#include "TGLParametric.h"

#include "TGLSurfacePainter.h"
#include "TGLTF3Painter.h"

#include "TGLRnrCtx.h"
#include "TGLIncludes.h"

/** \class TGLParametricEquationGL
\ingroup opengl
GL-renderer wrapper for TGLParametricEquation.
This allows rendering of parametric-equations in standard GL viewer.
*/

ClassImp(TGLParametricEquationGL);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TGLParametricEquationGL::TGLParametricEquationGL() : TGLPlot3D(), fM(0)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TGLParametricEquationGL::~TGLParametricEquationGL()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Set model object.

Bool_t TGLParametricEquationGL::SetModel(TObject* obj, const Option_t* opt)
{
   fM = SetModelDynCast<TGLParametricEquation>(obj);

   SetPainter( new TGLParametricPlot(fM, 0) );
   TString option(opt);
   fPlotPainter->AddOption(option);
   fPlotPainter->InitGeometry();

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Setup bounding-box.

void TGLParametricEquationGL::SetBBox()
{
   fBoundingBox.Set(fPlotPainter->RefBackBox().Get3DBox());
}

////////////////////////////////////////////////////////////////////////////////
/// Render the object.

void TGLParametricEquationGL::DirectDraw(TGLRnrCtx& /*rnrCtx*/) const
{
   fPlotPainter->RefBackBox().FindFrontPoint();
   glPushAttrib(GL_ENABLE_BIT | GL_LIGHTING_BIT);
   glEnable(GL_NORMALIZE);
   fPlotPainter->InitGL();
   fPlotPainter->DrawPlot();
   glPopAttrib();
}
