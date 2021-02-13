// @(#)root/gl:$Id$
// Author:  Matevz Tadel, Jun 2007

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TF2GL.h"

#include <TF2.h>
#include <TF3.h>
#include <TH2.h>

#include "TGLSurfacePainter.h"
#include "TGLTF3Painter.h"
#include "TGLAxisPainter.h"

#include "TGLRnrCtx.h"

#include "TGLIncludes.h"

/** \class TF2GL
\ingroup opengl
GL renderer for TF2.
TGLPlotPainter is used internally.
*/

ClassImp(TF2GL);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TF2GL::TF2GL() : TGLPlot3D(), fM(0), fH(0)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TF2GL::~TF2GL()
{
   delete fH;
}

////////////////////////////////////////////////////////////////////////////////
/// Set model object.

Bool_t TF2GL::SetModel(TObject* obj, const Option_t* opt)
{
   TString option(opt);
   option.ToLower();

   fM = SetModelDynCast<TF2>(obj);

   fH = (TH2*) fM->CreateHistogram();
   if (!fH) return kFALSE;

   fH->GetZaxis()->SetLimits(fH->GetMinimum(), fH->GetMaximum());

   if (dynamic_cast<TF3*>(fM))
      SetPainter( new TGLTF3Painter((TF3*)fM, fH, 0, &fCoord) );
   else
      SetPainter( new TGLSurfacePainter(fH, 0, &fCoord) );

   if (option.Index("sph") != kNPOS)
      fCoord.SetCoordType(kGLSpherical);
   else if (option.Index("pol") != kNPOS)
      fCoord.SetCoordType(kGLPolar);
   else if (option.Index("cyl") != kNPOS)
      fCoord.SetCoordType(kGLCylindrical);

   fPlotPainter->AddOption(option);
   fPlotPainter->InitGeometry();

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Setup bounding-box.

void TF2GL::SetBBox()
{
   fBoundingBox.Set(fPlotPainter->RefBackBox().Get3DBox());
}

////////////////////////////////////////////////////////////////////////////////
/// Render the object.

void TF2GL::DirectDraw(TGLRnrCtx & rnrCtx) const
{
   fPlotPainter->RefBackBox().FindFrontPoint();

   glPushAttrib(GL_ENABLE_BIT | GL_LIGHTING_BIT);

   glEnable(GL_NORMALIZE);
   glDisable(GL_COLOR_MATERIAL);

   fPlotPainter->InitGL();
   fPlotPainter->DrawPlot();

   glDisable(GL_CULL_FACE);
   glPopAttrib();

   // Axes
   TGLAxisPainterBox axe_painter;
   axe_painter.SetUseAxisColors(kFALSE);
   axe_painter.SetFontMode(TGLFont::kPixmap);
   axe_painter.PlotStandard(rnrCtx, fH, fBoundingBox);
}
