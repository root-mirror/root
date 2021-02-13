// @(#)root/gl:$Id$
// Author:  Matevz Tadel, Jun 2007

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TH2GL.h"
#include "TH2.h"
#include "TAxis.h"

#include "TGLSurfacePainter.h"
#include "TGLHistPainter.h"
#include "TGLLegoPainter.h"
#include "TGLBoxPainter.h"
#include "TGLTF3Painter.h"
#include "TGLAxisPainter.h"
#include "TGLCamera.h"

#include "TGLRnrCtx.h"

#include "TGLIncludes.h"

/** \class TH2GL
\ingroup opengl
Rendering of TH2 and derived classes.
Interface to plot-painters also used for gl-in-pad.
*/

ClassImp(TH2GL);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TH2GL::TH2GL() :
   TGLPlot3D(), fM(0)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TH2GL::~TH2GL()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Set model object.

Bool_t TH2GL::SetModel(TObject* obj, const Option_t* opt)
{
   TString option(opt);
   option.ToLower();

   fM = SetModelDynCast<TH2>(obj);

   // Plot type
   if (option.Index("surf") != kNPOS)
      SetPainter( new TGLSurfacePainter(fM, 0, &fCoord) );
   else
      SetPainter( new TGLLegoPainter(fM, 0, &fCoord) );

   if (option.Index("sph") != kNPOS)
      fCoord.SetCoordType(kGLSpherical);
   else if (option.Index("pol") != kNPOS)
      fCoord.SetCoordType(kGLPolar);
   else if (option.Index("cyl") != kNPOS)
      fCoord.SetCoordType(kGLCylindrical);

   fPlotPainter->AddOption(option);

   Ssiz_t pos = option.Index("fb");
   if (pos != kNPOS) {
      option.Remove(pos, 2);
      fPlotPainter->SetDrawFrontBox(kFALSE);
   }

   pos = option.Index("bb");
   if (pos != kNPOS)
      fPlotPainter->SetDrawBackBox(kFALSE);

   pos = option.Index("a");
   if (pos != kNPOS)
      fPlotPainter->SetDrawAxes(kFALSE);

   fPlotPainter->InitGeometry();

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Setup bounding-box.

void TH2GL::SetBBox()
{
   fBoundingBox.Set(fPlotPainter->RefBackBox().Get3DBox());
}

////////////////////////////////////////////////////////////////////////////////
/// Render the object.

void TH2GL::DirectDraw(TGLRnrCtx & rnrCtx) const
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
   if (fPlotPainter->GetDrawAxes()) {
      TGLAxisPainterBox axe_painter;
      axe_painter.SetUseAxisColors(kFALSE);
      axe_painter.SetFontMode(TGLFont::kPixmap);
      axe_painter.PlotStandard(rnrCtx, fM, fBoundingBox);
   }
}
