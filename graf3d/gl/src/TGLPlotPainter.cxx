// @(#)root/gl:$Id$
// Author:  Timur Pocheptsov  14/06/2006

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/
#include <cstdio>

#include "TVirtualPad.h"
#include "TVirtualPS.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TError.h"
#include "TColor.h"
#include "TAxis.h"
#include "TMath.h"
#include "TList.h"
#include "TH2Poly.h"
#include "TH1.h"
#include "TH3.h"
#include "TF3.h"
#include "TROOT.h"
#include "TVirtualMutex.h"

#include "TGLPlotPainter.h"
#include "TGLPlotCamera.h"
#include "TGLIncludes.h"
#include "TGLAdapter.h"
#include "TGLOutput.h"
#include "TGLUtil.h"
#include "TGL5D.h"
#include "gl2ps.h"

/** \class TGLPlotPainter
\ingroup opengl
Base class for plot-painters that provide GL rendering of various
2D and 3D histograms, functions and parametric surfaces.
*/

ClassImp(TGLPlotPainter);

////////////////////////////////////////////////////////////////////////////////
///TGLPlotPainter's ctor.

TGLPlotPainter::TGLPlotPainter(TH1 *hist, TGLPlotCamera *camera, TGLPlotCoordinates *coord,
                               Bool_t xoy, Bool_t xoz, Bool_t yoz)
                  : fPadColor(0),
                    fPhysicalShapeColor(0),
                    fPadPhi(45.),
                    fPadTheta(0.),
                    fHist(hist),
                    fXAxis(hist->GetXaxis()),
                    fYAxis(hist->GetYaxis()),
                    fZAxis(hist->GetZaxis()),
                    fCoord(coord),
                    fCamera(camera),
                    fUpdateSelection(kTRUE),
                    fSelectionPass(kFALSE),
                    fSelectedPart(0),
                    fXOZSectionPos(0.),
                    fYOZSectionPos(0.),
                    fXOYSectionPos(0.),
                    fBackBox(xoy, xoz, yoz),
                    fBoxCut(&fBackBox),
                    fHighColor(kFALSE),
                    fSelectionBase(kTrueColorSelectionBase),
                    fDrawPalette(kFALSE),
                    fDrawAxes(kTRUE)
{
   if (gPad) {
      fPadPhi   = gPad->GetPhi();
      fPadTheta = gPad->GetTheta();
   }
}

////////////////////////////////////////////////////////////////////////////////
///TGLPlotPainter's ctor.

TGLPlotPainter::TGLPlotPainter(TGL5DDataSet *data, TGLPlotCamera *camera, TGLPlotCoordinates *coord)
                  : fPadColor(0),
                    fPhysicalShapeColor(0),
                    fPadPhi(45.),
                    fPadTheta(0.),
                    fHist(0),
                    fXAxis(data->GetXAxis()),
                    fYAxis(data->GetYAxis()),
                    fZAxis(data->GetZAxis()),
                    fCoord(coord),
                    fCamera(camera),
                    fUpdateSelection(kTRUE),
                    fSelectionPass(kFALSE),
                    fSelectedPart(0),
                    fXOZSectionPos(0.),
                    fYOZSectionPos(0.),
                    fXOYSectionPos(0.),
                    fBackBox(kFALSE, kFALSE, kFALSE),
                    fBoxCut(&fBackBox),
                    fHighColor(kFALSE),
                    fSelectionBase(kTrueColorSelectionBase),
                    fDrawPalette(kFALSE),
                    fDrawAxes(kTRUE)
{
   if (gPad) {
      fPadPhi   = gPad->GetPhi();
      fPadTheta = gPad->GetTheta();
   }
}

////////////////////////////////////////////////////////////////////////////////
///TGLPlotPainter's ctor.

TGLPlotPainter::TGLPlotPainter(TGLPlotCamera *camera)
                  : fPadColor(0),
                    fPhysicalShapeColor(0),
                    fPadPhi(45.),
                    fPadTheta(0.),
                    fHist(0),
                    fXAxis(0),
                    fYAxis(0),
                    fZAxis(0),
                    fCoord(0),
                    fCamera(camera),
                    fUpdateSelection(kTRUE),
                    fSelectionPass(kFALSE),
                    fSelectedPart(0),
                    fXOZSectionPos(0.),
                    fYOZSectionPos(0.),
                    fXOYSectionPos(0.),
                    fBackBox(kFALSE, kFALSE, kFALSE),
                    fBoxCut(&fBackBox),
                    fHighColor(kFALSE),
                    fSelectionBase(kTrueColorSelectionBase),
                    fDrawPalette(kFALSE),
                    fDrawAxes(kTRUE)
{
   if (gPad) {
      fPadPhi   = gPad->GetPhi();
      fPadTheta = gPad->GetTheta();
   }
}

////////////////////////////////////////////////////////////////////////////////
///Draw lego/surf/whatever you can.

void TGLPlotPainter::Paint()
{
   R__LOCKGUARD(gROOTMutex);

   fHighColor = kFALSE;
   fSelectionBase = fHighColor ? kHighColorSelectionBase : kTrueColorSelectionBase;

   int vp[4] = {};
   glGetIntegerv(GL_VIEWPORT, vp);

   //GL pad painter does not use depth test,
   //so, switch it on now.
   glDepthMask(GL_TRUE);//[0
   //
   InitGL();
   //Save material/light properties in a stack.
   glPushAttrib(GL_LIGHTING_BIT);

   //Save projection and modelview matrix, used by glpad.
   SaveProjectionMatrix();
   SaveModelviewMatrix();

   //glOrtho etc.
   fCamera->SetCamera();
   //
   glClear(GL_DEPTH_BUFFER_BIT);
   //Set light.
   const Float_t pos[] = {0.f, 0.f, 0.f, 1.f};
   glLightfv(GL_LIGHT0, GL_POSITION, pos);
   //Set transformation - shift and rotate the scene.
   fCamera->Apply(fPadPhi, fPadTheta);
   fBackBox.FindFrontPoint();

   if (gVirtualPS)
      PrintPlot();

   DrawPlot();
   //Restore material properties from stack.
   glPopAttrib();
   //
   DeInitGL();//Disable/enable, what concrete plot painter enabled/disabled

   //Restore projection and modelview matrices.
   RestoreProjectionMatrix();
   RestoreModelviewMatrix();

   glViewport(vp[0], vp[1], vp[2], vp[3]);
   //GL pad painter does not use depth test, so,
   //switch it off now.
   glDepthMask(GL_FALSE);//0]

   if (fCoord && fCoord->GetCoordType() == kGLCartesian && fDrawAxes) {

      Bool_t old = gPad->TestBit(TGraph::kClipFrame);
      if (!old)
         gPad->SetBit(TGraph::kClipFrame);

      //Viewport on retina is bigger than real pixel coordinates in
      //a pad, scale it back.
      TGLUtil::InitializeIfNeeded();
      Float_t scale = TGLUtil::GetScreenScalingFactor();
      if (scale < 1.f)//Just ignore this.
         scale = 1.f;

      const Int_t viewport[] = {Int_t(fCamera->GetX() / scale),
                                Int_t(fCamera->GetY() / scale),
                                Int_t(fCamera->GetWidth() / scale),
                                Int_t(fCamera->GetHeight() / scale)};
      Rgl::DrawAxes(fBackBox.GetFrontPoint(), viewport, fBackBox.Get2DBox(), fCoord, fXAxis, fYAxis, fZAxis);
      if (fDrawPalette)
         DrawPaletteAxis();

      if (!old)
         gPad->ResetBit(TGraph::kClipFrame);
   } else if(fDrawPalette)
      DrawPaletteAxis();

}

////////////////////////////////////////////////////////////////////////////////
/// Generate PS using gl2ps

void TGLPlotPainter::PrintPlot()const
{
   using namespace std;

   TGLOutput::StartEmbeddedPS();

   FILE *output = fopen(gVirtualPS->GetName(), "a");
   if (!output) {
      Error("TGLPlotPainter::PrintPlot", "Could not (re)open ps file for GL output");
      //As soon as we started embedded ps, we have to close it before exiting.
      TGLOutput::CloseEmbeddedPS();
      return;
   }

   Int_t gl2psFormat = GL2PS_EPS;
   Int_t gl2psSort   = GL2PS_BSP_SORT;
   Int_t buffsize    = 0;
   Int_t state       = GL2PS_OVERFLOW;
   GLint gl2psoption = GL2PS_USE_CURRENT_VIEWPORT |
                       GL2PS_SILENT               |
                       GL2PS_BEST_ROOT            |
                       GL2PS_OCCLUSION_CULL       |
                       0;

   while (state == GL2PS_OVERFLOW) {
      buffsize += 1024*1024;
      gl2psBeginPage ("ROOT Scene Graph", "ROOT", NULL,
                      gl2psFormat, gl2psSort, gl2psoption,
                      GL_RGBA, 0, NULL,0, 0, 0,
                      buffsize, output, NULL);
      DrawPlot();
      state = gl2psEndPage();
   }

   fclose(output);
   TGLOutput::CloseEmbeddedPS();
   glFlush();
}

////////////////////////////////////////////////////////////////////////////////
///Read color buffer content to find selected object

Bool_t TGLPlotPainter::PlotSelected(Int_t px, Int_t py)
{
   if (fUpdateSelection) {
      //Save projection and modelview matrix, used by glpad.
      glMatrixMode(GL_PROJECTION);//[1
      glPushMatrix();
      glMatrixMode(GL_MODELVIEW);//[2
      glPushMatrix();

      fSelectionPass = kTRUE;
      fCamera->SetCamera();

      glDepthMask(GL_TRUE);
      glClearColor(0.f, 0.f, 0.f, 0.f);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

      fCamera->Apply(fPadPhi, fPadTheta);
      DrawPlot();

      glFinish();
      //fSelection.ReadColorBuffer(fCamera->GetWidth(), fCamera->GetHeight());

      fSelection.ReadColorBuffer(fCamera->GetX(), fCamera->GetY(), fCamera->GetWidth(), fCamera->GetHeight());
      fSelectionPass   = kFALSE;
      fUpdateSelection = kFALSE;

      glDepthMask(GL_FALSE);
      glDisable(GL_DEPTH_TEST);

      //Restore projection and modelview matrices.
      glMatrixMode(GL_PROJECTION);//1]
      glPopMatrix();
      glMatrixMode(GL_MODELVIEW);//2]
      glPopMatrix();
   }

   px -= Int_t(gPad->GetXlowNDC() * gPad->GetWw());
   py -= Int_t(gPad->GetWh() - gPad->YtoAbsPixel(gPad->GetY1()));

   //Pixel coords can be affected by scaling factor on retina displays.
   TGLUtil::InitializeIfNeeded();
   const Float_t scale = TGLUtil::GetScreenScalingFactor();

   if (scale > 1.f) {
      px *= scale;
      py *= scale;
   }

   //py = fCamera->GetHeight() - py;
   //Y is a number of a row, x - column.
   std::swap(px, py);
   Int_t newSelected(Rgl::ColorToObjectID(fSelection.GetPixelColor(px, py), fHighColor));

   if (newSelected != fSelectedPart) {
      //New object was selected (or surface deselected) - re-paint.
      fSelectedPart = newSelected;
      gPad->Update();
   }

   return fSelectedPart ? kTRUE : kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
///Used in a pad.

void TGLPlotPainter::SetPadColor(const TColor *c)
{
   fPadColor = c;
}

////////////////////////////////////////////////////////////////////////////////
///Set plot's back box color.

void TGLPlotPainter::SetFrameColor(const TColor *c)
{
   fBackBox.SetFrameColor(c);
}

////////////////////////////////////////////////////////////////////////////////
///Selection must be updated.

void TGLPlotPainter::InvalidateSelection()
{
   fUpdateSelection = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
///Get pad color.

const TColor *TGLPlotPainter::GetPadColor()const
{
   return fPadColor;
}

////////////////////////////////////////////////////////////////////////////////
///Create dynamic profile using selected plane

void TGLPlotPainter::MoveSection(Int_t px, Int_t py)
{
   //Coordinates are expected to be fixed for retina!

   const TGLVertex3 *frame = fBackBox.Get3DBox();
   const Int_t frontPoint  = fBackBox.GetFrontPoint();

   if (fSelectedPart == 1) {
      fXOYSectionPos = frame[0].Z();
      fSelectedPart = 6;
   } else if (fSelectedPart == 2) {
      if (frontPoint == 2) {
         fXOZSectionPos = frame[0].Y();
         fSelectedPart = 4;
      } else if (!frontPoint) {
         fXOZSectionPos = frame[2].Y();
         fSelectedPart = 4;
      } else if (frontPoint == 1) {
         fYOZSectionPos = frame[0].X();
         fSelectedPart = 5;
      } else if (frontPoint == 3) {
         fYOZSectionPos = frame[1].X();
         fSelectedPart = 5;
      }
   } else if (fSelectedPart == 3) {
      if (frontPoint == 2) {
         fYOZSectionPos = frame[0].X();
         fSelectedPart = 5;
      } else if (!frontPoint) {
         fYOZSectionPos = frame[1].X();
         fSelectedPart = 5;
      } else if (frontPoint == 1) {
         fXOZSectionPos = frame[2].Y();
         fSelectedPart = 4;
      } else if (frontPoint == 3) {
         fXOZSectionPos = frame[0].Y();
         fSelectedPart = 4;
      }
   }

   Double_t mv[16] = {0.};
   glGetDoublev(GL_MODELVIEW_MATRIX, mv);
   Double_t pr[16] = {0.};
   glGetDoublev(GL_PROJECTION_MATRIX, pr);
   Int_t vp[4] = {0};
   glGetIntegerv(GL_VIEWPORT, vp);
   Double_t winVertex[3] = {0.};

   if (fSelectedPart == 6)
      gluProject(0., 0., fXOYSectionPos, mv, pr, vp, &winVertex[0], &winVertex[1], &winVertex[2]);
   else
      gluProject(fSelectedPart == 5 ? fYOZSectionPos : 0.,
                 fSelectedPart == 4 ? fXOZSectionPos : 0.,
                 0., mv, pr, vp,
                 &winVertex[0], &winVertex[1], &winVertex[2]);
   winVertex[0] += px - fMousePosition.fX;
   winVertex[1] += py - fMousePosition.fY;
   Double_t newPoint[3] = {0.};
   gluUnProject(winVertex[0], winVertex[1], winVertex[2], mv, pr, vp,
                newPoint, newPoint + 1, newPoint + 2);

   if (fSelectedPart == 4)
      fXOZSectionPos = newPoint[1];
   else if (fSelectedPart == 5)
      fYOZSectionPos = newPoint[0];
   else
      fXOYSectionPos = newPoint[2];
}

////////////////////////////////////////////////////////////////////////////////
///Draw sections (if any).

void TGLPlotPainter::DrawSections()const
{
   const TGLVertex3 *frame = fBackBox.Get3DBox();

   if (fXOZSectionPos > frame[0].Y()) {
      if (fXOZSectionPos > frame[2].Y())
         fXOZSectionPos = frame[2].Y();
      const TGLVertex3 v1(frame[0].X(), fXOZSectionPos, frame[0].Z());
      const TGLVertex3 v2(frame[4].X(), fXOZSectionPos, frame[4].Z());
      const TGLVertex3 v3(frame[5].X(), fXOZSectionPos, frame[5].Z());
      const TGLVertex3 v4(frame[1].X(), fXOZSectionPos, frame[1].Z());

      if (fSelectionPass)
         Rgl::ObjectIDToColor(4, fHighColor);
      else if (fSelectedPart == 4)
         glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, Rgl::gBlueEmission);

      glEnable(GL_POLYGON_OFFSET_FILL);
      glPolygonOffset(1.f, 1.f);
      Rgl::DrawQuadFilled(v1, v2, v3, v4, TGLVector3(0., 1., 0.));
      glDisable(GL_POLYGON_OFFSET_FILL);
      //Zlevels here.
      if (!fSelectionPass) {
         if (fSelectedPart == 4)
            glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, Rgl::gNullEmission);
         const TGLDisableGuard lightGuard(GL_LIGHTING);
         const TGLEnableGuard  blendGuard(GL_BLEND);
         const TGLEnableGuard  lineSmooth(GL_LINE_SMOOTH);
         glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
         glDepthMask(GL_FALSE);
         DrawSectionXOZ();
         //Draw z-levels
         const TGLEnableGuard stippleGuard(GL_LINE_STIPPLE);//[1-1]
         const UShort_t stipple = 0x5555;
         glLineStipple(1, stipple);

         glColor3d(0., 0., 0.);
         glBegin(GL_LINES);
         for (UInt_t i = 0; i < fZLevels.size(); ++i) {
            glVertex3d(fBackBox.Get3DBox()[1].X(), fXOZSectionPos, fZLevels[i]);
            glVertex3d(fBackBox.Get3DBox()[0].X(), fXOZSectionPos, fZLevels[i]);
         }
         glEnd();
         glDepthMask(GL_TRUE);
      }
   }

   if (fYOZSectionPos > frame[0].X()) {
      if (fYOZSectionPos > frame[1].X())
         fYOZSectionPos = frame[1].X();
      TGLVertex3 v1(fYOZSectionPos, frame[0].Y(), frame[0].Z());
      TGLVertex3 v2(fYOZSectionPos, frame[3].Y(), frame[3].Z());
      TGLVertex3 v3(fYOZSectionPos, frame[7].Y(), frame[7].Z());
      TGLVertex3 v4(fYOZSectionPos, frame[4].Y(), frame[4].Z());

      if (fSelectionPass) {
         Rgl::ObjectIDToColor(5, fHighColor);
      } else if (fSelectedPart == 5)
         glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, Rgl::gBlueEmission);

      glEnable(GL_POLYGON_OFFSET_FILL);
      glPolygonOffset(1.f, 1.f);
      Rgl::DrawQuadFilled(v1, v2, v3, v4, TGLVector3(1., 0., 0.));
      glDisable(GL_POLYGON_OFFSET_FILL);

      if (!fSelectionPass) {
         if (fSelectedPart == 5)
            glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, Rgl::gNullEmission);
         const TGLDisableGuard lightHuard(GL_LIGHTING);
         const TGLEnableGuard blendGuard(GL_BLEND);
         const TGLEnableGuard lineSmooth(GL_LINE_SMOOTH);
         glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
         glDepthMask(GL_FALSE);
         DrawSectionYOZ();
         //Draw z-levels
         const TGLEnableGuard stippleGuard(GL_LINE_STIPPLE);//[1-1]
         glLineStipple(1, 0x5555);

         glColor3d(0., 0., 0.);
         glBegin(GL_LINES);
         for (UInt_t i = 0; i < fZLevels.size(); ++i) {
            glVertex3d(fYOZSectionPos, fBackBox.Get3DBox()[3].Y(), fZLevels[i]);
            glVertex3d(fYOZSectionPos, fBackBox.Get3DBox()[0].Y(), fZLevels[i]);
         }
         glEnd();
         glDepthMask(GL_TRUE);
      }
   }

   if (fXOYSectionPos > frame[0].Z()) {
      if (fXOYSectionPos > frame[4].Z())
         fXOYSectionPos = frame[4].Z();
      TGLVertex3 v1(frame[0].X(), frame[0].Y(), fXOYSectionPos);
      TGLVertex3 v2(frame[1].X(), frame[1].Y(), fXOYSectionPos);
      TGLVertex3 v3(frame[2].X(), frame[2].Y(), fXOYSectionPos);
      TGLVertex3 v4(frame[3].X(), frame[3].Y(), fXOYSectionPos);

      if (fSelectionPass) {
         Rgl::ObjectIDToColor(6, fHighColor);
      } else if (fSelectedPart == 6)
         glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, Rgl::gBlueEmission);

      glEnable(GL_POLYGON_OFFSET_FILL);
      glPolygonOffset(1.f, 1.f);
      //if (fSelectionPass || fSelectedPart == 6)
      Rgl::DrawQuadFilled(v1, v2, v3, v4, TGLVector3(0., 0., 1.));
      glDisable(GL_POLYGON_OFFSET_FILL);

      if (!fSelectionPass) {
         if (fSelectedPart == 6)
            glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, Rgl::gNullEmission);
         const TGLDisableGuard lightGuard(GL_LIGHTING);
         const TGLEnableGuard blendGuard(GL_BLEND);
         const TGLEnableGuard lineSmooth(GL_LINE_SMOOTH);
         glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
         glDepthMask(GL_FALSE);
         DrawSectionXOY();
         glDepthMask(GL_TRUE);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////

void TGLPlotPainter::ClearBuffers()const
{
/*
   // Clear buffer.
   Float_t rgb[3] = {1.f, 1.f, 1.f};
   if (const TColor *color = GetPadColor())
      color->GetRGB(rgb[0], rgb[1], rgb[2]);
   glClearColor(rgb[0], rgb[1], rgb[2], 1.);
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
   */
}

////////////////////////////////////////////////////////////////////////////////
///Draw. Palette. Axis.

void TGLPlotPainter::DrawPaletteAxis()const
{
}

////////////////////////////////////////////////////////////////////////////////

void TGLPlotPainter::SaveModelviewMatrix()const
{
   glMatrixMode(GL_MODELVIEW);
   glPushMatrix();
}

////////////////////////////////////////////////////////////////////////////////

void TGLPlotPainter::SaveProjectionMatrix()const
{
   glMatrixMode(GL_PROJECTION);
   glPushMatrix();
}

////////////////////////////////////////////////////////////////////////////////

void TGLPlotPainter::RestoreModelviewMatrix()const
{
   glMatrixMode(GL_MODELVIEW);
   glPopMatrix();
}

////////////////////////////////////////////////////////////////////////////////

void TGLPlotPainter::RestoreProjectionMatrix()const
{
   glMatrixMode(GL_PROJECTION);
   glPopMatrix();
}

/** \class TGLPlotCoordinates
\ingroup opengl
Helper class for plot-painters holding information about axis
ranges, numbers of bins and flags if certain axis is logarithmic.
*/

ClassImp(TGLPlotCoordinates);

////////////////////////////////////////////////////////////////////////////////
///Constructor.

TGLPlotCoordinates::TGLPlotCoordinates()
                        : fCoordType(kGLCartesian),
                          fXScale(1.),
                          fYScale(1.),
                          fZScale(1.),
                          fXLog(kFALSE),
                          fYLog(kFALSE),
                          fZLog(kFALSE),
                          fModified(kFALSE),
                          fFactor(1.)
{
}

////////////////////////////////////////////////////////////////////////////////
///Destructor.

TGLPlotCoordinates::~TGLPlotCoordinates()
{
}

////////////////////////////////////////////////////////////////////////////////
///If coord type was changed, plot must reset sections (if any),
///set fModified.

void TGLPlotCoordinates::SetCoordType(EGLCoordType type)
{
   if (fCoordType != type) {
      fModified = kTRUE;
      fCoordType = type;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get coordinates type.

EGLCoordType TGLPlotCoordinates::GetCoordType()const
{
   return fCoordType;
}

////////////////////////////////////////////////////////////////////////////////
///If log changed, sections must be reset,
///set fModified.

void TGLPlotCoordinates::SetXLog(Bool_t xLog)
{
   if (fXLog != xLog) {
      fXLog = xLog;
      fModified = kTRUE;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get X log.

Bool_t TGLPlotCoordinates::GetXLog()const
{
   return fXLog;
}

////////////////////////////////////////////////////////////////////////////////
///If log changed, sections must be reset,
///set fModified.

void TGLPlotCoordinates::SetYLog(Bool_t yLog)
{
   if (fYLog != yLog) {
      fYLog = yLog;
      fModified = kTRUE;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get Y log.

Bool_t TGLPlotCoordinates::GetYLog()const
{
   return fYLog;
}

////////////////////////////////////////////////////////////////////////////////
///If log changed, sections must be reset,
///set fModified.

void TGLPlotCoordinates::SetZLog(Bool_t zLog)
{
   if (fZLog != zLog) {
      fZLog = zLog;
      fModified = kTRUE;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get Z log.

Bool_t TGLPlotCoordinates::GetZLog()const
{
   return fZLog;
}

////////////////////////////////////////////////////////////////////////////////
/// Reset modified.

void TGLPlotCoordinates::ResetModified()
{
   fModified = !fModified;//kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Modified.

Bool_t TGLPlotCoordinates::Modified()const
{
   return fModified;
}

////////////////////////////////////////////////////////////////////////////////
///Set bin ranges, ranges.

Bool_t TGLPlotCoordinates::SetRanges(const TH1 *hist, Bool_t errors, Bool_t zBins)
{
   switch (fCoordType) {
   case kGLPolar:
      return SetRangesPolar(hist);
   case kGLCylindrical:
      return SetRangesCylindrical(hist);
   case kGLSpherical:
      return SetRangesSpherical(hist);
   case kGLCartesian:
   default:
      return SetRangesCartesian(hist, errors, zBins);
   }
}

////////////////////////////////////////////////////////////////////////////////
///Number of X bins.

Int_t TGLPlotCoordinates::GetNXBins()const
{
   return fXBins.second - fXBins.first + 1;
}

////////////////////////////////////////////////////////////////////////////////
///Number of Y bins.

Int_t TGLPlotCoordinates::GetNYBins()const
{
   return fYBins.second - fYBins.first + 1;
}

////////////////////////////////////////////////////////////////////////////////
///Number of Z bins.

Int_t TGLPlotCoordinates::GetNZBins()const
{
   return fZBins.second - fZBins.first + 1;
}

////////////////////////////////////////////////////////////////////////////////
///X bins range.

const Rgl::BinRange_t &TGLPlotCoordinates::GetXBins()const
{
   return fXBins;
}

////////////////////////////////////////////////////////////////////////////////
///Y bins range.

const Rgl::BinRange_t &TGLPlotCoordinates::GetYBins()const
{
   return fYBins;
}

////////////////////////////////////////////////////////////////////////////////
///Z bins range.

const Rgl::BinRange_t &TGLPlotCoordinates::GetZBins()const
{
   return fZBins;
}

////////////////////////////////////////////////////////////////////////////////
///X range.

const Rgl::Range_t &TGLPlotCoordinates::GetXRange()const
{
   return fXRange;
}

////////////////////////////////////////////////////////////////////////////////
///X length.

Double_t TGLPlotCoordinates::GetXLength()const
{
   return fXRange.second - fXRange.first;
}

////////////////////////////////////////////////////////////////////////////////
///Y range.

const Rgl::Range_t &TGLPlotCoordinates::GetYRange()const
{
   return fYRange;
}

////////////////////////////////////////////////////////////////////////////////
///Y length.

Double_t TGLPlotCoordinates::GetYLength()const
{
   return fYRange.second - fYRange.first;
}


////////////////////////////////////////////////////////////////////////////////
///Z range.

const Rgl::Range_t &TGLPlotCoordinates::GetZRange()const
{
   return fZRange;
}

////////////////////////////////////////////////////////////////////////////////
///Z length.

Double_t TGLPlotCoordinates::GetZLength()const
{
   return fZRange.second - fZRange.first;
}


////////////////////////////////////////////////////////////////////////////////
///Scaled range.

const Rgl::Range_t &TGLPlotCoordinates::GetXRangeScaled()const
{
   return fXRangeScaled;
}

////////////////////////////////////////////////////////////////////////////////
///Scaled range.

const Rgl::Range_t &TGLPlotCoordinates::GetYRangeScaled()const
{
   return fYRangeScaled;
}

////////////////////////////////////////////////////////////////////////////////
///Scaled range.

const Rgl::Range_t &TGLPlotCoordinates::GetZRangeScaled()const
{
   return fZRangeScaled;
}

////////////////////////////////////////////////////////////////////////////////
/// Get factor.

Double_t TGLPlotCoordinates::GetFactor()const
{
   return fFactor;
}

namespace {

Bool_t FindAxisRange(const TAxis *axis, Bool_t log, Rgl::BinRange_t &bins, Rgl::Range_t &range);
Bool_t FindAxisRange(const TH1 *hist, Bool_t logZ, const Rgl::BinRange_t &xBins,
                     const Rgl::BinRange_t &yBins, Rgl::Range_t &zRange,
                     Double_t &factor, Bool_t errors);

Bool_t FindAxisRange(TH2Poly *hist, Bool_t zLog, Rgl::Range_t &zRange);

}

////////////////////////////////////////////////////////////////////////////////
///Set bin ranges, ranges, etc.

Bool_t TGLPlotCoordinates::SetRangesCartesian(const TH1 *hist, Bool_t errors, Bool_t zAsBins)
{
   Rgl::BinRange_t xBins;
   Rgl::Range_t    xRange;
   if (!FindAxisRange(hist->GetXaxis(), fXLog, xBins, xRange)) {
      Error("TGLPlotCoordinates::SetRangesCartesian", "Cannot set X axis to log scale");
      return kFALSE;
   }

   Rgl::BinRange_t yBins;
   Rgl::Range_t    yRange;
   if (!FindAxisRange(hist->GetYaxis(), fYLog, yBins, yRange)) {
      Error("TGLPlotCoordinates::SetRangesCartesian", "Cannot set Y axis to log scale");
      return kFALSE;
   }

   Rgl::BinRange_t zBins;
   Rgl::Range_t zRange;
   Double_t factor = 1.;

   if (zAsBins) {
      if (!FindAxisRange(hist->GetZaxis(), fZLog, zBins, zRange)) {
         Error("TGLPlotCoordinates::SetRangesCartesian", "Cannot set Z axis to log scale");
         return kFALSE;
      }
   } else if (!FindAxisRange(hist, fZLog, xBins, yBins, zRange, factor, errors)) {
      Error("TGLPlotCoordinates::SetRangesCartesian",
            "Log scale is requested for Z, but maximum less or equal 0. (%f)", zRange.second);
      return kFALSE;
   }

   //Finds the maximum dimension and adjust scale coefficients
   const Double_t x = xRange.second - xRange.first;
   const Double_t y = yRange.second - yRange.first;
   const Double_t z = zRange.second - zRange.first;

   if (!x || !y || !z) {
      Error("TGLPlotCoordinates::SetRangesCartesian", "Zero axis range.");
      return kFALSE;
   }

   if (xRange != fXRange || yRange != fYRange || zRange != fZRange ||
       xBins != fXBins || yBins != fYBins || zBins != fZBins || fFactor != factor)
   {
      fModified = kTRUE;
   }

   fXRange = xRange, fXBins = xBins, fYRange = yRange, fYBins = yBins, fZRange = zRange, fZBins = zBins;
   fFactor = factor;
   /*
   const Double_t maxDim = TMath::Max(TMath::Max(x, y), z);
   fXScale = maxDim / x;
   fYScale = maxDim / y;
   fZScale = maxDim / z;
   */
   fXScale = 1. / x;
   fYScale = 1. / y;
   fZScale = 1. / z;

   fXRangeScaled.first = fXRange.first * fXScale, fXRangeScaled.second = fXRange.second * fXScale;
   fYRangeScaled.first = fYRange.first * fYScale, fYRangeScaled.second = fYRange.second * fYScale;
   fZRangeScaled.first = fZRange.first * fZScale, fZRangeScaled.second = fZRange.second * fZScale;

   return kTRUE;
}

//
////////////////////////////////////////////////////////////////////////////////
///Set bin ranges, ranges, etc.

Bool_t TGLPlotCoordinates::SetRanges(TH2Poly *hist)
{
   Rgl::BinRange_t xBins;
   Rgl::Range_t    xRange;
   FindAxisRange(hist->GetXaxis(), kFALSE, xBins, xRange);//kFALSE == never logarithmic.

   Rgl::BinRange_t yBins;
   Rgl::Range_t    yRange;
   FindAxisRange(hist->GetYaxis(), kFALSE, yBins, yRange);//kFALSE == never logarithmic.

   Rgl::BinRange_t zBins;
   Rgl::Range_t zRange;
   Double_t factor = 1.;

   if (!FindAxisRange(hist, fZLog, zRange))
      return kFALSE;

   //Finds the maximum dimension and adjust scale coefficients
   const Double_t x = xRange.second - xRange.first;
   const Double_t y = yRange.second - yRange.first;
   const Double_t z = zRange.second - zRange.first;

   if (!x || !y || !z) {
      Error("TGLPlotCoordinates::SetRanges", "Zero axis range.");
      return kFALSE;
   }

   if (xRange != fXRange || yRange != fYRange || zRange != fZRange ||
       xBins != fXBins || yBins != fYBins || zBins != fZBins || fFactor != factor)
   {
      fModified = kTRUE;
   }

   fXRange = xRange, fXBins = xBins, fYRange = yRange, fYBins = yBins, fZRange = zRange, fZBins = zBins;
   fFactor = factor;

   fXScale = Rgl::gH2PolyScaleXY / x;
   fYScale = Rgl::gH2PolyScaleXY / y;
   fZScale = 1. / z;

   fXRangeScaled.first = fXRange.first * fXScale, fXRangeScaled.second = fXRange.second * fXScale;
   fYRangeScaled.first = fYRange.first * fYScale, fYRangeScaled.second = fYRange.second * fYScale;
   fZRangeScaled.first = fZRange.first * fZScale, fZRangeScaled.second = fZRange.second * fZScale;

   return kTRUE;
}
//

////////////////////////////////////////////////////////////////////////////////
///Set bin ranges, ranges, etc.

Bool_t TGLPlotCoordinates::SetRanges(const TAxis *xAxis, const TAxis *yAxis, const TAxis *zAxis)
{
   Rgl::BinRange_t xBins;
   Rgl::Range_t    xRange;

   FindAxisRange(xAxis, kFALSE, xBins, xRange);

   Rgl::BinRange_t yBins;
   Rgl::Range_t    yRange;

   FindAxisRange(yAxis, kFALSE, yBins, yRange);

   Rgl::BinRange_t zBins;
   Rgl::Range_t zRange;
   Double_t factor = 1.;

   FindAxisRange(zAxis, kFALSE, zBins, zRange);

   //Finds the maximum dimension and adjust scale coefficients
   const Double_t x = xRange.second - xRange.first;
   const Double_t y = yRange.second - yRange.first;
   const Double_t z = zRange.second - zRange.first;

   if (!x || !y || !z) {
      Error("TGLPlotCoordinates::SetRangesCartesian", "Zero axis range.");
      return kFALSE;
   }

   if (xRange != fXRange || yRange != fYRange || zRange != fZRange ||
       xBins != fXBins || yBins != fYBins || zBins != fZBins || fFactor != factor)
   {
      fModified = kTRUE;
   }

   fXRange = xRange, fXBins = xBins, fYRange = yRange, fYBins = yBins, fZRange = zRange, fZBins = zBins;
   fFactor = factor;

/*   const Double_t maxDim = TMath::Max(TMath::Max(x, y), z);
   fXScale = maxDim / x;
   fYScale = maxDim / y;
   fZScale = maxDim / z;*/

   fXScale = 1. / x;
   fYScale = 1. / y;
   fZScale = 1. / z;

   fXRangeScaled.first = fXRange.first * fXScale, fXRangeScaled.second = fXRange.second * fXScale;
   fYRangeScaled.first = fYRange.first * fYScale, fYRangeScaled.second = fYRange.second * fYScale;
   fZRangeScaled.first = fZRange.first * fZScale, fZRangeScaled.second = fZRange.second * fZScale;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
///Set bin ranges, ranges, etc.

Bool_t TGLPlotCoordinates::SetRangesPolar(const TH1 *hist)
{
   Rgl::BinRange_t xBins;
   Rgl::Range_t phiRange;
   const TAxis *xAxis = hist->GetXaxis();
   FindAxisRange(xAxis, kFALSE, xBins, phiRange);
   if (xBins.second - xBins.first + 1 > 360) {
      Error("TGLPlotCoordinates::SetRangesPolar", "To many PHI sectors");
      return kFALSE;
   }

   Rgl::BinRange_t yBins;
   Rgl::Range_t roRange;
   const TAxis *yAxis = hist->GetYaxis();
   FindAxisRange(yAxis, kFALSE, yBins, roRange);

   Rgl::Range_t zRange;
   Double_t factor = 1.;
   if (!FindAxisRange(hist, fZLog, xBins, yBins, zRange, factor, kFALSE))
   {
      Error("TGLPlotCoordinates::SetRangesPolar",
            "Log scale is requested for Z, but maximum less or equal 0. (%f)", zRange.second);
      return kFALSE;
   }

   const Double_t z = zRange.second - zRange.first;
   if (!z || !(phiRange.second - phiRange.first) || !(roRange.second - roRange.first)) {
      Error("TGLPlotCoordinates::SetRangesPolar", "Zero axis range.");
      return kFALSE;
   }

   if (phiRange != fXRange || roRange != fYRange || zRange != fZRange ||
       xBins != fXBins || yBins != fYBins || fFactor != factor)
   {
      fModified = kTRUE;
      fXRange = phiRange, fXBins = xBins;
      fYRange = roRange,  fYBins = yBins;
      fZRange = zRange;
      fFactor = factor;
   }

   //const Double_t maxDim = TMath::Max(2., z);
   fXScale = 0.5;//maxDim / 2.;
   fYScale = 0.5;//maxDim / 2.;
   fZScale = 1. / z;//maxDim / z;
   fXRangeScaled.first = -fXScale, fXRangeScaled.second = fXScale;
   fYRangeScaled.first = -fYScale, fYRangeScaled.second = fYScale;
   fZRangeScaled.first = fZRange.first * fZScale, fZRangeScaled.second = fZRange.second * fZScale;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set ranges cylindrical.

Bool_t TGLPlotCoordinates::SetRangesCylindrical(const TH1 *hist)
{
   Rgl::BinRange_t xBins, yBins;
   Rgl::Range_t angleRange, heightRange, radiusRange;
   const TAxis *xAxis = hist->GetXaxis();
   const TAxis *yAxis = hist->GetYaxis();
   Double_t factor = 1.;

   FindAxisRange(xAxis, kFALSE, xBins, angleRange);
   if (xBins.second - xBins.first + 1 > 360) {
      Error("TGLPlotCoordinates::SetRangesCylindrical", "To many PHI sectors");
      return kFALSE;
   }
   if (!FindAxisRange(yAxis, fYLog, yBins, heightRange)) {
      Error("TGLPlotCoordinates::SetRangesCylindrical", "Cannot set Y axis to log scale");
      return kFALSE;
   }
   FindAxisRange(hist, kFALSE, xBins, yBins, radiusRange, factor, kFALSE);

   const Double_t x = angleRange.second  - angleRange.first;
   const Double_t y = heightRange.second - heightRange.first;
   const Double_t z = radiusRange.second - radiusRange.first;

   if (!x || !y || !z) {
      Error("TGLPlotCoordinates::SetRangesCylindrical", "Zero axis range.");
      return kFALSE;
   }

   if (angleRange != fXRange  || heightRange != fYRange ||
       radiusRange != fZRange || xBins != fXBins ||
       yBins != fYBins || fFactor != factor)
   {
      fModified = kTRUE;
      fXRange = angleRange,  fXBins = xBins;
      fYRange = heightRange, fYBins = yBins;
      fZRange = radiusRange;
      fFactor = factor;
   }

   // const Double_t maxDim = TMath::Max(2., y);
   fXScale = 0.5;//maxDim / 2.;
   fYScale = 1. / y;//maxDim / y;
   fZScale = 0.5;//maxDim / 2.;
   fXRangeScaled.first = -fXScale, fXRangeScaled.second = fXScale;
   fYRangeScaled.first = fYRange.first * fYScale, fYRangeScaled.second = fYRange.second * fYScale;
   fZRangeScaled.first = -fZScale, fZRangeScaled.second = fZScale;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set ranges spherical.

Bool_t TGLPlotCoordinates::SetRangesSpherical(const TH1 *hist)
{
   Rgl::BinRange_t xBins;
   Rgl::Range_t phiRange;
   FindAxisRange(hist->GetXaxis(), kFALSE, xBins, phiRange);
   if (xBins.second - xBins.first + 1 > 360) {
      Error("TGLPlotCoordinates::SetRangesSpherical", "To many PHI sectors");
      return kFALSE;
   }

   Rgl::BinRange_t yBins;
   Rgl::Range_t thetaRange;
   FindAxisRange(hist->GetYaxis(), kFALSE, yBins, thetaRange);
   if (yBins.second - yBins.first + 1 > 180) {
      Error("TGLPlotCoordinates::SetRangesSpherical", "To many THETA sectors");
      return kFALSE;
   }

   Rgl::Range_t radiusRange;
   Double_t factor = 1.;
   FindAxisRange(hist, kFALSE, xBins, yBins, radiusRange, factor, kFALSE);

   if (xBins != fXBins || yBins != fYBins ||
       phiRange != fXRange || thetaRange != fYRange ||
       radiusRange != fZRange || fFactor != factor)
   {
      fModified = kTRUE;
      fXBins    = xBins;
      fYBins    = yBins;
      fXRange   = phiRange;
      fYRange   = thetaRange,
      fZRange   = radiusRange;
      fFactor   = factor;
   }

   fXScale = 0.5;
   fYScale = 0.5;
   fZScale = 0.5;
   fXRangeScaled.first = -fXScale, fXRangeScaled.second = fXScale;
   fYRangeScaled.first = -fYScale, fYRangeScaled.second = fYScale;
   fZRangeScaled.first = -fZScale, fZRangeScaled.second = fZScale;

   return kTRUE;
}

namespace {

////////////////////////////////////////////////////////////////////////////////
/// Find minimal bin width.

Double_t FindMinBinWidth(const TAxis *axis)
{
   Int_t currBin = axis->GetFirst();
   Double_t width = axis->GetBinWidth(currBin);

   if (!axis->IsVariableBinSize())//equal bins
      return width;

   ++currBin;
   //variable size bins
   for (const Int_t lastBin = axis->GetLast(); currBin <= lastBin; ++currBin)
      width = TMath::Min(width, axis->GetBinWidth(currBin));

   return width;
}

////////////////////////////////////////////////////////////////////////////////
///"Generic" function, can be used for X/Y/Z axis.
///[low edge of first ..... up edge of last]
///If log is true, at least up edge of last MUST be positive or function fails (1).
///If log is true and low edge is negative, try to find bin with positive low edge, bin number
///must be less or equal to last (2). If no such bin, function fails.
///When looking for a such bin, I'm trying to find value which is 0.01 of
///MINIMUM bin width (3) (if bins are equidimensional, first's bin width is OK).
///But even such lookup can fail, so, it's a stupid idea to have negative ranges
///and logarithmic scale :)

Bool_t FindAxisRange(const TAxis *axis, Bool_t log, Rgl::BinRange_t &bins, Rgl::Range_t &range)
{
   bins.first = axis->GetFirst(), bins.second = axis->GetLast();
   range.first = axis->GetBinLowEdge(bins.first), range.second = axis->GetBinUpEdge(bins.second);

   if (log) {
      if (range.second <= 0.)
         return kFALSE;//(1)

      range.second = TMath::Log10(range.second);

      if (range.first <= 0.) {//(2)
         Int_t bin = axis->FindFixBin(FindMinBinWidth(axis) * 0.01);//(3)
         //Overflow or something stupid.
         if (bin > bins.second)
            return kFALSE;

         if (axis->GetBinLowEdge(bin) <= 0.) {
            ++bin;
            if (bin > bins.second)//Again, something stupid.
               return kFALSE;
         }

         bins.first = bin;
         range.first = axis->GetBinLowEdge(bin);
      }

      range.first = TMath::Log10(range.first);
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
///First, look through hist to find minimum and maximum values.

Bool_t FindAxisRange(const TH1 *hist, Bool_t logZ, const Rgl::BinRange_t &xBins,
                     const Rgl::BinRange_t &yBins, Rgl::Range_t &zRange,
                     Double_t &factor, Bool_t errors)
{
   const Bool_t minimum = hist->GetMinimumStored() != -1111;
   const Bool_t maximum = hist->GetMaximumStored() != -1111;
   const Double_t margin = gStyle->GetHistTopMargin();

   zRange.second = hist->GetBinContent(hist->GetBin(xBins.first, yBins.first)), zRange.first = zRange.second;
   Double_t summ = 0.;

   for (Int_t i = xBins.first; i <= xBins.second; ++i) {
      for (Int_t j = yBins.first; j <= yBins.second; ++j) {
         Double_t val = hist->GetBinContent(hist->GetBin(i, j));
         if (val > 0. && errors)
            val = TMath::Max(val, val + hist->GetCellError(i, j));
         zRange.second = TMath::Max(val, zRange.second);
         zRange.first = TMath::Min(val, zRange.first);
         summ += val;
      }
   }

   if (hist->GetMaximumStored() != -1111)
      zRange.second = hist->GetMaximumStored();
   if (hist->GetMinimumStored() != -1111)
      zRange.first = hist->GetMinimumStored();

   if (logZ && zRange.second <= 0.)
      return kFALSE;//cannot setup logarithmic scale

   if (zRange.first >= zRange.second)
      zRange.first = 0.001 * zRange.second;

   factor = hist->GetNormFactor() > 0. ? hist->GetNormFactor() : summ;
   if (summ) factor /= summ;
   if (!factor) factor = 1.;
   if (factor < 0.)
      Warning("TGLPlotPainter::ExtractAxisZInfo",
              "Negative factor, negative ranges - possible incorrect behavior");

   zRange.second *= factor;
   zRange.first  *= factor;

   if (logZ) {
      if (zRange.first <= 0.)
         zRange.first = TMath::Min(1., 0.001 * zRange.second);
      zRange.first = TMath::Log10(zRange.first);
      if (!minimum)
         zRange.first += TMath::Log10(0.5);
      zRange.second = TMath::Log10(zRange.second);
      if (!maximum)
         zRange.second += TMath::Log10(2*(0.9/0.95));//This magic numbers are from THistPainter.
      return kTRUE;
   }

   if (!maximum)
      zRange.second += margin * (zRange.second - zRange.first);
   if (!minimum) {
      if (gStyle->GetHistMinimumZero())
         zRange.first >= 0 ? zRange.first = 0. : zRange.first -= margin * (zRange.second - zRange.first);
      else
         zRange.first >= 0 && zRange.first - margin * (zRange.second - zRange.first) <= 0 ?
               zRange.first = 0 : zRange.first -= margin * (zRange.second - zRange.first);
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
///First, look through hist to find minimum and maximum values.

Bool_t FindAxisRange(TH2Poly *hist, Bool_t logZ, Rgl::Range_t &zRange)
{
   TList *bins = hist->GetBins();
   if (!bins || !bins->GetEntries()) {
      Error("FindAxisRange", "TH2Poly returned empty list of bins");
      return kFALSE;
   }

   zRange.first  = hist->GetMinimum();
   zRange.second = hist->GetMaximum();

   if (zRange.first >= zRange.second)
      zRange.first = 0.001 * zRange.second;

   if (logZ) {
      if (zRange.second < 1e-20) {
         Error("FindAxisRange", "Failed to switch Z axis to logarithmic scale");
         return kFALSE;
      }

      if (zRange.first <= 0.)
         zRange.first = TMath::Min(1., 0.001 * zRange.second);

      zRange.first  = TMath::Log10(zRange.first);
      zRange.first += TMath::Log10(0.5);
      zRange.second = TMath::Log10(zRange.second);
      zRange.second += TMath::Log10(2 * (0.9 / 0.95));//These magic numbers come from THistPainter.

      return kTRUE;
   }

   const Double_t margin = gStyle->GetHistTopMargin();
   zRange.second += margin * (zRange.second - zRange.first);
   if (gStyle->GetHistMinimumZero())
      zRange.first >= 0 ? zRange.first = 0. : zRange.first -= margin * (zRange.second - zRange.first);
   else
      zRange.first >= 0 && zRange.first - margin * (zRange.second - zRange.first) <= 0 ?
         zRange.first = 0 : zRange.first -= margin * (zRange.second - zRange.first);

   return kTRUE;
}

}

/** \class TGLBoxCut
\ingroup opengl
Used by plot-painters to determine the area of the plot that is cut away.
*/

ClassImp(TGLBoxCut);

////////////////////////////////////////////////////////////////////////////////
///Constructor.

TGLBoxCut::TGLBoxCut(const TGLPlotBox *box)
               : fXLength(0.),
                 fYLength(0.),
                 fZLength(0.),
                 fPlotBox(box),
                 fActive(kFALSE),
                 fFactor(1.)
{
}

////////////////////////////////////////////////////////////////////////////////
///Destructor.

TGLBoxCut::~TGLBoxCut()
{
}

////////////////////////////////////////////////////////////////////////////////
///Turn the box cut on/off.
///If it's on, it will be placed in front point of a plot.

void TGLBoxCut::TurnOnOff()
{
   fActive = !fActive;

   if (fActive) {
      ResetBoxGeometry();
   }
}

////////////////////////////////////////////////////////////////////////////////
///Turn the box cut on/off.

void TGLBoxCut::SetActive(Bool_t a)
{
   if (a == fActive)
      return;
   TurnOnOff();
}

////////////////////////////////////////////////////////////////////////////////
///Set geometry using plot's back box.

void TGLBoxCut::ResetBoxGeometry()
{
   const Int_t frontPoint = fPlotBox->GetFrontPoint();
   const TGLVertex3 *box = fPlotBox->Get3DBox();
   const TGLVertex3 center((box[0].X() + box[1].X()) / 2, (box[0].Y() + box[2].Y()) / 2,
                           (box[0].Z() + box[4].Z()) / 2);
   fXLength = fFactor * (box[1].X() - box[0].X());
   fYLength = fFactor * (box[2].Y() - box[0].Y());
   fZLength = fFactor * (box[4].Z() - box[0].Z());

   switch(frontPoint){
   case 0:
      fCenter.X() = box[0].X();
      fCenter.Y() = box[0].Y();
      break;
   case 1:
      fCenter.X() = box[1].X();
      fCenter.Y() = box[0].Y();
      break;
   case 2:
      fCenter.X() = box[2].X();
      fCenter.Y() = box[2].Y();
      break;
   case 3:
      fCenter.X() = box[0].X();
      fCenter.Y() = box[2].Y();
      break;
   }

   fCenter.Z() = box[0].Z() * 0.5 + box[4].Z() * 0.5;
   AdjustBox();
}

////////////////////////////////////////////////////////////////////////////////
///Draw cut as a semi-transparent box.

void TGLBoxCut::DrawBox(Bool_t selectionPass, Int_t selected)const
{
   if (!selectionPass) {
      glDisable(GL_LIGHTING);
      glLineWidth(3.f);

      selected == TGLPlotPainter::kXAxis ? glColor3d(1., 1., 0.) : glColor3d(1., 0., 0.);
      glBegin(GL_LINES);
      glVertex3d(fXRange.first, (fYRange.first + fYRange.second) / 2, (fZRange.first + fZRange.second) / 2);
      glVertex3d(fXRange.second, (fYRange.first + fYRange.second) / 2, (fZRange.first + fZRange.second) / 2);
      glEnd();

      selected == TGLPlotPainter::kYAxis ? glColor3d(1., 1., 0.) : glColor3d(0., 1., 0.);
      glBegin(GL_LINES);
      glVertex3d((fXRange.first + fXRange.second) / 2, fYRange.first, (fZRange.first + fZRange.second) / 2);
      glVertex3d((fXRange.first + fXRange.second) / 2, fYRange.second, (fZRange.first + fZRange.second) / 2);
      glEnd();

      selected == TGLPlotPainter::kZAxis ? glColor3d(1., 1., 0.) : glColor3d(0., 0., 1.);
      glBegin(GL_LINES);
      glVertex3d((fXRange.first + fXRange.second) / 2, (fYRange.first + fYRange.second) / 2, fZRange.first);
      glVertex3d((fXRange.first + fXRange.second) / 2, (fYRange.first + fYRange.second) / 2, fZRange.second);
      glEnd();

      glLineWidth(1.f);
      glEnable(GL_LIGHTING);

      GLboolean oldBlendState = kFALSE;
      glGetBooleanv(GL_BLEND, &oldBlendState);

      if (!oldBlendState)
         glEnable(GL_BLEND);

      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


      const Float_t diffuseColor[] = {0.f, 0.f, 1.f, 0.1f};
      glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuseColor);

      Rgl::DrawBoxFront(fXRange.first, fXRange.second, fYRange.first, fYRange.second,
                        fZRange.first, fZRange.second, fPlotBox->GetFrontPoint());

      if (!oldBlendState)
         glDisable(GL_BLEND);
   } else {
      glLineWidth(5.f);
      Rgl::ObjectIDToColor(TGLPlotPainter::kXAxis, kFALSE);
      glBegin(GL_LINES);
      glVertex3d(fXRange.first, (fYRange.first + fYRange.second) / 2, (fZRange.first + fZRange.second) / 2);
      glVertex3d(fXRange.second, (fYRange.first + fYRange.second) / 2, (fZRange.first + fZRange.second) / 2);
      glEnd();

      Rgl::ObjectIDToColor(TGLPlotPainter::kYAxis, kFALSE);
      glBegin(GL_LINES);
      glVertex3d((fXRange.first + fXRange.second) / 2, fYRange.first, (fZRange.first + fZRange.second) / 2);
      glVertex3d((fXRange.first + fXRange.second) / 2, fYRange.second, (fZRange.first + fZRange.second) / 2);
      glEnd();

      Rgl::ObjectIDToColor(TGLPlotPainter::kZAxis, kFALSE);
      glBegin(GL_LINES);
      glVertex3d((fXRange.first + fXRange.second) / 2, (fYRange.first + fYRange.second) / 2, fZRange.first);
      glVertex3d((fXRange.first + fXRange.second) / 2, (fYRange.first + fYRange.second) / 2, fZRange.second);
      glEnd();
      glLineWidth(1.f);
   }
}

////////////////////////////////////////////////////////////////////////////////
///Start cut's movement

void TGLBoxCut::StartMovement(Int_t px, Int_t py)
{
   fMousePos.fX = px;
   fMousePos.fY = py;
}

////////////////////////////////////////////////////////////////////////////////
///Move box cut along selected direction.

void TGLBoxCut::MoveBox(Int_t px, Int_t py, Int_t axisID)
{
   Double_t mv[16] = {0.};
   glGetDoublev(GL_MODELVIEW_MATRIX, mv);
   Double_t pr[16] = {0.};
   glGetDoublev(GL_PROJECTION_MATRIX, pr);
   Int_t vp[4] = {0};
   glGetIntegerv(GL_VIEWPORT, vp);
   Double_t winVertex[3] = {0.};

   switch(axisID){
   case TGLPlotPainter::kXAxis :
      gluProject(fCenter.X(), 0., 0., mv, pr, vp, &winVertex[0], &winVertex[1], &winVertex[2]);
      break;
   case TGLPlotPainter::kYAxis :
      gluProject(0., fCenter.Y(), 0., mv, pr, vp, &winVertex[0], &winVertex[1], &winVertex[2]);
      break;
   case TGLPlotPainter::kZAxis :
      gluProject(0., 0., fCenter.Z(), mv, pr, vp, &winVertex[0], &winVertex[1], &winVertex[2]);
      break;
   }

   winVertex[0] += px - fMousePos.fX;
   winVertex[1] += py - fMousePos.fY;
   Double_t newPoint[3] = {0.};
   gluUnProject(winVertex[0], winVertex[1], winVertex[2], mv, pr, vp,
                newPoint, newPoint + 1, newPoint + 2);

   const TGLVertex3 *box = fPlotBox->Get3DBox();

   switch(axisID){
   case TGLPlotPainter::kXAxis :
      if (newPoint[0] >= box[1].X() + 0.4 * fXLength)
         break;
      if (newPoint[0] <= box[0].X() - 0.4 * fXLength)
         break;
      fCenter.X() = newPoint[0];
      break;
   case TGLPlotPainter::kYAxis :
      if (newPoint[1] >= box[2].Y() + 0.4 * fYLength)
         break;
      if (newPoint[1] <= box[0].Y() - 0.4 * fYLength)
         break;
      fCenter.Y() = newPoint[1];
      break;
   case TGLPlotPainter::kZAxis :
      if (newPoint[2] >= box[4].Z() + 0.4 * fZLength)
         break;
      if (newPoint[2] <= box[0].Z() - 0.4 * fZLength)
         break;
      fCenter.Z() = newPoint[2];
      break;
   }

   fMousePos.fX = px;
   fMousePos.fY = py;

   AdjustBox();
}

////////////////////////////////////////////////////////////////////////////////
///Box cut is limited by plot's sizes.

void TGLBoxCut::AdjustBox()
{
   const TGLVertex3 *box = fPlotBox->Get3DBox();

   fXRange.first  = fCenter.X() - fXLength / 2.;
   fXRange.second = fCenter.X() + fXLength / 2.;
   fYRange.first  = fCenter.Y() - fYLength / 2.;
   fYRange.second = fCenter.Y() + fYLength / 2.;
   fZRange.first  = fCenter.Z() - fZLength / 2.;
   fZRange.second = fCenter.Z() + fZLength / 2.;

   fXRange.first  = TMath::Max(fXRange.first,  box[0].X());
   fXRange.first  = TMath::Min(fXRange.first,  box[1].X());
   fXRange.second = TMath::Min(fXRange.second, box[1].X());
   fXRange.second = TMath::Max(fXRange.second, box[0].X());

   fYRange.first  = TMath::Max(fYRange.first,  box[0].Y());
   fYRange.first  = TMath::Min(fYRange.first,  box[2].Y());
   fYRange.second = TMath::Min(fYRange.second, box[2].Y());
   fYRange.second = TMath::Max(fYRange.second, box[0].Y());

   fZRange.first  = TMath::Max(fZRange.first,  box[0].Z());
   fZRange.first  = TMath::Min(fZRange.first,  box[4].Z());
   fZRange.second = TMath::Min(fZRange.second, box[4].Z());
   fZRange.second = TMath::Max(fZRange.second, box[0].Z());
}

////////////////////////////////////////////////////////////////////////////////
///Check, if box defined by xmin/xmax etc. is in cut.

Bool_t TGLBoxCut::IsInCut(Double_t xMin, Double_t xMax, Double_t yMin, Double_t yMax,
                          Double_t zMin, Double_t zMax)const
{
   if (((xMin >= fXRange.first && xMin < fXRange.second) || (xMax > fXRange.first && xMax <= fXRange.second)) &&
       ((yMin >= fYRange.first && yMin < fYRange.second) || (yMax > fYRange.first && yMax <= fYRange.second)) &&
       ((zMin >= fZRange.first && zMin < fZRange.second) || (zMax > fZRange.first && zMax <= fZRange.second)))
      return kTRUE;
   return kFALSE;
}


/** \class TGLTH3Slice
\ingroup opengl
A slice of a TH3.
*/

ClassImp(TGLTH3Slice);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TGLTH3Slice::TGLTH3Slice(const TString &name, const TH3 *hist, const TGLPlotCoordinates *coord,
                         const TGLPlotBox *box, ESliceAxis axis)
               : TNamed(name, name),
                 fAxisType(axis),
                 fAxis(0),
                 fCoord(coord),
                 fBox(box),
                 fSliceWidth(1),
                 fHist(hist),
                 fF3(0)
{
   fAxis = fAxisType == kXOZ ? fHist->GetYaxis() : fAxisType == kYOZ ? fHist->GetXaxis() : fHist->GetZaxis();
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TGLTH3Slice::TGLTH3Slice(const TString &name, const TH3 *hist, const TF3 *fun, const TGLPlotCoordinates *coord,
                         const TGLPlotBox *box, ESliceAxis axis)
               : TNamed(name, name),
                 fAxisType(axis),
                 fAxis(0),
                 fCoord(coord),
                 fBox(box),
                 fSliceWidth(1),
                 fHist(hist),
                 fF3(fun)
{
   fAxis = fAxisType == kXOZ ? fHist->GetYaxis() : fAxisType == kYOZ ? fHist->GetXaxis() : fHist->GetZaxis();
}

////////////////////////////////////////////////////////////////////////////////
/// Set Slice width.

void TGLTH3Slice::SetSliceWidth(Int_t width)
{
   if (width <= 0)
      return;

   if (fAxis->GetLast() - fAxis->GetFirst() + 1 <= width)
      fSliceWidth = fAxis->GetLast() - fAxis->GetFirst() + 1;
   else
      fSliceWidth = width;
}

////////////////////////////////////////////////////////////////////////////////
/// Draw slice.

void TGLTH3Slice::DrawSlice(Double_t pos)const
{
   Int_t bin = 0;
   for (Int_t i = fAxis->GetFirst(), e = fAxis->GetLast(); i <= e; ++i) {
      if (pos >= fAxis->GetBinLowEdge(i) && pos <= fAxis->GetBinUpEdge(i)) {
         bin = i;
         break;
      }
   }

   if (bin) {
      Int_t low = 1, up = 2;
      if (bin - fSliceWidth + 1 >= fAxis->GetFirst()) {
         low = bin - fSliceWidth + 1;
         up  = bin + 1;
      } else {
         low = fAxis->GetFirst();
         up  = bin + (fSliceWidth - (bin - fAxis->GetFirst() + 1)) + 1;
      }

      if (!fF3)
         FindMinMax(low, up);

      if (!PreparePalette())
         return;

      PrepareTexCoords(pos, low, up);

      fPalette.EnableTexture(GL_REPLACE);
      const TGLDisableGuard lightGuard(GL_LIGHTING);
      DrawSliceTextured(pos);
      fPalette.DisableTexture();
      //highlight bins in a slice.
      //DrawSliceFrame(low, up);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Find minimum and maximum for slice.

void TGLTH3Slice::FindMinMax(Int_t /*low*/, Int_t /*up*/)const
{
  /* fMinMax.first = 0.;

   switch (fAxisType) {
   case kXOZ:
      for (Int_t level = low; level < up; ++ level)
         fMinMax.first += fHist->GetBinContent(fCoord->GetFirstXBin(), level, fCoord->GetFirstZBin());
      fMinMax.second = fMinMax.first;
      for (Int_t j = fCoord->GetFirstZBin(), jt = 0, ej = fCoord->GetLastZBin(); j <= ej; ++j, ++jt) {
         for (Int_t i = fCoord->GetFirstXBin(), it = 0, ei = fCoord->GetLastXBin(); i <= ei; ++i, ++it) {
            Double_t val = 0.;
            for (Int_t level = low; level < up; ++ level)
               val += fHist->GetBinContent(i, level, j);
            fMinMax.second = TMath::Max(fMinMax.second, val);
            fMinMax.first = TMath::Min(fMinMax.first, val);
         }
      }
      break;
   case kYOZ:
      for (Int_t level = low; level < up; ++ level)
         fMinMax.first += fHist->GetBinContent(level, fCoord->GetFirstYBin(), fCoord->GetFirstZBin());
      fMinMax.second = fMinMax.first;
      for (Int_t j = fCoord->GetFirstZBin(), jt = 0, ej = fCoord->GetLastZBin(); j <= ej; ++j, ++jt) {
         for (Int_t i = fCoord->GetFirstYBin(), it = 0, ei = fCoord->GetLastYBin(); i <= ei; ++i, ++it) {
            Double_t val = 0.;
            for (Int_t level = low; level < up; ++ level)
               val += fHist->GetBinContent(level, i, j);
            fMinMax.second = TMath::Max(fMinMax.second, val);
            fMinMax.first = TMath::Min(fMinMax.first, val);
         }
      }
      break;
   case kXOY:
      for (Int_t level = low; level < up; ++ level)
         fMinMax.first += fHist->GetBinContent(fCoord->GetFirstXBin(), fCoord->GetFirstYBin(), level);
      fMinMax.second = fMinMax.first;
      for (Int_t i = fCoord->GetFirstXBin(), ir = 0, ei = fCoord->GetLastXBin(); i <= ei; ++i, ++ir) {
         for (Int_t j = fCoord->GetFirstYBin(), jr = 0, ej = fCoord->GetLastYBin(); j <= ej; ++j, ++jr) {
            Double_t val = 0.;
            for (Int_t level = low; level < up; ++ level)
               val += fHist->GetBinContent(i, j, level);
            fMinMax.second = TMath::Max(fMinMax.second, val);
            fMinMax.first = TMath::Min(fMinMax.first, val);
         }
      }
      break;
   }*/
}

////////////////////////////////////////////////////////////////////////////////
///Initialize color palette.

Bool_t TGLTH3Slice::PreparePalette()const
{
   UInt_t paletteSize = ((TH1 *)fHist)->GetContour();
   if (!paletteSize && !(paletteSize = gStyle->GetNumberContours()))
      paletteSize = 20;

   return fPalette.GeneratePalette(paletteSize, fMinMax);
}

////////////////////////////////////////////////////////////////////////////////
/// Prepare TexCoords.

void TGLTH3Slice::PrepareTexCoords(Double_t pos, Int_t low, Int_t up)const
{
   switch (fAxisType) {
   case kXOZ:
      fTexCoords.resize(fCoord->GetNXBins() * fCoord->GetNZBins());
      fTexCoords.SetRowLen(fCoord->GetNXBins());
      if (!fF3) {

         for (Int_t j = fCoord->GetFirstZBin(), jt = 0, ej = fCoord->GetLastZBin(); j <= ej; ++j, ++jt) {
            for (Int_t i = fCoord->GetFirstXBin(), it = 0, ei = fCoord->GetLastXBin(); i <= ei; ++i, ++it) {
               Double_t val = 0.;
               for (Int_t level = low; level < up; ++ level)
                  val += fHist->GetBinContent(i, level, j);
               fTexCoords[jt][it] = fPalette.GetTexCoord(val);
            }
         }
      } else {
         for (Int_t j = fCoord->GetFirstZBin(), jt = 0, ej = fCoord->GetLastZBin(); j <= ej; ++j, ++jt) {
            for (Int_t i = fCoord->GetFirstXBin(), it = 0, ei = fCoord->GetLastXBin(); i <= ei; ++i, ++it) {
               Double_t val = fF3->Eval(fHist->GetXaxis()->GetBinCenter(i), pos, fHist->GetZaxis()->GetBinCenter(j));
               if (val > fMinMax.second)
                  val = fMinMax.second;
               else if (val < fMinMax.first)
                  val = fMinMax.first;
               fTexCoords[jt][it] = fPalette.GetTexCoord(val);
            }
         }
      }
      break;
   case kYOZ:
      fTexCoords.resize(fCoord->GetNYBins() * fCoord->GetNZBins());
      fTexCoords.SetRowLen(fCoord->GetNYBins());
      if (!fF3) {
         for (Int_t j = fCoord->GetFirstZBin(), jt = 0, ej = fCoord->GetLastZBin(); j <= ej; ++j, ++jt) {
            for (Int_t i = fCoord->GetFirstYBin(), it = 0, ei = fCoord->GetLastYBin(); i <= ei; ++i, ++it) {
               Double_t val = 0.;
               for (Int_t level = low; level < up; ++ level)
                  val += fHist->GetBinContent(level, i, j);
               fTexCoords[jt][it] = fPalette.GetTexCoord(val);
            }
         }
      } else {
         for (Int_t j = fCoord->GetFirstZBin(), jt = 0, ej = fCoord->GetLastZBin(); j <= ej; ++j, ++jt) {
            for (Int_t i = fCoord->GetFirstXBin(), it = 0, ei = fCoord->GetLastXBin(); i <= ei; ++i, ++it) {
               Double_t val = fF3->Eval(pos, fHist->GetYaxis()->GetBinCenter(i), fHist->GetZaxis()->GetBinCenter(j));
               if (val > fMinMax.second)
                  val = fMinMax.second;
               else if (val < fMinMax.first)
                  val = fMinMax.first;
               fTexCoords[jt][it] = fPalette.GetTexCoord(val);
            }
         }
      }
      break;
   case kXOY:
      fTexCoords.resize(fCoord->GetNXBins() * fCoord->GetNYBins());
      fTexCoords.SetRowLen(fCoord->GetNYBins());
      if (!fF3) {
         for (Int_t i = fCoord->GetFirstXBin(), ir = 0, ei = fCoord->GetLastXBin(); i <= ei; ++i, ++ir) {
            for (Int_t j = fCoord->GetFirstYBin(), jr = 0, ej = fCoord->GetLastYBin(); j <= ej; ++j, ++jr) {
               Double_t val = 0.;
               for (Int_t level = low; level < up; ++ level)
                  val += fHist->GetBinContent(i, j, level);
               fTexCoords[ir][jr] = fPalette.GetTexCoord(val);
            }
         }
      } else {
         for (Int_t i = fCoord->GetFirstXBin(), it = 0, ei = fCoord->GetLastXBin(); i <= ei; ++i, ++it) {
            for (Int_t j = fCoord->GetFirstYBin(), jt = 0, ej = fCoord->GetLastYBin(); j <= ej; ++j, ++jt) {
               Double_t val = fF3->Eval(fHist->GetXaxis()->GetBinCenter(i), fHist->GetYaxis()->GetBinCenter(j), pos);
               if (val > fMinMax.second)
                  val = fMinMax.second;
               else if (val < fMinMax.first)
                  val = fMinMax.first;
               fTexCoords[it][jt] = fPalette.GetTexCoord(val);
            }
         }

      }
      break;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw slice textured.

void TGLTH3Slice::DrawSliceTextured(Double_t pos)const
{
   const Double_t xScale = fCoord->GetXScale();
   const Double_t yScale = fCoord->GetYScale();
   const Double_t zScale = fCoord->GetZScale();
   const TAxis *xA = fHist->GetXaxis();
   const TAxis *yA = fHist->GetYaxis();
   const TAxis *zA = fHist->GetZaxis();

   switch (fAxisType) {
   case kXOZ:
      pos *= yScale;
      for (Int_t j = fCoord->GetFirstZBin(), jt = 0, ej = fCoord->GetLastZBin(); j < ej; ++j, ++jt) {
         for (Int_t i = fCoord->GetFirstXBin(), it = 0, ei = fCoord->GetLastXBin(); i < ei; ++i, ++it) {
            const Double_t xMin = xA->GetBinCenter(i) * xScale;
            const Double_t xMax = xA->GetBinCenter(i + 1) * xScale;
            const Double_t zMin = zA->GetBinCenter(j) * zScale;
            const Double_t zMax = zA->GetBinCenter(j + 1) * zScale;
            glBegin(GL_POLYGON);
            glTexCoord1d(fTexCoords[jt][it]);
            glVertex3d(xMin, pos, zMin);
            glTexCoord1d(fTexCoords[jt + 1][it]);
            glVertex3d(xMin, pos, zMax);
            glTexCoord1d(fTexCoords[jt + 1][it + 1]);
            glVertex3d(xMax, pos, zMax);
            glTexCoord1d(fTexCoords[jt][it + 1]);
            glVertex3d(xMax, pos, zMin);
            glEnd();
         }
      }
      break;
   case kYOZ:
      pos *= xScale;
      for (Int_t j = fCoord->GetFirstZBin(), jt = 0, ej = fCoord->GetLastZBin(); j < ej; ++j, ++jt) {
         for (Int_t i = fCoord->GetFirstYBin(), it = 0, ei = fCoord->GetLastYBin(); i < ei; ++i, ++it) {
            const Double_t yMin = yA->GetBinCenter(i) * yScale;
            const Double_t yMax = yA->GetBinCenter(i + 1) * yScale;
            const Double_t zMin = zA->GetBinCenter(j) * zScale;
            const Double_t zMax = zA->GetBinCenter(j + 1) * zScale;
            glBegin(GL_POLYGON);
            glTexCoord1d(fTexCoords[jt][it]);
            glVertex3d(pos, yMin, zMin);
            glTexCoord1d(fTexCoords[jt][it + 1]);
            glVertex3d(pos, yMax, zMin);
            glTexCoord1d(fTexCoords[jt + 1][it + 1]);
            glVertex3d(pos, yMax, zMax);
            glTexCoord1d(fTexCoords[jt + 1][it]);
            glVertex3d(pos, yMin, zMax);
            glEnd();
         }
      }
      break;
   case kXOY:
      pos *= zScale;
      for (Int_t j = fCoord->GetFirstXBin(), jt = 0, ej = fCoord->GetLastXBin(); j < ej; ++j, ++jt) {
         for (Int_t i = fCoord->GetFirstYBin(), it = 0, ei = fCoord->GetLastYBin(); i < ei; ++i, ++it) {
            const Double_t xMin = xA->GetBinCenter(j) * xScale;
            const Double_t xMax = xA->GetBinCenter(j + 1) * xScale;
            const Double_t yMin = yA->GetBinCenter(i) * yScale;
            const Double_t yMax = yA->GetBinCenter(i + 1) * yScale;
            glBegin(GL_POLYGON);
            glTexCoord1d(fTexCoords[jt + 1][it]);
            glVertex3d(xMax, yMin, pos);
            glTexCoord1d(fTexCoords[jt + 1][it + 1]);
            glVertex3d(xMax, yMax, pos);
            glTexCoord1d(fTexCoords[jt][it + 1]);
            glVertex3d(xMin, yMax, pos);
            glTexCoord1d(fTexCoords[jt][it]);
            glVertex3d(xMin, yMin, pos);
            glEnd();
         }
      }
      break;
   }
}

namespace {

////////////////////////////////////////////////////////////////////////////////

void DrawBoxOutline(Double_t xMin, Double_t xMax, Double_t yMin,
                    Double_t yMax, Double_t zMin, Double_t zMax)
{
   glBegin(GL_LINE_LOOP);
   glVertex3d(xMin, yMin, zMin);
   glVertex3d(xMax, yMin, zMin);
   glVertex3d(xMax, yMax, zMin);
   glVertex3d(xMin, yMax, zMin);
   glEnd();

   glBegin(GL_LINE_LOOP);
   glVertex3d(xMin, yMin, zMax);
   glVertex3d(xMax, yMin, zMax);
   glVertex3d(xMax, yMax, zMax);
   glVertex3d(xMin, yMax, zMax);
   glEnd();

   glBegin(GL_LINES);
   glVertex3d(xMin, yMin, zMin);
   glVertex3d(xMin, yMin, zMax);
   glVertex3d(xMax, yMin, zMin);
   glVertex3d(xMax, yMin, zMax);
   glVertex3d(xMax, yMax, zMin);
   glVertex3d(xMax, yMax, zMax);
   glVertex3d(xMin, yMax, zMin);
   glVertex3d(xMin, yMax, zMax);
   glEnd();
}

}

////////////////////////////////////////////////////////////////////////////////
/// Draw slice frame.

void TGLTH3Slice::DrawSliceFrame(Int_t low, Int_t up)const
{
   glColor3d(1., 0., 0.);
   const TGLVertex3 *box = fBox->Get3DBox();

   switch (fAxisType) {
   case kXOZ:
      DrawBoxOutline(box[0].X(), box[1].X(),
                     fAxis->GetBinLowEdge(low) * fCoord->GetYScale(),
                     fAxis->GetBinUpEdge(up - 1) * fCoord->GetYScale(),
                     box[0].Z(), box[4].Z());
      break;
   case kYOZ:
      DrawBoxOutline(fAxis->GetBinLowEdge(low) * fCoord->GetXScale(),
                     fAxis->GetBinUpEdge(up - 1) * fCoord->GetXScale(),
                     box[0].Y(), box[2].Y(),
                     box[0].Z(), box[4].Z());
      break;
   case kXOY:
      DrawBoxOutline(box[0].X(), box[1].X(),
                     box[0].Y(), box[2].Y(),
                     fAxis->GetBinLowEdge(low) * fCoord->GetZScale(),
                     fAxis->GetBinUpEdge(up - 1) * fCoord->GetZScale());
      break;
   }
}

namespace Rgl {

////////////////////////////////////////////////////////////////////////////////

PlotTranslation::PlotTranslation(const TGLPlotPainter *painter)
                   : fPainter(painter)
{
   const TGLVertex3 *box = fPainter->fBackBox.Get3DBox();
   const Double_t center[] = {(box[0].X() + box[1].X()) / 2,
                              (box[0].Y() + box[2].Y()) / 2,
                              (box[0].Z() + box[4].Z()) / 2};

   fPainter->SaveModelviewMatrix();
   glTranslated(-center[0], -center[1], -center[2]);
}

////////////////////////////////////////////////////////////////////////////////

PlotTranslation::~PlotTranslation()
{
   fPainter->RestoreModelviewMatrix();
}

namespace
{

const Double_t lr = 0.85;
const Double_t rr = 0.9;

}

////////////////////////////////////////////////////////////////////////////////
///Draw. Palette.

void DrawPalette(const TGLPlotCamera * camera, const TGLLevelPalette & palette)
{
   const TGLDisableGuard light(GL_LIGHTING);
   const TGLDisableGuard depth(GL_DEPTH_TEST);
   const TGLEnableGuard blend(GL_BLEND);

   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   glOrtho(0, camera->GetWidth(), 0, camera->GetHeight(), -1., 1.);

   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();

   const Double_t leftX = camera->GetWidth() * lr;
   const Double_t rightX = camera->GetWidth() * rr;
   const Double_t margin = 0.1 * camera->GetHeight();
   const Double_t h = (camera->GetHeight() * 0.8) / palette.GetPaletteSize();

   for (Int_t i = 0, e = palette.GetPaletteSize(); i < e; ++i) {
      glBegin(GL_POLYGON);
      const UChar_t * color = palette.GetColour(i);
      glColor4ub(color[0], color[1], color[2], 150);
      glVertex2d(leftX, margin + i * h);
      glVertex2d(rightX, margin + i * h);
      glVertex2d(rightX, margin + (i + 1) * h);
      glVertex2d(leftX, margin + (i + 1) * h);
      glEnd();
   }

   const TGLEnableGuard  smoothGuard(GL_LINE_SMOOTH);
   glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
   glColor4d(0., 0., 0., 0.5);

   for (Int_t i = 0, e = palette.GetPaletteSize(); i < e; ++i) {
      glBegin(GL_LINE_LOOP);
      glVertex2d(leftX, margin + i * h);
      glVertex2d(rightX, margin + i * h);
      glVertex2d(rightX, margin + (i + 1) * h);
      glVertex2d(leftX, margin + (i + 1) * h);
      glEnd();
   }

}

////////////////////////////////////////////////////////////////////////////////
///Draw. Palette.

void DrawPalette(const TGLPlotCamera * camera, const TGLLevelPalette & palette,
                 const std::vector<Double_t> & levels)
{
   const TGLDisableGuard light(GL_LIGHTING);
   const TGLDisableGuard depth(GL_DEPTH_TEST);
   const TGLEnableGuard blend(GL_BLEND);

   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   glOrtho(0, camera->GetWidth(), 0, camera->GetHeight(), -1., 1.);

   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();

   const Double_t leftX = camera->GetWidth() * lr;
   const Double_t rightX = camera->GetWidth() * rr;
   const Double_t margin = 0.1 * camera->GetHeight();
   const Double_t h = (camera->GetHeight() * 0.8);
   const Double_t range = levels.back() - levels.front();

   const UChar_t opacity = 200;

   for (Int_t i = 0, e = palette.GetPaletteSize(); i < e; ++i) {
      const Double_t yMin = margin + (levels[i] - levels.front()) / range * h;
      const Double_t yMax = margin + (levels[i + 1] - levels.front()) / range * h;
      glBegin(GL_POLYGON);
      const UChar_t * color = palette.GetColour(i);
      glColor4ub(color[0], color[1], color[2], opacity);
      glVertex2d(leftX, yMin);
      glVertex2d(rightX, yMin);
      glVertex2d(rightX, yMax);
      glVertex2d(leftX, yMax);
      glEnd();
   }

   const TGLEnableGuard  smoothGuard(GL_LINE_SMOOTH);
   glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
   glColor4d(0., 0., 0., 0.5);

   for (Int_t i = 0, e = palette.GetPaletteSize(); i < e; ++i) {
      const Double_t yMin = (levels[i] - levels.front()) / range * h;
      const Double_t yMax = (levels[i + 1] - levels.front()) / range * h;

      glBegin(GL_LINE_LOOP);
      glVertex2d(leftX, margin + yMin);
      glVertex2d(rightX, margin + yMin);
      glVertex2d(rightX, margin + yMax);
      glVertex2d(leftX, margin + yMax);
      glEnd();
   }

}

////////////////////////////////////////////////////////////////////////////////

void DrawPaletteAxis(const TGLPlotCamera * camera, const Range_t & minMax, Bool_t logZ)
{
   UInt_t pixelW = camera->GetWidth();
   UInt_t pixelH = camera->GetHeight();

   TGLUtil::InitializeIfNeeded();
   const Float_t scale = TGLUtil::GetScreenScalingFactor();
   if (scale > 1.) {
      pixelW = UInt_t(pixelW / scale);
      pixelH = UInt_t(pixelH / scale);
   }

   const Double_t x = gPad->AbsPixeltoX(Int_t(gPad->GetXlowNDC() * gPad->GetWw() + rr * pixelW));
   const Double_t yMin = gPad->AbsPixeltoY(Int_t(gPad->GetWh() - (pixelH * 0.1
                                           + gPad->GetYlowNDC() * gPad->GetWh())));
   const Double_t yMax = gPad->AbsPixeltoY(Int_t(gPad->GetWh() - (pixelH * 0.9
                                           + gPad->GetYlowNDC() * gPad->GetWh())));

   Double_t zMin = minMax.first;
   Double_t zMax = minMax.second;

   if (logZ) {
      zMin = TMath::Power(10, zMin);
      zMax = TMath::Power(10, zMax);
   }

   //Now, some stupid magic, to force ROOT's painting machine work as I want, not as it wants.
   const Bool_t logX = gPad->GetLogx();
   gPad->SetLogx(kFALSE);
   const Bool_t logY = gPad->GetLogy();
   gPad->SetLogy(kFALSE);

   TGaxis axisPainter(x, yMin, x, yMax, zMin, zMax, 510, logZ ? "G" : "");
   axisPainter.Paint();

   gPad->SetLogx(logX);
   gPad->SetLogy(logY);
}

//Constant for TGLH2PolyPainter.
const Double_t gH2PolyScaleXY = 1.2;

}

