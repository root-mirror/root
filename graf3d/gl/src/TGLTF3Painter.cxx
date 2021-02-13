// @(#)root/gl:$Id$
// Author:  Timur Pocheptsov  31/08/2006

/*************************************************************************
 * Copyright (C) 1995-2006, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <typeinfo>

#include "KeySymbols.h"
#include "TVirtualX.h"
#include "Buttons.h"
#include "TString.h"
#include "TROOT.h"
#include "TColor.h"
#include "TMath.h"
#include "TH3.h"
#include "TF3.h"

#include "TGLMarchingCubes.h"
#include "TGLPlotCamera.h"
#include "TGLTF3Painter.h"
#include "TGLIncludes.h"

/** \class TGLTF3Painter
\ingroup opengl
Plot-painter for TF3 functions.
*/

ClassImp(TGLTF3Painter);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TGLTF3Painter::TGLTF3Painter(TF3 *fun, TH1 *hist, TGLPlotCamera *camera, TGLPlotCoordinates *coord)
                  : TGLPlotPainter(hist, camera, coord, kFALSE, kFALSE, kFALSE),
                    fStyle(kDefault),
                    fF3(fun),
                    fXOZSlice("XOZ", (TH3 *)hist, fun, coord, &fBackBox, TGLTH3Slice::kXOZ),
                    fYOZSlice("YOZ", (TH3 *)hist, fun, coord, &fBackBox, TGLTH3Slice::kYOZ),
                    fXOYSlice("XOY", (TH3 *)hist, fun, coord, &fBackBox, TGLTH3Slice::kXOY)
{
}

////////////////////////////////////////////////////////////////////////////////
///Coords for point on surface under cursor.

char *TGLTF3Painter::GetPlotInfo(Int_t /*px*/, Int_t /*py*/)
{
   static char mess[] = { "fun3" };
   return mess;
}

////////////////////////////////////////////////////////////////////////////////
///Create mesh.

Bool_t TGLTF3Painter::InitGeometry()
{
   fCoord->SetCoordType(kGLCartesian);

   if (!fCoord->SetRanges(fHist, kFALSE, kTRUE))
      return kFALSE;

   fBackBox.SetPlotBox(fCoord->GetXRangeScaled(), fCoord->GetYRangeScaled(), fCoord->GetZRangeScaled());
   if (fCamera) fCamera->SetViewVolume(fBackBox.Get3DBox());

   //Build mesh for TF3 surface
   fMesh.ClearMesh();

   Rgl::Mc::TMeshBuilder<TF3, Double_t> builder(kFALSE);//no averaged normals.
   //Set grid parameters.
   Rgl::Mc::TGridGeometry<Double_t> geom(fXAxis, fYAxis, fZAxis, fCoord->GetXScale(),
                                         fCoord->GetYScale(), fCoord->GetZScale(),
                                         Rgl::Mc::TGridGeometry<Double_t>::kBinEdge);

   builder.BuildMesh(fF3, geom, &fMesh, 0.);

   if (fCoord->Modified()) {
      fUpdateSelection = kTRUE;
      const TGLVertex3 &vertex = fBackBox.Get3DBox()[0];
      fXOZSectionPos = vertex.Y();
      fYOZSectionPos = vertex.X();
      fXOYSectionPos = vertex.Z();
      fCoord->ResetModified();
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
///User clicks right mouse button (in a pad).

void TGLTF3Painter::StartPan(Int_t px, Int_t py)
{
   fMousePosition.fX = px;
   fMousePosition.fY = fCamera->GetHeight() - py;
   fCamera->StartPan(px, py);
   fBoxCut.StartMovement(px, fCamera->GetHeight() - py);
}

////////////////////////////////////////////////////////////////////////////////
///User's moving mouse cursor, with middle mouse button pressed (for pad).
///Calculate 3d shift related to 2d mouse movement.
///Slicing is disabled (since somebody has broken it).

void TGLTF3Painter::Pan(Int_t px, Int_t py)
{
   if (fSelectedPart >= fSelectionBase) {//Pan camera.
      SaveModelviewMatrix();
      SaveProjectionMatrix();

      fCamera->SetCamera();
      fCamera->Apply(fPadPhi, fPadTheta);
      fCamera->Pan(px, py);

      RestoreProjectionMatrix();
      RestoreModelviewMatrix();
   } else if (fSelectedPart > 0) {
      //Convert py into bottom-top orientation.
      //Possibly, move box here
      py = fCamera->GetHeight() - py;

      SaveModelviewMatrix();
      SaveProjectionMatrix();

      fCamera->SetCamera();
      fCamera->Apply(fPadPhi, fPadTheta);

      if (!fHighColor) {
         if (fBoxCut.IsActive() && (fSelectedPart >= kXAxis && fSelectedPart <= kZAxis)) {
            fBoxCut.MoveBox(px, py, fSelectedPart);
         } else {
            //MoveSection(px, py);
         }
      } else {
         //MoveSection(px, py);
      }

      RestoreProjectionMatrix();
      RestoreModelviewMatrix();
   }

   fMousePosition.fX = px, fMousePosition.fY = py;
   fUpdateSelection = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
///No options for tf3

void TGLTF3Painter::AddOption(const TString &/*option*/)
{
}

////////////////////////////////////////////////////////////////////////////////
///Change color scheme.

void TGLTF3Painter::ProcessEvent(Int_t event, Int_t /*px*/, Int_t py)
{
   if (event == kKeyPress) {
      if (py == kKey_s || py == kKey_S) {
         fStyle < kMaple2 ? fStyle = ETF3Style(fStyle + 1) : fStyle = kDefault;
      } else if (py == kKey_c || py == kKey_C) {
         if (fHighColor)
            Info("ProcessEvent", "Cut box does not work in high color, please, switch to true color");
         else {
            fBoxCut.TurnOnOff();
            fUpdateSelection = kTRUE;
         }
      }
   } else if (event == kButton1Double && (fBoxCut.IsActive() || HasSections())) {
      if (fBoxCut.IsActive())
         fBoxCut.TurnOnOff();
      const TGLVertex3 *frame = fBackBox.Get3DBox();
      fXOZSectionPos = frame[0].Y();
      fYOZSectionPos = frame[0].X();
      fXOYSectionPos = frame[0].Z();

      if (!gVirtualX->IsCmdThread())
         gROOT->ProcessLineFast(Form("((TGLPlotPainter *)0x%lx)->Paint()", (ULong_t)this));
      else
         Paint();
   }
}

////////////////////////////////////////////////////////////////////////////////
///Initialize OpenGL state variables.

void TGLTF3Painter::InitGL() const
{
   glEnable(GL_LIGHTING);
   glEnable(GL_LIGHT0);
   glEnable(GL_DEPTH_TEST);
   glDisable(GL_CULL_FACE);
   glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
}

////////////////////////////////////////////////////////////////////////////////
///Initialize OpenGL state variables.

void TGLTF3Painter::DeInitGL() const
{
   glDisable(GL_LIGHTING);
   glDisable(GL_LIGHT0);
   glDisable(GL_DEPTH_TEST);
   glDisable(GL_CULL_FACE);
   glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
}

////////////////////////////////////////////////////////////////////////////////
///Draw triangles, no normals, no lighting.

void TGLTF3Painter::DrawToSelectionBuffer() const
{
   Rgl::ObjectIDToColor(fSelectionBase, fHighColor);

   if (!fBoxCut.IsActive())
      Rgl::DrawMesh(fMesh.fVerts, fMesh.fTris);
   else
      Rgl::DrawMesh(fMesh.fVerts, fMesh.fTris, fBoxCut);
}

////////////////////////////////////////////////////////////////////////////////
///Surface with material properties and lighting.

void TGLTF3Painter::DrawDefaultPlot() const
{
   if (HasSections()) {
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      glDepthMask(GL_FALSE);
   }

   SetSurfaceColor();

   if (!fBoxCut.IsActive()) {
      Rgl::DrawMesh(fMesh.fVerts, fMesh.fNorms, fMesh.fTris);
   } else {
      Rgl::DrawMesh(fMesh.fVerts, fMesh.fNorms, fMesh.fTris, fBoxCut);
   }

   if (HasSections()) {
      glDisable(GL_BLEND);
      glDepthMask(GL_TRUE);
   }
}

////////////////////////////////////////////////////////////////////////////////
///Colored surface, without lighting and
///material properties.

void TGLTF3Painter::DrawMaplePlot() const
{
   const TGLDisableGuard lightGuard(GL_LIGHTING);

   if (HasSections() && fStyle < kMaple2) {
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      glDepthMask(GL_FALSE);
   }

   if (fStyle == kMaple1) {//Shaded polygons and outlines.
      glEnable(GL_POLYGON_OFFSET_FILL);//[1
      glPolygonOffset(1.f, 1.f);
   } else if (fStyle == kMaple2)//Colored outlines only.
      glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);//[2

   if(!fBoxCut.IsActive())
      Rgl::DrawMapleMesh(fMesh.fVerts, fMesh.fNorms, fMesh.fTris);
   else
      Rgl::DrawMapleMesh(fMesh.fVerts, fMesh.fNorms, fMesh.fTris, fBoxCut);

   if (fStyle == kMaple1) {
      //Draw outlines.
      glDisable(GL_POLYGON_OFFSET_FILL);//1]
      glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);//[3
      glColor4d(0., 0., 0., 0.25);

      if(!fBoxCut.IsActive())
         Rgl::DrawMesh(fMesh.fVerts, fMesh.fTris);
      else
         Rgl::DrawMesh(fMesh.fVerts, fMesh.fTris, fBoxCut);

      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);//[3
   } else if (fStyle == kMaple2)
      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

   if (HasSections() && fStyle < kMaple2) {
      glDisable(GL_BLEND);
      glDepthMask(GL_TRUE);
   }
}

////////////////////////////////////////////////////////////////////////////////
///Draw mesh.

void TGLTF3Painter::DrawPlot() const
{
   //Shift plot to point of origin.
   const Rgl::PlotTranslation trGuard(this);

   fBackBox.DrawBox(fSelectedPart, fSelectionPass, fZLevels, fHighColor);
   DrawSections();

   if (fSelectionPass) {
      DrawToSelectionBuffer();
   } else if (fStyle == kDefault) {
      DrawDefaultPlot();
   } else {
      DrawMaplePlot();
   }

   if (fBoxCut.IsActive())
      fBoxCut.DrawBox(fSelectionPass, fSelectedPart);
}

////////////////////////////////////////////////////////////////////////////////
///Set color for surface.

void TGLTF3Painter::SetSurfaceColor() const
{
   Float_t diffColor[] = {0.8f, 0.8f, 0.8f, 0.15f};

   if (fF3->GetFillColor() != kWhite)
      if (const TColor *c = gROOT->GetColor(fF3->GetFillColor()))
         c->GetRGB(diffColor[0], diffColor[1], diffColor[2]);

   glMaterialfv(GL_BACK, GL_DIFFUSE, diffColor);
   diffColor[0] /= 2, diffColor[1] /= 2, diffColor[2] /= 2;
   glMaterialfv(GL_FRONT, GL_DIFFUSE, diffColor);
   const Float_t specColor[] = {1.f, 1.f, 1.f, 1.f};
   glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specColor);
   glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 70.f);
}

////////////////////////////////////////////////////////////////////////////////
///Any section exists.

Bool_t TGLTF3Painter::HasSections() const
{
   return fXOZSectionPos > fBackBox.Get3DBox()[0].Y() ||
          fYOZSectionPos > fBackBox.Get3DBox()[0].X() ||
          fXOYSectionPos > fBackBox.Get3DBox()[0].Z();
}

////////////////////////////////////////////////////////////////////////////////
/// Draw XOZ parallel section.

void TGLTF3Painter::DrawSectionXOZ() const
{
   if (fSelectionPass)
      return;
   fXOZSlice.DrawSlice(fXOZSectionPos / fCoord->GetYScale());
}

////////////////////////////////////////////////////////////////////////////////
/// Draw YOZ parallel section.

void TGLTF3Painter::DrawSectionYOZ() const
{
   if (fSelectionPass)
      return;
   fYOZSlice.DrawSlice(fYOZSectionPos / fCoord->GetXScale());
}

////////////////////////////////////////////////////////////////////////////////
/// Draw XOY parallel section.

void TGLTF3Painter::DrawSectionXOY() const
{
   if (fSelectionPass)
      return;
   fXOYSlice.DrawSlice(fXOYSectionPos / fCoord->GetZScale());
}


/** \class TGLIsoPainter
\ingroup opengl
"gliso" option for TH3.
*/

ClassImp(TGLIsoPainter);

////////////////////////////////////////////////////////////////////////////////
///Constructor.

TGLIsoPainter::TGLIsoPainter(TH1 *hist, TGLPlotCamera *camera, TGLPlotCoordinates *coord)
                  : TGLPlotPainter(hist, camera, coord, kFALSE, kFALSE, kFALSE),
                    fXOZSlice("XOZ", (TH3 *)hist, coord, &fBackBox, TGLTH3Slice::kXOZ),
                    fYOZSlice("YOZ", (TH3 *)hist, coord, &fBackBox, TGLTH3Slice::kYOZ),
                    fXOYSlice("XOY", (TH3 *)hist, coord, &fBackBox, TGLTH3Slice::kXOY),
                    fInit(kFALSE)
{
   if (hist->GetDimension() < 3)
      Error("TGLIsoPainter::TGLIsoPainter", "Wrong type of histogramm, must have 3 dimensions");
}

////////////////////////////////////////////////////////////////////////////////
///Return info for plot part under cursor.

char *TGLIsoPainter::GetPlotInfo(Int_t /*px*/, Int_t /*py*/)
{
   static char mess[] = { "iso" };
   return mess;
}

////////////////////////////////////////////////////////////////////////////////
///Initializes meshes for 3d iso contours.

Bool_t TGLIsoPainter::InitGeometry()
{
   if (fHist->GetDimension() < 3) {
      Error("TGLIsoPainter::TGLIsoPainter", "Wrong type of histogramm, must have 3 dimensions");
      return kFALSE;
   }

   //Create mesh.
   if (fInit)
      return kTRUE;

   //Only in cartesian.
   fCoord->SetCoordType(kGLCartesian);
   if (!fCoord->SetRanges(fHist, kFALSE, kTRUE))
      return kFALSE;

   fBackBox.SetPlotBox(fCoord->GetXRangeScaled(), fCoord->GetYRangeScaled(), fCoord->GetZRangeScaled());
   if (fCamera) fCamera->SetViewVolume(fBackBox.Get3DBox());

   //Move old meshes into the cache.
   if (!fIsos.empty())
      fCache.splice(fCache.begin(), fIsos);
   //Number of contours == number of iso surfaces.
   UInt_t nContours = fHist->GetContour();

   if (nContours > 1) {
      fColorLevels.resize(nContours);
      FindMinMax();

      if (fHist->TestBit(TH1::kUserContour)) {
         //There are user defined contours (iso-levels).
         for (UInt_t i = 0; i < nContours; ++i)
            fColorLevels[i] = fHist->GetContourLevelPad(i);
      } else {
         //Equidistant iso-surfaces.
         const Double_t isoStep = (fMinMax.second - fMinMax.first) / nContours;
         for (UInt_t i = 0; i < nContours; ++i)
            fColorLevels[i] = fMinMax.first + i * isoStep;
      }

      fPalette.GeneratePalette(nContours, fMinMax, kFALSE);
   } else {
      //Only one iso (ROOT's standard).
      fColorLevels.resize(nContours = 1);
      fColorLevels[0] = fHist->GetSumOfWeights() / (fHist->GetNbinsX() * fHist->GetNbinsY() * fHist->GetNbinsZ());
   }

   MeshIter_t firstMesh = fCache.begin();
   //Initialize meshes, trying to reuse mesh from
   //mesh cache.
   for (UInt_t i = 0; i < nContours; ++i) {
      if (firstMesh != fCache.end()) {
         //There is a mesh in a cache.
         SetMesh(*firstMesh, fColorLevels[i]);
         MeshIter_t next = firstMesh;
         ++next;
         fIsos.splice(fIsos.begin(), fCache, firstMesh);
         firstMesh = next;
      } else {
         //No meshes in a cache.
         //Create new one and _swap_ data (look at Mesh_t::Swap in a header)
         //between empty mesh in a list and this mesh
         //to avoid real copying.
         Mesh_t newMesh;
         SetMesh(newMesh, fColorLevels[i]);
         fIsos.push_back(fDummyMesh);
         fIsos.back().Swap(newMesh);
      }
   }

   if (fCoord->Modified()) {
      fUpdateSelection = kTRUE;
      fXOZSectionPos = fBackBox.Get3DBox()[0].Y();
      fYOZSectionPos = fBackBox.Get3DBox()[0].X();
      fXOYSectionPos = fBackBox.Get3DBox()[0].Z();
      fCoord->ResetModified();
   }

   //Avoid rebuilding the mesh.
   fInit = kTRUE;

   return kTRUE;

}

////////////////////////////////////////////////////////////////////////////////
///User clicks right mouse button (in a pad).

void TGLIsoPainter::StartPan(Int_t px, Int_t py)
{
   fMousePosition.fX = px;
   fMousePosition.fY = fCamera->GetHeight() - py;
   fCamera->StartPan(px, py);
   fBoxCut.StartMovement(px, fCamera->GetHeight() - py);
}

////////////////////////////////////////////////////////////////////////////////
///User's moving mouse cursor, with middle mouse button pressed (for pad).
///Calculate 3d shift related to 2d mouse movement.
///User's moving mouse cursor, with middle mouse button pressed (for pad).
///Calculate 3d shift related to 2d mouse movement.

void TGLIsoPainter::Pan(Int_t px, Int_t py)
{
   if (fSelectedPart >= fSelectionBase) {//Pan camera.
      SaveModelviewMatrix();
      SaveProjectionMatrix();

      fCamera->SetCamera();
      fCamera->Apply(fPadPhi, fPadTheta);
      fCamera->Pan(px, py);

      RestoreProjectionMatrix();
      RestoreModelviewMatrix();
   } else if (fSelectedPart > 0) {
      //Convert py into bottom-top orientation.
      //Possibly, move box here
      py = fCamera->GetHeight() - py;

      SaveModelviewMatrix();
      SaveProjectionMatrix();

      fCamera->SetCamera();
      fCamera->Apply(fPadPhi, fPadTheta);

      if (!fHighColor) {
         if (fBoxCut.IsActive() && (fSelectedPart >= kXAxis && fSelectedPart <= kZAxis)) {
            fBoxCut.MoveBox(px, py, fSelectedPart);
         } else {
            //MoveSection(px, py);
         }
      } else {
         //MoveSection(px, py);
      }

      RestoreProjectionMatrix();
      RestoreModelviewMatrix();

   }

   fMousePosition.fX = px, fMousePosition.fY = py;
   fUpdateSelection = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
///No additional options for TGLIsoPainter.

void TGLIsoPainter::AddOption(const TString &/*option*/)
{
}

////////////////////////////////////////////////////////////////////////////////
///Change color scheme.

void TGLIsoPainter::ProcessEvent(Int_t event, Int_t /*px*/, Int_t py)
{
   if (event == kKeyPress) {
      if (py == kKey_c || py == kKey_C) {
         if (fHighColor)
            Info("ProcessEvent", "Cut box does not work in high color, please, switch to true color");
         else {
            fBoxCut.TurnOnOff();
            fUpdateSelection = kTRUE;
         }
      }
   } else if (event == kButton1Double && (fBoxCut.IsActive() || HasSections())) {
      if (fBoxCut.IsActive())
         fBoxCut.TurnOnOff();
      const TGLVertex3 *frame = fBackBox.Get3DBox();
      fXOZSectionPos = frame[0].Y();
      fYOZSectionPos = frame[0].X();
      fXOYSectionPos = frame[0].Z();

      if (!gVirtualX->IsCmdThread())
         gROOT->ProcessLineFast(Form("((TGLPlotPainter *)0x%lx)->Paint()", (ULong_t)this));
      else
         Paint();
   }
}

////////////////////////////////////////////////////////////////////////////////
///Initialize OpenGL state variables.

void TGLIsoPainter::InitGL() const
{
   glEnable(GL_LIGHTING);
   glEnable(GL_LIGHT0);
   glEnable(GL_DEPTH_TEST);
   glDisable(GL_CULL_FACE);
   glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
}

////////////////////////////////////////////////////////////////////////////////
///Initialize OpenGL state variables.

void TGLIsoPainter::DeInitGL() const
{
   glDisable(GL_LIGHTING);
   glDisable(GL_LIGHT0);
   glDisable(GL_DEPTH_TEST);
   glDisable(GL_CULL_FACE);
   glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
}

////////////////////////////////////////////////////////////////////////////////
///Draw mesh.

void TGLIsoPainter::DrawPlot() const
{
   //Shift plot to point of origin.
   const Rgl::PlotTranslation trGuard(this);


   fBackBox.DrawBox(fSelectedPart, fSelectionPass, fZLevels, fHighColor);
   DrawSections();

   if (fIsos.size() != fColorLevels.size()) {
      Error("TGLIsoPainter::DrawPlot", "Non-equal number of levels and isos");
      return;
   }

   if (!fSelectionPass && HasSections()) {
      //Surface is semi-transparent during dynamic profiling.
      //Having several complex nested surfaces, it's not easy
      //(possible?) to implement correct and _efficient_ transparency
      //drawing. So, artefacts are possbile.
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      glDepthMask(GL_FALSE);
   }

   UInt_t colorInd = 0;
   ConstMeshIter_t iso = fIsos.begin();

   for (; iso != fIsos.end(); ++iso, ++colorInd)
      DrawMesh(*iso, colorInd);

   if (!fSelectionPass && HasSections()) {
      glDisable(GL_BLEND);
      glDepthMask(GL_TRUE);
   }

   if (fBoxCut.IsActive())
      fBoxCut.DrawBox(fSelectionPass, fSelectedPart);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw XOZ parallel section.

void TGLIsoPainter::DrawSectionXOZ() const
{
   if (fSelectionPass)
      return;
   fXOZSlice.DrawSlice(fXOZSectionPos / fCoord->GetYScale());
}

////////////////////////////////////////////////////////////////////////////////
/// Draw YOZ parallel section.

void TGLIsoPainter::DrawSectionYOZ() const
{
   if (fSelectionPass)
      return;
   fYOZSlice.DrawSlice(fYOZSectionPos / fCoord->GetXScale());
}

////////////////////////////////////////////////////////////////////////////////
/// Draw XOY parallel section.

void TGLIsoPainter::DrawSectionXOY() const
{
   if (fSelectionPass)
      return;
   fXOYSlice.DrawSlice(fXOYSectionPos / fCoord->GetZScale());
}

////////////////////////////////////////////////////////////////////////////////
///Any section exists.

Bool_t TGLIsoPainter::HasSections() const
{
   return fXOZSectionPos > fBackBox.Get3DBox()[0].Y() || fYOZSectionPos > fBackBox.Get3DBox()[0].X() ||
          fXOYSectionPos > fBackBox.Get3DBox()[0].Z();
}

////////////////////////////////////////////////////////////////////////////////
///Set color for surface.

void TGLIsoPainter::SetSurfaceColor(Int_t ind) const
{
   Float_t diffColor[] = {0.8f, 0.8f, 0.8f, 0.25f};

   if (fColorLevels.size() == 1) {
      if (fHist->GetFillColor() != kWhite)
         if (const TColor *c = gROOT->GetColor(fHist->GetFillColor()))
            c->GetRGB(diffColor[0], diffColor[1], diffColor[2]);
   } else {
      const UChar_t *color = fPalette.GetColour(ind);
      diffColor[0] = color[0] / 255.;
      diffColor[1] = color[1] / 255.;
      diffColor[2] = color[2] / 255.;
   }

   glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffColor);
   const Float_t specColor[] = {1.f, 1.f, 1.f, 1.f};
   glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specColor);
   diffColor[0] /= 3.5, diffColor[1] /= 3.5, diffColor[2] /= 3.5;
   glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, diffColor);
   glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 30.f);
}

////////////////////////////////////////////////////////////////////////////////
///Grid geometry.

void TGLIsoPainter::SetMesh(Mesh_t &m, Double_t isoValue)
{
   Rgl::Mc::TGridGeometry<Float_t> geom(fXAxis, fYAxis, fZAxis, fCoord->GetXScale(),
                                        fCoord->GetYScale(), fCoord->GetZScale());
   //Clear mesh if it was from cache.
   m.ClearMesh();
   //Select correct TMeshBuilder type.
   if (typeid(*fHist) == typeid(TH3C)) {
      Rgl::Mc::TMeshBuilder<TH3C, Float_t> builder(kTRUE);
      builder.BuildMesh(static_cast<TH3C *>(fHist), geom, &m, isoValue);
   } else if (typeid(*fHist) == typeid(TH3S)) {
      Rgl::Mc::TMeshBuilder<TH3S, Float_t> builder(kTRUE);
      builder.BuildMesh(static_cast<TH3S *>(fHist), geom, &m, isoValue);
   } else if (typeid(*fHist) == typeid(TH3I)) {
      Rgl::Mc::TMeshBuilder<TH3I, Float_t> builder(kTRUE);
      builder.BuildMesh(static_cast<TH3I *>(fHist), geom, &m, isoValue);
   } else if (typeid(*fHist) == typeid(TH3F)) {
      Rgl::Mc::TMeshBuilder<TH3F, Float_t> builder(kTRUE);
      builder.BuildMesh(static_cast<TH3F *>(fHist), geom, &m, isoValue);
   } else if (typeid(*fHist) == typeid(TH3D)) {
      Rgl::Mc::TMeshBuilder<TH3D, Float_t> builder(kTRUE);
      builder.BuildMesh(static_cast<TH3D *>(fHist), geom, &m, isoValue);
   }
}

////////////////////////////////////////////////////////////////////////////////
///Draw TF3 surface

void TGLIsoPainter::DrawMesh(const Mesh_t &m, Int_t level) const
{
   if (!fSelectionPass)
      SetSurfaceColor(level);

   if (!fBoxCut.IsActive()) {
      if (!fSelectionPass)
         Rgl::DrawMesh(m.fVerts, m.fNorms, m.fTris);
      else {
         Rgl::ObjectIDToColor(fSelectionBase, fHighColor);
         Rgl::DrawMesh(m.fVerts, m.fTris);
      }
   } else {
      if (!fSelectionPass)
         Rgl::DrawMesh(m.fVerts, m.fNorms, m.fTris, fBoxCut);
      else {
         Rgl::ObjectIDToColor(fSelectionBase, fHighColor);
         Rgl::DrawMesh(m.fVerts, m.fTris, fBoxCut);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
///Find max/min bin contents for TH3.

void TGLIsoPainter::FindMinMax()
{
   fMinMax.first  = fHist->GetBinContent(fXAxis->GetFirst(), fYAxis->GetFirst(), fZAxis->GetFirst());
   fMinMax.second = fMinMax.first;

   for (Int_t i = fXAxis->GetFirst(), ei = fXAxis->GetLast(); i <= ei; ++i) {
      for (Int_t j = fYAxis->GetFirst(), ej = fYAxis->GetLast(); j <= ej; ++j) {
         for (Int_t k = fZAxis->GetFirst(), ek = fZAxis->GetLast(); k <= ek; ++k) {
            const Double_t binContent = fHist->GetBinContent(i, j, k);
            fMinMax.first  = TMath::Min(binContent, fMinMax.first);
            fMinMax.second = TMath::Max(binContent, fMinMax.second);
         }
      }
   }
}
