// @(#)root/gl:$Id$
// Author:  Timur Pocheptsov  28/07/2009

/*************************************************************************
 * Copyright (C) 1995-2009, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "KeySymbols.h"
#include "TVirtualX.h"
#include "Buttons.h"
#include "TString.h"
#include "TError.h"
#include "TROOT.h"
#include "TMath.h"

#include "TGLPlotCamera.h"
#include "TGL5DPainter.h"
#include "TGLPadUtils.h"
#include "TGLIncludes.h"
#include "TGL5D.h"

/** \class TGL5DPainter
\ingroup opengl
TGL5DPainter implements "gl5d" option for TTree::Draw.
Data (4D) is visualized as a set of iso-surfaces. 5D.
*/

////////////////////////////////////////////////////////////////////////////////
///Constructor.

TGL5DPainter::TGL5DPainter(TGL5DDataSet *data, TGLPlotCamera *camera, TGLPlotCoordinates *coord)
         : TGLPlotPainter(data, camera, coord),
           fMeshBuilder(kTRUE),//kTRUE == average normals.
           fInit(kFALSE),
           fData(data),
           fShowSlider(kFALSE),
           fAlpha(0.4),
           fNContours(kNContours)
{
   if (fData->fV4IsString)
      fNContours = Int_t(fData->fV4MinMax.second) - Int_t(fData->fV4MinMax.first) + 1;
}

////////////////////////////////////////////////////////////////////////////////
///Try to add new iso-surface.
///If something goes wrong, return
///iterator to the end of fIsos.

TGL5DPainter::SurfIter_t TGL5DPainter::AddSurface(Double_t v4, Color_t ci,
                                                  Double_t iso, Double_t sigma,
                                                  Double_t range, Int_t lownps)
{
   fData->SelectPoints(v4, range);

   if (fData->SelectedSize() < size_type(lownps)) {
      Warning("TGL5DPainter::AddSurface", "Too little points: %d", Int_t(fData->SelectedSize()));
      return fIsos.end();//This is a valid iterator, but an invalid surface.
   } else {
      Info("TGL5DPainter::AddSurface", "Selected %d points", Int_t(fData->SelectedSize()));
   }

   fKDE.BuildModel(fData, sigma);//Prepare density estimator.

   Info("TGL5DPainter::AddSurface", "Building the mesh ...");
   //Prepare grid parameters.
   Rgl::Mc::TGridGeometry<Float_t> geom(fXAxis, fYAxis, fZAxis,
                                        fCoord->GetXScale(),
                                        fCoord->GetYScale(),
                                        fCoord->GetZScale());
   Mesh_t mesh;
   fMeshBuilder.SetGeometry(fData);
   //Build a new mesh.
   fMeshBuilder.BuildMesh(&fKDE, geom, &mesh, iso);

   Info("TGL5DPainter::AddSurface", "Mesh has %d vertices", Int_t(mesh.fVerts.size() / 3));

   if (!mesh.fVerts.size())//I do not need an empty mesh.
      return fIsos.end();
   //Add surface with empty mesh and swap meshes.
   fIsos.push_front(fDummy);

   fIsos.front().fMesh.Swap(mesh);
   fIsos.front().f4D = v4;
   fIsos.front().fRange = range;
   fIsos.front().fShowCloud = kFALSE;
   fIsos.front().fHide = kFALSE;
   fIsos.front().fColor = ci;

   //Predictions for the 5-th variable.
   //Not-implemented yet.
   return fIsos.begin();
}

////////////////////////////////////////////////////////////////////////////////
///Add new surface. Simplified version for ged.

void TGL5DPainter::AddSurface(Double_t v4)
{
   const Rgl::Range_t &v4R = fData->fV4MinMax;
   const Bool_t isString   = fData->fV4IsString;
   const Double_t rms  = TMath::RMS(fData->fNP, fData->fV4);  //RMS of the N points.
   const Double_t d    = isString ? (v4R.second - v4R.first) / (fNContours - 1)
                                  : 6 * rms / fNContours;
   //alpha is in [0.1, 0.5], 1e-3 -s good for strings.
   const Double_t range = isString ? 1e-3 : fAlpha * d;

   AddSurface(v4, 1, 0.125, 0.05, range);
}

////////////////////////////////////////////////////////////////////////////////
///Remove iso-surface.

void TGL5DPainter::RemoveSurface(SurfIter_t surf)
{
   if (surf == fIsos.end()) {
      Error("TGL5DPainter::RemoveSurface", "Invalid iterator, surface does not exist.");
      return;
   }

   fIsos.erase(surf);
}

////////////////////////////////////////////////////////////////////////////////
///Return info for plot part under cursor.

char *TGL5DPainter::GetPlotInfo(Int_t /*px*/, Int_t /*py*/)
{
   static char mess[] = {"gl5d"};
   return mess;
}

////////////////////////////////////////////////////////////////////////////////
///Create mesh.
///InitGeometry creates surfaces for auto-iso levels.
///Called the first time and each time number of auto-levels is
///reset via the editor.

Bool_t TGL5DPainter::InitGeometry()
{
   if (fInit)
      return kTRUE;
   //Only in cartesian.
   fCoord->SetCoordType(kGLCartesian);

   if (!fCoord->SetRanges(fXAxis, fYAxis, fZAxis))
      return kFALSE;

   fIsos.clear();

   fBackBox.SetPlotBox(fCoord->GetXRangeScaled(),
                       fCoord->GetYRangeScaled(),
                       fCoord->GetZRangeScaled());
   if (fCamera)
      fCamera->SetViewVolume(fBackBox.Get3DBox());

   const Rgl::Range_t &v4R = fData->fV4MinMax;
   const Bool_t isString   = fData->fV4IsString;

   //Rene's code to automatically find iso-levels.
   const Double_t mean = TMath::Mean(fData->fNP, fData->fV4); //mean value of the NP points.
   const Double_t rms  = TMath::RMS(fData->fNP, fData->fV4);  //RMS of the N points.
   const Double_t min  = isString ? v4R.first : mean - 3 * rms; //take a range +- 3*xrms
   const Double_t d    = isString ? (v4R.second - v4R.first) / (fNContours - 1)
                                  : 6 * rms / fNContours;
   //alpha is in [0.1, 0.5], 1e-3 -s good for strings.
   const Double_t range = isString ? 1e-3 : fAlpha * d;

   Info("InitGeometry", "min = %g, mean = %g, rms = %g, dx = %g", min, mean, rms, d);

   for (Int_t j = 0; j < fNContours; ++j) {
      const Double_t isoLevel = min + j * d;
      Info("TGL5DPainter::InitGeometry", "Iso-level %g, range is %g ...", isoLevel, range);
      const Color_t color = j * 6 + 1;
      AddSurface(isoLevel, color, 0.125, 0.05, range);
   }

   if (fIsos.size())
      fBoxCut.TurnOnOff();

   return fInit = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
///User clicks right mouse button (in a pad).

void TGL5DPainter::StartPan(Int_t px, Int_t py)
{
   fMousePosition.fX = px;
   fMousePosition.fY = fCamera->GetHeight() - py;
   fCamera->StartPan(px, py);
   fBoxCut.StartMovement(px, fCamera->GetHeight() - py);
}

////////////////////////////////////////////////////////////////////////////////
///Mouse events handler.

void TGL5DPainter::Pan(Int_t px, Int_t py)
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
      py = fCamera->GetHeight() - py;

      SaveModelviewMatrix();
      SaveProjectionMatrix();

      fCamera->SetCamera();
      fCamera->Apply(fPadPhi, fPadTheta);

      if (!fHighColor) {
         if (fBoxCut.IsActive() && (fSelectedPart >= kXAxis && fSelectedPart <= kZAxis)) {
            fBoxCut.MoveBox(px, py, fSelectedPart);
         }
      }

      RestoreProjectionMatrix();
      RestoreModelviewMatrix();
   }

   fMousePosition.fX = px, fMousePosition.fY = py;
   fUpdateSelection = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
///No additional options for TGL5DPainter.

void TGL5DPainter::AddOption(const TString &/*option*/)
{
}

////////////////////////////////////////////////////////////////////////////////

void TGL5DPainter::ProcessEvent(Int_t event, Int_t /*px*/, Int_t py)
{
   //Change color scheme.
   if (event == kKeyPress) {
      if (py == kKey_c || py == kKey_C) {
         if (fHighColor)
            Info("ProcessEvent", "Cut box does not work in high color, please, switch to true color");
         else {
            fBoxCut.TurnOnOff();
            fUpdateSelection = kTRUE;
         }
      }
   } else if (event == kButton1Double && fBoxCut.IsActive()) {
      if (fBoxCut.IsActive())
         fBoxCut.TurnOnOff();
      if (!gVirtualX->IsCmdThread())
         gROOT->ProcessLineFast(Form("((TGLPlotPainter *)0x%lx)->Paint()", (ULong_t)this));
      else
         Paint();
   }
}

////////////////////////////////////////////////////////////////////////////////
///Set selection range parameter.

void TGL5DPainter::SetAlpha(Double_t newVal)
{
   if (fAlpha != newVal && !fData->fV4IsString) {
      fAlpha = newVal;
      fInit = kFALSE;
      InitGeometry();
   }

   if (fData->fV4IsString)
      Warning("SetAlpha", "Alpha is not required for string data (your 4-th dimension is string).");
}

////////////////////////////////////////////////////////////////////////////////
///Set the number of predefined contours.

void TGL5DPainter::SetNContours(Int_t n)
{
   if (n <= 0) {
      Warning("SetNContours", "Bad number of contours: %d", n);
      return;
   }

   fNContours = n;
   fInit = kFALSE;
   InitGeometry();
}

////////////////////////////////////////////////////////////////////////////////
///No need to create or delete meshes,
///number of meshes (iso-levels) are
///the same, but meshes must be rebuilt
///in new ranges.
///Only in cartesian.

void TGL5DPainter::ResetGeometryRanges()
{
   fCoord->SetRanges(fXAxis, fYAxis, fZAxis);
   fBackBox.SetPlotBox(fCoord->GetXRangeScaled(),
                       fCoord->GetYRangeScaled(),
                       fCoord->GetZRangeScaled());
   if (fCamera)
      fCamera->SetViewVolume(fBackBox.Get3DBox());
   //Iterate through all surfaces and re-calculate them.
   for (SurfIter_t surf = fIsos.begin(); surf != fIsos.end(); ++surf) {
      fData->SelectPoints(surf->f4D, surf->fRange);
      fKDE.BuildModel(fData, 0.05);//0.05 is sigma, will be controlled via GUI.
      Info("TGL5DPainter::ResetGeometryRanges", "Building the mesh ...");
      //Prepare grid parameters.
      Rgl::Mc::TGridGeometry<Float_t> geom(fXAxis, fYAxis, fZAxis,
                                           fCoord->GetXScale(),
                                           fCoord->GetYScale(),
                                           fCoord->GetZScale());
      fMeshBuilder.SetGeometry(fData);
      Mesh_t &mesh = surf->fMesh;
      //Clear old data.
      mesh.fVerts.clear();
      mesh.fNorms.clear();
      mesh.fTris.clear();
      //Build new mesh.
      fMeshBuilder.BuildMesh(&fKDE, geom, &mesh, 0.125);//0.125 will be set via GUI.
      Info("TGL5DPainter::AddSurface", "Mesh has %d vertices", Int_t(mesh.fVerts.size() / 3));
   }

   fBoxCut.ResetBoxGeometry();
}

////////////////////////////////////////////////////////////////////////////////
///std::list::begin.

TGL5DPainter::SurfIter_t TGL5DPainter::SurfacesBegin()
{
   return fIsos.begin();
}

////////////////////////////////////////////////////////////////////////////////
///std::list::end.

TGL5DPainter::SurfIter_t TGL5DPainter::SurfacesEnd()
{
   return fIsos.end();
}

////////////////////////////////////////////////////////////////////////////////
///Initialize OpenGL state variables.

void TGL5DPainter::InitGL() const
{
   glEnable(GL_LIGHTING);
   glEnable(GL_LIGHT0);
   glEnable(GL_DEPTH_TEST);
   glDisable(GL_CULL_FACE);
   glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
}

////////////////////////////////////////////////////////////////////////////////
///Return some gl states to original values.

void TGL5DPainter::DeInitGL()const
{
   glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
   glDisable(GL_CULL_FACE);
   glDisable(GL_DEPTH_TEST);
   glDisable(GL_LIGHT0);
   glDisable(GL_LIGHTING);
}

////////////////////////////////////////////////////////////////////////////////
///Draw a set of meshes.

void TGL5DPainter::DrawPlot() const
{
   //Shift plot to point of origin.
   const Rgl::PlotTranslation trGuard(this);

   fBackBox.DrawBox(fSelectedPart, fSelectionPass, fZLevels, fHighColor);
   //
   if (!fIsos.size())
      DrawCloud();
   else {
      //Two passes. First, non-transparent surfaces.
      Bool_t needSecondPass = kFALSE;
      for (ConstSurfIter_t it = fIsos.begin(); it != fIsos.end(); ++it) {
         //
         if (it->fHide)
            continue;
         if (it->fAlpha != 100) {
            needSecondPass = kTRUE;
            continue;
         }
         if (!fSelectionPass)
            SetSurfaceColor(it);
         glEnable(GL_POLYGON_OFFSET_FILL);
         glPolygonOffset(1.f, 1.f);
         DrawMesh(it);
         glDisable(GL_POLYGON_OFFSET_FILL);

         if (!fSelectionPass && it->fHighlight) {
            const TGLDisableGuard lightGuard(GL_LIGHTING);
            const TGLEnableGuard  blendGuard(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            glColor4d(1., 0.4, 0., 0.5);
            DrawMesh(it);
         }
      }
      //Second pass - semi-transparent surfaces.
      if (needSecondPass) {
         const TGLEnableGuard  blendGuard(GL_BLEND);
         glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
         glDepthMask(GL_FALSE);
         for (ConstSurfIter_t it = fIsos.begin(); it != fIsos.end(); ++it) {
            //
            if (it->fAlpha == 100)
               continue;
            if (!fSelectionPass)
               SetSurfaceColor(it);

            glEnable(GL_POLYGON_OFFSET_FILL);
            glPolygonOffset(1.f, 1.f);
            DrawMesh(it);
            glDisable(GL_POLYGON_OFFSET_FILL);

            if (!fSelectionPass && it->fHighlight) {
               const TGLDisableGuard lightGuard(GL_LIGHTING);
               glColor4d(1., 0.4, 0., it->fAlpha / 150.);
               DrawMesh(it);
            }
         }
         glDepthMask(GL_TRUE);
      }
   }

   if (fBoxCut.IsActive())
      fBoxCut.DrawBox(fSelectionPass, fSelectedPart);
}

////////////////////////////////////////////////////////////////////////////////
///Set the color for iso-surface.

void TGL5DPainter::SetSurfaceColor(ConstSurfIter_t it)const
{
   Color_t ind = it->fColor;
   Float_t rgba[] = {0.f, 0.f, 0.f, static_cast<Float_t>(it->fAlpha / 100.)};
   Rgl::Pad::ExtractRGBA(ind, rgba);
   //Set color for surface.
   glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, rgba);
   const Float_t specColor[] = {1.f, 1.f, 1.f, 1.f};
   glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specColor);
   glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 20.f);
}

////////////////////////////////////////////////////////////////////////////////
///Draw full cloud of points.

void TGL5DPainter::DrawCloud()const
{
   const TGLDisableGuard light(GL_LIGHTING);
   const TGLDisableGuard depth(GL_DEPTH_TEST);

   glColor3d(0.4, 0., 1.);
   glPointSize(3.f);

   glBegin(GL_POINTS);

   const Double_t xs = fCoord->GetXScale();
   const Double_t ys = fCoord->GetYScale();
   const Double_t zs = fCoord->GetZScale();

   for (Int_t i = 0; i < fData->fNP; ++i)
      glVertex3d(fData->fV1[i] * xs, fData->fV2[i] * ys, fData->fV3[i] * zs);

   glEnd();

   glPointSize(1.f);
}

////////////////////////////////////////////////////////////////////////////////
///Draw cloud for selected iso-surface.

void TGL5DPainter::DrawSubCloud(Double_t v4, Double_t range, Color_t ci)const
{
   const TGLDisableGuard light(GL_LIGHTING);

   Float_t rgba[4] = {};
   Rgl::Pad::ExtractRGBA(ci, rgba);

   glColor3fv(rgba);
   glPointSize(3.f);

   glBegin(GL_POINTS);

   const Double_t xs = fCoord->GetXScale();
   const Double_t ys = fCoord->GetYScale();
   const Double_t zs = fCoord->GetZScale();

   for (Int_t i = 0; i < fData->fNP; ++i)
      if (TMath::Abs(fData->fV4[i] - v4) < range)
         glVertex3d(fData->fV1[i] * xs, fData->fV2[i] * ys, fData->fV3[i] * zs);

   glEnd();

   glPointSize(1.f);
}

////////////////////////////////////////////////////////////////////////////////
///Draw one iso-surface.

void TGL5DPainter::DrawMesh(ConstSurfIter_t surf)const
{
   const Mesh_t &m = surf->fMesh;

   if (!fBoxCut.IsActive()) {
      if (!fSelectionPass)
         Rgl::DrawMesh(m.fVerts, m.fNorms, m.fTris);
      else {
         Rgl::ObjectIDToColor(fSelectionBase, fHighColor);
         Rgl::DrawMesh(m.fVerts, m.fTris);
      }
   } else {
      if (!fSelectionPass) {
         Rgl::DrawMesh(m.fVerts, m.fNorms, m.fTris, fBoxCut);
      } else {
         Rgl::ObjectIDToColor(fSelectionBase, fHighColor);
         Rgl::DrawMesh(m.fVerts, m.fTris, fBoxCut);
      }
   }
}
