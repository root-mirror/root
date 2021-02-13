// @(#)root/gl:$Id$
// Author:  Timur Pocheptsov  31/08/2006

/*************************************************************************
 * Copyright (C) 1995-2006, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/
#include <cctype>

#include "KeySymbols.h"
#include "TVirtualX.h"
#include "Buttons.h"
#include "TString.h"
#include "TROOT.h"
#include "TMath.h"
#include "TColor.h"
#include "TH3.h"
#include "TVirtualMutex.h"

#include "TPolyMarker3D.h"
#include "TGLPlotCamera.h"
#include "TGLBoxPainter.h"
#include "TGLIncludes.h"

/** \class TGLBoxPainter
\ingroup opengl
Paints TH3 histograms by rendering variable-sized boxes matching the
bin contents.
*/

ClassImp(TGLBoxPainter);

////////////////////////////////////////////////////////////////////////////////
/// Normal constructor.

TGLBoxPainter::TGLBoxPainter(TH1 *hist, TGLPlotCamera *cam, TGLPlotCoordinates *coord)
                  : TGLPlotPainter(hist, cam, coord, kTRUE, kTRUE, kTRUE),
                    fXOZSlice("XOZ", (TH3 *)hist, coord, &fBackBox, TGLTH3Slice::kXOZ),
                    fYOZSlice("YOZ", (TH3 *)hist, coord, &fBackBox, TGLTH3Slice::kYOZ),
                    fXOYSlice("XOY", (TH3 *)hist, coord, &fBackBox, TGLTH3Slice::kXOY),
                    fType(kBox),
                    fPolymarker(0)
{
}


////////////////////////////////////////////////////////////////////////////////
/// Normal constructor.

TGLBoxPainter::TGLBoxPainter(TH1 *hist, TPolyMarker3D * pm,
                             TGLPlotCamera *cam, TGLPlotCoordinates *coord)
                  : TGLPlotPainter(hist, cam, coord, kFALSE, kFALSE, kFALSE),
                    fXOZSlice("XOZ", (TH3 *)hist, coord, &fBackBox, TGLTH3Slice::kXOZ),
                    fYOZSlice("YOZ", (TH3 *)hist, coord, &fBackBox, TGLTH3Slice::kYOZ),
                    fXOYSlice("XOY", (TH3 *)hist, coord, &fBackBox, TGLTH3Slice::kXOY),
                    fType(kBox),
                    fPolymarker(pm)
{
}

////////////////////////////////////////////////////////////////////////////////
///Show box info (i, j, k, binContent).

char *TGLBoxPainter::GetPlotInfo(Int_t, Int_t)
{
   fPlotInfo = "";

   if (fSelectedPart) {
      if (fSelectedPart < fSelectionBase) {
         if (fHist->Class())
            fPlotInfo += fHist->Class()->GetName();
         fPlotInfo += "::";
         fPlotInfo += fHist->GetName();
      } else if (!fHighColor){
         const Int_t arr2Dsize = fCoord->GetNYBins() * fCoord->GetNZBins();
         const Int_t binI = (fSelectedPart - fSelectionBase) / arr2Dsize + fCoord->GetFirstXBin();
         const Int_t binJ = (fSelectedPart - fSelectionBase) % arr2Dsize / fCoord->GetNZBins() + fCoord->GetFirstYBin();
         const Int_t binK = (fSelectedPart - fSelectionBase) % arr2Dsize % fCoord->GetNZBins() + fCoord->GetFirstZBin();

         fPlotInfo.Form("(binx = %d; biny = %d; binz = %d; binc = %f)", binI, binJ, binK,
                        fHist->GetBinContent(binI, binJ, binK));
      } else
         fPlotInfo = "Switch to true color mode to get correct info";
   }

   return (Char_t *)fPlotInfo.Data();
}

////////////////////////////////////////////////////////////////////////////////
///Set ranges, find min and max bin content.

Bool_t TGLBoxPainter::InitGeometry()
{
   fCoord->SetZLog(kFALSE);
   fCoord->SetYLog(kFALSE);
   fCoord->SetXLog(kFALSE);

   if (!fCoord->SetRanges(fHist, kFALSE, kTRUE))//kFALSE == drawErrors, kTRUE == zAsBins
      return kFALSE;

   fBackBox.SetPlotBox(fCoord->GetXRangeScaled(), fCoord->GetYRangeScaled(), fCoord->GetZRangeScaled());
   if(fCamera) fCamera->SetViewVolume(fBackBox.Get3DBox());

   fMinMaxVal.second  = fHist->GetBinContent(fCoord->GetFirstXBin(), fCoord->GetFirstYBin(), fCoord->GetFirstZBin());
   fMinMaxVal.first = fMinMaxVal.second;
   //Bad. You can up-date some bin value and get wrong picture.
   for (Int_t ir = fCoord->GetFirstXBin(); ir <= fCoord->GetLastXBin(); ++ir) {
      for (Int_t jr = fCoord->GetFirstYBin(); jr <= fCoord->GetLastYBin(); ++jr) {
         for (Int_t kr = fCoord->GetFirstZBin();  kr <= fCoord->GetLastZBin(); ++kr) {
            fMinMaxVal.second = TMath::Max(fMinMaxVal.second, fHist->GetBinContent(ir, jr, kr));
            fMinMaxVal.first = TMath::Min(fMinMaxVal.first, fHist->GetBinContent(ir, jr, kr));
         }
      }
   }

   fXOYSlice.SetMinMax(fMinMaxVal);
   fXOZSlice.SetMinMax(fMinMaxVal);
   fYOZSlice.SetMinMax(fMinMaxVal);

   if (fPolymarker) {
      const Double_t xScale = fCoord->GetXScale();
      const Double_t yScale = fCoord->GetYScale();
      const Double_t zScale = fCoord->GetZScale();

      fPMPoints.assign(fPolymarker->GetP(), fPolymarker->GetP() + fPolymarker->GetN() * 3);
      for (unsigned i = 0; i < fPMPoints.size(); i += 3) {
         fPMPoints[i] *= xScale;
         fPMPoints[i + 1] *= yScale;
         fPMPoints[i + 2] *= zScale;
      }
   }

   if (fCoord->Modified()) {
      fUpdateSelection = kTRUE;
      fXOZSectionPos = fBackBox.Get3DBox()[0].Y();
      fYOZSectionPos = fBackBox.Get3DBox()[0].X();
      fXOYSectionPos = fBackBox.Get3DBox()[0].Z();
      fCoord->ResetModified();
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// User clicks right mouse button (in a pad).

void TGLBoxPainter::StartPan(Int_t px, Int_t py)
{
   fMousePosition.fX = px;
   fMousePosition.fY = fCamera->GetHeight() - py;
   fCamera->StartPan(px, py);
   fBoxCut.StartMovement(px, fCamera->GetHeight() - py);
}

////////////////////////////////////////////////////////////////////////////////
/// User's moving mouse cursor, with middle mouse button pressed (for pad).
/// Calculate 3d shift related to 2d mouse movement.

void TGLBoxPainter::Pan(Int_t px, Int_t py)
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
         if (fBoxCut.IsActive() && (fSelectedPart >= kXAxis && fSelectedPart <= kZAxis))
            fBoxCut.MoveBox(px, py, fSelectedPart);
         else
            MoveSection(px, py);
      } else {
         MoveSection(px, py);
      }

      RestoreProjectionMatrix();
      RestoreModelviewMatrix();
   }

   fMousePosition.fX = px, fMousePosition.fY = py;
   fUpdateSelection = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Box1 == spheres.

void TGLBoxPainter::AddOption(const TString &option)
{
   using namespace std;//isdigit must be in std. But ...

   const Ssiz_t boxPos = option.Index("box");//"box" _already_ _exists_ in a string.
   if (boxPos + 3 < option.Length() && isdigit(option[boxPos + 3]))
      option[boxPos + 3] - '0' == 1 ? fType = kBox1 : fType = kBox;
   else
      fType = kBox;
   option.Index("z") == kNPOS ? fDrawPalette = kFALSE : fDrawPalette = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove sections.

void TGLBoxPainter::ProcessEvent(Int_t event, Int_t /*px*/, Int_t py)
{
   if (event == kButton1Double && (HasSections() || fBoxCut.IsActive())) {
      fXOZSectionPos = fBackBox.Get3DBox()[0].Y();
      fYOZSectionPos = fBackBox.Get3DBox()[0].X();
      fXOYSectionPos = fBackBox.Get3DBox()[0].Z();
      if (fBoxCut.IsActive())
         fBoxCut.TurnOnOff();
      if (!gVirtualX->IsCmdThread())
         gROOT->ProcessLineFast(Form("((TGLPlotPainter *)0x%lx)->Paint()", (ULong_t)this));
      else
         Paint();
   } else if (event == kKeyPress && (py == kKey_c || py == kKey_C)) {
      if (fHighColor)
         Info("ProcessEvent", "Switch to true color mode to use box cut");
      else {
         fBoxCut.TurnOnOff();
         fUpdateSelection = kTRUE;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize some gl state variables.

void TGLBoxPainter::InitGL()const
{
   glEnable(GL_DEPTH_TEST);
   glEnable(GL_LIGHTING);
   glEnable(GL_LIGHT0);
   //For box option back polygons are culled (but not for dynamic profiles).
   glEnable(GL_CULL_FACE);
   glCullFace(GL_BACK);

   glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
}

////////////////////////////////////////////////////////////////////////////////
///Return back some gl state variables.

void TGLBoxPainter::DeInitGL()const
{
   glDisable(GL_DEPTH_TEST);
   glDisable(GL_LIGHTING);
   glDisable(GL_LIGHT0);
   glDisable(GL_CULL_FACE);
   glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
}

namespace {

   /////////////////////////////////////////////////////////////////////////////
   ///

   void DrawMinusSigns(Double_t xMin, Double_t xMax, Double_t yMin, Double_t yMax,
                       Double_t zMin, Double_t zMax, Int_t fp, Bool_t onSphere, Bool_t transp)
   {
      const TGLDisableGuard depthTest(GL_DEPTH_TEST);
      const TGLDisableGuard cullFace(GL_CULL_FACE);

      const Double_t ratio  = onSphere ? 0.4 : 0.15;
      const Double_t leftX = xMin + ratio * (xMax - xMin), rightX = xMax - ratio * (xMax - xMin);
      const Double_t leftY = yMin + ratio * (yMax - yMin), rightY = yMax - ratio * (yMax - yMin);
      const Double_t lowZ = zMin / 2. + zMax / 2. - 0.1 * (zMax - zMin);
      const Double_t upZ = zMin / 2. + zMax / 2. + 0.1 * (zMax - zMin);


      const Double_t minusVerts[][3] = {{xMin, leftY, lowZ}, {xMin, leftY, upZ}, {xMin, rightY, upZ}, {xMin, rightY, lowZ},
                                        {leftX, yMin, lowZ}, {rightX, yMin, lowZ}, {rightX, yMin, upZ}, {leftX, yMin, upZ},
                                        {xMax, leftY, lowZ}, {xMax, rightY, lowZ}, {xMax, rightY, upZ}, {xMax, leftY, upZ},
                                        {rightX, yMax, lowZ}, {leftX, yMax, lowZ}, {leftX, yMax, upZ}, {rightX, yMax, upZ}};
      const Int_t minusQuads[][4] = {{0, 1, 2, 3}, {4, 5, 6, 7}, {8, 9, 10, 11}, {12, 13, 14, 15}};


      TGLDisableGuard light(GL_LIGHTING);
      glColor3d(1., 0., 0.);

      const Int_t    frontPlanes[][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};//Code duplication again :(
      const Int_t *verts = minusQuads[frontPlanes[fp][0]];

      glBegin(GL_POLYGON);
      glVertex3dv(minusVerts[verts[0]]);
      glVertex3dv(minusVerts[verts[1]]);
      glVertex3dv(minusVerts[verts[2]]);
      glVertex3dv(minusVerts[verts[3]]);
      glEnd();

      verts = minusQuads[frontPlanes[fp][1]];

      glBegin(GL_POLYGON);
      glVertex3dv(minusVerts[verts[0]]);
      glVertex3dv(minusVerts[verts[1]]);
      glVertex3dv(minusVerts[verts[2]]);
      glVertex3dv(minusVerts[verts[3]]);
      glEnd();

      const Float_t nullEmission[] = {0.f, 0.f, 0.f, 1.f};
      glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, nullEmission);
      glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, nullEmission);

      glColor4d(0., 0., 0., 0.25);
      glPolygonMode(GL_FRONT, GL_LINE);

      if (!transp) {
         glEnable(GL_BLEND);
         glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      }

      glEnable(GL_LINE_SMOOTH);
      glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

      verts = minusQuads[frontPlanes[fp][0]];

      glBegin(GL_POLYGON);
      glVertex3dv(minusVerts[verts[0]]);
      glVertex3dv(minusVerts[verts[1]]);
      glVertex3dv(minusVerts[verts[2]]);
      glVertex3dv(minusVerts[verts[3]]);
      glEnd();

      verts = minusQuads[frontPlanes[fp][1]];

      glBegin(GL_POLYGON);
      glVertex3dv(minusVerts[verts[0]]);
      glVertex3dv(minusVerts[verts[1]]);
      glVertex3dv(minusVerts[verts[2]]);
      glVertex3dv(minusVerts[verts[3]]);
      glEnd();

      glPolygonMode(GL_FRONT, GL_FILL);

      if (!transp)
         glDisable(GL_BLEND);
   }

}

////////////////////////////////////////////////////////////////////////////////

void TGLBoxPainter::DrawPlot()const
{
   if (fPolymarker)
      return DrawCloud();

   // Draw set of boxes (spheres)

   //Shift plot to point of origin.
   const Rgl::PlotTranslation trGuard(this);

   fBackBox.DrawBox(fSelectedPart, fSelectionPass, fZLevels, fHighColor);
   glDisable(GL_CULL_FACE);
   DrawSections();
   glEnable(GL_CULL_FACE);

   if (!fSelectionPass) {
      glEnable(GL_POLYGON_OFFSET_FILL);//[0
      glPolygonOffset(1.f, 1.f);
      SetPlotColor();
      if (HasSections()) {
         //Boxes are semi-transparent if we have any sections.
         glEnable(GL_BLEND);//[1
         glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      }
   }

   //Using front point, find the correct order to draw boxes from
   //back to front/from bottom to top (it's important only for semi-transparent boxes).
   const Int_t frontPoint = fBackBox.GetFrontPoint();
   Int_t irInit = fCoord->GetFirstXBin(), iInit = 0;
   const Int_t nX = fCoord->GetNXBins();
   Int_t jrInit = fCoord->GetFirstYBin(), jInit = 0;
   const Int_t nY = fCoord->GetNYBins();
   Int_t krInit = fCoord->GetFirstZBin(), kInit = 0;
   const Int_t nZ = fCoord->GetNZBins();

   const Int_t addI = frontPoint == 2 || frontPoint == 1 ? 1 : (iInit = nX - 1, irInit = fCoord->GetLastXBin(), -1);
   const Int_t addJ = frontPoint == 2 || frontPoint == 3 ? 1 : (jInit = nY - 1, jrInit = fCoord->GetLastYBin(), -1);
   const Int_t addK = fBackBox.Get2DBox()[frontPoint + 4].Y() < fBackBox.Get2DBox()[frontPoint].Y() ? 1
                     : (kInit = nZ - 1, krInit = fCoord->GetLastZBin(),-1);
   const Double_t xScale = fCoord->GetXScale();
   const Double_t yScale = fCoord->GetYScale();
   const Double_t zScale = fCoord->GetZScale();
   const TAxis   *xA = fXAxis;
   const TAxis   *yA = fYAxis;
   const TAxis   *zA = fZAxis;

   if (fSelectionPass && fHighColor)
      Rgl::ObjectIDToColor(fSelectionBase, fHighColor);//base + 1 == 7

   Double_t maxContent = TMath::Max(TMath::Abs(fMinMaxVal.first), TMath::Abs(fMinMaxVal.second));
   if(!maxContent)//bad, find better way to check zero.
      maxContent = 1.;

   Double_t wmin = TMath::Max(fHist->GetMinimum(),0.);
   Double_t wmax = TMath::Max(TMath::Abs(fHist->GetMaximum()),
                              TMath::Abs(fHist->GetMinimum()));
   Double_t binContent;

   for(Int_t ir = irInit, i = iInit; addI > 0 ? i < nX : i >= 0; ir += addI, i += addI) {
      for(Int_t jr = jrInit, j = jInit; addJ > 0 ? j < nY : j >= 0; jr += addJ, j += addJ) {
         for(Int_t kr = krInit, k = kInit; addK > 0 ? k < nZ : k >= 0; kr += addK, k += addK) {
            binContent = fHist->GetBinContent(ir, jr, kr);
            if (binContent < wmin) continue;
            if (binContent > wmax) binContent = wmax;

            const Double_t w = TMath::Power(TMath::Abs(binContent-wmin) / (wmax-wmin),1./3.);
            if (!w)
               continue;

            const Double_t xMin = xScale * (xA->GetBinLowEdge(ir) / 2 + xA->GetBinUpEdge(ir) / 2 - w * xA->GetBinWidth(ir) / 2);
            const Double_t xMax = xScale * (xA->GetBinLowEdge(ir) / 2 + xA->GetBinUpEdge(ir) / 2 + w * xA->GetBinWidth(ir) / 2);
            const Double_t yMin = yScale * (yA->GetBinLowEdge(jr) / 2 + yA->GetBinUpEdge(jr) / 2 - w * yA->GetBinWidth(jr) / 2);
            const Double_t yMax = yScale * (yA->GetBinLowEdge(jr) / 2 + yA->GetBinUpEdge(jr) / 2 + w * yA->GetBinWidth(jr) / 2);
            const Double_t zMin = zScale * (zA->GetBinLowEdge(kr) / 2 + zA->GetBinUpEdge(kr) / 2 - w * zA->GetBinWidth(kr) / 2);
            const Double_t zMax = zScale * (zA->GetBinLowEdge(kr) / 2 + zA->GetBinUpEdge(kr) / 2 + w * zA->GetBinWidth(kr) / 2);

            if (fBoxCut.IsActive() && fBoxCut.IsInCut(xMin, xMax, yMin, yMax, zMin, zMax))
               continue;

            const Int_t binID = fSelectionBase + i * fCoord->GetNZBins() * fCoord->GetNYBins() + j * fCoord->GetNZBins() + k;

            if (fSelectionPass && !fHighColor)
               Rgl::ObjectIDToColor(binID, fHighColor);
            else if(!fHighColor && fSelectedPart == binID)
               glMaterialfv(GL_FRONT, GL_EMISSION, Rgl::gOrangeEmission);

            if (fType == kBox) {
               Rgl::DrawBoxFront(xMin, xMax, yMin, yMax, zMin, zMax, frontPoint);
            } else {
               Rgl::DrawSphere(&fQuadric, xMin, xMax, yMin, yMax, zMin, zMax);
            }

            if (binContent < 0. && !fSelectionPass)
               DrawMinusSigns(xMin, xMax, yMin, yMax, zMin, zMax, frontPoint, fType != kBox, HasSections());

            if (!fSelectionPass && !fHighColor && fSelectedPart == binID)
               glMaterialfv(GL_FRONT, GL_EMISSION, Rgl::gNullEmission);
         }
      }
   }

   if (fBoxCut.IsActive())
      fBoxCut.DrawBox(fSelectionPass, fSelectedPart);

   if (!fSelectionPass && fType != kBox1) {
      glDisable(GL_POLYGON_OFFSET_FILL);//0]
      TGLDisableGuard lightGuard(GL_LIGHTING);//[2 - 2]
      glColor4d(0., 0., 0., 0.25);
      glPolygonMode(GL_FRONT, GL_LINE);//[3

      const TGLEnableGuard blendGuard(GL_BLEND);//[4-4] + 1]
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      const TGLEnableGuard smoothGuard(GL_LINE_SMOOTH);//[5-5]
      glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

      for(Int_t ir = irInit, i = iInit; addI > 0 ? i < nX : i >= 0; ir += addI, i += addI) {
         for(Int_t jr = jrInit, j = jInit; addJ > 0 ? j < nY : j >= 0; jr += addJ, j += addJ) {
            for(Int_t kr = krInit, k = kInit; addK > 0 ? k < nZ : k >= 0; kr += addK, k += addK) {
               binContent = fHist->GetBinContent(ir, jr, kr);
               if (binContent < wmin) continue;
               if (binContent > wmax) binContent = wmax;
               const Double_t w = TMath::Power(TMath::Abs(binContent-wmin) / (wmax-wmin),1./3.);
               if (!w)
                  continue;

               const Double_t xMin = xScale * (xA->GetBinLowEdge(ir) / 2 + xA->GetBinUpEdge(ir) / 2 - w * xA->GetBinWidth(ir) / 2);
               const Double_t xMax = xScale * (xA->GetBinLowEdge(ir) / 2 + xA->GetBinUpEdge(ir) / 2 + w * xA->GetBinWidth(ir) / 2);
               const Double_t yMin = yScale * (yA->GetBinLowEdge(jr) / 2 + yA->GetBinUpEdge(jr) / 2 - w * yA->GetBinWidth(jr) / 2);
               const Double_t yMax = yScale * (yA->GetBinLowEdge(jr) / 2 + yA->GetBinUpEdge(jr) / 2 + w * yA->GetBinWidth(jr) / 2);
               const Double_t zMin = zScale * (zA->GetBinLowEdge(kr) / 2 + zA->GetBinUpEdge(kr) / 2 - w * zA->GetBinWidth(kr) / 2);
               const Double_t zMax = zScale * (zA->GetBinLowEdge(kr) / 2 + zA->GetBinUpEdge(kr) / 2 + w * zA->GetBinWidth(kr) / 2);

               if (fBoxCut.IsActive() && fBoxCut.IsInCut(xMin, xMax, yMin, yMax, zMin, zMax))
                  continue;

               Rgl::DrawBoxFront(xMin, xMax, yMin, yMax, zMin, zMax, frontPoint);
            }
         }
      }

      glPolygonMode(GL_FRONT, GL_FILL);//3]
   }

   if (!fSelectionPass && fDrawPalette && HasSections())
      DrawPalette();
}

////////////////////////////////////////////////////////////////////////////////
///Draw a frame and a polymarker inside.

void TGLBoxPainter::DrawCloud()const
{
   //Shift plot to the point of origin.
   const Rgl::PlotTranslation trGuard(this);

   //Frame.
   fBackBox.DrawBox(fSelectedPart, fSelectionPass, fZLevels, fHighColor);

   if (fPhysicalShapeColor)
      glColor3fv(fPhysicalShapeColor);

   glDisable(GL_LIGHTING);

   const TGLVertex3 *bb = fBackBox.Get3DBox();
   const Double_t dX = (bb[1].X() - bb[0].X()) / 40.;
   const Double_t dY = (bb[3].Y() - bb[0].Y()) / 40.;
   const Double_t dZ = (bb[4].Z() - bb[0].Z()) / 40.;
   //Now, draw the cloud of points (polymarker) inside the frame.
   TGLUtil::RenderPolyMarkers(*fPolymarker, fPMPoints, dX, dY, dZ);

   glEnable(GL_LIGHTING);
}

////////////////////////////////////////////////////////////////////////////////
/// Set boxes color.

void TGLBoxPainter::SetPlotColor()const
{
   Float_t diffColor[] = {0.8f, 0.8f, 0.8f, 0.05f};

   if (fPhysicalShapeColor) {
      diffColor[0] = fPhysicalShapeColor[0];
      diffColor[1] = fPhysicalShapeColor[1];
      diffColor[2] = fPhysicalShapeColor[2];
   } else {
      if (fHist->GetFillColor() != kWhite)
         if (const TColor *c = gROOT->GetColor(fHist->GetFillColor()))
            c->GetRGB(diffColor[0], diffColor[1], diffColor[2]);
   }

   glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffColor);
   const Float_t specColor[] = {1.f, 1.f, 1.f, 1.f};
   glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specColor);
   glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 70.f);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw XOZ parallel section.

void TGLBoxPainter::DrawSectionXOZ()const
{
   if (fSelectionPass)
      return;
   fXOZSlice.DrawSlice(fXOZSectionPos / fCoord->GetYScale());
}

////////////////////////////////////////////////////////////////////////////////
/// Draw YOZ parallel section.

void TGLBoxPainter::DrawSectionYOZ()const
{
   if (fSelectionPass)
      return;
   fYOZSlice.DrawSlice(fYOZSectionPos / fCoord->GetXScale());
}

////////////////////////////////////////////////////////////////////////////////
/// Draw XOY parallel section.

void TGLBoxPainter::DrawSectionXOY()const
{
   if (fSelectionPass)
      return;
   fXOYSlice.DrawSlice(fXOYSectionPos / fCoord->GetZScale());
}

////////////////////////////////////////////////////////////////////////////////
/// Check, if any section exists.

Bool_t TGLBoxPainter::HasSections()const
{
   return fXOZSectionPos > fBackBox.Get3DBox()[0].Y() || fYOZSectionPos> fBackBox.Get3DBox()[0].X() ||
          fXOYSectionPos > fBackBox.Get3DBox()[0].Z();
}

////////////////////////////////////////////////////////////////////////////////
///Draw. Palette.
///Originally, fCamera was never null.
///It can be a null now because of gl-viewer.

void TGLBoxPainter::DrawPalette()const
{
   if (!fCamera) {
      //Thank you, gl-viewer!
      return;
   }

   const TGLLevelPalette * palette = 0;
   const TGLVertex3 *frame = fBackBox.Get3DBox();

   if (fXOZSectionPos > frame[0].Y())
      palette = &fXOZSlice.GetPalette();
   else if (fYOZSectionPos > frame[0].X())
      palette = &fYOZSlice.GetPalette();
   else if (fXOYSectionPos > frame[0].Z())
      palette = &fXOYSlice.GetPalette();

   if (!palette || !palette->GetPaletteSize()) {
      return;
   }

   Rgl::DrawPalette(fCamera, *palette);

   glFinish();

   fCamera->SetCamera();
   fCamera->Apply(fPadPhi, fPadTheta);
}

////////////////////////////////////////////////////////////////////////////////
///Draw. Palette. Axis.

void TGLBoxPainter::DrawPaletteAxis()const
{
   if (HasSections()) {
      gVirtualX->SetDrawMode(TVirtualX::kCopy);//TCanvas by default sets in kInverse
      Rgl::DrawPaletteAxis(fCamera, fMinMaxVal, fCoord->GetCoordType() == kGLCartesian ? fCoord->GetZLog() : kFALSE);
   }
}
