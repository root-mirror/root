// @(#)root/gl:$Id$
// Author:  Timur Pocheptsov  14/06/2006

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/
#include <algorithm>
#include <cctype>

#include "KeySymbols.h"
#include "TVirtualX.h"
#include "Buttons.h"
#include "TString.h"
#include "TColor.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TMath.h"
#include "TH1.h"

#include "TGLLegoPainter.h"
#include "TGLPlotCamera.h"
#include "TGLIncludes.h"

/** \class TGLLegoPainter
\ingroup opengl
Plot-painter implementing LEGO rendering of TH2 histograms in
cartesian, polar, cylindrical and spherical coordinates.
*/

ClassImp(TGLLegoPainter);

////////////////////////////////////////////////////////////////////////////////
///Ctor.

TGLLegoPainter::TGLLegoPainter(TH1 *hist, TGLPlotCamera *cam, TGLPlotCoordinates *coord)
                  : TGLPlotPainter(hist, cam, coord, kFALSE, kTRUE, kTRUE),
                    fLegoType(kColorSimple),
                    fMinZ(0.),
                    fDrawErrors(kFALSE)
{
}

////////////////////////////////////////////////////////////////////////////////
///Obtain bin's info (i, j, value).

char *TGLLegoPainter::GetPlotInfo(Int_t /*px*/, Int_t /*py*/)
{
   fBinInfo = "";

   if (fSelectedPart) {
      if (fSelectedPart < fSelectionBase) {
         if (fHist->Class())
            fBinInfo += fHist->Class()->GetName();
         fBinInfo += "::";
         fBinInfo += fHist->GetName();
      } else if (!fHighColor) {
         const Int_t binI = (fSelectedPart - fSelectionBase) / fCoord->GetNYBins() + fCoord->GetFirstXBin();
         const Int_t binJ = (fSelectedPart - fSelectionBase) % fCoord->GetNYBins() + fCoord->GetFirstYBin();
         fBinInfo.Form("(binx = %d; biny = %d; binc = %f)", binI, binJ,
                       fHist->GetBinContent(binI, binJ));
      } else
         fBinInfo = "Switch to true-color mode to obtain correct info";
   }

   return (Char_t *)fBinInfo.Data();
}

////////////////////////////////////////////////////////////////////////////////
///Select method.

Bool_t TGLLegoPainter::InitGeometry()
{
   Bool_t ret = kFALSE;
   switch (fCoord->GetCoordType()) {
   case kGLCartesian:
      ret = InitGeometryCartesian(); break;
   case kGLPolar:
      ret = InitGeometryPolar(); break;
   case kGLCylindrical:
      ret = InitGeometryCylindrical(); break;
   case kGLSpherical:
      ret = InitGeometrySpherical(); break;
   default:
      return kFALSE;
   }
   if (ret && fCamera) fCamera->SetViewVolume(fBackBox.Get3DBox());
   return ret;
}

////////////////////////////////////////////////////////////////////////////////
///Geometry for lego in cartesian coords.

Bool_t TGLLegoPainter::InitGeometryCartesian()
{
   if (!fCoord->SetRanges(fHist, fDrawErrors, kFALSE))
      return kFALSE;

   fBackBox.SetPlotBox(fCoord->GetXRangeScaled(), fCoord->GetYRangeScaled(), fCoord->GetZRangeScaled());

   //Find bin edges
   const Int_t nX = fCoord->GetNXBins();
   const Double_t barWidth = fHist->GetBarWidth(), barOffset = fHist->GetBarOffset();
   const TGLVertex3 *frame = fBackBox.Get3DBox();
   fXEdges.resize(nX);

   if (fCoord->GetXLog()) {
      for (Int_t i = 0, ir = fCoord->GetFirstXBin(); i < nX; ++i, ++ir) {
         const Double_t xWidth = fXAxis->GetBinWidth(ir);
         Double_t low = fXAxis->GetBinLowEdge(ir) + xWidth * barOffset;
         fXEdges[i].first  = TMath::Log10(low) * fCoord->GetXScale();
         fXEdges[i].second = TMath::Log10(low + xWidth * barWidth) * fCoord->GetXScale();
         if (fXEdges[i].second > frame[1].X())
            fXEdges[i].second = frame[1].X();
         if (fXEdges[i].first < frame[0].X())
            fXEdges[i].first = frame[0].X();
         if (fXEdges[i].second < frame[0].X())
            fXEdges[i].second = frame[0].X();
      }
   } else {
      for (Int_t i = 0, ir = fCoord->GetFirstXBin(); i < nX; ++i, ++ir) {
         const Double_t xWidth = fXAxis->GetBinWidth(ir);
         fXEdges[i].first  = (fXAxis->GetBinLowEdge(ir) + xWidth * barOffset) * fCoord->GetXScale();
         fXEdges[i].second = fXEdges[i].first + xWidth * barWidth * fCoord->GetXScale();
         if (fXEdges[i].second > frame[1].X())
            fXEdges[i].second = frame[1].X();
         if (fXEdges[i].first < frame[0].X())
            fXEdges[i].first = frame[0].X();
         if (fXEdges[i].second < frame[0].X())
            fXEdges[i].second = frame[0].X();
      }
   }

   const Int_t nY = fCoord->GetNYBins();
   fYEdges.resize(nY);

   if (fCoord->GetYLog()) {
      for (Int_t j = 0, jr = fCoord->GetFirstYBin(); j < nY; ++j, ++jr) {
         const Double_t yWidth = fYAxis->GetBinWidth(jr);
         Double_t low = fYAxis->GetBinLowEdge(jr) + yWidth * barOffset;
         fYEdges[j].first  = TMath::Log10(low) * fCoord->GetYScale();
         fYEdges[j].second = TMath::Log10(low + yWidth * barWidth) * fCoord->GetYScale();
         if (fYEdges[j].second > frame[2].Y())
            fYEdges[j].second = frame[2].Y();
         if (fYEdges[j].first < frame[0].Y())
            fYEdges[j].first = frame[0].Y();
         if (fYEdges[j].second < frame[0].Y())
            fYEdges[j].second = frame[0].Y();
      }
   } else {
      for (Int_t j = 0, jr = fCoord->GetFirstYBin(); j < nY; ++j, ++jr) {
         const Double_t yWidth = fYAxis->GetBinWidth(jr);
         fYEdges[j].first  = (fYAxis->GetBinLowEdge(jr) + yWidth * barOffset) * fCoord->GetYScale();
         fYEdges[j].second = fYEdges[j].first + yWidth * barWidth * fCoord->GetYScale();
         if (fYEdges[j].second > frame[2].Y())
            fYEdges[j].second = frame[2].Y();
         if (fYEdges[j].first < frame[0].Y())
            fYEdges[j].first = frame[0].Y();
         if (fYEdges[j].second < frame[0].Y())
            fYEdges[j].second = frame[0].Y();
      }
   }

   fMinZ = frame[0].Z();
   if (fMinZ < 0.)
      frame[4].Z() > 0. ? fMinZ = 0. : fMinZ = frame[4].Z();

   if (fCoord->Modified()) {
      fUpdateSelection = kTRUE;
      fXOZSectionPos = frame[0].Y();
      fYOZSectionPos = frame[0].X();
      fXOYSectionPos = frame[0].Z();
      fCoord->ResetModified();
      Rgl::SetZLevels(fZAxis, fCoord->GetZRange().first, fCoord->GetZRange().second, fCoord->GetZScale(), fZLevels);
   }

   fMinMaxVal.first  = fHist->GetBinContent(fCoord->GetFirstXBin(), fCoord->GetFirstYBin());
   fMinMaxVal.second = fMinMaxVal.first;

   for (Int_t i = fCoord->GetFirstXBin(), e = fCoord->GetLastXBin(); i <= e; ++i) {
      for (Int_t j = fCoord->GetFirstYBin(), e1 = fCoord->GetLastYBin(); j <= e1; ++j) {
         Double_t val = fHist->GetBinContent(i, j);
         fMinMaxVal.first  = TMath::Min(fMinMaxVal.first, val);
         fMinMaxVal.second = TMath::Max(fMinMaxVal.second, val);
      }
   }

   ClampZ(fMinMaxVal.first);
   ClampZ(fMinMaxVal.second);

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
///Geometry for lego in polar coords.

Bool_t TGLLegoPainter::InitGeometryPolar()
{
   if (!fCoord->SetRanges(fHist, kFALSE, kFALSE))
      return kFALSE;

   fBackBox.SetPlotBox(fCoord->GetXRangeScaled(), fCoord->GetYRangeScaled(), fCoord->GetZRangeScaled());

   if (fCoord->Modified()) {
      fUpdateSelection = kTRUE;
      fCoord->ResetModified();
   }

   const Int_t nY = fCoord->GetNYBins();//yBins.second - yBins.first + 1;
   fYEdges.resize(nY);

   for (Int_t j = 0, jr = fCoord->GetFirstYBin(); j < nY; ++j, ++jr) {
      fYEdges[j].first = ((fYAxis->GetBinLowEdge(jr)) - fCoord->GetYRange().first) /
                           fCoord->GetYLength() * fCoord->GetYScale();
      fYEdges[j].second = ((fYAxis->GetBinUpEdge(jr)) - fCoord->GetYRange().first) /
                           fCoord->GetYLength() * fCoord->GetYScale();
   }

   const Int_t nX = fCoord->GetNXBins();
   fCosSinTableX.resize(nX + 1);
   const Double_t fullAngle = fXAxis->GetXmax() - fXAxis->GetXmin();
   const Double_t phiLow = fXAxis->GetXmin();
   Double_t angle = 0;
   for (Int_t i = 0, ir = fCoord->GetFirstXBin(); i < nX; ++i, ++ir) {
      angle = (fXAxis->GetBinLowEdge(ir) - phiLow) / fullAngle * TMath::TwoPi();
      fCosSinTableX[i].first  = TMath::Cos(angle);
      fCosSinTableX[i].second = TMath::Sin(angle);
   }
   angle = (fXAxis->GetBinUpEdge(fCoord->GetLastXBin()) - phiLow) / fullAngle * TMath::TwoPi();
   fCosSinTableX[nX].first  = TMath::Cos(angle);
   fCosSinTableX[nX].second = TMath::Sin(angle);

   fMinZ = fBackBox.Get3DBox()[0].Z();
   if (fMinZ < 0.)
      fBackBox.Get3DBox()[4].Z() > 0. ? fMinZ = 0. : fMinZ = fBackBox.Get3DBox()[4].Z();

   fMinMaxVal.first  = fHist->GetBinContent(fCoord->GetFirstXBin(), fCoord->GetFirstYBin());
   fMinMaxVal.second = fMinMaxVal.first;

   for (Int_t i = fCoord->GetFirstXBin(), e = fCoord->GetLastXBin(); i <= e; ++i) {
      for (Int_t j = fCoord->GetFirstYBin(), e1 = fCoord->GetLastYBin(); j <= e1; ++j) {
         Double_t val = fHist->GetBinContent(i, j);
         fMinMaxVal.first  = TMath::Min(fMinMaxVal.first, val);
         fMinMaxVal.second = TMath::Max(fMinMaxVal.second, val);
      }
   }

   ClampZ(fMinMaxVal.first);
   ClampZ(fMinMaxVal.second);

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
///Geometry for lego in cylindrical coords.

Bool_t TGLLegoPainter::InitGeometryCylindrical()
{
   if (!fCoord->SetRanges(fHist, kFALSE, kFALSE))
      return kFALSE;

   fBackBox.SetPlotBox(fCoord->GetXRangeScaled(), fCoord->GetYRangeScaled(), fCoord->GetZRangeScaled());

   const Int_t nY = fCoord->GetNYBins();
   fYEdges.resize(nY);

   if (fCoord->GetYLog()) {
      for (Int_t j = 0, jr = fCoord->GetFirstYBin(); j < nY; ++j, ++jr) {
         fYEdges[j].first  = TMath::Log10(fYAxis->GetBinLowEdge(jr)) * fCoord->GetYScale();
         fYEdges[j].second = TMath::Log10(fYAxis->GetBinUpEdge(jr))  * fCoord->GetYScale();
      }
   } else {
      for (Int_t j = 0, jr = fCoord->GetFirstYBin(); j < nY; ++j, ++jr) {
         fYEdges[j].first  = fYAxis->GetBinLowEdge(jr) * fCoord->GetYScale();
         fYEdges[j].second = fYAxis->GetBinUpEdge(jr)  * fCoord->GetYScale();
      }
   }

   const Int_t nX = fCoord->GetNXBins();
   fCosSinTableX.resize(nX + 1);
   const Double_t fullAngle = fXAxis->GetXmax() - fXAxis->GetXmin();
   const Double_t phiLow    = fXAxis->GetXmin();
   Double_t angle = 0.;
   for (Int_t i = 0, ir = fCoord->GetFirstXBin(); i < nX; ++i, ++ir) {
      angle = (fXAxis->GetBinLowEdge(ir) - phiLow) / fullAngle * TMath::TwoPi();
      fCosSinTableX[i].first  = TMath::Cos(angle);
      fCosSinTableX[i].second = TMath::Sin(angle);
   }
   angle = (fXAxis->GetBinUpEdge(fCoord->GetLastXBin()) - phiLow) / fullAngle * TMath::TwoPi();
   fCosSinTableX[nX].first  = TMath::Cos(angle);
   fCosSinTableX[nX].second = TMath::Sin(angle);

   if (fCoord->Modified()) {
      fUpdateSelection = kTRUE;
      fCoord->ResetModified();
   }

   fMinZ = fCoord->GetZRange().first;
   if (fMinZ < 0.)
      fCoord->GetZRange().second > 0. ? fMinZ = 0. : fMinZ = fCoord->GetZRange().second;


   fMinMaxVal.first  = fHist->GetBinContent(fCoord->GetFirstXBin(), fCoord->GetFirstYBin());
   fMinMaxVal.second = fMinMaxVal.first;

   for (Int_t i = fCoord->GetFirstXBin(), e = fCoord->GetLastXBin(); i <= e; ++i) {
      for (Int_t j = fCoord->GetFirstYBin(), e1 = fCoord->GetLastYBin(); j <= e1; ++j) {
         Double_t val = fHist->GetBinContent(i, j);
         fMinMaxVal.first  = TMath::Min(fMinMaxVal.first, val);
         fMinMaxVal.second = TMath::Max(fMinMaxVal.second, val);
      }
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
///Geometry for lego in spherical coords.

Bool_t TGLLegoPainter::InitGeometrySpherical()
{
   if (!fCoord->SetRanges(fHist, kFALSE, kFALSE))
      return kFALSE;

   fBackBox.SetPlotBox(fCoord->GetXRangeScaled(), fCoord->GetYRangeScaled(), fCoord->GetZRangeScaled());

   const Int_t nY = fCoord->GetNYBins();
   fCosSinTableY.resize(nY + 1);
   const Double_t fullTheta = fYAxis->GetXmax() - fYAxis->GetXmin();
   const Double_t thetaLow  = fYAxis->GetXmin();
   Double_t angle = 0.;
   for (Int_t j = 0, jr = fCoord->GetFirstYBin(); j < nY; ++j, ++jr) {
      angle = (fYAxis->GetBinLowEdge(jr) - thetaLow) / fullTheta * TMath::Pi();
      fCosSinTableY[j].first = TMath::Cos(angle);
      fCosSinTableY[j].second = TMath::Sin(angle);
   }
   angle = (fYAxis->GetBinUpEdge(fCoord->GetLastYBin()) - thetaLow) / fullTheta * TMath::Pi();
   fCosSinTableY[nY].first = TMath::Cos(angle);
   fCosSinTableY[nY].second = TMath::Sin(angle);

   const Int_t nX = fCoord->GetNXBins();
   fCosSinTableX.resize(nX + 1);
   const Double_t fullPhi = fXAxis->GetXmax() - fXAxis->GetXmin();
   const Double_t phiLow  = fXAxis->GetXmin();

   for (Int_t i = 0, ir = fCoord->GetFirstXBin(); i < nX; ++i, ++ir) {
      angle = (fXAxis->GetBinLowEdge(ir) - phiLow) / fullPhi * TMath::TwoPi();
      fCosSinTableX[i].first  = TMath::Cos(angle);
      fCosSinTableX[i].second = TMath::Sin(angle);
   }

   angle = (fXAxis->GetBinUpEdge(fCoord->GetLastXBin()) - phiLow) / fullPhi * TMath::TwoPi();
   fCosSinTableX[nX].first  = TMath::Cos(angle);
   fCosSinTableX[nX].second = TMath::Sin(angle);

   fMinZ = fCoord->GetZRange().first;
   if (fMinZ < 0.)
      fCoord->GetZRange().second > 0. ? fMinZ = 0. : fMinZ = fCoord->GetZRange().second;

   fMinMaxVal.first  = fHist->GetBinContent(fCoord->GetFirstXBin(), fCoord->GetFirstYBin());
   fMinMaxVal.second = fMinMaxVal.first;

   for (Int_t i = fCoord->GetFirstXBin(), e = fCoord->GetLastXBin(); i <= e; ++i) {
      for (Int_t j = fCoord->GetFirstYBin(), e1 = fCoord->GetLastYBin(); j <= e1; ++j) {
         Double_t val = fHist->GetBinContent(i, j);
         fMinMaxVal.first  = TMath::Min(fMinMaxVal.first, val);
         fMinMaxVal.second = TMath::Max(fMinMaxVal.second, val);
      }
   }


   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
///User clicks on a lego with middle mouse button (middle for pad).

void TGLLegoPainter::StartPan(Int_t px, Int_t py)
{
   fMousePosition.fX = px;
   fMousePosition.fY = fCamera->GetHeight() - py;
   fCamera->StartPan(px, py);
   fBoxCut.StartMovement(px, py);
}

////////////////////////////////////////////////////////////////////////////////
///Move lego or section.

void TGLLegoPainter::Pan(Int_t px, Int_t py)
{
   if (fSelectedPart >= fSelectionBase || fSelectedPart == 1) {
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
         if (fBoxCut.IsActive() && (fSelectedPart >= kXAxis && fSelectedPart <= kZAxis))
            fBoxCut.MoveBox(px, py, fSelectedPart);
         else
            MoveSection(px, py);
      } else
         MoveSection(px, py);

      RestoreProjectionMatrix();
      RestoreModelviewMatrix();
   }

   fMousePosition.fX = px, fMousePosition.fY = py;
   fUpdateSelection = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
///Parse additional options.

void TGLLegoPainter::AddOption(const TString &option)
{
   using namespace std;
   const Ssiz_t legoPos = option.Index("lego");//"lego" _already_ _exists_ in a string.
   if (legoPos + 4 < option.Length() && isdigit(option[legoPos + 4])) {
      switch (option[legoPos + 4] - '0') {
      case 1:
         fLegoType = kColorSimple;
         break;
      case 2:
         fLegoType = kColorLevel;
         break;
      case 3:
         fLegoType = kCylindricBars;
         break;
      default:
         fLegoType = kColorSimple;
         break;
      }
   } else
      fLegoType = kColorSimple;
   //check 'e' option
   Ssiz_t ePos = option.Index("e");
   if (ePos == legoPos + 1)
      ePos = option.Index("e", legoPos + 4);
   fDrawErrors = ePos != kNPOS ? kTRUE : kFALSE;

   option.Index("z") == kNPOS ? fDrawPalette = kFALSE : fDrawPalette = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
///Initialize some gl state variables.

void TGLLegoPainter::InitGL()const
{
   glEnable(GL_DEPTH_TEST);
   glEnable(GL_LIGHTING);
   glEnable(GL_LIGHT0);
   //For lego, back polygons are culled (but not for sections).
   glEnable(GL_CULL_FACE);
   glCullFace(GL_BACK);

   glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
}

////////////////////////////////////////////////////////////////////////////////
///Return some gl states to original values.

void TGLLegoPainter::DeInitGL()const
{
   glDisable(GL_DEPTH_TEST);
   glDisable(GL_LIGHTING);
   glDisable(GL_LIGHT0);
   glDisable(GL_CULL_FACE);
   glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
}

////////////////////////////////////////////////////////////////////////////////
///Select method corresponding to coordinate system.

void TGLLegoPainter::DrawPlot()const
{
   //Shift plot to point of origin.
   const Rgl::PlotTranslation trGuard(this);

   switch (fCoord->GetCoordType()) {
   case kGLCartesian:
      return DrawLegoCartesian();
   case kGLPolar:
      return DrawLegoPolar();
   case kGLCylindrical:
      return DrawLegoCylindrical();
   case kGLSpherical:
      return DrawLegoSpherical();
   default:;
   }
}

////////////////////////////////////////////////////////////////////////////////
///Lego in cartesian system.

void TGLLegoPainter::DrawLegoCartesian()const
{
   if (fCoord->GetCoordType() == kGLCartesian) {
      fBackBox.DrawBox(fSelectedPart, fSelectionPass, fZLevels, fHighColor);
      const TGLDisableGuard cullGuard(GL_CULL_FACE);
      DrawSections();
   }

   //const TGLDisableGuard depthTest(GL_DEPTH_TEST); //[0-0]

   if (!fSelectionPass) {
      glEnable(GL_POLYGON_OFFSET_FILL);//[0
      glPolygonOffset(1.f, 1.f);
      SetLegoColor();
      if (fXOZSectionPos > fBackBox.Get3DBox()[0].Y() || fYOZSectionPos> fBackBox.Get3DBox()[0].X()) {
         //Lego is semi-transparent if we have any sections.
         glEnable(GL_BLEND);//[1
         glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      }
   }

   //Using front point, find the correct order to draw bars from
   //back to front (it's important only for semi-transparent lego).
   //Only in cartesian.
   const Int_t nX = fXEdges.size();
   const Int_t nY = fYEdges.size();
   const Int_t frontPoint = fBackBox.GetFrontPoint();
   Int_t iInit = 0, jInit = 0, irInit = fCoord->GetFirstXBin(), jrInit = fCoord->GetFirstYBin();
   const Int_t addI = frontPoint == 2 || frontPoint == 1 ? 1 : (iInit = nX - 1, irInit = fCoord->GetLastXBin(), -1);
   const Int_t addJ = frontPoint == 2 || frontPoint == 3 ? 1 : (jInit = nY - 1, jrInit = fCoord->GetLastYBin(), -1);

   if (fLegoType == kColorLevel && !fSelectionPass) {
      if (!PreparePalette()) {
         fLegoType = kColorSimple;
         fDrawPalette = kFALSE;
      } else
         fPalette.EnableTexture(GL_MODULATE);
   }

   if (fSelectionPass && fHighColor)
      Rgl::ObjectIDToColor(fSelectionBase, kTRUE);

   for(Int_t i = iInit, ir = irInit; addI > 0 ? i < nX : i >= 0; i += addI, ir += addI) {
      for(Int_t j = jInit, jr = jrInit; addJ > 0 ? j < nY : j >= 0; j += addJ, jr += addJ) {
         Double_t zMax = fHist->GetBinContent(ir, jr) * fCoord->GetFactor();
         if (!ClampZ(zMax))
            continue;

         const Int_t binID = fSelectionBase + i * fCoord->GetNYBins() + j;

         if (fSelectionPass && !fHighColor)
            Rgl::ObjectIDToColor(binID, kFALSE);
         else if(!fHighColor && fSelectedPart == binID)
            glMaterialfv(GL_FRONT, GL_EMISSION, Rgl::gOrangeEmission);

         if (fLegoType == kCylindricBars) {
            Rgl::DrawCylinder(&fQuadric, fXEdges[i].first, fXEdges[i].second, fYEdges[j].first,
                              fYEdges[j].second, fMinZ, zMax);
         } else if (fLegoType == kColorLevel && !fSelectionPass) {
            Rgl::DrawBoxFrontTextured(fXEdges[i].first, fXEdges[i].second, fYEdges[j].first,
                                      fYEdges[j].second, fMinZ, zMax, fPalette.GetTexCoord(fMinZ),
                                      fPalette.GetTexCoord(zMax), frontPoint);
         } else {
            Rgl::DrawBoxFront(fXEdges[i].first, fXEdges[i].second, fYEdges[j].first,
                              fYEdges[j].second, fMinZ, zMax, frontPoint);
         }

         if (!fHighColor && !fSelectionPass && fSelectedPart == binID)
            glMaterialfv(GL_FRONT, GL_EMISSION, Rgl::gNullEmission);
      }
   }

   if (fLegoType == kColorLevel && !fSelectionPass)
      fPalette.DisableTexture();

   //Draw outlines for non-cylindrical bars.
   if (!fSelectionPass) {
      glDisable(GL_POLYGON_OFFSET_FILL);//0]
      const TGLDisableGuard lightGuard(GL_LIGHTING);//[2 - 2]
      if (fXOZSectionPos <= fBackBox.Get3DBox()[0].Y() && fYOZSectionPos <= fBackBox.Get3DBox()[0].X())
         glColor3d(0., 0., 0.);
      else
         glColor4d(0., 0., 0., 0.4);
      glPolygonMode(GL_FRONT, GL_LINE);//[3

      const TGLEnableGuard blendGuard(GL_BLEND);//[4-4] + 1]
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      const TGLEnableGuard smoothGuard(GL_LINE_SMOOTH);//[5-5]
      glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

      for(Int_t i = iInit, ir = irInit; addI > 0 ? i < nX : i >= 0; i += addI, ir += addI) {
         for(Int_t j = jInit, jr = jrInit; addJ > 0 ? j < nY : j >= 0; j += addJ, jr += addJ) {
            Double_t zMax = fHist->GetBinContent(ir, jr) * fCoord->GetFactor();
            if (!ClampZ(zMax))
               continue;
            if (fLegoType != kCylindricBars) {
               Rgl::DrawBoxFront(fXEdges[i].first, fXEdges[i].second, fYEdges[j].first,
                                 fYEdges[j].second, fMinZ, zMax, frontPoint);
            }
            if (fDrawErrors && zMax > 0.) {
               Double_t errorZMax = (fHist->GetBinContent(ir, jr) + fHist->GetBinError(ir, jr)) * fCoord->GetFactor();
               ClampZ(errorZMax);
               Rgl::DrawError(fXEdges[i].first, fXEdges[i].second, fYEdges[j].first,
                              fYEdges[j].second, zMax, errorZMax);
            }
         }
      }

      glPolygonMode(GL_FRONT, GL_FILL);//3]
   }

   if(!fSelectionPass && fDrawPalette)
      DrawPalette();
}

////////////////////////////////////////////////////////////////////////////////
///Lego in polar system.

void TGLLegoPainter::DrawLegoPolar()const
{
   const Int_t nX = fCosSinTableX.size() - 1;
   const Int_t nY = fYEdges.size();

   if (!fSelectionPass) {
      SetLegoColor();
      glEnable(GL_POLYGON_OFFSET_FILL);
      glPolygonOffset(1.f, 1.f);
   }

   Double_t points[4][2] = {};

   if (fLegoType == kColorLevel && !fSelectionPass) {
      if (!PreparePalette()) {
         fLegoType = kColorSimple;
         fDrawPalette = kFALSE;
      } else
         fPalette.EnableTexture(GL_MODULATE);
   }

   if (fHighColor && fSelectionPass)
      Rgl::ObjectIDToColor(fSelectionBase, kTRUE);

   for(Int_t i = 0, ir = fCoord->GetFirstXBin(); i < nX; ++i, ++ir) {
      for(Int_t j = 0, jr = fCoord->GetFirstYBin(); j < nY; ++j, ++jr) {
         Double_t zMax = fHist->GetBinContent(ir, jr);
         if (!ClampZ(zMax))
            continue;
         points[0][0] = fYEdges[j].first  * fCosSinTableX[i].first;
         points[0][1] = fYEdges[j].first  * fCosSinTableX[i].second;
         points[1][0] = fYEdges[j].second * fCosSinTableX[i].first;
         points[1][1] = fYEdges[j].second * fCosSinTableX[i].second;
         points[2][0] = fYEdges[j].second * fCosSinTableX[i + 1].first;
         points[2][1] = fYEdges[j].second * fCosSinTableX[i + 1].second;
         points[3][0] = fYEdges[j].first  * fCosSinTableX[i + 1].first;
         points[3][1] = fYEdges[j].first  * fCosSinTableX[i + 1].second;

         const Int_t binID = fSelectionBase + i * fCoord->GetNYBins() + j;

         if (!fHighColor && fSelectionPass)
            Rgl::ObjectIDToColor(binID, kFALSE);
         else if(!fHighColor && fSelectedPart == binID)
            glMaterialfv(GL_FRONT, GL_EMISSION, Rgl::gOrangeEmission);

         if (fLegoType == kColorLevel && !fSelectionPass)
            Rgl::DrawTrapezoidTextured(points, fMinZ, zMax, fPalette.GetTexCoord(fMinZ),
                                       fPalette.GetTexCoord(zMax));
         else
            Rgl::DrawTrapezoid(points, fMinZ, zMax);

         if (!fHighColor && !fSelectionPass && fSelectedPart == binID)
            glMaterialfv(GL_FRONT, GL_EMISSION, Rgl::gNullEmission);
      }
   }

   if (fLegoType == kColorLevel && !fSelectionPass)
      fPalette.DisableTexture();

   //Draw outlines.
   if (!fSelectionPass) {
      glDisable(GL_POLYGON_OFFSET_FILL);//0]
      const TGLDisableGuard lightGuard(GL_LIGHTING);//[2-2]
      glColor3d(0., 0., 0.);
      glPolygonMode(GL_FRONT, GL_LINE);//[3
      const TGLEnableGuard blendGuard(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      const TGLEnableGuard smoothGuard(GL_LINE_SMOOTH);
      glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

      for(Int_t i = 0, ir = fCoord->GetFirstXBin(); i < nX; ++i, ++ir) {
         for(Int_t j = 0, jr = fCoord->GetFirstYBin(); j < nY; ++j, ++jr) {
            Double_t zMax = fHist->GetBinContent(ir, jr);
            if (!ClampZ(zMax))
               continue;
            points[0][0] = fYEdges[j].first  * fCosSinTableX[i].first;
            points[0][1] = fYEdges[j].first  * fCosSinTableX[i].second;
            points[1][0] = fYEdges[j].second * fCosSinTableX[i].first;
            points[1][1] = fYEdges[j].second * fCosSinTableX[i].second;
            points[2][0] = fYEdges[j].second * fCosSinTableX[i + 1].first;
            points[2][1] = fYEdges[j].second * fCosSinTableX[i + 1].second;
            points[3][0] = fYEdges[j].first  * fCosSinTableX[i + 1].first;
            points[3][1] = fYEdges[j].first  * fCosSinTableX[i + 1].second;
            Rgl::DrawTrapezoid(points, fMinZ, zMax, kFALSE);
         }
      }

      glPolygonMode(GL_FRONT, GL_FILL);//3]
   }

   if(!fSelectionPass && fDrawPalette)
      DrawPalette();
}

////////////////////////////////////////////////////////////////////////////////
///Lego in cylindrical system.

void TGLLegoPainter::DrawLegoCylindrical()const
{
   const Int_t nX = fCosSinTableX.size() - 1;
   const Int_t nY = fYEdges.size();
   Double_t legoR = gStyle->GetLegoInnerR();
   if (legoR > 1. || legoR < 0.)
      legoR = 0.5;
   const Double_t rRange = fCoord->GetZLength();

   if (!fSelectionPass) {
      SetLegoColor();
      glEnable(GL_POLYGON_OFFSET_FILL);
      glPolygonOffset(1.f, 1.f);
   }

   Double_t points[4][2] = {};
   const Double_t sc = (1 - legoR) * fCoord->GetXScale();
   legoR *= fCoord->GetXScale();

   if (fLegoType == kColorLevel && !fSelectionPass) {
      if (!PreparePalette()) {
         fLegoType = kColorSimple;
         fDrawPalette = kFALSE;
      } else
         fPalette.EnableTexture(GL_MODULATE);
   }

   if (fHighColor && fSelectionPass)
      Rgl::ObjectIDToColor(fSelectionBase, kTRUE);

   for(Int_t i = 0, ir = fCoord->GetFirstXBin(); i < nX; ++i, ++ir) {
      for(Int_t j = 0, jr = fCoord->GetFirstYBin(); j < nY; ++j, ++jr) {
         Double_t zMin = legoR + (fMinZ - fCoord->GetZRange().first) / rRange * sc;
         Double_t zMax = legoR + (fHist->GetBinContent(ir, jr) - fCoord->GetZRange().first) / rRange * sc;
         if (zMin > zMax)
            std::swap(zMin, zMax);

         points[0][0] = fCosSinTableX[i].first * zMin;
         points[0][1] = fCosSinTableX[i].second * zMin;
         points[1][0] = fCosSinTableX[i].first * zMax;
         points[1][1] = fCosSinTableX[i].second * zMax;
         points[2][0] = fCosSinTableX[i + 1].first * zMax;
         points[2][1] = fCosSinTableX[i + 1].second * zMax;
         points[3][0] = fCosSinTableX[i + 1].first * zMin;
         points[3][1] = fCosSinTableX[i + 1].second * zMin;

         const Int_t binID = fSelectionBase + i * fCoord->GetNYBins() + j;

         if (fSelectionPass && !fHighColor)
            Rgl::ObjectIDToColor(binID, kFALSE);
         else if(!fHighColor && fSelectedPart == binID)
            glMaterialfv(GL_FRONT, GL_EMISSION, Rgl::gOrangeEmission);

         if (fLegoType == kColorLevel && !fSelectionPass){
            Rgl::DrawTrapezoidTextured2(points, fYEdges[j].first, fYEdges[j].second,
                                        fPalette.GetTexCoord(fMinZ), fPalette.GetTexCoord(fHist->GetBinContent(ir, jr)));
         }else
            Rgl::DrawTrapezoid(points, fYEdges[j].first, fYEdges[j].second);

         if(!fSelectionPass && !fHighColor && fSelectedPart == binID)
            glMaterialfv(GL_FRONT, GL_EMISSION, Rgl::gNullEmission);
      }
   }

   if (fLegoType == kColorLevel && !fSelectionPass)
      fPalette.DisableTexture();

   //Draw outlines.
   if (!fSelectionPass) {
      glDisable(GL_POLYGON_OFFSET_FILL);//0]
      const TGLDisableGuard lightGuard(GL_LIGHTING);//[2-2]
      glColor3d(0., 0., 0.);
      glPolygonMode(GL_FRONT, GL_LINE);//[3

      const TGLEnableGuard blendGuard(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      const TGLEnableGuard smoothGuard(GL_LINE_SMOOTH);
      glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

      for(Int_t i = 0, ir = fCoord->GetFirstXBin(); i < nX; ++i, ++ir) {
         for(Int_t j = 0, jr = fCoord->GetFirstYBin(); j < nY; ++j, ++jr) {
            Double_t zMin = legoR + (fMinZ - fCoord->GetZRange().first) / rRange * sc;
            Double_t zMax = legoR + (fHist->GetBinContent(ir, jr) - fCoord->GetZRange().first) / rRange * sc;
            if (zMin > zMax)
               std::swap(zMin, zMax);

            points[0][0] = fCosSinTableX[i].first * zMin;
            points[0][1] = fCosSinTableX[i].second * zMin;
            points[1][0] = fCosSinTableX[i].first * zMax;
            points[1][1] = fCosSinTableX[i].second * zMax;
            points[2][0] = fCosSinTableX[i + 1].first * zMax;
            points[2][1] = fCosSinTableX[i + 1].second * zMax;
            points[3][0] = fCosSinTableX[i + 1].first * zMin;
            points[3][1] = fCosSinTableX[i + 1].second * zMin;
            Rgl::DrawTrapezoid(points, fYEdges[j].first, fYEdges[j].second);
         }
      }

      glPolygonMode(GL_FRONT, GL_FILL);//3]
   }

   if(!fSelectionPass && fDrawPalette)
      DrawPalette();
}

////////////////////////////////////////////////////////////////////////////////
///Lego in spherical system.

void TGLLegoPainter::DrawLegoSpherical()const
{
   const Int_t nX = fCosSinTableX.size() - 1;
   const Int_t nY = fCosSinTableY.size() - 1;
   const Double_t rRange = fCoord->GetZLength();
   Double_t legoR = gStyle->GetLegoInnerR();
   if (legoR > 1. || legoR < 0.)
      legoR = 0.5;

   if (!fSelectionPass) {
      SetLegoColor();
      glEnable(GL_POLYGON_OFFSET_FILL);
      glPolygonOffset(1.f, 1.f);
   }

   Double_t points[8][3] = {};
   const Double_t sc = 1 - legoR;

   if (fLegoType == kColorLevel && !fSelectionPass) {
      if (!PreparePalette()) {
         fLegoType = kColorSimple;
         fDrawPalette = kFALSE;
      } else
         fPalette.EnableTexture(GL_MODULATE);
   }

   if (fSelectionPass && fHighColor)
      Rgl::ObjectIDToColor(fSelectionBase, kTRUE);

   for(Int_t i = 0, ir = fCoord->GetFirstXBin(); i < nX; ++i, ++ir) {
      for(Int_t j = 0, jr = fCoord->GetFirstYBin(); j < nY; ++j, ++jr) {
         Double_t zMin = legoR + (fMinZ - fCoord->GetZRange().first) / rRange * sc;
         Double_t zMax = legoR + (fHist->GetBinContent(ir, jr) - fCoord->GetZRange().first) / rRange * sc;
         if (zMin > zMax)
            std::swap(zMin, zMax);

         points[4][0] = zMin * fCosSinTableY[j].second * fCosSinTableX[i].first;
         points[4][1] = zMin * fCosSinTableY[j].second * fCosSinTableX[i].second;
         points[4][2] = zMin * fCosSinTableY[j].first;
         points[5][0] = zMin * fCosSinTableY[j].second * fCosSinTableX[i + 1].first;
         points[5][1] = zMin * fCosSinTableY[j].second * fCosSinTableX[i + 1].second;
         points[5][2] = zMin * fCosSinTableY[j].first;
         points[6][0] = zMax * fCosSinTableY[j].second * fCosSinTableX[i + 1].first;
         points[6][1] = zMax * fCosSinTableY[j].second * fCosSinTableX[i + 1].second;
         points[6][2] = zMax * fCosSinTableY[j].first;
         points[7][0] = zMax * fCosSinTableY[j].second * fCosSinTableX[i].first;
         points[7][1] = zMax * fCosSinTableY[j].second * fCosSinTableX[i].second;
         points[7][2] = zMax * fCosSinTableY[j].first;
         points[0][0] = zMin * fCosSinTableY[j + 1].second * fCosSinTableX[i].first;
         points[0][1] = zMin * fCosSinTableY[j + 1].second * fCosSinTableX[i].second;
         points[0][2] = zMin * fCosSinTableY[j + 1].first;
         points[1][0] = zMin * fCosSinTableY[j + 1].second * fCosSinTableX[i + 1].first;
         points[1][1] = zMin * fCosSinTableY[j + 1].second * fCosSinTableX[i + 1].second;
         points[1][2] = zMin * fCosSinTableY[j + 1].first;
         points[2][0] = zMax * fCosSinTableY[j + 1].second * fCosSinTableX[i + 1].first;
         points[2][1] = zMax * fCosSinTableY[j + 1].second * fCosSinTableX[i + 1].second;
         points[2][2] = zMax * fCosSinTableY[j + 1].first;
         points[3][0] = zMax * fCosSinTableY[j + 1].second * fCosSinTableX[i].first;
         points[3][1] = zMax * fCosSinTableY[j + 1].second * fCosSinTableX[i].second;
         points[3][2] = zMax * fCosSinTableY[j + 1].first;

         const Int_t binID = fSelectionBase + i * fCoord->GetNYBins() + j;

         if (fSelectionPass && !fHighColor)
            Rgl::ObjectIDToColor(binID, kFALSE);
         else if(!fHighColor && fSelectedPart == binID)
            glMaterialfv(GL_FRONT, GL_EMISSION, Rgl::gOrangeEmission);
         if (fLegoType == kColorLevel && !fSelectionPass)
            Rgl::DrawTrapezoidTextured(points, fPalette.GetTexCoord(fMinZ),
                                       fPalette.GetTexCoord(fHist->GetBinContent(ir, jr)));
         else
            Rgl::DrawTrapezoid(points);

         if(!fHighColor && fSelectedPart == binID)
            glMaterialfv(GL_FRONT, GL_EMISSION, Rgl::gNullEmission);
      }
   }

   if (fLegoType == kColorLevel && !fSelectionPass)
      fPalette.DisableTexture();

   //Draw outlines.
   if (!fSelectionPass) {
      glDisable(GL_POLYGON_OFFSET_FILL);//0]
      const TGLDisableGuard lightGuard(GL_LIGHTING);//[2-2]
      glColor3d(0., 0., 0.);
      glPolygonMode(GL_FRONT, GL_LINE);//[3
      const TGLEnableGuard blendGuard(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      const TGLEnableGuard smoothGuard(GL_LINE_SMOOTH);
      glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

      for(Int_t i = 0, ir = fCoord->GetFirstXBin(); i < nX; ++i, ++ir) {
         for(Int_t j = 0, jr = fCoord->GetFirstYBin(); j < nY; ++j, ++jr) {
            Double_t zMin = legoR + (fMinZ - fCoord->GetZRange().first) / rRange * sc;
            Double_t zMax = legoR + (fHist->GetBinContent(ir, jr) - fCoord->GetZRange().first) / rRange * sc;
            if (zMin > zMax)
               std::swap(zMin, zMax);

            points[4][0] = zMin * fCosSinTableY[j].second * fCosSinTableX[i].first;
            points[4][1] = zMin * fCosSinTableY[j].second * fCosSinTableX[i].second;
            points[4][2] = zMin * fCosSinTableY[j].first;
            points[5][0] = zMin * fCosSinTableY[j].second * fCosSinTableX[i + 1].first;
            points[5][1] = zMin * fCosSinTableY[j].second * fCosSinTableX[i + 1].second;
            points[5][2] = zMin * fCosSinTableY[j].first;
            points[6][0] = zMax * fCosSinTableY[j].second * fCosSinTableX[i + 1].first;
            points[6][1] = zMax * fCosSinTableY[j].second * fCosSinTableX[i + 1].second;
            points[6][2] = zMax * fCosSinTableY[j].first;
            points[7][0] = zMax * fCosSinTableY[j].second * fCosSinTableX[i].first;
            points[7][1] = zMax * fCosSinTableY[j].second * fCosSinTableX[i].second;
            points[7][2] = zMax * fCosSinTableY[j].first;
            points[0][0] = zMin * fCosSinTableY[j + 1].second * fCosSinTableX[i].first;
            points[0][1] = zMin * fCosSinTableY[j + 1].second * fCosSinTableX[i].second;
            points[0][2] = zMin * fCosSinTableY[j + 1].first;
            points[1][0] = zMin * fCosSinTableY[j + 1].second * fCosSinTableX[i + 1].first;
            points[1][1] = zMin * fCosSinTableY[j + 1].second * fCosSinTableX[i + 1].second;
            points[1][2] = zMin * fCosSinTableY[j + 1].first;
            points[2][0] = zMax * fCosSinTableY[j + 1].second * fCosSinTableX[i + 1].first;
            points[2][1] = zMax * fCosSinTableY[j + 1].second * fCosSinTableX[i + 1].second;
            points[2][2] = zMax * fCosSinTableY[j + 1].first;
            points[3][0] = zMax * fCosSinTableY[j + 1].second * fCosSinTableX[i].first;
            points[3][1] = zMax * fCosSinTableY[j + 1].second * fCosSinTableX[i].second;
            points[3][2] = zMax * fCosSinTableY[j + 1].first;
            Rgl::DrawTrapezoid(points);
         }
      }

      glPolygonMode(GL_FRONT, GL_FILL);//3]
   }

   if(!fSelectionPass && fDrawPalette)
      DrawPalette();
}

////////////////////////////////////////////////////////////////////////////////
///Set lego's color.

void TGLLegoPainter::SetLegoColor()const
{
   Float_t diffColor[] = {0.8f, 0.8f, 0.8f, 0.15f};

   if (fLegoType != kColorLevel && fHist->GetFillColor() != kWhite)
      if (const TColor *c = gROOT->GetColor(fHist->GetFillColor()))
         c->GetRGB(diffColor[0], diffColor[1], diffColor[2]);

   glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffColor);
   const Float_t specColor[] = {1.f, 1.f, 1.f, 1.f};
   glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specColor);
   glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 70.f);
}

////////////////////////////////////////////////////////////////////////////////
///XOZ plane parallel section.

void TGLLegoPainter::DrawSectionXOZ()const
{
   Int_t binY = -1;

   for (Int_t i = 0, e = fYEdges.size(); i < e; ++i) {
      if (fYEdges[i].first <= fXOZSectionPos && fXOZSectionPos <= fYEdges[i].second) {
         binY = i;
         break;
      }
   }

   if (binY >= 0) {
      binY += fCoord->GetFirstYBin();
      glColor3d(1., 0., 0.);
      glLineWidth(3.f);
      //Draw 2d hist on the profile's plane.
      for (UInt_t i = 0, ir = fCoord->GetFirstXBin(), e = fXEdges.size(); i < e; ++i, ++ir) {
         Double_t zMax = fHist->GetBinContent(Int_t(ir), binY);
         if (!ClampZ(zMax))
            continue;

         glBegin(GL_LINE_LOOP);
         glVertex3d(fXEdges[i].first,  fXOZSectionPos, fMinZ);
         glVertex3d(fXEdges[i].first,  fXOZSectionPos, zMax);
         glVertex3d(fXEdges[i].second, fXOZSectionPos, zMax);
         glVertex3d(fXEdges[i].second, fXOZSectionPos, fMinZ);
         glEnd();
      }

      glLineWidth(1.f);
   }
}

////////////////////////////////////////////////////////////////////////////////
///YOZ plane parallel section.

void TGLLegoPainter::DrawSectionYOZ()const
{
   Int_t binX = -1;

   for (Int_t i = 0, e = fXEdges.size(); i < e; ++i) {
      if (fXEdges[i].first <= fYOZSectionPos && fYOZSectionPos <= fXEdges[i].second) {
         binX = i;
         break;
      }
   }

   if (binX >= 0) {
      binX += fCoord->GetFirstXBin();//fBinsX.first;
      glColor3d(1., 0., 0.);
      glLineWidth(3.f);
      //Draw 2d hist on the profile's plane.
      for (UInt_t i = 0, ir = fCoord->GetFirstYBin(), e = fYEdges.size(); i < e; ++i, ++ir) {
         Double_t zMax = fHist->GetBinContent(binX, ir);
         if (!ClampZ(zMax))
            continue;

         glBegin(GL_LINE_LOOP);
         glVertex3d(fYOZSectionPos, fYEdges[i].first,  fMinZ);
         glVertex3d(fYOZSectionPos, fYEdges[i].first,   zMax);
         glVertex3d(fYOZSectionPos, fYEdges[i].second,  zMax);
         glVertex3d(fYOZSectionPos, fYEdges[i].second, fMinZ);
         glEnd();
      }

      glLineWidth(1.f);
   }
}

////////////////////////////////////////////////////////////////////////////////
///Empty. No such sections for lego.

void TGLLegoPainter::DrawSectionXOY()const
{
}

////////////////////////////////////////////////////////////////////////////////
///Remove all sections and repaint.

void TGLLegoPainter::ProcessEvent(Int_t event, Int_t /*px*/, Int_t py)
{
   const TGLVertex3 *frame = fBackBox.Get3DBox();
   if (event == kButton1Double && (fXOZSectionPos > frame[0].Y() || fYOZSectionPos > frame[0].X())) {
      fXOZSectionPos = frame[0].Y();
      fYOZSectionPos = frame[0].X();
      if (fBoxCut.IsActive())
         fBoxCut.TurnOnOff();
      //gGLManager->PaintSingleObject(this);
      if (!gVirtualX->IsCmdThread())
         gROOT->ProcessLineFast(Form("((TGLPlotPainter *)0x%lx)->Paint()", (ULong_t)this));
      else
         Paint();
   } else if (event == kKeyPress && (py == kKey_c || py == kKey_C)) {
      Info("ProcessEvent", "Box cut does not exist for lego");
   }
}

////////////////////////////////////////////////////////////////////////////////
///Clamp z value.

Bool_t TGLLegoPainter::ClampZ(Double_t &zVal)const
{
   if (fCoord->GetZLog())
      if (zVal <= 0.)
         return kFALSE;
      else
         zVal = TMath::Log10(zVal) * fCoord->GetZScale();
   else
      zVal *= fCoord->GetZScale();

   const TGLVertex3 *frame = fBackBox.Get3DBox();

   if (zVal > frame[4].Z())
      zVal = frame[4].Z();
   else if (zVal < frame[0].Z())
      zVal = frame[0].Z();

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
///Initialize color palette.

Bool_t TGLLegoPainter::PreparePalette()const
{
   if(fMinMaxVal.first == fMinMaxVal.second)
      return kFALSE;//must be std::abs(fMinMaxVal.second - fMinMaxVal.first) < ...

   //User-defined contours are disabled, to be fixed in a future.
   if (fHist->TestBit(TH1::kUserContour))
      fHist->ResetBit(TH1::kUserContour);

   UInt_t paletteSize = gStyle->GetNumberContours();
   if (!paletteSize)
      paletteSize = 20;

   return fPalette.GeneratePalette(paletteSize, Rgl::Range_t(fMinZ, fMinMaxVal.second));
}

////////////////////////////////////////////////////////////////////////////////
///Draw. Palette.
///Originally, fCamera was never null.
///It can be a null now because of gl-viewer.

void TGLLegoPainter::DrawPalette()const
{
   if (!fCamera) {
      //Thank you, gl-viewer!
      return;
   }

   Rgl::DrawPalette(fCamera, fPalette);

   glFinish();

   fCamera->SetCamera();
   fCamera->Apply(fPadPhi, fPadTheta);
}

////////////////////////////////////////////////////////////////////////////////
///Draw. Palette. Axis.

void TGLLegoPainter::DrawPaletteAxis()const
{
   gVirtualX->SetDrawMode(TVirtualX::kCopy);//TCanvas by default sets in kInverse
   Rgl::DrawPaletteAxis(fCamera, fMinMaxVal, fCoord->GetCoordType() == kGLCartesian ? fCoord->GetZLog() : kFALSE);
}
