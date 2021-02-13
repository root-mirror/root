// @(#)root/eve:$Id$
// Author: Alja Mrak-Tadel 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TAxis.h"
#include "THLimitsFinder.h"

#include "TGLViewer.h"
#include "TGLIncludes.h"
#include "TGLPhysicalShape.h"
#include "TGLRnrCtx.h"
#include "TGLSelectRecord.h"
#include "TGLScene.h"
#include "TGLCamera.h"
#include "TGLUtil.h"
#include "TColor.h"
#include "TROOT.h"


#include "TEveCaloLegoGL.h"
#include "TEveCalo.h"
#include "TEveManager.h"
#include "TEveRGBAPalette.h"

#include <algorithm>

/** \class TEveCaloLegoGL
\ingroup TEve
OpenGL renderer class for TEveCaloLego.
*/

ClassImp(TEveCaloLegoGL);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TEveCaloLegoGL::TEveCaloLegoGL() :
   TGLObject(),

   fGridColor(-1),
   fFontColor(-1),

   fEtaAxis(0),
   fPhiAxis(0),
   fZAxis(0),
   fM(0),
   fDLCacheOK(kFALSE),
   fMaxVal(0),
   fValToPixel(0),
   fCurrentPixelsPerBin(0),
   fCells3D(kTRUE),
   fBinStep(-1)
{
   fDLCache = kFALSE;

   fEtaAxis = new TAxis();
   fPhiAxis = new TAxis();
   fZAxis   = new TAxis();

   fAxisPainter.SetFontMode(TGLFont::kPixmap);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TEveCaloLegoGL::~TEveCaloLegoGL()
{
   DLCachePurge();

   delete fEtaAxis;
   delete fPhiAxis;
   delete fZAxis;
}

////////////////////////////////////////////////////////////////////////////////
/// Set model object.

Bool_t TEveCaloLegoGL::SetModel(TObject* obj, const Option_t* /*opt*/)
{
   fM = SetModelDynCast<TEveCaloLego>(obj);
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set bounding box.

void TEveCaloLegoGL::SetBBox()
{
   SetAxisAlignedBBox(((TEveCaloLego*)fExternalObj)->AssertBBox());
}

////////////////////////////////////////////////////////////////////////////////
/// Drop all display-list definitions.

void TEveCaloLegoGL::DLCacheDrop()
{
   fDLCacheOK = kFALSE;
   for (SliceDLMap_i i = fDLMap.begin(); i != fDLMap.end(); ++i)
      i->second = 0;

   TGLObject::DLCacheDrop();
}

////////////////////////////////////////////////////////////////////////////////
/// Unregister all display-lists.

void TEveCaloLegoGL::DLCachePurge()
{
   // all lego cells
   fDLCacheOK = kFALSE;
   if (! fDLMap.empty()) {
      for (SliceDLMap_i i = fDLMap.begin(); i != fDLMap.end(); ++i) {
         if (i->second) {
            PurgeDLRange(i->second, 1);
            i->second = 0;
         }
      }
   }
   TGLObject::DLCachePurge();
}

////////////////////////////////////////////////////////////////////////////////
/// Draw an axis-aligned box using quads.
/// ~~~ {.cpp}
///    z
///    |
///    |
///    |________y
///   /  6-------7
///  /  /|      /|
/// x  5-------4 |
///    | 2-----|-3
///    |/      |/
///    1-------0
/// ~~~

void TEveCaloLegoGL::MakeQuad(Float_t x1, Float_t y1, Float_t z1,
      Float_t xw, Float_t yw, Float_t h) const
{
   Float_t x2 = x1 + xw;
   Float_t y2 = y1 + yw;
   Float_t z2 = z1 + h;

   if (x1 < fM->GetEtaMin()) x1 = fM->GetEtaMin();
   if (x2 > fM->GetEtaMax()) x2 = fM->GetEtaMax();

   if (y1 < fM->GetPhiMin()) y1 = fM->GetPhiMin();
   if (y2 > fM->GetPhiMax()) y2 = fM->GetPhiMax();

   glBegin(GL_QUADS);
   {
      // bottom 0123
      glNormal3f(0, 0, -1);
      glVertex3f(x2, y2, z1);
      glVertex3f(x2, y1, z1);
      glVertex3f(x1, y1, z1);
      glVertex3f(x1, y2, z1);
      // top 4765
      glNormal3f(0, 0, 1);
      glVertex3f(x2, y2, z2);
      glVertex3f(x1, y2, z2);
      glVertex3f(x1, y1, z2);
      glVertex3f(x2, y1, z2);

      // back 0451
      glNormal3f(1, 0, 0);
      glVertex3f(x2, y2, z1);
      glVertex3f(x2, y2, z2);
      glVertex3f(x2, y1, z2);
      glVertex3f(x2, y1, z1);
      // front 3267
      glNormal3f(-1, 0, 0);
      glVertex3f(x1, y2, z1);
      glVertex3f(x1, y1, z1);
      glVertex3f(x1, y1, z2);
      glVertex3f(x1, y2, z2);

      // left  0374
      glNormal3f(0, 1, 0);
      glVertex3f(x2, y2, z1);
      glVertex3f(x1, y2, z1);
      glVertex3f(x1, y2, z2);
      glVertex3f(x2, y2, z2);
      // right 1562
      glNormal3f(0, -1, 0);
      glVertex3f(x2, y1, z1);
      glVertex3f(x2, y1, z2);
      glVertex3f(x1, y1, z2);
      glVertex3f(x1, y1, z1);
   }
   glEnd();
}

////////////////////////////////////////////////////////////////////////////////
/// Create display-list that draws histogram bars for non-rebinned data.
/// It is used for filled and outline passes.

void TEveCaloLegoGL::Make3DDisplayList(TEveCaloData::vCellId_t& cellList, SliceDLMap_t& dlMap, Bool_t selection) const
{
   TEveCaloData::CellData_t cellData;
   Int_t   prevTower = 0;
   Float_t offset = 0;

   // ids in eta phi rng
   Int_t nSlices = fM->fData->GetNSlices();
   for (Int_t s = 0; s < nSlices; ++s)
   {
      if (dlMap.empty() || dlMap[s] == 0)
         dlMap[s] = glGenLists(1);

      glNewList(dlMap[s], GL_COMPILE);

      for (UInt_t i = 0; i < cellList.size(); ++i)
      {
         if (cellList[i].fSlice > s) continue;
         if (cellList[i].fTower != prevTower) {
            offset = 0;
            prevTower = cellList[i].fTower;
         }

         fM->fData->GetCellData(cellList[i], cellData);
         if (s == cellList[i].fSlice)
         {
            if (selection) glLoadName(i);

            WrapTwoPi(cellData.fPhiMin, cellData.fPhiMax);
            MakeQuad(cellData.EtaMin(), cellData.PhiMin(), offset,
                     cellData.EtaDelta(), cellData.PhiDelta(), cellData.Value(fM->fPlotEt));
         }
         offset += cellData.Value(fM->fPlotEt);
      }
      glEndList();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Create display-list that draws histogram bars for rebinned data.
/// It is used for filled and outline passes.

void TEveCaloLegoGL::Make3DDisplayListRebin(TEveCaloData::RebinData_t& rebinData, SliceDLMap_t& dlMap, Bool_t selection) const
{
   Int_t nSlices = fM->fData->GetNSlices();
   Float_t *vals;
   Float_t offset;
   Float_t y0, y1;

   for (Int_t s = 0; s < nSlices; ++s)
   {
      if (dlMap.empty() || dlMap[s] == 0)
         dlMap[s] = glGenLists(1);

      glNewList(dlMap[s], GL_COMPILE);

      for (Int_t i = 1; i <= fEtaAxis->GetNbins(); ++i)
      {
         for (Int_t j = 1; j <= fPhiAxis->GetNbins(); ++j)
         {
            const Int_t bin = (i)+(j)*(fEtaAxis->GetNbins()+2);

            if (rebinData.fBinData[bin] !=-1)
            {
               vals = rebinData.GetSliceVals(bin);
               offset =0;
               for (Int_t t = 0; t < s; ++t)
                  offset += vals[t];

               y0 = fPhiAxis->GetBinLowEdge(j);
               y1 = fPhiAxis->GetBinUpEdge(j);
               WrapTwoPi(y0, y1);

               if (selection) glLoadName(bin);

               MakeQuad(fEtaAxis->GetBinLowEdge(i), y0, offset,
                        fEtaAxis->GetBinWidth(i), y1-y0, vals[s]);
            }
         }
      }
      glEndList();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the axis 3D title position.

void TEveCaloLegoGL::SetAxis3DTitlePos(TGLRnrCtx &rnrCtx, Float_t x0, Float_t x1, Float_t y0, Float_t y1) const
{
   const GLdouble *pm = rnrCtx.RefCamera().RefLastNoPickProjM().CArr();
   GLdouble mm[16];
   GLint    vp[4];
   glGetDoublev(GL_MODELVIEW_MATRIX,  mm);
   glGetIntegerv(GL_VIEWPORT, vp);
   GLdouble projX[4], projY[4], projZ[4];

   GLdouble cornerX[4];
   GLdouble cornerY[4];
   cornerX[0] = x0; cornerY[0] = y0;
   cornerX[1] = x1; cornerY[1] = y0;
   cornerX[2] = x1; cornerY[2] = y1;
   cornerX[3] = x0; cornerY[3] = y1;

   gluProject(cornerX[0], cornerY[0], 0, mm, pm, vp, &projX[0], &projY[0], &projZ[0]);
   gluProject(cornerX[1], cornerY[1], 0, mm, pm, vp, &projX[1], &projY[1], &projZ[1]);
   gluProject(cornerX[2], cornerY[2], 0, mm, pm, vp, &projX[2], &projY[2], &projZ[2]);
   gluProject(cornerX[3], cornerY[3], 0, mm, pm, vp, &projX[3], &projY[3], &projZ[3]);


   // Z axis location (left most corner)
   //
   Int_t idxLeft = 0;
   Float_t xt = projX[0];
   for (Int_t i = 1; i < 4; ++i) {
      if (projX[i] < xt) {
         xt  = projX[i];
         idxLeft = i;
      }
   }
   fZAxisTitlePos.Set(cornerX[idxLeft], cornerY[idxLeft], 1.05 * fMaxVal);


   // XY axis location (closest to eye) first
   //
   Float_t zt = 1.f;
   Float_t zMin = 0.f;
   Int_t idxFront = 0;
   for (Int_t i = 0; i < 4; ++i) {
      if (projZ[i] < zt) {
         zt  = projZ[i];
         idxFront = i;
      }
      if (projZ[i] > zMin) zMin = projZ[i];
   }


   Int_t xyIdx = idxFront;
   if (zMin - zt < 1e-2) xyIdx = 0; // avoid flipping in front view


   switch (xyIdx) {
      case 0:
         fXAxisTitlePos.fX = x1;
         fXAxisTitlePos.fY = y0;
         fYAxisTitlePos.fX = x0;
         fYAxisTitlePos.fY = y1;
         break;
      case 1:
         fXAxisTitlePos.fX = x0;
         fXAxisTitlePos.fY = y0;
         fYAxisTitlePos.fX = x1;
         fYAxisTitlePos.fY = y1;
         break;
      case 2:
         fXAxisTitlePos.fX = x0;
         fXAxisTitlePos.fY = y1;
         fYAxisTitlePos.fX = x1;
         fYAxisTitlePos.fY = y0;
         break;
      case 3:
         fXAxisTitlePos.fX = x1;
         fXAxisTitlePos.fY = y1;
         fYAxisTitlePos.fX = x0;
         fYAxisTitlePos.fY = y0;
         break;
   }

   // move title 5% over the axis length
   Float_t off = 0.05;
   Float_t tOffX = (x1-x0) * off; if (fYAxisTitlePos.fX > x0) tOffX = -tOffX;
   Float_t tOffY = (y1-y0) * off; if (fXAxisTitlePos.fY > y0) tOffY = -tOffY;
   fXAxisTitlePos.fX += tOffX;
   fYAxisTitlePos.fY += tOffY;


   // frame box
   //
   if (fM->fBoxMode)
   {
      // get corner closest to eye excluding left corner
      Double_t zm = 1.f;
      Int_t idxDepthT = 0;
      for (Int_t i = 0; i < 4; ++i)
      {
         if (projZ[i] < zm && projZ[i] >= zt && i != idxFront )
         {
            zm  = projZ[i];
            idxDepthT = i;
         }
      }
      if (idxFront == idxLeft)  idxFront =idxDepthT;

      switch (idxFront)
      {
         case 0:
            fBackPlaneXConst[0].Set(x1, y0, 0); fBackPlaneXConst[1].Set(x1, y1, 0);
            fBackPlaneYConst[0].Set(x0, y1, 0); fBackPlaneYConst[1].Set(x1, y1, 0);
            break;
         case 1:
            fBackPlaneXConst[0].Set(x0, y0, 0); fBackPlaneXConst[1].Set(x0, y1, 0);
            fBackPlaneYConst[0].Set(x0, y1, 0); fBackPlaneYConst[1].Set(x1, y1, 0);
            break;
         case 2:
            fBackPlaneXConst[0].Set(x0, y0, 0); fBackPlaneXConst[1].Set(x0, y1, 0);
            fBackPlaneYConst[0].Set(x0, y0, 0); fBackPlaneYConst[1].Set(x1, y0, 0);
            break;
         case 3:
            fBackPlaneXConst[0].Set(x1, y0, 0); fBackPlaneXConst[1].Set(x1, y1, 0);
            fBackPlaneYConst[0].Set(x0, y0, 0); fBackPlaneYConst[1].Set(x1, y0, 0);
            break;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw z-axis and z-box at the appropriate grid corner-point including
/// tick-marks and labels.

void TEveCaloLegoGL::DrawAxis3D(TGLRnrCtx & rnrCtx) const
{
   // set font size first depending on size of projected axis

   TGLMatrix mm;
   GLdouble pm[16];
   glGetDoublev(GL_MODELVIEW_MATRIX, mm.Arr());
   glGetDoublev(GL_PROJECTION_MATRIX, pm);
   Int_t* vp = rnrCtx.RefCamera().RefViewport().CArr();
   GLdouble dn[3];
   GLdouble up[3];
   gluProject(fXAxisTitlePos.fX, fXAxisTitlePos.fY, fXAxisTitlePos.fZ, mm.Arr(), pm, vp, &up[0], &up[1], &up[2]);
   gluProject(fYAxisTitlePos.fX, fYAxisTitlePos.fY, fYAxisTitlePos.fZ, mm.Arr(), pm, vp, &dn[0], &dn[1], &dn[2]);
   Float_t len = TMath::Sqrt((up[0] - dn[0]) * (up[0] - dn[0])
                              + (up[1] - dn[1]) * (up[1] - dn[1])
                              + (up[2] - dn[2]) * (up[2] - dn[2]));
   len = TMath::Min(len, rnrCtx.RefCamera().RefViewport().Diagonal()*0.7f);
   len /= TMath::Sqrt2();

   TGLVertex3 worldRef(fZAxisTitlePos.fX, fZAxisTitlePos.fY, fZAxisTitlePos.fZ);
   fAxisPainter.RefTMOff(0) = rnrCtx.RefCamera().ViewportDeltaToWorld(worldRef, -len, 0,  &mm);
   fAxisPainter.SetLabelPixelFontSize(TMath::Nint(TMath::Max(len*fM->GetData()->GetEtaBins()->GetLabelSize(), 0.02f)));
   fAxisPainter.SetTitlePixelFontSize(TMath::Nint(TMath::Max(len*fM->GetData()->GetEtaBins()->GetLabelSize(), 0.02f)));

   Float_t tickLength = TMath::Max(fM->GetData()->GetEtaBins()->GetTickLength(), 0.02f);
   Float_t labelOffset = TMath::Max(fM->GetData()->GetEtaBins()->GetLabelOffset(), 0.02f);

   // Z axis
   //
   if (fM->fData->Empty() == kFALSE)
   {
      Int_t    ondiv;
      Double_t omin=0, omax=0, bw1;
      THLimitsFinder::Optimize(0, fMaxVal, fM->fNZSteps, omin, omax, ondiv, bw1);
      worldRef.Set(fZAxisTitlePos.fX, fZAxisTitlePos.fY, fZAxisTitlePos.fZ);
      TGLVector3 zto = rnrCtx.RefCamera().ViewportDeltaToWorld(worldRef, 0, fAxisPainter.GetLabelPixelFontSize(),  &mm);
      // check z axis title does not overalp with label
      if ( fZAxisTitlePos.fZ - omax <   fAxisPainter.GetLabelPixelFontSize())
         fZAxisTitlePos.fZ = omax + zto.Z();


      fZAxis->SetAxisColor(fGridColor);
      fZAxis->SetLabelColor(fFontColor);
      fZAxis->SetTitleColor(fFontColor);
      fZAxis->SetNdivisions(fM->fNZSteps*100 + 10);
      fZAxis->SetLimits(0, fMaxVal);
      fZAxis->SetTitle(fM->GetPlotEt() ? "Et[GeV]" : "E[GeV]");

      fAxisPainter.SetTMNDim(1);
      fAxisPainter.RefDir().Set(0., 0., 1.);
      fAxisPainter.SetLabelAlign(TGLFont::kRight, TGLFont::kCenterV);
      glPushMatrix();
      glTranslatef(fZAxisTitlePos.fX, fZAxisTitlePos.fY, 0);

      // tickmark vector = 10 pixels left
      fAxisPainter.RefTitlePos().Set(fAxisPainter.RefTMOff(0).X()*0.05,  fAxisPainter.RefTMOff(0).Y()*0.05, fZAxisTitlePos.fZ);
      fZAxis->SetLabelOffset(labelOffset);
      fZAxis->SetTickLength(tickLength);
      fAxisPainter.PaintAxis(rnrCtx, fZAxis);
      glPopMatrix();

      // draw box frame
      //
      if (fM->fBoxMode) {

         glPushAttrib(GL_ENABLE_BIT | GL_LINE_BIT);

         // box verticals
         TGLUtil::LineWidth(1);
         glBegin(GL_LINES);
         TGLUtil::Color(fGridColor);

         glVertex3f(fBackPlaneXConst[0].fX   ,fBackPlaneXConst[0].fY   ,0);
         glVertex3f(fBackPlaneXConst[0].fX   ,fBackPlaneXConst[0].fY   ,fMaxVal);
         glVertex3f(fBackPlaneXConst[1].fX   ,fBackPlaneXConst[1].fY   ,0);
         glVertex3f(fBackPlaneXConst[1].fX   ,fBackPlaneXConst[1].fY   ,fMaxVal);


         glVertex3f(fBackPlaneYConst[0].fX   ,fBackPlaneYConst[0].fY   ,0);
         glVertex3f(fBackPlaneYConst[0].fX   ,fBackPlaneYConst[0].fY   ,fMaxVal);
         glVertex3f(fBackPlaneYConst[1].fX   ,fBackPlaneYConst[1].fY   ,0);
         glVertex3f(fBackPlaneYConst[1].fX   ,fBackPlaneYConst[1].fY   ,fMaxVal);

         // box top
         glVertex3f(fBackPlaneXConst[0].fX   ,fBackPlaneXConst[0].fY   ,fMaxVal);
         glVertex3f(fBackPlaneXConst[1].fX   ,fBackPlaneXConst[1].fY   ,fMaxVal);
         glVertex3f(fBackPlaneYConst[0].fX   ,fBackPlaneYConst[0].fY   ,fMaxVal);
         glVertex3f(fBackPlaneYConst[1].fX   ,fBackPlaneYConst[1].fY   ,fMaxVal);

         glEnd();

         // box horizontals stippled
         glEnable(GL_LINE_STIPPLE);
         glLineStipple(1, 0x5555);
         glBegin(GL_LINES);
         Float_t hz  = bw1;
         for (Int_t i = 1; i <= ondiv; ++i, hz += bw1) {
            glVertex3f(fBackPlaneXConst[0].fX   ,fBackPlaneXConst[0].fY   ,hz);
            glVertex3f(fBackPlaneXConst[1].fX   ,fBackPlaneXConst[1].fY   ,hz);
            glVertex3f(fBackPlaneYConst[0].fX   ,fBackPlaneYConst[0].fY   ,hz);
            glVertex3f(fBackPlaneYConst[1].fX   ,fBackPlaneYConst[1].fY   ,hz);
         }
         glEnd();

         glPopAttrib();
      }
   }

   // XY Axis
   //

   Float_t yOff  = fM->GetPhiRng();
   if (fXAxisTitlePos.fY < fM->GetPhiMax()) yOff = -yOff;

   Float_t xOff  = fM->GetEtaRng();
   if (fYAxisTitlePos.fX < fM->GetEtaMax()) xOff = -xOff;

   TAxis ax;
   ax.SetAxisColor(fGridColor);
   ax.SetLabelColor(fFontColor);
   ax.SetTitleColor(fFontColor);
   ax.SetTitleFont(fM->GetData()->GetEtaBins()->GetTitleFont());
   ax.SetLabelOffset(labelOffset);
   ax.SetTickLength(tickLength);
   fAxisPainter.SetTMNDim(2);
   fAxisPainter.RefTMOff(1).Set(0, 0, -fMaxVal);
   fAxisPainter.SetLabelAlign(TGLFont::kCenterH, TGLFont::kBottom);

   // eta
   glPushMatrix();
   fAxisPainter.RefDir().Set(1, 0, 0);
   fAxisPainter.RefTMOff(0).Set(0, yOff, 0);
   glTranslatef(0, fXAxisTitlePos.fY, 0);

   ax.SetNdivisions(fM->GetData()->GetEtaBins()->GetNdivisions());
   ax.SetLimits(fM->GetEtaMin(), fM->GetEtaMax());
   ax.SetTitle(fM->GetData()->GetEtaBins()->GetTitle());
   fAxisPainter.RefTitlePos().Set(fXAxisTitlePos.fX, yOff*1.5*ax.GetTickLength(), -fMaxVal*ax.GetTickLength());
   fAxisPainter.PaintAxis(rnrCtx, &ax);
   glPopMatrix();

   // phi
   fAxisPainter.RefDir().Set(0, 1, 0);
   fAxisPainter.RefTMOff(0).Set(xOff, 0, 0);
   ax.SetNdivisions(fM->GetData()->GetPhiBins()->GetNdivisions());
   ax.SetLimits(fM->GetPhiMin(), fM->GetPhiMax());
   ax.SetTitle(fM->GetData()->GetPhiBins()->GetTitle());
   glPushMatrix();
   glTranslatef(fYAxisTitlePos.fX, 0, 0);
   fAxisPainter.RefTitlePos().Set( xOff*1.5*ax.GetTickLength(), fYAxisTitlePos.fY,  -fMaxVal*ax.GetTickLength());
   fAxisPainter.PaintAxis(rnrCtx, &ax);
   glPopMatrix();

} // DrawAxis3D

////////////////////////////////////////////////////////////////////////////////
/// Get scale for matrix.

void  TEveCaloLegoGL::GetScaleForMatrix(Float_t& sx, Float_t& sy, Float_t& sz) const
{
   Double_t em, eM, pm, pM;
   fM->fData->GetEtaLimits(em, eM);
   fM->fData->GetPhiLimits(pm, pM);
   Double_t unit = ((eM - em) < (pM - pm)) ? (eM - em) : (pM - pm);
   sx = (eM - em) / (fM->GetEtaRng() * unit);
   sy = (pM - pm) / (fM->GetPhiRng() * unit);

   sz = 1;
   if (fM->fScaleAbs)
   {
      sz = fM->GetMaxTowerH() / fM->fMaxValAbs;
   }
   else if (!fM->fData->Empty())
   {
      sz = fM->GetMaxTowerH() / fMaxVal;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw XY axis.

void TEveCaloLegoGL::DrawAxis2D(TGLRnrCtx & rnrCtx) const
{
   if (fM->GetData()->Empty())
      fAxisPainter.SetTMNDim(1);

   TGLCamera& cam  = rnrCtx.RefCamera();

   TAxis ax;
   ax.SetAxisColor(fGridColor);
   ax.SetLabelColor(fFontColor);
   ax.SetTitleColor(fFontColor);
   ax.SetTitleFont(fM->GetData()->GetEtaBins()->GetTitleFont());
   ax.SetTitleSize(TMath::Max(fM->GetData()->GetEtaBins()->GetTitleSize(), 0.02f));
   ax.SetLabelOffset(TMath::Max(fM->GetData()->GetEtaBins()->GetLabelOffset(), 0.02f));
   ax.SetTickLength(TMath::Max(fM->GetData()->GetEtaBins()->GetTickLength(), 0.05f));

   // set fonts
   fAxisPainter.SetAttAxis(&ax);

   // get projected length of diagonal to determine
   TGLMatrix mm;
   GLdouble pm[16];
   GLint    vp[4];
   glGetDoublev(GL_MODELVIEW_MATRIX, mm.Arr());
   glGetDoublev(GL_PROJECTION_MATRIX, pm);
   glGetIntegerv(GL_VIEWPORT, vp);

   GLdouble dn[3];
   GLdouble up[3];
   gluProject(fM->GetEtaMin(), fM->GetPhiMin(), 0, mm.Arr(), pm, vp, &dn[0], &dn[1], &dn[2]);
   gluProject(fM->GetEtaMax(), fM->GetPhiMax(), 0, mm.Arr(), pm, vp, &up[0], &up[1], &up[2]);
   Double_t len = TMath::Sqrt((up[0] - dn[0]) * (up[0] - dn[0])
                              + (up[1] - dn[1]) * (up[1] - dn[1])
                              + (up[2] - dn[2]) * (up[2] - dn[2]));

   // lock upper limit to of relative font size relative to viewport diagonal
   Double_t vpLimit = cam.RefViewport().Diagonal()*0.5/TMath::Sqrt2();
   len = TMath::Min(len, vpLimit);

   // eta
   fAxisPainter.SetLabelPixelFontSize(TMath::Nint(len*fM->GetData()->GetEtaBins()->GetLabelSize()));
   fAxisPainter.SetTitlePixelFontSize(TMath::Nint(len*fM->GetData()->GetEtaBins()->GetTitleSize()));
   ax.SetNdivisions(fM->GetData()->GetEtaBins()->GetNdivisions());
   ax.SetLimits(fM->GetEtaMin(), fM->GetEtaMax());
   ax.SetTitle(fM->GetData()->GetEtaBins()->GetTitle());
   fAxisPainter.RefTitlePos().Set(fM->GetEtaMax(), -fM->GetPhiRng()*(ax.GetTickLength()+ ax.GetLabelOffset()), 0 );
   fAxisPainter.RefDir().Set(1, 0, 0);

   Float_t tmOffFrustX = cam.FrustumPlane(TGLCamera::kRight).D() + cam.FrustumPlane(TGLCamera::kLeft).D();
   fAxisPainter.RefTMOff(0).Set(0,  -TMath::Min(fM->GetPhiRng(), tmOffFrustX), 0);
   fAxisPainter.SetLabelAlign(TGLFont::kCenterH, TGLFont::kBottom);

   glPushMatrix();
   glTranslatef(0, fM->GetPhiMin(), 0);
   fAxisPainter.PaintAxis(rnrCtx, &ax);
   glPopMatrix();

   // phi
   ax.SetNdivisions(fM->GetData()->GetPhiBins()->GetNdivisions());
   ax.SetLimits(fM->GetPhiMin(), fM->GetPhiMax());
   ax.SetTitle(fM->GetData()->GetPhiBins()->GetTitle());
   fAxisPainter.RefTitlePos().Set(-fM->GetEtaRng()*(ax.GetTickLength()+ ax.GetLabelOffset()), fM->GetPhiMax(), 0);
   fAxisPainter.RefDir().Set(0, 1, 0);
   Float_t tmOffFrustY = cam.FrustumPlane(TGLCamera::kTop).D() + cam.FrustumPlane(TGLCamera::kBottom).D();
   fAxisPainter.RefTMOff(0).Set(-TMath::Min(fM->GetEtaRng(), tmOffFrustY), 0, 0);
   fAxisPainter.SetLabelAlign(TGLFont::kRight, TGLFont::kCenterV);

   glPushMatrix();
   glTranslatef(fM->GetEtaMin(), 0, 0);
   fAxisPainter.PaintAxis(rnrCtx, &ax);
   glPopMatrix();

   fAxisPainter.SetTMNDim(2);
}

////////////////////////////////////////////////////////////////////////////////
/// Calculate view-dependent grid density.

Int_t TEveCaloLegoGL::GetGridStep(TGLRnrCtx &rnrCtx) const
{
   TGLCamera &camera = rnrCtx.RefCamera();
   Float_t l = -camera.FrustumPlane(TGLCamera::kLeft).D();
   Float_t r =  camera.FrustumPlane(TGLCamera::kRight).D();
   Float_t t =  camera.FrustumPlane(TGLCamera::kTop).D();
   Float_t b = -camera.FrustumPlane(TGLCamera::kBottom).D();
   Float_t frustD    = TMath::Hypot(r-l, t-b);

   GLint   vp[4]; glGetIntegerv(GL_VIEWPORT, vp);
   Float_t viewportD = TMath::Sqrt((vp[1] - vp[0]) * (vp[1] - vp[0]) + (vp[3] - vp[1]) * (vp[3] - vp[1]));
   Float_t deltaToViewport = viewportD/frustD;

   // average bin width
   GLdouble em, eM, pm, pM;
   fM->GetData()->GetEtaLimits(pm, pM);
   fM->GetData()->GetPhiLimits(em, eM);
   Int_t i0 = fM->fData->GetEtaBins()->FindBin(fM->GetEtaMin());
   Int_t i1 = fM->fData->GetEtaBins()->FindBin(fM->GetEtaMax());
   Int_t j0 = fM->fData->GetPhiBins()->FindBin(fM->GetPhiMin());
   Int_t j1 = fM->fData->GetPhiBins()->FindBin(fM->GetPhiMax());

   Float_t averageBinWidth = TMath::Hypot(eM - em, pM - pm)/TMath::Sqrt((i0 - i1) * (i0 - i1) + (j0 - j1) * (j0 - j1));
   Float_t ppb = deltaToViewport*averageBinWidth;

   Int_t ngroup = 1;
   if (fM->fAutoRebin && fM->fPixelsPerBin > ppb)
   {
      // limit rebin realtive to number of axis bins
      Int_t maxGroup = TMath::Min(fM->fData->GetEtaBins()->GetNbins(), fM->fData->GetPhiBins()->GetNbins())/4;
      if (maxGroup > 1) {
         ngroup = TMath::Nint(fM->fPixelsPerBin*0.5/ppb); // symetrical rebin factor 2
         if (ngroup  > maxGroup) ngroup = maxGroup;
      }
   }
   fCurrentPixelsPerBin = TMath::Nint(ppb);

   return ngroup;
}

////////////////////////////////////////////////////////////////////////////////
/// Rebin eta, phi axis.

void TEveCaloLegoGL::RebinAxis(TAxis *orig, TAxis *curr) const
{
   Double_t center = 0.5 * (orig->GetXmin() + orig->GetXmax());
   Int_t    idx0   = orig->FindBin(center);
   Double_t bc     = orig->GetBinCenter(idx0);
   if (bc > center) --idx0;

   Int_t nbR = TMath::FloorNint(idx0/fBinStep) + TMath::FloorNint((orig->GetNbins() - idx0)/fBinStep);
   Int_t off = idx0 - TMath::FloorNint(idx0/fBinStep)*fBinStep;
   std::vector<Double_t> bins(nbR + 1);
   for (Int_t i = 0; i <= nbR; ++i)
   {
      bins[i] = orig->GetBinUpEdge(off + i*fBinStep);
   }
   curr->Set(nbR, &bins[0]);
}

////////////////////////////////////////////////////////////////////////////////
/// Draw basic histogram components: x-y grid

void TEveCaloLegoGL::DrawHistBase(TGLRnrCtx &rnrCtx) const
{
   Float_t eta0 = fM->fEtaMin;
   Float_t eta1 = fM->fEtaMax;
   Float_t phi0 = fM->GetPhiMin();
   Float_t phi1 = fM->GetPhiMax();

   // XY grid
   //
   TGLUtil::Color(fGridColor);
   TGLUtil::LineWidth(1);
   glBegin(GL_LINES);
   glVertex2f(eta0, phi0);
   glVertex2f(eta0, phi1);
   glVertex2f(eta1, phi0);
   glVertex2f(eta1, phi1);

   glVertex2f(eta0, phi0);
   glVertex2f(eta1, phi0);
   glVertex2f(eta0, phi1);
   glVertex2f(eta1, phi1);

   // eta grid
   Float_t val;
   Int_t neb = fEtaAxis->GetNbins();
   for (Int_t i = 0; i<= neb; i++)
   {
      val = fEtaAxis->GetBinUpEdge(i);
      if (val > eta0 && val < eta1 )
      {
         glVertex2f(val, phi0);
         glVertex2f(val, phi1);
      }
   }

   // phi grid
   Int_t npb = fPhiAxis->GetNbins();
   for (Int_t i = 1; i <= npb; i++) {
       val = fPhiAxis->GetBinUpEdge(i);
      if (val > phi0 && val < phi1)
      {
         glVertex2f(eta0, val);
         glVertex2f(eta1, val);
      }
   }

   glEnd();

   // XYZ axes
   //
   glPushAttrib(GL_ENABLE_BIT | GL_LINE_BIT | GL_POLYGON_BIT);
   TGLUtil::LineWidth(2);
   if (fCells3D)
   {
      SetAxis3DTitlePos(rnrCtx, eta0, eta1, phi0, phi1);
      DrawAxis3D(rnrCtx);
   }
   else
   {
      DrawAxis2D(rnrCtx);
   }
   glPopAttrib();
}

////////////////////////////////////////////////////////////////////////////////
/// Render the calo lego-plot with OpenGL.

void TEveCaloLegoGL::DrawCells3D(TGLRnrCtx & rnrCtx) const
{
   // quads
   {
      for (SliceDLMap_i i = fDLMap.begin(); i != fDLMap.end(); ++i) {
         TGLUtil::ColorTransparency(fM->GetDataSliceColor(i->first), fM->GetData()->GetSliceTransparency(i->first));
         glLoadName(i->first);
         glPushName(0);
         glCallList(i->second);
         glPopName();
      }
   }
   // outlines
   {
      if (rnrCtx.SceneStyle() == TGLRnrCtx::kFill) {
         glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
         glDisable(GL_POLYGON_OFFSET_FILL);
         TGLUtil::Color(1);
         for (SliceDLMap_i i = fDLMap.begin(); i != fDLMap.end(); ++i)
            glCallList(i->second);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Prepare cells 2D data non-rebinned for drawing.

void TEveCaloLegoGL::PrepareCell2DData(TEveCaloData::vCellId_t& cellList, vCell2D_t& cells2D) const
{
   Int_t   max_energy_slice, cellID=0;
   Float_t sum, max_energy;

   TEveCaloData::vCellId_t::iterator currentCell = cellList.begin();
   TEveCaloData::vCellId_t::iterator nextCell    = currentCell;
   ++nextCell;

   while (true)
   {
      TEveCaloData::CellData_t currentCellData;
      TEveCaloData::CellData_t nextCellData;

      fM->fData->GetCellData(*currentCell, currentCellData);
      sum = max_energy = currentCellData.Value(fM->fPlotEt);
      max_energy_slice = currentCell->fSlice;
      while (nextCell != cellList.end() && currentCell->fTower == nextCell->fTower)
      {
         fM->fData->GetCellData(*nextCell, nextCellData);
         Float_t energy = nextCellData.Value(fM->fPlotEt);
         sum += energy;
         if (energy > max_energy)
         {
            max_energy       = energy;
            max_energy_slice = nextCell->fSlice;
         }
         ++nextCell;
         ++cellID;
      }

      WrapTwoPi(currentCellData.fPhiMin, currentCellData.fPhiMax);
      cells2D.push_back(Cell2D_t(cellID, sum, max_energy_slice));
      cells2D.back().SetGeom(currentCellData.fEtaMin, currentCellData.fEtaMax,
                             currentCellData.fPhiMin, currentCellData.fPhiMax);

      if (nextCell == cellList.end())
         break;

      currentCell = nextCell;
      ++nextCell;
      ++cellID;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Prepare cells 2D rebinned data for drawing.

void TEveCaloLegoGL::PrepareCell2DDataRebin(TEveCaloData::RebinData_t& rebinData, vCell2D_t& cells2D) const
{
   const Int_t nEta = fEtaAxis->GetNbins();
   const Int_t nPhi = fPhiAxis->GetNbins();
   std::vector<Float_t> vec;
   vec.assign((nEta + 2)*(nPhi + 2), 0.f);
   std::vector<Float_t> max_e;
   std::vector<Int_t>   max_e_slice;
   max_e.assign((nEta + 2) * (nPhi + 2), 0.f);
   max_e_slice.assign((nEta + 2) * (nPhi + 2), -1);

   for (UInt_t bin = 0; bin < rebinData.fBinData.size(); ++bin) {
      Float_t ssum = 0;
      if (rebinData.fBinData[bin] != -1) {
         Float_t *val = rebinData.GetSliceVals(bin);
         for (Int_t s = 0; s < rebinData.fNSlices; ++s) {
            ssum += val[s];
            if (val[s] > max_e[bin]) {
               max_e[bin]       = val[s];
               max_e_slice[bin] = s;
            }
         }
      }
      vec[bin] = ssum;
   }

   // smallest threshold
   Float_t threshold = fM->GetDataSliceThreshold(0);
   for (Int_t s = 1; s < fM->fData->GetNSlices(); ++s) {
      if (threshold > fM->GetDataSliceThreshold(s))
         threshold = fM->GetDataSliceThreshold(s);
   }

   // write cells
   for (Int_t i = 1; i <= fEtaAxis->GetNbins(); ++i) {
      for (Int_t j = 1; j <= fPhiAxis->GetNbins(); ++j) {
         const Int_t bin = j * (nEta + 2) + i;
         if (vec[bin] > threshold && rebinData.fBinData[bin] != -1) {
            cells2D.push_back(Cell2D_t(bin, vec[bin], max_e_slice[bin]));
            cells2D.back().SetGeom(fEtaAxis->GetBinLowEdge(i), fEtaAxis->GetBinUpEdge(i),
                                   fPhiAxis->GetBinLowEdge(j), fPhiAxis->GetBinUpEdge(j));
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw cells in top view.

void TEveCaloLegoGL::DrawCells2D(TGLRnrCtx &rnrCtx, vCell2D_t& cells2D) const
{
   Float_t bws    = -1; //smallest bin
   Float_t logMax = -1;

   Float_t baseOffset = fM->GetFixedHeightValIn2DMode()*fMaxVal;

   if (fM->f2DMode == TEveCaloLego::kValColor)
   {
      fM->AssertPalette();
      UChar_t col[4];

      for (vCell2D_i i = cells2D.begin(); i != cells2D.end(); ++i)
      {
         if (rnrCtx.SecSelection()) glLoadName(i->fId);
         glBegin(GL_POLYGON);
         fM->fPalette->ColorFromValue(TMath::FloorNint(i->fSumVal), col);
         col[3] = fM->GetData()->GetSliceTransparency(i->fMaxSlice);
         TGLUtil::Color4ubv(col);
         Float_t z = fM->GetHasFixedHeightIn2DMode() ? baseOffset : i->fSumVal;
         glVertex3f(i->fX0, i->fY0, z);
         glVertex3f(i->fX1, i->fY0, z);
         glVertex3f(i->fX1, i->fY1, z);
         glVertex3f(i->fX0, i->fY1, z);
         glEnd();
      }
   }
   else
   {
      Float_t maxv = 0;
      bws = 1e5;
      for (vCell2D_i i = fCells2D.begin(); i != fCells2D.end(); ++i)
      {
         if (i->MinSize() < bws)   bws  = i->MinSize();
         if (i->fSumVal   > maxv)  maxv = i->fSumVal;
      }
      bws   *= 0.5f;
      logMax = TMath::Log10(maxv + 1);
      fValToPixel = bws / logMax;

      if (rnrCtx.SecSelection())
      {
         // Special draw for name stack.

         for (vCell2D_i i = cells2D.begin(); i != cells2D.end(); ++i)
         {
            glLoadName(i->fMaxSlice);
            glPushName(i->fId);

            glBegin(GL_QUADS);
            Float_t z = fM->GetHasFixedHeightIn2DMode() ? baseOffset : i->fSumVal;
            glVertex3f(i->fX0, i->fY0, z);
            glVertex3f(i->fX1, i->fY0, z);
            glVertex3f(i->fX1, i->fY1, z);
            glVertex3f(i->fX0, i->fY1, z);
            glEnd();

            glPopName();
         }
      }
      else
      {
         // Optimised draw without name stack.

         if ( ! rnrCtx.Highlight())
         {
            glBegin(GL_POINTS);
            for (vCell2D_i i = cells2D.begin(); i != cells2D.end(); ++i)
            {
               TGLUtil::ColorTransparency(fM->fData->GetSliceColor(i->fMaxSlice), fM->fData->GetSliceTransparency(i->fMaxSlice));
               Float_t z = fM->GetHasFixedHeightIn2DMode() ? baseOffset : i->fSumVal;
               glVertex3f(i->X(), i->Y() , z);
            }
            glEnd();
         }

         glBegin(GL_QUADS);
         for (vCell2D_i i = cells2D.begin(); i != cells2D.end(); ++i)
         {
            TGLUtil::ColorTransparency(fM->fData->GetSliceColor(i->fMaxSlice), fM->fData->GetSliceTransparency(i->fMaxSlice));
            Float_t bw = fValToPixel*TMath::Log10(i->fSumVal+1);
            Float_t x = i->X();
            Float_t y = i->Y();
            Float_t z = fM->GetHasFixedHeightIn2DMode() ? baseOffset : i->fSumVal;
            glVertex3f(x - bw, y - bw, z);
            glVertex3f(x + bw, y - bw, z);
            glVertex3f(x + bw, y + bw, z);
            glVertex3f(x - bw, y + bw, z);
         }
         glEnd();

         if (fM->f2DMode == TEveCaloLego::kValSizeOutline)
         {
            glPushAttrib(GL_ENABLE_BIT | GL_POLYGON_BIT);
            Float_t z    = 0;
            Float_t zOff = fMaxVal*0.001 ; // avoid polygon stipling
            glBegin(GL_QUADS);
            for ( vCell2D_i i = cells2D.begin(); i != cells2D.end(); ++i) {
               Char_t transp = TMath::Min(100, 80 + fM->fData->GetSliceTransparency(i->fMaxSlice) / 5);
               TGLUtil::ColorTransparency(fM->fData->GetSliceColor(i->fMaxSlice), transp);
               z = fM->GetHasFixedHeightIn2DMode() ? baseOffset : i->fSumVal;
               z -=  zOff;
               glVertex3f(i->fX0, i->fY0, z);
               glVertex3f(i->fX1, i->fY0, z);
               glVertex3f(i->fX1, i->fY1, z);
               glVertex3f(i->fX0, i->fY1, z);
            }
            glEnd();

            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            glBegin(GL_QUADS);
            for ( vCell2D_i i = cells2D.begin(); i != cells2D.end(); ++i) {
               TGLUtil::ColorTransparency(fM->fData->GetSliceColor(i->fMaxSlice), 60);
               z = fM->GetHasFixedHeightIn2DMode() ? baseOffset : i->fSumVal;
               z +=  zOff;
               glVertex3f(i->fX0, i->fY0, z);
               glVertex3f(i->fX1, i->fY0, z);
               glVertex3f(i->fX1, i->fY1, z);
               glVertex3f(i->fX0, i->fY1, z);
            }
            glEnd();
            glPopAttrib();
         }
      }
   }

   // text
   if (fCurrentPixelsPerBin > fM->fDrawNumberCellPixels &&
       ! rnrCtx.Selection() && ! rnrCtx.Highlight())
   {
      TGLUtil::Color(rnrCtx.ColorSet().Markup().GetColorIndex());
      TGLFont font;
      rnrCtx.RegisterFontNoScale(fM->fCellPixelFontSize, "arial", TGLFont::kPixmap, font);
      const char* txt;
      for (vCell2D_i i = cells2D.begin(); i != cells2D.end(); ++i) {

         Float_t val = i->fSumVal;
         if (val > 10)
            txt = Form("%d", TMath::Nint(val));
         else if (val > 1 )
            txt = Form("%.1f", val);
         else if (val > 0.01 )
            txt = Form("%.2f", 0.01*TMath::Nint(val*100));
         else
            txt = Form("~1e%d", TMath::Nint(TMath::Log10(val)));

         font.Render(txt, i->X(), i->Y(), val*1.2, TGLFont::kCenterH, TGLFont::kCenterV);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw highligted cells.

void TEveCaloLegoGL::DrawHighlight(TGLRnrCtx& rnrCtx, const TGLPhysicalShape* /*pshp*/, Int_t /*lvl*/) const
{
   if (fM->fData->GetCellsSelected().empty() && fM->fData->GetCellsHighlighted().empty())
   {
      return;
   }

   // modelview matrix
   glPushMatrix();
   Float_t sx, sy, sz;
   GetScaleForMatrix(sx, sy, sz);
   glScalef(sx, sy, sz);
   glTranslatef(-fM->GetEta(), -fM->fPhi, 0);

   if (fCells3D)
   {
      glPushAttrib(GL_ENABLE_BIT | GL_LINE_BIT | GL_POLYGON_BIT);
      glDisable(GL_LIGHTING);
      glDisable(GL_CULL_FACE);
      glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
      TGLUtil::LineWidth(2);
   }

   TGLUtil::LockColor();
   if (!fM->fData->GetCellsHighlighted().empty())
   {
      glColor4ubv(rnrCtx.ColorSet().Selection(3).CArr());
      DrawSelectedCells(rnrCtx, fM->fData->GetCellsHighlighted());
   }
   if (!fM->fData->GetCellsSelected().empty())
   {
      glColor4ubv(rnrCtx.ColorSet().Selection(1).CArr());
      DrawSelectedCells(rnrCtx, fM->fData->GetCellsSelected());
   }
   TGLUtil::UnlockColor();

   if (fCells3D)
   {
      glPopAttrib();
   }

   glPopMatrix();
}


////////////////////////////////////////////////////////////////////////////////
/// Draw selected cells in highlight mode.

void TEveCaloLegoGL::DrawSelectedCells(TGLRnrCtx & rnrCtx, TEveCaloData::vCellId_t cellsSelectedInput) const
{
   // check eta&phi range of selected cells
   TEveCaloData::vCellId_t  cellsSelected;
   TEveCaloData::CellData_t cellData;
   for (TEveCaloData::vCellId_i i = cellsSelectedInput.begin(); i != cellsSelectedInput.end(); ++i)
   {
      fM->fData->GetCellData((*i), cellData);
      if (fM->CellInEtaPhiRng(cellData))
         cellsSelected.push_back(*i);
   }

   // prepare rebin for 2D or 3D if necessary
   TEveCaloData::RebinData_t rebinDataSelected;
   if (fBinStep > 1)
   {
      fM->fData->Rebin(fEtaAxis, fPhiAxis, cellsSelected, fM->fPlotEt, rebinDataSelected);
      if (fM->fNormalizeRebin) {
         Float_t scale = 1.f / (fBinStep * fBinStep);
         for (std::vector<Float_t>::iterator it = rebinDataSelected.fSliceData.begin();
              it != rebinDataSelected.fSliceData.end(); ++it)
            (*it) *= scale;
      }
   }

   if (fCells3D)
   {
      Float_t offset =  0;
      if (fBinStep == 1)
      {
         for (TEveCaloData::vCellId_i j = cellsSelected.begin(); j != cellsSelected.end(); ++j)
         {
            offset = 0;
            {
               Int_t   orig_slice = j->fSlice;
               for (Int_t s = 0; s < orig_slice; ++s)
               {
                  j->fSlice = s;
                  fM->fData->GetCellData(*j, cellData);
                  offset += cellData.Value(fM->fPlotEt);
               }
               j->fSlice = orig_slice;
            }
            fM->fData->GetCellData(*j, cellData);
            WrapTwoPi(cellData.fPhiMin, cellData.fPhiMax);
            MakeQuad(cellData.EtaMin(), cellData.PhiMin(), offset,
                     cellData.EtaDelta(), cellData.PhiDelta(), cellData.Value(fM->fPlotEt));
         }
      }
      else
      {
         Float_t *vals;
         Float_t *valsRef;
         Float_t  y0, y1;
         Int_t    nSlices = fM->fData->GetNSlices();
         for (Int_t i = 1; i <= fEtaAxis->GetNbins(); ++i)
         {
            for (Int_t j = 1; j <= fPhiAxis->GetNbins(); ++j)
            {
               const Int_t bin = (i)+(j)*(fEtaAxis->GetNbins()+2);
               if (rebinDataSelected.fBinData[bin] !=-1)
               {
                  offset  = 0;
                  vals    = rebinDataSelected.GetSliceVals(bin);
                  valsRef = fRebinData.GetSliceVals(bin);
                  for (Int_t s = 0; s < nSlices; ++s)
                  {
                     if (vals[s] > 0)
                     {
                        y0 = fPhiAxis->GetBinLowEdge(j);
                        y1 = fPhiAxis->GetBinUpEdge(j);
                        WrapTwoPi(y0, y1);
                        MakeQuad(fEtaAxis->GetBinLowEdge(i), y0, offset,
                                 fEtaAxis->GetBinWidth(i), y1-y0, vals[s]);
                     }
                     offset += valsRef[s];
                  }
               }
            }
         }
      }
   }
   else
   {
      vCell2D_t cells2DSelected;
      if (fBinStep == 1)
      {
         // but is confusing since top view does not draw all slices at same time
         TEveCaloData::vCellId_i j    = cellsSelectedInput.begin();
         TEveCaloData::vCellId_i jEnd = cellsSelectedInput.end();
         std::set<Int_t> towers;
         while (j != jEnd)
         {
            towers.insert(j->fTower);
            ++j;
         }
         for (vCell2D_i i = fCells2D.begin(); i != fCells2D.end(); ++i)
         {
            TEveCaloData::CellId_t cell = fM->fCellList[i->fId];
            // std::set<Int_t>::iterator ti = towers.find(cell.fTower);
            if (towers.find(cell.fTower) != towers.end())
            {
               cells2DSelected.push_back(*i);
            }
         }
      }
      else
      {
         PrepareCell2DDataRebin(rebinDataSelected, cells2DSelected);
      }
      DrawCells2D(rnrCtx, cells2DSelected);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw the object.

void TEveCaloLegoGL::DirectDraw(TGLRnrCtx & rnrCtx) const
{
   if (! fM->fData || ! fM->fData->GetEtaBins() || ! fM->fData->GetPhiBins())
      return;

   // projection type
   if (fM->fProjection == TEveCaloLego::kAuto)
      fCells3D = (!(rnrCtx.RefCamera().IsOrthographic() && rnrCtx.RefCamera().GetCamBase().GetBaseVec(1).Z()));
   else if (fM->fProjection == TEveCaloLego::k2D)
      fCells3D = kFALSE;
   else if (fM->fProjection == TEveCaloLego::k3D)
      fCells3D = kTRUE;

   // rebin axsis , check limits, fix TwoPi cycling
   Int_t new_bin_step = GetGridStep(rnrCtx);

   // rebin data
   if (fM->AssertCellIdCache() || fBinStep != new_bin_step)
   {
      fBinStep = new_bin_step;
      fDLCacheOK   = kFALSE;
      fRebinData.Clear();

      RebinAxis(fM->fData->GetEtaBins(), fEtaAxis);
      RebinAxis(fM->fData->GetPhiBins(), fPhiAxis);

      if (fBinStep > 1)
      {
         fM->fData->Rebin(fEtaAxis, fPhiAxis, fM->fCellList, fM->fPlotEt, fRebinData);

         fMaxVal = 0;
         for (UInt_t i = 0; i < fRebinData.fSliceData.size(); i += fRebinData.fNSlices)
         {
            Double_t sum = 0;
            for (Int_t s = 0; s < fRebinData.fNSlices; s++)
            {
               sum += fRebinData.fSliceData[i+s];
            }
            if (sum > fMaxVal) fMaxVal = sum;
         }

         if (fM->fNormalizeRebin)
         {
            Float_t scale = 1.f / (fBinStep * fBinStep);
            for (std::vector<Float_t>::iterator it = fRebinData.fSliceData.begin(); it != fRebinData.fSliceData.end();
                 ++it) {
               (*it) *= scale;
            }
            fMaxVal *= scale;
         }
      }
      else
      {
         fMaxVal = fM->GetMaxVal();
      }
   }

   // modelview matrix
   glPushMatrix();
   Float_t sx, sy, sz;
   GetScaleForMatrix(sx, sy, sz);
   glScalef(sx, sy, sz);
   glTranslatef(-fM->GetEta(), -fM->fPhi, 0);

   fFontColor = fM->fFontColor;
   fGridColor = fM->fGridColor;
   if (fGridColor < 0 || fFontColor < 0)
   {
      TColor* c1 = gROOT->GetColor(rnrCtx.ColorSet().Markup().GetColorIndex());
      TColor* c2 = gROOT->GetColor(rnrCtx.ColorSet().Background().GetColorIndex());
      Float_t f1, f2;
      if (fFontColor < 0) {
         f1 = 0.8; f2 = 0.2;
         fFontColor = TColor::GetColor(c1->GetRed()  *f1  + c2->GetRed()  *f2,
                                       c1->GetGreen()*f1  + c2->GetGreen()*f2,
                                       c1->GetBlue() *f1  + c2->GetBlue() *f2);
      }
      if (fGridColor < 0) {
         f1 = 0.3; f2 = 0.3;
         fGridColor = TColor::GetColor(c1->GetRed()  *f1  + c2->GetRed()  *f2,
                                       c1->GetGreen()*f1  + c2->GetGreen()*f2,
                                       c1->GetBlue() *f1  + c2->GetBlue() *f2);
      }
   }

   glPushAttrib(GL_ENABLE_BIT | GL_LINE_BIT | GL_POLYGON_BIT);
   TGLUtil::LineWidth(1);
   glEnable(GL_BLEND);
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
   if (!fM->fData->Empty())
   {
      glPushName(0);
      if (fCells3D)
      {
         if (fDLCacheOK == kFALSE)
         {
            if (fBinStep == 1)
               Make3DDisplayList(fM->fCellList, fDLMap, kTRUE);
            else
               Make3DDisplayListRebin(fRebinData, fDLMap, kTRUE);
            fDLCacheOK = kTRUE;
         }
         glEnable(GL_NORMALIZE);
         glEnable(GL_POLYGON_OFFSET_FILL);
         glPolygonOffset(0.8, 1);

         DrawCells3D(rnrCtx);
      }
      else
      {
         glDisable(GL_LIGHTING);

         fCells2D.clear();
         if (fBinStep == 1)
            PrepareCell2DData(fM->fCellList, fCells2D);
         else
            PrepareCell2DDataRebin(fRebinData, fCells2D);

         DrawCells2D(rnrCtx, fCells2D);
      }
      glPopName();
   }
   glPopAttrib();

   // draw histogram base
   if (rnrCtx.Selection() == kFALSE && rnrCtx.IsDrawPassFilled())
   {
      glPushAttrib(GL_ENABLE_BIT | GL_LINE_BIT | GL_POLYGON_BIT);
      glDisable(GL_LIGHTING);
      DrawHistBase(rnrCtx);
      if (fM->fDrawHPlane) {
         glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
         glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
         glDisable(GL_CULL_FACE);
         TGLUtil::ColorTransparency(fM->fPlaneColor, fM->fPlaneTransparency);
         Float_t zhp = fM->fHPlaneVal * fMaxVal;
         glBegin(GL_POLYGON);
         glVertex3f(fM->fEtaMin, fM->GetPhiMin(), zhp);
         glVertex3f(fM->fEtaMax, fM->GetPhiMin(), zhp);
         glVertex3f(fM->fEtaMax, fM->GetPhiMax(), zhp);
         glVertex3f(fM->fEtaMin, fM->GetPhiMax(), zhp);
         glEnd();
      }
      glPopAttrib();
   }

   glPopMatrix();
}

////////////////////////////////////////////////////////////////////////////////
/// Processes tower selection from TGLViewer.

void TEveCaloLegoGL::ProcessSelection(TGLRnrCtx & /*rnrCtx*/, TGLSelectRecord & rec)
{
   TEveCaloData::vCellId_t sel;
   if (rec.GetN() > 2)
   {
      Int_t slice = rec.GetItem(1);
      Int_t cell  = rec.GetItem(2);

      if (fBinStep == 1)
      {
         Int_t tower = fM->fCellList[cell].fTower;
         while (cell > 0 && tower == fM->fCellList[cell].fTower)
         {
            sel.push_back(fM->fCellList[cell]);
            if (fCells3D) break;
            --cell;
         }
      }
      else
      {
         if (cell > 0)
         {
            Int_t nEta   = fEtaAxis->GetNbins();
            Int_t phiBin = Int_t(cell/(nEta+2));
            Int_t etaBin = cell - phiBin*(nEta+2);
            TEveCaloData::vCellId_t sl;
            fM->fData->GetCellList(fEtaAxis->GetBinCenter(etaBin), fEtaAxis->GetBinWidth(etaBin),
                                   fPhiAxis->GetBinCenter(phiBin), fPhiAxis->GetBinWidth(phiBin),
                                   sl);

            for (TEveCaloData::vCellId_i it = sl.begin(); it != sl.end(); ++it)
            {
               if (fCells3D) {
                  if ((*it).fSlice == slice ) sel.push_back(*it);
               } else {
                  if ((*it).fSlice <= slice ) sel.push_back(*it);
               }
            }
         }
      }
   }
   fM->fData->ProcessSelection(sel, rec);
}
