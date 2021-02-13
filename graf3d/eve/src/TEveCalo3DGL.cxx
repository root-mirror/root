// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TEveCalo3DGL.h"
#include "TEveCalo.h"

#include "TMath.h"
#include "TAxis.h"

#include "TGLRnrCtx.h"
#include "TGLSelectRecord.h"
#include "TGLPhysicalShape.h"
#include "TGLIncludes.h"
#include "TGLUtil.h"
#include "TEveRGBAPalette.h"
#include "TEveUtil.h"

/** \class TEveCalo3DGL
\ingroup TEve
OpenGL renderer class for TEveCalo3D.
*/

ClassImp(TEveCalo3DGL);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TEveCalo3DGL::TEveCalo3DGL() :
   TGLObject(), fM(0)
{
   fMultiColor = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set model object.

Bool_t TEveCalo3DGL::SetModel(TObject* obj, const Option_t* /*opt*/)
{
   fM = SetModelDynCast<TEveCalo3D>(obj);
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set bounding box.

void TEveCalo3DGL::SetBBox()
{
   // !! This ok if master sub-classed from TAttBBox
   SetAxisAlignedBBox(((TEveCalo3D*)fExternalObj)->AssertBBox());
}

////////////////////////////////////////////////////////////////////////////////
/// Override from TGLObject.
/// To account for large point-sizes we modify the projection matrix
/// during selection and thus we need a direct draw.

Bool_t TEveCalo3DGL::ShouldDLCache(const TGLRnrCtx& rnrCtx) const
{
   if (rnrCtx.Highlight() || rnrCtx.Selection()) return kFALSE;
   return TGLObject::ShouldDLCache(rnrCtx);
}

////////////////////////////////////////////////////////////////////////////////
/// Calculate cross-product.

inline void TEveCalo3DGL::CrossProduct(const Float_t a[3], const Float_t b[3],
                                       const Float_t c[3], Float_t out[3]) const
{
   const Float_t v1[3] = { a[0] - c[0], a[1] - c[1], a[2] - c[2] };
   const Float_t v2[3] = { b[0] - c[0], b[1] - c[1], b[2] - c[2] };

   out[0] = v1[1] * v2[2] - v1[2] * v2[1];
   out[1] = v1[2] * v2[0] - v1[0] * v2[2];
   out[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

////////////////////////////////////////////////////////////////////////////////
/// Render end cap grid.

void TEveCalo3DGL::RenderGridEndCap() const
{
   using namespace TMath;

   Float_t  rB = fM->GetBarrelRadius();
   Double_t zEF = fM->GetForwardEndCapPos();
   Double_t zEB = fM->GetBackwardEndCapPos();

   Float_t etaMin = fM->GetEtaMin();
   Float_t etaMax = fM->GetEtaMax();
   Float_t transF  = fM->GetTransitionEtaForward();
   Float_t transB  = fM->GetTransitionEtaBackward();
   Float_t phiMin = fM->GetPhiMin();
   Float_t phiMax = fM->GetPhiMax();

   TAxis *ax = fM->GetData()->GetEtaBins();
   Int_t  nx = ax->GetNbins();
   TAxis *ay = fM->GetData()->GetPhiBins();
   Int_t  ny = ay->GetNbins();


   Float_t r, z, theta, phiU, phiL, eta;

   // eta slices
   for (Int_t i=0; i<=nx; ++i)
   {
      eta = ax->GetBinUpEdge(i);
      if (eta >= transF && (eta > etaMin && eta < etaMax))
      {
         theta = TEveCaloData::EtaToTheta(eta);
         r = Abs(zEF*Tan(theta));
         z = Sign(zEF, ax->GetBinLowEdge(i));
         for (Int_t j=1; j<=ny; ++j)
         {
            phiL = ay->GetBinLowEdge(j);
            phiU = ay->GetBinUpEdge(j);
            if (TEveUtil::IsU1IntervalContainedByMinMax(phiMin, phiMax, phiL, phiU))
            {
               glVertex3f(r*Cos(phiL), r*Sin(phiL), z);
               glVertex3f(r*Cos(phiU), r*Sin(phiU), z);
            }
         }
      } else if (eta <= transB && (eta > etaMin && eta < etaMax)) {
         theta = TEveCaloData::EtaToTheta(eta);
         r = Abs(zEB*Tan(theta));
         z = Sign(zEB, ax->GetBinLowEdge(i));
         for (Int_t j=1; j<=ny; ++j)
         {
            phiL = ay->GetBinLowEdge(j);
            phiU = ay->GetBinUpEdge(j);
            if (TEveUtil::IsU1IntervalContainedByMinMax(phiMin, phiMax, phiL, phiU))
            {
               glVertex3f(r*Cos(phiL), r*Sin(phiL), z);
               glVertex3f(r*Cos(phiU), r*Sin(phiU), z);
            }
         }
      }
   }

   Float_t r1, r2;
   // phi slices front
   if (etaMax > transF)
   {
      r1 = zEF*Tan(TEveCaloData::EtaToTheta(etaMax));
      if (etaMin < transF)
         r2 = rB;
      else
         r2 = zEF*Tan(TEveCaloData::EtaToTheta(etaMin));

      for (Int_t j=1; j<=ny; ++j)
      {
         phiL = ay->GetBinLowEdge(j);
         phiU = ay->GetBinUpEdge(j);
         if (TEveUtil::IsU1IntervalContainedByMinMax(phiMin, phiMax, phiL, phiU))
         {
            glVertex3f( r1*Cos(phiU), r1*Sin(phiU), zEF);
            glVertex3f( r2*Cos(phiU), r2*Sin(phiU), zEF);
            glVertex3f( r1*Cos(phiL), r1*Sin(phiL), zEF);
            glVertex3f( r2*Cos(phiL), r2*Sin(phiL), zEF);
         }
      }
   }

   // phi slices back
   if (etaMin < transB)
   {
      r1 = zEB*Tan(TEveCaloData::EtaToTheta(etaMin));
      if (etaMax > transB)
         r2 = rB;
      else
         r2 = zEB*Tan(TEveCaloData::EtaToTheta(etaMax));

      r1 = Abs(r1);
      r2 = Abs(r2);
      for (Int_t j=1; j<=ny; ++j)
      {
         phiL = ay->GetBinLowEdge(j);
         phiU = ay->GetBinUpEdge(j);
         if (TEveUtil::IsU1IntervalContainedByMinMax(phiMin, phiMax, phiL, phiU))
         {
            glVertex3f(r1*Cos(phiU), r1*Sin(phiU), zEB);
            glVertex3f(r2*Cos(phiU), r2*Sin(phiU), zEB);
            glVertex3f(r1*Cos(phiL), r1*Sin(phiL), zEB);
            glVertex3f(r2*Cos(phiL), r2*Sin(phiL), zEB);
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Render barrel grid.

void TEveCalo3DGL::RenderGridBarrel() const
{
   using namespace TMath;

   Float_t etaMin = fM->GetEtaMin();
   Float_t etaMax = fM->GetEtaMax();
   Float_t transF  = fM->GetTransitionEtaForward();
   Float_t transB  = fM->GetTransitionEtaBackward();
   Float_t phiMin = fM->GetPhiMin();
   Float_t phiMax = fM->GetPhiMax();

   Float_t rB = fM->GetBarrelRadius();
   TAxis *ax  = fM->GetData()->GetEtaBins();
   Int_t nx   = ax->GetNbins();
   TAxis *ay  = fM->GetData()->GetPhiBins();
   Int_t ny   = ay->GetNbins();

   Float_t z, theta, phiL, phiU, eta, x, y;

   // eta slices
   for(Int_t i=0; i<=nx; i++)
   {
      eta = ax->GetBinUpEdge(i);
      if (eta<=transF && eta>=transB && (etaMin < eta && eta < etaMax))
      {
         theta = TEveCaloData::EtaToTheta(eta);
         z  = rB/Tan(theta);
         for (Int_t j=1; j<=ny; j++)
         {
            phiU = ay->GetBinUpEdge(j);
            phiL = ay->GetBinLowEdge(j);
            if (TEveUtil::IsU1IntervalContainedByMinMax(phiMin, phiMax, phiL, phiU))
            {
               glVertex3f(rB*Cos(phiL), rB*Sin(phiL), z);
               glVertex3f(rB*Cos(phiU), rB*Sin(phiU), z);
            }
         }
      }
   }

   // phi slices
   Float_t zF, zB;

   if (etaMin > transB)
      zB = rB/Tan(TEveCaloData::EtaToTheta(etaMin));
   else
      zB = fM->GetBackwardEndCapPos();


   if (etaMax < transF)
      zF =  rB/Tan(TEveCaloData::EtaToTheta(etaMax));
   else
      zF = fM->GetForwardEndCapPos();

   for (Int_t j=1; j<=ny; j++)
   {
      phiU = ay->GetBinUpEdge(j);
      phiL = ay->GetBinLowEdge(j);
      if (TEveUtil::IsU1IntervalContainedByMinMax(phiMin, phiMax, phiL, phiU))
      {
         x = rB * Cos(phiL);
         y = rB * Sin(phiL);
         glVertex3f(x, y, zB);
         glVertex3f(x, y, zF);
         x = rB * Cos(phiU);
         y = rB * Sin(phiU);
         glVertex3f(x, y, zB);
         glVertex3f(x, y, zF);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Draw frame reading eta, phi axis.

void TEveCalo3DGL::RenderGrid(TGLRnrCtx & rnrCtx) const
{
   if (rnrCtx.Highlight() || rnrCtx.Selection() || rnrCtx.IsDrawPassOutlineLine()) return;

   Bool_t transparent_p = fM->fFrameTransparency > 0;

   if (transparent_p)
   {
      glPushAttrib(GL_ENABLE_BIT | GL_DEPTH_BUFFER_BIT);

      glDepthMask(GL_FALSE);
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

      TGLUtil::ColorTransparency(fM->fFrameColor, fM->fFrameTransparency);
   }

   TGLCapabilitySwitch lights_off(GL_LIGHTING, kFALSE);

   TGLUtil::LineWidth(fM->GetFrameWidth());
   glBegin(GL_LINES);

   Float_t etaMin = fM->GetEtaMin();
   Float_t etaMax = fM->GetEtaMax();

   Float_t transF  = fM->GetTransitionEtaForward();
   Float_t transB  = fM->GetTransitionEtaBackward();
   if (fM->GetRnrBarrelFrame() && (etaMin < transF && etaMax > transB))
   {
      RenderGridBarrel();
   }

   if (fM->GetRnrEndCapFrame() && (etaMax > transF || etaMin < transB))
   {
      RenderGridEndCap();
   }

   glEnd();

   if (transparent_p)
   {
      glPopAttrib();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Render box with given points.
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

void TEveCalo3DGL::RenderBox(const Float_t pnts[8]) const
{
   const Float_t *p = pnts;
   Float_t cross[3];

   // bottom: 0123
   glBegin(GL_POLYGON);
   CrossProduct(p+3, p+9, p, cross);
   glNormal3fv(cross);
   glVertex3fv(p);
   glVertex3fv(p+3);
   glVertex3fv(p+6);
   glVertex3fv(p+9);
   glEnd();
   // top:    7654
   glBegin(GL_POLYGON);
   CrossProduct(p+21, p+15, p+12, cross);
   glNormal3fv(cross);
   glVertex3fv(p+21);
   glVertex3fv(p+18);
   glVertex3fv(p+15);
   glVertex3fv(p+12);
   glEnd();
   // back:   0451
   glBegin(GL_POLYGON);
   CrossProduct(p+12, p+3, p, cross);
   glNormal3fv(cross);
   glVertex3fv(p);
   glVertex3fv(p+12);
   glVertex3fv(p+15);
   glVertex3fv(p+3);
   glEnd();
   //front :  3267
   glBegin(GL_POLYGON);
   CrossProduct(p+6, p+21, p+9, cross);
   glNormal3fv(cross);
   glVertex3fv(p+9);
   glVertex3fv(p+6);
   glVertex3fv(p+18);
   glVertex3fv(p+21);
   glEnd();
   // left:    0374
   glBegin(GL_POLYGON);
   CrossProduct(p+21, p, p+9, cross);
   glNormal3fv(cross);
   glVertex3fv(p);
   glVertex3fv(p+9);
   glVertex3fv(p+21);
   glVertex3fv(p+12);
   glEnd();
   // right:   1562
   glBegin(GL_POLYGON);
   CrossProduct(p+15, p+6, p+3, cross);
   glNormal3fv(cross);
   glVertex3fv(p+3);
   glVertex3fv(p+15);
   glVertex3fv(p+18);
   glVertex3fv(p+6);
   glEnd();
}

////////////////////////////////////////////////////////////////////////////////
/// Render barrel cell.

void TEveCalo3DGL::RenderBarrelCell(const TEveCaloData::CellGeom_t &cellData, Float_t towerH, Float_t& offset ) const
{
   using namespace TMath;

   Float_t r1 = fM->GetBarrelRadius() + offset;
   Float_t r2 = r1 + towerH*Sin(cellData.ThetaMin());
   Float_t z1In, z1Out, z2In, z2Out;

   z1In  = r1/Tan(cellData.ThetaMax());
   z1Out = r2/Tan(cellData.ThetaMax());
   z2In  = r1/Tan(cellData.ThetaMin());
   z2Out = r2/Tan(cellData.ThetaMin());

   Float_t cos1 = Cos(cellData.PhiMin());
   Float_t sin1 = Sin(cellData.PhiMin());
   Float_t cos2 = Cos(cellData.PhiMax());
   Float_t sin2 = Sin(cellData.PhiMax());

   Float_t box[24];
   Float_t* pnts = box;
   // 0
   pnts[0] = r1*cos2;
   pnts[1] = r1*sin2;
   pnts[2] = z1In;
   pnts += 3;
   // 1
   pnts[0] = r1*cos1;
   pnts[1] = r1*sin1;
   pnts[2] = z1In;
   pnts += 3;
   // 2
   pnts[0] = r1*cos1;
   pnts[1] = r1*sin1;
   pnts[2] = z2In;
   pnts += 3;
   // 3
   pnts[0] = r1*cos2;
   pnts[1] = r1*sin2;
   pnts[2] = z2In;
   pnts += 3;
   //---------------------------------------------------
   // 4
   pnts[0] = r2*cos2;
   pnts[1] = r2*sin2;
   pnts[2] = z1Out;
   pnts += 3;
   // 5
   pnts[0] = r2*cos1;
   pnts[1] = r2*sin1;
   pnts[2] = z1Out;
   pnts += 3;
   // 6
   pnts[0] = r2*cos1;
   pnts[1] = r2*sin1;
   pnts[2] = z2Out;
   pnts += 3;
   // 7
   pnts[0] = r2*cos2;
   pnts[1] = r2*sin2;
   pnts[2] = z2Out;

   RenderBox(box);

   offset += towerH*Sin(cellData.ThetaMin());

}// end RenderBarrelCell

////////////////////////////////////////////////////////////////////////////////
/// Render an endcap cell.

void TEveCalo3DGL::RenderEndCapCell(const TEveCaloData::CellGeom_t &cellData, Float_t towerH, Float_t& offset ) const
{
   using namespace TMath;
   Float_t z1, r1In, r1Out, z2, r2In, r2Out;

   z1    = (cellData.EtaMin()<0) ? fM->fEndCapPosB - offset : fM->fEndCapPosF + offset;
   z2    = z1 + TMath::Sign(towerH, cellData.EtaMin());

   r1In  = z1*Tan(cellData.ThetaMin());
   r2In  = z2*Tan(cellData.ThetaMin());
   r1Out = z1*Tan(cellData.ThetaMax());
   r2Out = z2*Tan(cellData.ThetaMax());

   Float_t cos2 = Cos(cellData.PhiMin());
   Float_t sin2 = Sin(cellData.PhiMin());
   Float_t cos1 = Cos(cellData.PhiMax());
   Float_t sin1 = Sin(cellData.PhiMax());

   Float_t box[24];
   Float_t* pnts = box;
   // 0
   pnts[0] = r1In*cos1;
   pnts[1] = r1In*sin1;
   pnts[2] = z1;
   pnts += 3;
   // 1
   pnts[0] = r1In*cos2;
   pnts[1] = r1In*sin2;
   pnts[2] = z1;
   pnts += 3;
   // 2
   pnts[0] = r2In*cos2;
   pnts[1] = r2In*sin2;
   pnts[2] = z2;
   pnts += 3;
   // 3
   pnts[0] = r2In*cos1;
   pnts[1] = r2In*sin1;
   pnts[2] = z2;
   pnts += 3;
   //---------------------------------------------------
   // 4
   pnts[0] = r1Out*cos1;
   pnts[1] = r1Out*sin1;
   pnts[2] = z1;
   pnts += 3;
   // 5
   pnts[0] = r1Out*cos2;
   pnts[1] = r1Out*sin2;
   pnts[2] = z1;
   pnts += 3;
   // 6
   pnts[0] = r2Out*cos2;
   pnts[1] = r2Out*sin2;
   pnts[2] = z2;
   pnts += 3;
   // 7
   pnts[0] = r2Out*cos1;
   pnts[1] = r2Out*sin1;
   pnts[2] = z2;

   RenderBox(box);

   offset += towerH;

} // end RenderEndCapCell

////////////////////////////////////////////////////////////////////////////////
/// GL rendering.

void TEveCalo3DGL::DirectDraw(TGLRnrCtx &rnrCtx) const
{
   if ( fM->GetValueIsColor())  fM->AssertPalette();

   // check if eta phi range has changed
   if (fM->fCellIdCacheOK == kFALSE)
      fM->BuildCellIdCache();

   glPushAttrib(GL_ENABLE_BIT | GL_LINE_BIT | GL_POLYGON_BIT);
   glEnable(GL_LIGHTING);
   glEnable(GL_NORMALIZE);
   glEnable(GL_BLEND);
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

   TEveCaloData::CellData_t cellData;
   Float_t towerH = 0;
   Int_t   tower = 0;
   Int_t   prevTower = -1;
   Float_t offset = 0;
   Int_t cellID = 0;

   if (rnrCtx.SecSelection()) glPushName(0);

   fOffset.assign(fM->fCellList.size(), 0);
   for (TEveCaloData::vCellId_i i = fM->fCellList.begin(); i != fM->fCellList.end(); ++i)
   {
      fM->fData->GetCellData((*i), cellData);
      tower = i->fTower;
      if (tower != prevTower)
      {
         offset = 0;
         prevTower = tower;
      }
      fOffset[cellID] = offset;
      fM->SetupColorHeight(cellData.Value(fM->fPlotEt), (*i).fSlice, towerH);

      if (rnrCtx.SecSelection()) glLoadName(cellID);

      if ((cellData.Eta() > 0 && cellData.Eta() < fM->GetTransitionEtaForward()) ||
          (cellData.Eta() < 0 && cellData.Eta() > fM->GetTransitionEtaBackward()))
      {
         RenderBarrelCell(cellData, towerH, offset);
      }
      else
      {
         RenderEndCapCell(cellData, towerH, offset);
      }
      ++cellID;
   }

   if (rnrCtx.SecSelection()) glPopName();

   RenderGrid(rnrCtx);

   glPopAttrib();
}

////////////////////////////////////////////////////////////////////////////////
/// Draw polygons in highlight mode.

void TEveCalo3DGL::DrawHighlight(TGLRnrCtx & rnrCtx, const TGLPhysicalShape* /*pshp*/, Int_t /*lvl*/) const
{
   if (fM->fData->GetCellsSelected().empty() && fM->fData->GetCellsHighlighted().empty())
   {
      return;
   }

   glPushAttrib(GL_ENABLE_BIT | GL_LINE_BIT | GL_POLYGON_BIT);
   glDisable(GL_LIGHTING);
   glDisable(GL_CULL_FACE);
   glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

   TGLUtil::LineWidth(2);
   TGLUtil::LockColor();

   if (!fM->fData->GetCellsHighlighted().empty())
   {
      glColor4ubv(rnrCtx.ColorSet().Selection(3).CArr());
      DrawSelectedCells(fM->fData->GetCellsHighlighted());
   }
   if (!fM->fData->GetCellsSelected().empty())
   {
      Float_t dr[2];
      glGetFloatv(GL_DEPTH_RANGE,dr);
      glColor4ubv(rnrCtx.ColorSet().Selection(1).CArr());
      glDepthRange(dr[0], 0.8*dr[1]);
      DrawSelectedCells(fM->fData->GetCellsSelected());
      glDepthRange(dr[0], dr[1]);
   }

   TGLUtil::UnlockColor();
   glPopAttrib();
}

////////////////////////////////////////////////////////////////////////////////

void TEveCalo3DGL::DrawSelectedCells(TEveCaloData::vCellId_t cells) const
{
   TEveCaloData::CellData_t cellData;
   Float_t towerH = 0;

   for (TEveCaloData::vCellId_i i = cells.begin(); i != cells.end(); i++)
   {
      fM->fData->GetCellData(*i, cellData);
      fM->SetupColorHeight(cellData.Value(fM->fPlotEt), (*i).fSlice, towerH);

      // find tower with offsets
      Float_t offset = 0;
      for (Int_t j = 0; j < (Int_t) fM->fCellList.size(); ++j)
      {
         if (fM->fCellList[j].fTower == i->fTower && fM->fCellList[j].fSlice == i->fSlice )
         {
            offset = fOffset[j];
            break;
         }
      }

      if (fM->CellInEtaPhiRng(cellData))
      {
         if (cellData.Eta() < fM->GetTransitionEtaForward() && cellData.Eta() > fM->GetTransitionEtaBackward())
            RenderBarrelCell(cellData, towerH, offset);
         else
            RenderEndCapCell(cellData, towerH, offset);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Processes tower selection.
/// Virtual function from TGLogicalShape. Called from TGLViewer.

void TEveCalo3DGL::ProcessSelection(TGLRnrCtx& /*rnrCtx*/, TGLSelectRecord& rec)
{
   TEveCaloData::vCellId_t sel;
   if (rec.GetN() > 1)
   {
      sel.push_back(fM->fCellList[rec.GetItem(1)]);
   }
   fM->fData->ProcessSelection(sel, rec);
}
