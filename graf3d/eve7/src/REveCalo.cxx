// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "ROOT/REveCalo.hxx"
#include "ROOT/REveCaloData.hxx"
#include "ROOT/REveProjections.hxx"
#include "ROOT/REveProjectionManager.hxx"
#include "ROOT/REveRGBAPalette.hxx"
//#include "ROOT/REveText.hxx"
#include "ROOT/REveRenderData.hxx"
#include "ROOT/REveTrans.hxx"

#include "TClass.h"
#include "TMathBase.h"
#include "TMath.h"
#include "TAxis.h"

// #include "TGLUtil.h"

#include <cassert>
#include <iostream>

using namespace ROOT::Experimental;

/** \class REveCaloViz
\ingroup REve
Base class for calorimeter data visualization.
See REveCalo2D and REveCalo3D for concrete implementations.
*/


////////////////////////////////////////////////////////////////////////////////

REveCaloViz::REveCaloViz(REveCaloData* data, const char* n, const char* t) :
   REveElement(),
   REveProjectable(),

   fData(0),
   fCellIdCacheOK(kFALSE),

   fEtaMin(-10),
   fEtaMax(10),

   fPhi(0.),
   fPhiOffset(TMath::Pi()),

   fAutoRange(kTRUE),

   fBarrelRadius(-1.f),
   fEndCapPosF(-1.f),
   fEndCapPosB(-1.f),

   fPlotEt(kTRUE),

   fMaxTowerH(100),
   fScaleAbs(kFALSE),
   fMaxValAbs(100),

   fValueIsColor(kFALSE),
   fPalette(0)
{
   // Constructor.

   fPickable = kTRUE;
   SetNameTitle(n, t);
   SetData(data);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

REveCaloViz::~REveCaloViz()
{
   if (fPalette) fPalette->DecRefCount();
}

////////////////////////////////////////////////////////////////////////////////
/// Get threshold for given slice.

Float_t REveCaloViz::GetDataSliceThreshold(Int_t slice) const
{
   return fData->RefSliceInfo(slice).fThreshold;
}

////////////////////////////////////////////////////////////////////////////////
/// Management of selection state and ownership of selected cell list
/// is done in REveCaloData. This is a reason selection is forwarded to it.

REveElement* REveCaloViz::ForwardSelection()
{
   return fData;
}

////////////////////////////////////////////////////////////////////////////////
/// Management of selection state and ownership of selected cell list
/// is done in REveCaloData. We still want GUI editor to display
/// concrete calo-viz object.

REveElement* REveCaloViz::ForwardEdit()
{
   return this;
}

////////////////////////////////////////////////////////////////////////////////
/// Set threshold for given slice.

void REveCaloViz::SetDataSliceThreshold(Int_t slice, Float_t val)
{
   fData->SetSliceThreshold(slice, val);
}

////////////////////////////////////////////////////////////////////////////////
/// Get slice color from data.

Color_t REveCaloViz::GetDataSliceColor(Int_t slice) const
{
   return fData->RefSliceInfo(slice).fColor;
}

////////////////////////////////////////////////////////////////////////////////
/// Set slice color in data.

void REveCaloViz::SetDataSliceColor(Int_t slice, Color_t col)
{
   fData->SetSliceColor(slice, col);
}

////////////////////////////////////////////////////////////////////////////////
/// Set eta range.

void REveCaloViz::SetEta(Float_t l, Float_t u)
{
   fEtaMin=l;
   fEtaMax=u;

   InvalidateCellIdCache();
}

////////////////////////////////////////////////////////////////////////////////
/// Set E/Et plot.

void REveCaloViz::SetPlotEt(Bool_t isEt)
{
   fPlotEt=isEt;
   if (fPalette)
      fPalette->SetLimits(0, TMath::CeilNint(GetMaxVal()));

   InvalidateCellIdCache();
}

////////////////////////////////////////////////////////////////////////////////

Float_t REveCaloViz::GetMaxVal() const
{
   // Get maximum plotted value.

   return fData->GetMaxVal(fPlotEt);

}

////////////////////////////////////////////////////////////////////////////////
/// Set phi range.

void REveCaloViz::SetPhiWithRng(Float_t phi, Float_t rng)
{
   using namespace TMath;

   fPhi = phi;
   fPhiOffset = rng;

   InvalidateCellIdCache();
}

////////////////////////////////////////////////////////////////////////////////
/// Get transition angle between barrel and end-cap cells, assuming fEndCapPosF = -fEndCapPosB.

Float_t REveCaloViz::GetTransitionTheta() const
{
   return TMath::ATan(fBarrelRadius/fEndCapPosF);
}

////////////////////////////////////////////////////////////////////////////////
/// Get transition eta between barrel and end-cap cells, assuming fEndCapPosF = -fEndCapPosB.

Float_t REveCaloViz::GetTransitionEta() const
{
   using namespace TMath;
   Float_t t = GetTransitionTheta()*0.5f;
   return -Log(Tan(t));
}

////////////////////////////////////////////////////////////////////////////////
/// Get transition angle between barrel and forward end-cap cells.

Float_t REveCaloViz::GetTransitionThetaForward() const
{
   return TMath::ATan(fBarrelRadius/fEndCapPosF);
}

////////////////////////////////////////////////////////////////////////////////
/// Get transition eta between barrel and forward end-cap cells.

Float_t REveCaloViz::GetTransitionEtaForward() const
{
   using namespace TMath;
   Float_t t = GetTransitionThetaForward()*0.5f;
   return -Log(Tan(t));
}

////////////////////////////////////////////////////////////////////////////////
/// Get transition angle between barrel and backward end-cap cells.

Float_t REveCaloViz::GetTransitionThetaBackward() const
{
   return TMath::ATan(fBarrelRadius/fEndCapPosB);
}

////////////////////////////////////////////////////////////////////////////////
/// Get transition eta between barrel and backward end-cap cells.

Float_t REveCaloViz::GetTransitionEtaBackward() const
{
   using namespace TMath;
   Float_t t = GetTransitionThetaBackward()*0.5f;
   //negative theta means negative eta
   return Log(-Tan(t));
}


////////////////////////////////////////////////////////////////////////////////
/// Set calorimeter event data.

void REveCaloViz::SetData(REveCaloData* data)
{

   if (data == fData) return;
   fData = data;
   //   SetSelectionMaster(data);
   if (fData)
   {
      fData->AddNiece(this);
      DataChanged();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Update setting and cache on data changed.
/// Called from REvecaloData::BroadcastDataChange()

void REveCaloViz::DataChanged()
{
   Double_t min, max, delta;

   fData->GetEtaLimits(min, max);
   if (fAutoRange) {
      fEtaMin = min;
      fEtaMax = max;
   } else {
      if (fEtaMin < min) fEtaMin = min;
      if (fEtaMax > max) fEtaMax = max;
   }

   fData->GetPhiLimits(min, max);
   delta = 0.5*(max - min);
   if (fAutoRange || fPhi < min || fPhi > max) {
      fPhi       = 0.5*(max + min);
      fPhiOffset = delta;
   } else {
      if (fPhiOffset > delta) fPhiOffset = delta;
   }

   if (fPalette)
   {
      Int_t hlimit = TMath::CeilNint(GetMaxVal());
      fPalette->SetLimits(0, hlimit);
      fPalette->SetMin(0);
      fPalette->SetMax(hlimit);
   }

   InvalidateCellIdCache();
}

////////////////////////////////////////////////////////////////////////////////
/// Assert cell id cache is ok.
/// Returns true if the cache has been updated.

Bool_t REveCaloViz::AssertCellIdCache() const
{
   REveCaloViz* cv = const_cast<REveCaloViz*>(this);
   if (!fCellIdCacheOK) {
      cv->BuildCellIdCache();
      return kTRUE;
   } else {
      return kFALSE;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns true if given cell is in the ceta phi range.

Bool_t REveCaloViz::CellInEtaPhiRng(REveCaloData::CellData_t& cellData) const
{
   if (cellData.EtaMin() >= fEtaMin && cellData.EtaMax() <= fEtaMax)
   {
      if (REveUtil::IsU1IntervalContainedByMinMax
          (fPhi-fPhiOffset, fPhi+fPhiOffset, cellData.PhiMin(), cellData.PhiMax()))
         return kTRUE;
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Assign parameters from given model.

void REveCaloViz::AssignCaloVizParameters(REveCaloViz* m)
{
   SetData(m->fData);

   fEtaMin    = m->fEtaMin;
   fEtaMax    = m->fEtaMax;

   fPhi       = m->fPhi;
   fPhiOffset = m->fPhiOffset;

   fBarrelRadius = m->fBarrelRadius;
   fEndCapPosF    = m->fEndCapPosF;
   fEndCapPosB    = m->fEndCapPosB;

   if (m->fPalette)
   {
      REveRGBAPalette& mp = * m->fPalette;
      if (fPalette) fPalette->DecRefCount();
      fPalette = new REveRGBAPalette(mp.GetMinVal(), mp.GetMaxVal(), mp.GetInterpolate());
      fPalette->SetDefaultColor(mp.GetDefaultColor());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set REveRGBAPalette object pointer.

void REveCaloViz::SetPalette(REveRGBAPalette* p)
{
   if ( fPalette == p) return;
   if (fPalette) fPalette->DecRefCount();
   fPalette = p;
   if (fPalette) fPalette->IncRefCount();
}

////////////////////////////////////////////////////////////////////////////////
/// Get transformation factor from E/Et to height

Float_t REveCaloViz::GetValToHeight() const
{
   if (fScaleAbs)
   {
      return fMaxTowerH/fMaxValAbs;
   }
   else
   {
      if (fData->Empty()) {
         assert(false);
       return 1;
      }
      return fMaxTowerH/fData->GetMaxVal(fPlotEt);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Make sure the REveRGBAPalette pointer is not null.
/// If it is not set, a new one is instantiated and the range is set
/// to current min/max signal values.

REveRGBAPalette* REveCaloViz::AssertPalette()
{
   if (fPalette == 0) {
      fPalette = new REveRGBAPalette;
      fPalette->SetDefaultColor((Color_t)4);

      Int_t hlimit = TMath::CeilNint(GetMaxVal());
      fPalette->SetLimits(0, hlimit);
      fPalette->SetMin(0);
      fPalette->SetMax(hlimit);

   }
   return fPalette;
}


////////////////////////////////////////////////////////////////////////////////
/// Virtual from REveProjectable, returns REveCalo2D class.

TClass* REveCaloViz::ProjectedClass(const REveProjection*) const
{
   return TClass::GetClass<REveCalo2D>();
}

////////////////////////////////////////////////////////////////////////////////
/// Set color and height for a given value and slice using slice color or REveRGBAPalette.

void REveCaloViz::SetupHeight(Float_t value, Int_t /*slice*/, Float_t& outH) const
{
   if (fValueIsColor)
   {
      outH = GetValToHeight()*fData->GetMaxVal(fPlotEt);
      assert("fValueIsColor" && false);
      // UChar_t c[4];
      // fPalette->ColorFromValue((Int_t)value, c);
      // c[3] = fData->GetSliceTransparency(slice);
      // TGLUtil::Color4ubv(c);
   }
   else
   {
      // TGLUtil::ColorTransparency(fData->GetSliceColor(slice), fData->GetSliceTransparency(slice));
      outH = GetValToHeight()*value;
   }
}
///////////////////////////////////////////////////////////////////////////////
/// Fill core part of JSON representation.

Int_t REveCaloViz::WriteCoreJson(nlohmann::json &j, Int_t rnr_offset)
{
   // The slice colors need to be streamed becuse at EveElement contruction time, streamed caloData
   // is not available. Maybe this is not necessary if EveElements have EveManager globaly available

   Int_t ret = REveElement::WriteCoreJson(j, rnr_offset);
   j["dataId"] = fData->GetElementId();
   j["sliceColors"] =  nlohmann::json::array();
   for (int i = 0; i < fData->GetNSlices(); ++i)
   {
      j["sliceColors"].push_back(fData->GetSliceColor(i));
      }
   j["fSecondarySelect"] = true;
   return ret;
}

/** \class REveCalo3D
\ingroup REve
Visualization of a calorimeter event data in 3D.
*/


////////////////////////////////////////////////////////////////////////////////
/// Constructor.

REveCalo3D::REveCalo3D(REveCaloData* d, const char* n, const char* t):
   REveCaloViz(d, n, t),

   fRnrEndCapFrame    (kTRUE),
   fRnrBarrelFrame    (kTRUE),
   fFrameWidth        (0.5),
   fFrameColor        (kGray+1),
   fFrameTransparency (80)
{
   fCanEditMainColor        = kTRUE;
   fCanEditMainTransparency = kTRUE;
   fMainColorPtr = &fFrameColor;
}


////////////////////////////////////////////////////////////////////////////////
/// Make endcap cell
//
void REveCalo3D::MakeEndCapCell(const REveCaloData::CellGeom_t &cellData, float towerH, Float_t& offset, float *pnts) const
{
   using namespace TMath;
   Float_t z1, r1In, r1Out, z2, r2In, r2Out;

   z1    = (cellData.EtaMin()<0) ? fEndCapPosB - offset : fEndCapPosF + offset;
   z2    = z1 + TMath::Sign(towerH, cellData.EtaMin());

   r1In  = z1*Tan(cellData.ThetaMin());
   r2In  = z2*Tan(cellData.ThetaMin());
   r1Out = z1*Tan(cellData.ThetaMax());
   r2Out = z2*Tan(cellData.ThetaMax());

   Float_t cos2 = Cos(cellData.PhiMin());
   Float_t sin2 = Sin(cellData.PhiMin());
   Float_t cos1 = Cos(cellData.PhiMax());
   Float_t sin1 = Sin(cellData.PhiMax());

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
   // 6s
   pnts[0] = r2Out*cos2;
   pnts[1] = r2Out*sin2;
   pnts[2] = z2;
   pnts += 3;
   // 7
   pnts[0] = r2Out*cos1;
   pnts[1] = r2Out*sin1;
   pnts[2] = z2;

   offset += towerH;
}

////////////////////////////////////////////////////////////////////////////////
/// Make endcap cell
//
void REveCalo3D::MakeBarrelCell(const REveCaloData::CellGeom_t &cellData, float towerH, Float_t& offset, float *pnts) const
{
   using namespace TMath;

   float r1 = GetBarrelRadius() + offset;
   float r2 = r1 + towerH*Sin(cellData.ThetaMin());
   float z1In, z1Out, z2In, z2Out;

   z1In  = r1/Tan(cellData.ThetaMax());
   z1Out = r2/Tan(cellData.ThetaMax());
   z2In  = r1/Tan(cellData.ThetaMin());
   z2Out = r2/Tan(cellData.ThetaMin());

   float cos1 = Cos(cellData.PhiMin());
   float sin1 = Sin(cellData.PhiMin());
   float cos2 = Cos(cellData.PhiMax());
   float sin2 = Sin(cellData.PhiMax());

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


   offset += towerH*Sin(cellData.ThetaMin());

}

////////////////////////////////////////////////////////////////////////////////
/// Crates 3D point array for rendering.

void REveCalo3D::BuildRenderData()
{
   AssertCellIdCache();
   if (fCellList.empty())
   return;

   REveCaloData::CellData_t cellData;
   Float_t towerH = 0;
   Int_t   tower = 0;
   Int_t   prevTower = -1;
   Float_t offset = 0;

   fRenderData = std::make_unique<REveRenderData>("makeCalo3D");
   float pnts[24];
   for (REveCaloData::vCellId_i i = fCellList.begin(); i != fCellList.end(); ++i)
   {
      fData->GetCellData((*i), cellData);
      tower = i->fTower;
      if (tower != prevTower)
      {
         offset = 0;
         prevTower = tower;
      }
      // fOffset[cellID] = offset; this is needed to be stored for selection

      SetupHeight(cellData.Value(fPlotEt), (*i).fSlice, towerH);

      if ((cellData.Eta() > 0 && cellData.Eta() < GetTransitionEtaForward()) ||
          (cellData.Eta() < 0 && cellData.Eta() > GetTransitionEtaBackward()))
      {
         MakeBarrelCell(cellData, towerH, offset, pnts);
      }
      else
      {
         MakeEndCapCell(cellData, towerH, offset, pnts);
      }
      /*
      printf(" REveCalo3D::BuildRenderData push box vertces -------------------------\n");
      for (int t = 0; t < 8; ++t)
      {
         printf("(%f %f %f)\n", pnts[t*3],  pnts[t*3+1], pnts[t*3+2] );
      }
      */
      fRenderData->PushV(pnts, 24);

      //      REveCaloData::SliceInfo_t& sliceInfo = fData->RefSliceInfo(i->fSlice);
      fRenderData->PushI( i->fSlice);
      fRenderData->PushI( i->fTower);
      fRenderData->PushN(cellData.Value(fPlotEt));
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Fill core part of JSON representation.

Int_t REveCalo3D::WriteCoreJson(nlohmann::json &j, Int_t rnr_offset)
{
   return REveCaloViz::WriteCoreJson(j, rnr_offset);
}

////////////////////////////////////////////////////////////////////////////////
/// Fill core part of JSON representation for selection.

void REveCalo3D::WriteCoreJsonSelection(nlohmann::json &j, REveCaloData::vCellId_t cells)
{
   // selection
   auto sarr = nlohmann::json::array();
   REveCaloData::CellData_t cellData;
   for (REveCaloData::vCellId_i i = cells.begin(); i != cells.end(); i++)
   {
      fData->GetCellData(*i, cellData);
      if (CellInEtaPhiRng(cellData))
      {
         nlohmann::json jsc;
         jsc["t"] = i->fTower;
         jsc["s"] = i->fSlice;
         jsc["f"] = i->fFraction;
         sarr.push_back(jsc);
      }
   }

   nlohmann::json rec = {};
   rec["caloVizId"] = GetElementId();
   rec["cells"] = sarr;

   j.push_back(rec);
}


////////////////////////////////////////////////////////////////////////////////
/// Build list of drawn cell IDs. See REveCalo3DGL::DirectDraw().

void REveCalo3D::BuildCellIdCache()
{
   fCellList.clear();

   fData->GetCellList(GetEta(), GetEtaRng(), GetPhi(), GetPhiRng(), fCellList);
   fCellIdCacheOK = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill bounding-box information of the base-class TAttBBox (virtual method).
/// If member 'REveFrameBox* fFrame' is set, frame's corners are used as bbox.

void REveCalo3D::ComputeBBox()
{
   BBoxInit();

   Float_t th = (fData) ? GetValToHeight() * fData->GetMaxVal(fPlotEt) : 0;

   fBBox[0] = -fBarrelRadius - th;
   fBBox[1] =  fBarrelRadius + th;
   fBBox[2] =  fBBox[0];
   fBBox[3] =  fBBox[1];
   fBBox[4] =  fEndCapPosB - th;
   fBBox[5] =  fEndCapPosF + th;
}

/** \class REveCalo2D
\ingroup REve
Visualization of a calorimeter event data in 2D.
*/

////////////////////////////////////////////////////////////////////////////////
/// Client selection callback

void REveCalo3D::NewTowerPicked(Int_t tower, Int_t slice, Int_t selectionId, bool multi)
{
   REveCaloData::CellId_t cell(tower, slice, 1.0f);
   REveCaloData::vCellId_t sel;

   sel.push_back(cell);
   fData->ProcessSelection(sel, selectionId, multi);
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor.

REveCalo2D::REveCalo2D(const char* n, const char* t):
   REveCaloViz(0, n, t),
   REveProjected(),
   fOldProjectionType(REveProjection::kPT_Unknown),
   fMaxESumBin( 0),
   fMaxEtSumBin(0)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

REveCalo2D::~REveCalo2D()
{
   REveCaloData::vCellId_t* cids;
   UInt_t n;

   // clear selected cell ids
   n = fCellListsSelected.size();
   for(UInt_t i = 0; i < n; ++i) {
      cids = fCellListsSelected[i];
      if (cids) {
         cids->clear(); delete cids;
      }
   }
   fCellListsSelected.clear();

   // clear all cell dds
   n = fCellLists.size();
   for(UInt_t i = 0; i < n; ++i) {
      cids = fCellLists[i];
      if (cids) {
         cids->clear(); delete cids;
      }
   }
   fCellLists.clear();
}

////////////////////////////////////////////////////////////////////////////////
/// This is virtual method from base-class REveProjected.

void REveCalo2D::UpdateProjection()
{
   if (fManager->GetProjection()->GetType() != fOldProjectionType)
   {
      fCellIdCacheOK=kFALSE;
      fOldProjectionType = fManager->GetProjection()->GetType();
   }
   ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Set projection manager and model object.

void REveCalo2D::SetProjection(REveProjectionManager* mng, REveProjectable* model)
{
   REveProjected::SetProjection(mng, model);
   REveCaloViz* viz = dynamic_cast<REveCaloViz*>(model);
   AssignCaloVizParameters(viz);
}

////////////////////////////////////////////////////////////////////////////////
/// Is current projection type RPhi

Bool_t REveCalo2D::IsRPhi() const
{
   return fManager->GetProjection()->GetType() == REveProjection::kPT_RPhi;
}

////////////////////////////////////////////////////////////////////////////////
/// Build lists of drawn cell IDs. See REveCalo2DGL::DirecDraw().

void REveCalo2D::BuildCellIdCache()
{
   // clear old cache
   for (vBinCells_i it = fCellLists.begin(); it != fCellLists.end(); it++)
   {
      if (*it)
      {
         (*it)->clear();
         delete *it;
      }
   }
   fCellLists.clear();
   fCellLists.push_back(0);

   REveProjection::EPType_e pt = fManager->GetProjection()->GetType();
   REveCaloData::vCellId_t* clv; // ids per phi bin in r-phi projection else ids per eta bins in rho-z projection

   Bool_t isRPhi = (pt == REveProjection::kPT_RPhi);

   const TAxis* axis = isRPhi ? fData->GetPhiBins() :  fData->GetEtaBins();
   Int_t nBins = axis->GetNbins();

   Float_t min, max;
   if (isRPhi)
   {
      min = GetPhiMin() - fData->GetEps();
      max = GetPhiMax() + fData->GetEps();
      for (Int_t ibin = 1; ibin <= nBins; ++ibin) {
         clv = 0;
         if ( REveUtil::IsU1IntervalOverlappingByMinMax
              (min, max, axis->GetBinLowEdge(ibin), axis->GetBinUpEdge(ibin)))
         {
            clv = new REveCaloData::vCellId_t();
            fData->GetCellList(GetEta(), GetEtaRng(), axis->GetBinCenter(ibin), axis->GetBinWidth(ibin), *clv);
            if (!clv->size()) {
               delete clv; clv = 0;
            }
         }
         fCellLists.push_back(clv);
      }
   }
   else
   {
      min = GetEtaMin() - fData->GetEps();
      max = GetEtaMax() + fData->GetEps();
      for (Int_t ibin = 1; ibin <= nBins; ++ibin) {
         clv = 0;
         Float_t low = axis->GetBinLowEdge(ibin);
         Float_t up = axis->GetBinUpEdge(ibin) ;
         if (low >= min && up <= max)
         {
            clv = new REveCaloData::vCellId_t();
            fData->GetCellList(axis->GetBinCenter(ibin), axis->GetBinWidth(ibin), fPhi, GetPhiRng(), *clv);
            if (!clv->size()) {
               delete clv; clv = 0;
            }
         }
         fCellLists.push_back(clv);
      }
   }

   // cache max bin sum for auto scale
   if (!fScaleAbs)
   {
      fMaxESumBin  = 0;
      fMaxEtSumBin = 0;
      Float_t sumE  = 0;
      Float_t sumEt = 0;
      REveCaloData::CellData_t  cellData;
      for (Int_t ibin = 1; ibin <= nBins; ++ibin) {
         REveCaloData::vCellId_t* cids = fCellLists[ibin];
         if (cids)
         {
            sumE = 0; sumEt = 0;
            for (REveCaloData::vCellId_i it = cids->begin(); it != cids->end(); it++)
            {
               fData->GetCellData(*it, cellData);
               sumE  += cellData.Value(kFALSE);
               sumEt += cellData.Value(kTRUE);
            }
            fMaxESumBin  = TMath::Max(fMaxESumBin,  sumE);
            fMaxEtSumBin = TMath::Max(fMaxEtSumBin, sumEt);
         }
      }
      ComputeBBox();
   }

   fCellIdCacheOK= kTRUE;
}

//////////////////////////////////////////////s//////////////////////////////////
/// Sort selected cells in eta or phi bins.

void REveCalo2D::CellSelectionChangedInternal(REveCaloData::vCellId_t& inputCells, std::vector<REveCaloData::vCellId_t*>& outputCellLists)
{
   Bool_t isRPhi = (fManager->GetProjection()->GetType() == REveProjection::kPT_RPhi);
   const TAxis* axis = isRPhi ? fData->GetPhiBins() :  fData->GetEtaBins();

   // clear old cache
   for (vBinCells_i it = outputCellLists.begin(); it != outputCellLists.end(); it++)
   {
      if (*it)
      {
         (*it)->clear();
         delete *it;
      }
   }
   outputCellLists.clear();
   UInt_t nBins = axis->GetNbins();
   outputCellLists.resize(nBins+1);
   for (UInt_t b = 0; b <= nBins; ++b)
      outputCellLists[b] = 0;

   for(UInt_t bin = 1; bin <= nBins; ++bin)
   {
      REveCaloData::vCellId_t* idsInBin = fCellLists[bin];
      if (!idsInBin)
         continue;

      for (REveCaloData::vCellId_i i = idsInBin->begin(); i != idsInBin->end(); i++)
      {
         for (REveCaloData::vCellId_i j = inputCells.begin(); j != inputCells.end(); j++)
         {
            if( (*i).fTower == (*j).fTower && (*i).fSlice == (*j).fSlice)
            {
               if (!outputCellLists[bin])
                  outputCellLists[bin] = new REveCaloData::vCellId_t();

               outputCellLists[bin]->emplace_back((*i).fTower, (*i).fSlice, (*j).fFraction);
            }
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set absolute scale in projected calorimeter.

void REveCalo2D::SetScaleAbs(Bool_t sa)
{
   REveCaloViz::SetScaleAbs(sa);
   BuildCellIdCache();
}

////////////////////////////////////////////////////////////////////////////////
/// Virtual function of REveCaloViz.
/// Get transformation factor from E/Et to height.

Float_t REveCalo2D::GetValToHeight() const
{
   AssertCellIdCache();

   if (fScaleAbs)
   {
      return fMaxTowerH/fMaxValAbs;
   }
   else
   {
      if (fData->Empty())
         return 1;

      if (fPlotEt)
         return fMaxTowerH/fMaxEtSumBin;
      else
         return fMaxTowerH/fMaxESumBin;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Fill bounding-box information of the base-class TAttBBox (virtual method).
/// If member 'REveFrameBox* fFrame' is set, frame's corners are used as bbox.

void REveCalo2D::ComputeBBox()
{
   BBoxZero();

   Float_t x, y, z;
   Float_t th = fMaxTowerH                                           ;
   Float_t r  = fBarrelRadius + th;

   x = r,  y = 0, z = 0;
   fManager->GetProjection()->ProjectPoint(x, y, z, fDepth);
   BBoxCheckPoint(x, y, z);
   x = -r, y = 0, z = 0;
   fManager->GetProjection()->ProjectPoint(x, y, z, fDepth);
   BBoxCheckPoint(x, y, z);

   x = 0, y = 0, z = fEndCapPosF + th;
   fManager->GetProjection()->ProjectPoint(x, y, z, fDepth);
   BBoxCheckPoint(x, y, z);
   x = 0, y = 0, z = fEndCapPosB - th;
   fManager->GetProjection()->ProjectPoint(x, y, z, fDepth);
   BBoxCheckPoint(x, y, z);

   x = 0, y = r,  z = 0;
   fManager->GetProjection()->ProjectPoint(x, y, z, fDepth);
   BBoxCheckPoint(x, y, z);
   x = 0, y = -r, z = 0;
   fManager->GetProjection()->ProjectPoint(x, y, z, fDepth);
   BBoxCheckPoint(x, y, z);
}

////////////////////////////////////////////////////////////////////////////////
/// Fill core part of JSON representation.

Int_t REveCalo2D::WriteCoreJson(nlohmann::json &j, Int_t rnr_offset)
{
   Int_t ret = REveCaloViz::WriteCoreJson(j, rnr_offset);
   j["isRPhi"] = IsRPhi();
   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill core part of JSON representation for selection.

void REveCalo2D::WriteCoreJsonSelection(nlohmann::json &j, REveCaloData::vCellId_t cells)
{
   static const REveException eh("REveCalo2D::WriteCoreJsonSelection ");
   auto sarr = nlohmann::json::array();

   // selection
   // auto cellLists = isSel ? fCellListsSelected : fCellListsHighlighted;
   std::vector<REveCaloData::vCellId_t*> cellLists;
   CellSelectionChangedInternal(cells, cellLists);

   if (IsRPhi()) {
      REveCaloData::CellData_t cellData;
      Int_t  nSlices  = fData->GetNSlices();
      Float_t *sliceVal    = new Float_t[nSlices];
      Float_t *sliceValRef = new Float_t[nSlices];
      UInt_t nPhiBins = fData->GetPhiBins()->GetNbins();
      for(UInt_t phiBin = 1; phiBin <= nPhiBins; ++phiBin)
      {
         if (cellLists[phiBin])
         {
            if (!fCellLists[phiBin]) {
               delete[] sliceVal;
               delete[] sliceValRef;
               throw eh + "selected cell not in cell list cache.";
            }

            // selected eta sum
            for (Int_t s=0; s<nSlices; ++s) sliceVal[s] = 0;
            REveCaloData::vCellId_t& cids = *(cellLists[phiBin]);
            for (REveCaloData::vCellId_i i=cids.begin(); i!=cids.end(); i++) {
               fData->GetCellData((*i), cellData);
               sliceVal[i->fSlice] += cellData.Value(fPlotEt)*(*i).fFraction;
            }
            // referenced eta sum
            for (Int_t s=0; s<nSlices; ++s) sliceValRef[s] = 0;
            REveCaloData::vCellId_t& cidsRef = *(fCellLists[phiBin]);
            for (REveCaloData::vCellId_i i=cidsRef.begin(); i!=cidsRef.end(); i++) {
               fData->GetCellData(*i, cellData);
               sliceValRef[i->fSlice] += cellData.Value(fPlotEt)*(*i).fFraction;
            }

            // write
            for (Int_t s = 0; s < nSlices; ++s)  {
               if (sliceVal[s] > 0)
               {
                  nlohmann::json jsc;
                  jsc["b"] = phiBin;
                  jsc["s"] = s;
                  jsc["f"] = sliceVal[s]/sliceValRef[s];
                  sarr.push_back(jsc);
               }
            }
         }
      }
   }
   else {
      TAxis* axis        = fData->GetEtaBins();
      UInt_t nEtaBins    = axis->GetNbins();
      Int_t  nSlices     = fData->GetNSlices();

      Float_t *sliceValsUp     = new Float_t[nSlices];
      Float_t *sliceValsLow    = new Float_t[nSlices];
      Float_t *sliceValsUpRef  = new Float_t[nSlices];
      Float_t *sliceValsLowRef = new Float_t[nSlices];

      Float_t  towerH, towerHRef, offUp, offLow;
      REveCaloData::CellData_t cellData;

      for (UInt_t etaBin = 1; etaBin <= nEtaBins; ++etaBin)
      {
         if (cellLists[etaBin])
         {
            if (!fCellLists[etaBin]) {
               delete[] sliceValsUp;
               delete[] sliceValsLow;
               delete[] sliceValsUpRef;
               delete[] sliceValsLowRef;
               throw(eh + "selected cell not in cell list cache.");
            }
            offUp = 0; offLow =0;
            // selected phi sum
            for (Int_t s = 0; s < nSlices; ++s) {
               sliceValsUp[s] = 0; sliceValsLow[s] = 0;
            }
            REveCaloData::vCellId_t& cids = *(cellLists[etaBin]);
            for (REveCaloData::vCellId_i i=cids.begin(); i!=cids.end(); i++) {
               fData->GetCellData(*i, cellData);
               if (cellData.IsUpperRho())
                  sliceValsUp [i->fSlice] += cellData.Value(fPlotEt)*(*i).fFraction;
               else
                  sliceValsLow[i->fSlice] += cellData.Value(fPlotEt)*(*i).fFraction;
            }

            // reference phi sum
            for (Int_t s = 0; s < nSlices; ++s)
            {
               sliceValsUpRef[s] = 0; sliceValsLowRef[s] = 0;
            }
            REveCaloData::vCellId_t& cidsRef = *(fCellLists[etaBin]);
            for (REveCaloData::vCellId_i i=cidsRef.begin(); i!=cidsRef.end(); i++)
            {
               fData->GetCellData(*i, cellData);
               if (cellData.IsUpperRho())
                  sliceValsUpRef [i->fSlice] += cellData.Value(fPlotEt)*(*i).fFraction;
               else
                  sliceValsLowRef[i->fSlice] += cellData.Value(fPlotEt)*(*i).fFraction;
            }

            for (Int_t s = 0; s < nSlices; ++s)
            {
               //  phi +
               SetupHeight(sliceValsUpRef[s], s, towerHRef);
               if (sliceValsUp[s] > 0) {
                  SetupHeight(sliceValsUp[s], s, towerH);
                  nlohmann::json jsc;
                  jsc["b"] = etaBin;
                  jsc["s"] = s;
                  jsc["f"] = sliceValsUp[s]/sliceValsUpRef[s];
                  sarr.push_back(jsc);
               }
               offUp += towerHRef;

               // phi -
               SetupHeight(sliceValsLowRef[s], s, towerHRef);
               if (sliceValsLow[s] > 0) {
                  SetupHeight(sliceValsLow[s], s, towerH);
                  nlohmann::json jsc;
                  jsc["b"] = Int_t(-etaBin);
                  jsc["s"] = s;
                  jsc["f"] = sliceValsLow[s]/sliceValsLowRef[s];
                  sarr.push_back(jsc);
               }
               offLow += towerHRef;
            } // slices
         } // if eta bin
      } //eta bins
   } // RhoZ

   nlohmann::json rec = {};
   rec["caloVizId"] = GetElementId();
   rec["cells"] = sarr;
   j.push_back(rec);
}

////////////////////////////////////////////////////////////////////////////////
/// Creates 2D point array for rendering.

void REveCalo2D::BuildRenderData()
{
   AssertCellIdCache();
   bool isEmpty = fData->Empty();

   for (vBinCells_i it = fCellLists.begin(); it != fCellLists.end(); ++it)
   {
      if ((*it) && (*it)->empty())
      {
         isEmpty = false;
         break;
      }
   }
   if (isEmpty) return;

   fRenderData = std::make_unique<REveRenderData>("makeCalo2D");

   if (IsRPhi())
      BuildRenderDataRPhi();
   else
      BuildRenderDataRhoZ();
}


////////////////////////////////////////////////////////////////////////////////
/// Creates 2D point array in RhoZ projection.

void REveCalo2D::BuildRenderDataRhoZ()
{
   Int_t nSlices = fData->GetNSlices();

   REveCaloData::CellData_t cellData;
   Float_t *sliceValsUp  = new Float_t[nSlices];
   Float_t *sliceValsLow = new Float_t[nSlices];
   Bool_t   isBarrel;
   Float_t  towerH;
   Float_t transEtaF = GetTransitionEtaForward();
   Float_t transEtaB = GetTransitionEtaBackward();

   TAxis* axis = fData->GetEtaBins();
   UInt_t nEta = axis->GetNbins();
   Float_t pnts[12];
   for (UInt_t etaBin = 1; etaBin <= nEta; ++etaBin)
   {
      if (fCellLists[etaBin] )
      {
         assert(fCellLists[etaBin]);
         Float_t etaMin = axis->GetBinLowEdge(etaBin);
         Float_t etaMax = axis->GetBinUpEdge(etaBin);
         Float_t thetaMin = REveCaloData::EtaToTheta(etaMax);
         Float_t thetaMax = REveCaloData::EtaToTheta(etaMin);
         // printf("----------------------------------------- eta(%f, %f)\n", etaMin, etaMax);

         // clear
         Float_t offUp  = 0;
         Float_t offLow = 0;
         for (Int_t s = 0; s < nSlices; ++s) {
            sliceValsUp [s] = 0;
            sliceValsLow[s] = 0;
         }
         // values
         REveCaloData::vCellId_t* cids = fCellLists[etaBin];
         for (REveCaloData::vCellId_i it = cids->begin(); it != cids->end(); ++it)
         {
            fData->GetCellData(*it, cellData);
            if (cellData.IsUpperRho())
               sliceValsUp [it->fSlice] += cellData.Value(fPlotEt)*(*it).fFraction;
            else
               sliceValsLow[it->fSlice] += cellData.Value(fPlotEt)*(*it).fFraction;
         }

         isBarrel = !(etaMax > 0 && etaMax > transEtaF) && !(etaMin < 0 && etaMin < transEtaB);
         for (Int_t s = 0; s < nSlices; ++s)
         {
            //  phi +
            if (sliceValsUp[s])
            {
               SetupHeight(sliceValsUp[s], s, towerH);
               MakeRhoZCell(thetaMin, thetaMax, offUp, isBarrel, kTRUE , towerH, pnts);
               offUp += towerH;
               fRenderData->PushV(pnts, 12);
               fRenderData->PushI( s);
               fRenderData->PushI(etaBin);
               fRenderData->PushN(sliceValsUp[s]);
            }
            // phi -
            if (sliceValsLow[s])
            {
               SetupHeight(sliceValsLow[s], s, towerH);
               MakeRhoZCell(thetaMin, thetaMax, offLow, isBarrel, kFALSE , towerH, pnts);
               offLow += towerH;
               fRenderData->PushV(pnts, 12);
               fRenderData->PushI( s);
               fRenderData->PushI(etaBin);
               fRenderData->PushN(sliceValsLow[s]);
            }

         }
      }
   }

   delete [] sliceValsUp;
   delete [] sliceValsLow;
}

////////////////////////////////////////////////////////////////////////////////
/// Get cell vertices in RhoZ projection.
///
void REveCalo2D::MakeRhoZCell(Float_t thetaMin, Float_t thetaMax,
                              Float_t& offset, Bool_t isBarrel,  Bool_t phiPlus, Float_t towerH, float *pntsOut) const
{
   using namespace TMath;

   Float_t sin1 = Sin(thetaMin);
   Float_t cos1 = Cos(thetaMin);
   Float_t sin2 = Sin(thetaMax);
   Float_t cos2 = Cos(thetaMax);

   Float_t pnts[8];
   if (isBarrel)
   {
      Float_t r1 = fBarrelRadius/Abs(Sin(0.5f*(thetaMin+thetaMax))) + offset;
      Float_t r2 = r1 + towerH;

      pnts[0] = r1*sin1; pnts[1] = r1*cos1;
      pnts[2] = r2*sin1; pnts[3] = r2*cos1;
      pnts[4] = r2*sin2; pnts[5] = r2*cos2;
      pnts[6] = r1*sin2; pnts[7] = r1*cos2;
   }
   else
   {
      // endcap
      Float_t zE = GetForwardEndCapPos();
      // uses a different theta definition than GetTransitionThetaBackward(), so we need a conversion
      Float_t transThetaB = REveCaloData::EtaToTheta(GetTransitionEtaBackward());
      if (thetaMax >= transThetaB)
         zE = Abs(GetBackwardEndCapPos());
      Float_t r1 = zE/Abs(Cos(0.5f*(thetaMin+thetaMax))) + offset;
      Float_t r2 = r1 + towerH;

      pnts[0] = r1*sin1; pnts[1] = r1*cos1;
      pnts[2] = r2*sin1; pnts[3] = r2*cos1;
      pnts[4] = r2*sin2; pnts[5] = r2*cos2;
      pnts[6] = r1*sin2; pnts[7] = r1*cos2;
   }


   Float_t x, y, z;
   for (Int_t i = 0; i < 4; ++i)
   {
      x = 0.f;
      y = phiPlus ? Abs(pnts[2*i]) : -Abs(pnts[2*i]);
      z = pnts[2*i+1];
      fManager->GetProjection()->ProjectPoint(x, y, z, fDepth);

      int j = phiPlus ? i : (3 -i);
      pntsOut[j*3] = x;
      pntsOut[j*3 + 1] = y;
      pntsOut[j*3 + 2] = z;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Creates 2D point array in RPhi() projection.

void REveCalo2D::BuildRenderDataRPhi()
{
   REveCaloData* data = fData;
   Int_t    nSlices  = data->GetNSlices();
   Float_t *sliceVal = new Float_t[nSlices];
   REveCaloData::CellData_t cellData;
   Float_t towerH;

   UInt_t nPhi = data->GetPhiBins()->GetNbins();
   TAxis* axis = data->GetPhiBins();
   float pnts[12];
   for(UInt_t phiBin = 1; phiBin <= nPhi; ++phiBin)
   {
      if (fCellLists[phiBin] )
      {
         // reset values
         Float_t off = 0;
         for (Int_t s=0; s<nSlices; ++s)
            sliceVal[s] = 0;

         // sum eta cells
         REveCaloData::vCellId_t* cids = fCellLists[phiBin];
         for (REveCaloData::vCellId_i it = cids->begin(); it != cids->end(); it++)
         {
            data->GetCellData(*it, cellData);
            sliceVal[(*it).fSlice] += cellData.Value(fPlotEt)*(*it).fFraction;
         }
         for (Int_t s = 0; s < nSlices; ++s)
         {
            SetupHeight(sliceVal[s], s, towerH);
            MakeRPhiCell(axis->GetBinLowEdge(phiBin), axis->GetBinUpEdge(phiBin), towerH, off, pnts);
            fRenderData->PushV(pnts, 12);
            fRenderData->PushI(s);
            fRenderData->PushI(phiBin);
            fRenderData->PushN(sliceVal[s]);
            off += towerH;
         }
      }
   }

   delete [] sliceVal;
}

////////////////////////////////////////////////////////////////////////////////
/// Calculate vertices for the calorimeter cell in RPhi projection.
/// Returns outside radius of the tower.

void REveCalo2D::MakeRPhiCell(Float_t phiMin, Float_t phiMax,
                                Float_t towerH, Float_t offset, float* pntsOut) const
{
   using namespace TMath;

   Float_t r1 = fBarrelRadius + offset;
   Float_t r2 = r1 + towerH;

   Float_t pnts[8];
   pnts[0] = r1*Cos(phiMin); pnts[1] = r1*Sin(phiMin);
   pnts[2] = r2*Cos(phiMin); pnts[3] = r2*Sin(phiMin);
   pnts[4] = r2*Cos(phiMax); pnts[5] = r2*Sin(phiMax);
   pnts[6] = r1*Cos(phiMax); pnts[7] = r1*Sin(phiMax);

   Float_t x, y, z;
   for (Int_t i = 0; i < 4; ++i)
   {
      pntsOut[i*3]   = pnts[2*i];
      pntsOut[i*3+1] = pnts[2*i+1];
      pntsOut[i*3+2] = 0.f;
      fManager->GetProjection()->ProjectPoint(x, y, z, fDepth);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Client callback

void REveCalo2D::NewBinPicked(Int_t bin, Int_t slice, Int_t selectionId, bool multi)
{
   bool is_upper = bin >= 0;
   bin = abs(bin);

   REveCaloData::vCellId_t sel;
   for (REveCaloData::vCellId_i it = fCellLists[bin]->begin(); it != fCellLists[bin]->end(); ++it)
   {
      if ((*it).fSlice == slice)
      {
         if (IsRPhi())
         {
            sel.push_back(*it);
         }
         else
         {
            REveCaloData::CellData_t cd;
            fData->GetCellData(*it, cd);
            if ((is_upper && cd.IsUpperRho()) || (!is_upper && !cd.IsUpperRho()))
               sel.push_back(*it);
         }
      }
   }
   fData->ProcessSelection(sel, selectionId, multi);
}


/** \class REveCaloLego
\ingroup REve
Visualization of calorimeter data as eta/phi histogram.
*/

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

REveCaloLego::REveCaloLego(REveCaloData* d, const char* n, const char* t):
   REveCaloViz(d, n, t),

   fFontColor(-1),
   fGridColor(-1),
   fPlaneColor(kRed-5),
   fPlaneTransparency(60),

   fNZSteps(6),
   fZAxisStep(0.f),

   fAutoRebin(kTRUE),

   fPixelsPerBin(12),
   fNormalizeRebin(kFALSE),

   fProjection(kAuto),
   f2DMode(kValSize),
   fBoxMode(kBack),

   fDrawHPlane(kFALSE),
   fHPlaneVal(0),

   fHasFixedHeightIn2DMode(kFALSE),
   fFixedHeightValIn2DMode(0.f),

   fDrawNumberCellPixels(18), // draw numbers on cell above 30 pixels
   fCellPixelFontSize(12) // size of cell fonts in pixels
{
   fMaxTowerH = 4;
   SetNameTitle("REveCaloLego", "REveCaloLego");
}

////////////////////////////////////////////////////////////////////////////////
// Set data.

void REveCaloLego::SetData(REveCaloData* data)
{
   REveCaloViz::SetData(data);
}

////////////////////////////////////////////////////////////////////////////////
/// Build list of drawn cell IDs. For more information see REveCaloLegoGL:DirectDraw().

void REveCaloLego::BuildCellIdCache()
{
   fCellList.clear();

   fData->GetCellList(GetEta(), GetEtaRng(), GetPhi(), GetPhiRng(), fCellList);
   fCellIdCacheOK = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill bounding-box information of the base-class TAttBBox (virtual method).
/// If member 'REveFrameBox* fFrame' is set, frame's corners are used as bbox.

void REveCaloLego::ComputeBBox()
{

   // fBBox = Float_t[6] X(min,max), Y(min,max), Z(min,max)

   BBoxZero();

   Float_t ex = 1.2; // 20% offset for axis labels

   Float_t a = 0.5*ex;

   fBBox[0] = -a;
   fBBox[1] =  a;
   fBBox[2] = -a;
   fBBox[3] =  a;

   // scaling is relative to shortest XY axis
   Double_t em, eM, pm, pM;
   fData->GetEtaLimits(em, eM);
   fData->GetPhiLimits(pm, pM);
   Double_t r = (eM-em)/(pM-pm);
   if (r<1)
   {
      fBBox[2] /= r;
      fBBox[3] /= r;
   }
   else
   {
      fBBox[0] *= r;
      fBBox[1] *= r;
   }

   fBBox[4] =  0;
   if (fScaleAbs && !fData->Empty())
      fBBox[5] = GetMaxVal()*GetValToHeight();
   else
      fBBox[5] = fMaxTowerH;
}

