// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "ROOT/REveCaloData.hxx"
#include "ROOT/REveCalo.hxx"
#include "ROOT/REveUtil.hxx"
#include "ROOT/REveManager.hxx"
#include "ROOT/REveSelection.hxx"


// #include "TGLSelectRecord.h"

#include "TAxis.h"
#include "THStack.h"
#include "TH2.h"
#include "TMath.h"
#include "TList.h"

#include <cassert>
#include <algorithm>
#include <set>

using namespace ROOT::Experimental;

/** \class  REveCaloData::CellGeom_t
\ingroup REve
Cell geometry inner structure.
*/

////////////////////////////////////////////////////////////////////////////////
/// Print member data.

void REveCaloData::CellGeom_t::Dump() const
{
   printf("%f, %f %f, %f \n", fEtaMin, fEtaMax, fPhiMin, fPhiMax);
}

////////////////////////////////////////////////////////////////////////////////

void REveCaloData::CellGeom_t::Configure(Float_t etaMin, Float_t etaMax, Float_t phiMin, Float_t phiMax)
{
   fEtaMin = etaMin;
   fEtaMax = etaMax;

   fPhiMin = phiMin;
   fPhiMax = phiMax;

   // Complain if phi is out of [-2*pi, 2*pi] range.
   if (fPhiMin < - TMath::TwoPi() || fPhiMin > TMath::TwoPi() ||
       fPhiMax < - TMath::TwoPi() || fPhiMax > TMath::TwoPi())
   {
     ::Error("REveCaloData::CellGeom_t::Configure", "phiMin and phiMax should be between -2*pi and 2*pi (min=%f, max=%f). RhoZ projection will be wrong.",
             fPhiMin, fPhiMax);
   }

   fThetaMin = EtaToTheta(fEtaMax);
   fThetaMax = EtaToTheta(fEtaMin);
}

/** \class REveCaloData::CellData_t
\ingroup REve
Cell data inner structure.
*/

////////////////////////////////////////////////////////////////////////////////
/// Return energy value associated with the cell, usually Et.
/// If isEt is false it is transformed into energy E.

Float_t REveCaloData::CellData_t::Value(Bool_t isEt) const
{
   if (isEt)
      return fValue;
   else
      return TMath::Abs(fValue/TMath::Sin(Theta()));
}

////////////////////////////////////////////////////////////////////////////////
/// Print member data.

void REveCaloData::CellData_t::Dump() const
{
   printf("%f, %f %f, %f \n", fEtaMin, fEtaMax, fPhiMin, fPhiMax);
}

////////////////////////////////////////////////////////////////////////////////

Float_t* REveCaloData::RebinData_t::GetSliceVals(Int_t bin)
{
   //   printf("get val vec bin %d size %d\n", bin, fBinData.size());
   if (fBinData[bin] == -1)
   {
      fBinData[bin] = fSliceData.size();

      for (Int_t i=0; i<fNSlices; i++)
         fSliceData.push_back(0.f);
   }

   return &fSliceData[fBinData[bin]];
}

/** \class REveCaloData
\ingroup REve
A central manager for calorimeter event data. It provides a list of
cells within requested phi and eta range.
*/


////////////////////////////////////////////////////////////////////////////////

REveCaloData::REveCaloData(const char* n, const char* t):
   REveElement(),

   fWrapTwoPi(kTRUE),

   fMaxValEt(0),
   fMaxValE(0),

   fEps(0)
{
   SetNameTitle(n,t);
   // Constructor.
}

////////////////////////////////////////////////////////////////////////////////
/// Process newly selected cells with given select-record.

void REveCaloData::ProcessSelection(vCellId_t& sel_cells, UInt_t selectionId, Bool_t multiple)
{
   if (fSelector)
   {
      fSelector->ProcessSelection(sel_cells, selectionId, multiple);
   }
   else
   {
      REveSelection* selection = dynamic_cast<REveSelection*> (ROOT::Experimental::gEve->FindElementById(selectionId));

      std::set<int> secondary_idcs;
      for (vCellId_i i = sel_cells.begin(); i != sel_cells.end(); ++i)
      {
         int id = (i->fSlice << 24) + i->fTower;
         secondary_idcs.insert(id);
      }
      selection->NewElementPicked(GetElementId(), multiple, true, secondary_idcs);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Populate set impSelSet with derived / dependant elements.
///

void REveCaloData::FillImpliedSelectedSet(Set_t& impSelSet)
{
   // printf("REveCaloData::FillImpliedSelectedSet\n");
   for (auto &n : fNieces)
   {
      impSelSet.insert(n);
   }
}
////////////////////////////////////////////////////////////////////////////////
/// Set threshold for given slice.

void REveCaloData::SetSliceThreshold(Int_t slice, Float_t val)
{
   fSliceInfos[slice].fThreshold = val;
   InvalidateUsersCellIdCache();
}

////////////////////////////////////////////////////////////////////////////////
/// Get threshold for given slice.

Float_t REveCaloData::GetSliceThreshold(Int_t slice) const
{
   return fSliceInfos[slice].fThreshold;
}

////////////////////////////////////////////////////////////////////////////////
/// Set color for given slice.

void REveCaloData::SetSliceColor(Int_t slice, Color_t col)
{
   fSliceInfos[slice].fColor = col;
   for (auto &c : fNieces)
   {
      c->AddStamp(REveElement::kCBObjProps);
   }
   AddStamp(REveElement::kCBObjProps);
}

////////////////////////////////////////////////////////////////////////////////
/// Get color for given slice.

Color_t REveCaloData::GetSliceColor(Int_t slice) const
{
   return fSliceInfos[slice].fColor;
}

////////////////////////////////////////////////////////////////////////////////
/// Set transparency for given slice.

void REveCaloData::SetSliceTransparency(Int_t slice, Char_t t)
{
   fSliceInfos[slice].fTransparency = t;
   for (auto &c : fNieces)
   {
      c->AddStamp(REveElement::kCBObjProps);
   }
   AddStamp(REveElement::kCBObjProps);
}

////////////////////////////////////////////////////////////////////////////////
/// Get transparency for given slice.

Char_t REveCaloData::GetSliceTransparency(Int_t slice) const
{
   return fSliceInfos[slice].fTransparency;
}

////////////////////////////////////////////////////////////////////////////////
/// Invalidate cell ids cache on back ptr references.

void REveCaloData::InvalidateUsersCellIdCache()
{
   REveCaloViz* calo;
   for (auto &c : fNieces)
   {
      calo = dynamic_cast<REveCaloViz*>(c);
      calo->InvalidateCellIdCache();
      calo->StampObjProps();
   }
   AddStamp(REveElement::kCBObjProps);
}

////////////////////////////////////////////////////////////////////////////////
/// Tell users (REveCaloViz instances using this data) that data
/// has changed and they should update the limits/scales etc.
/// This is done by calling REveCaloViz::DataChanged().

void REveCaloData::DataChanged()
{
   REveCaloViz* calo;
   for (auto &c : fNieces)
   {
      calo = dynamic_cast<REveCaloViz*>(c);
      calo->DataChanged();
      calo->StampObjProps();
   }
}


////////////////////////////////////////////////////////////////////////////////

Int_t REveCaloData::WriteCoreJson(nlohmann::json &j, Int_t rnr_offset)
{
   Int_t ret = REveElement::WriteCoreJson(j, rnr_offset);

   auto sarr = nlohmann::json::array();
   for (auto &s : fSliceInfos)
   {
      nlohmann::json slice = {};
      slice["name"]     = s.fName;
      slice["threshold"] = s.fThreshold;
      slice["color"]    = s.fColor;
      sarr.push_back(slice);
   }
   j["sliceInfos"] = sarr;
   return ret;
}
////////////////////////////////////////////////////////////////////////////////

void  REveCaloData::FillExtraSelectionData(nlohmann::json& j, const std::set<int>& secondary_idcs) const
{
   vCellId_t cells;

   if (fSelector)
   {
      fSelector->GetCellsFromSecondaryIndices(secondary_idcs, cells);
   }
   else
   {
      for (auto &id : secondary_idcs ) {

         int s = (id >> 24);
         int t = id & 0xffffff;
         REveCaloData::CellId_t cell(t, s, 1.0f);
         cells.push_back(cell);
      }
   }
   for (auto &c : fNieces)
      ((REveCaloViz*)c)->WriteCoreJsonSelection(j, cells);
}

////////////////////////////////////////////////////////////////////////////////

Float_t REveCaloData::EtaToTheta(Float_t eta)
{
   using namespace TMath;

   if (eta < 0)
      return Pi() - 2*ATan(Exp(- Abs(eta)));
   else
      return 2*ATan(Exp(- Abs(eta)));
}

////////////////////////////////////////////////////////////////////////////////
std::string REveCaloData::GetHighlightTooltip(const std::set<int>& secondary_idcs) const
{
   std::string s;
   CellData_t cellData;

   Bool_t single = secondary_idcs.size() == 1;
   Float_t sum = 0;


   for (auto &id : secondary_idcs ) {

      int slice = (id >> 24);
      int tower = id & 0xffffff;
      REveCaloData::CellId_t cell(tower, slice, 1.0f);
      GetCellData(cell, cellData);

      s += TString::Format("%s %.2f (%.3f, %.3f)",
                           fSliceInfos[slice].fName.Data(), cellData.fValue,
                           cellData.Eta(), cellData.Phi());

      if (single) return s;
      s += "\n";
      sum += cellData.fValue;
   }
   s += TString::Format("Sum = %.2f", sum);
   return s;
}

/** \class  REveCaloDataVec
\ingroup REve
Calo data for universal cell geometry.
*/

////////////////////////////////////////////////////////////////////////////////

REveCaloDataVec::REveCaloDataVec(Int_t nslices):
   REveCaloData(),

   fTower(0),
   fEtaMin( 1e3),
   fEtaMax(-1e3),
   fPhiMin( 1e3),
   fPhiMax(-1e3)
{
   // Constructor.

   fSliceInfos.assign(nslices, SliceInfo_t());

   fSliceVec.assign(nslices, std::vector<Float_t> ());
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

REveCaloDataVec::~REveCaloDataVec()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Add new slice.

Int_t REveCaloDataVec::AddSlice()
{
  fSliceInfos.push_back(SliceInfo_t());
  fSliceVec.push_back(std::vector<Float_t> ());
  fSliceVec.back().resize(fGeomVec.size(), 0.f);

  return fSliceInfos.size() - 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Add tower within eta/phi range.

Int_t REveCaloDataVec::AddTower(Float_t etaMin, Float_t etaMax, Float_t phiMin, Float_t phiMax)
{
   assert (etaMin < etaMax);
   assert (phiMin < phiMax);

   fGeomVec.push_back(CellGeom_t(etaMin, etaMax, phiMin, phiMax));

   for (vvFloat_i it=fSliceVec.begin(); it!=fSliceVec.end(); ++it)
      (*it).push_back(0);

   if (etaMin < fEtaMin) fEtaMin = etaMin;
   if (etaMax > fEtaMax) fEtaMax = etaMax;

   if (phiMin < fPhiMin) fPhiMin = phiMin;
   if (phiMax > fPhiMax) fPhiMax = phiMax;

   fTower = fGeomVec.size() - 1;
   return fTower;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill given slice in the current tower.

void REveCaloDataVec::FillSlice(Int_t slice, Float_t val)
{
   fSliceVec[slice][fTower] = val;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill given slice in a given tower.

void REveCaloDataVec::FillSlice(Int_t slice, Int_t tower, Float_t val)
{
   fSliceVec[slice][tower] = val;
}


////////////////////////////////////////////////////////////////////////////////
/// Get list of cell-ids for given eta/phi range.

void REveCaloDataVec::GetCellList(Float_t eta, Float_t etaD,
                                  Float_t phi, Float_t phiD,
                                  REveCaloData::vCellId_t &out) const
{
   using namespace TMath;

   Float_t etaMin = eta - etaD*0.5;
   Float_t etaMax = eta + etaD*0.5;

   Float_t phiMin = phi - phiD*0.5;
   Float_t phiMax = phi + phiD*0.5;

   Int_t nS = fSliceVec.size();

   Int_t tower = 0;
   Float_t fracx=0, fracy=0, frac;
   Float_t minQ, maxQ;

   for(vCellGeom_ci i=fGeomVec.begin(); i!=fGeomVec.end(); i++)
   {
      const CellGeom_t &cg = *i;
      fracx = REveUtil::GetFraction(etaMin, etaMax, cg.fEtaMin, cg.fEtaMax);
      if (fracx > 1e-3)
      {
         minQ = cg.fPhiMin;
         maxQ = cg.fPhiMax;

         if (fWrapTwoPi)
         {
            if (maxQ < phiMin)
            {
               minQ += TwoPi(); maxQ += TwoPi();
            }
            else if (minQ > phiMax)
            {
               minQ -= TwoPi(); maxQ -= TwoPi();
            }
         }

         if (maxQ >= phiMin && minQ <= phiMax)
         {
            fracy = REveUtil::GetFraction(phiMin, phiMax, minQ, maxQ);
            if (fracy > 1e-3)
            {
               frac = fracx*fracy;
               for (Int_t s=0; s<nS; s++)
               {
                  if (fSliceVec[s][tower] > fSliceInfos[s].fThreshold)
                     out.push_back(CellId_t(tower, s, frac));
               }
            }
         }
      }
      tower++;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Rebin cells.

void REveCaloDataVec::Rebin(TAxis* ax, TAxis* ay, vCellId_t &ids, Bool_t et, RebinData_t& rdata) const
{
   rdata.fNSlices = GetNSlices();
   rdata.fBinData.assign((ax->GetNbins()+2)*(ay->GetNbins()+2), -1);

   CellData_t cd;
   for (vCellId_i it = ids.begin(); it != ids.end(); ++it)
   {
      GetCellData(*it, cd);
      Int_t iMin = ax->FindBin(cd.EtaMin());
      Int_t iMax = ax->FindBin(cd.EtaMax());
      Int_t jMin = ay->FindBin(cd.PhiMin());
      Int_t jMax = ay->FindBin(cd.PhiMax());
      for (Int_t i = iMin; i <= iMax; ++i)
      {
         if (i < 0 || i > ax->GetNbins()) continue;
         for (Int_t j = jMin; j <= jMax; ++j)
         {
            if (j < 0 || j > ay->GetNbins()) continue;

            Double_t ratio = REveUtil::GetFraction(ax->GetBinLowEdge(i), ax->GetBinUpEdge(i), cd.EtaMin(), cd.EtaMax())
                           * REveUtil::GetFraction(ay->GetBinLowEdge(j), ay->GetBinUpEdge(j), cd.PhiMin(), cd.PhiMax());

            if (ratio > 1e-6f)
            {
               Float_t* slices = rdata.GetSliceVals(i + j*(ax->GetNbins()+2));
               slices[(*it).fSlice] += ratio * cd.Value(et);
            }
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get cell geometry and value from cell ID.

void REveCaloDataVec::GetCellData(const REveCaloData::CellId_t &id,
                                  REveCaloData::CellData_t& cellData) const
{
   cellData.CellGeom_t::operator=( fGeomVec[id.fTower] );
   cellData.fValue = fSliceVec[id.fSlice][id.fTower];
}

////////////////////////////////////////////////////////////////////////////////
/// Update limits and notify data users.

void REveCaloDataVec::DataChanged()
{
   using namespace TMath;

   // update max E/Et values

   fMaxValE = 0;
   fMaxValEt = 0;
   Float_t sum=0;
   //   printf("geom vec %d slices %d\n",fGeomVec.size(), fSliceVec.size() );

   for (UInt_t tw=0; tw<fGeomVec.size(); tw++)
   {
      sum=0;
      for (vvFloat_i it=fSliceVec.begin(); it!=fSliceVec.end(); ++it)
         sum += (*it)[tw];

      if (sum > fMaxValEt ) fMaxValEt=sum;

      sum /= Abs(Sin(EtaToTheta(fGeomVec[tw].Eta())));

      if (sum > fMaxValE) fMaxValE=sum;
   }

   REveCaloData::DataChanged();
}


////////////////////////////////////////////////////////////////////////////////
/// Set XY axis from cells geometry.

void  REveCaloDataVec::SetAxisFromBins(Double_t epsX, Double_t epsY)
{
   std::vector<Double_t> binX;
   std::vector<Double_t> binY;

   for(vCellGeom_ci i=fGeomVec.begin(); i!=fGeomVec.end(); i++)
   {
      const CellGeom_t &ch = *i;

      binX.push_back(ch.EtaMin());
      binX.push_back(ch.EtaMax());
      binY.push_back(ch.PhiMin());
      binY.push_back(ch.PhiMax());
   }

   std::sort(binX.begin(), binX.end());
   std::sort(binY.begin(), binY.end());

   Int_t cnt = 0;
   Double_t sum = 0;
   Double_t val;

   // X axis
   Double_t dx = binX.back() - binX.front();
   epsX *= dx;
   std::vector<Double_t> newX;
   newX.push_back(binX.front()); // underflow
   Int_t nX = binX.size()-1;
   for(Int_t i=0; i<nX; i++)
   {
      val = (sum +binX[i])/(cnt+1);
      if (binX[i+1] -val > epsX)
      {
         newX.push_back(val);
         cnt = 0;
         sum = 0;
      }
      else
      {
         sum += binX[i];
         cnt++;
      }
   }
   newX.push_back(binX.back()); // overflow

   // Y axis
   cnt = 0;
   sum = 0;
   std::vector<Double_t> newY;
   Double_t dy = binY.back() - binY.front();
   epsY *= dy;
   newY.push_back(binY.front());// underflow
   Int_t nY = binY.size()-1;
   for(Int_t i=0 ; i<nY; i++)
   {
      val = (sum +binY[i])/(cnt+1);
      if (binY[i+1] -val > epsY )
      {
         newY.push_back(val);
         cnt = 0;
         sum = 0;
      }
      else
      {
         sum += binY[i];
         cnt++;
      }

   }
   newY.push_back(binY.back()); // overflow
   fEtaAxis = std::make_unique<TAxis>(newX.size()-1, &newX[0]);
   fEtaAxis = std::make_unique<TAxis>(newY.size()-1, &newY[0]);
   fEtaAxis->SetNdivisions(510);
   fPhiAxis->SetNdivisions(510);
}

/** \class REveCaloDataHist
\ingroup REve
A central manager for calorimeter data of an event written in TH2F.
X axis is used for eta and Y axis for phi.
*/


////////////////////////////////////////////////////////////////////////////////
/// Constructor.

REveCaloDataHist::REveCaloDataHist():
   REveCaloData(),

   fHStack(0)
{
   fHStack = new THStack();
   fEps    = 1e-5;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

REveCaloDataHist::~REveCaloDataHist()
{
   delete fHStack;
}

////////////////////////////////////////////////////////////////////////////////
/// Update limits and notify data users.

void REveCaloDataHist::DataChanged()
{
   using namespace TMath;

   // update max E/Et values
   fMaxValE  = 0;
   fMaxValEt = 0;

   if (GetNSlices() < 1) return;

   TH2* hist = GetHist(0);
   for (Int_t ieta = 1; ieta <= fEtaAxis->GetNbins(); ++ieta)
   {
      Double_t eta = fEtaAxis->GetBinCenter(ieta); // conversion E/Et
      for (Int_t iphi = 1; iphi <= fPhiAxis->GetNbins(); ++iphi)
      {
         Double_t value = 0;
         for (Int_t i = 0; i < GetNSlices(); ++i)
         {
            hist = GetHist(i);
            Int_t bin = hist->GetBin(ieta, iphi);
            value += hist->GetBinContent(bin);
         }

         if (value > fMaxValEt ) fMaxValEt = value;

         value /= Abs(Sin(EtaToTheta(eta)));

         if (value > fMaxValE) fMaxValE = value;
      }
   }
   REveCaloData::DataChanged();
}

////////////////////////////////////////////////////////////////////////////////
/// Get list of cell IDs in given eta and phi range.

void REveCaloDataHist::GetCellList(Float_t eta, Float_t etaD,
                                   Float_t phi, Float_t phiD,
                                   REveCaloData::vCellId_t &out) const
{
   using namespace TMath;

   Float_t etaMin = eta - etaD*0.5 -fEps;
   Float_t etaMax = eta + etaD*0.5 +fEps;

   Float_t phiMin = phi - phiD*0.5 -fEps;
   Float_t phiMax = phi + phiD*0.5 +fEps;

   Int_t nEta = fEtaAxis->GetNbins();
   Int_t nPhi = fPhiAxis->GetNbins();
   Int_t nSlices = GetNSlices();

   Int_t bin  = 0;

   Bool_t accept;
   for (Int_t ieta = 1; ieta <= nEta; ++ieta)
   {
      if (fEtaAxis->GetBinLowEdge(ieta) >= etaMin && fEtaAxis->GetBinUpEdge(ieta) <= etaMax)
      {
         for (Int_t iphi = 1; iphi <= nPhi; ++iphi)
         {
            if (fWrapTwoPi )
            {
               accept = REveUtil::IsU1IntervalContainedByMinMax
                  (phiMin, phiMax, fPhiAxis->GetBinLowEdge(iphi), fPhiAxis->GetBinUpEdge(iphi));
            }
            else
            {
               accept = fPhiAxis->GetBinLowEdge(iphi) >= phiMin &&  fPhiAxis->GetBinUpEdge(iphi) <= phiMax &&
                  fPhiAxis->GetBinLowEdge(iphi) >= phiMin &&  fPhiAxis->GetBinUpEdge(iphi) <= phiMax;
            }

            if (accept)
            {
               for (Int_t s = 0; s < nSlices; ++s)
               {
                  TH2F *hist = GetHist(s);
                  bin = hist->GetBin(ieta, iphi);
                  if (hist->GetBinContent(bin) > fSliceInfos[s].fThreshold)
                     out.push_back(REveCaloData::CellId_t(bin, s));
               } // hist slices
            }
         } // phi bins
      }
   } // eta bins
}

////////////////////////////////////////////////////////////////////////////////
/// Rebin

void REveCaloDataHist::Rebin(TAxis* ax, TAxis* ay, REveCaloData::vCellId_t &ids, Bool_t et, RebinData_t &rdata) const
{
   rdata.fNSlices = GetNSlices();
   rdata.fBinData.assign((ax->GetNbins()+2)*(ay->GetNbins()+2), -1);
   REveCaloData::CellData_t cd;
   Float_t *val;
   Int_t i, j, w;
   Int_t binx, biny;
   Int_t bin;

   for (vCellId_i it=ids.begin(); it!=ids.end(); ++it)
   {
      GetCellData(*it, cd);
      GetHist(it->fSlice)->GetBinXYZ((*it).fTower, i, j, w);
      binx = ax->FindBin(fEtaAxis->GetBinCenter(i));
      biny = ay->FindBin(fPhiAxis->GetBinCenter(j));
      bin = biny*(ax->GetNbins()+2)+binx;
      val = rdata.GetSliceVals(bin);
      Double_t ratio = REveUtil::GetFraction(ax->GetBinLowEdge(binx), ax->GetBinUpEdge(binx), cd.EtaMin(), cd.EtaMax())
                     * REveUtil::GetFraction(ay->GetBinLowEdge(biny), ay->GetBinUpEdge(biny), cd.PhiMin(), cd.PhiMax());

      val[(*it).fSlice] += cd.Value(et)*ratio;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get cell geometry and value from cell ID.

void REveCaloDataHist::GetCellData(const REveCaloData::CellId_t &id,
                                   REveCaloData::CellData_t& cellData) const
{
   TH2F* hist = GetHist(id.fSlice);

   Int_t x, y, z;
   hist->GetBinXYZ(id.fTower, x, y, z);

   cellData.fValue =  hist->GetBinContent(id.fTower);
   cellData.Configure(hist->GetXaxis()->GetBinLowEdge(x),
                      hist->GetXaxis()->GetBinUpEdge(x),
                      hist->GetYaxis()->GetBinLowEdge(y),
                      hist->GetYaxis()->GetBinUpEdge(y));
}

////////////////////////////////////////////////////////////////////////////////
/// Add new slice to calo tower. Updates cached variables fMaxValE
/// and fMaxValEt
/// Return last index in the vector of slice infos.

Int_t REveCaloDataHist::AddHistogram(TH2F* hist)
{
   if (!fEtaAxis) {
      fEtaAxis.reset(new TAxis(*hist->GetXaxis()));
      fPhiAxis.reset(new TAxis(*hist->GetYaxis()));
   }
   fHStack->Add(hist);
   fSliceInfos.push_back(SliceInfo_t());
   fSliceInfos.back().fName  = hist->GetName();
   fSliceInfos.back().fColor = hist->GetLineColor();

   DataChanged();

   return fSliceInfos.size() - 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Get histogram in given slice.

TH2F* REveCaloDataHist::GetHist(Int_t slice) const
{
   assert(slice >= 0 && slice < fHStack->GetHists()->GetSize());
   return (TH2F*) fHStack->GetHists()->At(slice);
}

////////////////////////////////////////////////////////////////////////////////
/// Get eta limits.

void REveCaloDataHist::GetEtaLimits(Double_t &min, Double_t &max) const
{
   min = fEtaAxis->GetXmin();
   max = fEtaAxis->GetXmax();
}

////////////////////////////////////////////////////////////////////////////////
/// Get phi limits.

void REveCaloDataHist::GetPhiLimits(Double_t &min, Double_t &max) const
{
   min = fPhiAxis->GetXmin();
   max = fPhiAxis->GetXmax();
}

////////////////////////////////////////////////////////////////////////////////
/// Process selection. Called from REveCaloViz object

void REveCaloDataSelector::ProcessSelection(REveCaloData::vCellId_t& sel_cells, UInt_t selectionId, Bool_t multi)
{
   // only one slice can be user selected at once
   fActiveSlice  = sel_cells.front().fSlice;
   for (auto &si : fSliceSelectors)
   {
      if (si->GetSliceIndex() == fActiveSlice) {
         si->ProcessSelection(sel_cells, selectionId, multi);
         break;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// GetCellsFromSecondaryIndices used in implied selection

void REveCaloDataSelector::GetCellsFromSecondaryIndices(const std::set<int>& idcs, REveCaloData::vCellId_t& out)
{
   for (auto &si : fSliceSelectors)
   {
      if (si->GetSliceIndex() == fActiveSlice) {
         si->GetCellsFromSecondaryIndices(idcs, out);
         break;
      }
   }
}
