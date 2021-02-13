// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TEveCaloData
#define ROOT_TEveCaloData

#include <vector>
#include "TEveElement.h"

#include "TMath.h"

class TGLSelectRecord;

class TH2F;
class TAxis;
class THStack;

class TEveCaloData: public TEveElement,
                    public TNamed
{
public:
   struct SliceInfo_t
   {
      TString  fName;      // Name of the slice, eg. ECAL, HCAL.
      Float_t  fThreshold; // Only display towers with higher energy.
      Color_t  fColor;     // Color used to draw this longitudinal slice.
      Color_t  fTransparency; // Transparency used to draw this longitudinal slice.

      SliceInfo_t(): fName(""), fThreshold(0), fColor(kRed), fTransparency(0) {}

      virtual ~SliceInfo_t() {}

      void Setup(const char* name, Float_t threshold, Color_t col, Char_t transp = 101)
      {
         fName      = name;
         fThreshold = threshold;
         fColor     = col;
         if (transp <= 100) fTransparency = transp;
      };

      ClassDef(SliceInfo_t, 0); // Slice info for histogram stack.
   };

   typedef std::vector<SliceInfo_t>           vSliceInfo_t;
   typedef std::vector<SliceInfo_t>::iterator vSliceInfo_i;

   /**************************************************************************/

   struct CellId_t
   {
      // Cell ID inner structure.

      Int_t fTower;
      Int_t fSlice;

      Float_t fFraction;

      CellId_t(Int_t t, Int_t s, Float_t f=1.0f) : fTower(t), fSlice(s), fFraction(f) {}

      bool operator<(const CellId_t& o) const
      { return (fTower == o.fTower) ? fSlice < o.fSlice : fTower < o.fTower; }
   };

   struct CellGeom_t
   {
      // Cell geometry inner structure.

      Float_t fPhiMin;
      Float_t fPhiMax;
      Float_t fEtaMin;
      Float_t fEtaMax;

      Float_t fThetaMin; // cached
      Float_t fThetaMax; // cached

      CellGeom_t(): fPhiMin(0), fPhiMax(0), fEtaMin(0), fEtaMax(0), fThetaMin(0), fThetaMax(0) {}
      CellGeom_t(Float_t etaMin, Float_t etaMax, Float_t phiMin, Float_t phiMax) {Configure(etaMin, etaMax, phiMin, phiMax);}
      virtual ~CellGeom_t() {}

      void Configure(Float_t etaMin, Float_t etaMax, Float_t phiMin, Float_t phiMax);

      Float_t EtaMin()   const { return fEtaMin; }
      Float_t EtaMax()   const { return fEtaMax; }
      Float_t Eta()      const { return (fEtaMin+fEtaMax)*0.5f; }
      Float_t EtaDelta() const { return fEtaMax-fEtaMin; }

      Float_t PhiMin()   const { return fPhiMin; }
      Float_t PhiMax()   const { return fPhiMax; }
      Float_t Phi()      const { return (fPhiMin+fPhiMax)*0.5f; }
      Float_t PhiDelta() const { return fPhiMax-fPhiMin; }

      Float_t ThetaMin() const { return fThetaMin; }
      Float_t ThetaMax() const { return fThetaMax; }
      Float_t Theta() const { return (fThetaMax+fThetaMin)*0.5f; }
      Float_t ThetaDelta() const { return fThetaMax-fThetaMin; }

      Bool_t  IsUpperRho() const
      {
         const Float_t phi = Phi();
         return ((phi > 0 && phi <= TMath::Pi()) || phi < - TMath::Pi());
      }

      virtual void  Dump() const;
   };

   struct CellData_t : public CellGeom_t
   {
      // Cell data inner structure.

      Float_t fValue;

      CellData_t() : CellGeom_t(), fValue(0) {}
      virtual ~CellData_t() {}

      Float_t Value(Bool_t) const;
      virtual void Dump() const;
   };


   struct RebinData_t
   {
      Int_t fNSlices;

      std::vector<Float_t> fSliceData;
      std::vector<Int_t>   fBinData;

      Float_t* GetSliceVals(Int_t bin);

      void Clear()
      {
         fSliceData.clear();
         fBinData.clear();
      }
   };

   /**************************************************************************/

   typedef std::vector<CellId_t>               vCellId_t;
   typedef std::vector<CellId_t>::iterator     vCellId_i;

   typedef std::vector<CellGeom_t>             vCellGeom_t;
   typedef std::vector<CellGeom_t>::iterator   vCellGeom_i;
   typedef std::vector<CellGeom_t>::const_iterator   vCellGeom_ci;

private:
   TEveCaloData(const TEveCaloData&);            // Not implemented
   TEveCaloData& operator=(const TEveCaloData&); // Not implemented

protected:
   vSliceInfo_t fSliceInfos;

   TAxis*       fEtaAxis;
   TAxis*       fPhiAxis;

   Bool_t       fWrapTwoPi;

   Float_t      fMaxValEt; // cached
   Float_t      fMaxValE;  // cached

   Float_t      fEps;

   vCellId_t    fCellsSelected;
   vCellId_t    fCellsHighlighted;

public:
   TEveCaloData(const char* n="TEveCalData", const char* t="");
   virtual ~TEveCaloData() {}

   virtual void UnSelected();
   virtual void UnHighlighted();

   virtual TString GetHighlightTooltip();

   virtual void    FillImpliedSelectedSet(Set_t& impSelSet);

   virtual void    GetCellList(Float_t etaMin, Float_t etaMax,
                               Float_t phi,    Float_t phiRng,
                               vCellId_t &out) const = 0;

   vCellId_t&      GetCellsSelected()    { return fCellsSelected; }
   vCellId_t&      GetCellsHighlighted() { return fCellsHighlighted; }
   void            PrintCellsSelected();
   void            ProcessSelection(vCellId_t& sel_cells, TGLSelectRecord& rec);

   virtual void    Rebin(TAxis *ax, TAxis *ay, vCellId_t &in, Bool_t et, RebinData_t &out) const = 0;


   virtual void    GetCellData(const CellId_t &id, CellData_t& data) const = 0;

   virtual void    InvalidateUsersCellIdCache();
   virtual void    DataChanged();
   virtual void    CellSelectionChanged();

   Int_t           GetNSlices()    const { return fSliceInfos.size(); }
   SliceInfo_t&    RefSliceInfo(Int_t s) { return fSliceInfos[s]; }
   void            SetSliceThreshold(Int_t slice, Float_t threshold);
   Float_t         GetSliceThreshold(Int_t slice) const;
   void            SetSliceColor(Int_t slice, Color_t col);
   Color_t         GetSliceColor(Int_t slice) const;
   void            SetSliceTransparency(Int_t slice, Char_t t);
   Char_t          GetSliceTransparency(Int_t slice) const;

   virtual void    GetEtaLimits(Double_t &min, Double_t &max) const = 0;

   virtual void    GetPhiLimits(Double_t &min, Double_t &max) const = 0;

   virtual Float_t GetMaxVal(Bool_t et) const { return et ? fMaxValEt : fMaxValE; }
   Bool_t  Empty() const { return fMaxValEt < 1e-5; }

   virtual TAxis*  GetEtaBins()    const { return fEtaAxis; }
   virtual void    SetEtaBins(TAxis* ax) { fEtaAxis=ax; }

   virtual TAxis*  GetPhiBins()    const { return fPhiAxis; }
   virtual void    SetPhiBins(TAxis* ax) { fPhiAxis=ax; }

   virtual Float_t GetEps()      const { return fEps; }
   virtual void    SetEps(Float_t eps) { fEps=eps; }

   Bool_t   GetWrapTwoPi() const { return fWrapTwoPi; }
   void     SetWrapTwoPi(Bool_t w) { fWrapTwoPi=w; }

   static  Float_t EtaToTheta(Float_t eta);


   ClassDef(TEveCaloData, 0); // Manages calorimeter event data.
};

/**************************************************************************/
/**************************************************************************/

class TEveCaloDataVec: public TEveCaloData
{

private:
   TEveCaloDataVec(const TEveCaloDataVec&);            // Not implemented
   TEveCaloDataVec& operator=(const TEveCaloDataVec&); // Not implemented

protected:
   typedef std::vector<Float_t>               vFloat_t;
   typedef std::vector<Float_t>::iterator     vFloat_i;

   typedef std::vector<vFloat_t>              vvFloat_t;
   typedef std::vector<vFloat_t>::iterator    vvFloat_i;

   vvFloat_t   fSliceVec;
   vCellGeom_t fGeomVec;

   Int_t       fTower; // current tower

   Float_t     fEtaMin;
   Float_t     fEtaMax;

   Float_t     fPhiMin;
   Float_t     fPhiMax;

public:
   TEveCaloDataVec(Int_t nslices);
   virtual ~TEveCaloDataVec();

   Int_t AddSlice();
   Int_t AddTower(Float_t etaMin, Float_t etaMax, Float_t phiMin, Float_t phiMax);
   void  FillSlice(Int_t slice, Float_t value);
   void  FillSlice(Int_t slice, Int_t tower, Float_t value);

   Int_t GetNCells() { return fGeomVec.size(); }
   std::vector<Float_t>&  GetSliceVals(Int_t slice) { return fSliceVec[slice]; }
   std::vector<TEveCaloData::CellGeom_t>& GetCellGeom() { return fGeomVec; }

   virtual void GetCellList(Float_t etaMin, Float_t etaMax,
                            Float_t phi,    Float_t phiRng,
                            vCellId_t &out) const;

   virtual void Rebin(TAxis *ax, TAxis *ay, vCellId_t &in, Bool_t et, RebinData_t &out) const;

   virtual void GetCellData(const TEveCaloData::CellId_t &id, TEveCaloData::CellData_t& data) const;
   virtual void GetEtaLimits(Double_t &min, Double_t &max) const { min=fEtaMin, max=fEtaMax;}
   virtual void GetPhiLimits(Double_t &min, Double_t &max) const { min=fPhiMin; max=fPhiMax;}


   virtual void  DataChanged();
   void          SetAxisFromBins(Double_t epsX=0.001, Double_t epsY=0.001);

   ClassDef(TEveCaloDataVec, 0); // Manages calorimeter event data.
};

/**************************************************************************/
/**************************************************************************/

class TEveCaloDataHist: public TEveCaloData
{
private:
   TEveCaloDataHist(const TEveCaloDataHist&);            // Not implemented
   TEveCaloDataHist& operator=(const TEveCaloDataHist&); // Not implemented

protected:
   THStack*    fHStack;

public:
   TEveCaloDataHist();
   virtual ~TEveCaloDataHist();

   virtual void GetCellList( Float_t etaMin, Float_t etaMax,
                             Float_t phi, Float_t phiRng, vCellId_t &out) const;

   virtual void Rebin(TAxis *ax, TAxis *ay, vCellId_t &in, Bool_t et, RebinData_t &out) const;

   virtual void GetCellData(const TEveCaloData::CellId_t &id, TEveCaloData::CellData_t& data) const;

   virtual void GetEtaLimits(Double_t &min, Double_t &max) const;
   virtual void GetPhiLimits(Double_t &min, Double_t &max) const;


   virtual void DataChanged();

   THStack* GetStack() { return fHStack; }

   TH2F*    GetHist(Int_t slice) const;

   Int_t   AddHistogram(TH2F* hist);

   ClassDef(TEveCaloDataHist, 0); // Manages calorimeter TH2F event data.
};

#endif

