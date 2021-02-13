// @(#)root/eve7:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_REvePointSet
#define ROOT7_REvePointSet

#include <ROOT/REveElement.hxx>
#include <ROOT/REveProjectionBases.hxx>
#include <ROOT/REveVector.hxx>

#include <TAttMarker.h>
#include <TAttBBox.h>

#include <cassert>

namespace ROOT {
namespace Experimental {

////////////////////////////////////////////////////////////////////////////////
// REvePointSet
// Set of 3D points with same marker attributes;
// optionally each point can be assigned an
// external TRef or a number of integer indices.
////////////////////////////////////////////////////////////////////////////////

class REvePointSet : public REveElement,
                     public REveProjectable,
                     public TAttMarker,
                     public TAttBBox
{
   friend class REvePointSetArray;

private:
   REvePointSet &operator=(const REvePointSet &) = delete;

protected:
   std::vector<REveVector> fPoints;
   int                     fCapacity{0};
   int                     fSize{0};

public:
   REvePointSet(const std::string& name="", const std::string& title="", Int_t n_points = 0);
   REvePointSet(const REvePointSet &e);
   virtual ~REvePointSet();

   REvePointSet *CloneElement() const override { return new REvePointSet(*this); }

   virtual void ClonePoints(const REvePointSet &e);

   void  Reset(Int_t n_points = 0);
   Int_t GrowFor(Int_t n_points);

   int   SetNextPoint(float x, float y, float z);
   int   SetPoint(int n, float x, float y, float z);

   int   GetCapacity() const { return fCapacity; }
   int   GetSize()     const { return fSize;     }

         REveVector& RefPoint(int n)       { assert (n < fSize); return fPoints[n]; }
   const REveVector& RefPoint(int n) const { assert (n < fSize); return fPoints[n]; }

   void SetMarkerColor(Color_t col) override { SetMainColor(col); }
   void SetMarkerStyle(Style_t mstyle = 1) override;
   void SetMarkerSize(Size_t msize = 1) override;

   void CopyVizParams(const REveElement *el) override;
   void WriteVizParams(std::ostream &out, const TString &var) override;

   TClass* ProjectedClass(const REveProjection *p) const override;

   Int_t WriteCoreJson(nlohmann::json &j, Int_t rnr_offset) override;
   void  BuildRenderData()override;

   void ComputeBBox() override;

   void PointSelected(Int_t id); // *SIGNAL*
};

/******************************************************************************/
// REvePointSetArray
// Array of REvePointSet's filled via a common point-source; range of displayed REvePointSet's can be
// controlled, based on a separating quantity provided on fill-time by a user.
/******************************************************************************/

class REvePointSetArray : public REveElement,
                          public REveProjectable,
                          public TAttMarker
{
   REvePointSetArray(const REvePointSetArray &) = delete;
   REvePointSetArray &operator=(const REvePointSetArray &) = delete;

protected:
   REvePointSet **fBins{nullptr};       //  Pointers to subjugated REvePointSet's.
   Int_t          fDefPointSetCapacity; //  Default capacity of subjugated REvePointSet's.
   Int_t          fNBins;               //  Number of subjugated REvePointSet's.
   Int_t          fLastBin;             //! Index of the last filled REvePointSet.
   Double_t       fMin, fCurMin;        //  Overall and current minimum value of the separating quantity.
   Double_t       fMax, fCurMax;        //  Overall and current maximum value of the separating quantity.
   Double_t       fBinWidth;            //  Separating quantity bin-width.
   std::string    fQuantName;           //  Name of the separating quantity.

public:
   REvePointSetArray(const std::string &name = "REvePointSetArray", const std::string &title = "");
   virtual ~REvePointSetArray();

   void RemoveElementLocal(REveElement *el) override;
   void RemoveElementsLocal() override;

   void SetMarkerColor(Color_t tcolor = 1) override;
   void SetMarkerStyle(Style_t mstyle = 1) override;
   void SetMarkerSize(Size_t msize = 1) override;

   Int_t Size(Bool_t under = kFALSE, Bool_t over = kFALSE) const;

   void InitBins(const std::string& quant_name, Int_t nbins, Double_t min, Double_t max);
   Bool_t Fill(Double_t x, Double_t y, Double_t z, Double_t quant);
   void CloseBins();

   Int_t GetDefPointSetCapacity()  const { return fDefPointSetCapacity; }
   void  SetDefPointSetCapacity(Int_t c) { fDefPointSetCapacity = c; }

   Int_t GetNBins() const { return fNBins; }
   REvePointSet *GetBin(Int_t bin) const { return fBins[bin]; }

   Double_t GetMin() const { return fMin; }
   Double_t GetCurMin() const { return fCurMin; }
   Double_t GetMax() const { return fMax; }
   Double_t GetCurMax() const { return fCurMax; }

   void SetRange(Double_t min, Double_t max);
};

/******************************************************************************/
// REvePointSetProjected
// Projected copy of a REvePointSet.
/******************************************************************************/

class REvePointSetProjected : public REvePointSet,
                              public REveProjected
{
private:
   REvePointSetProjected(const REvePointSetProjected &) = delete;
   REvePointSetProjected &operator=(const REvePointSetProjected &) = delete;

protected:
   void SetDepthLocal(Float_t d) override;

public:
   REvePointSetProjected();
   virtual ~REvePointSetProjected() {}

   void SetProjection(REveProjectionManager *proj, REveProjectable *model) override;
   void UpdateProjection() override;
   REveElement *GetProjectedAsElement() override { return this; }

   void PointSelected(Int_t id);
};

} // namespace Experimental
} // namespace ROOT

#endif
