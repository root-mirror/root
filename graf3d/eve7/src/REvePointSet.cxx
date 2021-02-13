// @(#)root/eve7:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/REvePointSet.hxx>

#include <ROOT/REveManager.hxx>
#include <ROOT/REveProjectionManager.hxx>
#include <ROOT/REveTrans.hxx>
#include <ROOT/REveRenderData.hxx>

#include "TArrayI.h"
#include "TClass.h"


using namespace ROOT::Experimental;

/** \class REvePointSet
\ingroup REve
REvePointSet is a render-element holding a collection of 3D points with
optional per-point TRef and an arbitrary number of integer ids (to
be used for signal, volume-id, track-id, etc).

3D point representation is implemented in base-class TPolyMarker3D.
Per-point TRef is implemented in base-class TPointSet3D.

By using the REvePointSelector the points and integer ids can be
filled directly from a TTree holding the source data.
Setting of per-point TRef's is not supported.

REvePointSet is a REveProjectable: it can be projected by using the
REveProjectionManager class.
*/

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

REvePointSet::REvePointSet(const std::string& name, const std::string& title, Int_t n_points) :
   REveElement(name, title),
   TAttMarker(),
   TAttBBox()
{
   fMarkerStyle = 20;

   SetMainColorPtr(&fMarkerColor);

   Reset(n_points);

   // Override from REveElement.
   fPickable = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.

REvePointSet::REvePointSet(const REvePointSet& e) :
   REveElement(e),
   REveProjectable(e),
   TAttMarker(e),
   TAttBBox(e)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

REvePointSet::~REvePointSet()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Clone points and all point-related information from point-set 'e'.

void REvePointSet::ClonePoints(const REvePointSet& e)
{
   fPoints   = e.fPoints;
   fCapacity = e.fCapacity;
   fSize     = e.fSize;
}

////////////////////////////////////////////////////////////////////////////////
/// Drop all data and set-up the data structures to recive new data.
/// n_points   specifies the initial size of the array.

void REvePointSet::Reset(Int_t n_points)
{
   fPoints.resize(n_points);
   fCapacity = n_points;
   fSize     = 0;

   ResetBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Resizes internal array to allow additional n_points to be stored.
/// Returns the old size which is also the location where one can
/// start storing new data.
/// The caller is *obliged* to fill the new point slots.

Int_t REvePointSet::GrowFor(Int_t n_points)
{
   assert(n_points >= 0);

   Int_t old_size = fCapacity;
   Int_t new_size = old_size + n_points;

   fPoints.resize(new_size);
   fCapacity = new_size;

   return old_size;
}

int REvePointSet::SetNextPoint(float x, float y, float z)
{
   return SetPoint(fSize, x, y, z);
}

int REvePointSet::SetPoint(int n, float x, float y, float z)
{
   if (n >= fCapacity)
   {
      fCapacity = std::max(n + 1, 2*fCapacity);
      fPoints.resize(fCapacity);
   }
   fPoints[n].Set(x, y, z);
   if (n >= fSize)
   {
      fSize = n + 1;
   }
   return fSize;
}

////////////////////////////////////////////////////////////////////////////////
/// Set marker style, propagate to projecteds.

void REvePointSet::SetMarkerStyle(Style_t mstyle)
{
   for (auto &pi: fProjectedList)
   {
      REvePointSet* pt = dynamic_cast<REvePointSet *>(pi);
      if (pt)
      {
         pt->SetMarkerStyle(mstyle);
         pt->StampObjProps();
      }
   }
   TAttMarker::SetMarkerStyle(mstyle);
}

////////////////////////////////////////////////////////////////////////////////
/// Set marker size, propagate to projecteds.

void REvePointSet::SetMarkerSize(Size_t msize)
{
   for (auto &pi: fProjectedList)
   {
      REvePointSet* pt = dynamic_cast<REvePointSet *>(pi);
      if (pt)
      {
         pt->SetMarkerSize(msize);
         pt->StampObjProps();
      }
   }
   TAttMarker::SetMarkerSize(msize);
   StampObjProps();
}

////////////////////////////////////////////////////////////////////////////////
/// Copy visualization parameters from element el.

void REvePointSet::CopyVizParams(const REveElement* el)
{
   const REvePointSet* m = dynamic_cast<const REvePointSet*>(el);
   if (m)
   {
      TAttMarker::operator=(*m);
   }

   REveElement::CopyVizParams(el);
}

////////////////////////////////////////////////////////////////////////////////
/// Write visualization parameters.

void REvePointSet::WriteVizParams(std::ostream& out, const TString& var)
{
   REveElement::WriteVizParams(out, var);

   TAttMarker::SaveMarkerAttributes(out, var);
}

////////////////////////////////////////////////////////////////////////////////
/// Virtual from REveProjectable, returns REvePointSetProjected class.

TClass* REvePointSet::ProjectedClass(const REveProjection*) const
{
   return TClass::GetClass<REvePointSetProjected>();
}


Int_t REvePointSet::WriteCoreJson(nlohmann::json& j, Int_t rnr_offset)
{
   Int_t ret = REveElement::WriteCoreJson(j, rnr_offset);

   j["fMarkerSize"]  = GetMarkerSize();
   j["fMarkerColor"] = GetMarkerColor();

   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// Crates 3D point array for rendering.

void REvePointSet::BuildRenderData()
{
   if (fSize > 0)
   {
      fRenderData = std::make_unique<REveRenderData>("makeHit", 3*fSize);
      fRenderData->PushV(&fPoints[0].fX, 3*fSize);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Compute bounding box.

void REvePointSet::ComputeBBox()
{
   if (fSize > 0) {
      BBoxInit();
      for (auto &p : fPoints)
      {
         BBoxCheckPoint(p.fX, p.fY, p.fZ);
      }
   } else {
      BBoxZero();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Virtual method of base class TPointSet3D. The function call is
/// invoked with secondary selection in TPointSet3DGL.

void REvePointSet::PointSelected(Int_t /* id */)
{
   // Emit("PointSelected(Int_t)", id);
}


//==============================================================================
// REvePointSetArray
//==============================================================================

/** \class REvePointSetArray
\ingroup REve
An array of point-sets with each point-set playing a role of a bin
in a histogram. When a new point is added to a REvePointSetArray,
an additional separating quantity needs to be specified: it
determines into which REvePointSet (bin) the point will actually be
stored. Underflow and overflow bins are automatically created but
they are not drawn by default.

By using the REvePointSelector the points and the separating
quantities can be filled directly from a TTree holding the source
data.
Setting of per-point TRef's is not supported.

After the filling, the range of separating variable can be
controlled with a slider to choose a sub-set of PointSets that are
actually shown.
*/

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

REvePointSetArray::REvePointSetArray(const std::string& name,
                                     const std::string& title) :
   REveElement(name, title),

   fBins(nullptr), fDefPointSetCapacity(128), fNBins(0), fLastBin(-1),
   fMin(0), fCurMin(0), fMax(0), fCurMax(0),
   fBinWidth(0),
   fQuantName()
{
   SetMainColorPtr(&fMarkerColor);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor: deletes the fBins array. Actual removal of
/// elements done by REveElement.

REvePointSetArray::~REvePointSetArray()
{
   // printf("REvePointSetArray::~REvePointSetArray()\n");
   delete [] fBins; fBins = nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Virtual from REveElement, provide bin management.

void REvePointSetArray::RemoveElementLocal(REveElement* el)
{
   for (Int_t i=0; i<fNBins; ++i) {
      if (fBins[i] == el) {
         fBins[i] = nullptr;
         break;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Virtual from REveElement, provide bin management.

void REvePointSetArray::RemoveElementsLocal()
{
   delete [] fBins; fBins = nullptr; fLastBin = -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Set marker color, propagate to children.

void REvePointSetArray::SetMarkerColor(Color_t tcolor)
{
   for (auto & el : fChildren)
   {
      TAttMarker* m = dynamic_cast<TAttMarker*>(el);
      if (m && m->GetMarkerColor() == fMarkerColor)
         m->SetMarkerColor(tcolor);
   }
   TAttMarker::SetMarkerColor(tcolor);
}

////////////////////////////////////////////////////////////////////////////////
/// Set marker style, propagate to children.

void REvePointSetArray::SetMarkerStyle(Style_t mstyle)
{
   for (auto & el : fChildren)
   {
      TAttMarker* m = dynamic_cast<TAttMarker*>(el);
      if (m && m->GetMarkerStyle() == fMarkerStyle)
         m->SetMarkerStyle(mstyle);
   }
   TAttMarker::SetMarkerStyle(mstyle);
}

////////////////////////////////////////////////////////////////////////////////
/// Set marker size, propagate to children.

void REvePointSetArray::SetMarkerSize(Size_t msize)
{
   for (auto & el : fChildren)
   {
      TAttMarker* m = dynamic_cast<TAttMarker*>(el);
      if (m && m->GetMarkerSize() == fMarkerSize)
         m->SetMarkerSize(msize);
   }
   TAttMarker::SetMarkerSize(msize);
}

////////////////////////////////////////////////////////////////////////////////
/// Get the total number of filled points.
/// 'under' and 'over' flags specify if under/overflow channels
/// should be added to the sum.

Int_t REvePointSetArray::Size(Bool_t under, Bool_t over) const
{
   Int_t size = 0;
   const Int_t min = under ? 0 : 1;
   const Int_t max = over  ? fNBins : fNBins - 1;
   for (Int_t i = min; i < max; ++i)
   {
      if (fBins[i])
         size += fBins[i]->GetSize();
   }
   return size;
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize internal point-sets with given binning parameters.
/// The actual number of bins is nbins+2, bin 0 corresponding to
/// underflow and bin nbin+1 to owerflow pointset.

void REvePointSetArray::InitBins(const std::string& quant_name,
                                 Int_t nbins, Double_t min, Double_t max)
{
   static const REveException eh("REvePointSetArray::InitBins ");

   if (nbins < 1) throw eh + "nbins < 1.";
   if (min > max) throw eh + "min > max.";

   RemoveElements();

   fQuantName = quant_name;
   fNBins     = nbins + 2; // under/overflow
   fLastBin   = -1;
   fMin = fCurMin = min;
   fMax = fCurMax = max;
   fBinWidth  = (fMax - fMin)/(fNBins - 2);

   fBins = new REvePointSet* [fNBins];

   for (Int_t i = 0; i < fNBins; ++i)
   {
      fBins[i] = new REvePointSet
         (Form("Slice %d [%4.3lf, %4.3lf]", i, fMin + (i-1)*fBinWidth, fMin + i*fBinWidth),
          "",
          fDefPointSetCapacity);
      fBins[i]->SetMarkerColor(fMarkerColor);
      fBins[i]->SetMarkerStyle(fMarkerStyle);
      fBins[i]->SetMarkerSize(fMarkerSize);
      AddElement(fBins[i]);
   }

   fBins[0]->SetName("Underflow");
   fBins[0]->SetRnrSelf(kFALSE);

   fBins[fNBins-1]->SetName("Overflow");
   fBins[fNBins-1]->SetRnrSelf(kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// Add a new point. Appropriate point-set will be chosen based on
/// the value of the separating quantity 'quant'.
/// If the selected bin does not have an associated REvePointSet
/// the point is discarded and false is returned.

Bool_t REvePointSetArray::Fill(Double_t x, Double_t y, Double_t z, Double_t quant)
{
   fLastBin = TMath::FloorNint((quant - fMin)/fBinWidth) + 1;

   if (fLastBin < 0)
   {
      fLastBin = 0;
   }
   else if (fLastBin > fNBins - 1)
   {
      fLastBin = fNBins - 1;
   }

   if (fBins[fLastBin] != 0)
   {
      fBins[fLastBin]->SetNextPoint(x, y, z);
      return kTRUE;
   }
   else
   {
      return kFALSE;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Call this after all the points have been filled.
/// At this point we can calculate bounding-boxes of individual
/// point-sets.

void REvePointSetArray::CloseBins()
{
   for (Int_t i=0; i<fNBins; ++i)
   {
      if (fBins[i])
      {
         fBins[i]->SetTitle(Form("N=%d", fBins[i]->GetSize()));
         fBins[i]->ComputeBBox();
      }
   }
   fLastBin = -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Set active range of the separating quantity.
/// Appropriate point-sets are tagged for rendering.
/// Over/underflow point-sets are left as they were.

void REvePointSetArray::SetRange(Double_t min, Double_t max)
{
   using namespace TMath;

   fCurMin = min; fCurMax = max;
   Int_t  low_b = Max(0,        FloorNint((min-fMin)/fBinWidth)) + 1;
   Int_t high_b = Min(fNBins-2, CeilNint ((max-fMin)/fBinWidth));

   for (Int_t i = 1; i < fNBins - 1; ++i)
   {
      if (fBins[i])
         fBins[i]->SetRnrSelf(i>=low_b && i<=high_b);
   }
}

/** \class REvePointSetProjected
\ingroup REve
Projected copy of a REvePointSet.
*/

////////////////////////////////////////////////////////////////////////////////
/// Default contructor.

REvePointSetProjected::REvePointSetProjected() :
   REvePointSet  (),
   REveProjected ()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Set projection manager and projection model.
/// Virtual from REveProjected.

void REvePointSetProjected::SetProjection(REveProjectionManager* proj,
                                          REveProjectable* model)
{
   REveProjected::SetProjection(proj, model);
   CopyVizParams(dynamic_cast<REveElement*>(model));
}

////////////////////////////////////////////////////////////////////////////////
/// Set depth (z-coordinate) of the projected points.

void REvePointSetProjected::SetDepthLocal(Float_t d)
{
   SetDepthCommon(d, this, fBBox);

   // XXXX rewrite

   Int_t    n = fSize;
   Float_t *p = & fPoints[0].fZ;
   for (Int_t i = 0; i < n; ++i, p+=3)
      *p = fDepth;
}

////////////////////////////////////////////////////////////////////////////////
/// Re-apply the projection.
/// Virtual from REveProjected.

void REvePointSetProjected::UpdateProjection()
{
   REveProjection &proj = * fManager->GetProjection();
   REvePointSet   &ps   = * dynamic_cast<REvePointSet*>(fProjectable);
   REveTrans      *tr   =   ps.PtrMainTrans(kFALSE);

   // XXXX rewrite

   Int_t n = ps.GetSize();
   Reset(n);
   fSize = n;
   const Float_t *o = & ps.RefPoint(0).fX;
         Float_t *p = & fPoints[0].fX;
   for (Int_t i = 0; i < n; ++i, o+=3, p+=3)
   {
      proj.ProjectPointfv(tr, o, p, fDepth);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Virtual method of base class TPointSet3D.
/// Forward to projectable.

void REvePointSetProjected::PointSelected(Int_t id)
{
   REvePointSet *ps = dynamic_cast<REvePointSet*>(fProjectable);
   ps->PointSelected(id);
}
