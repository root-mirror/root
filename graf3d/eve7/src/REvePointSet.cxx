// @(#)root/eve7:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
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

#include "TColor.h"
#include "TArrayI.h"

#include "json.hpp"


using namespace ROOT::Experimental;
namespace REX = ROOT::Experimental;

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

REvePointSet::REvePointSet(Int_t n_points, ETreeVarType_e tv_type) :
   REveElement(),
   TPointSet3D(n_points),
   REvePointSelectorConsumer(tv_type),
   REveProjectable(),

   fTitle          (),
   fIntIds         (nullptr),
   fIntIdsPerPoint (0)
{
   fMarkerStyle = 20;
   SetMainColorPtr(&fMarkerColor);

   // Override from REveElement.
   fPickable = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

REvePointSet::REvePointSet(const char* name, Int_t n_points, ETreeVarType_e tv_type) :
   REveElement(),
   TPointSet3D(n_points),
   REvePointSelectorConsumer(tv_type),
   REveProjectable(),

   fTitle          (),
   fIntIds         (nullptr),
   fIntIdsPerPoint (0)
{
   fMarkerStyle = 20;
   SetName(name);
   SetMainColorPtr(&fMarkerColor);

   // Override from REveElement.
   fPickable = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.

REvePointSet::REvePointSet(const REvePointSet& e) :
   REveElement(e),
   TPointSet3D(e),
   REvePointSelectorConsumer(e),
   REveProjectable(),

   fTitle          (e.fTitle),
   fIntIds         (e.fIntIds ? new TArrayI(*e.fIntIds) : nullptr),
   fIntIdsPerPoint (e.fIntIdsPerPoint)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

REvePointSet::~REvePointSet()
{
   delete fIntIds;
}

////////////////////////////////////////////////////////////////////////////////
/// Clone points and all point-related information from point-set 'e'.

void REvePointSet::ClonePoints(const REvePointSet& e)
{
   // TPolyMarker3D
   delete [] fP;
   fN = e.fN;
   if (fN > 0)
   {
      const Int_t nn = 3 * e.fN;
      fP = new Float_t [nn];
      for (Int_t i = 0; i < nn; i++) fP[i] = e.fP[i];
   } else {
      fP = 0;
   }
   fLastPoint = e.fLastPoint;

   // TPointSet3D
   CopyIds(e);

   // REvePointSet
   delete fIntIds;
   fIntIds         = e.fIntIds ? new TArrayI(*e.fIntIds) : nullptr;
   fIntIdsPerPoint = e.fIntIdsPerPoint;
}

////////////////////////////////////////////////////////////////////////////////
/// Drop all data and set-up the data structures to recive new data.
/// n_points   specifies the initial size of the arrays.
/// n_int_ids  specifies the number of integer ids per point.

void REvePointSet::Reset(Int_t n_points, Int_t n_int_ids)
{
   delete [] fP; fP = nullptr;
   fN = n_points;
   if (fN) {
      fP = new Float_t [3*fN];
      memset(fP, 0, 3*fN*sizeof(Float_t));
   }
   fLastPoint = -1;
   ClearIds();
   delete fIntIds; fIntIds = nullptr;
   fIntIdsPerPoint = n_int_ids;
   if (fIntIdsPerPoint > 0) fIntIds = new TArrayI(fIntIdsPerPoint*fN);
   ResetBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Resizes internal array to allow additional n_points to be stored.
/// Returns the old size which is also the location where one can
/// start storing new data.
/// The caller is *obliged* to fill the new point slots.

Int_t REvePointSet::GrowFor(Int_t n_points)
{
   Int_t old_size = Size();
   Int_t new_size = old_size + n_points;
   SetPoint(new_size - 1, 0, 0, 0);
   if (fIntIds)
      fIntIds->Set(fIntIdsPerPoint * new_size);
   return old_size;
}

////////////////////////////////////////////////////////////////////////////////
/// Assert that size of IntId array is compatible with the size of
/// the point array.

inline void REvePointSet::AssertIntIdsSize()
{
   Int_t exp_size = GetN()*fIntIdsPerPoint;
   if (fIntIds->GetSize() < exp_size)
      fIntIds->Set(exp_size);
}

////////////////////////////////////////////////////////////////////////////////
/// Return a pointer to integer ids of point with index p.
/// Existence of integer id array is checked, 0 is returned if it
/// does not exist.
/// Validity of p is *not* checked.

Int_t* REvePointSet::GetPointIntIds(Int_t p) const
{
   if (fIntIds)
      return fIntIds->GetArray() + p*fIntIdsPerPoint;
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Return i-th integer id of point with index p.
/// Existence of integer id array is checked, kMinInt is returned if
/// it does not exist.
/// Validity of p and i is *not* checked.

Int_t REvePointSet::GetPointIntId(Int_t p, Int_t i) const
{
   if (fIntIds)
      return * (fIntIds->GetArray() + p*fIntIdsPerPoint + i);
   return kMinInt;
}

////////////////////////////////////////////////////////////////////////////////
/// Set integer ids for the last point that was registered (most
/// probably via TPolyMarker3D::SetNextPoint(x,y,z)).

void REvePointSet::SetPointIntIds(Int_t* ids)
{
   SetPointIntIds(fLastPoint, ids);
}

////////////////////////////////////////////////////////////////////////////////
/// Set integer ids for point with index n.

void REvePointSet::SetPointIntIds(Int_t n, Int_t* ids)
{
   if (!fIntIds) return;
   AssertIntIdsSize();
   Int_t* x = fIntIds->GetArray() + n*fIntIdsPerPoint;
   for (Int_t i=0; i<fIntIdsPerPoint; ++i)
      x[i] = ids[i];
}

////////////////////////////////////////////////////////////////////////////////
/// Set marker style, propagate to projecteds.

void REvePointSet::SetMarkerStyle(Style_t mstyle)
{
   static const REveException eh("REvePointSet::SetMarkerStyle ");

   std::list<REveProjected*>::iterator pi = fProjectedList.begin();
   while (pi != fProjectedList.end())
   {
      REvePointSet* pt = dynamic_cast<REvePointSet*>(*pi);
      if (pt)
      {
         pt->SetMarkerStyle(mstyle);
         pt->StampObjProps();
      }
      ++pi;
   }
   TAttMarker::SetMarkerStyle(mstyle);
}

////////////////////////////////////////////////////////////////////////////////
/// Set marker size, propagate to projecteds.

void REvePointSet::SetMarkerSize(Size_t msize)
{
   static const REveException eh("REvePointSet::SetMarkerSize ");

   std::list<REveProjected*>::iterator pi = fProjectedList.begin();
   while (pi != fProjectedList.end())
   {
      REvePointSet* pt = dynamic_cast<REvePointSet*>(*pi);
      if (pt)
      {
         pt->SetMarkerSize(msize);
         pt->StampObjProps();
      }
      ++pi;
   }
   TAttMarker::SetMarkerSize(msize);
   StampObjProps();
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize point-set for new filling.
/// subIdNum gives the number of integer ids that can be assigned to
/// each point.

void REvePointSet::InitFill(Int_t subIdNum)
{
   if (subIdNum > 0) {
      fIntIdsPerPoint = subIdNum;
      if (!fIntIds)
         fIntIds = new TArrayI(fIntIdsPerPoint*GetN());
      else
         fIntIds->Set(fIntIdsPerPoint*GetN());
   } else {
      delete fIntIds; fIntIds = nullptr;
      fIntIdsPerPoint = 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Called from REvePointSelector when internal arrays of the tree-selector
/// are filled up and need to be processed.
/// Virtual from REvePointSelectorConsumer.

void REvePointSet::TakeAction(REvePointSelector* sel)
{
   static const REveException eh("REvePointSet::TakeAction ");

   if(sel == nullptr)
      throw(eh + "selector is <null>.");

   Int_t    n = sel->GetNfill();
   Int_t  beg = GrowFor(n);

   // printf("REvePointSet::TakeAction beg=%d n=%d size=%d nsubid=%d dim=%d\n",
   //        beg, n, Size(), sel->GetSubIdNum(), sel->GetDimension());

   Double_t *vx = sel->GetV1(), *vy = sel->GetV2(), *vz = sel->GetV3();
   Float_t  *p  = fP + 3*beg;

   switch(fSourceCS) {
      case kTVT_XYZ:
         while(n-- > 0) {
            p[0] = *vx; p[1] = *vy; p[2] = *vz;
            p += 3;
            ++vx; ++vy; ++vz;
         }
         break;
      case kTVT_RPhiZ:
         while(n-- > 0) {
            p[0] = *vx * TMath::Cos(*vy); p[1] = *vx * TMath::Sin(*vy); p[2] = *vz;
            p += 3;
            ++vx; ++vy; ++vz;
         }
         break;
      default:
         throw(eh + "unknown tree variable type.");
   }

   if (fIntIds) {
      Double_t** subarr = new Double_t* [fIntIdsPerPoint];
      for (Int_t i=0; i<fIntIdsPerPoint; ++i) {
         subarr[i] = sel->GetVal(sel->GetDimension() - fIntIdsPerPoint + i);
         if (subarr[i] == 0) {
            delete[] subarr;
            throw(eh + "sub-id array not available.");
         }
      }
      Int_t* ids = fIntIds->GetArray() + fIntIdsPerPoint*beg;
      n = sel->GetNfill();
      while (n-- > 0) {
         for (Int_t i=0; i<fIntIdsPerPoint; ++i) {
            ids[i] = TMath::Nint(*subarr[i]);
            ++subarr[i];
         }
         ids += fIntIdsPerPoint;
      }
      delete [] subarr;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Copy visualization parameters from element el.

void REvePointSet::CopyVizParams(const REveElement* el)
{
   const REvePointSet* m = dynamic_cast<const REvePointSet*>(el);
   if (m)
   {
      TAttMarker::operator=(*m);
      fOption = m->fOption;
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
   return REvePointSetProjected::Class();
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
   fRenderData = std::make_unique<REveRenderData>("makeHit", 3*fN);

   fRenderData->PushV(fP, 3*fN);
}

////////////////////////////////////////////////////////////////////////////////
/// Virtual method of base class TPointSet3D. The function call is
/// invoked with secondary selection in TPointSet3DGL.

void REvePointSet::PointSelected(Int_t id)
{
   // Emit("PointSelected(Int_t)", id);
   TPointSet3D::PointSelected(id);
}

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

REvePointSetArray::REvePointSetArray(const char* name,
                                     const char* title) :
   REveElement(),
   TNamed(name, title),

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
   static const REveException eh("REvePointSetArray::SetMarkerColor ");

   for (List_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
      TAttMarker* m = dynamic_cast<TAttMarker*>((*i)->GetObject(eh));
      if (m && m->GetMarkerColor() == fMarkerColor)
         m->SetMarkerColor(tcolor);
   }
   TAttMarker::SetMarkerColor(tcolor);
}

////////////////////////////////////////////////////////////////////////////////
/// Set marker style, propagate to children.

void REvePointSetArray::SetMarkerStyle(Style_t mstyle)
{
   static const REveException eh("REvePointSetArray::SetMarkerStyle ");

   for (List_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
      TAttMarker* m = dynamic_cast<TAttMarker*>((*i)->GetObject(eh));
      if (m && m->GetMarkerStyle() == fMarkerStyle)
         m->SetMarkerStyle(mstyle);
   }
   TAttMarker::SetMarkerStyle(mstyle);
}

////////////////////////////////////////////////////////////////////////////////
/// Set marker size, propagate to children.

void REvePointSetArray::SetMarkerSize(Size_t msize)
{
   static const REveException eh("REvePointSetArray::SetMarkerSize ");

   for (List_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
      TAttMarker* m = dynamic_cast<TAttMarker*>((*i)->GetObject(eh));
      if (m && m->GetMarkerSize() == fMarkerSize)
         m->SetMarkerSize(msize);
   }
   TAttMarker::SetMarkerSize(msize);
}

////////////////////////////////////////////////////////////////////////////////
/// Called from REvePointSelector when internal arrays of the tree-selector
/// are filled up and need to be processed.
/// Virtual from REvePointSelectorConsumer.

void REvePointSetArray::TakeAction(REvePointSelector* sel)
{
   static const REveException eh("REvePointSetArray::TakeAction ");

   if (sel == nullptr)
      throw eh + "selector is <null>.";

   Int_t n = sel->GetNfill();

   // printf("REvePointSetArray::TakeAction n=%d\n", n);

   Double_t *vx = sel->GetV1(), *vy = sel->GetV2(), *vz = sel->GetV3();
   Double_t *qq = sel->GetV4();

   if (qq == nullptr)
      throw eh + "requires 4-d varexp.";

   switch (fSourceCS)
   {
      case kTVT_XYZ:
      {
         while (n-- > 0)
         {
            Fill(*vx, *vy, *vz, *qq);
            ++vx; ++vy; ++vz; ++qq;
         }
         break;
      }
      case kTVT_RPhiZ:
      {
         while (n-- > 0)
         {
            Fill(*vx * TMath::Cos(*vy), *vx * TMath::Sin(*vy), *vz, *qq);
            ++vx; ++vy; ++vz; ++qq;
         }
         break;
      }
      default:
      {
         throw eh + "unknown tree variable type.";
      }
   }
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
         size += fBins[i]->Size();
   }
   return size;
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize internal point-sets with given binning parameters.
/// The actual number of bins is nbins+2, bin 0 corresponding to
/// underflow and bin nbin+1 to owerflow pointset.

void REvePointSetArray::InitBins(const char* quant_name,
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
/// Set external object id of the last added point.

void REvePointSetArray::SetPointId(TObject* id)
{
   if (fLastBin >= 0)
      fBins[fLastBin]->SetPointId(id);
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
         fBins[i]->SetTitle(Form("N=%d", fBins[i]->Size()));
         fBins[i]->ComputeBBox();
      }
   }
   fLastBin = -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Propagate id-object ownership to children.

void REvePointSetArray::SetOwnIds(Bool_t o)
{
   for (Int_t i=0; i<fNBins; ++i)
   {
      if (fBins[i])
         fBins[i]->SetOwnIds(o);
   }
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

   Int_t    n = Size();
   Float_t *p = GetP() + 2;
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

   Int_t n = ps.Size();
   Reset(n);
   fLastPoint = n - 1;
   Float_t *o = ps.GetP(), *p = GetP();
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
