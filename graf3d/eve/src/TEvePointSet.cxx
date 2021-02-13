// @(#)root/eve:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TEvePointSet.h"

#include "TEveManager.h"
#include "TEveProjectionManager.h"
#include "TEveTrans.h"

#include "TTree.h"
#include "TTreePlayer.h"
#include "TF3.h"

/** \class TEvePointSet
\ingroup TEve
TEvePointSet is a render-element holding a collection of 3D points with
optional per-point TRef and an arbitrary number of integer ids (to
be used for signal, volume-id, track-id, etc).

3D point representation is implemented in base-class TPolyMarker3D.
Per-point TRef is implemented in base-class TPointSet3D.

By using the TEvePointSelector the points and integer ids can be
filled directly from a TTree holding the source data.
Setting of per-point TRef's is not supported.

TEvePointSet is a TEveProjectable: it can be projected by using the
TEveProjectionManager class.
*/

ClassImp(TEvePointSet);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TEvePointSet::TEvePointSet(Int_t n_points, ETreeVarType_e tv_type) :
   TEveElement(),
   TPointSet3D(n_points),
   TEvePointSelectorConsumer(tv_type),
   TEveProjectable(),
   TQObject(),

   fTitle          (),
   fIntIds         (0),
   fIntIdsPerPoint (0)
{
   fMarkerStyle = 20;
   SetMainColorPtr(&fMarkerColor);

   // Override from TEveElement.
   fPickable = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TEvePointSet::TEvePointSet(const char* name, Int_t n_points, ETreeVarType_e tv_type) :
   TEveElement(),
   TPointSet3D(n_points),
   TEvePointSelectorConsumer(tv_type),
   TEveProjectable(),
   TQObject(),

   fTitle          (),
   fIntIds         (0),
   fIntIdsPerPoint (0)
{
   fMarkerStyle = 20;
   SetName(name);
   SetMainColorPtr(&fMarkerColor);

   // Override from TEveElement.
   fPickable = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.

TEvePointSet::TEvePointSet(const TEvePointSet& e) :
   TEveElement(e),
   TPointSet3D(e),
   TEvePointSelectorConsumer(e),
   TEveProjectable(),
   TQObject(),

   fTitle          (e.fTitle),
   fIntIds         (e.fIntIds ? new TArrayI(*e.fIntIds) : 0),
   fIntIdsPerPoint (e.fIntIdsPerPoint)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TEvePointSet::~TEvePointSet()
{
   delete fIntIds;
}

////////////////////////////////////////////////////////////////////////////////
/// Clone points and all point-related information from point-set 'e'.

void TEvePointSet::ClonePoints(const TEvePointSet& e)
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

   // TEvePointSet
   delete fIntIds;
   fIntIds         = e.fIntIds ? new TArrayI(*e.fIntIds) : 0;
   fIntIdsPerPoint = e.fIntIdsPerPoint;
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointset icon.

const TGPicture* TEvePointSet::GetListTreeIcon(Bool_t)
{
   return TEveElement::fgListTreeIcons[3];
}

////////////////////////////////////////////////////////////////////////////////
/// Drop all data and set-up the data structures to recive new data.
/// n_points   specifies the initial size of the arrays.
/// n_int_ids  specifies the number of integer ids per point.

void TEvePointSet::Reset(Int_t n_points, Int_t n_int_ids)
{
   delete [] fP; fP = 0;
   fN = n_points;
   if (fN) {
      fP = new Float_t [3*fN];
      memset(fP, 0, 3*fN*sizeof(Float_t));
   }
   fLastPoint = -1;
   ClearIds();
   delete fIntIds; fIntIds = 0;
   fIntIdsPerPoint = n_int_ids;
   if (fIntIdsPerPoint > 0) fIntIds = new TArrayI(fIntIdsPerPoint*fN);
   ResetBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Resizes internal array to allow additional n_points to be stored.
/// Returns the old size which is also the location where one can
/// start storing new data.
/// The caller is *obliged* to fill the new point slots.

Int_t TEvePointSet::GrowFor(Int_t n_points)
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

inline void TEvePointSet::AssertIntIdsSize()
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

Int_t* TEvePointSet::GetPointIntIds(Int_t p) const
{
   if (fIntIds)
      return fIntIds->GetArray() + p*fIntIdsPerPoint;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return i-th integer id of point with index p.
/// Existence of integer id array is checked, kMinInt is returned if
/// it does not exist.
/// Validity of p and i is *not* checked.

Int_t TEvePointSet::GetPointIntId(Int_t p, Int_t i) const
{
   if (fIntIds)
      return * (fIntIds->GetArray() + p*fIntIdsPerPoint + i);
   return kMinInt;
}

////////////////////////////////////////////////////////////////////////////////
/// Set integer ids for the last point that was registered (most
/// probably via TPolyMarker3D::SetNextPoint(x,y,z)).

void TEvePointSet::SetPointIntIds(Int_t* ids)
{
   SetPointIntIds(fLastPoint, ids);
}

////////////////////////////////////////////////////////////////////////////////
/// Set integer ids for point with index n.

void TEvePointSet::SetPointIntIds(Int_t n, Int_t* ids)
{
   if (!fIntIds) return;
   AssertIntIdsSize();
   Int_t* x = fIntIds->GetArray() + n*fIntIdsPerPoint;
   for (Int_t i=0; i<fIntIdsPerPoint; ++i)
      x[i] = ids[i];
}

////////////////////////////////////////////////////////////////////////////////
/// Set marker style, propagate to projecteds.

void TEvePointSet::SetMarkerStyle(Style_t mstyle)
{
   static const TEveException eh("TEvePointSet::SetMarkerStyle ");

   std::list<TEveProjected*>::iterator pi = fProjectedList.begin();
   while (pi != fProjectedList.end())
   {
      TEvePointSet* pt = dynamic_cast<TEvePointSet*>(*pi);
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

void TEvePointSet::SetMarkerSize(Size_t msize)
{
   static const TEveException eh("TEvePointSet::SetMarkerSize ");

   std::list<TEveProjected*>::iterator pi = fProjectedList.begin();
   while (pi != fProjectedList.end())
   {
      TEvePointSet* pt = dynamic_cast<TEvePointSet*>(*pi);
      if (pt)
      {
         pt->SetMarkerSize(msize);
         pt->StampObjProps();
      }
      ++pi;
   }
   TAttMarker::SetMarkerSize(msize);
}

////////////////////////////////////////////////////////////////////////////////
/// Paint point-set.

void TEvePointSet::Paint(Option_t*)
{
   PaintStandard(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize point-set for new filling.
/// subIdNum gives the number of integer ids that can be assigned to
/// each point.

void TEvePointSet::InitFill(Int_t subIdNum)
{
   if (subIdNum > 0) {
      fIntIdsPerPoint = subIdNum;
      if (!fIntIds)
         fIntIds = new TArrayI(fIntIdsPerPoint*GetN());
      else
         fIntIds->Set(fIntIdsPerPoint*GetN());
   } else {
      delete fIntIds; fIntIds = 0;
      fIntIdsPerPoint = 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Called from TEvePointSelector when internal arrays of the tree-selector
/// are filled up and need to be processed.
/// Virtual from TEvePointSelectorConsumer.

void TEvePointSet::TakeAction(TEvePointSelector* sel)
{
   static const TEveException eh("TEvePointSet::TakeAction ");

   if(sel == 0)
      throw(eh + "selector is <null>.");

   Int_t    n = sel->GetNfill();
   Int_t  beg = GrowFor(n);

   // printf("TEvePointSet::TakeAction beg=%d n=%d size=%d nsubid=%d dim=%d\n",
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

void TEvePointSet::CopyVizParams(const TEveElement* el)
{
   const TEvePointSet* m = dynamic_cast<const TEvePointSet*>(el);
   if (m)
   {
      TAttMarker::operator=(*m);
      fOption = m->fOption;
   }

   TEveElement::CopyVizParams(el);
}

////////////////////////////////////////////////////////////////////////////////
/// Write visualization parameters.

void TEvePointSet::WriteVizParams(std::ostream& out, const TString& var)
{
   TEveElement::WriteVizParams(out, var);

   TAttMarker::SaveMarkerAttributes(out, var);
}

////////////////////////////////////////////////////////////////////////////////
/// Virtual from TEveProjectable, returns TEvePointSetProjected class.

TClass* TEvePointSet::ProjectedClass(const TEveProjection*) const
{
   return TEvePointSetProjected::Class();
}

////////////////////////////////////////////////////////////////////////////////
/// Virtual method of base class TPointSet3D. The function call is
/// invoked with secondary selection in TPointSet3DGL.

void TEvePointSet::PointSelected(Int_t id)
{
   Emit("PointSelected(Int_t)", id);
   TPointSet3D::PointSelected(id);
}

//==============================================================================
/** \class TEvePointSetArray
\ingroup TEve
An array of point-sets with each point-set playing a role of a bin
in a histogram. When a new point is added to a TEvePointSetArray,
an additional separating quantity needs to be specified: it
determines into which TEvePointSet (bin) the point will actually be
stored. Underflow and overflow bins are automatically created but
they are not drawn by default.

By using the TEvePointSelector the points and the separating
quantities can be filled directly from a TTree holding the source
data.
Setting of per-point TRef's is not supported.

After the filling, the range of separating variable can be
controlled with a slider to choose a sub-set of PointSets that are
actually shown.
*/

ClassImp(TEvePointSetArray);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TEvePointSetArray::TEvePointSetArray(const char* name,
                                     const char* title) :
   TEveElement(),
   TNamed(name, title),

   fBins(0), fDefPointSetCapacity(128), fNBins(0), fLastBin(-1),
   fMin(0), fCurMin(0), fMax(0), fCurMax(0),
   fBinWidth(0),
   fQuantName()
{

   SetMainColorPtr(&fMarkerColor);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor: deletes the fBins array. Actual removal of
/// elements done by TEveElement.

TEvePointSetArray::~TEvePointSetArray()
{
   // printf("TEvePointSetArray::~TEvePointSetArray()\n");
   delete [] fBins; fBins = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Virtual from TEveElement, provide bin management.

void TEvePointSetArray::RemoveElementLocal(TEveElement* el)
{
   for (Int_t i=0; i<fNBins; ++i) {
      if (fBins[i] == el) {
         fBins[i] = 0;
         break;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Virtual from TEveElement, provide bin management.

void TEvePointSetArray::RemoveElementsLocal()
{
   delete [] fBins; fBins = 0; fLastBin = -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Set marker color, propagate to children.

void TEvePointSetArray::SetMarkerColor(Color_t tcolor)
{
   static const TEveException eh("TEvePointSetArray::SetMarkerColor ");

   for (List_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
      TAttMarker* m = dynamic_cast<TAttMarker*>((*i)->GetObject(eh));
      if (m && m->GetMarkerColor() == fMarkerColor)
         m->SetMarkerColor(tcolor);
   }
   TAttMarker::SetMarkerColor(tcolor);
}

////////////////////////////////////////////////////////////////////////////////
/// Set marker style, propagate to children.

void TEvePointSetArray::SetMarkerStyle(Style_t mstyle)
{
   static const TEveException eh("TEvePointSetArray::SetMarkerStyle ");

   for (List_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
      TAttMarker* m = dynamic_cast<TAttMarker*>((*i)->GetObject(eh));
      if (m && m->GetMarkerStyle() == fMarkerStyle)
         m->SetMarkerStyle(mstyle);
   }
   TAttMarker::SetMarkerStyle(mstyle);
}

////////////////////////////////////////////////////////////////////////////////
/// Set marker size, propagate to children.

void TEvePointSetArray::SetMarkerSize(Size_t msize)
{
   static const TEveException eh("TEvePointSetArray::SetMarkerSize ");

   for (List_i i=fChildren.begin(); i!=fChildren.end(); ++i) {
      TAttMarker* m = dynamic_cast<TAttMarker*>((*i)->GetObject(eh));
      if (m && m->GetMarkerSize() == fMarkerSize)
         m->SetMarkerSize(msize);
   }
   TAttMarker::SetMarkerSize(msize);
}

////////////////////////////////////////////////////////////////////////////////
/// Called from TEvePointSelector when internal arrays of the tree-selector
/// are filled up and need to be processed.
/// Virtual from TEvePointSelectorConsumer.

void TEvePointSetArray::TakeAction(TEvePointSelector* sel)
{
   static const TEveException eh("TEvePointSetArray::TakeAction ");

   if (sel == 0)
      throw eh + "selector is <null>.";

   Int_t n = sel->GetNfill();

   // printf("TEvePointSetArray::TakeAction n=%d\n", n);

   Double_t *vx = sel->GetV1(), *vy = sel->GetV2(), *vz = sel->GetV3();
   Double_t *qq = sel->GetV4();

   if (qq == 0)
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

Int_t TEvePointSetArray::Size(Bool_t under, Bool_t over) const
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

void TEvePointSetArray::InitBins(const char* quant_name,
                                 Int_t nbins, Double_t min, Double_t max)
{
   static const TEveException eh("TEvePointSetArray::InitBins ");

   if (nbins < 1) throw eh + "nbins < 1.";
   if (min > max) throw eh + "min > max.";

   RemoveElements();

   fQuantName = quant_name;
   fNBins     = nbins + 2; // under/overflow
   fLastBin   = -1;
   fMin = fCurMin = min;
   fMax = fCurMax = max;
   fBinWidth  = (fMax - fMin)/(fNBins - 2);

   fBins = new TEvePointSet* [fNBins];

   for (Int_t i = 0; i < fNBins; ++i)
   {
      fBins[i] = new TEvePointSet
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
/// If the selected bin does not have an associated TEvePointSet
/// the point is discarded and false is returned.

Bool_t TEvePointSetArray::Fill(Double_t x, Double_t y, Double_t z, Double_t quant)
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

void TEvePointSetArray::SetPointId(TObject* id)
{
   if (fLastBin >= 0)
      fBins[fLastBin]->SetPointId(id);
}

////////////////////////////////////////////////////////////////////////////////
/// Call this after all the points have been filled.
/// At this point we can calculate bounding-boxes of individual
/// point-sets.

void TEvePointSetArray::CloseBins()
{
   for (Int_t i=0; i<fNBins; ++i)
   {
      if (fBins[i] != 0)
      {
         fBins[i]->SetTitle(Form("N=%d", fBins[i]->Size()));
         fBins[i]->ComputeBBox();
      }
   }
   fLastBin = -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Propagate id-object ownership to children.

void TEvePointSetArray::SetOwnIds(Bool_t o)
{
   for (Int_t i=0; i<fNBins; ++i)
   {
      if (fBins[i] != 0)
         fBins[i]->SetOwnIds(o);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set active range of the separating quantity.
/// Appropriate point-sets are tagged for rendering.
/// Over/underflow point-sets are left as they were.

void TEvePointSetArray::SetRange(Double_t min, Double_t max)
{
   using namespace TMath;

   fCurMin = min; fCurMax = max;
   Int_t  low_b = Max(0,        FloorNint((min-fMin)/fBinWidth)) + 1;
   Int_t high_b = Min(fNBins-2, CeilNint ((max-fMin)/fBinWidth));

   for (Int_t i = 1; i < fNBins - 1; ++i)
   {
      if (fBins[i] != 0)
         fBins[i]->SetRnrSelf(i>=low_b && i<=high_b);
   }
}

/** \class TEvePointSetProjected
\ingroup TEve
Projected copy of a TEvePointSet.
*/
ClassImp(TEvePointSetProjected);

////////////////////////////////////////////////////////////////////////////////
/// Default contructor.

TEvePointSetProjected::TEvePointSetProjected() :
   TEvePointSet  (),
   TEveProjected ()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Set projection manager and projection model.
/// Virtual from TEveProjected.

void TEvePointSetProjected::SetProjection(TEveProjectionManager* proj,
                                          TEveProjectable* model)
{
   TEveProjected::SetProjection(proj, model);
   CopyVizParams(dynamic_cast<TEveElement*>(model));
}

////////////////////////////////////////////////////////////////////////////////
/// Set depth (z-coordinate) of the projected points.

void TEvePointSetProjected::SetDepthLocal(Float_t d)
{
   SetDepthCommon(d, this, fBBox);

   Int_t    n = Size();
   Float_t *p = GetP() + 2;
   for (Int_t i = 0; i < n; ++i, p+=3)
      *p = fDepth;
}

////////////////////////////////////////////////////////////////////////////////
/// Re-apply the projection.
/// Virtual from TEveProjected.

void TEvePointSetProjected::UpdateProjection()
{
   TEveProjection &proj = * fManager->GetProjection();
   TEvePointSet   &ps   = * dynamic_cast<TEvePointSet*>(fProjectable);
   TEveTrans      *tr   =   ps.PtrMainTrans(kFALSE);

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

void TEvePointSetProjected::PointSelected(Int_t id)
{
   TEvePointSet *ps = dynamic_cast<TEvePointSet*>(fProjectable);
   ps->PointSelected(id);
}
