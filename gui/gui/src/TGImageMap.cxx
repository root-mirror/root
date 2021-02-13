// @(#)root/gui:$Id$
// Author: Valeriy Onuchin & Fons Rademakers   18/10/2000

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TGImageMap (with TGRegion and TGRegionWithId help classes)           //
//                                                                      //
// A TGImageMap provides the functionality like a clickable image in    //
// a web browser with sensitive regions (MAP HTML tag).                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TGImageMap.h"
#include "TRefCnt.h"
#include "TGMenu.h"
#include "TGToolTip.h"
#include "TList.h"
#include "TArrayS.h"
#include "TVirtualX.h"


ClassImp(TGRegion);
ClassImp(TGRegionWithId);
ClassImpQ(TGImageMap)


TGRegionWithId  *gCurrentRegion; // current region

static TGRegion *gEmptyRegion = 0;
static Int_t gPointerX;  // current X mouse position
static Int_t gPointerY;  // current Y mouse position



class TGRegionData : public TRefCnt {

friend class TGRegion;

private:
   Region_t   fRgn;     // region handle
   Bool_t     fIsNull;  // true if null region

public:
   TGRegionData() { fRgn = 0; fIsNull = kTRUE; AddReference(); }
   ~TGRegionData() { }
   TGRegionData &operator=(const TGRegionData &r);
};

////////////////////////////////////////////////////////////////////////////////
/// Assignemnt of region data object.

TGRegionData &TGRegionData::operator=(const TGRegionData &r)
{
   if (this != &r) {
      fRefs   = r.fRefs;
      fRgn    = r.fRgn;
      fIsNull = r.fIsNull;
   }
   return *this;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/// Create a region object.

TGRegion::TGRegion()
{
   if (!gEmptyRegion)                      // avoid too many allocs
      gEmptyRegion = new TGRegion(kTRUE);

   fData = gEmptyRegion->fData;
   fData->AddReference();
}

////////////////////////////////////////////////////////////////////////////////
/// Create empty region.

TGRegion::TGRegion(Bool_t is_null)
{
   fData          = new TGRegionData;
   fData->fRgn    = gVirtualX->CreateRegion();
   fData->fIsNull = is_null;
}

////////////////////////////////////////////////////////////////////////////////
/// Create and initialize a region with a rectangle.

TGRegion::TGRegion(Int_t x, Int_t y, UInt_t w, UInt_t h, ERegionType)
{
   fData          = new TGRegionData;
   fData->fRgn    = gVirtualX->CreateRegion();
   fData->fIsNull = kFALSE;

   Rectangle_t xr;
   xr.fX      = (Short_t) x;
   xr.fY      = (Short_t) y;
   xr.fWidth  = (UShort_t) w;
   xr.fHeight = (UShort_t) h;
   gVirtualX->UnionRectWithRegion(&xr, fData->fRgn, fData->fRgn);
}

////////////////////////////////////////////////////////////////////////////////
/// Create and intialize a region with a polygon.

TGRegion::TGRegion(Int_t n, TPoint *points, Bool_t winding)
{
   fData            = new TGRegionData;
   fData->fIsNull   = kFALSE;
   Point_t *gpoints = new Point_t[n];

   for (int i = 0; i < n; i++) {
      gpoints[i].fX = (Short_t) points[i].GetX();
      gpoints[i].fY = (Short_t) points[i].GetY();
   }

   fData->fRgn = gVirtualX->PolygonRegion(gpoints, n, winding);
}

////////////////////////////////////////////////////////////////////////////////
/// Create and initialize a region with an X and a Y array of points.

TGRegion::TGRegion(const TArrayS &x, const TArrayS &y, Bool_t winding)
{
   fData          = new TGRegionData;
   fData->fIsNull = kFALSE;

   Int_t n = x.GetSize();
   if (n != y.GetSize()) {
      Error("TGRegion", "x and y arrays must have same length");
      return;
   }
   Point_t *gpoints = new Point_t[n];

   for (int i = 0; i < n; i++) {
      gpoints[i].fX = x.GetArray()[i];
      gpoints[i].fY = y.GetArray()[i];
   }

   fData->fRgn = gVirtualX->PolygonRegion(gpoints, n, winding);
}

////////////////////////////////////////////////////////////////////////////////
/// Create and initialize a region with an X and Y array of points.

TGRegion::TGRegion(Int_t n, Int_t *x, Int_t *y, Bool_t winding)
{
   fData          = new TGRegionData;
   fData->fIsNull = kFALSE;
   Point_t *gpoints = new Point_t[n];

   for (int i = 0; i < n; i++) {
      gpoints[i].fX = x[i];
      gpoints[i].fY = y[i];
   }

   fData->fRgn = gVirtualX->PolygonRegion(gpoints, n, winding);
}

////////////////////////////////////////////////////////////////////////////////
/// Region copy constructor.

TGRegion::TGRegion(const TGRegion &r) : TObject(r)
{
   fData = r.fData;
   fData->AddReference();
}

////////////////////////////////////////////////////////////////////////////////
/// Delete a region.

TGRegion::~TGRegion()
{
   if (fData->RemoveReference() <= 0) {
      gVirtualX->DestroyRegion(fData->fRgn);
      delete fData;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Region assignment operator.

TGRegion &TGRegion::operator=(const TGRegion &r)
{
   if (this != &r) {
      TObject::operator=(r);
      r.fData->AddReference();

      if (fData->RemoveReference() <= 0) {
         gVirtualX->DestroyRegion(fData->fRgn);
         delete fData;
      }
      fData = r.fData;
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy a region.

TGRegion TGRegion::CopyRegion() const
{
   TGRegion r(fData->fIsNull);
   gVirtualX->UnionRegion(fData->fRgn, r.fData->fRgn, r.fData->fRgn);
   return r;
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if region is not set.

Bool_t TGRegion::IsNull() const
{
   return fData->fIsNull;
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if region is empty.

Bool_t TGRegion::IsEmpty() const
{
   return fData->fIsNull || gVirtualX->EmptyRegion(fData->fRgn);
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if point p is contained in the region.

Bool_t TGRegion::Contains(const TPoint &p) const
{
   return gVirtualX->PointInRegion((Int_t)p.GetX(), (Int_t)p.GetY(), fData->fRgn);
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if point (x,y) is contained in the region.

Bool_t TGRegion::Contains(Int_t x, Int_t y) const
{
   return gVirtualX->PointInRegion(x, y, fData->fRgn);
}

////////////////////////////////////////////////////////////////////////////////
/// Return the union of this region with r.

TGRegion TGRegion::Unite(const TGRegion &r) const
{
   TGRegion result(kFALSE);
   gVirtualX->UnionRegion(fData->fRgn, r.fData->fRgn, result.fData->fRgn);
   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns a region which is the intersection of this region and r.

TGRegion TGRegion::Intersect(const TGRegion &r) const
{
   TGRegion result(kFALSE);
   gVirtualX->IntersectRegion(fData->fRgn, r.fData->fRgn, result.fData->fRgn);
   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns a region which is r subtracted from this region.

TGRegion TGRegion::Subtract(const TGRegion &r) const
{
   TGRegion result(kFALSE);
   gVirtualX->SubtractRegion(fData->fRgn, r.fData->fRgn, result.fData->fRgn);
   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns a region which is the difference between the union and
/// intersection this region and r.

TGRegion TGRegion::Eor(const TGRegion &r) const
{
   TGRegion result(kFALSE);
   gVirtualX->XorRegion(fData->fRgn, r.fData->fRgn, result.fData->fRgn);
   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Return dimension of region (widht, height).

TGDimension TGRegion::GetDimension() const
{
   Rectangle_t r = { 0, 0, 0, 0 };
   gVirtualX->GetRegionBox(fData->fRgn, &r);
   return TGDimension(r.fWidth, r.fHeight);
}

////////////////////////////////////////////////////////////////////////////////
/// Return position of region (x, y).

TGPosition TGRegion::GetPosition() const
{
   Rectangle_t r = { 0, 0, 0, 0 };
   gVirtualX->GetRegionBox(fData->fRgn, &r);
   return TGPosition(r.fX, r.fY);
}

////////////////////////////////////////////////////////////////////////////////
/// Region == operator.

Bool_t TGRegion::operator==(const TGRegion &r) const
{
   return fData == r.fData ?
             kTRUE : gVirtualX->EqualRegion(fData->fRgn, r.fData->fRgn);
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/// Create GUI region (with id and possible tooltip).

TGRegionWithId::TGRegionWithId() : TGRegion()
{
   fId    = 0;
   fTip   = 0;
   fPopup = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Create GUI region (with id and possible tooltip).

TGRegionWithId::TGRegionWithId(Int_t id, Int_t x, Int_t y,
                               UInt_t w, UInt_t h, ERegionType type) :
   TGRegion(x, y, w, h, type)
{
   fId    = id;
   fTip   = 0;
   fPopup = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Create GUI region (with id and possible tooltip).

TGRegionWithId::TGRegionWithId(Int_t id, Int_t n, TPoint *points,
                               Bool_t winding) :
   TGRegion(n, points, winding)
{
   fId    = id;
   fTip   = 0;
   fPopup = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.

TGRegionWithId::TGRegionWithId(const TGRegionWithId &reg) : TGRegion(reg)
{
   fId    = reg.GetId();
   fTip   = 0;
   fPopup = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy ctor which allows setting of new id.

TGRegionWithId::TGRegionWithId(const TGRegion &reg, Int_t id) :
   TGRegion(reg)
{
   fId    = id;
   fTip   = 0;
   fPopup = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Cleanup.

TGRegionWithId::~TGRegionWithId()
{
   delete fTip;
}

////////////////////////////////////////////////////////////////////////////////
/// Display popup menu associated with this region.

void TGRegionWithId::DisplayPopup()
{
   if (fPopup) fPopup->PlaceMenu(gPointerX, gPointerY, kTRUE, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Set tool tip text associated with this region. The delay is in
/// milliseconds (minimum 250). To remove tool tip call method with
/// text = 0.

void TGRegionWithId::SetToolTipText(const char *text, Long_t delayms,
                                    const TGFrame *frame)
{
   if (fTip) {
      delete fTip;
      fTip = 0;
   }

   if (text && strlen(text))
      fTip = new TGToolTip(gClient->GetDefaultRoot(), frame, text, delayms);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/// Create an image map widget.

TGImageMap::TGImageMap(const TGWindow *p, const TGPicture *pic) :
   TGPictureButton(p, pic)
{
   fCursorMouseOut  = kPointer;
   fCursorMouseOver = kHand;
   fListOfRegions   = new TList;
   fTrash           = new TList;
   fMainTip         = 0;
   fLastVisited     = 0;
   fNavMode = kNavRegions;

   SetDisabledPicture(fPic);
   SetState(kButtonDisabled);

   gVirtualX->GrabButton(fId, kAnyButton, kAnyModifier,
                         kButtonPressMask | kButtonReleaseMask |
                         kPointerMotionMask, kNone, kNone);

   AddInput(kKeyPressMask | kKeyReleaseMask | kPointerMotionMask |
            kStructureNotifyMask | kLeaveWindowMask);
   SetWindowName();
}

////////////////////////////////////////////////////////////////////////////////
/// Create an image map widget.

TGImageMap::TGImageMap(const TGWindow *p, const TString &pic) :
   TGPictureButton(p, pic.Data())
{
   fCursorMouseOut  = kPointer;
   fCursorMouseOver = kHand;
   fListOfRegions   = new TList;
   fTrash           = new TList;
   fMainTip         = 0;
   fLastVisited     = 0;
   fNavMode = kNavRegions;

   SetDisabledPicture(fPic);
   SetState(kButtonDisabled);

   gVirtualX->GrabButton(fId, kAnyButton, kAnyModifier,
                         kButtonPressMask | kButtonReleaseMask |
                         kPointerMotionMask, kNone, kNone);

   AddInput(kKeyPressMask | kKeyReleaseMask | kPointerMotionMask |
            kStructureNotifyMask | kLeaveWindowMask);
   SetWindowName();
}

////////////////////////////////////////////////////////////////////////////////
/// Cleanup image map widget.

TGImageMap::~TGImageMap()
{
   delete fMainTip;
   fTrash->Delete();
   delete fTrash;
   fListOfRegions->Delete();
   delete fListOfRegions;
}

////////////////////////////////////////////////////////////////////////////////
/// Add a region to the image map.

void TGImageMap::AddRegion(const TGRegion &region, Int_t id)
{
   fListOfRegions->Add(new TGRegionWithId(region, id));
}

////////////////////////////////////////////////////////////////////////////////
/// Create popoup menu or returns existing for regions with specified id.

TGPopupMenu *TGImageMap::CreatePopup(Int_t id)
{
   TIter next(fListOfRegions);
   TGRegionWithId *region;
   TGPopupMenu    *popup = 0;
   TGPopupMenu    *newpopup = 0;

   while ((region = (TGRegionWithId*)next())) {
      if (id == region->GetId()) {
         popup = region->GetPopup();
         if (!popup && !newpopup) {
            newpopup = new TGPopupMenu(this);
            fTrash->Add(newpopup);
         }
         if (newpopup) region->SetPopup(newpopup);
      }
   }
   return newpopup ? newpopup : popup;
}

////////////////////////////////////////////////////////////////////////////////
/// Return popup for regions with specified id.

TGPopupMenu *TGImageMap::GetPopup(Int_t id)
{
   TIter next(fListOfRegions);
   TGRegionWithId *region;

   while ((region = (TGRegionWithId*)next())) {
      if (id == region->GetId()) return region->GetPopup();
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle mouse motion events.

Bool_t TGImageMap::HandleMotion(Event_t *event)
{
   TIter next(fListOfRegions);
   TGRegionWithId *region;

   if (fNavMode != kNavRegions) return kTRUE;
   gPointerX = event->fX;
   gPointerY = event->fY;

   while ((region = (TGRegionWithId*)next())) {
      if (region->Contains(gPointerX, gPointerY)) {
         if (fLastVisited == region->GetId()) return kTRUE;
         if (fLastVisited) OnMouseOut(fLastVisited);
         fLastVisited = region->GetId();
         fTip = region->GetToolTipText();
         gCurrentRegion = region;
         OnMouseOver(fLastVisited);
         return kTRUE;
      }
   }

   if (fLastVisited) {
      OnMouseOut(fLastVisited);
      fTip = fMainTip;
   }
   fLastVisited = 0;  // main
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle double click events.

Bool_t TGImageMap::HandleDoubleClick(Event_t *event)
{
   TIter next(fListOfRegions);
   TGRegionWithId *region;

   if (fTip) fTip->Hide();
   if (event->fCode != kButton1 ) return kTRUE;
   if (fNavMode != kNavRegions) return kTRUE;

   gPointerX = event->fX;
   gPointerY = event->fY;

   while ((region = (TGRegionWithId*)next())) {
      if (region->Contains(gPointerX, gPointerY)) {
         DoubleClicked(region->GetId());
         gCurrentRegion = region;
         return kTRUE;
      }
   }
   DoubleClicked();
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Handle button events.

Bool_t TGImageMap::HandleButton(Event_t *event)
{
   TIter next(fListOfRegions);
   TGRegionWithId *region;
   TGPopupMenu *pop;

   if (fTip) fTip->Hide();
   if (fNavMode != kNavRegions) return kTRUE;

   gPointerX = event->fX;
   gPointerY = event->fY;

   while ((region = (TGRegionWithId*)next())) {
      if (region->Contains(gPointerX, gPointerY)) {
         gCurrentRegion = region;
         if (event->fType == kButtonPress) {
            if (event->fCode == kButton1 )
               RegionClicked(region->GetId());
            else if (event->fCode == kButton3 ) {
               pop = region->GetPopup();
               if (pop) pop->PlaceMenu(gPointerX, gPointerY, kTRUE, kTRUE);
            }
         }
         return kTRUE;
      }
   }
   if (event->fType == kButtonPress)
      Clicked();
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set tooltip text for main region.

void TGImageMap::SetToolTipText(const char *text, Long_t delayms)
{
   if (fMainTip) delete fMainTip;
   fMainTip = 0;

   if (text && strlen(text))
      fMainTip = new TGToolTip(fClient->GetDefaultRoot(), this, text, delayms);
}

////////////////////////////////////////////////////////////////////////////////
/// Set tooltip text for regions with specified id.

void TGImageMap::SetToolTipText(Int_t id, const char *text, Long_t delayms)
{
   TIter next(fListOfRegions);
   TGRegionWithId *region;

   while ((region = (TGRegionWithId*)next())) {
      if (id == region->GetId())
         region->SetToolTipText(text, delayms, this);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Handle when mouse moves over region id. Emits signal
/// OnMouseOver(Int_t).

void TGImageMap::OnMouseOver(Int_t id)
{
   if (fTip) fTip->Reset();
   if (fMainTip) fMainTip->Hide();
   gVirtualX->SetCursor(fId, gVirtualX->CreateCursor(fCursorMouseOver));
   Emit("OnMouseOver(Int_t)", id);
}

////////////////////////////////////////////////////////////////////////////////
/// Handle when mouse moves from region id. Emits signal
/// OnMouseOut(Int_t).

void TGImageMap::OnMouseOut(Int_t id)
{
   if(fTip) fTip->Hide();
   if(fMainTip) fMainTip->Reset();
   gVirtualX->SetCursor(fId,gVirtualX->CreateCursor(fCursorMouseOut));
   Emit("OnMouseOut(Int_t)",id);
}

////////////////////////////////////////////////////////////////////////////////
/// Handle when mouse was clicked on region id. Emits signal
/// RegionClicked(Int_t).

void TGImageMap::RegionClicked(Int_t id)
{
   Emit("RegionClicked(Int_t)",id);
}

////////////////////////////////////////////////////////////////////////////////
/// Handle when mouse is double clicked on main map. Emits signal
/// DoubleClicked().

void TGImageMap::DoubleClicked()
{
   Emit("DoubleClicked()");
}

////////////////////////////////////////////////////////////////////////////////
/// Handle when mouse is double clicked on region id. Emits signal
/// DoubleClicked(Int_t).

void TGImageMap::DoubleClicked(Int_t id)
{
   Emit("DoubleClicked(Int_t)",id);
}
