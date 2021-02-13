// @(#)root/eve:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TEveDigitSet.h"
#include "TEveManager.h"
#include "TEveTrans.h"

#include "TRefArray.h"


/** \class TEveDigitSet
\ingroup TEve
Base-class for storage of digit collections; provides
transformation matrix (TEveTrans), signal to color mapping
(TEveRGBAPalette) and visual grouping (TEveFrameBox).

Base-class for displaying a digit collection.
Provides common services for:
 - specifying signal / color per digit;
 - specifying object reference per digit;
 - controlling palette and thresholds (external object TEveRGBAPalette);
 - showing a frame around the digits (external object TEveFrameBox);
 - specifying transformation matrix for the whole collection;
   by data-member of class TEveTrans.

Use method DigitId(TObject* id) to assign additional identification
to the last created digit. By calling SetOwnIds(kTRUE) tje
digit-set becomes the owner of the assigned objects and deletes
them on destruction.
Note that TRef is used for referencing the objects and if you
instantiate the objects just to pass them to digit-set you should
also call  TProcessID::Get/SetObjectCount() at the beginning / end
of processing of an event. See documentation for class TRef, in
particular section 'ObjectNumber'.

If you use value-is-color mode and want to use transparency, set
the transparency to non-zero value so that GL-renderer will be
properly informed.

If you want to use single color for all elements call:
~~~ {.cpp}
   UseSingleColor()
~~~
Palette controls will not work in this case.

A pointer to a rectangle / box of class TEveFrameBox can be set via
~~~ {.cpp}
   void SetFrame(TEveFrameBox* b);
~~~
A single TEveFrameBox can be shared among several digit-sets (it is
reference-counted). The following flags affect how the frame-box will drawn
and used for selection and highlight:
~~~ {.cpp}
   Bool_t fSelectViaFrame;
   Bool_t fHighlightFrame;
~~~
TEveDigitSet is sub-classed from TEveSecondarySelectable -- this means
individual digits can be selected. By calling:
~~~ {.cpp}
   TEveSecondarySelectable::SetAlwaysSecSelect(kTRUE);
~~~
one can enforce immediate feedback (highlight, tooltip and select on normal
left-mouse click) on given digit-set.

See also:
~~~ {.cpp}
   TEveQuadSet: rectangle, hexagon or line per digit
   TEveBoxSet   a 3D box per digit
~~~
*/

ClassImp(TEveDigitSet);

////////////////////////////////////////////////////////////////////////////////

TEveDigitSet::TEveDigitSet(const char* n, const char* t) :
   TEveElement     (fColor),
   TNamed          (n, t),

   fDigitIds       (0),
   fDefaultValue   (kMinInt),
   fValueIsColor   (kFALSE),
   fSingleColor    (kFALSE),
   fAntiFlick      (kTRUE),
   fOwnIds         (kFALSE),
   fPlex           (),
   fLastDigit      (0),
   fLastIdx        (-1),

   fColor          (kWhite),
   fFrame          (0),
   fPalette        (0),
   fRenderMode     (kRM_AsIs),
   fSelectViaFrame (kFALSE),
   fHighlightFrame (kFALSE),
   fDisableLighting(kTRUE),
   fHistoButtons   (kTRUE),
   fEmitSignals    (kFALSE),
   fCallbackFoo    (0),
   fTooltipCBFoo   (0)
{
   // Constructor.

   fCanEditMainColor        = kTRUE;
   fCanEditMainTransparency = kTRUE;
   InitMainTrans();
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.
/// Unreference frame and palette. Destroy referenced objects if they
/// are owned by the TEveDigitSet.

TEveDigitSet::~TEveDigitSet()
{
   SetFrame(0);
   SetPalette(0);
   if (fOwnIds)
      ReleaseIds();
   delete fDigitIds;
}

////////////////////////////////////////////////////////////////////////////////
/// Protected method called whenever a new digit is added.

TEveDigitSet::DigitBase_t* TEveDigitSet::NewDigit()
{
   fLastIdx   = fPlex.Size();
   fLastDigit = new (fPlex.NewAtom()) DigitBase_t(fDefaultValue);
   return fLastDigit;
}

////////////////////////////////////////////////////////////////////////////////
/// Protected method. Release and delete the referenced objects, the
/// ownership is *NOT* checked.

void TEveDigitSet::ReleaseIds()
{
   if (fDigitIds)
   {
      const Int_t N = fDigitIds->GetSize();

      for (Int_t i = 0; i < N; ++i)
         delete fDigitIds->At(i);

      fDigitIds->Expand(0);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Instruct digit-set to use single color for its digits.
/// Call SetMainColor/Transparency to initialize it.

void TEveDigitSet::UseSingleColor()
{
   fSingleColor = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Override from TEveElement, forward to Frame.

void TEveDigitSet::SetMainColor(Color_t color)
{
   if (fSingleColor)
   {
      TEveElement::SetMainColor(color);
   }
   else if (fFrame)
   {
      fFrame->SetFrameColor(color);
      fFrame->StampBackPtrElements(kCBColorSelection);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Virtual function called when both fSelected is false and
/// fImpliedSelected is 0.

void TEveDigitSet::UnSelected()
{
   fSelectedSet.clear();
   TEveElement::UnSelected();
}

////////////////////////////////////////////////////////////////////////////////
/// Virtual function called when both fHighlighted is false and
/// fImpliedHighlighted is 0.

void TEveDigitSet::UnHighlighted()
{
   fHighlightedSet.clear();
   TEveElement::UnHighlighted();
}

////////////////////////////////////////////////////////////////////////////////
/// Return tooltip for highlighted element if always-sec-select is set.
/// Otherwise return the tooltip for this element.

TString TEveDigitSet::GetHighlightTooltip()
{
   if (fHighlightedSet.empty()) return "";

   if (GetAlwaysSecSelect())
   {
      if (fTooltipCBFoo)
      {
         return (fTooltipCBFoo)(this, *fHighlightedSet.begin());
      }
      else if (fDigitIds)
      {
         TObject *o = GetId(*fHighlightedSet.begin());
         if (o)
            return TString(o->GetName());
      }
      return TString::Format("%s; idx=%d", GetElementName(), *fHighlightedSet.begin());
   }
   else
   {
      return TEveElement::GetHighlightTooltip();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Instruct underlying memory allocator to regroup itself into a
/// contiguous memory chunk.

void TEveDigitSet::RefitPlex()
{
   fPlex.Refit();
}

////////////////////////////////////////////////////////////////////////////////
/// Iterate over the digits and determine min and max signal values.

void TEveDigitSet::ScanMinMaxValues(Int_t& min, Int_t& max)
{
   if (fValueIsColor || fPlex.Size() == 0)
   {
      min = max = 0;
      return;
   }

   min = kMaxInt;
   max = kMinInt;
   for (Int_t c=0; c<fPlex.VecSize(); ++c)
   {
      Char_t* a = fPlex.Chunk(c);
      Int_t   n = fPlex.NAtoms(c);
      while (n--)
      {
         Int_t v = ((DigitBase_t*)a)->fValue;
         if (v < min) min = v;
         if (v > max) max = v;
         a += fPlex.S();
      }
   }
   if (min == max)
      --min;
}

////////////////////////////////////////////////////////////////////////////////
/// Set current digit -- the one that will receive calls to
/// DigitValue/Color/Id/UserData() functions.
/// Note that various AddXyzz() functions set the current digit to the newly
/// added one.

void TEveDigitSet::SetCurrentDigit(Int_t idx)
{
   fLastIdx   = idx;
   fLastDigit = GetDigit(idx);
}

////////////////////////////////////////////////////////////////////////////////
/// Set signal value for the last digit added.

void TEveDigitSet::DigitValue(Int_t value)
{
   fLastDigit->fValue = value;
}

////////////////////////////////////////////////////////////////////////////////
/// Set color for the last digit added.

void TEveDigitSet::DigitColor(Color_t ci)
{
   TEveUtil::ColorFromIdx(ci, (UChar_t*) & fLastDigit->fValue, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Set color for the last digit added.

void TEveDigitSet::DigitColor(Color_t ci, Char_t transparency)
{
   TEveUtil::ColorFromIdx(ci, (UChar_t*) & fLastDigit->fValue, transparency);
}

////////////////////////////////////////////////////////////////////////////////
/// Set color for the last digit added.

void TEveDigitSet::DigitColor(UChar_t r, UChar_t g, UChar_t b, UChar_t a)
{
   UChar_t* x = (UChar_t*) & fLastDigit->fValue;
   x[0] = r; x[1] = g; x[2] = b; x[3] = a;
}

////////////////////////////////////////////////////////////////////////////////
/// Set color for the last digit added.

void TEveDigitSet::DigitColor(UChar_t* rgba)
{
   UChar_t* x = (UChar_t*) & fLastDigit->fValue;
   x[0] = rgba[0]; x[1] = rgba[1]; x[2] = rgba[2]; x[3] = rgba[3];
}

////////////////////////////////////////////////////////////////////////////////
/// Set external object reference for the last digit added.

void TEveDigitSet::DigitId(TObject* id)
{
   DigitId(fLastIdx, id);
}

////////////////////////////////////////////////////////////////////////////////
/// Set user-data for the last digit added.

void TEveDigitSet::DigitUserData(void* ud)
{
   fLastDigit->fUserData = ud;
}

////////////////////////////////////////////////////////////////////////////////
/// Set external object reference for digit n.

void TEveDigitSet::DigitId(Int_t n, TObject* id)
{
   if (!fDigitIds)
      fDigitIds = new TRefArray;

   if (fOwnIds && n < fDigitIds->GetSize() && fDigitIds->At(n))
      delete fDigitIds->At(n);

   fDigitIds->AddAtAndExpand(id, n);
}

////////////////////////////////////////////////////////////////////////////////
/// Set user-data for digit n.

void TEveDigitSet::DigitUserData(Int_t n, void* ud)
{
   GetDigit(n)->fUserData = ud;
}

////////////////////////////////////////////////////////////////////////////////
/// Return external TObject associated with digit n.

TObject* TEveDigitSet::GetId(Int_t n) const
{
   return fDigitIds ? fDigitIds->At(n) : 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Get user-data associated with digit n.

void* TEveDigitSet::GetUserData(Int_t n) const
{
   return GetDigit(n)->fUserData;
}

////////////////////////////////////////////////////////////////////////////////
/// Paint this object. Only direct rendering is supported.

void TEveDigitSet::Paint(Option_t*)
{
   PaintStandard(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Called from renderer when a digit with index idx is selected.
/// This is by-passed when always-secondary-select is active.

void TEveDigitSet::DigitSelected(Int_t idx)
{
   DigitBase_t *qb  = GetDigit(idx);
   TObject     *obj = GetId(idx);

   if (fCallbackFoo) {
      (fCallbackFoo)(this, idx, obj);
   }
   if (fEmitSignals) {
      SecSelected(this, idx);
   } else {
      printf("TEveDigitSet::DigitSelected idx=%d, value=%d, obj=0x%lx\n",
             idx, qb->fValue, (ULong_t)obj);
      if (obj)
         obj->Print();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Emit a SecSelected signal.
/// This is by-passed when always-secondary-select is active.

void TEveDigitSet::SecSelected(TEveDigitSet* qs, Int_t idx)
{
   Long_t args[2];
   args[0] = (Long_t) qs;
   args[1] = (Long_t) idx;

   Emit("SecSelected(TEveDigitSet*, Int_t)", args);
}

////////////////////////////////////////////////////////////////////////////////
/// Set TEveFrameBox pointer.

void TEveDigitSet::SetFrame(TEveFrameBox* b)
{
   if (fFrame == b) return;
   if (fFrame) fFrame->DecRefCount(this);
   fFrame = b;
   if (fFrame) {
      fFrame->IncRefCount(this);
      if (!fSingleColor) {
         SetMainColorPtr(fFrame->PtrFrameColor());
      }
   } else {
      SetMainColorPtr(&fColor);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set TEveRGBAPalette pointer.

void TEveDigitSet::SetPalette(TEveRGBAPalette* p)
{
   if (fPalette == p) return;
   if (fPalette) fPalette->DecRefCount();
   fPalette = p;
   if (fPalette) fPalette->IncRefCount();
}

////////////////////////////////////////////////////////////////////////////////
/// Make sure the TEveRGBAPalette pointer is not null.
/// If it is not set, a new one is instantiated and the range is set
/// to current min/max signal values.

TEveRGBAPalette* TEveDigitSet::AssertPalette()
{
   if (fPalette == 0) {
      fPalette = new TEveRGBAPalette;
      if (!fValueIsColor) {
         Int_t min, max;
         ScanMinMaxValues(min, max);
         fPalette->SetLimits(min, max);
         fPalette->SetMinMax(min, max);
      }
   }
   return fPalette;
}
