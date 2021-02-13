// @(#)root/eve7:$Id$
// Author: Matevz Tadel 2007, 2018

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/REveCompound.hxx>

#include "TClass.h"

using namespace ROOT::Experimental;

/** \class REveCompound
\ingroup REve
Description of REveCompound
*/

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

REveCompound::REveCompound(const std::string& n, const std::string& t,
                           Bool_t doColor, Bool_t doTransparency) :
   REveElement     (n, t),
   fCompoundOpen   (0),
   fDoColor        (doColor),
   fDoTransparency (doTransparency)
{
   if (fDoColor)
   {
      SetupDefaultColorAndTransparency(0, fDoColor, fDoTransparency);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// SetMainColor for the compound.
/// The color is also propagated to children with compound set to this
/// whose current color is the same as the old color.
///
/// The following CompoundSelectionColorBits have further influence:
///  - kCSCBApplyMainColorToAllChildren      - apply color to all children;
///  - kCSCBApplyMainColorToMatchingChildren - apply color to children who have
///                                            matching old color.

void REveCompound::SetMainColor(Color_t color)
{
   if (!fDoColor) return;

   Color_t old_color = GetMainColor();

   REveElement::SetMainColor(color);

   Bool_t color_all      = TestCSCBits(kCSCBApplyMainColorToAllChildren);
   Bool_t color_matching = TestCSCBits(kCSCBApplyMainColorToMatchingChildren);

   for (auto &c: fChildren)
   {
      if (color_all || (color_matching && c->GetMainColor() == old_color) ||
          (c->GetCompound() == this && c->GetMainColor() == old_color))
      {
         c->SetMainColor(color);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// SetMainTransparency for the compound.
/// The transparency is also propagated to children with compound set to this
/// whose current transparency is the same as the old transparency.
///
/// The following CompoundSelectionColorBits have further influence:
///  - kCSCBApplyMainTransparencyToAllChildren      - apply transparency to all children;
///  - kCSCBApplyMainTransparencyToMatchingChildren - apply transparency to children who have
///                                                   matching transparency.

void REveCompound::SetMainTransparency(Char_t t)
{
   if ( ! fDoTransparency) return;

   Char_t old_t = GetMainTransparency();

   REveElement::SetMainTransparency(t);

   Bool_t chg_all      = TestCSCBits(kCSCBApplyMainTransparencyToAllChildren);
   Bool_t chg_matching = TestCSCBits(kCSCBApplyMainTransparencyToMatchingChildren);

   for (auto &c: fChildren) {
      if (chg_all || (chg_matching && c->GetMainTransparency() == old_t) ||
          (c->GetCompound() == this && c->GetMainTransparency() == old_t))
      {
         c->SetMainTransparency(t);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Call base-class implementation.
/// If compound is open and compound of the new element is not set,
/// the el's compound is set to this.
/// You might also want to call RecheckImpliedSelections().

void REveCompound::AddElement(REveElement* el)
{
   REveElement::AddElement(el);
   if (IsCompoundOpen() && el->GetCompound() == nullptr)
      el->SetCompound(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Decompoundofy el, call base-class version.

void REveCompound::RemoveElementLocal(REveElement* el)
{
   if (el->GetCompound() == this)
      el->SetCompound(0);

   REveElement::RemoveElementLocal(el);
}

////////////////////////////////////////////////////////////////////////////////
/// Decompoundofy children, call base-class version.

void REveCompound::RemoveElementsLocal()
{
   for (auto &c: fChildren) {
      if (c->GetCompound() == this)
         c->SetCompound(nullptr);
   }

   REveElement::RemoveElementsLocal();
}

////////////////////////////////////////////////////////////////////////////////
/// Recurse on all children that are in this compound and
/// call the base-class version.
/// If SelectionColorBit kSCBImplySelectAllChildren is set, then all
/// children are added to the set.
///
/// Note that projected replicas of the compound will be added to
/// the set in base-class function that handles projectable.

void REveCompound::FillImpliedSelectedSet(Set_t& impSelSet)
{
   Bool_t select_all = TestCSCBits(kCSCBImplySelectAllChildren);

   for (auto &c: fChildren) {
      if (select_all || c->GetCompound() == this) {
         if (impSelSet.insert(c).second)
            c->FillImpliedSelectedSet(impSelSet);
      }
   }

   REveElement::FillImpliedSelectedSet(impSelSet);
}

////////////////////////////////////////////////////////////////////////////////
/// Virtual from REveProjectable, returns REveCompoundProjected class.

TClass* REveCompound::ProjectedClass(const REveProjection*) const
{
   return TClass::GetClass<REveCompoundProjected>();
}


/** \class REveCompoundProjected
\ingroup REve
Description of REveCompoundProjected
*/

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

REveCompoundProjected::REveCompoundProjected() :
   REveCompound  (),
   REveProjected ()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Revert back to the behaviour of REveElement as color
/// is propagated:
///  1. from projectable -> projected
///  2. from compound -> compound elements
/// and we do not need to do this twice for projected-compound-elements.

void REveCompoundProjected::SetMainColor(Color_t color)
{
   REveElement::SetMainColor(color);
}
