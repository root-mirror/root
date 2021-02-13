// @(#)root/eve7:$Id$
// Author: Matevz Tadel 2007

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/REveSelection.hxx>
#include <ROOT/REveProjectionBases.hxx>
#include <ROOT/REveCompound.hxx>
#include <ROOT/REveManager.hxx>
#include <ROOT/REveSecondarySelectable.hxx>

#include "TClass.h"
#include "TColor.h"

#include <iostream>

using namespace ROOT::Experimental;
namespace REX = ROOT::Experimental;

/** \class REveSelection
\ingroup REve
Make sure there is a SINGLE running REveSelection for each
selection type (select/highlight).
*/

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

REveSelection::REveSelection(const std::string& n, const std::string& t,
                             Color_t col_visible, Color_t col_hidden) :
   REveElement       (n, t),
   fVisibleEdgeColor (col_visible),
   fHiddenEdgeColor  (col_hidden),
   fActive           (kTRUE),
   fIsMaster         (kTRUE)
{
   // Managing complete selection state on element level.
   //
   // Method pointers for propagation of selected / implied selected state
   // to elements. This has to be done differently now -- and kept within
   // REveSelection.
   //
   // Also, see REveManager::PreDeleteElement. We might need some sort of
   // implied-selected-count after all (global, for all selections,
   // highlights) ... and traverse all selections if the element gets zapped.
   // Yup, we have it ...
   // XXXX but ... we can also go up to master and check there directly !!!!!

   AddPickToSelect(kPS_Master);
   AddPickToSelect(kPS_PableCompound);
   AddPickToSelect(kPS_Element);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

REveSelection::~REveSelection()
{
   DeactivateSelection();
   RemoveNieces();
}

////////////////////////////////////////////////////////////////////////////////
/// Set visible highlight color

void REveSelection::SetVisibleEdgeColorRGB(UChar_t r, UChar_t g, UChar_t b)
{
   fVisibleEdgeColor = TColor::GetColor(r, g, b);
   StampObjProps();
}

////////////////////////////////////////////////////////////////////////////////
/// Set hidden highlight color
void REveSelection::SetHiddenEdgeColorRGB(UChar_t r, UChar_t g, UChar_t b)
{
   fHiddenEdgeColor = TColor::GetColor(r, g, b);
   StampObjProps();
}

////////////////////////////////////////////////////////////////////////////////
/// Set to 'highlight' mode.

void REveSelection::SetHighlightMode()
{
   // Most importantly, this sets the pointers-to-function-members in
   // REveElement that are used to mark elements as (un)selected and
   // implied-(un)selected.

   fIsMaster     = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Select element indicated by the entry and fill its
/// implied-selected set.

void REveSelection::DoElementSelect(SelMap_i &entry)
{
   Set_t &imp_set = entry->second.f_implied;

   entry->first->FillImpliedSelectedSet(imp_set);

   auto i = imp_set.begin();
   while (i != imp_set.end())
   {
      if ((*i)->GetElementId() == 0)
      {
         if (gDebug > 0)
         {
            Info("REveSelection::DoElementSelect",
                 "Element '%s' [%s] with 0 id detected and removed.",
                 (*i)->GetCName(), (*i)->IsA()->GetName());
         }
         auto j = i++;
         imp_set.erase(j);
      }
      else
      {
         (*i)->IncImpliedSelected();
         ++i;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Deselect element indicated by the entry and clear its
/// implied-selected set.

void REveSelection::DoElementUnselect(SelMap_i &entry)
{
   Set_t &imp_set = entry->second.f_implied;

   for (auto &imp_el: imp_set) imp_el->DecImpliedSelected();

   imp_set.clear();
}

////////////////////////////////////////////////////////////////////////////////
/// Check if elemenet el is selected (not implied selected).

bool REveSelection::HasNiece(REveElement *el) const
{
   return fMap.find(el) != fMap.end();
}

////////////////////////////////////////////////////////////////////////////////
/// Check if any elements are selected.

bool REveSelection::HasNieces() const
{
   return ! fMap.empty();
}

////////////////////////////////////////////////////////////////////////////////
/// Pre-addition check. Deny addition if el is already selected.
/// Virtual from REveAunt.

bool REveSelection::AcceptNiece(REveElement* el)
{
   return el != this && fMap.find(el) == fMap.end() &&
          el->IsA()->InheritsFrom(TClass::GetClass<REveSelection>()) == kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Add an element into selection, virtual from REveAunt

void REveSelection::AddNieceInternal(REveElement* el)
{
   auto res = fMap.emplace(el, Record(el));
   if (fActive) {
      DoElementSelect(res.first);
      SelectionAdded(el);
   }
   StampObjPropsPreChk();
}

////////////////////////////////////////////////////////////////////////////////
/// Virtual from REveAunt.

void REveSelection::RemoveNieceInternal(REveElement* el)
{
   auto i = fMap.find(el);

   if (i != fMap.end())
   {
      if (fActive)
      {
         DoElementUnselect(i);
         SelectionRemoved(el);
      }
      fMap.erase(i);
      StampObjPropsPreChk();
   }
   else
   {
      Warning("REveSelection::RemoveNieceLocal", "element not found in map.");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Add an element into selection, virtual from REveAunt.
/// Overriden here just so that a signal can be emitted.

void REveSelection::RemoveNieces()
{
   if (IsEmpty()) return;

   for (auto i = fMap.begin(); i != fMap.end(); ++i)
   {
      i->first->RemoveAunt(this);
      if (fActive) DoElementUnselect(i);
   }
   fMap.clear();
   if (fActive) SelectionCleared();
   StampObjPropsPreChk();
}

////////////////////////////////////////////////////////////////////////////////
/// Remove element from all implied-selected sets.
///
/// This is called as part of the element destruction from
/// REveManager::PreDeleteElement() and should not be called
/// directly.

void REveSelection::RemoveImpliedSelected(REveElement *el)
{
   bool changed = false;

   for (auto &i : fMap)
   {
      auto j = i.second.f_implied.find(el);
      if (j != i.second.f_implied.end())
      {
         i.second.f_implied.erase(j);
         changed = true;
      }
   }

   if (changed) StampObjPropsPreChk();
}

////////////////////////////////////////////////////////////////////////////////
/// Recalculate implied-selected state for given selection entry.
/// Add new elements to implied-selected set and increase their
/// implied-selected count.

void REveSelection::RecheckImpliedSet(SelMap_i &smi)
{
   bool  changed = false;
   Set_t set;
   smi->first->FillImpliedSelectedSet(set);
   for (auto &i: set)
   {
      if (smi->second.f_implied.find(i) == smi->second.f_implied.end())
      {
         smi->second.f_implied.insert(i);
         i->IncImpliedSelected();
         changed = true;
      }
   }

   if (changed) StampObjPropsPreChk();
}

////////////////////////////////////////////////////////////////////////////////
/// If given element is selected or implied-selected within this
/// selection then recheck implied-set for given selection entry.

void REveSelection::RecheckImpliedSetForElement(REveElement *el)
{
   // Top-level selected.
   {
      auto i = fMap.find(el);
      if (i != fMap.end())
         RecheckImpliedSet(i);
   }

   // Implied selected (we can not tell if by this selection or some other),
   // then we need to loop over all.
   if (el->GetImpliedSelected() > 0)
   {
      for (auto i = fMap.begin(); i != fMap.end(); ++i)
      {
         if (i->second.f_implied.find(el) != i->second.f_implied.end())
            RecheckImpliedSet(i);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Emit SelectionAdded signal.

void REveSelection::SelectionAdded(REveElement* /*el*/)
{
   // XXXX
   // Emit("SelectionAdded(REveElement*)", (Long_t)el);
}

////////////////////////////////////////////////////////////////////////////////
/// Emit SelectionRemoved signal.

void REveSelection::SelectionRemoved(REveElement* /*el*/)
{
   // XXXX
   // Emit("SelectionRemoved(REveElement*)", (Long_t)el);
}

////////////////////////////////////////////////////////////////////////////////
/// Emit SelectionCleared signal.

void REveSelection::SelectionCleared()
{
   // XXXX
   // Emit("SelectionCleared()");
}

////////////////////////////////////////////////////////////////////////////////
/// Emit SelectionRepeated signal.

void REveSelection::SelectionRepeated(REveElement* /*el*/)
{
   // XXXX
   // Emit("SelectionRepeated(REveElement*)", (Long_t)el);
}

////////////////////////////////////////////////////////////////////////////////
/// Activate this selection.

void REveSelection::ActivateSelection()
{
   if (fActive) return;

   fActive = kTRUE;
   for (auto i = fMap.begin(); i != fMap.end(); ++i) {
      DoElementSelect(i);
      SelectionAdded(i->first);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Deactivate this selection.

void REveSelection::DeactivateSelection()
{
   if (!fActive) return;

   for (auto i = fMap.begin(); i != fMap.end(); ++i) {
      DoElementUnselect(i);
   }
   SelectionCleared();
   fActive = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Given element el that was picked or clicked by the user, find
/// the parent/ancestor element that should actually become the main
/// selected element according to current selection mode.

REveElement* REveSelection::MapPickedToSelected(REveElement* el)
{
   if (el == nullptr)
      return nullptr;

   for (int pick_to_select : fPickToSelect)
   {
      switch (pick_to_select)
      {
         case kPS_Ignore:
         {
            return nullptr;
         }
         case kPS_Element:
         {
            return el;
         }
         case kPS_Projectable:
         {
            REveProjected* pted = dynamic_cast<REveProjected*>(el);
            if (pted)
               return dynamic_cast<REveElement*>(pted->GetProjectable());
            break;
         }
         case kPS_Compound:
         {
            REveElement* cmpnd = el->GetCompound();
            if (cmpnd)
               return cmpnd;
            break;
         }
         case kPS_PableCompound:
         {
            REveProjected* pted = dynamic_cast<REveProjected*>(el);
            if (pted)
               el = dynamic_cast<REveElement*>(pted->GetProjectable());
            REveElement* cmpnd = el->GetCompound();
            if (cmpnd)
               return cmpnd;
            if (pted)
               return el;
            break;
         }
         case kPS_Master:
         {
            REveElement* mstr = el->GetSelectionMaster();
            if (mstr)
               return mstr;
            break;
         }
      }
   }

   return el;
}

////////////////////////////////////////////////////////////////////////////////
/// Called when user picks/clicks on an element. If multi is true,
/// the user is requiring a multiple selection (usually this is
/// associated with control-key being pressed at the time of pick
/// event).
/// XXXX Old interface, not used in EVE-7.

void REveSelection::UserPickedElement(REveElement* el, Bool_t multi)
{
   el = MapPickedToSelected(el);

   if (el || NotEmpty())
   {
      if ( ! multi)
         RemoveNieces();
      if (el)
      {
         if (HasNiece(el))
             RemoveNiece(el);
         else
            AddNiece(el);
      }
      StampObjProps();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Called when element selection is repeated.
/// XXXX Old interface, not used in EVE-7.

void REveSelection::UserRePickedElement(REveElement* el)
{
   el = MapPickedToSelected(el);
   if (el && HasNiece(el))
   {
      SelectionRepeated(el);
      StampObjProps();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Called when an element is unselected.
/// XXXX Old interface, not used in EVE-7.

void REveSelection::UserUnPickedElement(REveElement* el)
{
   el = MapPickedToSelected(el);
   if (el && HasNiece(el))
   {
      RemoveNiece(el);
      StampObjProps();
   }
}

//==============================================================================

void REveSelection::NewElementPicked(ElementId_t id, bool multi, bool secondary, const std::set<int>& in_secondary_idcs)
{
   static const REveException eh("REveSelection::NewElementPicked ");

   REveElement *pel = nullptr, *el = nullptr;

   // AMT the forth/last argument is optional and therefore need to be constant
   std::set<int> secondary_idcs = in_secondary_idcs;

   if (id > 0)
   {
     pel = REX::gEve->FindElementById(id);

     if ( ! pel) throw eh + "picked element id=" + id + " not found.";

     el = MapPickedToSelected(pel);

     if (el != pel) {
        REveSecondarySelectable* ss = dynamic_cast<REveSecondarySelectable*>(el);
        if (!secondary && ss) {
           secondary = true;
           secondary_idcs = ss->RefSelectedSet();
        }
     }
   }

   if (gDebug > 0) {
      std::string debug_secondary;
      if (secondary) {
         debug_secondary = " {";
         for (auto si : secondary_idcs) {
            debug_secondary.append(" ");
            debug_secondary.append(std::to_string(si));
         }
         debug_secondary.append(" }");
      }
      ::Info("REveSelection::NewElementPicked", "%p -> %p, multi: %d, secondary: %d  %s", pel, el, multi, secondary, debug_secondary.c_str());
   }

   Record *rec = find_record(el);

   bool changed = true;

   if (multi)
   {
      if (el)
      {
         if (rec)
         {
            assert(secondary == rec->is_secondary());
            if (secondary || rec->is_secondary())
            {
               std::set<int> dup;
               for (auto &ns :  secondary_idcs)
               {
                  int nsi = ns;
                  auto ir = rec->f_sec_idcs.insert(nsi);
                  if (!ir.second)
                     dup.insert(nsi);
               }

               // erase duplicates
               for (auto &dit :  dup)
                  rec->f_sec_idcs.erase(dit);

               secondary_idcs  = rec->f_sec_idcs;
               if (!secondary_idcs.empty()) {
                  AddNiece(el);
                  rec = find_record(el);
                  rec->f_is_sec   = true;
                  rec->f_sec_idcs = secondary_idcs;
               }
            }
            else
            {
               RemoveNiece(el);
            }
         }
         else
         {
            AddNiece(el);
            rec = find_record(el);
            rec->f_is_sec   = true;
            rec->f_sec_idcs = secondary_idcs;
         }
      }
      else
      {
         // Multiple selection with 0 element ... do nothing, I think.
         changed = false;
      }
   }
   else // single selection (not multi)
   {
      if (el)
      {
         if (rec)
         {
            if (secondary)
            {
                bool modified = (rec->f_sec_idcs != secondary_idcs);
               RemoveNieces();
               // re-adding is needed to refresh implied selected
               if (modified) {
                  AddNiece(el);
                  rec = find_record(el);
                  rec->f_is_sec   = true;
                  rec->f_sec_idcs = secondary_idcs;
               }
            }
            else
            {
               RemoveNiece(el);
            }
         }
         else
         {
            if (HasNieces()) RemoveNieces();
            AddNiece(el);
            if (secondary)
            {
               rec = find_record(el);
               rec->f_is_sec   = true;
               rec->f_sec_idcs = secondary_idcs;
            }
         }
      }
      else // Single selection with zero element --> clear selection.
      {
         if (HasNieces())
           RemoveNieces();
         else
           changed = false;
      }
   }

   if (changed)
     StampObjProps();
}

////////////////////////////////////////////////////////////////////////////////
/// Clear selection if not empty.

void REveSelection::ClearSelection()
{
   if (HasNieces())
   {
      RemoveNieces();
      StampObjProps();
   }
}

//==============================================================================

////////////////////////////////////////////////////////////////////////////////
/// Remove pointers to el from implied selected sets.

int REveSelection::RemoveImpliedSelectedReferencesTo(REveElement *el)
{
   int count = 0;

   for (auto &i : fMap)
   {
      auto j = i.second.f_implied.find(el);

      if (j != i.second.f_implied.end())
      {
         i.second.f_implied.erase(j);
         el->DecImpliedSelected();
         ++count;
      }
   }

   return count;
}

////////////////////////////////////////////////////////////////////////////////
/// Write core json. If rnr_offset negative, render data will not be written

Int_t REveSelection::WriteCoreJson(nlohmann::json &j, Int_t /* rnr_offset */)
{
   REveElement::WriteCoreJson(j, -1);

   j["fVisibleEdgeColor"] = fVisibleEdgeColor;
   j["fHiddenEdgeColor"]  = fHiddenEdgeColor;

   nlohmann::json sel_list = nlohmann::json::array();

   for (auto &i : fMap)
   {
      nlohmann::json rec = {}, imp = nlohmann::json::array(), sec = nlohmann::json::array();

      rec["primary"] = i.first->GetElementId();

      // XXX if not empty / f_is_sec is false ???
      for (auto &sec_id : i.second.f_sec_idcs)
         sec.push_back(sec_id);

      // XXX if not empty ???
      for (auto &imp_el : i.second.f_implied) {
         imp.push_back(imp_el->GetElementId());
         imp_el->FillExtraSelectionData(rec["extra"], sec);

      }
      rec["implied"]  = imp;


      if (i.first->RequiresExtraSelectionData()) {
         i.first->FillExtraSelectionData(rec["extra"], sec);
      }

      rec["sec_idcs"] = sec;

      // stream tooltip in highlight type
      if (!fIsMaster)
         rec["tooltip"] = i.first->GetHighlightTooltip(i.second.f_sec_idcs);

      sel_list.push_back(rec);
   }

   j["sel_list"] = sel_list;

   j["UT_PostStream"] = "UT_Selection_Refresh_State"; // XXXX to be canonized

   // std::cout << j.dump(4) << std::endl;

   return 0;
}
