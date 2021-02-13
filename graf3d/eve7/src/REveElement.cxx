// @(#)root/eve7:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007, 2018

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/REveElement.hxx>
#include <ROOT/REveUtil.hxx>
#include <ROOT/REveScene.hxx>
#include <ROOT/REveCompound.hxx>
#include <ROOT/REveTrans.hxx>
#include <ROOT/REveManager.hxx>
#include <ROOT/REveSelection.hxx>
#include <ROOT/REveProjectionBases.hxx>
#include <ROOT/REveProjectionManager.hxx>
#include <ROOT/REveRenderData.hxx>

#include "TGeoMatrix.h"

#include "TClass.h"
#include "TPRegexp.h"
#include "TROOT.h"
#include "TColor.h"

#include <cassert>


#include <algorithm>

using namespace ROOT::Experimental;
namespace REX = ROOT::Experimental;

/** \class REveElement
\ingroup REve
Base class for REveUtil visualization elements, providing hierarchy
management, rendering control and list-tree item management.

Class of acceptable children can be limited by setting the
fChildClass member.
*/

////////////////////////////////////////////////////////////////////////////////
/// Default constructor.

REveElement::REveElement(const std::string& name, const std::string& title) :
   fName                (name),
   fTitle               (title)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor. Does shallow copy.
/// For deep-cloning and children-cloning, see:
/// ~~~ {.cpp}
///   REveElement* CloneElementRecurse(Int_t level)
///   void         CloneChildrenRecurse(REveElement* dest, Int_t level)
/// ~~~
/// 'void* UserData' is NOT copied.
/// If the element is projectable, its projections are NOT copied.
///
/// Not implemented for most sub-classes, let us know.
/// Note that sub-classes of REveProjected are NOT and will NOT be copyable.

REveElement::REveElement(const REveElement& e) :
   fName                (e.fName),
   fTitle               (e.fTitle),
   fChildClass          (e.fChildClass),
   fVizTag              (e.fVizTag),
   fDestroyOnZeroRefCnt (e.fDestroyOnZeroRefCnt),
   fRnrSelf             (e.fRnrSelf),
   fRnrChildren         (e.fRnrChildren),
   fCanEditMainColor    (e.fCanEditMainColor),
   fCanEditMainTransparency(e.fCanEditMainTransparency),
   fCanEditMainTrans    (e.fCanEditMainTrans),
   fMainTransparency    (e.fMainTransparency),
   fPickable            (e.fPickable),
   fCSCBits             (e.fCSCBits)
{
   SetVizModel(e.fVizModel);
   // FIXME: from Sergey: one have to use other way to referencing main color
   if (e.fMainColorPtr)
      fMainColorPtr = (Color_t*)((const char*) this + ((const char*) e.fMainColorPtr - (const char*) &e));
   if (e.fMainTrans)
      fMainTrans = std::make_unique<REveTrans>(*e.fMainTrans.get());
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor. Do not call this method directly, either call Destroy() or
/// Annihilate(). See also DestroyElements() and AnnihilateElements() if you
/// need to delete all children of an element.

REveElement::~REveElement()
{
   if (fDestructing != kAnnihilate)
   {
      fDestructing = kStandard;
      RemoveElementsInternal();

      if (fMother) {
        fMother->RemoveElementLocal(this);
        fMother->fChildren.remove(this);
      }

      if (fScene) {
         fScene->SceneElementRemoved( fElementId);
      }

      for (auto &au : fAunts)
      {
         au->RemoveNieceInternal(this);
      }
   }
}

ElementId_t REveElement::get_mother_id() const
{
   return fMother ? fMother->GetElementId() : 0;
}

ElementId_t REveElement::get_scene_id() const
{
   return fScene ? fScene->GetElementId() : 0;
}

void REveElement::assign_element_id_recurisvely()
{
   assert(fElementId == 0);

   REX::gEve->AssignElementId(this);
   for (auto &c : fChildren)
      c->assign_element_id_recurisvely();
}

std::string REveElement::GetHighlightTooltip(const std::set<int>& iSet) const
{
   std::string res = fTitle;
   if (res.empty())
      res = fName;

   if (!iSet.empty())
      res = TString::Format("%s idx=%d", res.c_str(), *iSet.begin());

   return res;
}

void REveElement::assign_scene_recursively(REveScene* s)
{
   assert(fScene == nullptr);

   fScene = s;

   // XXX MT -- Why do we have fDestructing here? Can this really happen?
   //           If yes, shouldn't we block it in AddElement() already?
   if (fDestructing == kNone && fScene && fScene->IsAcceptingChanges())
   {
      StampElementAdded();
   }
   for (auto &c : fChildren)
      c->assign_scene_recursively(s);
}

////////////////////////////////////////////////////////////////////////////////
/// Called before the element is deleted, thus offering the last chance
/// to detach from acquired resources and from the framework itself.
/// Here the request is just passed to REveManager.
/// If you override it, make sure to call base-class version.

void REveElement::PreDeleteElement()
{
   if (fElementId != 0)
   {
      REX::gEve->PreDeleteElement(this);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Clone the element via copy constructor.
/// Should be implemented for all classes that require cloning support.

REveElement* REveElement::CloneElement() const
{
   return new REveElement(*this);
}

////////////////////////////////////////////////////////////////////////////////
/// Clone elements and recurse 'level' deep over children.
///  - If level ==  0, only the element itself is cloned (default).
///  - If level == -1, all the hierarchy is cloned.

REveElement* REveElement::CloneElementRecurse(Int_t level) const
{
   REveElement* el = CloneElement();
   if (level--)
   {
      CloneChildrenRecurse(el, level);
   }
   return el;
}

////////////////////////////////////////////////////////////////////////////////
/// Clone children and attach them to the dest element.
/// If level ==  0, only the direct descendants are cloned (default).
/// If level == -1, all the hierarchy is cloned.

void REveElement::CloneChildrenRecurse(REveElement* dest, Int_t level) const
{
   for (auto &c: fChildren)
      dest->AddElement(c->CloneElementRecurse(level));
}


//==============================================================================

////////////////////////////////////////////////////////////////////////////////
/// Set name of an element.

void REveElement::SetName(const std::string& name)
{
   fName = name;
   NameTitleChanged();
}

////////////////////////////////////////////////////////////////////////////////
/// Set title of an element.

void REveElement::SetTitle(const std::string& title)
{
   fTitle = title;
   NameTitleChanged();
}

////////////////////////////////////////////////////////////////////////////////
/// Set name and title of an element.
/// Here we attempt to cast the assigned object into TNamed and call
/// SetNameTitle() there.
/// If you override this call NameTitleChanged() from there.

void REveElement::SetNameTitle(const std::string& name, const std::string& title)
{
   fName  = name;
   fTitle = title;
   NameTitleChanged();
}

////////////////////////////////////////////////////////////////////////////////
/// Virtual function called when a name or title of the element has
/// been changed.
/// If you override this, call also the version of your direct base-class.

void REveElement::NameTitleChanged()
{
   // Should send out some message. Need a new stamp type?
}

////////////////////////////////////////////////////////////////////////////////
/// Set visualization-parameter model element.
/// Calling of this function from outside of EVE should in principle
/// be avoided as it can lead to dis-synchronization of viz-tag and
/// viz-model.

void REveElement::SetVizModel(REveElement* model)
{
   fVizModel = model;
}

////////////////////////////////////////////////////////////////////////////////
/// Find model element in VizDB that corresponds to previously
/// assigned fVizTag and set fVizModel accordingly.
/// If the tag is not found in VizDB, the old model-element is kept
/// and false is returned.

Bool_t REveElement::SetVizModelByTag()
{
   REveElement* model = REX::gEve->FindVizDBEntry(fVizTag);
   if (model)
   {
      SetVizModel(model);
      return kTRUE;
   }
   else
   {
      return kFALSE;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set the VizTag, find model-element from the VizDB and copy
/// visualization-parameters from it. If the model is not found and
/// fallback_tag is non-null, search for it is attempted as well.
/// For example: ApplyVizTag("TPC Clusters", "Clusters");
///
/// If the model-element can not be found a warning is printed and
/// false is returned.

Bool_t REveElement::ApplyVizTag(const TString& tag, const TString& fallback_tag)
{
   REveElement* model;

   if ((model = REX::gEve->FindVizDBEntry(tag)) != nullptr)
   {
      SetVizTag(tag);
   }
   else if ( ! fallback_tag.IsNull() && (model = REX::gEve->FindVizDBEntry(fallback_tag)) != nullptr)
   {
      SetVizTag(fallback_tag);
   }

   if (model)
   {
      SetVizModel(model);
      CopyVizParamsFromDB();
      return true;
   }
   Warning("REveElement::ApplyVizTag", "entry for tag '%s' not found in VizDB.", tag.Data());
   return false;
}

////////////////////////////////////////////////////////////////////////////////
/// Propagate visualization parameters to dependent elements.
///
/// MainColor is propagated independently in SetMainColor().
/// In this case, as fMainColor is a pointer to Color_t, it should
/// be set in TProperClass::CopyVizParams().
///
/// Render state is not propagated. Maybe it should be, at least optionally.

void REveElement::PropagateVizParamsToProjecteds()
{
   REveProjectable* pable = dynamic_cast<REveProjectable*>(this);
   if (pable && pable->HasProjecteds())
   {
      pable->PropagateVizParams();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Propagate visualization parameters from element el (defaulting
/// to this) to all children.
///
/// The primary use of this is for model-elements from
/// visualization-parameter database.

void REveElement::PropagateVizParamsToChildren(REveElement* el)
{
   if (!el) el = this;

   for (auto &c : fChildren)
   {
      c->CopyVizParams(el);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Copy visualization parameters from element el.
/// This method needs to be overriden by any class that introduces
/// new parameters.

void REveElement::CopyVizParams(const REveElement* el)
{
   fCanEditMainColor        = el->fCanEditMainColor;
   fCanEditMainTransparency = el->fCanEditMainTransparency;
   fMainTransparency        = el->fMainTransparency;
   if (fMainColorPtr == & fDefaultColor) fDefaultColor = el->GetMainColor();

   AddStamp(kCBColorSelection | kCBObjProps);
}

////////////////////////////////////////////////////////////////////////////////
/// Copy visualization parameters from the model-element fVizModel.
/// A warning is printed if the model-element fVizModel is not set.

void REveElement::CopyVizParamsFromDB()
{
   if (fVizModel)
   {
      CopyVizParams(fVizModel);
   }
   else
   {
      Warning("REveElement::CopyVizParamsFromDB", "VizModel has not been set.");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Save visualization parameters for this element with given tag.
///
/// This function creates the instantiation code, calls virtual
/// WriteVizParams() and, at the end, writes out the code for
/// registration of the model into the VizDB.

void REveElement::SaveVizParams(std::ostream& out, const TString& tag, const TString& var)
{
   static const REveException eh("REveElement::GetObject ");

   TString t = "   ";
   TString cls(IsA()->GetName());

   out << "\n";

   TString intro = " TAG='" + tag + "', CLASS='" + cls + "'";
   out << "   //" << intro << "\n";
   out << "   //" << TString('-', intro.Length()) << "\n";
   out << t << cls << "* " << var <<" = new " << cls << ";\n";

   WriteVizParams(out, var);

   out << t << "REX::gEve->InsertVizDBEntry(\"" << tag << "\", "<< var <<");\n";
}

////////////////////////////////////////////////////////////////////////////////
/// Write-out visual parameters for this object.
/// This is a virtual function and all sub-classes are required to
/// first call the base-element version.
/// The name of the element pointer is 'x%03d', due to cint limitations.
/// Three spaces should be used for indentation, same as in
/// SavePrimitive() methods.

void REveElement::WriteVizParams(std::ostream& out, const TString& var)
{
   TString t = "   " + var + "->";

   out << t << "SetElementName(\""  << fName  << "\");\n";
   out << t << "SetElementTitle(\"" << fTitle << "\");\n";
   out << t << "SetEditMainColor("  << fCanEditMainColor << ");\n";
   out << t << "SetEditMainTransparency(" << fCanEditMainTransparency << ");\n";
   out << t << "SetMainTransparency("     << fMainTransparency << ");\n";
}

////////////////////////////////////////////////////////////////////////////////
/// Set visual parameters for this object for given tag.

void REveElement::VizDB_Apply(const std::string& tag)
{
   if (ApplyVizTag(tag))
   {
      PropagateVizParamsToProjecteds();
      REX::gEve->Redraw3D();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Reset visual parameters for this object from VizDB.
/// The model object must be already set.

void REveElement::VizDB_Reapply()
{
   if (fVizModel)
   {
      CopyVizParamsFromDB();
      PropagateVizParamsToProjecteds();
      REX::gEve->Redraw3D();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Copy visual parameters from this element to viz-db model.
/// If update is set, all clients of the model will be updated to
/// the new value.
/// A warning is printed if the model-element fVizModel is not set.

void REveElement::VizDB_UpdateModel(Bool_t update)
{
   if (fVizModel)
   {
      fVizModel->CopyVizParams(this);
      if (update)
      {
         // XXX Back references from vizdb templates have been removed in Eve7.
         // XXX We could traverse all scenes and elementes and reset those that
         // XXX have a matching fVizModel. Or something.
         Error("VizDB_UpdateModel", "update from vizdb -> elements not implemented.");
         // fVizModel->PropagateVizParamsToElements(fVizModel);
         // REX::gEve->Redraw3D();
      }
   }
   else
   {
      Warning("VizDB_UpdateModel", "VizModel has not been set.");
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Create a replica of element and insert it into VizDB with given tag.
/// If replace is true an existing element with the same tag will be replaced.
/// If update is true, existing client of tag will be updated.

void REveElement::VizDB_Insert(const std::string& tag, Bool_t replace, Bool_t update)
{
   static const REveException eh("REveElement::GetObject ");

   TClass* cls     = IsA();
   REveElement* el = reinterpret_cast<REveElement*>(cls->New());
   if (!el) {
      Error("VizDB_Insert", "Creation of replica failed.");
      return;
   }
   el->CopyVizParams(this);
   Bool_t succ = REX::gEve->InsertVizDBEntry(tag, el, replace, update);
   if (succ && update)
      REX::gEve->Redraw3D();
}

////////////////////////////////////////////////////////////////////////////////
/// Add el into the list aunts.
///
/// Adding aunt is subordinate to adding a niece.
/// This is an internal function.

void REveElement::AddAunt(REveAunt *au)
{
   assert(au != nullptr);

   fAunts.emplace_back(au);
}

////////////////////////////////////////////////////////////////////////////////
/// Remove el from the list of aunts.
/// Removing aunt is subordinate to removing a niece.
/// This is an internal function.

void REveElement::RemoveAunt(REveAunt *au)
{
   assert(au != nullptr);

   fAunts.remove(au);
}

/******************************************************************************/

////////////////////////////////////////////////////////////////////////////////
/// Check external references to this and eventually auto-destruct
/// the render-element.

void REveElement::CheckReferenceCount(const std::string& from)
{
   if (fDestructing != kNone)
      return;

   if (fMother == nullptr && fDestroyOnZeroRefCnt && fDenyDestroy <= 0)
   {
      if (gDebug > 0)
         Info("REveElement::CheckReferenceCount", "(called from %s) auto-destructing '%s' on zero reference count.",
              from.c_str(), GetCName());

      PreDeleteElement();
      delete this;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return class for this element

TClass *REveElement::IsA() const
{
   return TClass::GetClass(typeid(*this), kTRUE, kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// Export render-element to CINT with variable name var_name.

void REveElement::ExportToCINT(const char *var_name)
{
   const char* cname = IsA()->GetName();
   gROOT->ProcessLine(TString::Format("%s* %s = (%s*)0x%lx;", cname, var_name, cname, (ULong_t)this));
}

////////////////////////////////////////////////////////////////////////////////
/// Set render state of this element, i.e. if it will be published
/// on next scene update pass.
/// Returns true if the state has changed.

Bool_t REveElement::SetRnrSelf(Bool_t rnr)
{
   if (SingleRnrState())
   {
      return SetRnrState(rnr);
   }

   if (rnr != fRnrSelf)
   {
      fRnrSelf = rnr;
      StampVisibility();
      PropagateRnrStateToProjecteds();

      return kTRUE;
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set render state of this element's children, i.e. if they will
/// be published on next scene update pass.
/// Returns true if the state has changed.

Bool_t REveElement::SetRnrChildren(Bool_t rnr)
{
   if (SingleRnrState())
   {
      return SetRnrState(rnr);
   }

   if (rnr != fRnrChildren)
   {
      fRnrChildren = rnr;
      StampVisibility();
      PropagateRnrStateToProjecteds();
      return kTRUE;
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set state for rendering of this element and its children.
/// Returns true if the state has changed.

Bool_t REveElement::SetRnrSelfChildren(Bool_t rnr_self, Bool_t rnr_children)
{
   if (SingleRnrState())
   {
      return SetRnrState(rnr_self);
   }

   if (fRnrSelf != rnr_self || fRnrChildren != rnr_children)
   {
      fRnrSelf     = rnr_self;
      fRnrChildren = rnr_children;
      StampVisibility();
      PropagateRnrStateToProjecteds();
      return kTRUE;
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set render state of this element and of its children to the same
/// value.
/// Returns true if the state has changed.

Bool_t REveElement::SetRnrState(Bool_t rnr)
{
   if (fRnrSelf != rnr || fRnrChildren != rnr)
   {
      fRnrSelf = fRnrChildren = rnr;
      StampVisibility();
      PropagateRnrStateToProjecteds();
      return kTRUE;
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Propagate render state to the projected replicas of this element.
/// Maybe this should be optional on REX::gEve/element level.

void REveElement::PropagateRnrStateToProjecteds()
{
   REveProjectable *pable = dynamic_cast<REveProjectable*>(this);
   if (pable && pable->HasProjecteds())
   {
      pable->PropagateRenderState(fRnrSelf, fRnrChildren);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set up element to use built-in main color and set flags allowing editing
/// of main color and transparency.

void REveElement::SetupDefaultColorAndTransparency(Color_t col, Bool_t can_edit_color, Bool_t can_edit_transparency)
{
   fMainColorPtr = & fDefaultColor;
   fDefaultColor = col;
   fCanEditMainColor = can_edit_color;
   fCanEditMainTransparency = can_edit_transparency;
}

////////////////////////////////////////////////////////////////////////////////
/// Set main color of the element.

void REveElement::SetMainColor(Color_t color)
{
   Color_t old_color = GetMainColor();

   if (fMainColorPtr)
   {
      *fMainColorPtr = color;
      StampColorSelection();
   }

   PropagateMainColorToProjecteds(color, old_color);
}

////////////////////////////////////////////////////////////////////////////////
/// Convert pixel to Color_t and call SetMainColor().

void REveElement::SetMainColorPixel(Pixel_t pixel)
{
   SetMainColor(TColor::GetColor(pixel));
}

////////////////////////////////////////////////////////////////////////////////
/// Convert RGB values to Color_t and call SetMainColor.

void REveElement::SetMainColorRGB(UChar_t r, UChar_t g, UChar_t b)
{
   SetMainColor(TColor::GetColor(r, g, b));
}

////////////////////////////////////////////////////////////////////////////////
/// Convert RGB values to Color_t and call SetMainColor.

void REveElement::SetMainColorRGB(Float_t r, Float_t g, Float_t b)
{
   SetMainColor(TColor::GetColor(r, g, b));
}

////////////////////////////////////////////////////////////////////////////////
/// Propagate color to projected elements.

void REveElement::PropagateMainColorToProjecteds(Color_t color, Color_t old_color)
{
   REveProjectable* pable = dynamic_cast<REveProjectable*>(this);
   if (pable && pable->HasProjecteds())
   {
      pable->PropagateMainColor(color, old_color);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set main-transparency.
/// Transparency is clamped to [0, 100].

void REveElement::SetMainTransparency(Char_t t)
{
   Char_t old_t = GetMainTransparency();

   if (t > 100) t = 100;
   fMainTransparency = t;
   StampColorSelection();

   PropagateMainTransparencyToProjecteds(t, old_t);
}

////////////////////////////////////////////////////////////////////////////////
/// Set main-transparency via float alpha variable.
/// Value of alpha is clamped t0 [0, 1].

void REveElement::SetMainAlpha(Float_t alpha)
{
   if (alpha < 0) alpha = 0;
   if (alpha > 1) alpha = 1;
   SetMainTransparency((Char_t) (100.0f*(1.0f - alpha)));
}

////////////////////////////////////////////////////////////////////////////////
/// Propagate transparency to projected elements.

void REveElement::PropagateMainTransparencyToProjecteds(Char_t t, Char_t old_t)
{
   REveProjectable* pable = dynamic_cast<REveProjectable*>(this);
   if (pable && pable->HasProjecteds())
   {
      pable->PropagateMainTransparency(t, old_t);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to main transformation. If 'create' flag is set (default)
/// it is created if not yet existing.

REveTrans *REveElement::PtrMainTrans(Bool_t create)
{
   if (!fMainTrans && create)
      InitMainTrans();

   return fMainTrans.get();
}

////////////////////////////////////////////////////////////////////////////////
/// Return reference to main transformation. It is created if not yet
/// existing.

REveTrans &REveElement::RefMainTrans()
{
   if (!fMainTrans)
      InitMainTrans();

   return *fMainTrans.get();
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize the main transformation to identity matrix.
/// If can_edit is true (default), the user will be able to edit the
/// transformation parameters via GUI.

void REveElement::InitMainTrans(Bool_t can_edit)
{
   if (fMainTrans)
      fMainTrans->UnitTrans();
   else
      fMainTrans = std::make_unique<REveTrans>();
   fCanEditMainTrans = can_edit;
}

////////////////////////////////////////////////////////////////////////////////
/// Destroy the main transformation matrix, it will always be taken
/// as identity. Editing of transformation parameters is disabled.

void REveElement::DestroyMainTrans()
{
   fMainTrans.reset(nullptr);
   fCanEditMainTrans = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Set transformation matrix from column-major array.

void REveElement::SetTransMatrix(Double_t* carr)
{
   RefMainTrans().SetFrom(carr);
}

////////////////////////////////////////////////////////////////////////////////
/// Set transformation matrix from TGeo's matrix.

void REveElement::SetTransMatrix(const TGeoMatrix& mat)
{
   RefMainTrans().SetFrom(mat);
}

////////////////////////////////////////////////////////////////////////////////
/// Check if el can be added to this element.
///
/// Here we make sure the new child is not equal to this and, if fChildClass
/// is set, that it is inherited from it.

Bool_t REveElement::AcceptElement(REveElement* el)
{
   if (el == this)
      return kFALSE;
   if (fChildClass && ! el->IsA()->InheritsFrom(fChildClass))
      return kFALSE;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Add el to the list of children.

void REveElement::AddElement(REveElement *el)
{
   static const REveException eh("REveElement::AddElement ");

   if (!el)                  throw eh + "called with nullptr argument.";
   if ( ! AcceptElement(el)) throw eh + Form("parent '%s' rejects '%s'.", GetCName(), el->GetCName());
   if (el->fElementId)       throw eh + "element already has an id.";
   // if (el->fScene)           throw eh + "element already has a Scene.";
   if (el->fMother)          throw eh + "element already has a Mother.";

   // XXX Implement reparent --> MoveElement() ????
   //     Actually, better to do new = old.Clone(), RemoveElement(old), AddElement(new);
   //     Or do magick with Inc/DecDenyDestroy().
   //     PITA with existing children !!!! Need to re-scene them.

   if (fElementId)             el->assign_element_id_recurisvely();
   if (fScene && ! el->fScene) el->assign_scene_recursively(fScene);

   el->fMother = this;

   fChildren.emplace_back(el);

   // XXXX This should be element added. Also, should be different for
   // "full (re)construction". Scenes should manage that and have
   // state like: none - constructing - clearing - nominal - updating.
   // I recon this means an element should have a ptr to its scene.
   //
   // ElementChanged();
}

////////////////////////////////////////////////////////////////////////////////
/// Remove el from the list of children.

void REveElement::RemoveElement(REveElement* el)
{
   static const REveException eh("REveElement::RemoveElement ");

   if (!el)                 throw eh + "called with nullptr argument.";
   if (el->fMother != this) throw eh + "this element is not mother of el.";

   RemoveElementLocal(el);

   el->fScene->SceneElementRemoved(fElementId);
   el->fMother = nullptr;
   el->fScene  = nullptr;

   el->CheckReferenceCount();

   fChildren.remove(el);

   // XXXX This should be ElementRemoved(). Also, think about recursion, deletion etc.
   // Also, this seems to be done above, in the call to fScene.
   //
   // ElementChanged();
}

////////////////////////////////////////////////////////////////////////////////
/// Perform additional local removal of el.
/// Called from RemoveElement() which does whole untangling.
/// Put into special function as framework-related handling of
/// element removal should really be common to all classes and
/// clearing of local structures happens in between removal
/// of list-tree-items and final removal.
/// If you override this, you should also override
/// RemoveElementsLocal().

void REveElement::RemoveElementLocal(REveElement* /*el*/)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Remove all elements. This assumes removing of all elements can
/// be done more efficiently then looping over them and removing one
/// by one. This protected function performs the removal on the
/// level of REveElement.

void REveElement::RemoveElementsInternal()
{
   RemoveElementsLocal();

   for (auto &c : fChildren)
   {
      c->fScene->SceneElementRemoved(c->fElementId);
      c->fMother = nullptr;
      c->fScene  = nullptr;

      c->CheckReferenceCount();
   }

   fChildren.clear();
}

////////////////////////////////////////////////////////////////////////////////
/// Remove all elements. This assumes removing of all elements can
/// be done more efficiently then looping over them and removing
/// them one by one.

void REveElement::RemoveElements()
{
   if (HasChildren())
   {
      RemoveElementsInternal();
      // ElementChanged();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Perform additional local removal of all elements.
/// See comment to RemoveElementLocal(REveElement*).

void REveElement::RemoveElementsLocal()
{
}

////////////////////////////////////////////////////////////////////////////////
/// If this is a projectable, loop over all projected replicas and
/// add the projected image of child 'el' there. This is supposed to
/// be called after you add a child to a projectable after it has
/// already been projected.
/// You might also want to call RecheckImpliedSelections() on this
/// element or 'el'.
///
/// If 'same_depth' flag is true, the same depth as for parent object
/// is used in every projection. Otherwise current depth of each
/// relevant projection-manager is used.

void REveElement::ProjectChild(REveElement* el, Bool_t same_depth)
{
   REveProjectable* pable = dynamic_cast<REveProjectable*>(this);
   if (pable && HasChild(el))
   {
      for (auto &pp: pable->RefProjecteds())
      {
         auto pmgr = pp->GetManager();
         Float_t cd = pmgr->GetCurrentDepth();
         if (same_depth) pmgr->SetCurrentDepth(pp->GetDepth());

         pmgr->SubImportElements(el, pp->GetProjectedAsElement());

         if (same_depth) pmgr->SetCurrentDepth(cd);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// If this is a projectable, loop over all projected replicas and
/// add the projected image of all children there. This is supposed
/// to be called after you destroy all children and then add new
/// ones after this element has already been projected.
/// You might also want to call RecheckImpliedSelections() on this
/// element.
///
/// If 'same_depth' flag is true, the same depth as for the
/// projected element is used in every projection. Otherwise current
/// depth of each relevant projection-manager is used.

void REveElement::ProjectAllChildren(Bool_t same_depth)
{
   REveProjectable* pable = dynamic_cast<REveProjectable*>(this);
   if (pable)
   {
      for (auto &pp: pable->RefProjecteds())
      {
         REveProjectionManager *pmgr = pp->GetManager();
         Float_t cd = pmgr->GetCurrentDepth();
         if (same_depth) pmgr->SetCurrentDepth(pp->GetDepth());

         pmgr->SubImportChildren(this, pp->GetProjectedAsElement());

         if (same_depth) pmgr->SetCurrentDepth(cd);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Check if element el is a child of this element.

Bool_t REveElement::HasChild(REveElement* el)
{
   return (std::find(fChildren.begin(), fChildren.end(), el) != fChildren.end());
}

////////////////////////////////////////////////////////////////////////////////
/// Find the first child with given name.  If cls is specified (non
/// 0), it is also checked.
///
/// Returns nullptr if not found.

REveElement *REveElement::FindChild(const TString&  name, const TClass* cls)
{
   for (auto &c: fChildren)
   {
      if (name.CompareTo(c->GetCName()) == 0)
      {
         if (!cls || c->IsA()->InheritsFrom(cls))
            return c;
      }
   }
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Find the first child whose name matches regexp. If cls is
/// specified (non 0), it is also checked.
///
/// Returns nullptr if not found.

REveElement* REveElement::FindChild(TPRegexp &regexp, const TClass *cls)
{
   for (auto &c: fChildren)
   {
      if (regexp.MatchB(c->GetName()))
      {
         if (!cls || c->IsA()->InheritsFrom(cls))
            return c;
      }
   }
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Find all children with given name and append them to matches
/// list. If class is specified (non 0), it is also checked.
///
/// Returns number of elements added to the list.

Int_t REveElement::FindChildren(List_t& matches,
                                const TString& name, const TClass* cls)
{
   Int_t count = 0;
   for (auto &c: fChildren)
   {
      if (name.CompareTo(c->GetCName()) == 0)
      {
         if (!cls || c->IsA()->InheritsFrom(cls))
         {
            matches.push_back(c);
            ++count;
         }
      }
   }
   return count;
}

////////////////////////////////////////////////////////////////////////////////
/// Find all children whose name matches regexp and append them to
/// matches list.
///
/// Returns number of elements added to the list.

Int_t REveElement::FindChildren(List_t &matches,
                                TPRegexp &regexp, const TClass *cls)
{
   Int_t count = 0;
   for (auto &c : fChildren)
   {
      if (regexp.MatchB(c->GetCName()))
      {
         if (!cls || c->IsA()->InheritsFrom(cls))
         {
            matches.push_back(c);
            ++count;
         }
      }
   }
   return count;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the first child element or 0 if the list is empty.

REveElement* REveElement::FirstChild() const
{
   return HasChildren() ? fChildren.front() : nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the last child element or 0 if the list is empty.

REveElement* REveElement::LastChild () const
{
   return HasChildren() ? fChildren.back() : nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Enable rendering of children and their list contents.
/// Arguments control how to set self/child rendering.

void REveElement::EnableListElements(Bool_t rnr_self,  Bool_t rnr_children)
{
   for (auto &c: fChildren)
      c->SetRnrSelfChildren(rnr_self, rnr_children);
}

////////////////////////////////////////////////////////////////////////////////
/// Disable rendering of children and their list contents.
/// Arguments control how to set self/child rendering.
///
/// Same as above function, but default arguments are different. This
/// is convenient for calls via context menu.

void REveElement::DisableListElements(Bool_t rnr_self,  Bool_t rnr_children)
{
   for (auto &c: fChildren)
      c->SetRnrSelfChildren(rnr_self, rnr_children);
}

////////////////////////////////////////////////////////////////////////////////
/// Protected member function called from REveElement::Annihilate().

void REveElement::AnnihilateRecursively()
{
   // projected  were already destroyed in REveElement::Anihilate(), now only clear its list
   REveProjectable* pable = dynamic_cast<REveProjectable*>(this);
   if (pable && pable->HasProjecteds())
      pable->ClearProjectedList();

   // same as REveElement::RemoveElementsInternal(), except parents are ignored
   RemoveElementsLocal();
   for (auto &c : fChildren)
      c->AnnihilateRecursively();

   fChildren.clear();

   fDestructing = kAnnihilate;
   PreDeleteElement();

   delete this;
}

////////////////////////////////////////////////////////////////////////////////
/// Optimized destruction without check of reference-count.
/// Parents are not notified about child destruction.
/// The method should only be used when an element does not have
/// more than one parent -- otherwise an exception is thrown.

void REveElement::Annihilate()
{
   static const REveException eh("REveElement::Annihilate ");

   fDestructing = kAnnihilate;

   // recursive annihilation of projecteds
   REveProjectable* pable = dynamic_cast<REveProjectable*>(this);
   if (pable && pable->HasProjecteds())
   {
      pable->AnnihilateProjecteds();
   }

   // detach from the parent
   if (fMother)
   {
      fMother->RemoveElement(this);
   }

   // XXXX wont the above already start off the destruction cascade ?????

   AnnihilateRecursively();

   // XXXX ????? Annihilate flag ???? Is it different than regular remove ????
   // REX::gEve->Redraw3D();
}

////////////////////////////////////////////////////////////////////////////////
/// Annihilate elements.

void REveElement::AnnihilateElements()
{
   while (!fChildren.empty())
   {
      auto c = fChildren.front();
      c->Annihilate();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Destroy this element. Throws an exception if deny-destroy is in force.
/// This method should be called instead of a destructor.
/// Note that an exception will be thrown if the element has been
/// protected against destruction with IncDenyDestroy().

void REveElement::Destroy()
{
   static const REveException eh("REveElement::Destroy ");

   if (fDenyDestroy > 0)
      throw eh + TString::Format("element '%s' (%s*) %p is protected against destruction.",
                                 GetCName(), IsA()->GetName(), this);

   PreDeleteElement();
   delete this;
   REX::gEve->Redraw3D();
}

////////////////////////////////////////////////////////////////////////////////
/// Destroy this element. Prints a warning if deny-destroy is in force.

void REveElement::DestroyOrWarn()
{
   try
   {
      Destroy();
   }
   catch (REveException &exc)
   {
      ::Warning("REveElement::DestroyOrWarn", "Error while destroy element %p : %s", this, exc.what());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Destroy all children of this element.

void REveElement::DestroyElements()
{
   while (HasChildren())
   {
      auto c = fChildren.front();
      if (c->fDenyDestroy <= 0)
      {
         try {
            c->Destroy();
         }
         catch (REveException &exc) {
            ::Warning("REveElement::DestroyElements", "element destruction failed: '%s'.", exc.what());
            RemoveElement(c);
         }
      }
      else
      {
         if (gDebug > 0)
           ::Info("REveElement::DestroyElements", "element '%s' is protected against destruction, removing locally.", c->GetCName());
         RemoveElement(c);
      }
   }

   REX::gEve->Redraw3D();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns state of flag determining if the element will be
/// destroyed when reference count reaches zero.
/// This is true by default.

Bool_t REveElement::GetDestroyOnZeroRefCnt() const
{
   return fDestroyOnZeroRefCnt;
}

////////////////////////////////////////////////////////////////////////////////
/// Sets the state of flag determining if the element will be
/// destroyed when reference count reaches zero.
/// This is true by default.

void REveElement::SetDestroyOnZeroRefCnt(Bool_t d)
{
   fDestroyOnZeroRefCnt = d;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the number of times deny-destroy has been requested on
/// the element.

Int_t REveElement::GetDenyDestroy() const
{
   return fDenyDestroy;
}

////////////////////////////////////////////////////////////////////////////////
/// Increases the deny-destroy count of the element.
/// Call this if you store an external pointer to the element.

void REveElement::IncDenyDestroy()
{
   ++fDenyDestroy;
}

////////////////////////////////////////////////////////////////////////////////
/// Decreases the deny-destroy count of the element.
/// Call this after releasing an external pointer to the element.

void REveElement::DecDenyDestroy()
{
   if (--fDenyDestroy <= 0)
      CheckReferenceCount("REveElement::DecDenyDestroy ");
}

////////////////////////////////////////////////////////////////////////////////
/// Set pickable state on the element and all its children.

void REveElement::SetPickableRecursively(Bool_t p)
{
   fPickable = p;
   for (auto &c: fChildren)
      c->SetPickableRecursively(p);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the master element - that is:
/// - master of projectable, if this is a projected;
/// - master of compound, if fCompound is set;
/// - master of mother, if kSCBTakeMotherAsMaster bit is set;
/// If non of the above is true, *this* is returned.

REveElement *REveElement::GetSelectionMaster()
{
   if (fSelectionMaster) return fSelectionMaster;

   REveProjected* proj = dynamic_cast<REveProjected*>(this);
   if (proj)
   {
      return dynamic_cast<REveElement*>(proj->GetProjectable())->GetSelectionMaster();
   }
   if (fCompound)
   {
      return fCompound->GetSelectionMaster();
   }
   if (TestCSCBits(kCSCBTakeMotherAsMaster) && fMother)
   {
      return fMother->GetSelectionMaster();
   }
   return this;
}

////////////////////////////////////////////////////////////////////////////////
/// Populate set impSelSet with derived / dependant elements.
///
/// If this is a REveProjectable, the projected replicas are added
/// to the set. Thus it does not have to be reimplemented for each
/// sub-class of REveProjected.
///
/// Note that this also takes care of projections of REveCompound
/// class, which is also a projectable.

void REveElement::FillImpliedSelectedSet(Set_t &impSelSet)
{
   REveProjectable* p = dynamic_cast<REveProjectable*>(this);
   if (p)
      p->AddProjectedsToSet(impSelSet);
}

////////////////////////////////////////////////////////////////////////////////
/// Call this if it is possible that implied-selection or highlight
/// has changed for this element or for implied-selection this
/// element is member of and you want to maintain consistent
/// selection state.
/// This can happen if you add elements into compounds in response
/// to user-interaction.

void REveElement::RecheckImpliedSelections()
{
   // XXXX MT 2019-01 --- RecheckImpliedSelections
   //
   // With removal of selection state from this class there might be some
   // corner cases requiring checking of implied-selected state in
   // selection/highlight objects.
   //
   // This could be done as part of begin / end changes on the EveManager level.
   //
   // See also those functions in TEveSelection.

   // if (fSelected || fImpliedSelected)
   //    REX::gEve->GetSelection()->RecheckImpliedSetForElement(this);

   // if (fHighlighted || fImpliedHighlighted)
   //    REX::gEve->GetHighlight()->RecheckImpliedSetForElement(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Add (bitwise or) given stamps to fChangeBits.
/// Register this element to REX::gEve as stamped.
/// This method is virtual so that sub-classes can add additional
/// actions. The base-class method should still be called (or replicated).

void REveElement::AddStamp(UChar_t bits)
{
   if (fDestructing == kNone && fScene && fScene->IsAcceptingChanges())
   {
      if (gDebug > 0)
         ::Info(Form("%s::AddStamp", GetCName()), "%d + (%d) -> %d", fChangeBits, bits, fChangeBits | bits);

      if (fChangeBits == 0)
      {
         fScene->SceneElementChanged(this);
      }

      fChangeBits |= bits;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Write transformation Matrix to render data
////////////////////////////////////////////////////////////////////////////////

void REveElement::BuildRenderData()
{
   if (fMainTrans.get())
   {
      fRenderData->SetMatrix(fMainTrans->Array());
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Convert Bool_t to string - kTRUE or kFALSE.
/// Needed in WriteVizParams().

const std::string& REveElement::ToString(Bool_t b)
{
   static const std::string true_str ("kTRUE");
   static const std::string false_str("kFALSE");

   return b ? true_str : false_str;
}

////////////////////////////////////////////////////////////////////////////////
/// Write core json. If rnr_offset is negative, render data shall not be
/// written.
/// Returns number of bytes written into binary render data.

Int_t REveElement::WriteCoreJson(nlohmann::json &j, Int_t rnr_offset)
{
   j["_typename"]  = IsA()->GetName();
   j["fName"]      = fName;
   j["fTitle"]     = fTitle;
   j["fElementId"] = GetElementId();
   j["fMotherId"]  = get_mother_id();
   j["fSceneId"]   = get_scene_id();
   j["fMasterId"]  = GetSelectionMaster()->GetElementId();

   j["fRnrSelf"]     = GetRnrSelf();
   j["fRnrChildren"] = GetRnrChildren();

   j["fMainColor"]        = GetMainColor();
   j["fMainTransparency"] = GetMainTransparency();
   j["fPickable"]         = fPickable;

   Int_t ret = 0;

   if (rnr_offset >= 0) {
      BuildRenderData();

      if (fRenderData) {
         nlohmann::json rd = {};

         rd["rnr_offset"] = rnr_offset;
         rd["rnr_func"] = fRenderData->GetRnrFunc();
         rd["vert_size"] = fRenderData->SizeV();
         rd["norm_size"] = fRenderData->SizeN();
         rd["index_size"] = fRenderData->SizeI();
         rd["trans_size"] = fRenderData->SizeT();

         j["render_data"] = rd;

         ret = fRenderData->GetBinarySize();
      }
   }

   return ret;
}
