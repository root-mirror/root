// @(#)root/eve:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TEveProjectionManager.h"
#include "TEveManager.h"
#include "TEveProjectionBases.h"
#include "TEveCompound.h"

#include "TBuffer3D.h"
#include "TBuffer3DTypes.h"
#include "TVirtualViewer3D.h"

#include "TClass.h"

#include <list>

/** \class TEveProjectionManager
\ingroup TEve
Manager class for steering of projections and managing projected objects.

Recursively projects TEveElement's and draws axis in the projected
scene. It enables to interactively set TEveProjection parameters
and updates projected scene accordingly.
*/

ClassImp(TEveProjectionManager);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TEveProjectionManager::TEveProjectionManager(TEveProjection::EPType_e type):
   TEveElementList("TEveProjectionManager",""),
   TAttBBox(),
   fProjection  (0),
   fCurrentDepth(0),
   fImportEmpty (kFALSE)
{
   for (Int_t i = 0; i < TEveProjection::kPT_End; ++i)
      fProjections[i] = 0;

   if (type != TEveProjection::kPT_Unknown)
      SetProjection(type);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.
/// Destroys also dependent elements.

TEveProjectionManager::~TEveProjectionManager()
{
   for (Int_t i = 0; i < TEveProjection::kPT_End; ++i)
   {
      delete fProjections[i];
   }
   while ( ! fDependentEls.empty())
   {
      fDependentEls.front()->Destroy();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Add el as dependent element.

void TEveProjectionManager::AddDependent(TEveElement* el)
{
   fDependentEls.push_back(el);
}

////////////////////////////////////////////////////////////////////////////////
/// Remove el as dependent element.

void TEveProjectionManager::RemoveDependent(TEveElement* el)
{
   fDependentEls.remove(el);
}

////////////////////////////////////////////////////////////////////////////////
/// Updates name to have consistent information with projection.

void TEveProjectionManager::UpdateName()
{
   if (fProjection->Is2D())
      SetName(Form ("%s (%3.1f)", fProjection->GetName(), fProjection->GetDistortion()*1000));
   else
      SetName(fProjection->GetName());
}

////////////////////////////////////////////////////////////////////////////////
/// Set projection type and distortion.

void TEveProjectionManager::SetProjection(TEveProjection::EPType_e type)
{
   static const TEveException eH("TEveProjectionManager::SetProjection ");

   if (fProjections[type] == 0)
   {
      switch (type)
      {
         case TEveProjection::kPT_RPhi:
         {
            fProjections[type] = new TEveRPhiProjection();
            break;
         }
         case TEveProjection::kPT_RhoZ:
         {
            fProjections[type] = new TEveRhoZProjection();
            break;
         }
         case TEveProjection::kPT_3D:
         {
            fProjections[type] = new TEve3DProjection();
            break;
         }
         default:
            throw eH + "projection type not valid.";
            break;
      }
   }

   if (fProjection && fProjection->Is2D() != fProjections[type]->Is2D())
   {
      throw eH + "switching between 2D and 3D projections not implemented.";
   }

   fProjection = fProjections[type];
   fProjection->SetCenter(fCenter);
   UpdateName();
}

////////////////////////////////////////////////////////////////////////////////
/// Set projection center and rebuild projected scene.

void TEveProjectionManager::SetCenter(Float_t x, Float_t y, Float_t z)
{
   fCenter.Set(x, y, z);
   fProjection->SetCenter(fCenter);
   ProjectChildren();
}

////////////////////////////////////////////////////////////////////////////////
/// React to element being pasted or dnd-ed.
/// Return true if redraw is needed (virtual method).

Bool_t TEveProjectionManager::HandleElementPaste(TEveElement* el)
{
   List_t::size_type n_children  = fChildren.size();
   ImportElements(el);
   return n_children != fChildren.size();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns true if element el should be imported.
///
/// Behaviour depends on the value of the fImportEmpty member:
///   false - el or any of its children must be projectable (default);
///   true  - always import.

Bool_t TEveProjectionManager::ShouldImport(TEveElement* el)
{
   if (fImportEmpty)
      return kTRUE;

   if (el->IsA() != TEveElementList::Class() && el->IsA()->InheritsFrom(TEveProjectable::Class()))
      return kTRUE;
   for (List_i i=el->BeginChildren(); i!=el->EndChildren(); ++i)
      if (ShouldImport(*i))
         return kTRUE;
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Update dependent elements' bounding box and mark scenes
/// containing element root or its children as requiring a repaint.

void TEveProjectionManager::UpdateDependentElsAndScenes(TEveElement* root)
{
   for (List_i i=fDependentEls.begin(); i!=fDependentEls.end(); ++i)
   {
      TAttBBox* bbox = dynamic_cast<TAttBBox*>(*i);
      if (bbox)
         bbox->ComputeBBox();
   }

   List_t scenes;
   root->CollectSceneParentsFromChildren(scenes, 0);
   gEve->ScenesChanged(scenes);
}

////////////////////////////////////////////////////////////////////////////////
/// If el is TEveProjectable add projected instance else add plain
/// TEveElementList to parent. Call the same function on el's
/// children.
///
/// Returns the projected replica of el. Can be 0, if el and none of
/// its children are projectable.

TEveElement* TEveProjectionManager::ImportElementsRecurse(TEveElement* el,
                                                          TEveElement* parent)
{
   static const TEveException eh("TEveProjectionManager::ImportElementsRecurse ");

   TEveElement *new_el = 0;

   if (ShouldImport(el))
   {
      TEveProjected   *new_pr = 0;
      TEveProjectable *pble   = dynamic_cast<TEveProjectable*>(el);
      if (pble)
      {
         new_el = (TEveElement*) pble->ProjectedClass(fProjection)->New();
         new_pr = dynamic_cast<TEveProjected*>(new_el);
         new_pr->SetProjection(this, pble);
         new_pr->SetDepth(fCurrentDepth);
      }
      else
      {
         new_el = new TEveElementList;
      }
      new_el->SetElementName (Form("%s [P]", el->GetElementName()));
      new_el->SetElementTitle(Form("Projected replica.\n%s", el->GetElementTitle()));
      new_el->SetRnrSelf     (el->GetRnrSelf());
      new_el->SetRnrChildren (el->GetRnrChildren());
      new_el->SetPickable    (el->IsPickable());
      parent->AddElement(new_el);

      TEveCompound *cmpnd    = dynamic_cast<TEveCompound*>(el);
      TEveCompound *cmpnd_pr = dynamic_cast<TEveCompound*>(new_el);
      for (List_i i=el->BeginChildren(); i!=el->EndChildren(); ++i)
      {
         TEveElement* child_pr = ImportElementsRecurse(*i, new_el);
         if (cmpnd && (*i)->GetCompound() == cmpnd)
            child_pr->SetCompound(cmpnd_pr);
      }
   }

   return new_el;
}

////////////////////////////////////////////////////////////////////////////////
/// Recursively import elements and apply projection to the newly
/// imported objects.
///
/// If ext_list is not 0 the new element is also added to the list.
/// This simplifies construction of complex views where projected
/// elements are distributed into several scenes for optimization of
/// updates and rendering.
///
/// Returns the projected replica of el. Can be 0, if el and none of
/// its children are projectable.

TEveElement* TEveProjectionManager::ImportElements(TEveElement* el,
                                                   TEveElement* ext_list)
{
   TEveElement* new_el = ImportElementsRecurse(el, this);
   if (new_el)
   {
      AssertBBox();
      ProjectChildrenRecurse(new_el);
      AssertBBoxExtents(0.1);
      StampTransBBox();

      UpdateDependentElsAndScenes(new_el);

      if (ext_list)
         ext_list->AddElement(new_el);
   }
   return new_el;
}

////////////////////////////////////////////////////////////////////////////////
/// Recursively import elements and apply projection to the newly
/// imported objects.
///
/// The proj_parent argument should be a projected replica of parent
/// of element 'el'. This allows to insert projected children of
/// a given element when they are added after the projection has
/// been already performed on the parent.
/// This is called from TEveElement::ProjectChild().
///
/// Returns the projected replica of el. Can be 0, if el and none of
/// its children are projectable.

TEveElement* TEveProjectionManager::SubImportElements(TEveElement* el,
                                                      TEveElement* proj_parent)
{
   TEveElement* new_el = ImportElementsRecurse(el, proj_parent);
   if (new_el)
   {
      AssertBBox();
      ProjectChildrenRecurse(new_el);
      AssertBBoxExtents(0.1);
      StampTransBBox();

      UpdateDependentElsAndScenes(new_el);
   }
   return new_el;
}

////////////////////////////////////////////////////////////////////////////////
/// Recursively import children elements of el and apply projection
/// to the newly imported objects.
///
/// The proj_parent argument should be a projected replica of
/// element 'el'. This allows to insert projected children of
/// a given element when they are added after the projection has
/// been already performed on the parent.
/// This is called from TEveElement::ProjectChild().
///
/// Returns the projected replica of el. Can be 0, if el and none of
/// its children are projectable.

Int_t TEveProjectionManager::SubImportChildren(TEveElement* el, TEveElement* proj_parent)
{
   List_t new_els;
   for (List_i i = el->BeginChildren(); i != el->EndChildren(); ++i)
   {
      TEveElement* new_el = ImportElementsRecurse(*i, proj_parent);
      if (new_el)
         new_els.push_back(new_el);
   }

   if ( ! new_els.empty())
   {
      AssertBBox();
      for (List_i i = new_els.begin(); i != new_els.end(); ++i)
      {
         ProjectChildrenRecurse(*i);
      }
      AssertBBoxExtents(0.1);
      StampTransBBox();

      UpdateDependentElsAndScenes(proj_parent);
   }
   return (Int_t) new_els.size();
}

////////////////////////////////////////////////////////////////////////////////
/// Project el (via TEveProjected::UpdateProjection()) and recurse
/// through el's children.
/// Bounding-box is updated along the recursion.

void TEveProjectionManager::ProjectChildrenRecurse(TEveElement* el)
{
   TEveProjected* pted = dynamic_cast<TEveProjected*>(el);
   if (pted)
   {
      pted->UpdateProjection();
      TAttBBox* bb = dynamic_cast<TAttBBox*>(pted);
      if (bb)
      {
         Float_t* b = bb->AssertBBox();
         BBoxCheckPoint(b[0], b[2], b[4]);
         BBoxCheckPoint(b[1], b[3], b[5]);
      }
      el->ElementChanged(kFALSE);
   }

   for (List_i i=el->BeginChildren(); i!=el->EndChildren(); ++i)
      ProjectChildrenRecurse(*i);
}

////////////////////////////////////////////////////////////////////////////////
/// Project all children recursively, update bounding-box and notify
/// TEveManger about the scenes that have been changed.

void TEveProjectionManager::ProjectChildren()
{
   BBoxInit();
   for (List_i i=BeginChildren(); i!=EndChildren(); ++i)
      ProjectChildrenRecurse(*i);
   AssertBBoxExtents(0.1);
   StampTransBBox();

   UpdateDependentElsAndScenes(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Virtual from TAttBBox; fill bounding-box information.
///
/// The bounding-box information is kept coherent during addition of
/// projected elements and projection parameter updates. This is
/// called only in case the manager has not been populated at all.

void TEveProjectionManager::ComputeBBox()
{
   static const TEveException eH("TEveProjectionManager::ComputeBBox ");

   if (HasChildren() == kFALSE) {
      BBoxZero();
      return;
   }

   BBoxInit();
}
