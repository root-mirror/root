// @(#)root/eve7:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2020

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/REveDataSimpleProxyBuilder.hxx>

#include <ROOT/REveDataCollection.hxx>
#include <ROOT/REveCompound.hxx>
#include <ROOT/REveScene.hxx>
#include <cassert>
#include <TClass.h>

using namespace ROOT::Experimental;

REveDataSimpleProxyBuilder::REveDataSimpleProxyBuilder()
{
}


REveDataSimpleProxyBuilder::~REveDataSimpleProxyBuilder()
{
}

void
REveDataSimpleProxyBuilder::Clean()
{
   for (auto &p: m_products)
   {
      if (p->m_elements)
      {
         REveElement *elms = p->m_elements;
         for (auto &c: elms->RefChildren())
            c->DestroyElements();
      }
   }

   CleanLocal();
}

//______________________________________________________________________________

void
REveDataSimpleProxyBuilder::Build(const REveDataCollection *collection,
                                  REveElement* product, const REveViewContext* vc)
{
   // printf("REveDataSimpleProxyBuilder::Build %s %d\n", collection->GetCName(), collection->GetNItems());
   auto size = collection->GetNItems();
   auto pIdx = product->RefChildren().begin();
   for (int index = 0; index < size; ++index)
   {
      const REveDataItem* di = Collection()->GetDataItem(index);
      REveElement *itemHolder = nullptr;

      if (index <  product->NumChildren())
      {
         itemHolder = *pIdx;
         itemHolder->SetRnrSelfChildren(true, true);
         ++pIdx;
      }
      else
      {
         itemHolder = CreateCompound(true, true);
         itemHolder->SetMainColor(collection->GetMainColor());
         itemHolder->SetName(Form("%s %d", collection->GetCName(), index));

         product->AddElement(itemHolder);
      }


      if (di->GetRnrSelf() && !di->GetFiltered())
      {
         Build(collection->GetDataPtr(index), index, itemHolder, vc);
      }
   }
}

void
REveDataSimpleProxyBuilder::BuildViewType(const REveDataCollection* collection,
                                          REveElement* product, const std::string& viewType, const REveViewContext* vc)
{
   auto size = collection->GetNItems();
   auto pIdx = product->RefChildren().begin();
   for (int index = 0; index < size; ++index)
   {
      auto di = Collection()->GetDataItem(index);
      REveElement* itemHolder = nullptr;

      if (index <  product->NumChildren())
      {
         itemHolder = *pIdx;
         itemHolder->SetRnrSelfChildren(true, true);
         ++pIdx;
      }
      else
      {
         itemHolder = CreateCompound(true, true);
         itemHolder->SetMainColor(collection->GetMainColor());
         itemHolder->SetName(Form("%s %d", collection->GetCName(), index));

         product->AddElement(itemHolder);
      }


      if (di->GetRnrSelf() && !di->GetFiltered())
      {
         BuildViewType(collection->GetDataPtr(index), index, itemHolder, viewType, vc);
      }
      }
}

//______________________________________________________________________________

namespace
{
   void applyColorAttrToChildren(REveElement* p) {
      for (auto &it: p->RefChildren())
      {
         REveElement* c = it;
         if (c->GetMainColor() != p->GetMainColor())
         {
            c->SetMainColor(p->GetMainColor());
            // printf("apply color %d to %s\n", p->GetMainColor(), c->GetCName());
         }
         applyColorAttrToChildren(c);
      }
   }
}

void
REveDataSimpleProxyBuilder::ModelChanges(const REveDataCollection::Ids_t& iIds, Product* p)
{
   // printf("REveDataSimple ProxyBuilderBase::ModelChanges >>>>> (%p)  %s \n", (void*)this, Collection()->GetCName());
   REveElement* elms = p->m_elements;
   assert(Collection() && static_cast<int>(Collection()->GetNItems()) <= elms->NumChildren() && "can not use default modelChanges implementation");

   for (auto itemIdx: iIds)
   {
      const REveDataItem* item = Collection()->GetDataItem(itemIdx);

      // printf("Edit compound for item index %d \n", itemIdx);
      // imitate FWInteractionList::modelChanges
      auto itElement = elms->RefChildren().begin();
      std::advance(itElement, itemIdx);
      REveElement* comp = *itElement;
      bool visible = ((!item->GetFiltered()) && item->GetRnrSelf()) && Collection()->GetRnrSelf();
      comp->SetRnrSelf(visible);
      comp->SetRnrChildren(visible);

      // printf("comapre %d %d\n", item->GetMainColor(), comp->GetMainColor());
      if (item->GetMainColor() != comp->GetMainColor()) {
         //printf("ffffffffffffffffffffffff set color to comp \n");
         comp->SetMainColor(item->GetMainColor());

      }
      applyColorAttrToChildren(comp);

      if (VisibilityModelChanges(itemIdx, comp, p->m_viewType, p->m_viewContext))
      {
         elms->ProjectChild(comp);
         // printf("---REveDataProxyBuilderBase project child\n");
      }
      else
      {
         LocalModelChanges(itemIdx, comp, p->m_viewContext);
      }
   }
}


REveCompound*
REveDataSimpleProxyBuilder::CreateCompound(bool set_color, bool propagate_color_to_all_children)
{
   REveCollectionCompound *c = new REveCollectionCompound(Collection());
   c->CSCImplySelectAllChildren();
   c->SetPickable(true);
   if (set_color)
   {
      c->SetMainColor(Collection()->GetMainColor());
      c->SetMainTransparency(Collection()->GetMainTransparency());
   }
   if (propagate_color_to_all_children)
   {
      c->CSCApplyMainColorToAllChildren();
      c->CSCApplyMainTransparencyToAllChildren();
   }
   else
   {
      c->CSCApplyMainColorToMatchingChildren();
      c->CSCApplyMainTransparencyToMatchingChildren();
   }
   return c;
}
//______________________________________________________________________________

bool
REveDataSimpleProxyBuilder::VisibilityModelChanges(int idx, REveElement* iCompound, const std::string& viewType, const REveViewContext* vc)
{
   const REveDataItem *item = Collection()->GetDataItem(idx);
   bool returnValue = false;
   if (item->GetVisible() && iCompound->NumChildren() == 0) {
      if (HaveSingleProduct())
         Build(Collection()->GetDataPtr(idx), idx, iCompound, vc);
      else
         BuildViewType(Collection()->GetDataPtr(idx), idx, iCompound, viewType, vc);
      returnValue = true;
   }
   return returnValue;
}

//______________________________________________________________________________

void
REveDataSimpleProxyBuilder::FillImpliedSelected(REveElement::Set_t& impSet, Product*p)
{
    REveElement* elms = p->m_elements;
    for (auto &s: Collection()->GetItemList()->RefSelectedSet()) {

      auto it = elms->RefChildren().begin();
      std::advance(it, s);
      REveElement* comp = *it;
      comp->FillImpliedSelectedSet(impSet);
   }
}

//==============================================================================

//==============================================================================

//==============================================================================
REveCollectionCompound::REveCollectionCompound(REveDataCollection* collection) : fCollection(collection)
{
   SetSelectionMaster(collection);
   // deny destroy ??
}

//______________________________________________________________________________

REveCollectionCompound::~REveCollectionCompound()
{
   // deny destroy ??
}

//______________________________________________________________________________


REveElement* REveCollectionCompound::GetSelectionMaster()
{
   if (!fCollection->GetScene()->IsAcceptingChanges()) return fCollection->GetItemList();
   fCollection->GetItemList()->RefSelectedSet().clear();

   auto m = GetMother();
   int idx = 0;
   for (auto &c : m->RefChildren()) {
      REveElement* ctest = c;
      if (ctest == this)
      {
         fCollection->GetItemList()->RefSelectedSet().insert(idx);
         break;
      }
      ++idx;
   }

   return fCollection->GetItemList();
}
