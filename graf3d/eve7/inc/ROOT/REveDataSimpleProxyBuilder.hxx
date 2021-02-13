// @(#)root/eve7:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel, 2019

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_REveDataProxySimpleBuilder
#define ROOT7_REveDataProxySimpleBuilder

#include <ROOT/REveDataProxyBuilderBase.hxx>

namespace ROOT {
namespace Experimental {

class REveDataCollection;
class REveElement;

class REveDataSimpleProxyBuilder : public REveDataProxyBuilderBase
{
public:
   REveDataSimpleProxyBuilder();
   virtual ~REveDataSimpleProxyBuilder();

protected:
   void Build(const REveDataCollection* iCollection, REveElement* product, const REveViewContext*) override;

   void BuildViewType(const REveDataCollection* iCollection, REveElement* product, const std::string& viewType, const REveViewContext*) override;

   // Called once for every item in collection, the void* points to the
   // item properly offset in memory.
   virtual void Build(const void* data, int index, REveElement* iCollectionHolder, const REveViewContext*) = 0;
   virtual void BuildViewType(const void* data, int index, REveElement* iCollectionHolder, const std::string& viewType, const REveViewContext*) = 0;

   void ModelChanges(const REveDataCollection::Ids_t& iIds, Product* p) override;
   void FillImpliedSelected(REveElement::Set_t& impSet, Product* p) override;
   void Clean() override; // Utility
   REveCompound* CreateCompound(bool set_color=true, bool propagate_color_to_all_children=false);


private:
   REveDataSimpleProxyBuilder(const REveDataSimpleProxyBuilder&); // stop default

   const REveDataSimpleProxyBuilder& operator=(const REveDataSimpleProxyBuilder&); // stop default

   bool VisibilityModelChanges(int idx, REveElement*, const std::string& viewType, const REveViewContext*) override;

};
//==============================================================================

class REveCollectionCompound : public REveCompound // ?? Should this be in as REveDataSimpleProxyBuilder.hxx ?????
{
private:
   REveDataCollection* fCollection {nullptr};

public:
   REveCollectionCompound(REveDataCollection* c);
   virtual ~REveCollectionCompound();
   //   Int_t WriteCoreJson(nlohmann::json &cj, Int_t rnr_offset) override;

   // virtual REveElement* GetSelectionMaster(const bool &secondary = false, const std::set<int>& secondary_idcs = {});
   virtual REveElement* GetSelectionMaster() override;
};


} // namespace Experimental
} // namespace ROOT

#endif
