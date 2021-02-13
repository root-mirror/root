// @(#)root/eve7:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel, 2018

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_REveDataProxyBuilderBase
#define ROOT7_REveDataProxyBuilderBase

#include <ROOT/REveElement.hxx>
#include <ROOT/REveCompound.hxx>
#include <ROOT/REveDataCollection.hxx>

namespace ROOT {
namespace Experimental {

class REveViewContext;
class REveTrackPropagator;

class REveDataProxyBuilderBase
{
public:
   struct Product
   {
      std::string              m_viewType;
      const REveViewContext   *m_viewContext{nullptr};
      REveCompound            *m_elements{nullptr};

      Product(std::string viewType, const REveViewContext* c);
      virtual ~Product();
   };

   // ---------- const member functions ---------------------

   const REveViewContext&    Context()    const;
   REveDataCollection* Collection() const { return m_collection;  }

   // ---------- constructor/destructor  ---------------------

   REveDataProxyBuilderBase();
   virtual ~REveDataProxyBuilderBase() {}

   virtual void SetCollection(REveDataCollection*);
   //   virtual void SetInteractionList(REveDataInteractionList*, const std::string&);

   virtual void CollectionBeingDestroyed(const REveDataCollection*);

   void Build();
   // virtual void Build(REveElement* product);

   REveElement* CreateProduct(const std::string& viewType, const REveViewContext*);
   //  void removePerViewProduct(const REveViewContext* vc);

   void FillImpliedSelected(REveElement::Set_t& impSet);
   void ModelChanges(const REveDataCollection::Ids_t&);
   void CollectionChanged(const REveDataCollection*);

   void SetupElement(REveElement* el, bool color = true);
   void SetupAddElement(REveElement* el, REveElement* parent,  bool set_color = true);

   bool GetHaveAWindow() const { return m_haveWindow; }
   void SetHaveAWindow(bool);

   // const member functions
   virtual bool HaveSingleProduct() const { return true; }

protected:
   // Override this if visibility changes can cause (re)-creation of proxies.
   // Returns true if new proxies were created.
   virtual bool VisibilityModelChanges(int idx, REveElement*, const std::string& viewType, const REveViewContext*);

   virtual void Build(const REveDataCollection* iItem, REveElement* product, const REveViewContext*);
   virtual void BuildViewType(const REveDataCollection* iItem, REveElement* product, const std::string& viewType, const REveViewContext*);

   virtual void ModelChanges(const REveDataCollection::Ids_t&, Product*) = 0;
   virtual void FillImpliedSelected( REveElement::Set_t& /*impSet*/, Product*) {};
   virtual void LocalModelChanges(int idx, REveElement* el, const REveViewContext* ctx);

   virtual void Clean();
   virtual void CleanLocal();

   std::vector<Product*> m_products;

private:
   REveDataCollection *m_collection{nullptr};

   float                 m_layer;
   bool                  m_haveWindow;
   bool                  m_modelsChanged;
};

} // namespace Experimental
} // namespace ROOT
#endif
