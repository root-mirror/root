// @(#)root/eve7:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2020

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_REveTableProxyBuilder
#define ROOT7_REveTableProxyBuilder

#include <ROOT/REveDataProxyBuilderBase.hxx>

namespace ROOT {
namespace Experimental {

class REveDataTable;
class REveTableInfo;

class REveTableProxyBuilder : public REveDataProxyBuilderBase
{
private:
   REveDataTable* fTable; // cached

protected:
   void Clean() override;

public:
   REveTableProxyBuilder();
   virtual ~REveTableProxyBuilder();

   virtual bool WillHandleInteraction() const { return true; }

   using REveDataProxyBuilderBase::ModelChanges;
   virtual void ModelChanges(const REveDataCollection::Ids_t&, REveDataProxyBuilderBase::Product* p) override;

   using REveDataProxyBuilderBase::Build;
   virtual void Build(const REveDataCollection* collection, REveElement* product, const REveViewContext* context) override;

   void SetCollection(REveDataCollection*) override;
   void ConfigChanged();
};

} // Experimental
} // ROOT

#endif
