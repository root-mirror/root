/*************************************************************************
 * Copyright (C) 1995-2017, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "ROOT/RDisplayItem.hxx"

#include "ROOT/RDrawable.hxx"

#include "TString.h"
#include "TObject.h"

using namespace ROOT::Experimental;

////////////////////////////////////////////////////////////////////////////
/// Assign id using arbitrary pointer value
/// Typically drawable pointer should be used here

void RDisplayItem::SetObjectIDAsPtr(const void *ptr)
{
   SetObjectID(ObjectIDFromPtr(ptr));
}

////////////////////////////////////////////////////////////////////////////
/// Build full id, including prefix and object index

void RDisplayItem::BuildFullId(const std::string &prefix)
{
   SetObjectID(prefix + std::to_string(GetIndex()) + "_" + GetObjectID());
}

////////////////////////////////////////////////////////////////////////////
/// Construct fillid using pointer value

std::string RDisplayItem::ObjectIDFromPtr(const void *ptr)
{
   auto hash = TString::Hash(&ptr, sizeof(ptr));
   return std::to_string(hash);
}

///////////////////////////////////////////////////////////
/// destructor

RDrawableDisplayItem::~RDrawableDisplayItem()
{
   if (fDrawable)
      fDrawable->OnDisplayItemDestroyed(this);
}


///////////////////////////////////////////////////////////
/// Constructor

RIndirectDisplayItem::RIndirectDisplayItem(const RDrawable &dr)
{
   fAttr = &dr.fAttr;
   fCssClass = &dr.fCssClass;
   fId = &dr.fId;
}


///////////////////////////////////////////////////////////
/// destructor
TObjectDisplayItem::~TObjectDisplayItem()
{
   if (fOwner) delete fObject;
}


