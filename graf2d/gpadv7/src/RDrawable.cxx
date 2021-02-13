/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "ROOT/RDrawable.hxx"
#include "ROOT/RDisplayItem.hxx"
#include "ROOT/RMenuItems.hxx"
#include "ROOT/RLogger.hxx"
#include "ROOT/RCanvas.hxx"

#include "TClass.h"
#include "TROOT.h"

#include <cassert>
#include <string>
#include <sstream>

using namespace ROOT::Experimental;

// destructor, pin vtable
RDrawable::~RDrawable() = default;

/////////////////////////////////////////////////////////////////////
// Fill context menu items for the ROOT class

void RDrawable::PopulateMenu(RMenuItems &items)
{
   auto isA = TClass::GetClass(typeid(*this));
   if (isA)
      items.PopulateObjectMenu(this, isA);
}

/////////////////////////////////////////////////////////////////////
// Execute command for the drawable

void RDrawable::Execute(const std::string &exec)
{
   auto isA = TClass::GetClass(typeid(*this));
   if (!isA) return;

   std::stringstream cmd;
   cmd << "((" << isA->GetName() << " *) " << std::hex << std::showbase << (size_t)this << ")->" << exec << ";";
   R__LOG_DEBUG(0, GPadLog()) << "RDrawable::Execute Obj " << this << " cmd " << exec;
   gROOT->ProcessLine(cmd.str().c_str());
}

/////////////////////////////////////////////////////////////////////////////
/// Preliminary method which checks if drawable matches with given selector
/// Following selector are allowed:
///   "type" or "#id" or ".class_name"
///  Here type is drawable kind like 'rect' or 'pad'
///       id is drawable identifier, specified with RDrawable::SetId() method
///       class_name is drawable class name, specified with RDrawable::SetCssClass() method

bool RDrawable::MatchSelector(const std::string &selector) const
{
   return (selector == fCssType) || (!fCssClass.empty() && (selector == std::string(".") + fCssClass)) || (!fId.empty() && (selector == std::string("#") + fId));
}

/////////////////////////////////////////////////////////////////////////////
/// Creates display item for drawable
/// By default item contains drawable data itself

std::unique_ptr<RDisplayItem> RDrawable::Display(const RDisplayContext &ctxt)
{
   if (GetVersion() > ctxt.GetLastVersion())
      return std::make_unique<RDrawableDisplayItem>(*this);

   return nullptr;
}
