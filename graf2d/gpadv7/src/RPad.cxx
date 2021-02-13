/*************************************************************************
 * Copyright (C) 1995-2017, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "ROOT/RPad.hxx"

#include "ROOT/RLogger.hxx"
#include <ROOT/RPadDisplayItem.hxx>
#include <ROOT/RCanvas.hxx>

#include <cassert>
#include <limits>

using namespace ROOT::Experimental;

/////////////////////////////////////////////////////////////////////////////////////////////////

ROOT::Experimental::RPad::~RPad() = default;

/////////////////////////////////////////////////////////////////////////////////////////////////
/// Create pad display item

std::unique_ptr<RDisplayItem> RPad::Display(const RDisplayContext &ctxt)
{
   auto paditem = std::make_unique<RPadDisplayItem>();

   RDisplayContext subctxt(ctxt.GetCanvas(), this, ctxt.GetLastVersion());

   DisplayPrimitives(*paditem.get(), subctxt);

   paditem->SetPadPosSize(&fPos, &fSize);

   return paditem;
}
