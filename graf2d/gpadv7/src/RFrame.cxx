/*************************************************************************
 * Copyright (C) 1995-2017, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "ROOT/RFrame.hxx"

#include "ROOT/RLogger.hxx"
#include "ROOT/RPadUserAxis.hxx"
#include "ROOT/RMenuItems.hxx"

#include "TROOT.h"

#include <cassert>
#include <sstream>

using namespace ROOT::Experimental;

////////////////////////////////////////////////////////////////////////////
/// Deprecated constructor, to be removed soon

RFrame::RFrame(std::vector<std::unique_ptr<RPadUserAxisBase>> &&coords) : RFrame()
{
   fUserCoord = std::move(coords);
}

////////////////////////////////////////////////////////////////////////////
/// Internal - extract range for specified axis

void RFrame::GetAxisRanges(unsigned ndim, const RAttrAxis &axis, RUserRanges &ranges) const
{
   if (axis.HasZoomMin())
      ranges.AssignMin(ndim, axis.GetZoomMin());

   if (axis.HasZoomMax())
      ranges.AssignMax(ndim, axis.GetZoomMax());
}

////////////////////////////////////////////////////////////////////////////
/// Internal - assign client zoomed range to specified axis

void RFrame::AssignZoomRange(unsigned ndim, RAttrAxis &axis, const RUserRanges &ranges)
{
   if (ranges.IsUnzoom(ndim)) {
      axis.ClearZoom();
   } else {
      if (ranges.HasMin(ndim))
         axis.SetZoomMin(ranges.GetMin(ndim));
      if (ranges.HasMax(ndim))
         axis.SetZoomMax(ranges.GetMax(ndim));
   }
}

////////////////////////////////////////////////////////////////////////////
/// Deprecated, to be removed soon

void RFrame::GrowToDimensions(size_t nDimensions)
{
   std::size_t oldSize = fUserCoord.size();
   if (oldSize >= nDimensions)
      return;
   fUserCoord.resize(nDimensions);
   for (std::size_t idx = oldSize; idx < nDimensions; ++idx)
      if (!fUserCoord[idx])
         fUserCoord[idx].reset(new RPadCartesianUserAxis);
}

////////////////////////////////////////////////////////////////////////////
/// Provide context menu items

void RFrame::PopulateMenu(RMenuItems & /* items */)
{
   // do not use online context menu for frame - make it fully client side
/*   auto is_x = items.GetSpecifier() == "x";
   auto is_y = items.GetSpecifier() == "y";

   if (is_x || is_y) {
      RAttrAxis &attr = is_x ? AttrX() : AttrY();
      std::string name = is_x ? "AttrX()" : "AttrY()";
      auto cl = TClass::GetClass<RAttrAxis>();
      auto log = attr.GetLog();
      items.AddChkMenuItem("Linear scale", "Set linear scale", !log, name + ".SetLog(0)", cl);
      items.AddChkMenuItem("Log scale", "Logarithmic scale", !log, name + ".SetLog(10)", cl);
   }
*/
}

////////////////////////////////////////////////////////////////////////////
/// Remember client range, can be used for drawing or stats box calculations

void RFrame::SetClientRanges(unsigned connid, const RUserRanges &ranges, bool ismainconn)
{
   if (ismainconn) {
      AssignZoomRange(0, AttrX(), ranges);
      AssignZoomRange(1, AttrY(), ranges);
      AssignZoomRange(2, AttrZ(), ranges);
   }

   if (fClientRanges.find(connid) == fClientRanges.end()) {
      RUserRanges ranges0;
      GetClientRanges(connid, ranges0);
      fClientRanges[connid] = ranges0;
   }

   fClientRanges[connid].Update(ranges);
}

////////////////////////////////////////////////////////////////////////////
/// Return ranges configured for the client

void RFrame::GetClientRanges(unsigned connid, RUserRanges &ranges)
{
   auto iter = fClientRanges.find(connid);

   if (iter != fClientRanges.end()) {
      ranges = iter->second;
   } else {
      GetAxisRanges(0, GetAttrX(), ranges);
      GetAxisRanges(1, GetAttrY(), ranges);
      GetAxisRanges(2, GetAttrZ(), ranges);
   }
}

