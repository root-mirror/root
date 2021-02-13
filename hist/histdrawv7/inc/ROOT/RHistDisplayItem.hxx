/// \file ROOT/RHistDisplayItem.h
/// \ingroup HistDraw ROOT7
/// \author Sergey Linev <s.linev@gsi.de>
/// \date 2020-06-25
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RHistDisplayItem
#define ROOT7_RHistDisplayItem

#include <ROOT/RDisplayItem.hxx>

#include <vector>

namespace ROOT {
namespace Experimental {

class RAxisBase;

class RHistDisplayItem : public RIndirectDisplayItem {
   std::vector<const RAxisBase *> fAxes;   ///< histogram axes, only temporary pointers
   std::vector<int> fIndicies;             ///< [left,right,step] for each axes
   std::vector<double> fBinContent;        ///< extracted bins values
   double fContMin{0.};                    ///< minimum content value
   double fContMinPos{0.};                 ///< minimum positive value
   double fContMax{0.};                    ///< maximum content value

public:
   RHistDisplayItem() = default;

   RHistDisplayItem(const RDrawable &dr);

   void AddAxis(const RAxisBase *axis, int left = -1, int right = -1, int step = 1)
   {
      fAxes.emplace_back(axis);
      fIndicies.emplace_back(left);
      fIndicies.emplace_back(right);
      fIndicies.emplace_back(step);
   }

   auto &GetBinContent() { return fBinContent; }

   void SetContentMinMax(double min, double minpos, double max)
   {
      fContMin = min;
      fContMinPos = minpos;
      fContMax = max;
   }
};

} // namespace Experimental
} // namespace ROOT

#endif
