/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TLeafProvider.hxx"

#include "TVirtualPad.h"

/** Provider for drawing of ROOT6 classes */

class TLeafDraw6Provider : public TLeafProvider<void> {
public:

   bool AddHist(TVirtualPad *pad, TH1 *hist, const std::string &opt)
   {
      if (!hist)
         return false;

      pad->GetListOfPrimitives()->Clear();

      pad->GetListOfPrimitives()->Add(hist, opt.c_str());

      return true;
   }

   TLeafDraw6Provider()
   {
      RegisterDraw6(TLeaf::Class(), [this](TVirtualPad *pad, std::unique_ptr<RHolder> &obj, const std::string &opt) -> bool {
         return AddHist(pad, DrawLeaf(obj), opt);
      });

      RegisterDraw6(TBranchElement::Class(), [this](TVirtualPad *pad, std::unique_ptr<RHolder> &obj, const std::string &opt) -> bool {
         return AddHist(pad, DrawBranchElement(obj), opt);
      });

      RegisterDraw6(TVirtualBranchBrowsable::Class(), [this](TVirtualPad *pad, std::unique_ptr<RHolder> &obj, const std::string &opt) -> bool {
         return AddHist(pad, DrawBranchBrowsable(obj), opt);
      });
   }

} newTLeafDraw6Provider;

