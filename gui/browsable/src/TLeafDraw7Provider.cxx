/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TLeafProvider.hxx"

#include <ROOT/RCanvas.hxx>
#include <ROOT/TObjectDrawable.hxx>

using namespace ROOT::Experimental;

/** Provider for drawing of ROOT7 classes */

class TLeafDraw7Provider : public TLeafProvider<void> {
public:
   bool AddHist(std::shared_ptr<RPadBase> &subpad, TH1 *hist, const std::string &opt)
   {
      if (!hist)
         return false;

      if (subpad->NumPrimitives() > 0) {
         subpad->Wipe();
         subpad->GetCanvas()->Modified();
         subpad->GetCanvas()->Update(true);
      }

      std::shared_ptr<TH1> shared;
      shared.reset(hist);

      subpad->Draw<ROOT::Experimental::TObjectDrawable>(shared, opt);
      subpad->GetCanvas()->Update(true);
      return true;
   }

   TLeafDraw7Provider()
   {
      RegisterDraw7(TLeaf::Class(), [this](std::shared_ptr<RPadBase> &subpad, std::unique_ptr<RHolder> &obj, const std::string &opt) -> bool {
         return AddHist(subpad, DrawLeaf(obj), opt);
      });

      RegisterDraw7(TBranchElement::Class(), [this](std::shared_ptr<RPadBase> &subpad, std::unique_ptr<RHolder> &obj, const std::string &opt) -> bool {
         return AddHist(subpad, DrawBranchElement(obj), opt);
      });

      RegisterDraw7(TVirtualBranchBrowsable::Class(), [this](std::shared_ptr<RPadBase> &subpad, std::unique_ptr<RHolder> &obj, const std::string &opt) -> bool {
         return AddHist(subpad, DrawBranchBrowsable(obj), opt);
      });

   }

} newTLeafDraw7Provider;

