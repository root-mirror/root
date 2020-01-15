/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/Browsable/RProvider.hxx>

#include "TObject.h"
#include <ROOT/RCanvas.hxx>
#include <ROOT/RObjectDrawable.hxx>

using namespace ROOT::Experimental;

/** Provider for drawing of ROOT7 classes */

class TObjectDraw7Provider : public Browsable::RProvider {
public:
   TObjectDraw7Provider()
   {
      RegisterDraw7(nullptr, [] (std::shared_ptr<RPadBase> &subpad, std::unique_ptr<Browsable::RHolder> &obj, const std::string &opt) -> bool {
         // here clear ownership is required
         // If it possible, TObject will be cloned by TObjectHolder
         auto tobj = obj->get_shared<TObject>();
         if (!tobj) return false;

         if (subpad->NumPrimitives() > 0) {
            subpad->Wipe();
            subpad->GetCanvas()->Modified();
            subpad->GetCanvas()->Update(true);
         }

         subpad->Draw<RObjectDrawable>(tobj, opt);
         return true;
      });

   }

} newTObjectDraw7Provider;
