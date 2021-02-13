/// \file RHistDraw7Provider.cxx
/// \ingroup rbrowser
/// \author Bertrand Bellenot <bertrand.bellenot@cern.ch>
/// \author Sergey Linev <S.Linev@gsi.de>
/// \date 2019-10-17
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback
/// is welcome!

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/Browsable/RProvider.hxx>
#include <ROOT/RHistDrawable.hxx>

#include <ROOT/RCanvas.hxx>

using namespace ROOT::Experimental::Browsable;

class RV7HistDrawProvider : public RProvider {
   template<class HistClass>
   void RegisterHistClass()
   {
      RegisterDraw7(TClass::GetClass<HistClass>(), [] (std::shared_ptr<ROOT::Experimental::RPadBase> &subpad, std::unique_ptr<RHolder> &obj, const std::string &) -> bool {
         auto hist = obj->get_shared<HistClass>();
         if (!hist) return false;

         if (subpad->NumPrimitives() > 0) {
            subpad->Wipe();
            subpad->GetCanvas()->Modified();
            subpad->GetCanvas()->Update(true);
         }

         subpad->Draw(hist);

         return true;
      });
   }

public:

   RV7HistDrawProvider()
   {
      RegisterHistClass<ROOT::Experimental::RH1D>();
      RegisterHistClass<ROOT::Experimental::RH2D>();
      RegisterHistClass<ROOT::Experimental::RH3D>();
   }

} newRV7HistDrawProvider;
