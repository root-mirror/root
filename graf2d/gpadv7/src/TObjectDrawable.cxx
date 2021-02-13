/*************************************************************************
 * Copyright (C) 1995-2017, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/TObjectDrawable.hxx>

#include <ROOT/RDisplayItem.hxx>
#include <ROOT/RLogger.hxx>
#include <ROOT/RMenuItems.hxx>

#include "TArrayI.h"
#include "TColor.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TROOT.h"
#include "TStyle.h"

#include <exception>
#include <sstream>
#include <iostream>

using namespace ROOT::Experimental;

TObjectDrawable::TObjectDrawable(EKind kind, bool persistent) : RDrawable("tobject"), fKind(kind)
{
   if (!persistent) return;

   fOpts = "persistent";
   fObj = CreateSpecials(kind);
}

////////////////////////////////////////////////////////////////////
/// Convert TColor to RGB string for using with SVG

const char *TObjectDrawable::GetColorCode(TColor *col)
{
   static TString code;

   if (col->GetAlpha() == 1)
      code.Form("rgb(%d,%d,%d)", (int) (255*col->GetRed()), (int) (255*col->GetGreen()), (int) (255*col->GetBlue()));
   else
      code.Form("rgba(%d,%d,%d,%5.3f)", (int) (255*col->GetRed()), (int) (255*col->GetGreen()), (int) (255*col->GetBlue()), col->GetAlpha());

   return code.Data();
}

////////////////////////////////////////////////////////////////////
/// Create instance of requested special object

std::unique_ptr<TObject> TObjectDrawable::CreateSpecials(int kind)
{
   switch (kind) {
      case kColors: {
         // convert list of colors into strings
         auto arr = std::make_unique<TObjArray>();
         arr->SetOwner(kTRUE);

         auto cols = gROOT->GetListOfColors();
         for (int n = 0; n <= cols->GetLast(); ++n) {
            auto col = dynamic_cast<TColor *>(cols->At(n));
            if (!col) continue;
            auto code = TString::Format("%d=%s", n, GetColorCode(col));
            arr->Add(new TObjString(code));
         }

         return arr;
      }
      case kStyle: {  // create copy of gStyle
         return std::make_unique<TStyle>(*gStyle);
      }
      case kPalette: { // copy color palette

         auto arr = std::make_unique<TObjArray>();
         arr->SetOwner(kTRUE);

         auto palette = TColor::GetPalette();
         for (int n = 0; n < palette.GetSize(); ++n) {
            auto col = gROOT->GetColor(palette[n]);
            arr->Add(new TObjString(GetColorCode(col)));
         }

         return arr;
      }
   }
   return nullptr;
}


////////////////////////////////////////////////////////////////////
/// Create display item which will be delivered to the client

std::unique_ptr<RDisplayItem> TObjectDrawable::Display(const RDisplayContext &ctxt)
{
   if (GetVersion() > ctxt.GetLastVersion()) {
      if ((fKind == kObject) || (fOpts == "persistent"))
         return std::make_unique<TObjectDisplayItem>(fKind, fObj.get(), fOpts);

      auto specials = CreateSpecials(fKind);
      return std::make_unique<TObjectDisplayItem>(fKind, specials.release(), fOpts, true);
   }

   return nullptr;
}

void TObjectDrawable::PopulateMenu(RMenuItems &items)
{
   // fill context menu items for the ROOT class
   if (fKind == kObject)
      items.PopulateObjectMenu(fObj.get(), fObj.get()->IsA());
}

void TObjectDrawable::Execute(const std::string &exec)
{
   if (fKind != kObject) return;

   TObject *obj = fObj.get();

   std::string sub, ex = exec;
   if (ex.compare(0, 6, "xaxis#") == 0) {
      ex.erase(0,6);
      ex.insert(0, "GetXaxis()->");
   } else if (ex.compare(0, 6, "yaxis#") == 0) {
      ex.erase(0,6);
      ex.insert(0, "GetYaxis()->");
   } else if (ex.compare(0, 6, "zaxis#") == 0) {
      ex.erase(0,6);
      ex.insert(0, "GetZaxis()->");
   }

   std::stringstream cmd;
   cmd << "((" << obj->ClassName() << "* ) " << std::hex << std::showbase << (size_t)obj << ")->" << ex << ";";
   std::cout << "TObjectDrawable::Execute Obj " << obj->GetName() << "Cmd " << cmd.str() << std::endl;
   gROOT->ProcessLine(cmd.str().c_str());
}
