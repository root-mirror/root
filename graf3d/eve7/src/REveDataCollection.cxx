// @(#)root/eve7:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007, 2018

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/REveDataCollection.hxx>
#include <ROOT/REveUtil.hxx>
#include <ROOT/REveSelection.hxx>
#include <ROOT/REveManager.hxx>

#include "TROOT.h"
#include "TMethod.h"
#include "TMethodArg.h"
#include "TColor.h"
#include "TClass.h"
#include "TList.h"
#include "TBaseClass.h"

#include <sstream>

using namespace ROOT::Experimental;


Color_t  REveDataCollection::fgDefaultColor  = kBlue;

//==============================================================================
//==============================================================================

REveDataItemList::REveDataItemList(const std::string& n, const std::string& t):
   REveElement(n,t)
{
   fAlwaysSecSelect = true;

   SetItemsChangeDelegate([&](REveDataItemList *collection, const REveDataCollection::Ids_t &ids) {
      REveDataItemList::DummyItemsChange(collection, ids);
   });

   SetFillImpliedSelectedDelegate([&](REveDataItemList *collection, REveElement::Set_t &impSelSet) {
      REveDataItemList::DummyFillImpliedSelected(collection, impSelSet);
   });

   SetupDefaultColorAndTransparency(REveDataCollection::fgDefaultColor, true, true);
}
//______________________________________________________________________________

void REveDataItemList::SetItemVisible(Int_t idx, Bool_t visible)
{
   fItems[idx]->SetRnrSelf(visible);
   ItemChanged(idx);
   StampObjProps();
}

//______________________________________________________________________________

void REveDataItemList::SetItemColorRGB(Int_t idx, UChar_t r, UChar_t g, UChar_t b)
{
   Color_t c = TColor::GetColor(r, g, b);
   fItems[idx]->SetMainColor(c);
   ItemChanged(idx);
   StampObjProps();
}
//______________________________________________________________________________

void REveDataItemList::ItemChanged(REveDataItem* iItem)
{
   int idx = 0;
   std::vector<int> ids;
   for (auto & chld : fItems)
   {
      if (chld == iItem) {
         ids.push_back(idx);
          fHandlerItemsChange( this , ids);
         return;
      }
      idx++;
   }
}

//______________________________________________________________________________

void REveDataItemList::ItemChanged(Int_t idx)
{
   std::vector<int> ids;
   ids.push_back(idx);
   fHandlerItemsChange( this , ids);
}

//______________________________________________________________________________

void REveDataItemList::FillImpliedSelectedSet( Set_t& impSelSet)
{
   /*
   printf("REveDataCollection::FillImpliedSelectedSet colecction setsize %zu\n",   RefSelectedSet().size());
   for (auto x : RefSelectedSet())
      printf("%d \n", x);
   */
   fHandlerFillImplied( this ,  impSelSet);
}

//______________________________________________________________________________


Int_t REveDataItemList::WriteCoreJson(nlohmann::json &j, Int_t rnr_offset)
{
   Int_t ret = REveElement::WriteCoreJson(j, rnr_offset);
   j["items"] =  nlohmann::json::array();
   for (auto & chld : fItems)
   {
      nlohmann::json i;
      i["fFiltered"] = chld->GetFiltered();
      i["fRnrSelf"] = chld->GetRnrSelf();
      i["fColor"] = chld->GetMainColor();
      j["items"].push_back(i);
   }

   return ret;
}

//______________________________________________________________________________

Bool_t REveDataItemList::SetRnrState(Bool_t iRnrSelf)
{
   Bool_t ret = REveElement::SetRnrState(iRnrSelf);
   std::vector<int> ids;

   for (size_t i = 0; i < fItems.size(); ++i ) {
      ids.push_back(i);
      fItems[i]->SetRnrSelf(fRnrSelf);
   }

   fHandlerItemsChange( this , ids);
   StampVisibility();
   StampObjProps();

   return ret;
}

//______________________________________________________________________________
void REveDataItemList::ProcessSelection(ElementId_t selectionId, bool multi, bool secondary, const std::set<int>& secondary_idcs)
{
   RefSelectedSet() = secondary_idcs;
   REveSelection* selection = (REveSelection*) ROOT::Experimental::gEve->FindElementById(selectionId);
   selection->NewElementPicked(GetElementId(), multi, secondary, secondary_idcs);
}

//______________________________________________________________________________
std::string REveDataItemList::GetHighlightTooltip(const std::set<int>& secondary_idcs) const
{
   if (secondary_idcs.empty()) {
    return GetName();
  }

   // print info for first selected index
   int idx = *secondary_idcs.begin();
   auto col = dynamic_cast<REveDataCollection*>(fMother);
   void* data = col->GetDataPtr(idx);
   std::string name = col->GetName();
   auto li = name.size();
   if (li && name[li-1] == 's')
   {
      name = name.substr(0, li-1);
   }

   std::string res;
   for (auto &z : secondary_idcs)
   {
      idx = z;
      data = col->GetDataPtr(idx);
      res +=  Form("%s %d",  name.c_str(), idx);
      for (auto &t : fTooltipExpressions) {
         std::string eval = t.fTooltipFunction.EvalExpr(data);
         res +=  Form("\n  %s = %s", t.fTooltipTitle.c_str(), eval.c_str());
      }
      res += "\n";
   }
   return res;
}

//______________________________________________________________________________
void REveDataItemList::AddTooltipExpression(const std::string &title, const std::string &expr)
{
   TTip tt;
   tt.fTooltipTitle = title;
   tt.fTooltipFunction.SetPrecision(2);
   auto col = dynamic_cast<REveDataCollection*>(fMother);
   auto icls = col->GetItemClass();
   tt.fTooltipFunction.SetExpressionAndType(expr, REveDataColumn::FT_Double, icls);
   fTooltipExpressions.push_back(tt);
}

//______________________________________________________________________________
void REveDataItemList::SetItemsChangeDelegate (ItemsChangeFunc_t handler_func)
{
   fHandlerItemsChange = handler_func;
}

//______________________________________________________________________________
void REveDataItemList::SetFillImpliedSelectedDelegate (FillImpliedSelectedFunc_t handler_func)
{
   fHandlerFillImplied = handler_func;
}

//______________________________________________________________________________
void REveDataItemList::DummyItemsChange(REveDataItemList*, const std::vector<int>&)
{
   if (gDebug) {
      printf("REveDataItemList::DummyItemsCahngeDelegate not implemented\n");
   }
}


//______________________________________________________________________________
void REveDataItemList::DummyFillImpliedSelected(REveDataItemList*, REveElement::Set_t&)
{
   if (gDebug) {
      printf("REveDataItemList::DummyFillImpliedSelectedDelegate not implemented\n");
   }
}

//==============================================================================
// REveDataCollection
//==============================================================================

REveDataCollection::REveDataCollection(const std::string& n, const std::string& t) :
   REveElement(n, t)
{
   std::string lname = n + "Items";
   fItemList = new REveDataItemList(lname.c_str());
   AddElement(fItemList);

   SetupDefaultColorAndTransparency(fgDefaultColor, true, true);
}

void REveDataCollection::AddItem(void *data_ptr, const std::string& /*n*/, const std::string& /*t*/)
{
   auto el = new REveDataItem(data_ptr, GetMainColor());
   fItemList->fItems.emplace_back(el);
}

//------------------------------------------------------------------------------

void REveDataCollection::SetFilterExpr(const TString& filter)
{
   static const REveException eh("REveDataCollection::SetFilterExpr ");

   if (!fItemClass) throw eh + "item class has to be set before the filter expression.";

   fFilterExpr = filter;

   std::stringstream s;
   s << "*((std::function<bool(" << fItemClass->GetName() << "*)>*)" << std::hex << std::showbase << (size_t)&fFilterFoo
     << ") = [](" << fItemClass->GetName() << "* p){" << fItemClass->GetName() << " &i=*p; return ("
     << fFilterExpr.Data() << "); }";

   // printf("%s\n", s.Data());
   try {
      gROOT->ProcessLine(s.str().c_str());
      // AMT I don't know why ApplyFilter call is separated
      ApplyFilter();
   }
   catch (const std::exception &exc)
   {
      std::cerr << "EveDataCollection::SetFilterExpr" << exc.what();
   }
}

void REveDataCollection::ApplyFilter()
{
   Ids_t ids;
   int idx = 0;
   for (auto &ii : fItemList->fItems)
   {
      bool res = fFilterFoo(ii->GetDataPtr());

      // printf("Item:%s -- filter result = %d\n", ii.fItemPtr->GetElementName(), res);

      ii->SetFiltered( ! res );

      ids.push_back(idx++);
   }
   StampObjProps();
   fItemList->StampObjProps();
   fItemList->fHandlerItemsChange( fItemList , ids);
}

//______________________________________________________________________________

void  REveDataCollection::StreamPublicMethods(nlohmann::json &j) const
{
   struct PubMethods
   {
      void FillJSON(TClass* c, nlohmann::json & arr)
      {
         TString  ctor = c->GetName(), dtor = "~";
         {
            int i = ctor.Last(':');
            if (i != kNPOS)
            {
               ctor.Replace(0, i + 1, "");
            }
            dtor += ctor;
         }

         TMethod *meth;
         TIter    next(c->GetListOfMethods());
         while ((meth = (TMethod*) next()))
         {
            // Filter out ctor, dtor, some ROOT stuff.
            {
               TString m(meth->GetName());
               if (m == ctor || m == dtor ||
                   m == "Class" || m == "Class_Name" || m == "Class_Version" || m == "Dictionary" || m == "IsA" ||
                   m == "DeclFileName" || m == "ImplFileName" || m == "DeclFileLine" || m == "ImplFileLine" ||
                   m == "Streamer" || m == "StreamerNVirtual" || m == "ShowMembers" ||
                   m == "CheckTObjectHashConsistency")
               {
                  continue;
               }
            }

            TString     ms;
            TMethodArg *ma;
            TIter       next_ma(meth->GetListOfMethodArgs());
            while ((ma = (TMethodArg*) next_ma()))
            {
               if ( ! ms.IsNull()) ms += ", ";

               ms += ma->GetTypeName();
               ms += " ";
               ms += ma->GetName();
            }
            char* entry = Form("i.%s(%s)",meth->GetName(),ms.Data());
            nlohmann::json jm ;
            jm["f"] = entry;
            jm["r"] = meth->GetReturnTypeName();
            jm["c"] = c->GetName();
            arr.push_back(jm);
         }
         {
            TBaseClass *base;
            TIter       blnext(c->GetListOfBases());
            while ((base = (TBaseClass*) blnext()))
            {
               FillJSON(base->GetClassPointer(), arr);
            }
         }
      }
   };
   j["fPublicFunctions"]  = nlohmann::json::array();
   PubMethods pm;
   pm.FillJSON(fItemClass, j["fPublicFunctions"]);
}

//______________________________________________________________________________

void REveDataCollection::SetMainColor(Color_t newv)
{
   int idx = 0;
   Ids_t ids;
   for (auto & chld : fItemList->fItems)
   {
      chld->SetMainColor(newv);
      ids.push_back(idx);
      idx++;
   }

   REveElement::SetMainColor(newv);
   for (auto & chld : fItemList->fItems)
   {
      chld->SetMainColor(newv);
   }
   fItemList->StampObjProps();
   fItemList->SetMainColor(newv);
   fItemList->fHandlerItemsChange( fItemList , ids);
}

//______________________________________________________________________________

Bool_t REveDataCollection::SetRnrState(Bool_t iRnrSelf)
{
   Bool_t ret = REveElement::SetRnrState(iRnrSelf);
   Ids_t ids;
   for (int i = 0; i < GetNItems(); ++i ) {
      ids.push_back(i);
      fItemList->fItems[i]->SetRnrSelf(fRnrSelf);
   }

   fItemList->StampObjProps();
   fItemList->fHandlerItemsChange( fItemList , ids);

   return ret;
}


//______________________________________________________________________________

Int_t REveDataCollection::WriteCoreJson(nlohmann::json &j, Int_t rnr_offset)
{
   Int_t ret = REveElement::WriteCoreJson(j, rnr_offset);
   j["fFilterExpr"] = fFilterExpr.Data();
   return ret;
}

