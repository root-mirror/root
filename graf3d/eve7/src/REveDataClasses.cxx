// @(#)root/eve7:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007, 2018

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include <ROOT/REveDataClasses.hxx>

#include "TROOT.h"
#include "TMethod.h"
#include "TMethodArg.h"

#include "json.hpp"
#include <sstream>


using namespace ROOT::Experimental;
namespace REX = ROOT::Experimental;


//==============================================================================
// REveDataCollection
//==============================================================================

REveDataCollection::REveDataCollection(const char* n, const char* t) :
   REveElementList(n, t)
{
   fChildClass = REveDataItem::Class();
}

void REveDataCollection::AddItem(void *data_ptr, const char *n, const char *t)
{
   auto el = new REveDataItem(n, t);
   AddElement(el);
   fItems.emplace_back(data_ptr, el);
}

//------------------------------------------------------------------------------

void REveDataCollection::SetFilterExpr(const TString& filter)
{
   static const REveException eh("REveDataCollection::SetFilterExpr ");

   if ( ! fItemClass) throw eh + "item class has to be set before the filter expression.";

   fFilterExpr = filter;

   std::stringstream s;
   s << "*((std::function<bool(" << fItemClass->GetName() << "*)>*)" << std::hex << std::showbase
     << (size_t)&fFilterFoo << ") = [](" << fItemClass->GetName() << "* p){" << fItemClass->GetName()
     << " &i=*p; return (" << fFilterExpr.Data() << "); }";

   // printf("%s\n", s.Data());
   try {
      gROOT->ProcessLine(s.str().c_str());
   }
   catch (const std::exception &exc)
   {
      std::cerr << "EveDataCollection::SetFilterExpr" << exc.what();
   }
}

void REveDataCollection::ApplyFilter()
{
   for (auto &ii : fItems)
   {
      bool res = fFilterFoo(ii.fDataPtr);

      // printf("Item:%s -- filter result = %d\n", ii.fItemPtr->GetElementName(), res);

      ii.fItemPtr->SetFiltered( ! res );
   }
}


Int_t REveDataCollection::WriteCoreJson(nlohmann::json &j, Int_t rnr_offset)
{
   Int_t ret = REveElement::WriteCoreJson(j, rnr_offset);
   j["fFilterExpr"] = fFilterExpr.Data();
   j["publicFunction"]  = nlohmann::json::array();

   TIter x( fItemClass->GetListOfAllPublicMethods());
   while (TObject *obj = x()) {
      // printf("func %s \n", obj->GetName());
      nlohmann::json m;


      TMethod* method = dynamic_cast<TMethod*>(obj);
      m["name"] = method->GetPrototype();
      j["publicFunction"].push_back(m);
   }

   return ret;
}
//==============================================================================
// REveDataItem
//==============================================================================

REveDataItem::REveDataItem(const char* n, const char* t) :
   REveElementList(n, t)
{
}

Int_t REveDataItem::WriteCoreJson(nlohmann::json &j, Int_t rnr_offset)
{
   Int_t ret = REveElement::WriteCoreJson(j, rnr_offset);
   j["fFiltered"] = fFiltered;
   return ret;
}

//==============================================================================
// REveDataTable
//==============================================================================

REveDataTable::REveDataTable(const char* n, const char* t) :
   REveElementList(n, t)
{
   fChildClass = REveDataColumn::Class();
}

void REveDataTable::PrintTable()
{
   Int_t Nit = fCollection->GetNItems();

   for (Int_t i = 0; i< Nit; ++i)
   {
      void         *data = fCollection->GetDataPtr(i);
      REveDataItem *item = fCollection->GetDataItem(i);

      printf("| %-20s |", item->GetElementName());

      for (auto & chld : fChildren)
      {
         auto clmn = dynamic_cast<REveDataColumn*>(chld);

         printf(" %10s |", clmn->EvalExpr(data).c_str());
      }
      printf("\n");
   }
}

Int_t REveDataTable::WriteCoreJson(nlohmann::json &j, Int_t rnr_offset)
{
   int ret = REveElement::WriteCoreJson(j, rnr_offset);
   Int_t Nit = fCollection->GetNItems();

   nlohmann::json jarr = nlohmann::json::array();

   for (Int_t i = 0; i< Nit; ++i)
   {
      void         *data = fCollection->GetDataPtr(i);
      nlohmann::json row;
      for (auto & chld : fChildren)
      {
         auto clmn = dynamic_cast<REveDataColumn*>(chld);
         row[chld->GetElementName()] = clmn->EvalExpr(data);
         // printf(" %10s |", clmn->EvalExpr(data).c_str());

      }
      jarr.push_back(row);
   }
   j["body"] = jarr;
   j["fCollectionId"] = fCollection->GetElementId();
   return ret;
}

void REveDataTable::AddNewColumn(const char* expr, const char* title, int prec)
{
   auto c = new REX::REveDataColumn(title);
   AddElement(c);

   c->SetExpressionAndType(expr, REX::REveDataColumn::FT_Double);
   c->SetPrecision(prec);

   StampObjProps();
}

//==============================================================================
// REveDataColumn
//==============================================================================

REveDataColumn::REveDataColumn(const char* n, const char* t) :
   REveElementList(n, t)
{
}

void REveDataColumn::SetExpressionAndType(const TString& expr, FieldType_e type)
{
   auto table = dynamic_cast<REveDataTable*>(fMother);
   auto coll = table->GetCollection();
   auto icls = coll->GetItemClass();

   fExpression = expr;
   fType       = type;

   const char *rtyp{nullptr};
   const void *fooptr{nullptr};

   switch (fType)
   {
      case FT_Double: rtyp = "double";      fooptr = &fDoubleFoo; break;
      case FT_Bool:   rtyp = "bool";        fooptr = &fBoolFoo;   break;
      case FT_String: rtyp = "std::string"; fooptr = &fStringFoo; break;
   }

   std::stringstream s;
   s << "*((std::function<" << rtyp << "(" << icls->GetName() << "*)>*)" << std::hex
     << std::showbase << (size_t)fooptr << ") = [](" << icls->GetName() << "* p){"
     << icls->GetName() << " &i=*p; return (" << fExpression.Data() << "); }";

   // printf("%s\n", s.Data());
   try {
      gROOT->ProcessLine(s.str().c_str());
   }
   catch (const std::exception &exc)
   {
      std::cerr << "REveDataColumn::SetExpressionAndType" << exc.what();
   }
}

void REveDataColumn::SetPrecision(Int_t prec)
{
   fPrecision = prec;
}

std::string REveDataColumn::EvalExpr(void *iptr)
{
   switch (fType)
   {
      case FT_Double:
      {
         TString ostr;
         ostr.Form("%.*f", fPrecision, fDoubleFoo(iptr));
         return ostr.Data();
      }
      case FT_Bool:
      {
         return fBoolFoo(iptr) ? fTrue : fFalse;
      }
      case FT_String:
      {
         return fStringFoo(iptr);
      }
   }
   return "XYZZ";
}
