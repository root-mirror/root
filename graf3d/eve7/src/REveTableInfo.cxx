// @(#)root/eve7:$Id$
// Authors: Matevz Tadel & Alja Mrak-Tadel: 2018

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TClass.h"
#include "TBaseClass.h"
#include "TROOT.h"
#include "TInterpreter.h"
#include "TMethod.h"
#include "TMethodArg.h"

#include <ROOT/REveTableInfo.hxx>
#include <ROOT/REveManager.hxx>

#include <sstream>

using namespace ROOT::Experimental;

REveTableViewInfo::REveTableViewInfo(const std::string &name, const std::string &title)
   : REveElement(name, title), fConfigChanged(false)
{
}

void REveTableViewInfo::SetDisplayedCollection(ElementId_t collectionId)
{
   fDisplayedCollection = collectionId;

   fConfigChanged = true;
   for (auto &it : fDelegates)
      it();

   fConfigChanged = false;
   StampObjProps();
}

void REveTableViewInfo::AddNewColumnToCurrentCollection(const std::string &expr, const std::string &title, int prec)
{
   if (!fDisplayedCollection)
      return;

   REveDataCollection *col = dynamic_cast<REveDataCollection *>(gEve->FindElementById(fDisplayedCollection));
   if (!col) {
      printf("REveTableViewInfo::AddNewColumnToCurrentCollection error: collection not found\n");
      return;
   }

   const char *rtyp = "void";
   auto icls = col->GetItemClass();
   std::function<void(void *)> fooptr;
   std::stringstream s;
   s << "*((std::function<" << rtyp << "(" << icls->GetName() << "*)>*)" << std::hex << std::showbase
     << (size_t)(&fooptr) << ") = [](" << icls->GetName() << "* p){" << icls->GetName() << " &i=*p; return (" << expr
     << "); }";

   int err;
   gROOT->ProcessLine(s.str().c_str(), &err);
   if (err != TInterpreter::kNoError) {
      std::cout << "REveTableViewInfo::AddNewColumnToCurrentCollection failed." << std::endl;
      return;
   }

   fConfigChanged = true;
   table(col->GetItemClass()->GetName()).column(title, prec, expr);

   for (auto &it : fDelegates)
      it();

   fConfigChanged = false;

   StampObjProps();
}

////////////////////////////////////////////////////////////////////////////////
/// Find column definitions for given class name.
//  Look for definition also in base classes
REveTableHandle::Entries_t &REveTableViewInfo::RefTableEntries(std::string cname)
{
   struct TableDictHelper {
      void fillPublicMethods(REveTableHandle::Entries_t &entries, TClass *c)
      {
         static size_t maxEnt = 3;
         TMethod *meth;
         TIter next(c->GetListOfAllPublicMethods());
         while ((meth = (TMethod *)next()) && entries.size() < maxEnt) {
            // take only methods without arguments
            if (!meth->GetListOfMethodArgs()->First()) {
               std::string mn = meth->GetName();
               std::string exp = "i." + mn + "()";

               TDataType *dt = gROOT->GetType(meth->GetReturnTypeName());
               if (dt) {
                  int t = dt->GetType();
                  if (t == EDataType::kInt_t || t == EDataType::kUInt_t || t == EDataType::kLong_t ||
                      t == EDataType::kULong_t || t == EDataType::kLong64_t || t == EDataType::kULong64_t ||
                      t == EDataType::kBool_t)
                     entries.push_back(REveTableEntry(mn, 0, exp));
                  else if (t == EDataType::kFloat_t || t == EDataType::kDouble_t || t == EDataType::kDouble32_t)
                     entries.push_back(REveTableEntry(mn, 3, exp));
               }
            }
         }

         // look in the base classes
         TBaseClass *base;
         TIter blnext(c->GetListOfBases());
         while ((base = (TBaseClass *)blnext())) {
            fillPublicMethods(entries, base->GetClassPointer());
         }
      }

      TClass *searchMatchInBaseClasses(TClass *c, REveTableHandle::Specs_t &specs)
      {
         TBaseClass *base;
         TIter blnext(c->GetListOfBases());
         while ((base = (TBaseClass *)blnext())) {
            auto bs = specs.find(base->GetName());
            if (bs != specs.end()) {
               return base->GetClassPointer();
            }
            return searchMatchInBaseClasses(base->GetClassPointer(), specs);
         }
         return nullptr;
      }
   };

   TableDictHelper helper;
   auto search = fSpecs.find(cname);
   if (search != fSpecs.end()) {
      return search->second;
   } else {
      TClass *b = helper.searchMatchInBaseClasses(TClass::GetClass(cname.c_str()), fSpecs);
      if (b) {
         return fSpecs[b->GetName()];
      }
   }

   // create new entry if not existing
   helper.fillPublicMethods(fSpecs[cname], TClass::GetClass(cname.c_str()));
   return fSpecs[cname];
}

////////////////////////////////////////////////////////////////////////////////
/// Fill core part of JSON representation.

Int_t REveTableViewInfo::WriteCoreJson(nlohmann::json &j, Int_t rnr_offset)
{
   auto ret = REveElement::WriteCoreJson(j, rnr_offset);
   j["fDisplayedCollection"] = fDisplayedCollection;
   return ret;
}
