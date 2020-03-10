/*************************************************************************
 * Copyright (C) 1995-2017, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "ROOT/RMenuItem.hxx"

#include "TROOT.h"
#include "TString.h"
#include "TClass.h"
#include "TList.h"
#include "TMethod.h"
#include "TMethodArg.h"
#include "TMethodCall.h"
#include "TBufferJSON.h"

//////////////////////////////////////////////////////////
/// Returns drawable id coded, fullid can include # with extra specifier

std::string ROOT::Experimental::RMenuItems::GetDrawableId()
{
   auto p = fId.find("#");
   if (p == std::string::npos)
      return fId;
   return fId.substr(0, p);
}


//////////////////////////////////////////////////////////
/// Returns specifier string after # in fullid

std::string ROOT::Experimental::RMenuItems::GetSpecifier()
{
   auto p = fId.find("#");
   if (p == std::string::npos)
      return "";
   return fId.substr(p+1);
}


//////////////////////////////////////////////////////////
/// Fill menu for provided object, using *MENU* as indicator in method comments

void ROOT::Experimental::RMenuItems::PopulateObjectMenu(void *obj, TClass *cl)
{
   fItems.clear();

   TList lst;
   cl->GetMenuItems(&lst);

   TIter iter(&lst);
   TMethod *m = nullptr;

   while ((m = (TMethod *)iter()) != nullptr) {

      if (m->IsMenuItem() == kMenuToggle) {
         TString getter;
         if (m->Getter() && strlen(m->Getter()) > 0) {
            getter = m->Getter();
         } else if (strncmp(m->GetName(), "Set", 3) == 0) {
            getter = TString(m->GetName())(3, strlen(m->GetName()) - 3);
            if (cl->GetMethodAllAny(TString("Has") + getter))
               getter = TString("Has") + getter;
            else if (cl->GetMethodAllAny(TString("Get") + getter))
               getter = TString("Get") + getter;
            else if (cl->GetMethodAllAny(TString("Is") + getter))
               getter = TString("Is") + getter;
            else
               getter = "";
         }

         if ((getter.Length() > 0) && cl->GetMethodAllAny(getter)) {
            // execute getter method to get current state of toggle item

            auto call = std::make_unique<TMethodCall>(cl, getter, "");

            if (call->ReturnType() == TMethodCall::kLong) {
               Long_t l(0);
               call->Execute(obj, l);

               AddChkMenuItem(m->GetName(), m->GetTitle(), l != 0, Form("%s(%s)", m->GetName(), (l != 0) ? "0" : "1"));

            } else {
               // Error("CheckModifiedFlag", "Cannot get toggle value with getter %s", getter.Data());
            }
         }
      } else {
         TList *args = m->GetListOfMethodArgs();

         if (!args || (args->GetSize() == 0)) {
            AddMenuItem(m->GetName(), m->GetTitle(), Form("%s()", m->GetName()));
         } else {
            auto item = std::make_unique<Detail::RArgsMenuItem>(m->GetName(), m->GetTitle());
            item->SetExec(Form("%s()", m->GetName()));

            TIter args_iter(args);
            TMethodArg *arg = nullptr;

            while ((arg = dynamic_cast<TMethodArg *>(args_iter())) != nullptr) {
               Detail::RMenuArgument menu_arg(arg->GetName(), arg->GetTitle(), arg->GetFullTypeName());
               if (arg->GetDefault()) menu_arg.SetDefault(arg->GetDefault());
               item->AddArg(menu_arg);
            }

            Add(std::move(item));
         }
      }
   }
}
