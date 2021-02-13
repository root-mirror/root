// Author:  Sergey Linev, GSI  29/06/2017

/*************************************************************************
 * Copyright (C) 1995-2021, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TWebMenuItem.h"

#include "TClass.h"
#include "TList.h"
#include "TMethod.h"
#include "TMethodArg.h"
#include "TMethodCall.h"

void TWebMenuItems::PopulateObjectMenu(void *obj, TClass *cl)
{
   fItems.clear();

   TList *lst = new TList;
   cl->GetMenuItems(lst);

   TIter iter(lst);
   TMethod *m = nullptr;

   Bool_t has_editor = kFALSE;

   TClass *last_class = nullptr;

   while ((m = (TMethod *)iter()) != nullptr) {

      Bool_t is_editor = kFALSE;

      if (strcmp(m->GetClass()->GetName(), "TH1") == 0) {
         if (strcmp(m->GetName(), "SetHighlight") == 0) continue;
         if (strcmp(m->GetName(), "DrawPanel") == 0) is_editor = kTRUE;
      } else if (strcmp(m->GetClass()->GetName(), "TGraph") == 0) {
         if (strcmp(m->GetName(), "SetHighlight") == 0) continue;
         if (strcmp(m->GetName(), "DrawPanel") == 0) is_editor = kTRUE;
      } else if (strcmp(m->GetClass()->GetName(), "TAttFill") == 0) {
         if (strcmp(m->GetName(), "SetFillAttributes") == 0) is_editor = kTRUE;
      } else if (strcmp(m->GetClass()->GetName(), "TAttLine") == 0) {
         if (strcmp(m->GetName(), "SetLineAttributes") == 0) is_editor = kTRUE;
      } else if (strcmp(m->GetClass()->GetName(), "TAttMarker") == 0) {
         if (strcmp(m->GetName(), "SetMarkerAttributes") == 0) is_editor = kTRUE;
      } else if (strcmp(m->GetClass()->GetName(), "TAttText") == 0) {
         if (strcmp(m->GetName(), "SetTextAttributes") == 0) is_editor = kTRUE;
      }

      if (is_editor) {
         if (!has_editor) {
            AddMenuItem("Editor", "Attributes editor", "Show:Editor", last_class ? last_class : m->GetClass());
            has_editor = kTRUE;
         }
         continue;
      }

      last_class = m->GetClass();

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

            TMethodCall *call = new TMethodCall(cl, getter, "");

            if (call->ReturnType() == TMethodCall::kLong) {
               Long_t l(0);
               call->Execute(obj, l);

               AddChkMenuItem(m->GetName(), m->GetTitle(), l != 0, Form("%s(%s)", m->GetName(), (l != 0) ? "0" : "1"), m->GetClass());

            } else {
               // Error("CheckModifiedFlag", "Cannot get toggle value with getter %s", getter.Data());
            }

            delete call;
         }
      } else {
         TList *args = m->GetListOfMethodArgs();

         if (!args || (args->GetSize() == 0)) {
            AddMenuItem(m->GetName(), m->GetTitle(), Form("%s()", m->GetName()), m->GetClass());
         } else {
            TWebArgsMenuItem *item = new TWebArgsMenuItem(m->GetName(), m->GetTitle());
            item->SetExec(Form("%s()", m->GetName()));
            if (m->GetClass()) item->SetClassName(m->GetClass()->GetName());

            TIter args_iter(args);
            TMethodArg *arg = nullptr;

            while ((arg = dynamic_cast<TMethodArg *>(args_iter())) != nullptr) {
               const char *dflt = arg->GetDefault();
               if (!dflt) dflt = "";
               item->GetArgs().emplace_back(arg->GetName(), arg->GetTitle(), arg->GetFullTypeName(), dflt);
            }

            Add(item);
         }
      }
   }

   delete lst;
}
