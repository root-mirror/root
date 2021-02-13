// @(#)root/core/meta:$Id$
// Author: Axel Naumann 2020-08-25

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TClingMemberIter.h"

#include "cling/Interpreter/Interpreter.h"
#include "cling/Interpreter/LookupHelper.h"

#include "clang/AST/DeclTemplate.h"

ClingMemberIterInternal::DCIter::DCIter(clang::DeclContext *DC, cling::Interpreter *interp) : fInterp(interp)
{
   cling::Interpreter::PushTransactionRAII RAII(fInterp);
   DC->collectAllContexts(fContexts);
   fDeclIter = fContexts[0]->decls_begin();
   // Skip initial empty decl contexts.
   while (IsValid() && fDeclIter == fContexts[fDCIdx]->decls_end()) {
      ++fDCIdx;
      if (fDCIdx < fContexts.size())
         fDeclIter = fContexts[fDCIdx]->decls_begin();
      else
         fDeclIter = fContexts.back()->decls_end();
   }
   AdvanceToFirstValidDecl();
}

bool ClingMemberIterInternal::DCIter::HandleInlineDeclContext()
{
   if (auto *NSD = llvm::dyn_cast<clang::NamespaceDecl>(*fDeclIter)) {
      if (NSD->isInlineNamespace() || NSD->isAnonymousNamespace()) {
         // Collect e.g. internal `__cling_N5xxx' inline namespaces; they will be traversed later
         // Top-most inline namespaces are folded into the iteration:
         fContexts.push_back(NSD);
         return true;
      }
   } else if (auto *ED = llvm::dyn_cast<clang::EnumDecl>(*fDeclIter)) {
      if (!ED->isScoped()) {
         // Inline enums folded into the iteration:
         fContexts.push_back(ED);
         return true;
      }
   } else if (auto *RD = llvm::dyn_cast<clang::RecordDecl>(*fDeclIter)) {
      if (RD->isAnonymousStructOrUnion()) {
         // Anonymous unions are folded into the iteration:
         fContexts.push_back(RD);
         return true;
      }
   }
   return false;
}

bool ClingMemberIterInternal::DCIter::AdvanceToFirstValidDecl()
{
   if (!IsValid())
      return false;

   while (HandleInlineDeclContext())
      if (!IterNext())
         return false;

   return true;
}

bool ClingMemberIterInternal::DCIter::IterNext()
{
   ++fDeclIter;
   while (fDeclIter == fContexts[fDCIdx]->decls_end()) {
      ++fDCIdx;
      if (fDCIdx == fContexts.size())
         return false;
      cling::Interpreter::PushTransactionRAII RAII(fInterp);
      fDeclIter = fContexts[fDCIdx]->decls_begin();
   }
   return true;
}

bool ClingMemberIterInternal::DCIter::Next()
{
   IterNext();
   return AdvanceToFirstValidDecl();
}

ClingMemberIterInternal::UsingDeclIter::UsingDeclIter(const clang::UsingDecl *UD, cling::Interpreter *interp)
   : fInterp(interp)
{
   cling::Interpreter::PushTransactionRAII RAII(interp);
   fUsingIterStack.push({UD});
}

bool ClingMemberIterInternal::UsingDeclIter::Next()
{
   ++Iter();
   while (true) {
      if (Iter() == End()) {
         // End of this UD's loop; continue iteration with parent.
         fUsingIterStack.pop();
         if (fUsingIterStack.empty())
            return false;
         ++Iter(); // parent was the UsingDecl we just finished, move to next.
         continue;
      }
      if (auto *UD = llvm::dyn_cast<clang::UsingDecl>(Iter()->getTargetDecl())) {
         if (UD->shadow_size()) {
            cling::Interpreter::PushTransactionRAII RAII(fInterp);
            fUsingIterStack.push({UD});
            // Continue with child.
         }
      } else {
         break;
      }
   };
   return true;
}

bool TClingMemberIter::Advance()
{
   fTemplateSpec = nullptr;
   do {
      const clang::Decl *D = Get();
      if (auto *UD = llvm::dyn_cast<clang::UsingDecl>(D)) {
         if (UD->shadow_size()) {
            assert(!fUsingDeclIter.IsValid() && "Expected UsingDecl to be already handled by UsingDeclIter!");
            AdvanceUnfiltered();
            fUsingDeclIter = ClingMemberIterInternal::UsingDeclIter(UD, fInterp);
            continue;
         }
      }
      if (auto *USD = llvm::dyn_cast<clang::UsingShadowDecl>(D)) {
         if (!ShouldSkip(USD))
            return true;
      } else if (auto *RTD = llvm::dyn_cast<clang::RedeclarableTemplateDecl>(D)) {
         if (const clang::Decl *DInst = InstantiateTemplateWithDefaults(RTD)) {
            fTemplateSpec = DInst;
            return true;
         }
      } else if (!ShouldSkip(D))
         return true;

      // Not interested in this one, continue.
      if (!AdvanceUnfiltered())
         break;
   } while (true);
   return false;
}
