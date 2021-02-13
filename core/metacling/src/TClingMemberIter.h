// @(#)root/core/meta:$Id$
// Author: Axel Naumann 2020-08-25

/*************************************************************************
 * Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TClingMemberIter
#define ROOT_TClingMemberIter

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TClingMemberIter                                                     //
//                                                                      //
// Iteration through function or data members,                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "clang/AST/DeclBase.h"
#include "clang/AST/DeclCXX.h"

#include "llvm/ADT/SmallVector.h"

#include <memory>
#include <stack>

namespace clang {
class Sema;
class RedeclarableDemplateDecl;
} // namespace clang

namespace cling {
class Interpreter;
class LookupHelper;
} // namespace cling

namespace ClingMemberIterInternal {
/// Iteration over collected DeclContexts
class DCIter {
   cling::Interpreter *fInterp = nullptr;
   llvm::SmallVector<clang::DeclContext *, 2> fContexts; // Set of DeclContext to iterate over.
   size_t fDCIdx = 0;                                    // Index of current element in fContexts.
   clang::DeclContext::decl_iterator fDeclIter;

   /// Do not return inline namespaces etc.
   bool AdvanceToFirstValidDecl();

   /// Return true if fDeclIter has triggered an inline or anonymous namespace / ...
   /// to be appended to fContexts.
   bool HandleInlineDeclContext();

   /// Increment the iterator, possibly moving to the next context.
   bool IterNext();

public:
   DCIter() = default;
   DCIter(clang::DeclContext *DC, cling::Interpreter *interp);

   bool Next();

   clang::Decl *operator->() const { return operator*(); }
   clang::Decl *operator*() const { return fInterp ? *fDeclIter : nullptr; }

   bool IsValid() const { return fInterp && fDCIdx < fContexts.size(); }
};

class UsingDeclIter {
   struct UsingDeclFrame {
      clang::UsingDecl::shadow_iterator fShadowedIter; // Shadowed decl iterator of fCurrentUsingDecl.
      clang::UsingDecl::shadow_iterator fShadowedEnd;  // Shadowed decl iterator of fCurrentUsingDecl.

      UsingDeclFrame(const clang::UsingDecl *UD) : fShadowedIter(UD->shadow_begin()), fShadowedEnd(UD->shadow_end()) {}
   };
   /// Handle nested UsingDecls:
   /// ```
   ///    struct A { void f(); void f(int); };
   ///    struct B: A { using A::f; void f(float); };
   ///    struct C: B { using B::f; void f(char); };
   /// ```
   /// will cause the UsingDecl iteration to be stacked,
   /// to recursively resolve all used decls for `f` in `C`.
   std::stack<UsingDeclFrame> fUsingIterStack;
   cling::Interpreter *fInterp = nullptr;

   clang::UsingDecl::shadow_iterator &Iter() { return fUsingIterStack.top().fShadowedIter; }
   const clang::UsingDecl::shadow_iterator &Iter() const { return fUsingIterStack.top().fShadowedIter; }
   const clang::UsingDecl::shadow_iterator &End() const { return fUsingIterStack.top().fShadowedEnd; }

public:
   UsingDeclIter() = default;
   UsingDeclIter(const clang::UsingDecl *UD, cling::Interpreter *interp);

   bool Next();

   clang::Decl *operator->() const { return operator*(); }

   clang::Decl *operator*() const { return IsValid() ? *Iter() : nullptr; }

   bool IsValid() const { return !fUsingIterStack.empty(); }
};

}; // namespace ClingMemberIterInternal

/// Iterate over all DeclT-s (or UsingShadowDecl-s pointing to DeclT-s) of a decl
/// context, skipping those for which DerivedT::ShouldSkip(const Decl*) returns `true`
/// when invoked with the Decl.
class TClingMemberIter {
private:
   cling::Interpreter *fInterp = nullptr;
   const clang::Decl *fTemplateSpec =
      nullptr; // an all-default-template-args member corresponding to the current iteration Decl.
   ClingMemberIterInternal::DCIter fDCIter;               // DeclContext iterator.
   ClingMemberIterInternal::UsingDeclIter fUsingDeclIter; // Iterating the shadowed decls of a using decl.

   const clang::Decl *GetDeclSlow() const;

   bool Advance();

protected:
   virtual clang::Decl *AdvanceUnfiltered()
   {
      if (fUsingDeclIter.IsValid() && fUsingDeclIter.Next())
         return *fUsingDeclIter;
      if (fDCIter.IsValid() && fDCIter.Next())
         return *fDCIter;
      return nullptr;
   }

   virtual const clang::Decl *InstantiateTemplateWithDefaults(const clang::RedeclarableTemplateDecl *TD) const
   {
      // TODO: add support for variable templates.
      return nullptr;
   }

   virtual bool ShouldSkip(const clang::Decl *D) const = 0;
   virtual bool ShouldSkip(const clang::UsingShadowDecl *USD) const { return ShouldSkip((clang::Decl *)USD); };

public:
   TClingMemberIter() = default;

   TClingMemberIter(cling::Interpreter *interp, clang::DeclContext *DC) : fInterp(interp), fDCIter(DC, interp) {}

   void Init()
   {
      // This function needs the derived class to be constructed (vtable and such),
      // and thus cannot be embedded in the constructor.
      if (IsValid())
         Advance();
   }

   /// Advance to next non-skipped; return false if no next decl exists.
   bool Next()
   {
      if (AdvanceUnfiltered())
         return Advance();
      return false;
   }

   virtual const clang::Decl *Get() const
   {
      if (fTemplateSpec)
         return fTemplateSpec;
      if (fUsingDeclIter.IsValid())
         return *fUsingDeclIter;
      if (fDCIter.IsValid())
         return *fDCIter;
      return nullptr;
   }

   const clang::Decl *operator->() const { return Get(); }

   const clang::Decl *operator*() const { return Get(); }

   virtual bool IsValid() const { return fInterp && (fDCIter.IsValid() || fUsingDeclIter.IsValid()); }

   cling::Interpreter *GetInterpreter() const { return fInterp; }
};

#endif // ROOT_TClingMethodInfo
