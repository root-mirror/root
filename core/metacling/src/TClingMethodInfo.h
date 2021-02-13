// @(#)root/core/meta:$Id$
// Author: Paul Russo   30/07/2012

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TClingMethodInfo
#define ROOT_TClingMethodInfo

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TClingMethodInfo                                                     //
//                                                                      //
// Emulation of the CINT MethodInfo class.                              //
//                                                                      //
// The CINT C++ interpreter provides an interface to metadata about     //
// a function through the MethodInfo class.  This class provides the    //
// same functionality, using an interface as close as possible to       //
// MethodInfo but the typedef metadata comes from the Clang C++         //
// compiler, not CINT.                                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TClingDeclInfo.h"
#include "TClingMemberIter.h"
#include "TDictionary.h"
#include "TString.h"

#include "llvm/ADT/SmallVector.h"

#include <string>

namespace clang {
   class Decl;
   class FunctionDecl;
}

namespace cling {
   class Interpreter;
}

namespace ROOT {
   namespace TMetaUtils {
      class TNormalizedCtxt;
   }
}

class TClingClassInfo;
class TClingTypeInfo;

/// Iterate over FunctionDecl and UsingShadowDecls of FunctionDecl, within a scope,
/// recursing through "transparent" scopes (see DCIter::HandleInlineDeclContext()).
class TClingCXXRecMethIter final: public TClingMemberIter {

class SpecFuncIter {
   using Vec_t = llvm::SmallVector<clang::CXXMethodDecl *,4>;
   Vec_t fDefDataSpecFuns; // Special functions materialized from DefinitionData.
   size_t fIDefDataSpecFuns = 0; // Current element in fDefDataSpecFuns.

public:
   SpecFuncIter() = default;
   SpecFuncIter(cling::Interpreter *interp, clang::DeclContext *DC,
                llvm::SmallVectorImpl<clang::CXXMethodDecl*> &&specFuncs);

   bool IsValid() const { return fIDefDataSpecFuns < fDefDataSpecFuns.size(); }

   bool Next() {
      ++fIDefDataSpecFuns;
      return IsValid();
   }

   clang::FunctionDecl *operator->() const {
      return operator*();
   }

   clang::FunctionDecl *operator*() const {
      return IsValid() ? fDefDataSpecFuns[fIDefDataSpecFuns] : nullptr;
   }
};

   SpecFuncIter fSpecFuncIter;

protected:
   const clang::Decl *
   InstantiateTemplateWithDefaults(const clang::RedeclarableTemplateDecl *TD) const final;

   bool ShouldSkip(const clang::Decl* FD) const final;
   bool ShouldSkip(const clang::UsingShadowDecl* USD) const final;

   clang::Decl *AdvanceUnfiltered() final {
      if (fSpecFuncIter.IsValid() && fSpecFuncIter.Next())
         return *fSpecFuncIter;
      return TClingMemberIter::AdvanceUnfiltered();
   }

public:
   TClingCXXRecMethIter() = default;
   TClingCXXRecMethIter(cling::Interpreter *interp, clang::DeclContext *DC,
                        llvm::SmallVectorImpl<clang::CXXMethodDecl*> &&specFuncs):
      TClingMemberIter(interp, DC), fSpecFuncIter(interp, DC, std::move(specFuncs)) {}

   const clang::Decl *Get() const final {
      if (!IsValid())
         return nullptr;

      if (fSpecFuncIter.IsValid())
         return *fSpecFuncIter;

      return TClingMemberIter::Get();
   }

   bool IsValid() const final {
      return fSpecFuncIter.IsValid() || TClingMemberIter::IsValid();
   }

};


class TClingMethodInfo final : public TClingDeclInfo {
private:
   cling::Interpreter                          *fInterp; // Cling interpreter, we do *not* own.
   bool                                         fFirstTime; // Flag for first time incrementing iterator, cint semantics are weird.
   TClingCXXRecMethIter                         fIter; // Our iterator.
   std::string                                  fTitle; // The meta info for the method.

   const clang::Decl* GetDeclSlow() const;

public:
   explicit TClingMethodInfo(cling::Interpreter *interp)
      : TClingDeclInfo(nullptr), fInterp(interp), fFirstTime(true), fTitle("") {}

   // Takes concrete decl and disables the iterator.
   TClingMethodInfo(cling::Interpreter *, const clang::Decl *);
   TClingMethodInfo(cling::Interpreter *, TClingClassInfo *);

   const clang::FunctionDecl                   *GetAsFunctionDecl() const;
   const clang::UsingShadowDecl                *GetAsUsingShadowDecl() const;

   /// Get the FunctionDecl, or if this represents a UsingShadowDecl, the underlying target FunctionDecl.
   const clang::FunctionDecl                   *GetTargetFunctionDecl() const;

   TDictionary::DeclId_t                        GetDeclId() const;
   cling::Interpreter                          *GetInterpreter() const { return fInterp; }
   void                                         CreateSignature(TString &signature) const;
   void                                         Init(const clang::FunctionDecl *);
   void                                        *InterfaceMethod(const ROOT::TMetaUtils::TNormalizedCtxt &normCtxt) const;

   const clang::Decl *GetDecl() const override {
     if (const clang::Decl* SingleDecl = TClingDeclInfo::GetDecl())
       return SingleDecl;
     return GetDeclSlow();
   }
   int                                          NArg() const;
   int                                          NDefaultArg() const;
   int                                          Next();
   long                                         Property() const;
   long                                         ExtraProperty() const;
   TClingTypeInfo                              *Type() const;
   std::string                                  GetMangledName() const;
   const char                                  *GetPrototype();
   const char                                  *Name() const override;
   const char                                  *TypeName() const;
   const char                                  *Title();
};

#endif // ROOT_TClingMethodInfo
