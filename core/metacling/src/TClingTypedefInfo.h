// @(#)root/core/meta:$Id$
// Author: Paul Russo   30/07/2012

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TClingTypedefInfo
#define ROOT_TClingTypedefInfo

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TClingTypedefInfo                                                    //
//                                                                      //
// Emulation of the CINT TypedefInfo class.                             //
//                                                                      //
// The CINT C++ interpreter provides an interface to metadata about     //
// a typedef through the TypedefInfo class.  This class provides the    //
// same functionality, using an interface as close as possible to       //
// TypedefInfo but the typedef metadata comes from the Clang C++        //
// compiler, not CINT.                                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TClingDeclInfo.h"

#include "cling/Interpreter/Interpreter.h"
#include "clang/AST/ASTContext.h"
#include "clang/AST/Decl.h"
#include "clang/Frontend/CompilerInstance.h"

#include <vector>
#include <string>

namespace ROOT {
   namespace TMetaUtils {
      class TNormalizedCtxt;
   }
}

class TClingTypedefInfo final : public TClingDeclInfo {

private:

   cling::Interpreter  *fInterp; // Cling interpreter, we do *not* own.
   bool                 fFirstTime; // We need to skip the first increment to support the cint Next() semantics.
   bool                 fDescend; // Flag for signaling the need to descend on this advancement.
   clang::DeclContext::decl_iterator fIter; // Current decl in scope.
   std::vector<clang::DeclContext::decl_iterator> fIterStack; // Recursion stack for traversing nested scopes.
   std::string          fTitle; // The meta info for the typedef.

public:

   explicit TClingTypedefInfo(cling::Interpreter *interp)
      : TClingDeclInfo(nullptr), fInterp(interp), fFirstTime(true), fDescend(false), fTitle("")
   {
      const clang::TranslationUnitDecl *TU = fInterp->getCI()->getASTContext().getTranslationUnitDecl();
      const clang::DeclContext *DC = llvm::cast<clang::DeclContext>(TU);
      cling::Interpreter::PushTransactionRAII RAII(fInterp);
      fIter = DC->decls_begin();
   }

   explicit TClingTypedefInfo(cling::Interpreter *, const char *);

   explicit TClingTypedefInfo(cling::Interpreter *, const clang::TypedefNameDecl *);

   void                 Init(const char *name);
   int                  InternalNext();
   int                  Next();
   long                 Property() const;
   int                  Size() const;
   const char          *TrueName(const ROOT::TMetaUtils::TNormalizedCtxt &normCtxt) const;
   const char          *Name() const override;
   const char          *Title();

};

#endif // ROOT_TClingTypedefInfo
