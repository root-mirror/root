// @(#)root/core/meta:$Id$
// Author: Paul Russo   30/07/2012

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TClingMethodArgInfo
Emulation of the CINT MethodInfo class.

The CINT C++ interpreter provides an interface to metadata about
the arguments to a function through the MethodArgInfo class.  This
class provides the same functionality, using an interface as close
as possible to MethodArgInfo but the typedef metadata comes from
the Clang C++ compiler, not CINT.
*/

#include "TClingMethodArgInfo.h"

#include "TDictionary.h"
#include "TClingMethodInfo.h"
#include "TClingTypeInfo.h"
#include "ThreadLocalStorage.h"

#include "cling/Interpreter/Interpreter.h"
#include "clang/AST/ASTContext.h"
#include "clang/AST/Decl.h"
#include "clang/AST/Expr.h"
#include "clang/AST/ExprCXX.h"
#include "clang/AST/PrettyPrinter.h"
#include "clang/AST/Type.h"
#include "clang/Sema/Sema.h"

#include "llvm/Support/Casting.h"
#include "llvm/Support/raw_ostream.h"

#include <string>

TClingMethodArgInfo::TClingMethodArgInfo(cling::Interpreter *interp, const TClingMethodInfo* mi) : TClingDeclInfo(mi->GetTargetFunctionDecl()), fInterp(interp), fIdx(-1) {}

bool TClingMethodArgInfo::IsValid() const
{
   // Calling the base class implementation is unsafe because we override the
   // GetDecl which it forwads to. That implementation depends on fIdx which is
   // used to determine validity.
   if (!fDecl)
      return false;

   auto FD = llvm::cast_or_null<clang::FunctionDecl>(TClingDeclInfo::GetDecl());
   int numParams = static_cast<int>(FD->getNumParams());
   return (fIdx > -1) && (fIdx < numParams);
}

int TClingMethodArgInfo::Next()
{
   ++fIdx;
   fNameCache.clear(); // invalidate the cache.
   return IsValid();
}

long TClingMethodArgInfo::Property() const
{
   if (!IsValid()) {
      return 0L;
   }
   long property = 0L;
   const clang::ParmVarDecl *pvd = GetDecl();
   if (pvd->hasDefaultArg() || pvd->hasInheritedDefaultArg()) {
      property |= kIsDefault;
   }
   clang::QualType qt = pvd->getOriginalType().getCanonicalType();

   return TClingDeclInfo::Property(property, qt);
}

const char *TClingMethodArgInfo::DefaultValue() const
{
   if (!IsValid()) {
      return 0;
   }
   const clang::ParmVarDecl *pvd = GetDecl();
   // Instantiate default arg if needed
   if (pvd->hasUninstantiatedDefaultArg()) {
      // Could deserialize / create instantiated decls.
      cling::Interpreter::PushTransactionRAII RAII(fInterp);
      auto fd = llvm::cast_or_null<clang::FunctionDecl>(TClingDeclInfo::GetDecl());
      fInterp->getSema().BuildCXXDefaultArgExpr(clang::SourceLocation(),
                                                const_cast<clang::FunctionDecl*>(fd),
                                                const_cast<clang::ParmVarDecl*>(pvd));
   }
   const clang::Expr *expr = 0;
   if (pvd->hasUninstantiatedDefaultArg()) {
      // We tried to instantiate it above; if we fail, use the uninstantiated one.
      expr = pvd->getUninstantiatedDefaultArg();
   } else {
      expr = pvd->getDefaultArg();
   }
   clang::ASTContext &context = pvd->getASTContext();
   clang::PrintingPolicy policy(context.getPrintingPolicy());
   TTHREAD_TLS_DECL( std::string, buf );
   buf.clear();
   llvm::raw_string_ostream out(buf);
   if (!expr) {
      // CINT returned NULL for non-defaulted args.
      return 0;
   }
   bool implicitInit = false;
   if (const clang::CXXConstructExpr *construct =
       llvm::dyn_cast<clang::CXXConstructExpr>(expr)) {
      implicitInit = (pvd->getInitStyle() == clang::VarDecl::CallInit) &&
         (construct->getNumArgs() == 0) &&
         !construct->isListInitialization();
   }
   if (!implicitInit) {
      if (pvd->getInitStyle() == clang::VarDecl::CallInit) {
         //out << "(";
      }
      else if (pvd->getInitStyle() == clang::VarDecl::CInit) {
         //out << " = ";
      }
      expr->printPretty(out, 0, policy, /*Indentation=*/0);
      if (pvd->getInitStyle() == clang::VarDecl::CallInit) {
         //out << ")";
      }
      out.flush();
   }
   return buf.c_str();
}

const TClingTypeInfo *TClingMethodArgInfo::Type() const
{
   TTHREAD_TLS_DECL_ARG( TClingTypeInfo, ti, fInterp);
   if (!IsValid()) {
      return &ti;
   }
   const clang::ParmVarDecl *pvd = GetDecl();
   clang::QualType qt = pvd->getOriginalType();
   ti.Init(qt);
   return &ti;
}

const char *TClingMethodArgInfo::TypeName() const
{
   if (!IsValid()) {
      return 0;
   }
   return Type()->Name();
}

