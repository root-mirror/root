// @(#)root/core/meta:$Id$
// Author: Paul Russo   30/07/2012

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TClingTypeInfo
Emulation of the CINT TypeInfo class.

The CINT C++ interpreter provides an interface to metadata about
a type through the TypeInfo class.  This class provides the same
functionality, using an interface as close as possible to TypeInfo
but the type metadata comes from the Clang C++ compiler, not CINT.
*/

#include "TClingTypeInfo.h"

#include "TDictionary.h"
#include "Rtypes.h" // for gDebug
#include "TClassEdit.h"
#include "TClingUtils.h"
#include "TInterpreter.h"
#include "ThreadLocalStorage.h"

#include "cling/Interpreter/Interpreter.h"
#include "cling/Interpreter/LookupHelper.h"
#include "cling/Utils/AST.h"

#include "clang/AST/ASTContext.h"
#include "clang/AST/DeclCXX.h"
#include "clang/AST/Type.h"
#include "clang/AST/PrettyPrinter.h"
#include "clang/Frontend/CompilerInstance.h"

#include <cstdio>
#include <string>

using namespace std;

////////////////////////////////////////////////////////////////////////////////

TClingTypeInfo::TClingTypeInfo(cling::Interpreter *interp, const char *name)
   : TClingDeclInfo(nullptr), fInterp(interp)
{
   Init(name);
}

////////////////////////////////////////////////////////////////////////////////

void TClingTypeInfo::Init(const char *name)
{
   fQualType = clang::QualType();
   if (gDebug > 0) {
      fprintf(stderr,
              "TClingTypeInfo::Init(name): looking up clang type: %s", name);
   }
   const cling::LookupHelper& lh = fInterp->getLookupHelper();
   clang::QualType QT = lh.findType(name,
                                    gDebug > 5 ? cling::LookupHelper::WithDiagnostics
                                    : cling::LookupHelper::NoDiagnostics);
   if (QT.isNull()) {
      if (gDebug > 0) {
         fprintf(stderr,
                 "TClingTypeInfo::Init(name): clang type not found: %s", name);
      }
      std::string buf = TClassEdit::InsertStd(name);
      QT = lh.findType(buf,
                       gDebug > 5 ? cling::LookupHelper::WithDiagnostics
                       : cling::LookupHelper::NoDiagnostics);
      if (QT.isNull()) {
         if (gDebug > 0) {
            fprintf(stderr,
                    "TClingTypeInfo::Init(name):  "
                    "clang type not found name: %s\n", buf.c_str());
         }
      }
      else {
         fQualType = QT;
         if (gDebug > 0) {
            fprintf(stderr,
                    "TClingTypeInfo::Init(name): found clang type name: %s\n",
                    buf.c_str());
         }
      }
   }
   else {
      fQualType = QT;
      if (gDebug > 0) {
         fprintf(stderr,
                 "TClingTypeInfo::Init(name): clang type found: %s\n", name);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////

const char *TClingTypeInfo::Name() const
{
   if (!IsValid()) {
      return "";
   }
   // Note: This *must* be static because we are returning a pointer inside it!
   TTHREAD_TLS_DECL( std::string, buf);
   buf.clear();

   R__LOCKGUARD(gInterpreterMutex);
   ROOT::TMetaUtils::GetFullyQualifiedTypeName(buf,fQualType,*fInterp);
   return buf.c_str();
}

////////////////////////////////////////////////////////////////////////////////

long TClingTypeInfo::Property() const
{
   if (!IsValid()) {
      return 0L;
   }
   long property = 0L;
   if (llvm::isa<clang::TypedefType>(*fQualType)) {
      property |= kIsTypedef;
   }
   clang::QualType QT = fQualType.getCanonicalType();
   property = TClingDeclInfo::Property(property, QT);
   const clang::TagType *tagQT = llvm::dyn_cast<clang::TagType>(QT.getTypePtr());
   if (tagQT) {
      // Note: Now we have class, enum, struct, union only.
      const clang::TagDecl *TD = llvm::dyn_cast<clang::TagDecl>(tagQT->getDecl());
      if (TD->isEnum()) {
         property |= kIsEnum;
      } else {
         // Note: Now we have class, struct, union only.
         const clang::CXXRecordDecl *CRD =
            llvm::dyn_cast<clang::CXXRecordDecl>(TD);
         if (CRD->isClass()) {
            property |= kIsClass;
         }
         else if (CRD->isStruct()) {
            property |= kIsStruct;
         }
         else if (CRD->isUnion()) {
            property |= kIsUnion;
         }
         // isAbstract can trigger deserialization
         cling::Interpreter::PushTransactionRAII RAII(fInterp);
         if (CRD->isThisDeclarationADefinition() && CRD->isAbstract()) {
            property |= kIsAbstract;
         }
      }
   }
   return property;
}

////////////////////////////////////////////////////////////////////////////////

int TClingTypeInfo::RefType() const
{
   if (!IsValid()) {
      return 0;
   }
   int cnt = 0;
   bool is_ref = false;
   clang::QualType QT = fQualType.getCanonicalType();
   while (1) {
      if (QT->isArrayType()) {
         QT = llvm::cast<clang::ArrayType>(QT)->getElementType();
         continue;
      }
      else if (QT->isReferenceType()) {
         is_ref = true;
         QT = llvm::cast<clang::ReferenceType>(QT)->getPointeeType();
         continue;
      }
      else if (QT->isPointerType()) {
         ++cnt;
         QT = llvm::cast<clang::PointerType>(QT)->getPointeeType();
         continue;
      }
      else if (QT->isMemberPointerType()) {
         QT = llvm::cast<clang::MemberPointerType>(QT)->getPointeeType();
         continue;
      }
      break;
   }
   int val = 0;
   if (cnt > 1) {
      val = cnt;
   }
   if (is_ref) {
      if (cnt < 2) {
         val = kParaReference;
      }
      else {
         val |= kParaRef;
      }
   }
   return val;
}

////////////////////////////////////////////////////////////////////////////////

int TClingTypeInfo::Size() const
{
   if (!IsValid()) {
      return 1;
   }
   if (fQualType->isDependentType()) {
      // Dependent on a template parameter, we do not know what it is yet.
      return 0;
   }
   if (const clang::RecordType *RT = fQualType->getAs<clang::RecordType>()) {
      if (!RT->getDecl()->getDefinition()) {
         // This is a forward-declared class.
         return 0;
      }
   }
   clang::ASTContext &Context = fInterp->getCI()->getASTContext();
   // Note: This is an int64_t.
   clang::CharUnits::QuantityType Quantity =
      Context.getTypeSizeInChars(fQualType).getQuantity();
   return static_cast<int>(Quantity);
}

////////////////////////////////////////////////////////////////////////////////
/// Return the normalized name of the type (i.e. fully qualified and without
/// the non-opaque typedefs.

const char *TClingTypeInfo::TrueName(const ROOT::TMetaUtils::TNormalizedCtxt &normCtxt) const
{
   if (!IsValid()) {
      return 0;
   }
   // Note: This *must* be static because we are returning a pointer inside it.
   TTHREAD_TLS_DECL( std::string, buf);
   buf.clear();

   ROOT::TMetaUtils::GetNormalizedName(buf,fQualType, *fInterp, normCtxt);

   return buf.c_str();
}

////////////////////////////////////////////////////////////////////////////////
/// Return the normalized name of the type (i.e. fully qualified and without
/// the non-opaque typedefs.

std::string TClingTypeInfo::NormalizedName(const ROOT::TMetaUtils::TNormalizedCtxt &normCtxt) const
{
   if (!IsValid()) {
      return "";
   }
   std::string buf;
   ROOT::TMetaUtils::GetNormalizedName(buf,fQualType, *fInterp, normCtxt);

   // in C++11 this will be efficient thanks to the move constructor.
   return buf;
}

