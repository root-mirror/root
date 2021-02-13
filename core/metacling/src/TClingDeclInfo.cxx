/// \file TClingDeclInfo.cxx
///
/// \brief The file contains a base class of TCling*Info classes.
///
/// \author Vassil Vassilev <vvasilev@cern.ch>
///
/// \date March, 2019
///
/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TClingDeclInfo.h"

#include "TDictionary.h"

#include "clang/AST/ASTContext.h"
#include "clang/AST/Decl.h"

using namespace clang;

// pin the vtable here.
TClingDeclInfo::~TClingDeclInfo() {}

const char* TClingDeclInfo::Name() const
{
   if (!IsValid())
      return 0;

   if (!fNameCache.empty())
      return fNameCache.c_str();

   const Decl* D = GetDecl();
   if (!isa<NamedDecl>(D))
      return 0;

   const NamedDecl* ND = cast<NamedDecl>(D);
   clang::PrintingPolicy policy(ND->getASTContext().getPrintingPolicy());
   llvm::raw_string_ostream stream(fNameCache);
   ND->getNameForDiagnostic(stream, policy, /*Qualified=*/false);
   stream.flush();
   return fNameCache.c_str();
}

long TClingDeclInfo::Property(long property, clang::QualType &qt) const
{
   if (!IsValid()) {
      return 0L;
   }
   if (qt.isConstQualified()) {
      property |= kIsConstant;
   }
   while (1) {
      if (qt->isArrayType()) {
         property |= kIsArray;
         qt = llvm::cast<clang::ArrayType>(qt)->getElementType();
         continue;
      }
      else if (qt->isReferenceType()) {
         property |= kIsReference;
         qt = llvm::cast<clang::ReferenceType>(qt)->getPointeeType();
         continue;
      }
      else if (qt->isPointerType()) {
         property |= kIsPointer;
         if (qt.isConstQualified()) {
            property |= kIsConstPointer;
         }
         qt = llvm::cast<clang::PointerType>(qt)->getPointeeType();
         continue;
      }
      else if (qt->isMemberPointerType()) {
         qt = llvm::cast<clang::MemberPointerType>(qt)->getPointeeType();
         continue;
      }
      break;
   }
   if (qt->isBuiltinType()) {
      property |= kIsFundamental;
   }
   if (qt.isConstQualified()) {
      property |= kIsConstant;
   }
   return property;
}
