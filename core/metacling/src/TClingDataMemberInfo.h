// @(#)root/core/meta:$Id$
// Author: Paul Russo   30/07/2012

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TClingDataMemberInfo
#define ROOT_TClingDataMemberInfo

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TClingDataMemberInfo                                                 //
//                                                                      //
// Emulation of the CINT DataMemberInfo class.                          //
//                                                                      //
// The CINT C++ interpreter provides an interface to metadata about     //
// the data members of a class through the DataMemberInfo class.  This  //
// class provides the same functionality, using an interface as close   //
// as possible to DataMemberInfo but the data member metadata comes     //
// from the Clang C++ compiler, not CINT.                               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TClingDeclInfo.h"
#include "TClingMemberIter.h"
#include "TDictionary.h"

#include <vector>
#include <string>

namespace clang {
   class Decl;
   class Type;
   class ValueDecl;
}

namespace ROOT {
   namespace TMetaUtils {
      class TNormalizedCtxt;
   }
}

class TClingClassInfo;

/// Iterate over VarDecl, FieldDecl, EnumConstantDecl, IndirectFieldDecl, and 
/// UsingShadowDecls thereof, within a scope, recursing through "transparent"
/// scopes (see DCIter::HandleInlineDeclContext()).
class TClingDataMemberIter final: public TClingMemberIter {
   TDictionary::EMemberSelection fSelection = TDictionary::EMemberSelection::kNoUsingDecls;

protected:
   // TODO:
   //const clang::Decl *
   //InstantiateTemplateWithDefaults(const clang::RedeclarableTemplateDecl *TD) const final;

   bool ShouldSkip(const clang::Decl* FD) const final;
   bool ShouldSkip(const clang::UsingShadowDecl* USD) const final;

public:
   TClingDataMemberIter() = default;
   TClingDataMemberIter(cling::Interpreter *interp, clang::DeclContext *DC, TDictionary::EMemberSelection selection)
      : TClingMemberIter(interp, DC), fSelection(selection)
   {
   }
};

class TClingDataMemberInfo final : public TClingDeclInfo {

private:

   cling::Interpreter    *fInterp;    // Cling interpreter, we do *not* own.
   TClingClassInfo       *fClassInfo = nullptr; // ClassInfo for the decl context, for X<Float16_t> vs X<float>.
   TClingDataMemberIter   fIter; // Current decl.
   std::string            fTitle; // The meta info for the member.
   bool                   fFirstTime = true; // We need to skip the first increment to support the cint Next() semantics.

   mutable std::string fIoType;
   mutable std::string fIoName;
   union {
      float fFloat;
      double fDouble;
      long fLong;
   } fConstInitVal; // Result of VarDecl::evaluateValue()
   inline void CheckForIoTypeAndName () const;

   // Invalidate the name caches.
   void ClearNames() {
      fNameCache.clear();
      fIoType.clear();
      fIoName.clear();
   }

public:

   explicit TClingDataMemberInfo(cling::Interpreter *interp)
      : TClingDeclInfo(nullptr), fInterp(interp) {}

   TClingDataMemberInfo(cling::Interpreter *interp, TClingClassInfo *ci, TDictionary::EMemberSelection selection);

   // Takes concrete decl and disables the iterator.
   // ValueDecl is the common base between enum constant, var decl and field decl
   TClingDataMemberInfo(cling::Interpreter *interp, const clang::ValueDecl *, TClingClassInfo *);

   typedef TDictionary::DeclId_t DeclId_t;

   int                ArrayDim() const;
   const clang::Decl *GetDecl() const override {
     if (const clang::Decl* SingleDecl = TClingDeclInfo::GetDecl())
       return SingleDecl;
     return *fIter;
   }
   const clang::ValueDecl       *GetAsValueDecl() const;
   const clang::UsingShadowDecl *GetAsUsingShadowDecl() const;

   /// Get the ValueDecl, or if this represents a UsingShadowDecl, the underlying target ValueDecl.
   const clang::ValueDecl       *GetTargetValueDecl() const;
   DeclId_t           GetDeclId() const;
   const clang::Type *GetClassAsType() const;
   int                MaxIndex(int dim) const;
   int                Next();
   long               Offset();
   long               Property() const;
   long               TypeProperty() const;
   int                TypeSize() const;
   const char        *TypeName() const;
   const char        *TypeTrueName(const ROOT::TMetaUtils::TNormalizedCtxt &normCtxt) const;
   const char        *Name() const override;
   const char        *Title();
   llvm::StringRef    ValidArrayIndex() const;

};

#endif // ROOT_TClingDataMemberInfo
