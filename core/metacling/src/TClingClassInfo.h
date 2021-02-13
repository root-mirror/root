// @(#)root/core/meta:$Id$
// Author: Paul Russo   30/07/2012

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TClingClassInfo
#define ROOT_TClingClassInfo

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TClingClassInfo                                                      //
//                                                                      //
// Emulation of the CINT ClassInfo class.                               //
//                                                                      //
// The CINT C++ interpreter provides an interface to metadata about     //
// a class through the ClassInfo class.  This class provides the same   //
// functionality, using an interface as close as possible to ClassInfo  //
// but the class metadata comes from the Clang C++ compiler, not CINT.  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TClingDeclInfo.h"
#include "TClingMethodInfo.h"
#include "TClingUtils.h"
#include "TDataType.h"
#include "TDictionary.h"

#include <vector>
#include <string>
#include <utility>

#include "llvm/ADT/DenseMap.h"

namespace cling {
   class Interpreter;
}

namespace clang {
   class CXXMethodDecl;
   class FunctionTemplateDecl;
   class ValueDecl;
}

namespace ROOT {
   namespace TMetaUtils {
      class TNormalizedCtxt;
   }
}

extern "C" typedef ptrdiff_t (*OffsetPtrFunc_t)(void*, bool);

class TClingClassInfo final : public TClingDeclInfo {

private:

   cling::Interpreter   *fInterp; // Cling interpreter, we do *not* own.
   bool                  fFirstTime: 1; // We need to skip the first increment to support the cint Next() semantics.
   bool                  fDescend : 1;  // Flag for signaling the need to descend on this advancement.
   bool                  fIterAll : 1;  // Flag whether iteration should be as complete as possible.
   bool                  fIsIter : 1;   // Flag whether this object was setup for iteration.
   clang::DeclContext::decl_iterator fIter; // Current decl in scope.
   const clang::Type    *fType; // Type representing the decl (conserves typedefs like Double32_t). (we do *not* own)
   std::vector<clang::DeclContext::decl_iterator> fIterStack; // Recursion stack for traversing nested scopes.
   std::string           fTitle; // The meta info for the class.
   std::string           fDeclFileName; // Name of the file where the underlying entity is declared.
   llvm::DenseMap<const clang::Decl*, std::pair<ptrdiff_t, OffsetPtrFunc_t> > fOffsetCache; // Functions already generated for offsets.

   explicit TClingClassInfo() = delete;
   TClingClassInfo &operator=(const TClingClassInfo &) = delete;
public: // Types

   enum EInheritanceMode {
      kInThisScope = 0,
      kWithInheritance = 1
   };

public:

   explicit TClingClassInfo(cling::Interpreter *, Bool_t all = kTRUE);
   explicit TClingClassInfo(cling::Interpreter *, const char *classname, bool intantiateTemplate = kTRUE);
   explicit TClingClassInfo(cling::Interpreter *, const clang::Type &);
   explicit TClingClassInfo(cling::Interpreter *, const clang::Decl *);
   void                 AddBaseOffsetFunction(const clang::Decl* decl, OffsetPtrFunc_t func) { fOffsetCache[decl] = std::make_pair(0L, func); }
   void                 AddBaseOffsetValue(const clang::Decl* decl, ptrdiff_t offset);
   long                 ClassProperty() const;
   void                 Delete(void *arena, const ROOT::TMetaUtils::TNormalizedCtxt &normCtxt) const;
   void                 DeleteArray(void *arena, bool dtorOnly, const ROOT::TMetaUtils::TNormalizedCtxt &normCtxt) const;
   void                 Destruct(void *arena, const ROOT::TMetaUtils::TNormalizedCtxt &normCtxt) const;
   const clang::ValueDecl *GetDataMember(const char *name) const;
   void SetDecl(const clang::Decl* D) {
     // FIXME: We should track down all sets and potentially avoid them.
     fDecl = D;
     fNameCache.clear(); // invalidate the cache.
   }
   TDictionary::DeclId_t   GetDeclId() const {
      if (!fDecl)
        return nullptr;
      return (const clang::Decl*)(fDecl->getCanonicalDecl());
   }
   const clang::FunctionTemplateDecl *GetFunctionTemplate(const char *fname) const;
   TClingMethodInfo     GetMethod(const char *fname) const;
   TClingMethodInfo     GetMethod(const char *fname, const char *proto,
                                  long *poffset, ROOT::EFunctionMatchMode mode = ROOT::kConversionMatch,
                                  EInheritanceMode imode = kWithInheritance) const;
   TClingMethodInfo     GetMethodWithArgs(const char *fname, const char *arglist,
                                  long *poffset, ROOT::EFunctionMatchMode mode = ROOT::kConversionMatch,
                                  EInheritanceMode imode = kWithInheritance) const;
   TClingMethodInfo     GetMethod(const char *fname, const char *proto, bool objectIsConst,
                                  long *poffset, ROOT::EFunctionMatchMode mode = ROOT::kConversionMatch,
                                  EInheritanceMode imode = kWithInheritance) const;
   TClingMethodInfo     GetMethodWithArgs(const char *fname, const char *arglist, bool objectIsConst,
                                  long *poffset, ROOT::EFunctionMatchMode mode = ROOT::kConversionMatch,
                                  EInheritanceMode imode = kWithInheritance) const;
   TClingMethodInfo     GetMethod(const char *fname, const llvm::SmallVectorImpl<clang::QualType> &proto,
                                  long *poffset, ROOT::EFunctionMatchMode mode = ROOT::kConversionMatch,
                                  EInheritanceMode imode = kWithInheritance) const;
   TClingMethodInfo     GetMethod(const char *fname, const llvm::SmallVectorImpl<clang::QualType> &proto, bool objectIsConst,
                                  long *poffset, ROOT::EFunctionMatchMode mode = ROOT::kConversionMatch,
                                  EInheritanceMode imode = kWithInheritance) const;
   int                  GetMethodNArg(const char *method, const char *proto, Bool_t objectIsConst, ROOT::EFunctionMatchMode mode = ROOT::kConversionMatch) const;
   long                 GetOffset(const clang::CXXMethodDecl* md) const;
   ptrdiff_t            GetBaseOffset(TClingClassInfo* toBase, void* address, bool isDerivedObject);
   const clang::Type   *GetType() const { return fType; } // Underlying representation with Double32_t
   std::vector<std::string> GetUsingNamespaces();
   ROOT::TMetaUtils::EIOCtorCategory HasDefaultConstructor(bool checkio = false, std::string *type_name = nullptr) const;
   bool                 HasMethod(const char *name) const;
   void                 Init(const char *name);
   void                 Init(const clang::Decl*);
   void                 Init(int tagnum);
   void                 Init(const clang::Type &);
   bool                 IsBase(const char *name) const;
   static bool          IsEnum(cling::Interpreter *interp, const char *name);
   bool                 IsScopedEnum() const;
   EDataType            GetUnderlyingType() const;
   bool                 IsLoaded() const;
   bool                 IsValidMethod(const char *method, const char *proto, Bool_t objectIsConst, long *offset, ROOT::EFunctionMatchMode mode = ROOT::kConversionMatch) const;
   int                  InternalNext();
   int                  Next();
   void                *New(const ROOT::TMetaUtils::TNormalizedCtxt &normCtxt) const;
   void                *New(int n, const ROOT::TMetaUtils::TNormalizedCtxt &normCtxt) const;
   void                *New(int n, void *arena, const ROOT::TMetaUtils::TNormalizedCtxt &normCtxt) const;
   void                *New(void *arena, const ROOT::TMetaUtils::TNormalizedCtxt &normCtxt) const;
   long                 Property() const;
   int                  RootFlag() const;
   int                  Size() const;
   long                 Tagnum() const;
   const char          *FileName();
   void                 FullName(std::string &output, const ROOT::TMetaUtils::TNormalizedCtxt &normCtxt) const;
   const char          *Title();
   const char          *TmpltName() const;

};

#endif // ROOT_TClingClassInfo
