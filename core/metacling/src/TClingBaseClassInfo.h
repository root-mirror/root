// @(#)root/core/meta:$Id$
// Author: Paul Russo   30/07/2012

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TClingBaseClassInfo
#define ROOT_TClingBaseClassInfo

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TClingBaseClassInfo                                                  //
//                                                                      //
// Emulation of the CINT BaseClassInfo class.                           //
//                                                                      //
// The CINT C++ interpreter provides an interface to metadata about     //
// the base classes of a class through the BaseClassInfo class.  This   //
// class provides the same functionality, using an interface as close   //
// as possible to BaseClassInfo but the base class metadata comes from  //
// the Clang C++ compiler, not CINT.                                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TClingClassInfo.h"

#include "clang/AST/DeclCXX.h"

#include <utility>
#include <vector>
#include <string>

namespace cling {
   class Interpreter;
}

class TClingClassInfo;

class TClingBaseClassInfo {

private:

   cling::Interpreter           *fInterp; // Cling interpreter, we do *not* own.
   TClingClassInfo              *fClassInfo; // Class we were intialized with, we own.
   bool                          fFirstTime; // Flag to provide Cint semantics for iterator advancement (not first time)
   bool                          fDescend; // Flag for signaling the need to descend on this advancement.
   const clang::Decl            *fDecl; // Current class whose bases we are iterating through, we do *not* own.
   clang::CXXRecordDecl::base_class_const_iterator fIter; // Current iterator.
   TClingClassInfo              *fBaseInfo; // Base class our iterator is currently pointing at, we own.
   std::vector<std::pair<std::pair<const clang::Decl*, clang::CXXRecordDecl::base_class_const_iterator>, long> > fIterStack; // Iterator stack.
   long                          fOffset; // Offset of the current base, fDecl, in the most-derived class.
   bool                          fClassInfoOwnership; // We created the fClassInfo and we need to delete it in the constructor.

public:

   ~TClingBaseClassInfo() {
      if (fClassInfoOwnership) delete fClassInfo;
      delete fBaseInfo;
   }

   TClingBaseClassInfo(cling::Interpreter*, TClingClassInfo*);
   TClingBaseClassInfo(cling::Interpreter*, TClingClassInfo* derived, TClingClassInfo* base);
   TClingBaseClassInfo(const TClingBaseClassInfo&);
   TClingBaseClassInfo& operator=(const TClingBaseClassInfo&);

   TClingClassInfo *GetBase() const;
   int           InternalNext(int onlyDirect);
   bool          IsValid() const;
   int           Next();
   int           Next(int onlyDirect);
   ptrdiff_t     Offset(void * address = 0, bool isDerivedObject = true) const;
   long          Property() const;
   long          Tagnum() const;
   void          FullName(std::string &output, const ROOT::TMetaUtils::TNormalizedCtxt &normCtxt) const;
   const char   *Name() const;
   const char   *TmpltName() const;

private:
   OffsetPtrFunc_t GenerateBaseOffsetFunction(TClingClassInfo* derivedClass, TClingClassInfo* targetClass, void* address, bool isDerivedObject) const;
};

#endif // ROOT_TClingBaseClassInfo
