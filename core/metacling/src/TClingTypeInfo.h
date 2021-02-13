// @(#)root/core/meta:$Id$
// Author: Paul Russo   30/07/2012

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TClingTypeInfo
#define ROOT_TClingTypeInfo

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TClingTypeInfo                                                       //
//                                                                      //
// Emulation of the CINT TypeInfo class.                                //
//                                                                      //
// The CINT C++ interpreter provides an interface to metadata about     //
// a type through the TypeInfo class.  This class provides the same     //
// functionality, using an interface as close as possible to TypeInfo   //
// but the type metadata comes from the Clang C++ compiler, not CINT.   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TClingDeclInfo.h"
#include <string>

#include "clang/AST/Type.h"

namespace cling {
class Interpreter;
}

namespace ROOT {
   namespace TMetaUtils {
      class TNormalizedCtxt;
   }
}

class TClingTypeInfo final : public TClingDeclInfo {

private:
   cling::Interpreter  *fInterp;    //Cling interpreter, we do *not* own.
   clang::QualType      fQualType;  //Clang qualified type we are querying.

public:

   explicit TClingTypeInfo(cling::Interpreter *interp)
      : TClingDeclInfo(nullptr), fInterp(interp) {}

   TClingTypeInfo(cling::Interpreter *interp, clang::QualType ty)
      : TClingDeclInfo(nullptr), fInterp(interp), fQualType(ty) {}

   TClingTypeInfo(cling::Interpreter *interp, const char *name);

   cling::Interpreter  *GetInterpreter() const { return fInterp; }

   clang::QualType      GetQualType() const { return fQualType; }

   void                 Init(const char *name); // Set type by name.
   void                 Init(clang::QualType ty) { fQualType = ty; }
   bool                 IsValid() const override { return !fQualType.isNull(); }
   const char          *Name() const override; // Get name of type.
   long                 Property() const; // Get properties of type.
   int                  RefType() const; // Get CINT reftype of type.
   int                  Size() const; // Get size in bytes of type.
   const char          *TrueName(const ROOT::TMetaUtils::TNormalizedCtxt &normCtxt) const; // Get name of type with no typedefs.
   std::string          NormalizedName(const ROOT::TMetaUtils::TNormalizedCtxt &normCtxt) const; // Get name of type with no typedefs.

};

#endif // ROOT_TClingTypeInfo
