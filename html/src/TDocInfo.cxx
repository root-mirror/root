// @(#)root/html:$Id$
// Author: Axel Naumann 2007-01-09

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TDocInfo.h"

#include "TClass.h"

//______________________________________________________________________________
//
// Caches class documentation information, like which module it belongs to,
// and whether THtml should generate documentation for the class.
//______________________________________________________________________________


ClassImp(TClassDocInfo);

const char* TClassDocInfo::GetName() const
{
   // Get the class name, or (UNKNOWN) is no TClass object was found.
   return fClass ? fClass->GetName() : "(UNKNOWN)";
}

ULong_t TClassDocInfo::Hash() const
{
   // Forward to TClass::Hash(), return -1 if no TClass object was found.
   return fClass ? fClass->Hash() : (ULong_t)-1;
}

Int_t TClassDocInfo::Compare(const TObject* obj) const
{
   // Compare two TClassDocInfo objects; used for sorting.
   return fClass ? fClass->Compare(obj) : obj < this;
}

//______________________________________________________________________________
//
// BEGIN_HTML
// <p>Represents modules of the documented product. Modules are sub-groups of
// sources of a product, which get separate index pages and user-provided
// documentation. For ROOT, a module is a sub-directory; it often corresponds
// to a library. TModuleDocInfo, for example, is part of the HTML module,
// which is documented <a href="./HTML_index.html">here</a>. The list of
// all modules is shown e.g. in the <a href="ClassIndex.html">class index</a>.</p>
// <p>A module's documentation is searched by combining its source directory
// (see <a href="#TModuleDocInfo:SetInputDir">SetInputDir()</a>) and the
// module documentation search path defined by
// <a href="./THtml.html#THtml:SetModuleDocPath">THtml::SetModuleDocPath()</a>;
// it defaults to "../doc", i.e. for a module's sources in "module/src" its
// documentation is searched in "module/doc".
// END_HTML
//______________________________________________________________________________


ClassImp(TModuleDocInfo);
