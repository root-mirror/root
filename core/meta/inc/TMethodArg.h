// @(#)root/meta:$Id$
// Author: Rene Brun   04/02/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TMethodArg
#define ROOT_TMethodArg


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMethodArg                                                           //
//                                                                      //
// Dictionary interface for a method argument.                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TDictionary.h"
#include <string>

#ifdef R__LESS_INCLUDES
class TDataMember;
#else
#include "TDataMember.h"
#endif

class TFunction;
class TMethod;

class TMethodArg : public TDictionary {

friend class TMethod;

private:
   TMethodArg(const TMethodArg&) = delete;
   TMethodArg& operator=(const TMethodArg&) = delete;

   MethodArgInfo_t   *fInfo;         //pointer to CINT method argument info
   TFunction         *fMethod;       //pointer to the method or global function
   TDataMember       *fDataMember;   //TDataMember pointed by this arg,to get values and options from.

public:
   TMethodArg(MethodArgInfo_t *info = nullptr, TFunction *method = nullptr);
   virtual       ~TMethodArg();
   const char    *GetDefault() const;
   TFunction     *GetMethod() const { return fMethod; }
   const char    *GetTypeName() const;
   const char    *GetFullTypeName() const;
   std::string    GetTypeNormalizedName() const;
   Long_t         Property() const;

   TDataMember   *GetDataMember() const;
   TList         *GetOptions() const;

   void           Update(MethodArgInfo_t *info);

   ClassDef(TMethodArg,0)  //Dictionary for a method argument
};

#endif

