// @(#)root/meta:$Id$
// Author: Rene Brun   04/02/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TMethodArg.h"
#include "TInterpreter.h"
#include "Strlen.h"
#include "TMethod.h"
#include "TMethodCall.h"
#include "TDataMember.h"

/** \class TMethodArg
Each ROOT method (see TMethod) has a linked list of its arguments.
This class describes one single method argument.
The method argument info is obtained via the CINT api.
See class TCling.

The method argument information is used a.o. in the TContextMenu
and THtml classes.
*/

ClassImp(TMethodArg);

////////////////////////////////////////////////////////////////////////////////
/// Default TMethodArg ctor. TMethodArgs are constructed in TFunction
/// via a call to TCling::CreateListOfMethodArgs().

TMethodArg::TMethodArg(MethodArgInfo_t *info, TFunction *method) : TDictionary()
{
   fDataMember = 0;
   fInfo       = info;
   fMethod     = method;
   if (fInfo) {
      SetName(gCling->MethodArgInfo_Name(fInfo));
      SetTitle(gCling->MethodArgInfo_TypeName(fInfo));
   }
}

////////////////////////////////////////////////////////////////////////////////
/// TMethodArg dtor deletes adopted CINT MethodArgInfo object.

TMethodArg::~TMethodArg()
{
   if (fInfo) gCling->MethodArgInfo_Delete(fInfo);
}

////////////////////////////////////////////////////////////////////////////////
/// Get default value of method argument.

const char *TMethodArg::GetDefault() const
{
   return gCling->MethodArgInfo_DefaultValue(fInfo);
}

////////////////////////////////////////////////////////////////////////////////
/// Get type of method argument, e.g.: "class TDirectory*" -> "TDirectory"
/// Result needs to be used or copied immediately.

const char *TMethodArg::GetTypeName() const
{
   return gCling->TypeName(gCling->MethodArgInfo_TypeName(fInfo));
}

////////////////////////////////////////////////////////////////////////////////
/// Get full type description of method argument, e.g.: "class TDirectory*".

const char *TMethodArg::GetFullTypeName() const
{
   return gCling->MethodArgInfo_TypeName(fInfo);
}

////////////////////////////////////////////////////////////////////////////////
/// Get the normalized name of the return type. A normalized name is fully
/// qualified and has all typedef desugared except for the 'special' typedef
/// which include Double32_t, Float16_t, [U]Long64_t and std::string.  It
/// also has std:: removed [This is subject to change].

std::string TMethodArg::GetTypeNormalizedName() const
{
   return gCling->MethodArgInfo_TypeNormalizedName(fInfo);
}

////////////////////////////////////////////////////////////////////////////////
/// Get property description word. For meaning of bits see EProperty.

Long_t TMethodArg::Property() const
{
   return gCling->MethodArgInfo_Property(fInfo);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns list of possible options - according to pointed datamember.
/// If there is no datamember field assigned to this methodarg - returns 0.

TList *TMethodArg::GetOptions() const
{
   return (TList*)(fDataMember ? fDataMember->GetOptions() : 0);
}

////////////////////////////////////////////////////////////////////////////////
/// Returns TDataMember pointed by this methodarg.
/// If you want to specify list of options or current value for your
/// MethodArg (i.e. it is used as initial values in argument-asking dialogs
/// popped up from context-meny),you can get this value from one of data
/// members of the class.
///
/// The only restriction is, that this DataMember object must have its
/// Getter/Setter methods set-up correctly - for details look at TDataMember.
/// To learn how to specify the data member to which the argument should
/// "point", look at TMethod. This is TMethod which sets up fDataMember,
/// so it could work correctly.

TDataMember *TMethodArg::GetDataMember() const
{
   return fDataMember;
}

////////////////////////////////////////////////////////////////////////////////
/// Update fInfo (to 0 for unloading and non-zero for reloading).
/// This takes ownership of 'info'

void TMethodArg::Update(MethodArgInfo_t *info)
{
   if (fInfo) gCling->MethodArgInfo_Delete(fInfo);
   fInfo = info;
   if (fInfo) {
      SetTitle(gCling->MethodArgInfo_TypeName(fInfo));
   }
}
