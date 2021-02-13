// @(#)root/meta:$Id$
// Author: Philippe Canal November 2013.

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TFunctionTemplate
Dictionary for function template
This class describes one single function template.
*/

#include "TFunctionTemplate.h"
#include "TInterpreter.h"
#include "TClass.h"

ClassImp(TFunctionTemplate);

////////////////////////////////////////////////////////////////////////////////
/// Default TFunctionTemplate ctor.

TFunctionTemplate::TFunctionTemplate(FuncTempInfo_t *info, TClass *cl) : TDictionary(),
   fInfo(info), fClass(cl)
{
   if (fInfo) {
      gCling->FuncTempInfo_Name(fInfo,fName);
      gCling->FuncTempInfo_Title(fInfo,fTitle);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Copy operator.

TFunctionTemplate::TFunctionTemplate(const TFunctionTemplate &orig) : TDictionary(orig)
{
   if (orig.fInfo) {
      fInfo = gCling->FuncTempInfo_FactoryCopy(orig.fInfo);
   } else
      fInfo = 0;
   fClass = orig.fClass;
}

////////////////////////////////////////////////////////////////////////////////
/// Assignment operator.

TFunctionTemplate& TFunctionTemplate::operator=(const TFunctionTemplate &rhs)
{
   if (this != &rhs) {
      gCling->FuncTempInfo_Delete(fInfo);
      if (rhs.fInfo) {
         fInfo = gCling->FuncTempInfo_FactoryCopy(rhs.fInfo);
         gCling->FuncTempInfo_Name(fInfo,fName);
         gCling->FuncTempInfo_Title(fInfo,fTitle);
      } else
         fInfo = 0;
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// TFunctionTemplate dtor deletes adopted CINT FuncTempInfo.

TFunctionTemplate::~TFunctionTemplate()
{
   gCling->FuncTempInfo_Delete(fInfo);
}

////////////////////////////////////////////////////////////////////////////////
/// Clone method.

TObject *TFunctionTemplate::Clone(const char *newname) const
{
   TNamed *newobj = new TFunctionTemplate(*this);
   if (newname && strlen(newname)) newobj->SetName(newname);
   return newobj;
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if this function template object is pointing to a currently
/// loaded function.  If a function is unloaded after the TFunction
/// is created, the TFunction will be set to be invalid.

Bool_t TFunctionTemplate::IsValid()
{
   // Register the transaction when checking the validity of the object.
   if (!fInfo && UpdateInterpreterStateMarker()) {
      // Only for global functions. For data member functions TMethod does it.
      DeclId_t newId = gInterpreter->GetFunction(0, fName);
      if (newId) {
         FuncTempInfo_t *info = gInterpreter->FuncTempInfo_Factory(newId);
         Update(info);
      }
      return newId != 0;
   }
   return fInfo != 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Number of function arguments.

UInt_t TFunctionTemplate::GetTemplateNargs() const
{
   return fInfo ? gCling->FuncTempInfo_TemplateNargs(fInfo) : 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Number of function optional (default) arguments.

UInt_t TFunctionTemplate::GetTemplateMinReqArgs() const
{
   // FIXME: when unload this is an over-estimate.
   return fInfo ? gCling->FuncTempInfo_TemplateMinReqArgs(fInfo) : 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Get property description word. For meaning of bits see EProperty.

Long_t TFunctionTemplate::Property() const
{
   return fInfo ? gCling->FuncTempInfo_Property(fInfo) : 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Get the properties not already defined in Property.See TDictionary's EFunctionProperty.

Long_t TFunctionTemplate::ExtraProperty() const
{
   return fInfo ? gCling->FuncTempInfo_ExtraProperty(fInfo) : 0;
}

////////////////////////////////////////////////////////////////////////////////

TDictionary::DeclId_t TFunctionTemplate::GetDeclId() const
{
   return gInterpreter->GetDeclId(fInfo);
}

////////////////////////////////////////////////////////////////////////////////
/// Update the TFunctionTemplate to reflect the new info.
///
/// This can be used to implement unloading (info == 0) and then reloading
/// (info being the 'new' decl address).

Bool_t TFunctionTemplate::Update(FuncTempInfo_t *info)
{
   if (info == 0) {
      if (fInfo) gCling->FuncTempInfo_Delete(fInfo);
      fInfo = 0;
      return kTRUE;
   } else {
      if (fInfo) gCling->FuncTempInfo_Delete(fInfo);
      fInfo = info;
      gCling->FuncTempInfo_Title(fInfo,fTitle);
      return kTRUE;
   }
}

