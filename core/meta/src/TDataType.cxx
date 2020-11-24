// @(#)root/meta:$Id$
// Author: Rene Brun   04/02/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons .               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TDataType
Basic data type descriptor (datatype information is obtained from
CINT). This class describes the attributes of type definitions
(typedef's). The TROOT class contains a list of all currently
defined types (accessible via TROOT::GetListOfTypes()).
*/

#include "TDataType.h"
#include "TInterpreter.h"
#include "TCollection.h"
#include "TVirtualMutex.h"
#include "ThreadLocalStorage.h"
#ifdef R__SOLARIS
#include <typeinfo>
#endif

ClassImp(TDataType);

TDataType* TDataType::fgBuiltins[kNumDataTypes] = {0};

////////////////////////////////////////////////////////////////////////////////
/// Default TDataType ctor. TDataTypes are constructed in TROOT via
/// a call to TCling::UpdateListOfTypes().

TDataType::TDataType(TypedefInfo_t *info) : TDictionary(),
   fTypeNameIdx(-1), fTypeNameLen(0)
{
   fInfo = info;

   if (fInfo) {
      R__LOCKGUARD(gInterpreterMutex);
      SetName(gCling->TypedefInfo_Name(fInfo));
      SetTitle(gCling->TypedefInfo_Title(fInfo));
      SetType(gCling->TypedefInfo_TrueName(fInfo));
      fProperty = gCling->TypedefInfo_Property(fInfo);
      fSize = gCling->TypedefInfo_Size(fInfo);
   } else {
      SetTitle("Builtin basic type");
      fProperty = 0;
      fSize = 0;
      fType = kNoType_t;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor for basic data types, like "char", "unsigned char", etc.

TDataType::TDataType(const char *typenam) : fInfo(0), fProperty(kIsFundamental),
   fTypeNameIdx(-1), fTypeNameLen(0)
{
   fInfo = 0;
   SetName(typenam);
   SetTitle("Builtin basic type");

   SetType(fName.Data());
}

////////////////////////////////////////////////////////////////////////////////
///copy constructor

TDataType::TDataType(const TDataType& dt) :
  TDictionary(dt),
  fInfo(gCling->TypedefInfo_FactoryCopy(dt.fInfo)),
  fSize(dt.fSize),
  fType(dt.fType),
  fProperty(dt.fProperty),
  fTrueName(dt.fTrueName),
  fTypeNameIdx(dt.fTypeNameIdx), fTypeNameLen(dt.fTypeNameLen)
{
}

////////////////////////////////////////////////////////////////////////////////
/// assignment operator

TDataType& TDataType::operator=(const TDataType& dt)
{
   if(this!=&dt) {
      TDictionary::operator=(dt);
      gCling->TypedefInfo_Delete(fInfo);
      fInfo=gCling->TypedefInfo_FactoryCopy(dt.fInfo);
      fSize=dt.fSize;
      fType=dt.fType;
      fProperty=dt.fProperty;
      fTrueName=dt.fTrueName;
      fTypeNameIdx=dt.fTypeNameIdx;
      fTypeNameLen=dt.fTypeNameLen;
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// TDataType dtor deletes adopted CINT TypedefInfo object.

TDataType::~TDataType()
{
   gCling->TypedefInfo_Delete(fInfo);
}

////////////////////////////////////////////////////////////////////////////////
/// Return the name of the type.

const char *TDataType::GetTypeName(EDataType type)
{
   switch (type) {
      case  1: return "Char_t";
      case  2: return "Short_t";
      case  3: return "Int_t";
      case  4: return "Long_t";
      case  5: return "Float_t";
      case  6: return "Int_t";
      case  7: return "char*";
      case  8: return "Double_t";
      case  9: return "Double32_t";
      case 11: return "UChar_t";
      case 12: return "UShort_t";
      case 13: return "UInt_t";
      case 14: return "ULong_t";
      case 15: return "UInt_t";
      case 16: return "Long64_t";
      case 17: return "ULong64_t";
      case 18: return "Bool_t";
      case 19: return "Float16_t";
      case kVoid_t: return "void";
      case kDataTypeAliasUnsigned_t: return "UInt_t";
      case kDataTypeAliasSignedChar_t: return "SignedChar_t";
      case kOther_t:  return "";
      case kNoType_t: return "";
      case kchar:     return "Char_t";
      default: return "";
   }
   return ""; // to silence compilers
}

////////////////////////////////////////////////////////////////////////////////
/// Get basic type of typedef, e,g.: "class TDirectory*" -> "TDirectory".
/// Result needs to be used or copied immediately.

TString TDataType::GetTypeName()
{
   if (fTypeNameLen) {
     return fTrueName(fTypeNameIdx, fTypeNameLen);
   }

   if (fInfo) {
      (const_cast<TDataType*>(this))->CheckInfo();
      TString typeName = gInterpreter->TypeName(fTrueName.Data());
      fTypeNameIdx = fTrueName.Index(typeName);
      if (fTypeNameIdx == -1) {
         Error("GetTypeName", "Cannot find type name %s in true name %s!",
               typeName.Data(), fTrueName.Data());
         return fName;
      }
      fTypeNameLen = typeName.Length();
      return fTrueName(fTypeNameIdx, fTypeNameLen);
   } else {
      if (fType != kOther_t) return fName.Data();
      return fTrueName;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get full type description of typedef, e,g.: "class TDirectory*".

const char *TDataType::GetFullTypeName() const
{
   if (fInfo) {
      (const_cast<TDataType*>(this))->CheckInfo();
      return fTrueName;
   } else {
     if (fType != kOther_t) return fName;
     return fTrueName;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set type id depending on name.

EDataType TDataType::GetType(const std::type_info &typeinfo)
{
   EDataType retType = kOther_t;

   if (!strcmp(typeid(unsigned int).name(), typeinfo.name())) {
      retType = kUInt_t;
   } else if (!strcmp(typeid(int).name(), typeinfo.name())) {
      retType = kInt_t;
   } else if (!strcmp(typeid(ULong_t).name(), typeinfo.name())) {
      retType = kULong_t;
   } else if (!strcmp(typeid(Long_t).name(), typeinfo.name())) {
      retType = kLong_t;
   } else if (!strcmp(typeid(ULong64_t).name(), typeinfo.name())) {
      retType = kULong64_t;
   } else if (!strcmp(typeid(Long64_t).name(), typeinfo.name())) {
      retType = kLong64_t;
   } else if (!strcmp(typeid(unsigned short).name(), typeinfo.name())) {
      retType = kUShort_t;
   } else if (!strcmp(typeid(short).name(), typeinfo.name())) {
      retType = kShort_t;
   } else if (!strcmp(typeid(unsigned char).name(), typeinfo.name())) {
      retType = kUChar_t;
   } else if (!strcmp(typeid(char).name(), typeinfo.name())) {
      retType = kChar_t;
   } else if (!strcmp(typeid(Bool_t).name(), typeinfo.name())) {
      retType = kBool_t;
   } else if (!strcmp(typeid(float).name(), typeinfo.name())) {
      retType = kFloat_t;
   } else if (!strcmp(typeid(Float16_t).name(), typeinfo.name())) {
      retType = kFloat16_t;
   } else if (!strcmp(typeid(double).name(), typeinfo.name())) {
      retType = kDouble_t;
   } else if (!strcmp(typeid(Double32_t).name(), typeinfo.name())) {
      retType = kDouble32_t;
   } else if (!strcmp(typeid(char*).name(), typeinfo.name())) {
      retType = kCharStar;
   } else if (!strcmp(typeid(signed char).name(), typeinfo.name())) {
      retType = kDataTypeAliasSignedChar_t;
   }
   return retType;
}

////////////////////////////////////////////////////////////////////////////////
/// Return string containing value in buffer formatted according to
/// the basic data type. The result needs to be used or copied immediately.

const char *TDataType::AsString(void *buf) const
{
   TTHREAD_TLS_DECL_ARG(TString, line ,81);
   const char *name;

   if (fInfo) {
      (const_cast<TDataType*>(this))->CheckInfo();
      name = fTrueName;
   } else {
      name = fName.Data();
   }

   line[0] = 0;
   if (!strcmp("unsigned int", name))
      line.Form( "%u", *(unsigned int *)buf);
   else if (!strcmp("unsigned", name))
      line.Form( "%u", *(unsigned int *)buf);
   else if (!strcmp("int", name))
      line.Form( "%d", *(int *)buf);
   else if (!strcmp("unsigned long", name))
      line.Form( "%lu", *(ULong_t *)buf);
   else if (!strcmp("long", name))
      line.Form( "%ld", *(Long_t *)buf);
   else if (!strcmp("unsigned long long", name))
      line.Form( "%llu", *(ULong64_t *)buf);
   else if (!strcmp("ULong64_t", name))
      line.Form( "%llu", *(ULong64_t *)buf);
   else if (!strcmp("long long", name))
      line.Form( "%lld", *(Long64_t *)buf);
   else if (!strcmp("Long64_t", name))
      line.Form( "%lld", *(Long64_t *)buf);
   else if (!strcmp("unsigned short", name))
      line.Form( "%hu", *(unsigned short *)buf);
   else if (!strcmp("short", name))
      line.Form( "%hd", *(short *)buf);
   else if (!strcmp("bool", name))
      line.Form( "%s", *(Bool_t *)buf ? "true" : "false");
   else if (!strcmp("unsigned char", name) || !strcmp("char", name) ) {
      line = *(char*)buf;
   } else if (!strcmp("float", name))
      line.Form( "%g", *(float *)buf);
   else if (!strcmp("double", name))
      line.Form( "%g", *(double *)buf);
   else if (!strcmp("char*", name))
      line.Form( "%s", *(char**)buf);

   return line;
}

////////////////////////////////////////////////////////////////////////////////
/// Get property description word. For meaning of bits see EProperty.

Long_t TDataType::Property() const
{
   if (fInfo) (const_cast<TDataType*>(this))->CheckInfo();
   return fProperty;
}

////////////////////////////////////////////////////////////////////////////////
/// Set type id depending on name.

void TDataType::SetType(const char *name)
{
   fTrueName = name;
   fType = kOther_t;
   fSize = 0;

   if (name==0) {
      return;
   } else if (!strcmp("unsigned int", name)) {
      fType = kUInt_t;
      fSize = sizeof(UInt_t);
   } else if (!strcmp("unsigned", name)) {
      fType = kUInt_t;
      fSize = sizeof(UInt_t);
   } else if (!strcmp("int", name)) {
      fType = kInt_t;
      fSize = sizeof(Int_t);
   } else if (!strcmp("unsigned long", name)) {
      fType = kULong_t;
      fSize = sizeof(ULong_t);
   } else if (!strcmp("long", name)) {
      fType = kLong_t;
      fSize = sizeof(Long_t);
   } else if (!strcmp("unsigned long long", name) || !strcmp("ULong64_t",name)) {
      fType = kULong64_t;
      fSize = sizeof(ULong64_t);
   } else if (!strcmp("long long", name) || !strcmp("Long64_t",name)) {
      fType = kLong64_t;
      fSize = sizeof(Long64_t);
   } else if (!strcmp("unsigned short", name)) {
      fType = kUShort_t;
      fSize = sizeof(UShort_t);
   } else if (!strcmp("short", name)) {
      fType = kShort_t;
      fSize = sizeof(Short_t);
   } else if (!strcmp("unsigned char", name)) {
      fType = kUChar_t;
      fSize = sizeof(UChar_t);
   } else if (!strcmp("char", name)) {
      fType = kChar_t;
      fSize = sizeof(Char_t);
   } else if (!strcmp("bool", name)) {
      fType = kBool_t;
      fSize = sizeof(Bool_t);
   } else if (!strcmp("float", name)) {
      fType = kFloat_t;
      fSize = sizeof(Float_t);
   } else if (!strcmp("double", name)) {
      fType = kDouble_t;
      fSize = sizeof(Double_t);
   } else if (!strcmp("signed char", name)) {
      fType = kChar_t; // kDataTypeAliasSignedChar_t;
      fSize = sizeof(Char_t);
   } else if (!strcmp("void", name)) {
      fType = kVoid_t;
      fSize = 0;
   }

   if (!strcmp("Float16_t", fName.Data())) {
      fType = kFloat16_t;
   }
   if (!strcmp("Double32_t", fName.Data())) {
      fType = kDouble32_t;
   }
   if (!strcmp("char*",fName.Data())) {
      fType = kCharStar;
   }
   // kCounter =  6, kBits     = 15
}

////////////////////////////////////////////////////////////////////////////////
/// Get size of basic typedef'ed type.

Int_t TDataType::Size() const
{
   if (fInfo) (const_cast<TDataType*>(this))->CheckInfo();
   return fSize;
}

////////////////////////////////////////////////////////////////////////////////
/// Refresh the underlying information.

void TDataType::CheckInfo()
{
   // This can be needed if the library defining this typedef was loaded after
   // another library and that this other library is unloaded (in which case
   // things can get renumbered inside CINT).

   if (!fInfo) return;

   // This intentionally cast the constness away so that
   // we can call CheckInfo from const data members.
   R__LOCKGUARD(gInterpreterMutex);

   if (!gCling->TypedefInfo_IsValid(fInfo) ||
       strcmp(gCling->TypedefInfo_Name(fInfo),fName.Data())!=0) {

      // The fInfo is invalid or does not
      // point to this typedef anymore, let's
      // refresh it

      gCling->TypedefInfo_Init(fInfo, fName.Data());

      if (!gCling->TypedefInfo_IsValid(fInfo)) return;

      SetTitle(gCling->TypedefInfo_Title(fInfo));
      SetType(gCling->TypedefInfo_TrueName(fInfo));
      fProperty = gCling->TypedefInfo_Property(fInfo);
      fSize = gCling->TypedefInfo_Size(fInfo);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Create the TDataType objects for builtins.

void TDataType::AddBuiltins(TCollection* types)
{
   if (fgBuiltins[kChar_t] == 0) {
      // Add also basic types (like a identity typedef "typedef int int")
      fgBuiltins[kChar_t] = new TDataType("char");
      fgBuiltins[kUChar_t] = new TDataType("unsigned char");
      fgBuiltins[kShort_t] = new TDataType("short");
      fgBuiltins[kUShort_t] = new TDataType("unsigned short");
      fgBuiltins[kInt_t] = new TDataType("int");
      fgBuiltins[kUInt_t] = new TDataType("unsigned int");
      fgBuiltins[kLong_t] = new TDataType("long");
      fgBuiltins[kULong_t] = new TDataType("unsigned long");
      fgBuiltins[kLong64_t] = new TDataType("long long");
      fgBuiltins[kULong64_t] = new TDataType("unsigned long long");
      fgBuiltins[kFloat_t] = new TDataType("float");
      fgBuiltins[kDouble_t] = new TDataType("double");
      fgBuiltins[kVoid_t] = new TDataType("void");
      fgBuiltins[kBool_t] = new TDataType("bool");
      fgBuiltins[kCharStar] = new TDataType("char*");

      fgBuiltins[kDataTypeAliasUnsigned_t] = new TDataType("unsigned");
      fgBuiltins[kDataTypeAliasSignedChar_t] = new TDataType("signed char");
   }

   for (Int_t i = 0; i < (Int_t)kNumDataTypes; ++i) {
      if (fgBuiltins[i]) types->Add(fgBuiltins[i]);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Given a EDataType type, get the TDataType* that represents it.

TDataType* TDataType::GetDataType(EDataType type)
{
   if (type == kOther_t || type >= kNumDataTypes) return 0;
   return fgBuiltins[(int)type];
}
