// @(#)root/xml:$Id$
// Author: Sergey Linev  10.05.2004

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//________________________________________________________________________
//
// Class TXMLSetup is used as storage of xml file settings
// This class is used in TXMLFile and in TXmlBuffer classes.
// Xml settings can be coded via a string in following format
//
//   "2xoo"
//    ||| \ .
//    || \ usage of name spaces.
//    | \ usage of DTD;
//     \ storage of TStreamerInfo objects in file;
//      layout of xml file (= 2 - specialized (default), = 3 - generic)
//
// For last three boolean parameters "x" means true, "o" - false
//
// Such string can be set as argument of TXMLFile constructor. In that
// case new TXMLFile with such parameters will be created.
// These settings automatically stored in xml file.

//________________________________________________________________________

#include "TXMLSetup.h"

#include "TROOT.h"
#include "TList.h"
#include "TClass.h"
#include "TStreamerElement.h"

#include <iostream>
#include <cstdlib>

ClassImp(TXMLSetup);

namespace xmlio {

const char *Root = "root";
const char *Setup = "setup";
const char *ClassVersion = "version";
const char *IOVersion = "version";
const char *OnlyVersion = "Version";
const char *Ptr = "ptr";
const char *Ref = "ref";
const char *Null = "null";
const char *IdBase = "id";
const char *Size = "size";
const char *Xmlobject = "XmlObject";
const char *Xmlkey = "XmlKey";
const char *Cycle = "cycle";
const char *XmlBlock = "XmlBlock";
const char *Zip = "zip";
const char *Object = "Object";
const char *ObjClass = "class";
const char *Class = "Class";
const char *Member = "Member";
const char *Item = "Item";
const char *Name = "name";
const char *Title = "title";
const char *CreateTm = "created";
const char *ModifyTm = "modified";
const char *ObjectUUID = "uuid";
const char *Type = "type";
const char *Value = "value";
const char *v = "v";
const char *cnt = "cnt";
const char *True = "true";
const char *False = "false";
const char *SInfos = "StreamerInfos";

const char *Array = "Array";
const char *Bool = "Bool_t";
const char *Char = "Char_t";
const char *Short = "Short_t";
const char *Int = "Int_t";
const char *Long = "Long_t";
const char *Long64 = "Long64_t";
const char *Float = "Float_t";
const char *Double = "Double_t";
const char *UChar = "UChar_t";
const char *UShort = "UShort_t";
const char *UInt = "UInt_t";
const char *ULong = "ULong_t";
const char *ULong64 = "ULong64_t";
const char *String = "string";
const char *CharStar = "CharStar";
};

TString TXMLSetup::fgNameSpaceBase = "http://root.cern.ch/root/htmldoc/";

////////////////////////////////////////////////////////////////////////////////
/// return default value for XML setup

TString TXMLSetup::DefaultXmlSetup()
{
   return TString("2xoo");
}

////////////////////////////////////////////////////////////////////////////////
/// set namespace base

void TXMLSetup::SetNameSpaceBase(const char *namespacebase)
{
   fgNameSpaceBase = namespacebase;
}

////////////////////////////////////////////////////////////////////////////////
/// creates TXMLSetup object getting values from string

TXMLSetup::TXMLSetup(const char *opt)
{
   ReadSetupFromStr(opt);
}

////////////////////////////////////////////////////////////////////////////////
/// copy constructor of TXMLSetup class

TXMLSetup::TXMLSetup(const TXMLSetup &src)
   : fXmlLayout(src.fXmlLayout), fStoreStreamerInfos(src.fStoreStreamerInfos), fUseDtd(src.fUseDtd),
     fUseNamespaces(src.fUseNamespaces)
{
}

////////////////////////////////////////////////////////////////////////////////
/// assign operator

TXMLSetup &TXMLSetup::operator=(const TXMLSetup &rhs)
{
   fXmlLayout = rhs.fXmlLayout;
   fStoreStreamerInfos = rhs.fStoreStreamerInfos;
   fUseDtd = rhs.fUseDtd;
   fUseNamespaces = rhs.fUseNamespaces;
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// return setup values as string

TString TXMLSetup::GetSetupAsString()
{
   char setupstr[10] = "2xxx";

   setupstr[0] = char(48 + fXmlLayout);
   setupstr[1] = fStoreStreamerInfos ? 'x' : 'o';
   setupstr[2] = fUseDtd ? 'x' : 'o';
   setupstr[3] = fUseNamespaces ? 'x' : 'o';

   return TString(setupstr);
}

////////////////////////////////////////////////////////////////////////////////
/// checks if string is valid setup

Bool_t TXMLSetup::IsValidXmlSetup(const char *setupstr)
{
   if (!setupstr || (strlen(setupstr) != 4))
      return kFALSE;
   TString str = setupstr;
   str.ToLower();
   if ((str[0] < 48) || (str[0] > 53))
      return kFALSE;
   for (int n = 1; n < 4; n++)
      if ((str[n] != 'o') && (str[n] != 'x'))
         return kFALSE;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// get values from string

Bool_t TXMLSetup::ReadSetupFromStr(const char *setupstr)
{
   if (!setupstr || (strlen(setupstr) < 4))
      return kFALSE;
   Int_t lay = EXMLLayout(setupstr[0] - 48);
   if (lay == kGeneralized)
      fXmlLayout = kGeneralized;
   else
      fXmlLayout = kSpecialized;

   fStoreStreamerInfos = setupstr[1] == 'x';
   fUseDtd = kFALSE;
   fUseNamespaces = setupstr[3] == 'x';
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// show setup values

void TXMLSetup::PrintSetup()
{
   std::cout << " *** Setup printout ***" << std::endl;
   std::cout << "Attribute mode = " << fXmlLayout << std::endl;
   std::cout << "Store streamer infos = " << (fStoreStreamerInfos ? "true" : "false") << std::endl;
   std::cout << "Use dtd = " << (fUseDtd ? "true" : "false") << std::endl;
   std::cout << "Use name spaces = " << (fUseNamespaces ? "true" : "false") << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// convert class name to exclude any special symbols like ':', '<' '>' ',' and spaces

const char *TXMLSetup::XmlConvertClassName(const char *clname)
{
   fStrBuf = clname;
   fStrBuf.ReplaceAll("<", "_");
   fStrBuf.ReplaceAll(">", "_");
   fStrBuf.ReplaceAll(",", "_");
   fStrBuf.ReplaceAll(" ", "_");
   fStrBuf.ReplaceAll(":", "_");
   return fStrBuf.Data();
}

////////////////////////////////////////////////////////////////////////////////
/// produce string which used as reference in class namespace definition

const char *TXMLSetup::XmlClassNameSpaceRef(const TClass *cl)
{
   TString clname = XmlConvertClassName(cl->GetName());
   fStrBuf = fgNameSpaceBase;
   fStrBuf += clname;
   if (fgNameSpaceBase == "http://root.cern.ch/root/htmldoc/")
      fStrBuf += ".html";
   return fStrBuf.Data();
}

////////////////////////////////////////////////////////////////////////////////
///  return converted name for TStreamerElement

const char *TXMLSetup::XmlGetElementName(const TStreamerElement *el)
{
   if (!el)
      return nullptr;
   if (!el->InheritsFrom(TStreamerSTL::Class()))
      return el->GetName();
   if (strcmp(el->GetName(), el->GetClassPointer()->GetName()) != 0)
      return el->GetName();
   return XmlConvertClassName(el->GetName());
}

////////////////////////////////////////////////////////////////////////////////
/// get item name for given element

const char *TXMLSetup::GetElItemName(TStreamerElement *el)
{
   if (!el)
      return nullptr;
   fStrBuf = el->GetName();
   fStrBuf += "_item";
   return fStrBuf.Data();
}

////////////////////////////////////////////////////////////////////////////////
/// define class for the converted class name, where
/// special symbols were replaced by '_'

TClass *TXMLSetup::XmlDefineClass(const char *xmlClassName)
{
   if (strchr(xmlClassName, '_') == 0)
      return TClass::GetClass(xmlClassName);

   TIter iter(gROOT->GetListOfClasses());
   TClass *cl = nullptr;
   while ((cl = (TClass *)iter()) != nullptr) {
      const char *name = XmlConvertClassName(cl->GetName());
      if (strcmp(xmlClassName, name) == 0)
         return cl;
   }
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// converts string to integer.
/// if error, returns default value

Int_t TXMLSetup::AtoI(const char *sbuf, Int_t def, const char *errinfo)
{
   if (sbuf)
      return atoi(sbuf);
   if (errinfo)
      std::cerr << "<Error in TXMLSetup::AtoI>" << errinfo << " not valid integer: sbuf <NULL>" << std::endl;
   return def;
}
