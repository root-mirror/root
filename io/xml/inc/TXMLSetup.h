// @(#)root/xml:$Id$
// Author: Sergey Linev  10.05.2004

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TXMLSetup
#define ROOT_TXMLSetup

#include "Rtypes.h"
#include "TString.h"

#ifdef Bool
#undef Bool
#endif
#ifdef True
#undef True
#endif
#ifdef False
#undef False
#endif

namespace xmlio {
extern const char *Root;
extern const char *Setup;
extern const char *ClassVersion;
extern const char *IOVersion;
extern const char *OnlyVersion;
extern const char *Ptr;
extern const char *Ref;
extern const char *Null;
extern const char *IdBase;
extern const char *Size;
extern const char *Xmlobject;
extern const char *Xmlkey;
extern const char *Cycle;
extern const char *XmlBlock;
extern const char *Zip;
extern const char *Object;
extern const char *ObjClass;
extern const char *Class;
extern const char *Member;
extern const char *Item;
extern const char *Name;
extern const char *Title;
extern const char *CreateTm;
extern const char *ModifyTm;
extern const char *ObjectUUID;
extern const char *Type;
extern const char *Value;
extern const char *v;
extern const char *cnt;
extern const char *True;
extern const char *False;
extern const char *SInfos;

extern const char *Array;
extern const char *Bool;
extern const char *Char;
extern const char *Short;
extern const char *Int;
extern const char *Long;
extern const char *Long64;
extern const char *Float;
extern const char *Double;
extern const char *UChar;
extern const char *UShort;
extern const char *UInt;
extern const char *ULong;
extern const char *ULong64;
extern const char *String;
extern const char *CharStar;
}

class TStreamerElement;

class TXMLSetup {
public:
   enum EXMLLayout { kSpecialized = 2, kGeneralized = 3 };

   TXMLSetup() = default;
   TXMLSetup(const char *opt);
   TXMLSetup(const TXMLSetup &src);
   virtual ~TXMLSetup() = default;

   TXMLSetup &operator=(const TXMLSetup &rhs);

   TString GetSetupAsString();

   void PrintSetup();

   EXMLLayout GetXmlLayout() const { return fXmlLayout; }
   Bool_t IsStoreStreamerInfos() const { return fStoreStreamerInfos; }
   Bool_t IsUseDtd() const { return fUseDtd; }
   Bool_t IsUseNamespaces() const { return fUseNamespaces; }

   virtual void SetXmlLayout(EXMLLayout layout) { fXmlLayout = layout; }
   virtual void SetStoreStreamerInfos(Bool_t iConvert = kTRUE) { fStoreStreamerInfos = iConvert; }
   virtual void SetUsedDtd(Bool_t use = kTRUE) { fUseDtd = use; }
   virtual void SetUseNamespaces(Bool_t iUseNamespaces = kTRUE) { fUseNamespaces = iUseNamespaces; }

   const char *XmlConvertClassName(const char *name);
   const char *XmlClassNameSpaceRef(const TClass *cl);
   const char *XmlGetElementName(const TStreamerElement *el);

   Int_t GetNextRefCounter() { return fRefCounter++; }

   static TString DefaultXmlSetup();
   static void SetNameSpaceBase(const char *namespacebase);

protected:
   TClass *XmlDefineClass(const char *xmlClassName);
   const char *GetElItemName(TStreamerElement *el);

   Bool_t IsValidXmlSetup(const char *setupstr);
   Bool_t ReadSetupFromStr(const char *setupstr);

   Int_t AtoI(const char *sbuf, Int_t def = 0, const char *errinfo = nullptr);

   EXMLLayout fXmlLayout{kSpecialized};
   Bool_t fStoreStreamerInfos{kTRUE};
   Bool_t fUseDtd{kFALSE};
   Bool_t fUseNamespaces{kFALSE};

   Int_t fRefCounter{0}; //!  counter , used to build id of xml references

   TString fStrBuf; //!  buffer, used in XmlDefineClass() function

   static TString fgNameSpaceBase;

   ClassDef(TXMLSetup, 1) // settings to be stored in XML files
};

#endif
