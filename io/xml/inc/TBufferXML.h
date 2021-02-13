// @(#)root/xml:$Id: d90d66e8fd2aa9daa4b05bcba9166aee1e2b2e7f $
// Author: Sergey Linev  10.05.2004

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TBufferXML
#define ROOT_TBufferXML

#include "Compression.h"
#include "TBufferText.h"
#include "TXMLSetup.h"
#include "TString.h"
#include "TXMLEngine.h"

#include <string>
#include <deque>
#include <memory>

class TExMap;
class TVirtualStreamerInfo;
class TStreamerInfo;
class TStreamerElement;
class TMemberStreamer;
class TXMLFile;
class TXMLStackObj;

class TBufferXML final : public TBufferText, public TXMLSetup {

   friend class TKeyXML;

public:
   TBufferXML(TBuffer::EMode mode);
   TBufferXML(TBuffer::EMode mode, TXMLFile *file);
   virtual ~TBufferXML();

   static TString ConvertToXML(const TObject *obj, Bool_t GenericLayout = kFALSE, Bool_t UseNamespaces = kFALSE);
   static TString
   ConvertToXML(const void *obj, const TClass *cl, Bool_t GenericLayout = kFALSE, Bool_t UseNamespaces = kFALSE);

   template <class T>
   static TString ToXML(const T *obj, Bool_t GenericLayout = kFALSE, Bool_t UseNamespaces = kFALSE)
   {
      return ConvertToXML(obj, TClass::GetClass<T>(), GenericLayout, UseNamespaces);
   }

   static TObject *ConvertFromXML(const char *str, Bool_t GenericLayout = kFALSE, Bool_t UseNamespaces = kFALSE);
   static void *ConvertFromXMLAny(const char *str, TClass **cl = nullptr, Bool_t GenericLayout = kFALSE,
                                  Bool_t UseNamespaces = kFALSE);

   template <class T>
   static Bool_t FromXML(T *&obj, const char *xml, Bool_t GenericLayout = kFALSE, Bool_t UseNamespaces = kFALSE)
   {
      if (obj)
         return kFALSE;
      obj = (T *)ConvertFromXMLChecked(xml, TClass::GetClass<T>(), GenericLayout, UseNamespaces);
      return obj != nullptr;
   }

   Int_t GetIOVersion() const { return fIOVersion; }
   void SetIOVersion(Int_t v) { fIOVersion = v; }

   // suppress class writing/reading

   TClass *ReadClass(const TClass *cl = nullptr, UInt_t *objTag = nullptr) final;
   void WriteClass(const TClass *cl) final;

   // redefined virtual functions of TBuffer

   Version_t ReadVersion(UInt_t *start = nullptr, UInt_t *bcnt = nullptr, const TClass *cl = nullptr) final;
   UInt_t WriteVersion(const TClass *cl, Bool_t useBcnt = kFALSE) final;

   void *ReadObjectAny(const TClass *clCast) final;
   void SkipObjectAny() final;

   void IncrementLevel(TVirtualStreamerInfo *) final;
   void SetStreamerElementNumber(TStreamerElement *elem, Int_t comp_type) final;
   void DecrementLevel(TVirtualStreamerInfo *) final;

   void ClassBegin(const TClass *, Version_t = -1) final;
   void ClassEnd(const TClass *) final;
   void ClassMember(const char *name, const char *typeName = nullptr, Int_t arrsize1 = -1, Int_t arrsize2 = -1) final;

   Int_t ReadArray(Bool_t *&b) final;
   Int_t ReadArray(Char_t *&c) final;
   Int_t ReadArray(UChar_t *&c) final;
   Int_t ReadArray(Short_t *&h) final;
   Int_t ReadArray(UShort_t *&h) final;
   Int_t ReadArray(Int_t *&i) final;
   Int_t ReadArray(UInt_t *&i) final;
   Int_t ReadArray(Long_t *&l) final;
   Int_t ReadArray(ULong_t *&l) final;
   Int_t ReadArray(Long64_t *&l) final;
   Int_t ReadArray(ULong64_t *&l) final;
   Int_t ReadArray(Float_t *&f) final;
   Int_t ReadArray(Double_t *&d) final;

   Int_t ReadStaticArray(Bool_t *b) final;
   Int_t ReadStaticArray(Char_t *c) final;
   Int_t ReadStaticArray(UChar_t *c) final;
   Int_t ReadStaticArray(Short_t *h) final;
   Int_t ReadStaticArray(UShort_t *h) final;
   Int_t ReadStaticArray(Int_t *i) final;
   Int_t ReadStaticArray(UInt_t *i) final;
   Int_t ReadStaticArray(Long_t *l) final;
   Int_t ReadStaticArray(ULong_t *l) final;
   Int_t ReadStaticArray(Long64_t *l) final;
   Int_t ReadStaticArray(ULong64_t *l) final;
   Int_t ReadStaticArray(Float_t *f) final;
   Int_t ReadStaticArray(Double_t *d) final;

   void ReadFastArray(Bool_t *b, Int_t n) final;
   void ReadFastArray(Char_t *c, Int_t n) final;
   void ReadFastArray(UChar_t *c, Int_t n) final;
   void ReadFastArray(Short_t *h, Int_t n) final;
   void ReadFastArray(UShort_t *h, Int_t n) final;
   void ReadFastArray(Int_t *i, Int_t n) final;
   void ReadFastArray(UInt_t *i, Int_t n) final;
   void ReadFastArray(Long_t *l, Int_t n) final;
   void ReadFastArray(ULong_t *l, Int_t n) final;
   void ReadFastArray(Long64_t *l, Int_t n) final;
   void ReadFastArray(ULong64_t *l, Int_t n) final;
   void ReadFastArray(Float_t *f, Int_t n) final;
   void ReadFastArray(Double_t *d, Int_t n) final;
   void ReadFastArrayString(Char_t *c, Int_t n) final;
   void ReadFastArray(void *start, const TClass *cl, Int_t n = 1, TMemberStreamer *s = nullptr,
                      const TClass *onFileClass = nullptr) final;
   void ReadFastArray(void **startp, const TClass *cl, Int_t n = 1, Bool_t isPreAlloc = kFALSE,
                      TMemberStreamer *s = nullptr, const TClass *onFileClass = nullptr) final;

   void WriteArray(const Bool_t *b, Int_t n) final;
   void WriteArray(const Char_t *c, Int_t n) final;
   void WriteArray(const UChar_t *c, Int_t n) final;
   void WriteArray(const Short_t *h, Int_t n) final;
   void WriteArray(const UShort_t *h, Int_t n) final;
   void WriteArray(const Int_t *i, Int_t n) final;
   void WriteArray(const UInt_t *i, Int_t n) final;
   void WriteArray(const Long_t *l, Int_t n) final;
   void WriteArray(const ULong_t *l, Int_t n) final;
   void WriteArray(const Long64_t *l, Int_t n) final;
   void WriteArray(const ULong64_t *l, Int_t n) final;
   void WriteArray(const Float_t *f, Int_t n) final;
   void WriteArray(const Double_t *d, Int_t n) final;

   void WriteFastArray(const Bool_t *b, Int_t n) final;
   void WriteFastArray(const Char_t *c, Int_t n) final;
   void WriteFastArray(const UChar_t *c, Int_t n) final;
   void WriteFastArray(const Short_t *h, Int_t n) final;
   void WriteFastArray(const UShort_t *h, Int_t n) final;
   void WriteFastArray(const Int_t *i, Int_t n) final;
   void WriteFastArray(const UInt_t *i, Int_t n) final;
   void WriteFastArray(const Long_t *l, Int_t n) final;
   void WriteFastArray(const ULong_t *l, Int_t n) final;
   void WriteFastArray(const Long64_t *l, Int_t n) final;
   void WriteFastArray(const ULong64_t *l, Int_t n) final;
   void WriteFastArray(const Float_t *f, Int_t n) final;
   void WriteFastArray(const Double_t *d, Int_t n) final;
   void WriteFastArrayString(const Char_t *c, Int_t n) final;
   void WriteFastArray(void *start, const TClass *cl, Int_t n = 1, TMemberStreamer *s = nullptr) final;
   Int_t WriteFastArray(void **startp, const TClass *cl, Int_t n = 1, Bool_t isPreAlloc = kFALSE,
                        TMemberStreamer *s = nullptr) final;

   void StreamObject(void *obj, const TClass *cl, const TClass *onFileClass = nullptr) final;
   using TBufferText::StreamObject;

   void ReadBool(Bool_t &b) final;
   void ReadChar(Char_t &c) final;
   void ReadUChar(UChar_t &c) final;
   void ReadShort(Short_t &s) final;
   void ReadUShort(UShort_t &s) final;
   void ReadInt(Int_t &i) final;
   void ReadUInt(UInt_t &i) final;
   void ReadLong(Long_t &l) final;
   void ReadULong(ULong_t &l) final;
   void ReadLong64(Long64_t &l) final;
   void ReadULong64(ULong64_t &l) final;
   void ReadFloat(Float_t &f) final;
   void ReadDouble(Double_t &d) final;
   void ReadCharP(Char_t *c) final;
   void ReadTString(TString &s) final;
   void ReadStdString(std::string *s) final;
   using TBuffer::ReadStdString;
   void ReadCharStar(char *&s) final;

   void WriteBool(Bool_t b) final;
   void WriteChar(Char_t c) final;
   void WriteUChar(UChar_t c) final;
   void WriteShort(Short_t s) final;
   void WriteUShort(UShort_t s) final;
   void WriteInt(Int_t i) final;
   void WriteUInt(UInt_t i) final;
   void WriteLong(Long_t l) final;
   void WriteULong(ULong_t l) final;
   void WriteLong64(Long64_t l) final;
   void WriteULong64(ULong64_t l) final;
   void WriteFloat(Float_t f) final;
   void WriteDouble(Double_t d) final;
   void WriteCharP(const Char_t *c) final;
   void WriteTString(const TString &s) final;
   void WriteStdString(const std::string *s) final;
   using TBuffer::WriteStdString;
   void WriteCharStar(char *s) final;

   TVirtualStreamerInfo *GetInfo() final;

protected:
   // redefined protected virtual functions

   void WriteObjectClass(const void *actualObjStart, const TClass *actualClass, Bool_t cacheReuse) final;

   // end redefined protected virtual functions

   static void *ConvertFromXMLChecked(const char *xml, const TClass *expectedClass, Bool_t GenericLayout = kFALSE,
                                      Bool_t UseNamespaces = kFALSE);

   TXMLFile *XmlFile();

   Int_t GetCompressionAlgorithm() const;
   Int_t GetCompressionLevel() const;
   Int_t GetCompressionSettings() const;
   void SetCompressionAlgorithm(Int_t algorithm = ROOT::RCompressionSetting::EAlgorithm::kUseGlobal);
   void SetCompressionLevel(Int_t level = ROOT::RCompressionSetting::ELevel::kUseMin);
   void SetCompressionSettings(Int_t settings = ROOT::RCompressionSetting::EDefaults::kUseCompiledDefault);
   void SetXML(TXMLEngine *xml) { fXML = xml; }

   void XmlWriteBlock(XMLNodePointer_t node);
   XMLNodePointer_t XmlWriteAny(const void *obj, const TClass *cl);

   void XmlReadBlock(XMLNodePointer_t node);
   void *XmlReadAny(XMLNodePointer_t node, void *obj, TClass **cl);

   TXMLStackObj *PushStack(XMLNodePointer_t current, Bool_t simple = kFALSE);
   TXMLStackObj *PopStack();
   void ShiftStack(const char *info = nullptr);

   XMLNodePointer_t StackNode();
   TXMLStackObj *Stack(UInt_t depth = 0)
   {
      return (depth < fStack.size()) ? (depth ? fStack[fStack.size() - depth - 1].get() : fStack.back().get()) : nullptr;
   }

   void WorkWithClass(TStreamerInfo *info, const TClass *cl = nullptr);
   void WorkWithElement(TStreamerElement *elem, Int_t comp_type);
   Bool_t VerifyNode(XMLNodePointer_t node, const char *name, const char *errinfo = nullptr);
   Bool_t VerifyStackNode(const char *name, const char *errinfo = nullptr);

   Bool_t VerifyAttr(XMLNodePointer_t node, const char *name, const char *value, const char *errinfo = nullptr);
   Bool_t VerifyStackAttr(const char *name, const char *value, const char *errinfo = nullptr);

   Bool_t ProcessPointer(const void *ptr, XMLNodePointer_t node);
   Bool_t ExtractPointer(XMLNodePointer_t node, void *&ptr, TClass *&cl);
   void ExtractReference(XMLNodePointer_t node, const void *ptr, const TClass *cl);

   XMLNodePointer_t CreateItemNode(const char *name);
   Bool_t VerifyItemNode(const char *name, const char *errinfo = nullptr);

   void CreateElemNode(const TStreamerElement *elem);
   Bool_t VerifyElemNode(const TStreamerElement *elem);

   void PerformPreProcessing(const TStreamerElement *elem, XMLNodePointer_t elemnode);
   void PerformPostProcessing();

   XMLNodePointer_t XmlWriteBasic(Char_t value);
   XMLNodePointer_t XmlWriteBasic(Short_t value);
   XMLNodePointer_t XmlWriteBasic(Int_t value);
   XMLNodePointer_t XmlWriteBasic(Long_t value);
   XMLNodePointer_t XmlWriteBasic(Long64_t value);
   XMLNodePointer_t XmlWriteBasic(Float_t value);
   XMLNodePointer_t XmlWriteBasic(Double_t value);
   XMLNodePointer_t XmlWriteBasic(Bool_t value);
   XMLNodePointer_t XmlWriteBasic(UChar_t value);
   XMLNodePointer_t XmlWriteBasic(UShort_t value);
   XMLNodePointer_t XmlWriteBasic(UInt_t value);
   XMLNodePointer_t XmlWriteBasic(ULong_t value);
   XMLNodePointer_t XmlWriteBasic(ULong64_t value);
   XMLNodePointer_t XmlWriteValue(const char *value, const char *name);

   void XmlReadBasic(Char_t &value);
   void XmlReadBasic(Short_t &value);
   void XmlReadBasic(Int_t &value);
   void XmlReadBasic(Long_t &value);
   void XmlReadBasic(Long64_t &value);
   void XmlReadBasic(Float_t &value);
   void XmlReadBasic(Double_t &value);
   void XmlReadBasic(Bool_t &value);
   void XmlReadBasic(UChar_t &value);
   void XmlReadBasic(UShort_t &value);
   void XmlReadBasic(UInt_t &value);
   void XmlReadBasic(ULong_t &value);
   void XmlReadBasic(ULong64_t &value);
   const char *XmlReadValue(const char *name);

   template <typename T>
   R__ALWAYS_INLINE void XmlReadArrayContent(T *arr, Int_t arrsize);

   template <typename T>
   R__ALWAYS_INLINE Int_t XmlReadArray(T *&arr, bool is_static = false);

   template <typename T>
   R__ALWAYS_INLINE void XmlReadFastArray(T *arr, Int_t n);

   template <typename T>
   R__ALWAYS_INLINE void XmlWriteArrayContent(const T *arr, Int_t arrsize);

   template <typename T>
   R__ALWAYS_INLINE void XmlWriteArray(const T *arr, Int_t arrsize);

   template <typename T>
   R__ALWAYS_INLINE void XmlWriteFastArray(const T *arr, Int_t n);

   XMLNodePointer_t XmlWriteObject(const void *obj, const TClass *objClass, Bool_t cacheReuse);
   void *XmlReadObject(void *obj, TClass **cl = nullptr);

   void BeforeIOoperation();
   void CheckVersionBuf();

   TXMLEngine *fXML{nullptr};                        ///<! instance of TXMLEngine for working with XML structures
   std::deque<std::unique_ptr<TXMLStackObj>> fStack; ///<! Stack of processed objects
   Version_t fVersionBuf{-111};                      ///<! Current version buffer
   TString fValueBuf;                                ///<! Current value buffer
   Int_t fErrorFlag{0};                              ///<! Error flag
   Bool_t fCanUseCompact{kFALSE};                    ///<! Flag indicate that basic type (like Int_t) can be placed in the same tag
   TClass *fExpectedBaseClass{nullptr};              ///<! Pointer to class, which should be stored as parent of current
   Int_t fCompressLevel{0};                          ///<! Compression level and algorithm
   Int_t fIOVersion{3};                              ///<! Indicates format of ROOT xml file

   ClassDefOverride(TBufferXML, 0); // a specialized TBuffer to read/write to XML files
};

//______________________________________________________________________________
inline Int_t TBufferXML::GetCompressionAlgorithm() const
{
   return (fCompressLevel < 0) ? -1 : fCompressLevel / 100;
}

//______________________________________________________________________________
inline Int_t TBufferXML::GetCompressionLevel() const
{
   return (fCompressLevel < 0) ? -1 : fCompressLevel % 100;
}

//______________________________________________________________________________
inline Int_t TBufferXML::GetCompressionSettings() const
{
   return (fCompressLevel < 0) ? -1 : fCompressLevel;
}

#endif
