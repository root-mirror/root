// $Id$
// Author: Sergey Linev  4.03.2014

/*************************************************************************
 * Copyright (C) 1995-2004, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TBufferJSON
#define ROOT_TBufferJSON

#include "TBufferText.h"
#include "TString.h"

#include <deque>

class TVirtualStreamerInfo;
class TStreamerInfo;
class TStreamerElement;
class TMemberStreamer;
class TDataMember;
class TJSONStackObj;

class TBufferJSON : public TBufferText {

public:
   TBufferJSON(TBuffer::EMode mode = TBuffer::kWrite);
   virtual ~TBufferJSON();

   void SetCompact(int level);

   static TString ConvertToJSON(const TObject *obj, Int_t compact = 0, const char *member_name = nullptr);
   static TString
   ConvertToJSON(const void *obj, const TClass *cl, Int_t compact = 0, const char *member_name = nullptr);
   static TString ConvertToJSON(const void *obj, TDataMember *member, Int_t compact = 0, Int_t arraylen = -1);

   static Int_t ExportToFile(const char *filename, const TObject *obj, const char *option = nullptr);
   static Int_t ExportToFile(const char *filename, const void *obj, const TClass *cl, const char *option = nullptr);

   static TObject *ConvertFromJSON(const char *str);
   static void *ConvertFromJSONAny(const char *str, TClass **cl = nullptr);

   template <class T>
   static TString ToJSON(const T *obj, Int_t compact = 0, const char *member_name = nullptr)
   {
      return ConvertToJSON(obj, TBuffer::GetClass(typeid(T)), compact, member_name);
   }

   template <class T>
   static Bool_t FromJSON(T *&obj, const char *json)
   {
      if (obj)
         return kFALSE;
      obj = (T *)ConvertFromJSONChecked(json, TBuffer::GetClass(typeid(T)));
      return obj != nullptr;
   }

   // suppress class writing/reading

   virtual TClass *ReadClass(const TClass *cl = nullptr, UInt_t *objTag = nullptr);
   virtual void WriteClass(const TClass *cl);

   // redefined virtual functions of TBuffer

   virtual Version_t ReadVersion(UInt_t *start = nullptr, UInt_t *bcnt = nullptr, const TClass *cl = nullptr);
   virtual UInt_t WriteVersion(const TClass *cl, Bool_t useBcnt = kFALSE);

   virtual void *ReadObjectAny(const TClass *clCast);
   virtual void SkipObjectAny();

   // these methods used in streamer info to indicate currently streamed element,
   virtual void IncrementLevel(TVirtualStreamerInfo *);
   virtual void SetStreamerElementNumber(TStreamerElement *elem, Int_t comp_type);
   virtual void DecrementLevel(TVirtualStreamerInfo *);

   virtual void ClassBegin(const TClass *, Version_t = -1);
   virtual void ClassEnd(const TClass *);
   virtual void ClassMember(const char *name, const char *typeName = nullptr, Int_t arrsize1 = -1, Int_t arrsize2 = -1);

   virtual Int_t ReadArray(Bool_t *&b);
   virtual Int_t ReadArray(Char_t *&c);
   virtual Int_t ReadArray(UChar_t *&c);
   virtual Int_t ReadArray(Short_t *&h);
   virtual Int_t ReadArray(UShort_t *&h);
   virtual Int_t ReadArray(Int_t *&i);
   virtual Int_t ReadArray(UInt_t *&i);
   virtual Int_t ReadArray(Long_t *&l);
   virtual Int_t ReadArray(ULong_t *&l);
   virtual Int_t ReadArray(Long64_t *&l);
   virtual Int_t ReadArray(ULong64_t *&l);
   virtual Int_t ReadArray(Float_t *&f);
   virtual Int_t ReadArray(Double_t *&d);

   virtual Int_t ReadStaticArray(Bool_t *b);
   virtual Int_t ReadStaticArray(Char_t *c);
   virtual Int_t ReadStaticArray(UChar_t *c);
   virtual Int_t ReadStaticArray(Short_t *h);
   virtual Int_t ReadStaticArray(UShort_t *h);
   virtual Int_t ReadStaticArray(Int_t *i);
   virtual Int_t ReadStaticArray(UInt_t *i);
   virtual Int_t ReadStaticArray(Long_t *l);
   virtual Int_t ReadStaticArray(ULong_t *l);
   virtual Int_t ReadStaticArray(Long64_t *l);
   virtual Int_t ReadStaticArray(ULong64_t *l);
   virtual Int_t ReadStaticArray(Float_t *f);
   virtual Int_t ReadStaticArray(Double_t *d);

   virtual void ReadFastArray(Bool_t *b, Int_t n);
   virtual void ReadFastArray(Char_t *c, Int_t n);
   virtual void ReadFastArrayString(Char_t *c, Int_t n);
   virtual void ReadFastArray(UChar_t *c, Int_t n);
   virtual void ReadFastArray(Short_t *h, Int_t n);
   virtual void ReadFastArray(UShort_t *h, Int_t n);
   virtual void ReadFastArray(Int_t *i, Int_t n);
   virtual void ReadFastArray(UInt_t *i, Int_t n);
   virtual void ReadFastArray(Long_t *l, Int_t n);
   virtual void ReadFastArray(ULong_t *l, Int_t n);
   virtual void ReadFastArray(Long64_t *l, Int_t n);
   virtual void ReadFastArray(ULong64_t *l, Int_t n);
   virtual void ReadFastArray(Float_t *f, Int_t n);
   virtual void ReadFastArray(Double_t *d, Int_t n);
   virtual void ReadFastArray(void *start, const TClass *cl, Int_t n = 1, TMemberStreamer *s = nullptr,
                              const TClass *onFileClass = nullptr);
   virtual void ReadFastArray(void **startp, const TClass *cl, Int_t n = 1, Bool_t isPreAlloc = kFALSE,
                              TMemberStreamer *s = nullptr, const TClass *onFileClass = nullptr);

   virtual void WriteArray(const Bool_t *b, Int_t n);
   virtual void WriteArray(const Char_t *c, Int_t n);
   virtual void WriteArray(const UChar_t *c, Int_t n);
   virtual void WriteArray(const Short_t *h, Int_t n);
   virtual void WriteArray(const UShort_t *h, Int_t n);
   virtual void WriteArray(const Int_t *i, Int_t n);
   virtual void WriteArray(const UInt_t *i, Int_t n);
   virtual void WriteArray(const Long_t *l, Int_t n);
   virtual void WriteArray(const ULong_t *l, Int_t n);
   virtual void WriteArray(const Long64_t *l, Int_t n);
   virtual void WriteArray(const ULong64_t *l, Int_t n);
   virtual void WriteArray(const Float_t *f, Int_t n);
   virtual void WriteArray(const Double_t *d, Int_t n);

   virtual void WriteFastArray(const Bool_t *b, Int_t n);
   virtual void WriteFastArray(const Char_t *c, Int_t n);
   virtual void WriteFastArrayString(const Char_t *c, Int_t n);
   virtual void WriteFastArray(const UChar_t *c, Int_t n);
   virtual void WriteFastArray(const Short_t *h, Int_t n);
   virtual void WriteFastArray(const UShort_t *h, Int_t n);
   virtual void WriteFastArray(const Int_t *i, Int_t n);
   virtual void WriteFastArray(const UInt_t *i, Int_t n);
   virtual void WriteFastArray(const Long_t *l, Int_t n);
   virtual void WriteFastArray(const ULong_t *l, Int_t n);
   virtual void WriteFastArray(const Long64_t *l, Int_t n);
   virtual void WriteFastArray(const ULong64_t *l, Int_t n);
   virtual void WriteFastArray(const Float_t *f, Int_t n);
   virtual void WriteFastArray(const Double_t *d, Int_t n);
   virtual void WriteFastArray(void *start, const TClass *cl, Int_t n = 1, TMemberStreamer *s = nullptr);
   virtual Int_t WriteFastArray(void **startp, const TClass *cl, Int_t n = 1, Bool_t isPreAlloc = kFALSE,
                                TMemberStreamer *s = nullptr);

   virtual void StreamObject(void *obj, const TClass *cl, const TClass *onFileClass = nullptr);
   using TBufferText::StreamObject;

   virtual void ReadBool(Bool_t &b);
   virtual void ReadChar(Char_t &c);
   virtual void ReadUChar(UChar_t &c);
   virtual void ReadShort(Short_t &s);
   virtual void ReadUShort(UShort_t &s);
   virtual void ReadInt(Int_t &i);
   virtual void ReadUInt(UInt_t &i);
   virtual void ReadLong(Long_t &l);
   virtual void ReadULong(ULong_t &l);
   virtual void ReadLong64(Long64_t &l);
   virtual void ReadULong64(ULong64_t &l);
   virtual void ReadFloat(Float_t &f);
   virtual void ReadDouble(Double_t &d);
   virtual void ReadCharP(Char_t *c);
   virtual void ReadTString(TString &s);
   virtual void ReadStdString(std::string *s);
   using TBuffer::ReadStdString;
   virtual void ReadCharStar(char *&s);

   virtual void WriteBool(Bool_t b);
   virtual void WriteChar(Char_t c);
   virtual void WriteUChar(UChar_t c);
   virtual void WriteShort(Short_t s);
   virtual void WriteUShort(UShort_t s);
   virtual void WriteInt(Int_t i);
   virtual void WriteUInt(UInt_t i);
   virtual void WriteLong(Long_t l);
   virtual void WriteULong(ULong_t l);
   virtual void WriteLong64(Long64_t l);
   virtual void WriteULong64(ULong64_t l);
   virtual void WriteFloat(Float_t f);
   virtual void WriteDouble(Double_t d);
   virtual void WriteCharP(const Char_t *c);
   virtual void WriteTString(const TString &s);
   virtual void WriteStdString(const std::string *s);
   using TBuffer::WriteStdString;
   virtual void WriteCharStar(char *s);

   virtual TVirtualStreamerInfo *GetInfo();

   // end of redefined virtual functions

   virtual void ReadBaseClass(void *start, TStreamerBase *elem);

protected:
   // redefined protected virtual functions

   virtual void WriteObjectClass(const void *actualObjStart, const TClass *actualClass, Bool_t cacheReuse);

   // end redefined protected virtual functions

   static void *ConvertFromJSONChecked(const char *str, const TClass *expectedClass);

   TString JsonWriteMember(const void *ptr, TDataMember *member, TClass *memberClass, Int_t arraylen);

   TJSONStackObj *PushStack(Int_t inclevel = 0, void *readnode = nullptr);
   TJSONStackObj *PopStack();
   TJSONStackObj *Stack() { return fStack.back(); }

   void WorkWithClass(TStreamerInfo *info, const TClass *cl = nullptr);
   void WorkWithElement(TStreamerElement *elem, Int_t);

   void JsonDisablePostprocessing();
   Int_t JsonSpecialClass(const TClass *cl) const;

   void JsonStartElement(const TStreamerElement *elem, const TClass *base_class = nullptr);

   void PerformPostProcessing(TJSONStackObj *stack, const TClass *obj_cl = nullptr);

   void JsonWriteBasic(Char_t value);
   void JsonWriteBasic(Short_t value);
   void JsonWriteBasic(Int_t value);
   void JsonWriteBasic(Long_t value);
   void JsonWriteBasic(Long64_t value);
   void JsonWriteBasic(Float_t value);
   void JsonWriteBasic(Double_t value);
   void JsonWriteBasic(Bool_t value);
   void JsonWriteBasic(UChar_t value);
   void JsonWriteBasic(UShort_t value);
   void JsonWriteBasic(UInt_t value);
   void JsonWriteBasic(ULong_t value);
   void JsonWriteBasic(ULong64_t value);

   void JsonWriteConstChar(const char *value, Int_t len = -1, const char * /*typname*/ = nullptr);

   void JsonWriteObject(const void *obj, const TClass *objClass, Bool_t check_map = kTRUE);

   void JsonWriteCollection(TCollection *obj, const TClass *objClass);

   void JsonReadCollection(TCollection *obj, const TClass *objClass);

   void JsonReadTObjectMembers(TObject *obj, void *node = nullptr);

   void *JsonReadObject(void *obj, const TClass *objClass = nullptr, TClass **readClass = nullptr);

   void AppendOutput(const char *line0, const char *line1 = nullptr);

   void JsonPushValue();

   template <typename T>
   R__ALWAYS_INLINE void JsonWriteArrayCompress(const T *vname, Int_t arrsize, const char *typname);

   template <typename T>
   R__ALWAYS_INLINE void JsonReadBasic(T &value);

   template <typename T>
   R__ALWAYS_INLINE Int_t JsonReadArray(T *value);

   template <typename T>
   R__ALWAYS_INLINE void JsonReadFastArray(T *arr, Int_t arrsize, bool asstring = false);

   template <typename T>
   R__ALWAYS_INLINE void JsonWriteFastArray(const T *arr, Int_t arrsize, const char *typname,
                                            void (TBufferJSON::*method)(const T *, Int_t, const char *));

   TString fOutBuffer;                 ///<!  main output buffer for json code
   TString *fOutput;                   ///<!  current output buffer for json code
   TString fValue;                     ///<!  buffer for current value
   unsigned fJsonrCnt;                 ///<!  counter for all objects, used for referencing
   std::deque<TJSONStackObj *> fStack; ///<!  hierarchy of currently streamed element
   Int_t fCompact;     ///<!  0 - no any compression, 1 - no spaces in the begin, 2 - no new lines, 3 - no spaces at all
   TString fSemicolon; ///<!  depending from compression level, " : " or ":"
   TString fArraySepar;    ///<!  depending from compression level, ", " or ","
   TString fNumericLocale; ///<!  stored value of setlocale(LC_NUMERIC), which should be recovered at the end

   ClassDef(TBufferJSON, 1) // a specialized TBuffer to only write objects into JSON format
};

#endif
