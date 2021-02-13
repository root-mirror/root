// @(#)root/base:$Id$
// Author: Rene Brun, Philippe Canal, Fons Rademakers   04/05/96

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TBuffer
#define ROOT_TBuffer


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TBuffer                                                              //
//                                                                      //
// Buffer base class used for serializing objects.                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TDataType.h"
#include "TObject.h"
#include "TClass.h"
#include "Bytes.h"

#include <vector>
#include <string>

class TVirtualStreamerInfo;
class TStreamerElement;
class TString;
class TProcessID;
class TClonesArray;
class TRefTable;
class TVirtualArray;
namespace TStreamerInfoActions {
   class TActionSequence;
}

class TBuffer : public TObject {

protected:
   typedef std::vector<TVirtualArray*> CacheList_t;

   Bool_t           fMode;          //Read or write mode
   Int_t            fVersion;       //Buffer format version
   Int_t            fBufSize;       //Size of buffer
   char            *fBuffer;        //Buffer used to store objects
   char            *fBufCur;        //Current position in buffer
   char            *fBufMax;        //End of buffer
   TObject         *fParent;        //Pointer to parent object owning this buffer
   ReAllocCharFun_t fReAllocFunc;   //! Realloc function to be used when extending the buffer.
   CacheList_t      fCacheStack;    //Stack of pointers to the cache where to temporarily store the value of 'missing' data members

   // Default ctor
   TBuffer() : TObject(), fMode(0), fVersion(0), fBufSize(0), fBuffer(0),
     fBufCur(0), fBufMax(0), fParent(0), fReAllocFunc(0), fCacheStack(0,(TVirtualArray*)0) {}

   // TBuffer objects cannot be copied or assigned
   TBuffer(const TBuffer &) = delete;
   void operator=(const TBuffer &) = delete;

   Int_t Read(const char *name) { return TObject::Read(name); }
   Int_t Write(const char *name, Int_t opt, Int_t bufs)
                              { return TObject::Write(name, opt, bufs); }
   Int_t Write(const char *name, Int_t opt, Int_t bufs) const
                              { return TObject::Write(name, opt, bufs); }

public:
   enum EMode { kRead = 0, kWrite = 1 };
   enum EStatusBits {
     kIsOwner = BIT(16), //if set TBuffer owns fBuffer
     kCannotHandleMemberWiseStreaming = BIT(17) //if set TClonesArray should not use member wise streaming
   };
   enum { kInitialSize = 1024, kMinimalSize = 128 };

   TBuffer(EMode mode);
   TBuffer(EMode mode, Int_t bufsiz);
   TBuffer(EMode mode, Int_t bufsiz, void *buf, Bool_t adopt = kTRUE, ReAllocCharFun_t reallocfunc = 0);
   virtual ~TBuffer();

   Int_t    GetBufferVersion() const { return fVersion; }
   Bool_t   IsReading() const { return (fMode & kWrite) == 0; }
   Bool_t   IsWriting() const { return (fMode & kWrite) != 0; }
   void     SetReadMode();
   void     SetWriteMode();
   void     SetBuffer(void *buf, UInt_t bufsiz = 0, Bool_t adopt = kTRUE, ReAllocCharFun_t reallocfunc = 0);
   ReAllocCharFun_t GetReAllocFunc() const;
   void     SetReAllocFunc(ReAllocCharFun_t reallocfunc = 0);
   void     SetBufferOffset(Int_t offset = 0) { fBufCur = fBuffer+offset; }
   void     SetParent(TObject *parent);
   TObject *GetParent()  const;
   char    *Buffer()     const { return fBuffer; }
   char    *GetCurrent() const { return fBufCur; }
   Int_t    BufferSize() const { return fBufSize; }
   void     DetachBuffer() { fBuffer = 0; }
   Int_t    Length()     const { return (Int_t)(fBufCur - fBuffer); }
   void     Expand(Int_t newsize, Bool_t copy = kTRUE);  // expand buffer to newsize
   void     AutoExpand(Int_t size_needed);  // expand buffer to newsize
   Bool_t   ByteSwapBuffer(Long64_t n, EDataType type);  // Byte-swap N primitive-elements in the buffer

   virtual Bool_t     CheckObject(const TObject *obj) = 0;
   virtual Bool_t     CheckObject(const void *obj, const TClass *ptrClass) = 0;

   virtual Int_t      ReadBuf(void *buf, Int_t max) = 0;
   virtual void       WriteBuf(const void *buf, Int_t max) = 0;

   virtual char      *ReadString(char *s, Int_t max) = 0;
   virtual void       WriteString(const char *s) = 0;

   virtual Int_t      GetVersionOwner() const  = 0;
   virtual Int_t      GetMapCount() const  = 0;
   virtual void       GetMappedObject(UInt_t tag, void* &ptr, TClass* &ClassPtr) const = 0;
   virtual void       MapObject(const TObject *obj, UInt_t offset = 1) = 0;
   virtual void       MapObject(const void *obj, const TClass *cl, UInt_t offset = 1) = 0;
   virtual void       Reset() = 0;
   virtual void       InitMap() = 0;
   virtual void       ResetMap() = 0;
   virtual void       SetReadParam(Int_t mapsize) = 0;
   virtual void       SetWriteParam(Int_t mapsize) = 0;

   virtual Int_t      CheckByteCount(UInt_t startpos, UInt_t bcnt, const TClass *clss) = 0;
   virtual Int_t      CheckByteCount(UInt_t startpos, UInt_t bcnt, const char *classname) = 0;
   virtual void       SetByteCount(UInt_t cntpos, Bool_t packInVersion = kFALSE)= 0;

   virtual void       SkipVersion(const TClass *cl = 0) = 0;
   virtual Version_t  ReadVersion(UInt_t *start = 0, UInt_t *bcnt = 0, const TClass *cl = 0) = 0;
   virtual Version_t  ReadVersionNoCheckSum(UInt_t *start = 0, UInt_t *bcnt = 0) = 0;
   virtual Version_t  ReadVersionForMemberWise(const TClass *cl = 0) = 0;
   virtual UInt_t     WriteVersion(const TClass *cl, Bool_t useBcnt = kFALSE) = 0;
   virtual UInt_t     WriteVersionMemberWise(const TClass *cl, Bool_t useBcnt = kFALSE) = 0;

   virtual void      *ReadObjectAny(const TClass* cast) = 0;
   virtual void       SkipObjectAny() = 0;

   virtual void       TagStreamerInfo(TVirtualStreamerInfo* info) = 0;
   virtual void       IncrementLevel(TVirtualStreamerInfo* info) = 0;
   virtual void       SetStreamerElementNumber(TStreamerElement *elem, Int_t comp_type) = 0;
   virtual void       DecrementLevel(TVirtualStreamerInfo*) = 0;

   virtual void       ClassBegin(const TClass*, Version_t = -1) = 0;
   virtual void       ClassEnd(const TClass*) = 0;
   virtual void       ClassMember(const char*, const char* = 0, Int_t = -1, Int_t = -1) = 0;
   virtual TVirtualStreamerInfo *GetInfo() = 0;

   virtual TVirtualArray *PeekDataCache() const;
   virtual TVirtualArray *PopDataCache();
   virtual void           PushDataCache(TVirtualArray *);

   virtual TClass    *ReadClass(const TClass *cl = 0, UInt_t *objTag = 0) = 0;
   virtual void       WriteClass(const TClass *cl) = 0;

   virtual TObject   *ReadObject(const TClass *cl) = 0;
   virtual void       WriteObject(const TObject *obj, Bool_t cacheReuse) = 0;

   template <class T> Int_t WriteObject(const T *objptr, Bool_t cacheReuse = kTRUE);

   virtual Int_t      WriteObjectAny(const void *obj, const TClass *ptrClass, Bool_t cacheReuse = kTRUE) = 0;

   virtual UShort_t   GetPidOffset() const  = 0;
   virtual void       SetPidOffset(UShort_t offset) = 0;
   virtual Int_t      GetBufferDisplacement() const  = 0;
   virtual void       SetBufferDisplacement() = 0;
   virtual void       SetBufferDisplacement(Int_t skipped) = 0;

   // basic types and arrays of basic types
   virtual   void     ReadFloat16 (Float_t *f, TStreamerElement *ele=0) = 0;
   virtual   void     WriteFloat16(Float_t *f, TStreamerElement *ele=0) = 0;
   virtual   void     ReadDouble32 (Double_t *d, TStreamerElement *ele=0) = 0;
   virtual   void     WriteDouble32(Double_t *d, TStreamerElement *ele=0) = 0;
   virtual   void     ReadWithFactor(Float_t *ptr, Double_t factor, Double_t minvalue) = 0;
   virtual   void     ReadWithNbits(Float_t *ptr, Int_t nbits) = 0;
   virtual   void     ReadWithFactor(Double_t *ptr, Double_t factor, Double_t minvalue) = 0;
   virtual   void     ReadWithNbits(Double_t *ptr, Int_t nbits) = 0;

   virtual   Int_t    ReadArray(Bool_t    *&b) = 0;
   virtual   Int_t    ReadArray(Char_t    *&c) = 0;
   virtual   Int_t    ReadArray(UChar_t   *&c) = 0;
   virtual   Int_t    ReadArray(Short_t   *&h) = 0;
   virtual   Int_t    ReadArray(UShort_t  *&h) = 0;
   virtual   Int_t    ReadArray(Int_t     *&i) = 0;
   virtual   Int_t    ReadArray(UInt_t    *&i) = 0;
   virtual   Int_t    ReadArray(Long_t    *&l) = 0;
   virtual   Int_t    ReadArray(ULong_t   *&l) = 0;
   virtual   Int_t    ReadArray(Long64_t  *&l) = 0;
   virtual   Int_t    ReadArray(ULong64_t *&l) = 0;
   virtual   Int_t    ReadArray(Float_t   *&f) = 0;
   virtual   Int_t    ReadArray(Double_t  *&d) = 0;
   virtual   Int_t    ReadArrayFloat16(Float_t *&f, TStreamerElement *ele=0) = 0;
   virtual   Int_t    ReadArrayDouble32(Double_t *&d, TStreamerElement *ele=0) = 0;

   virtual   Int_t    ReadStaticArray(Bool_t    *b) = 0;
   virtual   Int_t    ReadStaticArray(Char_t    *c) = 0;
   virtual   Int_t    ReadStaticArray(UChar_t   *c) = 0;
   virtual   Int_t    ReadStaticArray(Short_t   *h) = 0;
   virtual   Int_t    ReadStaticArray(UShort_t  *h) = 0;
   virtual   Int_t    ReadStaticArray(Int_t     *i) = 0;
   virtual   Int_t    ReadStaticArray(UInt_t    *i) = 0;
   virtual   Int_t    ReadStaticArray(Long_t    *l) = 0;
   virtual   Int_t    ReadStaticArray(ULong_t   *l) = 0;
   virtual   Int_t    ReadStaticArray(Long64_t  *l) = 0;
   virtual   Int_t    ReadStaticArray(ULong64_t *l) = 0;
   virtual   Int_t    ReadStaticArray(Float_t   *f) = 0;
   virtual   Int_t    ReadStaticArray(Double_t  *d) = 0;
   virtual   Int_t    ReadStaticArrayFloat16(Float_t  *f, TStreamerElement *ele=0) = 0;
   virtual   Int_t    ReadStaticArrayDouble32(Double_t  *d, TStreamerElement *ele=0) = 0;

   virtual   void     ReadFastArray(Bool_t    *b, Int_t n) = 0;
   virtual   void     ReadFastArray(Char_t    *c, Int_t n) = 0;
   virtual   void     ReadFastArrayString(Char_t *c, Int_t n) = 0;
   virtual   void     ReadFastArray(UChar_t   *c, Int_t n) = 0;
   virtual   void     ReadFastArray(Short_t   *h, Int_t n) = 0;
   virtual   void     ReadFastArray(UShort_t  *h, Int_t n) = 0;
   virtual   void     ReadFastArray(Int_t     *i, Int_t n) = 0;
   virtual   void     ReadFastArray(UInt_t    *i, Int_t n) = 0;
   virtual   void     ReadFastArray(Long_t    *l, Int_t n) = 0;
   virtual   void     ReadFastArray(ULong_t   *l, Int_t n) = 0;
   virtual   void     ReadFastArray(Long64_t  *l, Int_t n) = 0;
   virtual   void     ReadFastArray(ULong64_t *l, Int_t n) = 0;
   virtual   void     ReadFastArray(Float_t   *f, Int_t n) = 0;
   virtual   void     ReadFastArray(Double_t  *d, Int_t n) = 0;
   virtual   void     ReadFastArrayFloat16(Float_t  *f, Int_t n, TStreamerElement *ele=0) = 0;
   virtual   void     ReadFastArrayDouble32(Double_t  *d, Int_t n, TStreamerElement *ele=0) = 0;
   virtual   void     ReadFastArrayWithFactor(Float_t *ptr, Int_t n, Double_t factor, Double_t minvalue) = 0;
   virtual   void     ReadFastArrayWithNbits(Float_t *ptr, Int_t n, Int_t nbits) = 0;
   virtual   void     ReadFastArrayWithFactor(Double_t *ptr, Int_t n, Double_t factor, Double_t minvalue) = 0;
   virtual   void     ReadFastArrayWithNbits(Double_t *ptr, Int_t n, Int_t nbits) = 0;
   virtual   void     ReadFastArray(void  *start , const TClass *cl, Int_t n=1, TMemberStreamer *s=0, const TClass *onFileClass=0) = 0;
   virtual   void     ReadFastArray(void **startp, const TClass *cl, Int_t n=1, Bool_t isPreAlloc=kFALSE, TMemberStreamer *s=0, const TClass *onFileClass=0) = 0;

   virtual   void     WriteArray(const Bool_t    *b, Int_t n) = 0;
   virtual   void     WriteArray(const Char_t    *c, Int_t n) = 0;
   virtual   void     WriteArray(const UChar_t   *c, Int_t n) = 0;
   virtual   void     WriteArray(const Short_t   *h, Int_t n) = 0;
   virtual   void     WriteArray(const UShort_t  *h, Int_t n) = 0;
   virtual   void     WriteArray(const Int_t     *i, Int_t n) = 0;
   virtual   void     WriteArray(const UInt_t    *i, Int_t n) = 0;
   virtual   void     WriteArray(const Long_t    *l, Int_t n) = 0;
   virtual   void     WriteArray(const ULong_t   *l, Int_t n) = 0;
   virtual   void     WriteArray(const Long64_t  *l, Int_t n) = 0;
   virtual   void     WriteArray(const ULong64_t *l, Int_t n) = 0;
   virtual   void     WriteArray(const Float_t   *f, Int_t n) = 0;
   virtual   void     WriteArray(const Double_t  *d, Int_t n) = 0;
   virtual   void     WriteArrayFloat16(const Float_t  *f, Int_t n, TStreamerElement *ele=0) = 0;
   virtual   void     WriteArrayDouble32(const Double_t  *d, Int_t n, TStreamerElement *ele=0) = 0;

   virtual   void     WriteFastArray(const Bool_t    *b, Int_t n) = 0;
   virtual   void     WriteFastArray(const Char_t    *c, Int_t n) = 0;
   virtual   void     WriteFastArrayString(const Char_t    *c, Int_t n) = 0;
   virtual   void     WriteFastArray(const UChar_t   *c, Int_t n) = 0;
   virtual   void     WriteFastArray(const Short_t   *h, Int_t n) = 0;
   virtual   void     WriteFastArray(const UShort_t  *h, Int_t n) = 0;
   virtual   void     WriteFastArray(const Int_t     *i, Int_t n) = 0;
   virtual   void     WriteFastArray(const UInt_t    *i, Int_t n) = 0;
   virtual   void     WriteFastArray(const Long_t    *l, Int_t n) = 0;
   virtual   void     WriteFastArray(const ULong_t   *l, Int_t n) = 0;
   virtual   void     WriteFastArray(const Long64_t  *l, Int_t n) = 0;
   virtual   void     WriteFastArray(const ULong64_t *l, Int_t n) = 0;
   virtual   void     WriteFastArray(const Float_t   *f, Int_t n) = 0;
   virtual   void     WriteFastArray(const Double_t  *d, Int_t n) = 0;
   virtual   void     WriteFastArrayFloat16(const Float_t  *f, Int_t n, TStreamerElement *ele=0) = 0;
   virtual   void     WriteFastArrayDouble32(const Double_t  *d, Int_t n, TStreamerElement *ele=0) = 0;
   virtual   void     WriteFastArray(void  *start,  const TClass *cl, Int_t n=1, TMemberStreamer *s=0) = 0;
   virtual   Int_t    WriteFastArray(void **startp, const TClass *cl, Int_t n=1, Bool_t isPreAlloc=kFALSE, TMemberStreamer *s=0) = 0;

   virtual   void     StreamObject(void *obj, const std::type_info &typeinfo, const TClass* onFileClass = 0 ) = 0;
   virtual   void     StreamObject(void *obj, const char *className, const TClass* onFileClass = 0 ) = 0;
   virtual   void     StreamObject(void *obj, const TClass *cl, const TClass* onFileClass = 0 ) = 0;
   virtual   void     StreamObject(TObject *obj) = 0;

   virtual   void     ReadBool(Bool_t       &b) = 0;
   virtual   void     ReadChar(Char_t       &c) = 0;
   virtual   void     ReadUChar(UChar_t     &c) = 0;
   virtual   void     ReadShort(Short_t     &s) = 0;
   virtual   void     ReadUShort(UShort_t   &s) = 0;
   virtual   void     ReadInt(Int_t         &i) = 0;
   virtual   void     ReadUInt(UInt_t       &i) = 0;
   virtual   void     ReadLong(Long_t       &l) = 0;
   virtual   void     ReadULong(ULong_t     &l) = 0;
   virtual   void     ReadLong64(Long64_t   &l) = 0;
   virtual   void     ReadULong64(ULong64_t &l) = 0;
   virtual   void     ReadFloat(Float_t     &f) = 0;
   virtual   void     ReadDouble(Double_t   &d) = 0;
   virtual   void     ReadCharP(Char_t      *c) = 0;
   virtual   void     ReadTString(TString   &s) = 0;
   virtual   void     ReadStdString(std::string *s) = 0;
   virtual   void     ReadCharStar(char* &s) = 0;

   virtual inline void ReadStdString(std::string &s) { ReadStdString(&s); }

   virtual   void     WriteBool(Bool_t       b) = 0;
   virtual   void     WriteChar(Char_t       c) = 0;
   virtual   void     WriteUChar(UChar_t     c) = 0;
   virtual   void     WriteShort(Short_t     s) = 0;
   virtual   void     WriteUShort(UShort_t   s) = 0;
   virtual   void     WriteInt(Int_t         i) = 0;
   virtual   void     WriteUInt(UInt_t       i) = 0;
   virtual   void     WriteLong(Long_t       l) = 0;
   virtual   void     WriteULong(ULong_t     l) = 0;
   virtual   void     WriteLong64(Long64_t   l) = 0;
   virtual   void     WriteULong64(ULong64_t l) = 0;
   virtual   void     WriteFloat(Float_t     f) = 0;
   virtual   void     WriteDouble(Double_t   d) = 0;
   virtual   void     WriteCharP(const Char_t *c) = 0;
   virtual   void     WriteTString(const TString &s) = 0;
   virtual   void     WriteStdString(const std::string *s) = 0;
   virtual   void     WriteCharStar(char *s) = 0;

   virtual inline void WriteStdString(std::string &s) { WriteStdString(&s); }

   // Special basic ROOT objects and collections
   virtual   TProcessID *GetLastProcessID(TRefTable *reftable) const = 0;
   virtual   UInt_t      GetTRefExecId() = 0;
   virtual   TProcessID *ReadProcessID(UShort_t pidf) = 0;
   virtual   UShort_t    WriteProcessID(TProcessID *pid) = 0;

   // Utilities for TStreamerInfo
   virtual   void     ForceWriteInfo(TVirtualStreamerInfo *info, Bool_t force) = 0;
   virtual   void     ForceWriteInfoClones(TClonesArray *a) = 0;
   virtual   Int_t    ReadClones (TClonesArray *a, Int_t nobjects, Version_t objvers) = 0;
   virtual   Int_t    WriteClones(TClonesArray *a, Int_t nobjects) = 0;

   // Utilities for TClass
   virtual   Int_t    ReadClassEmulated(const TClass *cl, void *object, const TClass *onfile_class = 0) = 0;
   virtual   Int_t    ReadClassBuffer(const TClass *cl, void *pointer, const TClass *onfile_class = 0) = 0;
   virtual   Int_t    ReadClassBuffer(const TClass *cl, void *pointer, Int_t version, UInt_t start, UInt_t count, const TClass *onfile_class = 0) = 0;
   virtual   Int_t    WriteClassBuffer(const TClass *cl, void *pointer) = 0;

   // Utilites to streamer using sequences.
   virtual Int_t ApplySequence(const TStreamerInfoActions::TActionSequence &sequence, void *object) = 0;
   virtual Int_t ApplySequenceVecPtr(const TStreamerInfoActions::TActionSequence &sequence, void *start_collection, void *end_collection) = 0;
   virtual Int_t ApplySequence(const TStreamerInfoActions::TActionSequence &sequence, void *start_collection, void *end_collection) = 0;

   static TClass *GetClass(const std::type_info &typeinfo);
   static TClass *GetClass(const char *className);

   ClassDef(TBuffer,0)  //Buffer base class used for serializing objects
};

//---------------------- TBuffer default external operators --------------------

inline TBuffer &operator>>(TBuffer &buf, Bool_t &b)   { buf.ReadBool(b);   return buf; }
inline TBuffer &operator>>(TBuffer &buf, Char_t &c)   { buf.ReadChar(c);   return buf; }
inline TBuffer &operator>>(TBuffer &buf, UChar_t &c)  { buf.ReadUChar(c);  return buf; }
inline TBuffer &operator>>(TBuffer &buf, Short_t &s)  { buf.ReadShort(s);  return buf; }
inline TBuffer &operator>>(TBuffer &buf, UShort_t &s) { buf.ReadUShort(s); return buf; }
inline TBuffer &operator>>(TBuffer &buf, Int_t &i)    { buf.ReadInt(i);    return buf; }
inline TBuffer &operator>>(TBuffer &buf, UInt_t &i)   { buf.ReadUInt(i);   return buf; }
inline TBuffer &operator>>(TBuffer &buf, Long_t &l)   { buf.ReadLong(l);   return buf; }
inline TBuffer &operator>>(TBuffer &buf, ULong_t &l)  { buf.ReadULong(l);  return buf; }
inline TBuffer &operator>>(TBuffer &buf, Long64_t &l) { buf.ReadLong64(l); return buf; }
inline TBuffer &operator>>(TBuffer &buf, ULong64_t &l){ buf.ReadULong64(l);return buf; }
inline TBuffer &operator>>(TBuffer &buf, Float_t &f)  { buf.ReadFloat(f);  return buf; }
inline TBuffer &operator>>(TBuffer &buf, Double_t &d) { buf.ReadDouble(d); return buf; }
inline TBuffer &operator>>(TBuffer &buf, Char_t *c)   { buf.ReadCharP(c);  return buf; }
inline TBuffer &operator>>(TBuffer &buf, TString &s)  { buf.ReadTString(s);return buf; }

inline TBuffer &operator<<(TBuffer &buf, Bool_t b)   { buf.WriteBool(b);   return buf; }
inline TBuffer &operator<<(TBuffer &buf, Char_t c)   { buf.WriteChar(c);   return buf; }
inline TBuffer &operator<<(TBuffer &buf, UChar_t c)  { buf.WriteUChar(c);  return buf; }
inline TBuffer &operator<<(TBuffer &buf, Short_t s)  { buf.WriteShort(s);  return buf; }
inline TBuffer &operator<<(TBuffer &buf, UShort_t s) { buf.WriteUShort(s); return buf; }
inline TBuffer &operator<<(TBuffer &buf, Int_t i)    { buf.WriteInt(i);    return buf; }
inline TBuffer &operator<<(TBuffer &buf, UInt_t i)   { buf.WriteUInt(i);   return buf; }
inline TBuffer &operator<<(TBuffer &buf, Long_t l)   { buf.WriteLong(l);   return buf; }
inline TBuffer &operator<<(TBuffer &buf, ULong_t l)  { buf.WriteULong(l);  return buf; }
inline TBuffer &operator<<(TBuffer &buf, Long64_t l) { buf.WriteLong64(l); return buf; }
inline TBuffer &operator<<(TBuffer &buf, ULong64_t l){ buf.WriteULong64(l);return buf; }
inline TBuffer &operator<<(TBuffer &buf, Float_t f)  { buf.WriteFloat(f);  return buf; }
inline TBuffer &operator<<(TBuffer &buf, Double_t d) { buf.WriteDouble(d); return buf; }
inline TBuffer &operator<<(TBuffer &buf, const Char_t *c)  { buf.WriteCharP(c);  return buf; }
inline TBuffer &operator<<(TBuffer &buf, const TString &s) { buf.WriteTString(s);return buf; }

#ifndef __CINT__

#if defined(R__SOLARIS) && defined(R__GNU)
#include <typeinfo>
#endif

template <class Tmpl> TBuffer &operator>>(TBuffer &buf, Tmpl *&obj)
{
   // Read TObject derived classes from a TBuffer. Need to provide
   // custom version for non-TObject derived classes.

   // This operator has to be a templated and/or automatically
   // generated if we want to be able to check the type of the
   // incoming object. I.e. a operator>>(TBuffer &buf, TObject *&)
   // would not be sufficient to pass the information 'which class do we want'
   // since the pointer could be zero (so typeid(*obj) is not usable).

   auto cl = TClass::GetClass<Tmpl>();
   obj = (Tmpl *) ( (void*) buf.ReadObjectAny(cl) );
   return buf;
}

template <class Tmpl> TBuffer &operator<<(TBuffer &buf, const Tmpl *obj)
{
   auto cl = (obj) ? TClass::GetClass<Tmpl>() : nullptr;
   buf.WriteObjectAny(obj, cl);
   return buf;
}
#else
template <class Tmpl> TBuffer &operator>>(TBuffer &buf, Tmpl *&obj);
template <class Tmpl> TBuffer &operator<<(TBuffer &buf, Tmpl *&obj);
#endif

template <class T>
inline Int_t TBuffer::WriteObject(const T *objptr, Bool_t cacheReuse)
{
   auto cl = (objptr) ? TClass::GetClass<T>() : nullptr;
   return WriteObjectAny(objptr, cl, cacheReuse);
}

#endif // ROOT_TBuffer
