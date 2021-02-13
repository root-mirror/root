// @(#)root/io:$Id: 697641b2b52ed3d97bb5bde0fb5d2ff4a2f6c24f $
// Author: Rene Brun   17/01/2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TBufferFile
#define ROOT_TBufferFile


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TBufferFile                                                          //
//                                                                      //
// The concrete implementation of TBuffer for writing/reading to/from a //
// ROOT file or socket.                                                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TBufferIO.h"
#include "Bytes.h"

#include <vector>
#include <string>

#ifdef R__OLDHPACC
namespace std {
   using ::string;
   using ::vector;
}
#endif

class TVirtualStreamerInfo;
class TStreamerInfo;
class TStreamerElement;
class TClass;
class TVirtualArray;
namespace TStreamerInfoActions {
   class TActionSequence;
}

class TBufferFile : public TBufferIO {

protected:
   typedef std::vector<TStreamerInfo*> InfoList_t;

   TStreamerInfo  *fInfo{nullptr};  ///< Pointer to TStreamerInfo object writing/reading the buffer
   InfoList_t      fInfoStack;     ///< Stack of pointers to the TStreamerInfos

   // Default ctor
   TBufferFile() {} // NOLINT: not allowed to use = default because of TObject::kIsOnHeap detection, see ROOT-10300

   // TBuffer objects cannot be copied or assigned
   TBufferFile(const TBufferFile &) = delete;       ///<  not implemented
   void operator=(const TBufferFile &) = delete;    ///<  not implemented

   Int_t  CheckByteCount(UInt_t startpos, UInt_t bcnt, const TClass *clss, const char* classname);
   void  CheckCount(UInt_t offset) override;
   UInt_t CheckObject(UInt_t offset, const TClass *cl, Bool_t readClass = kFALSE);

   void  WriteObjectClass(const void *actualObjStart, const TClass *actualClass, Bool_t cacheReuse) override;

public:
   enum { kStreamedMemberWise = BIT(14) }; //added to version number to know if a collection has been stored member-wise

   TBufferFile(TBuffer::EMode mode);
   TBufferFile(TBuffer::EMode mode, Int_t bufsiz);
   TBufferFile(TBuffer::EMode mode, Int_t bufsiz, void *buf, Bool_t adopt = kTRUE, ReAllocCharFun_t reallocfunc = nullptr);
   virtual ~TBufferFile();

   Int_t      CheckByteCount(UInt_t startpos, UInt_t bcnt, const TClass *clss) override;
   Int_t      CheckByteCount(UInt_t startpos, UInt_t bcnt, const char *classname) override;
   void       SetByteCount(UInt_t cntpos, Bool_t packInVersion = kFALSE) override;

   void       SkipVersion(const TClass *cl = nullptr) override;
   Version_t  ReadVersion(UInt_t *start = nullptr, UInt_t *bcnt = nullptr, const TClass *cl = nullptr) override;
   Version_t  ReadVersionNoCheckSum(UInt_t *start = nullptr, UInt_t *bcnt = nullptr) override;
   Version_t  ReadVersionForMemberWise(const TClass *cl = nullptr) override;
   UInt_t     WriteVersion(const TClass *cl, Bool_t useBcnt = kFALSE) override;
   UInt_t     WriteVersionMemberWise(const TClass *cl, Bool_t useBcnt = kFALSE) override;

   void      *ReadObjectAny(const TClass* cast) override;
   void       SkipObjectAny() override;

   void       IncrementLevel(TVirtualStreamerInfo* info) override;
   void       SetStreamerElementNumber(TStreamerElement*,Int_t) override {}
   void       DecrementLevel(TVirtualStreamerInfo*) override;
   TVirtualStreamerInfo  *GetInfo() override { return (TVirtualStreamerInfo*)fInfo; }
   void       ClassBegin(const TClass*, Version_t = -1) override {}
   void       ClassEnd(const TClass*) override {}
   void       ClassMember(const char*, const char* = 0, Int_t = -1, Int_t = -1) override {}

   Int_t      ReadBuf(void *buf, Int_t max) override;
   void       WriteBuf(const void *buf, Int_t max) override;

   char      *ReadString(char *s, Int_t max) override;
   void       WriteString(const char *s) override;

   TClass    *ReadClass(const TClass *cl = nullptr, UInt_t *objTag = nullptr) override;
   void       WriteClass(const TClass *cl) override;

   TObject   *ReadObject(const TClass *cl) override;

   using TBufferIO::CheckObject;

   // basic types and arrays of basic types
   void     ReadFloat16 (Float_t *f, TStreamerElement *ele = nullptr) override;
   void     WriteFloat16(Float_t *f, TStreamerElement *ele = nullptr) override;
   void     ReadDouble32 (Double_t *d, TStreamerElement *ele = nullptr) override;
   void     WriteDouble32(Double_t *d, TStreamerElement *ele = nullptr) override;
   void     ReadWithFactor(Float_t *ptr, Double_t factor, Double_t minvalue) override;
   void     ReadWithNbits(Float_t *ptr, Int_t nbits) override;
   void     ReadWithFactor(Double_t *ptr, Double_t factor, Double_t minvalue) override;
   void     ReadWithNbits(Double_t *ptr, Int_t nbits) override;

   Int_t    ReadArray(Bool_t    *&b) override;
   Int_t    ReadArray(Char_t    *&c) override;
   Int_t    ReadArray(UChar_t   *&c) override;
   Int_t    ReadArray(Short_t   *&h) override;
   Int_t    ReadArray(UShort_t  *&h) override;
   Int_t    ReadArray(Int_t     *&i) override;
   Int_t    ReadArray(UInt_t    *&i) override;
   Int_t    ReadArray(Long_t    *&l) override;
   Int_t    ReadArray(ULong_t   *&l) override;
   Int_t    ReadArray(Long64_t  *&l) override;
   Int_t    ReadArray(ULong64_t *&l) override;
   Int_t    ReadArray(Float_t   *&f) override;
   Int_t    ReadArray(Double_t  *&d) override;
   Int_t    ReadArrayFloat16(Float_t  *&f, TStreamerElement *ele = nullptr) override;
   Int_t    ReadArrayDouble32(Double_t  *&d, TStreamerElement *ele = nullptr) override;

   Int_t    ReadStaticArray(Bool_t    *b) override;
   Int_t    ReadStaticArray(Char_t    *c) override;
   Int_t    ReadStaticArray(UChar_t   *c) override;
   Int_t    ReadStaticArray(Short_t   *h) override;
   Int_t    ReadStaticArray(UShort_t  *h) override;
   Int_t    ReadStaticArray(Int_t     *i) override;
   Int_t    ReadStaticArray(UInt_t    *i) override;
   Int_t    ReadStaticArray(Long_t    *l) override;
   Int_t    ReadStaticArray(ULong_t   *l) override;
   Int_t    ReadStaticArray(Long64_t  *l) override;
   Int_t    ReadStaticArray(ULong64_t *l) override;
   Int_t    ReadStaticArray(Float_t   *f) override;
   Int_t    ReadStaticArray(Double_t  *d) override;
   Int_t    ReadStaticArrayFloat16(Float_t  *f, TStreamerElement *ele = nullptr) override;
   Int_t    ReadStaticArrayDouble32(Double_t  *d, TStreamerElement *ele = nullptr) override;

   void     ReadFastArray(Bool_t    *b, Int_t n) override;
   void     ReadFastArray(Char_t    *c, Int_t n) override;
   void     ReadFastArrayString(Char_t    *c, Int_t n) override;
   void     ReadFastArray(UChar_t   *c, Int_t n) override;
   void     ReadFastArray(Short_t   *h, Int_t n) override;
   void     ReadFastArray(UShort_t  *h, Int_t n) override;
   void     ReadFastArray(Int_t     *i, Int_t n) override;
   void     ReadFastArray(UInt_t    *i, Int_t n) override;
   void     ReadFastArray(Long_t    *l, Int_t n) override;
   void     ReadFastArray(ULong_t   *l, Int_t n) override;
   void     ReadFastArray(Long64_t  *l, Int_t n) override;
   void     ReadFastArray(ULong64_t *l, Int_t n) override;
   void     ReadFastArray(Float_t   *f, Int_t n) override;
   void     ReadFastArray(Double_t  *d, Int_t n) override;
   void     ReadFastArrayFloat16(Float_t  *f, Int_t n, TStreamerElement *ele = nullptr) override;
   void     ReadFastArrayDouble32(Double_t  *d, Int_t n, TStreamerElement *ele = nullptr) override;
   void     ReadFastArrayWithFactor(Float_t *ptr, Int_t n, Double_t factor, Double_t minvalue)  override;
   void     ReadFastArrayWithNbits(Float_t *ptr, Int_t n, Int_t nbits) override;
   void     ReadFastArrayWithFactor(Double_t *ptr, Int_t n, Double_t factor, Double_t minvalue) override;
   void     ReadFastArrayWithNbits(Double_t *ptr, Int_t n, Int_t nbits)  override;
   void     ReadFastArray(void  *start , const TClass *cl, Int_t n=1, TMemberStreamer *s = nullptr, const TClass* onFileClass = nullptr) override;
   void     ReadFastArray(void **startp, const TClass *cl, Int_t n=1, Bool_t isPreAlloc=kFALSE, TMemberStreamer *s = nullptr, const TClass* onFileClass = nullptr) override;

   void     WriteArray(const Bool_t    *b, Int_t n) override;
   void     WriteArray(const Char_t    *c, Int_t n) override;
   void     WriteArray(const UChar_t   *c, Int_t n) override;
   void     WriteArray(const Short_t   *h, Int_t n) override;
   void     WriteArray(const UShort_t  *h, Int_t n) override;
   void     WriteArray(const Int_t     *i, Int_t n) override;
   void     WriteArray(const UInt_t    *i, Int_t n) override;
   void     WriteArray(const Long_t    *l, Int_t n) override;
   void     WriteArray(const ULong_t   *l, Int_t n) override;
   void     WriteArray(const Long64_t  *l, Int_t n) override;
   void     WriteArray(const ULong64_t *l, Int_t n) override;
   void     WriteArray(const Float_t   *f, Int_t n) override;
   void     WriteArray(const Double_t  *d, Int_t n) override;
   void     WriteArrayFloat16(const Float_t  *f, Int_t n, TStreamerElement *ele = nullptr) override;
   void     WriteArrayDouble32(const Double_t  *d, Int_t n, TStreamerElement *ele = nullptr) override;

   void     WriteFastArray(const Bool_t    *b, Int_t n) override;
   void     WriteFastArray(const Char_t    *c, Int_t n) override;
   void     WriteFastArrayString(const Char_t    *c, Int_t n) override;
   void     WriteFastArray(const UChar_t   *c, Int_t n) override;
   void     WriteFastArray(const Short_t   *h, Int_t n) override;
   void     WriteFastArray(const UShort_t  *h, Int_t n) override;
   void     WriteFastArray(const Int_t     *i, Int_t n) override;
   void     WriteFastArray(const UInt_t    *i, Int_t n) override;
   void     WriteFastArray(const Long_t    *l, Int_t n) override;
   void     WriteFastArray(const ULong_t   *l, Int_t n) override;
   void     WriteFastArray(const Long64_t  *l, Int_t n) override;
   void     WriteFastArray(const ULong64_t *l, Int_t n) override;
   void     WriteFastArray(const Float_t   *f, Int_t n) override;
   void     WriteFastArray(const Double_t  *d, Int_t n) override;
   void     WriteFastArrayFloat16(const Float_t  *f, Int_t n, TStreamerElement *ele = nullptr) override;
   void     WriteFastArrayDouble32(const Double_t  *d, Int_t n, TStreamerElement *ele = nullptr) override;
   void     WriteFastArray(void  *start,  const TClass *cl, Int_t n=1, TMemberStreamer *s = nullptr) override;
   Int_t    WriteFastArray(void **startp, const TClass *cl, Int_t n=1, Bool_t isPreAlloc=kFALSE, TMemberStreamer *s = nullptr) override;

   void     StreamObject(void *obj, const std::type_info &typeinfo, const TClass* onFileClass = nullptr) override;
   void     StreamObject(void *obj, const char *className, const TClass* onFileClass = nullptr) override;
   void     StreamObject(void *obj, const TClass *cl, const TClass* onFileClass = nullptr) override;
   void     StreamObject(TObject *obj) override;

   void     ReadBool(Bool_t       &b) override;
   void     ReadChar(Char_t       &c) override;
   void     ReadUChar(UChar_t     &c) override;
   void     ReadShort(Short_t     &s) override;
   void     ReadUShort(UShort_t   &s) override;
   void     ReadInt(Int_t         &i) override;
   void     ReadUInt(UInt_t       &i) override;
   void     ReadLong(Long_t       &l) override;
   void     ReadULong(ULong_t     &l) override;
   void     ReadLong64(Long64_t   &l) override;
   void     ReadULong64(ULong64_t &l) override;
   void     ReadFloat(Float_t     &f) override;
   void     ReadDouble(Double_t   &d) override;
   void     ReadCharP(Char_t      *c) override;
   void     ReadTString(TString   &s) override;
   void     ReadStdString(std::string *s) override;
   using    TBuffer::ReadStdString;
   void     ReadCharStar(char* &s) override;

   void     WriteBool(Bool_t       b) override;
   void     WriteChar(Char_t       c) override;
   void     WriteUChar(UChar_t     c) override;
   void     WriteShort(Short_t     s) override;
   void     WriteUShort(UShort_t   s) override;
   void     WriteInt(Int_t         i) override;
   void     WriteUInt(UInt_t       i) override;
   void     WriteLong(Long_t       l) override;
   void     WriteULong(ULong_t     l) override;
   void     WriteLong64(Long64_t   l) override;
   void     WriteULong64(ULong64_t l) override;
   void     WriteFloat(Float_t     f) override;
   void     WriteDouble(Double_t   d) override;
   void     WriteCharP(const Char_t *c) override;
   void     WriteTString(const TString &s) override;
   using    TBuffer::WriteStdString;
   void     WriteStdString(const std::string *s) override;
   void     WriteCharStar(char *s) override;

   // Utilities for TClass
   Int_t  ReadClassEmulated(const TClass *cl, void *object, const TClass *onfile_class) override;
   Int_t  ReadClassBuffer(const TClass *cl, void *pointer, const TClass *onfile_class) override;
   Int_t  ReadClassBuffer(const TClass *cl, void *pointer, Int_t version, UInt_t start, UInt_t count, const TClass *onfile_class) override;
   Int_t  WriteClassBuffer(const TClass *cl, void *pointer) override;

   // Utilities to streamer using sequences.
   Int_t ApplySequence(const TStreamerInfoActions::TActionSequence &sequence, void *object) override;
   Int_t ApplySequenceVecPtr(const TStreamerInfoActions::TActionSequence &sequence, void *start_collection, void *end_collection) override;
   Int_t ApplySequence(const TStreamerInfoActions::TActionSequence &sequence, void *start_collection, void *end_collection) override;

   ClassDefOverride(TBufferFile,0)  //concrete implementation of TBuffer for writing/reading to/from a ROOT file or socket.
};


//---------------------- TBufferFile inlines ---------------------------------------

//______________________________________________________________________________
inline void TBufferFile::WriteBool(Bool_t b)
{
   if (fBufCur + sizeof(UChar_t) > fBufMax) AutoExpand(fBufSize+sizeof(UChar_t));
   tobuf(fBufCur, b);
}

//______________________________________________________________________________
inline void TBufferFile::WriteChar(Char_t c)
{
   if (fBufCur + sizeof(Char_t) > fBufMax) AutoExpand(fBufSize+sizeof(Char_t));
   tobuf(fBufCur, c);
}

//______________________________________________________________________________
inline void TBufferFile::WriteUChar(UChar_t c)
{
   if (fBufCur + sizeof(UChar_t) > fBufMax) AutoExpand(fBufSize+sizeof(UChar_t));
   tobuf(fBufCur, (Char_t)c);
}

//______________________________________________________________________________
inline void TBufferFile::WriteShort(Short_t h)
{
   if (fBufCur + sizeof(Short_t) > fBufMax) AutoExpand(fBufSize+sizeof(Short_t));
   tobuf(fBufCur, h);
}

//______________________________________________________________________________
inline void TBufferFile::WriteUShort(UShort_t h)
{
   if (fBufCur + sizeof(UShort_t) > fBufMax) AutoExpand(fBufSize+sizeof(UShort_t));
   tobuf(fBufCur, (Short_t)h);
}

//______________________________________________________________________________
inline void TBufferFile::WriteInt(Int_t i)
{
   if (fBufCur + sizeof(Int_t) > fBufMax) AutoExpand(fBufSize+sizeof(Int_t));
   tobuf(fBufCur, i);
}

//______________________________________________________________________________
inline void TBufferFile::WriteUInt(UInt_t i)
{
   if (fBufCur + sizeof(UInt_t) > fBufMax) AutoExpand(fBufSize+sizeof(UInt_t));
   tobuf(fBufCur, (Int_t)i);
}

//______________________________________________________________________________
inline void TBufferFile::WriteLong(Long_t l)
{
   if (fBufCur + sizeof(Long_t) > fBufMax) AutoExpand(fBufSize+sizeof(Long_t));
   tobuf(fBufCur, l);
}

//______________________________________________________________________________
inline void TBufferFile::WriteULong(ULong_t l)
{
   if (fBufCur + sizeof(ULong_t) > fBufMax) AutoExpand(fBufSize+sizeof(ULong_t));
   tobuf(fBufCur, (Long_t)l);
}

//______________________________________________________________________________
inline void TBufferFile::WriteLong64(Long64_t ll)
{
   if (fBufCur + sizeof(Long64_t) > fBufMax) AutoExpand(fBufSize+sizeof(Long64_t));
   tobuf(fBufCur, ll);
}

//______________________________________________________________________________
inline void TBufferFile::WriteULong64(ULong64_t ll)
{
   if (fBufCur + sizeof(ULong64_t) > fBufMax) AutoExpand(fBufSize+sizeof(ULong64_t));
   tobuf(fBufCur, (Long64_t)ll);
}

//______________________________________________________________________________
inline void TBufferFile::WriteFloat(Float_t f)
{
   if (fBufCur + sizeof(Float_t) > fBufMax) AutoExpand(fBufSize+sizeof(Float_t));
   tobuf(fBufCur, f);
}

//______________________________________________________________________________
inline void TBufferFile::WriteDouble(Double_t d)
{
   if (fBufCur + sizeof(Double_t) > fBufMax) AutoExpand(fBufSize+sizeof(Double_t));
   tobuf(fBufCur, d);
}

//______________________________________________________________________________
inline void TBufferFile::WriteCharP(const Char_t *c)
{
   WriteString(c);
}

//______________________________________________________________________________
inline void TBufferFile::ReadBool(Bool_t &b)
{
   frombuf(fBufCur, &b);
}

//______________________________________________________________________________
inline void TBufferFile::ReadChar(Char_t &c)
{
   frombuf(fBufCur, &c);
}

//______________________________________________________________________________
inline void TBufferFile::ReadUChar(UChar_t &c)
{
   TBufferFile::ReadChar((Char_t &)c);
}

//______________________________________________________________________________
inline void TBufferFile::ReadShort(Short_t &h)
{
   frombuf(fBufCur, &h);
}

//______________________________________________________________________________
inline void TBufferFile::ReadUShort(UShort_t &h)
{
   TBufferFile::ReadShort((Short_t &)h);
}

//______________________________________________________________________________
inline void TBufferFile::ReadInt(Int_t &i)
{
   frombuf(fBufCur, &i);
}

//______________________________________________________________________________
inline void TBufferFile::ReadUInt(UInt_t &i)
{
   TBufferFile::ReadInt((Int_t &)i);
}


// in implementation file because special case with old version
//______________________________________________________________________________
//inline void TBufferFile::ReadLong(Long_t &ll)
//{
//   frombuf(fBufCur, &ll);
//}

//______________________________________________________________________________
inline void TBufferFile::ReadULong(ULong_t &ll)
{
   TBufferFile::ReadLong((Long_t&)ll);
}


//______________________________________________________________________________
inline void TBufferFile::ReadLong64(Long64_t &ll)
{
   frombuf(fBufCur, &ll);
}

//______________________________________________________________________________
inline void TBufferFile::ReadULong64(ULong64_t &ll)
{
   TBufferFile::ReadLong64((Long64_t &)ll);
}

//______________________________________________________________________________
inline void TBufferFile::ReadFloat(Float_t &f)
{
   frombuf(fBufCur, &f);
}

//______________________________________________________________________________
inline void TBufferFile::ReadDouble(Double_t &d)
{
   frombuf(fBufCur, &d);
}

//______________________________________________________________________________
inline void TBufferFile::ReadCharP(Char_t *c)
{
   ReadString(c, -1);
}

//______________________________________________________________________________
inline Int_t TBufferFile::ReadArray(UChar_t *&c)
   {  return TBufferFile::ReadArray((Char_t *&)c); }
//______________________________________________________________________________
inline Int_t TBufferFile::ReadArray(UShort_t *&h)
   {  return TBufferFile::ReadArray((Short_t *&)h); }
//______________________________________________________________________________
inline Int_t TBufferFile::ReadArray(UInt_t *&i)
   {  return TBufferFile::ReadArray((Int_t *&)i); }
//______________________________________________________________________________
inline Int_t TBufferFile::ReadArray(ULong_t *&l)
   {  return TBufferFile::ReadArray((Long_t *&)l); }
//______________________________________________________________________________
inline Int_t TBufferFile::ReadArray(ULong64_t *&ll)
   {  return TBufferFile::ReadArray((Long64_t *&)ll); }

//______________________________________________________________________________
inline Int_t TBufferFile::ReadStaticArray(UChar_t *c)
   {  return TBufferFile::ReadStaticArray((Char_t *)c); }
//______________________________________________________________________________
inline Int_t TBufferFile::ReadStaticArray(UShort_t *h)
   {  return TBufferFile::ReadStaticArray((Short_t *)h); }
//______________________________________________________________________________
inline Int_t TBufferFile::ReadStaticArray(UInt_t *i)
   {  return TBufferFile::ReadStaticArray((Int_t *)i); }
//______________________________________________________________________________
inline Int_t TBufferFile::ReadStaticArray(ULong_t *l)
   {  return TBufferFile::ReadStaticArray((Long_t *)l); }
//______________________________________________________________________________
inline Int_t TBufferFile::ReadStaticArray(ULong64_t *ll)
   {  return TBufferFile::ReadStaticArray((Long64_t *)ll); }

//______________________________________________________________________________
inline void TBufferFile::ReadFastArray(UChar_t *c, Int_t n)
   {        TBufferFile::ReadFastArray((Char_t *)c, n); }
//______________________________________________________________________________
inline void TBufferFile::ReadFastArray(UShort_t *h, Int_t n)
   {        TBufferFile::ReadFastArray((Short_t *)h, n); }
//______________________________________________________________________________
inline void TBufferFile::ReadFastArray(UInt_t *i, Int_t n)
   {        TBufferFile::ReadFastArray((Int_t *)i, n); }
//______________________________________________________________________________
inline void TBufferFile::ReadFastArray(ULong_t *l, Int_t n)
   {        TBufferFile::ReadFastArray((Long_t *)l, n); }
//______________________________________________________________________________
inline void TBufferFile::ReadFastArray(ULong64_t *ll, Int_t n)
   {        TBufferFile::ReadFastArray((Long64_t *)ll, n); }

//______________________________________________________________________________
inline void TBufferFile::WriteArray(const UChar_t *c, Int_t n)
   {        TBufferFile::WriteArray((const Char_t *)c, n); }
//______________________________________________________________________________
inline void TBufferFile::WriteArray(const UShort_t *h, Int_t n)
   {        TBufferFile::WriteArray((const Short_t *)h, n); }
//______________________________________________________________________________
inline void TBufferFile::WriteArray(const UInt_t *i, Int_t n)
   {        TBufferFile::WriteArray((const Int_t *)i, n); }
//______________________________________________________________________________
inline void TBufferFile::WriteArray(const ULong64_t *ll, Int_t n)
   {        TBufferFile::WriteArray((const Long64_t *)ll, n); }

//______________________________________________________________________________
inline void TBufferFile::WriteFastArray(const UChar_t *c, Int_t n)
   {        TBufferFile::WriteFastArray((const Char_t *)c, n); }
//______________________________________________________________________________
inline void TBufferFile::WriteFastArray(const UShort_t *h, Int_t n)
   {        TBufferFile::WriteFastArray((const Short_t *)h, n); }
//______________________________________________________________________________
inline void TBufferFile::WriteFastArray(const UInt_t *i, Int_t n)
   {        TBufferFile::WriteFastArray((const Int_t *)i, n); }
//______________________________________________________________________________
inline void TBufferFile::WriteFastArray(const ULong64_t *ll, Int_t n)
   {        TBufferFile::WriteFastArray((const Long64_t *)ll, n); }

#endif
