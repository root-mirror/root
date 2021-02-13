// @(#)root/sql
// Author: Sergey Linev  20/11/2005

#ifndef ROOT_TBufferSQL2
#define ROOT_TBufferSQL2

#include "TBufferText.h"
#include "TString.h"
#include <string>

class TMap;
class TExMap;
class TVirtualStreamerInfo;
class TStreamerElement;
class TObjArray;
class TMemberStreamer;

class TSQLStackObj;
class TSQLServer;
class TSQLResult;
class TSQLRow;
class TSQLFile;
class TSQLStructure;
class TSQLObjectData;
class TSQLClassInfo;

class TBufferSQL2 final : public TBufferText {

   friend class TSQLStructure;

protected:
   TSQLFile *fSQL;               ///<!   instance of TSQLFile
   Int_t fIOVersion;             ///<!   I/O version from TSQLFile
   TSQLStructure *fStructure;    ///<!   structures, created by object storing
   TSQLStructure *fStk;          ///<!   pointer on current active structure (stack head)
   TString fReadBuffer;          ///<!   Buffer for read value
   Int_t fErrorFlag;             ///<!   Error id value
   Int_t fCompressLevel;         ///<!   compress level used to minimize size of data in database
   Int_t fReadVersionBuffer;     ///<!   buffer, used to by ReadVersion method
   Long64_t fObjIdCounter;       ///<!   counter of objects id
   Bool_t fIgnoreVerification;   ///<!   ignore verification of names
   TSQLObjectData *fCurrentData; ///<!
   TObjArray *fObjectsInfos;     ///<!   array of objects info for selected key
   Long64_t fFirstObjId;         ///<!   id of first object to be read from the database
   Long64_t fLastObjId;          ///<!   id of last object correspond to this key
   TMap *fPoolsMap;              ///<!   map of pools with data from different tables

   // TBufferSQL2 objects cannot be copied or assigned
   TBufferSQL2(const TBufferSQL2 &);    // not implemented
   void operator=(const TBufferSQL2 &); // not implemented

   TBufferSQL2();

   // redefined protected virtual functions

   void WriteObjectClass(const void *actualObjStart, const TClass *actualClass, Bool_t cacheReuse) final;

   // end redefined protected virtual functions

   TSQLStructure *PushStack();
   TSQLStructure *PopStack();
   TSQLStructure *Stack(Int_t depth = 0);

   void WorkWithClass(const char *classname, Version_t classversion);
   void WorkWithElement(TStreamerElement *elem, Int_t comp_type);

   Int_t SqlReadArraySize();
   Bool_t SqlObjectInfo(Long64_t objid, TString &clname, Version_t &version);
   TSQLObjectData *SqlObjectData(Long64_t objid, TSQLClassInfo *sqlinfo);

   Bool_t SqlWriteBasic(Char_t value);
   Bool_t SqlWriteBasic(Short_t value);
   Bool_t SqlWriteBasic(Int_t value);
   Bool_t SqlWriteBasic(Long_t value);
   Bool_t SqlWriteBasic(Long64_t value);
   Bool_t SqlWriteBasic(Float_t value);
   Bool_t SqlWriteBasic(Double_t value);
   Bool_t SqlWriteBasic(Bool_t value);
   Bool_t SqlWriteBasic(UChar_t value);
   Bool_t SqlWriteBasic(UShort_t value);
   Bool_t SqlWriteBasic(UInt_t value);
   Bool_t SqlWriteBasic(ULong_t value);
   Bool_t SqlWriteBasic(ULong64_t value);
   Bool_t SqlWriteValue(const char *value, const char *tname);

   void SqlReadBasic(Char_t &value);
   void SqlReadBasic(Short_t &value);
   void SqlReadBasic(Int_t &value);
   void SqlReadBasic(Long_t &value);
   void SqlReadBasic(Long64_t &value);
   void SqlReadBasic(Float_t &value);
   void SqlReadBasic(Double_t &value);
   void SqlReadBasic(Bool_t &value);
   void SqlReadBasic(UChar_t &value);
   void SqlReadBasic(UShort_t &value);
   void SqlReadBasic(UInt_t &value);
   void SqlReadBasic(ULong_t &value);
   void SqlReadBasic(ULong64_t &value);
   const char *SqlReadValue(const char *tname);
   const char *SqlReadCharStarValue();

   Int_t SqlWriteObject(const void *obj, const TClass *objClass, Bool_t cacheReuse, TMemberStreamer *streamer = nullptr,
                        Int_t streamer_index = 0);
   void *SqlReadObject(void *obj, TClass **cl = nullptr, TMemberStreamer *streamer = nullptr, Int_t streamer_index = 0,
                       const TClass *onFileClass = nullptr);
   void *SqlReadObjectDirect(void *obj, TClass **cl, Long64_t objid, TMemberStreamer *streamer = nullptr,
                             Int_t streamer_index = 0, const TClass *onFileClass = nullptr);

   void StreamObjectExtra(void *obj, TMemberStreamer *streamer, const TClass *cl, Int_t n = 0,
                          const TClass *onFileClass = nullptr);

   template <typename T>
   R__ALWAYS_INLINE void SqlReadArrayContent(T *arr, Int_t arrsize, Bool_t withsize);

   template <typename T>
   R__ALWAYS_INLINE Int_t SqlReadArray(T *&arr, Bool_t is_static = kFALSE);

   template <typename T>
   R__ALWAYS_INLINE void SqlReadFastArray(T *arr, Int_t arrsize);

   template <typename T>
   R__ALWAYS_INLINE void SqlWriteArray(T *arr, Int_t arrsize, Bool_t withsize = kFALSE);

public:
   TBufferSQL2(TBuffer::EMode mode, TSQLFile *file = nullptr);
   virtual ~TBufferSQL2();

   void SetCompressionLevel(int level) { fCompressLevel = level; }

   TSQLStructure *GetStructure() const { return fStructure; }

   Int_t GetErrorFlag() const { return fErrorFlag; }

   void SetIgnoreVerification() { fIgnoreVerification = kTRUE; }

   TSQLStructure *SqlWriteAny(const void *obj, const TClass *cl, Long64_t objid);

   void *SqlReadAny(Long64_t keyid, Long64_t objid, TClass **cl, void *obj = nullptr);

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
   void ClassMember(const char *name, const char *typeName = 0, Int_t arrsize1 = -1, Int_t arrsize2 = -1) final;

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

   // end of redefined virtual functions

   ClassDefOverride(TBufferSQL2, 0); // a specialized TBuffer to convert data to SQL statements or read data from SQL tables
};

#endif
