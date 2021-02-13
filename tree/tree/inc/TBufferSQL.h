// @(#)root/tree:$Id$
// Author: Philippe Canal 2005

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TBufferSQL
#define ROOT_TBufferSQL

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TBufferSQL                                                           //
//                                                                      //
// Implement TBuffer for a SQL backend                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TBufferFile.h"

#include <vector>
#include <string>

class TSQLRow;

class TBufferSQL final : public TBufferFile {

private:
   std::vector<Int_t>::const_iterator fIter;

   std::vector<Int_t>  *fColumnVec{nullptr};   //!
   TString             *fInsertQuery{nullptr}; //!
   TSQLRow            **fRowPtr{nullptr};      //!

   // TBuffer objects cannot be copied or assigned
   TBufferSQL(const TBufferSQL &);        // not implemented
   void operator=(const TBufferSQL &);    // not implemented

public:
   TBufferSQL();
   TBufferSQL(TBuffer::EMode mode, std::vector<Int_t> *vc, TString *insert_query, TSQLRow **rowPtr);
   TBufferSQL(TBuffer::EMode mode, Int_t bufsiz, std::vector<Int_t> *vc, TString *insert_query, TSQLRow **rowPtr);
   TBufferSQL(TBuffer::EMode mode, Int_t bufsiz, std::vector<Int_t> *vc, TString *insert_query, TSQLRow **rowPtr,void *buf, Bool_t adopt = kTRUE);
   ~TBufferSQL();

   void ResetOffset();

   void     ReadBool(Bool_t       &b) final;
   void     ReadChar(Char_t       &c) final;
   void     ReadUChar(UChar_t     &c) final;
   void     ReadShort(Short_t     &s) final;
   void     ReadUShort(UShort_t   &s) final;
   void     ReadInt(Int_t         &i) final;
   void     ReadUInt(UInt_t       &i) final;
   void     ReadLong(Long_t       &l) final;
   void     ReadULong(ULong_t     &l) final;
   void     ReadLong64(Long64_t   &l) final;
   void     ReadULong64(ULong64_t &l) final;
   void     ReadFloat(Float_t     &f) final;
   void     ReadDouble(Double_t   &d) final;
   void     ReadCharP(Char_t      *c) final;
   void     ReadTString(TString   &s) final;
   void     ReadStdString(std::string *s) final;
   using    TBuffer::ReadStdString;
   void     ReadCharStar(char* &s) final;

   void     WriteBool(Bool_t       b) final;
   void     WriteChar(Char_t       c) final;
   void     WriteUChar(UChar_t     c) final;
   void     WriteShort(Short_t     s) final;
   void     WriteUShort(UShort_t   s) final;
   void     WriteInt(Int_t         i) final;
   void     WriteUInt(UInt_t       i) final;
   void     WriteLong(Long_t       l) final;
   void     WriteULong(ULong_t     l) final;
   void     WriteLong64(Long64_t   l) final;
   void     WriteULong64(ULong64_t l) final;
   void     WriteFloat(Float_t     f) final;
   void     WriteDouble(Double_t   d) final;
   void     WriteCharP(const Char_t *c) final;
   void     WriteTString(const TString  &s) final;
   void     WriteStdString(const std::string *s) final;
   using    TBuffer::WriteStdString;
   void     WriteCharStar(char *s) final;

   void     WriteFastArray(const Bool_t    *b, Int_t n) final;
   void     WriteFastArray(const Char_t    *c, Int_t n) final;
   void     WriteFastArrayString(const Char_t   *c, Int_t n) final;
   void     WriteFastArray(const UChar_t   *c, Int_t n) final;
   void     WriteFastArray(const Short_t   *h, Int_t n) final;
   void     WriteFastArray(const UShort_t  *h, Int_t n) final;
   void     WriteFastArray(const Int_t     *i, Int_t n) final;
   void     WriteFastArray(const UInt_t    *i, Int_t n) final;
   void     WriteFastArray(const Long_t    *l, Int_t n) final;
   void     WriteFastArray(const ULong_t   *l, Int_t n) final;
   void     WriteFastArray(const Long64_t  *l, Int_t n) final;
   void     WriteFastArray(const ULong64_t *l, Int_t n) final;
   void     WriteFastArray(const Float_t   *f, Int_t n) final;
   void     WriteFastArray(const Double_t  *d, Int_t n) final;
   void     WriteFastArray(void  *start,  const TClass *cl, Int_t n=1, TMemberStreamer *s=nullptr) final;
   Int_t    WriteFastArray(void **startp, const TClass *cl, Int_t n=1, Bool_t isPreAlloc=kFALSE, TMemberStreamer *s=nullptr) final;

   void     ReadFastArray(Bool_t    *, Int_t ) final;
   void     ReadFastArray(Char_t    *, Int_t ) final;
   void     ReadFastArrayString(Char_t   *, Int_t ) final;
   void     ReadFastArray(UChar_t   *, Int_t ) final;
   void     ReadFastArray(Short_t   *, Int_t ) final;
   void     ReadFastArray(UShort_t  *, Int_t ) final;
   void     ReadFastArray(Int_t     *, Int_t ) final;
   void     ReadFastArray(UInt_t    *, Int_t ) final;
   void     ReadFastArray(Long_t    *, Int_t ) final;
   void     ReadFastArray(ULong_t   *, Int_t ) final;
   void     ReadFastArray(Long64_t  *, Int_t ) final;
   void     ReadFastArray(ULong64_t *, Int_t ) final;
   void     ReadFastArray(Float_t   *, Int_t ) final;
   void     ReadFastArray(Double_t  *, Int_t ) final;
   void     ReadFastArrayFloat16(Float_t  *f, Int_t n, TStreamerElement *ele=nullptr) final;
   void     ReadFastArrayDouble32(Double_t  *d, Int_t n, TStreamerElement *ele=nullptr) final;
   void     ReadFastArrayWithFactor(Float_t *ptr, Int_t n, Double_t factor, Double_t minvalue)  final;
   void     ReadFastArrayWithNbits(Float_t *ptr, Int_t n, Int_t nbits) final;
   void     ReadFastArrayWithFactor(Double_t *ptr, Int_t n, Double_t factor, Double_t minvalue) final;
   void     ReadFastArrayWithNbits(Double_t *ptr, Int_t n, Int_t nbits)  final;
   void     ReadFastArray(void  *, const TClass *, Int_t n=1, TMemberStreamer *s=nullptr, const TClass *onFileClass=nullptr) final;
   void     ReadFastArray(void **, const TClass *, Int_t n=1, Bool_t isPreAlloc=kFALSE, TMemberStreamer *s=nullptr, const TClass *onFileClass=nullptr) final;

   ClassDefOverride(TBufferSQL, 0); // Implementation of TBuffer to load and write to a SQL database

};

#endif


