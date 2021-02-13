// @(#)root/tree:$Id$
// Author: Philippe Canal and al. 08/2004

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TBufferSQL
\ingroup tree
Implement TBuffer for a SQL backend.
*/

#include "TError.h"

#include "TBasketSQL.h"
#include "TBufferSQL.h"
#include "TSQLResult.h"
#include "TSQLRow.h"

#include <cstdio>
#include <cstdlib>
#include <iostream>

ClassImp(TBufferSQL);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TBufferSQL::TBufferSQL(TBuffer::EMode mode, std::vector<Int_t> *vc,
                       TString *insert_query, TSQLRow ** r) :
   TBufferFile(mode),
   fColumnVec(vc), fInsertQuery(insert_query), fRowPtr(r)
{
   fIter = fColumnVec->begin();
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TBufferSQL::TBufferSQL(TBuffer::EMode mode, Int_t bufsiz, std::vector<Int_t> *vc,
                       TString *insert_query, TSQLRow ** r) :
   TBufferFile(mode,bufsiz),
   fColumnVec(vc), fInsertQuery(insert_query), fRowPtr(r)
{
   fIter = fColumnVec->begin();
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TBufferSQL::TBufferSQL(TBuffer::EMode mode, Int_t bufsiz, std::vector<Int_t> *vc,
                       TString *insert_query, TSQLRow ** r,
                       void *buf, Bool_t adopt) :
   TBufferFile(mode,bufsiz,buf,adopt),
   fColumnVec(vc), fInsertQuery(insert_query), fRowPtr(r)
{
   fIter = fColumnVec->begin();
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TBufferSQL::TBufferSQL() : TBufferFile(), fColumnVec(0),fInsertQuery(0),fRowPtr(0)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TBufferSQL::~TBufferSQL()
{
   delete fColumnVec;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator>>

void TBufferSQL::ReadBool(Bool_t &b)
{
   b = (Bool_t)atoi((*fRowPtr)->GetField(*fIter));

   if (fIter != fColumnVec->end()) ++fIter;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator>>

void TBufferSQL::ReadChar(Char_t &c)
{
   c = (Char_t)atoi((*fRowPtr)->GetField(*fIter));

   if (fIter != fColumnVec->end()) ++fIter;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator>>

void TBufferSQL::ReadShort(Short_t &h)
{
   h = (Short_t)atoi((*fRowPtr)->GetField(*fIter));

   if (fIter != fColumnVec->end()) ++fIter;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator>>

void TBufferSQL::ReadInt(Int_t &i)
{
   i = atoi((*fRowPtr)->GetField(*fIter));

   if (fIter != fColumnVec->end()) ++fIter;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator>>

void TBufferSQL::ReadFloat(Float_t &f)
{
   f = atof((*fRowPtr)->GetField(*fIter));

   if (fIter != fColumnVec->end()) ++fIter;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator>>

void TBufferSQL::ReadLong(Long_t &l)
{
   l = atol((*fRowPtr)->GetField(*fIter));

   if (fIter != fColumnVec->end()) ++fIter;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator>>

void TBufferSQL::ReadDouble(Double_t &d)
{
   d = atof((*fRowPtr)->GetField(*fIter));

   if (fIter != fColumnVec->end()) ++fIter;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator<<

void TBufferSQL::WriteBool(Bool_t    b)
{
   (*fInsertQuery) += b;
   (*fInsertQuery) += ",";
   if (fIter != fColumnVec->end()) ++fIter;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator<<

void TBufferSQL::WriteChar(Char_t    c)
{
   (*fInsertQuery) += c;
   (*fInsertQuery) += ",";
   if (fIter != fColumnVec->end()) ++fIter;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator<<

void TBufferSQL::WriteShort(Short_t   h)
{
   (*fInsertQuery) += h;
   (*fInsertQuery) += ",";
   if (fIter != fColumnVec->end()) ++fIter;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator<<

void TBufferSQL::WriteInt(Int_t     i)
{
   (*fInsertQuery) += i;
   (*fInsertQuery) += ",";
   if (fIter != fColumnVec->end()) ++fIter;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator<<

void TBufferSQL::WriteLong(Long_t    l)
{
   (*fInsertQuery) += l;
   (*fInsertQuery) += ",";
   if (fIter != fColumnVec->end()) ++fIter;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator<<

void TBufferSQL::WriteFloat(Float_t   f)
{
   (*fInsertQuery) += f;
   (*fInsertQuery) += ",";
   if (fIter != fColumnVec->end()) ++fIter;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator<<

void TBufferSQL::WriteDouble(Double_t  d)
{
   (*fInsertQuery) += d;
   (*fInsertQuery) += ",";
   if (fIter != fColumnVec->end()) ++fIter;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator>>

void TBufferSQL::ReadUChar(UChar_t& uc)
{
   uc = (UChar_t)atoi((*fRowPtr)->GetField(*fIter));

   if (fIter != fColumnVec->end()) ++fIter;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator>>

void TBufferSQL::ReadUShort(UShort_t& us)
{
   us = (UShort_t)atoi((*fRowPtr)->GetField(*fIter));

   if (fIter != fColumnVec->end()) ++fIter;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator>>

void TBufferSQL::ReadUInt(UInt_t& ui)
{
   TString val = (*fRowPtr)->GetField(*fIter);
   Int_t code = sscanf(val.Data(), "%u",&ui);
   if(code == 0) Error("operator>>(UInt_t&)","Error reading UInt_t");

   if (fIter != fColumnVec->end()) ++fIter;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator>>

void TBufferSQL::ReadULong(ULong_t& ul)
{
   TString val = (*fRowPtr)->GetField(*fIter);
   Int_t code = sscanf(val.Data(), "%lu",&ul);
   if(code == 0) Error("operator>>(ULong_t&)","Error reading ULong_t");

   if (fIter != fColumnVec->end()) ++fIter;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator>>

void TBufferSQL::ReadLong64(Long64_t &ll)
{
   TString val = (*fRowPtr)->GetField(*fIter);
   Int_t code = sscanf(val.Data(), "%lld",&ll);
   if(code == 0) Error("operator>>(ULong_t&)","Error reading Long64_t");

   if (fIter != fColumnVec->end()) ++fIter;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator>>

void TBufferSQL::ReadULong64(ULong64_t &ull)
{
   TString val = (*fRowPtr)->GetField(*fIter);
   Int_t code = sscanf(val.Data(), "%llu",&ull);
   if(code == 0) Error("operator>>(ULong_t&)","Error reading ULong64_t");

   if (fIter != fColumnVec->end()) ++fIter;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator>>

void TBufferSQL::ReadCharP(Char_t *str)
{
   strcpy(str,(*fRowPtr)->GetField(*fIter));  // Legacy interface, we have no way to know the user's buffer size ....
   if (fIter != fColumnVec->end()) ++fIter;
}

////////////////////////////////////////////////////////////////////////////////
/// Read a TString

void TBufferSQL::ReadTString(TString   &s)
{
   s = (*fRowPtr)->GetField(*fIter);
   if (fIter != fColumnVec->end()) ++fIter;
}

////////////////////////////////////////////////////////////////////////////////
/// Write a TString

void TBufferSQL::WriteTString(const TString   &s)
{
   (*fInsertQuery) += s;
   (*fInsertQuery) += ",";
   if (fIter != fColumnVec->end()) ++fIter;
}



////////////////////////////////////////////////////////////////////////////////
/// Read a std::string

void TBufferSQL::ReadStdString(std::string *s)
{
   TBufferFile::ReadStdString(s);
}

////////////////////////////////////////////////////////////////////////////////
/// Write a std::string

void TBufferSQL::WriteStdString(const std::string *s)
{
   TBufferFile::WriteStdString(s);
}

////////////////////////////////////////////////////////////////////////////////
/// Read a char* string

void TBufferSQL::ReadCharStar(char* &s)
{
   TBufferFile::ReadCharStar(s);
}

////////////////////////////////////////////////////////////////////////////////
/// Write a char* string

void TBufferSQL::WriteCharStar(char *s)
{
   TBufferFile::WriteCharStar(s);
}


// Method to send to database.

////////////////////////////////////////////////////////////////////////////////
/// Operator<<

void TBufferSQL::WriteUChar(UChar_t uc)
{
   (*fInsertQuery) += uc;
   (*fInsertQuery) += ",";
   ++fIter;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator<<

void TBufferSQL::WriteUShort(UShort_t us)
{
   (*fInsertQuery) += us;
   (*fInsertQuery) += ",";
   ++fIter;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator<<

void TBufferSQL::WriteUInt(UInt_t ui)
{
   (*fInsertQuery) += ui;
   (*fInsertQuery) += ",";
   ++fIter;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator<<

void TBufferSQL::WriteULong(ULong_t ul)
{
   (*fInsertQuery) += ul;
   (*fInsertQuery) += ",";
   ++fIter;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator<<

void TBufferSQL::WriteLong64(Long64_t ll)
{
   (*fInsertQuery) += ll;
   (*fInsertQuery) += ",";
   ++fIter;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator<<

void TBufferSQL::WriteULong64(ULong64_t ull)
{
   (*fInsertQuery) += ull;
   (*fInsertQuery) += ",";
   ++fIter;
}

////////////////////////////////////////////////////////////////////////////////
/// Operator<<

void TBufferSQL::WriteCharP(const Char_t *str)
{
   (*fInsertQuery) += "\"";
   (*fInsertQuery) += str;
   (*fInsertQuery) += "\",";
   ++fIter;
}

////////////////////////////////////////////////////////////////////////////////
/// WriteFastArray SQL implementation.

void TBufferSQL::WriteFastArray(const Bool_t *b, Int_t n)
{
   for(int i=0; i<n; ++i) {
      (*fInsertQuery) += b[i];
      (*fInsertQuery) += ",";
      ++fIter;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// WriteFastArray SQL implementation.

void TBufferSQL::WriteFastArray(const Char_t *c, Int_t n)
{
   for(int i=0; i<n; ++i) {
      (*fInsertQuery) += (Short_t)c[i];
      (*fInsertQuery) += ",";
      ++fIter;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// WriteFastArray SQL implementation.

void TBufferSQL::WriteFastArrayString(const Char_t *c, Int_t /* n */)
{
   (*fInsertQuery) += "\"";
   (*fInsertQuery) += c;
   (*fInsertQuery) += "\",";
   ++fIter;
}

////////////////////////////////////////////////////////////////////////////////
/// WriteFastArray SQL implementation.

void TBufferSQL::WriteFastArray(const UChar_t *uc, Int_t n)
{
   for(int i=0; i<n; ++i) {
      (*fInsertQuery) += uc[i];
      (*fInsertQuery) += ",";
      ++fIter;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// WriteFastArray SQL implementation.

void TBufferSQL::WriteFastArray(const Short_t *h, Int_t n)
{
   for(int i=0; i<n; ++i) {
      (*fInsertQuery) += h[i];
      (*fInsertQuery) += ",";
      ++fIter;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// WriteFastArray SQL implementation.

void TBufferSQL::WriteFastArray(const UShort_t *us, Int_t n)
{
   for(int i=0; i<n; ++i) {
      (*fInsertQuery) += us[i];
      (*fInsertQuery) += ",";
      ++fIter;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// WriteFastArray SQL implementation.

void TBufferSQL::WriteFastArray(const Int_t     *ii, Int_t n)
{
    //   std::cerr << "Column: " <<*fIter << "   i:" << *ii << std::endl;
   for(int i=0; i<n; ++i) {
      (*fInsertQuery) += ii[i];
      (*fInsertQuery) += ",";
      ++fIter;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// WriteFastArray SQL implementation.

void TBufferSQL::WriteFastArray(const UInt_t *ui, Int_t n)
{
   for(int i=0; i<n; ++i) {
      (*fInsertQuery) += ui[i];
      (*fInsertQuery) += ",";
      ++fIter;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// WriteFastArray SQL implementation.

void TBufferSQL::WriteFastArray(const Long_t    *l, Int_t n)
{
   for(int i=0; i<n; ++i) {
      (*fInsertQuery)+= l[i];
      (*fInsertQuery)+= ",";
      ++fIter;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// WriteFastArray SQL implementation.

void TBufferSQL::WriteFastArray(const ULong_t   *ul, Int_t n)
{
   for(int i=0; i<n; ++i) {
      (*fInsertQuery) += ul[i];
      (*fInsertQuery) += ",";
      ++fIter;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// WriteFastArray SQL implementation.

void TBufferSQL::WriteFastArray(const Long64_t  *l, Int_t n)
{
   for(int i=0; i<n; ++i) {
      (*fInsertQuery) += l[i];
      (*fInsertQuery) += ",";
      ++fIter;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// WriteFastArray SQL implementation.

void TBufferSQL::WriteFastArray(const ULong64_t *ul, Int_t n)
{
   for(int i=0; i<n; ++i) {
      (*fInsertQuery) += ul[i];
      (*fInsertQuery) += ",";
      ++fIter;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// WriteFastArray SQL implementation.

void TBufferSQL::WriteFastArray(const Float_t   *f, Int_t n)
{
   for(int i=0; i<n; ++i) {
      (*fInsertQuery) += f[i];
      (*fInsertQuery) += ",";
      ++fIter;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// WriteFastArray SQL implementation.

void TBufferSQL::WriteFastArray(const Double_t  *d, Int_t n)
{
   for(int i=0; i<n; ++i) {
      (*fInsertQuery) += d[i];
      (*fInsertQuery )+= ",";
      ++fIter;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// WriteFastArray SQL implementation.

void TBufferSQL::WriteFastArray(void*, const TClass*, Int_t, TMemberStreamer *)
{
   Fatal("WriteFastArray(void*, const TClass*, Int_t, TMemberStreamer *)","Not implemented yet");
}

////////////////////////////////////////////////////////////////////////////////
/// WriteFastArray SQL implementation.

Int_t TBufferSQL::WriteFastArray(void **, const TClass*, Int_t, Bool_t, TMemberStreamer*)
{
   Fatal("WriteFastArray(void **, const TClass*, Int_t, Bool_t, TMemberStreamer*)","Not implemented yet");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// ReadFastArray SQL implementation.

void TBufferSQL::ReadFastArray(Bool_t *b, Int_t n)
{
   for(int i=0; i<n; ++i) {
      b[i] = (Bool_t)atoi((*fRowPtr)->GetField(*fIter));
      ++fIter;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// ReadFastArray SQL implementation.

void TBufferSQL::ReadFastArray(Char_t *c, Int_t n)
{
   for(int i=0; i<n; ++i) {
      c[i] = (Char_t)atoi((*fRowPtr)->GetField(*fIter));
      ++fIter;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// ReadFastArray SQL implementation.

void TBufferSQL::ReadFastArrayString(Char_t *c, Int_t /* n */)
{
   strcpy(c,((*fRowPtr)->GetField(*fIter)));
   ++fIter;
}

////////////////////////////////////////////////////////////////////////////////
/// ReadFastArray SQL implementation.

void TBufferSQL::ReadFastArray(UChar_t *uc, Int_t n)
{
   for(int i=0; i<n; ++i) {
      uc[i] = (UChar_t)atoi((*fRowPtr)->GetField(*fIter));
      ++fIter;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// ReadFastArray SQL implementation.

void TBufferSQL::ReadFastArray(Short_t *s, Int_t n)
{
   for(int i=0; i<n; ++i) {
      s[i] = (Short_t)atoi((*fRowPtr)->GetField(*fIter));
      ++fIter;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// ReadFastArray SQL implementation.

void TBufferSQL::ReadFastArray(UShort_t *us, Int_t n)
{
   for(int i=0; i<n; ++i) {
      us[i] = (UShort_t)atoi((*fRowPtr)->GetField(*fIter));
      ++fIter;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// ReadFastArray SQL implementation.

void     TBufferSQL::ReadFastArray(Int_t *in, Int_t n)
{
   for(int i=0; i<n; ++i) {
      in[i] = atoi((*fRowPtr)->GetField(*fIter));
      ++fIter;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// ReadFastArray SQL implementation.

void     TBufferSQL::ReadFastArray(UInt_t *ui, Int_t n)
{
   for(int i=0; i<n; ++i) {
      ui[i] = atoi((*fRowPtr)->GetField(*fIter));
      ++fIter;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// ReadFastArray SQL implementation.

void TBufferSQL::ReadFastArray(Long_t *l, Int_t n)
{
   for(int i=0; i<n; ++i) {
      l[i] = atol((*fRowPtr)->GetField(*fIter));
      ++fIter;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// ReadFastArray SQL implementation.

void TBufferSQL::ReadFastArray(ULong_t   *ul, Int_t n)
{
   for(int i=0; i<n; ++i) {
      (*this) >> ul[i];
   }
}

////////////////////////////////////////////////////////////////////////////////
/// ReadFastArray SQL implementation.

void TBufferSQL::ReadFastArray(Long64_t  *ll, Int_t n)
{
   for(int i=0; i<n; ++i) {
      (*this) >> ll[i];
   }
}

////////////////////////////////////////////////////////////////////////////////
/// ReadFastArray SQL implementation.

void TBufferSQL::ReadFastArray(ULong64_t *ull, Int_t n)
{
   for(int i=0; i<n; ++i) {
      (*this) >> ull[i];
   }
}

////////////////////////////////////////////////////////////////////////////////
/// ReadFastArray SQL implementation.

void TBufferSQL::ReadFastArray(Float_t   *f, Int_t n)
{
   for(int i=0; i<n; ++i) {
      f[i] = atof((*fRowPtr)->GetField(*fIter));
      ++fIter;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// ReadFastArray SQL implementation.

void TBufferSQL::ReadFastArray(Double_t *d, Int_t n)
{
   for(int i=0; i<n; ++i) {
      d[i] = atof((*fRowPtr)->GetField(*fIter));
      ++fIter;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// ReadFastArray SQL implementation.

void     TBufferSQL::ReadFastArrayFloat16(Float_t  *, Int_t , TStreamerElement *)
{
   Fatal("ReadFastArrayFloat16(Float_t  *, Int_t , TStreamerElement *)","Not implemented yet");
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Float16_t from buffer

void TBufferSQL::ReadFastArrayWithFactor(Float_t  *, Int_t , Double_t /* factor */, Double_t /* minvalue */)
{
   Fatal("ReadFastArrayWithFactor(Float_t  *, Int_t, Double_t, Double_t)","Not implemented yet");
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Float16_t from buffer

void TBufferSQL::ReadFastArrayWithNbits(Float_t  *, Int_t , Int_t /*nbits*/)
{
   Fatal("ReadFastArrayWithNbits(Float_t  *, Int_t , Int_t )","Not implemented yet");
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Double32_t from buffer

void TBufferSQL::ReadFastArrayWithFactor(Double_t  *, Int_t , Double_t /* factor */, Double_t /* minvalue */)
{
   Fatal("ReadFastArrayWithFactor(Double_t  *, Int_t, Double_t, Double_t)","Not implemented yet");
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of Double32_t from buffer

void TBufferSQL::ReadFastArrayWithNbits(Double_t  *, Int_t , Int_t /*nbits*/)
{
   Fatal("ReadFastArrayWithNbits(Double_t  *, Int_t , Int_t )","Not implemented yet");
}

////////////////////////////////////////////////////////////////////////////////
/// ReadFastArray SQL implementation.

void     TBufferSQL::ReadFastArrayDouble32(Double_t  *, Int_t , TStreamerElement *)
{
   Fatal("ReadFastArrayDouble32(Double_t  *, Int_t , TStreamerElement *)","Not implemented yet");
}

////////////////////////////////////////////////////////////////////////////////
/// ReadFastArray SQL implementation.

void     TBufferSQL::ReadFastArray(void  *, const TClass *, Int_t, TMemberStreamer *, const TClass *)
{
   Fatal("ReadFastArray(void  *, const TClass *, Int_t, TMemberStreamer *, const TClass *)","Not implemented yet");
}

////////////////////////////////////////////////////////////////////////////////
/// ReadFastArray SQL implementation.

void     TBufferSQL::ReadFastArray(void **, const TClass *, Int_t, Bool_t, TMemberStreamer *, const TClass *)
{
   Fatal("ReadFastArray(void **, const TClass *, Int_t, Bool_t, TMemberStreamer *, const TClass *)","Not implemented yet");
}

////////////////////////////////////////////////////////////////////////////////
/// Reset Offset.

void TBufferSQL::ResetOffset()
{
   fIter = fColumnVec->begin();
}

#if 0
////////////////////////////////////////////////////////////////////////////////

void TBufferSQL::insert_test(const char* dsn, const char* usr,
                             const char* pwd, const TString& tblname)
{
   TString str;
   TString select = "select * from ";
   TString sql;
   TSQLStatement* stmt;
   sql = select + "ins";

   con = gSQLDriverManager->GetConnection(dsn,usr,pwd);

   if(!con)
      printf("\n\n\nConnection NOT Successful\n\n\n");
   else
      printf("\n\n\nConnection Sucessful\n\n\n");

   stmt = con->CreateStatement(0, odbc::ResultSet::CONCUR_READ_ONLY);

   ptr = stmt->ExecuteQuery(sql.Data());
   if(!ptr) printf("No recorSet found!");

   ptr->Next();
   ptr->MoveToInsertRow();
   std::cerr << "IsAfterLast(): " << ptr->IsAfterLast() << std::endl;
   ptr->UpdateInt(1, 5555);
   ptr->InsertRow();
   con->Commit();

   ptr1 = stmt->ExecuteQuery(sql.Data());

}
#endif
