// @(#)root/io:$Id$
// Author: Rene Brun 17/01/2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/**
\file TBufferFile.cxx
\class TBufferFile
\ingroup IO

The concrete implementation of TBuffer for writing/reading to/from a ROOT file or socket.
*/

#include <string.h>
#include <typeinfo>
#include <string>

#include "TFile.h"
#include "TBufferFile.h"
#include "TExMap.h"
#include "TClass.h"
#include "TStorage.h"
#include "TError.h"
#include "TStreamer.h"
#include "TStreamerInfo.h"
#include "TStreamerElement.h"
#include "TSchemaRuleSet.h"
#include "TStreamerInfoActions.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"

#if (defined(__linux) || defined(__APPLE__)) && defined(__i386__) && \
     defined(__GNUC__)
#define USE_BSWAPCPY
#endif

#ifdef USE_BSWAPCPY
#include "Bswapcpy.h"
#endif


const UInt_t kNewClassTag       = 0xFFFFFFFF;
const UInt_t kClassMask         = 0x80000000;  // OR the class index with this
const UInt_t kByteCountMask     = 0x40000000;  // OR the byte count with this
const UInt_t kMaxMapCount       = 0x3FFFFFFE;  // last valid fMapCount and byte count
const Version_t kByteCountVMask = 0x4000;      // OR the version byte count with this
const Version_t kMaxVersion     = 0x3FFF;      // highest possible version number
const Int_t  kMapOffset         = 2;   // first 2 map entries are taken by null obj and self obj


ClassImp(TBufferFile);

////////////////////////////////////////////////////////////////////////////////
/// Thread-safe check on StreamerInfos of a TClass

static inline bool Class_Has_StreamerInfo(const TClass* cl)
{
   // NOTE: we do not need a R__LOCKGUARD2 since we know the
   //   mutex is available since the TClass constructor will make
   //   it
   R__LOCKGUARD(gInterpreterMutex);
   return cl->GetStreamerInfos()->GetLast()>1;
}

////////////////////////////////////////////////////////////////////////////////
/// Create an I/O buffer object. Mode should be either TBuffer::kRead or
/// TBuffer::kWrite. By default the I/O buffer has a size of
/// TBuffer::kInitialSize (1024) bytes.

TBufferFile::TBufferFile(TBuffer::EMode mode)
            :TBufferIO(mode),
             fInfo(nullptr), fInfoStack()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Create an I/O buffer object. Mode should be either TBuffer::kRead or
/// TBuffer::kWrite.

TBufferFile::TBufferFile(TBuffer::EMode mode, Int_t bufsiz)
            :TBufferIO(mode,bufsiz),
             fInfo(nullptr), fInfoStack()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Create an I/O buffer object.
/// Mode should be either TBuffer::kRead or
/// TBuffer::kWrite. By default the I/O buffer has a size of
/// TBuffer::kInitialSize (1024) bytes. An external buffer can be passed
/// to TBuffer via the buf argument. By default this buffer will be adopted
/// unless adopt is false.
/// If the new buffer is <b>not</b> adopted and no memory allocation routine
/// is provided, a Fatal error will be issued if the Buffer attempts to
/// expand.

TBufferFile::TBufferFile(TBuffer::EMode mode, Int_t bufsiz, void *buf, Bool_t adopt, ReAllocCharFun_t reallocfunc) :
   TBufferIO(mode,bufsiz,buf,adopt,reallocfunc),
   fInfo(nullptr), fInfoStack()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Delete an I/O buffer object.

TBufferFile::~TBufferFile()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Increment level.

void TBufferFile::IncrementLevel(TVirtualStreamerInfo* info)
{
   fInfoStack.push_back(fInfo);
   fInfo = (TStreamerInfo*)info;
}

////////////////////////////////////////////////////////////////////////////////
/// Decrement level.

void TBufferFile::DecrementLevel(TVirtualStreamerInfo* /*info*/)
{
   fInfo = fInfoStack.back();
   fInfoStack.pop_back();
}

////////////////////////////////////////////////////////////////////////////////
/// Handle old file formats.
///
/// Files written with versions older than 3.00/06 had a non-portable
/// implementation of Long_t/ULong_t. These types should not have been
/// used at all. However, because some users had already written many
/// files with these types we provide this dirty patch for "backward
/// compatibility"

static void frombufOld(char *&buf, Long_t *x)
{
#ifdef R__BYTESWAP
#ifdef R__B64
   char *sw = (char *)x;
   sw[0] = buf[7];
   sw[1] = buf[6];
   sw[2] = buf[5];
   sw[3] = buf[4];
   sw[4] = buf[3];
   sw[5] = buf[2];
   sw[6] = buf[1];
   sw[7] = buf[0];
#else
   char *sw = (char *)x;
   sw[0] = buf[3];
   sw[1] = buf[2];
   sw[2] = buf[1];
   sw[3] = buf[0];
#endif
#else
   memcpy(x, buf, sizeof(Long_t));
#endif
   buf += sizeof(Long_t);
}

////////////////////////////////////////////////////////////////////////////////
/// Read Long from TBuffer.

void TBufferFile::ReadLong(Long_t &l)
{
   TFile *file = (TFile*)fParent;
   if (file && file->GetVersion() < 30006) {
      frombufOld(fBufCur, &l);
   } else {
      frombuf(fBufCur, &l);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Read TString from TBuffer.

void TBufferFile::ReadTString(TString &s)
{
   Int_t   nbig;
   UChar_t nwh;
   *this >> nwh;
   if (nwh == 0) {
      s.UnLink();
      s.Zero();
   } else {
      if (nwh == 255)
         *this >> nbig;
      else
         nbig = nwh;

      s.Clobber(nbig);
      char *data = s.GetPointer();
      data[nbig] = 0;
      s.SetSize(nbig);
      ReadFastArray(data, nbig);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Write TString to TBuffer.

void TBufferFile::WriteTString(const TString &s)
{
   Int_t nbig = s.Length();
   UChar_t nwh;
   if (nbig > 254) {
      nwh = 255;
      *this << nwh;
      *this << nbig;
   } else {
      nwh = UChar_t(nbig);
      *this << nwh;
   }
   const char *data = s.GetPointer();
   WriteFastArray(data, nbig);
}

////////////////////////////////////////////////////////////////////////////////
/// Read std::string from TBuffer.

void TBufferFile::ReadStdString(std::string *obj)
{
   if (obj == nullptr) {
      Error("TBufferFile::ReadStdString","The std::string address is nullptr but should not");
      return;
   }
   Int_t   nbig;
   UChar_t nwh;
   *this >> nwh;
   if (nwh == 0)  {
       obj->clear();
   } else {
      if( obj->size() ) {
         // Insure that the underlying data storage is not shared
         (*obj)[0] = '\0';
      }
      if (nwh == 255)  {
         *this >> nbig;
         obj->resize(nbig,'\0');
         ReadFastArray((char*)obj->data(),nbig);
      }
      else  {
         obj->resize(nwh,'\0');
         ReadFastArray((char*)obj->data(),nwh);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Write std::string to TBuffer.

void TBufferFile::WriteStdString(const std::string *obj)
{
   if (obj==0) {
      *this << (UChar_t)0;
      WriteFastArray("",0);
      return;
   }

   UChar_t nwh;
   Int_t nbig = obj->length();
   if (nbig > 254) {
      nwh = 255;
      *this << nwh;
      *this << nbig;
   } else {
      nwh = UChar_t(nbig);
      *this << nwh;
   }
   WriteFastArray(obj->data(),nbig);
}

////////////////////////////////////////////////////////////////////////////////
/// Read char* from TBuffer.

void TBufferFile::ReadCharStar(char* &s)
{
   delete [] s;
   s = nullptr;

   Int_t nch;
   *this >> nch;
   if (nch > 0) {
      s = new char[nch+1];
      ReadFastArray(s, nch);
      s[nch] = 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Write char* into TBuffer.

void TBufferFile::WriteCharStar(char *s)
{
   Int_t nch = 0;
   if (s) {
      nch = strlen(s);
      *this  << nch;
      WriteFastArray(s,nch);
   } else {
      *this << nch;
   }

}

////////////////////////////////////////////////////////////////////////////////
/// Set byte count at position cntpos in the buffer. Generate warning if
/// count larger than kMaxMapCount. The count is excluded its own size.

void TBufferFile::SetByteCount(UInt_t cntpos, Bool_t packInVersion)
{
   UInt_t cnt = UInt_t(fBufCur - fBuffer) - cntpos - sizeof(UInt_t);
   char  *buf = (char *)(fBuffer + cntpos);

   // if true, pack byte count in two consecutive shorts, so it can
   // be read by ReadVersion()
   if (packInVersion) {
      union {
         UInt_t    cnt;
         Version_t vers[2];
      } v;
      v.cnt = cnt;
#ifdef R__BYTESWAP
      tobuf(buf, Version_t(v.vers[1] | kByteCountVMask));
      tobuf(buf, v.vers[0]);
#else
      tobuf(buf, Version_t(v.vers[0] | kByteCountVMask));
      tobuf(buf, v.vers[1]);
#endif
   } else
      tobuf(buf, cnt | kByteCountMask);

   if (cnt >= kMaxMapCount) {
      Error("WriteByteCount", "bytecount too large (more than %d)", kMaxMapCount);
      // exception
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Check byte count with current buffer position. They should
/// match. If not print warning and position buffer in correct
/// place determined by the byte count. Startpos is position of
/// first byte where the byte count is written in buffer.
/// Returns 0 if everything is ok, otherwise the bytecount offset
/// (< 0 when read too little, >0 when read too much).

Int_t TBufferFile::CheckByteCount(UInt_t startpos, UInt_t bcnt, const TClass *clss, const char *classname)
{
   if (!bcnt) return 0;

   Int_t  offset = 0;

   Long_t endpos = Long_t(fBuffer) + startpos + bcnt + sizeof(UInt_t);

   if (Long_t(fBufCur) != endpos) {
      offset = Int_t(Long_t(fBufCur) - endpos);

      const char *name = clss ? clss->GetName() : classname ? classname : 0;

      if (name) {
         if (offset < 0) {
            Error("CheckByteCount", "object of class %s read too few bytes: %d instead of %d",
                  name,bcnt+offset,bcnt);
         }
         if (offset > 0) {
            Error("CheckByteCount", "object of class %s read too many bytes: %d instead of %d",
                  name,bcnt+offset,bcnt);
            if (fParent)
               Warning("CheckByteCount","%s::Streamer() not in sync with data on file %s, fix Streamer()",
                       name, fParent->GetName());
            else
               Warning("CheckByteCount","%s::Streamer() not in sync with data, fix Streamer()",
                       name);
         }
      }
      if ( ((char *)endpos) > fBufMax ) {
         offset = fBufMax-fBufCur;
         Error("CheckByteCount",
               "Byte count probably corrupted around buffer position %d:\n\t%d for a possible maximum of %d",
               startpos, bcnt, offset);
         fBufCur = fBufMax;

      } else {

         fBufCur = (char *) endpos;

      }
   }
   return offset;
}

////////////////////////////////////////////////////////////////////////////////
/// Check byte count with current buffer position. They should
/// match. If not print warning and position buffer in correct
/// place determined by the byte count. Startpos is position of
/// first byte where the byte count is written in buffer.
/// Returns 0 if everything is ok, otherwise the bytecount offset
/// (< 0 when read too little, >0 when read too much).

Int_t TBufferFile::CheckByteCount(UInt_t startpos, UInt_t bcnt, const TClass *clss)
{
   if (!bcnt) return 0;
   return CheckByteCount( startpos, bcnt, clss, nullptr);
}

////////////////////////////////////////////////////////////////////////////////
/// Check byte count with current buffer position. They should
/// match. If not print warning and position buffer in correct
/// place determined by the byte count. Startpos is position of
/// first byte where the byte count is written in buffer.
/// Returns 0 if everything is ok, otherwise the bytecount offset
/// (< 0 when read too little, >0 when read too much).

Int_t TBufferFile::CheckByteCount(UInt_t startpos, UInt_t bcnt, const char *classname)
{
   if (!bcnt) return 0;
   return CheckByteCount( startpos, bcnt, nullptr, classname);
}

////////////////////////////////////////////////////////////////////////////////
/// Read a Float16_t from the buffer,
/// see comments about Float16_t encoding at TBufferFile::WriteFloat16().

void TBufferFile::ReadFloat16(Float_t *f, TStreamerElement *ele)
{
   if (ele && ele->GetFactor() != 0) {
      ReadWithFactor(f, ele->GetFactor(), ele->GetXmin());
   } else {
      Int_t nbits = 0;
      if (ele) nbits = (Int_t)ele->GetXmin();
      if (!nbits) nbits = 12;
      ReadWithNbits(f, nbits);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Read a Double32_t from the buffer,
/// see comments about Double32_t encoding at TBufferFile::WriteDouble32().

void TBufferFile::ReadDouble32(Double_t *d, TStreamerElement *ele)
{
   if (ele && ele->GetFactor() != 0) {
      ReadWithFactor(d, ele->GetFactor(), ele->GetXmin());
   } else {
      Int_t nbits = 0;
      if (ele) nbits = (Int_t)ele->GetXmin();
      if (!nbits) {
         //we read a float and convert it to double
         Float_t afloat;
         *this >> afloat;
         d[0] = (Double_t)afloat;
      } else {
         ReadWithNbits(d, nbits);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Read a Float16_t from the buffer when the factor and minimum value have been specified
/// see comments about Double32_t encoding at TBufferFile::WriteDouble32().

void TBufferFile::ReadWithFactor(Float_t *ptr, Double_t factor, Double_t minvalue)
{
   //a range was specified. We read an integer and convert it back to a double.
   UInt_t aint;
   frombuf(this->fBufCur,&aint);
   ptr[0] = (Float_t)(aint/factor + minvalue);
}

////////////////////////////////////////////////////////////////////////////////
/// Read a Float16_t from the buffer when the number of bits is specified (explicitly or not)
/// see comments about Float16_t encoding at TBufferFile::WriteFloat16().

void TBufferFile::ReadWithNbits(Float_t *ptr, Int_t nbits)
{
   //we read the exponent and the truncated mantissa of the float
   //and rebuild the float.
   union {
      Float_t fFloatValue;
      Int_t   fIntValue;
   } temp;
   UChar_t  theExp;
   UShort_t theMan;
   frombuf(this->fBufCur,&theExp);
   frombuf(this->fBufCur,&theMan);
   temp.fIntValue = theExp;
   temp.fIntValue <<= 23;
   temp.fIntValue |= (theMan & ((1<<(nbits+1))-1)) <<(23-nbits);
   if(1<<(nbits+1) & theMan) temp.fFloatValue = -temp.fFloatValue;
   ptr[0] = temp.fFloatValue;
}

////////////////////////////////////////////////////////////////////////////////
/// Read a Double32_t from the buffer when the factor and minimum value have been specified
/// see comments about Double32_t encoding at TBufferFile::WriteDouble32().

void TBufferFile::ReadWithFactor(Double_t *ptr, Double_t factor, Double_t minvalue)
{
   //a range was specified. We read an integer and convert it back to a double.
   UInt_t aint;
   frombuf(this->fBufCur,&aint);
   ptr[0] = (Double_t)(aint/factor + minvalue);
}

////////////////////////////////////////////////////////////////////////////////
/// Read a Double32_t from the buffer when the number of bits is specified (explicitly or not)
/// see comments about Double32_t encoding at TBufferFile::WriteDouble32().

void TBufferFile::ReadWithNbits(Double_t *ptr, Int_t nbits)
{
   //we read the exponent and the truncated mantissa of the float
   //and rebuild the float.
   union {
      Float_t fFloatValue;
      Int_t   fIntValue;
   } temp;
   UChar_t  theExp;
   UShort_t theMan;
   frombuf(this->fBufCur,&theExp);
   frombuf(this->fBufCur,&theMan);
   temp.fIntValue = theExp;
   temp.fIntValue <<= 23;
   temp.fIntValue |= (theMan & ((1<<(nbits+1))-1)) <<(23-nbits);
   if(1<<(nbits+1) & theMan) temp.fFloatValue = -temp.fFloatValue;
   ptr[0] = (Double_t)temp.fFloatValue;
}

////////////////////////////////////////////////////////////////////////////////
/// Write a Float16_t to the buffer.
///
/// The following cases are supported for streaming a Float16_t type
/// depending on the range declaration in the comment field of the data member:
/// Case | Example |
/// -----|---------|
///  A   | Float16_t     fNormal; |
///  B   | Float16_t     fTemperature; //[0,100]|
///  C   | Float16_t     fCharge;      //[-1,1,2]|
///  D   | Float16_t     fVertex[3];   //[-30,30,10]|
///  E   | Float16_t     fChi2;        //[0,0,6]|
///  F   | Int_t         fNsp;<br>Float16_t*    fPointValue;   //[fNsp][0,3]|
///
///   - In case A fNormal is converted from a Float_t to a Float_t with mantissa truncated to 12 bits
///   - In case B fTemperature is converted to a 32 bit unsigned integer
///   - In case C fCharge is converted to a 2 bits unsigned integer
///   - In case D the array elements of fVertex are converted to an unsigned 10 bits integer
///   - In case E fChi2 is converted to a Float_t with truncated precision at 6 bits
///   - In case F the fNsp elements of array fPointvalue are converted to an unsigned 32 bit integer
/// Note that the range specifier must follow the dimension specifier.
/// Case B has more precision (9 to 10 significative digits than case A (6 to 7 digits).
///
/// The range specifier has the general format: [xmin,xmax] or [xmin,xmax,nbits]
///   - [0,1];
///   - [-10,100];
///   - [-pi,pi], [-pi/2,pi/4],[-2pi,2*pi]
///   - [-10,100,16]
///   - [0,0,8]
/// if nbits is not specified, or nbits <2 or nbits>16 it is set to 16. If
/// (xmin==0 and xmax==0 and nbits <=14) the float word will have
/// its mantissa truncated to nbits significative bits.
///
/// ## IMPORTANT NOTE
/// ### NOTE 1
/// Lets assume an original variable float x:
/// When using the format [0,0,8] (ie range not specified) you get the best
/// relative precision when storing and reading back the truncated x, say xt.
/// The variance of (x-xt)/x will be better than when specifying a range
/// for the same number of bits. However the precision relative to the
/// range (x-xt)/(xmax-xmin) will be worst, and vice-versa.
/// The format [0,0,8] is also interesting when the range of x is infinite
/// or unknown.
///
/// ### NOTE 2
/// It is important to understand the difference with the meaning of nbits
///   - in case of [-1,1,nbits], nbits is the total number of bits used to make
/// the conversion from a float to an integer
///   - in case of [0,0,nbits], nbits is the number of bits used for the mantissa
///
///  See example of use of the Float16_t data type in tutorial double32.C
///  \image html tbufferfile_double32.gif

void TBufferFile::WriteFloat16(Float_t *f, TStreamerElement *ele)
{

   if (ele && ele->GetFactor() != 0) {
      //A range is specified. We normalize the double to the range and
      //convert it to an integer using a scaling factor that is a function of nbits.
      //see TStreamerElement::GetRange.
      Double_t x = f[0];
      Double_t xmin = ele->GetXmin();
      Double_t xmax = ele->GetXmax();
      if (x < xmin) x = xmin;
      if (x > xmax) x = xmax;
      UInt_t aint = UInt_t(0.5+ele->GetFactor()*(x-xmin)); *this << aint;
   } else {
      Int_t nbits = 0;
      //number of bits stored in fXmin (see TStreamerElement::GetRange)
      if (ele) nbits = (Int_t)ele->GetXmin();
      if (!nbits) nbits = 12;
      //a range is not specified, but nbits is.
      //In this case we truncate the mantissa to nbits and we stream
      //the exponent as a UChar_t and the mantissa as a UShort_t.
      union {
         Float_t fFloatValue;
         Int_t   fIntValue;
      };
      fFloatValue = f[0];
      UChar_t  theExp = (UChar_t)(0x000000ff & ((fIntValue<<1)>>24));
      UShort_t theMan = ((1<<(nbits+1))-1) & (fIntValue>>(23-nbits-1));
      theMan++;
      theMan = theMan>>1;
      if (theMan&1<<nbits) theMan = (1<<nbits) - 1;
      if (fFloatValue < 0) theMan |= 1<<(nbits+1);
      *this << theExp;
      *this << theMan;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Write a Double32_t to the buffer.
///
/// The following cases are supported for streaming a Double32_t type
/// depending on the range declaration in the comment field of the data member:
/// Case | Example |
/// -----|---------|
///  A   | Double32_t     fNormal; |
///  B   | Double32_t     fTemperature; //[0,100]|
///  C   | Double32_t     fCharge;      //[-1,1,2]|
///  D   | Double32_t     fVertex[3];   //[-30,30,10]|
///  E   | Double32_t     fChi2;        //[0,0,6]|
///  F   | Int_t         fNsp;<br>Double32_t*    fPointValue;   //[fNsp][0,3]|
///
/// In case A fNormal is converted from a Double_t to a Float_t
/// In case B fTemperature is converted to a 32 bit unsigned integer
/// In case C fCharge is converted to a 2 bits unsigned integer
/// In case D the array elements of fVertex are converted to an unsigned 10 bits integer
/// In case E fChi2 is converted to a Float_t with mantissa truncated precision at 6 bits
/// In case F the fNsp elements of array fPointvalue are converted to an unsigned 32 bit integer
///           Note that the range specifier must follow the dimension specifier.
/// Case B has more precision (9 to 10 significative digits than case A (6 to 7 digits).
/// See TBufferFile::WriteFloat16 for more information.
///
///  see example of use of the Double32_t data type in tutorial double32.C
///  \image html tbufferfile_double32.gif

void TBufferFile::WriteDouble32(Double_t *d, TStreamerElement *ele)
{

   if (ele && ele->GetFactor() != 0) {
      //A range is specified. We normalize the double to the range and
      //convert it to an integer using a scaling factor that is a function of nbits.
      //see TStreamerElement::GetRange.
      Double_t x = d[0];
      Double_t xmin = ele->GetXmin();
      Double_t xmax = ele->GetXmax();
      if (x < xmin) x = xmin;
      if (x > xmax) x = xmax;
      UInt_t aint = UInt_t(0.5+ele->GetFactor()*(x-xmin)); *this << aint;
   } else {
      Int_t nbits = 0;
      //number of bits stored in fXmin (see TStreamerElement::GetRange)
      if (ele) nbits = (Int_t)ele->GetXmin();
      if (!nbits) {
         //if no range and no bits specified, we convert from double to float
         Float_t afloat = (Float_t)d[0];
         *this << afloat;
      } else {
         //a range is not specified, but nbits is.
         //In this case we truncate the mantissa to nbits and we stream
         //the exponent as a UChar_t and the mantissa as a UShort_t.
         union {
            Float_t fFloatValue;
            Int_t   fIntValue;
         };
         fFloatValue = (Float_t)d[0];
         UChar_t  theExp = (UChar_t)(0x000000ff & ((fIntValue<<1)>>24));
         UShort_t theMan = ((1<<(nbits+1))-1) & (fIntValue>>(23-nbits-1)) ;
         theMan++;
         theMan = theMan>>1;
         if (theMan&1<<nbits) theMan = (1<<nbits)-1 ;
         if (fFloatValue < 0) theMan |= 1<<(nbits+1);
         *this << theExp;
         *this << theMan;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of bools from the I/O buffer. Returns the number of
/// bools read. If argument is a 0 pointer then space will be
/// allocated for the array.

Int_t TBufferFile::ReadArray(Bool_t *&b)
{
   R__ASSERT(IsReading());

   Int_t n;
   *this >> n;

   if (n <= 0 || n > fBufSize) return 0;

   if (!b) b = new Bool_t[n];

   if (sizeof(Bool_t) > 1) {
      for (int i = 0; i < n; i++)
         frombuf(fBufCur, &b[i]);
   } else {
      Int_t l = sizeof(Bool_t)*n;
      memcpy(b, fBufCur, l);
      fBufCur += l;
   }

   return n;
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of characters from the I/O buffer. Returns the number of
/// characters read. If argument is a 0 pointer then space will be
/// allocated for the array.

Int_t TBufferFile::ReadArray(Char_t *&c)
{
   R__ASSERT(IsReading());

   Int_t n;
   *this >> n;
   Int_t l = sizeof(Char_t)*n;

   if (l <= 0 || l > fBufSize) return 0;

   if (!c) c = new Char_t[n];

   memcpy(c, fBufCur, l);
   fBufCur += l;

   return n;
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of shorts from the I/O buffer. Returns the number of shorts
/// read. If argument is a 0 pointer then space will be allocated for the
/// array.

Int_t TBufferFile::ReadArray(Short_t *&h)
{
   R__ASSERT(IsReading());

   Int_t n;
   *this >> n;
   Int_t l = sizeof(Short_t)*n;

   if (l <= 0 || l > fBufSize) return 0;

   if (!h) h = new Short_t[n];

#ifdef R__BYTESWAP
# ifdef USE_BSWAPCPY
   bswapcpy16(h, fBufCur, n);
   fBufCur += l;
# else
   for (int i = 0; i < n; i++)
      frombuf(fBufCur, &h[i]);
# endif
#else
   memcpy(h, fBufCur, l);
   fBufCur += l;
#endif

   return n;
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of ints from the I/O buffer. Returns the number of ints
/// read. If argument is a 0 pointer then space will be allocated for the
/// array.

Int_t TBufferFile::ReadArray(Int_t *&ii)
{
   R__ASSERT(IsReading());

   Int_t n;
   *this >> n;
   Int_t l = sizeof(Int_t)*n;

   if (l <= 0 || l > fBufSize) return 0;

   if (!ii) ii = new Int_t[n];

#ifdef R__BYTESWAP
# ifdef USE_BSWAPCPY
   bswapcpy32(ii, fBufCur, n);
   fBufCur += l;
# else
   for (int i = 0; i < n; i++)
      frombuf(fBufCur, &ii[i]);
# endif
#else
   memcpy(ii, fBufCur, l);
   fBufCur += l;
#endif

   return n;
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of longs from the I/O buffer. Returns the number of longs
/// read. If argument is a 0 pointer then space will be allocated for the
/// array.

Int_t TBufferFile::ReadArray(Long_t *&ll)
{
   R__ASSERT(IsReading());

   Int_t n;
   *this >> n;
   Int_t l = sizeof(Long_t)*n;

   if (l <= 0 || l > fBufSize) return 0;

   if (!ll) ll = new Long_t[n];

   TFile *file = (TFile*)fParent;
   if (file && file->GetVersion() < 30006) {
      for (int i = 0; i < n; i++) frombufOld(fBufCur, &ll[i]);
   } else {
      for (int i = 0; i < n; i++) frombuf(fBufCur, &ll[i]);
   }
   return n;
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of long longs from the I/O buffer. Returns the number of
/// long longs read. If argument is a 0 pointer then space will be
/// allocated for the array.

Int_t TBufferFile::ReadArray(Long64_t *&ll)
{
   R__ASSERT(IsReading());

   Int_t n;
   *this >> n;
   Int_t l = sizeof(Long64_t)*n;

   if (l <= 0 || l > fBufSize) return 0;

   if (!ll) ll = new Long64_t[n];

#ifdef R__BYTESWAP
   for (int i = 0; i < n; i++)
      frombuf(fBufCur, &ll[i]);
#else
   memcpy(ll, fBufCur, l);
   fBufCur += l;
#endif

   return n;
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of floats from the I/O buffer. Returns the number of floats
/// read. If argument is a 0 pointer then space will be allocated for the
/// array.

Int_t TBufferFile::ReadArray(Float_t *&f)
{
   R__ASSERT(IsReading());

   Int_t n;
   *this >> n;
   Int_t l = sizeof(Float_t)*n;

   if (l <= 0 || l > fBufSize) return 0;

   if (!f) f = new Float_t[n];

#ifdef R__BYTESWAP
# ifdef USE_BSWAPCPY
   bswapcpy32(f, fBufCur, n);
   fBufCur += l;
# else
   for (int i = 0; i < n; i++)
      frombuf(fBufCur, &f[i]);
# endif
#else
   memcpy(f, fBufCur, l);
   fBufCur += l;
#endif

   return n;
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of doubles from the I/O buffer. Returns the number of doubles
/// read. If argument is a 0 pointer then space will be allocated for the
/// array.

Int_t TBufferFile::ReadArray(Double_t *&d)
{
   R__ASSERT(IsReading());

   Int_t n;
   *this >> n;
   Int_t l = sizeof(Double_t)*n;

   if (l <= 0 || l > fBufSize) return 0;

   if (!d) d = new Double_t[n];

#ifdef R__BYTESWAP
   for (int i = 0; i < n; i++)
      frombuf(fBufCur, &d[i]);
#else
   memcpy(d, fBufCur, l);
   fBufCur += l;
#endif

   return n;
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of floats (written as truncated float) from the I/O buffer.
/// Returns the number of floats read.
/// If argument is a 0 pointer then space will be allocated for the array.
/// see comments about Float16_t encoding at TBufferFile::WriteFloat16

Int_t TBufferFile::ReadArrayFloat16(Float_t *&f, TStreamerElement *ele)
{
   R__ASSERT(IsReading());

   Int_t n;
   *this >> n;

   if (n <= 0 || 3*n > fBufSize) return 0;

   if (!f) f = new Float_t[n];

   ReadFastArrayFloat16(f,n,ele);

   return n;
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of doubles (written as float) from the I/O buffer.
/// Returns the number of doubles read.
/// If argument is a 0 pointer then space will be allocated for the array.
/// see comments about Double32_t encoding at TBufferFile::WriteDouble32

Int_t TBufferFile::ReadArrayDouble32(Double_t *&d, TStreamerElement *ele)
{
   R__ASSERT(IsReading());

   Int_t n;
   *this >> n;

   if (n <= 0 || 3*n > fBufSize) return 0;

   if (!d) d = new Double_t[n];

   ReadFastArrayDouble32(d,n,ele);

   return n;
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of bools from the I/O buffer. Returns the number of bools
/// read.

Int_t TBufferFile::ReadStaticArray(Bool_t *b)
{
   R__ASSERT(IsReading());

   Int_t n;
   *this >> n;

   if (n <= 0 || n > fBufSize) return 0;

   if (!b) return 0;

   if (sizeof(Bool_t) > 1) {
      for (int i = 0; i < n; i++)
         frombuf(fBufCur, &b[i]);
   } else {
      Int_t l = sizeof(Bool_t)*n;
      memcpy(b, fBufCur, l);
      fBufCur += l;
   }

   return n;
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of characters from the I/O buffer. Returns the number of
/// characters read.

Int_t TBufferFile::ReadStaticArray(Char_t *c)
{
   R__ASSERT(IsReading());

   Int_t n;
   *this >> n;
   Int_t l = sizeof(Char_t)*n;

   if (l <= 0 || l > fBufSize) return 0;

   if (!c) return 0;

   memcpy(c, fBufCur, l);
   fBufCur += l;

   return n;
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of shorts from the I/O buffer. Returns the number of shorts
/// read.

Int_t TBufferFile::ReadStaticArray(Short_t *h)
{
   R__ASSERT(IsReading());

   Int_t n;
   *this >> n;
   Int_t l = sizeof(Short_t)*n;

   if (l <= 0 || l > fBufSize) return 0;

   if (!h) return 0;

#ifdef R__BYTESWAP
# ifdef USE_BSWAPCPY
   bswapcpy16(h, fBufCur, n);
   fBufCur += l;
# else
   for (int i = 0; i < n; i++)
      frombuf(fBufCur, &h[i]);
# endif
#else
   memcpy(h, fBufCur, l);
   fBufCur += l;
#endif

   return n;
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of ints from the I/O buffer. Returns the number of ints
/// read.

Int_t TBufferFile::ReadStaticArray(Int_t *ii)
{
   R__ASSERT(IsReading());

   Int_t n;
   *this >> n;
   Int_t l = sizeof(Int_t)*n;

   if (l <= 0 || l > fBufSize) return 0;

   if (!ii) return 0;

#ifdef R__BYTESWAP
# ifdef USE_BSWAPCPY
   bswapcpy32(ii, fBufCur, n);
   fBufCur += sizeof(Int_t)*n;
# else
   for (int i = 0; i < n; i++)
      frombuf(fBufCur, &ii[i]);
# endif
#else
   memcpy(ii, fBufCur, l);
   fBufCur += l;
#endif

   return n;
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of longs from the I/O buffer. Returns the number of longs
/// read.

Int_t TBufferFile::ReadStaticArray(Long_t *ll)
{
   R__ASSERT(IsReading());

   Int_t n;
   *this >> n;
   Int_t l = sizeof(Long_t)*n;

   if (l <= 0 || l > fBufSize) return 0;

   if (!ll) return 0;

   TFile *file = (TFile*)fParent;
   if (file && file->GetVersion() < 30006) {
      for (int i = 0; i < n; i++) frombufOld(fBufCur, &ll[i]);
   } else {
      for (int i = 0; i < n; i++) frombuf(fBufCur, &ll[i]);
   }
   return n;
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of long longs from the I/O buffer. Returns the number of
/// long longs read.

Int_t TBufferFile::ReadStaticArray(Long64_t *ll)
{
   R__ASSERT(IsReading());

   Int_t n;
   *this >> n;
   Int_t l = sizeof(Long64_t)*n;

   if (l <= 0 || l > fBufSize) return 0;

   if (!ll) return 0;

#ifdef R__BYTESWAP
   for (int i = 0; i < n; i++)
      frombuf(fBufCur, &ll[i]);
#else
   memcpy(ll, fBufCur, l);
   fBufCur += l;
#endif

   return n;
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of floats from the I/O buffer. Returns the number of floats
/// read.

Int_t TBufferFile::ReadStaticArray(Float_t *f)
{
   R__ASSERT(IsReading());

   Int_t n;
   *this >> n;
   Int_t l = sizeof(Float_t)*n;

   if (n <= 0 || l > fBufSize) return 0;

   if (!f) return 0;

#ifdef R__BYTESWAP
# ifdef USE_BSWAPCPY
   bswapcpy32(f, fBufCur, n);
   fBufCur += sizeof(Float_t)*n;
# else
   for (int i = 0; i < n; i++)
      frombuf(fBufCur, &f[i]);
# endif
#else
   memcpy(f, fBufCur, l);
   fBufCur += l;
#endif

   return n;
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of doubles from the I/O buffer. Returns the number of doubles
/// read.

Int_t TBufferFile::ReadStaticArray(Double_t *d)
{
   R__ASSERT(IsReading());

   Int_t n;
   *this >> n;
   Int_t l = sizeof(Double_t)*n;

   if (n <= 0 || l > fBufSize) return 0;

   if (!d) return 0;

#ifdef R__BYTESWAP
   for (int i = 0; i < n; i++)
      frombuf(fBufCur, &d[i]);
#else
   memcpy(d, fBufCur, l);
   fBufCur += l;
#endif

   return n;
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of floats (written as truncated float) from the I/O buffer.
/// Returns the number of floats read.
/// see comments about Float16_t encoding at TBufferFile::WriteFloat16

Int_t TBufferFile::ReadStaticArrayFloat16(Float_t *f, TStreamerElement *ele)
{
   R__ASSERT(IsReading());

   Int_t n;
   *this >> n;

   if (n <= 0 || 3*n > fBufSize) return 0;

   if (!f) return 0;

   ReadFastArrayFloat16(f,n,ele);

   return n;
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of doubles (written as float) from the I/O buffer.
/// Returns the number of doubles read.
/// see comments about Double32_t encoding at TBufferFile::WriteDouble32

Int_t TBufferFile::ReadStaticArrayDouble32(Double_t *d, TStreamerElement *ele)
{
   R__ASSERT(IsReading());

   Int_t n;
   *this >> n;

   if (n <= 0 || 3*n > fBufSize) return 0;

   if (!d) return 0;

   ReadFastArrayDouble32(d,n,ele);

   return n;
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of n bools from the I/O buffer.

void TBufferFile::ReadFastArray(Bool_t *b, Int_t n)
{
   if (n <= 0 || n > fBufSize) return;

   if (sizeof(Bool_t) > 1) {
      for (int i = 0; i < n; i++)
         frombuf(fBufCur, &b[i]);
   } else {
      Int_t l = sizeof(Bool_t)*n;
      memcpy(b, fBufCur, l);
      fBufCur += l;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of n characters from the I/O buffer.

void TBufferFile::ReadFastArray(Char_t *c, Int_t n)
{
   if (n <= 0 || n > fBufSize) return;

   Int_t l = sizeof(Char_t)*n;
   memcpy(c, fBufCur, l);
   fBufCur += l;
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of n characters from the I/O buffer.

void TBufferFile::ReadFastArrayString(Char_t *c, Int_t n)
{
   Int_t len;
   UChar_t lenchar;
   *this >> lenchar;
   if (lenchar < 255) {
      len = lenchar;
   } else {
      *this >> len;
   }
   if (len) {
      if (len <= 0 || len > fBufSize) return;
      Int_t blen = len;
      if (len >= n) len = n-1;

      Int_t l = sizeof(Char_t)*len;
      memcpy(c, fBufCur, l);
      fBufCur += blen;

      c[len] = 0;
   } else {
      c[0] = 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of n shorts from the I/O buffer.

void TBufferFile::ReadFastArray(Short_t *h, Int_t n)
{
   Int_t l = sizeof(Short_t)*n;
   if (n <= 0 || l > fBufSize) return;

#ifdef R__BYTESWAP
# ifdef USE_BSWAPCPY
   bswapcpy16(h, fBufCur, n);
   fBufCur += sizeof(Short_t)*n;
# else
   for (int i = 0; i < n; i++)
      frombuf(fBufCur, &h[i]);
# endif
#else
   memcpy(h, fBufCur, l);
   fBufCur += l;
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of n ints from the I/O buffer.

void TBufferFile::ReadFastArray(Int_t *ii, Int_t n)
{
   Int_t l = sizeof(Int_t)*n;
   if (l <= 0 || l > fBufSize) return;

#ifdef R__BYTESWAP
# ifdef USE_BSWAPCPY
   bswapcpy32(ii, fBufCur, n);
   fBufCur += sizeof(Int_t)*n;
# else
   for (int i = 0; i < n; i++)
      frombuf(fBufCur, &ii[i]);
# endif
#else
   memcpy(ii, fBufCur, l);
   fBufCur += l;
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of n longs from the I/O buffer.

void TBufferFile::ReadFastArray(Long_t *ll, Int_t n)
{
   Int_t l = sizeof(Long_t)*n;
   if (l <= 0 || l > fBufSize) return;

   TFile *file = (TFile*)fParent;
   if (file && file->GetVersion() < 30006) {
      for (int i = 0; i < n; i++) frombufOld(fBufCur, &ll[i]);
   } else {
      for (int i = 0; i < n; i++) frombuf(fBufCur, &ll[i]);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of n long longs from the I/O buffer.

void TBufferFile::ReadFastArray(Long64_t *ll, Int_t n)
{
   Int_t l = sizeof(Long64_t)*n;
   if (l <= 0 || l > fBufSize) return;

#ifdef R__BYTESWAP
   for (int i = 0; i < n; i++)
      frombuf(fBufCur, &ll[i]);
#else
   memcpy(ll, fBufCur, l);
   fBufCur += l;
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of n floats from the I/O buffer.

void TBufferFile::ReadFastArray(Float_t *f, Int_t n)
{
   Int_t l = sizeof(Float_t)*n;
   if (l <= 0 || l > fBufSize) return;

#ifdef R__BYTESWAP
# ifdef USE_BSWAPCPY
   bswapcpy32(f, fBufCur, n);
   fBufCur += sizeof(Float_t)*n;
# else
   for (int i = 0; i < n; i++)
      frombuf(fBufCur, &f[i]);
# endif
#else
   memcpy(f, fBufCur, l);
   fBufCur += l;
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of n doubles from the I/O buffer.

void TBufferFile::ReadFastArray(Double_t *d, Int_t n)
{
   Int_t l = sizeof(Double_t)*n;
   if (l <= 0 || l > fBufSize) return;

#ifdef R__BYTESWAP
   for (int i = 0; i < n; i++)
      frombuf(fBufCur, &d[i]);
#else
   memcpy(d, fBufCur, l);
   fBufCur += l;
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of n floats (written as truncated float) from the I/O buffer.
/// see comments about Float16_t encoding at TBufferFile::WriteFloat16

void TBufferFile::ReadFastArrayFloat16(Float_t *f, Int_t n, TStreamerElement *ele)
{
   if (n <= 0 || 3*n > fBufSize) return;

   if (ele && ele->GetFactor() != 0) {
      //a range was specified. We read an integer and convert it back to a float
      Double_t xmin = ele->GetXmin();
      Double_t factor = ele->GetFactor();
      for (int j=0;j < n; j++) {
         UInt_t aint; *this >> aint; f[j] = (Float_t)(aint/factor + xmin);
      }
   } else {
      Int_t i;
      Int_t nbits = 0;
      if (ele) nbits = (Int_t)ele->GetXmin();
      if (!nbits) nbits = 12;
      //we read the exponent and the truncated mantissa of the float
      //and rebuild the new float.
      union {
         Float_t fFloatValue;
         Int_t   fIntValue;
      };
      UChar_t  theExp;
      UShort_t theMan;
      for (i = 0; i < n; i++) {
         *this >> theExp;
         *this >> theMan;
         fIntValue = theExp;
         fIntValue <<= 23;
         fIntValue |= (theMan & ((1<<(nbits+1))-1)) <<(23-nbits);
         if(1<<(nbits+1) & theMan) fFloatValue = -fFloatValue;
         f[i] = fFloatValue;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of n floats (written as truncated float) from the I/O buffer.
/// see comments about Float16_t encoding at TBufferFile::WriteFloat16

void TBufferFile::ReadFastArrayWithFactor(Float_t *ptr, Int_t n, Double_t factor, Double_t minvalue)
{
   if (n <= 0 || 3*n > fBufSize) return;

   //a range was specified. We read an integer and convert it back to a float
   for (int j=0;j < n; j++) {
      UInt_t aint; *this >> aint; ptr[j] = (Float_t)(aint/factor + minvalue);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of n floats (written as truncated float) from the I/O buffer.
/// see comments about Float16_t encoding at TBufferFile::WriteFloat16

void TBufferFile::ReadFastArrayWithNbits(Float_t *ptr, Int_t n, Int_t nbits)
{
   if (n <= 0 || 3*n > fBufSize) return;

   if (!nbits) nbits = 12;
   //we read the exponent and the truncated mantissa of the float
   //and rebuild the new float.
   union {
      Float_t fFloatValue;
      Int_t   fIntValue;
   };
   UChar_t  theExp;
   UShort_t theMan;
   for (Int_t i = 0; i < n; i++) {
      *this >> theExp;
      *this >> theMan;
      fIntValue = theExp;
      fIntValue <<= 23;
      fIntValue |= (theMan & ((1<<(nbits+1))-1)) <<(23-nbits);
      if(1<<(nbits+1) & theMan) fFloatValue = -fFloatValue;
      ptr[i] = fFloatValue;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of n doubles (written as float) from the I/O buffer.
/// see comments about Double32_t encoding at TBufferFile::WriteDouble32

void TBufferFile::ReadFastArrayDouble32(Double_t *d, Int_t n, TStreamerElement *ele)
{
   if (n <= 0 || 3*n > fBufSize) return;

   if (ele && ele->GetFactor() != 0) {
      //a range was specified. We read an integer and convert it back to a double.
      Double_t xmin = ele->GetXmin();
      Double_t factor = ele->GetFactor();
      for (int j=0;j < n; j++) {
         UInt_t aint; *this >> aint; d[j] = (Double_t)(aint/factor + xmin);
      }
   } else {
      Int_t i;
      Int_t nbits = 0;
      if (ele) nbits = (Int_t)ele->GetXmin();
      if (!nbits) {
         //we read a float and convert it to double
         Float_t afloat;
         for (i = 0; i < n; i++) {
            *this >> afloat;
            d[i] = (Double_t)afloat;
         }
      } else {
         //we read the exponent and the truncated mantissa of the float
         //and rebuild the double.
         union {
            Float_t fFloatValue;
            Int_t   fIntValue;
         };
         UChar_t  theExp;
         UShort_t theMan;
         for (i = 0; i < n; i++) {
            *this >> theExp;
            *this >> theMan;
            fIntValue = theExp;
            fIntValue <<= 23;
            fIntValue |= (theMan & ((1<<(nbits+1))-1)) <<(23-nbits);
            if (1<<(nbits+1) & theMan) fFloatValue = -fFloatValue;
            d[i] = (Double_t)fFloatValue;
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of n doubles (written as float) from the I/O buffer.
/// see comments about Double32_t encoding at TBufferFile::WriteDouble32

void TBufferFile::ReadFastArrayWithFactor(Double_t *d, Int_t n, Double_t factor, Double_t minvalue)
{
   if (n <= 0 || 3*n > fBufSize) return;

   //a range was specified. We read an integer and convert it back to a double.
   for (int j=0;j < n; j++) {
      UInt_t aint; *this >> aint; d[j] = (Double_t)(aint/factor + minvalue);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Read array of n doubles (written as float) from the I/O buffer.
/// see comments about Double32_t encoding at TBufferFile::WriteDouble32

void TBufferFile::ReadFastArrayWithNbits(Double_t *d, Int_t n, Int_t nbits)
{
   if (n <= 0 || 3*n > fBufSize) return;

   if (!nbits) {
      //we read a float and convert it to double
      Float_t afloat;
      for (Int_t i = 0; i < n; i++) {
         *this >> afloat;
         d[i] = (Double_t)afloat;
      }
   } else {
      //we read the exponent and the truncated mantissa of the float
      //and rebuild the double.
      union {
         Float_t fFloatValue;
         Int_t   fIntValue;
      };
      UChar_t  theExp;
      UShort_t theMan;
      for (Int_t i = 0; i < n; i++) {
         *this >> theExp;
         *this >> theMan;
         fIntValue = theExp;
         fIntValue <<= 23;
         fIntValue |= (theMan & ((1<<(nbits+1))-1)) <<(23-nbits);
         if (1<<(nbits+1) & theMan) fFloatValue = -fFloatValue;
         d[i] = (Double_t)fFloatValue;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Read an array of 'n' objects from the I/O buffer.
/// Stores the objects read starting at the address 'start'.
/// The objects in the array are assume to be of class 'cl'.

void TBufferFile::ReadFastArray(void  *start, const TClass *cl, Int_t n,
                                TMemberStreamer *streamer, const TClass* onFileClass )
{
   if (streamer) {
      streamer->SetOnFileClass(onFileClass);
      (*streamer)(*this,start,0);
      return;
   }

   int objectSize = cl->Size();
   char *obj = (char*)start;
   char *end = obj + n*objectSize;

   for(; obj<end; obj+=objectSize) ((TClass*)cl)->Streamer(obj,*this, onFileClass);
}

////////////////////////////////////////////////////////////////////////////////
/// Read an array of 'n' objects from the I/O buffer.
///
/// The objects read are stored starting at the address '*start'
/// The objects in the array are assumed to be of class 'cl' or a derived class.
/// 'mode' indicates whether the data member is marked with '->'

void TBufferFile::ReadFastArray(void **start, const TClass *cl, Int_t n,
                                Bool_t isPreAlloc, TMemberStreamer *streamer, const TClass* onFileClass)
{
   // if isPreAlloc is true (data member has a ->) we can assume that the pointer (*start)
   // is never 0.

   if (streamer) {
      if (isPreAlloc) {
         for (Int_t j=0;j<n;j++) {
            if (!start[j]) start[j] = cl->New();
         }
      }
      streamer->SetOnFileClass(onFileClass);
      (*streamer)(*this,(void*)start,0);
      return;
   }

   if (!isPreAlloc) {

      for (Int_t j=0; j<n; j++){
         //delete the object or collection
         void *old = start[j];
         start[j] = ReadObjectAny(cl);
         if (old && old!=start[j] &&
             TStreamerInfo::CanDelete()
             // There are some cases where the user may set up a pointer in the (default)
             // constructor but not mark this pointer as transient.  Sometime the value
             // of this pointer is the address of one of the object with just created
             // and the following delete would result in the deletion (possibly of the
             // top level object we are goint to return!).
             // Eventhough this is a user error, we could prevent the crash by simply
             // adding:
             // && !CheckObject(start[j],cl)
             // However this can increase the read time significantly (10% in the case
             // of one TLine pointer in the test/Track and run ./Event 200 0 0 20 30000
             //
             // If ReadObjectAny returned the same value as we previous had, this means
             // that when writing this object (start[j] had already been written and
             // is indeed pointing to the same object as the object the user set up
             // in the default constructor).
             ) {
            ((TClass*)cl)->Destructor(old,kFALSE); // call delete and desctructor
         }
      }

   } else {
      //case //-> in comment

      for (Int_t j=0; j<n; j++){
         if (!start[j]) start[j] = ((TClass*)cl)->New();
         ((TClass*)cl)->Streamer(start[j],*this,onFileClass);
      }

   }
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of n bools into the I/O buffer.

void TBufferFile::WriteArray(const Bool_t *b, Int_t n)
{
   R__ASSERT(IsWriting());

   *this << n;

   if (n <= 0) return;

   R__ASSERT(b);

   Int_t l = sizeof(UChar_t)*n;
   if (fBufCur + l > fBufMax) AutoExpand(fBufSize+l);

   if (sizeof(Bool_t) > 1) {
      for (int i = 0; i < n; i++)
         tobuf(fBufCur, b[i]);
   } else {
      memcpy(fBufCur, b, l);
      fBufCur += l;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of n characters into the I/O buffer.

void TBufferFile::WriteArray(const Char_t *c, Int_t n)
{
   R__ASSERT(IsWriting());

   *this << n;

   if (n <= 0) return;

   R__ASSERT(c);

   Int_t l = sizeof(Char_t)*n;
   if (fBufCur + l > fBufMax) AutoExpand(fBufSize+l);

   memcpy(fBufCur, c, l);
   fBufCur += l;
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of n shorts into the I/O buffer.

void TBufferFile::WriteArray(const Short_t *h, Int_t n)
{
   R__ASSERT(IsWriting());

   *this << n;

   if (n <= 0) return;

   R__ASSERT(h);

   Int_t l = sizeof(Short_t)*n;
   if (fBufCur + l > fBufMax) AutoExpand(fBufSize+l);

#ifdef R__BYTESWAP
# ifdef USE_BSWAPCPY
   bswapcpy16(fBufCur, h, n);
   fBufCur += l;
# else
   for (int i = 0; i < n; i++)
      tobuf(fBufCur, h[i]);
# endif
#else
   memcpy(fBufCur, h, l);
   fBufCur += l;
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of n ints into the I/O buffer.

void TBufferFile::WriteArray(const Int_t *ii, Int_t n)
{
   R__ASSERT(IsWriting());

   *this << n;

   if (n <= 0) return;

   R__ASSERT(ii);

   Int_t l = sizeof(Int_t)*n;
   if (fBufCur + l > fBufMax) AutoExpand(fBufSize+l);

#ifdef R__BYTESWAP
# ifdef USE_BSWAPCPY
   bswapcpy32(fBufCur, ii, n);
   fBufCur += l;
# else
   for (int i = 0; i < n; i++)
      tobuf(fBufCur, ii[i]);
# endif
#else
   memcpy(fBufCur, ii, l);
   fBufCur += l;
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of n longs into the I/O buffer.

void TBufferFile::WriteArray(const Long_t *ll, Int_t n)
{
   R__ASSERT(IsWriting());

   *this << n;

   if (n <= 0) return;

   R__ASSERT(ll);

   Int_t l = 8*n;
   if (fBufCur + l > fBufMax) AutoExpand(fBufSize+l);
   for (int i = 0; i < n; i++) tobuf(fBufCur, ll[i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of n unsigned longs into the I/O buffer.
/// This is an explicit case for unsigned longs since signed longs
/// have a special tobuf().

void TBufferFile::WriteArray(const ULong_t *ll, Int_t n)
{
   R__ASSERT(IsWriting());

   *this << n;

   if (n <= 0) return;

   R__ASSERT(ll);

   Int_t l = 8*n;
   if (fBufCur + l > fBufMax) AutoExpand(fBufSize+l);
   for (int i = 0; i < n; i++) tobuf(fBufCur, ll[i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of n long longs into the I/O buffer.

void TBufferFile::WriteArray(const Long64_t *ll, Int_t n)
{
   R__ASSERT(IsWriting());

   *this << n;

   if (n <= 0) return;

   R__ASSERT(ll);

   Int_t l = sizeof(Long64_t)*n;
   if (fBufCur + l > fBufMax) AutoExpand(fBufSize+l);

#ifdef R__BYTESWAP
   for (int i = 0; i < n; i++)
      tobuf(fBufCur, ll[i]);
#else
   memcpy(fBufCur, ll, l);
   fBufCur += l;
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of n floats into the I/O buffer.

void TBufferFile::WriteArray(const Float_t *f, Int_t n)
{
   R__ASSERT(IsWriting());

   *this << n;

   if (n <= 0) return;

   R__ASSERT(f);

   Int_t l = sizeof(Float_t)*n;
   if (fBufCur + l > fBufMax) AutoExpand(fBufSize+l);

#ifdef R__BYTESWAP
# ifdef USE_BSWAPCPY
   bswapcpy32(fBufCur, f, n);
   fBufCur += l;
# else
   for (int i = 0; i < n; i++)
      tobuf(fBufCur, f[i]);
# endif
#else
   memcpy(fBufCur, f, l);
   fBufCur += l;
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of n doubles into the I/O buffer.

void TBufferFile::WriteArray(const Double_t *d, Int_t n)
{
   R__ASSERT(IsWriting());

   *this << n;

   if (n <= 0) return;

   R__ASSERT(d);

   Int_t l = sizeof(Double_t)*n;
   if (fBufCur + l > fBufMax) AutoExpand(fBufSize+l);

#ifdef R__BYTESWAP
   for (int i = 0; i < n; i++)
      tobuf(fBufCur, d[i]);
#else
   memcpy(fBufCur, d, l);
   fBufCur += l;
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of n floats (as truncated float) into the I/O buffer.
/// see comments about Float16_t encoding at TBufferFile::WriteFloat16

void TBufferFile::WriteArrayFloat16(const Float_t *f, Int_t n, TStreamerElement *ele)
{
   R__ASSERT(IsWriting());

   *this << n;

   if (n <= 0) return;

   R__ASSERT(f);

   Int_t l = sizeof(Float_t)*n;
   if (fBufCur + l > fBufMax) AutoExpand(fBufSize+l);

   WriteFastArrayFloat16(f,n,ele);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of n doubles (as float) into the I/O buffer.
/// see comments about Double32_t encoding at TBufferFile::WriteDouble32

void TBufferFile::WriteArrayDouble32(const Double_t *d, Int_t n, TStreamerElement *ele)
{
   R__ASSERT(IsWriting());

   *this << n;

   if (n <= 0) return;

   R__ASSERT(d);

   Int_t l = sizeof(Float_t)*n;
   if (fBufCur + l > fBufMax) AutoExpand(fBufSize+l);

   WriteFastArrayDouble32(d,n,ele);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of n bools into the I/O buffer.

void TBufferFile::WriteFastArray(const Bool_t *b, Int_t n)
{
   if (n <= 0) return;

   Int_t l = sizeof(UChar_t)*n;
   if (fBufCur + l > fBufMax) AutoExpand(fBufSize+l);

   if (sizeof(Bool_t) > 1) {
      for (int i = 0; i < n; i++)
         tobuf(fBufCur, b[i]);
   } else {
      memcpy(fBufCur, b, l);
      fBufCur += l;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of n characters into the I/O buffer.

void TBufferFile::WriteFastArray(const Char_t *c, Int_t n)
{
   if (n <= 0) return;

   Int_t l = sizeof(Char_t)*n;
   if (fBufCur + l > fBufMax) AutoExpand(fBufSize+l);

   memcpy(fBufCur, c, l);
   fBufCur += l;
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of n characters into the I/O buffer.

void TBufferFile::WriteFastArrayString(const Char_t *c, Int_t n)
{
   if (n < 255) {
      *this << (UChar_t)n;
   } else {
      *this << (UChar_t)255;
      *this << n;
   }

   if (n <= 0) return;

   Int_t l = sizeof(Char_t)*n;
   if (fBufCur + l > fBufMax) AutoExpand(fBufSize+l);

   memcpy(fBufCur, c, l);
   fBufCur += l;
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of n shorts into the I/O buffer.

void TBufferFile::WriteFastArray(const Short_t *h, Int_t n)
{
   if (n <= 0) return;

   Int_t l = sizeof(Short_t)*n;
   if (fBufCur + l > fBufMax) AutoExpand(fBufSize+l);

#ifdef R__BYTESWAP
# ifdef USE_BSWAPCPY
   bswapcpy16(fBufCur, h, n);
   fBufCur += l;
# else
   for (int i = 0; i < n; i++)
      tobuf(fBufCur, h[i]);
# endif
#else
   memcpy(fBufCur, h, l);
   fBufCur += l;
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of n ints into the I/O buffer.

void TBufferFile::WriteFastArray(const Int_t *ii, Int_t n)
{
   if (n <= 0) return;

   Int_t l = sizeof(Int_t)*n;
   if (fBufCur + l > fBufMax) AutoExpand(fBufSize+l);

#ifdef R__BYTESWAP
# ifdef USE_BSWAPCPY
   bswapcpy32(fBufCur, ii, n);
   fBufCur += l;
# else
   for (int i = 0; i < n; i++)
      tobuf(fBufCur, ii[i]);
# endif
#else
   memcpy(fBufCur, ii, l);
   fBufCur += l;
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of n longs into the I/O buffer.

void TBufferFile::WriteFastArray(const Long_t *ll, Int_t n)
{
   if (n <= 0) return;

   Int_t l = 8*n;
   if (fBufCur + l > fBufMax) AutoExpand(fBufSize+l);

   for (int i = 0; i < n; i++) tobuf(fBufCur, ll[i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of n unsigned longs into the I/O buffer.
/// This is an explicit case for unsigned longs since signed longs
/// have a special tobuf().

void TBufferFile::WriteFastArray(const ULong_t *ll, Int_t n)
{
   if (n <= 0) return;

   Int_t l = 8*n;
   if (fBufCur + l > fBufMax) AutoExpand(fBufSize+l);

   for (int i = 0; i < n; i++) tobuf(fBufCur, ll[i]);
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of n long longs into the I/O buffer.

void TBufferFile::WriteFastArray(const Long64_t *ll, Int_t n)
{
   if (n <= 0) return;

   Int_t l = sizeof(Long64_t)*n;
   if (fBufCur + l > fBufMax) AutoExpand(fBufSize+l);

#ifdef R__BYTESWAP
   for (int i = 0; i < n; i++)
      tobuf(fBufCur, ll[i]);
#else
   memcpy(fBufCur, ll, l);
   fBufCur += l;
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of n floats into the I/O buffer.

void TBufferFile::WriteFastArray(const Float_t *f, Int_t n)
{
   if (n <= 0) return;

   Int_t l = sizeof(Float_t)*n;
   if (fBufCur + l > fBufMax) AutoExpand(fBufSize+l);

#ifdef R__BYTESWAP
# ifdef USE_BSWAPCPY
   bswapcpy32(fBufCur, f, n);
   fBufCur += l;
# else
   for (int i = 0; i < n; i++)
      tobuf(fBufCur, f[i]);
# endif
#else
   memcpy(fBufCur, f, l);
   fBufCur += l;
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of n doubles into the I/O buffer.

void TBufferFile::WriteFastArray(const Double_t *d, Int_t n)
{
   if (n <= 0) return;

   Int_t l = sizeof(Double_t)*n;
   if (fBufCur + l > fBufMax) AutoExpand(fBufSize+l);

#ifdef R__BYTESWAP
   for (int i = 0; i < n; i++)
      tobuf(fBufCur, d[i]);
#else
   memcpy(fBufCur, d, l);
   fBufCur += l;
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of n floats (as truncated float) into the I/O buffer.
/// see comments about Float16_t encoding at TBufferFile::WriteFloat16

void TBufferFile::WriteFastArrayFloat16(const Float_t *f, Int_t n, TStreamerElement *ele)
{
   if (n <= 0) return;

   Int_t l = sizeof(Float_t)*n;
   if (fBufCur + l > fBufMax) AutoExpand(fBufSize+l);

   if (ele && ele->GetFactor()) {
      //A range is specified. We normalize the float to the range and
      //convert it to an integer using a scaling factor that is a function of nbits.
      //see TStreamerElement::GetRange.
      Double_t factor = ele->GetFactor();
      Double_t xmin = ele->GetXmin();
      Double_t xmax = ele->GetXmax();
      for (int j = 0; j < n; j++) {
         Float_t x = f[j];
         if (x < xmin) x = xmin;
         if (x > xmax) x = xmax;
         UInt_t aint = UInt_t(0.5+factor*(x-xmin)); *this << aint;
      }
   } else {
      Int_t nbits = 0;
      //number of bits stored in fXmin (see TStreamerElement::GetRange)
      if (ele) nbits = (Int_t)ele->GetXmin();
      if (!nbits) nbits = 12;
      Int_t i;
      //a range is not specified, but nbits is.
      //In this case we truncate the mantissa to nbits and we stream
      //the exponent as a UChar_t and the mantissa as a UShort_t.
      union {
         Float_t fFloatValue;
         Int_t   fIntValue;
      };
      for (i = 0; i < n; i++) {
         fFloatValue = f[i];
         UChar_t  theExp = (UChar_t)(0x000000ff & ((fIntValue<<1)>>24));
         UShort_t theMan = ((1<<(nbits+1))-1) & (fIntValue>>(23-nbits-1));
         theMan++;
         theMan = theMan>>1;
         if (theMan&1<<nbits) theMan = (1<<nbits) - 1;
         if (fFloatValue < 0) theMan |= 1<<(nbits+1);
         *this << theExp;
         *this << theMan;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Write array of n doubles (as float) into the I/O buffer.
/// see comments about Double32_t encoding at TBufferFile::WriteDouble32

void TBufferFile::WriteFastArrayDouble32(const Double_t *d, Int_t n, TStreamerElement *ele)
{
   if (n <= 0) return;

   Int_t l = sizeof(Float_t)*n;
   if (fBufCur + l > fBufMax) AutoExpand(fBufSize+l);

   if (ele && ele->GetFactor()) {
      //A range is specified. We normalize the double to the range and
      //convert it to an integer using a scaling factor that is a function of nbits.
      //see TStreamerElement::GetRange.
      Double_t factor = ele->GetFactor();
      Double_t xmin = ele->GetXmin();
      Double_t xmax = ele->GetXmax();
      for (int j = 0; j < n; j++) {
         Double_t x = d[j];
         if (x < xmin) x = xmin;
         if (x > xmax) x = xmax;
         UInt_t aint = UInt_t(0.5+factor*(x-xmin)); *this << aint;
      }
   } else {
      Int_t nbits = 0;
      //number of bits stored in fXmin (see TStreamerElement::GetRange)
      if (ele) nbits = (Int_t)ele->GetXmin();
      Int_t i;
      if (!nbits) {
         //if no range and no bits specified, we convert from double to float
         for (i = 0; i < n; i++) {
            Float_t afloat = (Float_t)d[i];
            *this << afloat;
         }
      } else {
         //a range is not specified, but nbits is.
         //In this case we truncate the mantissa to nbits and we stream
         //the exponent as a UChar_t and the mantissa as a UShort_t.
         union {
            Float_t fFloatValue;
            Int_t   fIntValue;
         };
         for (i = 0; i < n; i++) {
            fFloatValue = (Float_t)d[i];
            UChar_t  theExp = (UChar_t)(0x000000ff & ((fIntValue<<1)>>24));
            UShort_t theMan = ((1<<(nbits+1))-1) & (fIntValue>>(23-nbits-1));
            theMan++;
            theMan = theMan>>1;
            if(theMan&1<<nbits) theMan = (1<<nbits) - 1;
            if (fFloatValue < 0) theMan |= 1<<(nbits+1);
            *this << theExp;
            *this << theMan;
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Write an array of object starting at the address 'start' and of length 'n'
/// the objects in the array are assumed to be of class 'cl'

void TBufferFile::WriteFastArray(void  *start, const TClass *cl, Int_t n,
                                 TMemberStreamer *streamer)
{
   if (streamer) {
      (*streamer)(*this, start, 0);
      return;
   }

   char *obj = (char*)start;
   if (!n) n=1;
   int size = cl->Size();

   for(Int_t j=0; j<n; j++,obj+=size) {
      ((TClass*)cl)->Streamer(obj,*this);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Write an array of object starting at the address '*start' and of length 'n'
/// the objects in the array are of class 'cl'
/// 'isPreAlloc' indicates whether the data member is marked with '->'
/// Return:
///   - 0: success
///   - 2: truncated success (i.e actual class is missing. Only ptrClass saved.)

Int_t TBufferFile::WriteFastArray(void **start, const TClass *cl, Int_t n,
                                  Bool_t isPreAlloc, TMemberStreamer *streamer)
{
   // if isPreAlloc is true (data member has a ->) we can assume that the pointer
   // is never 0.

   if (streamer) {
      (*streamer)(*this,(void*)start,0);
      return 0;
   }

   int strInfo = 0;

   Int_t res = 0;

   if (!isPreAlloc) {

      for (Int_t j=0;j<n;j++) {
         //must write StreamerInfo if pointer is null
         if (!strInfo && !start[j]) {
            if (cl->Property() & kIsAbstract) {
               // Do not try to generate the StreamerInfo for an abstract class
            } else {
               TStreamerInfo *info = (TStreamerInfo*)((TClass*)cl)->GetStreamerInfo();
               ForceWriteInfo(info,kFALSE);
            }
         }
         strInfo = 2003;
         res |= WriteObjectAny(start[j],cl);
      }

   } else {
      //case //-> in comment

      for (Int_t j=0;j<n;j++) {
         if (!start[j]) start[j] = ((TClass*)cl)->New();
         ((TClass*)cl)->Streamer(start[j],*this);
      }

   }
   return res;
}

////////////////////////////////////////////////////////////////////////////////
/// Read object from I/O buffer. clReq is NOT used.
///
/// The value returned is the address of the actual start in memory of
/// the object. Note that if the actual class of the object does not
/// inherit first from TObject, the type of the pointer is NOT 'TObject*'.
/// [More accurately, the class needs to start with the TObject part, for
/// the pointer to be a real TObject*].
/// We recommend using ReadObjectAny instead of ReadObject

TObject *TBufferFile::ReadObject(const TClass * /*clReq*/)
{
   return (TObject*) ReadObjectAny(nullptr);
}

////////////////////////////////////////////////////////////////////////////////
/// Skip any kind of object from buffer

void TBufferFile::SkipObjectAny()
{
   UInt_t start, count;
   ReadVersion(&start, &count);
   SetBufferOffset(start+count+sizeof(UInt_t));
}

////////////////////////////////////////////////////////////////////////////////
/// Read object from I/O buffer.
///
/// A typical use for this function is:
///
///     MyClass *ptr = (MyClass*)b.ReadObjectAny(MyClass::Class());
///
/// I.e. clCast should point to a TClass object describing the class pointed
/// to by your pointer.
/// In case of multiple inheritance, the return value might not be the
/// real beginning of the object in memory.  You will need to use a
/// dynamic_cast later if you need to retrieve it.

void *TBufferFile::ReadObjectAny(const TClass *clCast)
{
   R__ASSERT(IsReading());

   // make sure fMap is initialized
   InitMap();

   // before reading object save start position
   UInt_t startpos = UInt_t(fBufCur-fBuffer);

   // attempt to load next object as TClass clCast
   UInt_t tag;       // either tag or byte count
   TClass *clRef = ReadClass(clCast, &tag);
   TClass *clOnfile = nullptr;
   Int_t baseOffset = 0;
   if (clRef && (clRef!=(TClass*)(-1)) && clCast) {
      //baseOffset will be -1 if clRef does not inherit from clCast.
      baseOffset = clRef->GetBaseClassOffset(clCast);
      if (baseOffset == -1) {
         // The 2 classes are unrelated, maybe there is a converter between the 2.

         if (!clCast->GetSchemaRules() ||
             !clCast->GetSchemaRules()->HasRuleWithSourceClass(clRef->GetName()))
         {
            // There is no converter
            Error("ReadObject", "got object of wrong class! requested %s but got %s",
                  clCast->GetName(), clRef->GetName());

            CheckByteCount(startpos, tag, (TClass *)nullptr); // avoid mis-leading byte count error message
            return 0; // We better return at this point
         }
         baseOffset = 0; // For now we do not support requesting from a class that is the base of one of the class for which there is transformation to ....

         Info("ReadObjectAny","Using Converter StreamerInfo from %s to %s",clRef->GetName(),clCast->GetName());
         clRef = const_cast<TClass*>(clCast);

      }
      if (clCast->GetState() > TClass::kEmulated && clRef->GetState() <= TClass::kEmulated) {
         //we cannot mix a compiled class with an emulated class in the inheritance
         Error("ReadObject", "trying to read an emulated class (%s) to store in a compiled pointer (%s)",
               clRef->GetName(),clCast->GetName());
         CheckByteCount(startpos, tag, (TClass *)nullptr); // avoid mis-leading byte count error message
         return 0;
      }
   }

   // check if object has not already been read
   // (this can only happen when called via CheckObject())
   char *obj;
   if (fVersion > 0) {
      obj = (char *) (Long_t)fMap->GetValue(startpos+kMapOffset);
      if (obj == (void*) -1) obj = nullptr;
      if (obj) {
         CheckByteCount(startpos, tag, (TClass *)nullptr);
         return (obj + baseOffset);
      }
   }

   // unknown class, skip to next object and return 0 obj
   if (clRef == (TClass*) -1) {
      if (fBufCur >= fBufMax) return 0;
      if (fVersion > 0)
         MapObject((TObject*) -1, startpos+kMapOffset);
      else
         MapObject((void*)nullptr, nullptr, fMapCount);
      CheckByteCount(startpos, tag, (TClass *)nullptr);
      return 0;
   }

   if (!clRef) {

      // got a reference to an already read object
      if (fVersion > 0) {
         tag += fDisplacement;
         tag = CheckObject(tag, clCast);
      } else {
         if (tag > (UInt_t)fMap->GetSize()) {
            Error("ReadObject", "object tag too large, I/O buffer corrupted");
            return 0;
            // exception
         }
      }
      obj = (char *) (Long_t)fMap->GetValue(tag);
      clRef = (TClass*) (Long_t)fClassMap->GetValue(tag);

      if (clRef && (clRef!=(TClass*)(-1)) && clCast) {
         //baseOffset will be -1 if clRef does not inherit from clCast.
         baseOffset = clRef->GetBaseClassOffset(clCast);
         if (baseOffset == -1) {
            Error("ReadObject", "Got object of wrong class (Got %s while expecting %s)",
                  clRef->GetName(),clCast->GetName());
            // exception
            baseOffset = 0;
         }
      }

      // There used to be a warning printed here when:
      //   obj && isTObject && !((TObject*)obj)->IsA()->InheritsFrom(clReq)
      // however isTObject was based on clReq (now clCast).
      // If the test was to fail, then it is as likely that the object is not a TObject
      // and then we have a potential core dump.
      // At this point (missing clRef), we do NOT have enough information to really
      // answer the question: is the object read of the type I requested.

   } else {

      // allocate a new object based on the class found
      obj = (char*)clRef->New();
      if (!obj) {
         Error("ReadObject", "could not create object of class %s",
               clRef->GetName());
         // exception
         return 0;
      }

      // add to fMap before reading rest of object
      if (fVersion > 0)
         MapObject(obj, clRef, startpos+kMapOffset);
      else
         MapObject(obj, clRef, fMapCount);

      // let the object read itself
      clRef->Streamer( obj, *this, clOnfile );

      CheckByteCount(startpos, tag, clRef);
   }

   return obj+baseOffset;
}

////////////////////////////////////////////////////////////////////////////////
/// Write object to I/O buffer.
/// This function assumes that the value of 'actualObjectStart' is the actual start of
/// the object of class 'actualClass'
/// If 'cacheReuse' is true (default) upon seeing an object address a second time,
/// we record the offset where its was written the first time rather than streaming
/// the object a second time.
/// If 'cacheReuse' is false, we always stream the object.  This allows the (re)use
/// of temporary object to store different data in the same buffer.

void TBufferFile::WriteObjectClass(const void *actualObjectStart, const TClass *actualClass, Bool_t cacheReuse)
{
   R__ASSERT(IsWriting());

   if (!actualObjectStart) {

      // save kNullTag to represent NULL pointer
      *this << (UInt_t) kNullTag;

   } else {

      // make sure fMap is initialized
      InitMap();

      ULong_t idx;
      UInt_t slot;
      ULong_t hash = Void_Hash(actualObjectStart);

      if ((idx = (ULong_t)fMap->GetValue(hash, (Long_t)actualObjectStart, slot)) != 0) {

         // truncation is OK the value we did put in the map is an 30-bit offset
         // and not a pointer
         UInt_t objIdx = UInt_t(idx);

         // save index of already stored object
         *this << objIdx;

      } else {

         // A warning to let the user know it will need to change the class code
         // to  be able to read this back.
         if (!actualClass->HasDefaultConstructor(kTRUE)) {
            Warning("WriteObjectAny", "since %s has no public constructor\n"
               "\twhich can be called without argument, objects of this class\n"
               "\tcan not be read with the current library. You will need to\n"
               "\tadd a default constructor before attempting to read it.",
               actualClass->GetName());
         }

         // reserve space for leading byte count
         UInt_t cntpos = UInt_t(fBufCur-fBuffer);
         fBufCur += sizeof(UInt_t);

         // write class of object first
         Int_t mapsize = fMap->Capacity(); // The slot depends on the capacity and WriteClass might induce an increase.
         WriteClass(actualClass);

         if (cacheReuse) {
            // add to map before writing rest of object (to handle self reference)
            // (+kMapOffset so it's != kNullTag)
            //MapObject(actualObjectStart, actualClass, cntpos+kMapOffset);
            UInt_t offset = cntpos+kMapOffset;
            if (mapsize == fMap->Capacity()) {
               fMap->AddAt(slot, hash, (Long_t)actualObjectStart, offset);
            } else {
               // The slot depends on the capacity and WriteClass has induced an increase.
               fMap->Add(hash, (Long_t)actualObjectStart, offset);
            }
            // No need to keep track of the class in write mode
            // fClassMap->Add(hash, (Long_t)obj, (Long_t)((TObject*)obj)->IsA());
            fMapCount++;
         }

         ((TClass*)actualClass)->Streamer((void*)actualObjectStart,*this);

         // write byte count
         SetByteCount(cntpos);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Read class definition from I/O buffer.
///
/// \param[in] clReq Can be used to cross check if the actually read object is of the requested class.
/// \param[in] objTag Set in case the object is a reference to an already read object.

TClass *TBufferFile::ReadClass(const TClass *clReq, UInt_t *objTag)
{
   R__ASSERT(IsReading());

   // read byte count and/or tag (older files don't have byte count)
   TClass *cl;
   if (fBufCur < fBuffer || fBufCur > fBufMax) {
      fBufCur = fBufMax;
      cl = (TClass*)-1;
      return cl;
   }
   UInt_t bcnt, tag, startpos = 0;
   *this >> bcnt;
   if (!(bcnt & kByteCountMask) || bcnt == kNewClassTag) {
      tag  = bcnt;
      bcnt = 0;
   } else {
      fVersion = 1;
      startpos = UInt_t(fBufCur-fBuffer);
      *this >> tag;
   }

   // in case tag is object tag return tag
   if (!(tag & kClassMask)) {
      if (objTag) *objTag = tag;
      return 0;
   }

   if (tag == kNewClassTag) {

      // got a new class description followed by a new object
      // (class can be 0 if class dictionary is not found, in that
      // case object of this class must be skipped)
      cl = TClass::Load(*this);

      // add class to fMap for later reference
      if (fVersion > 0) {
         // check if class was already read
         TClass *cl1 = (TClass *)(Long_t)fMap->GetValue(startpos+kMapOffset);
         if (cl1 != cl)
            MapObject(cl ? cl : (TObject*) -1, startpos+kMapOffset);
      } else
         MapObject(cl, fMapCount);

   } else {

      // got a tag to an already seen class
      UInt_t clTag = (tag & ~kClassMask);

      if (fVersion > 0) {
         clTag += fDisplacement;
         clTag = CheckObject(clTag, clReq, kTRUE);
      } else {
         if (clTag == 0 || clTag > (UInt_t)fMap->GetSize()) {
            Error("ReadClass", "illegal class tag=%d (0<tag<=%d), I/O buffer corrupted",
                  clTag, fMap->GetSize());
            // exception
         }
      }

      // class can be 0 if dictionary was not found
      cl = (TClass *)(Long_t)fMap->GetValue(clTag);
   }

   if (cl && clReq &&
       (!cl->InheritsFrom(clReq) &&
        !(clReq->GetSchemaRules() &&
          clReq->GetSchemaRules()->HasRuleWithSourceClass(cl->GetName()) )
        ) ) {
      Error("ReadClass", "The on-file class is \"'%s\" which is not compatible with the requested class: \"%s\"",
            cl->GetName(), clReq->GetName());
      // exception
   }

   // return bytecount in objTag
   if (objTag) *objTag = (bcnt & ~kByteCountMask);

   // case of unknown class
   if (!cl) cl = (TClass*)-1;

   return cl;
}

////////////////////////////////////////////////////////////////////////////////
/// Write class description to I/O buffer.

void TBufferFile::WriteClass(const TClass *cl)
{
   R__ASSERT(IsWriting());

   ULong_t idx;
   ULong_t hash = Void_Hash(cl);
   UInt_t slot;

   if ((idx = (ULong_t)fMap->GetValue(hash, (Long_t)cl,slot)) != 0) {

      // truncation is OK the value we did put in the map is an 30-bit offset
      // and not a pointer
      UInt_t clIdx = UInt_t(idx);

      // save index of already stored class
      *this << (clIdx | kClassMask);

   } else {

      // offset in buffer where class info is written
      UInt_t offset = UInt_t(fBufCur-fBuffer);

      // save new class tag
      *this << kNewClassTag;

      // write class name
      cl->Store(*this);

      // store new class reference in fMap (+kMapOffset so it's != kNullTag)
      CheckCount(offset+kMapOffset);
      fMap->AddAt(slot, hash, (Long_t)cl, offset+kMapOffset);
      fMapCount++;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Skip class version from I/O buffer.

void TBufferFile::SkipVersion(const TClass *cl)
{
   Version_t version;

   // not interested in byte count
   frombuf(this->fBufCur,&version);

   // if this is a byte count, then skip next short and read version
   if (version & kByteCountVMask) {
      frombuf(this->fBufCur,&version);
      frombuf(this->fBufCur,&version);
   }

   if (cl && cl->GetClassVersion() != 0  && version<=1) {
      if (version <= 0)  {
         UInt_t checksum = 0;
         //*this >> checksum;
         frombuf(this->fBufCur,&checksum);
         TStreamerInfo *vinfo = (TStreamerInfo*)cl->FindStreamerInfo(checksum);
         if (vinfo) {
            return;
         } else {
            // There are some cases (for example when the buffer was stored outside of
            // a ROOT file) where we do not have a TStreamerInfo.  If the checksum is
            // the one from the current class, we can still assume that we can read
            // the data so let use it.
            if (checksum==cl->GetCheckSum() || cl->MatchLegacyCheckSum(checksum)) {
               version = cl->GetClassVersion();
            } else {
               if (fParent) {
                  Error("SkipVersion", "Could not find the StreamerInfo with a checksum of %d for the class \"%s\" in %s.",
                        checksum, cl->GetName(), ((TFile*)fParent)->GetName());
               } else {
                  Error("SkipVersion", "Could not find the StreamerInfo with a checksum of %d for the class \"%s\" (buffer with no parent)",
                        checksum, cl->GetName());
               }
               return;
            }
         }
      }  else if (version == 1 && fParent && ((TFile*)fParent)->GetVersion()<40000 ) {
         // We could have a file created using a Foreign class before
         // the introduction of the CheckSum.  We need to check
         if ((!cl->IsLoaded() || cl->IsForeign()) &&
             Class_Has_StreamerInfo(cl) ) {

            const TList *list = ((TFile*)fParent)->GetStreamerInfoCache();
            const TStreamerInfo *local = list ? (TStreamerInfo*)list->FindObject(cl->GetName()) : 0;
            if ( local )  {
               UInt_t checksum = local->GetCheckSum();
               TStreamerInfo *vinfo = (TStreamerInfo*)cl->FindStreamerInfo(checksum);
               if (vinfo) {
                  version = vinfo->GetClassVersion();
               } else {
                  Error("SkipVersion", "Could not find the StreamerInfo with a checksum of %d for the class \"%s\" in %s.",
                        checksum, cl->GetName(), ((TFile*)fParent)->GetName());
                  return;
               }
            }
            else  {
               Error("SkipVersion", "Class %s not known to file %s.",
                     cl->GetName(), ((TFile*)fParent)->GetName());
               version = 0;
            }
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Read class version from I/O buffer.

Version_t TBufferFile::ReadVersion(UInt_t *startpos, UInt_t *bcnt, const TClass *cl)
{
   Version_t version;

   if (startpos) {
      // before reading object save start position
      *startpos = UInt_t(fBufCur-fBuffer);
   }

   // read byte count (older files don't have byte count)
   // byte count is packed in two individual shorts, this to be
   // backward compatible with old files that have at this location
   // only a single short (i.e. the version)
   union {
      UInt_t     cnt;
      Version_t  vers[2];
   } v;
#ifdef R__BYTESWAP
   frombuf(this->fBufCur,&v.vers[1]);
   frombuf(this->fBufCur,&v.vers[0]);
#else
   frombuf(this->fBufCur,&v.vers[0]);
   frombuf(this->fBufCur,&v.vers[1]);
#endif

   // no bytecount, backup and read version
   if (!(v.cnt & kByteCountMask)) {
      fBufCur -= sizeof(UInt_t);
      v.cnt = 0;
   }
   if (bcnt) *bcnt = (v.cnt & ~kByteCountMask);
   frombuf(this->fBufCur,&version);

   if (version<=1) {
      if (version <= 0)  {
         if (cl) {
            if (cl->GetClassVersion() != 0
                // If v.cnt < 6 then we have a class with a version that used to be zero and so there is no checksum.
                && (v.cnt && v.cnt >= 6)
                ) {
               UInt_t checksum = 0;
               //*this >> checksum;
               frombuf(this->fBufCur,&checksum);
               TStreamerInfo *vinfo = (TStreamerInfo*)cl->FindStreamerInfo(checksum);
               if (vinfo) {
                  return vinfo->TStreamerInfo::GetClassVersion(); // Try to get inlining.
               } else {
                  // There are some cases (for example when the buffer was stored outside of
                  // a ROOT file) where we do not have a TStreamerInfo.  If the checksum is
                  // the one from the current class, we can still assume that we can read
                  // the data so let use it.
                  if (checksum==cl->GetCheckSum() || cl->MatchLegacyCheckSum(checksum)) {
                     version = cl->GetClassVersion();
                  } else {
                     if (fParent) {
                        Error("ReadVersion", "Could not find the StreamerInfo with a checksum of 0x%x for the class \"%s\" in %s.",
                              checksum, cl->GetName(), ((TFile*)fParent)->GetName());
                     } else {
                        Error("ReadVersion", "Could not find the StreamerInfo with a checksum of 0x%x for the class \"%s\" (buffer with no parent)",
                              checksum, cl->GetName());
                     }
                     return 0;
                  }
               }
            }
         } else { // of if (cl) {
            UInt_t checksum = 0;
            //*this >> checksum;
            // If *bcnt < 6 then we have a class with 'just' version zero and no checksum
            if (v.cnt && v.cnt >= 6)
               frombuf(this->fBufCur,&checksum);
         }
      }  else if (version == 1 && fParent && ((TFile*)fParent)->GetVersion()<40000 && cl && cl->GetClassVersion() != 0) {
         // We could have a file created using a Foreign class before
         // the introduction of the CheckSum.  We need to check
         if ((!cl->IsLoaded() || cl->IsForeign()) &&
             Class_Has_StreamerInfo(cl) ) {

            const TList *list = ((TFile*)fParent)->GetStreamerInfoCache();
            const TStreamerInfo *local = list ? (TStreamerInfo*)list->FindObject(cl->GetName()) : 0;
            if ( local )  {
               UInt_t checksum = local->GetCheckSum();
               TStreamerInfo *vinfo = (TStreamerInfo*)cl->FindStreamerInfo(checksum);
               if (vinfo) {
                  version = vinfo->GetClassVersion();
               } else {
                  Error("ReadVersion", "Could not find the StreamerInfo with a checksum of 0x%x for the class \"%s\" in %s.",
                        checksum, cl->GetName(), ((TFile*)fParent)->GetName());
                  return 0;
               }
            }
            else  {
               Error("ReadVersion", "Class %s not known to file %s.",
                 cl->GetName(), ((TFile*)fParent)->GetName());
               version = 0;
            }
         }
      }
   }
   return version;
}

////////////////////////////////////////////////////////////////////////////////
/// Read class version from I/O buffer, when the caller knows for sure that
/// there is no checksum written/involved.

Version_t TBufferFile::ReadVersionNoCheckSum(UInt_t *startpos, UInt_t *bcnt)
{
   Version_t version;

   if (startpos) {
      // before reading object save start position
      *startpos = UInt_t(fBufCur-fBuffer);
   }

   // read byte count (older files don't have byte count)
   // byte count is packed in two individual shorts, this to be
   // backward compatible with old files that have at this location
   // only a single short (i.e. the version)
   union {
      UInt_t     cnt;
      Version_t  vers[2];
   } v;
#ifdef R__BYTESWAP
   frombuf(this->fBufCur,&v.vers[1]);
   frombuf(this->fBufCur,&v.vers[0]);
#else
   frombuf(this->fBufCur,&v.vers[0]);
   frombuf(this->fBufCur,&v.vers[1]);
#endif

   // no bytecount, backup and read version
   if (!(v.cnt & kByteCountMask)) {
      fBufCur -= sizeof(UInt_t);
      v.cnt = 0;
   }
   if (bcnt) *bcnt = (v.cnt & ~kByteCountMask);
   frombuf(this->fBufCur,&version);

   return version;
}

////////////////////////////////////////////////////////////////////////////////
/// Read class version from I/O buffer
///
/// To be used when streaming out member-wise streamed collection where we do not
/// care (not save) about the byte count and can safely ignore missing streamerInfo
/// (since they usually indicate empty collections).

Version_t TBufferFile::ReadVersionForMemberWise(const TClass *cl)
{
   Version_t version;

   // not interested in byte count
   frombuf(this->fBufCur,&version);

   if (version<=1) {
      if (version <= 0)  {
         if (cl) {
            if (cl->GetClassVersion() != 0) {
               UInt_t checksum = 0;
               frombuf(this->fBufCur,&checksum);
               TStreamerInfo *vinfo = (TStreamerInfo*)cl->FindStreamerInfo(checksum);
               if (vinfo) {
                  return vinfo->TStreamerInfo::GetClassVersion(); // Try to get inlining.
               } else {
                  // There are some cases (for example when the buffer was stored outside of
                  // a ROOT file) where we do not have a TStreamerInfo.  If the checksum is
                  // the one from the current class, we can still assume that we can read
                  // the data so let use it.
                  if (checksum==cl->GetCheckSum() || cl->MatchLegacyCheckSum(checksum)) {
                     version = cl->GetClassVersion();
                  } else {
                     // If we can not find the streamerInfo this means that
                     // we do not actually need it (the collection is always empty
                     // in this file), so no need to issue a warning.
                     return 0;
                  }
               }
            }
         } else { // of if (cl) {
            UInt_t checksum = 0;
            frombuf(this->fBufCur,&checksum);
         }
      }  else if (version == 1 && fParent && ((TFile*)fParent)->GetVersion()<40000 && cl && cl->GetClassVersion() != 0) {
         // We could have a file created using a Foreign class before
         // the introduction of the CheckSum.  We need to check
         if ((!cl->IsLoaded() || cl->IsForeign()) && Class_Has_StreamerInfo(cl) ) {

            const TList *list = ((TFile*)fParent)->GetStreamerInfoCache();
            const TStreamerInfo *local = list ? (TStreamerInfo*)list->FindObject(cl->GetName()) : 0;
            if ( local )  {
               UInt_t checksum = local->GetCheckSum();
               TStreamerInfo *vinfo = (TStreamerInfo*)cl->FindStreamerInfo(checksum);
               if (vinfo) {
                  version = vinfo->GetClassVersion();
               } else {
                  // If we can not find the streamerInfo this means that
                  // we do not actually need it (the collection is always empty
                  // in this file), so no need to issue a warning.
                  return 0;
               }
            }
            else  {
               Error("ReadVersion", "Class %s not known to file %s.",
                 cl->GetName(), ((TFile*)fParent)->GetName());
               version = 0;
            }
         }
      }
   }
   return version;
}

////////////////////////////////////////////////////////////////////////////////
/// Write class version to I/O buffer.

UInt_t TBufferFile::WriteVersion(const TClass *cl, Bool_t useBcnt)
{
   UInt_t cntpos = 0;
   if (useBcnt) {
      // reserve space for leading byte count
      cntpos   = UInt_t(fBufCur-fBuffer);
      fBufCur += sizeof(UInt_t);
   }

   Version_t version = cl->GetClassVersion();
   if (version<=1 && cl->IsForeign()) {
      *this << Version_t(0);
      *this << cl->GetCheckSum();
   } else {
      if (version > kMaxVersion) {
         Error("WriteVersion", "version number cannot be larger than %hd)",
               kMaxVersion);
         version = kMaxVersion;
      }
      *this <<version;
   }

   // return position where to store possible byte count
   return cntpos;
}

////////////////////////////////////////////////////////////////////////////////
/// Write class version to I/O buffer after setting the kStreamedMemberWise
/// bit in the version number.

UInt_t TBufferFile::WriteVersionMemberWise(const TClass *cl, Bool_t useBcnt)
{
   UInt_t cntpos = 0;
   if (useBcnt) {
      // reserve space for leading byte count
      cntpos   = UInt_t(fBufCur-fBuffer);
      fBufCur += sizeof(UInt_t);
   }

   Version_t version = cl->GetClassVersion();
   if (version<=1 && cl->IsForeign()) {
      Error("WriteVersionMemberWise", "Member-wise streaming of foreign collection not yet implemented!");
      *this << Version_t(0);
      *this << cl->GetCheckSum();
   } else {
      if (version > kMaxVersion) {
         Error("WriteVersionMemberWise", "version number cannot be larger than %hd)",
               kMaxVersion);
         version = kMaxVersion;
      }
      version |= kStreamedMemberWise;
      *this <<version;
   }

   // return position where to store possible byte count
   return cntpos;
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object given its C++ typeinfo information.

void TBufferFile::StreamObject(void *obj, const std::type_info &typeinfo, const TClass* onFileClass )
{
   TClass *cl = TClass::GetClass(typeinfo);
   if (cl) cl->Streamer(obj, *this, (TClass*)onFileClass );
   else Warning("StreamObject","No TClass for the type %s is available, the object was not read.", typeinfo.name());
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object given the name of its actual class.

void TBufferFile::StreamObject(void *obj, const char *className, const TClass* onFileClass)
{
   TClass *cl = TClass::GetClass(className);
   if (cl) cl->Streamer(obj, *this, (TClass*)onFileClass );
   else Warning("StreamObject","No TClass for the type %s is available, the object was not read.", className);
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object given a pointer to its actual class.

void TBufferFile::StreamObject(void *obj, const TClass *cl, const TClass* onFileClass )
{
   ((TClass*)cl)->Streamer(obj, *this, (TClass*)onFileClass );
}

////////////////////////////////////////////////////////////////////////////////
/// Stream an object inheriting from TObject using its streamer.

void TBufferFile::StreamObject(TObject *obj)
{
   obj->Streamer(*this);
}

////////////////////////////////////////////////////////////////////////////////
/// Check if offset is not too large (< kMaxMapCount) when writing.

void TBufferFile::CheckCount(UInt_t offset)
{
   if (IsWriting()) {
      if (offset >= kMaxMapCount) {
         Error("CheckCount", "buffer offset too large (larger than %d)", kMaxMapCount);
         // exception
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Check for object in the read map. If the object is 0 it still has to be
/// read. Try to read it from the buffer starting at location offset. If the
/// object is -1 then it really does not exist and we return 0. If the object
/// exists just return the offset.

UInt_t TBufferFile::CheckObject(UInt_t offset, const TClass *cl, Bool_t readClass)
{
   // in position 0 we always have the reference to the null object
   if (!offset) return offset;

   Long_t cli;

   if (readClass) {
      if ((cli = fMap->GetValue(offset)) == 0) {
         // No class found at this location in map. It might have been skipped
         // as part of a skipped object. Try to explicitly read the class.

         // save fBufCur and set to place specified by offset (-kMapOffset-sizeof(bytecount))
         char *bufsav = fBufCur;
         fBufCur = (char *)(fBuffer + offset-kMapOffset-sizeof(UInt_t));

         TClass *c = ReadClass(cl);
         if (c == (TClass*) -1) {
            // mark class as really not available
            fMap->Remove(offset);
            fMap->Add(offset, -1);
            offset = 0;
            if (cl)
               Warning("CheckObject", "reference to unavailable class %s,"
                       " pointers of this type will be 0", cl->GetName());
            else
               Warning("CheckObject", "reference to an unavailable class,"
                       " pointers of that type will be 0");
         }

         fBufCur = bufsav;

      } else if (cli == -1) {

         // class really does not exist
         return 0;
      }

   } else {

      if ((cli = fMap->GetValue(offset)) == 0) {
         // No object found at this location in map. It might have been skipped
         // as part of a skipped object. Try to explicitly read the object.

         // save fBufCur and set to place specified by offset (-kMapOffset)
         char *bufsav = fBufCur;
         fBufCur = (char *)(fBuffer + offset-kMapOffset);

         TObject *obj = ReadObject(cl);
         if (!obj) {
            // mark object as really not available
            fMap->Remove(offset);
            fMap->Add(offset, -1);
            Warning("CheckObject", "reference to object of unavailable class %s, offset=%d"
                    " pointer will be 0", cl ? cl->GetName() : "TObject",offset);
            offset = 0;
         }

         fBufCur = bufsav;

      } else if (cli == -1) {

         // object really does not exist
         return 0;
      }

   }

   return offset;
}


////////////////////////////////////////////////////////////////////////////////
/// Read max bytes from the I/O buffer into buf. The function returns
/// the actual number of bytes read.

Int_t TBufferFile::ReadBuf(void *buf, Int_t max)
{
   R__ASSERT(IsReading());

   if (max == 0) return 0;

   Int_t n = TMath::Min(max, (Int_t)(fBufMax - fBufCur));

   memcpy(buf, fBufCur, n);
   fBufCur += n;

   return n;
}

////////////////////////////////////////////////////////////////////////////////
/// Write max bytes from buf into the I/O buffer.

void TBufferFile::WriteBuf(const void *buf, Int_t max)
{
   R__ASSERT(IsWriting());

   if (max == 0) return;

   if (fBufCur + max > fBufMax) AutoExpand(fBufSize+max); // a more precise request would be: fBufSize + max - (fBufMax - fBufCur)

   memcpy(fBufCur, buf, max);
   fBufCur += max;
}

////////////////////////////////////////////////////////////////////////////////
/// Read string from I/O buffer. String is read till 0 character is
/// found or till max-1 characters are read (i.e. string s has max
/// bytes allocated). If max = -1 no check on number of character is
/// made, reading continues till 0 character is found.

char *TBufferFile::ReadString(char *s, Int_t max)
{
   R__ASSERT(IsReading());

   char  ch;
   Int_t nr = 0;

   if (max == -1) max = kMaxInt;

   while (nr < max-1) {

      *this >> ch;

      // stop when 0 read
      if (ch == 0) break;

      s[nr++] = ch;
   }

   s[nr] = 0;
   return s;
}

////////////////////////////////////////////////////////////////////////////////
/// Write string to I/O buffer. Writes string upto and including the
/// terminating 0.

void TBufferFile::WriteString(const char *s)
{
   WriteBuf(s, (strlen(s)+1)*sizeof(char));
}

////////////////////////////////////////////////////////////////////////////////
/// Read emulated class.

Int_t TBufferFile::ReadClassEmulated(const TClass *cl, void *object, const TClass *onFileClass)
{
   UInt_t start,count;
   //We assume that the class was written with a standard streamer
   //We attempt to recover if a version count was not written
   Version_t v = ReadVersion(&start,&count);

   if (count) {
      TStreamerInfo *sinfo = nullptr;
      if( onFileClass ) {
         sinfo = (TStreamerInfo*)cl->GetConversionStreamerInfo( onFileClass, v );
         if( !sinfo )
            return 0;
      }

      sinfo = (TStreamerInfo*)cl->GetStreamerInfo(v);
      ApplySequence(*(sinfo->GetReadObjectWiseActions()), object);
      if (sinfo->IsRecovered()) count=0;
      CheckByteCount(start,count,cl);
   } else {
      SetBufferOffset(start);
      TStreamerInfo *sinfo = ((TStreamerInfo*)cl->GetStreamerInfo());
      ApplySequence(*(sinfo->GetReadObjectWiseActions()), object);
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Deserialize information from a buffer into an object.
///
/// Note: This function is called by the xxx::Streamer() functions in
/// rootcint-generated dictionaries.
/// This function assumes that the class version and the byte count
/// information have been read.
///
/// \param[in] version The version number of the class
/// \param[in] start   The starting position in the buffer b
/// \param[in] count   The number of bytes for this object in the buffer
///

Int_t TBufferFile::ReadClassBuffer(const TClass *cl, void *pointer, Int_t version, UInt_t start, UInt_t count, const TClass *onFileClass)
{

   //---------------------------------------------------------------------------
   // The ondisk class has been specified so get foreign streamer info
   /////////////////////////////////////////////////////////////////////////////

   TStreamerInfo *sinfo = nullptr;
   if( onFileClass ) {
      sinfo = (TStreamerInfo*)cl->GetConversionStreamerInfo( onFileClass, version );
      if( !sinfo ) {
         Error("ReadClassBuffer",
               "Could not find the right streamer info to convert %s version %d into a %s, object skipped at offset %d",
               onFileClass->GetName(), version, cl->GetName(), Length() );
         CheckByteCount(start, count, onFileClass);
         return 0;
      }
   }
   //---------------------------------------------------------------------------
   // Get local streamer info
   /////////////////////////////////////////////////////////////////////////////
   /// The StreamerInfo should exist at this point.

   else {
      R__READ_LOCKGUARD(ROOT::gCoreMutex);
      auto infos = cl->GetStreamerInfos();
      auto ninfos = infos->GetSize();
      if (version < -1 || version >= ninfos) {
         Error("ReadClassBuffer", "class: %s, attempting to access a wrong version: %d, object skipped at offset %d",
               cl->GetName(), version, Length() );
         CheckByteCount(start, count, cl);
         return 0;
      }
      sinfo = (TStreamerInfo*)infos->At(version);
      if (sinfo == nullptr) {
         // Unless the data is coming via a socket connection from with schema evolution
         // (tracking) was not enabled.  So let's create the StreamerInfo if it is the
         // one for the current version, otherwise let's complain ...
         // We could also get here if there old class version was '1' and the new class version is higher than 1
         // AND the checksum is the same.
         R__WRITE_LOCKGUARD(ROOT::gCoreMutex);
         // check if another thread took care of this already
         sinfo = (TStreamerInfo*)cl->GetStreamerInfos()->At(version);
         if (sinfo == nullptr) {
            if ( version == cl->GetClassVersion() || version == 1 ) {
               const_cast<TClass*>(cl)->BuildRealData(pointer);
               // This creation is alright since we just checked within the
               // current 'locked' section.
               sinfo = new TStreamerInfo(const_cast<TClass*>(cl));
               const_cast<TClass*>(cl)->RegisterStreamerInfo(sinfo);
               if (gDebug > 0) Info("ReadClassBuffer", "Creating StreamerInfo for class: %s, version: %d", cl->GetName(), version);
               sinfo->Build();
            } else if (version==0) {
               // When the object was written the class was version zero, so
               // there is no StreamerInfo to be found.
               // Check that the buffer position corresponds to the byte count.
               CheckByteCount(start, count, cl);
               return 0;
            } else {
               Error("ReadClassBuffer", "Could not find the StreamerInfo for version %d of the class %s, object skipped at offset %d",
                     version, cl->GetName(), Length() );
               CheckByteCount(start, count, cl);
               return 0;
            }
         }
      } else if (!sinfo->IsCompiled()) {  // Note this read is protected by the above lock.
         // Streamer info has not been compiled, but exists.
         // Therefore it was read in from a file and we have to do schema evolution.
         R__WRITE_LOCKGUARD(ROOT::gCoreMutex);
         // check if another thread took care of this already
         if (!sinfo->IsCompiled()) {
            const_cast<TClass*>(cl)->BuildRealData(pointer);
            sinfo->BuildOld();
         }
      }
   }

   // Deserialize the object.
   ApplySequence(*(sinfo->GetReadObjectWiseActions()), (char*)pointer);
   if (sinfo->IsRecovered()) count=0;

   // Check that the buffer position corresponds to the byte count.
   CheckByteCount(start, count, cl);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Deserialize information from a buffer into an object.
///
/// Note: This function is called by the xxx::Streamer()
/// functions in rootcint-generated dictionaries.
///

Int_t TBufferFile::ReadClassBuffer(const TClass *cl, void *pointer, const TClass *onFileClass)
{
   // Read the class version from the buffer.
   UInt_t R__s = 0; // Start of object.
   UInt_t R__c = 0; // Count of bytes.
   Version_t version;

   if( onFileClass )
      version = ReadVersion(&R__s, &R__c, onFileClass);
   else
      version = ReadVersion(&R__s, &R__c, cl);

   Bool_t v2file = kFALSE;
   TFile *file = (TFile*)GetParent();
   if (file && file->GetVersion() < 30000) {
      version = -1; //This is old file
      v2file = kTRUE;
   }

   //---------------------------------------------------------------------------
   // The ondisk class has been specified so get foreign streamer info
   /////////////////////////////////////////////////////////////////////////////

   TStreamerInfo *sinfo = nullptr;
   if( onFileClass ) {
      sinfo = (TStreamerInfo*)cl->GetConversionStreamerInfo( onFileClass, version );
      if( !sinfo ) {
         Error("ReadClassBuffer",
               "Could not find the right streamer info to convert %s version %d into a %s, object skipped at offset %d",
               onFileClass->GetName(), version, cl->GetName(), Length() );
         CheckByteCount(R__s, R__c, onFileClass);
         return 0;
      }
   }
   //---------------------------------------------------------------------------
   // Get local streamer info
   /////////////////////////////////////////////////////////////////////////////
   /// The StreamerInfo should exist at this point.

   else {
      TStreamerInfo *guess = (TStreamerInfo*)cl->GetLastReadInfo();
      if (guess && guess->GetClassVersion() == version) {
         sinfo = guess;
      } else {
         // The last one is not the one we are looking for.
         {
            R__LOCKGUARD(gInterpreterMutex);

            const TObjArray *infos = cl->GetStreamerInfos();
            Int_t infocapacity = infos->Capacity();
            if (infocapacity) {
               if (version < -1 || version >= infocapacity) {
                  Error("ReadClassBuffer","class: %s, attempting to access a wrong version: %d, object skipped at offset %d",
                        cl->GetName(), version, Length());
                  CheckByteCount(R__s, R__c, cl);
                  return 0;
               }
               sinfo = (TStreamerInfo*) infos->UncheckedAt(version);
               if (sinfo) {
                  if (!sinfo->IsCompiled())
                  {
                     // Streamer info has not been compiled, but exists.
                     // Therefore it was read in from a file and we have to do schema evolution?
                     R__LOCKGUARD(gInterpreterMutex);
                     const_cast<TClass*>(cl)->BuildRealData(pointer);
                     sinfo->BuildOld();
                  }
                  // If the compilation succeeded, remember this StreamerInfo.
                  // const_cast okay because of the lock on gInterpreterMutex.
                  if (sinfo->IsCompiled()) const_cast<TClass*>(cl)->SetLastReadInfo(sinfo);
               }
            }
         }

         if (sinfo == nullptr) {
            // Unless the data is coming via a socket connection from with schema evolution
            // (tracking) was not enabled.  So let's create the StreamerInfo if it is the
            // one for the current version, otherwise let's complain ...
            // We could also get here when reading a file prior to the introduction of StreamerInfo.
            // We could also get here if there old class version was '1' and the new class version is higher than 1
            // AND the checksum is the same.
            if (v2file || version == cl->GetClassVersion() || version == 1 ) {
               R__LOCKGUARD(gInterpreterMutex);

               // We need to check if another thread did not get here first
               // and did the StreamerInfo creation already.
               auto infos = cl->GetStreamerInfos();
               auto ninfos = infos->GetSize();
               if (!(version < -1 || version >= ninfos)) {
                  sinfo = (TStreamerInfo *) infos->At(version);
               }
               if (!sinfo) {
                  const_cast<TClass *>(cl)->BuildRealData(pointer);
                  sinfo = new TStreamerInfo(const_cast<TClass *>(cl));
                  sinfo->SetClassVersion(version);
                  const_cast<TClass *>(cl)->RegisterStreamerInfo(sinfo);
                  if (gDebug > 0)
                     Info("ReadClassBuffer", "Creating StreamerInfo for class: %s, version: %d",
                           cl->GetName(), version);
                  if (v2file) {
                     sinfo->Build();             // Get the elements.
                     sinfo->Clear("build");      // Undo compilation.
                     sinfo->BuildEmulated(file); // Fix the types and redo compilation.
                  } else {
                     sinfo->Build();
                  }
               }
            } else if (version==0) {
               // When the object was written the class was version zero, so
               // there is no StreamerInfo to be found.
               // Check that the buffer position corresponds to the byte count.
               CheckByteCount(R__s, R__c, cl);
               return 0;
            } else {
               Error( "ReadClassBuffer", "Could not find the StreamerInfo for version %d of the class %s, object skipped at offset %d",
                     version, cl->GetName(), Length() );
               CheckByteCount(R__s, R__c, cl);
               return 0;
            }
         }
      }
   }

   //deserialize the object
   ApplySequence(*(sinfo->GetReadObjectWiseActions()), (char*)pointer );
   if (sinfo->TStreamerInfo::IsRecovered()) R__c=0; // 'TStreamerInfo::' avoids going via a virtual function.

   // Check that the buffer position corresponds to the byte count.
   CheckByteCount(R__s, R__c, cl);

   if (gDebug > 2) Info("ReadClassBuffer", "For class: %s has read %d bytes", cl->GetName(), R__c);

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Function called by the Streamer functions to serialize object at p
/// to buffer b. The optional argument info may be specified to give an
/// alternative StreamerInfo instead of using the default StreamerInfo
/// automatically built from the class definition.
/// For more information, see class TStreamerInfo.

Int_t TBufferFile::WriteClassBuffer(const TClass *cl, void *pointer)
{
   //build the StreamerInfo if first time for the class
   TStreamerInfo *sinfo = (TStreamerInfo*)const_cast<TClass*>(cl)->GetCurrentStreamerInfo();
   if (sinfo == nullptr) {
      //Have to be sure between the check and the taking of the lock if the current streamer has changed
      R__LOCKGUARD(gInterpreterMutex);
      sinfo = (TStreamerInfo*)const_cast<TClass*>(cl)->GetCurrentStreamerInfo();
      if (sinfo == nullptr) {
         const_cast<TClass*>(cl)->BuildRealData(pointer);
         sinfo = new TStreamerInfo(const_cast<TClass*>(cl));
         const_cast<TClass*>(cl)->SetCurrentStreamerInfo(sinfo);
         const_cast<TClass*>(cl)->RegisterStreamerInfo(sinfo);
         if (gDebug > 0) Info("WritedClassBuffer", "Creating StreamerInfo for class: %s, version: %d",cl->GetName(),cl->GetClassVersion());
         sinfo->Build();
      }
   } else if (!sinfo->IsCompiled()) {
      R__LOCKGUARD(gInterpreterMutex);
      // Redo the test in case we have been victim of a data race on fIsCompiled.
      if (!sinfo->IsCompiled()) {
         const_cast<TClass*>(cl)->BuildRealData(pointer);
         sinfo->BuildOld();
      }
   }

   //write the class version number and reserve space for the byte count
   UInt_t R__c = WriteVersion(cl, kTRUE);

   //NOTE: In the future Philippe wants this to happen via a custom action
   TagStreamerInfo(sinfo);
   ApplySequence(*(sinfo->GetWriteObjectWiseActions()), (char*)pointer);

   //write the byte count at the start of the buffer
   SetByteCount(R__c, kTRUE);

   if (gDebug > 2) Info("WritedClassBuffer", "For class: %s version %d has written %d bytes",cl->GetName(),cl->GetClassVersion(),UInt_t(fBufCur - fBuffer) - R__c - (UInt_t)sizeof(UInt_t));
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Read one collection of objects from the buffer using the StreamerInfoLoopAction.
/// The collection needs to be a split TClonesArray or a split vector of pointers.

Int_t TBufferFile::ApplySequence(const TStreamerInfoActions::TActionSequence &sequence, void *obj)
{
   if (gDebug) {
      //loop on all active members
      TStreamerInfoActions::ActionContainer_t::const_iterator end = sequence.fActions.end();
      for(TStreamerInfoActions::ActionContainer_t::const_iterator iter = sequence.fActions.begin();
          iter != end;
          ++iter) {
         (*iter).PrintDebug(*this,obj);
         (*iter)(*this,obj);
      }

   } else {
      //loop on all active members
      TStreamerInfoActions::ActionContainer_t::const_iterator end = sequence.fActions.end();
      for(TStreamerInfoActions::ActionContainer_t::const_iterator iter = sequence.fActions.begin();
          iter != end;
          ++iter) {
         (*iter)(*this,obj);
      }
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Read one collection of objects from the buffer using the StreamerInfoLoopAction.
/// The collection needs to be a split TClonesArray or a split vector of pointers.

Int_t TBufferFile::ApplySequenceVecPtr(const TStreamerInfoActions::TActionSequence &sequence, void *start_collection, void *end_collection)
{
   if (gDebug) {
      //loop on all active members
      TStreamerInfoActions::ActionContainer_t::const_iterator end = sequence.fActions.end();
      for(TStreamerInfoActions::ActionContainer_t::const_iterator iter = sequence.fActions.begin();
          iter != end;
          ++iter) {
         if (!start_collection || start_collection == end_collection)
            (*iter).PrintDebug(*this, nullptr);  // Warning: This limits us to TClonesArray and vector of pointers.
         else
            (*iter).PrintDebug(*this, *(char**)start_collection);  // Warning: This limits us to TClonesArray and vector of pointers.
         (*iter)(*this, start_collection, end_collection);
      }

   } else {
      //loop on all active members
      TStreamerInfoActions::ActionContainer_t::const_iterator end = sequence.fActions.end();
      for(TStreamerInfoActions::ActionContainer_t::const_iterator iter = sequence.fActions.begin();
          iter != end;
          ++iter) {
         (*iter)(*this,start_collection,end_collection);
      }
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Read one collection of objects from the buffer using the StreamerInfoLoopAction.

Int_t TBufferFile::ApplySequence(const TStreamerInfoActions::TActionSequence &sequence, void *start_collection, void *end_collection)
{
   TStreamerInfoActions::TLoopConfiguration *loopconfig = sequence.fLoopConfig;
   if (gDebug) {

      // Get the address of the first item for the PrintDebug.
      // (Performance is not essential here since we are going to print to
      // the screen anyway).
      void *arr0 = start_collection ? loopconfig->GetFirstAddress(start_collection,end_collection) : 0;
      // loop on all active members
      TStreamerInfoActions::ActionContainer_t::const_iterator end = sequence.fActions.end();
      for(TStreamerInfoActions::ActionContainer_t::const_iterator iter = sequence.fActions.begin();
          iter != end;
          ++iter) {
         (*iter).PrintDebug(*this,arr0);
         (*iter)(*this,start_collection,end_collection,loopconfig);
      }

   } else {
      //loop on all active members
      TStreamerInfoActions::ActionContainer_t::const_iterator end = sequence.fActions.end();
      for(TStreamerInfoActions::ActionContainer_t::const_iterator iter = sequence.fActions.begin();
          iter != end;
          ++iter) {
         (*iter)(*this,start_collection,end_collection,loopconfig);
      }
   }

   return 0;
}
