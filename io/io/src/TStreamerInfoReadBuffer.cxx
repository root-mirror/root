// @(#)root/io:$Id$
// Author: Rene Brun   12/10/2000

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TBuffer.h"
#include "TFile.h"
#include "TClass.h"
#include "TBufferFile.h"
#include "TClonesArray.h"
#include "TError.h"
#include "TRef.h"
#include "TProcessID.h"
#include "TStreamer.h"
#include "TStreamerElement.h"
#include "TStreamerInfo.h"
#include "TVirtualCollectionProxy.h"
#include "TContainerConverters.h"
#include "TVirtualArray.h"
#include "TVirtualObject.h"
#include "ThreadLocalStorage.h"

// GetCurrentElement.
// Currently only used by TRef::Streamer.

static TStreamerElement* &CurrentElement()
{
   //Pointer to current TStreamerElement
   //Thread local storage.

   TTHREAD_TLS(TStreamerElement*) fgElement(0);

   return fgElement;
}

////////////////////////////////////////////////////////////////////////////////
///static function returning a pointer to the current TStreamerElement
///fgElement points to the current TStreamerElement being read in ReadBuffer

TStreamerElement *TStreamerInfo::GetCurrentElement()
{
   return CurrentElement();
}

//==========CPP macros

#define DOLOOP for(Int_t k=0; k<narr; ++k)

#define ReadBasicTypeElem(name,index)           \
   {                                            \
      name *x=(name*)(arr[index]+ioffset);      \
      b >> *x;                                  \
   }

#define ReadBasicType(name)                     \
   {                                            \
      ReadBasicTypeElem(name,0);                \
   }

#define ReadBasicTypeLoop(name)                                \
   {                                                           \
      for(Int_t k=0; k<narr; ++k) ReadBasicTypeElem(name,k);   \
   }

#define ReadBasicArrayElem(name,index)          \
   {                                            \
      name *x=(name*)(arr[index]+ioffset);      \
      b.ReadFastArray(x,compinfo[i]->fLength);            \
   }

#define ReadBasicArray(name)                     \
   {                                             \
      ReadBasicArrayElem(name,0);                \
   }

#define ReadBasicArrayLoop(name)                               \
   {                                                           \
      for(Int_t k=0; k<narr; ++k) ReadBasicArrayElem(name,k)   \
   }

#define ReadBasicPointerElem(name,index)        \
   {                                            \
      Char_t isArray;                           \
      b >> isArray;                             \
      Int_t *l = (Int_t*)(arr[index]+imethod);  \
      if (*l < 0 || *l > b.BufferSize()) continue;     \
      name **f = (name**)(arr[index]+ioffset);  \
      int j;                                    \
      if (isArray) for(j=0;j<compinfo[i]->fLength;j++) {  \
         delete [] f[j];                        \
         f[j] = 0; if (*l <=0) continue;        \
         f[j] = new name[*l];                   \
         b.ReadFastArray(f[j],*l);              \
      }                                         \
      else for(j=0;j<compinfo[i]->fLength;j++) {  \
        delete [] f[j];                        \
        f[j] = 0;                              \
     }                                         \
   }

#define ReadBasicPointer(name)                  \
   {                                            \
      const int imethod = compinfo[i]->fMethod+eoffset;   \
      ReadBasicPointerElem(name,0);             \
   }

#define ReadBasicPointerLoop(name)              \
   {                                            \
      int imethod = compinfo[i]->fMethod+eoffset;         \
      for(int k=0; k<narr; ++k) {               \
         ReadBasicPointerElem(name,k);          \
      }                                         \
   }

#define SkipCBasicType(name)                    \
   {                                            \
      name dummy;                               \
      DOLOOP{ b >> dummy; }                     \
      break;                                    \
   }

#define SkipCFloat16(name)                       \
   {                                             \
      name dummy;                                \
      DOLOOP { b.ReadFloat16(&dummy,aElement); } \
      break;                                     \
   }

#define SkipCDouble32(name)                      \
   {                                             \
      name dummy;                                \
      DOLOOP { b.ReadDouble32(&dummy,aElement); }\
      break;                                     \
   }

#define SkipCBasicArray(name,ReadArrayFunc)              \
    {                                                    \
      name* readbuf = new name[compinfo->fLength];              \
      DOLOOP {                                           \
          b.ReadArrayFunc(readbuf, compinfo->fLength);          \
      }                                                  \
      delete[] readbuf;                                  \
      break;                                             \
    }

#define SkipCBasicPointer(name,ReadArrayFunc)                             \
   {                                                                      \
      Int_t addCounter = -111;                                            \
      if ((imethod>0) && (compinfo->fMethod>0)) addCounter = -1;                 \
      if((addCounter<-1) && (aElement!=0) && (aElement->IsA()==TStreamerBasicPointer::Class())) { \
         TStreamerElement* elemCounter = (TStreamerElement*) thisVar->GetElements()->FindObject(((TStreamerBasicPointer*)aElement)->GetCountName()); \
         if (elemCounter) addCounter = elemCounter->GetTObjectOffset();   \
      }                                                                   \
      if (addCounter>=-1) {                                               \
         int len = aElement->GetArrayDim()?aElement->GetArrayLength():1;  \
         Char_t isArray;                                                  \
         DOLOOP {                                                         \
            b >> isArray;                                                 \
            char *arr_k = arr[k];                                         \
            Int_t *l = (addCounter==-1 && arr_k) ? (Int_t*)(arr_k+imethod) : &addCounter;  \
            if (*l>0) {                                                   \
               name* readbuf = new name[*l];                              \
               for (int j=0;j<len;j++)                                    \
                  b.ReadArrayFunc(readbuf, *l);                           \
               delete[] readbuf;                                          \
            }                                                             \
         }                                                                \
      }                                                                   \
      break;                                                              \
   }

////////////////////////////////////////////////////////////////////////////////
/// Skip an element.

template <class T>
Int_t TStreamerInfo::ReadBufferSkip(TBuffer &b, const T &arr, const TCompInfo *compinfo, Int_t kase,
                                    TStreamerElement *aElement, Int_t narr,
                                    Int_t eoffset)
{
   TStreamerInfo* thisVar = this;

   //  Skip elements in a TClonesArray

   TClass* cle = compinfo->fClass;

   Int_t imethod = compinfo->fMethod+eoffset;

   switch (kase) {

      // skip basic types
      case TStreamerInfo::kSkip + TStreamerInfo::kBool:      SkipCBasicType(Bool_t);
      case TStreamerInfo::kSkip + TStreamerInfo::kChar:      SkipCBasicType(Char_t);
      case TStreamerInfo::kSkip + TStreamerInfo::kShort:     SkipCBasicType(Short_t);
      case TStreamerInfo::kSkip + TStreamerInfo::kInt:       SkipCBasicType(Int_t);
      case TStreamerInfo::kSkip + TStreamerInfo::kLong:      SkipCBasicType(Long_t);
      case TStreamerInfo::kSkip + TStreamerInfo::kLong64:    SkipCBasicType(Long64_t);
      case TStreamerInfo::kSkip + TStreamerInfo::kFloat:     SkipCBasicType(Float_t);
      case TStreamerInfo::kSkip + TStreamerInfo::kFloat16:   SkipCFloat16(Float_t);
      case TStreamerInfo::kSkip + TStreamerInfo::kDouble:    SkipCBasicType(Double_t);
      case TStreamerInfo::kSkip + TStreamerInfo::kDouble32:  SkipCDouble32(Double32_t)
      case TStreamerInfo::kSkip + TStreamerInfo::kUChar:     SkipCBasicType(UChar_t);
      case TStreamerInfo::kSkip + TStreamerInfo::kUShort:    SkipCBasicType(UShort_t);
      case TStreamerInfo::kSkip + TStreamerInfo::kUInt:      SkipCBasicType(UInt_t);
      case TStreamerInfo::kSkip + TStreamerInfo::kULong:     SkipCBasicType(ULong_t);
      case TStreamerInfo::kSkip + TStreamerInfo::kULong64:   SkipCBasicType(ULong64_t);
      case TStreamerInfo::kSkip + TStreamerInfo::kBits:      {
         UInt_t dummy;
         DOLOOP{
            b >> dummy;
            if ((dummy & kIsReferenced) != 0) {
               UShort_t pidf;
               b >> pidf;
            }
         }
         break;
      }

         // skip array of basic types  array[8]
      case TStreamerInfo::kSkipL + TStreamerInfo::kBool:     SkipCBasicArray(Bool_t,ReadFastArray);
      case TStreamerInfo::kSkipL + TStreamerInfo::kChar:     SkipCBasicArray(Char_t,ReadFastArray);
      case TStreamerInfo::kSkipL + TStreamerInfo::kShort:    SkipCBasicArray(Short_t,ReadFastArray);
      case TStreamerInfo::kSkipL + TStreamerInfo::kInt:      SkipCBasicArray(Int_t,ReadFastArray);
      case TStreamerInfo::kSkipL + TStreamerInfo::kLong:     SkipCBasicArray(Long_t,ReadFastArray);
      case TStreamerInfo::kSkipL + TStreamerInfo::kLong64:   SkipCBasicArray(Long64_t,ReadFastArray);
      case TStreamerInfo::kSkipL + TStreamerInfo::kFloat16:  SkipCBasicArray(Float_t,ReadFastArrayFloat16);
      case TStreamerInfo::kSkipL + TStreamerInfo::kFloat:    SkipCBasicArray(Float_t,ReadFastArray);
      case TStreamerInfo::kSkipL + TStreamerInfo::kDouble32: SkipCBasicArray(Double_t,ReadFastArrayDouble32)
      case TStreamerInfo::kSkipL + TStreamerInfo::kDouble:   SkipCBasicArray(Double_t,ReadFastArray);
      case TStreamerInfo::kSkipL + TStreamerInfo::kUChar:    SkipCBasicArray(UChar_t,ReadFastArray);
      case TStreamerInfo::kSkipL + TStreamerInfo::kUShort:   SkipCBasicArray(UShort_t,ReadFastArray);
      case TStreamerInfo::kSkipL + TStreamerInfo::kUInt:     SkipCBasicArray(UInt_t,ReadFastArray);
      case TStreamerInfo::kSkipL + TStreamerInfo::kULong:    SkipCBasicArray(ULong_t,ReadFastArray);
      case TStreamerInfo::kSkipL + TStreamerInfo::kULong64:  SkipCBasicArray(ULong64_t,ReadFastArray);

   // skip pointer to an array of basic types  array[n]
      case TStreamerInfo::kSkipP + TStreamerInfo::kBool:     SkipCBasicPointer(Bool_t,ReadFastArray);
      case TStreamerInfo::kSkipP + TStreamerInfo::kChar:     SkipCBasicPointer(Char_t,ReadFastArray);
      case TStreamerInfo::kSkipP + TStreamerInfo::kShort:    SkipCBasicPointer(Short_t,ReadFastArray);
      case TStreamerInfo::kSkipP + TStreamerInfo::kInt:      SkipCBasicPointer(Int_t,ReadFastArray);
      case TStreamerInfo::kSkipP + TStreamerInfo::kLong:     SkipCBasicPointer(Long_t,ReadFastArray);
      case TStreamerInfo::kSkipP + TStreamerInfo::kLong64:   SkipCBasicPointer(Long64_t,ReadFastArray);
      case TStreamerInfo::kSkipP + TStreamerInfo::kFloat:    SkipCBasicPointer(Float_t,ReadFastArray);
      case TStreamerInfo::kSkipP + TStreamerInfo::kFloat16:  SkipCBasicPointer(Float_t,ReadFastArrayFloat16);
      case TStreamerInfo::kSkipP + TStreamerInfo::kDouble:   SkipCBasicPointer(Double_t,ReadFastArray);
      case TStreamerInfo::kSkipP + TStreamerInfo::kDouble32: SkipCBasicPointer(Double_t,ReadFastArrayDouble32)
      case TStreamerInfo::kSkipP + TStreamerInfo::kUChar:    SkipCBasicPointer(UChar_t,ReadFastArray);
      case TStreamerInfo::kSkipP + TStreamerInfo::kUShort:   SkipCBasicPointer(UShort_t,ReadFastArray);
      case TStreamerInfo::kSkipP + TStreamerInfo::kUInt:     SkipCBasicPointer(UInt_t,ReadFastArray);
      case TStreamerInfo::kSkipP + TStreamerInfo::kULong:    SkipCBasicPointer(ULong_t,ReadFastArray);
      case TStreamerInfo::kSkipP + TStreamerInfo::kULong64:  SkipCBasicPointer(ULong64_t,ReadFastArray);

      // skip char*
      case TStreamerInfo::kSkip + TStreamerInfo::kCharStar: {
         DOLOOP {
            Int_t nch; b >> nch;
            if (nch>0) {
               char* readbuf = new char[nch];
               b.ReadFastArray(readbuf,nch);
               delete[] readbuf;
            }
         }
         break;
      }

      // skip Class*   derived from TObject
      case TStreamerInfo::kSkip + TStreamerInfo::kObjectP: {
         DOLOOP{
            for (Int_t j=0;j<compinfo->fLength;j++) {
               b.SkipObjectAny();
            }
         }
         break;
      }

      // skip array counter //[n]
      case TStreamerInfo::kSkip + TStreamerInfo::kCounter: {
         DOLOOP {
            Int_t dummy; b >> dummy;
            aElement->SetTObjectOffset(dummy);
         }
         break;
      }

      // skip Class *  derived from TObject with comment field  //->
      // skip Class    derived from TObject
      case TStreamerInfo::kSkip + TStreamerInfo::kObjectp:
      case TStreamerInfo::kSkip + TStreamerInfo::kObject:  {
         if (cle == TRef::Class()) {
            TRef refjunk;
            DOLOOP{ refjunk.Streamer(b);}
         } else {
            DOLOOP{
               b.SkipObjectAny();
            }
         }
         break;
      }

      // skip Special case for TString, TObject, TNamed
      case TStreamerInfo::kSkip + TStreamerInfo::kTString: {
         TString s;
         DOLOOP {
            s.Streamer(b);
         }
         break;
      }
      case TStreamerInfo::kSkip + TStreamerInfo::kTObject: {
         TObject x;
         DOLOOP {
            x.Streamer(b);
         }
         break;
      }
      case TStreamerInfo::kSkip + TStreamerInfo::kTNamed:  {
         TNamed n;
         DOLOOP {
            n.Streamer(b);
         }
      break;
      }

      // skip Class *  not derived from TObject with comment field  //->
      case TStreamerInfo::kSkip + TStreamerInfo::kAnyp: {
         DOLOOP {
            b.SkipObjectAny();
         }
         break;
      }

      // skip Class*   not derived from TObject
      case TStreamerInfo::kSkip + TStreamerInfo::kAnyP: {
         DOLOOP {
            for (Int_t j=0;j<compinfo->fLength;j++) {
               b.SkipObjectAny();
            }
         }
         break;
      }

      // skip Any Class not derived from TObject
      case TStreamerInfo::kSkip + TStreamerInfo::kAny:    {
         DOLOOP {
            b.SkipObjectAny();
         }
         break;
      }

      // skip Any Class not derived from TObject
      case TStreamerInfo::kSkip + TStreamerInfo::kSTLp:
      case TStreamerInfo::kSkip + TStreamerInfo::kSTLp + TStreamerInfo::kOffsetL:
      case TStreamerInfo::kSkip + TStreamerInfo::kSTL:     {
         if (fOldVersion<3) return 0;
         b.SkipObjectAny();
         break;
      }

      // skip Base Class
      case TStreamerInfo::kSkip + TStreamerInfo::kBase:    {
         DOLOOP {
            b.SkipObjectAny();
         }
         break;
      }

      case TStreamerInfo::kSkip + TStreamerInfo::kStreamLoop:
      case TStreamerInfo::kSkip + TStreamerInfo::kStreamer: {
         DOLOOP {
            b.SkipObjectAny();
            }
         break;
      }
      default:
         //Error("ReadBufferClones","The element type %d is not supported yet\n",compinfo->fType);
         return -1;
   }
   return 0;
}

#define ConvCBasicType(name,stream)                                       \
   {                                                                      \
      DOLOOP {                                                            \
         name u;                                                          \
         stream;                                                          \
         switch(compinfo->fNewType) {                                            \
            case TStreamerInfo::kBool:    {Bool_t   *x=(Bool_t*)(arr[k]+ioffset);   *x = (Bool_t)u;   break;} \
            case TStreamerInfo::kChar:    {Char_t   *x=(Char_t*)(arr[k]+ioffset);   *x = (Char_t)u;   break;} \
            case TStreamerInfo::kShort:   {Short_t  *x=(Short_t*)(arr[k]+ioffset);  *x = (Short_t)u;  break;} \
            case TStreamerInfo::kInt:     {Int_t    *x=(Int_t*)(arr[k]+ioffset);    *x = (Int_t)u;    break;} \
            case TStreamerInfo::kLong:    {Long_t   *x=(Long_t*)(arr[k]+ioffset);   *x = (Long_t)u;   break;} \
            case TStreamerInfo::kLong64:  {Long64_t *x=(Long64_t*)(arr[k]+ioffset); *x = (Long64_t)u; break;} \
            case TStreamerInfo::kFloat:   {Float_t  *x=(Float_t*)(arr[k]+ioffset);  *x = (Float_t)u;  break;} \
            case TStreamerInfo::kFloat16: {Float_t  *x=(Float_t*)(arr[k]+ioffset);  *x = (Float_t)u;  break;} \
            case TStreamerInfo::kDouble:  {Double_t *x=(Double_t*)(arr[k]+ioffset); *x = (Double_t)u; break;} \
            case TStreamerInfo::kDouble32:{Double_t *x=(Double_t*)(arr[k]+ioffset); *x = (Double_t)u; break;} \
            case TStreamerInfo::kUChar:   {UChar_t  *x=(UChar_t*)(arr[k]+ioffset);  *x = (UChar_t)u;  break;} \
            case TStreamerInfo::kUShort:  {UShort_t *x=(UShort_t*)(arr[k]+ioffset); *x = (UShort_t)u; break;} \
            case TStreamerInfo::kUInt:    {UInt_t   *x=(UInt_t*)(arr[k]+ioffset);   *x = (UInt_t)u;   break;} \
            case TStreamerInfo::kULong:   {ULong_t  *x=(ULong_t*)(arr[k]+ioffset);  *x = (ULong_t)u;  break;} \
            case TStreamerInfo::kULong64: {ULong64_t*x=(ULong64_t*)(arr[k]+ioffset);*x = (ULong64_t)u;break;} \
         }                                                                \
      } break;                                                            \
   }

#define ConvCBasicArrayTo(newtype)                                        \
   {                                                                      \
      newtype *f=(newtype*)(arr[k]+ioffset);                              \
      for (j=0;j<len;j++) f[j] = (newtype)readbuf[j];                     \
      break;                                                              \
   }

#define ConvCBasicArray(name,ReadArrayFunc)                               \
   {                                                                      \
      int j, len = compinfo->fLength;                                            \
      name* readbuf = new name[len];                                      \
      int newtype = compinfo->fNewType%20;                                       \
      DOLOOP {                                                            \
          b.ReadArrayFunc(readbuf, len);                                  \
          switch(newtype) {                                               \
             case TStreamerInfo::kBool:     ConvCBasicArrayTo(Bool_t);    \
             case TStreamerInfo::kChar:     ConvCBasicArrayTo(Char_t);    \
             case TStreamerInfo::kShort:    ConvCBasicArrayTo(Short_t);   \
             case TStreamerInfo::kInt:      ConvCBasicArrayTo(Int_t);     \
             case TStreamerInfo::kLong:     ConvCBasicArrayTo(Long_t);    \
             case TStreamerInfo::kLong64:   ConvCBasicArrayTo(Long64_t);  \
             case TStreamerInfo::kFloat:    ConvCBasicArrayTo(Float_t);   \
             case TStreamerInfo::kFloat16:  ConvCBasicArrayTo(Float_t);   \
             case TStreamerInfo::kDouble:   ConvCBasicArrayTo(Double_t);  \
             case TStreamerInfo::kDouble32: ConvCBasicArrayTo(Double_t);  \
             case TStreamerInfo::kUChar:    ConvCBasicArrayTo(UChar_t);   \
             case TStreamerInfo::kUShort:   ConvCBasicArrayTo(UShort_t);  \
             case TStreamerInfo::kUInt:     ConvCBasicArrayTo(UInt_t);    \
             case TStreamerInfo::kULong:    ConvCBasicArrayTo(ULong_t);   \
             case TStreamerInfo::kULong64:  ConvCBasicArrayTo(ULong64_t); \
          }                                                               \
      }                                                                   \
      delete[] readbuf;                                                   \
      break;                                                              \
   }

#define ConvCBasicPointerToOutOfRange(newtype,ReadArrayFunc)              \
   {                                                                      \
     newtype **f=(newtype**)(arr[k]+ioffset);                             \
     for (j=0;j<len;j++) {                                                \
       delete [] f[j];                                                    \
       f[j] = 0;                                                          \
     }                                                                    \
     break;                                                               \
   }

#define ConvCBasicPointerTo(newtype,ReadArrayFunc)                        \
   {                                                                      \
     newtype **f=(newtype**)(arr[k]+ioffset);                             \
     for (j=0;j<len;j++) {                                                \
       delete [] f[j];                                                   \
       f[j] = new newtype[*l];                                            \
       newtype *af = f[j];                                                \
       b.ReadArrayFunc(readbuf, *l);                                      \
       for (jj=0;jj<*l;jj++) af[jj] = (newtype)readbuf[jj];               \
     }                                                                    \
     break;                                                               \
   }

#define ConvCBasicPointer(name,ReadArrayFunc)                                           \
   {                                                                      \
      Char_t isArray;                                                     \
      int j, jj, len = aElement->GetArrayDim()?aElement->GetArrayLength():1; \
      name* readbuf = 0;                                                  \
      int newtype = compinfo->fNewType %20;                                      \
      Int_t imethod = compinfo->fMethod+eoffset;                                 \
      DOLOOP {                                                            \
         b >> isArray;                                                    \
         Int_t *l = (Int_t*)(arr[k]+imethod);                             \
         if (*l>0 && *l < b.BufferSize()) {                               \
            readbuf = new name[*l];                                       \
            switch(newtype) {                                             \
               case TStreamerInfo::kBool:     ConvCBasicPointerTo(Bool_t,ReadArrayFunc);   \
               case TStreamerInfo::kChar:     ConvCBasicPointerTo(Char_t,ReadArrayFunc);   \
               case TStreamerInfo::kShort:    ConvCBasicPointerTo(Short_t,ReadArrayFunc);  \
               case TStreamerInfo::kInt:      ConvCBasicPointerTo(Int_t,ReadArrayFunc);    \
               case TStreamerInfo::kLong:     ConvCBasicPointerTo(Long_t,ReadArrayFunc);   \
               case TStreamerInfo::kLong64:   ConvCBasicPointerTo(Long64_t,ReadArrayFunc); \
               case TStreamerInfo::kFloat:    ConvCBasicPointerTo(Float_t,ReadArrayFunc);  \
               case TStreamerInfo::kFloat16:  ConvCBasicPointerTo(Float_t,ReadArrayFunc);  \
               case TStreamerInfo::kDouble:   ConvCBasicPointerTo(Double_t,ReadArrayFunc); \
               case TStreamerInfo::kDouble32: ConvCBasicPointerTo(Double_t,ReadArrayFunc); \
               case TStreamerInfo::kUChar:    ConvCBasicPointerTo(UChar_t,ReadArrayFunc);  \
               case TStreamerInfo::kUShort:   ConvCBasicPointerTo(UShort_t,ReadArrayFunc); \
               case TStreamerInfo::kUInt:     ConvCBasicPointerTo(UInt_t,ReadArrayFunc);   \
               case TStreamerInfo::kULong:    ConvCBasicPointerTo(ULong_t,ReadArrayFunc);  \
               case TStreamerInfo::kULong64:  ConvCBasicPointerTo(ULong64_t,ReadArrayFunc); \
            }                                                             \
            delete[] readbuf;                                             \
         } else {                                                         \
            switch(newtype) {                                             \
               case TStreamerInfo::kBool:     ConvCBasicPointerToOutOfRange(Bool_t,ReadArrayFunc);   \
               case TStreamerInfo::kChar:     ConvCBasicPointerToOutOfRange(Char_t,ReadArrayFunc);   \
               case TStreamerInfo::kShort:    ConvCBasicPointerToOutOfRange(Short_t,ReadArrayFunc);  \
               case TStreamerInfo::kInt:      ConvCBasicPointerToOutOfRange(Int_t,ReadArrayFunc);    \
               case TStreamerInfo::kLong:     ConvCBasicPointerToOutOfRange(Long_t,ReadArrayFunc);   \
               case TStreamerInfo::kLong64:   ConvCBasicPointerToOutOfRange(Long64_t,ReadArrayFunc); \
               case TStreamerInfo::kFloat:    ConvCBasicPointerToOutOfRange(Float_t,ReadArrayFunc);  \
               case TStreamerInfo::kFloat16:  ConvCBasicPointerToOutOfRange(Float_t,ReadArrayFunc);  \
               case TStreamerInfo::kDouble:   ConvCBasicPointerToOutOfRange(Double_t,ReadArrayFunc); \
               case TStreamerInfo::kDouble32: ConvCBasicPointerToOutOfRange(Double_t,ReadArrayFunc); \
               case TStreamerInfo::kUChar:    ConvCBasicPointerToOutOfRange(UChar_t,ReadArrayFunc);  \
               case TStreamerInfo::kUShort:   ConvCBasicPointerToOutOfRange(UShort_t,ReadArrayFunc); \
               case TStreamerInfo::kUInt:     ConvCBasicPointerToOutOfRange(UInt_t,ReadArrayFunc);   \
               case TStreamerInfo::kULong:    ConvCBasicPointerToOutOfRange(ULong_t,ReadArrayFunc);  \
               case TStreamerInfo::kULong64:  ConvCBasicPointerToOutOfRange(ULong64_t,ReadArrayFunc); \
               }                                                          \
         }                                                                \
         readbuf = 0;                                                     \
      } break;                                                            \
   }

////////////////////////////////////////////////////////////////////////////////
/// Handle Artificial StreamerElement

template <class T>
Int_t TStreamerInfo::ReadBufferArtificial(TBuffer &b, const T &arr,
                                          TStreamerElement *aElement, Int_t narr,
                                          Int_t eoffset)
{
   TStreamerArtificial *artElement = (TStreamerArtificial*)aElement;
   ROOT::TSchemaRule::ReadRawFuncPtr_t rawfunc = artElement->GetReadRawFunc();

   if (rawfunc) {
      for(Int_t k=0; k<narr; ++k) {
         rawfunc( arr[k], b ); // Intentionally pass the object, so that the member can be set from other members.
      }
      return 0;
   }

   ROOT::TSchemaRule::ReadFuncPtr_t readfunc = artElement->GetReadFunc();
   // Process the result
   if (readfunc) {
      TVirtualObject obj(0);
      TVirtualArray *objarr = ((TBufferFile&)b).PeekDataCache();
      if (objarr) {
         obj.fClass = objarr->fClass;

         for(Int_t k=0; k<narr; ++k) {
            obj.fObject = objarr->GetObjectAt(k);
            readfunc(arr[k]+eoffset, &obj);
         }
         obj.fObject = 0; // Prevent auto deletion
      } else {
         for(Int_t k=0; k<narr; ++k) {
            readfunc(arr[k]+eoffset, &obj);
         }
      }
      return 0;
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
///  Convert elements of a TClonesArray

template <class T>
Int_t TStreamerInfo::ReadBufferConv(TBuffer &b, const T &arr,  const TCompInfo *compinfo, Int_t kase,
                                    TStreamerElement *aElement, Int_t narr,
                                    Int_t eoffset)
{
   Int_t ioffset = eoffset+compinfo->fOffset;

   switch (kase) {

      // convert basic types
      case TStreamerInfo::kConv + TStreamerInfo::kBool:    ConvCBasicType(Bool_t,b >> u);
      case TStreamerInfo::kConv + TStreamerInfo::kChar:    ConvCBasicType(Char_t,b >> u);
      case TStreamerInfo::kConv + TStreamerInfo::kShort:   ConvCBasicType(Short_t,b >> u);
      case TStreamerInfo::kConv + TStreamerInfo::kInt:     ConvCBasicType(Int_t,b >> u);
      case TStreamerInfo::kConv + TStreamerInfo::kLong:    if (compinfo->fNewType==TStreamerInfo::kLong64 || compinfo->fNewType==TStreamerInfo::kULong64) {
                                                              ConvCBasicType(Long64_t,b >> u);
                                                           } else {
                                                              ConvCBasicType(Long_t,b >> u);
                                                           }
      case TStreamerInfo::kConv + TStreamerInfo::kLong64:  ConvCBasicType(Long64_t,b >> u);
      case TStreamerInfo::kConv + TStreamerInfo::kFloat:   ConvCBasicType(Float_t,b >> u);
      case TStreamerInfo::kConv + TStreamerInfo::kFloat16: ConvCBasicType(Float_t,b.ReadFloat16(&u,aElement));
      case TStreamerInfo::kConv + TStreamerInfo::kDouble:  ConvCBasicType(Double_t,b >> u);
      case TStreamerInfo::kConv + TStreamerInfo::kDouble32:ConvCBasicType(Double_t,b.ReadDouble32(&u,aElement));
      case TStreamerInfo::kConv + TStreamerInfo::kUChar:   ConvCBasicType(UChar_t,b >> u);
      case TStreamerInfo::kConv + TStreamerInfo::kUShort:  ConvCBasicType(UShort_t,b >> u);
      case TStreamerInfo::kConv + TStreamerInfo::kUInt:    ConvCBasicType(UInt_t,b >> u);
      case TStreamerInfo::kConv + TStreamerInfo::kULong:   if (compinfo->fNewType==TStreamerInfo::kLong64 || compinfo->fNewType==TStreamerInfo::kULong64) {
#if defined(_MSC_VER) && (_MSC_VER <= 1200)
                                                              ConvCBasicType(Long64_t,b >> u);
#else
                                                              ConvCBasicType(ULong64_t,b >> u);
#endif
                                                           } else {
                                                              ConvCBasicType(ULong_t,b >> u);
                                                           }
#if defined(_MSC_VER) && (_MSC_VER <= 1200)
      case TStreamerInfo::kConv + TStreamerInfo::kULong64: ConvCBasicType(Long64_t,b >> u)
#else
      case TStreamerInfo::kConv + TStreamerInfo::kULong64: ConvCBasicType(ULong64_t,b >> u)
#endif
      case TStreamerInfo::kConv + TStreamerInfo::kBits:  {
         DOLOOP {
            UInt_t u;
            b >> u;
            if ((u & kIsReferenced) != 0) {
               UShort_t pidf;
               b >> pidf;
               pidf += b.GetPidOffset();
               TProcessID *pid = b.ReadProcessID(pidf);
               if (pid!=0) {
                  TObject *obj = (TObject*)(arr[k]+eoffset);
                  UInt_t gpid = pid->GetUniqueID();
                  UInt_t uid;
                  if (gpid>=0xff) {
                     uid = obj->GetUniqueID() | 0xff000000;
                  } else {
                     uid = ( obj->GetUniqueID() & 0xffffff) + (gpid<<24);
                  }
                  obj->SetUniqueID(uid);
                  pid->PutObjectWithID(obj);
               }
            }
            switch(compinfo->fNewType) {
               case TStreamerInfo::kBool:    {Bool_t   *x=(Bool_t*)(arr[k]+ioffset);   *x = (Bool_t)u;   break;}
               case TStreamerInfo::kChar:    {Char_t   *x=(Char_t*)(arr[k]+ioffset);   *x = (Char_t)u;   break;}
               case TStreamerInfo::kShort:   {Short_t  *x=(Short_t*)(arr[k]+ioffset);  *x = (Short_t)u;  break;}
               case TStreamerInfo::kInt:     {Int_t    *x=(Int_t*)(arr[k]+ioffset);    *x = (Int_t)u;    break;}
               case TStreamerInfo::kLong:    {Long_t   *x=(Long_t*)(arr[k]+ioffset);   *x = (Long_t)u;   break;}
               case TStreamerInfo::kLong64:  {Long64_t *x=(Long64_t*)(arr[k]+ioffset); *x = (Long64_t)u; break;}
               case TStreamerInfo::kFloat:   {Float_t  *x=(Float_t*)(arr[k]+ioffset);  *x = (Float_t)u;  break;}
               case TStreamerInfo::kFloat16: {Float_t  *x=(Float_t*)(arr[k]+ioffset);  *x = (Float_t)u;  break;}
               case TStreamerInfo::kDouble:  {Double_t *x=(Double_t*)(arr[k]+ioffset); *x = (Double_t)u; break;}
               case TStreamerInfo::kDouble32:{Double_t *x=(Double_t*)(arr[k]+ioffset); *x = (Double_t)u; break;}
               case TStreamerInfo::kUChar:   {UChar_t  *x=(UChar_t*)(arr[k]+ioffset);  *x = (UChar_t)u;  break;}
               case TStreamerInfo::kUShort:  {UShort_t *x=(UShort_t*)(arr[k]+ioffset); *x = (UShort_t)u; break;}
               case TStreamerInfo::kUInt:    {UInt_t   *x=(UInt_t*)(arr[k]+ioffset);   *x = (UInt_t)u;   break;}
               case TStreamerInfo::kULong:   {ULong_t  *x=(ULong_t*)(arr[k]+ioffset);  *x = (ULong_t)u;  break;}
               case TStreamerInfo::kULong64: {ULong64_t*x=(ULong64_t*)(arr[k]+ioffset);*x = (ULong64_t)u;break;}
            }
         } break;
      }

         // convert array of basic types  array[8]
      case TStreamerInfo::kConvL + TStreamerInfo::kBool:    ConvCBasicArray(Bool_t,ReadFastArray);
      case TStreamerInfo::kConvL + TStreamerInfo::kChar:    ConvCBasicArray(Char_t,ReadFastArray);
      case TStreamerInfo::kConvL + TStreamerInfo::kShort:   ConvCBasicArray(Short_t,ReadFastArray);
      case TStreamerInfo::kConvL + TStreamerInfo::kInt:     ConvCBasicArray(Int_t,ReadFastArray);
      case TStreamerInfo::kConvL + TStreamerInfo::kLong:
         if (compinfo->fNewType==TStreamerInfo::kLong64 || compinfo->fNewType==TStreamerInfo::kULong64) {
            ConvCBasicArray(Long64_t,ReadFastArray);
         } else {
            ConvCBasicArray(Long_t,ReadFastArray);
         }
      case TStreamerInfo::kConvL + TStreamerInfo::kLong64:  ConvCBasicArray(Long64_t,ReadFastArray);
      case TStreamerInfo::kConvL + TStreamerInfo::kFloat:   ConvCBasicArray(Float_t,ReadFastArray);
      case TStreamerInfo::kConvL + TStreamerInfo::kFloat16: ConvCBasicArray(Float_t,ReadFastArrayFloat16);
      case TStreamerInfo::kConvL + TStreamerInfo::kDouble:  ConvCBasicArray(Double_t,ReadFastArray);
      case TStreamerInfo::kConvL + TStreamerInfo::kDouble32:ConvCBasicArray(Double_t,ReadFastArrayDouble32);
      case TStreamerInfo::kConvL + TStreamerInfo::kUChar:   ConvCBasicArray(UChar_t,ReadFastArray);
      case TStreamerInfo::kConvL + TStreamerInfo::kUShort:  ConvCBasicArray(UShort_t,ReadFastArray);
      case TStreamerInfo::kConvL + TStreamerInfo::kUInt:    ConvCBasicArray(UInt_t,ReadFastArray);
      case TStreamerInfo::kConvL + TStreamerInfo::kULong:
         if (compinfo->fNewType==TStreamerInfo::kLong64 || compinfo->fNewType==TStreamerInfo::kULong64) {
#if defined(_MSC_VER) && (_MSC_VER <= 1200)
            ConvCBasicArray(Long64_t,ReadFastArray)
#else
            ConvCBasicArray(ULong64_t,ReadFastArray)
#endif
         } else {
            ConvCBasicArray(ULong_t,ReadFastArray);
         }
#if defined(_MSC_VER) && (_MSC_VER <= 1200)
      case TStreamerInfo::kConvL + TStreamerInfo::kULong64: ConvCBasicArray(Long64_t,ReadFastArray)
#else
      case TStreamerInfo::kConvL + TStreamerInfo::kULong64: ConvCBasicArray(ULong64_t,ReadFastArray)
#endif

   // convert pointer to an array of basic types  array[n]
      case TStreamerInfo::kConvP + TStreamerInfo::kBool:    ConvCBasicPointer(Bool_t,ReadFastArray);
      case TStreamerInfo::kConvP + TStreamerInfo::kChar:    ConvCBasicPointer(Char_t,ReadFastArray);
      case TStreamerInfo::kConvP + TStreamerInfo::kShort:   ConvCBasicPointer(Short_t,ReadFastArray);
      case TStreamerInfo::kConvP + TStreamerInfo::kInt:     ConvCBasicPointer(Int_t,ReadFastArray);
      case TStreamerInfo::kConvP + TStreamerInfo::kLong:
         if (compinfo->fNewType==TStreamerInfo::kLong64 || compinfo->fNewType==TStreamerInfo::kULong64) {
            ConvCBasicPointer(Long64_t,ReadFastArray);
         } else {
            ConvCBasicPointer(Long_t,ReadFastArray);
         }
      case TStreamerInfo::kConvP + TStreamerInfo::kLong64:  ConvCBasicPointer(Long64_t,ReadFastArray);
      case TStreamerInfo::kConvP + TStreamerInfo::kFloat:   ConvCBasicPointer(Float_t,ReadFastArray);
      case TStreamerInfo::kConvP + TStreamerInfo::kFloat16: ConvCBasicPointer(Float_t,ReadFastArrayFloat16);
      case TStreamerInfo::kConvP + TStreamerInfo::kDouble:  ConvCBasicPointer(Double_t,ReadFastArray);
      case TStreamerInfo::kConvP + TStreamerInfo::kDouble32:ConvCBasicPointer(Double_t,ReadFastArrayDouble32);
      case TStreamerInfo::kConvP + TStreamerInfo::kUChar:   ConvCBasicPointer(UChar_t,ReadFastArray);
      case TStreamerInfo::kConvP + TStreamerInfo::kUShort:  ConvCBasicPointer(UShort_t,ReadFastArray);
      case TStreamerInfo::kConvP + TStreamerInfo::kUInt:    ConvCBasicPointer(UInt_t,ReadFastArray);
      case TStreamerInfo::kConvP + TStreamerInfo::kULong:
         if (compinfo->fNewType==TStreamerInfo::kLong64 || compinfo->fNewType==TStreamerInfo::kULong64) {
#if defined(_MSC_VER) && (_MSC_VER <= 1200)
            ConvCBasicPointer(Long64_t,ReadFastArray)
#else
            ConvCBasicPointer(ULong64_t,ReadFastArray)
#endif
         } else {
            ConvCBasicPointer(ULong_t,ReadFastArray);
         }
#if defined(_MSC_VER) && (_MSC_VER <= 1200)
      case TStreamerInfo::kConvP + TStreamerInfo::kULong64: ConvCBasicPointer(Long64_t,ReadFastArray)
#else
      case TStreamerInfo::kConvP + TStreamerInfo::kULong64: ConvCBasicPointer(ULong64_t,ReadFastArray)
#endif

      default:
         // Warning("ReadBufferConv","The element type %d is not supported yet",compinfo->fType);
         return -1;

   }

   return 0;
}

// Helper function for TStreamerInfo::ReadBuffer
namespace {
   template <class T> Bool_t R__TestUseCache(TStreamerElement *element)
   {
      return element->TestBit(TStreamerElement::kCache);
   }

   template <> Bool_t R__TestUseCache<TVirtualArray>(TStreamerElement*)
   {
      // We are already using the cache, no need to recurse one more time.
      return kFALSE;
   }
}

////////////////////////////////////////////////////////////////////////////////
///  Deserialize information from buffer b into object at pointer
///  if (arrayMode & 1) ptr is a pointer to array of pointers to the objects
///  otherwise it is a pointer to a pointer to a single object.
///  This also means that T is of a type such that arr[i] is a pointer to an
///  object.  Currently the only anticipated instantiation are for T==char**
///  and T==TVirtualCollectionProxy

template <class T>
Int_t TStreamerInfo::ReadBuffer(TBuffer &b, const T &arr,
                                TCompInfo *const*const compinfo, Int_t first, Int_t last,
                                Int_t narr, Int_t eoffset, Int_t arrayMode)
{
   TStreamerInfo *thisVar = this;
   Bool_t needIncrement = !( arrayMode & 2 );
   arrayMode = arrayMode & (~2);

   if (needIncrement) b.IncrementLevel(thisVar);

   //loop on all active members

   // In order to speed up the case where the object being written is
   // not in a collection (i.e. arrayMode is false), we actually
   // duplicate the code for the elementary types using this typeOffset.
   static const int kHaveLoop = 1024;
   const Int_t typeOffset = arrayMode ? kHaveLoop : 0;

   TClass     *cle      = 0;
   TClass     *newCle   = 0;
   TMemberStreamer *pstreamer=0;
   Int_t isPreAlloc = 0;
   for (Int_t i=first;i<last;i++) {
      TStreamerElement * aElement  = (TStreamerElement*)compinfo[i]->fElem;
      CurrentElement() = aElement;

      if (needIncrement) b.SetStreamerElementNumber(aElement,compinfo[i]->fType);

      if (aElement->TestBit(TStreamerElement::kWrite)) continue;

      if (R__TestUseCache<T>(aElement)) {
         Int_t bufpos = b.Length();
         if (((TBufferFile&)b).PeekDataCache()==0) {
            Warning("ReadBuffer","Skipping %s::%s because the cache is missing.",thisVar->GetName(),aElement->GetName());
            thisVar->ReadBufferSkip(b,arr,compinfo[i],compinfo[i]->fType+TStreamerInfo::kSkip,aElement,narr,eoffset);
         } else {
            if (gDebug > 1) {
               printf("ReadBuffer, class:%s, name=%s, fType[%d]=%d,"
                  " %s, bufpos=%d, arr=%p, eoffset=%d, Redirect=%p\n",
                  fClass->GetName(),aElement->GetName(),i,compinfo[i]->fType,
                  aElement->ClassName(),b.Length(),arr[0], eoffset,((TBufferFile&)b).PeekDataCache()->GetObjectAt(0));
            }
            thisVar->ReadBuffer(b,*((TBufferFile&)b).PeekDataCache(),compinfo,i,i+1,narr,eoffset, arrayMode);
         }
         if (aElement->TestBit(TStreamerElement::kRepeat)) { b.SetBufferOffset(bufpos); }
         continue;
      }
      const Int_t ioffset = compinfo[i]->fOffset+eoffset;

      if (gDebug > 1) {
         printf("ReadBuffer, class:%s, name=%s, fType[%d]=%d,"
                " %s, bufpos=%d, arr=%p, offset=%d\n",
                fClass->GetName(),aElement->GetName(),i,compinfo[i]->fType,
                aElement->ClassName(),b.Length(),arr[0], ioffset);
      }

      Int_t kase = compinfo[i]->fType;

      switch (kase + typeOffset) {

         // read basic types
         case TStreamerInfo::kBool:               ReadBasicType(Bool_t);    continue;
         case TStreamerInfo::kChar:               ReadBasicType(Char_t);    continue;
         case TStreamerInfo::kShort:              ReadBasicType(Short_t);   continue;
         case TStreamerInfo::kInt:                ReadBasicType(Int_t);     continue;
         case TStreamerInfo::kLong:               ReadBasicType(Long_t);    continue;
         case TStreamerInfo::kLong64:             ReadBasicType(Long64_t);  continue;
         case TStreamerInfo::kFloat:              ReadBasicType(Float_t);   continue;
         case TStreamerInfo::kDouble:             ReadBasicType(Double_t);  continue;
         case TStreamerInfo::kUChar:              ReadBasicType(UChar_t);   continue;
         case TStreamerInfo::kUShort:             ReadBasicType(UShort_t);  continue;
         case TStreamerInfo::kUInt:               ReadBasicType(UInt_t);    continue;
         case TStreamerInfo::kULong:              ReadBasicType(ULong_t);   continue;
         case TStreamerInfo::kULong64:            ReadBasicType(ULong64_t); continue;
         case TStreamerInfo::kFloat16: {
            Float_t *x=(Float_t*)(arr[0]+ioffset);
            b.ReadFloat16(x,aElement);
            continue;
         }
         case TStreamerInfo::kDouble32: {
            Double_t *x=(Double_t*)(arr[0]+ioffset);
            b.ReadDouble32(x,aElement);
            continue;
         }

         case TStreamerInfo::kBool   + kHaveLoop: ReadBasicTypeLoop(Bool_t);    continue;
         case TStreamerInfo::kChar   + kHaveLoop: ReadBasicTypeLoop(Char_t);    continue;
         case TStreamerInfo::kShort  + kHaveLoop: ReadBasicTypeLoop(Short_t);   continue;
         case TStreamerInfo::kInt    + kHaveLoop: ReadBasicTypeLoop(Int_t);     continue;
         case TStreamerInfo::kLong   + kHaveLoop: ReadBasicTypeLoop(Long_t);    continue;
         case TStreamerInfo::kLong64 + kHaveLoop: ReadBasicTypeLoop(Long64_t);  continue;
         case TStreamerInfo::kFloat  + kHaveLoop: ReadBasicTypeLoop(Float_t);   continue;
         case TStreamerInfo::kDouble + kHaveLoop: ReadBasicTypeLoop(Double_t);  continue;
         case TStreamerInfo::kUChar  + kHaveLoop: ReadBasicTypeLoop(UChar_t);   continue;
         case TStreamerInfo::kUShort + kHaveLoop: ReadBasicTypeLoop(UShort_t);  continue;
         case TStreamerInfo::kUInt   + kHaveLoop: ReadBasicTypeLoop(UInt_t);    continue;
         case TStreamerInfo::kULong  + kHaveLoop: ReadBasicTypeLoop(ULong_t);   continue;
         case TStreamerInfo::kULong64+ kHaveLoop: ReadBasicTypeLoop(ULong64_t); continue;
         case TStreamerInfo::kFloat16 + kHaveLoop: {
            for(Int_t k=0; k<narr; ++k) {
               Float_t *x=(Float_t*)(arr[k]+ioffset);
               b.ReadFloat16(x,aElement);
            }
            continue;
         }
         case TStreamerInfo::kDouble32 + kHaveLoop: {
            for(Int_t k=0; k<narr; ++k) {
               Double_t *x=(Double_t*)(arr[k]+ioffset);
               b.ReadDouble32(x,aElement);
            }
            continue;
         }

         // read array of basic types  like array[8]
         case TStreamerInfo::kOffsetL + TStreamerInfo::kBool:   ReadBasicArray(Bool_t);    continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kChar:   ReadBasicArray(Char_t);    continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kShort:  ReadBasicArray(Short_t);   continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kInt:    ReadBasicArray(Int_t);     continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kLong:   ReadBasicArray(Long_t);    continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kLong64: ReadBasicArray(Long64_t);  continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kFloat:  ReadBasicArray(Float_t);   continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kDouble: ReadBasicArray(Double_t);  continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kUChar:  ReadBasicArray(UChar_t);   continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kUShort: ReadBasicArray(UShort_t);  continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kUInt:   ReadBasicArray(UInt_t);    continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kULong:  ReadBasicArray(ULong_t);   continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kULong64:ReadBasicArray(ULong64_t); continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kFloat16: {
            b.ReadFastArrayFloat16((Float_t*)(arr[0]+ioffset),compinfo[i]->fLength,aElement);
            continue;
         }
         case TStreamerInfo::kOffsetL + TStreamerInfo::kDouble32: {
            b.ReadFastArrayDouble32((Double_t*)(arr[0]+ioffset),compinfo[i]->fLength,aElement);
            continue;
         }

         case TStreamerInfo::kOffsetL + TStreamerInfo::kBool    + kHaveLoop: ReadBasicArrayLoop(Bool_t);    continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kChar    + kHaveLoop: ReadBasicArrayLoop(Char_t);    continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kShort   + kHaveLoop: ReadBasicArrayLoop(Short_t);   continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kInt     + kHaveLoop: ReadBasicArrayLoop(Int_t);     continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kLong    + kHaveLoop: ReadBasicArrayLoop(Long_t);    continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kLong64  + kHaveLoop: ReadBasicArrayLoop(Long64_t);  continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kFloat   + kHaveLoop: ReadBasicArrayLoop(Float_t);   continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kDouble  + kHaveLoop: ReadBasicArrayLoop(Double_t);  continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kUChar   + kHaveLoop: ReadBasicArrayLoop(UChar_t);   continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kUShort  + kHaveLoop: ReadBasicArrayLoop(UShort_t);  continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kUInt    + kHaveLoop: ReadBasicArrayLoop(UInt_t);    continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kULong   + kHaveLoop: ReadBasicArrayLoop(ULong_t);   continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kULong64 + kHaveLoop: ReadBasicArrayLoop(ULong64_t); continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kFloat16 + kHaveLoop: {
            for(Int_t k=0; k<narr; ++k) {
               b.ReadFastArrayFloat16((Float_t*)(arr[k]+ioffset),compinfo[i]->fLength,aElement);
            }
            continue;
         }
         case TStreamerInfo::kOffsetL + TStreamerInfo::kDouble32+ kHaveLoop: {
            for(Int_t k=0; k<narr; ++k) {
               b.ReadFastArrayDouble32((Double_t*)(arr[k]+ioffset),compinfo[i]->fLength,aElement);
            }
            continue;
         }

         // read pointer to an array of basic types  array[n]
         case TStreamerInfo::kOffsetP + TStreamerInfo::kBool:   ReadBasicPointer(Bool_t);  continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kChar:   ReadBasicPointer(Char_t);  continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kShort:  ReadBasicPointer(Short_t);  continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kInt:    ReadBasicPointer(Int_t);  continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kLong:   ReadBasicPointer(Long_t);  continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kLong64: ReadBasicPointer(Long64_t);  continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kFloat:  ReadBasicPointer(Float_t);  continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kDouble: ReadBasicPointer(Double_t);  continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kUChar:  ReadBasicPointer(UChar_t);  continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kUShort: ReadBasicPointer(UShort_t);  continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kUInt:   ReadBasicPointer(UInt_t);  continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kULong:  ReadBasicPointer(ULong_t);  continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kULong64:ReadBasicPointer(ULong64_t);  continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kFloat16: {
            Char_t isArray;
            b >> isArray;
            const int imethod = compinfo[i]->fMethod+eoffset;
            Int_t *l = (Int_t*)(arr[0]+imethod);
            Float_t **f = (Float_t**)(arr[0]+ioffset);
            int j;
            for(j=0;j<compinfo[i]->fLength;j++) {
               delete [] f[j];
               f[j] = 0; if (*l <=0) continue;
               f[j] = new Float_t[*l];
               b.ReadFastArrayFloat16(f[j],*l,aElement);
            }
            continue;
         }
         case TStreamerInfo::kOffsetP + TStreamerInfo::kDouble32: {
            Char_t isArray;
            b >> isArray;
            const int imethod = compinfo[i]->fMethod+eoffset;
            Int_t *l = (Int_t*)(arr[0]+imethod);
            Double_t **f = (Double_t**)(arr[0]+ioffset);
            int j;
            for(j=0;j<compinfo[i]->fLength;j++) {
               delete [] f[j];
               f[j] = 0; if (*l <=0) continue;
               f[j] = new Double_t[*l];
               b.ReadFastArrayDouble32(f[j],*l,aElement);
            }
            continue;
         }

         case TStreamerInfo::kOffsetP + TStreamerInfo::kBool    + kHaveLoop: ReadBasicPointerLoop(Bool_t);    continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kChar    + kHaveLoop: ReadBasicPointerLoop(Char_t);    continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kShort   + kHaveLoop: ReadBasicPointerLoop(Short_t);   continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kInt     + kHaveLoop: ReadBasicPointerLoop(Int_t);     continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kLong    + kHaveLoop: ReadBasicPointerLoop(Long_t);    continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kLong64  + kHaveLoop: ReadBasicPointerLoop(Long64_t);  continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kFloat   + kHaveLoop: ReadBasicPointerLoop(Float_t);   continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kDouble  + kHaveLoop: ReadBasicPointerLoop(Double_t);  continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kUChar   + kHaveLoop: ReadBasicPointerLoop(UChar_t);   continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kUShort  + kHaveLoop: ReadBasicPointerLoop(UShort_t);  continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kUInt    + kHaveLoop: ReadBasicPointerLoop(UInt_t);    continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kULong   + kHaveLoop: ReadBasicPointerLoop(ULong_t);   continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kULong64 + kHaveLoop: ReadBasicPointerLoop(ULong64_t); continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kFloat16 + kHaveLoop: {
            const int imethod = compinfo[i]->fMethod+eoffset;
            for(Int_t k=0; k<narr; ++k) {
               Char_t isArray;
               b >> isArray;
               Int_t *l = (Int_t*)(arr[k]+imethod);
               Float_t **f = (Float_t**)(arr[k]+ioffset);
               int j;
               for(j=0;j<compinfo[i]->fLength;j++) {
                  delete [] f[j];
                  f[j] = 0; if (*l <=0) continue;
                  f[j] = new Float_t[*l];
                  b.ReadFastArrayFloat16(f[j],*l,aElement);
               }
            }
            continue;
         }
         case TStreamerInfo::kOffsetP + TStreamerInfo::kDouble32+ kHaveLoop: {
            const int imethod = compinfo[i]->fMethod+eoffset;
            for(Int_t k=0; k<narr; ++k) {
               Char_t isArray;
               b >> isArray;
               Int_t *l = (Int_t*)(arr[k]+imethod);
               Double_t **f = (Double_t**)(arr[k]+ioffset);
               int j;
               for(j=0;j<compinfo[i]->fLength;j++) {
                  delete [] f[j];
                  f[j] = 0; if (*l <=0) continue;
                  f[j] = new Double_t[*l];
                  b.ReadFastArrayDouble32(f[j],*l,aElement);
               }
            }
            continue;
         }
      }

      switch (kase) {

         // char*
         case TStreamerInfo::kCharStar: {
            DOLOOP {
               char **f = (char**)(arr[k]+ioffset);
               b.ReadCharStar(*f);
            }
         }
         continue;

         // special case for TObject::fBits in case of a referenced object
         case TStreamerInfo::kBits: {
            DOLOOP {
               UInt_t *x=(UInt_t*)(arr[k]+ioffset); b >> *x;
               if ((*x & kIsReferenced) != 0) {
                  UShort_t pidf;
                  b >> pidf;
                  pidf += b.GetPidOffset();
                  TProcessID *pid = b.ReadProcessID(pidf);
                  if (pid!=0) {
                     TObject *obj = (TObject*)(arr[k]+eoffset);
                     UInt_t gpid = pid->GetUniqueID();
                     UInt_t uid;
                     if (gpid>=0xff) {
                        uid = obj->GetUniqueID() | 0xff000000;
                     } else {
                        uid = ( obj->GetUniqueID() & 0xffffff) + (gpid<<24);
                     }
                     obj->SetUniqueID(uid);
                     pid->PutObjectWithID(obj);
                  }
               }
            }
         }
         continue;

         // array counter //[n]
         case TStreamerInfo::kCounter: {
            DOLOOP {
               Int_t *x=(Int_t*)(arr[k]+ioffset);
               b >> *x;
            }
         }
         continue;


         // Special case for TString, TObject, TNamed
         case TStreamerInfo::kTString: { DOLOOP { ((TString*)(arr[k]+ioffset))->Streamer(b);         } } continue;
         case TStreamerInfo::kTObject: { DOLOOP { ((TObject*)(arr[k]+ioffset))->TObject::Streamer(b);} } continue;
         case TStreamerInfo::kTNamed:  { DOLOOP { ((TNamed*) (arr[k]+ioffset))->TNamed::Streamer(b) ;} } continue;

      }

   SWIT:
      isPreAlloc= 0;
      cle       = compinfo[i]->fClass;
      newCle    = compinfo[i]->fNewClass;
      pstreamer = compinfo[i]->fStreamer;

      switch (kase) {

         case TStreamerInfo::kAnyp:    // Class*  not derived from TObject with    comment field //->
         case TStreamerInfo::kAnyp+TStreamerInfo::kOffsetL:
         case TStreamerInfo::kObjectp: // Class*      derived from TObject with    comment field  //->
         case TStreamerInfo::kObjectp+TStreamerInfo::kOffsetL:
            isPreAlloc = 1;
            // Intentional fallthrough now that isPreAlloc is set.
         case TStreamerInfo::kObjectP: // Class* derived from TObject with no comment field NOTE: Re-added by Phil
         case TStreamerInfo::kObjectP+TStreamerInfo::kOffsetL:
         case TStreamerInfo::kAnyP:    // Class*  not derived from TObject with no comment field NOTE:: Re-added by Phil
         case TStreamerInfo::kAnyP+TStreamerInfo::kOffsetL: {
            DOLOOP {
               b.ReadFastArray((void**)(arr[k]+ioffset),cle,compinfo[i]->fLength,isPreAlloc,pstreamer);
            }
         }
         continue;

//        case TStreamerInfo::kSTLvarp:           // Variable size array of STL containers.
//             {
//                TMemberStreamer *pstreamer = compinfo[i]->fStreamer;
//                TClass *cl                 = compinfo[i]->fClass;
//                ROOT::NewArrFunc_t arraynew = cl->GetNewArray();
//                ROOT::DelArrFunc_t arraydel = cl->GetDeleteArray();
//                UInt_t start,count;
//                // Version_t v =
//                b.ReadVersion(&start, &count, cle);
//                if (pstreamer == 0) {
//                   Int_t size = cl->Size();
//                   Int_t imethod = compinfo[i]->fMethod+eoffset;
//                   DOLOOP {
//                      char **contp = (char**)(arr[k]+ioffset);
//                      const Int_t *counter = (Int_t*)(arr[k]+imethod);
//                      const Int_t sublen = (*counter);

//                      for(int j=0;j<compinfo[i]->fLength;++j) {
//                         if (arraydel) arraydel(contp[j]);
//                         contp[j] = 0;
//                         if (sublen<=0) continue;
//                         if (arraynew) {
//                            contp[j] = (char*)arraynew(sublen, 0);
//                            char *cont = contp[j];
//                            for(int k=0;k<sublen;++k) {
//                               cl->Streamer( cont, b );
//                               cont += size;
//                            }
//                         } else {
//                            // Can't create an array of object
//                            Error("ReadBuffer","The element %s::%s type %d (%s) can be read because of the class does not have access to new %s[..]\n",
//                                  GetName(),aElement->GetFullName(),kase,aElement->GetTypeName(),GetName());
//                            void *cont = cl->New();
//                            for(int k=0;k<sublen;++k) {
//                               cl->Streamer( cont, b );
//                            }
//                         }
//                      }
//                   }
//                } else {
//                   DOLOOP{(*pstreamer)(b,arr[k]+ioffset,compinfo[i]->fLength);}
//                }
//                b.CheckByteCount(start,count,aElement->GetFullName());
//             }
//             continue;

         case TStreamerInfo::kSTLp:            // Pointer to Container with no virtual table (stl) and no comment
         case TStreamerInfo::kSTLp + TStreamerInfo::kOffsetL: // array of pointers to Container with no virtual table (stl) and no comment
            {
               UInt_t start,count;
               Version_t vers = b.ReadVersion(&start, &count, cle);

               if ( vers & TBufferFile::kStreamedMemberWise ) {
                  // Collection was saved member-wise

                  vers &= ~( TBufferFile::kStreamedMemberWise );

                  TClass *newClass = aElement->GetNewClass();
                  TClass *oldClass = aElement->GetClassPointer();
                  if( vers < 9 && newClass && newClass!=oldClass ) {
                     Error( "ReadBuffer", "Unfortunately, version %d of TStreamerInfo (used in %s) did not record enough information to convert a %s into a %s.",
                           vers, b.GetParent() ? b.GetParent()->GetName() : "memory/socket", oldClass->GetName(), newClass->GetName() );
                     continue;
                  }

                  Version_t vClVersion = 0; // For vers less than 9, we have to use the current version.
                  if( vers >= 9 ) {
                     vClVersion = b.ReadVersionForMemberWise( cle->GetCollectionProxy()->GetValueClass() );
                  }

                  TVirtualCollectionProxy *newProxy = (newClass ? newClass->GetCollectionProxy() : 0);
                  TVirtualCollectionProxy *oldProxy = oldClass->GetCollectionProxy();
                  TStreamerInfo *subinfo = 0;

                  if( newProxy ) {
                     // coverity[dereference] oldProxy->GetValueClass() can not be null since this was streamed memberwise.
                     subinfo = (TStreamerInfo*)newProxy->GetValueClass()->GetConversionStreamerInfo( oldProxy->GetValueClass(), vClVersion );
                  } else {
                     subinfo = (TStreamerInfo*)oldProxy->GetValueClass()->GetStreamerInfo( vClVersion );
                     newProxy = oldProxy;
                  }
                  if (subinfo) {
                     DOLOOP {
                        void* env;
                        void **contp = (void**)(arr[k]+ioffset);
                        int j;
                        for(j=0;j<compinfo[i]->fLength;j++) {
                           void *cont = contp[j];
                           if (cont==0) {
                              contp[j] = cle->New();
                              cont = contp[j];
                           }
                           TVirtualCollectionProxy::TPushPop helper( newProxy, cont );
                           Int_t nobjects;
                           b >> nobjects;
                           env = newProxy->Allocate(nobjects,true);
                           subinfo->ReadBufferSTL(b,newProxy,nobjects,/* offset */ 0, vers>=7 );
                           newProxy->Commit(env);
                        }
                     }
                  }
                  b.CheckByteCount(start,count,aElement->GetFullName());
                  continue;
               }
               if (pstreamer == 0) {
                  DOLOOP {
                     void **contp = (void**)(arr[k]+ioffset);
                     int j;
                     for(j=0;j<compinfo[i]->fLength;j++) {
                        void *cont = contp[j];
                        if (cont==0) {
                           // int R__n;
                           // b >> R__n;
                           // b.SetOffset(b.GetOffset()-4); // rewind to the start of the int
                           // if (R__n) continue;
                           contp[j] = cle->New();
                           cont = contp[j];
                        }
                        cle->Streamer( cont, b );
                     }
                  }
               } else {
                  DOLOOP {(*pstreamer)(b,arr[k]+ioffset,compinfo[i]->fLength);}
               }
               b.CheckByteCount(start,count,aElement->GetFullName());
            }
            continue;

         case TStreamerInfo::kSTL:                // Container with no virtual table (stl) and no comment
         case TStreamerInfo::kSTL + TStreamerInfo::kOffsetL:     // array of Container with no virtual table (stl) and no comment
            {
               UInt_t start, count;
               Version_t vers = b.ReadVersion(&start, &count, cle);

               if ( vers & TBufferFile::kStreamedMemberWise ) {
                  // Collection was saved member-wise
                  vers &= ~( TBufferFile::kStreamedMemberWise );

                  TClass *newClass = aElement->GetNewClass();
                  TClass *oldClass = aElement->GetClassPointer();

                  if( vers < 8 && newClass && newClass!=oldClass ) {
                     Error( "ReadBuffer", "Unfortunately, version %d of TStreamerInfo (used in %s) did not record enough information to convert a %s into a %s.",
                           vers, b.GetParent() ? b.GetParent()->GetName() : "memory/socket", oldClass->GetName(), newClass->GetName() );
                     continue;
                  }
                  TVirtualCollectionProxy *oldProxy = oldClass->GetCollectionProxy();
                  TClass *valueClass = oldProxy ? oldProxy->GetValueClass() : 0;
                  Version_t vClVersion = 0; // For vers less than 8, we have to use the current version.
                  if( vers >= 8 ) {
                     vClVersion = b.ReadVersionForMemberWise( valueClass );
                  }

                  if (valueClass == 0) {
                     // MemberWise streaming applies to only collection of classes, and hence
                     // valueClass can only be null if we are reading without the original library
                     // and the collection is always empty,
                     // So let's skip the rest (which requires the StreamerInfo of the valueClass ... which we do not have)

                     b.SetBufferOffset(start+count+sizeof(UInt_t));
                     continue;
                  }

                  TVirtualCollectionProxy *newProxy = (newClass ? newClass->GetCollectionProxy() : 0);
                  TStreamerInfo *subinfo = 0;

                  if( newProxy ) {
                     // coverity[dereference] oldProxy->GetValueClass() can not be null since this was streamed memberwise.
                     subinfo = (TStreamerInfo*)newProxy->GetValueClass()->GetConversionStreamerInfo( oldProxy->GetValueClass(), vClVersion );
                  } else {
                     subinfo = (TStreamerInfo*)valueClass->GetStreamerInfo( vClVersion );
                     newProxy = oldProxy;
                  }
                  if (subinfo) {
                     DOLOOP {
                        int objectSize = cle->Size();
                        char *obj = arr[k]+ioffset;
                        char *end = obj + compinfo[i]->fLength*objectSize;

                        for(; obj<end; obj+=objectSize) {
                           TVirtualCollectionProxy::TPushPop helper( newProxy, obj );
                           Int_t nobjects;
                           b >> nobjects;
                           void* env = newProxy->Allocate(nobjects,true);
                           subinfo->ReadBufferSTL(b,newProxy,nobjects,/* offset */ 0, vers >= 7);
                           newProxy->Commit(env);
                        }
                     }
                  }
                  b.CheckByteCount(start,count,aElement->GetTypeName());
                  continue;
               }
               if (fOldVersion<3){   // case of old TStreamerInfo
                  //  Backward compatibility. Some TStreamerElement's where without
                  //  Streamer but were not removed from element list
                  if (aElement->IsBase() && aElement->IsA()!=TStreamerBase::Class()) {
                     b.SetBufferOffset(start);  //there is no byte count
                  } else if (vers==0) {
                     b.SetBufferOffset(start);  //there is no byte count
                  }
               }
               if (pstreamer == 0) {
                  if( !newCle ) {
                     newCle = cle;
                     cle = 0;
                  }
                  DOLOOP {
                     b.ReadFastArray((void*)(arr[k]+ioffset),newCle,compinfo[i]->fLength,(TMemberStreamer*)0, cle );
                  }
               } else {
                  DOLOOP {(*pstreamer)(b,arr[k]+ioffset,compinfo[i]->fLength);}
               }
               b.CheckByteCount(start,count,aElement->GetTypeName());
            }
            continue;

         case TStreamerInfo::kObject: // Class derived from TObject
            if (cle->IsStartingWithTObject() && cle->GetState() > TClass::kEmulated) {
               DOLOOP {((TObject*)(arr[k]+ioffset))->Streamer(b);}
               continue; // intentionally inside the if statement.
                      // if the class does not start with its TObject part (or does
                      // not have one), we use the generic case.
            }
         case TStreamerInfo::kAny:    // Class not derived from TObject
            if (pstreamer) {
               DOLOOP {(*pstreamer)(b,arr[k]+ioffset,0);}
            } else {
               if( newCle )
                  DOLOOP { newCle->Streamer( arr[k]+ioffset, b, cle ); }
               else
                  DOLOOP { cle->Streamer(arr[k]+ioffset,b);}
            }
            continue;

         case TStreamerInfo::kObject+TStreamerInfo::kOffsetL:  {
            TFile *file = (TFile*)b.GetParent();
            if (file && file->GetVersion() < 30208) {
               // For older ROOT file we use a totally different case to treat
               // this situation, so we change 'kase' and restart.
               kase = TStreamerInfo::kStreamer;
               goto SWIT;
            }
            // there is intentionally no break/continue statement here.
            // For newer ROOT file, we always use the generic case for kOffsetL(s)
         }

         case TStreamerInfo::kAny+TStreamerInfo::kOffsetL: {
            DOLOOP {
               b.ReadFastArray((void*)(arr[k]+ioffset),cle,compinfo[i]->fLength,pstreamer);
            }
            continue;
         }

         // Base Class
         case TStreamerInfo::kBase:
            if (!(arrayMode&1)) {
               if(pstreamer)  {kase = TStreamerInfo::kStreamer; goto SWIT;}
               DOLOOP { ((TStreamerBase*)aElement)->ReadBuffer(b,arr[k]);}
            } else {
               // FIXME: Rather than relying on the StreamerElement to
               // contain the base class version information we should
               // embed it in the bytestream even in the member-wise case.
               // For now rely, on the StreamerElement:
               TStreamerInfo *binfo = ((TStreamerInfo*)((TStreamerBase*)aElement)->GetBaseStreamerInfo());
               binfo->ReadBuffer(b,arr,binfo->fCompFull,0,binfo->fNfulldata,narr,ioffset,arrayMode);
            }
            continue;

         case TStreamerInfo::kOffsetL + TStreamerInfo::kTString:
         case TStreamerInfo::kOffsetL + TStreamerInfo::kTObject:
         case TStreamerInfo::kOffsetL + TStreamerInfo::kTNamed:
         {
            //  Backward compatibility. Some TStreamerElement's where without
            //  Streamer but were not removed from element list
            UInt_t start,count;
            Version_t v = b.ReadVersion(&start, &count, cle);
            if (fOldVersion<3){   // case of old TStreamerInfo
               if (count<= 0    || v   !=  fOldVersion) {
                  b.SetBufferOffset(start);
                  continue;
               }
            }
            DOLOOP {
               b.ReadFastArray((void*)(arr[k]+ioffset),cle,compinfo[i]->fLength,pstreamer);
            }
            b.CheckByteCount(start,count,aElement->GetFullName());
            continue;
         }


         case TStreamerInfo::kStreamer:{
            //  Backward compatibility. Some TStreamerElement's where without
            //  Streamer but were not removed from element list
            UInt_t start,count;
            Version_t v = b.ReadVersion(&start, &count, cle);
            if (fOldVersion<3){   // case of old TStreamerInfo
               if (aElement->IsBase() && aElement->IsA()!=TStreamerBase::Class()) {
                  b.SetBufferOffset(start);  //it was no byte count
               } else if (kase == TStreamerInfo::kSTL || kase == TStreamerInfo::kSTL+TStreamerInfo::kOffsetL ||
                          count<= 0    || v   !=  fOldVersion) {
                  b.SetBufferOffset(start);
                  continue;
               }
            }
            if (pstreamer == 0) {
               Error("ReadBuffer","Streamer for %s is null\n",aElement->GetName());
               if (gDebug > 0) {
                  aElement->ls(); continue;
               }
            } else {
               DOLOOP {(*pstreamer)(b,arr[k]+ioffset,compinfo[i]->fLength);}
            }
            b.CheckByteCount(start,count,aElement->GetFullName());
         }
         continue;

         case TStreamerInfo::kStreamLoop:
            // -- A pointer to a varying-length array of objects.
            // MyClass* ary; //[n]
            // -- Or a pointer to a varying-length array of pointers to objects.
            // MyClass** ary; //[n]
         case TStreamerInfo::kStreamLoop + TStreamerInfo::kOffsetL:
            // -- An array of pointers to a varying-length array of objects.
            // MyClass* ary[d]; //[n]
            // -- Or an array of pointers to a varying-length array of pointers to objects.
            // MyClass** ary[d]; //[n]
         {
            // Get the class of the data member.
            TClass* cl = compinfo[i]->fClass;
            // Which are we, an array of objects or an array of pointers to objects?
            Bool_t isPtrPtr = (strstr(aElement->GetTypeName(), "**") != 0);
            // Check for a private streamer.
            if (pstreamer) {
               // -- We have a private streamer.
               // Read the class version and byte count from the buffer.
               UInt_t start = 0;
               UInt_t count = 0;
               b.ReadVersion(&start, &count, cl);
               // Loop over the entries in the clones array or the STL container.
               for (Int_t k = 0; k < narr; ++k) {
                  Int_t* counter = (Int_t*) (arr[k] /*entry pointer*/ + eoffset /*entry offset*/ + compinfo[i]->fMethod /*counter offset*/);
                  // And call the private streamer, passing it the buffer, the object, and the counter.
                  (*pstreamer)(b, arr[k] /*entry pointer*/ + ioffset /*object offset*/, *counter);
               }
               b.CheckByteCount(start, count, aElement->GetFullName());
               // We are done, next streamer element.
               continue;
            }
            // At this point we do *not* have a private streamer.
            // Get the version of the file we are reading from.
            TFile* file = (TFile*) b.GetParent();
            // By default assume the file version is the newest.
            Int_t fileVersion = kMaxInt;
            if (file) {
               fileVersion = file->GetVersion();
            }
            // Read the class version and byte count from the buffer.
            UInt_t start = 0;
            UInt_t count = 0;
            b.ReadVersion(&start, &count, cl);
            if (fileVersion > 51508) {
               // -- Newer versions allow polymorphic pointers.
               // Loop over the entries in the clones array or the STL container.
               for (Int_t k = 0; k < narr; ++k) {
                  // Get the counter for the varying length array.
                  Int_t vlen = *((Int_t*) (arr[k] /*entry pointer*/ + eoffset /*entry offset*/ + compinfo[i]->fMethod /*counter offset*/));
                  //Int_t realLen;
                  //b >> realLen;
                  //if (realLen != vlen) {
                  //   fprintf(stderr, "read vlen: %d  realLen: %s\n", vlen, realLen);
                  //}
                  // Get a pointer to the array of pointers.
                  char** pp = (char**) (arr[k] /*entry pointer*/ + ioffset /*object offset*/);
                  if (!pp) {
                     continue;
                  }
                  // Loop over each element of the array of pointers to varying-length arrays.
                  for (Int_t ndx = 0; ndx < compinfo[i]->fLength; ++ndx) {
                     //if (!pp[ndx]) {
                        // -- We do not have a pointer to a varying-length array.
                        //Error("ReadBuffer", "The pointer to element %s::%s type %d (%s) is null\n", thisVar->GetName(), aElement->GetFullName(), compinfo[i]->fType, aElement->GetTypeName());
                        //continue;
                     //}
                     // Delete any memory at pp[ndx].
                     if (!isPtrPtr) {
                        cl->DeleteArray(pp[ndx]);
                        pp[ndx] = 0;
                     } else {
                        // Using vlen is wrong here because it has already
                        // been overwritten with the value needed to read
                        // the current record.  Fixing this will require
                        // doing a pass over the object at the beginning
                        // of the I/O and releasing all the buffer memory
                        // for varying length arrays before we overwrite
                        // the counter values.
                        //
                        // For now we will just leak memory, just as we
                        // have always done in the past.  Fix this.
                        //
                        //char** r = (char**) pp[ndx];
                        //if (r) {
                        //   for (Int_t v = 0; v < vlen; ++v) {
                        //      cl->Destructor(r[v]);
                        //      r[v] = 0;
                        //   }
                        //}
                        delete[] pp[ndx];
                        pp[ndx] = 0;
                     }
                     if (!vlen) {
                        continue;
                     }
                     // Note: We now have pp[ndx] is null.
                     // Allocate memory to read into.
                     if (!isPtrPtr) {
                        // -- We are a varying-length array of objects.
                        // Note: Polymorphism is not allowed here.
                        // Allocate a new array of objects to read into.
                        pp[ndx] = (char*) cl->NewArray(vlen);
                        if (!pp[ndx]) {
                           Error("ReadBuffer", "Memory allocation failed!\n");
                           continue;
                        }
                     } else {
                        // -- We are a varying-length array of pointers to objects.
                        // Note: The object pointers are allowed to be polymorphic.
                        // Allocate a new array of pointers to objects to read into.
                        pp[ndx] = (char*) new char*[vlen];
                        if (!pp[ndx]) {
                           Error("ReadBuffer", "Memory allocation failed!\n");
                           continue;
                        }
                        // And set each pointer to null.
                        memset(pp[ndx], 0, vlen * sizeof(char*)); // This is the right size we really have a char**: pp[ndx] = (char*) new char*[vlen];
                     }
                     if (!isPtrPtr) {
                        // -- We are a varying-length array of objects.
                        b.ReadFastArray(pp[ndx], cl, vlen, 0);
                     }
                     else {
                        // -- We are a varying-length array of object pointers.
                        b.ReadFastArray((void**) pp[ndx], cl, vlen, kFALSE, 0);
                     } // isPtrPtr
                  } // ndx
               } // k
            }
            else {
               // -- Older versions do *not* allow polymorphic pointers.
               // Loop over the entries in the clones array or the STL container.
               for (Int_t k = 0; k < narr; ++k) {
                  // Get the counter for the varying length array.
                  Int_t vlen = *((Int_t*) (arr[k] /*entry pointer*/ + eoffset /*entry offset*/ + compinfo[i]->fMethod /*counter offset*/));
                  //Int_t realLen;
                  //b >> realLen;
                  //if (realLen != vlen) {
                  //   fprintf(stderr, "read vlen: %d  realLen: %s\n", vlen, realLen);
                  //}
                  // Get a pointer to the array of pointers.
                  char** pp = (char**) (arr[k] /*entry pointer*/ + ioffset /*object offset*/);
                  if (!pp) {
                     continue;
                  }
                  // Loop over each element of the array of pointers to varying-length arrays.
                  for (Int_t ndx = 0; ndx < compinfo[i]->fLength; ++ndx) {
                     //if (!pp[ndx]) {
                        // -- We do not have a pointer to a varying-length array.
                        //Error("ReadBuffer", "The pointer to element %s::%s type %d (%s) is null\n", thisVar->GetName(), aElement->GetFullName(), compinfo[i]->fType, aElement->GetTypeName());
                        //continue;
                     //}
                     // Delete any memory at pp[ndx].
                     if (!isPtrPtr) {
                        cl->DeleteArray(pp[ndx]);
                        pp[ndx] = 0;
                     } else {
                        // Using vlen is wrong here because it has already
                        // been overwritten with the value needed to read
                        // the current record.  Fixing this will require
                        // doing a pass over the object at the beginning
                        // of the I/O and releasing all the buffer memory
                        // for varying length arrays before we overwrite
                        // the counter values.
                        //
                        // For now we will just leak memory, just as we
                        // have always done in the past.  Fix this.
                        //
                        //char** r = (char**) pp[ndx];
                        //if (r) {
                        //   for (Int_t v = 0; v < vlen; ++v) {
                        //      cl->Destructor(r[v]);
                        //      r[v] = 0;
                        //   }
                        //}
                        delete[] pp[ndx];
                        pp[ndx] = 0;
                     }
                     if (!vlen) {
                        continue;
                     }
                     // Note: We now have pp[ndx] is null.
                     // Allocate memory to read into.
                     if (!isPtrPtr) {
                        // -- We are a varying-length array of objects.
                        // Note: Polymorphism is not allowed here.
                        // Allocate a new array of objects to read into.
                        pp[ndx] = (char*) cl->NewArray(vlen);
                        if (!pp[ndx]) {
                           Error("ReadBuffer", "Memory allocation failed!\n");
                           continue;
                        }
                     } else {
                        // -- We are a varying-length array of pointers to objects.
                        // Note: The object pointers are allowed to be polymorphic.
                        // Allocate a new array of pointers to objects to read into.
                        pp[ndx] = (char*) new char*[vlen];
                        if (!pp[ndx]) {
                           Error("ReadBuffer", "Memory allocation failed!\n");
                           continue;
                        }
                        // And set each pointer to null.
                        memset(pp[ndx], 0, vlen * sizeof(char*)); // This is the right size we really have a char**: pp[ndx] = (char*) new char*[vlen];
                     }
                     if (!isPtrPtr) {
                        // -- We are a varying-length array of objects.
                        // Loop over the elements of the varying length array.
                        for (Int_t v = 0; v < vlen; ++v) {
                           // Read the object from the buffer.
                           cl->Streamer(pp[ndx] + (v * cl->Size()), b);
                        } // v
                     }
                     else {
                        // -- We are a varying-length array of object pointers.
                        // Get a pointer to the object pointer array.
                        char** r = (char**) pp[ndx];
                        // Loop over the elements of the varying length array.
                        for (Int_t v = 0; v < vlen; ++v) {
                           // Allocate an object to read into.
                           r[v] = (char*) cl->New();
                           if (!r[v]) {
                              // Do not print a second error message here.
                              //Error("ReadBuffer", "Memory allocation failed!\n");
                              continue;
                           }
                           // Read the object from the buffer.
                           cl->Streamer(r[v], b);
                        } // v
                     } // isPtrPtr
                  } // ndx
               } // k
            } // fileVersion
            b.CheckByteCount(start, count, aElement->GetFullName());
            continue;
         }

         case TStreamerInfo::kCacheNew:
            ((TBufferFile&)b).PushDataCache( new TVirtualArray( aElement->GetClassPointer(), narr ) );
            continue;
         case TStreamerInfo::kCacheDelete:
            delete ((TBufferFile&)b).PopDataCache();
            continue;

         case -1:
            // -- Skip an ignored TObject base class.
            continue;

         default: {
            int ans = -1;

            if (TStreamerInfo::kCache <= kase && kase < TStreamerInfo::kArtificial) {

               //T &cache_add = *(T*)b.PeekDataCacheArray();
               R__ASSERT(kFALSE); // cache_add);

               // thisVar->ReadBuffer(b,cache_addr,i,kase-TStreamerInfo::kCache,aElement,narr,eoffset)

               continue;
            }

            if (kase >= TStreamerInfo::kConv)
               ans = thisVar->ReadBufferConv(b,arr,compinfo[i],kase,aElement,narr,eoffset);
            if (ans==0) continue;

            if (kase >= TStreamerInfo::kSkip)
               ans = thisVar->ReadBufferSkip(b,arr,compinfo[i],kase,aElement,narr,eoffset);
            if (ans==0) continue;

            if (kase >= TStreamerInfo::kArtificial) {
               ans = thisVar->ReadBufferArtificial(b,arr,aElement,narr,eoffset);
            }
            if (ans==0) continue;
         }
         if (aElement)
            Error("ReadBuffer","The element %s::%s type %d (%s) is not supported yet\n",
                  thisVar->GetName(),aElement->GetFullName(),kase,aElement->GetTypeName());
         else
            Error("ReadBuffer","The TStreamerElement for %s %d is missing!\n",
                  thisVar->GetName(),i);

         continue;
      }
   }
   if (needIncrement) b.DecrementLevel(thisVar);
   return 0;
}

template Int_t TStreamerInfo::ReadBufferSkip<char**>(TBuffer &b, char** const &arr, const TCompInfo *compinfo, Int_t kase,
                                       TStreamerElement *aElement, Int_t narr,
                                       Int_t eoffset);
template Int_t TStreamerInfo::ReadBufferSkip<TVirtualCollectionProxy>(TBuffer &b, const TVirtualCollectionProxy &arr, const TCompInfo *compinfo, Int_t kase,
                                       TStreamerElement *aElement, Int_t narr,
                                       Int_t eoffset);
template Int_t TStreamerInfo::ReadBufferSkip<TVirtualArray>(TBuffer &b, const TVirtualArray &arr, const TCompInfo *compinfo, Int_t kase,
                                       TStreamerElement *aElement, Int_t narr,
                                       Int_t eoffset);

template Int_t TStreamerInfo::ReadBufferConv<char**>(TBuffer &b, char** const &arr,  const TCompInfo *compinfo, Int_t kase,
                                       TStreamerElement *aElement, Int_t narr,
                                       Int_t eoffset);
template Int_t TStreamerInfo::ReadBufferConv<TVirtualCollectionProxy>(TBuffer &b, const TVirtualCollectionProxy &arr,  const TCompInfo *compinfo, Int_t kase,
                                       TStreamerElement *aElement, Int_t narr,
                                       Int_t eoffset);
template Int_t TStreamerInfo::ReadBufferConv<TVirtualArray>(TBuffer &b, const TVirtualArray &arr, const TCompInfo *compinfo, Int_t kase,
                                       TStreamerElement *aElement, Int_t narr,
                                       Int_t eoffset);

template Int_t TStreamerInfo::ReadBufferArtificial<char**>(TBuffer &b, char** const &arr,
                                             TStreamerElement *aElement, Int_t narr,
                                             Int_t eoffset);
template Int_t TStreamerInfo::ReadBufferArtificial<TVirtualCollectionProxy>(TBuffer &b, const TVirtualCollectionProxy &arr,
                                             TStreamerElement *aElement, Int_t narr,
                                             Int_t eoffset);
template Int_t TStreamerInfo::ReadBufferArtificial<TVirtualArray>(TBuffer &b, const TVirtualArray &arr,
                                             TStreamerElement *aElement, Int_t narr,
                                             Int_t eoffset);

template Int_t TStreamerInfo::ReadBuffer<char**>(TBuffer &b, char** const &arr,
                                                 TCompInfo *const*const compinfo, Int_t first, Int_t last,
                                                 Int_t narr, Int_t eoffset, Int_t arrayMode);
template Int_t TStreamerInfo::ReadBuffer<TVirtualCollectionProxy>(TBuffer &b, const TVirtualCollectionProxy &arr,
                                                                  TCompInfo *const*const compinfo, Int_t first, Int_t last,
                                                                  Int_t narr, Int_t eoffset, Int_t arrayMode);
template Int_t TStreamerInfo::ReadBuffer<TVirtualArray>(TBuffer &b, const TVirtualArray &arr,
                                                        TCompInfo *const*const compinfo, Int_t first, Int_t last,
                                                        Int_t narr, Int_t eoffset, Int_t arrayMode);

////////////////////////////////////////////////////////////////////////////////
///  The STL vector/list is deserialized from the buffer b

Int_t TStreamerInfo::ReadBufferSTL(TBuffer &b, TVirtualCollectionProxy *cont,
                                   Int_t nc, Int_t eoffset, Bool_t v7 /* = kTRUE */)
{
   if (!nc && v7) return 0; // in version 6 of TStreamerInfo and below, we were calling ReadBuffer for empty collection.
   int ret = ReadBuffer(b, *cont,fCompFull,0,fNfulldata,nc,eoffset,1);
   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// Read for TClonesArray.
/// Note: This is no longer used.

Int_t TStreamerInfo::ReadBufferClones(TBuffer &b, TClonesArray *clones,
                                      Int_t nc, Int_t first, Int_t eoffset)
{
   char **arr = (char **)clones->GetObjectRef(0);
   return ReadBuffer(b,arr,fCompFull,first==-1?0:first,first==-1?fNfulldata:first+1,nc,eoffset,1);
}
