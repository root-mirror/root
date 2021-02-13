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
#include "TClass.h"
#include "TClonesArray.h"
#include "TError.h"
#include "TProcessID.h"
#include "TStreamer.h"
#include "TStreamerElement.h"
#include "TStreamerInfo.h"
#include "TVirtualCollectionProxy.h"
#include "TRefTable.h"
#include "TFile.h"

#include "TVirtualArray.h"
#include "TBufferFile.h"
#include "TInterpreter.h"

//==========CPP macros

#define DOLOOP for(int k=0; k<narr; ++k)

#define WriteBasicTypeElem(name,index)          \
   {                                            \
      name *x=(name*)(arr[index]+ioffset);      \
      b << *x;                                  \
   }

#define WriteBasicType(name)                    \
   {                                            \
      WriteBasicTypeElem(name,0);               \
   }

#define WriteBasicTypeLoop(name)                            \
   {                                                        \
      for(int k=0; k<narr; ++k) WriteBasicTypeElem(name,k); \
   }

#define WriteBasicArrayElem(name,index)         \
   {                                            \
      name *x=(name*)(arr[index]+ioffset);      \
      b.WriteFastArray(x,compinfo[i]->fLength);           \
   }

#define WriteBasicArray(name)                   \
   {                                            \
      WriteBasicArrayElem(name,0);              \
   }

#define WriteBasicArrayLoop(name)                              \
   {                                                           \
      for(int k=0; k<narr; ++k) WriteBasicArrayElem(name,k);   \
   }

#define WriteBasicPointerElem(name,index)       \
   {                                            \
      Int_t *l = (Int_t*)(arr[index]+imethod);  \
      name **f = (name**)(arr[index]+ioffset);  \
      name *af = *f;                            \
      if (af && *l)  b << Char_t(1);            \
      else          {b << Char_t(0); continue;} \
      int j;                                    \
      for(j=0;j<compinfo[i]->fLength;j++) {               \
         b.WriteFastArray(f[j],*l);             \
      }                                         \
   }

#define WriteBasicPointer(name)                 \
   {                                            \
      int imethod = compinfo[i]->fMethod+eoffset;         \
      WriteBasicPointerElem(name,0);            \
   }

#define WriteBasicPointerLoop(name)             \
   {                                            \
      int imethod = compinfo[i]->fMethod+eoffset;         \
      for(int k=0; k<narr; ++k) {               \
         WriteBasicPointerElem(name,k);         \
      }                                         \
   }

// Helper function for TStreamerInfo::WriteBuffer
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
///  The object at pointer is serialized to the buffer b
///  if (arrayMode & 1) ptr is a pointer to array of pointers to the objects
///  otherwise it is a pointer to a pointer to a single object.
///  This also means that T is of a type such that arr[i] is a pointer to an
///  object.  Currently the only anticipated instantiation are for T==char**
///  and T==TVirtualCollectionProxy

template <class T>
Int_t TStreamerInfo::WriteBufferAux(TBuffer &b, const T &arr,
                                    TCompInfo *const*const compinfo, Int_t first, Int_t last,
                                    Int_t narr, Int_t eoffset, Int_t arrayMode)
{
   Bool_t needIncrement = !( arrayMode & 2 );
   arrayMode = arrayMode & (~2);

   if (needIncrement) b.IncrementLevel(this);

   //mark this class as being used in the current file
   b.TagStreamerInfo(this);

   //============

   //loop on all active members
//   Int_t last;
//   if (first < 0) {first = 0; last = ninfo;}
//   else            last = first+1;

   // In order to speed up the case where the object being written is
   // not in a collection (i.e. arrayMode is false), we actually
   // duplicate the elementary types using this typeOffset.
   static const int kHaveLoop = 1024;
   const Int_t typeOffset = arrayMode ? kHaveLoop : 0;

   for (Int_t i=first;i<last;i++) {

      TStreamerElement *aElement = (TStreamerElement*)compinfo[i]->fElem;

      if (needIncrement) b.SetStreamerElementNumber(aElement,compinfo[i]->fType);

      Int_t ioffset = eoffset+compinfo[i]->fOffset;

      if (R__TestUseCache<T>(aElement)) {
         if (aElement->TestBit(TStreamerElement::kWrite)) {
            if (((TBufferFile&)b).PeekDataCache()==0) {
               Warning("WriteBuffer","Skipping %s::%s because the cache is missing.",GetName(),aElement->GetName());
            } else {
               if (gDebug > 1) {
                  printf("WriteBuffer, class:%s, name=%s, fType[%d]=%d,"
                         " %s, bufpos=%d, arr=%p, eoffset=%d, Redirect=%p\n",
                         fClass->GetName(),aElement->GetName(),i,compinfo[i]->fType,
                         aElement->ClassName(),b.Length(),arr[0], eoffset,((TBufferFile&)b).PeekDataCache()->GetObjectAt(0));
               }
               WriteBufferAux(b,*((TBufferFile&)b).PeekDataCache(),compinfo,i,i+1,narr,eoffset, arrayMode);
            }
            continue;
         } else {
            if (gDebug > 1) {
               printf("WriteBuffer, class:%s, name=%s, fType[%d]=%d,"
                      " %s, bufpos=%d, arr=%p, eoffset=%d, not a write rule, skipping.\n",
                      fClass->GetName(),aElement->GetName(),i,compinfo[i]->fType,
                      aElement->ClassName(),b.Length(),arr[0], eoffset);
            }
            // The rule was a cached element for a read, rule, the real offset is in the
            // next element (the one for the rule itself).
            if (aElement->TestBit(TStreamerElement::kRepeat)) continue;
            ioffset = eoffset+compinfo[i]->fOffset;
            continue;
         }
      }

      if (gDebug > 1) {
         printf("WriteBuffer, class:%s, name=%s, fType[%d]=%d, %s, "
               "bufpos=%d, arr=%p, offset=%d\n",
                fClass->GetName(),aElement->GetName(),i,compinfo[i]->fType,aElement->ClassName(),
                b.Length(),arr[0],ioffset);
      }

      switch (compinfo[i]->fType+typeOffset) {
         // In this switch we intentionally use 'continue' instead of
         // 'break' to avoid running the 2nd switch (see later in this
         // function).

         case TStreamerInfo::kBool:                WriteBasicType(Bool_t);    continue;
         case TStreamerInfo::kChar:                WriteBasicType(Char_t);    continue;
         case TStreamerInfo::kShort:               WriteBasicType(Short_t);   continue;
         case TStreamerInfo::kInt:                 WriteBasicType(Int_t);     continue;
         case TStreamerInfo::kLong:                WriteBasicType(Long_t);    continue;
         case TStreamerInfo::kLong64:              WriteBasicType(Long64_t);  continue;
         case TStreamerInfo::kFloat:               WriteBasicType(Float_t);   continue;
         case TStreamerInfo::kDouble:              WriteBasicType(Double_t);  continue;
         case TStreamerInfo::kUChar:               WriteBasicType(UChar_t);   continue;
         case TStreamerInfo::kUShort:              WriteBasicType(UShort_t);  continue;
         case TStreamerInfo::kUInt:                WriteBasicType(UInt_t);    continue;
         case TStreamerInfo::kULong:               WriteBasicType(ULong_t);   continue;
         case TStreamerInfo::kULong64:             WriteBasicType(ULong64_t); continue;
         case TStreamerInfo::kFloat16: {
            Float_t *x=(Float_t*)(arr[0]+ioffset);
            b.WriteFloat16(x,aElement);
            continue;
         }
         case TStreamerInfo::kDouble32: {
            Double_t *x=(Double_t*)(arr[0]+ioffset);
            b.WriteDouble32(x,aElement);
            continue;
         }

         case TStreamerInfo::kBool    + kHaveLoop: WriteBasicTypeLoop(Bool_t);    continue;
         case TStreamerInfo::kChar    + kHaveLoop: WriteBasicTypeLoop(Char_t);    continue;
         case TStreamerInfo::kShort   + kHaveLoop: WriteBasicTypeLoop(Short_t);   continue;
         case TStreamerInfo::kInt     + kHaveLoop: WriteBasicTypeLoop(Int_t);     continue;
         case TStreamerInfo::kLong    + kHaveLoop: WriteBasicTypeLoop(Long_t);    continue;
         case TStreamerInfo::kLong64  + kHaveLoop: WriteBasicTypeLoop(Long64_t);  continue;
         case TStreamerInfo::kFloat   + kHaveLoop: WriteBasicTypeLoop(Float_t);   continue;
         case TStreamerInfo::kDouble  + kHaveLoop: WriteBasicTypeLoop(Double_t);  continue;
         case TStreamerInfo::kUChar   + kHaveLoop: WriteBasicTypeLoop(UChar_t);   continue;
         case TStreamerInfo::kUShort  + kHaveLoop: WriteBasicTypeLoop(UShort_t);  continue;
         case TStreamerInfo::kUInt    + kHaveLoop: WriteBasicTypeLoop(UInt_t);    continue;
         case TStreamerInfo::kULong   + kHaveLoop: WriteBasicTypeLoop(ULong_t);   continue;
         case TStreamerInfo::kULong64 + kHaveLoop: WriteBasicTypeLoop(ULong64_t); continue;
         case TStreamerInfo::kFloat16+ kHaveLoop: {
            for(int k=0; k<narr; ++k) {
               Float_t *x=(Float_t*)(arr[k]+ioffset);
               b.WriteFloat16(x,aElement);
            }
            continue;
         }
         case TStreamerInfo::kDouble32+ kHaveLoop: {
            for(int k=0; k<narr; ++k) {
               Double_t *x=(Double_t*)(arr[k]+ioffset);
               b.WriteDouble32(x,aElement);
            }
            continue;
         }

         // write array of basic types  array[8]
         case TStreamerInfo::kOffsetL + TStreamerInfo::kBool:   WriteBasicArray(Bool_t);    continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kChar:   WriteBasicArray(Char_t);    continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kShort:  WriteBasicArray(Short_t);   continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kInt:    WriteBasicArray(Int_t);     continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kLong:   WriteBasicArray(Long_t);    continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kLong64: WriteBasicArray(Long64_t);  continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kFloat:  WriteBasicArray(Float_t);   continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kDouble: WriteBasicArray(Double_t);  continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kUChar:  WriteBasicArray(UChar_t);   continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kUShort: WriteBasicArray(UShort_t);  continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kUInt:   WriteBasicArray(UInt_t);    continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kULong:  WriteBasicArray(ULong_t);   continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kULong64:WriteBasicArray(ULong64_t); continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kFloat16: {
            b.WriteFastArrayFloat16((Float_t*)(arr[0]+ioffset),compinfo[i]->fLength,aElement);
            continue;
         }
         case TStreamerInfo::kOffsetL + TStreamerInfo::kDouble32: {
            b.WriteFastArrayDouble32((Double_t*)(arr[0]+ioffset),compinfo[i]->fLength,aElement);
            continue;
         }

         case TStreamerInfo::kOffsetL + TStreamerInfo::kBool    + kHaveLoop: WriteBasicArrayLoop(Bool_t);    continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kChar    + kHaveLoop: WriteBasicArrayLoop(Char_t);    continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kShort   + kHaveLoop: WriteBasicArrayLoop(Short_t);   continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kInt     + kHaveLoop: WriteBasicArrayLoop(Int_t);     continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kLong    + kHaveLoop: WriteBasicArrayLoop(Long_t);    continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kLong64  + kHaveLoop: WriteBasicArrayLoop(Long64_t);  continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kFloat   + kHaveLoop: WriteBasicArrayLoop(Float_t);   continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kDouble  + kHaveLoop: WriteBasicArrayLoop(Double_t);  continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kUChar   + kHaveLoop: WriteBasicArrayLoop(UChar_t);   continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kUShort  + kHaveLoop: WriteBasicArrayLoop(UShort_t);  continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kUInt    + kHaveLoop: WriteBasicArrayLoop(UInt_t);    continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kULong   + kHaveLoop: WriteBasicArrayLoop(ULong_t);   continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kULong64 + kHaveLoop: WriteBasicArrayLoop(ULong64_t); continue;
         case TStreamerInfo::kOffsetL + TStreamerInfo::kFloat16+ kHaveLoop: {
            for(int k=0; k<narr; ++k) {
               b.WriteFastArrayFloat16((Float_t*)(arr[k]+ioffset),compinfo[i]->fLength,aElement);
            }
            continue;
         }
         case TStreamerInfo::kOffsetL + TStreamerInfo::kDouble32+ kHaveLoop: {
            for(int k=0; k<narr; ++k) {
               b.WriteFastArrayDouble32((Double_t*)(arr[k]+ioffset),compinfo[i]->fLength,aElement);
            }
            continue;
         }

         // write pointer to an array of basic types  array[n]
         case TStreamerInfo::kOffsetP + TStreamerInfo::kBool:   WriteBasicPointer(Bool_t);    continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kChar:   WriteBasicPointer(Char_t);    continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kShort:  WriteBasicPointer(Short_t);   continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kInt:    WriteBasicPointer(Int_t);     continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kLong:   WriteBasicPointer(Long_t);    continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kLong64: WriteBasicPointer(Long64_t);  continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kFloat:  WriteBasicPointer(Float_t);   continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kDouble: WriteBasicPointer(Double_t);  continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kUChar:  WriteBasicPointer(UChar_t);   continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kUShort: WriteBasicPointer(UShort_t);  continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kUInt:   WriteBasicPointer(UInt_t);    continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kULong:  WriteBasicPointer(ULong_t);   continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kULong64:WriteBasicPointer(ULong64_t); continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kFloat16: {
            int imethod = compinfo[i]->fMethod+eoffset;
            Int_t *l = (Int_t*)(arr[0]+imethod);
            Float_t **f = (Float_t**)(arr[0]+ioffset);
            Float_t *af = *f;
            if (af && *l)  b << Char_t(1);
            else          {b << Char_t(0); continue;}
            int j;
            for(j=0;j<compinfo[i]->fLength;j++) {
               b.WriteFastArrayFloat16(f[j],*l,aElement);
            }
            continue;
         }
         case TStreamerInfo::kOffsetP + TStreamerInfo::kDouble32: {
            int imethod = compinfo[i]->fMethod+eoffset;
            Int_t *l = (Int_t*)(arr[0]+imethod);
            Double_t **f = (Double_t**)(arr[0]+ioffset);
            Double_t *af = *f;
            if (af && *l)  b << Char_t(1);
            else          {b << Char_t(0); continue;}
            int j;
            for(j=0;j<compinfo[i]->fLength;j++) {
               b.WriteFastArrayDouble32(f[j],*l,aElement);
            }
            continue;
         }

         case TStreamerInfo::kOffsetP + TStreamerInfo::kBool    + kHaveLoop: WriteBasicPointerLoop(Bool_t);    continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kChar    + kHaveLoop: WriteBasicPointerLoop(Char_t);    continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kShort   + kHaveLoop: WriteBasicPointerLoop(Short_t);   continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kInt     + kHaveLoop: WriteBasicPointerLoop(Int_t);     continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kLong    + kHaveLoop: WriteBasicPointerLoop(Long_t);    continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kLong64  + kHaveLoop: WriteBasicPointerLoop(Long64_t);  continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kFloat   + kHaveLoop: WriteBasicPointerLoop(Float_t);   continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kDouble  + kHaveLoop: WriteBasicPointerLoop(Double_t);  continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kUChar   + kHaveLoop: WriteBasicPointerLoop(UChar_t);   continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kUShort  + kHaveLoop: WriteBasicPointerLoop(UShort_t);  continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kUInt    + kHaveLoop: WriteBasicPointerLoop(UInt_t);    continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kULong   + kHaveLoop: WriteBasicPointerLoop(ULong_t);   continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kULong64 + kHaveLoop: WriteBasicPointerLoop(ULong64_t); continue;
         case TStreamerInfo::kOffsetP + TStreamerInfo::kFloat16+ kHaveLoop: {
            int imethod = compinfo[i]->fMethod+eoffset;
            for(int k=0; k<narr; ++k) {
               Int_t *l = (Int_t*)(arr[k]+imethod);
               Float_t **f = (Float_t**)(arr[k]+ioffset);
               Float_t *af = *f;
               if (af && *l)  b << Char_t(1);
               else          {b << Char_t(0); continue;}
               int j;
               for(j=0;j<compinfo[i]->fLength;j++) {
                  b.WriteFastArrayFloat16(f[j],*l,aElement);
               }
            }
            continue;
         }
         case TStreamerInfo::kOffsetP + TStreamerInfo::kDouble32+ kHaveLoop: {
            int imethod = compinfo[i]->fMethod+eoffset;
            for(int k=0; k<narr; ++k) {
               Int_t *l = (Int_t*)(arr[k]+imethod);
               Double_t **f = (Double_t**)(arr[k]+ioffset);
               Double_t *af = *f;
               if (af && *l)  b << Char_t(1);
               else          {b << Char_t(0); continue;}
               int j;
               for(j=0;j<compinfo[i]->fLength;j++) {
                  b.WriteFastArrayDouble32(f[j],*l,aElement);
               }
            }
            continue;
         }

         case TStreamerInfo::kCounter: {
            Int_t *x=(Int_t*)(arr[0]+ioffset);
            b << *x;
            if (i == last-1) {
               if (needIncrement) b.DecrementLevel(this);
               return x[0]; // info used by TBranchElement::FillLeaves
            }
            continue;
         }

         case TStreamerInfo::kCounter + kHaveLoop : {
            DOLOOP {
               Int_t *x=(Int_t*)(arr[k]+ioffset);
               b << x[0];
            }
            continue;
         }


      };
      Bool_t isPreAlloc = 0;

      switch (compinfo[i]->fType) {

         // char*
         case TStreamerInfo::kCharStar: { DOLOOP {
            char **f = (char**)(arr[k]+ioffset);
            b.WriteCharStar(*f);
         }
         continue; }

         // special case for TObject::fBits in case of a referenced object
         case TStreamerInfo::kBits: { DOLOOP {
            UInt_t *x=(UInt_t*)(arr[k]+ioffset); b << *x;
            if ((*x & kIsReferenced) != 0) {
               TObject *obj = (TObject*)(arr[k]+eoffset);
               TProcessID *pid = TProcessID::GetProcessWithUID(obj->GetUniqueID(),obj);
               TRefTable *table = TRefTable::GetRefTable();
               if(table) table->Add(obj->GetUniqueID(),pid);
               UShort_t pidf = b.WriteProcessID(pid);
               b << pidf;
            }
         }
         continue; }

         // Special case for TString, TObject, TNamed
         case TStreamerInfo::kTString: { DOLOOP{ ((TString*)(arr[k]+ioffset))->Streamer(b);         }; continue; }
         case TStreamerInfo::kTObject: { DOLOOP{ ((TObject*)(arr[k]+ioffset))->TObject::Streamer(b);}; continue; }
         case TStreamerInfo::kTNamed:  { DOLOOP{ ((TNamed *)(arr[k]+ioffset))->TNamed::Streamer(b); }; continue; }

         case TStreamerInfo::kAnyp:     // Class*   Class not derived from TObject and with comment field //->
         case TStreamerInfo::kAnyp    + TStreamerInfo::kOffsetL:

         case TStreamerInfo::kObjectp:  // Class *  Class     derived from TObject and with comment field //->
         case TStreamerInfo::kObjectp + TStreamerInfo::kOffsetL:

            isPreAlloc = kTRUE;

            // Intentional fallthrough now that isPreAlloc is set.
         case TStreamerInfo::kAnyP:         // Class*   Class not derived from TObject and no comment
         case TStreamerInfo::kAnyP    + TStreamerInfo::kOffsetL:

         case TStreamerInfo::kObjectP:  // Class*   Class derived from TObject
         case TStreamerInfo::kObjectP + TStreamerInfo::kOffsetL: {
            TClass *cl                 = compinfo[i]->fClass;
            TMemberStreamer *pstreamer = compinfo[i]->fStreamer;
            DOLOOP {
               Int_t res = b.WriteFastArray((void**)(arr[k]+ioffset),cl,compinfo[i]->fLength,isPreAlloc,pstreamer);
               if (res==2) {
                  Warning("WriteBuffer",
                          "The actual class of %s::%s is not available. Only the \"%s\" part will be written\n",
                          GetName(),aElement->GetName(),cl->GetName());
               }
            }
            continue;
         }

         case TStreamerInfo::kAnyPnoVT:     // Class*   Class not derived from TObject and no virtual table and no comment
         case TStreamerInfo::kAnyPnoVT    + TStreamerInfo::kOffsetL: {
            TClass *cl                 = compinfo[i]->fClass;
            TMemberStreamer *pstreamer = compinfo[i]->fStreamer;
            DOLOOP {
               void **f = (void**)(arr[k]+ioffset);
               int j;
               for(j=0;j<compinfo[i]->fLength;j++) {
                  void *af = f[j];
                  if (af)  b << Char_t(1);
                  else    {b << Char_t(0); continue;}
                  if (pstreamer) (*pstreamer)(b, af, 0);
                  else cl->Streamer( af, b );
               }
            }
            continue;
         }

//          case TStreamerInfo::kSTLvarp:           // Variable size array of STL containers.
//             {
//                TMemberStreamer *pstreamer = compinfo[i]->fStreamer;
//                TClass *cl                 = compinfo[i]->fClass;
//                UInt_t pos = b.WriteVersion(this->IsA(),kTRUE);
//                if (pstreamer == 0) {
//                   Int_t size = cl->Size();
//                   Int_t imethod = compinfo[i]->fMethod+eoffset;
//                   DOLOOP {
//                      char **contp = (char**)(arr[k]+ioffset);
//                      const Int_t *counter = (Int_t*)(arr[k]+imethod);
//                      const Int_t sublen = (*counter);

//                      for(int j=0;j<compinfo[i]->fLength;++j) {
//                         char *cont = contp[j];
//                         for(int k=0;k<sublen;++k) {
//                            cl->Streamer( cont, b );
//                            cont += size;
//                         }
//                      }
//                   }
//                } else {
//                   DOLOOP{(*pstreamer)(b,arr[k]+ioffset,compinfo[i]->fLength);}
//                }
//                b.SetByteCount(pos,kTRUE);
//             }
//             continue;


         case TStreamerInfo::kSTLp:                // Pointer to container with no virtual table (stl) and no comment
         case TStreamerInfo::kSTLp + TStreamerInfo::kOffsetL:     // array of pointers to container with no virtual table (stl) and no comment
            {
               TClass *cl                 = compinfo[i]->fClass;
               TMemberStreamer *pstreamer = compinfo[i]->fStreamer;
               TVirtualCollectionProxy *proxy = cl->GetCollectionProxy();
               TClass* vClass = proxy ? proxy->GetValueClass() : 0;

               if (!b.TestBit(TBuffer::kCannotHandleMemberWiseStreaming)
                   && proxy && vClass
                   && GetStreamMemberWise()
                   && cl->CanSplit()
                   && !(strspn(aElement->GetTitle(),"||") == 2)
                   && !(vClass->TestBit(TClass::kHasCustomStreamerMember)) ) {
                  // Let's save the collection member-wise.

                  UInt_t pos = b.WriteVersionMemberWise(this->IsA(),kTRUE);
                  b.WriteVersion( vClass, kFALSE );
                  DOLOOP {
                     char **contp = (char**)(arr[k]+ioffset);
                     for(int j=0;j<compinfo[i]->fLength;++j) {
                        char *cont = contp[j];
                        TVirtualCollectionProxy::TPushPop helper( proxy, cont );
                        Int_t nobjects = cont ? proxy->Size() : 0;
                        b << nobjects;
                        if (nobjects) {
                           auto actions = proxy->GetWriteMemberWiseActions();

                           char startbuf[TVirtualCollectionProxy::fgIteratorArenaSize];
                           char endbuf[TVirtualCollectionProxy::fgIteratorArenaSize];
                           void *begin = &(startbuf[0]);
                           void *end = &(endbuf[0]);
                           proxy->GetFunctionCreateIterators()(cont, &begin, &end, proxy);
                           // We can not get here with a split vector of pointer, so we can indeed assume
                           // that actions->fConfiguration != null.
                           b.ApplySequence(*actions, begin, end);
                           if (begin != &(startbuf[0])) {
                              // assert(end != endbuf);
                              proxy->GetFunctionDeleteTwoIterators()(begin,end);
                           }
                        }
                     }
                  }
                  b.SetByteCount(pos,kTRUE);
                  continue;
               }
               UInt_t pos = b.WriteVersion(this->IsA(),kTRUE);
               if (pstreamer == 0) {
                  DOLOOP {
                     char **contp = (char**)(arr[k]+ioffset);
                     for(int j=0;j<compinfo[i]->fLength;++j) {
                        char *cont = contp[j];
                        cl->Streamer( cont, b );
                     }
                  }
               } else {
                  DOLOOP{(*pstreamer)(b,arr[k]+ioffset,compinfo[i]->fLength);}
               }
               b.SetByteCount(pos,kTRUE);
            }
            continue;

         case TStreamerInfo::kSTL:             // container with no virtual table (stl) and no comment
         case TStreamerInfo::kSTL + TStreamerInfo::kOffsetL:  // array of containers with no virtual table (stl) and no comment
            {
               TClass *cl                 = compinfo[i]->fClass;
               TMemberStreamer *pstreamer = compinfo[i]->fStreamer;
               TVirtualCollectionProxy *proxy = cl->GetCollectionProxy();
               TClass* vClass = proxy ? proxy->GetValueClass() : 0;
               if (!b.TestBit(TBuffer::kCannotHandleMemberWiseStreaming)
                   && proxy && vClass
                   && GetStreamMemberWise() && cl->CanSplit()
                   && !(strspn(aElement->GetTitle(),"||") == 2)
                   && !(vClass->TestBit(TClass::kHasCustomStreamerMember)) ) {
                  // Let's save the collection in member-wise order.

                  UInt_t pos = b.WriteVersionMemberWise(this->IsA(),kTRUE);
                  b.WriteVersion( vClass, kFALSE );
                  DOLOOP {
                     char *obj = (char*)(arr[k]+ioffset);
                     Int_t n = compinfo[i]->fLength;
                     if (!n) n=1;
                     int size = cl->Size();

                     for(Int_t j=0; j<n; j++,obj+=size) {
                        TVirtualCollectionProxy::TPushPop helper( proxy, obj );
                        Int_t nobjects = proxy->Size();
                        b << nobjects;
                        if (nobjects) {
                           auto actions = proxy->GetWriteMemberWiseActions();

                           char startbuf[TVirtualCollectionProxy::fgIteratorArenaSize];
                           char endbuf[TVirtualCollectionProxy::fgIteratorArenaSize];
                           void *begin = &(startbuf[0]);
                           void *end = &(endbuf[0]);
                           proxy->GetFunctionCreateIterators()(obj, &begin, &end, proxy);
                           // We can not get here with a split vector of pointer, so we can indeed assume
                           // that actions->fConfiguration != null.
                           b.ApplySequence(*actions, begin, end);
                           if (begin != &(startbuf[0])) {
                              // assert(end != endbuf);
                              proxy->GetFunctionDeleteTwoIterators()(begin,end);
                           }
                        }
                     }
                  }
                  b.SetByteCount(pos,kTRUE);
                  continue;
               }
               UInt_t pos = b.WriteVersion(this->IsA(),kTRUE);
               if (pstreamer == 0) {
                  DOLOOP {
                     b.WriteFastArray((void*)(arr[k]+ioffset),cl,compinfo[i]->fLength,0);
                  }
               } else {
                  DOLOOP{(*pstreamer)(b,arr[k]+ioffset,compinfo[i]->fLength);}
               }
               b.SetByteCount(pos,kTRUE);

               continue;
            }

         case TStreamerInfo::kObject:   // Class      derived from TObject
         case TStreamerInfo::kAny:   // Class  NOT derived from TObject
         case TStreamerInfo::kOffsetL + TStreamerInfo::kObject:
         case TStreamerInfo::kAny     + TStreamerInfo::kOffsetL: {
            TClass *cl                 = compinfo[i]->fClass;
            TMemberStreamer *pstreamer = compinfo[i]->fStreamer;
            DOLOOP
               {b.WriteFastArray((void*)(arr[k]+ioffset),cl,compinfo[i]->fLength,pstreamer);}
            continue;
         }

         case TStreamerInfo::kOffsetL + TStreamerInfo::kTString:
         case TStreamerInfo::kOffsetL + TStreamerInfo::kTObject:
         case TStreamerInfo::kOffsetL + TStreamerInfo::kTNamed:
         {
            TMemberStreamer *pstreamer = compinfo[i]->fStreamer;
            TClass *cl                 = compinfo[i]->fClass;

            UInt_t pos = b.WriteVersion(this->IsA(),kTRUE);
            DOLOOP {b.WriteFastArray((void*)(arr[k]+ioffset),cl,compinfo[i]->fLength,pstreamer);}
            b.SetByteCount(pos,kTRUE);
            continue;
         }

            // Base Class
         case TStreamerInfo::kBase:
            if (!(arrayMode&1)) {
               TMemberStreamer *pstreamer = compinfo[i]->fStreamer;
               if(pstreamer) {
                  // See kStreamer case (similar code)
                  UInt_t pos = b.WriteVersion(this->IsA(),kTRUE);
                  DOLOOP{(*pstreamer)(b,arr[k]+ioffset,compinfo[i]->fLength);}
                  b.SetByteCount(pos,kTRUE);
               } else {
                  DOLOOP { ((TStreamerBase*)aElement)->WriteBuffer(b,arr[k]);}
               }
            } else {
               TClass *cl                 = compinfo[i]->fClass;
               TStreamerInfo *binfo = ((TStreamerInfo*)cl->GetStreamerInfo());
               binfo->WriteBufferAux(b,arr,binfo->fCompFull,0,binfo->fNfulldata,narr,ioffset,arrayMode);
            }
            continue;

         case TStreamerInfo::kStreamer:
         {
            TMemberStreamer *pstreamer = compinfo[i]->fStreamer;

            UInt_t pos = b.WriteVersion(this->IsA(),kTRUE);
            if (pstreamer == 0) {
               printf("ERROR, Streamer is null\n");
               aElement->ls();continue;
            } else {
               DOLOOP{(*pstreamer)(b,arr[k]+ioffset,compinfo[i]->fLength);}
            }
            b.SetByteCount(pos,kTRUE);
         }
         continue;

         case TStreamerInfo::kStreamLoop:
            // -- A pointer to a varying-length array of objects.
            // MyClass* ary; //[n]
            // -- Or a pointer to a varying-length array of pointers to objects.
            // MyClass** ary; //[n]
         case TStreamerInfo::kOffsetL + TStreamerInfo::kStreamLoop:
            // -- An array of pointers to a varying-length array of objects.
            // MyClass* ary[d]; //[n]
            // -- Or an array of pointers to a varying-length array of pointers to objects.
            // MyClass** ary[d]; //[n]
         {
            // Get the class of the data member.
            TClass* cl = compinfo[i]->fClass;
            // Get any private streamer which was set for the data member.
            TMemberStreamer* pstreamer = compinfo[i]->fStreamer;
            // Which are we, an array of objects or an array of pointers to objects?
            Bool_t isPtrPtr = (strstr(aElement->GetTypeName(), "**") != 0);
            if (pstreamer) {
               // -- We have a private streamer.
               UInt_t pos = b.WriteVersion(this->IsA(), kTRUE);
               // Loop over the entries in the clones array or the STL container.
               for (int k = 0; k < narr; ++k) {
                  // Get a pointer to the counter for the varying length array.
                  Int_t* counter = (Int_t*) (arr[k] /*entry pointer*/ + eoffset /*entry offset*/ + compinfo[i]->fMethod /*counter offset*/);
                  // And call the private streamer, passing it the buffer, the object, and the counter.
                  (*pstreamer)(b, arr[k] /*entry pointer*/ + ioffset /*object offset*/, *counter);
               }
               b.SetByteCount(pos, kTRUE);
               // We are done, next streamer element.
               continue;
            }
            // At this point we do *not* have a private streamer.
            // Get the version of the file we are writing to.
            TFile* file = (TFile*) b.GetParent();
            // By default assume the file version is the newest.
            Int_t fileVersion = kMaxInt;
            if (file) {
               fileVersion = file->GetVersion();
            }
            // Write the class version to the buffer.
            UInt_t pos = b.WriteVersion(this->IsA(), kTRUE);
            if (fileVersion > 51508) {
               // -- Newer versions allow polymorphic pointers to objects.
               // Loop over the entries in the clones array or the STL container.
               for (int k = 0; k < narr; ++k) {
                  // Get the counter for the varying length array.
                  Int_t vlen = *((Int_t*) (arr[k] /*entry pointer*/ + eoffset /*entry offset*/ + compinfo[i]->fMethod /*counter offset*/));
                  //b << vlen;
                  if (vlen) {
                     // Get a pointer to the array of pointers.
                     char** pp = (char**) (arr[k] /*entry pointer*/ + ioffset /*object offset*/);
                     // Loop over each element of the array of pointers to varying-length arrays.
                     for (Int_t ndx = 0; ndx < compinfo[i]->fLength; ++ndx) {
                        if (!pp[ndx]) {
                           // -- We do not have a pointer to a varying-length array.
                           Error("WriteBufferAux", "The pointer to element %s::%s type %d (%s) is null\n", GetName(), aElement->GetFullName(), compinfo[i]->fType, aElement->GetTypeName());
                           continue;
                        }
                        if (!isPtrPtr) {
                           // -- We are a varying-length array of objects.
                           // Write the entire array of objects to the buffer.
                           // Note: Polymorphism is not allowed here.
                           b.WriteFastArray(pp[ndx], cl, vlen, 0);
                        }
                        else {
                           // -- We are a varying-length array of pointers to objects.
                           // Write the entire array of object pointers to the buffer.
                           // Note: The object pointers are allowed to be polymorphic.
                           b.WriteFastArray((void**) pp[ndx], cl, vlen, kFALSE, 0);
                        } // isPtrPtr
                     } // ndx
                  } // vlen
               } // k
            }
            else {
               // -- Older versions do *not* allow polymorphic pointers to objects.
               // Loop over the entries in the clones array or the STL container.
               for (int k = 0; k < narr; ++k) {
                  // Get the counter for the varying length array.
                  Int_t vlen = *((Int_t*) (arr[k] /*entry pointer*/ + eoffset /*entry offset*/ + compinfo[i]->fMethod /*counter offset*/));
                  //b << vlen;
                  if (vlen) {
                     // Get a pointer to the array of pointers.
                     char** pp = (char**) (arr[k] /*entry pointer*/ + ioffset /*object offset*/);
                     // -- Older versions do *not* allow polymorphic pointers to objects.
                     // Loop over each element of the array of pointers to varying-length arrays.
                     for (Int_t ndx = 0; ndx < compinfo[i]->fLength; ++ndx) {
                        if (!pp[ndx]) {
                           // -- We do not have a pointer to a varying-length array.
                           Error("WriteBufferAux", "The pointer to element %s::%s type %d (%s) is null\n", GetName(), aElement->GetFullName(), compinfo[i]->fType, aElement->GetTypeName());
                           continue;
                        }
                        if (!isPtrPtr) {
                           // -- We are a varying-length array of objects.
                           // Loop over the elements of the varying length array.
                           for (Int_t v = 0; v < vlen; ++v) {
                              // Write the object to the buffer.
                              cl->Streamer(pp[ndx] + (v * cl->Size()), b);
                           } // v
                        }
                        else {
                           // -- We are a varying-length array of pointers to objects.
                           // Loop over the elements of the varying length array.
                           for (Int_t v = 0; v < vlen; ++v) {
                              // Get a pointer to the object pointer.
                              char** r = (char**) pp[ndx];
                              // Write the object to the buffer.
                              cl->Streamer(r[v], b);
                           } // v
                        } // isPtrPtr
                     } // ndx
                  } // vlen
               } // k
            } // fileVersion
            // Backpatch the byte count into the buffer.
            b.SetByteCount(pos, kTRUE);
            continue;
         }

         case TStreamerInfo::kCacheNew:
            ((TBufferFile&)b).PushDataCache( new TVirtualArray( aElement->GetClassPointer(), narr ) );
            continue;
         case TStreamerInfo::kCacheDelete:
            delete ((TBufferFile&)b).PopDataCache();
            continue;
         case TStreamerInfo::kArtificial:
#if 0
            ROOT::TSchemaRule::WriteFuncPtr_t writefunc = ((TStreamerArtificial*)aElement)->GetWriteFunc();
            if (writefunc) {
               DOLOOP( writefunc(arr[k]+eoffset, b) );
            }
#endif
            continue;
         case -1:
            // -- Skip an ignored TObject base class.
            continue;
         default:
            Error("WriteBuffer","The element %s::%s type %d (%s) is not supported yet\n",GetName(),aElement->GetFullName(),compinfo[i]->fType,aElement->GetTypeName());
            continue;
      }
   }

   if (needIncrement) b.DecrementLevel(this);

   return 0;
}

template Int_t TStreamerInfo::WriteBufferAux<char**>(TBuffer &b, char ** const &arr, TCompInfo *const*const compinfo, Int_t first, Int_t last, Int_t narr,Int_t eoffset,Int_t mode);

////////////////////////////////////////////////////////////////////////////////
/// Write for STL container. ('first' is an id between -1 and fNfulldata).

Int_t TStreamerInfo::WriteBufferSTL(TBuffer &b, TVirtualCollectionProxy *cont, Int_t nc)
{
   if (!nc) return 0;
   R__ASSERT((unsigned int)nc==cont->Size());


   int ret = WriteBufferAux(b,*cont,fCompFull,0,fNfulldata,nc,/* eoffset = */ 0,1);
   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// Write for STL container.  ('first' is an id between -1 and fNfulldata).
/// Note: This is no longer used.

Int_t TStreamerInfo::WriteBufferSTLPtrs(TBuffer &b, TVirtualCollectionProxy *cont, Int_t nc, Int_t first, Int_t eoffset )
{
   if (!nc) return 0;
   R__ASSERT((unsigned int)nc==cont->Size());
   int ret = WriteBufferAux(b, TPointerCollectionAdapter(cont),fCompFull,first==-1?0:first,first==-1?fNfulldata:first+1,nc,eoffset,1);
   return ret;
}


////////////////////////////////////////////////////////////////////////////////
/// General Write.  ('first' is an id between -1 and fNdata).
/// Note: This is no longer used.

Int_t TStreamerInfo::WriteBuffer(TBuffer &b, char *ipointer, Int_t first)
{
   return WriteBufferAux(b,&ipointer,fCompOpt,first==-1?0:first,first==-1?fNdata:first+1,1,0,0);
}

////////////////////////////////////////////////////////////////////////////////
/// Write for ClonesArray ('first' is an id between -1 and fNfulldata).
/// Note: This is no longer used.

Int_t TStreamerInfo::WriteBufferClones(TBuffer &b, TClonesArray *clones,
                                       Int_t nc, Int_t first, Int_t eoffset)
{
   char **arr = reinterpret_cast<char**>(clones->GetObjectRef(0));
   return WriteBufferAux(b,arr,fCompFull,first==-1?0:first,first==-1?fNfulldata:first+1,nc,eoffset,1);
}


