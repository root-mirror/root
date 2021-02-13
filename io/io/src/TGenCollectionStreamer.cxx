// @(#)root/io:$Id$
// Author: Markus Frank 28/10/04

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/*
\class TGenCollectionStreamer
\ingroup IO

Streamer around an arbitrary container, which implements basic
functionality and iteration.

In particular this is used to implement splitting and abstract
element access of any container. Access to compiled code is necessary
to implement the abstract iteration sequence and functionality like
size(), clear(), resize(). resize() may be a void operation.
**/

#include "TGenCollectionStreamer.h"
#include "TClassEdit.h"
#include "TError.h"
#include "TROOT.h"
#include "TStreamerInfo.h"
#include "TStreamerElement.h"
#include "TVirtualCollectionIterators.h"

TGenCollectionStreamer::TGenCollectionStreamer(const TGenCollectionStreamer& copy)
      : TGenCollectionProxy(copy), fReadBufferFunc(&TGenCollectionStreamer::ReadBufferDefault)
{
   // Build a Streamer for an emulated vector whose type is 'name'.
}

TGenCollectionStreamer::TGenCollectionStreamer(Info_t info, size_t iter_size)
      : TGenCollectionProxy(info, iter_size), fReadBufferFunc(&TGenCollectionStreamer::ReadBufferDefault)
{
   // Build a Streamer for a collection whose type is described by 'collectionClass'.
}

TGenCollectionStreamer::TGenCollectionStreamer(const ::ROOT::TCollectionProxyInfo &info, TClass *cl)
      : TGenCollectionProxy(info, cl), fReadBufferFunc(&TGenCollectionStreamer::ReadBufferDefault)
{
   // Build a Streamer for a collection whose type is described by 'collectionClass'.
}

TGenCollectionStreamer::~TGenCollectionStreamer()
{
   // Standard destructor.
}

TVirtualCollectionProxy* TGenCollectionStreamer::Generate() const
{
   // Virtual copy constructor.
   if (!fValue.load()) Initialize(kFALSE);
   return new TGenCollectionStreamer(*this);
}

template <typename T>
T* getaddress(TGenCollectionProxy::StreamHelper &itm);

template <>
bool* getaddress<bool>(TGenCollectionProxy::StreamHelper &itm)
{
   return &itm.boolean;
}

template <>
Char_t* getaddress<Char_t>(TGenCollectionProxy::StreamHelper &itm)
{
   return &itm.s_char;
}

template <>
Short_t* getaddress<Short_t>(TGenCollectionProxy::StreamHelper &itm)
{
   return &itm.s_short;
}

template <>
Int_t* getaddress<Int_t>(TGenCollectionProxy::StreamHelper &itm)
{
   return &itm.s_int;
}

template <>
Long_t* getaddress<Long_t>(TGenCollectionProxy::StreamHelper &itm)
{
   return &itm.s_long;
}

template <>
Long64_t* getaddress<Long64_t>(TGenCollectionProxy::StreamHelper &itm)
{
   return &itm.s_longlong;
}

template <>
Float_t* getaddress<Float_t>(TGenCollectionProxy::StreamHelper &itm)
{
   return &itm.flt;
}

template <>
Double_t* getaddress<Double_t>(TGenCollectionProxy::StreamHelper &itm)
{
   return &itm.dbl;
}

template <>
UChar_t* getaddress<UChar_t>(TGenCollectionProxy::StreamHelper &itm)
{
   return &itm.u_char;
}

template <>
UShort_t* getaddress<UShort_t>(TGenCollectionProxy::StreamHelper &itm)
{
   return &itm.u_short;
}

template <>
UInt_t* getaddress<UInt_t>(TGenCollectionProxy::StreamHelper &itm)
{
   return &itm.u_int;
}

template <>
ULong_t* getaddress<ULong_t>(TGenCollectionProxy::StreamHelper &itm)
{
   return &itm.u_long;
}

template <>
ULong64_t* getaddress<ULong64_t>(TGenCollectionProxy::StreamHelper &itm)
{
   return &itm.u_longlong;
}

template <typename From, typename To>
void ConvertArray(TGenCollectionProxy::StreamHelper *read, TGenCollectionProxy::StreamHelper *write, int nElements)
{
   From *r = getaddress<From>( *read );
   To *w = getaddress<To>( *write );
   for(int i = 0; i < nElements; ++i) {
      // getvalue<To>( write[i] ) = (To) getvalue<From>( read[i] );
      w[i] = (To)r[i];
   }
}

template <typename From>
void DispatchConvertArray(int writeType, TGenCollectionProxy::StreamHelper *read, TGenCollectionProxy::StreamHelper *write, int nElements)
{
   switch(writeType) {
      case kBool_t:
         ConvertArray<From,bool>(read,write,nElements);
         break;
      case kChar_t:
         ConvertArray<From,Char_t>(read,write,nElements);
         break;
      case kShort_t:
         ConvertArray<From,Short_t>(read,write,nElements);
         break;
      case kInt_t:
         ConvertArray<From,Int_t>(read,write,nElements);
         break;
      case kLong_t:
         ConvertArray<From,Long64_t>(read,write,nElements);
         break;
      case kLong64_t:
         ConvertArray<From,Long64_t>(read,write,nElements);
         break;
      case kFloat_t:
         ConvertArray<From,Float_t>(read,write,nElements);
         break;
      case kFloat16_t:
         ConvertArray<From,Float16_t>(read,write,nElements);
         break;
      case kDouble_t:
         ConvertArray<From,Double_t>(read,write,nElements);
         break;
      case kUChar_t:
         ConvertArray<From,UChar_t>(read,write,nElements);
         break;
      case kUShort_t:
         ConvertArray<From,UShort_t>(read,write,nElements);
         break;
      case kUInt_t:
         ConvertArray<From,UInt_t>(read,write,nElements);
         break;
      case kULong_t:
         ConvertArray<From,ULong_t>(read,write,nElements);
         break;
      case kULong64_t:
         ConvertArray<From,ULong64_t>(read,write,nElements);
         break;
      case kDouble32_t:
         ConvertArray<From,Double32_t>(read,write,nElements);
         break;
      case kchar:
      case kNoType_t:
      case kOther_t:
         Error("TGenCollectionStreamer", "fType %d is not supported yet!\n", writeType);
   }
}

void TGenCollectionStreamer::ReadPrimitives(int nElements, TBuffer &b, const TClass *onFileClass)
{
   // Primitive input streamer.
   size_t len = fValDiff * nElements;
   char   buffer[8096];
   Bool_t   feed = false;
   void*  memory = 0;
   StreamHelper* itmstore = 0;
   StreamHelper* itmconv = 0;
   fEnv->fSize = nElements;
   switch (fSTL_type)  {
      case ROOT::kSTLvector:
         if (fVal->fKind != kBool_t)  {
            fResize(fEnv->fObject,fEnv->fSize);
            fEnv->fIdx = 0;

            TVirtualVectorIterators iterators(fFunctionCreateIterators);
            iterators.CreateIterators(fEnv->fObject);
            itmstore = (StreamHelper*)iterators.fBegin;
            fEnv->fStart = itmstore;
            break;
         }
      default:
         feed = true;
         itmstore = (StreamHelper*)(len < sizeof(buffer) ? buffer : memory =::operator new(len));
         break;
   }
   fEnv->fStart = itmstore;

   StreamHelper *itmread;
   int readkind;
   if (onFileClass) {
      readkind = onFileClass->GetCollectionProxy()->GetType();
      itmconv = (StreamHelper*) ::operator new( nElements * onFileClass->GetCollectionProxy()->GetIncrement() );
      itmread = itmconv;
   } else {
      itmread = itmstore;
      readkind = fVal->fKind;
   }
   switch (readkind)   {
      case kBool_t:
         b.ReadFastArray(&itmread->boolean   , nElements);
         break;
      case kChar_t:
         b.ReadFastArray(&itmread->s_char    , nElements);
         break;
      case kShort_t:
         b.ReadFastArray(&itmread->s_short   , nElements);
         break;
      case kInt_t:
         b.ReadFastArray(&itmread->s_int     , nElements);
         break;
      case kLong_t:
         b.ReadFastArray(&itmread->s_long    , nElements);
         break;
      case kLong64_t:
         b.ReadFastArray(&itmread->s_longlong, nElements);
         break;
      case kFloat_t:
         b.ReadFastArray(&itmread->flt       , nElements);
         break;
      case kFloat16_t:
         b.ReadFastArrayFloat16(&itmread->flt, nElements);
         break;
      case kDouble_t:
         b.ReadFastArray(&itmread->dbl       , nElements);
         break;
      case kUChar_t:
         b.ReadFastArray(&itmread->u_char    , nElements);
         break;
      case kUShort_t:
         b.ReadFastArray(&itmread->u_short   , nElements);
         break;
      case kUInt_t:
         b.ReadFastArray(&itmread->u_int     , nElements);
         break;
      case kULong_t:
         b.ReadFastArray(&itmread->u_long    , nElements);
         break;
      case kULong64_t:
         b.ReadFastArray(&itmread->u_longlong, nElements);
         break;
      case kDouble32_t:
         b.ReadFastArrayDouble32(&itmread->dbl, nElements);
         break;
      case kchar:
      case kNoType_t:
      case kOther_t:
         Error("TGenCollectionStreamer", "fType %d is not supported yet!\n", readkind);
   }
   if (onFileClass) {
      switch (readkind)   {
         case kBool_t:
            DispatchConvertArray<bool>(fVal->fKind, itmread, itmstore, nElements);
            break;
         case kChar_t:
            DispatchConvertArray<Char_t>(fVal->fKind, itmread, itmstore, nElements);
            break;
         case kShort_t:
            DispatchConvertArray<Short_t>(fVal->fKind, itmread, itmstore, nElements);
            break;
         case kInt_t:
            DispatchConvertArray<Int_t>(fVal->fKind, itmread, itmstore, nElements);
            break;
         case kLong_t:
            DispatchConvertArray<Long_t>(fVal->fKind, itmread, itmstore, nElements);
            break;
         case kLong64_t:
            DispatchConvertArray<Long64_t>(fVal->fKind, itmread, itmstore, nElements);
            break;
         case kFloat_t:
            DispatchConvertArray<Float_t>(fVal->fKind, itmread, itmstore, nElements);
            break;
         case kFloat16_t:
            DispatchConvertArray<Float16_t>(fVal->fKind, itmread, itmstore, nElements);
            break;
         case kDouble_t:
            DispatchConvertArray<Double_t>(fVal->fKind, itmread, itmstore, nElements);
            break;
         case kUChar_t:
            DispatchConvertArray<UChar_t>(fVal->fKind, itmread, itmstore, nElements);
            break;
         case kUShort_t:
            DispatchConvertArray<UShort_t>(fVal->fKind, itmread, itmstore, nElements);
            break;
         case kUInt_t:
            DispatchConvertArray<UInt_t>(fVal->fKind, itmread, itmstore, nElements);
            break;
         case kULong_t:
            DispatchConvertArray<ULong_t>(fVal->fKind, itmread, itmstore, nElements);
            break;
         case kULong64_t:
            DispatchConvertArray<ULong64_t>(fVal->fKind, itmread, itmstore, nElements);
            break;
         case kDouble32_t:
            DispatchConvertArray<Double_t>(fVal->fKind, itmread, itmstore, nElements);
            break;
         case kchar:
         case kNoType_t:
         case kOther_t:
            Error("TGenCollectionStreamer", "fType %d is not supported yet!\n", readkind);
      }
      ::operator delete((void*)itmconv);
   }
   if (feed)  {      // need to feed in data...
      fEnv->fStart = fFeed(itmstore,fEnv->fObject,fEnv->fSize);
      if (memory)  {
         ::operator delete(memory);
      }
   }
}

void TGenCollectionStreamer::ReadObjects(int nElements, TBuffer &b, const TClass *onFileClass)
{
   // Object input streamer.
   Bool_t vsn3 = b.GetInfo() && b.GetInfo()->GetOldVersion() <= 3;
   size_t len = fValDiff * nElements;
   StreamHelper* itm = 0;
   char   buffer[8096];
   void*  memory = 0;

   TClass* onFileValClass = (onFileClass ? onFileClass->GetCollectionProxy()->GetValueClass() : 0);

   fEnv->fSize = nElements;
   switch (fSTL_type)  {
         // Simple case: contiguous memory. get address of first, then jump.
      case ROOT::kSTLvector:
#define DOLOOP(x) {int idx=0; while(idx<nElements) {StreamHelper* i=(StreamHelper*)(((char*)itm) + fValDiff*idx); { x ;} ++idx;} break;}
         fResize(fEnv->fObject,fEnv->fSize);
         fEnv->fIdx = 0;

         {
            TVirtualVectorIterators iterators(fFunctionCreateIterators);
            iterators.CreateIterators(fEnv->fObject);
            itm = (StreamHelper*)iterators.fBegin;
         }
         fEnv->fStart = itm;
         switch (fVal->fCase) {
            case kIsClass:
               DOLOOP(b.StreamObject(i, fVal->fType, onFileValClass ));
            case EProperty(kBIT_ISSTRING):
               DOLOOP(i->read_std_string(b));
            case EProperty(kIsPointer | kIsClass):
               DOLOOP(i->set(b.ReadObjectAny(fVal->fType)));
            case EProperty(kIsPointer | kBIT_ISSTRING):
               DOLOOP(i->read_std_string_pointer(b));
            case EProperty(kIsPointer | kBIT_ISTSTRING | kIsClass):
               DOLOOP(i->read_tstring_pointer(vsn3, b));
         }
#undef DOLOOP
         break;

         // No contiguous memory, but resize is possible
         // Hence accessing objects using At(i) should be not too much an overhead
      case ROOT::kSTLlist:
      case ROOT::kSTLforwardlist:
      case ROOT::kSTLdeque:
#define DOLOOP(x) {int idx=0; while(idx<nElements) {StreamHelper* i=(StreamHelper*)TGenCollectionProxy::At(idx); { x ;} ++idx;} break;}
         fResize(fEnv->fObject,fEnv->fSize);
         fEnv->fIdx = 0;
         fEnv->fStart = 0;
         switch (fVal->fCase) {
            case kIsClass:
               DOLOOP(b.StreamObject(i, fVal->fType, onFileValClass));
            case kBIT_ISSTRING:
               DOLOOP(i->read_std_string(b));
            case kIsPointer | kIsClass:
               DOLOOP(i->set(b.ReadObjectAny(fVal->fType)));
            case kIsPointer | kBIT_ISSTRING:
               DOLOOP(i->read_std_string_pointer(b));
            case kIsPointer | kBIT_ISTSTRING | kIsClass:
               DOLOOP(i->read_tstring_pointer(vsn3, b));
         }
#undef DOLOOP
         break;

         // Rather troublesome case: Objects can only be fed into the container
         // Once they are created. Need to take memory from stack or heap.
      case ROOT::kSTLmultiset:
      case ROOT::kSTLset:
      case ROOT::kSTLunorderedset:
      case ROOT::kSTLunorderedmultiset:
#define DOLOOP(x) {int idx=0; while(idx<nElements) {StreamHelper* i=(StreamHelper*)(((char*)itm) + fValDiff*idx); { x ;} ++idx;}}
         fEnv->fStart = itm = (StreamHelper*)(len < sizeof(buffer) ? buffer : memory =::operator new(len));
         fConstruct(itm,nElements);
         switch (fVal->fCase) {
            case kIsClass:
               DOLOOP(b.StreamObject(i, fVal->fType, onFileValClass));
               fFeed(fEnv->fStart,fEnv->fObject,fEnv->fSize);
               fDestruct(fEnv->fStart,fEnv->fSize);
               break;
            case EProperty(kBIT_ISSTRING):
               DOLOOP(i->read_std_string(b))
               fFeed(fEnv->fStart,fEnv->fObject,fEnv->fSize);
               fDestruct(fEnv->fStart,fEnv->fSize);
               break;
            case EProperty(kIsPointer | kIsClass):
               DOLOOP(i->set(b.ReadObjectAny(fVal->fType)));
               fFeed(fEnv->fStart,fEnv->fObject,fEnv->fSize);
               break;
            case EProperty(kIsPointer | kBIT_ISSTRING):
               DOLOOP(i->read_std_string_pointer(b))
               fFeed(fEnv->fStart,fEnv->fObject,fEnv->fSize);
               break;
            case EProperty(kIsPointer | kBIT_ISTSTRING | kIsClass):
               DOLOOP(i->read_tstring_pointer(vsn3, b));
               fFeed(fEnv->fStart,fEnv->fObject,fEnv->fSize);
               break;
         }
#undef DOLOOP
         break;
      default:
         break;
   }
   if (memory) {
      ::operator delete(memory);
   }
}

void TGenCollectionStreamer::ReadPairFromMap(int nElements, TBuffer &b)
{
   // Input streamer to convert a map into another collection

   Bool_t vsn3 = b.GetInfo() && b.GetInfo()->GetOldVersion() <= 3;
   size_t len = fValDiff * nElements;
   StreamHelper* itm = 0;
   char   buffer[8096];
   void*  memory = 0;

   TStreamerInfo *pinfo = (TStreamerInfo*)fVal->fType->GetStreamerInfo();
   R__ASSERT(pinfo);
   R__ASSERT(fVal->fCase == kIsClass);

   int nested = 0;
   std::vector<std::string> inside;
   TClassEdit::GetSplit(pinfo->GetName(), inside, nested);
   Value first(inside[1],kFALSE);
   Value second(inside[2],kFALSE);
   fValOffset = ((TStreamerElement*)pinfo->GetElements()->At(1))->GetOffset();

   fEnv->fSize = nElements;
   switch (fSTL_type)  {
         // Simple case: contiguous memory. get address of first, then jump.
      case ROOT::kSTLvector:
#define DOLOOP(x) {int idx=0; while(idx<nElements) {StreamHelper* i=(StreamHelper*)(((char*)itm) + fValDiff*idx); { x ;} ++idx;} break;}
         fResize(fEnv->fObject,fEnv->fSize);
         fEnv->fIdx = 0;

         {
            TVirtualVectorIterators iterators(fFunctionCreateIterators);
            iterators.CreateIterators(fEnv->fObject);
            itm = (StreamHelper*)iterators.fBegin;
         }
         fEnv->fStart = itm;
         switch (fVal->fCase) {
            case kIsClass:
               DOLOOP(
                  ReadMapHelper(i, &first, vsn3, b);
                  ReadMapHelper((StreamHelper*)(((char*)i) + fValOffset), &second, vsn3, b)
               );
         }
#undef DOLOOP
         break;

         // No contiguous memory, but resize is possible
         // Hence accessing objects using At(i) should be not too much an overhead
      case ROOT::kSTLlist:
      case ROOT::kSTLforwardlist:
      case ROOT::kSTLdeque:
#define DOLOOP(x) {int idx=0; while(idx<nElements) {StreamHelper* i=(StreamHelper*)TGenCollectionProxy::At(idx); { x ;} ++idx;} break;}
         fResize(fEnv->fObject,fEnv->fSize);
         fEnv->fIdx = 0;
         {
            TVirtualVectorIterators iterators(fFunctionCreateIterators);
            iterators.CreateIterators(fEnv->fObject);
            fEnv->fStart = iterators.fBegin;
         }
         switch (fVal->fCase) {
            case kIsClass:
               DOLOOP(
                  char **where = (char**)(void*) & i;
                  b.ApplySequence(*(pinfo->GetReadObjectWiseActions()), where);
               );
         }
#undef DOLOOP
         break;

         // Rather troublesome case: Objects can only be fed into the container
         // Once they are created. Need to take memory from stack or heap.
      case ROOT::kSTLmultiset:
      case ROOT::kSTLset:
      case ROOT::kSTLunorderedset:
      case ROOT::kSTLunorderedmultiset:
#define DOLOOP(x) {int idx=0; while(idx<nElements) {StreamHelper* i=(StreamHelper*)(((char*)itm) + fValDiff*idx); { x ;} ++idx;}}
         fEnv->fStart = itm = (StreamHelper*)(len < sizeof(buffer) ? buffer : memory =::operator new(len));
         fConstruct(itm,nElements);
         switch (fVal->fCase) {
            case kIsClass:
               DOLOOP(
                  char **where = (char**)(void*) & i;
                  b.ApplySequence(*(pinfo->GetReadObjectWiseActions()), where);
               );
               fFeed(fEnv->fStart,fEnv->fObject,fEnv->fSize);
               fDestruct(fEnv->fStart,fEnv->fSize);
               break;
         }
#undef DOLOOP
         break;
      default:
         break;
   }
   if (memory) {
      ::operator delete(memory);
   }
}


void TGenCollectionStreamer::ReadMapHelper(StreamHelper *i, Value *v, Bool_t vsn3,  TBuffer &b)
{
   // helper class to read std::map

   float f;

   switch (v->fCase) {
      case kIsFundamental:  // Only handle primitives this way
      case kIsEnum:
         switch (int(v->fKind))   {
            case kBool_t:
               b >> i->boolean;
               break;
            case kChar_t:
               b >> i->s_char;
               break;
            case kShort_t:
               b >> i->s_short;
               break;
            case kInt_t:
               b >> i->s_int;
               break;
            case kLong_t:
               b >> i->s_long;
               break;
            case kLong64_t:
               b >> i->s_longlong;
               break;
            case kFloat_t:
               b >> i->flt;
               break;
            case kFloat16_t:
               b >> f;
               i->flt = float(f);
               break;
            case kDouble_t:
               b >> i->dbl;
               break;
            case kUChar_t:
               b >> i->u_char;
               break;
            case kUShort_t:
               b >> i->u_short;
               break;
            case kUInt_t:
               b >> i->u_int;
               break;
            case kULong_t:
               b >> i->u_long;
               break;
            case kULong64_t:
               b >> i->u_longlong;
               break;
            case kDouble32_t:
               b >> f;
               i->dbl = double(f);
               break;
            case kchar:
            case kNoType_t:
            case kOther_t:
               Error("TGenCollectionStreamer", "fType %d is not supported yet!\n", v->fKind);
         }
         break;
      case kIsClass:
         b.StreamObject(i, v->fType);
         break;
      case kBIT_ISSTRING:
         i->read_std_string(b);
         break;
      case kIsPointer | kIsClass:
         i->set(b.ReadObjectAny(v->fType));
         break;
      case kIsPointer | kBIT_ISSTRING:
         i->read_std_string_pointer(b);
         break;
      case kIsPointer | kBIT_ISTSTRING | kIsClass:
         i->read_tstring_pointer(vsn3, b);
         break;
   }
}

template <typename To>
To readOneValue(TBuffer &b, int readtype) {
   TGenCollectionProxy::StreamHelper itm;
   TGenCollectionProxy::StreamHelper *i = &itm;
   switch (readtype)   {
   case kBool_t:
      b >> i->boolean;
      return (To)i->boolean;
      break;
   case kChar_t:
      b >> i->s_char;
      return (To)i->s_char;
      break;
   case kShort_t:
      b >> i->s_short;
      return (To)i->s_short;
      break;
   case kInt_t:
      b >> i->s_int;
      return (To)i->s_int;
      break;
   case kLong_t:
      b >> i->s_long;
      return (To)i->s_long;
      break;
   case kLong64_t:
      b >> i->s_longlong;
      return (To)i->s_longlong;
      break;
   case kFloat_t:
      b >> i->flt;
      return (To)i->flt;
      break;
   case kFloat16_t:
      b >> i->flt;
      return (To)i->flt;
      break;
   case kDouble_t:
      b >> i->dbl;
      return (To)i->dbl;
      break;
   case kUChar_t:
      b >> i->u_char;
      return (To)i->u_char;
      break;
   case kUShort_t:
      b >> i->u_short;
      return (To)i->u_short;
      break;
   case kUInt_t:
      b >> i->u_int;
      return (To)i->u_int;
      break;
   case kULong_t:
      b >> i->u_long;
      return (To)i->u_long;
      break;
   case kULong64_t:
      b >> i->u_longlong;
      return (To)i->u_longlong;
      break;
   case kDouble32_t: {
      float f;
      b >> f;
      i->dbl = double(f);
      return (To)i->dbl;
      break;
   }
   case kchar:
   case kNoType_t:
   case kOther_t:
      Error("TGenCollectionStreamer", "fType %d is not supported yet!\n", readtype);
   }
   return 0;
}


void TGenCollectionStreamer::ReadMap(int nElements, TBuffer &b, const TClass *onFileClass)
{
   // Map input streamer.
   Bool_t vsn3 = b.GetInfo() && b.GetInfo()->GetOldVersion() <= 3;
   size_t len = fValDiff * nElements;
   Value  *v;
   char buffer[8096], *addr, *temp;
   void* memory = 0;
   StreamHelper* i;
   float f;
   fEnv->fSize  = nElements;
   fEnv->fStart = (len < sizeof(buffer) ? buffer : memory =::operator new(len));
   addr = temp = (char*)fEnv->fStart;
   fConstruct(addr,nElements);

   int onFileValueKind[2];
   if (onFileClass) {
      TClass *onFileValueClass = onFileClass->GetCollectionProxy()->GetValueClass();
      TVirtualStreamerInfo *sourceInfo = onFileValueClass->GetStreamerInfo();
      onFileValueKind[0] = ((TStreamerElement*)sourceInfo->GetElements()->At(0))->GetType();
      onFileValueKind[1] = ((TStreamerElement*)sourceInfo->GetElements()->At(1))->GetType();
   }
   for (int loop, idx = 0; idx < nElements; ++idx)  {
      addr = temp + fValDiff * idx;
      v = fKey;
      for (loop = 0; loop < 2; loop++)  {
         i = (StreamHelper*)addr;
         switch (v->fCase) {
            case kIsFundamental:  // Only handle primitives this way
            case kIsEnum:
               if (onFileClass) {
                  int readtype = (int)(onFileValueKind[loop]);
                  switch (int(v->fKind))   {
                  case kBool_t:
                     i->boolean = readOneValue<bool>(b,readtype);
                     break;
                  case kChar_t:
                     i->s_char = readOneValue<Char_t>(b,readtype);
                     break;
                  case kShort_t:
                     i->s_short = readOneValue<Short_t>(b,readtype);
                     break;
                  case kInt_t:
                     i->s_int = readOneValue<Int_t>(b,readtype);
                     break;
                  case kLong_t:
                     i->s_long = readOneValue<Long_t>(b,readtype);
                     break;
                  case kLong64_t:
                     i->s_longlong = readOneValue<Long64_t>(b,readtype);
                     break;
                  case kFloat_t:
                     i->flt = readOneValue<Float_t>(b,readtype);
                     break;
                  case kFloat16_t:
                     i->flt = readOneValue<Float16_t>(b,readtype);
                     break;
                  case kDouble_t:
                     i->dbl = readOneValue<Double_t>(b,readtype);
                     break;
                  case kUChar_t:
                     i->u_char = readOneValue<UChar_t>(b,readtype);
                     break;
                  case kUShort_t:
                     i->u_short = readOneValue<UShort_t>(b,readtype);
                     break;
                  case kUInt_t:
                     i->u_int = readOneValue<UInt_t>(b,readtype);
                     break;
                  case kULong_t:
                     i->u_long = readOneValue<ULong_t>(b,readtype);
                     break;
                  case kULong64_t:
                     i->u_longlong = readOneValue<ULong64_t>(b,readtype);
                     break;
                  case kDouble32_t:
                     i->dbl = readOneValue<Double32_t>(b,readtype);
                     break;
                  case kchar:
                  case kNoType_t:
                  case kOther_t:
                     Error("TGenCollectionStreamer", "fType %d is not supported yet!\n", v->fKind);
                  }
               } else {
                  switch (int(v->fKind))   {
                  case kBool_t:
                     b >> i->boolean;
                     break;
                  case kChar_t:
                     b >> i->s_char;
                     break;
                  case kShort_t:
                     b >> i->s_short;
                     break;
                  case kInt_t:
                     b >> i->s_int;
                     break;
                  case kLong_t:
                     b >> i->s_long;
                     break;
                  case kLong64_t:
                     b >> i->s_longlong;
                     break;
                  case kFloat_t:
                     b >> i->flt;
                     break;
                  case kFloat16_t:
                     b >> f;
                     i->flt = float(f);
                     break;
                  case kDouble_t:
                     b >> i->dbl;
                     break;
                  case kUChar_t:
                     b >> i->u_char;
                     break;
                  case kUShort_t:
                     b >> i->u_short;
                     break;
                  case kUInt_t:
                     b >> i->u_int;
                     break;
                  case kULong_t:
                     b >> i->u_long;
                     break;
                  case kULong64_t:
                     b >> i->u_longlong;
                     break;
                  case kDouble32_t:
                     b >> f;
                     i->dbl = double(f);
                     break;
                  case kchar:
                  case kNoType_t:
                  case kOther_t:
                     Error("TGenCollectionStreamer", "fType %d is not supported yet!\n", v->fKind);
                  }
               }
               break;
            case kIsClass:
               b.StreamObject(i, v->fType);
               break;
            case EProperty(kBIT_ISSTRING):
               i->read_std_string(b);
               break;
            case EProperty(kIsPointer | kIsClass):
               i->set(b.ReadObjectAny(v->fType));
               break;
            case EProperty(kIsPointer | kBIT_ISSTRING):
               i->read_std_string_pointer(b);
               break;
            case EProperty(kIsPointer | kBIT_ISTSTRING | kIsClass):
               i->read_tstring_pointer(vsn3, b);
               break;
         }
         v = fVal;
         addr += fValOffset;
      }
   }
   fFeed(fEnv->fStart,fEnv->fObject,fEnv->fSize);
   fDestruct(fEnv->fStart,fEnv->fSize);
   if (memory) {
      ::operator delete(memory);
   }
}

void TGenCollectionStreamer::WritePrimitives(int nElements, TBuffer &b)
{
   // Primitive output streamer.
   size_t len = fValDiff * nElements;
   char   buffer[8192];
   void*  memory  = 0;
   StreamHelper* itm = 0;
   switch (fSTL_type)  {
      case ROOT::kSTLvector:
         if (fVal->fKind != kBool_t)  {
            itm = (StreamHelper*)(fEnv->fStart = fFirst.invoke(fEnv));
            break;
         }
      default:
         fEnv->fStart = itm = (StreamHelper*)(len < sizeof(buffer) ? buffer : memory =::operator new(len));
         fCollect(fEnv->fObject,itm);
         break;
   }
   switch (int(fVal->fKind))   {
      case kBool_t:
         b.WriteFastArray(&itm->boolean    , nElements);
         break;
      case kChar_t:
         b.WriteFastArray(&itm->s_char    , nElements);
         break;
      case kShort_t:
         b.WriteFastArray(&itm->s_short   , nElements);
         break;
      case kInt_t:
         b.WriteFastArray(&itm->s_int     , nElements);
         break;
      case kLong_t:
         b.WriteFastArray(&itm->s_long    , nElements);
         break;
      case kLong64_t:
         b.WriteFastArray(&itm->s_longlong, nElements);
         break;
      case kFloat_t:
         b.WriteFastArray(&itm->flt       , nElements);
         break;
      case kFloat16_t:
         b.WriteFastArrayFloat16(&itm->flt, nElements);
         break;
      case kDouble_t:
         b.WriteFastArray(&itm->dbl       , nElements);
         break;
      case kUChar_t:
         b.WriteFastArray(&itm->u_char    , nElements);
         break;
      case kUShort_t:
         b.WriteFastArray(&itm->u_short   , nElements);
         break;
      case kUInt_t:
         b.WriteFastArray(&itm->u_int     , nElements);
         break;
      case kULong_t:
         b.WriteFastArray(&itm->u_long    , nElements);
         break;
      case kULong64_t:
         b.WriteFastArray(&itm->u_longlong, nElements);
         break;
      case kDouble32_t:
         b.WriteFastArrayDouble32(&itm->dbl, nElements);
         break;
      case kchar:
      case kNoType_t:
      case kOther_t:
         Error("TGenCollectionStreamer", "fType %d is not supported yet!\n", fVal->fKind);
   }
   if (memory)  {
      ::operator delete(memory);
   }
}

void TGenCollectionStreamer::WriteObjects(int nElements, TBuffer &b)
{
   // Object output streamer.
   StreamHelper* itm = 0;
   switch (fSTL_type)  {
         // Simple case: contiguous memory. get address of first, then jump.
      case ROOT::kSTLvector:
#define DOLOOP(x) {int idx=0; while(idx<nElements) {StreamHelper* i=(StreamHelper*)(((char*)itm) + fValDiff*idx); { x ;} ++idx;} break;}
         itm = (StreamHelper*)fFirst.invoke(fEnv);
         switch (fVal->fCase) {
            case kIsClass:
               DOLOOP(b.StreamObject(i, fVal->fType));
               break;
            case kBIT_ISSTRING:
               DOLOOP(TString(i->c_str()).Streamer(b));
               break;
            case kIsPointer | kIsClass:
               DOLOOP(b.WriteObjectAny(i->ptr(), fVal->fType));
               break;
            case kBIT_ISSTRING | kIsPointer:
               DOLOOP(i->write_std_string_pointer(b));
               break;
            case kBIT_ISTSTRING | kIsClass | kIsPointer:
               DOLOOP(i->write_tstring_pointer(b));
               break;
         }
#undef DOLOOP
         break;

         // No contiguous memory, but resize is possible
         // Hence accessing objects using At(i) should be not too much an overhead
      case ROOT::kSTLlist:
      case ROOT::kSTLforwardlist:
      case ROOT::kSTLdeque:
      case ROOT::kSTLmultiset:
      case ROOT::kSTLset:
      case ROOT::kSTLunorderedset:
      case ROOT::kSTLunorderedmultiset:
#define DOLOOP(x) {int idx=0; while(idx<nElements) {StreamHelper* i=(StreamHelper*)TGenCollectionProxy::At(idx); { x ;} ++idx;} break;}
         switch (fVal->fCase) {
            case kIsClass:
               DOLOOP(b.StreamObject(i, fVal->fType));
            case kBIT_ISSTRING:
               DOLOOP(TString(i->c_str()).Streamer(b));
            case kIsPointer | kIsClass:
               DOLOOP(b.WriteObjectAny(i->ptr(), fVal->fType));
            case kBIT_ISSTRING | kIsPointer:
               DOLOOP(i->write_std_string_pointer(b));
            case kBIT_ISTSTRING | kIsClass | kIsPointer:
               DOLOOP(i->write_tstring_pointer(b));
         }
#undef DOLOOP
         break;
      default:
         break;
   }
}

void TGenCollectionStreamer::WriteMap(int nElements, TBuffer &b)
{
   // Map output streamer
   StreamHelper* i;
   Value  *v;

   for (int loop, idx = 0; idx < nElements; ++idx)  {
      char* addr = (char*)TGenCollectionProxy::At(idx);
      v = fKey;
      for (loop = 0; loop < 2; ++loop)  {
         i = (StreamHelper*)addr;
         switch (v->fCase) {
            case kIsFundamental:  // Only handle primitives this way
            case kIsEnum:
               switch (int(v->fKind))   {
                  case kBool_t:
                     b << i->boolean;
                     break;
                  case kChar_t:
                     b << i->s_char;
                     break;
                  case kShort_t:
                     b << i->s_short;
                     break;
                  case kInt_t:
                     b << i->s_int;
                     break;
                  case kLong_t:
                     b << i->s_long;
                     break;
                  case kLong64_t:
                     b << i->s_longlong;
                     break;
                  case kFloat_t:
                     b << i->flt;
                     break;
                  case kFloat16_t:
                     b << float(i->flt);
                     break;
                  case kDouble_t:
                     b << i->dbl;
                     break;
                  case kUChar_t:
                     b << i->u_char;
                     break;
                  case kUShort_t:
                     b << i->u_short;
                     break;
                  case kUInt_t:
                     b << i->u_int;
                     break;
                  case kULong_t:
                     b << i->u_long;
                     break;
                  case kULong64_t:
                     b << i->u_longlong;
                     break;
                  case kDouble32_t:
                     b << float(i->dbl);
                     break;
                  case kchar:
                  case kNoType_t:
                  case kOther_t:
                     Error("TGenCollectionStreamer", "fType %d is not supported yet!\n", v->fKind);
               }
               break;
            case kIsClass:
               b.StreamObject(i, v->fType);
               break;
            case kBIT_ISSTRING:
               TString(i->c_str()).Streamer(b);
               break;
            case kIsPointer | kIsClass:
               b.WriteObjectAny(i->ptr(), v->fType);
               break;
            case kBIT_ISSTRING | kIsPointer:
               i->write_std_string_pointer(b);
               break;
            case kBIT_ISTSTRING | kIsClass | kIsPointer:
               i->write_tstring_pointer(b);
               break;
         }
         addr += fValOffset;
         v = fVal;
      }
   }
}

template <typename From, typename To>
void TGenCollectionStreamer::ConvertBufferVectorPrimitives(TBuffer &b, void *obj, Int_t nElements)
{
   From *temp = new From[nElements];
   b.ReadFastArray(temp, nElements);
   std::vector<To> *const vec = (std::vector<To>*)(obj);
   for(Int_t ind = 0; ind < nElements; ++ind) {
      (*vec)[ind] = (To)temp[ind];
   }
   delete [] temp;
}

template <typename To>
void TGenCollectionStreamer::ConvertBufferVectorPrimitivesFloat16(TBuffer &b, void *obj, Int_t nElements)
{
   Float16_t *temp = new Float16_t[nElements];
   b.ReadFastArrayFloat16(temp, nElements);
   std::vector<To> *const vec = (std::vector<To>*)(obj);
   for(Int_t ind = 0; ind < nElements; ++ind) {
      (*vec)[ind] = (To)temp[ind];
   }
   delete [] temp;
}

template <typename To>
void TGenCollectionStreamer::ConvertBufferVectorPrimitivesDouble32(TBuffer &b, void *obj, Int_t nElements)
{
   Double32_t *temp = new Double32_t[nElements];
   b.ReadFastArrayDouble32(temp, nElements);
   std::vector<To> *const vec = (std::vector<To>*)(obj);
   for(Int_t ind = 0; ind < nElements; ++ind) {
      (*vec)[ind] = (To)temp[ind];
   }
   delete [] temp;
}

template <typename To>
void TGenCollectionStreamer::DispatchConvertBufferVectorPrimitives(TBuffer &b, void *obj, Int_t nElements, const TVirtualCollectionProxy *onFileProxy)
{
   switch ((TStreamerInfo::EReadWrite)onFileProxy->GetType()) {
      case TStreamerInfo::kBool:     ConvertBufferVectorPrimitives<Bool_t    ,To>(b,obj,nElements); break;
      case TStreamerInfo::kChar:     ConvertBufferVectorPrimitives<Char_t    ,To>(b,obj,nElements); break;
      case TStreamerInfo::kShort:    ConvertBufferVectorPrimitives<Short_t   ,To>(b,obj,nElements); break;
      case TStreamerInfo::kInt:      ConvertBufferVectorPrimitives<Int_t     ,To>(b,obj,nElements); break;
      case TStreamerInfo::kLong:     ConvertBufferVectorPrimitives<Long_t    ,To>(b,obj,nElements); break;
      case TStreamerInfo::kLong64:   ConvertBufferVectorPrimitives<Long64_t  ,To>(b,obj,nElements); break;
      case TStreamerInfo::kFloat:    ConvertBufferVectorPrimitives<Float_t   ,To>(b,obj,nElements); break;
      case TStreamerInfo::kFloat16:  ConvertBufferVectorPrimitives<Float16_t ,To>(b,obj,nElements); break;
      case TStreamerInfo::kDouble:   ConvertBufferVectorPrimitives<Double_t  ,To>(b,obj,nElements); break;
      case TStreamerInfo::kDouble32: ConvertBufferVectorPrimitives<Double32_t,To>(b,obj,nElements); break;
      case TStreamerInfo::kUChar:    ConvertBufferVectorPrimitives<UChar_t   ,To>(b,obj,nElements); break;
      case TStreamerInfo::kUShort:   ConvertBufferVectorPrimitives<UShort_t  ,To>(b,obj,nElements); break;
      case TStreamerInfo::kUInt:     ConvertBufferVectorPrimitives<UInt_t    ,To>(b,obj,nElements); break;
      case TStreamerInfo::kULong:    ConvertBufferVectorPrimitives<ULong_t   ,To>(b,obj,nElements); break;
      case TStreamerInfo::kULong64:  ConvertBufferVectorPrimitives<ULong64_t ,To>(b,obj,nElements); break;
     default: break;
   }
}

template <typename basictype>
void TGenCollectionStreamer::ReadBufferVectorPrimitives(TBuffer &b, void *obj, const TClass *onFileClass)
{
   int nElements = 0;
   b >> nElements;
   fResize(obj,nElements);

   if (onFileClass) {
      DispatchConvertBufferVectorPrimitives<basictype>(b,obj,nElements,onFileClass->GetCollectionProxy());
   } else {
      TVirtualVectorIterators iterators(fFunctionCreateIterators);
      iterators.CreateIterators(obj);
      b.ReadFastArray((basictype*)iterators.fBegin, nElements);
   }
}

void TGenCollectionStreamer::ReadBufferVectorPrimitivesFloat16(TBuffer &b, void *obj, const TClass *onFileClass)
{
   int nElements = 0;
   b >> nElements;
   fResize(obj,nElements);

   if (onFileClass) {
      DispatchConvertBufferVectorPrimitives<Float16_t>(b,obj,nElements,onFileClass->GetCollectionProxy());
   } else {
      TVirtualVectorIterators iterators(fFunctionCreateIterators);
      iterators.CreateIterators(obj);
      b.ReadFastArrayFloat16((Float16_t*)iterators.fBegin, nElements);
   }
}

void TGenCollectionStreamer::ReadBufferVectorPrimitivesDouble32(TBuffer &b, void *obj, const TClass *onFileClass)
{
   int nElements = 0;
   b >> nElements;
   fResize(obj,nElements);

   if (onFileClass) {
      DispatchConvertBufferVectorPrimitives<Double32_t>(b,obj,nElements,onFileClass->GetCollectionProxy());
   } else {
      TVirtualVectorIterators iterators(fFunctionCreateIterators);
      iterators.CreateIterators(obj);
      b.ReadFastArrayDouble32((Double32_t*)iterators.fBegin, nElements);
   }
}



void TGenCollectionStreamer::ReadBuffer(TBuffer &b, void *obj, const TClass *onFileClass)
{
   // Call the specialized function.  The first time this call ReadBufferDefault which
   // actually set to fReadBufferFunc to the 'right' specialized version.

   (this->*fReadBufferFunc)(b,obj,onFileClass);
}

void TGenCollectionStreamer::ReadBuffer(TBuffer &b, void *obj)
{
   // Call the specialized function.  The first time this call ReadBufferDefault which
   // actually set to fReadBufferFunc to the 'right' specialized version.

   (this->*fReadBufferFunc)(b,obj,0);
}

void TGenCollectionStreamer::ReadBufferDefault(TBuffer &b, void *obj, const TClass *onFileClass)
{

   fReadBufferFunc = &TGenCollectionStreamer::ReadBufferGeneric;

   // We will need this later, so let's make sure it is initialized.
   if ( !fValue.load() ) InitializeEx(kFALSE);
   if (!GetFunctionCreateIterators()) {
      Fatal("TGenCollectionStreamer::ReadBufferDefault","No CreateIterators function for %s",fName.c_str());
   }
   if (fSTL_type == ROOT::kSTLvector && ( fVal->fCase == kIsFundamental || fVal->fCase == kIsEnum ) )
   {
      // Only handle primitives this way
      switch (int(fVal->fKind))   {
         case kBool_t:
            // Nothing use generic for now
            break;
         case kChar_t:
            fReadBufferFunc = &TGenCollectionStreamer::ReadBufferVectorPrimitives<Char_t>;
            break;
         case kShort_t:
            fReadBufferFunc = &TGenCollectionStreamer::ReadBufferVectorPrimitives<Short_t>;
            break;
         case kInt_t:
            fReadBufferFunc = &TGenCollectionStreamer::ReadBufferVectorPrimitives<Int_t>;
            break;
         case kLong_t:
            fReadBufferFunc = &TGenCollectionStreamer::ReadBufferVectorPrimitives<Long_t>;
            break;
         case kLong64_t:
            fReadBufferFunc = &TGenCollectionStreamer::ReadBufferVectorPrimitives<Long64_t>;
            break;
         case kFloat_t:
            fReadBufferFunc = &TGenCollectionStreamer::ReadBufferVectorPrimitives<Float_t>;
            break;
         case kFloat16_t:
            fReadBufferFunc = &TGenCollectionStreamer::ReadBufferVectorPrimitivesFloat16;
            break;
         case kDouble_t:
            fReadBufferFunc = &TGenCollectionStreamer::ReadBufferVectorPrimitives<Double_t>;
            break;
         case kUChar_t:
            fReadBufferFunc = &TGenCollectionStreamer::ReadBufferVectorPrimitives<UChar_t>;
            break;
         case kUShort_t:
            fReadBufferFunc = &TGenCollectionStreamer::ReadBufferVectorPrimitives<UShort_t>;
            break;
         case kUInt_t:
            fReadBufferFunc = &TGenCollectionStreamer::ReadBufferVectorPrimitives<UInt_t>;
            break;
         case kULong_t:
            fReadBufferFunc = &TGenCollectionStreamer::ReadBufferVectorPrimitives<ULong_t>;
            break;
         case kULong64_t:
            fReadBufferFunc = &TGenCollectionStreamer::ReadBufferVectorPrimitives<ULong64_t>;
            break;
         case kDouble32_t:
            fReadBufferFunc = &TGenCollectionStreamer::ReadBufferVectorPrimitivesDouble32;
            break;
         case kchar:
         case kNoType_t:
         case kOther_t:
            // Nothing use the generic for now
            break;
      }
   }
   (this->*fReadBufferFunc)(b,obj,onFileClass);
}

void TGenCollectionStreamer::ReadBufferGeneric(TBuffer &b, void *obj, const TClass *onFileClass)
{
   TVirtualCollectionProxy::TPushPop env(this, obj);

   int nElements = 0;
   b >> nElements;

   if (nElements == 0) {
      if (obj) {
         TGenCollectionProxy::Clear("force");
      }
   } else if (nElements > 0)  {
      switch (fSTL_type)  {
         case ROOT::kSTLbitset:
            if (obj) {
               if (fProperties & kNeedDelete)   {
                  TGenCollectionProxy::Clear("force");
               }  else {
                  fClear.invoke(fEnv);
               }
            }
            ReadPrimitives(nElements, b, onFileClass);
            return;
         case ROOT::kSTLvector:
            if (obj) {
               if (fProperties & kNeedDelete)   {
                  TGenCollectionProxy::Clear("force");
               } // a resize will be called in ReadPrimitives/ReadObjects.
               else if (fVal->fKind == kBool_t) {
                  fClear.invoke(fEnv);
               }
            }
            switch (fVal->fCase) {
               case kIsFundamental:  // Only handle primitives this way
               case kIsEnum:
                  ReadPrimitives(nElements, b, onFileClass);
                  return;
               default:
                  ReadObjects(nElements, b, onFileClass);
                  return;
            }
            break;
         case ROOT::kSTLlist:
         case ROOT::kSTLforwardlist:
         case ROOT::kSTLdeque:
         case ROOT::kSTLmultiset:
         case ROOT::kSTLset:
         case ROOT::kSTLunorderedset:
         case ROOT::kSTLunorderedmultiset:
            if (obj) {
               if (fProperties & kNeedDelete)   {
                  TGenCollectionProxy::Clear("force");
               }  else {
                  fClear.invoke(fEnv);
               }
            }
            switch (fVal->fCase) {
               case kIsFundamental:  // Only handle primitives this way
               case kIsEnum:
                  ReadPrimitives(nElements, b, onFileClass);
                  return;
               default:
                  ReadObjects(nElements, b, onFileClass);
                  return;
            }
            break;
         case ROOT::kSTLmap:
         case ROOT::kSTLmultimap:
         case ROOT::kSTLunorderedmap:
         case ROOT::kSTLunorderedmultimap:
            if (obj) {
               if (fProperties & kNeedDelete)   {
                  TGenCollectionProxy::Clear("force");
               }  else {
                  fClear.invoke(fEnv);
               }
            }
            ReadMap(nElements, b, onFileClass);
            break;
      }
   }
}

void TGenCollectionStreamer::Streamer(TBuffer &b)
{
   // TClassStreamer IO overload.
   if (b.IsReading()) {    //Read mode
      int nElements = 0;
      b >> nElements;
      if (fEnv->fObject)   {
         TGenCollectionProxy::Clear("force");
      }
      if (nElements > 0)  {
         switch (fSTL_type)  {
            case ROOT::kSTLbitset:
               ReadPrimitives(nElements, b, fOnFileClass);
               return;
            case ROOT::kSTLvector:
            case ROOT::kSTLlist:
            case ROOT::kSTLdeque:
            case ROOT::kSTLmultiset:
            case ROOT::kSTLset:
            case ROOT::kSTLunorderedset:
            case ROOT::kSTLunorderedmultiset:
               switch (fVal->fCase) {
                  case kIsFundamental:  // Only handle primitives this way
                  case kIsEnum:
                     ReadPrimitives(nElements, b, fOnFileClass);
                     return;
                  default:
                     ReadObjects(nElements, b, fOnFileClass);
                     return;
               }
               break;
            case ROOT::kSTLmap:
            case ROOT::kSTLmultimap:
            case ROOT::kSTLunorderedmap:
            case ROOT::kSTLunorderedmultimap:
               ReadMap(nElements, b, fOnFileClass);
               break;
         }
      }
   } else {    // Write case
      int nElements = fEnv->fObject ? *(size_t*)fSize.invoke(fEnv) : 0;
      b << nElements;
      if (nElements > 0)  {
         switch (fSTL_type)  {
            case ROOT::kSTLbitset:
               WritePrimitives(nElements, b);
               return;
            case ROOT::kSTLvector:
            case ROOT::kSTLlist:
            case ROOT::kSTLforwardlist:
            case ROOT::kSTLdeque:
            case ROOT::kSTLmultiset:
            case ROOT::kSTLset:
            case ROOT::kSTLunorderedset:
            case ROOT::kSTLunorderedmultiset:
               switch (fVal->fCase) {
                  case kIsFundamental:  // Only handle primitives this way
                  case kIsEnum:
                     WritePrimitives(nElements, b);
                     return;
                  default:
                     WriteObjects(nElements, b);
                     return;
               }
               break;
            case ROOT::kSTLmap:
            case ROOT::kSTLmultimap:
            case ROOT::kSTLunorderedmap:
            case ROOT::kSTLunorderedmultimap:
               WriteMap(nElements, b);
               break;
         }
      }
   }
}

void TGenCollectionStreamer::StreamerAsMap(TBuffer &b)
{
   // TClassStreamer IO overload.
   if (b.IsReading()) {    //Read mode
      int nElements = 0;
      b >> nElements;
      if (fEnv->fObject)   {
         TGenCollectionProxy::Clear("force");
      }
      if (nElements > 0)  {
         switch (fSTL_type)  {
            case ROOT::kSTLmap:
            case ROOT::kSTLmultimap:
            case ROOT::kSTLunorderedmap:
            case ROOT::kSTLunorderedmultimap:
               ReadMap(nElements, b, fOnFileClass);
               break;
            case ROOT::kSTLvector:
            case ROOT::kSTLlist:
            case ROOT::kSTLforwardlist:
            case ROOT::kSTLdeque:
            case ROOT::kSTLmultiset:
            case ROOT::kSTLset:
            case ROOT::kSTLunorderedset:
            case ROOT::kSTLunorderedmultiset:{
                  ReadPairFromMap(nElements, b);
                  break;
               }
            default:
               break;
         }
      }
   } else {    // Write case
      Streamer(b);
   }
}
