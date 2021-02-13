// $Id$
// Author: Sergey Linev  21.12.2017

/*************************************************************************
 * Copyright (C) 1995-2017, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TBufferText
#define ROOT_TBufferText

#include "TBufferIO.h"

class TStreamerBase;
class TExMap;

class TBufferText : public TBufferIO {

protected:
   TBufferText();
   TBufferText(TBuffer::EMode mode, TObject *parent = nullptr);

public:
   virtual ~TBufferText();

   // virtual TBuffer methods, which are generic for all text-based streamers

   void StreamObject(void *obj, const std::type_info &typeinfo, const TClass *onFileClass = nullptr) override;
   void StreamObject(void *obj, const char *className, const TClass *onFileClass = nullptr) override;
   void StreamObject(TObject *obj) override;
   using TBuffer::StreamObject;

   Int_t ApplySequence(const TStreamerInfoActions::TActionSequence &sequence, void *object) final;
   Int_t ApplySequenceVecPtr(const TStreamerInfoActions::TActionSequence &sequence, void *start_collection,
                                     void *end_collection) final;
   Int_t
   ApplySequence(const TStreamerInfoActions::TActionSequence &sequence, void *start_collection, void *end_collection) final;

   void ReadFloat16(Float_t *f, TStreamerElement *ele = nullptr) final;
   void WriteFloat16(Float_t *f, TStreamerElement *ele = nullptr) final;
   void ReadDouble32(Double_t *d, TStreamerElement *ele = nullptr) final;
   void WriteDouble32(Double_t *d, TStreamerElement *ele = nullptr) final;
   void ReadWithFactor(Float_t *ptr, Double_t factor, Double_t minvalue) final;
   void ReadWithNbits(Float_t *ptr, Int_t nbits) final;
   void ReadWithFactor(Double_t *ptr, Double_t factor, Double_t minvalue) final;
   void ReadWithNbits(Double_t *ptr, Int_t nbits) final;

   Int_t ReadArrayFloat16(Float_t *&f, TStreamerElement *ele = nullptr) override;
   Int_t ReadArrayDouble32(Double_t *&d, TStreamerElement *ele = nullptr) override;

   Int_t ReadStaticArrayFloat16(Float_t *f, TStreamerElement *ele = nullptr) final;
   Int_t ReadStaticArrayDouble32(Double_t *d, TStreamerElement *ele = nullptr) final;

   void ReadFastArrayFloat16(Float_t *f, Int_t n, TStreamerElement *ele = nullptr) final;
   void ReadFastArrayDouble32(Double_t *d, Int_t n, TStreamerElement *ele = nullptr) final;
   void ReadFastArrayWithFactor(Float_t *ptr, Int_t n, Double_t factor, Double_t minvalue) final;
   void ReadFastArrayWithNbits(Float_t *ptr, Int_t n, Int_t nbits) final;
   void ReadFastArrayWithFactor(Double_t *ptr, Int_t n, Double_t factor, Double_t minvalue) final;
   void ReadFastArrayWithNbits(Double_t *ptr, Int_t n, Int_t nbits) final;

   void WriteArrayFloat16(const Float_t *f, Int_t n, TStreamerElement *ele = nullptr) final;
   void WriteArrayDouble32(const Double_t *d, Int_t n, TStreamerElement *ele = nullptr) final;

   void WriteFastArrayFloat16(const Float_t *d, Int_t n, TStreamerElement *ele = nullptr) final;
   void WriteFastArrayDouble32(const Double_t *d, Int_t n, TStreamerElement *ele = nullptr) final;

   // Utilities for TClass
   Int_t ReadClassBuffer(const TClass * /*cl*/, void * /*pointer*/, const TClass * /*onfile_class*/ = nullptr) override;
   Int_t ReadClassBuffer(const TClass * /*cl*/, void * /*pointer*/, Int_t /*version*/, UInt_t /*start*/,
                                 UInt_t /*count*/, const TClass * /*onfile_class*/ = nullptr) override;
   Int_t WriteClassBuffer(const TClass *cl, void *pointer) override;

   // virtual abstract TBuffer methods, which are not used in text streaming

   Int_t CheckByteCount(UInt_t /* startpos */, UInt_t /* bcnt */, const TClass * /* clss */) final { return 0; }
   Int_t CheckByteCount(UInt_t /* startpos */, UInt_t /* bcnt */, const char * /* classname */) final { return 0; }
   void SetByteCount(UInt_t /* cntpos */, Bool_t /* packInVersion */ = kFALSE) final {}
   void SkipVersion(const TClass *cl = nullptr) final;
   Version_t ReadVersionNoCheckSum(UInt_t *, UInt_t *) final { return 0; }

   Int_t ReadBuf(void * /*buf*/, Int_t /*max*/) final
   {
      Error("ReadBuf", "useless in text streamers");
      return 0;
   }
   void WriteBuf(const void * /*buf*/, Int_t /*max*/) final { Error("WriteBuf", "useless in text streamers"); }

   char *ReadString(char * /*s*/, Int_t /*max*/) final
   {
      Error("ReadString", "useless");
      return nullptr;
   }
   void WriteString(const char * /*s*/) final { Error("WriteString", "useless"); }

   Version_t ReadVersionForMemberWise(const TClass * /*cl*/ = nullptr) final
   {
      Error("ReadVersionForMemberWise", "not defined in text-based streamers");
      return 0;
   }
   UInt_t WriteVersionMemberWise(const TClass * /*cl*/, Bool_t /*useBcnt*/ = kFALSE) final
   {
      Error("WriteVersionMemberWise", "not defined in text-based streamers");
      return 0;
   }

   TObject *ReadObject(const TClass * /*cl*/) final
   {
      Error("ReadObject", "not yet implemented for text-based streamers");
      return nullptr;
   }

   // Utilities for TClass
   Int_t ReadClassEmulated(const TClass * /*cl*/, void * /*object*/, const TClass * /*onfile_class*/ = nullptr) final
   {
      Error("ReadClassEmulated", "not defined in text-based streamers");
      return 0;
   }

   virtual void WriteBaseClass(void *start, TStreamerBase *elem);

   virtual void ReadBaseClass(void *start, TStreamerBase *elem);

   static void SetFloatFormat(const char *fmt = "%e");
   static const char *GetFloatFormat();
   static void SetDoubleFormat(const char *fmt = "%.14e");
   static const char *GetDoubleFormat();

   static void CompactFloatString(char *buf, unsigned len);
   static const char *ConvertFloat(Float_t v, char *buf, unsigned len, Bool_t not_optimize = kFALSE);
   static const char *ConvertDouble(Double_t v, char *buf, unsigned len, Bool_t not_optimize = kFALSE);

protected:
   static const char *fgFloatFmt;  ///<!  printf argument for floats, either "%f" or "%e" or "%10f" and so on
   static const char *fgDoubleFmt; ///<!  printf argument for doubles, either "%f" or "%e" or "%10f" and so on

   ClassDefOverride(TBufferText, 0); // a TBuffer subclass for all text-based streamers
};

#endif
