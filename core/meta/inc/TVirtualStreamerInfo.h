// @(#)root/meta:$Id$
// Author: Rene Brun   05/02/2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TVirtualStreamerInfo
#define ROOT_TVirtualStreamerInfo


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TVirtualStreamerInfo   Abstract Interface class                      //
//                                                                      //
// Abstract Interface describing Streamer information for one class.    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TNamed.h"

#include "ESTLType.h"

class TFile;
class TClass;
class TObjArray;
class TStreamerElement;
class TStreamerBasicType;
class TVirtualCollectionProxy;
class TClassStreamer;
namespace ROOT { namespace Detail { class TCollectionProxyInfo; } }

class TVirtualStreamerInfo : public TNamed {

protected:
   Bool_t              fOptimized : 1;     //! true if the StreamerInfo has been optimized
   Bool_t              fIsBuilt : 1;       //! true if the StreamerInfo has been 'built' (i.e. has all the StreamerElements it should have)
   std::atomic<Bool_t> fIsCompiled;        //! true if the StreamerInfo has been compiled (i.e. fully built, ready to use for streaming).

protected:
   static  Bool_t    fgCanDelete;        //True if ReadBuffer can delete object
   static  Bool_t    fgOptimize;         //True if optimization on
   static  Bool_t    fgStreamMemberWise; //True if the collections are to be stream "member-wise" (when possible).
   static TVirtualStreamerInfo  *fgInfoFactory;

   TVirtualStreamerInfo(const TVirtualStreamerInfo& info);
   TVirtualStreamerInfo& operator=(const TVirtualStreamerInfo&);

   void  ResetIsCompiled() {
      fIsCompiled = kFALSE;
      ResetBit(kIsCompiled); /* for backward compatibility */
   }
   void  SetIsCompiled() {
      fIsCompiled = kTRUE;
      SetBit(kIsCompiled); /* for backward compatibility */
   }

public:

   //status bits
   enum EStatusBits {
      kCannotOptimize        = BIT(12),
      kIgnoreTObjectStreamer = BIT(13),  // eventhough BIT(13) is taken up by TObject (to preserve forward compatibility)
      kRecovered             = BIT(14),
      kNeedCheck             = BIT(15),
      kIsCompiled            = BIT(16),
      kBuildOldUsed          = BIT(17),
      kBuildRunning          = BIT(18)
   };

   enum class EStatusBitsDupExceptions {
      // This bit duplicates TObject::kInvalidObject. As the semantic of kIgnoreTObjectStreamer is a persistent,
      // we can not change its value without breaking forward compatibility.
      // Furthermore, TObject::kInvalidObject and its semantic is not (and should not be)
      // used in TVirtualStreamerInfo
      kIgnoreTObjectStreamer  = TVirtualStreamerInfo::kIgnoreTObjectStreamer,
   };

   enum EReadWrite {
      kBase        =  0,  kOffsetL = 20,  kOffsetP = 40,  kCounter =  6,  kCharStar = 7,
      kChar        =  1,  kShort   =  2,  kInt     =  3,  kLong    =  4,  kFloat    = 5,
      kDouble      =  8,  kDouble32=  9,
      kUChar       = 11,  kUShort  = 12,  kUInt    = 13,  kULong   = 14,  kBits     = 15,
      kLong64      = 16,  kULong64 = 17,  kBool    = 18,  kFloat16 = 19,
      kObject      = 61,  kAny     = 62,  kObjectp = 63,  kObjectP = 64,  kTString  = 65,
      kTObject     = 66,  kTNamed  = 67,  kAnyp    = 68,  kAnyP    = 69,  kAnyPnoVT = 70,
      kSTLp        = 71,
      kSkip        = 100, kSkipL = 120, kSkipP   = 140,
      kConv        = 200, kConvL = 220, kConvP   = 240,
      kSTL         = ROOT::kSTLany /* 300 */,
      kSTLstring   = ROOT::kSTLstring /* 365 */,
      kStreamer    = 500, kStreamLoop = 501,
      kCache       = 600,  // Cache the value in memory than is not part of the object but is accessible via a SchemaRule
      kArtificial  = 1000,
      kCacheNew    = 1001,
      kCacheDelete = 1002,
      kNeedObjectForVirtualBaseClass = 99997,
      kMissing     = 99999
   };

// Some comments about EReadWrite
// kBase    : base class element
// kOffsetL : fixed size array
// kOffsetP : pointer to object
// kCounter : counter for array size
// kCharStar: pointer to array of char
// kBits    : TObject::fBits in case of a referenced object
// kObject  : Class  derived from TObject
// kObjectp : Class* derived from TObject and with    comment field //->Class
// kObjectP : Class* derived from TObject and with NO comment field //->Class
// kAny     : Class  not derived from TObject
// kAnyp    : Class* not derived from TObject with    comment field //->Class
// kAnyP    : Class* not derived from TObject with NO comment field //->Class
// kAnyPnoVT: Class* not derived from TObject with NO comment field //->Class and Class has NO virtual table
// kSTLp    : Pointer to STL container.
// kTString : TString, special case
// kTObject : TObject, special case
// kTNamed  : TNamed , special case



   TVirtualStreamerInfo();
   TVirtualStreamerInfo(TClass * /*cl*/);
   virtual            ~TVirtualStreamerInfo();
   virtual void        Build(Bool_t isTransient = kFALSE) = 0;
   virtual void        BuildCheck(TFile *file = 0, Bool_t load = kTRUE) = 0;
   virtual void        BuildEmulated(TFile *file) = 0;
   virtual void        BuildOld() = 0;
   virtual Bool_t      BuildFor( const TClass *cl ) = 0;
   virtual void        CallShowMembers(const void* obj, TMemberInspector &insp, Bool_t isTransient) const = 0;
   virtual void        Clear(Option_t *) = 0;
   virtual Bool_t      CompareContent(TClass *cl,TVirtualStreamerInfo *info, Bool_t warn, Bool_t complete, TFile *file) = 0;
   virtual void        Compile() = 0;
   virtual void        ForceWriteInfo(TFile *file, Bool_t force=kFALSE) = 0;
   virtual Int_t       GenerateHeaderFile(const char *dirname, const TList *subClasses = 0, const TList *extrainfos = 0) = 0;
   virtual TClass     *GetActualClass(const void *obj) const = 0;
   virtual TClass     *GetClass() const  = 0;
   virtual UInt_t      GetCheckSum() const = 0;
   virtual Int_t       GetClassVersion() const = 0;
   virtual TStreamerElement *GetElem(Int_t id) const = 0;
   virtual TStreamerElement *GetElement(Int_t id) const = 0;
   virtual TObjArray  *GetElements() const = 0;
   virtual Int_t       GetOffset(const char *) const = 0;
   virtual Int_t       GetOffset(Int_t id) const = 0;
   virtual Int_t       GetElementOffset(Int_t id) const = 0;
   virtual Version_t   GetOldVersion() const = 0;
   virtual Int_t       GetOnFileClassVersion() const = 0;
   virtual Int_t       GetNumber()  const = 0;
   virtual Int_t       GetSize()    const = 0;
   virtual TStreamerElement *GetStreamerElement(const char*datamember, Int_t& offset) const = 0;
           Bool_t      IsBuilt() const { return fIsBuilt; }
           Bool_t      IsCompiled() const { return fIsCompiled; }
           Bool_t      IsOptimized() const { return fOptimized; }
           Int_t       IsRecovered() const { return TestBit(kRecovered); }
   virtual void        ls(Option_t *option="") const = 0;
   virtual TVirtualStreamerInfo *NewInfo(TClass *cl) = 0;
   virtual void       *New(void *obj = 0) = 0;
   virtual void       *NewArray(Long_t nElements, void* ary = 0) = 0;
   virtual void        Destructor(void* p, Bool_t dtorOnly = kFALSE) = 0;
   virtual void        DeleteArray(void* p, Bool_t dtorOnly = kFALSE) = 0;

   virtual void        SetCheckSum(UInt_t checksum) = 0;
   virtual void        SetClass(TClass *cl) = 0;
   virtual void        SetClassVersion(Int_t vers) = 0;
   static  Bool_t      SetStreamMemberWise(Bool_t enable = kTRUE);
   virtual void        TagFile(TFile *fFile) = 0;
   virtual void        Update(const TClass *oldClass, TClass *newClass) = 0;

   static const char         *GetElementCounterStart(const char *dmTitle);
   static TStreamerBasicType *GetElementCounter(const char *countName, TClass *cl);

   static Bool_t       CanOptimize();
   static Bool_t       GetStreamMemberWise();
   static void         Optimize(Bool_t opt=kTRUE);
   static Bool_t       CanDelete();
   static void         SetCanDelete(Bool_t opt=kTRUE);
   static void         SetFactory(TVirtualStreamerInfo *factory);

   // \brief Generate the TClass and TStreamerInfo for the requested pair.
   // This creates a TVirtualStreamerInfo for the pair and trigger the BuildCheck/Old to
   // provokes the creation of the corresponding TClass.  This relies on the dictionary for
   // std::pair<const int, int> to already exist (or the interpreter information being available)
   // as it is used as a template.
   virtual TVirtualStreamerInfo *GenerateInfoForPair(const std::string &pairclassname, bool silent, size_t hint_pair_offset, size_t hint_pair_size) = 0;
   virtual TVirtualStreamerInfo *GenerateInfoForPair(const std::string &firstname, const std::string &secondname, bool silent, size_t hint_pair_offset, size_t hint_pair_size) = 0;

   virtual TVirtualCollectionProxy *GenEmulatedProxy(const char* class_name, Bool_t silent) = 0;
   virtual TClassStreamer *GenEmulatedClassStreamer(const char* class_name, Bool_t silent) = 0;
   virtual TVirtualCollectionProxy *GenExplicitProxy( const ::ROOT::Detail::TCollectionProxyInfo &info, TClass *cl ) = 0;
   virtual TClassStreamer *GenExplicitClassStreamer( const ::ROOT::Detail::TCollectionProxyInfo &info, TClass *cl ) = 0;
   static TVirtualStreamerInfo *Factory();

   //WARNING this class version must be the same as TStreamerInfo
   ClassDef(TVirtualStreamerInfo,6)  //Abstract Interface describing Streamer information for one class
};

#endif
