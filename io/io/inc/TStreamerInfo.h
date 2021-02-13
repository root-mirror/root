// @(#)root/io:$Id$
// Author: Rene Brun   12/10/2000

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TStreamerInfo
#define ROOT_TStreamerInfo

#include <atomic>
#include <vector>

#include "TVirtualStreamerInfo.h"

#include "TVirtualCollectionProxy.h"

#include "TObjArray.h"

/**
\class TStreamerInfo
\ingroup IO
Describe Streamer information for one class version
*/

class TFile;
class TClass;
class TClonesArray;
class TDataMember;
class TMemberStreamer;
class TStreamerElement;
class TStreamerBasicType;
class TClassStreamer;
class TVirtualArray;
namespace ROOT { namespace Detail { class TCollectionProxyInfo; } }
namespace ROOT { class TSchemaRule; }

namespace TStreamerInfoActions { class TActionSequence; }

class TStreamerInfo : public TVirtualStreamerInfo {

   class TCompInfo {
   // Class used to cache information (see fComp)
   private:
      // TCompInfo(const TCompInfo&) = default;
      // TCompInfo& operator=(const TCompInfo&) = default;
   public:
      Int_t             fType;
      Int_t             fNewType;
      Int_t             fOffset;
      Int_t             fLength;
      TStreamerElement *fElem;     ///< Not Owned
      ULong_t           fMethod;
      TClass           *fClass;    ///< Not Owned
      TClass           *fNewClass; ///< Not Owned
      TString           fClassName;
      TMemberStreamer  *fStreamer; ///< Not Owned
      TCompInfo() : fType(-1), fNewType(0), fOffset(0), fLength(0), fElem(0), fMethod(0),
                    fClass(0), fNewClass(0), fClassName(), fStreamer(0) {};
      ~TCompInfo() {};
      void Update(const TClass *oldcl, TClass *newcl);
   };
   friend class TStreamerInfoActions::TActionSequence;

public:
   // make the opaque pointer public.
   typedef TCompInfo TCompInfo_t;

protected:
   //---------------------------------------------------------------------------
   // Adapter class used to handle streaming collection of pointers
   //---------------------------------------------------------------------------
   class TPointerCollectionAdapter
   {
   public:
      TPointerCollectionAdapter( TVirtualCollectionProxy *proxy ):
         fProxy( proxy ) {}

      char* operator[]( UInt_t idx ) const
      {
         char **el = (char**)fProxy->At(idx);
         return *el;
      }
   private:
      TVirtualCollectionProxy *fProxy;
   };

private:
   UInt_t            fCheckSum;          ///<Checksum of original class
   Int_t             fClassVersion;      ///<Class version identifier
   Int_t             fOnFileClassVersion;///<!Class version identifier as stored on file.
   Int_t             fNumber;            ///<!Unique identifier
   Int_t             fSize;              ///<!size of the persistent class
   Int_t             fNdata;             ///<!number of optimized elements
   Int_t             fNfulldata;         ///<!number of elements
   Int_t             fNslots;            ///<!total numbrer of slots in fComp.
   TCompInfo        *fComp;              ///<![fNslots with less than fElements->GetEntries()*1.5 used] Compiled info
   TCompInfo       **fCompOpt;           ///<![fNdata]
   TCompInfo       **fCompFull;          ///<![fElements->GetEntries()]
   TClass           *fClass;             ///<!pointer to class
   TObjArray        *fElements;          ///<Array of TStreamerElements
   Version_t         fOldVersion;        ///<! Version of the TStreamerInfo object read from the file
   Int_t             fNVirtualInfoLoc;   ///<! Number of virtual info location to update.
   ULong_t          *fVirtualInfoLoc;    ///<![fNVirtualInfoLoc] Location of the pointer to the TStreamerInfo inside the object (when emulated)
   TStreamerInfoActions::TActionSequence *fReadObjectWise;        ///<! List of read action resulting from the compilation.
   TStreamerInfoActions::TActionSequence *fReadMemberWise;        ///<! List of read action resulting from the compilation for use in member wise streaming.
   TStreamerInfoActions::TActionSequence *fReadMemberWiseVecPtr;  ///<! List of read action resulting from the compilation for use in member wise streaming.
   TStreamerInfoActions::TActionSequence *fReadText;              ///<! List of text read action resulting from the compilation, used for JSON.
   TStreamerInfoActions::TActionSequence *fWriteObjectWise;       ///<! List of write action resulting from the compilation.
   TStreamerInfoActions::TActionSequence *fWriteMemberWise;       ///<! List of write action resulting from the compilation for use in member wise streaming.
   TStreamerInfoActions::TActionSequence *fWriteMemberWiseVecPtr; ///<! List of write action resulting from the compilation for use in member wise streaming.
   TStreamerInfoActions::TActionSequence *fWriteText;             ///<! List of text write action resulting for the compilation, used for JSON.

   static std::atomic<Int_t>             fgCount;     ///<Number of TStreamerInfo instances

   template <typename T> static T GetTypedValueAux(Int_t type, void *ladd, int k, Int_t len);
   static void       PrintValueAux(char *ladd, Int_t atype, TStreamerElement * aElement, Int_t aleng, Int_t *count);

   UInt_t            GenerateIncludes(FILE *fp, char *inclist, const TList *extrainfos);
   void              GenerateDeclaration(FILE *fp, FILE *sfp, const TList *subClasses, Bool_t top = kTRUE);
   void              InsertArtificialElements(std::vector<const ROOT::TSchemaRule*> &rules);
   void              DestructorImpl(void* p, Bool_t dtorOnly);

private:
   TStreamerInfo(const TStreamerInfo&) = delete;            // TStreamerInfo are not copiable.  Not Implemented.
   TStreamerInfo& operator=(const TStreamerInfo&) = delete; // TStreamerInfo are not copiable.  Not Implemented.
   void AddReadAction(TStreamerInfoActions::TActionSequence *readSequence, Int_t index, TCompInfo *compinfo);
   void AddWriteAction(TStreamerInfoActions::TActionSequence *writeSequence, Int_t index, TCompInfo *compinfo);
   void AddReadTextAction(TStreamerInfoActions::TActionSequence *readSequence, Int_t index, TCompInfo *compinfo);
   void AddWriteTextAction(TStreamerInfoActions::TActionSequence *writeSequence, Int_t index, TCompInfo *compinfo);
   void AddReadMemberWiseVecPtrAction(TStreamerInfoActions::TActionSequence *readSequence, Int_t index, TCompInfo *compinfo);
   void AddWriteMemberWiseVecPtrAction(TStreamerInfoActions::TActionSequence *writeSequence, Int_t index, TCompInfo *compinfo);

public:

   /// Status bits
   /// See TVirtualStreamerInfo::EStatusBits for the values.

   /// EReadWrite Enumerator
   /// | Enum Constant | Description   |
   /// |-------------|--------------------|
   /// | kBase       | Base class element |
   /// | kOffsetL    | Fixed size array |
   /// | kOffsetP    | Pointer to object |
   /// | kCounter    | Counter for array size |
   /// | kCharStar   | Pointer to array of char |
   /// | kLegacyChar | Equal to TDataType's kchar |
   /// | kBits       | TObject::fBits in case of a referenced object |
   /// | kObject     | Class  derived from TObject, or for TStreamerSTL::fCtype non-pointer elements |
   /// | kObjectp    | Class* derived from TObject and with    comment field //->Class, or for TStreamerSTL::fCtype: pointer elements |
   /// | kObjectP    | Class* derived from TObject and with NO comment field //->Class |
   /// | kAny        | Class  not derived from TObject |
   /// | kAnyp       | Class* not derived from TObject with    comment field //->Class |
   /// | kAnyP       | Class* not derived from TObject with NO comment field //->Class |
   /// | kAnyPnoVT   | Class* not derived from TObject with NO comment field //->Class and Class has NO virtual table |
   /// | kSTLp       | Pointer to STL container |
   /// | kTString    | TString, special case |
   /// | kTObject    | TObject, special case |
   /// | kTNamed     | TNamed , special case |
   /// | kCache      | Cache the value in memory than is not part of the object but is accessible via a SchemaRule |
   enum EReadWrite {
      kBase        =  0,  kOffsetL = 20,  kOffsetP = 40,  kCounter =  6,  kCharStar = 7,
      kChar        =  1,  kShort   =  2,  kInt     =  3,  kLong    =  4,  kFloat    = 5,
      kDouble      =  8,  kDouble32=  9,
      kLegacyChar  = 10, /// Equal to TDataType's kchar
      kUChar       = 11,  kUShort  = 12,  kUInt    = 13,  kULong   = 14,  kBits     = 15,
      kLong64      = 16,  kULong64 = 17,  kBool    = 18,  kFloat16 = 19,
      kObject      = 61,  kAny     = 62,  kObjectp = 63,  kObjectP = 64,  kTString  = 65,
      kTObject     = 66,  kTNamed  = 67,  kAnyp    = 68,  kAnyP    = 69,  kAnyPnoVT = 70,
      kSTLp        = 71,
      kSkip        = 100, kSkipL = 120, kSkipP   = 140,
      kConv        = 200, kConvL = 220, kConvP   = 240,
      kSTL         = 300, kSTLstring = 365,
      kStreamer    = 500, kStreamLoop = 501,
      kCache       = 600,  /// Cache the value in memory than is not part of the object but is accessible via a SchemaRule
      kArtificial  = 1000,
      kCacheNew    = 1001,
      kCacheDelete = 1002,
      kNeedObjectForVirtualBaseClass = 99997,
      kMissing     = 99999
   };

   TStreamerInfo();
   TStreamerInfo(TClass *cl);
   virtual            ~TStreamerInfo();
   void                Build(Bool_t isTransient = kFALSE);
   void                BuildCheck(TFile *file = 0, Bool_t load = kTRUE);
   void                BuildEmulated(TFile *file);
   void                BuildOld();
   virtual Bool_t      BuildFor( const TClass *cl );
   void                CallShowMembers(const void* obj, TMemberInspector &insp, Bool_t isTransient) const;
   void                Clear(Option_t *);
   TObject            *Clone(const char *newname = "") const;
   Bool_t              CompareContent(TClass *cl,TVirtualStreamerInfo *info, Bool_t warn, Bool_t complete, TFile *file);
   void                Compile();
   void                ComputeSize();
   void                ForceWriteInfo(TFile *file, Bool_t force=kFALSE);
   Int_t               GenerateHeaderFile(const char *dirname, const TList *subClasses = 0, const TList *extrainfos = 0);
   TClass             *GetActualClass(const void *obj) const;
   TClass             *GetClass() const {return fClass;}
   UInt_t              GetCheckSum() const {return fCheckSum;}
   UInt_t              GetCheckSum(TClass::ECheckSum code) const;
   Int_t               GetClassVersion() const {return fClassVersion;}
   Int_t               GetDataMemberOffset(TDataMember *dm, TMemberStreamer *&streamer) const;
   TObjArray          *GetElements() const {return fElements;}
   TStreamerElement   *GetElem(Int_t id) const {return fComp[id].fElem;}  // Return the element for the list of optimized elements (max GetNdata())
   TStreamerElement   *GetElement(Int_t id) const {return (TStreamerElement*)fElements->At(id);} // Return the element for the complete list of elements (max GetElements()->GetEntries())
   Int_t               GetElementOffset(Int_t id) const {return fCompFull[id]->fOffset;}
   TStreamerInfoActions::TActionSequence *GetReadMemberWiseActions(Bool_t forCollection) { return forCollection ? fReadMemberWiseVecPtr : fReadMemberWise; }
   TStreamerInfoActions::TActionSequence *GetReadObjectWiseActions() { return fReadObjectWise; }
   TStreamerInfoActions::TActionSequence *GetReadTextActions() { return fReadText; }
   TStreamerInfoActions::TActionSequence *GetWriteMemberWiseActions(Bool_t forCollection) { return forCollection ? fWriteMemberWiseVecPtr : fWriteMemberWise; }
   TStreamerInfoActions::TActionSequence *GetWriteObjectWiseActions() { return fWriteObjectWise; }
   TStreamerInfoActions::TActionSequence *GetWriteTextActions() { return fWriteText; }
   Int_t               GetNdata()   const {return fNdata;}
   Int_t               GetNelement() const { return fElements->GetEntries(); }
   Int_t               GetNumber()  const {return fNumber;}
   Int_t               GetLength(Int_t id) const {return fComp[id].fLength;}
   ULong_t             GetMethod(Int_t id) const {return fComp[id].fMethod;}
   Int_t               GetNewType(Int_t id) const {return fComp[id].fNewType;}
   Int_t               GetOffset(const char *) const;
   Int_t               GetOffset(Int_t id) const {return fComp[id].fOffset;}
   Version_t           GetOldVersion() const {return fOldVersion;}
   Int_t               GetOnFileClassVersion() const {return fOnFileClassVersion;}
   Int_t               GetSize()    const;
   Int_t               GetSizeElements()    const;
   TStreamerElement   *GetStreamerElement(const char*datamember, Int_t& offset) const;
   TStreamerElement   *GetStreamerElementReal(Int_t i, Int_t j) const;
   Int_t               GetType(Int_t id)   const {return fComp[id].fType;}
   template <typename T> T GetTypedValue(char *pointer, Int_t i, Int_t j, Int_t len) const;
   template <typename T> T GetTypedValueClones(TClonesArray *clones, Int_t i, Int_t j, Int_t k, Int_t eoffset) const;
   template <typename T> T GetTypedValueSTL(TVirtualCollectionProxy *cont, Int_t i, Int_t j, Int_t k, Int_t eoffset) const;
   template <typename T> T GetTypedValueSTLP(TVirtualCollectionProxy *cont, Int_t i, Int_t j, Int_t k, Int_t eoffset) const;
   Double_t            GetValue(char *pointer, Int_t i, Int_t j, Int_t len) const { return GetTypedValue<Double_t>(pointer, i, j, len); }
   Double_t            GetValueClones(TClonesArray *clones, Int_t i, Int_t j, Int_t k, Int_t eoffset) const { return GetTypedValueClones<Double_t>(clones, i, j, k, eoffset); }
   Double_t            GetValueSTL(TVirtualCollectionProxy *cont, Int_t i, Int_t j, Int_t k, Int_t eoffset) const { return GetTypedValueSTL<Double_t>(cont, i, j, k, eoffset); }
   Double_t            GetValueSTLP(TVirtualCollectionProxy *cont, Int_t i, Int_t j, Int_t k, Int_t eoffset) const { return GetTypedValueSTLP<Double_t>(cont, i, j, k, eoffset); }
   void                ls(Option_t *option="") const;
   Bool_t              MatchLegacyCheckSum(UInt_t checksum) const;
   TVirtualStreamerInfo *NewInfo(TClass *cl) {return new TStreamerInfo(cl);}
   void               *New(void *obj = 0);
   void               *NewArray(Long_t nElements, void* ary = 0);
   void                Destructor(void* p, Bool_t dtorOnly = kFALSE);
   void                DeleteArray(void* p, Bool_t dtorOnly = kFALSE);
   void                PrintValue(const char *name, char *pointer, Int_t i, Int_t len, Int_t lenmax=1000) const;
   void                PrintValueClones(const char *name, TClonesArray *clones, Int_t i, Int_t eoffset, Int_t lenmax=1000) const;
   void                PrintValueSTL(const char *name, TVirtualCollectionProxy *cont, Int_t i, Int_t eoffset, Int_t lenmax=1000) const;

   template <class T>
   Int_t               ReadBuffer(TBuffer &b, const T &arrptr, TCompInfo *const*const compinfo, Int_t first, Int_t last, Int_t narr=1,Int_t eoffset=0,Int_t mode=0);
   template <class T>
   Int_t               ReadBufferSkip(TBuffer &b, const T &arrptr, const TCompInfo *compinfo,Int_t kase, TStreamerElement *aElement, Int_t narr, Int_t eoffset);
   template <class T>
   Int_t               ReadBufferConv(TBuffer &b, const T &arrptr, const TCompInfo *compinfo,Int_t kase, TStreamerElement *aElement, Int_t narr, Int_t eoffset);
   template <class T>
   Int_t               ReadBufferArtificial(TBuffer &b, const T &arrptr, TStreamerElement *aElement, Int_t narr, Int_t eoffset);

   Int_t               ReadBufferClones(TBuffer &b, TClonesArray *clones, Int_t nc, Int_t first, Int_t eoffset);
   Int_t               ReadBufferSTL(TBuffer &b, TVirtualCollectionProxy *cont, Int_t nc, Int_t eoffset, Bool_t v7 = kTRUE );
   void                SetCheckSum(UInt_t checksum) {fCheckSum = checksum;}
   void                SetClass(TClass *cl);
   void                SetClassVersion(Int_t vers) {fClassVersion=vers;}
   void                SetOnFileClassVersion(Int_t vers) {fOnFileClassVersion=vers;}
   void                TagFile(TFile *fFile);
private:
   // Try to remove those functions from the public interface.
   Int_t               WriteBuffer(TBuffer &b, char *pointer, Int_t first);
   Int_t               WriteBufferClones(TBuffer &b, TClonesArray *clones, Int_t nc, Int_t first, Int_t eoffset);
   Int_t               WriteBufferSTL   (TBuffer &b, TVirtualCollectionProxy *cont,   Int_t nc);
   Int_t               WriteBufferSTLPtrs( TBuffer &b, TVirtualCollectionProxy *cont, Int_t nc, Int_t first, Int_t eoffset);
public:
   virtual void        Update(const TClass *oldClass, TClass *newClass);

   // \brief Generate the TClass and TStreamerInfo for the requested pair.
   // This creates a TVirtualStreamerInfo for the pair and trigger the BuildCheck/Old to
   // provokes the creation of the corresponding TClass.  This relies on the dictionary for
   // std::pair<const int, int> to already exist (or the interpreter information being available)
   // as it is used as a template.
   virtual TVirtualStreamerInfo *GenerateInfoForPair(const std::string &pairclassname, bool silent, size_t hint_pair_offset, size_t hint_pair_size);
   virtual TVirtualStreamerInfo *GenerateInfoForPair(const std::string &firstname, const std::string &secondname, bool silent, size_t hint_pair_offset, size_t hint_pair_size);

   virtual TVirtualCollectionProxy *GenEmulatedProxy(const char* class_name, Bool_t silent);
   virtual TClassStreamer *GenEmulatedClassStreamer(const char* class_name, Bool_t silent);
   virtual TVirtualCollectionProxy *GenExplicitProxy( const ::ROOT::Detail::TCollectionProxyInfo &info, TClass *cl );
   virtual TClassStreamer *GenExplicitClassStreamer( const ::ROOT::Detail::TCollectionProxyInfo &info, TClass *cl );

   static TStreamerElement   *GetCurrentElement();

public:
   // For access by the StreamerInfoActions.
   template <class T>
   Int_t               WriteBufferAux      (TBuffer &b, const T &arr, TCompInfo *const*const compinfo, Int_t first, Int_t last, Int_t narr,Int_t eoffset,Int_t mode);

   //WARNING this class version must be the same as TVirtualStreamerInfo
   ClassDef(TStreamerInfo,9)  //Streamer information for one class version
};


#endif
