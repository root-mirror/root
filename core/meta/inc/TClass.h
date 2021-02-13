// @(#)root/meta:$Id$
// Author: Rene Brun   07/01/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TClass
#define ROOT_TClass

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TClass                                                               //
//                                                                      //
// Dictionary of a class.                                               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TDictionary.h"
#include "TString.h"

#ifdef R__LESS_INCLUDES
class TObjArray;
#else
#include "TObjArray.h"
// Not used in this header file; user code should #include this directly.
// #include "TObjString.h"
// #include "ThreadLocalStorage.h"
// #include <set>
#endif

#include <map>
#include <string>
#include <unordered_set>
#include <vector>
#include <atomic>

class TBaseClass;
class TBrowser;
class TDataMember;
class TCling;
class TMethod;
class TRealData;
class TBuffer;
class TVirtualStreamerInfo;
class TVirtualCollectionProxy;
class TMethodCall;
class TVirtualIsAProxy;
class TVirtualRefProxy;
class THashTable;
class TListOfFunctions;
class TListOfFunctionTemplates;
class TListOfDataMembers;
class TListOfEnums;
class TViewPubFunctions;
class TViewPubDataMembers;
class TFunctionTemplate;
class TProtoClass;

namespace ROOT {
   class TGenericClassInfo;
   class TMapTypeToTClass;
   class TMapDeclIdToTClass;
   namespace Detail {
      class TSchemaRuleSet;
      class TCollectionProxyInfo;
   }
   namespace Internal {
      class TCheckHashRecursiveRemoveConsistency;
   }
}

typedef ROOT::TMapTypeToTClass IdMap_t;
typedef ROOT::TMapDeclIdToTClass DeclIdMap_t;

class TClass : public TDictionary {

friend class TCling;
friend void ROOT::ResetClassVersion(TClass*, const char*, Short_t);
friend class ROOT::TGenericClassInfo;
friend class TProtoClass;
friend class ROOT::Internal::TCheckHashRecursiveRemoveConsistency;
friend class TStreamerInfo;

public:
   // TClass status bits
   enum EStatusBits {
      kReservedLoading = BIT(7), // Internal status bits, set and reset only during initialization

      kClassSaved  = BIT(12),
      kHasLocalHashMember = BIT(14),
      kIgnoreTObjectStreamer = BIT(15),
      kUnloaded    = BIT(16), // The library containing the dictionary for this class was
                              // loaded and has been unloaded from memory.
      kIsTObject = BIT(17),
      kIsForeign   = BIT(18),
      kIsEmulation = BIT(19), // Deprecated
      kStartWithTObject = BIT(20),  // see comments for IsStartingWithTObject()
      kWarned      = BIT(21),
      kHasNameMapNode = BIT(22),
      kHasCustomStreamerMember = BIT(23) // The class has a Streamer method and it is implemented by the user or an older (not StreamerInfo based) automatic streamer.
   };
   enum ENewType { kRealNew = 0, kClassNew, kDummyNew };
   enum ECheckSum {
      kCurrentCheckSum = 0,
      kNoEnum          = 1, // Used since v3.3
      kReflexNoComment = 2, // Up to v5.34.18 (has no range/comment and no typedef at all)
      kNoRange         = 3, // Up to v5.17
      kWithTypeDef     = 4, // Up to v5.34.18 and v5.99/06
      kReflex          = 5, // Up to v5.34.18 (has no typedef at all)
      kNoRangeCheck    = 6, // Up to v5.34.18 and v5.99/06
      kNoBaseCheckSum  = 7, // Up to v5.34.18 and v5.99/06
      kLatestCheckSum  = 8
   };

   // Describe the current state of the TClass itself.
   enum EState {
      kNoInfo,          // The state has not yet been initialized, i.e. the TClass
                        // was just created and/or there is no trace of it in the interpreter.
      kForwardDeclared, // The interpreted knows the entity is a class but that's it.
      kEmulated,        // The information about the class only comes from a TStreamerInfo
      kInterpreted,     // The class is described completely/only in the interpreter database.
      kHasTClassInit,   // The class has a TClass proper bootstrap coming from a run
                        // through rootcling/genreflex/TMetaUtils and the library
                        // containing this dictionary has been loaded in memory.
      kLoaded = kHasTClassInit,
      kNamespaceForMeta // Very transient state necessary to bootstrap namespace entries
                        // in ROOT Meta w/o interpreter information
   };

   // "Typed" pointer that recalls how TClass::New allocated the object.
   // It is returned by TClass:NewObject and should be passed to TClass::DeleteArray or TClass::Destructor
   // to delete the object.
   // It is also used in TVirtualCollectionProxy for the same reasons.
   struct ObjectPtr
   {
      void *fPtr = nullptr;

      TVirtualStreamerInfo *fAllocator = nullptr;

      ObjectPtr(void *ptr = nullptr, TVirtualStreamerInfo *allocator = nullptr) : fPtr(ptr), fAllocator(allocator) {}

      void *GetPtr() const { return fPtr; }

      TVirtualStreamerInfo *GetAllocator() const { return fAllocator; }

      operator bool() const { return fPtr != nullptr; }
   };

private:



   class TDeclNameRegistry {
      // A class which is used to collect decl names starting from normalised
      // names (typedef resolution is excluded here, just string manipulation
      // is performed). At the heart of the implementation, an unordered set.
   public:
      TDeclNameRegistry(Int_t verbLevel=0);
      void AddQualifiedName(const char *name);
      Bool_t HasDeclName(const char *name) const;
      ~TDeclNameRegistry();
   private:
      Int_t fVerbLevel=0;
      std::unordered_set<std::string> fClassNamesSet;
      mutable std::atomic_flag fSpinLock; // MSVC doesn't support = ATOMIC_FLAG_INIT;
   };

   class InsertTClassInRegistryRAII {
      // Trivial RAII used to insert names in the registry
      TClass::EState& fState;
      const char* fName;
      TDeclNameRegistry& fNoInfoOrEmuOrFwdDeclNameRegistry;
   public:
      InsertTClassInRegistryRAII(TClass::EState &state, const char *name, TDeclNameRegistry &emuRegistry);
      ~InsertTClassInRegistryRAII();
   };

   // TClass objects can be created as a result of opening a TFile (in which
   // they are in emulated mode) or as a result of loading the dictionary for
   // the corresponding class.   When a dictionary is loaded any pre-existing
   // emulated TClass is replaced by the one created/coming from the dictionary.
   // To have a reference that always point to the 'current' TClass object for
   // a given class, one should use a TClassRef.
   // TClassRef works by holding on to the fPersistentRef which is updated
   // atomically whenever a TClass is replaced.  During the replacement the
   // value of fPersistentRef is set to zero, leading the TClassRef to call
   // TClass::GetClass which is also locked by the replacement.   At the end
   // of the replacement, fPersistentRef points to the new TClass object.
   std::atomic<TClass**> fPersistentRef;//!Persistent address of pointer to this TClass object and its successors.

   typedef std::atomic<std::map<std::string, TObjArray*>*> ConvSIMap_t;

   mutable TObjArray  *fStreamerInfo;           //Array of TVirtualStreamerInfo
   mutable ConvSIMap_t fConversionStreamerInfo; //Array of the streamer infos derived from another class.
   TList              *fRealData;        //linked list for persistent members including base classes
   std::atomic<TList*> fBase;            //linked list for base classes
   std::atomic<TListOfDataMembers*> fData;            //linked list for data members; non-owning.
   std::atomic<TListOfDataMembers*> fUsingData;//linked list for data members pulled in through using decls.

   std::atomic<TListOfEnums*> fEnums;        //linked list for the enums
   TListOfFunctionTemplates  *fFuncTemplate; //linked list for function templates [Not public until implemented as active list]
   std::atomic<TListOfFunctions*> fMethod;   //linked list for methods

   TViewPubDataMembers*fAllPubData;     //all public data members (including from base classes)
   TViewPubFunctions *fAllPubMethod;    //all public methods (including from base classes)
   mutable TList     *fClassMenuList;   //list of class menu items

   const char        *fDeclFileName;    //name of class declaration file
   const char        *fImplFileName;    //name of class implementation file
   Short_t            fDeclFileLine;    //line of class declaration
   Short_t            fImplFileLine;    //line of class implementation
   UInt_t             fInstanceCount;   //number of instances of this class
   UInt_t             fOnHeap;          //number of instances on heap
   mutable std::atomic<UInt_t>  fCheckSum;        //checksum of data members and base classes
   TVirtualCollectionProxy *fCollectionProxy; //Collection interface
   Version_t          fClassVersion;    //Class version Identifier
   ClassInfo_t       *fClassInfo;       //pointer to CINT class info class
   TString            fContextMenuTitle;//context menu title
   const std::type_info *fTypeInfo;        //pointer to the C++ type information.
   ShowMembersFunc_t  fShowMembers;     //pointer to the class's ShowMembers function
   TClassStreamer    *fStreamer;        //pointer to streamer function
   TString            fSharedLibs;      //shared libraries containing class code

   TVirtualIsAProxy  *fIsA;             //!pointer to the class's IsA proxy.
   IsAGlobalFunc_t    fGlobalIsA;       //pointer to a global IsA function.
   mutable std::atomic<TMethodCall*> fIsAMethod;       //!saved info to call a IsA member function

   ROOT::MergeFunc_t   fMerge;          //pointer to a function implementing Merging objects of this class.
   ROOT::ResetAfterMergeFunc_t fResetAfterMerge; //pointer to a function implementing Merging objects of this class.
   ROOT::NewFunc_t     fNew;            //pointer to a function newing one object.
   ROOT::NewArrFunc_t  fNewArray;       //pointer to a function newing an array of objects.
   ROOT::DelFunc_t     fDelete;         //pointer to a function deleting one object.
   ROOT::DelArrFunc_t  fDeleteArray;    //pointer to a function deleting an array of objects.
   ROOT::DesFunc_t     fDestructor;     //pointer to a function call an object's destructor.
   ROOT::DirAutoAdd_t  fDirAutoAdd;     //pointer which implements the Directory Auto Add feature for this class.']'
   ClassStreamerFunc_t fStreamerFunc;   //Wrapper around this class custom Streamer member function.
   ClassConvStreamerFunc_t fConvStreamerFunc;   //Wrapper around this class custom conversion Streamer member function.
   Int_t               fSizeof;         //Sizeof the class.

   // Bit field
   Int_t fCanSplit : 3;          //!Indicates whether this class can be split or not. Values are -1, 0, 1, 2

   /// Indicates whether this class represents a pair and was not created from a dictionary nor interpreter info but has
   /// compiler compatible offset and size (and all the info is in the StreamerInfo per se)
   Bool_t fIsSyntheticPair : 1;  //!

   mutable std::atomic<Long_t> fProperty; //!Property See TClass::Property() for details
   mutable Long_t     fClassProperty;     //!C++ Property of the class (is abstract, has virtual table, etc.)

           // fHasRootPcmInfo needs to be atomic as long as GetListOfBases needs to modify it.
           std::atomic<Bool_t> fHasRootPcmInfo;      //!Whether info was loaded from a root pcm.
   mutable std::atomic<Bool_t> fCanLoadClassInfo;    //!Indicates whether the ClassInfo is supposed to be available.
   mutable std::atomic<Bool_t> fIsOffsetStreamerSet; //!saved remember if fOffsetStreamer has been set.
   mutable std::atomic<Bool_t> fVersionUsed;         //!Indicates whether GetClassVersion has been called

   enum class ERuntimeProperties : UChar_t {
      kNotInitialized = 0,
      kSet = BIT(0),
      // kInconsistent when kSet & !kConsistent.
      kConsistentHash = BIT(1)
   };
   friend bool operator&(UChar_t l, ERuntimeProperties r) {
      return l & static_cast<UChar_t>(r);
   }
   mutable std::atomic<UChar_t> fRuntimeProperties;    //! Properties that can only be evaluated at run-time

   mutable Long_t     fOffsetStreamer;  //!saved info to call Streamer
   Int_t              fStreamerType;    //!cached of the streaming method to use
   EState             fState;           //!Current 'state' of the class (Emulated,Interpreted,Loaded)
   mutable std::atomic<TVirtualStreamerInfo*>  fCurrentInfo;     //!cached current streamer info.
   mutable std::atomic<TVirtualStreamerInfo*>  fLastReadInfo;    //!cached streamer info used in the last read.
   TVirtualRefProxy  *fRefProxy;        //!Pointer to reference proxy if this class represents a reference
   ROOT::Detail::TSchemaRuleSet *fSchemaRules;  //! Schema evolution rules

   typedef void (*StreamerImpl_t)(const TClass* pThis, void *obj, TBuffer &b, const TClass *onfile_class);
#ifdef R__NO_ATOMIC_FUNCTION_POINTER
   mutable StreamerImpl_t fStreamerImpl; //! Pointer to the function implementing streaming for this class
#else
   mutable std::atomic<StreamerImpl_t> fStreamerImpl; //! Pointer to the function implementing streaming for this class
#endif

   Bool_t             CanSplitBaseAllow();
   TListOfFunctions  *GetMethodList();
   TMethod           *GetClassMethod(Long_t faddr);
   TMethod           *FindClassOrBaseMethodWithId(DeclId_t faddr);
   Int_t              GetBaseClassOffsetRecurse(const TClass *toBase);
   void Init(const char *name, Version_t cversion, const std::type_info *info,
             TVirtualIsAProxy *isa,
             const char *dfil, const char *ifil,
             Int_t dl, Int_t il,
             ClassInfo_t *classInfo,
             Bool_t silent);
   void ForceReload (TClass* oldcl);
   void LoadClassInfo() const;

   static TClass     *LoadClassDefault(const char *requestedname, Bool_t silent);
   static TClass     *LoadClassCustom(const char *requestedname, Bool_t silent);

   void               SetClassVersion(Version_t version);
   void               SetClassSize(Int_t sizof) { fSizeof = sizof; }
   TVirtualStreamerInfo* DetermineCurrentStreamerInfo();

   void SetStreamerImpl();

   void SetRuntimeProperties();

   // Various implementation for TClass::Stramer
   static void StreamerExternal(const TClass* pThis, void *object, TBuffer &b, const TClass *onfile_class);
   static void StreamerTObject(const TClass* pThis, void *object, TBuffer &b, const TClass *onfile_class);
   static void StreamerTObjectInitialized(const TClass* pThis, void *object, TBuffer &b, const TClass *onfile_class);
   static void StreamerTObjectEmulated(const TClass* pThis, void *object, TBuffer &b, const TClass *onfile_class);
   static void StreamerInstrumented(const TClass* pThis, void *object, TBuffer &b, const TClass *onfile_class);
   static void ConvStreamerInstrumented(const TClass* pThis, void *object, TBuffer &b, const TClass *onfile_class);
   static void StreamerStreamerInfo(const TClass* pThis, void *object, TBuffer &b, const TClass *onfile_class);
   static void StreamerDefault(const TClass* pThis, void *object, TBuffer &b, const TClass *onfile_class);

   static IdMap_t    *GetIdMap();       //Map from typeid to TClass pointer
   static DeclIdMap_t *GetDeclIdMap();  //Map from DeclId_t to TClass pointer
   static std::atomic<Int_t>     fgClassCount;  //provides unique id for a each class
                                                //stored in TObject::fUniqueID
   static TDeclNameRegistry fNoInfoOrEmuOrFwdDeclNameRegistry; // Store decl names of the forwardd and no info instances
   static Bool_t HasNoInfoOrEmuOrFwdDeclaredDecl(const char*);

   // Internal status bits, set and reset only during initialization and thus under the protection of the global lock.
   enum { kLoading = kReservedLoading, kUnloading = kReservedLoading };
   // Internal streamer type.
   enum EStreamerType {kDefault=0, kEmulatedStreamer=1, kTObject=2, kInstrumented=4, kForeign=8, kExternal=16};

   // These are the above-referenced hash tables.  (The pointers are null
   // if no entries have been made.)
   static THashTable* fgClassTypedefHash;

   TVirtualStreamerInfo     *GetStreamerInfoImpl(Int_t version, Bool_t silent) const;

   mutable TVirtualMutex *fOVRMutex = nullptr;
   typedef std::multimap<void*, Version_t> RepoCont_t;
   mutable RepoCont_t fObjectVersionRepository;

   void UnregisterAddressInRepository(const char *where, void *location, const TClass *what) const;
   void MoveAddressInRepository(const char *where, void *oldadd, void *newadd, const TClass *what) const;
   void RegisterAddressInRepository(const char *where, void *location, const TClass *what) const;

private:
   TClass(const TClass& tc) = delete;
   TClass& operator=(const TClass&) = delete;

   bool IsClassStructOrUnion() const { return Property() & (kIsClass|kIsStruct|kIsUnion); }
   TList *CreateListOfDataMembers(std::atomic<TListOfDataMembers*> &data, TDictionary::EMemberSelection selection, bool load);

protected:
   TVirtualStreamerInfo *FindStreamerInfo(TObjArray *arr, UInt_t checksum) const;
   void GetMissingDictionariesForBaseClasses(TCollection &result, TCollection &visited, bool recurse);
   void GetMissingDictionariesForMembers(TCollection &result, TCollection &visited, bool recurse);
   void GetMissingDictionariesWithRecursionCheck(TCollection &result, TCollection &visited, bool recurse);
   void GetMissingDictionariesForPairElements(TCollection &result, TCollection &visited, bool recurse);

public:
   TClass();
   TClass(const char *name, Bool_t silent = kFALSE);
   TClass(const char *name, Version_t cversion, Bool_t silent = kFALSE);
   TClass(const char *name, Version_t cversion, EState theState, Bool_t silent = kFALSE);
   TClass(ClassInfo_t *info, Version_t cversion,
          const char *dfil, const char *ifil = 0,
          Int_t dl = 0, Int_t il = 0, Bool_t silent = kFALSE);
   TClass(const char *name, Version_t cversion,
          const char *dfil, const char *ifil = 0,
          Int_t dl = 0, Int_t il = 0, Bool_t silent = kFALSE);
   TClass(const char *name, Version_t cversion,
          const std::type_info &info, TVirtualIsAProxy *isa,
          const char *dfil, const char *ifil,
          Int_t dl, Int_t il, Bool_t silent = kFALSE);
   virtual           ~TClass();

   void               AddInstance(Bool_t heap = kFALSE) { fInstanceCount++; if (heap) fOnHeap++; }
   void               AddImplFile(const char *filename, int line);
   static Bool_t      AddRule(const char *rule);
   static Int_t       ReadRules(const char *filename);
   static Int_t       ReadRules();
   void               AdoptSchemaRules( ROOT::Detail::TSchemaRuleSet *rules );
   virtual void       Browse(TBrowser *b);
   void               BuildRealData(void *pointer=0, Bool_t isTransient = kFALSE);
   void               BuildEmulatedRealData(const char *name, Long_t offset, TClass *cl, Bool_t isTransient = kFALSE);
   void               CalculateStreamerOffset() const;
   Bool_t             CallShowMembers(const void* obj, TMemberInspector &insp, Bool_t isTransient = kFALSE) const;
   Bool_t             CanSplit() const;
   Bool_t             CanIgnoreTObjectStreamer() { return TestBit(kIgnoreTObjectStreamer);}
   Long_t             ClassProperty() const;
   TObject           *Clone(const char *newname="") const;
   void               CopyCollectionProxy(const TVirtualCollectionProxy&);
   void               Draw(Option_t *option="");
   void               Dump() const { TDictionary::Dump(); }
   void               Dump(const void *obj, Bool_t noAddr = kFALSE) const;
   char              *EscapeChars(const char *text) const;
   TVirtualStreamerInfo     *FindStreamerInfo(UInt_t checksum, Bool_t isTransient = kFALSE) const;
   TVirtualStreamerInfo     *GetConversionStreamerInfo( const char* onfile_classname, Int_t version ) const;
   TVirtualStreamerInfo     *FindConversionStreamerInfo( const char* onfile_classname, UInt_t checksum ) const;
   TVirtualStreamerInfo     *GetConversionStreamerInfo( const TClass* onfile_cl, Int_t version ) const;
   TVirtualStreamerInfo     *FindConversionStreamerInfo( const TClass* onfile_cl, UInt_t checksum ) const;
   Bool_t             HasDataMemberInfo() const { return fIsSyntheticPair || fHasRootPcmInfo || HasInterpreterInfo(); }
   Bool_t             HasDefaultConstructor(Bool_t testio = kFALSE) const;
   Bool_t             HasInterpreterInfoInMemory() const { return 0 != fClassInfo; }
   Bool_t             HasInterpreterInfo() const { return fCanLoadClassInfo || fClassInfo; }
   UInt_t             GetCheckSum(ECheckSum code = kCurrentCheckSum) const;
   UInt_t             GetCheckSum(Bool_t &isvalid) const;
   UInt_t             GetCheckSum(ECheckSum code, Bool_t &isvalid) const;
   TVirtualCollectionProxy *GetCollectionProxy() const;
   TVirtualIsAProxy  *GetIsAProxy() const;
   TMethod           *GetClassMethod(const char *name, const char *params, Bool_t objectIsConst = kFALSE);
   TMethod           *GetClassMethodWithPrototype(const char *name, const char *proto,
                                                  Bool_t objectIsConst = kFALSE,
                                                  ROOT::EFunctionMatchMode mode = ROOT::kConversionMatch);
   Version_t          GetClassVersion() const {
      if (!fVersionUsed.load(std::memory_order_relaxed))
         fVersionUsed = kTRUE;
      return fClassVersion;
   }
   Int_t              GetClassSize() const { return Size(); }
   TDataMember       *GetDataMember(const char *datamember) const;
   Long_t             GetDataMemberOffset(const char *membername) const;
   const char        *GetDeclFileName() const;
   Short_t            GetDeclFileLine() const { return fDeclFileLine; }
   ROOT::DelFunc_t    GetDelete() const;
   ROOT::DesFunc_t    GetDestructor() const;
   ROOT::DelArrFunc_t GetDeleteArray() const;
   ClassInfo_t       *GetClassInfo() const {
      if (fCanLoadClassInfo && !TestBit(kLoading))
         LoadClassInfo();
      return fClassInfo;
   }
   const char        *GetContextMenuTitle() const { return fContextMenuTitle; }
   TVirtualStreamerInfo     *GetCurrentStreamerInfo() {
      auto current = fCurrentInfo.load(std::memory_order_relaxed);
      if (current) return current;
      else return DetermineCurrentStreamerInfo();
   }
   TVirtualStreamerInfo     *GetLastReadInfo() const { return fLastReadInfo; }
   void                      SetLastReadInfo(TVirtualStreamerInfo *info) { fLastReadInfo = info; }
   TList             *GetListOfDataMembers(Bool_t load = kTRUE);
   TList             *GetListOfUsingDataMembers(Bool_t load = kTRUE);
   TList             *GetListOfEnums(Bool_t load = kTRUE);
   TList             *GetListOfFunctionTemplates(Bool_t load = kTRUE);
   TList             *GetListOfBases();
   TList             *GetListOfMethods(Bool_t load = kTRUE);
   TCollection       *GetListOfMethodOverloads(const char* name) const;
   TList             *GetListOfRealData() const { return fRealData; }
   const TList       *GetListOfAllPublicMethods(Bool_t load = kTRUE);
   TList             *GetListOfAllPublicDataMembers(Bool_t load = kTRUE);
   const char        *GetImplFileName() const { return fImplFileName; }
   Short_t            GetImplFileLine() const { return fImplFileLine; }
   TClass            *GetActualClass(const void *object) const;
   TClass            *GetBaseClass(const char *classname);
   TClass            *GetBaseClass(const TClass *base);
   Int_t              GetBaseClassOffset(const TClass *toBase, void *address = 0, bool isDerivedObject = true);
   TClass            *GetBaseDataMember(const char *datamember);
   ROOT::ESTLType     GetCollectionType() const;
   ROOT::DirAutoAdd_t GetDirectoryAutoAdd() const;
   TFunctionTemplate *GetFunctionTemplate(const char *name);
   UInt_t             GetInstanceCount() const { return fInstanceCount; }
   UInt_t             GetHeapInstanceCount() const { return fOnHeap; }
   void               GetMenuItems(TList *listitems);
   TList             *GetMenuList() const;
   TMethod           *GetMethod(const char *method, const char *params, Bool_t objectIsConst = kFALSE);
   TMethod *GetMethodWithPrototype(const char *method, const char *proto, Bool_t objectIsConst = kFALSE,
                                   ROOT::EFunctionMatchMode mode = ROOT::kConversionMatch);
   TMethod           *GetMethodAny(const char *method);
   TMethod           *GetMethodAllAny(const char *method);
   Int_t              GetNdata();
   ROOT::MergeFunc_t  GetMerge() const;
   ROOT::ResetAfterMergeFunc_t  GetResetAfterMerge() const;
   ROOT::NewFunc_t    GetNew() const;
   ROOT::NewArrFunc_t GetNewArray() const;
   Int_t              GetNmethods();
   TClass      *const*GetPersistentRef() const { return fPersistentRef; }
   TRealData         *GetRealData(const char *name) const;
   TVirtualRefProxy  *GetReferenceProxy()  const   {  return fRefProxy; }
   const ROOT::Detail::TSchemaRuleSet *GetSchemaRules() const;
   ROOT::Detail::TSchemaRuleSet *GetSchemaRules(Bool_t create = kFALSE);
   const char        *GetSharedLibs();
   ShowMembersFunc_t  GetShowMembersWrapper() const { return fShowMembers; }
   EState             GetState() const { return fState; }
   TClassStreamer    *GetStreamer() const;
   ClassStreamerFunc_t GetStreamerFunc() const;
   ClassConvStreamerFunc_t GetConvStreamerFunc() const;
   const TObjArray          *GetStreamerInfos() const { return fStreamerInfo; }
   TVirtualStreamerInfo     *GetStreamerInfo(Int_t version=0, Bool_t isTransient = kFALSE) const;
   TVirtualStreamerInfo     *GetStreamerInfoAbstractEmulated(Int_t version=0) const;
   TVirtualStreamerInfo     *FindStreamerInfoAbstractEmulated(UInt_t checksum) const;
   const std::type_info     *GetTypeInfo() const { return fTypeInfo; };

   /// @brief Return 'true' if we can guarantee that if this class (or any class in
   /// this class inheritance hierarchy) overload TObject::Hash it also starts
   /// the RecursiveRemove process from its own destructor.
   Bool_t HasConsistentHashMember()
   {
      if (!fRuntimeProperties)
         SetRuntimeProperties();
      return fRuntimeProperties.load() & ERuntimeProperties::kConsistentHash;
   }
   Bool_t             HasDictionary() const;
   static Bool_t      HasDictionarySelection(const char* clname);
   Bool_t HasLocalHashMember() const;
   void               GetMissingDictionaries(THashTable& result, bool recurse = false);
   void               IgnoreTObjectStreamer(Bool_t ignore=kTRUE);
   Bool_t             InheritsFrom(const char *cl) const;
   Bool_t             InheritsFrom(const TClass *cl) const;
   void               InterpretedShowMembers(void* obj, TMemberInspector &insp, Bool_t isTransient);
   Bool_t             IsFolder() const { return kTRUE; }
   Bool_t             IsLoaded() const;
   Bool_t             IsForeign() const;
   Bool_t             IsStartingWithTObject() const;
   Bool_t             IsSyntheticPair() const { return fIsSyntheticPair; }
   Bool_t             IsVersioned() const { return !( GetClassVersion()<=1 && IsForeign() ); }
   Bool_t             IsTObject() const;
   static TClass     *LoadClass(const char *requestedname, Bool_t silent);
   void               ls(Option_t *opt="") const;
   void               MakeCustomMenuList();
   Bool_t             MatchLegacyCheckSum(UInt_t checksum) const;
   void               Move(void *arenaFrom, void *arenaTo) const;
   void              *New(ENewType defConstructor = kClassNew, Bool_t quiet = kFALSE) const;
   void              *New(void *arena, ENewType defConstructor = kClassNew) const;
   void              *NewArray(Long_t nElements, ENewType defConstructor = kClassNew) const;
   void              *NewArray(Long_t nElements, void *arena, ENewType defConstructor = kClassNew) const;
   ObjectPtr          NewObject(ENewType defConstructor = kClassNew, Bool_t quiet = kFALSE) const;
   ObjectPtr          NewObject(void *arena, ENewType defConstructor = kClassNew) const;
   ObjectPtr          NewObjectArray(Long_t nElements, ENewType defConstructor = kClassNew) const;
   ObjectPtr          NewObjectArray(Long_t nElements, void *arena, ENewType defConstructor = kClassNew) const;
   virtual void       PostLoadCheck();
   Long_t             Property() const;
   Int_t              ReadBuffer(TBuffer &b, void *pointer, Int_t version, UInt_t start, UInt_t count);
   Int_t              ReadBuffer(TBuffer &b, void *pointer);
   void               RegisterStreamerInfo(TVirtualStreamerInfo *info);
   void               RemoveStreamerInfo(Int_t slot);
   void               ReplaceWith(TClass *newcl) const;
   void               ResetCaches();
   void               ResetClassInfo(Long_t tagnum);
   void               ResetClassInfo();
   void               ResetInstanceCount() { fInstanceCount = fOnHeap = 0; }
   void               ResetMenuList();
   Int_t              Size() const;
   void               SetCanSplit(Int_t splitmode);
   void               SetCollectionProxy(const ROOT::Detail::TCollectionProxyInfo&);
   void               SetContextMenuTitle(const char *title);
   void               SetCurrentStreamerInfo(TVirtualStreamerInfo *info);
   void               SetGlobalIsA(IsAGlobalFunc_t);
   void               SetDeclFile(const char *name, int line) { fDeclFileName = name; fDeclFileLine = line; }
   void               SetDelete(ROOT::DelFunc_t deleteFunc);
   void               SetDeleteArray(ROOT::DelArrFunc_t deleteArrayFunc);
   void               SetDirectoryAutoAdd(ROOT::DirAutoAdd_t dirAutoAddFunc);
   void               SetDestructor(ROOT::DesFunc_t destructorFunc);
   void               SetImplFileName(const char *implFileName) { fImplFileName = implFileName; }
   void               SetMerge(ROOT::MergeFunc_t mergeFunc);
   void               SetResetAfterMerge(ROOT::ResetAfterMergeFunc_t resetFunc);
   void               SetNew(ROOT::NewFunc_t newFunc);
   void               SetNewArray(ROOT::NewArrFunc_t newArrayFunc);
   TVirtualStreamerInfo     *SetStreamerInfo(Int_t version, const char *info="");
   void               SetUnloaded();
   Int_t              WriteBuffer(TBuffer &b, void *pointer, const char *info="");

   void               AdoptReferenceProxy(TVirtualRefProxy* proxy);
   void               AdoptStreamer(TClassStreamer *strm);
   void               AdoptMemberStreamer(const char *name, TMemberStreamer *strm);
   void               SetMemberStreamer(const char *name, MemberStreamerFunc_t strm);
   void               SetStreamerFunc(ClassStreamerFunc_t strm);
   void               SetConvStreamerFunc(ClassConvStreamerFunc_t strm);

   // Function to retrieve the TClass object and dictionary function
   static void           AddClass(TClass *cl);
   static void           AddClassToDeclIdMap(TDictionary::DeclId_t id, TClass* cl);
   static void           RemoveClass(TClass *cl);
   static void           RemoveClassDeclId(TDictionary::DeclId_t id);
   static TClass        *GetClass(const char *name, Bool_t load = kTRUE, Bool_t silent = kFALSE);
   static TClass        *GetClass(const char *name, Bool_t load, Bool_t silent, size_t hint_pair_offset, size_t hint_pair_size);
   static TClass        *GetClass(const std::type_info &typeinfo, Bool_t load = kTRUE, Bool_t silent = kFALSE, size_t hint_pair_offset = 0, size_t hint_pair_size = 0);
   static TClass        *GetClass(ClassInfo_t *info, Bool_t load = kTRUE, Bool_t silent = kFALSE);
   template<typename T>
   static TClass        *GetClass(Bool_t load = kTRUE, Bool_t silent = kFALSE);
   static Bool_t         GetClass(DeclId_t id, std::vector<TClass*> &classes);
   static DictFuncPtr_t  GetDict (const char *cname);
   static DictFuncPtr_t  GetDict (const std::type_info &info);

   static Int_t       AutoBrowse(TObject *obj, TBrowser *browser);
   static ENewType    IsCallingNew();
   static TClass     *Load(TBuffer &b);
   void               Store(TBuffer &b) const;

   // Pseudo-method apply to the 'obj'. In particular those are used to
   // implement TObject like methods for non-TObject classes.
   Int_t              Browse(void *obj, TBrowser *b) const;
   void               DeleteArray(void *ary, Bool_t dtorOnly = kFALSE);
   void               DeleteArray(ObjectPtr ary, Bool_t dtorOnly = kFALSE);
   void               Destructor(void *obj, Bool_t dtorOnly = kFALSE);
   void               Destructor(ObjectPtr obj, Bool_t dtorOnly = kFALSE);
   void              *DynamicCast(const TClass *base, void *obj, Bool_t up = kTRUE);
   const void        *DynamicCast(const TClass *base, const void *obj, Bool_t up = kTRUE);
   Bool_t             IsFolder(void *obj) const;

   inline void        Streamer(void *obj, TBuffer &b, const TClass *onfile_class = 0) const
   {
      // Inline for performance, skipping one function call.
#ifdef R__NO_ATOMIC_FUNCTION_POINTER
      fStreamerImpl(this,obj,b,onfile_class);
#else
      auto t = fStreamerImpl.load();
      t(this,obj,b,onfile_class);
#endif
   }

   ClassDef(TClass,0)  //Dictionary containing class information
};

namespace ROOT {
namespace Internal {
template <typename T>
TClass *GetClassHelper(Bool_t, Bool_t, std::true_type)
{
   return T::Class();
}

template <typename T>
struct TClassGetClassHelper {
   static TClass *GetClass(Bool_t load, Bool_t silent) {
      return TClass::GetClass(typeid(T), load, silent);
   }
};

template <typename F, typename S>
struct TClassGetClassHelper<std::pair<F, S> > {
   static TClass *GetClass(Bool_t load, Bool_t silent) {
      std::pair<F, S> *p = nullptr;
      size_t hint_offset = ((char*)&(p->second)) - (char*)p;
      return TClass::GetClass(typeid(std::pair<F, S>), load, silent, hint_offset, sizeof(std::pair<F,S>));
   }
};

template <typename T>
TClass *GetClassHelper(Bool_t load, Bool_t silent, std::false_type)
{
   return TClassGetClassHelper<T>::GetClass(load, silent);
}

} // namespace Internal
} // namespace ROOT

template <typename T>
TClass *TClass::GetClass(Bool_t load, Bool_t silent)
{
   typename std::is_base_of<TObject, T>::type tag;
   return ROOT::Internal::GetClassHelper<T>(load, silent, tag);
}

namespace ROOT {

template <typename T> TClass *GetClass(T * /* dummy */)       { return TClass::GetClass<T>(); }
template <typename T> TClass *GetClass(const T * /* dummy */) { return TClass::GetClass<T>(); }

#ifndef R__NO_CLASS_TEMPLATE_SPECIALIZATION
   // This can only be used when the template overload resolution can distinguish between T* and T**
   template <typename T> TClass* GetClass(      T**       /* dummy */) { return TClass::GetClass<T>(); }
   template <typename T> TClass* GetClass(const T**       /* dummy */) { return TClass::GetClass<T>(); }
   template <typename T> TClass* GetClass(      T* const* /* dummy */) { return TClass::GetClass<T>(); }
   template <typename T> TClass* GetClass(const T* const* /* dummy */) { return TClass::GetClass<T>(); }
#endif

   extern TClass *CreateClass(const char *cname, Version_t id,
                              const char *dfil, const char *ifil,
                              Int_t dl, Int_t il);
}

#endif // ROOT_TClass
