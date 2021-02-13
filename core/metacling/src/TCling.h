// @(#)root/meta:$Id$
// Author: Axel Naumann, 2011-10-19

/*************************************************************************
 * Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


#ifndef ROOT_TCling
#define ROOT_TCling

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TCling                                                               //
//                                                                      //
// This class defines an interface to the cling C++ interpreter.        //
//                                                                      //
// Cling is a full ANSI compliant C++ interpreter based on              //
// clang/LLVM technology.                                               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TInterpreter.h"

#include <map>
#include <memory>
#include <set>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>
#include <utility>

#ifndef WIN32
#define TWin32SendClass char
#endif

namespace llvm {
   class GlobalValue;
   class StringRef;
}

namespace clang {
   class CXXRecordDecl;
   class Decl;
   class DeclContext;
   class EnumDecl;
   class FunctionDecl;
   class IdentifierIterator;
   class NamedDecl;
   class NamespaceDecl;
   class TagDecl;
   class Type;
   class QualType;
}
namespace cling {
   class Interpreter;
   class MetaProcessor;
   class Transaction;
   class Value;
}

class TClingCallbacks;
class TEnv;
class TFile;
class THashTable;
class TInterpreterValue;
class TMethod;
class TObjArray;
class TListOfDataMembers;
class TListOfFunctions;
class TListOfFunctionTemplates;
class TListOfEnums;

namespace ROOT {
   namespace TMetaUtils {
      class TNormalizedCtxt;
      class TClingLookupHelper;
   }
}

extern "C" {
   void TCling__UpdateListsOnCommitted(const cling::Transaction&,
                                       cling::Interpreter*);
   void TCling__UpdateListsOnUnloaded(const cling::Transaction&);
   void TCling__InvalidateGlobal(const clang::Decl*);
   void TCling__TransactionRollback(const cling::Transaction&);
   TObject* TCling__GetObjectAddress(const char *Name, void *&LookupCtx);
   const clang::Decl* TCling__GetObjectDecl(TObject *obj);
   void TCling__LibraryLoaded(const void* dyLibHandle,
                              const char* canonicalName);
   void TCling__LibraryUnloaded(const void* dyLibHandle,
                                const char* canonicalName);
   void TCling__RegisterRdictForLoadPCM(const std::string &pcmFileNameFullPath, llvm::StringRef *pcmContent);
}

class TCling final : public TInterpreter {
private: // Static Data Members

   static void* fgSetOfSpecials; // set of TObjects used in CINT variables

private: // Data Members

   Int_t           fMore;             // The brace indent level for the cint command line processor.
   Int_t           fExitCode;         // Value passed to exit() in interpreter.
   char            fPrompt[64];       // Command line prompt string.
   //cling::DictPosition fDictPos;          // dictionary context after initialization is complete.
   //cling::DictPosition fDictPosGlobals;   // dictionary context after ResetGlobals().
   TString         fSharedLibs;       // Shared libraries loaded by G__loadfile().
   Int_t           fGlobalsListSerial;// Last time we refreshed the ROOT list of globals.
   TString         fIncludePath;      // Interpreter include path.
   TString         fRootmapLoadPath;  // Dynamic load path for rootmap files.
   TEnv*           fMapfile;          // Association of classes to libraries.
   std::vector<std::string> fAutoLoadLibStorage; // A storage to return a const char* from GetClassSharedLibsForModule.
   std::map<size_t,std::vector<const char*>> fClassesHeadersMap; // Map of classes hashes and headers associated
   std::map<const cling::Transaction*,size_t> fTransactionHeadersMap; // Map which transaction contains which autoparse.
   std::set<size_t> fLookedUpClasses; // Set of classes for which headers were looked up already
   std::set<size_t> fPayloads; // Set of payloads
   std::set<const char*> fParsedPayloadsAddresses; // Set of payloads which were parsed
   std::hash<std::string> fStringHashFunction; // A simple hashing function
   std::unordered_set<const clang::NamespaceDecl*> fNSFromRootmaps;   // Collection of namespaces fwd declared in the rootmaps
   TObjArray*      fRootmapFiles;     // Loaded rootmap files.
   Bool_t          fLockProcessLine;  // True if ProcessLine should lock gInterpreterMutex.
   Bool_t          fCxxModulesEnabled;// True if C++ modules was enabled

   std::unique_ptr<cling::Interpreter>   fInterpreter;   // The interpreter.
   std::unique_ptr<cling::MetaProcessor> fMetaProcessor; // The metaprocessor.

   std::vector<cling::Value> *fTemporaries;    // Stack of temporaries
   ROOT::TMetaUtils::TNormalizedCtxt  *fNormalizedCtxt; // Which typedef to avoid stripping.
   ROOT::TMetaUtils::TClingLookupHelper *fLookupHelper; // lookup helper used by TClassEdit

   void*           fPrevLoadedDynLibInfo; // Internal info to mark the last loaded libray.
   std::vector<void*> fRegisterModuleDyLibs; // Stack of libraries currently running RegisterModule
   TClingCallbacks* fClingCallbacks; // cling::Interpreter owns it.
   struct CharPtrCmp_t {
      bool operator()(const char* a, const char *b) const {
         return strcmp(a, b) < 0;
      }
   };
   std::set<TClass*> fModTClasses;
   std::vector<std::pair<TClass*,DictFuncPtr_t> > fClassesToUpdate;
   void* fAutoLoadCallBack;
   ULong64_t fTransactionCount; // Cling counter for commited or unloaded transactions which changed the AST.
   std::vector<const char*> fCurExecutingMacros;

   typedef void* SpecialObjectLookupCtx_t;
   typedef std::unordered_map<std::string, TObject*> SpecialObjectMap_t;
   std::map<SpecialObjectLookupCtx_t, SpecialObjectMap_t> fSpecialObjectMaps;

   struct MutexStateAndRecurseCount {
      /// State of gCoreMutex when the first interpreter-related function was invoked.
      std::unique_ptr<ROOT::TVirtualRWMutex::State> fState;

      /// Interpreter-related functions will push the "entry" lock state to *this.
      /// Recursive calls will do that, too - but we must only forget about the lock
      /// state once this recursion count went to 0.
      Int_t fRecurseCount = 0;

      operator bool() const { return (bool)fState; }
   };

   struct MutexStateAndRecurseCountDelta {
      using StateDelta = ROOT::TVirtualRWMutex::StateDelta;

      MutexStateAndRecurseCount   fInitialState;
      std::unique_ptr<StateDelta> fDelta;
   };

   MutexStateAndRecurseCount fInitialMutex;

   DeclId_t GetDeclId(const llvm::GlobalValue *gv) const;

   static Int_t DeepAutoLoadImpl(const char *cls);
   static Int_t ShallowAutoLoadImpl(const char *cls);

   Bool_t fHeaderParsingOnDemand;
   Bool_t fIsAutoParsingSuspended;

   UInt_t AutoParseImplRecurse(const char *cls, bool topLevel);
   constexpr static const char* kNullArgv[] = {nullptr};

   bool fIsShuttingDown = false;

protected:
   Bool_t SetSuspendAutoParsing(Bool_t value);

public: // Public Interface

   virtual ~TCling();
   TCling(const char* name, const char* title, const char* const argv[]);
   TCling(const char* name, const char* title): TCling(name, title, kNullArgv) {}

   void    AddIncludePath(const char* path);
   void   *GetAutoLoadCallBack() const { return fAutoLoadCallBack; }
   void   *SetAutoLoadCallBack(void* cb) { void* prev = fAutoLoadCallBack; fAutoLoadCallBack = cb; return prev; }
   Int_t   AutoLoad(const char *classname, Bool_t knowDictNotLoaded = kFALSE);
   Int_t   AutoLoad(const std::type_info& typeinfo, Bool_t knowDictNotLoaded = kFALSE);
   Int_t   AutoParse(const char* cls);
   void*   LazyFunctionCreatorAutoload(const std::string& mangled_name);
   bool   LibraryLoadingFailed(const std::string&, const std::string&, bool, bool);
   Bool_t  IsAutoLoadNamespaceCandidate(const clang::NamespaceDecl* nsDecl);
   void    ClearFileBusy();
   void    ClearStack(); // Delete existing temporary values
   Bool_t  Declare(const char* code);
   void    EndOfLineAction();
   TClass *GetClass(const std::type_info& typeinfo, Bool_t load) const;
   Int_t   GetExitCode() const { return fExitCode; }
   TEnv*   GetMapfile() const { return fMapfile; }
   Int_t   GetMore() const { return fMore; }
   TClass *GenerateTClass(const char *classname, Bool_t emulation, Bool_t silent = kFALSE);
   TClass *GenerateTClass(ClassInfo_t *classinfo, Bool_t silent = kFALSE);
   Int_t   GenerateDictionary(const char* classes, const char* includes = "", const char* options = 0);
   char*   GetPrompt() { return fPrompt; }
   const char* GetSharedLibs();
   const char* GetClassSharedLibs(const char* cls);
   const char* GetSharedLibDeps(const char* lib, bool tryDyld = false);
   const char* GetIncludePath();
   virtual const char* GetSTLIncludePath() const;
   TObjArray*  GetRootMapFiles() const { return fRootmapFiles; }
   unsigned long long GetInterpreterStateMarker() const { return fTransactionCount;}
   virtual void Initialize();
   virtual void ShutDown();
   void    InspectMembers(TMemberInspector&, const void* obj, const TClass* cl, Bool_t isTransient);
   Bool_t  IsLoaded(const char* filename) const;
   Bool_t  IsLibraryLoaded(const char* libname) const;
   Bool_t  HasPCMForLibrary(const char *libname) const;
   Int_t   Load(const char* filenam, Bool_t system = kFALSE);
   void    LoadMacro(const char* filename, EErrorCode* error = 0);
   Int_t   LoadLibraryMap(const char* rootmapfile = 0);
   Int_t   RescanLibraryMap();
   Int_t   ReloadAllSharedLibraryMaps();
   Int_t   UnloadAllSharedLibraryMaps();
   Int_t   UnloadLibraryMap(const char* library);
   Long_t  ProcessLine(const char* line, EErrorCode* error = 0);
   Long_t  ProcessLineAsynch(const char* line, EErrorCode* error = 0);
   Long_t  ProcessLineSynch(const char* line, EErrorCode* error = 0);
   void    PrintIntro();
   bool    RegisterPrebuiltModulePath(const std::string& FullPath,
                                      const std::string& ModuleMapName = "module.modulemap") const;
   void    RegisterModule(const char* modulename,
                          const char** headers,
                          const char** includePaths,
                          const char* payloadCode,
                          const char* fwdDeclsCode,
                          void (*triggerFunc)(),
                          const FwdDeclArgsToKeepCollection_t& fwdDeclsArgToSkip,
                          const char** classesHeaders,
                          Bool_t lateRegistration = false,
                          Bool_t hasCxxModule = false);
   virtual void AddAvailableIndentifiers(TSeqCollection& Idents);
   void    RegisterTClassUpdate(TClass *oldcl,DictFuncPtr_t dict);
   void    UnRegisterTClassUpdate(const TClass *oldcl);

   Int_t   SetClassSharedLibs(const char *cls, const char *libs);
   void    SetGetline(const char * (*getlineFunc)(const char* prompt),
                      void (*histaddFunc)(const char* line));
   void    Reset();
   void    ResetAll();
   void    ResetGlobals();
   void    ResetGlobalVar(void* obj);
   void    RewindDictionary();
   Int_t   DeleteGlobal(void* obj);
   Int_t   DeleteVariable(const char *name);
   void    SaveContext();
   void    SaveGlobalsContext();
   void    UpdateListOfGlobals();
   void    UpdateListOfGlobalFunctions();
   void    UpdateListOfTypes();
   void    SetClassInfo(TClass* cl, Bool_t reload = kFALSE);

   ECheckClassInfo CheckClassInfo(const char *name, Bool_t autoload, Bool_t isClassOrNamespaceOnly = kFALSE);

   Bool_t  CheckClassTemplate(const char *name);
   Long_t  Calc(const char* line, EErrorCode* error = 0);
   void    CreateListOfBaseClasses(TClass* cl) const;
   void    CreateListOfDataMembers(TClass* cl) const;
   void    CreateListOfMethods(TClass* cl) const;
   void    CreateListOfMethodArgs(TFunction* m) const;
   void    UpdateListOfMethods(TClass* cl) const;
   void    UpdateListOfDataMembers(TClass* cl) const;

   virtual DeclId_t GetDataMember(ClassInfo_t *cl, const char *name) const;
   virtual DeclId_t GetDataMemberAtAddr(const void *addr) const;
   virtual DeclId_t GetDataMemberWithValue(const void *ptrvalue) const;
   virtual DeclId_t GetEnum(TClass *cl, const char *name) const;
   virtual TEnum*   CreateEnum(void *VD, TClass *cl) const;
   virtual void     UpdateEnumConstants(TEnum* enumObj, TClass* cl) const;
   virtual void     LoadEnums(TListOfEnums& cl) const;
   virtual std::string ToString(const char* type, void *obj);
   TString GetMangledName(TClass* cl, const char* method, const char* params, Bool_t objectIsConst = kFALSE);
   TString GetMangledNameWithPrototype(TClass* cl, const char* method, const char* proto, Bool_t objectIsConst = kFALSE, ROOT::EFunctionMatchMode mode = ROOT::kConversionMatch);
   void*   GetInterfaceMethod(TClass* cl, const char* method, const char* params, Bool_t objectIsConst = kFALSE);
   void*   GetInterfaceMethodWithPrototype(TClass* cl, const char* method, const char* proto, Bool_t objectIsConst = kFALSE, ROOT::EFunctionMatchMode mode = ROOT::kConversionMatch);
   DeclId_t GetFunction(ClassInfo_t *cl, const char *funcname);
   DeclId_t GetFunctionWithPrototype(ClassInfo_t *cl, const char* method, const char* proto, Bool_t objectIsConst = kFALSE, ROOT::EFunctionMatchMode mode = ROOT::kConversionMatch);
   DeclId_t GetFunctionWithValues(ClassInfo_t *cl, const char* method, const char* params, Bool_t objectIsConst = kFALSE);
   DeclId_t GetFunctionTemplate(ClassInfo_t *cl, const char *funcname);
   void     GetFunctionOverloads(ClassInfo_t *cl, const char *funcname, std::vector<DeclId_t>& res) const;
   virtual void     LoadFunctionTemplates(TClass* cl) const;

   virtual std::vector<std::string> GetUsingNamespaces(ClassInfo_t *cl) const;

   void    GetInterpreterTypeName(const char* name, std::string &output, Bool_t full = kFALSE);
   void    Execute(const char* function, const char* params, int* error = 0);
   void    Execute(TObject* obj, TClass* cl, const char* method, const char* params, int* error = 0);
   void    Execute(TObject* obj, TClass* cl, const char* method, const char* params, Bool_t objectIsConst, int* error = 0);
   void    Execute(TObject* obj, TClass* cl, TMethod* method, TObjArray* params, int* error = 0);
   void    ExecuteWithArgsAndReturn(TMethod* method, void* address, const void* args[] = 0, int nargs = 0, void* ret= 0) const;
   Long_t  ExecuteMacro(const char* filename, EErrorCode* error = 0);
   void    RecursiveRemove(TObject* obj);
   Bool_t  IsErrorMessagesEnabled() const;
   Bool_t  SetErrorMessages(Bool_t enable = kTRUE);
   Bool_t  IsProcessLineLocked() const {
      return fLockProcessLine;
   }
   void    SetProcessLineLock(Bool_t lock = kTRUE) {
      fLockProcessLine = lock;
   }
   const char* TypeName(const char* typeDesc);

   void     SnapshotMutexState(ROOT::TVirtualRWMutex* mtx);
   void     ForgetMutexState();

   void     ApplyToInterpreterMutex(void* delta);
   void    *RewindInterpreterMutex();

   static void  UpdateClassInfo(char* name, Long_t tagnum);
   static void  UpdateClassInfoWork(const char* name);
          void  RefreshClassInfo(TClass *cl, const clang::NamedDecl *def, bool alias);
          void  UpdateClassInfoWithDecl(const clang::NamedDecl* ND);
   static void  UpdateAllCanvases();

   // Misc
   virtual int    DisplayClass(FILE* fout, const char* name, int base, int start) const;
   virtual int    DisplayIncludePath(FILE* fout) const;
   virtual void*  FindSym(const char* entry) const;
   virtual void   GenericError(const char* error) const;
   virtual Long_t GetExecByteCode() const;
   virtual const char* GetTopLevelMacroName() const;
   virtual const char* GetCurrentMacroName() const;
   virtual int    GetSecurityError() const;
   virtual int    LoadFile(const char* path) const;
   virtual Bool_t LoadText(const char* text) const;
   virtual const char* MapCppName(const char*) const;
   virtual void   SetAlloclockfunc(void (*)()) const;
   virtual void   SetAllocunlockfunc(void (*)()) const;
   virtual int    SetClassAutoLoading(int) const;
   virtual int    SetClassAutoparsing(int) ;
           Bool_t IsAutoParsingSuspended() const { return fIsAutoParsingSuspended; }
   virtual void   SetErrmsgcallback(void* p) const;
   virtual void   SetTempLevel(int val) const;
   virtual int    UnloadFile(const char* path) const;

   void               CodeComplete(const std::string&, size_t&,
                                   std::vector<std::string>&);
   virtual int Evaluate(const char*, TInterpreterValue&);
   virtual std::unique_ptr<TInterpreterValue> MakeInterpreterValue() const;
   void               RegisterTemporary(const TInterpreterValue& value);
   void               RegisterTemporary(const cling::Value& value);
   const ROOT::TMetaUtils::TNormalizedCtxt& GetNormalizedContext() const {return *fNormalizedCtxt;};
   TObject* GetObjectAddress(const char *Name, void *&LookupCtx);


   // core/meta helper functions.
   virtual EReturnType MethodCallReturnType(TFunction *func) const;
   virtual void GetFunctionName(const clang::Decl *decl, std::string &name) const;
   virtual bool DiagnoseIfInterpreterException(const std::exception &e) const;

   // CallFunc interface
   virtual DeclId_t GetDeclId(CallFunc_t *info) const;
   virtual void   CallFunc_Delete(CallFunc_t* func) const;
   virtual void   CallFunc_Exec(CallFunc_t* func, void* address) const;
   virtual void   CallFunc_Exec(CallFunc_t* func, void* address, TInterpreterValue& val) const;
   virtual void   CallFunc_ExecWithReturn(CallFunc_t* func, void* address, void* ret) const;
   virtual void   CallFunc_ExecWithArgsAndReturn(CallFunc_t* func, void* address, const void* args[] = 0, int nargs = 0, void* ret = 0) const;
   virtual Long_t    CallFunc_ExecInt(CallFunc_t* func, void* address) const;
   virtual Long64_t  CallFunc_ExecInt64(CallFunc_t* func, void* address) const;
   virtual Double_t  CallFunc_ExecDouble(CallFunc_t* func, void* address) const;
   virtual CallFunc_t*   CallFunc_Factory() const;
   virtual CallFunc_t*   CallFunc_FactoryCopy(CallFunc_t* func) const;
   virtual MethodInfo_t* CallFunc_FactoryMethod(CallFunc_t* func) const;
   virtual void   CallFunc_IgnoreExtraArgs(CallFunc_t* func, bool ignore) const;
   virtual void   CallFunc_Init(CallFunc_t* func) const;
   virtual bool   CallFunc_IsValid(CallFunc_t* func) const;
   virtual CallFuncIFacePtr_t CallFunc_IFacePtr(CallFunc_t * func) const;
   virtual void   CallFunc_ResetArg(CallFunc_t* func) const;
   virtual void   CallFunc_SetArg(CallFunc_t* func, Long_t param) const;
   virtual void   CallFunc_SetArg(CallFunc_t* func, ULong_t param) const;
   virtual void   CallFunc_SetArg(CallFunc_t* func, Float_t param) const;
   virtual void   CallFunc_SetArg(CallFunc_t* func, Double_t param) const;
   virtual void   CallFunc_SetArg(CallFunc_t* func, Long64_t param) const;
   virtual void   CallFunc_SetArg(CallFunc_t* func, ULong64_t param) const;
   virtual void   CallFunc_SetArgArray(CallFunc_t* func, Long_t* paramArr, Int_t nparam) const;
   virtual void   CallFunc_SetArgs(CallFunc_t* func, const char* param) const;
   virtual void   CallFunc_SetFunc(CallFunc_t* func, ClassInfo_t* info, const char* method, const char* params, Long_t* Offset) const;
   virtual void   CallFunc_SetFunc(CallFunc_t* func, ClassInfo_t* info, const char* method, const char* params, bool objectIsConst, Long_t* Offset) const;
   virtual void   CallFunc_SetFunc(CallFunc_t* func, MethodInfo_t* info) const;
   virtual void   CallFunc_SetFuncProto(CallFunc_t* func, ClassInfo_t* info, const char* method, const char* proto, Long_t* Offset, ROOT::EFunctionMatchMode mode = ROOT::kConversionMatch) const;
   virtual void   CallFunc_SetFuncProto(CallFunc_t* func, ClassInfo_t* info, const char* method, const char* proto, bool objectIsConst, Long_t* Offset, ROOT::EFunctionMatchMode mode = ROOT::kConversionMatch) const;
   virtual void   CallFunc_SetFuncProto(CallFunc_t* func, ClassInfo_t* info, const char* method, const std::vector<TypeInfo_t*> &proto, Long_t* Offset, ROOT::EFunctionMatchMode mode = ROOT::kConversionMatch) const;
   virtual void   CallFunc_SetFuncProto(CallFunc_t* func, ClassInfo_t* info, const char* method, const std::vector<TypeInfo_t*> &proto, bool objectIsConst, Long_t* Offset, ROOT::EFunctionMatchMode mode = ROOT::kConversionMatch) const;

   virtual std::string CallFunc_GetWrapperCode(CallFunc_t *func) const;

   // ClassInfo interface
   virtual DeclId_t GetDeclId(ClassInfo_t *info) const;
   virtual Bool_t ClassInfo_Contains(ClassInfo_t *info, DeclId_t declid) const;
   virtual Long_t ClassInfo_ClassProperty(ClassInfo_t* info) const;
   virtual void   ClassInfo_Delete(ClassInfo_t* info) const;
   virtual void   ClassInfo_Delete(ClassInfo_t* info, void* arena) const;
   virtual void   ClassInfo_DeleteArray(ClassInfo_t* info, void* arena, bool dtorOnly) const;
   virtual void   ClassInfo_Destruct(ClassInfo_t* info, void* arena) const;
   virtual ClassInfo_t*  ClassInfo_Factory(Bool_t all = kTRUE) const;
   virtual ClassInfo_t*  ClassInfo_Factory(ClassInfo_t* cl) const;
   virtual ClassInfo_t*  ClassInfo_Factory(const char* name) const;
   virtual ClassInfo_t*  ClassInfo_Factory(DeclId_t declid) const;
   virtual Long_t   ClassInfo_GetBaseOffset(ClassInfo_t* fromDerived, ClassInfo_t* toBase, void * address, bool isDerivedObject) const;
   virtual int    ClassInfo_GetMethodNArg(ClassInfo_t* info, const char* method, const char* proto, Bool_t objectIsConst = false, ROOT::EFunctionMatchMode mode = ROOT::kConversionMatch) const;
   virtual bool   ClassInfo_HasDefaultConstructor(ClassInfo_t* info, Bool_t testio = kFALSE) const;
   virtual bool   ClassInfo_HasMethod(ClassInfo_t* info, const char* name) const;
   virtual void   ClassInfo_Init(ClassInfo_t* info, const char* funcname) const;
   virtual void   ClassInfo_Init(ClassInfo_t* info, int tagnum) const;
   virtual bool   ClassInfo_IsBase(ClassInfo_t* info, const char* name) const;
   virtual bool   ClassInfo_IsEnum(const char* name) const;
   virtual bool   ClassInfo_IsScopedEnum(ClassInfo_t* info) const;
   virtual EDataType ClassInfo_GetUnderlyingType(ClassInfo_t* info) const;
   virtual bool   ClassInfo_IsLoaded(ClassInfo_t* info) const;
   virtual bool   ClassInfo_IsValid(ClassInfo_t* info) const;
   virtual bool   ClassInfo_IsValidMethod(ClassInfo_t* info, const char* method, const char* proto, Long_t* offset, ROOT::EFunctionMatchMode /* mode */ = ROOT::kConversionMatch) const;
   virtual bool   ClassInfo_IsValidMethod(ClassInfo_t* info, const char* method, const char* proto, Bool_t objectIsConst, Long_t* offset, ROOT::EFunctionMatchMode /* mode */ = ROOT::kConversionMatch) const;
   virtual int    ClassInfo_Next(ClassInfo_t* info) const;
   virtual void*  ClassInfo_New(ClassInfo_t* info) const;
   virtual void*  ClassInfo_New(ClassInfo_t* info, int n) const;
   virtual void*  ClassInfo_New(ClassInfo_t* info, int n, void* arena) const;
   virtual void*  ClassInfo_New(ClassInfo_t* info, void* arena) const;
   virtual Long_t ClassInfo_Property(ClassInfo_t* info) const;
   virtual int    ClassInfo_Size(ClassInfo_t* info) const;
   virtual Long_t ClassInfo_Tagnum(ClassInfo_t* info) const;
   virtual const char* ClassInfo_FileName(ClassInfo_t* info) const;
   virtual const char* ClassInfo_FullName(ClassInfo_t* info) const;
   virtual const char* ClassInfo_Name(ClassInfo_t* info) const;
   virtual const char* ClassInfo_Title(ClassInfo_t* info) const;
   virtual const char* ClassInfo_TmpltName(ClassInfo_t* info) const;

   // BaseClassInfo interface
   virtual void   BaseClassInfo_Delete(BaseClassInfo_t* bcinfo) const;
   virtual BaseClassInfo_t*  BaseClassInfo_Factory(ClassInfo_t* info) const;
   virtual BaseClassInfo_t*  BaseClassInfo_Factory(ClassInfo_t* derived,
                                                   ClassInfo_t* base) const;
   virtual int    BaseClassInfo_Next(BaseClassInfo_t* bcinfo) const;
   virtual int    BaseClassInfo_Next(BaseClassInfo_t* bcinfo, int onlyDirect) const;
   virtual Long_t BaseClassInfo_Offset(BaseClassInfo_t* toBaseClassInfo, void * address, bool isDerivedObject) const;
   virtual Long_t BaseClassInfo_Property(BaseClassInfo_t* bcinfo) const;
   virtual Long_t BaseClassInfo_Tagnum(BaseClassInfo_t* bcinfo) const;
   virtual ClassInfo_t*BaseClassInfo_ClassInfo(BaseClassInfo_t * /* bcinfo */) const;
   virtual const char* BaseClassInfo_FullName(BaseClassInfo_t* bcinfo) const;
   virtual const char* BaseClassInfo_Name(BaseClassInfo_t* bcinfo) const;
   virtual const char* BaseClassInfo_TmpltName(BaseClassInfo_t* bcinfo) const;

   // DataMemberInfo interface
   virtual DeclId_t GetDeclId(DataMemberInfo_t *info) const;
   virtual int    DataMemberInfo_ArrayDim(DataMemberInfo_t* dminfo) const;
   virtual void   DataMemberInfo_Delete(DataMemberInfo_t* dminfo) const;
   virtual DataMemberInfo_t*  DataMemberInfo_Factory(ClassInfo_t* clinfo, TDictionary::EMemberSelection selection) const;
   virtual DataMemberInfo_t  *DataMemberInfo_Factory(DeclId_t declid, ClassInfo_t* clinfo) const;
   virtual DataMemberInfo_t*  DataMemberInfo_FactoryCopy(DataMemberInfo_t* dminfo) const;
   virtual bool   DataMemberInfo_IsValid(DataMemberInfo_t* dminfo) const;
   virtual int    DataMemberInfo_MaxIndex(DataMemberInfo_t* dminfo, Int_t dim) const;
   virtual int    DataMemberInfo_Next(DataMemberInfo_t* dminfo) const;
   virtual Long_t DataMemberInfo_Offset(DataMemberInfo_t* dminfo) const;
   virtual Long_t DataMemberInfo_Property(DataMemberInfo_t* dminfo) const;
   virtual Long_t DataMemberInfo_TypeProperty(DataMemberInfo_t* dminfo) const;
   virtual int    DataMemberInfo_TypeSize(DataMemberInfo_t* dminfo) const;
   virtual const char* DataMemberInfo_TypeName(DataMemberInfo_t* dminfo) const;
   virtual const char* DataMemberInfo_TypeTrueName(DataMemberInfo_t* dminfo) const;
   virtual const char* DataMemberInfo_Name(DataMemberInfo_t* dminfo) const;
   virtual const char* DataMemberInfo_Title(DataMemberInfo_t* dminfo) const;
   virtual const char* DataMemberInfo_ValidArrayIndex(DataMemberInfo_t* dminfo) const;
   virtual void SetDeclAttr(DeclId_t, const char* /* attribute */);


   // Function Template interface
   virtual DeclId_t GetDeclId(FuncTempInfo_t *info) const;
   virtual void   FuncTempInfo_Delete(FuncTempInfo_t * /* ft_info */) const;
   virtual FuncTempInfo_t  *FuncTempInfo_Factory(DeclId_t declid) const;
   virtual FuncTempInfo_t  *FuncTempInfo_FactoryCopy(FuncTempInfo_t * /* ft_info */) const;
   virtual Bool_t FuncTempInfo_IsValid(FuncTempInfo_t * /* ft_info */) const;
   virtual UInt_t FuncTempInfo_TemplateNargs(FuncTempInfo_t * /* ft_info */) const;
   virtual UInt_t FuncTempInfo_TemplateMinReqArgs(FuncTempInfo_t * /* ft_info */) const;
   virtual Long_t FuncTempInfo_Property(FuncTempInfo_t * /* ft_info */) const;
   virtual Long_t FuncTempInfo_ExtraProperty(FuncTempInfo_t * /* ft_info */) const;
   virtual void FuncTempInfo_Name(FuncTempInfo_t * /* ft_info */, TString& name) const;
   virtual void FuncTempInfo_Title(FuncTempInfo_t * /* ft_info */, TString& name) const;

   // MethodInfo interface
   virtual DeclId_t GetDeclId(MethodInfo_t *info) const;
   virtual void   MethodInfo_CreateSignature(MethodInfo_t* minfo, TString& signature) const;
   virtual void   MethodInfo_Delete(MethodInfo_t* minfo) const;
   virtual MethodInfo_t*  MethodInfo_Factory() const;
   virtual MethodInfo_t*  MethodInfo_Factory(ClassInfo_t *clinfo) const;
   virtual MethodInfo_t  *MethodInfo_Factory(DeclId_t declid) const;
   virtual MethodInfo_t*  MethodInfo_FactoryCopy(MethodInfo_t* minfo) const;
   virtual void*  MethodInfo_InterfaceMethod(MethodInfo_t* minfo) const;
   virtual bool   MethodInfo_IsValid(MethodInfo_t* minfo) const;
   virtual int    MethodInfo_NArg(MethodInfo_t* minfo) const;
   virtual int    MethodInfo_NDefaultArg(MethodInfo_t* minfo) const;
   virtual int    MethodInfo_Next(MethodInfo_t* minfo) const;
   virtual Long_t MethodInfo_Property(MethodInfo_t* minfo) const;
   virtual Long_t MethodInfo_ExtraProperty(MethodInfo_t* minfo) const;
   virtual TypeInfo_t*  MethodInfo_Type(MethodInfo_t* minfo) const;
   virtual EReturnType MethodInfo_MethodCallReturnType(MethodInfo_t* minfo) const;
   virtual const char* MethodInfo_GetMangledName(MethodInfo_t* minfo) const;
   virtual const char* MethodInfo_GetPrototype(MethodInfo_t* minfo) const;
   virtual const char* MethodInfo_Name(MethodInfo_t* minfo) const;
   virtual const char* MethodInfo_TypeName(MethodInfo_t* minfo) const;
   virtual std::string MethodInfo_TypeNormalizedName(MethodInfo_t* minfo) const;
   virtual const char* MethodInfo_Title(MethodInfo_t* minfo) const;

   // MethodArgInfo interface
   virtual void   MethodArgInfo_Delete(MethodArgInfo_t* marginfo) const;
   virtual MethodArgInfo_t*  MethodArgInfo_Factory() const;
   virtual MethodArgInfo_t*  MethodArgInfo_Factory(MethodInfo_t *minfo) const;
   virtual MethodArgInfo_t*  MethodArgInfo_FactoryCopy(MethodArgInfo_t* marginfo) const;
   virtual bool   MethodArgInfo_IsValid(MethodArgInfo_t* marginfo) const;
   virtual int    MethodArgInfo_Next(MethodArgInfo_t* marginfo) const;
   virtual Long_t MethodArgInfo_Property(MethodArgInfo_t* marginfo) const;
   virtual const char* MethodArgInfo_DefaultValue(MethodArgInfo_t* marginfo) const;
   virtual const char* MethodArgInfo_Name(MethodArgInfo_t* marginfo) const;
   virtual const char* MethodArgInfo_TypeName(MethodArgInfo_t* marginfo) const;
   virtual std::string MethodArgInfo_TypeNormalizedName(MethodArgInfo_t *marginfo) const;

   // TypeInfo interface
   virtual void   TypeInfo_Delete(TypeInfo_t* tinfo) const;
   virtual TypeInfo_t* TypeInfo_Factory() const;
   virtual TypeInfo_t *TypeInfo_Factory(const char* name) const;
   virtual TypeInfo_t* TypeInfo_FactoryCopy(TypeInfo_t* /* tinfo */) const;
   virtual void   TypeInfo_Init(TypeInfo_t* tinfo, const char* funcname) const;
   virtual bool   TypeInfo_IsValid(TypeInfo_t* tinfo) const;
   virtual const char* TypeInfo_Name(TypeInfo_t* /* info */) const;
   virtual Long_t TypeInfo_Property(TypeInfo_t* tinfo) const;
   virtual int    TypeInfo_RefType(TypeInfo_t* /* tinfo */) const;
   virtual int    TypeInfo_Size(TypeInfo_t* tinfo) const;
   virtual const char* TypeInfo_TrueName(TypeInfo_t* tinfo) const;

   // TypedefInfo interface
   virtual DeclId_t GetDeclId(TypedefInfo_t *info) const;
   virtual void   TypedefInfo_Delete(TypedefInfo_t* tinfo) const;
   virtual TypedefInfo_t*  TypedefInfo_Factory() const;
   virtual TypedefInfo_t*  TypedefInfo_Factory(const char* name) const;
   virtual TypedefInfo_t*  TypedefInfo_FactoryCopy(TypedefInfo_t* tinfo) const;
   virtual void   TypedefInfo_Init(TypedefInfo_t* tinfo, const char* name) const;
   virtual bool   TypedefInfo_IsValid(TypedefInfo_t* tinfo) const;
   virtual int    TypedefInfo_Next(TypedefInfo_t* tinfo) const;
   virtual Long_t TypedefInfo_Property(TypedefInfo_t* tinfo) const;
   virtual int    TypedefInfo_Size(TypedefInfo_t* tinfo) const;
   virtual const char* TypedefInfo_TrueName(TypedefInfo_t* tinfo) const;
   virtual const char* TypedefInfo_Name(TypedefInfo_t* tinfo) const;
   virtual const char* TypedefInfo_Title(TypedefInfo_t* tinfo) const;

   std::set<TClass*>& GetModTClasses() { return fModTClasses; }

   void HandleNewDecl(const void* DV, bool isDeserialized, std::set<TClass*>& modifiedClasses);
   void UpdateListsOnCommitted(const cling::Transaction &T);
   void UpdateListsOnUnloaded(const cling::Transaction &T);
   void InvalidateGlobal(const clang::Decl *D);
   void TransactionRollback(const cling::Transaction &T);
   void LibraryLoaded(const void* dyLibHandle, const char* canonicalName);
   void LibraryUnloaded(const void* dyLibHandle, const char* canonicalName);

private: // Private Utility Functions and Classes
   template <typename List, typename Object>
   static void RemoveAndInvalidateObject(List &L, Object *O) {
      // Invalidate stored information by setting the `xxxInfo_t' to nullptr.
      if (O && O->IsValid())
         L.Unload(O), O->Update(nullptr);
   }

   void InvalidateCachedDecl(const std::tuple<TListOfDataMembers*,
                                        TListOfFunctions*,
                                        TListOfFunctionTemplates*,
                                        TListOfEnums*> &Lists, const clang::Decl *D);

   class SuspendAutoLoadingRAII {
      TCling *fTCling = nullptr;
      bool fOldValue;

   public:
      SuspendAutoLoadingRAII(TCling *tcling) : fTCling(tcling) { fOldValue = fTCling->SetClassAutoLoading(false); }
      ~SuspendAutoLoadingRAII() { fTCling->SetClassAutoLoading(fOldValue); }
   };

   class TUniqueString {
   public:
      TUniqueString() = delete;
      TUniqueString(const TUniqueString &) = delete;
      TUniqueString(Long64_t size);
      const char *Data();
      bool Append(const std::string &str);
   private:
      std::string fContent;
      std::set<size_t> fLinesHashSet;
      std::hash<std::string> fHashFunc;
   };

   TCling();
   TCling(const TCling&) = delete;
   TCling& operator=(const TCling&) = delete;

   void Execute(TMethod*, TObjArray*, int* /*error*/ = nullptr)
   {
   }

   void UpdateListOfLoadedSharedLibraries();
   void RegisterLoadedSharedLibrary(const char* name);
   void AddFriendToClass(clang::FunctionDecl*, clang::CXXRecordDecl*) const;

   std::map<std::string, llvm::StringRef> fPendingRdicts;
   void RegisterRdictForLoadPCM(const std::string &pcmFileNameFullPath, llvm::StringRef *pcmContent);
   void LoadPCM(std::string pcmFileNameFullPath);
   void LoadPCMImpl(TFile &pcmFile);

   void InitRootmapFile(const char *name);
   int  ReadRootmapFile(const char *rootmapfile, TUniqueString* uniqueString = nullptr);
   Bool_t HandleNewTransaction(const cling::Transaction &T);
   bool IsClassAutoLoadingEnabled() const;
   void ProcessClassesToUpdate();
   cling::Interpreter *GetInterpreterImpl() const { return fInterpreter.get(); }
   cling::MetaProcessor *GetMetaProcessorImpl() const { return fMetaProcessor.get(); }

   friend void TCling__RegisterRdictForLoadPCM(const std::string &pcmFileNameFullPath, llvm::StringRef *pcmContent);
   friend cling::Interpreter* TCling__GetInterpreter();
};

#endif
