// @(#)root/meta:$Id$
// Author: Fons Rademakers   01/03/96

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


#ifndef ROOT_TInterpreter
#define ROOT_TInterpreter

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TInterpreter                                                         //
//                                                                      //
// This class defines an abstract interface to a generic command line   //
// interpreter.                                                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TDataType.h"
#include "TDictionary.h"
#include "TInterpreterValue.h"
#include "TNamed.h"
#include "TVirtualRWMutex.h"

#include <map>
#include <typeinfo>
#include <vector>
#include <string>
#include <utility>

class TClass;
class TEnv;
class TFunction;
class TMethod;
class TObjArray;
class TEnum;
class TListOfEnums;
class TSeqCollection;

R__EXTERN TVirtualMutex *gInterpreterMutex;

#if defined (_REENTRANT) || defined (WIN32)
# define R__LOCKGUARD_CLING(mutex)  ::ROOT::Internal::InterpreterMutexRegistrationRAII _R__UNIQUE_(R__guard)(mutex); { }
#else
# define R__LOCKGUARD_CLING(mutex)  (void)(mutex); { }
#endif

namespace ROOT {
namespace Internal {
struct InterpreterMutexRegistrationRAII {
   TLockGuard fLockGuard;
   InterpreterMutexRegistrationRAII(TVirtualMutex* mutex);
   ~InterpreterMutexRegistrationRAII();
};
}
}

class TInterpreter : public TNamed {

protected:
   virtual void Execute(TMethod *method, TObjArray *params, int *error = 0) = 0;
   virtual Bool_t SetSuspendAutoParsing(Bool_t value) = 0;

   friend class SuspendAutoParsing;

public:
   // See as in TSchemaType.h.
   typedef class std::map<std::string, std::string> MembersMap_t;

   enum EErrorCode {
      kNoError     = 0,
      kRecoverable = 1,
      kDangerous   = 2,
      kFatal       = 3,
      kProcessing  = 99
   };

   enum class EReturnType { kLong, kDouble, kString, kOther, kNoReturnType };

   struct CallFuncIFacePtr_t {
      enum EKind {
         kUninitialized,
         kGeneric,
         kCtor,
         kDtor
      };

      typedef void (*Generic_t)(void*, int, void**, void*);
      typedef void (*Ctor_t)(void**, void*, unsigned long);
      typedef void (*Dtor_t)(void*, unsigned long, int);

      CallFuncIFacePtr_t():
         fKind(kUninitialized), fGeneric(0) {}
      CallFuncIFacePtr_t(Generic_t func):
         fKind(kGeneric), fGeneric(func) {}
      CallFuncIFacePtr_t(Ctor_t func):
         fKind(kCtor), fCtor(func) {}
      CallFuncIFacePtr_t(Dtor_t func):
         fKind(kDtor), fDtor(func) {}

      EKind fKind;
      union {
         Generic_t fGeneric;
         Ctor_t fCtor;
         Dtor_t fDtor;
      };
   };

   class SuspendAutoParsing {
      TInterpreter *fInterp;
      Bool_t        fPrevious;
   public:
      SuspendAutoParsing(TInterpreter *where, Bool_t value = kTRUE) : fInterp(where), fPrevious(fInterp->SetSuspendAutoParsing(value)) {}
      ~SuspendAutoParsing() { fInterp->SetSuspendAutoParsing(fPrevious); }
   };
   virtual Bool_t IsAutoParsingSuspended() const = 0;

   class SuspendAutoLoadingRAII {
      TInterpreter *fInterp = nullptr;
      bool fOldValue;

   public:
      SuspendAutoLoadingRAII(TInterpreter *interp) : fInterp(interp)
      {
         fOldValue = fInterp->SetClassAutoLoading(false);
      }
      ~SuspendAutoLoadingRAII() { fInterp->SetClassAutoLoading(fOldValue); }
   };

   typedef int (*AutoLoadCallBack_t)(const char*);
   typedef std::vector<std::pair<std::string, int> > FwdDeclArgsToKeepCollection_t;

   TInterpreter() { }   // for Dictionary
   TInterpreter(const char *name, const char *title = "Generic Interpreter");
   virtual ~TInterpreter() { }

   virtual void     AddIncludePath(const char *path) = 0;
   virtual void    *SetAutoLoadCallBack(void* /*cb*/) { return 0; }
   virtual void    *GetAutoLoadCallBack() const { return 0; }
   virtual Int_t    AutoLoad(const char *classname, Bool_t knowDictNotLoaded = kFALSE) = 0;
   virtual Int_t    AutoLoad(const std::type_info& typeinfo, Bool_t knowDictNotLoaded = kFALSE) = 0;
   virtual Int_t    AutoParse(const char* cls) = 0;
   virtual void     ClearFileBusy() = 0;
   virtual void     ClearStack() = 0; // Delete existing temporary values
   virtual Bool_t   Declare(const char* code) = 0;
   virtual void     EndOfLineAction() = 0;
   virtual TClass  *GetClass(const std::type_info& typeinfo, Bool_t load) const = 0;
   virtual Int_t    GetExitCode() const = 0;
   virtual TEnv    *GetMapfile() const { return 0; }
   virtual Int_t    GetMore() const = 0;
   virtual TClass  *GenerateTClass(const char *classname, Bool_t emulation, Bool_t silent = kFALSE) = 0;
   virtual TClass  *GenerateTClass(ClassInfo_t *classinfo, Bool_t silent = kFALSE) = 0;
   virtual Int_t    GenerateDictionary(const char *classes, const char *includes = 0, const char *options = 0) = 0;
   virtual char    *GetPrompt() = 0;
   virtual const char *GetSharedLibs() = 0;
   virtual const char *GetClassSharedLibs(const char *cls) = 0;
   virtual const char *GetSharedLibDeps(const char *lib, bool tryDyld = false) = 0;
   virtual const char *GetIncludePath() = 0;
   virtual const char *GetSTLIncludePath() const { return ""; }
   virtual TObjArray  *GetRootMapFiles() const = 0;
   virtual void     Initialize() = 0;
   virtual void     ShutDown() = 0;
   virtual void     InspectMembers(TMemberInspector&, const void* obj, const TClass* cl, Bool_t isTransient) = 0;
   virtual Bool_t   IsLoaded(const char *filename) const = 0;
   virtual Bool_t   IsLibraryLoaded(const char *libname) const = 0;
   virtual Bool_t   HasPCMForLibrary(const char *libname) const = 0;
   virtual Int_t    Load(const char *filenam, Bool_t system = kFALSE) = 0;
   virtual void     LoadMacro(const char *filename, EErrorCode *error = 0) = 0;
   virtual Int_t    LoadLibraryMap(const char *rootmapfile = 0) = 0;
   virtual Int_t    RescanLibraryMap() = 0;
   virtual Int_t    ReloadAllSharedLibraryMaps() = 0;
   virtual Int_t    UnloadAllSharedLibraryMaps() = 0;
   virtual Int_t    UnloadLibraryMap(const char *library) = 0;
   virtual Long_t   ProcessLine(const char *line, EErrorCode *error = 0) = 0;
   virtual Long_t   ProcessLineSynch(const char *line, EErrorCode *error = 0) = 0;
   virtual void     PrintIntro() = 0;
   virtual bool     RegisterPrebuiltModulePath(const std::string& FullPath,
                                               const std::string& ModuleMapName = "module.modulemap") const = 0;
   virtual void     RegisterModule(const char* /*modulename*/,
                                   const char** /*headers*/,
                                   const char** /*includePaths*/,
                                   const char* /*payloadCode*/,
                                   const char* /*fwdDeclsCode*/,
                                   void (* /*triggerFunc*/)(),
                                   const FwdDeclArgsToKeepCollection_t& fwdDeclArgsToKeep,
                                   const char** classesHeaders,
                                   Bool_t lateRegistration = false,
                                   Bool_t hasCxxModule = false) = 0;
   virtual void     AddAvailableIndentifiers(TSeqCollection&) = 0;
   virtual void     RegisterTClassUpdate(TClass *oldcl,DictFuncPtr_t dict) = 0;
   virtual void     UnRegisterTClassUpdate(const TClass *oldcl) = 0;
   virtual Int_t    SetClassSharedLibs(const char *cls, const char *libs) = 0;
   virtual void     SetGetline(const char*(*getlineFunc)(const char* prompt),
                               void (*histaddFunc)(const char* line)) = 0;
   virtual void     Reset() = 0;
   virtual void     ResetAll() = 0;
   virtual void     ResetGlobals() = 0;
   virtual void     ResetGlobalVar(void *obj) = 0;
   virtual void     RewindDictionary() = 0;
   virtual Int_t    DeleteGlobal(void *obj) = 0;
   virtual Int_t    DeleteVariable(const char* name) = 0;
   virtual void     SaveContext() = 0;
   virtual void     SaveGlobalsContext() = 0;
   virtual void     UpdateListOfGlobals() = 0;
   virtual void     UpdateListOfGlobalFunctions() = 0;
   virtual void     UpdateListOfTypes() = 0;
   virtual void     SetClassInfo(TClass *cl, Bool_t reload = kFALSE) = 0;

   enum ECheckClassInfo {
      kUnknown = 0, // backward compatible with false
      kKnown = 1,
      kWithClassDefInline = 2
   };
   virtual ECheckClassInfo CheckClassInfo(const char *name, Bool_t autoload, Bool_t isClassOrNamespaceOnly = kFALSE) = 0;

   virtual Bool_t   CheckClassTemplate(const char *name) = 0;
   virtual Long_t   Calc(const char *line, EErrorCode* error = 0) = 0;
   virtual void     CreateListOfBaseClasses(TClass *cl) const = 0;
   virtual void     CreateListOfDataMembers(TClass *cl) const = 0;
   virtual void     CreateListOfMethods(TClass *cl) const = 0;
   virtual void     CreateListOfMethodArgs(TFunction *m) const = 0;
   virtual void     UpdateListOfMethods(TClass *cl) const = 0;
   virtual TString  GetMangledName(TClass *cl, const char *method, const char *params, Bool_t objectIsConst = kFALSE) = 0;
   virtual TString  GetMangledNameWithPrototype(TClass *cl, const char *method, const char *proto, Bool_t objectIsConst = kFALSE, ROOT::EFunctionMatchMode /* mode */ = ROOT::kConversionMatch) = 0;
   virtual void     GetInterpreterTypeName(const char *name, std::string &output, Bool_t full = kFALSE) = 0;
   virtual void    *GetInterfaceMethod(TClass *cl, const char *method, const char *params, Bool_t objectIsConst = kFALSE) = 0;
   virtual void    *GetInterfaceMethodWithPrototype(TClass *cl, const char *method, const char *proto, Bool_t objectIsConst = kFALSE, ROOT::EFunctionMatchMode /* mode */ = ROOT::kConversionMatch) = 0;
   virtual void     Execute(const char *function, const char *params, int *error = 0) = 0;
   virtual void     Execute(TObject *obj, TClass *cl, const char *method, const char *params, int *error = 0) = 0;
   virtual void     Execute(TObject *obj, TClass *cl, TMethod *method, TObjArray *params, int *error = 0) = 0;
   virtual void     ExecuteWithArgsAndReturn(TMethod *method, void* address, const void* args[] = 0, int /*nargs*/ = 0, void* ret= 0) const = 0;
   virtual Long_t   ExecuteMacro(const char *filename, EErrorCode *error = 0) = 0;
   virtual Bool_t   IsErrorMessagesEnabled() const = 0;
   virtual Bool_t   SetErrorMessages(Bool_t enable = kTRUE) = 0;
   virtual Bool_t   IsProcessLineLocked() const = 0;
   virtual void     SetProcessLineLock(Bool_t lock = kTRUE) = 0;
   virtual const char *TypeName(const char *s) = 0;
   virtual std::string ToString(const char *type, void *obj) = 0;

   virtual void     SnapshotMutexState(ROOT::TVirtualRWMutex* mtx) = 0;
   virtual void     ForgetMutexState() = 0;

   // All the functions below must be virtual with a dummy implementation
   // These functions are redefined in TCling.

   // Misc
   virtual int    DisplayClass(FILE * /* fout */,const char * /* name */,int /* base */,int /* start */) const {return 0;}
   virtual int    DisplayIncludePath(FILE * /* fout */) const {return 0;}
   virtual void  *FindSym(const char * /* entry */) const {return 0;}
   virtual void   GenericError(const char * /* error */) const {;}
   virtual Long_t GetExecByteCode() const {return 0;}
   virtual const char *GetTopLevelMacroName() const {return 0;};
   virtual const char *GetCurrentMacroName()  const {return 0;};
   virtual int    GetSecurityError() const{return 0;}
   virtual int    LoadFile(const char * /* path */) const {return 0;}
   virtual Bool_t LoadText(const char * /* text */) const {return kFALSE;}
   virtual const char *MapCppName(const char*) const {return 0;}
   virtual void   SetAlloclockfunc(void (*)()) const {;}
   virtual void   SetAllocunlockfunc(void (*)()) const {;}
   virtual int    SetClassAutoLoading(int) const {return 0;}
           int    SetClassAutoloading(int a) const { return SetClassAutoLoading(a); }  // Deprecated
   virtual int    SetClassAutoparsing(int) {return 0;};
   virtual void   SetErrmsgcallback(void * /* p */) const {;}
   virtual void   SetTempLevel(int /* val */) const {;}
   virtual int    UnloadFile(const char * /* path */) const {return 0;}

   /// The created temporary must be deleted by the caller.
   /// Deprecated! Please use MakeInterpreterValue().
   TInterpreterValue *CreateTemporary() const {
      return MakeInterpreterValue().release();
   }
   virtual std::unique_ptr<TInterpreterValue> MakeInterpreterValue() const { return 0; }
   virtual void   CodeComplete(const std::string&, size_t&,
                               std::vector<std::string>&) {;}
   virtual int Evaluate(const char*, TInterpreterValue&) {return 0;}

   // core/meta helper functions.
   virtual EReturnType MethodCallReturnType(TFunction *func) const = 0;
   virtual ULong64_t GetInterpreterStateMarker() const = 0;
   virtual bool DiagnoseIfInterpreterException(const std::exception &e) const = 0;

   typedef TDictionary::DeclId_t DeclId_t;
   virtual DeclId_t GetDeclId(CallFunc_t *info) const = 0;
   virtual DeclId_t GetDeclId(ClassInfo_t *info) const = 0;
   virtual DeclId_t GetDeclId(DataMemberInfo_t *info) const = 0;
   virtual DeclId_t GetDeclId(FuncTempInfo_t *info) const = 0;
   virtual DeclId_t GetDeclId(MethodInfo_t *info) const = 0;
   virtual DeclId_t GetDeclId(TypedefInfo_t *info) const = 0;

   virtual void SetDeclAttr(DeclId_t, const char* /* attribute */) = 0 ;

   virtual DeclId_t GetDataMember(ClassInfo_t *cl, const char *name) const = 0;
   virtual DeclId_t GetDataMemberAtAddr(const void *addr) const = 0;
   virtual DeclId_t GetDataMemberWithValue(const void *ptrvalue) const = 0;
   virtual DeclId_t GetEnum(TClass *cl, const char *name) const = 0;
   virtual TEnum*   CreateEnum(void *VD, TClass *cl) const = 0;
   virtual void     UpdateEnumConstants(TEnum* enumObj, TClass* cl) const = 0;
   virtual void     LoadEnums(TListOfEnums& cl) const = 0;
   virtual DeclId_t GetFunction(ClassInfo_t *cl, const char *funcname) = 0;
   virtual DeclId_t GetFunctionWithPrototype(ClassInfo_t *cl, const char* method, const char* proto, Bool_t objectIsConst = kFALSE, ROOT::EFunctionMatchMode mode = ROOT::kConversionMatch) = 0;
   virtual DeclId_t GetFunctionWithValues(ClassInfo_t *cl, const char* method, const char* params, Bool_t objectIsConst = kFALSE) = 0;
   virtual DeclId_t GetFunctionTemplate(ClassInfo_t *cl, const char *funcname) = 0;
   virtual void     GetFunctionOverloads(ClassInfo_t *cl, const char *funcname, std::vector<DeclId_t>& res) const = 0;
   virtual void     LoadFunctionTemplates(TClass* cl) const = 0;
   virtual std::vector<std::string> GetUsingNamespaces(ClassInfo_t *cl) const = 0;

   // CallFunc interface
   virtual void   CallFunc_Delete(CallFunc_t * /* func */) const {;}
   virtual void   CallFunc_Exec(CallFunc_t * /* func */, void * /* address */) const {;}
   virtual void   CallFunc_Exec(CallFunc_t * /* func */, void * /* address */, TInterpreterValue& /* val */) const {;}
   virtual void   CallFunc_ExecWithReturn(CallFunc_t * /* func */, void * /* address */, void * /* ret */) const {;}
   virtual void   CallFunc_ExecWithArgsAndReturn(CallFunc_t * /* func */, void * /* address */, const void* /* args */ [] = 0, int /*nargs*/ = 0, void * /* ret */ = 0) const {}
   virtual Long_t    CallFunc_ExecInt(CallFunc_t * /* func */, void * /* address */) const {return 0;}
   virtual Long64_t  CallFunc_ExecInt64(CallFunc_t * /* func */, void * /* address */) const {return 0;}
   virtual Double_t  CallFunc_ExecDouble(CallFunc_t * /* func */, void * /* address */) const {return 0;}
   virtual CallFunc_t   *CallFunc_Factory() const {return 0;}
   virtual CallFunc_t   *CallFunc_FactoryCopy(CallFunc_t * /* func */) const {return 0;}
   virtual MethodInfo_t *CallFunc_FactoryMethod(CallFunc_t * /* func */) const {return 0;}
   virtual void   CallFunc_IgnoreExtraArgs(CallFunc_t * /*func */, bool /*ignore*/) const {;}
   virtual void   CallFunc_Init(CallFunc_t * /* func */) const {;}
   virtual Bool_t CallFunc_IsValid(CallFunc_t * /* func */) const {return 0;}
   virtual CallFuncIFacePtr_t CallFunc_IFacePtr(CallFunc_t * /* func */) const {return CallFuncIFacePtr_t();}
   virtual void   CallFunc_ResetArg(CallFunc_t * /* func */) const {;}
   virtual void   CallFunc_SetArgArray(CallFunc_t * /* func */, Long_t * /* paramArr */, Int_t /* nparam */) const {;}
   virtual void   CallFunc_SetArgs(CallFunc_t * /* func */, const char * /* param */) const {;}

   virtual void   CallFunc_SetArg(CallFunc_t * /*func */, Long_t /* param */) const = 0;
   virtual void   CallFunc_SetArg(CallFunc_t * /*func */, ULong_t /* param */) const = 0;
   virtual void   CallFunc_SetArg(CallFunc_t * /* func */, Float_t /* param */) const = 0;
   virtual void   CallFunc_SetArg(CallFunc_t * /* func */, Double_t /* param */) const = 0;
   virtual void   CallFunc_SetArg(CallFunc_t * /* func */, Long64_t /* param */) const = 0;
   virtual void   CallFunc_SetArg(CallFunc_t * /* func */, ULong64_t /* param */) const = 0;

   void CallFunc_SetArg(CallFunc_t * func, Char_t param) const { CallFunc_SetArg(func,(Long_t)param); }
   void CallFunc_SetArg(CallFunc_t * func, Short_t param) const { CallFunc_SetArg(func,(Long_t)param); }
   void CallFunc_SetArg(CallFunc_t * func, Int_t param) const { CallFunc_SetArg(func,(Long_t)param); }

   void CallFunc_SetArg(CallFunc_t * func, UChar_t param) const { CallFunc_SetArg(func,(ULong_t)param); }
   void CallFunc_SetArg(CallFunc_t * func, UShort_t param) const { CallFunc_SetArg(func,(ULong_t)param); }
   void CallFunc_SetArg(CallFunc_t * func, UInt_t param) const { CallFunc_SetArg(func,(ULong_t)param); }

   template <typename T>
   void CallFunc_SetArgRef(CallFunc_t * func, T &param) const { CallFunc_SetArg(func,(ULong_t)&param); }

   void CallFunc_SetArg(CallFunc_t *func, void *arg)
   {
      CallFunc_SetArg(func,(Long_t) arg);
   }

   template <typename T>
   void CallFunc_SetArg(CallFunc_t *func, const T *arg)
   {
      CallFunc_SetArg(func,(Long_t) arg);
   }

   void CallFunc_SetArgImpl(CallFunc_t * /* func */)
   {
   }

   template <typename U>
   void CallFunc_SetArgImpl(CallFunc_t *func, const U& head)
   {
      CallFunc_SetArg(func, head);
   }

   template <typename U, typename... T>
   void CallFunc_SetArgImpl(CallFunc_t *func, const U& head, const T&... tail)
   {
      CallFunc_SetArg(func, head);
      CallFunc_SetArgImpl(func, tail...);
   }

   template <typename... T>
   void CallFunc_SetArguments(CallFunc_t *func, const T&... args)
   {
      R__LOCKGUARD(gInterpreterMutex);

      CallFunc_ResetArg(func);
      CallFunc_SetArgImpl(func,args...);
   }

   virtual void   CallFunc_SetFunc(CallFunc_t * /* func */, ClassInfo_t * /* info */, const char * /* method */, const char * /* params */, bool /* objectIsConst */, Long_t * /* Offset */) const {;}
   virtual void   CallFunc_SetFunc(CallFunc_t * /* func */, ClassInfo_t * /* info */, const char * /* method */, const char * /* params */, Long_t * /* Offset */) const {;}
   virtual void   CallFunc_SetFunc(CallFunc_t * /* func */, MethodInfo_t * /* info */) const {;}
   virtual void   CallFunc_SetFuncProto(CallFunc_t * /* func */, ClassInfo_t * /* info */, const char * /* method */, const char * /* proto */, Long_t * /* Offset */, ROOT::EFunctionMatchMode /* mode */ = ROOT::kConversionMatch) const {;}
   virtual void   CallFunc_SetFuncProto(CallFunc_t * /* func */, ClassInfo_t * /* info */, const char * /* method */, const char * /* proto */, bool /* objectIsConst */, Long_t * /* Offset */, ROOT::EFunctionMatchMode /* mode */ = ROOT::kConversionMatch) const {;}
   virtual void   CallFunc_SetFuncProto(CallFunc_t* func, ClassInfo_t* info, const char* method, const std::vector<TypeInfo_t*> &proto, Long_t* Offset, ROOT::EFunctionMatchMode mode = ROOT::kConversionMatch) const = 0;
   virtual void   CallFunc_SetFuncProto(CallFunc_t* func, ClassInfo_t* info, const char* method, const std::vector<TypeInfo_t*> &proto, bool objectIsConst, Long_t* Offset, ROOT::EFunctionMatchMode mode = ROOT::kConversionMatch) const = 0;

   virtual std::string CallFunc_GetWrapperCode(CallFunc_t *func) const = 0;

   // ClassInfo interface
   virtual Bool_t ClassInfo_Contains(ClassInfo_t *info, DeclId_t decl) const = 0;
   virtual Long_t ClassInfo_ClassProperty(ClassInfo_t * /* info */) const {return 0;}
   virtual void   ClassInfo_Delete(ClassInfo_t * /* info */) const {;}
   virtual void   ClassInfo_Delete(ClassInfo_t * /* info */, void * /* arena */) const {;}
   virtual void   ClassInfo_DeleteArray(ClassInfo_t * /* info */, void * /* arena */, bool /* dtorOnly */) const {;}
   virtual void   ClassInfo_Destruct(ClassInfo_t * /* info */, void * /* arena */) const {;}
   virtual ClassInfo_t  *ClassInfo_Factory(Bool_t /*all*/ = kTRUE) const = 0;
   virtual ClassInfo_t  *ClassInfo_Factory(ClassInfo_t * /* cl */) const = 0;
   virtual ClassInfo_t  *ClassInfo_Factory(const char * /* name */) const = 0;
   virtual ClassInfo_t  *ClassInfo_Factory(DeclId_t declid) const = 0;
   virtual Long_t   ClassInfo_GetBaseOffset(ClassInfo_t* /* fromDerived */,
                                            ClassInfo_t* /* toBase */, void* /* address */ = 0, bool /*isderived*/ = true) const {return 0;}
   virtual int    ClassInfo_GetMethodNArg(ClassInfo_t * /* info */, const char * /* method */,const char * /* proto */, Bool_t /* objectIsConst */ = false, ROOT::EFunctionMatchMode /* mode */ = ROOT::kConversionMatch) const {return 0;}
   virtual Bool_t ClassInfo_HasDefaultConstructor(ClassInfo_t * /* info */, Bool_t = kFALSE) const {return kFALSE;}
   virtual Bool_t ClassInfo_HasMethod(ClassInfo_t * /* info */, const char * /* name */) const {return kFALSE;}
   virtual void   ClassInfo_Init(ClassInfo_t * /* info */, const char * /* funcname */) const {;}
   virtual void   ClassInfo_Init(ClassInfo_t * /* info */, int /* tagnum */) const {;}
   virtual Bool_t ClassInfo_IsBase(ClassInfo_t * /* info */, const char * /* name */) const {return 0;}
   virtual Bool_t ClassInfo_IsEnum(const char * /* name */) const {return 0;}
   virtual Bool_t ClassInfo_IsScopedEnum(ClassInfo_t * /* info */) const {return 0;}
   virtual EDataType ClassInfo_GetUnderlyingType(ClassInfo_t * /* info */) const {return kNumDataTypes;}
   virtual Bool_t ClassInfo_IsLoaded(ClassInfo_t * /* info */) const {return 0;}
   virtual Bool_t ClassInfo_IsValid(ClassInfo_t * /* info */) const {return 0;}
   virtual Bool_t ClassInfo_IsValidMethod(ClassInfo_t * /* info */, const char * /* method */,const char * /* proto */, Long_t * /* offset */, ROOT::EFunctionMatchMode /* mode */ = ROOT::kConversionMatch) const {return 0;}
   virtual Bool_t ClassInfo_IsValidMethod(ClassInfo_t * /* info */, const char * /* method */,const char * /* proto */, Bool_t /* objectIsConst */, Long_t * /* offset */, ROOT::EFunctionMatchMode /* mode */ = ROOT::kConversionMatch) const {return 0;}
   virtual int    ClassInfo_Next(ClassInfo_t * /* info */) const {return 0;}
   virtual void  *ClassInfo_New(ClassInfo_t * /* info */) const {return 0;}
   virtual void  *ClassInfo_New(ClassInfo_t * /* info */, int /* n */) const {return 0;}
   virtual void  *ClassInfo_New(ClassInfo_t * /* info */, int /* n */, void * /* arena */) const {return 0;}
   virtual void  *ClassInfo_New(ClassInfo_t * /* info */, void * /* arena */) const {return 0;}
   virtual Long_t ClassInfo_Property(ClassInfo_t * /* info */) const {return 0;}
   virtual int    ClassInfo_Size(ClassInfo_t * /* info */) const {return 0;}
   virtual Long_t ClassInfo_Tagnum(ClassInfo_t * /* info */) const {return 0;}
   virtual const char *ClassInfo_FileName(ClassInfo_t * /* info */) const {return 0;}
   virtual const char *ClassInfo_FullName(ClassInfo_t * /* info */) const {return 0;}
   virtual const char *ClassInfo_Name(ClassInfo_t * /* info */) const {return 0;}
   virtual const char *ClassInfo_Title(ClassInfo_t * /* info */) const {return 0;}
   virtual const char *ClassInfo_TmpltName(ClassInfo_t * /* info */) const {return 0;}


   // BaseClassInfo interface
   virtual void   BaseClassInfo_Delete(BaseClassInfo_t * /* bcinfo */) const {;}
   virtual BaseClassInfo_t  *BaseClassInfo_Factory(ClassInfo_t * /* info */) const {return 0;}
   virtual BaseClassInfo_t  *BaseClassInfo_Factory(ClassInfo_t* /* derived */,
                                                   ClassInfo_t* /* base */) const {return 0;}
   virtual int    BaseClassInfo_Next(BaseClassInfo_t * /* bcinfo */) const {return 0;}
   virtual int    BaseClassInfo_Next(BaseClassInfo_t * /* bcinfo */, int  /* onlyDirect */) const {return 0;}
   virtual Long_t BaseClassInfo_Offset(BaseClassInfo_t * /* toBaseClassInfo */, void* /* address */ = 0 /*default for non-virtual*/, bool /*isderived*/ = true /*default for non-virtual*/) const {return 0;}
   virtual Long_t BaseClassInfo_Property(BaseClassInfo_t * /* bcinfo */) const {return 0;}
   virtual Long_t BaseClassInfo_Tagnum(BaseClassInfo_t * /* bcinfo */) const {return 0;}
   virtual ClassInfo_t*BaseClassInfo_ClassInfo(BaseClassInfo_t * /* bcinfo */) const = 0;
   virtual const char *BaseClassInfo_FullName(BaseClassInfo_t * /* bcinfo */) const {return 0;}
   virtual const char *BaseClassInfo_Name(BaseClassInfo_t * /* bcinfo */) const {return 0;}
   virtual const char *BaseClassInfo_TmpltName(BaseClassInfo_t * /* bcinfo */) const {return 0;}

   // DataMemberInfo interface
   virtual int    DataMemberInfo_ArrayDim(DataMemberInfo_t * /* dminfo */) const {return 0;}
   virtual void   DataMemberInfo_Delete(DataMemberInfo_t * /* dminfo */) const {;}
   virtual DataMemberInfo_t  *DataMemberInfo_Factory(ClassInfo_t * /* clinfo */, TDictionary::EMemberSelection /*selection*/) const {return 0;}
   virtual DataMemberInfo_t  *DataMemberInfo_Factory(DeclId_t declid, ClassInfo_t* clinfo) const = 0;
   virtual DataMemberInfo_t  *DataMemberInfo_FactoryCopy(DataMemberInfo_t * /* dminfo */) const {return 0;}
   virtual Bool_t DataMemberInfo_IsValid(DataMemberInfo_t * /* dminfo */) const {return 0;}
   virtual int    DataMemberInfo_MaxIndex(DataMemberInfo_t * /* dminfo */, Int_t  /* dim */) const {return 0;}
   virtual int    DataMemberInfo_Next(DataMemberInfo_t * /* dminfo */) const {return 0;}
   virtual Long_t DataMemberInfo_Offset(DataMemberInfo_t * /* dminfo */) const {return 0;}
   virtual Long_t DataMemberInfo_Property(DataMemberInfo_t * /* dminfo */) const {return 0;}
   virtual Long_t DataMemberInfo_TypeProperty(DataMemberInfo_t * /* dminfo */) const {return 0;}
   virtual int    DataMemberInfo_TypeSize(DataMemberInfo_t * /* dminfo */) const {return 0;}
   virtual const char *DataMemberInfo_TypeName(DataMemberInfo_t * /* dminfo */) const {return 0;}
   virtual const char *DataMemberInfo_TypeTrueName(DataMemberInfo_t * /* dminfo */) const {return 0;}
   virtual const char *DataMemberInfo_Name(DataMemberInfo_t * /* dminfo */) const {return 0;}
   virtual const char *DataMemberInfo_Title(DataMemberInfo_t * /* dminfo */) const {return 0;}
   virtual const char *DataMemberInfo_ValidArrayIndex(DataMemberInfo_t * /* dminfo */) const {return 0;}

   // Function Template interface
   virtual void   FuncTempInfo_Delete(FuncTempInfo_t * /* ft_info */) const = 0;
   virtual FuncTempInfo_t  *FuncTempInfo_Factory(DeclId_t declid) const = 0;
   virtual FuncTempInfo_t  *FuncTempInfo_FactoryCopy(FuncTempInfo_t * /* ft_info */) const = 0;
   virtual Bool_t FuncTempInfo_IsValid(FuncTempInfo_t * /* ft_info */) const = 0;
   virtual UInt_t FuncTempInfo_TemplateNargs(FuncTempInfo_t * /* ft_info */) const = 0;
   virtual UInt_t FuncTempInfo_TemplateMinReqArgs(FuncTempInfo_t * /* ft_info */) const = 0;
   virtual Long_t FuncTempInfo_Property(FuncTempInfo_t * /* ft_info */) const = 0;
   virtual Long_t FuncTempInfo_ExtraProperty(FuncTempInfo_t * /* ft_info */) const = 0;
   virtual void FuncTempInfo_Name(FuncTempInfo_t * /* ft_info */, TString &name) const = 0;
   virtual void FuncTempInfo_Title(FuncTempInfo_t * /* ft_info */, TString &title) const = 0;

   // MethodInfo interface
   virtual void   MethodInfo_CreateSignature(MethodInfo_t * /* minfo */, TString & /* signature */) const {;}
   virtual void   MethodInfo_Delete(MethodInfo_t * /* minfo */) const {;}
   virtual MethodInfo_t  *MethodInfo_Factory() const {return 0;}
   virtual MethodInfo_t  *MethodInfo_Factory(ClassInfo_t * /*clinfo*/) const {return 0;}
   virtual MethodInfo_t  *MethodInfo_Factory(DeclId_t declid) const = 0;
   virtual MethodInfo_t  *MethodInfo_FactoryCopy(MethodInfo_t * /* minfo */) const {return 0;}
   virtual void  *MethodInfo_InterfaceMethod(MethodInfo_t * /* minfo */) const {return 0;}
   virtual Bool_t MethodInfo_IsValid(MethodInfo_t * /* minfo */) const {return 0;}
   virtual int    MethodInfo_NArg(MethodInfo_t * /* minfo */) const {return 0;}
   virtual int    MethodInfo_NDefaultArg(MethodInfo_t * /* minfo */) const {return 0;}
   virtual int    MethodInfo_Next(MethodInfo_t * /* minfo */) const {return 0;}
   virtual Long_t MethodInfo_Property(MethodInfo_t * /* minfo */) const = 0;
   virtual Long_t MethodInfo_ExtraProperty(MethodInfo_t * /* minfo */) const = 0;
   virtual TypeInfo_t  *MethodInfo_Type(MethodInfo_t * /* minfo */) const {return 0;}
   virtual EReturnType MethodInfo_MethodCallReturnType(MethodInfo_t* minfo) const = 0;
   virtual const char *MethodInfo_GetMangledName(MethodInfo_t * /* minfo */) const {return 0;}
   virtual const char *MethodInfo_GetPrototype(MethodInfo_t * /* minfo */) const {return 0;}
   virtual const char *MethodInfo_Name(MethodInfo_t * /* minfo */) const {return 0;}
   virtual const char *MethodInfo_TypeName(MethodInfo_t * /* minfo */) const {return 0;}
   virtual std::string MethodInfo_TypeNormalizedName(MethodInfo_t * /* minfo */) const {return "";}
   virtual const char *MethodInfo_Title(MethodInfo_t * /* minfo */) const {return 0;}

   // MethodArgInfo interface
   virtual void   MethodArgInfo_Delete(MethodArgInfo_t * /* marginfo */) const {;}
   virtual MethodArgInfo_t  *MethodArgInfo_Factory() const {return 0;}
   virtual MethodArgInfo_t  *MethodArgInfo_Factory(MethodInfo_t * /*minfo*/) const {return 0;}
   virtual MethodArgInfo_t  *MethodArgInfo_FactoryCopy(MethodArgInfo_t * /* marginfo */) const {return 0;}
   virtual Bool_t MethodArgInfo_IsValid(MethodArgInfo_t * /* marginfo */) const {return 0;}
   virtual int    MethodArgInfo_Next(MethodArgInfo_t * /* marginfo */) const {return 0;}
   virtual Long_t MethodArgInfo_Property(MethodArgInfo_t * /* marginfo */) const {return 0;}
   virtual const char *MethodArgInfo_DefaultValue(MethodArgInfo_t * /* marginfo */) const {return 0;}
   virtual const char *MethodArgInfo_Name(MethodArgInfo_t * /* marginfo */) const {return 0;}
   virtual const char *MethodArgInfo_TypeName(MethodArgInfo_t * /* marginfo */) const {return 0;}
   virtual std::string MethodArgInfo_TypeNormalizedName(MethodArgInfo_t * /* marginfo */) const = 0;


   // TypeInfo interface
   virtual void    TypeInfo_Delete(TypeInfo_t * /* tinfo */) const {;}
   virtual TypeInfo_t *TypeInfo_Factory() const {return 0;}
   virtual TypeInfo_t *TypeInfo_Factory(const char* /* name */) const {return 0;}
   virtual TypeInfo_t *TypeInfo_FactoryCopy(TypeInfo_t * /* tinfo */) const {return 0;}
   virtual void   TypeInfo_Init(TypeInfo_t * /* tinfo */, const char * /* funcname */) const {;}
   virtual Bool_t TypeInfo_IsValid(TypeInfo_t * /* tinfo */) const {return 0;}
   virtual const char *TypeInfo_Name(TypeInfo_t * /* info */) const {return 0;}
   virtual Long_t TypeInfo_Property(TypeInfo_t * /* tinfo */) const {return 0;}
   virtual int    TypeInfo_RefType(TypeInfo_t * /* tinfo */) const {return 0;}
   virtual int    TypeInfo_Size(TypeInfo_t * /* tinfo */) const {return 0;}
   virtual const char *TypeInfo_TrueName(TypeInfo_t * /* tinfo */) const {return 0;}


   // TypedefInfo interface
   virtual void   TypedefInfo_Delete(TypedefInfo_t * /* tinfo */) const {;}
   virtual TypedefInfo_t  *TypedefInfo_Factory() const {return 0;}
   virtual TypedefInfo_t  *TypedefInfo_Factory(const char *) const {return 0;}
   virtual TypedefInfo_t  *TypedefInfo_FactoryCopy(TypedefInfo_t * /* tinfo */) const {return 0;}
   virtual void   TypedefInfo_Init(TypedefInfo_t * /* tinfo */, const char * /* funcname */) const {;}
   virtual Bool_t TypedefInfo_IsValid(TypedefInfo_t * /* tinfo */) const {return 0;}
   virtual int    TypedefInfo_Next(TypedefInfo_t* /*tinfo*/) const {return 0;}
   virtual Long_t TypedefInfo_Property(TypedefInfo_t * /* tinfo */) const {return 0;}
   virtual int    TypedefInfo_Size(TypedefInfo_t * /* tinfo */) const {return 0;}
   virtual const char *TypedefInfo_TrueName(TypedefInfo_t * /* tinfo */) const {return 0;}
   virtual const char *TypedefInfo_Name(TypedefInfo_t * /* tinfo */) const {return 0;}
   virtual const char *TypedefInfo_Title(TypedefInfo_t * /* tinfo */) const {return 0;}

   static TInterpreter *Instance();

   ClassDef(TInterpreter,0)  //ABC defining interface to generic interpreter
};


typedef TInterpreter *CreateInterpreter_t(void* shlibHandle, const char* argv[]);
typedef void *DestroyInterpreter_t(TInterpreter*);

#ifndef __CINT__
#define gInterpreter (TInterpreter::Instance())
R__EXTERN TInterpreter* gCling;
#endif

inline ROOT::Internal::InterpreterMutexRegistrationRAII::InterpreterMutexRegistrationRAII(TVirtualMutex* mutex):
   fLockGuard(mutex)
{
   if (gCoreMutex)
      ::gCling->SnapshotMutexState(gCoreMutex);
}
inline ROOT::Internal::InterpreterMutexRegistrationRAII::~InterpreterMutexRegistrationRAII()
{
   if (gCoreMutex)
      ::gCling->ForgetMutexState();
}

#endif
