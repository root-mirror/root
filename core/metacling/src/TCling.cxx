// @(#)root/meta:$Id$
// vim: sw=3 ts=3 expandtab foldmethod=indent

/*************************************************************************
 * Copyright (C) 1995-2012, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TCling

This class defines an interface to the cling C++ interpreter.

Cling is a full ANSI compliant C++-11 interpreter based on
clang/LLVM technology.
*/

#include "TCling.h"

#include "TClingBaseClassInfo.h"
#include "TClingCallFunc.h"
#include "TClingClassInfo.h"
#include "TClingDataMemberInfo.h"
#include "TClingMethodArgInfo.h"
#include "TClingMethodInfo.h"
#include "TClingTypedefInfo.h"
#include "TClingTypeInfo.h"
#include "TClingValue.h"

#include "TROOT.h"
#include "TApplication.h"
#include "TGlobal.h"
#include "TDataType.h"
#include "TClass.h"
#include "TClassEdit.h"
#include "TClassTable.h"
#include "TClingCallbacks.h"
#include "TBaseClass.h"
#include "TDataMember.h"
#include "TMemberInspector.h"
#include "TMethod.h"
#include "TMethodArg.h"
#include "TFunctionTemplate.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TString.h"
#include "THashList.h"
#include "TOrdCollection.h"
#include "TVirtualPad.h"
#include "TSystem.h"
#include "TVirtualMutex.h"
#include "TError.h"
#include "TEnv.h"
#include "TEnum.h"
#include "TEnumConstant.h"
#include "THashTable.h"
#include "RConversionRuleParser.h"
#include "RConfigure.h"
#include "compiledata.h"
#include "TClingUtils.h"
#include "TVirtualCollectionProxy.h"
#include "TVirtualStreamerInfo.h"
#include "TListOfDataMembers.h"
#include "TListOfEnums.h"
#include "TListOfEnumsWithLock.h"
#include "TListOfFunctions.h"
#include "TListOfFunctionTemplates.h"
#include "TProtoClass.h"
#include "TStreamerInfo.h" // This is here to avoid to use the plugin manager
#include "ThreadLocalStorage.h"
#include "TFile.h"
#include "TKey.h"
#include "ClingRAII.h"

#include "clang/AST/ASTContext.h"
#include "clang/AST/Decl.h"
#include "clang/AST/DeclarationName.h"
#include "clang/AST/GlobalDecl.h"
#include "clang/AST/RecordLayout.h"
#include "clang/AST/RecursiveASTVisitor.h"
#include "clang/AST/Type.h"
#include "clang/Basic/SourceLocation.h"
#include "clang/Basic/Specifiers.h"
#include "clang/Basic/TargetInfo.h"
#include "clang/CodeGen/ModuleBuilder.h"
#include "clang/Frontend/CompilerInstance.h"
#include "clang/Frontend/FrontendDiagnostic.h"
#include "clang/Lex/HeaderSearch.h"
#include "clang/Lex/Preprocessor.h"
#include "clang/Lex/PreprocessorOptions.h"
#include "clang/Sema/Lookup.h"
#include "clang/Sema/Sema.h"
#include "clang/Parse/Parser.h"

#include "cling/Interpreter/ClangInternalState.h"
#include "cling/Interpreter/DynamicLibraryManager.h"
#include "cling/Interpreter/Interpreter.h"
#include "cling/Interpreter/LookupHelper.h"
#include "cling/Interpreter/Value.h"
#include "cling/Interpreter/Transaction.h"
#include "cling/MetaProcessor/MetaProcessor.h"
#include "cling/Utils/AST.h"
#include "cling/Utils/ParserStateRAII.h"
#include "cling/Utils/SourceNormalization.h"
#include "cling/Interpreter/Exception.h"

#include "llvm/IR/GlobalValue.h"
#include "llvm/IR/Module.h"

#include "llvm/Support/DynamicLibrary.h"
#include "llvm/Support/raw_ostream.h"
#include "llvm/Support/Path.h"
#include "llvm/Support/Process.h"
#include "llvm/Object/ObjectFile.h"
#include "llvm/Support/FileSystem.h"
#include "llvm/Object/SymbolicFile.h"

#include <algorithm>
#include <iostream>
#include <cassert>
#include <map>
#include <set>
#include <stdexcept>
#include <stdint.h>
#include <fstream>
#include <sstream>
#include <string>
#include <tuple>
#include <typeinfo>
#include <unordered_map>
#include <utility>
#include <vector>
#include <functional>

#ifndef R__WIN32
#include <cxxabi.h>
#define R__DLLEXPORT __attribute__ ((visibility ("default")))
#include <sys/stat.h>
#endif
#include <limits.h>
#include <stdio.h>

#ifdef __APPLE__
#include <dlfcn.h>
#include <mach-o/dyld.h>
#include <mach-o/loader.h>
#endif // __APPLE__

#ifdef R__UNIX
#include <dlfcn.h>
#endif

#if defined(__CYGWIN__)
#include <sys/cygwin.h>
#define HMODULE void *
extern "C" {
   __declspec(dllimport) void * __stdcall GetCurrentProcess();
   __declspec(dllimport) bool __stdcall EnumProcessModules(void *, void **, unsigned long, unsigned long *);
   __declspec(dllimport) unsigned long __stdcall GetModuleFileNameExW(void *, void *, wchar_t *, unsigned long);
}
#endif

// Fragment copied from LLVM's raw_ostream.cpp
#if defined(_MSC_VER)
#ifndef STDIN_FILENO
# define STDIN_FILENO 0
#endif
#ifndef STDOUT_FILENO
# define STDOUT_FILENO 1
#endif
#ifndef STDERR_FILENO
# define STDERR_FILENO 2
#endif
#ifndef R__WIN32
//#if defined(HAVE_UNISTD_H)
# include <unistd.h>
//#endif
#else
#include "Windows4Root.h"
#include <Psapi.h>
#undef GetModuleFileName
#define RTLD_DEFAULT ((void *)::GetModuleHandle(NULL))
#define dlsym(library, function_name) ::GetProcAddress((HMODULE)library, function_name)
#define dlopen(library_name, flags) ::LoadLibraryA(library_name)
#define dlclose(library) ::FreeLibrary((HMODULE)library)
#define R__DLLEXPORT __declspec(dllexport)
#endif
#endif

//______________________________________________________________________________
// Infrastructure to detect and react to libCling being teared down.
//
namespace {
   class TCling_UnloadMarker {
   public:
      ~TCling_UnloadMarker() {
         if (ROOT::Internal::gROOTLocal) {
            ROOT::Internal::gROOTLocal->~TROOT();
         }
      }
   };
   static TCling_UnloadMarker gTClingUnloadMarker;
}



//______________________________________________________________________________
// These functions are helpers for debugging issues with non-LLVMDEV builds.
//
R__DLLEXPORT clang::DeclContext* TCling__DEBUG__getDeclContext(clang::Decl* D) {
   return D->getDeclContext();
}
R__DLLEXPORT clang::NamespaceDecl* TCling__DEBUG__DCtoNamespace(clang::DeclContext* DC) {
   return llvm::dyn_cast<clang::NamespaceDecl>(DC);
}
R__DLLEXPORT clang::RecordDecl* TCling__DEBUG__DCtoRecordDecl(clang::DeclContext* DC) {
   return llvm::dyn_cast<clang::RecordDecl>(DC);
}
R__DLLEXPORT void TCling__DEBUG__dump(clang::DeclContext* DC) {
   return DC->dumpDeclContext();
}
R__DLLEXPORT void TCling__DEBUG__dump(clang::Decl* D) {
   return D->dump();
}
R__DLLEXPORT void TCling__DEBUG__dump(clang::FunctionDecl* FD) {
   return FD->dump();
}
R__DLLEXPORT void TCling__DEBUG__decl_dump(void* D) {
   return ((clang::Decl*)D)->dump();
}
R__DLLEXPORT void TCling__DEBUG__printName(clang::Decl* D) {
   if (clang::NamedDecl* ND = llvm::dyn_cast<clang::NamedDecl>(D)) {
      std::string name;
      {
         llvm::raw_string_ostream OS(name);
         ND->getNameForDiagnostic(OS, D->getASTContext().getPrintingPolicy(),
                                  true /*Qualified*/);
      }
      printf("%s\n", name.c_str());
   }
}
//______________________________________________________________________________
// These functions are helpers for testing issues directly rather than
// relying on side effects.
// This is used for the test for ROOT-7462/ROOT-6070
R__DLLEXPORT bool TCling__TEST_isInvalidDecl(clang::Decl* D) {
   return D->isInvalidDecl();
}
R__DLLEXPORT bool TCling__TEST_isInvalidDecl(ClassInfo_t *input) {
   TClingClassInfo *info( (TClingClassInfo*) input);
   assert(info && info->IsValid());
   return info->GetDecl()->isInvalidDecl();
}

using namespace std;
using namespace clang;
using namespace ROOT;

namespace {
  static const std::string gInterpreterClassDef = R"ICF(
#undef ClassDef
#define ClassDef(name, id) \
_ClassDefInterp_(name,id,virtual,) \
static int DeclFileLine() { return __LINE__; }
#undef ClassDefNV
#define ClassDefNV(name, id) \
_ClassDefInterp_(name,id,,) \
static int DeclFileLine() { return __LINE__; }
#undef ClassDefOverride
#define ClassDefOverride(name, id) \
_ClassDefInterp_(name,id,,override) \
static int DeclFileLine() { return __LINE__; }
)ICF";

  static const std::string gNonInterpreterClassDef = R"ICF(
#define __ROOTCLING__ 1
#undef ClassDef
#define ClassDef(name,id) \
_ClassDefOutline_(name,id,virtual,) \
static int DeclFileLine() { return __LINE__; }
#undef ClassDefNV
#define ClassDefNV(name, id)\
_ClassDefOutline_(name,id,,)\
static int DeclFileLine() { return __LINE__; }
#undef ClassDefOverride
#define ClassDefOverride(name, id)\
_ClassDefOutline_(name,id,,override)\
static int DeclFileLine() { return __LINE__; }
)ICF";

// The macros below use ::Error, so let's ensure it is included
  static const std::string gClassDefInterpMacro = R"ICF(
#include "TError.h"

#define _ClassDefInterp_(name,id,virtual_keyword, overrd) \
private: \
public: \
   static TClass *Class() { static TClass* sIsA = 0; if (!sIsA) sIsA = TClass::GetClass(#name); return sIsA; } \
   static const char *Class_Name() { return #name; } \
   virtual_keyword Bool_t CheckTObjectHashConsistency() const overrd { return true; } \
   static Version_t Class_Version() { return id; } \
   static TClass *Dictionary() { return 0; } \
   virtual_keyword TClass *IsA() const overrd { return name::Class(); } \
   virtual_keyword void ShowMembers(TMemberInspector&insp) const overrd { ::ROOT::Class_ShowMembers(name::Class(), this, insp); } \
   virtual_keyword void Streamer(TBuffer&) overrd { ::Error("Streamer", "Cannot stream interpreted class."); } \
   void StreamerNVirtual(TBuffer&ClassDef_StreamerNVirtual_b) { name::Streamer(ClassDef_StreamerNVirtual_b); } \
   static const char *DeclFileName() { return __FILE__; } \
   static int ImplFileLine() { return 0; } \
   static const char *ImplFileName() { return __FILE__; }
)ICF";
}
R__EXTERN int optind;

// The functions are used to bridge cling/clang/llvm compiled with no-rtti and
// ROOT (which uses rtti)

//______________________________________________________________________________

// Class extracting recursively every Enum type defined for a class.
class EnumVisitor : public RecursiveASTVisitor<EnumVisitor> {
private:
   llvm::SmallVector<EnumDecl*,128> &fClassEnums;
public:
   EnumVisitor(llvm::SmallVector<EnumDecl*,128> &enums) : fClassEnums(enums)
   {}

   bool TraverseStmt(Stmt*) {
      // Don't descend into function bodies.
      return true;
   }

   bool shouldVisitTemplateInstantiations() const { return true; }

   bool TraverseClassTemplateDecl(ClassTemplateDecl*) {
      // Don't descend into templates (but only instances thereof).
      return true; // returning false will abort the in-depth traversal.
   }

   bool TraverseClassTemplatePartialSpecializationDecl(ClassTemplatePartialSpecializationDecl*) {
      // Don't descend into templates partial specialization (but only instances thereof).
      return true; // returning false will abort the in-depth traversal.
   }

   bool VisitEnumDecl(EnumDecl *TEnumD) {
      if (!TEnumD->getDeclContext()->isDependentContext())
         fClassEnums.push_back(TEnumD);
      return true; // returning false will abort the in-depth traversal.
   }
};

////////////////////////////////////////////////////////////////////////////////
/// Print a StackTrace!

extern "C"
void TCling__PrintStackTrace() {
   gSystem->StackTrace();
}

////////////////////////////////////////////////////////////////////////////////
/// Re-apply the lock count delta that TCling__ResetInterpreterMutex() caused.

extern "C" void TCling__RestoreInterpreterMutex(void *delta)
{
   ((TCling*)gCling)->ApplyToInterpreterMutex(delta);
}

////////////////////////////////////////////////////////////////////////////////
/// Lookup libraries in LD_LIBRARY_PATH and DYLD_LIBRARY_PATH with mangled_name,
/// which is extracted by error messages we get from callback from cling. Return true
/// when the missing library was autoloaded.

extern "C" bool TCling__LibraryLoadingFailed(const std::string& errmessage, const std::string& libStem, bool permanent, bool resolved)
{
   return ((TCling*)gCling)->LibraryLoadingFailed(errmessage, libStem, permanent, resolved);
}

////////////////////////////////////////////////////////////////////////////////
/// Reset the interpreter lock to the state it had before interpreter-related
/// calls happened.

extern "C" void *TCling__ResetInterpreterMutex()
{
   return ((TCling*)gCling)->RewindInterpreterMutex();
}

////////////////////////////////////////////////////////////////////////////////
/// Lock the interpreter.

extern "C" void *TCling__LockCompilationDuringUserCodeExecution()
{
   if (gInterpreterMutex) {
      gInterpreterMutex->Lock();
   }
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Unlock the interpreter.

extern "C" void TCling__UnlockCompilationDuringUserCodeExecution(void* /*state*/)
{
   if (gInterpreterMutex) {
      gInterpreterMutex->UnLock();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Update TClingClassInfo for a class (e.g. upon seeing a definition).

static void TCling__UpdateClassInfo(const NamedDecl* TD)
{
   static Bool_t entered = kFALSE;
   static vector<const NamedDecl*> updateList;
   Bool_t topLevel;

   if (entered) topLevel = kFALSE;
   else {
      entered = kTRUE;
      topLevel = kTRUE;
   }
   if (topLevel) {
      ((TCling*)gInterpreter)->UpdateClassInfoWithDecl(TD);
   } else {
      // If we are called indirectly from within another call to
      // TCling::UpdateClassInfo, we delay the update until the dictionary loading
      // is finished (i.e. when we return to the top level TCling::UpdateClassInfo).
      // This allows for the dictionary to be fully populated when we actually
      // update the TClass object.   The updating of the TClass sometimes
      // (STL containers and when there is an emulated class) forces the building
      // of the TClass object's real data (which needs the dictionary info).
      updateList.push_back(TD);
   }
   if (topLevel) {
      while (!updateList.empty()) {
         ((TCling*)gInterpreter)->UpdateClassInfoWithDecl(updateList.back());
         updateList.pop_back();
      }
      entered = kFALSE;
   }
}

void TCling::UpdateEnumConstants(TEnum* enumObj, TClass* cl) const {
   const clang::Decl* D = static_cast<const clang::Decl*>(enumObj->GetDeclId());
   if(const clang::EnumDecl* ED = dyn_cast<clang::EnumDecl>(D)) {
      // Add the constants to the enum type.
      for (EnumDecl::enumerator_iterator EDI = ED->enumerator_begin(),
                EDE = ED->enumerator_end(); EDI != EDE; ++EDI) {
         // Get name of the enum type.
         std::string constbuf;
         if (const NamedDecl* END = llvm::dyn_cast<NamedDecl>(*EDI)) {
            PrintingPolicy Policy((*EDI)->getASTContext().getPrintingPolicy());
            llvm::raw_string_ostream stream(constbuf);
            // Don't trigger fopen of the source file to count lines:
            Policy.AnonymousTagLocations = false;
            (END)->getNameForDiagnostic(stream, Policy, /*Qualified=*/false);
         }
         const char* constantName = constbuf.c_str();

         // Get value of the constant.
         Long64_t value;
         const llvm::APSInt valAPSInt = (*EDI)->getInitVal();
         if (valAPSInt.isSigned()) {
            value = valAPSInt.getSExtValue();
         } else {
            value = valAPSInt.getZExtValue();
         }

         // Create the TEnumConstant or update it if existing
         TEnumConstant* enumConstant = nullptr;
         TClingClassInfo* tcCInfo = (TClingClassInfo*)(cl ? cl->GetClassInfo() : 0);
         TClingDataMemberInfo* tcDmInfo = new TClingDataMemberInfo(fInterpreter, *EDI, tcCInfo);
         DataMemberInfo_t* dmInfo = (DataMemberInfo_t*) tcDmInfo;
         if (TObject* encAsTObj = enumObj->GetConstants()->FindObject(constantName)){
            ((TEnumConstant*)encAsTObj)->Update(dmInfo);
         } else {
            enumConstant = new TEnumConstant(dmInfo, constantName, value, enumObj);
         }

         // Add the global constants to the list of Globals.
         if (!cl) {
            TCollection* globals = gROOT->GetListOfGlobals(false);
            if (!globals->FindObject(constantName)) {
               globals->Add(enumConstant);
            }
         }
      }
   }
}

TEnum* TCling::CreateEnum(void *VD, TClass *cl) const
{
   // Handle new enum declaration for either global and nested enums.

   // Create the enum type.
   TEnum* enumType = 0;
   const clang::Decl* D = static_cast<const clang::Decl*>(VD);
   std::string buf;
   if (const EnumDecl* ED = llvm::dyn_cast<EnumDecl>(D)) {
      // Get name of the enum type.
      PrintingPolicy Policy(ED->getASTContext().getPrintingPolicy());
      llvm::raw_string_ostream stream(buf);
      // Don't trigger fopen of the source file to count lines:
      Policy.AnonymousTagLocations = false;
      ED->getNameForDiagnostic(stream, Policy, /*Qualified=*/false);
      // If the enum is unnamed we do not add it to the list of enums i.e unusable.
   }
   if (buf.empty()) {
      return 0;
   }
   const char* name = buf.c_str();
   enumType = new TEnum(name, VD, cl);
   UpdateEnumConstants(enumType, cl);

   return enumType;
}

void TCling::HandleNewDecl(const void* DV, bool isDeserialized, std::set<TClass*> &modifiedTClasses) {
   // Handle new declaration.
   // Record the modified class, struct and namespaces in 'modifiedTClasses'.

   const clang::Decl* D = static_cast<const clang::Decl*>(DV);

   if (!D->isCanonicalDecl() && !isa<clang::NamespaceDecl>(D)
       && !dyn_cast<clang::RecordDecl>(D)) return;

   if (isa<clang::FunctionDecl>(D->getDeclContext())
       || isa<clang::TagDecl>(D->getDeclContext()))
      return;

   // Don't list templates.
   if (const clang::CXXRecordDecl* RD = dyn_cast<clang::CXXRecordDecl>(D)) {
      if (RD->getDescribedClassTemplate())
         return;
   } else if (const clang::FunctionDecl* FD = dyn_cast<clang::FunctionDecl>(D)) {
      if (FD->getDescribedFunctionTemplate())
         return;
   }

   if (const RecordDecl *TD = dyn_cast<RecordDecl>(D)) {
      if (TD->isCanonicalDecl() || TD->isThisDeclarationADefinition())
         TCling__UpdateClassInfo(TD);
   }
   else if (const NamedDecl *ND = dyn_cast<NamedDecl>(D)) {

      if (const TagDecl *TD = dyn_cast<TagDecl>(D)) {
         // Mostly just for EnumDecl (the other TagDecl are handled
         // by the 'RecordDecl' if statement.
         TCling__UpdateClassInfo(TD);
      } else if (const NamespaceDecl* NSD = dyn_cast<NamespaceDecl>(D)) {
         TCling__UpdateClassInfo(NSD);
      }

      // We care about declarations on the global scope.
      if (!isa<TranslationUnitDecl>(ND->getDeclContext()))
         return;

      // Enums are lazyly created, thus we don not need to handle them here.
      if (isa<EnumDecl>(ND))
         return;

      // ROOT says that global is enum(lazylycreated)/var/field declared on the global
      // scope.
      if (!(isa<VarDecl>(ND)))
         return;

      // Skip if already in the list.
      if (gROOT->GetListOfGlobals()->FindObject(ND->getNameAsString().c_str()))
         return;

      // Put the global constants and global enums in the corresponding lists.
      gROOT->GetListOfGlobals()->Add(new TGlobal((DataMemberInfo_t *)
                                                new TClingDataMemberInfo(fInterpreter,
                                                                        cast<ValueDecl>(ND), 0)));
   }
}

extern "C"
void TCling__GetNormalizedContext(const ROOT::TMetaUtils::TNormalizedCtxt*& normCtxt)
{
   // We are sure in this context of the type of the interpreter
   normCtxt = &( (TCling*) gInterpreter)->GetNormalizedContext();
}

extern "C"
void TCling__UpdateListsOnCommitted(const cling::Transaction &T, cling::Interpreter*) {

   ((TCling*)gCling)->UpdateListsOnCommitted(T);
}

extern "C"
void TCling__UpdateListsOnUnloaded(const cling::Transaction &T) {

   ((TCling*)gCling)->UpdateListsOnUnloaded(T);
}

extern "C"
void TCling__TransactionRollback(const cling::Transaction &T) {

   ((TCling*)gCling)->TransactionRollback(T);
}

extern "C" void TCling__LibraryLoadedRTTI(const void* dyLibHandle,
                                          const char* canonicalName) {

   ((TCling*)gCling)->LibraryLoaded(dyLibHandle, canonicalName);
}

extern "C" void TCling__LibraryUnloadedRTTI(const void* dyLibHandle,
                                            const char* canonicalName) {

   ((TCling*)gCling)->LibraryUnloaded(dyLibHandle, canonicalName);
}


extern "C"
TObject* TCling__GetObjectAddress(const char *Name, void *&LookupCtx) {
   return ((TCling*)gCling)->GetObjectAddress(Name, LookupCtx);
}

extern "C" const Decl* TCling__GetObjectDecl(TObject *obj) {
   return ((TClingClassInfo*)obj->IsA()->GetClassInfo())->GetDecl();
}

extern "C" R__DLLEXPORT TInterpreter *CreateInterpreter(void* interpLibHandle,
                                                        const char* argv[])
{
   cling::DynamicLibraryManager::ExposeHiddenSharedLibrarySymbols(interpLibHandle);
   return new TCling("C++", "cling C++ Interpreter", argv);
}

extern "C" R__DLLEXPORT void DestroyInterpreter(TInterpreter *interp)
{
   delete interp;
}

// Load library containing specified class. Returns 0 in case of error
// and 1 in case if success.
extern "C" int TCling__AutoLoadCallback(const char* className)
{
   return ((TCling*)gCling)->AutoLoad(className);
}

extern "C" int TCling__AutoParseCallback(const char* className)
{
   return ((TCling*)gCling)->AutoParse(className);
}

extern "C" const char* TCling__GetClassSharedLibs(const char* className)
{
   return ((TCling*)gCling)->GetClassSharedLibs(className);
}

// // Returns 0 for failure 1 for success
// extern "C" int TCling__IsAutoLoadNamespaceCandidate(const char* name)
// {
//    return ((TCling*)gCling)->IsAutoLoadNamespaceCandidate(name);
// }

// Returns 0 for failure 1 for success
extern "C" int TCling__IsAutoLoadNamespaceCandidate(const clang::NamespaceDecl* nsDecl)
{
   return ((TCling*)gCling)->IsAutoLoadNamespaceCandidate(nsDecl);
}

extern "C" int TCling__CompileMacro(const char *fileName, const char *options)
{
   string file(fileName);
   string opt(options);
   return gSystem->CompileMacro(file.c_str(), opt.c_str());
}

extern "C" void TCling__SplitAclicMode(const char* fileName, string &mode,
                                       string &args, string &io, string &fname)
{
   string file(fileName);
   TString f, amode, arguments, aclicio;
   f = gSystem->SplitAclicMode(file.c_str(), amode, arguments, aclicio);
   mode = amode.Data(); args = arguments.Data();
   io = aclicio.Data(); fname = f.Data();
}

// Implemented in TClingCallbacks.
extern "C" void TCling__FindLoadedLibraries(std::vector<std::pair<uint32_t, std::string>> &sLibraries,
                                 std::vector<std::string> &sPaths,
                                 cling::Interpreter &interpreter, bool searchSystem);

//______________________________________________________________________________
//
//
//

#ifdef R__WIN32
extern "C" {
   char *__unDName(char *demangled, const char *mangled, int out_len,
                   void * (* pAlloc )(size_t), void (* pFree )(void *),
                   unsigned short int flags);
}
#endif

////////////////////////////////////////////////////////////////////////////////
/// Find a template decl within N nested namespaces, 0<=N<inf
/// Assumes 1 and only 1 template present and 1 and only 1 entity contained
/// by the namespace. Example: ns1::ns2::..::nsN::myTemplate
/// Returns nullptr in case of error

static clang::ClassTemplateDecl* FindTemplateInNamespace(clang::Decl* decl)
{
   using namespace clang;
   if (NamespaceDecl* nsd = llvm::dyn_cast<NamespaceDecl>(decl)){
      return FindTemplateInNamespace(*nsd->decls_begin());
   }

   if (ClassTemplateDecl* ctd = llvm::dyn_cast<ClassTemplateDecl>(decl)){
      return ctd;
   }

   return nullptr; // something went wrong.
}

////////////////////////////////////////////////////////////////////////////////
/// Autoload a library provided the mangled name of a missing symbol.

void* llvmLazyFunctionCreator(const std::string& mangled_name)
{
   return ((TCling*)gCling)->LazyFunctionCreatorAutoload(mangled_name);
}

//______________________________________________________________________________
//
//
//

int TCling_GenerateDictionary(const std::vector<std::string> &classes,
                              const std::vector<std::string> &headers,
                              const std::vector<std::string> &fwdDecls,
                              const std::vector<std::string> &unknown)
{
   //This function automatically creates the "LinkDef.h" file for templated
   //classes then executes CompileMacro on it.
   //The name of the file depends on the class name, and it's not generated again
   //if the file exist.
   if (classes.empty()) {
      return 0;
   }
   // Use the name of the first class as the main name.
   const std::string& className = classes[0];
   //(0) prepare file name
   TString fileName = "AutoDict_";
   std::string::const_iterator sIt;
   for (sIt = className.begin(); sIt != className.end(); ++sIt) {
      if (*sIt == '<' || *sIt == '>' ||
            *sIt == ' ' || *sIt == '*' ||
            *sIt == ',' || *sIt == '&' ||
            *sIt == ':') {
         fileName += '_';
      }
      else {
         fileName += *sIt;
      }
   }
   if (classes.size() > 1) {
      Int_t chk = 0;
      std::vector<std::string>::const_iterator it = classes.begin();
      while ((++it) != classes.end()) {
         for (UInt_t cursor = 0; cursor != it->length(); ++cursor) {
            chk = chk * 3 + it->at(cursor);
         }
      }
      fileName += TString::Format("_%u", chk);
   }
   fileName += ".cxx";
   if (gSystem->AccessPathName(fileName) != 0) {
      //file does not exist
      //(1) prepare file data
      // If STL, also request iterators' operators.
      // vector is special: we need to check whether
      // vector::iterator is a typedef to pointer or a
      // class.
      static const std::set<std::string> sSTLTypes {
         "vector","list","forward_list","deque","map","unordered_map","multimap",
         "unordered_multimap","set","unordered_set","multiset","unordered_multiset",
         "queue","priority_queue","stack","iterator"};
      std::vector<std::string>::const_iterator it;
      std::string fileContent("");
      for (it = headers.begin(); it != headers.end(); ++it) {
         fileContent += "#include \"" + *it + "\"\n";
      }
      for (it = unknown.begin(); it != unknown.end(); ++it) {
         TClass* cl = TClass::GetClass(it->c_str());
         if (cl && cl->GetDeclFileName()) {
            TString header(gSystem->BaseName(cl->GetDeclFileName()));
            TString dir(gSystem->DirName(cl->GetDeclFileName()));
            TString dirbase(gSystem->BaseName(dir));
            while (dirbase.Length() && dirbase != "."
                   && dirbase != "include" && dirbase != "inc"
                   && dirbase != "prec_stl") {
               gSystem->PrependPathName(dirbase, header);
               dir = gSystem->DirName(dir);
            }
            fileContent += TString("#include \"") + header + "\"\n";
         }
      }
      for (it = fwdDecls.begin(); it != fwdDecls.end(); ++it) {
         fileContent += "class " + *it + ";\n";
      }
      fileContent += "#ifdef __CINT__ \n";
      fileContent += "#pragma link C++ nestedclasses;\n";
      fileContent += "#pragma link C++ nestedtypedefs;\n";
      for (it = classes.begin(); it != classes.end(); ++it) {
         std::string n(*it);
         size_t posTemplate = n.find('<');
         std::set<std::string>::const_iterator iSTLType = sSTLTypes.end();
         if (posTemplate != std::string::npos) {
            n.erase(posTemplate, std::string::npos);
            if (n.compare(0, 5, "std::") == 0) {
               n.erase(0, 5);
            }
            iSTLType = sSTLTypes.find(n);
         }
         fileContent += "#pragma link C++ class ";
         fileContent +=    *it + "+;\n" ;
         fileContent += "#pragma link C++ class ";
         if (iSTLType != sSTLTypes.end()) {
            // STL class; we cannot (and don't need to) store iterators;
            // their shadow and the compiler's version don't agree. So
            // don't ask for the '+'
            fileContent +=    *it + "::*;\n" ;
         }
         else {
            // Not an STL class; we need to allow the I/O of contained
            // classes (now that we have a dictionary for them).
            fileContent +=    *it + "::*+;\n" ;
         }
      }
      fileContent += "#endif\n";
      //end(1)
      //(2) prepare the file
      FILE* filePointer;
      filePointer = fopen(fileName, "w");
      if (filePointer == NULL) {
         //can't open a file
         return 1;
      }
      //end(2)
      //write data into the file
      fprintf(filePointer, "%s", fileContent.c_str());
      fclose(filePointer);
   }
   //(3) checking if we can compile a macro, if not then cleaning
   Int_t oldErrorIgnoreLevel = gErrorIgnoreLevel;
   gErrorIgnoreLevel = kWarning; // no "Info: creating library..."
   Int_t ret = gSystem->CompileMacro(fileName, "k");
   gErrorIgnoreLevel = oldErrorIgnoreLevel;
   if (ret == 0) { //can't compile a macro
      return 2;
   }
   //end(3)
   return 0;
}

int TCling_GenerateDictionary(const std::string& className,
                              const std::vector<std::string> &headers,
                              const std::vector<std::string> &fwdDecls,
                              const std::vector<std::string> &unknown)
{
   //This function automatically creates the "LinkDef.h" file for templated
   //classes then executes CompileMacro on it.
   //The name of the file depends on the class name, and it's not generated again
   //if the file exist.
   std::vector<std::string> classes;
   classes.push_back(className);
   return TCling_GenerateDictionary(classes, headers, fwdDecls, unknown);
}

//______________________________________________________________________________
//
//
//

// It is a "fantom" method to synchronize user keyboard input
// and ROOT prompt line (for WIN32)
const char* fantomline = "TRint::EndOfLineAction();";

//______________________________________________________________________________
//
//
//

void* TCling::fgSetOfSpecials = 0;

//______________________________________________________________________________
//
// llvm error handler through exceptions; see also cling/UserInterface
//
namespace {
   // Handle fatal llvm errors by throwing an exception.
   // Yes, throwing exceptions in error handlers is bad.
   // Doing nothing is pretty terrible, too.
   void exceptionErrorHandler(void * /*user_data*/,
                              const std::string& reason,
                              bool /*gen_crash_diag*/) {
      throw std::runtime_error(std::string(">>> Interpreter compilation error:\n") + reason);
   }
}

//______________________________________________________________________________
//
//
//

////////////////////////////////////////////////////////////////////////////////

namespace{
   // An instance of this class causes the diagnostics of clang to be suppressed
   // during its lifetime
   class clangDiagSuppr {
   public:
      clangDiagSuppr(clang::DiagnosticsEngine& diag): fDiagEngine(diag){
         fOldDiagValue = fDiagEngine.getIgnoreAllWarnings();
         fDiagEngine.setIgnoreAllWarnings(true);
      }

      ~clangDiagSuppr() {
         fDiagEngine.setIgnoreAllWarnings(fOldDiagValue);
      }
   private:
      clang::DiagnosticsEngine& fDiagEngine;
      bool fOldDiagValue;
   };

}

////////////////////////////////////////////////////////////////////////////////
/// Allow calling autoparsing from TMetaUtils
bool TClingLookupHelper__AutoParse(const char *cname)
{
   return gCling->AutoParse(cname);
}

////////////////////////////////////////////////////////////////////////////////
/// Try hard to avoid looking up in the Cling database as this could enduce
/// an unwanted autoparsing.

bool TClingLookupHelper__ExistingTypeCheck(const std::string &tname,
                                           std::string &result)
{
   result.clear();

   unsigned long offset = 0;
   if (strncmp(tname.c_str(), "const ", 6) == 0) {
      offset = 6;
   }
   unsigned long end = tname.length();
   while( end && (tname[end-1]=='&' || tname[end-1]=='*' || tname[end-1]==']') ) {
      if ( tname[end-1]==']' ) {
         --end;
         while ( end && tname[end-1]!='[' ) --end;
      }
      --end;
   }
   std::string innerbuf;
   const char *inner;
   if (end != tname.length()) {
      innerbuf = tname.substr(offset,end-offset);
      inner = innerbuf.c_str();
   } else {
      inner = tname.c_str()+offset;
   }

   //if (strchr(tname.c_str(),'[')!=0) fprintf(stderr,"DEBUG: checking on %s vs %s %lu %lu\n",tname.c_str(),inner,offset,end);
   if (gROOT->GetListOfClasses()->FindObject(inner)
       || TClassTable::Check(inner,result) ) {
      // This is a known class.
      return true;
   }

   THashTable *typeTable = dynamic_cast<THashTable*>( gROOT->GetListOfTypes() );
   TDataType *type = (TDataType *)typeTable->THashTable::FindObject( inner );
   if (type) {
      // This is a raw type and an already loaded typedef.
      const char *newname = type->GetFullTypeName();
      if (type->GetType() == kLong64_t) {
         newname = "Long64_t";
      } else if (type->GetType() == kULong64_t) {
         newname = "ULong64_t";
      }
      if (strcmp(inner,newname) == 0) {
         return true;
      }
      if (offset) result = "const ";
      result += newname;
      if ( end != tname.length() ) {
         result += tname.substr(end,tname.length()-end);
      }
      if (result == tname) result.clear();
      return true;
   }

   // Check if the name is an enumerator
   const auto lastPos = TClassEdit::GetUnqualifiedName(inner);
   if (lastPos != inner)   // Main switch: case 1 - scoped enum, case 2 global enum
   {
      // We have a scope
      // All of this C gymnastic is to avoid allocations on the heap
      const auto enName = lastPos;
      const auto scopeNameSize = ((Long64_t)lastPos - (Long64_t)inner) / sizeof(decltype(*lastPos)) - 2;
      char *scopeName = new char[scopeNameSize + 1];
      strncpy(scopeName, inner, scopeNameSize);
      scopeName[scopeNameSize] = '\0';
      // Check if the scope is in the list of classes
      if (auto scope = static_cast<TClass *>(gROOT->GetListOfClasses()->FindObject(scopeName))) {
         auto enumTable = dynamic_cast<const THashList *>(scope->GetListOfEnums(false));
         if (enumTable && enumTable->THashList::FindObject(enName)) return true;
      }
      // It may still be in one of the loaded protoclasses
      else if (auto scope = static_cast<TProtoClass *>(gClassTable->GetProtoNorm(scopeName))) {
         auto listOfEnums = scope->GetListOfEnums();
         if (listOfEnums) { // it could be null: no enumerators in the protoclass
            auto enumTable = dynamic_cast<const THashList *>(listOfEnums);
            if (enumTable && enumTable->THashList::FindObject(enName)) return true;
         }
      }
      delete [] scopeName;
   } else
   {
      // We don't have any scope: this could only be a global enum
      auto enumTable = dynamic_cast<const THashList *>(gROOT->GetListOfEnums());
      if (enumTable && enumTable->THashList::FindObject(inner)) return true;
   }

   if (gCling->GetClassSharedLibs(inner))
   {
      // This is a class name.
      return true;
   }

   return false;
}

////////////////////////////////////////////////////////////////////////////////

TCling::TUniqueString::TUniqueString(Long64_t size)
{
   fContent.reserve(size);
}

////////////////////////////////////////////////////////////////////////////////

inline const char *TCling::TUniqueString::Data()
{
   return fContent.c_str();
}

////////////////////////////////////////////////////////////////////////////////
/// Append string to the storage if not added already.

inline bool TCling::TUniqueString::Append(const std::string& str)
{
   bool notPresent = fLinesHashSet.emplace(fHashFunc(str)).second;
   if (notPresent){
      fContent+=str;
   }
   return notPresent;
}

std::string TCling::ToString(const char* type, void* obj)
{
   return fInterpreter->toString(type, obj);
}

////////////////////////////////////////////////////////////////////////////////
///\returns true if the module was loaded.
static bool LoadModule(const std::string &ModuleName, cling::Interpreter &interp, bool Complain = true)
{
   if (interp.loadModule(ModuleName, Complain))
      return true;

   // When starting up ROOT, cling would load all modulemap files on the include
   // paths. However, in a ROOT session, it is very common to run aclic which
   // will invoke rootcling and possibly produce a modulemap and a module in
   // the current folder.
   //
   // Before failing, try loading the modulemap in the current folder and try
   // loading the requested module from it.
   Preprocessor &PP = interp.getCI()->getPreprocessor();
   FileManager& FM = PP.getFileManager();
   // FIXME: In a ROOT session we can add an include path (through .I /inc/path)
   // We should look for modulemap files there too.
   const DirectoryEntry *DE = FM.getDirectory(".");
   if (DE) {
      HeaderSearch& HS = PP.getHeaderSearchInfo();
      const FileEntry *FE = HS.lookupModuleMapFile(DE, /*IsFramework*/ false);
      if (FE && !gCling->IsLoaded("./module.modulemap")) {
         if (!HS.loadModuleMapFile(FE, /*IsSystem*/ false))
            return LoadModule(ModuleName, interp, Complain);
         Error("TCling::LoadModule", "Could not load modulemap in the current directory");
      }
   }

   if (Complain)
      Error("TCling::LoadModule", "Module %s not found!", ModuleName.c_str());

   return false;
}

////////////////////////////////////////////////////////////////////////////////
/// Loads the C++ modules that we require to run any ROOT program. This is just
/// supposed to make a C++ module from a modulemap available to the interpreter.
static void LoadModules(const std::vector<std::string> &modules, cling::Interpreter &interp)
{
   for (const auto &modName : modules)
      LoadModule(modName, interp);
}

static bool IsFromRootCling() {
  // rootcling also uses TCling for generating the dictionary ROOT files.
  const static bool foundSymbol = dlsym(RTLD_DEFAULT, "usedToIdentifyRootClingByDlSym");
  return foundSymbol;
}

static std::string GetModuleNameAsString(clang::Module *M, const clang::Preprocessor &PP)
{
   const HeaderSearchOptions &HSOpts = PP.getHeaderSearchInfo().getHeaderSearchOpts();

   std::string ModuleFileName;
   if (!HSOpts.PrebuiltModulePaths.empty())
      // Load the module from *only* in the prebuilt module path.
      ModuleFileName = PP.getHeaderSearchInfo().getModuleFileName(M->Name, /*ModuleMapPath*/"", /*UsePrebuiltPath*/ true);
   if (ModuleFileName.empty()) return "";

   std::string ModuleName = llvm::sys::path::filename(ModuleFileName);
   // Return stem of the filename
   return std::string(llvm::sys::path::stem(ModuleName));
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize the cling interpreter interface.
/// \param argv - array of arguments passed to the cling::Interpreter constructor
///               e.g. `-DFOO=bar`. The last element of the array must be `nullptr`.

TCling::TCling(const char *name, const char *title, const char* const argv[])
: TInterpreter(name, title), fGlobalsListSerial(-1), fInterpreter(0),
   fMetaProcessor(0), fNormalizedCtxt(0), fPrevLoadedDynLibInfo(0),
   fClingCallbacks(0), fAutoLoadCallBack(0),
   fTransactionCount(0), fHeaderParsingOnDemand(true), fIsAutoParsingSuspended(kFALSE)
{
   const bool fromRootCling = IsFromRootCling();

   fCxxModulesEnabled = false;
#ifdef R__USE_CXXMODULES
   fCxxModulesEnabled = true;
#endif

   llvm::install_fatal_error_handler(&exceptionErrorHandler);

   fTemporaries = new std::vector<cling::Value>();

   std::vector<std::string> clingArgsStorage;
   clingArgsStorage.push_back("cling4root");
   for (const char* const* arg = argv; *arg; ++arg)
      clingArgsStorage.push_back(*arg);

   // rootcling sets its arguments through TROOT::GetExtraInterpreterArgs().
   if (!fromRootCling) {
      ROOT::TMetaUtils::SetPathsForRelocatability(clingArgsStorage);

      // Add -I early so ASTReader can find the headers.
      std::string interpInclude(TROOT::GetEtcDir().Data());
      clingArgsStorage.push_back("-I" + interpInclude);

      // Add include path to etc/cling.
      clingArgsStorage.push_back("-I" + interpInclude + "/cling");

      // Add the root include directory and etc/ to list searched by default.
      clingArgsStorage.push_back(std::string(("-I" + TROOT::GetIncludeDir()).Data()));

      // Add the current path to the include path
      // TCling::AddIncludePath(".");

      // Attach the PCH (unless we have C++ modules enabled which provide the
      // same functionality).
      if (!fCxxModulesEnabled) {
         std::string pchFilename = interpInclude + "/allDict.cxx.pch";
         if (gSystem->Getenv("ROOT_PCH")) {
            pchFilename = gSystem->Getenv("ROOT_PCH");
         }

         clingArgsStorage.push_back("-include-pch");
         clingArgsStorage.push_back(pchFilename);
      }

      clingArgsStorage.push_back("-Wno-undefined-inline");
      clingArgsStorage.push_back("-fsigned-char");
   }

   // Process externally passed arguments if present.
   llvm::Optional<std::string> EnvOpt = llvm::sys::Process::GetEnv("EXTRA_CLING_ARGS");
   if (EnvOpt.hasValue()) {
      StringRef Env(*EnvOpt);
      while (!Env.empty()) {
         StringRef Arg;
         std::tie(Arg, Env) = Env.split(' ');
         clingArgsStorage.push_back(Arg.str());
      }
   }

   if (fCxxModulesEnabled) {
      clingArgsStorage.push_back("-modulemap_overlay=" + std::string(TROOT::GetIncludeDir().Data()));
   }

   // FIXME: This only will enable frontend timing reports.
   EnvOpt = llvm::sys::Process::GetEnv("ROOT_CLING_TIMING");
   if (EnvOpt.hasValue())
     clingArgsStorage.push_back("-ftime-report");

   std::vector<const char*> interpArgs;
   for (std::vector<std::string>::const_iterator iArg = clingArgsStorage.begin(),
           eArg = clingArgsStorage.end(); iArg != eArg; ++iArg)
      interpArgs.push_back(iArg->c_str());

   // Activate C++ modules support. If we are running within rootcling, it's up
   // to rootcling to set this flag depending on whether it wants to produce
   // C++ modules.
   TString vfsArg;
   if (fCxxModulesEnabled && !fromRootCling) {
      // We only set this flag, rest is done by the CIFactory.
      interpArgs.push_back("-fmodules");
      // We should never build modules during runtime, so let's enable the
      // module build remarks from clang to make it easier to spot when we do
      // this by accident.
      interpArgs.push_back("-Rmodule-build");
   }

#ifdef R__FAST_MATH
   // Same setting as in rootcling_impl.cxx.
   interpArgs.push_back("-ffast-math");
#endif

#ifdef R__EXTERN_LLVMDIR
   TString llvmResourceDir = R__EXTERN_LLVMDIR;
#else
   TString llvmResourceDir = TROOT::GetEtcDir() + "/cling";
#endif
   // Add statically injected extra arguments, usually coming from rootcling.
   for (const char** extraArgs = TROOT::GetExtraInterpreterArgs();
        extraArgs && *extraArgs; ++extraArgs) {
      if (!strcmp(*extraArgs, "-resource-dir")) {
         // Take the next arg as the llvm resource directory.
         llvmResourceDir = *(++extraArgs);
      } else {
         interpArgs.push_back(*extraArgs);
      }
   }

   for (const auto &arg: TROOT::AddExtraInterpreterArgs({})) {
      interpArgs.push_back(arg.c_str());
   }

   fInterpreter = new cling::Interpreter(interpArgs.size(),
                                         &(interpArgs[0]),
                                         llvmResourceDir);

   if (!fromRootCling) {
      fInterpreter->installLazyFunctionCreator(llvmLazyFunctionCreator);
   }

   // Don't check whether modules' files exist.
   fInterpreter->getCI()->getPreprocessorOpts().DisablePCHValidation = true;

   // Until we can disable autoloading during Sema::CorrectTypo() we have
   // to disable spell checking.
   fInterpreter->getCI()->getLangOpts().SpellChecking = false;


   // We need stream that doesn't close its file descriptor, thus we are not
   // using llvm::outs. Keeping file descriptor open we will be able to use
   // the results in pipes (Savannah #99234).
   static llvm::raw_fd_ostream fMPOuts (STDOUT_FILENO, /*ShouldClose*/false);
   fMetaProcessor = new cling::MetaProcessor(*fInterpreter, fMPOuts);

   if (fInterpreter->getCI()->getLangOpts().Modules) {
      // Setup core C++ modules if we have any to setup.

      // Load libc and stl first.
#ifdef R__MACOSX
      LoadModules({"Darwin", "std"}, *fInterpreter);
#else
      LoadModules({"libc", "stl"}, *fInterpreter);
#endif

      // Load core modules
      // This should be vector in order to be able to pass it to LoadModules
      std::vector<std::string> CoreModules = {"ROOT_Foundation_C","ROOT_Config",
         "ROOT_Foundation_Stage1_NoRTTI", "Core", "RIO"};
      // These modules contain global variables which conflict with users' code such as "PI".
      // FIXME: Reducing those will let us be less dependent on rootmap files
      static constexpr std::array<const char*, 3> ExcludeModules =
         { { "Rtools", "RSQLite", "RInterface"} }; 

      LoadModules(CoreModules, *fInterpreter);

      // Take this branch only from ROOT because we don't need to preload modules in rootcling
      if (!fromRootCling) {
         // Dynamically get all the modules and load them if they are not in core modules
         clang::CompilerInstance &CI = *fInterpreter->getCI();
         clang::ModuleMap &moduleMap = CI.getPreprocessor().getHeaderSearchInfo().getModuleMap();
         clang::Preprocessor &PP = CI.getPreprocessor();
         std::vector<std::string> ModulesPreloaded;

         for (auto I = moduleMap.module_begin(), E = moduleMap.module_end(); I != E; ++I) {
            clang::Module *M = I->second;
            assert(M);

            std::string ModuleName = GetModuleNameAsString(M, PP);
            if (!ModuleName.empty() &&
                  std::find(CoreModules.begin(), CoreModules.end(), ModuleName) == CoreModules.end()
                  && std::find(ExcludeModules.begin(), ExcludeModules.end(), ModuleName) == ExcludeModules.end()) {
               if (M->IsSystem && !M->IsMissingRequirement)
                  LoadModule(ModuleName, *fInterpreter);
               else if (!M->IsSystem && !M->IsMissingRequirement)
                  ModulesPreloaded.push_back(ModuleName);
            }
         }
         LoadModules(ModulesPreloaded, *fInterpreter);
      }

      // Check that the gROOT macro was exported by any core module.
      assert(fInterpreter->getMacro("gROOT") && "Couldn't load gROOT macro?");

      // C99 decided that it's a very good idea to name a macro `I` (the letter I).
      // This seems to screw up nearly all the template code out there as `I` is
      // common template parameter name and iterator variable name.
      // Let's follow the GCC recommendation and undefine `I` in case any of the
      // core modules have defined it:
      // https://www.gnu.org/software/libc/manual/html_node/Complex-Numbers.html
      fInterpreter->declare("#ifdef I\n #undef I\n #endif\n");
   }

   // For the list to also include string, we have to include it now.
   // rootcling does parts already if needed, e.g. genreflex does not want using
   // namespace std.
   if (fromRootCling) {
      fInterpreter->declare("#include \"RtypesCore.h\"\n"
                            "#include <string>\n"
                            "using std::string;\n"
                            "#include <cassert>\n");
   } else {
      fInterpreter->declare("#include \"Rtypes.h\"\n"
                            + gClassDefInterpMacro + "\n"
                            + gInterpreterClassDef + "\n"
                            "#undef ClassImp\n"
                            "#define ClassImp(X);\n"
                            "#include <string>\n"
                            "using namespace std;\n"
                            "#include <cassert>\n");
   }

   // We are now ready (enough is loaded) to init the list of opaque typedefs.
   fNormalizedCtxt = new ROOT::TMetaUtils::TNormalizedCtxt(fInterpreter->getLookupHelper());
   fLookupHelper = new ROOT::TMetaUtils::TClingLookupHelper(*fInterpreter, *fNormalizedCtxt, TClingLookupHelper__ExistingTypeCheck, TClingLookupHelper__AutoParse);
   TClassEdit::Init(fLookupHelper);

   // Initialize the cling interpreter interface.
   fMore      = 0;
   fPrompt[0] = 0;
   fMapfile   = 0;
//    fMapNamespaces   = 0;
   fRootmapFiles = 0;
   fLockProcessLine = kTRUE;

   fAllowLibLoad = !fromRootCling;
   // Disallow auto-parsing in rootcling
   fIsAutoParsingSuspended = fromRootCling;
   // Disable the autoloader until it is explicitly enabled.
   SetClassAutoloading(false);

   ResetAll();
#ifndef R__WIN32
   //optind = 1;  // make sure getopt() works in the main program
#endif // R__WIN32

   // Enable dynamic lookup
   if (!fromRootCling) {
      fInterpreter->enableDynamicLookup();
   }

   // Attach cling callbacks last; they might need TROOT::fInterpreter
   // and should thus not be triggered during the equivalent of
   // TROOT::fInterpreter = new TCling;
   std::unique_ptr<TClingCallbacks>
      clingCallbacks(new TClingCallbacks(fInterpreter));
   fClingCallbacks = clingCallbacks.get();
   fClingCallbacks->SetAutoParsingSuspended(fIsAutoParsingSuspended);
   fInterpreter->setCallbacks(std::move(clingCallbacks));

}


////////////////////////////////////////////////////////////////////////////////
/// Destroy the interpreter interface.

TCling::~TCling()
{
   delete fMapfile;
//    delete fMapNamespaces;
   delete fRootmapFiles;
   delete fMetaProcessor;
   delete fTemporaries;
   delete fNormalizedCtxt;
   delete fInterpreter;
   delete fLookupHelper;
   gCling = 0;
#if defined(R__MUST_REVISIT)
#if R__MUST_REVISIT(6,2)
   Warning("~TCling", "Interface not available yet.");
#ifdef R__COMPLETE_MEM_TERMINATION
   // remove all cling objects
#endif
#endif
#endif
   //--
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize the interpreter, once TROOT::fInterpreter is set.

void TCling::Initialize()
{
   fClingCallbacks->Initialize();
}

////////////////////////////////////////////////////////////////////////////////
/// Wrapper around dladdr (and friends)

static const char *FindLibraryName(void (*func)())
{
#if defined(__CYGWIN__) && defined(__GNUC__)
   return 0;
#elif defined(G__WIN32)
   MEMORY_BASIC_INFORMATION mbi;
   if (!VirtualQuery (func, &mbi, sizeof (mbi)))
   {
      return 0;
   }

   HMODULE hMod = (HMODULE) mbi.AllocationBase;
   TTHREAD_TLS_ARRAY(char, MAX_PATH, moduleName);

   if (!GetModuleFileNameA (hMod, moduleName, sizeof (moduleName)))
   {
      return 0;
   }
   return moduleName;
#else
   Dl_info info;
   if (dladdr((void*)func,&info)==0) {
      // Not in a known share library, let's give up
      return 0;
   } else {
      //fprintf(stdout,"Found address in %s\n",info.dli_fname);
      return info.dli_fname;
   }
#endif

}

////////////////////////////////////////////////////////////////////////////////
/// Helper to initialize TVirtualStreamerInfo's factor early.
/// Use static initialization to insure only one TStreamerInfo is created.
static bool R__InitStreamerInfoFactory()
{
   // Use lambda since SetFactory return void.
   auto setFactory = []() {
      TVirtualStreamerInfo::SetFactory(new TStreamerInfo());
      return kTRUE;
   };
   static bool doneFactory = setFactory();
   return doneFactory; // avoid unused variable warning.
}


////////////////////////////////////////////////////////////////////////////////
/// Tries to load a PCM; returns true on success.

bool TCling::LoadPCM(TString pcmFileName,
                     const char** headers,
                     void (*triggerFunc)()) const {
   // pcmFileName is an intentional copy; updated by FindFile() below.

   TString searchPath;

   if (triggerFunc) {
      const char *libraryName = FindLibraryName(triggerFunc);
      if (libraryName) {
         searchPath = llvm::sys::path::parent_path(libraryName);
#ifdef R__WIN32
         searchPath += ";";
#else
         searchPath += ":";
#endif
      }
   }
   // Note: if we know where the library is, we probably shouldn't even
   // look in other places.
   searchPath.Append( gSystem->GetDynamicPath() );

   if (!gSystem->FindFile(searchPath, pcmFileName))
      return kFALSE;

   // Prevent the ROOT-PCMs hitting this during auto-load during
   // JITting - which will cause recursive compilation.
   // Avoid to call the plugin manager at all.
   R__InitStreamerInfoFactory();

   if (gROOT->IsRootFile(pcmFileName)) {
      Int_t oldDebug = gDebug;
      if (gDebug > 5) {
         gDebug -= 5;
         ::Info("TCling::LoadPCM", "Loading ROOT PCM %s", pcmFileName.Data());
      } else {
         gDebug = 0;
      }

      TDirectory::TContext ctxt;

      TFile *pcmFile = new TFile(pcmFileName+"?filetype=pcm","READ");

      auto listOfKeys = pcmFile->GetListOfKeys();

      // This is an empty pcm
      if (
         listOfKeys &&
         (
            (listOfKeys->GetSize() == 0) || // Nothing here, or
            (
               (listOfKeys->GetSize() == 1) && // only one, and
               !strcmp(((TKey*)listOfKeys->At(0))->GetName(), "EMPTY") // name is EMPTY
            )
         )
      ) {
         delete pcmFile;
         gDebug = oldDebug;
         return kTRUE;
      }

      TObjArray *protoClasses;
      if (gDebug > 1)
            ::Info("TCling::LoadPCM","reading protoclasses for %s \n",pcmFileName.Data());

      pcmFile->GetObject("__ProtoClasses", protoClasses);

      if (protoClasses) {
         for (auto obj : *protoClasses) {
            TProtoClass * proto = (TProtoClass*)obj;
            TClassTable::Add(proto);
         }
         // Now that all TClass-es know how to set them up we can update
         // existing TClasses, which might cause the creation of e.g. TBaseClass
         // objects which in turn requires the creation of TClasses, that could
         // come from the PCH, but maybe later in the loop. Instead of resolving
         // a dependency graph the addition to the TClassTable above allows us
         // to create these dependent TClasses as needed below.
         for (auto proto : *protoClasses) {
            if (TClass* existingCl
                = (TClass*)gROOT->GetListOfClasses()->FindObject(proto->GetName())) {
               // We have an existing TClass object. It might be emulated
               // or interpreted; we now have more information available.
               // Make that available.
               if (existingCl->GetState() != TClass::kHasTClassInit) {
                  DictFuncPtr_t dict = gClassTable->GetDict(proto->GetName());
                  if (!dict) {
                     ::Error("TCling::LoadPCM", "Inconsistent TClassTable for %s",
                             proto->GetName());
                  } else {
                     // This will replace the existing TClass.
                     TClass *ncl = (*dict)();
                     if (ncl) ncl->PostLoadCheck();

                  }
               }
            }
         }

         protoClasses->Clear(); // Ownership was transfered to TClassTable.
         delete protoClasses;
      }

      TObjArray *dataTypes;
      pcmFile->GetObject("__Typedefs", dataTypes);
      if (dataTypes) {
         for (auto typedf: *dataTypes)
            gROOT->GetListOfTypes()->Add(typedf);
         dataTypes->Clear(); // Ownership was transfered to TListOfTypes.
         delete dataTypes;
      }

      TObjArray *enums;
      pcmFile->GetObject("__Enums", enums);
      if (enums) {
         // Cache the pointers
         auto listOfGlobals = gROOT->GetListOfGlobals();
         auto listOfEnums = dynamic_cast<THashList*>(gROOT->GetListOfEnums());
         // Loop on enums and then on enum constants
         for (auto selEnum: *enums){
            const char* enumScope = selEnum->GetTitle();
            const char* enumName = selEnum->GetName();
            if (strcmp(enumScope,"") == 0){
               // This is a global enum and is added to the
               // list of enums and its constants to the list of globals
               if (!listOfEnums->THashList::FindObject(enumName)){
                  ((TEnum*) selEnum)->SetClass(nullptr);
                  listOfEnums->Add(selEnum);
               }
               for (auto enumConstant: *static_cast<TEnum*>(selEnum)->GetConstants()){
                  if (!listOfGlobals->FindObject(enumConstant)){
                     listOfGlobals->Add(enumConstant);
                  }
               }
            }
            else {
               // This enum is in a namespace. A TClass entry is bootstrapped if
               // none exists yet and the enum is added to it
               TClass* nsTClassEntry = TClass::GetClass(enumScope);
               if (!nsTClassEntry){
                  nsTClassEntry = new TClass(enumScope,0,TClass::kNamespaceForMeta, true);
               }
               auto listOfEnums = nsTClassEntry->fEnums.load();
               if (!listOfEnums) {
                  if ( (kIsClass | kIsStruct | kIsUnion) & nsTClassEntry->Property() ) {
                     // For this case, the list will be immutable once constructed
                     // (i.e. in this case, by the end of this routine).
                     listOfEnums = nsTClassEntry->fEnums = new TListOfEnums(nsTClassEntry);
                  } else {
                     //namespaces can have enums added to them
                     listOfEnums = nsTClassEntry->fEnums = new TListOfEnumsWithLock(nsTClassEntry);
                  }
               }
               if (listOfEnums && !listOfEnums->THashList::FindObject(enumName)){
                  ((TEnum*) selEnum)->SetClass(nsTClassEntry);
                  listOfEnums->Add(selEnum);
               }
            }
         }
         enums->Clear();
         delete enums;
      }

      delete pcmFile;

      gDebug = oldDebug;
   } else {
      if (gDebug > 5)
         ::Info("TCling::LoadPCM", "Loading clang PCM %s", pcmFileName.Data());

   }
   return kTRUE;
}

//______________________________________________________________________________

namespace {
   using namespace clang;

   class ExtLexicalStorageAdder: public RecursiveASTVisitor<ExtLexicalStorageAdder>{
      // This class is to be considered an helper for autoparsing.
      // It visits the AST and marks all classes (in all of their redeclarations)
      // with the setHasExternalLexicalStorage method.
   public:
      bool VisitRecordDecl(clang::RecordDecl* rcd){
         if (gDebug > 2)
            Info("ExtLexicalStorageAdder",
                 "Adding external lexical storage to class %s",
                 rcd->getNameAsString().c_str());
         auto reDeclPtr = rcd->getMostRecentDecl();
         do {
            reDeclPtr->setHasExternalLexicalStorage();
         } while ((reDeclPtr = reDeclPtr->getPreviousDecl()));

         return false;
      }
   };


}

////////////////////////////////////////////////////////////////////////////////
/// List of dicts that have the PCM information already in the PCH.
static const std::unordered_set<std::string> gIgnoredPCMNames = {"libCore",
                                                                 "libRint",
                                                                 "libThread",
                                                                 "libRIO",
                                                                 "libImt",
                                                                 "libcomplexDict",
                                                                 "libdequeDict",
                                                                 "liblistDict",
                                                                 "libforward_listDict",
                                                                 "libvectorDict",
                                                                 "libmapDict",
                                                                 "libmultimap2Dict",
                                                                 "libmap2Dict",
                                                                 "libmultimapDict",
                                                                 "libsetDict",
                                                                 "libmultisetDict",
                                                                 "libunordered_setDict",
                                                                 "libunordered_multisetDict",
                                                                 "libunordered_mapDict",
                                                                 "libunordered_multimapDict",
                                                                 "libvalarrayDict",
                                                                 "G__GenVector32",
                                                                 "G__Smatrix32"};

////////////////////////////////////////////////////////////////////////////////
/// Inject the module named "modulename" into cling; load all headers.
/// headers is a 0-terminated array of header files to #include after
/// loading the module. The module is searched for in all $LD_LIBRARY_PATH
/// entries (or %PATH% on Windows).
/// This function gets called by the static initialization of dictionary
/// libraries.
/// The payload code is injected "as is" in the interpreter.
/// The value of 'triggerFunc' is used to find the shared library location.

void TCling::RegisterModule(const char* modulename,
                            const char** headers,
                            const char** includePaths,
                            const char* payloadCode,
                            const char* fwdDeclsCode,
                            void (*triggerFunc)(),
                            const FwdDeclArgsToKeepCollection_t& fwdDeclsArgToSkip,
                            const char** classesHeaders,
                            Bool_t lateRegistration /*=false*/,
                            Bool_t hasCxxModule /*=false*/)
{
   const bool fromRootCling = IsFromRootCling();
   // We need the dictionary initialization but we don't want to inject the
   // declarations into the interpreter, except for those we really need for
   // I/O; see rootcling.cxx after the call to TCling__GetInterpreter().
   if (fromRootCling) return;

   // When we cannot provide a module for the library we should enable header
   // parsing. This 'mixed' mode ensures gradual migration to modules.
   llvm::SaveAndRestore<bool> SaveHeaderParsing(fHeaderParsingOnDemand);
   fHeaderParsingOnDemand = !hasCxxModule;

   // Treat Aclic Libs in a special way. Do not delay the parsing.
   bool hasHeaderParsingOnDemand = fHeaderParsingOnDemand;
   bool isACLiC = strstr(modulename, "_ACLiC_dict") != nullptr;
   if (hasHeaderParsingOnDemand && isACLiC) {
      if (gDebug>1)
         Info("TCling::RegisterModule",
              "Header parsing on demand is active but this is an Aclic library. Disabling it for this library.");
      hasHeaderParsingOnDemand = false;
   }


   // Make sure we relookup symbols that were search for before we loaded
   // their autoparse information.  We could be more subtil and remove only
   // the failed one or only the one in this module, but for now this is
   // better than nothing.
   fLookedUpClasses.clear();

   // Make sure we do not set off autoloading or autoparsing during the
   // module registration!
   Int_t oldAutoloadValue = SetClassAutoloading(false);

   for (const char** inclPath = includePaths; *inclPath; ++inclPath) {
      TCling::AddIncludePath(*inclPath);
   }
   cling::Transaction* T = 0;
   // Put the template decls and the number of arguments to skip in the TNormalizedCtxt
   for (auto& fwdDeclArgToSkipPair : fwdDeclsArgToSkip){
      const std::string& fwdDecl = fwdDeclArgToSkipPair.first;
      const int nArgsToSkip = fwdDeclArgToSkipPair.second;
      auto compRes = fInterpreter->declare(fwdDecl.c_str(), &T);
      assert(cling::Interpreter::kSuccess == compRes &&
            "A fwd declaration could not be compiled");
      if (compRes!=cling::Interpreter::kSuccess){
         Warning("TCling::RegisterModule",
               "Problems in declaring string '%s' were encountered.",
               fwdDecl.c_str()) ;
         continue;
      }

      // Drill through namespaces recursively until the template is found
      if(ClassTemplateDecl* TD = FindTemplateInNamespace(T->getFirstDecl().getSingleDecl())){
         fNormalizedCtxt->AddTemplAndNargsToKeep(TD->getCanonicalDecl(), nArgsToSkip);
      }

   }

   // FIXME: Remove #define __ROOTCLING__ once PCMs are there.
   // This is used to give Sema the same view on ACLiC'ed files (which
   // are then #included through the dictionary) as rootcling had.
   TString code = gNonInterpreterClassDef;
   if (payloadCode)
      code += payloadCode;

   const char* dyLibName = nullptr;
   // If this call happens after dlopen has finished (i.e. late registration)
   // there is no need to dlopen the library recursively. See ROOT-8437 where
   // the dyLibName would correspond to the binary.
   if (!lateRegistration) {
      // We need to open the dictionary shared library, to resolve symbols
      // requested by the JIT from it: as the library is currently being dlopen'ed,
      // its symbols are not yet reachable from the process.
      // Recursive dlopen seems to work just fine.
      dyLibName = FindLibraryName(triggerFunc);
      if (dyLibName) {
         // We were able to determine the library name.
         void* dyLibHandle = dlopen(dyLibName, RTLD_LAZY | RTLD_GLOBAL);
         if (!dyLibHandle) {
#ifdef R__WIN32
            char dyLibError[1000];
            FormatMessageA(FORMAT_MESSAGE_FROM_SYSTEM, NULL, GetLastError(), MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
                           dyLibError, sizeof(dyLibError), NULL);
#else
            const char* dyLibError = dlerror();
            if (dyLibError)
#endif
            {
               if (gDebug > 0) {
                  ::Info("TCling::RegisterModule",
                         "Cannot open shared library %s for dictionary %s:\n  %s",
                         dyLibName, modulename, dyLibError);
               }
            }
            dyLibName = 0;
         } else {
            fRegisterModuleDyLibs.push_back(dyLibHandle);
         } // if (!dyLibHandle) .. else
      } // if (dyLibName)
   } // if (!lateRegistration)

   if (hasHeaderParsingOnDemand && fwdDeclsCode){
      // We now parse the forward declarations. All the classes are then modified
      // in order for them to have an external lexical storage.
      std::string fwdDeclsCodeLessEnums;
      {
         // Search for enum forward decls and only declare them if no
         // declaration exists yet.
         std::string fwdDeclsLine;
         std::istringstream fwdDeclsCodeStr(fwdDeclsCode);
         std::vector<std::string> scopes;
         while (std::getline(fwdDeclsCodeStr, fwdDeclsLine)) {
            const auto enumPos = fwdDeclsLine.find("enum  __attribute__((annotate(\"");
            // We check if the line contains a fwd declaration of an enum
            if (enumPos != std::string::npos) {
               // We clear the scopes which we may have carried from a previous iteration
               scopes.clear();
               // We check if the enum is not in a scope. If yes, save its name
               // and the names of the enclosing scopes.
               if (enumPos != 0) {
                  // it's enclosed in namespaces. We need to understand what they are
                  auto nsPos = fwdDeclsLine.find("namespace");
                  R__ASSERT(nsPos < enumPos && "Inconsistent enum and enclosing scope parsing!");
                  while (nsPos < enumPos && nsPos != std::string::npos) {
                     // we have a namespace, let's put it in the collection of scopes
                     const auto nsNameStart = nsPos + 10;
                     const auto nsNameEnd = fwdDeclsLine.find('{', nsNameStart);
                     const auto nsName = fwdDeclsLine.substr(nsNameStart, nsNameEnd - nsNameStart);
                     scopes.push_back(nsName);
                     nsPos = fwdDeclsLine.find("namespace", nsNameEnd);
                  }
               }
               clang::DeclContext* DC = 0;
               for (auto &&aScope: scopes) {
                  DC = cling::utils::Lookup::Namespace(&fInterpreter->getSema(), aScope.c_str(), DC);
                  if (!DC) {
                     // No decl context means we have to fwd declare the enum.
                     break;
                  }
               }
               if (scopes.empty() || DC) {
                  // We know the scope; let's look for the enum.
                  size_t posEnumName = fwdDeclsLine.find("\"))) ", 32);
                  R__ASSERT(posEnumName != std::string::npos && "Inconsistent enum fwd decl!");
                  posEnumName += 5; // skip "\"))) "
                  while (isspace(fwdDeclsLine[posEnumName]))
                     ++posEnumName;
                  size_t posEnumNameEnd = fwdDeclsLine.find(" : ", posEnumName);
                  R__ASSERT(posEnumNameEnd  != std::string::npos && "Inconsistent enum fwd decl (end)!");
                  while (isspace(fwdDeclsLine[posEnumNameEnd]))
                     --posEnumNameEnd;
                  // posEnumNameEnd now points to the last character of the name.

                  std::string enumName = fwdDeclsLine.substr(posEnumName,
                                                             posEnumNameEnd - posEnumName + 1);

                  if (clang::NamedDecl* enumDecl
                      = cling::utils::Lookup::Named(&fInterpreter->getSema(),
                                                    enumName.c_str(), DC)) {
                     // We have an existing enum decl (forward or definition);
                     // skip this.
                     R__ASSERT(llvm::dyn_cast<clang::EnumDecl>(enumDecl) && "not an enum decl!");
                     (void)enumDecl;
                     continue;
                  }
               }
            }

            fwdDeclsCodeLessEnums += fwdDeclsLine + "\n";
         }
      }

      if (fwdDeclsCodeLessEnums.size() != 0){ // Avoid the overhead if nothing is to be declared
         auto compRes = fInterpreter->declare(fwdDeclsCodeLessEnums, &T);
         assert(cling::Interpreter::kSuccess == compRes &&
               "The forward declarations could not be compiled");
         if (compRes!=cling::Interpreter::kSuccess){
            Warning("TCling::RegisterModule",
                  "Problems in compiling forward declarations for module %s: '%s'",
                  modulename, fwdDeclsCodeLessEnums.c_str()) ;
         }
         else if (T){
            // Loop over all decls in the transaction and go through them all
            // to mark them properly.
            // In order to do that, we first iterate over all the DelayedCallInfos
            // within the transaction. Then we loop over all Decls in the DeclGroupRef
            // contained in the DelayedCallInfos. For each decl, we traverse.
            ExtLexicalStorageAdder elsa;
            for (auto dciIt = T->decls_begin();dciIt!=T->decls_end();dciIt++){
               cling::Transaction::DelayCallInfo& dci = *dciIt;
               for(auto dit = dci.m_DGR.begin(); dit != dci.m_DGR.end(); ++dit) {
                  clang::Decl* declPtr = *dit;
                  elsa.TraverseDecl(declPtr);
               }
            }
         }
      }

      // Now we register all the headers necessary for the class
      // Typical format of the array:
      //    {"A", "classes.h", "@",
      //     "vector<A>", "vector", "@",
      //     "myClass", payloadCode, "@",
      //    nullptr};

      std::string temp;
      for (const char** classesHeader = classesHeaders; *classesHeader; ++classesHeader) {
         temp=*classesHeader;

         size_t theTemplateHash = 0;
         bool addTemplate = false;
         size_t posTemplate = temp.find('<');
         if (posTemplate != std::string::npos) {
            // Add an entry for the template itself.
            std::string templateName = temp.substr(0, posTemplate);
            theTemplateHash = fStringHashFunction(templateName);
            addTemplate = true;
         }
         size_t theHash = fStringHashFunction(temp);
         classesHeader++;
         for (const char** classesHeader_inner = classesHeader; 0!=strcmp(*classesHeader_inner,"@"); ++classesHeader_inner,++classesHeader){
            // This is done in order to distinguish headers from files and from the payloadCode
            if (payloadCode == *classesHeader_inner ){
               fPayloads.insert(theHash);
               if (addTemplate) fPayloads.insert(theTemplateHash);
            }
            if (gDebug > 2)
               Info("TCling::RegisterModule",
                     "Adding a header for %s", temp.c_str());
            fClassesHeadersMap[theHash].push_back(*classesHeader_inner);
            if (addTemplate) {
               if (fClassesHeadersMap.find(theTemplateHash) == fClassesHeadersMap.end()) {
                  fClassesHeadersMap[theTemplateHash].push_back(*classesHeader_inner);
               }
               addTemplate = false;
            }
         }
      }
   }

   if (gIgnoredPCMNames.find(modulename) == gIgnoredPCMNames.end()) {
      TString pcmFileName(ROOT::TMetaUtils::GetModuleFileName(modulename).c_str());
      if (!LoadPCM(pcmFileName, headers, triggerFunc)) {
         ::Error("TCling::RegisterModule", "cannot find dictionary module %s",
                 ROOT::TMetaUtils::GetModuleFileName(modulename).c_str());
      }
   }

   bool oldValue = false;
   if (fClingCallbacks)
     oldValue = SetClassAutoloading(false);

   clang::Sema &TheSema = fInterpreter->getSema();

   bool ModuleWasSuccessfullyLoaded = false;
   if (hasCxxModule) {
      std::string ModuleName = modulename;
      if (llvm::StringRef(modulename).startswith("lib"))
         ModuleName = llvm::StringRef(modulename).substr(3).str();

      // FIXME: We should only complain for modules which we know to exist. For example, we should not complain about
      // modules such as GenVector32 because it needs to fall back to GenVector.
      ModuleWasSuccessfullyLoaded = LoadModule(ModuleName, *fInterpreter, /*Complain=*/ false);
      if (!ModuleWasSuccessfullyLoaded) {
         // Only report if we found the module in the modulemap.
         clang::Preprocessor &PP = TheSema.getPreprocessor();
         clang::HeaderSearch &headerSearch = PP.getHeaderSearchInfo();
         clang::ModuleMap &moduleMap = headerSearch.getModuleMap();
         if (moduleMap.findModule(ModuleName))
            Info("TCling::RegisterModule", "Module %s in modulemap failed to load.", ModuleName.c_str());
      }
   }

   { // scope within which diagnostics are de-activated
   // For now we disable diagnostics because we saw them already at
   // dictionary generation time. That won't be an issue with the PCMs.

      clangDiagSuppr diagSuppr(TheSema.getDiagnostics());

#if defined(R__MUST_REVISIT)
#if R__MUST_REVISIT(6,2)
      Warning("TCling::RegisterModule","Diagnostics suppression should be gone by now.");
#endif
#endif

      if (!ModuleWasSuccessfullyLoaded && !hasHeaderParsingOnDemand){
         SuspendAutoParsing autoParseRaii(this);

         const cling::Transaction* watermark = fInterpreter->getLastTransaction();
         cling::Interpreter::CompilationResult compRes = fInterpreter->parseForModule(code.Data());
         if (isACLiC) {
            // Register an unload point.
            fMetaProcessor->registerUnloadPoint(watermark, headers[0]);
         }

         assert(cling::Interpreter::kSuccess == compRes &&
                        "Payload code of a dictionary could not be parsed correctly.");
         if (compRes!=cling::Interpreter::kSuccess) {
            Warning("TCling::RegisterModule",
                  "Problems declaring payload for module %s.", modulename) ;
         }
      }
   }

   // Now that all the header have been registered/compiled, let's
   // make sure to 'reset' the TClass that have a class init in this module
   // but already had their type information available (using information/header
   // loaded from other modules or from class rules).
   if (!ModuleWasSuccessfullyLoaded && !hasHeaderParsingOnDemand) {
      // This code is likely to be superseded by the similar code in LoadPCM,
      // and have been disabled, (inadvertently or awkwardly) by
      // commit 7903f09f3beea69e82ffba29f59fb2d656a4fd54 (Refactor the routines used for header parsing on demand)
      // whereas it seems that a more semantically correct conditional would have
      // been 'if this module does not have a rootpcm'.
      // Note: this need to be review when the clang pcm are being installed.
      //       #if defined(R__MUST_REVISIT)
      while (!fClassesToUpdate.empty()) {
         TClass *oldcl = fClassesToUpdate.back().first;
         if (oldcl->GetState() != TClass::kHasTClassInit) {
            // if (gDebug > 2) Info("RegisterModule", "Forcing TClass init for %s", oldcl->GetName());
            DictFuncPtr_t dict = fClassesToUpdate.back().second;
            fClassesToUpdate.pop_back();
            // Calling func could manipulate the list so, let maintain the list
            // then call the dictionary function.
            TClass *ncl = dict();
            if (ncl) ncl->PostLoadCheck();
         } else {
            fClassesToUpdate.pop_back();
         }
      }
   }

   if (fClingCallbacks)
     SetClassAutoloading(oldValue);

   if (!ModuleWasSuccessfullyLoaded && !hasHeaderParsingOnDemand) {
      // __ROOTCLING__ might be pulled in through PCH
      fInterpreter->declare("#ifdef __ROOTCLING__\n"
                            "#undef __ROOTCLING__\n"
                            + gInterpreterClassDef +
                            "#endif");
   }

   if (dyLibName) {
      void* dyLibHandle = fRegisterModuleDyLibs.back();
      fRegisterModuleDyLibs.pop_back();
      dlclose(dyLibHandle);
   }

   SetClassAutoloading(oldAutoloadValue);
}

////////////////////////////////////////////////////////////////////////////////
/// Register classes that already existed prior to their dictionary loading
/// and that already had a ClassInfo (and thus would not be refresh via
/// UpdateClassInfo.

void TCling::RegisterTClassUpdate(TClass *oldcl,DictFuncPtr_t dict)
{
   fClassesToUpdate.push_back(std::make_pair(oldcl,dict));
}

////////////////////////////////////////////////////////////////////////////////
/// If the dictionary is loaded, we can remove the class from the list
/// (otherwise the class might be loaded twice).

void TCling::UnRegisterTClassUpdate(const TClass *oldcl)
{
   typedef std::vector<std::pair<TClass*,DictFuncPtr_t> >::iterator iterator;
   iterator stop = fClassesToUpdate.end();
   for(iterator i = fClassesToUpdate.begin();
       i != stop;
       ++i)
   {
      if ( i->first == oldcl ) {
         fClassesToUpdate.erase(i);
         return;
      }
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Let cling process a command line.
///
/// If the command is executed and the error is 0, then the return value
/// is the int value corresponding to the result of the executed command
/// (float and double return values will be truncated).
///

// Method for handling the interpreter exceptions.
// the MetaProcessor is passing in as argument to teh function, because
// cling::Interpreter::CompilationResult is a nested class and it cannot be
// forward declared, thus this method cannot be a static member function
// of TCling.

static int HandleInterpreterException(cling::MetaProcessor* metaProcessor,
                                 const char* input_line,
                                 cling::Interpreter::CompilationResult& compRes,
                                 cling::Value* result)
{
   try {
      return metaProcessor->process(input_line, compRes, result);
   }
   catch (cling::InterpreterException& ex)
   {
      Error("HandleInterpreterException", "%s.\n%s", ex.what(), "Execution of your code was aborted.");
      ex.diagnose();
      compRes = cling::Interpreter::kFailure;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////

bool TCling::DiagnoseIfInterpreterException(const std::exception &e) const
{
   if (auto ie = dynamic_cast<const cling::InterpreterException*>(&e)) {
      ie->diagnose();
      return true;
   }
   return false;
}

////////////////////////////////////////////////////////////////////////////////

Long_t TCling::ProcessLine(const char* line, EErrorCode* error/*=0*/)
{
   // Copy the passed line, it comes from a static buffer in TApplication
   // which can be reentered through the Cling evaluation routines,
   // which would overwrite the static buffer and we would forget what we
   // were doing.
   //
   TString sLine(line);
   if (strstr(line,fantomline)) {
      // End-Of-Line action
      // See the comment (copied from above):
      // It is a "fantom" method to synchronize user keyboard input
      // and ROOT prompt line (for WIN32)
      // and is implemented by
      if (gApplication) {
         if (gApplication->IsCmdThread()) {
            R__LOCKGUARD(fLockProcessLine ? gInterpreterMutex : 0);
            gROOT->SetLineIsProcessing();

            UpdateAllCanvases();

            gROOT->SetLineHasBeenProcessed();
         }
      }
      return 0;
   }

   if (gGlobalMutex && !gInterpreterMutex && fLockProcessLine) {
      gGlobalMutex->Lock();
      if (!gInterpreterMutex)
         gInterpreterMutex = gGlobalMutex->Factory(kTRUE);
      gGlobalMutex->UnLock();
   }
   R__LOCKGUARD_CLING(fLockProcessLine ? gInterpreterMutex : 0);
   gROOT->SetLineIsProcessing();

   struct InterpreterFlagsRAII {
      cling::Interpreter* fInterpreter;
      bool fWasDynamicLookupEnabled;

      InterpreterFlagsRAII(cling::Interpreter* interp):
         fInterpreter(interp),
         fWasDynamicLookupEnabled(interp->isDynamicLookupEnabled())
      {
         fInterpreter->enableDynamicLookup(true);
      }
      ~InterpreterFlagsRAII() {
         fInterpreter->enableDynamicLookup(fWasDynamicLookupEnabled);
         gROOT->SetLineHasBeenProcessed();
      }
   } interpreterFlagsRAII(fInterpreter);

   // A non-zero returned value means the given line was
   // not a complete statement.
   int indent = 0;
   // This will hold the resulting value of the evaluation the given line.
   cling::Value result;
   cling::Interpreter::CompilationResult compRes = cling::Interpreter::kSuccess;
   if (!strncmp(sLine.Data(), ".L", 2) || !strncmp(sLine.Data(), ".x", 2) ||
       !strncmp(sLine.Data(), ".X", 2)) {
      // If there was a trailing "+", then CINT compiled the code above,
      // and we will need to strip the "+" before passing the line to cling.
      TString mod_line(sLine);
      TString aclicMode;
      TString arguments;
      TString io;
      TString fname = gSystem->SplitAclicMode(sLine.Data() + 3,
         aclicMode, arguments, io);
      if (aclicMode.Length()) {
         // Remove the leading '+'
         R__ASSERT(aclicMode[0]=='+' && "ACLiC mode must start with a +");
         aclicMode[0]='k';    // We always want to keep the .so around.
         if (aclicMode[1]=='+') {
            // We have a 2nd +
            aclicMode[1]='f'; // We want to force the recompilation.
         }
         if (!gSystem->CompileMacro(fname,aclicMode)) {
            // ACLiC failed.
            compRes = cling::Interpreter::kFailure;
         } else {
            if (strncmp(sLine.Data(), ".L", 2) != 0) {
               // if execution was requested.

               if (arguments.Length()==0) {
                  arguments = "()";
               }
               // We need to remove the extension.
               Ssiz_t ext = fname.Last('.');
               if (ext != kNPOS) {
                  fname.Remove(ext);
               }
               const char *function = gSystem->BaseName(fname);
               mod_line = function + arguments + io;
               indent = HandleInterpreterException(fMetaProcessor, mod_line, compRes, &result);
            }
         }
      } else {
         // not ACLiC
         size_t unnamedMacroOpenCurly;
         {
            std::string code;
            std::string codeline;
            std::ifstream in(fname);
            while (in) {
               std::getline(in, codeline);
               code += codeline + "\n";
            }
            unnamedMacroOpenCurly
              = cling::utils::isUnnamedMacro(code, fInterpreter->getCI()->getLangOpts());
         }

         fCurExecutingMacros.push_back(fname);
         if (unnamedMacroOpenCurly != std::string::npos) {
            compRes = fMetaProcessor->readInputFromFile(fname.Data(), &result,
                                                        unnamedMacroOpenCurly);
         } else {
            // No DynLookup for .x, .L of named macros.
            fInterpreter->enableDynamicLookup(false);
            indent = HandleInterpreterException(fMetaProcessor, mod_line, compRes, &result);
         }
         fCurExecutingMacros.pop_back();
      }
   } // .L / .X / .x
   else {
      if (0!=strncmp(sLine.Data(), ".autodict ",10) && sLine != ".autodict") {
         // explicitly ignore .autodict without having to support it
         // in cling.

         // Turn off autoparsing if this is an include directive
         bool isInclusionDirective = sLine.Contains("\n#include") || sLine.BeginsWith("#include");
         if (isInclusionDirective) {
            SuspendAutoParsing autoParseRaii(this);
            indent = HandleInterpreterException(fMetaProcessor, sLine, compRes, &result);
         } else {
            indent = HandleInterpreterException(fMetaProcessor, sLine, compRes, &result);
         }
      }
   }
   if (result.isValid())
      RegisterTemporary(result);
   if (indent) {
      if (error)
         *error = kProcessing;
      return 0;
   }
   if (error) {
      switch (compRes) {
      case cling::Interpreter::kSuccess: *error = kNoError; break;
      case cling::Interpreter::kFailure: *error = kRecoverable; break;
      case cling::Interpreter::kMoreInputExpected: *error = kProcessing; break;
      }
   }
   if (compRes == cling::Interpreter::kSuccess
       && result.isValid()
       && !result.isVoid())
   {
      return result.simplisticCastAs<long>();
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// No-op; see TRint instead.

void TCling::PrintIntro()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Add the given path to the list of directories in which the interpreter
/// looks for include files. Only one path item can be specified at a
/// time, i.e. "path1:path2" is NOT supported.

void TCling::AddIncludePath(const char *path)
{
   R__LOCKGUARD(gInterpreterMutex);
   // Favorite source of annoyance: gSystem->AddIncludePath() needs "-I",
   // gCling->AddIncludePath() does not! Work around that inconsistency:
   if (path[0] == '-' && path[1] == 'I')
      path += 2;

   fInterpreter->AddIncludePath(path);
}

////////////////////////////////////////////////////////////////////////////////
/// Visit all members over members, recursing over base classes.

void TCling::InspectMembers(TMemberInspector& insp, const void* obj,
                            const TClass* cl, Bool_t isTransient)
{
   if (insp.GetObjectValidity() == TMemberInspector::kUnset) {
      insp.SetObjectValidity(obj ? TMemberInspector::kValidObjectGiven
                             : TMemberInspector::kNoObjectGiven);
   }

   if (!cl || cl->GetCollectionProxy()) {
      // We do not need to investigate the content of the STL
      // collection, they are opaque to us (and details are
      // uninteresting).
      return;
   }

   static const TClassRef clRefString("std::string");
   if (clRefString == cl) {
      // We stream std::string without going through members..
      return;
   }

   if (TClassEdit::IsStdArray(cl->GetName())) {
      // We treat std arrays as C arrays
      return;
   }

   const char* cobj = (const char*) obj; // for ptr arithmetics

   // Treat the case of std::complex in a special manner. We want to enforce
   // the layout of a stl implementation independent class, which is the
   // complex as implemented in ROOT5.

   // A simple lambda to simplify the code
   auto inspInspect =  [&] (ptrdiff_t offset){
      insp.Inspect(const_cast<TClass*>(cl), insp.GetParent(), "_real", cobj, isTransient);
      insp.Inspect(const_cast<TClass*>(cl), insp.GetParent(), "_imag", cobj + offset, isTransient);
   };

   auto complexType = TClassEdit::GetComplexType(cl->GetName());
   switch(complexType) {
      case TClassEdit::EComplexType::kNone:
      {
        break;
      }
      case TClassEdit::EComplexType::kFloat:
      {
        inspInspect(sizeof(float));
        return;
      }
      case TClassEdit::EComplexType::kDouble:
      {
        inspInspect(sizeof(double));
        return;
      }
      case TClassEdit::EComplexType::kInt:
      {
        inspInspect(sizeof(int));
        return;
      }
      case TClassEdit::EComplexType::kLong:
      {
        inspInspect(sizeof(long));
        return;
      }
   }

   static clang::PrintingPolicy
      printPol(fInterpreter->getCI()->getLangOpts());
   if (printPol.Indentation) {
      // not yet initialized
      printPol.Indentation = 0;
      printPol.SuppressInitializers = true;
   }

   const char* clname = cl->GetName();
   // Printf("Inspecting class %s\n", clname);

   const clang::ASTContext& astContext = fInterpreter->getCI()->getASTContext();
   const clang::Decl *scopeDecl = 0;
   const clang::Type *recordType = 0;

   if (cl->GetClassInfo()) {
      TClingClassInfo * clingCI = (TClingClassInfo *)cl->GetClassInfo();
      scopeDecl = clingCI->GetDecl();
      recordType = clingCI->GetType();
   } else {
      const cling::LookupHelper& lh = fInterpreter->getLookupHelper();
      // Diags will complain about private classes:
      scopeDecl = lh.findScope(clname, cling::LookupHelper::NoDiagnostics,
                               &recordType);
   }
   if (!scopeDecl) {
      Error("InspectMembers", "Cannot find Decl for class %s", clname);
      return;
   }
   const clang::CXXRecordDecl* recordDecl
     = llvm::dyn_cast<const clang::CXXRecordDecl>(scopeDecl);
   if (!recordDecl) {
      Error("InspectMembers", "Cannot find Decl for class %s is not a CXXRecordDecl.", clname);
      return;
   }

   {
      // Force possible deserializations first. We need to have no pending
      // Transaction when passing control flow to the inspector below (ROOT-7779).
      cling::Interpreter::PushTransactionRAII deserRAII(fInterpreter);

      astContext.getASTRecordLayout(recordDecl);

      for (clang::RecordDecl::field_iterator iField = recordDecl->field_begin(),
              eField = recordDecl->field_end(); iField != eField; ++iField) {}
   }

   const clang::ASTRecordLayout& recLayout
      = astContext.getASTRecordLayout(recordDecl);

   // TVirtualCollectionProxy *proxy = cl->GetCollectionProxy();
   // if (proxy && ( proxy->GetProperties() & TVirtualCollectionProxy::kIsEmulated ) ) {
   //    Error("InspectMembers","The TClass for %s has an emulated proxy but we are looking at a compiled version of the collection!\n",
   //          cl->GetName());
   // }
   if (cl->Size() != recLayout.getSize().getQuantity()) {
      Error("InspectMembers","TClass and cling disagree on the size of the class %s, respectively %d %lld\n",
            cl->GetName(),cl->Size(),(Long64_t)recLayout.getSize().getQuantity());
   }

   unsigned iNField = 0;
   // iterate over fields
   // FieldDecls are non-static, else it would be a VarDecl.
   for (clang::RecordDecl::field_iterator iField = recordDecl->field_begin(),
        eField = recordDecl->field_end(); iField != eField;
        ++iField, ++iNField) {


      clang::QualType memberQT = iField->getType();
      if (recordType) {
         // if (we_need_to_do_the_subst_because_the_class_is_a_template_instance_of_double32_t)
         memberQT = ROOT::TMetaUtils::ReSubstTemplateArg(memberQT, recordType);
      }
      memberQT = cling::utils::Transform::GetPartiallyDesugaredType(astContext, memberQT, fNormalizedCtxt->GetConfig(), false /* fully qualify */);
      if (memberQT.isNull()) {
         std::string memberName;
         llvm::raw_string_ostream stream(memberName);
         // Don't trigger fopen of the source file to count lines:
         printPol.AnonymousTagLocations = false;
         iField->getNameForDiagnostic(stream, printPol, true /*fqi*/);
         stream.flush();
         Error("InspectMembers",
               "Cannot retrieve QualType for member %s while inspecting class %s",
               memberName.c_str(), clname);
         continue; // skip member
      }
      const clang::Type* memType = memberQT.getTypePtr();
      if (!memType) {
         std::string memberName;
         llvm::raw_string_ostream stream(memberName);
         // Don't trigger fopen of the source file to count lines:
         printPol.AnonymousTagLocations = false;
         iField->getNameForDiagnostic(stream, printPol, true /*fqi*/);
         stream.flush();
         Error("InspectMembers",
               "Cannot retrieve Type for member %s while inspecting class %s",
               memberName.c_str(), clname);
         continue; // skip member
      }

      const clang::Type* memNonPtrType = memType;
      Bool_t ispointer = false;
      if (memNonPtrType->isPointerType()) {
         ispointer = true;
         clang::QualType ptrQT
            = memNonPtrType->getAs<clang::PointerType>()->getPointeeType();
         if (recordType) {
            // if (we_need_to_do_the_subst_because_the_class_is_a_template_instance_of_double32_t)
            ptrQT = ROOT::TMetaUtils::ReSubstTemplateArg(ptrQT, recordType);
         }
         ptrQT = cling::utils::Transform::GetPartiallyDesugaredType(astContext, ptrQT, fNormalizedCtxt->GetConfig(), false /* fully qualify */);
         if (ptrQT.isNull()) {
            std::string memberName;
            llvm::raw_string_ostream stream(memberName);
            // Don't trigger fopen of the source file to count lines:
            printPol.AnonymousTagLocations = false;
            iField->getNameForDiagnostic(stream, printPol, true /*fqi*/);
            stream.flush();
            Error("InspectMembers",
                  "Cannot retrieve pointee Type for member %s while inspecting class %s",
                  memberName.c_str(), clname);
            continue; // skip member
         }
         memNonPtrType = ptrQT.getTypePtr();
      }

      // assemble array size(s): "[12][4][]"
      llvm::SmallString<8> arraySize;
      const clang::ArrayType* arrType = memNonPtrType->getAsArrayTypeUnsafe();
      unsigned arrLevel = 0;
      bool haveErrorDueToArray = false;
      while (arrType) {
         ++arrLevel;
         arraySize += '[';
         const clang::ConstantArrayType* constArrType =
         clang::dyn_cast<clang::ConstantArrayType>(arrType);
         if (constArrType) {
            constArrType->getSize().toStringUnsigned(arraySize);
         }
         arraySize += ']';
         clang::QualType subArrQT = arrType->getElementType();
         if (subArrQT.isNull()) {
            std::string memberName;
            llvm::raw_string_ostream stream(memberName);
            // Don't trigger fopen of the source file to count lines:
            printPol.AnonymousTagLocations = false;
            iField->getNameForDiagnostic(stream, printPol, true /*fqi*/);
            stream.flush();
            Error("InspectMembers",
                  "Cannot retrieve QualType for array level %d (i.e. element type of %s) for member %s while inspecting class %s",
                  arrLevel, subArrQT.getAsString(printPol).c_str(),
                  memberName.c_str(), clname);
            haveErrorDueToArray = true;
            break;
         }
         arrType = subArrQT.getTypePtr()->getAsArrayTypeUnsafe();
      }
      if (haveErrorDueToArray) {
         continue; // skip member
      }

      // construct member name
      std::string fieldName;
      if (memType->isPointerType()) {
         fieldName = "*";
      }

      // Check if this field has a custom ioname, if not, just use the one of the decl
      std::string ioname(iField->getName());
      ROOT::TMetaUtils::ExtractAttrPropertyFromName(**iField,"ioname",ioname);
      fieldName += ioname;
      fieldName += arraySize;

      // get member offset
      // NOTE currently we do not support bitfield and do not support
      // member that are not aligned on 'bit' boundaries.
      clang::CharUnits offset(astContext.toCharUnitsFromBits(recLayout.getFieldOffset(iNField)));
      ptrdiff_t fieldOffset = offset.getQuantity();

      // R__insp.Inspect(R__cl, R__insp.GetParent(), "fBits[2]", fBits);
      // R__insp.Inspect(R__cl, R__insp.GetParent(), "fName", &fName);
      // R__insp.InspectMember(fName, "fName.");
      // R__insp.Inspect(R__cl, R__insp.GetParent(), "*fClass", &fClass);

      // If the class has a custom streamer and the type of the filed is a
      // private enum, struct or class, skip it.
      if (!insp.IsTreatingNonAccessibleTypes()){
         auto iFiledQtype = iField->getType();
         if (auto tagDecl = iFiledQtype->getAsTagDecl()){
            auto declAccess = tagDecl->getAccess();
            if (declAccess == AS_private || declAccess == AS_protected) {
               continue;
            }
         }
      }

      insp.Inspect(const_cast<TClass*>(cl), insp.GetParent(), fieldName.c_str(), cobj + fieldOffset, isTransient);

      if (!ispointer) {
         const clang::CXXRecordDecl* fieldRecDecl = memNonPtrType->getAsCXXRecordDecl();
         if (fieldRecDecl) {
            // nested objects get an extra call to InspectMember
            // R__insp.InspectMember("FileStat_t", (void*)&fFileStat, "fFileStat.", false);
            std::string sFieldRecName;
            if (!ROOT::TMetaUtils::ExtractAttrPropertyFromName(*fieldRecDecl,"iotype",sFieldRecName)){
               ROOT::TMetaUtils::GetNormalizedName(sFieldRecName,
                                                   clang::QualType(memNonPtrType,0),
                                                   *fInterpreter,
                                                   *fNormalizedCtxt);
            }

            TDataMember* mbr = cl->GetDataMember(ioname.c_str());
            // if we can not find the member (which should not really happen),
            // let's consider it transient.
            Bool_t transient = isTransient || !mbr || !mbr->IsPersistent();

            insp.InspectMember(sFieldRecName.c_str(), cobj + fieldOffset,
                               (fieldName + '.').c_str(), transient);

         }
      }
   } // loop over fields

   // inspect bases
   // TNamed::ShowMembers(R__insp);
   unsigned iNBase = 0;
   for (clang::CXXRecordDecl::base_class_const_iterator iBase
        = recordDecl->bases_begin(), eBase = recordDecl->bases_end();
        iBase != eBase; ++iBase, ++iNBase) {
      clang::QualType baseQT = iBase->getType();
      if (baseQT.isNull()) {
         Error("InspectMembers",
               "Cannot find QualType for base number %d while inspecting class %s",
               iNBase, clname);
         continue;
      }
      const clang::CXXRecordDecl* baseDecl
         = baseQT->getAsCXXRecordDecl();
      if (!baseDecl) {
         Error("InspectMembers",
               "Cannot find CXXRecordDecl for base number %d while inspecting class %s",
               iNBase, clname);
         continue;
      }
      TClass* baseCl=nullptr;
      std::string sBaseName;
      // Try with the DeclId
      std::vector<TClass*> foundClasses;
      TClass::GetClass(static_cast<DeclId_t>(baseDecl), foundClasses);
      if (foundClasses.size()==1){
         baseCl=foundClasses[0];
      } else {
         // Try with the normalised Name, as a fallback
         if (!baseCl){
            ROOT::TMetaUtils::GetNormalizedName(sBaseName,
                                                baseQT,
                                                *fInterpreter,
                                                *fNormalizedCtxt);
            baseCl = TClass::GetClass(sBaseName.c_str());
         }
      }

      if (!baseCl){
         std::string qualNameForDiag;
         ROOT::TMetaUtils::GetQualifiedName(qualNameForDiag, *baseDecl);
         Error("InspectMembers",
               "Cannot find TClass for base class %s", qualNameForDiag.c_str() );
         continue;
      }

      int64_t baseOffset;
      if (iBase->isVirtual()) {
         if (insp.GetObjectValidity() == TMemberInspector::kNoObjectGiven) {
            if (!isTransient) {
               Error("InspectMembers",
                     "Base %s of class %s is virtual but no object provided",
                     sBaseName.c_str(), clname);
            }
            baseOffset = TVirtualStreamerInfo::kNeedObjectForVirtualBaseClass;
         } else {
            // We have an object to determine the vbase offset.
            TClingClassInfo* ci = (TClingClassInfo*)cl->GetClassInfo();
            TClingClassInfo* baseCi = (TClingClassInfo*)baseCl->GetClassInfo();
            if (ci && baseCi) {
               baseOffset = ci->GetBaseOffset(baseCi, const_cast<void*>(obj),
                                              true /*isDerivedObj*/);
               if (baseOffset == -1) {
                  Error("InspectMembers",
                        "Error calculating offset of virtual base %s of class %s",
                        sBaseName.c_str(), clname);
               }
            } else {
               Error("InspectMembers",
                     "Cannot calculate offset of virtual base %s of class %s",
                     sBaseName.c_str(), clname);
               continue;
            }
         }
      } else {
         baseOffset = recLayout.getBaseClassOffset(baseDecl).getQuantity();
      }
      // TOFIX: baseCl can be null here!
      if (baseCl->IsLoaded()) {
         // For loaded class, CallShowMember will (especially for TObject)
         // call the virtual ShowMember rather than the class specific version
         // resulting in an infinite recursion.
         InspectMembers(insp, cobj + baseOffset, baseCl, isTransient);
      } else {
         baseCl->CallShowMembers(cobj + baseOffset,
                                 insp, isTransient);
      }
   } // loop over bases
}

////////////////////////////////////////////////////////////////////////////////
/// Reset the interpreter internal state in case a previous action was not correctly
/// terminated.

void TCling::ClearFileBusy()
{
   // No-op there is not equivalent state (to be cleared) in Cling.
}

////////////////////////////////////////////////////////////////////////////////
/// Delete existing temporary values.

void TCling::ClearStack()
{
   // No-op for cling due to cling::Value.
}

////////////////////////////////////////////////////////////////////////////////
/// Declare code to the interpreter, without any of the interpreter actions
/// that could trigger a re-interpretation of the code. I.e. make cling
/// behave like a compiler: no dynamic lookup, no input wrapping for
/// subsequent execution, no automatic provision of declarations but just a
/// plain #include.
/// Returns true on success, false on failure.

bool TCling::Declare(const char* code)
{
   R__LOCKGUARD_CLING(gInterpreterMutex);

   int oldload = SetClassAutoloading(0);
   SuspendAutoParsing autoParseRaii(this);

   bool oldDynLookup = fInterpreter->isDynamicLookupEnabled();
   fInterpreter->enableDynamicLookup(false);
   bool oldRawInput = fInterpreter->isRawInputEnabled();
   fInterpreter->enableRawInput(true);

   Bool_t ret = LoadText(code);

   fInterpreter->enableRawInput(oldRawInput);
   fInterpreter->enableDynamicLookup(oldDynLookup);
   SetClassAutoloading(oldload);
   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// Enable the automatic loading of shared libraries when a class
/// is used that is stored in a not yet loaded library. Uses the
/// information stored in the class/library map (typically
/// $ROOTSYS/etc/system.rootmap).

void TCling::EnableAutoLoading()
{
   if (fAllowLibLoad) {
      LoadLibraryMap();
      SetClassAutoloading(true);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// It calls a "fantom" method to synchronize user keyboard input
/// and ROOT prompt line.

void TCling::EndOfLineAction()
{
   ProcessLineSynch(fantomline);
}

// This static function is a hop of TCling::IsLibraryLoaded, which is taking a lock and calling
// into this function. This is because we wanted to avoid a duplication in TCling::IsLoaded, which
// was already taking a lock.
static Bool_t s_IsLibraryLoaded(const char* libname, cling::Interpreter* fInterpreter)
{
   // Check shared library.
   TString tLibName(libname);
   if (gSystem->FindDynamicLibrary(tLibName, kTRUE))
      return fInterpreter->getDynamicLibraryManager()->isLibraryLoaded(tLibName.Data());
   return false;
}

Bool_t TCling::IsLibraryLoaded(const char* libname) const
{
   R__LOCKGUARD(gInterpreterMutex);
   return s_IsLibraryLoaded(libname, fInterpreter);
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if ROOT has cxxmodules pcm for a given library name.
// FIXME: We need to be able to support lazy loading of pcm generated by ACLiC.
Bool_t TCling::HasPCMForLibrary(const char *libname) const
{
   llvm::StringRef ModuleName(libname);
   ModuleName = llvm::sys::path::stem(ModuleName);
   ModuleName.consume_front("lib");

   clang::ModuleMap &moduleMap = fInterpreter->getCI()->getPreprocessor().getHeaderSearchInfo().getModuleMap();
   clang::Module *M = moduleMap.findModule(ModuleName);
   return M && !M->IsMissingRequirement && M->getASTFile();
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if the file has already been loaded by cint.
/// We will try in this order:
///   actual filename
///   filename as a path relative to
///            the include path
///            the shared library path

Bool_t TCling::IsLoaded(const char* filename) const
{
   R__LOCKGUARD(gInterpreterMutex);

   //FIXME: if we use llvm::sys::fs::make_absolute all this can go away. See
   // cling::DynamicLibraryManager.

   std::string file_name = filename;
   size_t at = std::string::npos;
   while ((at = file_name.find("/./")) != std::string::npos)
       file_name.replace(at, 3, "/");

   std::string filesStr = "";
   llvm::raw_string_ostream filesOS(filesStr);
   clang::SourceManager &SM = fInterpreter->getCI()->getSourceManager();
   cling::ClangInternalState::printIncludedFiles(filesOS, SM);
   filesOS.flush();

   llvm::SmallVector<llvm::StringRef, 100> files;
   llvm::StringRef(filesStr).split(files, "\n");

   std::set<std::string> fileMap;
   // Fill fileMap; return early on exact match.
   for (llvm::SmallVector<llvm::StringRef, 100>::const_iterator
           iF = files.begin(), iE = files.end(); iF != iE; ++iF) {
      if ((*iF) == file_name.c_str()) return kTRUE; // exact match
      fileMap.insert(*iF);
   }

   if (fileMap.empty()) return kFALSE;

   // Check MacroPath.
   TString sFilename(file_name.c_str());
   if (gSystem->FindFile(TROOT::GetMacroPath(), sFilename, kReadPermission)
       && fileMap.count(sFilename.Data())) {
      return kTRUE;
   }

   // Check IncludePath.
   TString incPath = gSystem->GetIncludePath(); // of the form -Idir1  -Idir2 -Idir3
   incPath.Append(":").Prepend(" "); // to match " -I" (note leading ' ')
   incPath.ReplaceAll(" -I", ":");      // of form :dir1 :dir2:dir3
   while (incPath.Index(" :") != -1) {
      incPath.ReplaceAll(" :", ":");
   }
   incPath.Prepend(".:");
   sFilename = file_name.c_str();
   if (gSystem->FindFile(incPath, sFilename, kReadPermission)
       && fileMap.count(sFilename.Data())) {
      return kTRUE;
   }

   // Check shared library.
   if (s_IsLibraryLoaded(file_name.c_str(), fInterpreter))
      return kTRUE;

   //FIXME: We must use the cling::Interpreter::lookupFileOrLibrary iface.
   const clang::DirectoryLookup *CurDir = 0;
   clang::Preprocessor &PP = fInterpreter->getCI()->getPreprocessor();
   clang::HeaderSearch &HS = PP.getHeaderSearchInfo();
   const clang::FileEntry *FE = HS.LookupFile(file_name.c_str(),
                                              clang::SourceLocation(),
                                              /*isAngled*/ false,
                                              /*FromDir*/ 0, CurDir,
                                              clang::ArrayRef<std::pair<const clang::FileEntry *,
                                                                        const clang::DirectoryEntry *>>(),
                                              /*SearchPath*/ 0,
                                              /*RelativePath*/ 0,
                                              /*RequestingModule*/ 0,
                                              /*SuggestedModule*/ 0,
                                              /*IsMapped*/ 0,
                                              /*SkipCache*/ false,
                                              /*BuildSystemModule*/ false,
                                              /*OpenFile*/ false,
                                              /*CacheFail*/ false);
   if (FE && FE->isValid()) {
      // check in the source manager if the file is actually loaded
      clang::SourceManager &SM = fInterpreter->getCI()->getSourceManager();
      // this works only with header (and source) files...
      clang::FileID FID = SM.translateFile(FE);
      if (!FID.isInvalid() && FID.getHashValue() == 0)
         return kFALSE;
      else {
         clang::SrcMgr::SLocEntry SLocE = SM.getSLocEntry(FID);
         if (SLocE.isFile() && SLocE.getFile().getContentCache()->getRawBuffer() == 0)
            return kFALSE;
         if (!FID.isInvalid())
            return kTRUE;
      }
      // ...then check shared library again, but with full path now
      sFilename = FE->getName();
      if (gSystem->FindDynamicLibrary(sFilename, kTRUE)
          && fileMap.count(sFilename.Data())) {
         return kTRUE;
      }
   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////

void TCling::UpdateListOfLoadedSharedLibraries()
{
#if defined(R__WIN32) || defined(__CYGWIN__)
   HMODULE hModules[1024];
   void *hProcess;
   unsigned long cbModules;
   unsigned int i;
   hProcess = (void *)::GetCurrentProcess();
   ::EnumProcessModules(hProcess, hModules, sizeof(hModules), &cbModules);
   // start at 1 to skip the executable itself
   for (i = 1; i < (cbModules / sizeof(void *)); i++) {
      static const int bufsize = 260;
      wchar_t winname[bufsize];
      char posixname[bufsize];
      ::GetModuleFileNameExW(hProcess, hModules[i], winname, bufsize);
#if defined(__CYGWIN__)
      cygwin_conv_path(CCP_WIN_W_TO_POSIX, winname, posixname, bufsize);
#else
      std::wstring wpath = winname;
      std::replace(wpath.begin(), wpath.end(), '\\', '/');
      string path(wpath.begin(), wpath.end());
      strncpy(posixname, path.c_str(), bufsize);
#endif
      if (!fSharedLibs.Contains(posixname)) {
         RegisterLoadedSharedLibrary(posixname);
      }
   }
#elif defined(R__MACOSX)
   // fPrevLoadedDynLibInfo stores the *next* image index to look at
   uint32_t imageIndex = (uint32_t) (size_t) fPrevLoadedDynLibInfo;

   while (const mach_header* mh = _dyld_get_image_header(imageIndex)) {
      // Skip non-dylibs
      if (mh->filetype == MH_DYLIB) {
         if (const char* imageName = _dyld_get_image_name(imageIndex)) {
            RegisterLoadedSharedLibrary(imageName);
         }
      }

      ++imageIndex;
   }
   fPrevLoadedDynLibInfo = (void*)(size_t)imageIndex;
#elif defined(R__LINUX)
   struct PointerNo4 {
      void* fSkip[3];
      void* fPtr;
   };
   struct LinkMap {
      void* fAddr;
      const char* fName;
      void* fLd;
      LinkMap* fNext;
      LinkMap* fPrev;
   };
   if (!fPrevLoadedDynLibInfo || fPrevLoadedDynLibInfo == (void*)(size_t)-1) {
      PointerNo4* procLinkMap = (PointerNo4*)dlopen(0,  RTLD_LAZY | RTLD_GLOBAL);
      // 4th pointer of 4th pointer is the linkmap.
      // See http://syprog.blogspot.fr/2011/12/listing-loaded-shared-objects-in-linux.html
      LinkMap* linkMap = (LinkMap*) ((PointerNo4*)procLinkMap->fPtr)->fPtr;
      RegisterLoadedSharedLibrary(linkMap->fName);
      fPrevLoadedDynLibInfo = linkMap;
      // reduce use count of link map structure:
      dlclose(procLinkMap);
   }

   LinkMap* iDyLib = (LinkMap*)fPrevLoadedDynLibInfo;
   while (iDyLib->fNext) {
      iDyLib = iDyLib->fNext;
      RegisterLoadedSharedLibrary(iDyLib->fName);
   }
   fPrevLoadedDynLibInfo = iDyLib;
#else
   Error("TCling::UpdateListOfLoadedSharedLibraries",
         "Platform not supported!");
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Register a new shared library name with the interpreter; add it to
/// fSharedLibs.

void TCling::RegisterLoadedSharedLibrary(const char* filename)
{
   // Ignore NULL filenames, aka "the process".
   if (!filename) return;

   // Tell the interpreter that this library is available; all libraries can be
   // used to resolve symbols.
   cling::DynamicLibraryManager* DLM = fInterpreter->getDynamicLibraryManager();
   if (!DLM->isLibraryLoaded(filename)) {
      DLM->loadLibrary(filename, true /*permanent*/);
   }

#if defined(R__MACOSX)
   // Check that this is not a system library
   if (!strncmp(filename, "/usr/lib/system/", 16)
       || !strncmp(filename, "/usr/lib/libc++", 15)
       || !strncmp(filename, "/System/Library/Frameworks/", 27)
       || !strncmp(filename, "/System/Library/PrivateFrameworks/", 34)
       || !strncmp(filename, "/System/Library/CoreServices/", 29)
       || !strcmp(filename, "cl_kernels") // yepp, no directory
       || strstr(filename, "/usr/lib/libSystem")
       || strstr(filename, "/usr/lib/libstdc++")
       || strstr(filename, "/usr/lib/libicucore")
       || strstr(filename, "/usr/lib/libbsm")
       || strstr(filename, "/usr/lib/libobjc")
       || strstr(filename, "/usr/lib/libresolv")
       || strstr(filename, "/usr/lib/libauto")
       || strstr(filename, "/usr/lib/libcups")
       || strstr(filename, "/usr/lib/libDiagnosticMessagesClient")
       || strstr(filename, "/usr/lib/liblangid")
       || strstr(filename, "/usr/lib/libCRFSuite")
       || strstr(filename, "/usr/lib/libpam")
       || strstr(filename, "/usr/lib/libOpenScriptingUtil")
       || strstr(filename, "/usr/lib/libextension"))
      return;
#elif defined(__CYGWIN__)
   // Check that this is not a system library
   static const int bufsize = 260;
   char posixwindir[bufsize];
   char *windir = getenv("WINDIR");
   if (windir)
      cygwin_conv_path(CCP_WIN_A_TO_POSIX, windir, posixwindir, bufsize);
   else
      snprintf(posixwindir, sizeof(posixwindir), "/Windows/");
   if (strstr(filename, posixwindir) ||
       strstr(filename, "/usr/bin/cyg"))
      return;
#elif defined(R__WIN32)
   if (strstr(filename, "/Windows/"))
      return;
#elif defined (R__LINUX)
   if (strstr(filename, "/ld-linux")
       || strstr(filename, "linux-gnu/")
       || strstr(filename, "/libstdc++.")
       || strstr(filename, "/libgcc")
       || strstr(filename, "/libc.")
       || strstr(filename, "/libdl.")
       || strstr(filename, "/libm."))
      return;
#endif
   // Update string of available libraries.
   if (!fSharedLibs.IsNull()) {
      fSharedLibs.Append(" ");
   }
   fSharedLibs.Append(filename);
}

////////////////////////////////////////////////////////////////////////////////
/// Load a library file in cling's memory.
/// if 'system' is true, the library is never unloaded.
/// Return 0 on success, -1 on failure.

Int_t TCling::Load(const char* filename, Bool_t system)
{
   if (!fAllowLibLoad) {
      Error("Load","Trying to load library (%s) from rootcling.",filename);
      return -1;
   }

   // Used to return 0 on success, 1 on duplicate, -1 on failure, -2 on "fatal".
   R__LOCKGUARD_CLING(gInterpreterMutex);
   cling::DynamicLibraryManager* DLM = fInterpreter->getDynamicLibraryManager();
   std::string canonLib = DLM->lookupLibrary(filename);
   cling::DynamicLibraryManager::LoadLibResult res
      = cling::DynamicLibraryManager::kLoadLibNotFound;
   if (!canonLib.empty()) {
      if (system)
         res = DLM->loadLibrary(filename, system);
      else {
         // For the non system libs, we'd like to be able to unload them.
         // FIXME: Here we lose the information about kLoadLibAlreadyLoaded case.
         cling::Interpreter::CompilationResult compRes;
         HandleInterpreterException(fMetaProcessor, Form(".L %s", canonLib.c_str()), compRes, /*cling::Value*/0);
         if (compRes == cling::Interpreter::kSuccess)
            res = cling::DynamicLibraryManager::kLoadLibSuccess;
      }
   }

   if (res == cling::DynamicLibraryManager::kLoadLibSuccess) {
      UpdateListOfLoadedSharedLibraries();
   }
   switch (res) {
   case cling::DynamicLibraryManager::kLoadLibSuccess: return 0;
   case cling::DynamicLibraryManager::kLoadLibAlreadyLoaded:  return 1;
   default: break;
   };
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
/// Load a macro file in cling's memory.

void TCling::LoadMacro(const char* filename, EErrorCode* error)
{
   ProcessLine(Form(".L %s", filename), error);
}

////////////////////////////////////////////////////////////////////////////////
/// Let cling process a command line asynch.

Long_t TCling::ProcessLineAsynch(const char* line, EErrorCode* error)
{
   return ProcessLine(line, error);
}

////////////////////////////////////////////////////////////////////////////////
/// Let cling process a command line synchronously, i.e we are waiting
/// it will be finished.

Long_t TCling::ProcessLineSynch(const char* line, EErrorCode* error)
{
   R__LOCKGUARD_CLING(fLockProcessLine ? gInterpreterMutex : 0);
   if (gApplication) {
      if (gApplication->IsCmdThread()) {
         return ProcessLine(line, error);
      }
      return 0;
   }
   return ProcessLine(line, error);
}

////////////////////////////////////////////////////////////////////////////////
/// Directly execute an executable statement (e.g. "func()", "3+5", etc.
/// however not declarations, like "Int_t x;").

Long_t TCling::Calc(const char* line, EErrorCode* error)
{
#ifdef R__WIN32
   // Test on ApplicationImp not being 0 is needed because only at end of
   // TApplication ctor the IsLineProcessing flag is set to 0, so before
   // we can not use it.
   if (gApplication && gApplication->GetApplicationImp()) {
      while (gROOT->IsLineProcessing() && !gApplication) {
         Warning("Calc", "waiting for cling thread to free");
         gSystem->Sleep(500);
      }
      gROOT->SetLineIsProcessing();
   }
#endif // R__WIN32
   R__LOCKGUARD_CLING(gInterpreterMutex);
   if (error) {
      *error = TInterpreter::kNoError;
   }
   cling::Value valRef;
   cling::Interpreter::CompilationResult cr = fInterpreter->evaluate(line, valRef);
   if (cr != cling::Interpreter::kSuccess) {
      // Failure in compilation.
      if (error) {
         // Note: Yes these codes are weird.
         *error = TInterpreter::kRecoverable;
      }
      return 0L;
   }
   if (!valRef.isValid()) {
      // Failure at runtime.
      if (error) {
         // Note: Yes these codes are weird.
         *error = TInterpreter::kDangerous;
      }
      return 0L;
   }

   if (valRef.isVoid()) {
      return 0;
   }

   RegisterTemporary(valRef);
#ifdef R__WIN32
   if (gApplication && gApplication->GetApplicationImp()) {
      gROOT->SetLineHasBeenProcessed();
   }
#endif // R__WIN32
   return valRef.simplisticCastAs<long>();
}

////////////////////////////////////////////////////////////////////////////////
/// Set a getline function to call when input is needed.

void TCling::SetGetline(const char * (*getlineFunc)(const char* prompt),
                                void (*histaddFunc)(const char* line))
{
   // If cling offers a replacement for G__pause(), it would need to
   // also offer a way to customize at least the history recording.

#if defined(R__MUST_REVISIT)
#if R__MUST_REVISIT(6,2)
   Warning("SetGetline","Cling should support the equivalent of SetGetlineFunc(getlineFunc, histaddFunc)");
#endif
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Helper function to increase the internal Cling count of transactions
/// that change the AST.

Bool_t TCling::HandleNewTransaction(const cling::Transaction &T)
{
   R__LOCKGUARD(gInterpreterMutex);

   if ((std::distance(T.decls_begin(), T.decls_end()) != 1)
      || T.deserialized_decls_begin() != T.deserialized_decls_end()
      || T.macros_begin() != T.macros_end()
      || ((!T.getFirstDecl().isNull()) && ((*T.getFirstDecl().begin()) != T.getWrapperFD()))) {
      fTransactionCount++;
      return true;
   }
   return false;
}

////////////////////////////////////////////////////////////////////////////////
/// Delete object from cling symbol table so it can not be used anymore.
/// cling objects are always on the heap.

void TCling::RecursiveRemove(TObject* obj)
{
   // NOTE: When replacing the mutex by a ReadWrite mutex, we **must**
   // put in place the Read/Write part here.  Keeping the write lock
   // here is 'catasptrophic' for scaling as it means that ALL calls
   // to RecursiveRemove will take the write lock and performance
   // of many threads trying to access the write lock at the same
   // time is relatively bad.
   R__READ_LOCKGUARD(ROOT::gCoreMutex);
   // Note that fgSetOfSpecials is supposed to be updated by TClingCallbacks::tryFindROOTSpecialInternal
   // (but isn't at the moment).
   if (obj->IsOnHeap() && fgSetOfSpecials && !((std::set<TObject*>*)fgSetOfSpecials)->empty()) {
      std::set<TObject*>::iterator iSpecial = ((std::set<TObject*>*)fgSetOfSpecials)->find(obj);
      if (iSpecial != ((std::set<TObject*>*)fgSetOfSpecials)->end()) {
         R__WRITE_LOCKGUARD(ROOT::gCoreMutex);
         DeleteGlobal(obj);
         ((std::set<TObject*>*)fgSetOfSpecials)->erase(iSpecial);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Pressing Ctrl+C should forward here. In the case where we have had
/// continuation requested we must reset it.

void TCling::Reset()
{
   fMetaProcessor->cancelContinuation();
   // Reset the Cling state to the state saved by the last call to
   // TCling::SaveContext().
#if defined(R__MUST_REVISIT)
#if R__MUST_REVISIT(6,2)
   R__LOCKGUARD(gInterpreterMutex);
   Warning("Reset","Cling should support the equivalent of scratch_upto(&fDictPos)");
#endif
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Reset the Cling state to its initial state.

void TCling::ResetAll()
{
#if defined(R__MUST_REVISIT)
#if R__MUST_REVISIT(6,2)
   R__LOCKGUARD(gInterpreterMutex);
   Warning("ResetAll","Cling should support the equivalent of complete reset (unload everything but the startup decls.");
#endif
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Reset in Cling the list of global variables to the state saved by the last
/// call to TCling::SaveGlobalsContext().
///
/// Note: Right now, all we do is run the global destructors.

void TCling::ResetGlobals()
{
   R__LOCKGUARD(gInterpreterMutex);
   // TODO:
   // Here we should iterate over the transactions (N-3) and revert.
   // N-3 because the first three internal to cling.

   fInterpreter->runAndRemoveStaticDestructors();
}

////////////////////////////////////////////////////////////////////////////////
/// Reset the Cling 'user' global objects/variables state to the state saved by the last
/// call to TCling::SaveGlobalsContext().

void TCling::ResetGlobalVar(void* obj)
{
#if defined(R__MUST_REVISIT)
#if R__MUST_REVISIT(6,2)
   R__LOCKGUARD(gInterpreterMutex);
   Warning("ResetGlobalVar","Cling should support the equivalent of resetglobalvar(obj)");
#endif
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Rewind Cling dictionary to the point where it was before executing
/// the current macro. This function is typically called after SEGV or
/// ctlr-C after doing a longjmp back to the prompt.

void TCling::RewindDictionary()
{
#if defined(R__MUST_REVISIT)
#if R__MUST_REVISIT(6,2)
   R__LOCKGUARD(gInterpreterMutex);
   Warning("RewindDictionary","Cling should provide a way to revert transaction similar to rewinddictionary()");
#endif
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Delete obj from Cling symbol table so it cannot be accessed anymore.
/// Returns 1 in case of success and 0 in case object was not in table.

Int_t TCling::DeleteGlobal(void* obj)
{
#if defined(R__MUST_REVISIT)
#if R__MUST_REVISIT(6,2)
   R__LOCKGUARD(gInterpreterMutex);
   Warning("DeleteGlobal","Cling should provide the equivalent of deleteglobal(obj), see also DeleteVariable.");
#endif
#endif
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Undeclare obj called name.
/// Returns 1 in case of success, 0 for failure.

Int_t TCling::DeleteVariable(const char* name)
{
#if defined(R__MUST_REVISIT)
#if R__MUST_REVISIT(6,2)
   Warning("DeleteVariable","should do more that just reseting the value to zero");
#endif
#endif

   R__LOCKGUARD(gInterpreterMutex);
   llvm::StringRef srName(name);
   const char* unscopedName = name;
   llvm::StringRef::size_type posScope = srName.rfind("::");
   const clang::DeclContext* declCtx = 0;
   if (posScope != llvm::StringRef::npos) {
      const cling::LookupHelper& lh = fInterpreter->getLookupHelper();
      const clang::Decl* scopeDecl
         = lh.findScope(srName.substr(0, posScope),
                        cling::LookupHelper::WithDiagnostics);
      if (!scopeDecl) {
         Error("DeleteVariable", "Cannot find enclosing scope for variable %s",
               name);
         return 0;
      }
      declCtx = llvm::dyn_cast<clang::DeclContext>(scopeDecl);
      if (!declCtx) {
         Error("DeleteVariable",
               "Enclosing scope for variable %s is not a declaration context",
               name);
         return 0;
      }
      unscopedName += posScope + 2;
   }
   // Could trigger deserialization of decls.
   cling::Interpreter::PushTransactionRAII RAII(fInterpreter);
   clang::NamedDecl* nVarDecl
      = cling::utils::Lookup::Named(&fInterpreter->getSema(), unscopedName, declCtx);
   if (!nVarDecl) {
      Error("DeleteVariable", "Unknown variable %s", name);
      return 0;
   }
   clang::VarDecl* varDecl = llvm::dyn_cast<clang::VarDecl>(nVarDecl);
   if (!varDecl) {
      Error("DeleteVariable", "Entity %s is not a variable", name);
      return 0;
   }

   clang::QualType qType = varDecl->getType();
   const clang::Type* type = qType->getUnqualifiedDesugaredType();
   // Cannot set a reference's address to nullptr; the JIT can place it
   // into read-only memory (ROOT-7100).
   if (type->isPointerType()) {
      int** ppInt = (int**)fInterpreter->getAddressOfGlobal(GlobalDecl(varDecl));
      // set pointer to invalid.
      if (ppInt) *ppInt = 0;
   }
   return 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Save the current Cling state.

void TCling::SaveContext()
{
#if defined(R__MUST_REVISIT)
#if R__MUST_REVISIT(6,2)
   R__LOCKGUARD(gInterpreterMutex);
   Warning("SaveContext","Cling should provide a way to record a state watermark similar to store_dictposition(&fDictPos)");
#endif
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Save the current Cling state of global objects.

void TCling::SaveGlobalsContext()
{
#if defined(R__MUST_REVISIT)
#if R__MUST_REVISIT(6,2)
   R__LOCKGUARD(gInterpreterMutex);
   Warning("SaveGlobalsContext","Cling should provide a way to record a watermark for the list of global variable similar to store_dictposition(&fDictPosGlobals)");
#endif
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// No op: see TClingCallbacks (used to update the list of globals)

void TCling::UpdateListOfGlobals()
{
}

////////////////////////////////////////////////////////////////////////////////
/// No op: see TClingCallbacks (used to update the list of global functions)

void TCling::UpdateListOfGlobalFunctions()
{
}

////////////////////////////////////////////////////////////////////////////////
/// No op: see TClingCallbacks (used to update the list of types)

void TCling::UpdateListOfTypes()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Check in what order the member of a tuple are layout.
enum class ETupleOrdering {
   kAscending,
   kDescending,
   kUnexpected
};

struct AlternateTupleIntDoubleAsc
{
   Int_t    _0;
   Double_t _1;
};

struct AlternateTupleIntDoubleDes
{
   Double_t _1;
   Int_t    _0;
};

static ETupleOrdering IsTupleAscending()
{
   std::tuple<int,double> value;
   AlternateTupleIntDoubleAsc asc;
   AlternateTupleIntDoubleDes des;

   size_t offset0 = ((char*)&(std::get<0>(value))) - ((char*)&value);
   size_t offset1 = ((char*)&(std::get<1>(value))) - ((char*)&value);

   size_t ascOffset0 = ((char*)&(asc._0)) - ((char*)&asc);
   size_t ascOffset1 = ((char*)&(asc._1)) - ((char*)&asc);

   size_t desOffset0 = ((char*)&(des._0)) - ((char*)&des);
   size_t desOffset1 = ((char*)&(des._1)) - ((char*)&des);

   if (offset0 == ascOffset0 && offset1 == ascOffset1) {
      return ETupleOrdering::kAscending;
   } else if (offset0 == desOffset0 && offset1 == desOffset1) {
      return ETupleOrdering::kDescending;
   } else {
      return ETupleOrdering::kUnexpected;
   }
}

std::string AlternateTuple(const char *classname)
{
   TClassEdit::TSplitType tupleContent(classname);
   std::string alternateName = "TEmulatedTuple";
   alternateName.append( classname + 5 );

   std::string guard_name;
   ROOT::TMetaUtils::GetCppName(guard_name,alternateName.c_str());
   std::ostringstream guard;
   guard << "ROOT_INTERNAL_TEmulated_";
   guard << guard_name;

   std::ostringstream alternateTuple;
   alternateTuple << "#ifndef " << guard.str() << "\n";
   alternateTuple << "#define " << guard.str() << "\n";
   alternateTuple << "namespace ROOT { namespace Internal {\n";
   alternateTuple << "template <class... Types> struct TEmulatedTuple;\n";
   alternateTuple << "template <> struct " << alternateName << " {\n";

   // This could also be a compile time choice ...
   switch(IsTupleAscending()) {
      case ETupleOrdering::kAscending: {
         unsigned int nMember = 0;
         auto iter = tupleContent.fElements.begin() + 1; // Skip the template name (tuple)
         auto theEnd = tupleContent.fElements.end() - 1; // skip the 'stars'.
         while (iter != theEnd) {
            alternateTuple << "   " << *iter << " _" << nMember << ";\n";
            ++iter;
            ++nMember;
         }
         break;
      }
      case ETupleOrdering::kDescending: {
         unsigned int nMember = tupleContent.fElements.size() - 3;
         auto iter = tupleContent.fElements.rbegin() + 1; // Skip the template name (tuple)
         auto theEnd = tupleContent.fElements.rend() - 1; // skip the 'stars'.
         while (iter != theEnd) {
            alternateTuple << "   " << *iter << " _" << nMember << ";\n";
            ++iter;
            --nMember;
         }
         break;
      }
      case ETupleOrdering::kUnexpected: {
         Fatal("TCling::SetClassInfo::AlternateTuple",
               "Layout of std::tuple on this platform is unexpected.");
         break;
      }
   }

   alternateTuple << "};\n";
   alternateTuple << "}}\n";
   alternateTuple << "#endif\n";
   if (!gCling->Declare(alternateTuple.str().c_str())) {
      Error("Load","Could not declare %s",alternateName.c_str());
      return "";
   }
   alternateName = "ROOT::Internal::" + alternateName;
   return alternateName;
}

////////////////////////////////////////////////////////////////////////////////
/// Set pointer to the TClingClassInfo in TClass.
/// If 'reload' is true, (attempt to) generate a new ClassInfo even if we
/// already have one.

void TCling::SetClassInfo(TClass* cl, Bool_t reload)
{
   R__LOCKGUARD(gInterpreterMutex);
   if (cl->fClassInfo && !reload) {
      return;
   }
   //Remove the decl_id from the DeclIdToTClass map
   TClingClassInfo* TClinginfo = (TClingClassInfo*) cl->fClassInfo;
   if (TClinginfo) {
      TClass::RemoveClassDeclId(TClinginfo->GetDeclId());
   }
   delete TClinginfo;
   cl->fClassInfo = 0;
   std::string name(cl->GetName());

   // Handle the special case of 'tuple' where we ignore the real implementation
   // details and just overlay a 'simpler'/'simplistic' version that is easy
   // for the I/O to understand and handle.
   if (!(fCxxModulesEnabled && IsFromRootCling()) && strncmp(cl->GetName(),"tuple<",strlen("tuple<"))==0) {

      name = AlternateTuple(cl->GetName());

   }

   TClingClassInfo* info = new TClingClassInfo(fInterpreter, name.c_str());
   if (!info->IsValid()) {
      if (cl->fState != TClass::kHasTClassInit) {
         if (cl->fStreamerInfo->GetEntries() != 0) {
            cl->fState = TClass::kEmulated;
         } else {
            cl->fState = TClass::kForwardDeclared;
         }
      }
      delete info;
      return;
   }
   cl->fClassInfo = (ClassInfo_t*)info; // Note: We are transferring ownership here.
   // In case a class contains an external enum, the enum will be seen as a
   // class. We must detect this special case and make the class a Zombie.
   // Here we assume that a class has at least one method.
   // We can NOT call TClass::Property from here, because this method
   // assumes that the TClass is well formed to do a lot of information
   // caching. The method SetClassInfo (i.e. here) is usually called during
   // the building phase of the TClass, hence it is NOT well formed yet.
   Bool_t zombieCandidate = kFALSE;
   if (
      info->IsValid() &&
      !(info->Property() & (kIsClass | kIsStruct | kIsNamespace))
   ) {
      zombieCandidate = kTRUE;
   }
   if (!info->IsLoaded()) {
      if (info->Property() & (kIsNamespace)) {
         // Namespaces can have info but no corresponding CINT dictionary
         // because they are auto-created if one of their contained
         // classes has a dictionary.
         zombieCandidate = kTRUE;
      }
      // this happens when no dictionary is available
      delete info;
      cl->fClassInfo = 0;
   }
   if (zombieCandidate && !cl->GetCollectionType()) {
      cl->MakeZombie();
   }
   // If we reach here, the info was valid (See early returns).
   if (cl->fState != TClass::kHasTClassInit) {
      if (cl->fClassInfo) {
         cl->fState = TClass::kInterpreted;
         cl->ResetBit(TClass::kIsEmulation);
      } else {
//         if (TClassEdit::IsSTLCont(cl->GetName()) {
//            There will be an emulated collection proxy, is that the same?
//            cl->fState = TClass::kEmulated;
//         } else {
         if (cl->fStreamerInfo->GetEntries() != 0) {
            cl->fState = TClass::kEmulated;
         } else {
            cl->fState = TClass::kForwardDeclared;
         }
//         }
      }
   }
   if (cl->fClassInfo) {
      TClass::AddClassToDeclIdMap(((TClingClassInfo*)cl->fClassInfo)->GetDeclId(), cl);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Checks if an entity with the specified name is defined in Cling.
/// Returns kUnknown if the entity is not defined.
/// Returns kWithClassDefInline if the entity exists and has a ClassDefInline
/// Returns kKnown if the entity is defined.
///
/// By default, structs, namespaces, classes, enums and unions are looked for.
/// If the flag isClassOrNamespaceOnly is true, classes, structs and
/// namespaces only are considered. I.e. if the name is an enum or a union,
/// the returned value is false.
///
/// In the case where the class is not loaded and belongs to a namespace
/// or is nested, looking for the full class name is outputting a lots of
/// (expected) error messages.  Currently the only way to avoid this is to
/// specifically check that each level of nesting is already loaded.
/// In case of templates the idea is that everything between the outer
/// '<' and '>' has to be skipped, e.g.: aap<pippo<noot>::klaas>::a_class

TInterpreter::ECheckClassInfo
TCling::CheckClassInfo(const char *name, Bool_t autoload, Bool_t isClassOrNamespaceOnly /* = kFALSE*/)
{
   R__LOCKGUARD(gInterpreterMutex);
   static const char *anonEnum = "anonymous enum ";
   static const int cmplen = strlen(anonEnum);

   if (0 == strncmp(name,anonEnum,cmplen)) {
      return kUnknown;
   }

   // Avoid the double search below in case the name is a fundamental type
   // or typedef to a fundamental type.
   THashTable *typeTable = dynamic_cast<THashTable*>( gROOT->GetListOfTypes() );
   TDataType *fundType = (TDataType *)typeTable->THashTable::FindObject( name );

   if (fundType && fundType->GetType() < TVirtualStreamerInfo::kObject
       && fundType->GetType() > 0) {
      // Fundamental type, no a class.
      return kUnknown;
   }

   // Migrated from within TClass::GetClass
   // If we want to know if a class or a namespace with this name exists in the
   // interpreter and this is an enum in the type system, before or after loading
   // according to the autoload function argument, return kUnknown.
   if (isClassOrNamespaceOnly && TEnum::GetEnum(name, autoload ? TEnum::kAutoload : TEnum::kNone))
      return kUnknown;

   const char *classname = name;

   int storeAutoload = SetClassAutoloading(autoload);

   // First we want to check whether the decl exist, but _without_
   // generating any template instantiation. However, the lookup
   // still will create a forward declaration of the class template instance
   // if it exist.  In this case, the return value of findScope will still
   // be zero but the type will be initialized.
   // Note in the corresponding code in ROOT 5, CINT was not instantiating
   // this forward declaration.
   const cling::LookupHelper& lh = fInterpreter->getLookupHelper();
   const clang::Type *type = 0;
   const clang::Decl *decl
      = lh.findScope(classname,
                     gDebug > 5 ? cling::LookupHelper::WithDiagnostics
                     : cling::LookupHelper::NoDiagnostics,
                     &type, /* intantiateTemplate= */ false );
   if (!decl) {
      std::string buf = TClassEdit::InsertStd(classname);
      decl = lh.findScope(buf,
                          gDebug > 5 ? cling::LookupHelper::WithDiagnostics
                          : cling::LookupHelper::NoDiagnostics,
                          &type,false);
   }

   if (type) {
      // If decl==0 and the type is valid, then we have a forward declaration.
      if (!decl) {
         // If we have a forward declaration for a class template instantiation,
         // we want to ignore it if it was produced/induced by the call to
         // findScope, however we can not distinguish those from the
         // instantiation induce by 'soft' use (and thus also induce by the
         // same underlying code paths)
         // ['soft' use = use not requiring a complete definition]
         // So to reduce the amount of disruption to the existing code we
         // would just ignore those for STL collection, for which we really
         // need to have the compiled collection proxy (and thus the TClass
         // bootstrap).
         clang::ClassTemplateSpecializationDecl *tmpltDecl =
            llvm::dyn_cast_or_null<clang::ClassTemplateSpecializationDecl>
               (type->getAsCXXRecordDecl());
         if (tmpltDecl && !tmpltDecl->getPointOfInstantiation().isValid()) {
            // Since the point of instantiation is invalid, we 'guess' that
            // the 'instantiation' of the forwarded type appended in
            // findscope.
            if (ROOT::TMetaUtils::IsSTLCont(*tmpltDecl)) {
               // For STL Collection we return kUnknown.
               SetClassAutoloading(storeAutoload);
               return kUnknown;
            }
         }
      }
      TClingClassInfo tci(fInterpreter, *type);
      if (!tci.IsValid()) {
         SetClassAutoloading(storeAutoload);
         return kUnknown;
      }
      auto propertiesMask = isClassOrNamespaceOnly ? kIsClass | kIsStruct | kIsNamespace :
                                                     kIsClass | kIsStruct | kIsNamespace | kIsEnum | kIsUnion;

      if (tci.Property() & propertiesMask) {
         bool hasClassDefInline = false;
         if (isClassOrNamespaceOnly) {
            // We do not need to check for ClassDefInline when this is called from
            // TClass::Init, we only do it for the call from TClass::GetClass.
            auto hasDictionary = tci.GetMethod("Dictionary", "", false, 0, ROOT::kExactMatch);
            auto implLineFunc = tci.GetMethod("ImplFileLine", "", false, 0, ROOT::kExactMatch);

            if (hasDictionary.IsValid() && implLineFunc.IsValid()) {
               int lineNumber = 0;
               bool success = false;
               std::tie(success, lineNumber) =
                  ROOT::TMetaUtils::GetTrivialIntegralReturnValue(implLineFunc.GetMethodDecl(), *fInterpreter);
               hasClassDefInline = success && (lineNumber == -1);
            }
         }

         // fprintf(stderr,"CheckClassInfo: %s had dict=%d  inline=%d\n",name,hasDictionary.IsValid()
         // , hasClassDefInline);

         // We are now sure that the entry is not in fact an autoload entry.
         SetClassAutoloading(storeAutoload);
         if (hasClassDefInline)
            return kWithClassDefInline;
         else
            return kKnown;
      } else {
         // We are now sure that the entry is not in fact an autoload entry.
         SetClassAutoloading(storeAutoload);
         return kUnknown;
      }
   }

   SetClassAutoloading(storeAutoload);
   if (decl)
      return kKnown;
   else
      return kUnknown;

   // Setting up iterator part of TClingTypedefInfo is too slow.
   // Copy the lookup code instead:
   /*
   TClingTypedefInfo t(fInterpreter, name);
   if (t.IsValid() && !(t.Property() & kIsFundamental)) {
      delete[] classname;
      SetClassAutoloading(storeAutoload);
      return kTRUE;
   }
   */

//   const clang::Decl *decl = lh.findScope(name);
//   if (!decl) {
//      std::string buf = TClassEdit::InsertStd(name);
//      decl = lh.findScope(buf);
//   }

//   SetClassAutoloading(storeAutoload);
//   return (decl);
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if there is a class template by the given name ...

Bool_t TCling::CheckClassTemplate(const char *name)
{
   const cling::LookupHelper& lh = fInterpreter->getLookupHelper();
   const clang::Decl *decl
      = lh.findClassTemplate(name,
                             gDebug > 5 ? cling::LookupHelper::WithDiagnostics
                             : cling::LookupHelper::NoDiagnostics);
   if (!decl) {
      std::string strname = "std::";
      strname += name;
      decl = lh.findClassTemplate(strname,
                                  gDebug > 5 ? cling::LookupHelper::WithDiagnostics
                                  : cling::LookupHelper::NoDiagnostics);
   }
   return 0 != decl;
}

////////////////////////////////////////////////////////////////////////////////
/// Create list of pointers to base class(es) for TClass cl.

void TCling::CreateListOfBaseClasses(TClass *cl) const
{
   R__LOCKGUARD(gInterpreterMutex);
   if (cl->fBase) {
      return;
   }
   TClingClassInfo *tci = (TClingClassInfo *)cl->GetClassInfo();
   if (!tci) return;
   TClingBaseClassInfo t(fInterpreter, tci);
   TList *listOfBase = new TList;
   while (t.Next()) {
      // if name cannot be obtained no use to put in list
      if (t.IsValid() && t.Name()) {
         TClingBaseClassInfo *a = new TClingBaseClassInfo(t);
         listOfBase->Add(new TBaseClass((BaseClassInfo_t *)a, cl));
      }
   }
   // Now that is complete, publish it.
   cl->fBase = listOfBase;
}

////////////////////////////////////////////////////////////////////////////////
/// Create list of pointers to enums for TClass cl.

void TCling::LoadEnums(TListOfEnums& enumList) const
{
   R__LOCKGUARD(gInterpreterMutex);

   const Decl * D;
   TClass* cl = enumList.GetClass();
   if (cl) {
      D = ((TClingClassInfo*)cl->GetClassInfo())->GetDecl();
   }
   else {
      D = fInterpreter->getCI()->getASTContext().getTranslationUnitDecl();
   }
   // Iterate on the decl of the class and get the enums.
   if (const clang::DeclContext* DC = dyn_cast<clang::DeclContext>(D)) {
      cling::Interpreter::PushTransactionRAII deserRAII(fInterpreter);
      // Collect all contexts of the namespace.
      llvm::SmallVector< DeclContext *, 4> allDeclContexts;
      const_cast< clang::DeclContext *>(DC)->collectAllContexts(allDeclContexts);
      for (llvm::SmallVector<DeclContext*, 4>::iterator declIter = allDeclContexts.begin(), declEnd = allDeclContexts.end();
           declIter != declEnd; ++declIter) {
         // Iterate on all decls for each context.
         for (clang::DeclContext::decl_iterator DI = (*declIter)->decls_begin(),
              DE = (*declIter)->decls_end(); DI != DE; ++DI) {
            if (const clang::EnumDecl* ED = dyn_cast<clang::EnumDecl>(*DI)) {
               // Get name of the enum type.
               std::string buf;
               PrintingPolicy Policy(ED->getASTContext().getPrintingPolicy());
               llvm::raw_string_ostream stream(buf);
               // Don't trigger fopen of the source file to count lines:
               Policy.AnonymousTagLocations = false;
               ED->getNameForDiagnostic(stream, Policy, /*Qualified=*/false);
               stream.flush();
               // If the enum is unnamed we do not add it to the list of enums i.e unusable.
               if (!buf.empty()) {
                  const char* name = buf.c_str();
                  // Add the enum to the list of loaded enums.
                  enumList.Get(ED, name);
               }
            }
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Create list of pointers to function templates for TClass cl.

void TCling::LoadFunctionTemplates(TClass* cl) const
{
   R__LOCKGUARD(gInterpreterMutex);

   const Decl * D;
   TListOfFunctionTemplates* funcTempList;
   if (cl) {
      D = ((TClingClassInfo*)cl->GetClassInfo())->GetDecl();
      funcTempList = (TListOfFunctionTemplates*)cl->GetListOfFunctionTemplates(false);
   }
   else {
      D = fInterpreter->getCI()->getASTContext().getTranslationUnitDecl();
      funcTempList = (TListOfFunctionTemplates*)gROOT->GetListOfFunctionTemplates();
   }
   // Iterate on the decl of the class and get the enums.
   if (const clang::DeclContext* DC = dyn_cast<clang::DeclContext>(D)) {
      cling::Interpreter::PushTransactionRAII deserRAII(fInterpreter);
      // Collect all contexts of the namespace.
      llvm::SmallVector< DeclContext *, 4> allDeclContexts;
      const_cast< clang::DeclContext *>(DC)->collectAllContexts(allDeclContexts);
      for (llvm::SmallVector<DeclContext*, 4>::iterator declIter = allDeclContexts.begin(),
           declEnd = allDeclContexts.end(); declIter != declEnd; ++declIter) {
         // Iterate on all decls for each context.
         for (clang::DeclContext::decl_iterator DI = (*declIter)->decls_begin(),
              DE = (*declIter)->decls_end(); DI != DE; ++DI) {
            if (const clang::FunctionTemplateDecl* FTD = dyn_cast<clang::FunctionTemplateDecl>(*DI)) {
                  funcTempList->Get(FTD);
            }
         }
      }
   }
}
////////////////////////////////////////////////////////////////////////////////
/// Create list of pointers to data members for TClass cl.
/// This is now a nop.  The creation and updating is handled in
/// TListOfDataMembers.

void TCling::CreateListOfDataMembers(TClass* cl) const
{
}

////////////////////////////////////////////////////////////////////////////////
/// Create list of pointers to methods for TClass cl.
/// This is now a nop.  The creation and updating is handled in
/// TListOfFunctions.

void TCling::CreateListOfMethods(TClass* cl) const
{
}

////////////////////////////////////////////////////////////////////////////////
/// Update the list of pointers to method for TClass cl
/// This is now a nop.  The creation and updating is handled in
/// TListOfFunctions.

void TCling::UpdateListOfMethods(TClass* cl) const
{
}

////////////////////////////////////////////////////////////////////////////////
/// Update the list of pointers to data members for TClass cl
/// This is now a nop.  The creation and updating is handled in
/// TListOfDataMembers.

void TCling::UpdateListOfDataMembers(TClass* cl) const
{
}

////////////////////////////////////////////////////////////////////////////////
/// Create list of pointers to method arguments for TMethod m.

void TCling::CreateListOfMethodArgs(TFunction* m) const
{
   R__LOCKGUARD(gInterpreterMutex);
   if (m->fMethodArgs) {
      return;
   }
   TList *arglist = new TList;
   TClingMethodArgInfo t(fInterpreter, (TClingMethodInfo*)m->fInfo);
   while (t.Next()) {
      if (t.IsValid()) {
         TClingMethodArgInfo* a = new TClingMethodArgInfo(t);
         arglist->Add(new TMethodArg((MethodArgInfo_t*)a, m));
      }
   }
   m->fMethodArgs = arglist;
}


////////////////////////////////////////////////////////////////////////////////
/// Generate a TClass for the given class.
/// Since the caller has already check the ClassInfo, let it give use the
/// result (via the value of emulation) rather than recalculate it.

TClass *TCling::GenerateTClass(const char *classname, Bool_t emulation, Bool_t silent /* = kFALSE */)
{
// For now the following line would lead to the (unwanted) instantiation
// of class template.  This could/would need to be resurrected only if
// we re-introduce so sort of automatic instantiation.   However this would
// have to include carefull look at the template parameter to avoid
// creating instance we can not really use (if the parameter are only forward
// declaration or do not have all the necessary interfaces).

   //   TClingClassInfo tci(fInterpreter, classname);
   //   if (1 || !tci.IsValid()) {

   Version_t version = 1;
   if (TClassEdit::IsSTLCont(classname)) {
      version = TClass::GetClass("TVirtualStreamerInfo")->GetClassVersion();
   }
   TClass *cl = new TClass(classname, version, silent);
   if (emulation) {
      cl->SetBit(TClass::kIsEmulation);
   } else {
      // Set the class version if the class is versioned.
      // Note that we cannot just call CLASS::Class_Version() as we might not have
      // an execution engine (when invoked from rootcling).

      // Do not call cl->GetClassVersion(), it has side effects!
      Version_t oldvers = cl->fClassVersion;
      if (oldvers == version && cl->GetClassInfo()) {
         // We have a version and it might need an update.
         Version_t newvers = oldvers;
         TClingClassInfo* cli = (TClingClassInfo*)cl->GetClassInfo();
         if (llvm::isa<clang::NamespaceDecl>(cli->GetDecl())) {
            // Namespaces don't have class versions.
            return cl;
         }
         TClingMethodInfo mi = cli->GetMethod("Class_Version", "", 0 /*poffset*/,
                                              ROOT::kExactMatch,
                                              TClingClassInfo::kInThisScope);
         if (!mi.IsValid()) {
            if (cl->TestBit(TClass::kIsTObject)) {
               Error("GenerateTClass",
                     "Cannot find %s::Class_Version()! Class version might be wrong.",
                     cl->GetName());
            }
            return cl;
         }
         newvers = ROOT::TMetaUtils::GetClassVersion(llvm::dyn_cast<clang::RecordDecl>(cli->GetDecl()),
                                                     *fInterpreter);
         if (newvers == -1) {
            // Didn't manage to determine the class version from the AST.
            // Use runtime instead.
            if ((mi.Property() & kIsStatic)
                && !fInterpreter->isInSyntaxOnlyMode()) {
               // This better be a static function.
               TClingCallFunc callfunc(fInterpreter, *fNormalizedCtxt);
               callfunc.SetFunc(&mi);
               newvers = callfunc.ExecInt(0);
            } else {
               Error("GenerateTClass",
                     "Cannot invoke %s::Class_Version()! Class version might be wrong.",
                     cl->GetName());
            }
         }
         if (newvers != oldvers) {
            cl->fClassVersion = newvers;
            cl->fStreamerInfo->Expand(newvers + 2 + 10);
         }
      }
   }

   return cl;

//   } else {
//      return GenerateTClass(&tci,silent);
//   }
}

#if 0
////////////////////////////////////////////////////////////////////////////////

static void GenerateTClass_GatherInnerIncludes(cling::Interpreter *interp, TString &includes,TClingClassInfo *info)
{
   includes += info->FileName();

   const clang::ClassTemplateSpecializationDecl *templateCl
      = llvm::dyn_cast<clang::ClassTemplateSpecializationDecl>(info->GetDecl());
   if (templateCl) {
      for(unsigned int i=0; i <  templateCl->getTemplateArgs().size(); ++i) {
          const clang::TemplateArgument &arg( templateCl->getTemplateArgs().get(i) );
          if (arg.getKind() == clang::TemplateArgument::Type) {
             const clang::Type *uType = ROOT::TMetaUtils::GetUnderlyingType( arg.getAsType() );

            if (!uType->isFundamentalType() && !uType->isEnumeralType()) {
               // We really need a header file.
               const clang::CXXRecordDecl *argdecl = uType->getAsCXXRecordDecl();
               if (argdecl) {
                  includes += ";";
                  TClingClassInfo subinfo(interp,*(argdecl->getASTContext().getRecordType(argdecl).getTypePtr()));
                  GenerateTClass_GatherInnerIncludes(interp, includes, &subinfo);
               } else {
                  std::string Result;
                  llvm::raw_string_ostream OS(Result);
                  arg.print(argdecl->getASTContext().getPrintingPolicy(),OS);
                  Warning("TCling::GenerateTClass","Missing header file for %s",OS.str().c_str());
               }
            }
          }
      }
   }
}
#endif

////////////////////////////////////////////////////////////////////////////////
/// Generate a TClass for the given class.

TClass *TCling::GenerateTClass(ClassInfo_t *classinfo, Bool_t silent /* = kFALSE */)
{
   TClingClassInfo *info = (TClingClassInfo*)classinfo;
   if (!info || !info->IsValid()) {
      Fatal("GenerateTClass","Requires a valid ClassInfo object");
      return 0;
   }
   // We are in the case where we have AST nodes for this class.
   TClass *cl = 0;
   std::string classname;
   info->FullName(classname,*fNormalizedCtxt); // Could we use Name()?
   if (TClassEdit::IsSTLCont(classname)) {
#if 0
      Info("GenerateTClass","Will (try to) generate the compiled TClass for %s.",classname.c_str());
      // We need to build up the list of required headers, by
      // looking at each template arguments.
      TString includes;
      GenerateTClass_GatherInnerIncludes(fInterpreter,includes,info);

      if (0 == GenerateDictionary(classname.c_str(),includes)) {
         // 0 means success.
         cl = TClass::LoadClass(classnam.c_str(), silent);
         if (cl == 0) {
            Error("GenerateTClass","Even though the dictionary generation for %s seemed successful we can't find the TClass bootstrap!",classname.c_str());
         }
      }
#endif
      if (cl == 0) {
         int version = TClass::GetClass("TVirtualStreamerInfo")->GetClassVersion();
         cl = new TClass(classinfo, version, 0, 0, -1, -1, silent);
         cl->SetBit(TClass::kIsEmulation);
      }
   } else {
      // For regular class, just create a TClass on the fly ...
      // Not quite useful yet, but that what CINT used to do anyway.
      cl = new TClass(classinfo, 1, 0, 0, -1, -1, silent);
   }
   // Add the new TClass to the map of declid and TClass*.
   if (cl) {
      TClass::AddClassToDeclIdMap(((TClingClassInfo*)classinfo)->GetDeclId(), cl);
   }
   return cl;
}

////////////////////////////////////////////////////////////////////////////////
/// Generate the dictionary for the C++ classes listed in the first
/// argument (in a semi-colon separated list).
/// 'includes' contains a semi-colon separated list of file to
/// #include in the dictionary.
/// For example:
/// ~~~ {.cpp}
///    gInterpreter->GenerateDictionary("vector<vector<float> >;list<vector<float> >","list;vector");
/// ~~~
/// or
/// ~~~ {.cpp}
///    gInterpreter->GenerateDictionary("myclass","myclass.h;myhelper.h");
/// ~~~

Int_t TCling::GenerateDictionary(const char* classes, const char* includes /* = "" */, const char* /* options  = 0 */)
{
   if (classes == 0 || classes[0] == 0) {
      Error("TCling::GenerateDictionary", "Cannot generate dictionary without passing classes.");
      return 0;
   }
   // Split the input list
   std::vector<std::string> listClasses;
   for (
      const char* current = classes, *prev = classes;
      *current != 0;
      ++current
   ) {
      if (*current == ';') {
         listClasses.push_back(std::string(prev, current - prev));
         prev = current + 1;
      }
      else if (*(current + 1) == 0) {
         listClasses.push_back(std::string(prev, current + 1 - prev));
         prev = current + 1;
      }
   }
   std::vector<std::string> listIncludes;
   if (!includes)
      includes = "";
   for (
      const char* current = includes, *prev = includes;
      *current != 0;
      ++current
   ) {
      if (*current == ';') {
         listIncludes.push_back(std::string(prev, current - prev));
         prev = current + 1;
      }
      else if (*(current + 1) == 0) {
         listIncludes.push_back(std::string(prev, current + 1 - prev));
         prev = current + 1;
      }
   }
   // Generate the temporary dictionary file
   return !TCling_GenerateDictionary(listClasses, listIncludes,
      std::vector<std::string>(), std::vector<std::string>());
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to cling Decl of global/static variable that is located
/// at the address given by addr.

TInterpreter::DeclId_t TCling::GetDataMember(ClassInfo_t *opaque_cl, const char *name) const
{
   R__LOCKGUARD(gInterpreterMutex);
   DeclId_t d;
   TClingClassInfo *cl = (TClingClassInfo*)opaque_cl;

   if (cl) {
      d = cl->GetDataMember(name);
      // We check if the decl of the data member has an annotation which indicates
      // an ioname.
      // In case this is true, if the name requested is not the ioname, we
      // return 0, as if the member did not exist. In some sense we override
      // the information in the TClassInfo instance, isolating the typesystem in
      // TClass from the one in the AST.
      if (const ValueDecl* decl = (const ValueDecl*) d){
         std::string ioName;
         bool hasIoName = ROOT::TMetaUtils::ExtractAttrPropertyFromName(*decl,"ioname",ioName);
         if (hasIoName && ioName != name) return 0;
      }
      return d;
   }
   // We are looking up for something on the TU scope.
   // FIXME: We do not want to go through TClingClassInfo(fInterpreter) because of redundant deserializations. That
   // interface will actually construct iterators and walk over the decls on the global scope. In would return the first
   // occurrence of a decl with the looked up name. However, that's not what C++ lookup would do: if we want to switch
   // to a more complete C++ lookup interface we need sift through the found names and pick up the declarations which
   // are only fulfilling ROOT's understanding for a Data Member.
   // FIXME: We should probably deprecate the TClingClassInfo(fInterpreter) interface and replace it withe something
   // similar as below.
   using namespace clang;
   Sema& SemaR = fInterpreter->getSema();
   DeclarationName DName = &SemaR.Context.Idents.get(name);

   LookupResult R(SemaR, DName, SourceLocation(), Sema::LookupOrdinaryName,
                  Sema::ForRedeclaration);

   // Could trigger deserialization of decls.
   cling::Interpreter::PushTransactionRAII RAII(fInterpreter);
   cling::utils::Lookup::Named(&SemaR, R);

   LookupResult::Filter F = R.makeFilter();
   // Filter the data-member looking decls.
   while (F.hasNext()) {
      NamedDecl *D = F.next();
      if (isa<VarDecl>(D) || isa<FieldDecl>(D) || isa<EnumConstantDecl>(D) ||
          isa<IndirectFieldDecl>(D))
         continue;
      F.erase();
   }
   F.done();

   if (R.isSingleResult())
      return R.getFoundDecl();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to cling Decl of global/static variable that is located
/// at the address given by addr.

TInterpreter::DeclId_t TCling::GetEnum(TClass *cl, const char *name) const
{
   R__LOCKGUARD(gInterpreterMutex);

   const clang::Decl* possibleEnum = 0;
   // FInd the context of the decl.
   if (cl) {
      TClingClassInfo *cci = (TClingClassInfo*)cl->GetClassInfo();
      if (cci) {
         const clang::DeclContext* dc = 0;
         if (const clang::Decl* D = cci->GetDecl()) {
            if (!(dc = dyn_cast<clang::NamespaceDecl>(D))) {
               dc = dyn_cast<clang::RecordDecl>(D);
            }
         }
         if (dc) {
            // If it is a data member enum.
            // Could trigger deserialization of decls.
            cling::Interpreter::PushTransactionRAII RAII(fInterpreter);
            possibleEnum = cling::utils::Lookup::Tag(&fInterpreter->getSema(), name, dc);
         } else {
            Error("TCling::GetEnum", "DeclContext not found for %s .\n", name);
         }
      }
   } else {
      // If it is a global enum.
      // Could trigger deserialization of decls.
      cling::Interpreter::PushTransactionRAII RAII(fInterpreter);
      possibleEnum = cling::utils::Lookup::Tag(&fInterpreter->getSema(), name);
   }
   if (possibleEnum && (possibleEnum != (clang::Decl*)-1)
       && isa<clang::EnumDecl>(possibleEnum)) {
      return possibleEnum;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to cling DeclId for a global value

TInterpreter::DeclId_t TCling::GetDeclId( const llvm::GlobalValue *gv ) const
{
   if (!gv) return 0;

   llvm::StringRef mangled_name = gv->getName();

   int err = 0;
   char* demangled_name_c = TClassEdit::DemangleName(mangled_name.str().c_str(), err);
   if (err) {
      if (err == -2) {
         // It might simply be an unmangled global name.
         DeclId_t d;
         TClingClassInfo gcl(fInterpreter);
         d = gcl.GetDataMember(mangled_name.str().c_str());
         return d;
      }
      return 0;
   }

   std::string scopename(demangled_name_c);
   free(demangled_name_c);

   //
   //  Separate out the class or namespace part of the
   //  function name.
   //
   std::string dataname;

   if (!strncmp(scopename.c_str(), "typeinfo for ", sizeof("typeinfo for ")-1)) {
      scopename.erase(0, sizeof("typeinfo for ")-1);
   } else if (!strncmp(scopename.c_str(), "vtable for ", sizeof("vtable for ")-1)) {
      scopename.erase(0, sizeof("vtable for ")-1);
   } else {
      // See if it is a function
      std::string::size_type pos = scopename.rfind('(');
      if (pos != std::string::npos) {
         return 0;
      }
      // Separate the scope and member name
      pos = scopename.rfind(':');
      if (pos != std::string::npos) {
         if ((pos != 0) && (scopename[pos-1] == ':')) {
            dataname = scopename.substr(pos+1);
            scopename.erase(pos-1);
         }
      } else {
         scopename.clear();
         dataname = scopename;
      }
   }
   //fprintf(stderr, "name: '%s'\n", name.c_str());
   // Now we have the class or namespace name, so do the lookup.


   DeclId_t d;
   if (scopename.size()) {
      TClingClassInfo cl(fInterpreter,scopename.c_str());
      d = cl.GetDataMember(dataname.c_str());
   }
   else {
      TClingClassInfo gcl(fInterpreter);
      d = gcl.GetDataMember(dataname.c_str());
   }
   return d;
}

////////////////////////////////////////////////////////////////////////////////
/// NOT IMPLEMENTED.

TInterpreter::DeclId_t TCling::GetDataMemberWithValue(const void *ptrvalue) const
{
   Error("GetDataMemberWithValue()", "not implemented");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to cling DeclId for a data member with a given name.

TInterpreter::DeclId_t TCling::GetDataMemberAtAddr(const void *addr) const
{
   // NOT IMPLEMENTED.
   Error("GetDataMemberAtAddr()", "not implemented");
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the cling mangled name for a method of a class with parameters
/// params (params is a string of actual arguments, not formal ones). If the
/// class is 0 the global function list will be searched.

TString TCling::GetMangledName(TClass* cl, const char* method,
                               const char* params, Bool_t objectIsConst /* = kFALSE */)
{
   R__LOCKGUARD(gInterpreterMutex);
   TClingCallFunc func(fInterpreter,*fNormalizedCtxt);
   if (cl) {
      Long_t offset;
      func.SetFunc((TClingClassInfo*)cl->GetClassInfo(), method, params, objectIsConst,
         &offset);
   }
   else {
      TClingClassInfo gcl(fInterpreter);
      Long_t offset;
      func.SetFunc(&gcl, method, params, &offset);
   }
   TClingMethodInfo* mi = (TClingMethodInfo*) func.FactoryMethod();
   if (!mi) return "";
   TString mangled_name( mi->GetMangledName() );
   delete mi;
   return mangled_name;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the cling mangled name for a method of a class with a certain
/// prototype, i.e. "char*,int,float". If the class is 0 the global function
/// list will be searched.

TString TCling::GetMangledNameWithPrototype(TClass* cl, const char* method,
                                            const char* proto, Bool_t objectIsConst /* = kFALSE */,
                                            EFunctionMatchMode mode /* = kConversionMatch */)
{
   R__LOCKGUARD(gInterpreterMutex);
   if (cl) {
      return ((TClingClassInfo*)cl->GetClassInfo())->
         GetMethod(method, proto, objectIsConst, 0 /*poffset*/, mode).GetMangledName();
   }
   TClingClassInfo gcl(fInterpreter);
   return gcl.GetMethod(method, proto, objectIsConst, 0 /*poffset*/, mode).GetMangledName();
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to cling interface function for a method of a class with
/// parameters params (params is a string of actual arguments, not formal
/// ones). If the class is 0 the global function list will be searched.

void* TCling::GetInterfaceMethod(TClass* cl, const char* method,
                                 const char* params, Bool_t objectIsConst /* = kFALSE */)
{
   R__LOCKGUARD(gInterpreterMutex);
   TClingCallFunc func(fInterpreter,*fNormalizedCtxt);
   if (cl) {
      Long_t offset;
      func.SetFunc((TClingClassInfo*)cl->GetClassInfo(), method, params, objectIsConst,
                   &offset);
   }
   else {
      TClingClassInfo gcl(fInterpreter);
      Long_t offset;
      func.SetFunc(&gcl, method, params, &offset);
   }
   return (void*) func.InterfaceMethod();
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to cling interface function for a method of a class with
/// a certain name.

TInterpreter::DeclId_t TCling::GetFunction(ClassInfo_t *opaque_cl, const char* method)
{
   R__LOCKGUARD(gInterpreterMutex);
   DeclId_t f;
   TClingClassInfo *cl = (TClingClassInfo*)opaque_cl;
   if (cl) {
      f = cl->GetMethod(method).GetDeclId();
   }
   else {
      TClingClassInfo gcl(fInterpreter);
      f = gcl.GetMethod(method).GetDeclId();
   }
   return f;

}

////////////////////////////////////////////////////////////////////////////////
/// Insert overloads of name in cl to res.

void TCling::GetFunctionOverloads(ClassInfo_t *cl, const char *funcname,
                                  std::vector<DeclId_t>& res) const
{
   clang::Sema& S = fInterpreter->getSema();
   clang::ASTContext& Ctx = S.Context;
   const clang::Decl* CtxDecl
      = cl ? (const clang::Decl*)((TClingClassInfo*)cl)->GetDeclId():
      Ctx.getTranslationUnitDecl();
   auto RecDecl = llvm::dyn_cast<const clang::RecordDecl>(CtxDecl);
   const clang::DeclContext* DeclCtx = RecDecl;

   if (!DeclCtx)
      DeclCtx = dyn_cast<clang::NamespaceDecl>(CtxDecl);
   if (!DeclCtx) return;

   clang::DeclarationName DName;
   // The DeclarationName is funcname, unless it's a ctor or dtor.
   // FIXME: or operator or conversion! See enum clang::DeclarationName::NameKind.

   if (RecDecl) {
      if (RecDecl->getNameAsString() == funcname) {
         clang::QualType QT = Ctx.getTypeDeclType(RecDecl);
         DName = Ctx.DeclarationNames.getCXXConstructorName(Ctx.getCanonicalType(QT));
      } else if (funcname[0] == '~' && RecDecl->getNameAsString() == funcname + 1) {
         clang::QualType QT = Ctx.getTypeDeclType(RecDecl);
         DName = Ctx.DeclarationNames.getCXXDestructorName(Ctx.getCanonicalType(QT));
      } else {
         DName = &Ctx.Idents.get(funcname);
      }
   } else {
      DName = &Ctx.Idents.get(funcname);
   }

   clang::LookupResult R(S, DName, clang::SourceLocation(),
                         Sema::LookupOrdinaryName, clang::Sema::ForRedeclaration);
   S.LookupQualifiedName(R, const_cast<DeclContext*>(DeclCtx));
   if (R.empty()) return;
   R.resolveKind();
   res.reserve(res.size() + (R.end() - R.begin()));
   for (clang::LookupResult::iterator IR = R.begin(), ER = R.end();
        IR != ER; ++IR) {
      if (const clang::FunctionDecl* FD
          = llvm::dyn_cast<const clang::FunctionDecl>(*IR)) {
         if (!FD->getDescribedFunctionTemplate()) {
            res.push_back(FD);
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to cling interface function for a method of a class with
/// a certain prototype, i.e. "char*,int,float". If the class is 0 the global
/// function list will be searched.

void* TCling::GetInterfaceMethodWithPrototype(TClass* cl, const char* method,
                                              const char* proto,
                                              Bool_t objectIsConst /* = kFALSE */,
                                              EFunctionMatchMode mode /* = kConversionMatch */)
{
   R__LOCKGUARD(gInterpreterMutex);
   void* f;
   if (cl) {
      f = ((TClingClassInfo*)cl->GetClassInfo())->
         GetMethod(method, proto, objectIsConst, 0 /*poffset*/, mode).InterfaceMethod(*fNormalizedCtxt);
   }
   else {
      TClingClassInfo gcl(fInterpreter);
      f = gcl.GetMethod(method, proto, objectIsConst, 0 /*poffset*/, mode).InterfaceMethod(*fNormalizedCtxt);
   }
   return f;
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to cling DeclId for a method of a class with
/// a certain prototype, i.e. "char*,int,float". If the class is 0 the global
/// function list will be searched.

TInterpreter::DeclId_t TCling::GetFunctionWithValues(ClassInfo_t *opaque_cl, const char* method,
                                                     const char* params,
                                                     Bool_t objectIsConst /* = kFALSE */)
{
   R__LOCKGUARD(gInterpreterMutex);
   DeclId_t f;
   TClingClassInfo *cl = (TClingClassInfo*)opaque_cl;
   if (cl) {
      f = cl->GetMethodWithArgs(method, params, objectIsConst, 0 /*poffset*/).GetDeclId();
   }
   else {
      TClingClassInfo gcl(fInterpreter);
      f = gcl.GetMethod(method, params, objectIsConst, 0 /*poffset*/).GetDeclId();
   }
   return f;
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to cling interface function for a method of a class with
/// a certain prototype, i.e. "char*,int,float". If the class is 0 the global
/// function list will be searched.

TInterpreter::DeclId_t TCling::GetFunctionWithPrototype(ClassInfo_t *opaque_cl, const char* method,
                                                        const char* proto,
                                                        Bool_t objectIsConst /* = kFALSE */,
                                                        EFunctionMatchMode mode /* = kConversionMatch */)
{
   R__LOCKGUARD(gInterpreterMutex);
   DeclId_t f;
   TClingClassInfo *cl = (TClingClassInfo*)opaque_cl;
   if (cl) {
      f = cl->GetMethod(method, proto, objectIsConst, 0 /*poffset*/, mode).GetDeclId();
   }
   else {
      TClingClassInfo gcl(fInterpreter);
      f = gcl.GetMethod(method, proto, objectIsConst, 0 /*poffset*/, mode).GetDeclId();
   }
   return f;
}

////////////////////////////////////////////////////////////////////////////////
/// Return pointer to cling interface function for a method of a class with
/// a certain name.

TInterpreter::DeclId_t TCling::GetFunctionTemplate(ClassInfo_t *opaque_cl, const char* name)
{
   R__LOCKGUARD(gInterpreterMutex);
   DeclId_t f;
   TClingClassInfo *cl = (TClingClassInfo*)opaque_cl;
   if (cl) {
      f = cl->GetFunctionTemplate(name);
   }
   else {
      TClingClassInfo gcl(fInterpreter);
      f = gcl.GetFunctionTemplate(name);
   }
   return f;

}

////////////////////////////////////////////////////////////////////////////////
/// The 'name' is known to the interpreter, this function returns
/// the internal version of this name (usually just resolving typedefs)
/// This is used in particular to synchronize between the name used
/// by rootcling and by the run-time environment (TClass)
/// Return 0 if the name is not known.

void TCling::GetInterpreterTypeName(const char* name, std::string &output, Bool_t full)
{
   output.clear();

   R__LOCKGUARD(gInterpreterMutex);

   TClingClassInfo cl(fInterpreter, name);
   if (!cl.IsValid()) {
      return ;
   }
   if (full) {
      cl.FullName(output,*fNormalizedCtxt);
      return;
   }
   // Well well well, for backward compatibility we need to act a bit too
   // much like CINT.
   TClassEdit::TSplitType splitname( cl.Name(), TClassEdit::kDropStd );
   splitname.ShortType(output, TClassEdit::kDropStd );

   return;
}

////////////////////////////////////////////////////////////////////////////////
/// Execute a global function with arguments params.
///
/// FIXME: The cint-based version of this code does not check if the
///        SetFunc() call works, and does not do any real checking
///        for errors from the Exec() call.  It did fetch the most
///        recent cint security error and return that in error, but
///        this does not really translate well to cling/clang.  We
///        should enhance these interfaces so that we can report
///        compilation and runtime errors properly.

void TCling::Execute(const char* function, const char* params, int* error)
{
   R__LOCKGUARD_CLING(gInterpreterMutex);
   if (error) {
      *error = TInterpreter::kNoError;
   }
   TClingClassInfo cl(fInterpreter);
   Long_t offset = 0L;
   TClingCallFunc func(fInterpreter,*fNormalizedCtxt);
   func.SetFunc(&cl, function, params, &offset);
   func.Exec(0);
}

////////////////////////////////////////////////////////////////////////////////
/// Execute a method from class cl with arguments params.
///
/// FIXME: The cint-based version of this code does not check if the
///        SetFunc() call works, and does not do any real checking
///        for errors from the Exec() call.  It did fetch the most
///        recent cint security error and return that in error, but
///        this does not really translate well to cling/clang.  We
///        should enhance these interfaces so that we can report
///        compilation and runtime errors properly.

void TCling::Execute(TObject* obj, TClass* cl, const char* method,
                     const char* params, Bool_t objectIsConst, int* error)
{
   R__LOCKGUARD_CLING(gInterpreterMutex);
   if (error) {
      *error = TInterpreter::kNoError;
   }
   // If the actual class of this object inherits 2nd (or more) from TObject,
   // 'obj' is unlikely to be the start of the object (as described by IsA()),
   // hence gInterpreter->Execute will improperly correct the offset.
   void* addr = cl->DynamicCast(TObject::Class(), obj, kFALSE);
   Long_t offset = 0L;
   TClingCallFunc func(fInterpreter,*fNormalizedCtxt);
   func.SetFunc((TClingClassInfo*)cl->GetClassInfo(), method, params, objectIsConst, &offset);
   void* address = (void*)((Long_t)addr + offset);
   func.Exec(address);
}

////////////////////////////////////////////////////////////////////////////////

void TCling::Execute(TObject* obj, TClass* cl, const char* method,
                    const char* params, int* error)
{
   Execute(obj,cl,method,params,false,error);
}

////////////////////////////////////////////////////////////////////////////////
/// Execute a method from class cl with the arguments in array params
/// (params[0] ... params[n] = array of TObjString parameters).
/// Convert the TObjArray array of TObjString parameters to a character
/// string of comma separated parameters.
/// The parameters of type 'char' are enclosed in double quotes and all
/// internal quotes are escaped.

void TCling::Execute(TObject* obj, TClass* cl, TMethod* method,
                     TObjArray* params, int* error)
{
   if (!method) {
      Error("Execute", "No method was defined");
      return;
   }
   TList* argList = method->GetListOfMethodArgs();
   // Check number of actual parameters against of expected formal ones

   Int_t nparms = argList->LastIndex() + 1;
   Int_t argc   = params ? params->GetEntries() : 0;

   if (argc > nparms) {
      Error("Execute","Too many parameters to call %s, got %d but expected at most %d.",method->GetName(),argc,nparms);
      return;
   }
   if (nparms != argc) {
     // Let's see if the 'missing' argument are all defaulted.
     // if nparms==0 then either we stopped earlier either argc is also zero and we can't reach here.
     assert(nparms > 0);

     TMethodArg *arg = (TMethodArg *) argList->At( 0 );
     if (arg && arg->GetDefault() && arg->GetDefault()[0]) {
        // There is a default value for the first missing
        // argument, so we are fine.
     } else {
        Int_t firstDefault = -1;
        for (Int_t i = 0; i < nparms; i ++) {
           arg = (TMethodArg *) argList->At( i );
           if (arg && arg->GetDefault() && arg->GetDefault()[0]) {
              firstDefault = i;
              break;
           }
        }
        if (firstDefault >= 0) {
           Error("Execute","Too few arguments to call %s, got only %d but expected at least %d and at most %d.",method->GetName(),argc,firstDefault,nparms);
        } else {
           Error("Execute","Too few arguments to call %s, got only %d but expected %d.",method->GetName(),argc,nparms);
        }
        return;
     }
   }

   const char* listpar = "";
   TString complete(10);
   if (params) {
      // Create a character string of parameters from TObjArray
      TIter next(params);
      for (Int_t i = 0; i < argc; i ++) {
         TMethodArg* arg = (TMethodArg*) argList->At(i);
         TClingTypeInfo type(fInterpreter, arg->GetFullTypeName());
         TObjString* nxtpar = (TObjString*) next();
         if (i) {
            complete += ',';
         }
         if (strstr(type.TrueName(*fNormalizedCtxt), "char")) {
            TString chpar('\"');
            chpar += (nxtpar->String()).ReplaceAll("\"", "\\\"");
            // At this point we have to check if string contains \\"
            // and apply some more sophisticated parser. Not implemented yet!
            complete += chpar;
            complete += '\"';
         }
         else {
            complete += nxtpar->String();
         }
      }
      listpar = complete.Data();
   }

   // And now execute it.
   R__LOCKGUARD_CLING(gInterpreterMutex);
   if (error) {
      *error = TInterpreter::kNoError;
   }
   // If the actual class of this object inherits 2nd (or more) from TObject,
   // 'obj' is unlikely to be the start of the object (as described by IsA()),
   // hence gInterpreter->Execute will improperly correct the offset.
   void* addr = cl->DynamicCast(TObject::Class(), obj, kFALSE);
   TClingCallFunc func(fInterpreter,*fNormalizedCtxt);
   TClingMethodInfo *minfo = (TClingMethodInfo*)method->fInfo;
   func.Init(*minfo);
   func.SetArgs(listpar);
   // Now calculate the 'this' pointer offset for the method
   // when starting from the class described by cl.
   const CXXMethodDecl * mdecl = dyn_cast<CXXMethodDecl>(minfo->GetMethodDecl());
   Long_t offset = ((TClingClassInfo*)cl->GetClassInfo())->GetOffset(mdecl);
   void* address = (void*)((Long_t)addr + offset);
   func.Exec(address);
}

////////////////////////////////////////////////////////////////////////////////

void TCling::ExecuteWithArgsAndReturn(TMethod* method, void* address,
                                      const void* args[] /*=0*/,
                                      int nargs /*=0*/,
                                      void* ret/*= 0*/) const
{
   if (!method) {
      Error("ExecuteWithArgsAndReturn", "No method was defined");
      return;
   }

   TClingMethodInfo* minfo = (TClingMethodInfo*) method->fInfo;
   TClingCallFunc func(*minfo,*fNormalizedCtxt);
   func.ExecWithArgsAndReturn(address, args, nargs, ret);
}

////////////////////////////////////////////////////////////////////////////////
/// Execute a cling macro.

Long_t TCling::ExecuteMacro(const char* filename, EErrorCode* error)
{
   R__LOCKGUARD_CLING(fLockProcessLine ? gInterpreterMutex : 0);
   fCurExecutingMacros.push_back(filename);
   Long_t result = TApplication::ExecuteFile(filename, (int*)error);
   fCurExecutingMacros.pop_back();
   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the file name of the current un-included interpreted file.
/// See the documentation for GetCurrentMacroName().

const char* TCling::GetTopLevelMacroName() const
{
   Warning("GetTopLevelMacroName", "Must change return type!");
   return fCurExecutingMacros.back();
}

////////////////////////////////////////////////////////////////////////////////
/// Return the file name of the currently interpreted file,
/// included or not. Example to illustrate the difference between
/// GetCurrentMacroName() and GetTopLevelMacroName():
/// ~~~ {.cpp}
///   void inclfile() {
///   std::cout << "In inclfile.C" << std::endl;
///   std::cout << "  TCling::GetCurrentMacroName() returns  " <<
///      TCling::GetCurrentMacroName() << std::endl;
///   std::cout << "  TCling::GetTopLevelMacroName() returns " <<
///      TCling::GetTopLevelMacroName() << std::endl;
///   }
/// ~~~
/// ~~~ {.cpp}
///   void mymacro() {
///   std::cout << "In mymacro.C" << std::endl;
///   std::cout << "  TCling::GetCurrentMacroName() returns  " <<
///      TCling::GetCurrentMacroName() << std::endl;
///   std::cout << "  TCling::GetTopLevelMacroName() returns " <<
///      TCling::GetTopLevelMacroName() << std::endl;
///   std::cout << "  Now calling inclfile..." << std::endl;
///   gInterpreter->ProcessLine(".x inclfile.C");;
///   }
/// ~~~
/// Running mymacro.C will print:
///
/// ~~~ {.cpp}
/// root [0] .x mymacro.C
/// ~~~
/// In mymacro.C
/// ~~~ {.cpp}
///   TCling::GetCurrentMacroName() returns  ./mymacro.C
///   TCling::GetTopLevelMacroName() returns ./mymacro.C
/// ~~~
///   Now calling inclfile...
/// In inclfile.h
/// ~~~ {.cpp}
///   TCling::GetCurrentMacroName() returns  inclfile.C
///   TCling::GetTopLevelMacroName() returns ./mymacro.C
/// ~~~

const char* TCling::GetCurrentMacroName() const
{
#if defined(R__MUST_REVISIT)
#if R__MUST_REVISIT(6,0)
   Warning("GetCurrentMacroName", "Must change return type!");
#endif
#endif
   return fCurExecutingMacros.back();
}

////////////////////////////////////////////////////////////////////////////////
/// Return the absolute type of typeDesc.
/// E.g.: typeDesc = "class TNamed**", returns "TNamed".
/// You need to use the result immediately before it is being overwritten.

const char* TCling::TypeName(const char* typeDesc)
{
   TTHREAD_TLS(char*) t = 0;
   TTHREAD_TLS(unsigned int) tlen = 0;

   unsigned int dlen = strlen(typeDesc);
   if (dlen > tlen) {
      delete[] t;
      t = new char[dlen + 1];
      tlen = dlen;
   }
   const char* s, *template_start;
   if (!strstr(typeDesc, "(*)(")) {
      s = strchr(typeDesc, ' ');
      template_start = strchr(typeDesc, '<');
      if (!strcmp(typeDesc, "long long")) {
         strlcpy(t, typeDesc, dlen + 1);
      }
      else if (!strncmp(typeDesc, "unsigned ", s + 1 - typeDesc)) {
         strlcpy(t, typeDesc, dlen + 1);
      }
      // s is the position of the second 'word' (if any)
      // except in the case of templates where there will be a space
      // just before any closing '>': eg.
      //    TObj<std::vector<UShort_t,__malloc_alloc_template<0> > >*
      else if (s && (template_start == 0 || (s < template_start))) {
         strlcpy(t, s + 1, dlen + 1);
      }
      else {
         strlcpy(t, typeDesc, dlen + 1);
      }
   }
   else {
      strlcpy(t, typeDesc, dlen + 1);
   }
   int l = strlen(t);
   while (l > 0 && (t[l - 1] == '*' || t[l - 1] == '&')) {
      t[--l] = 0;
   }
   return t;
}

static bool requiresRootMap(const char* rootmapfile, cling::Interpreter* interp)
{
   if (!rootmapfile || !*rootmapfile)
      return true;

   llvm::StringRef moduleName = llvm::sys::path::filename(rootmapfile);
   moduleName.consume_front("lib");
   moduleName.consume_back(".rootmap");

   Module *M = interp->getCI()->getPreprocessor().getHeaderSearchInfo().lookupModule(moduleName);

   return !(M && interp->getSema().isModuleVisible(M));
}

////////////////////////////////////////////////////////////////////////////////
/// Read and parse a rootmapfile in its new format, and return 0 in case of
/// success, -1 if the file has already been read, and -3 in case its format
/// is the old one (e.g. containing "Library.ClassName"), -4 in case of syntax
/// error.

int TCling::ReadRootmapFile(const char *rootmapfile, TUniqueString *uniqueString)
{
   // For "class ", "namespace ", "typedef ", "header ", "enum ", "var " respectively
   const std::map<char, unsigned int> keyLenMap = {{'c',6},{'n',10},{'t',8},{'h',7},{'e',5},{'v',4}};

   if (rootmapfile && *rootmapfile) {
      std::string rootmapfileNoBackslash(rootmapfile);
#ifdef _MSC_VER
      std::replace(rootmapfileNoBackslash.begin(), rootmapfileNoBackslash.end(), '\\', '/');
#endif
      // Add content of a specific rootmap file
      if (fRootmapFiles->FindObject(rootmapfileNoBackslash.c_str()))
         return -1;

      if (uniqueString)
         uniqueString->Append(std::string("\n#line 1 \"Forward declarations from ") + rootmapfileNoBackslash + "\"\n");

      std::ifstream file(rootmapfileNoBackslash);
      std::string line; line.reserve(200);
      std::string lib_name; line.reserve(100);
      bool newFormat=false;
      while (getline(file, line, '\n')) {
         if (!newFormat &&
             (line.compare(0, 8, "Library.") == 0 || line.compare(0, 8, "Declare.") == 0)) {
            file.close();
            return -3; // old format
         }
         newFormat=true;

         if (line.compare(0, 9, "{ decls }") == 0 && requiresRootMap(rootmapfile, fInterpreter)) {
            // forward declarations

            while (getline(file, line, '\n')) {
               if (line[0] == '[') break;
               if (!uniqueString) {
                  Error("ReadRootmapFile", "Cannot handle \"{ decls }\" sections in custom rootmap file %s",
                        rootmapfileNoBackslash.c_str());
                  return -4;
               }
               uniqueString->Append(line);
            }
         }
         const char firstChar=line[0];
         if (firstChar == '[') {
            // new section (library)
            auto brpos = line.find(']');
            if (brpos == string::npos) continue;
            lib_name = line.substr(1, brpos-1);
            size_t nspaces = 0;
            while( lib_name[nspaces] == ' ' ) ++nspaces;
            if (nspaces) lib_name.replace(0, nspaces, "");
            if (gDebug > 3) {
               TString lib_nameTstr(lib_name.c_str());
               TObjArray* tokens = lib_nameTstr.Tokenize(" ");
               const char* lib = ((TObjString *)tokens->At(0))->GetName();
               const char* wlib = gSystem->DynamicPathName(lib, kTRUE);
               if (wlib) {
                  Info("ReadRootmapFile", "new section for %s", lib_nameTstr.Data());
               }
               else {
                  Info("ReadRootmapFile", "section for %s (library does not exist)", lib_nameTstr.Data());
               }
               delete[] wlib;
               delete tokens;
            }
         }
         else {
            auto keyLenIt = keyLenMap.find(firstChar);
            if (keyLenIt == keyLenMap.end()) continue;
            unsigned int keyLen = keyLenIt->second;
            // Do not make a copy, just start after the key
            const char *keyname = line.c_str()+keyLen;
            if (gDebug > 6)
               Info("ReadRootmapFile", "class %s in %s", keyname, lib_name.c_str());
            TEnvRec* isThere = fMapfile->Lookup(keyname);
            if (isThere){
               if(lib_name != isThere->GetValue()){ // the same key for two different libs
                  if (firstChar == 'n') {
                     if (gDebug > 3)
                        Info("ReadRootmapFile", "namespace %s found in %s is already in %s",
                           keyname, lib_name.c_str(), isThere->GetValue());
                  } else if (firstChar == 'h'){ // it is a header: add the libname to the list of libs to be loaded.
                     lib_name+=" ";
                     lib_name+=isThere->GetValue();
                     fMapfile->SetValue(keyname, lib_name.c_str());
                  }
                  else if (!TClassEdit::IsSTLCont(keyname)) {
                     Warning("ReadRootmapFile", "%s %s found in %s is already in %s", line.substr(0, keyLen).c_str(),
                           keyname, lib_name.c_str(), isThere->GetValue());
                  }
               } else { // the same key for the same lib
                  if (gDebug > 3)
                        Info("ReadRootmapFile","Key %s was already defined for %s", keyname, lib_name.c_str());
               }

            } else {
               fMapfile->SetValue(keyname, lib_name.c_str());
            }
         }
      }
      file.close();
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a resource table and read the (possibly) three resource files, i.e
/// $ROOTSYS/etc/system<name> (or ROOTETCDIR/system<name>), $HOME/<name> and
/// ./<name>. ROOT always reads ".rootrc" (in TROOT::InitSystem()). You can
/// read additional user defined resource files by creating additional TEnv
/// objects. By setting the shell variable ROOTENV_NO_HOME=1 the reading of
/// the $HOME/<name> resource file will be skipped. This might be useful in
/// case the home directory resides on an automounted remote file system
/// and one wants to avoid the file system from being mounted.

void TCling::InitRootmapFile(const char *name)
{
   Bool_t ignore = fMapfile->IgnoreDuplicates(kFALSE);

   fMapfile->SetRcName(name);

   TString sname = "system";
   sname += name;
   char *s = gSystem->ConcatFileName(TROOT::GetEtcDir(), sname);

   Int_t ret = ReadRootmapFile(s);
   if (ret == -3) // old format
      fMapfile->ReadFile(s, kEnvGlobal);
   delete [] s;
   if (!gSystem->Getenv("ROOTENV_NO_HOME")) {
      s = gSystem->ConcatFileName(gSystem->HomeDirectory(), name);
      ret = ReadRootmapFile(s);
      if (ret == -3) // old format
         fMapfile->ReadFile(s, kEnvUser);
      delete [] s;
      if (strcmp(gSystem->HomeDirectory(), gSystem->WorkingDirectory())) {
         ret = ReadRootmapFile(name);
         if (ret == -3) // old format
            fMapfile->ReadFile(name, kEnvLocal);
      }
   } else {
      ret = ReadRootmapFile(name);
      if (ret == -3) // old format
         fMapfile->ReadFile(name, kEnvLocal);
   }
   fMapfile->IgnoreDuplicates(ignore);
}


namespace {
   using namespace clang;

   class ExtVisibleStorageAdder: public RecursiveASTVisitor<ExtVisibleStorageAdder>{
      // This class is to be considered an helper for autoloading.
      // It is a recursive visitor is used to inspect namespaces coming from
      // forward declarations in rootmaps and to set the external visible
      // storage flag for them.
   public:
      ExtVisibleStorageAdder(std::unordered_set<const NamespaceDecl*>& nsSet): fNSSet(nsSet) {};
      bool VisitNamespaceDecl(NamespaceDecl* nsDecl) {
         // We want to enable the external lookup for this namespace
         // because it may shadow the lookup of other names contained
         // in that namespace

         nsDecl->setHasExternalVisibleStorage();
         fNSSet.insert(nsDecl);
         return true;
      }
   private:
      std::unordered_set<const NamespaceDecl*>& fNSSet;

   };
}

////////////////////////////////////////////////////////////////////////////////
/// Load map between class and library. If rootmapfile is specified a
/// specific rootmap file can be added (typically used by ACLiC).
/// In case of error -1 is returned, 0 otherwise.
/// The interpreter uses this information to automatically load the shared
/// library for a class (autoload mechanism), see the AutoLoad() methods below.

Int_t TCling::LoadLibraryMap(const char* rootmapfile)
{
   R__LOCKGUARD(gInterpreterMutex);
   // open the [system].rootmap files
   if (!fMapfile) {
      fMapfile = new TEnv();
      fMapfile->IgnoreDuplicates(kTRUE);
//       fMapNamespaces = new THashTable();
//       fMapNamespaces->SetOwner();
      fRootmapFiles = new TObjArray;
      fRootmapFiles->SetOwner();
      InitRootmapFile(".rootmap");
   }
   bool needsRootMap = true;
   if (rootmapfile && *rootmapfile)
      needsRootMap = requiresRootMap(rootmapfile, fInterpreter);

   // Prepare a list of all forward declarations for cling
   // For some experiments it is easily as big as 500k characters. To be on the
   // safe side, we go for 1M.
   TUniqueString uniqueString(1048576);

   // Load all rootmap files in the dynamic load path ((DY)LD_LIBRARY_PATH, etc.).
   // A rootmap file must end with the string ".rootmap".
   TString ldpath = gSystem->GetDynamicPath();
   if (ldpath != fRootmapLoadPath) {
      fRootmapLoadPath = ldpath;
#ifdef WIN32
      TObjArray* paths = ldpath.Tokenize(";");
#else
      TObjArray* paths = ldpath.Tokenize(":");
#endif
      TString d;
      for (Int_t i = 0; i < paths->GetEntriesFast(); i++) {
         d = ((TObjString *)paths->At(i))->GetString();
         // check if directory already scanned
         Int_t skip = 0;
         for (Int_t j = 0; j < i; j++) {
            TString pd = ((TObjString *)paths->At(j))->GetString();
            if (pd == d) {
               skip++;
               break;
            }
         }
         if (!skip) {
            void* dirp = gSystem->OpenDirectory(d);
            if (dirp) {
               if (gDebug > 3) {
                  Info("LoadLibraryMap", "%s", d.Data());
               }
               const char* f1;
               while ((f1 = gSystem->GetDirEntry(dirp))) {
                  TString f = f1;
                  if (f.EndsWith(".rootmap")) {
                     TString p;
                     p = d + "/" + f;
                     if (!gSystem->AccessPathName(p, kReadPermission)) {
                        if (!fRootmapFiles->FindObject(f) && f != ".rootmap") {
                           if (gDebug > 4) {
                              Info("LoadLibraryMap", "   rootmap file: %s", p.Data());
                           }
                           Int_t ret;
                           if (needsRootMap)
                              ret = ReadRootmapFile(p,&uniqueString);
                           else
                              ret = ReadRootmapFile(p);

                           if (ret == 0)
                              fRootmapFiles->Add(new TNamed(gSystem->BaseName(f), p.Data()));
                           if (ret == -3) {
                              // old format
                              fMapfile->ReadFile(p, kEnvGlobal);
                              fRootmapFiles->Add(new TNamed(f, p));
                           }
                        }
                        // else {
                        //    fprintf(stderr,"Reject %s because %s is already there\n",p.Data(),f.Data());
                        //    fRootmapFiles->FindObject(f)->ls();
                        // }
                     }
                  }
                  if (f.BeginsWith("rootmap")) {
                     TString p;
                     p = d + "/" + f;
                     FileStat_t stat;
                     if (gSystem->GetPathInfo(p, stat) == 0 && R_ISREG(stat.fMode)) {
                        Warning("LoadLibraryMap", "please rename %s to end with \".rootmap\"", p.Data());
                     }
                  }
               }
            }
            gSystem->FreeDirectory(dirp);
         }
      }
      delete paths;
      if (fMapfile->GetTable() && !fMapfile->GetTable()->GetEntries()) {
         return -1;
      }
   }
   if (rootmapfile && *rootmapfile) {
      Int_t res;
      if (needsRootMap)
         res = ReadRootmapFile(rootmapfile, &uniqueString);
      else
         res = ReadRootmapFile(rootmapfile);
      if (res == 0) {
         //TString p = gSystem->ConcatFileName(gSystem->pwd(), rootmapfile);
         //fRootmapFiles->Add(new TNamed(gSystem->BaseName(rootmapfile), p.Data()));
         fRootmapFiles->Add(new TNamed(gSystem->BaseName(rootmapfile), rootmapfile));
      }
      else if (res == -3) {
         // old format
         Bool_t ignore = fMapfile->IgnoreDuplicates(kFALSE);
         fMapfile->ReadFile(rootmapfile, kEnvGlobal);
         fRootmapFiles->Add(new TNamed(gSystem->BaseName(rootmapfile), rootmapfile));
         fMapfile->IgnoreDuplicates(ignore);
      }
   }
   TEnvRec* rec;
   TIter next(fMapfile->GetTable());
   while ((rec = (TEnvRec*) next())) {
      TString cls = rec->GetName();
      if (!strncmp(cls.Data(), "Library.", 8) && cls.Length() > 8) {
         // get the first lib from the list of lib and dependent libs
         TString libs = rec->GetValue();
         if (libs == "") {
            continue;
         }
         TString delim(" ");
         TObjArray* tokens = libs.Tokenize(delim);
         const char* lib = ((TObjString*)tokens->At(0))->GetName();
         // convert "@@" to "::", we used "@@" because TEnv
         // considers "::" a terminator
         cls.Remove(0, 8);
         cls.ReplaceAll("@@", "::");
         // convert "-" to " ", since class names may have
         // blanks and TEnv considers a blank a terminator
         cls.ReplaceAll("-", " ");
         if (gDebug > 6) {
            const char* wlib = gSystem->DynamicPathName(lib, kTRUE);
            if (wlib) {
               Info("LoadLibraryMap", "class %s in %s", cls.Data(), wlib);
            }
            else {
               Info("LoadLibraryMap", "class %s in %s (library does not exist)", cls.Data(), lib);
            }
            delete[] wlib;
         }
         // Fill in the namespace candidate list
//          Ssiz_t last = cls.Last(':');
//          if (last != kNPOS) {
//             // Please note that the funny op overload does substring.
//             TString namespaceCand = cls(0, last - 1);
//             // This is a reference to a substring that lives in fMapfile
//             if (!fMapNamespaces->FindObject(namespaceCand.Data()))
//                fMapNamespaces->Add(new TNamed(namespaceCand.Data(), ""));
//          }
         delete tokens;
      }
      else if (!strncmp(cls.Data(), "Declare.", 8) && cls.Length() > 8) {
         cls.Remove(0, 8);
         // convert "-" to " ", since class names may have
         // blanks and TEnv considers a blank a terminator
         cls.ReplaceAll("-", " ");
         fInterpreter->declare(cls.Data());
      }
   }

   // Process the forward declarations collected
   cling::Transaction* T = nullptr;
   auto compRes= fInterpreter->declare(uniqueString.Data(), &T);
   assert(cling::Interpreter::kSuccess == compRes && "A declaration in a rootmap could not be compiled");

   if (compRes!=cling::Interpreter::kSuccess){
      Warning("LoadLibraryMap",
               "Problems in %s declaring '%s' were encountered.", rootmapfile, uniqueString.Data()) ;
   }

   if (T){
      ExtVisibleStorageAdder evsAdder(fNSFromRootmaps);
      for (auto declIt = T->decls_begin(); declIt < T->decls_end(); ++declIt) {
         if (declIt->m_DGR.isSingleDecl()) {
            if (Decl* D = declIt->m_DGR.getSingleDecl()) {
               if (NamespaceDecl* NSD = dyn_cast<NamespaceDecl>(D)) {
                  evsAdder.TraverseDecl(NSD);
               }
            }
         }
      }
   }

   // clear duplicates

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Scan again along the dynamic path for library maps. Entries for the loaded
/// shared libraries are unloaded first. This can be useful after reseting
/// the dynamic path through TSystem::SetDynamicPath()
/// In case of error -1 is returned, 0 otherwise.

Int_t TCling::RescanLibraryMap()
{
   UnloadAllSharedLibraryMaps();
   LoadLibraryMap();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Reload the library map entries coming from all the loaded shared libraries,
/// after first unloading the current ones.
/// In case of error -1 is returned, 0 otherwise.

Int_t TCling::ReloadAllSharedLibraryMaps()
{
   const TString sharedLibLStr = GetSharedLibs();
   const TObjArray* sharedLibL = sharedLibLStr.Tokenize(" ");
   const Int_t nrSharedLibs = sharedLibL->GetEntriesFast();
   for (Int_t ilib = 0; ilib < nrSharedLibs; ilib++) {
      const TString sharedLibStr = ((TObjString*)sharedLibL->At(ilib))->GetString();
      const  TString sharedLibBaseStr = gSystem->BaseName(sharedLibStr);
      const Int_t ret = UnloadLibraryMap(sharedLibBaseStr);
      if (ret < 0) {
         continue;
      }
      TString rootMapBaseStr = sharedLibBaseStr;
      if (sharedLibBaseStr.EndsWith(".dll")) {
         rootMapBaseStr.ReplaceAll(".dll", "");
      }
      else if (sharedLibBaseStr.EndsWith(".DLL")) {
         rootMapBaseStr.ReplaceAll(".DLL", "");
      }
      else if (sharedLibBaseStr.EndsWith(".so")) {
         rootMapBaseStr.ReplaceAll(".so", "");
      }
      else if (sharedLibBaseStr.EndsWith(".sl")) {
         rootMapBaseStr.ReplaceAll(".sl", "");
      }
      else if (sharedLibBaseStr.EndsWith(".dl")) {
         rootMapBaseStr.ReplaceAll(".dl", "");
      }
      else if (sharedLibBaseStr.EndsWith(".a")) {
         rootMapBaseStr.ReplaceAll(".a", "");
      }
      else {
         Error("ReloadAllSharedLibraryMaps", "Unknown library type %s", sharedLibBaseStr.Data());
         delete sharedLibL;
         return -1;
      }
      rootMapBaseStr += ".rootmap";
      const char* rootMap = gSystem->Which(gSystem->GetDynamicPath(), rootMapBaseStr);
      if (!rootMap) {
         Error("ReloadAllSharedLibraryMaps", "Could not find rootmap %s in path", rootMap);
         delete[] rootMap;
         delete sharedLibL;
         return -1;
      }
      const Int_t status = LoadLibraryMap(rootMap);
      if (status < 0) {
         Error("ReloadAllSharedLibraryMaps", "Error loading map %s", rootMap);
         delete[] rootMap;
         delete sharedLibL;
         return -1;
      }
      delete[] rootMap;
   }
   delete sharedLibL;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Unload the library map entries coming from all the loaded shared libraries.
/// Returns 0 if succesful

Int_t TCling::UnloadAllSharedLibraryMaps()
{
   const TString sharedLibLStr = GetSharedLibs();
   const TObjArray* sharedLibL = sharedLibLStr.Tokenize(" ");
   for (Int_t ilib = 0; ilib < sharedLibL->GetEntriesFast(); ilib++) {
      const TString sharedLibStr = ((TObjString*)sharedLibL->At(ilib))->GetString();
      const  TString sharedLibBaseStr = gSystem->BaseName(sharedLibStr);
      UnloadLibraryMap(sharedLibBaseStr);
   }
   delete sharedLibL;
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Unload library map entries coming from the specified library.
/// Returns -1 in case no entries for the specified library were found,
/// 0 otherwise.

Int_t TCling::UnloadLibraryMap(const char* library)
{
   if (!fMapfile || !library || !*library) {
      return 0;
   }
   TString libname(library);
   Ssiz_t idx = libname.Last('.');
   if (idx != kNPOS) {
      libname.Remove(idx);
   }
   size_t len = libname.Length();
   TEnvRec *rec;
   TIter next(fMapfile->GetTable());
   R__LOCKGUARD(gInterpreterMutex);
   Int_t ret = 0;
   while ((rec = (TEnvRec *) next())) {
      TString cls = rec->GetName();
      if (cls.Length() > 2) {
         // get the first lib from the list of lib and dependent libs
         TString libs = rec->GetValue();
         if (libs == "") {
            continue;
         }
         TString delim(" ");
         TObjArray* tokens = libs.Tokenize(delim);
         const char* lib = ((TObjString *)tokens->At(0))->GetName();
         if (!strncmp(cls.Data(), "Library.", 8) && cls.Length() > 8) {
            // convert "@@" to "::", we used "@@" because TEnv
            // considers "::" a terminator
            cls.Remove(0, 8);
            cls.ReplaceAll("@@", "::");
            // convert "-" to " ", since class names may have
            // blanks and TEnv considers a blank a terminator
            cls.ReplaceAll("-", " ");
         }
         if (!strncmp(lib, libname.Data(), len)) {
            if (fMapfile->GetTable()->Remove(rec) == 0) {
               Error("UnloadLibraryMap", "entry for <%s, %s> not found in library map table", cls.Data(), lib);
               ret = -1;
            }
         }
         delete tokens;
      }
   }
   if (ret >= 0) {
      TString library_rootmap(library);
      if (!library_rootmap.EndsWith(".rootmap"))
         library_rootmap.Append(".rootmap");
      TNamed* mfile = 0;
      while ((mfile = (TNamed *)fRootmapFiles->FindObject(library_rootmap))) {
         fRootmapFiles->Remove(mfile);
         delete mfile;
      }
      fRootmapFiles->Compress();
   }
   return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// Register the autoloading information for a class.
/// libs is a space separated list of libraries.

Int_t TCling::SetClassSharedLibs(const char *cls, const char *libs)
{
   if (!cls || !*cls)
      return 0;

   TString key = TString("Library.") + cls;
   // convert "::" to "@@", we used "@@" because TEnv
   // considers "::" a terminator
   key.ReplaceAll("::", "@@");
   // convert "-" to " ", since class names may have
   // blanks and TEnv considers a blank a terminator
   key.ReplaceAll(" ", "-");

   R__LOCKGUARD(gInterpreterMutex);
   if (!fMapfile) {
      fMapfile = new TEnv();
      fMapfile->IgnoreDuplicates(kTRUE);
//       fMapNamespaces = new THashTable();
//       fMapNamespaces->SetOwner();

      fRootmapFiles = new TObjArray;
      fRootmapFiles->SetOwner();

      InitRootmapFile(".rootmap");
   }
   //fMapfile->SetValue(key, libs);
   fMapfile->SetValue(cls, libs);
   return 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Demangle the name (from the typeinfo) and then request the class
/// via the usual name based interface (TClass::GetClass).

TClass *TCling::GetClass(const std::type_info& typeinfo, Bool_t load) const
{
   int err = 0;
   char* demangled_name = TClassEdit::DemangleTypeIdName(typeinfo, err);
   if (err) return 0;
   TClass* theClass = TClass::GetClass(demangled_name, load, kTRUE);
   free(demangled_name);
   return theClass;
}

////////////////////////////////////////////////////////////////////////////////
/// Load library containing the specified class. Returns 0 in case of error
/// and 1 in case if success.

Int_t TCling::AutoLoad(const std::type_info& typeinfo, Bool_t knowDictNotLoaded /* = kFALSE */)
{
   int err = 0;
   char* demangled_name_c = TClassEdit::DemangleTypeIdName(typeinfo, err);
   if (err) {
      return 0;
   }

   std::string demangled_name(demangled_name_c);
   free(demangled_name_c);

   // AutoLoad expects (because TClass::GetClass already prepares it that way) a
   // shortened name.
   TClassEdit::TSplitType splitname( demangled_name.c_str(), (TClassEdit::EModType)(TClassEdit::kLong64 | TClassEdit::kDropStd) );
   splitname.ShortType(demangled_name, TClassEdit::kDropStlDefault | TClassEdit::kDropStd);

   // No need to worry about typedef, they aren't any ... but there are
   // inlined namespaces ...

   Int_t result = AutoLoad(demangled_name.c_str());
   if (result == 0) {
      demangled_name = TClassEdit::GetLong64_Name(demangled_name);
      result = AutoLoad(demangled_name.c_str(), knowDictNotLoaded);
   }

   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Load library containing the specified class. Returns 0 in case of error
/// and 1 in case if success.

Int_t TCling::AutoLoad(const char *cls, Bool_t knowDictNotLoaded /* = kFALSE */)
{
   R__LOCKGUARD(gInterpreterMutex);

   if (!knowDictNotLoaded && gClassTable->GetDictNorm(cls)) {
      // The library is already loaded as the class's dictionary is known.
      // Return success.
      // Note: the name (cls) is expected to be normalized as it comes either
      // from a callbacks (that can/should calculate the normalized name from the
      // decl) or from TClass::GetClass (which does also calculate the normalized
      // name).
      return 1;
   }

   if (gDebug > 2) {
      Info("TCling::AutoLoad",
           "Trying to autoload for %s", cls);
   }

   Int_t status = 0;
   if (!gROOT || !gInterpreter || gROOT->TestBit(TObject::kInvalidObject)) {
      if (gDebug > 2) {
         Info("TCling::AutoLoad",
              "Disabled due to gROOT or gInterpreter being invalid/not ready (the class name is %s)", cls);
      }
      return status;
   }
   if (!fAllowLibLoad) {
      // Never load any library from rootcling/genreflex.
      if (gDebug > 2) {
         Info("TCling::AutoLoad",
              "Explicitly disabled (the class name is %s)", cls);
      }
      return 0;
   }
   // Prevent the recursion when the library dictionary are loaded.
   Int_t oldvalue = SetClassAutoloading(false);
   // Try using externally provided callback first.
   if (fAutoLoadCallBack) {
      int success = (*(AutoLoadCallBack_t)fAutoLoadCallBack)(cls);
      if (success) {
         SetClassAutoloading(oldvalue);
         return success;
      }
   }
   // lookup class to find list of dependent libraries
   TString deplibs = GetClassSharedLibs(cls);
   if (!deplibs.IsNull()) {
      TString delim(" ");
      TObjArray* tokens = deplibs.Tokenize(delim);
      for (Int_t i = (tokens->GetEntriesFast() - 1); i > 0; --i) {
         const char* deplib = ((TObjString*)tokens->At(i))->GetName();
         if (gROOT->LoadClass(cls, deplib) == 0) {
            if (gDebug > 0) {
               Info("TCling::AutoLoad",
                    "loaded dependent library %s for %s", deplib, cls);
            }
         }
         else {
            Error("TCling::AutoLoad",
                  "failure loading dependent library %s for %s",
                  deplib, cls);
         }
      }
      const char* lib = ((TObjString*)tokens->At(0))->GetName();
      if (lib && lib[0]) {
         if (gROOT->LoadClass(cls, lib) == 0) {
            if (gDebug > 0) {
               Info("TCling::AutoLoad",
                    "loaded library %s for %s", lib, cls);
            }
            status = 1;
         }
         else {
            Error("TCling::AutoLoad",
                  "failure loading library %s for %s", lib, cls);
         }
      }
      delete tokens;
   }

   SetClassAutoloading(oldvalue);
   return status;
}

////////////////////////////////////////////////////////////////////////////////
/// Parse the payload or header.

static cling::Interpreter::CompilationResult ExecAutoParse(const char *what,
                                                           Bool_t header,
                                                           cling::Interpreter *interpreter)
{
   std::string code = gNonInterpreterClassDef ;
   if (!header) {
      // This is the complete header file content and not the
      // name of a header.
      code += what;

   } else {
      code += ("#include \"");
      code += what;
      code += "\"\n";
   }
   code += ("#ifdef __ROOTCLING__\n"
            "#undef __ROOTCLING__\n"
            + gInterpreterClassDef +
            "#endif");

   cling::Interpreter::CompilationResult cr;
   {
      // scope within which diagnostics are de-activated
      // For now we disable diagnostics because we saw them already at
      // dictionary generation time. That won't be an issue with the PCMs.

      Sema &SemaR = interpreter->getSema();
      ROOT::Internal::ParsingStateRAII parsingStateRAII(interpreter->getParser(), SemaR);
      clangDiagSuppr diagSuppr(SemaR.getDiagnostics());

      #if defined(R__MUST_REVISIT)
      #if R__MUST_REVISIT(6,2)
      Warning("TCling::RegisterModule","Diagnostics suppression should be gone by now.");
      #endif
      #endif

      cr = interpreter->parseForModule(code);
   }
   return cr;
}

////////////////////////////////////////////////////////////////////////////////
/// Helper routine for TCling::AutoParse implementing the actual call to the
/// parser and looping over template parameters (if
/// any) and when they don't have a registered header to autoparse,
/// recurse over their template parameters.
///
/// Returns the number of header parsed.

UInt_t TCling::AutoParseImplRecurse(const char *cls, bool topLevel)
{
   // We assume the lock has already been taken.
   //    R__LOCKGUARD(gInterpreterMutex);

   Int_t nHheadersParsed = 0;
   unsigned long offset = 0;
   if (strncmp(cls, "const ", 6) == 0) {
      offset = 6;
   }

   // Loop on the possible autoparse keys
   bool skipFirstEntry = false;
   std::vector<std::string> autoparseKeys;
   if (strchr(cls, '<')) {
      int nestedLoc = 0;
      TClassEdit::GetSplit(cls + offset, autoparseKeys, nestedLoc, TClassEdit::kDropTrailStar);
      // Check if we can skip the name of the template in the autoparses
      // Take all the scopes one by one. If all of them are in the AST, we do not
      // need to autoparse for that particular template.
      if (!autoparseKeys.empty() && !autoparseKeys[0].empty()) {
         // autoparseKeys[0] is empty when the input is not a template instance.
         // The case strchr(cls, '<') != 0 but still not a template instance can
         // happens 'just' for string (GetSplit replaces the template by the short name
         // and then use that for thew splitting)
         TString templateName(autoparseKeys[0]);
         auto tokens = templateName.Tokenize("::");
         clang::NamedDecl* previousScopeAsNamedDecl = nullptr;
         clang::DeclContext* previousScopeAsContext = fInterpreter->getCI()->getASTContext().getTranslationUnitDecl();
         if (TClassEdit::IsStdClass(cls + offset))
            previousScopeAsContext = fInterpreter->getSema().getStdNamespace();
         auto nTokens = tokens->GetEntries();
         for (Int_t tk = 0; tk < nTokens; ++tk) {
            auto scopeObj = tokens->UncheckedAt(tk);
            auto scopeName = ((TObjString*) scopeObj)->String().Data();
            previousScopeAsNamedDecl = cling::utils::Lookup::Named(&fInterpreter->getSema(), scopeName, previousScopeAsContext);
            // Check if we have multiple nodes in the AST with this name
            if ((clang::NamedDecl*)-1 == previousScopeAsNamedDecl) break;
            previousScopeAsContext = llvm::dyn_cast_or_null<clang::DeclContext>(previousScopeAsNamedDecl);
            if (!previousScopeAsContext) break; // this is not a context
         }
         delete tokens;
         // Now, let's check if the last scope, the template, has a definition, i.e. it's not a fwd decl
         if ((clang::NamedDecl*)-1 != previousScopeAsNamedDecl) {
            if (auto templateDecl = llvm::dyn_cast_or_null<clang::ClassTemplateDecl>(previousScopeAsNamedDecl)) {
               if (auto templatedDecl = templateDecl->getTemplatedDecl()) {
                  skipFirstEntry = templatedDecl->hasDefinition();
               }
            }
         }

      }
   }
   if (topLevel) autoparseKeys.emplace_back(cls);

   for (const auto & apKeyStr : autoparseKeys) {
      if (skipFirstEntry) {
         skipFirstEntry=false;
         continue;
      }
      if (apKeyStr.empty()) continue;
      const char *apKey = apKeyStr.c_str();
      std::size_t normNameHash(fStringHashFunction(apKey));
      // If the class was not looked up
      if (gDebug > 1) {
         Info("TCling::AutoParse",
              "Starting autoparse for %s\n", apKey);
      }
      if (fLookedUpClasses.insert(normNameHash).second) {
         auto const &iter = fClassesHeadersMap.find(normNameHash);
         if (iter != fClassesHeadersMap.end()) {
            const cling::Transaction *T = fInterpreter->getCurrentTransaction();
            fTransactionHeadersMap.insert({T,normNameHash});
            auto const &hNamesPtrs = iter->second;
            if (gDebug > 1) {
               Info("TCling::AutoParse",
                    "We can proceed for %s. We have %s headers.", apKey, std::to_string(hNamesPtrs.size()).c_str());
            }
            for (auto & hName : hNamesPtrs) {
               if (fParsedPayloadsAddresses.count(hName) == 1) continue;
               if (0 != fPayloads.count(normNameHash)) {
                  float initRSSval=0.f, initVSIZEval=0.f;
                  (void) initRSSval; // Avoid unused var warning
                  (void) initVSIZEval;
                  if (gDebug > 0) {
                     Info("AutoParse",
                          "Parsing full payload for %s", apKey);
                     ProcInfo_t info;
                     gSystem->GetProcInfo(&info);
                     initRSSval = 1e-3*info.fMemResident;
                     initVSIZEval = 1e-3*info.fMemVirtual;
                  }
                  auto cRes = ExecAutoParse(hName, kFALSE, fInterpreter);
                  if (cRes != cling::Interpreter::kSuccess) {
                     if (hName[0] == '\n')
                        Error("AutoParse", "Error parsing payload code for class %s with content:\n%s", apKey, hName);
                  } else {
                     fParsedPayloadsAddresses.insert(hName);
                     nHheadersParsed++;
                     if (gDebug > 0){
                        ProcInfo_t info;
                        gSystem->GetProcInfo(&info);
                        float endRSSval = 1e-3*info.fMemResident;
                        float endVSIZEval = 1e-3*info.fMemVirtual;
                        Info("Autoparse", ">>> RSS key %s - before %.3f MB - after %.3f MB - delta %.3f MB", apKey, initRSSval, endRSSval, endRSSval-initRSSval);
                        Info("Autoparse", ">>> VSIZE key %s - before %.3f MB - after %.3f MB - delta %.3f MB", apKey, initVSIZEval, endVSIZEval, endVSIZEval-initVSIZEval);
                     }
                  }
               } else if (!IsLoaded(hName)) {
                  if (gDebug > 0) {
                     Info("AutoParse",
                          "Parsing single header %s", hName);
                  }
                  auto cRes = ExecAutoParse(hName, kTRUE, fInterpreter);
                  if (cRes != cling::Interpreter::kSuccess) {
                     Error("AutoParse", "Error parsing headerfile %s for class %s.", hName, apKey);
                  } else {
                     nHheadersParsed++;
                  }
               }
            }
         }
         else {
            // There is no header registered for this class, if this a
            // template, it will be instantiated if/when it is requested
            // and if we do no load/parse its components we might end up
            // not using an eventual specialization.
            if (strchr(apKey, '<')) {
               nHheadersParsed += AutoParseImplRecurse(apKey, false);
            }
         }
      }
   }

   return nHheadersParsed;

}

////////////////////////////////////////////////////////////////////////////////
/// Parse the headers relative to the class
/// Returns 1 in case of success, 0 in case of failure

Int_t TCling::AutoParse(const char *cls)
{
   R__LOCKGUARD(gInterpreterMutex);

   if (!fHeaderParsingOnDemand || fIsAutoParsingSuspended) {
      if (fClingCallbacks->IsAutoloadingEnabled()) {
         return AutoLoad(cls);
      } else {
         return 0;
      }
   }

   if (gDebug > 1) {
      Info("TCling::AutoParse",
           "Trying to autoparse for %s", cls);
   }

   // The catalogue of headers is in the dictionary
   if (fClingCallbacks->IsAutoloadingEnabled()
         && !gClassTable->GetDictNorm(cls)) {
      // Need RAII against recursive (dictionary payload) parsing (ROOT-8445).
      ROOT::Internal::ParsingStateRAII parsingStateRAII(fInterpreter->getParser(),
         fInterpreter->getSema());
      AutoLoad(cls, true /*knowDictNotLoaded*/);
   }

   // Prevent the recursion when the library dictionary are loaded.
   Int_t oldAutoloadValue = SetClassAutoloading(false);

   // No recursive header parsing on demand; we require headers to be standalone.
   SuspendAutoParsing autoParseRAII(this);

   Int_t nHheadersParsed = AutoParseImplRecurse(cls,/*topLevel=*/ true);

   if (nHheadersParsed != 0) {
      while (!fClassesToUpdate.empty()) {
         TClass *oldcl = fClassesToUpdate.back().first;
         if (oldcl->GetState() != TClass::kHasTClassInit) {
            // if (gDebug > 2) Info("RegisterModule", "Forcing TClass init for %s", oldcl->GetName());
            DictFuncPtr_t dict = fClassesToUpdate.back().second;
            fClassesToUpdate.pop_back();
            // Calling func could manipulate the list so, let maintain the list
            // then call the dictionary function.
            TClass *ncl = dict();
            if (ncl) ncl->PostLoadCheck();
         } else {
            fClassesToUpdate.pop_back();
         }
      }
   }

   SetClassAutoloading(oldAutoloadValue);

   return nHheadersParsed > 0 ? 1 : 0;
}

// This is a function which gets callback from cling when DynamicLibraryManager->loadLibrary failed for some reason.
// Try to solve the problem by autoloading. Return true when autoloading success, return
// false if not.
bool TCling::LibraryLoadingFailed(const std::string& errmessage, const std::string& libStem, bool permanent, bool resolved)
{
   StringRef errMsg(errmessage);
   if (errMsg.contains("undefined symbol: ")) {
   // This branch is taken when the callback was from DynamicLibraryManager::loadLibrary
      std::string mangled_name = std::string(errMsg.split("undefined symbol: ").second);
      void* res = ((TCling*)gCling)->LazyFunctionCreatorAutoload(mangled_name);
      cling::DynamicLibraryManager* DLM = fInterpreter->getDynamicLibraryManager();
      if (res && DLM && (DLM->loadLibrary(libStem, permanent, resolved) == cling::DynamicLibraryManager::kLoadLibSuccess))
        // Return success when LazyFunctionCreatorAutoload could find mangled_name
        return true;
   } else {
   // The callback is from IncrementalExecutor::diagnoseUnresolvedSymbols
      if ( ((TCling*)gCling)->LazyFunctionCreatorAutoload(errmessage))
         return true;
   }

   return false;
}

// This is a GNU implementation of hash used in bloom filter!
static uint32_t GNUHash(StringRef S) {
   uint32_t H = 5381;
   for (uint8_t C : S)
      H = (H << 5) + H + C;
   return H;
}

static StringRef GetGnuHashSection(llvm::object::ObjectFile *file) {
   for (auto S : file->sections()) {
      StringRef name;
      S.getName(name);
      if (name == ".gnu.hash") {
         StringRef content;
         S.getContents(content);
         return content;
      }
   }
   return "";
}

// Bloom filter. See https://blogs.oracle.com/solaris/gnu-hash-elf-sections-v2
// for detailed desctiption. In short, there is a .gnu.hash section in so files which contains
// bloomfilter hash value that we can compare with our mangled_name hash. This is a false positive
// probability data structure and enables us to skip libraries which doesn't contain mangled_name definition!
// PE and Mach-O files doesn't have .gnu.hash bloomfilter section, so this is a specific optimization for ELF.
// This is fine because performance critical part (data centers) are running on Linux :)
static bool LookupBloomFilter(llvm::object::ObjectFile *soFile, uint32_t hash) {
   const int bits = 64;

   StringRef contents = GetGnuHashSection(soFile);
   if (contents.size() < 16)
      // We need to search if the library doesn't have .gnu.hash section!
      return true;
   const char* hashContent = contents.data();

   // See https://flapenguin.me/2017/05/10/elf-lookup-dt-gnu-hash/ for .gnu.hash table layout.
   uint32_t maskWords = *reinterpret_cast<const uint32_t *>(hashContent + 8);
   uint32_t shift2 = *reinterpret_cast<const uint32_t *>(hashContent + 12);
   uint32_t hash2 = hash >> shift2;
   uint32_t n = (hash / bits) % maskWords;

   const char *bloomfilter = hashContent + 16;
   const char *hash_pos = bloomfilter + n*(bits/8); // * (Bits / 8)
   uint64_t word = *reinterpret_cast<const uint64_t *>(hash_pos);
   uint64_t bitmask = ( (1ULL << (hash % bits)) | (1ULL << (hash2 % bits)));
   return  (bitmask & word) == bitmask;
}

// Lookup for normal symbols
static bool LookupNormalSymbols(llvm::object::ObjectFile *RealSoFile, const std::string& mangled_name, const std::string& LibName)
{
   uint32_t hashedMangle = GNUHash(mangled_name);
   // Check Bloom filter. If false, it means that this library doesn't contain mangled_name defenition
   if (!LookupBloomFilter(RealSoFile, hashedMangle))
      return false;

   auto Symbols = RealSoFile->symbols();
   for (auto S : Symbols) {
      uint32_t Flags = S.getFlags();
      // DO NOT insert to table if symbol was undefined
      if (Flags & llvm::object::SymbolRef::SF_Undefined)
         continue;

      // Note, we are at last resort and loading library based on a weak
      // symbol is allowed. Otherwise, the JIT will issue an unresolved
      // symbol error.
      //
      // There are other weak symbol kinds (marked as 'V') to denote
      // typeinfo and vtables. It is unclear whether we should load such
      // libraries or from which library we should resolve the symbol.
      // We seem to not have a way to differentiate it from the symbol API.

      llvm::Expected<StringRef> SymNameErr = S.getName();
      if (!SymNameErr) {
         Warning("LookupNormalSymbols", "Failed to read symbol");
         continue;
      }

      if (SymNameErr.get() == mangled_name)
         return true;
   }

   return false;
}

static void* LazyFunctionCreatorAutoloadForModule(const std::string& mangled_name,
                                                  cling::Interpreter *fInterpreter) {
   using namespace llvm::object;
   using namespace llvm::sys::path;
   using namespace llvm::sys::fs;

   R__LOCKGUARD(gInterpreterMutex);

   static bool sFirstRun = true;
   static bool sFirstSystemLibrary = true;
   // sLibraies contains pair of sPaths[i] (eg. /home/foo/module) and library name (eg. libTMVA.so). The
   // reason why we're separating sLibraries and sPaths is that we have a lot of
   // dupulication in path, for example we have "/home/foo/module-release/lib/libFoo.so", "/home/../libBar.so", "/home/../lib.."
   // and it's waste of memory to store the full path.
   static std::vector< std::pair<uint32_t, std::string> > sLibraries;
   static std::vector<std::string> sPaths;

   // For system header autoloading
   static std::vector< std::pair<uint32_t, std::string> > sSysLibraries;
   static std::vector<std::string> sSysPaths;

   if (sFirstRun) {
      TCling__FindLoadedLibraries(sLibraries, sPaths, *fInterpreter, /* searchSystem */ false);
      sFirstRun = false;
   }

   // The JIT gives us a mangled name which has only one leading underscore on
   // all platforms, for instance _ZN8TRandom34RndmEv. However, on OSX the
   // linker stores this symbol as __ZN8TRandom34RndmEv (adding an extra _).
#ifdef R__MACOSX
   std::string name_in_so = "_" + mangled_name;
#else
   std::string name_in_so = mangled_name;
#endif


   // Iterate over files under this path. We want to get each ".so" files
   for (std::pair<uint32_t, std::string> &P : sLibraries) {
      llvm::SmallString<400> Vec(sPaths[P.first]);
      llvm::sys::path::append(Vec, StringRef(P.second));
      const std::string LibName = Vec.str();

      auto SoFile = ObjectFile::createObjectFile(LibName);
      if (!SoFile)
         continue;

      if (LookupNormalSymbols(SoFile.get().getBinary(), name_in_so, LibName)) {
         if (gSystem->Load(LibName.c_str(), "", false) < 0)
            Error("LazyFunctionCreatorAutoloadForModule", "Failed to load library %s", LibName.c_str());

         // We want to delete a loaded library from sLibraries cache, because sLibraries is
         // a vector of candidate libraries which might be loaded in the future.
         sLibraries.erase(std::remove(sLibraries.begin(), sLibraries.end(), P), sLibraries.end());
         void* addr = llvm::sys::DynamicLibrary::SearchForAddressOfSymbol(mangled_name.c_str());
         return addr;
      }
   }

   // Normal lookup failed! Fall back to system library
   if (sFirstSystemLibrary) {
      TCling__FindLoadedLibraries(sSysLibraries, sSysPaths, *fInterpreter, /* searchSystem */ true);
      sFirstSystemLibrary = false;
   }

   for (std::pair<uint32_t, std::string> &P : sSysLibraries) {
      llvm::SmallString<400> Vec(sSysPaths[P.first]);
      llvm::sys::path::append(Vec, StringRef(P.second));
      const std::string LibName = Vec.str();

      auto SoFile = ObjectFile::createObjectFile(LibName);
      if (!SoFile)
         continue;

      auto RealSoFile = SoFile.get().getBinary();

      if (LookupNormalSymbols(RealSoFile, name_in_so, LibName)) {
         if (gSystem->Load(LibName.c_str(), "", false) < 0)
            Error("LazyFunctionCreatorAutoloadForModule", "Failed to load library %s", LibName.c_str());

         sSysLibraries.erase(std::remove(sSysLibraries.begin(), sSysLibraries.end(), P), sSysLibraries.end());
         void* addr = llvm::sys::DynamicLibrary::SearchForAddressOfSymbol(mangled_name.c_str());
         return addr;
      }

      // Lookup for dynamic symbols
      auto sections = RealSoFile->sections();
      for (auto section : sections) {
         llvm::StringRef sectionName;
         section.getName(sectionName);

         // .dynstr contains string of dynamic symbols
         if (sectionName == ".dynstr") {
            llvm::StringRef dContents;
            section.getContents(dContents);
            // If this library contains mangled name
            if (dContents.contains(mangled_name)) {
               if (gSystem->Load(LibName.c_str(), "", false) < 0)
                  Error("LazyFunctionCreatorAutoloadForModule", "Failed to load library %s", LibName.c_str());

               // Delete a loaded library from sLibraries cache.
               sSysLibraries.erase(std::remove(sSysLibraries.begin(), sSysLibraries.end(), P), sSysLibraries.end());
               void* addr = llvm::sys::DynamicLibrary::SearchForAddressOfSymbol(mangled_name.c_str());
               return addr;
            }
         }
      }
   }

   // Lookup failed!!!!
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Autoload a library based on a missing symbol.

void* TCling::LazyFunctionCreatorAutoload(const std::string& mangled_name) {
   if (fCxxModulesEnabled)
      return LazyFunctionCreatorAutoloadForModule(mangled_name, fInterpreter);

   // First see whether the symbol is in the library that we are currently
   // loading. It will have access to the symbols of its dependent libraries,
   // thus checking "back()" is sufficient.
   if (!fRegisterModuleDyLibs.empty()) {
      if (void* addr = dlsym(fRegisterModuleDyLibs.back(),
                             mangled_name.c_str())) {
         return addr;
      }
   }

   int err = 0;
   char* demangled_name_c = TClassEdit::DemangleName(mangled_name.c_str(), err);
   if (err) {
      return 0;
   }

   std::string name(demangled_name_c);
   free(demangled_name_c);

   //fprintf(stderr, "demangled name: '%s'\n", demangled_name);
   //
   //  Separate out the class or namespace part of the
   //  function name.
   //

   std::string::size_type pos = name.find("__thiscall ");
   if (pos != std::string::npos) {
      name.erase(0, pos + sizeof("__thiscall ")-1);
   }
   pos = name.find("__cdecl ");
   if (pos != std::string::npos) {
      name.erase(0, pos + sizeof("__cdecl ")-1);
   }
   if (!strncmp(name.c_str(), "typeinfo for ", sizeof("typeinfo for ")-1)) {
      name.erase(0, sizeof("typeinfo for ")-1);
   } else if (!strncmp(name.c_str(), "vtable for ", sizeof("vtable for ")-1)) {
      name.erase(0, sizeof("vtable for ")-1);
   } else if (!strncmp(name.c_str(), "operator", sizeof("operator")-1)
              && !isalnum(name[sizeof("operator")])) {
     // operator...(A, B) - let's try with A!
     name.erase(0, sizeof("operator")-1);
     pos = name.rfind('(');
     if (pos != std::string::npos) {
       name.erase(0, pos + 1);
       pos = name.find(",");
       if (pos != std::string::npos) {
         // remove next arg up to end, leaving only the first argument type.
         name.erase(pos);
       }
       pos = name.rfind(" const");
       if (pos != std::string::npos) {
         name.erase(pos, strlen(" const"));
       }
       while (!name.empty() && strchr("&*", name.back()))
         name.erase(name.length() - 1);
     }
   } else {
      TClassEdit::FunctionSplitInfo fsi;
      TClassEdit::SplitFunction(name, fsi);
      name = fsi.fScopeName;
   }
   //fprintf(stderr, "name: '%s'\n", name.c_str());
   // Now we have the class or namespace name, so do the lookup.
   TString libs = GetClassSharedLibs(name.c_str());
   if (libs.IsNull()) {
      // Not found in the map, all done.
      return 0;
   }
   //fprintf(stderr, "library: %s\n", iter->second.c_str());
   // Now we have the name of the libraries to load, so load them.

   TString lib;
   Ssiz_t posLib = 0;
   while (libs.Tokenize(lib, posLib)) {
      if (gSystem->Load(lib, "", kFALSE /*system*/) < 0) {
         // The library load failed, all done.
         //fprintf(stderr, "load failed: %s\n", errmsg.c_str());
         return 0;
      }
   }

   //fprintf(stderr, "load succeeded.\n");
   // Get the address of the function being called.
   void* addr = llvm::sys::DynamicLibrary::SearchForAddressOfSymbol(mangled_name.c_str());
   //fprintf(stderr, "addr: %016lx\n", reinterpret_cast<unsigned long>(addr));
   return addr;
}

////////////////////////////////////////////////////////////////////////////////

Bool_t TCling::IsAutoLoadNamespaceCandidate(const char* name)
{
//    if (fMapNamespaces){
//       return fMapNamespaces->FindObject(name);
//    }
   return false;
}

////////////////////////////////////////////////////////////////////////////////

Bool_t TCling::IsAutoLoadNamespaceCandidate(const clang::NamespaceDecl* nsDecl)
{
   return fNSFromRootmaps.count(nsDecl) != 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Internal function. Inform a TClass about its new TagDecl or NamespaceDecl.

void TCling::UpdateClassInfoWithDecl(const void* vTD)
{
   const NamedDecl* ND = static_cast<const NamedDecl*>(vTD);
   const TagDecl* td = dyn_cast<TagDecl>(ND);
   std::string name;
   TagDecl* tdDef = 0;
   if (td) {
      tdDef = td->getDefinition();
      // Let's pass the decl to the TClass only if it has a definition.
      if (!tdDef) return;
      td = tdDef;
      ND = td;

      if (llvm::isa<clang::FunctionDecl>(td->getDeclContext())) {
         // Ignore declaration within a function.
         return;
      }
      clang::QualType type( td->getTypeForDecl(), 0 );

      auto declName=ND->getNameAsString();
      if (!TClass::HasNoInfoOrEmuOrFwdDeclaredDecl(declName.c_str())){
//          printf ("Impossible to find a TClassEntry in kNoInfo or kEmulated the decl of which would be called %s. Skip w/o building the normalized name.\n",declName );
         return;
      }

      ROOT::TMetaUtils::GetNormalizedName(name, type, *fInterpreter, *fNormalizedCtxt);
   } else {
      name = ND->getNameAsString();
   }

   // Supposedly we are being called while something is being
   // loaded ... let's now tell the autoloader to do the work
   // yet another time.
   int storedAutoloading = SetClassAutoloading(false);
   // FIXME: There can be more than one TClass for a single decl.
   // for example vector<double> and vector<Double32_t>
   TClass* cl = (TClass*)gROOT->GetListOfClasses()->FindObject(name.c_str());
   if (cl && GetModTClasses().find(cl) == GetModTClasses().end()) {
      TClingClassInfo* cci = ((TClingClassInfo*)cl->fClassInfo);
      if (cci) {
         // If we only had a forward declaration then update the
         // TClingClassInfo with the definition if we have it now.
         const TagDecl* tdOld = llvm::dyn_cast_or_null<TagDecl>(cci->GetDecl());
         if (!tdOld || (tdDef && tdDef != tdOld)) {
            cl->ResetCaches();
            TClass::RemoveClassDeclId(cci->GetDeclId());
            if (td) {
               // It's a tag decl, not a namespace decl.
               cci->Init(*cci->GetType());
               TClass::AddClassToDeclIdMap(cci->GetDeclId(), cl);
            }
         }
      } else if (!cl->TestBit(TClass::kLoading) && !cl->fHasRootPcmInfo) {
         cl->ResetCaches();
         // yes, this is almost a waste of time, but we do need to lookup
         // the 'type' corresponding to the TClass anyway in order to
         // preserve the opaque typedefs (Double32_t)
         cl->fClassInfo = (ClassInfo_t *)new TClingClassInfo(fInterpreter, cl->GetName());
         // We now need to update the state and bits.
         if (cl->fState != TClass::kHasTClassInit) {
            // if (!cl->fClassInfo->IsValid()) cl->fState = TClass::kForwardDeclared; else
            cl->fState = TClass::kInterpreted;
            cl->ResetBit(TClass::kIsEmulation);
         }
         TClass::AddClassToDeclIdMap(((TClingClassInfo*)(cl->fClassInfo))->GetDeclId(), cl);
      }
   }
   SetClassAutoloading(storedAutoloading);
}

////////////////////////////////////////////////////////////////////////////////
/// No op: see TClingCallbacks

void TCling::UpdateClassInfo(char* item, Long_t tagnum)
{
}

//______________________________________________________________________________
//FIXME: Factor out that function in TClass, because TClass does it already twice
void TCling::UpdateClassInfoWork(const char* item)
{
   // This is a no-op as part of the API.
   // TCling uses UpdateClassInfoWithDecl() instead.
}

////////////////////////////////////////////////////////////////////////////////
/// Update all canvases at end the terminal input command.

void TCling::UpdateAllCanvases()
{
   TIter next(gROOT->GetListOfCanvases());
   TVirtualPad* canvas;
   while ((canvas = (TVirtualPad*)next())) {
      canvas->Update();
   }
}

////////////////////////////////////////////////////////////////////////////////

void TCling::UpdateListsOnCommitted(const cling::Transaction &T) {
   std::set<TClass*> modifiedTClasses; // TClasses that require update after this transaction

   // If the transaction does not contain anything we can return earlier.
   if (!HandleNewTransaction(T)) return;

   bool isTUTransaction = false;
   if (!T.empty() && T.decls_begin() + 1 == T.decls_end() && !T.hasNestedTransactions()) {
      clang::Decl* FirstDecl = *(T.decls_begin()->m_DGR.begin());
      if (llvm::isa<clang::TranslationUnitDecl>(FirstDecl)) {
         // The is the first transaction, we have to expose to meta
         // what's already in the AST.
         isTUTransaction = true;
      }
   }

   std::set<const void*> TransactionDeclSet;
   if (!isTUTransaction && T.decls_end() - T.decls_begin()) {
      const clang::Decl* WrapperFD = T.getWrapperFD();
      for (cling::Transaction::const_iterator I = T.decls_begin(), E = T.decls_end();
          I != E; ++I) {
         if (I->m_Call != cling::Transaction::kCCIHandleTopLevelDecl
             && I->m_Call != cling::Transaction::kCCIHandleTagDeclDefinition)
            continue;

         for (DeclGroupRef::const_iterator DI = I->m_DGR.begin(),
                 DE = I->m_DGR.end(); DI != DE; ++DI) {
            if (*DI == WrapperFD)
               continue;
            TransactionDeclSet.insert(*DI);
            ((TCling*)gCling)->HandleNewDecl(*DI, false, modifiedTClasses);
         }
      }
   }

   // The above might trigger more decls to be deserialized.
   // Thus the iteration over the deserialized decls must be last.
   for (cling::Transaction::const_iterator I = T.deserialized_decls_begin(),
           E = T.deserialized_decls_end(); I != E; ++I) {
      for (DeclGroupRef::const_iterator DI = I->m_DGR.begin(),
              DE = I->m_DGR.end(); DI != DE; ++DI)
         if (TransactionDeclSet.find(*DI) == TransactionDeclSet.end()) {
            //FIXME: HandleNewDecl should take DeclGroupRef
            ((TCling*)gCling)->HandleNewDecl(*DI, /*isDeserialized*/true,
                                             modifiedTClasses);
         }
   }


   // When fully building the reflection info in TClass, a deserialization
   // could be triggered, which may result in request for building the
   // reflection info for the same TClass. This in turn will clear the caches
   // for the TClass in-flight and cause null ptr derefs.
   // FIXME: This is a quick fix, solving most of the issues. The actual
   // question is: Shouldn't TClass provide a lock mechanism on update or lock
   // itself until the update is done.
   //
   std::vector<TClass*> modifiedTClassesDiff(modifiedTClasses.size());
   std::vector<TClass*>::iterator it;
   it = set_difference(modifiedTClasses.begin(), modifiedTClasses.end(),
                       ((TCling*)gCling)->GetModTClasses().begin(),
                       ((TCling*)gCling)->GetModTClasses().end(),
                       modifiedTClassesDiff.begin());
   modifiedTClassesDiff.resize(it - modifiedTClassesDiff.begin());

   // Lock the TClass for updates
   ((TCling*)gCling)->GetModTClasses().insert(modifiedTClassesDiff.begin(),
                                              modifiedTClassesDiff.end());
   for (std::vector<TClass*>::const_iterator I = modifiedTClassesDiff.begin(),
           E = modifiedTClassesDiff.end(); I != E; ++I) {
      // Make sure the TClass has not been deleted.
      if (!gROOT->GetListOfClasses()->FindObject(*I)) {
         continue;
      }
      // Could trigger deserialization of decls.
      cling::Interpreter::PushTransactionRAII RAII(fInterpreter);
      // Unlock the TClass for updates
      ((TCling*)gCling)->GetModTClasses().erase(*I);

   }
}

////////////////////////////////////////////////////////////////////////////////
/// Helper function to go through the members of a class or namespace and unload them.

void TCling::UnloadClassMembers(TClass* cl, const clang::DeclContext* DC) {

   TDataMember* var = 0;
   TFunction* function = 0;
   TEnum* e = 0;
   TFunctionTemplate* functiontemplate = 0;
   TListOfDataMembers* datamembers = (TListOfDataMembers*)cl->GetListOfDataMembers();
   TListOfFunctions* functions = (TListOfFunctions*)cl->GetListOfMethods();
   TListOfEnums* enums = (TListOfEnums*)cl->GetListOfEnums();
   TListOfFunctionTemplates* functiontemplates = (TListOfFunctionTemplates*)cl->GetListOfFunctionTemplates();
   for (DeclContext::decl_iterator RI = DC->decls_begin(), RE = DC->decls_end(); RI != RE; ++RI) {
      if (isa<VarDecl>(*RI) || isa<FieldDecl>(*RI)) {
         const clang::ValueDecl* VD = dyn_cast<ValueDecl>(*RI);
         var = (TDataMember*)datamembers->FindObject(VD->getNameAsString().c_str());
         if (var) {
            // Unload the global by setting the DataMemberInfo_t to 0
            datamembers->Unload(var);
            var->Update(0);
         }
      } else if (const FunctionDecl* FD = dyn_cast<FunctionDecl>(*RI)) {
         function = (TFunction*)functions->FindObject(FD->getNameAsString().c_str());
         if (function) {
            functions->Unload(function);
            function->Update(0);
         }
      } else if (const EnumDecl* ED = dyn_cast<EnumDecl>(*RI)) {
         e = (TEnum*)enums->FindObject(ED->getNameAsString().c_str());
         if (e) {
            TIter iEnumConst(e->GetConstants());
            while (TEnumConstant* enumConst = (TEnumConstant*)iEnumConst()) {
               // Since the enum is already created and valid that ensures us that
               // we have the enum constants created as well.
               enumConst = (TEnumConstant*)datamembers->FindObject(enumConst->GetName());
               if (enumConst && enumConst->IsValid()) {
                  datamembers->Unload(enumConst);
                  enumConst->Update(0);
               }
            }
            enums->Unload(e);
            e->Update(0);
         }
      } else if (const FunctionTemplateDecl* FTD = dyn_cast<FunctionTemplateDecl>(*RI)) {
         functiontemplate = (TFunctionTemplate*)functiontemplates->FindObject(FTD->getNameAsString().c_str());
         if (functiontemplate) {
            functiontemplates->Unload(functiontemplate);
            functiontemplate->Update(0);
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////

void TCling::UpdateListsOnUnloaded(const cling::Transaction &T)
{
   HandleNewTransaction(T);

   // Unload the objects from the lists and update the objects' state.
   TListOfFunctions* functions = (TListOfFunctions*)gROOT->GetListOfGlobalFunctions();
   TListOfFunctionTemplates* functiontemplates = (TListOfFunctionTemplates*)gROOT->GetListOfFunctionTemplates();
   TListOfEnums* enums = (TListOfEnums*)gROOT->GetListOfEnums();
   TListOfDataMembers* globals = (TListOfDataMembers*)gROOT->GetListOfGlobals();
   cling::Transaction::const_nested_iterator iNested = T.nested_begin();
   for(cling::Transaction::const_iterator I = T.decls_begin(), E = T.decls_end();
       I != E; ++I) {
      if (I->m_Call == cling::Transaction::kCCIHandleVTable)
         continue;

      if (I->m_Call == cling::Transaction::kCCINone) {
         UpdateListsOnUnloaded(**iNested);
         ++iNested;
         continue;
      }

      for (DeclGroupRef::const_iterator DI = I->m_DGR.begin(),
              DE = I->m_DGR.end(); DI != DE; ++DI) {

         // Do not mark a decl as unloaded if we are going to keep it
         // (because it comes from the pch) ...
         if ( (*DI)->isFromASTFile() )
            continue;

         // Deal with global variables and global enum constants.
         if (isa<VarDecl>(*DI) || isa<EnumConstantDecl>(*DI)) {
            TObject *obj = globals->Find((TListOfDataMembers::DeclId_t)*DI);
            if (globals->GetClass()) {
               TDataMember* var = dynamic_cast<TDataMember*>(obj);
               if (var && var->IsValid()) {
                  // Unload the global by setting the DataMemberInfo_t to 0
                  globals->Unload(var);
                  var->Update(0);
               }
            } else {
               TGlobal *g = dynamic_cast<TGlobal*>(obj);
               if (g && g->IsValid()) {
                  // Unload the global by setting the DataMemberInfo_t to 0
                  globals->Unload(g);
                  g->Update(0);
               }
            }
         // Deal with global functions.
         } else if (const FunctionDecl* FD = dyn_cast<FunctionDecl>(*DI)) {
            TFunction* function = (TFunction*)functions->Find((TListOfFunctions::DeclId_t)FD);
            if (function && function->IsValid()) {
               functions->Unload(function);
               function->Update(0);
            }
         // Deal with global function templates.
         } else if (const FunctionTemplateDecl* FTD = dyn_cast<FunctionTemplateDecl>(*DI)) {
            TFunctionTemplate* functiontemplate = (TFunctionTemplate*)functiontemplates->FindObject(FTD->getNameAsString().c_str());
            if (functiontemplate) {
               functiontemplates->Unload(functiontemplate);
               functiontemplate->Update(0);
            }
         // Deal with global enums.
         } else if (const EnumDecl* ED = dyn_cast<EnumDecl>(*DI)) {
            if (TEnum* e = (TEnum*)enums->Find((TListOfEnums::DeclId_t)ED)) {
               globals = (TListOfDataMembers*)gROOT->GetListOfGlobals();
               TIter iEnumConst(e->GetConstants());
               while (TEnumConstant* enumConst = (TEnumConstant*)iEnumConst()) {
                  // Since the enum is already created and valid that ensures us that
                  // we have the enum constants created as well.
                  enumConst = (TEnumConstant*)globals->FindObject(enumConst->GetName());
                  if (enumConst) {
                     globals->Unload(enumConst);
                     enumConst->Update(0);
                  }
               }
               enums->Unload(e);
               e->Update(0);
            }
         // Deal with classes. Unload the class and the data members will be not accessible anymore
         // Cannot declare the members in a different declaration like redeclarable namespaces.
         } else if (const clang::RecordDecl* RD = dyn_cast<RecordDecl>(*DI)) {
            std::vector<TClass*> vectTClass;
            // Only update the TClass if the definition is being unloaded.
            if (RD->isCompleteDefinition()) {
               if (TClass::GetClass(RD, vectTClass)) {
                  for (std::vector<TClass*>::iterator CI = vectTClass.begin(), CE = vectTClass.end();
                       CI != CE; ++CI) {
                     UnloadClassMembers((*CI), RD);
                     (*CI)->ResetClassInfo();
                  }
               }
            }
         // Deal with namespaces. Unload the members of the current redeclaration only.
         } else if (const clang::NamespaceDecl* ND = dyn_cast<NamespaceDecl>(*DI)) {
            std::vector<TClass*> vectTClass;
            if (TClass::GetClass(ND->getCanonicalDecl(), vectTClass)) {
               for (std::vector<TClass*>::iterator CI = vectTClass.begin(), CE = vectTClass.end();
                    CI != CE; ++CI) {
                  UnloadClassMembers((*CI), ND);
                  if (ND->isOriginalNamespace()) {
                     (*CI)->ResetClassInfo();
                  }
               }
            }
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
// If an autoparse was done during a transaction and that it is rolled back,
// we need to make sure the next request for the same autoparse will be
// honored.
void TCling::TransactionRollback(const cling::Transaction &T) {
   auto const &triter = fTransactionHeadersMap.find(&T);
   if (triter != fTransactionHeadersMap.end()) {
      std::size_t normNameHash = triter->second;

      fLookedUpClasses.erase(normNameHash);

      auto const &iter = fClassesHeadersMap.find(normNameHash);
      if (iter != fClassesHeadersMap.end()) {
         auto const &hNamesPtrs = iter->second;
         for (auto &hName : hNamesPtrs) {
            if (gDebug > 0) {
               Info("TransactionRollback",
                    "Restoring ability to autoaparse: %s", hName);
            }
            fParsedPayloadsAddresses.erase(hName);
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////

void TCling::LibraryLoaded(const void* dyLibHandle, const char* canonicalName) {
// UpdateListOfLoadedSharedLibraries();
}

////////////////////////////////////////////////////////////////////////////////

void TCling::LibraryUnloaded(const void* dyLibHandle, const char* canonicalName) {
   fPrevLoadedDynLibInfo = 0;
   fSharedLibs = "";
}

////////////////////////////////////////////////////////////////////////////////
/// Return the list of shared libraries loaded into the process.

const char* TCling::GetSharedLibs()
{
   UpdateListOfLoadedSharedLibraries();
   return fSharedLibs;
}

////////////////////////////////////////////////////////////////////////////////
/// Get the list of shared libraries containing the code for class cls.
/// The first library in the list is the one containing the class, the
/// others are the libraries the first one depends on. Returns 0
/// in case the library is not found.

const char* TCling::GetClassSharedLibs(const char* cls)
{
   if (!cls || !*cls) {
      return 0;
   }
   // lookup class to find list of libraries
   if (fMapfile) {
      TEnvRec* libs_record = 0;
      libs_record = fMapfile->Lookup(cls);
      if (libs_record) {
         const char* libs = libs_record->GetValue();
         return (*libs) ? libs : 0;
      }
      else {
         // Try the old format...
         TString c = TString("Library.") + cls;
         // convert "::" to "@@", we used "@@" because TEnv
         // considers "::" a terminator
         c.ReplaceAll("::", "@@");
         // convert "-" to " ", since class names may have
         // blanks and TEnv considers a blank a terminator
         c.ReplaceAll(" ", "-");
         // Use TEnv::Lookup here as the rootmap file must start with Library.
         // and do not support using any stars (so we do not need to waste time
         // with the search made by TEnv::GetValue).
         TEnvRec* libs_record = 0;
         libs_record = fMapfile->Lookup(c);
         if (libs_record) {
            const char* libs = libs_record->GetValue();
            return (*libs) ? libs : 0;
         }
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Get the list a libraries on which the specified lib depends. The
/// returned string contains as first element the lib itself.
/// Returns 0 in case the lib does not exist or does not have
/// any dependencies.

const char* TCling::GetSharedLibDeps(const char* lib)
{
   if (!fMapfile || !lib || !lib[0]) {
      return 0;
   }
   TString libname(lib);
   Ssiz_t idx = libname.Last('.');
   if (idx != kNPOS) {
      libname.Remove(idx);
   }
   TEnvRec* rec;
   TIter next(fMapfile->GetTable());
   size_t len = libname.Length();
   while ((rec = (TEnvRec*) next())) {
      const char* libs = rec->GetValue();
      if (!strncmp(libs, libname.Data(), len) && strlen(libs) >= len
            && (!libs[len] || libs[len] == ' ' || libs[len] == '.')) {
         return libs;
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// If error messages are disabled, the interpreter should suppress its
/// failures and warning messages from stdout.

Bool_t TCling::IsErrorMessagesEnabled() const
{
#if defined(R__MUST_REVISIT)
#if R__MUST_REVISIT(6,2)
   Warning("IsErrorMessagesEnabled", "Interface not available yet.");
#endif
#endif
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// If error messages are disabled, the interpreter should suppress its
/// failures and warning messages from stdout. Return the previous state.

Bool_t TCling::SetErrorMessages(Bool_t enable)
{
#if defined(R__MUST_REVISIT)
#if R__MUST_REVISIT(6,2)
   Warning("SetErrorMessages", "Interface not available yet.");
#endif
#endif
   return TCling::IsErrorMessagesEnabled();
}

////////////////////////////////////////////////////////////////////////////////
/// Refresh the list of include paths known to the interpreter and return it
/// with -I prepended.

const char* TCling::GetIncludePath()
{
   R__LOCKGUARD(gInterpreterMutex);

   fIncludePath = "";

   llvm::SmallVector<std::string, 10> includePaths;//Why 10? Hell if I know.
   //false - no system header, true - with flags.
   fInterpreter->GetIncludePaths(includePaths, false, true);
   if (const size_t nPaths = includePaths.size()) {
      assert(!(nPaths & 1) && "GetIncludePath, number of paths and options is not equal");

      for (size_t i = 0; i < nPaths; i += 2) {
         if (i)
            fIncludePath.Append(' ');
         fIncludePath.Append(includePaths[i].c_str());

         if (includePaths[i] != "-I")
            fIncludePath.Append(' ');
         fIncludePath.Append('"');
         fIncludePath.Append(includePaths[i + 1], includePaths[i + 1].length());
         fIncludePath.Append('"');
      }
   }

   return fIncludePath;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the directory containing CINT's stl cintdlls.

const char* TCling::GetSTLIncludePath() const
{
   return "";
}

//______________________________________________________________________________
//                      M I S C
//______________________________________________________________________________

int TCling::DisplayClass(FILE* /*fout*/, const char* /*name*/, int /*base*/, int /*start*/) const
{
   // Interface to cling function
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Interface to cling function

int TCling::DisplayIncludePath(FILE *fout) const
{
   assert(fout != 0 && "DisplayIncludePath, 'fout' parameter is null");

   llvm::SmallVector<std::string, 10> includePaths;//Why 10? Hell if I know.
   //false - no system header, true - with flags.
   fInterpreter->GetIncludePaths(includePaths, false, true);
   if (const size_t nPaths = includePaths.size()) {
      assert(!(nPaths & 1) && "DisplayIncludePath, number of paths and options is not equal");

      std::string allIncludes("include path:");
      for (size_t i = 0; i < nPaths; i += 2) {
         allIncludes += ' ';
         allIncludes += includePaths[i];

         if (includePaths[i] != "-I")
            allIncludes += ' ';
         allIncludes += includePaths[i + 1];
      }

      fprintf(fout, "%s\n", allIncludes.c_str());
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Interface to cling function

void* TCling::FindSym(const char* entry) const
{
   return fInterpreter->getAddressOfGlobal(entry);
}

////////////////////////////////////////////////////////////////////////////////
/// Let the interpreter issue a generic error, and set its error state.

void TCling::GenericError(const char* error) const
{
#if defined(R__MUST_REVISIT)
#if R__MUST_REVISIT(6,2)
   Warning("GenericError","Interface not available yet.");
#endif
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// This routines used to return the address of the internal wrapper
/// function (of the interpreter) that was used to call *all* the
/// interpreted functions that were bytecode compiled (no longer
/// interpreted line by line).  In Cling, there is no such
/// wrapper function.
/// In practice this routines was use to decipher whether the
/// pointer returns by InterfaceMethod could be used to uniquely
/// represent the function.  In Cling if the function is in a
/// useable state (its compiled version is available), this is
/// always the case.
/// See TClass::GetMethod.

Long_t TCling::GetExecByteCode() const
{
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Interface to cling function

int TCling::GetSecurityError() const
{
#if defined(R__MUST_REVISIT)
#if R__MUST_REVISIT(6,2)
   Warning("GetSecurityError", "Interface not available yet.");
#endif
#endif
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Load a source file or library called path into the interpreter.

int TCling::LoadFile(const char* path) const
{
   cling::Interpreter::CompilationResult compRes;
   HandleInterpreterException(fMetaProcessor, TString::Format(".L %s", path), compRes, /*cling::Value*/0);
   return compRes == cling::Interpreter::kFailure;
}

////////////////////////////////////////////////////////////////////////////////
/// Load the declarations from text into the interpreter.
/// Note that this cannot be (top level) statements; text must contain
/// top level declarations.
/// Returns true on success, false on failure.

Bool_t TCling::LoadText(const char* text) const
{
   return (fInterpreter->declare(text) == cling::Interpreter::kSuccess);
}

////////////////////////////////////////////////////////////////////////////////
/// Interface to cling function

const char* TCling::MapCppName(const char* name) const
{
   TTHREAD_TLS_DECL(std::string,buffer);
   ROOT::TMetaUtils::GetCppName(buffer,name);
   return buffer.c_str();
}

////////////////////////////////////////////////////////////////////////////////
/// [Place holder for Mutex Lock]
/// Provide the interpreter with a way to
/// acquire a lock used to protect critical section
/// of its code (non-thread safe parts).

void TCling::SetAlloclockfunc(void (* /* p */ )()) const
{
   // nothing to do for now.
}

////////////////////////////////////////////////////////////////////////////////
/// [Place holder for Mutex Unlock] Provide the interpreter with a way to
/// release a lock used to protect critical section
/// of its code (non-thread safe parts).

void TCling::SetAllocunlockfunc(void (* /* p */ )()) const
{
   // nothing to do for now.
}

////////////////////////////////////////////////////////////////////////////////
/// Enable/Disable the Autoloading of libraries.
/// Returns the old value, i.e whether it was enabled or not.

int TCling::SetClassAutoloading(int autoload) const
{
   if (!autoload && !fClingCallbacks) return false;
   if (!fAllowLibLoad) return false;

   assert(fClingCallbacks && "We must have callbacks!");
   bool oldVal =  fClingCallbacks->IsAutoloadingEnabled();
   fClingCallbacks->SetAutoloadingEnabled(autoload);
   return oldVal;
}

////////////////////////////////////////////////////////////////////////////////
/// Enable/Disable the Autoparsing of headers.
/// Returns the old value, i.e whether it was enabled or not.

int TCling::SetClassAutoparsing(int autoparse)
{
   bool oldVal = fHeaderParsingOnDemand;
   fHeaderParsingOnDemand = autoparse;
   return oldVal;
}

////////////////////////////////////////////////////////////////////////////////
/// Suspend the Autoparsing of headers.
/// Returns the old value, i.e whether it was suspended or not.

Bool_t TCling::SetSuspendAutoParsing(Bool_t value) {
   Bool_t old = fIsAutoParsingSuspended;
   fIsAutoParsingSuspended = value;
   if (fClingCallbacks) fClingCallbacks->SetAutoParsingSuspended(value);
   return old;
}

////////////////////////////////////////////////////////////////////////////////
/// Set a callback to receive error messages.

void TCling::SetErrmsgcallback(void* p) const
{
#if defined(R__MUST_REVISIT)
#if R__MUST_REVISIT(6,2)
   Warning("SetErrmsgcallback", "Interface not available yet.");
#endif
#endif
}


////////////////////////////////////////////////////////////////////////////////
/// Create / close a scope for temporaries. No-op for cling; use
/// cling::Value instead.

void TCling::SetTempLevel(int val) const
{
}

////////////////////////////////////////////////////////////////////////////////

int TCling::UnloadFile(const char* path) const
{
   cling::DynamicLibraryManager* DLM = fInterpreter->getDynamicLibraryManager();
   std::string canonical = DLM->lookupLibrary(path);
   if (canonical.empty()) {
      canonical = path;
   }
   // Unload a shared library or a source file.
   cling::Interpreter::CompilationResult compRes;
   HandleInterpreterException(fMetaProcessor, Form(".U %s", canonical.c_str()), compRes, /*cling::Value*/0);
   return compRes == cling::Interpreter::kFailure;
}

std::unique_ptr<TInterpreterValue> TCling::MakeInterpreterValue() const {
   return std::unique_ptr<TInterpreterValue>(new TClingValue);
}

////////////////////////////////////////////////////////////////////////////////
/// The call to Cling's tab complition.

void TCling::CodeComplete(const std::string& line, size_t& cursor,
                          std::vector<std::string>& completions)
{
   fInterpreter->codeComplete(line, cursor, completions);
}

////////////////////////////////////////////////////////////////////////////////
/// Get the interpreter value corresponding to the statement.
int TCling::Evaluate(const char* code, TInterpreterValue& value)
{
   auto V = reinterpret_cast<cling::Value*>(value.GetValAddr());
   auto compRes = fInterpreter->evaluate(code, *V);
   return compRes!=cling::Interpreter::kSuccess ? 0 : 1 ;
}

////////////////////////////////////////////////////////////////////////////////

void TCling::RegisterTemporary(const TInterpreterValue& value)
{
   using namespace cling;
   const Value* V = reinterpret_cast<const Value*>(value.GetValAddr());
   RegisterTemporary(*V);
}

////////////////////////////////////////////////////////////////////////////////
/// Register value as a temporary, extending its lifetime to that of the
/// interpreter. This is needed for TCling's compatibility interfaces
/// returning long - the address of the temporary objects.
/// As such, "simple" types don't need to be stored; they are returned by
/// value; only pointers / references / objects need to be stored.

void TCling::RegisterTemporary(const cling::Value& value)
{
   if (value.isValid() && value.needsManagedAllocation()) {
      R__LOCKGUARD(gInterpreterMutex);
      fTemporaries->push_back(value);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// If the interpreter encounters Name, check whether that is an object ROOT
/// could retrieve. To not re-read objects from disk, cache the name/object
/// pair for a given LookupCtx.

TObject* TCling::GetObjectAddress(const char *Name, void *&LookupCtx)
{
   // The call to FindSpecialObject might induces any kind of use
   // of the interpreter ... (library loading, function calling, etc.)
   // ... and we _know_ we are in the middle of parsing, so let's make
   // sure to save the state and then restore it.

   if (gDirectory) {
      auto iSpecObjMap = fSpecialObjectMaps.find(gDirectory);
      if (iSpecObjMap != fSpecialObjectMaps.end()) {
         auto iSpecObj = iSpecObjMap->second.find(Name);
         if (iSpecObj != iSpecObjMap->second.end()) {
            LookupCtx = gDirectory;
            return iSpecObj->second;
         }
      }
   }

   // Save state of the PP
   Sema &SemaR = fInterpreter->getSema();
   ASTContext& C = SemaR.getASTContext();
   Preprocessor &PP = SemaR.getPreprocessor();
   Parser& P = const_cast<Parser&>(fInterpreter->getParser());
   Preprocessor::CleanupAndRestoreCacheRAII cleanupRAII(PP);
   Parser::ParserCurTokRestoreRAII savedCurToken(P);
   // After we have saved the token reset the current one to something which
   // is safe (semi colon usually means empty decl)
   Token& Tok = const_cast<Token&>(P.getCurToken());
   Tok.setKind(tok::semi);

   // We can't PushDeclContext, because we go up and the routine that pops
   // the DeclContext assumes that we drill down always.
   // We have to be on the global context. At that point we are in a
   // wrapper function so the parent context must be the global.
   Sema::ContextAndScopeRAII pushedDCAndS(SemaR, C.getTranslationUnitDecl(),
                                          SemaR.TUScope);

   TObject* specObj = gROOT->FindSpecialObject(Name, LookupCtx);
   if (specObj) {
      if (!LookupCtx) {
         Error("GetObjectAddress", "Got a special object without LookupCtx!");
      } else {
         fSpecialObjectMaps[LookupCtx][Name] = specObj;
      }
   }
   return specObj;
}

////////////////////////////////////////////////////////////////////////////////
/// Inject function as a friend into klass.
/// With function being f in void f() {new N::PrivKlass(); } this enables
/// I/O of non-public classes.

void TCling::AddFriendToClass(clang::FunctionDecl* function,
                              clang::CXXRecordDecl* klass) const
{
   using namespace clang;
   ASTContext& Ctx = klass->getASTContext();
   FriendDecl::FriendUnion friendUnion(function);
   // one dummy object for the source location
   SourceLocation sl;
   FriendDecl* friendDecl = FriendDecl::Create(Ctx, klass, sl, friendUnion, sl);
   klass->pushFriendDecl(friendDecl);
}

//______________________________________________________________________________
//
//  DeclId getter.
//

////////////////////////////////////////////////////////////////////////////////
/// Return a unique identifier of the declaration represented by the
/// CallFunc

TInterpreter::DeclId_t TCling::GetDeclId(CallFunc_t* func) const
{
   if (func) return ((TClingCallFunc*)func)->GetDecl()->getCanonicalDecl();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return a (almost) unique identifier of the declaration represented by the
/// ClassInfo.  In ROOT, this identifier can point to more than one TClass
/// when the underlying class is a template instance involving one of the
/// opaque typedef.

TInterpreter::DeclId_t TCling::GetDeclId(ClassInfo_t* cinfo) const
{
   if (cinfo) return ((TClingClassInfo*)cinfo)->GetDeclId();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return a unique identifier of the declaration represented by the
/// MethodInfo

TInterpreter::DeclId_t TCling::GetDeclId(DataMemberInfo_t* data) const
{
   if (data) return ((TClingDataMemberInfo*)data)->GetDeclId();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return a unique identifier of the declaration represented by the
/// MethodInfo

TInterpreter::DeclId_t TCling::GetDeclId(MethodInfo_t* method) const
{
   if (method) return ((TClingMethodInfo*)method)->GetDeclId();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return a unique identifier of the declaration represented by the
/// TypedefInfo

TInterpreter::DeclId_t TCling::GetDeclId(TypedefInfo_t* tinfo) const
{
   if (tinfo) return ((TClingTypedefInfo*)tinfo)->GetDecl()->getCanonicalDecl();
   return 0;
}

//______________________________________________________________________________
//
//  CallFunc interface
//

////////////////////////////////////////////////////////////////////////////////

void TCling::CallFunc_Delete(CallFunc_t* func) const
{
   delete (TClingCallFunc*) func;
}

////////////////////////////////////////////////////////////////////////////////

void TCling::CallFunc_Exec(CallFunc_t* func, void* address) const
{
   TClingCallFunc* f = (TClingCallFunc*) func;
   f->Exec(address);
}

////////////////////////////////////////////////////////////////////////////////

void TCling::CallFunc_Exec(CallFunc_t* func, void* address, TInterpreterValue& val) const
{
   TClingCallFunc* f = (TClingCallFunc*) func;
   f->Exec(address, &val);
}

////////////////////////////////////////////////////////////////////////////////

void TCling::CallFunc_ExecWithReturn(CallFunc_t* func, void* address, void* ret) const
{
   TClingCallFunc* f = (TClingCallFunc*) func;
   f->ExecWithReturn(address, ret);
}

////////////////////////////////////////////////////////////////////////////////

void TCling::CallFunc_ExecWithArgsAndReturn(CallFunc_t* func, void* address,
                                            const void* args[] /*=0*/,
                                            int nargs /*=0*/,
                                            void* ret/*=0*/) const
{
   TClingCallFunc* f = (TClingCallFunc*) func;
   f->ExecWithArgsAndReturn(address, args, nargs, ret);
}

////////////////////////////////////////////////////////////////////////////////

Long_t TCling::CallFunc_ExecInt(CallFunc_t* func, void* address) const
{
   TClingCallFunc* f = (TClingCallFunc*) func;
   return f->ExecInt(address);
}

////////////////////////////////////////////////////////////////////////////////

Long64_t TCling::CallFunc_ExecInt64(CallFunc_t* func, void* address) const
{
   TClingCallFunc* f = (TClingCallFunc*) func;
   return f->ExecInt64(address);
}

////////////////////////////////////////////////////////////////////////////////

Double_t TCling::CallFunc_ExecDouble(CallFunc_t* func, void* address) const
{
   TClingCallFunc* f = (TClingCallFunc*) func;
   return f->ExecDouble(address);
}

////////////////////////////////////////////////////////////////////////////////

CallFunc_t* TCling::CallFunc_Factory() const
{
   R__LOCKGUARD(gInterpreterMutex);
   return (CallFunc_t*) new TClingCallFunc(fInterpreter,*fNormalizedCtxt);
}

////////////////////////////////////////////////////////////////////////////////

CallFunc_t* TCling::CallFunc_FactoryCopy(CallFunc_t* func) const
{
   return (CallFunc_t*) new TClingCallFunc(*(TClingCallFunc*)func);
}

////////////////////////////////////////////////////////////////////////////////

MethodInfo_t* TCling::CallFunc_FactoryMethod(CallFunc_t* func) const
{
   TClingCallFunc* f = (TClingCallFunc*) func;
   return (MethodInfo_t*) f->FactoryMethod();
}

////////////////////////////////////////////////////////////////////////////////

void TCling::CallFunc_IgnoreExtraArgs(CallFunc_t* func, bool ignore) const
{
   TClingCallFunc* f = (TClingCallFunc*) func;
   f->IgnoreExtraArgs(ignore);
}

////////////////////////////////////////////////////////////////////////////////

void TCling::CallFunc_Init(CallFunc_t* func) const
{
   R__LOCKGUARD(gInterpreterMutex);
   TClingCallFunc* f = (TClingCallFunc*) func;
   f->Init();
}

////////////////////////////////////////////////////////////////////////////////

bool TCling::CallFunc_IsValid(CallFunc_t* func) const
{
   TClingCallFunc* f = (TClingCallFunc*) func;
   return f->IsValid();
}

////////////////////////////////////////////////////////////////////////////////

TInterpreter::CallFuncIFacePtr_t
TCling::CallFunc_IFacePtr(CallFunc_t * func) const
{
   TClingCallFunc* f = (TClingCallFunc*) func;
   return f->IFacePtr();
}

////////////////////////////////////////////////////////////////////////////////

void TCling::CallFunc_ResetArg(CallFunc_t* func) const
{
   TClingCallFunc* f = (TClingCallFunc*) func;
   f->ResetArg();
}

////////////////////////////////////////////////////////////////////////////////

void TCling::CallFunc_SetArg(CallFunc_t* func, Long_t param) const
{
   TClingCallFunc* f = (TClingCallFunc*) func;
   f->SetArg(param);
}

////////////////////////////////////////////////////////////////////////////////

void TCling::CallFunc_SetArg(CallFunc_t* func, ULong_t param) const
{
   TClingCallFunc* f = (TClingCallFunc*) func;
   f->SetArg(param);
}

////////////////////////////////////////////////////////////////////////////////

void TCling::CallFunc_SetArg(CallFunc_t* func, Float_t param) const
{
   TClingCallFunc* f = (TClingCallFunc*) func;
   f->SetArg(param);
}

////////////////////////////////////////////////////////////////////////////////

void TCling::CallFunc_SetArg(CallFunc_t* func, Double_t param) const
{
   TClingCallFunc* f = (TClingCallFunc*) func;
   f->SetArg(param);
}

////////////////////////////////////////////////////////////////////////////////

void TCling::CallFunc_SetArg(CallFunc_t* func, Long64_t param) const
{
   TClingCallFunc* f = (TClingCallFunc*) func;
   f->SetArg(param);
}

////////////////////////////////////////////////////////////////////////////////

void TCling::CallFunc_SetArg(CallFunc_t* func, ULong64_t param) const
{
   TClingCallFunc* f = (TClingCallFunc*) func;
   f->SetArg(param);
}

////////////////////////////////////////////////////////////////////////////////

void TCling::CallFunc_SetArgArray(CallFunc_t* func, Long_t* paramArr, Int_t nparam) const
{
   TClingCallFunc* f = (TClingCallFunc*) func;
   f->SetArgArray(paramArr, nparam);
}

////////////////////////////////////////////////////////////////////////////////

void TCling::CallFunc_SetArgs(CallFunc_t* func, const char* param) const
{
   TClingCallFunc* f = (TClingCallFunc*) func;
   f->SetArgs(param);
}

////////////////////////////////////////////////////////////////////////////////

void TCling::CallFunc_SetFunc(CallFunc_t* func, ClassInfo_t* info, const char* method, const char* params, Long_t* offset) const
{
   TClingCallFunc* f = (TClingCallFunc*) func;
   TClingClassInfo* ci = (TClingClassInfo*) info;
   f->SetFunc(ci, method, params, offset);
}

////////////////////////////////////////////////////////////////////////////////

void TCling::CallFunc_SetFunc(CallFunc_t* func, ClassInfo_t* info, const char* method, const char* params, bool objectIsConst, Long_t* offset) const
{
   TClingCallFunc* f = (TClingCallFunc*) func;
   TClingClassInfo* ci = (TClingClassInfo*) info;
   f->SetFunc(ci, method, params, objectIsConst, offset);
}
////////////////////////////////////////////////////////////////////////////////

void TCling::CallFunc_SetFunc(CallFunc_t* func, MethodInfo_t* info) const
{
   TClingCallFunc* f = (TClingCallFunc*) func;
   TClingMethodInfo* minfo = (TClingMethodInfo*) info;
   f->SetFunc(minfo);
}

////////////////////////////////////////////////////////////////////////////////
/// Interface to cling function

void TCling::CallFunc_SetFuncProto(CallFunc_t* func, ClassInfo_t* info, const char* method, const char* proto, Long_t* offset, EFunctionMatchMode mode /* = kConversionMatch */) const
{
   TClingCallFunc* f = (TClingCallFunc*) func;
   TClingClassInfo* ci = (TClingClassInfo*) info;
   f->SetFuncProto(ci, method, proto, offset, mode);
}

////////////////////////////////////////////////////////////////////////////////
/// Interface to cling function

void TCling::CallFunc_SetFuncProto(CallFunc_t* func, ClassInfo_t* info, const char* method, const char* proto, bool objectIsConst, Long_t* offset, EFunctionMatchMode mode /* = kConversionMatch */) const
{
   TClingCallFunc* f = (TClingCallFunc*) func;
   TClingClassInfo* ci = (TClingClassInfo*) info;
   f->SetFuncProto(ci, method, proto, objectIsConst, offset, mode);
}

////////////////////////////////////////////////////////////////////////////////
/// Interface to cling function

void TCling::CallFunc_SetFuncProto(CallFunc_t* func, ClassInfo_t* info, const char* method, const std::vector<TypeInfo_t*> &proto, Long_t* offset, EFunctionMatchMode mode /* = kConversionMatch */) const
{
   TClingCallFunc* f = (TClingCallFunc*) func;
   TClingClassInfo* ci = (TClingClassInfo*) info;
   llvm::SmallVector<clang::QualType, 4> funcProto;
   for (std::vector<TypeInfo_t*>::const_iterator iter = proto.begin(), end = proto.end();
        iter != end; ++iter) {
      funcProto.push_back( ((TClingTypeInfo*)(*iter))->GetQualType() );
   }
   f->SetFuncProto(ci, method, funcProto, offset, mode);
}

////////////////////////////////////////////////////////////////////////////////
/// Interface to cling function

void TCling::CallFunc_SetFuncProto(CallFunc_t* func, ClassInfo_t* info, const char* method, const std::vector<TypeInfo_t*> &proto, bool objectIsConst, Long_t* offset, EFunctionMatchMode mode /* = kConversionMatch */) const
{
   TClingCallFunc* f = (TClingCallFunc*) func;
   TClingClassInfo* ci = (TClingClassInfo*) info;
   llvm::SmallVector<clang::QualType, 4> funcProto;
   for (std::vector<TypeInfo_t*>::const_iterator iter = proto.begin(), end = proto.end();
        iter != end; ++iter) {
      funcProto.push_back( ((TClingTypeInfo*)(*iter))->GetQualType() );
   }
   f->SetFuncProto(ci, method, funcProto, objectIsConst, offset, mode);
}

std::string TCling::CallFunc_GetWrapperCode(CallFunc_t *func) const
{
   TClingCallFunc *f = (TClingCallFunc *)func;
   std::string wrapper_name;
   std::string wrapper;
   f->get_wrapper_code(wrapper_name, wrapper);
   return wrapper;
}

//______________________________________________________________________________
//
//  ClassInfo interface
//

////////////////////////////////////////////////////////////////////////////////
/// Return true if the entity pointed to by 'declid' is declared in
/// the context described by 'info'.  If info is null, look into the
/// global scope (translation unit scope).

Bool_t TCling::ClassInfo_Contains(ClassInfo_t *info, DeclId_t declid) const
{
   if (!declid) return kFALSE;

   const clang::Decl *scope;
   if (info) scope = ((TClingClassInfo*)info)->GetDecl();
   else scope = fInterpreter->getCI()->getASTContext().getTranslationUnitDecl();

   const clang::Decl *decl = reinterpret_cast<const clang::Decl*>(declid);
   const clang::DeclContext *ctxt = clang::Decl::castToDeclContext(scope);
   if (!decl || !ctxt) return kFALSE;
   if (decl->getDeclContext()->Equals(ctxt))
      return kTRUE;
   else if (decl->getDeclContext()->isTransparentContext() &&
            decl->getDeclContext()->getParent()->Equals(ctxt))
      return kTRUE;
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////

Long_t TCling::ClassInfo_ClassProperty(ClassInfo_t* cinfo) const
{
   TClingClassInfo* TClinginfo = (TClingClassInfo*) cinfo;
   return TClinginfo->ClassProperty();
}

////////////////////////////////////////////////////////////////////////////////

void TCling::ClassInfo_Delete(ClassInfo_t* cinfo) const
{
   delete (TClingClassInfo*) cinfo;
}

////////////////////////////////////////////////////////////////////////////////

void TCling::ClassInfo_Delete(ClassInfo_t* cinfo, void* arena) const
{
   TClingClassInfo* TClinginfo = (TClingClassInfo*) cinfo;
   TClinginfo->Delete(arena,*fNormalizedCtxt);
}

////////////////////////////////////////////////////////////////////////////////

void TCling::ClassInfo_DeleteArray(ClassInfo_t* cinfo, void* arena, bool dtorOnly) const
{
   TClingClassInfo* TClinginfo = (TClingClassInfo*) cinfo;
   TClinginfo->DeleteArray(arena, dtorOnly,*fNormalizedCtxt);
}

////////////////////////////////////////////////////////////////////////////////

void TCling::ClassInfo_Destruct(ClassInfo_t* cinfo, void* arena) const
{
   TClingClassInfo* TClinginfo = (TClingClassInfo*) cinfo;
   TClinginfo->Destruct(arena,*fNormalizedCtxt);
}

////////////////////////////////////////////////////////////////////////////////

ClassInfo_t* TCling::ClassInfo_Factory(Bool_t all) const
{
   R__LOCKGUARD(gInterpreterMutex);
   return (ClassInfo_t*) new TClingClassInfo(fInterpreter, all);
}

////////////////////////////////////////////////////////////////////////////////

ClassInfo_t* TCling::ClassInfo_Factory(ClassInfo_t* cinfo) const
{
   return (ClassInfo_t*) new TClingClassInfo(*(TClingClassInfo*)cinfo);
}

////////////////////////////////////////////////////////////////////////////////

ClassInfo_t* TCling::ClassInfo_Factory(const char* name) const
{
   R__LOCKGUARD(gInterpreterMutex);
   return (ClassInfo_t*) new TClingClassInfo(fInterpreter, name);
}

////////////////////////////////////////////////////////////////////////////////

int TCling::ClassInfo_GetMethodNArg(ClassInfo_t* cinfo, const char* method, const char* proto, Bool_t objectIsConst /* = false */, EFunctionMatchMode mode /* = kConversionMatch */) const
{
   TClingClassInfo* TClinginfo = (TClingClassInfo*) cinfo;
   return TClinginfo->GetMethodNArg(method, proto, objectIsConst, mode);
}

////////////////////////////////////////////////////////////////////////////////

bool TCling::ClassInfo_HasDefaultConstructor(ClassInfo_t* cinfo) const
{
   TClingClassInfo* TClinginfo = (TClingClassInfo*) cinfo;
   return TClinginfo->HasDefaultConstructor();
}

////////////////////////////////////////////////////////////////////////////////

bool TCling::ClassInfo_HasMethod(ClassInfo_t* cinfo, const char* name) const
{
   TClingClassInfo* TClinginfo = (TClingClassInfo*) cinfo;
   return TClinginfo->HasMethod(name);
}

////////////////////////////////////////////////////////////////////////////////

void TCling::ClassInfo_Init(ClassInfo_t* cinfo, const char* name) const
{
   R__LOCKGUARD(gInterpreterMutex);
   TClingClassInfo* TClinginfo = (TClingClassInfo*) cinfo;
   TClinginfo->Init(name);
}

////////////////////////////////////////////////////////////////////////////////

void TCling::ClassInfo_Init(ClassInfo_t* cinfo, int tagnum) const
{
   R__LOCKGUARD(gInterpreterMutex);
   TClingClassInfo* TClinginfo = (TClingClassInfo*) cinfo;
   TClinginfo->Init(tagnum);
}

////////////////////////////////////////////////////////////////////////////////

bool TCling::ClassInfo_IsBase(ClassInfo_t* cinfo, const char* name) const
{
   TClingClassInfo* TClinginfo = (TClingClassInfo*) cinfo;
   return TClinginfo->IsBase(name);
}

////////////////////////////////////////////////////////////////////////////////

bool TCling::ClassInfo_IsEnum(const char* name) const
{
   return TClingClassInfo::IsEnum(fInterpreter, name);
}

////////////////////////////////////////////////////////////////////////////////

bool TCling::ClassInfo_IsLoaded(ClassInfo_t* cinfo) const
{
   TClingClassInfo* TClinginfo = (TClingClassInfo*) cinfo;
   return TClinginfo->IsLoaded();
}

////////////////////////////////////////////////////////////////////////////////

bool TCling::ClassInfo_IsValid(ClassInfo_t* cinfo) const
{
   TClingClassInfo* TClinginfo = (TClingClassInfo*) cinfo;
   return TClinginfo->IsValid();
}

////////////////////////////////////////////////////////////////////////////////

bool TCling::ClassInfo_IsValidMethod(ClassInfo_t* cinfo, const char* method, const char* proto, Long_t* offset, EFunctionMatchMode mode /* = kConversionMatch */) const
{
   TClingClassInfo* TClinginfo = (TClingClassInfo*) cinfo;
   return TClinginfo->IsValidMethod(method, proto, false, offset, mode);
}

////////////////////////////////////////////////////////////////////////////////

bool TCling::ClassInfo_IsValidMethod(ClassInfo_t* cinfo, const char* method, const char* proto, Bool_t objectIsConst, Long_t* offset, EFunctionMatchMode mode /* = kConversionMatch */) const
{
   TClingClassInfo* TClinginfo = (TClingClassInfo*) cinfo;
   return TClinginfo->IsValidMethod(method, proto, objectIsConst, offset, mode);
}

////////////////////////////////////////////////////////////////////////////////

int TCling::ClassInfo_Next(ClassInfo_t* cinfo) const
{
   TClingClassInfo* TClinginfo = (TClingClassInfo*) cinfo;
   return TClinginfo->Next();
}

////////////////////////////////////////////////////////////////////////////////

void* TCling::ClassInfo_New(ClassInfo_t* cinfo) const
{
   TClingClassInfo* TClinginfo = (TClingClassInfo*) cinfo;
   return TClinginfo->New(*fNormalizedCtxt);
}

////////////////////////////////////////////////////////////////////////////////

void* TCling::ClassInfo_New(ClassInfo_t* cinfo, int n) const
{
   TClingClassInfo* TClinginfo = (TClingClassInfo*) cinfo;
   return TClinginfo->New(n,*fNormalizedCtxt);
}

////////////////////////////////////////////////////////////////////////////////

void* TCling::ClassInfo_New(ClassInfo_t* cinfo, int n, void* arena) const
{
   TClingClassInfo* TClinginfo = (TClingClassInfo*) cinfo;
   return TClinginfo->New(n, arena,*fNormalizedCtxt);
}

////////////////////////////////////////////////////////////////////////////////

void* TCling::ClassInfo_New(ClassInfo_t* cinfo, void* arena) const
{
   TClingClassInfo* TClinginfo = (TClingClassInfo*) cinfo;
   return TClinginfo->New(arena,*fNormalizedCtxt);
}

////////////////////////////////////////////////////////////////////////////////

Long_t TCling::ClassInfo_Property(ClassInfo_t* cinfo) const
{
   TClingClassInfo* TClinginfo = (TClingClassInfo*) cinfo;
   return TClinginfo->Property();
}

////////////////////////////////////////////////////////////////////////////////

int TCling::ClassInfo_Size(ClassInfo_t* cinfo) const
{
   TClingClassInfo* TClinginfo = (TClingClassInfo*) cinfo;
   return TClinginfo->Size();
}

////////////////////////////////////////////////////////////////////////////////

Long_t TCling::ClassInfo_Tagnum(ClassInfo_t* cinfo) const
{
   TClingClassInfo* TClinginfo = (TClingClassInfo*) cinfo;
   return TClinginfo->Tagnum();
}

////////////////////////////////////////////////////////////////////////////////

const char* TCling::ClassInfo_FileName(ClassInfo_t* cinfo) const
{
   TClingClassInfo* TClinginfo = (TClingClassInfo*) cinfo;
   return TClinginfo->FileName();
}

////////////////////////////////////////////////////////////////////////////////

const char* TCling::ClassInfo_FullName(ClassInfo_t* cinfo) const
{
   TClingClassInfo* TClinginfo = (TClingClassInfo*) cinfo;
   TTHREAD_TLS_DECL(std::string,output);
   TClinginfo->FullName(output,*fNormalizedCtxt);
   return output.c_str();
}

////////////////////////////////////////////////////////////////////////////////

const char* TCling::ClassInfo_Name(ClassInfo_t* cinfo) const
{
   TClingClassInfo* TClinginfo = (TClingClassInfo*) cinfo;
   return TClinginfo->Name();
}

////////////////////////////////////////////////////////////////////////////////

const char* TCling::ClassInfo_Title(ClassInfo_t* cinfo) const
{
   TClingClassInfo* TClinginfo = (TClingClassInfo*) cinfo;
   return TClinginfo->Title();
}

////////////////////////////////////////////////////////////////////////////////

const char* TCling::ClassInfo_TmpltName(ClassInfo_t* cinfo) const
{
   TClingClassInfo* TClinginfo = (TClingClassInfo*) cinfo;
   return TClinginfo->TmpltName();
}



//______________________________________________________________________________
//
//  BaseClassInfo interface
//

////////////////////////////////////////////////////////////////////////////////

void TCling::BaseClassInfo_Delete(BaseClassInfo_t* bcinfo) const
{
   delete(TClingBaseClassInfo*) bcinfo;
}

////////////////////////////////////////////////////////////////////////////////

BaseClassInfo_t* TCling::BaseClassInfo_Factory(ClassInfo_t* cinfo) const
{
   R__LOCKGUARD(gInterpreterMutex);
   TClingClassInfo* TClinginfo = (TClingClassInfo*) cinfo;
   return (BaseClassInfo_t*) new TClingBaseClassInfo(fInterpreter, TClinginfo);
}

////////////////////////////////////////////////////////////////////////////////

BaseClassInfo_t* TCling::BaseClassInfo_Factory(ClassInfo_t* derived,
   ClassInfo_t* base) const
{
   R__LOCKGUARD(gInterpreterMutex);
   TClingClassInfo* TClinginfo = (TClingClassInfo*) derived;
   TClingClassInfo* TClinginfoBase = (TClingClassInfo*) base;
   return (BaseClassInfo_t*) new TClingBaseClassInfo(fInterpreter, TClinginfo, TClinginfoBase);
}

////////////////////////////////////////////////////////////////////////////////

int TCling::BaseClassInfo_Next(BaseClassInfo_t* bcinfo) const
{
   TClingBaseClassInfo* TClinginfo = (TClingBaseClassInfo*) bcinfo;
   return TClinginfo->Next();
}

////////////////////////////////////////////////////////////////////////////////

int TCling::BaseClassInfo_Next(BaseClassInfo_t* bcinfo, int onlyDirect) const
{
   TClingBaseClassInfo* TClinginfo = (TClingBaseClassInfo*) bcinfo;
   return TClinginfo->Next(onlyDirect);
}

////////////////////////////////////////////////////////////////////////////////

Long_t TCling::BaseClassInfo_Offset(BaseClassInfo_t* toBaseClassInfo, void * address, bool isDerivedObject) const
{
   TClingBaseClassInfo* TClinginfo = (TClingBaseClassInfo*) toBaseClassInfo;
   return TClinginfo->Offset(address, isDerivedObject);
}

////////////////////////////////////////////////////////////////////////////////

Long_t TCling::ClassInfo_GetBaseOffset(ClassInfo_t* fromDerived, ClassInfo_t* toBase, void * address, bool isDerivedObject) const
{
   TClingClassInfo* TClinginfo = (TClingClassInfo*) fromDerived;
   TClingClassInfo* TClinginfoBase = (TClingClassInfo*) toBase;
   // Offset to the class itself.
   if (TClinginfo->GetDecl() == TClinginfoBase->GetDecl()) {
      return 0;
   }
   return TClinginfo->GetBaseOffset(TClinginfoBase, address, isDerivedObject);
}

////////////////////////////////////////////////////////////////////////////////

Long_t TCling::BaseClassInfo_Property(BaseClassInfo_t* bcinfo) const
{
   TClingBaseClassInfo* TClinginfo = (TClingBaseClassInfo*) bcinfo;
   return TClinginfo->Property();
}

////////////////////////////////////////////////////////////////////////////////

ClassInfo_t *TCling::BaseClassInfo_ClassInfo(BaseClassInfo_t *bcinfo) const
{
   TClingBaseClassInfo* TClinginfo = (TClingBaseClassInfo*) bcinfo;
   return (ClassInfo_t *)TClinginfo->GetBase();
}

////////////////////////////////////////////////////////////////////////////////

Long_t TCling::BaseClassInfo_Tagnum(BaseClassInfo_t* bcinfo) const
{
   TClingBaseClassInfo* TClinginfo = (TClingBaseClassInfo*) bcinfo;
   return TClinginfo->Tagnum();
}

////////////////////////////////////////////////////////////////////////////////

const char* TCling::BaseClassInfo_FullName(BaseClassInfo_t* bcinfo) const
{
   TClingBaseClassInfo* TClinginfo = (TClingBaseClassInfo*) bcinfo;
   TTHREAD_TLS_DECL(std::string,output);
   TClinginfo->FullName(output,*fNormalizedCtxt);
   return output.c_str();
}

////////////////////////////////////////////////////////////////////////////////

const char* TCling::BaseClassInfo_Name(BaseClassInfo_t* bcinfo) const
{
   TClingBaseClassInfo* TClinginfo = (TClingBaseClassInfo*) bcinfo;
   return TClinginfo->Name();
}

////////////////////////////////////////////////////////////////////////////////

const char* TCling::BaseClassInfo_TmpltName(BaseClassInfo_t* bcinfo) const
{
   TClingBaseClassInfo* TClinginfo = (TClingBaseClassInfo*) bcinfo;
   return TClinginfo->TmpltName();
}

//______________________________________________________________________________
//
//  DataMemberInfo interface
//

////////////////////////////////////////////////////////////////////////////////

int TCling::DataMemberInfo_ArrayDim(DataMemberInfo_t* dminfo) const
{
   TClingDataMemberInfo* TClinginfo = (TClingDataMemberInfo*) dminfo;
   return TClinginfo->ArrayDim();
}

////////////////////////////////////////////////////////////////////////////////

void TCling::DataMemberInfo_Delete(DataMemberInfo_t* dminfo) const
{
   delete(TClingDataMemberInfo*) dminfo;
}

////////////////////////////////////////////////////////////////////////////////

DataMemberInfo_t* TCling::DataMemberInfo_Factory(ClassInfo_t* clinfo /*= 0*/) const
{
   R__LOCKGUARD(gInterpreterMutex);
   TClingClassInfo* TClingclass_info = (TClingClassInfo*) clinfo;
   return (DataMemberInfo_t*) new TClingDataMemberInfo(fInterpreter, TClingclass_info);
}

////////////////////////////////////////////////////////////////////////////////

DataMemberInfo_t* TCling::DataMemberInfo_Factory(DeclId_t declid, ClassInfo_t* clinfo) const
{
   R__LOCKGUARD(gInterpreterMutex);
   const clang::Decl* decl = reinterpret_cast<const clang::Decl*>(declid);
   const clang::ValueDecl* vd = llvm::dyn_cast_or_null<clang::ValueDecl>(decl);
   return (DataMemberInfo_t*) new TClingDataMemberInfo(fInterpreter, vd, (TClingClassInfo*)clinfo);
}

////////////////////////////////////////////////////////////////////////////////

DataMemberInfo_t* TCling::DataMemberInfo_FactoryCopy(DataMemberInfo_t* dminfo) const
{
   TClingDataMemberInfo* TClinginfo = (TClingDataMemberInfo*) dminfo;
   return (DataMemberInfo_t*) new TClingDataMemberInfo(*TClinginfo);
}

////////////////////////////////////////////////////////////////////////////////

bool TCling::DataMemberInfo_IsValid(DataMemberInfo_t* dminfo) const
{
   TClingDataMemberInfo* TClinginfo = (TClingDataMemberInfo*) dminfo;
   return TClinginfo->IsValid();
}

////////////////////////////////////////////////////////////////////////////////

int TCling::DataMemberInfo_MaxIndex(DataMemberInfo_t* dminfo, Int_t dim) const
{
   TClingDataMemberInfo* TClinginfo = (TClingDataMemberInfo*) dminfo;
   return TClinginfo->MaxIndex(dim);
}

////////////////////////////////////////////////////////////////////////////////

int TCling::DataMemberInfo_Next(DataMemberInfo_t* dminfo) const
{
   TClingDataMemberInfo* TClinginfo = (TClingDataMemberInfo*) dminfo;
   return TClinginfo->Next();
}

////////////////////////////////////////////////////////////////////////////////

Long_t TCling::DataMemberInfo_Offset(DataMemberInfo_t* dminfo) const
{
   TClingDataMemberInfo* TClinginfo = (TClingDataMemberInfo*) dminfo;
   return TClinginfo->Offset();
}

////////////////////////////////////////////////////////////////////////////////

Long_t TCling::DataMemberInfo_Property(DataMemberInfo_t* dminfo) const
{
   TClingDataMemberInfo* TClinginfo = (TClingDataMemberInfo*) dminfo;
   return TClinginfo->Property();
}

////////////////////////////////////////////////////////////////////////////////

Long_t TCling::DataMemberInfo_TypeProperty(DataMemberInfo_t* dminfo) const
{
   TClingDataMemberInfo* TClinginfo = (TClingDataMemberInfo*) dminfo;
   return TClinginfo->TypeProperty();
}

////////////////////////////////////////////////////////////////////////////////

int TCling::DataMemberInfo_TypeSize(DataMemberInfo_t* dminfo) const
{
   TClingDataMemberInfo* TClinginfo = (TClingDataMemberInfo*) dminfo;
   return TClinginfo->TypeSize();
}

////////////////////////////////////////////////////////////////////////////////

const char* TCling::DataMemberInfo_TypeName(DataMemberInfo_t* dminfo) const
{
   TClingDataMemberInfo* TClinginfo = (TClingDataMemberInfo*) dminfo;
   return TClinginfo->TypeName();
}

////////////////////////////////////////////////////////////////////////////////

const char* TCling::DataMemberInfo_TypeTrueName(DataMemberInfo_t* dminfo) const
{
   TClingDataMemberInfo* TClinginfo = (TClingDataMemberInfo*) dminfo;
   return TClinginfo->TypeTrueName(*fNormalizedCtxt);
}

////////////////////////////////////////////////////////////////////////////////

const char* TCling::DataMemberInfo_Name(DataMemberInfo_t* dminfo) const
{
   TClingDataMemberInfo* TClinginfo = (TClingDataMemberInfo*) dminfo;
   return TClinginfo->Name();
}

////////////////////////////////////////////////////////////////////////////////

const char* TCling::DataMemberInfo_Title(DataMemberInfo_t* dminfo) const
{
   TClingDataMemberInfo* TClinginfo = (TClingDataMemberInfo*) dminfo;
   return TClinginfo->Title();
}

////////////////////////////////////////////////////////////////////////////////

const char* TCling::DataMemberInfo_ValidArrayIndex(DataMemberInfo_t* dminfo) const
{
   TTHREAD_TLS_DECL(std::string,result);

   TClingDataMemberInfo* TClinginfo = (TClingDataMemberInfo*) dminfo;
   result = TClinginfo->ValidArrayIndex().str();
   return result.c_str();
}

////////////////////////////////////////////////////////////////////////////////

void TCling::SetDeclAttr(DeclId_t declId, const char* attribute)
{
   Decl* decl = static_cast<Decl*>(const_cast<void*>(declId));
   ASTContext &C = decl->getASTContext();
   SourceRange commentRange; // this is a fake comment range
   decl->addAttr( new (C) AnnotateAttr( commentRange, C, attribute, 0 ) );
}

//______________________________________________________________________________
//
// Function Template interface
//

////////////////////////////////////////////////////////////////////////////////

static void ConstructorName(std::string &name, const clang::NamedDecl *decl,
                            cling::Interpreter &interp,
                            const ROOT::TMetaUtils::TNormalizedCtxt &normCtxt)
{
   const clang::TypeDecl* td = llvm::dyn_cast<clang::TypeDecl>(decl->getDeclContext());
   if (!td) return;

   clang::QualType qualType(td->getTypeForDecl(),0);
   ROOT::TMetaUtils::GetNormalizedName(name, qualType, interp, normCtxt);
   unsigned int level = 0;
   for(size_t cursor = name.length()-1; cursor != 0; --cursor) {
      if (name[cursor] == '>') ++level;
      else if (name[cursor] == '<' && level) --level;
      else if (level == 0 && name[cursor] == ':') {
         name.erase(0,cursor+1);
         break;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////

void TCling::GetFunctionName(const clang::FunctionDecl *decl, std::string &output) const
{
   output.clear();
   if (llvm::isa<clang::CXXConstructorDecl>(decl))
   {
      ConstructorName(output, decl, *fInterpreter, *fNormalizedCtxt);

   } else if (llvm::isa<clang::CXXDestructorDecl>(decl))
   {
      ConstructorName(output, decl, *fInterpreter, *fNormalizedCtxt);
      output.insert(output.begin(), '~');
   } else {
      llvm::raw_string_ostream stream(output);
      auto printPolicy = decl->getASTContext().getPrintingPolicy();
      // Don't trigger fopen of the source file to count lines:
      printPolicy.AnonymousTagLocations = false;
      decl->getNameForDiagnostic(stream, printPolicy, /*Qualified=*/false);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return a unique identifier of the declaration represented by the
/// FuncTempInfo

TInterpreter::DeclId_t TCling::GetDeclId(FuncTempInfo_t *info) const
{
   return (DeclId_t)info;
}

////////////////////////////////////////////////////////////////////////////////
/// Delete the FuncTempInfo_t

void   TCling::FuncTempInfo_Delete(FuncTempInfo_t * /* ft_info */) const
{
   // Currently the address of ft_info is actually the decl itself,
   // so we have nothing to do.
}

////////////////////////////////////////////////////////////////////////////////
/// Construct a FuncTempInfo_t

FuncTempInfo_t *TCling::FuncTempInfo_Factory(DeclId_t declid) const
{
   // Currently the address of ft_info is actually the decl itself,
   // so we have nothing to do.

   return (FuncTempInfo_t*)const_cast<void*>(declid);
}

////////////////////////////////////////////////////////////////////////////////
/// Construct a FuncTempInfo_t

FuncTempInfo_t *TCling::FuncTempInfo_FactoryCopy(FuncTempInfo_t *ft_info) const
{
   // Currently the address of ft_info is actually the decl itself,
   // so we have nothing to do.

   return (FuncTempInfo_t*)ft_info;
}

////////////////////////////////////////////////////////////////////////////////
/// Check validity of a FuncTempInfo_t

Bool_t TCling::FuncTempInfo_IsValid(FuncTempInfo_t *t_info) const
{
   // Currently the address of ft_info is actually the decl itself,
   // so we have nothing to do.

   return t_info != 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the maximum number of template arguments of the
/// function template described by ft_info.

UInt_t TCling::FuncTempInfo_TemplateNargs(FuncTempInfo_t *ft_info) const
{
   if (!ft_info) return 0;
   const clang::FunctionTemplateDecl *ft = (const clang::FunctionTemplateDecl*)ft_info;
   return ft->getTemplateParameters()->size();
}

////////////////////////////////////////////////////////////////////////////////
/// Return the number of required template arguments of the
/// function template described by ft_info.

UInt_t TCling::FuncTempInfo_TemplateMinReqArgs(FuncTempInfo_t *ft_info) const
{
   if (!ft_info) return 0;
   const clang::FunctionTemplateDecl *ft = (clang::FunctionTemplateDecl*)ft_info;
   return ft->getTemplateParameters()->getMinRequiredArguments();
}

////////////////////////////////////////////////////////////////////////////////
/// Return the property of the function template.

Long_t TCling::FuncTempInfo_Property(FuncTempInfo_t *ft_info) const
{
   if (!ft_info) return 0;

   long property = 0L;
   property |= kIsCompiled;

   const clang::FunctionTemplateDecl *ft = (clang::FunctionTemplateDecl*)ft_info;

   switch (ft->getAccess()) {
      case clang::AS_public:
         property |= kIsPublic;
         break;
      case clang::AS_protected:
         property |= kIsProtected;
         break;
      case clang::AS_private:
         property |= kIsPrivate;
         break;
      case clang::AS_none:
         if (ft->getDeclContext()->isNamespace())
            property |= kIsPublic;
         break;
      default:
         // IMPOSSIBLE
         break;
   }

   const clang::FunctionDecl *fd = ft->getTemplatedDecl();
   if (const clang::CXXMethodDecl *md =
       llvm::dyn_cast<clang::CXXMethodDecl>(fd)) {
      if (md->getTypeQualifiers() & clang::Qualifiers::Const) {
         property |= kIsConstant | kIsConstMethod;
      }
      if (md->isVirtual()) {
         property |= kIsVirtual;
      }
      if (md->isPure()) {
         property |= kIsPureVirtual;
      }
      if (const clang::CXXConstructorDecl *cd =
          llvm::dyn_cast<clang::CXXConstructorDecl>(md)) {
         if (cd->isExplicit()) {
            property |= kIsExplicit;
         }
      }
      else if (const clang::CXXConversionDecl *cd =
               llvm::dyn_cast<clang::CXXConversionDecl>(md)) {
         if (cd->isExplicit()) {
            property |= kIsExplicit;
         }
      }
   }
   return property;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the name of this function template.

void TCling::FuncTempInfo_Name(FuncTempInfo_t *ft_info, TString &output) const
{
   output.Clear();
   if (!ft_info) return;
   const clang::FunctionTemplateDecl *ft = (clang::FunctionTemplateDecl*)ft_info;
   std::string buf;
   GetFunctionName(ft->getTemplatedDecl(), buf);
   output = buf;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the comments associates with this function template.

void TCling::FuncTempInfo_Title(FuncTempInfo_t *ft_info, TString &output) const
{
   output.Clear();
   if (!ft_info) return;
   const clang::FunctionTemplateDecl *ft = (const clang::FunctionTemplateDecl*)ft_info;

   // Iterate over the redeclarations, we can have multiple definitions in the
   // redecl chain (came from merging of pcms).
   if (const RedeclarableTemplateDecl *AnnotFD
       = ROOT::TMetaUtils::GetAnnotatedRedeclarable((const RedeclarableTemplateDecl*)ft)) {
      if (AnnotateAttr *A = AnnotFD->getAttr<AnnotateAttr>()) {
         output = A->getAnnotation().str();
         return;
      }
   }
   if (!ft->isFromASTFile()) {
      // Try to get the comment from the header file if present
      // but not for decls from AST file, where rootcling would have
      // created an annotation
      output = ROOT::TMetaUtils::GetComment(*ft).str();
   }
}


//______________________________________________________________________________
//
//  MethodInfo interface
//

////////////////////////////////////////////////////////////////////////////////
/// Interface to cling function

void TCling::MethodInfo_Delete(MethodInfo_t* minfo) const
{
   delete(TClingMethodInfo*) minfo;
}

////////////////////////////////////////////////////////////////////////////////

void TCling::MethodInfo_CreateSignature(MethodInfo_t* minfo, TString& signature) const
{
   TClingMethodInfo* info = (TClingMethodInfo*) minfo;
   info->CreateSignature(signature);
}

////////////////////////////////////////////////////////////////////////////////

MethodInfo_t* TCling::MethodInfo_Factory() const
{
   R__LOCKGUARD(gInterpreterMutex);
   return (MethodInfo_t*) new TClingMethodInfo(fInterpreter);
}

////////////////////////////////////////////////////////////////////////////////

MethodInfo_t* TCling::MethodInfo_Factory(ClassInfo_t* clinfo) const
{
   R__LOCKGUARD(gInterpreterMutex);
   return (MethodInfo_t*) new TClingMethodInfo(fInterpreter, (TClingClassInfo*)clinfo);
}

////////////////////////////////////////////////////////////////////////////////

MethodInfo_t* TCling::MethodInfo_Factory(DeclId_t declid) const
{
   const clang::Decl* decl = reinterpret_cast<const clang::Decl*>(declid);
   R__LOCKGUARD(gInterpreterMutex);
   const clang::FunctionDecl* fd = llvm::dyn_cast_or_null<clang::FunctionDecl>(decl);
   return (MethodInfo_t*) new TClingMethodInfo(fInterpreter, fd);
}

////////////////////////////////////////////////////////////////////////////////

MethodInfo_t* TCling::MethodInfo_FactoryCopy(MethodInfo_t* minfo) const
{
   return (MethodInfo_t*) new TClingMethodInfo(*(TClingMethodInfo*)minfo);
}

////////////////////////////////////////////////////////////////////////////////

void* TCling::MethodInfo_InterfaceMethod(MethodInfo_t* minfo) const
{
   TClingMethodInfo* info = (TClingMethodInfo*) minfo;
   return info->InterfaceMethod(*fNormalizedCtxt);
}

////////////////////////////////////////////////////////////////////////////////

bool TCling::MethodInfo_IsValid(MethodInfo_t* minfo) const
{
   TClingMethodInfo* info = (TClingMethodInfo*) minfo;
   return info->IsValid();
}

////////////////////////////////////////////////////////////////////////////////

int TCling::MethodInfo_NArg(MethodInfo_t* minfo) const
{
   TClingMethodInfo* info = (TClingMethodInfo*) minfo;
   return info->NArg();
}

////////////////////////////////////////////////////////////////////////////////

int TCling::MethodInfo_NDefaultArg(MethodInfo_t* minfo) const
{
   TClingMethodInfo* info = (TClingMethodInfo*) minfo;
   return info->NDefaultArg();
}

////////////////////////////////////////////////////////////////////////////////

int TCling::MethodInfo_Next(MethodInfo_t* minfo) const
{
   TClingMethodInfo* info = (TClingMethodInfo*) minfo;
   return info->Next();
}

////////////////////////////////////////////////////////////////////////////////

Long_t TCling::MethodInfo_Property(MethodInfo_t* minfo) const
{
   TClingMethodInfo* info = (TClingMethodInfo*) minfo;
   return info->Property();
}

////////////////////////////////////////////////////////////////////////////////

Long_t TCling::MethodInfo_ExtraProperty(MethodInfo_t* minfo) const
{
   TClingMethodInfo* info = (TClingMethodInfo*) minfo;
   return info->ExtraProperty();
}

////////////////////////////////////////////////////////////////////////////////

TypeInfo_t* TCling::MethodInfo_Type(MethodInfo_t* minfo) const
{
   TClingMethodInfo* info = (TClingMethodInfo*) minfo;
   return (TypeInfo_t*)info->Type();
}

////////////////////////////////////////////////////////////////////////////////

const char* TCling::MethodInfo_GetMangledName(MethodInfo_t* minfo) const
{
   TClingMethodInfo* info = (TClingMethodInfo*) minfo;
   TTHREAD_TLS_DECL(TString, mangled_name);
   mangled_name = info->GetMangledName();
   return mangled_name;
}

////////////////////////////////////////////////////////////////////////////////

const char* TCling::MethodInfo_GetPrototype(MethodInfo_t* minfo) const
{
   TClingMethodInfo* info = (TClingMethodInfo*) minfo;
   return info->GetPrototype(*fNormalizedCtxt);
}

////////////////////////////////////////////////////////////////////////////////

const char* TCling::MethodInfo_Name(MethodInfo_t* minfo) const
{
   TClingMethodInfo* info = (TClingMethodInfo*) minfo;
   return info->Name(*fNormalizedCtxt);
}

////////////////////////////////////////////////////////////////////////////////

const char* TCling::MethodInfo_TypeName(MethodInfo_t* minfo) const
{
   TClingMethodInfo* info = (TClingMethodInfo*) minfo;
   return info->TypeName();
}

////////////////////////////////////////////////////////////////////////////////

std::string TCling::MethodInfo_TypeNormalizedName(MethodInfo_t* minfo) const
{
   TClingMethodInfo* info = (TClingMethodInfo*) minfo;
   if (info && info->IsValid())
      return info->Type()->NormalizedName(*fNormalizedCtxt);
   else
      return "";
}

////////////////////////////////////////////////////////////////////////////////

const char* TCling::MethodInfo_Title(MethodInfo_t* minfo) const
{
   TClingMethodInfo* info = (TClingMethodInfo*) minfo;
   return info->Title();
}

////////////////////////////////////////////////////////////////////////////////

auto TCling::MethodCallReturnType(TFunction *func) const -> EReturnType
{
   if (func) {
      return MethodInfo_MethodCallReturnType(func->fInfo);
   } else {
      return EReturnType::kOther;
   }
}

////////////////////////////////////////////////////////////////////////////////

auto TCling::MethodInfo_MethodCallReturnType(MethodInfo_t* minfo) const -> EReturnType
{
   TClingMethodInfo* info = (TClingMethodInfo*) minfo;
   if (info && info->IsValid()) {
      TClingTypeInfo *typeinfo = info->Type();
      clang::QualType QT( typeinfo->GetQualType().getCanonicalType() );
      if (QT->isEnumeralType()) {
         return EReturnType::kLong;
      } else if (QT->isPointerType()) {
         // Look for char*
         QT = llvm::cast<clang::PointerType>(QT)->getPointeeType();
         if ( QT->isCharType() ) {
            return EReturnType::kString;
         } else {
            return EReturnType::kOther;
         }
      } else if ( QT->isFloatingType() ) {
         int sz = typeinfo->Size();
         if (sz == 4 || sz == 8) {
            // Support only float and double.
            return EReturnType::kDouble;
         } else {
            return EReturnType::kOther;
         }
      } else if ( QT->isIntegerType() ) {
         int sz = typeinfo->Size();
         if (sz <= 8) {
            // Support only up to long long ... but
            // FIXME the TMethodCall::Execute only
            // return long (4 bytes) ...
            // The v5 implementation of TMethodCall::ReturnType
            // was not making the distinction so we let it go
            // as is for now, but we really need to upgrade
            // TMethodCall::Execute ...
            return EReturnType::kLong;
         } else {
            return EReturnType::kOther;
         }
      } else {
         return EReturnType::kOther;
      }
   } else {
      return EReturnType::kOther;
   }
}

//______________________________________________________________________________
//
//  MethodArgInfo interface
//

////////////////////////////////////////////////////////////////////////////////

void TCling::MethodArgInfo_Delete(MethodArgInfo_t* marginfo) const
{
   delete(TClingMethodArgInfo*) marginfo;
}

////////////////////////////////////////////////////////////////////////////////

MethodArgInfo_t* TCling::MethodArgInfo_Factory() const
{
   R__LOCKGUARD(gInterpreterMutex);
   return (MethodArgInfo_t*) new TClingMethodArgInfo(fInterpreter);
}

////////////////////////////////////////////////////////////////////////////////

MethodArgInfo_t* TCling::MethodArgInfo_Factory(MethodInfo_t *minfo) const
{
   R__LOCKGUARD(gInterpreterMutex);
   return (MethodArgInfo_t*) new TClingMethodArgInfo(fInterpreter, (TClingMethodInfo*)minfo);
}

////////////////////////////////////////////////////////////////////////////////

MethodArgInfo_t* TCling::MethodArgInfo_FactoryCopy(MethodArgInfo_t* marginfo) const
{
   return (MethodArgInfo_t*)
          new TClingMethodArgInfo(*(TClingMethodArgInfo*)marginfo);
}

////////////////////////////////////////////////////////////////////////////////

bool TCling::MethodArgInfo_IsValid(MethodArgInfo_t* marginfo) const
{
   TClingMethodArgInfo* info = (TClingMethodArgInfo*) marginfo;
   return info->IsValid();
}

////////////////////////////////////////////////////////////////////////////////

int TCling::MethodArgInfo_Next(MethodArgInfo_t* marginfo) const
{
   TClingMethodArgInfo* info = (TClingMethodArgInfo*) marginfo;
   return info->Next();
}

////////////////////////////////////////////////////////////////////////////////

Long_t TCling::MethodArgInfo_Property(MethodArgInfo_t* marginfo) const
{
   TClingMethodArgInfo* info = (TClingMethodArgInfo*) marginfo;
   return info->Property();
}

////////////////////////////////////////////////////////////////////////////////

const char* TCling::MethodArgInfo_DefaultValue(MethodArgInfo_t* marginfo) const
{
   TClingMethodArgInfo* info = (TClingMethodArgInfo*) marginfo;
   return info->DefaultValue();
}

////////////////////////////////////////////////////////////////////////////////

const char* TCling::MethodArgInfo_Name(MethodArgInfo_t* marginfo) const
{
   TClingMethodArgInfo* info = (TClingMethodArgInfo*) marginfo;
   return info->Name();
}

////////////////////////////////////////////////////////////////////////////////

const char* TCling::MethodArgInfo_TypeName(MethodArgInfo_t* marginfo) const
{
   TClingMethodArgInfo* info = (TClingMethodArgInfo*) marginfo;
   return info->TypeName();
}

////////////////////////////////////////////////////////////////////////////////

std::string TCling::MethodArgInfo_TypeNormalizedName(MethodArgInfo_t* marginfo) const
{
   TClingMethodArgInfo* info = (TClingMethodArgInfo*) marginfo;
   return info->Type()->NormalizedName(*fNormalizedCtxt);
}

//______________________________________________________________________________
//
//  TypeInfo interface
//

////////////////////////////////////////////////////////////////////////////////

void TCling::TypeInfo_Delete(TypeInfo_t* tinfo) const
{
   delete (TClingTypeInfo*) tinfo;
}

////////////////////////////////////////////////////////////////////////////////

TypeInfo_t* TCling::TypeInfo_Factory() const
{
   R__LOCKGUARD(gInterpreterMutex);
   return (TypeInfo_t*) new TClingTypeInfo(fInterpreter);
}

////////////////////////////////////////////////////////////////////////////////

TypeInfo_t* TCling::TypeInfo_Factory(const char *name) const
{
   R__LOCKGUARD(gInterpreterMutex);
   return (TypeInfo_t*) new TClingTypeInfo(fInterpreter, name);
}

////////////////////////////////////////////////////////////////////////////////

TypeInfo_t* TCling::TypeInfo_FactoryCopy(TypeInfo_t* tinfo) const
{
   return (TypeInfo_t*) new TClingTypeInfo(*(TClingTypeInfo*)tinfo);
}

////////////////////////////////////////////////////////////////////////////////

void TCling::TypeInfo_Init(TypeInfo_t* tinfo, const char* name) const
{
   R__LOCKGUARD(gInterpreterMutex);
   TClingTypeInfo* TClinginfo = (TClingTypeInfo*) tinfo;
   TClinginfo->Init(name);
}

////////////////////////////////////////////////////////////////////////////////

bool TCling::TypeInfo_IsValid(TypeInfo_t* tinfo) const
{
   TClingTypeInfo* TClinginfo = (TClingTypeInfo*) tinfo;
   return TClinginfo->IsValid();
}

////////////////////////////////////////////////////////////////////////////////

const char* TCling::TypeInfo_Name(TypeInfo_t* tinfo) const
{
   TClingTypeInfo* TClinginfo = (TClingTypeInfo*) tinfo;
   return TClinginfo->Name();
}

////////////////////////////////////////////////////////////////////////////////

Long_t TCling::TypeInfo_Property(TypeInfo_t* tinfo) const
{
   TClingTypeInfo* TClinginfo = (TClingTypeInfo*) tinfo;
   return TClinginfo->Property();
}

////////////////////////////////////////////////////////////////////////////////

int TCling::TypeInfo_RefType(TypeInfo_t* tinfo) const
{
   TClingTypeInfo* TClinginfo = (TClingTypeInfo*) tinfo;
   return TClinginfo->RefType();
}

////////////////////////////////////////////////////////////////////////////////

int TCling::TypeInfo_Size(TypeInfo_t* tinfo) const
{
   TClingTypeInfo* TClinginfo = (TClingTypeInfo*) tinfo;
   return TClinginfo->Size();
}

////////////////////////////////////////////////////////////////////////////////

const char* TCling::TypeInfo_TrueName(TypeInfo_t* tinfo) const
{
   TClingTypeInfo* TClinginfo = (TClingTypeInfo*) tinfo;
   return TClinginfo->TrueName(*fNormalizedCtxt);
}


//______________________________________________________________________________
//
//  TypedefInfo interface
//

////////////////////////////////////////////////////////////////////////////////

void TCling::TypedefInfo_Delete(TypedefInfo_t* tinfo) const
{
   delete(TClingTypedefInfo*) tinfo;
}

////////////////////////////////////////////////////////////////////////////////

TypedefInfo_t* TCling::TypedefInfo_Factory() const
{
   R__LOCKGUARD(gInterpreterMutex);
   return (TypedefInfo_t*) new TClingTypedefInfo(fInterpreter);
}

////////////////////////////////////////////////////////////////////////////////

TypedefInfo_t* TCling::TypedefInfo_Factory(const char *name) const
{
   R__LOCKGUARD(gInterpreterMutex);
   return (TypedefInfo_t*) new TClingTypedefInfo(fInterpreter, name);
}

////////////////////////////////////////////////////////////////////////////////

TypedefInfo_t* TCling::TypedefInfo_FactoryCopy(TypedefInfo_t* tinfo) const
{
   return (TypedefInfo_t*) new TClingTypedefInfo(*(TClingTypedefInfo*)tinfo);
}

////////////////////////////////////////////////////////////////////////////////

void TCling::TypedefInfo_Init(TypedefInfo_t* tinfo,
                              const char* name) const
{
   R__LOCKGUARD(gInterpreterMutex);
   TClingTypedefInfo* TClinginfo = (TClingTypedefInfo*) tinfo;
   TClinginfo->Init(name);
}

////////////////////////////////////////////////////////////////////////////////

bool TCling::TypedefInfo_IsValid(TypedefInfo_t* tinfo) const
{
   TClingTypedefInfo* TClinginfo = (TClingTypedefInfo*) tinfo;
   return TClinginfo->IsValid();
}

////////////////////////////////////////////////////////////////////////////////

Int_t TCling::TypedefInfo_Next(TypedefInfo_t* tinfo) const
{
   TClingTypedefInfo* TClinginfo = (TClingTypedefInfo*) tinfo;
   return TClinginfo->Next();
}

////////////////////////////////////////////////////////////////////////////////

Long_t TCling::TypedefInfo_Property(TypedefInfo_t* tinfo) const
{
   TClingTypedefInfo* TClinginfo = (TClingTypedefInfo*) tinfo;
   return TClinginfo->Property();
}

////////////////////////////////////////////////////////////////////////////////

int TCling::TypedefInfo_Size(TypedefInfo_t* tinfo) const
{
   TClingTypedefInfo* TClinginfo = (TClingTypedefInfo*) tinfo;
   return TClinginfo->Size();
}

////////////////////////////////////////////////////////////////////////////////

const char* TCling::TypedefInfo_TrueName(TypedefInfo_t* tinfo) const
{
   TClingTypedefInfo* TClinginfo = (TClingTypedefInfo*) tinfo;
   return TClinginfo->TrueName(*fNormalizedCtxt);
}

////////////////////////////////////////////////////////////////////////////////

const char* TCling::TypedefInfo_Name(TypedefInfo_t* tinfo) const
{
   TClingTypedefInfo* TClinginfo = (TClingTypedefInfo*) tinfo;
   return TClinginfo->Name();
}

////////////////////////////////////////////////////////////////////////////////

const char* TCling::TypedefInfo_Title(TypedefInfo_t* tinfo) const
{
   TClingTypedefInfo* TClinginfo = (TClingTypedefInfo*) tinfo;
   return TClinginfo->Title();
}

////////////////////////////////////////////////////////////////////////////////

void TCling::SnapshotMutexState(ROOT::TVirtualRWMutex* mtx)
{
   if (!fInitialMutex.back()) {
      if (fInitialMutex.back().fRecurseCount) {
         Error("SnapshotMutexState", "fRecurseCount != 0 even though initial mutex state is unset!");
      }
      fInitialMutex.back().fState = mtx->GetStateBefore();
   }
   // We will "forget" this lock once we backed out of all interpreter frames.
   // Here we are entering one, so ++.
   ++fInitialMutex.back().fRecurseCount;
}

////////////////////////////////////////////////////////////////////////////////

void TCling::ForgetMutexState()
{
   if (!fInitialMutex.back())
      return;
   if (fInitialMutex.back().fRecurseCount == 0) {
      Error("ForgetMutexState", "mutex state's recurse count already 0!");
   }
   else if (--fInitialMutex.back().fRecurseCount == 0) {
      // We have returned from all interpreter frames. Reset the initial lock state.
      fInitialMutex.back().fState.reset();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Re-apply the lock count delta that TCling__ResetInterpreterMutex() caused.

void TCling::ApplyToInterpreterMutex(void *delta)
{
   R__ASSERT(!fInitialMutex.empty() && "Inconsistent state of fInitialMutex!");
   if (gInterpreterMutex) {
      if (delta) {
         auto typedDelta = static_cast<TVirtualRWMutex::StateDelta *>(delta);
         std::unique_ptr<TVirtualRWMutex::StateDelta> uniqueP{typedDelta};
         gCoreMutex->Apply(std::move(uniqueP));
      }
   }
   fInitialMutex.pop_back();
}

////////////////////////////////////////////////////////////////////////////////
/// Reset the interpreter lock to the state it had before interpreter-related
/// calls happened.

void *TCling::RewindInterpreterMutex()
{
   if (fInitialMutex.back()) {
      std::unique_ptr<TVirtualRWMutex::StateDelta> uniqueP = gCoreMutex->Rewind(*fInitialMutex.back().fState);
      // Need to start a new recurse count.
      fInitialMutex.emplace_back();
      return uniqueP.release();
   }
   // Need to start a new recurse count.
   fInitialMutex.emplace_back();
   return nullptr;
}
