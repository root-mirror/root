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

#include "ROOT/FoundationUtils.hxx"

#include "TClingBaseClassInfo.h"
#include "TClingCallFunc.h"
#include "TClingClassInfo.h"
#include "TClingDataMemberInfo.h"
#include "TClingMethodArgInfo.h"
#include "TClingMethodInfo.h"
#include "TClingRdictModuleFileExtension.h"
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
#include "strlcpy.h"
#include "snprintf.h"
#include "TClingUtils.h"
#include "TVirtualCollectionProxy.h"
#include "TVirtualStreamerInfo.h"
#include "TListOfDataMembers.h"
#include "TListOfEnums.h"
#include "TListOfEnumsWithLock.h"
#include "TListOfFunctions.h"
#include "TListOfFunctionTemplates.h"
#include "TMemFile.h"
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
#include "clang/AST/DeclVisitor.h"
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
#include "clang/Parse/Parser.h"
#include "clang/Sema/Lookup.h"
#include "clang/Sema/Sema.h"
#include "clang/Serialization/ASTReader.h"
#include "clang/Serialization/GlobalModuleIndex.h"

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
#include "llvm/Object/ELFObjectFile.h"
#include "llvm/Object/ObjectFile.h"
#include "llvm/Object/SymbolicFile.h"
#include "llvm/Support/FileSystem.h"

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
#include <unordered_set>
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

#ifdef R__LINUX
# ifndef _GNU_SOURCE
#  define _GNU_SOURCE
# endif
# include <link.h> // dl_iterate_phdr()
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
         TClingDataMemberInfo* tcDmInfo = new TClingDataMemberInfo(GetInterpreterImpl(), *EDI, tcCInfo);
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
                                                 new TClingDataMemberInfo(GetInterpreterImpl(),
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
void TCling__InvalidateGlobal(const clang::Decl *D) {
   ((TCling*)gCling)->InvalidateGlobal(D);
}

extern "C"
void TCling__TransactionRollback(const cling::Transaction &T) {
   ((TCling*)gCling)->TransactionRollback(T);
}

extern "C" void TCling__LibraryLoadedRTTI(const void* dyLibHandle,
                                          const char* canonicalName) {
   ((TCling*)gCling)->LibraryLoaded(dyLibHandle, canonicalName);
}

extern "C" void TCling__RegisterRdictForLoadPCM(const std::string &pcmFileNameFullPath, llvm::StringRef *pcmContent)
{
   ((TCling *)gCling)->RegisterRdictForLoadPCM(pcmFileNameFullPath, pcmContent);
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
   auto tcling = new TCling("C++", "cling C++ Interpreter", argv);
   cling::DynamicLibraryManager::ExposeHiddenSharedLibrarySymbols(interpLibHandle);
   return tcling;
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
            TString header = gSystem->BaseName(cl->GetDeclFileName());
            TString dir = gSystem->GetDirName(cl->GetDeclFileName());
            TString dirbase(gSystem->BaseName(dir));
            while (dirbase.Length() && dirbase != "."
                   && dirbase != "include" && dirbase != "inc"
                   && dirbase != "prec_stl") {
               gSystem->PrependPathName(dirbase, header);
               dir = gSystem->GetDirName(dir);
            }
            fileContent += TString("#include \"") + header + "\"\n";
         }
      }
      for (it = fwdDecls.begin(); it != fwdDecls.end(); ++it) {
         fileContent += "class " + *it + ";\n";
      }
      fileContent += "#ifdef __CLING__ \n";
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
         if (enumTable && enumTable->THashList::FindObject(enName)) { delete [] scopeName; return true; }
      }
      // It may still be in one of the loaded protoclasses
      else if (auto scope = static_cast<TProtoClass *>(gClassTable->GetProtoNorm(scopeName))) {
         auto listOfEnums = scope->GetListOfEnums();
         if (listOfEnums) { // it could be null: no enumerators in the protoclass
            auto enumTable = dynamic_cast<const THashList *>(listOfEnums);
            if (enumTable && enumTable->THashList::FindObject(enName)) { delete [] scopeName; return true; }
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
static bool LoadModule(const std::string &ModuleName, cling::Interpreter &interp)
{
   // When starting up ROOT, cling would load all modulemap files on the include
   // paths. However, in a ROOT session, it is very common to run aclic which
   // will invoke rootcling and possibly produce a modulemap and a module in
   // the current folder.
   //
   // Before failing, try loading the modulemap in the current folder and try
   // loading the requested module from it.
   std::string currentDir = gSystem->WorkingDirectory();
   assert(!currentDir.empty());
   gCling->RegisterPrebuiltModulePath(currentDir);
   if (gDebug > 2)
      ::Info("TCling::__LoadModule", "Preloading module %s. \n",
             ModuleName.c_str());

   return interp.loadModule(ModuleName, /*Complain=*/true);
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

/// Checks if there is an ASTFile on disk for the given module \c M.
static bool HasASTFileOnDisk(clang::Module *M, const clang::Preprocessor &PP, std::string *FullFileName = nullptr)
{
   const HeaderSearchOptions &HSOpts = PP.getHeaderSearchInfo().getHeaderSearchOpts();

   std::string ModuleFileName;
   if (!HSOpts.PrebuiltModulePaths.empty())
      // Load the module from *only* in the prebuilt module path.
      ModuleFileName = PP.getHeaderSearchInfo().getPrebuiltModuleFileName(M->Name);
   if (FullFileName)
      *FullFileName = ModuleFileName;

   return !ModuleFileName.empty();
}

static bool HaveFullGlobalModuleIndex = false;
static GlobalModuleIndex *loadGlobalModuleIndex(cling::Interpreter &interp)
{
   CompilerInstance &CI = *interp.getCI();
   Preprocessor &PP = CI.getPreprocessor();
   auto ModuleManager = CI.getModuleManager();
   assert(ModuleManager);
   // StringRef ModuleIndexPath = HSI.getModuleCachePath();
   // HeaderSearch& HSI = PP.getHeaderSearchInfo();
   // HSI.setModuleCachePath(TROOT::GetLibDir().Data());
   std::string ModuleIndexPath = TROOT::GetLibDir().Data();
   if (ModuleIndexPath.empty())
      return nullptr;
   // Get an existing global index. This loads it if not already loaded.
   ModuleManager->resetForReload();
   ModuleManager->loadGlobalIndex();
   GlobalModuleIndex *GlobalIndex = ModuleManager->getGlobalIndex();

   // For finding modules needing to be imported for fixit messages,
   // we need to make the global index cover all modules, so we do that here.
   if (!GlobalIndex && !HaveFullGlobalModuleIndex) {
      ModuleMap &MMap = PP.getHeaderSearchInfo().getModuleMap();
      bool RecreateIndex = false;
      for (ModuleMap::module_iterator I = MMap.module_begin(), E = MMap.module_end(); I != E; ++I) {
         Module *TheModule = I->second;
         // We want the index only of the prebuilt modules.
         if (!HasASTFileOnDisk(TheModule, PP))
            continue;
         LoadModule(TheModule->Name, interp);
         RecreateIndex = true;
      }
      if (RecreateIndex) {
         cling::Interpreter::PushTransactionRAII deserRAII(&interp);
         clang::GlobalModuleIndex::UserDefinedInterestingIDs IDs;

         struct DefinitionFinder : public RecursiveASTVisitor<DefinitionFinder> {
            DefinitionFinder(clang::GlobalModuleIndex::UserDefinedInterestingIDs& IDs,
                             clang::TranslationUnitDecl* TU) : DefinitionIDs(IDs) {
               TraverseDecl(TU);
            }
            bool VisitNamedDecl(NamedDecl *ND) {
               if (!ND->isFromASTFile())
                  return true;
               if (!ND->getIdentifier())
                  return true;

               if (ND->getAccess() == AS_protected || ND->getAccess() == AS_private)
                  return true;

               if (TagDecl *TD = llvm::dyn_cast<TagDecl>(ND)) {
                  if (TD->isCompleteDefinition())
                     Register(TD);
               } else if (NamespaceDecl *NSD = llvm::dyn_cast<NamespaceDecl>(ND)) {
                  Register(NSD, /*AddSingleEntry=*/ false);
               }
               else if (TypedefNameDecl *TND = dyn_cast<TypedefNameDecl>(ND))
                  Register(TND);
               // FIXME: Add the rest...
               return true; // continue decending
            }
         private:
            clang::GlobalModuleIndex::UserDefinedInterestingIDs &DefinitionIDs;
            void Register(const NamedDecl* ND, bool AddSingleEntry = true) {
               assert(ND->isFromASTFile());
               // FIXME: All decls should have an owning module once rootcling
               // updates its generated decls from within the LookupHelper & co.
               if (!ND->hasOwningModule()) {
#ifndef NDEBUG
                  SourceManager &SM = ND->getASTContext().getSourceManager();
                  SourceLocation Loc = ND->getLocation();
                  const FileEntry *FE = SM.getFileEntryForID(SM.getFileID(Loc));
                  (void)FE;
                  assert(FE->getName().contains("input_line_"));
#endif
                  return;
               }

               Module *OwningModule = ND->getOwningModule()->getTopLevelModule();
               assert(OwningModule);
               assert(!ND->getName().empty() && "Empty name");
               if (AddSingleEntry && DefinitionIDs.count(ND->getName()))
                  return;
               // FIXME: The FileEntry in not stable to serialize.
               // FIXME: We might end up with many times with the same module.
               // FIXME: We might end up two modules containing a definition.
               // FIXME: What do we do if no definition is found.
               DefinitionIDs[ND->getName()].push_back(OwningModule->getASTFile());
            }
         };
         DefinitionFinder defFinder(IDs, CI.getASTContext().getTranslationUnitDecl());

         llvm::cantFail(GlobalModuleIndex::writeIndex(CI.getFileManager(),
                                                      CI.getPCHContainerReader(),
                                                      ModuleIndexPath,
                                                      &IDs));
         ModuleManager->resetForReload();
         ModuleManager->loadGlobalIndex();
         GlobalIndex = ModuleManager->getGlobalIndex();
      }
      HaveFullGlobalModuleIndex = true;
   }
   return GlobalIndex;
}

static void RegisterCxxModules(cling::Interpreter &clingInterp)
{
   if (!clingInterp.getCI()->getLangOpts().Modules)
      return;

   // Loading of a module might deserialize.
   cling::Interpreter::PushTransactionRAII deserRAII(&clingInterp);

   // Setup core C++ modules if we have any to setup.

   // Load libc and stl first.
   // Load vcruntime module for windows
#ifdef R__WIN32
   LoadModule("vcruntime", clingInterp);
   LoadModule("services", clingInterp);
#endif

#ifdef R__MACOSX
   LoadModule("Darwin", clingInterp);
#else
   LoadModule("libc", clingInterp);
#endif
   LoadModule("std", clingInterp);

   LoadModule("_Builtin_intrinsics", clingInterp);

   // Load core modules
   // This should be vector in order to be able to pass it to LoadModules
   std::vector<std::string> CoreModules = {"ROOT_Foundation_C",
                                           "ROOT_Config",
                                           "ROOT_Rtypes",
                                           "ROOT_Foundation_Stage1_NoRTTI",
                                           "Core",
                                           "Rint",
                                           "RIO"};

   LoadModules(CoreModules, clingInterp);

   // Take this branch only from ROOT because we don't need to preload modules in rootcling
   if (!IsFromRootCling()) {
      std::vector<std::string> CommonModules = {"MathCore"};
      LoadModules(CommonModules, clingInterp);

      // These modules should not be preloaded but they fix issues.
      // FIXME: Hist is not a core module but is very entangled to MathCore and
      // causes issues.
      std::vector<std::string> FIXMEModules = {"Hist"};
      LoadModules(FIXMEModules, clingInterp);

      clang::CompilerInstance &CI = *clingInterp.getCI();
      GlobalModuleIndex *GlobalIndex = nullptr;
      // Conservatively enable platform by platform.
      bool supportedPlatform =
#ifdef R__LINUX
         true
#elif defined(R__MACOSX)
         true
#else // Windows
         false
#endif
         ;
      // Allow forcefully enabling/disabling the GMI.
      llvm::Optional<std::string> envUseGMI = llvm::sys::Process::GetEnv("ROOT_USE_GMI");
      if (envUseGMI.hasValue()) {
         if (!envUseGMI->empty() && !ROOT::FoundationUtils::CanConvertEnvValueToBool(*envUseGMI))
            ::Warning("TCling__RegisterCxxModules",
                      "Cannot convert '%s' to bool, setting to false!",
                      envUseGMI->c_str());

         bool value = envUseGMI->empty() || ROOT::FoundationUtils::ConvertEnvValueToBool(*envUseGMI);

         if (supportedPlatform == value)
            ::Warning("TCling__RegisterCxxModules", "Global module index is%sused already!",
                     (value) ? " " :" not ");
         supportedPlatform = value;
      }

      if (supportedPlatform) {
         loadGlobalModuleIndex(clingInterp);
         // FIXME: The ASTReader still calls loadGlobalIndex and loads the file
         // We should investigate how to suppress it completely.
         GlobalIndex = CI.getModuleManager()->getGlobalIndex();
      }

      llvm::StringSet<> KnownModuleFileNames;
      if (GlobalIndex)
         GlobalIndex->getKnownModuleFileNames(KnownModuleFileNames);

      clang::Preprocessor &PP = CI.getPreprocessor();
      std::vector<std::string> PendingModules;
      PendingModules.reserve(256);
      ModuleMap &MMap = PP.getHeaderSearchInfo().getModuleMap();
      for (auto I = MMap.module_begin(), E = MMap.module_end(); I != E; ++I) {
         clang::Module *M = I->second;
         assert(M);

         // We want to load only already created modules.
         std::string FullASTFilePath;
         if (!HasASTFileOnDisk(M, PP, &FullASTFilePath))
            continue;

         if (GlobalIndex && KnownModuleFileNames.count(FullASTFilePath))
            continue;

         if (M->IsMissingRequirement)
            continue;

         if (GlobalIndex)
            LoadModule(M->Name, clingInterp);
         else {
            // FIXME: We may be able to remove those checks as cling::loadModule
            // checks if a module was alredy loaded.
            if (std::find(CoreModules.begin(), CoreModules.end(), M->Name) != CoreModules.end())
               continue; // This is a core module which was already loaded.

            // Load system modules now and delay the other modules after we have
            // loaded all system ones.
            if (M->IsSystem)
               LoadModule(M->Name, clingInterp);
            else
               PendingModules.push_back(M->Name);
         }
      }
      LoadModules(PendingModules, clingInterp);
   }

   // Check that the gROOT macro was exported by any core module.
   assert(clingInterp.getMacro("gROOT") && "Couldn't load gROOT macro?");

   // `ERROR` and `PI` are from loading R related modules, which conflict with
   // user's code.
   clingInterp.declare(R"CODE(
#ifdef PI
# undef PI
#endif
#ifdef ERROR
# undef ERROR
#endif
                       )CODE");
}

static void RegisterPreIncludedHeaders(cling::Interpreter &clingInterp)
{
   std::string PreIncludes;
   bool hasCxxModules = clingInterp.getCI()->getLangOpts().Modules;

   // For the list to also include string, we have to include it now.
   // rootcling does parts already if needed, e.g. genreflex does not want using
   // namespace std.
   if (IsFromRootCling()) {
      PreIncludes += "#include \"RtypesCore.h\"\n";
   } else {
      if (!hasCxxModules)
         PreIncludes += "#include \"Rtypes.h\"\n";

      PreIncludes += gClassDefInterpMacro + "\n"
                     + gInterpreterClassDef + "\n"
                     "#undef ClassImp\n"
                     "#define ClassImp(X);\n";
   }
   if (!hasCxxModules)
      PreIncludes += "#include <string>\n";

   // We must include it even when we have modules because it is marked as
   // textual in the modulemap due to the nature of the assert header.
#ifndef R__WIN32
   PreIncludes += "#include <cassert>\n";
#endif
   PreIncludes += "using namespace std;\n";
   clingInterp.declare(PreIncludes);
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize the cling interpreter interface.
/// \param argv - array of arguments passed to the cling::Interpreter constructor
///               e.g. `-DFOO=bar`. The last element of the array must be `nullptr`.

TCling::TCling(const char *name, const char *title, const char* const argv[])
: TInterpreter(name, title), fMore(0), fGlobalsListSerial(-1), fMapfile(nullptr),
  fRootmapFiles(nullptr), fLockProcessLine(true), fNormalizedCtxt(0),
  fPrevLoadedDynLibInfo(0), fClingCallbacks(0), fAutoLoadCallBack(0),
  fTransactionCount(0), fHeaderParsingOnDemand(true), fIsAutoParsingSuspended(kFALSE)
{
   fPrompt[0] = 0;
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

   auto GetEnvVarPath = [](const std::string &EnvVar,
                       std::vector<std::string> &Paths) {
      llvm::Optional<std::string> EnvOpt = llvm::sys::Process::GetEnv(EnvVar);
      if (EnvOpt.hasValue()) {
         StringRef Env(*EnvOpt);
         while (!Env.empty()) {
            StringRef Arg;
            std::tie(Arg, Env) = Env.split(ROOT::FoundationUtils::GetEnvPathSeparator());
            if (std::find(Paths.begin(), Paths.end(), Arg.str()) == Paths.end())
               Paths.push_back(Arg.str());
         }
      }
   };

   if (fCxxModulesEnabled) {
      std::vector<std::string> Paths;
      // ROOT usually knows better where its libraries are. This way we can
      // discover modules without having to should thisroot.sh and should fix
      // gnuinstall.
#ifdef R__WIN32
      Paths.push_back(TROOT::GetBinDir().Data());
#else
      Paths.push_back(TROOT::GetLibDir().Data());
#endif
      GetEnvVarPath("CLING_PREBUILT_MODULE_PATH", Paths);
      std::string EnvVarPath;
      for (const std::string& P : Paths)
         EnvVarPath += P + ROOT::FoundationUtils::GetEnvPathSeparator();
      // FIXME: We should make cling -fprebuilt-module-path work.
      gSystem->Setenv("CLING_PREBUILT_MODULE_PATH", EnvVarPath.c_str());
   }

   // FIXME: This only will enable frontend timing reports.
   EnvOpt = llvm::sys::Process::GetEnv("ROOT_CLING_TIMING");
   if (EnvOpt.hasValue())
     clingArgsStorage.push_back("-ftime-report");

   // Add the overlay file. Note that we cannot factor it out for both root
   // and rootcling because rootcling activates modules only if -cxxmodule
   // flag is passed.
   if (fCxxModulesEnabled && !fromRootCling) {
      // For now we prefer rootcling to enumerate explicitly its modulemaps.
      std::vector<std::string> ModuleMaps;
      std::string ModuleMapSuffix = ROOT::FoundationUtils::GetPathSeparator() + "module.modulemap";
      ModuleMaps.push_back(TROOT::GetIncludeDir().Data() + ModuleMapSuffix);
      GetEnvVarPath("CLING_MODULEMAP_FILES", ModuleMaps);

      std::string cwd = gSystem->WorkingDirectory();
      // Give highest precedence of the modulemap in the cwd if any.
      if (llvm::sys::fs::exists(cwd + ModuleMapSuffix))
         ModuleMaps.push_back(cwd + ModuleMapSuffix);

      for (const std::string& M : ModuleMaps)
         clingArgsStorage.push_back("-fmodule-map-file=" + M);

      std::string ModulesCachePath;
      EnvOpt = llvm::sys::Process::GetEnv("CLING_MODULES_CACHE_PATH");
      if (EnvOpt.hasValue()){
         StringRef Env(*EnvOpt);
         assert(llvm::sys::fs::exists(Env) && "Path does not exist!");
         ModulesCachePath = Env.str();
      } else {
         ModulesCachePath = TROOT::GetLibDir();
      }

      clingArgsStorage.push_back("-fmodules-cache-path=" + ModulesCachePath);
   }

   std::vector<const char*> interpArgs;
   for (std::vector<std::string>::const_iterator iArg = clingArgsStorage.begin(),
           eArg = clingArgsStorage.end(); iArg != eArg; ++iArg)
      interpArgs.push_back(iArg->c_str());

   // Activate C++ modules support. If we are running within rootcling, it's up
   // to rootcling to set this flag depending on whether it wants to produce
   // C++ modules.
   TString vfsArg;
   if (fCxxModulesEnabled) {
      if (!fromRootCling) {
         // We only set this flag, rest is done by the CIFactory.
         interpArgs.push_back("-fmodules");
         interpArgs.push_back("-fno-implicit-module-maps");
         // We should never build modules during runtime, so let's enable the
         // module build remarks from clang to make it easier to spot when we do
         // this by accident.
         interpArgs.push_back("-Rmodule-build");
      }
      // ROOT implements its AutoLoading upon module's link directives. We
      // generate module A { header "A.h" link "A.so" export * } where ROOT's
      // facilities use the link directive to dynamically load the relevant
      // library. So, we need to suppress clang's default autolink behavior.
      interpArgs.push_back("-fno-autolink");
   }

#ifdef R__FAST_MATH
   // Same setting as in rootcling_impl.cxx.
   interpArgs.push_back("-ffast-math");
#endif

   TString llvmResourceDir = TROOT::GetEtcDir() + "/cling";
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

   // Add the Rdict module file extension.
   cling::Interpreter::ModuleFileExtensions extensions;
   EnvOpt = llvm::sys::Process::GetEnv("ROOTDEBUG_RDICT");
   if (!EnvOpt.hasValue())
      extensions.push_back(std::make_shared<TClingRdictModuleFileExtension>());

   fInterpreter = llvm::make_unique<cling::Interpreter>(interpArgs.size(),
                                                        &(interpArgs[0]),
                                                        llvmResourceDir, extensions);

   // Don't check whether modules' files exist.
   fInterpreter->getCI()->getPreprocessorOpts().DisablePCHValidation = true;

   // Until we can disable AutoLoading during Sema::CorrectTypo() we have
   // to disable spell checking.
   fInterpreter->getCI()->getLangOpts().SpellChecking = false;

   // We need stream that doesn't close its file descriptor, thus we are not
   // using llvm::outs. Keeping file descriptor open we will be able to use
   // the results in pipes (Savannah #99234).
   static llvm::raw_fd_ostream fMPOuts (STDOUT_FILENO, /*ShouldClose*/false);
   fMetaProcessor = llvm::make_unique<cling::MetaProcessor>(*fInterpreter, fMPOuts);

   RegisterCxxModules(*fInterpreter);
   RegisterPreIncludedHeaders(*fInterpreter);

   // We are now ready (enough is loaded) to init the list of opaque typedefs.
   fNormalizedCtxt = new ROOT::TMetaUtils::TNormalizedCtxt(fInterpreter->getLookupHelper());
   fLookupHelper = new ROOT::TMetaUtils::TClingLookupHelper(*fInterpreter, *fNormalizedCtxt,
                                                            TClingLookupHelper__ExistingTypeCheck,
                                                            TClingLookupHelper__AutoParse,
                                                            &fIsShuttingDown);
   TClassEdit::Init(fLookupHelper);

   // Disallow auto-parsing in rootcling
   fIsAutoParsingSuspended = fromRootCling;

   ResetAll();

   // Enable dynamic lookup
   if (!fromRootCling) {
      fInterpreter->enableDynamicLookup();
   }

   // Enable ClinG's DefinitionShadower for ROOT.
   fInterpreter->getRuntimeOptions().AllowRedefinition = 1;

   // Attach cling callbacks last; they might need TROOT::fInterpreter
   // and should thus not be triggered during the equivalent of
   // TROOT::fInterpreter = new TCling;
   std::unique_ptr<TClingCallbacks>
      clingCallbacks(new TClingCallbacks(GetInterpreterImpl(), /*hasCodeGen*/ !fromRootCling));
   fClingCallbacks = clingCallbacks.get();
   fClingCallbacks->SetAutoParsingSuspended(fIsAutoParsingSuspended);
   fInterpreter->setCallbacks(std::move(clingCallbacks));

   if (!fromRootCling) {
      cling::DynamicLibraryManager& DLM = *fInterpreter->getDynamicLibraryManager();
      // Make sure cling looks into ROOT's libdir, even if not part of LD_LIBRARY_PATH
      // e.g. because of an RPATH build.
      DLM.addSearchPath(TROOT::GetLibDir().Data(), /*isUser=*/true,
                        /*prepend=*/true);
      auto ShouldPermanentlyIgnore = [](llvm::StringRef FileName) -> bool{
         llvm::StringRef stem = llvm::sys::path::stem(FileName);
         return stem.startswith("libNew") || stem.startswith("libcppyy_backend");
      };
      // Initialize the dyld for the llvmLazyFunctionCreator.
      DLM.initializeDyld(ShouldPermanentlyIgnore);
      fInterpreter->installLazyFunctionCreator(llvmLazyFunctionCreator);
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Destroy the interpreter interface.

TCling::~TCling()
{
   // ROOT's atexit functions require the interepreter to be available.
   // Run them before shutting down.
   if (!IsFromRootCling())
      GetInterpreterImpl()->runAtExitFuncs();
   fIsShuttingDown = true;
   delete fMapfile;
   delete fRootmapFiles;
   delete fTemporaries;
   delete fNormalizedCtxt;
   delete fLookupHelper;
   gCling = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize the interpreter, once TROOT::fInterpreter is set.

void TCling::Initialize()
{
   fClingCallbacks->Initialize();

   // We are set up. Enable ROOT's AutoLoading.
   if (IsFromRootCling())
      return;

   // Read the rules before enabling the auto loading to not inadvertently
   // load the libraries for the classes concerned even-though the user is
   // *not* using them.
   // Note this call must happen before the first call to LoadLibraryMap.
   assert(GetRootMapFiles() == 0 && "Must be called before LoadLibraryMap!");
   TClass::ReadRules(); // Read the default customization rules ...

   LoadLibraryMap();
   SetClassAutoLoading(true);
}

void TCling::ShutDown()
{
   fIsShuttingDown = true;
   ResetGlobals();
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
/// Register Rdict data for future loading by LoadPCM;

void TCling::RegisterRdictForLoadPCM(const std::string &pcmFileNameFullPath, llvm::StringRef *pcmContent)
{
   if (IsFromRootCling())
      return;

   if (llvm::sys::fs::exists(pcmFileNameFullPath)) {
      ::Error("TCling::RegisterRdictForLoadPCM", "Rdict '%s' is both in Module extension and in File system.", pcmFileNameFullPath.c_str());
      return;
   }

   // The pcmFileNameFullPath must be resolved already because we cannot resolve
   // a link to a non-existent file.
   fPendingRdicts[pcmFileNameFullPath] = *pcmContent;
}

////////////////////////////////////////////////////////////////////////////////
/// Tries to load a PCM from TFile; returns true on success.

void TCling::LoadPCMImpl(TFile &pcmFile)
{
   auto listOfKeys = pcmFile.GetListOfKeys();

   // This is an empty pcm
   if (listOfKeys && ((listOfKeys->GetSize() == 0) ||                           // Nothing here, or
                      ((listOfKeys->GetSize() == 1) &&                          // only one, and
                       !strcmp(((TKey *)listOfKeys->At(0))->GetName(), "EMPTY") // name is EMPTY
                       ))) {
      return;
   }

   TObjArray *protoClasses;
   if (gDebug > 1)
      ::Info("TCling::LoadPCMImpl", "reading protoclasses for %s \n", pcmFile.GetName());

   pcmFile.GetObject("__ProtoClasses", protoClasses);

   if (protoClasses) {
      for (auto obj : *protoClasses) {
         TProtoClass *proto = (TProtoClass *)obj;
         TClassTable::Add(proto);
      }
      // Now that all TClass-es know how to set them up we can update
      // existing TClasses, which might cause the creation of e.g. TBaseClass
      // objects which in turn requires the creation of TClasses, that could
      // come from the PCH, but maybe later in the loop. Instead of resolving
      // a dependency graph the addition to the TClassTable above allows us
      // to create these dependent TClasses as needed below.
      for (auto proto : *protoClasses) {
         if (TClass *existingCl = (TClass *)gROOT->GetListOfClasses()->FindObject(proto->GetName())) {
            // We have an existing TClass object. It might be emulated
            // or interpreted; we now have more information available.
            // Make that available.
            if (existingCl->GetState() != TClass::kHasTClassInit) {
               DictFuncPtr_t dict = gClassTable->GetDict(proto->GetName());
               if (!dict) {
                  ::Error("TCling::LoadPCM", "Inconsistent TClassTable for %s", proto->GetName());
               } else {
                  // This will replace the existing TClass.
                  TClass *ncl = (*dict)();
                  if (ncl)
                     ncl->PostLoadCheck();
               }
            }
         }
      }

      protoClasses->Clear(); // Ownership was transfered to TClassTable.
      delete protoClasses;
   }

   TObjArray *dataTypes;
   pcmFile.GetObject("__Typedefs", dataTypes);
   if (dataTypes) {
      for (auto typedf : *dataTypes)
         gROOT->GetListOfTypes()->Add(typedf);
      dataTypes->Clear(); // Ownership was transfered to TListOfTypes.
      delete dataTypes;
   }

   TObjArray *enums;
   pcmFile.GetObject("__Enums", enums);
   if (enums) {
      // Cache the pointers
      auto listOfGlobals = gROOT->GetListOfGlobals();
      auto listOfEnums = dynamic_cast<THashList *>(gROOT->GetListOfEnums());
      // Loop on enums and then on enum constants
      for (auto selEnum : *enums) {
         const char *enumScope = selEnum->GetTitle();
         const char *enumName = selEnum->GetName();
         if (strcmp(enumScope, "") == 0) {
            // This is a global enum and is added to the
            // list of enums and its constants to the list of globals
            if (!listOfEnums->THashList::FindObject(enumName)) {
               ((TEnum *)selEnum)->SetClass(nullptr);
               listOfEnums->Add(selEnum);
            }
            for (auto enumConstant : *static_cast<TEnum *>(selEnum)->GetConstants()) {
               if (!listOfGlobals->FindObject(enumConstant)) {
                  listOfGlobals->Add(enumConstant);
               }
            }
         } else {
            // This enum is in a namespace. A TClass entry is bootstrapped if
            // none exists yet and the enum is added to it
            TClass *nsTClassEntry = TClass::GetClass(enumScope);
            if (!nsTClassEntry) {
               nsTClassEntry = new TClass(enumScope, 0, TClass::kNamespaceForMeta, true);
            }
            auto listOfEnums = nsTClassEntry->fEnums.load();
            if (!listOfEnums) {
               if ((kIsClass | kIsStruct | kIsUnion) & nsTClassEntry->Property()) {
                  // For this case, the list will be immutable once constructed
                  // (i.e. in this case, by the end of this routine).
                  listOfEnums = nsTClassEntry->fEnums = new TListOfEnums(nsTClassEntry);
               } else {
                  // namespaces can have enums added to them
                  listOfEnums = nsTClassEntry->fEnums = new TListOfEnumsWithLock(nsTClassEntry);
               }
            }
            if (listOfEnums && !listOfEnums->THashList::FindObject(enumName)) {
               ((TEnum *)selEnum)->SetClass(nsTClassEntry);
               listOfEnums->Add(selEnum);
            }
         }
      }
      enums->Clear();
      delete enums;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Tries to load a rdict PCM, issues diagnostics if it fails.

void TCling::LoadPCM(std::string pcmFileNameFullPath)
{
   SuspendAutoLoadingRAII autoloadOff(this);
   SuspendAutoParsing autoparseOff(this);
   assert(!pcmFileNameFullPath.empty());
   assert(llvm::sys::path::is_absolute(pcmFileNameFullPath));

   // Easier to work with the ROOT interfaces.
   TString pcmFileName = pcmFileNameFullPath;

   // Prevent the ROOT-PCMs hitting this during auto-load during
   // JITting - which will cause recursive compilation.
   // Avoid to call the plugin manager at all.
   R__InitStreamerInfoFactory();

   TDirectory::TContext ctxt;
   llvm::SaveAndRestore<Int_t> SaveGDebug(gDebug);
   if (gDebug > 5) {
      gDebug -= 5;
      ::Info("TCling::LoadPCM", "Loading ROOT PCM %s", pcmFileName.Data());
   } else {
      gDebug = 0;
   }

   if (llvm::sys::fs::is_symlink_file(pcmFileNameFullPath))
      pcmFileNameFullPath = ROOT::TMetaUtils::GetRealPath(pcmFileNameFullPath);

   auto pendingRdict = fPendingRdicts.find(pcmFileNameFullPath);
   if (pendingRdict != fPendingRdicts.end()) {
      llvm::StringRef pcmContent = pendingRdict->second;
      TMemFile::ZeroCopyView_t range{pcmContent.data(), pcmContent.size()};
      std::string RDictFileOpts = pcmFileNameFullPath + "?filetype=pcm";
      TMemFile pcmMemFile(RDictFileOpts.c_str(), range);

      cling::Interpreter::PushTransactionRAII deserRAII(GetInterpreterImpl());
      LoadPCMImpl(pcmMemFile);
      fPendingRdicts.erase(pendingRdict);

      return;
   }

   if (!llvm::sys::fs::exists(pcmFileNameFullPath)) {
      ::Error("TCling::LoadPCM", "ROOT PCM %s file does not exist",
              pcmFileNameFullPath.data());
      if (!fPendingRdicts.empty())
         for (const auto &rdict : fPendingRdicts)
            ::Info("TCling::LoadPCM", "In-memory ROOT PCM candidate %s\n",
                   rdict.first.c_str());
      return;
   }

   if (!gROOT->IsRootFile(pcmFileName)) {
      Fatal("LoadPCM", "The file %s is not a ROOT as was expected\n", pcmFileName.Data());
      return;
   }
   TFile pcmFile(pcmFileName + "?filetype=pcm", "READ");
   LoadPCMImpl(pcmFile);
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
///\returns true if the module map was loaded, false on error or if the map was
///         already loaded.
bool TCling::RegisterPrebuiltModulePath(const std::string &FullPath,
                                        const std::string &ModuleMapName /*= "module.modulemap"*/) const
{
   assert(llvm::sys::path::is_absolute(FullPath));
   Preprocessor &PP = fInterpreter->getCI()->getPreprocessor();
   FileManager &FM = PP.getFileManager();
   // FIXME: In a ROOT session we can add an include path (through .I /inc/path)
   // We should look for modulemap files there too.
   const DirectoryEntry *DE = FM.getDirectory(FullPath);
   if (DE) {
      HeaderSearch &HS = PP.getHeaderSearchInfo();
      HeaderSearchOptions &HSOpts = HS.getHeaderSearchOpts();
      const auto &ModPaths = HSOpts.PrebuiltModulePaths;
      bool pathExists = std::find(ModPaths.begin(), ModPaths.end(), FullPath) != ModPaths.end();
      if (!pathExists)
         HSOpts.AddPrebuiltModulePath(FullPath);
      // We cannot use HS.lookupModuleMapFile(DE, /*IsFramework*/ false);
      // because its internal call to getFile has CacheFailure set to true.
      // In our case, modulemaps can appear any time due to ACLiC.
      // Code copied from HS.lookupModuleMapFile.
      llvm::SmallString<256> ModuleMapFileName(DE->getName());
      llvm::sys::path::append(ModuleMapFileName, ModuleMapName);
      const FileEntry *FE = FM.getFile(ModuleMapFileName, /*openFile*/ false,
                                       /*CacheFailure*/ false);
      if (FE) {
         if (!HS.loadModuleMapFile(FE, /*IsSystem*/ false))
            return true;
         Error("RegisterPrebuiltModulePath", "Could not load modulemap in %s", ModuleMapFileName.c_str());
      }
   }
   return false;
}

////////////////////////////////////////////////////////////////////////////////
/// List of dicts that have the PCM information already in the PCH.
static const std::unordered_set<std::string> gIgnoredPCMNames = {"libCore",
                                                                 "libRint",
                                                                 "libThread",
                                                                 "libRIO",
                                                                 "libImt",
                                                                 "libMultiProc",
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

static void PrintDlError(const char *dyLibName, const char *modulename)
{
#ifdef R__WIN32
   char dyLibError[1000];
   FormatMessageA(FORMAT_MESSAGE_FROM_SYSTEM, NULL, GetLastError(), MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
                  dyLibError, sizeof(dyLibError), NULL);
#else
   const char *dyLibError = dlerror();
#endif
   ::Error("TCling::RegisterModule", "Cannot open shared library %s for dictionary %s:\n  %s", dyLibName, modulename,
           (dyLibError) ? dyLibError : "");
}

////////////////////////////////////////////////////////////////////////////////
// Update all the TClass registered in fClassesToUpdate

void TCling::ProcessClassesToUpdate()
{
   while (!fClassesToUpdate.empty()) {
      TClass *oldcl = fClassesToUpdate.back().first;
      // If somehow the TClass has already been loaded (maybe it was registered several time),
      // we skip it.  Otherwise, the existing TClass is in mode kInterpreted, kEmulated or
      // maybe even kForwardDeclared and needs to replaced.
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

   // Make sure we do not set off AutoLoading or autoparsing during the
   // module registration!
   SuspendAutoLoadingRAII autoLoadOff(this);

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

   std::string dyLibName = cling::DynamicLibraryManager::getSymbolLocation(triggerFunc);
   assert(!llvm::sys::fs::is_symlink_file(dyLibName));

   if (dyLibName.empty()) {
      ::Error("TCling::RegisterModule", "Dictionary trigger function for %s not found", modulename);
      return;
   }

   // The triggerFunc may not be in a shared object but in an executable.
   bool isSharedLib = cling::DynamicLibraryManager::isSharedLibrary(dyLibName);

   bool wasDlopened = false;

   // If this call happens after dlopen has finished (i.e. late registration)
   // there is no need to dlopen the library recursively. See ROOT-8437 where
   // the dyLibName would correspond to the binary.
   if (!lateRegistration) {

      if (isSharedLib) {
         // We need to open the dictionary shared library, to resolve symbols
         // requested by the JIT from it: as the library is currently being dlopen'ed,
         // its symbols are not yet reachable from the process.
         // Recursive dlopen seems to work just fine.
         void* dyLibHandle = dlopen(dyLibName.c_str(), RTLD_LAZY | RTLD_GLOBAL);
         if (dyLibHandle) {
            fRegisterModuleDyLibs.push_back(dyLibHandle);
            wasDlopened = true;
         } else {
            PrintDlError(dyLibName.c_str(), modulename);
         }
      }
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

      if (!fwdDeclsCodeLessEnums.empty()){ // Avoid the overhead if nothing is to be declared
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

   clang::Sema &TheSema = fInterpreter->getSema();

   bool ModuleWasSuccessfullyLoaded = false;
   if (hasCxxModule) {
      std::string ModuleName = modulename;
      if (llvm::StringRef(modulename).startswith("lib"))
         ModuleName = llvm::StringRef(modulename).substr(3).str();

      // In case we are directly loading the library via gSystem->Load() without
      // specifying the relevant include paths we should try loading the
      // modulemap next to the library location.
      clang::Preprocessor &PP = TheSema.getPreprocessor();
      std::string ModuleMapName;
      if (isACLiC)
         ModuleMapName = ModuleName + ".modulemap";
      else
         ModuleMapName = "module.modulemap";
      RegisterPrebuiltModulePath(llvm::sys::path::parent_path(dyLibName),
                                 ModuleMapName);

      // FIXME: We should only complain for modules which we know to exist. For example, we should not complain about
      // modules such as GenVector32 because it needs to fall back to GenVector.
      cling::Interpreter::PushTransactionRAII deserRAII(GetInterpreterImpl());
      ModuleWasSuccessfullyLoaded = LoadModule(ModuleName, *fInterpreter);
      if (!ModuleWasSuccessfullyLoaded) {
         // Only report if we found the module in the modulemap.
         clang::HeaderSearch &headerSearch = PP.getHeaderSearchInfo();
         clang::ModuleMap &moduleMap = headerSearch.getModuleMap();
         if (moduleMap.findModule(ModuleName))
            Info("TCling::RegisterModule", "Module %s in modulemap failed to load.", ModuleName.c_str());
      }
   }

   if (gIgnoredPCMNames.find(modulename) == gIgnoredPCMNames.end()) {
      llvm::SmallString<256> pcmFileNameFullPath(dyLibName);
      // The path dyLibName might not be absolute. This can happen if dyLibName
      // is linked to an executable in the same folder.
      llvm::sys::fs::make_absolute(pcmFileNameFullPath);
      llvm::sys::path::remove_filename(pcmFileNameFullPath);
      llvm::sys::path::append(pcmFileNameFullPath,
                              ROOT::TMetaUtils::GetModuleFileName(modulename));
      LoadPCM(pcmFileNameFullPath.str().str());
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
   // loaded from other modules or from class rules or from opening a TFile
   // or from loading header in a way that did not provoke the loading of
   // the library we just loaded).
   ProcessClassesToUpdate();

   if (!ModuleWasSuccessfullyLoaded && !hasHeaderParsingOnDemand) {
      // __ROOTCLING__ might be pulled in through PCH
      fInterpreter->declare("#ifdef __ROOTCLING__\n"
                            "#undef __ROOTCLING__\n"
                            + gInterpreterClassDef +
                            "#endif");
   }

   if (wasDlopened) {
      assert(isSharedLib);
      void* dyLibHandle = fRegisterModuleDyLibs.back();
      fRegisterModuleDyLibs.pop_back();
      dlclose(dyLibHandle);
   }
}

void TCling::AddAvailableIndentifiers(TSeqCollection& Idents) {
   clang::CompilerInstance& CI = *GetInterpreterImpl()->getCI();
   ASTContext &C = CI.getASTContext();

   // Do not do anything if we have no global module index.
   // FIXME: This is mostly to real with false positives in the TTabCom
   // interface for non-modules.
   if (!fCxxModulesEnabled)
      return;

   if (IdentifierInfoLookup *External = C.Idents.getExternalIdentifierLookup()) {
      std::unique_ptr<IdentifierIterator> Iter(External->getIdentifiers());
      for (llvm::StringRef Ident = Iter->Next(); !Ident.empty(); Ident = Iter->Next()) {
         std::string I = Ident.str();
         if (!Idents.Contains(I.data()))
            Idents.Add(new TObjString(I.c_str()));
      }
   }
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

Longptr_t TCling::ProcessLine(const char* line, EErrorCode* error/*=0*/)
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
   } interpreterFlagsRAII(GetInterpreterImpl());

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
               indent = HandleInterpreterException(GetMetaProcessorImpl(), mod_line, compRes, &result);
            }
         }
      } else {
         // not ACLiC
         size_t unnamedMacroOpenCurly;
         {
            std::string code;
            std::string codeline;
            // Windows requires std::ifstream::binary to properly handle
            // CRLF and LF line endings
            std::ifstream in(fname, std::ifstream::binary);
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
            indent = HandleInterpreterException(GetMetaProcessorImpl(), mod_line, compRes, &result);
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
            indent = HandleInterpreterException(GetMetaProcessorImpl(), sLine, compRes, &result);
         } else {
            indent = HandleInterpreterException(GetMetaProcessorImpl(), sLine, compRes, &result);
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
      return result.simplisticCastAs<Longptr_t>();
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// No-op; see TRint instead.

void TCling::PrintIntro()
{
}

////////////////////////////////////////////////////////////////////////////////
/// \brief Add a directory to the list of directories in which the
///        interpreter looks for include files.
/// \param[in] path The path to the directory.
/// \note Only one path item can be specified at a time, i.e. "path1:path2" is
///       \b NOT supported.
/// \warning Only the path to the directory should be specified, without
///          prepending the \c -I prefix, i.e.
///          <tt>gCling->AddIncludePath("/path/to/my/includes")</tt>. If the
///          \c -I prefix is used it will be ignored.
void TCling::AddIncludePath(const char *path)
{
   R__LOCKGUARD(gInterpreterMutex);
   // Favorite source of annoyance: gSystem->AddIncludePath() needs "-I",
   // gCling->AddIncludePath() does not! Work around that inconsistency:
   if (path[0] == '-' && path[1] == 'I')
      path += 2;
   TString sPath(path);
   gSystem->ExpandPathName(sPath);
   fInterpreter->AddIncludePath(sPath.Data());
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
      cling::Interpreter::PushTransactionRAII deserRAII(GetInterpreterImpl());

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
         if (fieldRecDecl && !fieldRecDecl->isAnonymousStructOrUnion()) {
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

   SuspendAutoLoadingRAII autoLoadOff(this);
   SuspendAutoParsing autoParseRaii(this);

   bool oldDynLookup = fInterpreter->isDynamicLookupEnabled();
   fInterpreter->enableDynamicLookup(false);
   bool oldRawInput = fInterpreter->isRawInputEnabled();
   fInterpreter->enableRawInput(true);

   Bool_t ret = LoadText(code);

   fInterpreter->enableRawInput(oldRawInput);
   fInterpreter->enableDynamicLookup(oldDynLookup);
   return ret;
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
   return s_IsLibraryLoaded(libname, GetInterpreterImpl());
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if ROOT has cxxmodules pcm for a given library name.
// FIXME: We need to be able to support lazy loading of pcm generated by ACLiC.
Bool_t TCling::HasPCMForLibrary(const char *libname) const
{
   llvm::StringRef ModuleName(libname);
   ModuleName = llvm::sys::path::stem(ModuleName);
   ModuleName.consume_front("lib");

   // FIXME: In case when the modulemap is not yet loaded we will return the
   // wrong result. Consider a call to HasPCMForLibrary(../test/libEvent.so)
   // We will only load the modulemap for libEvent.so after we dlopen libEvent
   // which may happen after calling this interface. Maybe we should also check
   // if there is a Event.pcm file and a module.modulemap, load it and return
   // true.
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
   if (s_IsLibraryLoaded(file_name.c_str(), GetInterpreterImpl()))
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
                                              /*IsFrameworkFound*/ nullptr,
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


#if defined(R__MACOSX)

////////////////////////////////////////////////////////////////////////////////
/// Check if lib is in the dynamic linker cache, returns true if it is, and if so,
/// modifies the library file name parameter `lib` from `/usr/lib/libFOO.dylib`
/// to `-lFOO` such that it can be passed to the linker.
/// This is a unique feature of macOS 11.

static bool R__UpdateLibFileForLinking(TString &lib)
{
   const char *mapfile = nullptr;
#if __x86_64__
   mapfile = "/System/Library/dyld/dyld_shared_cache_x86_64.map";
#elif __arm64__
   mapfile = "/System/Library/dyld/dyld_shared_cache_arm64e.map";
#else
   #error unsupported architecture
#endif
   if (std::ifstream cacheMap{mapfile}) {
      std::string line;
      while (getline(cacheMap, line)) {
         if (line.find(lib) != std::string::npos) {
            lib.ReplaceAll("/usr/lib/lib","-l");
            lib.ReplaceAll(".dylib","");
            return true;
         }
      }
      return false;
   }
   return false;
}
#endif // R__MACOSX

#ifdef R__LINUX

////////////////////////////////////////////////////////////////////////////////
/// Callback for dl_iterate_phdr(), see `man dl_iterate_phdr`.
/// Collects opened libraries.

static int callback_for_dl_iterate_phdr(struct dl_phdr_info *info, size_t size, void *data)
{
   // This function is called through UpdateListOfLoadedSharedLibraries() which is locked.
   static std::unordered_set<decltype(info->dlpi_addr)> sKnownLoadedLibBaseAddrs;

   auto newLibs = static_cast<std::vector<std::string>*>(data);
   if (!sKnownLoadedLibBaseAddrs.count(info->dlpi_addr)) {
      // Skip \0, "", and kernel pseudo-libs linux-vdso.so.1 or linux-gate.so.1
      if (info->dlpi_name && info->dlpi_name[0]
          && strncmp(info->dlpi_name, "linux-vdso.so", 13)
          && strncmp(info->dlpi_name, "linux-vdso32.so", 15)
          && strncmp(info->dlpi_name, "linux-vdso64.so", 15)
          && strncmp(info->dlpi_name, "linux-gate.so", 13))
         newLibs->emplace_back(info->dlpi_name);
      sKnownLoadedLibBaseAddrs.insert(info->dlpi_addr);
   }
   // No matter what the doc says, return != 0 means "stop the iteration".
   return 0;
}

#endif // R__LINUX


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
   // fPrevLoadedDynLibInfo is unused on Linux.
   (void) fPrevLoadedDynLibInfo;

   std::vector<std::string> newLibs;
   dl_iterate_phdr(callback_for_dl_iterate_phdr, &newLibs);
   for (auto &&lib: newLibs)
      RegisterLoadedSharedLibrary(lib.c_str());
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
      DLM->loadLibrary(filename, true /*permanent*/, true /*resolved*/);
   }

#if defined(R__MACOSX)
   // Check that this is not a system library
   auto lenFilename = strlen(filename);
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
       || strstr(filename, "/usr/lib/libextension")
       || strstr(filename, "/usr/lib/libAudioToolboxUtility")
       || strstr(filename, "/usr/lib/liboah")
       || strstr(filename, "/usr/lib/libRosetta")
       // "cannot link directly with dylib/framework, your binary is not an allowed client of
       // /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/
       // SDKs/MacOSX.sdk/usr/lib/libAudioToolboxUtility.tbd for architecture x86_64
       || (lenFilename > 4 && !strcmp(filename + lenFilename - 4, ".tbd")))
      return;
   TString sFileName(filename);
   R__UpdateLibFileForLinking(sFileName);
   filename = sFileName.Data();
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
   assert(!IsFromRootCling() && "Trying to load library from rootcling!");

   // Used to return 0 on success, 1 on duplicate, -1 on failure, -2 on "fatal".
   R__LOCKGUARD_CLING(gInterpreterMutex);
   cling::DynamicLibraryManager* DLM = fInterpreter->getDynamicLibraryManager();
   std::string canonLib = DLM->lookupLibrary(filename);
   cling::DynamicLibraryManager::LoadLibResult res
      = cling::DynamicLibraryManager::kLoadLibNotFound;
   if (!canonLib.empty()) {
      if (system)
         res = DLM->loadLibrary(filename, system, true);
      else {
         // For the non system libs, we'd like to be able to unload them.
         // FIXME: Here we lose the information about kLoadLibAlreadyLoaded case.
         cling::Interpreter::CompilationResult compRes;
         HandleInterpreterException(GetMetaProcessorImpl(), Form(".L %s", canonLib.c_str()), compRes, /*cling::Value*/0);
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

Longptr_t TCling::ProcessLineAsynch(const char* line, EErrorCode* error)
{
   return ProcessLine(line, error);
}

////////////////////////////////////////////////////////////////////////////////
/// Let cling process a command line synchronously, i.e we are waiting
/// it will be finished.

Longptr_t TCling::ProcessLineSynch(const char* line, EErrorCode* error)
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

Longptr_t TCling::Calc(const char* line, EErrorCode* error)
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
   cling::Interpreter::CompilationResult cr = cling::Interpreter::kFailure;
   try {
      cr = fInterpreter->evaluate(line, valRef);
   }
   catch (cling::InterpreterException& ex)
   {
      Error("Calc", "%s.\n%s", ex.what(), "Evaluation of your expression was aborted.");
      ex.diagnose();
      cr = cling::Interpreter::kFailure;
   }

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
   return valRef.simplisticCastAs<Longptr_t>();
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
   cling::Interpreter::PushTransactionRAII RAII(GetInterpreterImpl());
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

static std::string AlternateTuple(const char *classname, const cling::LookupHelper& lh)
{
   TClassEdit::TSplitType tupleContent(classname);
   std::string alternateName = "TEmulatedTuple";
   alternateName.append( classname + 5 );

   std::string fullname = "ROOT::Internal::" + alternateName;
   if (lh.findScope(fullname, cling::LookupHelper::NoDiagnostics,
                    /*resultType*/nullptr, /* intantiateTemplate= */ false))
      return fullname;

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
   // We are shutting down, there is no point in reloading, it only triggers
   // redundant deserializations.
   if (fIsShuttingDown) {
      // Remove the decl_id from the DeclIdToTClass map
      if (cl->fClassInfo) {
         R__LOCKGUARD(gInterpreterMutex);
         TClingClassInfo* TClinginfo = (TClingClassInfo*) cl->fClassInfo;
         // Test again as another thread may have set fClassInfo to nullptr.
         if (TClinginfo) {
            TClass::RemoveClassDeclId(TClinginfo->GetDeclId());
         }
         delete TClinginfo;
         cl->fClassInfo = nullptr;
      }
      return;
   }

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
   if (strncmp(cl->GetName(),"tuple<",strlen("tuple<"))==0) {

      name = AlternateTuple(cl->GetName(), fInterpreter->getLookupHelper());

   }

   bool instantiateTemplate = !cl->TestBit(TClass::kUnloading);
   // FIXME: Rather than adding an option to the TClingClassInfo, we should consider combining code
   // that is currently in the caller (like SetUnloaded) that disable AutoLoading and AutoParsing and
   // code is in the callee (disabling template instantiation) and end up with a more explicit class:
   //      TClingClassInfoReadOnly.
   TClingClassInfo* info = new TClingClassInfo(GetInterpreterImpl(), name.c_str(), instantiateTemplate);
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

   if (fIsShuttingDown || 0 == strncmp(name, anonEnum, cmplen)) {
      return kUnknown;
   }

   // Do not turn on the AutoLoading if it is globally off.
   autoload = autoload && IsClassAutoLoadingEnabled();

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

   // RAII to suspend and restore auto-loading and auto-parsing based on some external conditions.
   class MaybeSuspendAutoLoadParse {
      int fStoreAutoLoad = 0;
      int fStoreAutoParse = 0;
      bool fSuspendedAutoParse = false;
   public:
      MaybeSuspendAutoLoadParse(int autoload) {
         fStoreAutoLoad = ((TCling*)gCling)->SetClassAutoLoading(autoload);
      }

      void SuspendAutoParsing() {
         fSuspendedAutoParse = true;
         fStoreAutoParse = ((TCling*)gCling)->SetSuspendAutoParsing(true);
      }

      ~MaybeSuspendAutoLoadParse() {
         if (fSuspendedAutoParse)
            ((TCling*)gCling)->SetSuspendAutoParsing(fStoreAutoParse);
         ((TCling*)gCling)->SetClassAutoLoading(fStoreAutoLoad);
      }
   };

   MaybeSuspendAutoLoadParse autoLoadParseRAII( autoload );
   if (TClassEdit::IsStdPair(classname) || TClassEdit::IsStdPairBase(classname))
      autoLoadParseRAII.SuspendAutoParsing();

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
               return kUnknown;
            }
         }
      }
      TClingClassInfo tci(GetInterpreterImpl(), *type);
      if (!tci.IsValid()) {
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
                  ROOT::TMetaUtils::GetTrivialIntegralReturnValue(implLineFunc.GetAsFunctionDecl(), *fInterpreter);
               hasClassDefInline = success && (lineNumber == -1);
            }
         }

         // fprintf(stderr,"CheckClassInfo: %s had dict=%d  inline=%d\n",name,hasDictionary.IsValid()
         // , hasClassDefInline);

         // We are now sure that the entry is not in fact an autoload entry.
         if (hasClassDefInline)
            return kWithClassDefInline;
         else
            return kKnown;
      } else {
         // We are now sure that the entry is not in fact an autoload entry.
         return kUnknown;
      }
   }

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
      return kTRUE;
   }
   */

//   const clang::Decl *decl = lh.findScope(name);
//   if (!decl) {
//      std::string buf = TClassEdit::InsertStd(name);
//      decl = lh.findScope(buf);
//   }

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
   TClingBaseClassInfo t(GetInterpreterImpl(), tci);
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
      cling::Interpreter::PushTransactionRAII deserRAII(GetInterpreterImpl());
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
      cling::Interpreter::PushTransactionRAII deserRAII(GetInterpreterImpl());
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
/// Get the scopes representing using declarations of namespace

std::vector<std::string> TCling::GetUsingNamespaces(ClassInfo_t *cl) const
{
   TClingClassInfo *ci = (TClingClassInfo*)cl;
   return ci->GetUsingNamespaces();
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
   TClingMethodArgInfo t(GetInterpreterImpl(), (TClingMethodInfo*)m->fInfo);
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
   R__LOCKGUARD(gInterpreterMutex);
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
         Version_t newvers = ROOT::TMetaUtils::GetClassVersion(llvm::dyn_cast<clang::RecordDecl>(cli->GetDecl()),
                                                     *fInterpreter);
         if (newvers == -1) {
            // Didn't manage to determine the class version from the AST.
            // Use runtime instead.
            if ((mi.Property() & kIsStatic)
                && !fInterpreter->isInSyntaxOnlyMode()) {
               // This better be a static function.
               TClingCallFunc callfunc(GetInterpreterImpl(), *fNormalizedCtxt);
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
                  Sema::ForExternalRedeclaration);

   // Could trigger deserialization of decls.
   cling::Interpreter::PushTransactionRAII RAII(GetInterpreterImpl());
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
            cling::Interpreter::PushTransactionRAII RAII(GetInterpreterImpl());
            possibleEnum = cling::utils::Lookup::Tag(&fInterpreter->getSema(), name, dc);
         } else {
            Error("TCling::GetEnum", "DeclContext not found for %s .\n", name);
         }
      }
   } else {
      // If it is a global enum.
      // Could trigger deserialization of decls.
      cling::Interpreter::PushTransactionRAII RAII(GetInterpreterImpl());
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
         TClingClassInfo gcl(GetInterpreterImpl());
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
      TClingClassInfo cl(GetInterpreterImpl(), scopename.c_str());
      d = cl.GetDataMember(dataname.c_str());
   }
   else {
      TClingClassInfo gcl(GetInterpreterImpl());
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
   TClingCallFunc func(GetInterpreterImpl(), *fNormalizedCtxt);
   if (cl) {
      Longptr_t offset;
      func.SetFunc((TClingClassInfo*)cl->GetClassInfo(), method, params, objectIsConst,
         &offset);
   }
   else {
      TClingClassInfo gcl(GetInterpreterImpl());
      Longptr_t offset;
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
   TClingClassInfo gcl(GetInterpreterImpl());
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
   TClingCallFunc func(GetInterpreterImpl(), *fNormalizedCtxt);
   if (cl) {
      Longptr_t offset;
      func.SetFunc((TClingClassInfo*)cl->GetClassInfo(), method, params, objectIsConst,
                   &offset);
   }
   else {
      TClingClassInfo gcl(GetInterpreterImpl());
      Longptr_t offset;
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
      TClingClassInfo gcl(GetInterpreterImpl());
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

   // NotForRedeclaration: we want to find names in inline namespaces etc.
   clang::LookupResult R(S, DName, clang::SourceLocation(),
                         Sema::LookupOrdinaryName, clang::Sema::NotForRedeclaration);
   R.suppressDiagnostics(); // else lookup with NotForRedeclaration will check access etc
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
      } else if (const auto *USD = llvm::dyn_cast<const clang::UsingShadowDecl>(*IR)) {
         // FIXME: multi-level using
         if (llvm::isa<clang::FunctionDecl>(USD->getTargetDecl())) {
            res.push_back(USD);
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
      TClingClassInfo gcl(GetInterpreterImpl());
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
      TClingClassInfo gcl(GetInterpreterImpl());
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
      TClingClassInfo gcl(GetInterpreterImpl());
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
      TClingClassInfo gcl(GetInterpreterImpl());
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

   TClingClassInfo cl(GetInterpreterImpl(), name);
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
   TClingClassInfo cl(GetInterpreterImpl());
   Longptr_t offset = 0L;
   TClingCallFunc func(GetInterpreterImpl(), *fNormalizedCtxt);
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
   Longptr_t offset = 0L;
   TClingCallFunc func(GetInterpreterImpl(), *fNormalizedCtxt);
   func.SetFunc((TClingClassInfo*)cl->GetClassInfo(), method, params, objectIsConst, &offset);
   void* address = (void*)((Longptr_t)addr + offset);
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
         TClingTypeInfo type(GetInterpreterImpl(), arg->GetFullTypeName());
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
   TClingCallFunc func(GetInterpreterImpl(), *fNormalizedCtxt);
   TClingMethodInfo *minfo = (TClingMethodInfo*)method->fInfo;
   func.Init(*minfo);
   func.SetArgs(listpar);
   // Now calculate the 'this' pointer offset for the method
   // when starting from the class described by cl.
   const CXXMethodDecl * mdecl = dyn_cast<CXXMethodDecl>(minfo->GetTargetFunctionDecl());
   Longptr_t offset = ((TClingClassInfo*)cl->GetClassInfo())->GetOffset(mdecl);
   void* address = (void*)((Longptr_t)addr + offset);
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

Longptr_t TCling::ExecuteMacro(const char* filename, EErrorCode* error)
{
   R__LOCKGUARD_CLING(fLockProcessLine ? gInterpreterMutex : 0);
   fCurExecutingMacros.push_back(filename);
   Longptr_t result = TApplication::ExecuteFile(filename, (int*)error);
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
   TTHREAD_TLS_DECL(std::string,t);

   if (!strstr(typeDesc, "(*)(")) {
      const char *s = strchr(typeDesc, ' ');
      const char *template_start = strchr(typeDesc, '<');
      if (!strcmp(typeDesc, "long long")) {
         t = typeDesc;
      }
      else if (!strncmp(typeDesc, "unsigned ", s + 1 - typeDesc)) {
         t = typeDesc;
      }
      // s is the position of the second 'word' (if any)
      // except in the case of templates where there will be a space
      // just before any closing '>': eg.
      //    TObj<std::vector<UShort_t,__malloc_alloc_template<0> > >*
      else if (s && (template_start == 0 || (s < template_start))) {
         t = s + 1;
      }
      else {
         t = typeDesc;
      }
   }
   else {
      t = typeDesc;
   }
   auto l = t.length();
   while (l > 0 && (t[l - 1] == '*' || t[l - 1] == '&'))
      --l;
   t.resize(l);
   return t.c_str(); // NOLINT
}

static bool requiresRootMap(const char* rootmapfile)
{
   assert(rootmapfile && *rootmapfile);

   llvm::StringRef libName = llvm::sys::path::filename(rootmapfile);
   libName.consume_back(".rootmap");

   return !gInterpreter->HasPCMForLibrary(libName.str().c_str());
}

////////////////////////////////////////////////////////////////////////////////
/// Read and parse a rootmapfile in its new format, and return 0 in case of
/// success, -1 if the file has already been read, and -3 in case its format
/// is the old one (e.g. containing "Library.ClassName"), -4 in case of syntax
/// error.

int TCling::ReadRootmapFile(const char *rootmapfile, TUniqueString *uniqueString)
{
   if (!(rootmapfile && *rootmapfile))
      return 0;

   if (!requiresRootMap(rootmapfile))
      return 0; // success

   // For "class ", "namespace ", "typedef ", "header ", "enum ", "var " respectively
   const std::map<char, unsigned int> keyLenMap = {{'c',6},{'n',10},{'t',8},{'h',7},{'e',5},{'v',4}};

   std::string rootmapfileNoBackslash(rootmapfile);
#ifdef _MSC_VER
   std::replace(rootmapfileNoBackslash.begin(), rootmapfileNoBackslash.end(), '\\', '/');
#endif
   // Add content of a specific rootmap file
   if (fRootmapFiles->FindObject(rootmapfileNoBackslash.c_str()))
      return -1;

   // Line 1 is `{ decls }`
   std::string lineDirective = std::string("\n#line 2 \"Forward declarations from ") + rootmapfileNoBackslash + "\"\n";

   std::ifstream file(rootmapfileNoBackslash);
   std::string line;
   line.reserve(200);
   std::string lib_name;
   line.reserve(100);
   bool newFormat = false;
   while (getline(file, line, '\n')) {
      if (!newFormat && (line.compare(0, 8, "Library.") == 0 || line.compare(0, 8, "Declare.") == 0)) {
         file.close();
         return -3; // old format
      }
      newFormat = true;

      if (line.compare(0, 9, "{ decls }") == 0) {
         // forward declarations

         while (getline(file, line, '\n')) {
            if (line[0] == '[')
               break;
            if (!uniqueString) {
               Error("ReadRootmapFile", "Cannot handle \"{ decls }\" sections in custom rootmap file %s",
                     rootmapfileNoBackslash.c_str());
               return -4;
            }
            if (!lineDirective.empty())
               uniqueString->Append(lineDirective);
            uniqueString->Append(line + '\n');
         }
      }
      const char firstChar = line[0];
      if (firstChar == '[') {
         // new section (library)
         auto brpos = line.find(']');
         if (brpos == string::npos)
            continue;
         lib_name = line.substr(1, brpos - 1);
         size_t nspaces = 0;
         while (lib_name[nspaces] == ' ')
            ++nspaces;
         if (nspaces)
            lib_name.replace(0, nspaces, "");
         if (gDebug > 3) {
            TString lib_nameTstr(lib_name.c_str());
            TObjArray *tokens = lib_nameTstr.Tokenize(" ");
            const char *lib = ((TObjString *)tokens->At(0))->GetName();
            const char *wlib = gSystem->DynamicPathName(lib, kTRUE);
            if (wlib) {
               Info("ReadRootmapFile", "new section for %s", lib_nameTstr.Data());
            } else {
               Info("ReadRootmapFile", "section for %s (library does not exist)", lib_nameTstr.Data());
            }
            delete[] wlib;
            delete tokens;
         }
      } else {
         auto keyLenIt = keyLenMap.find(firstChar);
         if (keyLenIt == keyLenMap.end())
            continue;
         unsigned int keyLen = keyLenIt->second;
         // Do not make a copy, just start after the key
         const char *keyname = line.c_str() + keyLen;
         if (gDebug > 6)
            Info("ReadRootmapFile", "class %s in %s", keyname, lib_name.c_str());
         TEnvRec *isThere = fMapfile->Lookup(keyname);
         if (isThere) {
            if (lib_name != isThere->GetValue()) { // the same key for two different libs
               if (firstChar == 'n') {
                  if (gDebug > 3)
                     Info("ReadRootmapFile", "namespace %s found in %s is already in %s", keyname, lib_name.c_str(),
                          isThere->GetValue());
               } else if (firstChar == 'h') { // it is a header: add the libname to the list of libs to be loaded.
                  lib_name += " ";
                  lib_name += isThere->GetValue();
                  fMapfile->SetValue(keyname, lib_name.c_str());
               } else if (!TClassEdit::IsSTLCont(keyname)) {
                  Warning("ReadRootmapFile", "%s %s found in %s is already in %s", line.substr(0, keyLen).c_str(),
                          keyname, lib_name.c_str(), isThere->GetValue());
               }
            } else { // the same key for the same lib
               if (gDebug > 3)
                  Info("ReadRootmapFile", "Key %s was already defined for %s", keyname, lib_name.c_str());
            }
         } else {
            fMapfile->SetValue(keyname, lib_name.c_str());
         }
      }
   }
   file.close();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Create a resource table and read the (possibly) three resource files,
/// i.e.\ $ROOTSYS/etc/system<name> (or ROOTETCDIR/system<name>), $HOME/<name>
/// and $PWD/<name>. ROOT always reads ".rootrc" (in TROOT::InitSystem()). You
/// can read additional user defined resource files by creating additional TEnv
/// objects. By setting the shell variable ROOTENV_NO_HOME=1 the reading of
/// the $HOME/<name> resource file will be skipped. This might be useful in
/// case the home directory resides on an automounted remote file system
/// and one wants to avoid the file system from being mounted.

void TCling::InitRootmapFile(const char *name)
{
   assert(requiresRootMap(name) && "We have a module!");

   if (!requiresRootMap(name))
      return;

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
      // This class is to be considered an helper for AutoLoading.
      // It is a recursive visitor is used to inspect namespaces and specializations
      // coming from forward declarations in rootmaps and to set the external visible
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
      bool VisitClassTemplateSpecializationDecl(ClassTemplateSpecializationDecl* specDecl) {
         // We want to enable the external lookup for this specialization
         // because we can provide a definition for it!
         if (specDecl->getTemplateSpecializationKind() == TSK_ExplicitSpecialization)
            //SpecSet.insert(specDecl);
            specDecl->setHasExternalLexicalStorage();

         // No need to recurse. On the contrary, recursing is actively harmful:
         // NOTE: must not recurse to prevent this visitor from triggering loading from
         // the external AST source (i.e. autoloading). This would be triggered right here,
         // before autoloading is even set up, as rootmap file parsing happens before that.
         // Even if autoloading is off and has no effect, triggering loading from external
         // AST source resets the flag setHasExternalLexicalStorage(), hiding this specialization
         // from subsequent autoloads!
         return false;
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
   if (rootmapfile && *rootmapfile && !requiresRootMap(rootmapfile))
      return 0;

   R__LOCKGUARD(gInterpreterMutex);

   // open the [system].rootmap files
   if (!fMapfile) {
      fMapfile = new TEnv();
      fMapfile->IgnoreDuplicates(kTRUE);
      fRootmapFiles = new TObjArray;
      fRootmapFiles->SetOwner();
      InitRootmapFile(".rootmap");
   }

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
                           Int_t ret = ReadRootmapFile(p, &uniqueString);

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
      Int_t res = ReadRootmapFile(rootmapfile, &uniqueString);
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

   if (T) {
      ExtVisibleStorageAdder evsAdder(fNSFromRootmaps);
      for (auto declIt = T->decls_begin(); declIt < T->decls_end(); ++declIt) {
         if (declIt->m_DGR.isSingleDecl()) {
            if (Decl* D = declIt->m_DGR.getSingleDecl()) {
               if (clang::isa<TagDecl>(D) || clang::isa<NamespaceDecl>(D)) {
                  evsAdder.TraverseDecl(D);
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
         Error("ReloadAllSharedLibraryMaps", "Could not find rootmap %s in path", rootMapBaseStr.Data());
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
/// Register the AutoLoading information for a class.
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
   assert(IsClassAutoLoadingEnabled() && "Calling when AutoLoading is off!");

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
// Get the list of 'published'/'known' library for the class and load them.
Int_t TCling::ShallowAutoLoadImpl(const char *cls)
{
   Int_t status = 0;

   // lookup class to find list of dependent libraries
   TString deplibs = gCling->GetClassSharedLibs(cls);
   if (!deplibs.IsNull()) {
      TString delim(" ");
      TObjArray* tokens = deplibs.Tokenize(delim);
      for (Int_t i = (tokens->GetEntriesFast() - 1); i > 0; --i) {
         const char* deplib = ((TObjString*)tokens->At(i))->GetName();
         if (gROOT->LoadClass(cls, deplib) == 0) {
            if (gDebug > 0) {
               gCling->Info("TCling::AutoLoad",
                            "loaded dependent library %s for %s", deplib, cls);
            }
         }
         else {
            gCling->Error("TCling::AutoLoad",
                          "failure loading dependent library %s for %s",
                          deplib, cls);
         }
      }
      const char* lib = ((TObjString*)tokens->At(0))->GetName();
      if (lib && lib[0]) {
         if (gROOT->LoadClass(cls, lib) == 0) {
            if (gDebug > 0) {
               gCling->Info("TCling::AutoLoad",
                            "loaded library %s for %s", lib, cls);
            }
            status = 1;
         }
         else {
            gCling->Error("TCling::AutoLoad",
                          "failure loading library %s for %s", lib, cls);
         }
      }
      delete tokens;
   }

   return status;
}

////////////////////////////////////////////////////////////////////////////////
// Iterate through the data member of the class (either through the TProtoClass
// or through Cling) and trigger, recursively, the loading the necessary libraries.
// \note `cls` is expected to be already normalized!
// \returns 1 on success.
Int_t TCling::DeepAutoLoadImpl(const char *cls, std::unordered_set<std::string> &visited,
                               bool nameIsNormalized)
{
   // Try to insert; if insertion failed because the entry existed, DeepAutoLoadImpl()
   // has previously (within the same call to `AutoLoad()`) tried to load this class
   // and we are done, whether success or not, as it won't work better now than before,
   // because there is no additional information now compared to before.
   if (!visited.insert(std::string(cls)).second)
      return 1;

   if (ShallowAutoLoadImpl(cls) == 0) {
      // If ShallowAutoLoadImpl() has an error, we have an error.
      return 0;
   }

   // Now look through the TProtoClass to load the required library/dictionary
   if (TProtoClass *proto = nameIsNormalized ? TClassTable::GetProtoNorm(cls) : TClassTable::GetProto(cls)) {
      for (auto element : proto->GetData()) {
         if (element->IsBasic())
            continue;
         const char *subtypename = element->GetTypeName();
         if (!TClassTable::GetDictNorm(subtypename)) {
            // Failure to load a dictionary is not (quite) a failure load
            // the top-level library.  If we return false here, then
            // we would end up in a situation where the library and thus
            // the dictionary is loaded for "cls" but the TClass is
            // not created and/or marked as unavailable (in case where
            // AutoLoad is called from TClass::GetClass).
            DeepAutoLoadImpl(subtypename, visited, true /*normalized*/);
         }
      }
      return 1;
   }

   // We found no TProtoClass for cls.
   auto classinfo = gInterpreter->ClassInfo_Factory(cls);
   if (classinfo && gInterpreter->ClassInfo_IsValid(classinfo)
         && !(gInterpreter->ClassInfo_Property(classinfo) & kIsEnum))
   {
      DataMemberInfo_t *memberinfo = gInterpreter->DataMemberInfo_Factory(classinfo, TDictionary::EMemberSelection::kNoUsingDecls);
      while (gInterpreter->DataMemberInfo_Next(memberinfo)) {
         if (gInterpreter->DataMemberInfo_TypeProperty(memberinfo) & ::kIsFundamental)
            continue;
         auto membertypename = TClassEdit::GetLong64_Name(gInterpreter->TypeName(gInterpreter->DataMemberInfo_TypeTrueName(memberinfo)));
         if (!TClassTable::GetDictNorm(membertypename.c_str())) {
            // Failure to load a dictionary is not (quite) a failure load
            // the top-level library.   See detailed comment in the TProtoClass
            // branch (above).
            (void)DeepAutoLoadImpl(membertypename.c_str(), visited, true /*normalized*/);
         }
      }
      gInterpreter->DataMemberInfo_Delete(memberinfo);
   }
   gInterpreter->ClassInfo_Delete(classinfo);
   return 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Load library containing the specified class. Returns 0 in case of error
/// and 1 in case if success.

Int_t TCling::AutoLoad(const char *cls, Bool_t knowDictNotLoaded /* = kFALSE */)
{
   // Prevent update to IsClassAutoloading between our check and our actions.
   R__READ_LOCKGUARD(ROOT::gCoreMutex);

   // TClass::GetClass explicitly calls gInterpreter->AutoLoad. When called from
   // rootcling (in *_rdict.pcm file generation) it is a no op.
   // FIXME: We should avoid calling autoload when we know we are not supposed
   // to and transform this check into an assert.
   if (!IsClassAutoLoadingEnabled()) {
      // Never load any library from rootcling/genreflex.
      if (gDebug > 2) {
         Info("TCling::AutoLoad", "Explicitly disabled (the class name is %s)", cls);
      }
      return 0;
   }

   assert(IsClassAutoLoadingEnabled() && "Calling when AutoLoading is off!");

   R__WRITE_LOCKGUARD(ROOT::gCoreMutex);

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

   if (!gROOT || !gInterpreter || gROOT->TestBit(TObject::kInvalidObject)) {
      if (gDebug > 2) {
         Info("TCling::AutoLoad",
              "Disabled due to gROOT or gInterpreter being invalid/not ready (the class name is %s)", cls);
      }
      return 0;
   }
   // Prevent the recursion when the library dictionary are loaded.
   SuspendAutoLoadingRAII autoLoadOff(this);
   // Try using externally provided callback first.
   if (fAutoLoadCallBack) {
      int success = (*(AutoLoadCallBack_t)fAutoLoadCallBack)(cls);
      if (success)
         return success;
   }

   // During the 'Deep' part of the search we will call GetClassSharedLibsForModule
   // (when module are enabled) which might end up calling AutoParsing but
   // that should only be for the cases where the library has no generated pcm
   // and in that case a rootmap should be available.
   // This avoids a very costly operation (for generally no gain) but reduce the
   // quality of the search (i.e. bad in case of library with no pcm and no rootmap
   // file).
   TInterpreter::SuspendAutoParsing autoParseRaii(this);
   std::unordered_set<std::string> visited;
   return DeepAutoLoadImpl(cls, visited, false /*normalized*/);
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
         auto nTokens = tokens->GetEntriesFast();
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
                  auto cRes = ExecAutoParse(hName, kFALSE, GetInterpreterImpl());
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
                  auto cRes = ExecAutoParse(hName, kTRUE, GetInterpreterImpl());
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
   if (llvm::StringRef(cls).contains("(lambda)"))
      return 0;

   if (!fHeaderParsingOnDemand || fIsAutoParsingSuspended) {
      if (fClingCallbacks->IsAutoLoadingEnabled()) {
         return AutoLoad(cls);
      } else {
         return 0;
      }
   }

   R__LOCKGUARD(gInterpreterMutex);

   if (gDebug > 1) {
      Info("TCling::AutoParse",
           "Trying to autoparse for %s", cls);
   }

   // The catalogue of headers is in the dictionary
   if (fClingCallbacks->IsAutoLoadingEnabled()
         && !gClassTable->GetDictNorm(cls)) {
      // Need RAII against recursive (dictionary payload) parsing (ROOT-8445).
      ROOT::Internal::ParsingStateRAII parsingStateRAII(fInterpreter->getParser(),
         fInterpreter->getSema());
      AutoLoad(cls, true /*knowDictNotLoaded*/);
   }

   // Prevent the recursion when the library dictionary are loaded.
   SuspendAutoLoadingRAII autoLoadOff(this);

   // No recursive header parsing on demand; we require headers to be standalone.
   SuspendAutoParsing autoParseRAII(this);

   Int_t nHheadersParsed = AutoParseImplRecurse(cls,/*topLevel=*/ true);

   ProcessClassesToUpdate();

   return nHheadersParsed > 0 ? 1 : 0;
}

// This is a function which gets callback from cling when DynamicLibraryManager->loadLibrary failed for some reason.
// Try to solve the problem by AutoLoading. Return true when AutoLoading success, return
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

static void* LazyFunctionCreatorAutoloadForModule(const std::string &mangled_name,
                                                  const cling::DynamicLibraryManager &DLM) {
   R__LOCKGUARD(gInterpreterMutex);

   auto LibLoader = [](const std::string& LibName) -> bool {
     if (gSystem->Load(LibName.c_str(), "", false) < 0) {
       Error("TCling__LazyFunctionCreatorAutoloadForModule",
             "Failed to load library %s", LibName.c_str());
       return false;
     }
     return true; //success.
   };

#ifdef R__MACOSX
   // The JIT gives us a mangled name which has only one leading underscore on
   // all platforms, for instance _ZN8TRandom34RndmEv. However, on OSX the
   // linker stores this symbol as __ZN8TRandom34RndmEv (adding an extra _).
   assert(!llvm::StringRef(mangled_name).startswith("__") && "Already added!");
   std::string libName = DLM.searchLibrariesForSymbol('_' + mangled_name,
                                                      /*searchSystem=*/ true);
#else
   std::string libName = DLM.searchLibrariesForSymbol(mangled_name,
                                                      /*searchSystem=*/ true);
#endif //R__MACOSX

   assert(!llvm::StringRef(libName).startswith("libNew") &&
          "We must not resolve symbols from libNew!");

   if (libName.empty())
      return nullptr;

   if (!LibLoader(libName))
      return nullptr;

   return llvm::sys::DynamicLibrary::SearchForAddressOfSymbol(mangled_name);

}

////////////////////////////////////////////////////////////////////////////////
/// Autoload a library based on a missing symbol.

void* TCling::LazyFunctionCreatorAutoload(const std::string& mangled_name) {
   if (fCxxModulesEnabled)
      return LazyFunctionCreatorAutoloadForModule(mangled_name,
                                                  *GetInterpreterImpl()->getDynamicLibraryManager());

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

Bool_t TCling::IsAutoLoadNamespaceCandidate(const clang::NamespaceDecl* nsDecl)
{
   return fNSFromRootmaps.count(nsDecl) != 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Internal function. Actually do the update of the ClassInfo when seeing
//  new TagDecl or NamespaceDecl.
void TCling::RefreshClassInfo(TClass *cl, const clang::NamedDecl *def, bool alias)
{

   TClingClassInfo *cci = ((TClingClassInfo *)cl->fClassInfo);
   if (cci) {
      // If we only had a forward declaration then update the
      // TClingClassInfo with the definition if we have it now.
      const NamedDecl *oldDef = llvm::dyn_cast_or_null<NamedDecl>(cci->GetDecl());
      if (!oldDef || (def && def != oldDef)) {
         cl->ResetCaches();
         TClass::RemoveClassDeclId(cci->GetDeclId());
         if (def) {
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
      if (!alias && def != nullptr)
         cl->fClassInfo = (ClassInfo_t *)new TClingClassInfo(GetInterpreterImpl(), def);
      else
         cl->fClassInfo = (ClassInfo_t *)new TClingClassInfo(GetInterpreterImpl(), cl->GetName());
      if (((TClingClassInfo *)cl->fClassInfo)->IsValid()) {
         // We now need to update the state and bits.
         if (cl->fState != TClass::kHasTClassInit) {
            // if (!cl->fClassInfo->IsValid()) cl->fState = TClass::kForwardDeclared; else
            cl->fState = TClass::kInterpreted;
            cl->ResetBit(TClass::kIsEmulation);
         }
         TClass::AddClassToDeclIdMap(((TClingClassInfo *)(cl->fClassInfo))->GetDeclId(), cl);
      } else {
         delete ((TClingClassInfo *)cl->fClassInfo);
         cl->fClassInfo = nullptr;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Internal function. Inform a TClass about its new TagDecl or NamespaceDecl.
void TCling::UpdateClassInfoWithDecl(const NamedDecl* ND)
{
   const TagDecl *td = dyn_cast<TagDecl>(ND);
   const NamespaceDecl *ns = dyn_cast<NamespaceDecl>(ND);
   const NamedDecl *canon = nullptr;

   std::string name;
   TagDecl* tdDef = 0;
   if (td) {
      canon = tdDef = td->getDefinition();
      // Let's pass the decl to the TClass only if it has a definition.
      if (!tdDef) return;

      if (!tdDef->isCompleteDefinition() || llvm::isa<clang::FunctionDecl>(tdDef->getDeclContext())) {
         // Ignore incomplete definition.
         // Ignore declaration within a function.
         return;
      }

      auto declName = tdDef->getNameAsString();
      // Check if we have registered the unqualified name into the list
      // of TClass that are in kNoInfo, kEmulated or kFwdDeclaredState.
      // Since this is used as heureutistic to avoid spurrious calls to GetNormalizedName
      // the unqualified name is sufficient (and the fully qualified name might be
      // 'wrong' if there is difference in spelling in the template paramters (for example)
      if (!TClass::HasNoInfoOrEmuOrFwdDeclaredDecl(declName.c_str())){
         // fprintf (stderr,"WARNING: Impossible to find a TClassEntry in kNoInfo or kEmulated the decl of which would be called %s. Skip w/o building the normalized name.\n",declName.c_str() );
         return;
      }

      clang::QualType type(tdDef->getTypeForDecl(), 0);
      ROOT::TMetaUtils::GetNormalizedName(name, type, *fInterpreter, *fNormalizedCtxt);
   } else if (ns) {
      canon = ns->getCanonicalDecl();
      name = ND->getQualifiedNameAsString();
   } else {
      name = ND->getQualifiedNameAsString();
   }

   // Supposedly we are being called while something is being
   // loaded ... let's now tell the autoloader to do the work
   // yet another time.
   SuspendAutoLoadingRAII autoLoadOff(this);
   // FIXME: There can be more than one TClass for a single decl.
   // for example vector<double> and vector<Double32_t>
   TClass* cl = (TClass*)gROOT->GetListOfClasses()->FindObject(name.c_str());
   if (cl && GetModTClasses().find(cl) == GetModTClasses().end()) {
      RefreshClassInfo(cl, canon, false);
   }
   // And here we should find the other 'aliases' (eg. vector<Double32_t>)
   // and update them too:
   // foreach(aliascl in gROOT->GetListOfClasses()->FindAliasesOf(name.c_str()))
   //    RefreshClassInfo(cl, tdDef, true);
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
      cling::Interpreter::PushTransactionRAII RAII(GetInterpreterImpl());
      // Unlock the TClass for updates
      ((TCling*)gCling)->GetModTClasses().erase(*I);

   }
}

///\brief Invalidate stored TCling state for declarations included in transaction `T'.
///
void TCling::UpdateListsOnUnloaded(const cling::Transaction &T)
{
   HandleNewTransaction(T);

   auto Lists = std::make_tuple((TListOfDataMembers *)gROOT->GetListOfGlobals(),
                                (TListOfFunctions *)gROOT->GetListOfGlobalFunctions(),
                                (TListOfFunctionTemplates *)gROOT->GetListOfFunctionTemplates(),
                                (TListOfEnums *)gROOT->GetListOfEnums());

   cling::Transaction::const_nested_iterator iNested = T.nested_begin();
   for (cling::Transaction::const_iterator I = T.decls_begin(), E = T.decls_end();
        I != E; ++I) {
      if (I->m_Call == cling::Transaction::kCCIHandleVTable)
         continue;
      if (I->m_Call == cling::Transaction::kCCINone) {
         UpdateListsOnUnloaded(*(*iNested));
         ++iNested;
         continue;
      }

      for (auto &D : I->m_DGR)
         InvalidateCachedDecl(Lists, D);
   }
}

///\brief Invalidate cached TCling information for the given global declaration.
///
void TCling::InvalidateGlobal(const Decl *D) {
   auto Lists = std::make_tuple((TListOfDataMembers *)gROOT->GetListOfGlobals(),
                                (TListOfFunctions *)gROOT->GetListOfGlobalFunctions(),
                                (TListOfFunctionTemplates *)gROOT->GetListOfFunctionTemplates(),
                                (TListOfEnums *)gROOT->GetListOfEnums());
   InvalidateCachedDecl(Lists, D);
}

///\brief Invalidate cached TCling information for the given declaration, and
/// removed it from the appropriate object list.
///\param[in] Lists - std::tuple<TListOfDataMembers&, TListOfFunctions&,
///                              TListOfFunctionTemplates&, TListOfEnums&>
///                   of pointers to the (global/class) object lists.
///\param[in] D - Decl to discard.
///
void TCling::InvalidateCachedDecl(const std::tuple<TListOfDataMembers*,
                                             TListOfFunctions*,
                                             TListOfFunctionTemplates*,
                                             TListOfEnums*> &Lists, const Decl *D) {
   if (D->isFromASTFile())  // `D' came from the PCH; ignore
      return;

   TListOfDataMembers &LODM = *(std::get<0>(Lists));
   TListOfFunctions &LOF = *(std::get<1>(Lists));
   TListOfFunctionTemplates &LOFT = *(std::get<2>(Lists));
   TListOfEnums &LOE = *(std::get<3>(Lists));

   if (isa<VarDecl>(D) || isa<FieldDecl>(D) || isa<EnumConstantDecl>(D)) {
      TObject *O = LODM.Find((TDictionary::DeclId_t)D);
      if (LODM.GetClass())
         RemoveAndInvalidateObject(LODM, static_cast<TDataMember *>(O));
      else
         RemoveAndInvalidateObject(LODM, static_cast<TGlobal *>(O));
   } else if (isa<FunctionDecl>(D)) {
      RemoveAndInvalidateObject(LOF, LOF.Find((TDictionary::DeclId_t)D));
   } else if (isa<FunctionTemplateDecl>(D)) {
      RemoveAndInvalidateObject(LOFT, LOFT.Get((TDictionary::DeclId_t)D));
   } else if (isa<EnumDecl>(D)) {
      TEnum *E = LOE.Find((TDictionary::DeclId_t)D);
      if (!E)
         return;

      // Try to invalidate enumerators (for unscoped enumerations).
      for (TIter I = E->GetConstants(); auto EC = (TEnumConstant *)I(); )
         RemoveAndInvalidateObject(LODM,
      	       	         	   (TEnumConstant *)LODM.FindObject(EC->GetName()));

      RemoveAndInvalidateObject(LOE, E);
   } else if (isa<RecordDecl>(D) || isa<NamespaceDecl>(D)) {
      if (isa<RecordDecl>(D) && !cast<RecordDecl>(D)->isCompleteDefinition())
         return;

      std::vector<TClass *> Classes;
      if (!TClass::GetClass(D->getCanonicalDecl(), Classes))
         return;
      for (auto &C : Classes) {
         auto Lists = std::make_tuple((TListOfDataMembers *)C->GetListOfDataMembers(),
                                      (TListOfFunctions *)C->GetListOfMethods(),
                                      (TListOfFunctionTemplates *)C->GetListOfFunctionTemplates(),
                                      (TListOfEnums *)C->GetListOfEnums());
         for (auto &I : cast<DeclContext>(D)->decls())
            InvalidateCachedDecl(Lists, I);

         // For NamespaceDecl (redeclarable), only invalidate this redecl.
         if (D->getKind() != Decl::Namespace
             || cast<NamespaceDecl>(D)->isOriginalNamespace())
            C->ResetClassInfo();
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
// R__LOCKGUARD_CLING(gInterpreterMutex);
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
   R__LOCKGUARD_CLING(gInterpreterMutex);
   UpdateListOfLoadedSharedLibraries();
   return fSharedLibs;
}

static std::string GetClassSharedLibsForModule(const char *cls, cling::LookupHelper &LH)
{
   if (!cls || !*cls)
      return {};

   using namespace clang;
   if (const Decl *D = LH.findScope(cls, cling::LookupHelper::NoDiagnostics,
                                    /*type*/ nullptr, /*instantiate*/ false)) {
      if (!D->isFromASTFile()) {
         if (gDebug > 5)
            Warning("GetClassSharedLibsForModule", "Decl found for %s is not part of a module", cls);
         return {};
      }
      class ModuleCollector : public ConstDeclVisitor<ModuleCollector> {
         llvm::DenseSet<Module *> &m_TopLevelModules;

      public:
         ModuleCollector(llvm::DenseSet<Module *> &TopLevelModules) : m_TopLevelModules(TopLevelModules) {}
         void Collect(const Decl *D) { Visit(D); }

         void VisitDecl(const Decl *D)
         {
            // FIXME: Such case is described ROOT-7765 where
            // ROOT_GENERATE_DICTIONARY does not contain the list of headers.
            // They are specified as #includes in the LinkDef file. This leads to
            // generation of incomplete modulemap files and this logic fails to
            // compute the corresponding module of D.
            // FIXME: If we want to support such a case, we should not rely on
            // the contents of the modulemap but mangle D and look it up in the
            // .so files.
            if (!D->hasOwningModule())
               return;
            if (Module *M = D->getOwningModule()->getTopLevelModule())
               m_TopLevelModules.insert(M);
         }

         void VisitTemplateArgument(const TemplateArgument &TA)
         {
            switch (TA.getKind()) {
            case TemplateArgument::Null:
            case TemplateArgument::Integral:
            case TemplateArgument::Pack:
            case TemplateArgument::NullPtr:
            case TemplateArgument::Expression:
            case TemplateArgument::Template:
            case TemplateArgument::TemplateExpansion: return;
            case TemplateArgument::Type:
               if (const TagType *TagTy = dyn_cast<TagType>(TA.getAsType()))
                  return Visit(TagTy->getDecl());
               return;
            case TemplateArgument::Declaration: return Visit(TA.getAsDecl());
            }
            llvm_unreachable("Invalid TemplateArgument::Kind!");
         }

         void VisitClassTemplateSpecializationDecl(const ClassTemplateSpecializationDecl *CTSD)
         {
            if (CTSD->getOwningModule())
               VisitDecl(CTSD);
            else
               VisitDecl(CTSD->getSpecializedTemplate());
            const TemplateArgumentList &ArgList = CTSD->getTemplateArgs();
            for (const TemplateArgument *Arg = ArgList.data(), *ArgEnd = Arg + ArgList.size(); Arg != ArgEnd; ++Arg) {
               VisitTemplateArgument(*Arg);
            }
         }
      };

      llvm::DenseSet<Module *> TopLevelModules;
      ModuleCollector m(TopLevelModules);
      m.Collect(D);
      std::string result;
      for (auto M : TopLevelModules) {
         // ROOT-unaware modules (i.e. not processed by rootcling) do not have a
         // link declaration.
         if (!M->LinkLibraries.size())
            continue;
         // We have preloaded the Core module thus libCore.so
         if (M->Name == "Core")
            continue;
         assert(M->LinkLibraries.size() == 1);
         if (!result.empty())
            result += ' ';
         result += M->LinkLibraries[0].Library;
      }
      return result;
   }
   return {};
}

////////////////////////////////////////////////////////////////////////////////
/// Get the list of shared libraries containing the code for class cls.
/// The first library in the list is the one containing the class, the
/// others are the libraries the first one depends on. Returns 0
/// in case the library is not found.

const char* TCling::GetClassSharedLibs(const char* cls)
{
   if (fCxxModulesEnabled) {
      llvm::StringRef className = cls;
      // If we get a class name containing lambda, we cannot parse it and we
      // can exit early.
      // FIXME: This works around a bug when we are instantiating a template
      // make_unique and the substitution fails. Seen in most of the dataframe
      // tests.
      if (className.contains("(lambda)"))
         return nullptr;
      // Limit the recursion which can be induced by GetClassSharedLibsForModule.
      SuspendAutoLoadingRAII AutoLoadingDisabled(this);
      cling::LookupHelper &LH = fInterpreter->getLookupHelper();
      std::string libs = GetClassSharedLibsForModule(cls, LH);
      if (!libs.empty()) {
         fAutoLoadLibStorage.push_back(libs);
         return fAutoLoadLibStorage.back().c_str();
      }
   }

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

/// This interface returns a list of dependent libraries in the form:
/// lib libA.so libB.so libC.so. The first library is the library we are
/// searching dependencies for.
/// Note: In order to speed up the search, we display the dependencies of the
/// libraries which are not yet loaded. For instance, if libB.so was already
/// loaded the list would contain: lib libA.so libC.so.
static std::string GetSharedLibImmediateDepsSlow(std::string lib,
                                                 cling::Interpreter *interp,
                                                 bool skipLoadedLibs = true)
{
   TString LibFullPath(lib);
   if (!llvm::sys::path::is_absolute(lib)) {
      if (!gSystem->FindDynamicLibrary(LibFullPath, /*quiet=*/true)) {
         Error("TCling__GetSharedLibImmediateDepsSlow", "Cannot find library '%s'", lib.c_str());
         return "";
      }
   } else {
      assert(llvm::sys::fs::exists(lib) && "Must exist!");
      lib = llvm::sys::path::filename(lib);
   }

   auto ObjF = llvm::object::ObjectFile::createObjectFile(LibFullPath.Data());
   if (!ObjF) {
      Warning("TCling__GetSharedLibImmediateDepsSlow", "Failed to read object file %s", lib.c_str());
      return "";
   }

   llvm::object::ObjectFile *BinObjFile = ObjF.get().getBinary();

   std::set<string> DedupSet;
   std::string Result = lib + ' ';
   for (const auto &S : BinObjFile->symbols()) {
      uint32_t Flags = S.getFlags();
      if (Flags & llvm::object::SymbolRef::SF_Undefined) {
         llvm::Expected<StringRef> SymNameErr = S.getName();
         if (!SymNameErr) {
            Warning("GetSharedLibDepsForModule", "Failed to read symbol");
            continue;
         }
         llvm::StringRef SymName = SymNameErr.get();
         if (SymName.empty())
            continue;

         if (BinObjFile->isELF()) {
            // Skip the symbols which are part of the C/C++ runtime and have a
            // fixed library version. See binutils ld VERSION. Those reside in
            // 'system' libraries, which we avoid in FindLibraryForSymbol.
            if (SymName.contains("@@GLIBCXX") || SymName.contains("@@CXXABI") ||
               SymName.contains("@@GLIBC") || SymName.contains("@@GCC"))
               continue;

            // Those are 'weak undefined' symbols produced by gcc. We can
            // ignore them.
            // FIXME: It is unclear whether we can ignore all weak undefined
            // symbols:
            // http://lists.llvm.org/pipermail/llvm-dev/2017-October/118177.html
            if (SymName == "_Jv_RegisterClasses" ||
               SymName == "_ITM_deregisterTMCloneTable" ||
               SymName == "_ITM_registerTMCloneTable")
               continue;
         }

// FIXME: this might really depend on MachO library format instead of R__MACOSX.
#ifdef R__MACOSX
         // MacOS symbols sometimes have an extra "_", see
         // https://developer.apple.com/library/archive/documentation/System/Conceptual/ManPages_iPhoneOS/man3/dlsym.3.html
         if (skipLoadedLibs && SymName[0] == '_'
             && llvm::sys::DynamicLibrary::SearchForAddressOfSymbol(SymName.drop_front()))
            continue;
#endif

         // If we can find the address of the symbol, we have loaded it. Skip.
         if (skipLoadedLibs && llvm::sys::DynamicLibrary::SearchForAddressOfSymbol(SymName))
            continue;

         R__LOCKGUARD(gInterpreterMutex);
         std::string found = interp->getDynamicLibraryManager()->searchLibrariesForSymbol(SymName, /*searchSystem*/false);
         // The expected output is just filename without the full path, which
         // is not very accurate, because our Dyld implementation might find
         // a match in location a/b/c.so and if we return just c.so ROOT might
         // resolve it to y/z/c.so and there we might not be ABI compatible.
         // FIXME: Teach the users of GetSharedLibDeps to work with full paths.
         if (!found.empty()) {
            std::string cand = llvm::sys::path::filename(found).str();
            if (!DedupSet.insert(cand).second)
               continue;

            Result += cand + ' ';
         }
      }
   }

   return Result;
}

static bool hasParsedRootmapForLibrary(llvm::StringRef lib)
{
   // Check if we have parsed a rootmap file.
   llvm::SmallString<256> rootmapName;
   if (!lib.startswith("lib"))
      rootmapName.append("lib");

   rootmapName.append(llvm::sys::path::filename(lib));
   llvm::sys::path::replace_extension(rootmapName, "rootmap");

   if (gCling->GetRootMapFiles()->FindObject(rootmapName.c_str()))
      return true;

   // Perform a last resort by dropping the lib prefix.
   llvm::StringRef rootmapNameNoLib = rootmapName.str();
   if (rootmapNameNoLib.consume_front("lib"))
      return gCling->GetRootMapFiles()->FindObject(rootmapNameNoLib.data());

   return false;
}

static bool hasPrecomputedLibraryDeps(llvm::StringRef lib)
{
   if (gCling->HasPCMForLibrary(lib.data()))
      return true;

   return hasParsedRootmapForLibrary(lib);
}

////////////////////////////////////////////////////////////////////////////////
/// Get the list a libraries on which the specified lib depends. The
/// returned string contains as first element the lib itself.
/// Returns 0 in case the lib does not exist or does not have
/// any dependencies. If useDyld is true, we iterate through all available
/// libraries and try to construct the dependency chain by resolving each
/// symbol.

const char* TCling::GetSharedLibDeps(const char* lib, bool useDyld/* = false*/)
{
   if (llvm::sys::path::is_absolute(lib) && !llvm::sys::fs::exists(lib))
      return nullptr;

   if (!hasParsedRootmapForLibrary(lib)) {
      llvm::SmallString<512> rootmapName(lib);
      llvm::sys::path::replace_extension(rootmapName, "rootmap");
      if (llvm::sys::fs::exists(rootmapName)) {
         if (gDebug > 0)
            Info("Load", "loading %s", rootmapName.c_str());
         gInterpreter->LoadLibraryMap(rootmapName.c_str());
      }
   }

   if (hasPrecomputedLibraryDeps(lib) && useDyld) {
      if (gDebug > 0)
         Warning("TCling::GetSharedLibDeps", "Precomputed dependencies available but scanning '%s'", lib);
   }

   if (useDyld) {
      std::string libs = GetSharedLibImmediateDepsSlow(lib, GetInterpreterImpl());
      if (!libs.empty()) {
         fAutoLoadLibStorage.push_back(libs);
         return fAutoLoadLibStorage.back().c_str();
      }
   }

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
   HandleInterpreterException(GetMetaProcessorImpl(), TString::Format(".L %s", path), compRes, /*cling::Value*/0);
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
   return buffer.c_str(); // NOLINT
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
/// Returns if class AutoLoading is currently enabled.

bool TCling::IsClassAutoLoadingEnabled() const
{
   if (IsFromRootCling())
      return false;
   if (!fClingCallbacks)
      return false;
   return fClingCallbacks->IsAutoLoadingEnabled();
}

////////////////////////////////////////////////////////////////////////////////
/// Enable/Disable the AutoLoading of libraries.
/// Returns the old value, i.e whether it was enabled or not.

int TCling::SetClassAutoLoading(int autoload) const
{
   // If no state change is required, exit early.
   // FIXME: In future we probably want to complain if we made a request which
   // was with the same state as before in order to catch programming errors.
   if ((bool) autoload == IsClassAutoLoadingEnabled())
      return autoload;

   assert(fClingCallbacks && "We must have callbacks!");
   bool oldVal = fClingCallbacks->IsAutoLoadingEnabled();
   fClingCallbacks->SetAutoLoadingEnabled(autoload);
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
   HandleInterpreterException(GetMetaProcessorImpl(), Form(".U %s", canonical.c_str()), compRes, /*cling::Value*/0);
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

Longptr_t TCling::CallFunc_ExecInt(CallFunc_t* func, void* address) const
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
   return (CallFunc_t*) new TClingCallFunc(GetInterpreterImpl(), *fNormalizedCtxt);
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

void TCling::CallFunc_SetArgArray(CallFunc_t* func, Longptr_t* paramArr, Int_t nparam) const
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

void TCling::CallFunc_SetFunc(CallFunc_t* func, ClassInfo_t* info, const char* method, const char* params, Longptr_t* offset) const
{
   TClingCallFunc* f = (TClingCallFunc*) func;
   TClingClassInfo* ci = (TClingClassInfo*) info;
   f->SetFunc(ci, method, params, offset);
}

////////////////////////////////////////////////////////////////////////////////

void TCling::CallFunc_SetFunc(CallFunc_t* func, ClassInfo_t* info, const char* method, const char* params, bool objectIsConst, Longptr_t* offset) const
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

void TCling::CallFunc_SetFuncProto(CallFunc_t* func, ClassInfo_t* info, const char* method, const char* proto, Longptr_t* offset, EFunctionMatchMode mode /* = kConversionMatch */) const
{
   TClingCallFunc* f = (TClingCallFunc*) func;
   TClingClassInfo* ci = (TClingClassInfo*) info;
   f->SetFuncProto(ci, method, proto, offset, mode);
}

////////////////////////////////////////////////////////////////////////////////
/// Interface to cling function

void TCling::CallFunc_SetFuncProto(CallFunc_t* func, ClassInfo_t* info, const char* method, const char* proto, bool objectIsConst, Longptr_t* offset, EFunctionMatchMode mode /* = kConversionMatch */) const
{
   TClingCallFunc* f = (TClingCallFunc*) func;
   TClingClassInfo* ci = (TClingClassInfo*) info;
   f->SetFuncProto(ci, method, proto, objectIsConst, offset, mode);
}

////////////////////////////////////////////////////////////////////////////////
/// Interface to cling function

void TCling::CallFunc_SetFuncProto(CallFunc_t* func, ClassInfo_t* info, const char* method, const std::vector<TypeInfo_t*> &proto, Longptr_t* offset, EFunctionMatchMode mode /* = kConversionMatch */) const
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

void TCling::CallFunc_SetFuncProto(CallFunc_t* func, ClassInfo_t* info, const char* method, const std::vector<TypeInfo_t*> &proto, bool objectIsConst, Longptr_t* offset, EFunctionMatchMode mode /* = kConversionMatch */) const
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
   if (!declid)
      return kFALSE;

   const clang::DeclContext *ctxt = nullptr;
   if (info) {
      ctxt = clang::Decl::castToDeclContext(((TClingClassInfo*)info)->GetDecl());
   } else {
      ctxt = fInterpreter->getCI()->getASTContext().getTranslationUnitDecl();
   }
   if (!ctxt)
      return kFALSE;

   const clang::Decl *decl = reinterpret_cast<const clang::Decl*>(declid);
   if (!decl)
      return kFALSE;

   const clang::DeclContext *declDC = decl->getDeclContext();
   // ClassInfo_t-s are always "spellable" scopes, never unnamed or inline ones.
   while (true) {
      if (declDC->isTransparentContext()) {
         declDC = declDC->getParent();
         continue;
      }
      if (const auto *declRD = llvm::dyn_cast<clang::RecordDecl>(declDC)) {
         if (declRD->isAnonymousStructOrUnion()) {
            declDC = declRD->getParent();
            continue;
         }
      }
      if (const auto *declNS = llvm::dyn_cast<clang::NamespaceDecl>(declDC)) {
         if (declNS->isAnonymousNamespace() || declNS->isInlineNamespace()) {
            declDC = declNS->getParent();
            continue;
         }
      }
      break;
   }

   return declDC->Equals(ctxt);
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
   return (ClassInfo_t*) new TClingClassInfo(GetInterpreterImpl(), all);
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
   return (ClassInfo_t*) new TClingClassInfo(GetInterpreterImpl(), name);
}

ClassInfo_t* TCling::ClassInfo_Factory(DeclId_t declid) const
{
   R__LOCKGUARD(gInterpreterMutex);
   return (ClassInfo_t*) new TClingClassInfo(GetInterpreterImpl(), (const clang::Decl*)declid);
}


////////////////////////////////////////////////////////////////////////////////

int TCling::ClassInfo_GetMethodNArg(ClassInfo_t* cinfo, const char* method, const char* proto, Bool_t objectIsConst /* = false */, EFunctionMatchMode mode /* = kConversionMatch */) const
{
   TClingClassInfo* TClinginfo = (TClingClassInfo*) cinfo;
   return TClinginfo->GetMethodNArg(method, proto, objectIsConst, mode);
}

////////////////////////////////////////////////////////////////////////////////

bool TCling::ClassInfo_HasDefaultConstructor(ClassInfo_t* cinfo, Bool_t testio) const
{
   TClingClassInfo *TClinginfo = (TClingClassInfo *) cinfo;
   return TClinginfo->HasDefaultConstructor(testio) != ROOT::TMetaUtils::EIOCtorCategory::kAbsent;
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
   return TClingClassInfo::IsEnum(GetInterpreterImpl(), name);
}

////////////////////////////////////////////////////////////////////////////////

Bool_t TCling::ClassInfo_IsScopedEnum(ClassInfo_t *info) const
{
   TClingClassInfo* TClinginfo = (TClingClassInfo*) info;
   return TClinginfo->IsScopedEnum();
}


////////////////////////////////////////////////////////////////////////////////

EDataType TCling::ClassInfo_GetUnderlyingType(ClassInfo_t* info) const
{
   TClingClassInfo* TClinginfo = (TClingClassInfo*) info;
   return TClinginfo->GetUnderlyingType();
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

bool TCling::ClassInfo_IsValidMethod(ClassInfo_t* cinfo, const char* method, const char* proto, Longptr_t* offset, EFunctionMatchMode mode /* = kConversionMatch */) const
{
   TClingClassInfo* TClinginfo = (TClingClassInfo*) cinfo;
   return TClinginfo->IsValidMethod(method, proto, false, offset, mode);
}

////////////////////////////////////////////////////////////////////////////////

bool TCling::ClassInfo_IsValidMethod(ClassInfo_t* cinfo, const char* method, const char* proto, Bool_t objectIsConst, Longptr_t* offset, EFunctionMatchMode mode /* = kConversionMatch */) const
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

Longptr_t TCling::ClassInfo_Tagnum(ClassInfo_t* cinfo) const
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
   return output.c_str(); // NOLINT
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
   return (BaseClassInfo_t*) new TClingBaseClassInfo(GetInterpreterImpl(), TClinginfo);
}

////////////////////////////////////////////////////////////////////////////////

BaseClassInfo_t* TCling::BaseClassInfo_Factory(ClassInfo_t* derived,
   ClassInfo_t* base) const
{
   R__LOCKGUARD(gInterpreterMutex);
   TClingClassInfo* TClinginfo = (TClingClassInfo*) derived;
   TClingClassInfo* TClinginfoBase = (TClingClassInfo*) base;
   return (BaseClassInfo_t*) new TClingBaseClassInfo(GetInterpreterImpl(), TClinginfo, TClinginfoBase);
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

Longptr_t TCling::BaseClassInfo_Offset(BaseClassInfo_t* toBaseClassInfo, void * address, bool isDerivedObject) const
{
   TClingBaseClassInfo* TClinginfo = (TClingBaseClassInfo*) toBaseClassInfo;
   return TClinginfo->Offset(address, isDerivedObject);
}

////////////////////////////////////////////////////////////////////////////////

Longptr_t TCling::ClassInfo_GetBaseOffset(ClassInfo_t* fromDerived, ClassInfo_t* toBase, void * address, bool isDerivedObject) const
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

Longptr_t TCling::BaseClassInfo_Tagnum(BaseClassInfo_t* bcinfo) const
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
   return output.c_str(); // NOLINT
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

DataMemberInfo_t* TCling::DataMemberInfo_Factory(ClassInfo_t* clinfo, TDictionary::EMemberSelection selection) const
{
   R__LOCKGUARD(gInterpreterMutex);
   TClingClassInfo* TClingclass_info = (TClingClassInfo*) clinfo;
   return (DataMemberInfo_t*) new TClingDataMemberInfo(GetInterpreterImpl(), TClingclass_info, selection);
}

////////////////////////////////////////////////////////////////////////////////

DataMemberInfo_t* TCling::DataMemberInfo_Factory(DeclId_t declid, ClassInfo_t* clinfo) const
{
   R__LOCKGUARD(gInterpreterMutex);
   const clang::Decl* decl = reinterpret_cast<const clang::Decl*>(declid);
   const clang::ValueDecl* vd = llvm::dyn_cast_or_null<clang::ValueDecl>(decl);
   return (DataMemberInfo_t*) new TClingDataMemberInfo(GetInterpreterImpl(), vd, (TClingClassInfo*)clinfo);
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

Longptr_t TCling::DataMemberInfo_Offset(DataMemberInfo_t* dminfo) const
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
   return result.c_str(); // NOLINT
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

static void ConstructorName(std::string &name, const clang::Decl *decl,
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

void TCling::GetFunctionName(const clang::Decl *decl, std::string &output) const
{
   output.clear();

   const auto *FD = llvm::dyn_cast<clang::FunctionDecl>(decl);
   if (const auto *USD = llvm::dyn_cast<clang::UsingShadowDecl>(decl)) {
      FD = llvm::dyn_cast<clang::FunctionDecl>(USD->getTargetDecl());
   }
   if (!FD) {
      Error("GetFunctionName", "NULL Decl!");
      return;
   }

   // For using-decls, show "Derived", not "Base", i.e. use the
   // name of the decl context of the UsingShadowDecl (aka `decl`)
   // not the name of FD's decl context.
   if (llvm::isa<clang::CXXConstructorDecl>(FD))
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
      FD->getNameForDiagnostic(stream, printPolicy, /*Qualified=*/false);
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
      if (md->getMethodQualifiers().hasConst()) {
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
/// Return the property not already defined in Property
/// See TDictionary's EFunctionProperty

Long_t TCling::FuncTempInfo_ExtraProperty(FuncTempInfo_t* ft_info) const
{
   if (!ft_info) return 0;

   long property = 0L;
   property |= kIsCompiled;

   const clang::FunctionTemplateDecl *ft = (clang::FunctionTemplateDecl*)ft_info;
   const clang::FunctionDecl *fd = ft->getTemplatedDecl();

   if (fd->isOverloadedOperator())
      property |= kIsOperator;
   if (llvm::isa<clang::CXXConversionDecl>(fd))
      property |= kIsConversion;
   if (llvm::isa<clang::CXXConstructorDecl>(fd))
      property |= kIsConstructor;
   if (llvm::isa<clang::CXXDestructorDecl>(fd))
      property |= kIsDestructor;
   if (fd->isInlined())
      property |= kIsInlined;
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
   return (MethodInfo_t*) new TClingMethodInfo(GetInterpreterImpl());
}

////////////////////////////////////////////////////////////////////////////////

MethodInfo_t* TCling::MethodInfo_Factory(ClassInfo_t* clinfo) const
{
   R__LOCKGUARD(gInterpreterMutex);
   return (MethodInfo_t*) new TClingMethodInfo(GetInterpreterImpl(), (TClingClassInfo*)clinfo);
}

////////////////////////////////////////////////////////////////////////////////

MethodInfo_t* TCling::MethodInfo_Factory(DeclId_t declid) const
{
   const clang::Decl* decl = reinterpret_cast<const clang::Decl*>(declid);
   R__LOCKGUARD(gInterpreterMutex);
   return (MethodInfo_t*) new TClingMethodInfo(GetInterpreterImpl(), decl);
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
   return info->GetPrototype();
}

////////////////////////////////////////////////////////////////////////////////

const char* TCling::MethodInfo_Name(MethodInfo_t* minfo) const
{
   TClingMethodInfo* info = (TClingMethodInfo*) minfo;
   return info->Name();
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
   return (MethodArgInfo_t*) new TClingMethodArgInfo(GetInterpreterImpl());
}

////////////////////////////////////////////////////////////////////////////////

MethodArgInfo_t* TCling::MethodArgInfo_Factory(MethodInfo_t *minfo) const
{
   R__LOCKGUARD(gInterpreterMutex);
   return (MethodArgInfo_t*) new TClingMethodArgInfo(GetInterpreterImpl(), (TClingMethodInfo*)minfo);
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
   return (TypeInfo_t*) new TClingTypeInfo(GetInterpreterImpl());
}

////////////////////////////////////////////////////////////////////////////////

TypeInfo_t* TCling::TypeInfo_Factory(const char *name) const
{
   R__LOCKGUARD(gInterpreterMutex);
   return (TypeInfo_t*) new TClingTypeInfo(GetInterpreterImpl(), name);
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
   return (TypedefInfo_t*) new TClingTypedefInfo(GetInterpreterImpl());
}

////////////////////////////////////////////////////////////////////////////////

TypedefInfo_t* TCling::TypedefInfo_Factory(const char *name) const
{
   R__LOCKGUARD(gInterpreterMutex);
   return (TypedefInfo_t*) new TClingTypedefInfo(GetInterpreterImpl(), name);
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
   if (!fInitialMutex) {
      if (fInitialMutex.fRecurseCount) {
         Error("SnapshotMutexState", "fRecurseCount != 0 even though initial mutex state is unset!");
      }
      fInitialMutex.fState = mtx->GetStateBefore();
   }
   // We will "forget" this lock once we backed out of all interpreter frames.
   // Here we are entering one, so ++.
   ++fInitialMutex.fRecurseCount;
}

////////////////////////////////////////////////////////////////////////////////

void TCling::ForgetMutexState()
{
   if (!fInitialMutex)
      return;
   if (fInitialMutex.fRecurseCount == 0) {
      Error("ForgetMutexState", "mutex state's recurse count already 0!");
   }
   else if (--fInitialMutex.fRecurseCount == 0) {
      // We have returned from all interpreter frames. Reset the initial lock state.
      fInitialMutex.fState.reset();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Re-apply the lock count delta that TCling__ResetInterpreterMutex() caused.

void TCling::ApplyToInterpreterMutex(void *delta)
{
   if (gInterpreterMutex) {
      if (delta) {
         auto typedDelta = static_cast<MutexStateAndRecurseCountDelta *>(delta);
         std::unique_ptr<MutexStateAndRecurseCountDelta> uniqueP{typedDelta};
         gCoreMutex->Apply(std::move(typedDelta->fDelta));
         // Now that we have the lock, update the global
         R__ASSERT(fInitialMutex.fRecurseCount == 0 && "Inconsistent state of fInitialMutex!  Another thread within Interpreter critical section.");
         std::swap(fInitialMutex, typedDelta->fInitialState);
      } else {
         // This case happens when EnableThreadSafety is first called from
         // the interpreter function we just handled.
         // Since thread safety was not enabled at the time we rewound, there was
         // no lock taken and even-though we should be locking the rest of this
         // interpreter handling/modifying code (since there might be threads in
         // flight), we can't because there would not be any lock guard to release the
         // locks
         if (fInitialMutex || fInitialMutex.fRecurseCount !=0)
            Error("ApplyToInterpreterMutex",
                 "After returning from user code that turned on thread safety support, we notice that fInitialMutex is already used ... "
                 "so the rest of this function/stack execution might have race condition (with the other thread that thinks it has exclusive access to the interpreter state.");
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Reset the interpreter lock to the state it had before interpreter-related
/// calls happened.

void *TCling::RewindInterpreterMutex()
{
   if (fInitialMutex) {
      // Need to start a new recurse count.
      std::unique_ptr<MutexStateAndRecurseCountDelta> uniqueP(new MutexStateAndRecurseCountDelta());
      std::swap(uniqueP->fInitialState, fInitialMutex);
      uniqueP->fDelta = gCoreMutex->Rewind(*uniqueP->fInitialState.fState);
      return uniqueP.release();
   }
   R__ASSERT(fInitialMutex.fRecurseCount == 0);
   return nullptr;
}
