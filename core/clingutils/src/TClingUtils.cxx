// @(#)root/metautils:$Id$
// Author: Paul Russo, 2009-10-06

/*************************************************************************
 * Copyright (C) 1995-2011, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//______________________________________________________________________________
//                                                                      //
// ROOT::TMetaUtils provides utility wrappers around                    //
// cling, the LLVM-based interpreter. It's an internal set of tools     //
// used by TCling and rootcling.                                        //
//                                                                      //
//______________________________________________________________________________
#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <unordered_set>

#include "RConfigure.h"
#include <ROOT/RConfig.hxx>
#include <ROOT/FoundationUtils.hxx>
#include "Rtypes.h"

#include "RStl.h"

#include "clang/AST/ASTContext.h"
#include "clang/AST/Attr.h"
#include "clang/AST/CXXInheritance.h"
#include "clang/AST/Decl.h"
#include "clang/AST/DeclTemplate.h"
#include "clang/AST/Mangle.h"
#include "clang/AST/Type.h"
#include "clang/AST/TypeVisitor.h"
#include "clang/Frontend/CompilerInstance.h"
#include "clang/Lex/HeaderSearch.h"
#include "clang/Lex/ModuleMap.h"
#include "clang/Lex/Preprocessor.h"
#include "clang/Lex/PreprocessorOptions.h"

#include "clang/Sema/Sema.h"
#include "clang/Sema/SemaDiagnostic.h"

#include "cling/Interpreter/LookupHelper.h"
#include "cling/Interpreter/Transaction.h"
#include "cling/Interpreter/Interpreter.h"
#include "cling/Utils/AST.h"

#include "llvm/Support/Path.h"
#include "llvm/Support/FileSystem.h"

// Intentionally access non-public header ...
#include "../../../interpreter/llvm/src/tools/clang/lib/Sema/HackForDefaultTemplateArg.h"

#include "TClingUtils.h"

#ifdef _WIN32
#define strncasecmp _strnicmp
#include <io.h>
#else
#include <unistd.h>
#endif // _WIN32

namespace ROOT {
namespace TMetaUtils {

std::string GetRealPath(const std::string &path)
{
   llvm::SmallString<256> result_path;
   llvm::sys::fs::real_path(path, result_path, /*expandTilde*/true);
   return result_path.str().str();
}


////////////////////////////////////////////////////////////////////////////////

class TNormalizedCtxtImpl {
   using DeclsCont_t = TNormalizedCtxt::Config_t::SkipCollection;
   using Config_t = TNormalizedCtxt::Config_t;
   using TypesCont_t = TNormalizedCtxt::TypesCont_t;
   using TemplPtrIntMap_t = TNormalizedCtxt::TemplPtrIntMap_t;
private:
   Config_t         fConfig;
   TypesCont_t      fTypeWithAlternative;
   static TemplPtrIntMap_t fTemplatePtrArgsToKeepMap;
public:
   TNormalizedCtxtImpl(const cling::LookupHelper &lh);

   const Config_t    &GetConfig() const { return fConfig; }
   const TypesCont_t &GetTypeWithAlternative() const { return fTypeWithAlternative; }
   void AddTemplAndNargsToKeep(const clang::ClassTemplateDecl* templ, unsigned int i);
   int GetNargsToKeep(const clang::ClassTemplateDecl* templ) const;
   const TemplPtrIntMap_t GetTemplNargsToKeepMap() const { return fTemplatePtrArgsToKeepMap; }
   void keepTypedef(const cling::LookupHelper &lh, const char* name,
                    bool replace = false);
};
}
}

namespace {

////////////////////////////////////////////////////////////////////////////////
/// Add default parameter to the scope if needed.

static clang::NestedNameSpecifier* AddDefaultParametersNNS(const clang::ASTContext& Ctx,
                                                           clang::NestedNameSpecifier* scope,
                                                           const cling::Interpreter &interpreter,
                                                           const ROOT::TMetaUtils::TNormalizedCtxt &normCtxt) {
   if (!scope) return 0;

   const clang::Type* scope_type = scope->getAsType();
   if (scope_type) {
      // this is not a namespace, so we might need to desugar
      clang::NestedNameSpecifier* outer_scope = scope->getPrefix();
      if (outer_scope) {
         outer_scope = AddDefaultParametersNNS(Ctx, outer_scope, interpreter, normCtxt);
      }

      clang::QualType addDefault =
         ROOT::TMetaUtils::AddDefaultParameters(clang::QualType(scope_type,0), interpreter, normCtxt );
      // NOTE: Should check whether the type has changed or not.
      if (addDefault.getTypePtr() != scope_type)
         return clang::NestedNameSpecifier::Create(Ctx,outer_scope,
                                                   false /* template keyword wanted */,
                                                   addDefault.getTypePtr());
   }
   return scope;
}

////////////////////////////////////////////////////////////////////////////////

static bool CheckDefinition(const clang::CXXRecordDecl *cl, const clang::CXXRecordDecl *context)
{
   if (!cl->hasDefinition()) {
      if (context) {
         ROOT::TMetaUtils::Error("CheckDefinition",
                                 "Missing definition for class %s, please #include its header in the header of %s\n",
                                 cl->getName().str().c_str(), context->getName().str().c_str());
      } else {
         ROOT::TMetaUtils::Error("CheckDefinition",
                                 "Missing definition for class %s\n",
                                 cl->getName().str().c_str());
      }
      return false;
   }
   return true;
}

////////////////////////////////////////////////////////////////////////////////
/// Check if 'scope' or any of its template parameter was substituted when
/// instantiating the class template instance and replace it with the
/// partially sugared types we have from 'instance'.

static clang::NestedNameSpecifier* ReSubstTemplateArgNNS(const clang::ASTContext &Ctxt,
                                                         clang::NestedNameSpecifier *scope,
                                                         const clang::Type *instance)
{
   if (!scope) return 0;

   const clang::Type* scope_type = scope->getAsType();
   if (scope_type) {
      clang::NestedNameSpecifier* outer_scope = scope->getPrefix();
      if (outer_scope) {
         outer_scope = ReSubstTemplateArgNNS(Ctxt, outer_scope, instance);
      }
      clang::QualType substScope =
         ROOT::TMetaUtils::ReSubstTemplateArg(clang::QualType(scope_type,0), instance);
      // NOTE: Should check whether the type has changed or not.
      scope = clang::NestedNameSpecifier::Create(Ctxt,outer_scope,
                                                 false /* template keyword wanted */,
                                                 substScope.getTypePtr());
   }
   return scope;
}

////////////////////////////////////////////////////////////////////////////////

static bool IsTypeInt(const clang::Type *type)
{
   const clang::BuiltinType * builtin = llvm::dyn_cast<clang::BuiltinType>(type->getCanonicalTypeInternal().getTypePtr());
   if (builtin) {
      return builtin->isInteger(); // builtin->getKind() == clang::BuiltinType::Int;
   } else {
      return false;
   }
}

////////////////////////////////////////////////////////////////////////////////

static bool IsFieldDeclInt(const clang::FieldDecl *field)
{
   return IsTypeInt(field->getType().getTypePtr());
}

////////////////////////////////////////////////////////////////////////////////
/// Return a data member name 'what' in the class described by 'cl' if any.

static const clang::FieldDecl *GetDataMemberFromAll(const clang::CXXRecordDecl &cl, llvm::StringRef what)
{
   for(clang::RecordDecl::field_iterator field_iter = cl.field_begin(), end = cl.field_end();
       field_iter != end;
       ++field_iter){
      if (field_iter->getName() == what) {
         return *field_iter;
      }
   }
   return 0;

}

////////////////////////////////////////////////////////////////////////////////

static bool CXXRecordDecl__FindOrdinaryMember(const clang::CXXBaseSpecifier *Specifier,
                                              clang::CXXBasePath &Path,
                                              const char *Name)
{
   clang::RecordDecl *BaseRecord = Specifier->getType()->getAs<clang::RecordType>()->getDecl();

   const clang::CXXRecordDecl *clxx = llvm::dyn_cast<clang::CXXRecordDecl>(BaseRecord);
   if (clxx == 0) return false;

   const clang::FieldDecl *found = GetDataMemberFromAll(*clxx,(const char*)Name);
   if (found) {
      // Humm, this is somewhat bad (well really bad), oh well.
      // Let's hope Paths never thinks it owns those (it should not as far as I can tell).
      clang::NamedDecl* NonConstFD = const_cast<clang::FieldDecl*>(found);
      clang::NamedDecl** BaseSpecFirstHack
      = reinterpret_cast<clang::NamedDecl**>(NonConstFD);
      Path.Decls = clang::DeclContextLookupResult(llvm::ArrayRef<clang::NamedDecl*>(BaseSpecFirstHack, 1));
      return true;
   }
   //
   // This is inspired from CXXInheritance.cpp:
   /*
    *      RecordDecl *BaseRecord =
    *        Specifier->getType()->castAs<RecordType>()->getDecl();
    *
    *  const unsigned IDNS = clang::Decl::IDNS_Ordinary | clang::Decl::IDNS_Tag | clang::Decl::IDNS_Member;
    *  clang::DeclarationName N = clang::DeclarationName::getFromOpaquePtr(Name);
    *  for (Path.Decls = BaseRecord->lookup(N);
    *       Path.Decls.first != Path.Decls.second;
    *       ++Path.Decls.first) {
    *     if ((*Path.Decls.first)->isInIdentifierNamespace(IDNS))
    *        return true;
    }
    */
   return false;

}

////////////////////////////////////////////////////////////////////////////////
/// Return a data member name 'what' in any of the base classes of the class described by 'cl' if any.

static const clang::FieldDecl *GetDataMemberFromAllParents(const clang::CXXRecordDecl &cl, const char *what)
{
   clang::CXXBasePaths Paths;
   Paths.setOrigin(const_cast<clang::CXXRecordDecl*>(&cl));
   if (cl.lookupInBases([=](const clang::CXXBaseSpecifier *Specifier, clang::CXXBasePath &Path) {
            return CXXRecordDecl__FindOrdinaryMember(Specifier, Path, what);}, Paths))
   {
      clang::CXXBasePaths::paths_iterator iter = Paths.begin();
      if (iter != Paths.end()) {
         // See CXXRecordDecl__FindOrdinaryMember, this is, well, awkward.
         const clang::FieldDecl *found = (clang::FieldDecl *)iter->Decls.data();
         return found;
      }
   }
   return 0;
}

static
cling::LookupHelper::DiagSetting ToLHDS(bool wantDiags) {
   return wantDiags
      ? cling::LookupHelper::WithDiagnostics
      : cling::LookupHelper::NoDiagnostics;
}

} // end of anonymous namespace


namespace ROOT {
namespace TMetaUtils {

////////////////////////////////////////////////////////////////////////////////
/// Add to the internal map the pointer of a template as key and the number of
/// template arguments to keep as value.

void TNormalizedCtxtImpl::AddTemplAndNargsToKeep(const clang::ClassTemplateDecl* templ,
                                             unsigned int i){
   if (!templ){
      Error("TNormalizedCtxt::AddTemplAndNargsToKeep",
            "Tring to specify a number of template arguments to keep for a null pointer. Exiting without assigning any value.\n");
      return;
   }

   const clang::ClassTemplateDecl* canTempl = templ->getCanonicalDecl();

   if(fTemplatePtrArgsToKeepMap.count(canTempl)==1 &&
      fTemplatePtrArgsToKeepMap[canTempl]!=(int)i){
      const std::string templateName (canTempl->getNameAsString());
      const std::string i_str (std::to_string(i));
      const std::string previousArgsToKeep(std::to_string(fTemplatePtrArgsToKeepMap[canTempl]));
      Error("TNormalizedCtxt::AddTemplAndNargsToKeep",
            "Tring to specify for template %s %s arguments to keep, while before this number was %s\n",
            canTempl->getNameAsString().c_str(),
            i_str.c_str(),
            previousArgsToKeep.c_str());
   }

   fTemplatePtrArgsToKeepMap[canTempl]=i;
}
////////////////////////////////////////////////////////////////////////////////
/// Get from the map the number of arguments to keep.
/// It uses the canonical decl of the template as key.
/// If not present, returns -1.

int TNormalizedCtxtImpl::GetNargsToKeep(const clang::ClassTemplateDecl* templ) const{
   const clang::ClassTemplateDecl* constTempl = templ->getCanonicalDecl();
   auto thePairPtr = fTemplatePtrArgsToKeepMap.find(constTempl);
   int nArgsToKeep = (thePairPtr != fTemplatePtrArgsToKeepMap.end() ) ? thePairPtr->second : -1;
   return nArgsToKeep;
}


////////////////////////////////////////////////////////////////////////////////

TNormalizedCtxt::TNormalizedCtxt(const cling::LookupHelper &lh):
   fImpl(new TNormalizedCtxtImpl(lh))
{}

TNormalizedCtxt::TNormalizedCtxt(const TNormalizedCtxt& other):
   fImpl(new TNormalizedCtxtImpl(*other.fImpl))
{}

TNormalizedCtxt::~TNormalizedCtxt() {
   delete fImpl;
}
const TNormalizedCtxt::Config_t &TNormalizedCtxt::GetConfig() const {
   return fImpl->GetConfig();
}
const TNormalizedCtxt::TypesCont_t &TNormalizedCtxt::GetTypeWithAlternative() const {
   return fImpl->GetTypeWithAlternative();
}
void TNormalizedCtxt::AddTemplAndNargsToKeep(const clang::ClassTemplateDecl* templ, unsigned int i)
{
   return fImpl->AddTemplAndNargsToKeep(templ, i);
}
int TNormalizedCtxt::GetNargsToKeep(const clang::ClassTemplateDecl* templ) const
{
   return fImpl->GetNargsToKeep(templ);
}
const TNormalizedCtxt::TemplPtrIntMap_t TNormalizedCtxt::GetTemplNargsToKeepMap() const {
   return fImpl->GetTemplNargsToKeepMap();
}
void TNormalizedCtxt::keepTypedef(const cling::LookupHelper &lh, const char* name,
                                  bool replace /*= false*/)
{
   return fImpl->keepTypedef(lh, name, replace);
}

std::string AnnotatedRecordDecl::BuildDemangledTypeInfo(const clang::RecordDecl *rDecl,
                                                        const std::string &normalizedName)
{
   // Types with strong typedefs must not be findable through demangled type names, or else
   // the demangled name will resolve to both sinblings double / Double32_t.
   if (normalizedName.find("Double32_t") != std::string::npos
       || normalizedName.find("Float16_t") != std::string::npos)
       return {};
   std::unique_ptr<clang::MangleContext> mangleCtx(rDecl->getASTContext().createMangleContext());
   std::string mangledName;
   {
      llvm::raw_string_ostream sstr(mangledName);
      if (const clang::TypeDecl* TD = llvm::dyn_cast<clang::TypeDecl>(rDecl)) {
         mangleCtx->mangleCXXRTTI(clang::QualType(TD->getTypeForDecl(), 0), sstr);
      }
   }
   if (!mangledName.empty()) {
      int errDemangle = 0;
#ifdef WIN32
      if (mangledName[0] == '\01')
         mangledName.erase(0, 1);
      char *demangledTIName = TClassEdit::DemangleName(mangledName.c_str(), errDemangle);
      if (!errDemangle && demangledTIName) {
         static const char typeinfoNameFor[] = " `RTTI Type Descriptor'";
         if (strstr(demangledTIName, typeinfoNameFor)) {
            std::string demangledName = demangledTIName;
            demangledName.erase(demangledName.end() - strlen(typeinfoNameFor), demangledName.end());
#else
      char* demangledTIName = TClassEdit::DemangleName(mangledName.c_str(), errDemangle);
      if (!errDemangle && demangledTIName) {
         static const char typeinfoNameFor[] = "typeinfo for ";
         if (!strncmp(demangledTIName, typeinfoNameFor, strlen(typeinfoNameFor))) {
            std::string demangledName = demangledTIName + strlen(typeinfoNameFor);
#endif
            free(demangledTIName);
            return demangledName;
         } else {
#ifdef WIN32
            ROOT::TMetaUtils::Error("AnnotatedRecordDecl::BuildDemangledTypeInfo",
                                    "Demangled typeinfo name '%s' does not contain `RTTI Type Descriptor'\n",
                                    demangledTIName);
#else
            ROOT::TMetaUtils::Error("AnnotatedRecordDecl::BuildDemangledTypeInfo",
                                    "Demangled typeinfo name '%s' does not start with 'typeinfo for'\n",
                                    demangledTIName);
#endif
         } // if demangled type_info starts with "typeinfo for "
      } // if demangling worked
      free(demangledTIName);
   } // if mangling worked
   return {};
}


////////////////////////////////////////////////////////////////////////////////
/// There is no requested type name.
/// Still let's normalized the actual name.

AnnotatedRecordDecl::AnnotatedRecordDecl(long index,
                                         const clang::RecordDecl *decl,
                                         bool rStreamerInfo,
                                         bool rNoStreamer,
                                         bool rRequestNoInputOperator,
                                         bool rRequestOnlyTClass,
                                         int rRequestedVersionNumber,
                                         const cling::Interpreter &interpreter,
                                         const TNormalizedCtxt &normCtxt) :
   fRuleIndex(index), fDecl(decl), fRequestStreamerInfo(rStreamerInfo), fRequestNoStreamer(rNoStreamer),
   fRequestNoInputOperator(rRequestNoInputOperator), fRequestOnlyTClass(rRequestOnlyTClass), fRequestedVersionNumber(rRequestedVersionNumber)
{
   TMetaUtils::GetNormalizedName(fNormalizedName, decl->getASTContext().getTypeDeclType(decl), interpreter,normCtxt);
   fDemangledTypeInfo = BuildDemangledTypeInfo(decl, fNormalizedName);
}

////////////////////////////////////////////////////////////////////////////////
/// Normalize the requested type name.

AnnotatedRecordDecl::AnnotatedRecordDecl(long index,
                                         const clang::Type *requestedType,
                                         const clang::RecordDecl *decl,
                                         const char *requestName,
                                         unsigned int nTemplateArgsToSkip,
                                         bool rStreamerInfo,
                                         bool rNoStreamer,
                                         bool rRequestNoInputOperator,
                                         bool rRequestOnlyTClass,
                                         int rRequestVersionNumber,
                                         const cling::Interpreter &interpreter,
                                         const TNormalizedCtxt &normCtxt) :
   fRuleIndex(index), fDecl(decl), fRequestedName(""), fRequestStreamerInfo(rStreamerInfo), fRequestNoStreamer(rNoStreamer),
   fRequestNoInputOperator(rRequestNoInputOperator), fRequestOnlyTClass(rRequestOnlyTClass), fRequestedVersionNumber(rRequestVersionNumber)
{
   // For comparison purposes.
   TClassEdit::TSplitType splitname1(requestName,(TClassEdit::EModType)(TClassEdit::kLong64 | TClassEdit::kDropStd));
   splitname1.ShortType(fRequestedName, 0);

   TMetaUtils::GetNormalizedName( fNormalizedName, clang::QualType(requestedType,0), interpreter, normCtxt);
   if ( 0!=TMetaUtils::RemoveTemplateArgsFromName( fNormalizedName, nTemplateArgsToSkip) ){
      ROOT::TMetaUtils::Warning("AnnotatedRecordDecl",
                                "Could not remove the requested template arguments.\n");
   }
   fDemangledTypeInfo = BuildDemangledTypeInfo(decl, fNormalizedName);
}

////////////////////////////////////////////////////////////////////////////////
/// Normalize the requested type name.

AnnotatedRecordDecl::AnnotatedRecordDecl(long index,
                                         const clang::Type *requestedType,
                                         const clang::RecordDecl *decl,
                                         const char *requestName,
                                         bool rStreamerInfo,
                                         bool rNoStreamer,
                                         bool rRequestNoInputOperator,
                                         bool rRequestOnlyTClass,
                                         int rRequestVersionNumber,
                                         const cling::Interpreter &interpreter,
                                         const TNormalizedCtxt &normCtxt) :
   fRuleIndex(index), fDecl(decl), fRequestedName(""), fRequestStreamerInfo(rStreamerInfo), fRequestNoStreamer(rNoStreamer),
   fRequestNoInputOperator(rRequestNoInputOperator), fRequestOnlyTClass(rRequestOnlyTClass), fRequestedVersionNumber(rRequestVersionNumber)
{
   // For comparison purposes.
   TClassEdit::TSplitType splitname1(requestName,(TClassEdit::EModType)(TClassEdit::kLong64 | TClassEdit::kDropStd));
   splitname1.ShortType(fRequestedName, 0);

   TMetaUtils::GetNormalizedName( fNormalizedName, clang::QualType(requestedType,0), interpreter, normCtxt);
   fDemangledTypeInfo = BuildDemangledTypeInfo(decl, fNormalizedName);
}

////////////////////////////////////////////////////////////////////////////////
/// Normalize the requested name.

AnnotatedRecordDecl::AnnotatedRecordDecl(long index,
                                         const clang::RecordDecl *decl,
                                         const char *requestName,
                                         bool rStreamerInfo,
                                         bool rNoStreamer,
                                         bool rRequestNoInputOperator,
                                         bool rRequestOnlyTClass,
                                         int rRequestVersionNumber,
                                         const cling::Interpreter &interpreter,
                                         const TNormalizedCtxt &normCtxt) :
   fRuleIndex(index), fDecl(decl), fRequestedName(""), fRequestStreamerInfo(rStreamerInfo), fRequestNoStreamer(rNoStreamer), fRequestNoInputOperator(rRequestNoInputOperator), fRequestOnlyTClass(rRequestOnlyTClass), fRequestedVersionNumber(rRequestVersionNumber)
{
   // const clang::ClassTemplateSpecializationDecl *tmplt_specialization = llvm::dyn_cast<clang::ClassTemplateSpecializationDecl> (decl);
   // if (tmplt_specialization) {
   //    tmplt_specialization->getTemplateArgs ().data()->print(decl->getASTContext().getPrintingPolicy(),llvm::outs());
   //    llvm::outs() << "\n";
   // }
   // const char *current = requestName;
   // Strips spaces and std::
   if (requestName && requestName[0]) {
      TClassEdit::TSplitType splitname(requestName,(TClassEdit::EModType)( TClassEdit::kLong64 | TClassEdit::kDropStd));
      splitname.ShortType( fRequestedName, TClassEdit::kLong64 | TClassEdit::kDropStd );

      fNormalizedName = fRequestedName;
   } else {
      TMetaUtils::GetNormalizedName( fNormalizedName, decl->getASTContext().getTypeDeclType(decl),interpreter,normCtxt);
   }
   fDemangledTypeInfo = BuildDemangledTypeInfo(decl, fNormalizedName);
}

////////////////////////////////////////////////////////////////////////////////

TClingLookupHelper::TClingLookupHelper(cling::Interpreter &interpreter,
                                       TNormalizedCtxt &normCtxt,
                                       ExistingTypeCheck_t existingTypeCheck,
                                       AutoParse_t autoParse,
                                       bool *shuttingDownPtr,
                                       const int* pgDebug /*= 0*/):
   fInterpreter(&interpreter),fNormalizedCtxt(&normCtxt),
   fExistingTypeCheck(existingTypeCheck),
   fAutoParse(autoParse),
   fInterpreterIsShuttingDownPtr(shuttingDownPtr),
   fPDebug(pgDebug)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Helper routine to ry hard to avoid looking up in the Cling database as
/// this could enduce an unwanted autoparsing.

bool TClingLookupHelper::ExistingTypeCheck(const std::string &tname,
                                           std::string &result)
{
   if (tname.empty()) return false;

   if (fExistingTypeCheck) return fExistingTypeCheck(tname,result);
   else return false;
}

////////////////////////////////////////////////////////////////////////////////

void TClingLookupHelper::GetPartiallyDesugaredName(std::string &nameLong)
{
   const cling::LookupHelper& lh = fInterpreter->getLookupHelper();
   clang::QualType t = lh.findType(nameLong, ToLHDS(WantDiags()));
   if (!t.isNull()) {
      clang::QualType dest = cling::utils::Transform::GetPartiallyDesugaredType(fInterpreter->getCI()->getASTContext(), t, fNormalizedCtxt->GetConfig(), true /* fully qualify */);
      if (!dest.isNull() && (dest != t)) {
         // getAsStringInternal() appends.
         nameLong.clear();
         dest.getAsStringInternal(nameLong, fInterpreter->getCI()->getASTContext().getPrintingPolicy());
      }
   }
}

////////////////////////////////////////////////////////////////////////////////

bool TClingLookupHelper::IsAlreadyPartiallyDesugaredName(const std::string &nondef,
                                                         const std::string &nameLong)
{
   const cling::LookupHelper& lh = fInterpreter->getLookupHelper();
   clang::QualType t = lh.findType(nondef.c_str(), ToLHDS(WantDiags()));
   if (!t.isNull()) {
      clang::QualType dest = cling::utils::Transform::GetPartiallyDesugaredType(fInterpreter->getCI()->getASTContext(), t, fNormalizedCtxt->GetConfig(), true /* fully qualify */);
      if (!dest.isNull() && (dest != t) &&
          nameLong == t.getAsString(fInterpreter->getCI()->getASTContext().getPrintingPolicy()))
         return true;
   }
   return false;
}

////////////////////////////////////////////////////////////////////////////////

bool TClingLookupHelper::IsDeclaredScope(const std::string &base, bool &isInlined)
{
   const cling::LookupHelper& lh = fInterpreter->getLookupHelper();
   const clang::Decl *scope = lh.findScope(base.c_str(), ToLHDS(WantDiags()), 0);

   if (!scope) {
      // the nesting namespace is not declared
      isInlined = false;
      return false;
   }
   const clang::NamespaceDecl *nsdecl = llvm::dyn_cast<clang::NamespaceDecl>(scope);
   isInlined = nsdecl && nsdecl->isInline();
   return true;
}

////////////////////////////////////////////////////////////////////////////////
/// We assume that we have a simple type:
/// [const] typename[*&][const]

bool TClingLookupHelper::GetPartiallyDesugaredNameWithScopeHandling(const std::string &tname,
                                                                    std::string &result,
                                                                    bool dropstd /* = true */)
{
   if (tname.empty()) return false;

   // Try hard to avoid looking up in the Cling database as this could enduce
   // an unwanted autoparsing.
   // Note: this is always done by the callers and thus is redundant.
   // Maybe replace with
   assert(! (fExistingTypeCheck && fExistingTypeCheck(tname,result)) );
   if (fExistingTypeCheck && fExistingTypeCheck(tname,result)) {
      return ! result.empty();
   }

   if (fAutoParse) fAutoParse(tname.c_str());

   // Since we already check via other means (TClassTable which is populated by
   // the dictonary loading, and the gROOT list of classes and enums, which are
   // populated via TProtoClass/Enum), we should be able to disable the autoloading
   // ... which requires access to libCore or libCling ...
   const cling::LookupHelper& lh = fInterpreter->getLookupHelper();
   clang::QualType t = lh.findType(tname.c_str(), ToLHDS(WantDiags()));
   // Technically we ought to try:
   //   if (t.isNull()) t =  lh.findType(TClassEdit::InsertStd(tname), ToLHDS(WantDiags()));
   // at least until the 'normalized name' contains the std:: prefix.

   if (!t.isNull()) {
      clang::QualType dest = TMetaUtils::GetNormalizedType(t, *fInterpreter, *fNormalizedCtxt);
      if (!dest.isNull() && dest != t) {
         // Since our input is not a template instance name, rather than going through the full
         // TMetaUtils::GetNormalizedName, we just do the 'strip leading std' and fix
         // white space.
         clang::PrintingPolicy policy(fInterpreter->getCI()->getASTContext().getPrintingPolicy());
         policy.SuppressTagKeyword = true; // Never get the class or struct keyword
         policy.SuppressScope = true;      // Force the scope to be coming from a clang::ElaboratedType.
         // The scope suppression is required for getting rid of the anonymous part of the name of a class defined in an anonymous namespace.
         // This gives us more control vs not using the clang::ElaboratedType and relying on the Policy.SuppressUnwrittenScope which would
         // strip both the anonymous and the inline namespace names (and we probably do not want the later to be suppressed).
         // getAsStringInternal() appends.
         result.clear();
         dest.getAsStringInternal(result, policy);
         // Strip the std::
         unsigned long offset = 0;
         if (strncmp(result.c_str(), "const ", 6) == 0) {
            offset = 6;
         }
         if (dropstd && strncmp(result.c_str()+offset, "std::", 5) == 0) {
            result.erase(offset,5);
         }
         for(unsigned int i = 1; i<result.length(); ++i) {
            if (result[i]=='s') {
               if (result[i-1]=='<' || result[i-1]==',' || result[i-1]==' ') {
                  if (dropstd && result.compare(i,5,"std::",5) == 0) {
                     result.erase(i,5);
                  }
               }
            }
            if (result[i]==' ') {
               if (result[i-1] == ',') {
                  result.erase(i,1);
                  --i;
               } else if ( (i+1) < result.length() &&
                          (result[i+1]=='*' || result[i+1]=='&' || result[i+1]=='[') ) {
                  result.erase(i,1);
                  --i;
               }
            }
         }

//         std::string alt;
//         TMetaUtils::GetNormalizedName(alt, dest, *fInterpreter, *fNormalizedCtxt);
//         if (alt != result) fprintf(stderr,"norm: %s vs result=%s\n",alt.c_str(),result.c_str());

         return true;
      }
   }
   return false;
}

////////////////////////////////////////////////////////////////////////////////
// TClassEdit will call this routine as soon as any of its static variable (used
// for caching) is destroyed.
void TClingLookupHelper::ShuttingDownSignal()
{
   if (fInterpreterIsShuttingDownPtr)
      *fInterpreterIsShuttingDownPtr = true;
}

   } // end namespace ROOT
} // end namespace TMetaUtils


////////////////////////////////////////////////////////////////////////////////
/// Insert the type with name into the collection of typedefs to keep.
/// if replace, replace occurrences of the canonical type by name.

void ROOT::TMetaUtils::TNormalizedCtxtImpl::keepTypedef(const cling::LookupHelper &lh,
                                                    const char* name,
                                                    bool replace /*=false*/) {
   clang::QualType toSkip = lh.findType(name, cling::LookupHelper::WithDiagnostics);
   if (const clang::Type* T = toSkip.getTypePtr()) {
      const clang::TypedefType *tt = llvm::dyn_cast<clang::TypedefType>(T);
      if (!tt) return;
      clang::Decl* D = tt->getDecl();
      fConfig.m_toSkip.insert(D);
      if (replace) {
         clang::QualType canon = toSkip->getCanonicalTypeInternal();
         fConfig.m_toReplace.insert(std::make_pair(canon.getTypePtr(),T));
      } else {
         fTypeWithAlternative.insert(T);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize the list of typedef to keep (i.e. make them opaque for normalization)
/// and the list of typedef whose semantic is different from their underlying type
/// (Double32_t and Float16_t).
/// This might be specific to an interpreter.

ROOT::TMetaUtils::TNormalizedCtxtImpl::TNormalizedCtxtImpl(const cling::LookupHelper &lh)
{
   keepTypedef(lh, "Double32_t");
   keepTypedef(lh, "Float16_t");
   keepTypedef(lh, "Long64_t", true);
   keepTypedef(lh, "ULong64_t", true);

   clang::QualType toSkip = lh.findType("string", cling::LookupHelper::WithDiagnostics);
   if (!toSkip.isNull()) {
      if (const clang::TypedefType* TT
          = llvm::dyn_cast_or_null<clang::TypedefType>(toSkip.getTypePtr()))
         fConfig.m_toSkip.insert(TT->getDecl());
   }
   toSkip = lh.findType("std::string", cling::LookupHelper::WithDiagnostics);
   if (!toSkip.isNull()) {
      if (const clang::TypedefType* TT
          = llvm::dyn_cast_or_null<clang::TypedefType>(toSkip.getTypePtr()))
         fConfig.m_toSkip.insert(TT->getDecl());

      clang::QualType canon = toSkip->getCanonicalTypeInternal();
      fConfig.m_toReplace.insert(std::make_pair(canon.getTypePtr(),toSkip.getTypePtr()));
   }
}

using TNCtxtFullQual = ROOT::TMetaUtils::TNormalizedCtxtImpl;
TNCtxtFullQual::TemplPtrIntMap_t TNCtxtFullQual::fTemplatePtrArgsToKeepMap=TNCtxtFullQual::TemplPtrIntMap_t{};

////////////////////////////////////////////////////////////////////////////////

inline bool IsTemplate(const clang::Decl &cl)
{
   return (cl.getKind() == clang::Decl::ClassTemplatePartialSpecialization
           || cl.getKind() == clang::Decl::ClassTemplateSpecialization);
}


////////////////////////////////////////////////////////////////////////////////

const clang::FunctionDecl* ROOT::TMetaUtils::ClassInfo__HasMethod(const clang::DeclContext *cl, const char* name,
                                                            const cling::Interpreter& interp)
{
   clang::Sema* S = &interp.getSema();
   const clang::NamedDecl* ND = cling::utils::Lookup::Named(S, name, cl);
   if (ND == (clang::NamedDecl*)-1)
      return (clang::FunctionDecl*)-1;
   return llvm::dyn_cast_or_null<clang::FunctionDecl>(ND);
}

////////////////////////////////////////////////////////////////////////////////
/// Return the scope corresponding to 'name' or std::'name'

const clang::CXXRecordDecl *
ROOT::TMetaUtils::ScopeSearch(const char *name, const cling::Interpreter &interp,
                              bool /*diagnose*/, const clang::Type** resultType)
{
   const cling::LookupHelper& lh = interp.getLookupHelper();
   // We have many bogus diagnostics if we allow diagnostics here. Suppress.
   // FIXME: silence them in the callers.
   const clang::CXXRecordDecl *result
      = llvm::dyn_cast_or_null<clang::CXXRecordDecl>
      (lh.findScope(name, cling::LookupHelper::NoDiagnostics, resultType));
   if (!result) {
      std::string std_name("std::");
      std_name += name;
      // We have many bogus diagnostics if we allow diagnostics here. Suppress.
      // FIXME: silence them in the callers.
      result = llvm::dyn_cast_or_null<clang::CXXRecordDecl>
         (lh.findScope(std_name, cling::LookupHelper::NoDiagnostics, resultType));
   }
   return result;
}


////////////////////////////////////////////////////////////////////////////////

bool ROOT::TMetaUtils::RequireCompleteType(const cling::Interpreter &interp, const clang::CXXRecordDecl *cl)
{
   clang::QualType qType(cl->getTypeForDecl(),0);
   return RequireCompleteType(interp,cl->getLocation(),qType);
}

////////////////////////////////////////////////////////////////////////////////

bool ROOT::TMetaUtils::RequireCompleteType(const cling::Interpreter &interp, clang::SourceLocation Loc, clang::QualType Type)
{
   clang::Sema& S = interp.getCI()->getSema();
   // Here we might not have an active transaction to handle
   // the caused instantiation decl.
   cling::Interpreter::PushTransactionRAII RAII(const_cast<cling::Interpreter*>(&interp));
   return S.RequireCompleteType(Loc, Type, clang::diag::err_incomplete_type);
}

////////////////////////////////////////////////////////////////////////////////

bool ROOT::TMetaUtils::IsBase(const clang::CXXRecordDecl *cl, const clang::CXXRecordDecl *base,
                              const clang::CXXRecordDecl *context, const cling::Interpreter &interp)
{
   if (!cl || !base) {
      return false;
   }

   if (!cl->getDefinition() || !cl->isCompleteDefinition()) {
      RequireCompleteType(interp,cl);
   }

   if (!CheckDefinition(cl, context) || !CheckDefinition(base, context)) {
      return false;
   }

   if (!base->hasDefinition()) {
      ROOT::TMetaUtils::Error("IsBase", "Missing definition for class %s\n", base->getName().str().c_str());
      return false;
   }
   return cl->isDerivedFrom(base);
}

////////////////////////////////////////////////////////////////////////////////

bool ROOT::TMetaUtils::IsBase(const clang::FieldDecl &m, const char* basename, const cling::Interpreter &interp)
{
   const clang::CXXRecordDecl* CRD = llvm::dyn_cast<clang::CXXRecordDecl>(ROOT::TMetaUtils::GetUnderlyingRecordDecl(m.getType()));
   if (!CRD) {
      return false;
   }

   const clang::NamedDecl *base
      = ScopeSearch(basename, interp, true /*diagnose*/, 0);

   if (base) {
      return IsBase(CRD, llvm::dyn_cast<clang::CXXRecordDecl>( base ),
                    llvm::dyn_cast<clang::CXXRecordDecl>(m.getDeclContext()),interp);
   }
   return false;
}

////////////////////////////////////////////////////////////////////////////////

int ROOT::TMetaUtils::ElementStreamer(std::ostream& finalString,
                                      const clang::NamedDecl &forcontext,
                                      const clang::QualType &qti,
                                      const char *R__t,int rwmode,
                                      const cling::Interpreter &interp,
                                      const char *tcl)
{
   static const clang::CXXRecordDecl *TObject_decl
      = ROOT::TMetaUtils::ScopeSearch("TObject", interp, true /*diag*/, 0);
   enum {
      kBIT_ISTOBJECT     = 0x10000000,
      kBIT_HASSTREAMER   = 0x20000000,
      kBIT_ISSTRING      = 0x40000000,

      kBIT_ISPOINTER     = 0x00001000,
      kBIT_ISFUNDAMENTAL = 0x00000020,
      kBIT_ISENUM        = 0x00000008
   };

   const clang::Type &ti( * qti.getTypePtr() );
   std::string tiName;
   ROOT::TMetaUtils::GetQualifiedName(tiName, clang::QualType(&ti,0), forcontext);

   std::string objType(ROOT::TMetaUtils::ShortTypeName(tiName.c_str()));

   const clang::Type *rawtype = ROOT::TMetaUtils::GetUnderlyingType(clang::QualType(&ti,0));
   std::string rawname;
   ROOT::TMetaUtils::GetQualifiedName(rawname, clang::QualType(rawtype,0), forcontext);

   clang::CXXRecordDecl *cxxtype = rawtype->getAsCXXRecordDecl() ;
   int isStre = cxxtype && ROOT::TMetaUtils::ClassInfo__HasMethod(cxxtype,"Streamer",interp);
   int isTObj = cxxtype && (IsBase(cxxtype,TObject_decl,nullptr,interp) || rawname == "TObject");

   long kase = 0;

   if (ti.isPointerType())           kase |= kBIT_ISPOINTER;
   if (rawtype->isFundamentalType()) kase |= kBIT_ISFUNDAMENTAL;
   if (rawtype->isEnumeralType())    kase |= kBIT_ISENUM;


   if (isTObj)              kase |= kBIT_ISTOBJECT;
   if (isStre)              kase |= kBIT_HASSTREAMER;
   if (tiName == "string")  kase |= kBIT_ISSTRING;
   if (tiName == "string*") kase |= kBIT_ISSTRING;


   if (tcl == 0) {
      tcl = " internal error in rootcling ";
   }
   //    if (strcmp(objType,"string")==0) RStl::Instance().GenerateTClassFor( "string", interp, normCtxt  );

   if (rwmode == 0) {  //Read mode

      if (R__t) finalString << "            " << tiName << " " << R__t << ";" << std::endl;
      switch (kase) {

      case kBIT_ISFUNDAMENTAL:
         if (!R__t)  return 0;
         finalString << "            R__b >> " << R__t << ";" << std::endl;
         break;

      case kBIT_ISPOINTER|kBIT_ISTOBJECT|kBIT_HASSTREAMER:
         if (!R__t)  return 1;
         finalString << "            " << R__t << " = (" << tiName << ")R__b.ReadObjectAny(" << tcl << ");"  << std::endl;
         break;

      case kBIT_ISENUM:
         if (!R__t)  return 0;
         //             fprintf(fp, "            R__b >> (Int_t&)%s;\n",R__t);
         // On some platforms enums and not 'Int_t' and casting to a reference to Int_t
         // induces the silent creation of a temporary which is 'filled' __instead of__
         // the desired enum.  So we need to take it one step at a time.
         finalString << "            Int_t readtemp;" << std::endl
                       << "            R__b >> readtemp;" << std::endl
                       << "            " << R__t << " = static_cast<" << tiName << ">(readtemp);" << std::endl;
         break;

      case kBIT_HASSTREAMER:
      case kBIT_HASSTREAMER|kBIT_ISTOBJECT:
         if (!R__t)  return 0;
         finalString << "            " << R__t << ".Streamer(R__b);" << std::endl;
         break;

      case kBIT_HASSTREAMER|kBIT_ISPOINTER:
         if (!R__t)  return 1;
         //fprintf(fp, "            fprintf(stderr,\"info is %%p %%d\\n\",R__b.GetInfo(),R__b.GetInfo()?R__b.GetInfo()->GetOldVersion():-1);\n");
         finalString << "            if (R__b.GetInfo() && R__b.GetInfo()->GetOldVersion()<=3) {" << std::endl;
         if (cxxtype && cxxtype->isAbstract()) {
            finalString << "               R__ASSERT(0);// " << objType << " is abstract. We assume that older file could not be produced using this streaming method." << std::endl;
         } else {
            finalString << "               " << R__t << " = new " << objType << ";" << std::endl
                          << "               " << R__t << "->Streamer(R__b);" << std::endl;
         }
         finalString << "            } else {" << std::endl
                       << "               " << R__t << " = (" << tiName << ")R__b.ReadObjectAny(" << tcl << ");" << std::endl
                       << "            }" << std::endl;
         break;

      case kBIT_ISSTRING:
         if (!R__t)  return 0;
         finalString << "            {TString R__str;" << std::endl
                       << "             R__str.Streamer(R__b);" << std::endl
                       << "             " << R__t << " = R__str.Data();}" << std::endl;
         break;

      case kBIT_ISSTRING|kBIT_ISPOINTER:
         if (!R__t)  return 0;
         finalString << "            {TString R__str;"  << std::endl
                       << "             R__str.Streamer(R__b);" << std::endl
                       << "             " << R__t << " = new string(R__str.Data());}" << std::endl;
         break;

      case kBIT_ISPOINTER:
         if (!R__t)  return 1;
         finalString << "            " << R__t << " = (" << tiName << ")R__b.ReadObjectAny(" << tcl << ");" << std::endl;
         break;

      default:
         if (!R__t) return 1;
         finalString << "            R__b.StreamObject(&" << R__t << "," << tcl << ");" << std::endl;
         break;
      }

   } else {     //Write case

      switch (kase) {

      case kBIT_ISFUNDAMENTAL:
      case kBIT_ISPOINTER|kBIT_ISTOBJECT|kBIT_HASSTREAMER:
         if (!R__t)  return 0;
         finalString << "            R__b << " << R__t << ";" << std::endl;
         break;

      case kBIT_ISENUM:
         if (!R__t)  return 0;
         finalString << "            {  void *ptr_enum = (void*)&" << R__t << ";\n";
         finalString << "               R__b >> *reinterpret_cast<Int_t*>(ptr_enum); }" << std::endl;
         break;

      case kBIT_HASSTREAMER:
      case kBIT_HASSTREAMER|kBIT_ISTOBJECT:
         if (!R__t)  return 0;
         finalString << "            ((" << objType << "&)" << R__t << ").Streamer(R__b);" << std::endl;
         break;

      case kBIT_HASSTREAMER|kBIT_ISPOINTER:
         if (!R__t)  return 1;
         finalString << "            R__b.WriteObjectAny(" << R__t << "," << tcl << ");" << std::endl;
         break;

      case kBIT_ISSTRING:
         if (!R__t)  return 0;
         finalString << "            {TString R__str(" << R__t << ".c_str());" << std::endl
                       << "             R__str.Streamer(R__b);};" << std::endl;
         break;

      case kBIT_ISSTRING|kBIT_ISPOINTER:
         if (!R__t)  return 0;
         finalString << "            {TString R__str(" << R__t << "->c_str());" << std::endl
                       << "             R__str.Streamer(R__b);}" << std::endl;
         break;

      case kBIT_ISPOINTER:
         if (!R__t)  return 1;
         finalString << "            R__b.WriteObjectAny(" << R__t << "," << tcl <<");" << std::endl;
         break;

      default:
         if (!R__t)  return 1;
         finalString << "            R__b.StreamObject((" << objType << "*)&" << R__t << "," << tcl << ");" << std::endl;
         break;
      }
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Checks if default constructor exists and accessible

bool ROOT::TMetaUtils::CheckDefaultConstructor(const clang::CXXRecordDecl* cl, const cling::Interpreter& interpreter)
{
   clang::CXXRecordDecl* ncCl = const_cast<clang::CXXRecordDecl*>(cl);

   // We may induce template instantiation
   cling::Interpreter::PushTransactionRAII clingRAII(const_cast<cling::Interpreter*>(&interpreter));

   if (auto* Ctor = interpreter.getCI()->getSema().LookupDefaultConstructor(ncCl)) {
      if (Ctor->getAccess() == clang::AS_public && !Ctor->isDeleted()) {
         return true;
      }
   }

   return false;
}


////////////////////////////////////////////////////////////////////////////////
/// Checks IO constructor - must be public and with specified argument

ROOT::TMetaUtils::EIOCtorCategory ROOT::TMetaUtils::CheckIOConstructor(const clang::CXXRecordDecl *cl,
                                                                       const char *typeOfArg,
                                                                       const clang::CXXRecordDecl *expectedArgType,
                                                                       const cling::Interpreter& interpreter)
{
   if (typeOfArg && !expectedArgType) {
      const cling::LookupHelper& lh = interpreter.getLookupHelper();
      // We can not use findScope since the type we are given are usually,
      // only forward declared (and findScope explicitly reject them).
      clang::QualType instanceType = lh.findType(typeOfArg, cling::LookupHelper::WithDiagnostics);
      if (!instanceType.isNull())
         expectedArgType = instanceType->getAsCXXRecordDecl();
   }

   if (!expectedArgType)
      return EIOCtorCategory::kAbsent;

   // FIXME: We should not iterate here. That costs memory!
   cling::Interpreter::PushTransactionRAII clingRAII(const_cast<cling::Interpreter*>(&interpreter));
   for (auto iter = cl->ctor_begin(), end = cl->ctor_end(); iter != end; ++iter)
      {
         if ((iter->getAccess() != clang::AS_public) || (iter->getNumParams() != 1))
            continue;

         // We can reach this constructor.
         clang::QualType argType((*iter->param_begin())->getType());
         argType = argType.getDesugaredType(cl->getASTContext());
         // Deal with pointers and references: ROOT-7723
         auto ioCtorCategory = EIOCtorCategory::kAbsent;
         if (argType->isPointerType()) {
            ioCtorCategory = EIOCtorCategory::kIOPtrType;
            argType = argType->getPointeeType();
         } else if (argType->isReferenceType()) {
            ioCtorCategory = EIOCtorCategory::kIORefType;
            argType = argType.getNonReferenceType();
         } else
            continue;

         argType = argType.getDesugaredType(cl->getASTContext());
         const clang::CXXRecordDecl *argDecl = argType->getAsCXXRecordDecl();
         if (argDecl) {
            if (argDecl->getCanonicalDecl() == expectedArgType->getCanonicalDecl()) {
               return ioCtorCategory;
            }
         } else {
            std::string realArg = argType.getAsString();
            std::string clarg("class ");
            clarg += typeOfArg;
            if (realArg == clarg)
               return ioCtorCategory;
         }
   } // for each constructor

   return EIOCtorCategory::kAbsent;
}


////////////////////////////////////////////////////////////////////////////////
/// Check if class has constructor of provided type - either default or with single argument

ROOT::TMetaUtils::EIOCtorCategory ROOT::TMetaUtils::CheckConstructor(const clang::CXXRecordDecl *cl,
                                                                     const RConstructorType &ioctortype,
                                                                     const cling::Interpreter& interpreter)
{
   const char *arg = ioctortype.GetName();

   if (!ioctortype.GetType() && (!arg || !arg[0])) {
      // We are looking for a constructor with zero non-default arguments.

      return CheckDefaultConstructor(cl, interpreter) ? EIOCtorCategory::kDefault : EIOCtorCategory::kAbsent;
   }

   return CheckIOConstructor(cl, arg, ioctortype.GetType(), interpreter);
}


////////////////////////////////////////////////////////////////////////////////

const clang::CXXMethodDecl *GetMethodWithProto(const clang::Decl* cinfo,
                                               const char *method, const char *proto,
                                               const cling::Interpreter &interp,
                                               bool diagnose)
{
   const clang::FunctionDecl* funcD
      = interp.getLookupHelper().findFunctionProto(cinfo, method, proto,
                                                   diagnose ? cling::LookupHelper::WithDiagnostics
                                                   : cling::LookupHelper::NoDiagnostics);
   if (funcD)
      return llvm::dyn_cast<const clang::CXXMethodDecl>(funcD);

   return nullptr;
}


////////////////////////////////////////////////////////////////////////////////

namespace ROOT {
   namespace TMetaUtils {
      RConstructorType::RConstructorType(const char *type_of_arg, const cling::Interpreter &interp) : fArgTypeName(type_of_arg),fArgType(nullptr)
      {
         const cling::LookupHelper& lh = interp.getLookupHelper();
         // We can not use findScope since the type we are given are usually,
         // only forward declared (and findScope explicitly reject them).
         clang::QualType instanceType = lh.findType(type_of_arg, cling::LookupHelper::WithDiagnostics);
         if (!instanceType.isNull())
            fArgType = instanceType->getAsCXXRecordDecl();
      }
      const char *RConstructorType::GetName() const { return fArgTypeName.c_str(); }
      const clang::CXXRecordDecl *RConstructorType::GetType() const { return fArgType; }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// return true if we can find an constructor calleable without any arguments
/// or with one the IOCtor special types.

bool ROOT::TMetaUtils::HasIOConstructor(const clang::CXXRecordDecl *cl,
                                        std::string& arg,
                                        const RConstructorTypes& ctorTypes,
                                        const cling::Interpreter &interp)
{
   if (cl->isAbstract()) return false;

   for (auto & ctorType : ctorTypes) {

      auto ioCtorCat = ROOT::TMetaUtils::CheckConstructor(cl, ctorType, interp);

      if (EIOCtorCategory::kAbsent == ioCtorCat)
         continue;

      std::string proto( ctorType.GetName() );
      bool defaultCtor = proto.empty();
      if (defaultCtor) {
         arg.clear();
      } else {
         // I/O constructors can take pointers or references to ctorTypes
        proto += " *";
        if (EIOCtorCategory::kIOPtrType == ioCtorCat) {
           arg = "( ("; //(MyType*)nullptr
        } else if (EIOCtorCategory::kIORefType == ioCtorCat) {
           arg = "( *("; //*(MyType*)nullptr
        }
        arg += proto;
        arg += ")nullptr )";
      }
      // Check for private operator new
      const clang::CXXMethodDecl *method
         = GetMethodWithProto(cl, "operator new", "size_t", interp,
                              cling::LookupHelper::NoDiagnostics);
      if (method && method->getAccess() != clang::AS_public) {
         // The non-public op new is not going to improve for other c'tors.
         return false;
      }

      // This one looks good!
      return true;
   }
   return false;
}

////////////////////////////////////////////////////////////////////////////////

bool ROOT::TMetaUtils::NeedDestructor(const clang::CXXRecordDecl *cl,
                                      const cling::Interpreter& interp)
{
   if (!cl) return false;

   if (cl->hasUserDeclaredDestructor()) {

      cling::Interpreter::PushTransactionRAII clingRAII(const_cast<cling::Interpreter*>(&interp));
      clang::CXXDestructorDecl *dest = cl->getDestructor();
      if (dest) {
         return (dest->getAccess() == clang::AS_public);
      } else {
         return true; // no destructor, so let's assume it means default?
      }
   }
   return true;
}

////////////////////////////////////////////////////////////////////////////////
/// Return true, if the function (defined by the name and prototype) exists and is public

bool ROOT::TMetaUtils::CheckPublicFuncWithProto(const clang::CXXRecordDecl *cl,
                                                const char *methodname,
                                                const char *proto,
                                                const cling::Interpreter &interp,
                                                bool diagnose)
{
   const clang::CXXMethodDecl *method
      = GetMethodWithProto(cl,methodname,proto, interp,
                           diagnose ? cling::LookupHelper::WithDiagnostics
                           : cling::LookupHelper::NoDiagnostics);
   return (method && method->getAccess() == clang::AS_public);
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if the class has a method DirectoryAutoAdd(TDirectory *)

bool ROOT::TMetaUtils::HasDirectoryAutoAdd(const clang::CXXRecordDecl *cl, const cling::Interpreter &interp)
{
   // Detect if the class has a DirectoryAutoAdd

   // Detect if the class or one of its parent has a DirectoryAutoAdd
   const char *proto = "TDirectory*";
   const char *name = "DirectoryAutoAdd";

   return CheckPublicFuncWithProto(cl,name,proto,interp, false /*diags*/);
}


////////////////////////////////////////////////////////////////////////////////
/// Return true if the class has a method Merge(TCollection*,TFileMergeInfo*)

bool ROOT::TMetaUtils::HasNewMerge(const clang::CXXRecordDecl *cl, const cling::Interpreter &interp)
{
   // Detect if the class has a 'new' Merge function.

   // Detect if the class or one of its parent has a DirectoryAutoAdd
   const char *proto = "TCollection*,TFileMergeInfo*";
   const char *name = "Merge";

   return CheckPublicFuncWithProto(cl,name,proto,interp, false /*diags*/);
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if the class has a method Merge(TCollection*)

bool ROOT::TMetaUtils::HasOldMerge(const clang::CXXRecordDecl *cl, const cling::Interpreter &interp)
{
   // Detect if the class has an old fashion Merge function.

   // Detect if the class or one of its parent has a DirectoryAutoAdd
   const char *proto = "TCollection*";
   const char *name = "Merge";

   return CheckPublicFuncWithProto(cl,name,proto, interp, false /*diags*/);
}


////////////////////////////////////////////////////////////////////////////////
/// Return true if the class has a method ResetAfterMerge(TFileMergeInfo *)

bool ROOT::TMetaUtils::HasResetAfterMerge(const clang::CXXRecordDecl *cl, const cling::Interpreter &interp)
{
   // Detect if the class has a 'new' Merge function.
   // bool hasMethod = cl.HasMethod("DirectoryAutoAdd");

   // Detect if the class or one of its parent has a DirectoryAutoAdd
   const char *proto = "TFileMergeInfo*";
   const char *name = "ResetAfterMerge";

   return CheckPublicFuncWithProto(cl,name,proto, interp, false /*diags*/);
}


////////////////////////////////////////////////////////////////////////////////
/// Return true if the class has a custom member function streamer.

bool ROOT::TMetaUtils::HasCustomStreamerMemberFunction(const AnnotatedRecordDecl &cl,
                                                       const clang::CXXRecordDecl* clxx,
                                                       const cling::Interpreter &interp,
                                                       const TNormalizedCtxt &normCtxt)
{
   static const char *proto = "TBuffer&";

   const clang::CXXMethodDecl *method
      = GetMethodWithProto(clxx,"Streamer",proto, interp,
                           cling::LookupHelper::NoDiagnostics);
   const clang::DeclContext *clxx_as_context = llvm::dyn_cast<clang::DeclContext>(clxx);

   return (method && method->getDeclContext() == clxx_as_context
           && ( cl.RequestNoStreamer() || !cl.RequestStreamerInfo()));
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if the class has a custom member function streamer.

bool ROOT::TMetaUtils::HasCustomConvStreamerMemberFunction(const AnnotatedRecordDecl &cl,
                                                           const clang::CXXRecordDecl* clxx,
                                                           const cling::Interpreter &interp,
                                                           const TNormalizedCtxt &normCtxt)
{
   static const char *proto = "TBuffer&,TClass*";

   const clang::CXXMethodDecl *method
      = GetMethodWithProto(clxx,"Streamer",proto, interp,
                        cling::LookupHelper::NoDiagnostics);
   const clang::DeclContext *clxx_as_context = llvm::dyn_cast<clang::DeclContext>(clxx);

   return (method && method->getDeclContext() == clxx_as_context
           && ( cl.RequestNoStreamer() || !cl.RequestStreamerInfo()));
}


////////////////////////////////////////////////////////////////////////////////
/// Main implementation relying on GetFullyQualifiedTypeName
/// All other GetQualifiedName functions leverage this one except the
/// one for namespaces.

void ROOT::TMetaUtils::GetQualifiedName(std::string &qual_name, const clang::QualType &type, const clang::NamedDecl &forcontext)
{
   ROOT::TMetaUtils::GetFullyQualifiedTypeName(qual_name, type, forcontext.getASTContext());
}

//----
std::string ROOT::TMetaUtils::GetQualifiedName(const clang::QualType &type, const clang::NamedDecl &forcontext)
{
   std::string result;
   ROOT::TMetaUtils::GetQualifiedName(result,
                                      type,
                                      forcontext);
   return result;
}


////////////////////////////////////////////////////////////////////////////////

void  ROOT::TMetaUtils::GetQualifiedName(std::string& qual_type, const clang::Type &type, const clang::NamedDecl &forcontext)
{
   clang::QualType qualType(&type,0);
   ROOT::TMetaUtils::GetQualifiedName(qual_type,
                                      qualType,
                                      forcontext);
}

//---
std::string ROOT::TMetaUtils::GetQualifiedName(const clang::Type &type, const clang::NamedDecl &forcontext)
{
   std::string result;
   ROOT::TMetaUtils::GetQualifiedName(result,
                                      type,
                                      forcontext);
   return result;
}

// //______________________________________________________________________________
// void ROOT::TMetaUtils::GetQualifiedName(std::string &qual_name, const clang::NamespaceDecl &cl)
// {
//    GetQualifiedName(qual_name,cl);
// }
//
// //----
// std::string ROOT::TMetaUtils::GetQualifiedName(const clang::NamespaceDecl &cl){
//    return GetQualifiedName(cl);
// }

////////////////////////////////////////////////////////////////////////////////
/// This implementation does not rely on GetFullyQualifiedTypeName

void ROOT::TMetaUtils::GetQualifiedName(std::string &qual_name, const clang::NamedDecl &cl)
{
   llvm::raw_string_ostream stream(qual_name);
   clang::PrintingPolicy policy( cl.getASTContext().getPrintingPolicy() );
   policy.SuppressTagKeyword = true; // Never get the class or struct keyword
   policy.SuppressUnwrittenScope = true; // Don't write the inline or anonymous namespace names.

   cl.getNameForDiagnostic(stream,policy,true);
   stream.flush(); // flush to string.

   if ( qual_name ==  "(anonymous " ) {
      size_t pos = qual_name.find(':');
      qual_name.erase(0,pos+2);
   }
}

//----
std::string ROOT::TMetaUtils::GetQualifiedName(const clang::NamedDecl &cl){
   std::string result;
   ROOT::TMetaUtils::GetQualifiedName(result, cl);
   return result;
}


////////////////////////////////////////////////////////////////////////////////

void ROOT::TMetaUtils::GetQualifiedName(std::string &qual_name, const clang::RecordDecl &recordDecl)
{
   const clang::Type* declType ( recordDecl.getTypeForDecl() );
   clang::QualType qualType(declType,0);
   ROOT::TMetaUtils::GetQualifiedName(qual_name,
                                      qualType,
                                      recordDecl);
}

//----
std::string ROOT::TMetaUtils::GetQualifiedName(const clang::RecordDecl &recordDecl)
{
   std::string result;
   ROOT::TMetaUtils::GetQualifiedName(result,recordDecl);
   return result;
}

////////////////////////////////////////////////////////////////////////////////

void ROOT::TMetaUtils::GetQualifiedName(std::string &qual_name, const ROOT::TMetaUtils::AnnotatedRecordDecl &annotated)
{
   ROOT::TMetaUtils::GetQualifiedName(qual_name, *annotated.GetRecordDecl());
}

//----
std::string ROOT::TMetaUtils::GetQualifiedName(const AnnotatedRecordDecl &annotated)
{
   std::string result;
   ROOT::TMetaUtils::GetQualifiedName(result, annotated);
   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Create the data member name-type map for given class

void ROOT::TMetaUtils::CreateNameTypeMap(const clang::CXXRecordDecl &cl, ROOT::MembersTypeMap_t& nameType)
{
   std::stringstream dims;
   std::string typenameStr;

   const clang::ASTContext& astContext =  cl.getASTContext();

   // Loop over the non static data member.
   for(clang::RecordDecl::field_iterator field_iter = cl.field_begin(), end = cl.field_end();
       field_iter != end;
       ++field_iter){
      // The CINT based code was filtering away static variables (they are not part of
      // the list starting with field_begin in clang), and const enums (which should
      // also not be part of this list).
      // It was also filtering out the 'G__virtualinfo' artificial member.

      typenameStr.clear();
      dims.str("");
      dims.clear();

      clang::QualType fieldType(field_iter->getType());
      if (fieldType->isConstantArrayType()) {
         const clang::ConstantArrayType *arrayType = llvm::dyn_cast<clang::ConstantArrayType>(fieldType.getTypePtr());
         while (arrayType) {
            dims << "[" << arrayType->getSize().getLimitedValue() << "]";
            fieldType = arrayType->getElementType();
            arrayType = llvm::dyn_cast<clang::ConstantArrayType>(arrayType->getArrayElementTypeNoTypeQual());
         }
      }

      GetFullyQualifiedTypeName(typenameStr, fieldType, astContext);
      nameType[field_iter->getName().str()] = ROOT::Internal::TSchemaType(typenameStr.c_str(),dims.str().c_str());
   }

   // And now the base classes
   // We also need to look at the base classes.
   for(clang::CXXRecordDecl::base_class_const_iterator iter = cl.bases_begin(), end = cl.bases_end();
       iter != end;
       ++iter){
      std::string basename( iter->getType()->getAsCXXRecordDecl()->getNameAsString() ); // Intentionally using only the unqualified name.
      nameType[basename] = ROOT::Internal::TSchemaType(basename.c_str(),"");
   }
}

////////////////////////////////////////////////////////////////////////////////

const clang::FunctionDecl *ROOT::TMetaUtils::GetFuncWithProto(const clang::Decl* cinfo,
                                                              const char *method,
                                                              const char *proto,
                                                              const cling::Interpreter &interp,
                                                              bool diagnose)
{
   return interp.getLookupHelper().findFunctionProto(cinfo, method, proto,
                                                     diagnose ? cling::LookupHelper::WithDiagnostics
                                                     : cling::LookupHelper::NoDiagnostics);
}

////////////////////////////////////////////////////////////////////////////////
/// It looks like the template specialization decl actually contains _less_ information
/// on the location of the code than the decl (in case where there is forward declaration,
/// that is what the specialization points to.
///
/// const clang::CXXRecordDecl* clxx = llvm::dyn_cast<clang::CXXRecordDecl>(decl);
/// if (clxx) {
///    switch(clxx->getTemplateSpecializationKind()) {
///       case clang::TSK_Undeclared:
///          // We want the default behavior
///          break;
///       case clang::TSK_ExplicitInstantiationDeclaration:
///       case clang::TSK_ExplicitInstantiationDefinition:
///       case clang::TSK_ImplicitInstantiation: {
///          // We want the location of the template declaration:
///          const clang::ClassTemplateSpecializationDecl *tmplt_specialization = llvm::dyn_cast<clang::ClassTemplateSpecializationDecl> (clxx);
///          if (tmplt_specialization) {
///             return GetLineNumber(const_cast< clang::ClassTemplateSpecializationDecl *>(tmplt_specialization)->getSpecializedTemplate());
///          }
///          break;
///       }
///       case clang::TSK_ExplicitSpecialization:
///          // We want the default behavior
///          break;
///       default:
///          break;
///    }
/// }

long ROOT::TMetaUtils::GetLineNumber(const clang::Decl *decl)
{
   clang::SourceLocation sourceLocation = decl->getLocation();
   clang::SourceManager& sourceManager = decl->getASTContext().getSourceManager();

  if (!sourceLocation.isValid() ) {
      return -1;
   }

   if (!sourceLocation.isFileID()) {
      sourceLocation = sourceManager.getExpansionRange(sourceLocation).second;
   }

   if (sourceLocation.isValid() && sourceLocation.isFileID()) {
      return sourceManager.getLineNumber(sourceManager.getFileID(sourceLocation),sourceManager.getFileOffset(sourceLocation));
   }
   else {
      return -1;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if the type is a Double32_t or Float16_t or
/// is a instance template that depends on Double32_t or Float16_t.

bool ROOT::TMetaUtils::hasOpaqueTypedef(clang::QualType instanceType, const ROOT::TMetaUtils::TNormalizedCtxt &normCtxt)
{
   while (llvm::isa<clang::PointerType>(instanceType.getTypePtr())
       || llvm::isa<clang::ReferenceType>(instanceType.getTypePtr()))
   {
      instanceType = instanceType->getPointeeType();
   }

   const clang::ElaboratedType* etype
      = llvm::dyn_cast<clang::ElaboratedType>(instanceType.getTypePtr());
   if (etype) {
      instanceType = clang::QualType(etype->getNamedType().getTypePtr(),0);
   }

   // There is no typedef to worried about, except for the opaque ones.

   // Technically we should probably used our own list with just
   // Double32_t and Float16_t
   if (normCtxt.GetTypeWithAlternative().count(instanceType.getTypePtr())) {
      return true;
   }


   bool result = false;
   const clang::CXXRecordDecl* clxx = instanceType->getAsCXXRecordDecl();
   if (clxx && clxx->getTemplateSpecializationKind() != clang::TSK_Undeclared) {
      // do the template thing.
      const clang::TemplateSpecializationType* TST
      = llvm::dyn_cast<const clang::TemplateSpecializationType>(instanceType.getTypePtr());
      if (TST==0) {
   //         std::string type_name;
   //         type_name =  GetQualifiedName( instanceType, *clxx );
   //         fprintf(stderr,"ERROR: Could not findS TST for %s\n",type_name.c_str());
         return false;
      }
      for(clang::TemplateSpecializationType::iterator
          I = TST->begin(), E = TST->end();
          I!=E; ++I)
      {
         if (I->getKind() == clang::TemplateArgument::Type) {
   //            std::string arg;
   //            arg = GetQualifiedName( I->getAsType(), *clxx );
   //            fprintf(stderr,"DEBUG: looking at %s\n", arg.c_str());
            result |= ROOT::TMetaUtils::hasOpaqueTypedef(I->getAsType(), normCtxt);
         }
      }
   }
   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if any of the argument is or contains a double32.

bool ROOT::TMetaUtils::hasOpaqueTypedef(const AnnotatedRecordDecl &cl,
                                        const cling::Interpreter &interp,
                                        const TNormalizedCtxt &normCtxt)
{
   const clang::CXXRecordDecl* clxx =  llvm::dyn_cast<clang::CXXRecordDecl>(cl.GetRecordDecl());
   if (clxx->getTemplateSpecializationKind() == clang::TSK_Undeclared) return 0;

   clang::QualType instanceType = interp.getLookupHelper().findType(cl.GetNormalizedName(),
                                                                    cling::LookupHelper::WithDiagnostics);
   if (instanceType.isNull()) {
      //Error(0,"Could not find the clang::Type for %s\n",cl.GetNormalizedName());
      return false;
   }

   return ROOT::TMetaUtils::hasOpaqueTypedef(instanceType, normCtxt);
}

////////////////////////////////////////////////////////////////////////////////
/// Extract attr string

int ROOT::TMetaUtils::extractAttrString(clang::Attr* attribute, std::string& attrString)
{
   clang::AnnotateAttr* annAttr = clang::dyn_cast<clang::AnnotateAttr>(attribute);
   if (!annAttr) {
      //TMetaUtils::Error(0,"Could not cast Attribute to AnnotatedAttribute\n");
      return 1;
   }
   attrString = annAttr->getAnnotation();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////

int ROOT::TMetaUtils::extractPropertyNameValFromString(const std::string attributeStr,std::string& attrName, std::string& attrValue)
{
   // if separator found, extract name and value
   size_t substrFound (attributeStr.find(propNames::separator));
   if (substrFound==std::string::npos) {
      //TMetaUtils::Error(0,"Could not find property name-value separator (%s)\n",ROOT::TMetaUtils::PropertyNameValSeparator.c_str());
      return 1;
   }
   size_t EndPart1 = attributeStr.find_first_of(propNames::separator)  ;
   attrName = attributeStr.substr(0, EndPart1);
   const int separatorLength(propNames::separator.size());
   attrValue = attributeStr.substr(EndPart1 + separatorLength);
   return 0;
}

////////////////////////////////////////////////////////////////////////////////

int ROOT::TMetaUtils::extractPropertyNameVal(clang::Attr* attribute, std::string& attrName, std::string& attrValue)
{
   std::string attrString;
   int ret = extractAttrString(attribute, attrString);
   if (0!=ret) return ret;
   return extractPropertyNameValFromString(attrString, attrName,attrValue);
}

////////////////////////////////////////////////////////////////////////////////
/// This routine counts on the "propName<separator>propValue" format

bool ROOT::TMetaUtils::ExtractAttrPropertyFromName(const clang::Decl& decl,
                                                   const std::string& propName,
                                                   std::string& propValue)
{
   for (clang::Decl::attr_iterator attrIt = decl.attr_begin();
    attrIt!=decl.attr_end();++attrIt){
      clang::AnnotateAttr* annAttr = clang::dyn_cast<clang::AnnotateAttr>(*attrIt);
      if (!annAttr) continue;

      llvm::StringRef attribute = annAttr->getAnnotation();
      std::pair<llvm::StringRef,llvm::StringRef> split = attribute.split(propNames::separator.c_str());
      if (split.first != propName.c_str()) continue;
      else {
         propValue = split.second;
         return true;
      }
   }
   return false;
}

////////////////////////////////////////////////////////////////////////////////
/// This routine counts on the "propName<separator>propValue" format

bool ROOT::TMetaUtils::ExtractAttrIntPropertyFromName(const clang::Decl& decl,
                                                      const std::string& propName,
                                                      int& propValue)
{
   for (clang::Decl::attr_iterator attrIt = decl.attr_begin();
    attrIt!=decl.attr_end();++attrIt){
      clang::AnnotateAttr* annAttr = clang::dyn_cast<clang::AnnotateAttr>(*attrIt);
      if (!annAttr) continue;

      llvm::StringRef attribute = annAttr->getAnnotation();
      std::pair<llvm::StringRef,llvm::StringRef> split = attribute.split(propNames::separator.c_str());
      if (split.first != propName.c_str()) continue;
      else {
         return split.second.getAsInteger(10,propValue);
      }
   }
   return false;
}

////////////////////////////////////////////////////////////////////////////////
/// FIXME: a function of 450+ lines!

void ROOT::TMetaUtils::WriteClassInit(std::ostream& finalString,
                                      const AnnotatedRecordDecl &cl,
                                      const clang::CXXRecordDecl *decl,
                                      const cling::Interpreter &interp,
                                      const TNormalizedCtxt &normCtxt,
                                      const RConstructorTypes& ctorTypes,
                                      bool& needCollectionProxy)
{
   std::string classname = TClassEdit::GetLong64_Name(cl.GetNormalizedName());

   std::string mappedname;
   ROOT::TMetaUtils::GetCppName(mappedname,classname.c_str());
   std::string csymbol = classname;
   std::string args;

   if ( ! TClassEdit::IsStdClass( classname.c_str() ) ) {

      // Prefix the full class name with '::' except for the STL
      // containers and std::string.  This is to request the
      // real class instead of the class in the namespace ROOT::Shadow
      csymbol.insert(0,"::");
   }

   int stl = TClassEdit::IsSTLCont(classname);
   bool bset = TClassEdit::IsSTLBitset(classname.c_str());

   bool isStd = TMetaUtils::IsStdClass(*decl);
   const cling::LookupHelper& lh = interp.getLookupHelper();
   bool isString = TMetaUtils::IsOfType(*decl,"std::string",lh);

   bool isStdNotString = isStd && !isString;

   finalString << "namespace ROOT {" << "\n";

   if (!ClassInfo__HasMethod(decl,"Dictionary",interp) || IsTemplate(*decl))
   {
      finalString << "   static TClass *" << mappedname.c_str() << "_Dictionary();\n"
                  << "   static void " << mappedname.c_str() << "_TClassManip(TClass*);\n";


   }

   if (HasIOConstructor(decl, args, ctorTypes, interp)) {
      finalString << "   static void *new_" << mappedname.c_str() << "(void *p = 0);" << "\n";

      if (args.size()==0 && NeedDestructor(decl, interp))
      {
         finalString << "   static void *newArray_";
         finalString << mappedname.c_str();
         finalString << "(Long_t size, void *p);";
         finalString << "\n";
      }
   }

   if (NeedDestructor(decl, interp)) {
      finalString << "   static void delete_" << mappedname.c_str() << "(void *p);" << "\n" << "   static void deleteArray_" << mappedname.c_str() << "(void *p);" << "\n" << "   static void destruct_" << mappedname.c_str() << "(void *p);" << "\n";
   }
   if (HasDirectoryAutoAdd(decl, interp)) {
      finalString << "   static void directoryAutoAdd_" << mappedname.c_str() << "(void *obj, TDirectory *dir);" << "\n";
   }
   if (HasCustomStreamerMemberFunction(cl, decl, interp, normCtxt)) {
      finalString << "   static void streamer_" << mappedname.c_str() << "(TBuffer &buf, void *obj);" << "\n";
   }
   if (HasCustomConvStreamerMemberFunction(cl, decl, interp, normCtxt)) {
      finalString << "   static void conv_streamer_" << mappedname.c_str() << "(TBuffer &buf, void *obj, const TClass*);" << "\n";
   }
   if (HasNewMerge(decl, interp) || HasOldMerge(decl, interp)) {
      finalString << "   static Long64_t merge_" << mappedname.c_str() << "(void *obj, TCollection *coll,TFileMergeInfo *info);" << "\n";
   }
   if (HasResetAfterMerge(decl, interp)) {
      finalString << "   static void reset_" << mappedname.c_str() << "(void *obj, TFileMergeInfo *info);" << "\n";
   }

   //--------------------------------------------------------------------------
   // Check if we have any schema evolution rules for this class
   /////////////////////////////////////////////////////////////////////////////

   ROOT::SchemaRuleClassMap_t::iterator rulesIt1 = ROOT::gReadRules.find( classname.c_str() );
   ROOT::SchemaRuleClassMap_t::iterator rulesIt2 = ROOT::gReadRawRules.find( classname.c_str() );

   ROOT::MembersTypeMap_t nameTypeMap;
   CreateNameTypeMap( *decl, nameTypeMap ); // here types for schema evo are written

   //--------------------------------------------------------------------------
   // Process the read rules
   /////////////////////////////////////////////////////////////////////////////

   if( rulesIt1 != ROOT::gReadRules.end() ) {
      int i = 0;
      finalString << "\n   // Schema evolution read functions\n";
      std::list<ROOT::SchemaRuleMap_t>::iterator rIt = rulesIt1->second.begin();
      while( rIt != rulesIt1->second.end() ) {

         //--------------------------------------------------------------------
         // Check if the rules refer to valid data members
         ///////////////////////////////////////////////////////////////////////

         std::string error_string;
         if( !HasValidDataMembers( *rIt, nameTypeMap, error_string ) ) {
            Warning(0, "%s", error_string.c_str());
            rIt = rulesIt1->second.erase(rIt);
            continue;
         }

         //---------------------------------------------------------------------
         // Write the conversion function if necessary
         ///////////////////////////////////////////////////////////////////////

         if( rIt->find( "code" ) != rIt->end() ) {
            WriteReadRuleFunc( *rIt, i++, mappedname, nameTypeMap, finalString );
         }
         ++rIt;
      }
   }




   //--------------------------------------------------------------------------
   // Process the read raw rules
   /////////////////////////////////////////////////////////////////////////////

   if( rulesIt2 != ROOT::gReadRawRules.end() ) {
      int i = 0;
      finalString << "\n   // Schema evolution read raw functions\n";
      std::list<ROOT::SchemaRuleMap_t>::iterator rIt = rulesIt2->second.begin();
      while( rIt != rulesIt2->second.end() ) {

         //--------------------------------------------------------------------
         // Check if the rules refer to valid data members
         ///////////////////////////////////////////////////////////////////////

         std::string error_string;
         if( !HasValidDataMembers( *rIt, nameTypeMap, error_string ) ) {
            Warning(0, "%s", error_string.c_str());
            rIt = rulesIt2->second.erase(rIt);
            continue;
         }

         //---------------------------------------------------------------------
         // Write the conversion function
         ///////////////////////////////////////////////////////////////////////

         if( rIt->find( "code" ) == rIt->end() )
            continue;

         WriteReadRawRuleFunc( *rIt, i++, mappedname, nameTypeMap, finalString );
         ++rIt;
      }
   }

   finalString << "\n" << "   // Function generating the singleton type initializer" << "\n";

   finalString << "   static TGenericClassInfo *GenerateInitInstanceLocal(const " << csymbol << "*)" << "\n" << "   {" << "\n";

   finalString << "      " << csymbol << " *ptr = 0;" << "\n";

   //fprintf(fp, "      static ::ROOT::ClassInfo< %s > \n",classname.c_str());
   if (ClassInfo__HasMethod(decl,"IsA",interp) ) {
      finalString << "      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< "  << csymbol << " >(0);" << "\n";
   }
   else {
      finalString << "      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(" << csymbol << "));" << "\n";
   }
   finalString << "      static ::ROOT::TGenericClassInfo " << "\n" << "         instance(\"" << classname.c_str() << "\", ";

   if (ClassInfo__HasMethod(decl,"Class_Version",interp)) {
      finalString << csymbol << "::Class_Version(), ";
   } else if (bset) {
      finalString << "2, "; // bitset 'version number'
   } else if (stl) {
      finalString << "-2, "; // "::TStreamerInfo::Class_Version(), ";
   } else if( cl.HasClassVersion() ) {
      finalString << cl.RequestedVersionNumber() << ", ";
   } else { // if (cl_input.RequestStreamerInfo()) {

      // Need to find out if the operator>> is actually defined for this class.
      static const char *versionFunc = "GetClassVersion";
      //      int ncha = strlen(classname.c_str())+strlen(versionFunc)+5;
      //      char *funcname= new char[ncha];
      //      snprintf(funcname,ncha,"%s<%s >",versionFunc,classname.c_str());
      std::string proto = classname + "*";
      const clang::Decl* ctxt = llvm::dyn_cast<clang::Decl>((*cl).getDeclContext());
      const clang::FunctionDecl *methodinfo
         = ROOT::TMetaUtils::GetFuncWithProto(ctxt, versionFunc, proto.c_str(),
                                              interp, cling::LookupHelper::NoDiagnostics);
      //      delete [] funcname;

      if (methodinfo &&
          ROOT::TMetaUtils::GetFileName(*methodinfo, interp).find("Rtypes.h") == llvm::StringRef::npos) {

         // GetClassVersion was defined in the header file.
         //fprintf(fp, "GetClassVersion((%s *)0x0), ",classname.c_str());
         finalString << "GetClassVersion< ";
         finalString << classname.c_str();
         finalString << " >(), ";
      }
      //static char temporary[1024];
      //sprintf(temporary,"GetClassVersion<%s>( (%s *) 0x0 )",classname.c_str(),classname.c_str());
      //fprintf(stderr,"DEBUG: %s has value %d\n",classname.c_str(),(int)G__int(G__calc(temporary)));
   }

   std::string filename = ROOT::TMetaUtils::GetFileName(*cl, interp);
   if (filename.length() > 0) {
      for (unsigned int i=0; i<filename.length(); i++) {
         if (filename[i]=='\\') filename[i]='/';
      }
   }
   finalString << "\"" << filename << "\", " << ROOT::TMetaUtils::GetLineNumber(cl)
      << "," << "\n" << "                  typeid(" << csymbol
      << "), ::ROOT::Internal::DefineBehavior(ptr, ptr)," << "\n" << "                  ";

   if (ClassInfo__HasMethod(decl,"Dictionary",interp) && !IsTemplate(*decl)) {
      finalString << "&" << csymbol << "::Dictionary, ";
   } else {
      finalString << "&" << mappedname << "_Dictionary, ";
   }

   enum {
      TClassTable__kHasCustomStreamerMember = 0x10 // See TClassTable.h
   };

   Int_t rootflag = cl.RootFlag();
   if (HasCustomStreamerMemberFunction(cl, decl, interp, normCtxt)) {
      rootflag = rootflag | TClassTable__kHasCustomStreamerMember;
   }
   finalString << "isa_proxy, " << rootflag << "," << "\n" << "                  sizeof(" << csymbol << ") );" << "\n";
   if (HasIOConstructor(decl, args, ctorTypes, interp)) {
      finalString << "      instance.SetNew(&new_" << mappedname.c_str() << ");" << "\n";
      if (args.size()==0 && NeedDestructor(decl, interp))
         finalString << "      instance.SetNewArray(&newArray_" << mappedname.c_str() << ");" << "\n";
   }
   if (NeedDestructor(decl, interp)) {
      finalString << "      instance.SetDelete(&delete_" << mappedname.c_str() << ");" << "\n" << "      instance.SetDeleteArray(&deleteArray_" << mappedname.c_str() << ");" << "\n" << "      instance.SetDestructor(&destruct_" << mappedname.c_str() << ");" << "\n";
   }
   if (HasDirectoryAutoAdd(decl, interp)) {
      finalString << "      instance.SetDirectoryAutoAdd(&directoryAutoAdd_" << mappedname.c_str() << ");" << "\n";
   }
   if (HasCustomStreamerMemberFunction(cl, decl, interp, normCtxt)) {
      // We have a custom member function streamer or an older (not StreamerInfo based) automatic streamer.
      finalString << "      instance.SetStreamerFunc(&streamer_" << mappedname.c_str() << ");" << "\n";
   }
   if (HasCustomConvStreamerMemberFunction(cl, decl, interp, normCtxt)) {
      // We have a custom member function streamer or an older (not StreamerInfo based) automatic streamer.
      finalString << "      instance.SetConvStreamerFunc(&conv_streamer_" << mappedname.c_str() << ");" << "\n";
   }
   if (HasNewMerge(decl, interp) || HasOldMerge(decl, interp)) {
      finalString << "      instance.SetMerge(&merge_" << mappedname.c_str() << ");" << "\n";
   }
   if (HasResetAfterMerge(decl, interp)) {
      finalString << "      instance.SetResetAfterMerge(&reset_" << mappedname.c_str() << ");" << "\n";
   }
   if (bset) {
      finalString << "      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::" << "Pushback" << "<Internal::TStdBitsetHelper< " << classname.c_str() << " > >()));" << "\n";

      needCollectionProxy = true;
   } else if (stl != 0 &&
              ((stl > 0 && stl<ROOT::kSTLend) || (stl < 0 && stl>-ROOT::kSTLend)) && // is an stl container
              (stl != ROOT::kSTLbitset && stl !=-ROOT::kSTLbitset) ){     // is no bitset
      int idx = classname.find("<");
      int stlType = (idx!=(int)std::string::npos) ? TClassEdit::STLKind(classname.substr(0,idx)) : 0;
      const char* methodTCP=0;
      switch(stlType)  {
         case ROOT::kSTLvector:
         case ROOT::kSTLlist:
         case ROOT::kSTLdeque:
            methodTCP="Pushback";
            break;
         case ROOT::kSTLforwardlist:
            methodTCP="Pushfront";
            break;
         case ROOT::kSTLmap:
         case ROOT::kSTLmultimap:
         case ROOT::kSTLunorderedmap:
         case ROOT::kSTLunorderedmultimap:
            methodTCP="MapInsert";
            break;
         case ROOT::kSTLset:
         case ROOT::kSTLmultiset:
         case ROOT::kSTLunorderedset:
         case ROOT::kSTLunorderedmultiset:
            methodTCP="Insert";
            break;
      }
      // FIXME Workaround: for the moment we do not generate coll proxies with unique ptrs since
      // they imply copies and therefore do not compile.
      auto classNameForIO = TClassEdit::GetNameForIO(classname);
      finalString << "      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::" << methodTCP << "< " << classNameForIO.c_str() << " >()));" << "\n";

      needCollectionProxy = true;
   }

   //---------------------------------------------------------------------------
   // Register Alternate spelling of the class name.
   /////////////////////////////////////////////////////////////////////////////

   if (cl.GetRequestedName()[0] && classname != cl.GetRequestedName()) {
      finalString << "\n" << "      ::ROOT::AddClassAlternate(\""
                  << classname << "\",\"" << cl.GetRequestedName() << "\");\n";
   }

   if (!cl.GetDemangledTypeInfo().empty()
         && cl.GetDemangledTypeInfo() != classname
         && cl.GetDemangledTypeInfo() != cl.GetRequestedName()) {
      finalString << "\n" << "      ::ROOT::AddClassAlternate(\""
                  << classname << "\",\"" << cl.GetDemangledTypeInfo() << "\");\n";

   }

   //---------------------------------------------------------------------------
   // Pass the schema evolution rules to TGenericClassInfo
   /////////////////////////////////////////////////////////////////////////////

   if( (rulesIt1 != ROOT::gReadRules.end() && rulesIt1->second.size()>0) || (rulesIt2 != ROOT::gReadRawRules.end()  && rulesIt2->second.size()>0) ) {
      finalString << "\n" << "      ROOT::Internal::TSchemaHelper* rule;" << "\n";
   }

   if( rulesIt1 != ROOT::gReadRules.end() ) {
      finalString << "\n" << "      // the io read rules" << "\n" << "      std::vector<ROOT::Internal::TSchemaHelper> readrules(" << rulesIt1->second.size() << ");" << "\n";
      ROOT::WriteSchemaList( rulesIt1->second, "readrules", finalString );
      finalString << "      instance.SetReadRules( readrules );" << "\n";
   }

   if( rulesIt2 != ROOT::gReadRawRules.end() ) {
      finalString << "\n" << "      // the io read raw rules" << "\n" << "      std::vector<ROOT::Internal::TSchemaHelper> readrawrules(" << rulesIt2->second.size() << ");" << "\n";
      ROOT::WriteSchemaList( rulesIt2->second, "readrawrules", finalString );
      finalString << "      instance.SetReadRawRules( readrawrules );" << "\n";
   }

   finalString << "      return &instance;" << "\n" << "   }" << "\n";

   if (!isStdNotString && !ROOT::TMetaUtils::hasOpaqueTypedef(cl, interp, normCtxt)) {
      // The GenerateInitInstance for STL are not unique and should not be externally accessible
      finalString << "   TGenericClassInfo *GenerateInitInstance(const " << csymbol << "*)" << "\n" << "   {\n      return GenerateInitInstanceLocal((" << csymbol << "*)0);\n   }" << "\n";
   }

   finalString << "   // Static variable to force the class initialization" << "\n";
   // must be one long line otherwise UseDummy does not work


   finalString << "   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const " << csymbol << "*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));" << "\n";

   if (!ClassInfo__HasMethod(decl,"Dictionary",interp) || IsTemplate(*decl)) {
      finalString <<  "\n" << "   // Dictionary for non-ClassDef classes" << "\n"
                  << "   static TClass *" << mappedname << "_Dictionary() {\n"
                  << "      TClass* theClass ="
                  << "::ROOT::GenerateInitInstanceLocal((const " << csymbol << "*)0x0)->GetClass();\n"
                  << "      " << mappedname << "_TClassManip(theClass);\n";
      finalString << "   return theClass;\n";
      finalString << "   }\n\n";

      // Now manipulate tclass in order to percolate the properties expressed as
      // annotations of the decls.
      std::string manipString;
      std::string attribute_s;
      std::string attrName, attrValue;
      // Class properties
      bool attrMapExtracted = false;
      if (decl->hasAttrs()){
         // Loop on the attributes
         for (clang::Decl::attr_iterator attrIt = decl->attr_begin();
              attrIt!=decl->attr_end();++attrIt){
            if ( 0!=ROOT::TMetaUtils::extractAttrString(*attrIt,attribute_s)){
               continue;
            }
            if (0!=ROOT::TMetaUtils::extractPropertyNameValFromString(attribute_s, attrName, attrValue)){
               continue;
            }
            if (attrName == "name" ||
                attrName == "pattern" ||
                attrName == "rootmap") continue;
            // A general property
            // 1) We need to create the property map (in the gen code)
            // 2) we need to take out the map (in the gen code)
            // 3) We need to bookkep the fact that the map is created and out (in this source)
            // 4) We fill the map (in the gen code)
            if (!attrMapExtracted){
               manipString+="      theClass->CreateAttributeMap();\n";
               manipString+="      TDictAttributeMap* attrMap( theClass->GetAttributeMap() );\n";
               attrMapExtracted=true;
            }
            manipString+="      attrMap->AddProperty(\""+attrName +"\",\""+attrValue+"\");\n";
         }
      } // end of class has properties

      // Member properties
      // Loop on declarations inside the class, including data members
      for(clang::CXXRecordDecl::decl_iterator internalDeclIt = decl->decls_begin();
          internalDeclIt != decl->decls_end(); ++internalDeclIt){
         if (!(!(*internalDeclIt)->isImplicit()
            && (clang::isa<clang::FieldDecl>(*internalDeclIt) ||
                clang::isa<clang::VarDecl>(*internalDeclIt)))) continue; // Check if it's a var or a field

         // Now let's check the attributes of the var/field
         if (!internalDeclIt->hasAttrs()) continue;

         attrMapExtracted = false;
         bool memberPtrCreated = false;

         for (clang::Decl::attr_iterator attrIt = internalDeclIt->attr_begin();
              attrIt!=internalDeclIt->attr_end();++attrIt){

            // Get the attribute as string
            if ( 0!=ROOT::TMetaUtils::extractAttrString(*attrIt,attribute_s)){
               continue;
            }

            // Check the name of the decl
            clang::NamedDecl* namedInternalDecl = clang::dyn_cast<clang::NamedDecl> (*internalDeclIt);
            if (!namedInternalDecl) {
               TMetaUtils::Error(0,"Cannot convert field declaration to clang::NamedDecl");
               continue;
            };
            const std::string memberName(namedInternalDecl->getName());
            const std::string cppMemberName="theMember_"+memberName;

            // Prepare a string to get the data member, it can be used later.
            const std::string dataMemberCreation= "      TDataMember* "+cppMemberName+" = theClass->GetDataMember(\""+memberName+"\");\n";

            // Let's now attack regular properties

            if (0!=ROOT::TMetaUtils::extractPropertyNameValFromString(attribute_s, attrName, attrValue)){
               continue;
            }

            // Skip these
            if (attrName == propNames::comment ||
                attrName == propNames::iotype ||
                attrName == propNames::ioname ) continue;

            if (!memberPtrCreated){
               manipString+=dataMemberCreation;
               memberPtrCreated=true;
            }

            if (!attrMapExtracted){
               manipString+="      "+cppMemberName+"->CreateAttributeMap();\n";
               manipString+="      TDictAttributeMap* memberAttrMap_"+memberName+"( theMember_"+memberName+"->GetAttributeMap() );\n";
            attrMapExtracted=true;
            }

            manipString+="      memberAttrMap_"+memberName+"->AddProperty(\""+attrName +"\",\""+attrValue+"\");\n";


         } // End loop on attributes
      } // End loop on internal declarations


      finalString << "   static void " << mappedname << "_TClassManip(TClass* " << (manipString.empty() ? "":"theClass") << "){\n"
                  << manipString
                  << "   }\n\n";
   } // End of !ClassInfo__HasMethod(decl,"Dictionary") || IsTemplate(*decl))

   finalString << "} // end of namespace ROOT" << "\n" << "\n";
}

////////////////////////////////////////////////////////////////////////////////
/// Return true if one of the class' enclosing scope is a namespace and
/// set fullname to the fully qualified name,
/// clsname to the name within a namespace
/// and nsname to the namespace fully qualified name.

bool ROOT::TMetaUtils::GetNameWithinNamespace(std::string &fullname,
                                                 std::string &clsname,
                                                 std::string &nsname,
                                                 const clang::CXXRecordDecl *cl)
{
   fullname.clear();
   nsname.clear();

   ROOT::TMetaUtils::GetQualifiedName(fullname,*cl);
   clsname = fullname;

   // Inline namespace are stripped from the normalized name, we need to
   // strip it from the prefix we want to remove.
   auto ctxt = cl->getEnclosingNamespaceContext();
   while(ctxt && ctxt!=cl && ctxt->isInlineNamespace()) {
      ctxt = ctxt->getParent();
   }
   if (ctxt) {
      const clang::NamedDecl *namedCtxt = llvm::dyn_cast<clang::NamedDecl>(ctxt);
      if (namedCtxt && namedCtxt!=cl) {
         const clang::NamespaceDecl *nsdecl = llvm::dyn_cast<clang::NamespaceDecl>(namedCtxt);
         if (nsdecl != 0 && !nsdecl->isAnonymousNamespace()) {
            ROOT::TMetaUtils::GetQualifiedName(nsname,*nsdecl);
            clsname.erase (0, nsname.size() + 2);
            return true;
         }
      }
   }
   return false;
}

////////////////////////////////////////////////////////////////////////////////

const clang::DeclContext *GetEnclosingSpace(const clang::RecordDecl &cl)
{
   const clang::DeclContext *ctxt = cl.getDeclContext();
   while(ctxt && !ctxt->isNamespace()) {
      ctxt = ctxt->getParent();
   }
   return ctxt;
}

////////////////////////////////////////////////////////////////////////////////
/// Write all the necessary opening part of the namespace and
/// return the number of closing brackets needed
/// For example for Space1::Space2
/// we write: namespace Space1 { namespace Space2 {
/// and return 2.

int ROOT::TMetaUtils::WriteNamespaceHeader(std::ostream &out, const clang::DeclContext *ctxt)
{
   int closing_brackets = 0;

   //fprintf(stderr,"DEBUG: in WriteNamespaceHeader for %s with %s\n",
   //    cl.Fullname(),namespace_obj.Fullname());
   if (ctxt && ctxt->isNamespace()) {
      closing_brackets = WriteNamespaceHeader(out,ctxt->getParent());
      for (int indent = 0; indent < closing_brackets; ++indent) {
         out << "   ";
      }
      const clang::NamespaceDecl *ns = llvm::dyn_cast<clang::NamespaceDecl>(ctxt);
      if (ns->isInline())
         out << "inline ";
      out << "namespace " << ns->getNameAsString() << " {" << std::endl;
      closing_brackets++;
   }

   return closing_brackets;
}

////////////////////////////////////////////////////////////////////////////////

int ROOT::TMetaUtils::WriteNamespaceHeader(std::ostream &out, const clang::RecordDecl *cl)
{
   return WriteNamespaceHeader(out, GetEnclosingSpace(*cl));
}

////////////////////////////////////////////////////////////////////////////////

bool ROOT::TMetaUtils::NeedTemplateKeyword(const clang::CXXRecordDecl *cl)
{
   clang::TemplateSpecializationKind kind = cl->getTemplateSpecializationKind();
   if (kind == clang::TSK_Undeclared ) {
      // Note a template;
      return false;
   } else if (kind == clang::TSK_ExplicitSpecialization) {
      // This is a specialized templated class
      return false;
   } else {
      // This is an automatically or explicitly instantiated templated class.
      return true;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// return true if we can find a custom operator new with placement

bool ROOT::TMetaUtils::HasCustomOperatorNewPlacement(const char *which, const clang::RecordDecl &cl, const cling::Interpreter &interp)
{
   const char *name = which;
   const char *proto = "size_t";
   const char *protoPlacement = "size_t,void*";

   // First search in the enclosing namespaces
   const clang::FunctionDecl *operatornew
      = ROOT::TMetaUtils::GetFuncWithProto(llvm::dyn_cast<clang::Decl>(cl.getDeclContext()),
                                           name, proto, interp,
                                           cling::LookupHelper::NoDiagnostics);
   const clang::FunctionDecl *operatornewPlacement
      = ROOT::TMetaUtils::GetFuncWithProto(llvm::dyn_cast<clang::Decl>(cl.getDeclContext()),
                                           name, protoPlacement, interp,
                                           cling::LookupHelper::NoDiagnostics);

   const clang::DeclContext *ctxtnew = 0;
   const clang::DeclContext *ctxtnewPlacement = 0;

   if (operatornew) {
      ctxtnew = operatornew->getParent();
   }
   if (operatornewPlacement) {
      ctxtnewPlacement = operatornewPlacement->getParent();
   }

   // Then in the class and base classes
   operatornew = ROOT::TMetaUtils::GetFuncWithProto(&cl, name, proto, interp,
                                                    false /*diags*/);
   operatornewPlacement
      = ROOT::TMetaUtils::GetFuncWithProto(&cl, name, protoPlacement, interp,
                                           false /*diags*/);

   if (operatornew) {
      ctxtnew = operatornew->getParent();
   }
   if (operatornewPlacement) {
      ctxtnewPlacement = operatornewPlacement->getParent();
   }

   if (ctxtnewPlacement == 0) {
      return false;
   }
   if (ctxtnew == 0) {
      // Only a new with placement, no hiding
      return true;
   }
   // Both are non zero
   if (ctxtnew == ctxtnewPlacement) {
      // Same declaration ctxt, no hiding
      return true;
   }
   const clang::CXXRecordDecl* clnew = llvm::dyn_cast<clang::CXXRecordDecl>(ctxtnew);
   const clang::CXXRecordDecl* clnewPlacement = llvm::dyn_cast<clang::CXXRecordDecl>(ctxtnewPlacement);
   if (clnew == 0 && clnewPlacement == 0) {
      // They are both in different namespaces, I am not sure of the rules.
      // we probably ought to find which one is closest ... for now bail
      // (because rootcling was also bailing on that).
      return true;
   }
   if (clnew != 0 && clnewPlacement == 0) {
      // operator new is class method hiding the outer scope operator new with placement.
      return false;
   }
   if (clnew == 0 && clnewPlacement != 0) {
      // operator new is a not class method and can not hide new with placement which is a method
      return true;
   }
   // Both are class methods
   if (clnew->isDerivedFrom(clnewPlacement)) {
      // operator new is in a more derived part of the hierarchy, it is hiding operator new with placement.
      return false;
   }
   // operator new with placement is in a more derived part of the hierarchy, it can't be hidden by operator new.
   return true;
}

////////////////////////////////////////////////////////////////////////////////
/// return true if we can find a custom operator new with placement

bool ROOT::TMetaUtils::HasCustomOperatorNewPlacement(const clang::RecordDecl &cl, const cling::Interpreter &interp)
{
   return HasCustomOperatorNewPlacement("operator new",cl, interp);
}

////////////////////////////////////////////////////////////////////////////////
/// return true if we can find a custom operator new with placement

bool ROOT::TMetaUtils::HasCustomOperatorNewArrayPlacement(const clang::RecordDecl &cl, const cling::Interpreter &interp)
{
   return HasCustomOperatorNewPlacement("operator new[]",cl, interp);
}

////////////////////////////////////////////////////////////////////////////////
/// std::string NormalizedName;
/// GetNormalizedName(NormalizedName, decl->getASTContext().getTypeDeclType(decl), interp, normCtxt);

void ROOT::TMetaUtils::WriteAuxFunctions(std::ostream& finalString,
                                         const AnnotatedRecordDecl &cl,
                                         const clang::CXXRecordDecl *decl,
                                         const cling::Interpreter &interp,
                                         const RConstructorTypes& ctorTypes,
                                         const TNormalizedCtxt &normCtxt)
{
   std::string classname = TClassEdit::GetLong64_Name(cl.GetNormalizedName());

   std::string mappedname;
   ROOT::TMetaUtils::GetCppName(mappedname,classname.c_str());

   // Write the functions that are need for the TGenericClassInfo.
   // This includes
   //    IsA
   //    operator new
   //    operator new[]
   //    operator delete
   //    operator delete[]

   ROOT::TMetaUtils::GetCppName(mappedname,classname.c_str());

   if ( ! TClassEdit::IsStdClass( classname.c_str() ) ) {

      // Prefix the full class name with '::' except for the STL
      // containers and std::string.  This is to request the
      // real class instead of the class in the namespace ROOT::Shadow
      classname.insert(0,"::");
   }

   finalString << "namespace ROOT {" << "\n";

   std::string args;
   if (HasIOConstructor(decl, args, ctorTypes, interp)) {
      // write the constructor wrapper only for concrete classes
      finalString << "   // Wrappers around operator new" << "\n";
      finalString << "   static void *new_" << mappedname.c_str() << "(void *p) {" << "\n" << "      return  p ? ";
      if (HasCustomOperatorNewPlacement(*decl, interp)) {
         finalString << "new(p) ";
         finalString << classname.c_str();
         finalString << args;
         finalString << " : ";
      } else {
         finalString << "::new((::ROOT::Internal::TOperatorNewHelper*)p) ";
         finalString << classname.c_str();
         finalString << args;
         finalString << " : ";
      }
      finalString << "new " << classname.c_str() << args << ";" << "\n";
      finalString << "   }" << "\n";

      if (args.size()==0 && NeedDestructor(decl, interp)) {
         // Can not can newArray if the destructor is not public.
         finalString << "   static void *newArray_";
         finalString << mappedname.c_str();
         finalString << "(Long_t nElements, void *p) {";
         finalString << "\n";
         finalString << "      return p ? ";
         if (HasCustomOperatorNewArrayPlacement(*decl, interp)) {
            finalString << "new(p) ";
            finalString << classname.c_str();
            finalString << "[nElements] : ";
         } else {
            finalString << "::new((::ROOT::Internal::TOperatorNewHelper*)p) ";
            finalString << classname.c_str();
            finalString << "[nElements] : ";
         }
         finalString << "new ";
         finalString << classname.c_str();
         finalString << "[nElements];";
         finalString << "\n";
         finalString << "   }";
         finalString << "\n";
      }
   }

   if (NeedDestructor(decl, interp)) {
      finalString << "   // Wrapper around operator delete" << "\n" << "   static void delete_" << mappedname.c_str() << "(void *p) {" << "\n" << "      delete ((" << classname.c_str() << "*)p);" << "\n" << "   }" << "\n" << "   static void deleteArray_" << mappedname.c_str() << "(void *p) {" << "\n" << "      delete [] ((" << classname.c_str() << "*)p);" << "\n" << "   }" << "\n" << "   static void destruct_" << mappedname.c_str() << "(void *p) {" << "\n" << "      typedef " << classname.c_str() << " current_t;" << "\n" << "      ((current_t*)p)->~current_t();" << "\n" << "   }" << "\n";
   }

   if (HasDirectoryAutoAdd(decl, interp)) {
       finalString << "   // Wrapper around the directory auto add." << "\n" << "   static void directoryAutoAdd_" << mappedname.c_str() << "(void *p, TDirectory *dir) {" << "\n" << "      ((" << classname.c_str() << "*)p)->DirectoryAutoAdd(dir);" << "\n" << "   }" << "\n";
   }

   if (HasCustomStreamerMemberFunction(cl, decl, interp, normCtxt)) {
      finalString << "   // Wrapper around a custom streamer member function." << "\n" << "   static void streamer_" << mappedname.c_str() << "(TBuffer &buf, void *obj) {" << "\n" << "      ((" << classname.c_str() << "*)obj)->" << classname.c_str() << "::Streamer(buf);" << "\n" << "   }" << "\n";
   }

   if (HasCustomConvStreamerMemberFunction(cl, decl, interp, normCtxt)) {
      finalString << "   // Wrapper around a custom streamer member function." << "\n" << "   static void conv_streamer_" << mappedname.c_str() << "(TBuffer &buf, void *obj, const TClass *onfile_class) {" << "\n" << "      ((" << classname.c_str() << "*)obj)->" << classname.c_str() << "::Streamer(buf,onfile_class);" << "\n" << "   }" << "\n";
   }

   if (HasNewMerge(decl, interp)) {
      finalString << "   // Wrapper around the merge function." << "\n" << "   static Long64_t merge_" << mappedname.c_str() << "(void *obj,TCollection *coll,TFileMergeInfo *info) {" << "\n" << "      return ((" << classname.c_str() << "*)obj)->Merge(coll,info);" << "\n" << "   }" << "\n";
   } else if (HasOldMerge(decl, interp)) {
      finalString << "   // Wrapper around the merge function." << "\n" << "   static Long64_t  merge_" << mappedname.c_str() << "(void *obj,TCollection *coll,TFileMergeInfo *) {" << "\n" << "      return ((" << classname.c_str() << "*)obj)->Merge(coll);" << "\n" << "   }" << "\n";
   }

   if (HasResetAfterMerge(decl, interp)) {
      finalString << "   // Wrapper around the Reset function." << "\n" << "   static void reset_" << mappedname.c_str() << "(void *obj,TFileMergeInfo *info) {" << "\n" << "      ((" << classname.c_str() << "*)obj)->ResetAfterMerge(info);" << "\n" << "   }" << "\n";
   }
   finalString << "} // end of namespace ROOT for class " << classname.c_str() << "\n" << "\n";
}

////////////////////////////////////////////////////////////////////////////////
/// Write interface function for STL members

void ROOT::TMetaUtils::WritePointersSTL(const AnnotatedRecordDecl &cl,
                                        const cling::Interpreter &interp,
                                        const TNormalizedCtxt &normCtxt)
{
   std::string a;
   std::string clName;
   TMetaUtils::GetCppName(clName, ROOT::TMetaUtils::GetFileName(*cl.GetRecordDecl(), interp).str().c_str());
   int version = ROOT::TMetaUtils::GetClassVersion(cl.GetRecordDecl(),interp);
   if (version == 0) return;
   if (version < 0 && !(cl.RequestStreamerInfo()) ) return;


   const clang::CXXRecordDecl *clxx = llvm::dyn_cast<clang::CXXRecordDecl>(cl.GetRecordDecl());
   if (clxx == 0) return;

   // We also need to look at the base classes.
   for(clang::CXXRecordDecl::base_class_const_iterator iter = clxx->bases_begin(), end = clxx->bases_end();
       iter != end;
       ++iter)
   {
      int k = ROOT::TMetaUtils::IsSTLContainer(*iter);
      if (k!=0) {
         Internal::RStl::Instance().GenerateTClassFor( iter->getType(), interp, normCtxt);
      }
   }

   // Loop over the non static data member.
   for(clang::RecordDecl::field_iterator field_iter = clxx->field_begin(), end = clxx->field_end();
       field_iter != end;
       ++field_iter)
   {
      std::string mTypename;
      ROOT::TMetaUtils::GetQualifiedName(mTypename, field_iter->getType(), *clxx);

      //member is a string
      {
         const char*shortTypeName = ROOT::TMetaUtils::ShortTypeName(mTypename.c_str());
         if (!strcmp(shortTypeName, "string")) {
            continue;
         }
      }

      if (!ROOT::TMetaUtils::IsStreamableObject(**field_iter, interp)) continue;

      int k = ROOT::TMetaUtils::IsSTLContainer( **field_iter );
      if (k!=0) {
         //          fprintf(stderr,"Add %s which is also",m.Type()->Name());
         //          fprintf(stderr," %s\n",R__TrueName(**field_iter) );
         clang::QualType utype(ROOT::TMetaUtils::GetUnderlyingType(field_iter->getType()),0);
         Internal::RStl::Instance().GenerateTClassFor(utype, interp, normCtxt);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// TrueName strips the typedefs and array dimensions.

std::string ROOT::TMetaUtils::TrueName(const clang::FieldDecl &m)
{
   const clang::Type *rawtype = m.getType()->getCanonicalTypeInternal().getTypePtr();
   if (rawtype->isArrayType()) {
      rawtype = rawtype->getBaseElementTypeUnsafe ();
   }

   std::string result;
   ROOT::TMetaUtils::GetQualifiedName(result, clang::QualType(rawtype,0), m);
   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the version number of the class or -1
/// if the function Class_Version does not exist.

int ROOT::TMetaUtils::GetClassVersion(const clang::RecordDecl *cl, const cling::Interpreter& interp)
{
   const clang::CXXRecordDecl* CRD = llvm::dyn_cast<clang::CXXRecordDecl>(cl);
   if (!CRD) {
      // Must be an enum or namespace.
      // FIXME: Make it work for a namespace!
      return -1;
   }
   const clang::FunctionDecl* funcCV = ROOT::TMetaUtils::ClassInfo__HasMethod(CRD,"Class_Version",interp);

   // if we have no Class_Info() return -1.
   if (!funcCV) return -1;

   // if we have many Class_Info() (?!) return 1.
   if (funcCV == (clang::FunctionDecl*)-1) return 1;

   return GetTrivialIntegralReturnValue(funcCV, interp).second;
}

////////////////////////////////////////////////////////////////////////////////
/// If the function contains 'just': return SomeValue;
/// this routine will extract this value and return it.
/// The first element is set to true we have the body of the function and it
/// is indeed a trivial function with just a return of a value.
/// The second element contains the value (or -1 is case of failure)

std::pair<bool, int>
ROOT::TMetaUtils::GetTrivialIntegralReturnValue(const clang::FunctionDecl *funcCV, const cling::Interpreter &interp)
{
   using res_t = std::pair<bool, int>;

   const clang::CompoundStmt* FuncBody
      = llvm::dyn_cast_or_null<clang::CompoundStmt>(funcCV->getBody());
   if (!FuncBody)
      return res_t{false, -1};
   if (FuncBody->size() != 1) {
      // This is a non-ClassDef(), complex function - it might depend on state
      // and thus we'll need the runtime and cannot determine the result
      // statically.
      return res_t{false, -1};
   }
   const clang::ReturnStmt* RetStmt
      = llvm::dyn_cast<clang::ReturnStmt>(FuncBody->body_back());
   if (!RetStmt)
      return res_t{false, -1};
   const clang::Expr* RetExpr = RetStmt->getRetValue();
   // ClassDef controls the content of Class_Version() but not the return
   // expression which is CPP expanded from what the user provided as second
   // ClassDef argument. It's usually just be an integer literal but it could
   // also be an enum or a variable template for all we know.
   // Go through ICE to be more general.
   llvm::APSInt RetRes;
   if (!RetExpr->isIntegerConstantExpr(RetRes, funcCV->getASTContext()))
      return res_t{false, -1};
   if (RetRes.isSigned()) {
      return res_t{true, (Version_t)RetRes.getSExtValue()};
   }
   // else
   return res_t{true, (Version_t)RetRes.getZExtValue()};
}

////////////////////////////////////////////////////////////////////////////////
/// Is this an STL container.

int ROOT::TMetaUtils::IsSTLContainer(const ROOT::TMetaUtils::AnnotatedRecordDecl &annotated)
{
   return TMetaUtils::IsSTLCont(*annotated.GetRecordDecl());
}

////////////////////////////////////////////////////////////////////////////////
/// Is this an STL container?

ROOT::ESTLType ROOT::TMetaUtils::IsSTLContainer(const clang::FieldDecl &m)
{
   clang::QualType type = m.getType();
   clang::RecordDecl *decl = ROOT::TMetaUtils::GetUnderlyingRecordDecl(type);

   if (decl) return TMetaUtils::IsSTLCont(*decl);
   else return ROOT::kNotSTL;
}

////////////////////////////////////////////////////////////////////////////////
/// Is this an STL container?

int ROOT::TMetaUtils::IsSTLContainer(const clang::CXXBaseSpecifier &base)
{
   clang::QualType type = base.getType();
   clang::RecordDecl *decl = ROOT::TMetaUtils::GetUnderlyingRecordDecl(type);

   if (decl) return TMetaUtils::IsSTLCont(*decl);
   else return ROOT::kNotSTL;
}

////////////////////////////////////////////////////////////////////////////////
/// Calls the given lambda on every header in the given module.
/// includeDirectlyUsedModules designates if the foreach should also loop over
/// the headers in all modules that are directly used via a `use` declaration
/// in the modulemap.
void ROOT::TMetaUtils::foreachHeaderInModule(const clang::Module &module,
                                             const std::function<void(const clang::Module::Header &)> &closure,
                                             bool includeDirectlyUsedModules)
{
   // Iterates over all headers in a module and calls the closure on each.

   // FIXME: We currently have to hardcode '4' to do this. Maybe we
   // will have a nicer way to do this in the future.
   // NOTE: This is on purpose '4', not '5' which is the size of the
   // vector. The last element is the list of excluded headers which we
   // obviously don't want to check here.
   const std::size_t publicHeaderIndex = 4;

   // Integrity check in case this array changes its size at some point.
   const std::size_t maxArrayLength = ((sizeof module.Headers) / (sizeof *module.Headers));
   static_assert(publicHeaderIndex + 1 == maxArrayLength,
                 "'Headers' has changed it's size, we need to update publicHeaderIndex");

   // Make a list of modules and submodules that we can check for headers.
   // We use a SetVector to prevent an infinite loop in unlikely case the
   // modules somehow are messed up and don't form a tree...
   llvm::SetVector<const clang::Module *> modules;
   modules.insert(&module);
   for (size_t i = 0; i < modules.size(); ++i) {
      const clang::Module *M = modules[i];
      for (const clang::Module *subModule : M->submodules())
         modules.insert(subModule);
   }

   for (const clang::Module *m : modules) {
      if (includeDirectlyUsedModules) {
         for (clang::Module *used : m->DirectUses) {
            foreachHeaderInModule(*used, closure, true);
         }
      }

      for (std::size_t i = 0; i < publicHeaderIndex; i++) {
         auto &headerList = m->Headers[i];
         for (const clang::Module::Header &moduleHeader : headerList) {
            closure(moduleHeader);
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return the absolute type of typeDesc.
/// E.g.: typeDesc = "class TNamed**", returns "TNamed".
/// we remove * and const keywords. (we do not want to remove & ).
/// You need to use the result immediately before it is being overwritten.

const char *ROOT::TMetaUtils::ShortTypeName(const char *typeDesc)
{
   static char t[4096];
   static const char* constwd = "const ";
   static const char* constwdend = "const";

   const char *s;
   char *p=t;
   int lev=0;
   for (s=typeDesc;*s;s++) {
      if (*s=='<') lev++;
      if (*s=='>') lev--;
      if (lev==0 && *s=='*') continue;
      if (lev==0 && (strncmp(constwd,s,strlen(constwd))==0
                     ||strcmp(constwdend,s)==0 ) ) {
         s+=strlen(constwd)-1; // -1 because the loop adds 1
         continue;
      }
      if (lev==0 && *s==' ' && *(s+1)!='*') { p = t; continue;}
      if (p - t > (long)sizeof(t)) {
         printf("ERROR (rootcling): type name too long for StortTypeName: %s\n",
                typeDesc);
         p[0] = 0;
         return t;
      }
      *p++ = *s;
   }
   p[0]=0;

   return t;
}

bool ROOT::TMetaUtils::IsStreamableObject(const clang::FieldDecl &m,
                                          const cling::Interpreter& interp)
{
   const char *comment = ROOT::TMetaUtils::GetComment( m ).data();

   // Transient
   if (comment[0] == '!') return false;

   clang::QualType type = m.getType();

   if (type->isReferenceType()) {
      // Reference can not be streamed.
      return false;
   }

   std::string mTypeName = type.getAsString(m.getASTContext().getPrintingPolicy());
   if (!strcmp(mTypeName.c_str(), "string") || !strcmp(mTypeName.c_str(), "string*")) {
      return true;
   }
   if (!strcmp(mTypeName.c_str(), "std::string") || !strcmp(mTypeName.c_str(), "std::string*")) {
      return true;
   }

   if (ROOT::TMetaUtils::IsSTLContainer(m)) {
      return true;
   }

   const clang::Type *rawtype = type.getTypePtr()->getBaseElementTypeUnsafe ();

   if (rawtype->isPointerType()) {
      //Get to the 'raw' type.
      clang::QualType pointee;
      while ( (pointee = rawtype->getPointeeType()) , pointee.getTypePtrOrNull() && pointee.getTypePtr() != rawtype)
      {
        rawtype = pointee.getTypePtr();
      }
   }

   if (rawtype->isFundamentalType() || rawtype->isEnumeralType()) {
      // not an ojbect.
      return false;
   }

   const clang::CXXRecordDecl *cxxdecl = rawtype->getAsCXXRecordDecl();
   if (cxxdecl && ROOT::TMetaUtils::ClassInfo__HasMethod(cxxdecl,"Streamer", interp)) {
      if (!(ROOT::TMetaUtils::ClassInfo__HasMethod(cxxdecl,"Class_Version", interp))) return true;
      int version = ROOT::TMetaUtils::GetClassVersion(cxxdecl,interp);
      if (version > 0) return true;
   }
   return false;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the absolute type of typeDesc.
/// E.g.: typeDesc = "class TNamed**", returns "TNamed".
/// we remove * and const keywords. (we do not want to remove & ).
/// You need to use the result immediately before it is being overwritten.

std::string ROOT::TMetaUtils::ShortTypeName(const clang::FieldDecl &m)
{
   const clang::Type *rawtype = m.getType().getTypePtr();

   //Get to the 'raw' type.
   clang::QualType pointee;
   while ( rawtype->isPointerType() && ((pointee = rawtype->getPointeeType()) , pointee.getTypePtrOrNull()) && pointee.getTypePtr() != rawtype)
   {
      rawtype = pointee.getTypePtr();
   }

   std::string result;
   ROOT::TMetaUtils::GetQualifiedName(result, clang::QualType(rawtype,0), m);
   return result;
}

////////////////////////////////////////////////////////////////////////////////

clang::RecordDecl *ROOT::TMetaUtils::GetUnderlyingRecordDecl(clang::QualType type)
{
   const clang::Type *rawtype = ROOT::TMetaUtils::GetUnderlyingType(type);

   if (rawtype->isFundamentalType() || rawtype->isEnumeralType()) {
      // not an object.
      return 0;
   }
   return rawtype->getAsCXXRecordDecl();
}

////////////////////////////////////////////////////////////////////////////////
/// Generate the code of the class
/// If the requestor is genreflex, request the new streamer format

void ROOT::TMetaUtils::WriteClassCode(CallWriteStreamer_t WriteStreamerFunc,
                                      const AnnotatedRecordDecl &cl,
                                      const cling::Interpreter &interp,
                                      const TNormalizedCtxt &normCtxt,
                                      std::ostream& dictStream,
                                      const RConstructorTypes& ctorTypes,
                                      bool isGenreflex=false)
{
   const clang::CXXRecordDecl* decl = llvm::dyn_cast<clang::CXXRecordDecl>(cl.GetRecordDecl());

   if (!decl  || !decl->isCompleteDefinition()) {
      return;
   }

   std::string fullname;
   ROOT::TMetaUtils::GetQualifiedName(fullname,cl);
   if (TClassEdit::IsSTLCont(fullname) ) {
      Internal::RStl::Instance().GenerateTClassFor(cl.GetNormalizedName(), llvm::dyn_cast<clang::CXXRecordDecl>(cl.GetRecordDecl()), interp, normCtxt);
      return;
   }

   if (ROOT::TMetaUtils::ClassInfo__HasMethod(cl,"Streamer",interp)) {
      // The !genreflex is there to prevent genreflex to select collections which are data members
      // This is to maintain the behaviour of ROOT5 and ROOT6 up to 6.07 included.
      if (cl.RootFlag() && !isGenreflex) ROOT::TMetaUtils::WritePointersSTL(cl, interp, normCtxt); // In particular this detect if the class has a version number.
      if (!(cl.RequestNoStreamer())) {
         (*WriteStreamerFunc)(cl, interp, normCtxt, dictStream, isGenreflex || cl.RequestStreamerInfo());
      } else
         ROOT::TMetaUtils::Info(0, "Class %s: Do not generate Streamer() [*** custom streamer ***]\n",fullname.c_str());
   } else {
      ROOT::TMetaUtils::Info(0, "Class %s: Streamer() not declared\n", fullname.c_str());

      // See comment above about the !isGenreflex
      if (cl.RequestStreamerInfo() && !isGenreflex) ROOT::TMetaUtils::WritePointersSTL(cl, interp, normCtxt);
   }
   ROOT::TMetaUtils::WriteAuxFunctions(dictStream, cl, decl, interp, ctorTypes, normCtxt);
}

////////////////////////////////////////////////////////////////////////////////
/// Add any unspecified template parameters to the class template instance,
/// mentioned anywhere in the type.
///
/// Note: this does not strip any typedef but could be merged with cling::utils::Transform::GetPartiallyDesugaredType
/// if we can safely replace TClassEdit::IsStd with a test on the declaring scope
/// and if we can resolve the fact that the added parameter do not take into account possible use/dependences on Double32_t
/// and if we decide that adding the default is the right long term solution or not.
/// Whether it is or not depend on the I/O on whether the default template argument might change or not
/// and whether they (should) affect the on disk layout (for STL containers, we do know they do not).

clang::QualType ROOT::TMetaUtils::AddDefaultParameters(clang::QualType instanceType,
                                                       const cling::Interpreter &interpreter,
                                                       const ROOT::TMetaUtils::TNormalizedCtxt &normCtxt)
{
   const clang::ASTContext& Ctx = interpreter.getCI()->getASTContext();

   clang::QualType originalType = instanceType;

   // In case of name* we need to strip the pointer first, add the default and attach
   // the pointer once again.
   if (llvm::isa<clang::PointerType>(instanceType.getTypePtr())) {
      // Get the qualifiers.
      clang::Qualifiers quals = instanceType.getQualifiers();
      clang::QualType newPointee = AddDefaultParameters(instanceType->getPointeeType(), interpreter, normCtxt);
      if (newPointee != instanceType->getPointeeType()) {
         instanceType = Ctx.getPointerType(newPointee);
         // Add back the qualifiers.
         instanceType = Ctx.getQualifiedType(instanceType, quals);
      }
      return instanceType;
   }

   // In case of Int_t& we need to strip the pointer first, desugar and attach
   // the pointer once again.
   if (llvm::isa<clang::ReferenceType>(instanceType.getTypePtr())) {
      // Get the qualifiers.
      bool isLValueRefTy = llvm::isa<clang::LValueReferenceType>(instanceType.getTypePtr());
      clang::Qualifiers quals = instanceType.getQualifiers();
      clang::QualType newPointee = AddDefaultParameters(instanceType->getPointeeType(), interpreter, normCtxt);

      if (newPointee != instanceType->getPointeeType()) {
         // Add the r- or l- value reference type back to the desugared one
         if (isLValueRefTy)
            instanceType = Ctx.getLValueReferenceType(newPointee);
         else
            instanceType = Ctx.getRValueReferenceType(newPointee);
         // Add back the qualifiers.
         instanceType = Ctx.getQualifiedType(instanceType, quals);
      }
      return instanceType;
   }

   // Treat the Scope.
   bool prefix_changed = false;
   clang::NestedNameSpecifier* prefix = 0;
   clang::Qualifiers prefix_qualifiers = instanceType.getLocalQualifiers();
   const clang::ElaboratedType* etype
      = llvm::dyn_cast<clang::ElaboratedType>(instanceType.getTypePtr());
   if (etype) {
      // We have to also handle the prefix.
      prefix = AddDefaultParametersNNS(Ctx, etype->getQualifier(), interpreter, normCtxt);
      prefix_changed = prefix != etype->getQualifier();
      instanceType = clang::QualType(etype->getNamedType().getTypePtr(),0);
   }

   // In case of template specializations iterate over the arguments and
   // add unspecified default parameter.

   const clang::TemplateSpecializationType* TST
      = llvm::dyn_cast<const clang::TemplateSpecializationType>(instanceType.getTypePtr());

   const clang::ClassTemplateSpecializationDecl* TSTdecl
      = llvm::dyn_cast_or_null<const clang::ClassTemplateSpecializationDecl>(instanceType.getTypePtr()->getAsCXXRecordDecl());

   // Don't add the default paramater onto std classes.
   // We really need this for __shared_ptr which add a enum constant value which
   // is spelled in its 'numeral' form and thus the resulting type name is
   // incorrect.  We also can used this for any of the STL collections where we
   // know we don't want the default argument.   For the other members of the
   // std namespace this is dubious (because TMetaUtils::GetNormalizedName would
   // not drop those defaults).  [I.e. the real test ought to be is std and
   // name is __shared_ptr or vector or list or set or etc.]
   bool isStdDropDefault = TSTdecl && IsStdDropDefaultClass(*TSTdecl);

   bool mightHaveChanged = false;
   if (TST && TSTdecl) {

      clang::Sema& S = interpreter.getCI()->getSema();
      clang::TemplateDecl *Template = TSTdecl->getSpecializedTemplate()->getMostRecentDecl();
      clang::TemplateParameterList *Params = Template->getTemplateParameters();
      clang::TemplateParameterList::iterator Param = Params->begin(); // , ParamEnd = Params->end();
      //llvm::SmallVectorImpl<TemplateArgument> Converted; // Need to contains the other arguments.
      // Converted seems to be the same as our 'desArgs'

      unsigned int dropDefault = normCtxt.GetConfig().DropDefaultArg(*Template);

      llvm::SmallVector<clang::TemplateArgument, 4> desArgs;
      unsigned int Idecl = 0, Edecl = TSTdecl->getTemplateArgs().size();
      unsigned int maxAddArg = TSTdecl->getTemplateArgs().size() - dropDefault;
      for(clang::TemplateSpecializationType::iterator
             I = TST->begin(), E = TST->end();
          Idecl != Edecl;
          I!=E ? ++I : 0, ++Idecl, ++Param) {

         if (I != E) {

            if (I->getKind() == clang::TemplateArgument::Template) {
               clang::TemplateName templateName = I->getAsTemplate();
               clang::TemplateDecl* templateDecl = templateName.getAsTemplateDecl();
               if (templateDecl) {
                  clang::DeclContext* declCtxt = templateDecl->getDeclContext();

                  if (declCtxt && !templateName.getAsQualifiedTemplateName()){
                     clang::NamespaceDecl* ns = clang::dyn_cast<clang::NamespaceDecl>(declCtxt);
                     clang::NestedNameSpecifier* nns;
                     if (ns) {
                        nns = cling::utils::TypeName::CreateNestedNameSpecifier(Ctx, ns);
                     } else if (clang::TagDecl* TD = llvm::dyn_cast<clang::TagDecl>(declCtxt)) {
                        nns = cling::utils::TypeName::CreateNestedNameSpecifier(Ctx,TD, false /*FullyQualified*/);
                     } else {
                        // TU scope
                        desArgs.push_back(*I);
                        continue;
                     }
                     clang::TemplateName templateNameWithNSS ( Ctx.getQualifiedTemplateName(nns, false, templateDecl) );
                     desArgs.push_back(clang::TemplateArgument(templateNameWithNSS));
                     mightHaveChanged = true;
                     continue;
                  }
               }
            }

            if (I->getKind() != clang::TemplateArgument::Type) {
               desArgs.push_back(*I);
               continue;
            }

            clang::QualType SubTy = I->getAsType();

            // Check if the type needs more desugaring and recurse.
            // (Originally this was limited to elaborated and templated type,
            // but we also need to do it for pointer and reference type
            // and who knows what, so do it always)
            clang::QualType newSubTy = AddDefaultParameters(SubTy,
                                                            interpreter,
                                                            normCtxt);
            if (SubTy != newSubTy) {
               mightHaveChanged = true;
               desArgs.push_back(clang::TemplateArgument(newSubTy));
            } else {
               desArgs.push_back(*I);
            }
            // Converted.push_back(TemplateArgument(ArgTypeForTemplate));
         } else if (!isStdDropDefault && Idecl < maxAddArg) {

            mightHaveChanged = true;

            const clang::TemplateArgument& templateArg
               = TSTdecl->getTemplateArgs().get(Idecl);
            if (templateArg.getKind() != clang::TemplateArgument::Type) {
               desArgs.push_back(templateArg);
               continue;
            }
            clang::QualType SubTy = templateArg.getAsType();

            clang::SourceLocation TemplateLoc = Template->getSourceRange ().getBegin(); //NOTE: not sure that this is the 'right' location.
            clang::SourceLocation RAngleLoc = TSTdecl->getSourceRange().getBegin(); // NOTE: most likely wrong, I think this is expecting the location of right angle

            clang::TemplateTypeParmDecl *TTP = llvm::dyn_cast<clang::TemplateTypeParmDecl>(*Param);
            {
               // We may induce template instantiation
               cling::Interpreter::PushTransactionRAII clingRAII(const_cast<cling::Interpreter*>(&interpreter));
               clang::sema::HackForDefaultTemplateArg raii;
               bool HasDefaultArgs;
               clang::TemplateArgumentLoc ArgType = S.SubstDefaultTemplateArgumentIfAvailable(
                                                                                              Template,
                                                                                              TemplateLoc,
                                                                                              RAngleLoc,
                                                                                              TTP,
                                                                                              desArgs,
                                                                                              HasDefaultArgs);
               // The substition can fail, in which case there would have been compilation
               // error printed on the screen.
               if (ArgType.getArgument().isNull()
                   || ArgType.getArgument().getKind() != clang::TemplateArgument::Type) {
                  ROOT::TMetaUtils::Error("ROOT::TMetaUtils::AddDefaultParameters",
                                          "Template parameter substitution failed for %s around %s\n",
                                          instanceType.getAsString().c_str(), SubTy.getAsString().c_str());
                  break;
               }
               clang::QualType BetterSubTy = ArgType.getArgument().getAsType();
               SubTy = cling::utils::Transform::GetPartiallyDesugaredType(Ctx,BetterSubTy,normCtxt.GetConfig(),/*fullyQualified=*/ true);
            }
            SubTy = AddDefaultParameters(SubTy,interpreter,normCtxt);
            desArgs.push_back(clang::TemplateArgument(SubTy));
         } else {
            // We are past the end of the list of specified arguements and we
            // do not want to add the default, no need to continue.
            break;
         }
      }

      // If we added default parameter, allocate new type in the AST.
      if (mightHaveChanged) {
         instanceType = Ctx.getTemplateSpecializationType(TST->getTemplateName(),
                                                          desArgs,
                                                          TST->getCanonicalTypeInternal());
      }
   }

   if (!prefix_changed && !mightHaveChanged) return originalType;
   if (prefix) {
      instanceType = Ctx.getElaboratedType(clang::ETK_None,prefix,instanceType);
      instanceType = Ctx.getQualifiedType(instanceType,prefix_qualifiers);
   }
   return instanceType;
}

////////////////////////////////////////////////////////////////////////////////
/// ValidArrayIndex return a static string (so use it or copy it immediatly, do not
/// call GrabIndex twice in the same expression) containing the size of the
/// array data member.
/// In case of error, or if the size is not specified, GrabIndex returns 0.
/// If errnum is not null, *errnum updated with the error number:
///   Cint::G__DataMemberInfo::G__VALID     : valid array index
///   Cint::G__DataMemberInfo::G__NOT_INT   : array index is not an int
///   Cint::G__DataMemberInfo::G__NOT_DEF   : index not defined before array
///                                          (this IS an error for streaming to disk)
///   Cint::G__DataMemberInfo::G__IS_PRIVATE: index exist in a parent class but is private
///   Cint::G__DataMemberInfo::G__UNKNOWN   : index is not known
/// If errstr is not null, *errstr is updated with the address of a static
///   string containing the part of the index with is invalid.

llvm::StringRef ROOT::TMetaUtils::DataMemberInfo__ValidArrayIndex(const clang::DeclaratorDecl &m, int *errnum, llvm::StringRef *errstr)
{
   llvm::StringRef title;

   // Try to get the comment either from the annotation or the header file if present
   if (clang::AnnotateAttr *A = m.getAttr<clang::AnnotateAttr>())
      title = A->getAnnotation();
   else
      // Try to get the comment from the header file if present
      title = ROOT::TMetaUtils::GetComment( m );

   // Let's see if the user provided us with some information
   // with the format: //[dimension] this is the dim of the array
   // dimension can be an arithmetical expression containing, literal integer,
   // the operator *,+ and - and data member of integral type.  In addition the
   // data members used for the size of the array need to be defined prior to
   // the array.

   if (errnum) *errnum = VALID;

   if (title.size() == 0 || (title[0] != '[')) return llvm::StringRef();
   size_t rightbracket = title.find(']');
   if (rightbracket == llvm::StringRef::npos) return llvm::StringRef();

   std::string working;
   llvm::StringRef indexvar(title.data()+1,rightbracket-1);

   // now we should have indexvar=dimension
   // Let's see if this is legal.
   // which means a combination of data member and digit separated by '*','+','-'
   // First we remove white spaces.
   unsigned int i;
   size_t indexvarlen = indexvar.size();
   for ( i=0; i<indexvarlen; i++) {
      if (!isspace(indexvar[i])) {
         working += indexvar[i];
      }
   };

   // Now we go through all indentifiers
   const char *tokenlist = "*+-";
   char *current = const_cast<char*>(working.c_str());
   current = strtok(current,tokenlist); // this method does not need to be reentrant

   while (current!=0) {
      // Check the token
      if (isdigit(current[0])) {
         for(i=0;i<strlen(current);i++) {
            if (!isdigit(current[i])) {
               // Error we only access integer.
               //NOTE: *** Need to print an error;
               //fprintf(stderr,"*** Datamember %s::%s: size of array (%s) is not an interger\n",
               //        member.MemberOf()->Name(), member.Name(), current);
               if (errstr) *errstr = current;
               if (errnum) *errnum = NOT_INT;
               return llvm::StringRef();
            }
         }
      } else { // current token is not a digit
         // first let's see if it is a data member:
         int found = 0;
         const clang::CXXRecordDecl *parent_clxx = llvm::dyn_cast<clang::CXXRecordDecl>(m.getDeclContext());
         const clang::FieldDecl *index1 = 0;
         if (parent_clxx)
            index1 = GetDataMemberFromAll(*parent_clxx, current );
         if ( index1 ) {
            if ( IsFieldDeclInt(index1) ) {
               found = 1;
               // Let's see if it has already been written down in the
               // Streamer.
               // Let's see if we already wrote it down in the
               // streamer.
               for(clang::RecordDecl::field_iterator field_iter = parent_clxx->field_begin(), end = parent_clxx->field_end();
                   field_iter != end;
                   ++field_iter)
               {
                  if ( field_iter->getNameAsString() == m.getNameAsString() ) {
                     // we reached the current data member before
                     // reaching the index so we have not written it yet!
                     //NOTE: *** Need to print an error;
                     //fprintf(stderr,"*** Datamember %s::%s: size of array (%s) has not been defined before the array \n",
                     //        member.MemberOf()->Name(), member.Name(), current);
                     if (errstr) *errstr = current;
                     if (errnum) *errnum = NOT_DEF;
                     return llvm::StringRef();
                  }
                  if ( field_iter->getNameAsString() == index1->getNameAsString() ) {
                     break;
                  }
               } // end of while (m_local.Next())
            } else {
               //NOTE: *** Need to print an error;
               //fprintf(stderr,"*** Datamember %s::%s: size of array (%s) is not int \n",
               //        member.MemberOf()->Name(), member.Name(), current);
               if (errstr) *errstr = current;
               if (errnum) *errnum = NOT_INT;
               return llvm::StringRef();
            }
         } else {
            // There is no variable by this name in this class, let see
            // the base classes!:
            index1 = GetDataMemberFromAllParents( *parent_clxx, current );
            if ( index1 ) {
               if ( IsFieldDeclInt(index1) ) {
                  found = 1;
               } else {
                  // We found a data member but it is the wrong type
                  //NOTE: *** Need to print an error;
                  //fprintf(stderr,"*** Datamember %s::%s: size of array (%s) is not int \n",
                  //  member.MemberOf()->Name(), member.Name(), current);
                  if (errnum) *errnum = NOT_INT;
                  if (errstr) *errstr = current;
                  //NOTE: *** Need to print an error;
                  //fprintf(stderr,"*** Datamember %s::%s: size of array (%s) is not int \n",
                  //  member.MemberOf()->Name(), member.Name(), current);
                  if (errnum) *errnum = NOT_INT;
                  if (errstr) *errstr = current;
                  return llvm::StringRef();
               }
               if ( found && (index1->getAccess() == clang::AS_private) ) {
                  //NOTE: *** Need to print an error;
                  //fprintf(stderr,"*** Datamember %s::%s: size of array (%s) is a private member of %s \n",
                  if (errstr) *errstr = current;
                  if (errnum) *errnum = IS_PRIVATE;
                  return llvm::StringRef();
               }
            }
            if (!found) {
               //NOTE: *** Need to print an error;
               //fprintf(stderr,"*** Datamember %s::%s: size of array (%s) is not known \n",
               //        member.MemberOf()->Name(), member.Name(), indexvar);
               if (errstr) *errstr = indexvar;
               if (errnum) *errnum = UNKNOWN;
               return llvm::StringRef();
            } // end of if not found
         } // end of if is a data member of the class
      } // end of if isdigit

      current = strtok(0,tokenlist);
   } // end of while loop on tokens

   return indexvar;

}

////////////////////////////////////////////////////////////////////////////////
/// Return (in the argument 'output') a mangled version of the C++ symbol/type (pass as 'input')
/// that can be used in C++ as a variable name.

void ROOT::TMetaUtils::GetCppName(std::string &out, const char *in)
{
   out.resize(strlen(in)*2);
   unsigned int i=0,j=0,c;
   while((c=in[i])) {
      if (out.capacity() < (j+3)) {
         out.resize(2*j+3);
      }
      switch(c) { // We resized the underlying buffer if needed
         case '+': strcpy(const_cast<char*>(out.data())+j,"pL"); j+=2; break;
         case '-': strcpy(const_cast<char*>(out.data())+j,"mI"); j+=2; break;
         case '*': strcpy(const_cast<char*>(out.data())+j,"mU"); j+=2; break;
         case '/': strcpy(const_cast<char*>(out.data())+j,"dI"); j+=2; break;
         case '&': strcpy(const_cast<char*>(out.data())+j,"aN"); j+=2; break;
         case '%': strcpy(const_cast<char*>(out.data())+j,"pE"); j+=2; break;
         case '|': strcpy(const_cast<char*>(out.data())+j,"oR"); j+=2; break;
         case '^': strcpy(const_cast<char*>(out.data())+j,"hA"); j+=2; break;
         case '>': strcpy(const_cast<char*>(out.data())+j,"gR"); j+=2; break;
         case '<': strcpy(const_cast<char*>(out.data())+j,"lE"); j+=2; break;
         case '=': strcpy(const_cast<char*>(out.data())+j,"eQ"); j+=2; break;
         case '~': strcpy(const_cast<char*>(out.data())+j,"wA"); j+=2; break;
         case '.': strcpy(const_cast<char*>(out.data())+j,"dO"); j+=2; break;
         case '(': strcpy(const_cast<char*>(out.data())+j,"oP"); j+=2; break;
         case ')': strcpy(const_cast<char*>(out.data())+j,"cP"); j+=2; break;
         case '[': strcpy(const_cast<char*>(out.data())+j,"oB"); j+=2; break;
         case ']': strcpy(const_cast<char*>(out.data())+j,"cB"); j+=2; break;
         case '!': strcpy(const_cast<char*>(out.data())+j,"nO"); j+=2; break;
         case ',': strcpy(const_cast<char*>(out.data())+j,"cO"); j+=2; break;
         case '$': strcpy(const_cast<char*>(out.data())+j,"dA"); j+=2; break;
         case ' ': strcpy(const_cast<char*>(out.data())+j,"sP"); j+=2; break;
         case ':': strcpy(const_cast<char*>(out.data())+j,"cL"); j+=2; break;
         case '"': strcpy(const_cast<char*>(out.data())+j,"dQ"); j+=2; break;
         case '@': strcpy(const_cast<char*>(out.data())+j,"aT"); j+=2; break;
         case '\'': strcpy(const_cast<char*>(out.data())+j,"sQ"); j+=2; break;
         case '\\': strcpy(const_cast<char*>(out.data())+j,"fI"); j+=2; break;
         default: out[j++]=c; break;
      }
      ++i;
   }
   out.resize(j);

   // Remove initial numbers if any
   std::size_t firstNonNumber = out.find_first_not_of("0123456789");
   out.replace(0,firstNonNumber,"");

   return;
}

static clang::SourceLocation
getFinalSpellingLoc(clang::SourceManager& sourceManager,
                    clang::SourceLocation sourceLoc) {
   // Follow macro expansion until we hit a source file.
   if (!sourceLoc.isFileID()) {
      return sourceManager.getExpansionRange(sourceLoc).second;
   }
   return sourceLoc;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the header file to be included to declare the Decl.

llvm::StringRef ROOT::TMetaUtils::GetFileName(const clang::Decl& decl,
                                              const cling::Interpreter& interp)
{
   // It looks like the template specialization decl actually contains _less_ information
   // on the location of the code than the decl (in case where there is forward declaration,
   // that is what the specialization points to).
   //
   // const clang::CXXRecordDecl* clxx = llvm::dyn_cast<clang::CXXRecordDecl>(decl);
   // if (clxx) {
   //    switch(clxx->getTemplateSpecializationKind()) {
   //       case clang::TSK_Undeclared:
   //          // We want the default behavior
   //          break;
   //       case clang::TSK_ExplicitInstantiationDeclaration:
   //       case clang::TSK_ExplicitInstantiationDefinition:
   //       case clang::TSK_ImplicitInstantiation: {
   //          // We want the location of the template declaration:
   //          const clang::ClassTemplateSpecializationDecl *tmplt_specialization = llvm::dyn_cast<clang::ClassTemplateSpecializationDecl> (clxx);
   //          if (tmplt_specialization) {
   //             // return GetFileName(const_cast< clang::ClassTemplateSpecializationDecl *>(tmplt_specialization)->getSpecializedTemplate());
   //          }
   //          break;
   //       }
   //       case clang::TSK_ExplicitSpecialization:
   //          // We want the default behavior
   //          break;
   //       default:
   //          break;
   //    }
   // }

   using namespace clang;
   SourceLocation headerLoc = decl.getLocation();

   static const char invalidFilename[] = "";
   if (!headerLoc.isValid()) return invalidFilename;

   HeaderSearch& HdrSearch = interp.getCI()->getPreprocessor().getHeaderSearchInfo();

   SourceManager& sourceManager = decl.getASTContext().getSourceManager();
   headerLoc = getFinalSpellingLoc(sourceManager, headerLoc);
   FileID headerFID = sourceManager.getFileID(headerLoc);
   SourceLocation includeLoc
      = getFinalSpellingLoc(sourceManager,
                            sourceManager.getIncludeLoc(headerFID));

   const FileEntry *headerFE = sourceManager.getFileEntryForID(headerFID);
   while (includeLoc.isValid() && sourceManager.isInSystemHeader(includeLoc)) {
      const DirectoryLookup *foundDir = 0;
      // use HeaderSearch on the basename, to make sure it takes a header from
      // the include path (e.g. not from /usr/include/bits/)
      assert(headerFE && "Couldn't find FileEntry from FID!");
      const FileEntry *FEhdr
         = HdrSearch.LookupFile(llvm::sys::path::filename(headerFE->getName()),
                                SourceLocation(),
                                true /*isAngled*/, 0/*FromDir*/, foundDir,
                                ArrayRef<std::pair<const FileEntry *, const DirectoryEntry *>>(),
                                0/*Searchpath*/, 0/*RelPath*/,
                                0/*IsMapped*/, 0/*RequestingModule*/, 0/*SuggestedModule*/,
                                false /*SkipCache*/,
                                false /*BuildSystemModule*/,
                                false /*OpenFile*/, true /*CacheFailures*/);
      if (FEhdr) break;
      headerFID = sourceManager.getFileID(includeLoc);
      headerFE = sourceManager.getFileEntryForID(headerFID);
      // If we have a system header in a module we can't just trace back the
      // original include with the preprocessor. But it should be enough if
      // we trace it back to the top-level system header that includes this
      // declaration.
      if (interp.getCI()->getLangOpts().Modules && !headerFE) {
         assert(decl.isFirstDecl() && "Couldn't trace back include from a decl"
                                      " that is not from an AST file");
         assert(StringRef(includeLoc.printToString(sourceManager)).startswith("<module-includes>"));
         break;
      }
      includeLoc = getFinalSpellingLoc(sourceManager,
                                       sourceManager.getIncludeLoc(headerFID));
   }

   if (!headerFE) return invalidFilename;
   llvm::StringRef headerFileName = headerFE->getName();

   // Now headerFID references the last valid system header or the original
   // user file.
   // Find out how to include it by matching file name to include paths.
   // We assume that the file "/A/B/C/D.h" can at some level be included as
   // "C/D.h". Be we cannot know whether that happens to be a different file
   // with the same name. Thus we first find the longest stem that can be
   // reached, say B/C/D.h. Then we find the shortest one, say C/D.h, that
   // points to the same file as the long version. If such a short version
   // exists it will be returned. If it doesn't the long version is returned.
   bool isAbsolute = llvm::sys::path::is_absolute(headerFileName);
   const FileEntry* FELong = 0;
   // Find the longest available match.
   for (llvm::sys::path::const_iterator
           IDir = llvm::sys::path::begin(headerFileName),
           EDir = llvm::sys::path::end(headerFileName);
        !FELong && IDir != EDir; ++IDir) {
      if (isAbsolute) {
         // skip "/" part
         isAbsolute = false;
         continue;
      }
      size_t lenTrailing = headerFileName.size() - (IDir->data() - headerFileName.data());
      llvm::StringRef trailingPart(IDir->data(), lenTrailing);
      assert(trailingPart.data() + trailingPart.size()
             == headerFileName.data() + headerFileName.size()
             && "Mismatched partitioning of file name!");
      const DirectoryLookup* FoundDir = 0;
      FELong = HdrSearch.LookupFile(trailingPart, SourceLocation(),
                                    true /*isAngled*/, 0/*FromDir*/, FoundDir,
                                    ArrayRef<std::pair<const FileEntry *, const DirectoryEntry *>>(),
                                    0/*IsMapped*/, 0/*Searchpath*/, 0/*RelPath*/,
                                    0/*RequestingModule*/, 0/*SuggestedModule*/);
   }

   if (!FELong) {
      // We did not find any file part in any search path.
      return invalidFilename;
   }

   // Iterates through path *parts* "C"; we need trailing parts "C/D.h"
   for (llvm::sys::path::reverse_iterator
           IDir = llvm::sys::path::rbegin(headerFileName),
           EDir = llvm::sys::path::rend(headerFileName);
        IDir != EDir; ++IDir) {
      size_t lenTrailing = headerFileName.size() - (IDir->data() - headerFileName.data());
      llvm::StringRef trailingPart(IDir->data(), lenTrailing);
      assert(trailingPart.data() + trailingPart.size()
             == headerFileName.data() + headerFileName.size()
             && "Mismatched partitioning of file name!");
      const DirectoryLookup* FoundDir = 0;
      // Can we find it, and is it the same file as the long version?
      // (or are we back to the previously found spelling, which is fine, too)
      if (HdrSearch.LookupFile(trailingPart, SourceLocation(),
                               true /*isAngled*/, 0/*FromDir*/, FoundDir,
                               ArrayRef<std::pair<const FileEntry *, const DirectoryEntry *>>(),
                               0/*IsMapped*/,
                               0/*Searchpath*/,
                               0/*RelPath*/,
                               0/*RequestingModule*/, 0 /*SuggestedModule*/) == FELong) {
         return trailingPart;
      }
   }

   return invalidFilename;
}

////////////////////////////////////////////////////////////////////////////////

void ROOT::TMetaUtils::GetFullyQualifiedTypeName(std::string &typenamestr,
                                                 const clang::QualType &qtype,
                                                 const clang::ASTContext &astContext)
{
   std::string fqname = cling::utils::TypeName::GetFullyQualifiedName(qtype, astContext);
   TClassEdit::TSplitType splitname(fqname.c_str(),
                                    (TClassEdit::EModType)(TClassEdit::kLong64 | TClassEdit::kDropStd | TClassEdit::kDropStlDefault | TClassEdit::kKeepOuterConst));
   splitname.ShortType(typenamestr,TClassEdit::kDropStd | TClassEdit::kDropStlDefault | TClassEdit::kKeepOuterConst);
}

////////////////////////////////////////////////////////////////////////////////

void ROOT::TMetaUtils::GetFullyQualifiedTypeName(std::string &typenamestr,
                                                 const clang::QualType &qtype,
                                                 const cling::Interpreter &interpreter)
{
   // We need this because GetFullyQualifiedTypeName is triggering deserialization
   // This calling the same name function GetFullyQualifiedTypeName, but this should stay here because
   // callee doesn't have an interpreter pointer
   cling::Interpreter::PushTransactionRAII RAII(const_cast<cling::Interpreter*>(&interpreter));

   GetFullyQualifiedTypeName(typenamestr,
                             qtype,
                             interpreter.getCI()->getASTContext());
}

////////////////////////////////////////////////////////////////////////////////
/// Get the template specialisation decl and template decl behind the qualtype
/// Returns true if successfully found, false otherwise

bool ROOT::TMetaUtils::QualType2Template(const clang::QualType& qt,
                                         clang::ClassTemplateDecl*& ctd,
                                         clang::ClassTemplateSpecializationDecl*& ctsd)
{
   using namespace clang;
   const Type* theType = qt.getTypePtr();
   if (!theType){
      ctd=nullptr;
      ctsd=nullptr;
      return false;
   }

   if (theType->isPointerType()) {
      return QualType2Template(theType->getPointeeType(), ctd, ctsd);
   }

   if (const RecordType* rType = llvm::dyn_cast<RecordType>(theType)) {
      ctsd = llvm::dyn_cast_or_null<ClassTemplateSpecializationDecl>(rType->getDecl());
      if (ctsd) {
         ctd = ctsd->getSpecializedTemplate();
         return true;
      }
   }

   if (const SubstTemplateTypeParmType* sttpType = llvm::dyn_cast<SubstTemplateTypeParmType>(theType)){
      return QualType2Template(sttpType->getReplacementType(), ctd, ctsd);
   }


   ctsd = llvm::dyn_cast_or_null<ClassTemplateSpecializationDecl>(qt->getAsCXXRecordDecl());
   if(ctsd){
      ctd = ctsd->getSpecializedTemplate();
      return true;
   }

   ctd=nullptr;
   ctsd=nullptr;
   return false;
}

////////////////////////////////////////////////////////////////////////////////
/// Extract from a qualtype the class template if this makes sense.
/// Retuns the ClassTemplateDecl or nullptr otherwise.

clang::ClassTemplateDecl* ROOT::TMetaUtils::QualType2ClassTemplateDecl(const clang::QualType& qt)
{
   using namespace clang;
   ClassTemplateSpecializationDecl* ctsd;
   ClassTemplateDecl* ctd;
   QualType2Template(qt,ctd,ctsd);
   return ctd;
}

////////////////////////////////////////////////////////////////////////////////
/// These manipulations are necessary because a template specialisation type
/// does not inherit from a record type (there is an asymmetry between
/// the decls and the types in the clang interface).
/// We may need therefore to step into the "Decl dimension" to then get back
/// to the "type dimension".

clang::TemplateName ROOT::TMetaUtils::ExtractTemplateNameFromQualType(const clang::QualType& qt)
{
   using namespace clang;
   TemplateName theTemplateName;

   const Type* theType = qt.getTypePtr();

   if (const TemplateSpecializationType* tst = llvm::dyn_cast_or_null<const TemplateSpecializationType>(theType)) {
      theTemplateName = tst->getTemplateName();
   } // We step into the decl dimension
   else if (ClassTemplateDecl* ctd = QualType2ClassTemplateDecl(qt)) {
      theTemplateName = TemplateName(ctd);
   }

   return theTemplateName;
}

////////////////////////////////////////////////////////////////////////////////

static bool areEqualTypes(const clang::TemplateArgument& tArg,
                   llvm::SmallVectorImpl<clang::TemplateArgument>& preceedingTArgs,
                   const clang::NamedDecl& tPar,
                   const cling::Interpreter& interp,
                   const ROOT::TMetaUtils::TNormalizedCtxt& normCtxt)
{
   using namespace ROOT::TMetaUtils;
   using namespace clang;

   // Check if this is a type for security
   TemplateTypeParmDecl* ttpdPtr = const_cast<TemplateTypeParmDecl*>(llvm::dyn_cast<TemplateTypeParmDecl>(&tPar));
   if (!ttpdPtr) return false;
   if (!ttpdPtr->hasDefaultArgument()) return false; // we should not be here in this case, but we protect us.

   // Try the fast solution
   QualType tParQualType = ttpdPtr->getDefaultArgument();
   const QualType tArgQualType = tArg.getAsType();

   // Now the equality tests for non template specialisations.

   // The easy cases:
   // template <class T=double> class A; or
   // template <class T=A<float>> class B;
   if (tParQualType.getTypePtr() == tArgQualType.getTypePtr()) return true;

   // Here the difficulty comes. We have to check if the argument is equal to its
   // default. We can do that bootstrapping an argument which has the default value
   // based on the preceeding arguments.
   // Basically we ask sema to give us the value of the argument given the template
   // of behind the parameter and the all the arguments.
   // So:

   // Take the template out of the parameter

   const clang::ElaboratedType* etype
      = llvm::dyn_cast<clang::ElaboratedType>(tParQualType.getTypePtr());
   while (etype) {
      tParQualType = clang::QualType(etype->getNamedType().getTypePtr(),0);
      etype = llvm::dyn_cast<clang::ElaboratedType>(tParQualType.getTypePtr());
   }

   const TemplateSpecializationType* tst =
            llvm::dyn_cast<TemplateSpecializationType>(tParQualType.getTypePtr());

   if(!tst) // nothing more to be tried. They are different indeed.
      return false;

   ClassTemplateSpecializationDecl* TSTdecl
      = llvm::dyn_cast_or_null<ClassTemplateSpecializationDecl>(tArgQualType->getAsCXXRecordDecl());

   if(!TSTdecl) // nothing more to be tried. They are different indeed.
      return false;

   TemplateDecl *Template = tst->getTemplateName().getAsTemplateDecl();

   // Take the template location
   SourceLocation TemplateLoc = Template->getSourceRange ().getBegin();

   // Get the position of the "<" (LA) of the specializaion
   SourceLocation LAngleLoc = TSTdecl->getSourceRange().getBegin();


   // Enclose in a scope for the RAII
   bool isEqual=false;
   TemplateArgument newArg = tArg;
   {
      clang::Sema& S = interp.getCI()->getSema();
      cling::Interpreter::PushTransactionRAII clingRAII(const_cast<cling::Interpreter*>(&interp));
      clang::sema::HackForDefaultTemplateArg raii; // Hic sunt leones
      bool HasDefaultArgs;
      TemplateArgumentLoc defTArgLoc = S.SubstDefaultTemplateArgumentIfAvailable(Template,
                                                                                 TemplateLoc,
                                                                                 LAngleLoc,
                                                                                 ttpdPtr,
                                                                                 preceedingTArgs,
                                                                                 HasDefaultArgs);
      // The substition can fail, in which case there would have been compilation
      // error printed on the screen.
      newArg = defTArgLoc.getArgument();
      if (newArg.isNull() ||
          newArg.getKind() != clang::TemplateArgument::Type) {
         ROOT::TMetaUtils::Error("areEqualTypes",
                                 "Template parameter substitution failed!");
      }

      ClassTemplateSpecializationDecl* nTSTdecl
         = llvm::dyn_cast_or_null<ClassTemplateSpecializationDecl>(newArg.getAsType()->getAsCXXRecordDecl());
//         std::cout << "nSTdecl is " << nTSTdecl << std::endl;

      isEqual =  nTSTdecl->getMostRecentDecl() == TSTdecl->getMostRecentDecl() ||
                 tParQualType.getTypePtr() == newArg.getAsType().getTypePtr();
   }


   return isEqual;
}


////////////////////////////////////////////////////////////////////////////////
///   std::cout << "Are equal values?\n";

static bool areEqualValues(const clang::TemplateArgument& tArg,
                           const clang::NamedDecl& tPar)
{
   using namespace clang;
   const NonTypeTemplateParmDecl* nttpdPtr = llvm::dyn_cast<NonTypeTemplateParmDecl>(&tPar);
   if (!nttpdPtr) return false;
   const NonTypeTemplateParmDecl& nttpd = *nttpdPtr;

   if (!nttpd.hasDefaultArgument())
      return false;

   // 64 bits wide and signed (non unsigned, that is why "false")
   llvm::APSInt defaultValueAPSInt(64, false);
   if (Expr* defArgExpr = nttpd.getDefaultArgument()) {
      const ASTContext& astCtxt = nttpdPtr->getASTContext();
      defArgExpr->isIntegerConstantExpr(defaultValueAPSInt, astCtxt);
   }

   const int value = tArg.getAsIntegral().getLimitedValue();

   //   std::cout << (value == defaultValueAPSInt ? "yes!":"no")  << std::endl;
   return  value == defaultValueAPSInt;
}

////////////////////////////////////////////////////////////////////////////////
/// Check if this NamedDecl is a template parameter with a default argument.
/// This is a single interface to treat both integral and type parameters.
/// Returns true if this is the case, false otherwise

static bool isTypeWithDefault(const clang::NamedDecl* nDecl)
{
   using namespace clang;
   if (!nDecl) return false;
   if (const TemplateTypeParmDecl* ttpd = llvm::dyn_cast<TemplateTypeParmDecl>(nDecl))
      return ttpd->hasDefaultArgument();
   if (const NonTypeTemplateParmDecl* nttpd = llvm::dyn_cast<NonTypeTemplateParmDecl>(nDecl))
      return nttpd->hasDefaultArgument();
   return false;

}

static void KeepNParams(clang::QualType& normalizedType,
                        const clang::QualType& vanillaType,
                        const cling::Interpreter& interp,
                        const ROOT::TMetaUtils::TNormalizedCtxt& normCtxt);

// Returns true if normTArg might have changed.
static bool RecurseKeepNParams(clang::TemplateArgument &normTArg,
                               const clang::TemplateArgument &tArg,
                               const cling::Interpreter& interp,
                               const ROOT::TMetaUtils::TNormalizedCtxt& normCtxt,
                               const clang::ASTContext& astCtxt)
{
   using namespace ROOT::TMetaUtils;
   using namespace clang;

   // Once we know there is no more default parameter, we can run through to the end
   // and/or recurse in the template parameter packs.

   // If this is a type,
   // we need first of all to recurse: this argument may need to be manipulated
   if (tArg.getKind() == clang::TemplateArgument::Type) {
      QualType thisNormQualType = normTArg.getAsType();
      QualType thisArgQualType = tArg.getAsType();
      KeepNParams(thisNormQualType,
                  thisArgQualType,
                  interp,
                  normCtxt);
      normTArg = TemplateArgument(thisNormQualType);
      return (thisNormQualType != thisArgQualType);
   } else if (normTArg.getKind() == clang::TemplateArgument::Pack) {
      assert( tArg.getKind() == clang::TemplateArgument::Pack );

      SmallVector<TemplateArgument, 2> desArgs;
      bool mightHaveChanged = true;
      for (auto I = normTArg.pack_begin(), E = normTArg.pack_end(),
           FI = tArg.pack_begin(), FE = tArg.pack_end();
           I != E && FI != FE; ++I, ++FI)
      {
         TemplateArgument pack_arg(*I);
         mightHaveChanged |= RecurseKeepNParams(pack_arg, *FI, interp, normCtxt, astCtxt);
         desArgs.push_back(pack_arg);
      }
      if (mightHaveChanged) {
         ASTContext &mutableCtx( const_cast<ASTContext&>(astCtxt) );
         normTArg = TemplateArgument::CreatePackCopy(mutableCtx, desArgs);
      }
      return mightHaveChanged;
   }
   return false;
}


////////////////////////////////////////////////////////////////////////////////
/// This function allows to manipulate the number of arguments in the type
/// of a template specialisation.

static void KeepNParams(clang::QualType& normalizedType,
                        const clang::QualType& vanillaType,
                        const cling::Interpreter& interp,
                        const ROOT::TMetaUtils::TNormalizedCtxt& normCtxt)
{
   using namespace ROOT::TMetaUtils;
   using namespace clang;

   // If this type has no template specialisation behind, we don't need to do
   // anything
   ClassTemplateSpecializationDecl* ctsd;
   ClassTemplateDecl* ctd;
   if (! QualType2Template(vanillaType, ctd,  ctsd)) return ;

   // Even if this is a template, if we don't keep any argument, return
   const int nArgsToKeep = normCtxt.GetNargsToKeep(ctd);

   // Important in case of early return: we must restore the original qualtype
   QualType originalNormalizedType = normalizedType;

   const ASTContext& astCtxt = ctsd->getASTContext();


   // In case of name* we need to strip the pointer first, add the default and attach
   // the pointer once again.
   if (llvm::isa<clang::PointerType>(normalizedType.getTypePtr())) {
      // Get the qualifiers.
      clang::Qualifiers quals = normalizedType.getQualifiers();
      auto valNormalizedType = normalizedType->getPointeeType();
      KeepNParams(valNormalizedType,vanillaType, interp, normCtxt);
      normalizedType = astCtxt.getPointerType(valNormalizedType);
      // Add back the qualifiers.
      normalizedType = astCtxt.getQualifiedType(normalizedType, quals);
      return;
   }

   // In case of Int_t& we need to strip the pointer first, desugar and attach
   // the pointer once again.
   if (llvm::isa<clang::ReferenceType>(normalizedType.getTypePtr())) {
      // Get the qualifiers.
      bool isLValueRefTy = llvm::isa<clang::LValueReferenceType>(normalizedType.getTypePtr());
      clang::Qualifiers quals = normalizedType.getQualifiers();
      auto valNormType = normalizedType->getPointeeType();
      KeepNParams(valNormType, vanillaType, interp, normCtxt);

      // Add the r- or l- value reference type back to the desugared one
      if (isLValueRefTy)
        normalizedType = astCtxt.getLValueReferenceType(valNormType);
      else
        normalizedType = astCtxt.getRValueReferenceType(valNormType);
      // Add back the qualifiers.
      normalizedType = astCtxt.getQualifiedType(normalizedType, quals);
      return;
   }

   // Treat the Scope (factorise the code out to reuse it in AddDefaultParameters)
   bool prefix_changed = false;
   clang::NestedNameSpecifier* prefix = nullptr;
   clang::Qualifiers prefix_qualifiers = normalizedType.getLocalQualifiers();
   const clang::ElaboratedType* etype
      = llvm::dyn_cast<clang::ElaboratedType>(normalizedType.getTypePtr());
   if (etype) {
      // We have to also handle the prefix.
      // TODO: we ought to be running KeepNParams
      prefix = AddDefaultParametersNNS(astCtxt, etype->getQualifier(), interp, normCtxt);
      prefix_changed = prefix != etype->getQualifier();
      normalizedType = clang::QualType(etype->getNamedType().getTypePtr(),0);
   }

   // The canonical decl does not necessarily have the template default arguments.
   // Need to walk through the redecl chain to find it (we know there will be no
   // inconsistencies, at least)
   const clang::ClassTemplateDecl* ctdWithDefaultArgs = ctd;
   for (const RedeclarableTemplateDecl* rd: ctdWithDefaultArgs->redecls()) {
      clang::TemplateParameterList* tpl = rd->getTemplateParameters();
      if (tpl->getMinRequiredArguments () < tpl->size())
         ctdWithDefaultArgs = llvm::dyn_cast<clang::ClassTemplateDecl>(rd);
   }
   TemplateParameterList* tParsPtr = ctdWithDefaultArgs->getTemplateParameters();
   const TemplateParameterList& tPars = *tParsPtr;
   const TemplateArgumentList& tArgs = ctsd->getTemplateArgs();

   // We extract the template name from the type
   TemplateName theTemplateName = ExtractTemplateNameFromQualType(normalizedType);
   if (theTemplateName.isNull()) {
      normalizedType=originalNormalizedType;
      return;
   }

   const TemplateSpecializationType* normalizedTst =
            llvm::dyn_cast<TemplateSpecializationType>(normalizedType.getTypePtr());
   if (!normalizedTst) {
      normalizedType=originalNormalizedType;
      return;
   }

   const clang::ClassTemplateSpecializationDecl* TSTdecl
      = llvm::dyn_cast_or_null<const clang::ClassTemplateSpecializationDecl>(normalizedType.getTypePtr()->getAsCXXRecordDecl());
   bool isStdDropDefault = TSTdecl && IsStdDropDefaultClass(*TSTdecl);

   // Loop over the template parameters and arguments recursively.
   // We go down the two lanes: the one of template parameters (decls) and the
   // one of template arguments (QualTypes) in parallel. The former are a
   // property of the template, independent of its instantiations.
   // The latter are a property of the instance itself.
   llvm::SmallVector<TemplateArgument, 4> argsToKeep;

   const int nArgs = tArgs.size();
   const int nNormArgs = normalizedTst->getNumArgs();

   bool mightHaveChanged = false;

   // becomes true when a parameter has a value equal to its default
   for (int formal = 0, inst = 0; formal != nArgs; ++formal, ++inst) {
      const NamedDecl* tParPtr = tPars.getParam(formal);
      if (!tParPtr) {
         Error("KeepNParams", "The parameter number %s is null.\n", formal);
         continue;
      }

      // Stop if the normalized TemplateSpecializationType has less arguments than
      // the one index is pointing at.
      // We piggy back on the AddDefaultParameters routine basically.
      if (formal == nNormArgs || inst == nNormArgs) break;

      const TemplateArgument& tArg = tArgs.get(formal);
      TemplateArgument normTArg(normalizedTst->getArgs()[inst]);

      bool shouldKeepArg = nArgsToKeep < 0 || inst < nArgsToKeep;
      if (isStdDropDefault) shouldKeepArg = false;

      // Nothing to do here: either this parameter has no default, or we have to keep it.
      // FIXME: Temporary measure to get Atlas started with this.
      // We put a hard cut on the number of template arguments to keep, w/o checking if
      // they are non default. This makes this feature UNUSABLE for cases like std::vector,
      // where 2 different entities would have the same name if an allocator different from
      // the default one is by chance used.
      if (!isTypeWithDefault(tParPtr) || shouldKeepArg) {
         if ( tParPtr->isTemplateParameterPack() ) {
            // This is the last template parameter in the template declaration
            // but it is signaling that there can be an arbitrary number of arguments
            // in the template instance.  So to avoid inadvertenly dropping those
            // arguments we just process all remaining argument and exit the main loop.
            for( ; inst != nNormArgs; ++inst) {
               normTArg = normalizedTst->getArgs()[inst];
               mightHaveChanged |= RecurseKeepNParams(normTArg, tArg, interp, normCtxt, astCtxt);
               argsToKeep.push_back(normTArg);
            }
            // Done.
            break;
         }
         mightHaveChanged |= RecurseKeepNParams(normTArg, tArg, interp, normCtxt, astCtxt);
         argsToKeep.push_back(normTArg);
         continue;
      } else {
         if (!isStdDropDefault) {
            // Here we should not break but rather check if the value is the default one.
            mightHaveChanged = true;
            break;
         }
         // For std, we want to check the default args values.
      }

      // Now, we keep it only if it not is equal to its default, expressed in the arg
      // Some gymnastic is needed to decide how to check for equality according to the
      // flavour of Type: templateType or Integer
      bool equal=false;
      auto argKind = tArg.getKind();
      if (argKind == clang::TemplateArgument::Type){
         // we need all the info
         equal = areEqualTypes(tArg, argsToKeep, *tParPtr, interp, normCtxt);
      } else if (argKind == clang::TemplateArgument::Integral){
         equal = areEqualValues(tArg, *tParPtr);
      }
      if (!equal) {
         mightHaveChanged |= RecurseKeepNParams(normTArg, tArg, interp, normCtxt, astCtxt);
         argsToKeep.push_back(normTArg);
      } else {
         mightHaveChanged = true;
      }


   } // of loop over parameters and arguments

   if (!prefix_changed && !mightHaveChanged) {
      normalizedType = originalNormalizedType;
      return;
   }

   // now, let's remanipulate our Qualtype
   if (mightHaveChanged) {
      Qualifiers qualifiers = normalizedType.getLocalQualifiers();
      normalizedType = astCtxt.getTemplateSpecializationType(theTemplateName,
                                                             argsToKeep,
                                                             normalizedType.getTypePtr()->getCanonicalTypeInternal());
      normalizedType = astCtxt.getQualifiedType(normalizedType, qualifiers);
   }

   // Here we have (prefix_changed==true || mightHaveChanged), in both case
   // we need to reconstruct the type.
   if (prefix) {
      normalizedType = astCtxt.getElaboratedType(clang::ETK_None,prefix,normalizedType);
      normalizedType = astCtxt.getQualifiedType(normalizedType,prefix_qualifiers);
   }

}

////////////////////////////////////////////////////////////////////////////////
/// Return the type normalized for ROOT,
/// keeping only the ROOT opaque typedef (Double32_t, etc.) and
/// adding default template argument for all types except those explicitly
/// requested to be drop by the user.
/// Default template for STL collections are not yet removed by this routine.

clang::QualType ROOT::TMetaUtils::GetNormalizedType(const clang::QualType &type, const cling::Interpreter &interpreter, const TNormalizedCtxt &normCtxt)
{
   clang::ASTContext &ctxt = interpreter.getCI()->getASTContext();

   // Modules can trigger deserialization.
   cling::Interpreter::PushTransactionRAII RAII(const_cast<cling::Interpreter*>(&interpreter));
   clang::QualType normalizedType = cling::utils::Transform::GetPartiallyDesugaredType(ctxt, type, normCtxt.GetConfig(), true /* fully qualify */);

   // Readd missing default template parameters
   normalizedType = ROOT::TMetaUtils::AddDefaultParameters(normalizedType, interpreter, normCtxt);

   // Get the number of arguments to keep in case they are not default.
   KeepNParams(normalizedType,type,interpreter,normCtxt);

   return normalizedType;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the type name normalized for ROOT,
/// keeping only the ROOT opaque typedef (Double32_t, etc.) and
/// adding default template argument for all types except the STL collections
/// where we remove the default template argument if any.
///
/// This routine might actually belong in the interpreter because
/// cache the clang::Type might be intepreter specific.

void ROOT::TMetaUtils::GetNormalizedName(std::string &norm_name, const clang::QualType &type, const cling::Interpreter &interpreter, const TNormalizedCtxt &normCtxt)
{
   if (type.isNull()) {
      norm_name = "";
      return;
   }

   clang::QualType normalizedType = GetNormalizedType(type,interpreter,normCtxt);

   clang::ASTContext &ctxt = interpreter.getCI()->getASTContext();
   clang::PrintingPolicy policy(ctxt.getPrintingPolicy());
   policy.SuppressTagKeyword = true; // Never get the class or struct keyword
   policy.SuppressScope = true;      // Force the scope to be coming from a clang::ElaboratedType.
   policy.AnonymousTagLocations = false; // Do not extract file name + line number for anonymous types.
   // The scope suppression is required for getting rid of the anonymous part of the name of a class defined in an anonymous namespace.
   // This gives us more control vs not using the clang::ElaboratedType and relying on the Policy.SuppressUnwrittenScope which would
   // strip both the anonymous and the inline namespace names (and we probably do not want the later to be suppressed).

   std::string normalizedNameStep1;

   // getAsStringInternal can trigger deserialization
   cling::Interpreter::PushTransactionRAII clingRAII(const_cast<cling::Interpreter*>(&interpreter));
   normalizedType.getAsStringInternal(normalizedNameStep1,policy);

   // Still remove the std:: and default template argument for STL container and
   // normalize the location and amount of white spaces.
   TClassEdit::TSplitType splitname(normalizedNameStep1.c_str(),(TClassEdit::EModType)(TClassEdit::kLong64 | TClassEdit::kDropStd | TClassEdit::kDropStlDefault | TClassEdit::kKeepOuterConst));
   splitname.ShortType(norm_name,TClassEdit::kDropStd | TClassEdit::kDropStlDefault );

   // The result of this routine is by definition a fully qualified name.  There is an implicit starting '::' at the beginning of the name.
   // Depending on how the user typed their code, in particular typedef declarations, we may end up with an explicit '::' being
   // part of the result string.  For consistency, we must remove it.
   if (norm_name.length()>2 && norm_name[0]==':' && norm_name[1]==':') {
      norm_name.erase(0,2);
   }

}

////////////////////////////////////////////////////////////////////////////////

void ROOT::TMetaUtils::GetNormalizedName(std::string &norm_name,
                                         const clang::TypeDecl* typeDecl,
                                         const cling::Interpreter &interpreter)
{
   ROOT::TMetaUtils::TNormalizedCtxt tNormCtxt(interpreter.getLookupHelper());
   const clang::Sema &sema = interpreter.getSema();
   clang::ASTContext& astCtxt = sema.getASTContext();
   clang::QualType qualType = astCtxt.getTypeDeclType(typeDecl);

   ROOT::TMetaUtils::GetNormalizedName(norm_name,
                                       qualType,
                                       interpreter,
                                       tNormCtxt);
}

////////////////////////////////////////////////////////////////////////////////
std::pair<std::string,clang::QualType>
ROOT::TMetaUtils::GetNameTypeForIO(const clang::QualType& thisType,
                                   const cling::Interpreter &interpreter,
                                   const TNormalizedCtxt &normCtxt,
                                   TClassEdit::EModType mode)
{
   std::string thisTypeName;
   GetNormalizedName(thisTypeName, thisType, interpreter, normCtxt );
   bool hasChanged;
   auto thisTypeNameForIO = TClassEdit::GetNameForIO(thisTypeName, mode, &hasChanged);
   if (!hasChanged) return std::make_pair(thisTypeName,thisType);

   if (hasChanged && ROOT::TMetaUtils::GetErrorIgnoreLevel() <= ROOT::TMetaUtils::kNote) {
      ROOT::TMetaUtils::Info("ROOT::TMetaUtils::GetTypeForIO",
        "Name changed from %s to %s\n", thisTypeName.c_str(), thisTypeNameForIO.c_str());
   }

   auto& lookupHelper = interpreter.getLookupHelper();

   const clang::Type* typePtrForIO;
   lookupHelper.findScope(thisTypeNameForIO,
                          cling::LookupHelper::DiagSetting::NoDiagnostics,
                          &typePtrForIO);

   // This should never happen
   if (!typePtrForIO) {
      ROOT::TMetaUtils::Fatal("ROOT::TMetaUtils::GetTypeForIO",
                              "Type not found: %s.",thisTypeNameForIO.c_str());
   }

   clang::QualType typeForIO(typePtrForIO,0);

   // Check if this is a class. Indeed it could well be a POD
   if (!typeForIO->isRecordType()) {
      return std::make_pair(thisTypeNameForIO,typeForIO);
   }

   auto thisDeclForIO = typeForIO->getAsCXXRecordDecl();
   if (!thisDeclForIO) {
      ROOT::TMetaUtils::Error("ROOT::TMetaUtils::GetTypeForIO",
       "The type for IO corresponding to %s is %s and it could not be found in the AST as class.\n", thisTypeName.c_str(), thisTypeNameForIO.c_str());
      return std::make_pair(thisTypeName,thisType);
   }

   return std::make_pair(thisTypeNameForIO,typeForIO);
}

////////////////////////////////////////////////////////////////////////////////

clang::QualType ROOT::TMetaUtils::GetTypeForIO(const clang::QualType& thisType,
                                               const cling::Interpreter &interpreter,
                                               const TNormalizedCtxt &normCtxt,
                                               TClassEdit::EModType mode)
{
   return GetNameTypeForIO(thisType, interpreter, normCtxt, mode).second;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the dictionary file name for a module

std::string ROOT::TMetaUtils::GetModuleFileName(const char* moduleName)
{
   std::string dictFileName(moduleName);
   dictFileName += "_rdict.pcm";
   return dictFileName;
}

int dumpDeclForAssert(const clang::Decl& D, const char* commentStart) {
   llvm::errs() << llvm::StringRef(commentStart, 80) << '\n';
   D.dump();
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns the comment (// striped away), annotating declaration in a meaningful
/// for ROOT IO way.
/// Takes optional out parameter clang::SourceLocation returning the source
/// location of the comment.
///
/// CXXMethodDecls, FieldDecls and TagDecls are annotated.
/// CXXMethodDecls declarations and FieldDecls are annotated as follows:
/// Eg. void f(); // comment1
///     int member; // comment2
/// Inline definitions of CXXMethodDecls after the closing } \n. Eg:
/// void f()
/// {...}  // comment3
/// TagDecls are annotated in the end of the ClassDef macro. Eg.
/// class MyClass {
/// ...
/// ClassDef(MyClass, 1) // comment4
///

llvm::StringRef ROOT::TMetaUtils::GetComment(const clang::Decl &decl, clang::SourceLocation *loc)
{
   clang::SourceManager& sourceManager = decl.getASTContext().getSourceManager();
   clang::SourceLocation sourceLocation = decl.getLocEnd();

   // If the location is a macro get the expansion location.
   sourceLocation = sourceManager.getExpansionRange(sourceLocation).second;
   // FIXME: We should optimize this routine instead making it do the wrong thing
   // returning an empty comment if the decl came from the AST.
   // In order to do that we need to: check if the decl has an attribute and
   // return the attribute content (including walking the redecl chain) and if
   // this is not the case we should try finding it in the header file.
   // This will allow us to move the implementation of TCling*Info::Title() in
   // TClingDeclInfo.
   if (!decl.hasOwningModule() && sourceManager.isLoadedSourceLocation(sourceLocation)) {
      // Do not touch disk for nodes coming from the PCH.
      return "";
   }

   bool invalid;
   const char *commentStart = sourceManager.getCharacterData(sourceLocation, &invalid);
   if (invalid)
      return "";

   bool skipToSemi = true;
   if (const clang::FunctionDecl* FD = clang::dyn_cast<clang::FunctionDecl>(&decl)) {
      if (FD->isImplicit()) {
         // Compiler generated function.
         return "";
      }
      if (FD->isExplicitlyDefaulted() || FD->isDeletedAsWritten()) {
         // ctorOrFunc() = xyz; with commentStart pointing somewhere into
         // ctorOrFunc.
         // We have to skipToSemi
      } else if (FD->doesThisDeclarationHaveABody()) {
         // commentStart is at body's '}'
         // But we might end up e.g. at the ')' of a CPP macro
         assert((decl.getLocEnd() != sourceLocation || *commentStart == '}'
                 || dumpDeclForAssert(*FD, commentStart))
                && "Expected macro or end of body at '}'");
         if (*commentStart) ++commentStart;

         // We might still have a ';'; skip the spaces and check.
         while (*commentStart && isspace(*commentStart)
                && *commentStart != '\n' && *commentStart != '\r') {
            ++commentStart;
         }
         if (*commentStart == ';') ++commentStart;

         skipToSemi = false;
      }
   } else if (const clang::EnumConstantDecl* ECD
              = clang::dyn_cast<clang::EnumConstantDecl>(&decl)) {
      // either "konstant = 12, //COMMENT" or "lastkonstant // COMMENT"
      if (ECD->getNextDeclInContext())
         while (*commentStart && *commentStart != ',' && *commentStart != '\r' && *commentStart != '\n')
            ++commentStart;
      // else commentStart already points to the end.

      skipToSemi = false;
   }

   if (skipToSemi) {
      while (*commentStart && *commentStart != ';' && *commentStart != '\r' && *commentStart != '\n')
         ++commentStart;
      if (*commentStart == ';') ++commentStart;
   }

   // Now skip the spaces until beginning of comments or EOL.
   while ( *commentStart && isspace(*commentStart)
           && *commentStart != '\n' && *commentStart != '\r') {
      ++commentStart;
   }

   if (commentStart[0] != '/' ||
       (commentStart[1] != '/' && commentStart[1] != '*')) {
      // not a comment
      return "";
   }

   // Treat by default c++ comments (+2) but also Doxygen comments (+4)
   //   Int_t fPx; ///< Some doxygen comment for persistent data.
   //   Int_t fPy; //!< Some doxygen comment for persistent data.
   //   Int_t fPz; /*!< Some doxygen comment for persistent data. */
   //   Int_t fPa; /**< Some doxygen comment for persistent data. */
   unsigned int skipChars = 2;
   if (commentStart[0] == '/' &&
       commentStart[1] == '/' &&
       (commentStart[2] == '/' || commentStart[2] == '!') &&
       commentStart[3] == '<') {
      skipChars = 4;
   } else if (commentStart[0] == '/' &&
              commentStart[1] == '*' &&
              (commentStart[2] == '*' || commentStart[2] == '!') &&
              commentStart[3] == '<') {
      skipChars = 4;
   }

   commentStart += skipChars;

   // Now skip the spaces after comment start until EOL.
   while ( *commentStart && isspace(*commentStart)
           && *commentStart != '\n' && *commentStart != '\r') {
      ++commentStart;
   }
   const char* commentEnd = commentStart;
   // Even for /* comments we only take the first line into account.
   while (*commentEnd && *commentEnd != '\n' && *commentEnd != '\r') {
      ++commentEnd;
   }

   // "Skip" (don't include) trailing space.
   // *commentEnd points behind comment end thus check commentEnd[-1]
   while (commentEnd > commentStart && isspace(commentEnd[-1])) {
      --commentEnd;
   }

   if (loc) {
      // Find the true beginning of a comment.
      unsigned offset = commentStart - sourceManager.getCharacterData(sourceLocation);
      *loc = sourceLocation.getLocWithOffset(offset - 1);
   }

   return llvm::StringRef(commentStart, commentEnd - commentStart);
}

////////////////////////////////////////////////////////////////////////////////
/// Return the class comment after the ClassDef:
/// class MyClass {
/// ...
/// ClassDef(MyClass, 1) // class comment
///

llvm::StringRef ROOT::TMetaUtils::GetClassComment(const clang::CXXRecordDecl &decl,
                                                  clang::SourceLocation *loc,
                                                  const cling::Interpreter &interpreter)
{
   using namespace clang;
   SourceLocation commentSLoc;
   llvm::StringRef comment;

   Sema& sema = interpreter.getCI()->getSema();

   const Decl* DeclFileLineDecl
      = interpreter.getLookupHelper().findFunctionProto(&decl, "DeclFileLine", "",
                                                        cling::LookupHelper::NoDiagnostics);
   if (!DeclFileLineDecl) return llvm::StringRef();

   // For now we allow only a special macro (ClassDef) to have meaningful comments
   SourceLocation maybeMacroLoc = DeclFileLineDecl->getLocation();
   bool isClassDefMacro = maybeMacroLoc.isMacroID() && sema.findMacroSpelling(maybeMacroLoc, "ClassDef");
   if (isClassDefMacro) {
      comment = ROOT::TMetaUtils::GetComment(*DeclFileLineDecl, &commentSLoc);
      if (comment.size()) {
         if (loc){
            *loc = commentSLoc;
         }
         return comment;
      }
   }
   return llvm::StringRef();
}

////////////////////////////////////////////////////////////////////////////////
/// Return the base/underlying type of a chain of array or pointers type.
/// Does not yet support the array and pointer part being intermixed.

const clang::Type *ROOT::TMetaUtils::GetUnderlyingType(clang::QualType type)
{
   const clang::Type *rawtype = type.getTypePtr();

   // NOTE: We probably meant isa<clang::ElaboratedType>
   if (rawtype->isElaboratedTypeSpecifier() ) {
      rawtype = rawtype->getCanonicalTypeInternal().getTypePtr();
   }
   if (rawtype->isArrayType()) {
      rawtype = type.getTypePtr()->getBaseElementTypeUnsafe ();
   }
   if (rawtype->isPointerType() || rawtype->isReferenceType() ) {
      //Get to the 'raw' type.
      clang::QualType pointee;
      while ( (pointee = rawtype->getPointeeType()) , pointee.getTypePtrOrNull() && pointee.getTypePtr() != rawtype)
      {
         rawtype = pointee.getTypePtr();

         if (rawtype->isElaboratedTypeSpecifier() ) {
            rawtype = rawtype->getCanonicalTypeInternal().getTypePtr();
         }
         if (rawtype->isArrayType()) {
            rawtype = rawtype->getBaseElementTypeUnsafe ();
         }
      }
   }
   if (rawtype->isArrayType()) {
      rawtype = rawtype->getBaseElementTypeUnsafe ();
   }
   return rawtype;
}

////////////////////////////////////////////////////////////////////////////////
/// Return true, if the decl is part of the std namespace.

bool ROOT::TMetaUtils::IsStdClass(const clang::RecordDecl &cl)
{
  return cling::utils::Analyze::IsStdClass(cl);
}

////////////////////////////////////////////////////////////////////////////////
/// Return true, if the decl is part of the std namespace and we want
/// its default parameter dropped.

bool ROOT::TMetaUtils::IsStdDropDefaultClass(const clang::RecordDecl &cl)
{
   // Might need to reduce it to shared_ptr and STL collection.s
   if (cling::utils::Analyze::IsStdClass(cl)) {
      static const char *names[] =
      {  "shared_ptr", "__shared_ptr",
         "vector", "list", "deque", "map", "multimap", "set", "multiset", "bitset"};
      llvm::StringRef clname(cl.getName());
      for(auto &&name : names) {
         if (clname == name) return true;
      }
   }
   return false;
}

////////////////////////////////////////////////////////////////////////////////
/// This is a recursive function

bool ROOT::TMetaUtils::MatchWithDeclOrAnyOfPrevious(const clang::CXXRecordDecl &cl,
                                                    const clang::CXXRecordDecl &currentCl)
{
   // We found it: let's return true
   if (&cl == &currentCl) return true;

   const  clang::CXXRecordDecl* previous = currentCl.getPreviousDecl();

   // There is no previous decl, so we cannot possibly find it
   if (NULL == previous){
      return false;
   }

   // We try to find it in the previous
   return ROOT::TMetaUtils::MatchWithDeclOrAnyOfPrevious(cl, *previous);

}

//______________________________________________________________________________

bool ROOT::TMetaUtils::IsOfType(const clang::CXXRecordDecl &cl, const std::string& typ, const cling::LookupHelper& lh)
{
   // Return true if the decl is of type.
   // A proper hashtable for caching results would be the ideal solution
   // 1) Only one lookup per type
   // 2) No string comparison
   // We may use a map which becomes an unordered map if c++11 is enabled?

   const clang::CXXRecordDecl *thisDecl =
      llvm::dyn_cast_or_null<clang::CXXRecordDecl>(lh.findScope(typ, cling::LookupHelper::WithDiagnostics));

   // this would be probably an assert given that this state is not reachable unless a mistake is somewhere
   if (! thisDecl){
      Error("IsOfType","Record decl of type %s not found in the AST.", typ.c_str());
      return false;
   }

   // Now loop on all previous decls to seek a match
   const clang::CXXRecordDecl *mostRecentDecl = thisDecl->getMostRecentDecl();
   bool matchFound = MatchWithDeclOrAnyOfPrevious (cl,*mostRecentDecl);

   return matchFound;
}

////////////////////////////////////////////////////////////////////////////////
///  type     : type name: vector<list<classA,allocator>,allocator>
///  result:    0          : not stl container
///             abs(result): code of container 1=vector,2=list,3=deque,4=map
///                           5=multimap,6=set,7=multiset

ROOT::ESTLType ROOT::TMetaUtils::IsSTLCont(const clang::RecordDecl &cl)
{
   // This routine could be enhanced to also support:
   //
   //  testAlloc: if true, we test allocator, if it is not default result is negative
   //  result:    0          : not stl container
   //             abs(result): code of container 1=vector,2=list,3=deque,4=map
   //                           5=multimap,6=set,7=multiset
   //             positive val: we have a vector or list with default allocator to any depth
   //                   like vector<list<vector<int>>>
   //             negative val: STL container other than vector or list, or non default allocator
   //                           For example: vector<deque<int>> has answer -1

   if (!IsStdClass(cl)) {
      return ROOT::kNotSTL;
   }

   return STLKind(cl.getName());
}

static bool hasSomeTypedefSomewhere(const clang::Type* T) {
  using namespace clang;
  struct SearchTypedef: public TypeVisitor<SearchTypedef, bool> {
    bool VisitTypedefType(const TypedefType* TD) {
      return true;
    }
    bool VisitArrayType(const ArrayType* AT) {
      return Visit(AT->getElementType().getTypePtr());
    }
    bool VisitDecltypeType(const DecltypeType* DT) {
      return Visit(DT->getUnderlyingType().getTypePtr());
    }
    bool VisitPointerType(const PointerType* PT) {
      return Visit(PT->getPointeeType().getTypePtr());
    }
    bool VisitReferenceType(const ReferenceType* RT) {
      return Visit(RT->getPointeeType().getTypePtr());
    }
    bool VisitSubstTemplateTypeParmType(const SubstTemplateTypeParmType* STST) {
      return Visit(STST->getReplacementType().getTypePtr());
    }
    bool VisitTemplateSpecializationType(const TemplateSpecializationType* TST) {
      for (int I = 0, N = TST->getNumArgs(); I < N; ++I) {
        const TemplateArgument& TA = TST->getArg(I);
        if (TA.getKind() == TemplateArgument::Type
            && Visit(TA.getAsType().getTypePtr()))
          return true;
      }
      return false;
    }
    bool VisitTemplateTypeParmType(const TemplateTypeParmType* TTPT) {
      return false; // shrug...
    }
    bool VisitTypeOfType(const TypeOfType* TOT) {
      return TOT->getUnderlyingType().getTypePtr();
    }
    bool VisitElaboratedType(const ElaboratedType* ET) {
      NestedNameSpecifier* NNS = ET->getQualifier();
      while (NNS) {
        if (NNS->getKind() == NestedNameSpecifier::TypeSpec) {
          if (Visit(NNS->getAsType()))
            return true;
        }
        NNS = NNS->getPrefix();
      }
      return Visit(ET->getNamedType().getTypePtr());
    }
  };

  SearchTypedef ST;
  return ST.Visit(T);
}

////////////////////////////////////////////////////////////////////////////////
/// Check if 'input' or any of its template parameter was substituted when
/// instantiating the class template instance and replace it with the
/// partially sugared types we have from 'instance'.

clang::QualType ROOT::TMetaUtils::ReSubstTemplateArg(clang::QualType input, const clang::Type *instance)
{
   if (!instance) return input;
   // if there is no typedef in instance then there is nothing guiding any
   // template parameter typedef replacement.
   if (!hasSomeTypedefSomewhere(instance))
     return input;

   using namespace llvm;
   using namespace clang;
   const clang::ASTContext &Ctxt = instance->getAsCXXRecordDecl()->getASTContext();

   // Treat scope (clang::ElaboratedType) if any.
   const clang::ElaboratedType* etype
      = llvm::dyn_cast<clang::ElaboratedType>(input.getTypePtr());
   if (etype) {
      // We have to also handle the prefix.

      clang::Qualifiers scope_qualifiers = input.getLocalQualifiers();
      assert(instance->getAsCXXRecordDecl()!=0 && "ReSubstTemplateArg only makes sense with a type representing a class.");

      clang::NestedNameSpecifier *scope = ReSubstTemplateArgNNS(Ctxt,etype->getQualifier(),instance);
      clang::QualType subTy = ReSubstTemplateArg(clang::QualType(etype->getNamedType().getTypePtr(),0),instance);

      if (scope) subTy = Ctxt.getElaboratedType(clang::ETK_None,scope,subTy);
      subTy = Ctxt.getQualifiedType(subTy,scope_qualifiers);
      return subTy;
   }

   QualType QT = input;

   // In case of Int_t* we need to strip the pointer first, ReSubst and attach
   // the pointer once again.
   if (isa<clang::PointerType>(QT.getTypePtr())) {
      // Get the qualifiers.
      Qualifiers quals = QT.getQualifiers();
      QualType nQT;
      nQT = ReSubstTemplateArg(QT->getPointeeType(),instance);
      if (nQT == QT->getPointeeType()) return QT;

      QT = Ctxt.getPointerType(nQT);
      // Add back the qualifiers.
      QT = Ctxt.getQualifiedType(QT, quals);
      return QT;
   }

   // In case of Int_t& we need to strip the pointer first, ReSubst and attach
   // the reference once again.
   if (isa<ReferenceType>(QT.getTypePtr())) {
      // Get the qualifiers.
      bool isLValueRefTy = isa<LValueReferenceType>(QT.getTypePtr());
      Qualifiers quals = QT.getQualifiers();
      QualType nQT;
      nQT = ReSubstTemplateArg(QT->getPointeeType(),instance);
      if (nQT == QT->getPointeeType()) return QT;

      // Add the r- or l-value reference type back to the desugared one.
      if (isLValueRefTy)
         QT = Ctxt.getLValueReferenceType(nQT);
      else
         QT = Ctxt.getRValueReferenceType(nQT);
      // Add back the qualifiers.
      QT = Ctxt.getQualifiedType(QT, quals);
      return QT;
   }

   // In case of Int_t[2] we need to strip the array first, ReSubst and attach
   // the array once again.
   if (isa<clang::ArrayType>(QT.getTypePtr())) {
      // Get the qualifiers.
      Qualifiers quals = QT.getQualifiers();

      if (const auto arr = dyn_cast<ConstantArrayType>(QT.getTypePtr())) {
         QualType newQT= ReSubstTemplateArg(arr->getElementType(),instance);

         if (newQT == arr->getElementType()) return QT;
         QT = Ctxt.getConstantArrayType (newQT,
                                        arr->getSize(),
                                        arr->getSizeModifier(),
                                        arr->getIndexTypeCVRQualifiers());

      } else if (const auto arr = dyn_cast<DependentSizedArrayType>(QT.getTypePtr())) {
         QualType newQT = ReSubstTemplateArg(arr->getElementType(),instance);

         if (newQT == QT) return QT;
         QT = Ctxt.getDependentSizedArrayType (newQT,
                                              arr->getSizeExpr(),
                                              arr->getSizeModifier(),
                                              arr->getIndexTypeCVRQualifiers(),
                                              arr->getBracketsRange());

      } else if (const auto arr = dyn_cast<IncompleteArrayType>(QT.getTypePtr())) {
         QualType newQT = ReSubstTemplateArg(arr->getElementType(),instance);

         if (newQT == arr->getElementType()) return QT;
         QT = Ctxt.getIncompleteArrayType (newQT,
                                          arr->getSizeModifier(),
                                          arr->getIndexTypeCVRQualifiers());

      } else if (const auto arr = dyn_cast<VariableArrayType>(QT.getTypePtr())) {
         QualType newQT = ReSubstTemplateArg(arr->getElementType(),instance);

         if (newQT == arr->getElementType()) return QT;
         QT = Ctxt.getVariableArrayType (newQT,
                                        arr->getSizeExpr(),
                                        arr->getSizeModifier(),
                                        arr->getIndexTypeCVRQualifiers(),
                                        arr->getBracketsRange());
      }

      // Add back the qualifiers.
      QT = Ctxt.getQualifiedType(QT, quals);
      return QT;
   }

   // If the instance is also an elaborated type, we need to skip
   etype = llvm::dyn_cast<clang::ElaboratedType>(instance);
   if (etype) {
      instance = etype->getNamedType().getTypePtr();
      if (!instance) return input;
   }

   const clang::TemplateSpecializationType* TST
      = llvm::dyn_cast<const clang::TemplateSpecializationType>(instance);

   if (!TST) return input;

   const clang::ClassTemplateSpecializationDecl* TSTdecl
      = llvm::dyn_cast_or_null<const clang::ClassTemplateSpecializationDecl>(instance->getAsCXXRecordDecl());

   const clang::SubstTemplateTypeParmType *substType
      = llvm::dyn_cast<clang::SubstTemplateTypeParmType>(input.getTypePtr());

   if (substType) {
      // Make sure it got replaced from this template
      const clang::ClassTemplateDecl *replacedCtxt = 0;

      const clang::DeclContext *replacedDeclCtxt = substType->getReplacedParameter()->getDecl()->getDeclContext();
      const clang::CXXRecordDecl *decl = llvm::dyn_cast<clang::CXXRecordDecl>(replacedDeclCtxt);
      unsigned int index = substType->getReplacedParameter()->getIndex();
      if (decl) {

         if (decl->getKind() == clang::Decl::ClassTemplatePartialSpecialization) {
            const clang::ClassTemplatePartialSpecializationDecl *spec = llvm::dyn_cast<clang::ClassTemplatePartialSpecializationDecl>(decl);

            unsigned int depth = substType->getReplacedParameter()->getDepth();

            const TemplateArgument *instanceArgs = spec->getTemplateArgs().data();
            unsigned int instanceNArgs = spec->getTemplateArgs().size();

            // Search for the 'right' replacement.

            for(unsigned int A = 0; A < instanceNArgs; ++A) {
               if (instanceArgs[A].getKind() == clang::TemplateArgument::Type) {
                  clang::QualType argQualType = instanceArgs[A].getAsType();

                  const clang::TemplateTypeParmType *replacementType;

                  replacementType = llvm::dyn_cast<clang::TemplateTypeParmType>(argQualType);

                  if (!replacementType) {
                     const clang::SubstTemplateTypeParmType *argType
                        = llvm::dyn_cast<clang::SubstTemplateTypeParmType>(argQualType);
                     if (argType) {
                        clang::QualType replacementQT = argType->getReplacementType();
                        replacementType = llvm::dyn_cast<clang::TemplateTypeParmType>(replacementQT);
                     }
                  }
                  if (replacementType &&
                      depth == replacementType->getDepth() &&
                      index == replacementType->getIndex() )
                  {
                     index = A;
                     break;
                  }
               }
            }
            replacedCtxt = spec->getSpecializedTemplate();
         } else {
            replacedCtxt = decl->getDescribedClassTemplate();
         }
      } else if (auto const declguide = llvm::dyn_cast<clang::CXXDeductionGuideDecl>(replacedDeclCtxt)) {
         replacedCtxt = llvm::dyn_cast<clang::ClassTemplateDecl>(declguide->getDeducedTemplate());
      } else if (auto const ctdecl = llvm::dyn_cast<clang::ClassTemplateDecl>(replacedDeclCtxt)) {
         replacedCtxt = ctdecl;
      } else {
         std::string astDump;
         llvm::raw_string_ostream ostream(astDump);
         instance->dump(ostream);
         ostream.flush();
         ROOT::TMetaUtils::Warning("ReSubstTemplateArg","Unexpected type of declaration context for template parameter: %s.\n\tThe responsible class is:\n\t%s\n",
                                   replacedDeclCtxt->getDeclKindName(), astDump.c_str());
         replacedCtxt = nullptr;
      }

      if ((replacedCtxt && replacedCtxt->getCanonicalDecl() == TSTdecl->getSpecializedTemplate()->getCanonicalDecl())
          || /* the following is likely just redundant */
          substType->getReplacedParameter()->getDecl()
          == TSTdecl->getSpecializedTemplate ()->getTemplateParameters()->getParam(index))
      {
         if ( index >= TST->getNumArgs() ) {
            // The argument replaced was a default template argument that is
            // being listed as part of the instance ...
            // so we probably don't really know how to spell it ... we would need to recreate it
            // (See AddDefaultParameters).
            return input;
         } else {
            return TST->getArg(index).getAsType();
         }
      }
   }
   // Maybe a class template instance, recurse and rebuild
   const clang::TemplateSpecializationType* inputTST
      = llvm::dyn_cast<const clang::TemplateSpecializationType>(input.getTypePtr());
   const clang::ASTContext& astCtxt = TSTdecl->getASTContext();

   if (inputTST) {
      bool mightHaveChanged = false;
      llvm::SmallVector<clang::TemplateArgument, 4> desArgs;
      for(clang::TemplateSpecializationType::iterator I = inputTST->begin(), E = inputTST->end();
          I != E; ++I) {
         if (I->getKind() != clang::TemplateArgument::Type) {
            desArgs.push_back(*I);
            continue;
         }

         clang::QualType SubTy = I->getAsType();
         // Check if the type needs more desugaring and recurse.
         if (llvm::isa<clang::SubstTemplateTypeParmType>(SubTy)
             || llvm::isa<clang::TemplateSpecializationType>(SubTy)) {
            clang::QualType newSubTy = ReSubstTemplateArg(SubTy,instance);
            mightHaveChanged = SubTy != newSubTy;
            if (!newSubTy.isNull()) {
               desArgs.push_back(clang::TemplateArgument(newSubTy));
            }
         } else
            desArgs.push_back(*I);
      }

      // If desugaring happened allocate new type in the AST.
      if (mightHaveChanged) {
         clang::Qualifiers qualifiers = input.getLocalQualifiers();
         input = astCtxt.getTemplateSpecializationType(inputTST->getTemplateName(),
                                                       desArgs,
                                                       inputTST->getCanonicalTypeInternal());
         input = astCtxt.getQualifiedType(input, qualifiers);
      }
   }

   return input;
}

////////////////////////////////////////////////////////////////////////////////
/// Remove the last n template arguments from the name

int ROOT::TMetaUtils::RemoveTemplateArgsFromName(std::string& name, unsigned int nArgsToRemove)
{
   if ( nArgsToRemove == 0 || name == "")
      return 0;

   // We proceed from the right to the left, counting commas which are not
   // enclosed by < >.
   const unsigned int length = name.length();
   unsigned int cur=0; // let's start beyond the first > from the right
   unsigned int nArgsRemoved=0;
   unsigned int nBraces=0;
   char c='@';
   while (nArgsRemoved!=nArgsToRemove && cur<length){
      c = name[cur];
      if (c == '<') nBraces++;
      if (c == '>') nBraces--;
      if (c == ',' && nBraces==1 /*So we are not in a sub-template*/) nArgsRemoved++;
      cur++;
   }
   cur--;
   name = name.substr(0,cur)+">";
   return 0;

}

////////////////////////////////////////////////////////////////////////////////
/// Converts STL container name to number. vector -> 1, etc..

ROOT::ESTLType ROOT::TMetaUtils::STLKind(const llvm::StringRef type)
{
   static const char *stls[] =                  //container names
      {"any","vector","list", "deque","map","multimap","set","multiset","bitset",
         "forward_list","unordered_set","unordered_multiset","unordered_map","unordered_multimap",0};
   static const ROOT::ESTLType values[] =
      {ROOT::kNotSTL, ROOT::kSTLvector,
       ROOT::kSTLlist, ROOT::kSTLdeque,
       ROOT::kSTLmap, ROOT::kSTLmultimap,
       ROOT::kSTLset, ROOT::kSTLmultiset,
       ROOT::kSTLbitset,
       ROOT::kSTLforwardlist,
       ROOT::kSTLunorderedset, ROOT::kSTLunorderedmultiset,
       ROOT::kSTLunorderedmap, ROOT::kSTLunorderedmultimap,
       ROOT::kNotSTL
      };
   //              kind of stl container
   for(int k=1;stls[k];k++) {if (type.equals(stls[k])) return values[k];}
   return ROOT::kNotSTL;
}

////////////////////////////////////////////////////////////////////////////////

const clang::TypedefNameDecl* ROOT::TMetaUtils::GetAnnotatedRedeclarable(const clang::TypedefNameDecl* TND)
{
   if (!TND)
      return 0;

   TND = TND->getMostRecentDecl();
   while (TND && !(TND->hasAttrs()))
      TND = TND->getPreviousDecl();

   return TND;
}

////////////////////////////////////////////////////////////////////////////////

const clang::TagDecl* ROOT::TMetaUtils::GetAnnotatedRedeclarable(const clang::TagDecl* TD)
{
   if (!TD)
      return 0;

   TD = TD->getMostRecentDecl();
   while (TD && !(TD->hasAttrs() && TD->isThisDeclarationADefinition()))
      TD = TD->getPreviousDecl();

   return TD;
}

////////////////////////////////////////////////////////////////////////////////
/// Extract the immediately outer namespace and then launch the recursion

void ROOT::TMetaUtils::ExtractEnclosingNameSpaces(const clang::Decl& decl,
                                                  std::list<std::pair<std::string,bool> >& enclosingNamespaces)
{
   const clang::DeclContext* enclosingNamespaceDeclCtxt = decl.getDeclContext();
   if (!enclosingNamespaceDeclCtxt) return;

   const clang::NamespaceDecl* enclosingNamespace =
             clang::dyn_cast<clang::NamespaceDecl>(enclosingNamespaceDeclCtxt);
   if (!enclosingNamespace) return;

   enclosingNamespaces.push_back(std::make_pair(enclosingNamespace->getNameAsString(),
                                                enclosingNamespace->isInline()));

   ExtractCtxtEnclosingNameSpaces(*enclosingNamespace, enclosingNamespaces);

}

////////////////////////////////////////////////////////////////////////////////
/// Extract enclosing namespaces recursively

void ROOT::TMetaUtils::ExtractCtxtEnclosingNameSpaces(const clang::DeclContext& ctxt,
                                                      std::list<std::pair<std::string,bool> >& enclosingNamespaces)
{
   const clang::DeclContext* enclosingNamespaceDeclCtxt = ctxt.getParent ();

   // If no parent is found, nothing more to be done
   if (!enclosingNamespaceDeclCtxt) {
      return;
   }

   // Check if the parent is a namespace (it could be a class for example)
   // if not, nothing to be done here
   const clang::NamespaceDecl* enclosingNamespace = clang::dyn_cast<clang::NamespaceDecl>(enclosingNamespaceDeclCtxt);
   if (!enclosingNamespace) return;

   // Add to the list of parent namespaces
   enclosingNamespaces.push_back(std::make_pair(enclosingNamespace->getNameAsString(),
                                                enclosingNamespace->isInline()));

   // here the recursion
   ExtractEnclosingNameSpaces(*enclosingNamespace, enclosingNamespaces);
}

////////////////////////////////////////////////////////////////////////////////
/// Extract the names and types of containing scopes.
/// Stop if a class is met and return its pointer.

const clang::RecordDecl* ROOT::TMetaUtils::ExtractEnclosingScopes(const clang::Decl& decl,
                                                                  std::list<std::pair<std::string,unsigned int> >& enclosingSc)
{
   const clang::DeclContext* enclosingDeclCtxt = decl.getDeclContext();
   if (!enclosingDeclCtxt) return 0;

   unsigned int scopeType;

   if (auto enclosingNamespacePtr =
             clang::dyn_cast<clang::NamespaceDecl>(enclosingDeclCtxt)){
      scopeType= enclosingNamespacePtr->isInline() ? 1 : 0; // inline or simple namespace
      enclosingSc.push_back(std::make_pair(enclosingNamespacePtr->getNameAsString(),scopeType));
      return ExtractEnclosingScopes(*enclosingNamespacePtr, enclosingSc);
   }

   if (auto enclosingClassPtr =
             clang::dyn_cast<clang::RecordDecl>(enclosingDeclCtxt)){
      return enclosingClassPtr;
   }

   return nullptr;

}

////////////////////////////////////////////////////////////////////////////////
/// Reimplementation of TSystem::ExpandPathName() that cannot be
/// used from TMetaUtils.

static void replaceEnvVars(const char* varname, std::string& txt)
{
   std::string::size_type beginVar = 0;
   std::string::size_type endVar = 0;
   while ((beginVar = txt.find('$', beginVar)) != std::string::npos
          && beginVar + 1 < txt.length()) {
      std::string::size_type beginVarName = beginVar + 1;
      std::string::size_type endVarName = std::string::npos;
      if (txt[beginVarName] == '(') {
         // "$(VARNAME)" style.
         endVarName = txt.find(')', beginVarName);
         ++beginVarName;
         if (endVarName == std::string::npos) {
            ROOT::TMetaUtils::Error(0, "Missing ')' for '$(' in $%s at %s\n",
                                    varname, txt.c_str() + beginVar);
            return;
         }
         endVar = endVarName + 1;
      } else {
         // "$VARNAME/..." style.
         beginVarName = beginVar + 1;
         endVarName = beginVarName;
         while (isalnum(txt[endVarName]) || txt[endVarName] == '_')
            ++endVarName;
         endVar = endVarName;
      }

      const char* val = getenv(txt.substr(beginVarName,
                                          endVarName - beginVarName).c_str());
      if (!val) val = "";

      txt.replace(beginVar, endVar - beginVar, val);
      int lenval = strlen(val);
      int delta = lenval - (endVar - beginVar); // these many extra chars,
      endVar += delta; // advance the end marker accordingly.

      // Look for the next one
      beginVar = endVar + 1;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Organise the parameters for cling in order to guarantee relocatability
/// It treats the gcc toolchain and the root include path
/// FIXME: enables relocatability for experiments' framework headers until PCMs
/// are available.

void ROOT::TMetaUtils::SetPathsForRelocatability(std::vector<std::string>& clingArgs )
{
   const char* envInclPath = getenv("ROOT_INCLUDE_PATH");

   if (!envInclPath)
      return;
   std::istringstream envInclPathsStream(envInclPath);
   std::string inclPath;
   while (std::getline(envInclPathsStream, inclPath, ':')) {
      // Can't use TSystem in here; re-implement TSystem::ExpandPathName().
      replaceEnvVars("ROOT_INCLUDE_PATH", inclPath);
      if (!inclPath.empty()) {
         clingArgs.push_back("-I");
         clingArgs.push_back(inclPath);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////

void ROOT::TMetaUtils::ReplaceAll(std::string& str, const std::string& from, const std::string& to,bool recurse)
{
   if(from.empty())
      return;
   size_t start_pos = 0;
   bool changed=true;
   while (changed){
      changed=false;
      start_pos = 0;
      while((start_pos = str.find(from, start_pos)) != std::string::npos) {
         str.replace(start_pos, from.length(), to);
         start_pos += to.length();
         if (recurse) changed = true;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return the separator suitable for this platform.
const std::string& ROOT::TMetaUtils::GetPathSeparator()
{
   return ROOT::FoundationUtils::GetPathSeparator();
}

////////////////////////////////////////////////////////////////////////////////

bool ROOT::TMetaUtils::EndsWith(const std::string &theString, const std::string &theSubstring)
{
   if (theString.size() < theSubstring.size()) return false;
   const unsigned int theSubstringSize = theSubstring.size();
   return 0 == theString.compare(theString.size() - theSubstringSize,
                                 theSubstringSize,
                                 theSubstring);
}

////////////////////////////////////////////////////////////////////////////////

bool ROOT::TMetaUtils::BeginsWith(const std::string &theString, const std::string &theSubstring)
{
   if (theString.size() < theSubstring.size()) return false;
   const unsigned int theSubstringSize = theSubstring.size();
   return 0 == theString.compare(0,
                                 theSubstringSize,
                                 theSubstring);
}



////////////////////////////////////////////////////////////////////////////////

bool ROOT::TMetaUtils::IsLinkdefFile(const char *filename)
{
   // Note, should change this into take llvm::StringRef.

   if ((strstr(filename, "LinkDef") || strstr(filename, "Linkdef") ||
         strstr(filename, "linkdef")) && strstr(filename, ".h")) {
      return true;
   }
   size_t len = strlen(filename);
   size_t linkdeflen = 9; /* strlen("linkdef.h") */
   if (len >= 9) {
      if (0 == strncasecmp(filename + (len - linkdeflen), "linkdef", linkdeflen - 2)
            && 0 == strcmp(filename + (len - 2), ".h")
         ) {
         return true;
      } else {
         return false;
      }
   } else {
      return false;
   }
}

////////////////////////////////////////////////////////////////////////////////

bool ROOT::TMetaUtils::IsHeaderName(const std::string &filename)
{
   return llvm::sys::path::extension(filename) == ".h" ||
          llvm::sys::path::extension(filename) == ".hh" ||
          llvm::sys::path::extension(filename) == ".hpp" ||
          llvm::sys::path::extension(filename) == ".H" ||
          llvm::sys::path::extension(filename) == ".h++" ||
          llvm::sys::path::extension(filename) == "hxx" ||
          llvm::sys::path::extension(filename) == "Hxx" ||
          llvm::sys::path::extension(filename) == "HXX";
}

////////////////////////////////////////////////////////////////////////////////

const std::string ROOT::TMetaUtils::AST2SourceTools::Decls2FwdDecls(const std::vector<const clang::Decl *> &decls, cling::Interpreter::IgnoreFilesFunc_t ignoreFiles, const cling::Interpreter &interp)
{
   clang::Sema &sema = interp.getSema();
   cling::Transaction theTransaction(sema);
   std::set<clang::Decl *> addedDecls;
   for (auto decl : decls) {
      // again waiting for cling
      clang::Decl *ncDecl = const_cast<clang::Decl *>(decl);
      theTransaction.append(ncDecl);
   }
   std::string newFwdDecl;
   llvm::raw_string_ostream llvmOstr(newFwdDecl);
   interp.forwardDeclare(theTransaction, sema.getPreprocessor(), sema.getASTContext(), llvmOstr, true, nullptr, ignoreFiles);
   llvmOstr.flush();
   return newFwdDecl;
}

////////////////////////////////////////////////////////////////////////////////
/// Take the namespaces which enclose the decl and put them around the
/// definition string.
/// For example, if the definition string is "myClass" which is enclosed by
/// the namespaces ns1 and ns2, one would get:
/// namespace ns2{ namespace ns1 { class myClass; } }

int ROOT::TMetaUtils::AST2SourceTools::EncloseInNamespaces(const clang::Decl& decl,
                                                           std::string& defString)
{
   auto rcd = EncloseInScopes(decl, defString);
   return rcd ? 1:0;
}

////////////////////////////////////////////////////////////////////////////////
/// Take the scopes which enclose the decl and put them around the
/// definition string.
/// If a class is encountered, bail out.

const clang::RecordDecl* ROOT::TMetaUtils::AST2SourceTools::EncloseInScopes(const clang::Decl& decl,
                                                                            std::string& defString)
{
   std::list<std::pair<std::string,unsigned int> > enclosingNamespaces;
   auto rcdPtr = ROOT::TMetaUtils::ExtractEnclosingScopes(decl,enclosingNamespaces);

   if (rcdPtr) return rcdPtr;

   // Check if we have enclosing namespaces
   static const std::string scopeType [] = {"namespace ", "inline namespace ", "class "};

   std::string scopeName;
   std::string scopeContent;
   unsigned int scopeIndex;
   for (auto const & encScope : enclosingNamespaces){
      scopeIndex = encScope.second;
      scopeName = encScope.first;
      scopeContent = " { " + defString + " }";
      defString = scopeType[scopeIndex] +
                  scopeName +
                  scopeContent;
   }
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Loop over the template parameters and build a string for template arguments
/// using the fully qualified name
/// There are different cases:
/// Case 1: a simple template parameter
///   E.g. template<typename T> class A;
/// Case 2: a non-type: either an integer or an enum
///   E.g. template<int I, Foo > class A; where Foo is enum Foo {red, blue};
/// 2 sub cases here:
///   SubCase 2.a: the parameter is an enum: bail out, cannot be treated.
///   SubCase 2.b: use the fully qualified name
/// Case 3: a TemplateTemplate argument
///   E.g. template <template <typename> class T> class container { };

int ROOT::TMetaUtils::AST2SourceTools::PrepareArgsForFwdDecl(std::string& templateArgs,
                          const clang::TemplateParameterList& tmplParamList,
                          const cling::Interpreter& interpreter)
{
   templateArgs="<";
   for (auto prmIt = tmplParamList.begin();
        prmIt != tmplParamList.end(); prmIt++){

      if (prmIt != tmplParamList.begin())
         templateArgs += ", ";

      auto nDecl = *prmIt;
      std::string typeName;

      // Case 1
      if (llvm::isa<clang::TemplateTypeParmDecl>(nDecl)){
         typeName = "typename ";
         if (nDecl->isParameterPack())
            typeName += "... ";
         typeName += (*prmIt)->getNameAsString();
      }
      // Case 2
      else if (auto nttpd = llvm::dyn_cast<clang::NonTypeTemplateParmDecl>(nDecl)){
         auto theType = nttpd->getType();
         // If this is an enum, use int as it is impossible to fwd declare and
         // this makes sense since it is not a type...
         if (theType.getAsString().find("enum") != std::string::npos){
            std::string astDump;
            llvm::raw_string_ostream ostream(astDump);
            nttpd->dump(ostream);
            ostream.flush();
            ROOT::TMetaUtils::Warning(0,"Forward declarations of templates with enums as template parameters. The responsible class is: %s\n", astDump.c_str());
            return 1;
         } else {
            ROOT::TMetaUtils::GetFullyQualifiedTypeName(typeName,
                                                        theType,
                                                        interpreter);
         }
      }
      // Case 3: TemplateTemplate argument
      else if (auto ttpd = llvm::dyn_cast<clang::TemplateTemplateParmDecl>(nDecl)){
         int retCode = FwdDeclFromTmplDecl(*ttpd,interpreter,typeName);
         if (retCode!=0){
            std::string astDump;
            llvm::raw_string_ostream ostream(astDump);
            ttpd->dump(ostream);
            ostream.flush();
            ROOT::TMetaUtils::Error(0,"Cannot reconstruct template template parameter forward declaration for %s\n", astDump.c_str());
            return 1;
         }
      }

      templateArgs += typeName;
   }

   templateArgs+=">";
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Convert a tmplt decl to its fwd decl

int ROOT::TMetaUtils::AST2SourceTools::FwdDeclFromTmplDecl(const clang::TemplateDecl& templDecl,
                                                           const cling::Interpreter& interpreter,
                                                           std::string& defString)
{
   std::string templatePrefixString;
   auto tmplParamList= templDecl.getTemplateParameters();
   if (!tmplParamList){ // Should never happen
      Error(0,
            "Cannot extract template parameter list for %s",
            templDecl.getNameAsString().c_str());
      return 1;
   }

   int retCode = PrepareArgsForFwdDecl(templatePrefixString,*tmplParamList,interpreter);
   if (retCode!=0){
      Warning(0,
               "Problems with arguments for forward declaration of class %s\n",
               templDecl.getNameAsString().c_str());
      return retCode;
   }
   templatePrefixString = "template " + templatePrefixString + " ";

   defString = templatePrefixString + "class ";
   if (templDecl.isParameterPack())
      defString += "... ";
   defString +=  templDecl.getNameAsString();
   if (llvm::isa<clang::TemplateTemplateParmDecl>(&templDecl)) {
      // When fwd declaring the template template arg of
      //   namespace N { template <template <class T> class C> class X; }
      // we don't need to put it into any namespace, and we want no trailing
      // ';'
      return 0;
   }
   defString += ';';
   return EncloseInNamespaces(templDecl, defString);
}

////////////////////////////////////////////////////////////////////////////////

static int TreatSingleTemplateArg(const clang::TemplateArgument& arg,
                                  std::string& argFwdDecl,
                                  const  cling::Interpreter& interpreter,
                                  bool acceptStl=false)
{
   using namespace ROOT::TMetaUtils::AST2SourceTools;

   // We do nothing in presence of ints, bools, templates.
   // We should probably in presence of templates though...
   if (clang::TemplateArgument::Type != arg.getKind()) return 0;

   auto argQualType = arg.getAsType();

   // Recursively remove all *
   while (llvm::isa<clang::PointerType>(argQualType.getTypePtr())) argQualType = argQualType->getPointeeType();

   auto argTypePtr = argQualType.getTypePtr();

   // Bail out on enums
   if (llvm::isa<clang::EnumType>(argTypePtr)){
      return 1;
   }

   // If this is a built-in, just return: fwd decl not necessary.
   if (llvm::isa<clang::BuiltinType>(argTypePtr)){
      return 0;
   }

   // Treat typedefs which are arguments
   if (auto tdTypePtr = llvm::dyn_cast<clang::TypedefType>(argTypePtr)) {
      FwdDeclFromTypeDefNameDecl(*tdTypePtr->getDecl(), interpreter, argFwdDecl);
      return 0;
   }

   if (auto argRecTypePtr = llvm::dyn_cast<clang::RecordType>(argTypePtr)){
      // Now we cannot but have a RecordType
      if (auto argRecDeclPtr = argRecTypePtr->getDecl()){
         FwdDeclFromRcdDecl(*argRecDeclPtr,interpreter,argFwdDecl,acceptStl);
      }
      return 0;
   }

   return 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Convert a rcd decl to its fwd decl
/// If this is a template specialisation, treat in the proper way.
/// If it is contained in a class, just fwd declare the class.

int ROOT::TMetaUtils::AST2SourceTools::FwdDeclFromRcdDecl(const clang::RecordDecl& recordDecl,
                                                          const cling::Interpreter& interpreter,
                                                          std::string& defString,
                                                          bool acceptStl)
{
   // Do not fwd declare the templates in the stl.
   if (ROOT::TMetaUtils::IsStdClass(recordDecl) && !acceptStl)
      return 0;

   // Do not fwd declare unnamed decls.
   if (!recordDecl.getIdentifier())
      return 0;

   // We may need to fwd declare the arguments of the template
   std::string argsFwdDecl;

   if (auto tmplSpecDeclPtr = llvm::dyn_cast<clang::ClassTemplateSpecializationDecl>(&recordDecl)){
      std::string argFwdDecl;
      if (GetErrorIgnoreLevel() == ROOT::TMetaUtils::kInfo)
         std::cout << "Class " << recordDecl.getNameAsString()
                   << " is a template specialisation. Treating its arguments.\n";
      for(auto arg : tmplSpecDeclPtr->getTemplateArgs().asArray()){
         int retCode = TreatSingleTemplateArg(arg, argFwdDecl, interpreter, acceptStl);
         if (GetErrorIgnoreLevel() == ROOT::TMetaUtils::kInfo){
            std::cout << " o Template argument ";
            if (retCode==0){
               std::cout << "successfully treated. Arg fwd decl: " << argFwdDecl << std::endl;
            } else {
               std::cout << "could not be treated. Abort fwd declaration generation.\n";
            }
         }

         if (retCode!=0){ // A sign we must bail out
            return retCode;
         }
         argsFwdDecl+=argFwdDecl;
      }

      if (acceptStl){
         defString=argsFwdDecl;
         return 0;
      }

      int retCode=0;
      if (auto tmplDeclPtr = tmplSpecDeclPtr->getSpecializedTemplate()){
         retCode = FwdDeclFromTmplDecl(*tmplDeclPtr,interpreter,defString);
      }
      defString = argsFwdDecl + "\n" + defString;
      return retCode;

   }

   defString = "class " + recordDecl.getNameAsString() + ";";
   const clang::RecordDecl* rcd =  EncloseInScopes(recordDecl, defString);

   if (rcd){
      FwdDeclFromRcdDecl(*rcd, interpreter,defString);
   }
   // Add a \n here to avoid long lines which contain duplications, for example (from MathCore):
   // namespace ROOT { namespace Math { class IBaseFunctionMultiDim; } }namespace ROOT { namespace Fit { template <typename FunType> class Chi2FCN; } }
   // namespace ROOT { namespace Math { class IGradientFunctionMultiDim; } }namespace ROOT { namespace Fit { template <typename FunType> class Chi2FCN; } }
   defString = argsFwdDecl + "\n" + defString;

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Extract "forward declaration" of a typedef.
/// If the typedef is contained in a class, just fwd declare the class.
/// If not, fwd declare the typedef and all the dependent typedefs and types if necessary.

int ROOT::TMetaUtils::AST2SourceTools::FwdDeclFromTypeDefNameDecl(const clang::TypedefNameDecl& tdnDecl,
                                                                  const cling::Interpreter& interpreter,
                                                                  std::string& fwdDeclString,
                                                                  std::unordered_set<std::string>* fwdDeclSetPtr)
{
   std::string buffer = tdnDecl.getNameAsString();
   std::string underlyingName;
   auto underlyingType = tdnDecl.getUnderlyingType().getCanonicalType();
   if (const clang::TagType* TT
       = llvm::dyn_cast<clang::TagType>(underlyingType.getTypePtr())) {
      if (clang::NamedDecl* ND = TT->getDecl()) {
         if (!ND->getIdentifier()) {
            // No fwd decl for unnamed underlying entities.
            return 0;
         }
      }
   }

   TNormalizedCtxt nCtxt(interpreter.getLookupHelper());
   ROOT::TMetaUtils::GetNormalizedName(underlyingName,
                                       underlyingType,
                                       interpreter,
                                       nCtxt);

   // Heuristic: avoid entities like myclass<myType1, myType2::xyz>
   if (underlyingName.find(">::") != std::string::npos)
      return 0;

   buffer="typedef "+underlyingName+" "+buffer+";";
   const clang::RecordDecl* rcd=EncloseInScopes(tdnDecl,buffer);
   if (rcd) {
      // We do not need the whole series of scopes, just the class.
      // It is enough to trigger an uncomplete type autoload/parse callback
      // for example: MyClass::blabla::otherNs::myTypedef
      return FwdDeclFromRcdDecl(*rcd, interpreter,fwdDeclString,fwdDeclSetPtr);
   }

   // Start Recursion if the underlying type is a TypedefNameDecl
   // Note: the simple cast w/o the getSingleStepDesugaredType call
   // does not work in case the typedef is in a namespace.
   auto& ctxt = tdnDecl.getASTContext();
   auto immediatelyUnderlyingType = underlyingType.getSingleStepDesugaredType(ctxt);

   if (auto underlyingTdnTypePtr = llvm::dyn_cast<clang::TypedefType>(immediatelyUnderlyingType.getTypePtr())){
      std::string tdnFwdDecl;
      auto underlyingTdnDeclPtr = underlyingTdnTypePtr->getDecl();
      FwdDeclFromTypeDefNameDecl(*underlyingTdnDeclPtr,
                                 interpreter,
                                 tdnFwdDecl,
                                 fwdDeclSetPtr);
      if (!fwdDeclSetPtr || fwdDeclSetPtr->insert(tdnFwdDecl).second)
         fwdDeclString+=tdnFwdDecl;
   } else if (auto CXXRcdDeclPtr = immediatelyUnderlyingType->getAsCXXRecordDecl()){
      std::string classFwdDecl;
      if (GetErrorIgnoreLevel() == ROOT::TMetaUtils::kInfo)
         std::cout << "Typedef " << tdnDecl.getNameAsString() << " hides a class: "
                   << CXXRcdDeclPtr->getNameAsString() << std::endl;
      int retCode = FwdDeclFromRcdDecl(*CXXRcdDeclPtr,
                                       interpreter,
                                       classFwdDecl,
                                       true /* acceptStl*/);
      if (retCode!=0){ // bail out
         return 0;
      }

      if (!fwdDeclSetPtr || fwdDeclSetPtr->insert(classFwdDecl).second)
         fwdDeclString+=classFwdDecl;
   }

   fwdDeclString+=buffer;

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Get the default value as string.
/// Limited at the moment to:
/// - Integers
/// - Booleans

int ROOT::TMetaUtils::AST2SourceTools::GetDefArg(const clang::ParmVarDecl& par,
                                                 std::string& valAsString,
                                                 const clang::PrintingPolicy& ppolicy)
{
   auto defArgExprPtr = par.getDefaultArg();
   auto& ctxt = par.getASTContext();
   if(!defArgExprPtr->isEvaluatable(ctxt)){
      return -1;
   }

   auto defArgType = par.getType();

   // The value is a boolean
   if (defArgType->isBooleanType()){
      bool result;
      defArgExprPtr->EvaluateAsBooleanCondition (result,ctxt);
      valAsString=std::to_string(result);
      return 0;
   }

   // The value is an integer
   if (defArgType->isIntegerType()){
      llvm::APSInt result;
      defArgExprPtr->EvaluateAsInt(result,ctxt);
      auto uintVal = *result.getRawData();
      if (result.isNegative()){
         long long int intVal=uintVal*-1;
         valAsString=std::to_string(intVal);
      } else {
         valAsString=std::to_string(uintVal);
      }

      return 0;
   }

   // The value is something else. We go for the generalised printer
   llvm::raw_string_ostream rso(valAsString);
   defArgExprPtr->printPretty(rso,nullptr,ppolicy);
   valAsString = rso.str();
   // We can be in presence of a string. Let's escape the characters properly.
   ROOT::TMetaUtils::ReplaceAll(valAsString,"\\\"","__TEMP__VAL__");
   ROOT::TMetaUtils::ReplaceAll(valAsString,"\"","\\\"");
   ROOT::TMetaUtils::ReplaceAll(valAsString,"__TEMP__VAL__","\\\"");

   return 0;
}

