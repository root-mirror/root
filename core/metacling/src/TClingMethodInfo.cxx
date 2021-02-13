// @(#)root/core/meta:$Id$
// Author: Paul Russo   30/07/2012

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TClingMethodInfo
Emulation of the CINT MethodInfo class.

The CINT C++ interpreter provides an interface to metadata about
a function through the MethodInfo class.  This class provides the
same functionality, using an interface as close as possible to
MethodInfo but the typedef metadata comes from the Clang C++
compiler, not CINT.
*/

#include "TClingMethodInfo.h"

#include "TClingCallFunc.h"
#include "TClingClassInfo.h"
#include "TClingMemberIter.h"
#include "TClingMethodArgInfo.h"
#include "TDictionary.h"
#include "TClingTypeInfo.h"
#include "TError.h"
#include "TClingUtils.h"
#include "TCling.h"
#include "ThreadLocalStorage.h"

#include "cling/Interpreter/Interpreter.h"
#include "cling/Interpreter/LookupHelper.h"
#include "cling/Utils/AST.h"

#include "clang/AST/ASTContext.h"
#include "clang/AST/CXXInheritance.h"
#include "clang/AST/Decl.h"
#include "clang/AST/DeclBase.h"
#include "clang/AST/DeclCXX.h"
#include "clang/AST/DeclTemplate.h"
#include "clang/AST/ExprCXX.h"
#include "clang/AST/GlobalDecl.h"
#include "clang/AST/Mangle.h"
#include "clang/AST/PrettyPrinter.h"
#include "clang/AST/Type.h"
#include "clang/Basic/IdentifierTable.h"
#include "clang/Sema/Lookup.h"
#include "clang/Sema/Sema.h"
#include "clang/Sema/Template.h"
#include "clang/Sema/TemplateDeduction.h"

#include "llvm/Support/Casting.h"
#include "llvm/Support/raw_ostream.h"

#include <algorithm>
#include <string>

using namespace clang;

TClingCXXRecMethIter::SpecFuncIter::SpecFuncIter(cling::Interpreter *interp, clang::DeclContext *DC,
                                                 llvm::SmallVectorImpl<clang::CXXMethodDecl *> &&specFuncs)
{
   auto *CXXRD = llvm::dyn_cast<CXXRecordDecl>(DC);
   if (!CXXRD)
      return;

   // Could trigger deserialization of decls.
   cling::Interpreter::PushTransactionRAII RAII(interp);

   auto emplaceSpecFunIfNeeded = [&](clang::CXXMethodDecl *D) {
      if (!D)
         return; // Handle "structor not found" case.

      if (std::find(CXXRD->decls_begin(), CXXRD->decls_end(), D) == CXXRD->decls_end()) {
         fDefDataSpecFuns.emplace_back(D);
      }
   };

   for (auto SpecFunc : specFuncs)
      emplaceSpecFunIfNeeded(SpecFunc);
}

bool TClingCXXRecMethIter::ShouldSkip(const clang::Decl *D) const
{
   if (const auto *FD = llvm::dyn_cast<clang::FunctionDecl>(D)) {
      if (FD->isDeleted())
         return true;
      if (const auto *RD = llvm::dyn_cast<clang::RecordDecl>(FD->getDeclContext())) {
         if (const auto *CXXMD = llvm::dyn_cast<clang::CXXMethodDecl>(FD)) {
            if (RD->isAnonymousStructOrUnion() &&
                GetInterpreter()->getSema().getSpecialMember(CXXMD) != clang::Sema::CXXInvalid) {
               // Do not enumerate special members of anonymous structs.
               return true;
            }
         }
      }
      return false;
   }
   return true;
}

bool TClingCXXRecMethIter::ShouldSkip(const clang::UsingShadowDecl *USD) const
{
   if (auto *FD = llvm::dyn_cast<clang::FunctionDecl>(USD->getTargetDecl())) {
      if (const auto *CXXMD = llvm::dyn_cast<clang::CXXMethodDecl>(FD)) {
         auto SpecMemKind = GetInterpreter()->getSema().getSpecialMember(CXXMD);
         if ((SpecMemKind == clang::Sema::CXXDefaultConstructor && CXXMD->getNumParams() == 0) ||
             ((SpecMemKind == clang::Sema::CXXCopyConstructor || SpecMemKind == clang::Sema::CXXMoveConstructor) &&
              CXXMD->getNumParams() == 1)) {
            // This is a special member pulled in through a using decl. Special
            // members of derived classes cannot be replaced; ignore this using decl,
            // and keep only the (still possibly compiler-generated) special member of the
            // derived class.
            // NOTE that e.g. `Klass(int = 0)` has SpecMemKind == clang::Sema::CXXDefaultConstructor,
            // yet this signature must be exposed, so check the argument count.
            return true;
         }
      }
      return ShouldSkip(FD);
   }
   // TODO: handle multi-level UsingShadowDecls.
   return true;
}

const clang::Decl *
TClingCXXRecMethIter::InstantiateTemplateWithDefaults(const clang::RedeclarableTemplateDecl *TD) const
{
   // Force instantiation if it doesn't exist yet, by looking it up.

   using namespace clang;

   cling::Interpreter *interp = GetInterpreter();
   Sema &S = interp->getSema();
   const cling::LookupHelper &LH = interp->getLookupHelper();

   if (!isa<FunctionTemplateDecl>(TD))
      return nullptr;

   auto templateParms = TD->getTemplateParameters();
   if (templateParms->containsUnexpandedParameterPack())
      return nullptr;

   if (templateParms->getMinRequiredArguments() > 0)
      return nullptr;

   if (templateParms->size() > 0) {
      NamedDecl *arg0 = *templateParms->begin();
      if (arg0->isTemplateParameterPack())
         return nullptr;
      if (auto TTP = dyn_cast<TemplateTypeParmDecl>(*templateParms->begin())) {
         if (!TTP->hasDefaultArgument())
            return nullptr;
      } else if (auto NTTP = dyn_cast<NonTypeTemplateParmDecl>(*templateParms->begin())) {
         if (!NTTP->hasDefaultArgument())
            return nullptr;
      } else {
         // TemplateTemplateParmDecl, pack
         return nullptr;
      }
   }

   const FunctionDecl *templatedDecl = llvm::dyn_cast<FunctionDecl>(TD->getTemplatedDecl());
   const Decl *declCtxDecl = dyn_cast<Decl>(TD->getDeclContext());

   // We have a function template
   //     template <class X = int, int i = 7> void func(int a0, X a1[i], X::type a2[i])
   // which has defaults for all its template parameters `X` and `i`. To
   // instantiate it we have to do a lookup, which in turn needs the function
   // argument types, e.g. `int[12]`.
   // If the function argument type is dependent (a1 and a2) we need to
   // substitute the types first, using the template arguments derived from the
   // template parameters' defaults.
   llvm::SmallVector<TemplateArgument, 8> defaultTemplateArgs(templateParms->size());
   for (int iParam = 0, nParams = templateParms->size(); iParam < nParams; ++iParam) {
      const NamedDecl *templateParm = templateParms->getParam(iParam);
      if (templateParm->isTemplateParameterPack()) {
         // shouldn't end up here
         assert(0 && "unexpected template parameter pack");
         return nullptr;
      }
      if (auto TTP = dyn_cast<TemplateTypeParmDecl>(templateParm)) {
         if (!TTP->hasDefaultArgument())
            return nullptr;
         defaultTemplateArgs[iParam] = TemplateArgument(TTP->getDefaultArgument());
      } else if (auto NTTP = dyn_cast<NonTypeTemplateParmDecl>(templateParm)) {
         if (!NTTP->hasDefaultArgument())
            return nullptr;
         defaultTemplateArgs[iParam] = TemplateArgument(NTTP->getDefaultArgument());
      } else if (auto TTP = dyn_cast<TemplateTemplateParmDecl>(templateParm)) {
         if (!TTP->hasDefaultArgument())
            return nullptr;
         defaultTemplateArgs[iParam] = TemplateArgument(TTP->getDefaultArgument().getArgument());
      } else {
         // shouldn't end up here
         assert(0 && "unexpected template parameter kind");
         return nullptr;
      }
   }

   cling::Interpreter::PushTransactionRAII RAII(interp);

   // Now substitute the dependent function parameter types given defaultTemplateArgs.
   llvm::SmallVector<QualType, 8> paramTypes;
   // Provide an instantiation context that suppresses errors:
   // DeducedTemplateArgumentSubstitution! (ROOT-8422)
   SmallVector<DeducedTemplateArgument, 4> DeducedArgs;
   sema::TemplateDeductionInfo Info{SourceLocation()};

   Sema::InstantiatingTemplate Inst(
      S, Info.getLocation(), const_cast<clang::FunctionTemplateDecl *>(llvm::dyn_cast<clang::FunctionTemplateDecl>(TD)),
      defaultTemplateArgs, Sema::CodeSynthesisContext::DeducedTemplateArgumentSubstitution, Info);

   // Collect the function arguments of the templated function, substituting
   // dependent types as possible.
   TemplateArgumentList templArgList(TemplateArgumentList::OnStack, defaultTemplateArgs);
   MultiLevelTemplateArgumentList MLTAL{templArgList};
   for (const clang::ParmVarDecl *param : templatedDecl->parameters()) {
      QualType paramType = param->getOriginalType();

      // If the function type is dependent, try to resolve it through the class's
      // template arguments. If that fails, skip this function.
      if (paramType->isDependentType()) {
         /*if (HasUnexpandedParameterPack(paramType, S)) {
            // We are not going to expand the pack here...
            Skip = true;
            break;
         }*/

         paramType = S.SubstType(paramType, MLTAL, SourceLocation(), templatedDecl->getDeclName());

         if (paramType.isNull() || paramType->isDependentType()) {
            // Even after resolving the types through the surrounding template
            // this argument type is still dependent: do not look it up.
            return nullptr;
         }
      }
      paramTypes.push_back(paramType);
   }

   return LH.findFunctionProto(declCtxDecl, TD->getNameAsString(), paramTypes, LH.NoDiagnostics,
                               templatedDecl->getType().isConstQualified());
}

TClingMethodInfo::TClingMethodInfo(cling::Interpreter *interp,
                                   TClingClassInfo *ci)
   : TClingDeclInfo(nullptr), fInterp(interp), fFirstTime(true), fTitle("")
{
   R__LOCKGUARD(gInterpreterMutex);

   if (!ci || !ci->IsValid()) {
      return;
   }
   clang::Decl *D = const_cast<clang::Decl *>(ci->GetDecl());
   auto *DC = llvm::dyn_cast<clang::DeclContext>(D);

   llvm::SmallVector<clang::CXXMethodDecl*, 8> SpecFuncs;

   if (auto *CXXRD = llvm::dyn_cast<CXXRecordDecl>(DC)) {
      // Initialize the CXXRecordDecl's special functions; could change the
      // DeclContext content!

      // Could trigger deserialization of decls.
      cling::Interpreter::PushTransactionRAII RAII(interp);

      auto &SemaRef = interp->getSema();
      SemaRef.ForceDeclarationOfImplicitMembers(CXXRD);

      // Assemble special functions (or FunctionTemplate-s) that are synthesized from DefinitionData but
      // won't be enumerated as part of decls_begin()/decls_end().
      for (clang::NamedDecl *ctor : SemaRef.LookupConstructors(CXXRD)) {
         // Filter out constructor templates, they are not functions we can iterate over:
         if (auto *CXXCD = llvm::dyn_cast<clang::CXXConstructorDecl>(ctor))
            SpecFuncs.emplace_back(CXXCD);
      }
      SpecFuncs.emplace_back(SemaRef.LookupCopyingAssignment(CXXRD, /*Quals*/ 0, /*RValueThis*/ false, 0 /*ThisQuals*/));
      SpecFuncs.emplace_back(SemaRef.LookupMovingAssignment(CXXRD, /*Quals*/ 0, /*RValueThis*/ false, 0 /*ThisQuals*/));
      SpecFuncs.emplace_back(SemaRef.LookupDestructor(CXXRD));
   }

   fIter = TClingCXXRecMethIter(interp, DC, std::move(SpecFuncs));
   fIter.Init();
}

TClingMethodInfo::TClingMethodInfo(cling::Interpreter *interp,
                                   const clang::Decl *D)
   : TClingDeclInfo(D), fInterp(interp), fFirstTime(true), fTitle("")
{
   if (!D)
      Error("TClingMethodInfo", "nullptr FunctionDecl passed!");
}

TDictionary::DeclId_t TClingMethodInfo::GetDeclId() const
{
   if (!IsValid()) {
      return TDictionary::DeclId_t();
   }
   if (auto *FD = GetAsFunctionDecl())
      return (const clang::Decl*)(FD->getCanonicalDecl());
   return (const clang::Decl*)(GetAsUsingShadowDecl()->getCanonicalDecl());
}

const clang::FunctionDecl *TClingMethodInfo::GetAsFunctionDecl() const
{
   return dyn_cast<FunctionDecl>(GetDecl());
}

const clang::UsingShadowDecl *TClingMethodInfo::GetAsUsingShadowDecl() const
{
   return dyn_cast<UsingShadowDecl>(GetDecl());
}

const clang::FunctionDecl *TClingMethodInfo::GetTargetFunctionDecl() const
{
   const Decl *D = GetDecl();
   do {
      if (auto FD = dyn_cast<FunctionDecl>(D))
         return FD;
   } while ((D = dyn_cast<UsingShadowDecl>(D)->getTargetDecl()));
   return nullptr;
}

void TClingMethodInfo::CreateSignature(TString &signature) const
{
   signature = "(";
   if (!IsValid()) {
      signature += ")";
      return;
   }

   R__LOCKGUARD(gInterpreterMutex);
   TClingMethodArgInfo arg(fInterp, this);

   int idx = 0;
   while (arg.Next()) {
      if (idx) {
         signature += ", ";
      }
      signature += arg.Type()->Name();
      if (arg.Name() && strlen(arg.Name())) {
         signature += " ";
         signature += arg.Name();
      }
      if (arg.DefaultValue()) {
         signature += " = ";
         signature += arg.DefaultValue();
      }
      ++idx;
   }
   auto decl = GetTargetFunctionDecl();
   if (decl && decl->isVariadic())
      signature += ",...";

   signature += ")";
}

void TClingMethodInfo::Init(const clang::FunctionDecl *decl)
{
   fFirstTime = true;
   fIter = {};
   fDecl = decl;
}

void *TClingMethodInfo::InterfaceMethod(const ROOT::TMetaUtils::TNormalizedCtxt &normCtxt) const
{
   if (!IsValid()) {
      return 0;
   }
   R__LOCKGUARD(gInterpreterMutex);
   TClingCallFunc cf(fInterp,normCtxt);
   cf.SetFunc(this);
   return cf.InterfaceMethod();
}

const clang::Decl* TClingMethodInfo::GetDeclSlow() const
{
   return *fIter;
}

int TClingMethodInfo::NArg() const
{
   if (!IsValid()) {
      return -1;
   }
   const clang::FunctionDecl *fd = GetTargetFunctionDecl();
   unsigned num_params = fd->getNumParams();
   // Truncate cast to fit cint interface.
   return static_cast<int>(num_params);
}

int TClingMethodInfo::NDefaultArg() const
{
   if (!IsValid()) {
      return -1;
   }
   const clang::FunctionDecl *fd = GetTargetFunctionDecl();
   unsigned num_params = fd->getNumParams();
   unsigned min_args = fd->getMinRequiredArguments();
   unsigned defaulted_params = num_params - min_args;
   // Truncate cast to fit cint interface.
   return static_cast<int>(defaulted_params);
}

/*
static bool HasUnexpandedParameterPack(clang::QualType QT, clang::Sema& S) {
   if (llvm::isa<PackExpansionType>(*QT)) {
      // We are not going to expand the pack here...
      return true;
   }
   SmallVector<UnexpandedParameterPack, 4> Unexpanded;
   S.collectUnexpandedParameterPacks (QT, Unexpanded);

   return !Unexpanded.empty();
}
 */

int TClingMethodInfo::Next()
{

   assert(!fDecl && "This is not an iterator!");

   fNameCache.clear(); // invalidate the cache.

   if (!fFirstTime && !fIter.IsValid()) {
      // Iterator is already invalid.
      return 0;
   }
   // Advance to the next decl.
   if (fFirstTime) {
      // The cint semantics are weird.
      fFirstTime = false;
   } else {
      fIter.Next();
   }
   return fIter.IsValid();
}

long TClingMethodInfo::Property() const
{
   if (!IsValid()) {
      return 0L;
   }
   long property = 0L;
   property |= kIsCompiled;

   // NOTE: this uses `GetDecl()`, to capture the access of the UsingShadowDecl,
   // which is defined in the derived class and might differ from the access of fd
   // in the base class.
   const Decl *declAccess = GetDecl();
   if (llvm::isa<UsingShadowDecl>(declAccess))
      property |= kIsUsing;

   const clang::FunctionDecl *fd = GetTargetFunctionDecl();
   clang::AccessSpecifier Access = clang::AS_public;
   if (!declAccess->getDeclContext()->isNamespace())
      Access = declAccess->getAccess();

   if ((property & kIsUsing) && llvm::isa<CXXConstructorDecl>(fd)) {
      Access = clang::AS_public;
      clang::CXXRecordDecl *typeCXXRD = llvm::cast<RecordType>(Type()->GetQualType())->getAsCXXRecordDecl();
      clang::CXXBasePaths basePaths;
      if (typeCXXRD->isDerivedFrom(llvm::dyn_cast<CXXRecordDecl>(fd->getDeclContext()), basePaths)) {
         // Access of the ctor is access of the base inheritance, and
         // cannot be overruled by the access of the using decl.

         for (auto el: basePaths) {
            if (el.Access > Access)
               Access = el.Access;
         }
      } else {
         Error("Property()", "UsingDecl of ctor not shadowing a base ctor!");
      }

      // But a private ctor stays private:
      if (fd->getAccess() > Access)
         Access = fd->getAccess();
   }
   switch (Access) {
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
         if (declAccess->getDeclContext()->isNamespace())
            property |= kIsPublic;
         break;
      default:
         // IMPOSSIBLE
         break;
   }

   if (fd->isConstexpr())
      property |= kIsConstexpr;
   if (fd->getStorageClass() == clang::SC_Static) {
      property |= kIsStatic;
   }
   clang::QualType qt = fd->getReturnType().getCanonicalType();

   property = TClingDeclInfo::Property(property, qt);

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

long TClingMethodInfo::ExtraProperty() const
{
   // Return the property not already defined in Property
   // See TDictionary's EFunctionProperty
   if (!IsValid()) {
      return 0L;
   }
   long property = 0;
   const clang::FunctionDecl *fd = GetTargetFunctionDecl();
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
   if (fd->getTemplatedKind() != clang::FunctionDecl::TK_NonTemplate)
      property |= kIsTemplateSpec;
   return property;
}

TClingTypeInfo *TClingMethodInfo::Type() const
{
   TTHREAD_TLS_DECL_ARG( TClingTypeInfo, ti, fInterp );
   if (!IsValid()) {
      ti.Init(clang::QualType());
      return &ti;
   }
   if (llvm::isa<clang::CXXConstructorDecl>(GetTargetFunctionDecl())) {
      // CINT claims that constructors return the class object.
      // For using-ctors of a base, claim that it "returns" the derived class.
      const clang::TypeDecl* ctorClass = llvm::dyn_cast_or_null<clang::TypeDecl>
         (GetDecl()->getDeclContext());
      if (!ctorClass) {
         Error("TClingMethodInfo::Type", "Cannot find DeclContext for constructor!");
      } else {
         clang::QualType qt(ctorClass->getTypeForDecl(), 0);
         ti.Init(qt);
      }
   } else {
      clang::QualType qt = GetTargetFunctionDecl()->getReturnType();
      ti.Init(qt);
   }
   return &ti;
}

std::string TClingMethodInfo::GetMangledName() const
{
   if (!IsValid()) {
      return "";
   }
   std::string mangled_name;
   mangled_name.clear();
   const FunctionDecl* D = GetTargetFunctionDecl();

   R__LOCKGUARD(gInterpreterMutex);
   cling::Interpreter::PushTransactionRAII RAII(fInterp);
   GlobalDecl GD;
   if (const CXXConstructorDecl* Ctor = dyn_cast<CXXConstructorDecl>(D))
     GD = GlobalDecl(Ctor, Ctor_Complete);
   else if (const CXXDestructorDecl* Dtor = dyn_cast<CXXDestructorDecl>(D))
     GD = GlobalDecl(Dtor, Dtor_Deleting);
   else
     GD = GlobalDecl(D);

   cling::utils::Analyze::maybeMangleDeclName(GD, mangled_name);
   return mangled_name;
}

const char *TClingMethodInfo::GetPrototype()
{
   if (!IsValid()) {
      return 0;
   }
   TTHREAD_TLS_DECL( std::string, buf );
   buf.clear();
   buf += Type()->Name();
   buf += ' ';
   const FunctionDecl *FD = GetTargetFunctionDecl();
   // Use the DeclContext of the decl, not of the target decl:
   // Used base functions should show as if they are part of the derived class,
   // e.g. `Derived Derived::Derived(int)`, not `Derived Base::Derived(int)`.
   if (const clang::TypeDecl *td = llvm::dyn_cast<clang::TypeDecl>(GetDecl()->getDeclContext())) {
      std::string name;
      clang::QualType qualType(td->getTypeForDecl(),0);
      ROOT::TMetaUtils::GetFullyQualifiedTypeName(name,qualType,*fInterp);
      buf += name;
      buf += "::";
   } else if (const clang::NamedDecl *nd = llvm::dyn_cast<clang::NamedDecl>(FD->getDeclContext())) {
      std::string name;
      clang::PrintingPolicy policy(FD->getASTContext().getPrintingPolicy());
      llvm::raw_string_ostream stream(name);
      nd->getNameForDiagnostic(stream, policy, /*Qualified=*/true);
      stream.flush();
      buf += name;
      buf += "::";
   }
   buf += Name();

   TString signature;
   CreateSignature(signature);
   buf += signature;

   if (const clang::CXXMethodDecl *md =
       llvm::dyn_cast<clang::CXXMethodDecl>(FD)) {
      if (md->getTypeQualifiers() & clang::Qualifiers::Const) {
         buf += " const";
      }
   }
   return buf.c_str();
}

const char *TClingMethodInfo::Name() const
{
   if (!IsValid()) {
      return 0;
   }
   if (!fNameCache.empty())
     return fNameCache.c_str();

   ((TCling*)gCling)->GetFunctionName(GetDecl(), fNameCache);
   return fNameCache.c_str();
}

const char *TClingMethodInfo::TypeName() const
{
   if (!IsValid()) {
      // FIXME: Cint does not check!
      return 0;
   }
   return Type()->Name();
}

const char *TClingMethodInfo::Title()
{
   if (!IsValid()) {
      return 0;
   }

   //NOTE: We can't use it as a cache due to the "thoughtful" self iterator
   //if (fTitle.size())
   //   return fTitle.c_str();

   // Try to get the comment either from the annotation or the header file if present

   // Iterate over the redeclarations, we can have multiple definitions in the
   // redecl chain (came from merging of pcms).
   const FunctionDecl *FD = GetTargetFunctionDecl();

   R__LOCKGUARD(gInterpreterMutex);

   // Could trigger deserialization of decls.
   cling::Interpreter::PushTransactionRAII RAII(fInterp);
   if (const FunctionDecl *AnnotFD
       = ROOT::TMetaUtils::GetAnnotatedRedeclarable(FD)) {
      if (AnnotateAttr *A = AnnotFD->getAttr<AnnotateAttr>()) {
         fTitle = A->getAnnotation().str();
         return fTitle.c_str();
      }
   }
   if (!FD->isFromASTFile()) {
      // Try to get the comment from the header file if present
      // but not for decls from AST file, where rootcling would have
      // created an annotation
      fTitle = ROOT::TMetaUtils::GetComment(*FD).str();
   }

   return fTitle.c_str();
}

