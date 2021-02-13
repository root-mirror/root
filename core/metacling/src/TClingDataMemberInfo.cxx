// @(#)root/core/meta:$Id$
// Author: Paul Russo   30/07/2012

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TClingDataMemberInfo

Emulation of the CINT DataMemberInfo class.

The CINT C++ interpreter provides an interface to metadata about
the data members of a class through the DataMemberInfo class.  This
class provides the same functionality, using an interface as close
as possible to DataMemberInfo but the data member metadata comes
from the Clang C++ compiler, not CINT.
*/

#include "TClingDataMemberInfo.h"

#include "TDictionary.h"
#include "TClingClassInfo.h"
#include "TClingTypeInfo.h"
#include "TClingUtils.h"
#include "TClassEdit.h"
#include "TError.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"

#include "cling/Interpreter/Interpreter.h"

#include "clang/AST/Attr.h"
#include "clang/AST/ASTContext.h"
#include "clang/AST/Decl.h"
#include "clang/AST/GlobalDecl.h"
#include "clang/AST/Expr.h"
#include "clang/AST/ExprCXX.h"
#include "clang/AST/PrettyPrinter.h"
#include "clang/AST/RecordLayout.h"
#include "clang/AST/Type.h"

#include "llvm/Support/Casting.h"
#include "llvm/Support/raw_ostream.h"
#include "llvm/ADT/APSInt.h"
#include "llvm/ADT/APFloat.h"

using namespace clang;

namespace {
   static bool IsRelevantKind(clang::Decl::Kind DK)
   {
      return DK == clang::Decl::Field || DK == clang::Decl::EnumConstant || DK == clang::Decl::Var;
   }
}

bool TClingDataMemberIter::ShouldSkip(const clang::Decl *D) const
{
   if (!TDictionary::WantsRegularMembers(fSelection))
      return true;

   if (const auto *ND = llvm::dyn_cast<NamedDecl>(D)) {
      // Skip unnamed declarations, e.g. in
      //    struct S {
      //       struct { int i; }
      //    };
      // the inner struct corresponds to an unnamed member variable,
      // where only `S::i` should be exposed.
      if (!ND->getIdentifier())
         return true;
   } else {
      // TClingDataMemberIter only cares about NamedDecls.
      return true;
   }

   return !IsRelevantKind(D->getKind());
}

bool TClingDataMemberIter::ShouldSkip(const clang::UsingShadowDecl *USD) const
{
   if (!TDictionary::WantsUsingDecls(fSelection))
      return true;

   if (auto *VD = llvm::dyn_cast<clang::ValueDecl>(USD->getTargetDecl())) {
      return !IsRelevantKind(VD->getKind());
   }

   // TODO: handle multi-level UsingShadowDecls.
   return true;
}

TClingDataMemberInfo::TClingDataMemberInfo(cling::Interpreter *interp,
                                           TClingClassInfo *ci,
                                           TDictionary::EMemberSelection selection)
: TClingDeclInfo(nullptr), fInterp(interp), fClassInfo(ci ? new TClingClassInfo(*ci) : new TClingClassInfo(interp))
{

   R__LOCKGUARD(gInterpreterMutex);

   if (!ci || !ci->IsValid()) {
      return;
   }

   auto *DC = llvm::dyn_cast<clang::DeclContext>(ci->GetDecl());

   fIter = TClingDataMemberIter(interp, DC, selection);
   fIter.Init();
}

TClingDataMemberInfo::TClingDataMemberInfo(cling::Interpreter *interp,
                                           const clang::ValueDecl *ValD,
                                           TClingClassInfo *ci)
: TClingDeclInfo(ValD), fInterp(interp), fClassInfo(ci ? new TClingClassInfo(*ci) : new TClingClassInfo(interp))
{

   using namespace llvm;
   const auto DC = ValD->getDeclContext();
   (void)DC;
   assert((ci || isa<TranslationUnitDecl>(DC) ||
          ((DC->isTransparentContext() || DC->isInlineNamespace()) && isa<TranslationUnitDecl>(DC->getParent()) ) ||
           isa<EnumConstantDecl>(ValD)) && "Not TU?");
   assert(IsRelevantKind(ValD->getKind()) &&
          "The decl should be either VarDecl or FieldDecl or EnumConstDecl");

}

void TClingDataMemberInfo::CheckForIoTypeAndName() const
{
   // Three cases:
   // 1) 00: none to be checked
   // 2) 01: type to be checked
   // 3) 10: none to be checked
   // 4) 11: both to be checked
   unsigned int code = fIoType.empty() + (int(fIoName.empty()) << 1);

   if (code == 0) return;

   const Decl* decl = GetTargetValueDecl();

   if (code == 3 || code == 2) ROOT::TMetaUtils::ExtractAttrPropertyFromName(*decl,"ioname",fIoName);
   if (code == 3 || code == 1) ROOT::TMetaUtils::ExtractAttrPropertyFromName(*decl,"iotype",fIoType);

}

TDictionary::DeclId_t TClingDataMemberInfo::GetDeclId() const
{
   if (!IsValid()) {
      return TDictionary::DeclId_t();
   }
   if (auto *VD = GetAsValueDecl())
      return (const clang::Decl*)(VD->getCanonicalDecl());
   return (const clang::Decl*)(GetAsUsingShadowDecl()->getCanonicalDecl());
}

const clang::ValueDecl *TClingDataMemberInfo::GetAsValueDecl() const
{
   return dyn_cast<ValueDecl>(GetDecl());
}

const clang::UsingShadowDecl *TClingDataMemberInfo::GetAsUsingShadowDecl() const
{
   return dyn_cast<UsingShadowDecl>(GetDecl());
}

const clang::ValueDecl *TClingDataMemberInfo::GetTargetValueDecl() const
{
   const Decl *D = GetDecl();
   do {
      if (auto VD = dyn_cast<ValueDecl>(D))
         return VD;
   } while ((D = dyn_cast<UsingShadowDecl>(D)->getTargetDecl()));
   return nullptr;
}

const clang::Type *TClingDataMemberInfo::GetClassAsType() const {
   return fClassInfo ? fClassInfo->GetType() : nullptr;
}

int TClingDataMemberInfo::ArrayDim() const
{
   if (!IsValid()) {
      return -1;
   }
   const clang::ValueDecl *VD = GetTargetValueDecl();
   // Sanity check the current data member.
   clang::Decl::Kind DK = VD->getKind();
   if (
       (DK != clang::Decl::Field) &&
       (DK != clang::Decl::Var) &&
       (DK != clang::Decl::EnumConstant)
       ) {
      // Error, was not a data member, variable, or enumerator.
      return -1;
   }
   if (DK == clang::Decl::EnumConstant) {
      // We know that an enumerator value does not have array type.
      return 0;
   }
   // To get this information we must count the number
   // of array type nodes in the canonical type chain.
   clang::QualType QT = VD->getType().getCanonicalType();
   int cnt = 0;
   while (1) {
      if (QT->isArrayType()) {
         ++cnt;
         QT = llvm::cast<clang::ArrayType>(QT)->getElementType();
         continue;
      }
      else if (QT->isReferenceType()) {
         QT = llvm::cast<clang::ReferenceType>(QT)->getPointeeType();
         continue;
      }
      else if (QT->isPointerType()) {
         QT = llvm::cast<clang::PointerType>(QT)->getPointeeType();
         continue;
      }
      else if (QT->isMemberPointerType()) {
         QT = llvm::cast<clang::MemberPointerType>(QT)->getPointeeType();
         continue;
      }
      break;
   }
   return cnt;
}

int TClingDataMemberInfo::MaxIndex(int dim) const
{
   if (!IsValid()) {
      return -1;
   }
   const clang::ValueDecl *VD = GetTargetValueDecl();
   // Sanity check the current data member.
   clang::Decl::Kind DK = GetDecl()->getKind();
   if (
       (DK != clang::Decl::Field) &&
       (DK != clang::Decl::Var) &&
       (DK != clang::Decl::EnumConstant)
       ) {
      // Error, was not a data member, variable, or enumerator.
      return -1;
   }
   if (DK == clang::Decl::EnumConstant) {
      // We know that an enumerator value does not have array type.
      return 0;
   }
   // To get this information we must count the number
   // of array type nodes in the canonical type chain.
   clang::QualType QT = VD->getType().getCanonicalType();
   int paran = ArrayDim();
   if ((dim < 0) || (dim >= paran)) {
      // Passed dimension is out of bounds.
      return -1;
   }
   int cnt = dim;
   int max = 0;
   while (1) {
      if (QT->isArrayType()) {
         if (cnt == 0) {
            if (const clang::ConstantArrayType *CAT =
                llvm::dyn_cast<clang::ConstantArrayType>(QT)
                ) {
               max = static_cast<int>(CAT->getSize().getZExtValue());
            }
            else if (llvm::dyn_cast<clang::IncompleteArrayType>(QT)) {
               max = INT_MAX;
            }
            else {
               max = -1;
            }
            break;
         }
         --cnt;
         QT = llvm::cast<clang::ArrayType>(QT)->getElementType();
         continue;
      }
      else if (QT->isReferenceType()) {
         QT = llvm::cast<clang::ReferenceType>(QT)->getPointeeType();
         continue;
      }
      else if (QT->isPointerType()) {
         QT = llvm::cast<clang::PointerType>(QT)->getPointeeType();
         continue;
      }
      else if (QT->isMemberPointerType()) {
         QT = llvm::cast<clang::MemberPointerType>(QT)->getPointeeType();
         continue;
      }
      break;
   }
   return max;
}

int TClingDataMemberInfo::Next()
{
   assert(!fDecl && "This is a single decl, not an iterator!");

   ClearNames();

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

long TClingDataMemberInfo::Offset()
{
   using namespace clang;

   if (!IsValid()) {
      return -1L;
   }

   const ValueDecl *D = GetTargetValueDecl();
   ASTContext& C = D->getASTContext();
   if (const FieldDecl *FldD = dyn_cast<FieldDecl>(D)) {
      // The current member is a non-static data member.
      const clang::RecordDecl *RD = FldD->getParent();
      const clang::ASTRecordLayout &Layout = C.getASTRecordLayout(RD);
      uint64_t bits = Layout.getFieldOffset(FldD->getFieldIndex());
      int64_t offset = C.toCharUnitsFromBits(bits).getQuantity();
      return static_cast<long>(offset);
   }
   else if (const VarDecl *VD = dyn_cast<VarDecl>(D)) {
      // Could trigger deserialization of decls, in particular in case
      // of constexpr, like:
      //   static constexpr Long64_t something = std::numeric_limits<Long64_t>::max();
      cling::Interpreter::PushTransactionRAII RAII(fInterp);

      if (long addr = reinterpret_cast<long>(fInterp->getAddressOfGlobal(GlobalDecl(VD))))
         return addr;
      auto evalStmt = VD->ensureEvaluatedStmt();
      if (evalStmt && evalStmt->Value) {
         if (const APValue* val = VD->evaluateValue()) {
            if (VD->getType()->isIntegralType(C)) {
               return reinterpret_cast<long>(val->getInt().getRawData());
            } else {
               // The VD stores the init value; its lifetime should the lifetime of
               // this offset.
               switch (val->getKind()) {
               case APValue::Int: {
                  if (val->getInt().isSigned())
                     fConstInitVal.fLong = (long)val->getInt().getSExtValue();
                  else
                     fConstInitVal.fLong = (long)val->getInt().getZExtValue();
                  return (long) &fConstInitVal.fLong;
               }
               case APValue::Float:
                  if (&val->getFloat().getSemantics()
                      == (const llvm::fltSemantics*)&llvm::APFloat::IEEEsingle()) {
                     fConstInitVal.fFloat = val->getFloat().convertToFloat();
                     return (long)&fConstInitVal.fFloat;
                  } else if (&val->getFloat().getSemantics()
                             == (const llvm::fltSemantics*) &llvm::APFloat::IEEEdouble()) {
                     fConstInitVal.fDouble = val->getFloat().convertToDouble();
                     return (long)&fConstInitVal.fDouble;
                  }
                  // else fall-through
               default:
                  ;// fall-through
               };
               // fall-through
            } // not integral type
         } // have an APValue
      } // have an initializing value
   }
   // FIXME: We have to explicitly check for not enum constant because the
   // implementation of getAddressOfGlobal relies on mangling the name and in
   // clang there is misbehaviour in MangleContext::shouldMangleDeclName.
   // enum constants are essentially numbers and don't get addresses. However
   // ROOT expects the address to the enum constant initializer to be returned.
   else if (const EnumConstantDecl *ECD = dyn_cast<EnumConstantDecl>(D))
      // The raw data is stored as a long long, so we need to find the 'long'
      // part.
#ifdef R__BYTESWAP
      // In this case at the beginning.
      return reinterpret_cast<long>(ECD->getInitVal().getRawData());
#else
      // In this case in the second part.
      return reinterpret_cast<long>(((char*)ECD->getInitVal().getRawData())+sizeof(long) );
#endif
   return -1L;
}

long TClingDataMemberInfo::Property() const
{
   if (!IsValid()) {
      return 0L;
   }
   long property = 0L;

   // If the declaration is public in a private nested struct, make the declaration
   // private nonetheless, as for outside access (e.g. ROOT I/O) it's private:
   // NOTE: this uses `GetDecl()`, to capture the access of the UsingShadowDecl,
   // which is defined in the derived class and might differ from the access of the decl
   // in the base class.
   // TODO: move this somewhere such that TClingMethodInfo can use this, too.
   const Decl *thisDecl = GetDecl();
   clang::AccessSpecifier strictestAccess = thisDecl->getAccess();
   const DeclContext *nonTransparentDC = thisDecl->getDeclContext();

   auto getParentAccessAndNonTransparentDC = [&]() {
      const Decl *declOrParent = thisDecl;
      for (const auto *Parent = declOrParent->getDeclContext(); !llvm::isa<TranslationUnitDecl>(Parent);
           Parent = declOrParent->getDeclContext()) {
         if (!Parent->isTransparentContext()) {
            if (const auto *RD = llvm::dyn_cast<clang::RecordDecl>(Parent)) {
               if (!RD->isAnonymousStructOrUnion()) {
                  nonTransparentDC = RD;
                  break;
               }
            } else {
               nonTransparentDC = Parent;
               break;
            }
         }

         declOrParent = llvm::dyn_cast<clang::Decl>(Parent);
         if (strictestAccess < declOrParent->getAccess()) {
            strictestAccess = declOrParent->getAccess();
         }
      }
   };

   getParentAccessAndNonTransparentDC();

   switch (strictestAccess) {
      case clang::AS_public:
         property |= kIsPublic;
         break;
      case clang::AS_protected:
         property |= kIsProtected;
         break;
      case clang::AS_private:
         property |= kIsPrivate;
         break;
      case clang::AS_none: //?
         property |= kIsPublic;
         break;
      default:
         // IMPOSSIBLE
         break;
   }
   if (llvm::isa<clang::UsingShadowDecl>(thisDecl))
      property |= kIsUsing;

   const clang::ValueDecl *vd = GetTargetValueDecl();
   if (const clang::VarDecl *vard = llvm::dyn_cast<clang::VarDecl>(vd)) {
      if (vard->isConstexpr())
         property |= kIsConstexpr;
      if (vard->getStorageClass() == clang::SC_Static) {
         property |= kIsStatic;
      } else if (nonTransparentDC->isNamespace()) {
         // Data members of a namespace are global variable which were
         // considered to be 'static' in the CINT (and thus ROOT) scheme.
         property |= kIsStatic;
      }
   } else if (llvm::isa<clang::EnumConstantDecl>(vd)) {
      // Enumeration constant are considered to be 'static' data member in
      // the CINT (and thus ROOT) scheme.
      property |= kIsStatic;
   }
   clang::QualType qt = vd->getType();
   if (llvm::isa<clang::TypedefType>(qt)) {
      property |= kIsTypedef;
   }
   qt = qt.getCanonicalType();
   property = TClingDeclInfo::Property(property, qt);
   const clang::TagType *tt = qt->getAs<clang::TagType>();
   if (tt) {
      // tt->getDecl() might deserialize.
      cling::Interpreter::PushTransactionRAII RAII(fInterp);
      const clang::TagDecl *td = tt->getDecl();
      if (td->isClass()) {
         property |= kIsClass;
      }
      else if (td->isStruct()) {
         property |= kIsStruct;
      }
      else if (td->isUnion()) {
         property |= kIsUnion;
      }
      else if (td->isEnum()) {
         property |= kIsEnum;
      }
   }

   if (const auto *RD = llvm::dyn_cast<RecordDecl>(thisDecl->getDeclContext())) {
      if (RD->isUnion())
         property |= kIsUnionMember;
   }
   // We can't be a namespace, can we?
   //   if (dc->isNamespace() && !dc->isTranslationUnit()) {
   //      property |= kIsNamespace;
   //   }
   return property;
}

long TClingDataMemberInfo::TypeProperty() const
{
   if (!IsValid()) {
      return 0L;
   }
   const clang::ValueDecl *vd = GetTargetValueDecl();
   clang::QualType qt = vd->getType();
   return TClingTypeInfo(fInterp, qt).Property();
}

int TClingDataMemberInfo::TypeSize() const
{
   if (!IsValid()) {
      return -1;
   }

   const clang::ValueDecl *vd = GetTargetValueDecl();
   // Sanity check the current data member.
   clang::Decl::Kind dk = vd->getKind();
   if ((dk != clang::Decl::Field) && (dk != clang::Decl::Var) &&
       (dk != clang::Decl::EnumConstant)) {
      // Error, was not a data member, variable, or enumerator.
      return -1;
   }
   clang::QualType qt = vd->getType();
   if (qt->isIncompleteType()) {
      // We cannot determine the size of forward-declared types.
      return -1;
   }
   clang::ASTContext &context = GetDecl()->getASTContext();
   // Truncate cast to fit to cint interface.
   return static_cast<int>(context.getTypeSizeInChars(qt).getQuantity());
}

const char *TClingDataMemberInfo::TypeName() const
{
   if (!IsValid()) {
      return 0;
   }

   CheckForIoTypeAndName();
   if (!fIoType.empty()) return fIoType.c_str();

   // Note: This must be static because we return a pointer inside it!
   static std::string buf;
   buf.clear();
   const clang::ValueDecl *vd = GetTargetValueDecl();
   clang::QualType vdType = vd->getType();
   // In CINT's version, the type name returns did *not* include any array
   // information, ROOT's existing code depends on it.
   while (vdType->isArrayType()) {
      vdType = GetDecl()->getASTContext().getQualifiedType(vdType->getBaseElementTypeUnsafe(),vdType.getQualifiers());
   }

   // if (we_need_to_do_the_subst_because_the_class_is_a_template_instance_of_double32_t)
   vdType = ROOT::TMetaUtils::ReSubstTemplateArg(vdType, GetClassAsType() );

   ROOT::TMetaUtils::GetFullyQualifiedTypeName(buf, vdType, *fInterp);

   return buf.c_str();
}

const char *TClingDataMemberInfo::TypeTrueName(const ROOT::TMetaUtils::TNormalizedCtxt &normCtxt) const
{
   if (!IsValid()) {
      return 0;
   }

   CheckForIoTypeAndName();
   if (!fIoType.empty()) return fIoType.c_str();

   // Note: This must be static because we return a pointer inside it!
   static std::string buf;
   buf.clear();
   const clang::ValueDecl *vd = GetTargetValueDecl();
   // if (we_need_to_do_the_subst_because_the_class_is_a_template_instance_of_double32_t)
   clang::QualType vdType = ROOT::TMetaUtils::ReSubstTemplateArg(vd->getType(), GetClassAsType());

   ROOT::TMetaUtils::GetNormalizedName(buf, vdType, *fInterp, normCtxt);

   // In CINT's version, the type name returns did *not* include any array
   // information, ROOT's existing code depends on it.
   // This might become part of the implementation of GetNormalizedName.
   while (buf.length() && buf[buf.length()-1] == ']') {
      size_t last = buf.rfind('['); // if this is not the bracket we are looking, the type is malformed.
      if (last != std::string::npos) {
         buf.erase(last);
      }
   }
   return buf.c_str();
}

const char *TClingDataMemberInfo::Name() const
{
   if (!IsValid()) {
      return 0;
   }

   CheckForIoTypeAndName();
   if (!fIoName.empty()) return fIoName.c_str();

   return TClingDeclInfo::Name();
}

const char *TClingDataMemberInfo::Title()
{
   if (!IsValid()) {
      return 0;
   }

   //NOTE: We can't use it as a cache due to the "thoughtful" self iterator
   //if (fTitle.size())
   //   return fTitle.c_str();

   bool titleFound=false;
   // Try to get the comment either from the annotation or the header file if present
   std::string attribute_s;
   const Decl* decl = GetTargetValueDecl();
   for (Decl::attr_iterator attrIt = decl->attr_begin();
        attrIt!=decl->attr_end() && !titleFound ;++attrIt){
      if (0 == ROOT::TMetaUtils::extractAttrString(*attrIt, attribute_s) &&
          attribute_s.find(ROOT::TMetaUtils::propNames::separator) == std::string::npos){
         fTitle = attribute_s;
         titleFound=true;
      }
   }

   if (!titleFound && !decl->isFromASTFile()) {
      // Try to get the comment from the header file if present
      // but not for decls from AST file, where rootcling would have
      // created an annotation
      fTitle = ROOT::TMetaUtils::GetComment(*GetTargetValueDecl()).str();
   }

   return fTitle.c_str();
}

// ValidArrayIndex return a static string (so use it or copy it immediately, do not
// call GrabIndex twice in the same expression) containing the size of the
// array data member.
llvm::StringRef TClingDataMemberInfo::ValidArrayIndex() const
{
   if (!IsValid()) {
      return llvm::StringRef();
   }
   const clang::DeclaratorDecl *FD = llvm::dyn_cast<clang::DeclaratorDecl>(GetTargetValueDecl());
   if (FD) return ROOT::TMetaUtils::DataMemberInfo__ValidArrayIndex(*FD);
   else return llvm::StringRef();
}

