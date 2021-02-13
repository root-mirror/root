// @(#)root/core/meta:$Id$
// Author: Paul Russo   30/07/2012

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TClingBaseClassInfo

Emulation of the CINT BaseClassInfo class.

The CINT C++ interpreter provides an interface to metadata about
the base classes of a class through the BaseClassInfo class.  This
class provides the same functionality, using an interface as close
as possible to BaseClassInfo but the base class metadata comes from
the Clang C++ compiler, not CINT.
*/

#include "TClingBaseClassInfo.h"

#include "TClingClassInfo.h"
#include "TDictionary.h"
#include "TClingUtils.h"

#include "TError.h"

#include "cling/Interpreter/Interpreter.h"
#include "cling/Interpreter/Transaction.h"


#include "clang/AST/ASTContext.h"
#include "clang/AST/Decl.h"
#include "clang/AST/DeclCXX.h"
#include "clang/AST/PrettyPrinter.h"
#include "clang/AST/RecordLayout.h"
#include "clang/AST/Type.h"
#include "clang/AST/CXXInheritance.h"


#include "llvm/Support/Casting.h"
#include "llvm/Support/raw_ostream.h"
#include "llvm/ExecutionEngine/ExecutionEngine.h"
#include "llvm/IR/Module.h"

#include <string>
#include <sstream>

using namespace llvm;
using namespace clang;
using namespace std;

TClingBaseClassInfo::TClingBaseClassInfo(cling::Interpreter* interp,
                                         TClingClassInfo* ci)
   : fInterp(interp), fClassInfo(0), fFirstTime(true), fDescend(false),
     fDecl(0), fIter(0), fBaseInfo(0), fOffset(0L), fClassInfoOwnership(true)
{
   // Constructs a base class iterator on ci; ci == 0 means global scope (which
   // is meaningless). The derived class info passed in as ci is copied.
   if (!ci) {
      fClassInfo = new TClingClassInfo(interp);
      return;
   }
   fClassInfo = new TClingClassInfo(*ci);
   if (!fClassInfo->GetDecl()) {
      return;
   }
   const clang::CXXRecordDecl* CRD =
      llvm::dyn_cast<clang::CXXRecordDecl>(fClassInfo->GetDecl());
   if (!CRD) {
      // We were initialized with something that is not a class.
      // FIXME: We should prevent this from happening!
      return;
   }
   fDecl = CRD;
   {
      // In particular if the base are templated, this might deserialize.
      cling::Interpreter::PushTransactionRAII RAII(fInterp);
      fIter = CRD->bases_begin();
   }
}

TClingBaseClassInfo::TClingBaseClassInfo(cling::Interpreter* interp,
                                         TClingClassInfo* derived,
                                         TClingClassInfo* base)
   : fInterp(interp), fClassInfo(0), fFirstTime(true), fDescend(false),
     fDecl(0), fIter(0), fBaseInfo(0), fOffset(0L), fClassInfoOwnership(false)
{
   // Constructs a single base class base (no iterator) of derived; derived must be != 0.
   // The derived class info is referenced during the lifetime of the TClingBaseClassInfo.
   if (!derived->GetDecl()) {
      return;
   }
   const clang::CXXRecordDecl* CRD =
      llvm::dyn_cast<clang::CXXRecordDecl>(derived->GetDecl());
   const clang::CXXRecordDecl* BaseCRD =
      llvm::dyn_cast<clang::CXXRecordDecl>(base->GetDecl());
   if (!CRD || !BaseCRD) {
      // We were initialized with something that is not a class.
      // FIXME: We should prevent this from happening!
      return;
   }

   fClassInfo = derived;
   fDecl = CRD;
   //CRD->isDerivedFrom(BaseCRD, Paths);
   // Check that base derives from derived.
   clang::CXXBasePaths Paths;

   // CXXRecordDecl::isDerivedFrom can trigger deserialization.
   cling::Interpreter::PushTransactionRAII RAII(fInterp);

   if (!CRD->isDerivedFrom(BaseCRD, Paths)) {
      //Not valid fBaseInfo = 0.
      return;
   }

   fBaseInfo = new TClingClassInfo(*base);
   fIter = CRD->bases_end();
}

TClingBaseClassInfo::TClingBaseClassInfo(const TClingBaseClassInfo& rhs)
   : fInterp(rhs.fInterp), fClassInfo(0), fFirstTime(rhs.fFirstTime),
     fDescend(rhs.fDescend), fDecl(rhs.fDecl), fIter(rhs.fIter), fBaseInfo(0),
     fIterStack(rhs.fIterStack), fOffset(rhs.fOffset), fClassInfoOwnership(true)
{
   // Copies a base class info including the base and derived class infos.
   fClassInfo = new TClingClassInfo(*rhs.fClassInfo);
   fBaseInfo = new TClingClassInfo(*rhs.fBaseInfo);
}

TClingBaseClassInfo& TClingBaseClassInfo::operator=(
   const TClingBaseClassInfo& rhs)
{
   if (this != &rhs) {
      fInterp = rhs.fInterp;
      if (fClassInfoOwnership)
         delete fClassInfo;
      fClassInfo = new TClingClassInfo(*rhs.fClassInfo);
      fFirstTime = rhs.fFirstTime;
      fDescend = rhs.fDescend;
      fDecl = rhs.fDecl;
      fIter = rhs.fIter;
      delete fBaseInfo;
      fBaseInfo = new TClingClassInfo(*rhs.fBaseInfo);
      fIterStack = rhs.fIterStack;
      fOffset = rhs.fOffset;
      fClassInfoOwnership = true;
   }
   return *this;
}

TClingClassInfo *TClingBaseClassInfo::GetBase() const
{
   if (!IsValid()) {
      return 0;
   }
   return fBaseInfo;
}

OffsetPtrFunc_t
TClingBaseClassInfo::GenerateBaseOffsetFunction(TClingClassInfo * fromDerivedClass,
                                                TClingClassInfo* toBaseClass,
                                                void* address, bool isDerivedObject) const
{
   // Generate a function at run-time that would calculate the offset
   // from the parameter derived class to the parameter toBase class for the
   // address.

   // rootcling can trigger this, too, and without CodeGen we cannot use any
   // offset calculation function.
   if (fInterp->isInSyntaxOnlyMode())
      return 0;

   // Get the dedcls for the two classes.
   const clang::RecordDecl* fromDerivedDecl
      = dyn_cast<clang::RecordDecl>(fromDerivedClass->GetDecl());
   if (!fromDerivedDecl) {
      ::Error("TClingBaseClassInfo::GenerateBaseOffsetFunction",
            "Offset of non-class %s is ill-defined!", fromDerivedClass->Name());
      return 0;
   }
   const clang::RecordDecl* toBaseDecl
      = dyn_cast<clang::RecordDecl>(toBaseClass->GetDecl());
   if (!toBaseDecl) {
      ::Error("TClingBaseClassInfo::GenerateBaseOffsetFunction",
            "Offset of non-class %s is ill-defined!", toBaseClass->Name());
      return 0;
   }

   // Make the wrapper name.
   string wrapper_name;
   {
      ostringstream buf;
      buf << "h" << fromDerivedDecl;
      buf << '_';
      buf << "h" << toBaseDecl;
      wrapper_name = buf.str();
   }
   string code;
   // Check whether the function was already generated.
   if (!fInterp->getAddressOfGlobal(wrapper_name)) {
      // Get the class or namespace name.
      string fromDerivedClassName;
      clang::QualType QTDerived(fromDerivedClass->GetType(), 0);
      ROOT::TMetaUtils::GetFullyQualifiedTypeName(fromDerivedClassName,
                                                     QTDerived, *fInterp);
      string toBase_class_name;
      clang::QualType QTtoBase(toBaseClass->GetType(), 0);
      ROOT::TMetaUtils::GetFullyQualifiedTypeName(toBase_class_name,
                                                  QTtoBase, *fInterp);
      //  Write the wrapper code.
      llvm::raw_string_ostream buf(code);
      buf << "extern \"C\" long " + wrapper_name + "(void* address, bool isDerivedObject) {\n"
      // If the object is not derived, will downcast to toBase first.
          << "  " << fromDerivedClassName << " *fromDerived;"
          << "  if (isDerivedObject) {"
          << "    fromDerived = (" << fromDerivedClassName << "*)address;\n"
          << "  } else {\n"
          << "    fromDerived = dynamic_cast<" << fromDerivedClassName << "*>((" << toBase_class_name << "*)address);\n"
          << "  }\n"
          << "  if (!fromDerived) {\n"
          << "    return -1; \n"
          << "  }\n"
          << "  " << toBase_class_name << " *toBase = fromDerived;\n"
          << "  return ((long)toBase - (long)fromDerived);\n}\n";
   }

   // If we have a GV then compileFunction will use it; empty code is enough.
   void* f = fInterp->compileFunction(wrapper_name, code, true /*ifUnique*/,
                                      false /*withAccessControl*/);
   if (!f) {
      ::Error("TClingBaseClassInfo::GenerateBaseOffsetFunction",
            "Compilation failed!");
      return 0;
   }

   return (OffsetPtrFunc_t) f;
}

bool TClingBaseClassInfo::IsValid() const
{
   return
      // inited with a valid class, and
      fClassInfo->IsValid() &&
      // the base class we are iterating over is valid, and
      fDecl &&
      // our current base has a TClingClassInfo, and
      fBaseInfo &&
      // our current base is a valid class
      fBaseInfo->IsValid();
}

int TClingBaseClassInfo::InternalNext(int onlyDirect)
{
   // Exit early if the iterator is already invalid.
   if (!fDecl || !fIter ||
         (fIter == llvm::dyn_cast<clang::CXXRecordDecl>(fDecl)->bases_end())) {
      return 0;
   }

   // Advance to the next valid base.
   while (1) {
      // Advance the iterator.
      if (fFirstTime) {
         // The cint semantics are strange.
         fFirstTime = false;
      }
      else if (!onlyDirect && fDescend) {
         // We previously processed a base class which itself has bases,
         // now we process the bases of that base class.

         // At least getASTRecordLayout() might deserialize.
         cling::Interpreter::PushTransactionRAII RAII(fInterp);
         fDescend = false;
         const clang::RecordType *Ty = fIter->getType()->
                                       getAs<clang::RecordType>();
         // Note: We made sure this would work when we selected the
         //       base for processing.
         clang::CXXRecordDecl *Base = llvm::cast<clang::CXXRecordDecl>(
                                         Ty->getDecl()->getDefinition());
         clang::ASTContext &Context = Base->getASTContext();
         const clang::RecordDecl *RD = llvm::dyn_cast<clang::RecordDecl>(fDecl);
         const clang::ASTRecordLayout &Layout = Context.getASTRecordLayout(RD);
         int64_t offset = Layout.getBaseClassOffset(Base).getQuantity();
         fOffset += static_cast<long>(offset);
         fIterStack.push_back(std::make_pair(std::make_pair(fDecl, fIter),
                                             static_cast<long>(offset)));
         fDecl = Base;
         fIter = Base->bases_begin();
      }
      else {
         // Simple case, move on to the next base class specifier.
         ++fIter;
      }
      // Fix it if we went past the end.
      while (
         (fIter == llvm::dyn_cast<clang::CXXRecordDecl>(fDecl)->bases_end()) &&
         fIterStack.size()
      ) {
         // All done with this base class.
         fDecl = fIterStack.back().first.first;
         fIter = fIterStack.back().first.second;
         fOffset -= fIterStack.back().second;
         fIterStack.pop_back();
         ++fIter;
      }
      // Check for final termination.
      if (fIter == llvm::dyn_cast<clang::CXXRecordDecl>(fDecl)->bases_end()) {
         // We have reached the end of the direct bases, all done.
         delete fBaseInfo;
         fBaseInfo = 0;
         // Iterator is now invalid.
         return 0;
      }
      // Check if current base class is a dependent type, that is, an
      // uninstantiated template class.
      const clang::TagType *Ty = fIter->getType()->getAs<clang::TagType>();
      if (!Ty) {
         // A dependent type (uninstantiated template), skip it.
         continue;
      }
      // Check if current base class has a definition.
      const clang::CXXRecordDecl *Base =
         llvm::cast_or_null<clang::CXXRecordDecl>(Ty->getDecl()->
               getDefinition());
      if (!Base) {
         // No definition yet (just forward declared), skip it.
         continue;
      }
      // Now that we are going to return this base, check to see if
      // we need to examine its bases next call.
      if (!onlyDirect && Base->getNumBases()) {
         fDescend = true;
      }
      // Update info for this base class.
      delete fBaseInfo;
      clang::QualType bType = ROOT::TMetaUtils::ReSubstTemplateArg(fIter->getType(),fClassInfo->GetType());
      fBaseInfo = new TClingClassInfo(fInterp, *bType);
      // Iterator is now valid.
      return 1;
   }
}

int TClingBaseClassInfo::Next(int onlyDirect)
{
   return InternalNext(onlyDirect);
}

int TClingBaseClassInfo::Next()
{
   return Next(1);
}

// This function is updating original one on http://clang.llvm.org/doxygen/CGExprCXX_8cpp_source.html#l01647
// To fit the needs.
static clang::CharUnits computeOffsetHint(clang::ASTContext &Context,
                                          const clang::CXXRecordDecl *Src,
                                          const clang::CXXRecordDecl *Dst,
                                          cling::Interpreter* interp)
{
   clang::CXXBasePaths Paths(/*FindAmbiguities=*/true, /*RecordPaths=*/true,
                      /*DetectVirtual=*/false);

   // If Dst is not derived from Src we can skip the whole computation below and
   // return that Src is not a public base of Dst.  Record all inheritance paths.
   if (!Dst->isDerivedFrom(Src, Paths))
     return clang::CharUnits::fromQuantity(-2);

   unsigned NumPublicPaths = 0;
   clang::CharUnits Offset;

   // Now walk all possible inheritance paths.
   for (clang::CXXBasePaths::paths_iterator I = Paths.begin(), E = Paths.end();
        I != E; ++I) {

     ++NumPublicPaths;

     for (clang::CXXBasePath::iterator J = I->begin(), JE = I->end(); J != JE; ++J) {
       // If the path contains a virtual base class we can't give any hint.
       // -1: no hint.
       if (J->Base->isVirtual())
         return clang::CharUnits::fromQuantity(-1);

       if (NumPublicPaths > 1) // Won't use offsets, skip computation.
         continue;

       // Accumulate the base class offsets.
       cling::Interpreter::PushTransactionRAII RAII(interp);
       const clang::ASTRecordLayout &L = Context.getASTRecordLayout(J->Class);
       Offset += L.getBaseClassOffset(J->Base->getType()->getAsCXXRecordDecl());
     }
   }

   // -2: Src is not a public base of Dst.
   if (NumPublicPaths == 0)
     return clang::CharUnits::fromQuantity(-2);

   // -3: Src is a multiple public base type but never a virtual base type.
   if (NumPublicPaths > 1)
     return clang::CharUnits::fromQuantity(-3);

   // Otherwise, the Src type is a unique public nonvirtual base type of Dst.
   // Return the offset of Src from the origin of Dst.
   return Offset;
 }

ptrdiff_t TClingBaseClassInfo::Offset(void * address, bool isDerivedObject) const
{
   // Compute the offset of the derived class to the base class.

   if (!IsValid()) {
      return -1;
   }
   // Check if current base class has a definition.
   const clang::CXXRecordDecl* Base =
      llvm::cast_or_null<clang::CXXRecordDecl>(fBaseInfo->GetDecl());
   if (!Base) {
      // No definition yet (just forward declared), invalid.
      return -1;
   }
   // If the base class has no virtual inheritance.
   if (!(Property() & kIsVirtualBase)) {
      clang::ASTContext& Context = Base->getASTContext();
      const clang::CXXRecordDecl* RD = llvm::dyn_cast<clang::CXXRecordDecl>(fDecl);
      if (!RD) {
         // No RecordDecl for the class.
         return -1;
      }
      long clang_val = computeOffsetHint(Context, Base, RD, fInterp).getQuantity();
      if (clang_val == -2 || clang_val == -3) {
         TString baseName;
         TString derivedName;
         {
            // Need TNormalizedCtxt otherwise...
            // Note: should we really be issuing a message here?  Shouldn't
            // the caller check and issue the message?
            std::string buf;
            PrintingPolicy Policy(fBaseInfo->GetDecl()->getASTContext().
                                  getPrintingPolicy());
            llvm::raw_string_ostream stream(buf);
            ((const clang::NamedDecl*)fBaseInfo->GetDecl())
               ->getNameForDiagnostic(stream, Policy, /*Qualified=*/true);
            stream.flush();
            baseName = buf;

            buf.clear();
            ((const clang::NamedDecl*)fClassInfo->GetDecl())
               ->getNameForDiagnostic(stream, Policy, /*Qualified=*/true);
            stream.flush();
            derivedName = buf;
         }
         if (clang_val == -2) {
            ::Error("TClingBaseClassInfo::Offset",
                  "The class %s does not derive from the base %s.",
                  derivedName.Data(), baseName.Data());
         } else {
            // clang_val == -3
            ::Error("TClingBaseClassInfo::Offset",
                  "There are multiple paths from derived class %s to base class %s.",
                  derivedName.Data(), baseName.Data());
         }
         clang_val = -1;
      }
      fClassInfo->AddBaseOffsetValue(fBaseInfo->GetDecl(), clang_val);
      return clang_val;
   }
   // Verify the address of the instantiated object
   if (!address) {
      ::Error("TClingBaseClassInfo::Offset", "The address of the object for virtual base offset calculation is not valid.");
      return -1;
   }

   // Virtual inheritance case
   OffsetPtrFunc_t executableFunc = GenerateBaseOffsetFunction(fClassInfo, fBaseInfo, address, isDerivedObject);
   if (executableFunc) {
      fClassInfo->AddBaseOffsetFunction(fBaseInfo->GetDecl(), executableFunc);
      return (*executableFunc)(address, isDerivedObject);
   }

   return -1;
}


long TClingBaseClassInfo::Property() const
{
   if (!IsValid()) {
      return 0L;
   }
   long property = 0L;

   if (fDecl == fClassInfo->GetDecl()) {
      property |= kIsDirectInherit;
   }

   const clang::CXXRecordDecl* CRD
      = llvm::dyn_cast<CXXRecordDecl>(fDecl);
   const clang::CXXRecordDecl* BaseCRD
      = llvm::dyn_cast<CXXRecordDecl>(fBaseInfo->GetDecl());
   if (!CRD || !BaseCRD) {
      ::Error("TClingBaseClassInfo::Property",
            "The derived class or the base class do not have a CXXRecordDecl.");
      return property;
   }

   clang::CXXBasePaths Paths(/*FindAmbiguities=*/false, /*RecordPaths=*/true,
                             /*DetectVirtual=*/true);
   if (!CRD->isDerivedFrom(BaseCRD, Paths)) {
      // Error really unexpected here, because construction / iteration guarantees
      //inheritance;
      ::Error("TClingBaseClassInfo", "Class not derived from given base.");
   }
   if (Paths.getDetectedVirtual()) {
      property |= kIsVirtualBase;
   }

   clang::AccessSpecifier AS = clang::AS_public;
   // Derived: public Mid; Mid : protected Base: Derived inherits protected Base?
   for (clang::CXXBasePaths::const_paths_iterator IB = Paths.begin(), EB = Paths.end();
        AS != clang::AS_private && IB != EB; ++IB) {
      switch (IB->Access) {
         // keep AS unchanged?
         case clang::AS_public: break;
         case clang::AS_protected: AS = clang::AS_protected; break;
         case clang::AS_private: AS = clang::AS_private; break;
         case clang::AS_none: break;
      }
   }
   switch (AS) {
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
         // IMPOSSIBLE
         break;
   }
   return property;
}

long TClingBaseClassInfo::Tagnum() const
{
   if (!IsValid()) {
      return -1L;
   }
   return fBaseInfo->Tagnum();
}

void TClingBaseClassInfo::FullName(std::string &output, const ROOT::TMetaUtils::TNormalizedCtxt &normCtxt) const
{
   if (!IsValid()) {
      output.clear();
      return;
   }
   fBaseInfo->FullName(output,normCtxt);
}

const char* TClingBaseClassInfo::Name() const
{
   if (!IsValid()) {
      return 0;
   }
   return fBaseInfo->Name();
}

const char* TClingBaseClassInfo::TmpltName() const
{
   if (!IsValid()) {
      return 0;
   }
   return fBaseInfo->TmpltName();
}

