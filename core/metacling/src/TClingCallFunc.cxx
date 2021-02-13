// root/core/meta
// vim: sw=3
// Author: Paul Russo   30/07/2012
// Author: Vassil Vassilev   9/02/2013

/*************************************************************************
 * Copyright (C) 1995-2013, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TClingCallFunc
Emulation of the CINT CallFunc class.

The CINT C++ interpreter provides an interface for calling
functions through the generated wrappers in dictionaries with
the CallFunc class. This class provides the same functionality,
using an interface as close as possible to CallFunc but the
function metadata and calling service comes from the Cling
C++ interpreter and the Clang C++ compiler, not CINT.
*/

#include "TClingCallFunc.h"

#include "TClingClassInfo.h"
#include "TClingMethodInfo.h"
#include "TInterpreterValue.h"
#include "TClingUtils.h"

#include "TError.h"
#include "TCling.h"

#include "cling/Interpreter/CompilationOptions.h"
#include "cling/Interpreter/Interpreter.h"
#include "cling/Interpreter/LookupHelper.h"
#include "cling/Interpreter/Transaction.h"
#include "cling/Interpreter/Value.h"
#include "cling/Utils/AST.h"

#include "clang/AST/ASTContext.h"
#include "clang/AST/Decl.h"
#include "clang/AST/DeclCXX.h"
#include "clang/AST/GlobalDecl.h"
#include "clang/AST/PrettyPrinter.h"
#include "clang/AST/RecordLayout.h"
#include "clang/AST/Type.h"
#include "clang/Frontend/CompilerInstance.h"
#include "clang/Lex/Preprocessor.h"
#include "clang/Sema/Sema.h"
#include "clang/Sema/Lookup.h"

#include "llvm/ADT/APInt.h"
#include "llvm/ExecutionEngine/ExecutionEngine.h"
#include "llvm/ExecutionEngine/GenericValue.h"
#include "llvm/Support/Casting.h"
#include "llvm/Support/raw_ostream.h"
#include "llvm/IR/LLVMContext.h"
#include "llvm/IR/DerivedTypes.h"
#include "llvm/IR/Function.h"
#include "llvm/IR/GlobalValue.h"
#include "llvm/IR/Module.h"
#include "llvm/IR/Type.h"

#include "clang/Sema/SemaInternal.h"

#include <map>
#include <string>
#include <sstream>

using namespace ROOT;
using namespace llvm;
using namespace clang;
using namespace std;

static unsigned long long gWrapperSerial = 0LL;
static const string kIndentString("   ");

static map<const Decl *, void *> gWrapperStore;
static map<const Decl *, void *> gCtorWrapperStore;
static map<const Decl *, void *> gDtorWrapperStore;

static
inline
void
indent(ostringstream &buf, int indent_level)
{
   for (int i = 0; i < indent_level; ++i) {
      buf << kIndentString;
   }
}

static
void
EvaluateExpr(cling::Interpreter &interp, const Expr *E, cling::Value &V)
{
   // Evaluate an Expr* and return its cling::Value
   ASTContext &C = interp.getCI()->getASTContext();
   APSInt res;
   if (E->EvaluateAsInt(res, C, /*AllowSideEffects*/Expr::SE_NoSideEffects)) {
      // IntTy or maybe better E->getType()?
      V = cling::Value(C.IntTy, interp);
      // We must use the correct signedness otherwise the zero extension
      // fails if the actual type is strictly less than long long.
      if (res.isSigned())
        V.getLL() = res.getSExtValue();
      else
        V.getULL() = res.getZExtValue();
      return;
   }
   // TODO: Build a wrapper around the expression to avoid decompilation and
   // compilation and other string operations.
   PrintingPolicy Policy(C.getPrintingPolicy());
   Policy.SuppressTagKeyword = true;
   Policy.SuppressUnwrittenScope = false;
   Policy.SuppressInitializers = false;
   Policy.AnonymousTagLocations = false;
   string buf;
   raw_string_ostream out(buf);
   E->printPretty(out, /*Helper=*/0, Policy, /*Indentation=*/0);
   out << ';'; // no value printing
   out.flush();
   // Evaluate() will set V to invalid if evaluation fails.
   interp.evaluate(buf, V);
}

namespace {
   template <typename returnType>
   returnType sv_to(const cling::Value &val)
   {
      QualType QT = val.getType().getCanonicalType();
      if (const BuiltinType *BT =
               dyn_cast<BuiltinType>(&*QT)) {
         //
         //  WARNING!!!
         //
         //  This switch is organized in order-of-declaration
         //  so that the produced assembly code is optimal.
         //  Do not reorder!
         //
         switch (BT->getKind()) {
            case BuiltinType::Void:
               // CINT used to expect a result of 0.
               return (returnType) 0;
               break;
               //
               //  Unsigned Types
               //
            case BuiltinType::Bool:
            case BuiltinType::Char_U: // char on targets where it is unsigned
            case BuiltinType::UChar:
               return (returnType) val.getULL();
               break;

            case BuiltinType::WChar_U:
               // wchar_t on targets where it is unsigned
               // The standard doesn't allow to specify signednedd of wchar_t
               // thus this maps simply to wchar_t.
               return (returnType)(wchar_t) val.getULL();
               break;

            case BuiltinType::Char16:
            case BuiltinType::Char32:
            case BuiltinType::UShort:
            case BuiltinType::UInt:
            case BuiltinType::ULong:
            case BuiltinType::ULongLong:
               return (returnType) val.getULL();
               break;

            case BuiltinType::UInt128:
               // __uint128_t
               break;

               //
               //  Signed Types
               //
            case BuiltinType::Char_S: // char on targets where it is signed
            case BuiltinType::SChar:
               return (returnType) val.getLL();
               break;

            case BuiltinType::WChar_S:
               // wchar_t on targets where it is signed
               // The standard doesn't allow to specify signednedd of wchar_t
               // thus this maps simply to wchar_t.
               return (returnType)(wchar_t) val.getLL();
               break;

            case BuiltinType::Short:
            case BuiltinType::Int:
            case BuiltinType::Long:
            case BuiltinType::LongLong:
               return (returnType) val.getLL();
               break;

            case BuiltinType::Int128:
               break;

            case BuiltinType::Half:
               // half in OpenCL, __fp16 in ARM NEON
               break;

            case BuiltinType::Float:
               return (returnType) val.getFloat();
               break;
            case BuiltinType::Double:
               return (returnType) val.getDouble();
               break;
            case BuiltinType::LongDouble:
               return (returnType) val.getLongDouble();
               break;

            case BuiltinType::NullPtr:
               return (returnType) 0;
               break;

            default:
               break;
         }
      }
      if (QT->isPointerType() || QT->isArrayType() || QT->isRecordType() ||
            QT->isReferenceType()) {
         return (returnType)(long) val.getPtr();
      }
      if (const EnumType *ET = dyn_cast<EnumType>(&*QT)) {
         if (ET->getDecl()->getIntegerType()->hasSignedIntegerRepresentation())
            return (returnType) val.getLL();
         else
            return (returnType) val.getULL();
      }
      if (QT->isMemberPointerType()) {
         const MemberPointerType *MPT = QT->getAs<MemberPointerType>();
         if (MPT->isMemberDataPointer()) {
            return (returnType)(ptrdiff_t)val.getPtr();
         }
         return (returnType)(long) val.getPtr();
      }
      ::Error("TClingCallFunc::sv_to", "Invalid Type!");
      QT->dump();
      return 0;
   }

   static
   long long sv_to_long_long(const cling::Value &val)
   {
      return sv_to<long long>(val);
   }
   static
   unsigned long long sv_to_ulong_long(const cling::Value &val)
   {
      return sv_to<unsigned long long>(val);
   }

} // unnamed namespace.

size_t TClingCallFunc::CalculateMinRequiredArguments()
{
   // This function is non-const to use caching overload of GetDecl()!
   return GetDecl()->getMinRequiredArguments();
}

void *TClingCallFunc::compile_wrapper(const string &wrapper_name, const string &wrapper,
                                      bool withAccessControl/*=true*/)
{
   return fInterp->compileFunction(wrapper_name, wrapper, false /*ifUnique*/,
                                   withAccessControl);
}

void TClingCallFunc::collect_type_info(QualType &QT, ostringstream &typedefbuf, std::ostringstream &callbuf,
                                       string &type_name, EReferenceType &refType, bool &isPointer, int indent_level,
                                       bool forArgument)
{
   //
   //  Collect information about type type of a function parameter
   //  needed for building the wrapper function.
   //
   const FunctionDecl *FD = GetDecl();
   PrintingPolicy Policy(FD->getASTContext().getPrintingPolicy());
   refType = kNotReference;
   if (QT->isRecordType() && forArgument) {
      ROOT::TMetaUtils::GetNormalizedName(type_name, QT, *fInterp, fNormCtxt);
      return;
   }
   if (QT->isFunctionPointerType()) {
      string fp_typedef_name;
      {
         ostringstream nm;
         nm << "FP" << gWrapperSerial++;
         type_name = nm.str();
         raw_string_ostream OS(fp_typedef_name);
         QT.print(OS, Policy, type_name);
         OS.flush();
      }
      for (int i = 0; i < indent_level; ++i) {
         typedefbuf << kIndentString;
      }
      typedefbuf << "typedef " << fp_typedef_name << ";\n";
      return;
   } else if (QT->isMemberPointerType()) {
      string mp_typedef_name;
      {
         ostringstream nm;
         nm << "MP" << gWrapperSerial++;
         type_name = nm.str();
         raw_string_ostream OS(mp_typedef_name);
         QT.print(OS, Policy, type_name);
         OS.flush();
      }
      for (int i = 0; i < indent_level; ++i) {
         typedefbuf << kIndentString;
      }
      typedefbuf << "typedef " << mp_typedef_name << ";\n";
      return;
   } else if (QT->isPointerType()) {
      isPointer = true;
      QT = cast<clang::PointerType>(QT)->getPointeeType();
   } else if (QT->isReferenceType()) {
      if (QT->isRValueReferenceType()) refType = kRValueReference;
      else refType = kLValueReference;
      QT = cast<ReferenceType>(QT)->getPointeeType();
   }
   // Fall through for the array type to deal with reference/pointer ro array type.
   if (QT->isArrayType()) {
      string ar_typedef_name;
      {
         ostringstream ar;
         ar << "AR" << gWrapperSerial++;
         type_name = ar.str();
         raw_string_ostream OS(ar_typedef_name);
         QT.print(OS, Policy, type_name);
         OS.flush();
      }
      for (int i = 0; i < indent_level; ++i) {
         typedefbuf << kIndentString;
      }
      typedefbuf << "typedef " << ar_typedef_name << ";\n";
      return;
   }
   ROOT::TMetaUtils::GetNormalizedName(type_name, QT, *fInterp, fNormCtxt);
}

void TClingCallFunc::make_narg_ctor(const unsigned N, ostringstream &typedefbuf,
                                    ostringstream &callbuf, const string &class_name,
                                    int indent_level)
{
   // Make a code string that follows this pattern:
   //
   // new ClassName(args...)
   //
   const FunctionDecl *FD = GetDecl();

   callbuf << "new " << class_name << "(";
   for (unsigned i = 0U; i < N; ++i) {
      const ParmVarDecl *PVD = FD->getParamDecl(i);
      QualType Ty = PVD->getType();
      QualType QT = Ty.getCanonicalType();
      string type_name;
      EReferenceType refType = kNotReference;
      bool isPointer = false;
      collect_type_info(QT, typedefbuf, callbuf, type_name,
                        refType, isPointer, indent_level, true);
      if (i) {
         callbuf << ',';
         if (i % 2) {
            callbuf << ' ';
         } else {
            callbuf << "\n";
            for (int j = 0; j <= indent_level; ++j) {
               callbuf << kIndentString;
            }
         }
      }
      if (refType != kNotReference) {
         callbuf << "(" << type_name.c_str() <<
                 (refType == kLValueReference ? "&" : "&&") << ")*(" << type_name.c_str() << "*)args["
                 << i << "]";
      } else if (isPointer) {
         callbuf << "*(" << type_name.c_str() << "**)args["
                 << i << "]";
      } else {
         callbuf << "*(" << type_name.c_str() << "*)args[" << i << "]";
      }
   }
   callbuf << ")";
}

void TClingCallFunc::make_narg_call(const std::string &return_type, const unsigned N, ostringstream &typedefbuf,
                                    ostringstream &callbuf, const string &class_name, int indent_level)
{
   //
   // Make a code string that follows this pattern:
   //
   // ((<class>*)obj)-><method>(*(<arg-i-type>*)args[i], ...)
   //
   const FunctionDecl *FD = GetDecl();

   // Sometimes it's necessary that we cast the function we want to call first
   // to its explicit function type before calling it. This is supposed to prevent
   // that we accidentially ending up in a function that is not the one we're
   // supposed to call here (e.g. because the C++ function lookup decides to take
   // another function that better fits).
   // This method has some problems, e.g. when we call a function with default
   // arguments and we don't provide all arguments, we would fail with this pattern.
   // Same applies with member methods which seem to cause parse failures even when
   // we supply the object parameter.
   // Therefore we only use it in cases where we know it works and set this variable
   // to true when we do.
   bool ShouldCastFunction = !isa<CXXMethodDecl>(FD) && N == FD->getNumParams();
   if (ShouldCastFunction) {
      callbuf << "(";
      callbuf << "(";
      callbuf << return_type << " (&)";
      {
         callbuf << "(";
         for (unsigned i = 0U; i < N; ++i) {
            if (i) {
               callbuf << ',';
               if (i % 2) {
                  callbuf << ' ';
               } else {
                  callbuf << "\n";
                  for (int j = 0; j <= indent_level; ++j) {
                     callbuf << kIndentString;
                  }
               }
            }
            const ParmVarDecl *PVD = FD->getParamDecl(i);
            QualType Ty = PVD->getType();
            QualType QT = Ty.getCanonicalType();
            std::string arg_type;
            ROOT::TMetaUtils::GetNormalizedName(arg_type, QT, *fInterp, fNormCtxt);
            callbuf << arg_type;
         }
         if (FD->isVariadic())
            callbuf << ", ...";
         callbuf << ")";
      }

      callbuf << ")";
   }

   if (const CXXMethodDecl *MD = dyn_cast<CXXMethodDecl>(FD)) {
      // This is a class, struct, or union member.
      if (MD->isConst())
         callbuf << "((const " << class_name << "*)obj)->";
      else
         callbuf << "((" << class_name << "*)obj)->";
   } else if (const NamedDecl *ND =
                 dyn_cast<NamedDecl>(GetDeclContext())) {
      // This is a namespace member.
      (void) ND;
      callbuf << class_name << "::";
   }
   //   callbuf << fMethod->Name() << "(";
   {
      std::string name;
      {
         llvm::raw_string_ostream stream(name);
         FD->getNameForDiagnostic(stream, FD->getASTContext().getPrintingPolicy(), /*Qualified=*/false);
      }
      callbuf << name;
   }
   if (ShouldCastFunction) callbuf << ")";

   callbuf << "(";
   for (unsigned i = 0U; i < N; ++i) {
      const ParmVarDecl *PVD = FD->getParamDecl(i);
      QualType Ty = PVD->getType();
      QualType QT = Ty.getCanonicalType();
      string type_name;
      EReferenceType refType = kNotReference;
      bool isPointer = false;
      collect_type_info(QT, typedefbuf, callbuf, type_name, refType, isPointer, indent_level, true);

      if (i) {
         callbuf << ',';
         if (i % 2) {
            callbuf << ' ';
         } else {
            callbuf << "\n";
            for (int j = 0; j <= indent_level; ++j) {
               callbuf << kIndentString;
            }
         }
      }

      if (refType != kNotReference) {
         callbuf << "(" << type_name.c_str() <<
                 (refType == kLValueReference ? "&" : "&&") << ")*(" << type_name.c_str() << "*)args["
                 << i << "]";
      } else if (isPointer) {
         callbuf << "*(" << type_name.c_str() << "**)args["
                 << i << "]";
      } else {
         // pointer falls back to non-pointer case; the argument preserves
         // the "pointerness" (i.e. doesn't reference the value).
         callbuf << "*(" << type_name.c_str() << "*)args[" << i << "]";
      }
   }
   callbuf << ")";
}

void TClingCallFunc::make_narg_ctor_with_return(const unsigned N, const string &class_name,
      ostringstream &buf, int indent_level)
{
   // Make a code string that follows this pattern:
   //
   // if (ret) {
   //    (*(ClassName**)ret) = new ClassName(args...);
   // }
   // else {
   //    new ClassName(args...);
   // }
   //
   for (int i = 0; i < indent_level; ++i) {
      buf << kIndentString;
   }
   buf << "if (ret) {\n";
   ++indent_level;
   {
      ostringstream typedefbuf;
      ostringstream callbuf;
      //
      //  Write the return value assignment part.
      //
      for (int i = 0; i < indent_level; ++i) {
         callbuf << kIndentString;
      }
      callbuf << "(*(" << class_name << "**)ret) = ";
      //
      //  Write the actual new expression.
      //
      make_narg_ctor(N, typedefbuf, callbuf, class_name, indent_level);
      //
      //  End the new expression statement.
      //
      callbuf << ";\n";
      for (int i = 0; i < indent_level; ++i) {
         callbuf << kIndentString;
      }
      callbuf << "return;\n";
      //
      //  Output the whole new expression and return statement.
      //
      buf << typedefbuf.str() << callbuf.str();
   }
   --indent_level;
   for (int i = 0; i < indent_level; ++i) {
      buf << kIndentString;
   }
   buf << "}\n";
   for (int i = 0; i < indent_level; ++i) {
      buf << kIndentString;
   }
   buf << "else {\n";
   ++indent_level;
   {
      ostringstream typedefbuf;
      ostringstream callbuf;
      for (int i = 0; i < indent_level; ++i) {
         callbuf << kIndentString;
      }
      make_narg_ctor(N, typedefbuf, callbuf, class_name, indent_level);
      callbuf << ";\n";
      for (int i = 0; i < indent_level; ++i) {
         callbuf << kIndentString;
      }
      callbuf << "return;\n";
      buf << typedefbuf.str() << callbuf.str();
   }
   --indent_level;
   for (int i = 0; i < indent_level; ++i) {
      buf << kIndentString;
   }
   buf << "}\n";
}

///////////////////////////////////////////////////////////////////////////////
// Returns the DeclContext corresponding to fMethod's Decl.
// \Note that this might be a FunctionDecl or a UsingShadowDecl; we use the
// DeclContext of the UsingShadowDecl e.g. for constructing a derived class
// object, even if invoking a function made available by a using declaration
// of a constructor of a base class (ROOT-11010).

const clang::DeclContext *TClingCallFunc::GetDeclContext() const {
   return fMethod->GetDecl()->getDeclContext();
}

int TClingCallFunc::get_wrapper_code(std::string &wrapper_name, std::string &wrapper)
{
   const FunctionDecl *FD = GetDecl();
   assert(FD && "generate_wrapper called without a function decl!");
   ASTContext &Context = FD->getASTContext();
   PrintingPolicy Policy(Context.getPrintingPolicy());
   //
   //  Get the class or namespace name.
   //
   string class_name;
   const clang::DeclContext *DC = GetDeclContext();
   if (const TypeDecl *TD = dyn_cast<TypeDecl>(DC)) {
      // This is a class, struct, or union member.
      QualType QT(TD->getTypeForDecl(), 0);
      ROOT::TMetaUtils::GetNormalizedName(class_name, QT, *fInterp, fNormCtxt);
   } else if (const NamedDecl *ND = dyn_cast<NamedDecl>(DC)) {
      // This is a namespace member.
      raw_string_ostream stream(class_name);
      ND->getNameForDiagnostic(stream, Policy, /*Qualified=*/true);
      stream.flush();
   }
   //
   //  Check to make sure that we can
   //  instantiate and codegen this function.
   //
   bool needInstantiation = false;
   const FunctionDecl *Definition = 0;
   if (!FD->isDefined(Definition)) {
      FunctionDecl::TemplatedKind TK = FD->getTemplatedKind();
      switch (TK) {
      case FunctionDecl::TK_NonTemplate: {
         // Ordinary function, not a template specialization.
         // Note: This might be ok, the body might be defined
         //       in a library, and all we have seen is the
         //       header file.
         //::Error("TClingCallFunc::make_wrapper",
         //      "Cannot make wrapper for a function which is "
         //      "declared but not defined!");
         // return 0;
      } break;
      case FunctionDecl::TK_FunctionTemplate: {
         // This decl is actually a function template,
         // not a function at all.
         ::Error("TClingCallFunc::make_wrapper", "Cannot make wrapper for a function template!");
         return 0;
      } break;
      case FunctionDecl::TK_MemberSpecialization: {
         // This function is the result of instantiating an ordinary
         // member function of a class template, or of instantiating
         // an ordinary member function of a class member of a class
         // template, or of specializing a member function template
         // of a class template, or of specializing a member function
         // template of a class member of a class template.
         if (!FD->isTemplateInstantiation()) {
            // We are either TSK_Undeclared or
            // TSK_ExplicitSpecialization.
            // Note: This might be ok, the body might be defined
            //       in a library, and all we have seen is the
            //       header file.
            //::Error("TClingCallFunc::make_wrapper",
            //      "Cannot make wrapper for a function template "
            //      "explicit specialization which is declared "
            //      "but not defined!");
            // return 0;
            break;
         }
         const FunctionDecl *Pattern = FD->getTemplateInstantiationPattern();
         if (!Pattern) {
            ::Error("TClingCallFunc::make_wrapper", "Cannot make wrapper for a member function "
                                                    "instantiation with no pattern!");
            return 0;
         }
         FunctionDecl::TemplatedKind PTK = Pattern->getTemplatedKind();
         TemplateSpecializationKind PTSK = Pattern->getTemplateSpecializationKind();
         if (
            // The pattern is an ordinary member function.
            (PTK == FunctionDecl::TK_NonTemplate) ||
            // The pattern is an explicit specialization, and
            // so is not a template.
            ((PTK != FunctionDecl::TK_FunctionTemplate) &&
             ((PTSK == TSK_Undeclared) || (PTSK == TSK_ExplicitSpecialization)))) {
            // Note: This might be ok, the body might be defined
            //       in a library, and all we have seen is the
            //       header file.
            break;
         } else if (!Pattern->hasBody()) {
            ::Error("TClingCallFunc::make_wrapper", "Cannot make wrapper for a member function "
                                                    "instantiation with no body!");
            return 0;
         }
         if (FD->isImplicitlyInstantiable()) {
            needInstantiation = true;
         }
      } break;
      case FunctionDecl::TK_FunctionTemplateSpecialization: {
         // This function is the result of instantiating a function
         // template or possibly an explicit specialization of a
         // function template.  Could be a namespace scope function or a
         // member function.
         if (!FD->isTemplateInstantiation()) {
            // We are either TSK_Undeclared or
            // TSK_ExplicitSpecialization.
            // Note: This might be ok, the body might be defined
            //       in a library, and all we have seen is the
            //       header file.
            //::Error("TClingCallFunc::make_wrapper",
            //      "Cannot make wrapper for a function template "
            //      "explicit specialization which is declared "
            //      "but not defined!");
            // return 0;
            break;
         }
         const FunctionDecl *Pattern = FD->getTemplateInstantiationPattern();
         if (!Pattern) {
            ::Error("TClingCallFunc::make_wrapper", "Cannot make wrapper for a function template"
                                                    "instantiation with no pattern!");
            return 0;
         }
         FunctionDecl::TemplatedKind PTK = Pattern->getTemplatedKind();
         TemplateSpecializationKind PTSK = Pattern->getTemplateSpecializationKind();
         if (
            // The pattern is an ordinary member function.
            (PTK == FunctionDecl::TK_NonTemplate) ||
            // The pattern is an explicit specialization, and
            // so is not a template.
            ((PTK != FunctionDecl::TK_FunctionTemplate) &&
             ((PTSK == TSK_Undeclared) || (PTSK == TSK_ExplicitSpecialization)))) {
            // Note: This might be ok, the body might be defined
            //       in a library, and all we have seen is the
            //       header file.
            break;
         }
         if (!Pattern->hasBody()) {
            ::Error("TClingCallFunc::make_wrapper", "Cannot make wrapper for a function template"
                                                    "instantiation with no body!");
            return 0;
         }
         if (FD->isImplicitlyInstantiable()) {
            needInstantiation = true;
         }
      } break;
      case FunctionDecl::TK_DependentFunctionTemplateSpecialization: {
         // This function is the result of instantiating or
         // specializing a  member function of a class template,
         // or a member function of a class member of a class template,
         // or a member function template of a class template, or a
         // member function template of a class member of a class
         // template where at least some part of the function is
         // dependent on a template argument.
         if (!FD->isTemplateInstantiation()) {
            // We are either TSK_Undeclared or
            // TSK_ExplicitSpecialization.
            // Note: This might be ok, the body might be defined
            //       in a library, and all we have seen is the
            //       header file.
            //::Error("TClingCallFunc::make_wrapper",
            //      "Cannot make wrapper for a dependent function "
            //      "template explicit specialization which is declared "
            //      "but not defined!");
            // return 0;
            break;
         }
         const FunctionDecl *Pattern = FD->getTemplateInstantiationPattern();
         if (!Pattern) {
            ::Error("TClingCallFunc::make_wrapper", "Cannot make wrapper for a dependent function template"
                                                    "instantiation with no pattern!");
            return 0;
         }
         FunctionDecl::TemplatedKind PTK = Pattern->getTemplatedKind();
         TemplateSpecializationKind PTSK = Pattern->getTemplateSpecializationKind();
         if (
            // The pattern is an ordinary member function.
            (PTK == FunctionDecl::TK_NonTemplate) ||
            // The pattern is an explicit specialization, and
            // so is not a template.
            ((PTK != FunctionDecl::TK_FunctionTemplate) &&
             ((PTSK == TSK_Undeclared) || (PTSK == TSK_ExplicitSpecialization)))) {
            // Note: This might be ok, the body might be defined
            //       in a library, and all we have seen is the
            //       header file.
            break;
         }
         if (!Pattern->hasBody()) {
            ::Error("TClingCallFunc::make_wrapper", "Cannot make wrapper for a dependent function template"
                                                    "instantiation with no body!");
            return 0;
         }
         if (FD->isImplicitlyInstantiable()) {
            needInstantiation = true;
         }
      } break;
      default: {
         // Will only happen if clang implementation changes.
         // Protect ourselves in case that happens.
         ::Error("TClingCallFunc::make_wrapper", "Unhandled template kind!");
         return 0;
      } break;
      }
      // We do not set needInstantiation to true in these cases:
      //
      // isInvalidDecl()
      // TSK_Undeclared
      // TSK_ExplicitInstantiationDefinition
      // TSK_ExplicitSpecialization && !getClassScopeSpecializationPattern()
      // TSK_ExplicitInstantiationDeclaration &&
      //    getTemplateInstantiationPattern() &&
      //    PatternDecl->hasBody() &&
      //    !PatternDecl->isInlined()
      //
      // Set it true in these cases:
      //
      // TSK_ImplicitInstantiation
      // TSK_ExplicitInstantiationDeclaration && (!getPatternDecl() ||
      //    !PatternDecl->hasBody() || PatternDecl->isInlined())
      //
   }
   if (needInstantiation) {
      clang::FunctionDecl *FDmod = const_cast<clang::FunctionDecl *>(FD);
      clang::Sema &S = fInterp->getSema();
      // Could trigger deserialization of decls.
      cling::Interpreter::PushTransactionRAII RAII(fInterp);
      S.InstantiateFunctionDefinition(SourceLocation(), FDmod,
                                      /*Recursive=*/true,
                                      /*DefinitionRequired=*/true);
      if (!FD->isDefined(Definition)) {
         ::Error("TClingCallFunc::make_wrapper", "Failed to force template instantiation!");
         return 0;
      }
   }
   if (Definition) {
      FunctionDecl::TemplatedKind TK = Definition->getTemplatedKind();
      switch (TK) {
      case FunctionDecl::TK_NonTemplate: {
         // Ordinary function, not a template specialization.
         if (Definition->isDeleted()) {
            ::Error("TClingCallFunc::make_wrapper", "Cannot make wrapper for a deleted function!");
            return 0;
         } else if (Definition->isLateTemplateParsed()) {
            ::Error("TClingCallFunc::make_wrapper", "Cannot make wrapper for a late template parsed "
                                                    "function!");
            return 0;
         }
         // else if (Definition->isDefaulted()) {
         //   // Might not have a body, but we can still use it.
         //}
         // else {
         //   // Has a body.
         //}
      } break;
      case FunctionDecl::TK_FunctionTemplate: {
         // This decl is actually a function template,
         // not a function at all.
         ::Error("TClingCallFunc::make_wrapper", "Cannot make wrapper for a function template!");
         return 0;
      } break;
      case FunctionDecl::TK_MemberSpecialization: {
         // This function is the result of instantiating an ordinary
         // member function of a class template or of a member class
         // of a class template.
         if (Definition->isDeleted()) {
            ::Error("TClingCallFunc::make_wrapper", "Cannot make wrapper for a deleted member function "
                                                    "of a specialization!");
            return 0;
         } else if (Definition->isLateTemplateParsed()) {
            ::Error("TClingCallFunc::make_wrapper", "Cannot make wrapper for a late template parsed "
                                                    "member function of a specialization!");
            return 0;
         }
         // else if (Definition->isDefaulted()) {
         //   // Might not have a body, but we can still use it.
         //}
         // else {
         //   // Has a body.
         //}
      } break;
      case FunctionDecl::TK_FunctionTemplateSpecialization: {
         // This function is the result of instantiating a function
         // template or possibly an explicit specialization of a
         // function template.  Could be a namespace scope function or a
         // member function.
         if (Definition->isDeleted()) {
            ::Error("TClingCallFunc::make_wrapper", "Cannot make wrapper for a deleted function "
                                                    "template specialization!");
            return 0;
         } else if (Definition->isLateTemplateParsed()) {
            ::Error("TClingCallFunc::make_wrapper", "Cannot make wrapper for a late template parsed "
                                                    "function template specialization!");
            return 0;
         }
         // else if (Definition->isDefaulted()) {
         //   // Might not have a body, but we can still use it.
         //}
         // else {
         //   // Has a body.
         //}
      } break;
      case FunctionDecl::TK_DependentFunctionTemplateSpecialization: {
         // This function is the result of instantiating or
         // specializing a  member function of a class template,
         // or a member function of a class member of a class template,
         // or a member function template of a class template, or a
         // member function template of a class member of a class
         // template where at least some part of the function is
         // dependent on a template argument.
         if (Definition->isDeleted()) {
            ::Error("TClingCallFunc::make_wrapper", "Cannot make wrapper for a deleted dependent function "
                                                    "template specialization!");
            return 0;
         } else if (Definition->isLateTemplateParsed()) {
            ::Error("TClingCallFunc::make_wrapper", "Cannot make wrapper for a late template parsed "
                                                    "dependent function template specialization!");
            return 0;
         }
         // else if (Definition->isDefaulted()) {
         //   // Might not have a body, but we can still use it.
         //}
         // else {
         //   // Has a body.
         //}
      } break;
      default: {
         // Will only happen if clang implementation changes.
         // Protect ourselves in case that happens.
         ::Error("TClingCallFunc::make_wrapper", "Unhandled template kind!");
         return 0;
      } break;
      }
   }
   unsigned min_args = GetMinRequiredArguments();
   unsigned num_params = FD->getNumParams();
   //
   //  Make the wrapper name.
   //
   {
      ostringstream buf;
      buf << "__cf";
      // const NamedDecl* ND = dyn_cast<NamedDecl>(FD);
      // string mn;
      // fInterp->maybeMangleDeclName(ND, mn);
      // buf << '_' << mn;
      buf << '_' << gWrapperSerial++;
      wrapper_name = buf.str();
   }
   //
   //  Write the wrapper code.
   // FIXME: this should be synthesized into the AST!
   //
   int indent_level = 0;
   ostringstream buf;
   buf << "#pragma clang diagnostic push\n"
          "#pragma clang diagnostic ignored \"-Wformat-security\"\n"
          "__attribute__((used)) "
          "extern \"C\" void ";
   buf << wrapper_name;
   buf << "(void* obj, int nargs, void** args, void* ret)\n"
          "{\n";
   ++indent_level;
   if (min_args == num_params) {
      // No parameters with defaults.
      make_narg_call_with_return(num_params, class_name, buf, indent_level);
   } else {
      // We need one function call clause compiled for every
      // possible number of arguments per call.
      for (unsigned N = min_args; N <= num_params; ++N) {
         for (int i = 0; i < indent_level; ++i) {
            buf << kIndentString;
         }
         buf << "if (nargs == " << N << ") {\n";
         ++indent_level;
         make_narg_call_with_return(N, class_name, buf, indent_level);
         --indent_level;
         for (int i = 0; i < indent_level; ++i) {
            buf << kIndentString;
         }
         buf << "}\n";
      }
   }
   --indent_level;
   buf << "}\n"
          "#pragma clang diagnostic pop";
   wrapper = buf.str();
   return 1;
}

void TClingCallFunc::make_narg_call_with_return(const unsigned N, const string &class_name,
      ostringstream &buf, int indent_level)
{
   // Make a code string that follows this pattern:
   //
   // if (ret) {
   //    new (ret) (return_type) ((class_name*)obj)->func(args...);
   // }
   // else {
   //    ((class_name*)obj)->func(args...);
   // }
   //
   const FunctionDecl *FD = GetDecl();
   if (const CXXConstructorDecl *CD = dyn_cast<CXXConstructorDecl>(FD)) {
      if (N <= 1 && llvm::isa<UsingShadowDecl>(GetFunctionOrShadowDecl())) {
         auto SpecMemKind = fInterp->getSema().getSpecialMember(CD);
         if ((N == 0 && SpecMemKind == clang::Sema::CXXDefaultConstructor) ||
             (N == 1 &&
              (SpecMemKind == clang::Sema::CXXCopyConstructor || SpecMemKind == clang::Sema::CXXMoveConstructor))) {
            // Using declarations cannot inject special members; do not call them
            // as such. This might happen by using `Base(Base&, int = 12)`, which
            // is fine to be called as `Derived d(someBase, 42)` but not as
            // copy constructor of `Derived`.
            return;
         }
      }
      make_narg_ctor_with_return(N, class_name, buf, indent_level);
      return;
   }
   QualType QT = FD->getReturnType().getCanonicalType();
   if (QT->isVoidType()) {
      ostringstream typedefbuf;
      ostringstream callbuf;
      for (int i = 0; i < indent_level; ++i) {
         callbuf << kIndentString;
      }
      make_narg_call("void", N, typedefbuf, callbuf, class_name, indent_level);
      callbuf << ";\n";
      for (int i = 0; i < indent_level; ++i) {
         callbuf << kIndentString;
      }
      callbuf << "return;\n";
      buf << typedefbuf.str() << callbuf.str();
   } else {
      for (int i = 0; i < indent_level; ++i) {
         buf << kIndentString;
      }

      string type_name;
      EReferenceType refType = kNotReference;
      bool isPointer = false;

      buf << "if (ret) {\n";
      ++indent_level;
      {
         ostringstream typedefbuf;
         ostringstream callbuf;
         //
         //  Write the placement part of the placement new.
         //
         for (int i = 0; i < indent_level; ++i) {
            callbuf << kIndentString;
         }
         callbuf << "new (ret) ";
         collect_type_info(QT, typedefbuf, callbuf, type_name,
                           refType, isPointer, indent_level, false);
         //
         //  Write the type part of the placement new.
         //
         callbuf << "(" << type_name.c_str();
         if (refType != kNotReference) {
            callbuf << "*) (&";
            type_name += "&";
         } else if (isPointer) {
            callbuf << "*) (";
            type_name += "*";
         } else {
            callbuf << ") (";
         }
         //
         //  Write the actual function call.
         //
         make_narg_call(type_name, N, typedefbuf, callbuf, class_name, indent_level);
         //
         //  End the placement new.
         //
         callbuf << ");\n";
         for (int i = 0; i < indent_level; ++i) {
            callbuf << kIndentString;
         }
         callbuf << "return;\n";
         //
         //  Output the whole placement new expression and return statement.
         //
         buf << typedefbuf.str() << callbuf.str();
      }
      --indent_level;
      for (int i = 0; i < indent_level; ++i) {
         buf << kIndentString;
      }
      buf << "}\n";
      for (int i = 0; i < indent_level; ++i) {
         buf << kIndentString;
      }
      buf << "else {\n";
      ++indent_level;
      {
         ostringstream typedefbuf;
         ostringstream callbuf;
         for (int i = 0; i < indent_level; ++i) {
            callbuf << kIndentString;
         }
         make_narg_call(type_name, N, typedefbuf, callbuf, class_name, indent_level);
         callbuf << ";\n";
         for (int i = 0; i < indent_level; ++i) {
            callbuf << kIndentString;
         }
         callbuf << "return;\n";
         buf << typedefbuf.str() << callbuf.str();
      }
      --indent_level;
      for (int i = 0; i < indent_level; ++i) {
         buf << kIndentString;
      }
      buf << "}\n";
   }
}

tcling_callfunc_Wrapper_t TClingCallFunc::make_wrapper()
{
   R__LOCKGUARD_CLING(gInterpreterMutex);

   const Decl *D = GetFunctionOrShadowDecl();
   string wrapper_name;
   string wrapper;

   if (get_wrapper_code(wrapper_name, wrapper) == 0) return 0;

   //fprintf(stderr, "%s\n", wrapper.c_str());
   //
   //  Compile the wrapper code.
   //
   void *F = compile_wrapper(wrapper_name, wrapper);
   if (F) {
      gWrapperStore.insert(make_pair(D, F));
   } else {
      ::Error("TClingCallFunc::make_wrapper",
            "Failed to compile\n  ==== SOURCE BEGIN ====\n%s\n  ==== SOURCE END ====",
            wrapper.c_str());
   }
   return (tcling_callfunc_Wrapper_t)F;
}

tcling_callfunc_ctor_Wrapper_t TClingCallFunc::make_ctor_wrapper(const TClingClassInfo *info,
      ROOT::TMetaUtils::EIOCtorCategory kind, const std::string &type_name)
{
   // Make a code string that follows this pattern:
   //
   // void
   // unique_wrapper_ddd(void** ret, void* arena, unsigned long nary)
   // {
   //    if (!arena) {
   //       if (!nary) {
   //          *ret = new ClassName;
   //       }
   //       else {
   //          *ret = new ClassName[nary];
   //       }
   //    }
   //    else {
   //       if (!nary) {
   //          *ret = new (arena) ClassName;
   //       }
   //       else {
   //          *ret = new (arena) ClassName[nary];
   //       }
   //    }
   // }
   //
   // When I/O constructor used:
   //
   // void
   // unique_wrapper_ddd(void** ret, void* arena, unsigned long nary)
   // {
   //    if (!arena) {
   //       if (!nary) {
   //          *ret = new ClassName((TRootIOCtor*)nullptr);
   //       }
   //       else {
   //          char *buf = malloc(nary * sizeof(ClassName));
   //          for (int k=0;k<nary;++k)
   //             new (buf + k * sizeof(ClassName)) ClassName((TRootIOCtor*)nullptr);
   //          *ret = buf;
   //       }
   //    }
   //    else {
   //       if (!nary) {
   //          *ret = new (arena) ClassName((TRootIOCtor*)nullptr);
   //       }
   //       else {
   //          for (int k=0;k<nary;++k)
   //             new ((char *) arena + k * sizeof(ClassName)) ClassName((TRootIOCtor*)nullptr);
   //          *ret = arena;
   //       }
   //    }
   // }
   //
   //
   // Note:
   //
   // If the class is of POD type, the form:
   //
   //    new ClassName;
   //
   // does not initialize the object at all, and the form:
   //
   //    new ClassName();
   //
   // default-initializes the object.
   //
   // We are using the form without parentheses because that is what
   // CINT did.
   //
   //--
   ASTContext &Context = info->GetDecl()->getASTContext();
   PrintingPolicy Policy(Context.getPrintingPolicy());
   Policy.SuppressTagKeyword = true;
   Policy.SuppressUnwrittenScope = true;
   //
   //  Get the class or namespace name.
   //
   string class_name;
   if (const TypeDecl *TD = dyn_cast<TypeDecl>(info->GetDecl())) {
      // This is a class, struct, or union member.
      QualType QT(TD->getTypeForDecl(), 0);
      ROOT::TMetaUtils::GetNormalizedName(class_name, QT, *fInterp, fNormCtxt);
   } else if (const NamedDecl *ND = dyn_cast<NamedDecl>(info->GetDecl())) {
      // This is a namespace member.
      raw_string_ostream stream(class_name);
      ND->getNameForDiagnostic(stream, Policy, /*Qualified=*/true);
      stream.flush();
   }


   //
   //  Make the wrapper name.
   //
   string wrapper_name;
   {
      ostringstream buf;
      buf << "__ctor";
      //const NamedDecl* ND = dyn_cast<NamedDecl>(FD);
      //string mn;
      //fInterp->maybeMangleDeclName(ND, mn);
      //buf << '_dtor_' << mn;
      buf << '_' << gWrapperSerial++;
      wrapper_name = buf.str();
   }

   string constr_arg;
   if (kind == ROOT::TMetaUtils::EIOCtorCategory::kIOPtrType)
      constr_arg = string("((") + type_name + "*)nullptr)";
   else if (kind == ROOT::TMetaUtils::EIOCtorCategory::kIORefType)
      constr_arg = string("(*((") + type_name + "*)arena))";

   //
   //  Write the wrapper code.
   //
   int indent_level = 0;
   ostringstream buf;
   buf << "__attribute__((used)) ";
   buf << "extern \"C\" void ";
   buf << wrapper_name;
   buf << "(void** ret, void* arena, unsigned long nary)\n";
   buf << "{\n";

   //    if (!arena) {
   //       if (!nary) {
   //          *ret = new ClassName;
   //       }
   //       else {
   //          *ret = new ClassName[nary];
   //       }
   //    }
   indent(buf, ++indent_level);
   buf << "if (!arena) {\n";
   indent(buf, ++indent_level);
   buf << "if (!nary) {\n";
   indent(buf, ++indent_level);
   buf << "*ret = new " << class_name << constr_arg << ";\n";
   indent(buf, --indent_level);
   buf << "}\n";
   indent(buf, indent_level);
   buf << "else {\n";
   indent(buf, ++indent_level);
   if (constr_arg.empty()) {
      buf << "*ret = new " << class_name << "[nary];\n";
   } else {
      buf << "char *buf = (char *) malloc(nary * sizeof(" << class_name << "));\n";
      indent(buf, indent_level);
      buf << "for (int k=0;k<nary;++k)\n";
      indent(buf, ++indent_level);
      buf << "new (buf + k * sizeof(" << class_name << ")) " << class_name << constr_arg << ";\n";
      indent(buf, --indent_level);
      buf << "*ret = buf;\n";
   }
   indent(buf, --indent_level);
   buf << "}\n";
   indent(buf, --indent_level);
   buf << "}\n";
   //    else {
   //       if (!nary) {
   //          *ret = new (arena) ClassName;
   //       }
   //       else {
   //          *ret = new (arena) ClassName[nary];
   //       }
   //    }
   indent(buf, indent_level);
   buf << "else {\n";
   indent(buf, ++indent_level);
   buf << "if (!nary) {\n";
   indent(buf, ++indent_level);
   buf << "*ret = new (arena) " << class_name << constr_arg << ";\n";
   indent(buf, --indent_level);
   buf << "}\n";
   indent(buf, indent_level);
   buf << "else {\n";
   indent(buf, ++indent_level);
   if (constr_arg.empty()) {
      buf << "*ret = new (arena) " << class_name << "[nary];\n";
   } else {
      buf << "for (int k=0;k<nary;++k)\n";
      indent(buf, ++indent_level);
      buf << "new ((char *) arena + k * sizeof(" << class_name << ")) " << class_name << constr_arg << ";\n";
      indent(buf, --indent_level);
      buf << "*ret = arena;\n";
   }
   indent(buf, --indent_level);
   buf << "}\n";
   indent(buf, --indent_level);
   buf << "}\n";
   // End wrapper.
   --indent_level;
   buf << "}\n";
   // Done.
   string wrapper(buf.str());
   //fprintf(stderr, "%s\n", wrapper.c_str());
   //
   //  Compile the wrapper code.
   //
   void *F = compile_wrapper(wrapper_name, wrapper,
                             /*withAccessControl=*/false);
   if (F) {
      gCtorWrapperStore.insert(make_pair(info->GetDecl(), F));
   } else {
      ::Error("TClingCallFunc::make_ctor_wrapper",
            "Failed to compile\n  ==== SOURCE BEGIN ====\n%s\n  ==== SOURCE END ====",
            wrapper.c_str());
   }
   return (tcling_callfunc_ctor_Wrapper_t)F;
}

tcling_callfunc_dtor_Wrapper_t
TClingCallFunc::make_dtor_wrapper(const TClingClassInfo *info)
{
   // Make a code string that follows this pattern:
   //
   // void
   // unique_wrapper_ddd(void* obj, unsigned long nary, int withFree)
   // {
   //    if (withFree) {
   //       if (!nary) {
   //          delete (ClassName*) obj;
   //       }
   //       else {
   //          delete[] (ClassName*) obj;
   //       }
   //    }
   //    else {
   //       typedef ClassName DtorName;
   //       if (!nary) {
   //          ((ClassName*)obj)->~DtorName();
   //       }
   //       else {
   //          for (unsigned long i = nary - 1; i > -1; --i) {
   //             (((ClassName*)obj)+i)->~DtorName();
   //          }
   //       }
   //    }
   // }
   //
   //--
   ASTContext &Context = info->GetDecl()->getASTContext();
   PrintingPolicy Policy(Context.getPrintingPolicy());
   Policy.SuppressTagKeyword = true;
   Policy.SuppressUnwrittenScope = true;
   //
   //  Get the class or namespace name.
   //
   string class_name;
   if (const TypeDecl *TD = dyn_cast<TypeDecl>(info->GetDecl())) {
      // This is a class, struct, or union member.
      QualType QT(TD->getTypeForDecl(), 0);
      ROOT::TMetaUtils::GetNormalizedName(class_name, QT, *fInterp, fNormCtxt);
   } else if (const NamedDecl *ND = dyn_cast<NamedDecl>(info->GetDecl())) {
      // This is a namespace member.
      raw_string_ostream stream(class_name);
      ND->getNameForDiagnostic(stream, Policy, /*Qualified=*/true);
      stream.flush();
   }
   //
   //  Make the wrapper name.
   //
   string wrapper_name;
   {
      ostringstream buf;
      buf << "__dtor";
      //const NamedDecl* ND = dyn_cast<NamedDecl>(FD);
      //string mn;
      //fInterp->maybeMangleDeclName(ND, mn);
      //buf << '_dtor_' << mn;
      buf << '_' << gWrapperSerial++;
      wrapper_name = buf.str();
   }
   //
   //  Write the wrapper code.
   //
   int indent_level = 0;
   ostringstream buf;
   buf << "__attribute__((used)) ";
   buf << "extern \"C\" void ";
   buf << wrapper_name;
   buf << "(void* obj, unsigned long nary, int withFree)\n";
   buf << "{\n";
   //    if (withFree) {
   //       if (!nary) {
   //          delete (ClassName*) obj;
   //       }
   //       else {
   //          delete[] (ClassName*) obj;
   //       }
   //    }
   ++indent_level;
   indent(buf, indent_level);
   buf << "if (withFree) {\n";
   ++indent_level;
   indent(buf, indent_level);
   buf << "if (!nary) {\n";
   ++indent_level;
   indent(buf, indent_level);
   buf << "delete (" << class_name << "*) obj;\n";
   --indent_level;
   indent(buf, indent_level);
   buf << "}\n";
   indent(buf, indent_level);
   buf << "else {\n";
   ++indent_level;
   indent(buf, indent_level);
   buf << "delete[] (" << class_name << "*) obj;\n";
   --indent_level;
   indent(buf, indent_level);
   buf << "}\n";
   --indent_level;
   indent(buf, indent_level);
   buf << "}\n";
   //    else {
   //       typedef ClassName Nm;
   //       if (!nary) {
   //          ((Nm*)obj)->~Nm();
   //       }
   //       else {
   //          for (unsigned long i = nary - 1; i > -1; --i) {
   //             (((Nm*)obj)+i)->~Nm();
   //          }
   //       }
   //    }
   indent(buf, indent_level);
   buf << "else {\n";
   ++indent_level;
   indent(buf, indent_level);
   buf << "typedef " << class_name << " Nm;\n";
   buf << "if (!nary) {\n";
   ++indent_level;
   indent(buf, indent_level);
   buf << "((Nm*)obj)->~Nm();\n";
   --indent_level;
   indent(buf, indent_level);
   buf << "}\n";
   indent(buf, indent_level);
   buf << "else {\n";
   ++indent_level;
   indent(buf, indent_level);
   buf << "do {\n";
   ++indent_level;
   indent(buf, indent_level);
   buf << "(((Nm*)obj)+(--nary))->~Nm();\n";
   --indent_level;
   indent(buf, indent_level);
   buf << "} while (nary);\n";
   --indent_level;
   indent(buf, indent_level);
   buf << "}\n";
   --indent_level;
   indent(buf, indent_level);
   buf << "}\n";
   // End wrapper.
   --indent_level;
   buf << "}\n";
   // Done.
   string wrapper(buf.str());
   //fprintf(stderr, "%s\n", wrapper.c_str());
   //
   //  Compile the wrapper code.
   //
   void *F = compile_wrapper(wrapper_name, wrapper,
                             /*withAccessControl=*/false);
   if (F) {
      gDtorWrapperStore.insert(make_pair(info->GetDecl(), F));
   } else {
      ::Error("TClingCallFunc::make_dtor_wrapper",
            "Failed to compile\n  ==== SOURCE BEGIN ====\n%s\n  ==== SOURCE END ====",
            wrapper.c_str());
   }

   return (tcling_callfunc_dtor_Wrapper_t)F;
}

class ValHolder {
public:
   union {
      long double ldbl;
      double dbl;
      float flt;
      //__uint128_t ui128;
      //__int128_t i128;
      unsigned long long ull;
      long long ll;
      unsigned long ul;
      long l;
      unsigned int ui;
      int i;
      unsigned short us;
      short s;
      //char32_t c32;
      //char16_t c16;
      //unsigned wchar_t uwc; - non-standard
      wchar_t wc;
      unsigned char uc;
      signed char sc;
      char c;
      bool b;
      void *vp;
   } u;
};

void TClingCallFunc::exec(void *address, void *ret)
{
   SmallVector<ValHolder, 8> vh_ary;
   SmallVector<void *, 8> vp_ary;

   unsigned num_args = fArgVals.size();
   {
      R__LOCKGUARD_CLING(gInterpreterMutex);

      const FunctionDecl *FD = GetDecl();

      //
      //  Convert the arguments from cling::Value to their
      //  actual type and store them in a holder for passing to the
      //  wrapper function by pointer to value.
      //
      unsigned num_params = FD->getNumParams();

      if (num_args < GetMinRequiredArguments()) {
         ::Error("TClingCallFunc::exec",
               "Not enough arguments provided for %s (%d instead of the minimum %d)",
               fMethod->Name(),
               num_args, (int)GetMinRequiredArguments());
         return;
      }
      if (address == 0 && dyn_cast<CXXMethodDecl>(FD)
          && !(dyn_cast<CXXMethodDecl>(FD))->isStatic()
          && !dyn_cast<CXXConstructorDecl>(FD)) {
         ::Error("TClingCallFunc::exec",
               "The method %s is called without an object.",
               fMethod->Name());
         return;
      }
      vh_ary.reserve(num_args);
      vp_ary.reserve(num_args);
      for (unsigned i = 0U; i < num_args; ++i) {
         QualType Ty;
         if (i < num_params) {
            const ParmVarDecl *PVD = FD->getParamDecl(i);
            Ty = PVD->getType();
         } else {
            Ty = fArgVals[i].getType();
         }
         QualType QT = Ty.getCanonicalType();
         if (const BuiltinType *BT =
                    dyn_cast<BuiltinType>(&*QT)) {
            //
            //  WARNING!!!
            //
            //  This switch is organized in order-of-declaration
            //  so that the produced assembly code is optimal.
            //  Do not reorder!
            //
            switch (BT->getKind()) {
                  //
                  //  Builtin Types
                  //
               case BuiltinType::Void: {
                     // void
                     ::Error("TClingCallFunc::exec(void*)",
                           "Invalid type 'Void'!");
                     return;
                  }
                  break;
                  //
                  //  Unsigned Types
                  //
               case BuiltinType::Bool: {
                     // bool
                     ValHolder vh;
                     vh.u.b = (bool) sv_to_ulong_long(fArgVals[i]);
                     vh_ary.push_back(vh);
                     vp_ary.push_back(&vh_ary.back());
                  }
                  break;
               case BuiltinType::Char_U: {
                     // char on targets where it is unsigned
                     ValHolder vh;
                     vh.u.c = (char) sv_to_ulong_long(fArgVals[i]);
                     vh_ary.push_back(vh);
                     vp_ary.push_back(&vh_ary.back());
                  }
                  break;
               case BuiltinType::UChar: {
                     // unsigned char
                     ValHolder vh;
                     vh.u.uc = (unsigned char) sv_to_ulong_long(fArgVals[i]);
                     vh_ary.push_back(vh);
                     vp_ary.push_back(&vh_ary.back());
                  }
                  break;
               case BuiltinType::WChar_U: {
                     // wchar_t on targets where it is unsigned.
                     // The standard doesn't allow to specify signednedd of wchar_t
                     // thus this maps simply to wchar_t.
                     ValHolder vh;
                     vh.u.wc = (wchar_t) sv_to_ulong_long(fArgVals[i]);
                     vh_ary.push_back(vh);
                     vp_ary.push_back(&vh_ary.back());
                  }
                  break;
               case BuiltinType::Char16: {
                     // char16_t
                     //ValHolder vh;
                     //vh.u.c16 = (char16_t) sv_to_ulong_long(fArgVals[i]);
                     //vh_ary.push_back(vh);
                     //vp_ary.push_back(&vh_ary.back());
                  }
                  break;
               case BuiltinType::Char32: {
                     // char32_t
                     //ValHolder vh;
                     //vh.u.c32 = (char32_t) sv_to_ulong_long(fArgVals[i]);
                     //vh_ary.push_back(vh);
                     //vp_ary.push_back(&vh_ary.back());
                  }
                  break;
               case BuiltinType::UShort: {
                     // unsigned short
                     ValHolder vh;
                     vh.u.us = (unsigned short) sv_to_ulong_long(fArgVals[i]);
                     vh_ary.push_back(vh);
                     vp_ary.push_back(&vh_ary.back());
                  }
                  break;
               case BuiltinType::UInt: {
                     // unsigned int
                     ValHolder vh;
                     vh.u.ui = (unsigned int) sv_to_ulong_long(fArgVals[i]);
                     vh_ary.push_back(vh);
                     vp_ary.push_back(&vh_ary.back());
                  }
                  break;
               case BuiltinType::ULong: {
                     // unsigned long
                     ValHolder vh;
                     vh.u.ul = (unsigned long) sv_to_ulong_long(fArgVals[i]);
                     vh_ary.push_back(vh);
                     vp_ary.push_back(&vh_ary.back());
                  }
                  break;
               case BuiltinType::ULongLong: {
                     // unsigned long long
                     ValHolder vh;
                     vh.u.ull = (unsigned long long) sv_to_ulong_long(fArgVals[i]);
                     vh_ary.push_back(vh);
                     vp_ary.push_back(&vh_ary.back());
                  }
                  break;
               case BuiltinType::UInt128: {
                     // __uint128_t
                  }
                  break;
                  //
                  //  Signed Types
                  //
                  //
                  //  Signed Types
                  //
               case BuiltinType::Char_S: {
                     // char on targets where it is signed
                     ValHolder vh;
                     vh.u.c = (char) sv_to_long_long(fArgVals[i]);
                     vh_ary.push_back(vh);
                     vp_ary.push_back(&vh_ary.back());
                  }
                  break;
               case BuiltinType::SChar: {
                     // signed char
                     ValHolder vh;
                     vh.u.sc = (signed char) sv_to_long_long(fArgVals[i]);
                     vh_ary.push_back(vh);
                     vp_ary.push_back(&vh_ary.back());
                  }
                  break;
               case BuiltinType::WChar_S: {
                     // wchar_t on targets where it is signed.
                     // The standard doesn't allow to specify signednedd of wchar_t
                     // thus this maps simply to wchar_t.
                     ValHolder vh;
                     vh.u.wc = (wchar_t) sv_to_long_long(fArgVals[i]);
                     vh_ary.push_back(vh);
                     vp_ary.push_back(&vh_ary.back());
                  }
                  break;
               case BuiltinType::Short: {
                     // short
                     ValHolder vh;
                     vh.u.s = (short) sv_to_long_long(fArgVals[i]);
                     vh_ary.push_back(vh);
                     vp_ary.push_back(&vh_ary.back());
                  }
                  break;
               case BuiltinType::Int: {
                     // int
                     ValHolder vh;
                     vh.u.i = (int) sv_to_long_long(fArgVals[i]);
                     vh_ary.push_back(vh);
                     vp_ary.push_back(&vh_ary.back());
                  }
                  break;
               case BuiltinType::Long: {
                     // long
                     ValHolder vh;
                     vh.u.l = (long) sv_to_long_long(fArgVals[i]);
                     vh_ary.push_back(vh);
                     vp_ary.push_back(&vh_ary.back());
                  }
                  break;
               case BuiltinType::LongLong: {
                     // long long
                     ValHolder vh;
                     vh.u.ll = (long long) sv_to_long_long(fArgVals[i]);
                     vh_ary.push_back(vh);
                     vp_ary.push_back(&vh_ary.back());
                  }
                  break;
               case BuiltinType::Int128: {
                     // __int128_t
                     ::Error("TClingCallFunc::exec(void*)",
                           "Invalid type 'Int128'!");
                     return;
                  }
                  break;
               case BuiltinType::Half: {
                     // half in OpenCL, __fp16 in ARM NEON
                     ::Error("TClingCallFunc::exec(void*)",
                           "Invalid type 'Half'!");
                     return;
                  }
                  break;
               case BuiltinType::Float: {
                     // float
                     ValHolder vh;
                     vh.u.flt = sv_to<float>(fArgVals[i]);
                     vh_ary.push_back(vh);
                     vp_ary.push_back(&vh_ary.back());
                  }
                  break;
               case BuiltinType::Double: {
                     // double
                     ValHolder vh;
                     vh.u.dbl = sv_to<double>(fArgVals[i]);
                     vh_ary.push_back(vh);
                     vp_ary.push_back(&vh_ary.back());
                  }
                  break;
               case BuiltinType::LongDouble: {
                     // long double
                     ValHolder vh;
                     vh.u.ldbl = sv_to<long double>(fArgVals[i]);
                     vh_ary.push_back(vh);
                     vp_ary.push_back(&vh_ary.back());
                  }
                  break;
                  //
                  //  Language-Specific Types
                  //
               case BuiltinType::NullPtr: {
                     // C++11 nullptr
                     ValHolder vh;
                     vh.u.vp = fArgVals[i].getPtr();
                     vh_ary.push_back(vh);
                     vp_ary.push_back(&vh_ary.back());
                  }
                  break;
               default: {
                     // There should be no others.  This is here in case
                     // this changes in the future.
                     ::Error("TClingCallFunc::exec(void*)",
                           "Unhandled builtin type!");
                     QT->dump();
                     return;
                  }
                  break;
            }
         } else if (QT->isReferenceType()) {
            // the argument is already a pointer value (point to the same thing
            // as the reference.
            vp_ary.push_back((void *) sv_to_ulong_long(fArgVals[i]));
         } else if (QT->isPointerType() || QT->isArrayType()) {
            ValHolder vh;
            vh.u.vp = (void *) sv_to_ulong_long(fArgVals[i]);
            vh_ary.push_back(vh);
            vp_ary.push_back(&vh_ary.back());
         } else if (QT->isRecordType()) {
            // the argument is already a pointer value (pointing to object passed
            // by value).
            vp_ary.push_back((void *) sv_to_ulong_long(fArgVals[i]));
         } else if (const EnumType *ET =
                    dyn_cast<EnumType>(&*QT)) {
            // Note: We may need to worry about the underlying type
            //       of the enum here.
            (void) ET;
            ValHolder vh;
            vh.u.i = (int) sv_to_long_long(fArgVals[i]);
            vh_ary.push_back(vh);
            vp_ary.push_back(&vh_ary.back());
         } else if (QT->isMemberPointerType()) {
            ValHolder vh;
            vh.u.vp = (void *) sv_to_ulong_long(fArgVals[i]);
            vh_ary.push_back(vh);
            vp_ary.push_back(&vh_ary.back());
         } else {
            ::Error("TClingCallFunc::exec(void*)",
                  "Invalid type (unrecognized)!");
            QT->dump();
            return;
         }
      }
   } // End of scope holding the lock
   (*fWrapper)(address, (int)num_args, (void **)vp_ary.data(), ret);
}

template <typename T>
void TClingCallFunc::execWithLL(void *address, cling::Value *val)
{
   T ret; // leave uninit for valgrind's sake!
   exec(address, &ret);
   val->getLL() = ret;
}

template <typename T>
void TClingCallFunc::execWithULL(void *address, cling::Value *val)
{
   T ret; // leave uninit for valgrind's sake!
   exec(address, &ret);
   val->getULL() = ret;
}

// Handle integral types.
template <class T>
TClingCallFunc::ExecWithRetFunc_t TClingCallFunc::InitRetAndExecIntegral(QualType QT, cling::Value &ret)
{
   ret = cling::Value::Create<T>(QT.getAsOpaquePtr(), *fInterp);
   static_assert(std::is_integral<T>::value, "Must be called with integral T");
   if (std::is_signed<T>::value)
      return [this](void* address, cling::Value& ret) { execWithLL<T>(address, &ret); };
   else
      return [this](void* address, cling::Value& ret) { execWithULL<T>(address, &ret); };
}

// Handle builtin types.
TClingCallFunc::ExecWithRetFunc_t
TClingCallFunc::InitRetAndExecBuiltin(QualType QT, const clang::BuiltinType *BT, cling::Value &ret) {
   switch (BT->getKind()) {
      case BuiltinType::Void: {
         ret = cling::Value::Create<void>(QT.getAsOpaquePtr(), *fInterp);
         return [this](void* address, cling::Value& ret) { exec(address, 0); };
         break;
      }

         //
         //  Unsigned Types
         //
      case BuiltinType::Bool:
         return InitRetAndExecIntegral<bool>(QT, ret);
         break;
      case BuiltinType::Char_U: // char on targets where it is unsigned
      case BuiltinType::UChar:
         return InitRetAndExecIntegral<char>(QT, ret);
         break;
      case BuiltinType::WChar_U:
         return InitRetAndExecIntegral<wchar_t>(QT, ret);
         break;
      case BuiltinType::Char16:
         ::Error("TClingCallFunc::InitRetAndExecBuiltin",
               "Invalid type 'char16_t'!");
         return {};
         break;
      case BuiltinType::Char32:
         ::Error("TClingCallFunc::InitRetAndExecBuiltin",
               "Invalid type 'char32_t'!");
         return {};
         break;
      case BuiltinType::UShort:
         return InitRetAndExecIntegral<unsigned short>(QT, ret);
         break;
      case BuiltinType::UInt:
         return InitRetAndExecIntegral<unsigned int>(QT, ret);
         break;
      case BuiltinType::ULong:
         return InitRetAndExecIntegral<unsigned long>(QT, ret);
         break;
      case BuiltinType::ULongLong:
         return InitRetAndExecIntegral<unsigned long long>(QT, ret);
         break;
      case BuiltinType::UInt128: {
            ::Error("TClingCallFunc::InitRetAndExecBuiltin",
                  "Invalid type '__uint128_t'!");
            return {};
         }
         break;

         //
         //  Signed Types
         //
      case BuiltinType::Char_S: // char on targets where it is signed
      case BuiltinType::SChar:
         return InitRetAndExecIntegral<signed char>(QT, ret);
         break;
      case BuiltinType::WChar_S:
         // wchar_t on targets where it is signed.
         // The standard doesn't allow to specify signednedd of wchar_t
         // thus this maps simply to wchar_t.
         return InitRetAndExecIntegral<wchar_t>(QT, ret);
         break;
      case BuiltinType::Short:
         return InitRetAndExecIntegral<short>(QT, ret);
         break;
      case BuiltinType::Int:
         return InitRetAndExecIntegral<int>(QT, ret);
         break;
      case BuiltinType::Long:
         return InitRetAndExecIntegral<long>(QT, ret);
         break;
      case BuiltinType::LongLong:
         return InitRetAndExecIntegral<long long>(QT, ret);
         break;
      case BuiltinType::Int128:
         ::Error("TClingCallFunc::InitRetAndExecBuiltin",
               "Invalid type '__int128_t'!");
         return {};
         break;
      case BuiltinType::Half:
         // half in OpenCL, __fp16 in ARM NEON
         ::Error("TClingCallFunc::InitRetAndExecBuiltin",
               "Invalid type 'Half'!");
         return {};
         break;
      case BuiltinType::Float: {
         ret = cling::Value::Create<float>(QT.getAsOpaquePtr(), *fInterp);
         return [this](void* address, cling::Value& ret) { exec(address, &ret.getFloat()); };
         break;
      }
      case BuiltinType::Double: {
         ret = cling::Value::Create<double>(QT.getAsOpaquePtr(), *fInterp);
         return [this](void* address, cling::Value& ret) { exec(address, &ret.getDouble()); };
         break;
      }
      case BuiltinType::LongDouble: {
         ret = cling::Value::Create<long double>(QT.getAsOpaquePtr(), *fInterp);
         return [this](void* address, cling::Value& ret) { exec(address, &ret.getLongDouble()); };
         break;
      }
         //
         //  Language-Specific Types
         //
      case BuiltinType::NullPtr:
         // C++11 nullptr
         ::Error("TClingCallFunc::InitRetAndExecBuiltin",
               "Invalid type 'nullptr'!");
         return {};
         break;
      default:
         break;
   }
   return {};
}


TClingCallFunc::ExecWithRetFunc_t
TClingCallFunc::InitRetAndExecNoCtor(clang::QualType QT, cling::Value &ret) {
   if (QT->isReferenceType()) {
      ret = cling::Value(QT, *fInterp);
      return [this](void* address, cling::Value& ret) { exec(address, &ret.getPtr()); };
   } else if (QT->isMemberPointerType()) {
      const MemberPointerType *MPT = QT->getAs<MemberPointerType>();
      if (MPT->isMemberDataPointer()) {
         // A member data pointer is a actually a struct with one
         // member of ptrdiff_t, the offset from the base of the object
         // storage to the storage for the designated data member.
         // But that's not relevant: we use it as a non-builtin, allocated
         // type.
         ret = cling::Value(QT, *fInterp);
         return [this](void* address, cling::Value& ret) { exec(address, ret.getPtr()); };
      }
      // We are a function member pointer.
      ret = cling::Value(QT, *fInterp);
      return [this](void* address, cling::Value& ret) { exec(address, &ret.getPtr()); };
   } else if (QT->isPointerType() || QT->isArrayType()) {
      // Note: ArrayType is an illegal function return value type.
      ret = cling::Value::Create<void*>(QT.getAsOpaquePtr(), *fInterp);
      return [this](void* address, cling::Value& ret) { exec(address, &ret.getPtr()); };
   } else if (QT->isRecordType()) {
      ret = cling::Value(QT, *fInterp);
      return [this](void* address, cling::Value& ret) { exec(address, ret.getPtr()); };
   } else if (const EnumType *ET = dyn_cast<EnumType>(&*QT)) {
      // Note: We may need to worry about the underlying type
      //       of the enum here.
      (void) ET;
      ret = cling::Value(QT, *fInterp);
      return [this](void* address, cling::Value& ret) { execWithLL<int>(address, &ret); };
   } else if (const BuiltinType *BT = dyn_cast<BuiltinType>(&*QT)) {
      return InitRetAndExecBuiltin(QT, BT, ret);
   }
   ::Error("TClingCallFunc::exec_with_valref_return",
         "Unrecognized return type!");
   QT->dump();
   return {};
}

TClingCallFunc::ExecWithRetFunc_t
TClingCallFunc::InitRetAndExec(const clang::FunctionDecl *FD, cling::Value &ret) {
   if (llvm::isa<CXXConstructorDecl>(FD)) {
      ASTContext &Context = FD->getASTContext();
      const TypeDecl *TD = dyn_cast<TypeDecl>(GetDeclContext());
      QualType ClassTy(TD->getTypeForDecl(), 0);
      QualType QT = Context.getLValueReferenceType(ClassTy);
      ret = cling::Value(QT, *fInterp);
      // Store the new()'ed address in getPtr()
      return [this](void* address, cling::Value& ret) { exec(address, &ret.getPtr()); };
   } else {
      QualType QT = FD->getReturnType().getCanonicalType();
      return InitRetAndExecNoCtor(QT, ret);
   }
}

void TClingCallFunc::exec_with_valref_return(void *address, cling::Value *ret)
{
   if (!ret) {
      exec(address, 0);
      return;
   }
   std::function<void(void*, cling::Value&)> execFunc;

   /* Release lock during user function execution*/
   {
      R__LOCKGUARD_CLING(gInterpreterMutex);
      execFunc = InitRetAndExec(GetDecl(), *ret);
   }

   if (execFunc)
      execFunc(address, *ret);
   return;
}

void TClingCallFunc::EvaluateArgList(const string &ArgList)
{
   R__LOCKGUARD_CLING(gInterpreterMutex);

   SmallVector<Expr *, 4> exprs;
   fInterp->getLookupHelper().findArgList(ArgList, exprs,
                                          gDebug > 5 ? cling::LookupHelper::WithDiagnostics
                                          : cling::LookupHelper::NoDiagnostics);
   for (SmallVectorImpl<Expr *>::const_iterator I = exprs.begin(),
         E = exprs.end(); I != E; ++I) {
      cling::Value val;
      EvaluateExpr(*fInterp, *I, val);
      if (!val.isValid()) {
         // Bad expression, all done.
         ::Error("TClingCallFunc::EvaluateArgList",
               "Bad expression in parameter %d of '%s'!",
               (int)(I - exprs.begin()),
               ArgList.c_str());
         return;
      }
      fArgVals.push_back(val);
   }
}

void TClingCallFunc::Exec(void *address, TInterpreterValue *interpVal/*=0*/)
{
   IFacePtr();
   if (!fWrapper) {
      ::Error("TClingCallFunc::Exec(address, interpVal)",
            "Called with no wrapper, not implemented!");
      return;
   }
   if (!interpVal) {
      exec(address, 0);
      return;
   }
   cling::Value *val = reinterpret_cast<cling::Value *>(interpVal->GetValAddr());
   exec_with_valref_return(address, val);
}

template <typename T>
T TClingCallFunc::ExecT(void *address)
{
   IFacePtr();
   if (!fWrapper) {
      ::Error("TClingCallFunc::ExecT",
            "Called with no wrapper, not implemented!");
      return 0;
   }
   cling::Value ret;
   exec_with_valref_return(address, &ret);
   if (!ret.isValid()) {
      // Sometimes we are called on a function returning void!
      return 0;
   }

   if (fReturnIsRecordType)
      ((TCling *)gCling)->RegisterTemporary(ret);
   return sv_to<T>(ret);
}

Long_t TClingCallFunc::ExecInt(void *address)
{
   return ExecT<long>(address);
}

long long TClingCallFunc::ExecInt64(void *address)
{
   return ExecT<long long>(address);
}

double TClingCallFunc::ExecDouble(void *address)
{
   return ExecT<double>(address);
}

void TClingCallFunc::ExecWithArgsAndReturn(void *address, const void *args[] /*= 0*/,
      int nargs /*= 0*/, void *ret/*= 0*/)
{
   IFacePtr();
   if (!fWrapper) {
      ::Error("TClingCallFunc::ExecWithArgsAndReturn(address, args, ret)",
            "Called with no wrapper, not implemented!");
      return;
   }
   (*fWrapper)(address, nargs, const_cast<void **>(args), ret);
}

void TClingCallFunc::ExecWithReturn(void *address, void *ret/*= 0*/)
{
   IFacePtr();
   if (!fWrapper) {
      ::Error("TClingCallFunc::ExecWithReturn(address, ret)",
            "Called with no wrapper, not implemented!");
      return;
   }
   exec(address, ret);
}

void *TClingCallFunc::ExecDefaultConstructor(const TClingClassInfo *info,
                                             ROOT::TMetaUtils::EIOCtorCategory kind,
                                             const std::string &type_name,
                                             void *address /*=0*/, unsigned long nary /*= 0UL*/)
{
   if (!info->IsValid()) {
      ::Error("TClingCallFunc::ExecDefaultConstructor", "Invalid class info!");
      return nullptr;
   }
   tcling_callfunc_ctor_Wrapper_t wrapper = nullptr;
   {
      R__LOCKGUARD_CLING(gInterpreterMutex);
      auto D = info->GetDecl();
      //if (!info->HasDefaultConstructor()) {
      //   // FIXME: We might have a ROOT ioctor, we might
      //   //        have to check for that here.
      //   ::Error("TClingCallFunc::ExecDefaultConstructor",
      //         "Class has no default constructor: %s",
      //         info->Name());
      //   return 0;
      //}
      auto I = gCtorWrapperStore.find(D);
      if (I != gCtorWrapperStore.end()) {
         wrapper = (tcling_callfunc_ctor_Wrapper_t) I->second;
      } else {
         wrapper = make_ctor_wrapper(info, kind, type_name);
      }
   }
   if (!wrapper) {
      ::Error("TClingCallFunc::ExecDefaultConstructor",
            "Called with no wrapper, not implemented!");
      return nullptr;
   }
   void *obj = 0;
   (*wrapper)(&obj, address, nary);
   return obj;
}

void TClingCallFunc::ExecDestructor(const TClingClassInfo *info, void *address /*=0*/,
                                    unsigned long nary /*= 0UL*/, bool withFree /*= true*/)
{
   if (!info->IsValid()) {
      ::Error("TClingCallFunc::ExecDestructor", "Invalid class info!");
      return;
   }

   tcling_callfunc_dtor_Wrapper_t wrapper = 0;
   {
      R__LOCKGUARD_CLING(gInterpreterMutex);
      const Decl *D = info->GetDecl();
      map<const Decl *, void *>::iterator I = gDtorWrapperStore.find(D);
      if (I != gDtorWrapperStore.end()) {
         wrapper = (tcling_callfunc_dtor_Wrapper_t) I->second;
      } else {
         wrapper = make_dtor_wrapper(info);
      }
   }
   if (!wrapper) {
      ::Error("TClingCallFunc::ExecDestructor",
            "Called with no wrapper, not implemented!");
      return;
   }
   (*wrapper)(address, nary, withFree);
}

TClingMethodInfo *
TClingCallFunc::FactoryMethod() const
{
   return new TClingMethodInfo(*fMethod);
}

void TClingCallFunc::Init()
{
   fMethod.reset();
   fWrapper = 0;
   fDecl = nullptr;
   fMinRequiredArguments = -1;
   ResetArg();
}

void TClingCallFunc::Init(const TClingMethodInfo &minfo)
{
   Init();
   fMethod = std::unique_ptr<TClingMethodInfo>(new TClingMethodInfo(minfo));
}

void TClingCallFunc::Init(std::unique_ptr<TClingMethodInfo> minfo)
{
   Init();
   fMethod = std::move(minfo);
}

void *TClingCallFunc::InterfaceMethod()
{
   if (!IsValid()) {
      return 0;
   }
   if (!fWrapper) {
      const Decl *decl = GetFunctionOrShadowDecl();

      R__LOCKGUARD_CLING(gInterpreterMutex);
      map<const Decl *, void *>::iterator I = gWrapperStore.find(decl);
      if (I != gWrapperStore.end()) {
         fWrapper = (tcling_callfunc_Wrapper_t) I->second;
      } else {
         fWrapper = make_wrapper();
      }
   }
   return (void *)fWrapper;
}

bool TClingCallFunc::IsValid() const
{
   if (!fMethod) {
      return false;
   }
   return fMethod->IsValid();
}

TInterpreter::CallFuncIFacePtr_t TClingCallFunc::IFacePtr()
{
   if (!IsValid()) {
      ::Error("TClingCallFunc::IFacePtr(kind)",
            "Attempt to get interface while invalid.");
      return TInterpreter::CallFuncIFacePtr_t();
   }
   if (!fWrapper) {
      const Decl *decl = GetFunctionOrShadowDecl();

      R__LOCKGUARD_CLING(gInterpreterMutex);
      map<const Decl *, void *>::iterator I = gWrapperStore.find(decl);
      if (I != gWrapperStore.end()) {
         fWrapper = (tcling_callfunc_Wrapper_t) I->second;
      } else {
         fWrapper = make_wrapper();
      }

      fReturnIsRecordType = GetDecl()->getReturnType().getCanonicalType()->isRecordType();
   }
   return TInterpreter::CallFuncIFacePtr_t(fWrapper);
}


void TClingCallFunc::ResetArg()
{
   fArgVals.clear();
}

void TClingCallFunc::SetArg(unsigned long param)
{
   const ASTContext &C = fInterp->getCI()->getASTContext();
   fArgVals.push_back(cling::Value(C.UnsignedLongTy, *fInterp));
   fArgVals.back().getLL() = param;
}

void TClingCallFunc::SetArg(long param)
{
   const ASTContext &C = fInterp->getCI()->getASTContext();
   fArgVals.push_back(cling::Value(C.LongTy, *fInterp));
   fArgVals.back().getLL() = param;
}

void TClingCallFunc::SetArg(float param)
{
   const ASTContext &C = fInterp->getCI()->getASTContext();
   fArgVals.push_back(cling::Value(C.FloatTy, *fInterp));
   fArgVals.back().getFloat() = param;
}

void TClingCallFunc::SetArg(double param)
{
   const ASTContext &C = fInterp->getCI()->getASTContext();
   fArgVals.push_back(cling::Value(C.DoubleTy, *fInterp));
   fArgVals.back().getDouble() = param;
}

void TClingCallFunc::SetArg(long long param)
{
   const ASTContext &C = fInterp->getCI()->getASTContext();
   fArgVals.push_back(cling::Value(C.LongLongTy, *fInterp));
   fArgVals.back().getLL() = param;
}

void TClingCallFunc::SetArg(unsigned long long param)
{
   const ASTContext &C = fInterp->getCI()->getASTContext();
   fArgVals.push_back(cling::Value(C.UnsignedLongLongTy, *fInterp));
   fArgVals.back().getULL() = param;
}

void TClingCallFunc::SetArgArray(long *paramArr, int nparam)
{
   ResetArg();
   for (int i = 0; i < nparam; ++i) {
      SetArg(paramArr[i]);
   }
}

void TClingCallFunc::SetArgs(const char *params)
{
   ResetArg();
   EvaluateArgList(params);
}

void TClingCallFunc::SetFunc(const TClingClassInfo *info, const char *method, const char *arglist,
                             long *poffset)
{
   SetFunc(info, method, arglist, false, poffset);
}

void TClingCallFunc::SetFunc(const TClingClassInfo *info, const char *method, const char *arglist,
                             bool objectIsConst, long *poffset)
{
   Init(std::unique_ptr<TClingMethodInfo>(new TClingMethodInfo(fInterp)));
   if (poffset) {
      *poffset = 0L;
   }
   ResetArg();
   if (!info->IsValid()) {
      ::Error("TClingCallFunc::SetFunc", "Class info is invalid!");
      return;
   }
   if (!strcmp(arglist, ")")) {
      // CINT accepted a single right paren as meaning no arguments.
      arglist = "";
   }
   *fMethod = info->GetMethodWithArgs(method, arglist, objectIsConst, poffset);
   if (!fMethod->IsValid()) {
      //::Error("TClingCallFunc::SetFunc", "Could not find method %s(%s)", method,
      //      arglist);
      return;
   }
   // FIXME: The arglist was already parsed by the lookup, we should
   //        enhance the lookup to return the resulting expression
   //        list so we do not need to parse it again here.
   EvaluateArgList(arglist);
}

void TClingCallFunc::SetFunc(const TClingMethodInfo *info)
{
   Init(std::unique_ptr<TClingMethodInfo>(new TClingMethodInfo(*info)));
   ResetArg();
   if (!fMethod->IsValid()) {
      return;
   }
}

void TClingCallFunc::SetFuncProto(const TClingClassInfo *info, const char *method,
                                  const char *proto, long *poffset,
                                  EFunctionMatchMode mode/*=kConversionMatch*/)
{
   SetFuncProto(info, method, proto, false, poffset, mode);
}

void TClingCallFunc::SetFuncProto(const TClingClassInfo *info, const char *method,
                                  const char *proto, bool objectIsConst, long *poffset,
                                  EFunctionMatchMode mode/*=kConversionMatch*/)
{
   Init(std::unique_ptr<TClingMethodInfo>(new TClingMethodInfo(fInterp)));
   if (poffset) {
      *poffset = 0L;
   }
   ResetArg();
   if (!info->IsValid()) {
      ::Error("TClingCallFunc::SetFuncProto", "Class info is invalid!");
      return;
   }
   *fMethod = info->GetMethod(method, proto, objectIsConst, poffset, mode);
   if (!fMethod->IsValid()) {
      //::Error("TClingCallFunc::SetFuncProto", "Could not find method %s(%s)",
      //      method, proto);
      return;
   }
}

void TClingCallFunc::SetFuncProto(const TClingClassInfo *info, const char *method,
                                  const llvm::SmallVectorImpl<clang::QualType> &proto, long *poffset,
                                  EFunctionMatchMode mode/*=kConversionMatch*/)
{
   SetFuncProto(info, method, proto, false, poffset, mode);
}

void TClingCallFunc::SetFuncProto(const TClingClassInfo *info, const char *method,
                                  const llvm::SmallVectorImpl<clang::QualType> &proto,
                                  bool objectIsConst, long *poffset,
                                  EFunctionMatchMode mode/*=kConversionMatch*/)
{
   Init(std::unique_ptr<TClingMethodInfo>(new TClingMethodInfo(fInterp)));
   if (poffset) {
      *poffset = 0L;
   }
   ResetArg();
   if (!info->IsValid()) {
      ::Error("TClingCallFunc::SetFuncProto", "Class info is invalid!");
      return;
   }
   *fMethod = info->GetMethod(method, proto, objectIsConst, poffset, mode);
   if (!fMethod->IsValid()) {
      //::Error("TClingCallFunc::SetFuncProto", "Could not find method %s(%s)",
      //      method, proto);
      return;
   }
}

