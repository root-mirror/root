//------------------------------------------------------------------------------
// CLING - the C++ LLVM-based InterpreterG :)
//
// This file is dual-licensed: you can choose to license it under the University
// of Illinois Open Source License or the GNU Lesser General Public License. See
// LICENSE.TXT for details.
//------------------------------------------------------------------------------

// RUN: cat %s | %cling -I%p | FileCheck %s

// This test covers cling::utils::GetFullyQualifiedName() interfaces.
// https://github.com/root-project/root/issues/7955


#include "cling/Interpreter/Interpreter.h"
#include "cling/Interpreter/LookupHelper.h"
#include "cling/Utils/AST.h"

#include "clang/AST/DeclCXX.h"
#include "clang/AST/DeclTemplate.h"

#include "llvm/Support/Casting.h"

.rawInput 1

template <class T>
struct Template1 {
  using type1 = int;
  void f0(type1);
};

Template1<float> xf1;

.rawInput 0

const auto& lookup = gCling->getLookupHelper();
using namespace clang;

auto *T1s = llvm::dyn_cast<ClassTemplateSpecializationDecl>(lookup.findScope("Template1<short>", cling::LookupHelper::WithDiagnostics));
// Find first function decl:
auto iFunc = std::find_if(T1s->decls_begin(), T1s->decls_end(),
  [](const Decl* D) {if (auto FD = llvm::dyn_cast<CXXMethodDecl>(D)) return FD->isUserProvided(); return false;});
using cling::utils::TypeName::GetFullyQualifiedName;
GetFullyQualifiedName(llvm::dyn_cast<FunctionDecl>(*iFunc)->getParamDecl(0)->getType(), T1s->getASTContext(), T1s) // CHECK: (std::string) "Template1<short>::type1"
