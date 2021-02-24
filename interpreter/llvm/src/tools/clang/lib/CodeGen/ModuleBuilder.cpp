//===--- ModuleBuilder.cpp - Emit LLVM Code from ASTs ---------------------===//
//
// Part of the LLVM Project, under the Apache License v2.0 with LLVM Exceptions.
// See https://llvm.org/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//
//
// This builds an AST and converts it to LLVM Code.
//
//===----------------------------------------------------------------------===//

#include "clang/CodeGen/ModuleBuilder.h"
#include "CGDebugInfo.h"
#include "CGCXXABI.h"
#include "CodeGenModule.h"
#include "clang/AST/ASTContext.h"
#include "clang/AST/DeclObjC.h"
#include "clang/AST/Expr.h"
#include "clang/Basic/CodeGenOptions.h"
#include "clang/Basic/Diagnostic.h"
#include "clang/Basic/TargetInfo.h"
#include "llvm/ADT/StringRef.h"
#include "llvm/IR/DataLayout.h"
#include "llvm/IR/LLVMContext.h"
#include "llvm/IR/Module.h"
#include <memory>

using namespace clang;
using namespace CodeGen;

namespace {
  struct CXXABICtxSwapper: clang::CodeGen::CGCXXABI {
    void SwapCtx(clang::CodeGen::CGCXXABI &other) {
      std::swap(MangleCtx, ((CXXABICtxSwapper&)other).MangleCtx);
    }
  };
}

namespace clang {
  class CodeGeneratorImpl : public CodeGenerator {
    DiagnosticsEngine &Diags;
    ASTContext *Ctx;
    const HeaderSearchOptions &HeaderSearchOpts; // Only used for debug info.
    const PreprocessorOptions &PreprocessorOpts; // Only used for debug info.
    CodeGenOptions CodeGenOpts;  // Intentionally copied in.

    unsigned HandlingTopLevelDecls;

    /// Use this when emitting decls to block re-entrant decl emission. It will
    /// emit all deferred decls on scope exit. Set EmitDeferred to false if decl
    /// emission must be deferred longer, like at the end of a tag definition.
    struct HandlingTopLevelDeclRAII {
      CodeGeneratorImpl &Self;
      bool EmitDeferred;
      HandlingTopLevelDeclRAII(CodeGeneratorImpl &Self,
                               bool EmitDeferred = true)
          : Self(Self), EmitDeferred(EmitDeferred) {
        ++Self.HandlingTopLevelDecls;
      }
      ~HandlingTopLevelDeclRAII() {
        unsigned Level = --Self.HandlingTopLevelDecls;
        if (Level == 0 && EmitDeferred)
          Self.EmitDeferredDecls();
      }
    };

    CoverageSourceInfo *CoverageInfo;

  protected:
    std::unique_ptr<llvm::Module> M;
    std::unique_ptr<CodeGen::CodeGenModule> Builder;

  private:
    SmallVector<FunctionDecl *, 8> DeferredInlineMemberFuncDefs;

  public:
    CodeGeneratorImpl(DiagnosticsEngine &diags, llvm::StringRef ModuleName,
                      const HeaderSearchOptions &HSO,
                      const PreprocessorOptions &PPO, const CodeGenOptions &CGO,
                      llvm::LLVMContext &C,
                      CoverageSourceInfo *CoverageInfo = nullptr)
        : Diags(diags), Ctx(nullptr), HeaderSearchOpts(HSO),
          PreprocessorOpts(PPO), CodeGenOpts(CGO), HandlingTopLevelDecls(0),
          CoverageInfo(CoverageInfo), M(new llvm::Module(ModuleName, C)) {
      C.setDiscardValueNames(CGO.DiscardValueNames);
    }

    ~CodeGeneratorImpl() override {
      // There should normally not be any leftover inline method definitions.
      assert(DeferredInlineMemberFuncDefs.empty() ||
             Diags.hasErrorOccurred());
    }

    CodeGenModule &CGM() {
      return *Builder;
    }

    llvm::Module *GetModule() {
      return M.get();
    }

    CGDebugInfo *getCGDebugInfo() {
      return Builder->getModuleDebugInfo();
    }

    llvm::Module *ReleaseModule() {
      // Remove pending etc decls in case of error; the asserts in StartModule()
      // will rightfully be confused otherwise, as none of the decls were
      // emitted.
      if (Diags.hasErrorOccurred())
        Builder->clear();
      return M.release();
    }

    const Decl *GetDeclForMangledName(StringRef MangledName) {
      GlobalDecl Result;
      if (!Builder->lookupRepresentativeDecl(MangledName, Result))
        return nullptr;
      const Decl *D = Result.getCanonicalDecl().getDecl();
      if (auto FD = dyn_cast<FunctionDecl>(D)) {
        if (FD->hasBody(FD))
          return FD;
      } else if (auto TD = dyn_cast<TagDecl>(D)) {
        if (auto Def = TD->getDefinition())
          return Def;
      }
      return D;
    }

    llvm::Constant *GetAddrOfGlobal(GlobalDecl global, bool isForDefinition) {
      return Builder->GetAddrOfGlobal(global, ForDefinition_t(isForDefinition));
    }

    llvm::Module *StartModule(llvm::StringRef ModuleName,
                              llvm::LLVMContext& C,
                              const CodeGenOptions& CGO) {
      assert(!M && "Replacing existing Module?");

      std::unique_ptr<CodeGen::CodeGenModule> OldBuilder;
      OldBuilder.swap(Builder);
      CodeGenOpts = CGO;
      M.reset(new llvm::Module(ModuleName, C));
      Initialize(*Ctx);

      assert(OldBuilder->DeferredDeclsToEmit.empty()
             && "Should have emitted all decls deferred to emit.");
      assert(Builder->DeferredDecls.empty()
             && "Newly created module should not have deferred decls");
      Builder->DeferredDecls.swap(OldBuilder->DeferredDecls);

      assert(OldBuilder->EmittedDeferredDecls.empty()
             && "Still have (unmerged) EmittedDeferredDecls deferred decls");

      assert(Builder->DeferredVTables.empty()
             && "Newly created module should not have deferred vtables");
      Builder->DeferredVTables.swap(OldBuilder->DeferredVTables);

      assert(Builder->Manglings.empty()
             && "Newly created module should not have manglings");
      // Calls swap() internally, *also* swapping the Allocator object which is
      // essential to keep the storage!
      Builder->Manglings = std::move(OldBuilder->Manglings);

      assert(Builder->WeakRefReferences.empty()
             && "Newly created module should not have weakRefRefs");
      Builder->WeakRefReferences.swap(OldBuilder->WeakRefReferences);

      ((CXXABICtxSwapper&)*Builder->ABI).SwapCtx(*OldBuilder->ABI);
      Builder->TBAA.swap(OldBuilder->TBAA);

      return M.get();
    }

    void print(llvm::raw_ostream& out) {
      out << "\n\nCodeGen:\n";
      //llvm::SmallPtrSet<llvm::GlobalValue*, 10> WeakRefReferences;
      out << " WeakRefReferences (llvm::SmallPtrSet<llvm::GlobalValue*, 10>)\n";
      for(auto I = Builder->WeakRefReferences.begin(),
            E = Builder->WeakRefReferences.end(); I != E; ++I) {
        (*I)->print(out);
        out << "\n";
      }

      //llvm::StringMap<GlobalDecl> DeferredDecls;
      out << " DeferredDecls (llvm::StringMap<GlobalDecl>)\n";
      for(auto I = Builder->DeferredDecls.begin(),
            E = Builder->DeferredDecls.end(); I != E; ++I) {
        out << I->first.str().c_str();
        I->second.getDecl()->print(out);
        out << "\n";
      }

      //std::vector<DeferredGlobal> DeferredDeclsToEmit;
      out << " DeferredDeclsToEmit (std::vector<DeferredGlobal>)\n";
      for(auto I = Builder->DeferredDeclsToEmit.begin(),
            E = Builder->DeferredDeclsToEmit.end(); I != E; ++I) {
        I->GD.getDecl()->print(out);
        I->GV->print(out);
        out << "\n";
      }

      //std::vector<GlobalDecl> Aliases;
      out << " Aliases (std::vector<GlobalDecl>)\n";
      for(auto I = Builder->Aliases.begin(),
            E = Builder->Aliases.end(); I != E; ++I) {
        I->getDecl()->print(out);
        out << "\n";
      }
      //typedef llvm::StringMap<llvm::TrackingVH<llvm::Constant> >
      // ReplacementsTy;
      //ReplacementsTy Replacements;
      out
        << " Replacements (llvm::StringMap<llvm::TrackingVH<llvm::Constant>>\n";
      for(auto I = Builder->Replacements.begin(),
            E = Builder->Replacements.end(); I != E; ++I) {
        out << I->getKey().str().c_str();
        (*I->getValue()).print(out);
        out << "\n";
      }

      //std::vector<const CXXRecordDecl*> DeferredVTables;
      out << " DeferredVTables (std::vector<const CXXRecordDecl*>\n";
      for(auto I = Builder->DeferredVTables.begin(),
            E = Builder->DeferredVTables.end(); I != E; ++I) {
        (*I)->print(out);
        out << "\n";
      }

      //std::vector<llvm::WeakVH> LLVMUsed;
      out << " LLVMUsed (std::vector<llvm::WeakVH> >\n";
      for(auto I = Builder->LLVMUsed.begin(),
            E = Builder->LLVMUsed.end(); I != E; ++I) {
        (*I)->print(out);
        out << "\n";
      }

      // typedef std::vector<std::pair<llvm::Constant*, int> > CtorList;
      //CtorList GlobalCtors;
      out << " GlobalCtors (std::vector<std::pair<llvm::Constant*, int> >\n";
      for(auto I = Builder->GlobalCtors.begin(),
            E = Builder->GlobalCtors.end(); I != E; ++I) {
        out << I->Initializer << " : " << I->AssociatedData;
        out << "\n";
      }

      //CtorList GlobalDtors;
      out << " GlobalDtors (std::vector<std::pair<llvm::Constant*, int> >\n";
      for(auto I = Builder->GlobalDtors.begin(),
            E = Builder->GlobalDtors.end(); I != E; ++I) {
        out << I->Initializer << " : " << I->AssociatedData;
        out << "\n";
      }

      //llvm::DenseMap<GlobalDecl, StringRef> MangledDeclNames;
      //std::vector<llvm::Constant*> Annotations;
      //llvm::StringMap<llvm::Constant*> AnnotationStrings;
      //llvm::StringMap<llvm::Constant*> CFConstantStringMap;
      //llvm::StringMap<llvm::GlobalVariable*> ConstantStringMap;
      out << " ConstantStringMap (llvm::DenseMap<llvm::Constant *, "
             "llvm::GlobalVariable *>)\n";
      for(auto I = Builder->ConstantStringMap.begin(),
            E = Builder->ConstantStringMap.end(); I != E; ++I) {
        I->first->print(out);
        I->second->print(out);
        out << "\n";
      }

      //llvm::DenseMap<const Decl*, llvm::Constant *> StaticLocalDeclMap;
      //llvm::DenseMap<const Decl*, llvm::GlobalVariable*> StaticLocalDeclGuardMap;
      //llvm::DenseMap<const Expr*, llvm::Constant *> MaterializedGlobalTemporaryMap;
      //llvm::DenseMap<QualType, llvm::Constant *> AtomicSetterHelperFnMap;
      //llvm::DenseMap<QualType, llvm::Constant *> AtomicGetterHelperFnMap;
      //llvm::DenseMap<QualType, llvm::Constant *> TypeDescriptorMap;
      //StaticExternCMap StaticExternCValues;
      //std::vector<std::pair<const VarDecl *, llvm::GlobalVariable *> >
      // CXXThreadLocals;
      //std::vector<llvm::Constant*> CXXThreadLocalInits;
      //std::vector<llvm::Constant*> CXXGlobalInits;
      //llvm::DenseMap<const Decl*, unsigned> DelayedCXXInitPosition;
      //SmallVector<GlobalInitData, 8> PrioritizedCXXGlobalInits;
      //std::vector<std::pair<llvm::WeakVH,llvm::Constant*> > CXXGlobalDtors;
      //llvm::SetVector<clang::Module *> ImportedModules;
      //SmallVector<llvm::Value *, 16> LinkerOptionsMetadata;
      //
      out.flush();
    }

    llvm::Module *StartModule(llvm::StringRef ModuleName,
                              llvm::LLVMContext &C) {
      assert(!M && "Replacing existing Module?");
      M.reset(new llvm::Module(ModuleName, C));
      Initialize(*Ctx);
      return M.get();
    }

    void forgetGlobal(llvm::GlobalValue* GV) {
      for (auto I = Builder->ConstantStringMap.begin(),
            E = Builder->ConstantStringMap.end(); I != E; ++I) {
        if (I->second == GV) {
          Builder->ConstantStringMap.erase(I);
          break;
        }
      }
    }

    void forgetDecl(const GlobalDecl& GD, llvm::StringRef MangledName) {
      Builder->DeferredDecls.erase(MangledName);
      Builder->Manglings.erase(MangledName);
    }

    void Initialize(ASTContext &Context) override {
      Ctx = &Context;

      M->setTargetTriple(Ctx->getTargetInfo().getTriple().getTriple());
      M->setDataLayout(Ctx->getTargetInfo().getDataLayout());
      const auto &SDKVersion = Ctx->getTargetInfo().getSDKVersion();
      if (!SDKVersion.empty())
        M->setSDKVersion(SDKVersion);
      Builder.reset(new CodeGen::CodeGenModule(Context, HeaderSearchOpts,
                                               PreprocessorOpts, CodeGenOpts,
                                               *M, Diags, CoverageInfo));

      for (auto &&Lib : CodeGenOpts.DependentLibraries)
        Builder->AddDependentLib(Lib);
      for (auto &&Opt : CodeGenOpts.LinkerOptions)
        Builder->AppendLinkerOptions(Opt);
    }

    void HandleCXXStaticMemberVarInstantiation(VarDecl *VD) override {
      if (Diags.hasErrorOccurred())
        return;

      Builder->HandleCXXStaticMemberVarInstantiation(VD);
    }

    bool HandleTopLevelDecl(DeclGroupRef DG) override {
      if (Diags.hasErrorOccurred())
        return true;

      HandlingTopLevelDeclRAII HandlingDecl(*this);

      // Make sure to emit all elements of a Decl.
      for (DeclGroupRef::iterator I = DG.begin(), E = DG.end(); I != E; ++I)
        Builder->EmitTopLevelDecl(*I);

      return true;
    }

    void EmitDeferredDecls() {
      if (DeferredInlineMemberFuncDefs.empty())
        return;

      // Emit any deferred inline method definitions. Note that more deferred
      // methods may be added during this loop, since ASTConsumer callbacks
      // can be invoked if AST inspection results in declarations being added.
      HandlingTopLevelDeclRAII HandlingDecl(*this);
      for (unsigned I = 0; I != DeferredInlineMemberFuncDefs.size(); ++I)
        Builder->EmitTopLevelDecl(DeferredInlineMemberFuncDefs[I]);
      DeferredInlineMemberFuncDefs.clear();
    }

    void HandleInlineFunctionDefinition(FunctionDecl *D) override {
      if (Diags.hasErrorOccurred())
        return;

      assert(D->doesThisDeclarationHaveABody());

      // We may want to emit this definition. However, that decision might be
      // based on computing the linkage, and we have to defer that in case we
      // are inside of something that will change the method's final linkage,
      // e.g.
      //   typedef struct {
      //     void bar();
      //     void foo() { bar(); }
      //   } A;
      DeferredInlineMemberFuncDefs.push_back(D);

      // Provide some coverage mapping even for methods that aren't emitted.
      // Don't do this for templated classes though, as they may not be
      // instantiable.
      if (!D->getLexicalDeclContext()->isDependentContext())
        Builder->AddDeferredUnusedCoverageMapping(D);
    }

    /// HandleTagDeclDefinition - This callback is invoked each time a TagDecl
    /// to (e.g. struct, union, enum, class) is completed. This allows the
    /// client hack on the type, which can occur at any point in the file
    /// (because these can be defined in declspecs).
    void HandleTagDeclDefinition(TagDecl *D) override {
      if (Diags.hasErrorOccurred())
        return;

      // Don't allow re-entrant calls to CodeGen triggered by PCH
      // deserialization to emit deferred decls.
      HandlingTopLevelDeclRAII HandlingDecl(*this, /*EmitDeferred=*/false);

      Builder->UpdateCompletedType(D);

      // For MSVC compatibility, treat declarations of static data members with
      // inline initializers as definitions.
      if (Ctx->getTargetInfo().getCXXABI().isMicrosoft()) {
        for (Decl *Member : D->decls()) {
          if (VarDecl *VD = dyn_cast<VarDecl>(Member)) {
            if (Ctx->isMSStaticDataMemberInlineDefinition(VD) &&
                Ctx->DeclMustBeEmitted(VD)) {
              Builder->EmitGlobal(VD);
            }
          }
        }
      }
      // For OpenMP emit declare reduction functions, if required.
      if (Ctx->getLangOpts().OpenMP) {
        for (Decl *Member : D->decls()) {
          if (auto *DRD = dyn_cast<OMPDeclareReductionDecl>(Member)) {
            if (Ctx->DeclMustBeEmitted(DRD))
              Builder->EmitGlobal(DRD);
          }
        }
      }
    }

    void HandleTagDeclRequiredDefinition(const TagDecl *D) override {
      if (Diags.hasErrorOccurred())
        return;

      // Don't allow re-entrant calls to CodeGen triggered by PCH
      // deserialization to emit deferred decls.
      HandlingTopLevelDeclRAII HandlingDecl(*this, /*EmitDeferred=*/false);

      if (CodeGen::CGDebugInfo *DI = Builder->getModuleDebugInfo())
        if (const RecordDecl *RD = dyn_cast<RecordDecl>(D))
          DI->completeRequiredType(RD);
    }

    void HandleTranslationUnit(ASTContext &Ctx) override {
      // Release the Builder when there is no error.
      if (!Diags.hasErrorOccurred() && Builder)
        Builder->Release();

      // If there are errors before or when releasing the Builder, reset
      // the module to stop here before invoking the backend.
      if (Diags.hasErrorOccurred()) {
        if (Builder)
          Builder->clear();
        M.reset();
        return;
      }
    }

    void AssignInheritanceModel(CXXRecordDecl *RD) override {
      if (Diags.hasErrorOccurred())
        return;

      Builder->RefreshTypeCacheForClass(RD);
    }

    void CompleteTentativeDefinition(VarDecl *D) override {
      if (Diags.hasErrorOccurred())
        return;

      Builder->EmitTentativeDefinition(D);
    }

    void HandleVTable(CXXRecordDecl *RD) override {
      if (Diags.hasErrorOccurred())
        return;

      Builder->EmitVTable(RD);
    }
  };
}

void CodeGenerator::anchor() { }

CodeGenModule &CodeGenerator::CGM() {
  return static_cast<CodeGeneratorImpl*>(this)->CGM();
}

llvm::Module *CodeGenerator::GetModule() {
  return static_cast<CodeGeneratorImpl*>(this)->GetModule();
}

llvm::Module *CodeGenerator::ReleaseModule() {
  return static_cast<CodeGeneratorImpl*>(this)->ReleaseModule();
}

CGDebugInfo *CodeGenerator::getCGDebugInfo() {
  return static_cast<CodeGeneratorImpl*>(this)->getCGDebugInfo();
}

const Decl *CodeGenerator::GetDeclForMangledName(llvm::StringRef name) {
  return static_cast<CodeGeneratorImpl*>(this)->GetDeclForMangledName(name);
}

llvm::Constant *CodeGenerator::GetAddrOfGlobal(GlobalDecl global,
                                               bool isForDefinition) {
  return static_cast<CodeGeneratorImpl*>(this)
           ->GetAddrOfGlobal(global, isForDefinition);
}

void CodeGenerator::print(llvm::raw_ostream& out) {
  static_cast<CodeGeneratorImpl*>(this)->print(out);
}

llvm::Module *CodeGenerator::StartModule(llvm::StringRef ModuleName,
                                         llvm::LLVMContext &C) {
  return static_cast<CodeGeneratorImpl*>(this)->StartModule(ModuleName, C);
}

void CodeGenerator::forgetGlobal(llvm::GlobalValue* GV) {
  static_cast<CodeGeneratorImpl*>(this)->forgetGlobal(GV);
}

void CodeGenerator::forgetDecl(const GlobalDecl& GD,
                               llvm::StringRef MangledName) {
  static_cast<CodeGeneratorImpl*>(this)->forgetDecl(GD, MangledName);
}


llvm::Module *CodeGenerator::StartModule(llvm::StringRef ModuleName,
                                         llvm::LLVMContext& C,
                                         const CodeGenOptions& CGO) {
  return static_cast<CodeGeneratorImpl*>(this)->StartModule(ModuleName, C, CGO);
}

CodeGenerator *clang::CreateLLVMCodeGen(
    DiagnosticsEngine &Diags, llvm::StringRef ModuleName,
    const HeaderSearchOptions &HeaderSearchOpts,
    const PreprocessorOptions &PreprocessorOpts, const CodeGenOptions &CGO,
    llvm::LLVMContext &C, CoverageSourceInfo *CoverageInfo) {
  return new CodeGeneratorImpl(Diags, ModuleName, HeaderSearchOpts,
                               PreprocessorOpts, CGO, C, CoverageInfo);
}
