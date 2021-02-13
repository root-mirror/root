// @(#)root/core/meta:$Id$
// Author: Vassil Vassilev   7/10/2012

/*************************************************************************
 * Copyright (C) 1995-2013, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "cling/Interpreter/InterpreterCallbacks.h"

#include <stack>
#include <string>

namespace clang {
   class Decl;
   class DeclarationName;
   class LookupResult;
   class NamespaceDecl;
   class Scope;
   class TagDecl;
   class Token;
   class FileEntry;
   class Module;
}

namespace cling {
   class Interpreter;
   class Transaction;
}

namespace llvm {
   class StringRef;
}

// The callbacks are used to update the list of globals in ROOT.
//
class TClingCallbacks : public cling::InterpreterCallbacks {
private:
   void *fLastLookupCtx = nullptr;
   clang::NamespaceDecl *fROOTSpecialNamespace = nullptr;
   bool fFirstRun = true;
   bool fIsAutoLoading = false;
   bool fIsAutoLoadingRecursively = false;
   bool fIsAutoParsingSuspended = false;
   bool fPPOldFlag = false;
   bool fPPChanged = false;
   bool fIsCodeGening = false;
   bool fIsLoadingModule = false;
   llvm::DenseMap<llvm::StringRef, clang::DeclarationName> m_LoadedModuleFiles;

public:
   TClingCallbacks(cling::Interpreter* interp, bool hasCodeGen);

   ~TClingCallbacks();

   void Initialize();

   void SetAutoLoadingEnabled(bool val = true) { fIsAutoLoading = val; }
   bool IsAutoLoadingEnabled() { return fIsAutoLoading; }

   void SetAutoParsingSuspended(bool val = true) { fIsAutoParsingSuspended = val; }
   bool IsAutoParsingSuspended() { return fIsAutoParsingSuspended; }

   bool LibraryLoadingFailed(const std::string &, const std::string &, bool, bool) override;

   void InclusionDirective(clang::SourceLocation /*HashLoc*/, const clang::Token & /*IncludeTok*/,
                           llvm::StringRef FileName, bool /*IsAngled*/, clang::CharSourceRange /*FilenameRange*/,
                           const clang::FileEntry * /*File*/, llvm::StringRef /*SearchPath*/,
                           llvm::StringRef /*RelativePath*/, const clang::Module * /*Imported*/) override;

   // Preprocessor callbacks used to handle special cases like for example:
   // #include "myMacro.C+"
   //
   bool FileNotFound(llvm::StringRef FileName, llvm::SmallVectorImpl<char> &RecoveryPath) override;

   bool LookupObject(clang::LookupResult &R, clang::Scope *S) override;
   bool LookupObject(const clang::DeclContext *DC, clang::DeclarationName Name) override;
   bool LookupObject(clang::TagDecl *Tag) override;

   // The callback is used to update the list of globals in ROOT.
   //
   void TransactionCommitted(const cling::Transaction &T) override;

   // The callback is used to inform ROOT when cling started code generation.
   //
   void TransactionCodeGenStarted(const cling::Transaction &T) override
   {
      assert(!fIsCodeGening);
      fIsCodeGening = true;
   }

   // The callback is used to inform ROOT when cling finished code generation.
   //
   void TransactionCodeGenFinished(const cling::Transaction &T) override
   {
      assert(fIsCodeGening);
      fIsCodeGening = false;
   }

   // The callback is used to update the list of globals in ROOT.
   //
   void TransactionUnloaded(const cling::Transaction &T) override;

   // The callback is used to clear the autoparsing caches.
   //
   void TransactionRollback(const cling::Transaction &T) override;

   ///\brief A previous definition has been shadowed; invalidate TCling' stored
   /// data about the old (global) decl.
   void DefinitionShadowed(const clang::NamedDecl *D) override;

   // Used to inform client about a new decl read by the ASTReader.
   //
   void DeclDeserialized(const clang::Decl *D) override;

   void LibraryLoaded(const void *dyLibHandle, llvm::StringRef canonicalName) override;
   void LibraryUnloaded(const void *dyLibHandle, llvm::StringRef canonicalName) override;

   void PrintStackTrace() override;

   void *EnteringUserCode() override;
   void ReturnedFromUserCode(void *stateInfo) override;
   void *LockCompilationDuringUserCodeExecution() override;
   void UnlockCompilationDuringUserCodeExecution(void *StateInfo) override;

private:
   bool tryAutoParseInternal(llvm::StringRef Name, clang::LookupResult &R,
                            clang::Scope *S, const clang::FileEntry* FE = 0);
   bool tryFindROOTSpecialInternal(clang::LookupResult &R, clang::Scope *S);
   bool tryResolveAtRuntimeInternal(clang::LookupResult &R, clang::Scope *S);
   bool shouldResolveAtRuntime(clang::LookupResult &R, clang::Scope *S);
   bool tryInjectImplicitAutoKeyword(clang::LookupResult &R, clang::Scope *S);
   bool findInGlobalModuleIndex(clang::DeclarationName Name, bool loadFirstMatchOnly = true);
};
