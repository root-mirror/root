// @(#)root/core/meta:$Id$
// Author: Vassil Vassilev   7/10/2012

/*************************************************************************
 * Copyright (C) 1995-2012, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TClingCallbacks.h"

#include "cling/Interpreter/Interpreter.h"
#include "cling/Interpreter/InterpreterCallbacks.h"
#include "cling/Interpreter/Transaction.h"
#include "cling/Utils/AST.h"

#include "clang/AST/ASTContext.h"
#include "clang/AST/DeclBase.h"
#include "clang/Lex/Preprocessor.h"
#include "clang/Parse/Parser.h"
#include "clang/Sema/Lookup.h"
#include "clang/Sema/Scope.h"

#include "clang/Frontend/CompilerInstance.h"
#include "clang/Lex/PPCallbacks.h"
#include "llvm/Support/FileSystem.h"

using namespace clang;
using namespace cling;

class TObject;

// Functions used to forward calls from code compiled with no-rtti to code 
// compiled with rtti.
extern "C" {
   void TCling__UpdateListsOnCommitted(const cling::Transaction&, Interpreter*);
   void TCling__UpdateListsOnUnloaded(const cling::Transaction&); 
   TObject* TCling__GetObjectAddress(const char *Name, void *&LookupCtx);
   Decl* TCling__GetObjectDecl(TObject *obj);
   int TCling__AutoLoadCallback(const char* className);
   void TCling__UpdateListsOnDeclDeserialized(const clang::Decl*);
   int TCling__CompileMacro(const char *fileName, const char *options);
   std::string TCling__SplitAclicMode(const char* fileName, std::string &mode,
                                      std::string &args, std::string &io);
}

// Preprocessor callbacks used to handle special cases
// like for example: #include "myMacro.C+"
//
class TPPClingCallbacks : public PPCallbacks {
private:
   cling::Interpreter *fInterpreter;
   Preprocessor *fPreprocessor;
   bool fOldFlag;
   bool fChanged;
public:
   TPPClingCallbacks(cling::Interpreter *inter, Preprocessor *PP) :
         fInterpreter(inter), fPreprocessor(PP), fOldFlag(false),
         fChanged(false) { }
   ~TPPClingCallbacks() { }

   virtual bool FileNotFound(llvm::StringRef FileName,
                             llvm::SmallVectorImpl<char> &RecoveryPath) {
      // Method called via Callbacks->FileNotFound(Filename, RecoveryPath)
      // in Preprocessor::HandleIncludeDirective(), initally allowing to
      // change the include path, and allowing us to compile code via ACLiC
      // when specifying #include "myfile.C+", and suppressing the preprocessor
      // error message:
      // input_line_23:1:10: fatal error: 'myfile.C+' file not found
      if ((!fPreprocessor) || (!fInterpreter))
         return false;

      // remove any trailing "\n
      std::string filename(FileName.str().substr(0,
                           FileName.str().find_last_of('"')));
      std::string mode, arguments, io;
      // extract the filename and ACliC mode
      std::string fname = TCling__SplitAclicMode(filename.c_str(), mode,
                                                 arguments, io);
      if (mode.length() > 0) {
         if (llvm::sys::fs::exists(fname)) {
            // format the CompileMacro() option string
            std::string options = "k";
            if (mode.find("++") != std::string::npos) options += "f";
            if (mode.find("g")  != std::string::npos) options += "g";
            if (mode.find("O")  != std::string::npos) options += "O";

            // Save state of the preprocessor
            Preprocessor::CleanupAndRestoreCacheRAII cleanupRAII(*fPreprocessor);
            Parser& P = const_cast<Parser&>(fInterpreter->getParser());
            Parser::ParserCurTokRestoreRAII savedCurToken(P);

            // After we have saved the token reset the current one to
            // something which is safe (semi colon usually means empty decl)
            Token& Tok = const_cast<Token&>(P.getCurToken());
            Tok.setKind(tok::semi);

            int retcode = TCling__CompileMacro(fname.c_str(), options.c_str());
            if (retcode) {
               // complation was successful, let's remember the original
               // preprocessor "include not found" error suppression flag
               if (!fChanged)
                  fOldFlag = fPreprocessor->GetSuppressIncludeNotFoundError();
               fPreprocessor->SetSuppressIncludeNotFoundError(true);
               fChanged = true;
            }
            return true;
         }
      }
      if (fChanged) {
         // restore the original preprocessor "include not found" error
         // suppression flag
         fPreprocessor->SetSuppressIncludeNotFoundError(fOldFlag);
         fChanged = false;
      }
      return false;
   }

};

TClingCallbacks::TClingCallbacks(cling::Interpreter* interp) 
   : InterpreterCallbacks(interp),
     fLastLookupCtx(0), fROOTSpecialNamespace(0), fFirstRun(true), 
     fIsAutoloading(false), fIsAutoloadingRecursively(false) {
   Transaction* T = 0;
   m_Interpreter->declare("namespace __ROOT_SpecialObjects{}", &T);
   fROOTSpecialNamespace = dyn_cast<NamespaceDecl>(T->getFirstDecl().getSingleDecl());
   // Add PreProcessor callback implementing FileNotFound(Filename, RecoveryPath)
   // in order to properly handle #include "myMacro.C+"
   Preprocessor &PP = interp->getCI()->getPreprocessor();
   PP.addPPCallbacks(new TPPClingCallbacks(interp, &PP));
}

//pin the vtable here
TClingCallbacks::~TClingCallbacks() {}

// On a failed lookup we have to try to more things before issuing an error.
// The symbol might need to be loaded by ROOT's autoloading mechanism or
// it might be a ROOT special object. 
// 
// Try those first and if still failing issue the diagnostics.
//
// returns true when a declaration is found and no error should be emitted.
//
bool TClingCallbacks::LookupObject(LookupResult &R, Scope *S) {

   if (tryAutoloadInternal(R, S))
      return true; // happiness.

   // If the autoload wasn't successful try ROOT specials.
   if (tryFindROOTSpecialInternal(R, S))
      return true;

   // For backward-compatibility with CINT we must support stmts like:
   // x = 4; y = new MyClass();
   // I.e we should "inject" a C++11 auto keyword in front of "x" and "y"
   // This has to have higher precedence than the dynamic scopes. It is claimed
   // that if one assigns to a name and the lookup of that name fails if *must*
   // auto keyword must be injected and the stmt evaluation must not be delayed
   // until runtime.
   // For now supported only at the prompt.
   if (tryInjectImplicitAutoKeyword(R, S)) {
      return true;
   }

   // Finally try to resolve this name as a dynamic name, i.e delay its 
   // resolution for runtime.
   return tryResolveAtRuntimeInternal(R, S);
}

// The symbol might be defined in the ROOT class autoloading map so we have to
// try to autoload it first and do secondary lookup to try to find it.
//
// returns true when a declaration is found and no error should be emitted.
//
bool TClingCallbacks::tryAutoloadInternal(LookupResult &R, Scope *S) {
   Sema &SemaR = m_Interpreter->getSema();
   DeclarationName Name = R.getLookupName();

   // Try to autoload first if autoloading is enabled
   if (IsAutoloadingEnabled()) {
     // Avoid tail chasing.
     if (fIsAutoloadingRecursively)
       return false;

     // We should try autoload only for special lookup failures. 
     Sema::LookupNameKind kind = R.getLookupKind();
     if (!(kind == Sema::LookupTagName || kind == Sema::LookupOrdinaryName)) 
        return false;

     fIsAutoloadingRecursively = true;

     bool lookupSuccess = false;
     if (getenv("ROOT_MODULES")) {
        if (TCling__AutoLoadCallback(Name.getAsString().c_str())) {
           lookupSuccess = SemaR.LookupName(R, S);
        }
     }
     else {
        // Save state of the PP
        ASTContext& C = SemaR.getASTContext();
        Preprocessor &PP = SemaR.getPreprocessor();
        Preprocessor::CleanupAndRestoreCacheRAII cleanupRAII(PP);
        Parser& P = const_cast<Parser&>(m_Interpreter->getParser());
        Parser::ParserCurTokRestoreRAII savedCurToken(P);
        // After we have saved the token reset the current one to something which 
        // is safe (semi colon usually means empty decl)
        Token& Tok = const_cast<Token&>(P.getCurToken());
        Tok.setKind(tok::semi);

        bool oldSuppressDiags = SemaR.getDiagnostics().getSuppressAllDiagnostics();
        SemaR.getDiagnostics().setSuppressAllDiagnostics();

        // We can't PushDeclContext, because we go up and the routine that pops 
        // the DeclContext assumes that we drill down always.
        // We have to be on the global context. At that point we are in a 
        // wrapper function so the parent context must be the global.
        Sema::ContextAndScopeRAII pushedDCAndS(SemaR, C.getTranslationUnitDecl(), 
                                               SemaR.TUScope);

        if (TCling__AutoLoadCallback(Name.getAsString().c_str())) {
            pushedDCAndS.pop();
           cleanupRAII.pop();
           lookupSuccess = SemaR.LookupName(R, S);
        }
 
        SemaR.getDiagnostics().setSuppressAllDiagnostics(oldSuppressDiags);
     }

     fIsAutoloadingRecursively = false;
   
     if (lookupSuccess)
       return true;
   }

   return false;
}

// If cling cannot find a name it should ask ROOT before it issues an error.
// If ROOT knows the name then it has to create a new variable with that name
// and type in dedicated for that namespace (eg. __ROOT_SpecialObjects).
// For example if the interpreter is looking for h in h-Draw(), this routine
// will create
// namespace __ROOT_SpecialObjects {
//   THist* h = (THist*) the_address;
// }
//
// Later if h is called again it again won't be found by the standart lookup
// because it is in our hidden namespace (nobody should do using namespace 
// __ROOT_SpecialObjects). It caches the variable declarations and their
// last address. If the newly found decl with the same name (h) has different
// address than the cached one it goes directly at the address and updates it.
//
// returns true when declaration is found and no error should be emitted.
//
bool TClingCallbacks::tryFindROOTSpecialInternal(LookupResult &R, Scope *S) {
   // User must be able to redefine the names that come from a file.
   if (R.isForRedeclaration())
      return false;
   // If there is a result abort.
   if (!R.empty())
      return false;

   Sema &SemaR = m_Interpreter->getSema();
   ASTContext& C = SemaR.getASTContext();
   Preprocessor &PP = SemaR.getPreprocessor();
   DeclContext *CurDC = SemaR.CurContext;
   DeclarationName Name = R.getLookupName();

   // Make sure that the failed lookup comes from a function body.
   if(!CurDC || !CurDC->isFunctionOrMethod())
      return false;

   // Save state of the PP, because TCling__GetObjectAddress may induce nested
   // lookup.
   Preprocessor::CleanupAndRestoreCacheRAII cleanupPPRAII(PP);
   TObject *obj = TCling__GetObjectAddress(Name.getAsString().c_str(), 
                                           fLastLookupCtx);
   cleanupPPRAII.pop(); // force restroing the cache

   if (obj) {

#if defined(R__MUST_REVISIT)
#if R__MUST_REVISIT(6,2)
      // Register the address in TCling::fgSetOfSpecials
      // to speed-up the execution of TCling::RecursiveRemove when 
      // the object is not a special.
      // See http://root.cern.ch/viewvc/trunk/core/meta/src/TCint.cxx?view=log#rev18109
      if (!fgSetOfSpecials) {
         fgSetOfSpecials = new std::set<TObject*>;
      }
      ((std::set<TObject*>*)fgSetOfSpecials)->insert((TObject*)*obj);
#endif
#endif

     VarDecl *VD = cast_or_null<VarDecl>(utils::Lookup::Named(&SemaR, Name, 
                                                        fROOTSpecialNamespace));
      if (VD) {
         //TODO: Check for same types.

         TObject **address = (TObject**)m_Interpreter->getAddressOfGlobal(VD);
         // Since code was generated already we cannot rely on the initializer 
         // of the decl in the AST, however we will update that init so that it
         // will be easier while debugging.
         CStyleCastExpr *CStyleCast = cast<CStyleCastExpr>(VD->getInit());
         Expr* newInit = utils::Synthesize::IntegerLiteralExpr(C, (uint64_t)obj);
         CStyleCast->setSubExpr(newInit);

         // The actual update happens here, directly in memory.
         *address = obj;
      }
      else {
         // Save state of the PP
         Preprocessor::CleanupAndRestoreCacheRAII cleanupRAII(PP);

         const Decl *TD = TCling__GetObjectDecl(obj);
         // We will declare the variable as pointer.
         QualType QT = C.getPointerType(C.getTypeDeclType(cast<TypeDecl>(TD)));
         
         VD = VarDecl::Create(C, fROOTSpecialNamespace, SourceLocation(), 
                              SourceLocation(), Name.getAsIdentifierInfo(), QT,
                              /*TypeSourceInfo*/0, SC_None);
         // Build an initializer
         Expr* Init 
           = utils::Synthesize::CStyleCastPtrExpr(&SemaR, QT, (uint64_t)obj);
         // Register the decl in our hidden special namespace
         VD->setInit(Init);
         fROOTSpecialNamespace->addDecl(VD);

         cling::CompilationOptions CO;
         CO.DeclarationExtraction = 0;
         CO.ValuePrinting = CompilationOptions::VPDisabled;
         CO.ResultEvaluation = 0;
         CO.DynamicScoping = 0;
         CO.Debug = 0;
         CO.CodeGeneration = 1;

         cling::Transaction T(CO, VD->getASTContext());
         T.append(VD);
         T.setState(cling::Transaction::kCompleted);

         m_Interpreter->emitAllDecls(&T);
         assert(T.getState() == Transaction::kCommitted
                && "Compilation should never fail!");
      }
      assert(VD && "Cannot be null!");
      R.addDecl(VD);
      return true;
   }

   return false;
}

bool TClingCallbacks::tryResolveAtRuntimeInternal(LookupResult &R, Scope *S) {
   if (!shouldResolveAtRuntime(R, S))
      return false;

   DeclarationName Name = R.getLookupName();
   IdentifierInfo* II = Name.getAsIdentifierInfo();
   SourceLocation Loc = R.getNameLoc();
   Sema& SemaRef = R.getSema();
   ASTContext& C = SemaRef.getASTContext();
   DeclContext* DC = C.getTranslationUnitDecl();
   assert(DC && "Must not be null.");
   VarDecl* Result = VarDecl::Create(C, DC, Loc, Loc, II, C.DependentTy,
                                     /*TypeSourceInfo*/0, SC_None);

   // Annotate the decl to give a hint in cling. FIXME: Current implementation
   // is a gross hack, because TClingCallbacks shouldn't know about 
   // EvaluateTSynthesizer at all!
    
   SourceRange invalidRange;
   Result->addAttr(new (C) AnnotateAttr(invalidRange, C, "__ResolveAtRuntime"));
   if (Result) {
      // Here we have the scope but we cannot do Sema::PushDeclContext, because
      // on pop it will try to go one level up, which we don't want.
      Sema::ContextRAII pushedDC(SemaRef, DC);
      R.addDecl(Result);
      SemaRef.PushOnScopeChains(Result, SemaRef.TUScope, /*Add to ctx*/true);
      // Say that we can handle the situation. Clang should try to recover
      return true;
   }
   // We cannot handle the situation. Give up
   return false;
}

bool TClingCallbacks::shouldResolveAtRuntime(LookupResult& R, Scope* S) {
   if (m_IsRuntime)
     return false;

   if (R.getLookupKind() != Sema::LookupOrdinaryName) 
      return false;

   if (R.isForRedeclaration()) 
      return false;

   if (!R.empty())
      return false;

   // FIXME: Figure out better way to handle:
   // C++ [basic.lookup.classref]p1:
   //   In a class member access expression (5.2.5), if the . or -> token is
   //   immediately followed by an identifier followed by a <, the
   //   identifier must be looked up to determine whether the < is the
   //   beginning of a template argument list (14.2) or a less-than operator.
   //   The identifier is first looked up in the class of the object
   //   expression. If the identifier is not found, it is then looked up in
   //   the context of the entire postfix-expression and shall name a class
   //   or function template.
   //
   // We want to ignore object(.|->)member<template>
   if (R.getSema().PP.LookAhead(0).getKind() == tok::less)
      // TODO: check for . or -> in the cached token stream
      return false;
   
   for (Scope* DepScope = S; DepScope; DepScope = DepScope->getParent()) {
      if (DeclContext* Ctx = static_cast<DeclContext*>(DepScope->getEntity())) {
         if (!Ctx->isDependentContext())
            // For now we support only the prompt.
            if (isa<FunctionDecl>(Ctx))
               return true;
      }
   }

   return false;
}

bool TClingCallbacks::tryInjectImplicitAutoKeyword(LookupResult &R, Scope *S) {
   // Make sure that the failed lookup comes the prompt. Currently, we support
   // only the prompt.

   // Should be disabled with the dynamic scopes.
   if (m_IsRuntime)
      return false;

   if (R.isForRedeclaration()) 
      return false;

   if (R.getLookupKind() != Sema::LookupOrdinaryName) 
      return false;

   if (!isa<FunctionDecl>(R.getSema().CurContext))
      return false;

   Sema& SemaRef = R.getSema();
   ASTContext& C = SemaRef.getASTContext();
   DeclContext* DC = C.getTranslationUnitDecl();
   assert(DC && "Must not be null.");


   Preprocessor& PP = R.getSema().getPreprocessor();
   //Preprocessor::CleanupAndRestoreCacheRAII cleanupRAII(PP);
   //PP.EnableBacktrackAtThisPos();
   if (PP.LookAhead(0).isNot(tok::equal)) {
      //PP.Backtrack();
      return false;
   }
   //PP.CommitBacktrackedTokens();
   //cleanupRAII.pop();
   DeclarationName Name = R.getLookupName();
   IdentifierInfo* II = Name.getAsIdentifierInfo();
   SourceLocation Loc = R.getNameLoc();
   VarDecl* Result = VarDecl::Create(C, DC, Loc, Loc, II, 
                                     C.getAutoType(QualType()),
                                     /*TypeSourceInfo*/0, SC_None);

   // Annotate the decl to give a hint in cling. FIXME: Current implementation
   // is a gross hack, because TClingCallbacks shouldn't know about 
   // EvaluateTSynthesizer at all!
    
   SourceRange invalidRange;
   Result->addAttr(new (C) AnnotateAttr(invalidRange, C, "__Auto"));

   if (Result) {
      R.addDecl(Result);
      // Here we have the scope but we cannot do Sema::PushDeclContext, because
      // on pop it will try to go one level up, which we don't want.
      Sema::ContextRAII pushedDC(SemaRef, DC);
      SemaRef.PushOnScopeChains(Result, SemaRef.TUScope, /*Add to ctx*/true);
      // Say that we can handle the situation. Clang should try to recover
      return true;
   }
   // We cannot handle the situation. Give up.
   return false;
}

void TClingCallbacks::Initialize(const ASTContext& Ctx) {
   // Replay existing decls from the AST.
   if (fFirstRun) {
      // Before setting up the callbacks register what cling have seen during init.
      cling::Transaction TPrev((cling::CompilationOptions(), Ctx));
      TPrev.append(Ctx.getTranslationUnitDecl());
      TCling__UpdateListsOnCommitted(TPrev, m_Interpreter);

      fFirstRun = false;
   }
}

// The callback is used to update the list of globals in ROOT.
//
void TClingCallbacks::TransactionCommitted(const Transaction &T) {
   // Even empty transactions must go through; any transaction even empty
   // will flush the deserialized decls into Meta.
   //if (!T.size())
   //   return;
   if (fFirstRun && T.size())
      Initialize(T.getASTContext());

   TCling__UpdateListsOnCommitted(T, m_Interpreter);
}

// The callback is used to update the list of globals in ROOT.
//
void TClingCallbacks::TransactionUnloaded(const Transaction &T) {
   if (!T.size())
      return;

   TCling__UpdateListsOnUnloaded(T);
}

void TClingCallbacks::DeclDeserialized(const clang::Decl* D) {
   TCling__UpdateListsOnDeclDeserialized(D);
}
