//------------------------------------------------------------------------------
// CLING - the C++ LLVM-based InterpreterG :)
// author:  Baozeng Ding <sploving1@gmail.com>
// author:  Vassil Vassilev <vasil.georgiev.vasilev@cern.ch>
//
// This file is dual-licensed: you can choose to license it under the University
// of Illinois Open Source License or the GNU Lesser General Public License. See
// LICENSE.TXT for details.
//------------------------------------------------------------------------------

#include "NullDerefProtectionTransformer.h"

#include "cling/Interpreter/Transaction.h"
#include "cling/Utils/AST.h"

#include "clang/AST/ASTContext.h"
#include "clang/AST/Decl.h"
#include "clang/AST/Mangle.h"
#include "clang/AST/StmtVisitor.h"
#include "clang/Basic/SourceLocation.h"
#include "clang/Sema/Lookup.h"

#include <bitset>
#include <map>

using namespace clang;

namespace cling {
  NullDerefProtectionTransformer::NullDerefProtectionTransformer(clang::Sema* S)
    : WrapperTransformer(S) {
  }

  NullDerefProtectionTransformer::~NullDerefProtectionTransformer()
  { }

  // Copied from clad - the clang/opencl autodiff project
  class NodeContext {
  public:
  private:
    typedef llvm::SmallVector<clang::Stmt*, 2> Statements;
    Statements m_Stmts;
  private:
    NodeContext() {};
  public:
    NodeContext(clang::Stmt* s) { m_Stmts.push_back(s); }
    NodeContext(clang::Stmt* s0, clang::Stmt* s1) {
      m_Stmts.push_back(s0);
      m_Stmts.push_back(s1);
    }

    bool isSingleStmt() const { return m_Stmts.size() == 1; }

    clang::Stmt* getStmt() {
      assert(isSingleStmt() && "Cannot get multiple stmts.");
      return m_Stmts.front();
    }
    const clang::Stmt* getStmt() const { return getStmt(); }
    const Statements& getStmts() const {
      return m_Stmts;
    }

    CompoundStmt* wrapInCompoundStmt(clang::ASTContext& C) const {
      assert(!isSingleStmt() && "Must be more than 1");
      llvm::ArrayRef<Stmt*> stmts
        = llvm::makeArrayRef(m_Stmts.data(), m_Stmts.size());
      clang::SourceLocation noLoc;
      return new (C) clang::CompoundStmt(C, stmts, noLoc, noLoc);
    }

    clang::Expr* getExpr() {
      assert(llvm::isa<clang::Expr>(getStmt()) && "Must be an expression.");
      return llvm::cast<clang::Expr>(getStmt());
    }
    const clang::Expr* getExpr() const {
      return getExpr();
    }

    void prepend(clang::Stmt* S) {
      m_Stmts.insert(m_Stmts.begin(), S);
    }

    void append(clang::Stmt* S) {
      m_Stmts.push_back(S);
    }
  };

  class IfStmtInjector : public StmtVisitor<IfStmtInjector, NodeContext> {
  private:
    Sema& m_Sema;
    typedef std::map<clang::FunctionDecl*, std::bitset<32> > decl_map_t;
    std::map<clang::FunctionDecl*, std::bitset<32> > m_NonNullArgIndexs;

    ///\brief Needed for the AST transformations, owned by Sema.
    ///
    ASTContext& m_Context;

    ///\brief cling_runtime_internal_throwIfInvalidPointer cache.
    ///
    LookupResult* m_LookupResult;

  public:
    IfStmtInjector(Sema& S) : m_Sema(S), m_Context(S.getASTContext()),
    m_LookupResult(0) {}
    CompoundStmt* Inject(CompoundStmt* CS) {
      NodeContext result = VisitCompoundStmt(CS);
      return cast<CompoundStmt>(result.getStmt());
    }

    NodeContext VisitStmt(Stmt* S) {
      return NodeContext(S);
    }

    NodeContext VisitCompoundStmt(CompoundStmt* CS) {
      ASTContext& C = m_Sema.getASTContext();
      llvm::SmallVector<Stmt*, 16> stmts;
      for (CompoundStmt::body_iterator I = CS->body_begin(), E = CS->body_end();
           I != E; ++I) {
        NodeContext nc = Visit(*I);
        if (nc.isSingleStmt())
          stmts.push_back(nc.getStmt());
        else
          stmts.append(nc.getStmts().begin(), nc.getStmts().end());
      }

      llvm::ArrayRef<Stmt*> stmtsRef(stmts.data(), stmts.size());
      CompoundStmt* newCS = new (C) CompoundStmt(C, stmtsRef,
                                                 CS->getLBracLoc(),
                                                 CS->getRBracLoc());
      return NodeContext(newCS);
    }

    NodeContext VisitCastExpr(CastExpr* CE) {
      NodeContext result = Visit(CE->getSubExpr());
      return result;
    }

    NodeContext VisitUnaryOperator(UnaryOperator* UnOp) {
      NodeContext result(UnOp);
      if (UnOp->getOpcode() == UO_Deref) {
        result = SynthesizeCheck(UnOp->getLocStart(),
                                 UnOp->getSubExpr());
      }
      return result;
    }

    NodeContext VisitMemberExpr(MemberExpr* ME) {
      NodeContext result(ME);
      if (ME->isArrow()) {
        result = SynthesizeCheck(ME->getLocStart(),
                                       ME->getBase()->IgnoreImplicit());
      }
      return result;
    }

    NodeContext VisitCallExpr(CallExpr* CE) {
      FunctionDecl* FDecl = CE->getDirectCallee();
      NodeContext result(CE);
      if (FDecl && isDeclCandidate(FDecl)) {
        decl_map_t::const_iterator it = m_NonNullArgIndexs.find(FDecl);
        const std::bitset<32>& ArgIndexs = it->second;
        Sema::ContextRAII pushedDC(m_Sema, FDecl);
        for (int index = 0; index < 32; ++index) {
          if (ArgIndexs.test(index)) {
            // Get the argument with the nonnull attribute.
            Expr* Arg = CE->getArg(index);
            result = SynthesizeCheck(Arg->getLocStart(), Arg);
          }
        }
      }
      return result;
    }

  private:
    Stmt* SynthesizeCheck(SourceLocation Loc, Expr* Arg) {
      assert(Arg && "Cannot call with Arg=0");

      if(!m_LookupResult)
        FindAndCacheRuntimeLookupResult();

      Expr* VoidSemaArg = utils::Synthesize::CStyleCastPtrExpr(&m_Sema,
                                                            m_Context.VoidPtrTy,
                                                            (uint64_t)&m_Sema);

      Expr* VoidExprArg = utils::Synthesize::CStyleCastPtrExpr(&m_Sema,
                                                          m_Context.VoidPtrTy,
                                                          (uint64_t)Arg);

      Expr *args[] = {VoidSemaArg, VoidExprArg, Arg};

      Scope* S = m_Sema.getScopeForContext(m_Sema.CurContext);

      CXXScopeSpec CSS;
      Expr* unresolvedLookup
        = m_Sema.BuildDeclarationNameExpr(CSS, *m_LookupResult,
                                         /*ADL*/ false).get();

      Expr* call = m_Sema.ActOnCallExpr(S, unresolvedLookup, Loc,
                                        args, Loc).get();

      TypeSourceInfo* TSI
              = m_Context.getTrivialTypeSourceInfo(Arg->getType(), Loc);
      Expr* castExpr = m_Sema.BuildCStyleCastExpr(Loc, TSI, Loc, call).get();

      return castExpr;
    }

    bool isDeclCandidate(FunctionDecl * FDecl) {
      if (m_NonNullArgIndexs.count(FDecl))
        return true;

      std::bitset<32> ArgIndexs;
      for (specific_attr_iterator<NonNullAttr>
             I = FDecl->specific_attr_begin<NonNullAttr>(),
             E = FDecl->specific_attr_end<NonNullAttr>(); I != E; ++I) {

        NonNullAttr *NonNull = *I;
        for (NonNullAttr::args_iterator i = NonNull->args_begin(),
               e = NonNull->args_end(); i != e; ++i) {
          ArgIndexs.set(*i);
        }
      }

      if (ArgIndexs.any()) {
        m_NonNullArgIndexs.insert(std::make_pair(FDecl, ArgIndexs));
        return true;
      }
      return false;
    }

    void FindAndCacheRuntimeLookupResult() {
      assert(!m_LookupResult && "Called multiple times!?");

      DeclarationName Name
        = &m_Context.Idents.get("cling_runtime_internal_throwIfInvalidPointer");

      SourceLocation noLoc;
      m_LookupResult = new LookupResult(m_Sema, Name, noLoc,
                                        Sema::LookupOrdinaryName,
                                        Sema::ForRedeclaration);
      m_Sema.LookupQualifiedName(*m_LookupResult,
                                 m_Context.getTranslationUnitDecl());
      assert(!m_LookupResult->empty() &&
              "cling_runtime_internal_throwIfInvalidPointer");
    }
  };

  ASTTransformer::Result
  NullDerefProtectionTransformer::Transform(clang::Decl* D) {
    FunctionDecl* FD = dyn_cast<FunctionDecl>(D);
    if (!FD || FD->isFromASTFile())
      return Result(D, true);

    CompoundStmt* CS = dyn_cast_or_null<CompoundStmt>(FD->getBody());
    if (!CS)
      return Result(D, true);

    IfStmtInjector injector(*m_Sema);
    CompoundStmt* newCS = injector.Inject(CS);
    FD->setBody(newCS);
    return Result(FD, true);
  }
} // end namespace cling
