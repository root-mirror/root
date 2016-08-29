//--------------------------------------------------------------------*- C++ -*-
// CLING - the C++ LLVM-based InterpreterG :)
// author:  Vassil Vassilev <vvasilev@cern.ch>
//
// This file is dual-licensed: you can choose to license it under the University
// of Illinois Open Source License or the GNU Lesser General Public License. See
// LICENSE.TXT for details.
//------------------------------------------------------------------------------

#ifndef CLING_DECL_UNLOADER
#define CLING_DECL_UNLOADER

#include "cling/Interpreter/Transaction.h"

#include "clang/AST/DeclVisitor.h"

#ifndef CLING_CONSTEXPR_FUNC
# if defined(_MSC_VER) && (_MSC_VER < 1900)
#  define CLING_CONSTEXPR_FUNC
# else
#  define CLING_CONSTEXPR_FUNC constexpr
# endif
#endif

namespace clang {

  class CodeGenerator;
  class GlobalDecl;

  using cling::Transaction;

  ///\brief The class does the actual work of removing a declaration and
  /// resetting the internal structures of the compiler
  ///
  class DeclUnloader : public DeclVisitor<DeclUnloader, bool> {
  private:
    typedef llvm::DenseSet<FileID> FileIDs;

    ///\brief The Sema object being unloaded (contains the AST as well).
    ///
    Sema* m_Sema;

    ///\brief The clang code generator, being recovered.
    ///
    clang::CodeGenerator* m_CodeGen;

    ///\brief The current transaction being unloaded.
    ///
    const Transaction* m_CurTransaction;

    ///\brief Unloaded declaration contains a SourceLocation, representing a
    /// place in the file where it was seen. Clang caches that file and even if
    /// a declaration is removed and the file is edited we hit the cached entry.
    /// This ADT keeps track of the files from which the unloaded declarations
    /// came from so that in the end they could be removed from clang's cache.
    ///
    FileIDs m_FilesToUncache;

  public:
    DeclUnloader(Sema* S, clang::CodeGenerator* CG, const Transaction* T)
      : m_Sema(S), m_CodeGen(CG), m_CurTransaction(T) { }
    ~DeclUnloader();

    ///\brief Interface with nice name, forwarding to Visit.
    ///
    ///\param[in] D - The declaration to forward.
    ///\returns true on success.
    ///
    bool UnloadDecl(Decl* D) { return Visit(D); }

    ///\brief If it falls back in the base class just remove the declaration
    /// only from the declaration context.
    /// @param[in] D - The declaration to be removed.
    ///
    ///\returns true on success.
    ///
    bool VisitDecl(Decl* D);

    ///\brief Removes the declaration from the lookup chains and from the
    /// declaration context.
    /// @param[in] ND - The declaration to be removed.
    ///
    ///\returns true on success.
    ///
    bool VisitNamedDecl(NamedDecl* ND);

    ///\brief Removes the declaration from Sema's unused decl registry
    /// @param[in] DD - The declaration to be removed.
    ///
    ///\returns true on success.
    ///
    bool VisitDeclaratorDecl(DeclaratorDecl* DD);

    ///\brief Removes a using shadow declaration, created in the cases:
    ///\code
    /// namespace A {
    ///   void foo();
    /// }
    /// namespace B {
    ///   using A::foo; // <- a UsingDecl
    ///                 // Also creates a UsingShadowDecl for A::foo() in B
    /// }
    ///\endcode
    ///\param[in] USD - The declaration to be removed.
    ///
    ///\returns true on success.
    ///
    bool VisitUsingShadowDecl(UsingShadowDecl* USD);

    ///\brief Removes a typedef name decls. A base class for TypedefDecls and
    /// TypeAliasDecls.
    ///\param[in] TND - The declaration to be removed.
    ///
    ///\returns true on success.
    ///
    bool VisitTypedefNameDecl(TypedefNameDecl* TND);

    ///\brief Removes the declaration from the lookup chains and from the
    /// declaration context and it rebuilds the redeclaration chain.
    /// @param[in] VD - The declaration to be removed.
    ///
    ///\returns true on success.
    ///
    bool VisitVarDecl(VarDecl* VD);

    ///\brief Removes the declaration from the lookup chains and from the
    /// declaration context and it rebuilds the redeclaration chain.
    /// @param[in] FD - The declaration to be removed.
    ///
    ///\returns true on success.
    ///
    bool VisitFunctionDecl(FunctionDecl* FD);

    ///\brief Specialize the removal of constructors due to the fact the we need
    /// the constructor type (aka CXXCtorType). The information is located in
    /// the CXXConstructExpr of usually VarDecls.
    /// See clang::CodeGen::CodeGenFunction::EmitCXXConstructExpr.
    ///
    /// What we will do instead is to brute-force and try to remove from the
    /// llvm::Module all ctors of this class with all the types.
    ///
    ///\param[in] CXXCtor - The declaration to be removed.
    ///
    ///\returns true on success.
    ///
    bool VisitCXXConstructorDecl(CXXConstructorDecl* CXXCtor);

    ///\brief Specialize the removal of destructors due to the fact the we need
    /// the to erase the dtor decl and the deleting operator.
    ///
    /// We will brute-force and try to remove from the llvm::Module all dtors of
    /// this class with all the types.
    ///
    ///\param[in] CXXDtor - The declaration to be removed.
    ///
    ///\returns true on success.
    ///
    bool VisitCXXDestructorDecl(CXXDestructorDecl* CXXDtor);

    ///\brief Removes the DeclCotnext and its decls.
    /// @param[in] DC - The declaration to be removed.
    ///
    ///\returns true on success.
    ///
    bool VisitDeclContext(DeclContext* DC);

    ///\brief Removes the namespace.
    /// @param[in] NSD - The declaration to be removed.
    ///
    ///\returns true on success.
    ///
    bool VisitNamespaceDecl(NamespaceDecl* NSD);

    ///\brief Removes a Tag (class/union/struct/enum). Most of the other
    /// containers fall back into that case.
    /// @param[in] TD - The declaration to be removed.
    ///
    ///\returns true on success.
    ///
    bool VisitTagDecl(TagDecl* TD);

    ///\brief Removes a RecordDecl. We shouldn't remove the implicit class
    /// declaration.
    ///\param[in] RD - The declaration to be removed.
    ///
    ///\returns true on success.
    ///
    bool VisitRecordDecl(RecordDecl* RD);

    ///\brief Remove the macro from the Preprocessor.
    /// @param[in] MD - The MacroDirectiveInfo containing the IdentifierInfo and
    ///                MacroDirective to forward.
    ///
    ///\returns true on success.
    ///
    bool VisitMacro(const Transaction::MacroDirectiveInfo MD);

    ///@name Templates
    ///@{

    ///\brief Removes template from the redecl chain. Templates are
    /// redeclarables also.
    /// @param[in] R - The declaration to be removed.
    ///
    ///\returns true on success.
    ///
    bool VisitRedeclarableTemplateDecl(RedeclarableTemplateDecl* R);


    ///\brief Removes the declaration clang's internal structures. This case
    /// looks very much to VisitFunctionDecl, but FunctionTemplateDecl doesn't
    /// derive from FunctionDecl and thus we need to handle it 'by hand'.
    /// @param[in] FTD - The declaration to be removed.
    ///
    ///\returns true on success.
    ///
    bool VisitFunctionTemplateDecl(FunctionTemplateDecl* FTD);

    ///\brief Removes a class template declaration from clang's internal
    /// structures.
    /// @param[in] CTD - The declaration to be removed.
    ///
    ///\returns true on success.
    ///
    bool VisitClassTemplateDecl(ClassTemplateDecl* CTD);

    ///\brief Removes a class template specialization declaration from clang's
    /// internal structures.
    /// @param[in] CTSD - The declaration to be removed.
    ///
    ///\returns true on success.
    ///
    bool VisitClassTemplateSpecializationDecl(ClassTemplateSpecializationDecl*
                                              CTSD);

    ///@}

    void MaybeRemoveDeclFromModule(GlobalDecl& GD) const;

    /// @name Helpers
    /// @{

    ///\brief Interface with nice name, forwarding to Visit.
    ///
    ///\param[in] MD - The MacroDirectiveInfo containing the IdentifierInfo and
    ///                MacroDirective to forward.
    ///\returns true on success.
    ///
    bool UnloadMacro(Transaction::MacroDirectiveInfo MD) {
      return VisitMacro(MD);
    }
    /// @}

  private:
    ///\brief Function that collects the files which we must reread from disk.
    ///
    /// For example: We must uncache the cached include, which brought a
    /// declaration or a macro directive definition in the AST.
    ///\param[in] Loc - The source location of the unloaded declaration.
    ///
    void CollectFilesToUncache(SourceLocation Loc);

    CLING_CONSTEXPR_FUNC static bool isDefinition(void*) { return false; }
    static bool isDefinition(TagDecl* R);

    static void resetDefinitionData(void*) {
      llvm_unreachable("resetDefinitionData on non-cxx record declaration");
    }

    static void resetDefinitionData(TagDecl *decl);

    template<typename DeclT>
    static void removeRedeclFromChain(DeclT* R);

    static void removeRedeclFromChain(void*) {
      llvm_unreachable("setLatestDeclImpl on non-redeclarable declaration");
    }

    template <typename T>
    bool VisitRedeclarable(clang::Redeclarable<T>* R, DeclContext* DC);
  };
} // namespace clang

namespace cling {
  using clang::DeclUnloader;

  /// \brief Unload a Decl from the AST, but not from CodeGen or Module.
  inline bool UnloadDecl(clang::Sema* S, clang::Decl* D) {
    DeclUnloader Unloader(S, nullptr, nullptr);
    return Unloader.UnloadDecl(D);
  }

  /// \brief Unload a Decl from the AST and CodeGen, but not from the Module.
  inline bool UnloadDecl(clang::Sema* S, clang::CodeGenerator* CG, clang::Decl* D) {
    DeclUnloader Unloader(S, CG, nullptr);
    return Unloader.UnloadDecl(D);
  }
} // namespace cling

#endif // CLING_DECL_UNLOADER
