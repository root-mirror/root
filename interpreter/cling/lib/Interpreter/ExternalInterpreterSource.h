//--------------------------------------------------------------------*- C++ -*-
// CLING - the C++ LLVM-based InterpreterG :)
// author:  Elisavet Sakellari <elisavet.sakellari@cern.ch>
//
// This file is dual-licensed: you can choose to license it under the University
// of Illinois Open Source License or the GNU Lesser General Public License. See
// LICENSE.TXT for details.
//------------------------------------------------------------------------------

#ifndef CLING_EXTERNAL_INTERPRETER_SOURCE
#define CLING_EXTERNAL_INTERPRETER_SOURCE

#include "clang/AST/ExternalASTSource.h"

#include <string>
#include <map>

namespace clang {
  class ASTContext;
  class ASTImporter;
  class Decl;
  class DeclContext;
  class DeclarationName;
  class ExternalASTSource;
  class NamedDecl;
  class Sema;
}

namespace cling {
  class Interpreter;
}

namespace cling {

    class ExternalInterpreterSource : public clang::ExternalASTSource {

      private:
        const cling::Interpreter *m_ParentInterpreter;
        cling::Interpreter *m_ChildInterpreter;

        clang::Sema *m_Sema;

        ///\brief We keep a mapping between the imported DeclContexts
        /// and the original ones from of the first Interpreter.
        /// Key: imported DeclContext
        /// Value: original DeclContext
        ///
        std::map<const clang::DeclContext *, clang::DeclContext *> m_ImportedDeclContexts;

        ///\brief A map for all the imported Decls (Contexts)
        /// according to their names.
        /// Key: Name of the Decl(Context) as a string.
        /// Value: The DeclarationName of this Decl(Context) is the one
        /// that comes from the first Interpreter.
        ///
        std::map <clang::DeclarationName, clang::DeclarationName > m_ImportedDecls;

      public:
        ExternalInterpreterSource(const cling::Interpreter *parent,
                                  cling::Interpreter *child);
        virtual ~ExternalInterpreterSource();

        void completeVisibleDeclsMap(const clang::DeclContext *DC) override;

        bool FindExternalVisibleDeclsByName(
                              const clang::DeclContext *childCurrentDeclContext,
                              clang::DeclarationName childDeclName) override;

        void InitializeSema(clang::Sema &S) { m_Sema = &S; }

        void ForgetSema() { m_Sema = nullptr; }

        bool Import(clang::DeclContext::lookup_result lookupResult,
                    clang::ASTContext &parentASTContext,
                    clang::ASTContext &childASTContext,
                    const clang::DeclContext *childCurrentDeclContext,
                    clang::DeclarationName &childDeclName,
                    clang::DeclarationName &parentDeclName);

        void ImportDeclContext(clang::DeclContext *declContextToImport,
                               clang::ASTImporter &importer,
                               clang::DeclarationName &childDeclName,
                               clang::DeclarationName &parentDeclName,
                               const clang::DeclContext *childCurrentDeclContext);

        void ImportDecl(clang::Decl *declToImport,
                        clang::ASTImporter &importer,
                        clang::DeclarationName &childDeclName,
                        clang::DeclarationName &parentDeclName,
                        const clang::DeclContext *childCurrentDeclContext);

        void addToImportedDecls(clang::DeclarationName child,
                                clang::DeclarationName parent) {
          m_ImportedDecls[child] = parent;
        }

        void addToImportedDeclContexts(clang::DeclContext *child,
                              clang::DeclContext *parent) {
          m_ImportedDeclContexts[child] = parent;
        }
    };
} // end namespace cling

#endif //CLING_ASTIMPORTSOURCE_H
